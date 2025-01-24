pub mod align_processor;
pub mod bam_writer;
pub mod params;
pub mod mapping_ext;
use std::{
    cmp,
    collections::HashMap,
    ops::{Deref, DerefMut},
    thread,
};

use crossbeam::channel::{Receiver, Sender};
use gskits::{dna::reverse_complement, ds::ReadInfo, fastx_reader::fasta_reader::FastaFileReader};
use mapping_ext::MappingExt;
use minimap2::{Aligner, Built, Mapping};
use params::{
    AlignParams, IndexParams, InputFilterParams, MapParams, OupParams, TOverrideAlignerParam,
};
use rust_htslib::bam::{
    record::{Aux, AuxArray, Cigar, CigarString},
    Read,
};

pub type BamRecord = rust_htslib::bam::record::Record;
pub type BamWriter = rust_htslib::bam::Writer;
pub type BamReader = rust_htslib::bam::Reader;

pub struct AlignResult {
    pub records: Vec<BamRecord>,
}

pub struct NoMemLeakAligner(pub Aligner<Built>);
impl Deref for NoMemLeakAligner {
    type Target = Aligner<Built>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for NoMemLeakAligner {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Drop for NoMemLeakAligner {
    fn drop(&mut self) {
        // let idx = self.idx.take().unwrap();
        // // unsafe {
        // //     mm_idx_destroy(*idx.as_ref());
        // // }
    }
}

impl From<Aligner<Built>> for NoMemLeakAligner {
    fn from(value: Aligner<Built>) -> Self {
        Self(value)
    }
}

pub fn build_aligner(
    preset: &str,
    index_args: &IndexParams,
    map_args: &MapParams,
    align_args: &AlignParams,
    oup_args: &OupParams,
    targets: &Vec<ReadInfo>,
) -> Vec<NoMemLeakAligner> {
    let aligners = thread::scope(|s| {
        let mut handles = vec![];
        for target in targets {
            let hd = s.spawn(|| {
                let aligner = Aligner::builder();

                let mut aligner = match preset {
                    "map-ont" => aligner.map_ont(),
                    "map-pb" => aligner.map_pb(),
                    "map-hifi" => aligner.map_hifi(),

                    pre => panic!("not implemented yet {}", pre),
                };

                index_args.modify_aligner(&mut aligner);
                map_args.modify_aligner(&mut aligner);
                align_args.modify_aligner(&mut aligner);
                oup_args.modify_aligner(&mut aligner);

                let mut aligner = aligner
                    .with_index_threads(4)
                    .with_cigar()
                    .with_sam_out()
                    .with_sam_hit_only()
                    .with_seq_and_id(target.seq.as_bytes(), target.name.as_bytes())
                    .unwrap();

                // https://github.com/lh3/minimap2/blob/618d33515e5853c4576d5a3d126fdcda28f0e8a4/options.c#L43
                aligner.mapopt.best_n = 5; // top best_n chains are subjected to DP alignment
                // aligner.mapopt.pri_ratio = 0.1; // secondary-to-primary score ratio to output secondary mappings


                aligner
            });
            handles.push(hd);
        }

        handles
            .into_iter()
            .map(|hd| NoMemLeakAligner(hd.join().unwrap()))
            .collect::<Vec<_>>()
    });

    aligners
}

pub fn query_seq_sender(
    filenames: &Vec<String>,
    sender: Sender<ReadInfo>,
    input_filter_params: &InputFilterParams,
) {
    let mut file_idx = 0;
    for filename in filenames {
        let qname_suffix = if filenames.len() == 1 {
            None
        } else {
            Some(format!("__{}", file_idx))
        };
        if filename.ends_with("fa") || filename.ends_with("fasta") || filename.ends_with("fna") {
            let fasta_reader = FastaFileReader::new(filename.clone());

            for mut record in fasta_reader {
                if let Some(suffix) = &qname_suffix {
                    record.name.push_str(suffix);
                }
                sender.send(record).unwrap();
            }
        } else if filename.ends_with("bam") {
            let mut bam_h = BamReader::from_path(filename).unwrap();
            bam_h.set_threads(4).unwrap();
            for record in bam_h.records() {
                let record = record.unwrap();
                if input_filter_params.valid(&record) {
                    sender
                        .send(ReadInfo::from_bam_record(
                            &record,
                            qname_suffix.as_ref().map(|v| v.as_str()),
                        ))
                        .unwrap();
                }
            }
        } else if filename.ends_with("fq") || filename.ends_with("fastq") {
            let fa_iter = FastaFileReader::new(filename.to_string());
            for mut record in fa_iter {
                if let Some(suffix) = &qname_suffix {
                    record.seq.push_str(suffix);
                }
                sender.send(record).unwrap();
            }
        } else {
            panic!(
                "invalid file format {}. bam/fa/fasta/fna supported",
                filename
            );
        }

        file_idx += 1;
    }
}

pub fn align_worker(
    query_record_recv: Receiver<gskits::ds::ReadInfo>,
    align_res_sender: Sender<AlignResult>,
    aligners: &Vec<NoMemLeakAligner>,
    target_idx: &HashMap<String, (usize, usize)>,
    oup_params: &OupParams,
) {
    for query_record in query_record_recv {
        let hits = align_single_query_to_targets(&query_record, aligners);
        let records = hits2records(&hits, &query_record, target_idx);

        if oup_params.discard_multi_align_reads && records.len() > 1 {
            continue;
        }

        if records.is_empty() {
            continue;
        }

        align_res_sender.send(AlignResult { records }).unwrap();
    }
}

pub fn align_single_query_to_targets(
    query_record: &gskits::ds::ReadInfo,
    aligners: &Vec<NoMemLeakAligner>,
) -> Vec<Mapping> {
    let mut all_hits = vec![];

    for aligner in aligners {
        for hit in aligner
            .map(
                query_record.seq.as_bytes(),
                true,
                true,
                None,
                Some(&[67108864]), // 67108864 eqx
                Some(query_record.name.as_bytes()),
            )
            .unwrap()
        {
            if hit.alignment.is_none() {
                // filter unmapped
                continue;
            }

            all_hits.push(hit);
        }
    }

    // set_primary_alignment_according_2_align_score(&mut all_hits);
    set_primary_supp_alignment_according_2_align_score(&mut all_hits);
    all_hits
    // all_hits
    //     .iter()
    //     .map(|hit| build_bam_record_from_mapping(hit, query_record, target_idx))
    //     .collect()
    // set_primary_alignment(&mut align_records);
}

pub fn hits2records(
    hits: &Vec<Mapping>,
    query_record: &gskits::ds::ReadInfo,
    target_idx: &HashMap<String, (usize, usize)>,
) -> Vec<BamRecord> {
    hits.iter()
        .map(|hit| build_bam_record_from_mapping(hit, query_record, target_idx))
        .collect()
}

pub fn build_bam_record_from_mapping(
    hit: &minimap2::Mapping,
    query_record: &gskits::ds::ReadInfo,
    target_idx: &HashMap<String, (usize, usize)>,
) -> BamRecord {
    // println!("{:?}", hit);

    let mut bam_record = BamRecord::new();

    let mut seq = &query_record.seq;
    let mut is_rev = false;

    let rev_seq = match hit.strand {
        minimap2::Strand::Forward => None,
        minimap2::Strand::Reverse => {
            is_rev = true;
            Some(reverse_complement(seq))
        }
    };

    if is_rev {
        seq = rev_seq.as_ref().unwrap();
        bam_record.set_reverse();
    }

    // hit.query_start, hit.query_end 是相对于原始 query 而言的(即 未 reverse 的 query 而言)
    let aln_info = hit.alignment.as_ref().unwrap();
    let cigar_str = convert_mapping_cigar_to_record_cigar(
        aln_info.cigar.as_ref().unwrap(),
        hit.query_start as usize,
        hit.query_end as usize,
        seq.len(),
        is_rev,
    );

    let qual = if let Some(ref qual_) = query_record.qual {
        if is_rev {
            qual_.iter().copied().rev().collect()
        } else {
            qual_.clone()
        }
    } else {
        vec![255; seq.len()]
    };

    bam_record.set(
        query_record.name.as_bytes(),
        Some(&cigar_str),
        seq.as_bytes(),
        &qual,
    );
    if is_rev {
        bam_record.set_reverse();
    }

    // reference start
    bam_record.set_pos(hit.target_start as i64);
    bam_record.set_mpos(-1);
    // bam_record.set_mpos(mpos);
    // bam_record.reference_end()
    bam_record.set_mapq(hit.mapq as u8);

    bam_record.set_tid(
        target_idx
            .get(hit.target_name.as_ref().unwrap().as_str())
            .unwrap()
            .0 as i32,
    );
    bam_record.set_mtid(-1);

    if hit.is_primary {
        bam_record.unset_secondary();
        bam_record.unset_supplementary();
    } else {
        bam_record.set_secondary();
        bam_record.set_supplementary();
    }

    if hit.is_supplementary {
        bam_record.unset_secondary();
        bam_record.set_supplementary();
    } else {
        bam_record.unset_supplementary();
    }
    bam_record.unset_unmapped();

    if let Some(cs) = aln_info.cs.as_ref() {
        bam_record
            .push_aux(b"cs", rust_htslib::bam::record::Aux::String(cs))
            .unwrap();
    }

    if let Some(md) = aln_info.md.as_ref() {
        bam_record
            .push_aux(b"md", rust_htslib::bam::record::Aux::String(md))
            .unwrap();
    }

    if let Some(np_) = query_record.np {
        bam_record.push_aux(b"np", Aux::U16(np_ as u16)).unwrap();
    }

    if let Some(ch_) = query_record.ch {
        bam_record.push_aux(b"ch", Aux::U32(ch_ as u32)).unwrap();
    }

    if let Some(rq_) = query_record.rq {
        bam_record.push_aux(b"rq", Aux::Float(rq_)).unwrap();
    }

    if let Some(dw_) = &query_record.dw {
        if is_rev {
            let dw_ = dw_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"dw", Aux::ArrayU8(AuxArray::from(&dw_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"dw", Aux::ArrayU8(AuxArray::from(dw_)))
                .unwrap();
        }
    }

    if let Some(ar_) = &query_record.ar {
        if is_rev {
            let ar_ = ar_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"ar", Aux::ArrayU8(AuxArray::from(&ar_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"ar", Aux::ArrayU8(AuxArray::from(ar_)))
                .unwrap();
        }
    }

    if let Some(cr_) = &query_record.cr {
        if is_rev {
            let cr_ = cr_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"cr", Aux::ArrayU8(AuxArray::from(&cr_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"cr", Aux::ArrayU8(AuxArray::from(cr_)))
                .unwrap();
        }
    }

    if let Some(be_) = &query_record.be {
        bam_record
            .push_aux(b"be", Aux::ArrayU32(AuxArray::from(be_)))
            .unwrap();
    }

    bam_record
}

/// if reverse alignment, the cigar will be reversed!
pub fn convert_mapping_cigar_to_record_cigar(
    mapping_cigar: &[(u32, u8)],
    mut query_start: usize,
    mut query_end: usize,
    query_len: usize,
    is_rev: bool,
) -> CigarString {
    let mut cigar_str = CigarString(vec![]);

    if is_rev {
        (query_start, query_end) = (query_len - query_end, query_len - query_start);
    }

    if query_start > 0 {
        cigar_str.push(Cigar::SoftClip(query_start as u32));
    }

    mapping_cigar.iter().for_each(|(cnt, op)| {
        let cnt = *cnt;
        let cur_cigar = match *op {
            0 => Cigar::Match(cnt),
            1 => Cigar::Ins(cnt),
            2 => Cigar::Del(cnt),
            3 => Cigar::RefSkip(cnt),
            4 => Cigar::SoftClip(cnt),
            5 => Cigar::HardClip(cnt),
            6 => Cigar::Pad(cnt),
            7 => Cigar::Equal(cnt),
            8 => Cigar::Diff(cnt),
            v => panic!("invalid cigar op :{}", v),
        };
        cigar_str.push(cur_cigar);
    });

    if query_end != query_len {
        cigar_str.push(Cigar::SoftClip((query_len - query_end) as u32));
    }
    cigar_str
}

/// {"target_name": (idx, length)}
pub fn targets_to_targetsidx(targets: &Vec<ReadInfo>) -> HashMap<String, (usize, usize)> {
    let mut target2idx = HashMap::new();
    targets.iter().enumerate().for_each(|(idx, target)| {
        target2idx.insert(target.name.clone(), (idx, target.seq.len()));
    });
    target2idx
}

pub fn set_primary_alignment_according_2_align_score(hits: &mut Vec<Mapping>) {
    if hits.is_empty() {
        return;
    }
    hits.sort_by_key(|v| {
        -v.alignment.as_ref().unwrap().alignment_score.unwrap()
            - if v.is_primary { 1_000_000_000 } else { 0 }
    });

    assert!(hits.first_mut().unwrap().is_primary); // assertion failed
    assert!(!hits.first_mut().unwrap().is_supplementary);

    let primary_hit = hits.first().unwrap();
    let (primary_qstart, primary_qend) = (primary_hit.query_start, primary_hit.query_end);

    hits.iter_mut().skip(1).for_each(|hit| {
        if hit.is_primary {
            let ratio = ovlp_ratio(primary_qstart, primary_qend, hit.query_start, hit.query_end);
            hit.is_primary = false;
            hit.is_supplementary = ratio <= 0.5;
        }
    });
}

/// 重新设置 primary 和 supplementary。其主要在 多个参考基因序列场景使用
/// 由于对于多个参考序列场景，每个参考序列都构建了一个aligner，所以我们会得到之多 N 个 primary。N表示参考序列个数。
/// primary设置逻辑：primary中，比对分数最高的作为最终的primary, 其余的primary基于 query 的 ovlp 来确定自己是 supplementary 还是 secondary
/// supplementary设置逻辑：所有hits重新排列，按照 primary 在第一位，supplementary 按照 比对得分位列2~？。secondary排到最后
///     依次判断 当前 range 和 之前的 supplementary range的 ovlp。如果 ovlp>0.5 当前 supp 就降级为 secondary
pub fn set_primary_supp_alignment_according_2_align_score(hits: &mut Vec<Mapping>) {
    if hits.is_empty() {
        return;
    }
    hits.sort_by_key(|v| {
        -v.alignment.as_ref().unwrap().alignment_score.unwrap()
            - if v.is_primary { 1_000_000_000 } else { 0 }
    });

    assert!(hits.first_mut().unwrap().is_primary); // assertion failed
    assert!(!hits.first_mut().unwrap().is_supplementary);

    let primary_hit = hits.first().unwrap();
    let (primary_qstart, primary_qend) = (primary_hit.query_start, primary_hit.query_end);

    hits.iter_mut().skip(1).for_each(|hit| {
        if hit.is_primary {
            let ratio = ovlp_ratio(primary_qstart, primary_qend, hit.query_start, hit.query_end);
            hit.is_primary = false;
            hit.is_supplementary = ratio <= 0.5;
        }
    });

    hits.sort_by_key(|v| {
        -v.alignment.as_ref().unwrap().alignment_score.unwrap()
            - if v.is_primary {
                1_000_000_000
            } else if v.is_supplementary {
                100_000_000
            } else {
                0
            }
    });

    let mut visited_ranges = vec![(primary_qstart, primary_qend)];
    hits.iter_mut().skip(1).for_each(|v| {
        if v.is_supplementary {
            let (start, end) = (v.query_start, v.query_end);

            let cnt = visited_ranges
                .iter()
                .map(|&(range_s, range_e)| ovlp_ratio(range_s, range_e, start, end))
                .filter(|&ovlp| ovlp > 0.5)
                .count();
            if cnt > 0 {
                v.is_supplementary = false;
            } else {
                visited_ranges.push((start, end));
            }
        }
    });
}

pub fn mapping2str(hit: &Mapping) -> String {
    format!(
        "qstart:{}, qend:{}, primary:{}, supp:{}, identity:{}, score:{:?}",
        hit.query_start,
        hit.query_end,
        hit.is_primary,
        hit.is_supplementary,
        MappingExt(hit).identity(),
        hit.alignment.as_ref().unwrap().alignment_score
    )
} 

#[allow(unused)]
fn get_query_start_end(hit: &Mapping, query_len: i32) -> (i32, i32) {
    match hit.strand {
        minimap2::Strand::Forward => (hit.query_start, hit.query_end),
        minimap2::Strand::Reverse => (query_len - hit.query_end, query_len - hit.query_start),
    }
}

fn ovlp_ratio(s1: i32, e1: i32, s2: i32, e2: i32) -> f32 {
    if s2 >= e1 || s1 >= e2 {
        return 0.0;
    }

    let s = cmp::max(s1, s2);
    let e = cmp::min(e1, e2);

    let ovlp_len = e - s;

    let min_len = cmp::min(e1 - s1, e2 - s2);

    if min_len == 0 {
        return 0.0;
    }

    ovlp_len as f32 / min_len as f32
}

#[cfg(test)]
mod tests {
    use std::time::{self, Duration};

    use gskits::fastx_reader::read_fastx;
    use minimap2::ffi::mm_idx_destroy;
    use rust_htslib::bam::ext::BamRecordExtensions;

    use super::*;

    #[test]
    fn test_align_single_query_to_target() {
        let ref_file = "test_data/GCF_000005845.2_ASM584v2_genomic.fna";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);
        let aligners = build_aligner(
            "map-ont",
            &IndexParams::default(),
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
        );
        let target2idx = targets_to_targetsidx(&targets);

        let query_file = "test_data/ecoli_query.fa";
        let query_filter_iter =
            gskits::fastx_reader::fasta_reader::FastaFileReader::new(query_file.to_string());
        for qr in query_filter_iter {
            let hits = align_single_query_to_targets(&qr, &aligners);
            let records = hits2records(&hits, &qr, &target2idx);
            for record in records {
                assert_eq!(record.reference_start(), 720);
                assert_eq!(record.reference_end(), 1920)
            }
        }
    }

    #[test]
    fn test_set_primary_alignment() {
        let mut r1 = BamRecord::new();
        r1.unset_secondary();
        let mut cigar_str = CigarString(vec![]);
        cigar_str.push(Cigar::Equal(3));
        cigar_str.push(Cigar::Diff(1));
        r1.set(b"1", Some(&cigar_str), b"AACG", &vec![255; 4]);

        let mut r2 = BamRecord::new();
        r2.unset_secondary();
        let mut cigar_str = CigarString(vec![]);
        cigar_str.push(Cigar::Equal(4));
        r2.set(b"2", Some(&cigar_str), b"AACT", &vec![255; 4]);

        let mut r3 = BamRecord::new();
        r3.set_secondary();
        let mut cigar_str = CigarString(vec![]);
        cigar_str.push(Cigar::Equal(2));
        cigar_str.push(Cigar::Diff(2));
        r3.set(b"3", Some(&cigar_str), b"AAGC", &vec![255; 4]);

        let mut records = vec![r1, r2, r3];

        // set_primary_alignment(&mut records);

        // assert_eq!(records[0].is_secondary(), true);
        // assert_eq!(records[1].is_secondary(), false);
        // assert_eq!(records[2].is_secondary(), true);
    }

    #[test]
    fn build_aligner_memory_leak() {
        for _ in 0..10000 {
            thread::sleep(Duration::from_millis(1));
            let aligner = Aligner::builder().map_ont();
            let aligner = aligner
            .with_index_threads(1)
            .with_cigar()
            .with_sam_out()
            .with_sam_hit_only()
            .with_seq_and_id(b"ACGGTAGAGAGGAAGAAGAAGGAATAGCGGACTTGTGTATTTTATCGTCATTCGTGGTTATCATATAGTTTATTGATTTGAAGACTACGTAAGTAATTTGAGGACTGATTAAAATTTTCTTTTTTAGCTTAGAGTCAATTAAAGAGGGCAAAATTTTCTCAAAAGACCATGGTGCATATGACGATAGCTTTAGTAGTATGGATTGGGCTCTTCTTTCATGGATGTTATTCAGAAGGAGTGATATATCGAGGTGTTTGAAACACCAGCGACACCAGAAGGCTGTGGATGTTAAATCGTAGAACCTATAGACGAGTTCTAAAATATACTTTGGGGTTTTCAGCGATGCAAAA",  b"ref")
            .unwrap();
            let aligner = NoMemLeakAligner(aligner);
        }
    }

    #[test]
    fn test_align_single_query_to_target2() {
        let ref_file = "/data/ccs_data/MG1655.fa";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
        );

        let aligner = &mut aligners[0];

        // aligner.mapopt.min_cnt = 2;
        aligner.mapopt.best_n = 100000;
        // aligner.mapopt.min_dp_max = 10; // min dp score
        // aligner.mapopt.pri_ratio = 0.5; // min dp score
        aligner.idxopt;
        println!("{:?}", aligner.mapopt);
        println!("{:?}", aligner.idxopt);

        // tot len
        let seq = b"AAAAGAGAGAGATAGGGGGCGACGAGCGCCAAACGGGCGGCAATTATAGGGATTTCATCGCCTGATACCAGTCGAATAGCGTTGCCGCGCGCTCAGGATGTTAATTGGTTGACAGAAGAATTCCCGGTGGCAAAATTACGTTGAAGATCAGTTTTATACAAGGTAAAAAATGTTATACGCAGTTGCGCAATTATCGCCTTTACGTCACTTTATGAGCATTCGCATATAAAATGTAAAACTTTTTGTAACTAGCATAAAACACAGAAACGAATACCTGGCGCTGGTCTTGCGATAAAGCAGGTAATGAGCAAACAACACAGCATGTATTAATTGCCCTGCCCACCCGCTGCTTCCACCTGGTCCAGTTTAAGGTTAGTCCTCGTTTACTTTACCCTTTTCTCGCTGAGCTTTCGCAAGTTTGGCACCCTCTCGCCCCACGCCTGTGGTTCCGACGGTCCACTGATGGTGGCGTTTATCGCCATGCCGGGGCGCAATGTGCCGGGAAATGCGCTGAGCTGTTCGCTGGAAATATCGCCGCATCCATCCTGCTTTTTTCCACCAGCCGCTGAACATGACCTGGACGACATCAATATTGTTGAAGCCGTGGTCGGGCAGTGCATCTCTCGCTCAAACAACAACGGAGGAGGAGAAAAAAGAGACAGATGCACTGCCCCGACCACGGCTTCAACACATATTGATGGTCGTCCAGGTCATTTCCACGAGTGGTGGAACAAAGCAGGAGGATGCGGCGTATTTCCAGGAACAGCTCACGCAATTCCCGGCACATGGCGCCCGGCAATGGCGAATAAAACGCCAACCATCATGAATGGACGTCGGAACCACATGGGCGATTGGTGCCAAACTGCGAAAGCTCAAGCGAGAAAAGGGTAAAGAAAAACGAGACTTATACCTAAATGACCAGGGAATGCAGCGGGTGGGGCAGGGCAATTAATTACATGGCTGTGATTGTTGCTCATACCGCTTTATCCCAAGACAGGTCGCCAGTATTCGTTTCTGTGTTTAATGCTATACAAAAGTTTATTTACATTTATATGCGAAATTGCTCATAAGTGAGGGTAAAGGCGATAATTTGCGCAACTGCGTAATAACATTTTTTTACCTTGTATAAAACTGATCAACGTAAAGTTGCCACCGGGATCTTCGTCAAAATTAACTCTGACGCGCGGCAACGCTTATTCGACTGGTATCAGGCGGATGAAATCCTATAAATTGCCGCGTTTGGCGCTTCCGTCGCCCCCTATCTCTCTCAAACAACAACGGAGAGGAGGAAAAAAGAGAGAAGATAGGGGGCACAGAAGCGCCCAACGCGCAATTAATAGGGATTCATCCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCATGATACCAGTCGAATAGCGTTGCGCGCGCTCAGAGTTACATTGTTGACGAAGAATTCCCGGGTGGCAAATTACGTTGATCAGTTTTATAACAAGGTAAAAAGTTATACGCAGTTGCGCAAATTATCCGCCTTACGTCACTTTAATGAGCAATTCCGCATAATAAAAATGTAAAACTTTGTAACTAGCATAAACACAGAAACGAATACTGGGGCGACCTGGTCTTGCGGAATAAGCGGTAATGAGCAAAAATCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGCTTCACCGGTCAGCGTTTAGGTTTAAGTCTCGTTTAATCTTTACCCTTTCTCGTTGAGCGAGCTTTCGAGTTTGCACCCAACTCGTCCCACTGTGGTTCCCGACGTCCATCATGTGGGCGTTTTAAAATCGCGCATGCCGGCGCATGTGGCCGGGAATTGCCTGAGCTGTTACGCTGGGAAATATCGCCGCATCCATCCTGGCTTTTTTTCCACAGCTACGCTGAACAGACCTGGACGACCACCATCAATATGTTGAAGCCGTGGTCGGGCGTGCATCTCTCTCAAACAACAACGGAGGAGGGAGGAAAAAAGAGAGAGATGCACTGCCCCGACCACGGTTTTTTGTCAACATATTGAAGGTCGTCCAGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAATCCCGGCCACAGCGCCCGGCAATGGCGATTAAAACGCCACCATCATGATGGAACGTCGGGAACCACAGTGGGGCGTGTTGGGTGCAAACTCGAAAAGCTCAAGCAGAAAAGGGTAAAGATAACACGAGACAAACCAAACTGACCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGTTTATCCGCAAAACAGGTCGCCAGTATTCGTTTCCTGTGTTTATGCTAGTACAAAAGTTTTACATTTTAATATGCGAATTGCTCATAAAGTGACGTAAGGCGGAAATTTGCGCAACTGCGTAATAACATTTTTTTAACTTGTATAAAACTGTTCAACGTAATTTGCCACGGAATTCTTCGTCAACAATTAACTACCTGAGCGCGCCAACGCTATTCGACTGGTATCAGGCGAATGAATCCCTATAATTGCCGCGTTTGGCGCTTCGTCCGCCCCCTATCTCTCTCAAACAACAACAGGAGGGGAGGAAAAAAGAGAGAGATAGGGTGGCGACGAAGCGCCAAACGCGGCAATTATAGGGATTTCATCCGCCTGATACCAGTCGAATAGCGTTGCCGCGCGCTCAGAGTTAATTGTTTGACGAAGAATTCCCGGTGGCAAATTACGTTGATCAGTTTTTATACAAGGTAAAAAAATGTTATACGCAGTTGCGCAAATTATCCGCTTTACAGTCCACTTTTAGAGCAATTCGCATAATAAAAGGTAAAACTTTTGTACTAGCATAAGAACACACACAGAAACAATACTGGCTGACCTGGTCTTGCGGATAAAGCGGTAATGAGCAAACATCACGCATGTATTAATTTGCCCCTCCCCACCGCTGCTTCAACCGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTAGTACCCTTTTCTCGCTTGAGCTTTGCAGTTTGGCACCCAACTCCCATGTGTTTCCGCGACGTCCATCATGATGGTGGCGTTTTATCGCATGCCGGGCGCATGTGGCGGGAATTGCGCTGAGCTGTTCGCTGGGACATATCGCCGGCATCCTCCCTGCTTTTTTCCCCAGCTCGCGCTGAACATGACCTGGACGACCAAGTCAATATTGTTGAGCCGTGGTCGGGTGCAAGTGCATCTCCTCAACAACAACGAGGAGGGGAAAAAAAAAGATGAGAGATGCACTGCCCCCGACCACGGCTTCACACAATATTGAGGGTCAGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGATGGAGCGCGATATTTCCCACGACAGCTCAGCGAAACAAAAAAAAAAAATTCCGCCACATGCGCCCGGCATGGCGATAAAACGCCACCTCCAGATGGACGCGGGAAACCACAGTGGGGCGAGTTGGGTGCCAAACTGCGAAAGTCAAGCAGAAAAGGGTAAAGAAAACGAGCTAAACTAACTGCACAGGTGAAGCAGCGGGTGGGCAGGCATAATACATGTGTGATTGTTTGGTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCATTACCGCTTTATCCCAAGACCAAGGTCGCCAGTATTCGTTTCTGGGTTTATGCTAGTACAAAAGTTTACATTTTATATGCGAATTGCTCCATAAAGTGGCCCCCCCCCCCCGTAAAGGCGGATAATTTGCGCAACTGCGTATAACATTTTTTACCTTGTATAAAACTGATCAAGTAATTTGCCACCGGAATTCTTTCGTCAACAATAACTCTGGCGCGCGGCAACGCTATTCGACTGGTATCAGGCGGAGAAATCCTAATAATTGCCGCGTTTGGCGCTTCGTCGCCCCTATCTCTCTCCAAACAAACAACGAGGAGGAGGAAAAAGAGAGAGATAGGGGCGACGAAGCGCCAAACGCGGCAAATTATAGGGATTCATCCGCCTGATAACCAGGTCGAATAGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCAGCGTTGCCCGCGCTCAATTTAATGTTGACGAAGAATTCCCGGTGGCAAATTAACACGTTGAATCAGTTTTATACAAGGTAAAAAATTTATAGCAGTTGCGCAATTAATCCGCCTTACGTCACTTATGAGCAATTCGCATATAAAATGCTAAAACTTTTTGTACTAGCATAAAACACAGAAACGAATACTGGCGACCTGGTCTCGCGGATAAAGCGGTAGATGAGCAAAACAATCACAGCACTGTATTAATGCCTGCCCCACCGCTGCTTCACCGAGTCAGTTTAGGTTTATAGTCTCGGTTATCTTTACCCTTTTCGCTTGAGCTTTCGCAGTTTGGCACCACACTGCCCATGTGGTTCCGCGTCCATCCATGATGGTGGCGTTTTTATCGCCATGCCGGGAGCGCATGTGGCCGGAATTGCGCTGAGCTGTTCGCTGGAAATATCGCCGCATCCGATCCTGCTTTTTTCCACCAGCTCGCTGAACATGACCGGACGACCATCAATATTGTGAAGCCGTGGTCGGGGCAGTGCATCTCTCTCAACAAACAAACAACAGGAGGAGAGAAAAAAGAGAAAGCATGCCCCGACCACGCTTCAACAATATTGAGTGGTCGTCCGGTCATGGTTCAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGCGATATTTCCCAGCGGAACACTCAGCGCAATTCCGGCCAAACAGCGCCCGGCATGGCTGATAAAACGCCACCATCATGATGGACGTGGGAACCACAAGTGGGGCGATTGGTGCCAACTGCGAAAGCTCAAGCGAGAAAGGGTAAAAGATAAACGATGACCTAAACCTAAACTGACCAAGTGAACAGCGGGTGGGGCAGGCAATTAATACATGTGTGATTGTTTGCCATTACCGCTTTATCCCAAGACCAGGTCGCCAGTATTCGATTTCTGTGTTTATGCTAGTACTAAAGTTTTACATTTATATGCAGAATTGCTCATAAAGTGAGTAAGGTCGGATAATTTGCGCACACTGCGTATAACATTTTTTACCTTGTATAAAACGATCCAACGTAATTTGCCACCGGGAATCTCGTCAACAATTAACTCTGAGCGCGCGGCAACGCTATTCGACTGTATCAGGCGGATGAAATCCCTATAATTGGCCGCGTTTGGCGCTTCGTCGCCCCCATCTCTCTCAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATAGGGGGCGACGAAGCGCCAAACGCGGCATATAGGGATTTCATCCGCCTGATACCAGTCCGAAAGCGTTGCCGCGCGCTCAGAGTAAGTTGACGAAGAATTCCCGTGTGGGCACAATAACGTGACAGTTTTATACAAGGTAAAAAATGTTATACAGCAGTTGCGCAAATTTCCCCTTTACGTCACTCACTTTATGAGCAATTCGCATATAAAATGTAAACCCTTTTGTACTAAGCATAAACAGCAGAAACGAATACCTGCGATACCTGGTCTTGCGGATAAACGGTAATGAGCACAACTCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGTCTTCAACCTGGTCATGTTTAGGTTTATCTCGTTTATCTTACCTTTCTCGCTTGAGCTTTCGCATATTGGCACCCAATCGCCCCACGTGGTTCCCGACGTCCATCATGAGGTGGGCGTTTTATCCCATGCCGGGCGGCAATGTGGCCGGAATTGCGCTGAGCTGTTCCGCTGGGAAAATCGCCGCATCCATCTAGCTTTTTTTCCACAGGCTCGCTGAACATGATACCTGGACGACCATCAATATTGTTGAAGCGGGTCGGCAGTGCATCTTCTCAAACAACAACCGTAGGAGGAGGAAAAAAGAGAGAGATGCACTGCCCCGACCAACGGCTTCACAATATTGATGGTCGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGTCGGCGATATTTCCCAGCGAACAGCTCAGCGCAATTCCCGGCCACATGCGCCCGGCATGCGATATAAAGCCACATCATGATGGACGTCGGAACCACATGGAGCGAGTTGGGTGCCAAACTGTGGAAAGCTCAGCGACAAAAGGGTAAAGATAAACGAGACTAAACCTAAACTGACCAGTGAAGCAGCGGGTGGGCAGGGCAATTAATACATGCTGTGATTGGTTTGCTCATACCGCTTACGCAAGACCAGGCCGCAGTATTCTTTCTGTGGTTTATGCTAGTACAAAAGTTTTACATTTTATATGCGAATTGCTCAAAAGTGACGTAAAGCGGATAATTTGCGCAACTCGTATAACATTTTTACCTTGGTATAAAACTGATCAACGTAATTGCCACCGGGAATCTTCGTCACAATTAAGACCTGAGCGCGCGGCAACGTATTCGACTGGTATCCAGCGGATGCAAATCCCTATAATTGCCGCGTTGGCGCTTCGTCGCCCCATTCTCTCAAACAACAACGGAAGGAGAGGAACAAAGAGAGAGCATAGGGGGGCGGGGGCGAACGAAGGCAAACAGCGGCAATTATAGAGGATTTCACTCCGCCTATACCAGTCGAATAGCGTTGCCGCGCGCCAGAGTTAATTGTGACGAAGAATTCCCGTGGCAAATTAACGATGATCAGTTTTATACAGGGTAAAAAATGTTATTACAGCAGTTGCGCAAATTATCCGCCTTACGTCACTTTATGCAGCAATTCGCATATAAAATGTAAAACATTTTGTACTAGCATAAACACAGAAACGAATACTGGCGACCTGGTCTTGCTCGGATAAAGCGGTTAATGAGCAAACAATCACCAGCAGTATTAATTGCCCTGCCCCCACCGCTGCCTTCACCTGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTCTCGCTTTGAGCTTCGCAGTTTGGCACCCAACTCGCCCCCACTGTGGTTCCCGACTCCCATCATGATGGTGGCGTTTTATCGCCATGCGGGCGCATGTTGGCCGGGAATGCGCTGACTGTTCGCTGGGAATATCGCCGCATCATCCTGCTTTTTTCACCAGCTCGCTGAACATGACCTGGACGACATAATAGTTGAAGCGTGTCGGCAGTGCATCTCTCTCAAAACAACAACGAGGAGGAGGAAAAAAGAGAAGATGCCTGCCCCCGCACAGGTTAACAATATGAGGTCGTCCAGGTCGATGTTCAGCGGACTGGTGGGAAAAAAAAAAAAAGCAGGATGGAATGCGCATATTTCCCAACAGCTCAGCGCAATTCCCGGCCCACAATGCGCCCGGCATGGCGATAAAACGCCACCACAGAGGACGTCGGAACCACAGTGGGGCGAGTTGGGTGCCAAATCTCGAAAGCCAGCGAGAAAAGGGTAAAGATAAAACGAGACTAAACCTAAACTGACCAGGTGAAGCAGCCGGGTGGGGCAGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCCGGCTTTATCCCAGACCAGGTCGCCAGTATTCGTTTCTGTGTTTATGCTAGTACAAAAGTTTAACATTTTATAGCGATTGTTCATAAAGTGACGTAAAGGCGGATAATTTGCCAACTGCGTATAACGATTTTTACCTTGTATAAAACTGATCAACGTATTTGCACGGGAATTCTTCGTCCCCAATTGAACTCTGAGCGCGCGCACGCTATTCGACCTGGTATCAGGCGGATGAAATCCGTATAATTGCCCGTTTGCGCTTCGTCGCCCCTATCTCTCTCAAACAACAACGGAAGGAGAGGAAAAAAAGAAGAGATAGGGGGCGACGAAGCGCAAACGCGCAATTTATAGGGATTTCTCCAGCCTGATACCAGTCGACATAGCGTTGCCGCGCGCTCAGAGTTAATGTTGACGAAGAATTCCCGGGTGCAAATTACGTTGATCAGTTTTTTTTTATACAAGTAAAAAATGTTATACGCAGTGGCAAATTTCCGCCATTTACGTCACTTTATGAGCAATTCGCATAAAAAATGTAAAGTTGCATTGTAAAATCTTTTGTACAGCATAAACACAGAAACGAATAGCTGGCGAACCTGGTCTGCGGATAAGAGCGGTAATGAGCAACACACAAGCATGGTATAATTGCCCTGCCCCACCCGCTGCTCACCTGGTCGTTTAGGTTAGTTCGTTTACTTTAACCCTTTCTCGCTTGAGCTTCGCAGTTTTGGCACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGGGCGTTTTACGCCATGCCCGGGCGAGTGGCCGGGAATTGCGCTGAGCTGTTCGCTGGGAATATCGCCGCACCATCCTGCTTTTTTCCACCAGGCTCCGCTGGAACATGCCTGGCGACCAATCATATTGTTGAAGCCGTGGTCGGGGCAAGTGCATCTCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGCATGCCCGACCACGCTTAACAATATTGATGGTCGTCCAGGTCATGTTCAGCGGAGCTGGTGGAAAAAAGCAGGATTGGATGCGGCGATATTTCCCAGCGAACAGCTCAGCGCAATTCCGCCACATGCGCCGGCATGGCGATAAAACGCACCATCAATGTGGACGTCGCTGAACCGAAGTGGGGCGAGTTGGTGCCAAACTGCGGAAAGCTCAAGGAGAAAAGGGTAAAGATAAACGAACAAACTAAATGACCAGGTGAAGCACGGTGGGGCCAGGGAATTAATACAGCTGTGATTGTTTGTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTCATTACCGCTTTATCCGAAGACCAGGTCGCCATATTCGTTTCTGTGGTTTAGCTAGTACAAAGTTTTACATTTTATATGCGAATTGCTCATAAAGTGACGTAAAGGCGAATAATTTGCGCACTGCAGTATAACATTTTTTACTGTATAAAACTGATCAACGTAATTTGCCACCGGGATTCTTCGTCAACAATTTACTCTGAGCGCGCGGCACGCTATTCGACTGGTAATCAGGCGGATGAAATCCCTATAATTGTCCCCCCCGGGGGGCGTTGGCGCTCGCGCCCCCTATCCTTCTCAACAACAACGGAGAGAGGAAAAAGAGAGAGATAGGGGGGGGCGACAAGCGCCAAACGCGGCAATTATAGGGATTTCCATCCGACCTGATACCATCGAATAGCGTTTGCCGCGCGCTCAGAGTGTAATTGTTGACAGAATTCCCGGTGGCAAATTACGTTGATCAGTTTATACAAAGGTAAAAATGTTATACGCAGTTGCGCAAATTATCCGCCTTTAGTCCTTATGAGCAATTCGCATAAAACCTAAAATGTAAAACTTTAGTACTTGCAATAAACACAGAAACAATACTTGGCGACCTGGTCCCTCTGCGGATAAAGCGGTAAATGAGCAAACAACACAGCATGTATTAATTGCCCTGCCCCAACCCGCTGCTTCACTGGTCAGTTTAGGTTTAGTCTGTTTTATCTTTACCCACTTTTCTCGCTGAGCTTCCGCAGTTTGGCACCCAACTCGCCCACTGTGTTCCCGACGTCCATCATGATGGTGGCGTTTATCGCCATGCCGGGCGCATGTGGCCGGAATTTGCGCTGAGCTGTTTGCTGGGGAAATATCGCTGCATCCATCCTGCTTTTTTCCACGAAAAAAACAGCTCGCAGAAATGACCTGGACGACCATCATATTAGTTGAGCCGTGGTCGGGGCATGCATCTCTCTCAAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGCACTGCCCCGAACACGGCTTCAACAATATTGATGGTCGTCCAGGTCCCCATGTTCAAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGTGCGATATTTCCCATCCCCCCCCCGAGCGAACAGCTCAGGCAATTCCCGGTCCACATGCGCCCGGCATGGCGGATAAAACGCCACCCGATCATGATGGACGTTCCAGGGAACCACATGGGGCGAGTTGGTGCCAAACTGCGAAAGCTCAAGGCGAGAAAAAGGGTAAAGATAAACGGACAAACCTAAACTGTCCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGGCTGTATTGTTTGCCATTAGCCCCCCCCCCCCCCCCCATTCCGCTTTTATCCCGCAAGACCAGGTCGCCATATTCGTTCTGTGTTTAATCTAGTACACAAGTTTACATTTATAGCTAATGCGAATTGCTCATAAGTGACGTAAAGGCGGATAAATTCGCAACTGCGTATAACATTTTTTACTTGTATAAAATGATCAACGTAATTTGCCACCGGGAATCTTCGCAACAATAACTCTGAGCGCGCGGCACGCTATTCGACTGGTATCAGGCGGTGAAATCCCTAAATTGCGCGTTTGGCGCTTCGTCGGCCTTCCTCTCAAACAACACGAGGAGGAGGAAAAAAAAAAAAGAGAGAGATAGGGGGCGACGAAGCGCAAACGCGCAATTATAGGGATTTATCCGCCTGATACCAGTCGATAGCAGTTGCCGCGTCGCTCAGAGTTAATTGTTGACAAGAATTCCGCGGTGGCAAATTACGTGATCAGTTTTATACAAGAGTAAAAAAGTTTATACGCAGTTGCGCAAATTATCCGCCTTTACAGTCACTTATGAGCAATTCGCATATAAATGTAAAACTTTTGTACTAGCATAAACACAGAACACGAATACTGGCGACCTGGTCTTGCGGAAAAGCGGTAATGAGCAAACAAATCACAGCATGTATTATTGCCCTGCCCCACCCGCTGCTTCAAACCTGGTCAGTTTAGGTTTAGTCTCGTTTATCTTACCTTTTGCTCGTTGAGCTTTCCAGTTTGGCACCCAACTCGCCCACTGTGGTTCCCGACGTCCATCATGATGGTGGCTGTTTTATCGCCCATGCCGGGCGCATGTGGCCGGAATTGCGCTGAGCTGTTCGCTGGGAAACTATCGCCGCATCCATCCGCTTTTTTCCACCATGCTCGCTGAACATGACCTGGACGCACCATTCATATTGTTGAAGCCGTGGTGGGGCATGCATCTCTCTCAAACAACAACGGAGGAAGAGGAGAAAAAGGAGAGAGATGCACTGCCCGACCACGGCTTCAACAATATTGAGGTCGTCAGGTCATGTTCAGCGAGCTGGTGGAAAAAAGCAGGTGGCAGGATGGGCGGCCGATATTGCCGAGCGAACAGCTCAGCGCAAATTCCGGCACATGCGCCCGGCATGGCGAATAAAACGTGCACCATCATGATGGACAGCGGGAACCACAGTGGGGCGAGTTTTGGGTGGCCAACTGCGAAAGCTCAACGAGAAAAGGGTAAAGATAAACGAGACTAAACTAAACTGACAGGTGGAAGCAGCGGTGGGCAGGCAATTAATACAGCTGTGATTGTTTGCTCATTACCGCTTTTCCGCAAGACCAGGTCGCAAGTATTCGTTTCTGTGTTTATGCTAGTAACAAAAGTTTTACTTTTATATGCGAAAATTGCTCAAAAGTGACAGTAAAAAGGCGGATAATTTTGCGCACATGTCGTATAACATTTTTACCTTGTATAAAACTGACAACGTAAATTTGCCACCAGGGAATTCTTCGTCCAACAATTTAACTCTGAGCGCGCGGCAACGCTATTCGACTGTATCAGGCGGATGAAATCCCTATAATTGCCGGCTGATTGGCTGCTTCGTCCGCCCCCTACGTCCTCAAACAAAACGGAGGAGGAGGAAAAAAAGAGAGATAAGGGGCGACAAGCGCAAAACCGGCAAATTAAGGGATTTTCATCCGCCTGATACCAGTCAGGAATAGCGTTGCCGAGCGCTCAGAGTTAATTGTGACAAGAATTCCGGTGGCACAAATTATTGTCAGTTTATACAAGGTAAAAAATTTATACGCAGTTGCGCAAATTATCCGCCTATTACGTCACCTTTAGACAATTCGCATAAAAATGTAAAACTTTTGTACTAAGCATAAACACAGAAACGAATACTGGCGACTTGGTCTTGCGGATAAAGCGTAATGAGCAAACAATCACACATGTATTAATTGCCCTGCCCCAGCCGCTGCTTCACCAGGTCCAGTTTAGGTTTAGTCTCGTTATCTTTACCTTTTCTGCTTGAGCTTTCGCAGTTTGGCACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGTGCGTTTTATCGCATGCCGGGCGCATGGCCTGGAATTGCGCTGAGCTGTTCAGCTGGGAATATCGCCGCATCCACCTGCTTTTTTCCCACAGCTCGCTGAACAGACCTGGAACGACCATCATATTGTTAGAAAGCCGTGGTCGGGGCCGTGCATCTCTCTCAAACCAACAACGGAGGACAGAAAAAAAGAGAATGCACTGCCCGCCACAGGCTTCAACAATATTGATGGTCGTCCAGGTCATGTTTCAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGCGATATTTCCCCAGCGAACAGCTCAGCGCAATTCCCGGCCAATGCGCCCGGCATGGCGATAAAACGGCCAACCATCATGATGGACGTCGGGACCCATGGGGGCGAGTTGGGTGCCAAACTGCAAGCTCAAGCAGAAGGGTAAAGATAAGATAAACGAGACTAAACCTAAACTTGACCAATGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTAGCTTTAATCCGCAAGACCAGGTCGCCAGTATTCGTTTCTGTGTTTATGCTAGTACAAAAGTTTTATGACATTTTTATATGCGAATTGCTCAAAAGTGACGTAAAGGGGATAAATTTGCGCAACTGGCGTATAACATTTTTTACTTGGTATAAAGCCTAAAACTGATCAACCGTAATTTGCCACCGGGAATTCTTCGTCAACAATTAAACTCTGAGCGCGCGGCAACGCTATTCCGACTGTAATCAGCGATGAAATCCCTATAATTGCCCGCGTTTGCCGCTTCTCGCCCTATCCTCTCAAACAACAACGGAGGAGGGGAAAAAAGAGAGAGATAGGGGCGACGAAGGCCGCCAAACGCGGCAATTATAGGGATTTTCATCCGCCTGATACCAGTCGAATAGCGTTTGCCCGCGCGCTCAGATTAATTGTTGACGAAGAATTCCCGGTTGGCAAATTACTCTGATGCATGCTTTTAGACAAGGTAAAAAGATGTTATACGCAGTTGCGCAAATATCCGCCTTTACGTCACTTTATGAGCAATTCGCATCATAAAATGTAAAACTTTGTACTAGCATAAACAAGAAACGAATACTGGCGAACCTGGTCTGCGGATAAAGCGGTAATGACAAACAATCACAGCATGTATTAATTGCCCTGCCCCACCGCTGCTTCACCGGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTTCTCGCTTGAGCTATTCGCAGTTTGGCACCCTTCTCGCCCACAGAGAAGAGGTGGTTGCCCGACGTCCCATCATGATGGTGGCTTTTATCGCCATGCCGGGCGCATGGGCCGAATAGCGCTGAGCTGTTCGCTGGGGAATATCGCCCGCATCCAATCCTGCTTTTTTCCACGAGCTCGCTGAACATGACCTGGACGACCACAATATTGTTGAAGCCGTGGGTCGGGGGCCAGTGCGATCTCTCTCAACAAAAACGGAGGAGGAGGAGAAAAAAGAGAGGAGCACTGCCCCGACCACGGCTTCACAATATTTGATGGTCGTCCAGGTCATGTCAGCAGAGCTGGTGGAAAAAAGCGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAAATTCCCGGCCAACATGCGCCCGGCAATGCGATAAACCGCCACCAATCATGATGGACGTCGGGGACCACAGTGGGGCGAGTTGGGTGCCAAAACTGCGAAAGCTCAAGCGAGAAAAGGGTAAAGATAAACGAGACAAACCTAAACTGGACCAGGTGAAGCAGCGTGGACAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGCTTTATCCGCAAGACAGGTCGCCAGTAATCGTTTCTTGTGTTTAGCTAGTACAAAAGTTTTTATTTTATATGCGAATTGCTCATAATGTGACTGTAAAAGCGGATAATTTGCGCAACTGCGTATAACATTTTTTACCTTGTATGAAAAACTGATCAACGTAATTTGCCAACCGGGAATCTTCGTCAACAAATAACTCTGAAGCAGCGCGCAACGGCTATTCCGACTGGTATCAGGGCGGCATGAAATCCCTATAATTGCCGCGTTTGGCGCTCGCGCCCCCTATCTCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGCAGATAGGGGCGAGAAGCGCCAAACGCGGCAATTATAGGGATTTCATCCGCCTGATACCAGTGAAATAGCGTTGCGCGCGCTCAGAGTAATTTTGAGAAGAATTCCCGGTGGCAATTACGTGATCAGTTTTATACAAGGTAAAATGTTATACGCAGTTGCCAAATTATCCGCCTTTACGTCACTTTTATGGCAATTCGCATATACAAATGTAAAACTTTTGTACTAAGCATAAACACAGAAACGAATAGCTGGCCGACCTGGTCCTTTGCGGATAAGGGGAAGCGGTAAATGAGCAAACATCACAGCATGTATTAATTGGCCCTGCCCCCCGCTGCTTCACCGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTTCTCTGAGCTTTTCGCAGTTTGGACCCAAATACTCGCCCCACTTGTGGTTCCCGACGTCCATCGATGATGGTGGCGTTTTTATCGCCAATGCCGGGCGCATGTGGCCAGGGAAATTGCGCTGAGCTGTTCGCTGGGAATATCCGCCGCCAATCCTCCTGCTTTTTTTCCACCAGCTCGCTGAACATGACCTCGGACGCCATCAATATTGTTGAAAGCCGTGGTCGGGGCAGTGCATCTCTCTCAAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATGCTCTGCCCACCACGGGCTTCAAAATATGATGGTCCGTCCAGGTCATGTTCAGCGAGCTGGTG";
        // 2048 ~ 2648
        // let seq = b"GCACTGCCCCGACCACGGTTTTTTGTCAACATATTGAAGGTCGTCCAGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAATCCCGGCCACAGCGCCCGGCAATGGCGATTAAAACGCCACCATCATGATGGAACGTCGGGAACCACAGTGGGGCGTGTTGGGTGCAAACTCGAAAAGCTCAAGCAGAAAAGGGTAAAGATAACACGAGACAAACCAAACTGACCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGTTTATCCGCAAAACAGGTCGCCAGTATTCGTTTCCTGTGTTTATGCTAGTACAAAAGTTTTACATTTTAATATGCGAATTGCTCATAAAGTGACGTAAGGCGGAAATTTGCGCAACTGCGTAATAACATTTTTTTAACTTGTATAAAACTGTTCAACGTAATTTGCCACGGAATTCTTCGTCAACAATTAACTACCTGAGCGCGCCAACGCTATTCGACTGGTATCAGGCGAATGAATCCCTATAATTGCCGCGTTTGGCGCTTCGTCCGCCCCC";
        
        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score
            )
        });
        
        // 3353~3662
        let seq = b"CACTGCCCCCGACCACGGCTTCACACAATATTGAGGGTCAGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGATGGAGCGCGATATTTCCCACGACAGCTCAGCGAAACAAAAAAAAAAAATTCCGCCACATGCGCCCGGCATGGCGATAAAACGCCACCTCCAGATGGACGCGGGAAACCACAGTGGGGCGAGTTGGGTGCCAAACTGCGAAAGTCAAGCAGAAAAGGGTAAAGAAAACGAGCTAAACTAACTGCACAGGTGAAGCAGCGGGTGGGCAGGCATAATACATGTGTGATTGTTTG";
        
        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score
            )
        });
    }


    
    #[test]
    fn test_align_single_query_to_target3() {
        let ref_file = "/data/ccs_data/MG1655.fa";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
        );

        let aligner = &mut aligners[0];

        // aligner.mapopt.min_cnt = 2;
        aligner.mapopt.best_n = 500;
        // aligner.mapopt.min_dp_max = 10; // min dp score
        // aligner.mapopt.pri_ratio = 0.5; // min dp score
        aligner.idxopt;
        println!("{:?}", aligner.mapopt);
        println!("{:?}", aligner.idxopt);

        // tot len
        let seq = b"AAAAGAGAGAACGAAAAAAAGAGAGAGATGAATCCAGCGAAAACCATGTACGGCTTAAATTCAGCCACCTTATTATCGCCAGAATCCCTTTTCGCACGAACTGTCCTTTGATGCGAAACAAAATTTGGCGCATTCACTGTGCATACCGGCCGGTCAAGCAGCAGGTAACAATGCGGAGCGGCGCGCCTGCCGTCCATGATCGGGATTTCCGGCGCCGTTAGACTTCGAAAACAATGTTATCGATGCCCAAAGCCCCGGAGAGCAGCATTGAGTGCCATCGGTTGAAATCCGTACATCTCTCGTGACCAGACACGTACAGAGCATGGTATCACGCACAGAATTTGGCATTCGGCCGGGGAAATCTACCGTGATTCAAAGTCGGTGCGACGATAGATGACCCGGTGTTGCCGGCGCAGGGCGTAACGTCAGGGTGACTTTTTGCCGGTATGTAAACCACAAACCCAGTCGCCTGAACAGAACGTTTATGGTCCTTGTATTGATCATCGTATTATCTCAGCCAATTACCTATCCAACCGAAGTGTACTATACATTCGGCCGCAGTTTTAGCACAAAAGAGCCTCGAAACCCAAATTCCAGCAATTCTTATTCAGCTTGCTTACGAGGAATCTGGGAATCCAGATAATCCGGCTCTTCCATTTGCGCGCATTGTCCAATTCACGACTTTAGCTTCCGCTTCTGCTCCCTGGGTCAGCGGAGCCCATCCCATGCTGCTGTAGCGATCCATCACTGGCTGCTGAACCTGCTTATTGGTCACCAGATCTTCAGAACGAACCATGCTGAATAGCCCTGCCATTGCCGCATCCAGTCCGCTCAAGGATGTTCCCATAATACGCCCCCCCCCCCCCCAGTTACGCATCGACCTTGCTTCCATCAGGCTGACCAACGGGGCTGGAAATACGTTTCGGGTTGCTTCTTTCATCAGGCCAGACGTGACCCGTGCATCACCGCAAATCAAACGGTCTGCTGGGTAAATCCTGTTTGTAAGGAGCATGACGCCAGCATTTGTAATTTGCATGATCGGTAACCTGGGCATGATTCATAAACACCACTGCAAATTTTTGTGTCGTGCCTGGTCTACTAGTCGTAAAAATTGATCGCGGACAAATTTCGCCAGCATCTCTCTCAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGTGTCTGGGCGAATATTTCCGCGTTCAATTTTTTACGACATAGTAGCTCCCAGGCACGACAGCAAAAATTTGCAGTGGCTGTTTATGAATCTGCAGGTTTTACCGTCTGCAAAATTACAATGCTGGCGTCATGCTCGTACAAATCAGGATTTACCAGCGAGACGTTTTGATTTAGCGGTGATCGCCACGGGTCGACGTCTGCTGATGAAGAAAGAAAGCAACCCGAACGTATTTTTCTCCAGCCCGTGTCAGGCCTGATGGAAGCAAAGGTCGATGCGTGTAACGTGGGTATTATGGGAACATCCTTGAGCGGACTGGATGCGGCAATGGCAGGCTATTCACATGGTTCGTTCATTGAAGATCTGGTGACCAATAAGCAGGTTCAGCAGCCATGGATGGATCGCTACATGCAGCATGGGATGGCTCCGCTGAAGGAGCAAAGCCGGTGCTAAAGTCGTGAATGACAATGCGCCGCAAACTGCGAAGAGCCGGATTCCTGGATATCCCTGCATTCCTGCGTAAGCCAACTGATTAAAGAATTGACTGGAATTTGGGTTCGAGGCTCTTTGTGCTAACTGGCCCCCGAATGTATAGTACACTTCGGTTGGATAGTAATTGGCGCAGATATTTCATGATCAAACAAAGACACTTAAACGTATCGTTCAGGCGTCGGGTGTCGGTTTACATACCGGCAAGAAAGTCACCCGCGATTACGCCTGCGCCGGCCAACACCGGGGTCATCTAACGTTCGCACCGACTTGAATCCACCGGTTAGATTTCCCGGCCCGATCCAAAATCGTGCGTGATACCATGCGCTGTACGTGTCTGGTCAACGGCAGAGTACGGATTTCAACCCGTAGAGCATACCTCAATGCTGCTCTCGCGGGCTTGGGCAATCGATACATTGTTATCGAATTAACGCGCGGAAATCCCGATCATGACGGCAGCGCCGCCCGTTTGTATTACTCTGCTTGACGCCGGTATCGACGAGTTGAACTGCGCAAAAAATTGTCCGCATCAAAGAGACTGTTCTGTCGTTGATGGCGTAAGTGGGCTGATTTAAGCCGTACAATGGTTTTTCGCTGGATTCATCTCCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGATTAAAGAGAGAGATGAAATCCAGCGAAAACCATGTAGGCTTAACATTCAGCCACTATCGCCTCTTCGACAGAACAGTCTCTTTGATGCGAACAATTTTTTGGCGCAGTTCAACTCCGTCGATACCGCGTCCAAAGCAGCAGGTATACAAATACAGAGGCGGCGCTGCGTCATGATCGGGATTTTCCGGCGCGTTAACTTTCGATAACAAGTTATCGTGCCCAAGCCCCGGAGAGCAGCATTGAGGTGCTCTACGGTTGATAATCCTACATCATGCTCCGTTGACCAGAAACACGTACAGACATGGTATCACGCACAAGATTTGGCATCGGCGGGAAATCTACCTCGTGGATTCAACGTCGGTGCGACGTAGGTGTCCCCGGTGTTGGCCGCGCAGGGTCACGGCAGGTGACTTTCTTGCCGGTATGTAAACCGACACCCGTCCGCCTGAACGATACGTTTAAGTGTCCTTTGTTTGATCATCGTATTATCTCGCCAATTACCTATCCAACCGAAGTGTCTAATACATTCGGCGGGCATTTAGCACAAAAGAGCCGAAACCCAAATTCCCATCAATTCTTAAATCAGCTTGCTGACGCAGAATGCTGGATATCCAGAATCCGGGCCTTTCGCAGTTTGCGCCGCATTGTCATCAAGACTTAGCACGGCTTCGCTCCTGGGTCAGCGGAAGCCATCCCCATGCTGCTGTAGCGATCCTCACTGGCTGCTGAACCTGCTTATTGTCACCAAGATTCCTTCTTTGAACCGAACCATGCTGAATAGCCATCCATTGCCGCATCCAGTCCGCTCAGGAGTTCCCAAAATACCCACGTTACACGCATCGACCTTTTGCTCCACAGGCCTGACCACGGGCTGGGAAAATACGTTCGGTTGCTTCTTCTTCACAGGCCAGACGTGACCCCGTGGCGATCACCGCTAAATCAAACGTCTCGCTGGTAATCCTGATTTGCAGCGAGCATGACGCCAGCATTTGTAATTTGCAGATCGTAACTGGCATGATTCATAAACAGCCACTGCAAATTTGCTGTCGTGCCTGGTCTACAGTCGTAAAAATTGATCGCGGAAATGATTTCGCCCAGCATCCTCTCCTCAAACAACAACGGATGGAGGAGGAAAATAACAGGAAAAAAGAGAGAGATGCTGGGCGAATATTTGGATCAATTTTTTACGACCTAGTAGACCAGGCACGACAGCAAAATTGCAGTGCTGTTTATGAATCATGCCAGTTACGGACTGCAAATTACAAAGCTGGCGTCATGCTCGCTACAAATCAGGTTTACCCAGCGAGGATTTTGATTTATAGCGTGATCGCCACGGGTCACGTCTGGCCTGATGAAGAAGAAGCACACCGAACGTATTTTTCCCAGCGTGGTCAGCTGATGGAAGCAAAGGTCGAATGCAGTAGTAACGTGGGTATTATGGGAACACCTAGCGGACTGGATGCGGCAAGCAGTGGCTATTCAGCATGGTTTCGTTCATTGAAGATTGGTGACCAATAGCAGGTTCAGCAGCCAGTGAGGAGCGCACCAGCAGCATGGATGGCTCCGCTGACCCAGGAGCAGAAGCCCGGTTGCTAAATCGTGAATGACAATGCGCCGCAACTGCGGAAAGAGCCGGATTAGATCTGGATATCCCAGCATCCTGCGTAAGCAAGGCTGATTAAGAATTGACTTGGAATTGGGTTTCGAGGCTCTTTGTGCTAAACTGGGCCCGACCCGCCGAATGTATAGTACACTTCGTGGATAGGAATTTTGGGCGAGATAGATAGTGACCTATACGATGATTTCAAACAAAGGACACTTAAAACGTATCGTTCAGGCGTCGGGTGTCGGTTTACAATACCGGCAAGAAAGTCACCTGACGTTACGCCCTGCGCGGCCACACCGGGTCATCTATCGTCCCACCGACTGAATCCACCTAGATTGTCCCGGCCGATGCCAAATCTGTGCGTGATACCTGTTAAGCTCTTGTAACGTGTCTGGTCAACGAGCATGAGTACAGGATTTCAACCGTAAGAGCTCCTCAATGCTGCTCTCGCGGGCTTGGGCATCATAAATTTGTTATCGAATTAACGACGCCGGAATCCCGATTCATGGACGGCAGCGCCCGCTCCGTTTGTATACCTGCTGCTGGACGCCGGTATCGACGGAGTAGAAGCTGCGCCAAAAAATTTGTTCGCATCAAAGAGACTGTTCGTGGGTCCGAAGATGGCGATAATGGGCTGAATTTAAGCCGTACAATGGTTTTTTCGCGGATTCATCTCTCTCCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGAAATCCAGCGAAAAACCATTGTACGGCTTAATTCAGCCACTTTATCGCCATCTCAGACACGAACAGTCTCTTTTGAGAACAAATTTTTTGGCGCAGTTTCAACTCGTCCGATACCGGCGTCAAGCAAGCTGGTAATACATAGGTGCGGCGCTGCCGTCAAGATCGGGATTTCCGGCGTCGTTAACTGCGATAACAATGTTACGATGTCCCAAGCCGCGGAGAGCAGCATTGAGGTGCTCTAAACGGTTGAAAATCCGTAACATCATGCTCGTTTGACCAGACACGTAAACAGAGCATGGTATGCCACGACAGATTTGGCATCGGCGGGGGGGGGGGGGGGGGGGGGGGGGGGAAATCTACCGTGGATTCAAGTCGGTGCGACGAAGATGACCCCCGGTGTTGGCCGGCGCAGGGCGTAAACGTCAGGTGAACTTTCTGCCGGTATGTAAACCGACCACCCGTCGCTGAACGATACGTTTAGTGTCCTTGTTGATCATCGTAATTATCTCGCCAATTACCTATTCAACCGAAGTGTACTAGTACATTCGGCTGGCCAAGTTTAGCAAAAGAGCCTCCGAAACCCAAATCCAGTCAATTCTTAATCACTTGCTTACGCAGGTAGCTGGGATATCCAGATAATCCGCTCTTTTCGCAGTTTGCGCGCATTGTCATTCCACGACTTTAGCACCGGCTTCTGCTCCGGGTCAGCGGAGCCATCCCATGGCTGCTGGTAGCGATCCATCACTGCTGCTGAACCTGCTTATTGGTCACCAGATCTTCAATGAACGAACCATGCTGAAATAGCCACTGCCATTGCCGCATCCAGTCCGCTCAAGGATGTTCCCATATATCCCACGACACGCATCGACCTTTGCTTCCTCAGGCCTGACCACGGCTGGGAAAAATACGTTCGGGTTGCTTCCTTCTTCATCAGGCCAGACAGTGACCCGTGGGTTCCCCACCGCTAAATCAAGTCTCCTGGGTAATCTGATTTGTAAGCGAGCATGACGCCAGCATTTGTATTTTGCAGATCGGTAACCTGGCATGATTCAAAAGCCCTTGCAAATTTTTGCTGTCGTGCCTGGTCTACTAGTCGTAAAAAATTGAATCCGGAAATATTCGCCACATCTCTCTCACAACAACAACGAGGAGGAGGAAAAAAGGAGAGAGAGCTGGGCGAATAATTTCCGCGATCAATTTTTACGACTAGTAGTCCAGGCACGACGCAAAAATTTTGCAGTGCGGTTTATGAATCATGCCAGGTACGATCTGCAAATAAATGCGGCGTCATGCTCGCTACAAATGCGAGGATTTTACCCAGCGAACGTTTGATTTAGCGGGTTCGCCCTCGGGTCACGTGCTGCCTGATGAAGAAGAAGCAACCCGAACGGTATTTTCCCAGCCCTGGTCAAAGGCTGATGGAAGCAAAGTCGATGCGTGTAAACGTGGTAAAATTATGGGAACATCCTTGAGCGGACTGAGATGCGGCAATGGCAGTGGGCTATTCAGCATGGTTTCGTTCATTGAAGATCTGGTGACCCATAAGCAGGTCAGCAGCCAGTGATGGATCGTCTACCAGCAAGCATGGGAATGGCTCCGCTGACCCAAGGAGCAGAGCCGTTGCTAAAGCGTGAATGACATTGCGCCGCAAAACTGCGAAAGAGCCGGATTATCTGGATATCCCAGGCATTCCTGCGTAAGCAAGCTGATTAAATTGACGGAATTTGGGTTTCGAGGTCTTTGTGCTAAACTGGCCCCGCGAATGATAGTACACTTCGGTGGATAGGTAAATTTGGCGAGATAATACGATGATCAAACAAAGGACATTAAACGTATCGTTCAGGCGCGGGTGTCGGTTTAACATACCGGCAAGAAGTCACCTGACGTTAGCCTGCGCCGGCCAACAACCGGGGTCATCATCGTCGCACCACTGAATCACACCGTGTAGATTTCCCGGCCCGATGCCAATCTGTGCAGTGATACCTGCTCTGTACGTGTCTGGATCAACGGCATGATGTAACGGATATGCAACGTAGATAGAGCAGCCTCAATGCTGCTCTTCGCGGGCTTGGGCATCATAACATTGTTATCGTGGTTAACGCGCCCCGGAATCCCGATCATGGACATAGGCAGCGCCGCTCCGTTTGTATACCTGCTTGCTTGACGCCGGTATCGCGAGTTTGAACTGCGCCAAAAAATTTGTTCGCATCAAGAGACTTGTTCGTGTCGAAGATGCGATAAAGTGGGCTGAATTTAAGCCGTACAATGGTTTTCTCGCTGGATTTCATCTCTCTCAAACAACAACGGAGGAGAGGAAAAAAGAGAGGATGAATCCCAGCGAAAAACCATTTGTACGGCTTAAATCAGCCCACTTATCGCCATCTTCGACACGAACAGTCTCTTTGATGCGAACAATTTTTTGCGCAAGTTCAACTCGTCATACCGCGTCAAGCATGCAGGTATACAAACGGAGCGCGCTGCCGTCCATGATCGGGATTTCCGGCGCGTTAACTTCGATAACAATGTTATCGATGCCCAAGCCCGCGGAGAGCAGCATTGAAGGTGCTCTACGGTGAAATCCGTACATCAGCTCGTTGACCAACACGTACAGAAGCATGGTATCACGCACAGATTTGGCATCGGCCGGAATCTACCGGTGATTCAAGTCGGTGCGAAGATGAATGACCCCCGGTGTTGGCCGGCGGCAGGGCGTAACGTCAGGGTGACTTCTTGGCGGTATGTAAACCCACAGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGTCGCCGAACGATACGTTTAAGGTCTTTTGTTTGATCATCGATTATCTCCCAAATTACCTATCCAACCCGAAGTGTACCTATACATTCGGCGGCCGTTTAGCACAAAGAGCCCGAAACCCAATTCCAGTCCATTCTTATCAGCTTGCTTACCCAGGAATCTGGATATCCAGAAATCCGGCTCTTTCCCGAGTTTGGCGCATTGTCATTCACGACTTTAGCAACCGGCTTCTGCTCCTGGGTCCAGCGAGGCCCTTCCATGCTGCTGTAGCGATCCATCACCTGGCTGCTGAACCTGCTTATGGTCCCACCAGATCTTCAATGAACGAACCATGCTGGAATACCACTGCCTTGCCGCATCAGTCCGCTCAAAGGATGTTCCCATAATTCCCACAGTTACAGCATCGACTTGCTCCATCAAGCCTGACCACGGCTGGGAAATACGTTCAGGTTGCTTTTCTTCATGCAAGGCCCAGAACCGTGACCGTGGCGATCACCGCTAAAATCAAACGTCTCGCTGGGTAAATCCTGATTTGTAGCGAGCATGACGCCAGCATTTGTAATTTGCAGATCGGTAACTGGATGTTCATAAACACACTGCAAATTTTTGCTGTCGTGCCAGGTCTACTAGTCGTAAAAATTGATCGCGGAATATTCGCCCAGCATCTCTCTCAAACAACAACAGAAGGAGGAGGAAAAAAGAGAGAGATGCTGGGCGAATAATTTTGCCCGATCAATTTTTTACGACTAGTAAGAGCCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCAGCACGCAGCAAAAATTTGCAGTGGCTAGTTTATTGAATCATGCCAGGTTTCCGAACTGCAAATTACAAAATGCTGGGTCAGCGCTACAAAATCAGGTTTTAACCCAGCGAGACGTTTTGATTTAGCGGTGATCGCCACGGTCTCGTCTGGCCCTGATGAAAAGAAGCAACCCGAACGTATTTTCCAGCCCGTGGTCAGGCCTGATGGAAGCAACAGGTCGATGGCGGTACACGTGGGTATTATGGGAACATCCTTGAGCGGACTGTGCGGCAATGCAGTGCTATTCAGCATGGTTCGTCCATTGAAGATCTGGGACCAATAGCTGGTTCAGCATGCCAGTGATGGATCGCTACCAGCAGATGGATGGCTCCCGTACCCAGGAGCAGAAGCCGGTTGCAAAGTCGTGGAATGACAATGCGCCGCAAACTGCGAAGAGCCGGAAATCTGGATAATCCCCCAGCATTCCTGCATAAGCAAGCTGATTAAGAATTGACTGGAATTTGGGTTCAGGCTCTTGTGCTAAACTGTAGCCCCGCCGATGATAGTCACTTCCGTTGGCTAGGTAAATTTGGCAGATAATCGATGATCAAACAAAGGACACTAAACGTATCGTTCAGTGGGCGACGTGTCGGTTTACATACCGCAGAAAAGTCCACCCTGACGTTACGCCCTGGCGCCGCAACACCGGGTCATCTAACGTCGCACCGACTTGAAATCCACCGGAGATTCCGGGCCAGATGCCAAATCTGTGCGTGGATACACAAGCTCTGTACGTGGTCTGGTCAACGACATGATGTACGGATTTCAACCGTAGAGCACCTCAATGCTGCTCTCGCGGGCTTGGGCATCATAAACATTGTTATCGAAGTTAACGCGGCCGGAAATCCGATCATGGCGGCAGCGCCGCTCCGTTTGATACCTGTGCTGCTTGACGCCGTATCGACGAGGTTGAACTGCGCCAAAAAATTTGTTCCGCATCAAAGAGAACTGTTCAGTGTCGAAGATGGCGAAAGTGGGCTGAATTTAAGCCCGACAATGGTTTTCGCTGGATTTCATCTCTCTCAAACACAACGGAGGAGGAGGAAAAAAGAGAGAGATGAAATCCAGCGAAAAACATTGTACGGCTTTAAATTCAGCCCACCTTATGCATCTTCGACACGAACAGTCCTCTTGATGCGCAAATTTTTGGCGCAGTTCACTCGTCATACCGGCGTCAAGCAGTGTATCAAATATCGGAGCAGGCCTGCCGCCATGATCGGGATTTCCGGCGCGTTAACTTGATAACAATGAATCCGAATGCCCAAGCCCGCGAGAGCGCATTGTGGTGCTCTACGTGTTTGAAATCCCGTAACATCATGCTCGTTGACCGACAGACAGAGCCATGGTATCACGCACAGATTTGGCATCGGCCGGGAAATCTACAGGTGGATCAAGTCGGGCGACCGATAGATGACCAGTGTTTGGCCGCGCCAGGCGTAAACGTCAGGGTGACTTTCGCCGGTATGTAATCCGACACCCGTCCCCGCCTGAACGATACGTTTAAGTGTCCCTTTTGTTTTGATCATCGTATTACTCGCCAAATTACCTATCCAACCGAAATGGTACTATACATTGGCGGCCAGTTTTAGCACAAAGAGCCTCGAAACAAATTCCCGTTATTCCTTAATCAGCTTGCTACGCAGGATGCTGGGATATCCAGAAATCGGCTCTTTCCAGTTTGCGGCGCATTGTCCATTCACGGACTTTAGCAACCGCTTCTGCTCCTGGGTCAGCGGTGCCATCCCTGCCTGCTGTAGCGATCCATCACTGGCTGCTGAACCTGCTTTATTTGGTCACCAGATCTTAATGACGAACCATGCTGAATAGCCACCTGCCTTGCCGCTCCAGTCCGCTCAAGGTGTTCCCATAATACACTTACACGCATCGACCTTTGCTTCACAGGCCTGACCACGGGCTGGGACAAATACGTTCGGGTTGCTTCGCTCTTCATCAGGCCAGACGTGACCGGCGATCAACCGCTAAATCAAACGTCTCGCTGGGTAAATCGCTGATTGTAGCGAGCAGACGCCAGCATTTGTAATTTGCAGATCGGTAAACCTGGCATGATCATAAACACCACTGCAAATTTTTGCTGTCGTGCCTGTCTACTAAGTCGAAATTGATCGCGAATATTCGCCCAGCATTCTCTCCAACAACAACGGAGGAGAGGCAAAAAGGAGAGAGATGCTGGGCAATATTTCCCGCGATCAATTTTTACGACTAGTAGACAGGGGCACCGACAGCAAAAAAATTTGCAGTGGCTGTTTATGAATCATGCCAGTTACGATCTGCAAATTACAAATGCTGGCGTCATGCTCGCCTACAAATCAGGATTTACCCAGCGAACGTTTGATTTAGCGGTGATGCGCCCGGTGTCACCGTCTGGCCTGAGAGAAGAAGCAACCCGAACGTATTTTTCCCAGCCCGTGGCAGGCCTGTGGAAGCAAAGGTCGATGCGTGAACCAGTGGGTATTATGGTGACATCCTTTGAGCGGACTGATGCCGCAATGGCAGTGGCTATTCACATGGTACGTTCATTGAAGAACTGTGACCAATAAAAAGCCAGGTTCAGCAGCCAGTGATGGATCGCTACCAGCAGAGGGATGGCTCCGCTGACCCAGGAGCAGAAGCCGGTTGCTAAAGTCTGATGACAAGGCCGCAAAACTGCGAAAGGAGCCGGATATCTGGATATCCCACATCCTGGCGTAGCAAGGCTGATTAAGAATTGCTGAATTTGGTTCGAGGCTCTTTGTGCTAAACTGGCCGCCGAAATGTATAGTACATTGTTGATAGGTAAATTTTTGGCGAGATAATACGATGATCAAACAAAGGACACTTAAACAGTACGTTCAGGCGACGGTGTCGGTTTACATACCGCAAGAAAGTCACCCTGACGTTACGCCCTGCGCCAGCCAACACGGGGTCATCTATCTCCACGACCTTGATCCGACCGGTAGATTTCCGGCGATGCCAAATCTGTGCGTGATACATGCTCTGTACGTGTCTGGTAACGAGCATGATGTAGGTTTTCAACCGTAGATCACCTCAATGCTGCTCTCGCCGGGCTTGGGCATCGATAACATTGGTTATCAAAGTTAACCGCCGGGAATCCGTCTGGACGGCAGCGCGCTCGTTTGTATACCTGCTGCTTGACCCGTACGAGAGTTGAACTGCCCAAAAAATTTGTTCGCATCAAAGAGACTGTTCGTGTCGAAGATGGCGATAAGTGGGCTGAATTTAAAGCCGTACACATGGGTTTTTCGCTGGATTTCATCTCTCTCAACAAAACGGAGGAGAGGAAAAAAGAGAAGATGAAAATCCAGCGAAAAAACCATTGTCGGCTTAAATCAGCCCACTTATCGCCATCTTCGACACGAACAGTCTCTTTGATGCGAACAAATTTTTTGGCGCAGTTCAACTCGTCGAATACCGGCGTCCCAAGCACAGGTATACAAACGGAGCGGCGCTGCCGTCCTGGATCGGGATTTCCGCGCGTTAACTTCGTAACAATGTTATTCGAAAGCCCAAAGCCCGCGAGAGCAGCATTGAGGTGCTCTACGGTTTGAAACCCGTACATCCATGCCTCGTTGACCAGACACGTTACAGGAGCATGTATCACGCACAGATTTGGCTCGGCCGGGAAATCTACCGGTGATTCAAGTCGGTGTCGACGATAGATGACCCCGGTGTTGGCCAGGCGCCAGGCGTAACGTCAGTGACTTTCTTGCCGGTATGTAAACCCGACCCCACCCTCGCCTGAACGATAGCGTTTAAGTGTCCTTTGTTGATCATCGTATTATCTCGCCAATTACCTATCCAAACCGAAGTGTACTATCACATTCGGGGCCAAGTTTAGCACAAGAGCTCGAAACCCAAAGTTCCAGTCATTCTAATCAGCTTGCTTACGCAGGAATGGCTGGGATATCCAGATAATCCAGCTCTTTCGCAGTTTGCGCGCATTGTCATTCACGACTTTAGCAACCGGCTTCTGCTCCTGGGTCAGCGAGCCATCCCATGCTGCTGGAGGATCCATCACTGCCTGCTGAACCTGCTTATTGGTCACCAGTTCTTCAATGAACGAACCATGCTGCAATAGCACTGCCATTGCCGCATCCAGTCGCTCAAGATGTTCCATAATACCCACGTTACGACGCATGACCTTGCTTCCATCAGGCCTGACCACGGCTGGGAAAATACGTTCGGGTTGCTTCTTCTTCATCAGGCCAGACAGTGGACCCGTGGCGAGCACCGCTAACAAACTCTCGCTGGGAAATCCTGATTTGTAGCGAGCATGACGCAGCATTTTAATTTGCAGAGTCGGTAACCTGGCATGATGCAAAACAGCCACTGCAAATTTTTGCTGCGTGGCCTGGTCACTAGTCGTAAAAATTGATCGCCGGAATAAATATTCCGCCCAGCATCCTCTCAACAACAACGGAGGAGGAGGAAAAAAGAAAGAGATGCTGGGCGAATATTCCGCGATCATTTTACGACTAGTTGTCCTGGCACGACAGCAAAAAATTTGCAGTGCTGTTTATGTTGCTTGCCAGGTTACCAACCTCTGCAATTTACAAAGCTGGCGTCAGCTCGCTACAAAATCAGGATTTACCCAGCGAGACGTTTGATTACGCGGGATCGCCACGGTCACGTCTCCTGATGAAGAAGAAGCAACCCGAACGTATTTTCCACAGCCGTGTCAGGCCTGATGGAAGCAAAGGCGATGCGTGTACGTGGGTATTATGGGAACATCCTTGAGCGGACTGGATGCGGCAATGGCAGTGCTATTCACATGCGTTCATGAATGCATCCTGGTGACCAAATAAGCAGTTCATGCAGCCAGTGATGATCGCTACCAGCAGCTTGGGATGGCTCCGCTGACCAGGACAGAGCCGTTGCTAAAGTCGTGATGACAATGCGCCGCAAACCGCGAAAGCCGGATTATCTGGATATCCCAGCATTCCCTGCCGTAATGCAAGCTGATTAAGAATTGCTGGATTTGGGTTCGAGGCTCTTGTGCTAAACTGGCCGCCCCGAGTATAGTATATTTCGGTTGGATAGGTAATTTGCGAGATAATAACTTGATCAAACCAAAGGACACTTAAACGTATCGTTCAGGCGACGGGTGTGCGGTTACATACCGCAAGAACAGTCACCCTGACGTTACGCCCTGCGCCGGCCACACCGGGGTCACTATCGTCGCACACTTGATAATCCACCGTAGAATTTCCCCGGCGAGCCAAATTGTGCGTATACCATGCTCTGGGTAAACGTGTCTGTCAACGGCATGATGTACGATTTCAACCGTAGAGCACCTCAATGCTGCTCTGCGGGCTTGGCATCGATAACATTGTTACAAAGTTAACGCGCCGGAAATCCCGTCACTGGACCGCAGCGCCGTCTCCGTTTTGATACCTGCTGCTTGACGCAGGTATCGACGAGATGAACTGGCGCCAAAAAATTGTTCGCATCAAAGAGACCTGTTCGTGTCGAAGATGGCGTAAGTGGGCTGAATTTATGCCGTACAATGGTTTTTCGCTGGATTTCACCTCTCTCAATCAAACAACAAGGAGGAGGAGGAAAAAAGAGAGAGATGAAACCCCAGCGAAAAACCATGTACGCTAAATTTCAGCCCACTTTCGCCATCTTCGACACGAACAGTCTCTTGATGCGAACAAATTTTTTGCGCAGTCAACTCGTCGTACCGGCGTCAAGCAGCAGGTATACCAAACGGTGCGCGCTTGCCGTCCATGATCGGGATTTCCGGGCCGCGTAACTTTCGTAAACATGTTATCGAGCCCAAGCCCGCGAGAGCGCATTGAGGTGCTCCTTACGGTTAAATCCGTACATCATGCTCGTGACCAGACACGTACAGAGCATGTATCACGCACAGAATTTGGCAATCCGGCCCCCGGGAAATCTACCGCTGGATTCAAGTCGGTGCGACGATAGTGACCCCGTGTTGGCCGGGCCGCCAGGGCGTAACGTCAAGGGGGTCCTTTTTCTTGCCCCCCGGTATGTAAAACCCGACACCCGTCCCCCCGCCCCTGAACGAAACCCCGTTTAAGTGGTTGTTTGGGAACAATCGTAATTTCTCCCATTTACTAAATCCCCCCCGAAACGGGGGGGGGTTTTAAACAA";
        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score
            )
        });
        
        // // 3353~3662
        // let seq = b"CACTGCCCCCGACCACGGCTTCACACAATATTGAGGGTCAGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGATGGAGCGCGATATTTCCCACGACAGCTCAGCGAAACAAAAAAAAAAAATTCCGCCACATGCGCCCGGCATGGCGATAAAACGCCACCTCCAGATGGACGCGGGAAACCACAGTGGGGCGAGTTGGGTGCCAAACTGCGAAAGTCAAGCAGAAAAGGGTAAAGAAAACGAGCTAAACTAACTGCACAGGTGAAGCAGCGGGTGGGCAGGCATAATACATGTGTGATTGTTTG";
        
        // let mut mapping = aligner
        //     .map(
        //         seq,
        //         false,
        //         false,
        //         None,
        //         Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
        //         Some(b"t"),
        //     )
        //     .unwrap();
        // mapping.sort_by_key(|hit| hit.query_start);
        // mapping.iter().for_each(|hit| {
        //     println!(
        //         "qstart:{}, qend:{}, primary:{}, supp:{}, identity:{}, score:{:?}",
        //         hit.query_start,
        //         hit.query_end,
        //         hit.is_primary,
        //         hit.is_supplementary,
        //         MappingExt(hit).identity(),
        //         hit.alignment.as_ref().unwrap().alignment_score
        //     )
        // });
// 
    }

    // #[test]
    // fn test_align_single_query_to_target2() {
    //     let ref_file =
    //         "/data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.fasta";
    //     let fa_iter = FastaFileReader::new(ref_file.to_string());
    //     let targets = read_fastx(fa_iter);
    //     let aligners = build_aligner(
    //         "map-ont",
    //         &IndexParams::default(),
    //         &MapParams::default(),
    //         &AlignParams::default(),
    //         &OupParams::default(),
    //         &targets,
    //     );

    //     let mut bam_reader = rust_htslib::bam::Reader::from_path(
    //         "/data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.bam",
    //     )
    //     .unwrap();
    //     bam_reader.set_threads(10).unwrap();

    //     let aligner = &aligners[0];

    //     for record in bam_reader.records() {
    //         let record = record.unwrap();
    //         let record_ext = BamRecordExt::new(&record);
    //         let seq = record_ext.get_seq();
    //         for hit in aligner
    //             .map(seq.as_bytes(), false, false, None, None, None)
    //             .unwrap()
    //         {
    //             println!("{:?}", hit);
    //         }
    //         break;
    //     }
    // }

    #[test]
    fn test_align_single_query_to_target4() {
        let ref_file = "test_data/MG1655.fa";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
        );

        let aligner = &mut aligners[0];

        // aligner.mapopt.min_cnt = 2;
        aligner.mapopt.best_n = 100000;
        // aligner.mapopt.min_dp_max = 10; // min dp score
        // aligner.mapopt.pri_ratio = 0.5; // min dp score
        aligner.idxopt;
        println!("{:?}", aligner.mapopt);
        println!("{:?}", aligner.idxopt);

        // tot len
        let seq = b"AAAAGAGAGAGATAGGGGGCGACGAGCGCCAAACGGGCGGCAATTATAGGGATTTCATCGCCTGATACCAGTCGAATAGCGTTGCCGCGCGCTCAGGATGTTAATTGGTTGACAGAAGAATTCCCGGTGGCAAAATTACGTTGAAGATCAGTTTTATACAAGGTAAAAAATGTTATACGCAGTTGCGCAATTATCGCCTTTACGTCACTTTATGAGCATTCGCATATAAAATGTAAAACTTTTTGTAACTAGCATAAAACACAGAAACGAATACCTGGCGCTGGTCTTGCGATAAAGCAGGTAATGAGCAAACAACACAGCATGTATTAATTGCCCTGCCCACCCGCTGCTTCCACCTGGTCCAGTTTAAGGTTAGTCCTCGTTTACTTTACCCTTTTCTCGCTGAGCTTTCGCAAGTTTGGCACCCTCTCGCCCCACGCCTGTGGTTCCGACGGTCCACTGATGGTGGCGTTTATCGCCATGCCGGGGCGCAATGTGCCGGGAAATGCGCTGAGCTGTTCGCTGGAAATATCGCCGCATCCATCCTGCTTTTTTCCACCAGCCGCTGAACATGACCTGGACGACATCAATATTGTTGAAGCCGTGGTCGGGCAGTGCATCTCTCGCTCAAACAACAACGGAGGAGGAGAAAAAAGAGACAGATGCACTGCCCCGACCACGGCTTCAACACATATTGATGGTCGTCCAGGTCATTTCCACGAGTGGTGGAACAAAGCAGGAGGATGCGGCGTATTTCCAGGAACAGCTCACGCAATTCCCGGCACATGGCGCCCGGCAATGGCGAATAAAACGCCAACCATCATGAATGGACGTCGGAACCACATGGGCGATTGGTGCCAAACTGCGAAAGCTCAAGCGAGAAAAGGGTAAAGAAAAACGAGACTTATACCTAAATGACCAGGGAATGCAGCGGGTGGGGCAGGGCAATTAATTACATGGCTGTGATTGTTGCTCATACCGCTTTATCCCAAGACAGGTCGCCAGTATTCGTTTCTGTGTTTAATGCTATACAAAAGTTTATTTACATTTATATGCGAAATTGCTCATAAGTGAGGGTAAAGGCGATAATTTGCGCAACTGCGTAATAACATTTTTTTACCTTGTATAAAACTGATCAACGTAAAGTTGCCACCGGGATCTTCGTCAAAATTAACTCTGACGCGCGGCAACGCTTATTCGACTGGTATCAGGCGGATGAAATCCTATAAATTGCCGCGTTTGGCGCTTCCGTCGCCCCCTATCTCTCTCAAACAACAACGGAGAGGAGGAAAAAAGAGAGAAGATAGGGGGCACAGAAGCGCCCAACGCGCAATTAATAGGGATTCATCCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCATGATACCAGTCGAATAGCGTTGCGCGCGCTCAGAGTTACATTGTTGACGAAGAATTCCCGGGTGGCAAATTACGTTGATCAGTTTTATAACAAGGTAAAAAGTTATACGCAGTTGCGCAAATTATCCGCCTTACGTCACTTTAATGAGCAATTCCGCATAATAAAAATGTAAAACTTTGTAACTAGCATAAACACAGAAACGAATACTGGGGCGACCTGGTCTTGCGGAATAAGCGGTAATGAGCAAAAATCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGCTTCACCGGTCAGCGTTTAGGTTTAAGTCTCGTTTAATCTTTACCCTTTCTCGTTGAGCGAGCTTTCGAGTTTGCACCCAACTCGTCCCACTGTGGTTCCCGACGTCCATCATGTGGGCGTTTTAAAATCGCGCATGCCGGCGCATGTGGCCGGGAATTGCCTGAGCTGTTACGCTGGGAAATATCGCCGCATCCATCCTGGCTTTTTTTCCACAGCTACGCTGAACAGACCTGGACGACCACCATCAATATGTTGAAGCCGTGGTCGGGCGTGCATCTCTCTCAAACAACAACGGAGGAGGGAGGAAAAAAGAGAGAGATGCACTGCCCCGACCACGGTTTTTTGTCAACATATTGAAGGTCGTCCAGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAATCCCGGCCACAGCGCCCGGCAATGGCGATTAAAACGCCACCATCATGATGGAACGTCGGGAACCACAGTGGGGCGTGTTGGGTGCAAACTCGAAAAGCTCAAGCAGAAAAGGGTAAAGATAACACGAGACAAACCAAACTGACCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGTTTATCCGCAAAACAGGTCGCCAGTATTCGTTTCCTGTGTTTATGCTAGTACAAAAGTTTTACATTTTAATATGCGAATTGCTCATAAAGTGACGTAAGGCGGAAATTTGCGCAACTGCGTAATAACATTTTTTTAACTTGTATAAAACTGTTCAACGTAATTTGCCACGGAATTCTTCGTCAACAATTAACTACCTGAGCGCGCCAACGCTATTCGACTGGTATCAGGCGAATGAATCCCTATAATTGCCGCGTTTGGCGCTTCGTCCGCCCCCTATCTCTCTCAAACAACAACAGGAGGGGAGGAAAAAAGAGAGAGATAGGGTGGCGACGAAGCGCCAAACGCGGCAATTATAGGGATTTCATCCGCCTGATACCAGTCGAATAGCGTTGCCGCGCGCTCAGAGTTAATTGTTTGACGAAGAATTCCCGGTGGCAAATTACGTTGATCAGTTTTTATACAAGGTAAAAAAATGTTATACGCAGTTGCGCAAATTATCCGCTTTACAGTCCACTTTTAGAGCAATTCGCATAATAAAAGGTAAAACTTTTGTACTAGCATAAGAACACACACAGAAACAATACTGGCTGACCTGGTCTTGCGGATAAAGCGGTAATGAGCAAACATCACGCATGTATTAATTTGCCCCTCCCCACCGCTGCTTCAACCGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTAGTACCCTTTTCTCGCTTGAGCTTTGCAGTTTGGCACCCAACTCCCATGTGTTTCCGCGACGTCCATCATGATGGTGGCGTTTTATCGCATGCCGGGCGCATGTGGCGGGAATTGCGCTGAGCTGTTCGCTGGGACATATCGCCGGCATCCTCCCTGCTTTTTTCCCCAGCTCGCGCTGAACATGACCTGGACGACCAAGTCAATATTGTTGAGCCGTGGTCGGGTGCAAGTGCATCTCCTCAACAACAACGAGGAGGGGAAAAAAAAAGATGAGAGATGCACTGCCCCCGACCACGGCTTCACACAATATTGAGGGTCAGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGATGGAGCGCGATATTTCCCACGACAGCTCAGCGAAACAAAAAAAAAAAATTCCGCCACATGCGCCCGGCATGGCGATAAAACGCCACCTCCAGATGGACGCGGGAAACCACAGTGGGGCGAGTTGGGTGCCAAACTGCGAAAGTCAAGCAGAAAAGGGTAAAGAAAACGAGCTAAACTAACTGCACAGGTGAAGCAGCGGGTGGGCAGGCATAATACATGTGTGATTGTTTGGTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCATTACCGCTTTATCCCAAGACCAAGGTCGCCAGTATTCGTTTCTGGGTTTATGCTAGTACAAAAGTTTACATTTTATATGCGAATTGCTCCATAAAGTGGCCCCCCCCCCCCGTAAAGGCGGATAATTTGCGCAACTGCGTATAACATTTTTTACCTTGTATAAAACTGATCAAGTAATTTGCCACCGGAATTCTTTCGTCAACAATAACTCTGGCGCGCGGCAACGCTATTCGACTGGTATCAGGCGGAGAAATCCTAATAATTGCCGCGTTTGGCGCTTCGTCGCCCCTATCTCTCTCCAAACAAACAACGAGGAGGAGGAAAAAGAGAGAGATAGGGGCGACGAAGCGCCAAACGCGGCAAATTATAGGGATTCATCCGCCTGATAACCAGGTCGAATAGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCAGCGTTGCCCGCGCTCAATTTAATGTTGACGAAGAATTCCCGGTGGCAAATTAACACGTTGAATCAGTTTTATACAAGGTAAAAAATTTATAGCAGTTGCGCAATTAATCCGCCTTACGTCACTTATGAGCAATTCGCATATAAAATGCTAAAACTTTTTGTACTAGCATAAAACACAGAAACGAATACTGGCGACCTGGTCTCGCGGATAAAGCGGTAGATGAGCAAAACAATCACAGCACTGTATTAATGCCTGCCCCACCGCTGCTTCACCGAGTCAGTTTAGGTTTATAGTCTCGGTTATCTTTACCCTTTTCGCTTGAGCTTTCGCAGTTTGGCACCACACTGCCCATGTGGTTCCGCGTCCATCCATGATGGTGGCGTTTTTATCGCCATGCCGGGAGCGCATGTGGCCGGAATTGCGCTGAGCTGTTCGCTGGAAATATCGCCGCATCCGATCCTGCTTTTTTCCACCAGCTCGCTGAACATGACCGGACGACCATCAATATTGTGAAGCCGTGGTCGGGGCAGTGCATCTCTCTCAACAAACAAACAACAGGAGGAGAGAAAAAAGAGAAAGCATGCCCCGACCACGCTTCAACAATATTGAGTGGTCGTCCGGTCATGGTTCAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGCGATATTTCCCAGCGGAACACTCAGCGCAATTCCGGCCAAACAGCGCCCGGCATGGCTGATAAAACGCCACCATCATGATGGACGTGGGAACCACAAGTGGGGCGATTGGTGCCAACTGCGAAAGCTCAAGCGAGAAAGGGTAAAAGATAAACGATGACCTAAACCTAAACTGACCAAGTGAACAGCGGGTGGGGCAGGCAATTAATACATGTGTGATTGTTTGCCATTACCGCTTTATCCCAAGACCAGGTCGCCAGTATTCGATTTCTGTGTTTATGCTAGTACTAAAGTTTTACATTTATATGCAGAATTGCTCATAAAGTGAGTAAGGTCGGATAATTTGCGCACACTGCGTATAACATTTTTTACCTTGTATAAAACGATCCAACGTAATTTGCCACCGGGAATCTCGTCAACAATTAACTCTGAGCGCGCGGCAACGCTATTCGACTGTATCAGGCGGATGAAATCCCTATAATTGGCCGCGTTTGGCGCTTCGTCGCCCCCATCTCTCTCAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATAGGGGGCGACGAAGCGCCAAACGCGGCATATAGGGATTTCATCCGCCTGATACCAGTCCGAAAGCGTTGCCGCGCGCTCAGAGTAAGTTGACGAAGAATTCCCGTGTGGGCACAATAACGTGACAGTTTTATACAAGGTAAAAAATGTTATACAGCAGTTGCGCAAATTTCCCCTTTACGTCACTCACTTTATGAGCAATTCGCATATAAAATGTAAACCCTTTTGTACTAAGCATAAACAGCAGAAACGAATACCTGCGATACCTGGTCTTGCGGATAAACGGTAATGAGCACAACTCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGTCTTCAACCTGGTCATGTTTAGGTTTATCTCGTTTATCTTACCTTTCTCGCTTGAGCTTTCGCATATTGGCACCCAATCGCCCCACGTGGTTCCCGACGTCCATCATGAGGTGGGCGTTTTATCCCATGCCGGGCGGCAATGTGGCCGGAATTGCGCTGAGCTGTTCCGCTGGGAAAATCGCCGCATCCATCTAGCTTTTTTTCCACAGGCTCGCTGAACATGATACCTGGACGACCATCAATATTGTTGAAGCGGGTCGGCAGTGCATCTTCTCAAACAACAACCGTAGGAGGAGGAAAAAAGAGAGAGATGCACTGCCCCGACCAACGGCTTCACAATATTGATGGTCGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGTCGGCGATATTTCCCAGCGAACAGCTCAGCGCAATTCCCGGCCACATGCGCCCGGCATGCGATATAAAGCCACATCATGATGGACGTCGGAACCACATGGAGCGAGTTGGGTGCCAAACTGTGGAAAGCTCAGCGACAAAAGGGTAAAGATAAACGAGACTAAACCTAAACTGACCAGTGAAGCAGCGGGTGGGCAGGGCAATTAATACATGCTGTGATTGGTTTGCTCATACCGCTTACGCAAGACCAGGCCGCAGTATTCTTTCTGTGGTTTATGCTAGTACAAAAGTTTTACATTTTATATGCGAATTGCTCAAAAGTGACGTAAAGCGGATAATTTGCGCAACTCGTATAACATTTTTACCTTGGTATAAAACTGATCAACGTAATTGCCACCGGGAATCTTCGTCACAATTAAGACCTGAGCGCGCGGCAACGTATTCGACTGGTATCCAGCGGATGCAAATCCCTATAATTGCCGCGTTGGCGCTTCGTCGCCCCATTCTCTCAAACAACAACGGAAGGAGAGGAACAAAGAGAGAGCATAGGGGGGCGGGGGCGAACGAAGGCAAACAGCGGCAATTATAGAGGATTTCACTCCGCCTATACCAGTCGAATAGCGTTGCCGCGCGCCAGAGTTAATTGTGACGAAGAATTCCCGTGGCAAATTAACGATGATCAGTTTTATACAGGGTAAAAAATGTTATTACAGCAGTTGCGCAAATTATCCGCCTTACGTCACTTTATGCAGCAATTCGCATATAAAATGTAAAACATTTTGTACTAGCATAAACACAGAAACGAATACTGGCGACCTGGTCTTGCTCGGATAAAGCGGTTAATGAGCAAACAATCACCAGCAGTATTAATTGCCCTGCCCCCACCGCTGCCTTCACCTGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTCTCGCTTTGAGCTTCGCAGTTTGGCACCCAACTCGCCCCCACTGTGGTTCCCGACTCCCATCATGATGGTGGCGTTTTATCGCCATGCGGGCGCATGTTGGCCGGGAATGCGCTGACTGTTCGCTGGGAATATCGCCGCATCATCCTGCTTTTTTCACCAGCTCGCTGAACATGACCTGGACGACATAATAGTTGAAGCGTGTCGGCAGTGCATCTCTCTCAAAACAACAACGAGGAGGAGGAAAAAAGAGAAGATGCCTGCCCCCGCACAGGTTAACAATATGAGGTCGTCCAGGTCGATGTTCAGCGGACTGGTGGGAAAAAAAAAAAAAGCAGGATGGAATGCGCATATTTCCCAACAGCTCAGCGCAATTCCCGGCCCACAATGCGCCCGGCATGGCGATAAAACGCCACCACAGAGGACGTCGGAACCACAGTGGGGCGAGTTGGGTGCCAAATCTCGAAAGCCAGCGAGAAAAGGGTAAAGATAAAACGAGACTAAACCTAAACTGACCAGGTGAAGCAGCCGGGTGGGGCAGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCCGGCTTTATCCCAGACCAGGTCGCCAGTATTCGTTTCTGTGTTTATGCTAGTACAAAAGTTTAACATTTTATAGCGATTGTTCATAAAGTGACGTAAAGGCGGATAATTTGCCAACTGCGTATAACGATTTTTACCTTGTATAAAACTGATCAACGTATTTGCACGGGAATTCTTCGTCCCCAATTGAACTCTGAGCGCGCGCACGCTATTCGACCTGGTATCAGGCGGATGAAATCCGTATAATTGCCCGTTTGCGCTTCGTCGCCCCTATCTCTCTCAAACAACAACGGAAGGAGAGGAAAAAAAGAAGAGATAGGGGGCGACGAAGCGCAAACGCGCAATTTATAGGGATTTCTCCAGCCTGATACCAGTCGACATAGCGTTGCCGCGCGCTCAGAGTTAATGTTGACGAAGAATTCCCGGGTGCAAATTACGTTGATCAGTTTTTTTTTATACAAGTAAAAAATGTTATACGCAGTGGCAAATTTCCGCCATTTACGTCACTTTATGAGCAATTCGCATAAAAAATGTAAAGTTGCATTGTAAAATCTTTTGTACAGCATAAACACAGAAACGAATAGCTGGCGAACCTGGTCTGCGGATAAGAGCGGTAATGAGCAACACACAAGCATGGTATAATTGCCCTGCCCCACCCGCTGCTCACCTGGTCGTTTAGGTTAGTTCGTTTACTTTAACCCTTTCTCGCTTGAGCTTCGCAGTTTTGGCACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGGGCGTTTTACGCCATGCCCGGGCGAGTGGCCGGGAATTGCGCTGAGCTGTTCGCTGGGAATATCGCCGCACCATCCTGCTTTTTTCCACCAGGCTCCGCTGGAACATGCCTGGCGACCAATCATATTGTTGAAGCCGTGGTCGGGGCAAGTGCATCTCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGCATGCCCGACCACGCTTAACAATATTGATGGTCGTCCAGGTCATGTTCAGCGGAGCTGGTGGAAAAAAGCAGGATTGGATGCGGCGATATTTCCCAGCGAACAGCTCAGCGCAATTCCGCCACATGCGCCGGCATGGCGATAAAACGCACCATCAATGTGGACGTCGCTGAACCGAAGTGGGGCGAGTTGGTGCCAAACTGCGGAAAGCTCAAGGAGAAAAGGGTAAAGATAAACGAACAAACTAAATGACCAGGTGAAGCACGGTGGGGCCAGGGAATTAATACAGCTGTGATTGTTTGTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTCATTACCGCTTTATCCGAAGACCAGGTCGCCATATTCGTTTCTGTGGTTTAGCTAGTACAAAGTTTTACATTTTATATGCGAATTGCTCATAAAGTGACGTAAAGGCGAATAATTTGCGCACTGCAGTATAACATTTTTTACTGTATAAAACTGATCAACGTAATTTGCCACCGGGATTCTTCGTCAACAATTTACTCTGAGCGCGCGGCACGCTATTCGACTGGTAATCAGGCGGATGAAATCCCTATAATTGTCCCCCCCGGGGGGCGTTGGCGCTCGCGCCCCCTATCCTTCTCAACAACAACGGAGAGAGGAAAAAGAGAGAGATAGGGGGGGGCGACAAGCGCCAAACGCGGCAATTATAGGGATTTCCATCCGACCTGATACCATCGAATAGCGTTTGCCGCGCGCTCAGAGTGTAATTGTTGACAGAATTCCCGGTGGCAAATTACGTTGATCAGTTTATACAAAGGTAAAAATGTTATACGCAGTTGCGCAAATTATCCGCCTTTAGTCCTTATGAGCAATTCGCATAAAACCTAAAATGTAAAACTTTAGTACTTGCAATAAACACAGAAACAATACTTGGCGACCTGGTCCCTCTGCGGATAAAGCGGTAAATGAGCAAACAACACAGCATGTATTAATTGCCCTGCCCCAACCCGCTGCTTCACTGGTCAGTTTAGGTTTAGTCTGTTTTATCTTTACCCACTTTTCTCGCTGAGCTTCCGCAGTTTGGCACCCAACTCGCCCACTGTGTTCCCGACGTCCATCATGATGGTGGCGTTTATCGCCATGCCGGGCGCATGTGGCCGGAATTTGCGCTGAGCTGTTTGCTGGGGAAATATCGCTGCATCCATCCTGCTTTTTTCCACGAAAAAAACAGCTCGCAGAAATGACCTGGACGACCATCATATTAGTTGAGCCGTGGTCGGGGCATGCATCTCTCTCAAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGCACTGCCCCGAACACGGCTTCAACAATATTGATGGTCGTCCAGGTCCCCATGTTCAAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGTGCGATATTTCCCATCCCCCCCCCGAGCGAACAGCTCAGGCAATTCCCGGTCCACATGCGCCCGGCATGGCGGATAAAACGCCACCCGATCATGATGGACGTTCCAGGGAACCACATGGGGCGAGTTGGTGCCAAACTGCGAAAGCTCAAGGCGAGAAAAAGGGTAAAGATAAACGGACAAACCTAAACTGTCCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGGCTGTATTGTTTGCCATTAGCCCCCCCCCCCCCCCCCATTCCGCTTTTATCCCGCAAGACCAGGTCGCCATATTCGTTCTGTGTTTAATCTAGTACACAAGTTTACATTTATAGCTAATGCGAATTGCTCATAAGTGACGTAAAGGCGGATAAATTCGCAACTGCGTATAACATTTTTTACTTGTATAAAATGATCAACGTAATTTGCCACCGGGAATCTTCGCAACAATAACTCTGAGCGCGCGGCACGCTATTCGACTGGTATCAGGCGGTGAAATCCCTAAATTGCGCGTTTGGCGCTTCGTCGGCCTTCCTCTCAAACAACACGAGGAGGAGGAAAAAAAAAAAAGAGAGAGATAGGGGGCGACGAAGCGCAAACGCGCAATTATAGGGATTTATCCGCCTGATACCAGTCGATAGCAGTTGCCGCGTCGCTCAGAGTTAATTGTTGACAAGAATTCCGCGGTGGCAAATTACGTGATCAGTTTTATACAAGAGTAAAAAAGTTTATACGCAGTTGCGCAAATTATCCGCCTTTACAGTCACTTATGAGCAATTCGCATATAAATGTAAAACTTTTGTACTAGCATAAACACAGAACACGAATACTGGCGACCTGGTCTTGCGGAAAAGCGGTAATGAGCAAACAAATCACAGCATGTATTATTGCCCTGCCCCACCCGCTGCTTCAAACCTGGTCAGTTTAGGTTTAGTCTCGTTTATCTTACCTTTTGCTCGTTGAGCTTTCCAGTTTGGCACCCAACTCGCCCACTGTGGTTCCCGACGTCCATCATGATGGTGGCTGTTTTATCGCCCATGCCGGGCGCATGTGGCCGGAATTGCGCTGAGCTGTTCGCTGGGAAACTATCGCCGCATCCATCCGCTTTTTTCCACCATGCTCGCTGAACATGACCTGGACGCACCATTCATATTGTTGAAGCCGTGGTGGGGCATGCATCTCTCTCAAACAACAACGGAGGAAGAGGAGAAAAAGGAGAGAGATGCACTGCCCGACCACGGCTTCAACAATATTGAGGTCGTCAGGTCATGTTCAGCGAGCTGGTGGAAAAAAGCAGGTGGCAGGATGGGCGGCCGATATTGCCGAGCGAACAGCTCAGCGCAAATTCCGGCACATGCGCCCGGCATGGCGAATAAAACGTGCACCATCATGATGGACAGCGGGAACCACAGTGGGGCGAGTTTTGGGTGGCCAACTGCGAAAGCTCAACGAGAAAAGGGTAAAGATAAACGAGACTAAACTAAACTGACAGGTGGAAGCAGCGGTGGGCAGGCAATTAATACAGCTGTGATTGTTTGCTCATTACCGCTTTTCCGCAAGACCAGGTCGCAAGTATTCGTTTCTGTGTTTATGCTAGTAACAAAAGTTTTACTTTTATATGCGAAAATTGCTCAAAAGTGACAGTAAAAAGGCGGATAATTTTGCGCACATGTCGTATAACATTTTTACCTTGTATAAAACTGACAACGTAAATTTGCCACCAGGGAATTCTTCGTCCAACAATTTAACTCTGAGCGCGCGGCAACGCTATTCGACTGTATCAGGCGGATGAAATCCCTATAATTGCCGGCTGATTGGCTGCTTCGTCCGCCCCCTACGTCCTCAAACAAAACGGAGGAGGAGGAAAAAAAGAGAGATAAGGGGCGACAAGCGCAAAACCGGCAAATTAAGGGATTTTCATCCGCCTGATACCAGTCAGGAATAGCGTTGCCGAGCGCTCAGAGTTAATTGTGACAAGAATTCCGGTGGCACAAATTATTGTCAGTTTATACAAGGTAAAAAATTTATACGCAGTTGCGCAAATTATCCGCCTATTACGTCACCTTTAGACAATTCGCATAAAAATGTAAAACTTTTGTACTAAGCATAAACACAGAAACGAATACTGGCGACTTGGTCTTGCGGATAAAGCGTAATGAGCAAACAATCACACATGTATTAATTGCCCTGCCCCAGCCGCTGCTTCACCAGGTCCAGTTTAGGTTTAGTCTCGTTATCTTTACCTTTTCTGCTTGAGCTTTCGCAGTTTGGCACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGTGCGTTTTATCGCATGCCGGGCGCATGGCCTGGAATTGCGCTGAGCTGTTCAGCTGGGAATATCGCCGCATCCACCTGCTTTTTTCCCACAGCTCGCTGAACAGACCTGGAACGACCATCATATTGTTAGAAAGCCGTGGTCGGGGCCGTGCATCTCTCTCAAACCAACAACGGAGGACAGAAAAAAAGAGAATGCACTGCCCGCCACAGGCTTCAACAATATTGATGGTCGTCCAGGTCATGTTTCAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGCGATATTTCCCCAGCGAACAGCTCAGCGCAATTCCCGGCCAATGCGCCCGGCATGGCGATAAAACGGCCAACCATCATGATGGACGTCGGGACCCATGGGGGCGAGTTGGGTGCCAAACTGCAAGCTCAAGCAGAAGGGTAAAGATAAGATAAACGAGACTAAACCTAAACTTGACCAATGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTAGCTTTAATCCGCAAGACCAGGTCGCCAGTATTCGTTTCTGTGTTTATGCTAGTACAAAAGTTTTATGACATTTTTATATGCGAATTGCTCAAAAGTGACGTAAAGGGGATAAATTTGCGCAACTGGCGTATAACATTTTTTACTTGGTATAAAGCCTAAAACTGATCAACCGTAATTTGCCACCGGGAATTCTTCGTCAACAATTAAACTCTGAGCGCGCGGCAACGCTATTCCGACTGTAATCAGCGATGAAATCCCTATAATTGCCCGCGTTTGCCGCTTCTCGCCCTATCCTCTCAAACAACAACGGAGGAGGGGAAAAAAGAGAGAGATAGGGGCGACGAAGGCCGCCAAACGCGGCAATTATAGGGATTTTCATCCGCCTGATACCAGTCGAATAGCGTTTGCCCGCGCGCTCAGATTAATTGTTGACGAAGAATTCCCGGTTGGCAAATTACTCTGATGCATGCTTTTAGACAAGGTAAAAAGATGTTATACGCAGTTGCGCAAATATCCGCCTTTACGTCACTTTATGAGCAATTCGCATCATAAAATGTAAAACTTTGTACTAGCATAAACAAGAAACGAATACTGGCGAACCTGGTCTGCGGATAAAGCGGTAATGACAAACAATCACAGCATGTATTAATTGCCCTGCCCCACCGCTGCTTCACCGGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTTCTCGCTTGAGCTATTCGCAGTTTGGCACCCTTCTCGCCCACAGAGAAGAGGTGGTTGCCCGACGTCCCATCATGATGGTGGCTTTTATCGCCATGCCGGGCGCATGGGCCGAATAGCGCTGAGCTGTTCGCTGGGGAATATCGCCCGCATCCAATCCTGCTTTTTTCCACGAGCTCGCTGAACATGACCTGGACGACCACAATATTGTTGAAGCCGTGGGTCGGGGGCCAGTGCGATCTCTCTCAACAAAAACGGAGGAGGAGGAGAAAAAAGAGAGGAGCACTGCCCCGACCACGGCTTCACAATATTTGATGGTCGTCCAGGTCATGTCAGCAGAGCTGGTGGAAAAAAGCGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAAATTCCCGGCCAACATGCGCCCGGCAATGCGATAAACCGCCACCAATCATGATGGACGTCGGGGACCACAGTGGGGCGAGTTGGGTGCCAAAACTGCGAAAGCTCAAGCGAGAAAAGGGTAAAGATAAACGAGACAAACCTAAACTGGACCAGGTGAAGCAGCGTGGACAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGCTTTATCCGCAAGACAGGTCGCCAGTAATCGTTTCTTGTGTTTAGCTAGTACAAAAGTTTTTATTTTATATGCGAATTGCTCATAATGTGACTGTAAAAGCGGATAATTTGCGCAACTGCGTATAACATTTTTTACCTTGTATGAAAAACTGATCAACGTAATTTGCCAACCGGGAATCTTCGTCAACAAATAACTCTGAAGCAGCGCGCAACGGCTATTCCGACTGGTATCAGGGCGGCATGAAATCCCTATAATTGCCGCGTTTGGCGCTCGCGCCCCCTATCTCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGCAGATAGGGGCGAGAAGCGCCAAACGCGGCAATTATAGGGATTTCATCCGCCTGATACCAGTGAAATAGCGTTGCGCGCGCTCAGAGTAATTTTGAGAAGAATTCCCGGTGGCAATTACGTGATCAGTTTTATACAAGGTAAAATGTTATACGCAGTTGCCAAATTATCCGCCTTTACGTCACTTTTATGGCAATTCGCATATACAAATGTAAAACTTTTGTACTAAGCATAAACACAGAAACGAATAGCTGGCCGACCTGGTCCTTTGCGGATAAGGGGAAGCGGTAAATGAGCAAACATCACAGCATGTATTAATTGGCCCTGCCCCCCGCTGCTTCACCGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTTCTCTGAGCTTTTCGCAGTTTGGACCCAAATACTCGCCCCACTTGTGGTTCCCGACGTCCATCGATGATGGTGGCGTTTTTATCGCCAATGCCGGGCGCATGTGGCCAGGGAAATTGCGCTGAGCTGTTCGCTGGGAATATCCGCCGCCAATCCTCCTGCTTTTTTTCCACCAGCTCGCTGAACATGACCTCGGACGCCATCAATATTGTTGAAAGCCGTGGTCGGGGCAGTGCATCTCTCTCAAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATGCTCTGCCCACCACGGGCTTCAAAATATGATGGTCCGTCCAGGTCATGTTCAGCGAGCTGGTG";

        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}, strand:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score,
                hit.strand
            )
        });

        // aligner.mapopt.min_cnt = 2;
        aligner.mapopt.best_n = 10;
        // aligner.mapopt.min_dp_max = 10; // min dp score
        aligner.mapopt.pri_ratio = 0.5; // min dp score
        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();

        println!("---------------------------------------------------");
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}. strand:{:?}, qlen:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score,
                hit.strand,
                hit.query_len
            )
        });
        

    }
}
