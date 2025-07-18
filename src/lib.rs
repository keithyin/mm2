pub mod align_processor;
pub mod bam_writer;
pub mod cigar_cvt;
pub mod mapping_ext;
pub mod params;
use std::{
    cmp,
    collections::HashMap,
    ops::{Deref, DerefMut},
    thread,
};

use crossbeam::channel::{Receiver, Sender};
use gskits::{
    dna::reverse_complement,
    ds::ReadInfo,
    fastx_reader::{fasta_reader::FastaFileReader, fastq_reader::FastqReader},
    phreq::phreq_list_2_quality,
};
use mapping_ext::MappingExt;
use minimap2::{Aligner, Built, Mapping};
use params::{
    AlignParams, IndexParams, InputFilterParams, MapParams, OupParams, TOverrideAlignerParam,
};
use rust_htslib::bam::{
    record::{Aux, AuxArray, Cigar, CigarString},
    Read,
};

pub use gskits;
pub use minimap2;
pub mod aligned_pairs;

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
    threads: usize,
) -> Vec<NoMemLeakAligner> {
    let per_target_threads = (threads / targets.len()).max(1);

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
                    .with_index_threads(per_target_threads)
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

/// does not work now
pub fn build_aligner_v2(
    preset: &str,
    index_args: &IndexParams,
    map_args: &MapParams,
    align_args: &AlignParams,
    oup_args: &OupParams,
    targets: &Vec<ReadInfo>,
) -> NoMemLeakAligner {
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

    let (seqs, ids): (Vec<Vec<u8>>, Vec<Vec<u8>>) = targets
        .iter()
        .map(|read_info| {
            (
                read_info.seq.as_bytes().to_vec(),
                read_info.name.as_bytes().to_vec(),
            )
        })
        .unzip();

    let mut aligner = aligner
        .with_index_threads(4)
        .with_cigar()
        .with_sam_out()
        .with_sam_hit_only()
        .with_seqs_and_ids(&seqs, &ids)
        .unwrap();

    // https://github.com/lh3/minimap2/blob/618d33515e5853c4576d5a3d126fdcda28f0e8a4/options.c#L43
    aligner.mapopt.best_n = 5; // top best_n chains are subjected to DP alignment
                               // aligner.mapopt.pri_ratio = 0.1; // secondary-to-primary score ratio to output secondary mappings

    NoMemLeakAligner(aligner)
}

pub fn query_seq_sender(
    filenames: &Vec<String>,
    sender: Sender<ReadInfo>,
    input_filter_params: &InputFilterParams,
    oup_params: &OupParams,
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
                            &oup_params.pass_through_tags,
                        ))
                        .unwrap();
                }
            }
        } else if filename.ends_with("fq") || filename.ends_with("fastq") {
            let fa_iter = FastqReader::new(filename.to_string());
            for mut record in fa_iter {
                if let Some(suffix) = &qname_suffix {
                    record.seq.push_str(suffix);
                }
                let rq = phreq_list_2_quality(record.qual.as_ref().unwrap()).unwrap_or(0.);
                record.rq = Some(rq);
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
            Some(String::from_utf8(reverse_complement(seq.as_bytes())).unwrap())
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

    if let Some(nn_) = &query_record.nn {
        if is_rev {
            let nn_ = nn_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"nn", Aux::ArrayU8(AuxArray::from(&nn_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"nn", Aux::ArrayU8(AuxArray::from(nn_)))
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

    cigar_str.extend(cigar_cvt::mapping_cigar2htslib_cigar_str(mapping_cigar).0);

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
    use std::time::Duration;

    use bio::alignment::pairwise::Scoring;
    use gskits::fastx_reader::read_fastx;
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
            10,
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
            10,
        );

        let aligner = &mut aligners[0];

        // aligner.mapopt.best_n = 100;
        // aligner.mapopt.pri_ratio = 0.1; // min dp score

        // tot len
        let seq = b"AAAAGAGAGAGATAGGGGGCGACGAGCGCCAAACGGGCGGCAATTATAGGGATTTCATCGCCTGATACCAGTCGAATAGCGTTGCCGCGCGCTCAGGATGTTAATTGGTTGACAGAAGAATTCCCGGTGGCAAAATTACGTTGAAGATCAGTTTTATACAAGGTAAAAAATGTTATACGCAGTTGCGCAATTATCGCCTTTACGTCACTTTATGAGCATTCGCATATAAAATGTAAAACTTTTTGTAACTAGCATAAAACACAGAAACGAATACCTGGCGCTGGTCTTGCGATAAAGCAGGTAATGAGCAAACAACACAGCATGTATTAATTGCCCTGCCCACCCGCTGCTTCCACCTGGTCCAGTTTAAGGTTAGTCCTCGTTTACTTTACCCTTTTCTCGCTGAGCTTTCGCAAGTTTGGCACCCTCTCGCCCCACGCCTGTGGTTCCGACGGTCCACTGATGGTGGCGTTTATCGCCATGCCGGGGCGCAATGTGCCGGGAAATGCGCTGAGCTGTTCGCTGGAAATATCGCCGCATCCATCCTGCTTTTTTCCACCAGCCGCTGAACATGACCTGGACGACATCAATATTGTTGAAGCCGTGGTCGGGCAGTGCATCTCTCGCTCAAACAACAACGGAGGAGGAGAAAAAAGAGACAGATGCACTGCCCCGACCACGGCTTCAACACATATTGATGGTCGTCCAGGTCATTTCCACGAGTGGTGGAACAAAGCAGGAGGATGCGGCGTATTTCCAGGAACAGCTCACGCAATTCCCGGCACATGGCGCCCGGCAATGGCGAATAAAACGCCAACCATCATGAATGGACGTCGGAACCACATGGGCGATTGGTGCCAAACTGCGAAAGCTCAAGCGAGAAAAGGGTAAAGAAAAACGAGACTTATACCTAAATGACCAGGGAATGCAGCGGGTGGGGCAGGGCAATTAATTACATGGCTGTGATTGTTGCTCATACCGCTTTATCCCAAGACAGGTCGCCAGTATTCGTTTCTGTGTTTAATGCTATACAAAAGTTTATTTACATTTATATGCGAAATTGCTCATAAGTGAGGGTAAAGGCGATAATTTGCGCAACTGCGTAATAACATTTTTTTACCTTGTATAAAACTGATCAACGTAAAGTTGCCACCGGGATCTTCGTCAAAATTAACTCTGACGCGCGGCAACGCTTATTCGACTGGTATCAGGCGGATGAAATCCTATAAATTGCCGCGTTTGGCGCTTCCGTCGCCCCCTATCTCTCTCAAACAACAACGGAGAGGAGGAAAAAAGAGAGAAGATAGGGGGCACAGAAGCGCCCAACGCGCAATTAATAGGGATTCATCCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCATGATACCAGTCGAATAGCGTTGCGCGCGCTCAGAGTTACATTGTTGACGAAGAATTCCCGGGTGGCAAATTACGTTGATCAGTTTTATAACAAGGTAAAAAGTTATACGCAGTTGCGCAAATTATCCGCCTTACGTCACTTTAATGAGCAATTCCGCATAATAAAAATGTAAAACTTTGTAACTAGCATAAACACAGAAACGAATACTGGGGCGACCTGGTCTTGCGGAATAAGCGGTAATGAGCAAAAATCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGCTTCACCGGTCAGCGTTTAGGTTTAAGTCTCGTTTAATCTTTACCCTTTCTCGTTGAGCGAGCTTTCGAGTTTGCACCCAACTCGTCCCACTGTGGTTCCCGACGTCCATCATGTGGGCGTTTTAAAATCGCGCATGCCGGCGCATGTGGCCGGGAATTGCCTGAGCTGTTACGCTGGGAAATATCGCCGCATCCATCCTGGCTTTTTTTCCACAGCTACGCTGAACAGACCTGGACGACCACCATCAATATGTTGAAGCCGTGGTCGGGCGTGCATCTCTCTCAAACAACAACGGAGGAGGGAGGAAAAAAGAGAGAGATGCACTGCCCCGACCACGGTTTTTTGTCAACATATTGAAGGTCGTCCAGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAATCCCGGCCACAGCGCCCGGCAATGGCGATTAAAACGCCACCATCATGATGGAACGTCGGGAACCACAGTGGGGCGTGTTGGGTGCAAACTCGAAAAGCTCAAGCAGAAAAGGGTAAAGATAACACGAGACAAACCAAACTGACCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGTTTATCCGCAAAACAGGTCGCCAGTATTCGTTTCCTGTGTTTATGCTAGTACAAAAGTTTTACATTTTAATATGCGAATTGCTCATAAAGTGACGTAAGGCGGAAATTTGCGCAACTGCGTAATAACATTTTTTTAACTTGTATAAAACTGTTCAACGTAATTTGCCACGGAATTCTTCGTCAACAATTAACTACCTGAGCGCGCCAACGCTATTCGACTGGTATCAGGCGAATGAATCCCTATAATTGCCGCGTTTGGCGCTTCGTCCGCCCCCTATCTCTCTCAAACAACAACAGGAGGGGAGGAAAAAAGAGAGAGATAGGGTGGCGACGAAGCGCCAAACGCGGCAATTATAGGGATTTCATCCGCCTGATACCAGTCGAATAGCGTTGCCGCGCGCTCAGAGTTAATTGTTTGACGAAGAATTCCCGGTGGCAAATTACGTTGATCAGTTTTTATACAAGGTAAAAAAATGTTATACGCAGTTGCGCAAATTATCCGCTTTACAGTCCACTTTTAGAGCAATTCGCATAATAAAAGGTAAAACTTTTGTACTAGCATAAGAACACACACAGAAACAATACTGGCTGACCTGGTCTTGCGGATAAAGCGGTAATGAGCAAACATCACGCATGTATTAATTTGCCCCTCCCCACCGCTGCTTCAACCGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTAGTACCCTTTTCTCGCTTGAGCTTTGCAGTTTGGCACCCAACTCCCATGTGTTTCCGCGACGTCCATCATGATGGTGGCGTTTTATCGCATGCCGGGCGCATGTGGCGGGAATTGCGCTGAGCTGTTCGCTGGGACATATCGCCGGCATCCTCCCTGCTTTTTTCCCCAGCTCGCGCTGAACATGACCTGGACGACCAAGTCAATATTGTTGAGCCGTGGTCGGGTGCAAGTGCATCTCCTCAACAACAACGAGGAGGGGAAAAAAAAAGATGAGAGATGCACTGCCCCCGACCACGGCTTCACACAATATTGAGGGTCAGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGATGGAGCGCGATATTTCCCACGACAGCTCAGCGAAACAAAAAAAAAAAATTCCGCCACATGCGCCCGGCATGGCGATAAAACGCCACCTCCAGATGGACGCGGGAAACCACAGTGGGGCGAGTTGGGTGCCAAACTGCGAAAGTCAAGCAGAAAAGGGTAAAGAAAACGAGCTAAACTAACTGCACAGGTGAAGCAGCGGGTGGGCAGGCATAATACATGTGTGATTGTTTGGTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCATTACCGCTTTATCCCAAGACCAAGGTCGCCAGTATTCGTTTCTGGGTTTATGCTAGTACAAAAGTTTACATTTTATATGCGAATTGCTCCATAAAGTGGCCCCCCCCCCCCGTAAAGGCGGATAATTTGCGCAACTGCGTATAACATTTTTTACCTTGTATAAAACTGATCAAGTAATTTGCCACCGGAATTCTTTCGTCAACAATAACTCTGGCGCGCGGCAACGCTATTCGACTGGTATCAGGCGGAGAAATCCTAATAATTGCCGCGTTTGGCGCTTCGTCGCCCCTATCTCTCTCCAAACAAACAACGAGGAGGAGGAAAAAGAGAGAGATAGGGGCGACGAAGCGCCAAACGCGGCAAATTATAGGGATTCATCCGCCTGATAACCAGGTCGAATAGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCAGCGTTGCCCGCGCTCAATTTAATGTTGACGAAGAATTCCCGGTGGCAAATTAACACGTTGAATCAGTTTTATACAAGGTAAAAAATTTATAGCAGTTGCGCAATTAATCCGCCTTACGTCACTTATGAGCAATTCGCATATAAAATGCTAAAACTTTTTGTACTAGCATAAAACACAGAAACGAATACTGGCGACCTGGTCTCGCGGATAAAGCGGTAGATGAGCAAAACAATCACAGCACTGTATTAATGCCTGCCCCACCGCTGCTTCACCGAGTCAGTTTAGGTTTATAGTCTCGGTTATCTTTACCCTTTTCGCTTGAGCTTTCGCAGTTTGGCACCACACTGCCCATGTGGTTCCGCGTCCATCCATGATGGTGGCGTTTTTATCGCCATGCCGGGAGCGCATGTGGCCGGAATTGCGCTGAGCTGTTCGCTGGAAATATCGCCGCATCCGATCCTGCTTTTTTCCACCAGCTCGCTGAACATGACCGGACGACCATCAATATTGTGAAGCCGTGGTCGGGGCAGTGCATCTCTCTCAACAAACAAACAACAGGAGGAGAGAAAAAAGAGAAAGCATGCCCCGACCACGCTTCAACAATATTGAGTGGTCGTCCGGTCATGGTTCAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGCGATATTTCCCAGCGGAACACTCAGCGCAATTCCGGCCAAACAGCGCCCGGCATGGCTGATAAAACGCCACCATCATGATGGACGTGGGAACCACAAGTGGGGCGATTGGTGCCAACTGCGAAAGCTCAAGCGAGAAAGGGTAAAAGATAAACGATGACCTAAACCTAAACTGACCAAGTGAACAGCGGGTGGGGCAGGCAATTAATACATGTGTGATTGTTTGCCATTACCGCTTTATCCCAAGACCAGGTCGCCAGTATTCGATTTCTGTGTTTATGCTAGTACTAAAGTTTTACATTTATATGCAGAATTGCTCATAAAGTGAGTAAGGTCGGATAATTTGCGCACACTGCGTATAACATTTTTTACCTTGTATAAAACGATCCAACGTAATTTGCCACCGGGAATCTCGTCAACAATTAACTCTGAGCGCGCGGCAACGCTATTCGACTGTATCAGGCGGATGAAATCCCTATAATTGGCCGCGTTTGGCGCTTCGTCGCCCCCATCTCTCTCAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATAGGGGGCGACGAAGCGCCAAACGCGGCATATAGGGATTTCATCCGCCTGATACCAGTCCGAAAGCGTTGCCGCGCGCTCAGAGTAAGTTGACGAAGAATTCCCGTGTGGGCACAATAACGTGACAGTTTTATACAAGGTAAAAAATGTTATACAGCAGTTGCGCAAATTTCCCCTTTACGTCACTCACTTTATGAGCAATTCGCATATAAAATGTAAACCCTTTTGTACTAAGCATAAACAGCAGAAACGAATACCTGCGATACCTGGTCTTGCGGATAAACGGTAATGAGCACAACTCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGTCTTCAACCTGGTCATGTTTAGGTTTATCTCGTTTATCTTACCTTTCTCGCTTGAGCTTTCGCATATTGGCACCCAATCGCCCCACGTGGTTCCCGACGTCCATCATGAGGTGGGCGTTTTATCCCATGCCGGGCGGCAATGTGGCCGGAATTGCGCTGAGCTGTTCCGCTGGGAAAATCGCCGCATCCATCTAGCTTTTTTTCCACAGGCTCGCTGAACATGATACCTGGACGACCATCAATATTGTTGAAGCGGGTCGGCAGTGCATCTTCTCAAACAACAACCGTAGGAGGAGGAAAAAAGAGAGAGATGCACTGCCCCGACCAACGGCTTCACAATATTGATGGTCGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGTCGGCGATATTTCCCAGCGAACAGCTCAGCGCAATTCCCGGCCACATGCGCCCGGCATGCGATATAAAGCCACATCATGATGGACGTCGGAACCACATGGAGCGAGTTGGGTGCCAAACTGTGGAAAGCTCAGCGACAAAAGGGTAAAGATAAACGAGACTAAACCTAAACTGACCAGTGAAGCAGCGGGTGGGCAGGGCAATTAATACATGCTGTGATTGGTTTGCTCATACCGCTTACGCAAGACCAGGCCGCAGTATTCTTTCTGTGGTTTATGCTAGTACAAAAGTTTTACATTTTATATGCGAATTGCTCAAAAGTGACGTAAAGCGGATAATTTGCGCAACTCGTATAACATTTTTACCTTGGTATAAAACTGATCAACGTAATTGCCACCGGGAATCTTCGTCACAATTAAGACCTGAGCGCGCGGCAACGTATTCGACTGGTATCCAGCGGATGCAAATCCCTATAATTGCCGCGTTGGCGCTTCGTCGCCCCATTCTCTCAAACAACAACGGAAGGAGAGGAACAAAGAGAGAGCATAGGGGGGCGGGGGCGAACGAAGGCAAACAGCGGCAATTATAGAGGATTTCACTCCGCCTATACCAGTCGAATAGCGTTGCCGCGCGCCAGAGTTAATTGTGACGAAGAATTCCCGTGGCAAATTAACGATGATCAGTTTTATACAGGGTAAAAAATGTTATTACAGCAGTTGCGCAAATTATCCGCCTTACGTCACTTTATGCAGCAATTCGCATATAAAATGTAAAACATTTTGTACTAGCATAAACACAGAAACGAATACTGGCGACCTGGTCTTGCTCGGATAAAGCGGTTAATGAGCAAACAATCACCAGCAGTATTAATTGCCCTGCCCCCACCGCTGCCTTCACCTGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTCTCGCTTTGAGCTTCGCAGTTTGGCACCCAACTCGCCCCCACTGTGGTTCCCGACTCCCATCATGATGGTGGCGTTTTATCGCCATGCGGGCGCATGTTGGCCGGGAATGCGCTGACTGTTCGCTGGGAATATCGCCGCATCATCCTGCTTTTTTCACCAGCTCGCTGAACATGACCTGGACGACATAATAGTTGAAGCGTGTCGGCAGTGCATCTCTCTCAAAACAACAACGAGGAGGAGGAAAAAAGAGAAGATGCCTGCCCCCGCACAGGTTAACAATATGAGGTCGTCCAGGTCGATGTTCAGCGGACTGGTGGGAAAAAAAAAAAAAGCAGGATGGAATGCGCATATTTCCCAACAGCTCAGCGCAATTCCCGGCCCACAATGCGCCCGGCATGGCGATAAAACGCCACCACAGAGGACGTCGGAACCACAGTGGGGCGAGTTGGGTGCCAAATCTCGAAAGCCAGCGAGAAAAGGGTAAAGATAAAACGAGACTAAACCTAAACTGACCAGGTGAAGCAGCCGGGTGGGGCAGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCCGGCTTTATCCCAGACCAGGTCGCCAGTATTCGTTTCTGTGTTTATGCTAGTACAAAAGTTTAACATTTTATAGCGATTGTTCATAAAGTGACGTAAAGGCGGATAATTTGCCAACTGCGTATAACGATTTTTACCTTGTATAAAACTGATCAACGTATTTGCACGGGAATTCTTCGTCCCCAATTGAACTCTGAGCGCGCGCACGCTATTCGACCTGGTATCAGGCGGATGAAATCCGTATAATTGCCCGTTTGCGCTTCGTCGCCCCTATCTCTCTCAAACAACAACGGAAGGAGAGGAAAAAAAGAAGAGATAGGGGGCGACGAAGCGCAAACGCGCAATTTATAGGGATTTCTCCAGCCTGATACCAGTCGACATAGCGTTGCCGCGCGCTCAGAGTTAATGTTGACGAAGAATTCCCGGGTGCAAATTACGTTGATCAGTTTTTTTTTATACAAGTAAAAAATGTTATACGCAGTGGCAAATTTCCGCCATTTACGTCACTTTATGAGCAATTCGCATAAAAAATGTAAAGTTGCATTGTAAAATCTTTTGTACAGCATAAACACAGAAACGAATAGCTGGCGAACCTGGTCTGCGGATAAGAGCGGTAATGAGCAACACACAAGCATGGTATAATTGCCCTGCCCCACCCGCTGCTCACCTGGTCGTTTAGGTTAGTTCGTTTACTTTAACCCTTTCTCGCTTGAGCTTCGCAGTTTTGGCACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGGGCGTTTTACGCCATGCCCGGGCGAGTGGCCGGGAATTGCGCTGAGCTGTTCGCTGGGAATATCGCCGCACCATCCTGCTTTTTTCCACCAGGCTCCGCTGGAACATGCCTGGCGACCAATCATATTGTTGAAGCCGTGGTCGGGGCAAGTGCATCTCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGCATGCCCGACCACGCTTAACAATATTGATGGTCGTCCAGGTCATGTTCAGCGGAGCTGGTGGAAAAAAGCAGGATTGGATGCGGCGATATTTCCCAGCGAACAGCTCAGCGCAATTCCGCCACATGCGCCGGCATGGCGATAAAACGCACCATCAATGTGGACGTCGCTGAACCGAAGTGGGGCGAGTTGGTGCCAAACTGCGGAAAGCTCAAGGAGAAAAGGGTAAAGATAAACGAACAAACTAAATGACCAGGTGAAGCACGGTGGGGCCAGGGAATTAATACAGCTGTGATTGTTTGTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTCATTACCGCTTTATCCGAAGACCAGGTCGCCATATTCGTTTCTGTGGTTTAGCTAGTACAAAGTTTTACATTTTATATGCGAATTGCTCATAAAGTGACGTAAAGGCGAATAATTTGCGCACTGCAGTATAACATTTTTTACTGTATAAAACTGATCAACGTAATTTGCCACCGGGATTCTTCGTCAACAATTTACTCTGAGCGCGCGGCACGCTATTCGACTGGTAATCAGGCGGATGAAATCCCTATAATTGTCCCCCCCGGGGGGCGTTGGCGCTCGCGCCCCCTATCCTTCTCAACAACAACGGAGAGAGGAAAAAGAGAGAGATAGGGGGGGGCGACAAGCGCCAAACGCGGCAATTATAGGGATTTCCATCCGACCTGATACCATCGAATAGCGTTTGCCGCGCGCTCAGAGTGTAATTGTTGACAGAATTCCCGGTGGCAAATTACGTTGATCAGTTTATACAAAGGTAAAAATGTTATACGCAGTTGCGCAAATTATCCGCCTTTAGTCCTTATGAGCAATTCGCATAAAACCTAAAATGTAAAACTTTAGTACTTGCAATAAACACAGAAACAATACTTGGCGACCTGGTCCCTCTGCGGATAAAGCGGTAAATGAGCAAACAACACAGCATGTATTAATTGCCCTGCCCCAACCCGCTGCTTCACTGGTCAGTTTAGGTTTAGTCTGTTTTATCTTTACCCACTTTTCTCGCTGAGCTTCCGCAGTTTGGCACCCAACTCGCCCACTGTGTTCCCGACGTCCATCATGATGGTGGCGTTTATCGCCATGCCGGGCGCATGTGGCCGGAATTTGCGCTGAGCTGTTTGCTGGGGAAATATCGCTGCATCCATCCTGCTTTTTTCCACGAAAAAAACAGCTCGCAGAAATGACCTGGACGACCATCATATTAGTTGAGCCGTGGTCGGGGCATGCATCTCTCTCAAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGCACTGCCCCGAACACGGCTTCAACAATATTGATGGTCGTCCAGGTCCCCATGTTCAAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGTGCGATATTTCCCATCCCCCCCCCGAGCGAACAGCTCAGGCAATTCCCGGTCCACATGCGCCCGGCATGGCGGATAAAACGCCACCCGATCATGATGGACGTTCCAGGGAACCACATGGGGCGAGTTGGTGCCAAACTGCGAAAGCTCAAGGCGAGAAAAAGGGTAAAGATAAACGGACAAACCTAAACTGTCCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGGCTGTATTGTTTGCCATTAGCCCCCCCCCCCCCCCCCATTCCGCTTTTATCCCGCAAGACCAGGTCGCCATATTCGTTCTGTGTTTAATCTAGTACACAAGTTTACATTTATAGCTAATGCGAATTGCTCATAAGTGACGTAAAGGCGGATAAATTCGCAACTGCGTATAACATTTTTTACTTGTATAAAATGATCAACGTAATTTGCCACCGGGAATCTTCGCAACAATAACTCTGAGCGCGCGGCACGCTATTCGACTGGTATCAGGCGGTGAAATCCCTAAATTGCGCGTTTGGCGCTTCGTCGGCCTTCCTCTCAAACAACACGAGGAGGAGGAAAAAAAAAAAAGAGAGAGATAGGGGGCGACGAAGCGCAAACGCGCAATTATAGGGATTTATCCGCCTGATACCAGTCGATAGCAGTTGCCGCGTCGCTCAGAGTTAATTGTTGACAAGAATTCCGCGGTGGCAAATTACGTGATCAGTTTTATACAAGAGTAAAAAAGTTTATACGCAGTTGCGCAAATTATCCGCCTTTACAGTCACTTATGAGCAATTCGCATATAAATGTAAAACTTTTGTACTAGCATAAACACAGAACACGAATACTGGCGACCTGGTCTTGCGGAAAAGCGGTAATGAGCAAACAAATCACAGCATGTATTATTGCCCTGCCCCACCCGCTGCTTCAAACCTGGTCAGTTTAGGTTTAGTCTCGTTTATCTTACCTTTTGCTCGTTGAGCTTTCCAGTTTGGCACCCAACTCGCCCACTGTGGTTCCCGACGTCCATCATGATGGTGGCTGTTTTATCGCCCATGCCGGGCGCATGTGGCCGGAATTGCGCTGAGCTGTTCGCTGGGAAACTATCGCCGCATCCATCCGCTTTTTTCCACCATGCTCGCTGAACATGACCTGGACGCACCATTCATATTGTTGAAGCCGTGGTGGGGCATGCATCTCTCTCAAACAACAACGGAGGAAGAGGAGAAAAAGGAGAGAGATGCACTGCCCGACCACGGCTTCAACAATATTGAGGTCGTCAGGTCATGTTCAGCGAGCTGGTGGAAAAAAGCAGGTGGCAGGATGGGCGGCCGATATTGCCGAGCGAACAGCTCAGCGCAAATTCCGGCACATGCGCCCGGCATGGCGAATAAAACGTGCACCATCATGATGGACAGCGGGAACCACAGTGGGGCGAGTTTTGGGTGGCCAACTGCGAAAGCTCAACGAGAAAAGGGTAAAGATAAACGAGACTAAACTAAACTGACAGGTGGAAGCAGCGGTGGGCAGGCAATTAATACAGCTGTGATTGTTTGCTCATTACCGCTTTTCCGCAAGACCAGGTCGCAAGTATTCGTTTCTGTGTTTATGCTAGTAACAAAAGTTTTACTTTTATATGCGAAAATTGCTCAAAAGTGACAGTAAAAAGGCGGATAATTTTGCGCACATGTCGTATAACATTTTTACCTTGTATAAAACTGACAACGTAAATTTGCCACCAGGGAATTCTTCGTCCAACAATTTAACTCTGAGCGCGCGGCAACGCTATTCGACTGTATCAGGCGGATGAAATCCCTATAATTGCCGGCTGATTGGCTGCTTCGTCCGCCCCCTACGTCCTCAAACAAAACGGAGGAGGAGGAAAAAAAGAGAGATAAGGGGCGACAAGCGCAAAACCGGCAAATTAAGGGATTTTCATCCGCCTGATACCAGTCAGGAATAGCGTTGCCGAGCGCTCAGAGTTAATTGTGACAAGAATTCCGGTGGCACAAATTATTGTCAGTTTATACAAGGTAAAAAATTTATACGCAGTTGCGCAAATTATCCGCCTATTACGTCACCTTTAGACAATTCGCATAAAAATGTAAAACTTTTGTACTAAGCATAAACACAGAAACGAATACTGGCGACTTGGTCTTGCGGATAAAGCGTAATGAGCAAACAATCACACATGTATTAATTGCCCTGCCCCAGCCGCTGCTTCACCAGGTCCAGTTTAGGTTTAGTCTCGTTATCTTTACCTTTTCTGCTTGAGCTTTCGCAGTTTGGCACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGTGCGTTTTATCGCATGCCGGGCGCATGGCCTGGAATTGCGCTGAGCTGTTCAGCTGGGAATATCGCCGCATCCACCTGCTTTTTTCCCACAGCTCGCTGAACAGACCTGGAACGACCATCATATTGTTAGAAAGCCGTGGTCGGGGCCGTGCATCTCTCTCAAACCAACAACGGAGGACAGAAAAAAAGAGAATGCACTGCCCGCCACAGGCTTCAACAATATTGATGGTCGTCCAGGTCATGTTTCAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGCGATATTTCCCCAGCGAACAGCTCAGCGCAATTCCCGGCCAATGCGCCCGGCATGGCGATAAAACGGCCAACCATCATGATGGACGTCGGGACCCATGGGGGCGAGTTGGGTGCCAAACTGCAAGCTCAAGCAGAAGGGTAAAGATAAGATAAACGAGACTAAACCTAAACTTGACCAATGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTAGCTTTAATCCGCAAGACCAGGTCGCCAGTATTCGTTTCTGTGTTTATGCTAGTACAAAAGTTTTATGACATTTTTATATGCGAATTGCTCAAAAGTGACGTAAAGGGGATAAATTTGCGCAACTGGCGTATAACATTTTTTACTTGGTATAAAGCCTAAAACTGATCAACCGTAATTTGCCACCGGGAATTCTTCGTCAACAATTAAACTCTGAGCGCGCGGCAACGCTATTCCGACTGTAATCAGCGATGAAATCCCTATAATTGCCCGCGTTTGCCGCTTCTCGCCCTATCCTCTCAAACAACAACGGAGGAGGGGAAAAAAGAGAGAGATAGGGGCGACGAAGGCCGCCAAACGCGGCAATTATAGGGATTTTCATCCGCCTGATACCAGTCGAATAGCGTTTGCCCGCGCGCTCAGATTAATTGTTGACGAAGAATTCCCGGTTGGCAAATTACTCTGATGCATGCTTTTAGACAAGGTAAAAAGATGTTATACGCAGTTGCGCAAATATCCGCCTTTACGTCACTTTATGAGCAATTCGCATCATAAAATGTAAAACTTTGTACTAGCATAAACAAGAAACGAATACTGGCGAACCTGGTCTGCGGATAAAGCGGTAATGACAAACAATCACAGCATGTATTAATTGCCCTGCCCCACCGCTGCTTCACCGGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTTCTCGCTTGAGCTATTCGCAGTTTGGCACCCTTCTCGCCCACAGAGAAGAGGTGGTTGCCCGACGTCCCATCATGATGGTGGCTTTTATCGCCATGCCGGGCGCATGGGCCGAATAGCGCTGAGCTGTTCGCTGGGGAATATCGCCCGCATCCAATCCTGCTTTTTTCCACGAGCTCGCTGAACATGACCTGGACGACCACAATATTGTTGAAGCCGTGGGTCGGGGGCCAGTGCGATCTCTCTCAACAAAAACGGAGGAGGAGGAGAAAAAAGAGAGGAGCACTGCCCCGACCACGGCTTCACAATATTTGATGGTCGTCCAGGTCATGTCAGCAGAGCTGGTGGAAAAAAGCGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAAATTCCCGGCCAACATGCGCCCGGCAATGCGATAAACCGCCACCAATCATGATGGACGTCGGGGACCACAGTGGGGCGAGTTGGGTGCCAAAACTGCGAAAGCTCAAGCGAGAAAAGGGTAAAGATAAACGAGACAAACCTAAACTGGACCAGGTGAAGCAGCGTGGACAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGCTTTATCCGCAAGACAGGTCGCCAGTAATCGTTTCTTGTGTTTAGCTAGTACAAAAGTTTTTATTTTATATGCGAATTGCTCATAATGTGACTGTAAAAGCGGATAATTTGCGCAACTGCGTATAACATTTTTTACCTTGTATGAAAAACTGATCAACGTAATTTGCCAACCGGGAATCTTCGTCAACAAATAACTCTGAAGCAGCGCGCAACGGCTATTCCGACTGGTATCAGGGCGGCATGAAATCCCTATAATTGCCGCGTTTGGCGCTCGCGCCCCCTATCTCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGCAGATAGGGGCGAGAAGCGCCAAACGCGGCAATTATAGGGATTTCATCCGCCTGATACCAGTGAAATAGCGTTGCGCGCGCTCAGAGTAATTTTGAGAAGAATTCCCGGTGGCAATTACGTGATCAGTTTTATACAAGGTAAAATGTTATACGCAGTTGCCAAATTATCCGCCTTTACGTCACTTTTATGGCAATTCGCATATACAAATGTAAAACTTTTGTACTAAGCATAAACACAGAAACGAATAGCTGGCCGACCTGGTCCTTTGCGGATAAGGGGAAGCGGTAAATGAGCAAACATCACAGCATGTATTAATTGGCCCTGCCCCCCGCTGCTTCACCGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTTCTCTGAGCTTTTCGCAGTTTGGACCCAAATACTCGCCCCACTTGTGGTTCCCGACGTCCATCGATGATGGTGGCGTTTTTATCGCCAATGCCGGGCGCATGTGGCCAGGGAAATTGCGCTGAGCTGTTCGCTGGGAATATCCGCCGCCAATCCTCCTGCTTTTTTTCCACCAGCTCGCTGAACATGACCTCGGACGCCATCAATATTGTTGAAAGCCGTGGTCGGGGCAGTGCATCTCTCTCAAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATGCTCTGCCCACCACGGGCTTCAAAATATGATGGTCCGTCCAGGTCATGTTCAGCGAGCTGGTG";
        // 2048 ~ 2648
        // let seq = b"GCACTGCCCCGACCACGGTTTTTTGTCAACATATTGAAGGTCGTCCAGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAATCCCGGCCACAGCGCCCGGCAATGGCGATTAAAACGCCACCATCATGATGGAACGTCGGGAACCACAGTGGGGCGTGTTGGGTGCAAACTCGAAAAGCTCAAGCAGAAAAGGGTAAAGATAACACGAGACAAACCAAACTGACCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGTTTATCCGCAAAACAGGTCGCCAGTATTCGTTTCCTGTGTTTATGCTAGTACAAAAGTTTTACATTTTAATATGCGAATTGCTCATAAAGTGACGTAAGGCGGAAATTTGCGCAACTGCGTAATAACATTTTTTTAACTTGTATAAAACTGTTCAACGTAATTTGCCACGGAATTCTTCGTCAACAATTAACTACCTGAGCGCGCCAACGCTATTCGACTGGTATCAGGCGAATGAATCCCTATAATTGCCGCGTTTGGCGCTTCGTCCGCCCCC";
        println!("{:?}", aligner.mapopt);
        println!("{:?}", aligner.idxopt);
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
            );
            println!("{:?}", MappingExt(hit).aligned_2_str(targets[0].seq.as_bytes(), seq));
        });

        // aligner.mapopt.best_n = 100;
        // aligner.mapopt.pri_ratio = 0.1; // min dp score
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
            10,
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
            10,
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

    // #[test]
    // fn test_build_aligner_v2() {
    //     let targets = vec![
    //         ReadInfo::new_fa_record(
    //             "t1".to_string(),
    //             "AAAAAGCAGGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAAT".to_string(),
    //         ),
    //         ReadInfo::new_fa_record(
    //             "t2".to_string(),
    //             "CCCGGCCACAGCGCCCGGCAATGGCGATTAAAACGCC".to_string(),
    //         ),
    //     ];
    //     let mut aligner = build_aligner_v2(
    //         "map-ont",
    //         &IndexParams::default(),
    //         &MapParams::default(),
    //         &AlignParams::default(),
    //         &OupParams::default(),
    //         &targets,
    //     );

    //     // let aligner = &mut aligner[0];
    //     aligner.mapopt.best_n = 1;
    //     aligner.mapopt.q_occ_frac = 0.0;

    //     aligner.idxopt.k = 4;
    //     aligner.idxopt.w = 1;

    //     aligner.mapopt.min_cnt = 2;
    //     aligner.mapopt.min_dp_max = 10; // min dp score
    //     aligner.mapopt.min_chain_score = 10; // this is important for short insert
    //     aligner.mapopt.min_ksw_len = 0;

    //     let query = b"TTTCCCAGCGAACAGCTCAGCGCAATCCCGGCCACAGCGCCCGGCAA";
    //     for hit in aligner
    //         .map(query, false, false, None, None, Some(b"query"))
    //         .unwrap()
    //     {
    //         println!("{:?}", hit);
    //     }
    // }

    #[test]
    fn test_align_n() {
        let seq = b"AANNNNNNNNNNNNNAAAAAAAAANNNNCCCGTTT";
        let seq = b"AAAAAAAAA";
        let mut aligner = Aligner::builder()
            .map_ont()
            .with_cigar()
            .with_sam_hit_only()
            .with_sam_out()
            .with_seq_and_id(seq, b"ref")
            .unwrap();
        aligner.idxopt.k = 3;
        aligner.idxopt.w = 1;
        aligner.mapopt.q_occ_frac = 0.;
        aligner.mapopt.occ_dist = 0;

        aligner.mapopt.min_cnt = 2;
        aligner.mapopt.min_dp_max = 2; // min dp score
        aligner.mapopt.min_chain_score = 2; // this is important for short insert
        aligner.mapopt.min_ksw_len = 0;
        aligner.mapopt.mid_occ_frac = 0.;
        aligner.mapopt.min_mid_occ = 0;

        println!("aligner.mapopt:{:?}", aligner.mapopt);
        println!("--------------------");

        println!("aligner.idxopt:{:?}", aligner.idxopt);

        // b"AACGTCGTCGTCGTAAAAAAAAACGTGCCCGTTT",

        for hit in aligner
            .map(
                b"AAAAAAAAA",
                false,
                false,
                None,
                Some(&[67108864, 68719476736]),
                Some(b"q"),
            )
            .unwrap()
        {
            println!("hit:{hit:?}");
        }
        println!("hello");
    }

    #[test]
    fn test_bio_aign_n() {
        let scoring = Scoring {
            gap_open: -2,
            gap_extend: -1,
            match_fn: |a: u8, b: u8| {
                if a == 'N' as u8 {
                    0
                } else if a == b {
                    1i32
                } else {
                    -3i32
                }
            },
            match_scores: Some((1, -3)),
            xclip_prefix: 0,
            xclip_suffix: 0,
            yclip_prefix: 0,
            yclip_suffix: 0,
        };
        let x = b"NNNNAAAAAANNNNAAAAAAA";
        let y = b"ACGTAAAAAAGGACGGAAAAAAA";
        let mut aligner =
            bio::alignment::pairwise::Aligner::with_capacity_and_scoring(x.len(), y.len(), scoring);
        let alignment = aligner.custom(x, y);
        println!("{}", alignment.pretty(x, y, 80));
        println!("{alignment:?}");
    }

    #[test]
    fn test_bio() {
        let score = |a: u8, b: u8| if a == b { 2i32 } else { -4i32 };
        let gap_open = -4;
        let gap_extension = -2;

        let mut aligner = bio::alignment::pairwise::Aligner::new(gap_open, gap_extension, &score);

        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let alignment = aligner.local(x, y);
        println!("{:?}", alignment);
        println!("{}", alignment.pretty(x, y, 30));

        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let alignment = aligner.semiglobal(x, y);
        println!("{alignment:?}");
        println!("{}", alignment.pretty(x, y, 30));
    }

    #[test]
    fn test_align_short() {
        let seq = b"GGTAGCGTTACAAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATGCTCTGCCCACCACGGGCTTCAAAATATGATGGTCCGTCCAGGTCATGTTCAGCGAGCTGGTG";
        let qry = b"GGTAGCGTTTCAAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATGCTCTGCCCACCACGGGCTTCAAAATATGATGGTCCGTCCAGGTCATGTTCAGCGAGCTGGTG";
        let mut aligner = Aligner::builder()
            .map_ont()
            .with_cigar()
            .with_sam_hit_only()
            .with_sam_out()
            .with_seq_and_id(seq, b"ref")
            .unwrap();
        aligner.idxopt.k = 3;
        aligner.idxopt.w = 1;
        aligner.mapopt.q_occ_frac = 0.;
        aligner.mapopt.occ_dist = 0;
        aligner.mapopt.a = 2;
        aligner.mapopt.b = 5;
        aligner.mapopt.q = 2;
        aligner.mapopt.q2 = 24;
        aligner.mapopt.e = 1;
        aligner.mapopt.e2 = 0;

        aligner.mapopt.min_cnt = 2;
        aligner.mapopt.min_dp_max = 1; // min dp score
        aligner.mapopt.min_chain_score = 1; // this is important for short insert
        aligner.mapopt.min_ksw_len = 0;
        aligner.mapopt.mid_occ_frac = 0.;
        aligner.mapopt.min_mid_occ = 0;

        println!("aligner.mapopt:{:?}", aligner.mapopt);
        println!("--------------------");

        println!("aligner.idxopt:{:?}", aligner.idxopt);

        // b"AACGTCGTCGTCGTAAAAAAAAACGTGCCCGTTT",

        for hit in aligner
            .map(
                qry,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]),
                Some(b"q"),
            )
            .unwrap()
        {
            let mapping_ext = MappingExt(&hit);
            let aligned = mapping_ext.aligned_2_str(seq, qry);
            println!("{}\n{}", aligned.0, aligned.1);
        }
        println!("hello");
    }

    #[test]
    fn test_align_short2() {
        let target = b"AAACAACAACGGAGGAGGAGGAAAAAA";
        let qry = b"AAAAAAAACCCCCGGGGTTTTTGAGAGAGATATCTCTCTCAAAACCCCGGGGTTTTAAACAACCAACGAAGGAGGAGGAAAAAAAAAAAAAACCCCGGGGTTTTTGAGAAGATATCCTCTCAAAACCCCGAGGGTTTTAAACAACAACGGAGGAGGAGGTAAAAAAAAAAAAAAAACCCCGGGTTTTGGAGAGAAATATCTCTCTCAAACCCCCGGGGTTTTAAACAACAACGGAGGAGGAGGAAAAAAAAAAAGAAAAAACCCACGGGTTTGATGAGAGATATCCTCTCAACCCCGGGTTTTTAACAAACACGGAGGAGGAGAAAAAAAAAAAAAAAAAAAAACCCCGGGGTTTTGGAGAGAGAATTTCTTCTCTCACAACCCCGGGGTTTTAAACACAACGAGGAGGAGGAAAAAAAAAAAACCCCGGGGTTTTGGAGAGAGAATATCTCTCTCAAACCCCGGGGTTTTTAAACAACAACGGAGGAGGAGGAAAAAAAAAACCCCGGGGTTTTGAAGAGATATCTCTCTCAAAAACCCCGGGTTTTAAAACAACAACGGAGGGGAGGAAACAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAGAAAAAAAAAAAACCCGGGAGTTTTTGAGAGATATCTCTCTCAAAAACCCCGGGGTTTTAAACAACAACGGAGGGAGGAGGTGAAGAGGAGGAGGAAAAAAAGGAGGAGAAAAAAAAAAAAAACCCCGGGGTTTTGAGAGAGATATCTCTCTCCAAAACCCTGGTTTTTAAACAACAAACGGAGGAGGAGGAAAAAAAAAACCCGGGGGTTTTGAGAGCAGATATCTCTCTCAAAACCCCGGGGTTTAAACAAAAACAACAACGGAAGGAGGGGAAAAAAAAAAAAAAAAAAAAAGAGGAAAAAAAAAAGCAAAAAAAAACCCCGGGGTTTTGAGAGAGAGAAGATTCTCTCTCAAAACCCCGGGGTTTTAAACAACTACGAGGAGGAGGAAAAAAAAAAAATCCCAAAAAAAAACACCGGGGTTTTGAGAGAGATATCTCTCTCAAACCCCGGGGTTTTAAACAACAGACGGAGGAGGAGGAAAAAAAAAAATAAAAAAAAAAAGGGCAAAAAAAACCCCGGGGTTTTGAGAGAGAATATCTCTCTCAAAACGCCCGGGTTTTAAAACAACAACAGGAGGAGGGAAAAAGAAAAGGAGGAAAAAAAAAACCCCCCGGGGTTTTGAGAAATATCTCCTCAAAACCCCGGGGTTTTAAACAACAACGAGGAGGAGGAAAAAAAAAACCCCGGGGTTTTGAGAGGATATCTCTCTCAAAATATCAAAACCCGGGGTTTAAACAACACGGATAGGACGAGGAAAAAAAAAACCCTCGGGGTTTTGAGAGAGATATCTCTCTCAAACCCCGGGGTTTTAAACAACAACGGAGGAGGAGGGAAAAAGTAAAAAAACCCCGGGGTTTTGAGAGAGAAGAAGAGATATCACTCTCAAAAACCCCGGAGGTTTTACCAACAACAA";

        let mut aligner = Aligner::builder()
            .sr()
            .with_cigar()
            .with_sam_hit_only()
            .with_sam_out();
        aligner.idxopt.k = 5;
        aligner.idxopt.w = 3;
        let mut aligner = aligner.with_seq_and_id(target, b"ref").unwrap();

        // aligner.mapopt.q_occ_frac = 0.;
        // aligner.mapopt.occ_dist = 0;

        aligner.mapopt.min_cnt = 3;
        aligner.mapopt.min_dp_max = 30; // min dp score
        aligner.mapopt.min_chain_score = 5; // this is important for short insert
        aligner.mapopt.min_ksw_len = 0;
        aligner.mapopt.mid_occ_frac = 0.;
        aligner.mapopt.min_mid_occ = 0;

        println!("aligner.mapopt:{:?}", aligner.mapopt);
        println!("--------------------");

        println!("aligner.idxopt:{:?}", aligner.idxopt);

        // b"AACGTCGTCGTCGTAAAAAAAAACGTGCCCGTTT",

        for hit in aligner
            .map(
                qry,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]),
                Some(b"q"),
            )
            .unwrap()
        {
            let mapping_ext = MappingExt(&hit);
            let aligned = mapping_ext.aligned_2_str(target, qry);
            println!("{}\n{}", aligned.0, aligned.1);
            println!("hit:{hit:?}");
        }
        println!("hello");
    }

    #[test]
    fn test_align() {
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let target = b"ACAAAAAAATTAAAATAAAGAATTCCCGGCTTGGGGGCCAGAGTCCTCACCTCAT
ATGTAATAAGATGTGTGGATCTGTGACAGTTTCACATCCTCCTGCCTCACTGTCTCCATATAAAATAAAACCGTTTCCTGGATGACAGAGCATGAGTCATAAACAGATGGGAAAGCCATGCAAAAGT
AAGAGGACTCTTCGTTTTTATGACTGACGGATCCTCCGCAAAAATACTTGCACGTCTCTCATGTGTGTCAGGTTCTCTATTGCAGGCCTGAATGCATACGATTAACTTGCTTATCATTTGGAAGTAAAGAGAAACACTTTCTGTCTTAGAACTTTTTGTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAGAGAGAGAGAAGAGAGAGAGAGAGAGAGAGAGAGAAAGAGAGAGAGAAAGAGAGAGAGAGAGAGAGAGAGAGGAGAGAGAGAGAGAGAGATTTTGCTATTCTTTGCATAGCATTTGGTTATGGCCTATTATAAGACTCATTATAGTCAGGGAAATTATGGGGTGAGAGAAACACTCTTGTGGAAAGATGAGCTGAGGTTTATATAAAGAAATTGACTAGTTTATGGAGGCTTATAAGTTTTTAGAAGCTACTCCCATTTGAGACTGACTTAGCAGTGATTTTAAACTTTAACATCAAGAAGGAAGAGGAGGAGAAGGAAGAGGAAGAAGAGGAAGAAGAGAAGAGGAAAAAGAGGAAGAAGAGGAAGAGGAAGAAGAGGAAAAAGGAGAAAGACAGACACCAAATAATCAACTGTATTTCTGTTTGTATCCCTTGACTGTCTTCCAAGTTTCTGCTTCACTTTTGTTTCCATAGCAACATGCAGAGATAAAACAGGTAGTAAATATTAAGCAAAAGTTGAATGAATAAAAGGAGAATAAAAAGCCAAGGAAAATGCTTAAAAATCTTACCTGTGGCTCAAGCATCCAACTACCAGTTCCTAGAAGATAGCTTACAGAAGATGGTAACTGTATCCATCAGGATACTGTCAGCAAAATCCTAGGGAAAAATGTTTTTAGACAAACCAACTACATTTCTTGAACAAAAGTAATTGCAAGTGGACAAGAGCTAGAGGAGGATTGGAGAACAAAGGAGACACAGGAGACATAGACATCAGTTGCAAAGAGTACACTTTATCTGGATCCTGGTTTGTTGTTGTTGTTTGTTTGTTTTGGTTTTTTGTATGTTTTGATTATTTTATTTTATGTGCATTGGTGTTTTGCCTGTGTGTATGTATATGAACCATGGGATGACTGGTGCTGGCAAAAGTCAGAAGAGGGTCTCATCTCCCCCAGAACTGAGGTAACAGATGGTTGTGGGACCCTAACCTGGGTGCTCTGCAAGAGCAGAAAGTGCTCTCAACTGCTGAGCCAATCTATTCTCTCAGCCCCAGAATTCTGATTTTTTAAAAAAATTATGATATCATGAGACATTTAGAAATCTGAACATTGCCAAATTTGGTGGCATGAGTTATCAATTGTTTTGCTTGTAATGTTTACAAATAAAATGATATCATATCTAGACTTTGCCTCAAAGTCATACGGGAGGAGGACTTGGGCAACAATAGGATGGGCTGAGGTAGGGACTGAGTCTGGAGAGGAGTTCCATGCTCACACATAAGGCATCACCACAGTGGGCAAGTTAATGTGTGCTTTGAGCTTCTTACTGAAGTGTGGCAAACAGAAAGCCCTGAAAAGCCCCTGGTGCACAACAGGAACCATTATACAAATAAAAAGATGGGTTCTCCCTCAGTGGGAAAAGGTTTTCGTTTTACATAAGGATATGTGTGATTTCTTACCTGTGCTCAAATTACCAAGCAATCAGGAACTTACAATCTTTTTCAGAGTAGCTACCCTCTGAGGAGATGTATTGTTTTTTAACACTAGTTGTTGCAGGCAGTCTTTCAATGCAAAGCAAACTGGGAAAATGGTCTCCTGGCACAATCGCCAACATTAGTAATAGGTAGTGAAGCTTGTTCACGTCTACAAAGTTCAAAGTCCTACAAGTTAAGCTGAGTGTCCTAAGGTGCTTATCACTGCAGTTAAACAAGGCGTGGTAACGGAATGTGCACAGCACTCTGGGAAAAATTAGTTATGTCCGTTTGCACAATGTTGCTTACTTAGAACGTTGACTTGCCTCAACCTCACACTGTGTCCAAGCTTTTAACTGACAACTGTTTTCCCAAAACCCAAAGATTATGTGTAAAAGAACAATCCTTACTTGTATCTGACTTGCCCTCAGAAACCCAGGTAAGGCTGTGTTCTCCAGCCCAAGCCCTCCCGTTGTAGACTCTGAATGGTTGGCACAAAGCATGGGTTTAAATGACCAAATGGTGTTTTTATGTTACAGTTTAGTTTGTAAGTTATTTGGAAAGCACTTCTTACCAGTAACCAGCCGCTTTGCCAGCTAACAGCTAAGTGCTTTTCAGAGATTCCTAATGTCTGTCAAGAAGAGTGTGGCTGCTCTGCTTCACTAATGAGGAAAGCAAGGCTCAGAGAAATGGGTGGCTCACCAGCCTGTAATGGCAGGCAGAGCAGGCAGTCAGGAGTGGGCTGTCTGGATTTATTGTTTTCTCCTAGTGTGCACGTGTGCATACGTGTGTGTGTGTGTGTGCATACGTGTGTGTGTGTGTGTGTGTGTTCATGGGTGGATAGATACAAAATGTGTGCATGTGAAGAGCAAAACTTCAGATGCCATTCCTTAGTCACTATCCTAACTTTTTCTCTCTCTTGAAGTAGGGTTTCTCACTGGACTGGGCAGTTAGCCACTACAGCTGGACCAACTGACTATGGAACCCACACAACCCTCTATACATCATCCCACACACAACTGGCATTACACCTTGCCATCATGCTATCATGCTGTTCCCCCTCCCTATTTCATATAAATACATCCTTCCTTCCCTCCCTAATACTTCCTCCTCCTTTCACATTTGTTCTGGGTACCAGACTCAGGTCCCCACACTTGCAAGGCAAGCACTTTCTGGCTGAGCCATCTCTCCAGCCATTATCCTATAGCATTAATCATTTTCTCATTGACATTGTAAGCCCCATCAGAGAACTCGTTCATTAGTCTCCAGGGCTCCTTAACCCTAGTAGAGTTATATAAACAATGTCGCTCTCAACAAATCGTTCAAGGGCTGATGTCCTGCAACTTGTCTAACCTTCTCAAGCCATTGCTGCTTACCAGTGTTTTCCAGGCTCATATTTTGAAACTGACCAAGTGTGATTATCCAGATAACTATTTAAGGGAAATGTTTGAGAACAACGTCTGCATCTGGAGCGAGAGAATTACACTTAACAAGGGAAAGAGGAATGAAAACCAAACACTGTAAGCAAGGTGAGGGAAACTGAGATTGAAATATAAGGACCAATCGGTCAGGTTGATGCAGTTTCCCTGTGTAGGCCCTAGAGTAAGTGGAAAACCCAACCAACAATGGAGAAAATTCTAGAACACTGGAACCAAGTGGAACACACACACACACACGGGATGGTTGTAGAAGATAGATACTCTAGATATAAAAAGGCAGCCACACTCTCCCTGTTCCCCTGTCCCCAGCCTGTTGGGACCTGTTTGCCTTTCCACTGCTTAGAGGAACATTTTTACTGGTGTGGTCTTATCCGGATTGAGAGTGATGAGGTTCCAGTCAAACAATTGGACAAGCACGGTGCTAAACCATAGCTGATGGCTATCCCATGCCCACGAACGTTATATGTGTTATCTAATTTAATTCTGATACACACAAGCAAAGGACGGAGATGCAGAGACTGTGGCTTTGCAGTTACAAGGTACTTATTAAATTCTAGGCTCAGGTTTAGCACTCTTGCTGGGTTTTGCTAGGGGAGGATGGCCAGTAGGGTTGGGAGAGTGCTACACCACATCTTTCCCCCAGCACTCGGGTGATAAAGGTGGCAAGTGTGTAGAGTGTTCTGCAGGAGCCAGCTCTAGTCATGTCTAGAGAAAACATTTGGGGATTGCTGCTTTGTCTACCATCTGATTCTCTGTGTCTGGATTCTGTAGGACTCGCCAGCTCATTTATCCTGGTTCCAGTTTATGCACCCAGCCTGGTTGACTTTTCTGGAACATCGAGGTTACTGTAAACTCTAAGGGCAAGGGCCTAGGCTCCTGGCAAGTTTCCAGTATGAGATCTGCTAAGTACAAACTATTGGTGGGTTCATGATTATGGTCTCTCACTGAGAGAGCGTGATTATACACTACTATATCATAACATTCAGGAAGATGACCAGGGCATTGGGAGTTCCAGGGCAGCCTGGGCTACATAGCCAACTCTTGTCAAAAGAGAAAGATTTGAAAATTATAAATCAATTTGAAATTATTTTGTTACTAGTAGAAGTGAACTAAATTCTGGAATCTTTAAGCATGTATAATTTTAAAATCTATGCTTTATGTGTTCTTACCAAATGCAATTTTTTTCCTGTAGGTTTTCTTAAAAGGTGTTGGCATTTTTGTAAGAAGTGATTGTTAATAATACCTGGAGTTTGACTTGACTTTTTAAATAGCATACTTGTCTTCCTAAAACAAAACTGGTCCATCTTTTGTTAAAAAGAATAAAGTTAATTTCTGAAGGTTGTCCAATAATGTATCTACTTGTTAATTTTAAAAAGGGTGTCTTGACTCAGGCTGCTATTACAGAATTCCATGTTAGATGTTTTGTATATAATAGAAATTTATATCTCATACTTCTGGAGTCCAACCTACTAAAATGTTCACTTTATACTTTATCAAGTTAGCAGTGAATATATTTTTAACAGCTGATGTTCCCCCATGCCCCCAGCTTGCTTCCACACCTTGGATATTTTGTGGGCGAGAAGTCTGAAGCCAACATTTCAGAAGGGCTCCCTATTTAGCACAATGTCTTTACACAAGCATACAGGTGCTCTCATATGATAATATTAGCAGGCCTGGCAGGAATAGCAAAGACCAATCTTTATATAGTATTGTTTAGGCATTAAGAACTGTTTTAGGTATTTTCTGTGTGTTAGACCATTTAATCCTTTCAACAGCTCAAGTAAAGTATGAAACGGTCATCCTCAGCACTTCATTACAGAGGAACCGAGGAACGGAGCCTGAGAGGGTTTGCCTGAGACTGTGTAGCTAGTGGAAGAAAGATGAGCATTCTCAATGGTGATTACCAACTTTTGTCCCCTTTGAATGAAATTCGCCTTTAAACTATGGTGCTCAGGTCGCACCTGTAATTGGCCTGGAACATTTTGTTTTATTTTATTTTATTTCATTTTATTTTATTTATTTATAGTCTCAATGTCATGATATCCTGCCTCTGCCCCTGAGTGTCTGTCAGGCTATAGGTATAGAACACTATATCCACATTGAAAAGGTTTTTAATGCTCCTTCGGTGCGCTGAATGACTGTTGGGAATCGGATGCCTTAGCTTTTCTACCACAGACCTGTAAGTGGCTCTTGAGCTGACACCACACCTAGTTAAAGGAACAGACAACAGAAGGGGGATTCTTGGCAATAAAATGGTAAAAGCTGAATGGACATTGAGCTGTCTACAGGTTCTTGGGCAACTGACTCAGTTGCCATTTCTCCATCTGCAAAAACAGAGTTTAAAGTAAAAGATCTGCTGGGCTTTCTAGCTCTGTGAGTCCATGAACTTACTTCTATAAAATTCTCTGTAATCCGTGAGTAGATGACCCTTCAGGAGTTGTTTCATGTGGGAACCAGTTTCCTTTGGCCAACCTTCTTATCAACTGGTGTGTGCCTAAGGAAGGCAGATTGGTTTTTACCTTTGCTTTGGGTGGAACTCAATGAGGGAGGAGGAAACTGTGTGCGCGCACACAAATCATCCTGAAGGTAAGGGGACTGATTTGGAGGGATTATCAGCTGGTTCAGATTCCCCACGTCACTGCCCCCTCCAATGGCTGCGGTAATTGAAAGCTGACCTGCGTTCTGATTGCTGACAGCTGAAGGATTCCTGTTGCCAGGGAACCTCTGGCCGTGTGCCAAAGTAACGACCTGATTAATTATAATAATCCTTGTGCTGTTCCTGCCCCCATCGGCTCCCTCCTCTGAGTGGCGTCCAGCCTCAGGATCTCAGAATGGAACAGTAAGCACCCACTAAATGAACAATCCTGTGTCCTGTGGCAGAACACATGTCAAAAGACTTTCGCACCCCTCCGACCGCTCTGCTCAGTGACATACCTGTTGAAAACAATCATAAATTCAAAGTGGGTCTGACCAAGCCACACCGTGTTTCCACATAGAAGGGTCCAGAAACCTGGCCACACCGTGTTTCCACATAGAAGGGGTCCAGAAACCTGGCCACACCGTGTTTCCATATAGAAGGGTCCAGAAACCTGGCACACCGTGTTTCCACATAGAAGGGTCCAGAACCTGGCCACACCGTGTTTCCACATAGAAGGGTCCCAGAAACCTGGCCACACCGTGTTTCCACATAGAAGGGTCCAGAAACCTGGCCACACCGTGTTTCCACATAGAAGGGTCCAGAAACCTGTTCTTCTCTCCTCCACTCTGACGCATCTTGAGATTACCCCAAATCTGGAGACTCCACGATAAAAACAACATATTCTGCCATGTCCCAAAAGCTTCTGGTAACTTGGGACTGACTGACTTAAACTAGTAATCGTGCTTGGGATTTTAGTTTGGAGGGTTACCTAGAATTCTGACACCCTTCTGAAATTCAGGTACCAAGTGACAAGATTTCTGGACCCTTCTATGTGGAAACACGGTGTGGCCAGGTTCTGGACCCTTCTATGTGGAAACACGGTGTGGCCAGGTTTCTGGACCCTTCTATGTGGAAACACGGTGTGGCCAGGTTTCTGGACCCTTCTATATGGAAACACGGTGTGGCCAGGTTTCTGGACCCTTCTATGTGGAAACACGGTGTGGCCAGGTTTCTGGACCCTTCTATGTGGAAACACGGGTGTGGCTGGTCAGACCCAATTTGAATTTATGATTGTTTTTAACAGGTATGTCACTGAGCAGAGCGGTCGGAGGGGGTGCGAAAGTCTTTTGACATGTGTTCCAACAGGACACAGGGATTGTTCATTTAGTGGGTGCTTACTGTTCCATTCTGAGATCCTGAGGCTGGACGCACTCAGAGGAGGGAGCCGATGGGGCAGGAACAGCACAAAGGAGGGATAATTAATCAGGTCGTTACTTGGCACACGGCCAGAGGTTCCCTGGCAACAGGAATCTTAGCTGTCAGCATCAGAACGCAGGCAGCTTTCAATTACCGCAGCCATTGGAGGGGGCAGTGACGTGGGGAATCTGAACCAGCTGATAATCCCTCCAAATCAGTCCCCTTACCTTCAGGATGATTTGTGTGCGCGCACACAGTTTCCTCCTCCCTCATTGAGTTCCACCCAAAGCAAAGTAAAAACCAATCTGCCTTCCTTAGGCACACACCAGTTGATAAGAAGGTTGGCCAAAGGAAACTGGTTCCCACATGAAAAAATCTGAAGGGTCATCTACTCACGGATTACAGAGAATTTTATAGAAGTAAGTTCATGGACTCACAGAGCTAGAAAGCGCAGCAGATCTTTTACTTAAACTCTGTTTTTGCAGATGGAGAAATGGCAACTGAGTCAGTTGCCCAGAACCAAATTGTAGACAGCTCAATGTCCATTCAGCTTTTACCATTTTATTGCCAAGAATCCCCCTTCTGTTGTCTGTTCCTTTAACTAGGTGTGTGTCAGCTCAAGAGCCACTACAGGTCTGTGGTAGAAAAGCTAAGGCATCCGATTCCCAACAGTCATTCAGCGCACCGAAGGAGCATTAAAAACCTTTTCAATGTGGATATAGTGTTCTATACCTATAGCTTGACAGACACTCAGGGGCAGAGGCAGGAGGATCATGACATTGAGACTATAAATAAAAAATAAAATAAATAAAATAGAAATAAACAAAATGTTCCAGGCCAATTACAGGTGCGACCTGAGCACCAGAGTTTAAGGCGAATTTCATTCAAAGGGGACAAAGTTGGGAATCACCATTGAGAATGCTCATCTTTTCTTCCACTAGCTACACAGTCTCAGCAAACCATCTCAGGCTCCGTTCCTCGTTCTCTGTAATGAAGTGCTGAGGATGACCGTTTCATACTTTACTTGAGCTGTTGAAAGGATTAAATGGTCTAACACACAGAAAATACCTAAAACAGTTCTTAATGCCTAAACAATACTATATAAAGATGGTCTGTTGCTATTCCTGCCAGGCCTGCTAATATTATCATATGAGAGCACCTGGTATGCGTGTGTAAAGACATTGTGCTAAATAGGGAGCCCTTCTGAAATGTTGGCTTCAGACTTCTCGCCCAAAAAATATCCAAGGTGTGGAAGCAAGCTGGGGGCAGGGGGAACATCAGCTGTTAAAAATAGATTCACTGCTAACTTGATAAAGGAGAAAGGGAACATTTTAGTAGGTTGGACTCCAGAAGTATGAGAGATAAATTTCTATTATATACAAAACATCTAACATGGAATTCTGTAATAGCAGCCTGAGTCAAGACACCCGTTTTTAAAATTAACAAGTAGATACATTATTGGACAACCTTCAGAAATTAACTTTATTCTTTTTAACAAAAGAGGGACCAGTTTTGTTTTAGGAAGACAAGTATGCTATTTAAAAAGTCAAGTCAAACTCCAGGGATTATTAACAAGCACTTCTTACAAAAATGCCAACACCTTTTAGAAAACCTACAGGAAAAAAATGTGCATTTGGTAAGAACACATAAAGCATAGATTTTAAAATTATACATGATTAAAGATTCCAGAATTTAGTTCACTTCTACTAGTAACAAAATAATTTCAAATTGATTTATAATTTTCAAATCTTTCTCTTTTTGACAGAGTGGCTATGTAGCCAGGCTGCCCTGGAACTCCCCATGCCCTGGTCATCTTCCTGAATTTAGGAATATAGAGTGTATAATCACGCTCTCTCAGTGAGAGACCATAATCATGAACCCACCAATGAGTTTGTACTTAGCAGATCTCCATACTGGAAACTTGCCAGGAGCCTAGGCCCTTGCCCTTAGAGTTTACAGTAACCTCGATGTTCCAGAAAAGTCAACCAGGCTGGGTGTAAATATACTGGAACCAGGATAAATGAGCTGGCGAGTCCTACAGAATCCAGACACAGAGAAGCAGATGGTAGACAAAGCAGCAATCCCCAAATGTTTTCTCTAGACATGACTAGAGCTGGCTCTGCAGAACACTCTACACACTTGCCACCTTTATCACCCGAGTGCTGGGGGAAAGATGTGTGGTAGCACTCTCCCAACCCTACTGGCCATCCTCCCTAGCAAAACCCAGCAAGAGTGCTAAACCTAGCCTAGAATTTAATAAGTAACCTTGTAACTGCAAAGCCACAGTCTCTGCATCTCCGTCCCTTTGCTTGTGTGTATCAGAATTAATTAGATAACACATATAACGTTCGTGGCATGGGATAGCCATCAGCTATGGTTTAGCACCGTGCTTGTCCAATGTTTGACTGGAACCTCATCACTCTCAATCCGGAATAAGACCACACCAGGAAAATGTTTCTTCTAAGGCAGTGGAAAGGCAAACAGGTCCCAACAGGCTGGGACAGGGGAACAGGGAGAGGTGTGGCTGCCTTTTTATATCTAGAGAGTCTATCTTCTACAACCATCCCGTGTGTGTGTGTGTGTTCCACTTGGTTCCAGTGTTCTAGAATTTTCTCCATTGTTGGTTGGGTTTTCCACTTACTTCTAGGGCCTACACAGGGAAACTGCATCAACCTGACCAGATTGGTCCTTATATTTCATCTCAGTTTCCCTCACCTTGTCTTACAGTGTTTGGTTTCATTCCTCTTTACCTTGTTAAGTGTAATTCTCTGCCTCCAGATGCAGACGTGTTCTCAAACATTTCCTTAAATAGTTATCGGATAATCACACTTGGTCAGTTTCAAAATATGAGCTGAAAACACTGGTAAGCAGCATGGCTTGAGAAGGTTAGACAAGTTGCAGGACATCAGCCTTGAACGATTTGTTGAGAGCGACAGTGTTTATATAACTCTACTAGGGTTAAGGAGCCCTGGAGACTAATGAACGAGTTCTCTGATGGGGCTTACAATGTCAATGAGAAAATGATTAATGCTATAGGATAATGGCTGGGAGAGATGGCTCAGCCAGAAAAGTGCTTGCCTTGCAAGTGTGGGGACCTGAGTCTGGTACCCAGAACAAATGGTGAAAGGAAGGAGGAAGGAAGGGAGGGAGGGAGGGAGGGAAGGAAGGAGGGAGGAGGGAGGAAATAGGGAGGGGGACAGCATGATAGCATGATGGCAAGGTGTAATGCCAGTTGTGTGTGGGATGATGTATAGAGGGTTGTGTGGGTTCCATAGTCAGTTGGTCCAGCTGTAGTGGCTAACTGCCAGTCCAGTGAGAAACCCTACTTCAAGAGAGAGAAAAGGTGGATAGTGACTAAGGAATGGCAATCTGAAGTTTTTGCATCTTCACATGCACACATTTGTATCTATCCACCCATGAACACACACACACACACACACGTATGCACACACAACACACACGTATGCACACGTGCACACTAGGAGAAAACAATAAATTCCAGACAGCCCACTCCTGACTGCCTGCTCTGCCTGCCACTTACAGGCTGGTGAGCCACCCATTTCTCTGAGCCTTGCTTTCCTCATTAGTGAAGCAAGAGCAGCCACACTCTTCTTGACAGACATTAGGAATCTCTGAAAAGCACTTAGCTGTTAGATGGCAAGCGGCTGGTTACTGGTAAGAAGTGCTTTCCAAATAACTTACAAACTAAACTGTAACATAAAACACCATTTGGTCATTTAAACCCATGCTTTGTGCCAACCAATTCAGAGTCTACAACGGAGGGCTTGGGCGAGAGAACACAGCCTTACCTGGGTTTCTGAGGGCAAGTCAGATACAAGTAAGGATTGTTCTTTTACACATAACTTTGGGTTTTGGGAAAACAGTTGTCAGTTAAAAGATTGGACACAGTGTGAGGTTGAGGCAAGTCAACGTTCTAAGTAAGCAACATTGTGCAAACGGAACATAACTAATTTTTCCCGAGTGCTTGCACATTCCGTTACCCACGCCTTGTTTTAACTGCAGTGATAAGCACCTTAGGACACTCAGCTTAACTTGTAGGACTTTGAACTTTGTAGACGTGAACAAGCTTCACTACCTATTACTAATGTTGGCGATTGTGCCAGGAGACCATTTTCCCAGTTGCTTTCATTGAAAGATGCCTGCAACAACTAGTGTTAAAAACAATACATCTCCTCAGAGGGTAGGCTACTCTGAAAAAGATGTAAGTTCCTGATTGGCTTGGTAATTTGAGCACAGGTAAGAAATCACACATATCCTTATGTAAAACGAAAACCTTTTCCCACTGAGGGGAGAACCCATCTTTTATTTGTATAATGGTTCTGTTGTGCACCAGGGGCTTTTTCAGGGCTTTCTGTTTGCCACACTTCAGTAAGAAGCTCAAAGCACACATTAACTTGCCCACTGTGGTGATGCCTTATGTGTGAGCATGGAACTCCTCTCCAGACTCAGTCCCTACCTCAGCCATCTATTGTTGCCCAAGTCCTATCCCGTATGACTTTGAGGCAAAGTCTAGATATGATATCATTTTATTTGTAAACATTACAGCAAACAATTGATAACTCATGCCACCAAATTGGCAATGTTCAGATTTATAAATGTCTCATGATATCATAATTTTTTTAAAAAATCAGAATTCTGGGGCTGGAGAGATAGATTGGCTCAGCAGTTGAGAGCACTTTCTGCTCTTGCAGAGCACCAGGTTAGGGTCCCACAACCATCTGTTACCTCAGTTCTGGGGGAGATGAGACCCTCTTCTGACTTTTGCCAGCACCAGGCATCCCATGGTTCATATACATACACACAGGCAAAACACCAATGCACATAAAAATAAATAATCAAAACATACAAAAAACCAAAACAAACAAACAACACAACAAACCAGGATCCAGATAAAGTGTACTCTTTGCAACTGATGTATATGTCTCCTGTGTATACTTTGTTCTCCAATCCTCCTCTAGCTCTTGTCCACTTGCAATTACTTTGTTCAAGAAATGTAGTTGCTTGTCTAAAAACATTTTCCCTAGATTTTGCTGACAGTATCCTGATGGATACAGTTAACATCTTCTGTAGCTATCTTCTAGGAATGGTAGTTGGATGCTTGAGCCACAGGTAAGATTTTAAGCATTTCCTTGGCTTTTTATTATAATTTTTTCATTCAACTTTTGCTTAATATTTACTACATGTTTATCTCTGCATGTTGCTATGGAAACAAAAGTGAAGCAGAAACTTGGAAGACAGTCAAAGGGATACAAACAGAAATACAGTTGATGTATTTGGTGTATGTATTTATACTTAATTTAATATTATTCCTCTTTTTATTTTTCCTCTTCTTCCTCTTATTATTATTCCTCCTTCTTGATGTTAAAGTTTAAAATCACTGCTAAGGCAGTCTCAAATGGGAGTAGCTTCTAAAAACTTATAAGCCTCATAAAACTAGTCAATTTATTTATATAAACATCAGCTCATCTTTTCCACAAGAGTGTTTATCTCACCCATAATTTCCTGACTATAATGAGTCTTATAATAGGCCATAACCAAATGCTATGCAAAGAATAGCAAAATCTCTTTCTCTCTCTCTCTATATATATCTCTCTATATCTTATATTATTTCTATATATCTATCTCTATCTCTCTATTATCTCTCTCACACACACACACAACACACACAACACACACAACAAAAAGTTCTAAGACAGAAAGTTGTTTCTCTTTACTTCCAAATGATAAGCAAGTTAATCGTATGCATTCAGGCCTGCAATAGAGAACCTGACACACATGAGAGACGTGCAAGTATTTTGCGGAGGATCCGTCAGTCATAAAAACGAAGAGTCTCTTACTTTTGCATGGCTTTCCATCTGTTTATGACTCATGCTCTGTCATCAGGAAACGGTTTTATTTATATGGAGACAGTGAGCAGGAGGATGTGAAACTGTCCACAGATCCACACAGCTTATTACATATGAGGTGAGGCCTCTTGCCCCCAAGCCGGGAATTCTTTATTTAATGTT";

        let aligner = Aligner::builder()
            .map_ont()
            .with_index_threads(4)
            .with_cigar()
            .with_sam_out()
            .with_sam_hit_only()
            .with_seq(target)
            .unwrap();

        // aligner.mapopt.min_cnt = 2;
        // aligner.mapopt.best_n = 500;
        // aligner.mapopt.min_dp_max = 10; // min dp score
        // aligner.mapopt.pri_ratio = 0.5; // min dp score
        println!("{:?}", aligner.mapopt);
        println!("{:?}", aligner.idxopt);

        // tot len
        let seq = b"CAAAGAGGAGAGTTAAACATAAACCAAAATTAAATAAGAATCACGGCTTTGGGGGCCAAGAGTCCTTCCACCTCATATGTAATAAGATGTGTGGATGCTGTGGACATTTCACCACCTCCTGCCCACTGTCACATAATAAAATAAACCGTTTCCTGGATGACAGAGACATAGTATAAACAGATGAAAGCCATGCAAAAGTAAGAGACTCTTCGTTTTTATGACTGACGGATCATCCGCAAAAATACTTGCAACGTCTCTCATTTGTCAGGTCTCTATTGCAGGACCTGAATGCATAAAGATAGACTTGCGTTATCGTTTGGAAGTAAAGAGAAACACTTTTCGGTCTTAGAACTTTTTGTTGTGTGTGTGTGTGTGGTGTGTGTGTGTGTGTGAGAAGAGAGAAGAGAGAGAGAAGAGAGAGAGAGAGAAAGAGAGAGAGACGAGGAGAGGAAGAGAGAGAGGAGAGAGAGAGAGAGAGATTTTGCTATTCCCTTTGCATAGCATTTGGTTATGGCCTATTATAAGAGTATATTAATAGTCAGGGAAATTATGGGGTGAGAGAAACATAATTGTGGAAAGATGAGCTGAGGTTTATAAAGAAATTGACTAGTTTATGGAGGCTATAAAAGTTTTAGAAAATAATACCATTGAGACTGACTTAGCAAGTGATTTTAAGATTTTAAACTTTAACATCAAGAAGGAGAGGAGGGAAGAGGAAGAAGAGGAAGAAAGAAGAAGAGAGAACAAGAGGAAAAGAGGAAGAGGTAGAGAGCGGAAAACGGAGAAAGATGAAAGGAGAAAGACAAGACAAAAAATAATTCCAACTGTTTTCTGTTGTATCCCTTGACTGTATTCCAAGTTTCTGCTTCACTTTTGTTTCCATAGCAACATGCAGAGTAAAACAGGAAGTAAATATTAAGCAAAAGTTGAATGATAAAAAAGAATAAAAAGCCAAGGAAAATGCTTAAAAATCTACTGTGCTCAAGCATCCAACTACCAGTTCCTATGAGATAGCTTACAGAAGATGGTAAATGTATCCATCAGGAAAATGGTTAACAGCAAAAATCCTAGGGGAAAAATGTTTTTAGACAAGAAAATACATTTCTTGACAAAAGTAATTGGCAGTGGCAAGAGATAGAGGAGGATTGGAGAACAAAGGAGACACGAGACATAGATGACATATCAGTTGCAAAGAGACAATTTATCTGGATCCTGGTTTGTTGTTGTTGTTTGTTTGTTTTGGTTTTTTGTATGTTTTGATTATTTTATTTTGATGTGCATTGGGTGGTTTTGCCTGTGGTGTATGTATATGAACCCATGGGATGAATGGTGGTGGCAAAAAGTCAGAAGAGGGTCCTCATCTCCCCCAGAAACTGAGGTAACAGATGGTTGTGGGCCCTAACCTGGGTGCTCTGTCAAGAAGCAGAAAGTGCTCTCAACGCGGAGCCAATCTATTATCCAGCCCCAGACATTCTGTTTTTAAAAAATTATGATATCATTGAGGACATTTAGAAATCTGAACATTGCCAAATTTGGTGGCATGAGTTACAATTGTTTTGCTTGTAATGTTGTAAAAATAAATATTACCATATCTAGACTTTGCATCAAGTCAACGGGAGGAGGACTTGGGCAACATAGGATGGGCTGAGGGTAGGGCTGAGTCTGGAGAGAGTTTCCATGATCACCATAAGGAATCACCCAGTGGGCAAGTTACATGTGTGTATTTGAGCTTCTTACTGAAGTGTGGCAAACAGAAAGCCTGGAAAAGCCCCTGGTGCACAACAGGAAACCATATACAAATAAAAAGAGGGTTCTCCCTCAGTGGGAAAAGGTTTTCGTTTTTTTACATAAGGTATGTGTGATTTGTTTACCTGTGCTCAATTACCAGCAATCAGGAATCTTCAATCTTTTTTCAGAGTAGATTACTCTGAGGAGATGTATTGTTTTTTAACATAGTTGTTGCAGGCAGTCTTTCAATGCAAAGCAAACTGGAAAATGGTATCCGCACAATCGCCAACATAGTAATAGGTAGTGAAGCTTGTTTCGACGTCTACAAAGTTCAAAGTCCTACAAGTTAAGCTGAGTGTATAAGTGCTTACACTGCAGTTAAACATGGGCGTGGGTAACGGAATGTGCACACACTCTGGGAAAAATTAGTTATGTAAAGTTGAACAATGTGATTCTTAGAACGTTGACTTGATCAACCTCACAGTATGGTCCAAGCTTTTAACTGACAACTGTTTTCCCAAAACCCAAAAAGATTATGTGTAAAAGAACAATCCACTTGTATTCTGCTTGCCCTCAGAAAACAAGGTAAGGTCTGTGATTCTCCAGCCAAGCCCTCCGTTGTAGACTCGATGGTTGGCACAAAGCAGGGTTTAAATGACCCAAATGGTGTTTTTATGTTACAGTTTAGTTTGAAGTTTTTGGAAAGCACTCTTACCAGTAACAAGCCCAGATTTGCCAGCTACAGCTAAGTGCTTTTCAAGATTCCTAATGTATGGTCAAGAAGAGTGTGGCTGCTCTCTGCTTCACTATGAAGGAAAGCAAGGCTCAGAGAAATGGGTGGCTCCCAGCCTGTAGGGCAGGCAGAGCAGCAGTCAGGAGTGGGCTGTCTGGATTATTGTTTTCCTAGTGTGCCGTGTGAATACGTGTGTGTGTGTGTGTGGAATACGGTGGTGGTTGTGGTGTGTTCATGGGTGGGATAGCTAAAATGTGTGCATGTGAAGAGCAAACTTCAGATGCCATTAATTTAGTCACATAACCTTTTTCTCCTCTTGAAGTAGGGTCCCTTCTCACTGGACTGGCAGTTAGCCACTACAGCTGACCAACTGATAAGGAACCCACACACCCTCTATACATCATCCCCACACACAACTGGGGGCATTACACCTTGCCATCATGCTTAAATGCTGTTCAACCTAAATATTTAATGAATAAATAATAATTAATTCCGAACAAAAAGTAAATAGATAACTTCCCTCCCTTCTCATACTTCCTTCACATTTTTCGGGTACAGATCAGGCACCACACTTCAAGAAAGCACTTTCTAGGATGAGCCAATATAAAAGCATTATCCTATAGCATTAAATCATTTTCTCATTGCAATTGTAAGAACCATCAGAGAACCTCGGTTCATTAGTCTCCAGGGCTCCTTAACCATAGTAGAGTTATATAAACAATGTCGTATCAACAATTCGTTCAAGGGCTGATGTCCTGCAACTTGTCTAACCTTCTCAAGCATGATGCTTACCAGTGTTTTCCAGGCTCATATTGAAACTGCCAAAGTGTGATTATCCAGATAACTATTTAAGGGAAATGTTTGAGAACACGTCTGCATTCTGAGGCAGAGAGAATATACCACCTTAACAAGGGAAGAGGAATGAAAAAAAAACTGTAGCAAGGTGAGGGGGAACGAGATTGAAATATCGAGGACCAATCGTCAGGTTGATGCAGTTTAAATGATGTAGGCCATAGAGTAAGTGGAAAACCCAACCAAAATGGAGAAATGTATAGAACACTGGAGACCAAGTGGAACACACACACACACAACGGGATGGTTGTAGAAGATAGAATCTCTAGATATAGAAATGGCAGCCCACCTCTCCCTGTTCCACCTGTCCCCAAGAAGTTGGACCTGTGACCTTGAAATTTCCACTGCTTAGAGGAACATTTTCAGGGTGGTCCTTATCCGGATTGAGAGAGAGAGTTCCAGATGAAACAATTGGACAAGTCACGGTGCTAAACCATAGCTGATGTATCCCATGCCCACGAATATATTGTAATCAATTTAATTCTGATACACACAAGCAAAAGGGACGGGAGCAGAGATGTGAGCTTTTGCAGTTACAAGTACTTATTAAATTCTAGGCTCGACGTTTAGCACTATTGCTGGGGGTTTTGCTAGGGAGGATGGCCAGTAGGTGGTTTGGGAGGTGCTACCACACACTTTCCCCCCAGCACTCGGGTGGATAAGGTGGCAATGTGTAGAGTGTTCTGCAGGAGCCAGCTCTAGTCATGTATTAGAGAAAACATTTGGGGATTGCTGCTTTGTTACCATATGTTGATTTCTGTGTCTGGATTCTTGTAGGACTCGCCAGATCCATTTACCCTGGTTCCAGTTTTATGCAAAGCCTGGTTGACTTTTCTGAGAACATCGAGGTTACTGTAAAATCTAGGGCAAGGGCCTAGGCTCCTGGCAGAGTTTCAGTATGAGTATGCTAAGTACAAACTAATTGGTGGGTTATGATATGGTCTCTCACTGAGAAGCAGTGATCTACAAATACTATATCATAACATCAGGAAGATGACCAGGGCAGGGAGTTCCAGGGCAGCCTGGGAGGTGGCTACATAAGCCAACTCTTGTCAAAAAAGATTTGAAAATTATAAATCAATTTGAAAATTATTTTGTTACTAGTAGAAGTGACTAAATTCTGGATCTTTAAGCATGTATAATTTTAGTTTAAAATATATGCTTTATGTGTTCTTACCAAATGCCATTTTTTTCCTGTAGGGTTTCTTAAAAGGTGTTGGCATTTTTGTAAGAAGTGCTGTTAAATAATAAATGGAGTTTGACTTGACTTTTTTAAATAGAATAATTTTATTCCTAAAACAAACCTTGGTCCATCTTTTGTAAAAAGAATAAAGTTAATTTCTGAAGGTTGTCCAATAATGTATCTAATTGTTAATTTTTAAAAAAGGGTGGTCTTGACTCAGCCTGCTATTACAGAATAAATTAGATGTTTTGTATATAATAAACAGTTATATCTCATACTTCTGGAGTCCACTACTAAAATGTTCCCTTTCTCCTTTATCAAGTACAGTATCTATTTTTATAACAGCTGATGCCCCCTGCCCCCAGCTGCTTCCCACAAACTTGGATATTTTGGTGGTGCGAGAAGTCATGAAGCCAACATTTCAAAGGGCTCCCATTTAGCACATGTCTTTACACACGCATACAGGTGCTCTCATATGCATAATATTAGAAAGGCCTGGCAGAATAGCAACGAAGACCAATATTTATAGTATTGTTTAGGAATTAAGTGTTTTAGGTATTTTCTGTGTGTTTAGAACCATTTAAATCCTTTCAAAGCTCAAGTAAAGATGAAACGGCATCCTCAGCACTTCATACAAGGAAGAACCGAGGAACGGAGCCTGAAGGGGTTTGCATGGAGACTGTGTAAGCTAGTGGAAGAAAAATAGCATTCCAATGGTGATCCCAACTTTTGTCCCCTTTGAAGGAAGATTCGCCTTTTAAAGCTCTGGTGCTCAGGTCAAGCACCTGTAATTGGCCGTGAACATTTTGTTTATTTTATTTTATTTCATTTTATTTTATTTATTTATAGTATCAAGAATGATATCCTGCCTCCTAGCCCCTGAGTGTCTGTTCAAGGCTATAGGTATAGAAACACTGATATCCACTTTGAAAGGTTTTTAATGCTCCTTCGGTGCGCTGATGACTGTTGGGAATCGGATGCCTAGCTTTTCTCACAGACCTGTAATGGCTCTGAGCTGATACCACACCAGTAAGGAACAACCAACAGAAGGGGGATCTTGGAAATAAAATGGTAAAAGCTGAATGGACATTGAGCTGTCTACAGGTTTCTTGGGCAAAACTGACTCAGTTGCCATTTCCTCATTCTGCAAAAACAGAGTTTAAAGTAAAAGATATGCTGGCGCTTCTAGCTCTGTGAGTCCATGAACTTATCTTCTATAAAATTCTCTTGTAATTCCGGGAAGTAGATGAAAATTAAGGAGTTGTTTCATGTGGAACCAGTTTCCTTTCCAAGACTTCTTATCAACTGGTGTGTGCCTAAGAAGGCAGATGGTTTTTACCTTTGATTTGGGTGCACCAACCTCAATGAGGGAGGAGGAAACGTGTGCGCGCACACAAATACTCCTGAAAGGTAAGGGGACTGATTTGAGGGATTACAGGTGGGTTCATGATTCCCCACGCCACTGCCCCCTCGCAATGGCTGCGGTAATTGAAAAGCCGCCTGCGTTCTATTGATTGACAGCTGAAGGATCCGTTGCCAGGGAACCTCTGGCCGTGTGCCTAAGGTAACGACCTGTGATTAATTATAACTACTTGTGCTGTCATGCACCCCTCGCTCCCTCCCTCTGAGTGGCGTCCCAGTAATCGGATATCAAATGGAACAGTAAAAATGCCCACACTAAATGAACAATACCTGTGCCTGTGGCAGAACACATGTCAAACAGAAAAAAAAAGACTTTTCGCACCCCCTCCGACAGATCTGCTCATGAATACCTGTTGAAAACAATCATAAATTCAAAGTGGGCTGACCAAACACACACGTTTCCACATAGAAGGGTCCAGAGAAGTAGAGGGTCCAGAAGCCGCAACGTGTTCCACATAGAAGGGTCCCAGCAAAACCTGGCCACACCGTGTTTCCATATAGAAGGTCCAGAAACGGCCACCCGTGTTACACCATAGAAGGTCCGAAACCGGCACACACCCGTGTTTCCACATAGAAGGGTCCAGAAACCTGGCCACACCGTGGTTTCACATAGAAGGTCCAGAAAACCTGGCCAACACCGTGGTTTCCACATAGAAGGGTCCAAAACCTGTTCTTCTCCATCCACTCTGACGCTGCTTGAGATTACCCAAATCTGGAGACTCCACGATAAAAACAACATTCTGCCATGTCCCAAAAGCTTCTGTACTTGGGACTGACTGACTTAAACTAGTAATAGTGCTGGATTTTAGTTTGGAGGGTTACCTAGAATTCTGGACCCCTTCTGAAATTCAGGTACCAAGTGACAAGATTTCTGGACCCTTCCCTATGTGGAAACACGGGTGTGGCCAGGTTCTGACCCTTCTATGGGAAACACGGTGTGGCAGGTTTCTGAAACTTCTATGTGGAAACACGGTGTGGGCAAGGTTTCTGGACATCATATGAACACGGTTGTGCCAGGTTTCTGGACCTCTATGTGGAAACACGGTGTGCCAGAGTTCTGGACTCTTCTATGTGGAACACGGGTGTGGCTGGTCAGACCCCTTGAATTTATGATTGTTTTCAAAGTATGTCACTGATCAGAGCGGTCGAGGGGAGTGCGAAAGCTTTTGACATGTGTTGCAAACAGACACAGGGATTGTTCATTAGTGGGCTTACTGTTCCATTCTGAGATCTGAGGTGAAGCCACTCAGAGGAGGGAGCCGATGGGCAGGAACAGCACAAAGGAGGATAATTAATCAGGTCGTTAGCTTGGCACACGGCCAGAGGTTCCCTGGCACAGATAATTCAGCTGTCAGCAATCAGAACAAGCAGGCAGCTTTCAATTACCGCAGCCATTGGGAGGGGGCAGGACGTGGGGAATCTGACACCAGCTGATAATCCCTCCAAATCAGTCCCTTACCTTCGGATGGATTTGTGATGCGCGCACCAAGTTTCCTCCTCCCTCATTGAGTCCACCCAAAGCAAAGTAAAAACCAGATCTGCCTTCCTTAGGCCACAACCAAGTTGATAAAAGGTTGGCCAAAGGAACTGGTTCCCACAGAACACTCCTGAGGGTCACTACTCACGATTACAGAGAATTTTATAGAAGTAAGTTCATGGACTCACAGAGCTAGAAAGCCCGCAGATCTTTACTTTAACTCTGTTTTTTTCGATGGAGAAATGAATGAGTCAGTTGCCCAGAAAAGATTGGTAAGACAGCGCAATCCAATTCAGCTTTTACCATTTTATTGCCAAGAATACACTTGTGTCTGGTTCCTTTAACTAGATGTGTGTCAGCTCAAGAGCCACTACAGGTCTGTGGGTAGAAAAGCTAGGACAGATTCCCAACAGTCATCAGCGCACCGAAGGAGCATTAAAAACCTTTTCAATGTGGAATAGTGTTCTATACTAAGCATGACAGACACTCAGGGGAGAGGCAGGAAGGATCAGACTTGAGACTATAACATAAAAAAAATAAAATAAATAAAATAAATAAACAAAATGTTCCAGGCCAATTACAGGTGCGACCTGAGCACCAGAGCTTAAGGCGAATTTCATTCAAAGGGGACAAAGTTGGGAATCACCAATTGAGAATGCTCATCTTTTCTTCCACTAGCTACGACAGTCTCAGGCAAACATATCAGGCTCCGTTCCTCGTTCCTCTGTAAGAAGTGCTGAGGATGACCGTTTCATAACTTTACTTGAGCTGTTGAAAGGATTAAATGGGTCTACACACAGAATACCTAAAACAGTTCTAATGCCTTAAACAATACTATAAAAGATGGTCTGTTGCTATTCCTGCCAGGCCTGCTAATATTAGTCATATGAGAGCACGCTGGTATGCAGTGTGTAAAGACATTGTGCTAAATAGGGAGCGCTTCATGAAATGTGGCTTCGACTTCTCGCCCAACAAAATATCCAGGTGTGGAAGCAAGCTGGGGGGCAGGGGGACATCGCTGTTAAAATAGATCTCACTGCAACTTTGATAAAAGAAAGGGAACATTTTAGTAGGTTGGACTCCAGAAGTATGAGAGATAAATTCTTTATATACAAAACCATCTAACAATGGATTCTGTAATAGTCAGGCCTGAGTCAAGACACCTTTTAAAATTACAAGTAGATACATTATTGGCAACCTTCAAAATTAACTTTAATCTTTTTAACAAAAGAAGGGGACCAGTTTTGTTTTAGGAAGACAAGTATGCTATTTAAAAAGTCCAAGTCAACTCCAGTGGATTATTGACAAGCACTTTCTTACAAAAATGCCAACACCTTTTAGAAAACCACAGGAAAAAAATTGTGCATTTGGTAGACACATAAAGCATAGATTTTAAAATTATAATTGATTAAAGATTCCAGATTTAGTCCACTTATCTAGTAACAAACCTAATTTCAACATTGATTTATAATTTTAAAATAATTTACTCTTTTTTGACAAGAGTGGCTATGTAGCCCAGGCTGCCCTGGAACTAAATCCAAATAAATTCGAAACTCACCCCATGCCCTGGTCCATCTTCCTGAATTTTAGGAAATAGGAGTGTATAACACGCTCTCTCAGTGAGAGACTAATCATGAACCCACAATGAGTTTGTACTTTAGCAGTCTCATACTGGAACTTGCCAGAGCCTAGGCCCTTGCCCTTAGAGTTTACAAGTACCTCGATGTTCCAGAAAAGTCAACCAGGCTGGGTGCATAAACTGGAACCAGGATAAATGAGCCTGGCGAGTCCTACAGAATCCAACACAGAGAAGAAAGGATGGTAGACAAAGCAGCAATCCCCAAAGTTTTCTCTAGCATATCCCCGAGTGCTGGGGAAGATGTGTGGTAGCACTCTCCCACCCTACTGGCCATCCTCCCTACAAAACCACAAGAGTGCTAAACCTGAGCCTAGAATTAATAAGTACCTTGTAACTGCAAAGCCACAGTCTCTGCATATCCGGTCCCTTTGCTTGTGTGTACAAGAATTAAATTAGATAACCACATATAACGTTCGTGGCATTGGGATAGCCATCACTGATGGTTTAGCACCGTGCTTGTCCATTGTTTGACTGGAACCTCATCCACCCAATCCGGATACAGACCACACCCAGGAAAAGGTACTCTTAAGGCAGTGGAAAGGCAAACAGGTCCCAACAGGCTGGGACAGAGGGAACAGAGGAGAGGGTGGCTGCCTTTTTATATCTAGAAGTCTATCTTCTACAACCATCCCGTGTGTGTGTGTGTGACACTTGGTTCAGTTTCATAGATTTTCTCCCATTTTGGTTGGGTTTTCCCATTATAACTATAGGGCCATAACAGGGAAAACACTGCATCAACCACTGATTAACAGATGGTCCTTATATTTTCATCTCATTTTCCCTCACTTGTCTTAACAGTGTTTGGTTTCATTACTCTTTTCCCTTGTTAAGTGTAATTCTCTGTCCTCCGATGCAGAGTGTTATCAAAACATTTCCCTTAAATAGTTATCGGATTACCACACTTGGTCAGGTTTCCAAAATTATGAGGCTGAAACCACTGGTAAAGCAGCATGGCTTGAGAAGGTTGAGAAGTTGGCAGGACATCACCTTGAACGATTGTTGAAGCGACAGTGTTTATATAACTCACTAGGTTAAGGAGCCCTGGAGACTAATTGAACGAGTTCTCTGATGGGGCTTTTACAATGTCATGAGAAAATATTAATGCTAATAGGAAAGGCTAGGGAGATGATCCAGCCGAAAGTGATTCCTTGCAAGTGTGGGGACCTGATCTGGTACCCCAGAACAATGGTGAAAGGAAGGAGGAAGGAAGGGAGGAGGGAGGGAGGGAAGGAAAGGAGGGAGGAGGGAGGAAATGAGGGAGGGGGAACAGCATGATGCATGATGGCAAGTGTAATGCCAGTTGTGTGTGGGATGATGTATAGAGGGTTGTGTGGGTTCCATAGTCAGTTGGTCCAGCTGTAGTGGCTAACTGGCCAGTCCAGTGAGAAAACCCTACTCAAGAGAGAGAAAAAGGTGGATAGTGACTAAGGAGATGGCAATCTGAAGTTTTTTGCTCTTCGACATGCACACATTGATCTATCCACCCAGACACACACACACACACACACACATGCACCACAACACACACACGTAATGCAACACGTGCACACTAGGAGAAAACAATAAATTCCAGAACACCCACTACTGACTGCCGCTCTGCCTGCCACTTACAGGCTGGTGAGCCACCCATTTCTCTGAGCCTTGCTTTCCTCATTTAGTGAAGCAGAGAGCAGCACACACTATTCTTGACAGACATTAAAGGAATCTCTGAAAAGACTTAGCTGTTAGATGGCAGCGGCTGGTTACTGGTAAGAGTGTGCTTTCCAAATAACTTACAAACTAAAACTGTAACCTAAAAACACCATTTGGTCATTTACAACCCATGCTTTGTGACCAACCAATTCAGAGTCTACAAACGGAGGCTTGGGCGGAGAACACAGCCTTACGGTTTCTGAGGGCAAGTCAGATACAAGTAAGGATTGTTCTTTTTACACATACCTTTGCGTTTTTGGTGAAAACAGTTGTCAGTTAAAGCTTGACACAGTGTGAGGTTGAGGCACAGGTCAACGTTCTAAGTAAGCAACATTGTGCAACGGACCATAACTAATTTTTCCCGGAGTGCTGTGCACATTCGTTACCACACCTTGTTTAACTGGACAGTGATAAGCACCTTAGGACATCAGCTTACTTGTAGGACTTTGAACTTTGTAGACAGTGAACAAGCTTCACTACCTATTCTAATGTTGGCGATTGTGCCAGGAAGACCATTTTTCCCAGTTGCTTTCATTGAAAGATGCCTGCAACAACTATGTTAAAAACAATACATCCATCCAGAGGTAGGCTACTCTGAAAAAATGTAAGTTCCTGATTGTTTTTTTTGGTATTTGAGACAGGTAAGAAATCACACCATATCCTTATAAACAGAAAACCTTTCCCACTGAGGGGAGAACCCATCTTATTATTCGTATAATGGTTCCTGTTGTGCACCAGGGGCCTTTTTCAGGCTTTCTGTTTGCCACACTTCAAGTAAGAAGCTCAAACACACATTAACTTGGCCCACTGTGGTGATGCCTTATGTGTGAGCATGGAACTCCTCTCCAGACTCAGTCCTACCTCGCCACACTATTGTTGCCCAAGTCCTCCTCGGATGACTTTGAGGCAAAGTCTAGTATGATATCATTTTTATGTGTAACCCACCCCATTCAAGCAAACAATTGATAAATCATGCCACCAAATTGCAAAGGTTCAGATTTCTAAATGTCTCATGATATCATAATTTTTTAAAAAATCAGAATTCGGGCTGGAGAGATAGGATTGGCTCAGCAGTTGAGAGCACTTTCTGCTCTTTGCAGAGCACCAGGTTAGGGTCCCACAACCATCTGTTTACCTCAGTTCTGGGGGAGAGAGACCCTCTTCATGACTTTGCCAGCACAGGCATCCATGGTTCATATACATACACACGGCAAAAACACCAATGCACATAAAAAATAAAATAAATCAAAGACATACAAAAACAAAACAAACAAACAACACCAACAAACCAGGATCACCAGGATAAAGTGTACTCTTTGCAAATGATGTATATGTCTCCTGATGTATACTTTGTTCTCCAATCCCCTTAGCCTTGGTCCACTTTGCAATTACTTTGTTCAAGAATGTAGTTGCTTGTCTAAAAACATTTTCCCTGATTTTGCTGACATATCCTGATGGATACAGTTAACATCTTCTGTTAGATCTTCTAGGAACTGTAGTTGGATGCTTGAGCCACAAGGTAAGATTTTTAAGTCATTTCCTGGCTTTTTATTATAATTTTATTCATTCAACTTTTGCTTAATAGTTTACTACATGTTTTATCTCTGCAGTTTGCTATTGGAAACGAAAAGTGAAAGCAGAAACTTGGAAGACAGTCAAGGGATACAAACCAGAATACGTTGAGTATTTGGTGTAGTAATTTCTCTTTAATTCTTCTTCCTCTTCCTCTTCTTCCCTTTTTAGTATTCTTAATCTTATTATGTCCTTCCCTCCATATCCTTCTTGATGTTAAAGTTTAAATCACTGCTAAGGCAGTCTCATGGGAGTAGCTTCTAAAAACTTATAAGCCTATAAAATAGCAATTTTCTTTGATATAAACCTCAGCTCATCTTTTCCACAGGTGTTTCTCTCACCCATAATTTCCCTGACTATAATGAGTCTTATAATAGGGCCATAACCAAATGATATGCAAAGAATAGCAAAATCTCGCTTCTCTCCTCTCTCTATTATATCTCTCTCTCGCTTTCTATCTCTCTTTCTCTCTCTCTATCGTCTCTCTCTCTCTCCTCTCTCTCACACACACACACACACACACACACACACACACAAACAAAAAGTTCTGAGACAGAAAGTTGTTTATCTTTTACTTTAAAATGATAAGCAAGTTTAACGTATGCATTCAGGCCTGCAATAGGAACATGACACAATGAGAAGACGTGCAAGTATTTTGCGGGGATCCGTCAGTCATAAAAACGAGAGTCTCTTACTTTTGCATGGCTTTCCATCTGTTATGGACTCAGCTCTGCATCAGGAAACGGTTATTGTATATGGAGACAGTAGCAGGAGAATGTGAAACTGTCCACCAGACCACACAGCTTATTACAATTGAGTGAGGACTCTGACCCCAAGCCGGGAATTCTTTATATAATTT";
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
}

/*

sr
aligner.mapopt:mm_mapopt_t { flag: 1077997644, seed: 11, sdust_thres: 0, max_qlen: 0, bw: 100, bw_long: 100, max_gap: 100, max_gap_ref: -1, max_frag_len: 800, max_chain_skip: 25, max_chain_iter: 5000, min_cnt: 2, min_chain_score: 1, chain_gap_scale: 0.8, chain_skip_scale: 0.0, rmq_size_cap: 100000, rmq_inner_dist: 1000, rmq_rescue_size: 1000, rmq_rescue_ratio: 0.1, mask_level: 0.5, mask_len: 2147483647, pri_ratio: 0.5, best_n: 20, alt_drop: 0.15, a: 2, b: 8, q: 12, e: 2, q2: 24, e2: 1, transition: 0, sc_ambi: 1, noncan: 0, junc_bonus: 0, zdrop: 100, zdrop_inv: 100, end_bonus: 10, min_dp_max: 1, min_ksw_len: 0, anchor_ext_len: 20, anchor_ext_shift: 6, max_clip_ratio: 1.0, rank_min_len: 500, rank_frac: 0.9, pe_ori: 1, pe_bonus: 33, mid_occ_frac: 0.0, q_occ_frac: 0.0, min_mid_occ: 0, max_mid_occ: 1000000, mid_occ: 1000, max_occ: 5000, max_max_occ: 4095, occ_dist: 0, mini_batch_size: 50000000, max_sw_mat: 100000000, cap_kalloc: 500000000, split_prefix: 0x0 }

aligner.mapopt:mm_mapopt_t { flag: 1077997644, seed: 11, sdust_thres: 0, max_qlen: 0, bw: 100, bw_long: 100, max_gap: 100, max_gap_ref: -1, max_frag_len: 800, max_chain_skip: 25, max_chain_iter: 5000, min_cnt: 2, min_chain_score: 1, chain_gap_scale: 0.8, chain_skip_scale: 0.0, rmq_size_cap: 100000, rmq_inner_dist: 1000, rmq_rescue_size: 1000, rmq_rescue_ratio: 0.1, mask_level: 0.5, mask_len: 2147483647, pri_ratio: 0.5, best_n: 20, alt_drop: 0.15, a: 2, b: 8, q: 12, e: 2, q2: 24, e2: 1, transition: 0, sc_ambi: 1, noncan: 0, junc_bonus: 0, zdrop: 100, zdrop_inv: 100, end_bonus: 10, min_dp_max: 1, min_ksw_len: 0, anchor_ext_len: 20, anchor_ext_shift: 6, max_clip_ratio: 1.0, rank_min_len: 500, rank_frac: 0.9, pe_ori: 1, pe_bonus: 33, mid_occ_frac: 0.0, q_occ_frac: 0.0, min_mid_occ: 0, max_mid_occ: 1000000, mid_occ: 1000, max_occ: 5000, max_max_occ: 4095, occ_dist: 0, mini_batch_size: 50000000, max_sw_mat: 100000000, cap_kalloc: 500000000, split_prefix: 0x0 }
ont
aligner.mapopt:mm_mapopt_t { flag: 1073741900, seed: 11, sdust_thres: 0, max_qlen: 0, bw: 500, bw_long: 20000, max_gap: 5000, max_gap_ref: -1, max_frag_len: 0, max_chain_skip: 25, max_chain_iter: 5000, min_cnt: 2, min_chain_score: 1, chain_gap_scale: 0.8, chain_skip_scale: 0.0, rmq_size_cap: 100000, rmq_inner_dist: 1000, rmq_rescue_size: 1000, rmq_rescue_ratio: 0.1, mask_level: 0.5, mask_len: 2147483647, pri_ratio: 0.8, best_n: 1, alt_drop: 0.15, a: 2, b: 4, q: 4, e: 2, q2: 24, e2: 1, transition: 0, sc_ambi: 1, noncan: 0, junc_bonus: 0, zdrop: 400, zdrop_inv: 200, end_bonus: -1, min_dp_max: 1, min_ksw_len: 0, anchor_ext_len: 20, anchor_ext_shift: 6, max_clip_ratio: 1.0, rank_min_len: 500, rank_frac: 0.9, pe_ori: 0, pe_bonus: 33, mid_occ_frac: 0.0, q_occ_frac: 0.0, min_mid_occ: 0, max_mid_occ: 1000000, mid_occ: 1000, max_occ: 0, max_max_occ: 4095, occ_dist: 0, mini_batch_size: 500000000, max_sw_mat: 100000000, cap_kalloc: 500000000, split_prefix: 0x0 }


*/
