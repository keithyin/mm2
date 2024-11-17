pub mod cli;
pub mod dna;
pub mod fille_reader;
pub mod bam_record_ext;
use std::{collections::HashMap, thread};

use cli::{AlignArgs, IndexArgs, MapArgs, OupArgs, TOverrideAlignerParam};
use crossbeam::channel::{Receiver, Sender};
use fille_reader::{FastaFileReader, QueryRecord};
use minimap2::Aligner;
use rust_htslib::bam::{header::HeaderRecord, record::{Cigar, CigarString}, Header, Read};

type BamRecord = rust_htslib::bam::record::Record;
type BamWriter = rust_htslib::bam::Writer;
type BamReader = rust_htslib::bam::Reader;

pub struct SingleQueryAlignResult {
    pub records: Vec<BamRecord>,
}

pub fn build_aligner(
    preset: &str,
    index_args: &IndexArgs,
    map_args: &MapArgs,
    align_args: &AlignArgs,
    oup_args: &OupArgs,
    targets: &Vec<QueryRecord>,
) -> Vec<Aligner> {
    let aligners = thread::scope(|s| {
        let mut handles = vec![];
        for target in targets {
            let hd = s.spawn(|| {
                let mut aligner = Aligner::builder();

                aligner = match preset {
                    "map-ont" => aligner.map_ont(),
                    "map-pb" => aligner.map_pb(),
                    "map-hifi" => aligner.map_hifi(),

                    pre => panic!("not implemented yet {}", pre),
                };

                index_args.modify_aligner(&mut aligner);
                map_args.modify_aligner(&mut aligner);
                align_args.modify_aligner(&mut aligner);
                oup_args.modify_aligner(&mut aligner);

                aligner = aligner
                    .with_index_threads(4)
                    .with_cigar()
                    .with_sam_out()
                    .with_sam_hit_only()
                    .with_seq_and_id(target.sequence.as_bytes(), target.qname.as_bytes())
                    .unwrap();

                aligner
            });
            handles.push(hd);
        }

        handles
            .into_iter()
            .map(|hd| hd.join().unwrap())
            .collect::<Vec<Aligner>>()
    });

    aligners
}

pub fn query_seq_sender(filenames: &Vec<String>, sender: Sender<QueryRecord>) {
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
                    record.qname.push_str(suffix);
                }
                sender.send(record).unwrap();
            }
        } else if filename.ends_with("bam") {
            let mut bam_h = BamReader::from_path(filename).unwrap();
            bam_h.set_threads(4).unwrap();
            for record in bam_h.records() {
                let record = record.unwrap();
                sender
                    .send(QueryRecord::from_bam_record(
                        &record,
                        qname_suffix.as_ref().map(|v| v.as_str()),
                    ))
                    .unwrap();
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
    query_record_recv: Receiver<QueryRecord>,
    align_res_sender: Sender<SingleQueryAlignResult>,
    aligners: &Vec<Aligner>,
    target_idx: &HashMap<String, (usize, usize)>,
) {
    for query_record in query_record_recv {
        let records = align_single_query_to_targets(&query_record, aligners, target_idx);
        align_res_sender
            .send(SingleQueryAlignResult { records })
            .unwrap();
    }
}

pub fn align_single_query_to_targets(
    query_record: &QueryRecord,
    aligners: &Vec<Aligner>,
    target_idx: &HashMap<String, (usize, usize)>,
) -> Vec<BamRecord> {
    let mut align_records = vec![];
    for aligner in aligners {
        for hit in aligner
            .map(
                query_record.sequence.as_bytes(),
                true,
                true,
                None,
                Some(&[67108864]), // 67108864 eqx
            )
            .unwrap()
        {
            if hit.alignment.is_none() {
                continue;
            }

            let record = build_bam_record_from_mapping(&hit, query_record, target_idx);
            align_records.push(record);
        }
    }

    align_records
}

/// targetidx: &HashMap<target_name, (idx, target_len)>
pub fn write_bam_worker(
    recv: Receiver<SingleQueryAlignResult>,
    target_idx: &HashMap<String, (usize, usize)>,
    o_path: &str,
) {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.5");
    hd.push_tag(b"SO", "unknown");
    header.push_record(&hd);

    let mut targets = target_idx.iter().collect::<Vec<_>>();
    targets.sort_by_key(|tup| tup.1 .0);

    for target in targets {
        let mut hd = HeaderRecord::new(b"SQ");
        hd.push_tag(b"SN", &target.0);
        hd.push_tag(b"LN", target.1 .1);
        header.push_record(&hd);
    }

    let mut writer = BamWriter::from_path(o_path, &header, rust_htslib::bam::Format::Bam).unwrap();
    writer.set_threads(4).unwrap();
    for align_res in recv {
        for record in align_res.records {
            writer.write(&record).unwrap();
        }
    }
}

pub fn build_bam_record_from_mapping(
    hit: &minimap2::Mapping,
    query_record: &QueryRecord,
    target_idx: &HashMap<String, (usize, usize)>,
) -> BamRecord {

    // println!("{:?}", hit);

    let mut bam_record = BamRecord::new();

    let mut seq = &query_record.sequence;
    let rev_seq = match hit.strand {
        minimap2::Strand::Forward => None,
        minimap2::Strand::Reverse => Some(dna::reverse_complement(seq)),
    };
    if rev_seq.is_some() {
        seq = rev_seq.as_ref().unwrap();
        bam_record.set_reverse();
    }

    let aln_info = hit.alignment.as_ref().unwrap();
    let cigar_str = convert_mapping_cigar_to_record_cigar(
        aln_info.cigar.as_ref().unwrap(),
        hit.query_start as usize,
        hit.query_end as usize,
        seq.len(),
    );

    bam_record.set(
        query_record.qname.as_bytes(),
        Some(&cigar_str),
        seq.as_bytes(),
        &vec![255; seq.len()],
    );

    // reference start
    bam_record.set_pos(hit.target_start as i64);
    bam_record.set_mpos(-1);
    // bam_record.set_mpos(mpos);
    // bam_record.reference_end()
    bam_record.set_mapq(hit.mapq as u8);

    bam_record.set_tid(target_idx.get(hit.target_name.as_ref().unwrap()).unwrap().0 as i32);
    bam_record.set_mtid(-1);

    if !hit.is_primary {
        bam_record.set_secondary();
    } else {
        bam_record.unset_secondary();
    }

    if hit.is_supplementary {
        bam_record.set_supplementary();
    } else {
        bam_record.unset_supplementary();
    }
    bam_record.unset_unmapped();

    bam_record
        .push_aux(
            b"cs",
            rust_htslib::bam::record::Aux::String(aln_info.cs.as_ref().unwrap()),
        )
        .unwrap();

    bam_record
        .push_aux(
            b"md",
            rust_htslib::bam::record::Aux::String(aln_info.md.as_ref().unwrap()),
        )
        .unwrap();

    bam_record
}

fn convert_mapping_cigar_to_record_cigar(
    mapping_cigar: &[(u32, u8)],
    query_start: usize,
    query_end: usize,
    query_len: usize,
) -> CigarString {
    let mut cigar_str = CigarString(vec![]);

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
pub fn targets_to_targetsidx(targets: &Vec<QueryRecord>) -> HashMap<String, (usize, usize)> {
    let mut target2idx = HashMap::new();
    targets.iter().enumerate().for_each(|(idx, target)| {
        target2idx.insert(target.qname.clone(), (idx, target.sequence.len()));
    });
    target2idx
}

#[cfg(test)]
mod tests {
    use bam_record_ext::record2str;
    use fille_reader::read_fasta;

    use super::*;

    #[test]
    fn test_align_single_query_to_target() {
        let ref_file = "test_data/GCF_000005845.2_ASM584v2_genomic.fna";
        let targets = read_fasta(ref_file).unwrap();
        let aligners = build_aligner(
            "map-ont",
            &IndexArgs::default(),
            &MapArgs::default(),
            &AlignArgs::default(),
            &OupArgs::default(),
            &targets,
        );
        let target2idx = targets_to_targetsidx(&targets);

        let query_file = "test_data/ecoli_query.fa";
        let query_filter_iter = FastaFileReader::new(query_file.to_string());
        for qr in query_filter_iter {
            let records = align_single_query_to_targets(&qr, &aligners, &target2idx);
            for record in records {
                println!("{:?}", record2str(&record));
            }
        }
    }
}
