pub mod cli;
pub mod dna;
pub mod fille_reader;
use std::{collections::HashMap, thread};

use cli::{AlignArgs, IndexArgs, MapArgs, OupArgs, TOverrideAlignerParam};
use crossbeam::channel::{Receiver, Sender};
use fille_reader::{FastaFileReader, QueryRecord};
use minimap2::Aligner;
use rust_htslib::bam::{
    header::HeaderRecord,
    record::{Cigar, CigarString},
    Header, Read,
};

type BamRecord = rust_htslib::bam::record::Record;
type BamWriter = rust_htslib::bam::Writer;
type BamReader = rust_htslib::bam::Reader;

pub struct AlignResult {
    record: BamRecord,
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

pub fn index_ref_file() {}

pub fn query2ref_align(
    query_files: &Vec<&str>,
    ref_file: Option<&str>,
    indexed_ref_file: Option<&str>,
    oup_filename: &str,
    mut aligner: Aligner,
    threads: usize,
) {
    if ref_file.is_none() && indexed_ref_file.is_none() {
        panic!("ref_file and indexed_ref_file cannot all be none");
    }

    if indexed_ref_file.is_some() {
    } else {
    }
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
            let mut bam_h = rust_htslib::bam::Reader::from_path(filename).unwrap();
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
            panic!("invalid file format {}. bam/fa/fasta/fna supported", filename);
        }

        file_idx += 1;
    }
}

pub fn align_worker(
    query_record_recv: Receiver<QueryRecord>,
    align_res_sender: Sender<AlignResult>,
    aligners: &Vec<Aligner>,
) {
    for query_record in query_record_recv {
        align_single_query_to_targets(&query_record, aligners);
    }
}

pub fn align_single_query_to_targets(query_record: &QueryRecord, aligners: &Vec<Aligner>) {
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
            println!("{:?}", hit);
        }
    }
}

/// targetidx: &HashMap<target_name, (idx, target_len)>
pub fn write_bam_worker(
    recv: Receiver<AlignResult>,
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
        writer.write(&align_res.record).unwrap();
    }
}

pub fn build_bam_record_from_mapping(
    hit: &minimap2::Mapping,
    query_record: &QueryRecord,
    target_idx: &HashMap<String, (usize, usize)>,
) -> BamRecord {
    let mut bam_record = BamRecord::new();

    let mut seq = &query_record.sequence;
    let rev_seq = match hit.strand {
        minimap2::Strand::Forward => None,
        minimap2::Strand::Reverse => Some(dna::reverse_complement(seq)),
    };
    if rev_seq.is_some() {
        seq = rev_seq.as_ref().unwrap();
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

    bam_record.set_tid(target_idx.get(hit.target_name.as_ref().unwrap()).unwrap().0 as i32);

    if !hit.is_primary {
        bam_record.set_secondary();
    }
    if hit.is_supplementary {
        bam_record.set_supplementary();
    }

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        println!("hello world");
    }

    #[test]
    fn hello() {}
}
