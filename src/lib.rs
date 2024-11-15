pub mod cli;
pub mod dna;
pub mod fille_reader;
use std::thread;

use cli::{AlignArgs, IndexArgs, MapArgs};
use crossbeam::channel::{Receiver, Sender};
use fille_reader::{FastaFileReader, QueryRecord};
use minimap2::{Aligner, Preset};
use rust_htslib::bam::{
    record::{Cigar, CigarString},
    Read,
};

type BamRecord = rust_htslib::bam::record::Record;

fn build_aligner(
    preset: &str,
    index_args: &IndexArgs,
    map_args: &MapArgs,
    align_args: &AlignArgs,

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
    for filename in filenames {
        if filename.ends_with("fa") || filename.ends_with("fasta") {
            let fasta_reader = FastaFileReader::new(filename.clone());

            for record in fasta_reader {
                sender.send(record).unwrap();
            }
        } else if filename.ends_with("bam") {
            let mut bam_h = rust_htslib::bam::Reader::from_path(filename).unwrap();
            bam_h.set_threads(4).unwrap();
            for record in bam_h.records() {
                let record = record.unwrap();
                sender.send(QueryRecord::from_bam_record(&record)).unwrap();
            }
        }
    }
}

pub fn aligner(query_record_recv: Receiver<QueryRecord>, aligners: &Vec<Aligner>) {
    for query_record in query_record_recv {
        aligner_core(&query_record, aligners);
    }
}

pub fn aligner_core(query_record: &QueryRecord, aligners: &Vec<Aligner>) {
    
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

pub fn writer() {}

pub fn build_bam_record_from_mapping(
    hit: &minimap2::Mapping,
    query_record: &QueryRecord,
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
    let cigar_str = mapping_cigar2_record_cigar(
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

fn mapping_cigar2_record_cigar(
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
}
