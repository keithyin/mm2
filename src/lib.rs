pub mod cli;
pub mod fille_reader;
use std::thread;

use cli::{AlignArgs, IndexArgs, MapArgs};
use crossbeam::channel::{Receiver, Sender};
use fille_reader::{FastaFileReader, QueryRecord};
use minimap2::{Aligner, Preset};
use rust_htslib::bam::Read;

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
            {}
        }
    }
}

pub fn writer() {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        println!("hello world");
    }
}
