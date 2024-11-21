pub mod cli;
pub mod dna;
pub mod fille_reader;

use std::{collections::HashMap, fs, io::BufReader, thread};

use gskits::{gsbam::bam_record_ext::BamRecordExt, pbar};
use cli::{AlignArgs, IndexArgs, MapArgs, OupArgs, TOverrideAlignerParam};
use crossbeam::channel::{Receiver, Sender};
use fille_reader::{FastaFileReader, FastqReaderIter, FastqRecord, QueryRecord};
use minimap2::Aligner;
use gskits::pbar::DEFAULT_INTERVAL;
use rust_htslib::bam::{
    header::HeaderRecord,
    record::{Aux, Cigar, CigarString},
    Header, Read,
};
use gskits::utils::command_line_str;

pub type BamRecord = rust_htslib::bam::record::Record;
pub type BamWriter = rust_htslib::bam::Writer;
pub type BamReader = rust_htslib::bam::Reader;

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
        } else if filename.ends_with("fq") || filename.ends_with("fastq") {
            let fastq_file = fs::File::open(&filename).unwrap();
            let mut buf_reader = BufReader::new(fastq_file);
            let fastq_iter = FastqReaderIter::new(&mut buf_reader);

            for FastqRecord(qname, seq, _) in fastq_iter {
                let mut record = QueryRecord {
                    qname: qname,
                    sequence: seq,
                    ch: None,
                    np: None,
                };

                if let Some(suffix) = &qname_suffix {
                    record.qname.push_str(suffix);
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

    set_primary_alignment(&mut align_records);

    align_records
}

/// targetidx: &HashMap<target_name, (idx, target_len)>
pub fn write_bam_worker(
    recv: Receiver<SingleQueryAlignResult>,
    target_idx: &HashMap<String, (usize, usize)>,
    o_path: &str,
    oup_args: &OupArgs,
    enable_pb: bool,
) {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.5");
    hd.push_tag(b"SO", "unknown");
    header.push_record(&hd);

    let mut hd = HeaderRecord::new(b"PG");
    hd.push_tag(b"ID", "gsmm2")
        .push_tag(b"PN", "gsmm2")
        .push_tag(b"CL", &command_line_str())
        .push_tag(b"VN", &env!("CARGO_PKG_VERSION"))
        ;
    header.push_record(&hd);


    let mut targets = target_idx.iter().collect::<Vec<_>>();
    targets.sort_by_key(|tup| tup.1 .0);

    for target in targets {
        let mut hd = HeaderRecord::new(b"SQ");
        hd.push_tag(b"SN", &target.0);
        hd.push_tag(b"LN", target.1 .1);
        header.push_record(&hd);
    }

    let pb = if enable_pb {
        Some(pbar::get_spin_pb(
            format!("writing alignment result"),
            DEFAULT_INTERVAL,
        ))
    } else {
        None
    };

    let mut writer = BamWriter::from_path(o_path, &header, rust_htslib::bam::Format::Bam).unwrap();
    writer.set_threads(4).unwrap();

    for align_res in recv {
        pb.as_ref().unwrap().inc(1);
        for mut record in align_res.records {
            if oup_args.valid(&record) {
                let record_ext = BamRecordExt::new(&record);
                let iy = record_ext.compute_identity();
                let ec = record_ext.compute_query_coverage();

                if iy < oup_args.oup_identity_threshold {
                    continue;
                }

                if ec < oup_args.oup_coverage_threshold {
                    continue;
                }

                record.push_aux(b"iy", Aux::Float(iy)).unwrap();
                record.push_aux(b"ec", Aux::Float(ec)).unwrap();

                writer.write(&record).unwrap();
            }
        }
    }
    pb.as_ref().unwrap().finish();
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
        rev_seq.is_some()
    );

    bam_record.set(
        query_record.qname.as_bytes(),
        Some(&cigar_str),
        seq.as_bytes(),
        &vec![255; seq.len()],
    );

    match hit.strand {
        minimap2::Strand::Reverse => bam_record.set_reverse(),
        _ => {}
    }

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
        bam_record.set_supplementary();
    } else {
        bam_record.unset_secondary();
        bam_record.unset_supplementary();
    }

    if hit.is_supplementary {
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
        bam_record.push_aux(b"ch", Aux::U16(ch_ as u16)).unwrap();

    }

    bam_record
}

pub fn convert_mapping_cigar_to_record_cigar(
    mapping_cigar: &[(u32, u8)],
    mut query_start: usize,
    mut query_end: usize,
    query_len: usize,
    is_rev: bool
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
pub fn targets_to_targetsidx(targets: &Vec<QueryRecord>) -> HashMap<String, (usize, usize)> {
    let mut target2idx = HashMap::new();
    targets.iter().enumerate().for_each(|(idx, target)| {
        target2idx.insert(target.qname.clone(), (idx, target.sequence.len()));
    });
    target2idx
}

pub fn set_primary_alignment(records: &mut Vec<BamRecord>) {
    let mut primary_reocrds = records
        .iter_mut()
        .filter(|record| !record.is_secondary())
        .collect::<Vec<_>>();
    if primary_reocrds.len() == 1 {
        return;
    }

    primary_reocrds.sort_by_key(|record| {
        let matched = record
            .cigar()
            .iter()
            .map(|cigar| match *cigar {
                Cigar::Equal(n) | Cigar::Match(n) => n as i64,
                _ => 0,
            })
            .reduce(|a, b| (a + b))
            .unwrap();
        -matched
    });

    primary_reocrds.iter_mut().skip(1).for_each(|record| {
        record.set_secondary();
        // record.set_supplementary();
    });
}

#[cfg(test)]
mod tests {
    use gskits::gsbam::bam_record_ext::record2str;
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

        set_primary_alignment(&mut records);

        for record in &records {
            println!("primary: {}", !record.is_secondary());
        }
    }
}
