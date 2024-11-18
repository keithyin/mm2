use std::fmt::Display;

use rust_htslib::bam::{ext::BamRecordExtensions, record::Cigar, Record};

use crate::BamRecord;

pub fn record2str(record: &Record) -> String {
    format!(
        "refstart:{}, refend:{}, seq:{}",
        record.reference_start(),
        record.reference_end(),
        String::from_utf8(record.seq().as_bytes()).unwrap()
    )
}

pub struct BamRecordExt<'a> {
    bam_record: &'a BamRecord,
    qname: Option<String>,
    seq: Option<String>,
}

impl<'a> BamRecordExt<'a> {
    pub fn new(bam_record: &'a BamRecord) -> Self {
        BamRecordExt {
            bam_record,
            qname: None,
            seq: None,
        }
    }

    /// eq / (eq + diff + ins + del)
    pub fn compute_identity(&self) -> f32 {
        let mut aligned_span = 0;
        let mut matched = 0;
        self.bam_record
            .cigar()
            .iter()
            .for_each(|cigar| match *cigar {
                Cigar::Equal(n) => {
                    matched += n;
                    aligned_span += n;
                }
                Cigar::Diff(n) => aligned_span += n,
                Cigar::Del(n) | Cigar::Ins(n) => aligned_span += n,
                Cigar::Match(_) => panic!("must use eqx cigar"),
                _ => {}
            });

        aligned_span = if aligned_span > 0 { aligned_span } else { 1 };

        matched as f32 / aligned_span as f32
    }

    /// (eq + diff + ins) / query_len
    pub fn compute_query_coverage(&self) -> f32 {
        let mut seq_len = self.bam_record.seq_len_from_cigar(true);
        let aligned_span = self
            .bam_record
            .cigar()
            .iter()
            .map(|cigar| match *cigar {
                Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Match(n) | Cigar::Ins(n) => n,
                _ => 0,
            })
            .reduce(|acc, b| acc + b)
            .unwrap_or(0);

        seq_len = if seq_len > 0 { seq_len } else { 1 };

        aligned_span as f32 / seq_len as f32
    }

    pub fn get_seq_cached(&mut self) -> &str {
        if self.seq.is_none() {
            self.seq = Some(self.get_seq());
        }

        self.seq.as_ref().unwrap()
    }

    pub fn get_qname_cached(&mut self) -> &str {
        if self.qname.is_none() {
            self.qname = Some(self.get_qname());
        }
        self.qname.as_ref().unwrap()
    }

    pub fn get_seq(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.bam_record.seq().as_bytes()) }
    }

    pub fn get_qname(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.bam_record.qname().to_vec()) }
    }


    pub fn get_ch(&self) -> Option<usize> {
        None
    }

    pub fn get_np(&self) -> Option<usize> {
        None
    }
    
}

impl<'a> Display for BamRecordExt<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "refstart:{}, refend:{}, seq:{}",
            self.bam_record.reference_start(),
            self.bam_record.reference_end(),
            self.get_seq()
        )
    }
}
