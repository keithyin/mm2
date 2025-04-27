use std::{
    fmt::Display,
    ops::{Deref, DerefMut},
};

use bio::bio_types::strand::ReqStrand;
use gskits::{
    dna::reverse_complement,
    gsbam::bam_record_ext::{BamRecord, BamRecordExt},
};
use minimap2::{Mapping, Strand};
use rust_htslib::bam::ext::BamRecordExtensions;

use crate::convert_mapping_cigar_to_record_cigar;

pub mod align_metric;

pub struct TseqAndRecord {
    pub ori_rstart: usize, // 用来存储 当前比对片段 对应 原始比对区域的 rstart
    pub ori_rend: usize,   // 用来存储 当前比对片段 对应 原始比对区域的 rend
    pub ori_qstart: usize, // 用来存储 当前比对片段 对应 原始比对区域的 qstart
    pub ori_qend: usize,   // 用来存储 当前比对片段 对应 原始比对区域的 qend
    pub tseq: String,      // 仅保留了比对区域的 reference sequence
    pub record: BamRecord, // 仅保留了 比对区域的 比对信息。等价于将 比对区域的 ref 和 query 都截取下来，然后重新比对之后的结果
}
impl TseqAndRecord {
    pub fn new(
        ori_rstart: usize,
        ori_rend: usize,
        ori_qstart: usize,
        ori_qend: usize,
        tseq: String,
        record: BamRecord,
    ) -> Self {
        match record.strand() {
            ReqStrand::Reverse => panic!("only forward strand is supported"),
            ReqStrand::Forward => (),
        }

        Self {
            ori_rstart,
            ori_rend,
            ori_qstart,
            ori_qend,
            tseq,
            record,
        }
    }

    pub fn alignment_pair(&self, qstart: Option<usize>, qend: Option<usize>) -> (String, String) {
        let record_ext = BamRecordExt::new(&self.record);
        let rstart = record_ext.reference_start() as i64;
        let rend = record_ext.reference_end() as i64;
        let qstart = qstart
            .map(|v| v - self.ori_qstart)
            .unwrap_or(record_ext.query_alignment_start()) as i64;
        let qend = qend
            .map(|v| v - self.ori_qstart)
            .unwrap_or(record_ext.query_alignment_end()) as i64;

        let mut r_cursor = None;
        let mut q_cursor = None;
        let qseq = record_ext.get_seq();
        let qseq = qseq.as_bytes();
        let rseq = self.tseq.as_bytes();

        let mut aligned_ref = String::new();
        let mut aligned_qry = String::new();
        for [qpos, rpos] in self.record.aligned_pairs_full() {
            if rpos.is_some() {
                r_cursor = rpos;
            }
            if qpos.is_some() {
                q_cursor = qpos;
            }

            if r_cursor.unwrap_or(0) >= rend || q_cursor.unwrap_or(0) >= qend {
                break;
            }

            if r_cursor.unwrap_or(0) >= rstart && q_cursor.unwrap_or(0) >= qstart {
                if let Some(qpos_) = qpos.map(|v| v as usize) {
                    aligned_qry.push(qseq[qpos_] as char);
                } else {
                    aligned_qry.push('-');
                }

                if let Some(rpos_) = rpos.map(|v| v as usize) {
                    aligned_ref.push(rseq[rpos_] as char);
                } else {
                    aligned_ref.push('-');
                }
            }
        }

        (aligned_ref, aligned_qry)
    }

    pub fn alignment_str(&self, qstart: Option<usize>, qend: Option<usize>) -> String {
        let (aligned_ref, aligned_qry) = self.alignment_pair(qstart, qend);

        let mut oup = String::new();
        oup.push_str(&format!(
            "ref:{}-{}\nqry:{}-{}\n",
            self.ori_rstart, self.ori_rend, self.ori_qstart, self.ori_qend
        ));
        let stepsize = 50;
        let numstep = (aligned_ref.len() - 1 + stepsize) / stepsize;
        (0..numstep).into_iter().for_each(|idx| {
            let start_pos = idx * stepsize;
            let end_pos = (start_pos + stepsize).min(aligned_ref.len());
            oup.push_str(&format!("ref:{}\n", &aligned_ref[start_pos..end_pos]));
            oup.push_str(&format!("qry:{}\n", &aligned_qry[start_pos..end_pos]));
            oup.push_str("\n");
        });
        oup
    }
}

/// 从 Mapping 中直接抽取出来的信息，使用该类作为一个中间存储
#[derive(Debug)]
pub struct AlignInfo {
    pub rstart: usize,
    pub rend: usize,
    pub qstart: usize,
    pub qend: usize,
    pub fwd: bool,
    pub primary: bool,
    pub cigar: Vec<(u32, u8)>,
}

impl Display for AlignInfo {
    /// 1:1000->1000:2000+ means query span 1~1000 forward aligned to ref span 1000:2000
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}->{}:{}{}",
            self.qstart,
            self.qend,
            self.rstart,
            self.rend,
            if self.fwd { "+" } else { "-" }
        )
    }
}

impl AlignInfo {
    pub fn fwd_record(&self, target_seq: &str, query_seq: &str) -> TseqAndRecord {
        let aligned_target_seq = if self.fwd {
            target_seq[self.rstart..self.rend].to_string()
        } else {
            String::from_utf8(reverse_complement(target_seq[self.rstart..self.rend].as_bytes())).unwrap()
        };

        let query_seq = query_seq[self.qstart..self.qend].to_string();
        let mut record = BamRecord::new();

        let mapping_cigar = if self.fwd {
            self.cigar.clone()
        } else {
            self.cigar.iter().rev().copied().collect::<Vec<_>>()
        };

        let cigar_str = convert_mapping_cigar_to_record_cigar(
            &mapping_cigar,
            0,
            query_seq.len(),
            query_seq.len(),
            false,
        );
        record.set_pos(0);

        record.set(
            b"-",
            Some(&cigar_str),
            query_seq.as_bytes(),
            &vec![255; query_seq.len()],
        );
        let (ori_rstart, ori_rend) = if self.fwd {
            (self.rstart, self.rend)
        } else {
            (self.rend, self.rstart)
        };
        TseqAndRecord::new(
            ori_rstart,
            ori_rend,
            self.qstart,
            self.qend,
            aligned_target_seq,
            record,
        )
    }
}

impl From<Mapping> for AlignInfo {
    fn from(value: Mapping) -> Self {
        let aln = value.alignment.unwrap();
        let fwd = match value.strand {
            Strand::Forward => true,
            Strand::Reverse => false,
        };
        Self {
            rstart: value.target_start as usize,
            rend: value.target_end as usize,
            qstart: value.query_start as usize,
            qend: value.query_end as usize,
            fwd: fwd,
            primary: value.is_primary,
            cigar: aln.cigar.unwrap(),
        }
    }
}

/// single query 的比对可能存在多个比对结果。1. 不同片段比对到ref的不同位置。2. 不同片段比对到了ref的相同位置。etc...
#[derive(Debug)]
pub struct SingleQueryAlignInfo(pub Vec<AlignInfo>);
impl Deref for SingleQueryAlignInfo {
    type Target = Vec<AlignInfo>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for SingleQueryAlignInfo {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
