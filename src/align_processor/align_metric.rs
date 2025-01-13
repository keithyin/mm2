use gskits::{dna::SEQ_NT4_TABLE, ds::region::Regions};
use lazy_static::lazy_static;
use std::fmt::Display;

use super::{AlignInfo, SingleQueryAlignInfo};

#[derive(Debug)]
pub struct Metric {
    qname: String,
    rname: String,
    qlen: usize,
    align_infos: SingleQueryAlignInfo,

    matched_cnt: [usize; 4],
    mismatched_cnt: [usize; 4],
    homodel_cnt: [usize; 4],
    non_homodel_cnt: [usize; 4],
    homoins_cnt: [usize; 4],
    non_homoins_cnt: [usize; 4],

    num_segs: usize,
    q_ovlp: usize,
    q_ovlp_ratio: f32,
    r_ovlp_ratio: f32,
    covlen: usize,
    primary_covlen: usize,
    ori_align_info: String,
    ori_q_gaps: String,
    merged_qry_span: String,
}
impl Metric {
    pub fn new(qname: String, qlen: usize, rname: String) -> Self {
        Self {
            qname,
            rname,
            qlen,
            align_infos: SingleQueryAlignInfo(vec![]),
            matched_cnt: [0; 4],
            mismatched_cnt: [0; 4],
            homodel_cnt: [0; 4],
            non_homodel_cnt: [0; 4],
            homoins_cnt: [0; 4],
            non_homoins_cnt: [0; 4],
            q_ovlp: 0,
            q_ovlp_ratio: 0.0,
            r_ovlp_ratio: 0.0,
            num_segs: 0,
            covlen: 0,
            primary_covlen: 0,
            ori_align_info: "".to_string(),
            ori_q_gaps: "".to_string(),
            merged_qry_span: "".to_string(),
        }
    }

    pub fn add_align_info(&mut self, align_info: AlignInfo) {
        self.align_infos.push(align_info);
    }

    pub fn compute_metric(&mut self, target_seq: &str, query_seq: &str) {
        if self.align_infos.is_empty() {
            return;
        }

        self.align_infos.iter().for_each(|v| {
            if v.primary {
                self.primary_covlen = v.qend - v.qstart;
            }
        });

        let mut tseq_and_records = self
            .align_infos
            .iter()
            .map(|v| v.fwd_record(target_seq, query_seq))
            .collect::<Vec<_>>();

        tseq_and_records.sort_by_key(|v| v.ori_qstart);

        let qstart_ends = tseq_and_records
            .iter()
            .map(|v| (v.ori_qstart, v.ori_qend))
            .collect::<Vec<(usize, usize)>>();

        let rstart_ends = tseq_and_records
            .iter()
            .map(|v| (v.ori_rstart.min(v.ori_rend), v.ori_rend.max(v.ori_rstart)))
            .collect::<Vec<(usize, usize)>>();

        let qstart_ends_region: Regions = (&qstart_ends).into();
        self.q_ovlp_ratio = qstart_ends_region.ovlp_ratio();

        self.ori_q_gaps = qstart_ends_region
            .gaps(Some(0), Some(self.qlen))
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");

        let rstart_ends_region: Regions = (&rstart_ends).into();
        self.r_ovlp_ratio = rstart_ends_region.ovlp_ratio();

        let mut truncated_qstarts_ends = vec![qstart_ends[0].clone()];
        let mut ovlp_len = 0;
        qstart_ends
            .iter()
            .skip(1)
            .for_each(|&(mut cur_start, cur_end)| {
                let (pre_start, pre_end) = truncated_qstarts_ends.last_mut().unwrap();
                if cur_start >= *pre_end {
                    truncated_qstarts_ends.push((cur_start, cur_end));
                } else {
                    ovlp_len += *pre_end - cur_start;

                    if (cur_end - cur_start) > (*pre_end - *pre_start) {
                        *pre_end = cur_start;
                    } else {
                        cur_start = *pre_end;
                    }

                    truncated_qstarts_ends.push((cur_start, cur_end));
                }
            });

        self.num_segs = truncated_qstarts_ends.len();
        let coverlen = truncated_qstarts_ends
            .iter()
            .map(|&(s, e)| e - s)
            .sum::<usize>();
        self.covlen = coverlen;
        tseq_and_records
            .iter()
            .zip(truncated_qstarts_ends.iter())
            .for_each(|(info, se)| {
                let align_pair = info.alignment_pair(Some(se.0), Some(se.1));
                fill_metric_from_align_pair(align_pair.0.as_bytes(), align_pair.1.as_bytes(), self);
            });

        let cov2 = self.matched_cnt.iter().copied().sum::<usize>()
            + self.mismatched_cnt.iter().copied().sum::<usize>()
            + self.homoins_cnt.iter().copied().sum::<usize>()
            + self.non_homoins_cnt.iter().copied().sum::<usize>();

        assert_eq!(cov2, coverlen);

        self.q_ovlp = ovlp_len;

        let mut ori_align_info = self
            .align_infos
            .iter()
            .map(|v| (v.qstart, format!("{}", v)))
            .collect::<Vec<_>>();
        ori_align_info.sort_by_key(|v| v.0);
        self.ori_align_info = ori_align_info
            .iter()
            .map(|v| v.1.as_str())
            .collect::<Vec<_>>()
            .join(";");
        // println!("{}", self.ori_align_info);
        self.merged_qry_span = truncated_qstarts_ends
            .iter()
            .map(|v| format!("{}:{}", v.0, v.1))
            .collect::<Vec<_>>()
            .join(";");
    }

    pub fn matched(&self) -> usize {
        self.matched_cnt.iter().copied().sum::<usize>()
    }

    pub fn mismatched(&self) -> usize {
        self.mismatched_cnt.iter().copied().sum::<usize>()
    }

    pub fn homoins(&self) -> usize {
        self.homoins_cnt.iter().copied().sum::<usize>()
    }

    pub fn non_homoins(&self) -> usize {
        self.non_homoins_cnt.iter().copied().sum::<usize>()
    }

    pub fn homodel(&self) -> usize {
        self.homodel_cnt.iter().copied().sum::<usize>()
    }

    pub fn non_homodel(&self) -> usize {
        self.non_homodel_cnt.iter().copied().sum::<usize>()
    }

    pub fn aligned_span(&self) -> usize {
        self.matched()
            + self.mismatched()
            + self.homoins()
            + self.non_homoins()
            + self.homodel()
            + self.non_homodel()
    }
}

lazy_static! {
    pub static ref METRIC_CSV_HEADER: Vec<String> = {
        let mut csv_header = vec![
            "qname".to_string(),
            "rname".to_string(),
            "qlen".to_string(),
            "segs".to_string(),
            "covlen".to_string(),
            "primaryCovlen".to_string(),
            "queryCoverage".to_string(),
            "identity".to_string(),
            "oriAlignInfo".to_string(),
            "oriQGaps".to_string(),
            "qOvlp".to_string(),
            "qOvlpRatio".to_string(),
            "rOvlpRatio".to_string(),
            "mergedQrySpan".to_string(),
            "match".to_string(),
            "misMatch".to_string(),
            "ins".to_string(),
            "homoIns".to_string(),
            "del".to_string(),
            "homoDel".to_string(),
        ];
        for base in ["A", "C", "G", "T"] {
            csv_header.push(format!("match-{}", base));
            csv_header.push(format!("misMatch-{}", base));
            csv_header.push(format!("ins-{}", base));
            csv_header.push(format!("homoIns-{}", base));
            csv_header.push(format!("del-{}", base));
            csv_header.push(format!("homoDel-{}", base));
        }
        csv_header
    };
}
impl Display for Metric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut res_str = String::new();
        res_str.push_str(&self.qname);
        res_str.push_str("\t");

        res_str.push_str(&self.rname);
        res_str.push_str("\t");

        res_str.push_str(&format!("{}\t", self.qlen));
        res_str.push_str(&format!("{}\t", self.num_segs));
        res_str.push_str(&format!("{}\t", self.covlen));
        res_str.push_str(&format!("{}\t", self.primary_covlen));
        res_str.push_str(&format!("{:.6}\t", self.covlen as f64 / self.qlen as f64));
        res_str.push_str(&format!(
            "{:.6}\t",
            self.matched() as f64 / self.aligned_span() as f64
        ));

        res_str.push_str(&self.ori_align_info);
        res_str.push_str("\t");
        res_str.push_str(&self.ori_q_gaps);
        res_str.push_str("\t");

        res_str.push_str(&format!("{}\t", self.q_ovlp));
        res_str.push_str(&format!("{:.6}\t", self.q_ovlp_ratio));
        res_str.push_str(&format!("{:.6}\t", self.r_ovlp_ratio));

        res_str.push_str(self.merged_qry_span.as_str());
        res_str.push_str("\t");

        res_str.push_str(&format!("{}\t", self.matched()));
        res_str.push_str(&format!("{}\t", self.mismatched()));
        res_str.push_str(&format!("{}\t", self.non_homoins()));
        res_str.push_str(&format!("{}\t", self.homoins()));
        res_str.push_str(&format!("{}\t", self.non_homodel()));
        res_str.push_str(&format!("{}\t", self.homodel()));

        for idx in 0..4 {
            res_str.push_str(&format!("{}\t", self.matched_cnt[idx]));
            res_str.push_str(&format!("{}\t", self.mismatched_cnt[idx]));
            res_str.push_str(&format!("{}\t", self.non_homoins_cnt[idx]));
            res_str.push_str(&format!("{}\t", self.homoins_cnt[idx]));
            res_str.push_str(&format!("{}\t", self.non_homodel_cnt[idx]));
            res_str.push_str(&format!("{}\t", self.homodel_cnt[idx]));
        }

        res_str.pop();
        write!(f, "{}", res_str)
    }
}

pub fn fill_metric_from_align_pair(aligned_ref: &[u8], aligned_qry: &[u8], metric: &mut Metric) {
    assert_eq!(aligned_ref.len(), aligned_qry.len());

    let mut ins_bases_idx = vec![];
    let mut del_bases_idx = vec![];

    for idx in 0..aligned_qry.len() {
        let qbase = aligned_qry[idx];
        let rbase = aligned_ref[idx];

        match (qbase as char, rbase as char) {
            ('-', '-') => panic!(""),

            ('-', _rbase) => {
                // deletion
                del_bases_idx.push(idx);
                fill_insertion_info(&ins_bases_idx, aligned_ref, aligned_qry, metric);
                ins_bases_idx.clear();
            }

            (_qbase, '-') => {
                // insertion
                ins_bases_idx.push(idx);
                fill_deletion_info_v2(&del_bases_idx, aligned_ref, metric);
                del_bases_idx.clear();
            }

            (qbase, rbase) => {
                if rbase == qbase {
                    metric.matched_cnt[SEQ_NT4_TABLE[qbase as usize] as usize] += 1;
                } else {
                    metric.mismatched_cnt[SEQ_NT4_TABLE[qbase as usize] as usize] += 1;
                }

                fill_deletion_info_v2(&del_bases_idx, aligned_ref, metric);
                fill_insertion_info(&ins_bases_idx, aligned_ref, aligned_qry, metric);

                del_bases_idx.clear();
                ins_bases_idx.clear();
            }
        }
    }

    fill_deletion_info_v2(&del_bases_idx, aligned_ref, metric);
    fill_insertion_info(&ins_bases_idx, aligned_ref, aligned_qry, metric);
}



/// homo 判定逻辑：1. del区域base需一致，2. del区域的 base 和其 左/右 临近的非 del 区域的base一致
/// ref:CGGGG  CCCCG  CCG CGG
/// sbr:C---G  C---G  C-G C-G
#[allow(unused)]
pub fn fill_deletion_info(del_bases_idx: &Vec<usize>, aligned_ref: &[u8], metric: &mut Metric) {
    if !del_bases_idx.is_empty() {
        let pre = get_pre_rbase(aligned_ref, del_bases_idx.first().copied().unwrap());
        let next = get_next_rbase(aligned_ref, del_bases_idx.last().copied().unwrap());

        let mut homo = false;
        let del_bases = del_bases_idx
            .iter()
            .map(|&idx| aligned_ref[idx])
            .collect::<Vec<_>>();
        if let Some(pre) = pre {
            homo |= del_bases.iter().filter(|&&b| b == pre).count() == del_bases.len();
        }
        if let Some(next) = next {
            homo |= del_bases.iter().filter(|&&b| b == next).count() == del_bases.len();
        }
        if homo {
            del_bases
                .iter()
                .for_each(|&v| metric.homodel_cnt[SEQ_NT4_TABLE[v as usize] as usize] += 1);
        } else {
            del_bases
                .iter()
                .for_each(|&v| metric.non_homodel_cnt[SEQ_NT4_TABLE[v as usize] as usize] += 1);
        }
    }
}

/// homo 判定逻辑：1. del 位置的 的 ref base 和其 左/右 ref-base 其一一致就是 homo del
/// ref:CGGGG  CCCCG  CCG CGG CGGAAC
/// sbr:C---G  C---G  C-G C-G C----C
#[allow(unused)]
pub fn fill_deletion_info_v2(del_bases_idx: &Vec<usize>, aligned_ref: &[u8], metric: &mut Metric) {
    if !del_bases_idx.is_empty() {

        del_bases_idx.iter().for_each(|&base_idx| {

            let pre_ref_base = get_pre_rbase(aligned_ref, base_idx);
            let next_ref_base = get_next_rbase(aligned_ref, base_idx);
            let cur_base = aligned_ref[base_idx];
            let mut homo = false;
            if let Some(pre_ref_base) = pre_ref_base {
                homo |= cur_base == pre_ref_base;
            }

            if let Some(next_ref_base) = next_ref_base {
                homo |= cur_base == next_ref_base;
            }

            if homo {
                metric.homodel_cnt[SEQ_NT4_TABLE[cur_base as usize] as usize] += 1;
            } else {
                metric.non_homodel_cnt[SEQ_NT4_TABLE[cur_base as usize] as usize] += 1;
            }

        });
    }
}


/// homo 判定逻辑：1. ins区域的base必须一致, 2. ins区域的base 和比对区域的上一个或者下一个base相同, 以下情况会被定义为homo
///     ref:C---G  C---G  C--G
///     sbr:CGGGG  CCCCG  ACCG
pub fn fill_insertion_info(
    ins_bases_idx: &Vec<usize>,
    aligned_ref: &[u8],
    aligned_qry: &[u8],
    metric: &mut Metric,
) {
    if !ins_bases_idx.is_empty() {
        let pre = get_pre_rbase(aligned_ref, ins_bases_idx.first().copied().unwrap());
        let next = get_next_rbase(aligned_ref, ins_bases_idx.last().copied().unwrap());
        let mut homo = false;
        let ins_bases = ins_bases_idx
            .iter()
            .map(|&idx| aligned_qry[idx])
            .collect::<Vec<_>>();
        if let Some(pre) = pre {
            homo |= ins_bases.iter().filter(|&&b| b == pre).count() == ins_bases.len();
        }
        if let Some(next) = next {
            homo |= ins_bases.iter().filter(|&&b| b == next).count() == ins_bases.len();
        }
        if homo {
            ins_bases
                .iter()
                .for_each(|&v| metric.homoins_cnt[SEQ_NT4_TABLE[v as usize] as usize] += 1);
        } else {
            ins_bases
                .iter()
                .for_each(|&v| metric.non_homoins_cnt[SEQ_NT4_TABLE[v as usize] as usize] += 1);
        }
    }
}

pub fn get_next_rbase(aligned_ref: &[u8], idx: usize) -> Option<u8> {
    let gap = '-' as u8;

    let mut idx_cursor_for_next = idx;

    let next = loop {
        if (idx_cursor_for_next + 1) >= aligned_ref.len() {
            break None;
        }

        idx_cursor_for_next += 1;
        if aligned_ref[idx_cursor_for_next] != gap {
            break Some(aligned_ref[idx_cursor_for_next]);
        }
    };
    next
}

pub fn get_pre_rbase(aligned_ref: &[u8], idx: usize) -> Option<u8> {
    let gap = '-' as u8;
    let mut idx_cursor_for_pre = idx;

    let pre = loop {
        if idx_cursor_for_pre == 0 {
            break None;
        }
        idx_cursor_for_pre -= 1;

        if aligned_ref[idx_cursor_for_pre] != gap {
            break Some(aligned_ref[idx_cursor_for_pre]);
        }
    };
    pre
}

#[cfg(test)]
mod test {
    use crate::align_processor::align_metric::{fill_metric_from_align_pair, Metric};

    #[test]
    fn test_fill_align_info_from_align_pair() {
        let aligned_ref = b"C-GGG-G";
        let aligned_qry = b"GG---CA";
        let mut metric = Metric::new("test".to_string(), 4, "test".to_string());
        fill_metric_from_align_pair(aligned_ref, aligned_qry, &mut metric);
        println!("{:?}", metric);
    }
}
