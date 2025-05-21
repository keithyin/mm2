use std::fmt::Display;

use super::{AlignInfo, TMetric};
use lazy_static::lazy_static;
#[derive(Debug, Clone, Copy, Default)]
pub struct TimeSpan {
    pub ins: usize,
    pub del: usize,
    pub eq: usize,
    pub diff: usize,
}

impl Display for TimeSpan {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}\t{}\t{}\t{}", self.eq, self.diff, self.ins, self.del)
    }
}

pub struct TimeErrMetric {
    qname: String,
    rname: String,
    base_secs: Vec<f64>,
    align_infos: Vec<AlignInfo>,
    time_spans: Vec<TimeSpan>,
    tot_time: usize, // minutes,
    interval: usize, //secs
}

impl TimeErrMetric {
    pub fn new(qname: String, rname: String, base_secs: Vec<f64>, interval: Option<usize>) -> Self {
        let tot_minutes = if let Some(tot_secs) = base_secs.last() {
            (*tot_secs / 60.).ceil()
        } else {
            0.0
        };
        Self {
            qname,
            rname,
            base_secs,
            align_infos: vec![],
            time_spans: vec![],
            tot_time: tot_minutes as usize,
            interval: interval.unwrap_or(600),
        }
    }

    fn secs2buckets(&self, secs: f64) -> usize {
        ((secs - 1.0).max(0.0) / self.interval as f64).floor() as usize
    }
}

impl Display for TimeErrMetric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut final_res_items = vec![];
        self.time_spans.iter().enumerate().for_each(|(idx, tspan)| {
            let mut result_items = vec![];

            result_items.push(self.qname.clone());
            result_items.push(self.rname.clone());
            result_items.push(format!("{}", idx * self.interval / 60));
            result_items.push(format!("{}", tspan));

            final_res_items.push(result_items.join("\t"));
        });

        if final_res_items.len() == 0 {
            final_res_items = vec![format!(
                "{}\t{}\t0\t{}",
                self.qname.clone(),
                self.rname.clone(),
                TimeSpan::default()
            )];
        }

        write!(f, "{}", final_res_items.join("\n"))
    }
}

lazy_static! {
    pub static ref METRIC_CSV_HEADER: Vec<String> = {
        let csv_header = vec![
            "qname".to_string(),
            "rname".to_string(),
            "minutes".to_string(),
            "eq".to_string(),
            "diff".to_string(),
            "ins".to_string(),
            "del".to_string(),
        ];
        csv_header
    };
}

/*

0 => Cigar::Match(cnt),
1 => Cigar::Ins(cnt),
2 => Cigar::Del(cnt),
3 => Cigar::RefSkip(cnt),
4 => Cigar::SoftClip(cnt),
5 => Cigar::HardClip(cnt),
6 => Cigar::Pad(cnt),
7 => Cigar::Equal(cnt),
8 => Cigar::Diff(cnt),
*/
impl TMetric for TimeErrMetric {
    fn add_align_info(&mut self, align_info: super::AlignInfo) {
        self.align_infos.push(align_info);
    }

    fn compute_metric(&mut self, _target_seq: &str, _query_seq: &str) {
        let num_buckets = (self.tot_time as f64 * 60. / self.interval as f64).ceil() as usize;
        let mut time_spans = vec![TimeSpan::default(); num_buckets];

        self.align_infos.iter().for_each(|align_info| {
            let mut qcursor = align_info.qstart;
            let cigar_iter: Box<dyn Iterator<Item = &(u32, u8)>> = if align_info.fwd {
                Box::new(align_info.cigar.iter())
            } else {
                Box::new(align_info.cigar.iter().rev())
            };

            cigar_iter.for_each(|&(cnt, op)| {
                // println!("{}, {}", cnt, op);
                let cnt = cnt as usize;
                (0..cnt).into_iter().for_each(|idx| {
                    let base_idx = match op {
                        1 | 7 | 8 => qcursor + idx,
                        _ => qcursor,
                    };

                    // println!("base_idx:{}", base_idx);

                    let bucket_idx = self.secs2buckets(self.base_secs[base_idx]);
                    let time_span = &mut time_spans[bucket_idx];
                    match op {
                        1 => time_span.ins += 1,
                        2 => time_span.del += 1,
                        7 => time_span.eq += 1,
                        8 => time_span.diff += 1,
                        0 => panic!("eqx needed"),
                        _ => (),
                    }
                });
                match op {
                    1 | 7 | 8 => qcursor += cnt,
                    _ => (),
                }
            });
        });

        self.time_spans = time_spans;
    }
}

#[cfg(test)]
mod test {
    use crate::align_processor::{AlignInfo, TMetric};

    use super::TimeErrMetric;

    #[test]
    fn test_time_err_metric() {
        let base_secs = vec![200, 400, 600, 800, 1000, 1200, 1400, 1600]
            .into_iter()
            .map(|v| v as f64)
            .collect();
        let mut metric = TimeErrMetric::new(
            "query".to_string(),
            "target".to_string(),
            base_secs,
            Some(600),
        );

        let align_info = AlignInfo {
            rstart: 0,
            rend: 100,
            qstart: 0,
            qend: 4,
            fwd: true,
            primary: true,
            cigar: vec![(1, 7), (1, 8), (1, 1), (1, 2)],
        };
        metric.add_align_info(align_info);

        let align_info = AlignInfo {
            rstart: 0,
            rend: 100,
            qstart: 3,
            qend: 8,
            fwd: false,
            primary: true,
            cigar: vec![(2, 7), (2, 8), (1, 1), (1, 2)],
        };
        metric.add_align_info(align_info);

        metric.compute_metric("", "");
        println!("{}", metric);
    }
}
