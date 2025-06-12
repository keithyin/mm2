use crate::{
    align_processor::{AlignInfo, TMetric},
    aligned_pairs::{AlignOp, TAlignedPairs},
};
use clap::builder::Str;
use lazy_static::lazy_static;
use std::{collections::HashMap, fmt::Display, sync::Arc};

lazy_static! {
    pub static ref METRIC_CSV_HEADER: Vec<String> = {
        let csv_header = vec![
            "qname".to_string(),
            "rname".to_string(),
            "motif".to_string(),
            "eq".to_string(),
            "diff".to_string(),
            "ins".to_string(),
            "del".to_string(),
        ];
        csv_header
    };
}

#[derive(Debug, Default)]
pub struct AlignCounter {
    pub eq: usize,
    pub diff: usize,
    pub ins: usize,
    pub del: usize,
}

impl AlignCounter {
    pub fn update(&mut self, align_op: AlignOp) {
        match align_op {
            AlignOp::Equal(n) => self.eq += n as usize,
            AlignOp::Diff(n) => self.diff += n as usize,
            AlignOp::Ins(n) => self.ins += n as usize,
            AlignOp::Del(n) => self.del += n as usize,
            _ => panic!("not a valid align op: {:?}", align_op),
        }
    }
}

pub struct HpTrMetric {
    qname: String,
    tname: Arc<String>,
    region2motif: Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>,
    align_infos: Vec<AlignInfo>,
    //Arc<String> for motif string. [fwd_counter, rev_counter]
    metric_core: HashMap<Arc<String>, [AlignCounter; 2]>,
}

impl HpTrMetric {
    pub fn new(
        qname: String,
        tname: Arc<String>,
        region2motif: Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>,
    ) -> Self {
        Self {
            qname,
            tname,
            region2motif,
            align_infos: vec![],
            metric_core: HashMap::new(),
        }
    }
}

impl Display for HpTrMetric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut result_items = vec![];
        for (key, fwd_rev_counters) in &self.metric_core {
            for (idx, counter) in fwd_rev_counters.iter().enumerate() {
                let mut innner_items = vec![];
                innner_items.push(self.qname.clone());
                innner_items.push(self.tname.as_ref().clone());
                innner_items.push(if idx == 0 {
                    key.as_ref().clone()
                } else {
                    reverse_complement_motif(&key)
                });

                innner_items.push(counter.eq.to_string());
                innner_items.push(counter.diff.to_string());
                innner_items.push(counter.ins.to_string());
                innner_items.push(counter.del.to_string());

                result_items.push(innner_items.join("\t"));
            }
        }
        write!(f, "{}", result_items.join("\n"))
    }
}

impl TMetric for HpTrMetric {
    fn add_align_info(&mut self, align_info: AlignInfo) {
        self.align_infos.push(align_info);
    }

    fn compute_metric(&mut self, _target_seq: &str, query_seq: &str) {
        for align_info in &self.align_infos {
            let mut pre_ins_cnt = 0;
            let target_start = align_info.target_start();
            let target_end = align_info.target_end();
            let strand_idx = if align_info.is_reverse() { 1 } else { 0 };

            for (_qpos, rpos, align_op) in align_info.aligned_pairs(query_seq.len()) {
                if let Some(rpos) = rpos {
                    if let Some(motifs) = self.region2motif.get(&rpos) {
                        motifs
                            .iter()
                            .filter(|motif| motif.0 .0 >= target_start && motif.0 .1 <= target_end)
                            .for_each(|motif| {
                                self.metric_core.entry(motif.1.clone()).or_default()[strand_idx]
                                    .update(align_op);

                                self.metric_core.entry(motif.1.clone()).or_default()[strand_idx]
                                    .update(AlignOp::Ins(pre_ins_cnt));
                            });
                    }
                }

                if rpos.is_none() {
                    pre_ins_cnt += 1;
                } else {
                    pre_ins_cnt = 0;
                }
            }
        }
    }
}

fn reverse_complement(dna: &str) -> String {
    dna.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c,
        })
        .collect()
}

fn reverse_complement_motif(s: &str) -> String {
    // 用正则提取括号里的序列
    let re = regex::Regex::new(r"\((?P<seq>[ACGT]+)\)").unwrap();

    re.replace_all(s, |caps: &regex::Captures| {
        let seq = &caps["seq"];
        let rc = reverse_complement(seq);
        format!("({})", rc)
    })
    .into_owned()
}
