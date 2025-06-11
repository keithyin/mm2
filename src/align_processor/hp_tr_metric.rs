use std::{fmt::Display, rc::Rc, sync::Arc};

use hp_tr_finder::intervaltree;
use rayon::vec;

use crate::align_processor::{AlignInfo, TMetric};

pub struct HpTrMetric {
    qname: String,
    tname: Arc<String>,
    region2motif: Arc<intervaltree::IntervalTree<usize, Arc<String>>>,
    align_infos: Vec<AlignInfo>,
}

impl HpTrMetric {
    pub fn new(
        qname: String,
        tname: Arc<String>,
        region2motif: Arc<intervaltree::IntervalTree<usize, Arc<String>>>,
    ) -> Self {
        Self {
            qname,
            tname,
            region2motif,
            align_infos: vec![],
        }
    }
}

impl Display for HpTrMetric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "")
    }
}

impl TMetric for HpTrMetric {
    fn add_align_info(&mut self, align_info: AlignInfo) {
        self.align_infos.push(align_info);
    }

    fn compute_metric(&mut self, target_seq: &str, query_seq: &str) {

        

    }
}
