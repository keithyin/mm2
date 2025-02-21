use std::collections::HashSet;

use gskits::{
    gsbam::bam_record_ext::{BamRecord, BamRecordExt},
    utils::Range,
};
use minimap2::{Aligner, PresetSet};

pub trait TOverrideAlignerParam {
    fn modify_aligner(&self, aligner: &mut Aligner<PresetSet>);
}

#[derive(Default)]
pub struct InputFilterParams {
    pub np_range: Option<Range<i32>>,
    pub rq_range: Option<Range<f32>>,
}

impl InputFilterParams {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_np_range(mut self, np_range_str: &str) -> Self {
        self.np_range = Some(Range::<i32>::new(np_range_str));
        self
    }

    pub fn set_rq_range(mut self, rq_range_str: &str) -> Self {
        self.rq_range = Some(Range::<f32>::new(rq_range_str));
        self
    }

    pub fn valid(&self, record: &BamRecord) -> bool {
        let record_ext = BamRecordExt::new(record);

        if let Some(np_range_) = &self.np_range {
            if let Some(np) = record_ext.get_np() {
                if !np_range_.within_range(np as i32) {
                    return false;
                }
            }
        }

        if let Some(rq_range_) = &self.rq_range {
            if let Some(rq) = record_ext.get_rq() {
                if !rq_range_.within_range(rq) {
                    return false;
                }
            }
        }

        true
    }
}

#[derive(Debug, Clone, Default)]
pub struct AlignParams {
    pub matching_score: Option<i32>,
    pub mismatch_penalty: Option<i32>,
    pub gap_open_penalty: Option<String>,
    pub gap_extension_penalty: Option<String>,
}

impl AlignParams {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_m_score(mut self, m_score: i32) -> Self {
        self.matching_score = Some(m_score);
        self
    }

    pub fn set_mm_score(mut self, mm_score: i32) -> Self {
        self.mismatch_penalty = Some(mm_score);
        self
    }

    pub fn set_gap_open_penalty(mut self, go: String) -> Self {
        self.gap_open_penalty = Some(go);
        self
    }

    pub fn set_gap_extension_penalty(mut self, ge: String) -> Self {
        self.gap_extension_penalty = Some(ge);
        self
    }
}

impl TOverrideAlignerParam for AlignParams {
    fn modify_aligner(&self, aligner: &mut Aligner<PresetSet>) {
        if let Some(m) = self.matching_score {
            aligner.mapopt.a = m;
        }

        if let Some(mm) = self.mismatch_penalty {
            aligner.mapopt.b = mm;
        }

        if let Some(gap_o) = &self.gap_open_penalty {
            if gap_o.contains(",") {
                let (o1, o2) = gap_o.rsplit_once(",").unwrap();
                aligner.mapopt.q = o1.parse::<i32>().unwrap();
                aligner.mapopt.q2 = o2.parse::<i32>().unwrap();
            } else {
                aligner.mapopt.q = gap_o.parse::<i32>().unwrap();
            }
        }

        if let Some(gap_e) = &self.gap_extension_penalty {
            if gap_e.contains(",") {
                let (e1, e2) = gap_e.rsplit_once(",").unwrap();
                aligner.mapopt.e = e1.parse::<i32>().unwrap();
                aligner.mapopt.e2 = e2.parse::<i32>().unwrap();
            } else {
                aligner.mapopt.e = gap_e.parse::<i32>().unwrap();
            }
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct MapParams {}

impl TOverrideAlignerParam for MapParams {
    fn modify_aligner(&self, _aligner: &mut Aligner<PresetSet>) {}
}

#[derive(Debug, Clone, Copy, Default)]
pub struct IndexParams {
    pub kmer: Option<usize>,
    pub wins: Option<usize>,
}

impl IndexParams {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_kmer(mut self, kmer: usize) -> Self {
        self.kmer = Some(kmer);
        self
    }

    pub fn set_wins(mut self, wins: usize) -> Self {
        self.wins = Some(wins);
        self
    }
}

impl TOverrideAlignerParam for IndexParams {
    fn modify_aligner(&self, aligner: &mut Aligner<PresetSet>) {
        if let Some(k) = self.kmer {
            aligner.idxopt.k = k as i16;
        }
        if let Some(w) = self.wins {
            aligner.idxopt.w = w as i16;
        }
    }
}

#[derive(Debug, Clone)]
pub struct OupParams {
    pub discard_secondary: bool,
    pub discard_supplementary: bool,
    pub oup_identity_threshold: f32,
    pub oup_coverage_threshold: f32,
    pub discard_multi_align_reads: bool,
    pub pass_through_tags: HashSet<String>, 
}

impl Default for OupParams {
    fn default() -> Self {
        Self {
            discard_secondary: false,
            discard_supplementary: false,
            oup_identity_threshold: -1.0,
            oup_coverage_threshold: -1.0,
            discard_multi_align_reads: false,
            pass_through_tags: HashSet::new(),
        }
    }
}

impl OupParams {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_discard_secondary(mut self, discard_secondary: bool) -> Self {
        if discard_secondary {
            assert!(self.discard_multi_align_reads == false);
        }
        self.discard_secondary = discard_secondary;
        self
    }

    pub fn set_discard_supplementary(mut self, discard_supplementary: bool) -> Self {
        if discard_supplementary {
            assert!(self.discard_multi_align_reads == false);
        }
        self.discard_supplementary = discard_supplementary;
        self
    }
    pub fn set_oup_identity_threshold(mut self, oup_identity_threshold: f32) -> Self {
        self.oup_identity_threshold = oup_identity_threshold;
        self
    }
    pub fn set_oup_coverage_threshold(mut self, oup_coverage_threshold: f32) -> Self {
        self.oup_coverage_threshold = oup_coverage_threshold;
        self
    }
    pub fn set_discard_multi_align_reads(mut self, discard_multi_align_reads: bool) -> Self {
        if discard_multi_align_reads {
            assert!(self.discard_secondary == false);
            assert!(self.discard_supplementary == false);
        }

        self.discard_multi_align_reads = discard_multi_align_reads;
        self
    }

    pub fn set_pass_through_tags(mut self, tags: Option<&String>) -> Self {
        if let Some(tags) = tags {
            self.pass_through_tags = tags.trim()
                .split(",")
                .map(|v| v.to_string())
                .collect::<HashSet<_>>();
        }

        self
    }

    pub fn valid(&self, record: &BamRecord) -> bool {
        if self.discard_secondary && record.is_secondary() {
            return false;
        }

        if self.discard_supplementary && record.is_supplementary() {
            return false;
        }

        let record_ext = BamRecordExt::new(&record);
        let iy = record_ext.compute_identity();
        let ec = record_ext.compute_query_coverage();

        if iy < self.oup_identity_threshold {
            return false;
        }

        if ec < self.oup_coverage_threshold {
            return false;
        }

        true
    }
}

impl TOverrideAlignerParam for OupParams {
    fn modify_aligner(&self, aligner: &mut Aligner<PresetSet>) {
        if self.discard_secondary {
            aligner.mapopt.flag &= !0x1000000000;
        }
    }
}
