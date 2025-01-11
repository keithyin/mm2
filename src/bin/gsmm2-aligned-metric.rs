use std::collections::HashMap;
use std::fmt::Display;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::ops::{Deref, DerefMut};
use std::{path, thread};

use bio::bio_types::strand::ReqStrand;
use crossbeam::channel::{Receiver, Sender};
use gskits::dna::{reverse_complement, SEQ_NT4_TABLE};
use gskits::fastx_reader::fasta_reader::FastaFileReader;
use gskits::fastx_reader::read_fastx;
use gskits::gsbam::bam_record_ext::{BamRecord, BamRecordExt};
use gskits::pbar::{self, DEFAULT_INTERVAL};
use minimap2::{Mapping, Strand};
use mm2::params::{AlignParams, IndexParams, InputFilterParams, MapParams, OupParams};
use mm2::{align_single_query_to_targets, convert_mapping_cigar_to_record_cigar, NoMemLeakAligner};
use mm2::{build_aligner, query_seq_sender};
use rust_htslib::bam::ext::BamRecordExtensions;

use std::str::FromStr;

use clap::{self, Args, Parser};

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(long="preset", default_value_t=String::from_str("map-ont").unwrap(), 
        help="read https://lh3.github.io/minimap2/minimap2.html for more details")]
    pub preset: String,

    #[command(flatten)]
    pub metric_args: MetricArgs,
}

#[derive(Debug, Args, Clone, Copy, Default)]
pub struct IndexArgs {
    #[arg(long, help = "minimizer kmer")]
    pub kmer: Option<usize>,

    #[arg(long, help = "minimizer window size")]
    pub wins: Option<usize>,
}

impl IndexArgs {
    fn to_index_params(&self) -> IndexParams {
        let mut param = IndexParams::new();
        param = if let Some(kmer) = self.kmer {
            param.set_kmer(kmer)
        } else {
            param
        };

        param = if let Some(wins) = self.wins {
            param.set_wins(wins)
        } else {
            param
        };

        param
    }
}

#[derive(Debug, Args, Clone)]
pub struct MetricArgs {
    #[command(flatten)]
    pub io_args: IoArgs,

    #[command(flatten)]
    pub index_args: IndexArgs,

    #[command(flatten)]
    pub map_args: MapArgs,

    #[command(flatten)]
    pub align_args: AlignArgs,

    #[command(flatten)]
    pub oup_args: OupArgs,
}

impl MetricArgs {
    pub fn get_oup_file(&self) -> String {
        if self.oup_args.discard_supplementary {
            return format!("{}/metric_noSupp.csv", self.io_args.get_oup_dir());
        }

        if self.oup_args.discard_multi_mapping_reads {
            return format!("{}/metric_noMar.csv", self.io_args.get_oup_dir());
        }

        format!("{}/metric.csv", self.io_args.get_oup_dir())
    }
}

#[derive(Debug, Args, Clone)]
pub struct IoArgs {
    #[arg(
        short = 'q',
        help = "query file paths, 
    if multiple query_filepath are provided, 
    their query_names will be rewriten to ___0, ___1, and so on, 
    based on the order of the filenames. valid file format bam/fa/fq",
        required = true
    )]
    pub query: Vec<String>,

    /// target.fasta
    #[arg(long = "target", short = 't')]
    pub target: String,
    #[arg(long = "oupdir", help = "oupdir")]
    pub oupdir: Option<String>,

    #[arg(
        long = "np-range",
        help = "1:3,5,7:9 means [[1, 3], [5, 5], [7, 9]]. only valid for bam input that contains np field"
    )]
    pub np_range: Option<String>,

    #[arg(
        long = "rq-range",
        help = "0.9:1.1 means 0.9<=rq<=1.1. only valid for bam input that contains rq field"
    )]
    pub rq_range: Option<String>,
}

impl IoArgs {
    pub fn get_oup_dir(&self) -> String {
        let oupdir = if let Some(oup) = &self.oupdir {
            oup.to_string()
        } else {
            self.query[0].rsplit_once(".").unwrap().0.to_string()
        };

        if !path::Path::new(&oupdir).exists() {
            fs::create_dir_all(&oupdir).unwrap()
        }
        oupdir
    }

    pub fn to_input_filter_params(&self) -> InputFilterParams {
        let mut param = InputFilterParams::new();
        param = if let Some(ref np_range_str) = self.np_range {
            param.set_np_range(np_range_str)
        } else {
            param
        };

        param = if let Some(ref rq_range_str) = self.rq_range {
            param.set_rq_range(rq_range_str)
        } else {
            param
        };

        param
    }
}

#[derive(Debug, Args, Clone, Default)]
pub struct MapArgs {}

impl MapArgs {
    pub fn to_map_params(&self) -> MapParams {
        MapParams::default()
    }
}

#[derive(Debug, Args, Clone, Default)]
pub struct AlignArgs {
    #[arg(short = 'm', help = "matching_score>=0, recommend 2")]
    matching_score: Option<i32>,

    #[arg(short = 'M', help = "mismatch_penalty >=0, recommend 4")]
    mismatch_penalty: Option<i32>,

    #[arg(short = 'o', help = "gap_open_penalty >=0, recommend 4,24")]
    gap_open_penalty: Option<String>,

    #[arg(short = 'e', help = "gap_extension_penalty >=0, recommend 2,1")]
    gap_extension_penalty: Option<String>,
}

impl AlignArgs {
    pub fn to_align_params(&self) -> AlignParams {
        let mut param = AlignParams::new();
        param = if let Some(ms) = self.matching_score {
            param.set_m_score(ms)
        } else {
            param
        };

        param = if let Some(mms) = self.mismatch_penalty {
            param.set_mm_score(mms)
        } else {
            param
        };

        param = if let Some(ref go) = self.gap_open_penalty {
            param.set_gap_open_penalty(go.to_string())
        } else {
            param
        };

        param = if let Some(ref ge) = self.gap_extension_penalty {
            param.set_gap_extension_penalty(ge.to_string())
        } else {
            param
        };

        param
    }
}

#[derive(Debug, Args, Clone, Default)]
pub struct OupArgs {
    #[arg(
        long = "noSupp",
        help = "discard supplementary alignment, only keep the primary alignemnt for metric"
    )]
    pub discard_supplementary: bool,

    #[arg(
        long = "noMar",
        help = "discard multi-mapping reads, only keep the unique mapping reads for metric"
    )]
    pub discard_multi_mapping_reads: bool,
}

impl OupArgs {
    fn to_oup_params(&self) -> OupParams {
        let mut param = OupParams::new();
        param = param
            .set_discard_secondary(true)
            .set_discard_supplementary(self.discard_supplementary)
            .set_oup_identity_threshold(-1.0)
            .set_oup_coverage_threshold(-1.0)
            .set_discard_multi_align_reads(self.discard_multi_mapping_reads);
        param
    }
}

fn metric_entrance(preset: &str, tot_threads: Option<usize>, args: &MetricArgs) {
    let tot_threads = tot_threads.unwrap_or(num_cpus::get());
    assert!(tot_threads >= 10, "at least 10 threads are needed");

    let target_filename = &args.io_args.target;

    let fa_iter = FastaFileReader::new(target_filename.to_string());
    let targets = read_fastx(fa_iter);
    let targetname2seq = targets
        .iter()
        .map(|v| (v.name.clone(), v.seq.clone()))
        .collect::<HashMap<_, _>>();

    let index_params = args.index_args.to_index_params();
    let map_params = args.map_args.to_map_params();
    let align_params = args.align_args.to_align_params();
    let oup_params = args.oup_args.to_oup_params();
    let inp_filter_params = args.io_args.to_input_filter_params();

    let aligners = build_aligner(
        preset,
        &index_params,
        &map_params,
        &align_params,
        &oup_params,
        &targets,
    );

    thread::scope(|s| {
        let aligners = &aligners;
        let targetname2seq = &targetname2seq;
        let inp_filter_params = &inp_filter_params;
        let (qs_sender, qs_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            query_seq_sender(&args.io_args.query, qs_sender, inp_filter_params);
        });

        let align_threads = tot_threads - 4;
        let (metric_sender, metric_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..align_threads {
            let qs_recv_ = qs_recv.clone();
            let metric_sender_ = metric_sender.clone();
            s.spawn(move || {
                metric_worker(
                    qs_recv_,
                    metric_sender_,
                    aligners,
                    targetname2seq,
                    &oup_params,
                )
            });
        }
        drop(qs_recv);
        drop(metric_sender);

        let oup_filename = args.get_oup_file();
        dump_metric_worker(metric_recv, &oup_filename, true);
    });
}

fn main() {
    let args = Cli::parse();

    let align_threads = args.threads.clone();

    let time_fmt = time::format_description::parse(
        "[year]-[month padding:zero]-[day padding:zero]#[hour]:[minute]:[second]",
    )
    .unwrap();

    let time_offset =
        time::UtcOffset::current_local_offset().unwrap_or_else(|_| time::UtcOffset::UTC);
    let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt.clone());
    tracing_subscriber::fmt::fmt().with_timer(timer).init();

    metric_entrance(&args.preset, align_threads, &args.metric_args);
}

pub fn metric_worker(
    query_record_recv: Receiver<gskits::ds::ReadInfo>,
    metric_sender: Sender<Metric>,
    aligners: &Vec<NoMemLeakAligner>,
    targetname2seq: &HashMap<String, String>,
    oup_params: &OupParams,
) {
    for query_record in query_record_recv {
        let metric = compute_metric(&query_record, aligners, targetname2seq, oup_params);
        metric_sender.send(metric).unwrap();
    }
}

pub fn dump_metric_worker(metric_recv: Receiver<Metric>, fname: &str, enable_pb: bool) {
    let mut writer =
        BufWriter::new(File::create(fname).expect(&format!("create file error: {}", fname)));
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

    writeln!(&mut writer, "{}", csv_header.join("\t")).unwrap();

    let pb = if enable_pb {
        Some(pbar::get_spin_pb(
            format!("gsmm2-aligned-metric: writing metric to {}", fname),
            DEFAULT_INTERVAL,
        ))
    } else {
        None
    };

    for metric in metric_recv {
        pb.as_ref().map(|v| v.inc(1));
        writeln!(&mut writer, "{}", metric).unwrap();
    }

    pb.as_ref().map(|v| v.finish());
}

pub fn compute_metric(
    query_record: &gskits::ds::ReadInfo,
    aligners: &Vec<NoMemLeakAligner>,
    targetname2seq: &HashMap<String, String>,
    oup_params: &OupParams,
) -> Metric {
    let hits = align_single_query_to_targets(&query_record, aligners);

    if hits.is_empty() || (hits.len() > 0 && oup_params.discard_multi_align_reads) {
        return Metric::new(
            query_record.name.clone(),
            query_record.seq.len(),
            "".to_string(),
        );
    }

    let hits = hits
        .into_iter()
        .filter(|v| v.is_primary || (v.is_supplementary && !oup_params.discard_supplementary))
        .collect::<Vec<_>>();

    let target_name = hits[0].target_name.as_ref().unwrap().as_ref().clone();
    let target_seq = targetname2seq.get(&target_name).unwrap();
    let mut metric = Metric::new(
        query_record.name.clone(),
        query_record.seq.len(),
        if hits.is_empty() {
            "".to_string()
        } else {
            target_name.clone()
        },
    );

    hits.into_iter()
        .for_each(|hit| metric.add_align_info(hit.into()));

    metric.compute_metric(target_seq, &query_record.seq);

    metric
}

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
                fill_align_info_from_align_pair(
                    align_pair.0.as_bytes(),
                    align_pair.1.as_bytes(),
                    self,
                );
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

/*

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
        "homoDel".to_string()
    ];
    for base in ["A", "C", "G", "T"] {
        csv_header.push(format!("match-{}", base));
        csv_header.push(format!("misMatch-{}", base));
        csv_header.push(format!("ins-{}", base));
        csv_header.push(format!("homoIns-{}", base));
        csv_header.push(format!("del-{}", base));
        csv_header.push(format!("homoDel-{}", base));
    }
*/
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

/// 从 Mapping 中直接抽取出来的信息，使用该类作为一个中间存储
#[derive(Debug)]
pub struct AlignInfo {
    rstart: usize,
    rend: usize,
    qstart: usize,
    qend: usize,
    fwd: bool,
    primary: bool,
    cigar: Vec<(u32, u8)>,
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
            reverse_complement(&target_seq[self.rstart..self.rend])
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

#[derive(Debug)]
pub struct SingleQueryAlignInfo(Vec<AlignInfo>);
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

pub struct TseqAndRecord {
    ori_rstart: usize, // 用来存储 当前比对片段 对应 原始比对区域的 rstart
    ori_rend: usize,   // 用来存储 当前比对片段 对应 原始比对区域的 rend
    ori_qstart: usize, // 用来存储 当前比对片段 对应 原始比对区域的 qstart
    ori_qend: usize,   // 用来存储 当前比对片段 对应 原始比对区域的 qend
    tseq: String,      // 仅保留了比对区域的 reference sequence
    record: BamRecord, // 仅保留了 比对区域的 比对信息。等价于将 比对区域的 ref 和 query 都截取下来，然后重新比对之后的结果
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

pub fn fill_align_info_from_align_pair(
    aligned_ref: &[u8],
    aligned_qry: &[u8],
    metric: &mut Metric,
) {
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
                fill_deletion_info(&del_bases_idx, aligned_ref, metric);
                del_bases_idx.clear();
            }

            (qbase, rbase) => {
                if rbase == qbase {
                    metric.matched_cnt[SEQ_NT4_TABLE[qbase as usize] as usize] += 1;
                } else {
                    metric.mismatched_cnt[SEQ_NT4_TABLE[qbase as usize] as usize] += 1;
                }

                fill_deletion_info(&del_bases_idx, aligned_ref, metric);
                fill_insertion_info(&ins_bases_idx, aligned_ref, aligned_qry, metric);

                del_bases_idx.clear();
                ins_bases_idx.clear();
            }
        }
    }

    fill_deletion_info(&del_bases_idx, aligned_ref, metric);
    fill_insertion_info(&ins_bases_idx, aligned_ref, aligned_qry, metric);
}

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

#[derive(Debug, Clone)]
struct Region {
    start: usize,
    end: usize,
}

impl Region {
    fn length(&self) -> usize {
        self.end - self.start
    }
}

struct Regions(Vec<Region>);

impl Deref for Regions {
    type Target = Vec<Region>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Regions {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<&Vec<(usize, usize)>> for Regions {
    fn from(value: &Vec<(usize, usize)>) -> Self {
        let regions = value
            .iter()
            .map(|v| Region {
                start: v.0,
                end: v.1,
            })
            .collect::<Vec<_>>();
        Self(regions)
    }
}

impl Regions {
    pub fn total_length(&self) -> usize {
        self.iter().map(|v| v.length()).sum()
    }

    pub fn merge_regions(&self) -> Self {
        let mut regions = self.0.clone();
        regions.sort_by_key(|v| v.start);
        let mut merged_regions: Vec<Region> = Vec::new();

        for region in regions {
            if let Some(last) = merged_regions.last_mut() {
                if region.start <= last.end {
                    last.end = last.end.max(region.end);
                } else {
                    merged_regions.push(region);
                }
            } else {
                merged_regions.push(region);
            }
        }
        Self(merged_regions)
    }

    pub fn ovlp_length(&self) -> usize {
        let mut total_overlap = 0;
        let mut events = Vec::new();

        for region in &self.0 {
            events.push((region.start, 1));
            events.push((region.end, -1));
        }

        events.sort();

        let mut current_depth = 0;
        let mut prev_pos = events[0].0;

        for (pos, change) in events {
            if current_depth > 1 {
                total_overlap += pos - prev_pos;
            }
            current_depth += change;
            prev_pos = pos;
        }

        total_overlap
    }

    pub fn ovlp_ratio(&self) -> f32 {
        let merged_regions = self.merge_regions();
        let len_after_merge = merged_regions.total_length();
        let ovlp_len = self.ovlp_length();
        if len_after_merge == 0 {
            0.0
        } else {
            ovlp_len as f32 / len_after_merge as f32
        }
    }

    /// the region must be sorted for use this function
    pub fn gaps(&self, init_pos: Option<usize>, end_pos: Option<usize>) -> Vec<i64> {
        if self.is_empty() {
            return vec![];
        }
        let mut gaps = vec![];

        let mut pre_pos = init_pos;
        self.iter().for_each(|region| {
            if let Some(pre_pos) = pre_pos {
                let gap = region.start as i64 - pre_pos as i64;
                gaps.push(gap);
            }
            pre_pos = Some(region.end);
        });

        if let Some(end_pos) = end_pos {
            let gap = end_pos as i64 - pre_pos.unwrap() as i64;
            gaps.push(gap);
        }

        gaps
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashMap;

    use gskits::{
        ds::ReadInfo,
        fastx_reader::{fasta_reader::FastaFileReader, read_fastx},
    };
    use mm2::{
        build_aligner,
        params::{AlignParams, IndexParams, MapParams, OupParams},
    };
    use rust_htslib::bcf::synced::pairing::SOME;

    use crate::{compute_metric, Regions};

    #[test]
    fn test_compute_metric() {
        let ref_file = "/data/ccs_data/MG1655.fa";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);

        let targetname2seq = targets
            .iter()
            .map(|v| (v.name.clone(), v.seq.clone()))
            .collect::<HashMap<String, String>>();

        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
        );

        //fwd
        let seq = b"GAGAGAAGAGATGTTCTACTGGTGTGCCGATGGTGGCGGATTAGTCCTCGCTCAATCTGGGGTCAGGCGGTGATGGTCTATTGCTATCAATTATACAACAATTAAACAACAACCGGCGAAAGTGATGCAACGGCAGACCAACATCACTGCAAGCTTTACGCGAACGAGCCAGACATTGCTGACGACTCTGGCAGTGCATAGATGACAATAAGTTCGACTGGTTACAAAACCCTTGGGGCTTTTAGAGCAACGAGACACGGCAATGTTGCACCGTTTGCTGCATGAATTGAATTGAACAAAAATATCACCAATAAAAAACGCCTTATGTAAATTTTTCAGCTTTTCATTCTGACTGCAACGGGGCATATGTCTCTGTGTAGATTAAAAAAAGAGTGTCTGATGCAGTTCTGAACTGGTTACCTGCCGTGAGTATAATTAATAAAATTTTATTGACTTGGTCACTAAATACTTTAACCAATATAGGCATAGCGCGACAGACAGATAAAATTACAGAGTACACAACCAATCCCATGGAAACGCATTAGCCCACCATTACCCCACCATCACCATTACCAACAGGTAAAAACGGTGCGGGCTGATCGCGCTATCAGGAAACACAGAAAAAGCCCGCACCTGACATGCGGGCTTTTTTTTTTCGACCAAAGGTAATATACGAGGTAACAACCATGCGAGTGTTTGAAGTTCGGCGGGTTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTTGTGCCGATATTCTGGAAAGCAATGCAGCAGGGGCAGGTGGCCACCGTCTCCTCTGACCCCGCCAAAATCACCACACCTGGTGCGATGATTGAAAAACCATTAGCGGCAAGGAGCTTTAAACCCAATATCAAGCGATGCCGAACGTATTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCAGCCGGGTTTCACCGCTGGCCAATTGAAAACTTTCGTTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACCTGCGCTGATTTGCCGTGCGAAGAAATGTCGATCAGCCATTATATCTCTC";

        // rev
        let seq = b"TATTCGATTGGAAATGCCCTTGACCAGGTAATTCAGTCTCATCACGAACCATGAGCGTACCTGGTGCTTGAGGATTTCCGGATATATTTTAATCAGGCAAGGGACGATCTGAACTGGCGGATGGGGGTAAATGGTGGCGGGGTTGAAGAACTTTAGCGCGAAGTAGGAAAGCTCCTCGCTTCGTTGCTCTAAAAAAGCCCCCAGGCGTTGTTTGTACCAGTCGACAAGTTTTAAATGTCATCTGCCACGCCAGAGTCGTCAGCAATGTCATGGCTCGTTCGCGTAAAAGCTTTGACAGTTGATGTTGGTCTGCCCGTTGCATCACTTTTCGGCCGGTTGTTGTATTAAGTTGCTAAATTGATAAGCAATAGATCACCGCCTGCCCAGATTGAGCGAAGGTAATCCGCCACCATCGGCAACACAGTAATGACGTCAGCCAGCCAAACGCTAACTCTTCGTTTAGTCAACCCGGAATCCTTCGCGACCCACCCAGCGCGGCATGGCCATCCATGAAAGATTTTTCCTCTAACAGCGCACCAGTTCAACTGGCGTGGCGTAAGTCAATGATATTTCGCCCGACTGCGCGCAGTGGTGGCGACAAGTGAAATCGACATCGTGTAACGATTCA";
        let query_record =
            ReadInfo::new_fa_record("tt".to_string(), String::from_utf8(seq.to_vec()).unwrap());
        compute_metric(
            &query_record,
            &aligners,
            &targetname2seq,
            &OupParams::default(),
        );
    }

    #[test]
    fn test_calcute_overlap_metrics() {
        let regions = vec![(1_usize, 5_usize), (1, 5)];
        let regions: Regions = (&regions).into();
        let res = regions.ovlp_ratio();
        println!("{:?}", res);

        let regions = vec![(1_usize, 5_usize), (5, 10)];
        let regions: Regions = (&regions).into();
        let res = regions.ovlp_ratio();
        println!("{:?}", res);
    }

    #[test]
    fn test_calcute_gaps() {
        let regions = vec![(1_usize, 5_usize), (1, 5)];
        let regions: Regions = (&regions).into();
        let res = regions.gaps(Some(0), Some(100));
        println!("{:?}", res);

        let regions = vec![(1_usize, 5_usize), (5, 10)];
        let regions: Regions = (&regions).into();
        let res = regions.gaps(Some(0), Some(100));

        println!("{:?}", res);

        let regions = vec![(1_usize, 5_usize), (5, 10)];
        let regions: Regions = (&regions).into();
        let res = regions.gaps(None, None);

        println!("{:?}", res);

        let regions = vec![(1_usize, 5_usize), (7, 10)];
        let regions: Regions = (&regions).into();
        let res = regions.gaps(None, None);

        println!("{:?}", res);
    }
}
