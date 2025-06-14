use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::sync::Arc;
use std::{path, thread};

use crossbeam::channel::{Receiver, Sender};
use gskits::fastx_reader::fasta_reader::FastaFileReader;
use gskits::fastx_reader::read_fastx;
use gskits::pbar::{self, DEFAULT_INTERVAL};
use hp_tr_finder::{all_seq_hp_tr_finder, UnitAndRepeats};
use mm2::align_processor::hp_tr_metric::HpTrMetric;
use mm2::align_processor::time_err_metric::TimeErrMetric;
use mm2::align_processor::TMetric;
use mm2::params::{AlignParams, IndexParams, InputFilterParams, MapParams, OupParams};
use mm2::{align_single_query_to_targets, mapping2str, NoMemLeakAligner};
use mm2::{build_aligner, query_seq_sender};

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

    #[arg(
        long = "short-aln",
        help = "if the 30 < query_len OR target <200, use short aln mode"
    )]
    pub short_aln: bool,

    #[command(flatten)]
    pub metric_args: MetricArgs,
}
impl Cli {
    pub fn post_process_param(&mut self) {
        if self.short_aln {
            if self.metric_args.index_args.kmer.is_none() {
                self.metric_args.index_args.kmer = Some(5);
            }
            if self.metric_args.index_args.wins.is_none() {
                self.metric_args.index_args.wins = Some(3);
            }

            if self.metric_args.align_args.min_cnt.is_none() {
                self.metric_args.align_args.min_cnt = Some(3);
            }
            if self.metric_args.align_args.min_dp_max.is_none() {
                self.metric_args.align_args.min_dp_max = Some(20);
            }

            if self.metric_args.align_args.min_chain_score.is_none() {
                self.metric_args.align_args.min_chain_score = Some(5);
            }

            if self.metric_args.align_args.min_ksw_len.is_none() {
                self.metric_args.align_args.min_ksw_len = Some(0);
            }
        }
    }
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
        self.io_args.get_oup_filename()
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
    #[arg(long = "out", help = "output filename")]
    pub outfn: Option<String>,

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

    #[arg(long = "qname-suffix", help = "suffix for query name")]
    pub qname_suffix: Option<String>,
}

impl IoArgs {
    pub fn get_oup_filename(&self) -> String {
        let out_fn = if let Some(oup) = &self.outfn {
            oup.to_string()
        } else {
            format!(
                "{}.gsmm2-hp-tr-fact.csv",
                self.query[0].rsplit_once(".").unwrap().0
            )
        };
        if let Some(dir) = path::Path::new(&out_fn).to_path_buf().parent() {
            if !dir.exists() {
                fs::create_dir_all(dir).unwrap();
            }
        }

        out_fn
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

        param.qname_suffix = self.qname_suffix.clone();

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

    #[arg(long = "min-cnt", help = "min_cnt")]
    pub min_cnt: Option<i32>,
    #[arg(long = "min-dp-max", help = "min dp max")]
    pub min_dp_max: Option<i32>,
    #[arg(long = "min-chain-score", help = "min chain score")]
    pub min_chain_score: Option<i32>,
    #[arg(long = "min-ksw-len", help = "min ksw len")]
    pub min_ksw_len: Option<i32>,
}

impl AlignArgs {
    pub fn to_align_params(&self) -> AlignParams {
        let mut param = AlignParams::new();
        param.matching_score = self.matching_score;
        param.mismatch_penalty = self.mismatch_penalty;

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

        param.min_cnt = self.min_cnt;
        param.min_dp_max = self.min_dp_max;
        param.min_chain_score = self.min_chain_score;
        param.min_ksw_len = self.min_ksw_len;

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
            .set_discard_multi_align_reads(self.discard_multi_mapping_reads)
            .set_pass_through_tags(Some(&"dw,ar".to_string()));
        param
    }
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

fn metric_entrance(preset: &str, tot_threads: Option<usize>, args: &MetricArgs) {
    let tot_threads = tot_threads.unwrap_or(num_cpus::get());
    assert!(tot_threads >= 10, "at least 10 threads are needed");

    let target_filename = &args.io_args.target;

    let fa_iter = FastaFileReader::new(target_filename.to_string());
    let targets = read_fastx(fa_iter);
    let targetname2seq = targets
        .iter()
        .map(|v| (Arc::new(v.name.clone()), v.seq.clone()))
        .collect::<HashMap<_, _>>();

    let all_regs = vec![
        UnitAndRepeats::new(1, 3).build_finder_regrex(),
        UnitAndRepeats::new(2, 3).build_finder_regrex(),
        UnitAndRepeats::new(3, 3).build_finder_regrex(),
        UnitAndRepeats::new(4, 3).build_finder_regrex(),
    ];

    let region2motif: HashMap<
        Arc<String>,
        Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>,
    > = all_seq_hp_tr_finder(&all_regs, &targetname2seq)
        .into_iter()
        .map(|(seqname, region2motif)| (seqname, Arc::new(region2motif.flatten())))
        .collect::<HashMap<_, _>>();

    let index_params = args.index_args.to_index_params();
    let map_params = args.map_args.to_map_params();
    let align_params = args.align_args.to_align_params();
    let oup_params = args.oup_args.to_oup_params();
    let inp_filter_params = args.io_args.to_input_filter_params();

    let mut aligners = build_aligner(
        preset,
        &index_params,
        &map_params,
        &align_params,
        &oup_params,
        &targets,
        tot_threads,
    );
    aligners.iter_mut().for_each(|aligner| {
        aligner.mapopt.best_n = 10000;
        aligner.mapopt.pri_ratio = 0.2;
        aligner.mapopt.max_gap = 100;
    });

    assert!(
        args.io_args
            .query
            .iter()
            .filter(|name| !name.ends_with(".bam"))
            .count()
            == 0,
        "bam file needed for this tool. but got {:?}",
        args.io_args.query
    );

    thread::scope(|s| {
        let aligners = &aligners;
        let targetname2seq = &targetname2seq;
        let inp_filter_params = &inp_filter_params;
        let oup_params = &oup_params;
        let region2motif = &region2motif;
        let (qs_sender, qs_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            query_seq_sender(
                &args.io_args.query,
                qs_sender,
                inp_filter_params,
                oup_params,
            );
        });

        let align_threads = tot_threads - 4;
        let (metric_sender, metric_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..align_threads {
            let qs_recv_ = qs_recv.clone();
            let metric_sender_ = metric_sender.clone();
            s.spawn(move || {
                hp_tr_metric_worker(
                    qs_recv_,
                    metric_sender_,
                    aligners,
                    targetname2seq,
                    region2motif,
                    oup_params,
                )
            });
        }
        drop(qs_recv);
        drop(metric_sender);

        let oup_filename = args.get_oup_file();
        dump_metric_worker(metric_recv, &oup_filename, true);
    });
}

pub fn hp_tr_metric_worker(
    query_record_recv: Receiver<gskits::ds::ReadInfo>,
    metric_sender: Sender<Box<dyn TMetric>>,
    aligners: &Vec<NoMemLeakAligner>,
    targetname2seq: &HashMap<Arc<String>, String>,
    region2motifs: &HashMap<Arc<String>, Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>>,
    oup_params: &OupParams,
) {
    for query_record in query_record_recv {
        let metric = build_hp_tr_metric(
            &query_record,
            aligners,
            targetname2seq,
            region2motifs,
            oup_params,
            false,
        );
        metric_sender.send(metric).unwrap();
    }
}

pub fn dump_metric_worker(metric_recv: Receiver<Box<dyn TMetric>>, fname: &str, enable_pb: bool) {
    let mut writer =
        BufWriter::new(File::create(fname).expect(&format!("create file error: {}", fname)));

    writeln!(
        &mut writer,
        "{}",
        mm2::align_processor::hp_tr_metric::METRIC_CSV_HEADER.join("\t")
    )
    .unwrap();

    let pb = if enable_pb {
        Some(pbar::get_spin_pb(
            format!("gsmm2-hp-tr-metric: writing metric to {}", fname),
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

pub fn build_hp_tr_metric(
    query_record: &gskits::ds::ReadInfo,
    aligners: &Vec<NoMemLeakAligner>,
    targetname2seq: &HashMap<Arc<String>, String>,
    region2motifs: &HashMap<Arc<String>, Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>>,
    oup_params: &OupParams,
    debug: bool,
) -> Box<dyn TMetric> {
    let hits = align_single_query_to_targets(&query_record, aligners);

    if hits.is_empty() || (hits.len() > 0 && oup_params.discard_multi_align_reads) {
        return Box::new(TimeErrMetric::new(
            query_record.name.clone(),
            "".to_string(),
            vec![],
            Some(600),
        ));
    }

    let mut hits = hits
        .into_iter()
        .filter(|v| v.is_primary || (v.is_supplementary && !oup_params.discard_supplementary))
        .collect::<Vec<_>>();

    if debug {
        hits.sort_by_key(|hit| hit.query_start);
        hits.iter().for_each(|hit| println!("{}", mapping2str(hit)));
    }

    let target_name = hits[0].target_name.as_ref().unwrap().clone();

    hits = hits
        .into_iter()
        .filter(|hit| hit.target_name.as_ref().unwrap() == &target_name)
        .collect::<Vec<_>>();

    let target_seq = targetname2seq.get(&target_name).unwrap();

    let target_seq_region2motif = region2motifs.get(&target_name).unwrap();

    let mut metric = HpTrMetric::new(
        query_record.name.clone(),
        if hits.is_empty() {
            Arc::new("".to_string())
        } else {
            target_name.clone()
        },
        target_seq_region2motif.clone(),
    );

    hits.into_iter()
        .for_each(|hit| metric.add_align_info(hit.into()));

    metric.compute_metric(target_seq, &query_record.seq);

    Box::new(metric)
}

#[cfg(test)]
mod test {
    use std::{collections::HashMap, sync::Arc};

    use gskits::{
        ds::ReadInfo,
        fastx_reader::{fasta_reader::FastaFileReader, read_fastx},
    };
    use hp_tr_finder::{all_seq_hp_tr_finder, intervaltree, UnitAndRepeats};
    use mm2::{
        build_aligner,
        params::{AlignParams, IndexParams, MapParams, OupParams},
    };

    use crate::build_hp_tr_metric;

    #[test]
    fn test_compute_metric() {
        let ref_file = "test_data/MG1655.fa";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);

        let all_regs = vec![
            UnitAndRepeats::new(1, 3).build_finder_regrex(),
            UnitAndRepeats::new(2, 3).build_finder_regrex(),
            UnitAndRepeats::new(3, 3).build_finder_regrex(),
            UnitAndRepeats::new(4, 3).build_finder_regrex(),
        ];

        let targetname2seq = targets
            .iter()
            .map(|v| (Arc::new(v.name.clone()), v.seq.clone()))
            .collect::<HashMap<_, _>>();

        let region2motif: HashMap<
            Arc<String>,
            Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>,
        > = all_seq_hp_tr_finder(&all_regs, &targetname2seq)
            .into_iter()
            .map(|(seqname, region2motif)| (seqname, Arc::new(region2motif.flatten())))
            .collect::<HashMap<_, _>>();

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
            10,
        );

        //fwd
        let seq = b"GAGAGAAGAGATGTTCTACTGGTGTGCCGATGGTGGCGGATTAGTCCTCGCTCAATCTGGGGTCAGGCGGTGATGGTCTATTGCTATCAATTATACAACAATTAAACAACAACCGGCGAAAGTGATGCAACGGCAGACCAACATCACTGCAAGCTTTACGCGAACGAGCCAGACATTGCTGACGACTCTGGCAGTGCATAGATGACAATAAGTTCGACTGGTTACAAAACCCTTGGGGCTTTTAGAGCAACGAGACACGGCAATGTTGCACCGTTTGCTGCATGAATTGAATTGAACAAAAATATCACCAATAAAAAACGCCTTATGTAAATTTTTCAGCTTTTCATTCTGACTGCAACGGGGCATATGTCTCTGTGTAGATTAAAAAAAGAGTGTCTGATGCAGTTCTGAACTGGTTACCTGCCGTGAGTATAATTAATAAAATTTTATTGACTTGGTCACTAAATACTTTAACCAATATAGGCATAGCGCGACAGACAGATAAAATTACAGAGTACACAACCAATCCCATGGAAACGCATTAGCCCACCATTACCCCACCATCACCATTACCAACAGGTAAAAACGGTGCGGGCTGATCGCGCTATCAGGAAACACAGAAAAAGCCCGCACCTGACATGCGGGCTTTTTTTTTTCGACCAAAGGTAATATACGAGGTAACAACCATGCGAGTGTTTGAAGTTCGGCGGGTTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTTGTGCCGATATTCTGGAAAGCAATGCAGCAGGGGCAGGTGGCCACCGTCTCCTCTGACCCCGCCAAAATCACCACACCTGGTGCGATGATTGAAAAACCATTAGCGGCAAGGAGCTTTAAACCCAATATCAAGCGATGCCGAACGTATTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCAGCCGGGTTTCACCGCTGGCCAATTGAAAACTTTCGTTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACCTGCGCTGATTTGCCGTGCGAAGAAATGTCGATCAGCCATTATATCTCTC";

        // rev
        let seq = b"TATTCGATTGGAAATGCCCTTGACCAGGTAATTCAGTCTCATCACGAACCATGAGCGTACCTGGTGCTTGAGGATTTCCGGATATATTTTAATCAGGCAAGGGACGATCTGAACTGGCGGATGGGGGTAAATGGTGGCGGGGTTGAAGAACTTTAGCGCGAAGTAGGAAAGCTCCTCGCTTCGTTGCTCTAAAAAAGCCCCCAGGCGTTGTTTGTACCAGTCGACAAGTTTTAAATGTCATCTGCCACGCCAGAGTCGTCAGCAATGTCATGGCTCGTTCGCGTAAAAGCTTTGACAGTTGATGTTGGTCTGCCCGTTGCATCACTTTTCGGCCGGTTGTTGTATTAAGTTGCTAAATTGATAAGCAATAGATCACCGCCTGCCCAGATTGAGCGAAGGTAATCCGCCACCATCGGCAACACAGTAATGACGTCAGCCAGCCAAACGCTAACTCTTCGTTTAGTCAACCCGGAATCCTTCGCGACCCACCCAGCGCGGCATGGCCATCCATGAAAGATTTTTCCTCTAACAGCGCACCAGTTCAACTGGCGTGGCGTAAGTCAATGATATTTCGCCCGACTGCGCGCAGTGGTGGCGACAAGTGAAATCGACATCGTGTAACGATTCA";
        let query_record =
            ReadInfo::new_fa_record("tt".to_string(), String::from_utf8(seq.to_vec()).unwrap());
        build_hp_tr_metric(
            &query_record,
            &aligners,
            &targetname2seq,
            &region2motif,
            &OupParams::default(),
            false,
        );
    }

    #[test]
    fn test_compute_metric2() {
        let ref_file = "test_data/MG1655.fa";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);

        let all_regs = vec![
            UnitAndRepeats::new(1, 3).build_finder_regrex(),
            UnitAndRepeats::new(2, 3).build_finder_regrex(),
            UnitAndRepeats::new(3, 3).build_finder_regrex(),
            UnitAndRepeats::new(4, 3).build_finder_regrex(),
        ];

        let targetname2seq = targets
            .iter()
            .map(|v| (Arc::new(v.name.clone()), v.seq.clone()))
            .collect::<HashMap<_, _>>();

        let region2motif: HashMap<
            Arc<String>,
            Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>,
        > = all_seq_hp_tr_finder(&all_regs, &targetname2seq)
            .into_iter()
            .map(|(seqname, region2motif)| (seqname, Arc::new(region2motif.flatten())))
            .collect::<HashMap<_, _>>();
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
            10,
        );

        // // aligner.mapopt.min_cnt = 2;
        // aligner.mapopt.best_n = 100000;
        // // aligner.mapopt.min_dp_max = 10; // min dp score
        // aligner.mapopt.pri_ratio = 0.5; // min dp score
        // aligner.idxopt;
        println!("{:?}", aligners[0].mapopt);
        println!("{:?}", aligners[0].idxopt);

        aligners[0].mapopt.best_n = 100000;
        aligners[0].mapopt.pri_ratio = 0.1;

        let seq = b"AAAAGAGAGAACGAAAAAAAGAGAGAGATGAATCCAGCGAAAACCATGTACGGCTTAAATTCAGCCACCTTATTATCGCCAGAATCCCTTTTCGCACGAACTGTCCTTTGATGCGAAACAAAATTTGGCGCATTCACTGTGCATACCGGCCGGTCAAGCAGCAGGTAACAATGCGGAGCGGCGCGCCTGCCGTCCATGATCGGGATTTCCGGCGCCGTTAGACTTCGAAAACAATGTTATCGATGCCCAAAGCCCCGGAGAGCAGCATTGAGTGCCATCGGTTGAAATCCGTACATCTCTCGTGACCAGACACGTACAGAGCATGGTATCACGCACAGAATTTGGCATTCGGCCGGGGAAATCTACCGTGATTCAAAGTCGGTGCGACGATAGATGACCCGGTGTTGCCGGCGCAGGGCGTAACGTCAGGGTGACTTTTTGCCGGTATGTAAACCACAAACCCAGTCGCCTGAACAGAACGTTTATGGTCCTTGTATTGATCATCGTATTATCTCAGCCAATTACCTATCCAACCGAAGTGTACTATACATTCGGCCGCAGTTTTAGCACAAAAGAGCCTCGAAACCCAAATTCCAGCAATTCTTATTCAGCTTGCTTACGAGGAATCTGGGAATCCAGATAATCCGGCTCTTCCATTTGCGCGCATTGTCCAATTCACGACTTTAGCTTCCGCTTCTGCTCCCTGGGTCAGCGGAGCCCATCCCATGCTGCTGTAGCGATCCATCACTGGCTGCTGAACCTGCTTATTGGTCACCAGATCTTCAGAACGAACCATGCTGAATAGCCCTGCCATTGCCGCATCCAGTCCGCTCAAGGATGTTCCCATAATACGCCCCCCCCCCCCCCAGTTACGCATCGACCTTGCTTCCATCAGGCTGACCAACGGGGCTGGAAATACGTTTCGGGTTGCTTCTTTCATCAGGCCAGACGTGACCCGTGCATCACCGCAAATCAAACGGTCTGCTGGGTAAATCCTGTTTGTAAGGAGCATGACGCCAGCATTTGTAATTTGCATGATCGGTAACCTGGGCATGATTCATAAACACCACTGCAAATTTTTGTGTCGTGCCTGGTCTACTAGTCGTAAAAATTGATCGCGGACAAATTTCGCCAGCATCTCTCTCAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGTGTCTGGGCGAATATTTCCGCGTTCAATTTTTTACGACATAGTAGCTCCCAGGCACGACAGCAAAAATTTGCAGTGGCTGTTTATGAATCTGCAGGTTTTACCGTCTGCAAAATTACAATGCTGGCGTCATGCTCGTACAAATCAGGATTTACCAGCGAGACGTTTTGATTTAGCGGTGATCGCCACGGGTCGACGTCTGCTGATGAAGAAAGAAAGCAACCCGAACGTATTTTTCTCCAGCCCGTGTCAGGCCTGATGGAAGCAAAGGTCGATGCGTGTAACGTGGGTATTATGGGAACATCCTTGAGCGGACTGGATGCGGCAATGGCAGGCTATTCACATGGTTCGTTCATTGAAGATCTGGTGACCAATAAGCAGGTTCAGCAGCCATGGATGGATCGCTACATGCAGCATGGGATGGCTCCGCTGAAGGAGCAAAGCCGGTGCTAAAGTCGTGAATGACAATGCGCCGCAAACTGCGAAGAGCCGGATTCCTGGATATCCCTGCATTCCTGCGTAAGCCAACTGATTAAAGAATTGACTGGAATTTGGGTTCGAGGCTCTTTGTGCTAACTGGCCCCCGAATGTATAGTACACTTCGGTTGGATAGTAATTGGCGCAGATATTTCATGATCAAACAAAGACACTTAAACGTATCGTTCAGGCGTCGGGTGTCGGTTTACATACCGGCAAGAAAGTCACCCGCGATTACGCCTGCGCCGGCCAACACCGGGGTCATCTAACGTTCGCACCGACTTGAATCCACCGGTTAGATTTCCCGGCCCGATCCAAAATCGTGCGTGATACCATGCGCTGTACGTGTCTGGTCAACGGCAGAGTACGGATTTCAACCCGTAGAGCATACCTCAATGCTGCTCTCGCGGGCTTGGGCAATCGATACATTGTTATCGAATTAACGCGCGGAAATCCCGATCATGACGGCAGCGCCGCCCGTTTGTATTACTCTGCTTGACGCCGGTATCGACGAGTTGAACTGCGCAAAAAATTGTCCGCATCAAAGAGACTGTTCTGTCGTTGATGGCGTAAGTGGGCTGATTTAAGCCGTACAATGGTTTTTCGCTGGATTCATCTCCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGATTAAAGAGAGAGATGAAATCCAGCGAAAACCATGTAGGCTTAACATTCAGCCACTATCGCCTCTTCGACAGAACAGTCTCTTTGATGCGAACAATTTTTTGGCGCAGTTCAACTCCGTCGATACCGCGTCCAAAGCAGCAGGTATACAAATACAGAGGCGGCGCTGCGTCATGATCGGGATTTTCCGGCGCGTTAACTTTCGATAACAAGTTATCGTGCCCAAGCCCCGGAGAGCAGCATTGAGGTGCTCTACGGTTGATAATCCTACATCATGCTCCGTTGACCAGAAACACGTACAGACATGGTATCACGCACAAGATTTGGCATCGGCGGGAAATCTACCTCGTGGATTCAACGTCGGTGCGACGTAGGTGTCCCCGGTGTTGGCCGCGCAGGGTCACGGCAGGTGACTTTCTTGCCGGTATGTAAACCGACACCCGTCCGCCTGAACGATACGTTTAAGTGTCCTTTGTTTGATCATCGTATTATCTCGCCAATTACCTATCCAACCGAAGTGTCTAATACATTCGGCGGGCATTTAGCACAAAAGAGCCGAAACCCAAATTCCCATCAATTCTTAAATCAGCTTGCTGACGCAGAATGCTGGATATCCAGAATCCGGGCCTTTCGCAGTTTGCGCCGCATTGTCATCAAGACTTAGCACGGCTTCGCTCCTGGGTCAGCGGAAGCCATCCCCATGCTGCTGTAGCGATCCTCACTGGCTGCTGAACCTGCTTATTGTCACCAAGATTCCTTCTTTGAACCGAACCATGCTGAATAGCCATCCATTGCCGCATCCAGTCCGCTCAGGAGTTCCCAAAATACCCACGTTACACGCATCGACCTTTTGCTCCACAGGCCTGACCACGGGCTGGGAAAATACGTTCGGTTGCTTCTTCTTCACAGGCCAGACGTGACCCCGTGGCGATCACCGCTAAATCAAACGTCTCGCTGGTAATCCTGATTTGCAGCGAGCATGACGCCAGCATTTGTAATTTGCAGATCGTAACTGGCATGATTCATAAACAGCCACTGCAAATTTGCTGTCGTGCCTGGTCTACAGTCGTAAAAATTGATCGCGGAAATGATTTCGCCCAGCATCCTCTCCTCAAACAACAACGGATGGAGGAGGAAAATAACAGGAAAAAAGAGAGAGATGCTGGGCGAATATTTGGATCAATTTTTTACGACCTAGTAGACCAGGCACGACAGCAAAATTGCAGTGCTGTTTATGAATCATGCCAGTTACGGACTGCAAATTACAAAGCTGGCGTCATGCTCGCTACAAATCAGGTTTACCCAGCGAGGATTTTGATTTATAGCGTGATCGCCACGGGTCACGTCTGGCCTGATGAAGAAGAAGCACACCGAACGTATTTTTCCCAGCGTGGTCAGCTGATGGAAGCAAAGGTCGAATGCAGTAGTAACGTGGGTATTATGGGAACACCTAGCGGACTGGATGCGGCAAGCAGTGGCTATTCAGCATGGTTTCGTTCATTGAAGATTGGTGACCAATAGCAGGTTCAGCAGCCAGTGAGGAGCGCACCAGCAGCATGGATGGCTCCGCTGACCCAGGAGCAGAAGCCCGGTTGCTAAATCGTGAATGACAATGCGCCGCAACTGCGGAAAGAGCCGGATTAGATCTGGATATCCCAGCATCCTGCGTAAGCAAGGCTGATTAAGAATTGACTTGGAATTGGGTTTCGAGGCTCTTTGTGCTAAACTGGGCCCGACCCGCCGAATGTATAGTACACTTCGTGGATAGGAATTTTGGGCGAGATAGATAGTGACCTATACGATGATTTCAAACAAAGGACACTTAAAACGTATCGTTCAGGCGTCGGGTGTCGGTTTACAATACCGGCAAGAAAGTCACCTGACGTTACGCCCTGCGCGGCCACACCGGGTCATCTATCGTCCCACCGACTGAATCCACCTAGATTGTCCCGGCCGATGCCAAATCTGTGCGTGATACCTGTTAAGCTCTTGTAACGTGTCTGGTCAACGAGCATGAGTACAGGATTTCAACCGTAAGAGCTCCTCAATGCTGCTCTCGCGGGCTTGGGCATCATAAATTTGTTATCGAATTAACGACGCCGGAATCCCGATTCATGGACGGCAGCGCCCGCTCCGTTTGTATACCTGCTGCTGGACGCCGGTATCGACGGAGTAGAAGCTGCGCCAAAAAATTTGTTCGCATCAAAGAGACTGTTCGTGGGTCCGAAGATGGCGATAATGGGCTGAATTTAAGCCGTACAATGGTTTTTTCGCGGATTCATCTCTCTCCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGAAATCCAGCGAAAAACCATTGTACGGCTTAATTCAGCCACTTTATCGCCATCTCAGACACGAACAGTCTCTTTTGAGAACAAATTTTTTGGCGCAGTTTCAACTCGTCCGATACCGGCGTCAAGCAAGCTGGTAATACATAGGTGCGGCGCTGCCGTCAAGATCGGGATTTCCGGCGTCGTTAACTGCGATAACAATGTTACGATGTCCCAAGCCGCGGAGAGCAGCATTGAGGTGCTCTAAACGGTTGAAAATCCGTAACATCATGCTCGTTTGACCAGACACGTAAACAGAGCATGGTATGCCACGACAGATTTGGCATCGGCGGGGGGGGGGGGGGGGGGGGGGGGGGGAAATCTACCGTGGATTCAAGTCGGTGCGACGAAGATGACCCCCGGTGTTGGCCGGCGCAGGGCGTAAACGTCAGGTGAACTTTCTGCCGGTATGTAAACCGACCACCCGTCGCTGAACGATACGTTTAGTGTCCTTGTTGATCATCGTAATTATCTCGCCAATTACCTATTCAACCGAAGTGTACTAGTACATTCGGCTGGCCAAGTTTAGCAAAAGAGCCTCCGAAACCCAAATCCAGTCAATTCTTAATCACTTGCTTACGCAGGTAGCTGGGATATCCAGATAATCCGCTCTTTTCGCAGTTTGCGCGCATTGTCATTCCACGACTTTAGCACCGGCTTCTGCTCCGGGTCAGCGGAGCCATCCCATGGCTGCTGGTAGCGATCCATCACTGCTGCTGAACCTGCTTATTGGTCACCAGATCTTCAATGAACGAACCATGCTGAAATAGCCACTGCCATTGCCGCATCCAGTCCGCTCAAGGATGTTCCCATATATCCCACGACACGCATCGACCTTTGCTTCCTCAGGCCTGACCACGGCTGGGAAAAATACGTTCGGGTTGCTTCCTTCTTCATCAGGCCAGACAGTGACCCGTGGGTTCCCCACCGCTAAATCAAGTCTCCTGGGTAATCTGATTTGTAAGCGAGCATGACGCCAGCATTTGTATTTTGCAGATCGGTAACCTGGCATGATTCAAAAGCCCTTGCAAATTTTTGCTGTCGTGCCTGGTCTACTAGTCGTAAAAAATTGAATCCGGAAATATTCGCCACATCTCTCTCACAACAACAACGAGGAGGAGGAAAAAAGGAGAGAGAGCTGGGCGAATAATTTCCGCGATCAATTTTTACGACTAGTAGTCCAGGCACGACGCAAAAATTTTGCAGTGCGGTTTATGAATCATGCCAGGTACGATCTGCAAATAAATGCGGCGTCATGCTCGCTACAAATGCGAGGATTTTACCCAGCGAACGTTTGATTTAGCGGGTTCGCCCTCGGGTCACGTGCTGCCTGATGAAGAAGAAGCAACCCGAACGGTATTTTCCCAGCCCTGGTCAAAGGCTGATGGAAGCAAAGTCGATGCGTGTAAACGTGGTAAAATTATGGGAACATCCTTGAGCGGACTGAGATGCGGCAATGGCAGTGGGCTATTCAGCATGGTTTCGTTCATTGAAGATCTGGTGACCCATAAGCAGGTCAGCAGCCAGTGATGGATCGTCTACCAGCAAGCATGGGAATGGCTCCGCTGACCCAAGGAGCAGAGCCGTTGCTAAAGCGTGAATGACATTGCGCCGCAAAACTGCGAAAGAGCCGGATTATCTGGATATCCCAGGCATTCCTGCGTAAGCAAGCTGATTAAATTGACGGAATTTGGGTTTCGAGGTCTTTGTGCTAAACTGGCCCCGCGAATGATAGTACACTTCGGTGGATAGGTAAATTTGGCGAGATAATACGATGATCAAACAAAGGACATTAAACGTATCGTTCAGGCGCGGGTGTCGGTTTAACATACCGGCAAGAAGTCACCTGACGTTAGCCTGCGCCGGCCAACAACCGGGGTCATCATCGTCGCACCACTGAATCACACCGTGTAGATTTCCCGGCCCGATGCCAATCTGTGCAGTGATACCTGCTCTGTACGTGTCTGGATCAACGGCATGATGTAACGGATATGCAACGTAGATAGAGCAGCCTCAATGCTGCTCTTCGCGGGCTTGGGCATCATAACATTGTTATCGTGGTTAACGCGCCCCGGAATCCCGATCATGGACATAGGCAGCGCCGCTCCGTTTGTATACCTGCTTGCTTGACGCCGGTATCGCGAGTTTGAACTGCGCCAAAAAATTTGTTCGCATCAAGAGACTTGTTCGTGTCGAAGATGCGATAAAGTGGGCTGAATTTAAGCCGTACAATGGTTTTCTCGCTGGATTTCATCTCTCTCAAACAACAACGGAGGAGAGGAAAAAAGAGAGGATGAATCCCAGCGAAAAACCATTTGTACGGCTTAAATCAGCCCACTTATCGCCATCTTCGACACGAACAGTCTCTTTGATGCGAACAATTTTTTGCGCAAGTTCAACTCGTCATACCGCGTCAAGCATGCAGGTATACAAACGGAGCGCGCTGCCGTCCATGATCGGGATTTCCGGCGCGTTAACTTCGATAACAATGTTATCGATGCCCAAGCCCGCGGAGAGCAGCATTGAAGGTGCTCTACGGTGAAATCCGTACATCAGCTCGTTGACCAACACGTACAGAAGCATGGTATCACGCACAGATTTGGCATCGGCCGGAATCTACCGGTGATTCAAGTCGGTGCGAAGATGAATGACCCCCGGTGTTGGCCGGCGGCAGGGCGTAACGTCAGGGTGACTTCTTGGCGGTATGTAAACCCACAGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGTCGCCGAACGATACGTTTAAGGTCTTTTGTTTGATCATCGATTATCTCCCAAATTACCTATCCAACCCGAAGTGTACCTATACATTCGGCGGCCGTTTAGCACAAAGAGCCCGAAACCCAATTCCAGTCCATTCTTATCAGCTTGCTTACCCAGGAATCTGGATATCCAGAAATCCGGCTCTTTCCCGAGTTTGGCGCATTGTCATTCACGACTTTAGCAACCGGCTTCTGCTCCTGGGTCCAGCGAGGCCCTTCCATGCTGCTGTAGCGATCCATCACCTGGCTGCTGAACCTGCTTATGGTCCCACCAGATCTTCAATGAACGAACCATGCTGGAATACCACTGCCTTGCCGCATCAGTCCGCTCAAAGGATGTTCCCATAATTCCCACAGTTACAGCATCGACTTGCTCCATCAAGCCTGACCACGGCTGGGAAATACGTTCAGGTTGCTTTTCTTCATGCAAGGCCCAGAACCGTGACCGTGGCGATCACCGCTAAAATCAAACGTCTCGCTGGGTAAATCCTGATTTGTAGCGAGCATGACGCCAGCATTTGTAATTTGCAGATCGGTAACTGGATGTTCATAAACACACTGCAAATTTTTGCTGTCGTGCCAGGTCTACTAGTCGTAAAAATTGATCGCGGAATATTCGCCCAGCATCTCTCTCAAACAACAACAGAAGGAGGAGGAAAAAAGAGAGAGATGCTGGGCGAATAATTTTGCCCGATCAATTTTTTACGACTAGTAAGAGCCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCAGCACGCAGCAAAAATTTGCAGTGGCTAGTTTATTGAATCATGCCAGGTTTCCGAACTGCAAATTACAAAATGCTGGGTCAGCGCTACAAAATCAGGTTTTAACCCAGCGAGACGTTTTGATTTAGCGGTGATCGCCACGGTCTCGTCTGGCCCTGATGAAAAGAAGCAACCCGAACGTATTTTCCAGCCCGTGGTCAGGCCTGATGGAAGCAACAGGTCGATGGCGGTACACGTGGGTATTATGGGAACATCCTTGAGCGGACTGTGCGGCAATGCAGTGCTATTCAGCATGGTTCGTCCATTGAAGATCTGGGACCAATAGCTGGTTCAGCATGCCAGTGATGGATCGCTACCAGCAGATGGATGGCTCCCGTACCCAGGAGCAGAAGCCGGTTGCAAAGTCGTGGAATGACAATGCGCCGCAAACTGCGAAGAGCCGGAAATCTGGATAATCCCCCAGCATTCCTGCATAAGCAAGCTGATTAAGAATTGACTGGAATTTGGGTTCAGGCTCTTGTGCTAAACTGTAGCCCCGCCGATGATAGTCACTTCCGTTGGCTAGGTAAATTTGGCAGATAATCGATGATCAAACAAAGGACACTAAACGTATCGTTCAGTGGGCGACGTGTCGGTTTACATACCGCAGAAAAGTCCACCCTGACGTTACGCCCTGGCGCCGCAACACCGGGTCATCTAACGTCGCACCGACTTGAAATCCACCGGAGATTCCGGGCCAGATGCCAAATCTGTGCGTGGATACACAAGCTCTGTACGTGGTCTGGTCAACGACATGATGTACGGATTTCAACCGTAGAGCACCTCAATGCTGCTCTCGCGGGCTTGGGCATCATAAACATTGTTATCGAAGTTAACGCGGCCGGAAATCCGATCATGGCGGCAGCGCCGCTCCGTTTGATACCTGTGCTGCTTGACGCCGTATCGACGAGGTTGAACTGCGCCAAAAAATTTGTTCCGCATCAAAGAGAACTGTTCAGTGTCGAAGATGGCGAAAGTGGGCTGAATTTAAGCCCGACAATGGTTTTCGCTGGATTTCATCTCTCTCAAACACAACGGAGGAGGAGGAAAAAAGAGAGAGATGAAATCCAGCGAAAAACATTGTACGGCTTTAAATTCAGCCCACCTTATGCATCTTCGACACGAACAGTCCTCTTGATGCGCAAATTTTTGGCGCAGTTCACTCGTCATACCGGCGTCAAGCAGTGTATCAAATATCGGAGCAGGCCTGCCGCCATGATCGGGATTTCCGGCGCGTTAACTTGATAACAATGAATCCGAATGCCCAAGCCCGCGAGAGCGCATTGTGGTGCTCTACGTGTTTGAAATCCCGTAACATCATGCTCGTTGACCGACAGACAGAGCCATGGTATCACGCACAGATTTGGCATCGGCCGGGAAATCTACAGGTGGATCAAGTCGGGCGACCGATAGATGACCAGTGTTTGGCCGCGCCAGGCGTAAACGTCAGGGTGACTTTCGCCGGTATGTAATCCGACACCCGTCCCCGCCTGAACGATACGTTTAAGTGTCCCTTTTGTTTTGATCATCGTATTACTCGCCAAATTACCTATCCAACCGAAATGGTACTATACATTGGCGGCCAGTTTTAGCACAAAGAGCCTCGAAACAAATTCCCGTTATTCCTTAATCAGCTTGCTACGCAGGATGCTGGGATATCCAGAAATCGGCTCTTTCCAGTTTGCGGCGCATTGTCCATTCACGGACTTTAGCAACCGCTTCTGCTCCTGGGTCAGCGGTGCCATCCCTGCCTGCTGTAGCGATCCATCACTGGCTGCTGAACCTGCTTTATTTGGTCACCAGATCTTAATGACGAACCATGCTGAATAGCCACCTGCCTTGCCGCTCCAGTCCGCTCAAGGTGTTCCCATAATACACTTACACGCATCGACCTTTGCTTCACAGGCCTGACCACGGGCTGGGACAAATACGTTCGGGTTGCTTCGCTCTTCATCAGGCCAGACGTGACCGGCGATCAACCGCTAAATCAAACGTCTCGCTGGGTAAATCGCTGATTGTAGCGAGCAGACGCCAGCATTTGTAATTTGCAGATCGGTAAACCTGGCATGATCATAAACACCACTGCAAATTTTTGCTGTCGTGCCTGTCTACTAAGTCGAAATTGATCGCGAATATTCGCCCAGCATTCTCTCCAACAACAACGGAGGAGAGGCAAAAAGGAGAGAGATGCTGGGCAATATTTCCCGCGATCAATTTTTACGACTAGTAGACAGGGGCACCGACAGCAAAAAAATTTGCAGTGGCTGTTTATGAATCATGCCAGTTACGATCTGCAAATTACAAATGCTGGCGTCATGCTCGCCTACAAATCAGGATTTACCCAGCGAACGTTTGATTTAGCGGTGATGCGCCCGGTGTCACCGTCTGGCCTGAGAGAAGAAGCAACCCGAACGTATTTTTCCCAGCCCGTGGCAGGCCTGTGGAAGCAAAGGTCGATGCGTGAACCAGTGGGTATTATGGTGACATCCTTTGAGCGGACTGATGCCGCAATGGCAGTGGCTATTCACATGGTACGTTCATTGAAGAACTGTGACCAATAAAAAGCCAGGTTCAGCAGCCAGTGATGGATCGCTACCAGCAGAGGGATGGCTCCGCTGACCCAGGAGCAGAAGCCGGTTGCTAAAGTCTGATGACAAGGCCGCAAAACTGCGAAAGGAGCCGGATATCTGGATATCCCACATCCTGGCGTAGCAAGGCTGATTAAGAATTGCTGAATTTGGTTCGAGGCTCTTTGTGCTAAACTGGCCGCCGAAATGTATAGTACATTGTTGATAGGTAAATTTTTGGCGAGATAATACGATGATCAAACAAAGGACACTTAAACAGTACGTTCAGGCGACGGTGTCGGTTTACATACCGCAAGAAAGTCACCCTGACGTTACGCCCTGCGCCAGCCAACACGGGGTCATCTATCTCCACGACCTTGATCCGACCGGTAGATTTCCGGCGATGCCAAATCTGTGCGTGATACATGCTCTGTACGTGTCTGGTAACGAGCATGATGTAGGTTTTCAACCGTAGATCACCTCAATGCTGCTCTCGCCGGGCTTGGGCATCGATAACATTGGTTATCAAAGTTAACCGCCGGGAATCCGTCTGGACGGCAGCGCGCTCGTTTGTATACCTGCTGCTTGACCCGTACGAGAGTTGAACTGCCCAAAAAATTTGTTCGCATCAAAGAGACTGTTCGTGTCGAAGATGGCGATAAGTGGGCTGAATTTAAAGCCGTACACATGGGTTTTTCGCTGGATTTCATCTCTCTCAACAAAACGGAGGAGAGGAAAAAAGAGAAGATGAAAATCCAGCGAAAAAACCATTGTCGGCTTAAATCAGCCCACTTATCGCCATCTTCGACACGAACAGTCTCTTTGATGCGAACAAATTTTTTGGCGCAGTTCAACTCGTCGAATACCGGCGTCCCAAGCACAGGTATACAAACGGAGCGGCGCTGCCGTCCTGGATCGGGATTTCCGCGCGTTAACTTCGTAACAATGTTATTCGAAAGCCCAAAGCCCGCGAGAGCAGCATTGAGGTGCTCTACGGTTTGAAACCCGTACATCCATGCCTCGTTGACCAGACACGTTACAGGAGCATGTATCACGCACAGATTTGGCTCGGCCGGGAAATCTACCGGTGATTCAAGTCGGTGTCGACGATAGATGACCCCGGTGTTGGCCAGGCGCCAGGCGTAACGTCAGTGACTTTCTTGCCGGTATGTAAACCCGACCCCACCCTCGCCTGAACGATAGCGTTTAAGTGTCCTTTGTTGATCATCGTATTATCTCGCCAATTACCTATCCAAACCGAAGTGTACTATCACATTCGGGGCCAAGTTTAGCACAAGAGCTCGAAACCCAAAGTTCCAGTCATTCTAATCAGCTTGCTTACGCAGGAATGGCTGGGATATCCAGATAATCCAGCTCTTTCGCAGTTTGCGCGCATTGTCATTCACGACTTTAGCAACCGGCTTCTGCTCCTGGGTCAGCGAGCCATCCCATGCTGCTGGAGGATCCATCACTGCCTGCTGAACCTGCTTATTGGTCACCAGTTCTTCAATGAACGAACCATGCTGCAATAGCACTGCCATTGCCGCATCCAGTCGCTCAAGATGTTCCATAATACCCACGTTACGACGCATGACCTTGCTTCCATCAGGCCTGACCACGGCTGGGAAAATACGTTCGGGTTGCTTCTTCTTCATCAGGCCAGACAGTGGACCCGTGGCGAGCACCGCTAACAAACTCTCGCTGGGAAATCCTGATTTGTAGCGAGCATGACGCAGCATTTTAATTTGCAGAGTCGGTAACCTGGCATGATGCAAAACAGCCACTGCAAATTTTTGCTGCGTGGCCTGGTCACTAGTCGTAAAAATTGATCGCCGGAATAAATATTCCGCCCAGCATCCTCTCAACAACAACGGAGGAGGAGGAAAAAAGAAAGAGATGCTGGGCGAATATTCCGCGATCATTTTACGACTAGTTGTCCTGGCACGACAGCAAAAAATTTGCAGTGCTGTTTATGTTGCTTGCCAGGTTACCAACCTCTGCAATTTACAAAGCTGGCGTCAGCTCGCTACAAAATCAGGATTTACCCAGCGAGACGTTTGATTACGCGGGATCGCCACGGTCACGTCTCCTGATGAAGAAGAAGCAACCCGAACGTATTTTCCACAGCCGTGTCAGGCCTGATGGAAGCAAAGGCGATGCGTGTACGTGGGTATTATGGGAACATCCTTGAGCGGACTGGATGCGGCAATGGCAGTGCTATTCACATGCGTTCATGAATGCATCCTGGTGACCAAATAAGCAGTTCATGCAGCCAGTGATGATCGCTACCAGCAGCTTGGGATGGCTCCGCTGACCAGGACAGAGCCGTTGCTAAAGTCGTGATGACAATGCGCCGCAAACCGCGAAAGCCGGATTATCTGGATATCCCAGCATTCCCTGCCGTAATGCAAGCTGATTAAGAATTGCTGGATTTGGGTTCGAGGCTCTTGTGCTAAACTGGCCGCCCCGAGTATAGTATATTTCGGTTGGATAGGTAATTTGCGAGATAATAACTTGATCAAACCAAAGGACACTTAAACGTATCGTTCAGGCGACGGGTGTGCGGTTACATACCGCAAGAACAGTCACCCTGACGTTACGCCCTGCGCCGGCCACACCGGGGTCACTATCGTCGCACACTTGATAATCCACCGTAGAATTTCCCCGGCGAGCCAAATTGTGCGTATACCATGCTCTGGGTAAACGTGTCTGTCAACGGCATGATGTACGATTTCAACCGTAGAGCACCTCAATGCTGCTCTGCGGGCTTGGCATCGATAACATTGTTACAAAGTTAACGCGCCGGAAATCCCGTCACTGGACCGCAGCGCCGTCTCCGTTTTGATACCTGCTGCTTGACGCAGGTATCGACGAGATGAACTGGCGCCAAAAAATTGTTCGCATCAAAGAGACCTGTTCGTGTCGAAGATGGCGTAAGTGGGCTGAATTTATGCCGTACAATGGTTTTTCGCTGGATTTCACCTCTCTCAATCAAACAACAAGGAGGAGGAGGAAAAAAGAGAGAGATGAAACCCCAGCGAAAAACCATGTACGCTAAATTTCAGCCCACTTTCGCCATCTTCGACACGAACAGTCTCTTGATGCGAACAAATTTTTTGCGCAGTCAACTCGTCGTACCGGCGTCAAGCAGCAGGTATACCAAACGGTGCGCGCTTGCCGTCCATGATCGGGATTTCCGGGCCGCGTAACTTTCGTAAACATGTTATCGAGCCCAAGCCCGCGAGAGCGCATTGAGGTGCTCCTTACGGTTAAATCCGTACATCATGCTCGTGACCAGACACGTACAGAGCATGTATCACGCACAGAATTTGGCAATCCGGCCCCCGGGAAATCTACCGCTGGATTCAAGTCGGTGCGACGATAGTGACCCCGTGTTGGCCGGGCCGCCAGGGCGTAACGTCAAGGGGGTCCTTTTTCTTGCCCCCCGGTATGTAAAACCCGACACCCGTCCCCCCGCCCCTGAACGAAACCCCGTTTAAGTGGTTGTTTGGGAACAATCGTAATTTCTCCCATTTACTAAATCCCCCCCGAAACGGGGGGGGGTTTTAAACAA";
        let query_record =
            ReadInfo::new_fa_record("tt".to_string(), String::from_utf8(seq.to_vec()).unwrap());
        build_hp_tr_metric(
            &query_record,
            &aligners,
            &targetname2seq,
            &region2motif,
            &OupParams::default(),
            true,
        );
    }
}
