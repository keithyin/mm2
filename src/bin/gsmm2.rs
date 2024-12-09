use std::thread;

use gskits::fastx_reader::read_fastx;
use gskits::{
    fastx_reader::fasta_reader::FastaFileReader,
    samtools::{samtools_bai, sort_by_coordinates},
};
use mm2::params::{AlignParams, IndexParams, InputFilterParams, MapParams, OupParams};
use mm2::{
    align_worker, bam_writer::write_bam_worker, build_aligner, query_seq_sender,
    targets_to_targetsidx,
};

use std::str::FromStr;

use clap::{self, Args, Parser, Subcommand};

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(long="preset", default_value_t=String::from_str("map-ont").unwrap(), 
        help="read https://lh3.github.io/minimap2/minimap2.html for more details")]
    pub preset: String,

    #[command(subcommand)]
    pub commands: Commands,
}

#[derive(Debug, Subcommand, Clone)]
pub enum Commands {
    Index(IndexArgs),
    Align(ReadsToRefAlignArgs),
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
pub struct ReadsToRefAlignArgs {
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
    #[arg(long = "target", short = 't', group = "target_group")]
    pub target: Option<String>,
    /// mmi file
    #[arg(long = "indexedTarget", group = "target_group")]
    pub indexed_target: Option<String>,

    #[arg(
        short = 'p',
        help = "output a file named ${p}.bam",
        required = true,
        requires = "target_group"
    )]
    pub prefix: String,

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
    pub fn get_oup_path(&self) -> String {
        format!("{}.bam", self.prefix)
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
    #[arg(long = "noSeco", help = "discard secondary alignment")]
    pub discard_secondary: bool,

    #[arg(long = "noSupp", help = "discard supplementary alignment")]
    pub discard_supplementary: bool,

    #[arg(long = "noMar", help = "discard multiple alignment reads")]
    pub discard_multi_mapping_reads: bool,

    #[arg(long="oupIyT", default_value_t=-1.0, help="remove the record from the result bam file when the identity < identity_threshold")]
    pub oup_identity_threshold: f32,

    #[arg(long="oupCovT", default_value_t=-1.0, help="remove the record from the result bam file when the coverage < coverage_threshold")]
    pub oup_coverage_threshold: f32,
}

impl OupArgs {
    fn to_oup_params(&self) -> OupParams {
        let mut param = OupParams::new();
        param = param
            .set_discard_secondary(self.discard_secondary)
            .set_discard_supplementary(self.discard_supplementary)
            .set_oup_identity_threshold(self.oup_identity_threshold)
            .set_oup_coverage_threshold(self.oup_coverage_threshold)
            .set_discard_multi_align_reads(self.discard_multi_mapping_reads);
        param
    }
}

fn alignment(preset: &str, tot_threads: Option<usize>, args: &ReadsToRefAlignArgs) {
    let tot_threads = tot_threads.unwrap_or(num_cpus::get());
    assert!(tot_threads >= 10, "at least 10 threads are needed");

    let target_filename = args
        .io_args
        .target
        .as_ref()
        .expect("target need to be provided");

    let fa_iter = FastaFileReader::new(target_filename.to_string());
    let targets = read_fastx(fa_iter);
    let target2idx = targets_to_targetsidx(&targets);

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

    /*
    1. query_seq_sender
    2. align
    3. write to bam
    */
    thread::scope(|s| {
        let aligners = &aligners;
        let target2idx = &target2idx;
        let inp_filter_params = &inp_filter_params;
        let (qs_sender, qs_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            query_seq_sender(&args.io_args.query, qs_sender, inp_filter_params);
        });

        let align_threads = tot_threads - 8;
        let (align_res_sender, align_res_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..align_threads {
            let qs_recv_ = qs_recv.clone();
            let align_res_sender_ = align_res_sender.clone();
            s.spawn(move || {
                align_worker(
                    qs_recv_,
                    align_res_sender_,
                    aligners,
                    target2idx,
                    &oup_params,
                )
            });
        }
        drop(qs_recv);
        drop(align_res_sender);

        write_bam_worker(
            align_res_recv,
            target2idx,
            &args.io_args.get_oup_path(),
            &oup_params,
            "gsmm2",
            env!("CARGO_PKG_VERSION"),
            true,
        );
    });
    sort_by_coordinates(&args.io_args.get_oup_path(), Some(tot_threads));
    samtools_bai(&args.io_args.get_oup_path(), true, Some(tot_threads)).unwrap();
}

fn main() {
    let args = Cli::parse();

    let preset = &args.preset;
    let align_threads = args.threads.clone();

    match args.commands {
        Commands::Align(ref args) => {
            alignment(preset, align_threads, args);
        }
        _ => panic!("not implemented yet"),
    }
}
