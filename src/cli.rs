
use std::str::FromStr;

use clap::{self, Parser, Subcommand, Args};

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(long="threads")]
    pub threads: Option<usize>,
    
    #[arg(long="preset", default_value_t=String::from_str("map-ont").unwrap())]
    pub preset: String,



    #[command(subcommand)]
    pub commands: Commands

}

#[derive(Debug, Subcommand, Clone)]
pub enum Commands {
    Index(IndexArgs),
    R2R(ReadsToRefAlignArgs),
    S2S(Sbreads2SmcAlignArgs)
}

#[derive(Debug, Args, Clone, Copy)]
pub struct IndexArgs {
    #[arg(long, help="minimizer kmer")]
    kmer: Option<usize>,

    #[arg(long, help="minimizer window size")]
    wins: Option<usize>,
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
    pub oup_args: OupArgs

}

#[derive(Debug, Args, Clone)]
pub struct Sbreads2SmcAlignArgs {
    #[command(flatten)]
    io_args: IoArgs,

    #[command(flatten)]
    index_args: IndexArgs,

    #[command(flatten)]
    map_args: MapArgs,

    #[command(flatten)]
    align_args: AlignArgs,

    #[command(flatten)]
    oup_args: OupArgs
}


#[derive(Debug, Args, Clone)]
pub struct IoArgs {

    #[arg(short='q', help="query file paths, 
    if multiple query_filepath are provided, 
    their query_names will be rewriten to ___0, ___1, and so on, 
    based on the order of the filenames", required=true)]
    pub query: Vec<String>,

    #[arg(long="target", group="target")]
    pub target: Option<String>,
    #[arg(long="indexedTarget", group="target")]
    pub indexed_target: Option<String>,

    #[arg(short='p', help="output a file named ${p}.bam")]
    pub prefix: String
}

impl IoArgs {
    pub fn get_oup_path(&self) -> String {
        format!("{}.bam", self.prefix)
    }
}

#[derive(Debug, Args, Clone)]
pub struct MapArgs {
    
}

#[derive(Debug, Args, Clone)]
pub struct AlignArgs {
    #[arg(short='m', help="matching_score>=0, recommend 2")]
    matching_score: Option<i32>,

    #[arg(short='M', help="mismatch_penalty >=0, recommend 4")]
    mismatch_penalty: Option<i32>,

    #[arg(short='o', help="gap_open_penalty >=0, recommend 4,24")]
    gap_open_penalty: Option<String>,

    #[arg(short='e', help="gap_extension_penalty >=0, recommend 2,1")]
    gap_extension_penalty: Option<String>,
}


#[derive(Debug, Args, Clone)]
pub struct OupArgs {

    #[arg(long="noSeco", help="discard secondary alignment")]
    discard_secondary: bool,

    #[arg(long="noSupp", help="discard supplementary alignment")]
    discard_supplementary: bool,

    #[arg(long="noMD", help="do not add MD to the final bam record")]
    no_md: bool,

    #[arg(long="noBai", help="do not generate bai file")]
    no_bai: bool,

    #[arg(long="noExtraInfo", help="don't append ch/np/iy... to result bam")]
    no_extra_info: bool,

    #[arg(long="oupIyT", default_value_t=0.0, help="remove the record from the result bam file when the identity < identity_threshold")]
    oup_identity_threshold: f32,

    #[arg(long="oupCovT", default_value_t=0.0, help="remove the record from the result bam file when the coverage < coverage_threshold")]
    oup_coverage_threshold: f32,
}