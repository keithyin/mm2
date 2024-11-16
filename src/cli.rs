
use std::str::FromStr;

use clap::{self, Parser, Subcommand, Args};
use minimap2::Aligner;

pub trait TOverrideAlignerParam {
    
    fn modify_aligner(&self, aligner: &mut Aligner);

}

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
    pub kmer: Option<usize>,

    #[arg(long, help="minimizer window size")]
    pub wins: Option<usize>,
}

impl TOverrideAlignerParam for IndexArgs {
    fn modify_aligner(&self, aligner: &mut Aligner) {
        if let Some(k) = self.kmer {
            aligner.idxopt.k = k as i16;
        }
        if let Some(w) = self.wins {
            aligner.idxopt.w = w as i16;
        }
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

    #[arg(short='p', help="output a file named ${p}.bam", required=true)]
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

impl TOverrideAlignerParam for  MapArgs {
    fn modify_aligner(&self, aligner: &mut Aligner) {
        
    }
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

impl TOverrideAlignerParam for AlignArgs {
    fn modify_aligner(&self, aligner: &mut Aligner) {
        if let Some(m) = self.matching_score {
            aligner.mapopt.a = m;
        }

        if let Some(mm) = self.mismatch_penalty {
            aligner.mapopt.b = mm;
        }

        if let Some(gap_o) = &self.gap_open_penalty {
            if gap_o.contains(",") {
                let (o1, o2) = gap_o.rsplit_once(".").unwrap();
                aligner.mapopt.q = o1.parse::<i32>().unwrap();
                aligner.mapopt.q2 = o2.parse::<i32>().unwrap();
            } else {
                aligner.mapopt.q = gap_o.parse::<i32>().unwrap();
            }
        }

        if let Some(gap_e) = &self.gap_extension_penalty {
            if gap_e.contains(",") {
                let (e1, e2) = gap_e.rsplit_once(".").unwrap();
                aligner.mapopt.e = e1.parse::<i32>().unwrap();
                aligner.mapopt.e2 = e2.parse::<i32>().unwrap();
            } else {
                aligner.mapopt.e = gap_e.parse::<i32>().unwrap();
            }
        }
    }
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

impl TOverrideAlignerParam for OupArgs {
    fn modify_aligner(&self, aligner: &mut Aligner) {
        if self.discard_secondary {
            aligner.mapopt.flag &=  !0x1000000000;
        }
    }
}