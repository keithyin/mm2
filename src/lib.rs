pub mod cli;
use cli::{AlignArgs, IndexArgs, MapArgs};
use minimap2::{Aligner, Preset};

fn build_aligner(preset: &str, index_args: &IndexArgs, map_args: &MapArgs, align_args: &AlignArgs) -> Aligner{

    let mut aligner = Aligner::builder();
    aligner = match preset {
        "map-ont" => aligner.map_ont(),
        ps => panic!("invalid preset {}", ps)
    };


    aligner

}


pub fn index_ref_file() {}

pub fn query2ref_align(
    query_files: &Vec<&str>,
    ref_file: Option<&str>,
    indexed_ref_file: Option<&str>,
    oup_filename: &str,
    mut aligner: Aligner,
    threads: usize
) {
    if ref_file.is_none() && indexed_ref_file.is_none() {
        panic!("ref_file and indexed_ref_file cannot all be none");
    }

    if indexed_ref_file.is_some() {
        
    } else {
        // aligner.with_index_threads(threads).with_seq(seq)
    }

    


}

pub fn subreads2smc_align() {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        println!("hello world");
    }
}
