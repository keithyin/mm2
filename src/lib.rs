pub mod cli;
mod samtools;
use cli::{AlignArgs, IndexArgs, MapArgs};
use minimap2::{Aligner, Preset};
use rust_htslib::bam::Read;

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

pub fn subreads2smc_align(subreads_bam: &str, smc_bam: &str) {

    let sorted_subreads_bam = samtools::sort(subreads_bam, "ch");
    let sorted_smc_bam = samtools::sort(smc_bam, "ch");

    let mut smc_bam_reader = rust_htslib::bam::Reader::from_path(&sorted_smc_bam).unwrap();
    smc_bam_reader.set_threads(5).unwrap();
    let mut subreads_bam_reader = rust_htslib::bam::Reader::from_path(&sorted_subreads_bam).unwrap();
    subreads_bam_reader.set_threads(5).unwrap();

    
    let mut smc_records = smc_bam_reader.records();
    let mut subreads_records = subreads_bam_reader.records();

    loop {
        
        let smc_record = smc_records.next();
        if smc_record.is_none() {
            break;
        }

        let smc_record = smc_record.unwrap().unwrap();


        let mut sbr_record = None;

        loop {
            if sbr_record.is_none() {
                sbr_record = subreads_records.next();
            }

            let sbr = sbr_record.unwrap().unwrap();
            

            sbr_record = subreads_records.next();
        }

    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        println!("hello world");
    }
}
