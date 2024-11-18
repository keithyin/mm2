
use std::process::Command;

use std::{path, fs, process};

pub fn sort(bam_file: &str, tag: &str) -> String {

    let out_bam = format!("{}.sort_by_{}.bam", bam_file.rsplit_once(".").unwrap().0, tag);

    let mut cmd = Command::new("samtools");
    cmd.args(["sort", "-n", "-t", tag, "-o", out_bam.as_str()]);

    let oup = cmd.output().expect(&format!("sort {} error", bam_file));
    if !oup.status.success() {
        panic!("sort {} error. {}", bam_file, String::from_utf8(oup.stderr).unwrap());
    }

    out_bam

}

pub fn samtools_bai(bam_file: &str, force: bool) -> anyhow::Result<()> {
    let res_filepath = format!("{}.bai", bam_file);

    if force {
        if path::Path::new(&res_filepath).exists() {
            fs::remove_file(&res_filepath).expect(&format!("remove {} error", res_filepath));
        }
    }

    let mut index_cmd = process::Command::new("samtools");
    
    index_cmd.args([
        "index",
        "-@", "40",
        bam_file,
        &res_filepath
    ]);

    if let Ok(res) = index_cmd.status() {
        if !res.success() {
            tracing::error!("{} ::> {:?}", "RUN CMD ERROR", index_cmd);
            panic!("");
        }
    } else {
        tracing::error!("{} ::> {:?}", "RUN CMD ERROR", index_cmd);
        panic!("");
    }

    Ok(())
}