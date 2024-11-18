use std::process::Command;

use std::{fs, path, process};

pub fn sort(bam_file: &str, tag: &str, threads: Option<usize>) -> String {
    let threads = threads.unwrap_or(num_cpus::get_physical() / 2);
    let threads = if threads > 1 { threads } else { 1 };

    let out_bam = format!(
        "{}.sort_by_{}.bam",
        bam_file.rsplit_once(".").unwrap().0,
        tag
    );

    let mut cmd = Command::new("samtools");
    cmd.args([
        "sort",
        "-@",
        threads.to_string().as_str(),
        "-n",
        "-t",
        tag,
        "-o",
        out_bam.as_str(),
    ]);

    let oup = cmd.output().expect(&format!("sort {} error", bam_file));
    if !oup.status.success() {
        panic!(
            "sort {} error. {}",
            bam_file,
            String::from_utf8(oup.stderr).unwrap()
        );
    }

    out_bam
}

pub fn samtools_bai(bam_file: &str, force: bool, threads: Option<usize>) -> anyhow::Result<()> {
    let res_filepath = format!("{}.bai", bam_file);

    let threads = threads.unwrap_or(num_cpus::get_physical() / 2);
    let threads = if threads > 1 { threads } else { 1 };

    if force {
        if path::Path::new(&res_filepath).exists() {
            fs::remove_file(&res_filepath).expect(&format!("remove {} error", res_filepath));
        }
    }

    let mut index_cmd = process::Command::new("samtools");

    index_cmd.args([
        "index",
        "-@",
        threads.to_string().as_str(),
        bam_file,
        &res_filepath,
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

pub fn sort_by_coordinates(bam_file: &str, threads: Option<usize>) {
    let threads = threads.unwrap_or(num_cpus::get_physical() / 2);
    let threads = if threads > 1 { threads } else { 1 };

    let out_bam = format!("{}.tmp", bam_file);

    let mut cmd = Command::new("samtools");
    cmd.args([
        "sort",
        "-o",
        out_bam.as_str(),
        "-@",
        threads.to_string().as_str(),
        bam_file,
    ]);

    let oup = cmd.output().expect(&format!("sort {} error", bam_file));
    if !oup.status.success() {
        panic!(
            "sort {} error. {}",
            bam_file,
            String::from_utf8(oup.stderr).unwrap()
        );
    }

    fs::rename(out_bam, bam_file).unwrap();
}
