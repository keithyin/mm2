
use std::process::Command;



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