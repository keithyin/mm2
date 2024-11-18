use rust_htslib::bam::{ext::BamRecordExtensions, Record};

pub fn record2str(record: &Record) -> String {
    format!(
        "refstart:{}, refend:{}, seq:{}",
        record.reference_start(),
        record.reference_end(),
        String::from_utf8(record.seq().as_bytes()).unwrap()
    )
}
