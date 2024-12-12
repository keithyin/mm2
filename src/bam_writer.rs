use std::collections::HashMap;

use crossbeam::channel::Receiver;
use gskits::{
    gsbam::bam_record_ext::{BamRecord, BamRecordExt, BamWriter},
    pbar::{self, DEFAULT_INTERVAL},
    utils::command_line_str,
};
use rust_htslib::bam::{header::HeaderRecord, record::Aux, Header};

use crate::{params::OupParams, AlignResult};

pub struct BamOupArgs {
    pub iy_t: f32,
    pub ec_t: f32,
    pub no_seco: bool,
    pub no_supp: bool,
    pub no_ma: bool
}
impl BamOupArgs {

    pub fn new(iy_t: f32, ec_t: f32, no_seco: bool, no_supp: bool, no_ma: bool) -> Self {
        Self{
            iy_t,
            ec_t,
            no_seco,
            no_supp,
            no_ma
        }
    }

    pub fn valid(&self, record: &BamRecord) -> bool {
        if self.no_seco && record.is_secondary() {
            return false;
        }

        if self.no_supp && record.is_supplementary() {
            return false;
        }

        let record_ext = BamRecordExt::new(&record);
        let iy = record_ext.compute_identity();
        let ec = record_ext.compute_query_coverage();

        if iy < self.iy_t {
            return false;
        }

        if ec < self.ec_t {
            return false;
        }

        true
    }
}


/// targetidx: &HashMap<target_name, (idx, target_len)>
pub fn write_bam_worker(
    recv: Receiver<AlignResult>,
    target_idx: &HashMap<String, (usize, usize)>,
    o_path: &str,
    oup_args: &OupParams,
    pg_name: &str,
    version: &str,
    enable_pb: bool,

) {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.5");
    hd.push_tag(b"SO", "unknown");
    header.push_record(&hd);

    let mut hd = HeaderRecord::new(b"PG");
    hd.push_tag(b"ID", pg_name)
        .push_tag(b"PN", pg_name)
        .push_tag(b"CL", &command_line_str())
        .push_tag(b"VN", version);
    header.push_record(&hd);

    let mut targets = target_idx.iter().collect::<Vec<_>>();
    targets.sort_by_key(|tup| tup.1 .0);

    for target in targets {
        let mut hd = HeaderRecord::new(b"SQ");
        hd.push_tag(b"SN", &target.0);
        hd.push_tag(b"LN", target.1 .1);
        header.push_record(&hd);
    }

    let pb = if enable_pb {
        Some(pbar::get_spin_pb(
            format!("{}: writing alignment result to {}", pg_name, o_path),
            DEFAULT_INTERVAL,
        ))
    } else {
        None
    };

    let mut writer = BamWriter::from_path(o_path, &header, rust_htslib::bam::Format::Bam).unwrap();
    writer.set_threads(4).unwrap();

    for align_res in recv {
        pb.as_ref().unwrap().inc(1);
        for mut record in align_res.records {
            if oup_args.valid(&record) {
                let record_ext = BamRecordExt::new(&record);
                let iy = record_ext.compute_identity();
                let ec = record_ext.compute_query_coverage();

                record.push_aux(b"iy", Aux::Float(iy)).unwrap();
                record.push_aux(b"ec", Aux::Float(ec)).unwrap();

                writer.write(&record).unwrap();
            }
        }
    }
    pb.as_ref().unwrap().finish();
}
