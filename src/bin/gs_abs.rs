
use clap::{self, Parser};
use gskits::{fastx_reader::fasta_reader::FastaFileReader, file_reader::{bed_reader::BedInfo, vcf_reader::VcfInfo}, gsbam::bam_record_ext::BamRecordExt, pbar};

#[derive(Parser)]
pub struct Cli {
    pub reffasta: String,
    pub aligned_bam: String,

    #[arg(long="hcregions")]
    pub hcregions: Option<String>,
    #[arg(long="hcvariants")]
    pub hcvariants: Option<String>,
    
    pub chrom: Option<String>,

}


use std::{collections::HashMap, fs, io::{BufReader, BufWriter, Write}, sync::Arc, thread};
use crossbeam;
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Cigar, Read};


struct RecordReplica {
    ch: usize,
    q_len: usize,
    passes: usize,
    rq: Option<f32>,
    cigars: Vec<bam::record::Cigar>,
    seq: String,
    align_ref_start: usize,
    align_ref_end: usize
}

impl RecordReplica {
    pub fn from_record(record: &bam::Record) -> Self {
        let cigars = record.cigar().take().0;
        let record_ext = BamRecordExt::new(record);
        Self { 
            ch: record_ext.get_ch().unwrap(), 
            q_len: record_ext.get_seq().len(), 
            passes: record_ext.get_np().unwrap_or(1), 
            rq: record_ext.get_rq(), 
            cigars: cigars,
            seq: record_ext.get_seq(),
            align_ref_start: record.reference_start() as usize,
            align_ref_end: record.reference_end() as usize
        }
    }

    pub fn get_ch(&self) -> usize {
        self.ch
    }

    pub fn get_seq_len(&self) -> usize {
        self.q_len
    }

    pub fn get_np(&self) -> usize {
        self.passes
    }

    pub fn get_rq(&self) -> Option<f32> {
        self.rq
    }

    pub fn get_cigars(&self) -> &Vec<bam::record::Cigar> {
        &self.cigars
    }

    pub fn get_seq(&self) -> &str {
        &self.seq
    }

    pub fn get_aligned_ref_start(&self) -> usize {
        self.align_ref_start
    }

    pub fn get_aligned_ref_end(&self) -> usize {
        self.align_ref_end
    }
}

pub struct FastaFile {
    ref_name2seq: HashMap<String, String>
}

impl FastaFile {
    
    pub fn new(filepath: &str) -> Self {
        let reader = FastaFileReader::new(filepath.to_string());
        Self { ref_name2seq: reader.into_iter().map(|record| (record.name, record.seq)).collect::<HashMap<String, String>>() }
    }

    pub fn get_ref_seq(&self, refname: &str) -> Option<&str> {
        self.ref_name2seq.get(refname).and_then(|v| Some(v.as_str()))
    }

    pub fn get_ref_name2seq(&self) -> &HashMap<String, String> {
        &self.ref_name2seq
    }

}

#[derive(Debug)]
struct Stat {
    ch: usize,
    passes: usize,
    q_len: usize,
    m: usize,
    mm: usize,
    non_h_del: usize,
    h_del: usize,
    non_h_ins: usize,
    h_ins: usize,
    rq: Option<f32>,
    ignore_bps: usize,
    ref_start: usize,
    ref_end: usize,
    mm_ref_positions: Vec<usize>,
    ins_ref_positions: Vec<usize>,
    del_ref_positions: Vec<usize>
}

impl Stat {
    fn new(ch: usize, q_len: usize, passes: usize, ref_start: usize, ref_end: usize, rq: Option<f32>) -> Self{
        Self { 
            ch: ch, 
            passes: passes,
            q_len: q_len, 
            m: 0, 
            mm: 0, 
            non_h_del: 0, 
            h_del: 0, 
            non_h_ins: 0, 
            h_ins: 0,
            rq: rq,
            ignore_bps: 0,
            ref_start: ref_start,
            ref_end: ref_end,
            mm_ref_positions: vec![],
            ins_ref_positions: vec![],
            del_ref_positions: vec![]
        }
    }
    
    fn ins_bp(&self) -> usize {
        self.non_h_ins + self.h_ins
    }

    fn del_bp(&self) -> usize {
        self.h_del + self.non_h_del
    }

    fn align_span(&self) -> usize {
        self.ins_bp() + self.del_bp() + self.m + self.mm
    }

    fn concordance(&self) -> f32 {
        (self.m as f32) / self.align_span() as f32
    }

    // qv is capped by 60
    fn concordance_qv(&self) -> f32 {
        let mut err_rate = 1.0 - self.concordance();
        err_rate = if err_rate < 1e-6 {1e-6} else {err_rate};
        -10. * err_rate.log10()
    }

    // qv coverage
    fn query_converage(&self) -> f32 {
        (self.ins_bp() + self.m + self.mm) as f32 / self.q_len as f32
    }
}


fn get_hc_regions(region_file: Option<&str>) -> Option<BedInfo>{
    if let Some(region_file) = region_file {
        let reader = fs::File::open(region_file).unwrap();
        let mut buf_reader = BufReader::new(reader);
        Some(BedInfo::new(&mut buf_reader))
    } else {
        None
    }
    
}


fn get_hcvariants(vcf_file: Option<&str>) -> Option<VcfInfo>{
    if let Some(vcf_file) = vcf_file {
        let reader = fs::File::open(vcf_file).unwrap();
        let mut buf_reader = BufReader::new(reader);
        Some(VcfInfo::new(&mut buf_reader))
    } else {
        None
    }
}


fn do_m_mm_stat(counter: &mut usize, start_pos: usize, n: usize, ref_name: &str, hc_regions: &Option<BedInfo>, hc_variants: &Option<VcfInfo>) -> (usize, Vec<usize>){
    let mut positions = vec![];
    let mut ignore_bps = 0_usize;
    if hc_regions.is_none() && hc_variants.is_none(){
        *counter += n as usize;
    } else {
        for shift in 0..n {
            let cur_pos = start_pos + shift;
            if valid_point(ref_name, cur_pos, &hc_regions, &hc_variants) {
                *counter += 1;
                positions.push(cur_pos);
            } else {
                ignore_bps += 1;
            }
        }
    }
    (ignore_bps, positions)
}


fn valid_insertion(ref_pos: usize, ref_name: &str, hc_regions: &Option<BedInfo>, hc_variants: &Option<VcfInfo>) -> bool {
    let mut valid = true;
    if let Some(hc_regions) = hc_regions {
        valid &= hc_regions.point_within_region(ref_name, ref_pos);
    }

    if let Some(hc_variants) = hc_variants {
        valid &= !hc_variants.point_hit(ref_name, ref_pos);
    }

    if ref_pos > 0 {
        if let Some(hc_regions) = hc_regions {
            valid &= hc_regions.point_within_region(ref_name, ref_pos - 1);
        }
    
        if let Some(hc_variants) = hc_variants {
            valid &= !hc_variants.point_hit(ref_name, ref_pos - 1);
        }
    }

    valid
}


fn valid_point(ref_name: &str, position: usize, hc_regions: &Option<BedInfo>, hc_variants: &Option<VcfInfo>) -> bool {

    let mut valid = true;
    if let Some(hc_regions) = hc_regions {
        valid &= hc_regions.point_within_region(ref_name, position);
    }

    if let Some(hc_variants) = hc_variants {
        valid &= !hc_variants.point_hit(ref_name, position);
    }
    
    valid  

}


fn stat_record_core(record: RecordReplica, ref_name: &str, references: &HashMap<String, String>, hc_regions: &Option<BedInfo>, hc_variants: &Option<VcfInfo>) -> Stat {
    let ref_seq = references.get(ref_name).expect(&format!("refname:'{}' not found. valid refnames are:'{:?}'", ref_name, references.keys()));
    let query_seq = record.get_seq();

    let mut rpos_cursor = record.get_aligned_ref_start();
    let mut qpos_cursor = 0_usize;
    let mut stat = Stat::new(record.get_ch(), record.get_seq_len(), record.get_np(),record.get_aligned_ref_start(), record.get_aligned_ref_end(), record.get_rq());
    // println!("cigar:{:?}", record.get_cigars());
    for cig in record.get_cigars(){
        match cig {
            &Cigar::Equal(n) => {
                let (ignore_bps, _) = do_m_mm_stat(&mut stat.m, rpos_cursor, n as usize, ref_name, &hc_regions, &hc_variants);
                stat.ignore_bps += ignore_bps;
            },
            &Cigar::Diff(n) => {
                let (ignore_bps, mm_ref_pos) = do_m_mm_stat(&mut stat.mm, rpos_cursor, n as usize, ref_name, &hc_regions, &hc_variants);
                stat.ignore_bps += ignore_bps;
                stat.mm_ref_positions.extend(mm_ref_pos.into_iter());
            },
            &Cigar::Match(_) => panic!("Match not valid"),
            _ => (),
        }

        match cig {
            &Cigar::Ins(n) => {
                if valid_insertion(rpos_cursor, ref_name, &hc_regions, &hc_variants) {
                    let n = n as usize;
                    let mut is_hp = false;
                    let cur_bp_seq = &query_seq[qpos_cursor..qpos_cursor+n];
                    if rpos_cursor > 0 {
                        let cur_rseq = vec![&ref_seq.as_str()[rpos_cursor-1..rpos_cursor]; n].join("");
                        assert_eq!(cur_bp_seq.len(), cur_rseq.len());
                        is_hp |= cur_rseq.as_str() == cur_bp_seq;
                    }
                    if rpos_cursor < ref_seq.len() {
                        let cur_rseq = vec![&ref_seq.as_str()[rpos_cursor..rpos_cursor+1]; n].join("");
                        assert_eq!(cur_bp_seq.len(), cur_rseq.len());
                        is_hp |= cur_rseq.as_str() == cur_bp_seq;
                    }

                    if is_hp {
                        stat.h_ins += n;
                    } else {
                        stat.non_h_ins += n;
                    }
                    stat.ins_ref_positions.push(rpos_cursor);

                } else {
                    stat.ignore_bps += n as usize;
                }
            },
            &Cigar::Del(n) => {
                let n = n as usize;
                for shift in 0..n {
                    let cur_r_pos = rpos_cursor + shift;
                    if valid_point(ref_name, cur_r_pos, &hc_regions, &hc_variants) {
                        let mut is_hp = false;
                        let cur_ref_bp = ref_seq.as_bytes()[cur_r_pos];
                        if cur_r_pos > 0 {
                            is_hp |= cur_ref_bp == ref_seq.as_bytes()[cur_r_pos-1];
                        }
                        if cur_r_pos < ref_seq.len() {
                            is_hp |= cur_ref_bp == ref_seq.as_bytes()[cur_r_pos + 1];
                        }
                        if is_hp {
                            stat.h_del += 1;
                        } else {
                            stat.non_h_del += 1;
                        }
                        stat.del_ref_positions.push(cur_r_pos);
                    } else {
                        stat.ignore_bps += 1;
                    }
                }
            },

            _ => (),
        }

        match cig {
            &Cigar::Equal(n) | &Cigar::Diff(n)=> {
                rpos_cursor += n as usize;
                qpos_cursor += n as usize;
            },
            &Cigar::Del(n) | &Cigar::RefSkip(n)=> rpos_cursor += n as usize,
            &Cigar::Ins(n) | &Cigar::SoftClip(n) => qpos_cursor += n as usize,
            _ => panic!("? cigar:{:?}. {}", cig, cig),
        }
    }

    stat
}


pub fn bam_concordance(args: &Cli) -> anyhow::Result<()> {
    let bed_filepath = args.hcregions.clone();
    let bed_thread = thread::spawn(
        move || get_hc_regions(bed_filepath.as_ref().and_then(|v| Some(v.as_str()))));

    let vcf_filepath = args.hcvariants.clone();
    let vcf_thread = thread::spawn(
        move || get_hcvariants(vcf_filepath.as_ref().and_then(|v| Some(v.as_str()))));
    
    
    let hc_regions = bed_thread.join().unwrap();
    let hc_variants = vcf_thread.join().unwrap();
    let reffasta = FastaFile::new(&args.reffasta);
    
    
    thread::scope(|thread_scope| {
        let hc_regions = &hc_regions;
        let hc_variants = &hc_variants;
        let reffasta = &reffasta;

        let (record_sender, record_receiver) = crossbeam::channel::bounded(200);
        thread_scope.spawn(move|| {
            let mut tid2refname = HashMap::new();
            let mut reader = bam::Reader::from_path(&args.aligned_bam).unwrap();
            reader.set_threads(10).unwrap();
            let header = bam::Header::from_template(reader.header());
            let header_view = bam::HeaderView::from_header(&header);
            for record in reader.records() {
                let record = record.unwrap();
                let tid = record.tid() as u32;
                if !tid2refname.contains_key(&tid) {
                    tid2refname.insert(tid, Arc::new(String::from_utf8(header_view.tid2name(tid).to_vec()).unwrap()));
                }

                record_sender.send(
                    (RecordReplica::from_record(&record), 
                        tid2refname.get(&tid).unwrap().clone(),
                        reffasta,
                        hc_regions,
                        hc_variants
                    )
                    ).unwrap();
            }
        });

        let stat_threads = 8 as usize;
        let (stat_sender, stat_receiver) = crossbeam::channel::bounded(200);
        for _ in 0..stat_threads {
            let record_receiver_ = record_receiver.clone();
            let stat_sender_ = stat_sender.clone();
            thread_scope.spawn(move|| {
                for info in record_receiver_ {
                    let stat = stat_record_core(info.0, info.1.as_str(), reffasta.get_ref_name2seq(), hc_regions, hc_variants);
                    stat_sender_.send(stat).unwrap();
                }
            });
        }
        drop(record_receiver);
        drop(stat_sender);

        // write result
        let csv_filepath = format!("{}.metric.csv", args.aligned_bam.rsplit_once(".").unwrap().0);
        let csv_writer = fs::File::create(&csv_filepath).unwrap();
        let mut buf_csv_writer = BufWriter::new(csv_writer);

        let pbar = pbar::get_spin_pb(format!("Processing {}", args.aligned_bam), pbar::DEFAULT_INTERVAL);

        writeln!(&mut buf_csv_writer, 
            "channel_id\treadLengthBp\tsubreadPasses\t\
            queryConverage\tpredictedConcordance\t\
            concordance\tconcordanceQv\tmatchBp\tmismatchBp\tnonHpInsertionBp\t\
            nonHpDeletionBp\thpInsertionBp\thpDeletionBp\tignoreBp\t\
            refStart\trefEnd\tmmRefPositions\tinsRefPositions\tdelRefPositions").unwrap();
        for stat in stat_receiver {
            writeln!(&mut buf_csv_writer, "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t\"{16}\"\t\"{17}\"\t\"{18}\"", 
                stat.ch, stat.q_len, stat.passes, 
                stat.query_converage(), stat.rq.unwrap_or(0.), stat.concordance(), 
                stat.concordance_qv(), stat.m, stat.mm, stat.non_h_ins, stat.non_h_del, 
                stat.h_ins, stat.h_del, stat.ignore_bps,
                stat.ref_start, stat.ref_end, stat.mm_ref_positions.iter().map(|v| v.to_string()).collect::<Vec<String>>().join(","),
                stat.ins_ref_positions.iter().map(|v| v.to_string()).collect::<Vec<String>>().join(","),
                stat.del_ref_positions.iter().map(|v| v.to_string()).collect::<Vec<String>>().join(","),
            ).unwrap();
            pbar.inc(1);
        }
        pbar.finish();

    });

    Ok(())
}


fn main() {
    let args = Cli::parse();

    bam_concordance(&args).unwrap();
    
}

#[cfg(test)]
mod test {
    use std::{fs, io::BufReader};

    use gskits::file_reader::vcf_reader::VcfInfo;
    use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};


    use crate::{Cli, stat_record_core, FastaFile, RecordReplica};

    use super::bam_concordance;
    #[test]
    fn test_stat_record_core() {
        let aligned_bam_path = "/data/ccs_data/ccs_eval2024q3/Mitochondria/subread/9948_subreads.smc_all_reads.aligned.bam";
        let mut bam_reader = bam::Reader::from_path(aligned_bam_path).unwrap();
        let mut cnt = 0;
        for (i, record) in bam_reader.records().enumerate() {
            let record = record.unwrap();
            if record.reference_start() == 0 {
                continue;
            }
            cnt += 1;
            println!("{:?}, ref_start:{}", record.cigar(), record.reference_start());
            if cnt > 10 {
                break;
            }
        }

    }

    #[test]
    fn test_logic_comp() {
        let mut a = true;
        a &= true;
        println!("{}", a);
        a &= false;
        println!("{}", a);
    }


    #[test]
    fn test_stat_record_core2() {
        let aligned_bam_path = "/data/ccs_data/ccs_eval2024q3/Ludaopei/subread/Sample3_subreads.aligned.bam";
        let ref_fasta = "/data/ccs_data/ccs_eval2024q3/Ludaopei/ref/Sample3.fa";
        let mut vcf_bufreader = BufReader::new(fs::File::open("/data/ccs_data/ccs_eval2024q3/Ludaopei/ref/Sample3.vcf").unwrap());
        let hc_variants = Some(VcfInfo::new(&mut vcf_bufreader));
        let ref_fasta = FastaFile::new(ref_fasta);
        let mut bam_reader = bam::Reader::from_path(aligned_bam_path).unwrap();
        let mut cnt = 0;
        for (i, record) in bam_reader.records().enumerate() {
            let record = record.unwrap();
            println!("{:?}, ref_start:{}, ref_end:{}", record.cigar(), record.reference_start(), record.reference_end());
            let record = RecordReplica::from_record(&record);
            let stat_info = stat_record_core(record, "BCR-ABL1-P210-e14a2", ref_fasta.get_ref_name2seq(), &None, &hc_variants);
            println!("{:?}", stat_info);
            cnt += 1;
            if cnt > 10 {
                break;
            }
        }
    }

    #[test]
    fn test_bam_concordance() {
        let reffasta = "/data/ccs_data/ccs_eval2024q3/Ludaopei/ref/Sample3.fa".to_string();
        let aligned_bam = "/data/ccs_data/ccs_eval2024q3/Ludaopei/subread/Sample3_subreads.aligned.bam".to_string();
        let hcvariants = "/data/ccs_data/ccs_eval2024q3/Ludaopei/ref/Sample3.vcf".to_string();
        let args = &Cli { reffasta: reffasta, aligned_bam: aligned_bam, hcregions: None, hcvariants: Some(hcvariants), chrom: None };
        bam_concordance(args).unwrap();
    }

    #[test]
    fn test_print() {

        println!("aa\
        bb");

    }


}



/*
docker run --rm -v`pwd`:/data -v/data/adapters:/adapters \
    192.168.3.44:5000/algo/adacus:release1.0.0_smc4.1.0_adapter_remover1.0.6 \
    adapter_remover_multi \
    -i /data/B4-20240812_Sync_Y0004_09_H01_Run0002_called.bam \
    -o /data/B4-20240812_Sync_Y0004_09_H01_Run0002_called.adapter.bam \
    -a /adapters/b4.fa \
    -ia /adapters/b4_init.fa \
    -mer 0.3 \
    -imer 0.2 \
    -mp 50
 */

