use rust_htslib::bam::record::{Cigar, CigarString};

pub fn mapping_cigar2htslib_cigar_str(mapping_cigar: &[(u32, u8)]) -> CigarString {
    let cigar_str = mapping_cigar
        .iter()
        .map(|&(cnt, op)| {
            let cur_cigar = match op {
                0 => Cigar::Match(cnt),
                1 => Cigar::Ins(cnt),
                2 => Cigar::Del(cnt),
                3 => Cigar::RefSkip(cnt),
                4 => Cigar::SoftClip(cnt),
                5 => Cigar::HardClip(cnt),
                6 => Cigar::Pad(cnt),
                7 => Cigar::Equal(cnt),
                8 => Cigar::Diff(cnt),
                v => panic!("invalid cigar op :{}", v),
            };
            cur_cigar
        })
        .collect::<Vec<_>>();
    CigarString(cigar_str)
}

pub fn mapping_cigar2str(mapping_cigar: &[(u32, u8)]) -> String {
    let cigar_str = mapping_cigar
        .iter()
        .map(|&(cnt, op)| {
            let cur_cigar = match op {
                0 => format!("{}M", cnt),
                1 => format!("{}I", cnt),
                2 => format!("{}D", cnt),
                3 => format!("{}N", cnt),
                4 => format!("{}S", cnt),
                5 => format!("{}H", cnt),
                6 => format!("{}P", cnt),
                7 => format!("{}=", cnt),
                8 => format!("{}X", cnt),
                v => panic!("invalid cigar op :{}", v),
            };
            cur_cigar
        })
        .collect::<Vec<_>>();
    cigar_str.join("")
}
