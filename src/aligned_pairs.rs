use crate::align_processor::AlignInfo;

/*

0 => Cigar::Match(cnt),
                1 => Cigar::Ins(cnt),
                2 => Cigar::Del(cnt),
                3 => Cigar::RefSkip(cnt),
                4 => Cigar::SoftClip(cnt),
                5 => Cigar::HardClip(cnt),
                6 => Cigar::Pad(cnt),
                7 => Cigar::Equal(cnt),
                8 => Cigar::Diff(cnt),

*/

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AlignOp {
    Match(u32),
    Ins(u32),
    Del(u32),
    RefSkip(u32),
    SoftClip(u32),
    HardClip(u32),
    Pad(u32),
    Equal(u32),
    Diff(u32),
}

impl From<u8> for AlignOp {
    fn from(value: u8) -> Self {
        match value {
            0 => AlignOp::Match(1),
            1 => AlignOp::Ins(1),
            2 => AlignOp::Del(1),
            3 => AlignOp::RefSkip(1),
            4 => AlignOp::SoftClip(1),
            5 => AlignOp::HardClip(1),
            6 => AlignOp::Pad(1),
            7 => AlignOp::Equal(1),
            8 => AlignOp::Diff(1),
            _ => panic!("invalid code. {}", value),
        }
    }
}

pub trait TAlignedPairs {
    fn target_start(&self) -> usize;
    fn target_end(&self) -> usize;
    fn query_start(&self) -> usize; // 原始序列的开始
    fn query_end(&self) -> usize; // 原始序列的结束
    fn is_reverse(&self) -> bool;
    fn cigars(&self) -> &[(u32, u8)]; // 如果是 reverse，那么是 reverse 后的 cigar。

    fn aligned_pairs(
        &self,
        query_len: usize,
    ) -> impl Iterator<Item = (Option<usize>, Option<usize>, AlignOp)> {
        let mut query_cursor = if self.is_reverse() {
            query_len - self.query_end()
        } else {
            self.query_start()
        };

        let mut target_cursor = self.target_start();

        self.cigars().iter().copied().flat_map(move |(cnt, op)| {
            let cnt = cnt as usize;
            match op {
                0 => panic!("eqx required"),
                1 => {
                    // ins
                    let old_query_cursor = query_cursor;
                    query_cursor += cnt;
                    Box::new(
                        (0..cnt)
                            .into_iter()
                            .map(move |shift| (Some(old_query_cursor + shift), None, op.into())),
                    )
                        as Box<dyn Iterator<Item = (Option<usize>, Option<usize>, AlignOp)>>
                }
                2 => {
                    // del
                    let old_target_cursor = target_cursor;
                    target_cursor += cnt;
                    Box::new(
                        (0..cnt)
                            .into_iter()
                            .map(move |shift| (None, Some(old_target_cursor + shift), op.into())),
                    )
                        as Box<dyn Iterator<Item = (Option<usize>, Option<usize>, AlignOp)>>
                }
                7 | 8 => {
                    let old_query_cursor = query_cursor;
                    query_cursor += cnt;
                    let old_target_cursor = target_cursor;
                    target_cursor += cnt;

                    Box::new((0..cnt).into_iter().map(move |shift| {
                        (
                            Some(old_query_cursor + shift),
                            Some(old_target_cursor + shift),
                            op.into(),
                        )
                    }))
                        as Box<dyn Iterator<Item = (Option<usize>, Option<usize>, AlignOp)>>
                }
                op => panic!("not a valid op:{}", op),
            }
        })
    }
}

impl TAlignedPairs for AlignInfo {
    fn cigars(&self) -> &[(u32, u8)] {
        &self.cigar
    }

    fn is_reverse(&self) -> bool {
        !self.fwd
    }

    fn query_end(&self) -> usize {
        self.qend
    }
    fn query_start(&self) -> usize {
        self.qstart
    }
    fn target_end(&self) -> usize {
        self.rend
    }
    fn target_start(&self) -> usize {
        self.rstart
    }
}

#[cfg(test)]
mod test {
    use crate::{align_processor::AlignInfo, aligned_pairs::TAlignedPairs};

    #[test]
    fn test_align_info() {
        let align_info = AlignInfo {
            rstart: 0,
            rend: 10,
            qstart: 1,
            qend: 12,
            fwd: true,
            primary: true,
            cigar: vec![(4, 7), (2, 1), (4, 7), (1, 2), (1, 7)],
        };
        for align_info in align_info.aligned_pairs(12) {
            println!("align_info:{:?}", align_info);
        }
    }
}
