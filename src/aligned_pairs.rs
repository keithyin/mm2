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
pub trait TAlignedPairs {
    fn target_start(&self) -> usize;
    fn target_end(&self) -> usize;
    fn query_start(&self) -> usize; // 原始序列的开始
    fn query_end(&self) -> usize; // 原始序列的结束
    fn is_reverse(&self) -> bool;
    fn query_len(&self) -> usize;
    fn cigars(&self) -> &[(u32, u8)]; // 如果是 reverse，那么是 reverse 后的 cigar。

    fn aligned_pairs(&self) -> impl Iterator<Item = (Option<usize>, Option<usize>)> {
        let mut query_cursor = if self.is_reverse() {
            self.query_len() - self.query_end()
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
                            .map(move |shift| (Some(old_query_cursor + shift), None)),
                    )
                        as Box<dyn Iterator<Item = (Option<usize>, Option<usize>)>>
                }
                2 => {
                    // del
                    let old_target_cursor = target_cursor;
                    target_cursor += cnt;
                    Box::new(
                        (0..cnt)
                            .into_iter()
                            .map(move |shift| (None, Some(old_target_cursor + shift))),
                    )
                        as Box<dyn Iterator<Item = (Option<usize>, Option<usize>)>>
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
                        )
                    }))
                        as Box<dyn Iterator<Item = (Option<usize>, Option<usize>)>>
                }
                op => panic!("not a valid op:{}", op),
            }
        })
    }
}

// impl TAlignedPairs for AlignInfo {
//     fn cigars(&self) -> &[(u32, u8)] {}
// }
