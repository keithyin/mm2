use minimap2::Mapping;

pub struct MappingExt<'a>(pub &'a Mapping);

impl<'a> MappingExt<'a> {
    pub fn identity(&self) -> f32 {
        if let Some(aln) = &self.0.alignment {
            let matched = aln
                .cigar
                .as_ref()
                .unwrap()
                .iter()
                .filter(|(_, op)| *op == 7)
                .map(|(cnt, _)| *cnt)
                .sum::<u32>();

            let tot = aln
                .cigar
                .as_ref()
                .unwrap()
                .iter()
                .map(|(cnt, _)| *cnt)
                .sum::<u32>();
            matched as f32 / tot as f32
        } else {
            0.0
        }
    }

    pub fn query_coverage(&self) -> f32 {
        let qlen = self.0.query_len.unwrap().get();
        let qstart = self.0.query_start;
        let qend = self.0.query_end;

        (qend - qstart) as f32 / qlen as f32
    }
}
