use std::ops::Deref;

use gskits::dna::reverse_complement;
use minimap2::{Mapping, Strand};

pub struct MappingExt<'a>(pub &'a Mapping);

impl<'a> Deref for MappingExt<'a> {
    type Target = Mapping;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a> MappingExt<'a> {
    pub fn identity(&self) -> f32 {
        if let Some(aln) = &self.0.alignment {
            let matched = self.num_matched();

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

    fn num_matched(&self) -> u32 {
        if let Some(aln) = &self.0.alignment {
            let matched = aln
                .cigar
                .as_ref()
                .unwrap()
                .iter()
                .filter(|&&(_, op)| {
                    assert!(op != 0, "eqx required");
                    op == 7
                })
                .map(|(cnt, _)| *cnt)
                .sum::<u32>();
            matched
        } else {
            0
        }
    }

    pub fn identity_without_long_indel(&self, indel_threshold: u32) -> f32 {
        let matched = self.num_matched();
        if matched == 0 {
            0.0
        } else {
            let tot = self
                .0
                .alignment
                .as_ref()
                .unwrap()
                .cigar
                .as_ref()
                .unwrap()
                .iter()
                .filter(|&&(cnt, op)| {
                    if op == 1 || op == 2 {
                        cnt < indel_threshold
                    } else {
                        true
                    }
                })
                .map(|&(cnt, _)| cnt)
                .sum::<u32>();

            matched as f32 / tot as f32
        }
    }

    pub fn identity_gap_compressed(&self) -> f32 {
        let matched = self.num_matched();
        if matched == 0 {
            0.0
        } else {
            let tot = self
                .0
                .alignment
                .as_ref()
                .unwrap()
                .cigar
                .as_ref()
                .unwrap()
                .iter()
                .map(|&(cnt, op)| if op == 1 || op == 2 { 1 } else { cnt })
                .sum::<u32>();

            matched as f32 / tot as f32
        }
    }

    pub fn query_coverage(&self) -> f32 {
        let qlen = self.0.query_len.unwrap().get();
        let qstart = self.0.query_start;
        let qend = self.0.query_end;
        (qend - qstart) as f32 / qlen as f32
    }

    pub fn target_coverage(&self) -> f32 {
        let t_len = self.0.target_len;
        let tstart = self.0.target_start;
        let tend = self.0.target_end;
        (tend - tstart) as f32 / t_len as f32
    }

    pub fn is_rev(&self) -> bool {
        match self.strand {
            Strand::Forward => false,
            Strand::Reverse => true,
        }
    }

    pub fn aligned_2_str(&self, target: &[u8], query: &[u8]) -> (String, String) {
        let sbr_target = &target[self.target_start as usize..self.target_end as usize];
        let sbr_target = String::from_utf8(sbr_target.to_vec()).unwrap();
        let sbr_target = sbr_target.as_str();
        let sbr_query = &query[self.query_start as usize..self.query_end as usize];
        let rev = self.is_rev();
        let sbr_query = if !rev {
            String::from_utf8(sbr_query.to_vec()).unwrap()
        } else {
            reverse_complement(String::from_utf8(sbr_query.to_vec()).unwrap().as_str())
        };

        let sbr_query = sbr_query.as_str();

        let mut aligned_target = String::new();
        let mut aligned_query = String::new();

        let mut target_start = 0;
        let mut query_start = 0;

        self.alignment
            .as_ref()
            .unwrap()
            .cigar
            .as_ref()
            .unwrap()
            .iter()
            .for_each(|&(cnt, code)| {
                let cnt = cnt as usize;
                match code {
                    0 | 7 | 8 => {
                        // match

                        aligned_target.push_str(&sbr_target[target_start..target_start + cnt]);
                        aligned_query.push_str(&sbr_query[query_start..query_start + cnt]);
                        target_start += cnt;
                        query_start += cnt;
                    }

                    1 => {
                        // ins
                        let ins = vec!["-"; cnt as usize].join("");

                        aligned_target.push_str(&ins);
                        aligned_query.push_str(&sbr_query[query_start..query_start + cnt]);
                        query_start += cnt;
                    }

                    2 => {
                        // del
                        let del = vec!["-"; cnt as usize].join("");

                        aligned_target.push_str(&sbr_target[target_start..target_start + cnt]);
                        aligned_query.push_str(&del);
                        target_start += cnt;
                    }

                    _ => {
                        panic!("invalid code. {}", code);
                    }
                }
            });

        (aligned_target, aligned_query)
    }
}

#[cfg(test)]
mod test {
    use minimap2::Aligner;

    use crate::mapping_ext::MappingExt;

    #[test]
    fn test_mapping_ext() {
        let mut aligner = Aligner::builder()
            .map_ont()
            .with_cigar()
            .with_sam_out()
            .with_sam_hit_only();

        aligner.mapopt.best_n = 1;
        aligner.mapopt.q_occ_frac = 0.0;

        // this is all map-ont default
        // aligner.mapopt.zdrop = 400;
        // aligner.mapopt.zdrop_inv = 200;
        // aligner.mapopt.bw = 500;
        // aligner.mapopt.end_bonus = 0;

        aligner.idxopt.k = 4;
        aligner.idxopt.w = 1;

        aligner.mapopt.min_cnt = 2;
        aligner.mapopt.min_dp_max = 10; // min dp score
        aligner.mapopt.min_chain_score = 10; // this is important for short insert
        aligner.mapopt.min_ksw_len = 0;

        /*
        AAATCCCCCCCG-------AAACGGGGGGGGGTTTTAAACAA
        AAATCCCCCCCGTTTTTTTAAACGGGGGGGGGTTTTAAACAA
        */
        let refseq = b"AAATCCCCCCCGAAACGGGGGGGGGTTTTAAACAA";
        let query_seq = b"AAATCCCCCCCGTTTTTTTAAACGGGGGGGGGTTTTAAACAA";
        let aligner = aligner.with_seq_and_id(refseq, b"ref").unwrap();

        for hit in aligner
            .map(
                query_seq,
                false,
                false,
                None,
                Some(&[67108864]),
                Some(b"query"),
            )
            .unwrap()
        {
            let hit_ext = MappingExt(&hit);
            let identity1 = hit_ext.identity_without_long_indel(5);
            assert!((identity1 - 1.0).abs() < 1e-3, "identity1:{}", identity1);
            let identity2 = hit_ext.identity();
            assert!((identity2 - 35.0 / 42.0) < 1e-3, "identity2:{}", identity2);

            let identity3 = hit_ext.identity_gap_compressed();
            assert!((identity3 - 35.0 / 36.0) < 1e-3, "identity3:{}", identity3);

            println!("{:?}", hit);

            println!("{:?}", hit_ext.aligned_2_str(refseq, query_seq));

            break;
        }

        /*
        AAATCCCCCCCGAAACGGGGGGGGGTTTTAAACAA
        AAATCCCCCCCGAAAC---------TTTTAAACAA
        */
        let query_seq = b"AAATCCCCCCCGTTTTTTTAAACTTTTAAACAA";
        for hit in aligner
            .map(
                query_seq,
                false,
                false,
                None,
                Some(&[67108864]),
                Some(b"query"),
            )
            .unwrap()
        {
            let hit_ext = MappingExt(&hit);
            let identity1 = hit_ext.identity_without_long_indel(5);
            assert!((identity1 - 1.0).abs() < 1e-3, "identity1:{}", identity1);
            let identity2 = hit_ext.identity();
            assert!((identity2 - 26.0 / 35.0) < 1e-3, "identity2:{}", identity2);

            let identity3 = hit_ext.identity_gap_compressed();
            assert!((identity3 - 26.0 / 27.0) < 1e-3, "identity3:{}", identity3);
            println!("{:?}", hit);
            break;
        }

        let query_seq = b"TTGTTTAAAACCCCCCCCCGTTTAAAAAAACGGGGGGGATTT";
        for hit in aligner
            .map(
                query_seq,
                false,
                false,
                None,
                Some(&[67108864]),
                Some(b"query"),
            )
            .unwrap()
        {
            let hit_ext = MappingExt(&hit);
            println!("{:?}", hit_ext.aligned_2_str(refseq, query_seq));

            break;
        }
    }
}
