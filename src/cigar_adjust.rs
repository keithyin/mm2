use crate::mapping_ext::MappingExt;

const GAP: u8 = '-' as u8;

fn shift(aligned_target: &mut [u8], aligned_seq: &mut [u8], cur_idx: usize) {
    for idx in cur_idx..(aligned_target.len() - 1) {
        let target_cur_base = aligned_target[idx];
        if target_cur_base != GAP
            && aligned_target[idx + 1] == GAP
            && target_cur_base == aligned_seq[idx]
            && target_cur_base == aligned_seq[idx + 1]
        {
            aligned_target.swap(idx, idx + 1);
        } else {
            break;
        }
    }
}

fn compute_cigar_state(target_aligned_base: u8, query_aligned_base: u8) -> u8 {
    if target_aligned_base == GAP {
        return 1; // ins
    }
    if query_aligned_base == GAP {
        return 2; // del
    }

    if query_aligned_base == target_aligned_base {
        return 7; // eq
    }

    if query_aligned_base != target_aligned_base {
        return 8; // mismatch
    }

    return 0;
}

/// Vec<(u32, u8)>
fn aligned_pair_seqs_2_cigar(aligned_target: &[u8], aligned_seq: &[u8]) -> Vec<(u32, u8)> {
    let mut cigars = vec![];
    let mut cur_state = (
        1_u32,
        compute_cigar_state(aligned_target[0], aligned_seq[0]),
    );

    for idx in 1..aligned_target.len() {
        let cigar_state = compute_cigar_state(aligned_target[idx], aligned_seq[idx]);
        if cigar_state == cur_state.1 {
            cur_state.0 += 1;
        } else {
            cigars.push(cur_state);
            cur_state = (1, cigar_state);
        }
    }
    cigars.push(cur_state);
    cigars
}

pub fn cigar_adjust_poly_gap_left_align(
    hit: &mut minimap2::Mapping,
    target_seq: &[u8],
    query_seq: &[u8],
) {
    // hit.alignment.unwrap().cigar
    let hit_ext = MappingExt(hit);

    let (mut aligned_target, mut aligned_seq) = hit_ext.aligned_2_str(target_seq, query_seq);

    let aligned_target = unsafe { aligned_target.as_bytes_mut() };
    let aligned_seq = unsafe { aligned_seq.as_bytes_mut() };

    for cur_idx in (0..(aligned_target.len() - 1)).rev() {
        shift(aligned_target, aligned_seq, cur_idx);
    }

    // println!("{:?}", hit.alignment.as_ref().unwrap().cigar);
    let new_cigars = aligned_pair_seqs_2_cigar(aligned_target, aligned_seq);
    hit.alignment.as_mut().unwrap().cigar = Some(new_cigars);
}

#[cfg(test)]
mod test {
    use minimap2::Aligner;

    use crate::{
        cigar_adjust::cigar_adjust_poly_gap_left_align, mapping_ext::MappingExt, visualization,
    };

    #[test]
    fn test_cigar_adjust_poly_gap_left_align_1() {
        let target_seq = "GGGGGGGGGGGGGGGGGCCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATTCCCCCGGTCCACCTGGTCTCTCCACCTTAGGTACAAAGACGGTTCGGCACTGCCGTAGAATTTCGGGTATTTCGCCTCGTGCCATCCATGCATTGAACATTTCCGCCTTCAAGTGCACAGGAACCGCACGCCACTGACCCGAACGTATACCGTCCGGGCCCGGCGAAGTTCGCCAGTCAAAGCGGGAGGCCTTGATCTCTTCCACCGAGATCGGCTTCCACAGCTGGGTGTAGTCGCGATTAGACCCAATTCCTGCTATCCCCAAGCACCGCGTCCAGTACTGCGGAATTCTCCAGCGTACCGTCGGCGCAGATAAATCCGCGCTGTCGTGCATCAGGGGGGCAGCAAGCCAACAGCCTCCGGG";
        let query_seq = "CCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATTCCCCCCGGGTCCACCTGGTCTCTCCACCTTAGGTACAAAGACGGTCGGCACTGCCGTAGAATTTCGGGTATTTCGCCTCGTGCCATCCATGCATTGAACATTTCCGCCTTCAAGTGCACAGGAACCGCACGCCACTGACCCGAACGTATACCGTCCGGGCCCGGCGAAGTTCGCCAGTCAAAGCGGGAGGCCTTGATCTCTTCCACCGAGATCGGCTTCCACAGCTGGGTGTAGTCGCGATTAGACCCAATTCCTGCTATCCCCAAGCACCGCGTCCAGTACTGCGGAATTCTCCAGCGTACCGTCGGGCGCAGATAAATCCGCGCTGTCGTGCATCAGGGGGGCAGCAAGCCAACAGCCTCGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

        let aligner = Aligner::builder()
            .map_ont()
            .with_cigar()
            .with_sam_out()
            .with_sam_hit_only();

        let aligner = aligner
            .with_seq_and_id(target_seq.as_bytes(), b"ref")
            .unwrap();

        for mut hit in aligner
            .map(
                query_seq.as_bytes(),
                false,
                false,
                None,
                Some(&[67108864, 68719476736]),
                Some(b"q"),
            )
            .unwrap()
        {
            visualization::blast_like_alignment(query_seq, target_seq, &hit, None);

            cigar_adjust_poly_gap_left_align(&mut hit, target_seq.as_bytes(), query_seq.as_bytes());

            let hit_ext = MappingExt(&hit);

            let (aligned_target, aligned_seq) =
                hit_ext.aligned_2_str(target_seq.as_bytes(), query_seq.as_bytes());

            assert_eq!(&aligned_target[..100], "CCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATT-CCCCC-GGTCCACCTGGTCTCTCCACCTTAGGTACAAAGACG");
            assert_eq!(&aligned_seq[..100],    "CCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATTCCCCCCGGGTCCACCTGGTCTCTCCACCTTAGGTACAAAGACG");

            visualization::blast_like_alignment(query_seq, target_seq, &hit, None);
        }
    }

    #[test]
    fn test_cigar_adjust_poly_gap_left_align_2() {
        let target_seq = "GGGGGGGGGGGGGGGGGCCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATTCCCCCGGTCCACCTGGTCTCTCCACCTTAGGTACAAAGACGGTTCGGCACTGCCGTAGAATTTCGGGTATTTCGCCTCGTGCCATCCATGCATTGAACATTTCCGCCTTCAAGTGCACAGGAACCGCACGCCACTGACCCGAACGTATACCGTCCGGGCCCGGCGAAGTTCGCCAGTCAAAGCGGGAGGCCTTGATCTCTTCCACCGAGATCGGCTTCCACAGCTGGGTGTAGTCGCGATTAGACCCAATTCCTGCTATCCCCAAGCACCGCGTCCAGTACTGCGGAATTCTCCAGCGTACCGTCGGCGCAGATAAATCCGCGCTGTCGTGCATCAGGGGGGCAGCAAGCCAACAGCCTCCGGG";
        let query_seq = "CCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATTCCCCCCCCGGGGTCCACCTGGTCTCTCCACCTTAGGTACAAAGACGGTCGGCACTGCCGTAGAATTTCGGGTATTTCGCCTCGTGCCATCCATGCATTGAACATTTCCGCCTTCAAGTGCACAGGAACCGCACGCCACTGACCCGAACGTATACCGTCCGGGCCCGGCGAAGTTCGCCAGTCAAAGCGGGAGGCCTTGATCTCTTCCACCGAGATCGGCTTCCACAGCTGGGTGTAGTCGCGATTAGACCCAATTCCTGCTATCCCCAAGCACCGCGTCCAGTACTGCGGAATTCTCCAGCGTACCGTCGGGCGCAGATAAATCCGCGCTGTCGTGCATCAGGGGGGCAGCAAGCCAACAGCCTCGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

        let aligner = Aligner::builder()
            .map_ont()
            .with_cigar()
            .with_sam_out()
            .with_sam_hit_only();

        let aligner = aligner
            .with_seq_and_id(target_seq.as_bytes(), b"ref")
            .unwrap();

        for mut hit in aligner
            .map(
                query_seq.as_bytes(),
                false,
                false,
                None,
                Some(&[67108864, 68719476736]),
                Some(b"q"),
            )
            .unwrap()
        {
            visualization::blast_like_alignment(query_seq, target_seq, &hit, None);

            cigar_adjust_poly_gap_left_align(&mut hit, target_seq.as_bytes(), query_seq.as_bytes());

            let hit_ext = MappingExt(&hit);

            let (aligned_target, aligned_seq) =
                hit_ext.aligned_2_str(target_seq.as_bytes(), query_seq.as_bytes());

            assert_eq!(&aligned_target[..100], "CCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATT---CCCCC--GGTCCACCTGGTCTCTCCACCTTAGGTACAAAG");
            assert_eq!(&aligned_seq[..100],    "CCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATTCCCCCCCCGGGGTCCACCTGGTCTCTCCACCTTAGGTACAAAG");

            visualization::blast_like_alignment(query_seq, target_seq, &hit, None);
        }
    }

    #[test]
    fn test_cigar_adjust_poly_gap_left_align_3() {
        let target_seq = "ACAGATCTTTCCCTCGGGCATTCTCAAGACGGATCCCCATTTCCAGACGATAAGGCTGCATTAAATCGAGCGGGCGGAGTACGCCATACAAGCCGGAAGCATTCGCAAATGCTGTTGGCAAAATCGAAATCGTCTTCGCTGAAGGTTTCGGCCTGCAAGCCGGTGTAGACATCACCTTTAAACGCCAGAATCGCCTGGCGGGCATTCGCCGGCGTGAAATCTGGCTGCCAGTCATGAAAGCGAGCGGCGTTGATACCCGCCAGTTTGTCGCTGATGCGCATCAGCGTGCTAATCTGCGGAGGCGTCAGTTTCCGCGCCTCATGGATCAACTGCTGGGAATTGTCTAACAGCTCCGGCAGCGTATAGCGCGTGGTGGTCAACGGGCTTTGGTAATCAAGCGTTTTCGCAGGTGAAATAAGAATCAGCATATCCAGTCCTTGCAGGAAATTTATGCCGACTTTAGCAAAAAATGAGAATGAGTTGATCGATAGTTGTGATTACTCCTGCGAAACATCATCCCACGCGTCCGGAGAAAGCTGGCGACCGATATCCGGATAACGCAATGGATCAAACACCGGGCGCACGCCGAGTTTACGCTGGCGTAGATAATCACTGGCAATGGTATGAACCACAGGCGAGAGCAGTAAAATGGCGGTCAAATTGGTAATAGCCATGCAGGCCATTATGATATCTGCCAGTTGCCACATCAGCGGAAGGCTTAGCAAGGTGCCGCCGATGACCGTTGCGAAGGTGCAGATCCGCAAACACCAGATCGCTTTAGGGTTGTTCAGGCGTAAAAAGAAGAGATTGTTTTCGGCATAAATGTAGTTGGCAACGATGGAGCTGAAGGCAAACAGAATAACCCACAAGGGTAACAAACTCAGCACCCCAGGAACCCATTAGCACCCGCATCGCCTTCTGGATAAGCTGAATACCTTCCAGCGGCATGTAGGTTGTGCCGTTACCCGCCAGTAATATCAGCATGGCGCTTGCCGTACAGATGACCAGGGTGTCGATAAAAATGCCAATCATCTGGACAATCCCCTTGCGCTGCCGGATGCGGAGGCCAGGACGCCGCTGCCGCTGCCGCGTTTGGCGTCGAACCCATTCCCGCCTCATTGGAAAACATACTGCGCTGAAAACCGTTAGTAATCGCCTGGCTTAAGGTATATCCGCCGCGCCGCCTGCCGCTTCCTGCCAGCCAAAAGCACTCTCAAAAATAGACCAAATGACGTGGGGAAGTTGCCCGATATTCATTACGCAAATTACCAGGCTGGTCAGTACCCAGATTATCGCCATCAACGGGACAAAGCCCTGCATGAGCCGGGCGACGCCATGAAGACGCGAGTGATTGCCAGCAGAGTAAAGACAGCGAGAATAATGCCTGTCACCAGCGGGGGAAAATCAAAAGAAAAACTCAGGGCGCGGGCAACGGCGTTCGCTTGAACTCCGCTGAAAATTATGCCATAGGCGATGAGCAAAAAGACGGCGAACAGAACGCCCATCCAGCGCATCCCCAGCCCGCGCGCCATATACCATGCCGGTCCGCCACGAAACTGCCCATTGACGTCACGTTCTTTATAAAGTTGTGCCAGAGAACATTCGGCAAACGAGGTCGCCATGCCGATAAACGCGGCAACCCACATCCAAAAGACGGCTCCAGGTCCACCGCGCGGTAATAGCCAGCGCAACGCCGGCCAGGTTGCCGCTACCCACGCGCGCCGCAAGACTGGTACACAATGACTGAAATGAGGTTAAACCGCCTGGCTGTGGATGAATGCTATTTTTAAGACTTTTGCCAAAC";
        let query_seq = "GTTTGGCAAAAGTCTTAAAAATAGCATTCATCCACACCCAGCGGTTTACCCATTCAGTCATTGTGTACCAGTCTTGCGGCGCGGCTTGGGTAGCGGCAACCTGGCCGGCGTTGGCGCGGCTATTACCGCGGTGGCCTGTGAGCCGTCTTTTGATGTGGGATTGCCGGCGTTTATCGCATGGCGACCTCGTTTGCCAAGTTCTCTGGCACAACTTTTATAAAGAACGTGACGTCCAATGGGCAGTTTCGTGGGCGGACCGCATGGTATATGGCGCGCGGGCTGGGGGATGCGCTGGATGGGCGTTCATGTTCGCCGTCTTTTTGCTCATCGCCTATGGCATAATTTTCAGCGAGTTCAGCGAACGCCGTTGCCGCGCCCTGAGTTTCTTTTGATTTTCCCCGCTGGTGACAGCATTATTCTCGATGTCTTTACTCTGCTGAATCACATCGCGTCTTCATGTGGCGTCGCCCGAGTCCTCATGCAGGGCTTTGTCCCGTTGATGGAGCTAATCTGGGTACTGACCAGATCTGGTAATTTGCGTATAGAATATGGTGCAACTTCCCCACGCCATTTGGTCTATTTTTTGGATGTTTTTTGGGCTGGCAGGAAGCGGCAGGCGGCGCGGCGGGATATACCTTAAAGCAGGCGATTACTAACGGTTTTTTTCAGCCAGTATGTTTTTCCAATGAAGCGGGAGGGTTTCGAGCCAAACGCGCAGCGGCGTGGGGGGGGGCGGCTCCCTCGCCTCCGCATCCGCAGCGCAAAGGATTGTCCAGAAGATTGCATTTTATCGACACCCTGGTCATCTGTACAGGAAGCGCCGCTGATATTACGCGGTAACGGCACAACCTACAGCGCTGAGGTATTCAGCTTATCCAGAAGGCGATGCGGGTGCTAATGGGTTCCGTGTGGTGCTGAGTTTGTTACCCTTGTGGTTATCTTTTGCCTTCACCCATCGTGCCAACTTACATTGATGCGAAAAAATTCTTCTTTTTTACGCCTGAAACAACCTAAAGCGATCTGGTGTTTGCAAAGGATCTGCACCTTTCGCAACAGGTCATCGCGGCACCTTGCTAAGCCTTCCGCGATGTGGCAACTGGCAGATATCATAATGGCTGCATGGCTATACCGATTTGACGCCATTTTACTGCTCTCGCCTTGTGGTTCTACCATTGCCAAGTGATTATCTACAGCCAGCGTAAACTCGGCGTGCGCCGGTGTTTTGATCATTGCGTTATCCCGGAATCGGTCGCCCAGCTTTCTCCGGACGCGTGGGATGATGTTTCGCAGGTAATCACAATCAACTCCATTCTCATTTTTGCCTAAAGTCCGGCCATAATTTTCCTGCAAGGACTGGATATGCTGATTTTAATTTCACCTGCGAAACGCTTGATTACCAAAAGCCTTGACCATCCACGCCGTTATACGCTGCCGAGCTGTTAGAAATTCGCAGCAGTTGATCCCAGAGGTCCGGAAACTGACGCCTCCCGCAGATTAGCACGTTTTTGGGATGCGCTAGCGACAACTGGCGGGGTATCACGCCGCTCGCTTTCATTGGACTGGCAGCCAGATTTCACGCCGGCGATGCCCGCCGGCGATTCTGGCGTTAAAGGTGATGTCTACCACCAGGCTTGCAGCCAAACCTTCAGGAAGACGATTTCGATTTTGCCCAACAGCATTTGCGAATGCTTCCGCTTGGTAGGCTTACTCCGCCCTCCGATTTAATGCCTTTATCGTCCTGGAAATGGGATCCGGTCTTGAGAATGCCCGAGGGAAAGATCTGT";

        let mut aligner = Aligner::builder()
            .map_ont()
            .with_cigar()
            .with_sam_out()
            .with_sam_hit_only();
        // aligner.mapopt.a =
        aligner.idxopt.set_hpc();
        aligner.mapopt.a = 4;
        aligner.mapopt.b = 8;
        aligner.mapopt.q = 4;
        aligner.mapopt.q2 = 48;
        aligner.mapopt.e = 2;
        aligner.mapopt.e2 = 1;

        let aligner = aligner
            .with_seq_and_id(target_seq.as_bytes(), b"ref")
            .unwrap();

        for mut hit in aligner
            .map(
                query_seq.as_bytes(),
                false,
                false,
                None,
                Some(&[67108864, 68719476736]),
                Some(b"q"),
            )
            .unwrap()
        {
            visualization::blast_like_alignment(query_seq, target_seq, &hit, None);

            // cigar_adjust_poly_gap_left_align(&mut hit, target_seq.as_bytes(), query_seq.as_bytes());

            // let hit_ext = MappingExt(&hit);

            // let (aligned_target, aligned_seq) =
            //     hit_ext.aligned_2_str(target_seq.as_bytes(), query_seq.as_bytes());

            // assert_eq!(&aligned_target[..100], "CCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATT---CCCCC--GGTCCACCTGGTCTCTCCACCTTAGGTACAAAG");
            // assert_eq!(&aligned_seq[..100],    "CCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATTCCCCCCCCGGGGTCCACCTGGTCTCTCCACCTTAGGTACAAAG");

            // visualization::blast_like_alignment(query_seq, target_seq, &hit, None);
        }
    }
}
