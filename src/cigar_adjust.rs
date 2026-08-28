use crate::{mapping_ext::MappingExt, minimap2};

const GAP: u8 = '-' as u8;

/// 左移 target 上的 gap（即 query 相对 target 的插入），返回是否发生了移动。
/// 前提：被跳过的 target 碱基与 gap 前 query 两列碱基相同（同聚物/等价区）。
/// 交换前后两条串投影掉 gap 后不变，是等价变换；gap 只会左移，保证终止。
fn shift(aligned_target: &mut [u8], aligned_seq: &[u8], cur_idx: usize) -> bool {
    let mut shifted = false;
    for idx in cur_idx..(aligned_target.len() - 1) {
        let target_cur_base = aligned_target[idx];
        if target_cur_base != GAP
            && aligned_target[idx + 1] == GAP
            && target_cur_base == aligned_seq[idx]
            && target_cur_base == aligned_seq[idx + 1]
        {
            aligned_target.swap(idx, idx + 1);
            shifted = true;
        } else {
            break;
        }
    }
    shifted
}

/// 与 shift 对称：左移 query 上的 gap（即 query 相对 target 的删除），返回是否发生了移动。
fn shift_del(aligned_target: &[u8], aligned_seq: &mut [u8], cur_idx: usize) -> bool {
    let mut shifted = false;
    for idx in cur_idx..(aligned_seq.len() - 1) {
        let seq_cur_base = aligned_seq[idx];
        if seq_cur_base != GAP
            && aligned_seq[idx + 1] == GAP
            && seq_cur_base == aligned_target[idx]
            && seq_cur_base == aligned_target[idx + 1]
        {
            aligned_seq.swap(idx, idx + 1);
            shifted = true;
        } else {
            break;
        }
    }
    shifted
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

    // 插入和删除的左移会互相解锁（一侧的交换可能让另一侧在相邻列获得可移动条件），
    // 因此迭代到不动点为止。每次交换都使 gap 严格左移，总和单调递减，必然终止。
    loop {
        let mut changed = false;
        for cur_idx in (0..(aligned_target.len() - 1)).rev() {
            changed |= shift(aligned_target, &aligned_seq, cur_idx);
            changed |= shift_del(&aligned_target, aligned_seq, cur_idx);
        }
        if !changed {
            break;
        }
    }

    // println!("{:?}", hit.alignment.as_ref().unwrap().cigar);
    let new_cigars = aligned_pair_seqs_2_cigar(aligned_target, aligned_seq);
    hit.alignment.as_mut().unwrap().cigar = Some(new_cigars);
}

#[cfg(test)]
mod tests {
    use super::*;
    use minimap2::{Alignment, Mapping, Strand};

    /// 构造一个 Forward 链、覆盖 target[..tend]/query[..qend] 的 hit
    fn make_hit(cigar: Vec<(u32, u8)>, tend: usize, qend: usize) -> Mapping {
        Mapping {
            strand: Strand::Forward,
            target_start: 0,
            target_end: tend as i32,
            query_start: 0,
            query_end: qend as i32,
            alignment: Some(Alignment {
                nm: 0,
                cigar: Some(cigar),
                cigar_str: None,
                md: None,
                cs: None,
                alignment_score: None,
            }),
            ..Default::default()
        }
    }

    fn adjusted_cigar(
        cigar: Vec<(u32, u8)>,
        target_seq: &[u8],
        query_seq: &[u8],
    ) -> Vec<(u32, u8)> {
        let mut hit = make_hit(cigar, target_seq.len(), query_seq.len());
        cigar_adjust_poly_gap_left_align(&mut hit, target_seq, query_seq);
        hit.alignment.unwrap().cigar.unwrap()
    }

    #[test]
    fn test_deletion_left_align_in_homopolymer() {
        // 删除了 6 个连续 A 中的 2 个（target=T+6A, query=T+4A），
        // minimap2 倾向于把删除放在右侧（3M 2D 2M）。
        // 期望规范化为：紧跟 T 之后的两个 A 被删掉 => 1M 2D 4M
        let target = b"TAAAAAA";
        let query = b"TAAAA";
        let cigar = vec![(3, 0), (2, 2), (2, 0)]; // 3M 2D 2M

        assert_eq!(
            adjusted_cigar(cigar, target, query),
            vec![(1, 7), (2, 2), (4, 7)]
        );
    }

    #[test]
    fn test_insertion_left_align_regression() {
        // 在同聚物中插入一个 A，右对齐摆放（4M 1I 1M），应保持原有左移能力：
        // 规范化为插入紧跟 T 之后 => 1M 1I 4M
        let target = b"TAAAA";
        let query = b"TAAAAA";
        let cigar = vec![(4, 0), (1, 1), (1, 0)]; // 4M 1I 1M

        assert_eq!(
            adjusted_cigar(cigar, target, query),
            vec![(1, 7), (1, 1), (4, 7)]
        );
    }

    #[test]
    fn test_interleaved_indel_requires_fixpoint() {
        // 相邻的删除 + 插入：初始列状态 t=A A - C / s=A - A A，
        // 即 del@col1 与 ins@col2 相邻。删除左移把 col2 的 '-' 让给 query 碱基后，
        // 插入才能继续左移——单次扫描无法收敛，必须迭代到不动点。
        // 最终列状态 t=A - A C / s=- A A A => 1D 1I 1= 1X
        let target = b"AAC";
        let query = b"AAA";
        // 列操作码: =(7) D(2) I(1) X(8)
        let cigar = vec![(1, 7), (1, 2), (1, 1), (1, 8)];

        assert_eq!(
            adjusted_cigar(cigar, target, query),
            vec![(1, 2), (1, 1), (1, 7), (1, 8)]
        );
    }
}
