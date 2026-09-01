use std::collections::HashMap;
use std::sync::Arc;

use hp_tr_finder::intervaltree::IntervalTree;
use hp_tr_finder::{single_seq_hp_tr_finder, UnitAndRepeats};
use lazy_static::lazy_static;
use minimap2::{Mapping, Strand};
use regex::Regex;

use crate::cigar_adjust::{aligned_pair_seqs_2_cigar, GAP};
use crate::mapping_ext::MappingExt;

lazy_static! {
    /// HP/TR 检测集: unit 长度 1~4、重复数 >=3。与 gsmm2-metric 的 `HP_TR_REG` 以及
    /// `bin/gsmm2-hp-tr-metric.rs` 使用的检测集一致
    pub static ref HP_TR_REG: Vec<HashMap<String, Regex>> = {
        vec![
            UnitAndRepeats::new(1, 3).build_finder_regrex(),
            UnitAndRepeats::new(2, 3).build_finder_regrex(),
            UnitAndRepeats::new(3, 3).build_finder_regrex(),
            UnitAndRepeats::new(4, 3).build_finder_regrex(),
        ]
    };
}

/// 某个 target 上的 motif 区间集合。区间为半开的 `[start, end)`(`query_point` 的覆盖判断
/// 即 `start <= pos < end`), value 是 motif 串, 如 `(AC)5`, 只用于日志
pub type HpTrRegions = IntervalTree<usize, Arc<String>>;

/// 扫描 target 序列, 建 HP/TR motif 索引。
///
/// 代价是 `HP_TR_REG` 里约 340 条正则各扫一遍 target, 属于建 index 时的一次性开销,
/// 由 `build_aligner` 的 per-target 线程并行分摊
pub fn build_hp_tr_regions(target_seq: &str) -> Arc<HpTrRegions> {
    // motif 串驻留表, 只在单个 target 内共享(并行构建时避免共享可变)
    let mut match_patterns: HashMap<String, Arc<String>> = HashMap::new();
    let region2motif = single_seq_hp_tr_finder(&HP_TR_REG, &mut match_patterns, target_seq);

    // 用 interval tree 而不是 gsmm2-metric 的 flatten() HashMap: 后者每个被覆盖的碱基都要
    // 存一条 (位置 -> motifs) 记录, 人基因组尺度内存扛不住, 而这里只需要点查询
    Arc::new(region2motif.to_interval_search_tree())
}

/// 覆盖 `pos` 的 motif 中最靠右的右端点
fn covering_max_end(regions: &HpTrRegions, pos: usize) -> Option<usize> {
    regions.query_point(pos).map(|e| e.range.end).max()
}

/// 覆盖 `pos` 的 motif 中最靠左的左端点
fn covering_min_start(regions: &HpTrRegions, pos: usize) -> Option<usize> {
    regions.query_point(pos).map(|e| e.range.start).min()
}

/// 收紧比对区间: 左端点若落在 HP/TR 内部就推到覆盖它的 motif 右端之外, 右端点同理往左拉到
/// motif 左端之外。要求 `end > start`。
///
/// 返回 `None` 表示收紧后区间塌空(整段比对都在重复区内部)。
///
/// 与 gsmm2-metric 的差别: metric 是单步收缩(`hp_tr_metric_v2.rs:117-131`), 收缩后的端点
/// 仍可能落在另一个 motif 里(motif 会嵌套, 如 `ACAACAAC[A]AAA` 上一个位置同时属于 `(ACA)3`
/// 和 `(A)4`; 也会首尾相接), 这里迭代到不动点。
pub fn shrink_interval(regions: &HpTrRegions, start: usize, end: usize) -> Option<(usize, usize)> {
    let mut s = start;
    // 覆盖 s 的区间必有 end > s, 故 s 严格递增, 迭代必然终止
    while let Some(motif_end) = covering_max_end(regions, s) {
        if motif_end <= s {
            break;
        }
        s = motif_end;
        if s >= end {
            return None;
        }
    }

    let mut e = end;
    // 覆盖 e-1 的区间必有 start < e, 故 e 严格递减, 迭代必然终止
    while e > s {
        let motif_start = match covering_min_start(regions, e - 1) {
            Some(motif_start) if motif_start < e => motif_start,
            _ => break,
        };
        e = motif_start;
    }

    if e <= s {
        return None;
    }

    Some((s, e))
}

/// 把比对末端落在 HP/TR 内部的 hit 收紧: 按 motif 表收缩 target 区间, 裁掉收缩区间外的比对列,
/// 重写 cigar 与 target/query 坐标。被裁掉的 read 碱基由
/// `crate::convert_mapping_cigar_to_record_cigar` 依据新的 query_start/query_end 补成
/// soft-clip, 不需要在这里显式写进 cigar。
///
/// 不重算 alignment_score / mapq / nm: 裁剪后的解与原 DP 解不再严格对应, 属于 clipping 的
/// 常规处理。cs / md 会按新的比对列重新生成, 否则它们会描述已经被裁掉的区域。
///
/// 返回 `Some(true)` 表示发生了改动, `Some(false)` 表示无需改动, `None` 表示收紧会让比对
/// 塌空或退化到没有 ref/query 跨度, 只能保持原样
fn shrink_hit_core(
    hit: &mut Mapping,
    target_seq: &[u8],
    query_seq: &[u8],
    regions: &HpTrRegions,
) -> Option<bool> {
    let (t_start, t_end) = (hit.target_start as usize, hit.target_end as usize);
    if t_end <= t_start {
        return Some(false);
    }

    let (new_start, new_end) = shrink_interval(regions, t_start, t_end)?;
    if (new_start, new_end) == (t_start, t_end) {
        return Some(false);
    }
    // 只往内收
    debug_assert!(new_start >= t_start && new_end <= t_end);

    let is_rev = matches!(hit.strand, Strand::Reverse);
    let (aligned_target, aligned_query) = MappingExt(hit).aligned_2_str(target_seq, query_seq);

    // 同 cigar_adjust: minimap2 生成 eqx 时比的是 2-bit 编码的碱基, 与大小写无关(软遮蔽的
    // reference 往往是小写的), 这里统一大写, 否则相同的碱基会被当成 mismatch
    let at: Vec<u8> = aligned_target
        .as_bytes()
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let aq: Vec<u8> = aligned_query
        .as_bytes()
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    if at.is_empty() || at.len() != aq.len() {
        return Some(false);
    }

    // 每列记录它开始时的 target 游标: 有 ref 碱基的列是该碱基的位置, 插入列是其后第一个 ref
    // 位点(BAM 里插入报在该位点之前)。游标单调不减, 所以保留的列必然是一个连续段
    let mut tpos = Vec::with_capacity(at.len());
    let mut qpos = Vec::with_capacity(at.len());
    let (mut tcur, mut qcur) = (t_start, 0usize);
    for idx in 0..at.len() {
        tpos.push(tcur);
        if at[idx] == GAP {
            // ins: 只有 query 碱基
            qpos.push(Some(qcur));
            qcur += 1;
        } else if aq[idx] == GAP {
            // del: 只有 ref 碱基
            qpos.push(None);
            tcur += 1;
        } else {
            qpos.push(Some(qcur));
            qcur += 1;
            tcur += 1;
        }
    }

    let c0 = match tpos.iter().position(|&p| p >= new_start && p < new_end) {
        Some(c0) => c0,
        None => return None, // 保留区间里没有可比对的列, 等价于塌空
    };
    let c1 = match tpos.iter().rposition(|&p| p >= new_start && p < new_end) {
        Some(c1) => c1 + 1,
        None => return None,
    };

    // 新的 target 边界: 保留列中第一个/最后一个消耗 ref 碱基的列
    let mut ref_span = (new_start, new_start);
    let mut j_lo = usize::MAX;
    let mut j_hi = 0usize;
    for idx in c0..c1 {
        if at[idx] != GAP {
            ref_span.0 = ref_span.0.min(tpos[idx]);
            ref_span.1 = ref_span.1.max(tpos[idx] + 1);
        }
        if let Some(j) = qpos[idx] {
            j_lo = j_lo.min(j);
            j_hi = j_hi.max(j + 1);
        }
    }
    // 全是插入列(0 ref 跨度)或全是删除列(0 query 跨度)都属于退化, 不裁剪
    if ref_span.1 <= ref_span.0 || j_lo >= j_hi {
        return None;
    }

    let (query_start, query_end) = if !is_rev {
        // aligned_query 的下标 j 就是 read[qstart..qend] 内的偏移
        (
            hit.query_start as usize + j_lo,
            hit.query_start as usize + j_hi,
        )
    } else {
        // reverse 时 cigar 走的是 revcomp(read[qstart..qend]), 下标 j 对应原始坐标
        // query_end - 1 - j, 故保留段 [j_lo, j_hi) 对应 [query_end - j_hi, query_end - j_lo)
        (hit.query_end as usize - j_hi, hit.query_end as usize - j_lo)
    };
    debug_assert!(query_start >= hit.query_start as usize);
    debug_assert!(query_end <= hit.query_end as usize);
    debug_assert!(query_start < query_end);

    // 只在原 hit 确实带了 cs / md 时才重写, 不凭空给 record 加 tag
    let (has_cs, has_md) = match hit.alignment.as_ref() {
        Some(aln) => (aln.cs.is_some(), aln.md.is_some()),
        // 没有 alignment 的 hit 本来就在 map_query_to_target 里被过滤掉了
        None => return Some(false),
    };
    let (cs, md) = rebuild_cs_md(&at[c0..c1], &aq[c0..c1]);

    let aln = hit.alignment.as_mut().unwrap();
    aln.cigar = Some(aligned_pair_seqs_2_cigar(&at[c0..c1], &aq[c0..c1]));
    // 本工具不消费 cigar_str, 但留着会与裁剪后的 cigar 矛盾
    aln.cigar_str = None;
    if has_cs {
        aln.cs = Some(cs);
    }
    if has_md {
        aln.md = Some(md);
    }

    hit.target_start = ref_span.0 as i32;
    hit.target_end = ref_span.1 as i32;
    hit.query_start = query_start as i32;
    hit.query_end = query_end as i32;

    Some(true)
}

/// `shrink_hit_core` 的对外入口: 塌空/退化的情况保持原样并记一条 warn
pub fn shrink_hit(hit: &mut Mapping, target_seq: &[u8], query_seq: &[u8], regions: &HpTrRegions) {
    let (t_start, t_end) = (hit.target_start, hit.target_end);
    let (q_start, q_end) = (hit.query_start, hit.query_end);
    // 日志要在可变借用结束后才能用, 这里先拷出来
    let qname = hit
        .query_name
        .as_ref()
        .map(|v| v.to_string())
        .unwrap_or_default();
    let target_name = hit
        .target_name
        .as_ref()
        .map(|v| v.to_string())
        .unwrap_or_default();
    match shrink_hit_core(hit, target_seq, query_seq, regions) {
        Some(true) => tracing::debug!(
            "hp-tr-shrink: qname:{}, target:{}, target region {}-{} -> {}-{}, query region {}-{} -> {}-{}",
            qname,
            target_name,
            t_start,
            t_end,
            hit.target_start,
            hit.target_end,
            q_start,
            q_end,
            hit.query_start,
            hit.query_end
        ),
        Some(false) => {}
        // 整段比对都在重复区内部, 收紧会塌成 0 长。0 长的对齐会让 reference_end == pos,
        // 排序/建索引都会出问题, 保持原样
        None => tracing::warn!(
            "hp-tr-shrink: whole alignment inside hp/tr region, skipped. qname:{}, target:{}, target region {}-{}",
            qname,
            target_name,
            t_start,
            t_end
        ),
    }
}

/// A/C/G/T(不区分大小写)原样大写返回, 其余(N、ONT 的 '*' 等)归一成 N。
/// minimap2 是把碱基转成 2-bit 编码后比较的, 非 ACGT 一律编成 4(=N), 这里保持一致
fn norm_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' => b.to_ascii_uppercase(),
        _ => b'N',
    }
}

/// 由裁剪后的比对列重建 cs / md, 与 minimap2 `format.c` 的
/// `write_cs_ds_core`(no_iden=1, 即 `mm_gen_cs` 最后一个参数传 true) 和 `write_MD_core`
/// 逐字一致:
///
/// - 两条列串都在 reference forward 方向(`MappingExt::aligned_2_str` 已对 reverse 的 query
///   片段做过 reverse complement), 与 minimap2 取 tseq/qseq 的方式相同
/// - 匹配段写 `:n`(不是 `=seq`), 每个匹配段结束就 flush, 所以 indel 前后各一个 `:` token
/// - substitution 逐碱基一个 `*ab`, 相邻 substitution 不会合并
/// - insertion 写成 `+seq`、deletion 写成 `-seq`, 碱基一律小写
/// - md 在每次 mismatch 与每次 deletion 之前都打印累计的匹配距离(相邻 mismatch 之间因此会
///   出现 `0`), ref 碱基大写, insertion 不产生任何内容
fn rebuild_cs_md(aligned_target: &[u8], aligned_query: &[u8]) -> (String, String) {
    assert_eq!(aligned_target.len(), aligned_query.len());

    let mut cs = String::new();
    let mut md = String::new();

    let mut match_run = 0usize; // cs 里待 flush 的匹配长度
    let mut md_dist = 0usize; // md 里待 flush 的匹配距离
    let mut idx = 0usize;
    while idx < aligned_target.len() {
        let t = aligned_target[idx];
        let q = aligned_query[idx];

        if t != GAP && q == GAP {
            // deletion: 连续的 del 列合成一个 `-` token
            let mut del = String::new();
            while idx < aligned_target.len()
                && aligned_target[idx] != GAP
                && aligned_query[idx] == GAP
            {
                del.push(norm_base(aligned_target[idx]) as char);
                idx += 1;
            }
            if match_run > 0 {
                cs.push_str(&format!(":{}", match_run));
                match_run = 0;
            }
            cs.push_str(&format!("-{}", del.to_lowercase()));
            md.push_str(&format!("{}^{}", md_dist, del));
            md_dist = 0;
            continue;
        }

        if t == GAP {
            // insertion: 连续的 ins 列合成一个 `+` token, 不影响 md
            let mut ins = String::new();
            while idx < aligned_target.len() && aligned_target[idx] == GAP {
                ins.push(norm_base(aligned_query[idx]) as char);
                idx += 1;
            }
            if match_run > 0 {
                cs.push_str(&format!(":{}", match_run));
                match_run = 0;
            }
            cs.push_str(&format!("+{}", ins.to_lowercase()));
            continue;
        }

        let (t, q) = (norm_base(t), norm_base(q));
        if t == q {
            match_run += 1;
            md_dist += 1;
        } else {
            if match_run > 0 {
                cs.push_str(&format!(":{}", match_run));
                match_run = 0;
            }
            cs.push_str(&format!(
                "*{}{}",
                (t as char).to_lowercase(),
                (q as char).to_lowercase()
            ));
            md.push_str(&format!("{}{}", md_dist, t as char));
            md_dist = 0;
        }
        idx += 1;
    }

    if match_run > 0 {
        cs.push_str(&format!(":{}", match_run));
    }
    if md_dist > 0 {
        md.push_str(&format!("{}", md_dist));
    }

    (cs, md)
}

#[cfg(test)]
mod test {
    use std::sync::Arc;

    use minimap2::{Alignment, Mapping, Strand};

    use super::{build_hp_tr_regions, rebuild_cs_md, shrink_hit, shrink_interval, HpTrRegions};

    /// 测试用 target, 84bp:
    /// `[0,20)` 非重复前缀, `[20,32)` = (AC)6, `[32,52)` 非重复中段, `[52,64)` = 12A 同聚物,
    /// `[64,84)` 非重复后缀。前后缀不能用 `"CGAT".repeat(n)`——那本身就是一条 (CGAT)n 串联重复
    fn target_seq() -> String {
        [
            "AGGTCCTAGATCCGATGACA",
            "ACACACACACAC",
            "GTTCAGGCTCAAGCTTGACG",
            "AAAAAAAAAAAA",
            "TTCGGATCAGGTCGATCGAA",
        ]
        .concat()
    }

    /// (AC)6 在 target 上的半开区间
    const STR: (usize, usize) = (20, 32);
    /// 12A 同聚物在 target 上的半开区间
    const HP: (usize, usize) = (52, 64);

    fn regions() -> Arc<HpTrRegions> {
        build_hp_tr_regions(&target_seq())
    }

    /// 手工构造一个只带 motif `[s, e)` 的索引, 便于精确验证收缩逻辑
    fn regions_of(motifs: &[(usize, usize)]) -> HpTrRegions {
        motifs
            .iter()
            .map(|(s, e)| (*s..*e, Arc::new(format!("motif{}-{}", s, e))))
            .collect()
    }

    /// motif 覆盖到的位置集合, 用于校验 fixture
    fn covered_positions(regions: &HpTrRegions, len: usize) -> Vec<usize> {
        (0..len)
            .filter(|&pos| regions.query_point(pos).next().is_some())
            .collect()
    }

    #[test]
    fn test_build_hp_tr_regions_detects_motifs() {
        // fixture 里只有 (AC)6 与 12A 两处会被判定成 HP/TR, 前后缀必须干净
        assert_eq!(
            covered_positions(&regions(), target_seq().len()),
            (STR.0..STR.1).chain(HP.0..HP.1).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_shrink_interval_single_motif() {
        let regions = regions_of(&[(100, 200)]);

        // 左端点落在 motif 内部 -> 推到 motif 右端
        assert_eq!(shrink_interval(&regions, 150, 300), Some((200, 300)));
        // 右端点(rend-1 是最后一个碱基)落在 motif 内部 -> 拉到 motif 左端
        assert_eq!(shrink_interval(&regions, 50, 150), Some((50, 100)));
        // 两端都落在 motif 内部, 整个 motif 被掏空
        assert_eq!(shrink_interval(&regions, 120, 180), None);
        // 整段在 motif 内 -> 塌空
        assert_eq!(shrink_interval(&regions, 120, 130), None);
        // motif 完整落在比对区间内部, 两端干净 -> 不动
        assert_eq!(shrink_interval(&regions, 50, 300), Some((50, 300)));
        // 端点正好压在 motif 边界上不算被覆盖(半开区间)
        assert_eq!(shrink_interval(&regions, 200, 300), Some((200, 300)));
        assert_eq!(shrink_interval(&regions, 50, 100), Some((50, 100)));
    }

    #[test]
    fn test_shrink_interval_nested_motif() {
        // 一个位置同时属于 (ACA)3 与 (A)4 那类嵌套 motif: 收缩必须取最靠外的边界
        let regions = regions_of(&[(100, 120), (110, 130)]);
        assert_eq!(shrink_interval(&regions, 105, 200), Some((130, 200)));
        assert_eq!(shrink_interval(&regions, 50, 125), Some((50, 100)));
    }

    #[test]
    fn test_shrink_interval_fixpoint() {
        // 首尾相接的两段 motif: 单步收缩会停在 [100,120) 的右端, 而该点仍被 [120,140) 覆盖
        let regions = regions_of(&[(100, 120), (120, 140)]);
        assert_eq!(shrink_interval(&regions, 105, 200), Some((140, 200)));

        let regions = regions_of(&[(60, 80), (80, 100)]);
        assert_eq!(shrink_interval(&regions, 50, 90), Some((50, 60)));
    }

    /// 构造一个 Forward 链的 hit。target/query 坐标都是 minimap2 的约定: query 坐标在
    /// **原始 read 方向**上, reverse 链的 cigar 走的是 revcomp(read[qstart..qend])
    fn make_hit(
        cigar: Vec<(u32, u8)>,
        target_start: i32,
        target_end: i32,
        query_start: i32,
        query_end: i32,
    ) -> Mapping {
        Mapping {
            strand: Strand::Forward,
            target_start,
            target_end,
            query_start,
            query_end,
            alignment: Some(Alignment {
                nm: 0,
                cigar: Some(cigar),
                cigar_str: None,
                md: Some(String::new()),
                cs: Some(String::new()),
                alignment_score: Some(1),
            }),
            ..Default::default()
        }
    }

    fn revcomp(seq: &str) -> String {
        let v: Vec<u8> = seq
            .bytes()
            .rev()
            .map(|b| match b.to_ascii_uppercase() {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                b => b,
            })
            .collect();
        String::from_utf8(v).unwrap()
    }

    fn aln_of(hit: &Mapping) -> &minimap2::Alignment {
        hit.alignment.as_ref().unwrap()
    }

    #[test]
    fn test_shrink_hit_tail_in_str() {
        // read 与 [0,26) 完全匹配, 末端停在 (AC)6 内部
        let target = target_seq();
        let read = target[..26].to_string();

        let mut hit = make_hit(vec![(26, 7)], 0, 26, 0, 26);
        shrink_hit(&mut hit, target.as_bytes(), read.as_bytes(), &regions());

        assert_eq!((hit.target_start, hit.target_end), (0, STR.0 as i32));
        assert_eq!((hit.query_start, hit.query_end), (0, 20));
        assert_eq!(aln_of(&hit).cigar, Some(vec![(20, 7)]));
        assert_eq!(aln_of(&hit).cs.as_deref(), Some(":20"));
        assert_eq!(aln_of(&hit).md.as_deref(), Some("20"));
    }

    #[test]
    fn test_shrink_hit_head_in_hp() {
        // read 从 12A 同聚物([52,64))中间的第 56 位开始, 之后是干净的非重复后缀
        let target = target_seq();
        let read = target[56..76].to_string();
        assert_eq!(target[56..64], "A".repeat(8));

        let mut hit = make_hit(vec![(20, 7)], 56, 76, 0, 20);
        shrink_hit(&mut hit, target.as_bytes(), read.as_bytes(), &regions());

        // 左端点被推到同聚物右端之外
        assert_eq!((hit.target_start, hit.target_end), (HP.1 as i32, 76));
        // 被裁掉的是 read 最前面的 8 个 A(它们会变成 soft-clip)
        assert_eq!((hit.query_start, hit.query_end), (8, 20));
        assert_eq!(aln_of(&hit).cigar, Some(vec![(12, 7)]));
    }

    #[test]
    fn test_shrink_hit_untouched() {
        // 两端都干净, (AC)6 与 12A 完整落在比对区间中间: 不改动
        let target = target_seq();
        let read = target.clone();

        let mut hit = make_hit(vec![(84, 7)], 0, 84, 0, 84);
        shrink_hit(&mut hit, target.as_bytes(), read.as_bytes(), &regions());

        assert_eq!((hit.target_start, hit.target_end), (0, 84));
        assert_eq!((hit.query_start, hit.query_end), (0, 84));
        assert_eq!(aln_of(&hit).cigar, Some(vec![(84, 7)]));
        // 没被裁过的 hit 不碰 cs/md
        assert_eq!(aln_of(&hit).cs.as_deref(), Some(""));
        assert_eq!(aln_of(&hit).md.as_deref(), Some(""));
    }

    #[test]
    fn test_shrink_hit_reverse() {
        // reverse 链: cigar 走 revcomp(read[qstart..qend]), 校验 query 坐标换回原始方向的结果
        let target = target_seq();
        let read = revcomp(&target[10..26]);

        let mut hit = make_hit(vec![(16, 7)], 10, 26, 0, 16);
        hit.strand = Strand::Reverse;
        shrink_hit(&mut hit, target.as_bytes(), read.as_bytes(), &regions());

        assert_eq!((hit.target_start, hit.target_end), (10, STR.0 as i32));
        // 落在重复区里的是 read 的 5' 端(原始坐标 [0,6)), 故 query_start 被推高
        assert_eq!((hit.query_start, hit.query_end), (6, 16));
        assert_eq!(aln_of(&hit).cigar, Some(vec![(10, 7)]));
        assert_eq!(aln_of(&hit).cs.as_deref(), Some(":10"));
    }

    #[test]
    fn test_shrink_hit_whole_alignment_in_repeat() {
        // 整段比对都在 (AC)6 内部 -> 收紧会塌空, 保持原样
        let target = target_seq();
        let read = target[22..28].to_string();

        let mut hit = make_hit(vec![(6, 7)], 22, 28, 0, 6);
        shrink_hit(&mut hit, target.as_bytes(), read.as_bytes(), &regions());

        assert_eq!((hit.target_start, hit.target_end), (22, 28));
        assert_eq!((hit.query_start, hit.query_end), (0, 6));
        assert_eq!(aln_of(&hit).cigar, Some(vec![(6, 7)]));
    }

    #[test]
    fn test_shrink_hit_keeps_indels_and_rebuilds_cs_md() {
        // target[10..23) = TCC GATG ACA ACA
        // read 带一个 2bp 插入(在 10+3 之后)、一个 1bp 删除(target[15])、一个错配
        // (target[16]=G -> read C), 末端停在 (AC)6 里
        let target = target_seq();
        let read = format!(
            "{}GG{}C{}",
            &target[10..13],
            &target[13..15],
            &target[17..23]
        );
        assert_eq!(read.as_bytes().len(), 14);

        let mut hit = make_hit(
            vec![(3, 7), (2, 1), (2, 7), (1, 2), (1, 8), (6, 7)],
            10,
            23,
            0,
            14,
        );
        shrink_hit(&mut hit, target.as_bytes(), read.as_bytes(), &regions());

        assert_eq!((hit.target_start, hit.target_end), (10, 20));
        assert_eq!((hit.query_start, hit.query_end), (0, 11));
        assert_eq!(
            aln_of(&hit).cigar,
            Some(vec![(3, 7), (2, 1), (2, 7), (1, 2), (1, 8), (3, 7)])
        );
        assert_eq!(aln_of(&hit).cs.as_deref(), Some(":3+gg:2-t*gc:3"));
        assert_eq!(aln_of(&hit).md.as_deref(), Some("5^T0G3"));
    }

    #[test]
    fn test_rebuild_cs_md_matches_minimap2() {
        // 期望值取自 minimap2 format.c 的 write_cs_ds_core(no_iden=1) / write_MD_core,
        // 并用 `gsmm2 align` 在真实数据上核对过形态(`:20*gt*cg:42` / `20G0C42` 这类)
        let (cs, md) = rebuild_cs_md(b"ACGT", b"ACGT");
        assert_eq!((cs.as_str(), md.as_str()), (":4", "4"));

        // 单碱基错配: 匹配段先 flush 成 `:n`, 再写 `*<ref><query>`
        let (cs, md) = rebuild_cs_md(b"ACGTACGT", b"ACGTACCT");
        assert_eq!((cs.as_str(), md.as_str()), (":6*gc:1", "6G1"));

        // 相邻错配不合并, md 每次错配前都打距离, 于是出现 0; 首尾各有 1 个匹配位点
        let (cs, md) = rebuild_cs_md(b"ACGT", b"AATT");
        assert_eq!((cs.as_str(), md.as_str()), (":1*ca*gt:1", "1C0G1"));

        // 插入: `+ttt` 单独成 token, md 完全不受影响(7 个 ref 碱基全匹配)
        let t = b"ACGCA---TG";
        let q = b"ACGCATTTTG";
        let (cs, md) = rebuild_cs_md(t, q);
        assert_eq!((cs.as_str(), md.as_str()), (":5+ttt:2", "7"));

        // 删除: cs `-cgg`, md `^CGG`, 距离照常跨过去累计
        let t = b"ACCGGTA";
        let q = b"AC---TA";
        let (cs, md) = rebuild_cs_md(t, q);
        assert_eq!((cs.as_str(), md.as_str()), (":2-cgg:2", "2^CGG2"));

        // 开头就是删除: 没有匹配段可 flush, cs 直接以 `-` 开头; md 无条件打距离 -> `0^`
        let (cs, md) = rebuild_cs_md(b"ACGT", b"--GT");
        assert_eq!((cs.as_str(), md.as_str()), ("-ac:2", "0^AC2"));

        // 非 ACGT 一律按 N 处理(2-bit 编码), N vs N 算匹配
        let (cs, md) = rebuild_cs_md(b"ACNT", b"ACNT");
        assert_eq!((cs.as_str(), md.as_str()), (":4", "4"));
        // ONT 数据里 read 的 '*' 也归一成 N
        let (cs, md) = rebuild_cs_md(b"ACGT", b"AC*T");
        assert_eq!((cs.as_str(), md.as_str()), (":2*gn:1", "2G1"));

        // 软遮蔽的小写 reference 不影响判定(minimap2 比的是 2-bit 编码)
        let (cs, md) = rebuild_cs_md(b"acgt", b"ACGT");
        assert_eq!((cs.as_str(), md.as_str()), (":4", "4"));
    }
}
