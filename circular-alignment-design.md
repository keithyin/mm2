# gsmm2 circular（环状参考）比对方案设计

状态：**方案已定（v4，2026-09-23 三次 review 后修订），未动代码**。
修订要点与依据见 §8、§9。v4 修掉的是三处会让实现卡住的硬伤（fold 的入参/出参载体、promote 会踩 `src/lib.rs:832`
的断言、§4.3.5 与 §4.5 的 noMar 口径冲突），并补上 minimap2 内部阈值的影响与测试的独立 oracle。

## 1. 目标

给 gsmm2 增加 circular 参考序列的比对能力。基本实现方式：

- 比对阶段：用 `seq + seq`（2L）建 index，照常比对（chaining / DP / 全部后处理都在 2L 坐标系上完成）。
- 输出阶段：结果折回原始 `seq`（L）的坐标系；跨 origin（seam）的比对拆成两段分别输出；复制带来的"孪生 hit"在这里被去掉。
- 正交性：现有所有功能（query-forward、polyNGapLeftAlign、hpTrShrinkAlnRegion、noSeco/noSupp/noMar、oupCovT/oupIyT 等）与 circular 组合均可用，且 circular 关闭时行为与现状完全一致（`fold_hits` 为空操作、过滤口径不变）。

为什么要复制序列：minimap2 2.28（`minimap2-sys-0.1.21+minimap2.2.28/minimap2/minimap.h`）没有任何 circular 语义
（`MM_F_*` / `MM_I_*` 里最接近的只有 `MM_I_HPC`/`MM_I_NO_SEQ`/`MM_I_NO_NAME`），复制序列 + 折回坐标是唯一可行路线，
也是社区对环状参考的通行做法（minimap2 官方没有 circular 模式，本文只把它当背景、不依赖任何外部结论）。
§8 的实验证明它确实给出了 linear 参考给不出的正确比对。

## 2. 已定决策

| 决策点 | 结论 |
|---|---|
| primary / supplementary 归属 | 跨 seam 拆出的两段中，**折回后 pos=1 的那段（B 段）继承父 hit 的 primary/supp 状态**（保留原 score / mapq），**尾部段（A 段）为 supplementary**。两个 record 打配对 tag：`XT:Z:CIRC` + `XG:i:<group>`（group 在 fold 阶段分配，见 §4.3.2），供下游把两段拼回。 |
| 两段的 alignment_score | **两段都沿用父 hit 的总分**，不重算 DP 分。理由：该值只参与 `set_primary_supp_alignment_according_2_align_score` 的排序，按列重算的"局部 DP 分"既没有 minimap2 end-bonus/zdrop 的语义，又要重新引入一套打分参数 |
| 折回规则 | **只有一条规则：target 坐标对 L 取模，跨 seam 的区间在 seam 处按列切一刀**（§4.3.1）。**不存在"丢弃 part2 的 hit"这种独立规则**——"保留 part1 / 丢弃 part2"与"把 part2 折回 part1 开头"是同一个操作，写成两条会互相打架 |
| 孪生 hit | **必须在 fold 阶段去重**（§4.3.2）。key = `(strand, target_name, 折回后 ts, 折回后 te, 折回后 qs, 折回后 qe)`，**必须含 query 区间**——`read span > L` 时同一条 hit 折回后的 A / B 两段会落在同一个 target 区间而 query 区间不同，只按 target 区间去重会静默吃掉一段真比对（§4.3.2、§4.3.5）。不去重会输出重复 record，并直接破坏 noMar / oupCovT |
| `read span > L`（绕圈） | 折回后可能产出**同 target 区间、不同 query 区间的多段，且这些段都是真的**（不是孪生、不能去重）；`s > 2L - c` 时根本没有单条表示，minimap2 会给出多条 query 区间不相交的 chain，**它们也全是真 supplementary**（§4.3.5、§4.7） |
| noMar / oupCovT / oupIyT 过滤时机 | 口径前移到 **hit 层面**：noMar 按**去重后的 group 数**判断（CIRC 对算一条），oupCovT / oupIyT 按**父 hit** 的 coverage / identity 判断，两段共享同一个值（§4.5）。执行顺序维持现状不变，非 circular 路径行为不变 |
| CLI 开关形式 | 全局 `--circular`（bool），所有 target 都当 circular。混合 FASTA（线性 + 环状）暂不支持 |
| mapq | 保持 minimap2 原值（**已知退化**：只要 read 在 2L 上存在第二条等价的 chain —— 不跨 origin 的 read 全都在此列，`s > L` 的 read 也会 —— mapq 就掉到 0，见 §4.6）。**这是本方案影响面最大的未决点**（§7 另有若干待定小项） |
| **数据载体（v4）** | `Mapping` 是 minimap2 crate 的类型（`minimap2-0.1.23/src/lib.rs:198`），字段固定，装不下 fold 需要的东西。**`TargetCtx` 加 `circular_len: Option<usize>`**（`None` = 非 circular → `fold_hits` 空操作；L 不能用 `target_seq.len() / 2` 反推，变体 B 下会错），**hits 载体升级成 `Vec<(Mapping, TargetCtx, HitMeta)>`**（`HitMeta { group: u32, circular: bool, parent_iy: f32, parent_ec: f32 }`），`align_single_query_to_targets*` 额外返回 `group_cnt` 供 noMar 使用（§4.3） |
| **近似重复（v4）** | 去重以精确坐标为主。**平移 / 截断出来的等价段来自另一个 DP 窗口，坐标可能差 1–2bp**，精确 key 会漏；漏掉的那条 record 会变成一个多余的 group，被 `--noMar` 误杀（§4.3.2、§4.5）。fold 里再加一条**带容差的包含式去重**（δ ≤ 3bp），且不得误并 `x500_c100` 那种"同 target 区间、不同 query 区间"的段 |

## 3. 现状管线（相关部分）

`src/bin/gsmm2.rs::alignment()`：

1. 读 target fasta → `Vec<ReadInfo>`。
2. 若 `--query-forward`，`targets_to_query_forward_targets`（`src/lib.rs:761`）展开成 `${name}___fwd` / `${name}___rev`。
3. `targets_to_targetsidx`（`src/lib.rs:742`）建 `name → (tid, len)`，用于 BAM header 的 SQ 记录（SN/LN）与 record.tid。
4. `build_aligner`（`src/lib.rs:124`）为每条展开后的 target 建一个 minimap2 index；`NoMemLeakAligner` 同时持有该 target 的 `target_seq`（供 cigar 后处理读碱基）与 `hp_tr_regions`（在同一序列上扫描 motif）。
5. `align_worker`（`src/lib.rs:306`）：`map_query_to_target`（map + DP）→ `hp_tr_shrink` → `poly_gap_left_align`（均通过 `target_seq[target_start..target_end]` 按 hit 坐标读 target 碱基）→ `set_primary_supp_alignment_according_2_align_score`（基于 query range overlap）→ `records.len() > 1` 的 noMar 过滤（`src/lib.rs:322`）。
6. `hits2records` → `build_bam_record_from_mapping`（`src/lib.rs:544`）：strand 决定 query 是否 RC；`convert_mapping_cigar_to_record_cigar` 加首尾 softclip；`pos = hit.target_start`；cs / md 原样透传。
7. `write_bam_worker`（`src/bam_writer.rs:57`）按 `targetsidx` 写 SQ header 后写 record；**写出前用 `OupParams::valid`（`src/params.rs:306`，由 `src/bam_writer.rs:105` 调用）做 noSeco / noSupp / oupIyT / oupCovT 过滤，并现算 `iy` / `ec` 两个 aux tag**（`src/bam_writer.rs:109-112`）。注意 `BamOupArgs`（`src/bam_writer.rs:13`）是没有任何调用点的死代码。

关键点：**index 建完之后，下游全部处于"建 index 所用序列"的坐标系**。circular 方案就是利用这一点，在入口处换序列、出口处折回坐标。

## 4. 方案：入口展开，出口折回

与现有 `query_forward` 模式同构：circular 是一个"target 展开 + 出口折回坐标"的步骤，中间段（index 构建、chaining、DP、全部后处理）保持在 2L 坐标系不动。

### 4.1 流程

```
targets ──[1. circular: seq → seq+seq]──[2. query-forward: fwd/rev 展开]──▶ build_aligner(2L 序列, 不变)
                                                                                     │
align_worker: map ──▶ hp_tr_shrink ──▶ poly_gap_left_align          （2L 坐标系, 全部不变）
            ──▶ [3. fold_hits: 孪生去重 + 分配 group + 折回 L / 跨 seam 拆段]   （新增, circular 关闭时为空操作）
            ──▶ set_primary_supp_alignment_according_2_align_score（对拆段后的 hits 跑）
            ──▶ hits2records（每段独立建 record；softclip 按段自己的 [qs, qe) 加；
                             iy / ec 用父 hit 的值；打 XT / XG 配对 tag）
            ──▶ noMar 过滤（按 fold 阶段分配好的 group 计数）
            ──▶ write_bam_worker: noSeco / noSupp + oupIyT / oupCovT（用 record 上的 iy / ec）
```

过滤的**执行顺序与现状完全一致**（noMar 在 align_worker、iy/ec 在 writer），改的只是**口径**（hit/group 层面而非 record 层面）——见 §4.5。

### 4.2 circular 展开

`targets_to_circular_targets(targets)`（新函数，样式照 `targets_to_query_forward_targets`）：

- 每条 target 的 `seq` 换成 `seq + seq`，**名字不变**（不需要后缀：同一 target 只有一份，且名字要原样出现在输出 BAM 里）。
- 执行顺序：circular 展开**在 query-forward 展开之前**。
- `targets_to_targetsidx` 在 circular 模式下记录**原始长度 L**（展开前的 `seq.len()`），BAM header 的 LN 因此仍是 L，`write_bam_worker` 不需要改（LN 本来就来自 targetsidx）。

与 query-forward 组合时每个 target 展开为两条链（执行顺序无关，因为 RC(S+S) = RC(S)+RC(S)）：

- `t___fwd`：序列 `S+S`，长 2L
- `t___rev`：序列 `RC(S+S)`，长 2L

**代价**：index 序列、`NoMemLeakAligner.target_seq`（`src/lib.rs:54`，索引序列的又一份拷贝）与 HP/TR motif 扫描
（~340 条正则）全部按 2L 走 → **内存与建索引时间约 2×**。变体 B 说的"省一半"指的就是这部分。

> **变体 B（仅当 §4.6 决定要救 mapq 时才采用）**：只复制侧翼，`S + S[0..K]`（K ≥ 最大 read 的 target span）。
> 孪生只出现在距 origin K 以内的 read（占比 ~K/L），染色体尺度下从 ~100% 降到 ~1%，还省一半 index 内存与构建时间；
> 代价是多一个 K 参数（CLI 给保守默认值 + 超限告警，或先扫一遍 query 文件取 max length），
> 折回从"两段"变成"三段"（`p < K` → `p + L - K`；`K ≤ p < K+L` → `p - K`；`p ≥ K+L` → `p - K - L`），
> 而且小 target（质粒 / 6kb）下 K ≥ L 会退化成 `S+S`、收益为零。默认采用 `S+S`（简单、无参数、无静默截断风险）。

### 4.3 fold_hits

在 `align_single_query_to_targets` / `align_single_query_to_targets_query_forward` 的 hp_tr_shrink、poly_gap_left_align **之后**、`set_primary_supp_alignment_according_2_align_score` **之前**调用。circular 关闭时直接原样返回。

**插入位置有硬约束**：折回来要重建 cigar / cs / md，必须读 2L 序列（列串重建），所以它得放在
`(Mapping, TargetCtx)` **还成对**的那一层（和 `shrink_hp_tr_aln_region_hits` / `poly_gap_left_align_hits` 并列），
不能等 `src/lib.rs:397` / `src/lib.rs:486` 把 ctx 丢掉之后再跑；**也必须在 query-forward 的 `pk_target_strands`
（`src/lib.rs:466`）之后**——pk 用的是 2L 坐标系上 minimap2 给的 primary（§8.1 实测这个 primary 常常是 part2 那份孪生），
fold 提前会把 pk 的输入改掉。签名变成 `&mut Vec<(Mapping, TargetCtx, HitMeta)>`（要 push 新 hit，也要把
group / 父 hit 的 `iy`、`ec` 带出去）。

**L 与开关的来源**：`TargetCtx`（`src/lib.rs:62`）现在只有 `target_seq` + `hp_tr_regions`，L 和开关两样都没有。给它加
`circular_len: Option<usize>`，由 `build_aligner`（`src/lib.rs:124`，`NoMemLeakAligner` 要同步加这个字段）在建 index 时写入：
`Some(L)` = 该 target 按 circular 处理，`None` = 空操作。**不要用 `target_seq.len() / 2` 反推**，变体 B 下是错的。

**出参**：group 数要传到 `align_worker`（noMar）、`HitMeta` 要传到 `hits2records` / `build_bam_record_from_mapping`
（`iy` / `ec` / `XT` / `XG`），所以 `align_single_query_to_targets*` 的返回从 `Vec<Mapping>` 变成
`(Vec<Mapping>, usize)`（或等价的 struct）；`hits2records` / `build_bam_record_from_mapping` 多带一个 `&HitMeta`
（`src/lib.rs:945/1046/1359/1421` 四个测试调用点要跟着改）。

#### 4.3.1 折回规则（target 坐标，L 为原始长度）

**只有一条规则**：把 target 坐标对 L 取模；跨 seam 的区间在 seam 处按列切一刀。下面三种形态是同一条规则的三个分支：

对每个 hit（`target_start < target_end`）：

1. **`target_end ≤ L`**（整段在 part1）：原样保留（`p mod L = p`）。
2. **`target_start ≥ L`**（整段在 part2）：两端各减 L。
3. **跨 seam：`target_start < L < target_end`**：拆成两段，见 §4.3.3。

> **不要把它写成"丢弃 part2 的 hit"**（形态 2 的相反规则）。两者在"孪生必然存在"时结果相同，但那个前提不是设计可以依赖的不变量：
> - **孪生的 primary / secondary 归属不确定**：实测 minimap2 常常把 **part2 那份**标成 primary（§8.1 的 `nc100_c50` / `plain`），按位置丢 part2 会丢掉 primary；
> - "孪生必然存在"依赖 minimap2 的 chaining 与 `pri_ratio` / `best_n` 剪枝，全是启发式；折回 + 去重是无条件的，不需要这个假设；
> - 位置规则无法表达形态 3 的 A / B 切分，两条规则并存时必然互相矛盾。

query-forward 组合时不需要镜像换算：fwd 链 `S+S[p mod L] == S[p mod L]`；rev 链同理 `RC(S+S)[p mod L] == RC(S)[p mod L]`，
两条链同构、seam 都在 target 坐标 L（`RC(S+S) = RC(S)+RC(S)`，物理上就是同一个接缝），
折回后 record 落在 `${name}___rev` 的同一 0-based 位置、reverse strand。
**注意这只对 target 坐标成立；query 侧的切点需要换算，见 §4.3.3。**

#### 4.3.2 孪生去重（新增，必须做）

复制序列后，**任何不跨 origin 的 read 在 2L 上都有两个完全等价的位置**（圈上坐标 c 处，以及 c+L 处，只要 c+span ≤ L）。
minimap2 会对其中一个标 primary、另一个标 secondary，两条折回后 **(target 区间, query 区间) 完全重合** → 不去重就会输出一条多余的重复 record。
跨 origin 的 read 还会多出一条被 2L 末端截断的复制品，折回后同样与父 hit 的某一段重合（末端截断的那条对上 A 段，
从另一端平移出来的那条对上 B 段——实测 `x100_c190` 报出的是后者，`10S90=`@0 与 B 段 `[0,90)`/`q[10,100)` 完全重合）。

去重规则（在折回之后、primary/supp 判定之前）：

- key = `(strand, target_name, 折回后的 target_start, 折回后的 target_end, 折回后的 query_start, 折回后的 query_end)`；
  **query 区间这一维是必须的**：`read span > L` 时同一条 hit 折回后的 A / B 两段可以落在同一个 target 区间而
  query 区间不同（§4.3.5 的 `x500_c100`），只按 target 区间去重会静默吃掉一段真比对。query 用**原始 query 方向**
  的值（与 `set_primary_supp_alignment_according_2_align_score` 同一口径，reverse 链由 §4.3.4 的共享 helper 换算）；
- 保留优先级：**`is_primary` > `is_supplementary` > query span 大 > score 高**。
  - 不要用"总是保留第一份拷贝"——实验里 minimap2 恰恰可能把第二份（part2 那份）标成 primary（§8.1）；
  - `is_supplementary` 排在 score 前面，是为了不依赖"补充段的父 hit 分数一定更高"这个未必成立的前提：
    同 key 的两段必然来自不同的父 hit，而 supplementary 表示该段携带独立的 query 区间，比分数更能说明谁该留下；
- 去重后若整个 hits 里一条 primary 都不剩（被丢掉的那条恰好是唯一的 primary），把保留者提升为 primary：
  `set_primary_supp_alignment_according_2_align_score` 只会降级 primary、不会提升，且 `src/lib.rs:831` 有
  `assert!(hits.first().unwrap().is_primary)`。
  **promote 必须同时置 `is_primary = true` 和 `is_supplementary = false`**：只置前者的话，`src/lib.rs:832` 的
  `assert!(!hits.first_mut().unwrap().is_supplementary)` 会直接 panic。可达路径：父 hit 是 primary，跨 seam 拆出的
  B 段因退化被丢（跨 seam 的一整段 D），只剩带 supplementary 位的 A 段被提升；
- **同时分配 group id**：一次 fold 里拆出的两段共享父 hit 的 group；被去重丢掉的孪生，其 group 合并到保留者。
  group 是 noMar / oupCovT 的口径基础（§4.5），也写成 record 上的 `XG:i:<group>`。
- **包含式去重（v4 新增，兜 §4.3.5 的近似重复）**：上面的精确 key 覆盖不到"平移 / 截断出来的拷贝与父 hit 的某一段差
  1–2bp"的情形。这类段的 query 区间必然是父段的前缀或后缀、target 区间几乎相同，所以再补一条包含式规则：按优先级
  排序后，若已保留的 Y 与 X 同 strand、同 target，且 X 的 query 区间（±δ）被 Y 包含、target 区间两端差 ≤ δ（δ = 3 够用），
  就丢弃 X（group 并入 Y）。**这条不得误并 `x500_c100` 那种"同 target 区间、不同 query 区间"的段**——它们的 query
  区间互不包含，天然安全。δ 是"依赖 minimap2 剪枝"与"多一条低分 record 打破 noMar / oupCovT 口径"之间的取舍。

#### 4.3.3 跨 seam 拆段

对 `target_start < L < target_end` 的 hit，在 2L 坐标系里按**列**切一刀，得到两条 hit，再把 B 段整体减 L 折回：

**切列规则**（精确）：取列串（§4.3.4），按每列起始 ref 坐标建游标 `tpos`（插入列的 `tpos` = 其后第一个 ref 位点，
沿用 `src/hp_tr_shrink.rs:139-159` 的约定），然后

```
c_A  = 最小的 i, 使得 at[i] != GAP 且 tpos[i] >= L      // 第一条吃掉 seam 处/之后 ref 碱基的列
A 段 = 列区间 [0, c_A)，target span [target_start, L)    // supplementary
B 段 = 列区间 [c_A, len)，target span [L, target_end)，折回后 [0, target_end - L)   // 继承父 hit 的状态
```

这条规则自动给出想要的行为：

- seam 处的**插入列**（`tpos == L` 且 `at[i] == GAP`）落在 `[0, c_A)` 里 → 归 A 段（尾部 I 比首部 I 更被下游容忍）；
- 跨 seam 的**删除**会按 ref 碱基被切开：A 段尾部留 `D(左半)`、B 段首部留 `D(右半)`。
  **不能把整段 D 归 A**：那样 A 的 cigar 会多消耗 ref 碱基，与 `target_end = L` 矛盾；
  若把 `target_end` 也设成 L+5，record 就超出 `LN = L` 的参考边界（非法 BAM，排序/索引也会出问题）；
- 切在 M/X 列中间时，两段各得一段更短的 `=`/`X`，query 与 target 同时被切开，天然自洽。

**两条可证明的不变量（直接写成 `debug_assert`，不必再写防御分支）**：case 3 下 `t_start < L` ⇒ 消耗 ref 碱基
`t_start` 的那一列坐标 < L ⇒ 它排在 `c_A` 之前 ⇒ `c_A ≥ 1` ⇒ **A 段必然含至少一个消耗 ref 的列**；B 段必含消耗
ref 碱基 `L` 的那一列 ⇒ **A / B 两段的 ref 跨度都 ≥ 1**。也就是说下面要处理的退化**只会发生在 query 侧**（整段 D）。

**首 / 尾 D 的形态（v4 补）**：跨 seam 的 D 被切开后，A 段以 `D(n)` 收尾、B 段以 `D(n)` 开头（B 段 `pos = 0`）。
按 SAM 规范 POS 是"第一个消耗 ref 碱基的 op"的位置，所以这是合法的；但多数下游默认首尾是 M/S/=/X，要按 §6.14
显式断言形态，并在文档里写明"保留、不转换成 softclip"。

**退化段**：若 A 段（或 B 段）的 query span 为 0（例如 alignment 开头就是一整段跨 seam 的 D），
不能产出 0 query span 的 record，直接只输出另一段；若被迫丢掉的是 B 段，A 段提升为 primary
（**提升时按 §4.3.2 同时把 `is_supplementary` 清掉，否则 `src/lib.rs:832` 的断言会炸**）。
（`src/hp_tr_shrink.rs:184` 已经在拒绝同类退化，复用它的判断；§4.3.3 的不变量保证退化只可能出在 query 侧。）

**两段的 flag**（`set_primary_supp_*` 之前）：

- B 段：原样继承父 hit 的 `is_primary` / `is_supplementary`；
- A 段：`is_primary = false`；`is_supplementary = 父 hit 的 is_primary || 父 hit 的 is_supplementary`
  （父 hit 是 primary 或 supplementary → A 是 supplementary；父 hit 是 secondary → A 也是 secondary）。
  这样"一个 read 不会因为拆段冒出两个 primary"。

**primary 归属在 A / B 大小悬殊时不直观**：规则固定让"折回后 pos=1 的 B 段"继承父 hit 的状态，所以当 B 段远小于
A 段时（`read span > L` 下常见，实测 `x300_c150` 的 primary 拆出的 B 段只有 50bp、read 长 300bp），该 read 的
primary record 只覆盖 read 的一小段。这与 §4.5 的"两段共享父 hit 的 ec"是自洽的，但下游若按 flag 取 primary
会拿到一个碎片。是否改成"折回后 query 覆盖最大的那段当 primary"见 §7。

**query 侧坐标必须换算**：列串的下标 j 是 **cigar 方向**的偏移，而 `Mapping.query_start/query_end`
是**原始 query 方向**的（`src/lib.rs:716` 的 `convert_mapping_cigar_to_record_cigar` 会再翻一次）。
非 query-forward 模式下、reverse 链的跨 seam hit，切点 `qm` 必须按 `query_end - j_hi` / `query_end - j_lo` 换算，
否则 softclip 与拆段位置全错。**这条由 §4.3.4 的共享 helper 负责，不要自己算。**

#### 4.3.4 实现路线：复用列区间重写，不要手写 cs / md 手术

`src/hp_tr_shrink.rs:102` 的 `shrink_hit_core` 做的正是"取 `aligned_2_str` → 建每列 `tpos`/`qpos` 游标 →
取一个列区间 `[c0, c1)` → 重建 cigar / cs / md / target 与 query 坐标"，而 circular 的拆段就是同一件事取**两个列区间**。
建议把它的尾巴（列区间 → 重写 hit）抽成共享 helper（例如 `rewrite_hit_from_columns(hit, &at, &aq, c0, c1)`），
`shrink_hit_core` 与 `fold_hits` 共用：

- `MappingExt::aligned_2_str`（`src/mapping_ext.rs:122`）拿列串（cigar 方向，reverse 时 query 已是 RC）；
- `aligned_pair_seqs_2_cigar`（`src/cigar_adjust.rs:65`，已是 `pub(crate)`）还原 cigar；
- `rebuild_cs_md`（`src/hp_tr_shrink.rs:293`，改成 `pub(crate)`）重建 cs / md。

这样 §7 里"sub-cs 的 `*`/`-` 归属"这类问题直接消失（列串里每个 ref 碱基是一列，D 自然被正确切开），
reverse 链的 query 坐标换算也已经在那里解决过（`src/hp_tr_shrink.rs:189-199`）。重建后照旧把
`aln.cigar_str = None`（留着会与裁剪后的 cigar 矛盾，`src/hp_tr_shrink.rs:215`）。CS/MD 只在父 hit 本来带 tag 时重写。

> 不要复用 `src/aligned_pairs.rs:55` 的 `TAlignedPairs::aligned_pairs`：它 reverse 时游标从
> `query_len - query_end` 起**递增**，而正确映射是 `query_end - 1 - j` 递减，等于把上面的陷阱重新引入一遍。

#### 4.3.5 `read span > L` 时折回后的段数

折回本身是逐 hit 的机械操作，但**"这一段是不是孪生"在 `s > L` 时会变得不直观**：

| 情形 | 2L 上的 hit | 折回后 |
|---|---|---|
| `s ≤ L` 且不跨 origin | 2 条（孪生对） | 各自 1 段，落点完全重合 → 去重成 1 条 |
| 跨 origin 且 `s ≤ 2L - c`（位置唯一） | 1 条 | 2 段（A / B） |
| `s > 2L - c`（无单条表示） | **≥ 2 条，query 区间不相交** | 每条各自折回，段数相加 |
| `s ≥ L`（绕圈） | 可能 1 条 hit 就覆盖整段 2L | **A / B 两段的 target 区间可以相同而 query 区间不同** |

两条必须写进实现的不变量：

1. **折回后落在同一个 target 区间、query 区间不同的多段都是真的**（read 绕了不止一圈），**不能去重**——这正是
   §4.3.2 的 key 必须含 query 区间的唯一原因。实测 L=200、`circle[100,500)`（2.5 圈）：第一段 hit
   `q[100,500) ↔ r[0,400)` 在 `r = L` 处切开后得到 A = `[0,200)` / `q[100,300)` 与 B = `[0,200)` / `q[300,500)`，
   两段 target 区间相同、query 区间不同，且都是 read 绕圈的真实经过。
2. **`s > 2L - c` 时 minimap2 给出的多条 chain 是"真 supplementary"**（query 区间不相交、彼此不重叠），不是复制
   出来的冗余，**一条都不能丢**。§4.7 只论证了可表示性，这里补上"不可表示时会退化成什么"。实测 L=200、
   `circle[150,300)`（`2L - c = 250 < 300`）得到 `50S250=` @0（primary）与 `250=50S` @150（secondary），
   折回后 4 段去重成 3 段，恰好是 read 的真实分解。

可以放心去重的"额外 supplementary"只有一种：**平移 / 截断出来的等价段**——在另一侧平移出来的拷贝（`p - L`，起点为负
被 `0` 截掉头部）对上 **B 段**（实测 L=200 的 `10S90=`@0 就是这种，与 B 段 `[0,90)`/`q[10,100)` 逐位重合）；
被 2L 末端截断的复制品（`p + L`，超出 2L 的尾巴被切掉）对上 **A 段**。两者都折回后与父 hit 的某一段重合，
§4.3.2 的 key 能命中。但注意它们是**另一个 DP 窗口**，只在理想情况下才逐位相等；read 末端带 indel / mismatch 时
两边可能差 1–2bp，精确坐标去重会漏出一条近似重复的低分 record —— **v4 起由 §4.3.2 的包含式去重兜住**（阈值 δ）。
实测这种低分链基本上先被 `pri_ratio 0.8` 剪掉（`ins_seam` 只报出一条），但那属于"依赖 minimap2 行为"，
不能当作去重的替代，要用测试钉死（§6.15）。

### 4.4 输出坐标与 header

- `record.pos` 为 0-based ∈ `[0, L)`；1-based 即 `[1, L]`，跨 origin read 的 B 段 `pos = 1`。
- header 的 SQ `LN` 用**原始长度 L**（§4.2）。
- `hit.target_name` 不变，record 的 tid 照旧从 `target_idx` 取（`src/lib.rs:607`）。
- 不变式：`cigar 的 ref 跨度 == target_end - target_start`，且两段之和 == 父 hit 的 ref 跨度；
  任何 record 的 `reference_end() <= L`（§4.3.3 的切列规则保证）。

### 4.5 过滤语义（noMar / oupCovT / oupIyT）

现状与问题：

- **noMar**：`records.len() > 1`（`src/lib.rs:322`）。拆段会把 CIRC 对变成 2 条 record 被误杀；
  而"改成按 fold-back 前的 hit 数判断"**同样不成立**——孪生也是 2 条 hit，等于把绝大多数不跨 origin 的 read 全丢掉。
  正确口径是**去重后的 group 数**：`group_cnt > 1` 才丢。非 circular 时 group 数 == hit 数 == record 数，行为不变。
  更准确地说：**group 数是 linear 路径下 record 数的对应物**——每条真实比对（孪生合并后）算一条。两者的差别只出现在
  跨 origin 的 read 上：linear 会把它算成 2 条（两端各自 chaining）而丢掉，circular 算 1 条而保留，这正是本功能要的。
  实现上 group 是 fold 阶段分配的（§4.3.2），一条 hit 拆出的 A / B 两段以及被合并掉的孪生都归同一个 group。
  **已知限制（v4）**：这条对齐依赖"孪生被精确去重"。§4.3.5 的近似重复一旦漏掉一条，group 数就比 linear 的 record 数
  多 1，`--noMar` 会把整条 read 丢掉（`--oupCovT` 通常不会：近似重复的父 ec 更低，会先被它自己丢掉）。要么用
  §4.3.2 的包含式去重兜住，要么把这条写进 `--noMar` 的帮助文本。
- **oupCovT / oupIyT**：现在在 writer 里按**单条 record** 现算（`src/params.rs:306`，`ec` 来自
  `BamRecordExt::compute_query_coverage()`）。跨 origin 的 read 拆成 400 + 600 两段后，两条 record 的 ec 都 ≈ 0.4 / 0.6，
  `--oupCovT 0.9` 会把**两条一起丢掉**，整条 read 消失；写进 aux 的 `iy` / `ec` 也会是 0.4 / 0.6，误导下游。

改法（口径前移、顺序不变）：

- 在 hit 层面算 `iy = MappingExt::identity()`、`ec = MappingExt::query_coverage()`，**CIRC 的两段共享父 hit 的值**；
  这两个函数与 gskits 的 record 级 `compute_identity` / `compute_query_coverage` 是同一套公式
  （matched / Σ(=,X,I,D)、query-consuming / read 全长），无拆段时逐位相等，所以**非 circular 路径行为不变**（要加回归测试）。
- 具体落地：`build_bam_record_from_mapping` 建 record 时就把父 hit 的 `iy` / `ec` 写进 aux，
  `write_bam_worker` 改为**按 record 上的这两个值过滤**（不再现算），noSeco / noSupp 仍按 record flag 过滤。
  这样 `src/bam_writer.rs` **必须改动**（原 §5 的"预计无改动"不成立）。

组合语义（写进文档，避免误用）：

- `--noSupp` + `--circular`：跨 origin read 的 A 段被删，只剩 B 段 —— 下游看到的是一条"半截对齐"（后半段是 softclip）。
  这是"预期行为"还是应当告警，见 §7；
  **tag 会自相矛盾（v4）**：此时 record 的 `ec` 仍是父 hit 的值（= 整条 read 的覆盖），而它自己的 cigar 只覆盖一部分
  （例：`ec = 1.0` 的 record 实际只放了 40% 的 read）。下游若把 `ec` 当"这条 record 放了多少 read"就会算错——要在文档 /
  帮助文本里写明 **`iy` / `ec` 是父级（整条 read）语义，segment 级请从 cigar 自己算**，或者额外打一个 segment 级 tag；
- `--circular` + `--query-forward`：pk 用的是 2L 坐标上的 primary，strand 与位置无关、结论不变，但 fold 必须排在 pk
  之后（§4.3）；
- `--noMar` 与 `--noSeco` / `--noSupp` 在 `src/params.rs:284-292` 是断言互斥的，`--circular` 不改变这一点。

### 4.6 mapq 语义（待拍板）

**mapq 在 chaining 阶段就算好了，一张 2L 的 index 会把"次优 chain"的位置填满**：只要 read 在 2L 上存在第二条与之
等价（或部分等价、分数接近）的 chain，mapq 就掉到 0。

> **v4 更正"输出层救不回来"的依据**：minimap2 的 `mm_set_parent` 在 `mm_select_sub` **之前**跑，所以 `subsc`（mapq 用的
> 次优分数）与 `pri_ratio` 无关；而 gsmm2 本来也没暴露这些开关（`best_n = 5` 硬编码在 `src/lib.rs:164`，`pri_ratio`
> 保持默认，`-N` / `--secondary` / `-p` 无 CLI 入口）。`--noSeco` 同样不影响 minimap2 侧：`OupParams::modify_aligner`
> （`src/params.rs:331`）清的是 `MM_F_SECONDARY_SEQ (0x1000000000)`——该位在 wrapper 里从未被置位，`&=` 本来也是 no-op——
> 而不是 `MM_F_NO_PRINT_2ND (0x4000)`。结论仍然成立，但**依据是代码路径，不是 §8 的实测**（§8 只跑了默认参数）。

退化的判据**不是"是否跨 origin"**，而是"2L 上有没有第二条等价的 chain"：

| 情形 | mapq | 实测 |
|---|---|---|
| 不跨 origin（必有孪生） | 0 | L=1000 `plain` / `chim_b`、L=200 `nc100_c50` 全部 0 |
| 跨 origin 且在 2L 上位置唯一 | **保持高分** | L=200 `x150_c150` 60、`x100_c190` 60；L=1000 `ins_seam` 56、`chim_a` 59 |
| `s > L` 且存在重复 / 部分重复 chain | 0 | L=200 `x300_c150` 的 primary mapq 0 |

所以"只有跨 origin、位置唯一的 read 保持高分"成立，影响面是**几乎全部 read**；而 §2 又明确要"保留原 hit 的 mapq"。

| 选项 | 说明 | 代价 |
|---|---|---|
| **1. 接受并文档化（默认）** | `--circular` 下 mapq 语义退化，下游不得按 mapq 过滤 | 若 gsda 下游按 mapq 过滤，这是静默的严重回归 |
| **2. 只复制侧翼（§4.2 变体 B）** | 孪生只影响距 origin K 以内的 read，占比 ~K/L | 多一个 K 参数；小 target 下退化成 `S+S`，收益为零 |
| 3. 两趟比对 | 第一趟 linear 拿正确 mapq，只对疑似跨 origin 的 read 用 2L index 补一遍 | 需要启发式候选判定 + 跨趟分数合并，与"一条 read 一次 map 全部 hit 一起判定"的现有架构不合 |
| **4. 折回后重算（v4 新增）** | fold 阶段已经能区分"竞争 chain 是自己的拷贝（偏移恰好 ±L）"与"真重复（偏移 ≠ L）"，所以可以**只对前者的 hit，按去重后的 chain 集合用 minimap2 的 `mm_set_mapq` 公式重算**，把 0 抬回合理值 | 要把 `mm_set_mapq` 抄一遍（`end_bonus` / `min_chain_sc` / anchor 数都要对），是近似而非等价；判错就是给真重复的 read 发了高分 mapq |

选 1、2 还是 4 取决于 gsda 下游是否消费 mapq；**这是本方案影响面最大的未决点**（§7 另有若干小项待定）。

### 4.7 适用域与限制

设 read 的 target span 为 `s`、在圈上的起点为 `c`，则在 `S+S`（2L）中存在连续表示 **当且仅当 `s ≤ 2L - c`**：

- `s ≤ L`：任意 `c` 都可表示（**安全区**）；
- `L < s ≤ 2L`：只有一部分起点可表示（§8 里 300bp read / 200bp 环恰好落在 `c = 100`，被成功表示为单条 300=）；
- `s > 2L`：**完全无法表示**，会被拆成多条 query 区间不相交的 chain。

质粒 + ONT 长 read 很容易落在第二、三档。设计要求：

- 文档写明"安全区是 read span ≤ L"；
- **`s > L` 时不是"静默截断"**：minimap2 退化出的多条 chain 都是真比对（§4.3.5），折回后必须全部保留，
  它们之间不能互相去重，也不能按"落在 part2"丢掉；用测试固定行为（§6）。
- **target 本身是 k 折串联重复**（周期 `L/k`）：`S+S` 里有 `2k` 份等价拷贝，折回 + 去重会把它们并成 1 条。这是 mod-L
  坐标语义的必然结果（圈上它们本来就是同一个位置），但与 linear 参考下报 `2k` 条的行为不同，要写明；
- **极短 target**：`L < k`（map-ont 下 k = 15）时 `seq+seq` 仍能建 index，但没有 minimizer、不会有 hit；target 为空串
  时 `with_seq_and_id` 的 `.unwrap()`（`src/lib.rs:160`）会 panic。两条都要有明确行为（至少别静默）。

**对 minimap2 内部阈值的影响（v4 补，§4.8 也记了一行）**：`with_seq_and_id` / `with_seqs_and_ids` 会把 `mapopt.mid_occ`
无条件钉成 **1000**（`minimap2-0.1.23/src/lib.rs:866`），且它**不随 target 长度缩放**；`max_max_occ = 4095`、
`occ_dist = 500` 来自默认值（`minimap2-sys-0.1.21/minimap2/options.c:37-38`，map-ont 不覆盖它们；只有 `short` 预设设
`max_occ`，gsmm2 用不到）。于是生效的是 `mm_seed_select` 的"高频 seed 稀疏化"（`seed.c:53-92`）：**index 命中数
> 1000 的 seed 若在 query 上连续成段，要么整段被掩掉、要么按 `occ_dist = 500` 只留一小撮；命中数 > 4095 的一律掩掉**。
复制序列让每个 minimizer 的 index 命中数翻倍，所以线性下命中数 ∈ (500, 1000] 的 seed 在 circular 下会新进入这条路径。
质粒尺度（命中数 ≤ 几十）完全没有影响；大环状基因组的高重复区可能轻微掉敏感度——真要在那种 target 上跑，先做一次
linear / circular 的 hit 数对照。`q_occ_frac` 不受影响（它数的是 **query** 内部的自重复，`seed.c:5-24`）。

### 4.8 正交性说明

| 功能 | 与 circular 的交互 | 需要的处理 |
|---|---|---|
| polyNGapLeftAlign | 在 2L 坐标系上操作，同聚物跨 seam 时 2L 序列已含完整同聚物（首尾拼接），左移反而比 linear 更正确。fold-back 在左移之后，坐标自然正确 | 无 |
| hpTrShrinkAlnRegion | motif 表扫在 S+S 上，跨 seam 的 motif 天然存在；shrink **只裁比对两端、不拆段**，fold-back 在其后，顺序天然正确 | 无 |
| noSeco / noSupp | 作用在 record 层面；noSupp 会把 CIRC 对的 A 段删掉（预期行为：只剩 B 段，注意 §4.5 的组合语义） | 无 |
| oupCovT / oupIyT / noMar | 见 §4.5，口径前移到 hit / group 层面 | **writer 与 align_worker 都要改**（§5） |
| query-forward | 见 §4.2 / §4.3.1，target 折回同构、无镜像换算；**query 侧切点需要换算**（§4.3.3） | 无（换算在共享 helper 里） |
| primary/supp 判定 | `set_primary_supp_alignment_according_2_align_score` 对拆段后的 hits 运行。两段的 query range 是 **`[qs,qm)` 与 `[qm,qe)`，完全不相交（ovlp = 0）**，所以 A 段在第二阶段（`src/lib.rs:856-872`）不会被降级成 secondary；B 段（primary）与其它 target 的 primary 之间照旧按 ovlp 竞争 | 无 |
| minimap2 内部阈值（v4 新增） | 复制使每个 minimizer 的 index 命中数翻倍，而 `mid_occ` 被钉死在 1000、`max_max_occ = 4095`、`occ_dist = 500` 都不随 target 缩放（§4.7）；`best_n = 5` 的名额还会被孪生占掉一个 | **无代码改动，但要记录**：质粒尺度无影响，大环状基因组的高重复区可能轻微掉敏感度 |
| mapq | 见 §4.6 | 待拍板 |

## 5. 改动文件清单

| 文件 | 内容 |
|---|---|
| `src/params.rs` | `MapParams` 加 `circular: bool` + `set_circular`（metric 子命令用 `MapParams::default()`，天然关闭） |
| `src/lib.rs` | `targets_to_circular_targets`；**`TargetCtx` 加 `circular_len: Option<usize>`，`NoMemLeakAligner` 同步、由 `build_aligner` 写入**；hits 载体升级成 `(Mapping, TargetCtx, HitMeta)`；`align_single_query_to_targets*` 在 post-processing 之后、primary 判定之前（且 query-forward 的 pk 之后）插 `fold_hits`，并**返回 `group_cnt`**；`hits2records` / `build_bam_record_from_mapping` 带 `HitMeta`（**5 个调用点要跟着改**）；`align_worker` 的 noMar 改为按 group 数；`build_bam_record_from_mapping` 写 `iy` / `ec`（父 hit 值）与 `XT` / `XG` tag |
| 新 `src/circular.rs` | `fold_hits`（折回 + 孪生去重 + 容差去重 + group 分配 + 拆段）+ 单测 |
| `src/hp_tr_shrink.rs` | 把 `shrink_hit_core` 的"列区间 → 重写 hit"抽成共享 helper；`rebuild_cs_md` 改 `pub(crate)` |
| `src/bin/gsmm2.rs` | `--circular` arg；circular 展开调用（放在 query-forward 展开之前）；`targets_to_targetsidx` 在 circular 模式下记录原始长度 |
| `src/bam_writer.rs` | **必须改**：oupIyT / oupCovT 改为按 record 上已有的 `iy` / `ec` 过滤，不再现算（`OupParams::valid` 要读 aux；**缺 tag 时回退到现算**，别让其它调用路径 panic）；顺带删掉死代码 `BamOupArgs` |
| **明确不改**：`src/align_processor/*`、`gsmm2-aligned-metric` / `gsmm2-hp-tr-metric` / `gsmm2-time-err` | 这三个子命令自己调 `align_single_query_to_targets`（如 `src/bin/gsmm2-time-err.rs:452`），而 `AlignInfo`（`src/align_processor/mod.rs:200`）把 hit 的 target 坐标直接当参考坐标用。**本次不支持 circular**，文档里要写明"它们的坐标假设是 linear target" |

## 6. 测试计划

1. **fold_hits 单测**：手工构造四类 hit（不跨 seam / 整体在第二份 / 跨 seam / 孪生对），验证折回坐标、拆段后的 cigar /
   cs / md、切列规则（seam 落在 I 列、落在 M 列中间、落在 D 列中间）、primary 归属、group 分配。
2. **孪生去重**：read 完全落在某一份拷贝内 → **只输出 1 条 record**；重复位点在 2L 上被 minimap2 标成 primary 时
   也只能留一条（实测 `nc100_c50` / `plain` 的 **part2 那份**才是 primary，去重后 primary 必须仍在、
   §4.3.2 的"提升为 primary"分支要生效、`src/lib.rs:831-832` 的断言都不能炸）。
3. **端到端**：200bp circular target + 跨 origin 的 read（后 100bp + 前 100bp）→ 一条 CIRC 对（pos=1 的 B 段 + 尾部 A 段），
   两段列数和 == read 长度；同一 read 比 linear 参考作对照（应出现大 softclip 或低分 hit）。
4. **`read span > L`（绕圈）**——这一组是 §4.3.5 的反例测试，L=200：
   - `circle[100,500)`（2.5 圈）→ 折回后有**两段落在同一个 target 区间 `[0,200)`、query 区间分别为 `[100,300)`
     与 `[300,500)`，两段都必须输出**。把去重 key 里的 query 维去掉，这条就会漏掉后者；
   - `circle[150,300)`（`2L - c = 250 < 300`，无单条表示）→ 折回后 3 段、无重复；primary 是 50bp 的碎片
     （当前 primary 归属规则的行为，先钉死，是否改见 §7）；
   - `circle[100,400)`（正好 2 圈）→ 3 段，与 read 的真实分解逐段一致。
5. **跨 seam 带 indel**：L=1000、`circle[950,1050)` 中间插 3bp → 一条 `100=3I100=` 跨 seam 的 hit，
   拆段后 seam 处的插入列归 A 段（§4.3.3 的切列规则），A / B 两段 query 之和 == read 长度。
6. **跨 seam 带 D**：断言 `reference_end() <= LN`、`cigar ref 跨度 == target_end - target_start`、两段 ref 跨度之和 == 父 hit。
7. **reverse 链**：非 query-forward 模式下、reverse 链跨 seam 的 read，校验 softclip 与 cs / md（切点换算）。
   实测 fwd / rev 两条链的折回结果应逐位镜像（§8 的 `plain_rc` / `chim_b_rc`）。
8. **退化段**：seam 切出 0 query span 的 A 段 / B 段（§4.3.3 的不变量保证退化只可能在 query 侧——0 ref span 的 A / B
   段不可能出现，这条要断言）；以及被 2L 末端截断成 0 ref span 的"孪生"（它根本不会有 alignment，应已在
   `map_query_to_target` 的 `hit.alignment.is_some()` 过滤里被丢掉）。
9. **过滤组合**：`--circular --noMar`（跨 origin 与不跨 origin 的 read **都必须存活**）、`--circular --oupCovT 0.9`
   （跨 origin 的 read 必须存活、且两段的 `ec` tag 是父 hit 的值）、`--noSupp` / `--noSeco`、`--query-forward`。
10. **适用域**：`s > 2L - c`（无单条表示，折回后必须是多条真 hit 且互不去重）与 `s > 2L` 的行为固定。
11. **回归**：现有全部测试不回归（circular 默认关闭，`fold_hits` 为空操作；iy / ec 口径对非拆段 hit 逐位相等——
    直接比对输出 BAM 里这两个 tag 的**字节**，不是只比重算公式）。
12. **rotate oracle（v4 新增，最重要的独立校验）**：把 circular target 旋转一个随机偏移 `r` 后按 **linear** 比对
    （`--circular` 关），取那些不跨新 origin 的 read，结果坐标 `− r (mod L)` 应与 circular 输出逐位一致（正负链各测一遍）。
    手工构造 hit 只能验证想得到的形态，这条能系统性覆盖折回 + 切列 + reverse 换算；可做成 N 个随机 read × M 个随机 `r` 的 fuzz。
13. **promote 路径（v4）**：构造"父 hit 是 primary + 跨 seam 的 B 段退化成 0 query span"（跨 seam 的一整段 D），
    断言 A 段被提升为 `is_primary` **且 `is_supplementary == false`**，`src/lib.rs:831-832` 两个断言都不炸。
14. **BAM 形态（v4）**：跨 seam 的 D 切开后，断言 B 段 cigar 以 `D(n)` 开头、A 段以 `D(n)` 结尾；对整份输出跑
    `samtools quickcheck` + `sort` / `index`（覆盖首 / 尾 D 与 `reference_end() <= LN`）。
15. **tag 与 group 的一致性（v4）**：CIRC 对的两条 record `XG` 相同、`XT` 都是 `CIRC`；不跨 seam 的 record 不带这两个
    tag；同一 read 的 group 数与预期一致（例如 `x500_c100` 是 2）。

## 7. 待定的小细节

- **`read span > L` 时 primary 归给哪一段**：现在是"折回后 pos=1 的 B 段固定继承父 hit 状态"（§4.3.3），
  在 A / B 大小悬殊时会让 primary record 只覆盖 read 的一小段（实测 `x300_c150`：50bp / 300bp）。
  是否改成"折回后 query 覆盖最大的那段当 primary"。
- `--noSupp` + `--circular` 导致跨 origin read 只剩半条：只写文档，还是打 warn？
- CIRC 配对 tag 的最终键名与内容（暂定 `XT:Z:CIRC` + `XG:i:<group>`，`XG` 同时是 group 计数与配对键；
  是否需要在段上额外存 partner 的 pos 由下游需求定）。
- 变体 B（侧翼复制）的 K 从哪里来：CLI 参数 + 保守默认值（并告警超限），还是先扫一遍 query 文件取 max length。
- 变体 B 若采用，`--circular` 与 `S+S` 是否要同时保留（例如小 target 强制走 `S+S`）。
- 配对 tag 的键名要不要换：`XT` / `XG` 是本地保留 tag（SAMtags 规定 X/Y/Z 开头留给本地使用），但别的流程里
  `XT:A` 有 unique / repeat 的既有约定，落地前核一下 gsda 上下游有没有撞名。
- `--circular` + 预先建好的 `.mmi`：`--indexedTarget`（`src/bin/gsmm2.rs:101`）目前没被 `alignment()` 使用，
  将来若要支持，L 得从 index 里恢复（另存 meta，或退回 `seq.len()/2`）。
- §4.3.2 的容差去重：δ 取多少、要不要做成隐藏参数，还是干脆不做只写进 `--noMar` 的已知限制。

## 8. 附：设计依据（实验记录，2026-09-23）

用 `target/release/gsmm2`（`--threads 12`，默认 `map-ont`，等比 pbsim 风格的随机序列，读段为 target 的精确子串）
验证了方案的两个前提与两个副作用。target 为单条序列，`S+S` 由 fasta 直接拼出：

| 用例 | linear 参考 | `seq+seq` |
|---|---|---|
| L=5000、1kb read 跨 origin（`S[4600..]+S[..600]`） | `400S600=` @1 + `400=600S` @4601（被硬拆成两段各自 chaining） | **单条 `1000=` @4601，mapq 60** ✅ |
| L=5000、1kb read 全在内部（`S[2000..3000]`） | `1000=` @2001，mapq 60 | `1000=` @2001，**mapq 0** + 一条重复 record ❌（孪生；注意 minimap2 把**第二份拷贝**标成了 primary） |
| L=200、100bp read 全在内部 | `100=` @51，mapq 60 | `100=` @51 与 @251（折回后重合），mapq 均 0 ❌ |
| L=200、100bp read（`S[190..]+S[..90]`） | `10S90=` @1 | `100=` @191 + `10S90=` @1（后者折回后与 B 段完全重合）❌ |
| L=200、300bp read（1.5 圈） | `100S200=` @1 + `100=200S` @101 | **单条 `300=` @101** ✅（`s ≤ 2L - c` 的边界用例） |

推论：

1. 复制序列确实把跨 origin 的 read 从"两段各自 chaining"变成"一条正确的完整比对"，设计前提成立；
2. 不跨 origin 的 read 实测必然出现孪生 hit（→ §4.3.2 去重、§4.5 的 group 口径）。方案不依赖"必然"这个前提，
   折回 + 去重在孪生缺席时也不会出错（§4.3.1）；
3. mapq 会掉到 0，且是 chaining 内部行为（`subsc` 在 `mm_set_parent` 里定，与 `pri_ratio` 无关），输出层参数救不回来
   （→ §4.6，依据是代码路径而不是这组实验；这里只验证了默认参数下的退化）；
4. `s ≤ 2L - c` 的可表示性推导与实测一致（→ §4.7）。

### 8.1 v3 补充：边界用例（L=200 / L=1000，2026-09-23）

同一台机、同一个 binary，随机序列（seed 7 / seed 11），read 是 target 的精确子串或只带一处 indel。
命名 `x<len>_c<start>` = 长 `<len>`、圈上起点 `<start>`；`nc` = 不跨 origin；`plain` / `chim_*` / `ins_seam`
是 L=1000 的用例。POS 已转成 0-based。

**L=200（target = 2L = 400）**

| read | 2L 上的输出 | 折回后 |
|---|---|---|
| `nc100_c50` | `100=`@50 (256, mapq 0) + `100=`@250 (**primary**, mapq 0) | 两段落点重合 → 1 条 |
| `x100_c190` | `10S90=`@0 (256) + `100=`@190 (**primary**, 60) | B 段 `[0,90)` 与 @0 那条重合 → 2 条 |
| `x150_c150` | `150=`@150 (primary, 60) | A `[150,200)` / B `[0,150)` → 2 条 |
| `x300_c150` | `50S250=`@0 (primary, **mapq 0**) + `250=50S`@150 (256) | 4 段去重成 3 段；primary 是 50bp 碎片 |
| `x400_c100` | `100S300=`@0 (primary, 60) + `300=100S`@100 (256, 60) | 3 段，与 read 的真实分解逐段一致 |
| `x500_c100` | `100S400=`@0 (primary, 60) + `300=200S`@100 (256, 60) | **A = `[0,200)`/`q[100,300)` 与 B = `[0,200)`/`q[300,500)`：同 target 区间、不同 query** |

**L=1000（target = 2L = 2000）**

| read | 2L 上的输出 |
|---|---|
| `plain` = `circle[300,500)` | `200=`@300 (256, mapq 0) + `200=`@1300 (**primary**, mapq 0) |
| `chim_a` = `circle[900,1000)` + `circle[0,300)` | 单条 `400=`@900 (primary, **59**)，本身就是跨 seam 的连续表示 |
| `chim_b` = `circle[100,300)` + `circle[880,980)` | `200=100S`@100 (primary, mapq 0) + `200=100S`@1100 (256, 0)；第二段（100bp）被 `pri_ratio` 剪掉 |
| `ins_seam` = `circle[950,1050)` 中间插 3bp | 单条 `100=3I100=`@950 (primary, **56**)；折回后 A = `50=3I`、B = `150=` |
| `*_rc`（上面几条的 reverse complement） | 与 fwd 逐位镜像（pos 相同 + 负链 flag），折回坐标一致 |

这几组数据支撑的结论：

1. **孪生对的 primary 常在 part2 那份**（`nc100_c50`、`plain`）→ "按位置丢 part2"会丢错，必须折回 + 去重（§4.3.1、§4.3.2）；
2. **mapq 退化的判据是"2L 上有没有第二条等价 chain"**，不是"是否跨 origin"：`x300_c150` 的 primary 也是 0（§4.6）；
3. **`s > L` 会产生同 target 区间、不同 query 区间的多段真比对**（`x500_c100`）→ 去重 key 必须含 query（§4.3.2、§4.3.5）；
4. **`s > 2L - c` 无单条表示时退化成多条真 chain**（`x300_c150`），不是"静默截断"（§4.7）。

## 9. 修订记录

- v4（2026-09-23，三次 review 后修订；核对过 minimap2 2.28 的 C 源、gskits 0.15.12、rust-htslib 0.49）：
  - **补数据载体（§2、§4.3、§5）**：`TargetCtx` 加 `circular_len`，hits 载体升级成 `(Mapping, TargetCtx, HitMeta)`，
    `align_single_query_to_targets*` 要额外返回 group 数——原稿的 `&mut Vec<(Mapping, TargetCtx)>` 根本没法把
    group / 父 hit 的 `iy`、`ec` 带出去，而 `Mapping` 的字段是 minimap2 crate 固定的，装不下；
  - **修 promote 的断言陷阱（§4.3.2）**：提升为 primary 时必须同时清 `is_supplementary`，否则 `src/lib.rs:832` 直接
    panic；并补出可达路径（父 hit 是 primary + B 段退化成 0 query span）；
  - **补近似重复的容差去重（§2、§4.3.2、§4.5）**：精确 key 漏掉的那条会多出一个 group，`--noMar` 会误杀整条 read；
    §4.5"与 linear 判定完全对齐"的措辞随之改写并标为已知限制；
  - **补 minimap2 内部阈值的影响（§4.7、§4.8）**：`mid_occ` 被 wrapper 钉死 1000、`max_max_occ = 4095`、`occ_dist = 500`
    都不随 target 缩放，复制使命中数翻倍 → 命中数 ∈ (500, 1000] 的 seed 新进入"高频 seed 稀疏化"路径（质粒尺度无影响）；
    `q_occ_frac` 是 query 侧、不受影响；
  - **§4.3.3 补两条可证明的不变量**（A / B 的 ref 跨度都 ≥ 1，退化只可能在 query 侧）与首 / 尾 D 的形态说明；
  - **§4.6 更正"输出层救不回来"的依据**（`mm_set_parent` 在 `mm_select_sub` 之前；gsmm2 没暴露 `-N` / `-p`；
    `--noSeco` 清错了 flag），并补第 4 个选项（折回后按去重后的 chain 集合重算 mapq）；
  - **§4.3 补顺序约束**（fold 必须在 query-forward 的 pk 之后）与 `NoMemLeakAligner` 的字段同步；
  - **§4.2 补代价（≈2× 内存与建索引时间）**；§4.7 补"target 为 k 折串联重复""极短 target"两个边界；
  - **§5 补 `TargetCtx` / `NoMemLeakAligner` / `hits2records` 的签名改动，并声明 metric 子命令不在范围内**；
  - **§6 补 rotate oracle、promote 路径、BAM 形态、tag 一致性四组测试**；§7 补 tag 撞名、`.mmi`、δ 三个待定项；
  - 新增 §10：记录相邻的三处既有问题（`--noSeco` 清错 flag、fastq 后缀拼进 seq、两个未使用的 CLI 参数）。

- v3（2026-09-23，二次 review 后修订；依据见 §8.1 的实测）：
  - **折回规则统一成一条**（§4.3.1）：删掉"丢弃 part2 的 hit"这个与形态 2 相反的表述，改为"target 坐标 mod L"。
    位置过滤只有在"孪生必然存在"时才等价于去重，而实测孪生对的 primary 常是 part2 那份，且该前提依赖 minimap2
    的剪枝启发式；
  - **去重 key 补上 query 区间**（§4.3.2）：`s > L` 时一条 hit 折回后的 A / B 两段可以落在同一个 target 区间而
    query 区间不同（`x500_c100`），只按 target 区间去重会静默吃掉一段真比对；去重优先级补 `is_supplementary`；
  - 新增 §4.3.5：`read span > L` 时折回的段数、"同 target 区间的多段都是真的"、`s > 2L - c` 退化成多条真
    supplementary、以及被 2L 末端截断的复制品在带 indel 时坐标可能差 1–2bp 的容错问题；
  - 更正 §4.7 的"静默截断"说法：`s > 2L` 是退化成多条 query 区间不相交的 chain，它们全是真比对，不能丢；
  - §4.5 明确 group 数 ≡ 非 circular 路径下该 read 的 record 数，noMar 与 linear 判定对齐；
  - §4.6 把 mapq 退化的判据从"是否跨 origin"改成"2L 上有没有第二条等价 chain"，并补实测：`s > L` 的 read
    同样会掉到 0；
  - §4.3.3 补注"primary 固定归 B 段"在 A / B 大小悬殊时会产出碎片化的 primary record；
  - §7 新增 primary 归属的待定项；§6 按 §8.1 扩充反例测试；新增 §8.1 实测数据。

- v2（2026-09-23，review 后修订）：
  - 新增"孪生去重 + group 分配"（§4.3.2）——原方案完全没提，会导致重复 record，且原定的 noMar 改法
    （按 fold-back 前的 hit 数）会把绝大多数不跨 origin 的 read 丢掉；
  - 新增 mapq 语义与决策表（§4.6）——复制序列对 mapq 的副作用原方案未提；
  - 更正 oupCovT / oupIyT 的过滤位置（在 `src/params.rs:306` 与 `src/bam_writer.rs:105`，不在 `align_worker`），
    并把 §5 里"`bam_writer.rs` 预计无改动"改成"必须改"（§4.5）；
  - 更正切列规则：原"I / D 列归入第一段"对 D 会产出非法 record（`reference_end > LN`），改为按 ref 碱基在 seam 处切开（§4.3.3）；
  - 更正 §4.1 / §4.5 中"两段 ovlp = 1.0"的说法：两段 query range 不相交（ovlp = 0），A 段保持 supplementary 的机制相反（§4.8）；
  - 明确实现路线：复用 `shrink_hit_core` 的"列区间 → 重写 hit"（含 reverse 链 query 坐标换算），
    并指出不要用 `src/aligned_pairs.rs:55` 的列遍历（§4.3.4）；
  - 两段的 `alignment_score` 改为沿用父 hit（原"按列重算 DP 分"无法与 minimap2 语义对齐）（§2）；
  - 明确 `fold_hits` 必须插在 `(Mapping, TargetCtx)` 成对的那一层，签名 `&mut Vec<(Mapping, TargetCtx)>`（§4.3）；
  - 新增适用域与限制（§4.7）、组合语义（§4.5）、实验依据（§8）、扩充测试计划（§6）。

## 10. 附：review 时发现的相邻既有问题（不属于本方案）

三次 review 里核对代码时顺带发现的三处，都不影响本方案的结论，但改到同一片代码时会碰到，单列出来：

1. `OupParams::modify_aligner`（`src/params.rs:331-336`）本意是让 `--noSeco` 影响 minimap2，但清的是
   `MM_F_SECONDARY_SEQ (0x1000000000)`——该位在 minimap2 crate 里从未被置位，`&=` 本身对未置位的位也是 no-op——
   而不是 `MM_F_NO_PRINT_2ND (0x4000)`。**结果是 `--noSeco` 只靠 writer 层的 record 过滤生效，minimap2 侧（含 mapq）
   完全没变**。§4.6 的结论不依赖它修好，但两处描述要对得上。
2. `query_seq_sender`（`src/lib.rs:286`）多文件 **fastq** 输入时把文件名后缀拼到了 `record.seq` 上（fasta 分支拼的是
   `record.name`），会把 `__1` 追加进序列、污染比对。
3. `--indexedTarget`（`src/bin/gsmm2.rs:101`）与 `--displayAlignment`（`src/bin/gsmm2.rs:126`）在 `alignment()` 里
   都没有被使用。
