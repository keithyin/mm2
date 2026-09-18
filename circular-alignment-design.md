# gsmm2 circular（环状参考）比对方案设计

状态：**方案已定，未动代码**。

> 注意：本文件暂存在 `projects/` 下，因为 `gsda/third_party/mm2/` 目录属 root、
> 仅 `bioinfo` 用户可写，当前会话（ubuntu 用户）无写入权限。
> 移入项目目录：`sudo cp /home/ubuntu/rocyin/projects/circular-alignment-design.md /home/ubuntu/rocyin/projects/gsda/third_party/mm2/`

## 1. 目标

给 gsmm2 增加 circular 参考序列的比对能力。基本实现方式：

- 比对阶段：用 `seq + seq`（2L）建 index，照常比对（chaining / DP / 全部后处理都在 2L 坐标系上完成）。
- 输出阶段：结果折回原始 `seq`（L）的坐标系；跨 origin（seam）的比对拆成两段分别输出。
- 正交性：现有所有功能（query-forward、polyNGapLeftAlign、hpTrShrinkAlnRegion、noSeco/noSupp/noMar、oupCovT/oupIyT 等）与 circular 组合均可用，且 circular 关闭时行为与现状完全一致（`fold_hits` 为空操作）。

## 2. 已定决策

| 决策点 | 结论 |
|---|---|
| primary / supplementary 归属 | 跨 seam 拆出的两段中，**pos=1 段为 primary**（保留原 hit 的 score / mapq），**尾部段为 supplementary**。两个 record 打配对 tag：`XT:A:Z:CIRC` + `XC:Z:<query_name>`，供下游把两段拼回。 |
| oupCovT / noMar 过滤时机 | **fold-back 之前**、hit 层面判断。跨 origin 的 hit 视为一条完整比对，coverage 按整条 read 计算，避免拆段后被误杀；`--noMar` 同理按 fold-back 前的 hit 数判断。 |
| CLI 开关形式 | 全局 `--circular`（bool），所有 target 都当 circular。 |

## 3. 现状管线（相关部分）

`src/bin/gsmm2.rs::alignment()`：

1. 读 target fasta → `Vec<ReadInfo>`。
2. 若 `--query-forward`，`targets_to_query_forward_targets`（`src/lib.rs:761`）展开成 `${name}___fwd` / `${name}___rev`。
3. `targets_to_targetsidx`（`src/lib.rs:742`）建 `name → (tid, len)`，用于 BAM header 的 SQ 记录（SN/LN）与 record.tid。
4. `build_aligner`（`src/lib.rs:124`）为每条展开后的 target 建一个 minimap2 index；`NoMemLeakAligner` 同时持有该 target 的 `target_seq`（供 cigar 后处理读碱基）与 `hp_tr_regions`（在同一序列上扫描 motif）。
5. `align_worker`（`src/lib.rs:306`）：`map_query_to_target`（map + DP）→ `hp_tr_shrink` → `poly_gap_left_align`（均通过 `target_seq[target_start..target_end]` 按 hit 坐标读 target 碱基）→ `set_primary_supp_alignment_according_2_align_score`（基于 query range overlap）。
6. `hits2records` → `build_bam_record_from_mapping`（`src/lib.rs:544`）：strand 决定 query 是否 RC；`convert_mapping_cigar_to_record_cigar` 加首尾 softclip；`pos = hit.target_start`；cs / md 原样透传。
7. `write_bam_worker`（`src/bam_writer.rs:57`）按 `targetsidx` 写 SQ header 后写 record。

关键点：**index 建完之后，下游全部处于"建 index 所用序列"的坐标系**。circular 方案就是利用这一点，在入口处换序列、出口处折回坐标。

## 4. 方案：入口展开，出口折回

与现有 `query_forward` 模式同构：circular 是一个"target 展开 + 出口折回坐标"的步骤，中间段（index 构建、chaining、DP、全部后处理）保持在 2L 坐标系不动。

### 4.1 流程

```
targets ──[1. circular: seq → seq+seq]──[2. query-forward: fwd/rev 展开]──▶ build_aligner(2L 序列, 不变)
                                                                                     │
align_worker: map ──▶ hp_tr_shrink ──▶ poly_gap_left_align          （2L 坐标系, 全部不变）
            ──▶ [3. fold_hits: 折回 L / 跨 seam 拆段]                （新增, circular 关闭时为空操作）
            ──▶ set_primary_supp_alignment_according_2_align_score（对拆段后的 hits 跑,
                                                        拆出的两段 ovlp=1.0, 不会互相降级）
            ──▶ noMar / oupCovT 过滤                                （hit 层面, 拆段前）
            ──▶ hits2records（每段独立建 record, softclip 按段自己的 [qs, qe) 加）
```

### 4.2 circular 展开

`targets_to_circular_targets(targets)`（新函数，样式照 `targets_to_query_forward_targets`）：

- 每条 target 的 `seq` 换成 `seq + seq`，**名字不变**（不需要后缀：同一 target 只有一份，且名字要原样出现在输出 BAM 里）。
- 执行顺序：circular 展开**在 query-forward 展开之前**。

与 query-forward 组合时每个 target 展开为两条链（执行顺序无关，因为 RC(S+S) = RC(S)+RC(S)）：

- `t___fwd`：序列 `S+S`，长 2L
- `t___rev`：序列 `RC(S+S)`，长 2L

### 4.3 fold_hits

在 `align_single_query_to_targets` / `align_single_query_to_targets_query_forward` 的后处理（hp_tr_shrink、poly_gap_left_align）之后、`set_primary_supp_alignment_according_2_align_score` 之前调用。circular 关闭时直接原样返回。

对每个 hit（`target_start < target_end`，L 为原始长度）：

1. **`target_end ≤ L`**：原样保留。
2. **`target_start ≥ L`**：两端各减 L。
3. **跨 seam：`target_start < L < target_end`**：拆成两个 hit。
   - 切点 `qm` 由逐列走 cigar 确定：从比对首列开始累计 target 消耗量（M / X / D 列各消耗 1，I 列不消耗），累计到 `L - target_start` 的那个列就是切点；I / D 列归入第一段（seam 处 target 的最后一列）。
   - 段 A `[target_start, L)` × `[qs, qm)`：**supplementary**，`alignment_score` = 自己覆盖列的 DP 分。
   - 段 B `[0, target_end - L)` × `[qm, qe)`：**primary**，保留原 hit 的总分与 mapq。
   - 两段的 sub-cigar / sub-cs / sub-md 按列重建（复用 `cigar_cvt.rs` 的列级转换逻辑与 `mapping_ext::aligned_2_str` 的列遍历逻辑；cs 的 `*` / `-` 字符按列重建）。

#### 与 query-forward 组合时不需要镜像换算

- fwd 链：`S+S[p mod L] == S[p mod L]`，直接 `p mod L` 折回。
- rev 链：hit 的 pos 本来就是 `RC(S+S)` 坐标系（现有 `___rev` 机制即如此，`build_bam_record_from_mapping` 直接以 `hit.target_start` 作 pos）。而 `RC(S+S)[p mod L] == RC(S)[p mod L]` 同样成立，所以折回 L 坐标系后 record 落在 `t___rev` 的同一 0-based 位置、reverse strand。**两条链同构，都是 `p mod L`，无任何额外换算。**
- 跨 origin 的 split 点也不变：两条链都用「target 坐标 = L」作 seam（fwd 链是 S 的 origin，rev 链是 RC 坐标下的 origin）。

### 4.4 输出坐标与 header

- `record.pos` 为 0-based ∈ `[0, L)`；1-based 即 `[1, L]`，跨 origin read 的第二段 `pos = 1`。
- header 的 SQ `LN` 用**原始长度 L**：`targets_to_targetsidx` 在 circular 模式下记录原始长度（展开前的 seq.len()），`write_bam_worker` 不需要改（LN 本来就来自 targetsidx）。

### 4.5 正交性说明

| 功能 | 与 circular 的交互 | 需要的处理 |
|---|---|---|
| polyNGapLeftAlign | 在 2L 坐标系上操作，同聚物跨 seam 时 2L 序列已含完整同聚物（首尾拼接），左移反而比 linear 更正确。fold-back 在左移之后，坐标自然正确 | 无 |
| hpTrShrinkAlnRegion | motif 表扫在 S+S 上，跨 seam 的 motif 天然存在；shrink 可能把跨 seam 的 hit 切成两段（各带 softclip），fold-back 在其后，顺序天然正确 | 无 |
| noSeco / noSupp | 作用在 record 层面；noSupp 会把 CIRC 对的 supp 段删掉（预期行为：只剩 primary 段） | 无 |
| oupCovT / noMar | 见 §2，均在 fold-back 前的 hit 层面判断 | **过滤点前移**：现状 `noMar` 在 `align_worker` 里按 `records.len() > 1`（src/lib.rs:322）判断，拆段会让 CIRC 对变成 2 条 record 被误杀，需改为按 fold-back 前的 hit 数判断。该改动对非 circular 路径无行为变化 |
| query-forward | 见 §4.2 / §4.3，fold-back 同构，无镜像换算 | 无 |
| primary/supp 判定 | `set_primary_supp_alignment_according_2_align_score` 对拆段后的 hits 运行；CIRC 两段 query range 完全重叠（ovlp=1.0 > 0.5），supp 段不会与 primary 之外的其它 supp 产生冲突判定 | 无 |

## 5. 改动文件清单（预告）

| 文件 | 内容 |
|---|---|
| `src/params.rs` | `MapParams` 加 `circular: bool` + `set_circular` |
| `src/lib.rs` | `targets_to_circular_targets`；`align_single_query_to_targets*` 在 post-processing 后、primary 判定前插 `fold_hits`；`align_worker` 的 noMar 过滤改为 fold-back 前的 hit 数 |
| 新 `src/circular.rs` | `fold_hits` + 按列拆 cigar / cs / md 的列级工具 + 单测 |
| `src/bin/gsmm2.rs` | `--circular` arg；circular 展开调用（放在 query-forward 展开之前）；`targets_to_targetsidx` 在 circular 模式下记录原始长度 |
| `src/bam_writer.rs` | 预计无改动（LN 来自 targetsidx） |

## 6. 测试计划

1. **fold_hits 单测**：手工构造三种 hit（不跨 seam / 整体在第二份 / 跨 seam），验证折回坐标、拆段后的 cigar / cs / md、primary 归属、score 分配。
2. **端到端**：200bp circular target，构造跨 origin 的 read（后 100bp + 前 100bp）→ 应出一条 CIRC 对（pos=1 primary + 尾部 supp），两段列数和 == read 长度；同一 read 比 linear 参考作对照（应出现大 gap 或两个低分 hit）。
3. **正交矩阵**：circular × { query-forward（fwd / rev 两条 read）、polyNGapLeftAlign（同聚物跨 seam）、hpTrShrinkAlnRegion、noSupp、noSeco、noMar、oupCovT } 各跑一个 case。
4. **回归**：现有全部测试不回归（circular 默认关闭，`fold_hits` 为空操作）。

## 7. 实现时待定的小细节

- noMar 由 `records.len() > 1` 改为 fold-back 前的 hit 数判断，具体在 `align_worker` 的哪一层做（`align_single_query_to_targets*` 返回后、`hits2records` 前）。
- sub-cs 重建时 `*`（I 列）与 `-`（D 列）在切点列的归属（已定：归第一段，实现时核对 minimap2 cs 字符语义）。
- CIRC 配对 tag 的最终键名（暂定 `XT:Z:CIRC` + `XC:Z:<qname>`；是否需要在 supp 段额外存 partner 的 pos 由下游需求定）。
