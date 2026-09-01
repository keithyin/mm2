pub mod align_processor;
pub mod bam_writer;
pub mod cigar_cvt;
pub mod mapping_ext;
pub mod params;
use std::{
    cmp,
    collections::HashMap,
    ops::{Deref, DerefMut},
    sync::Arc,
    thread,
};

use cigar_adjust::cigar_adjust_poly_gap_left_align;
use crossbeam::channel::{Receiver, Sender};
use gskits::{
    dna::reverse_complement,
    ds::ReadInfo,
    fastx_reader::{fasta_reader::FastaFileReader, fastq_reader::FastqReader},
    phreq::phreq_list_2_quality,
};
use hp_tr_shrink::HpTrRegions;
use mapping_ext::MappingExt;
use minimap2::{Aligner, Built, Mapping};
use params::{
    AlignParams, IndexParams, InputFilterParams, MapParams, OupParams, TOverrideAlignerParam,
};
use rust_htslib::bam::{
    record::{Aux, AuxArray, Cigar, CigarString},
    Read,
};

pub use gskits;
pub use minimap2;
pub mod aligned_pairs;
pub mod cigar_adjust;
pub mod hp_tr_shrink;
pub mod visualization;

pub type BamRecord = rust_htslib::bam::record::Record;
pub type BamWriter = rust_htslib::bam::Writer;
pub type BamReader = rust_htslib::bam::Reader;

pub struct AlignResult {
    pub records: Vec<BamRecord>,
}

/// 单个 target 的 aligner 包装。除了 aligner 本身, 还保存该 aligner 的 index 所对应的
/// target 序列, 以便对 hits 的 cigar 做后处理时需要读取 target 上的碱基
/// (见 `cigar_adjust::cigar_adjust_poly_gap_left_align`), 以及该 target 上的 HP/TR motif 表
/// (见 `hp_tr_shrink::shrink_hit`, 只在 --hpTrShrinkAlnRegion 打开时构建)
pub struct NoMemLeakAligner {
    pub aligner: Aligner<Built>,
    pub target_seq: Vec<u8>,
    pub hp_tr_regions: Option<Arc<HpTrRegions>>,
}

/// 一个 hit 所属 target 的后处理上下文。`target_seq` 与 `hp_tr_regions` 都在该 aligner
/// 实际建 index 的那条序列上取, 因此与 hit 的 target_start/target_end 同坐标系
/// (--query-forward 展开出的 `___rev` 链也是如此, 不需要额外换算)
#[derive(Clone, Copy)]
pub struct TargetCtx<'a> {
    pub target_seq: &'a [u8],
    pub hp_tr_regions: Option<&'a HpTrRegions>,
}

impl Default for TargetCtx<'_> {
    fn default() -> Self {
        Self {
            target_seq: &[],
            hp_tr_regions: None,
        }
    }
}

impl NoMemLeakAligner {
    pub fn new(aligner: Aligner<Built>, target_seq: &[u8]) -> Self {
        Self {
            aligner,
            target_seq: target_seq.to_vec(),
            hp_tr_regions: None,
        }
    }

    pub fn with_hp_tr_regions(mut self, regions: Option<Arc<HpTrRegions>>) -> Self {
        self.hp_tr_regions = regions;
        self
    }

    pub fn target_seq(&self) -> &[u8] {
        &self.target_seq
    }

    pub fn target_ctx(&self) -> TargetCtx<'_> {
        TargetCtx {
            target_seq: &self.target_seq,
            hp_tr_regions: self.hp_tr_regions.as_deref(),
        }
    }
}

impl Deref for NoMemLeakAligner {
    type Target = Aligner<Built>;
    fn deref(&self) -> &Self::Target {
        &self.aligner
    }
}

impl DerefMut for NoMemLeakAligner {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.aligner
    }
}

impl Drop for NoMemLeakAligner {
    fn drop(&mut self) {
        // let idx = self.aligner.idx.take().unwrap();
        // // unsafe {
        // //     mm_idx_destroy(*idx.as_ref());
        // // }
    }
}

pub fn build_aligner(
    preset: &str,
    index_args: &IndexParams,
    map_args: &MapParams,
    align_args: &AlignParams,
    oup_args: &OupParams,
    targets: &Vec<ReadInfo>,
    threads: usize,
) -> Vec<NoMemLeakAligner> {
    let per_target_threads = (threads / targets.len()).max(1);
    let hp_tr_shrink = map_args.hp_tr_shrink_aln_region;

    let aligners = thread::scope(|s| {
        let mut handles = vec![];
        for target in targets {
            let hd = s.spawn(|| {
                let aligner = Aligner::builder();

                let mut aligner = match preset {
                    "map-ont" => aligner.map_ont(),
                    "map-pb" => aligner.map_pb(),
                    "map-hifi" => aligner.map_hifi(),

                    pre => panic!("not implemented yet {}", pre),
                };

                index_args.modify_aligner(&mut aligner);
                map_args.modify_aligner(&mut aligner);
                align_args.modify_aligner(&mut aligner);
                oup_args.modify_aligner(&mut aligner);

                let mut aligner = aligner
                    .with_index_threads(per_target_threads)
                    .with_cigar()
                    .with_sam_out()
                    .with_sam_hit_only()
                    .with_seq_and_id(target.seq.as_bytes(), target.name.as_bytes())
                    .unwrap();

                // https://github.com/lh3/minimap2/blob/618d33515e5853c4576d5a3d126fdcda28f0e8a4/options.c#L43
                aligner.mapopt.best_n = 5; // top best_n chains are subjected to DP alignment
                                           // aligner.mapopt.pri_ratio = 0.1; // secondary-to-primary score ratio to output secondary mappings

                // motif 表在 aligner 实际索引的那条序列上扫描, 与 hit 的 target 坐标同系。
                // 扫描代价是每条 target 一次性的(约 340 条正则各过一遍序列), 由这里的
                // per-target 线程并行分摊
                let regions = if hp_tr_shrink {
                    Some(hp_tr_shrink::build_hp_tr_regions(target.seq.as_str()))
                } else {
                    None
                };

                (aligner, regions)
            });
            handles.push(hd);
        }

        handles
            .into_iter()
            // spawn 顺序与 targets 一致, 按下标配对即可拿到 aligner 对应的 target 序列
            .zip(targets.iter())
            .map(|(hd, target)| {
                let (aligner, regions) = hd.join().unwrap();
                NoMemLeakAligner::new(aligner, target.seq.as_bytes()).with_hp_tr_regions(regions)
            })
            .collect::<Vec<_>>()
    });

    aligners
}

/// does not work now
pub fn build_aligner_v2(
    preset: &str,
    index_args: &IndexParams,
    map_args: &MapParams,
    align_args: &AlignParams,
    oup_args: &OupParams,
    targets: &Vec<ReadInfo>,
) -> NoMemLeakAligner {
    let aligner = Aligner::builder();

    let mut aligner = match preset {
        "map-ont" => aligner.map_ont(),
        "map-pb" => aligner.map_pb(),
        "map-hifi" => aligner.map_hifi(),

        pre => panic!("not implemented yet {}", pre),
    };

    index_args.modify_aligner(&mut aligner);
    map_args.modify_aligner(&mut aligner);
    align_args.modify_aligner(&mut aligner);
    oup_args.modify_aligner(&mut aligner);

    let (seqs, ids): (Vec<Vec<u8>>, Vec<Vec<u8>>) = targets
        .iter()
        .map(|read_info| {
            (
                read_info.seq.as_bytes().to_vec(),
                read_info.name.as_bytes().to_vec(),
            )
        })
        .unzip();

    let mut aligner = aligner
        .with_index_threads(4)
        .with_cigar()
        .with_sam_out()
        .with_sam_hit_only()
        .with_seqs_and_ids(&seqs, &ids)
        .unwrap();

    // https://github.com/lh3/minimap2/blob/618d33515e5853c4576d5a3d126fdcda28f0e8a4/options.c#L43
    aligner.mapopt.best_n = 5; // top best_n chains are subjected to DP alignment
                               // aligner.mapopt.pri_ratio = 0.1; // secondary-to-primary score ratio to output secondary mappings

    // 一个 aligner 对应多个 target, 无法给出单一的 target_seq, 依赖它的 cigar 后处理不可用
    NoMemLeakAligner::new(aligner, &[])
}

pub fn query_seq_sender(
    filenames: &Vec<String>,
    sender: Sender<ReadInfo>,
    input_filter_params: &InputFilterParams,
    oup_params: &OupParams,
) {
    let mut file_idx = 0;
    for filename in filenames {
        let qname_suffix = if filenames.len() == 1 {
            None
        } else {
            Some(format!("__{}", file_idx))
        };
        if filename.ends_with("fa") || filename.ends_with("fasta") || filename.ends_with("fna") {
            let fasta_reader = FastaFileReader::new(filename.clone());

            for mut record in fasta_reader {
                if let Some(suffix) = &qname_suffix {
                    record.name.push_str(suffix);
                }
                sender.send(record).unwrap();
            }
        } else if filename.ends_with("bam") {
            let mut bam_h = BamReader::from_path(filename).unwrap();
            bam_h.set_threads(4).unwrap();
            for record in bam_h.records() {
                let record = record.unwrap();
                if input_filter_params.valid(&record) {
                    sender
                        .send(ReadInfo::from_bam_record(
                            &record,
                            qname_suffix.as_ref().map(|v| v.as_str()),
                            &oup_params.pass_through_tags,
                        ))
                        .unwrap();
                }
            }
        } else if filename.ends_with("fq") || filename.ends_with("fastq") {
            let fa_iter = FastqReader::new(filename.to_string());
            for mut record in fa_iter {
                if let Some(suffix) = &qname_suffix {
                    record.seq.push_str(suffix);
                }
                let rq = phreq_list_2_quality(record.qual.as_ref().unwrap()).unwrap_or(0.);

                if input_filter_params.rq_valid(rq) {
                    record.rq = Some(rq);
                    sender.send(record).unwrap();
                }
            }
        } else {
            panic!(
                "invalid file format {}. bam/fa/fasta/fna supported",
                filename
            );
        }

        file_idx += 1;
    }
}

pub fn align_worker(
    query_record_recv: Receiver<gskits::ds::ReadInfo>,
    align_res_sender: Sender<AlignResult>,
    aligners: &Vec<NoMemLeakAligner>,
    target_idx: &HashMap<String, (usize, usize)>,
    oup_params: &OupParams,
    map_params: &MapParams,
) {
    for query_record in query_record_recv {
        let hits = if map_params.query_forward {
            align_single_query_to_targets_query_forward(&query_record, aligners, map_params)
        } else {
            align_single_query_to_targets(&query_record, aligners, map_params)
        };
        let records = hits2records(&hits, &query_record, target_idx);

        if oup_params.discard_multi_align_reads && records.len() > 1 {
            continue;
        }

        if records.is_empty() {
            continue;
        }

        align_res_sender.send(AlignResult { records }).unwrap();
    }
}

/// 将 query 比对到单个 aligner(即单个 target) 上, 过滤掉没有 alignment 的 hit
fn map_query_to_target(
    aligner: &NoMemLeakAligner,
    query_record: &gskits::ds::ReadInfo,
) -> Vec<Mapping> {
    aligner
        .map(
            query_record.seq.as_bytes(),
            true,
            true,
            None,
            Some(&[67108864]), // 67108864 eqx
            Some(query_record.name.as_bytes()),
        )
        .unwrap()
        .into_iter()
        .filter(|hit| hit.alignment.is_some()) // filter unmapped
        .collect()
}

/// 对 hits 做多聚物区 gap 左移规范化。query 传入原始方向的序列
/// (`aligned_2_str` 会按 hit 的 strand 自行 reverse complement)
fn poly_gap_left_align_hits(hits: &mut [(Mapping, TargetCtx)], query_seq: &[u8]) {
    for (hit, ctx) in hits.iter_mut() {
        cigar_adjust_poly_gap_left_align(hit, ctx.target_seq, query_seq);
    }
}

/// 收紧 hits 的比对区域: 比对末端落在 HP/TR 内部的那一段被裁掉。未开启
/// --hpTrShrinkAlnRegion 时 target_ctx.hp_tr_regions 为 None, 这里什么也不做
fn shrink_hp_tr_aln_region_hits(hits: &mut [(Mapping, TargetCtx)], query_seq: &[u8]) {
    for (hit, ctx) in hits.iter_mut() {
        if let Some(regions) = ctx.hp_tr_regions {
            hp_tr_shrink::shrink_hit(hit, ctx.target_seq, query_seq, regions);
        }
    }
}

pub fn align_single_query_to_targets(
    query_record: &gskits::ds::ReadInfo,
    aligners: &Vec<NoMemLeakAligner>,
    map_params: &MapParams,
) -> Vec<Mapping> {
    // hit 与其所属 target 的后处理上下文成对保存: gap 左移与比对区域收紧都要读 target 上的碱基
    let mut all_hits: Vec<(Mapping, TargetCtx)> = vec![];

    for aligner in aligners {
        all_hits.extend(
            map_query_to_target(aligner, query_record)
                .into_iter()
                .map(|hit| (hit, aligner.target_ctx())),
        );
    }

    // 先收紧窗口再左移 gap: 反过来的话 gap 可能被左移到随后被裁掉的列上
    if map_params.hp_tr_shrink_aln_region {
        shrink_hp_tr_aln_region_hits(&mut all_hits, query_record.seq.as_bytes());
    }

    if map_params.poly_n_gap_left_align {
        poly_gap_left_align_hits(&mut all_hits, query_record.seq.as_bytes());
    }

    let mut all_hits: Vec<Mapping> = all_hits.into_iter().map(|(hit, _)| hit).collect();

    // set_primary_alignment_according_2_align_score(&mut all_hits);
    set_primary_supp_alignment_according_2_align_score(&mut all_hits);
    all_hits
}

/// query_forward 模式下的比对。
///
/// aligners 需要由 `targets_to_query_forward_targets` 产出的 target 构建, 即每个真实的
/// target 都有 `${name}___fwd` / `${name}___rev` 两份。对同一个 target 的这两份结果做
/// 二选一, 只保留 primary alignment 的 query strand 为 forward 的那份, 从而保证输出的
/// record 中 query 始终处于原始方向(不会被 reverse complement)。
///
/// - 两份都是 forward: 取 primary score 高的那份
/// - 两份都是 reverse(或两侧都没有 hit): 该 target 的结果被过滤掉
///
/// 胜出方的 secondary/supplementary hits 原样保留。
pub fn align_single_query_to_targets_query_forward(
    query_record: &gskits::ds::ReadInfo,
    aligners: &Vec<NoMemLeakAligner>,
    map_params: &MapParams,
) -> Vec<Mapping> {
    #[derive(Default)]
    struct TargetHits<'a> {
        fwd: Vec<Mapping>,
        rev: Vec<Mapping>,
        // 每一侧的 hits 都来自唯一一个 aligner(即唯一一个 target), 记录其上下文供
        // gap 左移与比对区域收紧使用
        fwd_ctx: TargetCtx<'a>,
        rev_ctx: TargetCtx<'a>,
    }

    // 按 target 首次出现的顺序排列, 保证结果可复现
    let mut grouped: Vec<(String, TargetHits)> = vec![];
    let mut name2pos: HashMap<String, usize> = HashMap::new();

    for aligner in aligners {
        let mut hits = map_query_to_target(aligner, query_record);
        if hits.is_empty() {
            // 该 target 上没有 hit, 无需参与 pk
            continue;
        }

        let (base_name, is_fwd) = split_target_strand_name(
            hits.first()
                .and_then(|hit| hit.target_name.as_ref())
                .map(|name| name.as_str())
                .unwrap_or(""),
        );

        let pos = *name2pos.entry(base_name.clone()).or_insert_with(|| {
            grouped.push((base_name, TargetHits::default()));
            grouped.len() - 1
        });
        let target_hits = &mut grouped[pos].1;
        if is_fwd {
            target_hits.fwd.append(&mut hits);
            target_hits.fwd_ctx = aligner.target_ctx();
        } else {
            target_hits.rev.append(&mut hits);
            target_hits.rev_ctx = aligner.target_ctx();
        }
    }

    // hit 与其所属 target 的后处理上下文成对保存
    let mut all_hits: Vec<(Mapping, TargetCtx)> = vec![];
    for (_, target_hits) in grouped {
        // 两侧都可能有 hit, pk 之后再取用胜出方的上下文, 避免调整被丢弃的 hits
        let keep_fwd = pk_target_strands(&target_hits.fwd, &target_hits.rev);
        let (hits, ctx) = match keep_fwd {
            Some(true) => (target_hits.fwd, target_hits.fwd_ctx),
            Some(false) => (target_hits.rev, target_hits.rev_ctx),
            // 两份 index 上的 primary 都不在 query forward 链上, 无法保持 query 方向不变
            None => continue,
        };

        all_hits.extend(hits.into_iter().map(|hit| (hit, ctx)));
    }

    // 先收紧窗口再左移 gap, 理由同 align_single_query_to_targets
    if map_params.hp_tr_shrink_aln_region {
        shrink_hp_tr_aln_region_hits(&mut all_hits, query_record.seq.as_bytes());
    }

    if map_params.poly_n_gap_left_align {
        poly_gap_left_align_hits(&mut all_hits, query_record.seq.as_bytes());
    }

    let mut all_hits: Vec<Mapping> = all_hits.into_iter().map(|(hit, _)| hit).collect();

    set_primary_supp_alignment_according_2_align_score(&mut all_hits);
    all_hits
}

/// 将 `${name}___fwd` / `${name}___rev` 拆成 (name, 是否为正向链)。
/// 不带后缀的 target 名视为正向链, 即没有对手参与 pk。
fn split_target_strand_name(target_name: &str) -> (String, bool) {
    if let Some(base_name) = target_name.strip_suffix(QUERY_FORWARD_FWD_SUFFIX) {
        (base_name.to_string(), true)
    } else if let Some(base_name) = target_name.strip_suffix(QUERY_FORWARD_REV_SUFFIX) {
        (base_name.to_string(), false)
    } else {
        (target_name.to_string(), true)
    }
}

/// 对同一 target 的 fwd / rev 两组 hits 二选一。
/// `Some(true)` 取 fwd 组, `Some(false)` 取 rev 组, `None` 表示两侧的 primary 都不在
/// query forward 链上(该 target 的结果需要被过滤掉)。
fn pk_target_strands(fwd: &[Mapping], rev: &[Mapping]) -> Option<bool> {
    let (fwd_primary, rev_primary) = (primary_hit(fwd), primary_hit(rev));

    match (
        fwd_primary.map(|hit| hit.strand),
        rev_primary.map(|hit| hit.strand),
    ) {
        (Some(minimap2::Strand::Forward), Some(minimap2::Strand::Forward)) => {
            Some(alignment_score(fwd_primary.unwrap()) >= alignment_score(rev_primary.unwrap()))
        }
        (Some(minimap2::Strand::Forward), _) => Some(true),
        (_, Some(minimap2::Strand::Forward)) => Some(false),
        _ => None,
    }
}

/// 单个 target 的 hits 中的 primary, 没有 primary 时退化为 score 最高的 hit
fn primary_hit(hits: &[Mapping]) -> Option<&Mapping> {
    hits.iter()
        .find(|hit| hit.is_primary)
        .or_else(|| hits.iter().max_by_key(|hit| alignment_score(hit)))
}

fn alignment_score(hit: &Mapping) -> i32 {
    hit.alignment.as_ref().unwrap().alignment_score.unwrap()
}

pub fn hits2records(
    hits: &Vec<Mapping>,
    query_record: &gskits::ds::ReadInfo,
    target_idx: &HashMap<String, (usize, usize)>,
) -> Vec<BamRecord> {
    hits.iter()
        .map(|hit| build_bam_record_from_mapping(hit, query_record, target_idx))
        .collect()
}

pub fn build_bam_record_from_mapping(
    hit: &minimap2::Mapping,
    query_record: &gskits::ds::ReadInfo,
    target_idx: &HashMap<String, (usize, usize)>,
) -> BamRecord {
    // println!("{:?}", hit);

    let mut bam_record = BamRecord::new();

    let mut seq = &query_record.seq;
    let mut is_rev = false;

    let rev_seq = match hit.strand {
        minimap2::Strand::Forward => None,
        minimap2::Strand::Reverse => {
            is_rev = true;
            Some(String::from_utf8(reverse_complement(seq.as_bytes())).unwrap())
        }
    };

    if is_rev {
        seq = rev_seq.as_ref().unwrap();
        bam_record.set_reverse();
    }

    // hit.query_start, hit.query_end 是相对于原始 query 而言的(即 未 reverse 的 query 而言)
    let aln_info = hit.alignment.as_ref().unwrap();
    let cigar_str = convert_mapping_cigar_to_record_cigar(
        aln_info.cigar.as_ref().unwrap(),
        hit.query_start as usize,
        hit.query_end as usize,
        seq.len(),
        is_rev,
    );

    let qual = if let Some(ref qual_) = query_record.qual {
        if is_rev {
            qual_.iter().copied().rev().collect()
        } else {
            qual_.clone()
        }
    } else {
        vec![255; seq.len()]
    };

    bam_record.set(
        query_record.name.as_bytes(),
        Some(&cigar_str),
        seq.as_bytes(),
        &qual,
    );
    if is_rev {
        bam_record.set_reverse();
    }

    // reference start
    bam_record.set_pos(hit.target_start as i64);
    bam_record.set_mpos(-1);
    // bam_record.set_mpos(mpos);
    // bam_record.reference_end()
    bam_record.set_mapq(hit.mapq as u8);

    bam_record.set_tid(
        target_idx
            .get(hit.target_name.as_ref().unwrap().as_str())
            .unwrap()
            .0 as i32,
    );
    bam_record.set_mtid(-1);

    if hit.is_primary {
        bam_record.unset_secondary();
        bam_record.unset_supplementary();
    } else {
        bam_record.set_secondary();
        bam_record.set_supplementary();
    }

    if hit.is_supplementary {
        bam_record.unset_secondary();
        bam_record.set_supplementary();
    } else {
        bam_record.unset_supplementary();
    }
    bam_record.unset_unmapped();

    if let Some(cs) = aln_info.cs.as_ref() {
        bam_record
            .push_aux(b"cs", rust_htslib::bam::record::Aux::String(cs))
            .unwrap();
    }

    if let Some(md) = aln_info.md.as_ref() {
        bam_record
            .push_aux(b"md", rust_htslib::bam::record::Aux::String(md))
            .unwrap();
    }

    if let Some(np_) = query_record.np {
        bam_record.push_aux(b"np", Aux::U16(np_ as u16)).unwrap();
    }

    if let Some(ch_) = query_record.ch {
        bam_record.push_aux(b"ch", Aux::U32(ch_ as u32)).unwrap();
    }

    if let Some(rq_) = query_record.rq {
        bam_record.push_aux(b"rq", Aux::Float(rq_)).unwrap();
    }

    if let Some(dw_) = &query_record.dw {
        if is_rev {
            let dw_ = dw_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"dw", Aux::ArrayU8(AuxArray::from(&dw_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"dw", Aux::ArrayU8(AuxArray::from(dw_)))
                .unwrap();
        }
    }

    if let Some(ar_) = &query_record.ar {
        if is_rev {
            let ar_ = ar_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"ar", Aux::ArrayU8(AuxArray::from(&ar_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"ar", Aux::ArrayU8(AuxArray::from(ar_)))
                .unwrap();
        }
    }

    if let Some(cr_) = &query_record.cr {
        if is_rev {
            let cr_ = cr_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"cr", Aux::ArrayU8(AuxArray::from(&cr_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"cr", Aux::ArrayU8(AuxArray::from(cr_)))
                .unwrap();
        }
    }

    if let Some(nn_) = &query_record.nn {
        if is_rev {
            let nn_ = nn_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"nn", Aux::ArrayU8(AuxArray::from(&nn_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"nn", Aux::ArrayU8(AuxArray::from(nn_)))
                .unwrap();
        }
    }

    if let Some(be_) = &query_record.be {
        bam_record
            .push_aux(b"be", Aux::ArrayU32(AuxArray::from(be_)))
            .unwrap();
    }

    bam_record
}

/// if reverse alignment, the cigar will be reversed!
pub fn convert_mapping_cigar_to_record_cigar(
    mapping_cigar: &[(u32, u8)],
    mut query_start: usize,
    mut query_end: usize,
    query_len: usize,
    is_rev: bool,
) -> CigarString {
    let mut cigar_str = CigarString(vec![]);

    if is_rev {
        (query_start, query_end) = (query_len - query_end, query_len - query_start);
    }

    if query_start > 0 {
        cigar_str.push(Cigar::SoftClip(query_start as u32));
    }

    cigar_str.extend(cigar_cvt::mapping_cigar2htslib_cigar_str(mapping_cigar).0);

    if query_end != query_len {
        cigar_str.push(Cigar::SoftClip((query_len - query_end) as u32));
    }
    cigar_str
}

/// {"target_name": (idx, length)}
pub fn targets_to_targetsidx(targets: &Vec<ReadInfo>) -> HashMap<String, (usize, usize)> {
    let mut target2idx = HashMap::new();
    targets.iter().enumerate().for_each(|(idx, target)| {
        target2idx.insert(target.name.clone(), (idx, target.seq.len()));
    });
    target2idx
}

/// query_forward 模式下, 正向链 target 的名字后缀
pub const QUERY_FORWARD_FWD_SUFFIX: &str = "___fwd";
/// query_forward 模式下, 反向链 target 的名字后缀
pub const QUERY_FORWARD_REV_SUFFIX: &str = "___rev";

/// 将每个 target 展开为两份: `${name}___fwd`(原始序列) 和 `${name}___rev`
/// (reverse complement), 两份交错排列。
///
/// 用于 query_forward 模式: 同一个 target 的两条链都建了 index, 比对时总可以选到
/// query 处于 forward 链的那份, 详见
/// `align_single_query_to_targets_query_forward`。
pub fn targets_to_query_forward_targets(targets: &[ReadInfo]) -> Vec<ReadInfo> {
    let mut res = Vec::with_capacity(targets.len() * 2);

    for target in targets {
        res.push(ReadInfo {
            name: format!("{}{}", target.name, QUERY_FORWARD_FWD_SUFFIX),
            seq: target.seq.clone(),
            ..Default::default()
        });

        res.push(ReadInfo {
            name: format!("{}{}", target.name, QUERY_FORWARD_REV_SUFFIX),
            seq: rev_complement_for_index(target.seq.as_str()),
            ..Default::default()
        });
    }

    res
}

/// reverse complement 一个 target, 使得结果可以直接用于构建 minimap2 index。
///
/// gskits::dna::reverse_complement 会把 ACGTN/acgtn/-/* 之外的字符补成 NUL, 而 index
/// 构建接口用 CStr 传序列, 遇到 NUL 会 panic, 这里统一将其替换成 N。
fn rev_complement_for_index(seq: &str) -> String {
    let mut seq = reverse_complement(seq.as_bytes());
    seq.iter_mut()
        .for_each(|b| *b = if *b == 0 { b'N' } else { *b });

    String::from_utf8(seq).unwrap()
}

pub fn set_primary_alignment_according_2_align_score(hits: &mut Vec<Mapping>) {
    if hits.is_empty() {
        return;
    }
    hits.sort_by_key(|v| {
        -v.alignment.as_ref().unwrap().alignment_score.unwrap()
            - if v.is_primary { 1_000_000_000 } else { 0 }
    });

    assert!(hits.first_mut().unwrap().is_primary); // assertion failed
    assert!(!hits.first_mut().unwrap().is_supplementary);

    let primary_hit = hits.first().unwrap();
    let (primary_qstart, primary_qend) = (primary_hit.query_start, primary_hit.query_end);

    hits.iter_mut().skip(1).for_each(|hit| {
        if hit.is_primary {
            let ratio = ovlp_ratio(primary_qstart, primary_qend, hit.query_start, hit.query_end);
            hit.is_primary = false;
            hit.is_supplementary = ratio <= 0.5;
        }
    });
}

/// 重新设置 primary 和 supplementary。其主要在 多个参考基因序列场景使用
/// 由于对于多个参考序列场景，每个参考序列都构建了一个aligner，所以我们会得到之多 N 个 primary。N表示参考序列个数。
/// primary设置逻辑：primary中，比对分数最高的作为最终的primary, 其余的primary基于 query 的 ovlp 来确定自己是 supplementary 还是 secondary
/// supplementary设置逻辑：所有hits重新排列，按照 primary 在第一位，supplementary 按照 比对得分位列2~？。secondary排到最后
///     依次判断 当前 range 和 之前的 supplementary range的 ovlp。如果 ovlp>0.5 当前 supp 就降级为 secondary
pub fn set_primary_supp_alignment_according_2_align_score(hits: &mut Vec<Mapping>) {
    if hits.is_empty() {
        return;
    }
    hits.sort_by_key(|v| {
        -v.alignment.as_ref().unwrap().alignment_score.unwrap()
            - if v.is_primary { 1_000_000_000 } else { 0 }
    });

    assert!(hits.first_mut().unwrap().is_primary); // assertion failed
    assert!(!hits.first_mut().unwrap().is_supplementary);

    let primary_hit = hits.first().unwrap();
    let (primary_qstart, primary_qend) = (primary_hit.query_start, primary_hit.query_end);

    hits.iter_mut().skip(1).for_each(|hit| {
        if hit.is_primary {
            let ratio = ovlp_ratio(primary_qstart, primary_qend, hit.query_start, hit.query_end);
            hit.is_primary = false;
            hit.is_supplementary = ratio <= 0.5;
        }
    });

    hits.sort_by_key(|v| {
        -v.alignment.as_ref().unwrap().alignment_score.unwrap()
            - if v.is_primary {
                1_000_000_000
            } else if v.is_supplementary {
                100_000_000
            } else {
                0
            }
    });

    let mut visited_ranges = vec![(primary_qstart, primary_qend)];
    hits.iter_mut().skip(1).for_each(|v| {
        if v.is_supplementary {
            let (start, end) = (v.query_start, v.query_end);

            let cnt = visited_ranges
                .iter()
                .map(|&(range_s, range_e)| ovlp_ratio(range_s, range_e, start, end))
                .filter(|&ovlp| ovlp > 0.5)
                .count();
            if cnt > 0 {
                v.is_supplementary = false;
            } else {
                visited_ranges.push((start, end));
            }
        }
    });
}

pub fn mapping2str(hit: &Mapping) -> String {
    format!(
        "qstart:{}, qend:{}, primary:{}, supp:{}, identity:{}, score:{:?}",
        hit.query_start,
        hit.query_end,
        hit.is_primary,
        hit.is_supplementary,
        MappingExt(hit).identity(),
        hit.alignment.as_ref().unwrap().alignment_score
    )
}

#[allow(unused)]
fn get_query_start_end(hit: &Mapping, query_len: i32) -> (i32, i32) {
    match hit.strand {
        minimap2::Strand::Forward => (hit.query_start, hit.query_end),
        minimap2::Strand::Reverse => (query_len - hit.query_end, query_len - hit.query_start),
    }
}

fn ovlp_ratio(s1: i32, e1: i32, s2: i32, e2: i32) -> f32 {
    if s2 >= e1 || s1 >= e2 {
        return 0.0;
    }

    let s = cmp::max(s1, s2);
    let e = cmp::min(e1, e2);

    let ovlp_len = e - s;

    let min_len = cmp::min(e1 - s1, e2 - s2);

    if min_len == 0 {
        return 0.0;
    }

    ovlp_len as f32 / min_len as f32
}

#[cfg(test)]
mod tests {
    use std::time::Duration;

    use bio::alignment::pairwise::Scoring;
    use gskits::fastx_reader::read_fastx;
    use rust_htslib::bam::ext::BamRecordExtensions;

    use super::*;

    #[test]
    fn test_align_single_query_to_target() {
        let ref_file = "test_data/GCF_000005845.2_ASM584v2_genomic.fna";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);
        let aligners = build_aligner(
            "map-ont",
            &IndexParams::default(),
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
            10,
        );
        let target2idx = targets_to_targetsidx(&targets);

        let query_file = "test_data/ecoli_query.fa";
        let query_filter_iter =
            gskits::fastx_reader::fasta_reader::FastaFileReader::new(query_file.to_string());
        for qr in query_filter_iter {
            let hits = align_single_query_to_targets(&qr, &aligners, &MapParams::default());
            let records = hits2records(&hits, &qr, &target2idx);
            for record in records {
                assert_eq!(record.reference_start(), 720);
                assert_eq!(record.reference_end(), 1920)
            }
        }
    }

    #[test]
    fn test_targets_to_query_forward_targets() {
        let targets = vec![ReadInfo::new_fa_record(
            "t".to_string(),
            "ACGTAAACNRT".to_string(),
        )];

        let qf_targets = targets_to_query_forward_targets(&targets);

        assert_eq!(qf_targets.len(), 2);
        assert_eq!(qf_targets[0].name, "t___fwd");
        assert_eq!(qf_targets[0].seq, "ACGTAAACNRT");
        assert_eq!(qf_targets[1].name, "t___rev");
        // R 不是 ACGTN, reverse_complement 会补成 NUL(CStr 不允许), 这里替换成了 N
        assert_eq!(qf_targets[1].seq, "ANNGTTTACGT");
    }

    #[test]
    fn test_align_single_query_to_targets_query_forward() {
        let targets = read_fastx(FastaFileReader::new("test_data/ch11-6000.fa".to_string()));
        let qf_targets = targets_to_query_forward_targets(&targets);
        assert_eq!(qf_targets.len(), targets.len() * 2);

        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &qf_targets,
            10,
        );
        let target2idx = targets_to_targetsidx(&qf_targets);
        assert_eq!(aligners.len(), qf_targets.len());

        let queries = read_fastx(FastaFileReader::new(
            "test_data/query_of_ch11.fa".to_string(),
        ));
        assert!(!queries.is_empty());

        for query in queries {
            let rev_query = ReadInfo {
                name: query.name.clone(),
                seq: String::from_utf8(reverse_complement(query.seq.as_bytes())).unwrap(),
                ..Default::default()
            };

            // 同一条 read 的正反两条链, query_forward 模式下结果应当等价
            let hits = align_single_query_to_targets_query_forward(
                &query,
                &aligners,
                &MapParams::default(),
            );
            let rev_hits = align_single_query_to_targets_query_forward(
                &rev_query,
                &aligners,
                &MapParams::default(),
            );
            assert!(!hits.is_empty() && !rev_hits.is_empty());

            let primary = primary_hit(&hits).unwrap();
            let rev_primary = primary_hit(&rev_hits).unwrap();

            // primary 始终落在 query forward 链上, 即 read 不会被 reverse complement;
            // 被选中的是 target 的哪一份链, 取决于 read 自身的方向, 两条 read 必然互补
            assert_eq!(primary.strand, minimap2::Strand::Forward);
            assert_eq!(rev_primary.strand, minimap2::Strand::Forward);

            let (base_name, primary_on_fwd_target) =
                split_target_strand_name(primary.target_name.as_ref().unwrap().as_str());
            let (rev_base_name, rev_primary_on_fwd_target) =
                split_target_strand_name(rev_primary.target_name.as_ref().unwrap().as_str());
            assert_eq!(base_name, rev_base_name);
            assert_ne!(primary_on_fwd_target, rev_primary_on_fwd_target);

            // 两份 target 上的 hit 互为镜像: fwd 的 [s, e) <-> rev 的 [L - e, L - s)
            assert_eq!(primary.target_len, rev_primary.target_len);
            let target_len = primary.target_len;
            let (fwd_hit, rev_hit) = if primary_on_fwd_target {
                (primary, rev_primary)
            } else {
                (rev_primary, primary)
            };
            assert_eq!(fwd_hit.target_start, target_len - rev_hit.target_end);
            assert_eq!(fwd_hit.target_end, target_len - rev_hit.target_start);
            assert_eq!(alignment_score(primary), alignment_score(rev_primary));

            // record 中 query 始终没有被 reverse complement
            for (hits, query_record) in [(&hits, &query), (&rev_hits, &rev_query)] {
                let records = hits2records(hits, query_record, &target2idx);
                assert!(!records.is_empty());
                for record in records {
                    assert!(!record.is_reverse());
                    assert_eq!(
                        record.seq().as_bytes().to_vec(),
                        query_record.seq.as_bytes().to_vec()
                    );
                }
            }
        }
    }

    #[test]
    fn test_align_single_query_to_targets_query_forward_no_hit() {
        // 一条完全比对不上的 read, 两侧都没有 primary, 结果被过滤
        let targets = vec![ReadInfo::new_fa_record("t".to_string(), "A".repeat(2000))];
        let qf_targets = targets_to_query_forward_targets(&targets);

        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);
        let aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &qf_targets,
            10,
        );

        let query = ReadInfo::new_fa_record(
            "q".to_string(),
            format!("{}{}", "CGAT".repeat(500), "GCTA".repeat(500)),
        );
        assert!(align_single_query_to_targets_query_forward(
            &query,
            &aligners,
            &MapParams::default()
        )
        .is_empty());
    }

    /// 一个 30A 的同聚物 run, 前后缀都不含 AA(前缀以 T 结尾, 后缀以 G 开头),
    /// 因此 run 的边界是唯一的
    fn poly_gap_left_align_target() -> String {
        format!(
            "{}{}{}",
            "CGAT".repeat(15),
            "A".repeat(30),
            "GCTA".repeat(15)
        )
    }

    /// 与 `poly_gap_left_align_target` 等长, 但碱基完全不同
    fn poly_gap_left_align_decoy_target() -> String {
        format!(
            "{}{}{}",
            "GCTA".repeat(15),
            "C".repeat(30),
            "CGAT".repeat(15)
        )
    }

    fn poly_gap_left_align_targets() -> Vec<ReadInfo> {
        vec![ReadInfo::new_fa_record(
            "t".to_string(),
            poly_gap_left_align_target(),
        )]
    }

    /// 上述 target 删掉 run 最左侧 5 个 A 之后的 read
    fn poly_gap_left_align_query() -> ReadInfo {
        ReadInfo::new_fa_record(
            "q".to_string(),
            format!(
                "{}{}{}",
                "CGAT".repeat(15),
                "A".repeat(25),
                "GCTA".repeat(15)
            ),
        )
    }

    fn build_poly_gap_left_align_aligners(targets: &Vec<ReadInfo>) -> Vec<NoMemLeakAligner> {
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            targets,
            10,
        )
    }

    fn primary_cigar(hit: &Mapping) -> Vec<(u32, u8)> {
        hit.alignment
            .as_ref()
            .unwrap()
            .cigar
            .as_ref()
            .unwrap()
            .clone()
    }

    /// --polyNGapLeftAlign: 5A 的删除被摆到同聚物 run 的最左侧, 即紧跟 60bp 的前缀
    const LEFT_ALIGNED_CIGAR: [(u32, u8); 3] = [(60, 7), (5, 2), (85, 7)];

    #[test]
    fn test_align_single_query_to_targets_poly_gap_left_align() {
        let targets = vec![
            // 与真 target 等长但碱基完全不同: 一旦 gap 左移用错了 target 序列,
            // 重建出的 cigar 会在前后缀上出现大量 mismatch(op 8)
            ReadInfo::new_fa_record("decoy".to_string(), poly_gap_left_align_decoy_target()),
            ReadInfo::new_fa_record("t".to_string(), poly_gap_left_align_target()),
        ];
        let aligners = build_poly_gap_left_align_aligners(&targets);
        let query = poly_gap_left_align_query();

        let adjusted_hits = align_single_query_to_targets(
            &query,
            &aligners,
            &MapParams::default().set_poly_n_gap_left_align(true),
        );
        let adjusted = adjusted_hits.iter().find(|hit| hit.is_primary).unwrap();
        assert_eq!(adjusted.target_name.as_ref().unwrap().as_str(), "t");
        assert_eq!(adjusted.target_start, 0);
        assert_eq!(adjusted.query_start, 0);
        assert_eq!(primary_cigar(adjusted), LEFT_ALIGNED_CIGAR);

        // 开关只改写 gap 的位置, 不影响 primary/supp 判定与比对分数。
        // minimap2 在这组数据上本来就输出左对齐的 gap, 因此两侧 cigar 一致
        // (非左对齐的输入由 cigar_adjust 的单测覆盖)
        let plain_hits = align_single_query_to_targets(&query, &aligners, &MapParams::default());
        assert_eq!(plain_hits.len(), adjusted_hits.len());
        for (plain, adjusted) in plain_hits.iter().zip(adjusted_hits.iter()) {
            assert_eq!(plain.is_primary, adjusted.is_primary);
            assert_eq!(plain.is_supplementary, adjusted.is_supplementary);
            assert_eq!(alignment_score(plain), alignment_score(adjusted));
            assert_eq!(primary_cigar(plain), primary_cigar(adjusted));
        }
    }

    #[test]
    fn test_align_single_query_to_targets_query_forward_poly_gap_left_align() {
        let qf_targets = targets_to_query_forward_targets(&poly_gap_left_align_targets());
        let aligners = build_poly_gap_left_align_aligners(&qf_targets);
        let query = poly_gap_left_align_query();

        let hits = align_single_query_to_targets_query_forward(
            &query,
            &aligners,
            &MapParams::default().set_poly_n_gap_left_align(true),
        );
        assert_eq!(hits.len(), 1);
        let hit = hits.first().unwrap();
        assert_eq!(hit.strand, minimap2::Strand::Forward);
        assert_eq!(hit.target_name.as_ref().unwrap().as_str(), "t___fwd");
        assert_eq!(hit.target_start, 0);
        assert_eq!(primary_cigar(hit), LEFT_ALIGNED_CIGAR);
    }

    /// HP/TR 比对区域收紧(--hpTrShrinkAlnRegion)端到端用例的 target。
    /// 上面只有两处重复: `(AC)6` 占 [60,72), 21A 同聚物占 [102,123); 其余碱基是刻意挑出来
    /// 的非重复序列(不能像 `poly_gap_left_align_target` 那样用 `"CGAT".repeat(n)`,
    /// 那本身就是一条 (CGAT)n 串联重复, 会被整个当成 motif)
    fn hp_tr_shrink_target() -> String {
        [
            "CCGTAATGCCTTCCTAACAGAGTTCGAACTCGTGTTGTCGAGCGACGGAATTAGATCAGT",
            "ACACACACACAC",
            "TAATGGCAGAACTGGCAGGCTTAGTCGTGG",
            "AAAAAAAAAAAAAAAAAAAA",
            "AGATGATCAGTGGTAAGGTGGCGCGGTAACGCGCCTAAGGC",
        ]
        .concat()
    }

    const HP_TR_STR: (usize, usize) = (60, 72);
    /// 同聚物实际是 4th 段的 20 个 A 加上后缀首字符的 A, 共 21 个
    const HP_TR_HP: (usize, usize) = (102, 123);

    fn build_hp_tr_shrink_aligners(
        map_params: &MapParams,
    ) -> (Vec<NoMemLeakAligner>, Vec<ReadInfo>) {
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let targets = vec![ReadInfo::new_fa_record(
            "t".to_string(),
            hp_tr_shrink_target(),
        )];
        let targets: Vec<ReadInfo> = if map_params.query_forward {
            targets_to_query_forward_targets(&targets)
        } else {
            targets.into_iter().collect()
        };

        let aligners = build_aligner(
            "map-ont",
            &index_params,
            map_params,
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
            10,
        );
        // motif 表只在开关打开时构建, 关闭时不付扫描代价
        assert_eq!(
            aligners.iter().all(|a| a.hp_tr_regions.is_some()),
            map_params.hp_tr_shrink_aln_region
        );

        (aligners, targets)
    }

    fn hp_tr_shrink_align(query: &ReadInfo, shrink: bool) -> Option<Mapping> {
        let map_params = MapParams::default().set_hp_tr_shrink_aln_region(shrink);
        let (aligners, _) = build_hp_tr_shrink_aligners(&map_params);
        align_single_query_to_targets(query, &aligners, &map_params)
            .into_iter()
            .find(|hit| hit.is_primary)
    }

    #[test]
    fn test_align_single_query_to_targets_hp_tr_shrink() {
        let target = hp_tr_shrink_target();
        assert_eq!(target.len(), 163);
        // fixture 上只有这两处会被扫成 motif
        assert_eq!(&target[HP_TR_STR.0..HP_TR_STR.1], "ACACACACACAC");
        assert_eq!(&target[HP_TR_HP.0..HP_TR_HP.1], "A".repeat(21));

        // 比对末端停在 (AC)6 内部: 拉到重复区左端之外, 多出来的 read 碱基变 soft-clip
        let query = ReadInfo::new_fa_record("q".to_string(), target[..70].to_string());
        let plain = hp_tr_shrink_align(&query, false).expect("no hit");
        assert_eq!((plain.target_start, plain.target_end), (0, 70));
        assert_eq!(primary_cigar(&plain), vec![(70, 7)]);

        let shrunk = hp_tr_shrink_align(&query, true).expect("no hit");
        assert_eq!((shrunk.target_start, shrunk.target_end), (0, 60));
        assert_eq!((shrunk.query_start, shrunk.query_end), (0, 60));
        assert_eq!(primary_cigar(&shrunk), vec![(60, 7)]);
        // cs / md 跟着裁剪后的列重建, 不再描述被裁掉的 10 个碱基
        let aln = shrunk.alignment.as_ref().unwrap();
        assert_eq!(aln.cs.as_deref(), Some(":60"));
        assert_eq!(aln.md.as_deref(), Some("60"));

        // 两端都落在重复区里: 只留下中间那段干净的比对
        let query = ReadInfo::new_fa_record("q".to_string(), target[64..120].to_string());
        let plain = hp_tr_shrink_align(&query, false).expect("no hit");
        assert_eq!((plain.target_start, plain.target_end), (64, 123));
        let shrunk = hp_tr_shrink_align(&query, true).expect("no hit");
        assert_eq!(
            (shrunk.target_start, shrunk.target_end),
            (HP_TR_STR.1 as i32, HP_TR_HP.0 as i32)
        );
        assert_eq!((shrunk.query_start, shrunk.query_end), (8, 38));
        assert_eq!(primary_cigar(&shrunk), vec![(30, 7)]);

        // 比对始端停在同聚物内部
        let query = ReadInfo::new_fa_record("q".to_string(), target[110..162].to_string());
        let shrunk = hp_tr_shrink_align(&query, true).expect("no hit");
        assert_eq!((shrunk.target_start, shrunk.target_end), (123, 162));
        assert_eq!((shrunk.query_start, shrunk.query_end), (13, 52));
        assert_eq!(primary_cigar(&shrunk), vec![(39, 7)]);

        // 两端都干净(重复区完整落在比对区间中间): 开关不产生任何改动
        let query = ReadInfo::new_fa_record("q".to_string(), target.clone());
        let plain = hp_tr_shrink_align(&query, false).expect("no hit");
        let shrunk = hp_tr_shrink_align(&query, true).expect("no hit");
        assert_eq!(
            (
                plain.target_start,
                plain.target_end,
                plain.query_start,
                plain.query_end
            ),
            (0, 163, 0, 163)
        );
        assert_eq!(
            (
                shrunk.target_start,
                shrunk.target_end,
                shrunk.query_start,
                shrunk.query_end
            ),
            (
                plain.target_start,
                plain.target_end,
                plain.query_start,
                plain.query_end
            )
        );
        assert_eq!(primary_cigar(&shrunk), primary_cigar(&plain));
    }

    /// 收紧后的 hit 落到 record 上: 被裁掉的 read 碱基必须是 soft-clip, 且 pos/reference_end
    /// 与 cigar 覆盖的 ref 跨度自洽
    #[test]
    fn test_hp_tr_shrink_record() {
        let target = hp_tr_shrink_target();
        let query = ReadInfo::new_fa_record("q".to_string(), target[..70].to_string());

        let map_params = MapParams::default().set_hp_tr_shrink_aln_region(true);
        let (aligners, targets) = build_hp_tr_shrink_aligners(&map_params);
        let hits = align_single_query_to_targets(&query, &aligners, &map_params);
        let target2idx = targets_to_targetsidx(&targets);
        let records = hits2records(&hits, &query, &target2idx);
        assert_eq!(records.len(), 1);

        let record = &records[0];
        assert_eq!(record.reference_start(), 0);
        assert_eq!(record.reference_end(), HP_TR_STR.0 as i64);
        // 被裁掉的 10 个 read 碱基落在 soft-clip 里, cigar 覆盖的 ref 跨度 = ref_end - pos
        assert_eq!(
            record.cigar().iter().copied().collect::<Vec<_>>(),
            vec![Cigar::Equal(60), Cigar::SoftClip(10)]
        );
        assert_eq!(record.seq_len() as usize, 70);
    }

    /// --query-forward 与 --hpTrShrinkAlnRegion 组合: motif 表建在每个 aligner 自己索引的
    /// 那条序列上, `___rev` 链上的坐标天然一致, 不需要额外换算
    #[test]
    fn test_align_single_query_to_targets_query_forward_hp_tr_shrink() {
        let target = hp_tr_shrink_target();

        let map_params = MapParams::default()
            .set_query_forward(true)
            .set_hp_tr_shrink_aln_region(true);
        let (aligners, targets) = build_hp_tr_shrink_aligners(&map_params);
        assert_eq!(aligners.len(), targets.len());

        // 原始方向的 read, 两端都落在重复区里
        let query = ReadInfo::new_fa_record("q".to_string(), target[64..120].to_string());
        let hits = align_single_query_to_targets_query_forward(&query, &aligners, &map_params);
        let hit = hits.iter().find(|hit| hit.is_primary).expect("no hit");
        assert_eq!(hit.strand, minimap2::Strand::Forward);
        assert_eq!(hit.target_name.as_ref().unwrap().as_str(), "t___fwd");
        assert_eq!((hit.target_start, hit.target_end), (72, 102));
        assert_eq!((hit.query_start, hit.query_end), (8, 38));

        // read 反向时命中的是 `___rev` 链, 收紧后的区间是同一条 target 的镜像:
        // fwd 的 [s, e) <-> rev 的 [L - e, L - s)
        let rev_query = ReadInfo {
            name: "q".to_string(),
            seq: String::from_utf8(reverse_complement(query.seq.as_bytes())).unwrap(),
            ..Default::default()
        };
        let hits = align_single_query_to_targets_query_forward(&rev_query, &aligners, &map_params);
        let hit = hits.iter().find(|hit| hit.is_primary).expect("no hit");
        assert_eq!(hit.target_name.as_ref().unwrap().as_str(), "t___rev");
        assert_eq!(
            (hit.target_start, hit.target_end),
            ((target.len() - 102) as i32, (target.len() - 72) as i32)
        );

        // record 里 query 始终没有被 reverse complement
        let target2idx = targets_to_targetsidx(&targets);
        for (hits, query_record) in [
            (
                align_single_query_to_targets_query_forward(&query, &aligners, &map_params),
                &query,
            ),
            (
                align_single_query_to_targets_query_forward(&rev_query, &aligners, &map_params),
                &rev_query,
            ),
        ] {
            let records = hits2records(&hits, query_record, &target2idx);
            assert!(!records.is_empty());
            for record in records {
                assert!(!record.is_reverse());
                assert_eq!(
                    record.seq().as_bytes().to_vec(),
                    query_record.seq.as_bytes().to_vec()
                );
            }
        }
    }

    #[test]
    fn test_set_primary_alignment() {
        let mut r1 = BamRecord::new();
        r1.unset_secondary();
        let mut cigar_str = CigarString(vec![]);
        cigar_str.push(Cigar::Equal(3));
        cigar_str.push(Cigar::Diff(1));
        r1.set(b"1", Some(&cigar_str), b"AACG", &vec![255; 4]);

        let mut r2 = BamRecord::new();
        r2.unset_secondary();
        let mut cigar_str = CigarString(vec![]);
        cigar_str.push(Cigar::Equal(4));
        r2.set(b"2", Some(&cigar_str), b"AACT", &vec![255; 4]);

        let mut r3 = BamRecord::new();
        r3.set_secondary();
        let mut cigar_str = CigarString(vec![]);
        cigar_str.push(Cigar::Equal(2));
        cigar_str.push(Cigar::Diff(2));
        r3.set(b"3", Some(&cigar_str), b"AAGC", &vec![255; 4]);

        let mut records = vec![r1, r2, r3];

        // set_primary_alignment(&mut records);

        // assert_eq!(records[0].is_secondary(), true);
        // assert_eq!(records[1].is_secondary(), false);
        // assert_eq!(records[2].is_secondary(), true);
    }

    #[test]
    fn build_aligner_memory_leak() {
        for _ in 0..10000 {
            thread::sleep(Duration::from_millis(1));
            let aligner = Aligner::builder().map_ont();
            let aligner = aligner
            .with_index_threads(1)
            .with_cigar()
            .with_sam_out()
            .with_sam_hit_only()
            .with_seq_and_id(b"ACGGTAGAGAGGAAGAAGAAGGAATAGCGGACTTGTGTATTTTATCGTCATTCGTGGTTATCATATAGTTTATTGATTTGAAGACTACGTAAGTAATTTGAGGACTGATTAAAATTTTCTTTTTTAGCTTAGAGTCAATTAAAGAGGGCAAAATTTTCTCAAAAGACCATGGTGCATATGACGATAGCTTTAGTAGTATGGATTGGGCTCTTCTTTCATGGATGTTATTCAGAAGGAGTGATATATCGAGGTGTTTGAAACACCAGCGACACCAGAAGGCTGTGGATGTTAAATCGTAGAACCTATAGACGAGTTCTAAAATATACTTTGGGGTTTTCAGCGATGCAAAA",  b"ref")
            .unwrap();
            let aligner = NoMemLeakAligner::new(aligner, &[]);
        }
    }

    #[test]
    fn test_align_single_query_to_target2() {
        let ref_file = "test_data/MG1655.fa";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
            10,
        );

        let aligner = &mut aligners[0];

        // aligner.mapopt.best_n = 100;
        // aligner.mapopt.pri_ratio = 0.1; // min dp score

        // tot len
        let seq = b"AAAAGAGAGAGATAGGGGGCGACGAGCGCCAAACGGGCGGCAATTATAGGGATTTCATCGCCTGATACCAGTCGAATAGCGTTGCCGCGCGCTCAGGATGTTAATTGGTTGACAGAAGAATTCCCGGTGGCAAAATTACGTTGAAGATCAGTTTTATACAAGGTAAAAAATGTTATACGCAGTTGCGCAATTATCGCCTTTACGTCACTTTATGAGCATTCGCATATAAAATGTAAAACTTTTTGTAACTAGCATAAAACACAGAAACGAATACCTGGCGCTGGTCTTGCGATAAAGCAGGTAATGAGCAAACAACACAGCATGTATTAATTGCCCTGCCCACCCGCTGCTTCCACCTGGTCCAGTTTAAGGTTAGTCCTCGTTTACTTTACCCTTTTCTCGCTGAGCTTTCGCAAGTTTGGCACCCTCTCGCCCCACGCCTGTGGTTCCGACGGTCCACTGATGGTGGCGTTTATCGCCATGCCGGGGCGCAATGTGCCGGGAAATGCGCTGAGCTGTTCGCTGGAAATATCGCCGCATCCATCCTGCTTTTTTCCACCAGCCGCTGAACATGACCTGGACGACATCAATATTGTTGAAGCCGTGGTCGGGCAGTGCATCTCTCGCTCAAACAACAACGGAGGAGGAGAAAAAAGAGACAGATGCACTGCCCCGACCACGGCTTCAACACATATTGATGGTCGTCCAGGTCATTTCCACGAGTGGTGGAACAAAGCAGGAGGATGCGGCGTATTTCCAGGAACAGCTCACGCAATTCCCGGCACATGGCGCCCGGCAATGGCGAATAAAACGCCAACCATCATGAATGGACGTCGGAACCACATGGGCGATTGGTGCCAAACTGCGAAAGCTCAAGCGAGAAAAGGGTAAAGAAAAACGAGACTTATACCTAAATGACCAGGGAATGCAGCGGGTGGGGCAGGGCAATTAATTACATGGCTGTGATTGTTGCTCATACCGCTTTATCCCAAGACAGGTCGCCAGTATTCGTTTCTGTGTTTAATGCTATACAAAAGTTTATTTACATTTATATGCGAAATTGCTCATAAGTGAGGGTAAAGGCGATAATTTGCGCAACTGCGTAATAACATTTTTTTACCTTGTATAAAACTGATCAACGTAAAGTTGCCACCGGGATCTTCGTCAAAATTAACTCTGACGCGCGGCAACGCTTATTCGACTGGTATCAGGCGGATGAAATCCTATAAATTGCCGCGTTTGGCGCTTCCGTCGCCCCCTATCTCTCTCAAACAACAACGGAGAGGAGGAAAAAAGAGAGAAGATAGGGGGCACAGAAGCGCCCAACGCGCAATTAATAGGGATTCATCCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCATGATACCAGTCGAATAGCGTTGCGCGCGCTCAGAGTTACATTGTTGACGAAGAATTCCCGGGTGGCAAATTACGTTGATCAGTTTTATAACAAGGTAAAAAGTTATACGCAGTTGCGCAAATTATCCGCCTTACGTCACTTTAATGAGCAATTCCGCATAATAAAAATGTAAAACTTTGTAACTAGCATAAACACAGAAACGAATACTGGGGCGACCTGGTCTTGCGGAATAAGCGGTAATGAGCAAAAATCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGCTTCACCGGTCAGCGTTTAGGTTTAAGTCTCGTTTAATCTTTACCCTTTCTCGTTGAGCGAGCTTTCGAGTTTGCACCCAACTCGTCCCACTGTGGTTCCCGACGTCCATCATGTGGGCGTTTTAAAATCGCGCATGCCGGCGCATGTGGCCGGGAATTGCCTGAGCTGTTACGCTGGGAAATATCGCCGCATCCATCCTGGCTTTTTTTCCACAGCTACGCTGAACAGACCTGGACGACCACCATCAATATGTTGAAGCCGTGGTCGGGCGTGCATCTCTCTCAAACAACAACGGAGGAGGGAGGAAAAAAGAGAGAGATGCACTGCCCCGACCACGGTTTTTTGTCAACATATTGAAGGTCGTCCAGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAATCCCGGCCACAGCGCCCGGCAATGGCGATTAAAACGCCACCATCATGATGGAACGTCGGGAACCACAGTGGGGCGTGTTGGGTGCAAACTCGAAAAGCTCAAGCAGAAAAGGGTAAAGATAACACGAGACAAACCAAACTGACCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGTTTATCCGCAAAACAGGTCGCCAGTATTCGTTTCCTGTGTTTATGCTAGTACAAAAGTTTTACATTTTAATATGCGAATTGCTCATAAAGTGACGTAAGGCGGAAATTTGCGCAACTGCGTAATAACATTTTTTTAACTTGTATAAAACTGTTCAACGTAATTTGCCACGGAATTCTTCGTCAACAATTAACTACCTGAGCGCGCCAACGCTATTCGACTGGTATCAGGCGAATGAATCCCTATAATTGCCGCGTTTGGCGCTTCGTCCGCCCCCTATCTCTCTCAAACAACAACAGGAGGGGAGGAAAAAAGAGAGAGATAGGGTGGCGACGAAGCGCCAAACGCGGCAATTATAGGGATTTCATCCGCCTGATACCAGTCGAATAGCGTTGCCGCGCGCTCAGAGTTAATTGTTTGACGAAGAATTCCCGGTGGCAAATTACGTTGATCAGTTTTTATACAAGGTAAAAAAATGTTATACGCAGTTGCGCAAATTATCCGCTTTACAGTCCACTTTTAGAGCAATTCGCATAATAAAAGGTAAAACTTTTGTACTAGCATAAGAACACACACAGAAACAATACTGGCTGACCTGGTCTTGCGGATAAAGCGGTAATGAGCAAACATCACGCATGTATTAATTTGCCCCTCCCCACCGCTGCTTCAACCGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTAGTACCCTTTTCTCGCTTGAGCTTTGCAGTTTGGCACCCAACTCCCATGTGTTTCCGCGACGTCCATCATGATGGTGGCGTTTTATCGCATGCCGGGCGCATGTGGCGGGAATTGCGCTGAGCTGTTCGCTGGGACATATCGCCGGCATCCTCCCTGCTTTTTTCCCCAGCTCGCGCTGAACATGACCTGGACGACCAAGTCAATATTGTTGAGCCGTGGTCGGGTGCAAGTGCATCTCCTCAACAACAACGAGGAGGGGAAAAAAAAAGATGAGAGATGCACTGCCCCCGACCACGGCTTCACACAATATTGAGGGTCAGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGATGGAGCGCGATATTTCCCACGACAGCTCAGCGAAACAAAAAAAAAAAATTCCGCCACATGCGCCCGGCATGGCGATAAAACGCCACCTCCAGATGGACGCGGGAAACCACAGTGGGGCGAGTTGGGTGCCAAACTGCGAAAGTCAAGCAGAAAAGGGTAAAGAAAACGAGCTAAACTAACTGCACAGGTGAAGCAGCGGGTGGGCAGGCATAATACATGTGTGATTGTTTGGTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCATTACCGCTTTATCCCAAGACCAAGGTCGCCAGTATTCGTTTCTGGGTTTATGCTAGTACAAAAGTTTACATTTTATATGCGAATTGCTCCATAAAGTGGCCCCCCCCCCCCGTAAAGGCGGATAATTTGCGCAACTGCGTATAACATTTTTTACCTTGTATAAAACTGATCAAGTAATTTGCCACCGGAATTCTTTCGTCAACAATAACTCTGGCGCGCGGCAACGCTATTCGACTGGTATCAGGCGGAGAAATCCTAATAATTGCCGCGTTTGGCGCTTCGTCGCCCCTATCTCTCTCCAAACAAACAACGAGGAGGAGGAAAAAGAGAGAGATAGGGGCGACGAAGCGCCAAACGCGGCAAATTATAGGGATTCATCCGCCTGATAACCAGGTCGAATAGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCAGCGTTGCCCGCGCTCAATTTAATGTTGACGAAGAATTCCCGGTGGCAAATTAACACGTTGAATCAGTTTTATACAAGGTAAAAAATTTATAGCAGTTGCGCAATTAATCCGCCTTACGTCACTTATGAGCAATTCGCATATAAAATGCTAAAACTTTTTGTACTAGCATAAAACACAGAAACGAATACTGGCGACCTGGTCTCGCGGATAAAGCGGTAGATGAGCAAAACAATCACAGCACTGTATTAATGCCTGCCCCACCGCTGCTTCACCGAGTCAGTTTAGGTTTATAGTCTCGGTTATCTTTACCCTTTTCGCTTGAGCTTTCGCAGTTTGGCACCACACTGCCCATGTGGTTCCGCGTCCATCCATGATGGTGGCGTTTTTATCGCCATGCCGGGAGCGCATGTGGCCGGAATTGCGCTGAGCTGTTCGCTGGAAATATCGCCGCATCCGATCCTGCTTTTTTCCACCAGCTCGCTGAACATGACCGGACGACCATCAATATTGTGAAGCCGTGGTCGGGGCAGTGCATCTCTCTCAACAAACAAACAACAGGAGGAGAGAAAAAAGAGAAAGCATGCCCCGACCACGCTTCAACAATATTGAGTGGTCGTCCGGTCATGGTTCAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGCGATATTTCCCAGCGGAACACTCAGCGCAATTCCGGCCAAACAGCGCCCGGCATGGCTGATAAAACGCCACCATCATGATGGACGTGGGAACCACAAGTGGGGCGATTGGTGCCAACTGCGAAAGCTCAAGCGAGAAAGGGTAAAAGATAAACGATGACCTAAACCTAAACTGACCAAGTGAACAGCGGGTGGGGCAGGCAATTAATACATGTGTGATTGTTTGCCATTACCGCTTTATCCCAAGACCAGGTCGCCAGTATTCGATTTCTGTGTTTATGCTAGTACTAAAGTTTTACATTTATATGCAGAATTGCTCATAAAGTGAGTAAGGTCGGATAATTTGCGCACACTGCGTATAACATTTTTTACCTTGTATAAAACGATCCAACGTAATTTGCCACCGGGAATCTCGTCAACAATTAACTCTGAGCGCGCGGCAACGCTATTCGACTGTATCAGGCGGATGAAATCCCTATAATTGGCCGCGTTTGGCGCTTCGTCGCCCCCATCTCTCTCAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATAGGGGGCGACGAAGCGCCAAACGCGGCATATAGGGATTTCATCCGCCTGATACCAGTCCGAAAGCGTTGCCGCGCGCTCAGAGTAAGTTGACGAAGAATTCCCGTGTGGGCACAATAACGTGACAGTTTTATACAAGGTAAAAAATGTTATACAGCAGTTGCGCAAATTTCCCCTTTACGTCACTCACTTTATGAGCAATTCGCATATAAAATGTAAACCCTTTTGTACTAAGCATAAACAGCAGAAACGAATACCTGCGATACCTGGTCTTGCGGATAAACGGTAATGAGCACAACTCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGTCTTCAACCTGGTCATGTTTAGGTTTATCTCGTTTATCTTACCTTTCTCGCTTGAGCTTTCGCATATTGGCACCCAATCGCCCCACGTGGTTCCCGACGTCCATCATGAGGTGGGCGTTTTATCCCATGCCGGGCGGCAATGTGGCCGGAATTGCGCTGAGCTGTTCCGCTGGGAAAATCGCCGCATCCATCTAGCTTTTTTTCCACAGGCTCGCTGAACATGATACCTGGACGACCATCAATATTGTTGAAGCGGGTCGGCAGTGCATCTTCTCAAACAACAACCGTAGGAGGAGGAAAAAAGAGAGAGATGCACTGCCCCGACCAACGGCTTCACAATATTGATGGTCGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGTCGGCGATATTTCCCAGCGAACAGCTCAGCGCAATTCCCGGCCACATGCGCCCGGCATGCGATATAAAGCCACATCATGATGGACGTCGGAACCACATGGAGCGAGTTGGGTGCCAAACTGTGGAAAGCTCAGCGACAAAAGGGTAAAGATAAACGAGACTAAACCTAAACTGACCAGTGAAGCAGCGGGTGGGCAGGGCAATTAATACATGCTGTGATTGGTTTGCTCATACCGCTTACGCAAGACCAGGCCGCAGTATTCTTTCTGTGGTTTATGCTAGTACAAAAGTTTTACATTTTATATGCGAATTGCTCAAAAGTGACGTAAAGCGGATAATTTGCGCAACTCGTATAACATTTTTACCTTGGTATAAAACTGATCAACGTAATTGCCACCGGGAATCTTCGTCACAATTAAGACCTGAGCGCGCGGCAACGTATTCGACTGGTATCCAGCGGATGCAAATCCCTATAATTGCCGCGTTGGCGCTTCGTCGCCCCATTCTCTCAAACAACAACGGAAGGAGAGGAACAAAGAGAGAGCATAGGGGGGCGGGGGCGAACGAAGGCAAACAGCGGCAATTATAGAGGATTTCACTCCGCCTATACCAGTCGAATAGCGTTGCCGCGCGCCAGAGTTAATTGTGACGAAGAATTCCCGTGGCAAATTAACGATGATCAGTTTTATACAGGGTAAAAAATGTTATTACAGCAGTTGCGCAAATTATCCGCCTTACGTCACTTTATGCAGCAATTCGCATATAAAATGTAAAACATTTTGTACTAGCATAAACACAGAAACGAATACTGGCGACCTGGTCTTGCTCGGATAAAGCGGTTAATGAGCAAACAATCACCAGCAGTATTAATTGCCCTGCCCCCACCGCTGCCTTCACCTGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTCTCGCTTTGAGCTTCGCAGTTTGGCACCCAACTCGCCCCCACTGTGGTTCCCGACTCCCATCATGATGGTGGCGTTTTATCGCCATGCGGGCGCATGTTGGCCGGGAATGCGCTGACTGTTCGCTGGGAATATCGCCGCATCATCCTGCTTTTTTCACCAGCTCGCTGAACATGACCTGGACGACATAATAGTTGAAGCGTGTCGGCAGTGCATCTCTCTCAAAACAACAACGAGGAGGAGGAAAAAAGAGAAGATGCCTGCCCCCGCACAGGTTAACAATATGAGGTCGTCCAGGTCGATGTTCAGCGGACTGGTGGGAAAAAAAAAAAAAGCAGGATGGAATGCGCATATTTCCCAACAGCTCAGCGCAATTCCCGGCCCACAATGCGCCCGGCATGGCGATAAAACGCCACCACAGAGGACGTCGGAACCACAGTGGGGCGAGTTGGGTGCCAAATCTCGAAAGCCAGCGAGAAAAGGGTAAAGATAAAACGAGACTAAACCTAAACTGACCAGGTGAAGCAGCCGGGTGGGGCAGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCCGGCTTTATCCCAGACCAGGTCGCCAGTATTCGTTTCTGTGTTTATGCTAGTACAAAAGTTTAACATTTTATAGCGATTGTTCATAAAGTGACGTAAAGGCGGATAATTTGCCAACTGCGTATAACGATTTTTACCTTGTATAAAACTGATCAACGTATTTGCACGGGAATTCTTCGTCCCCAATTGAACTCTGAGCGCGCGCACGCTATTCGACCTGGTATCAGGCGGATGAAATCCGTATAATTGCCCGTTTGCGCTTCGTCGCCCCTATCTCTCTCAAACAACAACGGAAGGAGAGGAAAAAAAGAAGAGATAGGGGGCGACGAAGCGCAAACGCGCAATTTATAGGGATTTCTCCAGCCTGATACCAGTCGACATAGCGTTGCCGCGCGCTCAGAGTTAATGTTGACGAAGAATTCCCGGGTGCAAATTACGTTGATCAGTTTTTTTTTATACAAGTAAAAAATGTTATACGCAGTGGCAAATTTCCGCCATTTACGTCACTTTATGAGCAATTCGCATAAAAAATGTAAAGTTGCATTGTAAAATCTTTTGTACAGCATAAACACAGAAACGAATAGCTGGCGAACCTGGTCTGCGGATAAGAGCGGTAATGAGCAACACACAAGCATGGTATAATTGCCCTGCCCCACCCGCTGCTCACCTGGTCGTTTAGGTTAGTTCGTTTACTTTAACCCTTTCTCGCTTGAGCTTCGCAGTTTTGGCACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGGGCGTTTTACGCCATGCCCGGGCGAGTGGCCGGGAATTGCGCTGAGCTGTTCGCTGGGAATATCGCCGCACCATCCTGCTTTTTTCCACCAGGCTCCGCTGGAACATGCCTGGCGACCAATCATATTGTTGAAGCCGTGGTCGGGGCAAGTGCATCTCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGCATGCCCGACCACGCTTAACAATATTGATGGTCGTCCAGGTCATGTTCAGCGGAGCTGGTGGAAAAAAGCAGGATTGGATGCGGCGATATTTCCCAGCGAACAGCTCAGCGCAATTCCGCCACATGCGCCGGCATGGCGATAAAACGCACCATCAATGTGGACGTCGCTGAACCGAAGTGGGGCGAGTTGGTGCCAAACTGCGGAAAGCTCAAGGAGAAAAGGGTAAAGATAAACGAACAAACTAAATGACCAGGTGAAGCACGGTGGGGCCAGGGAATTAATACAGCTGTGATTGTTTGTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTCATTACCGCTTTATCCGAAGACCAGGTCGCCATATTCGTTTCTGTGGTTTAGCTAGTACAAAGTTTTACATTTTATATGCGAATTGCTCATAAAGTGACGTAAAGGCGAATAATTTGCGCACTGCAGTATAACATTTTTTACTGTATAAAACTGATCAACGTAATTTGCCACCGGGATTCTTCGTCAACAATTTACTCTGAGCGCGCGGCACGCTATTCGACTGGTAATCAGGCGGATGAAATCCCTATAATTGTCCCCCCCGGGGGGCGTTGGCGCTCGCGCCCCCTATCCTTCTCAACAACAACGGAGAGAGGAAAAAGAGAGAGATAGGGGGGGGCGACAAGCGCCAAACGCGGCAATTATAGGGATTTCCATCCGACCTGATACCATCGAATAGCGTTTGCCGCGCGCTCAGAGTGTAATTGTTGACAGAATTCCCGGTGGCAAATTACGTTGATCAGTTTATACAAAGGTAAAAATGTTATACGCAGTTGCGCAAATTATCCGCCTTTAGTCCTTATGAGCAATTCGCATAAAACCTAAAATGTAAAACTTTAGTACTTGCAATAAACACAGAAACAATACTTGGCGACCTGGTCCCTCTGCGGATAAAGCGGTAAATGAGCAAACAACACAGCATGTATTAATTGCCCTGCCCCAACCCGCTGCTTCACTGGTCAGTTTAGGTTTAGTCTGTTTTATCTTTACCCACTTTTCTCGCTGAGCTTCCGCAGTTTGGCACCCAACTCGCCCACTGTGTTCCCGACGTCCATCATGATGGTGGCGTTTATCGCCATGCCGGGCGCATGTGGCCGGAATTTGCGCTGAGCTGTTTGCTGGGGAAATATCGCTGCATCCATCCTGCTTTTTTCCACGAAAAAAACAGCTCGCAGAAATGACCTGGACGACCATCATATTAGTTGAGCCGTGGTCGGGGCATGCATCTCTCTCAAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGCACTGCCCCGAACACGGCTTCAACAATATTGATGGTCGTCCAGGTCCCCATGTTCAAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGTGCGATATTTCCCATCCCCCCCCCGAGCGAACAGCTCAGGCAATTCCCGGTCCACATGCGCCCGGCATGGCGGATAAAACGCCACCCGATCATGATGGACGTTCCAGGGAACCACATGGGGCGAGTTGGTGCCAAACTGCGAAAGCTCAAGGCGAGAAAAAGGGTAAAGATAAACGGACAAACCTAAACTGTCCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGGCTGTATTGTTTGCCATTAGCCCCCCCCCCCCCCCCCATTCCGCTTTTATCCCGCAAGACCAGGTCGCCATATTCGTTCTGTGTTTAATCTAGTACACAAGTTTACATTTATAGCTAATGCGAATTGCTCATAAGTGACGTAAAGGCGGATAAATTCGCAACTGCGTATAACATTTTTTACTTGTATAAAATGATCAACGTAATTTGCCACCGGGAATCTTCGCAACAATAACTCTGAGCGCGCGGCACGCTATTCGACTGGTATCAGGCGGTGAAATCCCTAAATTGCGCGTTTGGCGCTTCGTCGGCCTTCCTCTCAAACAACACGAGGAGGAGGAAAAAAAAAAAAGAGAGAGATAGGGGGCGACGAAGCGCAAACGCGCAATTATAGGGATTTATCCGCCTGATACCAGTCGATAGCAGTTGCCGCGTCGCTCAGAGTTAATTGTTGACAAGAATTCCGCGGTGGCAAATTACGTGATCAGTTTTATACAAGAGTAAAAAAGTTTATACGCAGTTGCGCAAATTATCCGCCTTTACAGTCACTTATGAGCAATTCGCATATAAATGTAAAACTTTTGTACTAGCATAAACACAGAACACGAATACTGGCGACCTGGTCTTGCGGAAAAGCGGTAATGAGCAAACAAATCACAGCATGTATTATTGCCCTGCCCCACCCGCTGCTTCAAACCTGGTCAGTTTAGGTTTAGTCTCGTTTATCTTACCTTTTGCTCGTTGAGCTTTCCAGTTTGGCACCCAACTCGCCCACTGTGGTTCCCGACGTCCATCATGATGGTGGCTGTTTTATCGCCCATGCCGGGCGCATGTGGCCGGAATTGCGCTGAGCTGTTCGCTGGGAAACTATCGCCGCATCCATCCGCTTTTTTCCACCATGCTCGCTGAACATGACCTGGACGCACCATTCATATTGTTGAAGCCGTGGTGGGGCATGCATCTCTCTCAAACAACAACGGAGGAAGAGGAGAAAAAGGAGAGAGATGCACTGCCCGACCACGGCTTCAACAATATTGAGGTCGTCAGGTCATGTTCAGCGAGCTGGTGGAAAAAAGCAGGTGGCAGGATGGGCGGCCGATATTGCCGAGCGAACAGCTCAGCGCAAATTCCGGCACATGCGCCCGGCATGGCGAATAAAACGTGCACCATCATGATGGACAGCGGGAACCACAGTGGGGCGAGTTTTGGGTGGCCAACTGCGAAAGCTCAACGAGAAAAGGGTAAAGATAAACGAGACTAAACTAAACTGACAGGTGGAAGCAGCGGTGGGCAGGCAATTAATACAGCTGTGATTGTTTGCTCATTACCGCTTTTCCGCAAGACCAGGTCGCAAGTATTCGTTTCTGTGTTTATGCTAGTAACAAAAGTTTTACTTTTATATGCGAAAATTGCTCAAAAGTGACAGTAAAAAGGCGGATAATTTTGCGCACATGTCGTATAACATTTTTACCTTGTATAAAACTGACAACGTAAATTTGCCACCAGGGAATTCTTCGTCCAACAATTTAACTCTGAGCGCGCGGCAACGCTATTCGACTGTATCAGGCGGATGAAATCCCTATAATTGCCGGCTGATTGGCTGCTTCGTCCGCCCCCTACGTCCTCAAACAAAACGGAGGAGGAGGAAAAAAAGAGAGATAAGGGGCGACAAGCGCAAAACCGGCAAATTAAGGGATTTTCATCCGCCTGATACCAGTCAGGAATAGCGTTGCCGAGCGCTCAGAGTTAATTGTGACAAGAATTCCGGTGGCACAAATTATTGTCAGTTTATACAAGGTAAAAAATTTATACGCAGTTGCGCAAATTATCCGCCTATTACGTCACCTTTAGACAATTCGCATAAAAATGTAAAACTTTTGTACTAAGCATAAACACAGAAACGAATACTGGCGACTTGGTCTTGCGGATAAAGCGTAATGAGCAAACAATCACACATGTATTAATTGCCCTGCCCCAGCCGCTGCTTCACCAGGTCCAGTTTAGGTTTAGTCTCGTTATCTTTACCTTTTCTGCTTGAGCTTTCGCAGTTTGGCACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGTGCGTTTTATCGCATGCCGGGCGCATGGCCTGGAATTGCGCTGAGCTGTTCAGCTGGGAATATCGCCGCATCCACCTGCTTTTTTCCCACAGCTCGCTGAACAGACCTGGAACGACCATCATATTGTTAGAAAGCCGTGGTCGGGGCCGTGCATCTCTCTCAAACCAACAACGGAGGACAGAAAAAAAGAGAATGCACTGCCCGCCACAGGCTTCAACAATATTGATGGTCGTCCAGGTCATGTTTCAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGCGATATTTCCCCAGCGAACAGCTCAGCGCAATTCCCGGCCAATGCGCCCGGCATGGCGATAAAACGGCCAACCATCATGATGGACGTCGGGACCCATGGGGGCGAGTTGGGTGCCAAACTGCAAGCTCAAGCAGAAGGGTAAAGATAAGATAAACGAGACTAAACCTAAACTTGACCAATGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTAGCTTTAATCCGCAAGACCAGGTCGCCAGTATTCGTTTCTGTGTTTATGCTAGTACAAAAGTTTTATGACATTTTTATATGCGAATTGCTCAAAAGTGACGTAAAGGGGATAAATTTGCGCAACTGGCGTATAACATTTTTTACTTGGTATAAAGCCTAAAACTGATCAACCGTAATTTGCCACCGGGAATTCTTCGTCAACAATTAAACTCTGAGCGCGCGGCAACGCTATTCCGACTGTAATCAGCGATGAAATCCCTATAATTGCCCGCGTTTGCCGCTTCTCGCCCTATCCTCTCAAACAACAACGGAGGAGGGGAAAAAAGAGAGAGATAGGGGCGACGAAGGCCGCCAAACGCGGCAATTATAGGGATTTTCATCCGCCTGATACCAGTCGAATAGCGTTTGCCCGCGCGCTCAGATTAATTGTTGACGAAGAATTCCCGGTTGGCAAATTACTCTGATGCATGCTTTTAGACAAGGTAAAAAGATGTTATACGCAGTTGCGCAAATATCCGCCTTTACGTCACTTTATGAGCAATTCGCATCATAAAATGTAAAACTTTGTACTAGCATAAACAAGAAACGAATACTGGCGAACCTGGTCTGCGGATAAAGCGGTAATGACAAACAATCACAGCATGTATTAATTGCCCTGCCCCACCGCTGCTTCACCGGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTTCTCGCTTGAGCTATTCGCAGTTTGGCACCCTTCTCGCCCACAGAGAAGAGGTGGTTGCCCGACGTCCCATCATGATGGTGGCTTTTATCGCCATGCCGGGCGCATGGGCCGAATAGCGCTGAGCTGTTCGCTGGGGAATATCGCCCGCATCCAATCCTGCTTTTTTCCACGAGCTCGCTGAACATGACCTGGACGACCACAATATTGTTGAAGCCGTGGGTCGGGGGCCAGTGCGATCTCTCTCAACAAAAACGGAGGAGGAGGAGAAAAAAGAGAGGAGCACTGCCCCGACCACGGCTTCACAATATTTGATGGTCGTCCAGGTCATGTCAGCAGAGCTGGTGGAAAAAAGCGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAAATTCCCGGCCAACATGCGCCCGGCAATGCGATAAACCGCCACCAATCATGATGGACGTCGGGGACCACAGTGGGGCGAGTTGGGTGCCAAAACTGCGAAAGCTCAAGCGAGAAAAGGGTAAAGATAAACGAGACAAACCTAAACTGGACCAGGTGAAGCAGCGTGGACAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGCTTTATCCGCAAGACAGGTCGCCAGTAATCGTTTCTTGTGTTTAGCTAGTACAAAAGTTTTTATTTTATATGCGAATTGCTCATAATGTGACTGTAAAAGCGGATAATTTGCGCAACTGCGTATAACATTTTTTACCTTGTATGAAAAACTGATCAACGTAATTTGCCAACCGGGAATCTTCGTCAACAAATAACTCTGAAGCAGCGCGCAACGGCTATTCCGACTGGTATCAGGGCGGCATGAAATCCCTATAATTGCCGCGTTTGGCGCTCGCGCCCCCTATCTCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGCAGATAGGGGCGAGAAGCGCCAAACGCGGCAATTATAGGGATTTCATCCGCCTGATACCAGTGAAATAGCGTTGCGCGCGCTCAGAGTAATTTTGAGAAGAATTCCCGGTGGCAATTACGTGATCAGTTTTATACAAGGTAAAATGTTATACGCAGTTGCCAAATTATCCGCCTTTACGTCACTTTTATGGCAATTCGCATATACAAATGTAAAACTTTTGTACTAAGCATAAACACAGAAACGAATAGCTGGCCGACCTGGTCCTTTGCGGATAAGGGGAAGCGGTAAATGAGCAAACATCACAGCATGTATTAATTGGCCCTGCCCCCCGCTGCTTCACCGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTTCTCTGAGCTTTTCGCAGTTTGGACCCAAATACTCGCCCCACTTGTGGTTCCCGACGTCCATCGATGATGGTGGCGTTTTTATCGCCAATGCCGGGCGCATGTGGCCAGGGAAATTGCGCTGAGCTGTTCGCTGGGAATATCCGCCGCCAATCCTCCTGCTTTTTTTCCACCAGCTCGCTGAACATGACCTCGGACGCCATCAATATTGTTGAAAGCCGTGGTCGGGGCAGTGCATCTCTCTCAAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATGCTCTGCCCACCACGGGCTTCAAAATATGATGGTCCGTCCAGGTCATGTTCAGCGAGCTGGTG";
        // 2048 ~ 2648
        // let seq = b"GCACTGCCCCGACCACGGTTTTTTGTCAACATATTGAAGGTCGTCCAGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAATCCCGGCCACAGCGCCCGGCAATGGCGATTAAAACGCCACCATCATGATGGAACGTCGGGAACCACAGTGGGGCGTGTTGGGTGCAAACTCGAAAAGCTCAAGCAGAAAAGGGTAAAGATAACACGAGACAAACCAAACTGACCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGTTTATCCGCAAAACAGGTCGCCAGTATTCGTTTCCTGTGTTTATGCTAGTACAAAAGTTTTACATTTTAATATGCGAATTGCTCATAAAGTGACGTAAGGCGGAAATTTGCGCAACTGCGTAATAACATTTTTTTAACTTGTATAAAACTGTTCAACGTAATTTGCCACGGAATTCTTCGTCAACAATTAACTACCTGAGCGCGCCAACGCTATTCGACTGGTATCAGGCGAATGAATCCCTATAATTGCCGCGTTTGGCGCTTCGTCCGCCCCC";
        println!("{:?}", aligner.mapopt);
        println!("{:?}", aligner.idxopt);
        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score
            );
            println!("{:?}", MappingExt(hit).aligned_2_str(targets[0].seq.as_bytes(), seq));
        });

        // aligner.mapopt.best_n = 100;
        // aligner.mapopt.pri_ratio = 0.1; // min dp score
        // 3353~3662
        let seq = b"CACTGCCCCCGACCACGGCTTCACACAATATTGAGGGTCAGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGATGGAGCGCGATATTTCCCACGACAGCTCAGCGAAACAAAAAAAAAAAATTCCGCCACATGCGCCCGGCATGGCGATAAAACGCCACCTCCAGATGGACGCGGGAAACCACAGTGGGGCGAGTTGGGTGCCAAACTGCGAAAGTCAAGCAGAAAAGGGTAAAGAAAACGAGCTAAACTAACTGCACAGGTGAAGCAGCGGGTGGGCAGGCATAATACATGTGTGATTGTTTG";

        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score
            )
        });
    }

    #[test]
    fn test_align_single_query_to_target3() {
        let ref_file = "test_data/MG1655.fa";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
            10,
        );

        let aligner = &mut aligners[0];

        // aligner.mapopt.min_cnt = 2;
        aligner.mapopt.best_n = 500;
        // aligner.mapopt.min_dp_max = 10; // min dp score
        // aligner.mapopt.pri_ratio = 0.5; // min dp score
        aligner.idxopt;
        println!("{:?}", aligner.mapopt);
        println!("{:?}", aligner.idxopt);

        // tot len
        let seq = b"AAAAGAGAGAACGAAAAAAAGAGAGAGATGAATCCAGCGAAAACCATGTACGGCTTAAATTCAGCCACCTTATTATCGCCAGAATCCCTTTTCGCACGAACTGTCCTTTGATGCGAAACAAAATTTGGCGCATTCACTGTGCATACCGGCCGGTCAAGCAGCAGGTAACAATGCGGAGCGGCGCGCCTGCCGTCCATGATCGGGATTTCCGGCGCCGTTAGACTTCGAAAACAATGTTATCGATGCCCAAAGCCCCGGAGAGCAGCATTGAGTGCCATCGGTTGAAATCCGTACATCTCTCGTGACCAGACACGTACAGAGCATGGTATCACGCACAGAATTTGGCATTCGGCCGGGGAAATCTACCGTGATTCAAAGTCGGTGCGACGATAGATGACCCGGTGTTGCCGGCGCAGGGCGTAACGTCAGGGTGACTTTTTGCCGGTATGTAAACCACAAACCCAGTCGCCTGAACAGAACGTTTATGGTCCTTGTATTGATCATCGTATTATCTCAGCCAATTACCTATCCAACCGAAGTGTACTATACATTCGGCCGCAGTTTTAGCACAAAAGAGCCTCGAAACCCAAATTCCAGCAATTCTTATTCAGCTTGCTTACGAGGAATCTGGGAATCCAGATAATCCGGCTCTTCCATTTGCGCGCATTGTCCAATTCACGACTTTAGCTTCCGCTTCTGCTCCCTGGGTCAGCGGAGCCCATCCCATGCTGCTGTAGCGATCCATCACTGGCTGCTGAACCTGCTTATTGGTCACCAGATCTTCAGAACGAACCATGCTGAATAGCCCTGCCATTGCCGCATCCAGTCCGCTCAAGGATGTTCCCATAATACGCCCCCCCCCCCCCCAGTTACGCATCGACCTTGCTTCCATCAGGCTGACCAACGGGGCTGGAAATACGTTTCGGGTTGCTTCTTTCATCAGGCCAGACGTGACCCGTGCATCACCGCAAATCAAACGGTCTGCTGGGTAAATCCTGTTTGTAAGGAGCATGACGCCAGCATTTGTAATTTGCATGATCGGTAACCTGGGCATGATTCATAAACACCACTGCAAATTTTTGTGTCGTGCCTGGTCTACTAGTCGTAAAAATTGATCGCGGACAAATTTCGCCAGCATCTCTCTCAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGTGTCTGGGCGAATATTTCCGCGTTCAATTTTTTACGACATAGTAGCTCCCAGGCACGACAGCAAAAATTTGCAGTGGCTGTTTATGAATCTGCAGGTTTTACCGTCTGCAAAATTACAATGCTGGCGTCATGCTCGTACAAATCAGGATTTACCAGCGAGACGTTTTGATTTAGCGGTGATCGCCACGGGTCGACGTCTGCTGATGAAGAAAGAAAGCAACCCGAACGTATTTTTCTCCAGCCCGTGTCAGGCCTGATGGAAGCAAAGGTCGATGCGTGTAACGTGGGTATTATGGGAACATCCTTGAGCGGACTGGATGCGGCAATGGCAGGCTATTCACATGGTTCGTTCATTGAAGATCTGGTGACCAATAAGCAGGTTCAGCAGCCATGGATGGATCGCTACATGCAGCATGGGATGGCTCCGCTGAAGGAGCAAAGCCGGTGCTAAAGTCGTGAATGACAATGCGCCGCAAACTGCGAAGAGCCGGATTCCTGGATATCCCTGCATTCCTGCGTAAGCCAACTGATTAAAGAATTGACTGGAATTTGGGTTCGAGGCTCTTTGTGCTAACTGGCCCCCGAATGTATAGTACACTTCGGTTGGATAGTAATTGGCGCAGATATTTCATGATCAAACAAAGACACTTAAACGTATCGTTCAGGCGTCGGGTGTCGGTTTACATACCGGCAAGAAAGTCACCCGCGATTACGCCTGCGCCGGCCAACACCGGGGTCATCTAACGTTCGCACCGACTTGAATCCACCGGTTAGATTTCCCGGCCCGATCCAAAATCGTGCGTGATACCATGCGCTGTACGTGTCTGGTCAACGGCAGAGTACGGATTTCAACCCGTAGAGCATACCTCAATGCTGCTCTCGCGGGCTTGGGCAATCGATACATTGTTATCGAATTAACGCGCGGAAATCCCGATCATGACGGCAGCGCCGCCCGTTTGTATTACTCTGCTTGACGCCGGTATCGACGAGTTGAACTGCGCAAAAAATTGTCCGCATCAAAGAGACTGTTCTGTCGTTGATGGCGTAAGTGGGCTGATTTAAGCCGTACAATGGTTTTTCGCTGGATTCATCTCCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGATTAAAGAGAGAGATGAAATCCAGCGAAAACCATGTAGGCTTAACATTCAGCCACTATCGCCTCTTCGACAGAACAGTCTCTTTGATGCGAACAATTTTTTGGCGCAGTTCAACTCCGTCGATACCGCGTCCAAAGCAGCAGGTATACAAATACAGAGGCGGCGCTGCGTCATGATCGGGATTTTCCGGCGCGTTAACTTTCGATAACAAGTTATCGTGCCCAAGCCCCGGAGAGCAGCATTGAGGTGCTCTACGGTTGATAATCCTACATCATGCTCCGTTGACCAGAAACACGTACAGACATGGTATCACGCACAAGATTTGGCATCGGCGGGAAATCTACCTCGTGGATTCAACGTCGGTGCGACGTAGGTGTCCCCGGTGTTGGCCGCGCAGGGTCACGGCAGGTGACTTTCTTGCCGGTATGTAAACCGACACCCGTCCGCCTGAACGATACGTTTAAGTGTCCTTTGTTTGATCATCGTATTATCTCGCCAATTACCTATCCAACCGAAGTGTCTAATACATTCGGCGGGCATTTAGCACAAAAGAGCCGAAACCCAAATTCCCATCAATTCTTAAATCAGCTTGCTGACGCAGAATGCTGGATATCCAGAATCCGGGCCTTTCGCAGTTTGCGCCGCATTGTCATCAAGACTTAGCACGGCTTCGCTCCTGGGTCAGCGGAAGCCATCCCCATGCTGCTGTAGCGATCCTCACTGGCTGCTGAACCTGCTTATTGTCACCAAGATTCCTTCTTTGAACCGAACCATGCTGAATAGCCATCCATTGCCGCATCCAGTCCGCTCAGGAGTTCCCAAAATACCCACGTTACACGCATCGACCTTTTGCTCCACAGGCCTGACCACGGGCTGGGAAAATACGTTCGGTTGCTTCTTCTTCACAGGCCAGACGTGACCCCGTGGCGATCACCGCTAAATCAAACGTCTCGCTGGTAATCCTGATTTGCAGCGAGCATGACGCCAGCATTTGTAATTTGCAGATCGTAACTGGCATGATTCATAAACAGCCACTGCAAATTTGCTGTCGTGCCTGGTCTACAGTCGTAAAAATTGATCGCGGAAATGATTTCGCCCAGCATCCTCTCCTCAAACAACAACGGATGGAGGAGGAAAATAACAGGAAAAAAGAGAGAGATGCTGGGCGAATATTTGGATCAATTTTTTACGACCTAGTAGACCAGGCACGACAGCAAAATTGCAGTGCTGTTTATGAATCATGCCAGTTACGGACTGCAAATTACAAAGCTGGCGTCATGCTCGCTACAAATCAGGTTTACCCAGCGAGGATTTTGATTTATAGCGTGATCGCCACGGGTCACGTCTGGCCTGATGAAGAAGAAGCACACCGAACGTATTTTTCCCAGCGTGGTCAGCTGATGGAAGCAAAGGTCGAATGCAGTAGTAACGTGGGTATTATGGGAACACCTAGCGGACTGGATGCGGCAAGCAGTGGCTATTCAGCATGGTTTCGTTCATTGAAGATTGGTGACCAATAGCAGGTTCAGCAGCCAGTGAGGAGCGCACCAGCAGCATGGATGGCTCCGCTGACCCAGGAGCAGAAGCCCGGTTGCTAAATCGTGAATGACAATGCGCCGCAACTGCGGAAAGAGCCGGATTAGATCTGGATATCCCAGCATCCTGCGTAAGCAAGGCTGATTAAGAATTGACTTGGAATTGGGTTTCGAGGCTCTTTGTGCTAAACTGGGCCCGACCCGCCGAATGTATAGTACACTTCGTGGATAGGAATTTTGGGCGAGATAGATAGTGACCTATACGATGATTTCAAACAAAGGACACTTAAAACGTATCGTTCAGGCGTCGGGTGTCGGTTTACAATACCGGCAAGAAAGTCACCTGACGTTACGCCCTGCGCGGCCACACCGGGTCATCTATCGTCCCACCGACTGAATCCACCTAGATTGTCCCGGCCGATGCCAAATCTGTGCGTGATACCTGTTAAGCTCTTGTAACGTGTCTGGTCAACGAGCATGAGTACAGGATTTCAACCGTAAGAGCTCCTCAATGCTGCTCTCGCGGGCTTGGGCATCATAAATTTGTTATCGAATTAACGACGCCGGAATCCCGATTCATGGACGGCAGCGCCCGCTCCGTTTGTATACCTGCTGCTGGACGCCGGTATCGACGGAGTAGAAGCTGCGCCAAAAAATTTGTTCGCATCAAAGAGACTGTTCGTGGGTCCGAAGATGGCGATAATGGGCTGAATTTAAGCCGTACAATGGTTTTTTCGCGGATTCATCTCTCTCCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGAAATCCAGCGAAAAACCATTGTACGGCTTAATTCAGCCACTTTATCGCCATCTCAGACACGAACAGTCTCTTTTGAGAACAAATTTTTTGGCGCAGTTTCAACTCGTCCGATACCGGCGTCAAGCAAGCTGGTAATACATAGGTGCGGCGCTGCCGTCAAGATCGGGATTTCCGGCGTCGTTAACTGCGATAACAATGTTACGATGTCCCAAGCCGCGGAGAGCAGCATTGAGGTGCTCTAAACGGTTGAAAATCCGTAACATCATGCTCGTTTGACCAGACACGTAAACAGAGCATGGTATGCCACGACAGATTTGGCATCGGCGGGGGGGGGGGGGGGGGGGGGGGGGGGAAATCTACCGTGGATTCAAGTCGGTGCGACGAAGATGACCCCCGGTGTTGGCCGGCGCAGGGCGTAAACGTCAGGTGAACTTTCTGCCGGTATGTAAACCGACCACCCGTCGCTGAACGATACGTTTAGTGTCCTTGTTGATCATCGTAATTATCTCGCCAATTACCTATTCAACCGAAGTGTACTAGTACATTCGGCTGGCCAAGTTTAGCAAAAGAGCCTCCGAAACCCAAATCCAGTCAATTCTTAATCACTTGCTTACGCAGGTAGCTGGGATATCCAGATAATCCGCTCTTTTCGCAGTTTGCGCGCATTGTCATTCCACGACTTTAGCACCGGCTTCTGCTCCGGGTCAGCGGAGCCATCCCATGGCTGCTGGTAGCGATCCATCACTGCTGCTGAACCTGCTTATTGGTCACCAGATCTTCAATGAACGAACCATGCTGAAATAGCCACTGCCATTGCCGCATCCAGTCCGCTCAAGGATGTTCCCATATATCCCACGACACGCATCGACCTTTGCTTCCTCAGGCCTGACCACGGCTGGGAAAAATACGTTCGGGTTGCTTCCTTCTTCATCAGGCCAGACAGTGACCCGTGGGTTCCCCACCGCTAAATCAAGTCTCCTGGGTAATCTGATTTGTAAGCGAGCATGACGCCAGCATTTGTATTTTGCAGATCGGTAACCTGGCATGATTCAAAAGCCCTTGCAAATTTTTGCTGTCGTGCCTGGTCTACTAGTCGTAAAAAATTGAATCCGGAAATATTCGCCACATCTCTCTCACAACAACAACGAGGAGGAGGAAAAAAGGAGAGAGAGCTGGGCGAATAATTTCCGCGATCAATTTTTACGACTAGTAGTCCAGGCACGACGCAAAAATTTTGCAGTGCGGTTTATGAATCATGCCAGGTACGATCTGCAAATAAATGCGGCGTCATGCTCGCTACAAATGCGAGGATTTTACCCAGCGAACGTTTGATTTAGCGGGTTCGCCCTCGGGTCACGTGCTGCCTGATGAAGAAGAAGCAACCCGAACGGTATTTTCCCAGCCCTGGTCAAAGGCTGATGGAAGCAAAGTCGATGCGTGTAAACGTGGTAAAATTATGGGAACATCCTTGAGCGGACTGAGATGCGGCAATGGCAGTGGGCTATTCAGCATGGTTTCGTTCATTGAAGATCTGGTGACCCATAAGCAGGTCAGCAGCCAGTGATGGATCGTCTACCAGCAAGCATGGGAATGGCTCCGCTGACCCAAGGAGCAGAGCCGTTGCTAAAGCGTGAATGACATTGCGCCGCAAAACTGCGAAAGAGCCGGATTATCTGGATATCCCAGGCATTCCTGCGTAAGCAAGCTGATTAAATTGACGGAATTTGGGTTTCGAGGTCTTTGTGCTAAACTGGCCCCGCGAATGATAGTACACTTCGGTGGATAGGTAAATTTGGCGAGATAATACGATGATCAAACAAAGGACATTAAACGTATCGTTCAGGCGCGGGTGTCGGTTTAACATACCGGCAAGAAGTCACCTGACGTTAGCCTGCGCCGGCCAACAACCGGGGTCATCATCGTCGCACCACTGAATCACACCGTGTAGATTTCCCGGCCCGATGCCAATCTGTGCAGTGATACCTGCTCTGTACGTGTCTGGATCAACGGCATGATGTAACGGATATGCAACGTAGATAGAGCAGCCTCAATGCTGCTCTTCGCGGGCTTGGGCATCATAACATTGTTATCGTGGTTAACGCGCCCCGGAATCCCGATCATGGACATAGGCAGCGCCGCTCCGTTTGTATACCTGCTTGCTTGACGCCGGTATCGCGAGTTTGAACTGCGCCAAAAAATTTGTTCGCATCAAGAGACTTGTTCGTGTCGAAGATGCGATAAAGTGGGCTGAATTTAAGCCGTACAATGGTTTTCTCGCTGGATTTCATCTCTCTCAAACAACAACGGAGGAGAGGAAAAAAGAGAGGATGAATCCCAGCGAAAAACCATTTGTACGGCTTAAATCAGCCCACTTATCGCCATCTTCGACACGAACAGTCTCTTTGATGCGAACAATTTTTTGCGCAAGTTCAACTCGTCATACCGCGTCAAGCATGCAGGTATACAAACGGAGCGCGCTGCCGTCCATGATCGGGATTTCCGGCGCGTTAACTTCGATAACAATGTTATCGATGCCCAAGCCCGCGGAGAGCAGCATTGAAGGTGCTCTACGGTGAAATCCGTACATCAGCTCGTTGACCAACACGTACAGAAGCATGGTATCACGCACAGATTTGGCATCGGCCGGAATCTACCGGTGATTCAAGTCGGTGCGAAGATGAATGACCCCCGGTGTTGGCCGGCGGCAGGGCGTAACGTCAGGGTGACTTCTTGGCGGTATGTAAACCCACAGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGTCGCCGAACGATACGTTTAAGGTCTTTTGTTTGATCATCGATTATCTCCCAAATTACCTATCCAACCCGAAGTGTACCTATACATTCGGCGGCCGTTTAGCACAAAGAGCCCGAAACCCAATTCCAGTCCATTCTTATCAGCTTGCTTACCCAGGAATCTGGATATCCAGAAATCCGGCTCTTTCCCGAGTTTGGCGCATTGTCATTCACGACTTTAGCAACCGGCTTCTGCTCCTGGGTCCAGCGAGGCCCTTCCATGCTGCTGTAGCGATCCATCACCTGGCTGCTGAACCTGCTTATGGTCCCACCAGATCTTCAATGAACGAACCATGCTGGAATACCACTGCCTTGCCGCATCAGTCCGCTCAAAGGATGTTCCCATAATTCCCACAGTTACAGCATCGACTTGCTCCATCAAGCCTGACCACGGCTGGGAAATACGTTCAGGTTGCTTTTCTTCATGCAAGGCCCAGAACCGTGACCGTGGCGATCACCGCTAAAATCAAACGTCTCGCTGGGTAAATCCTGATTTGTAGCGAGCATGACGCCAGCATTTGTAATTTGCAGATCGGTAACTGGATGTTCATAAACACACTGCAAATTTTTGCTGTCGTGCCAGGTCTACTAGTCGTAAAAATTGATCGCGGAATATTCGCCCAGCATCTCTCTCAAACAACAACAGAAGGAGGAGGAAAAAAGAGAGAGATGCTGGGCGAATAATTTTGCCCGATCAATTTTTTACGACTAGTAAGAGCCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCAGCACGCAGCAAAAATTTGCAGTGGCTAGTTTATTGAATCATGCCAGGTTTCCGAACTGCAAATTACAAAATGCTGGGTCAGCGCTACAAAATCAGGTTTTAACCCAGCGAGACGTTTTGATTTAGCGGTGATCGCCACGGTCTCGTCTGGCCCTGATGAAAAGAAGCAACCCGAACGTATTTTCCAGCCCGTGGTCAGGCCTGATGGAAGCAACAGGTCGATGGCGGTACACGTGGGTATTATGGGAACATCCTTGAGCGGACTGTGCGGCAATGCAGTGCTATTCAGCATGGTTCGTCCATTGAAGATCTGGGACCAATAGCTGGTTCAGCATGCCAGTGATGGATCGCTACCAGCAGATGGATGGCTCCCGTACCCAGGAGCAGAAGCCGGTTGCAAAGTCGTGGAATGACAATGCGCCGCAAACTGCGAAGAGCCGGAAATCTGGATAATCCCCCAGCATTCCTGCATAAGCAAGCTGATTAAGAATTGACTGGAATTTGGGTTCAGGCTCTTGTGCTAAACTGTAGCCCCGCCGATGATAGTCACTTCCGTTGGCTAGGTAAATTTGGCAGATAATCGATGATCAAACAAAGGACACTAAACGTATCGTTCAGTGGGCGACGTGTCGGTTTACATACCGCAGAAAAGTCCACCCTGACGTTACGCCCTGGCGCCGCAACACCGGGTCATCTAACGTCGCACCGACTTGAAATCCACCGGAGATTCCGGGCCAGATGCCAAATCTGTGCGTGGATACACAAGCTCTGTACGTGGTCTGGTCAACGACATGATGTACGGATTTCAACCGTAGAGCACCTCAATGCTGCTCTCGCGGGCTTGGGCATCATAAACATTGTTATCGAAGTTAACGCGGCCGGAAATCCGATCATGGCGGCAGCGCCGCTCCGTTTGATACCTGTGCTGCTTGACGCCGTATCGACGAGGTTGAACTGCGCCAAAAAATTTGTTCCGCATCAAAGAGAACTGTTCAGTGTCGAAGATGGCGAAAGTGGGCTGAATTTAAGCCCGACAATGGTTTTCGCTGGATTTCATCTCTCTCAAACACAACGGAGGAGGAGGAAAAAAGAGAGAGATGAAATCCAGCGAAAAACATTGTACGGCTTTAAATTCAGCCCACCTTATGCATCTTCGACACGAACAGTCCTCTTGATGCGCAAATTTTTGGCGCAGTTCACTCGTCATACCGGCGTCAAGCAGTGTATCAAATATCGGAGCAGGCCTGCCGCCATGATCGGGATTTCCGGCGCGTTAACTTGATAACAATGAATCCGAATGCCCAAGCCCGCGAGAGCGCATTGTGGTGCTCTACGTGTTTGAAATCCCGTAACATCATGCTCGTTGACCGACAGACAGAGCCATGGTATCACGCACAGATTTGGCATCGGCCGGGAAATCTACAGGTGGATCAAGTCGGGCGACCGATAGATGACCAGTGTTTGGCCGCGCCAGGCGTAAACGTCAGGGTGACTTTCGCCGGTATGTAATCCGACACCCGTCCCCGCCTGAACGATACGTTTAAGTGTCCCTTTTGTTTTGATCATCGTATTACTCGCCAAATTACCTATCCAACCGAAATGGTACTATACATTGGCGGCCAGTTTTAGCACAAAGAGCCTCGAAACAAATTCCCGTTATTCCTTAATCAGCTTGCTACGCAGGATGCTGGGATATCCAGAAATCGGCTCTTTCCAGTTTGCGGCGCATTGTCCATTCACGGACTTTAGCAACCGCTTCTGCTCCTGGGTCAGCGGTGCCATCCCTGCCTGCTGTAGCGATCCATCACTGGCTGCTGAACCTGCTTTATTTGGTCACCAGATCTTAATGACGAACCATGCTGAATAGCCACCTGCCTTGCCGCTCCAGTCCGCTCAAGGTGTTCCCATAATACACTTACACGCATCGACCTTTGCTTCACAGGCCTGACCACGGGCTGGGACAAATACGTTCGGGTTGCTTCGCTCTTCATCAGGCCAGACGTGACCGGCGATCAACCGCTAAATCAAACGTCTCGCTGGGTAAATCGCTGATTGTAGCGAGCAGACGCCAGCATTTGTAATTTGCAGATCGGTAAACCTGGCATGATCATAAACACCACTGCAAATTTTTGCTGTCGTGCCTGTCTACTAAGTCGAAATTGATCGCGAATATTCGCCCAGCATTCTCTCCAACAACAACGGAGGAGAGGCAAAAAGGAGAGAGATGCTGGGCAATATTTCCCGCGATCAATTTTTACGACTAGTAGACAGGGGCACCGACAGCAAAAAAATTTGCAGTGGCTGTTTATGAATCATGCCAGTTACGATCTGCAAATTACAAATGCTGGCGTCATGCTCGCCTACAAATCAGGATTTACCCAGCGAACGTTTGATTTAGCGGTGATGCGCCCGGTGTCACCGTCTGGCCTGAGAGAAGAAGCAACCCGAACGTATTTTTCCCAGCCCGTGGCAGGCCTGTGGAAGCAAAGGTCGATGCGTGAACCAGTGGGTATTATGGTGACATCCTTTGAGCGGACTGATGCCGCAATGGCAGTGGCTATTCACATGGTACGTTCATTGAAGAACTGTGACCAATAAAAAGCCAGGTTCAGCAGCCAGTGATGGATCGCTACCAGCAGAGGGATGGCTCCGCTGACCCAGGAGCAGAAGCCGGTTGCTAAAGTCTGATGACAAGGCCGCAAAACTGCGAAAGGAGCCGGATATCTGGATATCCCACATCCTGGCGTAGCAAGGCTGATTAAGAATTGCTGAATTTGGTTCGAGGCTCTTTGTGCTAAACTGGCCGCCGAAATGTATAGTACATTGTTGATAGGTAAATTTTTGGCGAGATAATACGATGATCAAACAAAGGACACTTAAACAGTACGTTCAGGCGACGGTGTCGGTTTACATACCGCAAGAAAGTCACCCTGACGTTACGCCCTGCGCCAGCCAACACGGGGTCATCTATCTCCACGACCTTGATCCGACCGGTAGATTTCCGGCGATGCCAAATCTGTGCGTGATACATGCTCTGTACGTGTCTGGTAACGAGCATGATGTAGGTTTTCAACCGTAGATCACCTCAATGCTGCTCTCGCCGGGCTTGGGCATCGATAACATTGGTTATCAAAGTTAACCGCCGGGAATCCGTCTGGACGGCAGCGCGCTCGTTTGTATACCTGCTGCTTGACCCGTACGAGAGTTGAACTGCCCAAAAAATTTGTTCGCATCAAAGAGACTGTTCGTGTCGAAGATGGCGATAAGTGGGCTGAATTTAAAGCCGTACACATGGGTTTTTCGCTGGATTTCATCTCTCTCAACAAAACGGAGGAGAGGAAAAAAGAGAAGATGAAAATCCAGCGAAAAAACCATTGTCGGCTTAAATCAGCCCACTTATCGCCATCTTCGACACGAACAGTCTCTTTGATGCGAACAAATTTTTTGGCGCAGTTCAACTCGTCGAATACCGGCGTCCCAAGCACAGGTATACAAACGGAGCGGCGCTGCCGTCCTGGATCGGGATTTCCGCGCGTTAACTTCGTAACAATGTTATTCGAAAGCCCAAAGCCCGCGAGAGCAGCATTGAGGTGCTCTACGGTTTGAAACCCGTACATCCATGCCTCGTTGACCAGACACGTTACAGGAGCATGTATCACGCACAGATTTGGCTCGGCCGGGAAATCTACCGGTGATTCAAGTCGGTGTCGACGATAGATGACCCCGGTGTTGGCCAGGCGCCAGGCGTAACGTCAGTGACTTTCTTGCCGGTATGTAAACCCGACCCCACCCTCGCCTGAACGATAGCGTTTAAGTGTCCTTTGTTGATCATCGTATTATCTCGCCAATTACCTATCCAAACCGAAGTGTACTATCACATTCGGGGCCAAGTTTAGCACAAGAGCTCGAAACCCAAAGTTCCAGTCATTCTAATCAGCTTGCTTACGCAGGAATGGCTGGGATATCCAGATAATCCAGCTCTTTCGCAGTTTGCGCGCATTGTCATTCACGACTTTAGCAACCGGCTTCTGCTCCTGGGTCAGCGAGCCATCCCATGCTGCTGGAGGATCCATCACTGCCTGCTGAACCTGCTTATTGGTCACCAGTTCTTCAATGAACGAACCATGCTGCAATAGCACTGCCATTGCCGCATCCAGTCGCTCAAGATGTTCCATAATACCCACGTTACGACGCATGACCTTGCTTCCATCAGGCCTGACCACGGCTGGGAAAATACGTTCGGGTTGCTTCTTCTTCATCAGGCCAGACAGTGGACCCGTGGCGAGCACCGCTAACAAACTCTCGCTGGGAAATCCTGATTTGTAGCGAGCATGACGCAGCATTTTAATTTGCAGAGTCGGTAACCTGGCATGATGCAAAACAGCCACTGCAAATTTTTGCTGCGTGGCCTGGTCACTAGTCGTAAAAATTGATCGCCGGAATAAATATTCCGCCCAGCATCCTCTCAACAACAACGGAGGAGGAGGAAAAAAGAAAGAGATGCTGGGCGAATATTCCGCGATCATTTTACGACTAGTTGTCCTGGCACGACAGCAAAAAATTTGCAGTGCTGTTTATGTTGCTTGCCAGGTTACCAACCTCTGCAATTTACAAAGCTGGCGTCAGCTCGCTACAAAATCAGGATTTACCCAGCGAGACGTTTGATTACGCGGGATCGCCACGGTCACGTCTCCTGATGAAGAAGAAGCAACCCGAACGTATTTTCCACAGCCGTGTCAGGCCTGATGGAAGCAAAGGCGATGCGTGTACGTGGGTATTATGGGAACATCCTTGAGCGGACTGGATGCGGCAATGGCAGTGCTATTCACATGCGTTCATGAATGCATCCTGGTGACCAAATAAGCAGTTCATGCAGCCAGTGATGATCGCTACCAGCAGCTTGGGATGGCTCCGCTGACCAGGACAGAGCCGTTGCTAAAGTCGTGATGACAATGCGCCGCAAACCGCGAAAGCCGGATTATCTGGATATCCCAGCATTCCCTGCCGTAATGCAAGCTGATTAAGAATTGCTGGATTTGGGTTCGAGGCTCTTGTGCTAAACTGGCCGCCCCGAGTATAGTATATTTCGGTTGGATAGGTAATTTGCGAGATAATAACTTGATCAAACCAAAGGACACTTAAACGTATCGTTCAGGCGACGGGTGTGCGGTTACATACCGCAAGAACAGTCACCCTGACGTTACGCCCTGCGCCGGCCACACCGGGGTCACTATCGTCGCACACTTGATAATCCACCGTAGAATTTCCCCGGCGAGCCAAATTGTGCGTATACCATGCTCTGGGTAAACGTGTCTGTCAACGGCATGATGTACGATTTCAACCGTAGAGCACCTCAATGCTGCTCTGCGGGCTTGGCATCGATAACATTGTTACAAAGTTAACGCGCCGGAAATCCCGTCACTGGACCGCAGCGCCGTCTCCGTTTTGATACCTGCTGCTTGACGCAGGTATCGACGAGATGAACTGGCGCCAAAAAATTGTTCGCATCAAAGAGACCTGTTCGTGTCGAAGATGGCGTAAGTGGGCTGAATTTATGCCGTACAATGGTTTTTCGCTGGATTTCACCTCTCTCAATCAAACAACAAGGAGGAGGAGGAAAAAAGAGAGAGATGAAACCCCAGCGAAAAACCATGTACGCTAAATTTCAGCCCACTTTCGCCATCTTCGACACGAACAGTCTCTTGATGCGAACAAATTTTTTGCGCAGTCAACTCGTCGTACCGGCGTCAAGCAGCAGGTATACCAAACGGTGCGCGCTTGCCGTCCATGATCGGGATTTCCGGGCCGCGTAACTTTCGTAAACATGTTATCGAGCCCAAGCCCGCGAGAGCGCATTGAGGTGCTCCTTACGGTTAAATCCGTACATCATGCTCGTGACCAGACACGTACAGAGCATGTATCACGCACAGAATTTGGCAATCCGGCCCCCGGGAAATCTACCGCTGGATTCAAGTCGGTGCGACGATAGTGACCCCGTGTTGGCCGGGCCGCCAGGGCGTAACGTCAAGGGGGTCCTTTTTCTTGCCCCCCGGTATGTAAAACCCGACACCCGTCCCCCCGCCCCTGAACGAAACCCCGTTTAAGTGGTTGTTTGGGAACAATCGTAATTTCTCCCATTTACTAAATCCCCCCCGAAACGGGGGGGGGTTTTAAACAA";
        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}, strand:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score,
                hit.strand
            )
        });

        // // 3353~3662
        // let seq = b"CACTGCCCCCGACCACGGCTTCACACAATATTGAGGGTCAGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGATGGAGCGCGATATTTCCCACGACAGCTCAGCGAAACAAAAAAAAAAAATTCCGCCACATGCGCCCGGCATGGCGATAAAACGCCACCTCCAGATGGACGCGGGAAACCACAGTGGGGCGAGTTGGGTGCCAAACTGCGAAAGTCAAGCAGAAAAGGGTAAAGAAAACGAGCTAAACTAACTGCACAGGTGAAGCAGCGGGTGGGCAGGCATAATACATGTGTGATTGTTTG";

        // let mut mapping = aligner
        //     .map(
        //         seq,
        //         false,
        //         false,
        //         None,
        //         Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
        //         Some(b"t"),
        //     )
        //     .unwrap();
        // mapping.sort_by_key(|hit| hit.query_start);
        // mapping.iter().for_each(|hit| {
        //     println!(
        //         "qstart:{}, qend:{}, primary:{}, supp:{}, identity:{}, score:{:?}",
        //         hit.query_start,
        //         hit.query_end,
        //         hit.is_primary,
        //         hit.is_supplementary,
        //         MappingExt(hit).identity(),
        //         hit.alignment.as_ref().unwrap().alignment_score
        //     )
        // });
        //
    }

    // #[test]
    // fn test_align_single_query_to_target2() {
    //     let ref_file =
    //         "/data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.fasta";
    //     let fa_iter = FastaFileReader::new(ref_file.to_string());
    //     let targets = read_fastx(fa_iter);
    //     let aligners = build_aligner(
    //         "map-ont",
    //         &IndexParams::default(),
    //         &MapParams::default(),
    //         &AlignParams::default(),
    //         &OupParams::default(),
    //         &targets,
    //     );

    //     let mut bam_reader = rust_htslib::bam::Reader::from_path(
    //         "/data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.bam",
    //     )
    //     .unwrap();
    //     bam_reader.set_threads(10).unwrap();

    //     let aligner = &aligners[0];

    //     for record in bam_reader.records() {
    //         let record = record.unwrap();
    //         let record_ext = BamRecordExt::new(&record);
    //         let seq = record_ext.get_seq();
    //         for hit in aligner
    //             .map(seq.as_bytes(), false, false, None, None, None)
    //             .unwrap()
    //         {
    //             println!("{:?}", hit);
    //         }
    //         break;
    //     }
    // }

    #[test]
    fn test_align_single_query_to_target4() {
        let ref_file = "test_data/MG1655.fa";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &AlignParams::default(),
            &OupParams::default(),
            &targets,
            10,
        );

        let aligner = &mut aligners[0];

        // aligner.mapopt.min_cnt = 2;
        aligner.mapopt.best_n = 100000;
        // aligner.mapopt.min_dp_max = 10; // min dp score
        // aligner.mapopt.pri_ratio = 0.5; // min dp score
        aligner.idxopt;
        println!("{:?}", aligner.mapopt);
        println!("{:?}", aligner.idxopt);

        // tot len
        let seq = b"AAAAGAGAGAGATAGGGGGCGACGAGCGCCAAACGGGCGGCAATTATAGGGATTTCATCGCCTGATACCAGTCGAATAGCGTTGCCGCGCGCTCAGGATGTTAATTGGTTGACAGAAGAATTCCCGGTGGCAAAATTACGTTGAAGATCAGTTTTATACAAGGTAAAAAATGTTATACGCAGTTGCGCAATTATCGCCTTTACGTCACTTTATGAGCATTCGCATATAAAATGTAAAACTTTTTGTAACTAGCATAAAACACAGAAACGAATACCTGGCGCTGGTCTTGCGATAAAGCAGGTAATGAGCAAACAACACAGCATGTATTAATTGCCCTGCCCACCCGCTGCTTCCACCTGGTCCAGTTTAAGGTTAGTCCTCGTTTACTTTACCCTTTTCTCGCTGAGCTTTCGCAAGTTTGGCACCCTCTCGCCCCACGCCTGTGGTTCCGACGGTCCACTGATGGTGGCGTTTATCGCCATGCCGGGGCGCAATGTGCCGGGAAATGCGCTGAGCTGTTCGCTGGAAATATCGCCGCATCCATCCTGCTTTTTTCCACCAGCCGCTGAACATGACCTGGACGACATCAATATTGTTGAAGCCGTGGTCGGGCAGTGCATCTCTCGCTCAAACAACAACGGAGGAGGAGAAAAAAGAGACAGATGCACTGCCCCGACCACGGCTTCAACACATATTGATGGTCGTCCAGGTCATTTCCACGAGTGGTGGAACAAAGCAGGAGGATGCGGCGTATTTCCAGGAACAGCTCACGCAATTCCCGGCACATGGCGCCCGGCAATGGCGAATAAAACGCCAACCATCATGAATGGACGTCGGAACCACATGGGCGATTGGTGCCAAACTGCGAAAGCTCAAGCGAGAAAAGGGTAAAGAAAAACGAGACTTATACCTAAATGACCAGGGAATGCAGCGGGTGGGGCAGGGCAATTAATTACATGGCTGTGATTGTTGCTCATACCGCTTTATCCCAAGACAGGTCGCCAGTATTCGTTTCTGTGTTTAATGCTATACAAAAGTTTATTTACATTTATATGCGAAATTGCTCATAAGTGAGGGTAAAGGCGATAATTTGCGCAACTGCGTAATAACATTTTTTTACCTTGTATAAAACTGATCAACGTAAAGTTGCCACCGGGATCTTCGTCAAAATTAACTCTGACGCGCGGCAACGCTTATTCGACTGGTATCAGGCGGATGAAATCCTATAAATTGCCGCGTTTGGCGCTTCCGTCGCCCCCTATCTCTCTCAAACAACAACGGAGAGGAGGAAAAAAGAGAGAAGATAGGGGGCACAGAAGCGCCCAACGCGCAATTAATAGGGATTCATCCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCATGATACCAGTCGAATAGCGTTGCGCGCGCTCAGAGTTACATTGTTGACGAAGAATTCCCGGGTGGCAAATTACGTTGATCAGTTTTATAACAAGGTAAAAAGTTATACGCAGTTGCGCAAATTATCCGCCTTACGTCACTTTAATGAGCAATTCCGCATAATAAAAATGTAAAACTTTGTAACTAGCATAAACACAGAAACGAATACTGGGGCGACCTGGTCTTGCGGAATAAGCGGTAATGAGCAAAAATCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGCTTCACCGGTCAGCGTTTAGGTTTAAGTCTCGTTTAATCTTTACCCTTTCTCGTTGAGCGAGCTTTCGAGTTTGCACCCAACTCGTCCCACTGTGGTTCCCGACGTCCATCATGTGGGCGTTTTAAAATCGCGCATGCCGGCGCATGTGGCCGGGAATTGCCTGAGCTGTTACGCTGGGAAATATCGCCGCATCCATCCTGGCTTTTTTTCCACAGCTACGCTGAACAGACCTGGACGACCACCATCAATATGTTGAAGCCGTGGTCGGGCGTGCATCTCTCTCAAACAACAACGGAGGAGGGAGGAAAAAAGAGAGAGATGCACTGCCCCGACCACGGTTTTTTGTCAACATATTGAAGGTCGTCCAGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAATCCCGGCCACAGCGCCCGGCAATGGCGATTAAAACGCCACCATCATGATGGAACGTCGGGAACCACAGTGGGGCGTGTTGGGTGCAAACTCGAAAAGCTCAAGCAGAAAAGGGTAAAGATAACACGAGACAAACCAAACTGACCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGTTTATCCGCAAAACAGGTCGCCAGTATTCGTTTCCTGTGTTTATGCTAGTACAAAAGTTTTACATTTTAATATGCGAATTGCTCATAAAGTGACGTAAGGCGGAAATTTGCGCAACTGCGTAATAACATTTTTTTAACTTGTATAAAACTGTTCAACGTAATTTGCCACGGAATTCTTCGTCAACAATTAACTACCTGAGCGCGCCAACGCTATTCGACTGGTATCAGGCGAATGAATCCCTATAATTGCCGCGTTTGGCGCTTCGTCCGCCCCCTATCTCTCTCAAACAACAACAGGAGGGGAGGAAAAAAGAGAGAGATAGGGTGGCGACGAAGCGCCAAACGCGGCAATTATAGGGATTTCATCCGCCTGATACCAGTCGAATAGCGTTGCCGCGCGCTCAGAGTTAATTGTTTGACGAAGAATTCCCGGTGGCAAATTACGTTGATCAGTTTTTATACAAGGTAAAAAAATGTTATACGCAGTTGCGCAAATTATCCGCTTTACAGTCCACTTTTAGAGCAATTCGCATAATAAAAGGTAAAACTTTTGTACTAGCATAAGAACACACACAGAAACAATACTGGCTGACCTGGTCTTGCGGATAAAGCGGTAATGAGCAAACATCACGCATGTATTAATTTGCCCCTCCCCACCGCTGCTTCAACCGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTAGTACCCTTTTCTCGCTTGAGCTTTGCAGTTTGGCACCCAACTCCCATGTGTTTCCGCGACGTCCATCATGATGGTGGCGTTTTATCGCATGCCGGGCGCATGTGGCGGGAATTGCGCTGAGCTGTTCGCTGGGACATATCGCCGGCATCCTCCCTGCTTTTTTCCCCAGCTCGCGCTGAACATGACCTGGACGACCAAGTCAATATTGTTGAGCCGTGGTCGGGTGCAAGTGCATCTCCTCAACAACAACGAGGAGGGGAAAAAAAAAGATGAGAGATGCACTGCCCCCGACCACGGCTTCACACAATATTGAGGGTCAGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGATGGAGCGCGATATTTCCCACGACAGCTCAGCGAAACAAAAAAAAAAAATTCCGCCACATGCGCCCGGCATGGCGATAAAACGCCACCTCCAGATGGACGCGGGAAACCACAGTGGGGCGAGTTGGGTGCCAAACTGCGAAAGTCAAGCAGAAAAGGGTAAAGAAAACGAGCTAAACTAACTGCACAGGTGAAGCAGCGGGTGGGCAGGCATAATACATGTGTGATTGTTTGGTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCATTACCGCTTTATCCCAAGACCAAGGTCGCCAGTATTCGTTTCTGGGTTTATGCTAGTACAAAAGTTTACATTTTATATGCGAATTGCTCCATAAAGTGGCCCCCCCCCCCCGTAAAGGCGGATAATTTGCGCAACTGCGTATAACATTTTTTACCTTGTATAAAACTGATCAAGTAATTTGCCACCGGAATTCTTTCGTCAACAATAACTCTGGCGCGCGGCAACGCTATTCGACTGGTATCAGGCGGAGAAATCCTAATAATTGCCGCGTTTGGCGCTTCGTCGCCCCTATCTCTCTCCAAACAAACAACGAGGAGGAGGAAAAAGAGAGAGATAGGGGCGACGAAGCGCCAAACGCGGCAAATTATAGGGATTCATCCGCCTGATAACCAGGTCGAATAGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCAGCGTTGCCCGCGCTCAATTTAATGTTGACGAAGAATTCCCGGTGGCAAATTAACACGTTGAATCAGTTTTATACAAGGTAAAAAATTTATAGCAGTTGCGCAATTAATCCGCCTTACGTCACTTATGAGCAATTCGCATATAAAATGCTAAAACTTTTTGTACTAGCATAAAACACAGAAACGAATACTGGCGACCTGGTCTCGCGGATAAAGCGGTAGATGAGCAAAACAATCACAGCACTGTATTAATGCCTGCCCCACCGCTGCTTCACCGAGTCAGTTTAGGTTTATAGTCTCGGTTATCTTTACCCTTTTCGCTTGAGCTTTCGCAGTTTGGCACCACACTGCCCATGTGGTTCCGCGTCCATCCATGATGGTGGCGTTTTTATCGCCATGCCGGGAGCGCATGTGGCCGGAATTGCGCTGAGCTGTTCGCTGGAAATATCGCCGCATCCGATCCTGCTTTTTTCCACCAGCTCGCTGAACATGACCGGACGACCATCAATATTGTGAAGCCGTGGTCGGGGCAGTGCATCTCTCTCAACAAACAAACAACAGGAGGAGAGAAAAAAGAGAAAGCATGCCCCGACCACGCTTCAACAATATTGAGTGGTCGTCCGGTCATGGTTCAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGCGATATTTCCCAGCGGAACACTCAGCGCAATTCCGGCCAAACAGCGCCCGGCATGGCTGATAAAACGCCACCATCATGATGGACGTGGGAACCACAAGTGGGGCGATTGGTGCCAACTGCGAAAGCTCAAGCGAGAAAGGGTAAAAGATAAACGATGACCTAAACCTAAACTGACCAAGTGAACAGCGGGTGGGGCAGGCAATTAATACATGTGTGATTGTTTGCCATTACCGCTTTATCCCAAGACCAGGTCGCCAGTATTCGATTTCTGTGTTTATGCTAGTACTAAAGTTTTACATTTATATGCAGAATTGCTCATAAAGTGAGTAAGGTCGGATAATTTGCGCACACTGCGTATAACATTTTTTACCTTGTATAAAACGATCCAACGTAATTTGCCACCGGGAATCTCGTCAACAATTAACTCTGAGCGCGCGGCAACGCTATTCGACTGTATCAGGCGGATGAAATCCCTATAATTGGCCGCGTTTGGCGCTTCGTCGCCCCCATCTCTCTCAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATAGGGGGCGACGAAGCGCCAAACGCGGCATATAGGGATTTCATCCGCCTGATACCAGTCCGAAAGCGTTGCCGCGCGCTCAGAGTAAGTTGACGAAGAATTCCCGTGTGGGCACAATAACGTGACAGTTTTATACAAGGTAAAAAATGTTATACAGCAGTTGCGCAAATTTCCCCTTTACGTCACTCACTTTATGAGCAATTCGCATATAAAATGTAAACCCTTTTGTACTAAGCATAAACAGCAGAAACGAATACCTGCGATACCTGGTCTTGCGGATAAACGGTAATGAGCACAACTCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGTCTTCAACCTGGTCATGTTTAGGTTTATCTCGTTTATCTTACCTTTCTCGCTTGAGCTTTCGCATATTGGCACCCAATCGCCCCACGTGGTTCCCGACGTCCATCATGAGGTGGGCGTTTTATCCCATGCCGGGCGGCAATGTGGCCGGAATTGCGCTGAGCTGTTCCGCTGGGAAAATCGCCGCATCCATCTAGCTTTTTTTCCACAGGCTCGCTGAACATGATACCTGGACGACCATCAATATTGTTGAAGCGGGTCGGCAGTGCATCTTCTCAAACAACAACCGTAGGAGGAGGAAAAAAGAGAGAGATGCACTGCCCCGACCAACGGCTTCACAATATTGATGGTCGTCCAGGTCATGTTCAGCGAGCTGGTGGAAAAAGCAGGATGGATGTCGGCGATATTTCCCAGCGAACAGCTCAGCGCAATTCCCGGCCACATGCGCCCGGCATGCGATATAAAGCCACATCATGATGGACGTCGGAACCACATGGAGCGAGTTGGGTGCCAAACTGTGGAAAGCTCAGCGACAAAAGGGTAAAGATAAACGAGACTAAACCTAAACTGACCAGTGAAGCAGCGGGTGGGCAGGGCAATTAATACATGCTGTGATTGGTTTGCTCATACCGCTTACGCAAGACCAGGCCGCAGTATTCTTTCTGTGGTTTATGCTAGTACAAAAGTTTTACATTTTATATGCGAATTGCTCAAAAGTGACGTAAAGCGGATAATTTGCGCAACTCGTATAACATTTTTACCTTGGTATAAAACTGATCAACGTAATTGCCACCGGGAATCTTCGTCACAATTAAGACCTGAGCGCGCGGCAACGTATTCGACTGGTATCCAGCGGATGCAAATCCCTATAATTGCCGCGTTGGCGCTTCGTCGCCCCATTCTCTCAAACAACAACGGAAGGAGAGGAACAAAGAGAGAGCATAGGGGGGCGGGGGCGAACGAAGGCAAACAGCGGCAATTATAGAGGATTTCACTCCGCCTATACCAGTCGAATAGCGTTGCCGCGCGCCAGAGTTAATTGTGACGAAGAATTCCCGTGGCAAATTAACGATGATCAGTTTTATACAGGGTAAAAAATGTTATTACAGCAGTTGCGCAAATTATCCGCCTTACGTCACTTTATGCAGCAATTCGCATATAAAATGTAAAACATTTTGTACTAGCATAAACACAGAAACGAATACTGGCGACCTGGTCTTGCTCGGATAAAGCGGTTAATGAGCAAACAATCACCAGCAGTATTAATTGCCCTGCCCCCACCGCTGCCTTCACCTGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTCTCGCTTTGAGCTTCGCAGTTTGGCACCCAACTCGCCCCCACTGTGGTTCCCGACTCCCATCATGATGGTGGCGTTTTATCGCCATGCGGGCGCATGTTGGCCGGGAATGCGCTGACTGTTCGCTGGGAATATCGCCGCATCATCCTGCTTTTTTCACCAGCTCGCTGAACATGACCTGGACGACATAATAGTTGAAGCGTGTCGGCAGTGCATCTCTCTCAAAACAACAACGAGGAGGAGGAAAAAAGAGAAGATGCCTGCCCCCGCACAGGTTAACAATATGAGGTCGTCCAGGTCGATGTTCAGCGGACTGGTGGGAAAAAAAAAAAAAGCAGGATGGAATGCGCATATTTCCCAACAGCTCAGCGCAATTCCCGGCCCACAATGCGCCCGGCATGGCGATAAAACGCCACCACAGAGGACGTCGGAACCACAGTGGGGCGAGTTGGGTGCCAAATCTCGAAAGCCAGCGAGAAAAGGGTAAAGATAAAACGAGACTAAACCTAAACTGACCAGGTGAAGCAGCCGGGTGGGGCAGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCCGGCTTTATCCCAGACCAGGTCGCCAGTATTCGTTTCTGTGTTTATGCTAGTACAAAAGTTTAACATTTTATAGCGATTGTTCATAAAGTGACGTAAAGGCGGATAATTTGCCAACTGCGTATAACGATTTTTACCTTGTATAAAACTGATCAACGTATTTGCACGGGAATTCTTCGTCCCCAATTGAACTCTGAGCGCGCGCACGCTATTCGACCTGGTATCAGGCGGATGAAATCCGTATAATTGCCCGTTTGCGCTTCGTCGCCCCTATCTCTCTCAAACAACAACGGAAGGAGAGGAAAAAAAGAAGAGATAGGGGGCGACGAAGCGCAAACGCGCAATTTATAGGGATTTCTCCAGCCTGATACCAGTCGACATAGCGTTGCCGCGCGCTCAGAGTTAATGTTGACGAAGAATTCCCGGGTGCAAATTACGTTGATCAGTTTTTTTTTATACAAGTAAAAAATGTTATACGCAGTGGCAAATTTCCGCCATTTACGTCACTTTATGAGCAATTCGCATAAAAAATGTAAAGTTGCATTGTAAAATCTTTTGTACAGCATAAACACAGAAACGAATAGCTGGCGAACCTGGTCTGCGGATAAGAGCGGTAATGAGCAACACACAAGCATGGTATAATTGCCCTGCCCCACCCGCTGCTCACCTGGTCGTTTAGGTTAGTTCGTTTACTTTAACCCTTTCTCGCTTGAGCTTCGCAGTTTTGGCACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGGGCGTTTTACGCCATGCCCGGGCGAGTGGCCGGGAATTGCGCTGAGCTGTTCGCTGGGAATATCGCCGCACCATCCTGCTTTTTTCCACCAGGCTCCGCTGGAACATGCCTGGCGACCAATCATATTGTTGAAGCCGTGGTCGGGGCAAGTGCATCTCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGCATGCCCGACCACGCTTAACAATATTGATGGTCGTCCAGGTCATGTTCAGCGGAGCTGGTGGAAAAAAGCAGGATTGGATGCGGCGATATTTCCCAGCGAACAGCTCAGCGCAATTCCGCCACATGCGCCGGCATGGCGATAAAACGCACCATCAATGTGGACGTCGCTGAACCGAAGTGGGGCGAGTTGGTGCCAAACTGCGGAAAGCTCAAGGAGAAAAGGGTAAAGATAAACGAACAAACTAAATGACCAGGTGAAGCACGGTGGGGCCAGGGAATTAATACAGCTGTGATTGTTTGTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTCATTACCGCTTTATCCGAAGACCAGGTCGCCATATTCGTTTCTGTGGTTTAGCTAGTACAAAGTTTTACATTTTATATGCGAATTGCTCATAAAGTGACGTAAAGGCGAATAATTTGCGCACTGCAGTATAACATTTTTTACTGTATAAAACTGATCAACGTAATTTGCCACCGGGATTCTTCGTCAACAATTTACTCTGAGCGCGCGGCACGCTATTCGACTGGTAATCAGGCGGATGAAATCCCTATAATTGTCCCCCCCGGGGGGCGTTGGCGCTCGCGCCCCCTATCCTTCTCAACAACAACGGAGAGAGGAAAAAGAGAGAGATAGGGGGGGGCGACAAGCGCCAAACGCGGCAATTATAGGGATTTCCATCCGACCTGATACCATCGAATAGCGTTTGCCGCGCGCTCAGAGTGTAATTGTTGACAGAATTCCCGGTGGCAAATTACGTTGATCAGTTTATACAAAGGTAAAAATGTTATACGCAGTTGCGCAAATTATCCGCCTTTAGTCCTTATGAGCAATTCGCATAAAACCTAAAATGTAAAACTTTAGTACTTGCAATAAACACAGAAACAATACTTGGCGACCTGGTCCCTCTGCGGATAAAGCGGTAAATGAGCAAACAACACAGCATGTATTAATTGCCCTGCCCCAACCCGCTGCTTCACTGGTCAGTTTAGGTTTAGTCTGTTTTATCTTTACCCACTTTTCTCGCTGAGCTTCCGCAGTTTGGCACCCAACTCGCCCACTGTGTTCCCGACGTCCATCATGATGGTGGCGTTTATCGCCATGCCGGGCGCATGTGGCCGGAATTTGCGCTGAGCTGTTTGCTGGGGAAATATCGCTGCATCCATCCTGCTTTTTTCCACGAAAAAAACAGCTCGCAGAAATGACCTGGACGACCATCATATTAGTTGAGCCGTGGTCGGGGCATGCATCTCTCTCAAAACAACAACGGAGGAGGAGGAAAAAAGAGAGAGATGCACTGCCCCGAACACGGCTTCAACAATATTGATGGTCGTCCAGGTCCCCATGTTCAAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGTGCGATATTTCCCATCCCCCCCCCGAGCGAACAGCTCAGGCAATTCCCGGTCCACATGCGCCCGGCATGGCGGATAAAACGCCACCCGATCATGATGGACGTTCCAGGGAACCACATGGGGCGAGTTGGTGCCAAACTGCGAAAGCTCAAGGCGAGAAAAAGGGTAAAGATAAACGGACAAACCTAAACTGTCCAGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGGCTGTATTGTTTGCCATTAGCCCCCCCCCCCCCCCCCATTCCGCTTTTATCCCGCAAGACCAGGTCGCCATATTCGTTCTGTGTTTAATCTAGTACACAAGTTTACATTTATAGCTAATGCGAATTGCTCATAAGTGACGTAAAGGCGGATAAATTCGCAACTGCGTATAACATTTTTTACTTGTATAAAATGATCAACGTAATTTGCCACCGGGAATCTTCGCAACAATAACTCTGAGCGCGCGGCACGCTATTCGACTGGTATCAGGCGGTGAAATCCCTAAATTGCGCGTTTGGCGCTTCGTCGGCCTTCCTCTCAAACAACACGAGGAGGAGGAAAAAAAAAAAAGAGAGAGATAGGGGGCGACGAAGCGCAAACGCGCAATTATAGGGATTTATCCGCCTGATACCAGTCGATAGCAGTTGCCGCGTCGCTCAGAGTTAATTGTTGACAAGAATTCCGCGGTGGCAAATTACGTGATCAGTTTTATACAAGAGTAAAAAAGTTTATACGCAGTTGCGCAAATTATCCGCCTTTACAGTCACTTATGAGCAATTCGCATATAAATGTAAAACTTTTGTACTAGCATAAACACAGAACACGAATACTGGCGACCTGGTCTTGCGGAAAAGCGGTAATGAGCAAACAAATCACAGCATGTATTATTGCCCTGCCCCACCCGCTGCTTCAAACCTGGTCAGTTTAGGTTTAGTCTCGTTTATCTTACCTTTTGCTCGTTGAGCTTTCCAGTTTGGCACCCAACTCGCCCACTGTGGTTCCCGACGTCCATCATGATGGTGGCTGTTTTATCGCCCATGCCGGGCGCATGTGGCCGGAATTGCGCTGAGCTGTTCGCTGGGAAACTATCGCCGCATCCATCCGCTTTTTTCCACCATGCTCGCTGAACATGACCTGGACGCACCATTCATATTGTTGAAGCCGTGGTGGGGCATGCATCTCTCTCAAACAACAACGGAGGAAGAGGAGAAAAAGGAGAGAGATGCACTGCCCGACCACGGCTTCAACAATATTGAGGTCGTCAGGTCATGTTCAGCGAGCTGGTGGAAAAAAGCAGGTGGCAGGATGGGCGGCCGATATTGCCGAGCGAACAGCTCAGCGCAAATTCCGGCACATGCGCCCGGCATGGCGAATAAAACGTGCACCATCATGATGGACAGCGGGAACCACAGTGGGGCGAGTTTTGGGTGGCCAACTGCGAAAGCTCAACGAGAAAAGGGTAAAGATAAACGAGACTAAACTAAACTGACAGGTGGAAGCAGCGGTGGGCAGGCAATTAATACAGCTGTGATTGTTTGCTCATTACCGCTTTTCCGCAAGACCAGGTCGCAAGTATTCGTTTCTGTGTTTATGCTAGTAACAAAAGTTTTACTTTTATATGCGAAAATTGCTCAAAAGTGACAGTAAAAAGGCGGATAATTTTGCGCACATGTCGTATAACATTTTTACCTTGTATAAAACTGACAACGTAAATTTGCCACCAGGGAATTCTTCGTCCAACAATTTAACTCTGAGCGCGCGGCAACGCTATTCGACTGTATCAGGCGGATGAAATCCCTATAATTGCCGGCTGATTGGCTGCTTCGTCCGCCCCCTACGTCCTCAAACAAAACGGAGGAGGAGGAAAAAAAGAGAGATAAGGGGCGACAAGCGCAAAACCGGCAAATTAAGGGATTTTCATCCGCCTGATACCAGTCAGGAATAGCGTTGCCGAGCGCTCAGAGTTAATTGTGACAAGAATTCCGGTGGCACAAATTATTGTCAGTTTATACAAGGTAAAAAATTTATACGCAGTTGCGCAAATTATCCGCCTATTACGTCACCTTTAGACAATTCGCATAAAAATGTAAAACTTTTGTACTAAGCATAAACACAGAAACGAATACTGGCGACTTGGTCTTGCGGATAAAGCGTAATGAGCAAACAATCACACATGTATTAATTGCCCTGCCCCAGCCGCTGCTTCACCAGGTCCAGTTTAGGTTTAGTCTCGTTATCTTTACCTTTTCTGCTTGAGCTTTCGCAGTTTGGCACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGTGCGTTTTATCGCATGCCGGGCGCATGGCCTGGAATTGCGCTGAGCTGTTCAGCTGGGAATATCGCCGCATCCACCTGCTTTTTTCCCACAGCTCGCTGAACAGACCTGGAACGACCATCATATTGTTAGAAAGCCGTGGTCGGGGCCGTGCATCTCTCTCAAACCAACAACGGAGGACAGAAAAAAAGAGAATGCACTGCCCGCCACAGGCTTCAACAATATTGATGGTCGTCCAGGTCATGTTTCAGCGAGCTGGTGGAAAAAAGCAGGATGGATGCGGCGATATTTCCCCAGCGAACAGCTCAGCGCAATTCCCGGCCAATGCGCCCGGCATGGCGATAAAACGGCCAACCATCATGATGGACGTCGGGACCCATGGGGGCGAGTTGGGTGCCAAACTGCAAGCTCAAGCAGAAGGGTAAAGATAAGATAAACGAGACTAAACCTAAACTTGACCAATGGTGAAGCAGCGGGTGGGGCAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTAGCTTTAATCCGCAAGACCAGGTCGCCAGTATTCGTTTCTGTGTTTATGCTAGTACAAAAGTTTTATGACATTTTTATATGCGAATTGCTCAAAAGTGACGTAAAGGGGATAAATTTGCGCAACTGGCGTATAACATTTTTTACTTGGTATAAAGCCTAAAACTGATCAACCGTAATTTGCCACCGGGAATTCTTCGTCAACAATTAAACTCTGAGCGCGCGGCAACGCTATTCCGACTGTAATCAGCGATGAAATCCCTATAATTGCCCGCGTTTGCCGCTTCTCGCCCTATCCTCTCAAACAACAACGGAGGAGGGGAAAAAAGAGAGAGATAGGGGCGACGAAGGCCGCCAAACGCGGCAATTATAGGGATTTTCATCCGCCTGATACCAGTCGAATAGCGTTTGCCCGCGCGCTCAGATTAATTGTTGACGAAGAATTCCCGGTTGGCAAATTACTCTGATGCATGCTTTTAGACAAGGTAAAAAGATGTTATACGCAGTTGCGCAAATATCCGCCTTTACGTCACTTTATGAGCAATTCGCATCATAAAATGTAAAACTTTGTACTAGCATAAACAAGAAACGAATACTGGCGAACCTGGTCTGCGGATAAAGCGGTAATGACAAACAATCACAGCATGTATTAATTGCCCTGCCCCACCGCTGCTTCACCGGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTTCTCGCTTGAGCTATTCGCAGTTTGGCACCCTTCTCGCCCACAGAGAAGAGGTGGTTGCCCGACGTCCCATCATGATGGTGGCTTTTATCGCCATGCCGGGCGCATGGGCCGAATAGCGCTGAGCTGTTCGCTGGGGAATATCGCCCGCATCCAATCCTGCTTTTTTCCACGAGCTCGCTGAACATGACCTGGACGACCACAATATTGTTGAAGCCGTGGGTCGGGGGCCAGTGCGATCTCTCTCAACAAAAACGGAGGAGGAGGAGAAAAAAGAGAGGAGCACTGCCCCGACCACGGCTTCACAATATTTGATGGTCGTCCAGGTCATGTCAGCAGAGCTGGTGGAAAAAAGCGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAAATTCCCGGCCAACATGCGCCCGGCAATGCGATAAACCGCCACCAATCATGATGGACGTCGGGGACCACAGTGGGGCGAGTTGGGTGCCAAAACTGCGAAAGCTCAAGCGAGAAAAGGGTAAAGATAAACGAGACAAACCTAAACTGGACCAGGTGAAGCAGCGTGGACAGGGCAATTAATACATGCTGTGATTGTTTGCTCATTACCGCTTTATCCGCAAGACAGGTCGCCAGTAATCGTTTCTTGTGTTTAGCTAGTACAAAAGTTTTTATTTTATATGCGAATTGCTCATAATGTGACTGTAAAAGCGGATAATTTGCGCAACTGCGTATAACATTTTTTACCTTGTATGAAAAACTGATCAACGTAATTTGCCAACCGGGAATCTTCGTCAACAAATAACTCTGAAGCAGCGCGCAACGGCTATTCCGACTGGTATCAGGGCGGCATGAAATCCCTATAATTGCCGCGTTTGGCGCTCGCGCCCCCTATCTCTCTCAAACAACAACGGAGGAGGAGGAAAAAAGAGAGCAGATAGGGGCGAGAAGCGCCAAACGCGGCAATTATAGGGATTTCATCCGCCTGATACCAGTGAAATAGCGTTGCGCGCGCTCAGAGTAATTTTGAGAAGAATTCCCGGTGGCAATTACGTGATCAGTTTTATACAAGGTAAAATGTTATACGCAGTTGCCAAATTATCCGCCTTTACGTCACTTTTATGGCAATTCGCATATACAAATGTAAAACTTTTGTACTAAGCATAAACACAGAAACGAATAGCTGGCCGACCTGGTCCTTTGCGGATAAGGGGAAGCGGTAAATGAGCAAACATCACAGCATGTATTAATTGGCCCTGCCCCCCGCTGCTTCACCGGTCAGTTTAGGTTTAGTCTCGTTTATCTTTACCCTTTTCTCTGAGCTTTTCGCAGTTTGGACCCAAATACTCGCCCCACTTGTGGTTCCCGACGTCCATCGATGATGGTGGCGTTTTTATCGCCAATGCCGGGCGCATGTGGCCAGGGAAATTGCGCTGAGCTGTTCGCTGGGAATATCCGCCGCCAATCCTCCTGCTTTTTTTCCACCAGCTCGCTGAACATGACCTCGGACGCCATCAATATTGTTGAAAGCCGTGGTCGGGGCAGTGCATCTCTCTCAAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATGCTCTGCCCACCACGGGCTTCAAAATATGATGGTCCGTCCAGGTCATGTTCAGCGAGCTGGTG";

        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}, strand:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score,
                hit.strand
            )
        });

        // aligner.mapopt.min_cnt = 2;
        aligner.mapopt.best_n = 10;
        // aligner.mapopt.min_dp_max = 10; // min dp score
        aligner.mapopt.pri_ratio = 0.5; // min dp score
        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();

        println!("---------------------------------------------------");
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}. strand:{:?}, qlen:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score,
                hit.strand,
                hit.query_len
            )
        });
    }

    // #[test]
    // fn test_build_aligner_v2() {
    //     let targets = vec![
    //         ReadInfo::new_fa_record(
    //             "t1".to_string(),
    //             "AAAAAGCAGGATGGATGCGGTCGATATTTCCCAGCGAACAGCTCAGCGCAAT".to_string(),
    //         ),
    //         ReadInfo::new_fa_record(
    //             "t2".to_string(),
    //             "CCCGGCCACAGCGCCCGGCAATGGCGATTAAAACGCC".to_string(),
    //         ),
    //     ];
    //     let mut aligner = build_aligner_v2(
    //         "map-ont",
    //         &IndexParams::default(),
    //         &MapParams::default(),
    //         &AlignParams::default(),
    //         &OupParams::default(),
    //         &targets,
    //     );

    //     // let aligner = &mut aligner[0];
    //     aligner.mapopt.best_n = 1;
    //     aligner.mapopt.q_occ_frac = 0.0;

    //     aligner.idxopt.k = 4;
    //     aligner.idxopt.w = 1;

    //     aligner.mapopt.min_cnt = 2;
    //     aligner.mapopt.min_dp_max = 10; // min dp score
    //     aligner.mapopt.min_chain_score = 10; // this is important for short insert
    //     aligner.mapopt.min_ksw_len = 0;

    //     let query = b"TTTCCCAGCGAACAGCTCAGCGCAATCCCGGCCACAGCGCCCGGCAA";
    //     for hit in aligner
    //         .map(query, false, false, None, None, Some(b"query"))
    //         .unwrap()
    //     {
    //         println!("{:?}", hit);
    //     }
    // }

    #[test]
    fn test_align_n() {
        let seq = b"AANNNNNNNNNNNNNAAAAAAAAANNNNCCCGTTT";
        let seq = b"AAAAAAAAA";
        let mut aligner = Aligner::builder()
            .map_ont()
            .with_cigar()
            .with_sam_hit_only()
            .with_sam_out()
            .with_seq_and_id(seq, b"ref")
            .unwrap();
        aligner.idxopt.k = 3;
        aligner.idxopt.w = 1;
        aligner.mapopt.q_occ_frac = 0.;
        aligner.mapopt.occ_dist = 0;

        aligner.mapopt.min_cnt = 2;
        aligner.mapopt.min_dp_max = 2; // min dp score
        aligner.mapopt.min_chain_score = 2; // this is important for short insert
        aligner.mapopt.min_ksw_len = 0;
        aligner.mapopt.mid_occ_frac = 0.;
        aligner.mapopt.min_mid_occ = 0;

        println!("aligner.mapopt:{:?}", aligner.mapopt);
        println!("--------------------");

        println!("aligner.idxopt:{:?}", aligner.idxopt);

        // b"AACGTCGTCGTCGTAAAAAAAAACGTGCCCGTTT",

        for hit in aligner
            .map(
                b"AAAAAAAAA",
                false,
                false,
                None,
                Some(&[67108864, 68719476736]),
                Some(b"q"),
            )
            .unwrap()
        {
            println!("hit:{hit:?}");
        }
        println!("hello");
    }

    #[test]
    fn test_bio_aign_n() {
        let scoring = Scoring {
            gap_open: -2,
            gap_extend: -1,
            match_fn: |a: u8, b: u8| {
                if a == 'N' as u8 {
                    0
                } else if a == b {
                    1i32
                } else {
                    -3i32
                }
            },
            match_scores: Some((1, -3)),
            xclip_prefix: 0,
            xclip_suffix: 0,
            yclip_prefix: 0,
            yclip_suffix: 0,
        };
        let x = b"NNNNAAAAAANNNNAAAAAAA";
        let y = b"ACGTAAAAAAGGACGGAAAAAAA";
        let mut aligner =
            bio::alignment::pairwise::Aligner::with_capacity_and_scoring(x.len(), y.len(), scoring);
        let alignment = aligner.custom(x, y);
        println!("{}", alignment.pretty(x, y, 80));
        println!("{alignment:?}");
    }

    #[test]
    fn test_bio() {
        let score = |a: u8, b: u8| if a == b { 2i32 } else { -4i32 };
        let gap_open = -4;
        let gap_extension = -2;

        let mut aligner = bio::alignment::pairwise::Aligner::new(gap_open, gap_extension, &score);

        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let alignment = aligner.local(x, y);
        println!("{:?}", alignment);
        println!("{}", alignment.pretty(x, y, 30));

        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let alignment = aligner.semiglobal(x, y);
        println!("{alignment:?}");
        println!("{}", alignment.pretty(x, y, 30));
    }

    #[test]
    fn test_align_short() {
        let seq = b"GGTAGCGTTACAAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATGCTCTGCCCACCACGGGCTTCAAAATATGATGGTCCGTCCAGGTCATGTTCAGCGAGCTGGTG";
        let qry = b"GGTAGCGTTTCAAACAACAAGGAGGAGGAGGAAAAAAGAGAAGATGCTCTGCCCACCACGGGCTTCAAAATATGATGGTCCGTCCAGGTCATGTTCAGCGAGCTGGTG";
        let mut aligner = Aligner::builder()
            .map_ont()
            .with_cigar()
            .with_sam_hit_only()
            .with_sam_out()
            .with_seq_and_id(seq, b"ref")
            .unwrap();
        aligner.idxopt.k = 3;
        aligner.idxopt.w = 1;
        aligner.mapopt.q_occ_frac = 0.;
        aligner.mapopt.occ_dist = 0;
        aligner.mapopt.a = 2;
        aligner.mapopt.b = 5;
        aligner.mapopt.q = 2;
        aligner.mapopt.q2 = 24;
        aligner.mapopt.e = 1;
        aligner.mapopt.e2 = 0;

        aligner.mapopt.min_cnt = 2;
        aligner.mapopt.min_dp_max = 1; // min dp score
        aligner.mapopt.min_chain_score = 1; // this is important for short insert
        aligner.mapopt.min_ksw_len = 0;
        aligner.mapopt.mid_occ_frac = 0.;
        aligner.mapopt.min_mid_occ = 0;

        println!("aligner.mapopt:{:?}", aligner.mapopt);
        println!("--------------------");

        println!("aligner.idxopt:{:?}", aligner.idxopt);

        // b"AACGTCGTCGTCGTAAAAAAAAACGTGCCCGTTT",

        for hit in aligner
            .map(
                qry,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]),
                Some(b"q"),
            )
            .unwrap()
        {
            let mapping_ext = MappingExt(&hit);
            let aligned = mapping_ext.aligned_2_str(seq, qry);
            println!("{}\n{}", aligned.0, aligned.1);
        }
        println!("hello");
    }

    #[test]
    fn test_align_short2() {
        let target = b"AAACAACAACGGAGGAGGAGGAAAAAA";
        let qry = b"AAAAAAAACCCCCGGGGTTTTTGAGAGAGATATCTCTCTCAAAACCCCGGGGTTTTAAACAACCAACGAAGGAGGAGGAAAAAAAAAAAAAACCCCGGGGTTTTTGAGAAGATATCCTCTCAAAACCCCGAGGGTTTTAAACAACAACGGAGGAGGAGGTAAAAAAAAAAAAAAAACCCCGGGTTTTGGAGAGAAATATCTCTCTCAAACCCCCGGGGTTTTAAACAACAACGGAGGAGGAGGAAAAAAAAAAAGAAAAAACCCACGGGTTTGATGAGAGATATCCTCTCAACCCCGGGTTTTTAACAAACACGGAGGAGGAGAAAAAAAAAAAAAAAAAAAAACCCCGGGGTTTTGGAGAGAGAATTTCTTCTCTCACAACCCCGGGGTTTTAAACACAACGAGGAGGAGGAAAAAAAAAAAACCCCGGGGTTTTGGAGAGAGAATATCTCTCTCAAACCCCGGGGTTTTTAAACAACAACGGAGGAGGAGGAAAAAAAAAACCCCGGGGTTTTGAAGAGATATCTCTCTCAAAAACCCCGGGTTTTAAAACAACAACGGAGGGGAGGAAACAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAGAAAAAAAAAAAACCCGGGAGTTTTTGAGAGATATCTCTCTCAAAAACCCCGGGGTTTTAAACAACAACGGAGGGAGGAGGTGAAGAGGAGGAGGAAAAAAAGGAGGAGAAAAAAAAAAAAAACCCCGGGGTTTTGAGAGAGATATCTCTCTCCAAAACCCTGGTTTTTAAACAACAAACGGAGGAGGAGGAAAAAAAAAACCCGGGGGTTTTGAGAGCAGATATCTCTCTCAAAACCCCGGGGTTTAAACAAAAACAACAACGGAAGGAGGGGAAAAAAAAAAAAAAAAAAAAAGAGGAAAAAAAAAAGCAAAAAAAAACCCCGGGGTTTTGAGAGAGAGAAGATTCTCTCTCAAAACCCCGGGGTTTTAAACAACTACGAGGAGGAGGAAAAAAAAAAAATCCCAAAAAAAAACACCGGGGTTTTGAGAGAGATATCTCTCTCAAACCCCGGGGTTTTAAACAACAGACGGAGGAGGAGGAAAAAAAAAAATAAAAAAAAAAAGGGCAAAAAAAACCCCGGGGTTTTGAGAGAGAATATCTCTCTCAAAACGCCCGGGTTTTAAAACAACAACAGGAGGAGGGAAAAAGAAAAGGAGGAAAAAAAAAACCCCCCGGGGTTTTGAGAAATATCTCCTCAAAACCCCGGGGTTTTAAACAACAACGAGGAGGAGGAAAAAAAAAACCCCGGGGTTTTGAGAGGATATCTCTCTCAAAATATCAAAACCCGGGGTTTAAACAACACGGATAGGACGAGGAAAAAAAAAACCCTCGGGGTTTTGAGAGAGATATCTCTCTCAAACCCCGGGGTTTTAAACAACAACGGAGGAGGAGGGAAAAAGTAAAAAAACCCCGGGGTTTTGAGAGAGAAGAAGAGATATCACTCTCAAAAACCCCGGAGGTTTTACCAACAACAA";

        let mut aligner = Aligner::builder()
            .sr()
            .with_cigar()
            .with_sam_hit_only()
            .with_sam_out();
        aligner.idxopt.k = 5;
        aligner.idxopt.w = 3;
        let mut aligner = aligner.with_seq_and_id(target, b"ref").unwrap();

        // aligner.mapopt.q_occ_frac = 0.;
        // aligner.mapopt.occ_dist = 0;

        aligner.mapopt.min_cnt = 3;
        aligner.mapopt.min_dp_max = 30; // min dp score
        aligner.mapopt.min_chain_score = 5; // this is important for short insert
        aligner.mapopt.min_ksw_len = 0;
        aligner.mapopt.mid_occ_frac = 0.;
        aligner.mapopt.min_mid_occ = 0;

        println!("aligner.mapopt:{:?}", aligner.mapopt);
        println!("--------------------");

        println!("aligner.idxopt:{:?}", aligner.idxopt);

        // b"AACGTCGTCGTCGTAAAAAAAAACGTGCCCGTTT",

        for hit in aligner
            .map(
                qry,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]),
                Some(b"q"),
            )
            .unwrap()
        {
            let mapping_ext = MappingExt(&hit);
            let aligned = mapping_ext.aligned_2_str(target, qry);
            println!("{}\n{}", aligned.0, aligned.1);
            println!("hit:{hit:?}");
        }
        println!("hello");
    }

    #[test]
    fn test_align() {
        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let target = b"ACAAAAAAATTAAAATAAAGAATTCCCGGCTTGGGGGCCAGAGTCCTCACCTCAT
ATGTAATAAGATGTGTGGATCTGTGACAGTTTCACATCCTCCTGCCTCACTGTCTCCATATAAAATAAAACCGTTTCCTGGATGACAGAGCATGAGTCATAAACAGATGGGAAAGCCATGCAAAAGT
AAGAGGACTCTTCGTTTTTATGACTGACGGATCCTCCGCAAAAATACTTGCACGTCTCTCATGTGTGTCAGGTTCTCTATTGCAGGCCTGAATGCATACGATTAACTTGCTTATCATTTGGAAGTAAAGAGAAACACTTTCTGTCTTAGAACTTTTTGTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAGAGAGAGAGAAGAGAGAGAGAGAGAGAGAGAGAGAAAGAGAGAGAGAAAGAGAGAGAGAGAGAGAGAGAGAGGAGAGAGAGAGAGAGAGATTTTGCTATTCTTTGCATAGCATTTGGTTATGGCCTATTATAAGACTCATTATAGTCAGGGAAATTATGGGGTGAGAGAAACACTCTTGTGGAAAGATGAGCTGAGGTTTATATAAAGAAATTGACTAGTTTATGGAGGCTTATAAGTTTTTAGAAGCTACTCCCATTTGAGACTGACTTAGCAGTGATTTTAAACTTTAACATCAAGAAGGAAGAGGAGGAGAAGGAAGAGGAAGAAGAGGAAGAAGAGAAGAGGAAAAAGAGGAAGAAGAGGAAGAGGAAGAAGAGGAAAAAGGAGAAAGACAGACACCAAATAATCAACTGTATTTCTGTTTGTATCCCTTGACTGTCTTCCAAGTTTCTGCTTCACTTTTGTTTCCATAGCAACATGCAGAGATAAAACAGGTAGTAAATATTAAGCAAAAGTTGAATGAATAAAAGGAGAATAAAAAGCCAAGGAAAATGCTTAAAAATCTTACCTGTGGCTCAAGCATCCAACTACCAGTTCCTAGAAGATAGCTTACAGAAGATGGTAACTGTATCCATCAGGATACTGTCAGCAAAATCCTAGGGAAAAATGTTTTTAGACAAACCAACTACATTTCTTGAACAAAAGTAATTGCAAGTGGACAAGAGCTAGAGGAGGATTGGAGAACAAAGGAGACACAGGAGACATAGACATCAGTTGCAAAGAGTACACTTTATCTGGATCCTGGTTTGTTGTTGTTGTTTGTTTGTTTTGGTTTTTTGTATGTTTTGATTATTTTATTTTATGTGCATTGGTGTTTTGCCTGTGTGTATGTATATGAACCATGGGATGACTGGTGCTGGCAAAAGTCAGAAGAGGGTCTCATCTCCCCCAGAACTGAGGTAACAGATGGTTGTGGGACCCTAACCTGGGTGCTCTGCAAGAGCAGAAAGTGCTCTCAACTGCTGAGCCAATCTATTCTCTCAGCCCCAGAATTCTGATTTTTTAAAAAAATTATGATATCATGAGACATTTAGAAATCTGAACATTGCCAAATTTGGTGGCATGAGTTATCAATTGTTTTGCTTGTAATGTTTACAAATAAAATGATATCATATCTAGACTTTGCCTCAAAGTCATACGGGAGGAGGACTTGGGCAACAATAGGATGGGCTGAGGTAGGGACTGAGTCTGGAGAGGAGTTCCATGCTCACACATAAGGCATCACCACAGTGGGCAAGTTAATGTGTGCTTTGAGCTTCTTACTGAAGTGTGGCAAACAGAAAGCCCTGAAAAGCCCCTGGTGCACAACAGGAACCATTATACAAATAAAAAGATGGGTTCTCCCTCAGTGGGAAAAGGTTTTCGTTTTACATAAGGATATGTGTGATTTCTTACCTGTGCTCAAATTACCAAGCAATCAGGAACTTACAATCTTTTTCAGAGTAGCTACCCTCTGAGGAGATGTATTGTTTTTTAACACTAGTTGTTGCAGGCAGTCTTTCAATGCAAAGCAAACTGGGAAAATGGTCTCCTGGCACAATCGCCAACATTAGTAATAGGTAGTGAAGCTTGTTCACGTCTACAAAGTTCAAAGTCCTACAAGTTAAGCTGAGTGTCCTAAGGTGCTTATCACTGCAGTTAAACAAGGCGTGGTAACGGAATGTGCACAGCACTCTGGGAAAAATTAGTTATGTCCGTTTGCACAATGTTGCTTACTTAGAACGTTGACTTGCCTCAACCTCACACTGTGTCCAAGCTTTTAACTGACAACTGTTTTCCCAAAACCCAAAGATTATGTGTAAAAGAACAATCCTTACTTGTATCTGACTTGCCCTCAGAAACCCAGGTAAGGCTGTGTTCTCCAGCCCAAGCCCTCCCGTTGTAGACTCTGAATGGTTGGCACAAAGCATGGGTTTAAATGACCAAATGGTGTTTTTATGTTACAGTTTAGTTTGTAAGTTATTTGGAAAGCACTTCTTACCAGTAACCAGCCGCTTTGCCAGCTAACAGCTAAGTGCTTTTCAGAGATTCCTAATGTCTGTCAAGAAGAGTGTGGCTGCTCTGCTTCACTAATGAGGAAAGCAAGGCTCAGAGAAATGGGTGGCTCACCAGCCTGTAATGGCAGGCAGAGCAGGCAGTCAGGAGTGGGCTGTCTGGATTTATTGTTTTCTCCTAGTGTGCACGTGTGCATACGTGTGTGTGTGTGTGTGCATACGTGTGTGTGTGTGTGTGTGTGTTCATGGGTGGATAGATACAAAATGTGTGCATGTGAAGAGCAAAACTTCAGATGCCATTCCTTAGTCACTATCCTAACTTTTTCTCTCTCTTGAAGTAGGGTTTCTCACTGGACTGGGCAGTTAGCCACTACAGCTGGACCAACTGACTATGGAACCCACACAACCCTCTATACATCATCCCACACACAACTGGCATTACACCTTGCCATCATGCTATCATGCTGTTCCCCCTCCCTATTTCATATAAATACATCCTTCCTTCCCTCCCTAATACTTCCTCCTCCTTTCACATTTGTTCTGGGTACCAGACTCAGGTCCCCACACTTGCAAGGCAAGCACTTTCTGGCTGAGCCATCTCTCCAGCCATTATCCTATAGCATTAATCATTTTCTCATTGACATTGTAAGCCCCATCAGAGAACTCGTTCATTAGTCTCCAGGGCTCCTTAACCCTAGTAGAGTTATATAAACAATGTCGCTCTCAACAAATCGTTCAAGGGCTGATGTCCTGCAACTTGTCTAACCTTCTCAAGCCATTGCTGCTTACCAGTGTTTTCCAGGCTCATATTTTGAAACTGACCAAGTGTGATTATCCAGATAACTATTTAAGGGAAATGTTTGAGAACAACGTCTGCATCTGGAGCGAGAGAATTACACTTAACAAGGGAAAGAGGAATGAAAACCAAACACTGTAAGCAAGGTGAGGGAAACTGAGATTGAAATATAAGGACCAATCGGTCAGGTTGATGCAGTTTCCCTGTGTAGGCCCTAGAGTAAGTGGAAAACCCAACCAACAATGGAGAAAATTCTAGAACACTGGAACCAAGTGGAACACACACACACACACGGGATGGTTGTAGAAGATAGATACTCTAGATATAAAAAGGCAGCCACACTCTCCCTGTTCCCCTGTCCCCAGCCTGTTGGGACCTGTTTGCCTTTCCACTGCTTAGAGGAACATTTTTACTGGTGTGGTCTTATCCGGATTGAGAGTGATGAGGTTCCAGTCAAACAATTGGACAAGCACGGTGCTAAACCATAGCTGATGGCTATCCCATGCCCACGAACGTTATATGTGTTATCTAATTTAATTCTGATACACACAAGCAAAGGACGGAGATGCAGAGACTGTGGCTTTGCAGTTACAAGGTACTTATTAAATTCTAGGCTCAGGTTTAGCACTCTTGCTGGGTTTTGCTAGGGGAGGATGGCCAGTAGGGTTGGGAGAGTGCTACACCACATCTTTCCCCCAGCACTCGGGTGATAAAGGTGGCAAGTGTGTAGAGTGTTCTGCAGGAGCCAGCTCTAGTCATGTCTAGAGAAAACATTTGGGGATTGCTGCTTTGTCTACCATCTGATTCTCTGTGTCTGGATTCTGTAGGACTCGCCAGCTCATTTATCCTGGTTCCAGTTTATGCACCCAGCCTGGTTGACTTTTCTGGAACATCGAGGTTACTGTAAACTCTAAGGGCAAGGGCCTAGGCTCCTGGCAAGTTTCCAGTATGAGATCTGCTAAGTACAAACTATTGGTGGGTTCATGATTATGGTCTCTCACTGAGAGAGCGTGATTATACACTACTATATCATAACATTCAGGAAGATGACCAGGGCATTGGGAGTTCCAGGGCAGCCTGGGCTACATAGCCAACTCTTGTCAAAAGAGAAAGATTTGAAAATTATAAATCAATTTGAAATTATTTTGTTACTAGTAGAAGTGAACTAAATTCTGGAATCTTTAAGCATGTATAATTTTAAAATCTATGCTTTATGTGTTCTTACCAAATGCAATTTTTTTCCTGTAGGTTTTCTTAAAAGGTGTTGGCATTTTTGTAAGAAGTGATTGTTAATAATACCTGGAGTTTGACTTGACTTTTTAAATAGCATACTTGTCTTCCTAAAACAAAACTGGTCCATCTTTTGTTAAAAAGAATAAAGTTAATTTCTGAAGGTTGTCCAATAATGTATCTACTTGTTAATTTTAAAAAGGGTGTCTTGACTCAGGCTGCTATTACAGAATTCCATGTTAGATGTTTTGTATATAATAGAAATTTATATCTCATACTTCTGGAGTCCAACCTACTAAAATGTTCACTTTATACTTTATCAAGTTAGCAGTGAATATATTTTTAACAGCTGATGTTCCCCCATGCCCCCAGCTTGCTTCCACACCTTGGATATTTTGTGGGCGAGAAGTCTGAAGCCAACATTTCAGAAGGGCTCCCTATTTAGCACAATGTCTTTACACAAGCATACAGGTGCTCTCATATGATAATATTAGCAGGCCTGGCAGGAATAGCAAAGACCAATCTTTATATAGTATTGTTTAGGCATTAAGAACTGTTTTAGGTATTTTCTGTGTGTTAGACCATTTAATCCTTTCAACAGCTCAAGTAAAGTATGAAACGGTCATCCTCAGCACTTCATTACAGAGGAACCGAGGAACGGAGCCTGAGAGGGTTTGCCTGAGACTGTGTAGCTAGTGGAAGAAAGATGAGCATTCTCAATGGTGATTACCAACTTTTGTCCCCTTTGAATGAAATTCGCCTTTAAACTATGGTGCTCAGGTCGCACCTGTAATTGGCCTGGAACATTTTGTTTTATTTTATTTTATTTCATTTTATTTTATTTATTTATAGTCTCAATGTCATGATATCCTGCCTCTGCCCCTGAGTGTCTGTCAGGCTATAGGTATAGAACACTATATCCACATTGAAAAGGTTTTTAATGCTCCTTCGGTGCGCTGAATGACTGTTGGGAATCGGATGCCTTAGCTTTTCTACCACAGACCTGTAAGTGGCTCTTGAGCTGACACCACACCTAGTTAAAGGAACAGACAACAGAAGGGGGATTCTTGGCAATAAAATGGTAAAAGCTGAATGGACATTGAGCTGTCTACAGGTTCTTGGGCAACTGACTCAGTTGCCATTTCTCCATCTGCAAAAACAGAGTTTAAAGTAAAAGATCTGCTGGGCTTTCTAGCTCTGTGAGTCCATGAACTTACTTCTATAAAATTCTCTGTAATCCGTGAGTAGATGACCCTTCAGGAGTTGTTTCATGTGGGAACCAGTTTCCTTTGGCCAACCTTCTTATCAACTGGTGTGTGCCTAAGGAAGGCAGATTGGTTTTTACCTTTGCTTTGGGTGGAACTCAATGAGGGAGGAGGAAACTGTGTGCGCGCACACAAATCATCCTGAAGGTAAGGGGACTGATTTGGAGGGATTATCAGCTGGTTCAGATTCCCCACGTCACTGCCCCCTCCAATGGCTGCGGTAATTGAAAGCTGACCTGCGTTCTGATTGCTGACAGCTGAAGGATTCCTGTTGCCAGGGAACCTCTGGCCGTGTGCCAAAGTAACGACCTGATTAATTATAATAATCCTTGTGCTGTTCCTGCCCCCATCGGCTCCCTCCTCTGAGTGGCGTCCAGCCTCAGGATCTCAGAATGGAACAGTAAGCACCCACTAAATGAACAATCCTGTGTCCTGTGGCAGAACACATGTCAAAAGACTTTCGCACCCCTCCGACCGCTCTGCTCAGTGACATACCTGTTGAAAACAATCATAAATTCAAAGTGGGTCTGACCAAGCCACACCGTGTTTCCACATAGAAGGGTCCAGAAACCTGGCCACACCGTGTTTCCACATAGAAGGGGTCCAGAAACCTGGCCACACCGTGTTTCCATATAGAAGGGTCCAGAAACCTGGCACACCGTGTTTCCACATAGAAGGGTCCAGAACCTGGCCACACCGTGTTTCCACATAGAAGGGTCCCAGAAACCTGGCCACACCGTGTTTCCACATAGAAGGGTCCAGAAACCTGGCCACACCGTGTTTCCACATAGAAGGGTCCAGAAACCTGTTCTTCTCTCCTCCACTCTGACGCATCTTGAGATTACCCCAAATCTGGAGACTCCACGATAAAAACAACATATTCTGCCATGTCCCAAAAGCTTCTGGTAACTTGGGACTGACTGACTTAAACTAGTAATCGTGCTTGGGATTTTAGTTTGGAGGGTTACCTAGAATTCTGACACCCTTCTGAAATTCAGGTACCAAGTGACAAGATTTCTGGACCCTTCTATGTGGAAACACGGTGTGGCCAGGTTCTGGACCCTTCTATGTGGAAACACGGTGTGGCCAGGTTTCTGGACCCTTCTATGTGGAAACACGGTGTGGCCAGGTTTCTGGACCCTTCTATATGGAAACACGGTGTGGCCAGGTTTCTGGACCCTTCTATGTGGAAACACGGTGTGGCCAGGTTTCTGGACCCTTCTATGTGGAAACACGGGTGTGGCTGGTCAGACCCAATTTGAATTTATGATTGTTTTTAACAGGTATGTCACTGAGCAGAGCGGTCGGAGGGGGTGCGAAAGTCTTTTGACATGTGTTCCAACAGGACACAGGGATTGTTCATTTAGTGGGTGCTTACTGTTCCATTCTGAGATCCTGAGGCTGGACGCACTCAGAGGAGGGAGCCGATGGGGCAGGAACAGCACAAAGGAGGGATAATTAATCAGGTCGTTACTTGGCACACGGCCAGAGGTTCCCTGGCAACAGGAATCTTAGCTGTCAGCATCAGAACGCAGGCAGCTTTCAATTACCGCAGCCATTGGAGGGGGCAGTGACGTGGGGAATCTGAACCAGCTGATAATCCCTCCAAATCAGTCCCCTTACCTTCAGGATGATTTGTGTGCGCGCACACAGTTTCCTCCTCCCTCATTGAGTTCCACCCAAAGCAAAGTAAAAACCAATCTGCCTTCCTTAGGCACACACCAGTTGATAAGAAGGTTGGCCAAAGGAAACTGGTTCCCACATGAAAAAATCTGAAGGGTCATCTACTCACGGATTACAGAGAATTTTATAGAAGTAAGTTCATGGACTCACAGAGCTAGAAAGCGCAGCAGATCTTTTACTTAAACTCTGTTTTTGCAGATGGAGAAATGGCAACTGAGTCAGTTGCCCAGAACCAAATTGTAGACAGCTCAATGTCCATTCAGCTTTTACCATTTTATTGCCAAGAATCCCCCTTCTGTTGTCTGTTCCTTTAACTAGGTGTGTGTCAGCTCAAGAGCCACTACAGGTCTGTGGTAGAAAAGCTAAGGCATCCGATTCCCAACAGTCATTCAGCGCACCGAAGGAGCATTAAAAACCTTTTCAATGTGGATATAGTGTTCTATACCTATAGCTTGACAGACACTCAGGGGCAGAGGCAGGAGGATCATGACATTGAGACTATAAATAAAAAATAAAATAAATAAAATAGAAATAAACAAAATGTTCCAGGCCAATTACAGGTGCGACCTGAGCACCAGAGTTTAAGGCGAATTTCATTCAAAGGGGACAAAGTTGGGAATCACCATTGAGAATGCTCATCTTTTCTTCCACTAGCTACACAGTCTCAGCAAACCATCTCAGGCTCCGTTCCTCGTTCTCTGTAATGAAGTGCTGAGGATGACCGTTTCATACTTTACTTGAGCTGTTGAAAGGATTAAATGGTCTAACACACAGAAAATACCTAAAACAGTTCTTAATGCCTAAACAATACTATATAAAGATGGTCTGTTGCTATTCCTGCCAGGCCTGCTAATATTATCATATGAGAGCACCTGGTATGCGTGTGTAAAGACATTGTGCTAAATAGGGAGCCCTTCTGAAATGTTGGCTTCAGACTTCTCGCCCAAAAAATATCCAAGGTGTGGAAGCAAGCTGGGGGCAGGGGGAACATCAGCTGTTAAAAATAGATTCACTGCTAACTTGATAAAGGAGAAAGGGAACATTTTAGTAGGTTGGACTCCAGAAGTATGAGAGATAAATTTCTATTATATACAAAACATCTAACATGGAATTCTGTAATAGCAGCCTGAGTCAAGACACCCGTTTTTAAAATTAACAAGTAGATACATTATTGGACAACCTTCAGAAATTAACTTTATTCTTTTTAACAAAAGAGGGACCAGTTTTGTTTTAGGAAGACAAGTATGCTATTTAAAAAGTCAAGTCAAACTCCAGGGATTATTAACAAGCACTTCTTACAAAAATGCCAACACCTTTTAGAAAACCTACAGGAAAAAAATGTGCATTTGGTAAGAACACATAAAGCATAGATTTTAAAATTATACATGATTAAAGATTCCAGAATTTAGTTCACTTCTACTAGTAACAAAATAATTTCAAATTGATTTATAATTTTCAAATCTTTCTCTTTTTGACAGAGTGGCTATGTAGCCAGGCTGCCCTGGAACTCCCCATGCCCTGGTCATCTTCCTGAATTTAGGAATATAGAGTGTATAATCACGCTCTCTCAGTGAGAGACCATAATCATGAACCCACCAATGAGTTTGTACTTAGCAGATCTCCATACTGGAAACTTGCCAGGAGCCTAGGCCCTTGCCCTTAGAGTTTACAGTAACCTCGATGTTCCAGAAAAGTCAACCAGGCTGGGTGTAAATATACTGGAACCAGGATAAATGAGCTGGCGAGTCCTACAGAATCCAGACACAGAGAAGCAGATGGTAGACAAAGCAGCAATCCCCAAATGTTTTCTCTAGACATGACTAGAGCTGGCTCTGCAGAACACTCTACACACTTGCCACCTTTATCACCCGAGTGCTGGGGGAAAGATGTGTGGTAGCACTCTCCCAACCCTACTGGCCATCCTCCCTAGCAAAACCCAGCAAGAGTGCTAAACCTAGCCTAGAATTTAATAAGTAACCTTGTAACTGCAAAGCCACAGTCTCTGCATCTCCGTCCCTTTGCTTGTGTGTATCAGAATTAATTAGATAACACATATAACGTTCGTGGCATGGGATAGCCATCAGCTATGGTTTAGCACCGTGCTTGTCCAATGTTTGACTGGAACCTCATCACTCTCAATCCGGAATAAGACCACACCAGGAAAATGTTTCTTCTAAGGCAGTGGAAAGGCAAACAGGTCCCAACAGGCTGGGACAGGGGAACAGGGAGAGGTGTGGCTGCCTTTTTATATCTAGAGAGTCTATCTTCTACAACCATCCCGTGTGTGTGTGTGTGTTCCACTTGGTTCCAGTGTTCTAGAATTTTCTCCATTGTTGGTTGGGTTTTCCACTTACTTCTAGGGCCTACACAGGGAAACTGCATCAACCTGACCAGATTGGTCCTTATATTTCATCTCAGTTTCCCTCACCTTGTCTTACAGTGTTTGGTTTCATTCCTCTTTACCTTGTTAAGTGTAATTCTCTGCCTCCAGATGCAGACGTGTTCTCAAACATTTCCTTAAATAGTTATCGGATAATCACACTTGGTCAGTTTCAAAATATGAGCTGAAAACACTGGTAAGCAGCATGGCTTGAGAAGGTTAGACAAGTTGCAGGACATCAGCCTTGAACGATTTGTTGAGAGCGACAGTGTTTATATAACTCTACTAGGGTTAAGGAGCCCTGGAGACTAATGAACGAGTTCTCTGATGGGGCTTACAATGTCAATGAGAAAATGATTAATGCTATAGGATAATGGCTGGGAGAGATGGCTCAGCCAGAAAAGTGCTTGCCTTGCAAGTGTGGGGACCTGAGTCTGGTACCCAGAACAAATGGTGAAAGGAAGGAGGAAGGAAGGGAGGGAGGGAGGGAGGGAAGGAAGGAGGGAGGAGGGAGGAAATAGGGAGGGGGACAGCATGATAGCATGATGGCAAGGTGTAATGCCAGTTGTGTGTGGGATGATGTATAGAGGGTTGTGTGGGTTCCATAGTCAGTTGGTCCAGCTGTAGTGGCTAACTGCCAGTCCAGTGAGAAACCCTACTTCAAGAGAGAGAAAAGGTGGATAGTGACTAAGGAATGGCAATCTGAAGTTTTTGCATCTTCACATGCACACATTTGTATCTATCCACCCATGAACACACACACACACACACACGTATGCACACACAACACACACGTATGCACACGTGCACACTAGGAGAAAACAATAAATTCCAGACAGCCCACTCCTGACTGCCTGCTCTGCCTGCCACTTACAGGCTGGTGAGCCACCCATTTCTCTGAGCCTTGCTTTCCTCATTAGTGAAGCAAGAGCAGCCACACTCTTCTTGACAGACATTAGGAATCTCTGAAAAGCACTTAGCTGTTAGATGGCAAGCGGCTGGTTACTGGTAAGAAGTGCTTTCCAAATAACTTACAAACTAAACTGTAACATAAAACACCATTTGGTCATTTAAACCCATGCTTTGTGCCAACCAATTCAGAGTCTACAACGGAGGGCTTGGGCGAGAGAACACAGCCTTACCTGGGTTTCTGAGGGCAAGTCAGATACAAGTAAGGATTGTTCTTTTACACATAACTTTGGGTTTTGGGAAAACAGTTGTCAGTTAAAAGATTGGACACAGTGTGAGGTTGAGGCAAGTCAACGTTCTAAGTAAGCAACATTGTGCAAACGGAACATAACTAATTTTTCCCGAGTGCTTGCACATTCCGTTACCCACGCCTTGTTTTAACTGCAGTGATAAGCACCTTAGGACACTCAGCTTAACTTGTAGGACTTTGAACTTTGTAGACGTGAACAAGCTTCACTACCTATTACTAATGTTGGCGATTGTGCCAGGAGACCATTTTCCCAGTTGCTTTCATTGAAAGATGCCTGCAACAACTAGTGTTAAAAACAATACATCTCCTCAGAGGGTAGGCTACTCTGAAAAAGATGTAAGTTCCTGATTGGCTTGGTAATTTGAGCACAGGTAAGAAATCACACATATCCTTATGTAAAACGAAAACCTTTTCCCACTGAGGGGAGAACCCATCTTTTATTTGTATAATGGTTCTGTTGTGCACCAGGGGCTTTTTCAGGGCTTTCTGTTTGCCACACTTCAGTAAGAAGCTCAAAGCACACATTAACTTGCCCACTGTGGTGATGCCTTATGTGTGAGCATGGAACTCCTCTCCAGACTCAGTCCCTACCTCAGCCATCTATTGTTGCCCAAGTCCTATCCCGTATGACTTTGAGGCAAAGTCTAGATATGATATCATTTTATTTGTAAACATTACAGCAAACAATTGATAACTCATGCCACCAAATTGGCAATGTTCAGATTTATAAATGTCTCATGATATCATAATTTTTTTAAAAAATCAGAATTCTGGGGCTGGAGAGATAGATTGGCTCAGCAGTTGAGAGCACTTTCTGCTCTTGCAGAGCACCAGGTTAGGGTCCCACAACCATCTGTTACCTCAGTTCTGGGGGAGATGAGACCCTCTTCTGACTTTTGCCAGCACCAGGCATCCCATGGTTCATATACATACACACAGGCAAAACACCAATGCACATAAAAATAAATAATCAAAACATACAAAAAACCAAAACAAACAAACAACACAACAAACCAGGATCCAGATAAAGTGTACTCTTTGCAACTGATGTATATGTCTCCTGTGTATACTTTGTTCTCCAATCCTCCTCTAGCTCTTGTCCACTTGCAATTACTTTGTTCAAGAAATGTAGTTGCTTGTCTAAAAACATTTTCCCTAGATTTTGCTGACAGTATCCTGATGGATACAGTTAACATCTTCTGTAGCTATCTTCTAGGAATGGTAGTTGGATGCTTGAGCCACAGGTAAGATTTTAAGCATTTCCTTGGCTTTTTATTATAATTTTTTCATTCAACTTTTGCTTAATATTTACTACATGTTTATCTCTGCATGTTGCTATGGAAACAAAAGTGAAGCAGAAACTTGGAAGACAGTCAAAGGGATACAAACAGAAATACAGTTGATGTATTTGGTGTATGTATTTATACTTAATTTAATATTATTCCTCTTTTTATTTTTCCTCTTCTTCCTCTTATTATTATTCCTCCTTCTTGATGTTAAAGTTTAAAATCACTGCTAAGGCAGTCTCAAATGGGAGTAGCTTCTAAAAACTTATAAGCCTCATAAAACTAGTCAATTTATTTATATAAACATCAGCTCATCTTTTCCACAAGAGTGTTTATCTCACCCATAATTTCCTGACTATAATGAGTCTTATAATAGGCCATAACCAAATGCTATGCAAAGAATAGCAAAATCTCTTTCTCTCTCTCTCTATATATATCTCTCTATATCTTATATTATTTCTATATATCTATCTCTATCTCTCTATTATCTCTCTCACACACACACACAACACACACAACACACACAACAAAAAGTTCTAAGACAGAAAGTTGTTTCTCTTTACTTCCAAATGATAAGCAAGTTAATCGTATGCATTCAGGCCTGCAATAGAGAACCTGACACACATGAGAGACGTGCAAGTATTTTGCGGAGGATCCGTCAGTCATAAAAACGAAGAGTCTCTTACTTTTGCATGGCTTTCCATCTGTTTATGACTCATGCTCTGTCATCAGGAAACGGTTTTATTTATATGGAGACAGTGAGCAGGAGGATGTGAAACTGTCCACAGATCCACACAGCTTATTACATATGAGGTGAGGCCTCTTGCCCCCAAGCCGGGAATTCTTTATTTAATGTT";

        let aligner = Aligner::builder()
            .map_ont()
            .with_index_threads(4)
            .with_cigar()
            .with_sam_out()
            .with_sam_hit_only()
            .with_seq(target)
            .unwrap();

        // aligner.mapopt.min_cnt = 2;
        // aligner.mapopt.best_n = 500;
        // aligner.mapopt.min_dp_max = 10; // min dp score
        // aligner.mapopt.pri_ratio = 0.5; // min dp score
        println!("{:?}", aligner.mapopt);
        println!("{:?}", aligner.idxopt);

        // tot len
        let seq = b"CAAAGAGGAGAGTTAAACATAAACCAAAATTAAATAAGAATCACGGCTTTGGGGGCCAAGAGTCCTTCCACCTCATATGTAATAAGATGTGTGGATGCTGTGGACATTTCACCACCTCCTGCCCACTGTCACATAATAAAATAAACCGTTTCCTGGATGACAGAGACATAGTATAAACAGATGAAAGCCATGCAAAAGTAAGAGACTCTTCGTTTTTATGACTGACGGATCATCCGCAAAAATACTTGCAACGTCTCTCATTTGTCAGGTCTCTATTGCAGGACCTGAATGCATAAAGATAGACTTGCGTTATCGTTTGGAAGTAAAGAGAAACACTTTTCGGTCTTAGAACTTTTTGTTGTGTGTGTGTGTGTGGTGTGTGTGTGTGTGTGAGAAGAGAGAAGAGAGAGAGAAGAGAGAGAGAGAGAAAGAGAGAGAGACGAGGAGAGGAAGAGAGAGAGGAGAGAGAGAGAGAGAGATTTTGCTATTCCCTTTGCATAGCATTTGGTTATGGCCTATTATAAGAGTATATTAATAGTCAGGGAAATTATGGGGTGAGAGAAACATAATTGTGGAAAGATGAGCTGAGGTTTATAAAGAAATTGACTAGTTTATGGAGGCTATAAAAGTTTTAGAAAATAATACCATTGAGACTGACTTAGCAAGTGATTTTAAGATTTTAAACTTTAACATCAAGAAGGAGAGGAGGGAAGAGGAAGAAGAGGAAGAAAGAAGAAGAGAGAACAAGAGGAAAAGAGGAAGAGGTAGAGAGCGGAAAACGGAGAAAGATGAAAGGAGAAAGACAAGACAAAAAATAATTCCAACTGTTTTCTGTTGTATCCCTTGACTGTATTCCAAGTTTCTGCTTCACTTTTGTTTCCATAGCAACATGCAGAGTAAAACAGGAAGTAAATATTAAGCAAAAGTTGAATGATAAAAAAGAATAAAAAGCCAAGGAAAATGCTTAAAAATCTACTGTGCTCAAGCATCCAACTACCAGTTCCTATGAGATAGCTTACAGAAGATGGTAAATGTATCCATCAGGAAAATGGTTAACAGCAAAAATCCTAGGGGAAAAATGTTTTTAGACAAGAAAATACATTTCTTGACAAAAGTAATTGGCAGTGGCAAGAGATAGAGGAGGATTGGAGAACAAAGGAGACACGAGACATAGATGACATATCAGTTGCAAAGAGACAATTTATCTGGATCCTGGTTTGTTGTTGTTGTTTGTTTGTTTTGGTTTTTTGTATGTTTTGATTATTTTATTTTGATGTGCATTGGGTGGTTTTGCCTGTGGTGTATGTATATGAACCCATGGGATGAATGGTGGTGGCAAAAAGTCAGAAGAGGGTCCTCATCTCCCCCAGAAACTGAGGTAACAGATGGTTGTGGGCCCTAACCTGGGTGCTCTGTCAAGAAGCAGAAAGTGCTCTCAACGCGGAGCCAATCTATTATCCAGCCCCAGACATTCTGTTTTTAAAAAATTATGATATCATTGAGGACATTTAGAAATCTGAACATTGCCAAATTTGGTGGCATGAGTTACAATTGTTTTGCTTGTAATGTTGTAAAAATAAATATTACCATATCTAGACTTTGCATCAAGTCAACGGGAGGAGGACTTGGGCAACATAGGATGGGCTGAGGGTAGGGCTGAGTCTGGAGAGAGTTTCCATGATCACCATAAGGAATCACCCAGTGGGCAAGTTACATGTGTGTATTTGAGCTTCTTACTGAAGTGTGGCAAACAGAAAGCCTGGAAAAGCCCCTGGTGCACAACAGGAAACCATATACAAATAAAAAGAGGGTTCTCCCTCAGTGGGAAAAGGTTTTCGTTTTTTTACATAAGGTATGTGTGATTTGTTTACCTGTGCTCAATTACCAGCAATCAGGAATCTTCAATCTTTTTTCAGAGTAGATTACTCTGAGGAGATGTATTGTTTTTTAACATAGTTGTTGCAGGCAGTCTTTCAATGCAAAGCAAACTGGAAAATGGTATCCGCACAATCGCCAACATAGTAATAGGTAGTGAAGCTTGTTTCGACGTCTACAAAGTTCAAAGTCCTACAAGTTAAGCTGAGTGTATAAGTGCTTACACTGCAGTTAAACATGGGCGTGGGTAACGGAATGTGCACACACTCTGGGAAAAATTAGTTATGTAAAGTTGAACAATGTGATTCTTAGAACGTTGACTTGATCAACCTCACAGTATGGTCCAAGCTTTTAACTGACAACTGTTTTCCCAAAACCCAAAAAGATTATGTGTAAAAGAACAATCCACTTGTATTCTGCTTGCCCTCAGAAAACAAGGTAAGGTCTGTGATTCTCCAGCCAAGCCCTCCGTTGTAGACTCGATGGTTGGCACAAAGCAGGGTTTAAATGACCCAAATGGTGTTTTTATGTTACAGTTTAGTTTGAAGTTTTTGGAAAGCACTCTTACCAGTAACAAGCCCAGATTTGCCAGCTACAGCTAAGTGCTTTTCAAGATTCCTAATGTATGGTCAAGAAGAGTGTGGCTGCTCTCTGCTTCACTATGAAGGAAAGCAAGGCTCAGAGAAATGGGTGGCTCCCAGCCTGTAGGGCAGGCAGAGCAGCAGTCAGGAGTGGGCTGTCTGGATTATTGTTTTCCTAGTGTGCCGTGTGAATACGTGTGTGTGTGTGTGTGGAATACGGTGGTGGTTGTGGTGTGTTCATGGGTGGGATAGCTAAAATGTGTGCATGTGAAGAGCAAACTTCAGATGCCATTAATTTAGTCACATAACCTTTTTCTCCTCTTGAAGTAGGGTCCCTTCTCACTGGACTGGCAGTTAGCCACTACAGCTGACCAACTGATAAGGAACCCACACACCCTCTATACATCATCCCCACACACAACTGGGGGCATTACACCTTGCCATCATGCTTAAATGCTGTTCAACCTAAATATTTAATGAATAAATAATAATTAATTCCGAACAAAAAGTAAATAGATAACTTCCCTCCCTTCTCATACTTCCTTCACATTTTTCGGGTACAGATCAGGCACCACACTTCAAGAAAGCACTTTCTAGGATGAGCCAATATAAAAGCATTATCCTATAGCATTAAATCATTTTCTCATTGCAATTGTAAGAACCATCAGAGAACCTCGGTTCATTAGTCTCCAGGGCTCCTTAACCATAGTAGAGTTATATAAACAATGTCGTATCAACAATTCGTTCAAGGGCTGATGTCCTGCAACTTGTCTAACCTTCTCAAGCATGATGCTTACCAGTGTTTTCCAGGCTCATATTGAAACTGCCAAAGTGTGATTATCCAGATAACTATTTAAGGGAAATGTTTGAGAACACGTCTGCATTCTGAGGCAGAGAGAATATACCACCTTAACAAGGGAAGAGGAATGAAAAAAAAACTGTAGCAAGGTGAGGGGGAACGAGATTGAAATATCGAGGACCAATCGTCAGGTTGATGCAGTTTAAATGATGTAGGCCATAGAGTAAGTGGAAAACCCAACCAAAATGGAGAAATGTATAGAACACTGGAGACCAAGTGGAACACACACACACACAACGGGATGGTTGTAGAAGATAGAATCTCTAGATATAGAAATGGCAGCCCACCTCTCCCTGTTCCACCTGTCCCCAAGAAGTTGGACCTGTGACCTTGAAATTTCCACTGCTTAGAGGAACATTTTCAGGGTGGTCCTTATCCGGATTGAGAGAGAGAGTTCCAGATGAAACAATTGGACAAGTCACGGTGCTAAACCATAGCTGATGTATCCCATGCCCACGAATATATTGTAATCAATTTAATTCTGATACACACAAGCAAAAGGGACGGGAGCAGAGATGTGAGCTTTTGCAGTTACAAGTACTTATTAAATTCTAGGCTCGACGTTTAGCACTATTGCTGGGGGTTTTGCTAGGGAGGATGGCCAGTAGGTGGTTTGGGAGGTGCTACCACACACTTTCCCCCCAGCACTCGGGTGGATAAGGTGGCAATGTGTAGAGTGTTCTGCAGGAGCCAGCTCTAGTCATGTATTAGAGAAAACATTTGGGGATTGCTGCTTTGTTACCATATGTTGATTTCTGTGTCTGGATTCTTGTAGGACTCGCCAGATCCATTTACCCTGGTTCCAGTTTTATGCAAAGCCTGGTTGACTTTTCTGAGAACATCGAGGTTACTGTAAAATCTAGGGCAAGGGCCTAGGCTCCTGGCAGAGTTTCAGTATGAGTATGCTAAGTACAAACTAATTGGTGGGTTATGATATGGTCTCTCACTGAGAAGCAGTGATCTACAAATACTATATCATAACATCAGGAAGATGACCAGGGCAGGGAGTTCCAGGGCAGCCTGGGAGGTGGCTACATAAGCCAACTCTTGTCAAAAAAGATTTGAAAATTATAAATCAATTTGAAAATTATTTTGTTACTAGTAGAAGTGACTAAATTCTGGATCTTTAAGCATGTATAATTTTAGTTTAAAATATATGCTTTATGTGTTCTTACCAAATGCCATTTTTTTCCTGTAGGGTTTCTTAAAAGGTGTTGGCATTTTTGTAAGAAGTGCTGTTAAATAATAAATGGAGTTTGACTTGACTTTTTTAAATAGAATAATTTTATTCCTAAAACAAACCTTGGTCCATCTTTTGTAAAAAGAATAAAGTTAATTTCTGAAGGTTGTCCAATAATGTATCTAATTGTTAATTTTTAAAAAAGGGTGGTCTTGACTCAGCCTGCTATTACAGAATAAATTAGATGTTTTGTATATAATAAACAGTTATATCTCATACTTCTGGAGTCCACTACTAAAATGTTCCCTTTCTCCTTTATCAAGTACAGTATCTATTTTTATAACAGCTGATGCCCCCTGCCCCCAGCTGCTTCCCACAAACTTGGATATTTTGGTGGTGCGAGAAGTCATGAAGCCAACATTTCAAAGGGCTCCCATTTAGCACATGTCTTTACACACGCATACAGGTGCTCTCATATGCATAATATTAGAAAGGCCTGGCAGAATAGCAACGAAGACCAATATTTATAGTATTGTTTAGGAATTAAGTGTTTTAGGTATTTTCTGTGTGTTTAGAACCATTTAAATCCTTTCAAAGCTCAAGTAAAGATGAAACGGCATCCTCAGCACTTCATACAAGGAAGAACCGAGGAACGGAGCCTGAAGGGGTTTGCATGGAGACTGTGTAAGCTAGTGGAAGAAAAATAGCATTCCAATGGTGATCCCAACTTTTGTCCCCTTTGAAGGAAGATTCGCCTTTTAAAGCTCTGGTGCTCAGGTCAAGCACCTGTAATTGGCCGTGAACATTTTGTTTATTTTATTTTATTTCATTTTATTTTATTTATTTATAGTATCAAGAATGATATCCTGCCTCCTAGCCCCTGAGTGTCTGTTCAAGGCTATAGGTATAGAAACACTGATATCCACTTTGAAAGGTTTTTAATGCTCCTTCGGTGCGCTGATGACTGTTGGGAATCGGATGCCTAGCTTTTCTCACAGACCTGTAATGGCTCTGAGCTGATACCACACCAGTAAGGAACAACCAACAGAAGGGGGATCTTGGAAATAAAATGGTAAAAGCTGAATGGACATTGAGCTGTCTACAGGTTTCTTGGGCAAAACTGACTCAGTTGCCATTTCCTCATTCTGCAAAAACAGAGTTTAAAGTAAAAGATATGCTGGCGCTTCTAGCTCTGTGAGTCCATGAACTTATCTTCTATAAAATTCTCTTGTAATTCCGGGAAGTAGATGAAAATTAAGGAGTTGTTTCATGTGGAACCAGTTTCCTTTCCAAGACTTCTTATCAACTGGTGTGTGCCTAAGAAGGCAGATGGTTTTTACCTTTGATTTGGGTGCACCAACCTCAATGAGGGAGGAGGAAACGTGTGCGCGCACACAAATACTCCTGAAAGGTAAGGGGACTGATTTGAGGGATTACAGGTGGGTTCATGATTCCCCACGCCACTGCCCCCTCGCAATGGCTGCGGTAATTGAAAAGCCGCCTGCGTTCTATTGATTGACAGCTGAAGGATCCGTTGCCAGGGAACCTCTGGCCGTGTGCCTAAGGTAACGACCTGTGATTAATTATAACTACTTGTGCTGTCATGCACCCCTCGCTCCCTCCCTCTGAGTGGCGTCCCAGTAATCGGATATCAAATGGAACAGTAAAAATGCCCACACTAAATGAACAATACCTGTGCCTGTGGCAGAACACATGTCAAACAGAAAAAAAAAGACTTTTCGCACCCCCTCCGACAGATCTGCTCATGAATACCTGTTGAAAACAATCATAAATTCAAAGTGGGCTGACCAAACACACACGTTTCCACATAGAAGGGTCCAGAGAAGTAGAGGGTCCAGAAGCCGCAACGTGTTCCACATAGAAGGGTCCCAGCAAAACCTGGCCACACCGTGTTTCCATATAGAAGGTCCAGAAACGGCCACCCGTGTTACACCATAGAAGGTCCGAAACCGGCACACACCCGTGTTTCCACATAGAAGGGTCCAGAAACCTGGCCACACCGTGGTTTCACATAGAAGGTCCAGAAAACCTGGCCAACACCGTGGTTTCCACATAGAAGGGTCCAAAACCTGTTCTTCTCCATCCACTCTGACGCTGCTTGAGATTACCCAAATCTGGAGACTCCACGATAAAAACAACATTCTGCCATGTCCCAAAAGCTTCTGTACTTGGGACTGACTGACTTAAACTAGTAATAGTGCTGGATTTTAGTTTGGAGGGTTACCTAGAATTCTGGACCCCTTCTGAAATTCAGGTACCAAGTGACAAGATTTCTGGACCCTTCCCTATGTGGAAACACGGGTGTGGCCAGGTTCTGACCCTTCTATGGGAAACACGGTGTGGCAGGTTTCTGAAACTTCTATGTGGAAACACGGTGTGGGCAAGGTTTCTGGACATCATATGAACACGGTTGTGCCAGGTTTCTGGACCTCTATGTGGAAACACGGTGTGCCAGAGTTCTGGACTCTTCTATGTGGAACACGGGTGTGGCTGGTCAGACCCCTTGAATTTATGATTGTTTTCAAAGTATGTCACTGATCAGAGCGGTCGAGGGGAGTGCGAAAGCTTTTGACATGTGTTGCAAACAGACACAGGGATTGTTCATTAGTGGGCTTACTGTTCCATTCTGAGATCTGAGGTGAAGCCACTCAGAGGAGGGAGCCGATGGGCAGGAACAGCACAAAGGAGGATAATTAATCAGGTCGTTAGCTTGGCACACGGCCAGAGGTTCCCTGGCACAGATAATTCAGCTGTCAGCAATCAGAACAAGCAGGCAGCTTTCAATTACCGCAGCCATTGGGAGGGGGCAGGACGTGGGGAATCTGACACCAGCTGATAATCCCTCCAAATCAGTCCCTTACCTTCGGATGGATTTGTGATGCGCGCACCAAGTTTCCTCCTCCCTCATTGAGTCCACCCAAAGCAAAGTAAAAACCAGATCTGCCTTCCTTAGGCCACAACCAAGTTGATAAAAGGTTGGCCAAAGGAACTGGTTCCCACAGAACACTCCTGAGGGTCACTACTCACGATTACAGAGAATTTTATAGAAGTAAGTTCATGGACTCACAGAGCTAGAAAGCCCGCAGATCTTTACTTTAACTCTGTTTTTTTCGATGGAGAAATGAATGAGTCAGTTGCCCAGAAAAGATTGGTAAGACAGCGCAATCCAATTCAGCTTTTACCATTTTATTGCCAAGAATACACTTGTGTCTGGTTCCTTTAACTAGATGTGTGTCAGCTCAAGAGCCACTACAGGTCTGTGGGTAGAAAAGCTAGGACAGATTCCCAACAGTCATCAGCGCACCGAAGGAGCATTAAAAACCTTTTCAATGTGGAATAGTGTTCTATACTAAGCATGACAGACACTCAGGGGAGAGGCAGGAAGGATCAGACTTGAGACTATAACATAAAAAAAATAAAATAAATAAAATAAATAAACAAAATGTTCCAGGCCAATTACAGGTGCGACCTGAGCACCAGAGCTTAAGGCGAATTTCATTCAAAGGGGACAAAGTTGGGAATCACCAATTGAGAATGCTCATCTTTTCTTCCACTAGCTACGACAGTCTCAGGCAAACATATCAGGCTCCGTTCCTCGTTCCTCTGTAAGAAGTGCTGAGGATGACCGTTTCATAACTTTACTTGAGCTGTTGAAAGGATTAAATGGGTCTACACACAGAATACCTAAAACAGTTCTAATGCCTTAAACAATACTATAAAAGATGGTCTGTTGCTATTCCTGCCAGGCCTGCTAATATTAGTCATATGAGAGCACGCTGGTATGCAGTGTGTAAAGACATTGTGCTAAATAGGGAGCGCTTCATGAAATGTGGCTTCGACTTCTCGCCCAACAAAATATCCAGGTGTGGAAGCAAGCTGGGGGGCAGGGGGACATCGCTGTTAAAATAGATCTCACTGCAACTTTGATAAAAGAAAGGGAACATTTTAGTAGGTTGGACTCCAGAAGTATGAGAGATAAATTCTTTATATACAAAACCATCTAACAATGGATTCTGTAATAGTCAGGCCTGAGTCAAGACACCTTTTAAAATTACAAGTAGATACATTATTGGCAACCTTCAAAATTAACTTTAATCTTTTTAACAAAAGAAGGGGACCAGTTTTGTTTTAGGAAGACAAGTATGCTATTTAAAAAGTCCAAGTCAACTCCAGTGGATTATTGACAAGCACTTTCTTACAAAAATGCCAACACCTTTTAGAAAACCACAGGAAAAAAATTGTGCATTTGGTAGACACATAAAGCATAGATTTTAAAATTATAATTGATTAAAGATTCCAGATTTAGTCCACTTATCTAGTAACAAACCTAATTTCAACATTGATTTATAATTTTAAAATAATTTACTCTTTTTTGACAAGAGTGGCTATGTAGCCCAGGCTGCCCTGGAACTAAATCCAAATAAATTCGAAACTCACCCCATGCCCTGGTCCATCTTCCTGAATTTTAGGAAATAGGAGTGTATAACACGCTCTCTCAGTGAGAGACTAATCATGAACCCACAATGAGTTTGTACTTTAGCAGTCTCATACTGGAACTTGCCAGAGCCTAGGCCCTTGCCCTTAGAGTTTACAAGTACCTCGATGTTCCAGAAAAGTCAACCAGGCTGGGTGCATAAACTGGAACCAGGATAAATGAGCCTGGCGAGTCCTACAGAATCCAACACAGAGAAGAAAGGATGGTAGACAAAGCAGCAATCCCCAAAGTTTTCTCTAGCATATCCCCGAGTGCTGGGGAAGATGTGTGGTAGCACTCTCCCACCCTACTGGCCATCCTCCCTACAAAACCACAAGAGTGCTAAACCTGAGCCTAGAATTAATAAGTACCTTGTAACTGCAAAGCCACAGTCTCTGCATATCCGGTCCCTTTGCTTGTGTGTACAAGAATTAAATTAGATAACCACATATAACGTTCGTGGCATTGGGATAGCCATCACTGATGGTTTAGCACCGTGCTTGTCCATTGTTTGACTGGAACCTCATCCACCCAATCCGGATACAGACCACACCCAGGAAAAGGTACTCTTAAGGCAGTGGAAAGGCAAACAGGTCCCAACAGGCTGGGACAGAGGGAACAGAGGAGAGGGTGGCTGCCTTTTTATATCTAGAAGTCTATCTTCTACAACCATCCCGTGTGTGTGTGTGTGACACTTGGTTCAGTTTCATAGATTTTCTCCCATTTTGGTTGGGTTTTCCCATTATAACTATAGGGCCATAACAGGGAAAACACTGCATCAACCACTGATTAACAGATGGTCCTTATATTTTCATCTCATTTTCCCTCACTTGTCTTAACAGTGTTTGGTTTCATTACTCTTTTCCCTTGTTAAGTGTAATTCTCTGTCCTCCGATGCAGAGTGTTATCAAAACATTTCCCTTAAATAGTTATCGGATTACCACACTTGGTCAGGTTTCCAAAATTATGAGGCTGAAACCACTGGTAAAGCAGCATGGCTTGAGAAGGTTGAGAAGTTGGCAGGACATCACCTTGAACGATTGTTGAAGCGACAGTGTTTATATAACTCACTAGGTTAAGGAGCCCTGGAGACTAATTGAACGAGTTCTCTGATGGGGCTTTTACAATGTCATGAGAAAATATTAATGCTAATAGGAAAGGCTAGGGAGATGATCCAGCCGAAAGTGATTCCTTGCAAGTGTGGGGACCTGATCTGGTACCCCAGAACAATGGTGAAAGGAAGGAGGAAGGAAGGGAGGAGGGAGGGAGGGAAGGAAAGGAGGGAGGAGGGAGGAAATGAGGGAGGGGGAACAGCATGATGCATGATGGCAAGTGTAATGCCAGTTGTGTGTGGGATGATGTATAGAGGGTTGTGTGGGTTCCATAGTCAGTTGGTCCAGCTGTAGTGGCTAACTGGCCAGTCCAGTGAGAAAACCCTACTCAAGAGAGAGAAAAAGGTGGATAGTGACTAAGGAGATGGCAATCTGAAGTTTTTTGCTCTTCGACATGCACACATTGATCTATCCACCCAGACACACACACACACACACACACATGCACCACAACACACACACGTAATGCAACACGTGCACACTAGGAGAAAACAATAAATTCCAGAACACCCACTACTGACTGCCGCTCTGCCTGCCACTTACAGGCTGGTGAGCCACCCATTTCTCTGAGCCTTGCTTTCCTCATTTAGTGAAGCAGAGAGCAGCACACACTATTCTTGACAGACATTAAAGGAATCTCTGAAAAGACTTAGCTGTTAGATGGCAGCGGCTGGTTACTGGTAAGAGTGTGCTTTCCAAATAACTTACAAACTAAAACTGTAACCTAAAAACACCATTTGGTCATTTACAACCCATGCTTTGTGACCAACCAATTCAGAGTCTACAAACGGAGGCTTGGGCGGAGAACACAGCCTTACGGTTTCTGAGGGCAAGTCAGATACAAGTAAGGATTGTTCTTTTTACACATACCTTTGCGTTTTTGGTGAAAACAGTTGTCAGTTAAAGCTTGACACAGTGTGAGGTTGAGGCACAGGTCAACGTTCTAAGTAAGCAACATTGTGCAACGGACCATAACTAATTTTTCCCGGAGTGCTGTGCACATTCGTTACCACACCTTGTTTAACTGGACAGTGATAAGCACCTTAGGACATCAGCTTACTTGTAGGACTTTGAACTTTGTAGACAGTGAACAAGCTTCACTACCTATTCTAATGTTGGCGATTGTGCCAGGAAGACCATTTTTCCCAGTTGCTTTCATTGAAAGATGCCTGCAACAACTATGTTAAAAACAATACATCCATCCAGAGGTAGGCTACTCTGAAAAAATGTAAGTTCCTGATTGTTTTTTTTGGTATTTGAGACAGGTAAGAAATCACACCATATCCTTATAAACAGAAAACCTTTCCCACTGAGGGGAGAACCCATCTTATTATTCGTATAATGGTTCCTGTTGTGCACCAGGGGCCTTTTTCAGGCTTTCTGTTTGCCACACTTCAAGTAAGAAGCTCAAACACACATTAACTTGGCCCACTGTGGTGATGCCTTATGTGTGAGCATGGAACTCCTCTCCAGACTCAGTCCTACCTCGCCACACTATTGTTGCCCAAGTCCTCCTCGGATGACTTTGAGGCAAAGTCTAGTATGATATCATTTTTATGTGTAACCCACCCCATTCAAGCAAACAATTGATAAATCATGCCACCAAATTGCAAAGGTTCAGATTTCTAAATGTCTCATGATATCATAATTTTTTAAAAAATCAGAATTCGGGCTGGAGAGATAGGATTGGCTCAGCAGTTGAGAGCACTTTCTGCTCTTTGCAGAGCACCAGGTTAGGGTCCCACAACCATCTGTTTACCTCAGTTCTGGGGGAGAGAGACCCTCTTCATGACTTTGCCAGCACAGGCATCCATGGTTCATATACATACACACGGCAAAAACACCAATGCACATAAAAAATAAAATAAATCAAAGACATACAAAAACAAAACAAACAAACAACACCAACAAACCAGGATCACCAGGATAAAGTGTACTCTTTGCAAATGATGTATATGTCTCCTGATGTATACTTTGTTCTCCAATCCCCTTAGCCTTGGTCCACTTTGCAATTACTTTGTTCAAGAATGTAGTTGCTTGTCTAAAAACATTTTCCCTGATTTTGCTGACATATCCTGATGGATACAGTTAACATCTTCTGTTAGATCTTCTAGGAACTGTAGTTGGATGCTTGAGCCACAAGGTAAGATTTTTAAGTCATTTCCTGGCTTTTTATTATAATTTTATTCATTCAACTTTTGCTTAATAGTTTACTACATGTTTTATCTCTGCAGTTTGCTATTGGAAACGAAAAGTGAAAGCAGAAACTTGGAAGACAGTCAAGGGATACAAACCAGAATACGTTGAGTATTTGGTGTAGTAATTTCTCTTTAATTCTTCTTCCTCTTCCTCTTCTTCCCTTTTTAGTATTCTTAATCTTATTATGTCCTTCCCTCCATATCCTTCTTGATGTTAAAGTTTAAATCACTGCTAAGGCAGTCTCATGGGAGTAGCTTCTAAAAACTTATAAGCCTATAAAATAGCAATTTTCTTTGATATAAACCTCAGCTCATCTTTTCCACAGGTGTTTCTCTCACCCATAATTTCCCTGACTATAATGAGTCTTATAATAGGGCCATAACCAAATGATATGCAAAGAATAGCAAAATCTCGCTTCTCTCCTCTCTCTATTATATCTCTCTCTCGCTTTCTATCTCTCTTTCTCTCTCTCTATCGTCTCTCTCTCTCTCCTCTCTCTCACACACACACACACACACACACACACACACACAAACAAAAAGTTCTGAGACAGAAAGTTGTTTATCTTTTACTTTAAAATGATAAGCAAGTTTAACGTATGCATTCAGGCCTGCAATAGGAACATGACACAATGAGAAGACGTGCAAGTATTTTGCGGGGATCCGTCAGTCATAAAAACGAGAGTCTCTTACTTTTGCATGGCTTTCCATCTGTTATGGACTCAGCTCTGCATCAGGAAACGGTTATTGTATATGGAGACAGTAGCAGGAGAATGTGAAACTGTCCACCAGACCACACAGCTTATTACAATTGAGTGAGGACTCTGACCCCAAGCCGGGAATTCTTTATATAATTT";
        let mut mapping = aligner
            .map(
                seq,
                false,
                false,
                None,
                Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
                Some(b"t"),
            )
            .unwrap();
        mapping.sort_by_key(|hit| hit.query_start);
        mapping.iter().for_each(|hit| {
            println!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, identity:{}, score:{:?}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                MappingExt(hit).identity(),
                hit.alignment.as_ref().unwrap().alignment_score
            )
        });
    }
}

/*

sr
aligner.mapopt:mm_mapopt_t { flag: 1077997644, seed: 11, sdust_thres: 0, max_qlen: 0, bw: 100, bw_long: 100, max_gap: 100, max_gap_ref: -1, max_frag_len: 800, max_chain_skip: 25, max_chain_iter: 5000, min_cnt: 2, min_chain_score: 1, chain_gap_scale: 0.8, chain_skip_scale: 0.0, rmq_size_cap: 100000, rmq_inner_dist: 1000, rmq_rescue_size: 1000, rmq_rescue_ratio: 0.1, mask_level: 0.5, mask_len: 2147483647, pri_ratio: 0.5, best_n: 20, alt_drop: 0.15, a: 2, b: 8, q: 12, e: 2, q2: 24, e2: 1, transition: 0, sc_ambi: 1, noncan: 0, junc_bonus: 0, zdrop: 100, zdrop_inv: 100, end_bonus: 10, min_dp_max: 1, min_ksw_len: 0, anchor_ext_len: 20, anchor_ext_shift: 6, max_clip_ratio: 1.0, rank_min_len: 500, rank_frac: 0.9, pe_ori: 1, pe_bonus: 33, mid_occ_frac: 0.0, q_occ_frac: 0.0, min_mid_occ: 0, max_mid_occ: 1000000, mid_occ: 1000, max_occ: 5000, max_max_occ: 4095, occ_dist: 0, mini_batch_size: 50000000, max_sw_mat: 100000000, cap_kalloc: 500000000, split_prefix: 0x0 }

aligner.mapopt:mm_mapopt_t { flag: 1077997644, seed: 11, sdust_thres: 0, max_qlen: 0, bw: 100, bw_long: 100, max_gap: 100, max_gap_ref: -1, max_frag_len: 800, max_chain_skip: 25, max_chain_iter: 5000, min_cnt: 2, min_chain_score: 1, chain_gap_scale: 0.8, chain_skip_scale: 0.0, rmq_size_cap: 100000, rmq_inner_dist: 1000, rmq_rescue_size: 1000, rmq_rescue_ratio: 0.1, mask_level: 0.5, mask_len: 2147483647, pri_ratio: 0.5, best_n: 20, alt_drop: 0.15, a: 2, b: 8, q: 12, e: 2, q2: 24, e2: 1, transition: 0, sc_ambi: 1, noncan: 0, junc_bonus: 0, zdrop: 100, zdrop_inv: 100, end_bonus: 10, min_dp_max: 1, min_ksw_len: 0, anchor_ext_len: 20, anchor_ext_shift: 6, max_clip_ratio: 1.0, rank_min_len: 500, rank_frac: 0.9, pe_ori: 1, pe_bonus: 33, mid_occ_frac: 0.0, q_occ_frac: 0.0, min_mid_occ: 0, max_mid_occ: 1000000, mid_occ: 1000, max_occ: 5000, max_max_occ: 4095, occ_dist: 0, mini_batch_size: 50000000, max_sw_mat: 100000000, cap_kalloc: 500000000, split_prefix: 0x0 }
ont
aligner.mapopt:mm_mapopt_t { flag: 1073741900, seed: 11, sdust_thres: 0, max_qlen: 0, bw: 500, bw_long: 20000, max_gap: 5000, max_gap_ref: -1, max_frag_len: 0, max_chain_skip: 25, max_chain_iter: 5000, min_cnt: 2, min_chain_score: 1, chain_gap_scale: 0.8, chain_skip_scale: 0.0, rmq_size_cap: 100000, rmq_inner_dist: 1000, rmq_rescue_size: 1000, rmq_rescue_ratio: 0.1, mask_level: 0.5, mask_len: 2147483647, pri_ratio: 0.8, best_n: 1, alt_drop: 0.15, a: 2, b: 4, q: 4, e: 2, q2: 24, e2: 1, transition: 0, sc_ambi: 1, noncan: 0, junc_bonus: 0, zdrop: 400, zdrop_inv: 200, end_bonus: -1, min_dp_max: 1, min_ksw_len: 0, anchor_ext_len: 20, anchor_ext_shift: 6, max_clip_ratio: 1.0, rank_min_len: 500, rank_frac: 0.9, pe_ori: 0, pe_bonus: 33, mid_occ_frac: 0.0, q_occ_frac: 0.0, min_mid_occ: 0, max_mid_occ: 1000000, mid_occ: 1000, max_occ: 0, max_max_occ: 4095, occ_dist: 0, mini_batch_size: 500000000, max_sw_mat: 100000000, cap_kalloc: 500000000, split_prefix: 0x0 }


*/
