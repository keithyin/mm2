

requirements
* samtools : https://www.htslib.org/download/

```
# install samtools

cd samtools-1.x    # and similarly for bcftools and htslib
./configure --prefix=/usr/
make
make install
```



```
cargo install mm2

gsmm2 align -q query.fa --target target.fa -p query2target
```


# changelog

## 1.2.0

gsmm2 align

* new flag `--hpTrShrinkAlnRegion`: when an alignment end falls inside a homopolymer or a
  tandem repeat of the target (unit <= 4, repeats >= 3), the alignment region is shrunk so
  the end sits outside the repeat. the trimmed read bases become soft clips.
  A boundary placed inside a repeat is arbitrary -- every phase of the run is equivalent --
  so the bases there are not trustworthy and should not be reported as aligned.
  Migrated from gsmm2-metric's per-metric-target window shrink
  (`hp_tr_metric_v2.rs` / `hp_tr_metric_v3.rs`), without its second alignment pass and
  without the `get_target_substr` same-base trim.
  * alignment score, mapq and nm are left as they were. cs and md are rebuilt from the
    trimmed columns, byte-for-byte in the format minimap2's `write_cs_ds_core` /
    `write_MD_core` produce, so the record stays self-consistent
  * the repeat scan iterates to a fixpoint rather than taking one step like gsmm2-metric:
    motifs nest (`ACAACAAC[A]AAA` is covered by both `(ACA)3` and `(A)4`) and abut, so a
    single step can stop on a position that is still inside a repeat
  * if the whole alignment sits inside a repeat the shrink would collapse it to zero
    reference span, which breaks sorting/indexing, so the hit is left alone and a warn is
    logged
  * one interval tree per target, attached to `NoMemLeakAligner` next to `target_seq`,
    built on the sequence that aligner actually indexes -- so the `___rev` targets of
    `--query-forward` need no coordinate conversion. The scan is ~340 regexes over the
    target once, ~7 MB/s per target thread (0.7s for the 4.6 MB E. coli reference), and
    only runs when the flag is on
  * runs before `--polyNGapLeftAlign`, so a gap is never left-aligned onto a column that is
    then trimmed

## 0.23.0 

* more tags to pass through

## 0.22.0

* Pin minimap2 to version 0.1.23

## 0.21.0

gsmm2-aligned-metric

* new field: longIndel. only the abs(indel) >=10 will be recorded. 14,-31 means there are two long indels (ins(14), del(31))

## 0.20.2

gsmm2-aligned-metric

```rust
// 1 -> 1000, 0.8 -> 0.2
aligner.mapopt.best_n = 10000; // Output at most INT secondary alignments
aligner.mapopt.pri_ratio = 0.2; // Minimal secondary-to-primary score ratio to output secondary mappings
```


## 0.20.0

* homodel logic

## 0.18.1

* primaryCovlen
* qOvlpRatio
* rOvlpRatio

## 0.18.0

* gsmm2-aligned-metric.rs
* align_single_query_to_targets TO (align_single_query_to_targets + hits2records)

## 0.17.0

* set_primary_alignment_according_2_align_score -> set_primary_supp_alignment_according_2_align_score

## 0.11.0

* np_range & rq_range. 
* discard multiple alignment reads


## 0.9.0

* dump qual if the query file is bam or fastq. if query file is fa, the qual will be [255; query_len]

## 0.6.0

* if the query file is bam format, it will try to dump np & ch & rq to the result bam

## 0.5.0

* if the query file is bam format, it will try to dump np & ch to the result bam