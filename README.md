

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