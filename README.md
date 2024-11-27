

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

## 0.11.0

* np_range & rq_range. 
* discard multiple alignment reads


## 0.9.0

* dump qual if the query file is bam or fastq. if query file is fa, the qual will be [255; query_len]

## 0.6.0

* if the query file is bam format, it will try to dump np & ch & rq to the result bam

## 0.5.0

* if the query file is bam format, it will try to dump np & ch to the result bam