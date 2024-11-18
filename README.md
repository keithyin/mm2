

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