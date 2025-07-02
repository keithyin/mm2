build:
	cargo build --release

bai:
	cargo build --release
	cp target/release/gsmm2 /usr/bin/
	cp target/release/gsmm2-aligned-metric /usr/bin/
	cp target/release/gsmm2-time-err /usr/bin/
	cp target/release/bam-rw /usr/bin/
	cp target/release/gsmm2-hp-tr-metric /usr/bin/

install:
	cp target/release/gsmm2 /usr/bin/
	cp target/release/gsmm2-aligned-metric /usr/bin/
	cp target/release/gsmm2-time-err /usr/bin/
	cp target/release/bam-rw /usr/bin/
	cp target/release/gsmm2-hp-tr-metric /usr/bin/

clean:
	rm -rf target