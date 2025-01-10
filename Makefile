build:
	cargo build --release

install:
	cp target/release/gsmm2 /usr/bin/
	cp target/release/gsmm2-aligned-metric /usr/bin/

clean:
	rm -rf target