build:
	cargo build --release

install:
	cp target/release/gsmm2 /usr/bin/

clean:
	rm -rf target