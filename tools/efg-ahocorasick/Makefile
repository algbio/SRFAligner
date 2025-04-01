CPPFLAGS=-Ofast -march=native --std=c++20

all : efg-ahocorasick extractor

efg-ahocorasick : src/efg-ahocorasick.rs ../daachorse/target/release/libdaachorse.rlib
	rustc src/efg-ahocorasick.rs -C opt-level=3 \
		--extern daachorse=../daachorse/target/release/libdaachorse.rlib

../daachorse/target/release/libdaachorse.rlib : ../daachorse/.git
	cargo build --release --manifest-path=../daachorse/Cargo.toml

extractor : src/extractor.cpp src/efg.hpp src/efg-locate.hpp
	g++ $(CPPFLAGS) $(HEADERS) -o extractor \
	src/extractor.cpp \
