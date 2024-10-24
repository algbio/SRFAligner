CPPFLAGS=-Ofast -march=native --std=c++20 -pthread
#CPPFLAGS=-g -O0 --std=c++20 -pthread
HEADERS=-I ../sdsl-lite-v3/include -I ../concurrentqueue

all : chainx-block-graph 

chainx-block-graph : chainx-block-graph.cpp chainx-block-graph.hpp chaining.hpp efg.hpp command-line-parsing/cmdline.h command-line-parsing/cmdline.c
	g++ $(CPPFLAGS) $(HEADERS) \
	chainx-block-graph.cpp command-line-parsing/cmdline.c \
	-o chainx-block-graph

command-line-parsing/cmdline%c command-line-parsing/cmdline%h : command-line-parsing/config.ggo
	gengetopt \
		--input=./command-line-parsing/config.ggo \
		--output-dir=./command-line-parsing/ \
		--unnamed-opts

test : chainx-block-graph test/test.sh
	test/test.sh

.PHONY : clean all test

clean :
	rm -Rf chainx-block-graph command-line-parsing/cmdline.{c,h} test/output*
