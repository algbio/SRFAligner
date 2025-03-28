CPPFLAGS=-Ofast -march=native --std=c++20
#CPPFLAGS=-g -O0 --std=c++20

all : efg-gaf-splitter

efg-gaf-splitter : efg-gaf-splitter.cpp efg.hpp command-line-parsing/cmdline.h command-line-parsing/cmdline.c
	gcc $(CPPFLAGS) \
	efg-gaf-splitter.cpp command-line-parsing/cmdline.c \
	-o efg-gaf-splitter

# uncomment for development
#command-line-parsing/cmdline%c command-line-parsing/cmdline%h : command-line-parsing/config.ggo
#	gengetopt \
#		--input=./command-line-parsing/config.ggo \
#		--output-dir=./command-line-parsing/ \
#		--unnamed-opts

.PHONY : clean all cleanall

# uncomment for development
#clean :
#	rm -Rf command-line-parsing/cmdline.c command-line-parsing/cmdline.h
#cleanall :
#	rm -Rf efg-gaf-splitter command-line-parsing/cmdline.c command-line-parsing/cmdline.h
cleanall :
	rm -Rf efg-gaf-splitter
