CPPFLAGS=-Ofast -march=native --std=c++20
#CPPFLAGS=-O0 -g --std=c++20

all : efg-simplify

efg-simplify : efg-simplify.cpp command-line-parsing/cmdline.h command-line-parsing/cmdline.c
	g++ $(CPPFLAGS) \
	efg-simplify.cpp command-line-parsing/cmdline.c \
	-o efg-simplify

#command-line-parsing/cmdline%c command-line-parsing/cmdline%h : command-line-parsing/config.ggo
#	gengetopt \
#		--input=./command-line-parsing/config.ggo \
#		--output-dir=./command-line-parsing/ \
#		--unnamed-opts

.PHONY : clean all

clean :
	rm -Rf efg-simplify
#	rm -Rf efg-simplify command-line-parsing/cmdline.{c,h}
