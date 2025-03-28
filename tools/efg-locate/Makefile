CPPFLAGS=-Ofast -march=native --std=c++20 -pthread
#CPPFLAGS=-g -O0 --std=c++20 -pthread
HEADERS=-I ../sdsl-lite-v3/include -I ../concurrentqueue

all : efg-locate

efg-locate : efg-locate.cpp efg.hpp algo.cpp command-line-parsing/cmdline.h command-line-parsing/cmdline.c
	gcc $(CPPFLAGS) $(HEADERS) \
	efg-locate.cpp command-line-parsing/cmdline.c \
	-o efg-locate

# uncomment for development
#command-line-parsing/cmdline%c command-line-parsing/cmdline%h : command-line-parsing/config.ggo
#	gengetopt \
#		--input=./command-line-parsing/config.ggo \
#		--output-dir=./command-line-parsing/ \
#		--unnamed-opts

.PHONY : clean all cleanall

#clean :
#	rm -Rf command-line-parsing/cmdline.c command-line-parsing/cmdline.h
#cleanall :
#	rm -Rf efg-locate command-line-parsing/cmdline.c command-line-parsing/cmdline.h
cleanall :
	rm -Rf efg-locate
