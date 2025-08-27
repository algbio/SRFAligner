# adapted from https://github.com/jltsiren/gcsa2/blob/master/benchmark/Makefile
SDSL_DIR=../vgteam-sdsl-lite
include $(SDSL_DIR)/Make.helper

GCSA_DIR=../gcsa2
LIBRARY=$(GCSA_DIR)/lib/libgcsa2.a

# Multithreading with OpenMP and libstdc++ Parallel Mode.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) -lsdsl -ldivsufsort -ldivsufsort64
PROGRAMS=gcsa2-locate

# Apple Clang does not support OpenMP directly, so we need special handling.
ifeq ($(shell uname -s), Darwin)
    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        $(info The compiler is Apple Clang that needs libomp for OpenMP support.)

        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        # Find libomp installed by Homebrew or MacPorts.
        ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/include/omp.h ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp installed by Homebrew and linked to $(HOMEBREW_PREFIX).)
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/include
            LIBS += -L$(HOMEBREW_PREFIX)/lib
        else ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
            $(info Found a keg-only libomp installed by Homebrew at $(HOMEBREW_PREFIX)/opt/libomp.)
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
            LIBS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp installed by MacPorts at /opt/local.)
            PARALLEL_FLAGS += -I/opt/local/include/libomp
            LIBS += -L/opt/local/lib/libomp
        else
            $(error Could not find libomp. Please install it using Homebrew or MacPorts.)
        endif

        # We also need to link it.
        LIBS += -lomp
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(GCSA_DIR)/include -I$(INC_DIR) -I ../concurrentqueue
# uncomment for debugging
#MY_CXX_DEBUG_FLAGS= -std=c++17 -Wall -Wextra $(CODE_COVER)
#MY_CXX_OPT_DEBUG_FLAGS= -O0 -g -march=native
#CXX_FLAGS=$(MY_CXX_DEBUG_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_DEBUG_FLAGS) -I$(GCSA_DIR)/include -I$(INC_DIR) -I ../concurrentqueue

HEADERS=$(wildcard $(GCSA_DIR)/include/gcsa/*.h)

.PHONY: all clean test
all:gcsa2-locate

gcsa2-locate: gcsa2-locate.cpp gcsa2-locate.hpp command-line-parsing/cmdline.h command-line-parsing/cmdline.c $(LIBRARY) graph.hpp algo.hpp gafanchor.hpp
	$(MY_CXX) command-line-parsing/cmdline.c $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

# uncomment for development
#command-line-parsing/cmdline%c command-line-parsing/cmdline%h : command-line-parsing/config.ggo
#	gengetopt \
#		--input=./command-line-parsing/config.ggo \
#		--output-dir=./command-line-parsing/ \
#		--unnamed-opts

clean:
	rm -f $(PROGRAMS) *.o
