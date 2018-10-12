PROGRAM_NAME := ChiBenchmark
#CXX=g++ -lm -ffast-math -ftree-loop-vectorize 
CXX ?= mpiicpc


program_CXX_SRCS := $(wildcard *.cpp)
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}

#program_CXX_SRCS += $(wildcard ../../*.c) #Find C source files from additonal directories
program_C_OBJS := ${program_C_SRCS:.c=.o}
#
program_CU_SRCS := $(wildcard *.cu)
#program_CU_SRCS += $(wildcard ../../*.cu) #Find CUDA source files from additional directories
#program_CU_HEADERS := $(wildcard *.cuh) #Optional: Include .cuh files as dependencies
#program_CU_HEADERS += $(wildcard ../../*.cuh) #Find .cuh files from additional directories
program_CU_OBJS := ${program_CU_SRCS:.cu=.cuo}
#
program_INCLUDE_DIRS := . /usr/local/cuda/include/ #C++ Include directories
program_INCLUDE_DIRS += /usr/include/cfitsio/ 
program_INCLUDE_DIRS += /users/fgilles/Projects/Libs/cfitsio/include
program_INCLUDE_DIRS += /users/fgilles/bin/iaca-lin32/include
program_INCLUDE_DIRS += /usr/include
program_INCLUDE_DIRS += $(LENSTOOLHPC_ROOT)/src 
program_INCLUDE_DIRS += $(LENSTOOL_ROOT)/include 
program_INCLUDE_DIRS += $(GSL_ROOT)/include/ 
program_INCLUDE_DIRS += $(WCSTOOL_ROOT)
program_INCLUDE_DIRS += $(CFITSIO_ROOT)

#program_CU_INCLUDE_DIRS := /home/users/amclaugh/CUB/cub-1.3.2/ #CUDA Include directories
#
program_INCLUDE_LIBS := /usr/lib64/ #Include libraries
program_INCLUDE_LIBS += /users/fgilles/Projects/Libs/cfitsio/lib #Include libraries
program_INCLUDE_LIBS += $(LENSTOOLHPC_ROOT)/src 
program_INCLUDE_LIBS += $(LENSTOOL_ROOT)/include 
program_INCLUDE_LIBS += $(LENSTOOL_ROOT)/src 
program_INCLUDE_LIBS += $(LENSTOOL_ROOT)/liblt 
program_INCLUDE_LIBS += $(WCSTOOL_ROOT) 
program_INCLUDE_LIBS += $(CFITSIO_ROOT)
program_INCLUDE_LIBS += $(GSL_ROOT)/lib/ 
#
# Compiler flags
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CPPFLAGS += $(foreach includelib,$(program_INCLUDE_LIBS),-L$(includelib))
#CPPFLAGS += -D__WITH_LENSTOOL
CPPFLAGS += -D__WITH_MPI
CPPFLAGS += -D_double
CPPFLAGS += -fPIC -qopenmp -xHost -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -qopenmp -march=core-avx2 -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -llenstoolhpc -qopenmp -xHost -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -llenstoolhpc -qopenmp -axMIC-AVX512,CORE-AVX2 -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -qopt-prefetch-distance=64,8 -qopt-streaming-cache-evict=0 -llenstoolhpc -qopenmp -xMIC-AVX512 -g -O3 -std=c++0x -Wall -pedantic
LDFLAGS := -llenstoolhpc -llenstool -llt -lgsl -lgslcblas -lm -lwcs -lcfitsio


OBJECTS = $(program_CXX_OBJS) $(program_C_OBJS)

.PHONY: all clean distclean
#
all: $(PROGRAM_NAME) 
#
debug: CXXFLAGS = -g -O0 -std=c++0x -Wall -pedantic -DDEBUG $(EXTRA_FLAGS)
debug: $(PROGRAM_NAME)

$(PROGRAM_NAME): $(OBJECTS) 
	CC $(CPPFLAGS) -o $@ $(program_CXX_OBJS) $(program_C_OBJS) $(CUO_O_OBJECTS) $(LDFLAGS)

clean:
	@- $(RM) $(PROGRAM_NAME) $(OBJECTS) *~ *.o *.optrpt chi

distclean: clean


