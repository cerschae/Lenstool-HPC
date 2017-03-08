PROGRAM_NAME := ChiBenchmark
#CXX=g++ -lm -ffast-math -ftree-loop-vectorize 
CXX=icpc


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
program_INCLUDE_DIRS += $(CFITSIO_ROOT)/include
#program_INCLUDE_DIRS += /users/fgilles/bin/iaca-lin32/include
program_INCLUDE_DIRS += $(LENSTOOL_ROOT)/include
program_INCLUDE_DIRS += $(LENSTOOL_ROOT)/src
program_INCLUDE_DIRS += $(LENSTOOL_ROOT)/liblit

program_INCLUDE_DIRS += $(GSL_ROOT)
program_INCLUDE_DIRS += $(LENSTOOLHPC_ROOT)/src
#
#program_CU_INCLUDE_DIRS := /home/users/amclaugh/CUB/cub-1.3.2/ #CUDA Include directories
#
#program_INCLUDE_LIBS := ./ #Include libraries
program_INCLUDE_LIBS += $(CFITSIO_ROOT)/lib #Include libraries
program_INCLUDE_LIBS += $(LENSTOOL_ROOT)/src
program_INCLUDE_LIBS += $(LENSTOOL_ROOT)/liblt
program_INCLUDE_LIBS += $(LENSTOOL_ROOT)/include
program_INCLUDE_LIBS += $(LENSTOOLHPC_ROOT)/src/
program_INCLUDE_LIBS += $(GSL_ROOT)/lib
program_INCLUDE_LIBS += $(WCSTOOL_ROOT)
#
#
# Compiler flags
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CPPFLAGS += $(foreach includelib,$(program_INCLUDE_LIBS),-L$(includelib))
CPPFLAGS += -D__WITH_LENSTOOL
CPPFLAGS += -qopenmp -xHost -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -qopenmp -march=core-avx2 -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -llenstoolhpc -qopenmp -xHost -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -llenstoolhpc -qopenmp -axMIC-AVX512,CORE-AVX2 -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -qopt-prefetch-distance=64,8 -qopt-streaming-cache-evict=0 -llenstoolhpc -qopenmp -xMIC-AVX512 -g -O3 -std=c++0x -Wall -pedantic
LDFLAGS := -llenstoolhpc -lcfitsio -llenstool -llt -lgsl -lgslcblas -lm -lwcs

NVFLAGS := -O3 -g -G -arch=sm_35 -rdc=true -ccbin icpc -Xcompiler -qopenmp #rdc=true needed for separable compilation
NVFLAGS += -D__WITH_LENSTOOL 
NVFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
NVFLAGS += $(foreach includelib,$(program_INCLUDE_LIBS),-L$(includelib)) 

CUO_O_OBJECTS := ${program_CU_OBJS:.cuo=.cuo.o}


OBJECTS = $(program_CXX_OBJS) $(program_C_OBJS) $(program_CU_OBJS)

.PHONY: all clean distclean

all: $(PROGRAM_NAME) 
#
debug: CXXFLAGS = -g -O0 -std=c++0x -Wall -pedantic -DDEBUG $(EXTRA_FLAGS)
debug: $(PROGRAM_NAME)
%.cuo: %.cu %.cuh
	nvcc $(NVFLAGS) -o $@ -dc $<
#nvcc $(NVFLAGS)  -Xptxas -v,-abi=no -ptx -src-in-ptx -dc $<

$(PROGRAM_NAME): $(OBJECTS) 
	@ for cu_obj in $(program_CU_OBJS); \
	do        \
		mv $$cu_obj $$cu_obj.o; \
	done #append a .o suffix for nvcc
	nvcc $(NVFLAGS) -o $@ $(program_CXX_OBJS) $(program_C_OBJS) $(CUO_O_OBJECTS) $(LDFLAGS)
	@ for cu_obj in $(CUO_O_OBJECTS);  \
	do        \
		mv $$cu_obj $${cu_obj%.*} ;  \
	done #remove the .o for make

clean:
	@- $(RM) $(PROGRAM_NAME) $(OBJECTS) *~ *.o *.optrpt

distclean: clean