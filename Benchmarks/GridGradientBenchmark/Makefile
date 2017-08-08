all: ChiBenchmark_CPU ChiBenmark_GPU


#
#
PROGRAM_NAME_GPU := ChiBenchmark_GPU
PROGRAM_NAME := ChiBenchmark_CPU
#all: $(PROGRAM_NAME_GPU) $(PROGRAM_NAME)
all: $(PROGRAM_NAME) 
#CXX=g++ -lm -ffast-math -ftree-loop-vectorize 
CXX=xlc++
#
program_CXX_SRCS := $(wildcard *.cpp)
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
#
program_C_SRCS := $(wildcard *.c)
program_C_OBJS := ${program_C_SRCS:.c=.o}
#
program_CU_SRCS := $(wildcard *.cu)
program_CU_OBJS := ${program_CU_SRCS:.cu=.cu.o}
#
# include
#
program_INCLUDE_DIRS := . /usr/local/cuda/include/ #C++ Include directories
program_INCLUDE_DIRS += $(CFITSIO_ROOT)/include
program_INCLUDE_DIRS += $(LENSTOOL_ROOT)/include
program_INCLUDE_DIRS += $(GSL_ROOT)
program_INCLUDE_DIRS += $(LENSTOOLHPC_ROOT)/src
#
# libs
#
program_INCLUDE_LIBS += $(CFITSIO_ROOT)/lib #Include libraries
program_INCLUDE_LIBS += $(LENSTOOL_ROOT)/src/
program_INCLUDE_LIBS += $(LENSTOOL_ROOT)/liblt
program_INCLUDE_LIBS += $(LENSTOOLHPC_ROOT)/src/
program_INCLUDE_LIBS += $(GSL_ROOT)/lib
program_INCLUDE_LIBS += $(WCSTOOL_ROOT)
#
#
# Compiler flags
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CPPFLAGS += $(foreach includelib,$(program_INCLUDE_LIBS),-L$(includelib))
CPPFLAGS += -D__WITH_LENSTOOL 
#CPPFLAGS += -D__USE_GPU
CPPFLAGS += -qsmp -g -O3 -Wall -pedantic
CPPFLAGS += -pedantic  
CPPFLAGS += -fPIC
#CPPFLAGS += -no-vec
#CPPFLAGS += -qopenmp -march=core-avx2 -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -llenstoolhpc -qopenmp -xHost -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -llenstoolhpc -qopenmp -axMIC-AVX512,CORE-AVX2 -g -O3 -std=c++0x -Wall -pedantic
#CPPFLAGS += -qopt-prefetch-distance=64,8 -qopt-streaming-cache-evict=0 -llenstoolhpc -qopenmp -xMIC-AVX512 -g -O3 -std=c++0x -Wall -pedantic
LDFLAGS := -llenstoolhpc -lcfitsio -llenstool -llt -lgsl -lgslcblas -lm -lwcs

NVFLAGS := -O3 -g -G -arch=sm_35 -rdc=true -ccbin icpc -Xcompiler -qopenmp -D__WITH_LENSTOOL #rdc=true needed for separable compilation
#NVFLAGS += --maxrregcount=40 
NVFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
NVFLAGS += $(foreach includelib,$(program_INCLUDE_LIBS),-L$(includelib)) 

CUO_O_OBJECTS := ${program_CU_OBJS:.cuo=.cuo.o}
C_O_OBJECTS := ${program_C_OBJS:.c=.o}

GPU_OBJECTS = $(program_CXX_OBJS) $(program_C_OBJS) $(program_CU_OBJS)
CPU_OBJECTS = $(program_CXX_OBJS) $(program_C_OBJS)
#
%.cu.o: %.cu 
	nvcc $(NVFLAGS) -o $@ -c -rdc=true -std=c++11 $<
	#nvcc $(NVFLAGS) -Xptxas -dlcm=cg --ptxas-options=-v -o $@ -dc $<
#
%.cpp: %.cpp %.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(CPPFLAGS)
#
$(PROGRAM_NAME): $(CPU_OBJECTS) 
	xlc++ -o $@ $^ $(LDFLAGS) 

clean:
	@- $(RM) $(PROGRAM_NAME) $(OBJECTS) *~ *.o *.optrpt

distclean: clean

.PHONY: all clean distclean
