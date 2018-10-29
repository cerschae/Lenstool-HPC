PROGRAM_NAME := GradientBenchmark
#CXX=g++ -lm -ffast-math -ftree-loop-vectorize 
CXX=icpc

program_CXX_SRCS := $(wildcard *.cpp)
#program_CXX_SRCS += $(wildcard ../../*.cpp) #Find C++ source files from additonal directories
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}

program_C_SRCS := $(wildcard *.c)
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
program_INCLUDE_DIRS += $(GSL_ROOT)/include
program_INCLUDE_DIRS += $(LENSTOOLHPC_ROOT)/src
# 
program_CU_INCLUDE_DIRS := /home/users/amclaugh/CUB/cub-1.3.2/ #CUDA Include directories
#
program_INCLUDE_LIBS := ./ #Include libraries
program_INCLUDE_LIBS += $(CFITSIO_ROOT)/lib #Include libraries
program_INCLUDE_LIBS += $(LENSTOOL_ROOT)/src
program_INCLUDE_LIBS += $(LENSTOOL_ROOT)/liblt 
program_INCLUDE_LIBS += $(LENSTOOLHPC_ROOT)/src
program_INCLUDE_LIBS += $(GSL_ROOT)/lib 
program_INCLUDE_LIBS += $(WCSTOOL_ROOT) 
#
# Compiler flags
# 
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CPPFLAGS += $(foreach includelib,$(program_INCLUDE_LIBS),-L$(includelib))
#CPPFLAGS += -qopt-prefetch-distance=64,8 -qopt-streaming-cache-evict=0 -qopenmp -xMIC-AVX512 -pedantic -llenstoolhpc
#CPPFLAGS += -qopt-streaming-stores=always -qopenmp -xMIC-AVX512 -pedantic -llenstoolhpc
#CPPFLAGS += -ftz -qopenmp -xHost -llenstool -fma -pedantic -llenstoolhpc  -llt -lcfitsio -lwcs -lgsl -lgslcblas 
CPPFLAGS += -llenstool
CPPFLAGS += -ftz -qopenmp -pedantic -llenstoolhpc -llt -lcfitsio -lwcs -lgsl -lgslcblas 
CPPFLAGS += -xHost -fma
CPPFLAGS += -g -O3 
CPPFLAGS += -std=c++0x -Wall 
CPPFLAGS += -fPIC 
#CPPFLAGS += -D__WITH_LENSTOOL
include $(LENSTOOLHPC_ROOT)/src/precision.inc
CPPFLAGS += $(PRECISION)
#CPPFLAGS += -D__USE_GPU
#
GEN_SM35 := -gencode=arch=compute_35,code=\"sm_35,compute_35\" #Target CC 3.5, for example
NVFLAGS := -O3 -g -G -rdc=true #rdc=true needed for separable compilation
#NVFLAGS += $(GEN_SM35)
NVFLAGS += $(foreach includedir,$(program_CU_INCLUDE_DIRS),-I$(includedir))
#
CUO_O_OBJECTS := ${program_CU_OBJS:.cuo=.cuo.o}
#
OBJECTS = $(program_CU_OBJS) $(program_CXX_OBJS) $(program_C_OBJS)
#
.PHONY: all clean distclean
#
all: $(PROGRAM_NAME) 
#
debug: CPPFLAGS = -g -O0 -std=c++0x -Wall -pedantic -DDEBUG -llenstoolhpc -llenstool $(EXTRA_FLAGS)
debug: NVFLAGS = -O0 $(GEN_SM35) -g -G 
debug: NVFLAGS += $(foreach includedir,$(program_CU_INCLUDE_DIRS),-I$(includedir))
debug: $(PROGRAM_NAME)

# Rule for compilation of CUDA source (C++ source can be handled automatically)
%.cuo: %.cu %.cuh
	nvcc $(NVFLAGS) $(CPPFLAGS) -o $@ -dc $<

%.cpp: %.cpp %.h
	$(CXX) $(CPPFLAGS)  -o $@ $< $(CPPFLAGS) 

# This is pretty ugly...details below
# The program depends on both C++ and CUDA Objects, but storing CUDA objects as .o files results in circular dependency
# warnings from Make. However, nvcc requires that object files for linking end in the .o extension, else it will throw
# an error saying that it doesn't know how to handle the file. Using a non .o rule for Make and then renaming the file 
# to have the .o extension for nvcc won't suffice because the program will depend on the non .o file but the files in
# the directory after compilation will have a .o suffix. Thus, when one goes to recompile the code all files will be
# recompiled instead of just the ones that have been updated. 
#
# The solution below solves these issues by silently converting the non .o suffix needed by make to the .o suffix 
# required by nvcc, calling nvcc, and then converting back to the non .o suffix for future, dependency-based 
# compilation.
$(PROGRAM_NAME): $(OBJECTS) 
	$(CXX) -o $@ $(program_CXX_OBJS) $(CUO_O_OBJECTS) $(program_C_OBJS) $(CPPFLAGS) 





#$(PROGRAM_NAME): $(OBJECTS) 
#	@ for cu_obj in $(program_CU_OBJS); \
#	do				\
#		mv $$cu_obj $$cu_obj.o; \
#	done				#append a .o suffix for nvcc
#	nvcc $(NVFLAGS) $(CPPFLAGS) -o $@ $(program_CXX_OBJS) $(program_C_OBJS) $(CUO_O_OBJECTS)
#	@ for cu_obj in $(CUO_O_OBJECTS); 	\
#	do					\
#		mv $$cu_obj $${cu_obj%.*};	\
#	done				#remove the .o for make

clean:
	@- $(RM) $(PROGRAM_NAME) $(OBJECTS) *~ *.o *.optrpt 

distclean: clean

