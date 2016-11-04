PROGRAM_NAME := GradientTest


program_CXX_SRCS := $(wildcard *.cpp)
program_CXX_SRCS += $(wildcard ../GradientBenchmark/*.cpp) #Find C++ source files from additonal directories
program_CXX_SRCS := $(filter-out ../GradientBenchmark/main.cpp, $(program_CXX_SRCS))

program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}

#program_C_SRCS := $(wildcard ../GradientBenchmark/*.cpp)
#program_CXX_SRCS += $(wildcard ../../*.c) #Find C source files from additonal directories
#program_C_OBJS := ${program_C_SRCS:.c=.o}

#program_INCLUDE_DIRS += /usr/local/cuda/include/ #C++ Include directories
program_INCLUDE_DIRS := . /usr/include/cfitsio/ 
program_INCLUDE_DIRS += ../GradientBenchmark/
#program_CU_INCLUDE_DIRS := /home/users/amclaugh/CUB/cub-1.3.2/ #CUDA Include directories

program_INCLUDE_LIBS := /usr/lib64/ #Include libraries

# Compiler flags
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CPPFLAGS += $(foreach includelib,$(program_INCLUDE_LIBS),-L$(includelib) -lcfitsio)
CXXFLAGS += -g -O3 -std=c++0x -Wall -pedantic


OBJECTS = $(program_CXX_OBJS) $(program_C_OBJS)

.PHONY: all clean distclean

all: $(PROGRAM_NAME) 

$(PROGRAM_NAME): $(OBJECTS) 
	icpc $(CPPFLAGS) -o $@ $(program_CXX_OBJS) $(program_C_OBJS)
	

clean:
	rm *.o

