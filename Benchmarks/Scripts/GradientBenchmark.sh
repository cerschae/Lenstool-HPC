#!/bin/bash
echo "Starting Mapping Test"

#if recompiling lenstool is needed
#cd ../Libs/lenstool-6.8.1/
#make
#cd ../../lenstool-hpc

cd ../../src/
make -f Makefile.intel clean
make -f Makefile.intel
make -f Makefile.GPU.intel

cd ../Benchmarks/GradientBenchmark
make -f Makefile clean
make -f MakefileDouble
./GradientBenchmark

echo "Finish Double test"

