#!/bin/bash
echo "Starting Mapping Test"

#if recompiling lenstool is needed
#cd ../Libs/lenstool-6.8.1/
#make
#cd ../../lenstool-hpc

cd ../../src/
make clean
make


cd ../Benchmarks/GradientBenchmark
make clean
make 
./GradientBenchmark

echo "Finish Double test"

