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

cd ../Benchmarks/MappingBenchmark
make -f Makefile.GPU clean
make -f Makefile.GPU
rm -r tmp

#./Bayesmap_GPU m1931.par T
#./GridGradient_GPU ../ConfigFiles/TestResultParameter.par T
./GridGradient_GPU ../ConfigFiles/90Pot81.par T
#./ChiBenchmark_GPU ../ConfigFiles/MarkusBenchmark.par T
#./ChiBenchmark_GPU ../ConfigFiles/MarkusBenchmark1SIS.par T

#cd ../Benchmarks/GradientBenchmark
#make -f Makefile clean
#make -f MakefileDouble
#./GradientBenchmark

echo "Finish Double test"

