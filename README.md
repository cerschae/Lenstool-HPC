# Lenstool-HPC

Mass-modelling tool and fast Lens-map generation based on Lenstool using HPC techniques. 
Designed for CPU and GPU hardware clusters.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine 
for development and testing purposes.

### Prerequisites

* [Cuda toolkit](https://developer.nvidia.com/cuda-toolkit)
* [CuDNN](https://developer.nvidia.com/cudnn)
* Intel Compiler
* ([Lenstool](https://projets.lam.fr/projects/lenstool/wiki))
* Cfitsio library
* GSL library
* Wcstools library


### Installing

The Makefiles have not yet been modified to work with g++, only icpc and mpiicpc. To compile 
the lenstool-hpc library, run make in the main directory. 

To compile the executables, following environment variables have to be set:

``export LENSTOOLHPC_ROOT=/path/lenstool-hpc``
``export LD_LIBRARY_PATH+=:/path/lenstool-hpc/src``

``export CFITSIO_ROOT=/path/cfitsio``
``export WCSTOOL_ROOT=/path/libwcs``
``export GSL_ROOT=/path/gsl-2.3``
``export LD_LIBRARY_PATH+=:/path/gsl-2.3/lib``

To compile with LENSTOOL (necessary for test and benchmarks) the user has to set following 
environment variables with the path to corresponding libraries:
``export LENSTOOL_ROOT=/path/lenstool-6.8.1``

The compilation of the src creates two libraries: liblenstoolhpc and liblenstoolhpc_GPU (if 
nvcc was found) in the src folder. Compilation of the utils folder creates the executables 
for map generation. Compilation of the executables necessitates nvcc.

* Lenstool_HPC (Single map generation)
* Bayesmap_GPU (Multi map generation)





### Example

Generating a single map necessitates two input parameters: A Lenstool parameter file and a 
path to a folder for the results.

```
./Lenstool_HPC .parfile folder
```

Concrete example:

```
./Lenstool_HPC ../../Benchmarks/ConfigFiles/m1931.par Test
```
Generating multiple maps using a bayes.dat file needs: A Lenstool parameter file, a path to
 a folder for the results and the bayes.dat file in the same folder.

```
cd utils/maps
```

Concrete example:

```
./Bayesmap_GPU ../../Benchmarks/ConfigFiles/m1931.par Test
```

## Running the tests and Benchmarks

Tests with Lenstool-HPC compare the result with the corresponding Lenstool functions. A 
compiling Lenstool Library has to be provided. Tests can be run by launching the corresponding 
Script located in ``Benchmark/Scripts``. These Test double also as a Benchmark comparison between 
Lenstool and Lenstool-HPC.

* Gradient Benchmark: Testing gradient Computation test
* GridGradient Benchmark: Testing first order derivative of deflection potential computation over a grid 
* GridGradient2 Benchmark: Testing second order derivative of deflection potential computation over a grid 
* GridPotential Benchmark: Testing deflection potential computation over a grid 
* Bench*** : Benchmark over a HFF cluster


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details 
on our code of conduct, and the process for submitting pull requests to us.


## Authors

* **Christoph Schäfer** - [cerschae](https://github.com/cerschae)
* **Gilles Fourestey** [ursache](https://github.com/ursache)
* **Markus Rexroth** - [MarkusRe](https://github.com/MarkusRe)

## License

This project is licensed under the GPL3 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* SCITAS
* CSCS
