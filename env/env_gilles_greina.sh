#echo "\$BASH_SOURCE ${BASH_SOURCE[@]}"
#echo "$PWD/${BASH_SOURCE[@]}"
#echo dirname "$PWD/${BASH_SOURCE[@]}"
#echo $(realpath "$(dirname "$PWD/${BASH_SOURCE[@]}")")
MYPATH=$(realpath "$(dirname "$PWD/../${BASH_SOURCE[@]}")")
#echo $MYPATH
export LENSTOOL_ROOT=/scratch/fgilles/Projects/lenstool-6.8.1/
#export LENSTOOL_ROOT=/users/fgilles/Projects/lenstool-6.8.1/src/
#export LENSTOOLHPC_ROOT=/users/fgilles/GPU-Projects/lenstool-hpc-master/
export LENSTOOLHPC_ROOT=$MYPATH
export CFITSIO_ROOT=/scratch/fgilles/Projects/Libs/cfitsio/
export WCSTOOL_ROOT=/scratch/fgilles/Projects/Libs/wcstools-3.9.4/libwcs/
export GSL_ROOT=/scratch/fgilles/Projects/Libs/gsl-2.2/
#
export LD_LIBRARY_PATH=$LENSTOOLHPC_ROOT/src:$GSL_ROOT/lib:$LD_LIBRARY_PATH 
#export PATH=/cm/shared/apps/intel-2017/vtune_amplifier_xe/bin64/:$PATH
#
module load intel/compiler intel/mpi intel/mkl
export CUDA_ROOT=/usr/local/cuda/
export PATH=$CUDA_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_ROOT/lib64:$LD_LIBRARY_PATH
#
#module load cuda92
#
export CXX=mpiicpc
