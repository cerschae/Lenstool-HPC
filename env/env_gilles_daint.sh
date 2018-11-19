#echo "\$BASH_SOURCE ${BASH_SOURCE[@]}"
#echo "$PWD/${BASH_SOURCE[@]}"
#echo dirname "$PWD/${BASH_SOURCE[@]}"
#echo $(realpath "$(dirname "$PWD/${BASH_SOURCE[@]}")")
MYPATH=$(realpath "$(dirname "$PWD/../${BASH_SOURCE[@]}")")
echo $MYPATH
export LENSTOOL_ROOT=/scratch/snx3000/fgilles/Lenstool/lenstool-6.8.1/
#export LENSTOOL_ROOT=/users/fgilles/Projects/lenstool-6.8.1/src/
#export LENSTOOLHPC_ROOT=/users/fgilles/GPU-Projects/lenstool-hpc-master/
export LENSTOOLHPC_ROOT=$MYPATH
export CFITSIO_ROOT=/scratch/snx3000/fgilles/Lenstool/Libs/cfitsio/
export WCSTOOL_ROOT=/scratch/snx3000/fgilles/Lenstool/Libs/wcstools-3.9.4/libwcs/
export GSL_ROOT=/scratch/snx3000/fgilles/Lenstool/Libs/gsl-2.2/
export CUDAROOT=$CUDATOOLKIT_HOME
module swap PrgEnv-cray PrgEnv-intel

export LD_LIBRARY_PATH=$LENSTOOLHPC_ROOT/src:$LD_LIBRARY_PATH 
export LD_LIBRARY_PATH=$GSL_ROOT/lib:$LD_LIBRARY_PATH
#source /cm/shared/apps/intel-2017/compilers_and_libraries_2017/linux/mpi/intel64/bin/mpivars.sh
#module load intel/2017
#export PATH=/cm/shared/apps/intel-2017/vtune_amplifier_xe/bin64/:$PATH

#module load intel-compilers/2016.3.210 impi cuda80
#source /cm/shared/apps/INTEL/2016/impi/5.1.3.210/bin64/mpivars.sh

module load daint-gpu
module swap intel/18.0.2.199 intel/17.0.4.196
module load cudatoolkit/9.1.85_3.18-6.0.7.0_5.1__g2eb7c52


export CXX=CC
