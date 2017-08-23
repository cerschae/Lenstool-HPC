#echo "\$BASH_SOURCE ${BASH_SOURCE[@]}"
#echo "$PWD/${BASH_SOURCE[@]}"
#echo dirname "$PWD/${BASH_SOURCE[@]}"
#echo $(realpath "$(dirname "$PWD/${BASH_SOURCE[@]}")")
MYPATH=$(realpath "$(dirname "$PWD/${BASH_SOURCE[@]}")")
#echo $MYPATH
export LENSTOOL_ROOT=/users/fgilles/Projects/lenstool-6.8.1/
#export LENSTOOL_ROOT=/users/fgilles/Projects/lenstool-6.8.1/src/
#export LENSTOOLHPC_ROOT=/users/fgilles/GPU-Projects/lenstool-hpc-master/
export LENSTOOLHPC_ROOT=$MYPATH
export CFITSIO_ROOT=/users/fgilles/Projects/Libs/cfitsio/
export WCSTOOL_ROOT=/users/fgilles/Projects/Libs/wcstools-3.9.4/libwcs/
export GSL_ROOT=/users/fgilles/Projects/Libs/gsl-2.2/

export LD_LIBRARY_PATH=$LENSTOOLHPC_ROOT/src:$LD_LIBRARY_PATH 
module load intel-compilers/2016.3.210 cuda80
