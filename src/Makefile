#
# author: gilles.fourestey@epfl.ch
#

SHELL:=/bin/bash

arch = intel
prec = double
all: utils precision compile

#cd ../utils/
#make

utils:
	make -C ../utils
	sh ../utils/gpu.sh > ./gpusm.inc
	

precision:
	echo "PRECISION := -D_double" > precision.inc   
	if [ "$(prec)" == "single" ]; then \
		echo "PRECISION := -D_single" > precision.inc; \
	fi 
#fi

compile:
	@make -f Makefile.$(arch)
	if which nvcc > /dev/null; then \
		echo "nvcc detected, compiling the GPU lib"; \
                make -f Makefile.GPU.$(arch); \
	fi  
	

clean:
	make -C ../utils clean
	make -f Makefile.$(arch) clean
                
#make --silent -f Makefile.GPU.$(arch); 
