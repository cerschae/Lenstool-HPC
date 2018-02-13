#
# author: gilles.fourestey@epfl.ch
#

SHELL:=/bin/bash

arch = intel
prec = double
all: precision compile


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
                make --silent -f Makefile.GPU.$(arch); \
        fi  
	

clean:
	make -f Makefile.$(arch) clean
