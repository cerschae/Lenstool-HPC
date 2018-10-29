#
# author: christophernstrerne.schaefer@epfl.ch
#

SHELL:=/bin/bash

#LENSTOOL_DEBUG = true;

All: SRC Lenstool_HPC BayesMapGPU

SRC:
	$(MAKE) -C src/
Lenstool_HPC:
	$(MAKE) -C utils/exec
BayesMapGPU:
	$(MAKE) -C utils/maps

clean:
	+$(MAKE) clean -C src
	+$(MAKE) clean -C utils/exec
	+$(MAKE) clean -C utils/maps
	
#test: 
#	export LENSTOOL_DEBUG; \


