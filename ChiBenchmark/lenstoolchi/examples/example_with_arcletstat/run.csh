#!/bin/tcsh
# Script(!) for running the lenstool on the nfw and gnfw arcletstat examples
# Sand and Marshall (UCSB) 2007-11-30
#
# Sole argument is string 'nfw' or 'gnfw', eg

#   run.csh gnfw

# Results to compare agasint should be in "ref_?nfw" directories...


if ($#argv == 0) then
  echo "Usage: run.csh [nfw/gnfw]"
  goto FINISH
else if ($argv[1] == 'nfw') then  
  set parfile = nfw.par
  set subdir = results_nfw
else if ($argv[1] == 'gnfw') then 
  set parfile = gnfw.par
  set subdir = results_gnfw
else
  echo "ERROR: unrecognised run name '$argv[1]'"
  echo "Usage: run.csh [nfw/gnfw]"
  goto FINISH
endif

echo "Running lenstool on $parfile in subdirectory ./$subdir..."
set datafile = nfw_arclet.cat
  
mkdir -p $subdir
chdir $subdir

ln -sf ../$datafile .
ln -sf ../$parfile .

lenstool $argv[1] -n 

FINISH:
