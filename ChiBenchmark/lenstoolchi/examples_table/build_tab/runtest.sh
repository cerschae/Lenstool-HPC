#!/bin/sh
# Script(!) for running the lenstool_tab binary.
# Sand and Marshall (UCSB) 2007-11-29
#
# The output file is called gnfw.tab - only one of these is needed, and it
# should be kept in the lenstool/table_src directory (where lenstool will look 
# for it). 

\rm -f lenstool.tab
${LENSTOOL_DIR}/table_src/lenstool_tab tab.in

# After successful completion the stdout looks like:

# Hello - welcome to lenstool_tab
# Input file tab.in read in OK:
#  0.010000 1.700000 0.010000 0.010000 30.000000 1.010000
# Steps in alpha: 170
# Steps in x: 805
# Total steps: 136850
# Progress : alpha 1.700000/1.700000  x 29.811548/30.000000
# Total steps = 136850
# Alpha count = 170
# xcount = 805
# Table written in lenstool.tab
