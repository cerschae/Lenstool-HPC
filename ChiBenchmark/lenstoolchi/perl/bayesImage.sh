#!/bin/bash
# Call the bayesImage tool, combine all the image*.all files
# and plot then on DS9
#

NAMAX=5000              # Maximum number of images to plot
mult=bayesImage.all     # Name of the combined catalog

if [ -z $1 -o $1 == "-h" ]
then
	echo "Syntax : bayesImage.sh [OPTION] <.par>"
    echo "Available OPTIONS:"
    echo " -b : compute images of system barycenter [default]"
    echo " -i : compute counter images, image per image"
	exit 1
fi

# Analyse options
if [ ${1:0:1} == "-" ]
then
    mode=$1
    par=$2
else
    par=$1
fi

[ -e images/image0.all ] && rm images/image[0-9]*.all
trap '' INT  # Intercept the INT signal when CTRL-C is pressed
ok=0 && bayesImage $mode $par && ok=1

if [ $ok -eq 1 ] 
then
	[ -e $mult ] && rm $mult
	head -n1 images/image0.all > $mult
	cat images/image[0-9]*.all | grep -v "#" >> $mult
	if [ $(wc -l $mult | awk '{print $1}') -gt $NAMAX ] 
	then
		echo "WARNING: Too many images... catalog cut at $NAMAX images"
		head -n$((NAMAX+1)) $mult >> a
		mv a $mult
	fi
	pelli $mult red clean 100 point
#	rm image[0-9]*.all
fi
