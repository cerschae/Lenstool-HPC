#!/bin/bash
# Compute and display the counter images of the selected 
# ellipse in DS9
# 
# How it works:
# - select an ellipse with text in DS9
# - run the script
#
# Option: You can specify the source redshift by typing 
#  counterImage.sh -z<redshift> 
#
# and again
#  counterImage.sh 
# 
# to effectively compute the counter images.
# 
# Be aware :
# - Check that your best.par file is in the working
#   directory of DS9
# 
# EJ 09/2007
#

##############################################
# ADAPT THIS VALUE TO THE SIZE OF YOUR FIELD
##############################################
dmax=40   # size in arcsec of the champ section around the image


imcat=counter.cat
mod=best.par
sreg=s.reg
oreg=o.reg
zfile=z
wdir=$(xpaget ds9 cd)
z=1 #default redshift

cd $wdir

#############################################
# FIRST FUNCTION: set the image redshift with -z option
##############################################
#
# Write the redshift in z file and exit
if( [ "x${1:0:2}" = "x-z" ] )
then
	z=${1##-z}
	echo $z >| $zfile
	exit
fi

#############################################
# SECOND FUNCTION: compute counter images
##############################################
#
# Save the selected ellipses in DS9 format file $sreg, and count them
rm -f e.reg $sreg $oreg
xpaget ds9 regions selected > $sreg
nsel=$(grep "ellipse" $sreg | wc -l)
if( [ $nsel -lt 1 ] )
then
	echo "ERROR: you must select at least one ellipse"
	exit
fi
# Keep unrelated regions in $oreg
xpaset -p ds9 regions group counterImage delete
xpaset -p ds9 regions delete select
xpaget ds9 regions > $oreg
# Clean everything
xpaset -p ds9 regions delete all

# Create $imcat file (Input catalog of images in LT format)
# Write the header
rm -f $imcat
radec=$(awk '/\treference.*/ {print $3,$4}' $mod)
echo "#REFERENCE 3 $radec" > $imcat

# Convert DS9 format file $sreg to LT format file $imcat
if( cat $sreg | garc.pl STDIN $imcat > /dev/null )
then
	if( [ -e $zfile ] )
	then 
		z=$(cat $zfile)
	fi

	# detect loop format 'z0;z1;dz' in $zfile
	if( [ $(echo $z | awk -F';' '{print NF}') -eq 3 ] )
	then
		z0=$(echo $z | cut -d';' -f1)
		z1=$(echo $z | cut -d';' -f2)
		dz=$(echo $z | cut -d';' -f3)
	else
		z0=$(echo $z | cut -d';' -f1)
		z1=$z0
		dz=1.
	fi

	# Setup z0, z1 and dz to fit the for integer format
	s0=$(echo "scale($z0)" | bc)
	s1=$(echo "scale($z1)" | bc)
	ds=$(echo "scale($dz)" | bc)
	smax=$s0
	[ $s1 -gt $smax ] && smax=$s1
	[ $ds -gt $smax ] && smax=$ds
	z0=$(echo "scale=0; $z0*10^$smax/1." | bc)
	z1=$(echo "scale=0; $z1*10^$smax/1." | bc)
	dz=$(echo "scale=0; $dz*10^$smax/1." | bc)

	# Loop over the redshift range in $zfile
	for((i=$z0;i<=$z1;i+=$dz))
	do
		z=$(echo "scale=$smax; $i/(10^$smax)" | bc)
		echo "Redshift : $z"
		# set the redshift of the selected image in $imcat
		awk -v z=$z '/#/; ! /#/ {$7=z; print};' $imcat > a
	        mv a $imcat
	
		# modify the $mod file
		sed -i "s/\timage.*/\timage 3 $imcat/" $mod
	
		# Focus the champ section around the image
		xc=( $(awk '! /#/ {print $2}' $imcat) )
		yc=( $(awk '! /#/ {print $3}' $imcat) )
		xmin=$(echo "${xc[0]}-$dmax" | bc)
		xmax=$(echo "${xc[0]}+$dmax" | bc)
		ymin=$(echo "${yc[0]}-$dmax" | bc)
		ymax=$(echo "${yc[0]}+$dmax" | bc)
		sed -i "s/\txmin.*/\txmin  $xmin/" $mod
		sed -i "s/\txmax.*/\txmax  $xmax/" $mod
		sed -i "s/\tymin.*/\tymin  $ymin/" $mod
		sed -i "s/\tymax.*/\tymax  $ymax/" $mod
		echo -n "Run lenstool $mod..." 
		lenstool $mod -n  > /dev/null
		echo "done"
	
		# Convert image.all to e.reg and append counterImage tag
		pelli image.all green clean > /dev/null
		xpaset -p ds9 regions deleteall
		sed -i -e "s/$/ tag={counterImage}/" e.reg
		sed -i -e "s/text={\([0-z]*\)}/text={\1$z}/" e.reg

		# Plot the selected, unrelated and counter images back
		cat $sreg | xpaset ds9 regions
		xpaset -p ds9 regions select all
		cat $oreg | xpaset ds9 regions
		cat e.reg | xpaset ds9 regions
		[ $(($i+$dz)) -le $z1 ] && sleep 2
	done
fi
