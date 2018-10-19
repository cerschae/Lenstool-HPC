#!/bin/bash -l

#set -x
#../utils/gpu  | grep "capability" | cut -d " " -f11 #| paste -s -d, - > ./gpusm.inc
#./gpu  | grep "capability" | cut -d " " -f11 | sort | uniq  > ./.gpusm.inc
$LENSTOOLHPC_ROOT/utils/gpu  | grep "capability" | cut -d " " -f11 | sort | uniq > ./.gpusm.inc
while read line; 
	do val=$(echo $line | tr -d '.') 
	str="$str -gencode arch=compute_$val,code=sm_$val " 
done < ./.gpusm.inc
echo "GENCODE += $str"
rm -f ./.gpusm.inc
#sed -ir 's/,/,sm_/g' ./gpusm.inc
#cat ./.gpusm.inc | sed -r 's/=/=sm/g' > ./gpusm.inc
#cp ./gpusm.inc ./.gpusm.inc
#cat ./.gpusm.inc | sed -r 's/=sm/=compute/g' > ./gpusm.inc
#cat ./gpusm.inc | sed -r 's/,sm/,compute/g' > ./.gpusm.inc
#cat ./.gpusm.inc >> ./gpusm.inc
#sed -i '1s/^/SM:=sm/' ./gpusm.inc
#sed -i '2s/^/ARCH:=compute/' ./gpusm.inc
#cat ./gpusm.inc | sed -r 's/,SM=/,arch=/g' > ./gpusm.inc
#rm -f ./.gpusm.inc
