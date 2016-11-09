#!/bin/bash
# Convert the bayes.dat header into a ASCII header readable by TOPCAT
#

if( [ $1 -a -e $1 ] )
then
	bayes=$1
else
	bayes=bayes.dat
fi

root=${bayes%%.dat}
head=$(grep "^#" $bayes | tr -s ' ' '_' | tr '\n' ' ' | tr '#' ' ' | sed 's/^\s//' )
echo "#$head" > head
grep -v "^#" $bayes > body

cat head body > $root.top
rm head body

echo "Starting topcat..."
topcat -f ascii $root.top



