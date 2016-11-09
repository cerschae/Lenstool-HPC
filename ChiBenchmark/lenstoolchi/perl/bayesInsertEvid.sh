#!/bin/sh
# Insert a fake Evidence column into an old format bayes.dat file.
#

if( [ ! -e bayes.dat ] )
then
	echo "bayes.dat not found in current directory."
	exit
fi

if( [ ! -z "$(grep Evidence bayes.dat)" ] )
then
	echo "bayes.dat already in the new format."
	exit
fi

sed -i "/#Chi2/ a #Evidence" bayes.dat
sed -i "s/^[0-9]* /&1 /" bayes.dat
