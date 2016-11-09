#!/usr/bin/perl
# Remove the Evidence line and Evidence column in old bayes.dat
# files.
#

if( @ARGV > 0 )
{
	$file = $ARGV[0];
}
else
{
	$file = "bayes.dat";
}

open(IN, $file) || die "$file not found.";

while( $line = <IN> )
{
	if( $line =~ "#" )
	{
		print $line  if( $line ne "#Evidence\n" );
	} 
	else
	{
		@fld = split " ", $line;
		printf "%s\n", join( " ", @{ [ @fld[0..1], @fld[3..$#fld] ] } );
	}
}

