#!/usr/bin/perl
# Take catalog and convert the coordinates from absolute to relative

$PI = 3.141592653589793238462643;
$\ = "\n";

open(IN, $ARGV[0]);

if( $ARGV[1] )
{
	open($out, '>', $ARGV[1] );
}
else
{
	$out = *STDOUT;
}

LINE: while(<IN>)
{
	chop;
	if( /#REFERENCE/ )
	{
		($null,$iref,$ra,$dec) = split;
		printf $out "#REFERENCE 3 $ra $dec\n" if( $iref == 0);
		printf $out "#REFERENCE 0 $ra $dec\n" if( $iref == 3);
		next LINE;
	} 

	@fld = split;
	if( $iref == 0 )
	{
		if( $fld[1] < $ra + 180 ) {$fld[1] = $fld[1] - $ra;}
		else {$fld[1] = $fld[1] - 360 - $ra; }
		$fld[2] -= $dec;
		$fld[2] *= 3600.;
		$fld[1] *= -3600.*cos($dec*$PI/180.);
	}
	elsif( $iref == 3 )
	{
		$fld[1] /= -3600.*cos($dec*$PI/180.);
		$fld[1] += $ra;
		$fld[2] = $fld[2]/3600. + $dec
	}

	$fld[1] = sprintf( "%.7f", $fld[1]);
	$fld[2] = sprintf( "%.7f", $fld[2]);
	printf $out "%s\n", join(" ",@fld);
}
