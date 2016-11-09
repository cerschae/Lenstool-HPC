#!/usr/bin/perl
# From a ds9 FITS file, read a ds9 region, analyse it and write
# a contour file usable by cleanlens in lenstool
# If contFile is not provided, print on STDOUT
#
# syntax : gcontour.pl <.par> [<contFile>]
#

# Analyse arguments
$par=$ARGV[0];
$contFile = $ARGV[1];

if( ! @ARGV )
{
    print "Syntax: gcontour.pl <.par>\n";
    exit 1;
}

# Extract reference RA and DEC from .par
open( par, $par) || die "File not found $par\n";
@sref = grep {/reference/} <par>;
close par;

chop $sref[0];
($dump, $iref, $ra, $dec) = split " ", $sref[0];

# Extract the polygon coordinates from ds9 FITS file polygon region
$title = "ds9";
$title = $ENV{'DS9'} if( $ENV{'DS9'} );
system( "xpaset -p $title regions system wcs; xpaget $title regions selected > e.reg" );
open( reg, "e.reg" ) || die "File not found e.reg\n";
while( <reg> )
{
	if( /polygon\((.*)\)/ )
	{
		@coords = split /,/, $1;
	}
}
close reg;

# Write the contour file for cleanlens
if( $contFile )
{
	open( $contf, ">$contFile" );
}
else
{
	$contf = *STDOUT;
}

for( $i = 0; $i <= $#coords; $i+=2 )
{
	$x = $coords[$i] - $ra;
	$x *= -3600. * cos( $dec / 180. * 3.1415927 );
	$y = $coords[$i+1] - $dec;
	$y *= 3600.;
	printf $contf "%d  %.2f  %.2f\n", $i/2+1, $x, $y;
}
