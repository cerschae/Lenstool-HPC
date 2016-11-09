#!/usr/bin/perl
# Plot a cleanlens contour file on an open ds9 FITS file. 
# Contour coordinates are in relative arcsec to the reference frame in lenstool .par file
# Add the clean argument to delete all regions before plotting the polygon
#
# Syntax : pcontour.pl <.par> <contFile> [clean]
#

# Analyse arguments
$par=$ARGV[0];
$contFile=$ARGV[1];

if( scalar @ARGV < 2 )
{
    print "Syntax: pcontour.pl <.par> <contour file>\n";
    exit 1;
}

# Extract reference RA and DEC from .par
open( par, $par) || die "File not found $par\n";
@sref = grep {/reference/} <par>;
close par;

chop $sref[0];
($dump, $iref, $ra, $dec) = split " ",$sref[0];
print "REFERENCE in $par  $iref $ra $dec\n";

# Process contour file
open( contf, $contFile) || die "File not found $contFile\n";
while( <contf> )
{
	chop;
	($i, $x, $y) = split;
	print "$x $y\n";
	# $x and $y are in relative arcsec --> conversion to WCS
	$x = $x / -3600. / cos( $dec / 180. * 3.1415927 ) + $ra;
	$y = $y / 3600. + $dec;
	push @polyx, $x;
	push @polyy, $y;
}
close contf;

# Write the ds9 file
open( ds9, ">e.reg");
printf ds9 "fk5; polygon(%f,%f",$polyx[0],$polyy[0];
for( $i=1; $i <= $#polyx; $i++ )
{
	printf ds9 ",$polyx[$i],$polyy[$i]";
}
printf ds9 ")";
close ds9;

# Display e.reg in ds9
$title = "ds9";
$title = $ENV{'DS9'} if( $ENV{'DS9'} );
system( "xpaset -p $title regions delete all" ) if( $ARGV == "clean" );
system( "cat e.reg | xpaset $title regions" );
