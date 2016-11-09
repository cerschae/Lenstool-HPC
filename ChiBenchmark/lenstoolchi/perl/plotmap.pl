#!/usr/bin/perl

# Plot a map.iso file produced by lenstool.

# Input/Output : NONE
# Requirement : PGPLOT for perl

use PGPLOT; # Load PGPLOT module

$\ = "\n";

open(in, "map.iso") or die "ERROR : map.iso doesn't exists\n";

$n1 = $n2 = 10;
$i=0; $j=0;
while(<in>) 
{
	if( $_ =~ /Nrcut:(\d+)/ ) {
		$n1 = $1;
	}
	if( $_ =~ /Nsig:(\d+)/ ) {
		$n2 = $1;
	}
	if( $_ =~ /^([\d|\.]*) ([\d|\.]*) ([\d|\.]*)/ )
	{	
		($$rcut[$i][$j], $$sig[$i][$j], $$chi2[$i][$j]) = ($1, $2, $3);
#		printf "%d %d %f %f %f\n",  $i,$j,$$rcut[$i][$j], $$sig[$i][$j], $$chi2[$i][$j];
		$i++;

		if( $i == $n1 ) { 
			$i = 0; $j++;
		}
	}
}
printf "%d x %d points read\n", $n1, $n2;
close(in);

# Compute the chi2 min and max values
$chi2min = 1e20;
$chi2max = -1e20;
$i=0; $j=0;
while( $j < $n2 )
{
	if( $chi2min > $$chi2[$i][$j] ) {
		$chi2min = $$chi2[$i][$j];
		$imin = $i; $jmin = $j;
	}
	if( $chi2max < $$chi2[$i][$j] ) {
		$chi2max = $$chi2[$i][$j];
		$imax = $i; $jmax = $j;
	}
	$i++;

	if( $i == $n1 ) {
		$i = 0; $j++;
	}
}

printf "chi2min:%f chi2max:%f \n", $chi2min, $chi2max;
printf "best_rcut:%.2f best_sig:%.2f chi2:%.3f\n", $$rcut[$imin][$jmin], $$sig[$imin][$jmin], $chi2min;

# Set the min and max for rcut and sigma
$rcutmin = $$rcut[0][0];
$rcutmax = $$rcut[$n1-1][$n2-1];
$sigmin = $$sig[0][0];
$sigmax = $$sig[$n1-1][$n2-1];

# Define the transformation matrix array coordinates to (rcut, sigma) for pgcont
# rcutval = rcutmin + (rcutmax - rcutmin)/(n1-1) * (I-1)
# sigval = sigmin + (sigmax - sigmin)/(n2-1) * (J-1)
$tr[0] = $rcutmin;
$tr[2] = 0.;
$tr[1] = ( $rcutmax - $rcutmin ) / ($n1 - 1);
$tr[0] -= $tr[1];
$tr[3] = $sigmin;
$tr[5] = ( $sigmax - $sigmin ) / ($n2 - 1);
$tr[3] -= $tr[5];
$tr[4] = 0.;

# Set the 3 chi2 contour levels (3, 2, 1)sig
$alev[0] = $chi2min + 3 / 10 * ($chi2max - $chi2min);
$alev[1] = $chi2min + 2 / 10 * ($chi2max - $chi2min);
$alev[2] = $chi2min + 1 / 10 * ($chi2max - $chi2min);

# Begin plot section
pgbegin(0, "/xserve", 1, 1); # Open plot device
pgscf(2); # Set character font
pgslw(4); # Set line width
pgsch(1.6); # Set character height

# plot the viewport
pgenv( $rcutmin, $rcutmax, $sigmin, $sigmax, 0, 0 );
pglab( "rcut (kpc)", "sigma (km/s)", "(white, red, green) = (3sig, 2sig, 1sig)");
pggray( $chi2, $n1, $n2, 1, $n1, 1, $n2, $chi2min, $chi2max, \@tr);
pgsci(1);
pgbox('BCNTS',0.,0,'BCNTS',0.,0);

# plot the 3 contours with 3 different colors
pgsci(1);
pgcont( $chi2, $n1, $n2, 1, $n1, 1, $n2, $alev[0], 1, \@tr);
pgsci(2);
pgcont( $chi2, $n1, $n2, 1, $n1, 1, $n2, $alev[1], 1, \@tr);
pgsci(3);
pgcont( $chi2, $n1, $n2, 1, $n1, 1, $n2, $alev[2], 1, \@tr);

# Plot the min(blue) and max(red) chi2 position
pgsci(4);
pgpt( 1, $$rcut[$imin][$jmin], $$sig[$imin][$jmin], 17);
pgsci(2);
pgpt( 1, $$rcut[$imax][$jmax], $$sig[$imax][$jmax], 17);

# interactive : plot the chi2 value when click and quit when 'q'
$ch = "";
pgsci(1);
pgmtext( 'B', 2, -0.15, 0., "chi2 :");
while( $ch ne "q" ) {
	pgcurs( $x, $y, $ch);
	pgsci(0);
	pgmtext( 'B', 2, 0., 0., $$chi2[$j][$i] );
	$i = $x + $tr[1] / 2.; $j = $y + $tr[5] / 2.;
	$i = ($i - $rcutmin) * ($n1 - 1) / ($rcutmax - $rcutmin);
	$j = ($j - $sigmin) * ($n2 - 1) / ($sigmax - $sigmin);
	pgsci(1);
	pgmtext( 'B', 2, 0., 0., $$chi2[$j][$i] );
#	printf "%f %f %d %d\n", $x, $y, $i, $j;
}

pgend;
