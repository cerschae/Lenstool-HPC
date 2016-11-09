#!/usr/bin/perl
# Read all the chires*.dat in the current directory 
# and produce a chires_mean.dat file with the mean
# rms for each system.
# 

use warnings;
use strict;

sub readChires {
	my ($dir, @chires);
	opendir($dir, "chires") || die "Cannot open chires/ directory\n";
	@chires = grep { /chires[0-9]*\.dat/ } readdir($dir);
	closedir $dir;
	return @chires;
}

#---------------------------------------------------------------------- 
# MAIN PROGRAM
#---------------------------------------------------------------------- 

if( @ARGV < 1 )
{
	print "Syntax: bayesChires.pl <.par>\n";
	exit -1;
}

# Delete chires[0-9]*dat files in the current directory
#@chires = readChires();
#printf "INFO: Delete %d chires[0-9]*dat files\n", unlink @chires;

# Create chires[0-9]*dat files and process them
my $status = system("bayesChires $ARGV[0]");
die "ERROR: bayesChires error\n" if $status;

my @chires = readChires();
if( @chires eq 0 )
{
	print "ERROR: No chires/chires*.dat files found\n";
	exit 1;
}

# Read the chires*.dat files
my (@n, @id, @z, @narc, @chip, @chix, @chiy, @chia);
my (@n_cl, @x_cl, @y_cl, @chil);
my (@n_s, @id_s, @z_s, @chis, @gamma, @kappa, @g, @es1, @es2);
my (@rmss2, @rmsi2, @dx, @dy, @nwarn);
my ($i, $l, $l_cl, $l_s);
my ($secmult, $seccrit, $secwl);
my $nchires = scalar @chires;
for( @chires )
{
	$_ = "chires/".$_;
    my $ifh;
	open( $ifh,  $_ );
	$l = 0; # line index for multiple images in the chires.dat file
	$l_cl = 0; # line index for critical lines in the chires.dat file
	$l_s = 0; # line index for arclets in the chires.dat file
    $secmult = 0;
    $seccrit = 0;
    $secwl = 0;
LOOP:	while( <$ifh> )
	{
		chop;
        # Check the section limits
        if( /chi multiples/ )
        {
            $secmult = 1;
            my $nextline = <$ifh>;   # skip line N ID z ...
            next LOOP;
        }
        if( /chimul/ )
        {
            $secmult = 0;
            next LOOP;
        }
        if( /chi critical/ )
        {
            $seccrit = 1;
            next LOOP;
        }
        if( /chil/ )
        { 
            $seccrit = 0;
            next LOOP;
        }
        if( /chis arclets/ )
        {
            $secwl = 1;
            my $nextline = <$ifh>;   # skip line N ID z ...
            next LOOP;
        }
        if( /chis/ )
        { 
            $secwl = 0;
            next LOOP;
        }

        # Process a line in section chi multiple
		my @fld = split;
        if( $secmult )
        {
    		$n[$l] = $fld[0];
    		$id[$l] = $fld[1];
    		$z[$l] += $fld[2];
    		$narc[$l] = $fld[3]; # always the same number
    		$chip[$l] += $fld[4];
    		$chix[$l] += $fld[5];
    		$chiy[$l] += $fld[6];
    		$chia[$l] += $fld[7];
        	$rmss2[$l] += $fld[8]*$fld[8]; #total rmss
          	$rmsi2[$l] += $fld[9]*$fld[9]; #total rmsi
            if( $narc[$l] == 1 ) 
            {
                $dx[$l] += $fld[10];
                $dy[$l] += $fld[11];
            }
    		$nwarn[$l] += $fld[12];
    		$l++; # next system
        }

        # Process a line in section chi critical line
        if( $seccrit )
        {
            $n_cl[$l_cl] = $fld[0];
            $x_cl[$l_cl] = $fld[1];
            $y_cl[$l_cl] = $fld[2];
            $chil[$l_cl] += $fld[3];
            $l_cl++; # next critical line
        }

        # Process a line in section chi arclets
        if( $secwl )
        {
            $n_s[$l_s] = $fld[0];
            $id_s[$l_s] = $fld[1];
            $z_s[$l_s] = $fld[2];
            $chis[$l_s] += $fld[3];
            $gamma[$l_s] += $fld[4];
            $kappa[$l_s] += $fld[5];
            $g[$l_s] += $fld[6];
            $es1[$l_s] += $fld[7];
            $es2[$l_s] += $fld[8];
            $l_s++;  # next arclet
        }

	}
    close($ifh);
}
# Write the chires_mean.dat file
open(OUT, ">chires_mean.dat");

# Print summary for arclets
my $chistot = 0.;
if( $l_s > 0)
{
    print OUT "chis arclets\n";
    print OUT "     N       ID      z        chi2      gamma         1-k          g            es1          es2\n";
    for( $i = 0; $i < $l_s; $i ++ )
    {
        $chis[$i] /= $nchires;
        $gamma[$i] /= $nchires;
        $kappa[$i] /= $nchires;
        $g[$i] /= $nchires;
        $es1[$i] /= $nchires;
        $es2[$i] /= $nchires;
        $chistot += $chis[$i];
        printf OUT "%6ld   %6s    %.3lf    %6.4lf    %6.3le    %6.3le    %6.3le    %6.3le    %6.3le\n", $n_s[$i], $id_s[$i], $z_s[$i], $chis[$i], $gamma[$i], $kappa[$i], $g[$i], $es1[$i], $es2[$i];
    }
    printf OUT "chis               %.2f\n", $chistot;
    print OUT "\n";
}

# Strong lensing summary
my $chiptot = 0.;
if ( $l > 0 )
{
    print OUT "chi multiples\n";
    print OUT " N    ID    z   Narcs    chip    chix    chiy    chia   rmss     rmsi    dx       dy    nwarn\n";
    my ($chixtot, $chiytot, $chiatot);
    my ($rmss2tot, $rmsi2tot, $nwarntot, $nsys); 
    my ($sdx, $sdy);
    $chiptot=$chixtot=$chiytot=$chiatot=0.;
    $rmss2tot=0;
    $rmsi2tot=0;
    $nwarntot=0;
    $nsys=0;
    for( $i = 0; $i < $l; $i++ )
    {
        my ($rmss, $rmsi);
    	$rmss = $rmsi = 0.;
        $sdx = $sdy = "  N/A ";
    
    	my $name=$id[$i];
        $z[$i] /= $nchires;
    	$rmss = $rmss2[$i] / $nchires;
    	$chip[$i] /= $nchires;
    	$chix[$i] /= $nchires;
    	$chiy[$i] /= $nchires;
    	$chia[$i] /= $nchires;
    
    	# Be careful with non predicted images
        if( $nchires - $nwarn[$i] > 0 )
        {
    		$rmsi = $rmsi2[$i] / ($nchires - $nwarn[$i]);
        	if( $narc[$i] == 1 )
            {
                $sdx = sprintf("%6.2f", $dx[$i] / ($nchires - $nwarn[$i]));
                $sdy = sprintf("%6.2f", $dy[$i] / ($nchires - $nwarn[$i]));
        	}	
        }
    
    	$name = $name.'*' if( $nwarn[$i] > 0 );
    
        # print a line for an image
    	printf OUT "%2d %6s %.3f   %d   %7.2f  %6.2f  %6.2f  %6.2f  %6.3f  %6.2f  %s   %s  %d\n", $n[$i], $name, $z[$i], $narc[$i], $chip[$i],
    	       	$chix[$i], $chiy[$i], $chia[$i], sqrt($rmss), sqrt($rmsi), $sdx, $sdy, $nwarn[$i];
    
    	if( $narc[$i] != 1 )
    	{
    		$rmss2tot += $rmss;
    		$rmsi2tot += $rmsi;
    		$chiptot += $chip[$i];
    		$chixtot += $chix[$i];
    		$chiytot += $chiy[$i];
    		$chiatot += $chia[$i];
            $nwarntot += $nwarn[$i];
    		$nsys++;
    	}
    }
    $rmss2tot /= $nsys;
    $rmsi2tot /= $nsys;
    printf OUT "chimul                %7.2f  %6.2f  %6.2f  %6.2f  %6.3f  %6.2f    N/A      N/A   %d\n", $chiptot, $chixtot, $chiytot, $chiatot, sqrt($rmss2tot), sqrt($rmsi2tot),$nwarntot;
}

# Print summary for critical lines
my $chiltot = 0.;
if( $l_cl > 0)
{
    printf OUT "chi critical line\n";
    for( $i = 0; $i < $l_cl; $i ++ )
    {
        $chil[$i] /= $nchires;
        $chiltot += $chil[$i];
        printf OUT "%d  %.2f %.2f %.2f\n", $n_cl[$i], $x_cl[$i], $y_cl[$i], $chil[$i];
    }
    printf OUT "chil %.2f\n", $chiltot;
}

printf OUT "\nchitot           %.2f\n", $chiltot + $chiptot + $chistot;

close(OUT);

# Delete chires[0-9]*dat files
#@chires = readChires();
#printf "INFO: Delete %d chires[0-9]*dat files\n", unlink @chires;
