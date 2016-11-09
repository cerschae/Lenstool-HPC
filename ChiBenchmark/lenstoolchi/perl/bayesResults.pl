#! /usr/bin/perl
# Read the bayes.dat file and print on STDOUT the parameter mean values and
# standard deviation.
# 
# Syntax  : bayesResults.pl [bayes.dat] [scale=<arcsec_to_kpc>]
#

use Term::ANSIColor;
use IO::Handle;

# The 3 sigma errors corresponds to 3 times the std deviation.   

$PI = 3.1415927;

do "$ENV{'LENSTOOL_DIR'}/perl/median.pl";
do "$ENV{'LENSTOOL_DIR'}/perl/histogram.pl";
do "$ENV{'LENSTOOL_DIR'}/perl/rescaleBayes.pl";

$file = "bayes.dat";
$scale = 1;
while( @ARGV > 0 )
{
	if( $ARGV[0] =~ /scale=/ )
	{
		@fld = split /=/,$ARGV[0];
		$scale = $fld[1];
	} 
	else
	{
		$file = $ARGV[0];
	}
	shift @ARGV;
}

if( $scale != 1. )
{
	rescaleBayes($file, $scale);
	$file = "bayes_rescaled.dat";
}

open(bayes, $file) || die "ERROR : $file not found.\n";

$nVal=0;
$bestline=0;
$bestchi2=1e100;
$lhood=0;  # chi2 mode by default
while($line = <bayes>) 
{
	$line =~ s/\r$//;   # change DOS to UNIX fileformat
	chop($line);
	if( $line =~ /#/ )
	{
		push @param, $line;	
		$lhood=1 if( $line =~ /Lhood/ );
	} else
	{
		@fld = split / /, $line;
		
		# Find the best chi2
		$tmp = (($lhood==1)?-1:1)*$fld[1];
		if( $bestchi2 > $tmp )
		{
			$bestchi2 = $tmp; 
			$bestline = $nVal;
		}

		for( $i = 0; $i <= $#fld ; $i++ )
		{
			$$values[$i][$nVal] = $fld[$i];
		}
		$nVal++;
	}
}
close(bayes);
$nParam = $#param + 1;
printf "Read $nParam columns and $nVal lines\n";
printf "Param (nbin): <median>   <best>  <mode> <gausserr> <asymerr> (68%)\n";

for( $i = 1; $i < $nParam; $i++ )
{
	# Set the best value
	$best = $$values[$i][$bestline];

	# Compute the mean value
	$mean = 0;
	for( $j = 0 ; $j < $nVal ; $j++ )
	{
		$mean += $$values[$i][$j];
	}
	$mean /= $nVal;

	# Compute the median value
	for( $j = 0 ; $j < $nVal ; $j++ )
	{
		$list[$j] = $$values[$i][$j];
	}

    #$median = median($nVal, \@list );
	@slist = sort{ $a <=> $b } @list;
    $median = $slist[$nVal/2];

	# Compute asymmetric error bars
	# Compute the Freedman-Diaconis bin size
	$binsize = $slist[$nVal*0.75] - $slist[$nVal*0.25];
	$binsize *= 2.*$nVal**-0.3333;
	$binsize = 1 if( $binsize == 0);
	$nbin = ($slist[$#slist]-$slist[0])/$binsize;
	$nbin = 1 if( $nbin == 0);
	($x, $histo) = histo($nVal, $nbin, \@list);

	@shisto = sort{ $a <=> $b } @{$histo}; $hmax = $shisto[$#shisto];
	# find the index of the largest bin
	$hmaxid = 0;
	$hmaxid++ while( $hmax != $histo->[$hmaxid] );
	$mode = $x->[$hmaxid];
	# find the mode index in <@list>
	@slist = sort{ $a <=> $b } @list;
	$modeid = 0;
	$modeid++ while( $slist[$modeid] < $mode && $modeid < $nVal );
	$modeid--;
	# find 68% to the left of <modeid>
	$eminid = $modeid - $modeid*0.68;
	#$eminid-- while( $modeid - $eminid < $nVal*0.68 && $eminid > 0 );
	$emin = $mode - $slist[$eminid+1];
	# find 68% to the right of <modeid>
	$emaxid = $modeid + ($nVal - $modeid)*0.68;
	#$emaxid++ while( $emaxid - $modeid < $nVal*0.68 && $emaxid < $nVal );
	$emax = $slist[$emaxid-1] - $mode;
	
	# Compute the stddev value (bias corrected variance)
	$stddev = 0;
	for( $j = 0 ; $j < $nVal; $j++ )
	{
		$stddev += ($$values[$i][$j] - $mean)*($$values[$i][$j] - $mean);
	}
	$stddev /= $nVal - 1;
	$stddev = sqrt( $stddev );

	# Do not print Evidence line
	#next if( $param[$i] =~ "Evidence" );
	print color 'red' if( $hmaxid==0 || $hmaxid==int($nbin-1) );
	if( $median < 1e4 )
	{
		$line = "%s (%.0f): %.4f %.4f %.4f +-%.4f +%.4f -%.4f (68\%)";
	}
	else
	{
		$line = "%s (%.0f): %.4e %.4e %.4e +-%.4e +%.4e -%.4e (68\%)";
	}
	printf $line, $param[$i], $nbin, $median, $best, $mode, $stddev, $emax, $emin;
	print color 'reset';
	print "\n";
}

