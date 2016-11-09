#!/usr/bin/perl
# Compute the sum of the relative errors of all or part of the parameters
# in function of the number of samples
#
# Syntax : bayesError.pl [-r <i-j>] 
#
#
use Getopt::Std;

getopt('r');
@cols=split('-', $opt_r) if( $opt_r );

open(IN, "bayes.dat");
open(OUT,">bayesError.dat");

$samp=0;
LINE: while(<IN>)
{
	next LINE if( /^#/ );

	print "Process sample $samp\r";
	@fld=split(' ');
	@cols = (3,$#fld+1) if( ! @cols );
	$samp++;
	$sum=0.;
	# for all the parameters on this line
	for( $i = $cols[0]-1; $i < $cols[1]; $i++ )
	{
		$avg[$i]+=$fld[$i];
		$err[$i]+=$fld[$i]*$fld[$i];
		
		# mean and stddev for the $i parameter 
		$mean=$avg[$i]/$samp;
		$var=$err[$i]/$samp - $mean*$mean;

		#$stddev=sqrt($var);
		#printf OUT "%f ",$stddev/$mean;
		$sum+=$var/$mean/$mean;
	}
	printf OUT "%f\n",sqrt($sum);

}
