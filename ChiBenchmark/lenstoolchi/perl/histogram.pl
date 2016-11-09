# Algorithm that return the histogram of an array $value
#
# Origin : http://student.ulb.ac.be/~dgonze/SCRIPTS/PERL/histogram.pl
#
# Syntax  : 
#	histo(<number of points>, <number of bins>, <ptr to a list of points>)
#
# Call example :
#  ($index, $a) =  histo(100, 20, \@value); 
#  print @{$index}, @{$a};
#
# where @value is a list of data and @a contain the number of points
# in each bin.
#
sub histo {

	my ($npts, $nbin, $list) = @_;
	my ($i,$j)=(0,0);
	my $jmax;
	my @z=0;
	my @index,
	my @res=0;

	# find the min and max values
	@x=sort{$a <=> $b} @{$list} ; $min=$x[0]; $max=$x[$#x];
	# bin size
	$jmax=($max-$min)/$nbin;
	$i=0;
	$cmax=$min+$jmax;
	$index[0] = ($min + $cmax)/2.;

	for ($j=0;$j<$npts;$j++)
	{
		$z[$i]++;
		if($x[$j]>$cmax && $i < $nbin-1)
		{
			$z[$i]--;
			$i ++;
			$z[$i]++;
			$cmax+=$jmax;
			$index[$i]=$cmax - $jmax/2.;
		}
	}
	
	$res[0] = \@index;
	$res[1] = \@z;
	return @res;
}
	
