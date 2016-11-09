# Algorithm that return the median of an array $a.
#
# Algorithm : WIRTH.
#
# Call example :
#  print median(100, \@valuea); 
#
# where @value is a list of data
#
sub kth_smallest
{
	my ( $a, $n, $k) = @_;
	my ($i, $j, $l, $m);
	my ($x, $tmp);
	$l = 0; $m = $n -1;

	while( $l < $m )
	{
		$x = $a->[$k];
		$i = $l ;
		$j = $m ;
		do {
			$i++ while( $a->[$i] < $x );
			$j-- while( $x < $a->[$j] );
			if( $i <= $j )
			{
				$tmp = $a->[$i]; 
				$a->[$i] = $a->[$j];
				$a->[$j] = $tmp;
				$i++; $j--;
			}	
		} while( $i <= $j );
		$l = $i if( $j < $k );
		$m = $j if( $k < $i );
	}
	return $a->[$k];
}

sub median
{
	my ($n, $a) = @_;
    return 0 if( $n < 1 ); 
    if( $n % 2 ) {
        $k = ($n-1)/2;
        return kth_smallest($a, $n, $k);
    } else {
        $k = $n/2;
        return kth_smallest( $a, $n, $k);
    }
}
