#!/usr/bin/perl
# Rescale the bayes.dat file

sub rescaleBayes
{
	my ($file,$scale) = @_;
	my @fld, @col, $line;
	open(IN, $file);
	open(OUT, '>', "bayes_rescaled.dat");
	
	# Analyse which column as to be rescaled and rescale it
	$line = 0;
	while(<IN>)
	{
		chop;
		@fld = split;
		if( $_ =~ /#/ )
		{
			if( $_ =~ / rc/  || $_ =~ / Re/ || 
			    $_ =~ / rvir/ || $_ =~ / rs/ )
			{
				push @col, $line;
				$fld[3] = "(kpc)" if( $fld[3] eq "(arcsec)" );
				$fld[2] = "(kpc)" if( $fld[2] eq "(arcsec)" );
			} elsif( $_ =~ /theta/ )
			{
				push @th, $line;
				$fld[2] = "PA";
			}
		}
		else
		{
			for( $i = 0; $i < @col; $i++ )
			{
				$fld[$col[$i]] *= $scale;
			}
			for( $i = 0; $i < @th; $i++ )
			{
				$fld[$th[$i]] += 90;
			}
		}

		printf OUT "%s\n", join( " ", @fld);
		$line++;
	}
	close(IN);
	close(OUT);
}
