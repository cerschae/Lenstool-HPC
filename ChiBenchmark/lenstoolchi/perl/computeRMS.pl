#!/usr/bin/perl 
# Compute the RMS between 2 sets of multiple images
# from 2 files whose names are given in arguments
#
# file1 : list of predicted images with position relative in arcsec
# file2 : observed image with position relative in arcsec
#

# subroutine that converts an X and Y value read from DS9 according to 
# a reference declaration
# Input @xc and @yc arrays are modified
sub convertXY 
{
	my $iref, $ra, $dec, $pixel, $pixelx, $pxc, $pyc;
	($pxc, $pyc, $iref, $ra, $dec, $pixel, $pixelx) = @_;

	# Convert the coordinates to relative coordinates
	# if $iref == 0, nothing has to be done as it stays in WCS
	if( $iref == 1 or $iref == 3 )
	{
		for( $i=0; $i<= $#$pxc; $i++)
		{
			@$pxc[$i] -= $ra + 360 if( abs(@$pxc[$i] - $ra) > 1 );  #correction if ra=0
			@$pxc[$i] = (@$pxc[$i] - $ra)*$pixelx;
			@$pyc[$i] = (@$pyc[$i] - $dec)*$pixel;
		}
	}
	elsif( $iref == 2 )
	{
		for( $i=0; $i<= $#$pxc; $i++)
		{
			@$pxc[$i] -= $ra;
			@$pyc[$i] -= $dec;
		}
	}
} # end of convertXY subroutine

# Read a list of images from DS9, save them in a file
# and modify the two pointers to the ra dec in relative arcsec arrays
# Return the name of the last image
sub readDS9
{
	my $file;
	( $file, $pra, $pdec ) = @_;

	@ds9=`xpaget ds9 regions selected`;
	open( OUT, ">$file");
	foreach $line ( @ds9 )
	{
		if( $line =~ /ellipse\(([\d|\.]+),([\d|\.|\+|\-]+)/ )
		{
			chop;
			@fld = split( '[{|}]',$line);
			push @{$pra},$1;
			push @{$pdec},$2;
			print OUT "$1 $2\n";
		}
	}
	close(OUT);
	return $fld[1];
}

# Read a catalog of images and return 2 pointers of ra dec arrays
sub readFile
{
	my $file;
	($file,$pra,$pdec) = @_;

	open(IN,$file);
	while(<IN>)
	{
		chop;
		@fld = split, ' ';
		push @{$pra}, $fld[0];
		push @{$pdec}, $fld[1];
	}
	close(IN);
}

#----------------------------------------------

# Main program
# ---------------------------------------------
#

# Private declaration to the main prog
my @ra1,@dec1,@ra2,@dec2;

$file1='s1.cat';
$file2='s2.cat';


# Read DS9,  save into $file1 and exit
if( ! -e $file1 )
{
	print "INFO: Save $file1\n";
	readDS9($file1,\@ra1,\@dec1);
	exit 0;
}

# Read the second set of images and convert to relative coords
$name=readDS9($file2,\@ra2,\@dec2);
@ref=(3,$ra2[0],$dec2[0],3600,-3600*cos($dec2[0]/180.*3.1415926));
convertXY(\@ra2,\@dec2,@ref);


# Retreive the first set of images and convert it
readFile($file1,\@ra1,\@dec1) || die "ERROR: $file1 not found\n";
convertXY(\@ra1,\@dec1,@ref);

# Compute the RMS
for( $i=0; $i <= $#ra1; $i++ )
{
	$dx=$ra1[$i]-$ra2;
	$dy=$dec1[$i]-$dec2;
	$avg+=sqrt($dx*$dx+$dy*$dy);
	$err+=$dx*$dx+$dy*$dy;
}

$avg/=$i;
$var=$err/$i - $avg*$avg;
$rms=sqrt($var);

# Print the distance and error bars
print "$name $avg +-$rms\n";

system("rm $file1 $file2");
