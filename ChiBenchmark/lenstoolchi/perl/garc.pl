#!/usr/bin/perl
# Transform the ellipses in ds9 to a lenstool catalog. If no
# argument is given, write the lenstool catalog on STDOUT.
#
# If the argument is STDIN, read the ds9 region file from
# stdin.
#
# If an argument is given, the ellipses are written in the
# file and if the file already contains ellipses with the
# same identifiers, their position, size and orientations are
# updated.  Their redshifts and magnitudes are not modified.

# The ellipses in ds9 must have different identifiers.

$[ = 0;			# set array base to 1
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

# subroutine that return $ra and $dec according to $iref
sub getRADEC
{
	local $iref, $ra, $dec, $pixel, $pixelx;
	($iref, $ra, $dec, $pixel, $pixelx) = @_;

	#If we have the reference coordinates in sexagesimal form
	#convert them to degrees
	if ($iref == 1)
	{
		($hh,$mm,$ss)=split(':',$ra);
		($dd,$nn,$tt)=split(':',$dec);
		$sign=1;
		if (substr($dec,0,1) eq '-')
	        {
		        $sign=-1;
        		$dd=abs($dd);
		}

		$ra=($hh+$mm/60+$ss/3600)*15;
		$dec=$sign*($dd+$nn/60+$tt/3600);
		$pixel=3600;
		$pixelx=-3600*cos($dec/180*3.1415926);
		$ds9type = 'fk5';
	}
	elsif( $iref == 2 )
	{
		$ds9type = 'image';
	}
	else #( $iref == 3 or $iref == 0 ) default case
	{
		$pixel=3600;
		$pixelx=-3600*cos($dec/180*3.1415926);
		$ds9type = 'fk5';
	}
} # end of getRADEC subroutine

# subroutine that converts an X and Y value read from DS9 according to 
# a reference declaration
sub convertXY 
{
	local $iref, $ra, $dec, $pixel, $pixelx, $xc, $yc, @ref;
	($xc, $yc, $iref, $ra, $dec, $pixel, $pixelx) = @_;
	# Convert the coordinates to relative coordinates
	# if $iref == 0, nothing has to be done as it stays in WCS
	if( $iref == 1 or $iref == 3 )
	{
		$xc -= $ra + 360 if( abs($xc - $ra) > 1 );  #correction if ra=0
		$xc = ($xc - $ra)*$pixelx;
		$yc = ($yc - $dec)*$pixel;
	}
	elsif( $iref == 2 )
	{
		$xc -= $ra;
		$yc -= $dec;
	}
	return ($xc, $yc);
} # end of convertXY subroutine

#----------------------------------------------

#
# Beginning of the main script
#

# Get the right DS9 id
if( $ARGV[0] ne "STDIN" )
{
	if( ! $ENV{DS9} )
	{
		@xpa = `xpaget xpans`;
		@xpa = split(/ /,$xpa[0]);
		$title=$xpa[1];
	}
	else
	{
		$title=$ENV{DS9};
	}

	# Read the DS9 regions output
	system("xpaset -p $title regions system wcs");
	system("xpaget $title regions > e.reg");
	open $ds9 ,"e.reg";
}
else
{
	$ds9 = *STDIN;
	shift @ARGV;
}

#Read DS9 file and fill the %arcs hash map
while (<$ds9>)
{
	chop;
	if( $_ =~ /^#/ )
	{} elsif ( $_=~ /ellipse\((.*)\).*text={([\w|\d|\.]+)}/ )
	{
		if ( $2 ne "" )
		{
			@fld=split(/,/,$1);
			unshift(@fld,$2);
			$fld[6]=0.;	#redshift
			$fld[7]=0.;	#magnitude
			$arcs{$2}=[ @fld ]; # the key is the arc id
		}
	}
}

printf STDERR "%d regions read\n",scalar keys %arcs;
exit 1 if( scalar keys %arcs == 0 );

#if something has to be changed in ARGV[0]
#if ds9 regions have been read
if (scalar keys %arcs) 
{	
	#Default : No #REFERENCE keyword
	$iref = -1;

	#Read $ARGV[0] 
	if( @ARGV > 0 ) {
		open(in,"$ARGV[0]");
		while($line=<in>)
		{
			#Check if there is a reference keyword
			if( $line =~ /#REFERENCE/i )
			{
				($null, $iref, $ra, $dec) = split ' ', $line;
				getRADEC( $iref, $ra, $dec, $pixel, $pixelx );
			}
			@fld=split ' ',$line;
			# Lookup in the %arcs hash map for the same arc id
			if(exists($arcs{$fld[0]}))
			{
				#Update the old values with the new ds9 values
				@fldDs9=@{ $arcs{$fld[0]} };

				# update ra, dec, a, b and theta
				for($i=1;$i<6;$i++) 
				{
					$fld[$i]=$fldDs9[$i];
				}

				# save the new array
				$arcs{$fld[0]}=[ @fld ];
			} else
			{
				#Consider those lines as header
				push @header, $line;
			}
		}
		close(in);
	}
	
	#Switch the writting between ARGV[0] or stdout
	if( @ARGV > 0 ) {
		open($out,">$ARGV[0]");
	} else {
		$out = *STDOUT;
	}
	
	# Print the #REFERENCE keyword if it doesn't exist
	printf $out "#REFERENCE 0\n" if( $iref == -1 );
	
	# Print the header
	for( @header )
	{
		printf $out $_;
	}

	# Print the modified arcs values
	foreach $key ( sort keys %arcs ) 
	{
		@fld=@ { $arcs{$key} };	#retreive the arc array with its key

		# convert the center coordinates
		@ref = ($iref, $ra, $dec, $pixel, $pixelx);
		($fld[1], $fld[2]) = convertXY( $fld[1], $fld[2], @ref );

#		$fld[3]=substr($fld[3],0,-1);
#		$fld[4]=substr($fld[4],0,-1);
		printf $out "%s %.7f %.7f %.5f %.5f %6.1f %.3f %.4f\n", @fld;
	}
	close($out) if $out ne *STDOUT;
}
exit 0;
