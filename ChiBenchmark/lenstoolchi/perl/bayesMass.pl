#!/usr/bin/env perl
# Return a meanmass.dat file which contains the mass profile and
# the error bars
#
# syntax : bayesMass.pl <.par>
#

sub mydfits
{
    my ($filename, $regexp ) = @_;
    my $fh;
    open($fh, $filename) or die "ERROR: File $filename not found\n";
    $str = "";
    while( not $str =~ "END     " and not $str =~ $regexp )
    {
        sysread($fh, $str, 80);
    }
    close($fh);
    return $str;
}

use Getopt::Long;
do "$ENV{'LENSTOOL_DIR'}/perl/covariance.pl";

if( @ARGV < 1 )
{
	print "Syntax : bayesMass.pl [-r -xc <xc> -yc <yc>] <.par>\n";
	exit 1;
}

# Get options
my ($xc, $yc, $ring);
my $status = GetOptions("xc=f" => \$xc, "yc=f" => \$yc, "r" => \$ring);

# Create the tmp/mass??.fits files
#system("bayesMass $ARGV[0]");

# Get the properties of the fits files
#@fmass=<tmp/mass*.fits>; we cannot use globbing because of Unix glob limit
opendir(DIR,"./tmp") || die "Cannot open ./tmp directory : $!\n";
@fmass = grep { /fits$/ } readdir(DIR);
@fmass = @fmass[0..1000] if( scalar @fmass > 1000);
closedir DIR;

unless( defined $xc && defined $yc )
{
    $imsize=mydfits("tmp/$fmass[0]", "NAXIS1"); 
    $imsize=~s/.+=([\s|\w]+)/$1/;
    $imsize=$imsize/2+0.5;
    $xc = $imsize;
    $yc = $imsize;
}
$pixsize=mydfits("tmp/$fmass[0]", "CDELT2"); 
$pixsize=~s/.+=([\s|\w]+)/$1/;
$pixsize*=3600;
if( $pixsize == 0 )
{
    #$pixsize = 4.6875;
    $pixsize = 8. / 1.28; 
    print "WARNING: pixel scale not found in FITS header. Default $pixsize\"/pix\n";
}
open(IN, "para.out") || die "ERROR: para.out file not found\n";
foreach ( <IN> )
{
	$scale=$1 if( $_ =~ /Conversion.*==([\s|\d|\.]+)/ );
}
close(IN);

# Check intmass or ringmass
my $intsoft = "intmass";
if ( $ring )
{
    print "Use ringmass to compute radial density profile\n";
    $intsoft = "ringmass";
}

# Loop over the tmp/mass??.fits files
$nfiles = 0;
my @vm;   # vector of all the mass from all the files
foreach $f ( @fmass )
{
	$f="tmp/".$f;
	printf "Process file %d/%d : %s  dm/m:%.2f\r",$nfiles+1,$#fmass+1,$f,$erel;
	$f=~s/\.fits$//;

	if( ! -e "$f.dat" )
	{	
		# Prepare the intmass script
		open(IN, ">mass.int" );
	
		print IN << "EOT";
input
	filein $f.fits
	xc $xc
	yc $yc
	eps 0.
	theta 0.
	pixel $pixsize
	scale $scale
	end
output
	fileout $f.dat
	end
fini
EOT
		close(IN);

		if( system("$intsoft mass 2>1 /dev/null") != 0 )
		{
			print "ERROR: $intsoft produced error $?\n";
			exit 1;
		}
	}

	# Analyse the produced file
	open( IN2, "$f.dat") or die "ERROR: cannot open $f.dat\n";

	$i=0;
	$erel = 0.;  # cumulated relative error
	while( <IN2> )
	{
		($ra[$i],$rkpc[$i],$npix[$i],$m) = split(' ');

        push @vm, $m;
		$avg[$i]+=$m;
		$sdev[$i]+=$m*$m;

		# cumulated stats
		$a = $avg[$i]/($nfiles+1);
		$e = sqrt($sdev[$i]/($nfiles+1) - $a*$a);
		$erel+= $e*$e/$a/$a if( $a != 0);

		$i++;
	}
	close(IN2);
	$ni=$i;
	$nfiles++;
}
print "\n";

# Perform the statistics and write the result in meanmass.dat
#my $mmat = pdl( @vm )->reshape($ni, $nfiles);
#my $cov = covariance( $mmat );

open(OUT, ">meanmass.dat");
print OUT "# xc $xc yc $yc $ARGV[0]\n";
for( $i=0; $i<$ni; $i++ )
{
	$avg[$i]/=$nfiles;
	$sdev[$i] = sqrt($sdev[$i]/$nfiles - $avg[$i]*$avg[$i]);
    # $sdev[$i] = sqrt( $cov->slice(",$i")->sum() );
	print OUT "$ra[$i] $rkpc[$i] $npix[$i] $avg[$i] $sdev[$i]\n";
}
close OUT;
