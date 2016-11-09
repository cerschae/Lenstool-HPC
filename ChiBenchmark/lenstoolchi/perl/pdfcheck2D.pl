#!/usr/bin/perl

use Getopt::Long;

do "$ENV{'LENSTOOL_DIR'}/perl/rescaleBayes.pl";
$\="\n";
		
$file = "bayes.dat";
$scale = 1;
$noplot = 0; # default: display the PDF on screen
$i = 2; $j = 3; # default: start with parameter 2 and 3 

# Analyse command line arguments
if( @ARGV > 0 )
{
	GetOptions( 
		'scale=f' => \$scale,
		'l1=i' => \$i, 
		'l2=i' => \$j, 
		'n' => \$noplot,
		'h' => \$help);


	$file = $ARGV[0] if( $ARGV[0] );
}
# TODO test $i >= 2

# help
if( $help )
{
	print "Syntax : pdfcheck2D.pl [-scale=<kpcArcsec>] [-l1=<column>] [-l2=<column>] [-n] [-h] [<bayes.fits>]";
	exit;
}

if( ! ( -e $file ) )
{
    print "ERROR: $file not found.\n";
    exit;
}
		
# any scaling of the data?
#if( $scale != 1. )
#{
#	rescaleBayes($fnames, $scale);
#	$fnames = "bayes_rescaled.dat";
#}

# Get the xtitle and ytitle for the plots from bayes.dat
$fnames = "bayes.dat";
if( -e "names.dat" )
{
	$fnames = 'names.dat';
}


print "Read names from $fnames\n";
open(in, $fnames) || die "ERROR: $fnames not found\n";
@lines= grep { /^#/ } <in>;
chop @lines;
shift @lines; # remove #Nsamples
shift @lines; # remove #chi2
close(in);

# Start the loops
$rep="";
$i-=2;
$j-=2;
while( $rep ne "q" )
{
	@fld1= split /\(|\)/, $lines[$i];
	@fld2= split /\(|\)/, $lines[$j];
	$fld1[0] =~ s/ $//; $fld1[0] .= "_";
	$fld2[0] =~ s/ $//;
	$name = $fld1[0].$fld2[0];
	$name =~ s/#//g;
	$name =~ s/ ://g;
	$name =~ s/ /_/g;
	print "Process lines : $lines[$i], $lines[$j] --> $name.ps";
	goto LOOP if( -e $name.'.ps' );

	open(out,'>plot.inp');
	printf out "$file\n";	#file
	printf out "%d\n", $i+3;	#line index
	printf out "%s\n", $lines[$i];  #name
	printf out "%d\n", $j+3;	#line index
	printf out "%s\n", $lines[$j];  #name
	printf out "\n";   # no plot title
	printf out "100\n";  # nb of pixels
	printf out "n\n";		#change X range
	printf out "n\n";		#change Y range
	printf out "n\n";  		#automatic FWHM
	printf out "4\n"; 		#default: 4 pixels
	printf out "n\n";		#resmooth
	printf out "n\n"; 		#plot points
	printf out "n\n"; 		#plot ref points
	close(out);
	system('Histogram2D -n < plot.inp > /dev/null');
	system("mv histogram2D.ps $name.ps");

LOOP:
	if( ! $noplot )
	{
		system("gv $name.ps");
		printf "help:? next:[ENTER] previous:b quit:q ?? ";
		read STDIN, $rep, 1;
	} else
	{
		$rep = "\n";
	}

	if( $rep eq "\n" ) 
	{
		if( $j < $#lines )
		{
			$j++;
		}
		else
		{
			if( $i < $#lines-1 )
			{
				$i++;
			} else
			{
				$rep="q";
			}	
			$j = $i+1;
		}
	}

	if( $rep eq "?" )
	{
		read STDIN, $rep, 1;    # To empty the STDIN queue
		print "pdfcheck2D.pl help\n";
		print "----------------\n";
		print "[ENTER] : go forward\n";
		print "q : quit\n";
	}
}

system("rm -f plot.inp");
