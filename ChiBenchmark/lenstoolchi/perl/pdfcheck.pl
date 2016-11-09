#!/usr/bin/perl

use Getopt::Long;

do "$ENV{'LENSTOOL_DIR'}/perl/rescaleBayes.pl";
$\="\n";
		
$file = "bayes.dat";
$scale = 1;
$noplot = 0; # default: display the PDF on screen
$i = 2;  # default: start with parameter 2 #chi2

# Analyse command line arguments
if( @ARGV > 0 )
{
	GetOptions( 
		'scale=f' => \$scale,
		'l=i' => \$i, 
		'n' => \$noplot);

	$file = $ARGV[0];
}
		
# any scaling of the data?
if( $scale != 1. )
{
	rescaleBayes($file, $scale);
	$file = "bayes_rescaled.dat";
}

# Read the number of bins from the bayesResults.pl STDOUT
open(in, 'bayesResults.pl $file |');
@lines=<in>;
close(in);
shift @lines; # remove Read 20 columns 10010 lines
shift @lines; # remove Param(bin): median ... line[0] becomes #chi2

# Optionally, read the prior.dat 
if( -e 'prior.dat' )
{
	open(P, 'prior.dat');
	@prior = <P>;
	close(P);
}

$rep="";
$i-=2;
while( $rep ne "q" )
{
	chop $lines[$i];
	@fld=split /\(|\)/, $lines[$i];

	# Remove the leading color codes in the name
	$fld[0]=substr($fld[0],index($fld[0],"#"));

	open(out,'>plot.inp');
	printf out "$file\n";	#file
	printf out "%d\n", $i+2;	#line index
	if( $#fld > 5 )
	{	# parameter with unit
		printf out "%s)\n", join("(", @fld[0..$#fld-5]);	#name
	}
	else
	{	# parameter without unit
		printf out "%s\n", join("(", @fld[0..$#fld-4]);        #name
	}
	printf out "n\n";		#change range
	printf out "%d\n", $fld[$#fld-3];	#nbin
	printf out "n\n";		#rebin
	if( @prior )
	{
		printf out "y\n";		#overlay std dist
		@fldp =  split ' ', $prior[$i];
		printf out " $fldp[1]\n";
		printf out " $fldp[2] $fldp[3]\n";
		printf out "n\n";		# Try again?
	} else
	{
		printf out "n\n";		#overlay std dist
	}
	printf out "y\n";		#plot PS
	printf out "n\n"; 		#plot points
	printf out "n\n"; 		#plot ref points
	printf out "n\n"; 		#smooth
	close(out);
	system('Histogram -n < plot.inp > /dev/null');
	chop $fld[0];
	$fld[0] =~ s/#//;
	$fld[0] =~ s/ ://;
	$fld[0] =~ s/ /_/g;
	#print "mv pgplot.ps $fld[0].ps";
	system("mv pgplot.ps $fld[0].ps");

	if( ! $noplot )
	{
		system("gv $fld[0].ps");
		printf "help:? next:[ENTER] previous:b quit:q ?? ";
		read STDIN, $rep, 1;
	} else
	{
		$rep = "\n";
	}

	if( $rep eq "\n" ) 
	{
		if( $i < $#lines )
		{
			$i++;
		} else
		{
			$rep="q";
		}
	}
	if( $rep eq "b" )
	{
		read STDIN, $rep, 1;	# To empty the STDIN queue
		if( $i > 0 )
		{
			$i--;
		}
	}
	if( $rep eq "?" )
	{
		read STDIN, $rep, 1;    # To empty the STDIN queue
		print "pdfcheck.pl help\n";
		print "----------------\n";
		print "[ENTER] : go forward\n";
		print "b : go backward\n";
		print "q : quit\n";
	}
}

system("rm -f plot.inp");
