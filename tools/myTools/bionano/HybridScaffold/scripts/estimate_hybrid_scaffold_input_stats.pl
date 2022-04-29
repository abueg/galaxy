# $Id: estimate_hybrid_scaffold_input_stats.pl 7410 2018-02-23 23:00:29Z apang $

#!/usr/bin/perl

use strict;
use warnings;

my ($fasta, $cmap1, $cmap2, $out) = ("", "", "", "");
use Getopt::Long;

GetOptions	(
	"fasta=s"	=>	\$fasta,
	"cmap1=s"	=>	\$cmap1,
	"cmap2=s"	=>	\$cmap2,
	"o=s"		=>	\$out
); 

die "ERROR: Please enter in a FASTA file name\n" if ($fasta =~ /^$/);
die "ERROR: Please enter in a CMAP 1 file name\n" if ($cmap1 =~ /^$/);
die "ERROR: Please enter in a output fil ename\n" if ($out =~ /^$/);

# open a output file name
my $summaryOut = $out;
$summaryOut =~ s/\.\w+$/_summary.txt/;

open(OUT, ">$out") or die "ERROR: cannot write to $out: $!\n";
open(OUT2, ">$summaryOut") or die "ERROR: cannot write to $summaryOut: $!\n";

my $command = "perl estimate_fasta_stats.pl $fasta";
my $statistics = `$command`;
print OUT "$statistics\n";
my $summaryStatistics = '{"fasta":';
$summaryStatistics = $summaryStatistics.(convertToJson($statistics)).",";


$command = "perl estimate_cmap_stats.pl $cmap1";
$statistics = `$command`;
print OUT "$statistics\n";
$summaryStatistics = $summaryStatistics.'"cmap1":'.(convertToJson($statistics)).",";

if ($cmap2 !~ /^$/)	{
	$command = "perl estimate_cmap_stats.pl $cmap2";
	$statistics = `$command`;
	print OUT "$statistics\n";
	$summaryStatistics = $summaryStatistics.'"cmap2":'.(convertToJson($statistics)).",";
} # if cmap2 exists

$summaryStatistics =~ s/,$/}/;
print OUT2 "$summaryStatistics";

close OUT;
close OUT2;

sub convertToJson	{
	my ($string) = @_;

	# this subroutine converts the hybrid scaffold statistics output to json format

	$string =~ s/\s+=\s+/":"/g;
	$string =~ s/\s+$//;	# the very last whitespace (\n)
	$string =~ s/\n/","/g;	# change all \n to comma
	$string =~ s/\s+/_/g;	# change the rest of the whitespace to underscore

	# add curly brackets
	$string = '{"'.$string.'"}';

	return $string;
} # convertToJson
