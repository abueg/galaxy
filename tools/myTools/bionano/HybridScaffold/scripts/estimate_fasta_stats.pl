# $Id: estimate_fasta_stats.pl 7410 2018-02-23 23:00:29Z apang $

#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV) != 1 )	{
	die "ERROR: please input the FASTA file name (e.g. fasta.fa)\n"; 
} # if scalar
my ($fastaFile) = @ARGV;


my $seqLengthRef = readFasta($fastaFile);

my ($min, $max, $mean, $median, $n50value, $total) = getContigStat($seqLengthRef);

printContigStat($fastaFile, scalar(@$seqLengthRef), $min, $max, $mean, $median, $n50value, $total);


sub readFasta	{
	my ($file) = @_;
	
	my @seqLength = ();
	my ($curHeader, $curSeqLength) = ("", 0);
	open(IN, $file) or die "ERROR: cannot read in $file: $!\n";
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;
		if ($line =~ /^>/)	{
			# header line
			if ($curHeader !~ /^$/ && $curSeqLength != 0)	{
				# record
				push (@seqLength, $curSeqLength);
			} # if curHeader
			
			if ($line =~ /^>(\S+)\s*/)	{
				# read up to the first whitespace
				$curHeader = $1;
				$curSeqLength = 0;

				# replaces all characters that are not underscore/alphabet/numerals into underscore
				$curHeader =~ s/\W/_/g;
			} # if line			
		} else	{
			# data line
			if ($curHeader !~ /^$/)	{
				$curSeqLength += length($line);
			} # if curHeader
		} # if line
	} # while line
	
	# last entry
	if ($curHeader !~ /^$/ && $curSeqLength != 0)	{
		# record
		push(@seqLength, $curSeqLength);
	} # if curHeader
	
	close IN;
	
	return \@seqLength;
} # readFasta

sub getContigStat	{
	# give a contig sequence array reference, calculate its statistics
	my ($lenRef) = @_;

	my $numSeq = scalar(@$lenRef);
	
	if ($numSeq == 0)	{
		return ("NA", "NA", "NA", "NA", "NA", "NA");
	} # if scalar

	my ($min, $max, $mean, $median, $n50Value, $total) = (-1, -1, -1, -1, -1, -1);	
	@$lenRef = sort {$b <=> $a} @$lenRef;
	my $r = $numSeq % 2;
	my $m = ($numSeq - $r) / 2;
	if ($r != 0)	{
		# odd
		$median = $lenRef->[$m];
	} else	{
		# even
		$median = ($lenRef->[$m - 1] + $lenRef->[$m]) / 2;	
	} # if r
	$min = $lenRef->[-1];
	$max = $lenRef->[0];
	$total = 0.0;
	for (my $i = 0; $i < $numSeq; $i += 1)	{
		$total += $lenRef->[$i];
	} # for i
	$mean = $total / $numSeq;
	
	my $n50 = 0.0;
	$n50Value = 0.0;
	if ($numSeq == 1)	{
		$n50Value = $max;
		$n50 = $max;
	} else	{
		for (my $i = 0; $i < $numSeq; $i += 1)	{
			$n50 += $lenRef->[$i];
			if ($n50 >= $total / 2.0)	{
				$n50Value = $lenRef->[$i];
				last;
			} # if exceed or equal to half the total length
		} # for i
	} # if scalar
	return ($min, $max, $mean, $median, $n50Value, $total);
} # getContigStat

sub printContigStat	{
	my ($inFile, $numContig, $min, $max, $mean, $median, $n50value, $total) = @_;
	print "Name = $inFile\n";
	print "Count  = $numContig\n";
	if ($numContig == 0) {
		print "Min length (Mbp) = N/A \n";
		print "Median length (Mbp) = N/A \n";
		print "Mean length (Mbp) = N/A \n";
		print "N50 length (Mbp) = N/A \n";
		print "Max length (Mbp) = N/A \n";
		print "Total length (Mbp) = N/A \n";
	} else {
		print "Min length (Mbp) = " . (sprintf("%.3f", $min/1000000)) . "\n";
		print "Median length (Mbp) = " . (sprintf("%.3f", $median/1000000)) . "\n";
		print "Mean length (Mbp) = " . (sprintf("%.3f", $mean/1000000)) . "\n";
		print "N50 length (Mbp) = " . (sprintf("%.3f", $n50value/1000000)) . "\n";
		print "Max length (Mbp) = " . (sprintf("%.3f", $max/1000000)) . "\n";
		print "Total length (Mbp) = " . (sprintf("%.3f", $total/1000000)) . "\n";
	} # if numContig
} # printContigStat
