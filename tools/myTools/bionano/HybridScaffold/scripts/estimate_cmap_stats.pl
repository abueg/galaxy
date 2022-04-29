# $Id: estimate_cmap_stats.pl 7466 2018-03-12 17:13:39Z apang $

#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV) != 1)	{
	die "ERROR: please input the CMAP file name (e.g. hello.cmap)\n";
} # if scalar
my ($cmapFile) = @ARGV;

my ($contigLengthRef) = readCMap($cmapFile);
my ($min, $max, $mean, $median, $n50value, $total) = getContigStat($contigLengthRef);
printContigStat($cmapFile, scalar(@$contigLengthRef), $min, $max, $mean, $median, $n50value, $total);

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



##
# read a CMAP file
##
sub readCMap	{
	my ($cmap_file) = @_;
	# read cmap file (first time to see if there is the keyword "ChimQuality")
	open(IN, "$cmap_file") or die "ERROR: Unable to read in file $cmap_file: $!\n";
	my $foundChimQuality = 0;
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/)	{
			# header line
			if ($line =~ /#h\s+/){
				$foundChimQuality = 1 if ($line =~ /ChimQuality/);
				last;	# done
			}
		} # if line
	} # while line
	close IN;
	# then read it again, will read in according whether ChimQuality is there or not
	
	my $cmap={};
	my $c_cmapId = 0;
	my $numcontig = 0;
	my @contigLength = ();
	my @header = ();
	my @data_name = ();
	my @data_type;
	my $cmap_version = "";
	open(IN, "$cmap_file") or die "ERROR: Unable to read in file $cmap_file: $!\n";
	my $NUM_C_LIMITED = ($foundChimQuality == 1) ? (10) : (9);	# assumption: column 10 has ChimQuality, if it is present
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/)	{
			# header line
			if ($line =~ /^#h\s+/)	{
				# data name:
				my $headTag = "";	# the #h
				($headTag, @data_name) = split(/\s+/, $line);
				splice @data_name, $NUM_C_LIMITED;
				push(@header, join("\t", "$headTag $data_name[0]", @data_name[1..$#data_name]));
			} elsif ($line =~ /^#f\s+/)	{
				# data type:
				my $headTag = "";       # the #f
				($headTag, @data_type) = split(/\s+/, $line);
				splice @data_type, $NUM_C_LIMITED;
				push(@header, join("\t", "$headTag $data_type[0]", @data_type[1..$#data_type]));
			} elsif ($line =~ /^#\s+CMAP File Version:/)	{
				my (@tmp) = split(/\s+/, $line);
				$cmap_version = $tmp[-1];
				push(@header, $line);
			} else	{
				push(@header, $line);
			} # if 
			next;
		} # if header line
		
		# data lines
		if (! exists $cmap->{"headers"})	{
			# first data line, populate the relevant hash values
			$cmap->{"headers"} = \@header;
			$cmap->{"dataName"} = \@data_name;
			$cmap->{"dataType"} = \@data_type;
			$cmap->{"MAPFileVersion"} = $cmap_version;
		} # if exists
		
		my $numc = $NUM_C_LIMITED;
		
		# now store information of that data line
		$line =~ s/^\s*//;	$line =~ s/\s+/\t/g;
		my @d_items = split(/\t/, $line);
		if ($d_items[0] != $c_cmapId)	{
			# a new contig:
			$numcontig += 1;
			$c_cmapId=$d_items[0];
			# ContigLength:
			$cmap->{contigs}->{$c_cmapId}->{$data_name[1]} = $d_items[1];
			push(@contigLength, $d_items[1]);
			# NumSites:
			$cmap->{contigs}->{$c_cmapId}->{$data_name[2]} = $d_items[2];
			
			# the rest of the columns
			for (my $i = 3; $i < $numc; $i += 1)	{
				my $ad = [];
				push(@$ad, $d_items[$i]);
				$cmap->{contigs}->{$c_cmapId}->{$data_name[$i]} = $ad;
			} # for i
		} else	{
			# another site of the same contig:
			for (my $i = 3; $i < $numc; $i += 1)	{
				my $ad = $cmap->{contigs}->{$c_cmapId}->{$data_name[$i]};
				push(@$ad, $d_items[$i]);
				$cmap->{contigs}->{$c_cmapId}->{$data_name[$i]} = $ad;
			} # for i
		} # if same cmap id
	} # while line
	close IN;
	$cmap->{nContigs} = $numcontig;
	return (\@contigLength);
} # readCMap

