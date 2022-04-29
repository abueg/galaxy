# $Id: calc_scaffolded_seq_fasta_stats.pl 11051 2020-05-20 23:32:23Z apang $

#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

# This adds "${CURRENT_SCRIPT_PATH}/perl5/" direcory at run time to the @INC array
# This script sould sit one level above the additional Perl modules directory.
BEGIN {
        my $script_path = abs_path(dirname($0));
        my $module_path2 = abs_path($script_path . "/perl5");
        unshift @INC, $module_path2;
        my $lib4;
        if ($] >= 5.010000 && $] <= 5.011000) {
                $module_path2 = $module_path2."/5.10.1";
                $lib4 = $module_path2."/x86_64-linux-thread-multi";}
        elsif ($] >= 5.014000 && $] <= 5.015000) {
                $module_path2 = $module_path2."/5.10.1";
                $lib4 = $module_path2."/x86_64-linux-thread-multi";}
        elsif ($] >= 5.016000 && $] <= 5.017000) {
                $module_path2 = $module_path2."/5.10.1";
                $lib4 = $module_path2."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.018000 && $] <= 5.019000) {
		$module_path2 = $module_path2."/5.18.2";
		$lib4 = $module_path2."/x86_64-linux-thread-multi";}
	else {
		print "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X and 5.16.X and 5.18.X\n";
		die "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X and 5.16.X and 5.18.X\n"; 
                exit 1; }
        unshift @INC, $module_path2;
        unshift @INC, $lib4;
        #print "$0 library paths:\n\t"; print join("\n\t", @INC); print "\n";
}

use BNG::Utility;

if (scalar(@ARGV) != 2)	{
	dieLog("ERROR: Please enter in the original (or post-cut original) sequence FASTA file, and the FASTA file containing sequences not used in the hybrid scaffold (as generated by agp-fasta export)\n");
	exit 1;
} # if scalar
my ($oriSeqFile, $notScaffoldedSeqFile) = @ARGV;

# read in the original FASTA file
my $oriSeqRef = getSeq($oriSeqFile, 1);

# determine which of the original sequences were not scaffolded
my $notScaffoldedSeq = getSeq($notScaffoldedSeqFile, 0); 

# calculate the statistics of the remaining, thus scaffolded sequences
my @scaffoldedSeqLength = ();
foreach my $oriSeqName (keys %$oriSeqRef)	{
	my $seqName2 = $oriSeqName."_obj";	# the not scaffolded fasta would have the _obj suffix
	if (exists $notScaffoldedSeq->{$seqName2})	{
		# that sequence was not scaffolded, record
		$oriSeqRef->{$oriSeqName}{scaffolded} = 0;
	} # if exists
} # foreach oriSeqName
foreach my $oriSeqName (keys %$oriSeqRef)	{
	push(@scaffoldedSeqLength, $oriSeqRef->{$oriSeqName}{seqLength}) if ($oriSeqRef->{$oriSeqName}{scaffolded} == 1);
} # foreach oriSeqName
my $numSeq = scalar(@scaffoldedSeqLength);

my ($min, $max, $mean, $median, $n50value, $total) = getContigStat(@scaffoldedSeqLength);
printContigStat($numSeq, $min, $max, $mean, $median, $n50value, $total);



sub getSeq	{
	my ($file, $scaffoldedFlag) = @_;
	my %oriSeq = ();
	my ($curHeader, $curSeq) = ("", "");
	open(IN, $file) or dieLog "ERROR: Cannot open original sequence file $file: $!\n";
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		if ($line =~ /^>/)	{
			# header line
			if ($curHeader !~ /^$/ && $curSeq ne "")	{
				# record
				my $theSeqLength = length($curSeq);
				$oriSeq{$curHeader} = {seqLength => $theSeqLength, scaffolded => $scaffoldedFlag}; 
			} # if curHeader

			if ($line =~ /^>(\S+)\s*/)	{
				# read up to the first whitespace
				$curHeader = $1;
				$curSeq = "";

				# replaces all characters that are not underscore/alphabet/numerals into underscore
				$curHeader =~ s/\W/_/g;
			} # if line
		} else	{
			# sequence line
			if ($curHeader !~ /^$/)	{
				$curSeq .= $line;
			} # if curHeader
		} # if line
	} # while line
	# record the last line
	if ($curHeader !~ /^$/ && $curSeq ne "")	{
		# record
		my $theSeqLength = length($curSeq);
		$oriSeq{$curHeader} = {seqLength => $theSeqLength, scaffolded => $scaffoldedFlag};
	} # if curHeader
	
	close IN;
	
	return \%oriSeq;
} # getSeq