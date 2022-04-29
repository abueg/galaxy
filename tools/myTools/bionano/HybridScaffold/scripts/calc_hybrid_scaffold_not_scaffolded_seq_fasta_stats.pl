# $Id: calc_hybrid_scaffold_not_scaffolded_seq_fasta_stats.pl 7288 2018-01-02 23:18:12Z apang $

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
	dieLog("ERROR: please input the hybrid scaffold FASTA and the not-scaffolded sequence FASTA files\n"); 
	exit 1;
} # if scalar
my ($hybridScaffoldFastaFile, $notScaffoldedSeqFastaFile) = @ARGV;

my %seq = ();
my $seqRef = \%seq;
my $numSeq = 0;
my @seqLength = ();
my $seqLengthRef = \@seqLength;

($seqRef, $numSeq, $seqLengthRef) = readFasta($hybridScaffoldFastaFile, $seqRef, $numSeq, $seqLengthRef);
($seqRef, $numSeq, $seqLengthRef) = readFasta($notScaffoldedSeqFastaFile, $seqRef, $numSeq, $seqLengthRef);

my @theSeqLengths = @$seqLengthRef;
my ($min, $max, $mean, $median, $n50value, $total) = getContigStat(@seqLength);
printContigStat($numSeq, $min, $max, $mean, $median, $n50value, $total);

sub readFasta	{
	my ($file, $seqRef, $numSeq, $seqLengthRef) = @_;
	
	open(IN, "$file") or die "readFasta: cannot open $file: $!\n";
	my ($curHeader, $curSeq) = ("", "");
	while (my $line = <IN>)	{
		chomp $line;	
		$line =~ s/\r//g;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		if ($line =~ /^>/)	{
			# header line
			if ($curHeader !~ /^$/ && $curSeq ne "")	{
				# record
				$seqRef->{$curHeader} = $curSeq;
				$numSeq += 1;
				my $theSeqLength = length($curSeq);	
				push(@$seqLengthRef, $theSeqLength);
			} #curHeader
			if ($line =~ /^>(\S+)\s*/)	{
				# read up to the first whitespace
				$curHeader = $1;
				$curSeq = "";

				# replaces all characters that are not underscore/alphabet/numerals into underscore
				$curHeader =~ s/\W/_/g;
			} # if line

		} else	{
			if ($curHeader !~ /^$/)	{
				$curSeq .= $line;
			} # if curHeader
		} # if line
	} # while line
	if ($curHeader !~ /^$/ && $curSeq ne "")	{
		# get the last sequence
		$seqRef->{$curHeader} = $curSeq;
		$numSeq += 1;
		my $theSeqLength = length($curSeq);	
		push(@$seqLengthRef, $theSeqLength);
	} # if curHeader	
	close IN;
	return ($seqRef, $numSeq, $seqLengthRef);
} # readFasta
