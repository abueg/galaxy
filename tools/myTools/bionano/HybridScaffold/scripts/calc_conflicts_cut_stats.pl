# $Id: calc_conflicts_cut_stats.pl 7237 2017-12-15 21:03:32Z apang $

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
if (scalar(@ARGV) != 1 || $ARGV[0] !~ /\w+/)	{
	dieLog("ERROR: please input the conflicts cut status file name (e.g. conflicts_cut_status.txt)\n");
	exit 1;
} # if scalar

my ($conflictCutStatusFile) =  @ARGV;

my ($numCutOnSeq, $numCutOnMap) = (0, 0);
my ($numSeqCut, $numMapCut) = (0, 0);

my %nrSeqCut = ();
my %nrMapCut = ();

open(IN, $conflictCutStatusFile) or dieLog ("ERROR: Unable to read in conflict cut status file $conflictCutStatusFile: $!\n");
while (my $line = <IN>)	{
	chomp $line;
	$line =~ s/\r+//g;
	# white space
	$line =~ s/^\s+|\s+$//;
	
	next if ($line !~ /^\d+/);
	
	my @content = split(/\t/, $line);
	my ($refId, $refLeftBkptCutStatus, $refRightBkptCutStatus) = ($content[2], $content[6], $content[7]);
	my ($qryId, $qryLeftBkptCutStatus, $qryRightBkptCutStatus) = ($content[10], $content[14], $content[15]);
	
	if ($refLeftBkptCutStatus =~ /cut/i && $refId !~ /^-1$/)	{
		$numCutOnSeq += 1;
		$nrSeqCut{$refId} = 1;
	} 
	if ($refRightBkptCutStatus =~ /cut/i && $refId !~ /^-1$/)	{
		$numCutOnSeq += 1;
		$nrSeqCut{$refId} = 1;
	} 
	if ($qryLeftBkptCutStatus =~ /cut/i && $qryId !~ /^-1$/)	{
		$numCutOnMap += 1;
		$nrMapCut{$qryId} = 1;
	} 
	if ($qryRightBkptCutStatus =~ /cut/i && $qryId !~ /^-1$/)	{
		$numCutOnMap += 1;
		$nrMapCut{$qryId} = 1;
	} 
} # while line
close IN;

print "Number of conflict cuts made to Bionano maps: $numCutOnMap\n";
print "Number of conflict cuts made to NGS sequences: $numCutOnSeq\n";
print "Number of Bionano maps to be cut: ".(scalar(keys %nrMapCut))."\n";
print "Number of NGS sequences to be cut: ".(scalar(keys %nrSeqCut))."\n";
