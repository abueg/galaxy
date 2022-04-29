# $Id: calc_cmap_stats.pl 7288 2018-01-02 23:18:12Z apang $

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
my $cmap_file = $ARGV[0];

my ($cmap_data, $numContig, $contigLength_ref) = readCMap($cmap_file);
my @contigLength=@$contigLength_ref;

my ($min, $max, $mean, $median, $n50value, $total) = getContigStat(@contigLength);
printContigStat($numContig, $min, $max, $mean, $median, $n50value, $total);

=pod

=head1 NAME

calc_cmpa_stats.pl - Calculate stats about CMAPs

=head1 SYNOPSIS

calc_cmpa_stats.pl <CMAP_File>

Examples:

    calc_cmap_stat.pl Ler0_GM0327.cmap

=head1 DESCRIPTION

This script is used to caculate statitical data about a CMAP.

[1] "N contigs = 173"
[1] "Min contig length = 0.1436865 Mb"
[1] "Median contig length = 0.481299 Mb"
[1] "Mean contig length = 0.624446769364162 Mb"
[1] "Contig N50 = 0.75731195 Mb"
[1] "Max contig length = 2.683027 Mb"
[1] "Total contig length = 108.0292911 Mb"


=head1 ARGUMENTS

calc_cmap_stat.pl takes the following arguments:

=over

=item <CMAP_File>

A CMAP File with version 0.1

=back

=cut
