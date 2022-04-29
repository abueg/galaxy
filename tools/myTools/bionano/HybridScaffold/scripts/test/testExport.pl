use strict;
require "ExportAGP.pl";

print("abc\n");
my @a = [1,2,3,4,4,5];
my $cut_file = "/home/users3/apang/hello/pearl/retry_052020/output/assignAlignType/cut_conflicts/auto_cut_NGS_coord_translation.txt";
my $sorted_xmap_file = "/home/users3/apang/hello/pearl/retry_052020/output/agp_fasta/EXP_REFINEFINAL1_bppAdjust_cmap_PEARL_HiFicontigs_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap_sorted.xmap";
my $cut_map = readCutCoordFile($cut_file);
my $entry_count = (keys %{$cut_map});
$entry_count == 201 or die "entry count in cut coord file not correct";
print "Number of entries: ".$entry_count."\n";

my $maxId = getMaxKey($cut_map);
print "MaxId: ".$maxId."\n";

my $xmap = readXMap($sorted_xmap_file);
$xmap = computeGapLengths($xmap);
print "Total alignment parsed ".$xmap->{totalHits}."\n";

(my $adjusted_xmap, my $adjust_cut_map) = trimOverlapAlignment($xmap, $cut_map);
#my $adjust_cut_map = $cut_map;

my $new_cut_coordinate_file = '/home/users/jwang/workspace/codebase/HybridScaffoldThree_trunk/trunk/scripts/tmp/cut_coord_modified.txt';
printCutCoordinateFile($adjust_cut_map, $new_cut_coordinate_file);

open(my $fh, '>/home/users/jwang/workspace/codebase/HybridScaffoldThree_trunk/trunk/scripts/tmp/align_modified.xmap');
my $header = 1;
for(my $i = 0; $i < $xmap->{totalHits}; $i++){
	my @keys = keys%{$xmap->{hits}};
	if($header){
		print  $fh "#";
		print $fh (join "\t", @keys);
		$header=0;
		print $fh "\n";
	}
	foreach my $key(@keys){
		print $fh $xmap->{hits}->{$key}[$i];
		print $fh "\t";
	}
	print $fh "\n";
	#print join "\t", (values %{$entry});
	#print "\n";
}