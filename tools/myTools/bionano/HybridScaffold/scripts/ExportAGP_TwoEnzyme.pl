# $Id: ExportAGP.pl 4514 2016-02-01 18:26:31Z apang $
#!/usr/bin/perl -w

########################################################################
# File: ExportAGP.pl                                                   #
# Date: 9/23/2015                                                      #
# Purpose: Convert xmap output from HybridScaffold to agp format       #
#                                                                      #
# Author: Jian Wang                                                    #
# Email : jwang@bionano.com                                            #
# Affiliation: Research Department, BioNano Genomics Inc.              #
## Usage:                                                              #
#  
########################################################################

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
	my $lib2; my $lib4;
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
print "\nInfo: Running the command: $0 @ARGV\n";

use BNG::Utility;
use Getopt::Std;
use Tie::File;

#global varaibles

our $DEBUG = 0;
our $TYPE_GAP=0;
our $TYPE_CONTIG=1;
our $MAX_MEM = 2000000;
our $path_sep="/";


#enzyme cut site sequences
our %enzyme = (
	"BSPQI" => "GCTCTTC",
	"BBVCI" => "CCTCAGC",
	"BSMI"  => "GAATGC",
	"BSRDI" => "GCAATG",
	"BSECI" => "ATCGAT",
	"DLE1" => "CTTAAG",
	"BSSSI" => "CACGAG"
	# You can add more enzymes here ...
	
);

sub initCli;
sub sortXMap;
sub updateNGSFiles;
sub printUnUsedNGS;
sub printCutCoordinateFile;
sub trimOverlapAlignment;
sub readCutCoordFile;
sub printXmapFile;
sub cleanUp;


###################Main entry point of the export script########################

our $INPUTS = initCli() || die "";
our $sorted_xmap_file = sortXMap($INPUTS->{xmap_file}, $INPUTS->{out_dir});
our $xmap = readXMap($sorted_xmap_file);


#if need to trim overlapping ngs contig
print "checking ".(defined $INPUTS->{trim_overlap_contig})."\n";
if(defined $INPUTS->{trim_overlap_contig}){
	print "Trim option is on, updating alignment and cut coord file to cut out overlapping contig\n";
	my $cut_map;
	if($INPUTS->{cut_coord_file}){
		$cut_map = readCutCoordFile($INPUTS->{cut_coord_file});	
	}else{
		$cut_map = createCutCoordFile($xmap);
	}
	$xmap = computeGapLengths($xmap); 
	(my $adjusted_xmap, my $adjust_cut_map) = trimOverlapAlignment($xmap, $cut_map);
	my $new_cut_coordinate_file = $INPUTS->{out_dir}."/"."auto_cut_NGS_coord_translation_trimmed_overlap.txt";
	printCutCoordinateFile($adjust_cut_map, $new_cut_coordinate_file);
	$xmap = $adjusted_xmap;
	my $new_xmap_file = $INPUTS->{out_dir}."/"."Align_NGS_trimmed_overlap.xmap";
	printXmapFile($xmap, $new_xmap_file);
	$INPUTS->{cut_coord_file} = $new_cut_coordinate_file;
}


if($INPUTS->{cut_coord_file}){
	print "Detect cut-coord file, updating fasta and fasta name map file\n";
	my ($fasta_file, $ngs_namemap_file) = updateNGSFiles($INPUTS->{fasta_file}, 
														 $INPUTS->{ngs_namemap_file}, 
										                 $INPUTS->{cut_coord_file}, $INPUTS->{out_dir}); 
	$INPUTS->{fasta_file} = $fasta_file;
	$INPUTS->{ngs_namemap_file} = $ngs_namemap_file;										                 
}	

print("Processing fasta related files: ".$INPUTS->{fasta_file}."\t".$INPUTS->{ngs_namemap_file});
our $ngsMap = getNGSMap($INPUTS->{ngs_namemap_file});
our $seqMap = readFasta($INPUTS->{fasta_file}, $INPUTS->{out_dir});


our ($hybridCmap, $numcontig, $contigLength) = readCMap($INPUTS->{hybrid_cmap_file});
print("Number of contigs in hybrid $numcontig\t$contigLength");


processAlign($xmap, $INPUTS->{gap_out_file}, 
					$INPUTS->{agp_out_file}, 
					$INPUTS->{begin_end_file},
					$INPUTS->{paddingGapLen},
					$ngsMap);
printFasta($seqMap, $xmap, $INPUTS->{fasta_out_file}, $hybridCmap, $INPUTS->{cut_site}, $INPUTS->{paddingGap});
printFasta($seqMap, $xmap, $INPUTS->{NCBI_fasta_out_file}, $hybridCmap, $INPUTS->{NCBI_cut_site}, $INPUTS->{paddingGap});
printUnUsedNGS($INPUTS->{agp_out_file}, $ngsMap, $INPUTS->{begin_end_file}, $INPUTS->{Unused_fasta_out}, $seqMap);
cleanUp($INPUTS);

###################Utility functions ############################################
#get the file name from a path strings
sub getFileName{
	my $path = $_[0];
	my $sep = $_[1];
	if(not defined $sep){
		$sep = '/';
	}
	return substr($path, rindex($path, $sep)+1);
}

#strip file extension from file path
sub stripExtension{
	my $file = $_[0];
	my $sep = $_[1];
	if(not defined $sep){
		$sep = '.';
	}
	my $endInd = rindex($file, $sep);
	if($endInd < 0){
		$endInd = length($file);
	}
	return substr($file, 0, $endInd);
}

#obtain the basename from a path
sub getFileBasename{
	my $basename = `basename $_[0]`;
	chomp $basename;
	return $basename;
}


sub round{
	return int($_[0] + 0.5);
}

#creates array by repeating some values n times
sub rep{
	my $unit = $_[0];
	my $len = $_[1];
	my @arry=();
	for(my $i = 0; $i < $len; $i++){
		push(@arry, $unit);
	}
	return(\@arry);
}

###############End of Utility functions ###########################################################



#############Function for pre-processing inputs for the exporter ###################################

#Print help usage message
sub Usage{
	my $usage="Usage:                                                       
     ExportAGP.pl <-h> <-i xmap_file> <-o output_file>                   
       -h    : This help message                                        
       -i    : Input xmap file
       -m    : fasta to cmap name mapping file (i.e. generated from fa2cmap.pl)                                          
       -o    : Output directory   
       -s    : Input fasta file from NGS data
       -c    : CMap file for the hybrid contigs
       -e    : Names of the enzymes or the sequence pattern of the cut-site
       -r    : The output file when hybrid-scaffold was ran with -M option (i.e cut ngs contig/bionano map to resolve conflict)
       -p    : Length of gap to be inserted between overlapping NGF contigs in a super-scaffold contig\n";
    print $usage;
}



#sort xmap by refID and increasing refStartPos and decreasing refEndPos
#this simpifly checking if one contigs is completely embedded in another contigs
sub sortXMap{
	my $xmapFile = $_[0];
	my $out_dir = $_[1];
	my $headerCnt = `grep -c '^#' $xmapFile`;
	#print "header count: $headerCnt\n";
	chomp $headerCnt;
	$headerCnt = $headerCnt; 
	my $sorted_xmap = "$out_dir/".getFileBasename($xmapFile)."_sorted.xmap";
	print $sorted_xmap."\n";
	`head -n $headerCnt $xmapFile > $sorted_xmap`;  #export header to sorted xmap 
     $headerCnt = $headerCnt+1; #set to appropriate param for tail
	`tail -n+$headerCnt $xmapFile | sort -k3,3n -k6,6n -k7,7nr >> $sorted_xmap`;
	return($sorted_xmap);
}


#Parse command line options and extract input parameters
#Parameters are stored in the Hash %INPUTS
sub initCli(){
	my $opt_string = "hti:o:s:c:e:m:r:p:";
	my %options=();
	my %INPUTS =();
	if(!getopts("$opt_string", \%options)){
		print("ERROR: Invalid parameter(s)! Must have -i and -o options.\n");
		Usage();
		exit 1;
	}

	if(!defined $options{"i"} || !defined $options{"o"} || !defined $options{"m"}){
		print("ERROR: Missing parameter(s)! Try -h for more information.\n");
		Usage();
		exit 1;
	}
	
	$INPUTS{xmap_file} = $options{"i"};
	$INPUTS{out_dir} = $options{"o"};
	$INPUTS{agp_out_file} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file}))).".agp";
	$INPUTS{begin_end_file} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file})))."_trimHeadTailGap.coord";
	$INPUTS{fasta_out_file} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file}))).".fasta";
	$INPUTS{NCBI_fasta_out_file} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file})))."_NCBI.fasta";
	$INPUTS{Unused_fasta_out} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file})))."_NOT_SCAFFOLDED.fasta"; #fasta file for unused ngs contigs
	$INPUTS{gap_out_file} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file}))).".gap";
	#$auxiliary_file = stripExtension($out_file)."_trimmedcontig.txt";
	print "output agp file: $INPUTS{agp_out_file}\n";
	
	$INPUTS{ngs_namemap_file} = $options{"m"};
	
	if(defined $options{"r"}){
		$INPUTS{cut_coord_file} = $options{"r"};
	}
	
	if(defined $options{"p"}){
		$INPUTS{paddingGapLen} = $options{"p"};
		if($INPUTS{paddingGapLen} < 1){
			warn("Padding gap length must be at least 1, it is now $INPUTS{paddingGapLen}, using default value 13 instead");
			$INPUTS{paddingGapLen} = 13;
		}
	}else{
		$INPUTS{paddingGapLen} = 13;
	}
	print "Padding gap length between overlapping ngs contig: $INPUTS{paddingGapLen}\n";
	my $gap="";
	for(my $c=0; $c < $INPUTS{paddingGapLen}; $c++){
		$gap = $gap."N";
	}
	
	$INPUTS{paddingGap}=$gap;
	

	if(defined $options{"t"}){
		print "Trimmed option on\n";
		$INPUTS{trim_overlap_contig} = 1;
	}

	if(defined $options{"s"} || defined $options{"c"} || defined $options{"e"}){
		if(defined $options{"s"} && defined $options{"c"} && defined $options{"e"}){		
			$INPUTS{fasta_file} = $options{"s"};
			$INPUTS{hybrid_cmap_file} = $options{"c"};
			$INPUTS{cut_site} = $options{"e"};
			$INPUTS{fasta_out} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file}))).".fasta";
			#if we can find the name of the enzyme in table we used the known cut-site sequence
			#otherwise we assume the cut_site is input instead
			#cut_sites is in multi-color format: enzyme1 label1 enzyme2 label2
			my @tokens = split(/\s+/, $INPUTS{cut_site});
			#print "argument values for -e: $cut_site\n";
			#my $all_cap_name = uc $tokens[0];  #only use first enzyme for now 
                        #if($enzyme{$all_cap_name}){
			#	print "Found enzyme $INPUTS{cut_site} with cut site: ".($enzyme{$all_cap_name})."\n";
			#	$INPUTS{cut_site} = $enzyme{$all_cap_name};
			#}else{
			#	print "Cannot find enzyme $INPUTS{cut_site}, using input as cut-site sequence\n";
			#}
			$INPUTS{cut_site} = getEnzymeSeqMap(\@tokens, \%enzyme);
			$INPUTS{NCBI_cut_site} = ['N', 'N', 'N', 'N', 'N'];
		}else{
			print("ERROR: Missing -s or -c or -e options! Try -h for more information.\n");
			print("-s: ".$options{"s"}."\n");
			print("-c: ".$options{"c"}."\n");
			print("-e: ".$options{"e"}."\n");
			Usage();
			return 0;
		}
	}
	if(defined $options{"x"}){
	 
	}
	return \%INPUTS;
}

sub getEnzymeSeqMap{
    #enzymes should be in <enzyme name> <channel #> <enzyme name> <channel number> format
    my @enzymes_name = @{$_[0]};
    print "@enzymes_name\n";
    my $enzymes_table = $_[1];
    my @enzymes_seq = (0) x ($#enzymes_name+1)/2;
    for(my $i = 0; $i <= $#enzymes_name; $i=$i+2){
	my $name = uc $enzymes_name[$i];
	if($enzymes_table->{$name}){
	    print "Found enzyme ".$enzymes_name[$i]." with cut sites: ".($enzymes_table->{$enzymes_name[$i]})."\n"; 
	    $enzymes_seq[$enzymes_name[$i+1]] = $enzymes_table->{$name};
	}else{
	    print "Cannot find enzyme ".$name." using input as cut-site sequence\n";
	    $enzymes_seq[$enzymes_name[$i+1]] = $enzymes_name[$i];
	}
    }
    if($DEBUG){print "Sequence cut sites: @enzymes_seq\n";}
    return \@enzymes_seq;
}

#parse the ngs mapping file to create name mapping between contig id in cmap and xmap and its ngs name and other stats
sub getNGSMap{
	my $name_map_file = $_[0];
	my $cut_coord_file = $_[1];
	my %name_map=();
	
	open(NAMEMAP, $name_map_file) || dieLog "ERROR: cannot open sequence name map file $name_map_file";
	while(my $line = <NAMEMAP>){
		if($line =~ /^#/ || $line !~ /^\d/){
			next;
		}
		chomp $line;
		my @tokens = split(/\t/, $line);
		#we store the ngs stat as 1) name; 2) length; 3) a flag whether it is used scaffolding; 4) the begin coord 
		#of this contig in case it is split by the pipline into multiple smaller contigs
		my $name = (split(/\s+/, $tokens[1]))[0];
		#$name_map{$tokens[0]} = [$tokens[1], $tokens[2], 0, 1];
		$name_map{$tokens[0]} = [$name, $tokens[2], 0, 1];
		
	}
	if(defined $cut_coord_file && length($cut_coord_file) > 0){
		print "cut_coord_file is: $cut_coord_file\n";
		open(CUT_COORD, $cut_coord_file) || dieLog "ERROR: cannot open cut coordinate files $cut_coord_file\n";
		while(my $line = <CUT_COORD>){
			if($line =~ /^oldId/){ #skipping header
				next;
			}
			chomp $line;
			my @tokens = split(/\t/, $line);
			if(exists $name_map{$tokens[0]}){
				my @ngsstat = @{$name_map{$tokens[0]}};
				my $newlength = $tokens[2] - $tokens[1] + 1; 
				$name_map{$tokens[3]} = [$ngsstat[0], $newlength, 0, $tokens[1]+1];
			}else{
				warn "Warning map with Id: $tokens[0] is not recognized to be a valide ID, check $cut_coord_file or $name_map_file\n";
			}
		}
	}
	if($DEBUG){
		print "Total ngs named parsed: ".(keys %name_map)."\n";
		print "printing out ngs map file\n";
		my $count=1;
		foreach my $key(keys %name_map){
			print $key."\t";
			print join "\t", @{$name_map{$key}};
			print "\n";
			$count = $count + 1;
		}
	}	
	return \%name_map;
}


#read and parse ngs coordinate file
#return a table keyed by newId, which is needed
#for trimming of overlapping contigs
sub readCutCoordFile{
	my $cut_coord_file = $_[0];
	my $key_ind = 3; #this control which column to key the table 
	my %name_map=();
	if(defined $cut_coord_file && length($cut_coord_file) > 0){
		print "cut_coord_file is: $cut_coord_file\n";
		open(CUT_COORD, $cut_coord_file) || dieLog "ERROR: cannot open cut coordinate files $cut_coord_file\n";
		while(my $line = <CUT_COORD>){
			if($line =~ /^oldId/){ #skipping header
				next;
			}
			chomp $line;
			my @tokens = split(/\t/, $line);
			my %new_entry = (oldId => $tokens[0],
			                         oldStart => $tokens[1],
									 oldEnd => $tokens[2],
									 newId => $tokens[3],
									 newStart=> $tokens[4],
									 newEnd => $tokens[5]
			                        );
			$name_map{$tokens[$key_ind]} = \%new_entry;
		}
	}
	return \%name_map;
}


#create a dummy cut coordinate files when none is supplied
sub createCutCoordFile{
	my $xmap = $_[0];
	my %name_map=();
	my $i = 0;
	while($i < $xmap->{totalHits}){
		my $qryId = $xmap->{hits}->{QryContigID}[$i];
		my $qryLen = $xmap->{hits}->{QryLen}[$i];
		my %new_entry = (oldId => $qryId,
			                         oldStart => 0,
									 oldEnd => $qryLen - 1,
									 newId => $qryId,
									 newStart=> 0,
									 newEnd => $qryLen - 1
			                        );
		$name_map{$qryId} = \%new_entry;
		$i++;
	}
	return \%name_map;
}

#when hybrid-scaffold cut ngs contigs to resolve conflicts
#we need to generate new fasta files and ngs name map files for 
#the cutted ngs contigs
sub updateNGSFiles{
	my $fasta_file = $_[0];
	my $ngs_namemap_file = $_[1];
	my $cut_coord_file = $_[2];
	my $out_dir = $_[3];
	
	my $cutted_fasta_file = "$out_dir/".getFileBasename($fasta_file).".cutted.fasta";
	my $cutted_ngsNameMap_file = "$out_dir/".getFileBasename($ngs_namemap_file).".cutted.txt";
	
	my  %ngs_map = %{getNGSMap($ngs_namemap_file, $cut_coord_file)};
	my $fasta_accessor = readFasta($fasta_file, $out_dir);
	
	

	open(my $cutted_fasta, ">".$cutted_fasta_file) || dieLog "ERROR: cannot open file for writing: $cutted_fasta_file";
	open(my $cutted_ngs_nameMap, ">".$cutted_ngsNameMap_file) || dieLog "ERROR: cannot open file for writing: $cutted_ngsNameMap_file";

	print $cutted_ngs_nameMap "CompntId\tCompntName\tCompntLength\n";
	$DEBUG=1;
	foreach my $key(keys %ngs_map){
		my @ngscontig = @{$ngs_map{$key}};
		my $fasta_seq = $fasta_accessor->($ngscontig[0]);
		if(length($fasta_seq)- $ngscontig[1] < 2 ){ 
			print $cutted_fasta ">$ngscontig[0]\n";
			print $cutted_fasta $fasta_seq."\n";
			print $cutted_ngs_nameMap $key."\t".$ngscontig[0]."\t".$ngscontig[1]."\n";
			if($DEBUG){
				print "getting seq: $ngscontig[0]\n";
				print "total length: ".(length($fasta_seq))."\n";
				print (join "\t", @ngscontig);
				print "\n";
			} 			
		}else{
			my $new_ngs_name = $ngscontig[0]."_subseq_".$ngscontig[3].":".($ngscontig[3] + $ngscontig[1]-1);
			if($DEBUG){
				print $cutted_fasta ">$new_ngs_name\n";
				print "getting seq: $new_ngs_name\n";
				print "total length: ".(length($fasta_seq))."\n";
				print (join "\t", @ngscontig);
				print "\n";
			}
			print $cutted_fasta substr($fasta_seq, $ngscontig[3] - 1, $ngscontig[1])."\n";
			print $cutted_ngs_nameMap $key."\t".$new_ngs_name."\t".$ngscontig[1]."\n";
		}	
	}
	print "DONE updating NGS files\n";
	close($cutted_fasta);
	close($cutted_ngs_nameMap);
	return ($cutted_fasta_file, $cutted_ngsNameMap_file);
}



#Read in a fasta file, return a function pointer which can be used 
#to access the fasta sequence by NGS name
sub readFasta{
	#open FASTA, $_[0] || die "cannot open fasta file";
	my $fasta_file = $_[0];
	my $out_dir = $_[1];
	
	my ($fastaMap, $tmp_file) = processFastaFile($fasta_file, $out_dir);
	tie my @fasta_arry, 'Tie::File',  $tmp_file,  or dieLog "ERROR: cannot process fasta file";
	my %fasta_map = %{$fastaMap};
	
	my $accesser = sub{
		my $ID = $_[0];
		if(exists $fasta_map{$ID}){
			my $ind = $fasta_map{$ID};
			return $fasta_arry[$ind];
		}else{
			print "Warning: cannot find sequence for ID: $ID\n";
			return "";
		}
	};	
	return $accesser;
}

#This function remove extra white space in a fasta file, so the sequence remains in one line
#It also creates an index table of seq name to the which line in the file
#corresponding to that sequence
sub processFastaFile{
	my $tmp_file = getFileBasename($_[0]);
	my $out_dir = $_[1];
	$tmp_file =  "$out_dir/$tmp_file\.tmp.fa\.tmp";
	open FASTA, $_[0] || dieLog "ERROR: cannot open fasta file";
	open TMP, ">".$tmp_file || dieLog "ERROR: cannot access file system for temporary file";
	my %fastaTable=();
	my $curr="";
	my $currHeader="";
	my $numLines=1;
	while(<FASTA>){
		my $line = $_;
		chomp $line;
		if($line =~ /^>/){
			if(length($curr) > 0){
				$fastaTable{substr($currHeader,1)}=$numLines;
				print TMP $currHeader."\n";
				print TMP $curr."\n";
				$numLines=$numLines+2;
			}else{
				warn("Sequence  ".$currHeader." has zero length...\n")
			}
			my @tokens = split(/\s+/, $line);
			$currHeader = $tokens[0];
			$curr="";
		}else{
			$curr=$curr.$line;
		}
	}
	#taking care of the last seq contig
	$fastaTable{substr($currHeader,1)}=$numLines;	
	print TMP $currHeader."\n";
	print TMP $curr."\n";
	
	close(FASTA);
	close(TMP);
	return (\%fastaTable, $tmp_file)
}




##########################Functions for processing the hybrid-scaffold to NGS contigs alignment (i.e. the xmap file) ##########################

#This function read-in the alignments btw hybrid-scaffold contigs and the NGS contigs and compute gap length
#bewteen each consecutive pair of scaffolded NGS contigs and output this information into a file

sub processAlign{
	my $xmap = $_[0]; 
	my $gapLen_out_file = $_[1];
	my %alignMap = %{getAlignMap($xmap)};
	my $agp_out_file = $_[2];
	my $begin_end_file = $_[3];
	my $paddingGapLen = $_[4];
	my $ngs_map = $_[5];
	$xmap = mapNGSStat($xmap, $ngs_map);
	$xmap = computeGapLengths($xmap);
	printGapFile($xmap, $gapLen_out_file);
	printAGPFile($xmap, $agp_out_file, $begin_end_file, $paddingGapLen, $ngs_map);
}

#This modify the alignment xmap by trimming the overlapping NGS contig alignment 
#the "repeating" sequence in nearby contigs
#note we assume the xmap has already compute gap length
sub trimOverlapAlignment{
	my $xmap = $_[0];
	my %cut_coord_map = %{$_[1]};
	my $maxOldId = getMaxKey(\%cut_coord_map);
	my $counter = 1;
	my $i = 1;
	while($i < $xmap->{totalHits}){
		my $gapLen = $xmap->{hits}->{GapLen}[$i];
		my $adjustedGapLen = $xmap->{hits}->{AdjustedGapLen}[$i];
		my $embedded = $xmap->{hits}->{IsEmbedded}[$i];
		my $qryId = $xmap->{hits}->{QryContigID}[$i];
		#print "Qry: ".$qryId."\t".$embedded."\n";
		#print "stat: ".$gapLen."\t".$adjustedGapLen."\t".$embedded."\n";
		if((not $embedded) && $adjustedGapLen  < 0){
			my $qryId = $xmap->{hits}->{QryContigID}[$i];
			my %cut_coord_entry = %{$cut_coord_map{$qryId}};
			my $newId = $maxOldId + $counter;
			my %new_entry;
			#gap length is negative for overlapping, we only want the magnitude here
			$adjustedGapLen = -1*$adjustedGapLen; 

			#getting info from previous alignment
			my $adjustedPrevRefEnd = $xmap->{hits}->{PrevAdjustedRefEnd}[$i];
			
			#update reference
			#$xmap->{hits}->{RefStartPos}[$i] = $xmap->{hits}->{RefStartPos}[$i] + $adjustedGapLen + 1;
			my $PosAdjustLen = $adjustedPrevRefEnd - $xmap->{hits}->{RefStartPos}[$i] + 1;

			if($PosAdjustLen > 0){
				$xmap->{hits}->{RefStartPos}[$i] = $xmap->{hits}->{RefStartPos}[$i] + $PosAdjustLen;
			}

			#update query and modify coordinate and add entry in cut coordinate file 
			#print "updateing cut coordintate file\n";
			#note below we actually use a trick that actually cut the redundant part from the middle of the contig
			#to make the position all work out the alignment is not realy valid in a strict sense
			if($xmap->{hits}->{Orientation}[$i] eq '+'){	
				$xmap->{hits}->{QryLen}[$i] = $xmap->{hits}->{QryLen}[$i] - $adjustedGapLen;
				
				if($PosAdjustLen < 0){
					$xmap->{hits}->{QryStartPos}[$i] = -1*($PosAdjustLen-1); 
				}else{
					$xmap->{hits}->{QryStartPos}[$i] = 1;
				}
				$xmap->{hits}->{QryEndPos}[$i] = $xmap->{hits}->{QryEndPos}[$i] - $adjustedGapLen; 

				#create new entry in cut file for overlap/duplicate ngs contig
				%new_entry = (oldId => $cut_coord_entry{oldId},
			                  oldStart => $cut_coord_entry{oldStart},
							  oldEnd => $cut_coord_entry{oldStart} + $adjustedGapLen - 1,
							  newId => $newId,
							  newStart=> 0,
							  newEnd => $adjustedGapLen - 1
			                );

				#shift aligned contig in cut coord file
				$cut_coord_entry{oldStart} = $cut_coord_entry{oldStart} + $adjustedGapLen;
				$cut_coord_entry{newEnd} = $cut_coord_entry{newEnd} - $adjustedGapLen;
			}else{
				$xmap->{hits}->{QryLen}[$i] = $xmap->{hits}->{QryLen}[$i] - $adjustedGapLen;
				#$xmap->{hits}->{QryStartPos}[$i] = $xmap->{hits}->{QryStartPos}[$i] - $adjustedGapLen - 1;
				if($PosAdjustLen < 0){
					#$xmap->{hits}->{QryStartPos}[$i] = $xmap->{hits}->{QryLen}[$i] - $PosAdjustLen;
				}else{
					$xmap->{hits}->{QryStartPos}[$i] =  $xmap->{hits}->{QryLen}[$i];
					$xmap->{hits}->{EndPos}[$i] =  1;
				}

				%new_entry = (oldId => $cut_coord_entry{oldId},
			                  oldStart => $cut_coord_entry{oldEnd} - ($adjustedGapLen -1),
							  oldEnd => $cut_coord_entry{oldEnd},
							  newId => $newId,
							  newStart=> 0,
							  newEnd => $adjustedGapLen - 1
			                  );

				$cut_coord_entry{oldEnd} = $cut_coord_entry{oldEnd} - $adjustedGapLen;
				$cut_coord_entry{newEnd} = $cut_coord_entry{newEnd} - $adjustedGapLen;
			}
			$cut_coord_map{$qryId} = \%cut_coord_entry;
			$cut_coord_map{$newId} = \%new_entry;
			$counter++;
		}
		$i++;
	}
	return ($xmap, \%cut_coord_map);
}

#get the maximum key from a hash table
sub getMaxKey{
	my %table = %{$_[0]};
	my @keys = keys %table;
	my @sorted_keys = sort { $a <=> $b } @keys;
	return($sorted_keys[-1])
}

#This function maps the QryContigID in xmap to the ngs stats such as name and length 
#from the original ngs contigs 
sub mapNGSStat{
	my $xmap = $_[0];
	my $ngs_map = $_[1];
	$xmap->{hits}->{NGSName} = rep("", $xmap->{totalHits});
	$xmap->{hits}->{NGSLen} = rep(0, $xmap->{totalHits});
	$xmap->{hits}->{NGSStart} = rep(1, $xmap->{totalHits});
	for(my $i = 0; $i < $xmap->{totalHits}; $i++){
		my $ID = $xmap->{hits}->{QryContigID}[$i];
		if(exists $ngs_map->{$ID}){
			#print "NSG name: ".$ngs_map->{$ID}->[0]."\n"; 
			$xmap->{hits}->{NGSName}[$i] =$ngs_map->{$ID}->[0];
			$xmap->{hits}->{NGSLen}[$i] =$ngs_map->{$ID}->[1];
			$xmap->{hits}->{NGSStart}[$i] = $ngs_map->{$ID}->[3];
		}else{
			print "Warning: cannot find ngs name for QryContig ID: ".$ID.". Using stats from qry in Xmap instead.\n";
			#if we cannot find ngs for some reasons use stat from original xmap
			$xmap->{hits}->{NGSName}[$i] =$ngs_map->{hits}->{QryContigID}[$i];
			$xmap->{hits}->{NGSLen}[$i] =$ngs_map->{hits}->{QryLen}[$i];
		} 
	}
	return $xmap;
}


#This function group alignments that has the same RefContigID
#It return a map of row indices of all alignment has the same RefContigID
sub getAlignMap{
	my $xmap = $_[0];	
	my %alignMap = ();
	print "Total alignments: ".($xmap->{totalHits})."\n";
	for(my $i = 0; $i < $xmap->{totalHits}; $i++){
		my $hybridContigID = $xmap->{hits}->{RefContigID}[$i];
		#print "hybrid: $hybridContigID\n";
		my @alignmentInds=();
		if(exists $alignMap{$hybridContigID}){
			@alignmentInds = @{$alignMap{$hybridContigID}};
		}
		#print "i: $i\n";
		push(@alignmentInds, $i);
		#print "arry: @alignments\n";
		$alignMap{$hybridContigID}=\@alignmentInds;
	}
	return \%alignMap;
}


#check direction
sub isForwardDir{
	return $_[0] eq '+';
}	

#adjust the ref start/end positions of alignments taking into account
#un-aligned sequence and direction of the query contig
sub adjustRefStart{
	my($qryStart, $qryEnd, $refStart, $qryLen, $direction) = @_;
	return $refStart - beginGap($qryStart, $qryEnd, $qryLen, $direction);	
}

sub beginGap{
	my($qryStart, $qryEnd, $qryLen, $direction) = @_;
	
	if(isForwardDir($direction)){
		return ($qryStart - 1);
	}else{
		return ($qryLen - $qryStart);
	}
}

sub adjustRefEnd{
	my($qryStart, $qryEnd, $refEnd, $qryLen, $direction) = @_;
	#print "$qryStart\t$qryEnd\t$refEnd\t$qryLen\t$direction\n";
	return $refEnd + EndGap($qryStart, $qryEnd, $qryLen, $direction);
}

sub EndGap{
	my($qryStart, $qryEnd, $qryLen, $direction) = @_;

	#print "$qryStart\t$qryEnd\t$refEnd\t$qryLen\t$direction\n";
	if(isForwardDir($direction)){
		return ($qryLen - $qryEnd);
	}else{
		return ($qryEnd - 1);
	}
}

#compute gaps for a particular hybrid-scaffolded contig
#given the alignments of all ngs contigs to the hybrid
sub computeGapLengths{
	my $xmap = $_[0];
	#the xmap data table will be used as the main table to store data
	#we added three-more columns to store the gap length information
	$xmap->{hits}->{GapLen} = rep(0, $xmap->{totalHits});
	$xmap->{hits}->{AdjustedGapLen} = rep(0, $xmap->{totalHits}); 
	$xmap->{hits}->{IsEmbedded} = rep(0, $xmap->{totalHits});
	
	my $alignGapLen = 0;
	my $adjustedGapLen = 0;
	my $adjustedRefStart1=0;
	my $adjustedRefEnd1 = 0;
	my $adjustedRefStart2= 0;
	my $adjustedRefEnd2 = 0;
	
	my $currInd = 0;
	my $nextInd = 1;
	
	while($nextInd < $xmap->{totalHits}){
		if($DEBUG){
			print "GapInd $currInd\t$nextInd\n";
		}
		my $refID1 = $xmap->{hits}->{RefContigID}[$currInd];
		my $refID2 = $xmap->{hits}->{RefContigID}[$nextInd];
		if($refID1 ne $refID2){
			$alignGapLen = 0;
			$adjustedGapLen = 0;
			$currInd = $nextInd;
		}else{
			my $qryStartPos1 = $xmap->{hits}->{QryStartPos}[$currInd];
			my $qryStartPos2 = $xmap->{hits}->{QryStartPos}[$nextInd];
			my $qryEndPos1 = $xmap->{hits}->{QryEndPos}[$currInd];
			my $qryEndPos2 = $xmap->{hits}->{QryEndPos}[$nextInd];
		
			my $refStartPos1 = $xmap->{hits}->{RefStartPos}[$currInd];
			my $refStartPos2 = $xmap->{hits}->{RefStartPos}[$nextInd];
			my $refEndPos1 = $xmap->{hits}->{RefEndPos}[$currInd];
			my $refEndPos2 = $xmap->{hits}->{RefEndPos}[$nextInd];

			my $qryLen1 = $xmap->{hits}->{QryLen}[$currInd];
			my $qryLen2 = $xmap->{hits}->{QryLen}[$nextInd];
		
			my $direction1 = $xmap->{hits}->{Orientation}[$currInd];
			my $direction2 = $xmap->{hits}->{Orientation}[$nextInd];
			
			my $qryID1 = $xmap->{hits}->{QryContigID}[$currInd];
			my $qryID2 = $xmap->{hits}->{QryContigID}[$nextInd];
						
			$alignGapLen = round($refStartPos2) - round($refEndPos1) - 1;  #11/12/2015 remove plus one from gap length computation
			$adjustedRefStart1 = round(adjustRefStart($qryStartPos1, $qryEndPos1, $refStartPos1, $qryLen1, $direction1));
			$adjustedRefEnd1 = round(adjustRefEnd($qryStartPos1, $qryEndPos1, $refEndPos1, $qryLen1, $direction1));
			$adjustedRefStart2=round(adjustRefStart($qryStartPos2, $qryEndPos2, $refStartPos2, $qryLen2, $direction2));
			$adjustedRefEnd2=round(adjustRefEnd($qryStartPos2, $qryEndPos2, $refEndPos2, $qryLen2, $direction2));

			if($DEBUG){
				print "$qryStartPos1\t$qryEndPos1\t$qryLen1\t$refStartPos1\t$refEndPos1\t$direction1\t$qryID1\t$refID1\t$adjustedRefStart1\t$adjustedRefEnd1\n";
				print "$qryStartPos2\t$qryEndPos2\t$qryLen2\t$refStartPos2\t$refEndPos2\t$direction2\t$qryID2\t$refID2\t$adjustedRefStart2\t$adjustedRefEnd2\n";
			}
			
			$adjustedGapLen = $adjustedRefStart2 - $adjustedRefEnd1  - 1;
			
			if($DEBUG){print("Gap is ".$alignGapLen."\t".$adjustedGapLen."\n");}
			#detecting contig that is completely covered by the previous contig
			#if this is the case we do not advance the current contig index until the next contig that is not embedded
			if(-1*$adjustedGapLen  > $qryLen2){
				$xmap->{hits}->{IsEmbedded}[$nextInd] = 1;
			}else{	
				$currInd = $nextInd;
			}	
			#7/7/2016 note in this revision we store the gap length information in the contig at the end of the gap 
			#rather than beginning as compared to previous version to properly handle potential cases of  multiple embedded contigs
			$xmap->{hits}->{GapLen}[$nextInd] = $alignGapLen;
			$xmap->{hits}->{AdjustedGapLen}[$nextInd] = $adjustedGapLen;
			$xmap->{hits}->{AdjustedGapBegin}[$nextInd] = $adjustedRefEnd1 + 1;
			$xmap->{hits}->{AdjustedGapEnd}[$nextInd] = $adjustedRefStart2-1;

			#storing ref position about the alignment about the consecitve contigs, this informatioin is needed for trimming overlapping contigs
			$xmap->{hits}->{AdjustedRefStart}[$nextInd] = $adjustedRefStart2;
			$xmap->{hits}->{AdjustedRefEnd}[$nextInd] = $adjustedRefEnd2;
			$xmap->{hits}->{PrevAdjustedRefStart}[$nextInd] = $adjustedRefStart1;
			$xmap->{hits}->{PrevAdjustedRefEnd}[$nextInd] = $adjustedRefEnd1;
			$xmap->{hits}->{PrevRefStart}[$nextInd] = $refStartPos1;
			$xmap->{hits}->{PrevRefEnd}[$nextInd] = $refEndPos1;
		}
		
		$nextInd = $nextInd + 1;
	}
	#if($DEBUG){
		printXmapFile($xmap, $INPUTS->{out_dir}."/Align_NGS_trimmed_overlap_post_compute_gap.xmap");
	#}
	return $xmap;
}



###########################Functions for outputting agp/fasta or its auxilliary files###########################################


#print the gap informaton from alignments to a file
sub printGapFile{
	my $xmap = $_[0];
	my $out_put = $_[1];
	print "Gap length output file $out_put\n";
	#print "Hello\n";
	open(my $fh, ">".$out_put);
	print $fh "#NGSId1\tNGSId2\tSuperScaffoldId\tXmapGapLength\tAdjustedGapLength\tNGSLength1\tNGSLength2\n";
	my $numAlign = $xmap->{totalHits};
	my $prevInd = 0;
	my $currInd = 1;
	while($currInd < $numAlign){
		print "Index: $prevInd\t$currInd\n";
		if($currInd < $numAlign 
			&& $xmap->{hits}->{RefContigID}[$prevInd] == $xmap->{hits}->{RefContigID}[$currInd]){
			print $fh $xmap->{hits}->{NGSName}[$prevInd]."\t".
				$xmap->{hits}->{NGSName}[$currInd]."\t".
				$xmap->{hits}->{RefContigID}[$currInd]."\t".
				$xmap->{hits}->{GapLen}[$currInd]."\t".
				$xmap->{hits}->{AdjustedGapLen}[$currInd]."\t".
				$xmap->{hits}->{QryLen}[$prevInd]."\t".
				$xmap->{hits}->{QryLen}[$currInd]."\n";
			if(!$xmap->{hits}->{IsEmbedded}[$currInd]){
				$prevInd = $currInd;
			}	
		}else{
			print $fh $xmap->{hits}->{NGSName}[$prevInd]."\t".
				$xmap->{hits}->{NGSName}[$prevInd]."\t".
				$xmap->{hits}->{RefContigID}[$prevInd]."\t".
				#$xmap->{hits}->{GapLen}[$i]."\t".
				#$xmap->{hits}->{AdjustedGapLen}[$i-1]."\t".
				"0\t0\t".
				$xmap->{hits}->{QryLen}[$prevInd]."\t".
				$xmap->{hits}->{QryLen}[$prevInd]."\n";
				$prevInd = $currInd;
		}
		$currInd = $currInd + 1;
	}
}

#print header of AGP file
sub printAGPHeader{
my $fh = $_[0];
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $agp_header = "##agp-version\t2.0\n".
				"# Organism:  \n".
				"# Platform:     \n".
				"# Model:        \n".
				"# Enzyme(s):    \n".
				"# BioSample:    \n".
				"# BioProject:   \n".
             "# Obj_Name\tObj_Start\tObj_End\tPartNum\tCompnt_Type\tCompntId_GapLength\tCompntStart_GapType\tCompntEnd_Linkage\tOrientation_LinkageEvidence\n";
	print $fh $agp_header;
}


sub printAGPAuxHeader{
	my $fh = $_[0];
	my $header= "##agp-version\t2.0\n".
				"# Organism:   \n".
				"# Platform:     \n".
				"# Model:        \n".
				"# Enzyme(s):    \n".
				"# BioSample:    \n".
				"# BioProject:   \n".
				"Obj_Id\tHeadTrimmedLength\tTailTrimmedLength\n";
	print $fh $header;
}



#print the agp files
sub printAGPFile{
	my ($xmap, $out_put, $aux_out_file, $dummyGapLen, $ngs_map)= @_;
	print "AGP output file $out_put\n";
	open(my $fh, ">".$out_put);
	printAGPHeader($fh);
	open(my $aux_fh, ">".$aux_out_file);
	printAGPAuxHeader($aux_fh);
	my $numAlign = $xmap->{totalHits};
	my $currRefStart = 1;
	my $currRefEnd = 1;
	my $Ind = 1;
	
	#outputing ngs contigs
	my $currInd = 0;
	my $nextInd = 1;	
	while($currInd < $numAlign){
		#if contig is embbed inside another contig, skip this contig
		if($xmap->{hits}->{IsEmbedded}[$currInd]){
			#marking the ngs contig as embedded
			$ngs_map->{$xmap->{hits}->{QryContigID}[$currInd]}->[2] = -1;
			$currInd = $currInd + 1;
			next;
		}
		
		if($xmap->{hits}->{IsEmbedded}[$nextInd]){
			#marking the ngs contig as embeddedq
			$ngs_map->{$xmap->{hits}->{QryContigID}[$nextInd]}->[2] = -1;
			$nextInd = $nextInd + 1;
			next;
		}
					
		#printing begin gap for first contig
		#note first contig cannot be embedded so this condition is valid
		if($currInd == 0){
			my $begin_gap = $xmap->{hits}->{RefStartPos}[$currInd]-1;
			print $aux_fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\t$begin_gap\t";
		}
		#printing ngs sequence
		#print "Index $i\n";
		$currRefEnd = $currRefStart + round($xmap->{hits}->{NGSLen}[$currInd]) - 1;
		print $fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\t".
				"$currRefStart\t".	
				$currRefEnd."\t".
				"$Ind\tW\t".
				$xmap->{hits}->{NGSName}[$currInd]."\t".
				#"1\t".$xmap->{hits}->{QryLen}[$currInd]."\t".
				$xmap->{hits}->{NGSStart}[$currInd]."\t".($xmap->{hits}->{NGSStart}[$currInd]+$xmap->{hits}->{NGSLen}[$currInd]-1)."\t".
				$xmap->{hits}->{Orientation}[$currInd]."\n";
		$currRefStart = $currRefEnd + 1;
		#marking in the ngs table that this contig was scaffolded
		if(exists $ngs_map->{$xmap->{hits}->{QryContigID}[$currInd]}){
			$ngs_map->{$xmap->{hits}->{QryContigID}[$currInd]}->[2] = 1;
		}
		
		#printing gap		
		if($nextInd < $numAlign && $xmap->{hits}->{RefContigID}[$currInd] == $xmap->{hits}->{RefContigID}[$nextInd]){
			#my $gapLen = $xmap->{hits}->{GapLen}[$currInd];
			my $gapLen = $xmap->{hits}->{AdjustedGapLen}[$nextInd];
			
			#if two contigs begin/stop at the same label from alignment, we assume it is not overlap
			#if($gapLen == -1){
			#	$gapLen = 0;
			#}
			if($gapLen < 0){
				$gapLen = $dummyGapLen;
			}elsif($gapLen < $dummyGapLen + 10){
				$gapLen = $dummyGapLen + 10
			}
			$currRefEnd = $currRefStart + $gapLen - 1;
			if($gapLen > 0){
				$Ind++;
				print $fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\t".
					"$currRefStart\t".	
					"$currRefEnd\t".
					"$Ind\tN\t".$gapLen."\tscaffold\tyes\tmap\n";
			}
			$currRefStart = $currRefEnd + 1;		
		}else{								
			#handle beginning and end gap sequence
			#end gap for the previous hybrid contig
			my $end_gap = $xmap->{hits}->{RefLen}[$currInd] - $xmap->{hits}->{RefEndPos}[$currInd];
			#print "printing END gap: ".$xmap->{hits}->{RefLen}[$currInd]."\t". $xmap->{hits}->{RefEndPos}[$currInd]."\n";
			print $aux_fh round($end_gap)."\n"; 
			#begin gap for the current hybrid contig
			if($nextInd < $numAlign){ 
				my $begin_gap = $xmap->{hits}->{RefStartPos}[$nextInd]-1;
				print $aux_fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$nextInd]."\t".round($begin_gap)."\t";
			}
			$currRefStart = 1;
			$currRefEnd = 0;
			$Ind = 0;
		}
		$currInd = $nextInd;
		$nextInd = $nextInd + 1;
		$Ind++;
	}
	close($fh);
	close($aux_fh);
}


#convert the interval objects to fasta file format
sub printFasta{
	my ($fasta_map, $xmap, $out_file, $hybridCmap, $cut_site, $paddingGap) = @_;
	my $numAlign = $xmap->{totalHits};
	my $isBegin = 1;
	open(my $fh, ">".$out_file) || die "cannot open output fasta file $out_file\n";
	print "output fasta file $out_file\n";
	my $currInd = 0;
	my $nextInd = 1;
	
	my $minGap = $paddingGap; #min gap is padding + 10base
	for(my $i=0; $i < 10; $i++){
		$minGap = $minGap."N";
	} 
		
	while($currInd < $numAlign){
		#if contig is embbed inside another contig, skip this contig
		if($xmap->{hits}->{IsEmbedded}[$currInd]){
			$currInd = $currInd + 1;
			next;
		}
		
		if($xmap->{hits}->{IsEmbedded}[$nextInd]){
			$nextInd = $nextInd + 1;
			next;
		}

		my $numAlign = $xmap->{totalHits};
		#if contig is embbed inside another contig, skip this contig
		if($xmap->{hits}->{IsEmbedded}[$currInd]){
			next;
		}
		#handling the first line
		if($isBegin){
			print $fh ">Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\n";
			$isBegin=0;
		}
		
		#printing ngs sequence
		print "NGSname $xmap->{hits}->{NGSName}[$currInd]\n";
		my $ngsseq = $fasta_map->($xmap->{hits}->{NGSName}[$currInd]);
		if(1 | $DEBUG){print "NGS seq length ".length($ngsseq)." expected: ".$xmap->{hits}->{NGSLen}[$currInd]."\n"};
		if(isForwardDir($xmap->{hits}->{Orientation}[$currInd])){
			print $fh $ngsseq;
		}else{
			print $fh reverseComplement($ngsseq);
		}

		#printing gap		
		if($nextInd < $numAlign && $xmap->{hits}->{RefContigID}[$currInd] == $xmap->{hits}->{RefContigID}[$nextInd]){
			#my $gapBegin = round($xmap->{hits}->{RefEndPos}[$currInd]);
			#my $gapEnd = round($xmap->{hits}->{RefStartPos}[$nextInd])-1; #gap end at one b.p.before next start
			my $gapBegin = $xmap->{hits}->{AdjustedGapBegin}[$nextInd];
			my $gapEnd = $xmap->{hits}->{AdjustedGapEnd}[$nextInd];
			my $gaplen = $xmap->{hits}->{AdjustedGapLen}[$nextInd];
			
			my $RefContigID = $xmap->{hits}->{RefContigID}[$currInd];
			my $gapSeq = $paddingGap;
			
			if($gaplen >=0){
				if($gaplen < length($minGap)){
					$gapSeq = $minGap
				}else{
			    	#if($gaplen >= 0){  
					$gapSeq = getSeqFromCmap([$gapBegin, $gapEnd], $hybridCmap, $RefContigID, $cut_site);
				}
			}
			
			#when neighboring contigs begin and end at same position, we assume they are not overlap (this by our current definition translate to gaplen of -1)
			#if($gaplen == -1){
			#    $gapSeq = "";
			#}			
			if(1 | $DEBUG){print "Gap seq length: ".length($gapSeq)." expected ".($gapEnd -$gapBegin + 1)."\n"};
			print $fh $gapSeq;
		}else{								
			print $fh "\n";	
			if($nextInd < $numAlign){
				print $fh ">Super-Scaffold_".$xmap->{hits}->{RefContigID}[$nextInd]."\n";
			}
		}
		$currInd = $nextInd;
		$nextInd = $nextInd + 1;
		
	}
	close($fh);
}


#printing out a sequence accoriding to Cmap coordinate
#This is essentially a N-sequence with location of cut sites
#specified by the cut-site sequences of the enzyme/enzymes

sub getSeqFromCmap{
	my @interval = @{$_[0]}; #begin end coordinate of interval
	if($interval[0] > $interval[1]){
		print STDERR "Warning: Begin index is not smaller than end Index in interval: \[$interval[0], $interval[1]\]\n";
		return "";
	}
	my $cmap = $_[1]; #cmap pointer
	my $contigID = $_[2]; #contigID of the reference contig
	#print "hybrid map ID: ".$contigID."\n"; 
	my $contig =$cmap->{contigs}->{$contigID};
	if(!exists $cmap->{contigs}->{$contigID}){
	    print "Cannot find hybrid scaffold Id ".$contigID." cmap file\n";
	    return "";
	}
	
	my $cut_sites = $_[3]; #cut_sites table 
	my $seq="";
	(my $begin, my $end) = getCoordinatesSubset(\@interval, $contig);
	
	#print "hybrid interval $interval[0]\t$interval[1]\n";
	#print "begin $begin $end end\n";
	my $coord = $interval[0];
	for(my $i = $begin; $i < $end; $i++){
		#print "end pos: ".($contig->{Position}[$i])."\n";
		for(; $coord < $contig->{Position}[$i]; $coord++){
			$seq = $seq."N";
		}
		#print "pre-length ".length($seq)."\n";
		my $cut_site = $cut_sites->[$contig->{LabelChannel}[$i]];
		#print "cut site ".($cut_site)."\t at position ".$contig->{Position}[$i]."\n";
		if(!$cut_site){
		    $cut_site = "NNNNNN"; #we cannot lookup this cut_site jus use N to fill gaps
		}
		if($coord > $contig->{Position}[$i]){
		    print "overlap cut sites\n";
		    my $ind = $coord - $contig->{Position}[$i];
		    $cut_site = substr($cut_site, $ind);
		}		
		$seq=$seq.$cut_site;
		#print "length ".length($seq)."\n";
		$coord = $coord+length($cut_site);
	}
	for(; $coord <= $interval[1]; $coord++){
			$seq = $seq."N";
	}
	#adding nick site exceed gap length, need to trim the end, even if losing nicksite
	if($coord > $interval[1] + 1){
	    $seq = substr($seq, 0, length($seq) - ($coord-$interval[1]-1));
	}
	if($DEBUG){print "length of gap: ".(length($seq))."\n"};
	return $seq;
}

#return begin and end index (right exclusive) in a cmap contig within 
#a specified coordinate range note that to avoid scanning 
#the whole cmap data all the time, one can "batch lookup" 
#the coordinates in successive calls of this function on
#many non-overlapping intervals in increasing coordinate order
#the method will remember the last coordinate left-off using
#a static/state variable

{
	my $prevInd = 0;
	my $prevContig;
	sub getCoordinatesSubset{
		my @interval = @{$_[0]};
		my $contig = $_[1];
		my $begin = 0;
		my $end = 0;
		#resetting the remembered index
		if($DEBUG){print "previous stopping at $prevInd\t$prevContig\t$contig\t".($prevContig == $contig)."\n"};
		if(!(defined $prevInd) || !(defined $prevContig) || !($prevContig == $contig)){
			#print "resetting\n";
			$prevInd = 0;
		}
		#user can reset the index and decide to scan from a specified index
		if($#_ == 2){
			$prevInd = $_[2];
		}
		
		my @positions = @{$contig->{Position}};
		for(my $i = $prevInd; $i < $#positions; $i++){
			if($DEBUG){print "$interval[0]\t$interval[1]\t".$positions[$i]."\n"};
			if($positions[$i] > $interval[0] && $begin <= 0){
				$begin = $i;
			}
			
			if($positions[$i] > $interval[1]){
				$end = $i;
				$prevInd = $end-1;
				last;
			}
			
			if($i == $#positions && $begin > 0){
				$end = $i+1;
				$prevInd = $end;
			}
		}
		$prevContig = $contig;
		if($begin > $end){
			$begin = 0;
			$end = 0;
		}
		if($DEBUG){print "Interval $begin $end\n";}
		return ($begin, $end);		
	}
}


sub reverseComplement{
	my $str = $_[0];
	#print "forward: ".substr($str, length($str)-50, 50)."\n";
	my @match_forward = $str =~ /GCTCTTC/g;
	$str =~ tr/ATCGatcg/TAGCtagc/;
	$str = reverse($str);
	#my @match_reverse = $str =~ /GAAGAGC/g;
	#if($#match_forward != $#match_reverse){
	#	print "foward and reverse matches are not the same: $#match_forward\t$#match_reverse\n"
	#}
	#print "rcmp: ".substr(reverse($str), 0, 50)."\n";
	return $str;
}



#printing out the un-used ngs object as "singleton" to agp file as well as fasta
#file
#@TODO when only part of the ngs contig is not used only print out 
#part of the contigs
sub printUnUsedNGS(){
	my $agp_out_file = $_[0];
	my %ngs_map = %{$_[1]};
	my $aux_out_file = $_[2];
	my $fasta_out=$_[3]; 
	my $seqMap = $_[4];
	open(my $fh, ">>".$agp_out_file) || die "Cannot open agp output file";
	open(my $aux_fh, ">>".$aux_out_file) || die "Cannot open agp auxilliary file";
	open(my $fasta_fh, ">".$fasta_out) || die "Cannot open fasta output file";
	
	my $unUsedCount=0;
	foreach my $key(keys %ngs_map){
		my @ngs = @{$ngs_map{$key}};
		my $ngsStart = $ngs[3];
		my $ngsEnd = $ngsStart + $ngs[1] - 1;
		#print "key is: $key\t$ngs[2]\n";
		if($ngs[2] <= 0){
			print $fh "$ngs[0]_obj\t1\t$ngs[1]\t1\tW\t$ngs[0]\t1\t$ngs[1]\t\+\n";
			print $aux_fh "$ngs[0]_obj\t0\t0\n";
			#if($ngs[2] < 0){
				#print $fasta_fh ">".$ngs[0]."_overlap\n";	
			#}else{
			if(length($seqMap->($ngs[0])) > 0){				
				print $fasta_fh ">".$ngs[0]."_obj\n";
				print $fasta_fh $seqMap->($ngs[0])."\n";
			}else{
				warn("Sequence $ngs[0] has zero length, skipping it in output fasta\n");
			}
			$unUsedCount++;
		}
	}
	close($fh);
	close($aux_fh);
	close($fasta_fh);
	#print "Unused NGS: $unUsedCount\n";
}


#this printout the cut coordinate file 
sub printCutCoordinateFile(){
	my $adjust_cut_map = $_[0];
	my $output_file = $_[1];
	open(my $fh, ">".$output_file);
	my $header = 1;
	foreach my $key(sort { $a <=> $b } keys %{$adjust_cut_map}){
		my $entry = $adjust_cut_map->{$key};
		my @keys = ("oldId", "oldStart", "oldEnd", "newId", "newStart", "newEnd");
		if($header){
			print $fh (join "\t", @keys);
			$header=0;
			print $fh "\n";
		}
		my @values;
		foreach my $mykey(@keys){
			push @values, $entry->{$mykey};
		}
		print $fh (join "\t", @values);
		print $fh "\n";
	}
	close($fh);
}

sub printXmapFile(){
	my $xmap = $_[0];
	my $output_file = $_[1];
	my $header = 1;
	open(my $fh, ">".$output_file);
	for(my $i = 0; $i < $xmap->{totalHits}; $i++){
		my @keys = keys%{$xmap->{hits}};
		my @keys = ("XmapEntryID",  "QryContigID", "RefContigID", "QryStartPos", "QryEndPos", "RefStartPos", 
					"RefEndPos", "Orientation",  "Confidence", "QryLen", "RefLen", "AdjustedRefStart", "AdjustedRefEnd", 
					"PrevAdjustedRefStart", "PrevAdjustedRefEnd", "AdjustedGapLen", "IsEmbedded");
		if($header){
			print  $fh "#h ";
			print $fh (join "\t", @keys);
			$header=0;
			print $fh "\n";
		}
		my @values;
		foreach my $key(@keys){
			#print "key ".$key."\n";
			my $value = -1;
			if($xmap->{hits}->{$key}[$i]){
				$value = $xmap->{hits}->{$key}[$i];
			}
			push @values, $value;
		}
		print $fh (join "\t", @values);
		print $fh "\n";
		#print join "\t", (values %{$entry});
		#print "\n";
	}
	close($fh);
}

sub cleanUp(){
	my %INPUTS = %{$_[0]};
	my $dir = $INPUTS{out_dir};
	opendir(DIR, $dir);
	while(my $file = readdir(DIR)){
		if($file =~ /tmp\.fa\.tmp$/){
			print "Removing temp file: ". $file."\n";
			unlink $dir."/".$file;
		}
	}
	close(DIR);
}
