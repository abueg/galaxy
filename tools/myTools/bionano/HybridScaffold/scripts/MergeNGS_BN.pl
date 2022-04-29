# $Id: MergeNGS_BN.pl 10891 2020-04-21 18:28:08Z apang $

#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);
use Getopt::Long;


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
		$module_path2 = $module_path2."/5.14.4";
		$lib4 = $module_path2."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.016000 && $] <= 5.017000) {
		$module_path2 = $module_path2."/5.16.3";
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
	#print "$] version number\n";
	#print "$0 library paths:\n\t"; print join("\n\t", @INC); print "\n";
}

use File::Path qw(make_path);
use File::Copy;
use BNG::Utility;
use BNG::refAlignerRun;
use IPC::Open3;
use IO::Select;
use XML::Simple;

print "\nInfo: Running the command: $0 @ARGV\n";


# declare and read in command line args
my $outputDir = "";
my $refaligner = "";
my $merge_Tvalue = 1e-13;
my $scratchDir = "";
my $logFile = "";
my $ngs_cmap_fn = "";
my $bng_cmap_fn = "";
my $id_shift=100000;
my $readparameters = "";
my $xmlFile = "";

GetOptions      (
	"outputDir=s"		=>	\$outputDir,
	"refaligner=s"		=>	\$refaligner,
	"logFile=s"		=>	\$logFile,
	"ngs_cmap_fn=s"		=>	\$ngs_cmap_fn,
	"bng_cmap_fn=s"		=>	\$bng_cmap_fn,
	"id_shift=i"		=>	\$id_shift,
	"readparameters=s"	=>      \$readparameters,
	"xmlFile=s"		=>	\$xmlFile) or dieLog ("ERROR: MergeNGS_BN: error in command line arguments.\n");

$scratchDir = $outputDir;
$scratchDir .= "/" if ($scratchDir !~ /\/$/);
my $minsites = 0;
my $unMergeableMaxSites = 4;
my $mergeableMinSites = $unMergeableMaxSites + 1;


#####################
# Step 0. Initiation.
#####################
chdir $outputDir;		# change current directory to $outputDir
print "cwd=$outputDir\n";
make_path($scratchDir); 	# make sure $scratchDir exist
open(LOG, ">$logFile");
my $logFH = *LOG;

# 0.1 Read in the BioNanoAssembled Cmap files
# 0.2 Shift BioNano ContigId by offset to prevent ContigID collision,
my $bionano_idshift_cmap_fn = $scratchDir."bionano.idshift.cmap";
my ($oriBngHeadersRef, $oriBngCmapRef) = readWriteCmapShiftId($bng_cmap_fn, $bionano_idshift_cmap_fn, $id_shift);
logMessage($logFH, "Read and write BioNano contigs with ID shift of '" . $id_shift . "' completed and output to '" . $bionano_idshift_cmap_fn, "'.");
logMessage($logFH, "Step 0. Initiation Completed");


# 0.3 Separate the sequence entries into those with too few sites and those with normal number of sites
my $runContigExtract =
	new BNG::refAlignerRun({
		binPath			=>	$refaligner,
		i			=>	$ngs_cmap_fn,
		o			=>	$scratchDir."ngs_contig_minsites_$minsites"."_maxsites_$unMergeableMaxSites",
		merge			=>	1,
		minsites		=>	$minsites,
		maxsites		=>	$unMergeableMaxSites,
		f			=>	1,
		stdout			=>	1,
		stderr			=>	1
});
logMessage($logFH, $runContigExtract->getCMD());
my ($outResults1, $errResults1, $jobStatus1) = $runContigExtract->runCMD();
if ($jobStatus1 != 0)	{
	# print out error
	$errResults1 = "ERROR: 0.3 In separating sequence entries with too few sites, error was encountered with exit code $jobStatus1; out info: $outResults1; error info: $errResults1";
	logMessage($logFH, $errResults1);
	close($logFH);
	dieLog("$errResults1\n");
} # if jobStatus1

$runContigExtract =
	new BNG::refAlignerRun({
		binPath			=>	$refaligner,
		i			=>	$ngs_cmap_fn,
		o			=>	$scratchDir."ngs_contig_minsites_$mergeableMinSites",
		merge			=>	1,
		minsites		=>	$mergeableMinSites,
		f			=>	1,
		stdout			=>	1,
		stderr			=>	1
});
logMessage($logFH, $runContigExtract->getCMD());
($outResults1, $errResults1, $jobStatus1) = $runContigExtract->runCMD();
if ($jobStatus1 != 0)	{
	# print out error
	$errResults1 = "ERROR: 0.3 In extracting sequence entries with sufficient sites for merge, error was encountered with exit code $jobStatus1; out info: $outResults1; error info: $errResults1";
	logMessage($logFH, $errResults1);
	close($logFH);
	dieLog("$errResults1\n");
} # if jobStatus1
my $mergeable_ngs_cmap_fn = $scratchDir . "ngs_contig_minsites_$mergeableMinSites.cmap";        # only look at those NGS that are potentially mergeable

# 0.4 end-mask to mark the entity of Bionano maps and sequences so that pairmerge does not merge between Bionano and Bionano, sequence with sequence
my $masked_bionano_idshift_cmap_fn = $bionano_idshift_cmap_fn;
$masked_bionano_idshift_cmap_fn =~ s/\.cmap$/_entity_masked.cmap/;
markCmapEntity($bionano_idshift_cmap_fn, $masked_bionano_idshift_cmap_fn, $refaligner, "0x8");
$bionano_idshift_cmap_fn = $masked_bionano_idshift_cmap_fn;

my $masked_mergeable_ngs_cmap_fn = $mergeable_ngs_cmap_fn;
$masked_mergeable_ngs_cmap_fn =~ s/\.cmap$/_entity_masked.cmap/;
markCmapEntity($mergeable_ngs_cmap_fn, $masked_mergeable_ngs_cmap_fn, $refaligner, "0x10");
$mergeable_ngs_cmap_fn = $masked_mergeable_ngs_cmap_fn;

my @mrg_rounds_ids = ("A".."Z", "AA".."AZ", "BA".."BZ", "CA".."CZ", "DA".."DZ");

#####################
# Step 1. Merge
#####################
logMessage($logFH, "Step 1. Pair Merge Repeat between NGS\/BioNano started.");
my $mrg_output_prefix = $scratchDir."Mrg";

# parse the XML file
my $XML = new XML::Simple(KeyAttr=>[]);
print STDERR "###XML file=$xmlFile\n";
my $configRef = $XML->XMLin($xmlFile);

my %stageStack = ();
my @commandStack = ();
my $commandStackRef = \@commandStack;

# refAligner
push(@$commandStackRef, $refaligner);
# output
push(@$commandStackRef, "-o");
push(@$commandStackRef, $mrg_output_prefix);
# stdout stderr
push(@$commandStackRef, "-stdout");
push(@$commandStackRef, "-stderr");

# the -i
push(@$commandStackRef, "-i");
push(@$commandStackRef, $bionano_idshift_cmap_fn);
push(@$commandStackRef, "-i");
push(@$commandStackRef, $mergeable_ngs_cmap_fn);

# maxmem maxthreads
my $baseStage = "global";
$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
# minsites
push(@$commandStackRef, "-minsites");
push(@$commandStackRef, $minsites);

# all the parameters from XML file
$baseStage = "mergeNGS_BN";
$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);

# noise parameters
push(@$commandStackRef, "-readparameters");
push(@$commandStackRef, $readparameters);

# add the no same-entity merge criterion
push(@$commandStackRef, ("-TrimmedNoMerge", 9999, 0, "0x18", 4, 5));

# put the final modifications
my $commandStackString = join(" ", @$commandStackRef);
$commandStackString .= " ";
$commandStackString =~ s/-merge_Tvalue/-T/;
$commandStackString =~ s/-id_shift\s+\d+//;
$commandStackString =~ s/-max_merge_rounds\s+\d+//;
@$commandStackRef = split(/\s+/, $commandStackString);

# now the RefAligner call
print "\nRunning command: $commandStackString\n";
my ($outResults0, $errResults0) = runCommand($commandStackRef);

errCheckRefAligner($mrg_output_prefix.".stdout", "END of output", "Pairmerge completed.", "ERROR: Pairmerge cannot be completed.", 1);

@commandStack = ();
$commandStackRef = -1;
$commandStackString = "";

my ($numMrg, $this_mrg_pairs) = parsingFastMrgStdOut("$mrg_output_prefix.stdout");

if ($numMrg == 0)	{
	logMessage($logFH, "1.1 No merge was possible. Stop the merge program.");
	logMessage($logFH, "ERROR: 1.1. There was no merge possible between NGS and BioNano map.");
	close($logFH);
	dieLog("ERROR: 1.1. There are no merge possible between NGS and BioNano map. Please loosen parameters (e.g. -T or -pairmerge) and check input files.");
} # if numMrg

# keep record of which pairs of fragments were merged
writeAllFastMrgPairs($scratchDir."step1.merge.pairs.txt", $this_mrg_pairs, $numMrg, \@mrg_rounds_ids);
logMessage($logFH, "1. We successfully merged $numMrg pairs.");

# now figure out which contigs are hybrid scaffolds, and sequences, and BioNano
my ($usedBioNanoRef, $usedSeqRef, $hybridRef) = determineParticipants($this_mrg_pairs, $numMrg, $id_shift);

# now iterates the output directory to figure out which file is hybrid, sequence leftover or bionano leftover
my ($hybridCmapsRef, $seqLeftOverCmapsRef, $bioNanoLeftOverCmapsRef) = categorizeCmapFiles($scratchDir, "Mrg_contig", $id_shift, $usedBioNanoRef, $usedSeqRef, $hybridRef);
my ($numHybridFiles, $numSeqLeftOverFiles, $numBioNanoLeftOverFiles) = (scalar(keys %$hybridCmapsRef), scalar(keys %$seqLeftOverCmapsRef), scalar(keys %$bioNanoLeftOverCmapsRef));
logMessage($logFH, "1.1 There are $numHybridFiles hybrid scaffold cmap files, $numSeqLeftOverFiles sequence left over cmap files, $numBioNanoLeftOverFiles BioNano left over cmap files");

# now make a single hybrid cmap, a single sequence left over and a single BioNano left over cmap files
logMessage($logFH, "1.2 Concatenating into single hybrid, and single left over cmap files now");
# print a test file
printIdFile($scratchDir."step1.hybrid.id.txt", $hybridCmapsRef);
printIdFile($scratchDir."step1.BN.id.txt", $bioNanoLeftOverCmapsRef);
printIdFile($scratchDir."step1.NGS.id.txt", $seqLeftOverCmapsRef);

my $quickConcatCmap = new BNG::refAlignerRun({
	binPath			=>	$refaligner,
	"if"			=>	$scratchDir."step1.hybrid.id.txt",
	o			=>	$scratchDir."step1.hybrid",
	merge			=>	1,
	stdout			=>	1,
	f			=>	1
});
my $cmd = $quickConcatCmap->getCMD();
logMessage($logFH, $cmd);
my ($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
if ($job_status != 0)	{
	# print out error
	$errResults = "ERROR: 1.2 In concatenating a single hybrid step1.hybrid, error was encountered with exit code=$job_status; out info: $outResults; error info: $errResults";
	logMessage($logFH, $errResults);	
	close($logFH);
	dieLog("$errResults\n");
} # if job_status

$quickConcatCmap = new BNG::refAlignerRun({
	binPath 		=>	$refaligner,
	"if"			=>	$scratchDir."step1.hybrid.id.txt",
	o			=>	$scratchDir."step2.hybrid",
	merge			=>	1,
	stdout			=>	1,
	f			=>	1
});
$cmd = $quickConcatCmap->getCMD();
logMessage($logFH, $cmd);
($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
if ($job_status != 0)	{
	# print out error
	$errResults = "ERROR: 1.2 In concatenating a single hybrid step2.hybrid, error was encountered with exit code=$job_status; out info: $outResults; error info: $errResults";
	logMessage($logFH, $errResults);
	close($logFH);
	dieLog("$errResults\n");
} # if job_status

if ($numBioNanoLeftOverFiles != 0)	{
	$quickConcatCmap = new BNG::refAlignerRun({
		binPath			=>	$refaligner,
		"if"			=>	$scratchDir."step1.BN.id.txt",
		o			=>	$scratchDir."step1.BN.naive",
		merge			=>	1,
		stdout			=>	1,
		f			=>	1
	});
	$cmd = $quickConcatCmap->getCMD();
	logMessage($logFH, $cmd);
	($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
	if ($job_status != 0)	{
		# print out error
		$errResults = "ERROR: 1.2 In concatenating a single step1.BN.naive, error was encountered with exit code=$job_status; out info: $outResults; error info: $errResults";
		logMessage($logFH, $errResults);
		close($logFH);
		dieLog("$errResults\n");
	} # if job_status
} else	{
	# no genome map left over, just create an empty naive file
	open(OUT, ">$scratchDir"."step1.BN.naive.cmap") or dieLog ("ERROR: Cannot write to $scratchDir"."step1.BN.naive.cmap: $!\n");
	close OUT;
} # if numBioNanoLeftOverFiles

if ($numSeqLeftOverFiles != 0)	{
	$quickConcatCmap = new BNG::refAlignerRun({
		binPath			=>	$refaligner,
		"if"			=>	$scratchDir."step1.NGS.id.txt",
		o			=>	$scratchDir."step1.NGS.naive",
		merge			=>	1,
		stdout			=>	1,
		f			=>	1
	});
	$cmd = $quickConcatCmap->getCMD();
	logMessage($logFH, $cmd);
	($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
	if ($job_status != 0)	{
		# print out error
		$errResults = "ERROR: 1.2 In concatenating a single step1.NGS.naive, error was encountered with exit code=$job_status; out info: $outResults; error info: $errResults";
		logMessage($logFH, $errResults);
		close($logFH);
		dieLog("$errResults\n");
	} # if job_status
} else	{
	# no sequence left over, just create an empty naive file
	open(OUT, ">$scratchDir"."step1.NGS.naive.cmap") or dieLog ("ERROR: Cannot write to $scratchDir"."step1.NGS.naive.cmap: $!\n");
	close OUT;	
} # if numSeqLeftOverFiles

logMessage($logFH, "Step 1.2 file concatenation completed successfully\n");

close($logFH);
print "$0 finished successfully.\n\n";

### subroutines ###

sub readWriteCmapShiftId	{
	my ($file, $outFile, $idShift) = @_;
	
	open(IN, $file) or dieLog("ERROR: Unable to read in file $file: $!\n");
	open(OUT, ">$outFile") or dieLog("ERROR: Unable to write to $outFile: $!\n");
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;
		
		if ($line =~ /^#/)	{
			# header line
			print OUT "$line\n";
			
			next;
		} # if line
		
		# data line
		my @content = split(/\t/, $line);
		my ($cmapId) = ($content[0]);
		my $newId = $cmapId + $idShift;
		
		print OUT join("\t", $newId, @content[1..$#content])."\n";
		
	} # while line
	close OUT;
	close IN;
	
} # readCmapShiftId

# this subroutine marks the Bionano maps with a Mask bit, or marks the sequence with a different Mask bit 
sub markCmapEntity	{
	my ($file, $outFile, $refAligner, $maskBit) = @_;
	
	my $outFilePrefix = $outFile;
	$outFilePrefix =~ s/\.cmap$//;
	
	my @commandStack = ($refAligner, "-o", $outFilePrefix, "-stdout", "-stderr", "-i", $file, "-merge", "-f", "-orMask", $maskBit);
	my ($outResults, $errResults) = runCommand(\@commandStack);
	errCheckRefAligner("$outFilePrefix.stdout", "END of output", "Entity mask of $file to $outFile completed.", "ERROR: Entity mask of $file to $outFile cannot be completed.", 1);
	
} # markCmapEntity

sub printIdFile	{
	my ($file, $idsRef) = @_;
	open(OUT, ">$file") or dieLog("ERROR: printIdFile: cannot write to $file: $!\n");
	foreach my $id (sort numeric keys %$idsRef)	{
		print OUT "$idsRef->{$id}\n";	# print the file name
	} # foreach id
	close OUT;
} # printIdFile

sub numeric	{	$a	<=>	$b	}

sub categorizeCmapFiles	{
	my ($dir, $filePrefix, $idShift, $usedBioNanoRef, $usedSeqRef, $hybridRef) = @_;
	my %hybridCmaps = ();
	my %seqLeftOverCmaps = ();
	my %bioNanoLeftOverCmaps = ();
	opendir(DIR, $dir) or dieLog "ERROR: categorizeCmapFiles: cannot open dir $dir: $!\n";
	my @cmapFiles = grep {$_ =~ /^$filePrefix\d+\.cmap$/i} readdir DIR;
	closedir DIR;

	for (my $i = 0; $i < scalar(@cmapFiles); $i += 1)	{
		my $theId = $cmapFiles[$i];	$theId =~ s/^$filePrefix//;	$theId =~ s/\.cmap$//;
		if (exists $hybridRef->{$theId})	{
			# this file is a hybrid
			$hybridCmaps{$theId} = $dir."$cmapFiles[$i]";
		} else	{
			# this file is a left over, check whether it is a sequence or a BioNano
			if ($theId < $idShift)	{
				$seqLeftOverCmaps{$theId} = $dir."$cmapFiles[$i]";
			} else	{
				$bioNanoLeftOverCmaps{$theId} = $dir."$cmapFiles[$i]";
			} # if theId
		} # if exists
	} # for i	
	return (\%hybridCmaps, \%seqLeftOverCmaps, \%bioNanoLeftOverCmaps);
} # categorizeCmapFiles

sub determineParticipants	{
	my ($mrgPairsRef, $numMrg, $idShift) = @_;
	my %usedBioNano = ();
	my %usedSeq = ();
	my %hybrid = ();	# stores hybrid ids, but would not store those that are not present in the output directory in the end, BECAUSE they have been merged with an entity either with a different id or was totally encompassed

	for (my $i = 0; $i < $numMrg; $i += 1)	{
		my ($contig1, $contig2, $theHybrid) = ($mrgPairsRef->{ContigID1}[$i], $mrgPairsRef->{ContigID2}[$i], $mrgPairsRef->{ResultContigID}[$i]);
		if ($contig1 < $idShift)	{
			# a sequence
			$usedSeq{$contig1} = 1;
		} else	{
			# a genome map
			$usedBioNano{$contig1} = 1;
		} # if mrgPairsRef

		if ($contig2 < $idShift)	{
			# a sequence
			$usedSeq{$contig2} = 1;
		} else	{
			# a genome map
			$usedBioNano{$contig2} = 1;
		} # if mrgPairsRef

		$hybrid{$theHybrid} = 1;

		# now figure out if the particpant was a hybrid already AND that this merge resulted in an hybrid with a different id (either be a different id, or that the participant was completely encompassed)
		delete $hybrid{$contig1} if (exists $hybrid{$contig1} && $theHybrid != $contig1);
		delete $hybrid{$contig2} if (exists $hybrid{$contig2} && $theHybrid != $contig2);
	} # for i
	return (\%usedBioNano, \%usedSeq, \%hybrid);
} # determineParticipants

sub runCommand  {
	my ($argsRef) = @_;
	my $pid = open3(my $CMD_IN, my $out, my $err, @$argsRef);

	close($CMD_IN);

	my $outResults = "";
	my $errResults = "";
	my $sel = new IO::Select;
	$sel->add($out, $err);
	while(my @fhs = $sel->can_read) {
		foreach my $fh (@fhs) {
			my $line = <$fh>;
			unless(defined $line) {
				$sel->remove($fh);
				next;
			} # unless line
			if($fh == $out) {
				$outResults .= "$line";
			}elsif($fh == $err) {
				$errResults .= "$line";
			}else{
				dieLog ("ERROR: This should never execute!");
			} # if fh
		} # foreach fh
	} # while
	my $ret=waitpid ($pid, 0); # reap the exit code
	return ($outResults, $errResults);
} # runCommand

sub errCheckRefAligner  {
	my ($file, $completeString, $completeMsg, $dieMsg, $exitFlag) = @_;
	open(IN, "$file") or dieLog ("ERROR: Cannot open $file: $!\n");
	my $presence = 0;
	while (my $line = <IN>) {
		if ($line =~ /$completeString/) {
			# if the line contains the string that indicates successful completion
			$presence = 1;
		} # if line
	} # while line
	close IN;

	return $presence if ($exitFlag != 1);

	if ($presence == 0)     {
		dieLog ("ERROR: $dieMsg\n");
	} else  {
		print "$completeMsg\n";
	} # if presence
	return $presence;
} # errCheckRefAligner

sub parseConfig	{
	my ($configRef, $theStage, $commandStackRef) = @_;
	foreach my $flagIncludeKey (keys $configRef->{$theStage})	{
		if ($flagIncludeKey =~ /flag/)	{
			# inside flag section
			# check to see if there is just one flag, in which case under flag is a hash table; if more than one entry under flag, then an array
			if (ref($configRef->{$theStage}{$flagIncludeKey}) eq "ARRAY")	{
			# iterate the array of attributes
			for (my $i = 0; $i < scalar(@{$configRef->{$theStage}{$flagIncludeKey}}); $i += 1)	{
				my $tempAttr = "";
				my %tempVals = ();
				foreach my $flagKey (keys %{$configRef->{$theStage}{$flagIncludeKey}[$i]})	{
					# capture the attribute and values
					if ($flagKey =~ /attr/)	{
						# check if attr has "-" in front
						$tempAttr = ($configRef->{$theStage}{$flagIncludeKey}[$i]{$flagKey} !~ /^-/) ? ("-$configRef->{$theStage}{$flagIncludeKey}[$i]{$flagKey}") : ($configRef->{$theStage}{$flagIncludeKey}[$i]{$flagKey});
					} elsif ($flagKey =~ /val(\d+)/)	{
						$tempVals{$1} = $configRef->{$theStage}{$flagIncludeKey}[$i]{$flagKey};
					} else	{
						# no op, only want attr and val\d+
					} # if flagKeys
				} # foreach flagKeys
				# at the end of the XML line, now assign the sorted values to attr
				if ($tempAttr ne "")	{
					# make sure there is an attribute to store
					push(@$commandStackRef, $tempAttr);
					foreach my $vals(sort numeric keys %tempVals)	{
						push(@$commandStackRef, $tempVals{$vals});
					} # foreach vals
				} # if tempAttr
			} # for i
			} elsif (ref($configRef->{$theStage}{$flagIncludeKey}) eq "HASH")	{
				my $tempAttr = "";
				my %tempVals = ();
				foreach my $flagKey (keys %{$configRef->{$theStage}{$flagIncludeKey}})	{
					# capture the attribute and values
					if ($flagKey =~ /attr/)	{
						# check if attr has "-" in front
						$tempAttr = ($configRef->{$theStage}{$flagIncludeKey}{$flagKey} !~ /^-/) ? ("-$configRef->{$theStage}{$flagIncludeKey}{$flagKey}") : ($configRef->{$theStage}{$flagIncludeKey}{$flagKey});
					} elsif ($flagKey =~ /val(\d+)/)	{
						$tempVals{$1} = $configRef->{$theStage}{$flagIncludeKey}{$flagKey};
					} else	{
						# no op, only want attr and val\d+
					} # if flagKey
				} # foreach flagKey
				if ($tempAttr ne "")	{
					# make sure there is an attribute to store
					push(@$commandStackRef, $tempAttr);
					foreach my $vals (sort numeric keys %tempVals)	{
						push(@$commandStackRef, $tempVals{$vals});
					} # foreach vals
				} # if tempAttr
			} else	{
				# op operation, as it must be an array or hash under flag
			} # if under flag is of type array or hash
		} # if flag
		if ($flagIncludeKey =~ /include/)	{
			# inside include section
			# check if there are multiple includes, in which case it will be an array under include; single include, then a hash under include	
			if (ref ($configRef->{$theStage}{$flagIncludeKey}) eq "ARRAY")	{
				# iterate each include
				for (my $i = 0; $i < scalar(@{$configRef->{$theStage}{$flagIncludeKey}}); $i += 1)	{
					foreach my $includeKey (keys %{$configRef->{$theStage}{$flagIncludeKey}[$i]})	{
						if ($includeKey =~ /val(\d+)/)	{
							# find out which stage to include, and recursively parse that stage
							my $childStage = $configRef->{$theStage}{$flagIncludeKey}[$i]{$includeKey};
							$commandStackRef = parseConfig($configRef, $childStage, $commandStackRef);
						} # if includeKey
					} # foreach includeKey
				} # for i
			} elsif (ref ($configRef->{$theStage}{$flagIncludeKey}) eq "HASH")	{
				# only 1 stage to include
				foreach my $includeKey (keys %{$configRef->{$theStage}{$flagIncludeKey}})	{
					if ($includeKey =~ /val(\d+)/)	{
						# find out which stage to include, and recursively parse that stage
						my $childStage = $configRef->{$theStage}{$flagIncludeKey}{$includeKey};
						$commandStackRef = parseConfig($configRef, $childStage, $commandStackRef);
					} # if includeKey
				} # foreach includeKey
			} else	{
				# no op, as it must an array or hash under include
			} # if under include is of type array or hash type
		} # if include 
	} # foreach flagIncludKeys
	return $commandStackRef
} # parseConfig
