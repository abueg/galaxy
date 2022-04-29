# $Id: align_final_seq_two_passes.pl 7455 2018-03-08 18:12:27Z jwang $

#!/usr/bin/perl -w

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
use Getopt::Long;

print "\nInfo: Running the command $0 @ARGV\n";

# declear and read in command line args
my $outDir = "";
my $outFilePrefix = "";
my $refAligner = "";
my $xmlFile = "";
my $seqFile = "";
my $scaffoldFile = "";
my $logFile = "";

GetOptions	(
	"outDir=s"		=>	\$outDir,
	"outFilePrefix=s"	=>	\$outFilePrefix,
	"refAligner=s"		=>	\$refAligner,
	"xmlFile=s"		=>	\$xmlFile,
	"seqFile=s"		=>	\$seqFile,
	"scaffoldFile=s"	=>	\$scaffoldFile,
	"logFile=s"		=>	\$logFile
) or dieLog("ERROR: align_final_seq_two_passes: error in command line arguments.\n");

# final CMAP and XMAP file names
mkdir $outDir if (! -e $outDir);
my $outSeqCmapFile = "$outDir/$outFilePrefix"."_q.cmap";
my $outHybridCmapFile = "$outDir/$outFilePrefix"."_r.cmap";
my $outXmapFile = "$outDir/$outFilePrefix.xmap";

# intermediate file names
my $outFile1stPassPrefix = "$outFilePrefix"."_1st_pass";
my $outFile2ndPassPrefix = $outFilePrefix."_2nd_pass"; 
my $outFileNotFilteredPrefix = "$outFilePrefix"."_not_filtered";

open(LOG_FILE, ">$logFile") or dieLog("ERROR: align_final_two_passes: cannot open log file $logFile: $!\n");
my $logFileHandle = *LOG_FILE;

### 1st pass ###
### align the input sequences (post-conflict resolved) to the hybrid scaffold 
### use the sequence as the reference, and the hybrid scaffold as the query
logMessage($logFileHandle, "Running 1st pass alignment of sequence to hybrid scaffold...\n");
my @baseStages = ("global", "align_final_1st_pass");
my $baseStagesRef = \@baseStages;
runAlignment($refAligner, $scaffoldFile, $seqFile, "$outDir/$outFile1stPassPrefix", $baseStagesRef, $logFileHandle);

my %usedSeq1stPass = ();		my $usedSeq1stPassRef = \%usedSeq1stPass;
my %usedScaffold1stPass = ();		my $usedScaffold1stPassRef = \%usedScaffold1stPass;
my %usedSeq2ndPass = ();		my $usedSeq2ndPassRef = \%usedSeq2ndPass;
my %usedScaffold2ndPass = ();		my $usedScaffold2ndPassRef = \%usedScaffold2ndPass;

($usedSeq1stPassRef, $usedScaffold1stPassRef) = getUsedSeqScaffold("$outDir/$outFile1stPassPrefix.xmap", $usedSeq1stPassRef, $usedScaffold1stPassRef);
my (undef, $numTotalSeq, undef) = readCMap($seqFile);

logMessage($logFileHandle, "1st pass alignment of sequence to hybrid scaffold completed successfully\n");


if ($numTotalSeq == scalar(keys %$usedSeq1stPassRef))	{
	# if every sequence is aligned in 1st pass
	# then no 2nd pass is run, so combine, flip, output on just the 1st pass result
	combineSwitchReferenceQueryXmap("$outDir/$outFile1stPassPrefix.xmap", "", "$outDir/$outFileNotFilteredPrefix.xmap", $outHybridCmapFile, $outSeqCmapFile);
} else	{
	### 2nd pass ###
	# find the left over sequence and align using the leftoverSequence as the reference, and the hybrid scaffold as the query
	logMessage($logFileHandle, "Running 2nd pass alignment of sequence to hybrid scaffold...\n");
	my $notUsedSeqCmapFilePrefix = $outFile1stPassPrefix."_not_used_seq";
	findLeftOverSequences($seqFile, "$outDir/$notUsedSeqCmapFilePrefix", $usedSeq1stPassRef, $refAligner, $logFileHandle);
	$baseStagesRef->[1] = "align_final_2nd_pass";
	runAlignment($refAligner, $scaffoldFile, "$outDir/$notUsedSeqCmapFilePrefix.cmap", "$outDir/$outFile2ndPassPrefix", $baseStagesRef, $logFileHandle);

	($usedSeq2ndPassRef, $usedScaffold2ndPassRef) = getUsedSeqScaffold("$outDir/$outFile2ndPassPrefix.xmap", $usedSeq2ndPassRef, $usedScaffold2ndPassRef);

	if (scalar(keys %$usedSeq2ndPassRef) == 0)	{
		# if no sequence can be aligned in 2nd pass
		if (scalar(keys %$usedSeq1stPassRef) == 0)	{
			# if no sequence was aligned 1st pass, ERROR, nothing was aligned
			logMessage($logFileHandle, "ERROR: No alignment was possible in both align final 1st and 2nd pass between sequence and hybrid scaffold.");
			close($logFileHandle);
			dieLog("ERROR: No alignment was possible in both align final 1st and 2nd pass between sequence and hybrid scaffold.");
		} else	{
			# some sequence was aligned in 1st pass, but none in 2nd pass
			# so combine, flip, output on just the 1st pass result	
			combineSwitchReferenceQueryXmap("$outDir/$outFile1stPassPrefix.xmap", "", "$outDir/$outFileNotFilteredPrefix.xmap", $outHybridCmapFile, $outSeqCmapFile);
		} # if no sequences aligned in 1st pass
	} else	{
		# some sequence can be aligned in 2nd pass
		if (scalar(keys %$usedSeq1stPassRef) == 0)	{
			# if no sequence was aligned in 1st pass, but some in 2nd pass
			# so combine, flip, output just on the 2nd pass
			combineSwitchReferenceQueryXmap("", "$outDir/$outFile2ndPassPrefix.xmap", "$outDir/$outFileNotFilteredPrefix.xmap", $outHybridCmapFile, $outSeqCmapFile);
		} else	{
			# some sequence was aligned in 1st and 2nd pass
			# combine, flip, output, on both 1st and 2nd pass
			combineSwitchReferenceQueryXmap("$outDir/$outFile1stPassPrefix.xmap", "$outDir/$outFile2ndPassPrefix.xmap", "$outDir/$outFileNotFilteredPrefix.xmap", $outHybridCmapFile, $outSeqCmapFile);
		} # if no sequence aligned in 1st pass

	} # if no sequences aligned in 2nd pass
	
	logMessage($logFileHandle, "2nd pass alignment of sequence to hybrid scaffold completed successfully\n");
} # if numTotalSeq


### print the used sequences CMAP ###
logMessage($logFileHandle, "Running combination of 1st and 2nd pass alignment and switch between reference and query...\n");

### find best location for each sequence then output ###
logMessage($logFileHandle, "Running selection of best reference position...\n");
my ($bestAlignmentPairsHeadersRef, $bestAlignmentPairsRef) = getBestAlignmentPairsXmap("$outDir/$outFileNotFilteredPrefix.xmap");
printBestAlignmentPairsXmap($outXmapFile, $bestAlignmentPairsHeadersRef, $bestAlignmentPairsRef);
logMessage($logFileHandle, "Selection of best reference position finished sucessfully.");

my @alignmentSeqFiles = ("$outDir/$outFile1stPassPrefix"."_r.cmap", "$outDir/$outFile2ndPassPrefix"."_r.cmap");
my $usedAllSeqRef = ();
my $usedAllScaffoldRef = ();
($usedAllScaffoldRef, $usedAllSeqRef) = getUsedSeqScaffold($outXmapFile, $usedAllScaffoldRef, $usedAllSeqRef);
printUsedSeqScaffold(\@alignmentSeqFiles, $outSeqCmapFile, $usedAllSeqRef);
### print the used hybrid scaffold CMAP ###
my @alignmentScaffoldFiles = ("$outDir/$outFile1stPassPrefix"."_q.cmap", "$outDir/$outFile2ndPassPrefix"."_q.cmap");
printUsedSeqScaffold(\@alignmentScaffoldFiles, $outHybridCmapFile, 0);
logMessage($logFileHandle, "Combination of 1st and 2nd pass alignment and switch between reference and query finished successfully.");



logMessage($logFileHandle, "Two step alignment of sequence to hybrid scaffold completed successfully\n");

close($logFileHandle);
print "$0 finished successfully.\n\n";


### subroutines ###
sub printBestAlignmentPairsXmap	{
	my ($file, $headersRef, $bestDataRef) = @_;
	
	open(OUT, ">$file") or dieLog("ERROR: printBestAlignmentPairsXmap: cannot write to $file: $!\n");
	# header lines
	for (my $i = 0; $i < scalar(@$headersRef); $i += 1)	{
		print OUT "$headersRef->[$i]\n";
	} # for i
	# data lines
	foreach my $xId (sort numeric keys %$bestDataRef)	{
		print OUT "$bestDataRef->{$xId}[0]{line}\n";
	} # foreach xId
	close OUT;
} # printBestAlignmentPairsXmap

sub getBestAlignmentPairsXmap	{
	my ($file) = @_;
	my %data = ();
	my @headers = ();
	open(IN, $file) or dieLog ("getBestAlignmentPairsXmap: cannot read in $file: $!\n");
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;	# trim whitespace
		
		# header line
		if ($line =~ /^#/)	{
			push(@headers, $line);
			next;
		} # if line

		# data line
		my @content = split(/\t/, $line);
		my ($xId, $qryId, $refId, $confidence) = ($content[0], $content[1], $content[2], $content[8]);
		push(@{$data{"$qryId"}}, {confidence => $confidence, refId => $refId, qryId => $qryId, line => $line, xId => $xId});
	} # while line
	close IN;


	# sort the aligned reference query pairs by scores
	my $dataRef = sortReverseByHashKeyValue(\%data, "confidence");


	# return only the top alignment for each reference query pair
	my %bestAlignmentPairs = ();
	foreach my $id (keys %$dataRef)	{
		my $aRef = $dataRef->{$id}[0];	# only the top alignment
		push(@{$bestAlignmentPairs{$aRef->{xId}}}, {confidence => $aRef->{confidence}, refId => $aRef->{refId}, qryId => $aRef->{qryId}, line => $aRef->{line} });
	} # foreach idPair
	return (\@headers, \%bestAlignmentPairs);
} # getBestAlignmentPairsXmap

sub sortReverseByHashKeyValue	{
	my ($dataRef, $theKey) = @_;
	foreach my $id (keys %$dataRef)	{
		@{$dataRef->{$id}} = sort	{
			$b->{$theKey}	<=>	$a->{$theKey}
		} @{$dataRef->{$id}}
	} # foreach id
	return $dataRef;
} # sortByHashKeyValue

sub printUsedSeqScaffold	{
	my ($alignmentFilesRef, $outFile, $IDs) = @_;
	
	my %usedIds = ();
	my %usedHeaders = ();
	my %usedData = ();
	
	if($IDs != 0){
		my @ids =  (keys %{$IDs});
		print "Number of IDs: $#ids\n";
	}
	
	for (my $i = 0; $i < scalar(@$alignmentFilesRef); $i += 1)	{
		my $inFile = $alignmentFilesRef->[$i];
		next if (! -e $inFile);	# if not exist that file
		open(IN, $inFile) or dieLog("ERROR: printUsedSeqScaffold: cannot open $inFile: $!\n");
		while (my $line = <IN>)	{
			chomp $line;
			$line =~ s/\r+//g;
			$line =~ s/^\s+|\s+$//g;	# trim whitespace
			
			# header line, only keep the necessary lines
			if ($line =~ /^#/)	{
				if ($line =~ /# CMAP File Version:/ || $line =~ /# Label Channels:/ || $line =~ /# Nickase Recognition Site/ || $line =~ /# Values corresponding/ || $line =~ /#h/ || $line =~ /#f/)	{
					push(@{$usedHeaders{$i}}, $line);

				} elsif ($line =~ /(# Number of Consensus Maps:)/)	{
					### update this during output!!!
					push(@{$usedHeaders{$i}}, $line);
				} else	{
					# no operation
				} # if line
				next;		
			} # if line
			
			# data line, record the lines, and record which inFile this comes from
			my @content = split(/\t/, $line);
			my $theCmapId = $content[0];
			if($IDs == 0 || $IDs->{$theCmapId}){
				push(@{$usedData{$theCmapId}{$i}}, $line);
				
				# keep a record of which input cmap file (index) this entry comes from 
				if (! exists $usedIds{$theCmapId})	{
					$usedIds{$theCmapId} = $i;
				} # if exists
			}
		} # while line
		
		close IN;
		
	} # for i
	
	# now print out 
	open(OUT, ">$outFile") or dieLog("ERROR: printUsedSeqScaffold: cannot write to $outFile: $!\n");
	# header lines
	foreach my $inFileIndex (keys %usedHeaders)	{
		for (my $i = 0; $i < scalar(@{$usedHeaders{$inFileIndex}}); $i += 1)	{
			my $theHeader = $usedHeaders{$inFileIndex}[$i];
			
			if ($theHeader =~ /(# Number of Consensus Maps:)/)	{
				# update the number of maps
				print OUT $1."\t".(scalar keys %usedIds)."\n";
			} else	{
				print OUT "$theHeader\n";
			} # if theHeader
		} # for i
		last;	# just print one header
	} # foreach inFileIndex
	
	# data lines
	foreach my $theCmapId (sort numeric keys %usedIds)	{
		# print out the lines for each theCmapId (from the appropriate input cmap file)
		for (my $i = 0; $i < scalar(@{$usedData{$theCmapId}{$usedIds{$theCmapId}}}) ; $i += 1)	{
			my $theData = $usedData{$theCmapId}{$usedIds{$theCmapId}}[$i];
			print OUT "$theData\n";
		} # for i
	} # foreach theId
	close OUT;
	
} # printUsedSeqScaffold

sub combineSwitchReferenceQueryXmap	{
	my ($xmapFile1, $xmapFile2, $finalXmapFile, $finalReferenceCmapFile, $finalQueryCmapFile) = @_;

	my %xmapHeaders = ();
	my %xmapContent = ();
	my $xmapHeadersRef = \%xmapHeaders;	# two dimensional arrays
	my $xmapContentRef = \%xmapContent;

	my $numIteration = 0;
	if ($xmapFile1 ne "" && $xmapFile2 ne "")	{
		$numIteration = 2;
	} elsif (($xmapFile1 eq "" && $xmapFile2 ne "") || ($xmapFile1 ne "" && $xmapFile2 eq "") )	{
		$numIteration = 1;
	} else	{
		# no op, should not even be here
	} # if 

	foreach my $xmapFile ("$xmapFile1", "$xmapFile2")	{
		next if ($xmapFile eq "");
		# initialize 
		$xmapHeadersRef->{$xmapFile} = [];
		$xmapContentRef->{$xmapFile} = [];
		($xmapHeadersRef->{$xmapFile}, $xmapContentRef->{$xmapFile}) = readAlignFinalXmap($xmapFile, $xmapHeadersRef->{$xmapFile}, $xmapContentRef->{$xmapFile});
	} # foreach xmapFile


	# combine the xmap content
	my @outContent = ();
	my $outContentRef = \@outContent;
	$outContentRef = combineAlignFinalXmap($xmapContentRef, $outContentRef);

	# print new XMAP file
	my $theXmapFile = "";
	foreach my $xmapFile ("$xmapFile1", "$xmapFile2")	{
		next if ($xmapFile eq "");
		$theXmapFile = $xmapFile;
		last;
	} # foreach xmapFile
	printNewXmap($finalXmapFile, $xmapHeadersRef->{$theXmapFile}, $finalReferenceCmapFile, $finalQueryCmapFile, $outContentRef);
	
} # combineSwitchReferenceQueryXmap

sub printNewXmap	{
	my ($file, $headersRef, $refCmapFile, $qryCmapFile, $contentRef) = @_;
	open(OUT, ">$file") or dieLog("ERROR: printNewXmap: cannot write to $file: $!\n");
	# header lines
	for (my $i = 0; $i < scalar(@$headersRef); $i += 1)	{
		my $line = $headersRef->[$i];
		if ($line =~ /(# Reference Maps From:)/)	{
			$line = "$1\t$refCmapFile";
		} # if line
		if ($line =~ /(# Query Maps From:)/)	{
			$line = "$1\t$qryCmapFile";
		} # if line
		print OUT "$line\n";
	} # for i
	# data lines	
	
	for (my $i = 0; $i < scalar(@$contentRef); $i += 1)	{
		my $aRef = $contentRef->[$i];
		my $line = join("\t", $aRef->{xId}, $aRef->{qryId}, $aRef->{refId})."\t";
		$line .= join("\t", $aRef->{qryStart}, $aRef->{qryEnd}, $aRef->{refStart}, $aRef->{refEnd})."\t";
		$line .= join("\t", $aRef->{orientation}, $aRef->{confidence}, $aRef->{cigarString});

		# if a newer XMAP format
		if ($aRef->{labelChannel} != -1 && $aRef->{alignmentString} ne "")	{
			$line .= "\t";
			$line .= join("\t", $aRef->{qryLength}, $aRef->{refLength}, $aRef->{labelChannel}, $aRef->{alignmentString});
		} # if newer XMAP format

		# if a newer XMAP format
		if (! floatEqual($aRef->{mapWt}, -1.0, 6 ))	{
			$line .= "\t";
			$line .= join("\t", $aRef->{mapWt});
		} # if floatEqual

		print OUT "$line\n";
	} # for i

	close OUT;
} # printNewXmap

sub floatEqual	{
	my ($a, $b, $dp) = @_;
	return sprintf("%.${dp}g", $a) eq sprintf("%.${dp}g", $b);
} # floatEqual

sub combineAlignFinalXmap	{
	my ($contentRef, $outContentRef) = @_;

	my ($maxXmapId, $prevMaxXmapId) = (0, 0);
	foreach my $xmapFile (keys %$contentRef)	{
		for (my $j = 0; $j < scalar(@{$contentRef->{$xmapFile}}); $j += 1)	{
			my $aRef = $contentRef->{$xmapFile}[$j];
			# new XMAP id number
			$maxXmapId = ($maxXmapId < $aRef->{xId}) ? ($aRef->{xId}) : ($maxXmapId);
			my $newXId = $prevMaxXmapId + $aRef->{xId};
			# swap query and reference id
			my ($newQryId, $newRefId) = swap($aRef->{qryId}, $aRef->{refId});
			# swap the co-ordinates
			my ($newQryStart, $newQryEnd, $newRefStart, $newRefEnd) = swapPositions($aRef->{qryStart}, $aRef->{qryEnd}, $aRef->{refStart}, $aRef->{refEnd}, $aRef->{orientation});
			# cigar string
			my $cigarContentRef = findCigarMatchPositions('\d+M|\d+I|\d+D', $aRef->{cigarString});
			my $newCigarString = swapCigar($cigarContentRef, $aRef->{orientation});

			# if this is a newer XMAP format
			my ($newQryLength, $newRefLength, $newAlignmentString) = (-1.0, -1.0, -1, "");
			if ($aRef->{labelChannel} != -1 && $aRef->{alignmentString} ne "")	{
				# swap query and reference length
				($newQryLength, $newRefLength) = swap($aRef->{qryLength}, $aRef->{refLength});
				# alignment string
				my $alignmentContentRef = findAlignmentMatchPositions('\)\(', $aRef->{alignmentString});
				$newAlignmentString = swapAlignment($alignmentContentRef, $aRef->{orientation});
			} # if a newer XMAP format
			
			# store the new record
			push(@$outContentRef, {xId => $newXId, qryId => $newQryId, refId => $newRefId, qryStart => $newQryStart, qryEnd => $newQryEnd, refStart => $newRefStart, refEnd => $newRefEnd, orientation => $aRef->{orientation}, confidence => $aRef->{confidence}, cigarString => $newCigarString, qryLength => $newQryLength, refLength => $newRefLength, labelChannel => $aRef->{labelChannel}, alignmentString => $newAlignmentString, mapWt => $aRef->{mapWt}});	
		} # for j
		# update the prevMaxXmapId by the current maxXmapId
		$prevMaxXmapId = $maxXmapId;
	} # for i
	return $outContentRef;	
} # combineAlignFinalXmap

sub swapAlignment	{
	my ($alignmentContentRef, $orientation) = @_;
	my @newAlignmentContent = ();
	for (my $i = 0; $i < scalar(@$alignmentContentRef); $i += 1)	{
		# swap the ref and qry
		($alignmentContentRef->[$i]{refIndex}, $alignmentContentRef->[$i]{qryIndex}) = swap($alignmentContentRef->[$i]{refIndex}, $alignmentContentRef->[$i]{qryIndex});
	} # for i
	# if negative orientation, then reverse the order of the duplet
	if ($orientation =~ /^-$/)	{
		for (my $i = $#{$alignmentContentRef}; $i >= 0; $i -= 1)	{
			push(@newAlignmentContent, {refIndex => $alignmentContentRef->[$i]{refIndex}, qryIndex => $alignmentContentRef->[$i]{qryIndex}});
		} # for i
	} else	{
		for (my $i = 0; $i <= $#{$alignmentContentRef}; $i += 1)	{
			push(@newAlignmentContent, {refIndex => $alignmentContentRef->[$i]{refIndex}, qryIndex => $alignmentContentRef->[$i]{qryIndex}});
		} # for i
	} # orientation

	my $newAlignmentString = "";
	for (my $i = 0; $i < scalar(@newAlignmentContent); $i += 1)	{
		$newAlignmentString .= "$newAlignmentContent[$i]{refIndex},$newAlignmentContent[$i]{qryIndex}".")(";
	} # for i
	$newAlignmentString =~ s/\($//;
	$newAlignmentString = "(".$newAlignmentString;
	return $newAlignmentString;
} # swapAlignment

sub findAlignmentMatchPositions	{
	my ($regex, $string) =@_;
	my @locations = ();
	# first remove the leading ( and trailing )
	$string =~ s/^\(|\)$//g;
	# identify and parse the separator )(
	my @temp = split(/$regex/, $string);
	# split into dublets
	for (my $i = 0; $i < scalar(@temp); $i += 1)	{
		my ($refIndex, $qryIndex) = split(/,/, $temp[$i]);
		push(@locations, {refIndex => $refIndex, qryIndex => $qryIndex});
	} # for i
	return \@locations;
} # findAlignmentMatchPositions

sub swapCigar	{
	# this subroutine changes insertions to deletions and vice versa
	# if the alignment orientation is negative, then it would reverse the tokens
	my ($cigarContentRef, $orientation) = @_;
	my @newCigarContent = ();
	for (my $i = 0; $i < scalar(@$cigarContentRef); $i += 1)	{
		my $newToken = $cigarContentRef->[$i]{token};
		$newToken =~ tr/ID/DI/;
		$cigarContentRef->[$i]{token} = $newToken;
	} # for i
	if ($orientation =~ /^-$/)	{
		for (my $i = $#{$cigarContentRef}; $i >= 0; $i -= 1)	{
			push(@newCigarContent, {start => $cigarContentRef->[$i]{start}, end => $cigarContentRef->[$i]{end}, token => $cigarContentRef->[$i]{token}});
		} # for i
	} else	{
		for (my $i = 0; $i <= $#{$cigarContentRef}; $i += 1)	{
			push(@newCigarContent, {start => $cigarContentRef->[$i]{start}, end => $cigarContentRef->[$i]{end}, token => $cigarContentRef->[$i]{token}});
		} # for i
	} # orientation is negative

	my $newCigarString = "";
	for (my $i = 0; $i < scalar(@newCigarContent); $i += 1)	{
		$newCigarString .= $newCigarContent[$i]{token};
	} # for i

	return $newCigarString;
} # swapCigar

sub findCigarMatchPositions   {
	my ($regex, $string) = @_;
	my @locations = ();
	while ($string =~ /($regex)/g)  {
		push(@locations, {start => pos($string) - length($1), end => pos($string), token => $1});
	} # while string
	return \@locations;
} # findCigarMatchPositions

sub swapPositions	{
	my ($qryStart, $qryEnd, $refStart, $refEnd, $orientation) = @_;
	($qryStart, $refStart) = swap($qryStart, $refStart);
	($qryEnd, $refEnd) = swap($qryEnd, $refEnd);
	# reference start should always be less than reference end
	if ($orientation =~ /^-$/)	{
		# negative orientation alignment, needs to swap the start and end co-ordinates
		($qryStart, $qryEnd) = swap($qryStart, $qryEnd);
		($refStart, $refEnd) = swap($refStart, $refEnd);
	} # if orientation
	return ($qryStart, $qryEnd, $refStart, $refEnd, $orientation);
} # swapPositions

sub swap	{
	my ($data1, $data2) = @_;
	return ($data2, $data1);
} # swap

sub readAlignFinalXmap	{
	my ($file, $headersRef, $contentRef) = @_;
	
	open(IN, $file) or dieLog("ERROR: readAlignFinalXmap: cannot open $file: $!\n");
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;

		if ($line =~ /^#/)	{
			# header lines
			if ($line =~ /# XMAP File Version:/ || $line =~ /# Label Channels:/ || $line =~ /# Reference Maps From:/ || $line =~ /# Query Maps From:/ || $line =~ /#h/ || $line =~ /#f/)	{
				push(@$headersRef, $line);	
			} # if essential header line
			next;
		} # if header line
		
		# data line
		my @content = split(/\t/, $line);
		my ($xId, $qryId, $refId, $qryStart, $qryEnd, $refStart, $refEnd, $orientation, $confidence, $cigarString) = @content[0..9];
		my ($qryLength, $refLength, $labelChannel, $alignmentString) = (-1.0, -1.0, -1, "");
		if (scalar(@content) >= 14)	{
			($qryLength, $refLength, $labelChannel, $alignmentString) = @content[10..13];
		} # if newer xmap format		
		my ($mapWt) = (-1.0);
		if (scalar(@content) >= 15)	{
			($mapWt) = ($content[14]);
		} # if newer xmap format
	
		push(@$contentRef, {xId => $xId, qryId => $qryId, refId => $refId, qryStart => $qryStart, qryEnd => $qryEnd, refStart => $refStart, refEnd => $refEnd, orientation => $orientation, confidence => $confidence, cigarString => $cigarString, qryLength => $qryLength, refLength => $refLength, labelChannel => $labelChannel, alignmentString => $alignmentString, mapWt => $mapWt});	
	} # while line

	close IN;
	return ($headersRef, $contentRef);
} # readAlignFinalXmap

sub runAlignment	{
	my ($refAligner, $qryFile, $refFile, $outFilePrefix, $baseStagesRef, $logFileHandle) = @_;
	
	my $XML = new XML::Simple(KeyAttr=>[]);
	my $configRef = $XML->XMLin($xmlFile);
	my %stageStack = ();
	my @commandStack = ();
	my $commandStackRef = \@commandStack;

	push(@$commandStackRef, $refAligner);
	# output
	push(@$commandStackRef, ("-o", "$outFilePrefix"));
	# stdout stderr
	push(@$commandStackRef, ("-stdout", "-stderr"));
	# -i, -ref
	push(@$commandStackRef, ("-i", $qryFile, "-ref", $refFile));
	
	for (my $i = 0; $i < scalar(@$baseStagesRef); $i += 1)	{
		# maxmem maxthreads etc and other parameters
		my $baseStage = $baseStagesRef->[$i];
		$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);
	} # for i

	# command call
	print "\nRunning commmand: ".(join(" ", @$commandStackRef))."\n";
	my ($outResults0, $errResults0, $jobStatus0) = runCommand($commandStackRef);
	if ($jobStatus0 != 0)	{
		# print OUT ERROR
		$errResults0 .= "\nERROR: command ".(join(" ", @$commandStackRef))."\nerror was encountered with exit code $jobStatus0.";
		logMessage($logFileHandle, $errResults0);
		close($logFileHandle);
		dieLog("$errResults0\n");
	} # if jobStatus
	@commandStack = ();
} # runAlignment

sub findLeftOverSequences	{
	# this subroutine parses the 1st_pass xmap file, determine which sequence (reference) was aligned
	# compares with the original sequence (post-conflict resolved)
	# finds the sequences ids that were left out, writes them to a text file
	# call RefAligner to extract the left out seq ids from the original sequence
	my ($seqFile, $notUsedSeqCmapFilePrefix, $usedSeqRef, $refAligner, $logFileHandle) = @_;

	my $notUsedSeqRef = getNotUsedSeq($seqFile, $usedSeqRef);
	my $notUsedSeqIdFile = $notUsedSeqCmapFilePrefix."_id.txt";
	printIds($notUsedSeqIdFile, $notUsedSeqRef);

	# call RefAligner
	my ($outResults, $errResults, $jobStatus) = ("", "", 0);
	if (scalar(keys %$notUsedSeqRef) > 0)	{
		my $command = "$refAligner -f -merge -selectidf $notUsedSeqIdFile -i $seqFile -o $notUsedSeqCmapFilePrefix -minsites 0 -stdout -stderr";
		my @commandContent = split(/\s+/, $command);
		my ($outResults, $errResults, $jobStatus) = runCommand(\@commandContent);
		if ($jobStatus != 0)	{
			# print OUT ERROR
			$errResults = "ERROR: findLeftOverSequences ".(join(" ", @commandContent))."\nerror was encountered with exit code $jobStatus.";
			close($logFileHandle);
			dieLog("$errResults\n");
		} # if jobStatus
	} # if scalar
	
} # findLeftOverSequences

sub printIds	{
	my ($file, $idsRef) = @_;
	open(OUT, ">$file") or dieLog("ERROR: printIds: cannot write to $file: $!\n");
	foreach my $id (sort numeric keys %$idsRef)	{
		print OUT "$id\n";
	} # foreach id
	close OUT;
} # printIds

sub numeric	{	$a	<=>	$b	}

sub getNotUsedSeq	{
	my ($file, $usedSeqRef) = @_;
	my %notUsedSeq = ();

	open(IN, $file) or dieLog ("ERROR: getNotUsedSeq: cannot open original (post-conflict resolved) sequence file: $!\n");
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;
		next if ($line =~ /^#/);

		# data line, ehck the id and see if that ocntig has been aligned
		my @content = split(/\t/, $line);
		my $id = $content[0];
		if (! exists $usedSeqRef->{$id})	{
			$notUsedSeq{$id} = 1;
		} # if not exists	
	} # while line
	close IN;
	return \%notUsedSeq;
} # getNotUsedSeq

sub getUsedSeqScaffold	{
	my ($alignXmapFile, $usedSeqRef, $usedScaffoldRef) = @_;
	open(IN, $alignXmapFile) or dieLog("ERROR: getUsedSeqScaffold: cannot open alignment XMAP file $alignXmapFile: $!\n");
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;

		next if ($line =~ /^#/);

		# data lines
		# sequence is column 3
		my @content = split(/\t/, $line);
		my ($qryId, $refId) = ($content[1], $content[2]);
		$usedSeqRef->{$refId} = 1;
		$usedScaffoldRef->{$qryId} = 1;
	} # while line
	close IN;
	return ($usedSeqRef, $usedScaffoldRef);
} # getUsedSeqScaffold

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
	my $childExitStatus = $? >> 8;
	return ($outResults, $errResults, $childExitStatus);
} # runCommand

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
