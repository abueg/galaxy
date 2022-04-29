# $Id: align_final_bng.pl 7271 2017-12-28 04:08:17Z apang $

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
my $bngFile = "";
my $scaffoldFile = "";
my $logFile = "";

GetOptions	(
	"outDir=s"		=>	\$outDir,
	"outFilePrefix=s"	=>	\$outFilePrefix,
	"refAligner=s"		=>	\$refAligner,
	"xmlFile=s"		=>	\$xmlFile,
	"bngFile=s"		=>	\$bngFile,
	"scaffoldFile=s"	=>	\$scaffoldFile,
	"logFile=s"		=>	\$logFile
) or dieLog("ERROR: align_final_bng: erro in command arguments.\n");

# final CMAP and XMAjP file names
mkdir $outDir if (! -e $outDir);
my $outBngCmapFile = "$outDir/$outFilePrefix"."_q.cmap";
my $outHybridCmapFile = "$outDir/$outFilePrefix"."_r.cmap";
my $outXmapFile = "$outDir/$outFilePrefix.xmap";

open(LOG_FILE, ">$logFile") or dieLog("ERROR: align_final_bng: cannot open log file $logFile: $!\n");
my $logFileHandle = *LOG_FILE;

### align the input Bionano genome map (post-conflict resolved) to the hybrid scaffold ###
my @baseStages = ("global", "align_final_BNG");
my $baseStagesRef = \@baseStages;
logMessage($logFileHandle, "Running alignment of bionano genome map to hybrid scaffold...\n");
runAlignment($refAligner, $bngFile, $scaffoldFile, "$outDir/$outFilePrefix", $baseStagesRef, $logFileHandle);
logMessage($logFileHandle, "alignment of Bionano genome map to hybrid scaffold completed successfully\n");

close($logFile);
print "$0 finished successfully.\n\n";

### subroutines ###
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

sub numeric	{	$a	<=>	$b	}

