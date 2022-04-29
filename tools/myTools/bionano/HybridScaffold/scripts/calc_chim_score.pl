# $Id: calc_chim_score.pl 7281 2018-01-02 02:59:16Z apang $

#!/usr/bin/perl -w

# 1) needs hybrid scaffold pipine xml, needs auto noise errbin file, needs molecule bnx, needs the refinefinal cmap file
# 2) parse the xml to determine parameter
# 3) run RefAligner with refine 1 and chim quality flag
# 4) output the cmap file with chimeric quality score for each label

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
        #print "$0 library paths:\n\t"; print join("\n\t", @INC); print "\n";
}

use BNG::Utility;
use File::Copy qw(copy move);
use Getopt::Long;
use IPC::Open3;
use IO::Select;
use XML::Simple;
print "\nInfo: Running the command: $0 @ARGV\n";

my $refAligner = "";
my ($inMapDir, $inMapFile) = ("", "");
my ($outMapDir, $outMapFile) = ("", "");
my ($xmlDir, $xmlFile) = ("", "");
my ($noiseParamDir, $noiseParamFile) = ("", "");
my ($bnxDir, $bnxFile) = ("", "");
my $maxmem = 128;
my $helpFlag = "";

if (scalar(@ARGV) == 0)	{
	usage();
} # if scalar

GetOptions	(
	"help"				=>	\$helpFlag,
	"refAligner=s"			=>	\$refAligner,
	"inMapFile=s"			=>	\$inMapFile,
	"outMapFile=s"			=>	\$outMapFile,
	"xmlFile=s"			=>	\$xmlFile,
	"noiseParamFile=s"		=>	\$noiseParamFile,
	"bnxFile=s"			=>	\$bnxFile,
	"maxmem=i"			=>	\$maxmem
) or dieLog ("ERROR: calc_chim_score: error in command line arguments.\n ");

if ($helpFlag !~ /^$/)	{
	usage();
} # if helpFlag


($inMapDir, $inMapFile) = separateDirFile(abs_path(glob($inMapFile)));
($outMapDir, $outMapFile) = separateDirFile(abs_path(glob($outMapFile)));
($xmlDir, $xmlFile) = separateDirFile(abs_path(glob($xmlFile)));
($noiseParamDir, $noiseParamFile) = separateDirFile(abs_path(glob($noiseParamFile)));
($bnxDir, $bnxFile) = separateDirFile(abs_path(glob($bnxFile)));

# quick checks on the existence of the files
dieLog ("ERROR: chim_qual_calc: RefAligner $refAligner does not exist\n") if (! -e $refAligner);
dieLog ("ERROR: chim_qual_calc: genome map file $inMapDir/$inMapFile does not exist\n") if (! -e "$inMapDir/$inMapFile");
dieLog ("ERROR: chim_qual_calc: XML file $xmlDir/$xmlFile does not exist\n") if (! -e "$xmlDir/$xmlFile");
dieLog ("ERROR: chim_qual_calc: molecule bnx file $bnxDir/$bnxFile does not exist\n") if (! -e "$bnxDir/$bnxFile");
dieLog ("ERROR: chim_qual_calc: maximum memory specified=$maxmem is not valid\n") if ($maxmem !~ /\d+/);

mkdir $outMapDir if (! -e $outMapDir);

# parse the XML file
my $XML = new XML::Simple(KeyAttr=>[]);
my $configRef = $XML->XMLin("$xmlDir/$xmlFile");

my %stageStack = ();
my @commandStack = ();	
my $commandStackRef = \@commandStack;


# refAligner 
push(@$commandStackRef, $refAligner);
# output 
my $outMapFilePrefix = $outMapFile;	$outMapFilePrefix =~ s/\.cmap$//;
push(@$commandStackRef, "-o");
push(@$commandStackRef, "$outMapDir/$outMapFilePrefix");

# the -ref
push(@$commandStackRef, "-ref");
push(@$commandStackRef, "$inMapDir/$inMapFile");

# stdout stderr
push(@$commandStackRef, "-stdout");
push(@$commandStackRef, "-stderr");

# set the global variables
my $baseStage = "global";
$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);

# parse the refineFinal command from xml file and store the data in @$commandStackRef array
$baseStage = "refineFinal1";
$commandStackRef = parseConfig($configRef, $baseStage, $commandStackRef);

my $commandStackString = join(" ", @$commandStackRef);
$commandStackString .= " ";

# now change refine 3 to refine 1 - assume there must be refine 3
$commandStackString =~ s/-refine 3/-refine 1/;

# make sure we set EndTrim to 0 so that we don't change the genome map
if ($commandStackString =~ /-EndTrim \S+/)	{
	$commandStackString =~ s/-EndTrim \S+/-EndTrim 0.0/g;
} else	{
	$commandStackString .= "-EndTrim 0.0 ";
} # commandStackString EndTrim

# change maxmem
if ($commandStackString =~ /-maxmem \d+/)	{
	$commandStackString =~ s/-maxmem \d+/-maxmem $maxmem/;
} else	{
	$commandStackString .= "-maxmem $maxmem ";
} # commandStackString maxmem

# add chim quality score
if ($commandStackString !~ /-ChimQuality/)	{
	$commandStackString .= "-ChimQuality ";
} # if commandStackString ChimQuality

$commandStackString =~ s/\s+$//;
my @commandContent = split(/\s/, $commandStackString);

# now add the noise parameters
if ($noiseParamFile =~ /\.errbin$/ && -e "$noiseParamDir/$noiseParamFile" && -B "$noiseParamDir/$noiseParamFile")	{
	# an errbin binary file
	print "\nFor noise parameter: an .errbin file was provided\n";
	push(@commandContent, "-readparameters");
	push(@commandContent, "$noiseParamDir/$noiseParamFile");
} elsif ($noiseParamFile =~ /\.err$/ && -e "$noiseParamDir/$noiseParamFile" && -T "$noiseParamDir/$noiseParamFile")	{
	# not an errbin file, but the .err text file
	print "\nFor noise parameter: an .err file was provided\n";
	my ($errFileHeaderRef, $errFileContentRef) = parseErrFile($noiseParamDir, $noiseParamFile);
	# in case the .err file is very old
	push(@commandContent, "-FP");		if (exists $errFileHeaderRef->{FP})	{	push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{FP}]);	}	else	{	push(@commandContent, 1.5);	}
	push(@commandContent, "-FN");
	if (exists $errFileHeaderRef->{FN})	{
		push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{FN}]);
	} elsif (exists $errFileHeaderRef->{FNrate})	{
		push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{FNrate}]);
	} else	{
		push(@commandContent, 0.15);
	} # if FN
	push(@commandContent, "-sd");
	if (exists $errFileHeaderRef->{sd})	{
		push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{sd}]);
	} elsif (exists $errFileHeaderRef->{ScalingSD})	{
		push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{ScalingSD}]);
	} else	{
		push(@commandContent, 0.0);	
	} # if sd
	push(@commandContent, "-sf");	
	if (exists $errFileHeaderRef->{sf})	{
		push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{sf}]);
	} elsif (exists $errFileHeaderRef->{SiteSD})	{
		push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{SiteSD}]);
	} else	{
		push(@commandContent, 0.2);
	} # if sf
	push(@commandContent, "-sr");		if (exists $errFileHeaderRef->{sr}) {	push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{sr}]);	} else	{	push(@commandContent, 0.03);	}
	push(@commandContent, "-res");		if (exists $errFileHeaderRef->{res})	{	push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{res}]);	}	else	{	push(@commandContent, 3.3);	}
	push(@commandContent, "-bpp");		if (exists $errFileHeaderRef->{bpp})	{	push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{bpp}]);	}	else	{	push(@commandContent, 500);	}
	push(@commandContent, "-se");		if (exists $errFileHeaderRef->{se}) {	push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{se}]);	} else	{	push(@commandContent, 0.0);	}
	push(@commandContent, "-resSD");	if (exists $errFileHeaderRef->{resSD})	{	push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{resSD}]);	}	else	{	push(@commandContent, 0.75);	}
	push(@commandContent, "-mres");		if (exists $errFileHeaderRef->{mres})	{	push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{mres}]);	}	else	{	push(@commandContent, 0.9);	}
	push(@commandContent, "-mresSD");	if (exists $errFileHeaderRef->{mresSD})	{	push(@commandContent, $errFileContentRef->[$errFileHeaderRef->{mresSD}]);	}	else	{	push(@commandContent, 0.0);	}
} else	{
	# does not have any noise parameter file, just set default
	# -FP 1.5 -FN 0.15 -sd 0.0 -sf 0.2 -sr 0.03 -res 3.3 -bpp 500 -se 0.0 -resSD 0.75 -mres 0.9 -mresSD 0.0
	print "\nFor noise parameter: no noise parameter file was provided, thus using default -FP 1.5 -FN 0.15 -sd 0.0 -sf 0.2 -sr 0.03 -res 3.3 -bpp 500 -se 0.0 -resSD 0.75 -mres 0.9 -mresSD 0.0\n";
	push(@commandContent, "-FP");		push(@commandContent, 1.5);
	push(@commandContent, "-FN");		push(@commandContent, 0.15);
	push(@commandContent, "-sd");		push(@commandContent, 0.0);
	push(@commandContent, "-sf");		push(@commandContent, 0.2);
	push(@commandContent, "-sr");		push(@commandContent, 0.03);
	push(@commandContent, "-res");		push(@commandContent, 3.3);
	push(@commandContent, "-bpp");		push(@commandContent, 500);
	push(@commandContent, "-se");		push(@commandContent, 0.0);	
	push(@commandContent, "-resSD");	push(@commandContent, 0.75);
	push(@commandContent, "-mres");		push(@commandContent, 0.9);
	push(@commandContent, "-mresSD");	push(@commandContent, 0.0);
} # if errBinFile

# finally the input bnx file
push(@commandContent, "-i");
push(@commandContent, "$bnxDir/$bnxFile");

# now do the RefAligner call
print "\nRunning command: ".(join(" ", @commandContent))."\n";
my ($outResults, $errResults, $jobStatus) = ("", "", 0);
($outResults, $errResults, $jobStatus) = runCommand(\@commandContent);
if ($jobStatus != 0)	{
	# print out error
	$errResults = "ERROR: in generating chimeric quality score ".(join(" ", @commandContent))." with exit code $jobStatus; out info: $outResults; error info: $errResults";
	dieLog("$errResults\n");
} # if jobStatus
@commandContent = ();

# look at the output directory and look for $outMapFilePrefix_contig\d+.cmap files
getFileList($outMapDir, $outMapFilePrefix."_contig", $outMapDir, "file_list.txt");
# call RefAligner to merge those individual contig cmap files into one
@commandContent = ($refAligner, "-merge", "-f", "-if", "$outMapDir/file_list.txt", "-o", "$outMapDir/$outMapFilePrefix");
print "\nRunning command: ".(join(" ", @commandContent))."\n";
($outResults, $errResults, $jobStatus) = runCommand(\@commandContent);
if ($jobStatus != 0)	{
	# print out error
	$errResults = "ERROR: in generating chimeric quality score ".(join(" ", @commandContent))." with exit code $jobStatus; out info: $outResults; error info: $errResults";
	dieLog("$errResults\n");
} # if jobStatus
@commandContent = ();



## copy from _r.cmap to .cmap, as the refinefinal command will generate a genome map with chim score in _r.cmap
#copy ("$outMapDir/$outMapFilePrefix"."_r.cmap", "$outMapDir/$outMapFile") or dieLog "ERROR: calc_chim_score: cannot perform copy $outMapDir/$outMapFilePrefix"."_r.cmap to $outMapDir/$outMapFile: $!\n";

#print Dumper($config);

sub getFileList	{
	# this subroutine simply scales $dir for files with prefix $filePrefix
	# then it prints that list in to $outDir/$outFile
	my ($dir, $filePrefix, $outDir, $outFile) = @_;
	opendir(DIR, $dir) or dieLog("ERROR: getFileList: cannot open dir $dir: $!\n");
	my @files = grep {$_ =~ /^$filePrefix/} readdir DIR;
	closedir DIR;
	
	open(OUT, ">$outDir/$outFile") or dieLog("ERROR: getfileList: cannot write to $outDir/$outFile: $!\n");
	foreach my $file (sort @files)	{
		print OUT "$outDir/$file\n";
	} # foreach file
	close OUT;
} # getFileList

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

sub parseErrFile	{
	my ($dir, $file) = @_;
	my %header = ();	# stores the column number as value and the column name (e.g. FP) as key
	my @content = ();	# stores the last line value
	open(IN, "$dir/$file") or dieLog "ERROR: parseErrFile: cannot open $dir/$file: $!\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		# skip header lines
		next if ($line =~ /^#/);
		$line =~ s/^\s+//;	$line =~ s/\s+$//;	$line =~ s/\s+/\t/g;
		if ($line =~ /^Iteration/i)	{
			# the line contains the label of the columns
			my @headerContent = split(/\t/, $line);
			# push the label to %header 
			for (my $i = 0; $i < @headerContent; $i += 1)	{
				# strip off the trailing brackets and contents
				$headerContent[$i] =~ s/\(.+\)//;	# for here, don't need to differentiate Log10LR(/Maps) and log10LR(/GoodMaps)
				$header{$headerContent[$i]} = $i;
			} # for i
		} else	{
			# data line
			@content = split(/\t/, $line);
		} # if line	
	} # while line
	close IN;
	return (\%header, \@content);
} # parseErrFile

sub separateDirFile	{
	my ($fullFile) = @_;
	my ($dir, $file) = ("" ,"");
	if ($fullFile =~ /^$/)	{
		return ($dir, $file);
	} # if fullFile
	my $tempFile = $fullFile;
	my @tempContent = split(/\//, $tempFile);
	$file = $tempContent[$#tempContent];
	
	if ($fullFile =~ /^\// && scalar(@tempContent) == 2)	{
		# the file is in the "/" the very top directory
		# tempContent[0] is "", tempContent[1] is "$file"
		return ("/", $file);
	} # if fullFile
	$dir = $tempFile;
	$dir =~ s/$file$//;	
	$dir =~ s/\/$//;	# remove trailing /
	# if only a file with no path has been specified (meaning the file is in the current directory)
	# fullFile did not start with "/" either
	$dir = ($dir =~ /^$/) ? (".") : ($dir);	
	return ($dir, $file);
} # separateDirFile

sub usage	{
	print << "EOF";

Usage: perl calc_chim_score.pl <-help> <-refAligner RefAligner> <-inMapFile in_cmap_file> <-outMapFile out_cmap_file> <-xmlFile xml_file> <-noiseParamFile noise_parameter_file> <-bnxFile bnx_file> <maxmem maximum_memory>
	-help			: This help message
	-refAligner		: RefAligner program [required]
	-inMapFile		: Input genome map CMAP file to generate chimeric quality scores for [required]
	-outMapFile		: Output genome map CMAP file with chimeric quality scores [required]
	-xmlFile		: Input Hybrid Scaffold optArguments XML file [required]
	-noiseParamFile		: The .errbin or .err file generated during auto noise estimation in the de novo assembly [optional]
	-bnxFile		: The molecule BNX file used to build the input genome map in the de novo assembly [required]
	-maxmem			: Maximum memory allowed [optional]
EOF
	exit; 
} # usage
