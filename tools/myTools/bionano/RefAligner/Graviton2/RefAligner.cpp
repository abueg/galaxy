#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <math.h>
#include <ctype.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <fcntl.h>
#include <malloc.h>

#ifndef WIN32
#include <unistd.h>
#endif

#include <sys/stat.h>

using namespace std;


#include "globals.h"
#include "parameters.h"
#include "constants.h"
#include "hash.h"
#include "Ccontig.h"

#include "timers.h"

#define FORK_DEBUG 0 // to debug fork() & waitpid (currently does NOT work on MIC)

#if FORK_DEBUG
#include <sys/wait.h>
#endif

#if defined(_MSC_VER)
#define strtoll _strtoi64
#endif

#ifndef WIN32
#include <sys/types.h>

#endif


/* TODO :
  A1. Try -skip_dist 0.001, 0.01, 0.02, 0.05 0.1 in refineFinal to see if the speedup is worth the degradation in SV ROC 

  A2. Avoid spurious endoutliers during Extend & Split in high label density regions:
     A. Try MIS_VITERBI = 0 in refalign.cpp to raise scores of regions with mis-resolved labels
     B. Try -biaswtOutlier 1.0 (or value > 0.0)
     C. Try adding correct to -biaswt for misresolved regions : current misresolved regions are excessively penalized by -biaswt since the biaswt applied expected penalty does not reflect misresolved sites.
        Expected misresolve penalty for each contig interval Y[i..i+1] == Yi is approximately equal to Pr(Yi)log(Pr(Yi)) + (1-Pr(Yi)log(1-Pr(Yi))

  B. Figure out why refiner cannot fix consensus sizing using medians when false insertion and false deletions are present and all molecules crossing them have outliers.
      (contig 170 and contig907 in Final Refine of /mnt/bionf_tmp/apang/IPS/sample_3/human_ips3_49x_run1/contigs, and contig244,322,846,312 in .../sample_1/human_ips1_94x_run1/contigs/)
  1a. test -TrimNorm on DB assembly and optimize -CovTrim (1,2,3) -CovTrimLen (45,55,65kb) and -contigsplit 0.2 0.35(0.25,0.3,0.35) ID and -splitrev (0,1,2) AND -TrimNorm (0,1,2,3,4,5) AND -TrimNormChim (1,2,3,4,5)
      NOTE: CovTrimLen should be specified as multiple of Xlambda instead of absolute kb, and TRIMNORM_FIX = 0,1. Also -splitsite and -TrimOutlier
  1b. Exclude alignments with internal outliers in count of confirming molecules.
  1c. Add option -CovTrimCenter <width> to specify a width over which long molecules must be aligned along with -CovTrim & CovTrimLen on each side, otherwise a chimeric merge due to 2 paralog regions of size <width> can be assumed..
  1d. Test -TrimNormMin <2,3>
  1e. Test -TrimNormEnd <2,3,4,5>

  2. Test -PVres 2, -PVendoutlier and -AlignRes
   
  3. Try -endoutlier 1e-8 -outlier 1e-8 during error parameters estimation, to avoid classifying large sizing errors as outliers. Check if this would also help during extension (and possibly refinement B)
     Conversely when refining reference (for SV correction etc) may need to use use -outlier 1e-3 -endoutlier 1e-3 to encourage more molecules to align initially, but then use -endoutlierRef in reverse order to use more stringent pvalues during refinement (possibly with OUTLIER_REFINE==2 in refine.h)

  4. Add adaptive scaling option for refinement
  
  5. Add E/M estimation of res & resSD (use seperat list of site "triangles" to supplement the errors[] list of aligned intervals) 
 */

static int GenerateFasta = 0;// see option -merge -fasta
static long long IncID = 0;// see option -incID

const char *SVN_ID =  (char *)"$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/RefAligner.cpp 11543 2020-08-26 18:53:06Z tanantharaman $";

char *usage1 = (char *)"RefAligner \n"
  "-help : print this help message\n"
  "                      NOTE: Except for -i, -if, -ref and -reff, only the last instance of an option on the commandline has any effect\n"
  "-i [input.vfixx | input.cmap | input.bnx] : input maps from .vfixx, .cmap or .bnx file (Use more than one -i option to input multiple input files)\n"
  "-if input.filenames : input maps from all files listed in input.filenames (one per line), as if they had been input with multiple -i options\n"
  "                      If file names in input.filenames are relative (not starting with /) it is assumed to be relative to the directory in which input.filenames is located\n"
  "                      Use multiple -if options to input multiple sets of input files. \n"
  "  -bnxerror : If input map from .bnx file has sites beyond the end of the map, exit with an error message (Default is print WARNING and extend map ends)\n"
  "  -SNR : Require input .bnx file to have SNR information (one line per color) [Default : OFF]\n"
  "  -Intensity : Require input .bnx file to have Intensity information (one line per color) [Default : OFF]\n"
  "  -stitch : Require input .bnx file to have image stitch information (one line per color, starting with QX13 for BNX 1.0) [Default : OFF]\n"
  "  -stitchLoc : Require input .bnx file to have image stitch location information (one line per color, starting with QX14 for BNX 1.0) [Default : OFF]\n"
  "  -PSFWidth : Require input .bnx file to have PSF width information (one line per color, starting with QX15 for BNX 1.0) [Default : OFF]\n"
  "  -QXfilter : Delete all QX information except as specified (by -SNR -stitch -PSFWidth etc) from the input BNX [Default : OFF]\n"
  "                                        NOTE : Image Location information (QX16/QX26) is NOT deleted by -QXfilter\n"
  "  -QXerror <0/1> : If != 0 : exit with error if different BNX input files have different sets of QX codes [Default 0 : Discard non-common QX codes]\n"
  "  -bnxmergeerror <N> : If <N> = 1 : Don't allow BNX 1.x errors like missing Run Data or duplicate Run Data lines (<N> = 0 : just print warning) [Default 0]\n"
  "  -bnxBigid <N> : If <N> = 1 : Use 19 digit id encoding, if possible (If not possible due to field width overflows, use sequential ids starting at 1)\n"
  "                  If <N> = 0 : Avoid changing molecule ids (If duplicate ids in input(s), use sequential ids starting at 1) [Default 0]\n"
  "                  If any ids are changed, output <prefix>_renumbered.bnx showing the new and old ids\n"
  "  -bnxversion <V> : Specify BNX version (0.x or 1.x) required for input and used for output [Default : based on bnx input OR Version 0.1 if no bnx input]\n"
  "-ref [input.spots | input.cmap] : input reference map from .spots or .cmap file (Use more than one -ref option to merge multiple reference files)\n"
  "-reff ref.filenames : input reference maps from all files listed in ref.filesnames (one per line), as if they had been input as multiple -ref options\n"
  "                      If file names in ref.filenames are relative (not starting with /) it is assumed to be relative to the directory in which ref.filenames is located\n"
  "                      Use multiple -reff options to input multiple sets of input files\n"
  "-NFSdelay <N> : Maximum delay in seconds between file write (fclose) on one NFS client and file becoming visible on another NFS client (should match acdirmax + 1) [Default: 1]\n"
  "-Tracksites : track original index of labels in _q.cmap and _r.cmap back to input map file.\n"
  "              Label index values for 2-color data are interleaved index values (as used in CMAP files)\n"
  "-refmap <.map file> : override startloc/endloc of input.vfixx based on .map file\n"
  "-o [output prefix] : use specified prefix name for all output file names [It is an error not to specify this option. Use -o /dev/null to skip all file output]\n"
  "-output-veto-filter regex : regular expression to specify which files to not write. [Default None]\n"
  "-output-filter regex : regular expression to specify which files to write. [Default .*]\n"
  "-stdout [<file.stdout>] : Redirect stdout to file.stdout. Does NOT apply to error messages from this and any previous options on the command line. <file.stdout> defaults to <prefix>_id<ID>.stdout, if -o -id were specified previously\n"
  "-stderr [<file.stderr>] : Redirect stderr to file.stderr. Does NOT apply to error messages from this any any previous options on the command line. <file.stderr> defaults to <prefix>_id<ID>.stdout, if -o -id were specified previously\n"
  "-f | -force : force output file overwrites [Default: exit with error message if output file already exists]\n"
  "-colors <N> : <N> is the number of colors required in the input data : currently any value other than 1 or 2 is not supported. [Default : 1 or value in first -i or -readparameters input file header]\n"
  "-usecolor <N> : If -i input data has 2 colors only the <N> th color is used as if only that one color were present (but output BNX or _q.cmap will include both colors).\n"
  "-hashcolor <N> : If -i input data has 2 colors only the <N> th color is used by the hashtable. Should specify color with best site density. NOT required with -usecolor. [Default : OFF or 0]\n"
  "-bnx : output -merge output as .bnx file (instead of .cmap file). Also output _q.bnx instead of _q.cmap for .xmap output\n"
  "-minlen <L> : If L != 0 filter out all -i input maps that are smaller than L kb (But see -ScanScaling, which defers application of -minlen to until rescaled BNX output) [Default 0]\n"
  "-maxlen <L> : If L != 0 filter out all -i input maps that are larger than L kb (But see -ScanScaling, which defers application of -maxlen to until rescaled BNX output) [Default 0]\n"
  "-bppadjust <0/1> : If != 0 : Apply -minlen and -maxlen AFTER performing -bpp adjustment (instead of BEFORE) [Default 1] (Obsolete)\n"
  "-minsites <N> : filter out all -i input maps with fewer than N sites [Default 1]\n"
  "-maxsites <N> [ <M> ]: If N != 0 filter out all -i input maps with more than N sites [Default: OFF or 0]\n"
  "                NOTE : If -colors 2 and -usecolor is NOT applied, -minsites and -maxsites are based on the total number of sites on both colors\n"
  "                NOTE : If -minSNRestimate is specified, the N based filter is only applied on output of <prefix>_rescaled.bnx. In that case another larger value M can be applied right after -i input\n"
  "-maxSiteDensity <D> : If D != 0 : filter out all -i input maps with more than D sites per 100 kb [Default: OFF or 0]\n"
  "                NOTE : If -minSNRestimate is specified, this filter is only applied on output of <prefix>_rescaled.bnx.\n"
  "-minSiteDensity <D> : If D != 0 : filter out all -i input maps with less than D sites per 100kb [Default: OFF or 999]\n"
  "                NOTE : If -minSNRestimate is specified, this filter is only applied on output of <prefix>_rescaled.bnx.\n"
  "-maxEndSiteDensity <D> <L> : If D != 0 : filter out all -i input maps with more than D sites per 100kb for any end region of L kb or more [Default: OFF or 0]\n"
  "                NOTE : This filter is only applied on output of <prefix>_rescaled.bnx by autoNoise.\n"
  "-maxEnd <L> : Trim unlabeled molecule ends to no more than <L> kb [Default OFF]\n"
  "-minEnd <L> : Extend unlabeled molecule ends to at least <L> kb. Also applies to map fragments generated by -break or -contigsplit [Default 0.020 kb]\n"
  "-maxInterval <L> [ <ref> ]: If label interval is larger than <L> break molecule into 2 pieces. Applies only to BNX input unless -ref or -reff inputs are absent [Default OFF or -1 0]\n"
  "                            If <ref> is specified as 1, also applies to -ref CMAP input. CMAP maps are split into multiple maps. For BNX maps only largest piece is retained.\n"
  "-breakChiSq <PV> <F> [ <L> <Cmax> <S> <R> <I> <SR> <Cmin> <CR> <Ctrim> <N> <Lmax> <Fr> <Sd> <Ch>]: Break -i CMAP input map at intervals above I kb and below L kb IF any of the following applies:\n"
  "                        OutlierFrac value is above F (%) \n"
  "                     OR MolCov is below Cmax AND ChiSq score below PV AND MolSd above ExpSd + S kb + R times interval size AND MolSd above ExpSd * SR\n"
  "                     OR MolCov is below both Cmax and Cmin + CR * RangeRatio (ratio of larger map range to average label interval size) \n"
  "                            AND MolCov rises by at least factor Ctrim in both directions.\n"
  "                     OR FragileLeft+FragileRight is above Fr%\n"
  "                     OR SegDupLeft+SegDupRight is above Sd%\n"
  "                     If N > 0 : Instead of breaking the map, flag the qry interval involved (+ (N-1) neighboring intervals on each side) and suppress all deletion calls below Lmax kb\n"
  "                           that spans any flagged qry intervals AND ignore Ctrim (force its value to 1.0)\n"
  "                     If Ch > 0 : Replace MolCov with min(MolCov, ChimNorm * ChimQuality * Ch / 100)\n"
  "                     Requires that CMAP was previously refined with -OutlierQuality 2 to generate OutlierFrac and MolSd, ExpSd, MolCov and ChiSq values for each interval.\n"
  "                     NOTE: Requires MolSd etc in CMAP input (see -OutlierQuality): Can be used to ignore outliers caused by DNA knots in the molecules. [Default OFF or 0 100]\n"
  "-MaxIntensity <V> : If V != 0 : Filter out all -i input maps with AvgIntensity above V [Default 0]\n"
  "-minSNR <V> ... : Remove all -i input map sites with SNR below V (There is one V value per color) [Default 0.0]\n"
  "-maxPSFWidth <W> .. : Remove molecules with any label with PSFWidth larger than <W> bp (There is one W value per color). Ignored if no information on PSFWidth (QX15 or QX25) is present in BNX input. [Default OFF]\n"
  "-skipidf <file> : skip maps with ids that are in a file (one id per line)\n"
  "-selectid <ID1> ... <IDn> : select a subset of -i input maps based on the specfied set of map IDs specified\n"
  "-selectidf <file> : same as -selectid but the list of IDs are in a file (one id per line)\n"
  "-selectScanId <ID1> <IDn> : select a subset of the -i input maps based on matching a specific GlobalScanNumber value range ID1 .. IDn (NOTE: GlobalScanNumber is same as RunId in Saphyr BNX molecule header)\n"
  "-sort-sizeinc : Sort -i input maps in order of increasing size [Default OFF]\n"
  "-sort-sizedec : Sort -i input maps in order of decreasing size [Default OFF]\n"
  "-sort-sitesinc : Sort -i input maps in order of increasing number of sites [Default OFF]\n"
  "-sort-sitesdec : Sort -i input maps in order of decreasing number of sites [Default OFF]\n"
  "-sort-idinc : Sort -i input maps in order of increasing map id [Default OFF]\n"
  "-sort-runindexinc : Sort -i input maps in order of increasing RunIndex and ScanNumber. NOTE: RunIndex is the original order of the RunData header lines (ignoring duplicates, in order of merged BNX files) [Default OFF]\n"
  "-sort-runindexShuffle <MinLen> [Seed] : First sort -i input maps in order of increasing RunIndex and ScanNumber, then shuffle the order so the 1st map for each RunIndex and ScanNumber is followed by the 2nd map for each RunIndex and \n"
  "                        ScanNumber. If followed by -subset, this equalizes the number of samples for each RunIndex and ScanNumber. Molecules shorter than <MinLen> are placed last within each\n"
  "                        RunIndex/ScanNumber bucket : This causes -subset to pick molecules shorter than <MinLen> only if there are not enough longer molecules.\n"
  "                        If <Seed> is specified, all maps are ordered randomly before applying -sort-runindexShuffle\n"
  "-randomize [Seed] : Randomize order of maps [Default OFF]. Seed is the random number generator seed (integer) and defaults to 1\n"
  "-subset <start> <end> [<frac> <perCohort>]: select subset of -i input maps specified by an integer interval corresponding to maps being numbered sequentially from 1 ... (after -randomize or any sorting)\n"
  "                        If <end> is a frnactional value less than 1.0, that fraction of the total number of maps will be selected\n"      
  "                        If <frac> is specified and != 0, at least that fraction of the total number of maps will be selected\n"
  "                        If <perCohort> is specified and != 0, the number of maps selected will be at least <perCohort> times the number of Cohorts (RunIndex & ScanNumber values)\n"
  "-subsetGB <GB> <MinLen> [ <rnd> ] : select subset of -i input maps consisting of first N molecules, where N is selected so the total size of molecules longer than MinLen KiloBases exceeds <GB> GigaBases\n"
  "                          If <rnd> = 1 : The value of N is rounded up to include all molecules with the same RunIndex and ScanNumber as the last molecule (requires -sort-runindexinc) [Default rnd= 1]\n"
  "-finalsort-sitesdec : Re-sort -i input maps in order of decreasing number of sites (after applying -sort or -randomize and -subset). Helps with load balancing, but increases memory usage\n"
  "-hashgen [<Win> <HashScore> <SDmax> <SDrms> <RelErr> <OffsetKB> <NumErrI> <NumErrP> [<NumResP>]] : Generate hashtable <prefix>.hashbin for all input maps using the following parameters:\n"
  "        <Win> : window size in number of site intervals [Default 5]\n"
  "        <Hashscore> : minimum number of windows that must match (with similar Offset, see <OffsetKB>) [Default 3]\n"
  "        <SDmax> : maximum sizing error of each interval as multiple of -sf -sd -sr values [Default 2.2]\n"
  "        <SDrms> :  maximum RMS sizing error average over entire window (as multiple of -sf -sd -sr values) [Default 1.2]\n"
  "        <RelErr> : maximum relative sizing error of each interval (absolute value, used when larger than <SDmax> based value) [Default 0.05]\n"
  "        <OffsetKB> : maximum Offset error in kb <OffsetKB> of matching windows [Default 3.0 kb]\n"
  "        <NumErrI> : maximum number of extra sites per window for inserted maps (FPs in query maps with -i) [Default 1] max value = 2\n"
  "        <NumErrP> : maximum number of extra sites per window for probe maps (FNs for reference maps with -ref) [Default 1] max value = 2\n"
  "        <NumResP> : maximum number of sites that failed to resolve per window for probe maps (reference maps with -ref) [Default 0] max value = 4\n"
  "          NOTE: Program will terminate after generating hashtable and must be re-run with -hash <prefix>.hashbin <HashScore> instead of -hashgen ..., OR use -hash in the same command\n"
  "    -id <ID> : Output file names are modified by changing all output filename prefixes to <prefix>_id<refid> : Useful when the same command is run for different -ref maps and separate output files for each is required\n"
  "         NOTE : the <ID> value is not used, it just needs to be non-zero. Instead the the reference map id is used in the filename prefix.\n"
  "               Even if multiple -ref maps are used together, seperate output files will be used for each refid, provided each input -ref file had a single map & refid (see also -mapped-unsplit).\n"
  "    -hash    : the program will continue with alignment using the hashtable just generated, then delete the hashtable file (Must be specified after -id)\n"
  "    -insertThreads <M> : Number of threads to use during map insertion stage [Default 1]\n"
  "                         Using M > 1 may increase memory usage due to the creation of multiple hashtables. Use caution since -hashmaxmem or -maxmem are not applied during hashtable creation\n"
  "                         Reasonable values would be 4-6 on Xeon-Phi co-processors with 8Gb memory (16-24 on compute nodes with 64Gb or more memory)\n"
  "    -queryThreads <N> : Number of threads to use during hashtable query stage [Default : Same as -maxthreads]. Will be reduced if -hashmaxmem or -maxmem are exceeded\n"
  "    -subsetmaps : Output subset of query maps with HashTable matches in <prefix>_hash.bnx (requires -ref)\n"
  "                  A mapping from the final bnx mapid value to the original mapid value used in the hashtable file is output in <prefix>.hashmap\n"
  "    -defermerge : Defer Merge Sort of HashTable parts : just output a list of file names in <prefix>.hashbin that can be combined by Merge Sort later (Cannot be used with -hashtxt)\n"
  "    -hashtxt : Output a text version of the hashtable in <prefix>.hash (instead of the binary <prefix>.hashbin) : used for debugging\n"
  "    -Hash_Bits <HASH_BITS> : Override value of HASH_BITS in hash.h [Default 21]\n"
  "    -MHash_Bits <MHASH_BITS> : Override value of MHASH_BITS in hash.h [Default HASH_BITS - 7]\n"
  "    -HHash_Bits <HHASH_BITS> : Override value of HHASH_BITS in hash.h [Default HASH_BITS - 7]\n"
  "    -hashtxt : Output a text version of the hashtable in <prefix>.hash (instead of the binary <prefix>.hashbin) : used for debugging\n"
  "    -hashMaxCov : Suppress hashtable entries for any key with more than this many entries [Default OFF]\n"
  "    -hashmatrix : output matrix of matching intervals (including matches from same map) to <prefix>.matrix and exit (no hashtable will be generated). Works best for single reference map\n"
  "    -hashoffset <N> : Add the number of matching windows with offset within <N> times <OffsetKB>. Only used with -ref (Value fixed at 1 for pairwise alignment) [Default 2]\n"
  "    -hashmaxmem : Maximum memory in Gigabytes for the hashtable. Default same as -maxmem\n"
  "    -HSDrange <V> [ <M> ] : Set hashtable maximum un-resolved site interval to res + max(V * resSD, M * res) : Values of V below 2.0 will speed up hashtable but miss some alignments due to \n"
  "                            unresolved sites. M can be useful for SV alignment since resSD can be misleading for assemblies [Default : 2.0 0.0]\n"
  "    -hashMultiMatch <L> [ <N> [ <D> ] ]: Allow multiple hashtable matches per map pair, as long as relative offsets (for same orientation) differ by at least <L> times <OffsetKB> [Default : OFF]\n"
  "                                 If specified, <N> limits the number of hashtable matches per map pair (in each orientation) (0 means No Limit) [Default : OFF or 0]\n"
  "                                 If specified, <D> filters out hashtable matches that score worse than the best hashtable match for same map pair by <D> matching windows [Default : OFF or 0]\n"
  "    -hashbest <D> : For each query map only return those hashtable matches (across all ref maps) within <D> of best score [Default : OFF or -1]\n"
  "                    Use -hashbest 0 to return only the best scoring matches (fastest alignment speed). Not recommended with -BestRef 0\n"
  "    -HashMaxHits <N> [ <X> <G> ]: During pairwise alignment limit number of hashtable matches per ref map to <N> : <N> may be rounded up to include all matches with same hash score as Nth best one [Default : OFF or 0].\n"
  "                  If X & G are specified adjust N to correspond to X times the genome coverage of the map input given Genome size G (Mb) [Default: OFF or 0].\n"
  "    -hashGC <N> : Reclaim hashtable memory every N probe map labels, to reduce memory useage for large probe maps [Default : OFF or 0]\n"
  "    -hashT2 <0/1> : Use 2 threads per probe map (nested parallelism) [Default : OFF or 0]\n"
  "    -hashGrouped <S> <H> : Make sure at least <H> matched windows with same offset are within OffsetKB x (2 ^ S) on the inserted map, otherwise ignore match for that offset [Default : OFF or 31]\n"
  "    -hashkeys <0/1> : Use 1-bit table for presence of 30-bit keys in hashtable. Does not change results and is faster except for high coverage data (alignmol, pairwise) [Default : OFF or 0]\n"
  "    -hashSF <V> <R> : Force hashgen to change -sf value to max(V, R * sf) [Default : 0.10 0.5]\n" 
  "    -hashSD <V> <R> : Force hashgen to change -sd value to max(V, R * sd) [Default : 0.06 0.5]\n" 
  "    -hashSR <V> <R> : Force hashgen to change -sr value to max(V, R * sr) [Default : 0.025 0.5]\n" 
  "-hash [[input.hash | input.hashbin] <threshold>] : use hashtable file and corresponding hash score threshold for pairwise and -ref alignment. Implies -nosplit 2\n"
  "                                         NOTE: If not used with -hashgen, must specify a filename (eg input.hashbin)\n"
  "                                         NOTE: -hash with no args defaults to -hash <prefix>_id<id>.hashbin <HashScore>, where <HashScore> is the -hashgen score parameter\n"
  "    -hashdelta <N> [<inc> [<lim>]]: limit alignment range so difference in aligned sites (from left ends) are within <N> x 2 x <OffsetKB> (see -hashgen) [Default: OFF] (NOT implemented for pairwise)\n"
  "                            If <inc> is specified the value <N> may be incremented by <inc> if initial alignment uses too much of the original alignment range (leaving a margin of less then <inc> )\n"
  "                            <N>/2 and <inc> should be both set to at least deltaX + 1 to avoid missing outliers located inside the endoutlier region\n"
  "                            If specified <lim> limits the largest value to which <N> can be incremented (to limit runtime)\n"
  "                            If -hashMultiMatch <L> .. is used, <inc> should be set to at least <L> / 2\n"
  "    -hashrange <0/1> : If 0, base -hashdelta range on distance (multiple of average interval distance) \n"
  "                       If 1, base -hashdelta range on number of labels [Default : 0]\n"
  "    -hashScaleDelta <0/1/2/3> : If 0, -ScaleDelta is ignored by the hashtable, except that sd and sr may be adjusted upwards\n"
  "                              If >= 1, -ScaleDelta is used by hashtable to try all possible scalings [Default 1]\n"
  "                              If >= 2, -ScaleDelta only checks matchgroups with scaling that scored best in hashtable (without -MultiMatches, there will be only one matchgroup per map pair)\n"
  "                                        Note : When using multiple -M iterations, only the first sub-iterations will be restricted to the scaling that scored best in hashtable\n"
  "                              If >= 3, -ScaleDelta with -MultiMatches only retains matchgroups with the same scaling as the best matchgroup (can be used without -hash)\n"
  "-merge : merge one or more cmap or bnx input files into a single cmap or bnx output file, while assigning unique map IDs if the same ID is present on 2 input maps\n"
  "              Any option that modifies or filters maps or labels on input can be used, along with the following options that are only used with -merge : no alignments are done\n"
  "              If -ref input is used, no -i inputs can be used and the -ref input maps are converted to -i input maps\n"
  "    -id <ID> : If > 0 and only one map is present, change the id of the map to <ID> [Default OFF]\n"
  "    -range <left> <right> : If only one contig map is present, output only the range (in kb) from <left> to <right> [Default OFF]\n"
  "    -inversion <left> <right> : If only one contig map is present, invert the interval (in kb) from <left> to <right> [Default OFF]\n"
  "    -swap : Swap the 2 colors in the maps [Default OFF]\n"
  "    -subsetbin <bin> <numbins> : select a subset that divides the file into <numbins> consecutive sections of roughly equal size and output subset <bin> (after -randomize or any sorting)\n"
  "                             A <bin> value of 0, means output all <bin> values 1 .. <numbins> and the output <prefix> is changed for each <bin> value to <prefix>_<bin>_of_<numbins>\n"
  "    -split : split all input contigs into separate output files <prefix>_contig<N>.cmap (or .bnx with -bnx) each with 1 map/contig with id <N>\n"
  "    -break <INC> <ID1> <breakpt1a> ... <breakpt1z> <ID2> <breakpt2a> .. <breakpt2z> ... : break apart specified contigs (integer ids) at specified locations (float values, in kb)\n"
  "                                       Output new fragments as contigs with original ID incremented by <INC> for each new fragment derived from the same original contig\n"
  "                                       NOTE: If breaks occur close to a label the fragment ends will be extended to a minimum value of 0.020 kb (see -minEnd to change this value)\n"
  "    -winbreak <Len> <Shift> <N> : break apart each contig into fragments of length <Len> kb starting at the left end and shifted by <Shift> kb until the right end of each contig is reached\n"
  "                     The Kth fragment of input contig I will have id = I * <N> + K. N will be increased to at least the number of fragments of the largest input contig.\n"
  "                      Useful to break a reference into smaller overlapped pieces. To break into pieces of 2200kb with 200kb overlap use -winbreak 2200 2000 1000 -minsites 10\n"
  "    -setMask <M> : Set End Mask values of all maps to M. A value of 0, clears all End Mask values. A value of -1 means no change to End Mask values [Default: OFF or -1]\n"
  "    -orMask <M> : Modify End Mask values of all maps by OR'ing the value M into any previous value. A value of 0 means no change to End Mask values [Default: OFF or 0]\n"
  "    -forcemerge <MergeFile> : Perform pairwise merge of contig maps based on specification in MergeFile. Each line of MergeFile specifies one pair of input maps to merge as 6 integers:\n"
  "        ContigID1 Orientation1 LabelID1 ContigID2 Orientation2 LabelID2 (Orientation = 1 means reversed, 0 means normal)\n"
  "        This specifies a merged map to consist of map w/ ContigID1 in Orientation1 from left end to LabelID1 followed by map w/ ContigID2 in Orientatio2 from LabelID2 to right end\n"
  "        The output will consist of all input maps, but with each pair of merged maps replaced by a single map, whose ID will be the smaller of the two component map IDs\n"
  "   -fasta : output fasta file for each input contig, using base 'A' for all regions between labels. Contig <N> will produce <prefix>_contig<N>.fasta [Default: OFF]\n"
  "   -incID <N> : increment all -i input map IDs by <N> (a 63-bit value) [Default: OFF or 0]\n"
  "-FP <False Positives per 100kb, one per color>[Default 1.0]\n"
  "-FN <False Negative Probability, one per color>[Default 0.15]\n"
  "-sd <sizing error in root-kb, one per color>[Default 0.15]\n"
  "-sf <fixed sizing error in kb, one per color>[Default 0.10]\n"
  "-sr <relative sizing error, one per color>[Default 0.00]\n"
  "-se <resolution limited sizing error, one per color>[Default 0.00]\n"
  "-ma <mean> <sf> <sd> : 2-color misalignment values between sites of opposite colors : bias = <mean>, var = <sf>*<sf> + <sd>*<sd> * TrueDistance [Default : 0.0 1.0 0.10]\n"
  "-bpp <Scaling/500> : Example : -bpp 510 will rescale all -i input BNX or CMAP map sizes by 510/500 [Default : 500 for BNX or CMAP files]\n"
  "                    NOTE : This option no longer has anything to do with pixels and is named this way for historical reasons. It specifies a re-scaling factor relative to 500\n"
  "                    NOTE : If -pairsplit or -sv is specified -bpp is NOT applied to the first -i input map set (treated as reference map set)\n"
  "-NoBpp : Do NOT apply -bpp correction to -i input maps (includes bpp value from -readparameters file). Useful if -i input has already been scaled (eg by using -merge -bpp OR with _q.cmap)\n"
  "                    NOTE : This overrides the bpp values from either -bpp OR -readparameters with a value of 500, so it is not the same as NOT specifying bpp at all for .vfixx input files\n"
  "-ScaleDelta <stepsize> <Steps> : Try different per molecule scaling values of 1 + <stepsize> * (1, 2 ... <Steps>). Example: -ScaleDelta 0.05 1 will rescale molecules by 1 +- 0.05 [Default: OFF]\n"
  "                                 NOTE : Molecules are rescaled before estimating error parameters (see -M), so the same -ScaleDelta must be used for Assembly/Refinement as during AutoNoise\n"
  "                                 NOTE : not validated for pairwise alignment, except with -pairmerge\n"
  "       -ScaleDeltaBPP [ <N> ] : Estimate error parameters (see -M) based on original unscaled molecules : Recommended if -ScaleDelta will NOT be used for Assembly/Refinement but only during AutoNoise\n"
  "                                Also required if using -MapScale, so map scaling factor is based on original unscaled molecules\n"
  "                    If <N> is specified : Only the first N major iterations of -M will apply -ScaleDeltaBPP : requires use of -ScaleDelta in Assembly/Refinement, but allows -ScanScaling -resbias\n"
  "                                which will only be used in the first N major iterations of -M after which the rescaled BNX will be output before continuing without -ScaleDeltaBPP (NOT YET IMPLEMENTED)\n"
  "                    In the absense of -ScaleDeltaBPP average molecule scaling will be adjusted to equal 1 (NOT YET IMPLEMENTED)\n"
  "-MapScale [ <stepsize> <Steps> ] : Apply per map scaling correction before each -M iteration based on alignments from the previous iteration, with average scaling correction equal to 1 [Default OFF or 0]\n"
  "              Internal Outliers in the Alignment are excluded when computing the scaling factors\n"
  "              For the first iteration alignments are repeated after applying the per map scaling : the first alignments use less stringent -T -A -S -L thresholds and -BestRef 0\n"
  "              If <stepsize> <Steps> are specified, -ScaleDelta <stepsize> <Steps> -ScaleDeltaBPP is applied for the initial alignments only. Subsequent -M iterations will only use \n"
  "              -ScaleDelta if specified seperately. Molecules output by -mapped will be rescaled by -MapScale, while molecules output in <prefix>_rescaled.bnx by -ScanScaling will not reflect -MapScale\n"
  "              If -MapScale is used during error parameters estimation (see -M) sizing errors estimated will be smaller, hence -MapScale must be used subsequently during Assembly/Refinement\n"
  "-res <r> : r is Resolution in pixels, one value per color [Default 3.5] : NOTE : Depends on bpp value\n"
  "-resSD <rs> : rs is standard deviation of resolution in pixels, one value per color [Default 0.75] : NOTE : Depends on bpp value\n"
  "-mres <pixels> : Reduces resolution of -i input maps to specified multiple of 0.5 kb. [Default: 0.001]\n"
  "-mresSD <pixelsSD> : Varies the resolution of -i input maps by randomly varying -mres value by specified Gaussian SD. [Default 0]: introduces randomness in map resolution\n"
  "                    NOTE : -mres and -mresSD are applied after -bpp (including -bpp from -readparameters)\n"
  "-rres <pixels> : Reduces resolution of -ref input maps to specified number of multiples of 0.5 kb. Also limit resolution of refined consensus maps to this value [Default: 0]\n"
  "-local <pvalue> <fdr>[Default 0.0, 0.1] : Output chimeric molecules in <prefix>.chim (if -endoutlier is used, the larger <pvalue> is used)\n"
  "-nosplit <0 allows unlimited chimeric map splits (with -A sites) at outlier ends, 1 only allows one split, 2 allows no splits> [Default OFF or 2]\n"
  "-pairsplit <Psplit> : compute multiple local alignments during pairwise alignments by breaking alignments (and molecules) at internal outliers\n"
  "                      with specified Pvalue <Psplit> AND at outlier ends. Requires -endoutlier [Default 0.0]\n"
  "                      The first -i input file is assumed to contain reference maps and the remaining -i input files are assumed to contain query maps (see -first to modify this assumption)\n"
  "     -pairsplitExtend : When breaking molecules at internal outliers extend both pieces to include part of the outlier region (including any other outliers nearby) up to A - 1 sites (see -A)\n"
  "     -pairsplitMerge : Try to merge overlapped or adjacent alignments [Default: OFF]\n"
  "     -pairsplitMinOutlier <N> : Do NOT break alignments at internal outliers with less than <N> site intervals on both maps (intended to work with -indel) [Default OFF]\n"
  "     -pairsplitPV [<0/1>] : Prioritize alignments by LogPV (confidence) rather than LogLikelihood score [Default: 0 or OFF]\n"
  "     -PairwiseOutlierFix [<0/1>] : enable fixes in handling end outliers during pairwise alignment (recommended with -pairsplit or -sv) [Default : ON]\n"
  "     -Palindromic PvalueRatio : Mark alignment as Palindromic if the reverse alignment has Pvalue within this ratio of normal (-T) threshold [Default : 0.01]\n"
  "     -pairsplitRef [<0/1/2>] : If != 0 Assume first -i input is a in-silico Reference Map set [Default : 2]\n"
  "                              -bpp is not applied to these maps and (if >=2 ) no sizing errors are assumed for them.\n"
  "     -SeparateQueries <0/1> : With multiple query maps, process each query map separately from other query maps [Default : 1]\n"
  "-querysplit <Psplit> : Same as -pairsplit but assume 1st set of maps do not need to be split (ie have no chimeric junctions or duplicated regions)\n"
  "-first <N> : Perform pairwise alignments only between any one of the first N maps and one of the other maps.\n"
  "             -first -K takes N value from number of maps in first K input map files [Default OFF (-1 if -pairsplit)]\n"
  "    -firstPairs 1 :  Also include pairs when both maps are within first N maps\n"
  "    -firstPairs 2 :  Also include pairs when both maps are outside of first N maps\n"
  "-pairmerge [<minoverlap> <maxendoutlierFrac> <maxoutlierKB> <maxendoutlierKB> [<maxendoutlierExpandKB> <maxendoutlierN> <maxoutlierM>] : output pairwise merged maps in CMAP format when alignment satisfies -S -T -A -L thresholds and has a minimum overlap specified in kb [Default OFF or 0]\n"
  "             Note that the overlap normally includes internal outliers and averages the size for the two map intervals that are aligned : See -pairmergeExcludeOutliers to exclude outliers\n"
  "             If maxendoutlierFrac is specified, any misaligned regions caused by -endoutlier may not exceed the specified fraction of the overlap region [Default 0.20]\n"
  "             If maxoutlierKB is specified, any internal outliers caused by -outlier may not exceed the specified size difference in kb [Default OFF]\n"
  "                   Any internal outliers near the ends, within maxendoutlierExpandKB (defaults to maxendoutlierKB, if not specified) is treated as part of the endoutlier\n"
  "             If maxendoutlierKB is specified, any misaligned ends caused by -endoutlier may not exceed the specified size in kb [Default OFF]\n"
  "             If maxendoutlierN is specified, misaligned ends over maxendoutlierKB (or maxendoutlierFrac) are OK if there are less than maxendoutlierN labels in the endoutlier of shorter map.\n"
  "                              (NOTE: Does NOT count overlapped misaligned labels of larger map!) [Default OFF or 0]\n"
  "             If maxoutlierM is specified, any internal outliers may not have more than maxoutlierM misaligned labels even if the outlier size difference is below maxoutlierKB [Default OFF or 999]\n"
  "             Each input map must have distinct ID values (across all input files)\n"
  "             Unmerged maps are copied unchanged to PREFIX_contig<N>.cmap, where <N> is the ID of the map.\n"
  "             Merged map pairs (starting with highest Pvalue) are written to PREFIX_contig<M>.cmap, where <M> is the smaller \n"
  "             of the two map IDs (the higher map ID will be missing from the output files). The merged map transitions between the\n"
  "             original maps at the midpoint of the alignment. Conflicts between pairwise alignments are resolved in favor of the best Pvalue (Cannot use -ref) \n"
  "      -outlierBalanced <N> : For balanced outliers with more than <N> misaligned labels, replace the outlier size with the smallest inversion size (if larger) [Default: 1]\n"
  "      -pairmergeOutlierChiSq <PV> <F> [ <L> <Cmax> <S> <R> <I> <SR> <Cmin> <CR> <Ch> <CRt>]: Ignore outliers with size below L kb if smaller map interval (or any sub-interval) has CMAP OutlierFrac value above F % OR \n"
  "                                 MolCov below Cmax AND ChiSq score below PV AND MolSd above ExpSd + S kb + R times interval size AND MolSd above ExpSd * SR AND larger map interval size is at least I kb\n"
  "                                 Also ignore outliers with size below L kb AND larger map interval size of at least I kb IF MolCov of smaller map interval is below both Cmax and \n"
  "                                            Cmin + CR * RangeRatio (ratio of larger map interval size to average interval.\n"
  "                                            If CRt > 0 : If qry map and ref map have different number of unaligned labels, limit RangeRatio to no more than CRt\n"
  "                                 Requires that CMAP was refined with -OutlierQuality 2 to generate OutlierFrac and MolSd, ExpSd, MolCov and ChiSq values for each interval. \n"
  "                                 If Ch > 0 : Replace MolCov with min(MolCov, ChimNorm * ChimQuality * Ch / 100)\n"
  "                                 Can be used to ignore outliers caused by DNA knots in the molecules. [Default OFF or 0 100]\n"
  "      -preserverOutlier <L> <LE> : preserve outliers over L kb (and under <maxendoutlierKB>) in overlap region by cloning the part of the overlap discarded and bias the merge transition point to lie \n"
  "             on one side of all such outliers. Also preserve endoutliers over LE (and under <maxendoutlierKB> kb) by copying the endoutlier region along with the overlap region. \n"
  "             Useful with -pairmergeHmap to preserve outliers without requiring maxoutlierKB (-see pairmerge) to be set so small that extension merges are prevented [NOT YET IMPLEMENTED]\n"
  "             Useful with -SplitSegDup to encourage merging ends of contigs with endoutliers between L kb and SplitSegDup <E> value, without losing the endoutlier (alternate Allele).\n"
  "      -pairmergeOverlapped <L> <Margin> [ <S> <N> <E>]: Reduce minimum overlap required to <L> when one map completely overlaps another map (except by margin not exceeding <Margin>) [Default OFF or 0]\n"
  "                           Prefer to keep map that is larger by S kb OR has N extra labels (and is within E kb of other map) [Default OFF or 0 0.001 1.0 1 0.1]\n"
  "      -pairmergeExcludeOutliers <0/1> : Exclude internal outliers from the computation of the overlap, which must exceed <minoverlap> and is used by other suboptions [Default OFF or 0]\n"
  "      -preferNGS : when merging preserve as much of the NGS map (StdDev=0) as possible : transition based on minimizing average StdDev of result\n"
  "      -pairmergeRepeat : repeat merge operation until no more maps merge (faster than calling the basic command repeatedly)\n"
  "             Output file names (eg XMAP) will have M1, M2 etc attached to the -o prefix\n"
  "      -TrimmedNoMerge <A> <End> [<D> <S> <N>]: Don't merge any pair of maps if any overlapped alignement ends are marked as trimmed or non-extendable UNLESS either the aligned length is at least <A> kb\n"
  "             OR one map is fully overlapped by the other larger map (except for an extension of no more than <End> kb). \n"
  "             <A> should be set to around 1.5x the value used by -CloneLimitLen and -Truncate.  [Default OFF or 0]\n"
  "             Note that -TrimmedNoMerge is ignored for final merges with -pairmergeHmap except for ends of contigs broken (and marked) by -SplitSegDup.\n"
  "             Note that <A> value is ignored (assumed very large) if any overlapped end is marked as a chromosome end (or fragile site contig end).\n"
  "             If <D> is specified non-zero, it will allow merges with overlapped alignments if only one of the ends is marked, or both ends are marked with bitwise non-intersecting values, ignoring\n"
  "                     end mask bits not present in <D>. In addition if <S> is non-zero, any internal label with mask intersecting with <S> will be considered part of a SegDup region and merges will\n"
  "                     be disallowed unless the overlap includes at least <N> labels beyond the SegDup region on the side NOT being extend by the overlap\n"
  "                     By default (<D> = 0) a mark on either end disables merges. [Default : OFF or 0 0 0 0 0 0] NOTE : <D> and <S> are bitmasks & can be represented in Hex (eg 0x18 instead of 24)\n"
  "      -TrimmedNoMergeStrict <0/1> : Also dont't merge any pair of maps if both maps have any end marked as trimmed or non-extendable, even if the marked end is not overlapped by the alignment.\n"
  "             This is stricter than the normal -TrimmedNoMerge and overrides any merges by -TrimmedNoMerge except when one map is fully overlapped by the other larger map [Default: OFF or 0]\n"
  "      -Truncate <L> <N> <E> [<F> <M>] : With partial alignment between two maps of at least L kb and N aligned labels and no outliers larger than E and only one side having an endoutlier\n"
  "             truncate the shorter map from the end without the endoutlier to no more than L kb or N aligned labels and mark the trimmed end as non-extendable. It is assumed that the two maps\n"
  "             are alternate Alleles and that the extension of the shorter one can be inferred from the extension of the longer one : this assumption requires use of a small E of 15kb or\n"
  "             less and L at least 200kb. Works well with -extsplit to limit the size of the bubble Allele generated to L kb on both sides of the bubble and should use the same L value as\n"
  "             the -CloneLimitLen sub-option of -extsplit. Outliers are only required to be below E kb (& <= M misaligned labels) over the trimmed region + L kb (or N labels) [Default : OFF or 0 0 0 1 3]\n"
  "             If <F> is 1 then Don't truncate map at one end if it is the longer map at the other end\n"
  "      -BigMapsFirst <0/1> : Check pairwise alignments in descending order of Map size rather than just descending order of alignment confidence [Default : OFF or 0]\n"
  "                    Helpful after -extsplit with -Truncate so all split contigs are first truncated by aligning with the larger contigs they were split from.\n"
  "      -MaxPalindrome <L> <N> <SL> [<D>] : Check for end overlap alignments of any Contig Map with its inverted self with at least L kb and N labels and split the Contig at the midpoint of the alignment + SL kb\n"
  "                                    The aligned region must have no outliers or endoutliers larger than required by -pairmerge.\n"
  "                                    If D > 0 : The split ends are marked as non-extendable using flag value D : D defaults to 0x1\n"
  "                                    Useful for breaking Chimeric extensions caused by Palindromic regions over 60kb (eg -MaxPalindrome 400 40 60) [Default: OFF or 0]\n"
  "      -MaxInvPalindrome <L> <N> <SL> [<D>] : Check for end overlap alignments with 1 endoutlier of any Contig Map with its inverted self where the aligned region is at least L kb and N labels and one\n"
  "                                  endoutliers is at least as large as the aligned region + SL kb and the other endoutlier and internal outliers are no larger than permitted by -pairmerge.\n"
  "                                  Alternately there is no large endoutlier but one large internal outlier of at least SL kb that spans the center of the aligned region and each aligned side \n"
  "                                  is at least L kb and N labels long. The contig will be broken into two overlapping contigs, where the overlap region is the large internal outlier or the \n"
  "                                  excess/unique region in the larger end outlier. If D > 0, the broken ends are marked as non-extendable with flag value D : D defaults to 0x1\n"
  "                                  Useful for breaking Chimeric extensions caused by Heterozygous inversions over 60kb (eg -MaxInvPalindrome 200 20 60) . [Default: OFF or 0]\n"
  "      -SplitSegDup <L> <N> <E> <EN> <C> <LM> <NM> [ <MoutKB> ]: With partial alignment between two maps of at least L kb and N labels, with one or both sides having both endoutliers over E kb (& EN labels), and \n"
  "          both sides having at least one map with unaligned length of E kb (& EN labels), assume the aligned region is a SegDup, provided the total coverage of both maps in the overlap regions averages \n"
  "          C times the Genome average : split the maps (with segdup region included in both splits), but only keeping ends with <maxendoutlierKB> ( & <maxendoutlierN> lables) beyond the overlap region, \n"
  "          and mark the cut ends non-extendable. Instead if neither side has an endoutlier of at least E kb & EN labels AND overlap region exceeds LM kb and NM labels AND coverage does NOT exceed C times \n"
  "          Genome average, then merge the two maps but also save copies of each endoutlier over <maxendoutlierKB> ( & maxendoutlier<N> labels) with the overlap region as a new Contig. This is applied only\n"
  "          after no more regular merges are possible (but before -pairmergeHmap) and any subsquent merging of the split maps ends is disabled (as with -TrimmedNoMerge 1000000)\n"
  "          Requires either -pairmergeRepeat OR that pairmerge was previously run until no more map pairs can merge. [Default: OFF] [<C> coverage requirement OR merger NOT YET IMPLEMENTED]\n"
  "          Internal outliers above <maxoutlierKB> cannot be included in the partial alignment and alignments are adjusted (if needed) as if <maxendoutlierExpandKB> were unlimited (large), but only\n"
  "                   if original alignment already has one endoutlier over E kb (& EN labels).\n" 
  "          If -pairmergeExcludeOutliers is specified, the partial alignment must be above L kb and N labels without counting any outliers.\n"
  "          To avoid false splitting due to -extsplit side contigs E (and EN) should exceed the values for L (and N) used with -CloneLimitLen & -Truncate\n"
  "          -SplitSegDup is applied in a seperate stage (with multiple rounds of pairwise alignments) after processing all other -pairmerge merges but before processing -pairmergeHmap\n"
  "          If <MoutKB> is specified, its value replaces <maxoutlierKB> during the -SplitSegDup stage\n"
  "      -SegDupMask <M> : Mark overlap region of both maps for segdup alignments (see -SplitSegdup, but applied even if overlap is smaller than <L> or <N>) with this Mask value [Default : 0x4]\n"
  "                          Also used to specify that this Mask value should be propagated through subsequent refinement\n"
  "      -SplitSegDup_OneSided <K> <L> : If K == 0 : Ignore SegDup alignments (for -SplitSegDup AND -SegDupMask) with only one side having an endoutlier over E kb (& EN labels)\n" 
  "                                      If K == 1 : Allow SegDup alignments with only one side having an endoutlier over E kb (& EN labels) [Default : 1 0]\n"
  "                                      If K == 2 : Like K==1, but ignore SegDup alignments if end without endoutlier has been truncated (see -Truncate), unless endoutlier size exceeds L kb\n"
  "      -BreakEndOutlier <L> <N> <E> <EN> : With partial alignment between two maps of at least L kb (& N Labels) that cannot be merged due to endoutliers or outliers, each map is broken at the \n"
  "           junction of each endoutlier over E kb (& EN labels). Useful for mrg0 to break up chimeric contigs and (for roughdraft mrg0) to break up regions with different scalings [Default OFF]\n"
  "      -BreakOutlier <L> <N> <I> : With alignment betweeen two maps of at least L kb (& N Lables) and no endoutliers that cannot be merged due to large internal outliers over I kb, \n"
  "           break the smaller map at the largest internal outlier. Useful for roughdraft mrg0 to break up regions with different scalings [Default OFF]\n"
  "      -pairmergeHmap <L> <N> [ <LE> <LenE> ]: After no more map pairs can merge and after -SplitSegDup has been applied, re-examine all pairwise alignments and try to combine as many as possible\n"
  "             into Haplotype Maps and output them as .hmap contigs. First merge Bubbles (overlapped maps that could not be merged due to an internal outlier larger than <maxoutlierKB>)\n"
  "             The Bubble must overlap by at least L kb and N labels on both ends on either side of outliers larger than Truncate's <E> kb, with no outliers larger than <E> kb\n"
  "             Any merge that produces conflicting Het indels will be suppressed : currently any overlapping Het indels (even if identical) will prevent the merge\n"
  "             For end overlap merges with Bubbles additional restrictions apply due to the higher risk of chimeric merges : overlap must be at least <LE> kb and both maps must be at least <LenE> kb long.\n"
  "                   If not specified <LE> and <LenE> default to <minoverlap> and 0\n"
  "             If -SplitSegDup is specified, allow larger endoutliers up to Truncate's <E> for end overlaps, but also save copies of each endoutlier over <maxendoutlierKB> along with the overlap\n"
  "                   region as a new contig[NOT YET IMPLEMENTED].\n"
  "             Requires -Truncate (for 3rd arg <E>). Also requires either -pairmergeRepeat OR that pairmerge was previously run until no more map pairs can merge\n"
  "      -pairmergeXMAP : output XMAP files to display pairwise alignments (Also works without -pairmerge for any pairwise alignments)\n"
  "      -CIDrenumber <0/1> : renumber contig ids to smallest possible positive IDs (sequentially starting at 1) [Default: OFF or 0]\n"
  "-delta <N> : N = largest pairwise alignment interval in sites (maximum misaligned sites per interval is N-1) [Default 4]\n"
  "-deltaY <N> [<M>] : N = largest reference map alignment interval in labels [Default 6]\n"
  "                    If specified M is a smaller value that is only used before the last iteration of -M [Default OFF or same as N]\n"
  "-deltaX <N> [<M>] : N = largest query map alignment interval in labels [Default 4]\n"
  "                    If specified M is a smaller value that is only before the last iteration of -M [Default OFF or same as N]\n"
  "-outlierExtend <N> [ <Lim> ] : increase deltaX,deltaY as needed to enlarge internal outliers up to N sites beyond the current outlier boundaries  [Default 0]. [Requires -biaswt 0 with -ref]\n"
  "                   If Lim is specified, outliers are not enlarged beyond Lim sites [Default OFF]\n"
  "                   If two values of -deltaX or -deltaY are specified, -outlierExtend will only be used in the last iteration of -M\n"
  "-extend <N> : If N > 0 allow query map to extend beyond reference map (if N = 2, also allow refinement of extended reference)  [Default 0]\n"
  "              NOTE: refalign does not score (or Bias) the part of the query map extending beyond the reference, unlike refinement (which treats all such sites as FP).\n"
  "-thetaScale <F> : Scale (x/theta) score term in refalign alignment by F. Lower values increase sensitivity to FN labels, increasing likelihood of endoutliers and outliers [Default: 1.0]\n"
  "-thetaScaleRef <F> : Scale (x/theta) score term during refinement by F. Lower values increase sensitivity to FN labels, increasing likelihood of endoutliers and outliers [Default: 1.0]\n"
  "-outlier <Pvalue> [<StitchPV>]: outlier Pvalues per aligned interval. If -stitch, then <StitchPV> specifies a seperate Pvalue for intervals with an image stitch [Default 0.0003 0.01]\n"
  "    -outlierMax <delta> : Disallow outliers with sizing difference larger than <delta> in kb. NOTE : only used with -maptype 0 and helps speed. OFF if <delta> is greater than 1000.0 kb. [Default OFF]\n"
  "    -outlierLambda <lambda> : Decrease -outlier Pvalue by exp(-|x-y|/lambda), where x-y is sizing difference of the aligned interval in kb. OFF if <lambda> is greater than 1000.0 kb. [Default OFF]\n"
  "                              If maptype = 1 (and -outlierType1 0) : Instead decrease -outlier Pvalue by exp(-(x+y)/lambda)\n"
  "    -outlierType0 <W> : Modifies outlier penalty when maptype==0 [Default 0]\n"
  "                          If W=0, don't penalize misaligned sites inside an outlier (causes more intervals to be classified as outliers based on misaligned sites and sizing error)\n"
  "                          If W=1, Penalize misaligned sites, even inside an outlier (causes fewer intervals to be classified as outliers based only on sizing error)\n"
  "                               Values of W between 0 and 1 may be used to partially penalize misaligned sites inside and outlier\n"
  "    -outlierType1 <W> : Modifies outlier penalty when maptype==1 [Default 0]\n"
  "                          If W=0, don't penalize misaligned sites inside an outlier\n"
  "                          If W=1, Penalize misaligned sites, even inside an outlier\n"
  "                               Values of W between 0 and 1 may be used to partially penalize misaligned sites inside and outlier\n"
  "    -outlierNorm0 <0/1/2> : Modifies outlier penalty normalization factor when maptype==0 [Default 0]\n"
  "                          If >= 1, don't adjust outlier penalty normalization factor based on -outlierLambda value \n"
  "                          If >= 2, don't adjust outlier penalty normalization factor based on -outlierMax value \n"
  "    -outlierNorm1 <0/1> : Modified outlier penalty normalization factor when maptype==1 [Default 0]\n"
  "                          If 1, don't adjust outlier penalty normalization factor based on -outlierLambda value \n"
  "    -outlierBC : Apply BonFerroni correction to outlier score : makes outlier Likelihood more stringent by number of site intervals in outlier region of both maps [Default OFF]\n"
  "-endoutlier <pvalue> : end outlier(blemish or chimeric junction) prior probability per true positive site and -K bias includes unaligned end region of nanomap. Also used during refinement [Default 0]\n"
  "-maptype <n> : Specify if -i input (Query) map is Molecule (<n> = 0) or Assembly/Contig (<n> = 1). By default .vfixx or .bnx is assumed to be Molecules and .cmap Assembly/Contig\n"
  "                   Molecules are only permitted to have outliers with sizing error smaller than -outlierMax (representing knots and stitch errors)\n"
  "                   Assembly/Contig may have outliers of any size (representing SV insertions or deletions relative to the reference)\n"
  "-biaswt <WT> [ <PR> ] : Scale bias added to alignment score (log likelihood ratio) by WT. For WT=1 bias value is determined based on score expectation at a percentile specified by -K,-KL,-KF Bias [Default 1.0 0.0]\n"
  "               Use WT=0 to disable -K,-KL,-KF based bias (in that case alignment score is just the log likelihood ratio\n"
  "               Use <PR> == 1.0 to correct score expectation to include misresolved label penalties (OFF by default) NOT YET IMPLEMENTED\n"
  "               During refinement -biaswt is not used : Use -biaswtRefine to apply -biaswt during refinment\n"
  "-biaswtEnd <WT> : biaswt value to use for endoutlier regions of query maps, as a multiple of biaswt value. NOTE: Value is always 0 during refinement and pairwise alignment. [Default 1.0]\n"
  "               NOTE : biaswtEnd = 1 will not change which regions are classified as End outliers as -biaswt is changed, while biaswtEnd = 0 will classify more low scoring regions as End outliers\n"
  "-biaswtOutlier <WT> : biaswt value to use for outlier regions of query maps, as a multiple of biaswt value. [Default 1.0]\n"
  "               NOTE : biaswtOutlier = 1 will not change which regions are classified as outliers as -biaswt is changed, while biaswtOutlier = 0 will classify more low scoring regions as outliers\n"
  "-biaswtRefine <WT> <type> [ <PR> ] : Apply -biaswt during refinement (If -extend 2 then sites are added/deleted in the extension region before using -biaswt OR -LRbias, since initial score in the extension region will be very low)\n"
  "                            <type> = 0 means only use -biaswt <WT> when computing final coverage and Chimeric Quality score\n"
  "                            <type> = 1 means use as soon as -LRbias is applied, after initial refinement stages (including Haplotyping) \n"
  "                            <type> = 2 means use at all times, even during initial refinement stages [Default 0.0 0 0.0]\n"
  "                            Use <PR> != 0 to correct score expectation to include misresolved label penalties (OFF by default) NOT YET IMPLEMENTED\n"
  "                            NOTE : -biaswtRefine is NOT applied to endoutliers : use with -biaswtEnd 0.0 for consistency between refalign and refinement\n"
  "                            NOTE : -biaswtRefine IS applied to internal outliers, but modified by -biaswtOutlier\n"
  "-alignscore : output detailed pairwise alignment scores in .align[Default OFF]\n"
  "-SecondBest <MinSpacing> [ <Pvalue2> <score2> ] : output the 2nd best alignment score (with right end at least MinSpacing * QueryLength from best alignment) to .xmap\n"
  "                         Only 2nd best alignments above the specified <score2> and <Pvalue2> threshold are output [Default: Pvalue2 = 1, score2 = -1000]\n"
  "               NOTE: If used with -hashgen, requires -hashMultiMatch\n"
  "       -SinglelineXmap : The 2nd best match is output on the same line as the best match (non-standard xmap format) : disables use of <Pvalue2> & <score2> thresholds\n"

  "-MultiMatches <UniqueLabels> [ <T2> <S2> <A2> <L2> <E2> <I2> <M2>] : Output all alignments (in addition to best alignment), provided each alignment has at least the specified number of unique labels\n"
  "                         in either Query or Reference relative to every higher scoring alignment. A value of 0 turns off -MultiMatches [Default : OFF or 0]\n"
  "                         If specified <T2> and <S2> are a lower Pvalue and score threshold for alignments other than the best alignment [Default : same as previous -T and -S values]\n"
  "                         If specified <A2> is a lower minimum aligned sites threshold [Default : same as previous -A value]\n"
  "                         If specified <E2> is a lower maximum Endoutliers threshold [Default : same as previous -E value]\n"
  "                         If specified <I2> is a lower maximum internal outlier threshold [Default : same as previous -I <D> value]\n"
  "                         If <M2> is specified : M2 - 1 is the maximum number of internal outlier misaligned labels [Default : same as previous -I <M> value\n"
  "                         NOTE : If -T2 is specified and less stringent than -T, the alternate matchgroups are only retained if at least one matchgroup exceeds the -T threshold\n"
  "                         NOTE : The dependence on previous -A -E -I -T and -S values means subsequent -A -E -I -T -S values will not undo the effect of previous values on -MultiMatches\n"
  "    -MultiMatchesDelta <D> : Also output overlapped alignments that do not have enough unique labels, provided the offset between maps differs by at least <D> kb [Default : OFF or 0]\n"
  "    -MultiMatchesFilter <V> : V >= 1 : Filter out matchgroups below <Pvalue2> and <A2> thresholds early on, to avoid keeping higher scoring overlapped matchgroups that are below threshold\n"
  "                              V >= 2 : Earlier Filtering out of matchgroups that overlap with a higher scoring matchgroup for all except the rightmost <A2> - 1 aligned labels (faster) [Default : OFF or 0]\n"
  "    -MultiMatchesUniqueQuery <Q> <I> <T>: Filter out an alignment if another alignment with same query has higher confidence, and lower confidence alignment has fewer than Q unique aligments either end.\n"
  "                        If I >= 1 : Avoid filtering if I consecutive aligned qry labels of lower confidence alignment are NOT aligned in higher confidence alignment even if overlapped on both sides\n"
  "                                    Allows small inversions and translocations where the small alignment overlaps a large insertion outlier in the query of the larger alignment\n"
  "                        If T >= 1 : Trim lower confidence alignment where it overlaps the higher confidence alignment on the query : trim the alignment with lower confidence in the overlap region\n"
  "                        Applied before computing alignment weights with -BestRefWT 1. (I >= 1 : NOT YET IMPLEMENTED) [Default: OFF or 0 0 0]\n"
  "    -MinDupSites <N> : Minimum number of sites for non-inverted Duplication region (must be at least 1) [Default : 1]\n"
  "    -CutFlip <N> <Err> <Exp> <EndExp> [ <Pvalue3> <score3> <Pvalue4> <score4> <L> ] : Try to align query region of any internal outlier with <N> or more labels (on both query & ref) and relative error less than <Err>, \n"
  "                              in reverse orientation to the same query/ref region expanded by up to <Exp> labels on each side. The Flipped alignment found must satisfy at least the 2nd set of -MultiMatch thresholds \n"
  "                              (<Pvalue2> etc) OR <Pvalue3> & <score3> (if less stringent) AND the original (unflipped) alignment must satisfy the primary thresholds (-T etc) [Default : OFF or 0]\n"
  "                              If <EndExp> is non-zero, qry end outliers are also aligned in reverse orientation with either endoutlier in the ref, but ref is expanded in non-aligned direction by <EndExp> labels.\n"
  "                              However endoutlier inversions are suppressed if there is already another matchgroup extending into either the qry or ref endoutlier regions by <L> or more labels (Default same as <N>).\n"
  "                              For end outlier inversions, a more stringent <Pvalue4> & <score4> can be substituted for <Pvalue3> & <Score3> (see also -CutFlipEndBC) \n"
  "      -CutFlipBC <0/1/2>      If >= 1 : Apply local Bonferroni correction to Pvalue for end outlier CutFlip inversions matchgroups based on length of endoutlier vs minimum inversion size <N>\n"
  "                              If >= 2 : Also apply Bonferroni correction to internal outlier CutFlip inversion matchgroups based on length of internal outlier vs minimum inversion size <N>\n"
  "      -CutFlipSkip <M>        Save time by skipping -CutFlip for outliers with <M> or more labels on the Query map : These cases are expected to be handled by regular -RefSplit matchgroups [Default : OFF]\n"
  "      -CutFlipMerge <Psplit>  When two outliers are seperated by a aligned region with Alignment Score no better the -log(<Psplit>), merge the outliers before applying -CutFlip [Default : Same as -RefSplit <Psplit>]\n"
  "      -CutFlipFilter <MaxOverlap> [ <MaxMiss> <OverlapFilter> <FilterConf> ] : Filter out (ignore) CutFlip inversions with inverted aligned labels overlapping non-inversion aligned labels by more than\n"
  "                <MaxOverlap> fraction (ignoring palindromic inversion ends of smaller matchgroup) [Default 0.5]\n"
  "                Also Filter out CutFlip inversions if gap in either breakpoint has missaligned labels that exceed the fraction <MaxMiss> of the inverted aligned labels [Default 0.5 0.0 0.0 0.0]\n"
  "         If <OverlapFilter> != 0 : Before detecting Small SegDups and CutFlip Inversions, filter out matchgroups overlapped by larger matchgroups (same qry and ref) if query overlap fraction (in kb)\n"
  "           exceeds both <OverlapFilter> and <QryRatio> (from -OverlapFilter), but save the smaller matchgroup if one of the following three conditions apply :\n"
  "              1. Larger matchgroup has opposite orientation AND fully overlaps smaller matchgroup in reference AND confidence of larger matchgroup does NOT exceed that of smaller matchgroup by more \n"
  "                 than <FilterConf> AND aligned label overlap fraction in query does NOT exceed larger of 0.5 and <MaxOverlap> (ignoring palindromic inversion ends of smaller matchgroup).\n"
  "              2. Larger matchgroup has opposite orientation AND does NOT fully overlap smaller matchgroup in reference AND aligned label overlap fraction in query does NOT exceed <MaxOverlap> AND\n"
  "                 query kb overlap fraction does NOT exceeed <InvQryRatio> of -InversionOverlapFilter AND confidence of larger matchgroup does NOT exceed that of smaller matchgroup by more than\n"
  "                 <ConfDelta> of -InversionOverlapFilter.\n"
  "              3. Larger matchgroup has same orientation AND aligned label overlap fraction in query does NOT exceed larger of 0.5 and <MaxOverlap> of -SmallDupFilter.\n"
  "           See -OverlapFilter for more stringent matchgroup Filter applied after calling Small SegDups and CutFlip Inversions and before calling any other SVs.\n"
  "      -CutFlipInvDup <DupExp> <MinErr> : Also try to align query region of insertion outlier with <N> or more labels (query larger than reference) with the two neighboring regions of the reference (inverted)\n"
  "                               Expand query interval by <Exp> (see -CutFlip) and reference interval to <DupExp> labels more than query interval beyond outlier in either direction  [Default : OFF or 0]\n"
  "                               The insertion outlier must have a minimum value for (QrySize - RefSize)/(QrySize + RefSize) of <MinErr> (typical value 0.8) \n"
  "    -InversionNormalPV <PV> [<PV5>] : Normal (Non-Inverted) matchgroup's Pvalue must be better than <PV> to call an inversion or inverted duplicate [Default : OFF or 1 1e-15]\n"
  "    -OverlapFilter <QryRatio> <ConfDelta> <QryLabRatio> ]: Filter out smaller matchgroups that are overlapped on query by larger matchgroup by at least <QryRatio> fraction in kb AND <QryLabRatio>\n"
  "                             fraction in Labels AND larger matchgroup confidence exceeds smaller matchgroup by <ConfDelta> [Default : 0.7 0.0 0.0]\n"
  "                             Also if matchgroup is already used in an inversion or duplication call as smaller matchgroup it can only be used subsequently in non-indel SV calls as the larger martchgroup\n"
  "       -OverlapFilterSameRef <0/1> : Restrict -OverlapFilter to apply only when both matchgroups share the same reference contig. [Default: OFF or 0]\n"
  "       -InversionOverlapFilter <InvQryRatio> <InvConfDelta> [ <InvQryLabRatio>] : Modify -Overlapfilter by keeping matchgroups overlapped in query by larger matchgroups in opposite orientation (same ref) or\n"
  "                             another ref, unless Query overlap also exceeeds <InvQryRatio> fraction (in kb) AND <InvQryLabRatio> fraction (in Labels) AND larger matchgroup confidence exceeds that of \n"
  "                             target matchgroup by <InvConfDelta>. These saved overlapped matchgroups can only be used for inversion calls. [Default : OFF or 0]\n"
  "                             If the saved matchgroup is already used in an Inversion call as the smaller matchgroup, it can only be used subsequently for inverted duplication as larger matchgroup\n"
  "       -MergeMG <MaxGap> : If <MaxGap> != 0 : After applying -OverlapFilter merge nearest neighboring non-overlapping matchgroups (same qry & ref) with smallest total gap size in kb, provided\n"
  "                             the gap size in both ref and query does not exceed <MaxGap> kb. Then rerun -OverlapFilter. Finally undo merging of the matchgroups.\n"
  "                             This eliminates duplicative matchgroups in segdup regions and acts as a conservative parsimony filter [Default OFF or 0]\n"
  "       NOTE : The -OverlapFilter is applied AFTER checking for Small SegDups and CutFlip inversions (with smaller matchgroup fully overlapped in both query and reference).\n"
  "              See <OverlapFilter> option of -CutFlipFilter for less stringent matchgroup filter applied before checking for Small SegDups and CutFlip inversions.\n"
  "    -SmallDupFilter <MaxOverlap> <MinOverlap> <MaxMiss> <pvalueBC>\n"
  "                             Call small Duplications based on overlap of small matchgroup with internal outlier insertion (on query) of larger matchgroup in same orientation, but does not overlap reference internal outlier\n"
  "                             Filter out small Duplications if overlap on query includes aligned labels in the non gap region that exceed <MaxOverlap> times total aligned labels of the smaller matchgroup\n"
  "                             Filter out small Duplications if overlap on reference does NOT include aligned labels in the non gap region that exceed <MinOverlap> times total Aligned labels of smaller matchgroup\n"
  "                             Filter out small Duplications if the reference region includes aligned labels in the smaller matchgroup that do NOT correspond to aligned labels of the larger matchgroup and their number\n"
  "                                   exceeds <MaxMiss> times the number of aligned labels of the smaller matchgroup.\n"
  "                             Filter out small Duplications if the pvalue of the smaller matchgroup does not satisfy <pvalueBC> after applying Bonferroni correction based on number of labels seperating the two duplicate\n"
  "                                regions in the query. [Default : 0.2 0.8 0.2 1e-4]\n"
  "    -MATCHGROUP_TRIM <N>   : If N > 0 : Enable trimming (including complete deletion) of matchgroups that are fully overlapped in query :\n"
  "                             N = 1 : Only if fully overlapped on both query and reference AND both matchgroups have same reference contig\n"
  "                             N = 2 : Even if NOT overlapped on reference.\n"
  "                             N = 3 : Even if matchgroups have different reference contigs.  [Default: N = 3]\n"
  "                      NOTE : For N <= 2 even with different reference contigs, trimming still occurs if smaller MG overlaps outlier of larger MG : only the overlapped part is trimmed\n"
  "    -MATCHGROUP_PARTIAL_TRIM <N>  : If N > 0 : Enable trimming of matchgroups that are partially overlapped in query, trimming smaller matchgroup\n"
  "                                    N = 1 : Only if both matchgroups have same reference contig\n"
  "                                    N = 2 : Even if matchgroups have different reference contigs. [Default: N = 1]\n"
  "    -RefSplit <Psplit> <N> [ <Rsplit> ] : Split each alignment (but not maps) at internal regions with score < log(Psplit/BC) and at least N site intervals on either map. [Default : OFF or 0.0]\n"
  "                                         Internal regions can consist of any internal alignment interval with -ve score < log(Psplit) provided no subinterval with score > log(1/Rsplit) is included\n"
  "                                            and the first and last alignment subintervals have -ve scores. If not specified, Rsplit is same value as Psplit. \n"
  "                                            Identified internal regions with -ve score < log(Psplit) but less than N sites intervals are classified as new outliers : -indel will output these as Indels\n"
  "                                         NOTE : BC is the Bonferroni correction for Internal region size, if -outlierBC was specified\n"
  "               NOTE: If MultiMatch is used with -hashgen, must use -hashMultiMatch\n"
  "    -RefSplitStitch <0/1/2> : If >=1 : After applying -RefSplit, check if any remaining matchgroups can be stitched together (end to end). [Default : OFF or 0]\n"
  "                              If >=2 : Also saves partial matchgroups that do not satisfy -S -T -A -L thresholds until after trying to stitch them together, provided they have a +ve score AND were originally part of a larger\n"
  "                            matchgroup that satisfied the -S -T -A -L thresholds\n"
  "    -MultiMatchesRev <0/1> : Run -MultiMatches twice with reference in both orientations, since some matchgroups only show up in one orientation [Default : OFF or 0]\n"
  "    -MultiMatchesTotScore <N> <D> : If N >=1 : Resolve overlapped alignments to maximize the total score of unique parts by avoiding internal outliers, at the expense of losing the longer alignments\n"
  "                                    If N >=2 : Do this before applying -MultiMatchesDelta filter, in case alignments have multiple outliers whose sizes cancel each other so the ends have only a small offset.\n"
  "                                                 This takes extra time, so only do this if at least one outlier is <D> kb or larger (should match size of smallest tandem duplication to be detected) [Default OFF or 0]\n"
  "                                                 Reducing MultiMatchesDelta to <D> instead is not only slower (generating more shifted repeat regions) but will not guarantee finding duplications above <D> kb\n"
  "                                    If N >=3 : Modify Duplicate Alignment checking after -CutFlip alignments exist to reduce MultiMatchesDelta to <D> kb, when a CutFlip alignment is involved, AND require alignment\n"
  "                                                 offsets to differ by less than MultiMatchesDelta at both ends of the alignment before discarding the lower scoring alignment\n"
  "    -MultiMatchesDupes <S> <R> <I> : Eliminate overlap regions between matchgroups, provided score of shared region is at least S AND at least R times total score. [Default 5.0 0.5 5.0]\n"
  "                        If MultiMatchesTotScore, allow the higher scoring matchgroupto be broken, provided this improves the total score by at least I (vs breaking the lower scoring matchgroup)\n"
  "-partial <bin> <numbins> [skipmaps] : Obsolete\n"
  "-K <K> : Set log-likelihood Bias (sizing error allowance) multiplier. [Default 2.2]\n"
  "-KL <KL> : Set log-likelihood Bias (error allowance) per kb multipler. [Default 2.1]\n"
  "-KF <KF> : Set log-likelihood Bias (error allowance) per aligned interval multipler. [Default 2.1]\n"
  "-G <G> : Set misaligned penalty scaling factor OR Nguyen Gamma [Default 1.0]\n"
  "-Beta <B> : Set Beta value [Default 1.0]\n"
  "-KRef <K> : Set log-likelihood Bias (sizing error allowance) multiplier during refinement. [Default 2.2]\n"
  "-KLRef <KL> : Set log-likelihood Bias (error allowance) per kb multipler during refinement. [Default 2.1]\n"
  "-KFRef <KF> : Set log-likelihood Bias (error allowance) per aligned interval multipler during refinement. [Default 2.1]\n"
  "-S <score> : only use alignments with score above this threshold [Default 0.0]\n"
  "-T <Pvalue> : only use alignments with Pvalue below this threshold [Default 0.0001]\n"
  "-A <Aligned-Sites> : only use alignments with at least this many sites [Default 8]\n"
  "-L <Aligned-Length> : only use alignments with at least this length (in kb) aligned [Default 0]. NOTE: Does not yet apply to pairwise alignment\n"
  "-E <N> : only keep alignments with no more than N EndOutliers. See alternate option -Erefine. NOTE: Does not apply to pairwise alignment. [Default OFF or 2]\n"
  "-I <D> [<M>]: only use alignments that have no internal outliers \"larger than D kb OR with M or more misaligned labels\". NOTE: Does not apply to pairwise alignment. [Default OFF or 1000 1000]\n"
  "-F <F> : During autoNoise : only keep alignments where the fraction in kb that is aligned (excluding endoutliers) exceeds F (see also -minoverlap). [Default OFF or 0]\n"
  "-D <D> : During autoNoise : only keep alignments with aligned label density below D labels per 100kb. [Default OFF or 0]\n"
  "-TE <Pvalue> : less stringent -T threshold for alignments that extend beyond reference contig ends (requires -extend) [Default: same as -T]\n"
  "-AE <aligned-sites> : less stringent -A threshold for alignments that extend beyond reference contig ends (requires -extend) [Default: same as -A]\n"
  "-LE <alignment-length> : less stringent -L threshold for alignments that extend beyond reference contig ends (requires -extend) [Default: same as -L]\n"
  "-TS <Pvalue> : less stringent -T threshold for -extsplit alignments, both extending AND non-extending (requires -extsplit) [Default: same as -TE]\n"
  "-AS <aligned-sites> : less stringent -A threshold for -extsplit alignments, both extending AND non-extending  (requires -extsplit) [Default: same as -AE]\n"
  "-LS <alignment-length> : less stringent -L threshold for -extsplit alignments,both extending AND non-extending (requires -extsplit) [Default: same as -LE]\n"
  "-Erefine <N> <wt> [<L> <N> <T>]: similar to -E <N> but discarded maps are retained with weight multiplied by <wt> during refinement, but restored when computing -Haplotype, -ChimQuality, -contigsplit\n"
  "                    If L & N are specified, only EndOutliers of at least L kb and N unaligned labels (in both ref & qry) are counted (NOTE: small EndOutliers can happen due to sizing errors near end)\n"
  "                    If T == 1: weights are NOT restored when computing -Haplotype, -ChimQuality, -contigsplit (see also effect on -Irefine) [Default OFF or 2 0 0 0 0]\n"
  "                    NOTE that -extendWT with -BestRefWT also reduces weight of most maps with endoutlier in extension direction, but only if they have a better alignment elsewhere.\n"
  "-Irefine <D> <M> <wt> : Similar to -I <D>, but discarded maps are retained with weight multipled by wt during refinement, but restored when computing -Haplotype, ChimQuality, -contigsplit\n"
  "                    If -Erefine's T == 1: weights are NOT restored when computing -Haplotype, -ChimQuality, -contigsplit (see -Erefine) [NOT YET IMPLEMENTED][Default OFF or 1000 1]\n"
  "-M <R1> [<R2>] : Repeat reference alignments with error parameter updates <R1> times and save results in <prefix>.err [Default 1 1]\n"
  "                               If <R2> is specified perform <R1> x <R2>  parameter updates, redoing hashtable every <R1> updates : only useful if -hashgen is specified\n"
  "     -relerr <relerr>  : Terminate if last -M update did not change any parameter by more than this fraction [Default: 1e-5]\n"
  "     -Msave <N> : Instead of saving last set of error parameters in .errbin save the values after the <N>'th -M iteration.[Default: <R1>]\n"
  "                      If <R2> was specified, the <N>th iteration for the last hashtable iteration is saved\n"
  "     -outlierRate : output outlier and Endoutlier Rate in .err file : these are per label interval rates (NOT per molecule)\n"
  "     -perLabelScore : output score per label in .err file\n"
  "     -LabelDensityType <0/1> :  If 0, output reference based label density from aligned intervals\n"
  "                            If 1, output query based label density from aligned intervals [Default 0]\n"
  "  -resbias <xmax> <N> : Compute correction for sizing bias caused by resolution limits for small site intervals below <xmax> using <N> segment piecewise linear regression. [Default 0 0]\n"
  "                      NOTE : previously computed correction must be read in from an .errbin file with -readparameters\n"
  "  -ScanScaling <N> [ <minSNR> <minlen> <Min>]: Apply per scan scaling correction (only works with BNX 1.0 data). <N> specifies the number of iterations of scaling estimation per -M iteration [Default 0]\n"
  "                   Rescaled BNX data will be output to <prefix>_rescaled.bnx (NOTE: If -subset is used, all map will be rescaled and output, but only the subset will be used to estimate scaling)\n"
  "                   -minsites, -maxsites, -minlen and -maxlen will be applied just before output to <prefix>_rescale.bnx and -minlen,-maxlen are only applied with a 15% margin on input\n"
  "                   If <minSNR> and <minlen> are specified, -minSNR and -minlen value is changed to this value just before output of _rescaled.bnx (before applying -minsites, -maxsites, -minlen and -maxlen)\n"
  "                   A smaller <minSNR> or <minlen> value allows the _rescaled.bnx to be used as input to a subsequent parameter estimate (autoNoise) step with -minSNRestimate\n"
  "                   If Min > 0 is specified, cohort scaling is NOT corrected unless at least Min maps in that cohort aligned\n"
  "      -ScanFilter <M> <Rm> <Rb> <Re>: Filter out scans(cohorts) with less than M maps aligned OR with Mapping rate (by KB) below Rm times average Mapping rate (by KB)\n"
  "                      OR with bppSD above Rb times average value OR with Endoutlier % above Re tunes average value . [Default : OFF or 0]\n"
  "      -MaxCoverage <C> <L> : Reduce output coverage of <prefix>_rescale.bnx to no more then C (x -ref size) by increasing minlen up to L kb and then subsampling to reduce coverage\n"
  "  -resEstimate [<N>] : Estimate -res and -resSD error parameters [Default OFF]\n"
  "  -minSNRestimate <min> <max> [ <step> <N> <F>] : Estimate -minSNR values within range specified (For 2-color data same range is used for both colors) [Default OFF]\n"
  "                                   If specified, <step> is used as a stepsize to perform an initial global scan using no more than <N> molecules [Default OFF]\n"
  "                                   If specified, <F> is used to raise the initial value of minSNR so that molecule label interval is no more than <F> times the reference label interval size [Default 0]\n"
  "         -bppScan <Stepsize> <Steps> <MinMaps> : If Initial -bpp (possibly from -readparameters) results in fewer than <MinMaps> aligned maps (out of <N>), try alternate values:\n"
  "                    Try changing initial bpp value by <Stepsize> (in both directions) up to <Steps> times. Stop at the first value that produces <MinMaps> aligned maps. [Default : OFF]\n"
  "                    NOTE: -bppScan only applies with -minSNRestimate and is applied after estimating minSNR with Label density <F> but before applying the minSNR scan\n" 
  "                          -bppScan changes to bpp are undone if -ScaleDelta steps exceeds 3, so that bpp remains centered around original value\n"
  "  -MapRate <MapRate> [minT [maxS]] : If needed, adjust -T down (or -S up) so the specified MapRate is not exceeded. Works best with a lose -T (or and -S) value [Default : OFF]\n"
  "                     -BestRefPV 1 uses -T threshold, -BestRefPV 0 uses -S threshold\n"
  "                     If specified minT and maxS specify a limit on -T and -S adjustment : useful to handle unexpectedly good data sets [Default : OFF]\n"
  "  -Mfast <0,1 or 2> : speed up repeated alignment iterations (with -M) by assuming they don't change much from the first iteration (0 = none, 1 = conservative, 2 = agressive)[Default = 0]\n"
  "  -relaxBnds <0/1> : Remove -MinSR -MaxSR etc bounds or error parameters in the last -M iteration so final estimates are accurate but avoiding converging to low mapping rate [Default : OFF or 0]\n"
  "  -biaswtDefer <0/1> : Force -biaswt to 0 until final -M iteration so final Mapping rate estimate is accurate but avoiding converging to low mapping & error rates [Default : OFF or 0]\n"
  "-readparameters <file.errbin> : read in error parameters from specified file generated by a previous -M run.\n"
  "            NOTE: error values in <file.errbin> silently override any command line error values before or after this option (except for -mres -mresSD, unless -mres_override 0 is specified)\n"
  "-mres_override <0/1> : Value of 1 means -mres -mresSD commandline values override values from -readparameters [Default : 1]\n"
  "-maxmem <Gbytes> : avoid using more than specified amount of memory by reducing number of threads etc [ Default 128]\n"
  "    -maxmemIncrease <ResMem> : If SysMem > 0 : Increase maxmem to system memory size - ResMem (in GB) [Default OFF or 0]\n"
  "-maxvirtmem <Gbytes> : avoid using more than specified amount of virtual memory by reducing number of threads etc. 0 means unlimited [ Default : 3x -maxmem]\n"
;
  char *usage2 = (char *)
  "-refine <level> : refine the input reference map [Default OFF or 0]:\n"
  "                  level = 1 : Just apply -rres -EndTrim and -ChimQuality (fast, with no changes to map, unless -rres or -EndTrim are used)\n"
  "                  level = 2 : Add/delete sites and optimize site interval sizes(slow)\n"
  "                  level = 3 : Same as level=2 unless -extonly is used : in that cases keeps all maps for -ChimQuality and .xmap output, even for regions not refined by -extonly\n"
  "  -RefineOverwrite <0/1> : If 0, skip contigs which already exist in refined form; otherwise overwrite previously refined contig [Default 1 (0 with -Haplotype)]\n"
  "  -UseEndoutliers <0/1> : If 1, Use sites in endoutlier region of molecules to hypothesize where sites should be added in consensus map [Default 1]\n"
  "  -Mprobeval [ <res> <rel> ] : Specify <res> = Site interval accuracy (in kb) [Default 0.01]. Using 0.001 kb is over 2x slower than 0.01 kb.\n"
  "                               If specifield <rel> specifies the Site interval accuracy as a fraction (applies for intervals where that is less than <res>). [Default 0.01 0.001]\n"
  "  -RefineStep <R1> <R2> [ <R3> <R4> <N> <R5> <R6> <R7>]: Use R1 as initial sizing change step size during 1st refinement stage, and R2 during subsequent stages (or all stages of RefineA). [Default 0.10 0.20]\n"
  "                           Use R3 & R4 instead of R1 and R2 (respectively) after initial sizing change and limit number of steps (size values) to N. This is only done if previous round used more\n"
  "                           than 2 values AND best value was at range end OR either of neighboring interval changed by more than <R7>. [Default R3=0.10 R4=0.20 N=3 R7=0.10]\n"
  "                           Reducing R1 & R2 can improve accuracy of consensus map and ability to handle multimodal maxima at expense of larger run time\n"
  "                           Note that the number of steps (size values) with R3 & R4 is limited to N, so just reducing R3 & R4 trades off search space for precision, without changing run time\n"
  "                           The minimum step sizes with R1 or R3 are R5 (kb) and with R2 or R4 are R6 (kb) [Default 0.2 kb]\n"
  "                           With -MprobevalTest, the values are forced to 0.02,0.02,0.05,0.05,14,0.02,0.05,0.02 (along with many other changes)\n"
  "                           The values set by -MprobevalTest can be overridden by this option on the commandline\n"
  "  -RefineSizeRange <Low1> <High1> <Low2> <High2> : Use ratios <Low1> and <High1> during initizing sizing change as scan range relative to current interval size\n"
  "                            Use ratios <Low2> and <High2> during subsequent stages (or all stages of RefineA). [Default: 0.01 2.0 0.7 1.4]\n"
  "  -SiteMerge <0/1> : During refinement periodically try to merge all nearby labels, if that improves likelihood [Default : 1]\n"
  "  -RefineFix <0/1/2> : turn on adjustments and bug fixes in Refinement introduced after revision 4535 to support newer -Haplotype suboptions [Default: OFF or 0 unless -HapIndelMerge]\n"
  "                     This option is also recommended when using the refinement Switch options like -outlierTypeSwitch etc without -HapIndelMerge\n"
  "  -RefineErr <E1> <E2> <E3> : Allow reduced likelihood accuracy specified by log-likelihood errors of E1, E2 and E3 when computing endoutlier likelihood, total alignment likelihood and\n"
  "                              sub-alignment likelihood respectively. Increasing E1,E2 or E3 will speed up refinement with more predictable loss in accuracy by reducing deltaX, deltaY\n"
  "                              and RefineRange dynamically, compared to using fixed reduced values [Default : 1e-8 0 0]\n"
  "                              NOTE: Any non-zero value of E3 will incur a small run-time overhead\n"
  "  -TB <PV1> ... <PVn> : A sequence of Pvalues that is used to downweight maps (eg in the extension region) after preliminary refinement [Default : None](requires -extend)\n"
  "                        Maps are re-weighted based on actual map Pvalue : wt = TB/(TB + Pvalue)\n"
  "       [-TBmult <PVmult>] : Interpolate values between <PV1> and <PVn> in ratio steps of <PVmult> : Values between <PV1> and <PVn> are ignored:\n"
  "                            Example : -TB 1e-7 2e-12 -TBmult 0.1 is equivalent to -TB 1e-7 1e-8 1e-9 1e-10 1e-11 2e-12  [Default : OFF or 0]\n"
  "                        NOTE : Not used during Haplotype Refinement, instead the pre-refinement mapWT is used. See -LRbiasSwitch and -LRbias to adjust LR value threshold instead\n"
  "  -LRbias <bias> : A minimum LR (Likelihood Ratio) value per map. This value is added to the actual LR to reduce the influence of maps with low LR values. [Default 0]\n" 
  "                 If -extend 2 then sites are added/deleted before applying the -LRbias option, since the initial LR in the extension region will be very low\n"
  "                 Also reduces the final weight of maps that score below <bias> by wt = bias/(bias + LR) : Used in the coverage computation by -contigsplit -ChimQuality\n"
  "                 NOTE: if both -LRbias and -TB are specified final map weight is product of both map weights\n"
  "  -deltaXRef <N> : change deltaX during refinement to N. Useful since -outlierExtend does NOT apply during refinement [Default : deltaX]\n"
  "  -deltaYRef <N> : change deltaY during refinement to N. Useful since -outlierExtend does NOT apply during refinement [Default : deltaY]\n"
  "    -deltaExtXRef <M> : If M > deltaX, allow query map intervals to be as large as M, whenever reference map intervals are <= deltaY [Default : deltaX]\n"
  "    -deltaExtYRef <M> : If M > deltaY, allow reference map intervals to be as large as M, whenever query map intervals are <= deltaX [Default : deltaY]\n"
  "                      NOTE : Using larger values for deltaExtXRef, deltaExtYRef is faster than using the same values for deltaXRef,deltaYRef and \n"
  "                             can handle the same size Haplotype insertion or deletion (but not a Haplotype Inversion)\n"
  "    -deltaEndRef <U> : If U > 0 : allow non-outlier ends with up to deltaExtXRef (or deltaExtYRef) + U - 1) Labels, instead of deltaXRef,deltaYRef, during refinement [Default : 0]\n"
  "  -RefineRange <R> <R2> : R = largest change in map alignment (in labels) allowed during each Refinement step [Default 8]\n"
  "                        Reducing this value to 4 will dramatically speed up refinement but may sometimes produce a worse consensus map\n"
  "                        R2 is a larger value used only for initial alignment and periodic re-alignment during refinement (Using same value as <R> will be more stable) [Default 8]\n"
  "  -RefineXRange <0/1/2> : If 1, use extrapolation of alignments based on both RefineRange and deltaExtXRef (slightly more stable and slower refinement) [Default: OFF or 0]\n"
  "  -LocalRange <0/1> : If 1, Use label interval size of current molecules instead of global average value : will speed up refinement for high label density regions. [Default 0]\n"
  "  -LocalTheta <0/1> [MaxTheta] : If 1, Use local estimate of consensus label interval size instead of genome wide average value : will speed up refinement for high label density regions. \n"
  "                                 If MaxTheta is specified, reduce estimate of label interval size to no more than MaxTheta during non-Haplotype refinement [Default 0]\n"
  "  -ScanRange <X> : During non-hap refinement : Multiple of molecule label interval to scan ahead/back for improvements before confirming them (picking locally best improvements) [Default 3.0]\n"
  "  -HScanRange <X> : During haplotypeing refinement : Multiple of molecule label interval to scan ahead/back for improvements before confirming them (picking locally best improvements) [Default 3.0]\n"
  "  -RefineAlignExt <0/1> : If 1, force extrapolation of alignments to reach ends of molecules (slightly more stable and slower refinement) [Default: OFF or 0]\n"
  "  -RefineFillMap <0/1> : If 1, use older method for extrapolation of alignments during non-haplotype refinement [Default: ON or 1]\n"
  "  -RefineFillHmap <0/1> : If 1, use older method for extrapolation of alignments during haplotype refinement [Default: ON or 1]\n"
  "  -endoutlierRef <E1> <E2> ... <En> : -endoutlier Pvalue sequence used during successive phases of refinement (if -extend 2) to converge to the consensus estimate of the extension region\n"
  "                        This can help when the extension consists of a mixed population by allowing the lower scoring extensions to become outliers and stop contributing to the consensus\n"
  "                        Using at least two different values causes refinement to first try adding/deleting labels (recommended during extension refinement), instead of adjusting sizing first\n"
  "  -outlierRef <P1> <P2> ... <Pn> : -outlier Pvalue sequence used during successive phases of refinement (if -extend 2) to converge to the consensus estimate of the extension region\n"
  "                        Defaults to the same sequence as for -endoutlierRef. Length of sequence will be matched to -endoutlierRef by deleting the first few values OR adding values\n"
  "                         at the begining matching those of -endoutlierRef \n"
  "  -outlierTypeRef <0/1> : Modifies outlier penalty during Refinement [Default 1]\n"
  "                          If 0, don't penalize misaligned sites inside an outlier\n"
  "                          If 1, Penalize misaligned sites, even inside an outlier\n"
  "  -outlierNormRef <0/1/2> : Modifies outlier penalty normalization factor during Refinement [Default 1]\n"
  "                          If 1, don't adjust outlier penalty normalization factor based on -outlierLambda value \n"
  "                          If 2, don't adjust outlier penalty normalization factor based on -outlierMax value \n"
  "  -outlierLambdaSwitch [<Label> <Size>] : suppress outlierLambda (or force it to <Label>) when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "                                          restore outlierLambda (or force it to <Size>) when adjusting sizing (or checking for Haplotype indels) during refinement\n"
  "  -outlierTypeSwitch [<Label> <Size>] : restore outlierTypeRef to 1 (or force it to <Label>) when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "			                   change outlierTypeRef to 0 (or force it to <Size>) when adjusting sizing (or checking for Haplotype indels) during refinement\n"
  "  -outlierSwitch <Label> <Size> : change outlier Pvalue to <Label> when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "                                  change outlier Pvalue to <Size> when adjusting sizing (or check for Haplotype indels) during refinement\n"
  "                                  NOTE : May be modified by -RefineSwitchFade\n"
  "  -endoutlierSwitch <Label> <Size> : change endoutlier Pvalue to <Label> when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "                                     change endoutlier Pvalue to <Size> when adjusting sizing (or check for Haplotype indels) during refinement\n"
  "  -LRbiasSwitch <Label> <Size> : change LRbias value to <Label> when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "                                 change LRbias value to <Size> when adjusting sizing (or checking for Haplotype indels) during refinement\n"
  "  -MinSFSwitch <Label> <Size> [ <Liter> <Siter> ]: increase sf value to at least <Label>  when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "                                                   increase sf value to at least <Size> when adjusting sizing (or checking for Haplotype Indels) during refinement\n"  
  "                                If Liter > 0 : Limit sf adjustment when adding/deleting labels to the first Liter such iterations (otherwise see -RefineSwitchFade)\n"
  "                                If Siter > 0 : Limit sf adjustment when adjusting sizing to the first Siter such iterations (otherwise see -RefineSwitchFade)\n" 
  "  -MinSRSwitch <Label> <Size> [ <Liter> <Siter> ] : increase sr value to at least <Label>  when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "                                                    increase sr value to at least <Size> when adjusting sizing (or checking for Haplotype Indels) during refinement\n"  
  "                                If Liter > 0 : Limit sr adjustment when adding/deleting labels to the first Liter such iterations (otherwise see -RefineSwitchFade)\n"
  "                                If Siter > 0 : Limit sr adjustment when adjusting sizing to the first Siter such iterations (otherwise see -RefineSwitchFade)\n" 
  "  -ViterbiSwitch <Label> <Size> [ <Liter> <Siter> ] : Use Viterbi & Bayesian interpolated scoring with Bayesian weight <Label> when adding/deleting labels during Haplotype refinement\n"
  "                                                      Use Viterbi & Bayesian interpolated scoring with Bayesian weight <Size> when adjusting sizing during Haplotype refinement\n"
  "                                If Liter > 0 : Revert to Bayesian weight = 1 after Liter iterations of adding/deleting labels (also reverts eventually with -RefineSwitchFade)\n"
  "                                If Siter > 0 : Revert to Bayesian weight = 1 after Siter iterations of adjusting sizing (also reverts eventually with -RefineSwitchFade) [Default: OFF]\n" 
  "  -IndelScore <outlierLambda> <outlierType> <outlier> <outlierEnd> <LRbias> [ <MinSF> <MinSR> ]: Modify values of -outlierLambda, -outlierType, -outlier, -outlierEnd, -LRbias during final Haplotype Indel\n"
  "                                 scoring and during Indel deletion check. [Default OFF or 0,-1,-1,-1,0]\n"
  "  -endoutlierFinal <Pvalue> : Increase Pvalue to this value after refinement for final alignment and coverage profile computation : useful to reveal chimeric sites with -contigsplit [Default 0.0001]\n"
  "  -extonly [ <L> <L2>] : Only refine the extension region + a region of L kb at each end of original contig[Default OFF, 100kb] (requires -extend 2)\n"
  "                     If L2 < L : During initial (agressive) refinement (see -endoutlierRef) restrict refinement to extension region + a region of L2 kb at each end of original contig\n"
  "                     Use -refine 3 to keep alignment information of all maps in .xmap output (and -ChimQuality scores), even those maps that aligned to the non-extension region\n"
  "  -extsplit <C> <F> <N> <NL> <M> <L> <S> [ <LS> <E> <LE> ] <CntFile> : Clone part (see -CloneLimitLen for size) of contig on one side of any breakpoint site where at least max(C, F * coverage)\n"
  "                              endoutliers in the direction away from the cloned part are present. Don't count endoutliers shorter than L (kb) excluding unlabeled end. [Default: OFF or 0.0]\n"
  "                              Don't count endoutliers whose alignments ends more than N sites AND NL kb before target site or more than M sites after the target site.\n"
  "                              Only consider the largest cluster of endoutliers within max(S,N,M) labels (but if -extsplitSpacingFix >= 1, limit scan to LS times average label interval size, provided\n"
  "                                  peak has been verified OR -extsplitSpacingFix 2) : LS defaults to max(S,N,M).\n"
  "                              If -extsplitSpacingFix >= 1 : Only consider break point sites at least M+E-1 labels AND LE times average label interval size from original contig ends:\n"
  "                                  E defaults to S and LE defaults to LS.\n" 
  "                              File named <CntFile> contains a 63-bit integer that is added to a contig ID for each new split contig (same as with -contigsplit option). If <CntFile> is NOT an\n"
  "                              absolute pathname, it is located in the same folder as the -o output file prefix. If <CntFile> does not exist it will be created by looking for the largest\n"
  "                              ID in the group_manifest file in the stage 0 folder. The stage is infered from the -o output file prefix refineB1, extension1_<N> OR refineFinal1\n"
  "                              The -extsplit option is intended to identify and resolve heterozygous breakpoints as well as incorrectly merged segdup regions during extension refinement\n"
  "                              Requires -extend 2 (and optionally -extonly) and -E 2 (but can use -Erefine 0, if desired)\n"
  "                  NOTE : -TE and -T (or -T2 if -MultiMatches) should be set to same value since endoutliers are defined differently by -TE than -extsplit. If less stringent -TS is used both\n"
  "                         -TE and -T are changed to -TS for split contigs. Also recommend using -extsplitNE 1 (the default), so all molecules in split contigs are extending molecules\n"
  "      -extsplitVerifyPeak <0/1> : If 1 : Suppress -extsplit unless the endoutlier coverage values are decreasing on both sides of the breakpoint. [Default 1]\n"
  "      -extsplitSpacingFix <SF> : If SF >= 1 : Suppress -extsplit if within N+E labels (or within LE times average label interval) of contig end.\n"
  "                                              Also limit scanning range for detecting coverage peaks to LS times average label interval, provided peak has been verified.\n"
  "                                 If SF >= 2 : Limit scanning range to LS times average label interval even if peak has NOT been verified. [Default 1]\n"
  "      -splitcnt : Determine contig ids of contig splits by repeatedly adding the value in <CntFile> to the original contig id (this avoids relying on contig locking) [Default OFF]\n"
  "                  Otherwise the next contig id is determined by incrementing the value in <CntFile> by 1\n"
  "      -CloneLimitLen <L> <N> : Limit length of cloned part of contig to no more than L kb or N labels (whichever is greater) and mark the trimmed end as non-extendable [Default OFF]\n"
  "            -CloneTrim <0/1> : If 0 : Do NOT allow non-extension end of split contig to be trimmed by -EndTrim below original length.\n"
  "                               If 1 : Allow non-extension end of split contig to be trimed by -EndTrim below original length (this avoids formation of chimeric contigs) [Default: 1]\n"
  "      -extsplitMinWT <mWT> : Skip alignments with weight below <mWT> while computing endoutlier counts and coverage. Note that -extsplitMinWT is applied before -splitWT [Default 0.01]\n"
  "      -splitFilt <mWT> [ <N> <NL> <M> ] : Skip (end)outliers moved to splits if weight is below <mWT> OR whole alignments ends more than N sites AND NL kb before target site or more than M sites after \n"
  "                 target site. N,NL and M default to -extsplit values (use larger values to get more molecules without wasting time on false splits from larger -extsplit values) [Default OFF or 0]\n"
  "                 NOTE : This option (with mWT >= 0) is required to activate -extendWT even without -extsplit\n"
  "      -extendWT <L> <WT> <N> <NL> [<mWT> <mD> [ <mM> ] <rWT> <mDs> <mMs>]: Undo effect of -BestRefWT on input mapWT by raising it to at least WT for alignments that extend beyond the contig ends by at least L kb, not counting\n"
  "                  the unlabeled molecule ends, and excluding molecules with more than N labels AND more than NL kb in unaligned & overlapped endoutlier region before the contig end.\n"
  "                         For split contigs, max endoutlier misaligned size N and NL do not apply, but see -extsplit N & NL values instead\n"
  "                         For split contigs, minimum extension size L kb does not apply, but see -extsplit L value instead\n"    
  "                  If mWT is specified, do not modify mapWT if original value is below mWT.\n"
  "                  If mD is specified, do not increase mapWT if alignments has internal outlier size difference above mD kb OR at least mM misaligned labels: mD and mM defaults to 1000.\n"
  "                      If rWT < WT is specified, reduce mapWT to rWT if internal outlier size above mD (or >= mM  misaligned labels) was detected\n"
  "                      For split contigs, only remaining alignment is checked for internal outlier and replace mD by mDs and mM by mMs : mDs and mMs default to 1000.\n"
  "                  This option does not affect endoutlier counts and overage and hence deciding which contig splits happen.\n"
  "                  For example -extendWT 20.0 1.0 6 60.0 may increase rate of extension by allowing same molecule to be used on two ends extending towards\n"
  "                  each other, while still allowing -BestRefWT with -splitFilt 0.2 to suppress molecules aligned to the wrong Allele and avoid rediscovering splits in each extension round. [Default OFF]\n"
  "                  NOTE : This option applies to endoutliers moved to the split contig (see -extsplit) if -splitWT with sWT > 0 is used\n"
  "                  NOTE : This option applies even without -extsplit IF -splitFilt with mWT >= 0 is used (but see -extendWTend)\n"
  "                  NOTE : Does NOT apply to alignments filtered by -splitFilt, which is applied before -extendWT\n"
  "        -extendWTdup <0/1> :  If 0, allow only one alignment per molecule to have its mapWT raised by -extendWT [Default: 1]\n"
  "        -extendWTend <0/1> :  If 0, allow only split contig extension alignments to have their mapWT raised by -extendWT, otherwise applies to regular extension as well. [Default: 0]\n"
  "      -splitWT <sWT> : Apply -extendWT (see <L> <WT> <N> <NL> above) to endoutliers moved to the split contig (as extensions). Also increase weight of all alignments to <sWT> while computing endoutlier\n"
  "                  counts and coverage. Using WT = sWT = 1.0 will cause cause all splits to be rediscovered in each extension round, but avoids missing one of the possible phasings of breakpoints\n"
  "                  due to segdupes or het inversions. Using sWT = 0.01 will avoid rediscovering splits in each extension round. [Default OFF or 0]\n"
  "      -keepsplitmaps <0/1> : If 1, keep maps with endoutliers that are used in contig splits in the original contig, but at reduced weight based on -BestRefWT. [Default OFF or 0]\n"
  "      -extTrim <0/1> : If 1, trim back <N> + <M> labels at start of extension refinement of each split contig. [Default: OFF or 0]\n"
  "      -extsplitNE <0/1> [ <L> <E> ] : If 1, include non-endoutliers that align only to the split contig region. If 0, don't include non-endoutliers in split contigs AND reduce -extonly values to L \n"
  "                  (but at least N (or NL kb) + M label intervals, see -extsplit) and reduce -maxExtend by the modified -extonly value, but not less than E or half orig value[Default: 1 0.0 0.0]\n"
  "      -extsplitOutlier <LO> [ <M> <LS> <AS> [ <L1> <L2> <L3> ... <Ln>]] : In addition to endoutliers also count outliers of LO kb or more (size difference) OR with M or more misaligned labels\n"
  "                      with either outlier boundary matching the endoutlier cluster.  [Default: OFF or 0.0 1000]\n"
  "                  The equivalent alignment portion (to one side of outlier) must have at least LS kb and AS labels in aligned region.\n"
  "                  If specified, L1, L2 ... Ln forms a ladder of outlier sizes in kb starting at L0 (L0 < L1 < L2 < .. < Ln) : Outliers are clustered seperately in each pair of consecutive intervals\n"
  "                      starting with the interval L0 .. L2 and ending in the intervals L<n-2> .. Ln and above L<n-1>, with the last interval including endoutliers. There must be at least 2 values (n>=2).\n"
  "                      All internal outliers with M or more misaligned labels are treated as falling in the last interval including endoutliers.\n"
  "          -extsplitOutlierC <C0> <C1> ... <Cm> : Modify minimum number of outliers/endoutiers (-extsplit C) to C0 for all insertion intervals and C1 for deletion interval L0..L2, C2 for L1..L3 etc\n"
  "                                       The original -extsplit C value is still used for all intervals with endoutliers\n"
  "                                       Values Ck with k >= n are ignored. If m < n-1, Additional values C<m+1> .. C<n-1> equal to Cm are added\n"
  "  -refineMinWT <WT> : Skip alignments with weight (see -BestRefWT) below <WT>. Note that -refineMinWT is applied after -extsplit has removed molecules for split extension and -extendWT has modified \n"
  "                      the weights, hence use mWT value in -extendWT instead (For -extsplit can also use -extsplitMinWT & -splitFilt)\n"
  "                      Skipped alignments may be restored for final QualityScore computation (see -ChimQuality).  [Default OFF or 0.0]\n"
  "                      Even small values of WT (eg 0.02) can speed up processing significantly. Larger values (eg 0.2) can prevent already discovered alternate Alleles from being rediscovered in refineFinal1.\n"
  "     -refineWT <bWT> [ <N> <NL> <mD> <mM>] : Increase weight of all alignments that pass -refineMinWT to at least <bWT> during refinement (requires refineMinWT's WT > 0). [Default OFF or 0.0]\n"
  "                      If N, NL and mD are specified do NOT increase weight for alignments with more than N labels AND more than NL kb in any endoutlier or more than mD kb size difference\n"
  "                                OR >= mM unaligned labels in any internal outlier\n"
  "  -maxExtend <L> : limit refinement of each extension region to the first <L> kb (to speed up processing) [Default OFF or 0]\n"
  "  -minoverlap <R> : Ignore extending nanomaps with less than the specified ratio R of their length overlapping the original reference contig. Values above 0.6 prevent chimeric extensions.\n"
  "                    May prevent catching low allelelic fraction extensions or -extsplit. [Default: OFF or 0]\n"
  "  -cres <pixels> <N> [ <F> ] : Check refined consensus map for consecutive sites (2 or more) closer than <pixels> times 0.5 kb and combine some of them if there are not at least \n"
  "                                <N> molecules AND at least ( <F> x coverage ) molecules that have all of the sites [Default : 0 3 0.10]\n"
  "                                NOTE: Sites closer than rres * 0.5kb are always merged before applying -cres\n"
  "        -cresRefine <0/1/-1> : Add a size refinement step between each stage of -cres Site merging (otherwise do this only once after last -cres stage, which is faster) [Default : 1]\n"
  "                               A value of -1 skips any size refinement step, even after the last -cres stage\n"
  "        -cresMaxLPdrop <LP> : Don't merge sites if LP drops by more than <LP> [Default : OFF or -1]\n"
  "        -cresCheckDel <0/1> : When merging labels with -cres, check if deleting either Label scores better than merging them at the midpoint [Default : 0 or OFF]\n"
  "        -cresHapMaxLPdrop <L> <HLP> : Merge opposite Allele SNPs within L kb if LP drops by no more than HLP. Defaults to L = cres * 0.5 and HLP = 0\n"
  "  -id <ID> : Previously used to replace the <ID> value in the refined .cmap output (since .spots input had no ID). It also changed the output file name prefix to PREFIX_contig<ID>\n"
  "             Now with any <ID> value other than 0, the ID from the reference .cmap is used instead for both purposes. With multiple reference contigs, multiple refined .cmap outputs will get\n"
  "             their id from the contig id of the reference map regardless of the <ID> value used, but only a single .xmap (and _r.cmap, _q.cmap && .map) is output unless each reference contig originates\n"
  "             from a seperate -ref input file AND <ID> is non-zero (see option -mapped-unsplit to override the choice of single or multiple .xmap files).\n"
  "             The non-zero <ID> value is also used to add a _id<ID> suffix to the .err and .hashbin and _hash.bnx output files to make them unique. [Default: 0]\n"
  "  -RTHETA_FIX <0/1> : Set to 0 to better match endoutliers before and after refinement (1 reduces score of endoutliers during refinement) [Default : ON or  1]\n"
  "  -MinMapwt <wt> [ <Uwt> ]: In refined XMAP only show molecule maps with map weight above specified threshold, to filter out poorly aligned maps [Default 0.30 0.0]\n"
  "                    The map weight is the relative contribution of the molecule map to the refined consensus map between 0 and 1 (see -BestRefWT -LRbias and -TB for how weight is computed)\n"
  "                    If <Uwt> is specified the Unrefined (normal) XMAP and -mapped BNX output is also filtered for molecules with map weight above <Uwt>, based on -BestRefWT only\n"
  "  -contigsplit <R1> <R2> <CntFile> [<R3>]: split input reference map at internal site interval if coverage drops below <R1> x median coverage AND below <R2> x median trimmed coverage (used with -refine)[Default OFF]\n"
  "                                           <R3> is the minimum ratio of the trimmed coverage at nearby peak levels vs the split location (otherwise the contig is not split) [Default 1.5]\n"
  "                                           If -TrimNormMed <N> is specified, R1 is ignored and threshold is set to R2 x N\n"
  "      -splitrev [ <N> ] : Use new revision <N> of -contigsplit [Default : 0]: NOTE -splitrev without an arg defaults to -splitrev 1\n"
  "      -splitcnt : Determine contig ids of contig splits by repeatedly adding the value in <CntFile> to the original contig id (this avoids relying on contig locking) [Default OFF]\n"
  "      -splitsite : Apply -contigsplit only based on coverage at site locations (and not between site locations) [Default OFF]\n"
  "      -splitFragile <D> <R> <E> : Also break contig if Fragile Left and Right scores differ by at least <D> AND their ratio exceeds <R>. If <E> is 1, mark end with larger Fragile score as non-extendable (chromosome end).[Default OFF or 0]\n"
  "                                  Requires -ChimQuality 2. If -splitsite, the check is only applied to at site locations.\n"
  "  -CovTrim <sites> : minimum flanking aligned sites required to anchor a map on both sides of the target site to count that map towards ChimQuality (type N1) at the target site [Default 3]\n"
  "                   This value should be adjusted down if minlen is reduced (below 14 average site intervals) and can be adjusted up if minlen/minsites is increased (relative to average site interval)\n"
  "     -CovTrimLen <kb> : Similar to -CovTrim but expressed as the minimum kb to be aligned on both sides of the target site : must be satisfied in addition to -CovTrim [Default 0]\n"
    "     -TrimNorm [<N> [<L>]]: caused the trimmed coverage of each site to be normalized by the number of molecules available that extend by at least CovTrimLen to both sides of the site [Default OFF or -1.0]\n"
  "                      If <N> greater than 0 is specified, a pseudo-count of <N> is added to both the trimmed coverage and the normalization count\n"
  "                      If <L> greater than 0 is specified, use the largest normalization value within L kb (increases sensitivity, particularly when small segdups bridge chimeric contig)\n"
  "     -TrimNormChim <N> [ <M> ] : Try to suppresss counting chimeric molecules with -TrimNorm : don't count molecules whose alignment ends more than N sites (in consensus) before the target site [Default OFF or -1]\n"
  "             If <M> is specified, also don't count molecules whose alignment ends more than M sites (in consensus) after the target site, but alignment does NOT extend CovTrimLen kb AND CovTrim sites beyond target site [Default OFF or -1]\n"
  "             NOTE : N is the maximum number of false negative labels tolerated in alignments of non-chimeric molecules near the target site in the consensus\n"
  "             NOTE : M is the maximum number of false positive aligned labels tolerated in alignments of non-chimeric molecules near the target site in the consensus. Should be less than CovTrim\n"
  "     -TrimNormMin <N> [ <C> <M> <R>]: Modify -TrimNorm to require a ChimQuality coverage (type N1) of at least <N>, otherwise force ChimQuality to 0 [Default 0.0 0.0]\n"
  "                        If C > 0 is specified, only apply TrimNormMin if (actual) coverage is below C\n"
  "                        If M is specified don't apply near contig ends until actual coverage reaches at least M (M defaults to N)\n"
  "                        If R is specified don't apply unless N1 is below R times the highest value of N1 beyond the current label (applied in both directions)\n"
  "     -TrimNormEnd <N> : Modify -TrimNorm to exclude from trimmed coverage molecules with endoutliers with at least <N> unaligned sites, even if nowhere near the target site [Default OFF]\n"
  "     -TrimNormMed <N> : Modify -TrimNorm to use <N> instead of median of trimmed coverage to compute threshold using <R2> of -contigsplit (ignores <R1>)[Default OFF]\n"
  "     -TrimOutlier <0/1> [<D> <M>]: For molecule alignments with outliers treat each aligned region between outliers as a seperate molecule when computing N1 [Default OFF or 0 0 1000]\n"
  "                                   Ignore outlier unless sizing error is at least D kb OR there are at least M misaligned labels\n"
  "                                   If -outlierTypeRef 1 is specified, then any aligned interval with at least M misaligned labels is treated as an outlier\n"
  "     -ReplaceCov : causes the trimmed coverage to be output in place of the actual coverage as the consensus contig coverage [Default OFF]\n"
  "  -CovSegDup <endK> <endN> <L1> <L2> ... : Also output for each label counts of number of molecules that align across a central region of <L1> <L2> ... kb plus a flanking region of at least <endK> kb\n"
  "       AND <endN> labels. The results counts will be output in the CMAP as optional columsn with name SD<L1> SD<L2> etc. Note that the N1 count computed by -CovTrim <N> -CovTrimLen <K> is the same\n"
  "       as the count that would be computed by -CovSegDup <N> <K> 0. [Default OFF]\n"
  "  -FragileQuality <L> <FN> <FP> : Compute Left/Right Fragile Quality scores at each consensus (target) label by counting aligned molecules with Left/Right end within <L> kb and\n"
  "                                  no more than <FN> false negative labels before target label and no more than <FP> False Positive aligned labels beyond target label\n"
  "  -OutlierQuality <0/1/2> : Compute OutlierQuality at each consensus label interval as fraction of aligned molecules with internal outliers across that interval\n"
  "                            If 2 : Also compute sample sd, cov and ChiSq Pvalue from all non-outliers molecule alignments across each consensus label interval\n"
  "  -PVcontigsplit <Pvalue> <N1> <N2> <CntFile> : split input reference map at internal site interval based on count of molecules with alignment on both sides having Pvalue no more than specified <Pvalue> :\n"
  "                                               Anywhere this count is less than <N1> AND the count rises to at least <N2> on both sides the contig will be split. Example -PVcontigsplit 1e-6 3 6 IDfile [Default : Off]\n"
  "        -ReplaceCov : causes the Pvalue based coverage count to be output in place of the actual coverage as the consensus contig coverage\n"
  "  -MinSplitLen <len> : minimum length in kb of contig split generated by -contigsplit OR of complete refined map [Default 100]\n"
  "                         The new split contig output file name will be PREFIX_contig<N>_refined.cmap, where <N> is the integer value in the file <CntFile>, which is incremented by 1\n"
  "                         If no contig split has the minimum length an empty file PREFIX_contig<ID>_refined.cmap will be generated\n"
  "  -TrimFactor <TF> : Set the initial value of the outlier fraction used to estimate consensus interval sizes (not used with -MultiMode) [Default 0.26]\n"
  "  -MaxCov <C> : limit coverage at each consensus location to <C> by discarding maps with the lowest alignment Pvalue [Default OFF]\n"
  "        -MaxCovWt <0/1> : Use weighted coverage (see -BestRefWT) for -MaxCov. Slower but more accurate. [Default OFF or 0]\n"
  "  -FilterXmap <0/1> : If != 0 : Filter out molecule alignments in unrefined .xmap to match what was filtered out from refined .xmap due to -MaxCov, -TE, -minoverlap. -RepeatMask or -extonly\n"
  "                      Use 0 to see all molecule alignments in unrefined .xmap  [Default ON or 1]\n"
  "  -FilterXmapWT <0/1> : If != 0 : Filter molecule alignments in unrefined .xmap to match what was filtered out from refined .xmap due to -refineMinWT [Default OFF or 0]\n"
  "  -FilterXmapSplit <0/1> : If != 0 : Filter molecule alignments in unrefined .xmap to match what was filtered out from refined .xmap due to -keepsplitmaps 0 [Default OFF or 0]\n"
  "       NOTE1 : Molecule alignments with low weight may additionally be filtered out just before output from both refined and/or unrefined .xmap based on -MinMapwt.\n"
  "       NOTE2 : If -keepsplitmaps 1 is used : Molecule alignments moved to split contigs are retained in both unrefined and refined .xmap of the original contig.\n"
  "  -Site_Pen <score> [ <extscore> ] : Score penalty for each consensus site to balance FP and FN rate on consensus maps [Default 3.0 0.0]\n"
  "                                     If specified, <extscore> is a lower penalty that is applied during initial extension refinement, before reverting to <score> for final refinement [Default 3.0 0.0]\n"
  "  -skip_dist <D> [<Dext>] : Speed up refinement by not checking for new sites more often than every <D> (kb) [Default OFF or 0]\n"
  "                              If <Dext> is specified and larger than <D>, this value is used during initial extension refinement, before reverting to <D> for final refinement [Default 0.500]\n"
  "  -skip_score <S> [<Sext>] : Speed up adding/deleting sites during refinement by NOT rechecking sites that cause Likelihood Score to drop by <S> or more [Default 20.0]\n"
  "                              If <Sext> is specified and smaller than <S>, this value is used during initial extension refinement, before reverting to <S> for final refinement [Default 10.0]\n"
  "  -FN_Refine <FNmul> [<FNmulExt> <FNmin> <FNminExt> ] : Increase FN value during refinement by factor <FNmul>. If <FNmulext> is specified and larger than <FNmul>, this value is used \n"
  "                              during initial extension refinement, before reverting to <FNmul> [Default OFF]\n"
  "                              If <FNmin> is specified AND FNmul > 1, FN is increased to max(FNmin, FN*FNmul). Similarly <FNminExt> is the minimum absolute value of FN used during initial extension\n"
  "                              refinement if FNmulExt > 1. This can model variation in nicking efficiency to avoid having more FN than FP in consensus Map\n"
  "  -MinSF_Refine <Minsf> <MinsfExt> : Increase sf during refinement to at least <Minsf> and to at least <MinsfExt> during initial extension refinement [Default OFF or 0]\n"
  "  -MinSR_Refine <Minsr> <MinsrExt> : Increase sr during refinement to at least <Minsr> and to at least <MinsrExt> during initial extension refinement [Default OFF or 0]\n"
  "  -MprobevalTest : More thorough refinement (Typically 4x slower) [Default OFF]\n"
  "  -maxContigSiteDensity <D> : If D != 0 : Use -refine 1 -EndTrim 0 for refinment of contigs with Site Density greater than D per 100kb [Default OFF or 0]\n"
  "                              NOTE: Contig will be output unchanged except if -rres was used\n"
  "  -MinKB <L> : Try to keep labels at least <L> kb apart during refinement. Values between 0 and rres*0.5 can speed up refinement. <L> = 0 implies a value of rres*0.5 [Default 0.001 kb]\n"
  "-Haplotype <HapSitePvalue> <HapIndelPvalue> [ <PhasePvalue> ] : During refinement detect Heterozygous sites and indels (with Pvalue = <HapSitePvalue> & <HapIndelPvalue> respectively) and then try to phase them\n"
  "      (with Pvalue = <PhasePvalue> Default = 0.01)\n"
  "      Will output a haplotype format .hmap with both haplotypes and two .cmaps with one Allele each. -contigsplit is currently only effective during output of the seperate Alleles.\n"
  "      If the original contig id was N, the seperate Alleles will use 10*N+1 and 10*N+2 as id and the .hmap with both haplotypes will use 10*N as id.\n"
  "   -HapSiteRes <minKBadd> : Do not add Haplotype Site within <minKBadd> of another Site on same Allele\n"  
  "   -HapSiteResDelay <minKBadd2> <N> <L> <A> : Do not add(or create) Haplotype Site within <minKBadd2> of another Site on the same Allele during the first N iterations. If L = 2, also block Haplotype Site\n"
  "                                          within <minKBadd2> of another Site on the opposite Allele. Has little affect with L = 1 unless minKBadd2 > minKBadd [Default : OFF or 0 11 1 0]\n"
  "                                          If A = 1 : Only applies to adding a new Haplotype Site (and NOT also to change a non-Haplotype Site into a Haplotype Site) \n"
  "   -HapSiteUnphased <SNR> <Pvalue> : Filter out unphased Het site unless its mean SNR (if present) exceeds <SNR> AND Confidence is better than -Log10(<Pvalue>) [Defaults : OFF]\n"
  "   -HapMinCov <C> <R> <P> [ <S> <Pvalue> <RL> <PvalueL> <MC> <MM> <NN> ] : Filter out Het site or indel if the local total coverage of all molecules is below C OR\n"
  "                                  If the fraction of minor Allele coverage is below R (NOTE : The minor Allele refers to the Allele with lower coverage, so this fraction is 0.5 or less)\n"
  "                                  Molecules are assigned to one allele if that allele's relative probability exceeds P. P should be between 0.5 and 0.9 [Defaults : OFF]\n"
  "                                  If S is specified, the C & R cov filter is NOT applied to indels larger than S kb\n"
  "                                  If Pvalue is specified, filter is only applied to Het sites with Confidence worse than Pvalue\n"
  "                                  If RL is specified AND Indel is under 35kb AND the minor Allele is smaller in size THEN R is increased by RL * min(1, label interval in region of indel / 100kb)\n"
  "                                        Also log(Pvalue) is increased by the same proportion for Indels. This is intended to guard against DNA knots being misclassified as deletions\n"
  "                                  If PvalueL is specified, filter is only applied to Het Indels with Confidence worse than -log10(Pvalue) * (RL/R) (PvalueL defaults to same value as Pvalue)\n"
  "                                  If MC is specified, always filter out any Het Indel with minor Allele coverage below absolute value MC, except if it includes too many SNPs (see MM,NN) (Default 0.0)\n"
  "                                  If NN & MM are specified, never filter out any Het Indel that overlaps more than NN SNPs & MM SNPs / 100kb : this overrides even -HapMinCovPL [Default OFF or 99.0 999])\n"
  "                            NOTE : All coverages are sums of molecule alignment weights. If the LR value for Alleles A and B are LRa and LRb, the molecule weight = (LRa+LRb)/(LRa+LRb+2*LRbias)\n"
  "                                   Molecule weights are summed into one of 3 coverage values based on the value of LRa/(LRa+LRb) falling below 1-P (covB = allele B coverage), \n"
  "                                   above P (covC = allele A coverage) and neither (covC = unclassified coverage). Total coverge is covA+covB+covC. Minor Allele fraction is min(covA,covB)/(covA+covB)\n"
  "                                   Molecules that are junk or match neither of 2 alleles well, will have low weights and contribute little to any of the 3 coverage values\n"
  "                                    (see also -TB which can downweight molecule alignments further based on the alignment Pvalue)\n"
  "        -HapMinCovPL <P1> <P2> [<L> <CL>] : Override -HapMinCov action for Hap Indels in certain cases based on ChiSquare Pvalue for sizes of minor allele:\n"
  "                                       NEVER filter out Hap Indels IF ChiSquare Pvalue for map sizes of minor allele is above P1, \n"
  "                                               UNLESS (total coverage is below HapMinCov's C value OR minor Allele coverage is below CL value)\n"
  "                                       ALWAYS filter out HapIndels IF ChiSquare Pvalue for map sizes of minor Allele is below P2 [Default : OFF or 1.0 0.0 0 0]\n"
  "                                               However Het Indels with more than NN SNPs and MM SNPs/100kb are never filtered out (see -HapMinCov)\n"
  "                                               Also (see -MinorIsBigger) Hap Indels are not filtered if minor Allele is larger in size.\n"
  "                                  The ChiSquare Pvalue is computed 2 ways based on estimating sizing errors from either autoNoise estimate error parameters or local sizing error of major allele : \n"
  "                                      If L = 0 : Both versions of Pvalue must agree (both must exceed P1 or be below P2) before -HapMinCov action is overridden [Default L = 0]\n"
  "                                      If L = 1 : Only the Pvalue based on autoNoise estimate is used.\n"
  "                                      If L = 2 : Only the Pvalue based on major allele sizing error.\n"
  "                                  If CL > 0 is specified : the minor Allele coverage must be at least CL, otherwise the ChiSquare Pvalue is considered unreliable and no override of the\n"
  "                                      regular -HapMinCov action occurs [Default CL = 0]\n"
  "                                  NOTE : This is a more sensitive method to handle DNA knots and may be able to handle low frequency minor allele (eg with L = 1 and P1 = P2 = 0.01)\n"
  "        -MinorIsBigger <0/1> : If 1, Disables -HapMinCov action for Hap Indels if minor Allele is larger in size than major Allele\n"
  "                                               UNLESS (total coverage is below HapMinCov's C value OR minor Allele coverage is below HapMinCovPL's CL value)\n"
  "                                     Useful if DNA knots only cause shrinkage in size [Default: ON or 1]\n"
  "        -HapEndTrim <T> [ <N> <F> <R> <L>]  Also trim ends of Allele A as long as (covA + N) * (covA+covB+covC + 2N) / (covA+covB + 2N) - N is below T [Default : 0.0 1 0.5 0.2 99.9 [OFF]]\n"
  "                                    This method works better than -EndTrim <T>, which would approximately check if covA + covC/2 is below T and may fail to trim some ends with high covC values:\n"
  "                                                                This method distributes covC to covA & covB in a pro-rata manner\n"
  "                                    A small positive value of N is used to increase covA & CovB : this is needed to handle ends where covA = covB = 0 and covC is non-zero\n"
  "                                    Undo trimming if either Allele loses more than a fraction <F> of its original labels\n"
  "                                    If R > 0 : Increase T to at least R times the highest coverage over the next L kb\n"
  "   -HapThresh <HapSite> <HapIndel> [ <HapSiteWin> <Win> <HapIndelSum> <HapIndelDensity> <skipEndSNP> <MinHapCov> <MinHapRelCov>] : Output seperate contig .cmap for each allele + .hmap ONLY if at least <HapSite> Het Sites OR at least <HapIndel> Het Indels were found [Default: 1 1]\n"
  "                                     If <HapSiteWin> and <Win> are specified, also seperate alleles if at least <HapSiteWin> Het sites were found somewhere within a window of <Win> kb [Default 0 100.0]\n"
  "                                     If <HapIndelSum> is specified, sum of absolute Het Indel sizes must add up to at least <HapIndelSum>, otherwise treat same as no Het Indels found [Default 0]\n"
  "                                     If <HapIndelDensity> is specified, sum of absolute Het Indel sizes must add up to <HapIndelDensity> kb per Mb [Default 0]\n"
  "                                     If <skipEndSNP> is != 0, then don't count SNPs beyond the outermost non-SNP labels [Default 1]\n"
  "                                     If <MinHapCov> != 0, each Allele's average coverage based on maps assigned to Allele with 70% certainty must be at least MinHapcov [Default 0]\n"
  "                                     If <MinHapRelCov> != 0, each Allele's average coverage based on maps assigned to Allele with 70% certainty must be at least MinHapRelCov of total [Default 0]\n"
  "                                     NOTE : This is applied after -HapSiteUnphased and -HapMinCov have filtered out low quality Hap calls\n"
  "   -HapIndelMerge                  : Treat consecutive Hap Indels seperated Hap Sites as a single indel. Prevents large Haplotype Indels with labels on one Allele from being penaltized\n"
  "                                      multiple times by <HapIndelPvalue>. Forces -RefineFix option [Default: OFF]\n"
  "   -HapIndelSkip <N1> <N2> [ <E> <EM> ] : Skip expensive Haplotype Indel removal check N1 times after every removal check that failed to improve log likelihood by E + EM * maps or more\n"
  "                                     Skip expensive Haplotype best alignment update N2 times after every update that failed to improve log likelihood by E + EM * maps [Default: 0 0 1.0 0.0]\n"
  "   -HapMinLL <E1> <E2> <E3> <E4> <E5> <E6> : Terminate Haplotyping if log likelihood (LL) improvement during previous iteration is less than E1 for all interval size changes (respectively E2\n"
  "                      for Hap Indel and E3 for SNP changes). After -rres & -cres label merging, replace E1,E2,E3 by E4,E5,E6 (respectively) : typically E4 >= E1, E5 >= E2, E6 >= E3\n"
  "                      Increase these values to speed up Haplotyping at risk of incomplete convergence [Default : 1.0 1.0 0.1 2.0 1.0 1.0]\n"
  "                      NOTE: For backward compatibility, unless -HapLocalStop 1 is specified, E2 is replaced by E5 until stringency of Hap Indels is increased to final value\n"
  "   -HapFast <N> : N >= 1 : Speed up Haplotyping by recomputing likelihood only where consensus was changed (will still occassionally check full likelihood and backtrack if needed) [Default : OFF or 0]\n"
  "                  N >= 2 : Avoid global re-alignment of molecules to handle large changes on consensus\n"
  "                  N >= 3 : Avoid checking alternate ways of dividing Hap Indel when a new SNP is added inside a Hap Indel\n"
  "                  N >= 4 : Try making multiple SNP changes in a single step, backing off to verifying each change only if this no longer improves Likelihood\n"
  "                  N >= 5 : Try making multiple HapDelta changes in a single step, backing off to verify each change only if this no longer improves Likelihood\n"
  "   -HapLocalStop <N> [ <S> <S2>]: Terminate Haplotyping locally in large contigs and only work on regions that are still making progress in at least one Allele, based on -HapMinLL [Default : OFF or 0 1]\n"
  "                      Contigs are divided into regions of equal size (in kb) equal to 3x the average label interval size\n"
  "                        N >= 1 : Apply to interval size or Hap Indel changes\n"
  "                        N >= 2 : Also apply to Hap Site changes\n"
  "                        N >= 3 : Also skip reactivation of all regions when stringency of Hap Indels is increased\n"
  "                      Any type improvement in one region triggers continued checking of that region and the nearest S neighboring regions [Default: S = 1]\n"
  "                      If N==0 : Any Hap Delta improvement in one region triggers continued checking in nearst S2 neighboring regions [Default: S2 = 16]\n"
  "   -HapFilterLast <0/1> : If 0 : Apply -HapMinCov filters based on Allelic imbalance before  -rres & -cres based label merging [Default: 0]\n"
  "                          If 1 : Apply -HapMinCov filters based on Allelic imbalance last (after -rres & -cres based label merging) \n"
  "   -HapMergeSNPs <0/1> : If 1 : Try to merge nearby SNPs (with distance specified by -cres * 500bp) that are on opposite Alleles into a regular site [Default: OFF or 0]\n"
  "   -HapIndelHiRes <N> [ <M> ] : Number of iterations of High Resolution scan for Haplotype Indels. Increasing this value may find more Hap Indels at the expense of runtime. [Default 1]\n"
  "                                <M> : Initial Indel sizing ladder is a geometric progression with size range = 0.100kb ... 409.6kb, number of sizes = 1 + (24 << M) [Default M=2]\n"
  "   -RefineSwitchFade <0/1/2> : If >= 1, turn off the *Switch options during later stages of Haplotype Refinment (more stable and rapid convergence) [Default: OFF or 0]\n"
  "                               If 2, set -outlierType to value used when adding/deleting labels (instead of the -outlierTypeRef value)\n"
  "   -MaxDelta <H> : Largest Het Indel size in kb to searchfor during Haplotyping. Smaller values are faster. [Default 400.0 kb]\n"
  "   -HapMaxIter <N> <F> : Max Number of EM iterations per Haplotyping stage = N + F * TotalLabels. [Default 16 0.25]\n"
  "   -HapPhaseRev <N> <M> : Number of times all Haplotype phase reversals are checked is N (+ one final time if M == 1). [Default 2 1]\n"
  "   -HapCres <N> <S> : Number of times -rres -cres are applied. If S==1, don't resume Haplotyping after last -rres -cres is applied. [Default 2 0]\n"
  "   -HapFilterStop <S> : If S==1, Don't resume Haplotyping after -HapMinCov filters. [Default 0]\n"
  "   -FastSnpsStartup <N> : If > 0 : start Haplotyping with this iteration instead of 0 when SNPs/Indels are already present (must be <= 10) [Default OFF or 0]\n"
  "   -FastSnps <N> : Switch to Fast SNP detection after N iterations (See -HapFast with N >= 4). Should be at least 1 to avoid phasing errors [Default 1]\n"
  "   -FastIndel <N> : Switch to Fast Hap Indel detection after N iterations (See -HapFast with N >= 5). [Default 1]\n"
  "   -MaxDeltaIter <N> : Maximum number of iterations in a row with only sizing changes (no Hap Indel or SNPs) before terimination is forced [Default 5]\n"
  "   -MaxSwitchIter <N> : Maximum Hap Indel iterations before Switch options must be turned off AND transition to more stringent Hap Indel & SNP Pvalue starts [Default 100]\n"
  "   -SplitHmap <0/1> : By default (-SplitHmap 0) the Haplotype nature of the input Cmaps is ignored and the average of the 2 Alleles treated as a regular non-Haplotype Cmap [Default 0]\n"
  "                      With -SplitHmap 1 input consensus ref maps (see -ref) that are already Haplotyped are separated into 2 Allele maps (see -RefineHap for two ways to proceed). [Default 0]\n"
  "                      Also applied separately to input query maps (see -i ) if -maptype 1 is specified : has no effect on input from BNX files\n"
  "                      NOTE : haplotype information is preserved after splitting ref maps (see -RefineHmap) but discarded after splitting query maps\n"
  "   -RefineHmap <0/1> : By default (or with -RefineHap 0) input maps that are already Haplotyped are separated into 2 Allele maps and each one Refined + Haplotyped independently of the other. \n"
  "                      Specifying -RefineHap 1, will cause Haplotyped input maps to be preserved and Haplotyped starting from the previous Haplotype map. In this case non-Haplotype refinement\n"
  "                      will NOT be performed before Haplotyped refinement. Input Haplotype maps can be from a previous haplotype refinement OR from using -pairmergeHap [Default OFF or 0]\n"
  "                      Note that with -SplitHap 0 there is no effect from specifying -RefineHap.\n"
  "   -rresDefer <0/1> : Delay application of -rres until after first Haplotyping stage (this is always the case with -cres) [Default : 0 or OFF]\n"
  "   -cresCheckDel <0/1> : When merging labels with -cres, check if deleting either Label scores better than merging them at the midpoint [Default : 0 or OFF]\n"
  "   -refineDefer <0/1> : skip regular refinement and proceed directly to Haplotyping [Default : 0 or OFF]\n"
  "-EndTrim <C> [ <R> <N> <W> ] : Refined consensus map ends are trimmed back while coverage is below <C> [Default 7.9 0.0 0.0 0.0]\n"
  "                       If <R> is specified : Additionally trim unlabeled ends to R'th percentile in unlabeled lengths (in descending order : 0.20th is larger than 0.50th percentile)\n"
  "                       If <N> is specified : If at least N molecules have ends that matched the median (R=0.50) unlabeled end size (to within sizing error) mark end as non-extendable\n"
  "                                             AND increase R to 0.50 (if needed) to recompute trimmed unlabeled end size [Default OFF]\n"
  "                       If <W> is specified and non-zero : Only count molecules with map weight above W (instead of counting molecules displayed in refined XMAP, see -MinMapwt)\n"
  "    -EndTrimNoRefine <L> : Don't refine consensus ends more than <L> kb beyond initial estimate of consensus map ends [Default OFF] (NOT YET IMPLEMENTED)\n"
  "    -EndTrimCQ <0/1> : Trim ends to within -CovTrimLen of start of valid -ChimQuality score values [Default OFF or 0]\n"
  "    -EndTrimRatio <F> : Before applying -EndTrim <C>, trim back the last K labels at either end if highest label count for all of these K labels is <= F times the label cnt of the next label.\n"
  "                       Avoid F values larger then -EndTrim R. [Default: OFF or 0]\n"
  "    -EndTrimInternal <C2> <L2> : Trim back up to any internal region with coverage below C2, but not more than L2 kb beyond original ends [Default OFF or 0 0]\n"
  "-everb : display location of each misaligned site [Default OFF]\n"
  "-ErrDist <distance> : minimum distance in KB of FP or FN site from nearest true positive site[Default 2.0]\n"
  "-BppOutput <0/1> : If 1, output bpp value for each Molecule in .map file [Default 0]\n"
  "-mapped <fileprefix> : Output those input maps that mapped to the reference map (whose alignments are in the .xmap output) as BNX file <fileprefix>_contig<N>.bnx (Requires -ref)\n"
  "   -grouped <group_manifest> : output from -mapped argument is split into seperate files for each group of reference contigs listed in <groupfile> : \n"
  "                          The first line of <group_manifest> contains reference ids for group 1 output as <fileprefix>_group1.bnx etc\n"
  "   -centeronly <L> : output only those maps that do NOT align within <L> kb of the reference ends\n"
  "    NOTE : map output is always .bnx unless input files were .cmap AND (-maptype 0 OR -bnx) were NOT used : in that case .cmap files are output\n"
  "-unmapped <fileprefix> : Output those input maps that did NOT map to the reference map (whose alignments are NOT in the .xmap output) as BNX file <fileprefix>.bnx (Requires -ref). Also uses -mapped-unsplit\n"
  "    NOTE : map output is always .bnx unless input files were .cmap AND -maptype 0 or -bnx options were NOT used : in that case .cmap files are output\n"
  "-partmapped <fileprefix> : Output those input maps that only partly mapped to the reference map due to -endoutlier as BNX file <fileprefix>.bnx (Requires -ref). Also uses -mapped-unsplit\n"
  "    NOTE : map output is always .bnx unless input files were .cmap AND -maptype 0 or -bnx options were NOT used : in that case .cmap files are output\n"
  "-mapped-unsplit <0/1> : If 1 : Output from -mapped,-unmapped,-partmapped and .xmap output  will not split into separate files for each reference contig\n"
  "            [Default : 0 IFF number of -ref files matches the number of contigs AND -id was specified non-zero]\n"
  "-BNXnoRunData <0/1> : If 1 : Don't output RunData header lines (if present from BNX 1.x) to BNX output\n"
  "            Useful to reduce size of mapped BNX files. Downgrades BNX format from 1.x to 1.0 [Default 0]\n"
  "-Repeat <minlenKB> : Check each input map for repeats of minimum length <minlenKB> by aligning them with themselves with the specified minimum offset (in either orientation)\n"
  "-RepeatMask <MaxShift> [ <PvalueRatio> ] : Ignore alignments which can be shifted by 1 (OR more, up to MaxShift) sites without mean-square sizing error Pvalue decreasing by more than PvalueRatio [Default 0 0.5]\n"
  "                    NOTE : A smaller PvalueRatio will classify more alignments as Repeats, while PvalueRatio >= 1 is unlikely to classify any alignment as a Repeat.\n"
  "    -RepeatRec <LogPvalueRatio> <MinShiftRatio> [ <MaxShiftRatio> <MaxLabelDrop> <LogPvalueRatioPerLabel>] : Improved version that checks the alignment array to see if the best shifted alignments\n"
  "              have LogPvalue decreasing by a factor of no more than <LogPvalueRatio> [Default OFF or 1 0 3 -1 0] Recommend 0.6 0.6 1.4 1 0.9\n"
  "              The minimum shift in kb along the entire alignment must be at least <MinShiftRatio> times the rightmost alignment shift (in kb) based on the specified number of sites shifted\n"
  "              The maximum shift in kb along the entire alignment must not exceed <MaxShiftRatio> times the rightmost alignment shift (in kb) based on the specified number of sites shifted\n"
  "              If <MaxLabelDrop> is specified, the shifted alignment cannot lose more aligned labels than Shift + MaxLabelDrop\n"
  "              If <LogPvalueRatioPerLabel> is specified, the shifted alignment's LogPvalue per label interval must exceed <LogPvalueRatioPerLabel> times original value\n"
  "    NOTE: with -MultiMatch alternate alignments that fail the main score thresholds (eg -T) are not checked for repeats (to save time)\n"
  "-sv <L> : Detect structural variations (SV) based on xmap entries. Enables -pairsplit and -endoutlier using default values 1e-6 for both unless otherwise specified. '-sv' options below apply only when this option is used with L > 0. [Default OFF or 0].\n"
  "  -svMinConf <X> : include only SVs with raw confidence >= X (only if raw confidence is defined for that SV) [Default 3]\n"
  "  -svAlignConfByLen <X> : discard SVs based on alignments with confidence / length (kb) is < X [Default 0.07]\n"
  "  -svAlignConfByInterval <X> : discard SVs based on alignments with confidence / Label-Intervals is < X [Default 0.0]\n"
  "  -svConfFile <file> : use <file> to replace SV confidence scaling : converts raw confidence to PPV (Positive Predictive Value) based on a benchmark data. [Default values loaded at compile time]\n"
  "  -svConfMinsize <L> <NA> : use value <NA> for scaled confidence value when indel size is below <L> kb [Default 0.5 -1.0]\n"
  "  -svSideMinLen <kb> : For SV detection, minimum aligned length for each matchgroup [Default 50.0]\n"
  "  -svIndelMaxoverlap <kb> : Sets sv_indel_maxoverlap (Obsolete) [Default 50.0]\n"
  "  -svDelMaxQryOverlap <kb> : Set max query overlap allowed for deletions : larger overlaps are assumed to be FPs caused by repeat compressions during Assembly [Default: 110kb]\n"
  "  -svDelBigQryOverlap <kb> : If query overlap for deletion exceeds this value, reduce deletion confidence to minimum value (see -svMinConf) [Default: 60kb]\n"
  "  -svInsMaxGapOverlap <F> : Filter out Insertion, if 3rd MG overlaps gap in Qry by at least a fraction F. [Default: 0.5]\n"
  "  -svDelMaxGapOverlap <F> : Filter out Deletion, if 3rd MG overlaps gap in Ref by at least a fraction F. [Default: 0.5]\n"
  "  -svEndMinLen <kb> : For end-type SVs (unaligned ends of contigs), at least this length must be unaligned to be an SV [Default 50.0]\n"
  "  -svEndMinNLabels <N> : For end-type SVs (unaligned ends of contigs), at least this number of labels must be unaligned to be an SV [Default 5]\n"
    //"  -svIndelTinySize <kb> : Indels whose size is <= this (kb) are called deletion_tiny or insetion_tiny [default 1.0]\n" //with RefSplit we no longer want to use this, though it is still enabled
  //"  -svPalindromic <ConfidenceDelta> : Mark alignment as palindromic if best alignment in opposite orientation has confidence within <Delta> of overall best alignment. [Default 2.0]\n" //this not yet implemented
  "  -svTransMask <mask_file> : input mask file for filtering translocations: must provide chromosome (reference contig id) and coordinate for each breakpoint; translocation calls in 30 kb window are flagged \"common\" and filtered in Pipeline. Ignored unless -sv is used.\n"
  "  -svTransMaxGap <gap size in kb> : maximum separation of matchgroups on query allowed for translocations [default 500 kb]\n"
  "  -svTransMaxOverlapFrac <F> : maximum overlap fraction of matchgroups on query allowed for translocations [default 0.3]\n"
  "  -svTransMaxOverlapSize <L> : maximum overlap size (kb) of matchgroups on query allowed for translocations [Default 140 kb]\n"
  "  -svInvMaxOverlapSize <L> : maximum overlap size (kb) of matchgroups on query allowed for inversions. If undefine uses -OverlapFilter <QryRatio> instead. [Default: 140 kb]\n"
  "  -svInvConfMinInterval <L> <I> : When computing confidence of small inversions, count only aligned intervals above L kb.\n"
  "                                  Filter out small inversions with fewer than I aligned intervals above L kb. [Default: 1 2]\n"
  "  -svZygosity <0/1> : If 1, compute and output SV zygosity in .smap [Default ON or 1]\n"
  "  -svSize <0/1> : If 1, compute and output SV size in .smap [Default ON or 1]\n"
  "  -svFreq <0/1/2> : If >= 1, compute and output SV frequency fraction in .smap [Default OFF or 0]\n"
  "                    If >= 2, also output the total coverage and SV coverage used to compute the SV frequency.\n"
  "  -svTransOrient <0/1> : If 1, compute and output Translocation Orientation in .smap (+/+, +/-, -/+ or -/-) [Default OFF or 0]\n"
  "  -svDupSplit <maxqryoverlap> <maxqrygap> <maxend> <L> <N> [<T>]: Specify in labels the maximum query matchgroup overlap, gap and unaligned ends allowed for duplication_split calls : \n"
  "                 all except gap should be small. Any qry end below <maxend> unaligned labels can be extrapolated provided the nearest matchgroup is at least <L> kb AND <N> labels in length.\n"
  "                 Either both qry ends can be extrapolated or one of the qry ends can be extrapolated AND either T >= 3 OR one of the following conditions must hold:\n"
  "                         A. T >= 1 AND MGs overlap on reference\n"
  "                         B. T >= 2 AND MG that cannot be extrapolated must be at least 50% of the reference interval spanning both MGs\n"
  "                 NOTE : <L> and <N> should not exceed corresponding -CloneLimitLen & -pairmergeHmap args. [Default: 2 30 3 200 15 3]\n"
  "  -svDup_minRefOverlapSites <N> : Minimum overlap in sites between matchgroups in reference to call an Duplication [Default 4]\n"
  "  -svDup_maxQryGapMult <M> : Maximum Qry Gap size (in sites) as a multiple of the reference matchgroup overlap (in sites) to call a Duplication [Default 0.8]\n"
  "  -svDup_maxSpacing <L> : Maximum Gap (in kb) on Query between Duplication events (should correspond to max distance that can be phased) [Default OFF or 0]\n"
  "  -svInvDup_minRefOverlapFrac <F> : Minimum overlap fraction between matchgroups in sites in reference to call an Inverted Duplication [Default 0.5]\n"
  "                                    Also used to suppress inversions if overlap fraction in Sites is >= F for either reference and query: see -svInv_maxRefOverlapFrac or -svInv_maxQryOverlapFrac\n"
  "  -svInv_maxRefOverlapFrac <F> : Suppress inversions if overlap fraction between matchgroups in sites in reference is above F. Defaults to -svInvDup_minRefOverlapFrac\n"
  "  -svInv_maxQryOverlapFrac <F> : Suppress inversions if overlap fraction between matchgroups in sites in query is above F. Defaults to -svInvDup_minRefOverlapFrac\n"
  "  -svInvDup_maxQryGapMult <M> : Maximum Qry Gap size between matchgroups (in sites) as a Multiple of the reference matchgroup overlap (in sites) to call an Inverted Duplication [Default 0.8]\n"
  "  -svInvDup_maxSpacing <L> : Maximum Spacing(in kb) on Query between Inverted Duplication events (should correspond to max distance that can be phased reliably) [Default OFF or 0]\n"
  "  -svInv_TrimFrac <Min> <Max> : For Inversion calls with non-inverted matchgroup fully overlapped by inverted matchgroup on ref : the upper range of acceptable fraction of the inverted matchgroup\n"
  "                                that needs to be trimmed (or ignored) in the overlap region and beyond. If less than <Min> needs to be trimmed, rule out Inverted Duplication. If more than <Max> \n"
  "                                needs to be trimmed, rule out Inversion. In between those values, choose Inverted Duplication unless ruled out by -svIndDup_maxSpacing. [Default OFF or 0.0 to 1.0]\n"
  "  -svInv_minRefSites <N> : Minimum non-overlapped reference sites in inverted matchgroup of Inversion [Default 4]\n"
  "  -svInv_minQrySites <N> : Minimum non-overlapped reference sites in inverted matchgroup of Inversion [Default 0]\n"
  "  -svIndelConfBySize <0/1> : If 1, base Indel confidence on size difference only (ignoring misaligned labels) [Default OFF or 0]\n"
  "  -svSmallDelConfirm <L> : For Deletions in ref intervals below L kb, confirm deletion by including neighboring intervals [Default OFF or 0]\n"
  "      -svIndelRangeExpand <R> <S> <F> : For all Indels recheck sizing error after expanding alignment range if boundary label has misresolved labels extending over R kb or more,\n"
  "                                    OR the boundary labels on ref or qry are not seperated from next label in either direction by at least S kb. However if the range expansion does NOT change\n"
  "                                    the Indel size by at least the fraction F, then undo the expansion of the alignment range, but use the modified SVsize to adjust coordinates [Default 0.8 1.7 0.2]\n"
  "  -svGapOverlapRatio <R> : If SV gap region is overlapped in both ref and qry by another matchgroup by at least R times the smaller of the matchgroup's or gap's size, suppress the SV [Default: OFF or 0]\n"
  "  -svPrimaryOrientation <R> : Check if MG has a primary orientation based on total Qry sizes for all MGs with the same qryid and refid in each orientation : If total size for one orientation \n"
  "                          exceeds the total size for the other orientation by <R>, disallow Indels in orientation other than primary orientation or inversions with invalid orientation. [Default: OFF or 0]\n"
  "  -svRepeat <S> <N> <P> : Mark translocations & inversions as repeat IF there are at least N consecutive intervals that differ by at most S times SD error on either query or ref of either matchgroup\n"
  "                          If P > 0 : Also require pro-rata confidence of the non-repeat intervals (minus P times the number of repeat regions) of the matchgroup to drop below -T [Default: 1.7 2 1.5]\n"
  "  -svTransMaxOverlap <R> : Disallow translocation if 3rd MG overlaps one of the 2 trans MGs by at least a fraction R of the trans MG qry size and maps back to the other reference [Default: 0.7]\n"
  "                           If the 3rd MG confidence exceeds the minimum confidnece for translocation MGs, include overlap with the gap in the total overlap fraction.\n"
  "  -svTransMaxGapOverlap <G> : Disallow translocation if a 3rd MG overlaps the query gap of the 2 translocation MGs by at least a fraction G of the 3rd MG qry size [Default: 0.7]\n"
  "                              NOTE : 3rd MG confidence must exceed minimum confidence for translocation MGs.\n"
  "  -svMaxGap <G> [<GInv>]: Reclassify any 2-MG SV where the reference and query gap size differ by over G kb as an intra-chromosomal translocation\n"
  "                          If GInv is specified, treat 2-MG inversions with reference or query gap size over GInv as an intra-chromosomal translocation [Default : 5000 5000]\n"
  "-bed <filename.bed> : input bed file for filtering SVs (currently insertions and deletions only) for gaps (N base). Only type (4th column) gap are used.\n"
    //"-scaffold <scaffold file> : Must supply a .scaffold file. After alignment, checks for gaps spanned by query or reference maps and outputs .amap describing these. Must use -ref.\n" //this feature is obsolete
  "-maxthreads <N> : Limit number of threads to N (or OMP_NUM_THREADS if less). 0 means use OMP_NUM_THREADS. Default with -ref is OMP_NUM_THREADS, otherwise (pairwise) 1\n"
  "-RefineThreads <N> : Limit number of threads to N. Default is same as -maxthreads\n"
  "-TotalThreads [ <M> <R>] : pro-rate -maxmem and -hashmaxmem by N/M : Assumes -maxmem and -hashmaxmem specify the total usable memory on the machine. <M> defaults to total system threads. [Default OFF]\n"
  "                 NOTE : If N > OMP_NUM_THREADS, memory is pro-rated based on original N value, before N is reduced to OMP_NUM_THREADS\n"
  "                 If R is specified, adjust M to (M*OMP_NUM_THREADS)/R : effectively M/R specifies the amount of thread overloading, the exact values are immaterial.\n"
  "-boostThreads <T> [ <K> ] : After elapsed time of <T> seconds double RefineThreads (or use <K> threads, if specified) subject to original -maxmem and -hashmaxmem. [Default OFF]\n"
  "-BestRef <0/1> : 0 means allow Molecules to align to multiple -ref contigs, 1 means allow Molecules to align to best -ref contig only[Default 0]\n"
  "    -BestRefPV <0/1> : Use Pvalue (PV) instead of likelihood to compute weights of molecule alignments\n"
  "    -BestRefMargin <val> : Prefer earlier ref contigs over later ref contigs by this margin in score (or Log10(PV) [Default 0]\n"
  "    -BestRefOutlier <kb> : Only suppress alignments with internal outlier sizing error above <kb> [Default OFF]. NOTE: Ignored if -BestRefExt is specified\n"
  "    -BestRefExt <L> <PVdelta> : modify -BestRef to suppress only extension alignments IFF there is a center alignment with another -ref contig NOT within <L> kb of ends \n"
  "                                AND a -Log10(Pvalue) that exceeds that of the extension alignment by <Pvdelta> (can be -ve) [Default OFF]\n"
  "                                (NOTE: see -svcheck, which modifies how -BestRefExt works)\n"
  "-BestRefWT <0/1> : Alternative to -BestRef that only downweights molecule alignments based on relative likelihood ( = exp(score)) so total weights of a molecule add up to 1.0 [Default OFF or 0]\n"
  "                   The reduced weights are used during Refinement, including the computation of Quality Scores and coverage\n"
  "    -BestRefPV <0/1> [<C>]: modify -BestRef to pick best -ref contig based on Log10(PV) (confidence) rather than score [Default 1 9999]\n"
  "                            If C is specified, limit alignment confidence to no more than C. If C matches Log10 confidence of -T option, all alignments have equal weight\n" 
  "    -AddFragsWT <0/1>  : Add the alignment weights of all alignments of the same query map with the same ref map in same orientation, treating them as fragments of one large alignment\n"
  "                         This can be useful with -MultiMatches if molecules are aligning across large indels of single ref map in multiple matchgroups, eg with -extsplit\n"
  "                         Risks generating FP alternate Allele in regions with tandem repeats with alternate Allele having repeat compression relative to the main Allele\n"
  "-XMAPmapWT <0/1> : Output molecule alignment weights to .xmap as extra column at right end : Useful with -BestRefWT -LRbias -TB or -Haplotype. [Default OFF or 0]\n"
  "                   See also -MinMapwt to filter out alignments with weight below specified 2nd value\n"
  "-bestalignments <N> : Output only the best <N> alignments based on Pvalue (confidence) : requires -ref [Default OFF]\n"
  "-firstalignments <N> : Compute only the first <N> (or more) alignments above threshold : may output some additional alignments if multithreaded [Default OFF]\n" 
  "-PVendoutlier [<0/1/2>] : Apply a correction to Pvalue (-T) to reflect use of End outliers [Default OFF]\n"
  "                          -PVendoutlier 2 applies a more stringent correction when there are 2 End outliers (-PVendoutlier without a value defaults to -PVendoutlier 2)\n"
  "-PVres [ <N> ]: Apply a correction to Pvalue (-T) to reflect limited resolution (-res & -bpp). <N> = 2 uses a more aggressive correction. The effect is more significant with high site density [Default OFF]\n"
  "-AlignRes [ <range> ] : Apply more complete scoring of alignments with limited resolution (-res). Results in more alignments if sites are very dense [Default 1.0\n"
  "                        <range> limits range of sites that failed to resolve to (res + resSD*3) * <range>. A Value of 1.0 (or less) disables -AlignRes [Default 2.0]\n"
  "-EndScoreThresh <S1> <S2> : For pairwise alignment use a heuristic score thresholds of <S1> and <S2> as upper bounds of the most misaligned sites in the right ends [Default 10.0 2.0]\n"
  "-fastmode <N> : positive values turn on various heuristic optimizations [Default 0]\n"
  "-fastmodetolerance <S> : tolerance value for pair oracle [Default 4]\n"
  "-ChimQuality [ <N> [ <W> ] ]: output Chimeric Quality and other assembly quality scores for each site in consensus map [Default: OFF or 0 0]. The value of <N> defaults to 1.\n"
  "                      N >= 1 : Output just the Chimeric Quality score\n"
  "                      N >= 2 : Also output two Segdup scores : triggers near Left and Right end of Segdup region\n"
  "                               Also output two Fragile Site (or Chr End) scores : triggers near Left and Right end of chromosome or Fragile Site interval\n"
  "                      W >= 1 : Give full weight to maps that match both Alleles equally : other cases are handled by using twice the Allele probability, but not more than 1.0 as Allele weight (still uses BestRefWT & LRbias)\n"
  "                      W >= 2 : Also ignore -BestRefWT && -TB && -refineMinWT, only reducing weight based on LRbias (NOTE: Maps below -refineMinWT are still skipped during refinement and only used for Quality Scores)\n"
  "                      W >= 3 : Also ignore LRbias based weight (NOTE : W >= 2 also modifies Quality score computation without Haplotyping, while W==1 only has an effect with Haplotyping)\n"
  "-contigQual : output contig quality score (average pvalue of aligned molecules) for each site in consensus map [Default: OFF]\n"
  "-xmaplen : add query and reference map length columns to .xmap output [Default: off]\n"
  "-xmapchim <N> [ <L> ] : check for chimeric query maps that align with two or more reference maps having at least <N> unique sites in each reference map [Default: OFF or 0]\n"
  "   If <L> is specified, also check for two alignments to same reference map separated by at least <L> kb. Ignored unless -nosplit 1 -MultiMatches was used. [Default: Off or 0]\n"
  "-xmapUnique <N> : suppress xmap alignments for query maps, if there is a higher logPV alignment of the same query map that overlaps all but <N> sites of the smaller alignment [Default : OFF]\n"
  "                  Will also suppress -indel calls from overlapped alignment Query regions of the lower Pvalue alignment, even if more than <N> unique query labels are present in both alignments\n"
  "                  Ignores 2nd best alignments, if -SecondBest was specified\n"
  "                  No alignments are suppressed if -MultiMatch was specified (but still suppresses -indel calls from overlapped alignment Query regions)\n"
  "-nostat [<0/1>] : skip computation of alignment statistics [Default: OFF(0) unless -refine was specified]\n"
  "-CmapSNR : output complete SNR distribution for each site of consensus CMAP output (for condensed _r.cmap and refined .cmap or .hmap)\n"
  "-SNRpvalue <filename> : Use SNR to improve Pvalue of map alignments given the background SNR distribution in <filename> : Requires complete SNR distribution in input maps\n"
  "-SNRbackground <filename> <N> : Output a sub-sampled SNR distribution of all -i input .bnx data to <filename> with <N> sample values and exit\n"
  "-CovNorm <MapLambda> : Output Normalized coverage distribution (for condensed _r.cmap and refined .cmap) assuming map lengths are exponentially distributed with specified mean kB length\n"
  "                       NOTE : Average Map length = MapLambda + MinLen(cutoff) can be estimated from the complete BNX map set using -minlen -merge\n"
  "-XmapStatWrite <filename> : Compute statistics of -i input maps (after filtering) and save them in specified file\n"
  "-XmapStatRead <filename> : Read map statistics from specified filename and use them in place of statistics of actual -i input maps\n"
  "-readLoc <0/1>: Read simulated data location information from Version 0.1 BNX input file [Default : 1]\n"
  "-indel             : output .indel file (same format as .smap) of all outliers in the .xmap file (provided -ref was used) [Default : off]\n"
  "   -indelends      : Also output end outliers as end type in .indel file [Default : off]\n"
  "   -indelNoOverlap [ <N> [<Pv>]] : If the same query region aligns to two different reference regions, do not output any indels from the overlapped region of the query [Default: OFF or 0]\n"
  "                     NOTE: This is stricter than -xmapUnique which would only suppress the indels of the overlapped region from the lower scoring alignment\n"
  "                     The rationale is that overlapped matchgroups typically involve segdups that were incorrectly merged and hence may have resulted in a hybrid map that doesn't match either segdup\n"
  "                     Does not apply unless both matchgroups have confidence above the main -T threshold, hence does not apply to -CutFlip or -MultiMatch based secondary (lower confidence) matchgroups\n"
  "                         If <Pv> is specified, also does not apply unless both matchgroups have confidence above <Pv> [NOT YET IMPLEMENTED]\n"
  "                     N <= 1: Also applied when small matchgroup is inside a larger matchroup (eg inversion) which is not desirable since no segdup is involved\n"
  "                     N >= 2: Not applied  when small matchgroup is inside a larger matchroup (eg inversion), but if -xmapUnique <M> is specified, apply if M unique (non-overlapped) labels are present\n"
  "                     N >= 3: Not applied when both matchgroups have opposite orientation AND involve the same reference Map AND are overlapped on the reference Map by at least 1 label\n"
  "   -indelMinConf <C> : Filter out indels if neither Left nor Right flanking matchgroup confidence exceeds <C> : both matchgroups have original confidence above -T before -RefSplit [Default: OFF or 0]\n"
  "   -indelChiSq <PV> <F> [ <L> <Cmax> <S> <R> <I> <SR> <Cmin> <CR> <W> <Fr> <Sd> <Ch> <CRt>]: Ignore Deletions with size below L kb for ref map interval size above I kb IF any of the following applies:\n"
  "                     IF OutlierFrac value is above F (%d) \n"
  "                     OR MolCov is below Cmax AND ChiSq score below PV AND MolSd above ExpSd + S kb + R times interval size AND MolSd above ExpSd * SR\n"
  "                     OR MolCov is below both Cmax and Cmin + CR * RangeRatio (ratio of ref map range to average label interval size) \n"
  "                            If CRt > 0 : If qry map and ref map have different number of unaligned labels, limit RangeRatio to no more than CRt\n"
  "                     OR FragileLeft+FragileRight is above Fr%\n"
  "                     OR SegDupLeft+SegDupRight is above Sd%\n"
  "                     If W > 0 : Also suppress Deletions with size below L kb in W neighboring intervals [NOT YET IMPLEMENTED]\n"
  "                     If Ch > 0 : Replace MolCov with min(MolCov, ChimNorm * ChimQuality * Ch / 100)\n"
  "                     Requires that CMAP was refined with -OutlierQuality 2 to generate OutlierFrac and MolSd, ExpSd, MolCov and ChiSq values for each interval. \n"
  "                     NOTE: Requires MolSd etc in CMAP input (see -OutlierQuality): Can be used to ignore outliers caused by DNA knots in the molecules. [Default OFF or 0 100]\n"
  "-xmapNoQueryOverlap <L> : Suppress Xmap alignments where the Query region overlaps by <L> kB with a better Pvalue alignment (processed in descending order of log Pvalue) [Default: OFF]\n"
  "-xmapFilter <Pv1> <Pv2> <PVdelta> : Changes alignment Pvalue threshold to -T <Pv2/PVdelta> , then filters (removes) Xmap alignments in the following order:\n"
  "                                   1. Apply -xmapNoQueryOverlap (if present) and remember for each ref contig the highest confidence alignment removed.\n"
  "                                   2. If a Ref contig has multiple remaining alignments, keep them IFF the best PValue is better than that of the 2nd best alignment by at least\n"
  "                                                           a factor of PVdelta OR the best PValue is better than Pv1.\n"
  "                                   3. Discard all alignments other than the best one for each Ref contig.\n"
  "                                   4. If a remaining alignment's PValue is worse than an alignment of the same Ref contig removed in step 1, only keep it if PValue is better than Pv1\n"
  "                                   5. Remove any remaining alignment if its PValue is worse than Pv2\n"
  "               NOTE1 : 0 < Pv1 < Pv2 < PVdelta < 1\n"
  ".              NOTE2 : -xmapFilter is designed to be used with -xmapNoQueryOverlap for alignment of NGS scaffolds as Reference with BNG maps are Query\n"
  "-xmapUnique <N> : will suppress output of Xmap alignments (and corresponding indels) if the Query map region is overlapped by another alignment with better Pvalue\n"
  "                  However, if the overlap is partial with <N> or more non-overlapped labels, the alignment is not suppressed, but the indels in the overlap region are suppressed\n"
  "-svcheck <SVfilename> [<extension> <mergegap> [ <confidence>]] : check indel SVs in input .indel or .smap by aligning -i input maps against two version of each SV neighborhood of <extension> kb [Default Len = 300 (kb)]\n"
  "                                  All Indels with confidence below -svMinConf are ignored [Default 3.0]\n"
  "                                  If SVs are near each other (within <mergegap> kb), all possible combinations/merges of such SVs are also considered as seperate SVs\n"
  "                                  The consensus maps are computed from the query and reference files named in <SVfilename>\n"
  "                                  An SV is corrected if the new confidence values is below <confidence>. Use +ve values to bias in favor of reference (SV not present) [Default 0.0]\n"
  "                                  Outputs corrected versions of each query map as <prefix>_contig<NNN>.cmap along with updated PValues for SVs in <prefix>.indel (otherwise same as <SVfilename>)\n"
  "                                  If -BestRef and -BestRefExt <L> is specified, only suppress molecule alignment if better alignment does NOT overlap the same reference region by at least <L> kb\n"
  "                                  NOTE : Does not handle Haplotype Indels correctly and will give confidence around 0 for them. Use <confidence> of -20 or less to avoid removing them\n"
  "-floatcov [<0/1>] : output Coverage and Occurance columns in .cmap as float rather than int [Default: 1 (float)]\n"
  "-svcompare <file1.indel> <file2.indel> <overlapKb> <Fraction> [ <Conf> <DeltaMatch> <SiteMatch> <xmapOverlap> ] : compare 2 sets of indel SVs (must use same reference maps) and display SVs that fail to overlap by at least the specified kb or Fraction\n"
  "                                  This is a standalone command that requires no other parameters, except optionally, -o -stdout -stderr -T : -T can be used to filter out SVs below specified Pvalue\n"
  "                                  <overlapKb> and <Fraction> can be -ve to allow nearby SVs to be matched, provided no other SV is closer\n"
  "                                  <Conf> (if specified) will also display the total and unique SVs above specified confidence : Note that an SV is unique only if it does match any SV (including SVs below <Conf>)\n"
  "                                  <DeltaMatch> (if != 0) requires matching SVs to be the same type (deletion or insertion) with size agreeing better than a fraction of <DeltaMatch> times larger SV [Default : OFF or 0]\n"
  "                                  <SiteMatch> (if != 0) requires matching SVs to have number of inserted & deleted sites to differ by less than <SiteMatch> sites [Default : 0]\n"
  "                                  <XmapOverlap> (if != 0) : Categorize unique SVs as confirmed FP if the other Assembly has a contig that overlaps SV site by at least this amount on both sides\n"
  "     -svfile <filename>           Append summary statistics to specified file\n"
  "               NOTE: -svcompare must be the last option. This means -svfile and all other options must be specified before -svcompare\n"
  "-align_format <0/1> : 0 means use original text format, 1 means use binary format for output of .align files [Default : 0]\n"
  "-XmapInterleaved <0/1> : 1 means use interleaved label ids in XMAP for 2-color alignments (as in _q.cmap), 0 means label ids in XMAP are numbered separately for each color (as in _q.bnx) [Default : 1]\n"
  "-MinFP <V> : Force FP values to be at least <V> at all times [Default 0.01]\n"
  "-MaxFP <V> : Force FP values to be no more than <V> at all times [Default 4.0]\n"
  "-MinFN <V> : Force FN values to be at least <V> at all times [Default 0.001]\n"
  "-MaxFN <V> : Force FN values to be no more than <V> at all times [Default 0.001]\n"
  "-MinSF <V> : Force sf values to be at least <V> at all times [Default 0.001]\n"
  "-MaxSF <V> : Force sf values to be no more than <V> at all times [Default 1.00]\n"
  "-MinSD <V> : Force sd values to be at least <V> at all times [Default -0.50]\n"
  "-MinSD_Mult <F> : Limit -ve values of MinSD so sizing variance cannot be reduced by more than fraction F vs MinSD value of 0 [Default 0.99]\n"
  "-MaxSD <V> : Force sd values to be no more than <V> at all times [Default 0.50]\n"
  "-MinSR <V> : Force sr values to be at least <V> at all times [Default 0.00]\n"
  "-MaxSR <V> : Force sr values to be no more than <V> at all times [Default 0.50]\n"
  "-MinSE <V> : Force se values to be at least <V> at all times [Default 0.00]\n"
  "-MaxSE <V> : Force se values to be no more than <V> at all times [Default 0.00]\n"
  "-MinRes <V> : Force res value to be at least <V> at all times [Default 0.01]\n"
  "-MaxRes <V> : Force res value to be no more than <V> at all times [Default 9.99]\n"
  "-MinResSD <V> : Force res valueSD to be at least <V> at all times [Default 0.01]\n"
  "-MaxResSD <V> : Force res value to be no more than <V> at all times [Default 9.99]\n"
  "-Kmax <K> : Specify <K> as the maximum number of consecutive site intervals in reference that might (mis)resolve into a single site in query [Default 6]\n"
  "-RAmem <M> <T> [ <L> ]: Speed Memory tradeoff in refalign: [Default: 1 0]\n"
  "             M=0 : Maximize speed by doing a single large memory allocation : Can be much slower if memory starts swapping.\n"
  "             M=1 : Adaptive memory allocation, but memory only grows in steps 2x : much less likely to swap than M=0 and almost as fast.\n"
  "             M=2 .. 7 : Adaptive memory allocation, memory can grow and shrink in steps of (M+1)/M : progressively less memory usage compared to M=1, but more frequent reallocation reduces speed.\n"
  "             M <= 7 AND T > 0 : memory never shrinks, but instead releases resident pages once every T seconds. All threads are used, but if resident memory exceeds -maxmem, threads will\n"
  "                      release their resident pages and pause until resident memory is below 90% of -maxmem (however the last thread is never paused, even if resident memory exceeds -maxmem)\n"
  "                      If resident memory is below a fraction <L> of -maxmem, threads skip releasing resident pages every T seconds except when threads are paused. <L> defaults to 0.60 or 60% of -maxmem.\n"
  "                           For pairalign, when resident memory is below this fraction <L>, resident memory is only checked every 2.0 seconds instead of every 20msec (in each thread)\n"
  "             M=8 : Re-allocate memory for each alignment (Minimizes memory usage)\n"
  "     NOTE : With pairalign M= 0, 8 work the same way (but M=1..7 are treated same as M=0) : also supports M=8, T > 0. Typically M = 3, T = 30 works well for both refalign and pairalign\n"
  "-RAscoremem <T> : Another speed memory tradeoff in refalign: [Default: 1]\n"
  "             T=0 : Maximize speed by doing a single score table memory allocation (May swap if refinement uses lots of memory)\n"
  "             T=1 : Minimize memory used by score tables (May be slow if refining large numbers of small contigs)\n"
  "-QPmem <M> : Speed Memory tradeoff in qprobeval (see -RAmem for explanation of values, only M=0..7 are implemented, assuming T = 0)\n"
  "-MPmem <M> <T> : Speed Memory tradeoff in mprobeval (see -RAmem for explanation of values)\n"
  "-simpleRepeat : enable simple repeat detection on all -i inputs based on comparing site intervals (not using alignments). Will output .rmap describing repeats found, and repeatStats.txt file summarizing repeat statistics.\n"
  "  -simpleRepeatTolerance <X> : fraction which adjacent intervals are allowed to differ to be considered the same repeat element [Default: 0.1]\n"
  "  -simpleRepeatMinEle <X> : minimum number of repeat elements required [Default: 5]\n"
  "  -simpleRepeatStandalone : same as -simpleRepeat, but all alignments are skipped (and their output files are not produced), and -merge is disregarded.\n"
  "  -simpleRepeatFilter <0,1,2,3> : if non-zero, output filtered maps: 1 = output only maps which do not contain any repeats\n"
  "     2 = output only maps which contain any repeats, 3 = output all maps but with all repeats masked (ie, labels removed)\n"
  "     map output format is same as input (-i) format [Default: 0]\n"
  "-CoverageIntRcmap <0/1> : If 1 : Enable Interval based Coverage in _r.cmap output file [Default : 0 or OFF]\n"
;

static char buf[LINESIZ];

/* sort long long in increasing order */
static int LLInc(long long *p1, long long *p2)
{
  return (p1[0] > p2[0]) ? 1 : (p1[0] < p2[0]) ? -1 : 0;
}

/* sort doubles in increasing order */
static int doubleInc(double *p1, double *p2)
{
  return (p1[0] > p2[0]) ? 1 : (p1[0] < p2[0]) ? -1 : 0;
}

/* sort by increasing order of map size */
static int CmapSizeInc(register Cmap **p1, register Cmap **p2)
{
  register int M1 = p1[0]->numsite[0];
  register int M2 = p2[0]->numsite[0];
  register FLOAT *X1 = p1[0]->site[0];
  register FLOAT *X2 = p2[0]->site[0];
  register FLOAT size1 = X1[M1+1];
  register FLOAT size2 = X2[M2+1];
  return (size1 > size2) ? 1 : (size1 < size2) ? -1 : 0;
}

/* sort by decreasing order of map size */
static int CmapSizeDec(register Cmap **p1, register Cmap **p2)
{
  register int M1 = p1[0]->numsite[0];
  register int M2 = p2[0]->numsite[0];
  register FLOAT *X1 = p1[0]->site[0];
  register FLOAT *X2 = p2[0]->site[0];
  register FLOAT size1 = X1[M1+1];
  register FLOAT size2 = X2[M2+1];
  return (size1 < size2) ? 1 : (size1 > size2) ? -1 : 0;
}

/* sort by increasing order of number of sites */
static int CmapSitesInc(register Cmap **p1, register Cmap **p2)
{
  return p1[0]->numsite[0] - p2[0]->numsite[0];
}

/* sort by decreasing order of number of sites */
static int CmapSitesDec(register Cmap **p1, register Cmap **p2)
{
  return p2[0]->numsite[0] - p1[0]->numsite[0];
}

/* sort in increasing order of paired field */
static int CmapPairedInc(register Cmap **p1, register Cmap **p2)
{
  return p1[0]->paired - p2[0]->paired;
}

/* sort by increasing order of id */
static int CmapIdInc(register Cmap **p1, register Cmap **p2)
{
  return (p1[0]->id > p2[0]->id) ? 1 : (p1[0]->id < p2[0]->id) ? -1 : 0;
}

/* sort by increasing order of RunIndex and ScanNumber */
static int CmapRunIndexScanNumberInc(register Cmap **p1, register Cmap **p2)
{
  return (p1[0]->RunIndex > p2[0]->RunIndex) ? 1 : (p1[0]->RunIndex < p2[0]->RunIndex) ? -1 : 
    (p1[0]->ScanNumber > p2[0]->ScanNumber) ? 1 : (p1[0]->ScanNumber < p2[0]->ScanNumber) ? -1 : 
    (p1 > p2) ? 1 : (p1 < p2) ? -1 : 0;// preserves previous order in case of ties
}

/* sort by increasing order of RunIndex and ScanNumber, placing molecules under CmapSortMinLen kb last in each group with same RunIndex,ScanNumber*/
static int CmapRunIndexScanNumberIncMinLen(register Cmap **pp1, register Cmap **pp2)
{
  Cmap *p1 = pp1[0];
  Cmap *p2 = pp2[0];
  int M1 = p1->numsite[0];
  int M2 = p2->numsite[0];
  FLOAT len1 = p1->site[0][M1+1];
  FLOAT len2 = p2->site[0][M2+1];
  int small1 = (len1 > CmapSortMinLen) ? 0 : 1;
  int small2 = (len2 > CmapSortMinLen) ? 0 : 1;

  return (p1->RunIndex > p2->RunIndex) ? 1 : (p1->RunIndex < p2->RunIndex) ? -1 : 
    (p1->ScanNumber > p2->ScanNumber) ? 1 : (p1->ScanNumber < p2->ScanNumber) ? -1 : 
    (small1 > small2) ? 1 : (small1 < small2) ? -1 :
    (pp1 > pp2) ? 1 : (pp1 < pp2) ? -1 : 0;
}

static char *refmap_filename = 0;

extern double mtime();
extern double wtime();

extern int PmapLenDec(Cmap **p1, Cmap **p2);

double origSplitSegDup;
int PairMergeIDrenumbered = 0;

int PairMergeCmapOutput = 0;

int gargc;
char ** gargv;

int main(int argc, char **argv)
{
#if 0
#ifdef __AVX512F__
  printf("__AVX512F__ = %d\n",__AVX512F__);
  fflush(stdout);
#endif
#ifdef __AVX512BW__
  printf("__AVX512BW__ = %d\n",__AVX512BW__);
  fflush(stdout);
#endif
#ifdef __AVX512CD__
  printf("__AVX512CD__ = %d\n",__AVX512CD__);
  fflush(stdout);
#endif
#ifdef __AVX512DQ__
  printf("__AVX512DQ__ = %d\n",__AVX512DQ__);
  fflush(stdout);
#endif
#ifdef __AVX512VL__
  printf("__AVX512VL__ = %d\n",__AVX512VL__);
  fflush(stdout);
#endif
#ifdef __AVX2__
  printf("__AVX2__ = %d\n",__AVX2__);
#endif
#ifdef __AVX__
  printf("__AVX__ = %d\n",__AVX__);
#endif

  printf("USE_MIC= %d, USE_SSE= %d, USE_AVX= %d, USE_AVX512= %d\n", USE_MIC, USE_SSE, USE_AVX, USE_AVX512);
  fflush(stdout);
#endif

  int avx512 = 0, avx2 = 0, avx = 0, sse3 = 0;
#if __AVX512F__ && __AVX512BW__  && __AVX512CD__ && __AVX512DQ__ // && __AVX512VL__
  avx512 = 1;
#endif
#if __AVX2__
  avx2 = 1;
#endif  
#if __AVX__
  avx = 1;
#endif
#if __SSE3__
  sse3 = 1;
#endif

#if !USE_MIC
#ifndef __ARM_ARCH
  extern void check_cpu(int argc, char** argv,int sse3,int avx,int avx2,int avx512);
  check_cpu(argc,argv,sse3,avx,avx2,avx512);
#else
  printf("ARM target\n");
  fflush(stdout);
#endif // ifndef __ARM_ARCH
#endif // !USE_MIC

  START(main);
  (void) mtime();
  (void) wtime();

  if(VERB>=2){
    printf("Start of main():time=%0.6f,wall time=%0.6f\n",mtime(),wtime());
    dumpmemmap();
    fflush(stdout);
  }

  #ifdef WIN32
  _set_fmode(_O_BINARY);
  #endif

  assert(sizeof(int) == 4);
  assert(sizeof(size_t) == 8);
  assert(sizeof(long long) == 8);
  assert(sizeof(unsigned long long) == 8);

#if USE_MIC
  /* mark process as a memory hog for OOM killer in case we run out of memory : kill this process first */
  linux_tweak("/proc/self/oom_adj", "15");
  linux_tweak("/proc/self/oom_score_adj", "1000");
#else
  /* try to protect this process from OOM killer (might help with RESCALE, by encouraging OOM killer to go after rogue file transfer job first) */
  linux_tweak("/proc/self/oom_adj", "-17");
  linux_tweak("/proc/self/oom_score_adj", "-1000");
#endif
  
  if(!mallopt(M_MMAP_THRESHOLD, 128*1024)){
    printf("mallopt(M_MMAP_THRESHOLD,128*1024) failed\n");
    fflush(stdout);
  }
  if(!mallopt(M_MMAP_MAX, 256*1024)){
    printf("mallopt(M_MMAP_MAX,256*1024) failed\n");
    fflush(stdout);
  }

  pid = getpid();/* save process id : used to create unique file names */

  char *qt;

  assert(ENDEXTEND>0.0);

  gargc = argc;
  gargv = argv;

  /* initialize parameter arrays */
  for(int i= 0; i < MAXCOLOR; i++){
    res[i] = 3.5;
    resSD[i] = 0.75;

    FP[i] = 1.0;
    FN[i] = 0.15;
    SF[i] = 0.10;
    SD[i] = 0.15;
    SR[i] = 0.0;
    SE[i] = 0.0;
  }
  MA_mean = 0.0;
  MA_SF = 1.0;
  MA_SD = 0.10;

  int svlvl = 0; /*local for sv--will set globals below*/
  bool setoutlier = false, setextend = false, setbiaswt = false, seta = false, sett = false; /*local for seeing if outlier, extend, biaswt, A is set on command line*/
  char* scaf_filename = 0; /*no need for this to be global*/
  char* bed_filename = 0;
  char* mask_filename = 0;
  bool simpleRepeat = false, simpleRepeatStandalone = false;
  int simpleRepeatMinEle = 5, simpleRepeatFilter = 0;
  float simpleRepeatTolerance = 0.1;
  /* parse command line arguments : modifies global parameter variables */
  argv++,argc--;

  while(argc > 0){
    //    printf("argc=%d, argv[0]= %s\n",argc,argv[0]);
    //    fflush(stdout);

    if(!strcmp(argv[0],"-help")){
      printf("usage:\n%s%s\n",usage1, usage2);
      fflush(stdout);
      exit(0);
    }
    if(argv[0][0] != '-'){
      printf("invalid arg: \"%s\"(expecting an arg starting with '%c' (ASCII %d), instead of with '%c' (ASCII %d))\n",argv[0],'-','-',argv[0][0],argv[0][0]);
      if(argc > 1){
	printf("Subsequent args are:");
	for(int i = 1; i < argc; i++)
	  printf(" \"%s\"", argv[i]);
	printf("\n");
      } else
	printf("See last arg on line\n");
      printf("Use -help for list of options\n");
      fflush(stdout);
      exit(1);
    }

    switch(argv[0][1]){
    case 'v':
      if(!(strcmp(argv[0],"-version"))){
	printversion(stdout);
	fflush(stdout);
	exit(0);
      }
      break;

    case 'A':
      if(!strcmp(argv[0],"-AddFragsWT")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	AddFragsWT = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || AddFragsWT < 0 || AddFragsWT > 1){
	  printf("%s value is invalid (must be 0 or 1 ):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-AlignRes")){
	AlignRes = 1;
	if(argc >= 2 && !(argv[1][0]=='-' && !isdigit(argv[1][1]))){
	  AlignResMult = strtod(argv[1],&qt);
	  if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || AlignResMult < 1.0){
	    printf("%s value has invalid syntax or value (must be >= 1.0):%s\n",argv[0],argv[1]);
	    fflush(stdout); exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-A")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by positive number\n",argv[0]);
	  fflush(stdout); exit(1);
	}

	seta = true;
	AlignedSiteThreshold = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || AlignedSiteThreshold <= 0 ){
	  printf("%s value is invalid (must be positive integer):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-AE")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by positive number\n",argv[0]);
	  fflush(stdout); exit(1);
	}

	AlignedSiteThresholdE = strtol(argv[1],&qt,10);
	AEinit = 1;
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || AlignedSiteThresholdE <= 0 ){
	  printf("%s value is invalid (must be positive integer):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-AS")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by positive number\n",argv[0]);
	  fflush(stdout); exit(1);
	}

	AlignedSiteThresholdS = strtol(argv[1],&qt,10);
	ASinit = 1;
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || AlignedSiteThresholdS <= 0 ){
	  printf("%s value is invalid (must be positive integer):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;

    case 'B':
      if(!strcmp(argv[0],"-BigMapsFirst")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PairMergeBigMapsFirst = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PairMergeBigMapsFirst > 1){
	  printf("%s 1st arg is invalid : should be 0 or 1: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-BreakEndOutlier")){
	if(argc < 5 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-' || argv[4][0] == '-'){
	  printf("%s flag must be followed by 4 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	BreakEndOutlier = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st arg is invalid : should be non-negative number: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	BreakEndOutlierN = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd arg is invalid : should be non-negative number: %s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	BreakEndOutlierE = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd arg is invalid : should be non-negative number: %s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	BreakEndOutlierEN = strtol(argv[4],&qt,10);
	if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	  printf("%s 4th arg is invalid : should be non-negative number: %s\n",argv[0],argv[4]);
	  fflush(stdout);exit(1);
	}
	
	argv += 5;
	argc -= 5;
	continue;
      }

      if(!strcmp(argv[0],"-BreakOutlier")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s flag must be followed by 3 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	BreakOutlier = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st arg is invalid : should be non-negative number: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	BreakOutlierN = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd arg is invalid : should be non-negative number: %s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	BreakOutlierI = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd arg is invalid : should be non-negative number: %s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	argv += 4;
	argc -= 4;
	continue;
      }

      if(!strcmp(argv[0],"-BNXnoRunData")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	BNXnoRunData = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || BNXnoRunData < 0 || BNXnoRunData > 1){
	  printf("%s value is invalid (must be 0 or 1 ):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-BestRef")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	BestRef = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || BestRef < 0 || BestRef > 1){
	  printf("%s value is invalid (must be 0 or 1 ):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-BestRefWT")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	BestRefWT = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || BestRefWT < 0 || BestRefWT > 1){
	  printf("%s value is invalid (must be 0 or 1 ):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-BestRefExt")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	BestRefExt = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	if(argc >= 3){
	  BestRefExtMargin = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value is invalid:%s\n",argv[0],argv[2]);
	    fflush(stdout);
	    exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-BestRefOutlier")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	BestRefOutlier = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s numerical value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-BestRefMargin")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	BestRefMargin = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-BestRefPV")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by 1-2 non-negative numbers\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	BestRefPV = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || BestRefPV < 0 || BestRefPV > 1){
	  printf("%s 1st value is invalid (must be 0 or 1 ):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  BestRefPV_MaxConf = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value is invalid (must be non-negative value):%s\n",argv[0],argv[2]);
	    fflush(stdout); exit(1);
	  }

	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-BppOutput")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	scale_output = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || scale_output < 0 || scale_output > 1){
	  printf("%s value is invalid (must be 0 or 1 ):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-Beta")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);  exit(1);
	}
	Beta = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	Betainit = 1;
	argv += 2;
	argc -= 2;
	continue;
      }
      break; 

    case 'C':
      if(!strcmp(argv[0],"-CoverageIntRcmap")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	CoverageIntRcmap = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || CoverageIntRcmap > 1){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-CutFlip")){
	LogPvThreshold4 = LogPvThreshold3 = -1.0;
	ScoreThreshold4 = ScoreThreshold3 = -1.0;

	if(argc < 5 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-' || argv[4][0] == '-'){
	  printf("%s flag must be followed by 4-9 non-negative numbers\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	CutFlip = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	
	CutFlipMaxErr = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative integer):%s\n",argv[0],argv[2]);
	  fflush(stdout); exit(1);
	}

	CutFlipExpand = strtol(argv[3],&qt,10);
	if(qt==argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd value is invalid (must be non-negative integer):%s\n",argv[0],argv[3]);
	  fflush(stdout); exit(1);
	}

	CutFlipEndExpand = strtol(argv[4],&qt,10);
	if(qt==argv[4] || !(*qt==0 || isspace(*qt))){
	  printf("%s 4th value is invalid (must be non-negative integer):%s\n",argv[0],argv[4]);
	  fflush(stdout); exit(1);
	}

	CutFlipEndFilter = CutFlip;

	if(argc >= 6 && argv[5][0] != '-'){
	  double Pvalue = strtod(argv[5],&qt);
	  if(qt==argv[5] || !(*qt==0 || isspace(*qt)) || Pvalue > 1.0){
	    printf("%s 5th value is invalid (must be between 0 and 1):%s\n",argv[0],argv[5]);
	    fflush(stdout); exit(1);
	  }
	  LogPvThreshold3 = -log(Pvalue)/log(10.0);

	  if(argc >= 7 && !(argv[6][0] == '-' && !isdigit(argv[6][1]))){
	    ScoreThreshold3 = strtod(argv[6],&qt);
	    if(qt==argv[6] || !(*qt==0 || isspace(*qt))){
	      printf("%s 6thh value is invalid (must be number):%s\n",argv[0],argv[6]);
	      fflush(stdout); exit(1);
	    }

	    if(argc >= 8 && argv[7][0] != '-'){
	      double Pvalue = strtod(argv[7],&qt);
	      if(qt==argv[7] || !(*qt==0 || isspace(*qt)) || Pvalue > 1.0 ){
		printf("%s 7th value is invalid (must be between 0 and 1):%s\n",argv[0],argv[7]);
		fflush(stdout); exit(1);
	      }
	      LogPvThreshold4 = -log(Pvalue)/log(10.0);
	      LogPvThreshold4 = max(LogPvThreshold3,LogPvThreshold4);
	      
	      if(argc >= 9 && !(argv[8][0] == '-' && !isdigit(argv[8][1]))){
		ScoreThreshold4 = strtod(argv[8],&qt);
		if(qt==argv[8] || !(*qt==0 || isspace(*qt))){
		  printf("%s 8th value is invalid (must be number):%s\n",argv[0],argv[8]);
		  fflush(stdout); exit(1);
		}
		ScoreThreshold4 = max(ScoreThreshold3,ScoreThreshold4);

		if(argc >= 10 && argv[9][0] != '-'){
		  CutFlipEndFilter = strtol(argv[9],&qt,10);
		  if(qt==argv[9] || !(*qt==0 || isspace(*qt))){
		    printf("%s 9th value is invalid (must be a non-negative integer):%s\n",argv[0],argv[9]);
		    fflush(stdout); exit(1);
		  }

		  argv++;
		  argc--;
		}
		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 5;
	argc -= 5;
	continue;
      }
      if(!strcmp(argv[0],"-CutFlipBC")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0, 1 or 2\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CutFlipBC = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || CutFlipBC > 2){
	  printf("%s value is invalid (must be 0, 1 or 2):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-CutFlipSkip")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CutFlipSkip = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-CutFlipInvDup")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by a 2 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CutFlipDupExpand = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	CutFlipMinErr = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || CutFlipMinErr > 1.0){
	  printf("%s 2nd value is invalid (must be a non-negative number <= 1.0):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-CutFlipMerge")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CutFlipMerge = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || CutFlipMerge > 1.0){
	  printf("%s value is invalid (must be non-negative number between 0 and 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-CutFlipFilter")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1 to 4 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CutFlipMaxOverlap = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid (must be non-negative number):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  CutFlipMaxMiss = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value is invalid (must be non-negative number):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 4 && argv[3][0] != '-'){
	    CutFlipOverlapFilter = strtod(argv[3],&qt);
	    if(qt==argv[3] || !(*qt==0 || isspace(*qt))){
	      printf("%s 3rd value is invalid (must be non-negative number):%s\n",argv[0],argv[3]);
	      fflush(stdout);exit(1);
	    }
	    if(argc >= 5 && argv[4][0] != '-'){
	      CutFlipFilterConf = strtod(argv[4],&qt);
	      if(qt==argv[4] || !(*qt==0 || isspace(*qt))){
		printf("%s 4th value is invalid (must be non-negative number):%s\n",argv[0],argv[4]);
		fflush(stdout);exit(1);
	      }

	      argv++;
	      argc--;
	    }

	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-CIDrenumber")){
	if(argc < 2){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);  exit(1);
	}
	PairMergeIDrenumber = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || PairMergeIDrenumber < 0 || PairMergeIDrenumber > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-CloneTrim")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);  exit(1);
	}
	CloneTrim = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || CloneTrim > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-CloneLimitLen")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by a Length (in kb) and Labels (integer)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CloneLimitLen = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative number):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	CloneLimitN = strtol(argv[2],&qt,10);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-ChimQuality")){
	ChimQuality = 1;
	HapMapWT = 0;
	if(argc >= 2 && argv[1][0] != '-'){
	  ChimQuality = strtol(argv[1],&qt,10);
	  if(argv[1]==qt || !(*qt==0 || isspace(*qt))){
	    printf("%s value %s is not valid (must be non-negative integer)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 3 && argv[2][0] != '-'){
	    HapMapWT = strtol(argv[2],&qt,10);
	    if(argv[2]==qt || !(*qt==0 || isspace(*qt))){
	      printf("%s value %s is not valid (must be non-negative integer)\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-ChimQualityFix")){
	ChimQualityFix = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-CmapSNR")){
	CmapSNR = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-CovTrim")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	COVERAGE_TRIM = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || COVERAGE_TRIM <= 0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-CovTrimLen")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by one or more numbers\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	COVERAGE_TRIM_LEN = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || COVERAGE_TRIM_LEN < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	if(argc > 2 && argv[2][0] != '-'){/* parse alternate values */
	  mCNT = 0;
	  mCOVERAGE_TRIM_LEN[mCNT++] = COVERAGE_TRIM_LEN;
	  
	  int T = 2;
	  for(; argc > T && argv[T][0] != '-'; T++){
	    double val = strtod(argv[T], &qt);
	    if(qt== argv[T] || !(*qt==0 || isspace(*qt)) || val <= mCOVERAGE_TRIM_LEN[mCNT-1]){
	      printf("%s %d'th value is invalid (must be >= previous value %0.3f): %s\n",argv[0],T,mCOVERAGE_TRIM_LEN[mCNT-1],argv[T]);
	      fflush(stdout);exit(1);
	    }
	    if(mCNT >= MAX_CHIM){
	      if(T == mCNT + 1){
		printf("%s : number of values exceeds maximum of %d values (see MAX_CHIM in constants.h) : ignoring additionl values\n", argv[0], MAX_CHIM);
		fflush(stdout);
	      }
	      continue;
	    } 

	    mCOVERAGE_TRIM_LEN[mCNT++] = val;
	  }

	  argv += T-2;
	  argc -= T-2;
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-CovNorm")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	CovNorm = 1;
	CovLambda = strtod(argv[1],&qt);
	if(qt==argv[1] || CovLambda <= 0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'D':
      if(!strcmp(argv[0],"-D")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	AlignedLabelDensityThreshold = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid (must be number a non-negative number):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'E':
      if(!strcmp(argv[0],"-Erefine")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by two numbers\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	RefineEndOutlierThreshold = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || RefineEndOutlierThreshold < 0 ){
	  printf("%s 1st value is invalid (must be 0, 1 or 2):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	RefineEndOutlierWt = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || RefineEndOutlierWt < 0.0){
	  printf("%s 2nd value is invalid (must be non-negative value):%s\n",argv[0],argv[2]);
	  fflush(stdout); exit(1);
	}
	RefineEndOutlierLen = 0.0;
	if(argc >= 4 && argv[3][0] != '-'){
	  RefineEndOutlierLen = strtod(argv[3],&qt);
	  if(qt==argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value is invalid (must be non-negative size in kb):%s\n",argv[0],argv[3]);
	    fflush(stdout); exit(1);
	  }
	  if(argc >= 5 && argv[4][0] != '-'){
	    RefineEndOutlierN = strtol(argv[4],&qt,10);
	    if(qt==argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th value is invalid (must be non-negative integer):%s\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }
	    if(argc >= 6 && argv[5][0] != '-'){
	      RefineEndOutlierType = strtol(argv[5],&qt,10);
	      if(qt==argv[5] || !(*qt==0 || isspace(*qt)) || RefineEndOutlierType > 1){
		printf("%s 5th value is invalid (must be 0 or 1):%s\n",argv[0],argv[5]);
		fflush(stdout);exit(1);
	      }
	      
	      argv++;
	      argc--;
	    }

	    argv++;
	    argc--;
	  }

	  argv++;
	  argc--;
	}
	
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-E")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	AlignedEndOutlierThreshold = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || AlignedEndOutlierThreshold < 0 ){
	  printf("%s value is invalid (must be +ve integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-EndScoreThresh")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0]=='-'){
	  printf("%s flag must be followed by two numbers\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	EndScoreThresh = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value S1 is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	EndScoreThresh = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s value S2 is invalid:%s\n",argv[0],argv[2]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-ErrDist")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	err_dist = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || err_dist < 0.0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-EndTrimRatio")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by number from 0.0 to 1.0\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	EndTrimRatio = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || EndTrimRatio > 1.0 ){
	  printf("%s value is invalid (must be a number from 0.0 to 1.0):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-EndTrimCQ")){
	EndTrimCQ = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  EndTrimCQ = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || EndTrimCQ > 1){
	    printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	    fflush(stdout); exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-EndTrimInternal")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by 2 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	EndTrimCov2 = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid (must non-negative number)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	EndTrimLen2 = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid (must non-negative number)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-EndTrim")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);  exit(1);
	}
	EndTrimCov = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || EndTrimCov < 0){
	  printf("%s 1st value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  EndTrimFrac = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || EndTrimFrac >= 1.0){
	    printf("%s 2nd value is invalid (must be non-negative number below 1.0):%s\n",argv[0],argv[2]);
	    fflush(stdout); exit(1);
	  }

	  if(argc >= 4 && argv[3][0] != '-'){
	    EndTrimCnt = strtod(argv[3],&qt);

	    if(qt==argv[3] || !(*qt==0 || isspace(*qt))){
	      printf("%s 3rd value is invalid (must be non-negative number):%s\n",argv[0],argv[3]);
	      fflush(stdout); exit(1);
	    }

	    if(argc >= 5 && argv[4][0] != '-'){
	      EndTrimWt = strtod(argv[4],&qt);

	      if(qt==argv[4] || !(*qt==0 || isspace(*qt))){
		printf("%s 4th value is invalid (must be non-negative number):%s\n",argv[0],argv[4]);
		fflush(stdout); exit(1);
	      }

	      argv++;
	      argc--;
	    }

	    argv++;
	    argc--;
	  }

	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-EndTrimNoRefine")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	EndTrimNoRefine = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || EndTrimNoRefine < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'F':
      if(!strcmp(argv[0],"-F")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	AlignedFractionThreshold = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || AlignedFractionThreshold > 1.0 ){
	  printf("%s value is invalid (must be number from 0.0 to 1.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FilterXmap")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	FilterXmap = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || FilterXmap > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FilterXmapWT")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	FilterXmapWT = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || FilterXmapWT > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FilterXmapSplit")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	FilterXmapSplit = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || FilterXmapSplit > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FastSnpsStartup")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by positive integer\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	FastSnpsStartup = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || FastSnpsStartup > 10){
	  printf("%s 1st value is invalid (must be non-negative integer <= 10):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FastSnps")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative integer\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	FastSnps = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid: must be non-negative integer :%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FastIndel")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative integer\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	FastIndel = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid: must be non-negative integer :%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FragileQuality")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s must be followed by 3 non-negative values\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	FragileQualityEnd = strtod(argv[1],&qt);
	if(argv[1]==qt || !(*qt==0 || isspace(*qt)) || FragileQualityEnd <= 0.0){
	  printf("%s 1st arg must be +ve number:\n%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	FragileQualityFN = strtol(argv[2],&qt,10);
	if(argv[2]==qt || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd arg must be non-negative integer:\n%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	FragileQualityFP = strtol(argv[3],&qt,10);
	if(argv[3]==qt || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd arg must be non-negative integer:\n%s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-FN_Refine")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by 1 to 4 non-negative values\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	FN_Refine = strtod(argv[1],&qt);
	if(argv[1]==qt || !(*qt==0 || isspace(*qt)) || FN_Refine < 0.0){
	  printf("The %s value %s is not valid (must be non-negative value)\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  FN_Ext = strtod(argv[2],&qt);
	  if(argv[2]==qt || !(*qt==0 || isspace(*qt)) || FN_Ext < 0.0){
	    printf("The 2nd %s value of %s is not valid (must be non-negative value)\n",argv[0],argv[2]);
	    fflush(stdout);  exit(1);
	  }

	  if(argc >= 4 && argv[3][0] != '-'){
	    FNmin = strtod(argv[3],&qt);
	    if(argv[3]==qt || !(*qt==0 || isspace(*qt)) || FNmin < 0.0){
	      printf("The 3rd %s value of %s is not valid (must be non-negative value)\n",argv[0],argv[3]);
	      fflush(stdout);  exit(1);
	    }

	    if(argc >= 5 && argv[4][0] != '-'){
	      FNminExt = strtod(argv[4],&qt);
	      if(argv[4]==qt || !(*qt==0 || isspace(*qt)) || FNminExt < 0.0){
		printf("The 4th %s value of %s is not valid (must be non-negative value)\n",argv[0],argv[4]);
		fflush(stdout);  exit(1);
	      }

	      argv++;
	      argc--;
	    }

	    argv++;
	    argc--;
	  }

	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FP")){
	if(!colors){
	  printf("-colors must be specified before %s\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	if(argc < 2){
	  printf("%s must be followed by at least 1 value (up to one per color)\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	int c = 0;
	for(; c < MAXCOLOR; c++){
	  if(argc <= 1+c || (c >= 1 && argv[c+1][0] == '-' && argv[c+1][1] && !isdigit(argv[c+1][1]) && argv[c+1][1] != '.'))
	    break;
	  FP[c] = strtod(argv[c+1],&qt);
	  if(argv[c+1]==qt || !(*qt==0 || isspace(*qt)) || FP[c] < 0.0){
	    if(colors > 1)
	      printf("The %d. %s value %s is not valid (must be non-negative value)\n",c+1,argv[0],argv[c+1]);
	    else
	      printf("The %s value %s is not valid (must be non-negative value)\n",argv[0],argv[c+1]);
	    fflush(stdout);
	    exit(1);
	  }
	}
	for(int nc = c; nc < MAXCOLOR; nc++)
	  FP[nc] = FP[c-1];
	argv += 1+c;
	argc -= 1+c;
	continue;
      }
      if(!strcmp(argv[0],"-FN")){
	if(!colors){
	  printf("-colors must be specified before %s\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	if(argc < 2){
	  printf("%s must be followed by at least 1 value (up to one per color)\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	int c = 0;
	for(; c < MAXCOLOR; c++){
	  if(argc <= 1+c || (c >= 1 && argv[c+1][0] == '-' && argv[c+1][1] && !isdigit(argv[c+1][1]) && argv[c+1][1] != '.'))
	    break;
	  FN[c] = strtod(argv[c+1],&qt);
	  if(argv[c+1]==qt || !(*qt==0 || isspace(*qt)) || FN[c] < 0.0){
	    if(colors > 1)
	      printf("The %d. %s value %s is not valid\n",c+1,argv[0],argv[c+1]);
	    else
	      printf("The %s value %s is not valid\n",argv[0],argv[c+1]);
	    fflush(stdout);
	    exit(1);
	  }
	}
	for(int nc = c; nc < MAXCOLOR; nc++)
	  FN[nc] = FN[c-1];
	argv += 1+c;
	argc -= 1+c;
	continue;
      }
      break;
    case 'G':
      if(!strcmp(argv[0],"-G")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a positive number\n",argv[0]);
	  fflush(stdout);  exit(1);
	}
	Gbias = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || Gbias <= 0.0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	Gbiasinit = 1;
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'H':
      if(!strcmp(argv[0],"-HashMaxHits")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1 or 3 non-negative integers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	HashMaxHits = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid : should be non-negative integer: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 4 && argv[2][0] != '-' && argv[3][0] != '-'){
	  HashHitsCov = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value is invalid : should be non-negative integer: %s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  HashHitsGenSiz = strtol(argv[3],&qt,10);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value is invalid : should be non-negative integer: %s\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  
	  argv += 2;
	  argc -= 2;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-HScanRange")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	HScanRange = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid : should be non-negative number: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-HapMaxIter")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by two non-negative numbers\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	HapStageIter = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HapStageIter <= 0){
	  printf("%s 1st value is invalid (must be positive integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	HapStageFrac = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative fraction):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-HapPhaseRev")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by two non-negative integers\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	HREVERSAL = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	FINAL_HREVERSAL = strtol(argv[2],&qt,10);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-HapCres")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by two non-negative integers\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	HAP_MAX_DERES = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	DERES_STOP = strtol(argv[2],&qt,10);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-HapFilterStop")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by positive integer\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	FILTER_STOP = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }

      if(!strcmp(argv[0],"-HapFast")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by positive integer\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	HapFast = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-HapFilterLast")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative integer\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	HapFilterLast = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HapFilterLast > 1){
	  printf("%s 1st value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-HapMergeSNPs")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative integer\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	HapMergeSNPs = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HapMergeSNPs > 1){
	  printf("%s 1st value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-HapLocalStop")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative integer\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	HapLocalStop = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HapLocalStop > 3){
	  printf("%s 1st value is invalid (must be between 0 and 3):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  HapLocalScan = strtol(argv[2],&qt,10);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value is invalid (must be non-negative integer):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  if(argc >= 4 && argv[3][0] != '-'){
	    HapLocalScan2 = strtol(argv[3],&qt,10);
	    if(qt==argv[3] || !(*qt==0 || isspace(*qt))){
	      printf("%s 3rd value is invalid (must be non-negative integer):%s\n",argv[0],argv[3]);
	      fflush(stdout);exit(1);
	    }

	    argv++;
	    argc--;
	  }

	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-HapEndTrim")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1-5 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	HapEndTrim = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}

	HapEndTrimMin = 1.0;
	HapEndTrimMaxFrac = 0.50;
	HapEndTrimRatio = 0.20;
	HapEndTrimLen = 99.9;

	if(argc >= 3 && argv[2][0] != '-'){
	  HapEndTrimMin = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value is invalid (must be >= 0):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  if(argc >= 4 && argv[3][0] != '-'){
	    HapEndTrimMaxFrac = strtod(argv[3],&qt);
	    if(qt==argv[3] || !(*qt==0 || isspace(*qt))){
	      printf("%s 3rd value is invalid (must be >= 0):%s\n",argv[0],argv[3]);
	      fflush(stdout);exit(1);
	    }

	    if(argc >= 5 && argv[4][0] != '-'){
	      HapEndTrimRatio = strtod(argv[4],&qt);
	      if(qt==argv[4] || !(*qt==0 || isspace(*qt))){
		printf("%s 4th value is invalid (must be >= 0):%s\n",argv[0],argv[4]);
		fflush(stdout);exit(1);
	      }

	      if(argc >= 6 && argv[5][0] != '-'){
		HapEndTrimLen = strtod(argv[5],&qt);
		if(qt==argv[5] || !(*qt==0 || isspace(*qt))){
		  printf("%s 5th value is invalid (must be > 0):%s\n",argv[0],argv[5]);
		  fflush(stdout);exit(1);
		}

		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-HapMinLL")){
	if(argc < 5 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-' || argv[4][0] == '-'){
	  printf("%s flag must be followed by 4 or 5 positive numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	IntervalEps1 = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	HapIndelEps1 = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be > 0):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	HapLabelEps1 = strtod(argv[3],&qt);
	if(qt==argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd value is invalid:%s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	IntervalEps2 = strtod(argv[4],&qt);
	if(qt==argv[4] || !(*qt==0 || isspace(*qt))){
	  printf("%s 4th value is invalid:%s\n",argv[0],argv[4]);
	  fflush(stdout);exit(1);
	}
	HapIndelEps2 = HapIndelEps1;
	HapLabelEps2 = HapLabelEps1;
	if(argc >= 6 && argv[5][0] != '-'){
	  HapIndelEps2 = strtod(argv[5],&qt);
	  if(qt==argv[5] || !(*qt==0 || isspace(*qt))){
	    printf("%s 5th value is invalid:%s\n",argv[0],argv[5]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 7 && argv[6][0] != '-'){
	    HapLabelEps2 = strtod(argv[6],&qt);
	    if(qt==argv[6] || !(*qt==0 || isspace(*qt))){
	      printf("%s 6th value is invalid:%s\n",argv[0],argv[7]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 5;
	argc -= 5;
	continue;
      }
      if(!strcmp(argv[0],"-HapIndelHiRes")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by positive integer\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	HapIndelHiRes = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HapIndelHiRes <= 0){
	  printf("%s 1st value is invalid (must be positive integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(VERB>=2){
	  if(argc >= 3)
	    printf("argc=3:argv[0]=%s,argv[1]=%s,argv[2]=%s\n",argv[0],argv[1],argv[2]);
	  else
	    printf("argc=%d:argv[0]=%s,argv[1]=%s\n",argc,argv[0],argv[1]);
	  fflush(stdout);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  MAX_DELTA_OVERSAMPLE = strtol(argv[2],&qt,10);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || MAX_DELTA_OVERSAMPLE < 0){
	    printf("%s 2nd value is invalid (must be positive integer):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-HapIndelSkip")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by 2 non-negative integers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	HapIndelDelskip = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	HapUpdatemapSkip = strtol(argv[2],&qt,10);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be >= 0):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 4 && argv[3][0] != '-'){
	  HapSkipEps = strtod(argv[3],&qt);
	  if(qt==argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value is invalid:%s\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 5 && argv[4][0] != '-'){
	    HapSkipEpsPermap = strtod(argv[4],&qt);
	    if(qt==argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th value is invalid:%s\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-HSDrange")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1-2 non-negative numbers\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	HSDrange = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HSDrange < 0.0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	
	if(argc >= 3 && argv[2][0] != '-'){
	  HResMul = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value is invalid (must be >= 0):%s\n",argv[0],argv[2]);
	    fflush(stdout); exit(1);
	  }

	  argv++;
	  argc--;
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-HapIndelMerge")){
	HapIndelMerge = 1;
	RefineFix = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-HapThresh")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by two or more non-negative numbers\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	HapSiteMin = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HapSiteMin < 0){
	  printf("%s 1st value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}

	HapIndelMin = strtol(argv[2],&qt,10);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || HapIndelMin < 0){
	  printf("%s 2nd value is invalid (must be >= 0):%s\n",argv[0],argv[2]);
	  fflush(stdout); exit(1);
	}

	if(argc >= 5 && argv[3][0] != '-' && argv[4][0] != '-'){
	  HapSiteWin = strtol(argv[3], &qt, 10);
	  if(qt==argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value is invalid (must be integer >= 0):%s\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  
	  HapSiteWinSize = strtod(argv[4], &qt);
	  if(qt==argv[4] || !(*qt==0 || isspace(*qt)) || HapSiteWinSize <= 0.0){
	    printf("%s 4th value is invalid (must be number > 0):%s\n",argv[0],argv[4]);
	    fflush(stdout);exit(1);
	  }

	  if(argc >= 6 && argv[5][0] != '-'){
	    HapIndelMinSum = strtod(argv[5], &qt);
	    if(qt==argv[5] || !(*qt==0 || isspace(*qt)) || HapIndelMinSum < 0.0){
	      printf("%s 5th value is invalid (must be number >= 0):%s\n",argv[0],argv[5]);
	      fflush(stdout); exit(1);
	    }

	    if(argc >= 7 && argv[6][0] != '-'){
	      HapIndelMinDensity = strtod(argv[6],&qt);
	      if(qt==argv[6] || !(*qt==0 || isspace(*qt)) || HapIndelMinSum < 0.0){
		printf("%s 6th value is invalid (must be number >= 0):%s\n",argv[0],argv[6]);
		fflush(stdout); exit(1);
	      }
	      
	      if(argc >= 8 && argv[7][0] != '-'){
		HapIndelSkipEndSNP = strtol(argv[7],&qt,10);
		if(qt== argv[7] || !(*qt==0 || isspace(*qt)) || HapIndelSkipEndSNP > 1){
		  printf("%s 7th value is invalid (must be 0 or 1):%s\n",argv[0],argv[7]);
		  fflush(stdout); exit(1);
		}

		if(argc >= 9 && argv[8][0] != '-'){
		  MinHapCovThresh = strtod(argv[8],&qt);
		  if(qt== argv[8] || !(*qt==0 || isspace(*qt))){
		    printf("%s 8th value is invalid (must be non-zero number):%s\n",argv[0],argv[8]);
		    fflush(stdout); exit(1);
		  }

		  if(argc >= 10 && argv[9][0] != '-'){
		    MinHapRelCovThresh = strtod(argv[9],&qt);
		    if(qt == argv[9] || !(*qt==0 || isspace(*qt))){
		      printf("%s 9th value is invalid (must be non-zero number): %s\n",argv[0],argv[9]);
		      fflush(stdout); exit(1);
		    }

		    argv++;
		    argc--;
		  }

		  argv++;
		  argc--;
		}

		argv++;
		argc--;
	      }

	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }

	  argv += 2;
	  argc -= 2;
	}

	argv += 3;
	argc -= 3;
	continue;
      }

      if(!strcmp(argv[0],"-HapSiteUnphased")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by 2 non-negative numbers\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	HapSiteUnphasedSNR = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HapSiteUnphasedSNR < 0.0){
	  printf("%s 1st value is invalid (must be non-negative number):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	HapSiteUnphasedPvalue = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || HapSiteUnphasedPvalue < 0.0 || HapSiteUnphasedPvalue > 1.0){
	  printf("%s 2nd value is invalid (must be between 0 and 1):%s\n",argv[0],argv[2]);
	  fflush(stdout); exit(1);
	}
	
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-HapMinCovPL")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by 2 to 4 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	HapMinCovPL1 = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HapMinCovPL1 < 0.0 || HapMinCovPL1 > 1.0){
	  printf("%s 1st value is invalid (must be a fractional number from 0 to 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	HapMinCovPL2 = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || HapMinCovPL2 < 0.0 || HapMinCovPL2 > 1.0){
	  printf("%s 2nd value is invalid (must be a fractional number from 0 to 1):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	if(HapMinCovPL1 < HapMinCovPL2){
	  printf("%s 1st value %s cannot be smaller than 2nd value %s\n",argv[0],argv[1],argv[2]);
	  fflush(stdout);exit(1);
	}

	HapMinCovL = 0;
	HapMinCovCL = 0.0;

	if(argc >= 4 && argv[3][0] != '-'){
	  HapMinCovL = strtol(argv[3],&qt,10);
	  if(qt==argv[3] || !(*qt==0 || isspace(*qt)) || HapMinCovL > 2){
	    printf("%s 3rd value is invalid (must be 0, 1 or 2):%s\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }

	  if(argc >= 5 && argv[4][0] != '-'){
	    HapMinCovCL = strtod(argv[4],&qt);
	    if(qt==argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th value is invalid (must be non-negative coverager value):%s\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }
	    
	    argv++;
	    argc--;
	  }

	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-HapMinCov")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s flag must be followed by 3-7 non-negative numbers\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	HapMinCoverage = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HapMinCoverage < 0.0){
	  printf("%s 1st value is invalid (must be non-negative number):%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	HapAlleleCoverage = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || HapAlleleCoverage < 0.0 || HapAlleleCoverage > 1.0){
	  printf("%s 2nd value is invalid (must be between 0 and 1):%s\n",argv[0],argv[2]);
	  fflush(stdout);  exit(1);
	}
	HapAlleleProb = strtod(argv[3],&qt);
	if(qt==argv[3] || !(*qt==0 || isspace(*qt)) || HapAlleleProb < 0.0 || HapAlleleProb > 1.0){
	  printf("%s 3rd value is invalid (must be between 0 and 1):%s\n",argv[0],argv[3]);
	  fflush(stdout);  exit(1);
	}

	HapMinCovSize = 0.0;
	HapMinCovPvalue = 0.0;
	HapAlleleCovInt = 0.0;
	HapMinCovPvalueL = 0.0;
	HapMinCovAbs = 0.0;
	HapMinCovAbsM = 99.0;
	HapMinCovN = 999;

	if(argc >= 5 && argv[4][0] != '-'){
	  HapMinCovSize = strtod(argv[4],&qt);
	  if(qt==argv[4] || !(*qt==0 || isspace(*qt)) || HapMinCovSize < 0.0){
	    printf("%s 4th value is invalid (must be 0 or greater):%s\n",argv[0],argv[4]);
	    fflush(stdout); exit(1);
	  }
	  if(argc >= 6 && argv[5][0] != '-'){
	    HapMinCovPvalueL = HapMinCovPvalue = strtod(argv[5],&qt);
	    if(qt==argv[5] || !(*qt==0 || isspace(*qt)) || HapMinCovPvalue < 0.0){
	      printf("%s 5th value is invalid (must be 0 or greater):%s\n",argv[0],argv[5]);
	      fflush(stdout);
	      exit(1);
	    }
	    if(argc >= 7 && argv[6][0] != '-'){
	      HapAlleleCovInt = strtod(argv[6],&qt);
	      if(qt==argv[6] || !(*qt==0 || isspace(*qt)) || HapAlleleCovInt < 0.0){
		printf("%s 6th value is invalid (must be 0 or greater):%s\n",argv[0],argv[6]);
		fflush(stdout); exit(1);
	      }
	      if(argc >= 8 && argv[7][0] != '-'){
		HapMinCovPvalueL = strtod(argv[7],&qt);
		if(qt==argv[7] || !(*qt==0 || isspace(*qt)) || HapMinCovPvalueL < 0.0){
		  printf("%s 7th value is invalid (must be 0 or greater):%s\n",argv[0],argv[7]);
		  fflush(stdout); exit(1);
		}
		if(argc >= 9 && argv[8][0] != '-'){
		  HapMinCovAbs = strtod(argv[8],&qt);
		  if(qt==argv[8] || !(*qt==0 || isspace(*qt))){
		    printf("%s 8th value is invalid (must be coverage value greater than 0.0):%s\n",argv[0],argv[8]);
		    fflush(stdout); exit(1);
		  }
		  if(argc >= 10 && argv[9][0] != '-'){
		    HapMinCovAbsM = strtod(argv[9],&qt);
		    if(qt==argv[9] || !(*qt==0 || isspace(*qt))){
		      printf("%s 9th value is invlid (must non-negative number, Labels per 100kb):%s\n",argv[0],argv[9]);
		      fflush(stdout);exit(1);
		    }
		    if(argc >= 11 && argv[10][0] != '-'){
		      HapMinCovN = strtol(argv[10],&qt,10);
		      if(qt==argv[10] || !(*qt==0 || isspace(*qt))){
			printf("%s 10th value is invlid (must non-negative integer):%s\n",argv[0],argv[10]);
			fflush(stdout);exit(1);
		      }
		      argv++;
		      argc--;
		    }
		    argv++;
		    argc--;
		  }
		  argv++;
		  argc--;
		}
		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-HapSiteRes")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	HapSiteRes = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }

      if(!strcmp(argv[0],"-HapSiteResDelay")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s flag must be followed by a kb size and two positive integers\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	HapSiteResDelay = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be size in kb):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}

	HapSiteResN = strtol(argv[2],&qt,10);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || HapSiteResN <= 0){
	  printf("%s 2nd value is invalid (must be positive integer):%s\n",argv[0],argv[2]);
	  fflush(stdout); exit(1);
	}

	HapSiteResL = strtol(argv[3],&qt,10);
	if(qt==argv[3] || !(*qt==0 || isspace(*qt)) || HapSiteResL <= 0){
	  printf("%s 3rd value is invalid (must be positive integer):%s\n",argv[0],argv[3]);
	  fflush(stdout); exit(1);
	}

	if(argc >= 5 && argv[4][0] != '-'){
	  HapSiteResA = strtol(argv[4],&qt,10);
	  if(qt==argv[4] || !(*qt==0 || isspace(*qt)) || HapSiteResA > 1){
	    printf("%s 4th value is invalid (must be 0 or 1):%s\n",argv[0],argv[4]);
	    fflush(stdout);exit(1);
	  }

	  argv++;
	  argc--;
	}

	argv += 4;
	argc -= 4;
	continue;
      }

      if(!strcmp(argv[0],"-Haplotype")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be following by 2 non-negative numbers (Pvalues)\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	HapSitePvalue = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HapSitePvalue <= 0.0 || HapSitePvalue > 1.0){
	  printf("%s value %s is invalid (must be between 0 and 1)\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	HapIndelPvalue = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || HapIndelPvalue < 0.0 || HapIndelPvalue > 1.0){
	  printf("%s value %s is invalid (must be between 0 and 1)\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	if(argc > 3 && argv[3][0] != '-'){
	  PhasePvalue = strtod(argv[3],&qt);
	  if(qt==argv[3] || !(*qt==0 || isspace(*qt)) || PhasePvalue < 0.0 || PhasePvalue > 1.0){
	    printf("%s value %s is invalid (must be between 0 and 1)\n",argv[0],argv[1]);
	    fflush(stdout);
	    exit(1);
	  }
	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }

      if(!strcmp(argv[0],"-Hash_Bits")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	Hash_Bits = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || Hash_Bits <= 0){
	  printf("%s value is invalid (must be +ve integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }

      if(!strcmp(argv[0],"-HHash_Bits")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	HHash_Bits = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HHash_Bits <= 0){
	  printf("%s value is invalid (must be +ve integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }

      break;

    case 'I':
      if(!strcmp(argv[0],"-IndelScore")){
	if(argc < 6 || argv[1][0] == '-' || (argv[2][0] == '-' && !isdigit(argv[2][1])) || (argv[3][0] == '-' && !isdigit(argv[3][1])) || (argv[4][0] == '-' && !isdigit(argv[4][1])) || argv[5][0] == '-'){
	  printf("%s flag must be followed by 5 numbers (1st and 5th must be non-negative)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	IndelScore = 1;
	
	IndelDelOutlierLambda = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative size in kb):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	IndelDelOutlierType = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt == 0 || isspace(*qt)) || IndelDelOutlierType < -1 || IndelDelOutlierType > 1){
	  printf("%s 2nd value is invalid (must be -1, 0 or 1):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	IndelDelPoutlier = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt == 0 || isspace(*qt)) || IndelDelPoutlier < -1.0 || IndelDelPoutlier > 1.0){
	  printf("%s 3rd value is invalid (must be -1 or fraction from 0.0 to 1.0):%s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	IndelDelPoutlierEnd = strtod(argv[4],&qt);
	if(qt == argv[4] || !(*qt == 0 || isspace(*qt)) || IndelDelPoutlierEnd < -1.0 || IndelDelPoutlierEnd > 1.0){
	  printf("%s 4th value is invalid (must be -1 or fraction from 0.0 to 1.0):%s\n",argv[0],argv[4]);
	  fflush(stdout);exit(1);
	}

	IndelDelLRbias = strtod(argv[5],&qt);
	if(qt == argv[5] || !(*qt == 0 || isspace(*qt))){
	  printf("%s 5th value is invalid (must be non-negative number):%s\n",argv[0],argv[5]);
	  fflush(stdout);exit(1);	  
	}

	MinSFIndelDel = MinSRIndelDel = 0.0;

	if(argc >= 7 && argv[6][0] != '-'){
	  MinSFIndelDel = strtod(argv[6],&qt);
	  if(qt == argv[6] || !(*qt == 0 || isspace(*qt))){
	    printf("%s 6th value is ivnalid (must be non-negative number):%s\n",argv[0],argv[6]);
	    fflush(stdout);exit(1);
	  }

	  if(argc >= 8 && argv[7][0] != '-'){
	    MinSRIndelDel = strtod(argv[7],&qt);
	    if(qt == argv[7] || !(*qt == 0 || isspace(*qt))){
	      printf("%s 7th value is ivnalid (must be non-negative number):%s\n",argv[0],argv[7]);
	      fflush(stdout);exit(1);
	    }
	    
	    argv++;
	    argc--;
	  }

	  argv++;
	  argc--;
	}

	if(IndelDelOutlierLambda == 0.0 && IndelDelOutlierType < 0 && IndelDelPoutlier < 0.0 && IndelDelPoutlierEnd < 0.0 && IndelDelLRbias == 0.0 && MinSFIndelDel == 0.0 && MinSRIndelDel == 0.0)
	  IndelScore = 0;

	argv += 6;
	argc -= 6;
	continue;
      }
      if(!strcmp(argv[0],"-InversionNormalPV")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed 1 or 2 non-negative numbers\n", argv[0]);
	  fflush(stdout);  exit(1);
	}
	double Pvalue = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s 1st value has invalid syntax for float:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(Pvalue <= 0.0 || Pvalue > 1.0){
	  printf("%s 1st value of %0.6e is invalid (must be above 0 and <= 1):%s\n",argv[0], Pvalue, argv[1]);
	  fflush(stdout);exit(1);
	}
	InversionNormalLogPv = -log(Pvalue)/log(10.0);

	if(argc >= 3 && argv[2][0] != '-'){
	  Pvalue = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value has invalid syntax for float:%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  if(Pvalue <= 0.0 || Pvalue > 1.0){
	    printf("%s 2nd value of %0.6e is invalid (must be above 0 and <= 1):%s\n",argv[0], Pvalue, argv[1]);
	    fflush(stdout);exit(1);
	  }
	  InversionNormalLogPv5 = -log(Pvalue)/log(10.0);

	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-InversionOverlapFilter")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by 2 non-negative numbers:%s\n", argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	InvOverlapFilterQryRatio = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || InvOverlapFilterQryRatio > 1.0){
	  printf("%s 1st value is invalid (must be non-negative number <= 1.0:%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	InvOverlapFilterConfDelta = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative number):%s\n",argv[0],argv[2]);
	  fflush(stdout);  exit(1);
	}
	if(argc >= 4 && argv[3][0] != '-'){
	  InvOverlapFilterQryLabRatio = strtod(argv[3],&qt);
	  if(qt==argv[3] || !(*qt==0 || isspace(*qt)) || InvOverlapFilterQryLabRatio > 1.0){
	    printf("%s 3rd value is invalid (must be non-negative number <= 1.0:%s\n",argv[0],argv[3]);
	    fflush(stdout);  exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-Intensity")){
	MapIntensity = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-I")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1-2 non-negative numbers:%s\n", argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	AlignedOutlierThreshold = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative kb value):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}

	AlignedOutlierLabels = 1000;
	if(argc >= 3 && argv[2][0] != '-'){
	  AlignedOutlierLabels = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value is invalid (must be non-negative integer):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;

    case 'K':
      if(!strcmp(argv[0],"-Kmax")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	KMAX = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || KMAX <= 0 || KMAX > K_MAX){
	  printf("%s value is invalid (must be between 1 and %d):%s\n",argv[0],K_MAX,argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-K")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	Kbias = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || Kbias <= 0.0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	Kbiasinit = 1;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-KL")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	KLbias = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	KLbiasinit = 1;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-KF")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	KFbias = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || KFbias < 0.0){
	  printf("%s value is invalid (must be >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	KFbiasinit = 1;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-KRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	KbiasRef = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || KbiasRef <= 0.0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-KLRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);  exit(1);
	}
	KLbiasRef = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-KFRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);  exit(1);
	}
	KFbiasRef = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || KFbiasRef < 0.0){
	  printf("%s value is invalid (must be >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'L':
      if(!strcmp(argv[0],"-LabelDensityType")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	LabelDensityType = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || LabelDensityType > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-LocalRange")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	LocalRange = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || LocalRange > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-LocalTheta")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1 (and optional 2nd value)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	LocalTheta = strtol(argv[1],&qt,10);
	MaxLocalTheta = 0.0;
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || LocalTheta > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  MaxLocalTheta = strtod(argv[1],&qt);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value is invalid (must be non-negative number):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-L")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	AlignedLengthThreshold = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || AlignedLengthThreshold < 0 ){
	  printf("%s value is invalid (must be number >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-LE")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	AlignedLengthThresholdE = strtod(argv[1],&qt);
	LEinit = 1;
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || AlignedLengthThresholdE < 0 ){
	  printf("%s value is invalid (must be number >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-LS")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	AlignedLengthThresholdS = strtod(argv[1],&qt);
	LSinit = 1;
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || AlignedLengthThresholdS < 0 ){
	  printf("%s value is invalid (must be number >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-LRbias")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	LRbias = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || LRbias < 0.0){
	  printf("%s value is invalid (must be >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-LRbiasSwitch")){
	LRbiasSwitch = 1;
	if(!(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-')){
	  printf("%s must be followed by two positive numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	LRbiasLabel = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	LRbiasSize = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value of %s is invalid (must be positive)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      break;
    case 'M':
      if(!strcmp(argv[0],"-MultiMatchesDupes")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s must be followed by 3 non-negative integers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MultiMatchesDupesMinscore = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value of %s is invalid (must be non-negative number)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	MultiMatchesDupesMinscore = max(0.01,MultiMatchesDupesMinscore);// A value of 0 risks an exponential blowup

	MultiMatchesDupesMinscoreRel = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value of %s is invalid (must be non-negative number)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	
	MultiMatchesDupesMMPen = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd value of %s is invalid (must be non-negative number)\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-MultiMatchesUniqueQuery")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by 2 or 3 non-negative integers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MultiMatchesUniqueQuery = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value of %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	MultiMatchesUniqueQueryInv = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value of %s is invalid (must be non-negative integer)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	MultiMatchesUniqueQueryTrim = 0;
	if(argc >= 4 && argv[3][0] != '-'){
	  MultiMatchesUniqueQueryTrim = strtol(argv[3],&qt,10);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value of %s is invalid (must be non-negative integer)\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }

	  argv++;
	  argc--;
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-MATCHGROUP_TRIM")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0,1,2 or 3\n",argv[0]);
	  fflush(stdout);exit(1);
	}	
	smap_MATCHGROUP_TRIM = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || smap_MATCHGROUP_TRIM > 3){
	  printf("%s value %s is invalid : must be 0,1,2 or 3\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MATCHGROUP_PARTIAL_TRIM")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0,1 or 2\n",argv[0]);
	  fflush(stdout);exit(1);
	}	
	smap_MATCHGROUP_PARTIAL_TRIM = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || smap_MATCHGROUP_PARTIAL_TRIM > 2){
	  printf("%s value %s is invalid : must be 0,1 or 2\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinorIsBigger")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	checkMinorIsBigger = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || checkMinorIsBigger > 1){
	  printf("%s value %s is invalid : must be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxCoverage")){
	if(!(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-')){
	  printf("%s must be followed by two non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxCoverageAN = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value of %s is invalid (must be positive number)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	MaxCoverageLenAN = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value of %s is invalid (must be positive number)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-MinSFSwitch")){
	MinSFSwitch = 1;
	if(!(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-')){
	  printf("%s must be followed by two sf values in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinSFLabel = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	MinSFSize = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value of %s is invalid (must be positive)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	if(argc >= 4 && argv[3][0] != '-'){
	  MinSFLabelIter = strtol(argv[3],&qt,10);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd vale of %s is invalid (must be non-negative integer)\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  
	  if(argc >= 5 && argv[4][0] != '-'){
	    MinSFSizeIter = strtol(argv[4],&qt,10);
	    if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th vale of %s is invalid (must be non-negative integer)\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }
	    
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-MinSRSwitch")){
	MinSRSwitch = 1;
	if(!(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-')){
	  printf("%s must be followed by two sf values in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinSRLabel = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	MinSRSize = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value of %s is invalid (must be positive)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	if(argc >= 4 && argv[3][0] != '-'){
	  MinSRLabelIter = strtol(argv[3],&qt,10);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd vale of %s is invalid (must be non-negative integer)\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  
	  if(argc >= 5 && argv[4][0] != '-'){
	    MinSRSizeIter = strtol(argv[4],&qt,10);
	    if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th vale of %s is invalid (must be non-negative integer)\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }
	    
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }

      if(!strcmp(argv[0],"-MinSF_Refine")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by two non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	
	MinSF_Refine = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative value in kb):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	MinSF_Ext = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative value in kb):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-MinSR_Refine")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by two non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MinSR_Refine = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative value):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	MinSR_Ext = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative value):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-MaxSwitchIter")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	HapMaxSwitchIter = strtol(argv[1],&qt, 10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be integer >= 1)\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}	

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxDeltaIter")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	HapMaxDeltaIter = strtol(argv[1],&qt, 10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be integer >= 1)\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}	

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxDelta")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by positive size in kb\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxHapDelta = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxHapDelta > 350.0){
	  printf("%s 1st value is invalid (must be positive number <= 350.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	MaxHapDelta *= 0.5;// half of largest Het Indel size

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinKB")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number (max refined interval size in kb)\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	ZERO_MINKB = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative number)\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}	

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MergeMG")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number (max gap size in kb)\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	MergeMG = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative number)\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}	

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinDupSites")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	MinDupSites = strtol(argv[1],&qt, 10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || MinDupSites <= 0){
	  printf("%s value %s is invalid (must be integer >= 1)\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}	

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxPalindrome")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s flag must be followed by 3 or 4 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MaxPalindrome = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxPalindrome <= 0.0){
	  printf("%s : 1st value should be greater than 0: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	MaxPalindromeN = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || MaxPalindromeN <= 0){
	  printf("%s : 2nd value should be greater than 0: %s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	MaxPalindromeSL = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || MaxPalindromeSL <= 0.0){
	  printf("%s : 3rd value should be greater than 0: %s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	MaxPalindromeMask = END_NOEXT;
	if(argc >= 5 && argv[4][0] != '-'){
	  MaxPalindromeMask = strtoll(argv[4],&qt,0);
	  if(qt == argv[4] || !(*qt == 0 || isspace(*qt))){
	    printf("%s 4th value is invalid : must be non-negative 63-bit mask, typically 0x%x or 0x0: %s\n", argv[0],END_NOEXT,argv[4]);
	    fflush(stdout);exit(1);
	  }

	  argv++;
	  argc--;
	}
	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-MaxInvPalindrome")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s flag must be followed by 3 or 4 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MaxInvPalindrome = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxInvPalindrome <= 0.0){
	  printf("%s : 1st value should be greater than 0: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	MaxInvPalindromeN = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || MaxInvPalindromeN <= 0){
	  printf("%s : 2nd value should be greater than 0: %s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	MaxInvPalindromeSL = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || MaxInvPalindromeSL <= 0.0){
	  printf("%s : 3rd value should be greater than 0: %s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	MaxInvPalindromeMask = END_NOEXT;
	if(argc >= 5 && argv[4][0] != '-'){
	  MaxInvPalindromeMask = strtoll(argv[4],&qt,0);
	  if(qt == argv[4] || !(*qt == 0 || isspace(*qt))){
	    printf("%s 4th value is invalid : must be non-negative 63-bit mask, typically 0x%x or 0x0: %s\n",argv[0],END_NOEXT,argv[4]);
	    fflush(stdout);exit(1);
	  }
	  
	  argv++;
	  argc--;
	}
	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-MPmem")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MP_MIN_MEM = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MP_MIN_MEM < 0){
	  printf("%s %s : value should be 0 or greater\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  MP_MIN_TIME = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s %s : value should be 0 or greater\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MapScale")){
	MapScale = 1;
	if(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-'){
	  MapScaleDeltaSize = strtod(argv[1],&qt);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MapScaleDeltaSize < 0.0){
	    printf("%s 1st value of %s is invalid (must be >= 0.0)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }

	  MapScaleDelta = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value of %s is invalid (must be non-negative integer)\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  argv += 2;
	  argc -= 2;
	}
	if(MapScaleDeltaSize <= 0.0 || MapScaleDelta <= 0){
	  if(VERB && MapScale){
	    printf("WARNING : -MapScale disabled due to missing or zero <stepsize> and <Steps> values\n");
	    fflush(stdout);
	  }
	  MapScale = 0;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-MinMapwt")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	RefineMinMapwt = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || RefineMinMapwt < 0.0){
	  printf("%s value is invalid (must be >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  UnrefineMinMapwt = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt == 0 || isspace(*qt)) || UnrefineMinMapwt < 0.0){
	    printf("%s 2nd value is invalid (must be >= 0.0):%s\n",argv[0],argv[2]);
	    fflush(stdout);  exit(1);
	  }
	  argv++;
	  argc--;
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MultiMatchesDelta")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MultiMatchesDelta = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || MultiMatchesDelta < 0.0){
	  printf("%s value is invalid (must be >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MultiMatchesFilter")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0,1 or 2\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MultiMatchesFilter = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || MultiMatchesFilter > 2){
	  printf("%s value is invalid (must be 0,1 or 2):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MultiMatchesRev")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MultiMatchesRev = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || MultiMatchesRev > 1){
	  printf("%s value %s is invalid : must be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MultiMatchesTotScore")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by an integer and size in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MultiMatchesTotScore = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || MultiMatchesTotScore > 3){
	  printf("%s 1st arg %s invalid : must be 0, 1, 2 or 3\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	MultiMatchesTotDelta = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt == 0 || isspace(*qt))){
	  printf("%s 2nd arg %s invalid : must size in kb\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-MultiMatches")){
	MultiMatches = 1;
	LogPvThreshold2 = -1.0;
	ScoreThreshold2 = -1.0;
	AlignedSiteThreshold2 = -1;
	AlignedLengthThreshold2 = -1.0;
	AlignedEndOutlierThreshold2 = -1;
	AlignedOutlierThreshold2 = -1.0;
	AlignedOutlierLabels2 = -1;

	if(argc >= 2 && !(argv[1][0] == '-' && !isdigit(argv[1][1]))){
	  MultiMatches = strtol(argv[1],&qt,10);
	  if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MultiMatches < 0){
	    printf("%s 1st value of %s is invalid\n",argv[0],argv[1]);
	    fflush(stdout);
	    exit(1);
	  }
	  if(argc >= 3 && argv[2][0] != '-'){
	    double Pvalue = strtod(argv[2], &qt);
	    if(qt==argv[2] || !(*qt == 0 || isspace(*qt)) || Pvalue > 1.0){
	      printf("%s 2nd value of %s is invalid (must be between 0 and 1)\n",argv[0],argv[2]);
	      fflush(stdout); exit(1);
	    }
	    LogPvThreshold2 = -log(Pvalue)/log(10.0);
	    if(argc >= 4 && argv[3][0] != '-'){
	      ScoreThreshold2 = strtod(argv[3],&qt);
	      if(qt==argv[3] || !(*qt == 0 || isspace(*qt))){
		printf("%s 3rd value of %s is invalid (must be >= 0)\n",argv[0],argv[3]);
		fflush(stdout);	exit(1);
	      }
	      if(argc >= 5 && argv[4][0] != '-'){
		AlignedSiteThreshold2 = strtol(argv[4],&qt,10);
		if(qt==argv[4] || !(*qt == 0 || isspace(*qt))){
		  printf("%s 4th value of %s is invalid (must be >= 0)\n",argv[0],argv[4]);
		  fflush(stdout);exit(1);
		}
		if(argc >= 6 && argv[5][0] != '-'){
		  AlignedLengthThreshold2 = strtod(argv[5],&qt);
		  if(qt==argv[5] || !(*qt == 0 || isspace(*qt))){
		    printf("%s 5th value of %s is invalid (must be >= 0)\n",argv[0],argv[5]);
		    fflush(stdout);exit(1);
		  }
		  if(argc >= 7 && argv[6][0] != '-'){
		    AlignedEndOutlierThreshold2 = strtol(argv[6],&qt,10);
		    if(qt==argv[6] || !(*qt == 0 || isspace(*qt))){
		      printf("%s 6th value of %s is invalid (must be >= 0)\n",argv[0],argv[6]);
		      fflush(stdout);exit(1);
		    }
		    if(argc >= 8 && argv[7][0] != '-'){
		      AlignedOutlierThreshold2 = strtol(argv[7],&qt,10);
		      if(qt==argv[7] || !(*qt == 0 || isspace(*qt))){
			printf("%s 7th value of %s is invalid (must be >= 0)\n",argv[0],argv[7]);
			fflush(stdout);exit(1);
		      }
		      if(argc >= 9 && argv[8][0] != '-'){
			AlignedOutlierLabels2 = strtol(argv[8],&qt,10);
			if(qt==argv[8] || !(*qt == 0 || isspace(*qt))){
			  printf("%s 8th value of %s is invalid (must be non-negative integer)\n",argv[0],argv[8]);
			  fflush(stdout);exit(1);
			}

			argv++;
			argc--;
		      }
		      argv++;
		      argc--;
		    }
		    argv++;
		    argc--;
		  }
		  argv++;
		  argc--;
		}
		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	if(VERB>=2){
	  printf("MultiMatches=%d,S2=%0.4f,T2=%0.2f,A2=%d,L2=%0.3f,E2=%d,I2=%0.3f\n",
		 MultiMatches,ScoreThreshold2,LogPvThreshold2,AlignedSiteThreshold2,AlignedLengthThreshold2,AlignedEndOutlierThreshold2,AlignedOutlierThreshold2);
	  fflush(stdout);
	}

	continue;
      }
      if(!strcmp(argv[0],"-MaxIntensity")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MaxIntensity = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxIntensity < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MprobevalTest")){
	MprobevalTest = 1;
	RefineStep1 = min(RefineStep1, 0.02);
	RefineStep2 = min(RefineStep2, 0.02);
	RefineStep3 = min(RefineStep3, 0.05);
	RefineStep4 = min(RefineStep4, 0.05);
	RefineStepCnt = max(RefineStepCnt, 14);
        RefineStep5 = min(RefineStep5, 0.02);
        RefineStep6 = min(RefineStep6, 0.05);
	RefineStep7 = min(RefineStep7, 0.02);
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-MinSF")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MinSF = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MinSF < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinSR")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MinSR = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MinSR < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinSE")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MinSE = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MinSE < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinRes")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MinRes = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MinRes < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinResSD")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MinResSD = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MinResSD < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxRes")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MaxRes = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxRes < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxResSD")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MaxResSD = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxResSD < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinSD")){
	if(argc < 2){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MinSD = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid (must be a float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinSD_Mult")){
	if(argc < 2){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MinSD_Mult = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MinSD_Mult < 0.0 || MinSD_Mult > 1.0){
	  printf("%s value is invalid (must be a float between 0.0 and 1.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxSD")){
	if(argc < 2){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MaxSD = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid (must be a float):%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinFP")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);  exit(1);
	}
	MinFP = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MinFP < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxFP")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);  exit(1);
	}
	MaxFP = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxFP < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxFN")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MaxFN = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxFN < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinFN")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MinFN = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MinFN < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxSF")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MaxSF = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxSF < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxSR")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MaxSR = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxSR < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxSE")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MaxSE = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxSE < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MapRate")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MapRate = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MapRate < 0.0){
	  printf("%s value is invalid (must be non-negative):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  double minT = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt == 0 || isspace(*qt)) || minT <= 0.0){
	    printf("%s 2nd value is invalid : must be positive number:%s\n",argv[0],argv[2]);
	    fflush(stdout);
	    exit(1);
	  }
	  MaxLogPvThreshold = -log(minT)/log(10.0);
	  if(argc >= 4 && argv[3][0] != '-'){
	    MaxScoreThreshold = strtod(argv[3],&qt);
	    if(qt==argv[3] || !(*qt == 0 || isspace(*qt)) || MaxScoreThreshold < 0.0){
	      printf("%s 3rd value is invalid : must be non-negative number:%s\n",argv[0],argv[2]);
	      fflush(stdout);
	      exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MHash_Bits")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MHash_Bits = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MHash_Bits <= 0){
	  printf("%s value is invalid (must be +ve integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-Mfast")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	Mfast = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || Mfast < 0 || Mfast > 2){
	  printf("%s value %s is invalid (must be 0,1 or 2)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MultiMode")){
	MultiMode = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-Mprobeval")){
	Mprobeval = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  DELTA_RES = strtod(argv[1],&qt);
	  if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || DELTA_RES < 0.0){
	    printf("%s 1st value is invalid (must be >= 0.0):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 3 && argv[2][0] != '-'){
	    DELTA_REL = strtod(argv[2],&qt);
	    if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || DELTA_REL < 0.0){
	      printf("%s 2nd value is invalid (must be > 0.0):%s\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;

	if(DELTA_RES <= 0.0)
	  Mprobeval = 0;
	continue;
      }
      if(!strcmp(argv[0],"-MinSplitLen")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinSplitLen = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MinSplitLen <= 0.0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-Msave")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	SaveRepeat = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || SaveRepeat < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-M")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a positive number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefRepeats = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || RefRepeats < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	RefRepeats2 = 1;
	if(argc >= 3 && argv[2][0] != '-'){
	  RefRepeats2 = strtol(argv[2],&qt,10);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || RefRepeats2 <= 0){
	    printf("%s 2nd value is invalid (must be >= 1):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxCov")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxCov = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || MaxCov <= 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxCovWt")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxCovWt = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || MaxCovWt < 0 || MaxCovWt > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'N':
      if(!strcmp(argv[0],"-NFSdelay")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	NFSdelay = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value is invalid (must be non-negative number):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-NoBpp")){
	NoBpp = 1;
	argv++;
	argc--;
	continue;
      }
    case 'O':
      if(!strcmp(argv[0],"-OverlapFilterSameRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative integer:%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	OverlapFilterSameRef = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || OverlapFilterSameRef > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-OverlapFilter")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by 2 non-negative numbers:%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	OverlapFilterQryRatio = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || OverlapFilterQryRatio > 1.0){
	  printf("%s 1st value is invalid (must be non-negative number <= 1.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);  exit(1);
	}
	OverlapFilterConfDelta = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative number):%s\n",argv[0],argv[2]);
	  fflush(stdout);  exit(1);
	}
	OverlapFilterQryLabRatio = strtod(argv[3],&qt);
	if(qt==argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd value is invalid (must be non-negative number <= 1.0):%s\n",argv[0],argv[3]);
	  fflush(stdout);  exit(1);
	}

	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-OutlierQuality")){
	OutlierQuality = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  OutlierQuality = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt== 0 || isspace(*qt)) || OutlierQuality < 0 || OutlierQuality > 2){
	    printf("%s value is invalid (must be 0,1 or 2):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      break;
    case 'P':
      if(!strcmp(argv[0],"-PSFWidth")){
	MapPSFWidth = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-PairwiseOutlierFix")){
	if(argc < 2 || argv[1][0] == '-'){
	  PairwiseOutlierFix = 1;
	  argv++;
	  argc--;
	} else {
	  PairwiseOutlierFix = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PairwiseOutlierFix < 0 || PairwiseOutlierFix > 1){
	    printf("%s value of %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv += 2;
	  argc -= 2;
	}
	continue;
      }
      if(!strcmp(argv[0],"-PVendoutlier")){
	EndOutlierPV = 2;
	if(argc >= 2 && argv[1][0] != '-'){
	  EndOutlierPV = strtol(argv[1],&qt,10);
	  if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || EndOutlierPV > 2){
	    printf("%s value has invalid syntax (must be 0, 1 or 2):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-PVres")){
	PVres = 1;
	if(argc >= 2 && !(argv[1][0]=='-' && !isdigit(argv[1][1]))){
	  PVres = strtol(argv[1],&qt,10);
	  if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || PVres < 0){
	    printf("%s value has invalid syntax for non-negative int:%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-Palindromic")){
	double PVratio = 0.01;
	if(argc >= 2 && !(argv[1][0]=='-' && !isdigit(argv[1][1]))){
	  PVratio = strtod(argv[1],&qt);
	  if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || PVratio < 0.0){
	    printf("%s value has invalid syntax for non-negative number:%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	Palindromic_Delta = -log(PVratio)/log(10.0);
	argv++;
	argc--;
	continue;
      }
      break;
    case 'Q':
      if(!strcmp(argv[0],"-QPmem")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	QP_MIN_MEM = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || QP_MIN_MEM < 0){
	  printf("%s %s : value should be 0 or greater\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-QXfilter")){
	QXfilter = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-QXerror")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	QXerror = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || QXerror < 0 || QXerror > 1){
	  printf("%s value has invalid syntax (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'R':
      if(!strcmp(argv[0],"-RefineSizeRange")){
	if(argc < 5 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-' || argv[4][0] == '-'){
	  printf("%s flag must be following by 4 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	RefineMinRange1 = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineMinRange1 <= 0.0){
	  printf("%s 1st value invalid : should be positive ratio: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	RefineMaxRange1 = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || RefineMaxRange1 < RefineMinRange1){
	  printf("%s 2nd value invalid : should be ratio greater than 1st value: %s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	RefineMinRange2 = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || RefineMinRange2 <= 0.0){
	  printf("%s 3rd value invalid : should be positive ratio: %s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	RefineMaxRange2 = strtod(argv[4],&qt);
	if(qt == argv[4] || !(*qt==0 || isspace(*qt)) || RefineMaxRange2 < RefineMinRange2){
	  printf("%s 4th value invalid : should be ratio greater than 3rd value: %s\n",argv[0],argv[4]);
	  fflush(stdout);exit(1);
	}

	argv += 5;
	argc -= 5;
	continue;
      }
      if(!strcmp(argv[0],"-RAscoremem")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RAscoremem = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RAscoremem > 1){
	  printf("%s %s : value should be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RAmem")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RA_MIN_MEM = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RA_MIN_MEM < 0){
	  printf("%s %s : value should be 0 or greater\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  RA_MIN_TIME = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value should be 0 or greater:%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 4 && argv[3][0] != '-'){
	    RA_MAX_MEM = strtod(argv[3],&qt);
	    if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || RA_MAX_MEM > 1.0){
	      printf("%s 3rd value invalid (should be fraction between 0 and 1):%s\n",argv[0],argv[3]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefineHmap")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefineHmap = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineHmap < 0 || RefineHmap > 1){
	  printf("%s %s : value should be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefineXRange")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0, 1 or 2\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefineXRange = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineXRange > 2){
	  printf("%s value is invalid (must be 0 or 1 or 2):%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefineSwitchFade")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefineSwitchFade = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineSwitchFade > 2){
	  printf("%s value is invalid (must be 0,1 or 2):%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefineFillHmap")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefineFillHmap = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineFillHmap > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefineFillMap")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefineFillMap = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineFillMap > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefineAlignExt")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefineAlignExt = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineAlignExt > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefineOverwrite")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefineOverwrite = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineOverwrite > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefineFix")){
	RefineFix = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  RefineFix = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineFix < 0 || RefineFix > 2){
	    printf("%s %s : 2nd field must be an integer from 0 to 2\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-RefineErr")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s must be followed by 3 non-negative error rates\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefineEps1 = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineEps1 < 0.0){
	  printf("%s 1st value is invalid (must be non-negative error rate):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	RefineEps2 = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || RefineEps2 < 0.0){
	  printf("%s 2nd value is invalid (must be non-negative error rate):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	RefineEps3 = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || RefineEps3 < 0.0){
	  printf("%s 3rd value is invalid (must be non-negative error rate):%s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-RefineRange")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by two positive integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefineRange = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineRange <= 0){
	  printf("%s 1st value is invalid (must be positive integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	RangeSwitch = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || RangeSwitch <= 0){
	  printf("%s 2nd value is invalid (must be positive integer):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-RefineStep")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by two positive numbers\n",argv[0]);
	  exit(1);
	}
	RefineStep1 = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefineStep1 <= 0.0){
	  printf("%s 1st value is invalid (must be positive number):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	RefineStep2 = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || RefineStep2 <= 0.0){
	  printf("%s 2nd value is invalid (must be positive number):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	RefineStep3 = 0.10;
	RefineStep4 = 0.20;
	RefineStepCnt = 3;
	RefineStep5 = 0.20;
	RefineStep6 = 0.20;
	RefineStep7 = 0.10;

	if(argc >= 6 && argv[3][0] != '-' && argv[4][0] != '-' && argv[5][0] != '-'){
	  RefineStep3 = strtod(argv[3],&qt);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || RefineStep3 <= 0.0){
	    printf("%s 3rd value is invalid (must be positive number):%s\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }

	  RefineStep4 = strtod(argv[4],&qt);
	  if(qt == argv[4] || !(*qt==0 || isspace(*qt)) || RefineStep4 <= 0.0){
	    printf("%s 4th value is invalid (must be positive number):%s\n",argv[0],argv[4]);
	    fflush(stdout);exit(1);
	  }

	  RefineStepCnt = strtol(argv[5],&qt,10);
	  if(qt == argv[5] || !(*qt==0 || isspace(*qt)) || RefineStepCnt < 0){
	    printf("%s 5th value is invalid (must be non-negative integer):%s\n",argv[0],argv[5]);
	    fflush(stdout);exit(1);
	  }
	  
	  if(argc >= 8 && argv[6][0] != '-' && argv[7][0] != '-'){
	    RefineStep5 = strtod(argv[6],&qt);
	    if(qt == argv[6] || !(*qt==0 || isspace(*qt)) || RefineStep5 <= 0.0){
	      printf("%s 6th value is invalid (must be positive number):%s\n",argv[0],argv[6]);
	      fflush(stdout);exit(1);
	    }

	    RefineStep6 = strtod(argv[7],&qt);
	    if(qt == argv[7] || !(*qt==0 || isspace(*qt)) || RefineStep6 <= 0.0){
	      printf("%s 7th value is invalid (must be positive number):%s\n",argv[0],argv[7]);
	      fflush(stdout);exit(1);
	    }

	    if(argc >= 9 && argv[8][0] != '-'){
	      RefineStep7 = strtod(argv[8],&qt);
	      if(qt == argv[8] || !(*qt==0 || isspace(*qt)) || RefineStep7 <= 0.0){
		printf("%s 8th value is invalid (must be positive number):%s\n",argv[0],argv[8]);
		fflush(stdout);exit(1);
	      }
	      argv++;
	      argc--;
	    }

	    argv += 2;
	    argc -= 2;
	  }
	  argv += 3;
	  argc -= 3;
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-RTHETA_FIX")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a number (0 or 1)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RTHETA_FIX = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RTHETA_FIX < 0 || RTHETA_FIX > 1){
	  printf("%s value %s is invalid : must be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefineThreads")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	RefineThreads = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || RefineThreads <= 0){
	  printf("%s value %s is invalid : must be +ve integer\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefSplitStitch")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0,1 or 2\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	RefSplitStitch = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || RefSplitStitch > 2){
	  printf("%s value %s is invalid : must be 0,1 or 2\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RefSplit")){
	RefSplit = 1;
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by 2 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	Psplit = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || Psplit < 0.0){
	  printf("%s 1st value is invalid (must be non-negative fraction):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	RefSplitMinLabels = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt == 0 || isspace(*qt)) || RefSplitMinLabels < 0){
	  printf("%s 2nd value is invalid (must be non-negative integer):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 4 && argv[3][0] != '-'){
	  Rsplit = strtod(argv[3],&qt);
	  if(qt == argv[3] || !(*qt == 0 || isspace(*qt)) || Rsplit < 0.0){
	    printf("%s 3rd value is invalid (must be non-negative fraction):%s\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	if(Psplit <= 0.0)
	  RefSplit = 0;

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-ReplaceCov")){
	REPLACE_COV = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-RepeatRec")){
	RepeatRec = 1;
	if(argc < 3 || argv[1][0] == '-'){
	  argv++;
	  argc--;
	  continue;
	}
	RepeatLogPvalueRatio = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RepeatLogPvalueRatio < 0.0 || RepeatLogPvalueRatio > 1.0){
	  printf("%s value %s is invalid (LogPvalueRatio must be a number between 0.0 and 1.0)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	RepeatKbShiftRatio = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid size in kb\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 4 && argv[3][0] != '-'){
	  RepeatKbShiftRatio2 = strtod(argv[3],&qt);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value %s is invalid size in kb\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 5 && argv[4][0] != '-'){
	    RepeatMaxLabelDrop = strtol(argv[4],&qt,10);
	    if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th value %s is invalid integer\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);	      
	    }
	    if(argc >= 6 && argv[5][0] != '-'){
	      RepeatLogPvalueRatioPerLabel = strtod(argv[5],&qt);
	      if(qt == argv[5] || !(*qt==0 || isspace(*qt)) || RepeatLogPvalueRatioPerLabel < 0.0 || RepeatLogPvalueRatioPerLabel > 1.0){
		printf("%s value 6th %s is invalid (LogPvalueRatioPerLabel must be a number between 0.0 and 1.0)\n",argv[0],argv[5]);
		fflush(stdout);exit(1);
	      }
	      
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-Repeat")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number (KB)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RepeatShift = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RepeatShift < 0.0){
	  printf("%s value %s is invalid (minlenKB must be a non-negative)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-RepeatMask")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer (MaxShift)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RepeatMaxShift = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RepeatMaxShift < 0){
	  printf("%s %s is invalid (MaxShift must be a non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  RepeatPvalueRatio = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || RepeatPvalueRatio > 1.0){
	    printf("%s %s %s is invalid (PvalueRatio must be <= 1.0)\n",argv[0],argv[1],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'S':
      if(!strcmp(argv[0],"-ScanRange")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ScanRange = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid : should be non-negative number: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-SplitSegDup_OneSided")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by 2 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	SplitSegDup_OneSided = strtol(argv[1],&qt,0);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || SplitSegDup_OneSided > 2){
	  printf("%s 1st arg is invalid : should be 0, 1 or 2: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	SplitSegDup_OneSided_MaxLen = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd arg is invalid : should be length in kb: %s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }

      if(!strcmp(argv[0],"-SegDupMask")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative 63-bit mask value\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	SegDupMask = strtoll(argv[1],&qt,0);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st arg is invalid : should be non-negative 63-bit mask value: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-SiteMerge")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	SiteMerge = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || SiteMerge > 1){
	  printf("%s arg is invalid : must be 0 or 1: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-SplitSegDup")){
	if(argc < 8 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-' || argv[4][0] == '-' || argv[5][0] == '-' || argv[6][0] == '-' || argv[7][0] == '-'){
	  printf("%s flag must be followed by 7 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	SplitSegDup = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st arg is invalid : should be non-negative number: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	SplitSegDupN = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd arg is invalid : should be non-negative number: %s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	SplitSegDupE = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd arg is invalid : should be non-negative number: %s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	SplitSegDupEN = strtol(argv[4],&qt,10);
	if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	  printf("%s 4th arg is invalid : should be non-negative number: %s\n",argv[0],argv[4]);
	  fflush(stdout);exit(1);
	}
	
	SplitSegDupC = strtod(argv[5],&qt);
	if(qt == argv[5] || !(*qt==0 || isspace(*qt))){
	  printf("%s 5th arg is invalid : should be non-negative number: %s\n",argv[0],argv[5]);
	  fflush(stdout);exit(1);
	}

	SplitSegDupLM = strtod(argv[6],&qt);
	if(qt == argv[6] || !(*qt==0 || isspace(*qt))){
	  printf("%s 6th arg is invalid : should be non-negative number: %s\n",argv[0],argv[6]);
	  fflush(stdout);exit(1);
	}
	
	SplitSegDupNM = strtol(argv[7],&qt,10);
	if(qt == argv[7] || !(*qt==0 || isspace(*qt))){
	  printf("%s 7th arg is invalid : should be non-negative number: %s\n",argv[0],argv[7]);
	  fflush(stdout);exit(1);
	}

	SplitSegDupMaxOutlier = -1.0;

	if(argc >= 9 && argv[8][0] != '-'){
	  SplitSegDupMaxOutlier = strtod(argv[8],&qt);
	  if(qt == argv[8] || !(*qt==0 || isspace(*qt))){
	    printf("%s 8th arg is invalid : should be non-negative number: %s\n",argv[0],argv[8]);
	    fflush(stdout);exit(1);
	  }

	  argv++;
	  argc--;
	}

	argv += 8;
	argc -= 8;
	continue;
      }
      if(!strcmp(argv[0],"-SmallDupFilter")){
	if(argc < 5 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-' || argv[4][0] == '-'){
	  printf("%s flag must be followed by 4 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	
	SmallDupMaxOverlap = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st arg is invalid : should be non-negative number: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	SmallDupMinOverlap = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd arg is invalid : should be non-negative number: %s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	SmallDupMaxMiss = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd arg is invalid : should be non-negative number: %s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	double Pvalue = strtod(argv[4],&qt);
	if(qt == argv[4] || !(*qt==0 || isspace(*qt)) || Pvalue <= 0.0 || Pvalue > 1.0){
	  printf("%s 4th arg is invalid : should be positive number <= 1.0: %s\n",argv[0],argv[4]);
	  fflush(stdout);exit(1);
	}
	SmallDupConfBC = -log(Pvalue)/log(10.0);

	argv += 5;
	argc -= 5;
	continue;
      }
      if(!strcmp(argv[0],"-SinglelineXmap")){
	SinglelineXmap = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-SplitHmap")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	SplitHmap = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || SplitHmap < 0 || SplitHmap > 1){
	  printf("%s %s : value should be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-SeparateQueries")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PairsplitSeparateQueries = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PairsplitSeparateQueries < 0 || PairsplitSeparateQueries > 1){
	  printf("%s %s : value should be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-ScanFilter")){
	if(argc < 5 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-' || argv[4][0] == '-'){
	  printf("%s flag must be followed by 4 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	ScanFilter = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st arg %s is not a non-negative integer\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	ScanFilterRm = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s : 2nd arg %s is not a valid non-negative number\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	ScanFilterRb = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s : 3rd arg %s is not a valid non-negative number\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	ScanFilterRe = strtod(argv[4],&qt);
	if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	  printf("%s : 4th arg %s is not a valid non-negative number\n",argv[0],argv[4]);
	  fflush(stdout);exit(1);
	}

	argv += 5;
	argc -= 5;
	continue;
      }
      if(!strcmp(argv[0],"-ScanScaling")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1 to 4 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ScanCorrection = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st arg %s is not a non-negative integer\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  ScanCorrectMinSNR = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s : 2nd arg %s is not a valid non-negative minSNR value\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 4 && argv[3][0] != '-'){
	    ScanCorrectMinLen = strtod(argv[3],&qt);
	    if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	      printf("%s : 3rd arg %s is not a valid non-negative kb value\n",argv[0],argv[3]);
	      fflush(stdout);exit(1);
	    }
	    if(argc >= 5 && argv[4][0] != '-'){
	      ScanCorrectMinMaps = strtol(argv[4],&qt,10);
	      if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
		printf("%s : 4th arg %s is not a valid non-negative integer\n",argv[0],argv[4]);
		fflush(stdout);exit(1);
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-Site_Pen")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	Site_Pen = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  Site_Pen_Ext = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value %s is invalid\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-SNRpvalue")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	BackgroundSNR = strdup(argv[1]);
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-SNRbackground")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by filename and integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	BackgroundSNR = strdup(argv[1]);
	bgSNRcnt[0] = strtol(argv[2],&qt,10);	
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || bgSNRcnt[0] <= 0){
	  printf("%s 2nd value of %s is invalid (must be > 0)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-SNR")){
	MapSNR = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-ScaleDelta")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0]=='-'){
	  printf("%s must be followed by two +ve numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ScaleDeltaSize = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || ScaleDeltaSize <= 0.0){
	  printf("%s 1st value of %s is invalid (must be > 0.0)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	ScaleDelta = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || ScaleDelta < 0){
	  printf("%s 2nd value of %s is invalid (must be >= 0)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-ScaleDeltaBPP")){
	ScaleDeltaBPP = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-SecondBest")){
	SecondBest = 1;
	if(argc >= 2){
	  SecondSpacing = strtod(argv[1],&qt);
	  if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || SecondSpacing < 0.0){
	    printf("%s value of %s is invalid\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 3 && argv[2][0] != '-'){
	    double Pvalue = strtod(argv[2], &qt);
	    if(qt==argv[2] || !(*qt == 0 || isspace(*qt)) || Pvalue < 0.0 || Pvalue > 1.0){
	      printf("%s 2nd value of %s is invalid (must be between 0 and 1)\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }
	    LogPvThreshold2 = -log(Pvalue)/log(10.0);
	    if(argc >= 4 && argv[3][0] != '-'){
	      ScoreThreshold2 = strtod(argv[3],&qt);
	      if(qt==argv[3] || !(*qt == 0 || isspace(*qt))){
		printf("%s 3rd value of %s is invalid (must be >= 0)\n",argv[0],argv[3]);
		fflush(stdout);exit(1);
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-S")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(MultiMatches){
	  printf("WARNING:%s must be set before -MultiMatches\n",argv[0]);
	  fflush(stdout);
	}
	if(CutFlip){
	  printf("WARNING:%s must be set before -CutFlip\n",argv[0]);
	  fflush(stdout);
	}
	ScoreThreshold = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'T':
      if(!strcmp(argv[0],"-Truncate")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s must be followed by 3 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	Truncate = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid : must be a Length in kb\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	TruncateN = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt == 0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid : must be a non-negative integer number of Labels\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	TruncateMaxOutlier = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt == 0 || isspace(*qt))){
	  printf("%s 3rd value %s is invalid : must be a Length in kb\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	TruncateFix = 1;
	TruncateMaxOutlierM = 3;

	if(argc >= 5 && argv[4][0] != '-'){
	  TruncateFix = strtol(argv[4],&qt,10);
	  if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	    printf("%s 4th value %s is invalid : must be a Length in kb\n",argv[0],argv[4]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 6 && argv[5][0] != '-'){
	    TruncateMaxOutlierM = strtol(argv[5],&qt,10);
	    if(qt == argv[5] || !(*qt==0 || isspace(*qt))){
	      printf("%s 5th value %s is invalid : must be max number of misaligned labels for non-outlier intervals\n",argv[0],argv[5]);
	      fflush(stdout);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-TrimmedNoMerge")){
	TrimmedNoMerge = 9999.0;
	if(argc >= 2 && argv[1][0] != '-'){
	  TrimmedNoMerge = strtod(argv[1],&qt);
	  if(qt == argv[1] || !(*qt == 0 || isspace(*qt))){
	    printf("%s value %s is invalid : must be a length in kb\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 3 && argv[2][0] != '-'){
	    TrimmedNoMergeMaxEnd = strtod(argv[2],&qt);
	    if(qt == argv[2] || !(*qt == 0 || isspace(*qt))){
	      printf("%s value %s is invalid : must be a length in kb\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }
            if(argc >= 4 && argv[3][0] != '-'){
	      TrimmedNoMergeBoth = strtol(argv[3],&qt,0);
	      if(qt == argv[3] || !(*qt == 0 || isspace(*qt))){
		printf("%s 3rd value is invalid : must be 63-bit non-negative mask: %s\n", argv[0], argv[3]);
		fflush(stdout);exit(1);
	      }

	      TrimmedNoMergeSegDupMask = 0;
	      if(argc >= 5 && argv[4][0] != '-'){
		TrimmedNoMergeSegDupMask = strtol(argv[4],&qt,0);
		if(qt == argv[4] || !(*qt == 0 || isspace(*qt))){
		  printf("%s 4th value is invalid : must be 63-bit non-negative mask: %s\n", argv[0], argv[4]);
		  fflush(stdout);exit(1);
		}
		
		if(argc >= 6 && argv[5][0] != '-'){
		  TrimmedNoMergeSegDupFlank = strtol(argv[5],&qt,10);
		  if(qt == argv[5] || !(*qt == 0 || isspace(*qt))){
		    printf("%s 5th value is invalid : must be non-negative integer: %s\n", argv[0],argv[5]);
		    fflush(stdout);exit(1);
		  }

		  argv++;
		  argc--;
		}

		argv++;
		argc--;
	      }

	      argv++;
	      argc--;
            }

	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-Tracksites")){
	Tracksites = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-TotalThreads")){
	TotalThreads = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  TotalThreads = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt == 0 || isspace(*qt))){
	    printf("%s value %s is invalid : must be integer\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 3 && argv[2][0] != '-'){
	    TotalSlots = strtol(argv[2],&qt,10);
	    if(qt == argv[2] || !(*qt == 0 || isspace(*qt))){
	      printf("%s 2nd value %s is invalid : must be integer\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }

	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;

	continue;
      }
      if(!strcmp(argv[0],"-TrimmedNoMergeStrict")){
	TrimmedNoMergeStrict = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  TrimmedNoMergeStrict = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || TrimmedNoMergeStrict > 1){
	    printf("%s value %s is invalid : must be 0 or 1\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv += 2;
	  argc -= 2;
	} else {
	  argv++;
	  argc--;
	}

	continue;
      }
      if(!strcmp(argv[0],"-TrimNorm")){
	TrimNorm = 0.0;
	if(argc >= 2 && !(argv[1][0]=='-' && !isdigit(argv[1][1]))){
	  TrimNorm = strtod(argv[1],&qt);
	  if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	    printf("%s value has invalid syntax for non-negative number:%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 3 && argv[2][0] != '-'){
	    TrimNormLen = strtod(argv[2],&qt);
	    if(qt==argv[2] || !(*qt == 0 || isspace(*qt))){
	      printf("%s 2nd value has invalid syntax for non-negative number:%s\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-TrimNormMin")){
	TrimNormMin = 1.0;
	TrimNormMinCov = 0.0;
	TrimNormMinEnd = 1.0;
	TrimNormMinRel = 1.0;

	if(argc >= 2 && !(argv[1][0]=='-' && !isdigit(argv[1][1]))){
	  TrimNormMinEnd = TrimNormMin = strtod(argv[1],&qt);
	  if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || TrimNormMin < 0.0){
	    printf("%s 1st value has invalid syntax for non-negative number:%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 3 && !(argv[2][0]=='-' && !isdigit(argv[2][1]))){
	    TrimNormMinCov = strtod(argv[2],&qt);
	    if(qt==argv[2] || !(*qt == 0 || isspace(*qt)) || TrimNormMinCov < 0.0){
	      printf("%s 2nd value has invalid syntax for non-negative number:%s\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }
	    if(argc >= 4 && !(argv[3][0]=='-' && !isdigit(argv[3][1]))){
	      TrimNormMinEnd = strtod(argv[3],&qt);
	      if(qt==argv[3] || !(*qt == 0 || isspace(*qt)) || TrimNormMinEnd < 0.0){
		printf("%s 3rd value has invalid syntax for non-negative number:%s\n",argv[0],argv[3]);
		fflush(stdout);exit(1);
	      }
	      if(argc >= 5 && !(argv[4][0]=='-' && !isdigit(argv[4][1]))){
		TrimNormMinRel = strtod(argv[4],&qt);
		if(qt==argv[4] || !(*qt==0 || isspace(*qt)) || TrimNormMinRel < 0.0 || TrimNormMinRel > 1.0){
		  printf("%s 4th value is invalid : must be number between 0 and 1:%s\n",argv[0],argv[4]);
		  fflush(stdout);exit(1);
		}
		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-TrimOutlier")){
	TrimOutlier = 1;
	TrimOutlierKB = 0.0;
	TrimOutlierLabels = 1000;

	if(argc >= 2 && argv[1][0] != '-'){
	  TrimOutlier = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt== 0 || isspace(*qt)) || TrimOutlier > 1){
	    printf("%s 1st value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }

	  if(argc >= 3 && argv[2][0] != '-'){
	    TrimOutlierKB = strtod(argv[2],&qt);
	    if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	      printf("%s 2nd value is invalid (must be non-negative size in kb):%s\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }

	    if(argc >= 4 && argv[3][0] != '-'){
	      TrimOutlierLabels = strtol(argv[3],&qt,10);
	      if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
		printf("%s 3rd value is invalid (must be non-negative integer):%s\n",argv[0],argv[3]);
		fflush(stdout);exit(1);
	      }

	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-TrimmedNoMergeStrict")){
	TrimmedNoMergeStrict = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  TrimmedNoMergeStrict = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt== 0 || isspace(*qt))){
	    printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-TrimNormChim")){
	if(argc < 2){
	  printf("%s must be followed by 1-2 integers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	TrimNormChim = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value has invalid syntax for integer:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  TrimNormFP = strtol(argv[2],&qt,10);
	  if(qt==argv[2] || !(*qt == 0 || isspace(*qt))){
	    printf("%s 2nd value has invalid syntax for integer:%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-TrimNormMed")){
	if(argc < 2){
	  printf("%s must be followed by an integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	TrimNormMed = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value has invalid syntax integer:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-TrimNormEnd")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	TrimNormEnd = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || TrimNormEnd < 0){
	  printf("%s value has invalid syntax for non-negative int:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-TrimFactor")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	TrimFactor = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value has invalid syntax for float:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(TrimFactor < 0.0 || TrimFactor >= 0.50){
	  printf("%s value is invalid (must be above between 0 and 0.50):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-T")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(MultiMatches){
	  printf("WARNING:%s must be set before -MultiMatches\n",argv[0]);
	  fflush(stdout);
	}
	if(CutFlip){
	  printf("WARNING:%s must be set before -CutFlip\n",argv[0]);
	  fflush(stdout);
	}
	double Pvalue = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value has invalid syntax for float:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(Pvalue <= 0.0 || Pvalue > 1.0){
	  printf("%s value of %0.6e is invalid (must be above 0 and <= 1):%s\n",argv[0], Pvalue, argv[1]);
	  fflush(stdout);exit(1);
	}
	sett = true;
	LogPvThreshold = -log(Pvalue)/log(10.0);
	if(VERB>=2){
	  printf("LogPvThreshold = %0.2f\n",LogPvThreshold);
	  fflush(stdout);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-TE")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	double Pvalue = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value has invalid syntax for float:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(Pvalue <= 0.0 || Pvalue > 1.0){
	  printf("%s value of %0.6e is invalid (must be above 0 and <= 1):%s\n",argv[0],Pvalue,argv[1]);
	  fflush(stdout);exit(1);
	}
	LogPvThresholdTE = -log(Pvalue)/log(10.0);
	TEinit = 1;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-TS")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	double Pvalue = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value has invalid syntax for float:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(Pvalue <= 0.0 || Pvalue > 1.0){
	  printf("%s value of %0.6e is invalid (must be above 0 and <= 1):%s\n",argv[0],Pvalue,argv[1]);
	  fflush(stdout);exit(1);
	}
	LogPvThresholdTS = -log(Pvalue)/log(10.0);
	TSinit = 1;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-TBmult")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	TBmult = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || TBmult >= 1.0 || TBmult < 0.0){
	  printf("%s value has invalid syntax or value (must be between 0 and 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-TB")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by Refinement Pvalue sequence\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	TBcnt = 0;
        while(argc >=2+TBcnt && argv[1+TBcnt][0] != '-'){
	  if(TBcnt >= 16){
	    printf("-TB : only 16 values (or less) are supported: excess arg = %s\n", argv[1+TBcnt]);
	    fflush(stdout);exit(1);
	  }
	  double Pvalue = strtod(argv[1+TBcnt],&qt);
	  if(qt==argv[1+TBcnt] || !(*qt==0 || isspace(*qt)) || Pvalue <= 0.0 || Pvalue >= 1.0){
	    printf("%s : %d'th Pvalue is invalid (must be between 0 and 1):%s\n",argv[0],TBcnt+1,argv[1]);
	    fflush(stdout);exit(1);
	  }
	  LogPvThresholdTB[TBcnt++] = -log(Pvalue)/log(10.0);
	}
	TBinit = 1;
	if(VERB>=2){
	  printf("TBcnt=%d:TB=%0.2f .. %0.2f\n",TBcnt,LogPvThresholdTB[0],LogPvThresholdTB[TBcnt-1]);
	  fflush(stdout);
	}
	argv += 1+TBcnt;
	argc -= 1+TBcnt;
	continue;
      }
      break;
    case 'U':
      if(!strcmp(argv[0],"-UseEndoutliers")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	UseEndoutliers = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || UseEndoutliers < 0 || UseEndoutliers > 1){
	  printf("%s value is invalid (must be 0 or 1 ):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'V':
      if(!strcmp(argv[0],"-ViterbiSwitch")){
	ViterbiSwitch = 1;
	if(!(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-')){
	  printf("%s must be followed by two sf values in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ViterbiLabel = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value of %s is invalid (must be non-negative)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	ViterbiSize = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value of %s is invalid (must be non-negative)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	if(argc >= 4 && argv[3][0] != '-'){
	  ViterbiLabelIter = strtol(argv[3],&qt,10);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd vale of %s is invalid (must be non-negative integer)\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  
	  if(argc >= 5 && argv[4][0] != '-'){
	    ViterbiSizeIter = strtol(argv[4],&qt,10);
	    if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th vale of %s is invalid (must be non-negative integer)\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }
	    
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      break;
    case 'X':
      if(!strcmp(argv[0],"-XMAPmapWT")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	XMAPmapWT = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || XMAPmapWT < 0 || XMAPmapWT > 1){
	  printf("%s value is invalid (must be 0 or 1 ):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-XmapStatWrite")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	XmapStatWrite = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-XmapStatRead")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	XmapStatRead = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-XmapInterleaved")){
	if(argc<2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	XmapInterleaved = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || XmapInterleaved < 0 || XmapInterleaved > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'a':
      if(!strcmp(argv[0],"-align_format")){
	if(argc<2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	align_format = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || align_format < 0 || align_format > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-alignscore")){
	AlignScore = 1;
	argv++;
	argc--;
	continue;
      }
      break;
    case 'b':
      if(!strcmp(argv[0],"-breakChiSq")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by 2 to 16 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	BreakChiSq = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || BreakChiSq > 1.0){
	  printf("%s 1st value of %s is invalid (must be number between 0.0 and 1.0)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	BreakFrac = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt == 0 || isspace(*qt)) || BreakChiSq > 1.0){
	  printf("%s 2nd value of %s is invalid (must be number between 0.0 and 1.0)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	
	BreakMaxSize = 999.0;
	BreakMaxCov = 60.0;
	BreakMinSf = 0.300 * 1000.0;
	BreakMinSr = 0.02 * 1000.0;

	BreakMinInterval = 2.0;

	BreakSdRatio = 1.2;
	BreakMinCov = 3.0;
	BreakMinCovScale = 1.5;
	BreakMinCovRatio = 1.2;

	BreakFragile = 100.0;
	BreakSegDup = 100.0;
	BreakChimNorm = 0;

	if(argc >= 4 && argv[3][0] != '-'){
	  BreakMaxSize = strtod(argv[3],&qt);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value of %s is invalid (must be non-negative size in kb)\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 5 && argv[4][0] != '-'){
	    BreakMaxCov = strtod(argv[4],&qt);
	    if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th value of %s is invalid (must be non-negative number)\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }
	    if(argc >= 6 && argv[5][0] != '-'){
	      BreakMinSf = strtod(argv[5],&qt);
	      if(qt == argv[5] || !(*qt==0 || isspace(*qt))){
		printf("%s 5th value of %s is invalid (must be non-negative number)\n",argv[0],argv[5]);
		fflush(stdout);exit(1);
	      }
	      BreakMinSf *= 1000.0;/* Arg is in kb, MolSd is in bp */

	      if(argc >= 7 && argv[6][0] != '-'){
		BreakMinSr = strtod(argv[6],&qt);
		if(qt == argv[6] || !(*qt==0 || isspace(*qt))){
		  printf("%s 6th value of %s is invalid (must be non-negative number)\n",argv[0],argv[6]);
		  fflush(stdout);exit(1);
		}
		BreakMinSr *= 1000.0;/* Interval is in kb, MolSd is in bp */

		if(argc >= 8 && argv[7][0] != '-'){
		  BreakMinInterval = strtod(argv[7],&qt);
		  if(qt == argv[7] || !(*qt==0 || isspace(*qt))){
		    printf("%s 7th value of %s is invalid (must be non-negative number)\n",argv[0],argv[7]);
		    fflush(stdout);exit(1);
		  }

		  if(argc >= 9 && argv[8][0] != '-'){
		    BreakSdRatio = strtod(argv[8],&qt);
		    if(qt == argv[8] || !(*qt==0 || isspace(*qt))){
		      printf("%s 8th value of %s is invalid (must be non-negative number)\n",argv[0],argv[8]);
		      fflush(stdout);exit(1);
		    }

		    if(argc >= 10 && argv[9][0] != '-'){
		      BreakMinCov = strtod(argv[9],&qt);
		      if(qt == argv[9] || !(*qt==0 || isspace(*qt))){
			printf("%s 9th value of %s is invalid (must be non-negative number)\n",argv[0],argv[9]);
			fflush(stdout);exit(1);
		      }
		    
		      if(argc >= 11 && argv[10][0] != '-'){
			BreakMinCovScale = strtod(argv[10],&qt);
			if(qt == argv[10] || !(*qt==0 || isspace(*qt))){
			  printf("%s 10th value of %s is invalid (must be non-negative number)\n",argv[0],argv[10]);
			  fflush(stdout);exit(1);
			}
			
			if(argc >= 12 && argv[11][0] != '-'){
			  BreakMinCovRatio = strtod(argv[11],&qt);
			  if(qt == argv[11] || !(*qt==0 || isspace(*qt))){
			    printf("%s 11th value of %s is invalid (must be non-negative number)\n",argv[0],argv[11]);
			    fflush(stdout);exit(1);
			  }
			  
			  if(argc >= 13 && argv[12][0] != '-'){
			    NoBreak = strtol(argv[12],&qt,10);
			    if(qt == argv[12] || !(*qt==0 || isspace(*qt))){
			      printf("%s 12th value %s is invalid (must be non-negative integer)\n",argv[0],argv[12]);
			      fflush(stdout);exit(1);
			    }
			    if(argc >= 14 && argv[13][0] != '-'){
			      IndelMaxSize = strtod(argv[13],&qt);
			      if(qt == argv[13] || !(*qt==0 || isspace(*qt))){
				printf("%s 13th value %s is invalid (must be size in kb)\n",argv[0],argv[13]);
				fflush(stdout);exit(1);
			      }
			      if(argc >= 15 && argv[14][0] != '-'){
				BreakFragile = strtod(argv[14],&qt);
				if(qt == argv[14] || !(*qt==0 || isspace(*qt))){
				  printf("%s 14th value %s is invalid (must be non-negative number))\n",argv[0],argv[14]);
				  fflush(stdout);exit(1);
				}
				if(argc >= 16 && argv[15][0] != '-'){
				  BreakSegDup = strtod(argv[15],&qt);
				  if(qt == argv[15] || !(*qt==0 || isspace(*qt))){
				    printf("%s 15th value %s is invalid (must be non-negative number)\n",argv[0],argv[15]);
				    fflush(stdout);exit(1);
				  }
				  if(argc >= 17 && argv[16][0] != '-'){
				    BreakChimNorm = strtod(argv[16],&qt);
				    if(qt == argv[16] || !(*qt==0 || isspace(*qt))){
				      printf("%s 16th value %s is invalid (must be non-negative number)\n",argv[0],argv[16]);
				      fflush(stdout);exit(1);
				    }

				    argv++;
				    argc--;
				  }
				  argv++;
				  argc--;
				}
				argv++;
				argc--;
			      }
			      argv++;
			      argc--;
			    }
			    argv++;
			    argc--;
			  }
			  argv++;
			  argc--;
			}
			argv++;
			argc--;
		      }
		      argv++;
		      argc--;
		    }
		    argv++;
		    argc--;
		  }
		  argv++;
		  argc--;
		}
		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-boostThreads")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number (Time in seconds)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	BoostThreadsTime = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || BoostThreadsTime < 0){
	  printf("%s 1st value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	BoostThreads = (BoostThreadsTime > 0) ? 1 : 0;

	if(argc >= 3 && argv[2][0] != '-'){
	  BoostThreads = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt == 0 || isspace(*qt)) || BoostThreads <= 0){
	    printf("%s 2nd value %s is invalid : must be positive integer\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv += 3;
	  argc -= 3;
	} else {
	  argv += 2;
	  argc -= 2;
	}

	continue;
      }
      if(!strcmp(argv[0],"-break")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by at least 2 ints\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	BreakIDinc = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || BreakIDinc <= 0){
	  printf("%s value is invalid (must be +ve integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	/* count number of IDs (ints) in argv[2..argc-1] */
	BreakIDcnt = 0;
	for(int i = 2; i < argc; i++){
	  if(argv[i][0] == '-')
	    break;
	  if(strchr(argv[i],'.') == NULL)
	     BreakIDcnt++;
	}
	if(VERB>=2){
	  printf("Detected %d -break IDs on command line\n",BreakIDcnt);
	  fflush(stdout);
	}
	if(BreakID)
	  delete [] BreakID;
	BreakID = new long long[BreakIDcnt];
	if(BreakPointCnt)
	  delete [] BreakPointCnt;
	BreakPointCnt = new int[BreakIDcnt];
	if(BreakPoint)
	  delete [] BreakPoint;
	BreakPoint = new double*[BreakIDcnt];
	int b = 0, i = 2;
	for(; i < argc; ){
	  if(argv[i][0] == '-')
	    break;
	  if(strchr(argv[i],'.') != NULL){
	    printf("%s : %d'th arg (%s) must be integer (b=%d)\n",argv[0],i,argv[i],b);
	    fflush(stdout);exit(1);
	  }
	  if(DEBUG) assert(b < BreakIDcnt);
	  BreakID[b] = strtol(argv[i],&qt,10);
	  if(qt==argv[i] || !(*qt == 0 || isspace(*qt)) || BreakID[b] <= 0){
	    printf("%s %dth arg is invalid (must be +ve integer):%s\n",argv[0],i,argv[i]);
	    fflush(stdout);exit(1);
	  }
	  i++;

	  /* count number of floats till next int (or next option) */
	  BreakPointCnt[b] = 0;
	  for(int j = i; j < argc; j++, BreakPointCnt[b]++)
	    if(argv[j][0] == '-' || !strchr(argv[j],'.'))
	      break;
	  if(BreakPointCnt[b] <= 0){
	    printf("%s : ID[%d]=%lld is not followed by at least one float (breakpoint location):%s\n",argv[0],b,BreakID[b],argv[i]);
	    fflush(stdout);exit(1);
	  }
	  if(VERB>=2){
	    printf("BreakID[b=%d]=%lld, BreakPointCnt[%d]=%d, next argv[i=%d]=%s\n",b,BreakID[b],b,BreakPointCnt[b],i,argv[i]);
	    fflush(stdout);
	  }
	  BreakPoint[b] = new double[BreakPointCnt[b]];
	  int cnt = 0;
	  for(int j = i; j < i + BreakPointCnt[b]; j++, cnt++){
	    BreakPoint[b][cnt] = strtod(argv[j],&qt);
	    if(qt == argv[j] || !(*qt == 0 || isspace(*qt)) || BreakPoint[b][cnt] <= 0.0){
	      printf("%s : After ID[%d]=%lld breakpoint %s is invalid (must be +ve float)\n",argv[0],b,BreakID[b],argv[j]);
	      fflush(stdout);exit(1);
	    }
	    if(cnt > 0 && BreakPoint[b][cnt-1] >= BreakPoint[b][cnt]){
	      printf("%s : After ID[%d]=%lld breakpoints %s %s are invalid (must be in ascending order)\n",argv[0],b,BreakID[b],argv[j-1],argv[j]);
	      fflush(stdout);exit(1);
	    }
	  }

	  i += BreakPointCnt[b++]; 
	}
	if(DEBUG) assert(b == BreakIDcnt);
	argv += i;
	argc -= i;
	continue;
      }
      if(!strcmp(argv[0],"-bppadjust")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PostAdjust = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || PostAdjust < 0 || PostAdjust > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-bestalignments")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	BestAlignments = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || BestAlignments <= 0){
	  printf("%s value is invalid (must be between +ve integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-bnx")){
	CmapMergeBnx = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-bnxerror")){
	bnxerror = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-bnxmergeerror")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	bnxmergeerror = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || bnxmergeerror < 0){
	  printf("%s value is invalide (must be non-negative):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-bnxBigid")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	bnxBigid = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || bnxBigid < 0){
	  printf("%s value is invalide (must be non-negative):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-bnxversion")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a positive float\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	double Version = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || Version < 0.0){
	  printf("%s value is invalid (must be +ve):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	int majorVersion = (int)floor(Version);
	//	int minorVersion = ((int)floor(Version*10.0))%10;
	if(majorVersion==0)/* Version 0.x */
	  BNXVersion = 0;
	else if(majorVersion==1)/* Version 1.x */
	  BNXVersion = 1;
	else {
	  printf("%s value is not a valid BNX version (must be 0.x or 1.x):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(VERB>=2){
	  printf("BNXVersion=%d\n",BNXVersion);
	  fflush(stdout);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-biaswt")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative float\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	setbiaswt = true;
	biasWT = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || biasWT < 0.0 || biasWT > 1.0){
	  printf("%s value is invalid (must be between 0 and 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  PRbiasWT = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt == 0 || isspace(*qt)) || PRbiasWT < 0.0 || PRbiasWT > 1.0){
	    printf("%s value is invalid (must be between 0 and 1):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-biaswtEnd")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative float\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	biasWTend = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || biasWTend < 0.0 || biasWTend > 1.0){
	  printf("%s value is invalid (must be between 0 and 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-biaswtOutlier")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative float\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	biasWToutlier = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || biasWToutlier < 0.0 || biasWToutlier > 1.0){
	  printf("%s value is invalid (must be between 0 and 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-biaswtRefine")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative float\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	biasWTrefine = strtod(argv[1],&qt);

	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || biasWTrefine < 0.0 || biasWTrefine > 1.0){
	  printf("%s value is invalid (must be between 0 and 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	biasWTrefineType = 0;
	PRbiasWTrefine = 0.0;
	if(argc >= 3 && argv[2][0] != '-'){
	  biasWTrefineType = strtol(argv[2],&qt,10);
	  if(qt==argv[2] || !(*qt == 0 || isspace(*qt)) || biasWTrefineType < 0 || biasWTrefineType > 2){
	    printf("%s 2nd value is invalid (must be 0, 1 or 2):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 4 && argv[3][0] != '-'){
	    PRbiasWTrefine = strtod(argv[3],&qt);
	    if(qt==argv[3] || !(*qt == 0 || isspace(*qt)) || PRbiasWTrefine < 0.0 || PRbiasWTrefine > 1.0){
	      printf("%s value is invalid (must be between 0 and 1):%s\n",argv[0],argv[3]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	
	if(VERB>=2){
	  printf("biasWTrefine= %0.6f, biasWTrefineType= %d, PRbiasWTrefine= %0.6f\n",biasWTrefine,biasWTrefineType,PRbiasWTrefine);
	  fflush(stdout);
	}

	continue;
      }
      if(!strcmp(argv[0],"-bppScan")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s msut be followed by 3 positive numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	BppStepsiz = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || BppStepsiz <= 0.0){
	  printf("%s 1st value is invalid (must be a +ve size in kb):\n%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	BppSteps = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative integer)\n%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	BppMincnt = strtol(argv[3],&qt,10);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || BppMincnt <= 0){
	  printf("%s 3rd value is invalid (must be positive integer)\n%s\n", argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-bpp")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by resolution in pixels\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	bpp = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || bpp <= 0.0 ){
	  printf("%s value has invalid syntax or value:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	bpp *= 0.001;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-bed")){ 
	bed_filename = argv[1];
	//dobed = 1; //necessary?
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-biaswtDefer")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	biaswtDefer = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || biaswtDefer > 1){
	  printf("%s value %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }

      break;

    case 'c':
      if(!strcmp(argv[0],"-cresCheckDel")){
	cresFix = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  cresFix = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || cresFix > 1){
	    printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-cresRefine")){
	if(argc < 2 || (argv[1][0] == '-' && !isdigit(argv[1][1]))){
	  printf("%s flag must be followed by a number (0 or 1 or -1)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	cresRefine = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || cresRefine < -1 || cresRefine > 1){
	  printf("%s value %s is invalid : must be 0 or 1 or -1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-cresMaxLPdrop")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	cresMaxLPdrop = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid : must be float number\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-cresHapMaxLPdrop")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by two numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	cresHapL = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid : must be float number\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	cresHapMaxLPdrop = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid : must be float number\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-cres")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be floowed by 2 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	cres = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || cres < 0.0){
	  printf("%s %s %s : 1st value is invalid (must be positive number)\n",argv[0],argv[1],argv[2]);
	  fflush(stdout);exit(1);
	}

	cresN = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || cresN < 0){
	  printf("%s %s %s : 2nd value is invalid (must be positive integer)\n",argv[0],argv[1],argv[2]);
	  fflush(stdout);exit(1);
	}

	if(argc >= 4 && argv[3][0] != '-'){
	  cresF = strtod(argv[3],&qt);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || cresF < 0){
	    printf("%s %s %s %s :  3rd value is invalid (must be non-negative fraction)\n",argv[0],argv[1],argv[2],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-centeronly")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by distance in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CenterMapped = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || CenterMapped < 0.0){
	  printf("%s value has invalid syntax or value (must be non-negative number):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-contigsplit")){
	if(argc < 4 || argv[1][0]=='-'){
	  printf("%s must be followed by 2 fractions and filename (and an optional 3rd number)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ContigSplitRatio = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || ContigSplitRatio <= 0.0 || ContigSplitRatio >= 1.0){
	  printf("%s %s %s %s: R1 value of %s is invalid (must be between 0 and 1)\n",argv[0],argv[1],argv[2],argv[3],argv[1]);
	  fflush(stdout);exit(1);
	}
	ContigSplitRatio2 = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || ContigSplitRatio2 < 0.0 || ContigSplitRatio2 >= 1.0){
	  printf("%s %s %s %s : R2 value of %s is invalid (must be between 0 and 1)\n",argv[0],argv[1],argv[2],argv[3], argv[2]);
	  fflush(stdout);exit(1);
	}
	
	if(argv[3][0]=='-'){
	  printf("%s %s %s %s : CntFile value of %s is invalid (cannot begin with '-')\n",argv[0],argv[1],argv[2],argv[3],argv[3]);
	  fflush(stdout);exit(1);
	}

	ContigCntFile = argv[3];
	if(argc > 4 && argv[4][0] != '-'){
	  ContigSplitRatio3 = strtod(argv[4],&qt);
	  if(qt == argv[4] || !(*qt==0 || isspace(*qt)) || ContigSplitRatio3 <= 1.0){
	    printf("%s %s %s %s %s : R3 value of %s is invalid (must be above 1)\n", argv[0],argv[1],argv[2],argv[3],argv[4], argv[4]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}

	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-colors")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	int NColors = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || NColors <= 0){
	  printf("-colors value is invalid (must be >= 1):%s\n",argv[1]);
	  fflush(stdout);exit(1);
	}
	if(NColors > MAXCOLOR){
	  printf("-colors value %d is too large (MAXCOLOR=%d)\n",colors,MAXCOLOR);
	  fflush(stdout);exit(1);
	}
	if(DEBUG) assert(colors > 0);
	for(int c = colors; c < NColors; c++){
	  res[c] = res[0];
	  resSD[c] = resSD[0];
	  
	  FP[c] = FP[0];
	  FN[c] = FN[0];
	  SF[c] = SF[0];
	  SD[c] = SD[0];
	  SR[c] = SR[0];
	  SE[c] = SE[0];
	}
	colors = NColors;
	colors_set = 1;// NEW
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-contigQual")){
	noQuality = 0; //noQuality 0 means do quality
	argv++;
	argc--;
	continue;
      }
      break;
    case 'd':
      if(!strcmp(argv[0],"-defermerge")){
	DeferMerge = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-delta")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	DELTA = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || DELTA <= 0){
	  printf("%s value of %s is invalid (must be positive integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-deltaY")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	DELTA_Y = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || DELTA_Y <= 0 || DELTA_Y > DELTA_MAX){
	  printf("%s value of %s is invalid (must be positive integer <= %d)\n",argv[0],argv[1],DELTA_MAX);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  DELTA_Y2 = strtol(argv[2],&qt,10);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || DELTA_Y2 < 0 || DELTA_Y2 > DELTA_Y){
	    printf("%s value of %s is invalid (must be between 0 and %d)\n",argv[0],argv[2], DELTA_Y);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-deltaYRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	deltaYRef = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || deltaYRef < 0 || deltaYRef > DELTA_MAX){
	  printf("%s value of %s is invalid (must be non-negative integer <= %d)\n",argv[0],argv[1], DELTA_MAX);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-deltaExtYRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	deltaExtYRef = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || deltaExtYRef < 0 || deltaExtYRef > DELTA_MAX){
	  printf("%s value of %s is invalid (must be non-negative integer <= %d)\n",argv[0],argv[1], DELTA_MAX);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-deltaX")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	DELTA_X = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || DELTA_X <= 0 || DELTA_X > DELTA_MAX){
	  printf("%s value of %s is invalid (must be positive integer <= %d)\n",argv[0],argv[1], DELTA_MAX);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  DELTA_X2 = strtol(argv[2],&qt,10);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || DELTA_X2 < 0 || DELTA_X2 > DELTA_X){
	    printf("%s value of %s is invalid (must be between 0 and %d)\n",argv[0],argv[2], DELTA_X);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-deltaXRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	deltaXRef = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || deltaXRef < 0 || deltaXRef > DELTA_MAX){
	  printf("%s value of %s is invalid (must be non-negative integer <= %d)\n",argv[0],argv[1], DELTA_MAX);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-deltaExtXRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	deltaExtXRef = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || deltaExtXRef < 0 || deltaExtXRef > DELTA_MAX){
	  printf("%s value of %s is invalid (must be non-negative integer <= %d)\n",argv[0],argv[1], DELTA_MAX);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-deltaEndRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	deltaEndRef = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || deltaEndRef < 0 || deltaEndRef > DELTA_MAX){
	  printf("%s value of %s is invalid (must be non-negative integer <= %d)\n",argv[0],argv[1], DELTA_MAX);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'e':
      if(!strcmp(argv[0],"-extsplitOutlierC")){
	if(argc < 3 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1 to %d non-negative number(s)\n",argv[0],MAX_KB_BINS);
	  fflush(stdout); exit(1);
	}

	KCcnt = 0;

	extsplitOutlierC[KCcnt++] = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid : must be non-negative number\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	while(argc >= 2+KCcnt && argv[1+KCcnt][0] != '-'){
	  if(KCcnt >= MAX_KB_BINS){
	    printf("Too many args for %s option : Cannot exceed %d coverage values\n",argv[0],MAX_KB_BINS);
	    fflush(stdout);exit(1);
	  }

	  extsplitOutlierC[KCcnt] = strtod(argv[1+KCcnt],&qt);
	  if(qt == argv[1+KCcnt] || !(*qt==0 || isspace(*qt))){
	    printf("%s %d'th value %s is invalid : must be non-negative number\n",argv[0],1+KCcnt,argv[1+KCcnt]);
	    fflush(stdout);exit(1);
	  }
	  KCcnt++;
	}

	argv += 1+KCcnt;
	argc -= 1+KCcnt;

	continue;
      }
      if(!strcmp(argv[0],"-extsplitOutlier")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1 to %d non-negative number(s)\n",argv[0],4+MAX_KB_BINS);
	  fflush(stdout);exit(1);
	}
	extsplitOutlier = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid : must be non-negative size in kb\n",argv[0],argv[1]);
          fflush(stdout);exit(1);
	}
	extsplitOutlierLabels = 1000;
	extsplitOutlierLS = -1.0;
	extsplitOutlierAS = -1;
	KBcnt = 0;

	int shift = 0;
	if(argc >= 3 && argv[2][0] != '-' && !strchr(argv[2],'.')){
	  extsplitOutlierLabels = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value %s is invalid : must be non-negative integer (OR float with decimal point)\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  shift++;
	  argv++;
	  argc--;
	}

	if(argc >= 3 && argv[2][0] != '-'){
	  extsplitOutlierLS = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s %d'th value %s is invalid : must be non-negative size in kb\n",argv[-shift],shift+2,argv[2]);
	    fflush(stdout);exit(1);
	  }
	  
	  if(argc >= 4 && argv[3][0] != '-'){
	    extsplitOutlierAS = strtol(argv[3],&qt,10);
	    if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	      printf("%s %d'th value %s is invalid : must be non-negative integer\n",argv[-shift],shift+3,argv[3]);
	      fflush(stdout);exit(1);
	    }

	    KBcnt = 0;	    
	    extsplitOutlierKB[0] = extsplitOutlier;

	    while(argc >= 5+KBcnt && argv[4+KBcnt][0] != '-'){
	      if(KBcnt >= MAX_KB_BINS){
		printf("Too many args for %s option : Cannot exceed %d KB values (%d total args)\n",argv[-shift],MAX_KB_BINS,shift+3+MAX_KB_BINS);
		fflush(stdout);exit(1);
	      }
	      extsplitOutlierKB[KBcnt+1] = strtod(argv[4+KBcnt],&qt);
	      if(qt == argv[4+KBcnt] || !(*qt==0 || isspace(*qt))){
		printf("%s %d'th value %s is invalid : must be non-negative size in kb\n",argv[-shift],shift+4+KBcnt,argv[4+KBcnt]);
		fflush(stdout);exit(1);
	      }
	      if(extsplitOutlierKB[KBcnt+1] <= extsplitOutlier || extsplitOutlierKB[KBcnt+1] <= extsplitOutlierKB[KBcnt]){
		printf("%s %d'th value %s is too small : must exceed previous value %0.3f kb\n",argv[-shift],shift+4+KBcnt,argv[4+KBcnt], extsplitOutlierKB[KBcnt]);
		fflush(stdout);exit(1);
	      }

	      KBcnt++;
	    }

	    if(KBcnt == 1){
	      printf("%s must have at least 2 additional outlier sizes L1 .. Ln\n",argv[-shift]);
	      fflush(stdout);exit(1);
	    }
	    extsplitOutlierKB[KBcnt+1] = 9999.0;

	    argv += 1+KBcnt;
	    argc -= 1+KBcnt;
	  }
	  argv++;
	  argc--;
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-extsplitNE")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1 or 2 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	extsplitNEL = 0.0;
	extsplitNEE = 0.0;
	extsplitNE = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || extsplitNE > 1){
	  printf("%s 1st value %s is invalid : must be 0 or 1\n",argv[0],argv[1]);
          fflush(stdout);exit(1);
	}

	if(argc >= 3 && argv[2][0] != '-'){
	  extsplitNEL = strtod(argv[2],&qt);// WAS98 strtod(argv[1],&qt)
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value %s is invalid : must be size in kb\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  if(argc >= 4 && argv[3][0] != '-'){
	    extsplitNEE = strtod(argv[3],&qt);
	    if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	      printf("%s 3rd value %s is invalid : must be size in kb\n",argv[0],argv[3]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }

	  argv++;
	  argc--;
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-extsplitSpacingFix")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0, 1 or 2\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ExtsplitSpacingFix = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || ExtsplitSpacingFix > 2){
	  printf("%s 2nd value %s is invalid : must be 0,1 or 2\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-extsplitVerifyPeak")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ExtsplitVerifyPeak = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || ExtsplitVerifyPeak > 1){
	  printf("%s 2nd value %s is invalid : must be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-extTrim")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	extTrim = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || extTrim > 1){
	  printf("%s 2nd value %s is invalid : must be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-extsplitMinWT")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number <= 1.0\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	extsplitMinWT = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || extsplitMinWT > 1.0){
	  printf("%s 2nd value %s is invalid : must be non-negative number <= 1.0\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-endoutlierSwitch")){
	PoutlierEndSwitch = 1;
	if(!(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-')){
	  printf("%s must be followed by two numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PoutlierEndLabel = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PoutlierEndLabel <= 0.0){
	  printf("%s 1st value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	PoutlierEndSize = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || PoutlierEndSize <= 0.0){
	  printf("%s 2nd value of %s is invalid (must be positive)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-extonly")){
	extendonly = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  EndLen2 = EndLen = strtod(argv[1],&qt);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || EndLen < 0.0){
	    printf("%s 1st value of %s is invalid (Cannot be -ve)\n",argv[0], argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(EndLen <= 0.0)
	    extendonly = 0;
	  if(argc >= 3 && argv[2][0] != '-'){
	    EndLen2 = strtod(argv[2],&qt);
	    if(qt == argv[2] || !(*qt == 0 || isspace(*qt)) || EndLen2 < 0.0){
	      printf("%s 2nd value of %s is invalid (Cannot be -ve)\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }
	    EndLen2 = min(EndLen, EndLen2);
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-extsplit")){
	if(argc < 2 || argv[1][0] == '-'
	   || (argc > 2 && argv[2][0] != '-' && (argc < 8 || argv[3][0] == '-' || argv[4][0] == '-' || argv[5][0] == '-' || argv[6][0] == '-' || argv[7][0] == '-' || argv[8][0] == '-'))){
	  printf("%s must be followed by 0 OR 7-10 numbers and 1 filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	extendSplit = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid (non-negative number expected)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	if(argc > 2 && argv[2][0] != '-'){

	  extendFrac = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || extendFrac > 1.0){
	    printf("%s 2nd value %s is invalid (fractional number expected)\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  extendSplitFN = strtol(argv[3],&qt,10);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value %s is invalid (integer expected)\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }

	  extendSplitFNlen = strtod(argv[4],&qt);
	  if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	    printf("%s 4th value %s is invalid (size in kb expected)\n",argv[0],argv[4]);
	    fflush(stdout);exit(1);
	  }

	  extendSplitFP = strtol(argv[5],&qt,10);
	  if(qt == argv[5] || !(*qt==0 || isspace(*qt))){
	    printf("%s 5th value %s is invalid (integer expected)\n",argv[0],argv[5]);
	    fflush(stdout);exit(1);
	  }

	  extendSplitLen = strtod(argv[6],&qt);
	  if(qt == argv[6] || !(*qt==0 || isspace(*qt))){
	    printf("%s 6th value %s is invalid (integer expected)\n",argv[0],argv[6]);
	    fflush(stdout);exit(1);
	  }

	  extendSplitSpacing = strtol(argv[7],&qt,10);
	  if(qt == argv[7] || !(*qt==0 || isspace(*qt))){
	    printf("%s 7th value %s is invalid (integer expected)\n",argv[0],argv[7]);
	    fflush(stdout);exit(1);
	  }
		
	  extendSplitEndSpacingX = extendSplitSpacingX = max(extendSplitSpacing, max(extendSplitFP, extendSplitFN));
	  extendSplitEndSpacing = extendSplitSpacing;

	  int cnt = 0;/* number of optional args found */

	  if(argc >= 9 && argv[8][0] != '-' && isdigit(argv[8][0])){
	    cnt++;
	    extendSplitSpacingX = strtod(argv[8],&qt);
	    if(qt == argv[8] || !(*qt==0 || isspace(*qt))){
	      printf("%s 8th value %s is invalid (number or filename expected)\n",argv[0],argv[8]);
	      fflush(stdout);exit(1);
	    }
	    extendSplitEndSpacingX = extendSplitSpacingX;

	    if(argc >= 10 && argv[9][0] != '-' && isdigit(argv[9][0])){
	      cnt++;
	      extendSplitEndSpacing = strtol(argv[9],&qt,10);
	      if(qt == argv[9] || !(*qt==0 || isspace(*qt))){
		printf("%s 9th value %s is invalid (number or filename expected)\n",argv[0],argv[9]);
		fflush(stdout);exit(1);
	      }

	      if(argc >= 11 && argv[10][0] != '-' && isdigit(argv[10][0])){
		cnt++;
		extendSplitEndSpacingX = strtod(argv[10],&qt);
		if(qt == argv[10] || !(*qt==0 || isspace(*qt))){
		  printf("%s 10th value %s is invalid (number or filename expected)\n",argv[0],argv[10]);
		  fflush(stdout);exit(1);
		}
	      }
	    }
	  }

	  if(argc < 9+cnt){
	    printf("%s missing %d'th value (filename)\n",argv[0], 8+cnt);
	    fflush(stdout);exit(1);
	  }
	  if(argv[8+cnt][0] == '-' || isdigit(argv[8+cnt][0])){
	    printf("%s %d'th value %s is not a valid filename (cannot start with an integer or minus sign)\n",argv[0], 8+cnt, argv[8+cnt]);
	    fflush(stdout);exit(1);
	  }

	  ContigCntFile = argv[8+cnt];

	  argv += 7+cnt;
	  argc -= 7+cnt;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-everb")){
	EVERB = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-endoutlierFinal")){
	if(argc<2){
	  printf("%s flag must be followed by Pvalue\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PoutlierEndFinal = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PoutlierEndFinal < 0.0 || PoutlierEndFinal > 1.0){
	  printf("%s value  of %s is invalid (must be between 0 and 1)\n",argv[0], argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-endoutlier")){
	if(argc<2){
	  printf("%s flag must be followed by Pvalue\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PoutlierEnd = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PoutlierEnd < 0.0 || PoutlierEnd > 1.0){
	  printf("%s value  of %s is invalid (must be between 0 and 1)\n",argv[0], argv[1]);
	  fflush(stdout);exit(1);
	}
	PoutlierEnd_init = 1;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-endoutlierRef")){
	if(PoutlierEndCnt){
	  printf("WARNING : -endoutlierRef specified twice : ignoring first occurance with %d values\n", PoutlierEndCnt);
	  fflush(stdout);
	  PoutlierEndCnt = 0;
	}
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by Refinement end outlier probability value(s)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	while(argc >= 2+PoutlierEndCnt && argv[1+PoutlierEndCnt][0] != '-'){
	  if(PoutlierEndCnt >= 16){
	    printf("-endoutlierRef : only 16 values (or less) is supported: excess arg = %s\n", argv[1+PoutlierEndCnt]);
	    fflush(stdout);exit(1);
	  }
	  double Pvalue = strtod(argv[1+PoutlierEndCnt],&qt);
	  if(qt == argv[1+PoutlierEndCnt] || !(*qt==0 || isspace(*qt)) || Pvalue < (/*PoutlierEndCnt > 0 ? PoutlierEndRefine[PoutlierEndCnt-1] : */ 0.0) || Pvalue > 1.0){
	    //	    printf("%s : %d'th value  of %s is invalid (must be between 0 and 1 and in ascending order)\n",argv[0], PoutlierEndCnt+1, argv[1+PoutlierEndCnt]);
	    printf("%s : %d'th value  of %s is invalid (must be between 0 and 1)\n",argv[0], PoutlierEndCnt+1, argv[1+PoutlierEndCnt]);
	    fflush(stdout);exit(1);
	  }
	  PoutlierEndRefine[PoutlierEndCnt++] = Pvalue;
	}
	if(VERB>=2){
	  for(int t = 0; t < PoutlierEndCnt; t++)
	    printf("PoutlierEndRefine[%d]= %0.6e\n",t,PoutlierEndRefine[t]);
	  fflush(stdout);
	}
	argv += 1+PoutlierEndCnt;
	argc -= 1+PoutlierEndCnt;
	continue;
      }
      if(!strcmp(argv[0],"-extend")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by integer between 0 and 2\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	setextend = true;
	extend = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || extend < 0 || extend > 2){
	  printf("%s value of %s is invalid (must be between 0 and 2)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-extendWT")){
	if(argc < 5 || argv[1][0]=='-' || argv[2][0]=='-' || argv[3][0]=='-' || argv[4][0]=='-'){
	  printf("%s flag must be followed by 4 to 10 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	extendLen = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || extendLen < 0.0){
	  printf("%s 1st value of %s is invalid (must be non-negative length in kb\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	extendWT = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || extendWT < 0.0 || extendWT > 1.0){
	  printf("%s 2nd value of %s is invalid (must be WT between 0.0 and 1.0\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	extendFN = strtol(argv[3],&qt,10);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd value %s is invalid (integer expected)\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	extendFNlen = strtod(argv[4],&qt);
	if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	  printf("%s 4th value %s is invalid (size in kb expected)\n",argv[0],argv[4]);
	  fflush(stdout);exit(1);
	}

	extendMinWT = 0.0;
	extendMaxOutlierKBsplit = extendMaxOutlierKB = 1000.0;
	extendMaxOutlierLabelsSplit = extendMaxOutlierLabels = 1000;
	extendMaxWT = 1.0;

	if(argc >= 6 && argv[5][0] != '-'){
	  extendMinWT = strtod(argv[5],&qt);
	  if(qt == argv[5] || !(*qt==0 || isspace(*qt)) || extendMinWT > 1.0){
	    printf("%s 5th value %s is invalid (must be WT between 0.0 and 1.0)\n",argv[0],argv[5]);
	    fflush(stdout);exit(1);
	  }

	  if(argc >= 7 && argv[6][0] != '-'){
	    extendMaxOutlierKBsplit = extendMaxOutlierKB = strtod(argv[6],&qt);
	    if(qt == argv[6] || !(*qt==0 || isspace(*qt))){
	      printf("%s 6th value %s is invalid (must be Max outlier size in kb)\n",argv[0],argv[6]);
	      fflush(stdout);exit(1);
	    }

	    int shift = 0;
	    if(argc >= 8 && argv[7][0] != '-' && !strchr(argv[7],'.')){
	      extendMaxOutlierLabelsSplit = extendMaxOutlierLabels = strtol(argv[7],&qt,10);
	      if(qt == argv[7] || !(*qt==0 || isspace(*qt))){
		printf("%s 7th value %s is invalid (must be non-negative integer OR non-negative float with decimal point)\n",argv[0],argv[7]);
		fflush(stdout);exit(1);
	      }

	      shift++;

	      argv++;
	      argc--;
	    }

	    if(argc >= 8 && argv[7][0] != '-'){
	      extendMaxWT = strtod(argv[7],&qt);
	      if(qt == argv[7] || !(*qt==0 || isspace(*qt)) || extendMaxWT > 1.0){
		printf("%s %dth value %s is invalid (must be WT between 0.0 and 1.0)\n",argv[-shift],7+shift,argv[7]);
		fflush(stdout);exit(1);
	      }

	      if(argc >= 9 && argv[8][0] != '-'){
		extendMaxOutlierKBsplit = strtod(argv[8],&qt);
		if(qt == argv[8] || !(*qt==0 || isspace(*qt))){
		  printf("%s %dth value %s is invalid (must be Max outlier size in kb)\n",argv[-shift],8+shift,argv[8]);
		  fflush(stdout);exit(1);
		}
		
		if(argc >= 10 && argv[9][0] != '-'){
		  extendMaxOutlierLabelsSplit = strtol(argv[9],&qt,10);
		  if(qt == argv[9] || !(*qt==0 || isspace(*qt))){
		    printf("%s %dth value %s is invalid (must be non-negative integer)\n",argv[-shift],9+shift,argv[9]);
		    fflush(stdout);exit(1);
		  }

		  argv++;
		  argc--;
		}
		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 5;
	argc -= 5;
	continue;
      }
      if(!strcmp(argv[0],"-extendWTdup")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	extendWTdup = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || extendWTdup > 1){
	  printf("%s value %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-extendWTend")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	extendWTend = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || extendWTend > 1){
	  printf("%s value %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }

      break;
    
    case 'f':
      if(!strcmp(argv[0],"-fasta")){
	GenerateFasta = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-finalsort-sitesdec")){
	FinalSort = 4;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-forcemerge")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a Filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ForceMerge = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-floatcov")){
	if(argc < 2 || argv[1][0] == '-'){
	  floatcov = 1;
	  argv++;
	  argc--;
	} else {
	  floatcov = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || floatcov < 0 || floatcov > 1){
	    printf("%s value of %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv += 2;
	  argc -= 2;
	}
	continue;
      }
      if(!strcmp(argv[0],"-fastmode")){
	fast_mode= atoi(argv[1]);
	argv+=2;
	argc-=2;
	continue;
      }
      if(!strcmp(argv[0],"-fastmodetolerance")){
       fast_mode_tolerance= atof(argv[1]);
       argv+=2;
       argc-=2;
       continue;
      }
      if(!strcmp(argv[0],"-firstalignments")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	FirstAlignments = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || FirstAlignments <= 0){
	  printf("%s value is invalid (must be between +ve integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-f") || !strcmp(argv[0],"-force")){
	ForceOverwrite = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-firstPairs")){
	if(argc < 2){
	  printf("-first flag must be followed by a number\n");
	  fflush(stdout);exit(1);
	}
	CfirstPairs = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || CfirstPairs < 0 || CfirstPairs > 2){
	  printf("-first value is invalid (must be 0,1 or 2):%s\n", argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-first")){
	if(argc < 2){
	  printf("-first flag must be followed by a number\n");
	  fflush(stdout);exit(1);
	}
	Cfirst = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || Cfirst == 0){
	  printf("-first value is invalid (must be != 0):%s\n", argv[1]);
	  fflush(stdout);exit(1);
	}
	if(VERB>=2){
	  printf("Cfirst=%d : %s %s\n",Cfirst,argv[0],argv[1]);
	  fflush(stdout);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'g':
      if(!strcmp(argv[0],"-grouped")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	groupfile = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'h':
      if(!strcmp(argv[0],"-hashSF")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by 2 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	
	HashSF = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid (must be non-negative number)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	HashSFscale = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid (must be non-negative number)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-hashSD")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by 2 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	
	HashSD = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid (must be non-negative number)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	HashSDscale = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid (must be non-negative number)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-hashSR")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by 2 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	
	HashSR = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid (must be non-negative number)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	HashSRscale = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid (must be non-negative number)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-hashkeys")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	hashkeys = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-hashScaleDelta")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0, 1, 2 or 3\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	hashScaleDelta = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || hashScaleDelta > 3){
	  printf("%s value %s is invalid (must be 0, 1 , 2 or 3)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-hashGrouped")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by 2 non-negative integers\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	ShiftOffset1 = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value is invalid (must be non-negative integer): %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	GroupedScore = strtol(argv[2],&qt,10);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value is invalid (must be non-negative integer): %s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-hashT2")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	HashT2 = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || HashT2 > 1 ){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-hashGC")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative integer\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	HashGC = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout); exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-hashbest")){
	if(argc < 2 || (argv[1][0] == '-' && !isdigit(argv[1][1]))){
	  printf("%s must be following by integer (use -1 to disable)\n",argv[0]);
	  fflush(stdout); exit(1);
	}
	HashBest = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid (must be integer >= -1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-hashcolor")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	hashcolor = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || hashcolor < 0 || usecolor > MAXCOLOR){
	  printf("%s value is invalid (must be between 0 and %d):%s\n",argv[0],MAXCOLOR,argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-hashmaxmem")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by Maximum memory in Gigabytes\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	HMaxMem = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || HMaxMem <= 0){
	  printf("%s value has invalid syntax or value (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-hashmatrix")){
	HashMatrix = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-hashoffset")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	OffsetSpread = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || OffsetSpread < 0){
	  printf("%s value is invalid (must be between 0 and %d):%s\n",argv[0],OFFSET_SPREAD,argv[1]);
	  fflush(stdout);exit(1);
	}
	if(OffsetSpread > OFFSET_SPREAD){
	  printf("%s value reduced from %d to %d (see OFFSET_SPREAD in hash.h to increase max value)\n",argv[0],OffsetSpread,OFFSET_SPREAD);
	  fflush(stdout);
	  OffsetSpread = OFFSET_SPREAD;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-hashMultiMatch")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	HashMultiMatch = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || HashMultiMatch < 0){
	  printf("%s value is invalid (must be non-negative integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	HashMultiMatchMax = 0;
	HashMultiMatchMaxDelta = 0;

	if(argc >= 3 && argv[2][0] != '-'){
	  HashMultiMatchMax = strtol(argv[2],&qt,10);
	  if(qt==argv[2] || !(*qt == 0 || isspace(*qt)) || HashMultiMatchMax < 0){
	    printf("%s 2nd value is invalid (must be non-negative integer):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 4 && argv[3][0] != '-'){
	    HashMultiMatchMaxDelta = strtol(argv[3],&qt,10);
	    if(qt==argv[3] || !(*qt == 0 || isspace(*qt)) || HashMultiMatchMaxDelta < 0){
	      printf("%s 3rd value is invalid (must be non-negative integer):%s\n",argv[0],argv[3]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-hashMaxCov")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	hashMaxCov = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || hashMaxCov <= 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-hashtxt")){
	HashTxt = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-hashgen")){
	HashGen = 1;

	if(argc >= 2 && argv[1][0] != '-'){
	  if(argc < 9 || argv[2][0] == '-' || argv[3][0] == '-' || argv[4][0] == '-' || argv[5][0] == '-' || argv[6][0] == '-' || argv[7][0] == '-' || argv[8][0] == '-'){
	    printf("%s flag must be followed by 8 or 9 non-negative numbers (or none)\n",argv[0]);
	    fflush(stdout);exit(1);
	  }
	  int cnt = 1;

	  HashWin = strtol(argv[cnt],&qt,10);
	  if(qt == argv[cnt] || !(*qt==0 || isspace(*qt))){
	    printf("%s Win value of %s is invalid (must be >= 0)\n",argv[0],argv[cnt]);
	    fflush(stdout);exit(1);
	  }
	  cnt++;

	  HashScore = strtol(argv[cnt],&qt,10);
	  if(qt == argv[cnt] || !(*qt==0 || isspace(*qt)) || HashScore <= 0){
	    printf("%s HashScore value of %s is invalid (must be > 0)\n",argv[0],argv[cnt]);
	    fflush(stdout);exit(1);
	  }
	  cnt++;

	  SDmax = strtod(argv[cnt],&qt);
	  if(qt == argv[cnt] || !(*qt==0 || isspace(*qt)) || SDmax <= 0.0){
	    printf("%s SDmax value of %s is invalid (must be > 0)\n",argv[0],argv[cnt]);
	    fflush(stdout);exit(1);
	  }
	  cnt++;

	  SDrms = strtod(argv[cnt],&qt);
	  if(qt == argv[cnt] || !(*qt==0 || isspace(*qt)) || SDrms < 0.0){
	    printf("%s SDrms value of %s is invalid (must be >= 0)\n",argv[0],argv[cnt]);
	    fflush(stdout);exit(1);
	  }
	  cnt++;

	  RelErr = strtod(argv[cnt],&qt);
	  if(qt == argv[cnt] || !(*qt==0 || isspace(*qt)) || RelErr <= 0.0){
	    printf("%s SDrms value of %s is invalid (must be > 0)\n",argv[0],argv[cnt]);
	    fflush(stdout);exit(1);
	  }
	  cnt++;

	  OffsetKB = strtod(argv[cnt],&qt);
	  if(qt == argv[cnt] || !(*qt==0 || isspace(*qt)) || OffsetKB <= 0.0){
	    printf("%s OffsetKB value of %s is invalid (must be > 0)\n",argv[0],argv[cnt]);
	    fflush(stdout);exit(1);
	  }
	  cnt++;

	  NumErrI = strtol(argv[cnt],&qt,10);
	  if(qt == argv[cnt] || !(*qt==0 || isspace(*qt)) || NumErrI < 0){
	    printf("%s NumErrI value of %s is invalid (must be >= 0)\n", argv[0],argv[cnt]);
	    fflush(stdout);exit(1);
	  }
	  cnt++;

	  NumErrQ = strtol(argv[cnt],&qt,10);
	  if(qt == argv[cnt] || !(*qt==0 || isspace(*qt)) || NumErrQ < 0){
	    printf("%s NumErrQ value of %s is invalid (must be >= 0)\n", argv[0],argv[cnt]);
	    fflush(stdout);exit(1);
	  }
	  cnt++;
	  if(DEBUG) assert(cnt==9);
	  if(argc >= 10 && argv[cnt][0] != '-'){/* 9th arg  is NumResQ */
	    NumResQ = strtol(argv[cnt],&qt,10);
	    if(qt == argv[cnt] || !(*qt==0 || isspace(*qt)) || NumResQ < 0){
	      printf("%s NumResQ value of %s is invalid (must be >= 0)\n",argv[0],argv[cnt]);
	      fflush(stdout);exit(1);
	    }
	    cnt++;
	  }
	
	  argv += cnt;
	  argc -= cnt;
	} else {
	  argv++;
	  argc--;
	}
	
	if(HashWin <= 0)
	  HashGen = 0;

	continue;
      }
      if(!strcmp(argv[0],"-hash")){
	if(argc >= 3 && argv[1][0] != '-'){/* 2 args */
	  if(argv[2][0] == '-'){
	    printf("-hash must be following by 0 or 2 args\n%s %s %s\n", argv[0],argv[1],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  hash_filename = argv[1];
	  hash_threshold = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || hash_threshold < 0){
	    printf("hash threshold value of %s is invalid\n",argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv += 3;
	  argc -= 3;
	} else {
	  if(argc >=2 && argv[1][0] != '-'){
	    printf("-hash must be following by 0 or 2 args\n%s %s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  hash_filename = (char *)(-1);/* placeholder : will be set later, since -hashgen, -o and -id may occur later on commandline */
	  argv++;
	  argc--;
	}
	continue;
      }
      if(!strcmp(argv[0],"-hashsplit")){
	hashsplit = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-hashrange")){
	HASHRANGE = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  HASHRANGE = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || HASHRANGE < 0 || HASHRANGE > 1){
	    printf("%s float value of %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-hashdelta")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	hashdelta = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || hashdelta < 0.0){
	  printf("%s float value of %s is invalid\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	hashdeltaAdjust = 0.0;
	if(argc >= 3 && argv[2][0] != '-'){
	  hashdeltaAdjust = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || hashdeltaAdjust < 0.0){
	    printf("%s 2nd float value of %s is invalid\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  hashdeltaLim = 0.0;
	  if(argc >= 4 && argv[3][0] != '-'){
	    hashdeltaLim = strtod(argv[3],&qt);
	    if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || hashdeltaLim < 0.0){
	      printf("%s 3rd float value of %s is invalid\n",argv[0],argv[3]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'i':
      if(!strcmp(argv[0],"-incID")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative 63-bit integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	IncID = strtoll(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || IncID < 0){
	  printf("%s value %s is invalid (must be non-negative 63-bit integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-inversion")){
	if(argc<3){
	  printf("%s flag must be followed two numbers\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	InversionLeft = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || InversionLeft < 0.0){
	  printf("%s left value is invalid (must be >= 0):left=%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	InversionRight = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || InversionRight <= InversionLeft){
	  printf("%s right value is invalid (must be > left value):left=%s right=%s\n",argv[0],argv[1],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-indel")){
	InDel = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-indelends")){
	InDelEnds = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-indelMinConf")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	InDelMinConf = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid (must be non-negative number):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-indelChiSq")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by 2 to 15 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	IndelChiSq = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || IndelChiSq > 1.0){
	  printf("%s 1st value of %s is invalid (must be number between 0.0 and 1.0)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	IndelFrac = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt == 0 || isspace(*qt)) || IndelChiSq > 1.0){
	  printf("%s 2nd value of %s is invalid (must be number between 0.0 and 1.0)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	
	IndelMaxSize = 30.0;
	IndelMaxCov = 60.0;
	IndelMinSf = 0.300 * 1000.0;
	IndelMinSr = 0.02 * 1000.0;

	IndelMinInterval = 2.0;

	IndelSdRatio = 1.2;
	IndelMinCov = 3.0;
	IndelMinCovScale = 1.5;

	IndelWin = 0;
	IndelFragile = 999.0;
        IndelSegDup = 999.0;
        IndelChimNorm = 0.0;

	IndelMinCovRatio = 0.0;

	if(argc >= 4 && argv[3][0] != '-'){
	  IndelMaxSize = strtod(argv[3],&qt);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value of %s is invalid (must be non-negative size in kb)\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 5 && argv[4][0] != '-'){
	    IndelMaxCov = strtod(argv[4],&qt);
	    if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th value of %s is invalid (must be non-negative number)\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }
	    if(argc >= 6 && argv[5][0] != '-'){
	      IndelMinSf = strtod(argv[5],&qt);
	      if(qt == argv[5] || !(*qt==0 || isspace(*qt))){
		printf("%s 5th value of %s is invalid (must be non-negative number)\n",argv[0],argv[5]);
		fflush(stdout);exit(1);
	      }
	      IndelMinSf *= 1000.0;/* Arg is in kb, MolSd is in bp */

	      if(argc >= 7 && argv[6][0] != '-'){
		IndelMinSr = strtod(argv[6],&qt);
		if(qt == argv[6] || !(*qt==0 || isspace(*qt))){
		  printf("%s 6th value of %s is invalid (must be non-negative number)\n",argv[0],argv[6]);
		  fflush(stdout);exit(1);
		}
		IndelMinSr *= 1000.0;/* Interval is in kb, MolSd is in bp */

		if(argc >= 8 && argv[7][0] != '-'){
		  IndelMinInterval = strtod(argv[7],&qt);
		  if(qt == argv[7] || !(*qt==0 || isspace(*qt))){
		    printf("%s 7th value of %s is invalid (must be non-negative number)\n",argv[0],argv[7]);
		    fflush(stdout);exit(1);
		  }

		  if(argc >= 9 && argv[8][0] != '-'){
		    IndelSdRatio = strtod(argv[8],&qt);
		    if(qt == argv[8] || !(*qt==0 || isspace(*qt))){
		      printf("%s 8th value of %s is invalid (must be non-negative number)\n",argv[0],argv[8]);
		      fflush(stdout);exit(1);
		    }

		    if(argc >= 10 && argv[9][0] != '-'){
		      IndelMinCov = strtod(argv[9],&qt);
		      if(qt == argv[9] || !(*qt==0 || isspace(*qt))){
			printf("%s 9th value of %s is invalid (must be non-negative number)\n",argv[0],argv[9]);
			fflush(stdout);exit(1);
		      }
		    
		      if(argc >= 11 && argv[10][0] != '-'){
			IndelMinCovScale = strtod(argv[10],&qt);
			if(qt == argv[10] || !(*qt==0 || isspace(*qt))){
			  printf("%s 10th value of %s is invalid (must be non-negative number)\n",argv[0],argv[10]);
			  fflush(stdout);exit(1);
			}
			if(argc >= 12 && argv[11][0] != '-'){
			  IndelWin = strtol(argv[11],&qt,10);
			  if(qt == argv[11] || !(*qt==0 || isspace(*qt))){
			    printf("%s 11th value %s is invalid (must be non-negative integer)\n",argv[0],argv[11]);
			    fflush(stdout);exit(1);
			  }			  
			  if(argc >= 13 && argv[12][0] != '-'){
			    IndelFragile = strtod(argv[12],&qt);
			    if(qt == argv[12] || !(*qt==0 || isspace(*qt))){
			      printf("%s 12th value %s is invalid (must be non-negative number))\n",argv[0],argv[12]);
			      fflush(stdout);exit(1);
			    }
			    if(argc >= 14 && argv[13][0] != '-'){
			      IndelSegDup = strtod(argv[13],&qt);
			      if(qt == argv[13] || !(*qt==0 || isspace(*qt))){
				printf("%s 13th value %s is invalid (must be non-negative number)\n",argv[0],argv[13]);
				fflush(stdout);exit(1);
			      }
			      if(argc >= 15 && argv[14][0] != '-'){
				IndelChimNorm = strtod(argv[14],&qt);
				if(qt == argv[14] || !(*qt==0 || isspace(*qt))){
				  printf("%s 14th value %s is invalid (must be non-negative number)\n",argv[0],argv[14]);
				  fflush(stdout);exit(1);
				}
				if(argc >= 16 && argv[15][0] != '-'){
				  IndelMinCovRatio = strtod(argv[15],&qt);
				  if(qt == argv[15] || !(*qt==0 || isspace(*qt))){
				    printf("%s 15th value %s is invalid (must be non-negative number)\n",argv[0],argv[15]);
				    fflush(stdout);exit(1);
				  }

				  argv++;
				  argc--;
				}
				argv++;
				argc--;
			      }
			      argv++;
			      argc--;
			    }
			    argv++;
			    argc--;
			  }
			  argv++;
			  argc--;
			}
			argv++;
			argc--;
		      }
		      argv++;
		      argc--;
		    }
		    argv++;
		    argc--;
		  }
		  argv++;
		  argc--;
		}
		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-indelNoOverlap")){
	InDelNoOverlap = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  InDelNoOverlap = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || InDelNoOverlap > 3){
	    printf("%s 1st value has invalid syntax or value (must be 0, 1, 2 or 3):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 3 && argv[2][0] != '-'){
	    double Pvalue = strtod(argv[2],&qt);
	    if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || Pvalue >= 1.0 || Pvalue <= 0.0){
	      printf("%s 2nd value has invalid syntax : must be a Pvalue between 0.0 and 1.0: %s\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }
	    InDelNoOverlapT = -log(Pvalue)/log(10.0);

	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-insertThreads")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by number of threads\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	HashInsertThreads = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || HashInsertThreads <= 0){
	  printf("%s value has invalid syntax or value (must be >= 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-i")){
	if(argc<2){
	  printf("-i flag must be followed by .vfixx OR .cmap filename\n");
	  fflush(stdout);exit(1);
	}
	if(num_files >= MAXFILES){
	  printf("Too many input files : increase MAXFILES=%d\n", MAXFILES);
	  fflush(stdout);exit(1);
	}
	vfixx_filename[num_files++] = strdup(argv[1]);
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-if")){
	if(argc<2 || argv[1][0] == '-'){
	  printf("-if flag must be followed by the name of a file that contains input file names (one per line)\n");
	  fflush(stdout);exit(1);
	}

	errno = 0;
	FILE *fp = Fopen(argv[1], "r", 65);
	/*	int maxcnt = 65;
	for(int cnt = 0; cnt <= maxcnt && (fp = fopen(argv[1],"r")) == NULL; cnt++){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("Cannot read file %s (to get file names for %s), retry=%d/%d:errno=%d:%s\n",argv[1],argv[0],cnt,maxcnt,eno,err);
	  printf("Current Time= %s\n", currentDateTime());
	  fflush(stdout);
	  sleep(1);
	  errno = 0;
	  }*/
	if(fp == NULL)
	  exit(1);
	int orignumfiles = num_files;
	int linecnt = 1;
	for(;fgets(buf,LINESIZ,fp) != NULL; linecnt++){
	  int len = strlen(buf);
	  if(len <= 0) continue;
	  if(buf[0]=='#' || buf[0]=='\n' || buf[0]=='\r')continue;
	  if(len >= min(LINESIZ,PATH_MAX)-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	    printf("Line too long or not terminated in -if input file %s on line %d (len=%d,buf[len-1]=%c):\n%s\n",
		   argv[1],linecnt,len,buf[len-1],buf);
	    fflush(stdout);exit(1);
	  }
	  buf[--len] = '\0';
	  while(len > 0 && isspace(buf[len-1]))
	    buf[--len] = '\0';
	  if(num_files < MAXFILES){
	    if(VERB>=2){
	      printf("vfixx_filename[%d] = %s (from line %d in %s)\n",num_files,buf,linecnt,argv[1]);
	      fflush(stdout);
	    }
#ifndef WIN32
	    if(buf[0] != '/'){/* relative file names are interpreted relative to directory in which -if file is located */
	      char filename[BUFSIZ];
	      strcpy(filename, argv[1]);
	      char *pt = strrchr(filename,'/');
	      if(pt == NULL)
		vfixx_filename[num_files++] = strdup(buf);
	      else {
		sprintf(&pt[1],"%s",buf);
		vfixx_filename[num_files++] = strdup(filename);
	      }
	    } else
#endif
	      vfixx_filename[num_files++] = strdup(buf);
	  } else
	    num_files++;
	}
	fclose(fp);
	if(num_files > MAXFILES){
	  printf("Too many input files : increase MAXFILES=%d to at least %d\n",
		 MAXFILES, num_files);
	  fflush(stdout);exit(1);
	}
	if(num_files == orignumfiles){
	  printf("WARNING: input file list %s %s has no file names\n",argv[0],argv[1]);
	  fflush(stdout);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-id")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	startCMapID = CMapID = strtoll(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || CMapID <= 0){
	  printf("%s value has invalid syntax or value (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'k':
      if(!strcmp(argv[0],"-keepsplitmaps")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	keepsplitmaps = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || keepsplitmaps > 1){
	  printf("%s value %s is invalid (must be 0 or 1)\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'l':
      if(!strcmp(argv[0],"-local")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by Chimerism probability and False Discovery Rate\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	Ch = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || Ch < 0.0){
	  printf("chimerism probability is invalid (must be >= 0):%s\n",argv[1]);
	  fflush(stdout);exit(1);
	}
	ChFdr = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || ChFdr <= 0.0){
	  printf("chimerism False Discovery Rate is invalid (must be > 0):%s\n",argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      break;
    case 'm':
      if(!strcmp(argv[0],"-maxExtend")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxExtend = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || MaxExtend < 0.0){
	  printf("%s valuie %s is invalid : must be non-negative size in kb\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maxEnd")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MaxEnd = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || MaxEnd < 0.0){
	  printf("%s value %s is invalid : must be non-negative size in kb\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maxInterval")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1-2 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MaxInterval = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value %s is invalid : must be non-negative size in kb\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	refMaxInterval = 0;
	if(argc >= 3 && argv[2][0] != '-'){
	  refMaxInterval = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt == 0 || isspace(*qt)) || refMaxInterval > 1){
	    printf("%s 2nd value %s is invalid : must be 0 or 1\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  
	  argv++;
	  argc--;
	}


	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-minEnd")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MININTERVAL = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || MININTERVAL < 0.0){
	  printf("%s value %s is invalid : must be non-negative size in kb\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-minSNRestimate")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by 2 values\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	minSNRestimate = 1;
	SNRmin = strtod(argv[1],&qt);
	if(argv[1]==qt || !(*qt==0 || isspace(*qt)) || SNRmin < 0.0){
	  printf("%s : 1st value %s is not valid (must be non-negative number)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	SNRmax = strtod(argv[2],&qt);
	if(argv[2]==qt || !(*qt==0 || isspace(*qt)) || SNRmax < SNRmin){
	  printf("%s : 2nd value %s is not valid (must be >= 1st value %s)\n",argv[0],argv[2],argv[1]);
	  fflush(stdout);exit(1);
	}
	SNRstep = -1.0;
	if(argc >= 4 && argv[3][0] != '-'){
	  SNRstep = strtod(argv[3],&qt);
	  if(argv[3]==qt || !(*qt==0 || isspace(*qt)) || SNRstep < 0.0){
	    printf("%s : 3rd value %s is not valid (must be non-negative number)\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }

	  SNRmaxmaps = 5000;
	  if(argc >= 5 && argv[4][0] != '-'){
	    SNRmaxmaps = strtol(argv[4],&qt,10);
	    if(argv[4]==qt || !(*qt==0 || isspace(*qt)) || SNRmaxmaps < 0){
	      printf("%s : 4th value %s is not valid (must be +ve number)\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }	
	    if(argc >= 6 && argv[5][0] != '-'){
	      SNRminratio = strtod(argv[5],&qt);
	      if(argv[5]==qt || !(*qt==0 || isspace(*qt)) || SNRminratio < 0.0){
		printf("%s : 5th value %s is not valid (must be non-negative fraction)\n",argv[0],argv[5]);
		fflush(stdout);exit(1);
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  
	  argv++;
	  argc--;
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-minSNR")){
	if(argc < 2){
	  printf("%s must be followed by at least 1 value (up to one per color)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	int i = 0;
	for(; i < MAXCOLOR; i++){
	  if(argc <= 1+i || (i >= 1 && argv[i+1][0] == '-' && argv[i+1][1] && !isdigit(argv[i+1][1]) && argv[i+1][1] != '.'))
	    break;
	  minSNR[i] = strtod(argv[i+1],&qt);
	  if(argv[i+1]==qt || !(*qt==0 || isspace(*qt)) || minSNR[i] < 0.0){
	    if(i > 0)
	      printf("The %d'th %s value %s is not valid (must be non-negative number)\n",i+1,argv[0],argv[i+1]);
	    else
	      printf("The 1st %s value %s is not valid (must be non-negative number)\n",argv[0],argv[i+1]);
	    fflush(stdout);exit(1);
	  }
	}
	for(int nc = i; nc < MAXCOLOR; nc++)
	  minSNR[nc] = minSNR[i-1];
	argv += 1+i;
	argc -= 1+i;
	continue;
      }
      if(!strcmp(argv[0],"-maxPSFWidth")){
	if(!colors){
	  printf("-colors must be specified before %s\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(argc < 2){
	  printf("%s must be followed by at least 1 value (up to one per color)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	int i = 0;
	for(; i < MAXCOLOR; i++){
	  if(argc <= 1+i || (i >= 1 && argv[i+1][0] == '-' && argv[i+1][1] && !isdigit(argv[i+1][1]) && argv[i+1][1] != '.'))
	    break;
	  maxPSFWidth[i] = strtod(argv[i+1],&qt);
	  if(argv[i+1]==qt || !(*qt==0 || isspace(*qt)) || maxPSFWidth[i] < 0.0){
	    if(colors > 1)
	      printf("The %d. %s value %s is not valid (must be non-negative number)\n",i+1,argv[0],argv[i+1]);
	    else
	      printf("The %s value %s is not valid (must be non-negative number)\n",argv[0],argv[i+1]);
	    fflush(stdout);exit(1);
	  }
	}
	for(int nc = i; nc < MAXCOLOR; nc++)
	  maxPSFWidth[nc] = maxPSFWidth[i-1];
	argv += 1+i;
	argc -= 1+i;
	continue;
      }
      if(!strcmp(argv[0],"-maxthreads")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MaxThreads = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || MaxThreads <= 0){
	  printf("%s value %s is invalid : must be +ve integer\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maptype")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a number (0 or 1)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MapType = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MapType < 0 || MapType > 1){
	  printf("%s value %s is invalid : must be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-ma")){
	if(argc < 4){
	  printf("%s flag must be followed by 3 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(colors != 2){
	  printf("%s flag must be preceded by -colors 2\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	MA_mean = strtod(argv[1],&qt);	
	if(qt == argv[1] || (*qt && !isspace(*qt))){
	  printf("%s 1st value of %s has invalid syntax\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	MA_SF = strtod(argv[2],&qt);
	if(qt == argv[2] || (*qt && !isspace(*qt))){
	  printf("%s 2nd value of %s has invalid syntax\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	MA_SD = strtod(argv[3],&qt);
	if(qt == argv[3] || (*qt && !isspace(*qt))){
	  printf("%s 3rd value of %s has invalid syntax\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-mapped")){
	if(argc<2){
	  printf("-mapped flag must be followed by output filename prefix\n");
	  fflush(stdout);exit(1);
	}
	MappedPrefix = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-mapped-unsplit")){
	if(argc<2){
	  printf("-mapped-unsplit flag must be followed by 0 or 1\n");
	  fflush(stdout);exit(1);
	}
	MappedUnsplit = atoi(argv[1]);
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-merge")){
	CmapMerge = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-minoverlap")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative float\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinOverlapRatio = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt))){
	  printf("%s value has invalid syntax:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-minlen")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by minimum -i Map length in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinLen = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MinLen < 0.0){
	  printf("%s value has invalid syntax or value:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maxlen")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by maximum -i Map length in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxLen = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxLen < 0.0){
	  printf("%s value has invalid syntax or value:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-minsites")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by minimum -i Map sites\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinSites = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MinSites < 0){
	  printf("%s value has invalid syntax or value (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maxsites")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by maximum -i Map sites\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxSites = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxSites < 0){
	  printf("%s value has invalid syntax or value (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	MaxSites2 = 0;
	if(argc >= 3 && argv[2][0] != '-'){
	  MaxSites2 = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || MaxSites2 < 0){
	    printf("%s 2nd value has invalid syntax or value (must be >= 0):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maxEndSiteDensity")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by 2 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxEndSiteDensity = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxEndSiteDensity < 0.0){
	  printf("%s 1st value has invalid syntax or value (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	MaxEndSiteDensityLength = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || MaxEndSiteDensityLength < 0.0){
	  printf("%s 2nd value has invalid syntax or value (must be >= 0):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-maxSiteDensity")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by maximum -i Map sites per 100kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxSiteDensity = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxSiteDensity < 0){
	  printf("%s value has invalid syntax or value (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-minSiteDensity")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by minimum -i Map sites per 100kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinSiteDensity = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MinSiteDensity < 0){
	  printf("%s value has invalid syntax or value (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maxContigSiteDensity")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by maximum -i Map sites per 100kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxContigSiteDensity = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxContigSiteDensity < 0){
	  printf("%s value has invalid syntax or value (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maxmemIncrease")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by Reserve Memory in Gigabytes\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxMemIncrease = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxMemIncrease < 0){
	  printf("%s value has invalid syntax or value (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maxmem")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by Maximum memory in Gigabytes\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxMem = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxMem <= 0){
	  printf("%s value has invalid syntax or value (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maxvirtmem")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by Maximum virtual memory in Gigabytes (0 for unlimited)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxVirtMem = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxVirtMem < 0.0){
	  printf("%s value has invalid syntax or value (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-mres_override")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	mres_override = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || mres_override < 0 || mres_override > 1){
	  printf("%s value has invalid syntax (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-mres")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by resolution in multiples of 0.5 kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	mres_defined = 1;
	mres = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || mres < 0.0){
	  printf("-mres value has invalid syntax (must be non-negative number):%s\n",argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-mresSD")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by resolution SD in multiples of 0.5 kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	mresSD = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || mresSD < 0.0){
	  printf("%s value has invalid syntax:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'n':
      if(!strcmp(argv[0],"-nostat")){
	if(argc < 2 || argv[1][0] == '-'){
	  NoStat = 1;
	  argv++;
	  argc--;
	} else {
	  NoStat = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || NoStat < 0 || NoStat > 1){
	    printf("%s value of %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv += 2;
	  argc -= 2;
	}
	continue;
      }
      if(!strcmp(argv[0],"-nosplit")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by chimeric nosplit level (0 1 or 2)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	NoSplitSet = 1;
	NoSplit = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || NoSplit < 0 || NoSplit > 2){
	  printf("%s value of %s is invalid (must be between 0 and 2)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'o':
      if(!strcmp(argv[0],"-outlierBalanced")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PairMergeMaxOutlierBalanced = strtol(argv[1],&qt,0);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid : must be non-negative integer\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-orMask")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	orMask = strtol(argv[1],&qt,0);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid : must be non-negative number\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlierRef")){
	if(PoutlierCnt){
	  printf("WARNING : -outlierRef specified twice : ignoring first occurance with %d values\n", PoutlierCnt);
	  fflush(stdout);
	  PoutlierCnt = 0;
	}
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by Refinement outlier Pvalue(s)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	while(argc >= 2+PoutlierCnt && argv[1+PoutlierCnt][0] != '-'){
	  if(PoutlierCnt >= 16){
	    printf("-outlierRef : only 16 values (or less) is supported: excess arg = %s\n", argv[1+PoutlierCnt]);
	    fflush(stdout);exit(1);
	  }
	  double Pvalue = strtod(argv[1+PoutlierCnt],&qt);
	  if(qt == argv[1+PoutlierCnt] || !(*qt==0 || isspace(*qt)) || Pvalue < (PoutlierCnt > 0 ? PoutlierRefine[PoutlierCnt-1] : 0.0) || Pvalue > 1.0){
	    printf("%s : %d'th value  of %s is invalid (must be between 0 and 1 and in ascending order)\n",argv[0], PoutlierCnt+1, argv[1+PoutlierCnt]);
	    fflush(stdout);exit(1);
	  }
	  PoutlierRefine[PoutlierCnt++] = Pvalue;
	}
	if(VERB>=2){
	  for(int t = 0; t < PoutlierCnt; t++)
	    printf("PoutlierRefine[%d]= %0.6e\n",t,PoutlierRefine[t]);
	  fflush(stdout);
	}
	argv += 1 + PoutlierCnt;
	argc -= 1 + PoutlierCnt;
	continue;
      }
      if(!strcmp(argv[0],"-outlierRate")){
	outlierRate = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-outlierBC")){
	outlierBC = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-outlierMax")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	outlierMax = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || outlierMax <= 0.0){
	  printf("%s value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(!OUTLIER_MAX && outlierMax <= 1000.0){
	  printf("ERROR: %s %s not supported : enable OUTLIER_MAX in constants.h and recompile\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlierType0")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	OutlierType0 = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || OutlierType0 > 1.0){
	  printf("%s value of %s is invalid (must between 0 and 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlierType1")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	OutlierType1 = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || OutlierType1 > 1){
	  printf("%s value of %s is invalid (must between 0 and 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlierNorm0")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	outlierNorm0 = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || outlierNorm0 < 0 || outlierNorm0 > 2){
	  printf("%s value of %s is invalid (must be 0, 1 or 2)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlierNorm1")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	outlierNorm1 = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || outlierNorm1 < 0 || outlierNorm1 > 1){
	  printf("%s value of %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlierTypeRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	OutlierType = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || OutlierType < 0 || OutlierType > 1){
	  printf("%s value of %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlierNormRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	outlierNormRef = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || outlierNormRef < 0 || outlierNormRef > 2){
	  printf("%s value of %s is invalid (must be 0, 1 or 2)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlierTypeSwitch")){
	OutlierTypeSwitch = 1;
	OutlierTypeLabel = 1;
	OutlierTypeSize = 0;
	if(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-'){
	  OutlierTypeLabel = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || OutlierTypeLabel < 0){
	    printf("%s 1st value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  OutlierTypeSize = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || OutlierTypeSize < 0){
	    printf("%s 2nd value of %s is invalid (must be positive)\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv += 2;
	  argc -= 2;
	}
	if(VERB>=2){
	  printf("OutlierTypeSwitch=%d,OutlierTypeLabel=%d,OutlierTypeSize=%d\n",OutlierTypeSwitch,OutlierTypeLabel,OutlierTypeSize);
	  fflush(stdout);
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-outlierLambda")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a positive number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	outlierLambda = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || outlierLambda <= 0.0){
	  printf("%s value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlierLambdaRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a positive number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	outlierLambdaRef = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || outlierLambdaRef <= 0.0){
	  printf("%s value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlierLambdaSwitch")){
	outlierLambdaSwitch = 1;
	outlierLambdaLabel = 1e+30;
	outlierLambdaSize = 0.0;
	if(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-'){
	  outlierLambdaLabel = strtod(argv[1],&qt);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || outlierLambdaLabel <= 0.0){
	    printf("%s 1st value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  outlierLambdaSize = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || outlierLambdaSize <= 0.0){
	    printf("%s 2nd value of %s is invalid (must be positive)\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv += 2;
	  argc -= 2;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-outlierSwitch")){
	PoutlierSwitch = 1;
	if(!(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-')){
	  printf("%s must be followed by two numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PoutlierLabel = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PoutlierLabel <= 0.0){
	  printf("%s 1st value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	PoutlierSize = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || PoutlierSize <= 0.0){
	  printf("%s 2nd value of %s is invalid (must be positive)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-outlierExtend")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	outlierExtend = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || outlierExtend < 0){
	  printf("%s value of %s is invalid (must be between 0 and 2)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	outlierExtendLimit = 0;
	if(argc >= 3 && argv[2][0] != '-'){
	  outlierExtendLimit = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || outlierExtendLimit < 0){
	    printf("%s 2nd value of %s is invalid (must be positive integer)\n", argv[0], argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-o")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by output filename prefix\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	output_prefix = argv[1];
	if(strlen(output_prefix) > PATH_MAX - 100){
	  printf("output file name prefix must not exceed %d characters:\n%s %s\n",
		 PATH_MAX - 100,argv[0],argv[1]);
	  fflush(stdout);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-output-veto-filter")){
	if(argc<2){
	  printf("%s flag must be followed by output veto filter regular expression\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	output_veto_filter = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-output-filter")){
	if(argc<2){
	  printf("%s flag must be followed by output filter regular expression\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	output_filter = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-outlier")){
	if(argc<2){
	  printf("%s flag must be followed by outlier(blemish) probability per true positive site\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	Poutlier = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) ||  Poutlier < 0.0 || Poutlier > 1.0){
	  printf("%s value %s is invalid (must be between 0 and 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	setoutlier = true;
	if(argc >= 3 && argv[2][0] != '-'){
	  PoutlierS = strtod(argv[2],&qt);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) ||  PoutlierS < 0.0 || PoutlierS > 1.0){
	    printf("%s 2nd value %s is invalid (must be between 0 and 1)\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'p':
      if(!strcmp(argv[0],"-pairmergeOutlierChiSq")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s must be followed by 2 to 12 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	PairMergeOutlierChiSq = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || PairMergeOutlierChiSq > 1.0){
	  printf("%s 1st value of %s is invalid (must be number between 0.0 and 1.0)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	PairMergeOutlierFrac = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt == 0 || isspace(*qt)) || PairMergeOutlierChiSq > 1.0){
	  printf("%s 2nd value of %s is invalid (must be number between 0.0 and 1.0)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	
	PairMergeOutlierMaxSize = 35.0;
	PairMergeOutlierMaxCov = 60.0;
	PairMergeOutlierMinSf = 0.300 * 1000.0;
	PairMergeOutlierMinSr = 0.02 * 1000.0;

	PairMergeOutlierMinInterval = 2.0;

	PairMergeOutlierSdRatio = 1.2;
	PairMergeOutlierMinCov = 3.0;
	PairMergeOutlierMinCovScale = 1.5;

	PairMergeChimNorm = 0.0;
	PairMergeOutlierMinCovRatio = 0.0;

	if(argc >= 4 && argv[3][0] != '-'){
	  PairMergeOutlierMaxSize = strtod(argv[3],&qt);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value of %s is invalid (must be non-negative size in kb)\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 5 && argv[4][0] != '-'){
	    PairMergeOutlierMaxCov = strtod(argv[4],&qt);
	    if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th value of %s is invalid (must be non-negative number)\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }
	    if(argc >= 6 && argv[5][0] != '-'){
	      PairMergeOutlierMinSf = strtod(argv[5],&qt);
	      if(qt == argv[5] || !(*qt==0 || isspace(*qt))){
		printf("%s 5th value of %s is invalid (must be non-negative number)\n",argv[0],argv[5]);
		fflush(stdout);exit(1);
	      }
	      PairMergeOutlierMinSf *= 1000.0;/* Arg is in kb, MolSd is in bp */

	      if(argc >= 7 && argv[6][0] != '-'){
		PairMergeOutlierMinSr = strtod(argv[6],&qt);
		if(qt == argv[6] || !(*qt==0 || isspace(*qt))){
		  printf("%s 6th value of %s is invalid (must be non-negative number)\n",argv[0],argv[6]);
		  fflush(stdout);exit(1);
		}
		PairMergeOutlierMinSr *= 1000.0;/* Interval is in kb, MolSd is in bp */

		if(argc >= 8 && argv[7][0] != '-'){
		  PairMergeOutlierMinInterval = strtod(argv[7],&qt);
		  if(qt == argv[7] || !(*qt==0 || isspace(*qt))){
		    printf("%s 7th value of %s is invalid (must be non-negative number)\n",argv[0],argv[7]);
		    fflush(stdout);exit(1);
		  }

		  if(argc >= 9 && argv[8][0] != '-'){
		    PairMergeOutlierSdRatio = strtod(argv[8],&qt);
		    if(qt == argv[8] || !(*qt==0 || isspace(*qt))){
		      printf("%s 8th value of %s is invalid (must be non-negative number)\n",argv[0],argv[8]);
		      fflush(stdout);exit(1);
		    }

		    if(argc >= 10 && argv[9][0] != '-'){
		      PairMergeOutlierMinCov = strtod(argv[9],&qt);
		      if(qt == argv[9] || !(*qt==0 || isspace(*qt))){
			printf("%s 9th value of %s is invalid (must be non-negative number)\n",argv[0],argv[9]);
			fflush(stdout);exit(1);
		      }
		    
		      if(argc >= 11 && argv[10][0] != '-'){
			PairMergeOutlierMinCovScale = strtod(argv[10],&qt);
			if(qt == argv[10] || !(*qt==0 || isspace(*qt))){
			  printf("%s 10th value of %s is invalid (must be non-negative number)\n",argv[0],argv[10]);
			  fflush(stdout);exit(1);
			}
			if(argc >= 12 && argv[11][0] != '-'){
			  PairMergeChimNorm = strtod(argv[11],&qt);
			  if(qt == argv[11] || !(*qt==0 || isspace(*qt))){
			    printf("%s 11th value of %s is invalid (must be non-negative number)\n",argv[0],argv[11]);
			    fflush(stdout);exit(1);
			  }
			  if(argc >= 13 && argv[12][0] != '-'){
			    PairMergeOutlierMinCovRatio = strtod(argv[12],&qt);
			    if(qt == argv[12] || !(*qt==0 || isspace(*qt))){
			      printf("%s 12th value of %s is invalid (must be non-negative number)\n",argv[0],argv[12]);
			      fflush(stdout);exit(1);
			    }

			    argv++;
			    argc--;
			  }
			  argv++;
			  argc--;
			}
			argv++;
			argc--;
		      }
		      argv++;
		      argc--;
		    }
		    argv++;
		    argc--;
		  }
		  argv++;
		  argc--;
		}
		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-pairmergeOverlapped")){
	PairMergeOverlapped = -1.0;
	if(argc >= 2 && argv[1][0] != '-'){
	  PairMergeOverlapped = strtod(argv[1],&qt);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	    printf("%s 1st value of %s is invalid (must be min overlap length in kb)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }

	  if(argc >= 3 && argv[2][0] != '-'){
	    PairMergeOverlappedMargin = strtod(argv[2],&qt);
	    if(qt == argv[2] || !(*qt==0 || isspace(*qt))){	    
	      printf("%s 2nd value of %s is invalid (must be non-negative size in kb)\n",argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }

	    PairMergeOverlappedSmargin = 1.0;
	    PairMergeOverlappedNmargin = 1;
	    PairMergeOverlappedEmargin = 0.1;

	    if(argc >= 4 && argv[3][0] != '-'){
	      PairMergeOverlappedSmargin = strtod(argv[3],&qt);
	      if(qt == argv[3] || !(*qt==0 || isspace(*qt))){	    
		printf("%s 3rd value of %s is invalid (must be non-negative size in kb)\n",argv[0],argv[3]);
		fflush(stdout);exit(1);
	      }

	      if(argc >= 5 && argv[4][0] != '-'){
		PairMergeOverlappedNmargin = strtol(argv[4],&qt,10);
		if(qt == argv[4] || !(*qt==0 || isspace(*qt))){	    
		  printf("%s 4th value of %s is invalid (must be non-negative size in kb)\n",argv[0],argv[4]);
		  fflush(stdout);exit(1);
		}

		if(argc >= 6 && argv[5][0] != '-'){
		  PairMergeOverlappedSmargin = strtod(argv[5],&qt);
		  if(qt == argv[5] || !(*qt==0 || isspace(*qt))){	    
		    printf("%s 5th value of %s is invalid (must be non-negative size in kb)\n",argv[0],argv[5]);
		    fflush(stdout);exit(1);
		  }

		  argv++;
		  argc--;
		}
		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-pairmergeExcludeOutliers")){
	PairMergeExcludeOutliers = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  PairMergeExcludeOutliers = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PairMergeExcludeOutliers > 1){
	    printf("%s value of %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-pairmergeXMAP")){
	PairMergeXmap = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-pairmergeHmap")){
	PairMergeHmap = 1;
	TruncateHmap = -1.0;
	TruncateNHmap = -1;
	PairMergeEndHmap = MinMapLenHmap = 0.0;
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by 2 or 4 non-negative numbers\n", argv[0]);
	  fflush(stdout);exit(1);
	}

	TruncateHmap = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid (must be length in kb)\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	TruncateNHmap = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid (must be +ve integer)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	if(argc >= 4 && argv[3][0] != '-'){
	  PairMergeEndHmap = strtod(argv[3],&qt);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value %s is invalid (must be length in kb)\n", argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 5 && argv[4][0] != '-'){
	    MinMapLenHmap = strtod(argv[4],&qt);
	    if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("%s 4th value %s is invalid (must be length in kb)\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-pairmergeRepeat")){
	PairMergeRepeat = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-perLabelScore")){
	perLabelScore = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-preferNGS")){
	PairMergePreferNGS = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-pairsplitPV")){
	if(argc < 2 || argv[1][0] == '-'){
	  PairsplitPV = 1;
	  argv++;
	  argc--;
	} else {
	  PairsplitPV = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PairsplitPV < 0 || PairsplitPV > 1){
	    printf("%s value of %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv += 2;
	  argc -= 2;
	}
	continue;
      }
      if(!strcmp(argv[0],"-pairsplitRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  PairsplitRef = 1;
	  argv++;
	  argc--;
	} else {
	  PairsplitRef = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PairsplitRef < 0 || PairsplitRef > 1){
	    printf("%s value of %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv += 2;
	  argc -= 2;
	}
	continue;
      }
      if(!strcmp(argv[0],"-pairsplitExtend")){
	PairsplitExtend = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-pairsplitMinOutlier")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PairsplitMinOutlier = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || PairsplitMinOutlier < 0){
	  printf("%s value is invalid (must be non-negative):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-pairsplitMerge")){
	PairsplitMerge = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-partmapped")){
	if(argc < 2){
	  printf("-ref flag must be followed by output filename prefix\n");
	  fflush(stdout);exit(1);
	}
	PartMappedPrefix = argv[1];

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-pairmerge")){
	PairMerge = 0.001;
	if(argc >= 2 &&  argv[1][0] != '-'){
	  PairMerge = strtod(argv[1],&qt);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || PairMerge <= 0.0){
	    printf("%s value %s is invalid (must be > 0)\n", argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 3 && argv[2][0] != '-'){
	    PairMergeMaxEnd = strtod(argv[2],&qt);
	    if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || PairMergeMaxEnd < 0.0){
	      printf("%s 2nd value %s is invalid (must be >= 0)\n", argv[0],argv[2]);
	      fflush(stdout);exit(1);
	    }
	    if(argc >= 4 && argv[3][0] != '-'){
	      PairMergeMaxOutlier = strtod(argv[3],&qt);
	      if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || PairMergeMaxOutlier < 0.0){
		printf("%s 3rd value %s is invalid (must be >= 0)\n", argv[0],argv[3]);
		fflush(stdout);exit(1);
	      }
	      if(argc >= 5 && argv[4][0] != '-'){
		PairMergeMaxEndKB = strtod(argv[4],&qt);
		if(qt == argv[4] || !(*qt == 0 || isspace(*qt)) || PairMergeMaxEndKB < 0.0){
		  printf("%s 4th value %s is invalid (must be >= 0)\n",argv[0],argv[4]);
		  fflush(stdout);exit(1);
		}
		PairMergeMaxEndExpandKB = PairMergeMaxEndKB;
		PairMergeMaxEndN = 0;
		PairMergeMaxOutlierM = 999;

		if(argc >= 6 && argv[5][0] != '-'){
		  PairMergeMaxEndExpandKB = strtod(argv[5],&qt);
		  if(qt == argv[5] || !(*qt == 0 || isspace(*qt)) || PairMergeMaxEndExpandKB < 0.0){
		    printf("%s 5th value %s is invalid (must be >= 0)\n",argv[0],argv[5]);
		    fflush(stdout);exit(1);
		  }
		  if(argc >= 7 && argv[6][0] != '-'){
		    PairMergeMaxEndN = strtol(argv[6],&qt,10);
		    if(qt == argv[6] || !(*qt == 0 || isspace(*qt))){
		      printf("%s 6th value %s is invalid (must be >= 0)\n",argv[0],argv[6]);
		      fflush(stdout);exit(1);
		    }
		    if(argc >= 8 && argv[7][0] != '-'){
		      PairMergeMaxOutlierM = strtol(argv[7],&qt,10);
		      if(qt == argv[7] || !(*qt == 0 || isspace(*qt))){
			printf("%s 7th value %s is invalid (must be >= 0)\n",argv[0],argv[7]);
			fflush(stdout);exit(1);
		      }
		      argv++;
		      argc--;
		    }
		    argv++;
		    argc--;
		  }
		  argv++;
		  argc--;
		}
		argv++;
		argc--;
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-pairsplit")){
	if(argc<2){
	  printf("%s flag must be followed by local indel probability per true positive site\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PairSplit = 1;
	Psplit = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || Psplit < 0.0 || Psplit > 1.0){
	  printf("%s value %s is invalid (must be between 0 and 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-partial")){
	if(argc < 3 || argv[1][0]=='-' || argv[2][0] == '-'){
	  printf("%s flag must be following by 2 or 3 +ve integers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	bin = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || bin <= 0){
	  printf("%s bin value is invalid (must be +ve integer)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	numbins = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || numbins < bin || numbins < 2){
	  printf("%s numbins value is invalid (must be +ve integer and >= %d)\n",argv[0],max(2,bin));
	  fflush(stdout);exit(1);
	}
	if(argc >= 4 && argv[3][0] != '-'){
	  skipcnt = strtol(argv[3],&qt,10);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || skipcnt < 0){
	    printf("%s %s %s %s: skipcnt is invalid (must be >= 0)\n",argv[0],argv[1],argv[2],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      break;
    case 'q':
      if(!strcmp(argv[0],"-queryThreads")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by number of threads\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	HashQueryThreads = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || HashQueryThreads <= 0){
	  printf("%s value has invalid syntax or value (must be >= 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-querysplit")){
	if(argc<2){
	  printf("%s flag must be followed by local indel probability per true positive site\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	PairSplit = -1;
	Psplit = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || Psplit < 0.0 || Psplit > 1.0){
	  printf("%s value %s is invalid (must be between 0 and 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'r':
      if(!strcmp(argv[0],"-relaxBnds")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	relaxBnds = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || relaxBnds > 1){
	  printf("%s value %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-refineMinWT")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number <= 1.0\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	refineMinWT = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || refineMinWT > 1.0){
	  printf("%s 2nd value %s is invalid : must be non-negative number <= 1.0\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-refineWT")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1, 4 or 5 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	refineWT = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || refineWT > 1.0){
	  printf("%s 1st value %s is invalid : must be non-negative number <= 1.0\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 5 && argv[2][0] != '-' && argv[3][0] != '-' && argv[4][0] != '-'){
	  refineWT_N = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value %s is invalid : must be max endoutlier labels\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  refineWT_NL = strtod(argv[3],&qt);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value %s is invalid : must be max endoutlier size in kb\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }

	  refineWT_MaxOutlierKB = strtod(argv[4],&qt);
	  if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	    printf("%s 4th value %s is invalid : must be max outlier size (delta) in kb\n",argv[0],argv[4]);
	    fflush(stdout);exit(1);
	  }

	  refineWT_MaxOutlierLabels = 1000;

	  if(argc >= 6 && argv[5][0] != '-'){
	    refineWT_MaxOutlierLabels = strtol(argv[5],&qt,10);
	    if(qt == argv[5] || !(*qt==0 || isspace(*qt))){
	      printf("%s 5th value %s is invalid : must non-negative integer\n",argv[0],argv[5]);
	      fflush(stdout);exit(1);
	    }

	    argv++;
	    argc--;
	  }
	  argv += 3;
	  argc -= 3;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-rresDefer")){
	rresDefer = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  rresDefer = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || rresDefer > 1){
	    printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-refineDefer")){
	refineDefer = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  refineDefer = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || refineDefer > 1){
	    printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-resEstimate")){
	ResEstimate = 1;
	if(argc >= 2 && argv[1][0] != '-'){
	  ResEstimate = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || ResEstimate < 0){
	    printf("%s value has invalid syntax (must be non-negative integer):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-rres")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by resolution in pixels\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	rres = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || rres < 0.0){
	  printf("-mres value has invalid syntax (must be non-negative number):%s\n",argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-readLoc")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	HasLoc = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || HasLoc < 0 || HasLoc > 1){
	  printf("%s value has invalid syntax (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-resbias")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by cutoff size and Number of linear segments\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	maxresbias = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || maxresbias < 0.0){
	  printf("%s 1st arg of %s is invalid (must be >= 0.0)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	origResBins[0] = ResBins[0] = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || ResBins[0] <= 0 || ResBins[0] > RESBINS){
	  printf("%s 2nd arg of %s is invalid (must be between 1 and %d)\n",argv[0],argv[2],RESBINS);
	  fflush(stdout);exit(1);
	}
	for(int c = 1; c < MAXCOLOR; c++)
	  origResBins[c] = ResBins[c] = ResBins[0];
	/* set reasonable default values corresponding to no bias */
	for(int c = 0; c < MAXCOLOR; c++){
	  ResBins[c] = 1;
	  resbiasX[c][0] = mres * 0.5;
	  resbiasX[c][1] = maxresbias;
	  for(int Bin = 0; Bin <= ResBins[c]; Bin++)
	    resbias[c][Bin] = 0.0;
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-relerr")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a ratio\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	RefRepeatErr = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RefRepeatErr <= 0.0){
	  printf("%s value of %s is invalid (must be > 0)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-randomize")){
	CmapSort = 6;
	if(argc >= 2 && argv[1][0] != '-'){
	  RandSeed = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RandSeed <= 0){
	    printf("%s Random Seed value of %s is invalid (must be > 0)\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-readparameters")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by .errbin filename\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	parametersfile = strdup(argv[1]);
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-range")){
	if(argc<3){
	  printf("%s flag must be followed two numbers\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	RangeLeft = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || RangeLeft < 0.0){
	  printf("%s left value is invalid (must be >= 0):left=%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	RangeRight = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || RangeRight <= RangeLeft){
	  printf("%s right value is invalid (must be > left value):left=%s right=%s\n",argv[0],argv[1],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }

      if(!strcmp(argv[0],"-refine")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by refine level\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	Refine =  strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || Refine < 0 || Refine > 3){
	  printf("%s value  of %s is invalid (must be between 0 and 3)\n",argv[0], argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-reff")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("-reff flag must be followed by name of file that contains -ref filenames (one per line)\n");
	  fflush(stdout);exit(1);
	}
	if(DEBUG) assert(spots_filename[num_reffiles]==NULL);

	errno = 0;
	FILE *fp = Fopen(argv[1],"r",65);
	/*	int maxcnt = 65;
	for(int cnt = 0; cnt <= maxcnt && (fp = fopen(argv[1],"r")) == NULL; cnt++){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("Cannot read file %s (to get filename list for %s),retry=%d/%d:errno=%d:%s\n",argv[1],argv[0],cnt,maxcnt,eno,err);
	  printf("\t Current Time= %s\n", currentDateTime());
	  fflush(stdout);
	  sleep(1);// sleep 1 second so this error cannot happen multiple times within the same second
	  errno = 0;
	  }*/
	if(fp == NULL)
	  exit(1);

	int orignumfiles = num_reffiles;
	int linecnt = 1;
	for(;fgets(buf,LINESIZ,fp) != NULL; linecnt++){
	  int len = strlen(buf);
	  if(len <= 0) continue;
	  if(buf[0]=='#' || buf[0]=='\n' || buf[0]=='\r')continue;
	  if(len >= min(LINESIZ,PATH_MAX)-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	    printf("Line too long or not terminated in -reff input file %s on line %d (len=%d,buf[len-1]=%c):\n%s\n",
		   argv[1],linecnt,len,buf[len-1],buf);
	    fflush(stdout);exit(1);
	  }
	  if(buf[0]=='#' || buf[0]=='\n' || buf[0]=='\r')
	    continue;
	  buf[len-1] = '\0';
	  if(num_reffiles < MAXFILES){
#ifndef WIN32
	    if(buf[0] != '/'){/* relative files names are interpreted relative to directory in which -reff file is located */
	      char filename[BUFSIZ];
	      strcpy(filename, argv[1]);
	      char *pt = strrchr(filename,'/');
	      if(pt == NULL)
		spots_filename[num_reffiles++] = strdup(buf);
	      else {
		sprintf(&pt[1],"%s",buf);
		spots_filename[num_reffiles++] = strdup(filename);
	      }
	    } else
#endif
	      spots_filename[num_reffiles++] = strdup(buf);
	  } else
	    num_reffiles++;
	}
	fclose(fp);
	if(num_reffiles > MAXFILES){
	  printf("Too many input files : increase MAXFILES=%d to at least %d\n",
		 MAXFILES, num_reffiles);
	  fflush(stdout);exit(1);
	}
	if(num_reffiles == orignumfiles){
	  printf("ERROR: input file list %s %s has no file names\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	spots_filename[num_reffiles] = NULL;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-ref")){
	if(argc<2){
	  printf("-ref flag must be followed by filename\n");
	  fflush(stdout);exit(1);
	}
	if(num_reffiles >= MAXFILES){
	  printf("Too many input files : increase MAXFILES=%d to at least %d\n",
		 MAXFILES, num_reffiles+1);
	  fflush(stdout);exit(1);
	}
	spots_filename[num_reffiles++] = strdup(argv[1]);
	spots_filename[num_reffiles] = NULL;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-refmap")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by .map filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	refmap_filename = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-res")){
	if(!colors){
	  printf("-colors must be specified before %s\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(argc < 2){
	  printf("%s must be followed by at least 1 value (up to one per color)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	int c = 0;
	for(; c < MAXCOLOR; c++){
	  if(argc <= 1+c || (c >= 1 && argv[c+1][0] == '-' && argv[c+1][1] && !isdigit(argv[c+1][1]) && argv[c+1][1] != '.'))
	    break;
	  res[c] = strtod(argv[c+1],&qt);
	  if(qt == argv[c+1] || !(*qt==0 || isspace(*qt)) || res[c] < 0.0){
	    if(colors > 1)
	      printf("The %d'th %s value %s is not valie\n",c+1,argv[0],argv[c+1]);
	    else
	      printf("The %s value %s is not valid\n",argv[0],argv[c+1]);
	    fflush(stdout);exit(1);
	  }
	}
	for(int nc = c; nc < MAXCOLOR; nc++)
	  res[nc] = res[c-1];
	argv += 1+c;
	argc -= 1+c;
	continue;
      }
      if(!strcmp(argv[0],"-resSD")){
	if(!colors){
	  printf("-colors must be specified before %s\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(argc < 2){
	  printf("%s must be followed by at least 1 value (up to one per color)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	int c = 0;
	for(; c < MAXCOLOR; c++){
	  if(argc <= 1+c || (c >= 1 && argv[c+1][0] == '-' && argv[c+1][1] && !isdigit(argv[c+1][1]) && argv[c+1][1] != '.'))
	    break;
	  resSD[c] = strtod(argv[c+1],&qt);
	  if(qt == argv[c+1] || !(*qt==0 || isspace(*qt)) || res[c] < 0.0){
	    if(colors > 1)
	      printf("The %d'th %s value %s is not valie\n",c+1,argv[0],argv[c+1]);
	    else
	      printf("The %s value %s is not valid\n",argv[0],argv[c+1]);
	    fflush(stdout);exit(1);
	  }
	}
	for(int nc = c; nc < MAXCOLOR; nc++)
	  resSD[nc] = resSD[c-1];
	argv += 1+c;
	argc -= 1+c;
	continue;
      }
      break;
    case 's':
      if(!strcmp(argv[0],"-splitFilt")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1-4 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	splitFiltMinWT = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || splitFiltMinWT > 1.0){
	  printf("%s 1st value %s is invalid : must be non-negative fraction in the range 0.0 to 1.0\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  splitFiltFN = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value %s is invalid : must be non-negative integer\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 4 && argv[3][0] != '-'){
	    splitFiltFNlen = strtod(argv[3],&qt);
	    if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	      printf("%s 3rd value %s is invalid : must be non-negative size in kb\n",argv[0],argv[3]);
	      fflush(stdout);exit(1);
	    }
	    if(argc >= 5 && argv[4][0] != '-'){
	      splitFiltFP = strtol(argv[4],&qt,10);
	      if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
		printf("%s 4th value %s is invalid : must be non-negative integer\n",argv[0],argv[4]);
		fflush(stdout);exit(1);
	      }
	      argv++;
	      argc--;
	    }

	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svRepeat")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by 1-3 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_Repeat_Tolerance = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s must be followed by non-negative number: %s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  smap_Repeat_Minele = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || smap_Repeat_Minele <= 1){
	    printf("%s 2nd value must be a positive integer at least 2: %s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }	  
	  if(argc >= 4 && argv[3][0] != '-'){
	    smap_Repeat_ConfPenalty = strtod(argv[3],&qt);
	    if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	      printf("%s 3rd value must be a non-negative value: %s\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	    }	  
	    
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInv_TrimFrac")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be following by two non-negative fractions in the range 0.0 to 1.0\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_Inv_MinTrimFrac = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_Inv_MinTrimFrac > 1.0){
	  printf("%s 1st value %s is invalid : must be non-negative fraction in the range 0.0 to 1.0\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	smap_Inv_MaxTrimFrac = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || smap_Inv_MaxTrimFrac > 1.0){
	  printf("%s 2nd value %s is invalid : must be non-negative fraction in the range 0.0 to 1.0\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-svGapOverlapRatio")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be following by non-negative fraction between 0 and 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_GapOverlapRatio = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_GapOverlapRatio > 1.0){
	  printf("%s value %s is invalid : must be non-negative fraction between 0 and 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svPrimaryOrientation")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be following by non-negative ratio\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_PrimaryOrientation = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid : must be non-negative ratio\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svTransMaxOverlap")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be following by non-negative ratio\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_TransMaxOverlap = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid : must be non-negative ratio\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svTransMaxGapOverlap")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be following by non-negative ratio\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_TransMaxGapOverlap = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid : must be non-negative ratio\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svIndelConfBySize")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be following by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_IndelConfBySize = strtol(argv[1],&qt,0);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid : must be 0 or 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svIndelRangeExpand")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s flag must be following by three non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_IndelRangeExpandRes = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid : must be non-negative value in kb\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	smap_IndelRangeExpandSep = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid : must be non-negative value in kb\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	smap_IndelRangeMinSVchange = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || smap_IndelRangeMinSVchange >= 1.0){
	  printf("%s 3rd value %s is invalid : must be non-negative value less than 1.0\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-svSmallDelConfirm")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be following by non-negative size in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_SmallDelConfirm = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid : must be non-negative size in kb\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-setMask")){
	if(argc < 2 || (argv[1][0] == '-' && strcmp(argv[1],"-1"))){
	  printf("%s flag must be followed by non-negative number OR -1\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	setMask = strtol(argv[1],&qt,0);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid : must be non-negative number\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-splitFragile")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s flag must be followed by 3 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	SplitFragileD = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid : must be non-negative number\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	SplitFragileR = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || SplitFragileR <= 1.0){
	  printf("%s 3rd value %s is invalid : must be positive number > 1.0\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	SplitFragileE = strtol(argv[3],&qt,10);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || SplitFragileE > 1){
	  printf("%s 4th value %s is invalid : must be 0 or 1\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-splitWT")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number <= 1.0\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	splitWT = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || splitWT > 1.0){
	  printf("%s 2nd value %s is invalid : must be non-negative number <= 1.0\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-skip_dist")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	skip_dist = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  skip_dist_Ext = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value %s is invalid\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-skip_score")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	skip_score = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && argv[2][0] != '-'){
	  skip_score_Ext = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 3rd value %s is invalid\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-swap")){
	CmapSwap = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-svfile")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	SVfile = argv[1];
	if(VERB>=2){
	  printf("SVfile=%s\n",SVfile);
	  fflush(stdout);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svcompare")){
	if(argc < 5 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by 2 filenames and 2 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	svcompare1 = strdup(argv[1]);
	svcompare2 = strdup(argv[2]);
	SVcompareKB = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt))){	
	  printf("%s 3rd value is invalid:%s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	SVcompareFrac = strtod(argv[4],&qt);
	if(qt == argv[4] || !(*qt==0 || isspace(*qt))){	
	  printf("%s 4th value is invalid:%s\n",argv[0],argv[4]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 6 && argv[5][0] != '-'){
	  SVConf = strtod(argv[5],&qt);
	  if(qt == argv[5] || !(*qt==0 || isspace(*qt)) || SVConf < 0.0){
	    printf("%s 5th value is invalid (must be non-negative):%s\n",argv[0],argv[5]);
	    fflush(stdout);exit(1);
	  }

	  if(argc >= 7 && argv[6][0] != '-'){
	    SVDeltaMatch = strtod(argv[6],&qt);
	    if(qt == argv[6] || !(*qt==0 || isspace(*qt)) || SVDeltaMatch < 0.0){
	      printf("%s 6th value is invalid (must be non-negative):%s\n",argv[0],argv[6]);
	      fflush(stdout);exit(1);
	    }

	    if(argc >= 8 && argv[7][0] != '-'){
	      SVSiteMatch = strtol(argv[7],&qt,10);
	      if(qt == argv[7] || !(*qt==0 || isspace(*qt)) || SVSiteMatch < 0){
		printf("%s 7th value is invalid (must be non-negative integer):%s\n",argv[0],argv[7]);
		fflush(stdout);exit(1);
	      }

	      if(argc >= 9 && argv[8][0] != '-'){
		SVXmapOverlap = strtod(argv[8],&qt);
		if(qt == argv[8] || !(*qt==0 || isspace(*qt)) || SVXmapOverlap < 0.0){
		  printf("%s 8th value is invalid (must be non-negative):%s\n", argv[0], argv[8]);
		  fflush(stdout);exit(1);
		}

		argc--;
		argv++;
	      }

	      argc--;
	      argv++;
	    }

	    argc--;
	    argv++;
	  }

	  argc--;
	  argv++;
	}

	argc -= 5;
	argv += 5;

	if(argc > 0){
	  printf("ERROR: -svcompare must be last option on commandline, %d subsequent args:", argc);
	  for(int t = 0; t < argc; t++)
	    printf(" %s",argv[t]);
	  printf("\n");
	  fflush(stdout);
	  fflush(stdout);exit(1);
	}

	if(!stdout_file)
	  printversion(stdout);
	indel_compare(svcompare1,svcompare2,SVcompareKB,SVcompareFrac);
	exit(0);
      }
      if(!strcmp(argv[0],"-svcheck")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	svcheck = argv[1];
	if(argc >= 3 && argv[2][0] != '-'){
	  svcheck_confidence = 0.0;
	  svcheckKB = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || svcheckKB < 0.0){
	    printf("2nd %s value is invalid:%s (must be >= 0)\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  if(argc >= 4 && argv[3][0] != '-'){
	    svcheckMerge = strtod(argv[3],&qt);
	    if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || svcheckMerge < 0.0){
	      printf("3rd %s value is invalid:%s (must be >= 0)\n", argv[0], argv[3]);
	      fflush(stdout);exit(1);
	    }
	    if(argc >= 5 && !(argv[4][0] == '-' && !isdigit(argv[4][1]))){
	      svcheck_confidence = strtod(argv[4],&qt);
	      if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
		printf("4th %s value is invalid:%s\n",argv[0],argv[4]);
		fflush(stdout);exit(1);
	      }
	      argv++;
	      argc--;
	    }
	    argv++;
	    argc--;
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-subsetmaps")){
	SubsetMaps = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-splitcnt")){
	SplitCnt = -1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-splitsite")){
	SplitSite = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-splitrev")){
	if(argc < 2 || argv[1][0] == '-'){
	  SplitRev = 1;
	  argv++;
	  argc--;
	} else {
	  SplitRev = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || SplitRev < 0){
	    printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv += 2;
	  argc -= 2;
	}
	continue;
      }
      if(!strcmp(argv[0],"-stdout")){
	char filename[PATH_MAX];
	if(argc < 2 || argv[1][0] == '-'){/* 0 args */
	  if(!output_prefix){
	    printf("-o option must be specified before %s with no args\n",argv[0]);
	    fflush(stdout);exit(1);
	  }
	  if(!CMapID)
	    sprintf(filename,"%s.stdout", output_prefix);
	  else
	    sprintf(filename,"%s_id%lld.stdout", output_prefix, CMapID);
	  stdout_file = strdup(filename);
	  argv++;
	  argc--;
	} else { /* 1 arg */
	  stdout_file = strdup(argv[1]);
	  argv += 2;
	  argc -= 2;
	}
	
	printf("Appending stdout to %s\n",stdout_file);
	fflush(stdout);

	//	if(!(stderr_file && !strcmp(stderr_file,stdout_file)))
	//	  checkFile(stdout_file);

	if(freopen(stdout_file,"a",stdout)==NULL){
	  int myerrno = errno;
	  fprintf(stderr,"Unable to redirect stdout to %s: %s\n",stdout_file, strerror(myerrno));
	  fflush(stderr);exit(1);
	}
	if(!(stderr_file && !strcmp(stderr_file,stdout_file))){
	  printversion(stdout);

	  time_t my_time = time(NULL); 
	  printf("START TIME: %s", ctime(&my_time));   // ctime() used to give the present time 
	  fflush(stdout);
	}

	continue;
      }
      if(!strcmp(argv[0],"-stderr")){
	char filename[PATH_MAX];
	if(argc < 2 || argv[1][0] == '-'){/* 0 args */
	  if(!output_prefix){
	    printf("-o option must be specified before %s with no args\n",argv[0]);
	    fflush(stdout);exit(1);
	  }
	  if(!CMapID)
	    sprintf(filename,"%s.stdout", output_prefix);
	  else
	    sprintf(filename,"%s_id%lld.stdout", output_prefix, CMapID);
	  stderr_file = strdup(filename);
	  argv++;
	  argc--;
	} else { /* 1 args */
	  stderr_file = strdup(argv[1]);
	  argv += 2;
	  argc -= 2;
	}
	fprintf(stderr,"Appending stderr to %s\n",stderr_file);
	fflush(stderr);

	//	if(!(stdout_file && !strcmp(stderr_file,stdout_file)))
	//	  checkFile(stderr_file);

	if(freopen(stderr_file,"a",stderr)==NULL){
	  int myerrno = errno;
	  printf("Unable to redirect stderr to %s: %s\n",stderr_file, strerror(myerrno));
	  fflush(stdout);exit(1);
	}
	if(!(stdout_file && !strcmp(stderr_file,stdout_file))){
	  printversion(stderr);

	  time_t my_time = time(NULL); 
	  fprintf(stderr,"START TIME: %s", ctime(&my_time));   // ctime() used to give the present time 
	  fflush(stderr);
	}

	continue;
      }
      if(!strcmp(argv[0],"-stitch")){
	MapStitched = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-stitchLoc")){
	MapStitchLoc = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-skipidf")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by name of file containing molecule ids (one per line)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	FILE *fp;
	if((fp = fopen(argv[1],"r"))==NULL){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("Cannot read file %s:errno=%d:%s\n",argv[1],eno,err);
	  fflush(stdout);exit(1);
	}
	int SkipIDmax = 1000;
	SkipID = new long long[SkipIDmax];
	SkipIDcnt = 0;
	int linecnt = 1;
	for(;fgets(buf,BUFSIZ,fp) != NULL; linecnt++){
	  int len = strlen(buf);
	  if(len <= 0) continue;
	  if(buf[0]=='#' || buf[0] == '\n' || buf[0]== '\r') continue;
	  if(len >= BUFSIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	    printf("Line too long or not terminated in input file %s on line %d:\n%s\n",
		   argv[1],linecnt,buf);
	    fflush(stdout);exit(1);
	  }
	  buf[--len] = '\0';
	  while(len > 0 && isspace(buf[len-1]))
	    buf[--len] = '\0';
	  if(SkipIDcnt >= SkipIDmax){
	    SkipIDmax *= 2;
	    long long *newSkipID = new long long[SkipIDmax];
	    memcpy(newSkipID,SkipID,SkipIDcnt * sizeof(long long));
	    delete [] SkipID;
	    SkipID = newSkipID;
	  }
	  SkipID[SkipIDcnt] = strtoll(buf,&qt,10);
	  if(qt == buf || !(*qt==0 || isspace(*qt) || !strcmp(qt,"LL")) || SkipID[SkipIDcnt] < 0){
	    printf("-selectid value %s on line %d of %s has invalid long long integer syntax(qt=%s)\n",buf,linecnt,argv[1],qt);
	    fflush(stdout);exit(1);
	  }
	  SkipIDcnt++;
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-selectidf")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by name of file containing molecule ids (one per line)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(VERB>=2){
	  printf("Checking %s for molecule ids (one per line)\n",argv[1]);
	  fflush(stdout);
	}
	FILE *fp;
	if((fp = fopen(argv[1],"r"))==NULL){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("Cannot read file %s:errno=%d:%s\n",argv[1],eno,err);
	  fflush(stdout);exit(1);
	}
	SelectIDcnt = max(SelectIDcnt, 0LL);
	int linecnt = 1;
	for(;fgets(buf,BUFSIZ,fp) != NULL; linecnt++){
	  int len = strlen(buf);
	  if(VERB>=2){
            printf("%s Line %d (SelectIDcnt=%d,len=%d):%s\n\n",argv[1],linecnt,SelectIDcnt,len,buf);
	    fflush(stdout);
	  }
	  if(len <= 0) continue;
	  if(buf[0]=='#' || buf[0] == '\n' || buf[0]== '\r') continue;
	  if(len >= BUFSIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	    printf("Line too long or not terminated in input file %s on line %d:\n%s\n",
		   argv[1],linecnt,buf);
	    fflush(stdout);exit(1);
	  }
	  buf[--len] = '\0';
	  while(len > 0 && isspace(buf[len-1]))
	    buf[--len] = '\0';
	  if(SelectIDcnt >= SelectIDmax){
	    SelectIDmax = max(1024,SelectIDmax*2);
	    long long *newSelectID = new long long[SelectIDmax];
	    if(SelectIDcnt)
	      memcpy(newSelectID,SelectID,SelectIDcnt * sizeof(long long));
	    delete [] SelectID;
	    SelectID = newSelectID;
	  }
	  SelectID[SelectIDcnt] = strtoll(buf,&qt,10);
	  if(VERB>=2){
	    printf("Line %d: len=%d, SelectIDcnt=%d, strtoll= %lld, qt=%s,buf=%s\n",linecnt,len,SelectIDcnt,SelectID[SelectIDcnt],qt,buf);
	    fflush(stdout);
	  }
	  if(qt == buf ||  !(*qt==0 || (isspace(*qt) && *qt != '\r') || !strcmp(qt,"LL")) || SelectID[SelectIDcnt] < 0){
	    printf("-selectid value %s on line %d of %s has invalid long long integer syntax(qt=%s)\n",buf,linecnt,argv[1],qt);
	    fflush(stdout);exit(1);
	  }
	  SelectIDcnt++;
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-selectScanId")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by the RundId range (2 positive integers)\n",argv[0]);
	  fflush(stdout);
	}
	SelectScanId = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || SelectScanId <= 0){
	  printf("%s 1st value %s is invalid (must be positive integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	SelectScanId2 = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || SelectScanId2 < SelectScanId){
	  printf("%s 2nd value %s is invalid (must be positive integer >= 1st value %s)\n",argv[0],argv[2],argv[1]);
	  fflush(stdout);exit(1);
	}

	SelectScanId--;/* Cmap.h UniqueScanId is 0-based, while value in BNX is 1-based */
	SelectScanId2--;/* Cmap.h UniqueScanId is 0-based, while value in BNX is 1-based */

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-selectid")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a list of +ve map IDs\n",argv[0]);
	  fflush(stdout);
	}
	SelectIDcnt = 0;
	// NEW	int origSelectIDcnt = SelectIDcnt = max(SelectIDcnt, 0LL);
	while(1+SelectIDcnt < argc && argv[1+SelectIDcnt][0] != '-')
	  SelectIDcnt++;
	if(SelectIDcnt >= SelectIDmax){
	  SelectIDmax = max(SelectIDcnt,SelectIDmax*2 + 1024);
	  long long *newSelectID = new long long[SelectIDmax];
	  // NEW	  if(origSelectIDcnt)  memcpy(newSelectID,SelectID,origSelectIDcnt * sizeof(long long));
	  delete [] SelectID;
	  SelectID = newSelectID;
	}
	for(int i = 0; i < SelectIDcnt/* NEW WAS SelectIDcnt - origSelectIDcnt */; i++){
	  SelectID[i /* NEW WAS origSelectIDcnt + i */] = strtoll(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt) || !strcmp(qt,"LL")) || SelectID[i] < 0){
	    printf("-selectid value %s has invalid long long integer syntax (qt=%s)\n",argv[1],qt);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-subset")){
	if(argc < 3){
	  printf("%s flag must be followed by two to four numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CmapStart = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || CmapStart <= 0){
	  printf("1st %s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	CmapEnd = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || CmapEnd <= 0.0){
	  printf("2nd %s value is invalid:%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	CmapFrac = 0.0;
	perCohort = 0;

	if(argc >= 4 && argv[3][0] != '-'){
	  CmapFrac = strtod(argv[3],&qt);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || CmapFrac > 1.0){
	    printf("3rd %s value is invalid (must be non-negative fraction) :%s\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  
	  if(argc >= 5 && argv[4][0] != '-'){
	    perCohort = strtol(argv[4],&qt,10);
	    if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	      printf("4th %s value is invalid (must be non-negative integer) :%s\n",argv[0],argv[4]);
	      fflush(stdout);exit(1);
	    }	    

	    argv++;
	    argc--;
	  }

	  argv++;
	  argc--;
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-subsetGB")){
	if(argc < 3 || argv[1][0] == '-' || argv[2][0] == '-'){
	  printf("%s flag must be followed by two negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CmapGB = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || CmapGB <= 0){
	  printf("1st %s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	CmapMinLen = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || CmapMinLen < 0.0){
	  printf("2nd %s value is invalid:%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	CmapGBrnd = 1;
	if(argc >= 3 && argv[3][0] != '-'){
	  CmapGBrnd = strtol(argv[3],&qt,10);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || CmapMinLen < 0.0){
	    printf("3rd %s value is invalid:%s\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}

	CmapStart = -1;
	CmapEnd = -1.0;
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-subsetbin")){
	if(argc < 3){
	  printf("%s flag must be followed by two integers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CmapBin = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || CmapBin < 0){
	  printf("1st %s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	CmapNumBins = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || CmapNumBins < CmapBin){
	  printf("2nd %s value of %s is invalid (must be at least as large as 1st value %s)\n",argv[0],argv[2],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-split")){
	CmapSplit = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-sort-sizeinc")){
	CmapSort = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-sort-sizedec")){
	CmapSort = 2;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-sort-sitesinc")){
	CmapSort = 3;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-sort-sitesdec")){
	CmapSort = 4;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-sort-idinc")){
	CmapSort = 5;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-sort-runindexinc")){
	CmapSort = 7;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-sort-runindexShuffle")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by min length in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CmapSort = 8;
	CmapSortMinLen = strtod(argv[1],&qt); 
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || CmapSortMinLen <= 0.0){
	  printf("%s value %s is invalid : must be +ve length in kb\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	RandSeed = -1;/* do NOT randomize by default */

	if(argc >= 3 && argv[2][0] != '-'){
	  RandSeed = strtol(argv[2],&qt,10);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || RandSeed <= 0){
	    printf("%s Random Seed value of %s is invalid (must be > 0)\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  argv++;
	  argc--;
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-s") || !strcmp(argv[0],"-sd")){
	if(!colors){
	  printf("-colors must be specified before %s\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(argc < 2){
	  printf("%s must be followed by at least 1 value (up to one per color)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	int i = 0;
	for(; i < MAXCOLOR; i++){
	  if(argc <= 1+i || (i >= 1 && argv[i+1][0] == '-' && argv[i+1][1] && !isdigit(argv[i+1][1]) && argv[i+1][1] != '.'))
	    break;
	  SD[i] = strtod(argv[i+1],&qt);
	  if(argv[i+1]==qt || !(*qt==0 || isspace(*qt))){
	    if(colors > 1)
	      printf("The %d. %s value %s is not valid\n",i+1,argv[0],argv[i+1]);
	    else
	      printf("The %s value %s is not valid\n",argv[0],argv[i+1]);
	    fflush(stdout);exit(1);
	  }
	}
	for(int nc = i; nc < MAXCOLOR; nc++)
	  SD[nc] = SD[i-1];
	argv += 1+i;
	argc -= 1+i;

	continue;
      }
      if(!strcmp(argv[0],"-sf")){
	if(!colors){
	  printf("-colors must be specified before %s\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(argc < 2){
	  printf("%s must be followed by at least 1 value (up to one per color)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	int i = 0;
	for(; i < MAXCOLOR; i++){
	  if(argc <= 1+i || (i >= 1 && argv[i+1][0] == '-' && argv[i+1][1] && !isdigit(argv[i+1][1]) && argv[i+1][1] != '.'))
	    break;
	  SF[i] = strtod(argv[i+1],&qt);
	  if(argv[i+1]==qt || !(*qt==0 || isspace(*qt)) || SF[i] < 0.0){
	    if(colors > 1)
	      printf("The %d. %s value %s is not valid\n", i+1,argv[0],argv[i+1]);
	    else
	      printf("The %s value %s is not valid\n", argv[0],argv[i+1]);
	    fflush(stdout);exit(1);
	  }
	}
	for(int nc = i; nc < MAXCOLOR; nc++)
	  SF[nc] = SF[i-1];
	argv += 1+i;
	argc -= 1+i;

	continue;
      }
      if(!strcmp(argv[0],"-sr")){
	if(!colors){
	  printf("-colors must be specified before %s\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(argc < 2){
	  printf("%s must be followed by at least 1 value (up to one per color)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	int i = 0;
	for(; i < MAXCOLOR; i++){
	  if(argc <= 1+i || (i >= 1 && argv[i+1][0] == '-' && argv[i+1][1] && !isdigit(argv[i+1][1]) && argv[i+1][1] != '.'))
	    break;
	  SR[i] = strtod(argv[i+1],&qt);
	  if(argv[i+1]==qt || !(*qt==0 || isspace(*qt)) || SR[i] < 0.0){
	    if(colors > 1)
	      printf("The %d. %s value %s is not valid\n", i+1,argv[0],argv[i+1]);
	    else
	      printf("The %s value %s is not valid\n", argv[0],argv[i+1]);
	    fflush(stdout);exit(1);
	  }
	}
	for(int nc = i; nc < MAXCOLOR; nc++)
	  SR[nc] = SR[i-1];
	if(VERB>=2){
	  for(int t = 0; t < 2; t++)
	    printf("SR[%d]=%0.6f\n",t,SR[t]);
	  fflush(stdout);
	}
	argv += 1+i;
	argc -= 1+i;

	continue;
      }
      if(!strcmp(argv[0],"-se")){
	if(!colors){
	  printf("-colors must be specified before %s\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	if(argc < 2){
	  printf("%s must be followed by at least 1 value (up to one per color)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	int i = 0;
	for(; i < MAXCOLOR; i++){
	  if(argc <= 1+i || (i >= 1 && argv[i+1][0] == '-' && argv[i+1][1] && !isdigit(argv[i+1][1]) && argv[i+1][1] != '.'))
	    break;
	  SE[i] = strtod(argv[i+1],&qt);
	  if(argv[i+1]==qt || !(*qt==0 || isspace(*qt)) || SE[i] < 0.0){
	    if(colors > 1)
	      printf("The %d. %s value %s is not valid\n", i+1,argv[0],argv[i+1]);
	    else
	      printf("The %s value %s is not valid\n", argv[0],argv[i+1]);
	    fflush(stdout);exit(1);
	  }
	}
	for(int nc = i; nc < MAXCOLOR; nc++)
	  SE[nc] = SE[i-1];
	argv += 1+i;
	argc -= 1+i;

	continue;
      }
      if(!strcmp(argv[0],"-sv")){
	if(argc >= 2 && argv[1][0] != '-'){ /* this is optional */
	  svlvl = strtol(argv[1],&qt,10);
	  if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || svlvl < 0 || svlvl > 3){
	    printf("-sv arg must be 0, 1, 2, or 3--defaulting to 1.\n");
	    svlvl = 1;
	  }
	  argv++;
	  argc--;
	}
	else {
	  svlvl = 1;
	}
	dosvdetect = (svlvl > 0) ? 1 : 0; /* see parameters.cpp/h */
	argv++;
	argc--;
	continue;
      }  
      if(!strcmp(argv[0],"-svDup_maxQryGapMult")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_Dup_maxQryGapMult = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svDup_maxSpacing")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_Dup_maxSpacing = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative number)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInv_minRefSites")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_Inv_minRefSites = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInv_minQrySites")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_Inv_minQrySites = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }

      if(!strcmp(argv[0],"-svInvDup_maxSpacing")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_InvDup_maxSpacing = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svDup_minRefOverlapSites")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_Dup_minRefOverlapSites = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInvDup_maxQryGapMult")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_InvDup_maxQryGapMult = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInvDup_minRefOverlapFrac")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_InvDup_minRefOverlapFrac = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInv_maxRefOverlapFrac")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_Inv_maxRefOverlapFrac = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInv_maxQryOverlapFrac")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_Inv_maxQryOverlapFrac = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInvDup_maxSpacing")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_InvDup_maxSpacing = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svDup_maxSpacing")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_Dup_maxSpacing = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svDupSplit")){
	if(argc < 6 || argv[1][0]== '-' || argv[2][0] == '-' || argv[3][0] == '-' || argv[4][0] == '-' || argv[5][0] == '-'){
	  printf("%s must be followed by 5 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	smap_DupSplit_maxoverlap = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	smap_DupSplit_maxgap = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid (must be non-negative integer)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	smap_DupSplit_maxend = strtol(argv[3],&qt,10);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt))){
	  printf("%s 3rd value %s is invalid (must be non-negative integer)\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	smap_DupSplit_MinQryLen = strtod(argv[4],&qt);
	if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	  printf("%s 4th value %s is invalid (must be non-negative length in kb)\n",argv[0],argv[4]);
	  fflush(stdout);exit(1);
	}

	smap_DupSplit_MinQrySites = strtol(argv[5],&qt,10);
	if(qt == argv[5] || !(*qt==0 || isspace(*qt))){
	  printf("%s 5th value %s is invalid (must be non-negative integer)\n", argv[0], argv[5]);
	  fflush(stdout);exit(1);
	}

	if (argc >= 7 && argv[6][0] != '-'){
	  smap_DupSplit_singleEnd = strtol(argv[6],&qt,10);
	  if(qt == argv[6] || !(*qt==0 || isspace(*qt)) || smap_DupSplit_singleEnd > 3){
	    printf("%s 6th value %s is invalid (must be integer from 0 to 3)\n",argv[0],argv[6]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}

	argv += 6;
	argc -= 6;
	continue;
      }
      if(!strcmp(argv[0],"-svZygosity")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_zygosity = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_zygosity > 1){
	  printf("%s value %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svTransOrient")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smapTransOrient = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smapTransOrient > 1){
	  printf("%s value %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svSize")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smapSize = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smapSize > 1){
	  printf("%s value %s is invalid (must be 0 or 1)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svFreq")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by 0,1 or 2\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smapFreq = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smapFreq > 2){
	  printf("%s value %s is invalid (must be 0,1 or 2)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svMinConf")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by confidence\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_min_conf = strtod(argv[1],&qt);
	//require between 0 and 100 becuase larger than this is quite rare
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_min_conf < 0.0 || smap_min_conf > 100){
	  printf("%s value %s is invalid (must be between 0 and 100)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svAlignConfByLen")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by conf / length (kb)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_confbylen = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_confbylen < 0.0){
	  printf("%s value %s is invalid (must be non-negative)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svAlignConfByInterval")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by conf per aligned Label Interval\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_ConfByInterval = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_ConfByInterval < 0.0){
	  printf("%s value %s is invalid (must be non-negative)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svConfFile")){ 
	conf_filename = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svConfMinsize")){
	if(argc < 3 || argv[1][0]=='-'){
	  printf("%s must be followed by minimum length (kb) and PPV value\n",argv[0]);
	  fflush(stdout);
	}
	svConfMinsize = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || svConfMinsize < 0.0){
	  printf("%s 1st value %s is invalid (must be non-negative size in kb)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	svConfNAvalue = strtod(argv[1],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid (must be a number)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-svSideMinLen")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by length (kb)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_side_alignkb = strtod(argv[1],&qt);
	//require between 0 and 1e4 (kb) becuase larger than this doesn't make sense
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_side_alignkb < 0.0 || smap_side_alignkb > 1e4){
	  printf("%s value %s is invalid (must be between 0 and 1e4)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svMaxGap")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by 1 or 2 non-negative sizes in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_sv_sizekb = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s must be followed by non-negative length in kb:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	smap_inv_MaxGap = 0.0;
	if(argc >= 3 && argv[2][0] !='-'){
	  smap_inv_MaxGap = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value is invalid : must be non-negative length in kb: %s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }

	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svIndelMaxoverlap")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by non-negative length in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	sv_indel_maxoverlap = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be be non-negative size in kb)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svDelMaxQryOverlap")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by non-negative length in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_del_maxqryoverlap = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be be non-negative size in kb)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svDelBigQryOverlap")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by non-negative length in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_del_bigqryoverlap = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be be non-negative size in kb)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svDelMaxGapOverlap")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by non-negative fraction <= 1.0\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	svDelMaxGapOverlap = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || svDelMaxGapOverlap > 1.0){
	  printf("%s value %s is invalid (must be be non-negative fraction <= 1.0)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInsMaxGapOverlap")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by non-negative fraction <= 1.0\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	svInsMaxGapOverlap = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || svInsMaxGapOverlap > 1.0){
	  printf("%s value %s is invalid (must be be non-negative fraction <= 1.0)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svEndMinLen")){
	if(argc<2){
	  printf("%s must be followed by length (kb)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_min_end_lenkb = strtod(argv[1],&qt);
	//require between 0 and 1e4 (kb) becuase larger than this doesn't make sense
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_min_end_lenkb < 0.0 || smap_min_end_lenkb > 1e4){
	  printf("%s value %s is invalid (must be between 0 and 1e4)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svEndMinNLabels")){
	if(argc<2){
	  printf("%s must be followed by length (kb)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_min_end_nlabels = strtol(argv[1],&qt,10);
	//require between 1 and 1000 becuase 0 labels unaligned doesn't make sense
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_min_end_nlabels < 1 || smap_min_end_nlabels > 1000){
	  printf("%s value %s is invalid (must be between 1 and 1000)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svIndelTinySize")){
	if(argc<2){
	  printf("%s must be followed by length (kb)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_indel_tiny_maxsize = strtod(argv[1],&qt);
	//require between 0 and 100 becuase < 0 is impossible and > 100kb isn't really very tiny
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_indel_tiny_maxsize < 0 || smap_indel_tiny_maxsize > 100){
	  printf("%s value %s is invalid (must be between 0 and 100)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svOverlap")){
	if(argc < 2){
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	//	bool dosv_simpleoverlap = bool(strtol(argv[1],&qt,10));
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be integer; non-zero value interpreted as true)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	printf("WARNING: -svOverlap is obsolete and will be ignored. See -OverlapFilter instead\n");
	fflush(stdout);

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svTransMask")){ 
	mask_filename = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svTransMaxGap")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by length (kb)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_translocation_maxsep = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_translocation_maxsep < 0.0 || smap_translocation_maxsep > 1e3){
	  printf("%s value %s is invalid (must be between 0 and 1e3)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svTransMaxOverlapFrac")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by fraction\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_translocation_ovlfr = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || smap_translocation_ovlfr < 0.0 || smap_translocation_ovlfr > 0.5){
	  printf("%s value %s is invalid (must be between 0 and 0.5)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svTransMaxOverlapSize")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by length (kb)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_translocation_ovlmx = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInvMaxOverlapSize")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by length (kb)\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_InvMaxOverlapSize = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid (must be non-negative)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-svInvConfMinInterval")){
	if(argc < 3 || argv[1][0]=='-' || argv[2][0]=='-'){
	  printf("%s must be followed by 2 non-negative numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	smap_Inv_ConfMinIntervalSize = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 1st value %s is invalid (must be non-negative size in kb)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	smap_Inv_ConfMinNumIntervals = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value %s is invalid (must be non-negative integer)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}

	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-scaffold")){ 
	scaf_filename = argv[1];
	doscaffold = 1;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-simpleRepeat")){
	simpleRepeat = true;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-simpleRepeatStandalone")){
	simpleRepeatStandalone = true;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-simpleRepeatMinEle")){
	if(argc<2){
	  printf("%s must be followed by integer value\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	simpleRepeatMinEle = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || simpleRepeatMinEle < 2){
	  printf("%s value %s is invalid (must be integer >= 2)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }      
      if(!strcmp(argv[0],"-simpleRepeatFilter")){
	if(argc<2){
	  printf("%s must be followed by integer value\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	simpleRepeatFilter = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || simpleRepeatFilter < 0 || simpleRepeatFilter > 3){
	  printf("%s value %s is invalid (must be integer from 0 to 3)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }      
      if(!strcmp(argv[0],"-simpleRepeatTolerance")){
	if(argc<2){
	  printf("%s must be followed by integer value\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	simpleRepeatTolerance = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || simpleRepeatTolerance < 0 || simpleRepeatTolerance > .5){
	  printf("%s value %s is invalid (must be between 0 to 0.5)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }      
      break;
    case 't':
      if(!strcmp(argv[0],"-thetaScale")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a fraction in the range 0 .. 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	thetaScale = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || thetaScale > 1.0){
	  printf("%s value is invalid (must be in the range 0 .. 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-thetaScaleRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a fraction in the range 0 .. 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	thetaScaleRef = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || thetaScaleRef > 1.0){
	  printf("%s value is invalid (must be in the range 0 .. 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'u':
      if(!strcmp(argv[0],"-usecolor")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	usecolor = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || usecolor <= 0 || usecolor > MAXCOLOR){
	  printf("%s value is invalid (must be between 1 and %d):%s\n",argv[0],MAXCOLOR,argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-unmapped")){
	if(argc<2){
	  printf("-ref flag must be followed by output filename prefix\n");
	  fflush(stdout);exit(1);
	}
	UnMappedPrefix = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'w':
      if(!strcmp(argv[0],"-winbreak")){
	if(argc < 4 || argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-'){
	  printf("%s flag must be followed by 3 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	WinBreak = 1;

	WinBreakLen = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || WinBreakLen <= 0.0){
	  printf("%s 1st value is invalid (must be +ve number):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	
	WinBreakShift = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt == 0 || isspace(*qt)) || WinBreakLen <= 0.0){
	  printf("%s 2nd value is invalid (must be +ve number):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	
	WinBreak = strtol(argv[3],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || WinBreak <= 0){
	  printf("%s 3rd value is invalid (must be +ve integer):%s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}

	argv += 4;
	argc -= 4;
	continue;
      }
      break;

    case 'x':
      if(!strcmp(argv[0],"-xmapchim")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	XmapChim = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || XmapChim < 0){
	  printf("%s value %s is invalid (must be non-negative integer)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 3 && !(argv[2][0] == '-' && argv[2][1] && !isdigit(argv[2][1]))){
	  XmapChimSpacing = strtod(argv[2],&qt);
	  if(qt==argv[2] || !(*qt==0 || isspace(*qt))){
	    printf("%s 2nd value %s is invalid (must be non-negative value)\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-xmapNoQueryOverlap")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by non-negative length in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	XmapNoQueryOverlap = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-xmapFilter")){
	if(argc < 4 || argv[1][0]=='-' || argv[2][0]=='-' || argv[3][0]=='-'){
	  printf("%s must be followed by 3 positive numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	double Pvalue1 = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || Pvalue1 <= 0.0 || Pvalue1 > 1.0){
	  printf("%s 1st value %s is invalid : must lie between 0 and 1\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	double Pvalue2 = strtod(argv[2],&qt);
	if(qt==argv[2] || !(*qt==0 || isspace(*qt)) || Pvalue2 <= Pvalue1 || Pvalue2 > 1.0){
	  printf("%s 2nd value %s is invalid : must lie between 1st value and 1\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	double Pvalue3 = strtod(argv[3],&qt);
	if(qt==argv[3] || !(*qt==0 || isspace(*qt)) || Pvalue3 <= Pvalue2 || Pvalue3 > 1.0){
	  printf("%s 3rd value %s is invalid : must lie between 2nd value and 1\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	double Ilog10 = 1.0/log(10.0);
	XmapFilterDelta = -log(Pvalue3)*Ilog10;
	XmapFilterThresh1 = -log(Pvalue1)*Ilog10;
	XmapFilterThresh2 = -log(Pvalue2)*Ilog10;

	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-xmapUnique")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s must be followed by non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	XmapUnique = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value %s is invalid\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-xmaplen")){
	XmapLen = 1;
	argv++;
	argc--;
	continue;
      }
      break;
    default:
      break;
    }

    printf("unknown option:'%s' (argc = %d)\n",argv[0],argc);
    if(argc >= 3)
      printf("     Next 3 options are '%s' '%s' '%s'\n",argv[1],argv[2],argv[3]);
    else if(argc >= 2)
      printf("     Next 2 options are '%s' '%s'\n",argv[1],argv[2]);
    else if(argc >= 1)
      printf("     Next option is '%s'\n",argv[1]);
    printf("Use -help for list of options\n");
    fflush(stdout);exit(1);
  }

  if(!stdout_file){
    printversion(stdout);
    time_t my_time = time(NULL); 
    printf("START TIME: %s", ctime(&my_time));   // ctime() used to give the present time 
    fflush(stdout);
  }

  if(VERB){
    char *OMP_WAIT_POLICY = getenv("OMP_WAIT_POLICY");
    printf("OMP_WAIT_POLICY = %s\n",OMP_WAIT_POLICY);
    fflush(stdout);
  }

  if(KCcnt > 0){
    if(KCcnt > KBcnt){
      printf("-extsplitOutlierC: Ignoring last %d values %0.2f .. %0.2f\n",KCcnt-KBcnt,extsplitOutlierC[KBcnt],extsplitOutlierC[KCcnt-1]);
      fflush(stdout);

      KCcnt = KBcnt;
    } 
    if(KCcnt >= 2 && KCcnt < KBcnt){
      printf("-extsplitOutlierC: Padding %d additional values with %0.2f\n",KBcnt-KCcnt,extsplitOutlierC[KCcnt-1]);
      fflush(stdout);

      for(int B = KCcnt; B < KBcnt; B++)
	extsplitOutlierC[B] = extsplitOutlierC[KCcnt-1];

      KCcnt = KBcnt;
    }
  }

  if(cresMaxLPdrop > 0.0 && HapSitePvalue <= 0.0 && !cresFix){
    printf("WARNING : -cresMaxLPdrop %0.2f has no effect without -cresCheckDel or -Haplotype\n",  cresMaxLPdrop);
    fflush(stdout);
  }

  if(ShiftOffset1 > 31 || ShiftOffset1 <= 0){
    printf("Disabling -hashGrouped %d %d since first arg %d is not between 1 and 30\n",ShiftOffset1,GroupedScore,ShiftOffset1);
    fflush(stdout);
    ShiftOffset1 = 31;
  }
  if(GroupedScore <= 0){
    printf("Disabling -hashGrouped %d %d since 2nd arg is 0\n",ShiftOffset1,GroupedScore);
    fflush(stdout);
    ShiftOffset1 = 31;
  }
  if(ShiftOffset1 == 31)
    GroupedScore = HashScore;
  if(ShiftOffset1 <= 30 && GroupedScore > HashScore){
    printf("-hashGrouped %d %d : 2nd arg %d is larger than -hashgen score value of %d : reducing to -hashGrouped %d %d\n",ShiftOffset1,GroupedScore,GroupedScore,HashScore,ShiftOffset1,HashScore);
    fflush(stdout);
    GroupedScore = HashScore;
  }

  if(NoStat < 0)
    NoStat = (Refine ? 1 : 0);

  if(MapScale && !NoStat && ScanCorrection){
    printf("WARNING: -ScanScaling output <prefix>_rescaled.bnx will not reflect -MapScale, which must be applied again each time this BNX file is used\n");
    printf("         -MapScale can only be applied to maps that align to reference\n");
    fflush(stdout);
  }
  if(MapScale && Refine){
    printf("WARNING: -MapScale with -Refine will apply molecule rescaling a 2nd time if input BNX is based on -mapped output with -MapScale\n");
    fflush(stdout);
  }

  if(BestRefWT && BestRef){
    printf("WARNING:Cannot use both -BestRef and -BestRefWT : will ignore -BestRef\n");
    fflush(stdout);
    BestRef = 0;
  }
  if(BestRefWT && NoSplit != 2){
    printf("WARNING:-BestRefWT requires -nosplit 2 : changing to -nosplit 2\n");
    fflush(stdout);
    NoSplit = 2;
  }

  if(ChimQuality && !(TrimNorm >= 0 && COVERAGE_TRIM > 0 && COVERAGE_TRIM_LEN > 0.0)){
    printf("-ChimQuality requires -TrimNorm N -CovTrim M -CovTrimLen L (N >= 0, M > 0, L > 0)\n");
    printf("   Recommended values : -TrimNorm 0 -CovTrim 3 -CovTrimLen 55.0\n");
    fflush(stdout);exit(1);
  }

  if(TrimNormFP > COVERAGE_TRIM){
    printf("-TrimNormChim 2nd arg %d cannot exceed -CovTrim %d : TrimNormFP reduced to %d\n",TrimNormFP,COVERAGE_TRIM, COVERAGE_TRIM);
    fflush(stdout);
    TrimNormFP = COVERAGE_TRIM;
  }
  if(FragileQualityFP > COVERAGE_TRIM){
    printf("-FragileQuality 3rd arg %d cannot exceed -CovTrim %d : FragileQualityFP reduced to %d\n",FragileQualityFP,COVERAGE_TRIM, COVERAGE_TRIM);
    fflush(stdout);
    FragileQualityFP = COVERAGE_TRIM;
  }

  if(hashdelta && hashdelta < 2.0 * (DELTA_X + 1)){
    printf("WARNING: -hashdelta %0.1f %0.1f : values should be at least 2 * deltaX + 2 = %d\n",hashdelta,hashdeltaAdjust, 2*DELTA_X + 2);
    fflush(stdout);
  }
  if(hashdelta && hashdelta < 2.0 * (max(DELTA_X,outlierExtendLimit) + 1) && !hashdeltaAdjust){
    printf("WARNING: -hashdelta %0.1f %0.1f : 2nd value should be 10 or greater to support outlierExtend=%d\n",hashdelta,hashdeltaAdjust, outlierExtendLimit);
    fflush(stdout);
  }
											       
  deltaExtXRef = max(DELTA_X,max(deltaExtXRef,deltaXRef));
  deltaExtYRef = max(DELTA_Y,max(deltaExtYRef,deltaYRef));
  if(!outlierLambdaRef)
    outlierLambdaRef = outlierLambda;
  if(!outlierLambdaLabel)
    outlierLambdaLabel = outlierLambdaRef;
  if(!outlierLambdaSize)
    outlierLambdaSize = outlierLambdaRef;

  if(OutlierTypeLabel < 0)
    OutlierTypeLabel = OutlierType;
  if(OutlierTypeSize < 0)
    OutlierTypeSize = OutlierType;

  if(VERB>=2){
    printf("outlierLambda=%0.2e,outlierLambdaRef=%0.2e,outlierLambdaSwitch=%d,outlierLambdaLabel=%0.2e,outlierLambdaSize=%0.2e\n",
	   outlierLambda,outlierLambdaRef,outlierLambdaSwitch,outlierLambdaLabel,outlierLambdaSize);
    fflush(stdout);
  }

  if(0){/* exit with error message */
    fprintf(stdout,"Test : exiting after 60 seconds\n");
    fflush(stdout);
    sleep(60);
    fprintf(stderr,"END of output\n");
    fflush(stdout);exit(0);
  }

  if(parametersfile){
    Cparameters parameters;
    input_err(parametersfile, &parameters);
  }

  if(PoutlierEndCnt > 0){
    int origPoutlierEndCnt = PoutlierEndCnt;
    int origPoutlierCnt = PoutlierCnt;
    if(PoutlierCnt > PoutlierEndCnt){
      int delta = PoutlierCnt - PoutlierEndCnt;
      for(int i = delta; i < PoutlierCnt; i++)
	PoutlierRefine[i-delta] = PoutlierRefine[i];
      PoutlierCnt -= delta;
    } else if(PoutlierCnt < PoutlierEndCnt){
      int delta = PoutlierEndCnt - PoutlierCnt;
      for(int i = PoutlierEndCnt; --i >= delta;)
	PoutlierRefine[i] = PoutlierRefine[i-delta];
      for(int i = 0; i < delta; i++)
	PoutlierRefine[i] = PoutlierEndRefine[i];
    }
    if(VERB>=2){
      printf("PoutlierEndCnt = %d -> %d, PoutlierCnt = %d -> %d\n",origPoutlierEndCnt, PoutlierEndCnt, origPoutlierCnt, PoutlierCnt);
      fflush(stdout);
    }
  }

  if(minSNRestimate){
    // mapSNR = 1; 

    if(usecolor){
      minSNR[usecolor-1] = max(SNRmin, min(SNRmax,minSNR[usecolor-1]));
      if(VERB>=2){
	printf("minSNR[%d] -> %0.3f\n",usecolor-1,minSNR[usecolor-1]);
	fflush(stdout);
      }
    } else
      for(int c = 0; c < colors; c++){
	minSNR[c] = max(SNRmin, min(SNRmax,minSNR[c]));
	if(VERB>=2){
	  printf("minSNR[%d] -> %0.3f\n",c,minSNR[c]);
	  fflush(stdout);
	}
      }
  }

  if(colors >= 2 && usecolor >= 2){/* swap parameters for color "usecolor" with parameters for 1st color */
    if(VERB){
      printf("Swapping colors 1 and %d due to -usecolor %d\n",usecolor,usecolor);
      fflush(stdout);
    }
    swapparam(usecolor,NULL,0);
  }

  if(VERB>=2){
    printf("After parsing command line:num_files=%d,num_reffiles=%d\n",num_files,num_reffiles);
    fflush(stdout);
  }

  if(outlierExtend > 0 && outlierExtendLimit > 0 && outlierExtendLimit <= min(DELTA_X,DELTA_Y)){
    printf("WARNING: -outlierExtend %d %d has no effect with deltaX=%d deltaY=%d\n",outlierExtend,outlierExtendLimit,DELTA_X,DELTA_Y);
    fflush(stdout);
    outlierExtend = 0;
  }
  /*  if(outlierExtend > 0 && spots_filename[0] && !(biasWT <= 0.0)){
    printf("WARNING: -outlierExtend not validated with -ref and -biaswt %0.3f\n",biasWT);
    fflush(stdout);
    }*/
  if(PairMerge > 0.0 && PairMergeOverlapped > PairMerge){
    printf("WARNING: -pairmergeOverlapped %0.3f is inconsistent with -pairmerge %0.3f ... : changing to -pairmergeOverlapped %0.3f\n",
	   PairMergeOverlapped, PairMerge, PairMerge);
    fflush(stdout);

    PairMergeOverlapped = PairMerge;
  }

  if(PairMergeHmap){
    if(DEBUG) assert(TruncateHmap >= 0.0);
    if(DEBUG) assert(TruncateNHmap >= 0.0);
  }

  if(PairsplitMinOutlier > 0 && 
     ((outlierExtend==0) ? (PairsplitMinOutlier > DELTA) :
      outlierExtendLimit ? (PairsplitMinOutlier > max(DELTA,outlierExtendLimit)) : 0)){
    printf("-pairsplitMinOutlier %d cannot be used with -delta %d and -outlierExtend %d %d: pairsplitMinOutlier must be no larger than largest outlier size\n",
	   PairsplitMinOutlier, DELTA, outlierExtend, outlierExtendLimit);
    fflush(stdout);exit(1);
  }

  /*  if(PairMergeRepeat && Cfirst != 0 && CfirstPairs == 0){
    printf("-pairmergeRepeat not supported with -first, unless -firstPairs 1 or -firstPairs 2 is used: Either remove pairmergeRepeat or add -firstPairs\n");
    fflush(stdout);exit(1);
    }*/

  if(Refine && MapRate){
    printf("Cannot use -MapRate %0.3f with -refine %d\n", MapRate, Refine);
    fflush(stdout);exit(1);
  }

  if(HashMultiMatch && !spots_filename[0] && !CmapMerge){
    printf("-hashMultiMatch not supported for pairwise alignment : If aligning large consensus maps, try using -pairsplit without hashtable to find multiple matchgroups\n");
    fflush(stdout);exit(1);
  }

  if(RepeatRec && RepeatMaxShift > 0 && RepeatLogPvalueRatio > 0.0 && !extend){/* make sure extend >= 1 */
    printf("WARNING: -RepeatMask with -RepeatRec requires -extend 1 or -extend 2 (using -extend 1)\n");
    fflush(stdout);

    extend = 1;
  }

  //  if(bgSNRcnt[0] > 0)
  //    MapSNR = 1;/* force SNR values to be present in input files */

  MapSNROrig = MapSNR;
  MapStitchedOrig = MapStitched;
  MapIntensityOrig = MapIntensity;
  MapStitchLocOrig = MapStitchLoc;
  MapPSFWidthOrig = MapPSFWidth;

  if(groupfile && MappedUnsplit >= 0){
    printf("Cannot specify both -grouped and -mapped-unsplit\n");
    fflush(stdout);exit(1);
  }

  NumScaleFactor = 1+2*ScaleDelta;
  int NumScaleFactor2 = 1 + 2*MapScaleDelta;
  ScaleFactor = new double[max(NumScaleFactor,NumScaleFactor2)];
  ScaleFactor[0] = 1.0;
  for(int i = 1; i <= ScaleDelta; i++){
    double scale = 1.0 + i * ScaleDeltaSize;
    ScaleFactor[2*i-1] = scale;
    ScaleFactor[2*i] = 1.0/scale;
  }
  if(VERB>=2){
    for(int i = 0; i <= ScaleDelta*2; i++)
      printf("ScaleFactor[scaleID=%d]= %0.6f\n", i, ScaleFactor[i]);
    fflush(stdout);
  }

  long long MemSize = 0, SwapSize = 0, AvailableMem = 0;
  double GB = 1024.0 * 1024.0 * 1024.0;

#ifndef WIN32
  if(MaxMem > 0.0){

    /* call "free -m" to see how much actual swap space is available : reduce MaxVirtMem to this value + MaxMem */
#if USE_MIC
    getswap(MemSize,SwapSize);/* Total Mem, free Swap (and Available Mem) based on "free -m" */
    AvailableMem = MemSize;
#else
    getswap(MemSize,SwapSize, &AvailableMem);/* Total Mem, free Swap (and Available Mem) based on "free -m" */
#endif
      
    long long VmSize, VmRSS, VmSwap;
    getmem(VmSize,VmRSS,VmSwap);

    uid_t euid = geteuid();
    if(VERB/* HERE HERE >=2 */){
      if(MaxVirtMem != 0.0)
	 printf("euid= %d, TotalThreads=%d,%d total Mem= %0.4f(available= %0.4f), VmRSS= %0.4f, free Swap= %0.4f, VmSwap= %0.4f, MaxMem= %0.4f, MaxVirtMem= %0.4f, Margin= %0.4f GB\n",
		(int)euid, TotalThreads,TotalSlots, MemSize / GB, AvailableMem / GB, VmRSS / GB, SwapSize / GB, VmSwap / GB, MaxMem, MaxVirtMem, MaxMemIncrease);
      else
	 printf("euid= %d, TotalThreads=%d,%d total Mem= %0.4f(available= %0.4f), VmRSS= %0.4f, free Swap= %0.4f, VmSwap= %0.4f, MaxMem= %0.4f, Margin= %0.4f GB\n",
		(int)euid, TotalThreads,TotalSlots, MemSize / GB, AvailableMem / GB, VmRSS / GB, SwapSize / GB, VmSwap / GB, MaxMem, MaxMemIncrease);
      fflush(stdout);
    }

    //    MemSize = AvailableMem + VmRSS; // NOTE : this would exclude memory used by other RefAligner jobs
    SwapSize += VmSwap;

#if USE_MIC
    if(MemSize > 8.0 * GB){
      if(VERB){
	printf("Increasing -maxmem from %0.4f Gb to %0.4f Gb based on estimated system MemSiz= %0.4f Gb\n",
	       MaxMem, MaxMem * MemSize / (8.0 * GB), MemSize / GB);
	fflush(stdout);
      }
      MaxMem *= MemSize / (8.0 * GB);
    }
#else
    if(MemSize < (MaxMem + MaxMemIncrease) * GB){
      if(VERB){
	printf("Changed -maxmem from %0.4f Gb to %0.4f Gb based on estimated system memory= %0.4f, Margin = %0.4f GB\n",
          MaxMem, MemSize / GB - MaxMemIncrease, MemSize / GB, MaxMemIncrease);
	fflush(stdout);
      }
      if(MemSize < (MaxMem + MaxMemIncrease)*GB)
	MaxMem = MemSize / GB - MaxMemIncrease;
    }
    if(MaxMemIncrease > 0.0 && MemSize > (MaxMem + MaxMemIncrease) * GB) {
      if(VERB){
	printf("Changed -maxmem from %0.4f Gb to %0.4f Gb based on estimated system memory= %0.4f, Margin= %0.4f GB\n",
	       MaxMem, MemSize / GB - MaxMemIncrease, MemSize / GB, MaxMemIncrease);
	fflush(stdout);
      }
      MaxMem = MemSize / GB - MaxMemIncrease;
    }
#endif    

    if(MaxVirtMem < 0.0)
      MaxVirtMem = 3.0 * MaxMem;
    if(MaxVirtMem != 0.0){
      double MaxVirt = SwapSize / GB + MaxMem;
      if(VERB && MaxVirt != MaxVirtMem){
	printf("MemSize = %lld, SwapSize = %lld : Setting MaxVirtMem = %0.4f Gb (MaxMem= %0.4f Gb)\n",MemSize, SwapSize, MaxVirt, MaxMem);
	fflush(stdout);
      }
      MaxVirtMem = MaxVirt;
    }
    if(HMaxMem <= 0.0)
      HMaxMem = MaxMem;
    else if(HMaxMem > MaxMem){
      if(VERB){
	printf("Reducing -hashmaxmem from %0.4f to %0.4f Gb\n",HMaxMem, MaxMem);
	fflush(stdout);
      }
      HMaxMem = MaxMem;
    }
  }
#endif // ifndef WIN#2

  if(MaxThreads < 0)/* default initialization */
    MaxThreads = 0;

  int nthreads = 1;

#ifdef _OPENMP
  nthreads = omp_get_max_threads();
  if(MaxThreads <= 0)
    MaxThreads = nthreads;
  if(TotalThreads == 1)
    TotalThreads = nthreads;
  if(TotalSlots > 0) {
    int nTotalThreads = floor((TotalThreads * nthreads + TotalSlots*0.5) / TotalSlots);
    if(VERB){
      printf("Adjusting TotalThreads %d -> %d : omp_get_max_threads()= %d, TotalSlots= %d\n",TotalThreads,nTotalThreads,nthreads,TotalSlots);
      fflush(stdout);
    }
    TotalThreads = nTotalThreads;
  }
#endif

  if(MaxThreads < TotalThreads){
    if(VERB){
      printf("Reducing maxmem from %0.2f to %0.2f due to TotalThreads=%d,MaxThreads=%d\n",MaxMem, (MaxMem * MaxThreads)/TotalThreads, TotalThreads, MaxThreads);
      printf("Reducing maxvirtmem from %0.2f to %0.2f due to TotalThreads=%d,MaxThreads=%d\n",MaxVirtMem, (MaxVirtMem * MaxThreads)/TotalThreads, TotalThreads, MaxThreads);
      printf("Reducing hashmaxmem from %0.2f to %0.2f due to TotalThreads=%d,MaxThreads=%d\n",HMaxMem, (HMaxMem * MaxThreads)/TotalThreads, TotalThreads, MaxThreads);
      fflush(stdout);
      MaxMem = (MaxMem * MaxThreads) / TotalThreads;
      HMaxMem = (HMaxMem * MaxThreads) / TotalThreads;
      MaxVirtMem = (MaxVirtMem * MaxThreads) / TotalThreads;
    }
  }

#ifdef _OPENMP
  if(MaxThreads > nthreads){
    printf("WARNING: system has only %d threads : reducing -maxthreads from %d to %d\n",nthreads,MaxThreads,nthreads);
    fflush(stdout);
    MaxThreads = nthreads;
  }
  if(RefineThreads <= 0)
    RefineThreads = MaxThreads;
  if(BoostThreads){
    int origBoostThreads = BoostThreads;
    if(BoostThreads == 1)
      BoostThreads = nthreads;
    BoostThreads = min(nthreads, 2*RefineThreads);
    if(BoostThreads < min(nthreads, 2*origBoostThreads))
       BoostThreads = min(nthreads, 2*origBoostThreads);
    if(VERB && BoostThreads){
      printf("-boostThreads %d %d : number of threads will be increased to %d after %d seconds\n",
	     BoostThreadsTime,origBoostThreads,BoostThreads,BoostThreadsTime);
      fflush(stdout);
    }
  }
#endif

#if 0 // USE_MIC
  if(MaxMem > AvailableMem / GB){
    if(VERB){
      printf("Reducing MaxMem from %0.2f to %0.2f Gb due to Available Mem\n", MaxMem, AvailableMem / GB);
      fflush(stdout);
    }
    MaxMem = AvailableMem / GB;
  }
#endif

  int numthreads = 1;
  numthreads = MaxThreads;
  numthreads = max(1,min(numthreads, totalmaps/1024));

  int inputThreads = MaxThreads;
  if(Cfirst < 0 || PairSplit)
    inputThreads = 1;
  if(USE_MIC)
    inputThreads = min(16,inputThreads);

  if(colors == 2){
    SF[colors] = MA_SF;
    SD[colors] = MA_SD;
  }

  if(hashsplit && !hash_filename){
    printf("-hashsplit requires -hash filename\n");
    fflush(stdout);exit(1);
  }

  if(hashsplit && numbins <= 1){
    printf("-hashsplit requires -partial with at least 2 bins\n");
    fflush(stdout);exit(1);
  }

  if(MappedPrefix && !spots_filename[0]){
    printf("-mapped cannot be used without -ref (-mapped will be ignored)\n");
    fflush(stdout);
  }

  if(MultiMatches && !CutFlip){/* initialize LogPvThreshold3 etc */
    ScoreThreshold4 = ScoreThreshold3 = ScoreThreshold2;
    LogPvThreshold4 = LogPvThreshold3 = LogPvThreshold2;
  }

  if(XmapFilterThresh1 > 0.0){
    double origLogPvThreshold = LogPvThreshold;
    LogPvThreshold = XmapFilterThresh2 - XmapFilterDelta;
    if(VERB){
      printf("Changing -T confidence from %0.2f to %0.2f due to -xmapFilter (Final thresholds are %0.2f and %0.2f for unambigous and ambigous alignments respectively)\n", 
	     origLogPvThreshold,LogPvThreshold, XmapFilterThresh2, XmapFilterThresh1);
      fflush(stdout);
    }
  }

  if(RefRepeats > 1 && Refine >= 1){
    printf("-refine %d ignored due to -M %d\n",Refine, RefRepeats);
    fflush(stdout);
    Refine = 0;
  }
  if(RefRepeats > 1 && NoStat){
    printf("-nostat ignored due to -M %d\n",RefRepeats);
    fflush(stdout);
    NoStat = 0;
  }
  if(extendonly && !(Refine >= 2 && extend)){
    printf("-extendonly ignored without -refine 2 and -extend\n");
    fflush(stdout);
    extendonly = 0;
  }
  if(spots_filename[0]){
    for(int c = 0; c < colors; c++)
      if(resSD[c] <= 0.0){
	printf("-ref not supported with resSD[%d]= 0.0 : changed to 0.001\n", c+1);
	fflush(stdout);
	resSD[c] = 0.001;
      }
  }
  if(SecondBest && !spots_filename[0]){
    printf("-SecondBest not supported without -ref\n");
    fflush(stdout);exit(1);
  }
  if(MultiMatches){
    if(LogPvThreshold2 < -99.0)
      LogPvThreshold2 = LogPvThreshold;
    if(ScoreThreshold2 < -1.0e29)
      ScoreThreshold2 = ScoreThreshold;

    if(!MULTIMATCHES_FIX){
      if(MultiMatches > AlignedSiteThreshold){
	printf("WARNING: -MultiMatches %d is inconsistent with -A %d : changing to -A %d\n",MultiMatches, AlignedSiteThreshold,MultiMatches);
	fflush(stdout);
	AlignedSiteThreshold = MultiMatches;
      }
      if(MultiMatches > AlignedSiteThreshold2){
	printf("WARNING: -MultiMatches %d is inconsistent with -A2 %d : changing to -A2 %d\n",MultiMatches, AlignedSiteThreshold2,MultiMatches);
	fflush(stdout);
	AlignedSiteThreshold2 = MultiMatches;
      }
    }

    LogPvThreshold2 = min(LogPvThreshold2,LogPvThreshold);
    ScoreThreshold2 = min(ScoreThreshold2,ScoreThreshold);
    AlignedSiteThreshold2 = min(AlignedSiteThreshold2,AlignedSiteThreshold);
    AlignedLengthThreshold2 = min(AlignedLengthThreshold2,AlignedLengthThreshold);
    AlignedEndOutlierThreshold2 = max(AlignedEndOutlierThreshold2,AlignedEndOutlierThreshold);
    AlignedOutlierThreshold2 = max(AlignedOutlierThreshold2,AlignedOutlierThreshold);
    if(VERB>=2){
      printf("MultiMatches=%d,S2=%0.4f,T2=%0.2f,A2=%d,L2=%0.3f,E2=%d,I2=%0.3f\n",MultiMatches,ScoreThreshold2,LogPvThreshold2,AlignedSiteThreshold2,AlignedLengthThreshold2,AlignedEndOutlierThreshold2,AlignedOutlierThreshold2);
      fflush(stdout);
    }
  }
  if(!(MultiMatches || SecondBest)){
    LogPvThreshold2 = LogPvThreshold;
    ScoreThreshold2 = ScoreThreshold;
    AlignedSiteThreshold2 = AlignedSiteThreshold;
    AlignedLengthThreshold2 = AlignedLengthThreshold;
    AlignedEndOutlierThreshold2 = AlignedEndOutlierThreshold;
    AlignedOutlierThreshold2 = AlignedOutlierThreshold;
  }
  if(!TEinit)
    LogPvThresholdTE = LogPvThreshold; // WAS LogPvThreshold2;
  if(!AEinit)
    AlignedSiteThresholdE = AlignedSiteThreshold;
  if(!LEinit)
    AlignedLengthThresholdE = AlignedLengthThreshold;

  if(!TSinit)
    LogPvThresholdTS = LogPvThresholdTE;
  if(!ASinit)
    AlignedSiteThresholdS = AlignedSiteThresholdE;
  if(!LSinit)
    AlignedLengthThresholdS = AlignedLengthThresholdE;
    
  if(extsplitOutlier > 0.0 && extsplitOutlierLS < 0.0)
    extsplitOutlierLS = AlignedLengthThresholdS;
  if(extsplitOutlier > 0.0 && extsplitOutlierAS < 0.0)
    extsplitOutlierAS = AlignedSiteThresholdS;

  if(XmapChim && XmapChimSpacing && !MultiMatches){
    printf("WARNING: Ignoring chimeric alignments with same reference map : requires -MultiMatches\n");
    fflush(stdout);
    XmapChimSpacing = 0.0;
  }
  if(SecondBest && MultiMatches){
    printf("-SecondBest not supported with -MultiMatch\n");
    fflush(stdout);exit(1);
  }

  if(SecondBest){
    for(int c = 0; c < colors; c++)
      if(resSD[c] <= 0.0){
	printf("-SecondBest not supported with resSD[%d]= 0.0 : changed to 0.001\n", c+1);
	fflush(stdout);
	resSD[c] = 0.001;
      }
  }
  if(Refine){
    for(int c = 0; c < colors; c++)
      if(resSD[c] <= 0.0){
	printf("-refine not supported with resSD[%d]= 0.0 : changed to 0.001\n", c+1);
	fflush(stdout);
	resSD[c] = 0.001;
      }
  }
  if(extendonly && ContigSplitRatio > 0.0){
    printf("-contigsplit ignored due to -extonly\n");
    fflush(stdout);
    ContigSplitRatio = 0.0;
  }
  if(TBinit && !(Refine >= 2 && extend)){
    printf("-TB ignored without -refine 2 and -extend\n");
    fflush(stdout);
    TBinit = 0;
    TBcnt = 0;
  }
  if(TEinit && LogPvThresholdTE > LogPvThreshold + 0.001){
    printf("WARNING:-TE threshold is more stringent than -T threshold\n");
    fflush(stdout);exit(1);
  }
  if(TBcnt > 0){/* check that all -TB values are in order of increasing stringency */
    for(int i = 1; i < TBcnt; i++){
      if(LogPvThresholdTB[i] < LogPvThresholdTB[i-1]){
	printf("%d'th -TB value is less strigent than the previous value\n",i+1);
	fflush(stdout);exit(1);
      }
    }
  }

  if(maxresbias > 0.0 && !(maxresbias > mres * 0.5)){
    printf("-resbias cutoff must exceed -mres * 0.500 = %0.3f kb: ignoring -resbias\n", mres*0.5);
    fflush(stdout);
    maxresbias = 0.0;
  }

#if 0
  if(!parametersfile && !(maxresbias > mres * 0.5) && Refine){ /* Make sure mres >= rres */
    if(mres < rres){
      if(VERB){
	printf("Changing -mres %0.3f to -mres %0.3f\n",mres, rres);
	fflush(stdout);
      }
      mres = rres;
    }
  }
#endif

  double useval = 0; /* for scaffold and sv -- sv will be prioritized over scaffold for this if both */
  bool printext = false, printbw = false, printa = false, printt = false; /* report to user that extend, biaswt, A was set automatically */
  if(doscaffold) {
    useval = 1e-4; /* current recommendation is this for endoutlier and outlier--pairsplit not supported */
    if(!spots_filename[0]) {
      printf("ERROR: -scaffold is only supported for -ref aligning\n");
      fflush(stdout);exit(1);
    }
    if(dosvdetect || PairSplit) { /* currently, this isn't supported because sv is only pairwise, and scaffold is only ref */
      printf("ERROR: cannot do sv or pairsplit and scaffold simultaneously\n");
      fflush(stdout);exit(1);
    }
    if( !setextend ) {
      extend = 1;
      printext = true;
    }
    if( !setbiaswt ) {
      biasWT = 0;
      printbw = true;
    }
    if( !seta ) {
      AlignedSiteThreshold = 6; /* 6 sites */
      printa = true;
    }
    if( !sett ) {
      LogPvThreshold = -log(1e-6)/log(10.0); /* 10^-6 for 6 sites */
      printt = true;
    }
    input_scaffold(scaf_filename);
  }

  if( svlvl == 1 ) /* use svlvl to set pairsplit and endoutlier */
    useval = 1e-6;
  else if( svlvl == 2 )
    useval = 1e-5;
  else if( svlvl == 3 )
    useval = 1e-4;
  else if( svlvl != 0 ) { /* svlvl defaults to 0 now, anything else is an error */
    printf("-sv arg must be 1, 2, or 3--command line parsing error\n");
    fflush(stdout);exit(1);
  }

  if(LogPvThreshold2 < 0.0)
    LogPvThreshold2 = LogPvThreshold;
  if(ScoreThreshold2 < 0.0)
    ScoreThreshold2 = ScoreThreshold;
  if(AlignedSiteThreshold2 < 0)
    AlignedSiteThreshold2 = AlignedSiteThreshold;
  if(AlignedLengthThreshold2 < 0.0)
    AlignedLengthThreshold2 = AlignedLengthThreshold;
  if(AlignedEndOutlierThreshold2 < 0.0)
    AlignedEndOutlierThreshold2 = AlignedEndOutlierThreshold;
  if(AlignedOutlierThreshold2 < 0.0)
    AlignedOutlierThreshold2 = AlignedOutlierThreshold;
  if(AlignedOutlierLabels2 < 0)
    AlignedOutlierLabels2 = AlignedOutlierLabels;

  if(ScoreThreshold3 < 0.0)
    ScoreThreshold3 = ScoreThreshold2;
  if(LogPvThreshold3 < 0.0)
    LogPvThreshold3 = LogPvThreshold2;
  if(ScoreThreshold4 < 0.0)
    ScoreThreshold4 = ScoreThreshold3;
  if(LogPvThreshold4 < 0.0)
    LogPvThreshold4 = LogPvThreshold3;

  if(PairSplit && MultiMatches && RefSplit){
    printf("-RefSplit cannot be used with -PairSplit\n");
    fflush(stdout);exit(1);
  }
  if (MultiMatches) {
    if(BestRef){
      printf("-MultiMatches has no effect with -BestRef: Disabling MultiMatches\n");
      fflush(stdout);
      MultiMatches = 0;
    }
    if(NoSplit < 2){
      printf("-MultiMatches requires -nosplit 2\n");
      fflush(stdout);exit(1);
    }
    if(NumScaleFactor > 1 && hashScaleDelta < 2){
      printf("WARNING: -MultiMatches with -ScaleDelta not validated, unless -hashScaleDelta 2 is specified\n");
      fflush(stdout);// exit(1);
    }
    if(BestRefOutlier <= 999.0){
      printf("-MultiMatches with -BestRefOutlier not supported\n");
      fflush(stdout);exit(1);
    }
    if(PoutlierEnd <= 0.0){
      printf("-MultiMatches requires non-zero -endoutlier value\n");
      fflush(stdout);exit(1);
    }
  }
  if(RefSplit && !MultiMatches){
    printf("WARNING: -RefSplit not supported without -MultiMatches: Disabling RefSplit\n");
    fflush(stdout);
    RefSplit = 0;
  }

  bool printeo = false, printps = false, printo = false; /* these are just to shorten the stdout */
  if(dosvdetect && !(MultiMatches && RefSplit)) { /* need pairsplit for sv detect--turn it on with default 1e-6--print this */
    printf("-sv requires -RefSplit and -MultiMatches\n");
    fflush(stdout);exit(1);
  }

  /*  if(MultiMatches && RefSplit && Psplit < Poutlier){
    printf("WARNING: -RefSplit pvalue cannot be more stringent than -outlier pvalue: changing to -RefSplit %0.2e %d\n",Poutlier,RefSplitMinLabels);
    fflush(stdout);
    Psplit = Poutlier;
    }*/

  /* need -endoutlier for sv detect--turn it on with default 1e-6--print this */
  /* for scaffolding, it is recommended, but at 1e-4 */
  if((dosvdetect || doscaffold) && !PoutlierEnd_init) { 
    printeo = true;
    PoutlierEnd_init = 1;
    PoutlierEnd = useval;
  }

  if(dosvdetect && !setoutlier) { /* always set -outlier same as -pairsplit for -sv (unless you set it yourself) */
    printo = true;
    Poutlier = useval;
  } else if(dosvdetect && Poutlier < Psplit) { /* -outlier value _must_ be >= to -pairsplit value--set it equal if it's < */
    if(RefSplit){
      printf("Changing RefSplit value from %f to equal outlier Pvalue value (%f) because -sv was used\n", Psplit,Poutlier);
      fflush(stdout);
      Psplit = Poutlier;
    } else {
      printf("Changing outlier Pvalue from %f to equal pairsplit value (%f) because -sv was used\n", Poutlier,Psplit);
      fflush(stdout);
      Poutlier = Psplit;
    }
  } else if(doscaffold && !setoutlier) { /* for scaffold, use 1e-4 also as a baseline if not set */
    printo = true;
    Poutlier = useval;
  }

  /* print warning to user that options are used implicitly */
  if(printps || printeo || printo || printext || printbw || printa || printt) { 
    printf("Using: ");
    if(printps) /* pairsplit */
      printf("-pairsplit %f ", Psplit);
    if(printeo) /* endoutlier */
      printf("-endoutlier %f ", PoutlierEnd);
    if(printo) /* outlier */
      printf("-outlier %f ", Poutlier);
    if(printext) /* extend */
      printf("-extend %i ", extend);
    if(printbw) /* biaswt */
      printf("-biaswt %.1f ", biasWT);
    if(printa) /* A */
      printf("-A %i ", AlignedSiteThreshold);
    if(printt) /* T */
      printf("-T %f ", LogPvThreshold ); /* may as well just print the log, though I guess it's slightly misleading */
    if(dosvdetect)
      printf("because -sv used\n");
    else if(doscaffold)
      printf("because -scaffold used\n");
  }

  if(dosvdetect && mask_filename) //this is ignored if not -sv
    input_mask(mask_filename);
  else if(mask_filename) //warn user that this argument is ignored in this case
    printf("Warning: mask file (%s) is ignored because -sv is not enabled\n", mask_filename);

  if(conf_filename)
    input_confidence(conf_filename);

  if(outlierExtend > 0 && PoutlierEnd <= 0.0){
    printf("-outlierExtend %d ignored due to -outlier 0\n",outlierExtend);
    fflush(stdout);
    outlierExtend = 0;
  }
  if(Ch > 0.0 && ChFdr > 0.0 && PoutlierEnd > 0.0){
    printf("Using -endoutlier AND -local : RefAligner will use whichever is less stringent (refinement only uses -endoutlier)\n");
    fflush(stdout);
  }
  if(PairSplit && spots_filename[0]){
    printf("-pairsplit or -querysplit can only be used with pairwise alignment (no -ref option allowed)\n");
    fflush(stdout);exit(1);
  }
  if(PairSplit && !(PoutlierEnd > 0.0)){
    printf("Cannot perform pairwise multiple local alignments (-pairsplit) without -endoutlier\n");
    fflush(stdout);exit(1);
  }
  if(PairSplit && !( Cfirst == 0 && !spots_filename[0])){
    if(num_files < 2){
      printf("-pairsplit must use at least 2 input (-i) files\n");
      fflush(stdout);exit(1);
    } else if(num_files > 2){
      printf("WARNING:-pairsplit assumes first -i input file contains the reference maps (use -first to modify)\n");
      fflush(stdout);
    }
  }
  if(PairSplit && !Cfirst)
    Cfirst = -1;
  if(Cfirst && skipcnt > 0){
    printf("Cannot combine -first with -partial with skipcnt=%d\n",skipcnt);
    fflush(stdout);exit(1);
  }
  if(PairSplit && numbins > 1){
    printf("Cannot perform pairwise multiple location alignments with -partial\n");
    fflush(stdout);exit(1);
  }
  /*  if(num_files > 1 && spots_filename[0]){
    printf("Cannot use %d -i input files with -ref : use -merge to combine the -i inputs into a single file first\n",num_files);
    fflush(stdout);exit(1);
    }*/
  if(!num_files && !CmapMerge){
    printf("No map input file specified: only -merge supported with null input map set\n");
    fflush(stdout);exit(1);
  }
  if(!output_prefix){
    printf("-o option must be specified\n");
    fflush(stdout);exit(1);
  }

  if(HashGen){/* generate hash_prefix */
    char buf[PATH_MAX];
    sprintf(buf,"%s_p%lld_id%lld", output_prefix, (long long)pid, CMapID); 
    hash_prefix = strdup(buf);
    if(VERB){
      printf("hash_prefix = %s\n",hash_prefix);
      fflush(stdout);
    }
    if(hash_filename && hash_filename != (char *)(-1)){
      printf("WARNING: cannot use -hash to specify a hash_filename when using -hashgen : replacing %s with %s.hashbin\n",hash_filename,hash_prefix);
      fflush(stdout);
      hash_filename = (char *)(-1);
    }
  }

  if(hash_filename == (char *)(-1)){
    char filename[PATH_MAX];
    if(hash_prefix){ /* if -hashgen already set the hashtable file prefix, use it */;
      if(DEBUG) assert(HashGen);
      sprintf(filename,"%s.hashbin", hash_prefix);
    } else /* make up a default value */
      sprintf(filename,"%s_id%lld.hashbin", output_prefix, CMapID); 
    
    hash_filename = strdup(filename);
    hash_threshold = HashScore;

    if(VERB){
      printf("hash_filename = %s\n",hash_filename);
      fflush(stdout);
    }
  } 

  draft_prefix = output_prefix;
  if((CMapID > 0 && numrefmaps==1) && !CmapMerge){
    char buf[PATH_MAX];
    sprintf(buf,"%s_contig%lld",output_prefix,CMapID);
    output_prefix = strdup(buf);
  }

  if(extendSplit && MappedPrefix && extendLen > 0.0 && extendWT > 0.0){
    printf("WARNING: Cannot use -extsplit with BOTH -mapped AND -extendWT : ignoring -extsplit\n");
    fflush(stdout);

    extendSplit = 0;
  }
  if(MultiMatches && AddFragsWT && MappedPrefix){
    printf("WARNING: Cannot use -MultiMatches -BestRefWT AND -AddFragsWT with -mapped : ignoring -AddFragsWT\n");
    fflush(stdout);
    
    AddFragsWT = 0;
  }

  if(SplitCnt < 0 && ContigCntFile && (Refine || SplitHmap) && (extendSplit || ContigSplitRatio > 0.0 || SplitHmap)){
    SplitCnt = nextContigId(ContigCntFile);
    if(DEBUG) assert(SplitCnt > 0);
    if(VERB/* HERE >=2 */){
      printf("nextContigId(\"%s\") returned SplitCnt=%lld\n", ContigCntFile, SplitCnt);
      fflush(stdout);
    }
  }

  if(!colors){
    printf("-colors parameter missing\n");
    fflush(stdout);exit(1);
  }
  for(register int i=0;i<colors;i++){
    if(FP[i] < 0.0){
      printf("False Postive density FP[%d] = %0.6f is too small (must be > 0)\n", i, FP[i]);
      fflush(stdout);exit(1);
    }
    if(FN[i] < 0.0){
      printf("False Negative rate FN[%d] = %0.6f is too small (must be > 0)\n", i, FN[i]);
      fflush(stdout);exit(1);
    }
    if(FN[i] > 1.0){
      printf("False Negative rate FN[%d] = %0.6f is too large (must be < 1)\n", i, FN[i]);
      fflush(stdout);exit(1);
    }
    if(SF[i] < 0.0){
      printf("Sizing error SF[%d] = %0.6f is too small (must be > 0)\n", i, SF[i]);
      fflush(stdout);exit(1);
    }
    if(SR[i] < 0.0){
      printf("Sizing error SR[%d] = %0.6f is too small (must be > 0)\n", i, SR[i]);
      fflush(stdout);exit(1);
    }
  }
  if(Cfirst < 0 && -Cfirst > num_files){
    printf("-first %d not compatible with %d input files\n",Cfirst,num_files);
    fflush(stdout);exit(1);
  }

  if(PairSplit){
    if(Cfirst < 0){
      RefIndex = -Cfirst;
      QueryIndex = -Cfirst;
    } else {
      if(DEBUG)assert(Cfirst > 0);
      RefIndex = num_files;
      QueryIndex = 0;
    }
  }

// -bpp adjustment has been moved to mapcorrect() in input_vfixx.cpp
  
  if(svcheck && InDel){
    printf("WARNING: Cannot use -indel with -svcheck (since both write to the same .indel filename) : ignoring -indel\n");
    fflush(stdout);
    InDel = 0;
  }

  //  printf("Before reading input maps:Sleeping 60 seconds\n");
  //  fflush(stdout);
  //  sleep(60);

  if(NoBpp){
    printf("Changing bpp= %0.3f -> %0.3f due to -NoBpp\n",bpp * 1000.0, 500.0);
    fflush(stdout);
    bpp = 0.500;
  }
  PixelLen = bpp;

  if(ScanCorrection && (Refine || NoStat)){
    printf("WARNING:Cannot use -ScanScaling %d with -refine %d or -nostat : ignoring -ScanScaling\n", ScanCorrection, Refine);
    fflush(stdout);

    ScanCorrection = 0;
  }

  std::set<long long> ids;
  if(SelectIDcnt > 0){
    /* eliminate duplicates to avoid strange warnings later */
    qsort(SelectID,SelectIDcnt,sizeof(long long),(intcmp*)LLInc);
    
    int j = 1;
    for(int i = 1; i < SelectIDcnt; i++){
      if(SelectID[i] == SelectID[j-1]){
	if(VERB>=2){
	  printf("WARNING: duplicate ID= %lld in -selectid inputs\n",SelectID[i]);
	  fflush(stdout);
	}

	continue;
      }
      if(j < i)
	SelectID[j] = SelectID[i];
      j++;
    }
    if(VERB && j < SelectIDcnt){
      printf("WARNING:Removed %d duplicate IDs from -selectid(f), %d IDs remaining\n",SelectIDcnt-j, j);
      fflush(stdout);
    }
    SelectIDcnt = j;

    for(int i = 0; i < SelectIDcnt; i++)
      ids.insert(SelectID[i]);

    if(VERB/* HERE >=2 */){
      printf("Set of %d ids generated from -selectidf or -selectid\n",SelectIDcnt);
      fflush(stdout);
    }
  }

  int UniqueScanIdMax = 0;
  int UniqueScanIdMin = MASK(31);

  double startSNRmin = SNRmin;
  double startSNRmax = SNRmax;
  double startSNR = minSNR[0];

  int Hapnummaps = nummaps;
  int cmap_files = 0;/* number of -i input files that are .cmap or .hmap files */

  if(!spots_filename[0] && !PairSplit && !PairMerge && !RefineOverwrite){/* check if .align file already exists */
    char filename[PATH_MAX];
    FILE *fp;
    sprintf(filename,"%s.align",output_prefix);
    if(VERB>=2){
      printf("Checking existence of %s\n",filename);
      fflush(stdout);
    }

    if((fp = fopen(filename,"r")) != NULL){
      fclose(fp);
      if(VERB){
	printf("WARNING: Output file %s already exists : skipping pairwise alignment\n",filename);
	fflush(stdout);
      }
      goto Lfree;
    }

    if(VERB/* HERE >=2 */){
      printf("Could not find %s : proceeding with pairwise alignment\n",filename);
      fflush(stdout);
    }
  }

  /* read input maps (.bnx or .cmaps format, also supports .vfixx files) */
  if(VERB){
    if(num_files > 1)
      printf("Reading input maps from %s ... %s (%d files)\n",vfixx_filename[0],vfixx_filename[num_files-1],num_files);
    else if(num_files==1)
      printf("Reading input maps from %s\n",vfixx_filename[0]);
    fflush(stdout);
  }

  extern int merrorexit;
  merrorexit = 0;

  if(1){/* adjust number of threads to use for multithreading file input : only implemented for .cmap and .hmap files */
    for(int i = 0; i < num_files; i++){
      char *filename = vfixx_filename[i];
      size_t len = strlen(filename);
      
      if(len >= strlen(".cmap") && (!strcmp(&filename[len-strlen(".cmap")],".cmap") || !strcmp(&filename[len-strlen(".cmap")],".CMAP")))
	cmap_files++;
      
      if(len >= strlen(".hmap") && !strcmp(&filename[len-strlen(".hmap")],".hmap"))
	cmap_files++;
    }  
    if(cmap_files == num_files){// NEW305
      inputThreads = max(1,min(cmap_files,inputThreads));
      inputThreads = min(MaxThreads/2, inputThreads);
    } else
      inputThreads = 1;

    if(VERB/* HERE >=2 */ && inputThreads > 1){
      printf("Using %d threads for input of %d files (including %d CMAP files)\n",inputThreads, num_files,cmap_files);
      fflush(stdout);
    }
  }

  if(!(cmap_files == num_files || cmap_files == 0)){
    printf("ERROR: -i input files must be all be either .bnx or .cmap files: num_files=%d, cmap_files=%d\n",num_files,cmap_files);
    fflush(stdout);exit(1);
  }

  startmaps = 0;

  #pragma omp parallel num_threads(inputThreads) 
  {
    char *mbuf = new char[LINESIZ];

    #pragma omp for schedule(dynamic,1)
    for(int i = 0; i < num_files; i++){
      if(VERB>=2 && num_files > 1){
        printf("Reading input file %s (%d of %d):\n", vfixx_filename[i],i+1,num_files);
	fflush(stdout);
      }
      int origstartmaps = startmaps;
      int inputmaps = input_vfixx(i,nummaps,maxmaps,Gmap, (SelectIDcnt > 0) ? &ids : NULL, (inputThreads > 1),mbuf);

      #pragma omp critical
      {
        startmaps += inputmaps;

	if(VERB>=2){
          int tid = 0;
#ifdef _OPENMP
	  tid = omp_get_thread_num ();
#endif

	  printf("After reading input file %d (%s): %d maps read:(Cfirst=%d),MapSNR=%d,PixelLen=%0.8f,origPixelLen=%0.8f (tid=%d)\n", 
            i, vfixx_filename[i], startmaps-origstartmaps,Cfirst,MapSNR,PixelLen,origPixelLen,tid);
	  if(VERB>=3)
  	    for(int i = origstartmaps; i < startmaps; i++)
	      printf("  Gmap[%d]->id = %lld\n", i, Gmap[i]->id);
          fflush(stdout);
        }
      }

      if(Cfirst < 0 && i+1 == -Cfirst)
        Cfirst = startmaps;
      if(PairSplit && QueryIndex == 0 && RefIndex == num_files && origstartmaps < Cfirst && startmaps >= Cfirst){
        if(startmaps == Cfirst)
          RefIndex = QueryIndex = i+1;
	else {
          RefIndex = i+1;
	  QueryIndex = i;
        }
      }
    } // pragma omp for schedule(dynamic,1)

    delete [] mbuf;
  }// pragam omp parallel num_threads(inputThreads)

  if(DEBUG) assert(startmaps == nummaps);
  if(merrorexit){
    printf("ERROR:One or more input files %s ... had errors\n",vfixx_filename[0]);
    fflush(stdout);exit(1);
  }

  if(VERB>=2 && startmaps > 0){
    Cmap **maplist = new Cmap*[startmaps];

    double totlen = 0.0, totlen2 = 0.0;
    size_t nsites = 0, nsites2 = 0;
    for(int i= 0; i < startmaps; i++){
      Cmap *pmap = Gmap[i];
      maplist[i] = pmap;
      int M = pmap->numsite[0];
      totlen += pmap->site[0][M+1];
      nsites += M;
      if(colors >= 2){
	int M2 = pmap->numsite[1];
	totlen2 += pmap->site[1][M2 + 1];
	nsites2 += M2;
      }
    }
    qsort(maplist, startmaps, sizeof(Cmap*), (intcmp*)PmapLenDec);/* sort map pointers in descending order of length */
    double N50Sum= 0.0, N50Len = 0.0;
    for(int i = 0; i < startmaps; i++){
      Cmap *pmap = maplist[i];
      int M = pmap->numsite[0];
      if((N50Sum += pmap->site[0][M+1]) >= totlen * 0.5){
	N50Len = pmap->site[0][M+1];
	break;
      }
    }
    if(DEBUG && totlen > 0.0) assert(N50Len > 0.0);
	
    if(colors >= 2)
      printf("Maps=%d, sites=%lu,%lu, length= %0.3f kb (avg= %0.3f kb, label density1= %0.3f /100kb, label density2= %0.3f /100kb, N50= %0.4f Mb):wtime=%0.6f\n",
	     startmaps,nsites,nsites2,totlen, 0.5 * (totlen+totlen2)/startmaps, nsites*100.0/max(0.001,totlen), nsites2*100.0/max(0.001,totlen2), N50Len*0.001,wtime());
    else
      printf("Maps=%d, sites=%lu, length= %0.3f kb (avg= %0.3f kb, label density= %0.3f /100kb, N50= %0.4f Mb):wtime=%0.6f\n",
	     startmaps,nsites,totlen,totlen/startmaps, nsites*100.0/max(0.001,totlen), N50Len * 0.001,wtime());

    delete [] maplist;
  }

  if(inputThreads > 1){
    if(VERB/* HERE HERE >=2 */){
      printf("Completed input of %d maps from %d files (including %d CMAP files) using %d threads: calling check_cmap(): time= %0.6f(wall= %0.6f)\n",
	     nummaps,num_files,cmap_files,inputThreads,mtime(),wtime());
      fflush(stdout);
    }
    if(DEBUG) assert(cmap_files == num_files);

    check_cmap(0, NULL, 0, nummaps, maxmaps, Gmap,0);

    startmaps = nummaps;
  }
  if(VERB/* HERE HERE >=2 */){
    printf("Completed input of %d maps: time= %0.6f(wall= %0.6f)\n",nummaps,mtime(),wtime());
    fflush(stdout);
  }

  if(VERB>=2 && startmaps > 0){
    Cmap **maplist = new Cmap*[startmaps];

    double totlen = 0.0, totlen2 = 0.0;
    size_t nsites = 0, nsites2 = 0;
    for(int i= 0; i < startmaps; i++){
      Cmap *pmap = Gmap[i];
      maplist[i] = pmap;
      int M = pmap->numsite[0];
      totlen += pmap->site[0][M+1];
      nsites += M;
      if(colors >= 2){
	int M2 = pmap->numsite[1];
	totlen2 += pmap->site[1][M2 + 1];
	nsites2 += M2;
      }
    }
    qsort(maplist, startmaps, sizeof(Cmap*), (intcmp*)PmapLenDec);/* sort map pointers in descending order of length */
    double N50Sum= 0.0, N50Len = 0.0;
    for(int i = 0; i < startmaps; i++){
      Cmap *pmap = maplist[i];
      int M = pmap->numsite[0];
      if((N50Sum += pmap->site[0][M+1]) >= totlen * 0.5){
	N50Len = pmap->site[0][M+1];
	break;
      }
    }
    if(DEBUG && totlen > 0.0) assert(N50Len > 0.0);
	
    if(colors >= 2)
      printf("Maps=%d, sites=%lu,%lu, length= %0.3f kb (avg= %0.3f kb, label density1= %0.3f /100kb, label density2= %0.3f /100kb, N50= %0.4f Mb):wtime=%0.6f\n",
	     startmaps,nsites,nsites2,totlen, 0.5 * (totlen+totlen2)/startmaps, nsites*100.0/max(0.001,totlen), nsites2*100.0/max(0.001,totlen2), N50Len*0.001,wtime());
    else
      printf("Maps=%d, sites=%lu, length= %0.3f kb (avg= %0.3f kb, label density= %0.3f /100kb, N50= %0.4f Mb):wtime=%0.6f\n",
	     startmaps,nsites,totlen,totlen/startmaps, nsites*100.0/max(0.001,totlen), N50Len * 0.001,wtime());

    delete [] maplist;
  }

  if(CmapChimQuality >= 3 && BreakChiSq > 0.0){/* break input CMAPs based on -BreakChiSq */
    int orignummaps = nummaps;

    extern void mapbreakApply(Cmap **&Map, int &nummaps, int &maxmaps);
    mapbreakApply(Gmap, nummaps, maxmaps);
    
    if(nummaps > orignummaps){
      printf("Applied -breakChiSq options to break select input maps increasing number of maps from %d to %d\n", orignummaps,nummaps);
      fflush(stdout);
    }

    startmaps = nummaps;
  }


#if FORK_DEBUG // DEBUG code to test fork() & waitpid() : currently fopen of new files on parent fails on MIC after fork (errno = 2 or 90, as reported by scif_daemon.c), but works on child!
  fflush(stdout);
  fflush(stderr);

  int F = 1;
  char filename[PATH_MAX];

  sprintf(filename,"%s_test.cmap",output_prefix);
  printf("Before Fork: writing %s\n",filename);
  fflush(stdout);

  errno = 0;
  FILE *fp = fopen(filename,"w");
  if(fp == NULL){
    int eno = errno;
    printf("Unable to fopen %s for write:errno=%d:%s\n",filename,eno,strerror(eno));
    fflush(stdout);
  } else {
    fprintf(fp,"# Header of %s\n",filename);
    FFlush(fp,0);
    FILEclose(fp,0);

    printf("Finished writing %s\n",filename);
    fflush(stdout);
  }

  if(!CMapID)
    sprintf(filename,"%s_C%d.stdout",output_prefix,F);
  else
    sprintf(filename,"%s_id%lld_C%d.stdout",output_prefix,CMapID,F);
  if(VERB/* HERE HERE >=2 */){
    printf("Forking child : stdout will be in %s\n",filename);
    fflush(stdout);
  }

  pid_t pid = fork();
  if(pid == -1){/* fork failed */
    int myerrno = errno;
    fprintf(stderr,"Unable to fork %d'th child: %s\n",F,strerror(myerrno));
    fflush(stdout);exit(1);
  }
  
  if(!pid){/* child */
#if USE_MIC
    extern SCIF_IO_CONTEXT *main_context;
    main_context = NULL;// need to force reset of SCIF interface in child
#endif
    if(freopen(filename,"a",stdout)==NULL){
      int myerrno = errno;
      fprintf(stderr,"Unable to redirect stdout of %d'th child to %s: %s\n",F, filename, strerror(myerrno));
      fflush(stderr);
    }
    printversion(stdout);

    sprintf(filename,"%s_ctest.cmap",output_prefix);

    printf("Child job: writing %s\n",filename);
    fflush(stdout);
    
    errno = 0;
    FILE *fp = fopen(filename,"w");
    if(fp == NULL){
      int eno = errno;
      printf("Unable to fopen %s for write in child:errno=%d:%s\n",filename,eno,strerror(eno));
      fflush(stdout);
    } else {
      fprintf(fp,"# Header of %s\n",filename);
      FFlush(fp,0);
      FILEclose(fp,0);

      printf("Finished writing %s in child\n",filename);
      fflush(stdout);
    }
    printf("Child Debug Exit\n");
    fflush(stdout);
    exit(0);
  } else {// parent
    sprintf(filename,"%s_ptest.cmap",output_prefix);

    printf("Parent job: writing %s\n",filename);
    fflush(stdout);

    errno = 0;
    FILE *fp = fopen(filename,"w");
    if(fp == NULL){
      int eno = errno;
      printf("Unable to fopen %s for write in parent:errno=%d:%s\n",filename,eno,strerror(eno));
      fflush(stdout);
    } else {
      fprintf(fp,"# Header of %s\n",filename);
      FFlush(fp,0);
      FILEclose(fp,0);

      printf("Finished writing %s in parent\n",filename);
      fflush(stdout);
    }
    int status;
    pid_t pid2;

    printf("Waiting for child (pid= %d)\n",pid);
    fflush(stdout);

    while((pid2 = waitpid(pid, &status, 0)) == 0){
      printf("Waiting for child (pid= %d): sleeping 1 second after interrupt\n",pid);
      fflush(stdout);
      sleep(1);
    }
    if(pid2 < 0){
      int myerrno = errno;
      fprintf(stderr,"Error waiting for child (pid= %d): %s\n",pid,strerror(myerrno));
      fflush(stderr);exit(1);
    }
    if(WIFEXITED(status) && WEXITSTATUS(status) != 0){
      printf("Non-zero exit status %d from child (pid= %d)\n",WEXITSTATUS(status),pid);
      fflush(stdout);exit(1);
    }
    if(WIFSIGNALED(status)){
      printf("Child (pid= %d) killed with signal %d\n",pid, WTERMSIG(status));
      fflush(stdout);exit(1);
    }
    printf("Parent Debug Exit\n");
    fflush(stdout);
    exit(0);
  }

  exit(1);
#endif // DEBUG

  extern int RunIndexWarning, RunDataWarning;/* see input_vfixx.cpp */

  if(MappedPrefix && BNXVersion >= 1 && (RunIndexWarning || RunDataWarning || RunDataListLen <= 0)){
    printf("WARNING: Option -BestRefWT with -mapped %s requires BNX 1.2 input to output valid map weights (No map weights will be output)\n",  MappedPrefix);
    fflush(stdout);
  }

  //  printf("After reading input maps:Sleeping 60 seconds\n");
  //  fflush(stdout);
  //  sleep(60);

  if((!QXerror && (DEBUG>=2 || QXmismatch)) || QXfilter){
    if(MapSNR != MapSNROrig){/* check to make sure all maps have SNR information : if not reset MapSNR = 0 */
      if(DEBUG) assert(MapSNROrig==0);
      if(QXfilter)
	MapSNR = 0;
      else {
	for(int i = 0; i < startmaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(!pmap->SNR[0]){
	    if(VERB){
	      printf("WARNING: only a subset of input maps have SNR information : ignoring SNR\n");
	      fflush(stdout);
	    }
	    MapSNR = 0;
	    break;
	  }
	}
      }
    }
    if(MapStitched != MapStitchedOrig){
      if(DEBUG) assert(MapStitchedOrig==0);
      if(QXfilter)
	MapStitched = 0;
      else {
	for(int i = 0; i < startmaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(!pmap->stitch[0]){
	    if(VERB){
	      printf("WARNING: only a subset of input maps have stitch information : ignoring stitch\n");
	      fflush(stdout);
	    }
	    MapStitched = 0;
	    break;
	  }
	}
      }
    }
    if(MapIntensity != MapIntensityOrig){
      if(DEBUG) assert(MapIntensityOrig==0);
      if(QXfilter)
	MapIntensity = 0;
      else {
	for(int i = 0; i < startmaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(!pmap->Intensity[0]){
	    if(VERB){
	      printf("WARNING: only a subset of input maps have Intensity information : ignoring Intensity\n");
	      fflush(stdout);
	    }
	    MapIntensity = 0;
	    break;
	  }
	}
      }
    }
    if(VERB>=2){
      printf("MapStitchLoc=%d, MapStitchLocOrig=%d\n",MapStitchLoc,MapStitchLocOrig);
      fflush(stdout);
    }
    if(MapStitchLoc != MapStitchLocOrig){
      if(DEBUG) assert(MapStitchLocOrig==0);
      if(QXfilter)
	MapStitchLoc = 0;
      else {
	for(int i = 0; i < startmaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(!pmap->stitchLocation[0]){
	    if(VERB){
	      printf("WARNING: only a subset of input maps have stitch Location information : ignoring stitch Location\n");
	      fflush(stdout);
	    }
	    MapStitchLoc = 0;
	    break;
	  }
	}
      }
    }
    if(MapPSFWidth != MapPSFWidthOrig){/* check to make sure all maps have PSFWidth information : if not reset MapPSFWidth = 0 */
      if(DEBUG) assert(MapPSFWidthOrig==0);
      if(QXfilter)
	MapPSFWidth = 0;
      else {
	for(int i = 0; i < startmaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(!pmap->PSFWidth[0]){
	    if(VERB){
	      printf("WARNING: only a subset of input maps have PSFWidth information : ignoring PSFWidth\n");
	      fflush(stdout);
	    }
	    MapPSFWidth = 0;
	    break;
	  }
	}
      }
    }
    if(MapImageLoc != MapImageLocOrig){/* check to make sure all maps have PSFWidth information : if not reset MapImageLoc = 0 */
      if(DEBUG) assert(MapImageLocOrig==0);
      if(QXfilter)
	MapImageLoc = 0;
      else {
	for(int i = 0; i < startmaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(!pmap->ImageFOV[0]){
	    if(VERB){
	      printf("WARNING: only a subset of input maps have Image Location (FOV,X,Y) information : ignoring Image Location\n");
	      fflush(stdout);
	    }
	    MapImageLoc = 0;
	    break;
	  }
	}
      }
    }
  }

  if(MapType >= 0)
    maptype = MapType;

  if(minSNRestimate && !MapSNR){
    printf("WARNING: No SNR values in BNX input : ignoring -minSNRestimate\n");
    fflush(stdout);
    minSNRestimate = 0;
  }
  if(minSNRestimate && MaxSites && !ScanCorrection && !MaxSites2){
    printf("WARNING: -maxsites %d with -minSNRestimate is only applied AFTER estimating minSNR and has no effect without -ScanScaling.\n", MaxSites);
    printf(" Please use a 2nd arg for -maxsites to apply BEFORE estimating minSNR OR Try running RefAligner -merge -minsites %d -maxsites %d as a subsequent command\n", MinSites,MaxSites);
    fflush(stdout);
  }

  if(bgSNRcnt[0] > 0 && !MapSNR){
    printf("No SNR values in BNX input : -SNRbackground cannot be used without SNR values\n");
    fflush(stdout); exit(1);
  }

  if(outlierMax < 1e+10 && maptype==1){
    printf("-outlierMax %e cannot be used with -maptype 1 (or -i CMAP input, if -maptype is not specified)\n",outlierMax);
    fflush(stdout);exit(1);
  }

  if(!startmaps && !CmapMerge){
    if(num_files > 0){
      printf("WARNING: no maps in %d -i input files\n",num_files);
      fflush(stdout);
    } else {
      printf("No -i input files specified: only supported with -merge\n");
      fflush(stdout);exit(1);
    }
  }

  if(CmapMerge && startmaps && spots_filename[0]){
    printf("Cannot combine -i and -ref inputs with -merge\n");
    fflush(stdout);exit(1);
  }

  if(DEBUG && !spots_filename[0]) assert(Cfirst >= 0);
  if(DEBUG && PairSplit) assert(RefIndex > 0);
  if(DEBUG && PairSplit && !(QueryIndex < num_files)){
    printf("QueryIndex= %d, RefIndex= %d, num_files= %d, startmaps= %d, Cfirst= %d\n",QueryIndex,RefIndex,num_files,startmaps, Cfirst);
    fflush(stdout);
    assert(QueryIndex < num_files);
  }
  if(startmaps > 0){
    input_vfixx_duplicates();

    if(SelectScanId >= 0){
      if(SelectScanId2 >= UniqueScans){
	printf("ERROR: Input BNX has %d UniqueScanId's, -selectScanId %d %d is out of range\n",UniqueScans, SelectScanId+1,SelectScanId2+1);
	fflush(stdout);exit(1);
      }
      int j=0;
      for(int i = 0; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	int scan = pmap->UniqueScanId;
	if(DEBUG) assert(0 <= scan && scan < UniqueScans);
	if(!(SelectScanId <= scan && scan <= SelectScanId2))
	  continue;
	
	if(j < i){
	  Cmap *tmp = Gmap[j];
	  pmap->mapid = j;
	  Gmap[j] = pmap;
	  tmp->mapid = i;
	  Gmap[i] = tmp;
	}
	j++;
      }
      if(VERB && j != nummaps){
	printf("Reduced number of input maps from %d to %d due to -selectScanID %d %d\n",nummaps,j,SelectScanId+1,SelectScanId2+1);
	fflush(stdout);
      }
      if(VERB && j > 0){
	Cmap **maplist = new Cmap*[j];

	double totlen = 0.0, totlen2 = 0.0;
	size_t nsites = 0, nsites2 = 0;
	for(int i= 0; i < j; i++){
	  Cmap *pmap = Gmap[i];
	  maplist[i] = pmap;
	  int M = pmap->numsite[0];
	  totlen += pmap->site[0][M+1];
	  nsites += M;
	  if(colors >= 2){
	    int M2 = pmap->numsite[1];
	    totlen2 += pmap->site[1][M2 + 1];
	    nsites2 += M2;
	  }
	}
	qsort(maplist, j, sizeof(Cmap*), (intcmp*)PmapLenDec);/* sort map pointers in descending order of length */
	double N50Sum= 0.0, N50Len = 0.0;
	for(int i = 0; i < j; i++){
	  Cmap *pmap = maplist[i];
	  int M = pmap->numsite[0];
	  if((N50Sum += pmap->site[0][M+1]) >= totlen * 0.5){
	    N50Len = pmap->site[0][M+1];
	    break;
	  }
	}
	if(DEBUG && totlen > 0.0) assert(N50Len > 0.0);
	
	if(colors >= 2)
	  printf("After applying -selectScanID %d %d : maps=%d, sites=%lu,%lu, length= %0.3f kb (avg= %0.3f kb, label density1= %0.3f /100kb, label density2= %0.3f /100kb, N50= %0.4f Mb):wtime=%0.6f\n",
		 SelectScanId+1, SelectScanId2+1,j,nsites,nsites2,totlen, 0.5 * (totlen+totlen2)/j, nsites*100.0/max(0.001,totlen), nsites2*100.0/max(0.001,totlen2), N50Len*0.001,wtime());
	else
	  printf("After applying -selectScanID %d %d : maps=%d, sites=%lu, length= %0.3f kb (avg= %0.3f kb, label density= %0.3f /100kb, N50= %0.4f Mb):wtime=%0.6f\n",
		 SelectScanId+1, SelectScanId2+1,j,nsites,totlen,totlen/j, nsites*100.0/max(0.001,totlen), N50Len * 0.001,wtime());

	delete [] maplist;
      }
      
      nummaps = j;
    }
  }

  if(VERB>=2 && nummaps > 0){
    printf("After input_vfixx_duplicates:Gmap[0]->numsite[0]=M=%d,Gmap[0]->site[0][M]=%0.3f\n",
	   Gmap[0]->numsite[0],Gmap[0]->site[0][Gmap[0]->numsite[0]]);
    fflush(stdout);
  }
  if(VERB>=2){
    printf("After input_vfixx_duplicates:Cfirst=%d,nummaps=%d,startmaps=%d:\n",Cfirst,nummaps,startmaps);
    for(int i = 0; i < nummaps; i++)
      printf("  Gmap[%d]->id = %lld\n",i, Gmap[i]->id);
    fflush(stdout);
  }
  if(VERB>=2 && MapImageLoc){
    printf("After input_vfixx_duplicates:Cfirst=%d,nummaps=%d,startmaps=%d:\n",Cfirst,nummaps,startmaps);
    for(int i = 0; i < nummaps; i++){
      Cmap *pmap = Gmap[i];
      if(pmap->id != 54001)
	continue;
      printf("Gmap[%d]->id = %lld, numsite[0]=%d\n",i, pmap->id,pmap->numsite[0]);
      for(int J = 1; J <= pmap->numsite[0]; J++)
	printf("\t J=%d:ImageFOV[0][J]=%d,ImageX[0][J]=%0.2f,ImageY[0][J]=%0.2f\n",J,pmap->ImageFOV[0][J],pmap->ImageX[0][J],pmap->ImageY[0][J]);
      fflush(stdout);
    }
  }

  Hapnummaps = nummaps;
  HapQuerymaps = 0;
  if(SplitHmap && nummaps > 0 && maptype == 1){/* check if any of the input Query Cmaps are already Haplotyped : if so generate both Alleles */
    MaxQueryID = 0;
    for(int k = 0; k < Hapnummaps; k++)
      MaxQueryID = max(MaxQueryID, Gmap[k]->id);

    long long MQid = MaxQueryID;
    if(MaxQueryID < 1000000000000000LL){/* round MQid to nearest power of 10 below 1e15, which is around MAXINT64/9223 */
      MQid = 10;
      while(MQid < MaxQueryID)
	MQid *= 10;
    }

    if(VERB/* HERE >=2 */){
      if(MQid != MaxQueryID)
	printf("Using MaxQueryID=%lld -> %lld to shift ID of 2nd Allele of input haplotype Query maps\n",MaxQueryID,MQid);
      else
	printf("Using MaxQueryID=%lld to shift ID of 2nd Allele of input haplotype Query maps\n",MaxQueryID);
      fflush(stdout);
    }

    MaxQueryID = MQid;

    for(int k = 0; k < Hapnummaps; k++)
      if(Gmap[k]->contig){
	split_hapmap(k, Gmap[k]->id + MaxQueryID, Gmap, nummaps, maxmaps, Hapnummaps);
	if(DEBUG) assert(Gmap[nummaps-1]->id == Gmap[k]->id + MaxQueryID);
	if(DEBUG) assert(Gmap[nummaps-1]->Allele == k);
	if(DEBUG) assert(Gmap[k]->Allele == nummaps-1);
	if(DEBUG) assert(Gmap[k]->contig->contig[0].id == Gmap[k]->id);
	if(DEBUG) assert(Gmap[k]->contig->contig[1].id == Gmap[k]->id + MaxQueryID);
	HapQuerymaps++;

	/* discard haplotype information : it will never be needed and could cause problems when generating _q.cmap which cannot be a haplotype map */
	delete Gmap[k]->contig;
	Gmap[k]->contig = NULL;
	Gmap[k]->Allele = Gmap[nummaps-1]->Allele = -1;
      }
      
    if(VERB && nummaps > Hapnummaps){
      printf("Split %d input Haplotyped Query Cmaps into 2 Alleles each : nummaps = %d -> %d (MaxQueryID=%lld)\n",nummaps - Hapnummaps, Hapnummaps, nummaps, MaxQueryID);
      fflush(stdout);
    }
  }

  if(hashcolor > colors){
    printf("WARNING : -hashcolor %d not compatible with -colors %d : Ignoring -hashcolor\n",hashcolor,colors);
    hashcolor = 0;
  }
  if(colors >= 2 && usecolor && hashcolor){
    if(usecolor == hashcolor)
      hashcolor = 0;/* not needed */
    else
      hashcolor = 2;/* If -usecolor 2 -hashcolor 1 :  will need to use 2nd color after colors are swapped due to -usecolor 2 */
  }

  if(!CmapMerge && (cmap_files == num_files || cmap_files == 0) && colors==1 && !usecolor){ /* HERE : check why this segfaults when colors==2 OR when input is mixed BNX & CMAP */
    mapcompact(nummaps,maxmaps,Gmap);

    if(VERB){
      long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
      getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    

      printf("VmRSS= %0.4f, VmHWM= %0.4f : CPU time= %0.6f, wall time= %0.6f\n",VmRSS * 1e-9, VmHWM * 1e-9, mtime(), wtime());      
      fflush(stdout);
    }
  }

  if(VERB>=2 && MapImageLoc){
    printf("After mapcompact:Cfirst=%d,nummaps=%d,startmaps=%d:\n",Cfirst,nummaps,startmaps);
    for(int i = 0; i < nummaps; i++){
      Cmap *pmap = Gmap[i];
      if(pmap->id != 54001)
	continue;
      printf("Gmap[%d]->id = %lld, numsite[0]=%d\n",i, pmap->id,pmap->numsite[0]);
      for(int J = 1; J <= pmap->numsite[0]; J++)
	printf("\t J=%d:ImageFOV[0][J]=%d,ImageX[0][J]=%0.2f,ImageY[0][J]=%0.2f\n",J,pmap->ImageFOV[0][J],pmap->ImageX[0][J],pmap->ImageY[0][J]);
      fflush(stdout);
    }
  }

  // printf("After mapcompact:Sleeping 60 seconds\n");
  //  fflush(stdout);
  //  sleep(60);

  if(usecolor && (colors==2 || usecolor > 1)){
    if(colors == 1 && usecolor > 1){
      printf("ERROR: -usecolor %d is not valid with single color -i input maps\n",usecolor);
      fflush(stdout);exit(1);
    }
    if(DEBUG) assert(colors==2);
    if(DEBUG>=2){
      for(int i = 0; i < startmaps; i++){
	Cmap *pmap = Gmap[i];
	double len1 = pmap->site[0][pmap->numsite[0]+1];
	double len2 = pmap->site[1][pmap->numsite[1]+1];
	if(/* pmap->id == 147 || pmap->id == 14929 || */ (DEBUG && !(fabs(len1-len2) < 0.001))){
	  printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	  for(int c = 0; c < colors; c++)
	    printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.3f\n",
		   c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	  fflush(stdout);
	}
	assert(fabs(len1-len2) < 0.001);
      }
    }

    if(usecolor==2){/* swap 1st and 2nd color information */
      for(int i = 0; i < startmaps; i++)
	Gmap[i]->colorswap(usecolor);
      char *tmp = Nickase[0]; Nickase[0] = Nickase[1]; Nickase[1] = tmp;

      if(VERB/* HERE HERE >=2 */){
	printf("Swapped 1st and 2nd color for %d query maps due to -usecolor 2\n",startmaps);
	printf("Nickase[0]= %s\n",Nickase[0]);
	printf("Nickase[1]= %s\n",Nickase[1]);
	fflush(stdout);
      }
    }

    if(DEBUG>=2 || BIAS_TRACE){
      for(int i = 0; i < startmaps; i++){
	Cmap *pmap = Gmap[i];
	double len1 = pmap->site[0][pmap->numsite[0]+1];
	double len2 = pmap->site[1][pmap->numsite[1]+1];
	if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
	  printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	  for(int c = 0; c < colors; c++)
	    printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.6f\n",
		   c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	  fflush(stdout);
	  assert(fabs(len1-len2) < 0.001);
	}
	if(DEBUG>=2) assert(fabs(len1-len2) < 0.001);
      }
    }

    if(VERB/* HERE HERE >=2 */){
      printf("Changing colors from %d to 1 due to -usecolor %d, colors_set from %d to 1\n",colors,usecolor,colors_set);
      fflush(stdout);
    }
    colors = colors_set = 1;

  } else if(usecolor){/* redundant use of -usecolor 1 with only a single color */
    if(DEBUG) assert(colors==1);
    if(VERB/* HERE HERE >=2 */){
      printf("Ignoring redudant -usecolor %d due to colors=1\n",usecolor);
      fflush(stdout);
    }
    usecolor = 0;/* otherwise it will assume there were originally 2 colors and try to write them out when needed */
  }

  if((DEBUG>=2 || BIAS_TRACE) && (colors >= 2 || usecolor)){
    for(int i = 0; i < startmaps; i++){
      Cmap *pmap = Gmap[i];
      double len1 = pmap->site[0][pmap->numsite[0]+1];
      double len2 = pmap->site[1][pmap->numsite[1]+1];
      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG>=2 && !(fabs(len1-len2) < 0.001))){
	printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	for(int c = 0; c < max(2,colors); c++)
	  printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.6f\n",
		 c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	fflush(stdout);
	assert(fabs(len1-len2) < 0.001);
      }
    }
  }

  if(ScanCorrection > 0 && (BNXVersion != 1 || UniqueScans <= 0)){
    printf("-ScanScaling : ignored since input is not BNX version 1.1 (BNXVersion=%d,UniqueScans=%d)\n",BNXVersion,UniqueScans);
    fflush(stdout);
    ScanCorrection = 0;
  }
  /*  if(ScanCorrection && (CmapStart > 0 && NoSplit != 2)){
    printf("-ScanScaling not implemented with -subset and -nosplit %d\n",NoSplit);
    fflush(stdout);exit(1);
    } */
  if(HashGen && hash_filename && NoSplit <= 0){
    printf("hash table does not support -nosplit 0\n");
    fflush(stdout);exit(1);
  }

  if( bed_filename != 0 ) { //-bed argument was used on command line
    input_bed(bed_filename);
  }

  for(int c = 0; c < colors; c++)
    XmapCount[c] = 0;

  if((XmapStatWrite && !minSNRestimate) || PairMergeRepeat){/* compute statistics : XmapCount, XmapSites, XmapLenth */
    for(int c = 0; c < colors; c++){
      XmapCount[c] = XmapSites[c] = 0;
      XmapLength[c] = XmapGtheta[c] = 0.0;
      for(int i = 0; i < startmaps; i++){
	Cmap *pmap = Gmap[i];
	int M = pmap->numsite[c];
	FLOAT *X = pmap->site[c];
	XmapCount[c]++;
	XmapSites[c] += M;
	XmapLength[c] += X[M+1];
	for(int k= 1; k < M; k++)
	  XmapGtheta[c] += log(X[k+1]-X[k]);
      }
      XmapGtheta[c] /= max(1,XmapSites[c] - XmapCount[c]);
      XmapGtheta[c] = exp(XmapGtheta[c]);
    }
  }

  if(XmapStatWrite && !minSNRestimate){/* output XmapCount,XmapSites and XmapLength to specified file */
    FILE *fp;
    if(VERB){
      printf("Generating -i input map statistics in %s\n",XmapStatWrite);
      fflush(stdout);
    }
    if((fp = fopen(XmapStatWrite,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for writing:errno=%d:%s\n",XmapStatWrite,eno,err);
      fflush(stdout);exit(1);
    }

    /* write out commandline */
    printversion(fp);

    fprintf(fp,"XmapCount=");
    for(int c = 0; c < colors; c++)
      fprintf(fp," %d", XmapCount[c]);
    fprintf(fp,"\nXmapSites=");
    for(int c = 0; c < colors; c++)
      fprintf(fp," %lld", XmapSites[c]);
    fprintf(fp,"\nXmapLength=");
    for(int c = 0; c < colors; c++)
      fprintf(fp," %0.3f", XmapLength[c] /* NEW12 */ * origPixelLen / PixelLen);
    fprintf(fp,"\nXmapGtheta=");
    for(int c = 0; c < colors; c++)
      fprintf(fp," %0.6f", XmapGtheta[c] /* NEW12 */ * origPixelLen / PixelLen);
    fprintf(fp,"\n");
    fclose(fp);
  }

  if(XmapStatRead){/* read in XmapCount,XmapSites and XmapLength from specified file */
    if(colors < 1 && colors > 2){
      printf("colors=%d not supported for -XmapStatRead\n",colors);
      fflush(stdout);exit(1);
    }

    char *filename = XmapStatRead;
    FILE *fp;
    if(VERB){
      printf("Reading -i input map statistics from %s\n",filename);
      fflush(stdout);
    }
    if((fp = fopen(filename,"r"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for reading:errno=%d:%s\n",filename,eno,err);
      fflush(stdout);exit(1);
    }
    
    int linecnt = 1,cnt, found = 0;
    for(;fgets(buf,LINESIZ,fp) != NULL; linecnt++){
      int len = strlen(buf);
      if(len <= 0) continue;
      if(buf[0]=='#' || buf[0] == '\n' || buf[0]== '\r') continue;
      if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	printf("Line too long or not terminated in input file %s on line %d:\n%s\n",filename,linecnt,buf);
	fflush(stdout);exit(1);
      }
      if(VERB>=2){
	printf("Line %d in %s:%s\n",linecnt,filename,buf);
	fflush(stdout);
      }
      if(buf[0] == '#')/* comment lines are ignored */
	continue;

      /* not a comment line */
      switch(colors){
      case 2:
	switch(found){
	case 0:/* XmapCount */
	  if((cnt = sscanf(buf,"XmapCount= %d %d\n",&XmapCount[0],&XmapCount[1])) != 2){
	    printf("error reading line %d in %s : found %d integers (expected 2):\n%s\n",linecnt,filename,cnt,buf);
	    fflush(stdout);exit(1);
	  }
	  break;
	case 1:/* XmapSites */
	  if((cnt = sscanf(buf,"XmapSites= %lld %lld\n",&XmapSites[0],&XmapSites[1])) != 2){
	    printf("error reading line %d in %s : found %d integers (expected 2):\n%s\n",linecnt,filename,cnt,buf);
	    fflush(stdout);exit(1);
	  }
	  break;
	case 2:/* XmapLength */
	  if((cnt = sscanf(buf,"XmapLength= %lf %lf\n",&XmapLength[0],&XmapLength[1])) != 2){
	    printf("error reading line %d in %s : found %d floats (expected 2):\n%s\n",linecnt,filename,cnt,buf);
	    fflush(stdout);exit(1);
	  }
	  break;
	case 3:/* XmapGtheta */
	  if((cnt = sscanf(buf,"XmapGtheta= %lf %lf\n",&XmapGtheta[0],&XmapGtheta[1])) != 2){
	    printf("error reading line %d in %s : found %d floats (expected 2):\n%s\n",linecnt,filename,cnt,buf);
	    fflush(stdout);exit(1);
	  }
	  break;
	default:
	  printf("Unknown extra line %d in %s:\n%s\n",linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	break;
      default:/* 1 color */
	switch(found){
	case 0:/* XmapCount */
	  if((cnt = sscanf(buf,"XmapCount= %d\n",&XmapCount[0])) != 1){
	    printf("error reading line %d in %s : found %d integers (expected 1):\n%s\n",linecnt,filename,cnt,buf);
	    fflush(stdout);exit(1);
	  }
	  break;
	case 1:/* XmapSites */
	  if((cnt = sscanf(buf,"XmapSites= %lld\n",&XmapSites[0])) != 1){
	    printf("error reading line %d in %s : found %d integers (expected 1):\n%s\n",linecnt,filename,cnt,buf);
	    fflush(stdout);exit(1);
	  }
	  break;
	case 2:/* XmapLength */
	  if((cnt = sscanf(buf,"XmapLength= %lf\n",&XmapLength[0])) != 1){
	    printf("error reading line %d in %s : found %d floats (expected 1):\n%s\n",linecnt,filename,cnt,buf);
	    fflush(stdout);exit(1);
	  }
	  break;
	case 3:/* XmapGtheta */
	  if((cnt = sscanf(buf,"XmapGtheta= %lf\n",&XmapGtheta[0])) != 1){
	    printf("error reading line %d in %s : found %d floats (expected 1):\n%s\n",linecnt,filename,cnt,buf);
	    fflush(stdout);exit(1);
	  }
	  break;
	default:
	  printf("Unknown extra line %d in %s:\n%s\n",linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	break;
      }
      found++;
      continue;
    }
    (void)fclose(fp);
    if(found != 4){
      printf("WARNING:Missing lines in %s (only found %d non-comment lines, expected 4 lines)\n",filename,found);
      fflush(stdout);
    }
    if(VERB/* HERE HERE >=2 */){
      printf("XmapCount[0]= %d, XmapSites[0]= %lld, XmapLength[0]= %0.6f, XmapGtheta[0]= %0.6f\n",XmapCount[0],XmapSites[0],XmapLength[0],XmapGtheta[0]);
      fflush(stdout);
    }
  }

  if(bgSNRcnt[0] > 0 && !checkFile(BackgroundSNR)){/* output SNR background distribution with bgSNRcnt samples in file BackgroundSNR and exit */
    FILE *fp;
    if(VERB){
      printf("Generating %s (Background SNR distribution)\n",BackgroundSNR);
      fflush(stdout);
    }
    if((fp = fopen(BackgroundSNR,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for writing:errno=%d:%s\n",BackgroundSNR,eno,err);
      fflush(stdout);exit(1);
    }
    /* write out commandline */
    printversion(fp);

    /* output key header lines */
    fprintf(fp,"# SampleID\t Channel\t SNR\n");

    bgSNR[0] = new double[bgSNRcnt[0]];

    for(int c = 0; c < colors; c++){/* compute SNR background distribution for each channel seperately */
      /* count number of sites in all molecules */
      int totcnt = 0;
      for(int i = 0; i < startmaps; i++){
	Cmap *pmap = Gmap[i];
	if(!pmap->SNR[c])
	  continue;
	totcnt += pmap->numsite[c];
      }
      if(totcnt < bgSNRcnt[0]){
	printf("Only %d SNR values found in %d input maps for color=%d (at least %d SNR values needed)\n",totcnt,startmaps,c+1,bgSNRcnt[0]);
	fflush(stdout);exit(1);
      }
      double *SNRdist = new double[totcnt];
      int SNRcnt = 0;
      for(int i = 0; i < startmaps; i++){
	Cmap *pmap = Gmap[i];
	if(!pmap->SNR[c])
	  continue;
	for(int t = pmap->numsite[c]; t > 0; t--)
	  SNRdist[SNRcnt++] = pmap->SNR[c][t];
      }
      if(DEBUG)assert(SNRcnt == totcnt);

      /* sort all the SNR values */
      qsort(SNRdist,SNRcnt,sizeof(double),(intcmp*)doubleInc);
      
      /* draw bgSNRcnt samples uniformly distributed over all input values */
      for(int i = 0; i < bgSNRcnt[0]; i++){
	double quantile = (i+0.5)/bgSNRcnt[0];
	double samplepoint = quantile * SNRcnt - 0.5;
	int sample = (int)floor(samplepoint);
	if(sample < SNRcnt-1){
	  double wt = samplepoint - sample;
	  bgSNR[0][i] = exp((1.0-wt)*log(SNRdist[sample]) + wt*log(SNRdist[sample+1]));
	} else
	  bgSNR[0][i] = SNRdist[sample];
      }
      
      fprintf(fp,"# Number of Samples= %d out of %d for channel %d\n",bgSNRcnt[0],totcnt,c+1);

      /* write out the sampled SNR values */
      for(int i = 0; i < bgSNRcnt[0]; i++)
	fprintf(fp,"%10d\t%8d\t%0.6f\n",i+1,c+1,bgSNR[0][i]);

      /* free up memory */
      delete [] SNRdist;
    }/* c = 0 .. colors-1 */
    fclose(fp);

    /* free up memory */
    delete [] bgSNR[0];
    bgSNR[0] = NULL;
    bgSNRcnt[0] = 0;

    free(BackgroundSNR);
    BackgroundSNR = 0;
    exit(0);
  } 
  if(BackgroundSNR){/* read in background distribution from file */
    for(int c = 0; c < colors; c++)
      bgSNRcnt[c] = 0;
    char *filename = BackgroundSNR;
    FILE *fp;
    if((fp = fopen(filename,"r"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("Failed to open file %s:errno=%d:%s\n",filename,eno,err);
      fflush(stdout);exit(1);
    }
    int bgSNRmax[MAXCOLOR];

    for(int c = 0; c < colors; c++)
      bgSNRcnt[c] = bgSNRmax[c] = 0;

    int color = -1;
    int linecnt = 1;
    char *pt,*qt;

    for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++) {
      int len = strlen(buf);
      if(len <= 0) continue;
      if(buf[0]=='#' || buf[0] == '\n' || buf[0]== '\r') continue;
      if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	printf("Line too long or not terminated in input file %s on line %d:\n%s\n",filename,linecnt,buf);
	fflush(stdout);exit(1);
      }
      if(buf[0] == '#'){/* most comment lines are ignored except "Number of Samples=" */
	char *key = (char *)"# Number of Samples=";
	if((pt = strstr(buf,key))){
	  pt += strlen(key);
	  int cnt = strtol(pt,&qt,10);
	  if(qt == pt || !isspace(*qt)){
	    printf("Line %d in file %s has invalid Samples count:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  key = (char *)"for channel ";
	  if(!(pt = strstr(qt,key))){
	    printf("Could not find channel id on line %d in file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  pt += strlen(key);
	  color = strtol(pt,&qt,10);
	  if(qt == pt || !(*qt==0 || isspace(*qt)) || color <= 0){
	    printf("Line %d in file %s has invalid channel id:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  color--;
	  if(DEBUG) assert(color >= 0);
	  bgSNRmax[color] = cnt;
	  bgSNR[color] = new double[cnt];

	  if(VERB>=2){
	    printf("Expecting %d SNR values for Channel=%d (see line %d of %s)\n",bgSNRmax[color],color+1,linecnt,filename);
	    fflush(stdout);
	  }
	  continue;
	}
	if(VERB>=2){
	  printf("Skipping comment on line %d of %s:\n%s\n",linecnt,filename,buf);
	  fflush(stdout);
	}
	continue;/* skip all other comment lines */      	
      }

      if(color < 0){
	printf("SNR line encountered on line %d of %s before Number of Samples were specified\n",linecnt,filename);
	printf("color=%d\n",color);
	fflush(stdout);exit(1);
      }

      /* not a comment line */
      pt = buf;
      int index = strtol(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt)) || index <= 0 || index > bgSNRmax[color]){
	printf("Invalid SampleID value (must be between 1 and %d for color=%d) on line %d of %s:\n%s\n",bgSNRmax[color],color+1, linecnt,filename,buf);
	fflush(stdout);exit(1);
      }

      pt = qt;
      int c = strtol(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt)) || c != color+1){
	printf("Invalid Channel=%d (expected %d) on line %d of %s:\n%s\n",c,color+1,linecnt,filename,buf);
	fflush(stdout);exit(1);
      } 
      
      pt = qt;
      double SNR = strtod(pt,&qt);
      if(pt == qt || !(*qt==0 || isspace(*qt)) || SNR < 0.0){
	printf("Invalid SNR value on line %d of %s:\n%s\n",linecnt,filename,buf);
	fflush(stdout);exit(1);
      }
      bgSNR[color][bgSNRcnt[color]++] = SNR;
    }
    
    /* check if all colors have been read in */
    for(int c = 0; c < colors; c++){
      if(bgSNRmax[c] <= 0){
	printf("No SNR values for Channel=%d found in %s\n",c+1,filename);
	fflush(stdout);exit(1);
      }
      if(bgSNRcnt[c] < bgSNRmax[c]){
	printf("Found %d SNR values for Channel=%d in %s (expected %d values)\n",bgSNRcnt[c],c+1,filename,bgSNRmax[c]);
	fflush(stdout);exit(1);
      }
      if(VERB){
	printf("Read in %d SNR value for Channel=%d from %s\n",bgSNRcnt[c],c+1,filename);
	fflush(stdout);
      }
    }
  }

  if(DEBUG>=2 && (colors >= 2 || usecolor)){
    for(int i = 0; i < startmaps; i++){
      Cmap *pmap = Gmap[i];
      double len1 = pmap->site[0][pmap->numsite[0]+1];
      double len2 = pmap->site[1][pmap->numsite[1]+1];
      if(/* pmap->id == 147 || pmap->id == 14929 || */ (DEBUG && !(fabs(len1-len2) < 0.001))){
	printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	for(int c = 0; c < max(2,colors); c++)
	  printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.3f\n",
		 c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	fflush(stdout);
      }
      assert(fabs(len1-len2) < 0.001);
    }
  }

  if(DEBUG) assert(PixelLen > 0.0);
  if(DEBUG && num_files > 0) assert(origPixelLen > 0.0);
  if(origPixelLen <= 0.0) origPixelLen = 0.500;

  if(VERB>=2){
    printf("bpp=%0.6f, PixelLen=%0.6f, origPixelLen=%0.6f\n",bpp*1000.0,PixelLen*1000.0, origPixelLen*1000.0);
    fflush(stdout);
  }

  // -bpp adjustment has been moved to mapcorrect() in input_vfixx.cpp
  assert(fabs(PixelLen - bpp) < PixelLen*0.0001);

  if(BIAS_TRACE){/* print map id BIAS_TRACE_ID */
    for(int m = 0; m < startmaps; m++){
      Cmap *pmap = Gmap[m];
      if(pmap->id == BIAS_TRACE_ID){
	register FLOAT *Y =  pmap->site[0];
	register int N = pmap->numsite[0];
	printf("After bpp correction: id=%lld,N=%d:\n",pmap->id,N);
	for(int J = 1; J <= N+1; J++)
	  printf(" Y[%d]=%0.6f\n",J,Y[J]);
	fflush(stdout);
      }
    }
  }

  if(DEBUG>=2 && (colors >= 2 || usecolor)){
    for(int i = 0; i < startmaps; i++){
      Cmap *pmap = Gmap[i];
      double len1 = pmap->site[0][pmap->numsite[0]+1];
      double len2 = pmap->site[1][pmap->numsite[1]+1];
      if(/* pmap->id == 84 || pmap->id == 147 || pmap->id == 14929 || */ (DEBUG && !(fabs(len1-len2) < 0.001))){
	printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	for(int c = 0; c < max(2,colors); c++)
	  printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.3f\n",
		 c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	fflush(stdout);
      }
      assert(fabs(len1-len2) < 0.001);
    }
  }

  if(ResBins[0] > 0){ /* initialize rawsite[] and apply bias correction to site[] */
    BiasCorrect(Gmap,0,startmaps, 1);
    rawsitemaps = startmaps;
  }

  if(VERB>=2 && (colors >= 2 || usecolor)){
    printf("Called BiasCorrect()\n");
    for(int i = 0; i < startmaps; i++){
      Cmap *pmap = Gmap[i];
      double len1 = pmap->site[0][pmap->numsite[0]+1];
      double len2 = pmap->site[1][pmap->numsite[1]+1];
      if(/* pmap->id == 84 || pmap->id == 147 || pmap->id == 14929 || */ (DEBUG && !(fabs(len1-len2) < 0.001))){
	printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	for(int c = 0; c < max(2,colors); c++)
	  printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.3f\n",
		 c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	fflush(stdout);
      }
      assert(fabs(len1-len2) < 0.001);
    }
  }

  if(VERB>=2 && MapImageLoc){
    printf("After BiasCorrect:Cfirst=%d,nummaps=%d,startmaps=%d:\n",Cfirst,nummaps,startmaps);
    for(int i = 0; i < nummaps; i++){
      Cmap *pmap = Gmap[i];
      if(pmap->id != 54001)
	continue;
      printf("Gmap[%d]->id = %lld, numsite[0]=%d\n",i, pmap->id,pmap->numsite[0]);
      for(int J = 1; J <= pmap->numsite[0]; J++)
	printf("\t J=%d:ImageFOV[0][J]=%d,ImageX[0][J]=%0.2f,ImageY[0][J]=%0.2f\n",
	       J,pmap->ImageFOV[0][J],pmap->ImageX[0][J],pmap->ImageY[0][J]);
      fflush(stdout);
    }
  }

  if(BIAS_TRACE){/* print map id BIAS_TRACE_ID */
    for(int m = 0; m < startmaps; m++){
      Cmap *pmap = Gmap[m];
      if(pmap->id == BIAS_TRACE_ID){
	register FLOAT *Y =  pmap->site[0];
	register int N = pmap->numsite[0];
	printf("After bias correction: id=%lld,N=%d:\n",pmap->id,N);
	for(int J = 1; J <= N+1; J++)
	  printf(" Y[%d]=%0.6f\n",J,Y[J]);
	fflush(stdout);
      }
    }
  }

  if(VERB && num_files >= 1){
    Cmap **maplist = new Cmap*[nummaps];

    if(DEBUG) assert(vfixx_filename[0] != 0);
    double totlen = 0.0, totlen2 = 0.0;
    size_t nsites = 0, nsites2 = 0;
    for(int k = 0; k < nummaps; k++){
      Cmap *pmap = Gmap[k];
      maplist[k] = pmap;
      int M = pmap->numsite[0];
      totlen += pmap->site[0][M+1];
      nsites += M;
      if(colors >= 2){
	int M2 = pmap->numsite[1];
	totlen2 += pmap->site[1][M2 + 1];
	nsites2 += M2;
      }
    }
    qsort(maplist, nummaps, sizeof(Cmap*), (intcmp*)PmapLenDec);/* sort map pointers in descending order of length */
    double N50Sum= 0.0, N50Len = 0.0;
    for(int i = 0; i < nummaps; i++){
      Cmap *pmap = maplist[i];
      int M = pmap->numsite[0];
      if((N50Sum += pmap->site[0][M+1]) >= totlen * 0.5){
	N50Len = pmap->site[0][M+1];
	break;
      }
    }
    if(DEBUG && totlen > 0.0) assert(N50Len > 0.0);
    delete [] maplist;

    if(colors >= 2){
      if(num_files > 1)
	printf("Read %d Maps from %d files %s %s ... : total maps=%d, sites=%lu,%lu length=%0.3f kb (avg=%0.3f kb, label density1 = %0.3f /100kb, label density2 = %0.3f / 100kb, N50= %0.4f Mb):wall time=%0.6f\n",
	       startmaps, num_files, vfixx_filename[0],vfixx_filename[1], nummaps, nsites,nsites2, totlen, totlen/nummaps, nsites*100.0/max(0.001,totlen), nsites2*100.0/max(0.001,totlen2), N50Len*1e-3, wtime());
      else
	printf("Read %d Maps from 1 file %s : total maps=%d, sites=%lu,%lu length=%0.3f kb (avg=%0.3f kb, label density1 = %0.3f /100kb, label density2 = %0.3f / 100 kb, N50= %0.4f Mb):wall time=%0.6f\n",
	       startmaps,vfixx_filename[0],nummaps, nsites, nsites2, totlen, totlen/nummaps, nsites*100.0/max(0.001,totlen), nsites2*100.0/max(0.001,totlen2),N50Len*1e-3,wtime());
    } else {
      if(num_files > 1)
	printf("Read %d Maps from %d files %s %s ... : total maps=%d, sites=%lu, length=%0.3f kb (avg=%0.3f kb, label density = %0.3f /100kb,N50= %0.4f Mb):wall time=%0.6f\n",
	       startmaps, num_files, vfixx_filename[0],vfixx_filename[1], nummaps, nsites, totlen, totlen/nummaps, nsites*100.0/max(0.001,totlen), N50Len*1e-3, wtime());
      else
	printf("Read %d Maps from 1 file %s : total maps=%d, sites=%lu, length=%0.3f kb (avg=%0.3f kb, label density = %0.3f /100kb, N50= %0.4f Mb):wall time=%0.6f\n",
	       startmaps,vfixx_filename[0],nummaps, nsites, totlen, totlen/nummaps, nsites*100.0/max(0.001,totlen), N50Len * 1e-3, wtime());
    }
    //    printf("MapSNR=%d\n",MapSNR);
    fflush(stdout);
  }

  if(PairMerge > 0.0){ /* check that all the input maps are CMAPs (have sitecov and sitecnt) */
    if(spots_filename[0]){
      printf("Cannot use -pairmerge with -ref %s\n", spots_filename[0] );
      fflush(stdout);exit(1);
    }
    for(int i = 0; i < nummaps; i++){
      Cmap *pmap = Gmap[i];
      if(!pmap->sitecov[0] || !pmap->sitecnt[0] || !pmap->siteSD[0]){
	printf("Cannot use -pairmerge without CMAP input: map id=%lld from %s is not CMAP\n",
	       pmap->id,vfixx_filename[pmap->fileid]);
	fflush(stdout);exit(1);
      }
    }
  }

  if(!spots_filename[0]){
    if(BestAlignments > 0){
      printf("Must use -bestalignments with -ref\n");
      fflush(stdout);exit(1);
    }
  }

  if(RepeatShift > 0.0){
    if(spots_filename[0]){
      printf("-Repeat cannot be used with -ref\n");
      fflush(stdout);exit(1);
    }
    if(numbins > 1 && skipcnt > 0){
      printf("-Repeat cannot be used with -partial\n");
      fflush(stdout);exit(1);
    }
    if(PairSplit){
      printf("-Repeat cannot be used with -pairsplit\n");
      fflush(stdout);exit(1);
    }
  }

  if(spots_filename[0]){
    refPixelLen = PixelLen;
    for(int c = 0; c < colors; c++)
      if(res[c] <= 0.0){
	printf("res[%d] parameter must be > 0 (required to condense reference map data and specify data resolution) : changing to 0.001\n",c+1);
	fflush(stdout);
	res[c] = 0.001;
      }
  }

  if(IncID)
    for(int i = 0; i < nummaps; i++)
      Gmap[i]->id += IncID;

  if(CmapMerge){
    if(!nummaps && spots_filename[0]){
      int i = 0;
      merrorexit = 0;
      for(;spots_filename[i]!=NULL;i++) {
	if(VERB>=2)
	  printf( "Reading reference file \"%s\"\n", spots_filename[i]);
	(void)input_spots(spots_filename[i], numrefmaps, maxrefmaps, refmap, i);
      }
      if(merrorexit){
	printf("ERROR:One or more input files %s ... had errors\n",spots_filename[0]);
	fflush(stdout);exit(1);
      }

      if(numrefmaps <= 0){
	printf("ERROR: No reference maps in -ref input files %s ...\n", spots_filename[0]);
	fflush(stdout);exit(1);
      }
      if(VERB){
	double totlen = 0.0;
	size_t nsites = 0;
	for(int k = 0; k < numrefmaps; k++){
	  register Cmap *pmap = refmap[k];
	  int N = pmap->numsite[0];
	  totlen += pmap->site[0][N+1];
	  nsites += N;
	}

	printf("Read %d reference maps from %d input files : total maps=%d, sites=%lu, length=%0.3f kb (avg=%0.3f kb, label density= %0.3f /100kb): wall time= %0.6f\n",
	       numrefmaps,i,numrefmaps,nsites, totlen,totlen/numrefmaps, nsites * 100.0/totlen, wtime());
	fflush(stdout);
      }

      /* copy refmap[0..numrefmaps-1] to maps[] */
      maxmapalloc(numrefmaps, maxmaps, Gmap, 0, 1);
      nummaps = numrefmaps;
      /*      printf("maxmaps=%d,nummaps=%d\n",maxmaps,nummaps);
	      fflush(stdout);*/
      for(int i = 0; i < nummaps; i++){
	Cmap *tmp = Gmap[i];
	Gmap[i] = refmap[i];
	refmap[i] = tmp;
      }
      numrefmaps = 0;
    }

    if(VERB){
      printf("nummaps=%d,CmapStart=%d,CmapEnd= %0.3f,CmapBin=%d,CmapNumBins=%d\n",
	     nummaps,CmapStart,CmapEnd,CmapBin,CmapNumBins);
      fflush(stdout);
    }

    if(nummaps > 0){
      if(numrefmaps > 0){
	printf("Cannot combine -ref and -i with -merge\n");
	fflush(stdout);exit(1);
      }
      if(CMapID > 0){
	if(nummaps != 1){
	  printf("-id %lld invalid since there are multiple (%d) input maps\n",CMapID,nummaps);
	  fflush(stdout);exit(1);
	}
	Gmap[0]->id = CMapID;
      }

      if(RangeRight > RangeLeft){
	if(nummaps != 1){
	  printf("-range %0.3f %0.3f invalid since there are multiple (%d) input maps\n",RangeLeft,RangeRight,nummaps);
	  fflush(stdout);exit(1);
	}
	if(MapSNR){
	  printf("SNR values ignored due with -range\n");
	  fflush(stdout);

	  MapSNR = 0;
	}
	Gmap[0]->trim(RangeLeft,RangeRight);
      }
      if(InversionRight > InversionLeft){
	if(nummaps != 1){
	  printf("-inversion %0.3f %0.3f invalid since there are multiple (%d) input maps\n",InversionLeft,InversionRight,nummaps);
	  fflush(stdout);exit(1);
	}
	double len = Gmap[0]->site[0][Gmap[0]->numsite[0]+1];
	if(InversionLeft >= len){
	  printf("-inversion %0.3f %0.3f invalid since the map length is %0.3f is smaller than the left end of the region\n",InversionLeft,InversionRight,len);
	  fflush(stdout);exit(1);
	}
	if(InversionRight > len)
	  InversionRight = len;
	Gmap[0]->inversion(InversionLeft,InversionRight);
      }
      
      if(GenerateFasta){
	extern char *Nickase[MAXCOLOR];

        size_t max_elen = 0;
	for(int c = 0; c < colors; c++){
	  if(!Nickase[c] || !strcmp(Nickase[c],"unknown")){
	    printf("Nickase unknown for color=%d : cannot generate -fasta file\n",c+1);
	    fflush(stdout);exit(1);
	  }
	  max_elen = max(max_elen, sizeof(Nickase[c]));
	}

	for(int i = 0; i < nummaps; i++){
	  Cmap *pmap = Gmap[i];

	  char filename[BUFSIZ];
	  sprintf(filename,"%s_contig%lld.fasta",output_prefix,pmap->id);
	  FILE *fp = fopen(filename,"w");
	  if(fp == NULL){
	    printf("Cannot write file %s\n",filename);
	    fflush(stdout);exit(1);
	  }

	  int N = pmap->numsite[0];
	  double len = pmap->site[0][N+1] * 1000.0;
	  size_t ilen = floor(len + 0.5);

	  char *seq = new char[ilen + max_elen + 1];
	  for(size_t t = 0; t < ilen + max_elen; t++)
	    seq[t] = 'A';

	  size_t end = ilen + 1;

	  for(int c = 0; c < colors; c++){
	    int Nc = pmap->numsite[c];
	    size_t elen = strlen(Nickase[c]);
	    if(DEBUG) assert(elen <= max_elen);

	    for(int I = 1; I <= Nc; I++){
	      double site = pmap->site[c][I] * 1000.0;
	      int isite = floor(site + 0.5);
	      end = max(end, isite + elen);
	      for(size_t t = 0; t < elen; t++)
		seq[isite + t] = Nickase[c][t];
	    }
	  }
	  
	  if(DEBUG) assert(end <= ilen + max_elen);

	  size_t cnt = 0, lcnt = 0;

	  for(size_t t = 0; t < end; t += 256ul){
	    size_t blen = min(end - t, 256ul);
	    size_t r = fwrite(&seq[t],1,blen,fp);
	    if(r != blen){
	      printf("ERROR: wrote %lu/%lu bytes to %s starting at %lu\n",r,blen,filename,t);
	      fflush(stdout);exit(1);
	    }
	    fprintf(fp,"\n");
	    cnt += blen;
	    lcnt++;
	  }

	  if(VERB/* HERE HERE >=2 */){
	    printf("Wrote %lu bases in %lu lines to %s\n",cnt,lcnt,filename);
	    fflush(stdout);
	  }
	  (void)fclose(fp);
	  delete [] seq;
	}
      }

      if(ForceMerge){
	if(colors != 1){
	  printf("-forcemerge not yet implemented for more than 1 color (Label Channel): colors=%d\n", colors);
	  fflush(stdout);exit(1);
	}

	Cid2mapid *id2mapid = new Cid2mapid[nummaps];
	for(int i = nummaps; --i >= 0;){
	  Cid2mapid *p = &id2mapid[i];
	  Cmap *pmap = Gmap[i];
	  p->id = pmap->id;
	  p->mapid = i;
	  pmap->paired = 0;
	}
	qsort(id2mapid,nummaps,sizeof(Cid2mapid),(intcmp *)idinc);

	char *filename = ForceMerge;
	FILE *fp;
	if((fp = fopen(filename,"r"))==NULL){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("Cannot read file %s (to get -forcemerge specifications):errno=%d:%s\n",filename,eno,err);
	  fflush(stdout);exit(1);
	}
	int linecnt = 1, mergecnt = 0;
	for(;fgets(buf,BUFSIZ,fp) != NULL; linecnt++){
	  int len = strlen(buf);
	  if(len <= 0) continue;
	  if(buf[0]=='#' || buf[0]=='\n' || buf[0]=='\r')continue;
	  if(len >= BUFSIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	    printf("Line too long or not terminated in -if input file %s on line %d:\n%s\n",
		   argv[1],linecnt,buf);
	    fflush(stdout);exit(1);
	  }
	  buf[--len] = '\0';
	  while(len > 0 && isspace(buf[len-1]))
	    buf[--len] = '\0';
	  
	  char *pt = buf, *qt;
	  long long ID1 = strtoll(pt,&qt,10);
	  if(pt==qt || (*qt && !isspace(*qt)) || ID1 < 0){
	    printf("ERROR:Invalid ContigID2 on line %d of %s (must be +ve integer):\n%s\n",linecnt,filename,pt);
	    fflush(stdout);exit(1);
	  }
	  int index1 = findid(ID1,id2mapid,nummaps);/* translate ID1 to index value in Gmap[0..nummaps-1] */
	  if(index1 < 0){
	    printf("ERROR:ContigID1=%lld on line %d of %s not found in %d input maps from %s .. (%d files)\n",ID1,linecnt,filename, nummaps, vfixx_filename[0],num_files+num_reffiles);
	    fflush(stdout);exit(1);
	  }
	  Cmap *Ymap = Gmap[index1];
	  if(Ymap->paired){
	    printf("ERROR:ContigID1=%lld on line %d of %s: This map was previously merged based on line %d\n",ID1,linecnt,filename,Ymap->paired);
	    fflush(stdout);exit(1);
	  }

	  pt = qt;
	  int Orientation1 = strtol(pt,&qt,10);
	  if(pt==qt || (*qt && !isspace(*qt)) || !(Orientation1 ==0 || Orientation1 ==1)){
	    printf("ERROR:Orientation1 on line %d of %s invalid (must be 0 or 1):\n%s\n",linecnt,filename,pt);
	    fflush(stdout);exit(1);
	  }
	  
	  pt = qt;
	  int LabelID1 = strtol(pt,&qt,10);
	  if(pt==qt || (*qt && !isspace(*qt)) || LabelID1 <= 0 || LabelID1 > Ymap->numsite[0]){
	    printf("ERROR:LabelID1 on line %d of %s invalid (must be between 1 and %d for ContigID1=%lld)\n%s\n",linecnt,filename,Ymap->numsite[0],ID1,pt);
	    fflush(stdout);exit(1);
	  }

	  pt = qt;
	  long long ID2 = strtoll(pt,&qt,10);
	  if(pt==qt || (*qt && !isspace(*qt)) || ID2 < 0){
	    printf("ERROR:Invalid ContigID2 on line %d of %s (must be +ve integer):\n%s\n",linecnt,filename,pt);
	    fflush(stdout);exit(1);
	  }
	  if(ID1 == ID2){
	    printf("ERROR:Cannot merge same map with itself : ContigID1=%lld and ContigID2=%lld are same on line %d of %s\n",ID1,ID2,linecnt,filename);
	    fflush(stdout);exit(1);
	  }
	  int index2 = findid(ID2,id2mapid,nummaps);/* translate ID2 to index value in Gmap[0..nummaps-1] */
	  if(index1 < 0){
	    printf("ERROR:ContigID2=%lld on line %d of %s not found in %d input maps from %s .. (%d files)\n",ID2,linecnt,filename, nummaps, vfixx_filename[0],num_files+num_reffiles);
	    fflush(stdout);exit(1);
	  }
	  Cmap *Xmap = Gmap[index2];
	  if(Xmap->paired){
	    printf("ERROR:ContigID2=%lld on line %d of %s: This map was previously merged based on line %d\n",ID2,linecnt,filename,Xmap->paired);
	    fflush(stdout);exit(1);
	  }

	  pt = qt;
	  int Orientation2 = strtol(pt,&qt,10);
	  if(pt==qt || (*qt && !isspace(*qt)) || !(Orientation2 ==0 || Orientation2 ==1)){
	    printf("ERROR:Orientation2 on line %d of %s invalid (must be 0 or 1):\n%s\n",linecnt,filename,pt);
	    fflush(stdout);exit(1);
	  }
	  
	  pt = qt;
	  int LabelID2 = strtol(pt,&qt,10);
	  if(pt==qt || (*qt && !isspace(*qt)) || LabelID2 <= 0 || LabelID2 > Xmap->numsite[0]){
	    printf("ERROR:LabelID2 on line %d of %s invalid (must be between 1 and %d for ContigID2=%lld)\n%s\n",linecnt,filename,Xmap->numsite[0],ID2,pt);
	    fflush(stdout);exit(1);
	  }

	  mergecnt++;

	  Xmap->paired = linecnt;
	  Ymap->paired = linecnt;

	  /* merge Xmap and Ymap and overwrite the lower ID map and mark other map for deletion later with id= -1 */
	  Cmap *tmp = new Cmap(Ymap,Xmap,Orientation1,Orientation2,LabelID1,LabelID2);
	  tmp->id = min(ID1,ID2);
	  Cmap *newmap = (ID1 < ID2) ? Ymap : Xmap;
	  newmap->allfree();
	  *newmap = *tmp;
	  if(LEAKDEBUG){
	    tmp->init();
	    delete tmp;
	  }
	  if(ID1 < ID2)
	    Xmap->id = -1;
	  else
	    Ymap->id = -1;
        }

	/* remove maps with id == -1 (marked for deletion) */
	int j = 0;
	for(int i = 0; i < nummaps; i++){
	  if(Gmap[i]->id < 0)
	    continue;
	  if(j < i){
	    Cmap *tmp = Gmap[j];
	    Gmap[j] = Gmap[i];
	    Gmap[i] = tmp;
	  }
	  j++;
	}
	nummaps = j;
      }

      if(WinBreak > 0 && WinBreakShift > 0.0 && WinBreakLen > 0.0){
	/* first check for largest input contig size */
	double maxlen = 0.0;
	for(int i = 0; i < nummaps; i++){
	  Cmap *pmap = Gmap[i];
	  double len = pmap->site[0][pmap->numsite[0]+1];
	  maxlen = max(len,maxlen);
	}
	int maxPieces = (int)floor((maxlen - WinBreakLen)/WinBreakShift) + 2;
	if(maxPieces > WinBreak){
	  if(VERB){
	    printf("-winbreak : increasing ID multiple from %d to %d due to largest contig len= %0.3f kb\n",WinBreak,maxPieces,maxlen);
	    fflush(stdout);
	  }
	  WinBreak = maxPieces;
	}
	if(VERB/* HERE >=2 */){
	  printf("WinBreak=%d,WinBreakLen=%0.3f,WinBreakShift=%0.3f, maxPieces=%d\n",WinBreak,WinBreakLen,WinBreakShift,maxPieces);
	  fflush(stdout);
	}
	
	for(int i = nummaps; --i >= 0; ){
	  Cmap *pmap = Gmap[i];
	  double len = pmap->site[0][pmap->numsite[0]+1];

	  maxmapalloc(nummaps + WinBreak, maxmaps, Gmap, 0, 1);

	  double maxLeft = len - (WinBreakLen - WinBreakShift);

	  int PieceID = 2;
	  for(double Left = WinBreakShift; Left < maxLeft; Left += WinBreakShift, PieceID++){
	    double Right = min(len, Left + WinBreakLen);
	    if(DEBUG) assert(Right > Left);
	    
	    if(DEBUG>=2) assert(pmap->siteSD[0] && pmap->sitecov[0] && pmap->sitecnt[0]);

	    Cmap *newmap = Gmap[nummaps];
	    if(newmap){
	      newmap->allfree();
	      if(DEBUG>=2) assert(pmap->siteSD[0] && pmap->sitecov[0] && pmap->sitecnt[0]);
	      Cmap *tmp = new Cmap(pmap);
	      if(DEBUG>=2) assert(tmp->siteSD[0] && tmp->sitecov[0] && tmp->sitecnt[0]);
	      *newmap = *tmp;
	      if(DEBUG>=2) assert(newmap->siteSD[0] && newmap->sitecov[0] && newmap->sitecnt[0]);
	      if(LEAKDEBUG){
		tmp->init();
		delete tmp;
	      }
	      if(DEBUG>=2) assert(newmap->siteSD[0] && newmap->sitecov[0] && newmap->sitecnt[0]);
	    } else {
	      Gmap[nummaps] = newmap = new Cmap(pmap);
	      if(DEBUG>=2) assert(newmap->siteSD[0] && newmap->sitecov[0] && newmap->sitecnt[0]);
	    }
	    newmap->id = pmap->id * WinBreak + PieceID;
	    newmap->trim(Left,Right);
	    if(newmap->numsite[0] + (colors >= 2 ? newmap->numsite[1] : 0) >= MinSites)
	      nummaps++;
	  }
	  pmap->trim(0.0,min(len, WinBreakLen));
	  /*	  if(VERB){
	    printf("i=%d:pmap->id = %lld -> %lld, PieceID= %d\n",i, pmap->id, pmap->id * WinBreak + 1, PieceID);
	    fflush(stdout);
	    }*/
	  pmap->id = pmap->id * WinBreak + 1;
	  if(DEBUG) assert(nummaps <= maxmaps);
	  if(pmap->numsite[0] + (colors >= 2 ? pmap->numsite[1] : 0) < MinSites)
	    pmap->id = -1;
	}
	
	/* remove id = -1 */
	int j = 0;
	for(int i = 0; i < nummaps; i++){
	  if(Gmap[i]->id < 0)
	    continue;
	  if(j < i){
	    Cmap *tmp = Gmap[j];
	    Gmap[j] = Gmap[i];
	    Gmap[i] = tmp;
	  }
	  j++;
	}
	nummaps = j;

	qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapIdInc);
      }

      if(BreakIDcnt > 0){
	if(MapSNR){
	  printf("WARNING:SNR values ignored with -break\n");
	  fflush(stdout);

	  MapSNR = 0;
	}
	for(int b = 0; b < BreakIDcnt; b++){
	  long long id = BreakID[b];
	  int cnt = BreakPointCnt[b];
	  /* locate map with id */
	  int i = 0;
	  for(; i < nummaps; i++){
	    Cmap *pmap = Gmap[i];
	    if(pmap->id == id){
	      double Len = pmap->site[0][pmap->numsite[0]+1];
	      maxmapalloc(nummaps + cnt, maxmaps, Gmap, 0, 1);
	      for(int j = 1; j <= cnt; j++){
		double Left = BreakPoint[b][j-1];
		if(DEBUG) assert(0.0 < Left);
		if(Left >= Len){
		  printf("BreakPoint[%d][%d]=%0.3f is >= length %0.3f kb of contig with id=%lld\n",
			 b,j-1,Left,Len,id);
		  fflush(stdout);exit(1);
		}
		double Right = (j < cnt) ? BreakPoint[b][j] : Len;
		if(DEBUG) assert(Left < Right);
		if(Right > Len){
		  printf("BreakPoint[%d][%d]=%0.3f is >= Length %0.3f kb of Gmap[%d] with id=%lld\n",
			 b,j,Right,Len,i,id);
		  fflush(stdout);exit(1);
		}
		Cmap *newmap = Gmap[nummaps];
		if(newmap){
		  newmap->allfree();
		  Cmap *tmp = new Cmap(pmap);
		  *newmap = *tmp;
		  if(LEAKDEBUG){
		    tmp->init();
		    delete tmp;
		  }
		} else
		  Gmap[nummaps] = newmap = new Cmap(pmap);
		newmap->id = id + BreakIDinc * j;
		newmap->trim(Left,Right);
		nummaps++;
	      }
	      pmap->trim(0.0, BreakPoint[b][0]);

	      if(DEBUG) assert(nummaps <= maxmaps);
	      break;
	    }
	  }
	  if(i >= nummaps){
	    printf("WARNING for -break : id=%lld not found in input maps : skipping\n",id);
	    fflush(stdout);
	  }
	}
	if(LEAKDEBUG){
	  for(int b = 0; b < BreakIDcnt; b++)
	    delete [] BreakPoint[b];
	  delete [] BreakPoint; BreakPoint = NULL;
	  delete [] BreakPointCnt; BreakPointCnt = NULL;
	  delete [] BreakID; BreakID = NULL;
	}
      }

      if(CmapSwap){
	if(colors != 2){
	  printf("-swap cannot be used with %d colors (or -usecolor) (must be 2 colors)\n", colors);
	  fflush(stdout);exit(1);
	}
	
	for(int i = 0; i < startmaps; i++)
	  Gmap[i]->colorswap(2);
	char *tmp = Nickase[0]; Nickase[0] = Nickase[1]; Nickase[1] = tmp;
      }
    }
  }

  if(!CmapMerge){
    if(DEBUG && !SplitHmap) assert(startmaps == nummaps);
    if(DEBUG && SplitHmap) assert(startmaps + HapQuerymaps == nummaps);
    totalmaps = nummaps;
  }

  /* find largest UniqueScanId value */
  UniqueScanIdMax = 0;
  UniqueScanIdMin = MASK(31);
  for(int i = 0; i < nummaps; i++){
    UniqueScanIdMax = max(Gmap[i]->UniqueScanId, UniqueScanIdMax);
    UniqueScanIdMin = min(Gmap[i]->UniqueScanId, UniqueScanIdMin);
  }
  UniqueScanIdMax += 1;/* values run from UniqueScanIdMin .. UniqueScanIdMax - 1 */

  if(nummaps > 0){

    if(CmapSort){/* sort Gmap[0..nummaps-1] based on CmapSort value */
      switch(CmapSort){
      case 1:
	printf("Reordering %d maps by size in kb in ascending order\n",nummaps);
	qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapSizeInc);
	break;
      case 2:
	printf("Reordering %d maps by size in kb in descending order\n",nummaps);
	qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapSizeDec);
	break;
      case 3:
	printf("Reordering %d maps by number of sites in ascending order\n",nummaps);
	qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapSitesInc);
	break;
      case 4:
	printf("Reordering %d maps by number of sites in descending order\n",nummaps);
	qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapSitesDec);
	break;
      case 5:/* sort by id in ascending order */
	printf("Reordering %d maps by id in ascending order\n",nummaps);
	qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapIdInc);
	break;
      case 6:
	printf("Randomizing order of %d maps (seed=%d)\n",nummaps,RandSeed);
	srand(RandSeed);
	for(register int i = 0; i < nummaps; i++)
	  Gmap[i]->paired = rand();
	qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapPairedInc);
	break;
      case 7:
	printf("Reordering %d maps in order of RunIndex and ScanNumber\n",nummaps);
	qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapRunIndexScanNumberInc);
	break;
      case 8:
	if(RandSeed > 0){
	  //	  printf("Reordering %d maps by id in ascending order\n",nummaps);
	  //	  qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapIdInc);
	  printf("Randomizing order of %d maps (seed=%d)\n",nummaps,RandSeed);
	  srand(RandSeed);
	  for(register int i = 0; i < nummaps; i++)
	    Gmap[i]->paired = rand();
	  qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapPairedInc);

	  if(VERB>=2){
	    size_t numsites = 0, cnt = 0;
	    double totlen = 0.0;
	    for(int i = max(0,CmapStart - 1); i < (CmapEnd<=0 ? nummaps : CmapEnd); i++){
	      Cmap *pmap = Gmap[i];
	      int M = pmap->numsite[0];
	      double len = pmap->site[0][M+1];
	      numsites += M;
	      totlen += len;
	      cnt++;
	    }
	    if(cnt > 0){
	      printf("After randomizing with Rand Seed=%d: maps=%lu, sites=%lu, length=%0.3f kb (avg=%0.3f kb, label density= %0.3f /100kb)\n",
		     RandSeed, cnt,numsites,totlen,totlen/cnt,numsites * 100.0/totlen);
	      fflush(stdout);
	    }
	  }
	}

#if 1
	if(VERB){
	  printf("Reordering %d maps in order of RunIndex and ScanNumber placing molecules under %0.1f kb last\n", nummaps, CmapSortMinLen);
	  fflush(stdout);
	}
	qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapRunIndexScanNumberIncMinLen);
	
	if(VERB>=2){
	  size_t numsites = 0, cnt = 0;
	  double totlen = 0.0;
	  for(int i = max(0,CmapStart - 1); i < (CmapEnd<=0 ? nummaps : CmapEnd); i++){
	    Cmap *pmap = Gmap[i];
	    int M = pmap->numsite[0];
	    double len = pmap->site[0][M+1];
	    numsites += M;
	    totlen += len;
	    cnt++;
	  }
	  if(cnt > 0){
	    printf("After Sorting by RunIndex/ScanNumber: maps=%lu, sites=%lu, length=%0.3f kb (avg=%0.3f kb, label density= %0.3f /100kb)\n",
		   cnt,numsites,totlen,totlen/cnt,numsites * 100.0/totlen);
	    fflush(stdout);
	  }
	}
#endif

	if(1){
	  if(VERB){
	    printf("Shuffling %d maps across %d RunIndex/ScanNumber bins\n", nummaps, UniqueScanIdMax);
	    fflush(stdout);
	  }

	  /* now shuffle the maps : first locate start of each bin */
	  if(DEBUG && !(UniqueScanIdMin >= 0)){
	    printf("UniqueScanIdMin= %d, UniqueScanIdMax= %d\n",UniqueScanIdMin, UniqueScanIdMax-1);
	    fflush(stdout);
	    assert(UniqueScanIdMin >= 0);
	  }
	  int *UniqueScanIdStart = new int[UniqueScanIdMax];
	  int *UniqueScanIdEnd = new int[UniqueScanIdMax];
	  for(int i = 0; i < UniqueScanIdMax; i++){
	    UniqueScanIdStart[i] = 0;
	    UniqueScanIdEnd[i] = -1;
	  }
	  
	  for(int i = nummaps; --i >= 0; )
	    UniqueScanIdStart[Gmap[i]->UniqueScanId] = i;
	  for(int i = 0; i < nummaps; i++ )
	    UniqueScanIdEnd[Gmap[i]->UniqueScanId] = i;
	  
	  if(DEBUG>=2){
	    int totcnt = 0;
	    for(int i = 0; i < UniqueScanIdMax; i++){
	      if(UniqueScanIdEnd[i] >= 0){
		assert(UniqueScanIdStart[i] <= UniqueScanIdEnd[i] && UniqueScanIdEnd[i] < nummaps);
		totcnt += UniqueScanIdEnd[i] - UniqueScanIdStart[i] + 1;
	      }
	      if(VERB>=2){
		if(UniqueScanIdStart[i] <= UniqueScanIdEnd[i])
		  printf("i=%d/%d:UniqueScanIdStart[i]= %d, UniqueScanIdEnd[i]= %d (cnt=%d), totcnt -> %d, nummaps=%d : RunIndex=%d, ScanNumber=%d\n",
			 i, UniqueScanIdMax, UniqueScanIdStart[i], UniqueScanIdEnd[i], UniqueScanIdEnd[i] - UniqueScanIdStart[i] + 1, 
			 totcnt,nummaps, Gmap[UniqueScanIdStart[i]]->RunIndex, Gmap[UniqueScanIdStart[i]]->ScanNumber);
		else
		  printf("i=%d/%d:UniqueScanIdStart[i]= %d, UniqueScanIdEnd[i]= %d (cnt=%d), totcnt -> %d, nummaps=%d\n",
			 i, UniqueScanIdMax, UniqueScanIdStart[i], UniqueScanIdEnd[i], UniqueScanIdEnd[i] - UniqueScanIdStart[i] + 1, totcnt,nummaps);
		fflush(stdout);
	      }
	    }
	    if(totcnt != nummaps){
	      printf("totcnt=%d, nummaps=%d\n",totcnt,nummaps);
	      fflush(stdout);
	      assert(totcnt == nummaps);
	    }
	  }

	  Cmap **newmap = new Cmap *[nummaps];
	  int nextmap = 0, cnt = 0;
	  while(nextmap < nummaps){
	    cnt++;
	    /* pick next map from each bin and append it to newmap[0 .. nextmap-1] */
	    for(int i = 0; i < UniqueScanIdMax; i++){
	      if(UniqueScanIdStart[i] > UniqueScanIdEnd[i])
		continue;/* no more maps in this bucket */

	      if(DEBUG/* HERE >=2 */ && !(nextmap < nummaps)){
		printf("nextmap=%d,nummaps=%d: cnt=%d,i=%d,UniqueScanIdMax=%d\n",nextmap,nummaps,cnt,i,UniqueScanIdMax);
		fflush(stdout);
		assert(nextmap < nummaps);
	      }

	      newmap[nextmap++] = Gmap[UniqueScanIdStart[i]++];
	    }
	  }
	  
	  if(DEBUG) assert(nextmap == nummaps);
	  if(DEBUG>=2)
	    for(int i = 0; i < UniqueScanIdMax; i++)
	      assert(UniqueScanIdStart[i] == UniqueScanIdEnd[i] + 1);
	  
	  delete [] Gmap;
	  Gmap = newmap;
	
	  delete [] UniqueScanIdStart;
	  delete [] UniqueScanIdEnd;

	  if(VERB>=2){
	    size_t numsites = 0, cnt = 0;
	    double totlen = 0.0;
	    for(int i = max(0,CmapStart - 1); i < (CmapEnd<=0 ? nummaps : CmapEnd); i++){
	      Cmap *pmap = Gmap[i];
	      int M = pmap->numsite[0];
	      double len = pmap->site[0][M+1];
	      numsites += M;
	      totlen += len;
	      cnt++;
	    }
	    if(cnt > 0){
	      printf("After Shuffling by RunIndex/ScanNumber: maps=%lu, sites=%lu, length=%0.3f kb (avg=%0.3f kb, label density= %0.3f /100kb)\n",
		     cnt,numsites,totlen,totlen/cnt,numsites * 100.0/totlen);
	      fflush(stdout);
	    }
	  }
	}


	break;

      default:
	printf("CmapSort=%d: unknown sorting order\n",CmapSort);
	fflush(stdout);exit(1);
      }
      for(int i = 0; i < nummaps; i++)
	Gmap[i]->mapid = i;

      
      fflush(stdout);
    }

    if(DEBUG>=2 && (colors >= 2 || usecolor)){
      for(int i = 0; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	double len1 = pmap->site[0][pmap->numsite[0]+1];
	double len2 = pmap->site[1][pmap->numsite[1]+1];
	if(/* pmap->id == 147 || pmap->id == 14929 || */(DEBUG && !(fabs(len1-len2) < 0.001))){
	  printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	  for(int c = 0; c < max(2,colors); c++)
	    printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.3f\n",
		   c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	  fflush(stdout);
	}
	assert(fabs(len1-len2) < 0.001);
      }
    }

    if(CmapStart >= 1 && CmapEnd < 1.0){
      CmapEnd = floor(CmapEnd * nummaps + 0.5);
      if(DEBUG) assert(0.0 <= CmapEnd && CmapEnd <= (double)nummaps);
      if(CmapEnd < CmapStart)
	CmapEnd = CmapStart;
    }
    if(CmapFrac > 0.0){
      int Frac = floor(CmapFrac * nummaps + 0.5);
      if(VERB && Frac > CmapEnd){
	printf("Increased number of samples from %d to %d due to CmapFrac = %0.6f, nummaps= %d\n", (int)CmapEnd, Frac, CmapFrac, nummaps);
	fflush(stdout);
      }
      CmapEnd = max((int)CmapEnd, Frac);
    }
    if(perCohort > 0){
      int minMaps = perCohort * UniqueScanIdMax;
      if(VERB && minMaps > CmapEnd){
	printf("Increased number of samples from %d to %d due to %d Cohorts (perCohort=%d)\n", (int)CmapEnd, minMaps, UniqueScanIdMax, perCohort);
	fflush(stdout);
      }
      CmapEnd = max((int)CmapEnd, minMaps);
    }

    if(CmapGB > 0.0 && CmapStart < 0){
      CmapStart = 1;
      double GBLen = CmapGB * 1000000;
      double totLen = 0.0;
      int RunIndex = -1, ScanNumber = -1;
      if(CmapSort != 7 && CmapGBrnd)
	CmapGBrnd = 0;

      int i = 0;
      for(; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	double len = pmap->site[0][pmap->numsite[0]+1];
	if(len < CmapMinLen)
	  continue;
	totLen += len;
	if(totLen > GBLen){
	  if(VERB>=3){
	    printf("i=%d/%d: len= %0.3f, totLen = %0.3f, GBlen= %0.3f kb,CmapGBrnd=%d\n", i,nummaps,len, totLen, GBLen,CmapGBrnd);
	    fflush(stdout);
	  }
	  if(!CmapGBrnd)
	    break;
	  if(RunIndex < 0){/* include remaining maps with current RunIndex, ScanNumber */
	    RunIndex = pmap->RunIndex;
	    ScanNumber = pmap->ScanNumber;
	  } else if(pmap->RunIndex != RunIndex || pmap->ScanNumber != ScanNumber)
	    break;
	}
      }
      CmapEnd = min(nummaps,i+1);/* CmapStart and CmapEnd are 1 based */

      if(VERB){
	if(RunIndex >= 0 || CmapEnd < nummaps)
	  printf("Keeping only maps %d .. %d : totLen = %0.1f kb (not counting maps below %0.3f kb), up to RunIndex=%d,ScanNumber=%d\n",CmapStart,(int)CmapEnd, totLen,  CmapMinLen, RunIndex+1,ScanNumber);
	else
	  printf("Keeping all maps %d .. %d : totLen = %0.1f kb (not counting maps below %0.3f kb)\n",CmapStart,(int)CmapEnd, totLen,  CmapMinLen);

	printf("GBLen= %0.3f kb, totLen= %0.3f kb, CmapGBrnd=%d\n",GBLen, totLen, CmapGBrnd);
	fflush(stdout);
      }
    }

    if(CmapEnd > nummaps){
      if(nummaps < CmapStart){
	printf("-subset range %d .. %d is not consistent with number of maps=%d\n",CmapStart,(int)CmapEnd, nummaps);
	fflush(stdout);exit(1);
      }
      CmapEnd = nummaps;
    }

    if(SkipIDcnt > 0){
      if(CmapStart > 0 && CmapEnd >= CmapStart){
	printf("Cannot use -subset with -skipidf\n");
	fflush(stdout);exit(1);
      }	
      if(CmapNumBins > 0){
	printf("Cannot use -subsetbin with -selectid\n");
	fflush(stdout);exit(1);
      }	
      if(VERB){
	printf("Applying -skipidf with %d ids and %d maps\n",SkipIDcnt,nummaps);
	fflush(stdout);
      }
      qsort(SkipID,SkipIDcnt,sizeof(long long),(intcmp*)LLInc);
      qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapIdInc);
      if(VERB>=2){
	printf("Sorted Ids and maps\n");
	fflush(stdout);
      }

      int i = 0, j = 0;
      int skipcnt=0;
      for(int k = 0; k < SkipIDcnt; k++){
	if(DEBUG>=2 && k+1 < SkipIDcnt && !(SkipID[k] <= SkipID[k+1])){
	  printf("k=%d:SkipID[k]=%lld,SkipID[k+1]=%lld\n",k,SkipID[k],SkipID[k+1]);
	  fflush(stdout);
	  assert(SkipID[k] <= SkipID[k+1]);
	}
	if(k+1 < SkipIDcnt && SkipID[k] == SkipID[k+1])
	  continue;

	for(; i < nummaps; i++,j++){
	  if(DEBUG>=2 && i+1 < nummaps) assert(Gmap[i]->id < Gmap[i+1]->id);

	  if(SkipID[k] <= Gmap[i]->id)
	    break;
	  if(j < i){
	    Cmap *p = Gmap[i];
	    Gmap[i] = Gmap[j];
	    Gmap[j] = p;
	  }
	}
	if(VERB>=2){
	  printf("k=%d:SkipID[k]=%lld:i=%d,j=%d:Gmap[i]->id=%lld\n",k,SkipID[k],i,j,Gmap[i]->id);
	  fflush(stdout);
	}
	if(i >= nummaps)
	  break;
	if(SkipID[k] == Gmap[i]->id){
	  skipcnt++;
	  i++;
	}
	continue;
      }
      if(j < i){
	for(;i < nummaps; i++,j++){
	  Cmap *p = Gmap[i];
	  Gmap[i] = Gmap[j];
	  Gmap[j] = p;
	}	
      } else
	j = i = nummaps;
      if(VERB){
	printf("Reduced number of maps from %d to %d due to -skipidf (i=%d,skipcnt=%d)\n",nummaps,j,i,skipcnt);
	fflush(stdout);
      }
      nummaps = j;
      delete [] SkipID;
      SkipID = 0;
    }

    if(SelectIDcnt >= 0){
      if(CmapStart > 0 && CmapEnd >= CmapStart){
	printf("Cannot use -subset with -selectid\n");
	fflush(stdout);exit(1);
      }	
      if(CmapNumBins > 0){
	printf("Cannot use -subsetbin with -selectid\n");
	fflush(stdout);exit(1);
      }	
      if(VERB){
	printf("Applying -selectid with %d ids and %d maps\n",SelectIDcnt,nummaps);
	fflush(stdout);
      }
      //      qsort(SelectID,SelectIDcnt,sizeof(long long),(intcmp*)LLInc);
      qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapIdInc);
      if(VERB>=2){
	printf("Sorted maps by id\n");
	fflush(stdout);
      }

      int i = 0, j = 0;
      for(int k = 0; k < SelectIDcnt; k++){
	for(; i < nummaps; i++)
	  if(SelectID[k] <= Gmap[i]->id)
	    break;
	if(i >= nummaps || SelectID[k] != Gmap[i]->id){
	  if(i >= nummaps)
	    printf("WARNING:Cannot find -selectid %lld in input mapset(nummaps=%d): i=%d,j=%d,k=%d,Gmap[nummaps-1]->id=%lld\n",SelectID[k],nummaps,i,j,k,Gmap[nummaps-1]->id);
	  else if(i > 0 && i < nummaps - 1)
	    printf("WARNING:Cannot find -selectid %lld in input mapset(nummaps=%d): i=%d,j=%d,k=%d,Gmap[i-1,i,i+1]->id=%lld,%lld,%lld\n",SelectID[k],nummaps,i,j,k,Gmap[i-1]->id,Gmap[i]->id,Gmap[i+1]->id);
	  else
	    printf("WARNING:Cannot find -selectid %lld in input mapset(nummaps=%d): i=%d,j=%d,k=%d,Gmap[i]->id=%lld\n",SelectID[k],nummaps,i,j,k,Gmap[i]->id);
	  if(k > 0 && k < SelectIDcnt-1)
	    printf("\t SelectID[k-1,k,k+1] = %lld,%lld,%lld\n",SelectID[k-1],SelectID[k],SelectID[k+1]);
	  fflush(stdout);
	  continue;
	}
	if(j < i){
	  Cmap *p = Gmap[i];
	  Gmap[i] = Gmap[j];
	  Gmap[j] = p;
	}
	j++;
	i++;
      }
      if(VERB && j < nummaps){
	printf("Reduced number of maps from %d to %d due to -selectid\n",nummaps,j);
	fflush(stdout);
      }
      nummaps = j;
      if(LEAKDEBUG){
	delete [] SelectID;
	SelectID = NULL;
      }
    }
  }

  //  printf("CmapBin=%d,CmapNumBins=%d:nummaps=%d\n",CmapBin,CmapNumBins,nummaps);

  if(CmapNumBins > 0){/* divide nummaps into equal sized section and output one (or all) of them */
    if(!CmapMerge){
      printf("WARNING: -subsetbin requires -merge : will terminate after outputing subset files\n");
      fflush(stdout);
    }
    if(CmapStart > 0 && CmapEnd >= CmapStart){
      printf("Cannot use both -subset and -subsetbin\n");
      fflush(stdout);exit(1);
    }
    if(CmapSplit){
      printf("Cannot use both -subsetbin and -split\n");
      fflush(stdout);exit(1);
    }
    char prefix[PATH_MAX];
    if(CmapBin > 0)
      strcpy(prefix,output_prefix);
    int BinStart = (CmapBin <= 0) ? 1 : CmapBin;
    int BinEnd = (CmapBin <= 0) ? CmapNumBins : CmapBin;
    int binsize = nummaps/CmapNumBins;
    /*      if(usecolor==2)
	    for(int i = 0; i < nummaps; i++)
	    Gmap[i]->colorswap(2); */
      
    //    printf("CmapBin=%d,CmapNumBins=%d:BinStart=%d,BinEnd=%d,nummaps=%d,binsize=%d\n",CmapBin,CmapNumBins,BinStart,BinEnd,nummaps,binsize);

    // NOTE : -merge assumes -usecolor is used to select one of 2 input colors to output, so the output will be 1 color
    /*    if(usecolor)
      colors = 2;
    if(usecolor==2)
      for(int i = 0; i < nummaps; i++)
      Gmap[i]->colorswap(usecolor);*/

    for(int Bin = BinStart; Bin <= BinEnd; Bin++){
      CmapStart = 1 + (Bin-1)*binsize;
      int CmapEnd = (Bin == CmapNumBins) ? nummaps : Bin*binsize;
      if(CmapBin <= 0)
	sprintf(prefix,"%s_%d_of_%d",output_prefix,Bin,CmapNumBins);
      if(VERB){
	int sitecnt = 0,mapcnt = 0;
	double totlen = 0.0;
	for(int i = CmapStart-1; i < CmapEnd; i++){
	  Cmap *pmap = Gmap[i];
	  int M = pmap->numsite[0];
	  double len = pmap->site[0][M+1];
	  sitecnt += M;
	  totlen += len;
	  mapcnt++;
	}
	if(CmapStart < CmapEnd)
	  printf("subsetbin: Bin=%d : Keeping only maps %d .. %d (out of %d maps): ids = %lld .. %lld, %d sites(color 1),totlen=%0.3f kb(avg=%0.3f kb, label density = %0.3f /100kb)\n", 
		 Bin, CmapStart, CmapEnd, nummaps, Gmap[CmapStart-1]->id, Gmap[CmapEnd-1]->id, sitecnt, totlen, totlen/mapcnt, sitecnt*100.0/totlen);
	else
	  printf("setsetbin: Bin=%d : Keeping no maps (out of %d maps)\n", Bin, nummaps);
	fflush(stdout);
      }
      if(CmapMergeBnx)
	output_bnx(prefix,Gmap,CmapStart-1,CmapEnd-1, 0, 0,0,-1);/* output maps as a single .bnx file */
      else
	output_cmap(prefix,Gmap,CmapStart-1,CmapEnd-1);/* output maps as a single cmap file */
    }

    goto Lfree;
  }

 
  if(nummaps > 0){
    if(CmapSplit){
      if(!CmapMerge){
	printf("WARNING: -subsetbin requires -merge : will terminate after outputing each contig in seperate files\n");
	fflush(stdout);
      }
      // NOTE : -merge assumes -usecolor is used to select one of 2 input colors to output, so the output will be 1 color
      /*      if(usecolor)
	colors = 2;
      if(usecolor==2)
	for(int i = 0; i < nummaps; i++)
	Gmap[i]->colorswap(usecolor);*/

      char prefix[PATH_MAX];
      for(int i = 0; i < nummaps; i++){
	sprintf(prefix,"%s_contig%lld",output_prefix,Gmap[i]->id);
	if(CmapMergeBnx)
	  output_bnx(prefix,Gmap,i,i, 0, 0,0,-1);/* output maps as a single .bnx file */
	else
	  output_cmap(prefix,Gmap,i,i);/* output maps as a single cmap file */
      }
      
      goto Lfree;
    }

    if(CmapStart > 0){
      int maxcnt = (int)CmapEnd - CmapStart + 1;
      int cnt = maxcnt;
      if(CmapStart > 1 || ScanCorrection){
	cnt = 0;
	for(int i = CmapStart-1; i < startmaps; i++){	    /* swap Gmap[i] with Gmap[cnt] */
	  Cmap *p = Gmap[i];
	  int M = p->numsite[0];
	  double len = p->site[0][M+1];
	  if(ScanCorrection && (len < MinLen || (MaxLen > 0.0 && len > MaxLen)))
	    continue;
	  Gmap[i] = Gmap[cnt];
	  Gmap[cnt] = p;
	  cnt++;
	  if(cnt >= maxcnt)
	    break;
	}	    
	if(DEBUG && ScanCorrection) assert(cnt <= maxcnt);
	if(DEBUG && !ScanCorrection) assert(cnt == maxcnt);
      }
      
      if(VERB){
	int sitecnt = 0, i = CmapStart - 1;
	int mapcnt = 0;
	double totlen = 0.0;
	for(; i < startmaps; i++){
	  Cmap *p = Gmap[i];
	  int M = p->numsite[0];
	  double len = p->site[0][M+1];
	  if(ScanCorrection && (len < MinLen || (MaxLen > 0.0 && len > MaxLen)))
	    continue;
	  sitecnt += M;
	  mapcnt++;
	  totlen += p->site[0][M+1];
	  if(mapcnt >= (int)CmapEnd - CmapStart + 1)
	    break;
	}
	if(DEBUG && !ScanCorrection) assert(i+1 == (int)CmapEnd);

	if(mapcnt > 0) {
	  if(ScanCorrection)
	    printf("Keeping only %d maps from %d .. %d (out of %d maps, MinLen=%0.1f,MaxLen=%0.1f): ids = %lld .. %lld, %d sites(color 1),length=%0.3f kb(avg=%0.3f kb, label density = %0.3f /100kb)\n",
		   mapcnt, CmapStart, i+1, nummaps, MinLen, MaxLen, Gmap[CmapStart-1]->id, Gmap[(int)CmapEnd-1]->id, sitecnt,totlen,totlen/mapcnt,sitecnt*100.0/totlen);
	  else
	    printf("Keeping only maps from %d .. %d (out of %d maps): ids = %lld .. %lld, %d sites(color 1),length=%0.3f kb(avg=%0.3f kb, label density = %0.3f /100kb)\n",
		   CmapStart, (int)CmapEnd, nummaps, Gmap[CmapStart-1]->id, Gmap[(int)CmapEnd-1]->id, sitecnt, totlen,totlen/mapcnt,sitecnt*100.0/totlen);
	} else
	  printf("Keeping no maps (out of %d maps)\n",nummaps);
	fflush(stdout);
	if(DEBUG && !ScanCorrection) assert(mapcnt == (int)CmapEnd - CmapStart + 1);
	if(DEBUG && (CmapStart > 1 || ScanCorrection)) assert(mapcnt == cnt);
      }

      totalmaps = nummaps;/* original number of maps used by -ScanScaling */
      nummaps = cnt;
    }

    if(DEBUG>=2)
      for(int m = 0; m < nummaps; m++)
	Gmap[m]->mapid = m;
  }

  if(simpleRepeatStandalone) {
    if(CmapMerge)
      printf("Warning: skipping -merge because -simpleRepeatStandalone used; use -simpleRepeat instead to do both\n");
    else if(spots_filename[0])
      printf("Warning: skipping reference alignments because -simpleRepeatStandalone used; use -simpleRepeat instead to do both\n");
    else
      printf("Warning: skipping pairwise alignments because -simpleRepeatStandalone used; use -simpleRepeat instead to do both\n");
    findRepeatSimpleMaps(simpleRepeatFilter, simpleRepeatTolerance, simpleRepeatMinEle, false); //last arg is verbose
    goto Lfree;
  }

  if(CmapMerge){
	/*    if(usecolor==2)
      for(int i = 0; i < nummaps; i++)
      Gmap[i]->colorswap(2); */

    if(VERB/* HERE >=2 */ && nummaps > 0){
      Cmap **maplist = new Cmap*[nummaps];

      double totlen = 0.0, totlen2 = 0.0;
      size_t nsites = 0, nsites2 = 0;
      for(int i= 0; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	if(setMask >= 0){
	  for(int c = 0; c < colors; c++){
	    int N = pmap->numsite[c];
	    if(!pmap->Mask[c]){
	      pmap->Mask[c] = new size_t[N+2];
	      memset(pmap->Mask[c], 0, (pmap->numsite[c] + 2) * sizeof(size_t));
	    }
	    pmap->Mask[c][1] = setMask;
	    pmap->Mask[c][N+1] = setMask;
	  }
	}
	if(orMask > 0){
	  for(int c = 0; c < colors; c++){
	    int N = pmap->numsite[c];
	    if(!pmap->Mask[c]){
	      pmap->Mask[c] = new size_t[N+2];
	      memset(pmap->Mask[c], 0, (pmap->numsite[c] + 2) * sizeof(size_t));
	    }
	    pmap->Mask[c][1] |= orMask;
	    pmap->Mask[c][N+1] |= orMask;
	  }
	}
	maplist[i] = pmap;
	int M = pmap->numsite[0];
	totlen += pmap->site[0][M+1];
	nsites += M;
	if(colors >= 2){
	  int M2 = pmap->numsite[1];
	  totlen2 += pmap->site[1][M2 + 1];
	  nsites2 += M2;
	}
      }
      qsort(maplist, nummaps, sizeof(Cmap*), (intcmp*)PmapLenDec);/* sort map pointers in descending order of length */
      double N50Sum= 0.0, N50Len = 0.0;
      for(int i = 0; i < nummaps; i++){
	Cmap *pmap = maplist[i];
	int M = pmap->numsite[0];
	if((N50Sum += pmap->site[0][M+1]) >= totlen * 0.5){
	  N50Len = pmap->site[0][M+1];
	  break;
	}
      }
      if(DEBUG && totlen > 0.0) assert(N50Len > 0.0);
      delete [] maplist;

      if(colors >= 2)
	printf("Final maps=%d, sites=%lu, length= %0.3fkb,%0.3fkb (avg= %0.3f kb, label density1= %0.3f /100kb, label density2= %0.3f /100kb))\n",
	       nummaps,nsites,totlen,totlen2, 0.5 * (totlen+totlen2)/nummaps, nsites*100.0/max(0.001,totlen), nsites2*100.0/max(0.001,totlen2) );
      else
	printf("Final maps=%d, sites=%lu, length= %0.3fkb (avg= %0.3f kb, label density= %0.3f /100kb)\n",
	       nummaps,nsites,totlen,totlen/nummaps, nsites*100.0/max(0.001,totlen));
      fflush(stdout);
    }

    if(CmapMergeBnx)
      output_bnx(output_prefix,Gmap,0,nummaps-1, 0, NULL, NULL,-1);/* output all maps as a single .bnx file */
    else
      output_cmap(output_prefix,Gmap,0,nummaps-1);/* output all maps as a single cmap file */

    if(simpleRepeat)
      findRepeatSimpleMaps(simpleRepeatFilter, simpleRepeatTolerance, simpleRepeatMinEle, false); //last arg is verbose

    goto Lfree;
  }

  if(svcheck){/* construct reference maps in pairs from indels in SV file */
    if(spots_filename[0]){
      printf("Cannot use -svcheck %s with -ref %s ...\n",svcheck,spots_filename[0]);
      fflush(stdout);exit(1);
    }
    if(!nummaps){
      printf("No -i input maps : Cannot perform -svcheck\n");
      fflush(stdout);exit(1);
    }
    if(VERB>=2){
      printf("Calling input_svcheck: nummaps=%d, numrefmaps=%d\n",nummaps, numrefmaps);
      fflush(stdout);
    }
    input_svcheck(svcheck,svcheckKB,numrefmaps,maxrefmaps,refmap);
    if(!numrefmaps){
      printf("No indels in %s : Nothing to do\n",svcheck);
      goto Lfree;
    }
  }

  if(refmap_filename)
    input_refmap(refmap_filename);/* read in reference alignments to reset startloc/endloc */

  if(spots_filename[0]){
    /* read reference map (spots or CMAP format) */
    int i = 0;
    merrorexit = 0;
    for(; spots_filename[i] != NULL; i++){
      if(VERB>=2){
	printf("Reading reference file \"%s\"\n", spots_filename[i]);
	fflush(stdout);
      }
      (void)input_spots(spots_filename[i], numrefmaps, maxrefmaps, refmap, i);
    }
    if(merrorexit){
      printf("ERROR:One or more input files %s ... had errors\n",spots_filename[0]);
      fflush(stdout);exit(1);
    }
    if(numrefmaps <= 0){
      printf("ERROR: No reference maps read in from %s ... \n", spots_filename[0]);
      fflush(stdout);exit(1);
    }
    if(VERB){
      double totlen = 0.0;
      size_t nsites = 0;
      for(int k = 0; k < numrefmaps; k++){
	register Cmap *pmap = refmap[k];
	int N = pmap->numsite[0];
	totlen += pmap->site[0][N+1];
	nsites += N;
      }
      
      printf("Read %d reference maps from %d input files :total maps=%d, sites=%lu, length=%0.3f kb (avg=%0.3f kb, Label Density=%0.3f/100kb)\n",
	     numrefmaps,i,numrefmaps,nsites,totlen,totlen/numrefmaps,nsites * 100.0/totlen);
      fflush(stdout);
    }
    
#if 0
    if(!groupfile){ // cannot be used with group_manifest since it requires that all ref id's be accounted for
      /* filter out any maps with less than AlignedSiteThreshold labels or shorter than AlignedLengthThreshold */
      int j = 0;
      for(int refid = 0; refid < numrefmaps; refid++){
	Cmap *rmap = refmap[refid];
	int N = rmap->numsite[0];
	FLOAT *Y = rmap->site[0];
	if(N < max(1,AlignedSiteThreshold) || Y[N] <= Y[1] + max(res[0]*PixelLen, AlignedLengthThreshold))
	  continue;
	if(j < refid){
	  Cmap *tmp = refmap[j];
	  refmap[j] = refmap[refid];
	  refmap[refid] = tmp;
	}
	j++;
      }
      if(VERB && j < numrefmaps){
	printf("Reduced ref maps from %d to %d due to -A %d -L %0.3f \n",numrefmaps, j, AlignedSiteThreshold, AlignedLengthThreshold);
	fflush(stdout);
      }
      if(j > 0)
	numrefmaps = j;
    }
#endif

    if(FinalSort==4 && spots_filename[0]){/* sort both Gmap[0..nummaps-1] and refmap[0..numrefmaps-1] in descending order of number of labels for better load balancing */
      if(CmapSort != 4){
	if(VERB){
	  printf("sorting %d maps in descending order of labels for better load balancing\n", nummaps);
	  fflush(stdout);
	}
	qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapSitesDec);
      }
      
      // WAS  qsort(refmap,numrefmaps,sizeof(Cmap *),(intcmp*)CmapSitesDec);
    }

    /* need to remember original refmap[] order : also do this for Gmap[], so hashgen() can preserve it */
    for(int i = 0; i < nummaps; i++)
      Gmap[i]->mapid = i;
    for(int i = 0; i < numrefmaps; i++)
      refmap[i]->mapid = i;

    Hapnumrefmaps = numrefmaps;/* remember original number of refmaps : If RefineHap, only these will be refined (& Haplotyped) */
    Hapmaps = 0;

    if(SplitHmap){/* check if any of the input refmaps are already Haplotyped : if so generate both Alleles */
      /* HERE HERE : instead of using SplitCnt > 0 or MaxContigID to shift ID, use the normal Hap CMapID convention to change Allele1's id to N*10+1 and Allele2's id to N*10+2
	             (and N*10 for non-Hap input maps)
	             then fix input_group() to handle the reverse ID mapping when determining the group id and skip changing Allele2's id before calling input_group()
      */
      MaxContigID = 0;
      if(SplitCnt > 0){
	MaxContigID = SplitCnt;

	long long MaxID = 0;
	for(int k = 0; k < Hapnumrefmaps; k++)
	  MaxID = max(MaxID, refmap[k]->id);
	if(MaxID > MaxContigID){
	  printf("WARNING: largest input id = %lld is larger than MaxContigID=%lld (from ID file)\n", MaxID, MaxContigID);
	  MaxContigID = MaxID;
	  if(Refine){
	    printf("WARNING: Input of haplotyped map will have the 2nd Allele assigned an ID that is unique only for this RefAligner job\n");
	    printf("WARNING: Pre-refinement .xmap output filenames may conflict with those of other RefAligner jobs from the same stage\n");         
	    fflush(stdout);
	  }
	}
      } else {
	for(int k = 0; k < Hapnumrefmaps; k++)
	  MaxContigID = max(MaxContigID, refmap[k]->id);

	if(Refine){
	  printf("WARNING: Input of haplotyped map will have the 2nd Allele assigned an ID that is unique only for this RefAligner job\n");
	  printf("WARNING: Pre-refinement .xmap output filenames may conflict with those of other RefAligner jobs from the same stage : Use -contigsplit with largest contig id in ID file to avoid this\n");         
	  fflush(stdout);
	}
      }

      if(VERB/* HERE >=2 */){
	printf("Using MaxContigID=%lld to shift ID of 2nd Allele of input haplotype maps\n",MaxContigID);
	fflush(stdout);
      }

      for(int k = 0; k < Hapnumrefmaps; k++){
	int refid = k;
	Cmap *rmap = refmap[refid];
	if((DEBUG>=2 && rmap->contig && rmap->contig->HapSite[0]) || (VERB>=2 && (refid==3 || refid==2) && rmap->contig && rmap->contig->HapSite[0]) ){
	  int m = rmap->contig->numsite[0];/* real sites in input Hapmap */
	  int *HapSite = rmap->contig->HapSite[0];/* HapSite[0..m+1] for input Hapmap */
	  if((DEBUG && !(HapSite[0]==0 && HapSite[m+1]==0)) || (VERB>=2 && (refid==3 || refid==2))){
	    printf("refid=%d: m=%d,HapSite[0]=%d,HapSite[m+1]=%d\n",
		   k,m,HapSite[0],HapSite[m+1]);
	    fflush(stdout); 
	    assert(HapSite[0]==0 && HapSite[m+1]==0);
	  }
	}

	if(refmap[k]->contig){
	  split_hapmap(k, refmap[k]->id + MaxContigID, refmap, numrefmaps, maxrefmaps, Hapnumrefmaps);
	  if(DEBUG) assert(refmap[numrefmaps-1]->id == refmap[k]->id + MaxContigID);
	  if(DEBUG) assert(refmap[numrefmaps-1]->Allele == k);
	  if(DEBUG) assert(refmap[k]->Allele == numrefmaps-1);
	  if(DEBUG) assert(refmap[k]->contig->contig[0].id == refmap[k]->id);
	  if(DEBUG) assert(refmap[k]->contig->contig[1].id == refmap[k]->id + MaxContigID);
	  if(VERB>=2){
	    printf("refmap[k=%d]=%p (id=%lld) : split off refmap[%d]=%p (id=%lld) as 2nd Allele : refmap[k]->Allele= %d, numrefmaps=%d, MaxContigID= %lld)\n",
		   k,refmap[k], refmap[k]->id, numrefmaps-1, refmap[numrefmaps-1], refmap[numrefmaps-1]->id, refmap[k]->Allele, numrefmaps, MaxContigID);
	    fflush(stdout);
	  }

	  int refid = k;
	  Cmap *rmap = refmap[refid];
	  if((DEBUG>=2 && rmap->contig && rmap->contig->HapSite[0]) || (VERB>=2 && (refid==3 || refid==2) && rmap->contig && rmap->contig->HapSite[0]) ){
	    int m = rmap->contig->numsite[0];/* real sites in input Hapmap */
	    int *HapSite = rmap->contig->HapSite[0];/* HapSite[0..m+1] for input Hapmap */
	    if((DEBUG && !(HapSite[0]==0 && HapSite[m+1]==0)) || (VERB>=2 && (refid==3 || refid==2))){
	      printf("refid=%d: m=%d,HapSite[0]=%d,HapSite[m+1]=%d\n",
		     refid,m,HapSite[0],HapSite[m+1]);
	      fflush(stdout); 
	      assert(HapSite[0]==0 && HapSite[m+1]==0);
	    }
	  }

	  Hapmaps++;
	}

	if(VERB>=2){
	  printf("refmap[t=%d]= %p: id=%lld, Allele=%d, contig=%p\n", k, refmap[k], refmap[k]->id, refmap[k]->Allele, refmap[k]->contig);
	  if(refmap[k]->contig && refmap[k]->contig->contig)
	    printf("\t contig->contig[0].id= %lld, contig->contig[1].id= %lld\n", refmap[k]->contig->contig[0].id, refmap[k]->contig->contig[1].id);
	}
      }

      if(VERB && numrefmaps > Hapnumrefmaps){
	printf("Split %d input Haplotyped Cmaps into 2 Alleles each : numrefmaps = %d -> %d (MaxContigID=%lld)\n",numrefmaps - Hapnumrefmaps, Hapnumrefmaps, numrefmaps, MaxContigID);
	fflush(stdout);
      }

      for(int i = 0; i < numrefmaps; i++)
	refmap[i]->mapid = i;
    }
    if(!RefineHmap)
      Hapnumrefmaps = numrefmaps;/* refine (& Haplotype) all contigs, including (if SplitHmap) each Allele of input Haplotype maps (separately)  */
  }

  if(DEBUG>=3 && (colors >= 2 || usecolor)){
    for(int i = 0; i < startmaps; i++){
      Cmap *pmap = Gmap[i];
      double len1 = pmap->site[0][pmap->numsite[0]+1];
      double len2 = pmap->site[1][pmap->numsite[1]+1];
      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG>=2 && !(fabs(len1-len2) < 0.001))){
	printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	for(int c = 0; c < max(2,colors); c++)
	  printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.6f\n",
		 c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	fflush(stdout);
	assert(fabs(len1-len2) < 0.001);
      }
      
      len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
      len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG>=2 && !(fabs(len1-len2) < 0.001))){
	printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	for(int c = 0; c < max(2,colors); c++)
	  printf("  c=%d:pmap->orignumsite[c]=%d,pmap->origsite[c][pmap->orignumsite[c]+1] = %0.6f\n",
		 c,pmap->orignumsite[c],pmap->origsite[c][pmap->orignumsite[c]+1]);
	fflush(stdout);
	assert(fabs(len1-len2) < 0.001);
      }
    }
  }

  if(VERB>=2 && MapImageLoc){
    printf("Before calling refalign():Cfirst=%d,nummaps=%d,startmaps=%d:\n",Cfirst,nummaps,startmaps);
    for(int i = 0; i < nummaps; i++){
      Cmap *pmap = Gmap[i];
      if(pmap->id != 54001)
	continue;
      printf("Gmap[%d]->id = %lld, numsite[0]=%d\n",i, pmap->id,pmap->numsite[0]);
      for(int J = 1; J <= pmap->numsite[0]; J++)
	printf("\t J=%d:ImageFOV[0][J]=%d,ImageX[0][J]=%0.2f,ImageY[0][J]=%0.2f\n",
	       J,pmap->ImageFOV[0][J],pmap->ImageX[0][J],pmap->ImageY[0][J]);
      fflush(stdout);
    }
  }

  if(VERB>=2){
    printf("minSNRestimate=%d,minSNR=%0.4f..%0.4f,SNRminratio=%0.4f,nummaps=%d,startmaps=%d,totalmaps=%d\n",minSNRestimate,SNRmin,SNRmax,SNRminratio,nummaps,startmaps,totalmaps);
    fflush(stdout);
  }

  if(minSNRestimate && SNRminratio > 0.0){/* raise SNRmin so Xtheta >= SNRminratio * Ylambda */
    double wtime1 = wtime();

    double resKB = res[0] * PixelLen;
    double rawlambda = 0.0;
    int rawYsites = 0;
    for(int refid = 0;refid < numrefmaps;refid++){
      Cmap *rmap = refmap[refid];
      FLOAT *Y = rmap->site[0];
      int N = rmap->numsite[0];
      if(N < max(1,AlignedSiteThreshold) || Y[N] <= Y[1] + max(res[0]*PixelLen, AlignedLengthThreshold))
	continue;
      rawYsites++;
      rawlambda += Y[N+1];
      for(int i = 1; i < N; i++){
	double y = Y[i+1] - Y[i];
	if(y > resKB)
	  rawYsites++;
      }
    }

    if(rawYsites > 0.0){
      rawlambda /= rawYsites;

      int maxmaps = min(100000,totalmaps);

      if(VERB){
	printf("Using %d maps to estimate initial minSNR range (totalmaps=%d,nummaps=%d)\n",maxmaps,totalmaps,nummaps);
	fflush(stdout);
      }
      double rawthetaMin = -1.0;

      double SNR = SNRmin;
      for(; SNR <= SNRmax; SNR += 0.1){
	minSNR[0] = SNR;
	mapcorrectApply(Gmap, 0, maxmaps, maxmaps, MapSNR, 1, numthreads);
	if(ResBins[0] > 0){ /* (re)initialize rawsite[] and apply bias correction to site[] */
	  BiasCorrect(Gmap,0,maxmaps, 1);
	  rawsitemaps = maxmaps;
	}

	double rawtheta = 0.0;
	size_t rawXsites = 0;
	int LoopCnt = 0;

        #pragma omp parallel num_threads(numthreads)
	{
	  int myLoopCnt = 0;
	  double myrawtheta = 0.0;
	  size_t myrawXsites = 0;

          #pragma omp for schedule(static,1024)
	  for(int i= 0; i < maxmaps; i++){
	    if(OMP_DEBUG) myLoopCnt++;
	    Cmap *pmap = Gmap[i];
	    FLOAT *X = pmap->site[0];
	    int N = pmap->numsite[0];
	    myrawtheta += X[N+1];
	    myrawXsites += N;
	  }

	  #pragma omp critical
	  {
	    rawtheta += myrawtheta;
	    rawXsites += myrawXsites;
	    if(OMP_DEBUG) LoopCnt += myLoopCnt;
	  }
        }
	if(OMP_DEBUG) assert(LoopCnt == maxmaps);// check that omp is working on MIC

	rawtheta /= max(1ul,rawXsites);

	if(VERB/* HERE >=2 */){
	  printf("minSNR= %0.2f: Xtheta= %0.3f kb, Ylambda= %0.3f kb\n",minSNR[0], rawtheta, rawlambda);
	  fflush(stdout);
	}

	if(rawtheta >= SNRminratio * rawlambda){
	  if(rawtheta > SNRminratio * rawlambda * 1.10)
	    SNR -= 0.1;
	  if(!rawXsites)
	    SNRmax = SNR;
	  break;
        }

	if(rawthetaMin < 0.0)
	  rawthetaMin = rawtheta;
	else if(rawtheta <= rawthetaMin * 1.01)
	  SNRmin = SNR;
      }

      minSNR[0] = min(max(SNR,SNRmin),SNRmax);
      if(!HashGen && SNRstep > 0.0)
	SNRmin = max(SNRmin, SNR - 8.0);

      mapcorrectApply(Gmap, 0, totalmaps, maxmaps, MapSNR, 1, numthreads);
      if(ResBins[0] > 0){ /* (re)initialize rawsite[] and apply bias correction to site[] */
	BiasCorrect(Gmap,0,startmaps, 1);
	rawsitemaps = startmaps;
      }
      if(VERB && (SNRmin != startSNRmin || minSNR[0] != startSNR)){
	double wtime2 = wtime();
	printf("SNRmin = %0.3f -> %0.3f, SNRmax = %0.3f -> %0.3f, minSNR = %0.3f -> %0.3f (SNRminratio= %0.4f): wall time=%0.6f (cum=%0.6f)\n",
	  startSNRmin, SNRmin, startSNRmax, SNRmax, startSNR, minSNR[0], SNRminratio, wtime2-wtime1,wtime2);
	fflush(stdout);
      }
    }
  }

  if(spots_filename[0] || (svcheck && numrefmaps > 0)){
    extern int RAorigMultiMatches;
    extern int RAorigRefSplit;
    extern double RAorighashdeltaAdjust;
    extern int RAorigHashBest;
    extern double RAorigbiasWT;
    extern int RAorigHashScore;
    extern int RAorigNumResQ;
    extern int RAorigHashMultiMatch;
    extern int RAorigDeltaX;
    extern int RAorigDeltaY;
    extern int RAorigOutlierExtend;
    extern int RAorigOutlierBC;

    //    extern double origFP[2], origFN[2], origSF[2], origSD[2], origSR[2], origSE[2];

    RAorigMultiMatches = MultiMatches;
    RAorigRefSplit = RefSplit;
    RAorighashdeltaAdjust = hashdeltaAdjust;
    RAorigHashBest = HashBest;
    RAorigbiasWT = biasWT;
    RAorigHashScore = HashScore;
    RAorigNumResQ = NumResQ;
    if(VERB/* HERE HERE >=2 */){
      printf("RAorigNumResQ = %d\n",RAorigNumResQ);
      fflush(stdout);
    }
    RAorigHashMultiMatch = HashMultiMatch;
    RAorigDeltaX = DELTA_X;
    RAorigDeltaY = DELTA_Y;
    RAorigOutlierExtend = outlierExtend;
    RAorigOutlierBC = outlierBC;

    if(RefRepeats2 > 1){
      if(VERB){
	/* HERE HERE	if(MultiMatches)
		printf("Changing -MultiMatches from %d to 0\n",MultiMatches);*/
	if(RefSplit)
	  printf("Changing -RefSplit from %d to 0\n",RefSplit);
	if(hashdeltaAdjust > 0)
	  printf("Changing -hashdelta %0.1f %0.1f %0.1f to -hashdelta %0.1f\n",hashdelta,hashdeltaAdjust,hashdeltaLim,hashdelta);
	/*	if(HashBest > 0)
		printf("Changing -hashbest from %d to 0\n",HashBest); */
	if(biasWT > 0.0 && (maptype==1/*NEW29*/ || biaswtDefer))
	  printf("Changing -biaswt from %0.3f to 0\n", biasWT);
	if(HashScore && maptype==1)
	  printf("Changing HashScore (2nd arg of -hashgen) from %d to %d\n", HashScore, HashScore+10);
	/* HERE HERE 	if(NumResQ > (RefRepeats2 > 3 ? 1 : RefRepeats2 > 2 ? 2 : 3))// NEW147
	   printf("Changing NumResQ (last arg of -hashgen) from %d to %d\n", NumResQ, (RefRepeats2 > 3 ? 1 : RefRepeats2 > 2 ? 2 : 3));*/
	/* if(HashMultiMatch)
	   printf("Changing -hashMultiMatch from %d to 0\n", HashMultiMatch); */
	if(DELTA_X != 4)
	  printf("Changing -deltaX from %d to 4\n",DELTA_X);
	if(DELTA_Y != 6)
	  printf("Changing -deltaY from %d to 6\n",DELTA_Y);
	if(outlierExtend)
	  printf("Changing -outlierExtend from %d %d to 0\n",outlierExtend,outlierExtendLimit);
	if(outlierBC)
	  printf("Changing outlierBC from %d to 0\n",outlierBC);
	fflush(stdout);
      }
      // HERE HERE      MultiMatches = 0;
      RefSplit = 0;
      hashdeltaAdjust = 0;
      /*   if(HashBest > 0)
	   HashBest = 0; */
      if(biasWT > 0.0 && (maptype==1/*NEW29*/ || biaswtDefer))
      	biasWT = 0.0;
      if(HashScore && maptype==1)
	HashScore += 10;
      /* HERE HERE       if(NumResQ > (RefRepeats2 > 3 ? 1 : RefRepeats2 > 2 ? 2 : 3))// NEW147
	 NumResQ = (RefRepeats2 > 3) ? 1 : RefRepeats2 > 2 ? 2 : 3; */
      /* HashMultiMatch = 0; */
      DELTA_X = 4;
      DELTA_Y = 6;
      outlierExtend = 0;
      outlierBC = 0;
    } else if(RefRepeats > 1 && biaswtDefer && biasWT > 0.0){
      if(VERB){
	printf("Changing -biaswt from %0.3f to 0\n", biasWT);
	fflush(stdout);
      }
      biasWT = 0.0;
    }

    for(giter2 = 0; giter2 < RefRepeats2; giter2++){
      if(VERB && RefRepeats2 > 1){
	printf("Repeat %d of %d of -M %d\n", giter2+1,RefRepeats2, RefRepeats);
	//	printf("    SR[0]=%0.6f,SR[1]=%0.6f\n",SR[0],SR[1]);
	fflush(stdout);
      }

      if(giter2 > 0 && giter2 < RefRepeats2 - 1){// NEW147
	if(NumResQ < RAorigNumResQ - (RefRepeats2-1-giter2)){
	  if(VERB){
	    printf("giter=%d/%d: Changing -hashgen last arg (NumResQ) from %d to %d\n",giter2,RefRepeats2,NumResQ,RAorigNumResQ - (RefRepeats2-1-giter2));
	    fflush(stdout);
	  }
	  NumResQ = RAorigNumResQ - (RefRepeats2-1-giter2);
	}
      }

      if(giter2 == RefRepeats2 - 1){
	if(VERB){
	  if(MultiMatches != RAorigMultiMatches)
	    printf("giter2=%d/%d: Changing -MultiMatches from %d to %d\n",giter2, RefRepeats2, MultiMatches, RAorigMultiMatches);
	  if(RefSplit != RAorigRefSplit)
	    printf("giter2=%d/%d: Changing -RefSplit from %d to %d\n",giter2, RefRepeats2, RefSplit, RAorigRefSplit);
	  if(hashdeltaAdjust != RAorighashdeltaAdjust)
	    printf("giter2=%d/%d: Changing -hashdelta %0.1f to -hashdelta %0.1f %0.1f %0.1f\n",giter2, RefRepeats2, hashdelta,hashdelta,RAorighashdeltaAdjust,hashdeltaLim);
	  if(HashBest != RAorigHashBest)
	    printf("giter2=%d/%d: Changing -hashbest from %d to %d\n",giter2, RefRepeats2,HashBest,RAorigHashBest);
	  if(biasWT != RAorigbiasWT && (RefRepeats == 1 || !biaswtDefer))
	    printf("giter2=%d/%d: Changing -biaswt from %0.3f to %0.3f\n", giter2, RefRepeats2, biasWT, RAorigbiasWT);
	  if(HashScore != RAorigHashScore)
	    printf("giter2=%d/%d: Changing -hashgen 2nd arg (hashscore) from %d to %d\n", giter2,RefRepeats2,HashScore,RAorigHashScore);
	  if(NumResQ != RAorigNumResQ)
	    printf("giter=%d/%d: Changing -hashgen last arg (NumResQ) from %d to %d\n",giter2,RefRepeats2,NumResQ, RAorigNumResQ);
	  if(HashMultiMatch != RAorigHashMultiMatch)
	    printf("giter=%d/%d: Changing -hashMultiMatch from %d to %d\n", giter2, RefRepeats2,HashMultiMatch, RAorigHashMultiMatch);
	  if(DELTA_X != RAorigDeltaX)
	    printf("giter=%d/%d: Changing -deltaX from %d to %d\n",giter2,RefRepeats2,DELTA_X, RAorigDeltaX);
	  if(DELTA_Y != RAorigDeltaY)
	    printf("giter=%d/%d: Changing -deltaY from %d to %d\n",giter2,RefRepeats2,DELTA_Y, RAorigDeltaY);
	  if(outlierExtend != RAorigOutlierExtend)
	    printf("giter=%d/%d: Changing -outlierExtend from %d to %d %d\n",giter2,RefRepeats2,outlierExtend, RAorigOutlierExtend,outlierExtendLimit);
	  if(outlierBC != RAorigOutlierBC)
	    printf("giter=%d/%d: Changing -OutlierBC from %d to %d\n",giter2,RefRepeats2,outlierBC, RAorigOutlierBC);
	  fflush(stdout);
	}
	MultiMatches = RAorigMultiMatches;
	RefSplit = RAorigRefSplit;
	hashdeltaAdjust = RAorighashdeltaAdjust;
	HashBest = RAorigHashBest;
	if(RefRepeats == 1 || !biaswtDefer)
	  biasWT = RAorigbiasWT;
	HashScore = RAorigHashScore;
	NumResQ = RAorigNumResQ;
	HashMultiMatch = RAorigHashMultiMatch;
	DELTA_X = RAorigDeltaX;
	DELTA_Y = RAorigDeltaY;
	outlierExtend = RAorigOutlierExtend;
	outlierBC = RAorigOutlierBC;
      }

      /* save values before global scan of minSNR or bpp */
      //      int origNumResQ = NumResQ;
      //      int origHashBest = HashBest;

      int origstartmaps = startmaps;
      int orignummaps = nummaps;
      int origtotalmaps = totalmaps;

      int origHash_Bits = Hash_Bits;
      int origMHash_Bits = MHash_Bits;
      int origHHash_Bits = HHash_Bits;

      int origScaleDelta = ScaleDelta;

      if(minSNRestimate && SNRstep > 0.0 && !NoStat && colors==1 && SNRmaxmaps > 0 && nummaps > SNRmaxmaps){/* use a smaller sample size for minSNR[] scan */
	if(VERB){
          printf("Reducing number of maps from %d to %d for global scan of minSNR and disabling -ScaleDelta %0.3f %d\n",nummaps,SNRmaxmaps,ScaleDeltaSize,ScaleDelta);
	  fflush(stdout);
	}
	startmaps = totalmaps = nummaps = SNRmaxmaps;
	Hash_Bits = MHash_Bits = HHash_Bits = 0;
	ScaleDelta = 0;	NumScaleFactor = 1;
      }

      if(minSNRestimate && giter2 == 0 && giter2 < RefRepeats2 - 1){
	if(VERB && NumResQ > 1){
	  printf("giter2=%d/%d:Reducing NumResQ (9th arg of -hashgen) from %d to %d\n",giter2,RefRepeats2,NumResQ,1);
	  fflush(stdout);
	}
	NumResQ = min(1,NumResQ);

	if(HashBest >= 0){
	  if(VERB && HashBest > 0){
	    printf("giter2=%d/%d:Reducing -hashbest from %d to 0\n",giter2,RefRepeats2,HashBest);
	    fflush(stdout);
	  }
	  HashBest = 0;
	}
      }

      if(giter2 > 0){// redundant ?
	if(VERB){
	  if(NumResQ != RAorigNumResQ)
	    printf("giter=%d/%d: Restoring NumResQ (9th arg of -hashgen) from %d to %d\n",giter2,RefRepeats2,NumResQ,RAorigNumResQ);
	  if(HashBest != RAorigHashBest)
	    printf("giter=%d/%d: Restoring -hashbest from %d to %d\n",giter2,RefRepeats2,HashBest,RAorigHashBest);
	  fflush(stdout);
	}
	NumResQ = RAorigNumResQ;
	HashBest = RAorigHashBest;
      }

      if(HashGen){
	if(VERB>=2){
	  printf("Before calling hash_generate():nummaps=%d\n",nummaps);
	  fflush(stdout);
	  dumpmemmap();      
	}
	if(DEBUG) assert(hash_prefix != NULL);

	try {
	  hash_generate(&Gmap[0],nummaps, &refmap[0], numrefmaps, 0, hash_prefix);
	} catch (exception& e){
	  cout << e.what() << endl;
	  printf("hash_generate threw an exception\n");
	  dumpmemmap();
	  fflush(stdout);
	  assert(0);
	}
	if(!hash_filename)
	  goto Lfree;
      }

      if(DEBUG>=3 && (colors >= 2 || usecolor)){
	for(int i = 0; i < startmaps; i++){
	  Cmap *pmap = Gmap[i];
	  double len1 = pmap->site[0][pmap->numsite[0]+1];
	  double len2 = pmap->site[1][pmap->numsite[1]+1];
	  if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG>=2 && !(fabs(len1-len2) < 0.001))){
	    printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	    for(int c = 0; c < max(2,colors); c++)
	      printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.6f\n",
		     c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	    fflush(stdout);
	    assert(fabs(len1-len2) < 0.001);
	  }

	  len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	  len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	  if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG>=2 && !(fabs(len1-len2) < 0.001))){
	    printf("i=%d/%d:colors=%d:pmap->id=%lld:\n",i,startmaps,colors,pmap->id);
	    for(int c = 0; c < max(2,colors); c++)
	      printf("  c=%d:pmap->orignumsite[c]=%d,pmap->origsite[c][pmap->orignumsite[c]+1] = %0.6f\n",
		     c,pmap->orignumsite[c],pmap->origsite[c][pmap->orignumsite[c]+1]);
	    fflush(stdout);
	    assert(fabs(len1-len2) < 0.001);
	  }
	}
      }

      /* compute best alignment of each map with reference */
      try {
	if(minSNRestimate && SNRstep > 0.0 && !NoStat && colors==1){/* first try a few values of minSNR[] to find a good starting value */

	  double Sstep = (giter2 <= 0) ? SNRstep : SNRstep * ((giter2 <= 1) ? 0.3 : (giter2 <= 2) ? 0.1 : 0.05);
	  double Smin = (giter2 <= 0) ? SNRmin : max(SNRmin, minSNR[0] - ((giter2 >= 3) ? 3 : 5) * Sstep);
	  double Smax = (giter2 <= 0) ? SNRmax : min(SNRmax, minSNR[0] + ((giter2 >= 3) ? 3 : 5) * Sstep);
	  if(VERB){
	    printf("giter=%d:Performing global scan of minSNR from %0.3f to %0.3f in steps of %0.3f using %d maps\n",giter2,Smin,Smax,Sstep,nummaps);
	    fflush(stdout);
	  }

	  double origSNR = minSNR[0];
	  double bestSNR = minSNR[0];
	  double bestLL = -1e+300;

	  double startPixelLen = PixelLen;

	  if(BppSteps && giter2 == 0){/* first try initial minSNR to see if initial bpp needs to be adjusted */
	    if(VERB/* >=2 */){
	      printf("Performing global scan of bpp with minSNR[0] = %0.3f\n", minSNR[0]);
	      fflush(stdout);
	    }
	    double SNRtotLL,NPixelLen;
	    refalign_allpairs(nummaps,numrefmaps,&SNRtotLL, &NPixelLen);
	    
	    if(VERB/* HERE >=2 */){
	      printf("At bpp=%0.2f -> %0.2f : LL= %0.6f (BppMincnt= %d)\n", PixelLen * 1000.0, NPixelLen * 1000.0, SNRtotLL, BppMincnt);
	      fflush(stdout);
	    }

	    //	    double bestPixelLen = PixelLen;
	    double besttotLL = SNRtotLL;
	    double bestNPixelLen = NPixelLen;

	    if(SNRtotLL < BppMincnt){/* try alternate BPP values */
	      double startPixelLen = PixelLen;

	      for(int BppCnt = 1; BppCnt <= BppSteps; BppCnt++){
		for(int dir = -1; dir <= 1; dir += 2){
		  double nPixelLen = startPixelLen + BppStepsiz * BppCnt * dir * 0.001;
		  double C = nPixelLen / PixelLen;
		
		  if(C != 1.0){
		    if(VERB/* HERE >=2 */){
		      printf("    Rescaling all map site[] and rawsite[] values  by C=%0.12f (reflecting change in bpp to %0.2f)\n",C, C * PixelLen * 1000.0);
		      fflush(stdout);
		    }

		    int origcolors = colors;
		    if(usecolor)
		      colors = 2;

		    #pragma omp parallel for num_threads(numthreads) schedule(dynamic,256) if(numthreads > 1)
		    for(int i= 0; i < origtotalmaps; i++){
		      Cmap *pmap = Gmap[i];
		      if((DEBUG || (BIAS_TRACE && pmap->id==BIAS_TRACE_ID)) && colors>=2 && !(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-4)){
		        #pragma omp critical
		        {
		          printf("id=%lld:numsite[0,1]=%d,%d, pmap->site[0][numsite[0]+1]= %0.6f, pmap->site[1][numsite[1]+1]= %0.6f: Before rescaling by C= %0.12f\n",
		            pmap->id,pmap->numsite[0],pmap->numsite[1],pmap->site[0][pmap->numsite[0]+1],pmap->site[1][pmap->numsite[1]+1],C);
   		          fflush(stdout);
			  assert(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-4);	    
		        }
		      }
		      for(int c = 0; c < colors; c++){
			FLOAT *X = pmap->site[c];
			int M = pmap->numsite[c];
			for(int j= M+1; j > 0; j--)
			  X[j] *= C;
		      }
		      if((DEBUG || (BIAS_TRACE && pmap->id==BIAS_TRACE_ID)) && colors>=2 && !(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-4)){
		        #pragma omp critical
		        {
		          printf("id=%lld:numsite[0,1]=%d,%d, pmap->site[0][numsite[0]+1]= %0.6f, pmap->site[1][numsite[1]+1]= %0.6f: After rescaling by C= %0.12f\n",
		            pmap->id,pmap->numsite[0],pmap->numsite[1],pmap->site[0][pmap->numsite[0]+1],pmap->site[1][pmap->numsite[1]+1],C);
   		          fflush(stdout);
			  assert(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-4);	    
		        }
		      }
		    }

		    if(maxresbias > mres * 0.5 && rawsitemaps > 0){/* also rescale rawsite[]. NOTE : applied to all maps including -subset maps (if -ScanScaling) and split maps */

                      #pragma omp parallel for num_threads(numthreads) schedule(dynamic,256) if(numthreads > 1)
		      for(int i= 0; i < rawsitemaps; i++){
			Cmap *pmap = Gmap[i];
			if(DEBUG && colors>=2) assert(fabs(pmap->rawsite[0][pmap->numsite[0]+1] - pmap->rawsite[1][pmap->numsite[1]+1]) < 1e-4);
			
			for(int c = 0; c < colors; c++){
			  FLOAT *X = pmap->rawsite[c];
			  int M = pmap->numsite[c];
			  for(int j= M+1; j > 0; j--)
			    X[j] *= C;
			}

			if(DEBUG && colors>=2) assert(fabs(pmap->rawsite[0][pmap->numsite[0]+1] - pmap->rawsite[1][pmap->numsite[1]+1]) < 1e-4);
		      }

		      /* also scale resbias[c][] values */
		      for(int c = 0; c < colors; c++){
			for(int Bin = 0; Bin <= ResBins[c]; Bin++){
			  resbias[c][Bin] *= C;
			  resbiasX[c][Bin] *= C;
			}
		      }
		      maxresbias *= C;

		      if(maxresbias <= mres * 0.5)
			maxresbias = mres * 0.5 + 1e-6;/* Fudge : should never happen */
		      if(DEBUG>=2)
			for(int c = 0; c < colors; c++)
			  assert(resbiasX[c][ResBins[c]] <= maxresbias);
		    }

		    colors = origcolors;
		    
		    PixelLen *= C;
		  
		    if(HashGen){
		      if(VERB>=2){
			printf("Before calling hash_generate():nummaps=%d\n",nummaps);
			fflush(stdout);
			dumpmemmap();      
		      }
		    
		      if(DEBUG) assert(hash_prefix != NULL);

		      try {
			hash_generate(&Gmap[0],nummaps, &refmap[0], numrefmaps, 0, hash_prefix);
		      } catch (exception& e){
			cout << e.what() << endl;
			printf("hash_generate threw an exception\n");
			dumpmemmap();
			fflush(stdout);
			assert(0);
		      }
		    }
		  }

		  refalign_allpairs(nummaps,numrefmaps,&SNRtotLL, &NPixelLen);
		  if(SNRtotLL > besttotLL){
		    //		    bestPixelLen = PixelLen;
		    bestNPixelLen = NPixelLen;
		    besttotLL = SNRtotLL;
		  }
		  if(VERB/* HERE >=2 */){
		    printf("At bpp=%0.2f -> %0.2f : LL= %0.6f (best bpp=%0.2f, besttotLL= %0.6f)\n", PixelLen * 1000.0,NPixelLen * 1000.0, SNRtotLL, bestNPixelLen * 1000.0, besttotLL);
		    fflush(stdout);
		  }
		  if(besttotLL >= BppMincnt)
		    break;
		}// for dir = -1 .. 1
		if(besttotLL >= BppMincnt)
		  break;
	      }// BppCnt = 1 .. BppSteps

	      if(VERB && besttotLL < BppMincnt){
		printf("WARNING: BPP value scan failed to find a valid value : please try a larger scan range. Will use bpp= %0.2f\n", bestNPixelLen * 1000.0);
		fflush(stdout);
	      }

	      if(PixelLen != bestNPixelLen){
		double C = bestNPixelLen / PixelLen;

		if(C != 1.0){
		  if(VERB/* HERE >=2 */){
		    printf("    Rescaling all map site[] and rawsite[] values  by C=%0.12f (reflecting change in bpp to %0.2f) and computing initial LL\n",C, C * PixelLen * 1000.0);
		    //		    printf("\t origtotalmaps= %d, Gmap[%d]->site[0] = x%p\n",origtotalmaps, 1003815, Gmap[1003815]->site[0]);
		    fflush(stdout);
		  }

		  int origcolors = colors;
		  if(usecolor)
		    colors = 2;

                  #pragma omp parallel for num_threads(numthreads) schedule(dynamic,256) if(numthreads > 1)
		  for(int i= 0; i < origtotalmaps; i++){
		    Cmap *pmap = Gmap[i];
		    if(DEBUG && colors>=2) assert(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-4);
		    for(int c = 0; c < colors; c++){
		      FLOAT *X = pmap->site[c];
		      int M = pmap->numsite[c];
		      for(int j= M+1; j > 0; j--)
			X[j] *= C;
		    }
		    if(DEBUG && colors>=2) assert(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-4);	    
		  }

		  if(maxresbias > mres * 0.5 && rawsitemaps > 0){/* also rescale rawsite[]. NOTE : applied to all maps including -subset maps (if -ScanScaling) and split maps */

                    #pragma omp parallel for num_threads(numthreads) schedule(dynamic,256) if(numthreads > 1)
		    for(int i= 0; i < rawsitemaps; i++){
		      Cmap *pmap = Gmap[i];
		      if(DEBUG && colors>=2) assert(fabs(pmap->rawsite[0][pmap->numsite[0]+1] - pmap->rawsite[1][pmap->numsite[1]+1]) < 1e-4);
			
		      for(int c = 0; c < colors; c++){
			FLOAT *X = pmap->rawsite[c];
			int M = pmap->numsite[c];
			for(int j= M+1; j > 0; j--)
			  X[j] *= C;
		      }
		      
		      if(DEBUG && colors>=2) assert(fabs(pmap->rawsite[0][pmap->numsite[0]+1] - pmap->rawsite[1][pmap->numsite[1]+1]) < 1e-4);
		    }

		    /* also scale resbias[c][] values */
		    for(int c = 0; c < colors; c++){
		      for(int Bin = 0; Bin <= ResBins[c]; Bin++){
			resbias[c][Bin] *= C;
			resbiasX[c][Bin] *= C;
		      }
		    }
		    maxresbias *= C;
		    
		    if(maxresbias <= mres * 0.5)
		      maxresbias = mres * 0.5 + 1e-6;/* Fudge : should never happen */
		    if(DEBUG>=2)
		      for(int c = 0; c < colors; c++)
			assert(resbiasX[c][ResBins[c]] <= maxresbias);
		  } 

		  colors = origcolors;
		  
		  PixelLen *= C;

		  if(minSNRestimate && SNRminratio > 0.0){/* recompute SNRmin to account for new bpp value */
		    SNRmin = startSNRmin;
		    minSNR[0] = startSNR;

		    double rawthetaMin = -1.0;

		    double wtime1 = wtime();

		    double resKB = res[0] * PixelLen;
		    double rawlambda = 0.0;
		    int rawYsites = 0;
		    for(int refid = 0;refid < numrefmaps;refid++){
		      Cmap *rmap = refmap[refid];
		      FLOAT *Y = rmap->site[0];
		      int N = rmap->numsite[0];
		      rawYsites++;
		      rawlambda += Y[N+1];
		      for(int i = 1; i < N; i++){
			double y = Y[i+1] - Y[i];
			if(y > resKB)
			  rawYsites++;
		      }
		    }

		    if(rawYsites > 0.0){
		      rawlambda /= rawYsites;
		      double SNR = SNRmin;
		      for(; SNR <= SNRmax; SNR += 0.1){
			minSNR[0] = SNR;
			mapcorrectApply(Gmap, 0, totalmaps, maxmaps, MapSNR, 1, numthreads);
			if(ResBins[0] > 0){ /* (re)initialize rawsite[] and apply bias correction to site[] */
			  BiasCorrect(Gmap,0,startmaps, 1);
			  rawsitemaps = startmaps;
			}

			double rawtheta = 0.0;
			size_t rawXsites = 0;
			int LoopCnt = 0;

                        #pragma omp parallel num_threads(numthreads)
			{
			  double myrawtheta = 0.0;
			  size_t myrawXsites = 0;
			  int myLoopCnt = 0;

                          #pragma omp for schedule(static,1024)
			  for(int i= 0; i < totalmaps; i++){
			    if(OMP_DEBUG) myLoopCnt++;
			    Cmap *pmap = Gmap[i];
			    FLOAT *X = pmap->site[0];
			    int N = pmap->numsite[0];
			    myrawtheta += X[N+1];
			    myrawXsites += N;
			  }

                          #pragma omp critical
			  {
			    rawtheta += myrawtheta;
			    rawXsites += myrawXsites;
			    if(OMP_DEBUG) LoopCnt += myLoopCnt;
		          }
		        }
			if(OMP_DEBUG) assert(LoopCnt == totalmaps);// check that OMP is working correctly

			if(!rawXsites)
			  break;
			rawtheta /= rawXsites;
			if(VERB/* HERE >=2 */){
			  printf("minSNR= %0.2f: Xtheta= %0.3f kb, Ylambda= %0.3f kb\n",minSNR[0], rawtheta, rawlambda);
			  fflush(stdout);
		        }
			if(rawtheta >= SNRminratio * rawlambda)
			  break;

			if(rawthetaMin < 0.0)
			  rawthetaMin = rawtheta;
			else if(rawtheta <= rawthetaMin * 1.000001 /* WAS22 1.01 */)
			  SNRmin = SNR;
		      }

		      minSNR[0] = origSNR; // Preserve minSNR value used with bppScan // WAS max(SNR,startSNR);

		      if(!HashGen && SNRstep > 0.0)
			SNRmin = min(max(SNRmin,origSNR - 8.0), origSNR);

		      double origSNRmax = SNRmax;
		      SNRmax = max(SNRmax, max(SNRmin,origSNR));// NEW22
		      minSNR[0] = min(max(minSNR[0], SNRmin),SNRmax);// NEW22

		      mapcorrectApply(Gmap, 0, totalmaps, maxmaps, MapSNR, 1, numthreads);
		      if(ResBins[0] > 0){ /* (re)initialize rawsite[] and apply bias correction to site[] */
			BiasCorrect(Gmap,0,startmaps, 1);
			rawsitemaps = startmaps;
		      }
		      if(VERB && (SNRmin != startSNRmin || minSNR[0] != startSNR)){
			double wtime2 = wtime();
			printf("SNRmin = %0.3f -> %0.3f, SNRmax = %0.3f -> %0.3f, minSNR[0] = %0.3f -> %0.3f (SNRminratio= %0.4f): wall time=%0.6f (cum=%0.6f)\n",
			  startSNRmin, SNRmin, origSNR, origSNRmax, SNRmax, minSNR[0], SNRminratio, wtime2-wtime1,wtime2);
			fflush(stdout);
		      }
		    }

		    Sstep = (giter2 <= 0) ? SNRstep : SNRstep * ((giter2 <= 1) ? 0.3 : (giter2 <= 2) ? 0.1 : 0.05);
		    Smin = (giter2 <= 0) ? SNRmin : max(SNRmin, minSNR[0] - ((giter2 >= 3) ? 3 : 5) * Sstep);
		    Smax = (giter2 <= 0) ? SNRmax : min(SNRmax, minSNR[0] + ((giter2 >= 3) ? 3 : 5) * Sstep);
		    if(VERB){
		      printf("giter=%d:Performing global scan of minSNR from %0.3f to %0.3f in steps of %0.3f using %d maps\n",giter2,Smin,Smax,Sstep,nummaps);
		      fflush(stdout);
		    }

		    origSNR = minSNR[0];
		    bestSNR = minSNR[0];
  		  }// if(SNRestimate && SNRminratio > 0.0)

		  if(HashGen){
		    if(VERB>=2){
		      printf("Before calling hash_generate():nummaps=%d\n",nummaps);
		      fflush(stdout);
		      dumpmemmap();      
		    }
		    
		    if(DEBUG) assert(hash_prefix != NULL);

		    try {
		      hash_generate(&Gmap[0],nummaps, &refmap[0], numrefmaps, 0, hash_prefix);
		    } catch (exception& e){
		      cout << e.what() << endl;
		      printf("hash_generate threw an exception\n");
		      dumpmemmap();
		      fflush(stdout);
		      assert(0);
		    }
		  }
		}

		BppSteps = 0;// disable any subsequent use of -bppScan
		refalign_allpairs(nummaps,numrefmaps,&SNRtotLL);

		if(SNRtotLL > bestLL){
		  bestSNR = minSNR[0];
		  bestLL = SNRtotLL;
		}
		if(VERB){
		  printf("minSNR=%0.3f, SNRtotLL=%0.6f: bestSNR=%0.3f,bestLL=%0.6f\n",minSNR[0],SNRtotLL,bestSNR,bestLL);
		  fflush(stdout);
		}
	      }// if(PixelLen != bestNPixelLen)

	      if(DEBUG) assert(fabs(PixelLen - bestNPixelLen) < 1e-6);
   	    }/* if(SNRtotLL < BppMincnt) */ else {
	      BppSteps = 0;// disable any subsequent use of BppMincnt 
	      if(DEBUG && besttotLL >= BppMincnt && !(SNRtotLL >= besttotLL - 0.001)){
		printf("WARNING: besttotLL= %0.6f, BppMincnt= %d, SNRtotLL= %0.6f\n", besttotLL, BppMincnt, SNRtotLL);
		fflush(stdout);
		//	      assert( SNRtotLL >= besttotLL - 0.001);
	      }

	      bestLL = besttotLL;
	    }
          }// if(BppSteps ... )

	  for(double SNR = Smin; SNR - Sstep < Smax - 0.099; SNR += Sstep){
	    if(BppSteps && giter2 == 0 && fabs(SNR - origSNR) < 0.001)
	      continue;
	    if(BppSteps && giter2 == 0 && fabs(SNR - bestSNR) < 0.001)
	      continue;

	    minSNR[0] = SNR;
	    if(VERB/* >=2 */){
	      printf("minSNR[0] -> %0.3f\n", SNR);
	      fflush(stdout);
	    }

	    mapcorrectApply(Gmap, 0, totalmaps, maxmaps, MapSNR, 1, numthreads);
	    if(ResBins[0] > 0){ /* (re)initialize rawsite[] and apply bias correction to site[] */
	      BiasCorrect(Gmap,0,startmaps, 1);
	      rawsitemaps = startmaps;
	    }

	    double SNRtotLL;
	    refalign_allpairs(nummaps,numrefmaps,&SNRtotLL);

	    if(SNRtotLL > bestLL){
	      bestSNR = SNR;
	      bestLL = SNRtotLL;
	    }
	    if(VERB){
	      printf("minSNR=%0.3f, SNRtotLL=%0.6f: bestSNR=%0.3f,bestLL=%0.6f\n",SNR,SNRtotLL,bestSNR,bestLL);
	      fflush(stdout);
	    }
          } // for SNR = Smin .. Smax - 0.099
	  
	  // select the lower BND for minSNR estimate for current giter2
          double SNR = (giter2 == 0 && RefRepeats2 >= 3) ? bestSNR : (RefRepeats2 <= 2 || giter2 < RefRepeats2 - 2) ? max(bestSNR - Sstep, Smin) : bestSNR;

	  // restore startmaps,totalmaps, ScaleDelta 
          if(nummaps < orignummaps || bestSNR != origSNR || Hash_Bits != origHash_Bits || MHash_Bits != origMHash_Bits || HHash_Bits != origHHash_Bits || ScaleDelta != origScaleDelta){
	    nummaps = orignummaps;
	    startmaps = origstartmaps;
	    totalmaps = origtotalmaps;

	    ScaleDelta = origScaleDelta;NumScaleFactor = 1 + 2 * ScaleDelta;

	    if(VERB && ScaleDelta){
	      printf("Restored number of maps to %d and restored -ScaleDelta %0.3f %d\n",nummaps,ScaleDeltaSize,ScaleDelta);
	      fflush(stdout);
	    }

	    if(ScaleDelta > 3 && fabs(PixelLen - startPixelLen) >= 0.001){
	      if(VERB){
		printf("Restoring bpp = %0.2f -> %0.2f due to -ScaleDelta %0.3f %d\n",PixelLen * 1000.0, startPixelLen * 1000.0, ScaleDeltaSize,ScaleDelta);
		fflush(stdout);
	      }
	      PixelLen = startPixelLen;
            } // restore PixelLen : changes to maps will be made in mapcorrectApply() below

	    if(HashGen){
	      Hash_Bits = origHash_Bits;
	      MHash_Bits = origMHash_Bits;
	      HHash_Bits = origHHash_Bits;

	      // select the best minSNR value for hashtable
	      minSNR[0] = bestSNR; // (giter2 == 0 && RefRepeats2 >= 3) ? (bestSNR > origSNR ? min(bestSNR + Sstep, SNRmax) : max(bestSNR - Sstep , SNRmin)) : bestSNR;

	      if(VERB){
		printf("minSNR[0] -> %0.3f\n", minSNR[0]);
		fflush(stdout);
	      }
	      
	      mapcorrectApply(Gmap, 0, totalmaps, maxmaps, MapSNR, 1, numthreads);	  
	      if(ResBins[0] > 0){ /* (re)initialize rawsite[] and apply bias correction to site[] */
		BiasCorrect(Gmap,0,startmaps, 1);
		rawsitemaps = startmaps;
	      }

	      if(VERB>=2){
		printf("Before calling hash_generate():nummaps=%d\n",nummaps);
		fflush(stdout);
		dumpmemmap();      
	      }

	      if(DEBUG) assert(hash_prefix != NULL);

	      try {
		hash_generate(&Gmap[0],nummaps, &refmap[0], numrefmaps, 0, hash_prefix);
	      } catch (exception& e){
		cout << e.what() << endl;
		printf("hash_generate threw an exception\n");
		dumpmemmap();
		fflush(stdout);
		assert(0);
	      }
	      if(!hash_filename)
		goto Lfree;
	    }
          } // if(nummaps < orignummaps || bestSNR != origSNR || ScaleDelta != origScaleDelta ...)

	  if(VERB){
	    if(SNR != bestSNR)
	      printf("biasing minSNR[0] %0.3f -> %0.3f\n",bestSNR, SNR);
	    else
	      printf("minSNR[0] -> %0.3f (bestSNR= %0.3f)\n", SNR,bestSNR);
	    fflush(stdout);
	  }

	  minSNR[0] = SNR;

	  mapcorrectApply(Gmap, 0, totalmaps, maxmaps, MapSNR, 1, numthreads);	  
	  if(ResBins[0] > 0){ /* (re)initialize rawsite[] and apply bias correction to site[] */
	    BiasCorrect(Gmap,0,startmaps, 1);
	    rawsitemaps = startmaps;
	  }
        } // if(minSNRestimate ...

	if(VERB>=2){
	  printf("Before calling refalign_allpairs():\n");
	  fflush(stdout);
	  dumpmemmap();      
	}

	refalign_allpairs(nummaps, numrefmaps, NULL);
      } catch (exception& e){
	cout << e.what() << endl;
	printf("main():refalign_allpairs threw an exception\n");
	dumpmemmap();
	fflush(stdout);
	assert(0);
      }
    }

  } else { /* perform pairwise alignment */

    if(PairSplit)
      printf("Running pairwise multiple local alignment\n");
    else
      printf("No -ref file specified : running pairwise alignment\n");
    fflush(stdout);

    if(VERB/* HERE >=2 */){
      printf("PairSplit=%d, PairMerge=%0.1f, RefineOverwrite=%d\n", PairSplit,PairMerge,RefineOverwrite);
      fflush(stdout);
    }

    for(int c = 0; c < colors; c++)
      if(res[c] <= 0.0){
	printf("res[%d] parameter missing or zero (required to avoid optimistic Pvalue when intervals match exactly): using -res 3.5\n",c+1);
	res[c] = 3.5;
      }

    if(HashGen){
      if(DEBUG) assert(!spots_filename[0] && !(svcheck && numrefmaps > 0));
      try {
	if(Cfirst)
	  hash_generate(&Gmap[Cfirst],nummaps-Cfirst, &Gmap[0], Cfirst, 1, output_prefix);
	else
	  hash_generate(Gmap, nummaps, Gmap, nummaps, 1, output_prefix);
      } catch (exception& e){
	cout << e.what() << endl;
	printf("hash_generate threw an exception\n");
	dumpmemmap();
	fflush(stdout);
	assert(0);
      }
      if(!hash_filename)
	goto Lfree;
    }

    if(VERB>=2){
      printf("HashGen=%d,PairMerge= %0.3f, PairMergeRepeat=%d: completed hashtable call\n",HashGen,PairMerge,PairMergeRepeat);
      fflush(stdout);
    }

    if(PairMerge){

      if(SplitSegDupMaxOutlier <= 0.0)
	SplitSegDupMaxOutlier = PairMergeMaxOutlier;

      if(SplitSegDup > 0.0 && SplitSegDupE < PairMergeMaxEndKB){
	printf("WARNING: -SplitSegDup 3rd value %0.3f increased to %0.3f to match -pairmerge maxendoutlierKB value\n",
	       SplitSegDupE,PairMergeMaxEndKB);
	fflush(stdout);

	SplitSegDupE = PairMergeMaxEndKB;
      }
      if(SplitSegDup > 0.0 && SplitSegDupEN < PairMergeMaxEndN){
	printf("WARNING: -SplitSegDup 4th value %d increased to %d to match -pairmerge maxendoutlierN value\n",
	       SplitSegDupEN,PairMergeMaxEndN);
	fflush(stdout);
	
	SplitSegDupEN = PairMergeMaxEndN;
      }

      if(SplitSegDup > 0.0 && TrimmedNoMerge > 0.0){
	printf("WARNING: -SplitSegDup %0.3f used with -TrimmedNoMerge %0.3f may miss SegDups since contigs cannot merge opposite direction Extend & Split Contigs at a segdup\n",
	       SplitSegDup, TrimmedNoMerge);
	printf("\t Turning off -TrimmedNoMerge\n");
	fflush(stdout);
	
	TrimmedNoMerge = 0.0;
      }

      for(int i=0; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	pmap->paired = 1;/* set paired flag for all input maps (marks them as candidates for pairwise alignment) */
	pmap->centercloned = 0;/* reset previously merged flag for all input maps (marks if map was modified and hence orig* are no longer valid) */
	pmap->origid = Gmap[i]->id;
	if(VERB>=2 /* && pmap->id == 7585*/){
	  if(pmap->Mask[0])
	    printf("Gmap[%d]: id=%lld, N=%d, Mask[0][1]= %lx, Mask[0][N+1]= %lx\n",i, pmap->id, pmap->numsite[0], pmap->Mask[0][1], pmap->Mask[0][pmap->numsite[0]+1]);
	  else
	    printf("Gmap[%d]: id=%lld, N=%d, Mask[0] = %p\n",i, pmap->id, pmap->numsite[0], pmap->Mask[0]);
	  fflush(stdout);
	}
	if(!TrimmedNoMerge && SplitSegDup && pmap->Mask[0]){
	  /* change flagged ends from END_NOEXT to END_NOEXT2, so only newly created ends due to -SplitSegDup are flagged with END_NOEXT */
	  int N = pmap->numsite[0];
	  if(pmap->Mask[0][1] & END_NOEXT){
	    //	    if(DEBUG) assert(!(pmap->Mask[0][1] & END_NOEXT2));
	    pmap->Mask[0][1] |= END_NOEXT2;
	    pmap->Mask[0][1] &= ~END_NOEXT;
	    pmap->centercloned = 1;
	  }
	  if(pmap->Mask[0][N+1] & END_NOEXT){
	    //	    if(DEBUG) assert(!(pmap->Mask[0][N+1] & END_NOEXT2));
	    pmap->Mask[0][N+1] |= END_NOEXT2;
	    pmap->Mask[0][N+1] &= ~END_NOEXT;
	    pmap->centercloned = 1;
	  }
	}
	  
	if(Gmap[i]->contig){
	  delete Gmap[i]->contig;
	  Gmap[i]->contig = NULL;
	}
      }
    } // PairMerge

    char *orighash_filename = hash_filename;
    origSplitSegDup = SplitSegDup;
    int origPairMergeHmap = PairMergeHmap;
    double origPairMerge = PairMerge;
    double origPairMergeMaxOutlier = PairMergeMaxOutlier;
    double origPairMergeMaxEndExpandKB = PairMergeMaxEndExpandKB;
    //    double origPairMergeMaxEnd = PairMergeMaxEnd;
    //    double origPairMergeMaxEndN = PairMergeMaxEndN;

    if(SplitSegDup){
      printf("SplitSegDup = %0.3f -> 0\n",SplitSegDup);
      fflush(stdout);
      SplitSegDup = 0.0;
    }
    if(PairMergeHmap)
      PairMergeHmap = -1;

    int iter = 1, maxiter = 16 /* WAS88 99 */;
    for(; ; iter++){/* loop if PairMergeRepeat */
      if(DEBUG>=1+RELEASE) assert(!PairMergeIDrenumbered);

      if(DEBUG>=1+RELEASE && PairMerge) assert(!PairMergeCmapOutput);// make sure CMAPs have NOT been saved previously

      int progress = pairalign_allpairs(0, nummaps-1, iter, maxiter);
      if(VERB/* HERE HERE >=2 */){
	printf("iter=%d/%d:pairalign_allpairs(): progress=%d\n",iter,maxiter,progress);
	fflush(stdout);
      }
      if(!(PairMerge && PairMergeRepeat && progress))
	break;

      if(iter >= maxiter && !PairMergePreferNGS) { /* avoid infinite loop */
	printf("WARNING: exceeded %d iterations of PairMerge (terminating loop)\n",maxiter);
	fflush(stdout);
	break;
      }

      alignfree();/* reset alignment data structures */

      if(!(Cfirst && CfirstPairs==2))
	hash_filename = NULL;  /* disable use of hashtable for subsequent rounds */
      else if(hash_filename && HashGen){/* recompute hashtable */
#if 0 // HERE HERE Use Cfirst && CfirstPairs==2 (faster, but does not work yet)
	if(1){/* verify that Gmap[0..Cfirst-1]->id are all smaller than Gmap[Cfirst .. nummaps-1]->id : if not renumber ids of Gmap[Cfirst .. nummaps-1] */
	  long long id1 = Gmap[0]->id, id2 = Gmap[Cfirst]->id;
	  for(int i = 1; i < Cfirst; i++)
	    id1 = max(id1, Gmap[i]->id);
	  for(int i = Cfirst + 1; i < nummaps; i++)
	    id2 = min(id2, Gmap[i]->id);
	  
	  if(id1 >= id2){
	    printf("Gmap[0..%d] has largest id1= %lld, Gmap[%d..%d] has smallest id2= %lld : renumbering Gmap[%d..%d] so smallest id > id1\n",Cfirst-1,id1,Cfirst,nummaps-1,id2,Cfirst,nummaps-1);
	    fflush(stdout);
	    for(int i = Cfirst; i < nummaps; i++)
	      Gmap[i]->id = id1 + (long long)(i-Cfirst + 1);
	    //	    assert(id1 < id2);// NOTE : if this happens, fix by renumbering the  values of Gmap[Cfirst .. nummaps-1]->id to all be larger than id1
	  }
	}
	try {
	  hash_generate(&Gmap[Cfirst],nummaps-Cfirst, &Gmap[0], nummaps, 1, output_prefix);
	} catch (exception& e){
	  cout << e.what() << endl;
	  printf("hash_generate threw an exception\n");
	  dumpmemmap();
	  fflush(stdout);
	  assert(0);
	}
#else // use Cfirst == 0 to try to align all map pairs each time
	if(VERB){
	  printf("Reverting to Cfirst=0 and re-checking all map pairs: nummaps=%d\n",nummaps);
	  fflush(stdout);
	}
	Cfirst = CfirstPairs = 0;
	try {
	  hash_generate(Gmap,nummaps,Gmap,nummaps, 1, output_prefix);
	} catch (exception& e){
	  cout << e.what() << endl;
	  printf("hash_generate threw an exception\n");
	  dumpmemmap();
	  fflush(stdout);
	  assert(0);
	}
#endif
      }

    }
    
    hash_filename = orighash_filename;
    if(DEBUG) assert(!SplitSegDup);
    if(PairMerge && PairMergeRepeat && origSplitSegDup){
      if(VERB){
	printf("SplitSegDup = %0.3f -> %0.3f\n",SplitSegDup, origSplitSegDup);
	printf("TrimmedNoMerge = %0.1f -> %0.1f kb\n", TrimmedNoMerge, PAIRMERGE_HUGE);
	if(Truncate)
	  printf("Truncate = %0.1f -> 0.0\n",Truncate);
	if(PairMerge != origSplitSegDup)
	  printf("Min Overlap = %0.3f -> %0.3f kb\n", PairMerge, origSplitSegDup);
	if(PairMergeMaxOutlier != SplitSegDupMaxOutlier)
	  printf("PairMergeMaxOutlier = %0.3f -> %0.3f kb\n", PairMergeMaxOutlier, SplitSegDupMaxOutlier);
	printf("PairMergeMaxEndExpandKB = %0.1f -> %0.1f kb\n", PairMergeMaxEndExpandKB,PAIRMERGE_HUGE);
	fflush(stdout);
      }

      SplitSegDup = origSplitSegDup;

      TrimmedNoMerge = PAIRMERGE_HUGE;
      Truncate = 0.0;
      
      PairMerge = SplitSegDup;
      PairMergeMaxOutlier = SplitSegDupMaxOutlier;
      PairMergeMaxEndExpandKB = PAIRMERGE_HUGE;
      /*      PairMergeMaxEnd = PAIRMERGE_HUGE;
	      PairMergeMaxEndN = 9999999;*/

      if(DEBUG) assert(PairMerge && PairMergeRepeat && SplitSegDup);
    }

    if(PairMerge && PairMergeRepeat && SplitSegDup){/* need to call pairalign() with SplitSegDup starting with Cfirst = 0 to recheck all pairwise alignments */
      if(VERB/* HERE >=2 */){
        printf("Calling pairalign() again to check for -SplitSegDup\n");
        fflush(stdout);
      }

      alignfree();/* reset alignment data structures */

      if(HashGen && hash_filename){
	try {
	  hash_generate(Gmap, nummaps, Gmap, nummaps, 1, output_prefix);
	} catch (exception& e){
	  cout << e.what() << endl;
	  printf("hash_generate threw an exception\n");
	  dumpmemmap();
	  fflush(stdout);
	  assert(0);
	}
      }
      
      if(VERB>=2 && PairMergeIDrenumber){
	printf("After SplitSegDup hashtable sorted Gmap[0..%d]:\n",nummaps-1);
	for(int i = 0; i < nummaps; i++){
	  Cmap *pmap = Gmap[i];
	  int M = pmap->numsite[0];
	  printf("Gmap[%d]=%p: mapid=%d,id=%lld,paired=%d,len=%0.4f,M=%d\n",
		 i,pmap,pmap->mapid,pmap->id,pmap->paired,pmap->site[0][M+1],M);
	}
	fflush(stdout);
      }

      Cfirst = 0;

      /* change flagged ends from END_NOEXT to END_NOEXT2, so only newly created ends due to -SplitSegDup are flagged with END_NOEXT */
      for(int i = 0; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	if(pmap->Mask[0]){/* remove flagged ends, so only newly created ends due to -SplitSegDup are subject to TrimmedNoMerge (particularly during later -pairmergeHmap stage */
	  int N = pmap->numsite[0];
	  if(pmap->Mask[0][1] & END_NOEXT){
	    //	    if(DEBUG) assert(!(pmap->Mask[0][1] & END_NOEXT2));
	    pmap->Mask[0][1] |= END_NOEXT2;
	    pmap->Mask[0][1] &= ~END_NOEXT;
	    pmap->centercloned = 1;
	  }
	  if(pmap->Mask[0][N+1] & END_NOEXT){
	    //	    if(DEBUG) assert(!(pmap->Mask[0][N+1] & END_NOEXT2));
	    pmap->Mask[0][N+1] |= END_NOEXT2;
	    pmap->Mask[0][N+1] &= ~END_NOEXT;
	    pmap->centercloned = 1;
	  }
	}
      }

      // Cfirst==0 will re-enable -MaxPalindrome and -MaxInvPalindrome options (self alignment) for final breakup, but TrimmedNoMerge will prevent any re-merges

      int maxiter = 99;

      for(;++iter; ){
	if(DEBUG>=1+RELEASE) assert(!PairMergeIDrenumbered);

	if(DEBUG>=1+RELEASE && PairMerge) assert(!PairMergeCmapOutput);// make sure CMAPs have NOT been saved previously

	int progress = pairalign_allpairs(0, nummaps-1, iter, maxiter);
	if(!progress)
	  break;

	if(iter >= maxiter){	// Limit iterations to avoid infinite loop due to incorrectly handled corner cases
	  printf("WARNING: exceeded %d iterations of PairMerge with SplitSegDup : possible logic bug (terminating loop)\n",maxiter);
	  fflush(stdout);
	  break;
	}

	alignfree();/* reset alignment data structures */

	// subsequent iterations allow duplicate Cloned regions to be merged together (if fully overlapped, TrimmedNoMerge will not have any effect)
	// subsequent iterations are also required to catch all -SplitSegDup cases for large contigs that may have multiple SegDup regions


	/* disable use of hashtable for subsequent rounds */
	hash_filename = NULL;
      }
    }


    hash_filename = orighash_filename;

    if(VERB){
      if(PairMergeHmap != origPairMergeHmap)
	printf("PairMergeHmap = %d -> %d\n",PairMergeHmap, origPairMergeHmap);
      if(PairMerge != origPairMerge)
	printf("PairMerge -> %0.3f\n", origPairMerge);
      if(PairMergeMaxOutlier != origPairMergeMaxOutlier)
	printf("PairMergeMaxOutlier -> %0.3f\n", origPairMergeMaxOutlier);
      if(PairMergeMaxEndExpandKB != origPairMergeMaxEndExpandKB)
	printf("PairMergeMaxEndExpand -> %0.3f\n", origPairMergeMaxEndExpandKB);
      fflush(stdout);
    }

    PairMergeHmap = origPairMergeHmap;

    PairMerge = origPairMerge;
    PairMergeMaxOutlier = origPairMergeMaxOutlier;
    PairMergeMaxEndExpandKB = origPairMergeMaxEndExpandKB;
    //    PairMergeMaxEnd = origPairMergeMaxEnd;
    //    PairMergeMaxEndN = origPairMergeMaxEndN;

    if(PairMerge && PairMergeHmap){/* need to call pairalign() one more time without Cfirst */
      alignfree();/* reset alignment data structures */

      Cfirst = 0;
      
      if(HashGen && hash_filename){
	try {
	  hash_generate(Gmap, nummaps, Gmap, nummaps, 1, output_prefix);
	} catch (exception& e){
	  cout << e.what() << endl;
	  printf("hash_generate threw an exception\n");
	  dumpmemmap();
	  fflush(stdout);
	  assert(0);
	}
      }
      
      if(VERB>=2 && PairMergeIDrenumber){
	printf("After PairMergeHmap hashtable sorted Gmap[0..%d]:\n",nummaps-1);
	for(int i = 0; i < nummaps; i++){
	  Cmap *pmap = Gmap[i];
	  int M = pmap->numsite[0];
	  printf("Gmap[%d]=%p: mapid=%d,id=%lld,paired=%d,len=%0.4f,M=%d\n",
		 i,pmap,pmap->mapid,pmap->id,pmap->paired,pmap->site[0][M+1],M);
	}
	fflush(stdout);
      }

      if(DEBUG>=1+RELEASE) assert(!PairMergeIDrenumbered);

      iter++;

      if(DEBUG>=1+RELEASE && PairMerge) assert(!PairMergeCmapOutput);// make sure CMAPs have NOT been saved previously

      (void)pairalign_allpairs(0, nummaps-1, iter, 999999);
    } // PairMerge && PairMergeHmap

    if(DEBUG && PairMerge) assert(PairMergeCmapOutput);// make sure CMAPs were saved

#if 0 // files have already been output in pairalign(), so no point in changing maps now
    if(PairMerge && SplitSegDup){/* restore END_NOEXT flags saved as END_NOEXT2 */
      for(int i = 0; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	int N = pmap->numsite[0];
	if(pmap->Mask[0]){
	  if(pmap->Mask[0][1] & END_NOEXT2){
	    pmap->Mask[0][1] |= END_NOEXT;
	    pmap->Mask[0][1] &= ~END_NOEXT2;
	  }
	  if(pmap->Mask[0][N+1] & END_NOEXT2){
	    pmap->Mask[0][N+1] |= END_NOEXT;
	    pmap->Mask[0][N+1] &= ~END_NOEXT2;
	  }
	}
      }
    }
#endif


  } // pairwise alignment */

  if(XmapStatWrite && minSNRestimate){/* output XmapCount,XmapSites and XmapLength to specified file (after estimating minSNR) */
    FILE *fp;
    if(VERB){
      printf("Generating -i input map statistics in %s (startmaps=%d, minSNR[0]= %0.6f)\n",XmapStatWrite, startmaps, minSNR[0]);
      fflush(stdout);
    }
    if((fp = fopen(XmapStatWrite,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for writing:errno=%d:%s\n",XmapStatWrite,eno,err);
      fflush(stdout);exit(1);
    }

    /* write out commandline */
    printversion(fp);

    /* compute statistics */
    for(int c = 0; c < colors; c++){
      XmapCount[c] = XmapSites[c] = 0;
      XmapLength[c] = XmapGtheta[c] = 0.0;
      for(int i = 0; i < startmaps; i++){
	Cmap *pmap = Gmap[i];
	int M = pmap->numsite[c];
	FLOAT *X = pmap->site[c];
	XmapCount[c]++;
	XmapSites[c] += M;
	XmapLength[c] += X[M+1];
	for(int k= 1; k < M; k++)
	  XmapGtheta[c] += log(X[k+1]-X[k]);
      }
      XmapGtheta[c] /= max(1,XmapSites[c] - XmapCount[c]);
      XmapGtheta[c] = exp(XmapGtheta[c]);
    }
    fprintf(fp,"XmapCount=");
    for(int c = 0; c < colors; c++)
      fprintf(fp," %d", XmapCount[c]);
    fprintf(fp,"\nXmapSites=");
    for(int c = 0; c < colors; c++)
      fprintf(fp," %lld", XmapSites[c]);
    fprintf(fp,"\nXmapLength=");
    for(int c = 0; c < colors; c++)
      fprintf(fp," %0.3f", XmapLength[c] /* NEW12 */ * origPixelLen / PixelLen);
    fprintf(fp,"\nXmapGtheta=");
    for(int c = 0; c < colors; c++)
      fprintf(fp," %0.6f", XmapGtheta[c] /* NEW12 */ * origPixelLen / PixelLen);
    fprintf(fp,"\n");
    fclose(fp);
  }

  if(HashGen && !HashTxt && HASH_STREAM < 3){
    char filename[PATH_MAX];
    sprintf(filename,"%s_id%lld.hashbin",output_prefix, CMapID);
    errno = 0;
    if(unlink(filename)){
      int eno = errno;
      char *err = strerror(eno);
      printf("WARNING:unlink(%s) failed:errno=%d:%s\n",filename,eno,err);
      fflush(stdout);
    }
  }

  if(simpleRepeat)
    findRepeatSimpleMaps(simpleRepeatFilter, simpleRepeatTolerance, simpleRepeatMinEle, false); //last arg is verbose

 Lfree:  /* free up memory */

  if(VERB>=2){
    dumpmemmap();      
    fflush(stdout);
  }

  fflush(stderr);

  if(VERB){
    static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
    getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
    printf("VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb : CPU time= %0.6f, wall time= %0.6f\n",VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(), wtime());
    fflush(stdout);
  }

  if(LEAKDEBUG){
    if(usecolor)
      colors = 2;

    delete [] ScaleFactor;
    if(VERB/* HERE HERE >=2 */){
      size_t sites = 0;
      for(size_t i = 0; i < numaligns; i++)
	sites += alignment[i]->numpairs;
      printf("nummaps=%d,numrefmaps=%d,numaligns=%lu(sites=%lu),sizeof(Calign)=%lu:map=%p,refmap=%p,alignment=%p\n",nummaps,numrefmaps,numaligns,sites,sizeof(Calign),Gmap,refmap,alignment);
      fflush(stdout);
    }
    if(VERB){
      printf("Calling malloc_stats()\n");
      fflush(stdout);
      malloc_stats();
      fflush(stderr);
    }

    malloc_trim(0);
    if(VERB){
      static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
      getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
      printf("After malloc_trim(0):VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb : CPU time= %0.6f, wall time= %0.6f\n",VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(), wtime());
      fflush(stdout);
    }
    if(VERB){
      printf("Calling malloc_stats()\n");
      fflush(stdout);
      malloc_stats();
      fflush(stderr);
    }

    alignfree();
    if(VERB){
      static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
      getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
      printf("After alignfree():VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb : CPU time= %0.6f, wall time= %0.6f\n",VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(), wtime());
      fflush(stdout);
    }

    if(Gmap){
      if(DEBUG) assert(totalmaps <= maxmaps);
      for(int i = 0; i < totalmaps; i++){
	if(Gmap[i]){
	  Gmap[i]->allfree();
	  Gmap[i] = NULL;
	}
      }
      delete [] Gmap;
      Gmap = NULL;
    }
    if(refmap){
      if(DEBUG) assert(numrefmaps <= maxrefmaps);
      for(int i = 0; i < maxrefmaps; i++){
	if(refmap[i]){
	  refmap[i]->allfree();
	  refmap[i] = NULL;
	}
      }
      delete [] refmap;
      refmap = NULL;
    }

    extern char *Nickase[MAXCOLOR];/* see input_vfixx.cpp */
    for(int i = 0; i < colors;i++)
      if(Nickase[i]){
	/*      printf("Nickase[%d]=x%p\n",i,Nickase[i]);
		fflush(stdout);*/
	free(Nickase[i]);
	Nickase[i] = 0;
      }

    maxmap_free();
    if(stderr_file) free(stderr_file);
    if(stdout_file) free(stdout_file);

    if(RunDataListMax > 0){
      for(int i = 0; i < RunDataListLen; i++)
	free(RunDataList[i]);
      delete [] RunDataList;
    }
    for(int i = 0; spots_filename[i]!=NULL; i++)
      (void)free(spots_filename[i]);
    for(int i = 0; i < num_files; i++)
      free(vfixx_filename[i]);

    if(xmapentries){
      for(int i = 0; i < numxmapentries; i++)
	delete xmapentries[i];
      delete [] xmapentries;
      xmapentries = NULL;
    }
    if(usexmapentries){
      delete [] usexmapentries; 
      usexmapentries = NULL;
    }

    if(parametersfile) {
      free(parametersfile);
      parametersfile = NULL;
    }

    malloc_trim(0);
    if(VERB){
      static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
      getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
      printf("After malloc_trim(0):VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb : CPU time= %0.6f, wall time= %0.6f\n",VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(), wtime());
      fflush(stdout);
    }

    if(VERB){
      printf("Calling dumpmemmap()\n");
      dumpmemmap();    
      fflush(stdout);
    }
    if(VERB){
      printf("Calling malloc_stats()\n");
      fflush(stdout);
      malloc_stats();
      fflush(stderr);
    }

    if(VERB){
      static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
      getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
      printf("VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb : CPU time= %0.6f, wall time= %0.6f\n",VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(), wtime());
      fflush(stdout);
    }
  }

  STOP(main);
  printTimers(stdout);
  fflush(stdout);

  time_t my_time = time(NULL); 
  printf("END TIME: %s", ctime(&my_time));   // ctime() used to give the present time 

  printf("END of output\n");
  fflush(stdout);

#ifndef WIN32
  if(stdout_file){
    int fd = fileno(stdout);
    if(fd < 0){
      fprintf(stdout,"fileno(stdout) returned -1\n");
      fflush(stdout);
    } else {
      if(fsync(fd) < 0){
	char *err = strerror(errno);
	printf("fsync(fd=%d) returned -1:%s\n",fd,err);
	printf("STDOUT_FILENO=%d\n",STDOUT_FILENO);
	printf("CPU time= %0.6f, wall time= %0.6f\n",mtime(),wtime());
	printf("END of output\n");
	fflush(stdout);
      }
    }
  }
#endif

  if(stderr_file && !(stdout_file && !strcmp(stdout_file,stderr_file))){
    fprintf(stderr,"CPU time= %0.6f, wall time= %0.6f\n",mtime(),wtime());
    fprintf(stderr,"END of output\n");
    fflush(stderr);

#ifndef WIN32
    int fd = fileno(stderr);
    if(fd < 0){
      fprintf(stderr,"fileno(stderr) returned -1\n");
      fflush(stderr);
    } else {
      if(fsync(fd) < 0){
	char *err = strerror(errno);
	fprintf(stderr,"fsync(fd=%d) returned -1:%s\n",fd,err);
	fprintf(stderr,"STDERR_FILENO=%d\n",STDERR_FILENO);
	fprintf(stderr,"CPU time= %0.6f, wall time= %0.6f\n",mtime(),wtime());
	fprintf(stderr,"END of output\n");
	fflush(stderr);
      }
    }
#endif
  }

  exit(0);
  return 0;
}
