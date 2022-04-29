#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <math.h>
#include <ctype.h>
#include <fcntl.h>
#include <malloc.h>

#ifndef WIN32
#include <unistd.h>
#endif

#include "globals.h"

//#undef DEBUG
//#define DEBUG 2

#include "Assembler_parameters.h"
#include "Assembler_graph.h"
#include "Ccontig.h"

#if CALIGN_SMALL
//#define Calign CalignA
//#define Calign_block CalignA_block

//#define alignment alignmentA
#define align_block alignA_block
#define num_align_blocks num_alignA_blocks

//#define align_blockalloc alignA_blockalloc
#endif

/** \file Assembler.cpp \brief Main Routine and I/O for solving the de novo graph. 
*
* Significant solving done in calls to Assembler_graph.h and 
* Assembler_graph.cpp
*
* Minimal class structure
*/

extern void hashROC(char *);
extern void input_refmap(char *filename);
extern double mtime();
extern double wtime();

static char buf[LINESIZ];

Cmap **YYmap = 0,**XXmap= 0;/* required by IsRepeatRegion : defined below after reading maps */
int numY= 0,numX= 0;

extern Cid2mapid *id2mapid;

const char *SVN_ID = (char *)"$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Assembler.cpp 11529 2020-08-22 21:05:37Z tanantharaman $";

char *usage = (char *)"Assembler \n"
  "-a input.align : input alignments (Use more than one -a option to merge multiple alignment files)\n"
  "-af align.filenames : input alignments from all files listed in align.filenames (one per line), as if they had been input with multiple -a options\n"
  "-i [input.vfixx | input.cmap] : input map from .vfixx or .cmap file (Use more than one -i option to merge multiple input files)\n"
  "-if input.filenames : input maps from all files listed in input.filenames (one per line), as if they had been input with multiple -i options\n"
  "  -bnxerror : If input map from .bnx file has sites beyond the end of the map, exit with an error message (Default is print WARNING and extend map ends)\n"
  "  -SNR : Require input .bnx file to have SNR information (one line per color, starting with QX02) [Default : OFF]\n"
  "  -stitch : Require input .bnx file to have image stitch information (one line per color, starting with QX04) [Default : OFF]\n"
  "  -bnxmergeerror <N> : If <N> = 1 : Don't allow BNX 1.0 errors like missing Run Data or duplicate Run Data lines (<N> = 0 : just print warning) [Default 0]\n"
  "  -bnxBigid <N> : If <N> = 1 : Use 19 digit id encoding, if possible (If not possible due to field width overflows, use sequential ids starting at 1)\n"
  "                  If <N> = 0 : Avoid changing molecule ids (If duplicate ids in input(s), use sequential ids starting at 1) [Default 0]\n"
  "                  If any ids are changed, output <prefix>_renumbered.bnx showing the new ids\n"
  "-hashROC [input.hashbin | input.hash] [-hashtxt] : input Hashtable (.hashbin in binary format) and compute ROC curve (requires map and alignment input data) : No Assembly performed\n"
  "                                    If -hashtxt is specified, a text format version of the HashTable is generated from the binary .hashbin\n"
  "  -dumpalign : just dump the sorted alignments (use specified hashfile name as output prefix)\n"
  "  -first <N> : option used during pairwise alignment (needed by -hashROC to determine the total number of pairwise alignments involved)\n"
  "-ref input.cmap <provide true reference map for debugging>\n"
  "-reff ref.filenames : input reference maps from all files listed in ref.filesnames (one per lien), as if they had been input as multiple -ref options\n"
  "-refmap input.map <provide true reference alignments for debugging>\n"
  "-contigs input.contigs start end <provide draft contigs assembled previously and just run refinement on specified range of contigs>\n"
  "-contigs_format <0/1> : 0 means use original text format, 1 means use binary format for faster IO of .contigs file [Default : 0]\n"
  "-o [output prefix] : use specified prefix name for all output file names [Default : It is an error not to specify this option]\n"
  "-output-veto-filter regex : regular expression to specify which files to not write. [Default None]\n"
  "-output-filter regex : regular expression to specify which files to write. [Default .*]\n"
  "-stdout [file.stdout] : Redirect stdout to file.stdout. Does NOT apply to error messages from this and any previous options on the command line. <file.stdout> defaults to <prefix>.stdout, if -o was specified previously\n"
  "-stderr [file.stderr] : Redirect stderr to file.stderr. Does NOT apply to error messages from this any any previous options on the command line. <file.stderr> defaults to <prefix>.stdout, if -o was specified previously\n"
  "-f | -force : force output file overwrites [Default: error exit if output file already exists]\n"
  "-colors <N> : <N> is the number of colors required in the input data : currently any value other than 1 or 2 is not supported. [Default : 1 or value in first -i input file header]\n"
  "-usecolor <N> : If -i input data has 2 colors only the <N> th color is used as if only that one color were present. Any -ref input must then have 1 color\n"
  "-FP <False Positives per 100kb, one per color>[Default 1.0]\n"
  "-FN <False Negative Probability, one per color>[Default 0.15]\n"
  "-sd <scaling error in root-kb, one per color>[Default 0.15]\n"
  "-sf <fixed error in kb, one per color>[Default 0.05]\n"
  "-sr <relative sizing error, one per color>[Default 0.00]\n"
  "-sm <fixed error in kb for graph distances[Default 2.0]\n"
  "-bpp <Bases Per Pixel>[Default : from vfixx file]\n"
  "-res <resolution in pixels, one value per color>[Default 3.5] : Depends on bpp value\n"
  "-resSD <standard deviation of resolution in pixels, one value per color>[Default 0.75] : Depends on bpp value\n"
  "-mres <pixels> : Reduces resolution of -i input maps to specified multiple of 0.5 kb [Default: 0.001]\n"
  "-mresSD <pixelsSD> : Varies the resolution of -i input maps by randomly varying -mres value by specified Gaussian SD. [Default 0]: introduces randomness in map resolution\n"
  "                    NOTE : -mres and -mresSD are applied after -bpp\n"
  "-rres <pixels> : Limit resolution of refined consensus maps to this value [Default: 0]\n"
  "-readparameters <file.errbin> : read in error parameters from specified file generated by a previous -M run.\n"
  "            NOTE: error values in <file.errbin> silently override any command line error values before are after this option\n"
  "-MinFP <V> : Force FP values to be at least <V> at all times [Default 0.01]\n"
  "-MaxFP <V> : Force FP values to be no more than <V> at all times [Default 4.0]\n"
  "-MinFN <V> : Force FN values to be at least <V> at all times [Default 0.001]\n"
  "-MinSF <V> : Force sf values to be at least <V> at all times [Default 0.001]\n"
  "-MaxSF <V> : Force sf values to be no more than <V> at all times [Default 1.00]\n"
  "-MinSD <V> : Force sd values to be at least <V> at all times [Default -0.50]\n"
  "-MinSD_Mult <F> : Limit -ve values of MinSD so sizing variance cannot be reduced by more than fraction F vs MinSD value of 0 [Default 0.99]\n"
  "-MaxSD <V> : Force sd values to be no more than <V> at all times [Default 0.50]\n"
  "-MinSR <V> : Force sr values to be at least <V> at all times [Default 0.00]\n"
  "-MaxSR <V> : Force sr values to be no more than <V> at all times [Default 0.50]\n"
  "-MaxSE <V> : Force se values to be no more than <V> at all times [Default 0.00]\n"
  "-MinRes <V> : Force res value to be at least <V> at all times [Default 0.01]\n"
  "-MaxRes <V> : Force res value to be no more than <V> at all times [Default 9.99]\n"
  "-MinResSD <V> : Force res valueSD to be at least <V> at all times [Default 0.01]\n"
  "-MaxResSD <V> : Force res value to be no more than <V> at all times [Default 9.99]\n"
  "-Kmax <K> : Specify <K> as the maximum number of consecutive site intervals in reference that might (mis)resolve into a single site in query [Default 6]\n"
  "-deltaY <N> : N = largest reference map alignment interval [Default 6]\n"
  "-deltaX <N> : N = largest query map alignment interval [Default 4]\n"
  "-extend <N> : If N > 0 allow query map to extend beyond reference map (if N = 2, also allow refinement of extended reference). Only used with -refine  [Default 0]\n"
  "-draftsize <N> : Type of draft sizing computation (0 = old method, 1 = new method using -outlier to try to filter out outliers)\n"
  "-outlier <Pvalue> [<StitchPV>]: outlier(blemish) Pvalues per true positive site. If -stitch, then <StitchPV> specifies a seperate Pvalue for intervals with a image stitch [Default 0.0003 0.01]\n"
  "    -outlierMax <delta> : Disallow outliers with sizing difference larger than <delta> in kb [Default OFF]. NOTE : only works with -maptype 0 (not recommended with stitching errors)\n"
  "    -outlierLambda <lambda> : Decrease -outlier Pvalue by exp(-|x-y|/lambda), where x-y is sizing difference of the aligned interval in kb. OFF if <lambda> is greater than 1000.0 kb. [Default OFF]\n"
  "-endoutlier <end outlier(blemish or chimeric junction) prior probability per true positive site>[Default 0]\n"
  "-maptype <n> : Specify if -i input (Query) map is Molecule (<n> = 0) or Assembly/Contig (<n> = 1). By default .vfixx or .bnx is assumed to be Molecules and .cmap Assembly/Contig\n"
  "                   Molecules are only permitted to have outliers with Query interval smaller than Reference interval (representing DNA knots)\n"
  "                   Assembly/Contig may have outliers with Query interval smaller or larger than Reference interval (representing SV insertions or deletions)\n"
  "-alignscore : input detailed alignment scores from .align[Default OFF]\n"
  "-MinCov <min-coverage> : Find longest assembly paths with this minimum coverage [Default 2]\n"
  "-MinRatio <min-coverage-ratio> : Delete mismatched path if coverage is this much lower than other path [Default 0.20]\n"
  "-CycleRatio <min-coverage-ratio> : Delete weakest cycle edge if coverage is this much lower than average cycle edge [Default 0.20]\n"
  "-MinEdgeRatio <min-Edge-ratio> : Delete Edge if coverage is this much lower than another edge from same node in same direction [Default 0.20]\n"  
  "-ChimRatio <Chimeric-coverage-ratio> : Delete chimeric nodes with node coverage to forward/backward edge coverage below this value [Default 0.20]\n"
  "-S <score-threshold>[Default 0.0]: only use alignments with score above this threshold\n"
  "-T <Pvalue-threshold>[Default 0.0001]: only use alignments with Pvalue below this threshold\n"
  "-A <Aligned-Sites-threshold>[Default 4] : only use alignments with at least this many sites\n"
  "-L <Aligned-Length-threshold>[Default 0] : only use alignments with at least this length (in kb) aligned\n"
  "-minlen <L> : If L != 0 : filter out all -i input maps that are smaller than L kb [Default 0]\n"
  "-maxlen <L> : If L != 0 : filter out all -i input maps that are larger than L kb [Default 0]\n"
  "-minsites <N> : filter out all -i input maps with fewer than N sites [Default 1]\n"
  "-maxsites <N> : If N != 0 filter out all -i input maps with more than N sites [Default 0]\n"
  "                NOTE : If -colors 2 and -usecolor is NOT applied, -minsites and -maxsites are based on the total number of sites on both colors\n"
  "-maxSiteDensity <D> : If D != 0 filter out all -i input maps with more than D sites per 100 kb [Default: OFF or 0]\n"
  "-maxEnd <L> : Trim unlabeled molecule ends to no more than <L> kb [Default OFF]\n"
  "-minEnd <L> : Extend unlabeled molecule ends to at least <L> kb. Also applies to map fragments generated by -break or -contigsplit [Default 0.020 kb]\n"
  "-maxInterval <L> : If label interval is larger than <L> break molecule into 2 pieces [Default OFF]\n"
  "-MaxIntensity <V> : If V != 0 : Filter out all -i input maps with AvgIntensity above V [Default 0]\n"
  "-minSNR <V> ... : Remove all -i input map sites with SNR below V (There is One V value per color) [Default 0.0]\n"
  "-maxmem <Maximum Memory in Gbytes> [Default 256] : avoid using more than specified amount of memory by reducing number of threads (during Refinement). Only used during Assembly with -AlignmentFilter\n"
  "    -maxmemIncrease <ResMem> : If SysMem > 0 : Increase maxmem to system memory size - ResMem (in GB) [Default OFF or 0]\n"
  "-maxvirtmem <Gbytes> : avoid using more than specified amount of virtual memory by reducing number of threads etc. 0 means unlimited [ Default : 3x -maxmem]\n"
  "-refine <level> : refine each assembled contig consensus map (<level> = 1 : sizes only(slow), <level> = 2 : sizes and sites(very slow))[Default 0]\n"
  "  -RefineOverwrite <0/1> : If 0, skip contigs which already exist in refined form; otherwise overwrite previously refined contig [Default 1]\n"
  "  -MultiMode : Enable Multi Mode sizing optimization (slower) [Default OFF]\n"
  "  -RefineStep <R1> <R2> [ <R3> <R4> <N> <R5> <R6> <R7>]: Use R1 as initial sizing change step size during 1st refinement stage, and R2 during subsequent stages (or all stages of RefineA). [Default 0.10 0.20]\n"
  "                           Use R3 & R4 instead of R1 and R2 (respectively) after initial sizing change and limit number of steps (size values) to N. This is only done of previous round used more\n"
  "                           than 2 values AND best value was at range end OR either of neighboring interval changed by more than <R7> done if. [Default R3=0.10 R4=0.20 N=3 R7=0.10]\n"
  "                           Reducing R1 & R2 can improve accuracy of consensus map and ability to handle multimodal maxima at expense of larger run time\n"
  "                           Note that the number of steps (size values) with R3 & R4 is limited to N, so just reducing R3 & R4 trades off search space for precision, without changing run time\n"
  "                           The minimum step sizes with R1 or R3 are R5 (kb) and with R2 or R4 are R6 (kb) [Default 0.2 kb]\n"
  "  -RefineSizeRange <Low1> <High1> <Low2> <High2> : Use ratios <Low1> and <High1> during initizing sizing change as scan range relative to current interval size\n"
  "                            Use ratios <Low2> and <High2> during subsequent stages (or all stages of RefineA). [Default: 0.01 2.0 0.7 1.4]\n"
  "  -SiteMerge <0/1> : During refinement periodically try to merge all nearby labels, if that improves likelihood [Default : 1]\n"
  "  -RefineErr <E1> <E2> <E3> : Allow reduced likelihood accuracy specified by log-likelihood errors of E1, E2 and E3 when computing endoutlier likelihood, total alignment likelihood and\n"
  "                              sub-alignment likelihood respectively. Increasing E1,E2,E2 will speed up refinement with more predictable loss in accuracy by reducing deltaX, deltaY\n"
  "                              and RefineRange dynamically, compared to using fixed reduced values [Default : 1e-8 0 0]\n"
  "                              NOTE: Any non-zero value of E3 will incur a small run-time overhead\n"
  "  -LRbias <bias> : A minimum LR (Likelihood Ratio) value per map. This value is added to the actual LR to reduce the influence of maps with low LR values. [Default 0]\n" 
  "  -MinSplitLen <len> : minimum length in kb of refined map [Default 100 kb]\n"
  "                         If contig is below the minimum length an empty file PREFIX_contig<N>_refined.cmap will be generated\n"
  "  -deltaXRef <N> : change deltaX during refinement to N. Useful since -outlierExtend does NOT apply during refinement [Default : deltaX]\n"
  "  -deltaYRef <N> : change deltaY during refinement to N. Useful since -outlierExtend does NOT apply during refinement [Default : deltaY]\n"
  "    -deltaExtXRef <M> : If M > deltaX, allow query map intervals to be as large as M, whenever reference map intervals are <= deltaY [Default : deltaX]\n"
  "    -deltaExtYRef <M> : If M > deltaY, allow reference map intervals to be as large as M, whenever query map intervals are <= deltaX [Default : deltaY]\n"
  "                      NOTE : Using larger values for deltaExtXRef, deltaExtYRef is faster than using the same values for deltaXRef,deltaYRef but \n"
  "                             can handle the same size Haplotype insertion or deletion (but not a Haplotype Inversion)\n"
  "    -deltaEndRef <U> : If U > 0 : allow non-outlier ends with up to deltaExtXRef (or deltaExtYRef) + U - 1) Labels, instead of deltaXRef,deltaYRef, during refinement [Default : 0]\n"
  "  -RefineRange <R> <R2> : R = largest change in map alignment (in labels) allowed during each Refinement step [Default 8]\n"
  "                        Reducing this value to 4 will dramatically speed up refinement but may sometimes produce a worse consensus map\n"
  "                        R2 is a larger value used only for initial alignment and periodic re-alignment during refinement (Should at least equal deltaXRef) [Default 8]\n"
  "  -RefineXRange <0/1/2> : If 1, use extrapolation of alignments based on both RefineRange and deltaExtXRef (slightly more stable and slower refinement) [Default: OFF or 0]\n"
  "  -LocalRange <0/1> : If 1, Use label interval size of current molecules instead of global average value : will speed up refinement for high label density regions. [Default 0]\n"
  "  -LocalTheta <0/1> : If 1, Use local estimate of consensus label interval size instead of genome wide average value : will speed up refinement for high label density regions. [Default 0]\n"
  "  -RefineAlignExt <0/1> : If 1, force extrapolation of alignments to reach ends of molecules (slightly more stable and slower refinement) [Default: OFF or 0]\n"
  "  -RefineFillMap <0/1> : If 1, use older method for extrapolation of alignments during non-haplotype refinement [Default: ON or 1]\n"
  "  -endoutlierRef <E1> <E2> ... <En> : -endoutlier Pvalue sequence used during successive phases of refinement (if -extend 2) to converge to the consensus estimate of the extension region\n"
  "                        This can help when the extension consists of a mixed population by allowing the lower scoring extensions to become outliers and stop contributing to the consensus\n"
  "  -outlierRef <P1> <P2> ... <Pn> : -outlier Pvalue sequence used during successive phases of refinement (if -extend 2) to converge to the consensus estimate of the extension region\n"
  "                        Defaults to the same sequence as for -endoutlierRef. Length of sequence will be matched to -endoutlierRef by deleting the first few values OR adding values\n"
  "                         at the begining matching those of -endoutlierRef \n"
  "  -outlierTypeRef <0/1> : Modifies outlier penalty during Refinement [Default 1]\n"
  "                          If 0, don't penalize misaligned sites inside an outlier\n"
  "                          If 1, Penalize misaligned sites, even inside an outlier\n"
  "  -outlierLambdaSwitch [<Label> <Size>] : suppress outlierLambda (or force it to <Label>) when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "                                          restore outlierLambda (or force it to <Size>) when adjusting sizing (or checking for Haplotype indels) during refinement\n"
  "  -outlierTypeSwitch [<Label> <Size>] : restore outlierTypeRef to 1 (or force it to <Label>) when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "			                   change outlierTypeRef to 0 (or force it to <Size>) when adjusting sizing (or checking for Haplotype indels) during refinement\n"
  "  -outlierSwitch <Label> <Size> : change outlier Pvalue to <Label> when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "                                  change outlier Pvalue to <Size> when adjusting sizing (or check for Haplotype indels) during refinement\n"
  "  -endoutlierSwitch <Label> <Size> : change endoutlier Pvalue to <Label> when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "                                     change endoutlier Pvalue to <Size> when adjusting sizing (or check for Haplotype indels) during refinement\n"
  "  -LRbiasSwitch <Label> <Size> : change LRbias value to <Label> when adding/deleting labels (or checking for Haplotype labels) during refinement\n"
  "                                 change LRbias value to <Size> when adjusting sizing (or checking for Haplotype indels) during refinement\n"
  "  -SaveMapset : For each refined contig<N> generate a bnx file PREFIX_contig<N>.bnx with all the input maps used for contig<N> [Default OFF]\n"
  "  -Site_Pen <score> [ <extscore> ] : Score penalty for each consensus site to balance FP and FN rate on consensus maps [Default 3.0]\n"
  "                                     If specified, <extscore> is a lower penalty that is applied during initial extension refinement, before reverting to <score> for final refinement [Default 0.0]\n"
  "  -skip_dist <D> [<Dext>] : Speed up refinement by not checking for new sites more often than every <D> (kb) [Default OFF or 0]\n"
  "                              If <Dext> is specified and larger than <D>, this value is used during initial extension refinement, before reverting to <D> for final refinement [Default 0.500]\n"
  "  -skip_score <S> [<Sext>] : Speed up adding/deleting sites during refinement by NOT rechecking sites that cause Likelihood Score to drop by <S> or more [Default 20.0]\n"
  "                              If <Sext> is specified and smaller than <S>, this value is used during initial extension refinement, before reverting to <S> for final refinement [Default 10.0]\n"
  "  -FN_Refine <FNmul> [<FNmulExt> <FNmin> <FNminExt> ] : Increase FN value during refinement by factor <FNmul>. If <FNmulext> is specified and larger than <FNmul>, this value is used during initial extension refinement, before reverting to <FNmul> [Default OFF]\n"
  "                              If <FNmin> is specified, FN is increased to this absolute value if it is larger than FN times <FNmul>. Similarly <FNminExt> is the minimum absolute value of FN used during initial extension refinement.\n"
  "                              This is intended to model variation in nicking efficiency to avoid having more FN than FP in consensus Map\n"
  "  -MprobevalTest : More thorough refinement (Typically 4x slower) [Default OFF]\n"
  "  -maxContigSiteDensity <D> : If D != 0 : Use -refine 1 -EndTrim 0 for refinment of contigs with Site Density greater than D per 100kb [Default OFF or 0]\n"
  "                              NOTE: Contig will be output unchanged except if -rres was used (NOT YET VALIDATED for refineA)\n"
  "  -biaswtOutlier <WT> : biaswt value to use for outlier regions of query maps, as a multiple of biaswt value. [Default 1.0]\n"
  "  -biaswtRefine <WT> <type> : Apply -biaswt during refinement (If -extend 2 then sites are added/deleted in the extension region before using -biaswt OR -LRbias, since initial score in the extension region will be very low)\n"
  "                            <type> = 0 means only use biaswt when computing final coverage and Chimeric Quality scores, otherwise use as soon as -LRbias is applied [Default 0.0 0]\n"
  "                            NOTE : biaswtRefine is never applied to endoutliers : use with -biaswtEnd 0.0 for consistency between refalign and refinement\n"
  "  -RTHETA_FIX <0/1> : Set to 0 to better match endoutliers before and after refinement\n"
  "  -MinKB <L> : Try to keep labels at least <L> kb apart during refinement. Values between 0 and rres*0.5 can speed up refinement. <L> = 0 implies a value of rres*0.5 [Default 0 OR rres*0.5]\n"
  "-PVendoutlier : Apply a correction to Pvalue (-T) to reflect use of End outliers [Default OFF]\n"
  "-PVres [ <N> ]: Apply a correction to Pvalue (-T) to reflect limited resolution (-res & -bpp). <N> = 2 uses a more aggressive correction. The effect is more significant with high site density [Default OFF]\n"
  "-AlignRes [ <range> ] : Apply more complete scoring of alignments with limited resolution (-res). Results in more alignments if sites are very dense [Default OFF]\n"
  "                        <range> limits range of sites that failed to resolve to (res + resSD*3) * <range>. A Value of 1.0 (or less) disables -AlignRes [Default 2.0]\n"
  "-MinAvCov <min-average-coverage> : Output assembly paths only if average coverage is above this level [Default 4.0]\n"
  "-MinMaps <min-maps> : Output assembly paths only if the number of maps exceeds this number [Default 50]\n"
  "-MinContigLen <len> : Minimum assembly path length in kb (contig length in excess of average map length) [Default 50]\n"
  "-EndTrim <coverage below which consensus ends are trimmed>[Default 8]\n"
  "-TrimCoverage <graph coverage below which input maps are removed>[Default 2.99]\n"
  "-BulgeCoverage <gragh coverage above which weakest edge in inconsistent paths are not broken>[Default 3.9]\n"
  "-CycleCoverage <gragh coverage above which weakest cycle edge is not broken>[Default 19.9]\n"
  "-MaxRelCoverage <multiple> <absolute> [<absolute2>] : If input alignments for any map exceed this <multiple> of the average (or <absolute> level), delete any edges that rank lower than this on both nodes\n"
  "                             connecting them (in order of decreasing -Log(Pvalue)). If <absolute2> is specified (and lower) also delete edges where one node ranks below <absolute2> [Default 100 200 0]\n"
  "-MinCcntmax <N> : Set minimum number of edge confirmations required for at least 1 edge in each direction from each map (otherwise map is deleted)[Default 2]\n"
  "-dumpalign : Dump the alignments and maps remaining after applying -MaxRelCoverage\n"
  "-skipfilter : Assume input alignments are result of -dumpalign and skip -MaxRelCoverage [Default : Off]\n"
  "-align_format <0/1> : 0 means use original text format, 1 means use binary format for faster output of .align files [Default : 0]. NOTE : input of .align detects the format automatically\n"
  "-MaxSR <V> : Force sr values to be no more than <V> at all times [Default 0.50]\n"
  "-MaxSE <V> : Force se values to be no more than <V> at all times [Default 0.00]\n"
  "-PVchim <pvalue> <N> : If an alignment has N or more segments in the middle with score below log(pvalue), the alignment is considered chimeric and ignored (Requres -alignscore)[Default 0.001 3]\n"
  "-mapped-unsplit <0/1> : If 1 : Output from .xmap output  will not split into separate files for each reference contig [Default : 1]\n"
  "-RepeatMask <MaxShift> [ <PvalueRatio> ] : Ignore pairwise alignments which can be shifted by 1 or more sites (up to MaxShift) without decreasing mean-square sizing error Pvalue by more than PvalueRatio [Default 0 0.5]\n"
  "                    NOTE : A smaller PvalueRatio will classify more alignments as Repeats, while PvalueRatio >= 1 is unlikely to classify any alignment as a Repeat\n"
  "-MaxBulge <Len> : Largest bulge (bubble) size to break if lengths match OR one path has coverage below BulgeBreakMaxCov. [Default 2000.0 kb]\n" 
  "   NOTE : Assembler assumes matching lengths imply matching maps so there is a risk of losing alternate Haplotypes if this value is larger than 100kb, but default of 2000kb speeds up Assembly signficantly\n"
  "-FastBulge <N> : speed up alternate path resolution by skippign Dykstra shortest path check for paths with no more than <N> edges[Default 8]\n"
  "-FastPaths <1 or 2> : speed up final longest path (contigs) discovery (1 = conservative, 2 = agressive) : may find fewer (or suboptimial) contigs\n"
  "-maxthreads <N> : Limit number of threads to <N> (0 means use OMP_NUM_THREADS). [Default OMP_NUM_THREADS]\n"
  "-TotalThreads [ <M> <R>] : pro-rate -maxmem and -hashmaxmem by N/M : Assumes -maxmem and -hashmaxmem specify the total usable memory on the machine. <M> defaults to total system threads. [Default OFF]\n"
  "                 NOTE : If N > OMP_NUM_THREADS, memory is pro-rated based on original N value, before N is reduced to OMP_NUM_THREADS\n"
  "                 If R is specified, adjust M to (M*OMP_NUM_THREADS)/R : effectively M/R specifies the amount of thread overloading, the exact values are immaterial.\n"
  "-boostThreads <T> [ <K> ] : After elapsed time of <T> seconds double RefineThreads (or use <K> threads, if specified) subject to original -maxmem and -hashmaxmem. [Default OFF]\n"
  "-minoverlap <ratio> : Ignore pairwise alignments when overlap as fraction of smaller map is below specifed ratio [Default 0]\n"
  "-MaxCoverage <N> : Sets MaxCoverage : If coverage is above this level for a graph edge or node, it is considered reliable and mergeble with other reliable edges [Default 64]\n"
  "-FragilePreserve <0/1> : Try to avoid breaking edges near fragile sites (makes multithreading non-deterministic) [Default 1]\n"
  "-SideBranch <Len> : If <Len> != 0 : Output possible duplicate assembly paths with minimum (center to center) length of at least <Len> kb (and satisfying -MinMaps) [Default 0]\n"
  "  -SideBranchExt : Extend -SideBranch and alternate main paths in both directions up to the next branch point [Default OFF] (NOT YET IMPLEMENTED)\n"
  "-SideChain : Enables hiding side branches which helps simplify the graph : can be slow [Default : OFF]\n"
  "-CovNorm <MapLambda> : Output Normalized coverage distribution (for refined .cmap) assuming map lengths are exponentially distributed with specified mean kB length\n"
  "                       NOTE : Average Map length = MapLambda + MinLen(cutoff) can be estimated from the complete BNX map set using -minlen -merge\n"
  "-XmapStatRead <filename> : Read map statistics from specified filename and use them in place of statistics of actual -i input maps\n"
  "-readLoc : Read simulated data location information from Version 0.1 BNX input file [Default : off]\n"
  "-GraphWeight <PVbase> <MaxWeight> : Weight graph edges (alignments) based on -log10(PValue) and -T threshold : weight = min(MaxWeight,(pow(PVbase,log10(T)-log10(Pvalue))-1)/(PVbase-1)) [Default : OFF]\n"
  "-AlignmentFilter <N> <minlenDelta> <PvalueMult> <M>: If number of alignments per map exceeds <N>, increment -minlen by <minlenDelta> and multiply -T by <PvalueMult> [Default : OFF]\n"
  "                                                  Will also limit total alignments based on -maxmem to M (default 1 million alignments) per Gbyte\n"
  "                                                  Use M <= 1.0, to force checking actual memory usage (VmRss and VmSwap) instead and limiting VmRss+VmSwap to M * maxmem\n"
  "-RAmem <M> : Speed Memory tradeoff in refalign:\n"
  "             M=0 : Maximize speed by doing a single large memory allocation : Can be much slower if memory starts swapping.\n"
  "             M=1 : Adaptive memory allocation, but memory only grows : much less likely to swap than M=0 and almost as fast. [Default: 1]\n"
  "             M=2..7 : Adaptive memory allocation, memory can grow and shrink in steps of (M+1)/M : progressively less memory usage compared to M=1, but more frequent reallocation reduces speed.\n"
  "             M=8 : Re-allocate for each molecule (Minimizes memory usage)\n"
  "-QPmem <M> : Speed Memory tradeoff in qprobeval (see -RAmem for explanation of values)\n"
  "-MPmem <M> : Speed Memory tradeoff in mprobeval (see -RAmem for explanation of values)\n"
;

void SintDetailRA(double X, double Y, int m, int n, int J, int I, int K /* right */, int T /* left */, FLOAT *Yref, double &Bias, double &Pen, double &Gauss, double &PenSm, int verb)
{
  printf("SintDetailRA() called in Assembler : not implemented for Assembler binary\n");
  fflush(stdout);exit(1);
}

int verb=0;

int gargc;
char **gargv;

static void Assembler_cleanup()
{
  fflush(stderr);
  if(VERB){
    static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
    getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
    printf("VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb : CPU time= %0.6f, wall time= %0.6f\n",VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(), wtime());
  }
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
	fprintf(stdout,"fsync(fd) returned -1\n");
	fflush(stdout);
      }
    }
  }
#endif

  if(stderr_file && !(stdout_file && !strcmp(stdout_file,stderr_file))){
    printf("CPU time= %0.6f, wall time= %0.6f\n",mtime(),wtime());
    fprintf(stderr,"END of output\n");
    fflush(stderr);
#ifndef WIN32
    int fd = fileno(stderr);
    if(fd < 0){
      fprintf(stderr,"fileno(stderr) returned -1\n");
      fflush(stderr);
    } else {
      if(fsync(fd) < 0){
	fprintf(stderr,"fsync(fd) returned -1\n");
	fflush(stderr);
      }
    }
#endif
  }

  if(LEAKDEBUG) {/* free up memory */
    if(Gmap){
      for(int i = 0; i < maxmaps; i++){
	if(Gmap[i]){
	  Gmap[i]->allfree();
	  Gmap[i] = NULL;
	}
      }
      delete [] Gmap;
      Gmap = NULL;
    }
    free(output_prefix);
    alignfree();
    align_blockfree();

    extern char *Nickase[MAXCOLOR];/* see input_vfixx.cpp */
    for(register int i = 0; i < colors;i++)
      if(Nickase[i]){
	/*      printf("Nickase[%d]=x%p\n",i,Nickase[i]);
		fflush(stdout);*/
	free(Nickase[i]);
	Nickase[i] = 0;
      }

    graph_free();
    maxmap_free();
    for(int i = 0; spots_filename[i]!=NULL; i++)
      (void)free(spots_filename[i]);
    for(int i = 0; i < num_files; i++)
      free(vfixx_filename[i]);
    for(int i = 0; i < num_align_files; i++)
      free(align_filename[i]);
    if(stderr_file) free(stderr_file);
    if(stdout_file) free(stdout_file);
    
    if(VERB>=2){
      printf("Calling malloc_trim(0) : wtime=%0.6f\n",wtime());
      fflush(stdout);
    }
    malloc_trim(0);
  }
  if(VERB>=2){
    dumpmemmap();      
    fflush(stdout);
  }
}

/** Command line parser and launch assembler
 - parse command line, modify globals, input_bnx(),input_align()
 - graph_build() - initialize the graph object
 - Perform refining, edge removal in the graph with the following methods:
    + graph_nodetrim()
    + graph_edgetrim()
    + graph_chimtrim()
    + graph_redundancy()
    + graph_collapse()
    + graph_bulge()
    + and graph_prune()
 - Output the identified contigs with graph_components()
  */
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

  assert(sizeof(int) == 4);
  assert(sizeof(size_t) == 8);
  assert(sizeof(long long) == 8);
  assert(sizeof(unsigned long long) == 8);

  (void) mtime();
  (void) wtime();

  ZERO_MINKB = 0.0;// change default value for ZERO_MINKB

#if USE_MIC
  /* mark process as a memory hog for OOM killer in case we run out of memory : kill this process first */
  linux_tweak("/proc/self/oom_adj", "15");
  linux_tweak("/proc/self/oom_score_adj", "1000");
#else
  /* try to protect this process from OOM killer (might help with RESCALE, by encouraging OOM killer to go after rogue file transfer job first) */
  linux_tweak("/proc/self/oom_adj", "-17");
  linux_tweak("/proc/self/oom_score_adj", "-1000");
#endif

  #ifdef WIN32
  _set_fmode(_O_BINARY);
  #endif
  
  if(!mallopt(M_MMAP_THRESHOLD, 128*1024)){
    printf("malloc(M_MMAP_THRESHOLD,128*1024) failed\n");
    fflush(stdout);
  }
  if(!mallopt(M_MMAP_MAX, 256*1024)){
    printf("malloc(M_MMAP_MAX,256*1024) failed\n");
    fflush(stdout);
  }

  char *qt;

  gargc = argc;
  gargv = argv;

  /* initialize parameter arrays */
  for(register int i=0;i<MAXCOLOR;i++){
    FP[i] = 1.0;
    FN[i] = 0.15;
    SF[i] = 0.16;
    SD[i] = 0.16;
    SR[i] = 0.0;
  }
  MA_mean = 0.0;
  MA_SF = 1.0;
  MA_SD = 0.10;

  /* parse command line arguments : modifies global parameter variables */
  argv++,argc--;
  while(argc>0){
    if(!strcmp(argv[0],"-help")){
      printf("usage:\n%s\n",usage);
      fflush(stdout);exit(0);
    }
    if(argv[0][0] != '-'){
      printf("invalid arg: \"%s\"(expecting an arg starting with '-')\n",argv[0]);
      if(argc > 1){
	printf("Subsequent args are:");
	for(int i = 1; i < argc; i++)
	  printf(" \"%s\"", argv[i]);
	printf("\n");
      }
      printf("Use -help for list of options\n");
      fflush(stdout);exit(1);
    }

    switch(argv[0][1]){
    case 'v':
      if(!(strcmp(argv[0],"-version"))){
	printversion(stdout);
	fflush(stdout);exit(0);
      }
      break;
    case 'A':
      if(!strcmp(argv[0],"-AlignRes")){
	AlignRes = 1;
	if(argc >= 2 && !(argv[1][0]=='-' && !isdigit(argv[1][1]))){
	  AlignResMult = strtod(argv[1],&qt);
	  if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || AlignResMult < 1.0){
	    printf("%s value has invalid syntax or value (must be >= 1.0):%s\n",argv[0],argv[1]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-AlignmentFilter")){
	if(argc < 4 || argv[1][0]=='-' || argv[2][0]=='-' || argv[3][0]=='-'){
	  printf("%s flag must be followed by 3 numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	AlignmentFilter = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || AlignmentFilter <= 0){
	  printf("%s 1st value is invalid (must be +ve integer):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	minlenDelta = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || minlenDelta <= 0.0){
	  printf("%s 2nd value is invalid (must be +ve number):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	LogPvDelta = strtod(argv[3],&qt);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || LogPvDelta <= 0.0 || LogPvDelta > 1.0){
	  printf("%s 3rd value is invalid (must be +ve number <= 1.0):%s\n",argv[0],argv[3]);
	  fflush(stdout);exit(1);
	}
	LogPvDelta = -log(LogPvDelta)/log(10.0);

	if(argc >= 5 && argv[4][0] != '-'){
	  AlignmentsPerGb = strtod(argv[4],&qt);
	  if(qt == argv[4] || !(*qt==0 || isspace(*qt))){
	    printf("%s 4rd value is invalid (must be non-negative number):%s\n",argv[0],argv[4]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}

	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-A")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	AlignedSiteThreshold = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || AlignedSiteThreshold <= 0 ){
	  printf("%s value is invalid (must be integer > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'B':
      if(!strcmp(argv[0],"-BulgeCoverage")){
	if(argc < 2){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	BulgeBreakMaxCov = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || BulgeBreakMaxCov < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'C':
      if(!strcmp(argv[0],"-CycleRatio")){
	if(argc < 2){
	  printf("%s flag must be followed by a fractional number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CycleCoverageRatio = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || CycleCoverageRatio <= 0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-CycleCoverage")){
	if(argc < 2){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	CycleBreakMaxCov = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || CycleBreakMaxCov < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-ChimRatio")){
	if(argc < 2){
	  printf("%s flag must be followed by a fractional number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ChimRatio = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || ChimRatio <= 0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'E':
      if(!strcmp(argv[0],"-EndTrim")){
	if(argc < 2){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	EndTrimCov = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || EndTrimCov < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'F':
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
      if(!strcmp(argv[0],"-FragilePreserve")){
	if(argc < 2){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	FragilePreserve = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || FragilePreserve < 0 || FragilePreserve > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FastPaths")){
	if(argc < 2){
	  printf("%s flag must be followed by a fractional number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	FastPaths = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || FastPaths < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FastBulge")){
	if(argc < 2){
	  printf("%s flag must be followed by a fractional number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	BulgeFast = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || BulgeFast < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-FP")){
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
	  FP[c] = strtod(argv[c+1],&qt);
	  if(argv[c+1]==qt || !(*qt==0 || isspace(*qt)) || FP[c] < 0.0){
	    if(colors > 1)
	      printf("The %d. %s value %s is not valid (must be non-negative value)\n",c+1,argv[0],argv[c+1]);
	    else
	      printf("The %s value %s is not valid (must be non-negative value)\n",argv[0],argv[c+1]);
	    fflush(stdout);exit(1);
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
	  FN[c] = strtod(argv[c+1],&qt);
	  if(argv[c+1]==qt || !(*qt==0 || isspace(*qt)) || FN[c] < 0.0){
	    if(colors > 1)
	      printf("The %d. %s value %s is not valid\n",c+1,argv[0],argv[c+1]);
	    else
	      printf("The %s value %s is not valid\n",argv[0],argv[c+1]);
	    fflush(stdout);exit(1);
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
      if(!strcmp(argv[0],"-GraphWeight")){
	if(argc < 3 || argv[1][0]=='-' || argv[2][0]=='-'){
	  printf("%s flag must be followed by two numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	GraphWeight = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || !(*qt == 0 || isspace(*qt)) || GraphWeight < 0.0){
	  printf("%s value is invalid (must be >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	GraphMaxWeight = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || !(*qt == 0 || isspace(*qt)) || GraphMaxWeight < 1.0){
	  printf("%s value is invalid (must be >= 1.0):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
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
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || KMAX <= 0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'L':
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
	  printf("%s must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	LocalTheta = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || LocalTheta > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n", argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-LRbiasSwitch")){
	LRbiasSwitch = 1;
	if(!(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-')){
	  printf("%s must be followed by two numbers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	LRbiasLabel = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || LRbiasLabel <= 0.0){
	  printf("%s 1st value of %s is invalid (must be positive)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	LRbiasSize = strtod(argv[2],&qt);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || LRbiasSize <= 0.0){
	  printf("%s 2nd value of %s is invalid (must be positive)\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-L")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	AlignedLengthThreshold = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || AlignedLengthThreshold < 0 ){
	  printf("%s value is invalid (must be number >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-LRbias")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	LRbias = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || !(*qt == 0 || isspace(*qt)) || LRbias < 0.0){
	  printf("%s value is invalid (must be >= 0.0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'M':
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
	argv += 2;
	argc -= 2;
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
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinFP")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
	}
	MinFP = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MinFP < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);
	  exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxFP")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);
	  exit(1);
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
      if(!strcmp(argv[0],"-MaxIntensity")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxIntensity = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxIntensity < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxSR")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxSR = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxSR < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxSE")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxSE = strtod(argv[1],&qt);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || MaxSE < 0.0){
	  printf("%s value is invalid (must be a non-negative float):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MprobevalTest")){
	MprobevalTest = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-MaxCoverage")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by positive integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxCoverage = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxCoverage <= 0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinSplitLen")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinSplitLen = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MinSplitLen <= 0.0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
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
      if(!strcmp(argv[0],"-MinCcntmax")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinCcntmax = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MinCcntmax < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinAvCov")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ContigMinCov = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || ContigMinCov <= 0.0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinMaps")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ContigMinMaps = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || ContigMinMaps <= 0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinContigLen")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ContigMinLen = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || ContigMinLen < 0.0){
	  printf("%s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MaxRelCoverage")){
	if(argc < 3){
	  printf("%s flag must be followed by two integers\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MaxRelCoverage = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MaxRelCoverage <= 0){
	  printf("%s value is invalid (must be integer > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	MaxAbsCoverage = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || MaxAbsCoverage < 0){
	  printf("%s value is invalid (must be integer >= 0):%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	if(argc >= 4 && argv[3][0] != '-'){
	  MaxAbsCoverage2 = strtol(argv[3],&qt,10);
	  if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || MaxAbsCoverage2 < 0){
	    printf("%s value is invalid (must be integer >= 0):%s\n",argv[0],argv[3]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 3;
	argc -= 3;
	continue;
      }
      if(!strcmp(argv[0],"-MinCov")){
	if(argc < 2){
	  printf("%s flag must be followed by an integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinCoverage = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MinCoverage <= 0){
	  printf("%s value is invalid (must be integer > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinEdgeRatio")){
	if(argc < 2){
	  printf("%s flag must be followed by a fractional number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	EdgeCoverageRatio = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || EdgeCoverageRatio <= 0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-MinRatio")){
	if(argc < 2){
	  printf("%s flag must be followed by a fractional number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinCoverageRatio = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MinCoverageRatio <= 0){
	  printf("%s value is invalid (must be > 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'P':
      if(!strcmp(argv[0],"-PVendoutlier")){
	EndOutlierPV = 1;
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
      if(!strcmp(argv[0],"-PVchim")){
	if(argc < 3){
	  printf("%s flag must be followed by 1 float and 1 number\n", argv[0]);
	  fflush(stdout);exit(1);
	}
	Pchim = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || Pchim < 0.0 || Pchim > 1.0){
	  printf("%s pvalue is invalid : must be between 0 and 1:\n%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
        Pchim_minlen = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || Pchim_minlen <= 0){
	  printf("%s 2nd value is invalid : must be +ve integer:\n%s\n",argv[0],argv[2]);
	  fflush(stdout);exit(1);
	}
	argv += 3;
	argc -= 3;
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
      break;
    case 'R':
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
      if(!strcmp(argv[0],"-RepeatMask")){
	if(argc < 2){
	  printf("%s flag must be followed by integer (MaxShift)\n",argv[0]);
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
      if(!strcmp(argv[0],"-Site_Pen")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	Site_Pen = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s 2nd value of %s is invalid\n",argv[0],argv[2]);
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
      if(!strcmp(argv[0],"-SideChain")){
	SideChain = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-SideBranchExt")){
	SidebranchExt = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-SideBranch")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinSidebranch = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || MinSidebranch < 0){
	  printf("%s %s is invalid (must be a non-negative number)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	MaxSidebranch = 1.0e+30;
	if(argc >= 3 && argv[2][0] != '-'){
	  MaxSidebranch = strtod(argv[2],&qt);
	  if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || MaxSidebranch < 0){
	    printf("%s %s is invalid (must be a non-negative number)\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-SaveMapset")){
	SAVE_MAPSET = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-SNR")){
	MapSNR = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-S")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	ScoreThreshold = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value is invalid:%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'T':
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
      if(!strcmp(argv[0],"-TrimCoverage")){
	if(argc < 2){
	  printf("%s flag must be followed by a number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	TrimCoverage = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || TrimCoverage < 0){
	  printf("%s value is invalid (must be >= 0):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-T")){
	if(argc < 2){
	  printf("%s flag must be followed by number\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	double Pvalue = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || Pvalue <= 0.0 || Pvalue > 1.0){
	  printf("%s value is invalid (must be between 0 and 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	LogPvThreshold = -log(Pvalue)/log(10.0);
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'X':
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
      if(!strcmp(argv[0],"-a")){
	if(argc<2){
	  printf("-a flag must be followed by .align filename\n");
	  fflush(stdout);exit(1);
	}
	if(num_align_files >= MAXFILES){
	  printf("Too many input align files (-a) : increase MAXFILES=%d\n",MAXFILES);
	  fflush(stdout);exit(1);
	}
	align_filename[num_align_files++] = strdup(argv[1]);
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-af")){
	if(argc<2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by filelist filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	errno = 0;
	FILE *fp = NULL;
	for(int cnt = 0; cnt <= 65 && (fp = fopen(argv[1],"r")) == NULL; cnt++){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("Cannot read file %s (to get filename list), retry=%d:errno=%d:%s\n",argv[1],cnt,eno,err);
	  printf("Current Time= %s\n", currentDateTime());
	  fflush(stdout);
	  sleep(1);
	  errno = 0;
	}
	if(fp == NULL)
	  exit(1);
	char buf[BUFSIZ];
	int orignumfiles = num_align_files;
	int linecnt = 1;
	for(;fgets(buf,BUFSIZ,fp) != NULL; linecnt++){
	  int len = strlen(buf);
	  if(len >= BUFSIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	    printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",
		   buf[len-1],argv[1],linecnt,buf);
	    fflush(stdout);exit(1);
	  }
	  buf[--len] = '\0';
	  while(len > 0 && isspace(buf[len-1]))
	    buf[--len] = '\0';
	  if(num_align_files < MAXFILES)
	    align_filename[num_align_files++] = strdup(buf);
	}
	fclose(fp);
	if(orignumfiles + linecnt > MAXFILES){
	  printf("Too many input align files (-af) : increase MAXFILES=%d to at least %d\n",
		 MAXFILES, orignumfiles + linecnt);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'b':
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
	if(argc >= 3 && argv[2][0] != '-'){
	  biasWTrefineType = strtol(argv[2],&qt,10);
	  if(qt==argv[2] || !(*qt == 0 || isspace(*qt)) || biasWTrefineType < 0 || biasWTrefineType > 1){
	    printf("%s 2nd value is invalid (must be 0 or 1):%s\n",argv[0],argv[2]);
	    fflush(stdout);exit(1);
	  }
	  argv++;
	  argc--;
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-bnxmergeerror")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a non-negative integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	bnxmergeerror = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt == 0 || isspace(*qt)) || bnxmergeerror < 0){
	  printf("%s value is invalid (must be non-negative):%s\n",argv[0],argv[1]);
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
      if(!strcmp(argv[0],"-bnxerror")){
	bnxerror = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-bpp")){
	if(argc < 2 || argv[1][0]=='-'){
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
      break;
    case 'c':
      if(!strcmp(argv[0],"-contigs_format")){
	if(argc < 2 || argv[1][0]== '-'){
	  printf("%s flag must be followed by 0 or 1\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	contig_format = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || contig_format < 0 || contig_format > 1){
	  printf("%s value is invalid (must be 0 or 1):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-contigs")){
	if(argc < 4){
	  printf("%s must be followed by file name and contig index range\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	contig_filename = argv[1];
	contig_start = strtol(argv[2],&qt,10);
	if(qt == argv[2] || !(*qt==0 || isspace(*qt)) || contig_start <= 0){
	  printf("-contigs start value of %d is invalid:\n%s %s %s %s\n",contig_start,argv[0],argv[1],argv[2],argv[3]);
	  fflush(stdout);exit(1);
	}
	contig_end = strtol(argv[3],&qt,10);
	if(qt == argv[3] || !(*qt==0 || isspace(*qt)) || contig_end <  contig_start){
	  printf("-contigs end value of %d is invalid:\n%s %s %s %s\n",contig_end,argv[0],argv[1],argv[2],argv[3]);
	  fflush(stdout);exit(1);
	}
	contig_start--;
	contig_end--;
	argv += 4;
	argc -= 4;
	continue;
      }
      if(!strcmp(argv[0],"-colors")){
	if(argc < 2){
	  printf("-colors must be followed by a number\n");
	  fflush(stdout);exit(1);
	}
	colors = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || colors <= 0){
	  printf("-colors value is invalid (must be >= 1):%s\n",argv[1]);
	  fflush(stdout);exit(1);
	}
	if(colors > MAXCOLOR){
	  printf("-colors value %d is too large (MAXCOLOR=%d)\n",colors,MAXCOLOR);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'd':
      if(!strcmp(argv[0],"-dumpalign")){
	MaxCoverageDump = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-draftsize")){
	if(argc < 2 || argv[1][0]=='-'){
	  printf("%s flag must be following by integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	SizeType = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || SizeType < 0 || SizeType > 3){
	  printf("%s value of %s is invalid (must be 0,1,2 or 3)\n",argv[0],argv[1]);
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
      break;
    case 'e':
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
	  if(qt == argv[1+PoutlierEndCnt] || !(*qt==0 || isspace(*qt)) || Pvalue < (PoutlierEndCnt > 0 ? PoutlierEndRefine[PoutlierEndCnt-1] : 0.0) || Pvalue > 1.0){
	    printf("%s : %d'th value  of %s is invalid (must be between 0 and 1 and in ascending order)\n",argv[0], PoutlierEndCnt+1, argv[1+PoutlierEndCnt]);
	    fflush(stdout);exit(1);
	  }
	  PoutlierEndRefine[PoutlierEndCnt++] = Pvalue;
	}
	argv += 1+PoutlierEndCnt;
	argc -= 1+PoutlierEndCnt;
	continue;
      }
      if(!strcmp(argv[0],"-extend")){
	if(argc<2){
	  printf("%s flag must be followed by integer between 0 and 2\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	extend = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || extend < 0 || extend > 2){
	  printf("%s value of %s is invalid (must be between 0 and 2)\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-endoutlier")){
	if(argc<2){
	  printf("%s flag must be followed by end outlier(chimerism or blemish) probability per true positive site\n",argv[0]);
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
      if(!strcmp(argv[0],"-endsites")){
	if(argc < 2){
	  printf("-endsites must be followed by a number\n");
	  fflush(stdout);exit(1);
	}
	ENDSITES = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || ENDSITES < 1){
	  printf("-endsites value is invalid (must be >= 1):%s\n",argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'f':
      if(!strcmp(argv[0],"-f") || !strcmp(argv[0],"-force")){
	ForceOverwrite = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-first")){
	if(argc < 2){
	  printf("-first flag must be followed by number\n");
	  fflush(stdout);exit(1);
	}
	Cfirst = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || Cfirst == 0){
	  printf("-first value is invalid (must be != 0):%s\n", argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'h':
      if(!strcmp(argv[0],"-hashtxt")){
	HashTxt = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-hashROC")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s must be following by hashtable filename\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	hash_filename = argv[1];
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'i':
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
	  printf("-if flag must be followed by filelistname\n");
	  fflush(stdout);exit(1);
	}

	errno = 0;
	FILE *fp = NULL;
	for(int cnt = 0; cnt <= 65 && (fp = fopen(argv[1],"r")) == NULL; cnt++){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("Cannot read file %s (to get filename list),retry=%d:errno=%d:%s\n",argv[1],cnt,eno,err);
	  printf("Current Time= %s\n", currentDateTime());
	  fflush(stdout);
	  sleep(1);
	  errno = 0;
	}
	if(fp == NULL)
	  exit(1);
	char buf[BUFSIZ];
	int orignumfiles = num_files;
	int linecnt = 1;
	for(;fgets(buf,BUFSIZ,fp) != NULL; linecnt++){
	  int len = strlen(buf);
	  if(len >= BUFSIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	    printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",
		   buf[len-1],argv[1],linecnt,buf);
	    fflush(stdout);exit(1);
	  }
	  buf[--len] = '\0';
	  while(len > 0 && isspace(buf[len-1]))
	    buf[--len] = '\0';
	  if(num_files < MAXFILES){
#ifndef WIN32
	    if(buf[0] != '/'){/* relative files names are interpreted relative to directory in which -if file is located */
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
	  }
	}
	fclose(fp);
	if(orignumfiles + linecnt > MAXFILES){
	  printf("Too many input files : increase MAXFILES=%d to at least %d\n",
		 MAXFILES, orignumfiles + linecnt);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
    case 'm':
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
      if(!strcmp(argv[0],"-maxInterval")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by non-negative number\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MaxInterval = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || MaxInterval < 0.0){
	  printf("%s value %s is invalid : must be non-negative size in kb\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

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
      if(!strcmp(argv[0],"-minoverlap")){
	if(argc < 2){
	  printf("%s flag must be followed by float\n",argv[0]);
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
      if(!strcmp(argv[0],"-maxthreads")){
	if(argc < 2){
	  printf("%s flag must be followed integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}

	MaxThreads = strtol(argv[1],&qt,10);
	if(qt == argv[1] || !(*qt == 0 || isspace(*qt)) || MaxThreads < 0){
	  printf("%s value %s is invalid : must be integer\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}

	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-maptype")){
	if(argc < 2){
	  printf("%s flag must be followed by number (0 or 1)\n",argv[0]);
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
      if(!strcmp(argv[0],"-minlen")){
	if(argc < 2){
	  printf("%s flag must be followed by minimum -i Map length in kb\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	MinLen = strtod(argv[1],&qt);
	if(qt == argv[1] || !(*qt==0 || isspace(*qt))){
	  printf("%s value has invalid syntax:%s\n",argv[0],argv[1]);
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
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-minsites")){
	if(argc < 2 || argv[0][1] == '-'){
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
    case 'o':
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
	argv += 1 + PoutlierCnt;
	argc -= 1 + PoutlierCnt;
	continue;
      }
      if(!strcmp(argv[0],"-outlierTypeRef")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a positive number\n",argv[0]);
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
      if(!strcmp(argv[0],"-o")){
	if(argc<2){
	  printf("%s flag must be followed by output filename prefix\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	output_prefix = strdup(argv[1]);
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
	if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || Poutlier < 0.0 || Poutlier > 1.0){
	  printf("%s value  of %s is invalid (must be between 0 and 1)\n",argv[0], argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
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
	if(!OUTLIER_MAX && outlierMax < 1e10){
	  printf("ERROR: %s %s not supported : enable OUTLIER_MAX in constants.h and recompile\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
	continue;
      }
      break;
     case 'r':
       if(!strcmp(argv[0],"-rres")){
	 if(argc < 2 || argv[1][0] == '-'){
	   printf("%s flag must be followed by resolution in pixels\n",argv[0]);
	   fflush(stdout);exit(1);
	 }
	 rres = strtod(argv[1],&qt);
	 if(qt == argv[1] || !(*qt==0 || isspace(*qt)) || rres < 0.0){
	   printf("-rres value has invalid syntax (must be non-negative number):%s\n",argv[1]);
	   fflush(stdout);exit(1);
	 }
	 argv += 2;
	 argc -= 2;
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
       if(!strcmp(argv[0],"-readLoc")){
	 HasLoc = 1;
	 argv++;
	 argc--;
	 continue;
       }
       if(!strcmp(argv[0],"-refine")){
	if(argc<2){
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
       if(!strcmp(argv[0],"-ref")){
	if(argc<2){
	  printf("-ref flag must be followed by spots filename\n");
	  fflush(stdout);exit(1);
	}
	int i=0;
	while(spots_filename[i] != NULL)
	  i++;
	if(i >= MAXFILES){
	  printf("Too many input files : increase MAXFILES=%d to at least %d\n",
		 MAXFILES, i+1);
	  fflush(stdout);exit(1);
	}
	spots_filename[i] = strdup(argv[1]);
	spots_filename[i+1] = NULL;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-reff")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("-reff flag must be followed by filelistname\n");
	  fflush(stdout);exit(1);
	}
	int num_reffiles = 0;
	while(spots_filename[num_reffiles] != NULL)
	  num_reffiles++;

	FILE *fp;
	if((fp = fopen(argv[1],"r"))==NULL){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("Cannot read file %s (to get filename list for %s):errno=%d:%s\n",argv[1],argv[0],eno,err);
	  fflush(stdout);exit(1);
	}
	char buf[BUFSIZ];
	int orignumfiles = num_reffiles;
	int linecnt = 1;
	for(;fgets(buf,BUFSIZ,fp) != NULL; linecnt++){
	  int len = strlen(buf);
	  if(len >= BUFSIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	    printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",
		   buf[len-1],argv[1],linecnt,buf);
	    fflush(stdout);exit(1);
	  }
	  buf[len-1] = '\0';
	  if(num_reffiles < MAXFILES)
	    spots_filename[num_reffiles++] = strdup(buf);
	}
	fclose(fp);
	if(orignumfiles + linecnt > MAXFILES){
	  printf("Too many input files : increase MAXFILES=%d to at least %d\n",
		 MAXFILES, orignumfiles + linecnt);
	  fflush(stdout);exit(1);
	}
	spots_filename[num_reffiles] = NULL;
	argv += 2;
	argc -= 2;
	continue;
      }
      if(!strcmp(argv[0],"-refmap")){
	if(argc < 2){
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
      if(!strcmp(argv[0],"-skipfilter")){
	SkipFilter = 1;
	argv++;
	argc--;
	continue;
      }
      if(!strcmp(argv[0],"-stdout")){
	char filename[PATH_MAX];
	if(argc < 2 || argv[1][0] == '-'){/* 0 args */
	  if(!output_prefix){
	    printf("-o option must be specified before %s with no args\n",argv[0]);
	    fflush(stdout);exit(1);
	  }
	  if(contig_filename)
	    sprintf(filename,"%s_id%d.stdout", output_prefix,contig_start+1);
	  else
	    sprintf(filename,"%s.stdout", output_prefix);
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
	    printf("-o options must be specified before %s with no args\n",argv[0]);
	    fflush(stdout);exit(1);
	  }
	  if(contig_filename)
	    sprintf(filename,"%s_id%d.stdout", output_prefix, contig_start+1);
	  else
	    sprintf(filename,"%s.stdout", output_prefix);
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
      if(!strcmp(argv[0],"-sm")){
	if(argc < 2){
	  printf("%s must be followed by float\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	SM = strtod(argv[1],&qt);
	if(argv[1]==qt || SM < 0.0){
	  printf("The %s value %s is not valid\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
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
      break;
    case 'u':
      if(!strcmp(argv[0],"-usecolor")){
	if(argc < 2 || argv[1][0] == '-'){
	  printf("%s flag must be followed by a +ve integer\n",argv[0]);
	  fflush(stdout);exit(1);
	}
	usecolor = strtol(argv[1],&qt,10);
	if(qt==argv[1] || !(*qt==0 || isspace(*qt)) || usecolor <= 0 || usecolor > 2){
	  printf("%s value is invalid (must be 1 or 2):%s\n",argv[0],argv[1]);
	  fflush(stdout);exit(1);
	}
	argv += 2;
	argc -= 2;
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

  deltaExtXRef = max(DELTA_X,max(deltaExtXRef,deltaXRef));
  deltaExtYRef = max(DELTA_Y,max(deltaExtYRef,deltaYRef));
  if(!outlierLambdaRef)
    outlierLambdaRef = outlierLambda;
  if(!outlierLambdaLabel)
    outlierLambdaLabel = outlierLambdaRef;
  if(!outlierLambdaSize)
    outlierLambdaSize = outlierLambdaRef;
  if(!OutlierTypeLabel)
    OutlierTypeLabel = OutlierType;
  if(!OutlierTypeSize)
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

  if(AlignmentFilter > 0 && MinLen <= 0.0){
    printf("WARNING: -AlignmentFilter used without specifying -minlen\n");
    fflush(stdout);
  }
  if(AlignmentFilter > 0 && minlenDelta <= 0.0 && LogPvDelta <= 0.0){
    printf("Cannot use -AlignmentFilter without specifying a change in minlen or Pvalue\n");
    fflush(stdout);exit(1);
  }

  MapSNROrig = MapSNR;
  MapStitchedOrig = MapStitched;
  MapIntensityOrig = MapIntensity;
  MapStitchLocOrig = MapStitchLoc;

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

#ifndef WIN32
  if(MaxMem > 0.0){
    double GB = 1024.0 * 1024.0 * 1024.0;

    /* call "free -m" to see how much actual swap space is available : reduce MaxVirtMem to this value + MaxMem */
    long long MemSize, SwapSize, AvailableMem = 0;
#if USE_MIC
    getswap(MemSize,SwapSize);/* Total Mem, free Swap (and Available Mem) based on "free -m" */
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

    //    MemSize += VmRSS;
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
    if(MemSize < (MaxMem + MaxMemIncrease)*GB){
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
#endif

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

  if(!num_files){
    printf("No map input file specified\n");
    fflush(stdout);exit(1);
  }
  if(!output_prefix){
    printf("-o option must be specified\n");
    fflush(stdout);exit(1);
  }
  draft_prefix = output_prefix;

  if(!colors){
    printf("-colors parameter missing\n");
    fflush(stdout);exit(1);
  }
  if(Refine){
    for(int c = 0; c < colors; c++)
      if(resSD[c] <= 0.0){
	printf("-refine not supported with resSD[%d]= 0.0 : changed to 0.001\n", c+1);
	fflush(stdout);
	resSD[c] = 0.001;
      }
  }
  if(RepeatMaxShift > 0 && !AlignScore){
    printf("-RepeatMask cannot be used without -alignscore\n");
    fflush(stdout);exit(1);
  }
  
// -bpp adjustment has been moved to mapcorrect() in input_vfixx.cpp

  if(Cfirst < 0 && -Cfirst > num_files){
    printf("-first %d not compatible with %d input files\n",Cfirst,num_files);
    fflush(stdout);exit(1);
  }

  for(int c = 0; c < colors; c++)
    XmapCount[c] = 0;

  int bnxfile = 1;/* check if all input files are .bnx files */
  for(register int i = 0; i < num_files; i++){
    char *filename = vfixx_filename[i];
    size_t len = strlen(filename);
    if(len >= strlen(".bnx") && !strcmp(&filename[len-strlen(".bnx")],".bnx"))
      continue;
    bnxfile = 0;
  }

  PixelLen = bpp;

  if(!(AVOID_BNX && contig_filename && bnxfile)){
    /* read input maps (vfixx or .cmaps format) */
    startmaps = 0;
    for(register int i = 0; i < num_files; i++){
      startmaps += input_vfixx(i,nummaps,maxmaps,Gmap);
      if(Cfirst < 0 && i+1 == -Cfirst)
	Cfirst = startmaps;
    }

    if(DEBUG>=2){
      if(MapSNR != MapSNROrig){/* check to make sure all maps have SNR information : if not reset MapSNR = 0 */
	if(DEBUG) assert(MapSNROrig==0);
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
      if(MapStitched != MapStitchedOrig){
	if(DEBUG) assert(MapStitchedOrig==0);
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
      if(MapIntensity != MapIntensityOrig){
	if(DEBUG) assert(MapIntensityOrig==0);
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
    if(MapStitchLoc != MapStitchLocOrig){
      if(DEBUG) assert(MapStitchLocOrig==0);
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

    input_vfixx_duplicates();

    YYmap = XXmap = Gmap;
    numX = numY = nummaps;

    if(DEBUG) assert(PixelLen > 0.0);
    if(DEBUG) assert(origPixelLen > 0.0);

    if(VERB){
      if(num_files > 1)
	printf("Read %d Maps from %d files %s %s ... (total maps=%d)\n",startmaps, num_files, vfixx_filename[0],vfixx_filename[1], nummaps);
      else
	printf("Read %d Maps from %s (total maps=%d)\n",startmaps,vfixx_filename[0],nummaps);
      if(VERB>=2)
	dumpmemmap();
      fflush(stdout);
    }
  }

  if(!(AVOID_BNX && contig_filename && bnxfile)){      /* Handle -usecolor */
    if(usecolor && (colors==2 || usecolor > 1)){
      if(colors == 1 && usecolor > 1){
	printf("ERROR: -usecolor %d is not valid with single color -i input maps\n",usecolor);
	fflush(stdout);exit(1);
      }
      if(DEBUG) assert(colors==2);
      if(usecolor==2){/* swap 1st and 2nd color information */
	for(int i = 0; i < startmaps; i++)
	  Gmap[i]->colorswap(usecolor);
	char *tmp = Nickase[0]; Nickase[0] = Nickase[1]; Nickase[1] = tmp;
      }
      
      colors = colors_set = 1;
    } else if(usecolor){/* redundant use of -usecolor 1 with only a single color */
      if(DEBUG) assert(colors==1);
      usecolor = 0;/* otherwise it will assume there were originally 2 colors and try to write them out when needed */
    }
  }

  if(refmap_filename)
    input_refmap(refmap_filename);/* read in reference alignments for debugging */

  if(spots_filename[0] && !contig_filename){
    for(int c = 0; c < colors; c++)
      if(!res[c]){
	printf("res[%d] parameter must be > 0 (required to condense reference map data and specify data resolution) : changing to 0.001\n",c+1);
	fflush(stdout);
	res[c] = 0.001;
      }

    /* read reference map (spots format) */
    for(int i=0; spots_filename[i] != NULL; i++)
      (void)input_spots(spots_filename[i], numrefmaps, maxrefmaps, refmap, i);

    if(numrefmaps <= 0){
      printf("No reference maps read in\n");
      fflush(stdout);exit(1);
    }
  }

  if(contig_filename){/* just refine previously saved contigs and exit */
    if(!Refine){
      printf("-contigs requires -refine (>0)\n");
      fflush(stdout);exit(1);
    }
    if(spots_filename[0]){
      printf("reference map %s (and any others) ignored (due to -contigs)\n",spots_filename[0]);
      fflush(stdout);
    }
    if(num_align_files > 0){
      printf("%d pairwise alignment files ignored (due to -contigs)\n",num_align_files);
      fflush(stdout);
    }

    Ccontig *contig = 0;
    int numcontigs = 0;

    if(!RefineOverwrite){
      char filename[PATH_MAX];
      FILE *fp;


      long long i = contig_start;
      for(; i <= contig_end; i++){
	sprintf(filename,"%s_contig%lld.cmap",output_prefix,i+1);
	if(VERB/* HERE HERE >=2 */){
	  printf("Checking existence of %s\n",filename);
	  fflush(stdout);
	}

	if((fp = fopen(filename,"r")) != NULL){
	  fclose(fp);
	  if(VERB/* HERE HERE >=2 */){
	    printf("WARNING: Output file %s already exists\n",filename);
	    fflush(stdout);
	  }
	  continue;
	}
	break;
      }
      if(i > contig_end){/* All contigs exist */
	if(VERB){
	  printf("WARNING: All contigs N = %d .. %d exist as %s_contigN.cmap : skipping refineA\n",contig_start+1,contig_end+1,output_prefix);
	  fflush(stdout);
	}
	goto Lfree;
      }
    }

    input_contigs(contig,numcontigs,contig_filename);
    if(contig_end >= numcontigs){
      printf("contig index range reduced to %d .. %d since there are only %d total contigs\n",
	     contig_start,numcontigs-1,numcontigs);
      fflush(stdout);
      contig_end = numcontigs-1;
    }

    if(VERB>=2){
      printf("After reading contigs: contig_start=%d,contig_end=%d:\n",contig_start,contig_end);
      for(int i = 0; i < numcontigs; i++){
	Ccontig *p = &contig[i];
	if((max(0,contig_start-1) <= i && i <= contig_end+2) || p->nummaps > 0)
	  printf("contig[%d]=%p: mapid=%d, nummaps=%d, contig= %p\n",i, p, p->mapid, p->nummaps, p->contig);
      }
      fflush(stdout);
    }

    if(!(AVOID_BNX && bnxfile)){
      if(colors==1 && !usecolor)/* HERE : check why this segfaults when colors==2 */
	mapcompact(nummaps, maxmaps, Gmap);
    } else {/* read in only the bnx maps needed */
      /* assemble the set of map ids needed */
      std::set<long long> ids;
      for(int C = contig_start; C <= contig_end; C++){
	Ccontig *Ycontig = &contig[C];
	int MD = Ycontig->nummaps;
	Ccontig *Xcontig = Ycontig->contig;
	for(int m = 0; m < MD; m++){
	  if(VERB>=2){
	    printf("C=%d,m=%d/%d:Adding contig[C]->contig[m].id = %lld to idset\n",
		   C,m,MD,Xcontig[m].id);
	    fflush(stdout);
	  }
	  ids.insert(Xcontig[m].id);
	}
	if(VERB/* HERE >=2 */){
	  printf("C=%d:number of bnx ids to be read in = %lu\n",C,ids.size());
	  fflush(stdout);
	}
      }

      /* read input bnx maps selectively */
      extern int input_bnx(int fileid, int &nummaps, int &maxmaps, Cmap **&Map, std::set<long long> *ids);
      startmaps = 0;
      for(register int i = 0; i < num_files; i++)
	startmaps += input_bnx(i,nummaps,maxmaps,Gmap, &ids);
      
      input_vfixx_duplicates();

      if(VERB>=2){/* print map id 4155456 */
	for(int m = 0; m < nummaps; m++){
	  Cmap *pmap = Gmap[m];
	  if(pmap->id == 4155456LL && pmap->SNR[0]){
	    int N = pmap->numsite[0];
	    double *SNR = pmap->SNR[0];
	    FLOAT *Y =  pmap->site[0];
	    printf("After input_vfixx_duplicates: id=%lld,N=%d:\n",pmap->id,N);
	    for(int J = 1; J <= N+1; J++)
	      printf(" Y[%d]=%0.4f, SNR=%0.4f\n",J,Y[J], (J <= N) ? SNR[J] : -1.0);
	    fflush(stdout);
	  }
	}
      }

      YYmap = XXmap = Gmap;
      numX = numY = nummaps;

      if(DEBUG) assert(PixelLen > 0.0);
      if(DEBUG) assert(origPixelLen > 0.0);

#if 0 // -bpp adjustment has been moved to mapcorrect() in input_vfixx.cpp
      if(!PixelLen){
	printf("Default BasesPerPixel = 500 used\n");
	fflush(stdout);
	PixelLen = 0.500;/* default value when no -i input is used OR all -i inputs are .bnx files */
      }
      origPixelLen = PixelLen;

      if(bpp > 0.0 && bpp != PixelLen){
	register double C = bpp/PixelLen;
	for(register int i= 0; i < startmaps; i++){
	  register Cmap *pmap = Gmap[i];
	  for(register int c=0; c < colors; c++){
	    register FLOAT *X = pmap->site[c];
	    register int N = pmap->numsite[c];
	    for(register int j= N+1; j > 0; j--)
	      X[j] *= C;
	  }
	}
	if(VERB){
	  printf("Adjusted PixelLen from %0.2f to %0.2f(bp)\n",PixelLen*1000.0,PixelLen*C*1000.0);
	  fflush(stdout);
	}
	MinLen *= C;
	PixelLen *= C;
      }
#endif // old code

      assert(fabs(PixelLen - bpp) < PixelLen*0.0001);

      if(VERB){
	if(num_files > 1)
	  printf("Read %d Maps from %d files %s %s ... (total maps=%d)\n",startmaps, num_files, vfixx_filename[0],vfixx_filename[1], nummaps);
	else
	  printf("Read %d Maps from %s (total maps=%d)\n",startmaps,vfixx_filename[0],nummaps);
	fflush(stdout);
      }

      /* Handle -usecolor */
      if(usecolor && (colors==2 || usecolor > 1)){
	if(colors == 1 && usecolor > 1){
	  printf("ERROR: -usecolor %d is not valid with single color -i input maps\n",usecolor);
	  fflush(stdout);exit(1);
	}
	if(DEBUG) assert(colors==2);
	if(usecolor==2){/* swap 1st and 2nd color information */
	  for(int i = 0; i < startmaps; i++)
	    Gmap[i]->colorswap(usecolor);
	  char *tmp = Nickase[0]; Nickase[0] = Nickase[1]; Nickase[1] = tmp;
	}
      
	colors = colors_set = 1;
      } else if(usecolor){/* redundant use of -usecolor 1 with only a single color */
	if(DEBUG) assert(colors==1);
	usecolor = 0;/* otherwise it will assume there were originally 2 colors and try to write them out when needed */
      }

      if(colors==1 && !usecolor)/* HERE : check why this segfaults when colors==2 */
	mapcompact(nummaps,maxmaps,Gmap);

      if(VERB>=2){
	printf("Renumbering %d map ids in contig[%d..%d]: wall time=%0.6f\n",nummaps, contig_start,contig_end, wtime());
	fflush(stdout);
      }

      /* renumber the mapid's in the contigs */
      /* HERE : switch to using HashFindID() : see Assembler_input.cpp */
      delete [] id2mapid;
      id2mapid = new Cid2mapid[nummaps];
      for(int i = 0; i < nummaps; i++){
	register Cid2mapid *p = &id2mapid[i];
	p->id = Gmap[i]->id;
	p->mapid = i;
      }
      qsort(id2mapid,nummaps,sizeof(Cid2mapid),(intcmp *)idinc);

      for(int C = contig_start; C <= contig_end; C++){
	register Ccontig *Ycontig = &contig[C];
	int MD = Ycontig->nummaps;
	register Ccontig *Xcontig = Ycontig->contig;
	for(register int m = 0; m < MD; m++){
	  long long id = Xcontig[m].id;
	  int index = findid(id,id2mapid,nummaps);
	  if(index < 0){
	    printf("contig[%d].Xcontig[%d]:id=%lld not found in %d maps\n",C,m,id,nummaps);
	    fflush(stdout);exit(1);
	  }
	  Xcontig[m].mapid = index;
	}
      }

      if(VERB>=2){
	printf("Renumbered %d map ids in contig[%d..%d]: wall time=%0.6f\n",nummaps, contig_start,contig_end, wtime());
	fflush(stdout);
      }

      if(DEBUG){/* check that numsites match in contig[C].contig[m] and Gmap[mapid] */
	for(int C = contig_start; C <= contig_end; C++){
	  Ccontig *Ycontig = &contig[C];
	  int MD = Ycontig->nummaps;
	  Ccontig *Xcontig = Ycontig->contig;
	  for(int m = 0; m < MD; m++){
	    int mapid = Xcontig[m].mapid;
	    Cmap *pmap = Gmap[mapid];
	    if(DEBUG) assert(pmap->id == Xcontig[m].id);
	    for(int c = 0; c < colors; c++){
	      if(DEBUG) assert(0 <= Xcontig[m].trimL[c] && Xcontig[m].trimL[c] <= Xcontig[m].trimR[c]);
	      if(DEBUG) assert(Xcontig[m].numsite[c] == Xcontig[m].trimR[c] - Xcontig[m].trimL[c] - 1);
	      if(Xcontig[m].trimR[c] > pmap->numsite[c] + 1){
		printf(" mismatch between %s and %s:\n",contig_filename,vfixx_filename[pmap->fileid]);
		printf("  contig=%d, m=%d: mapid=%d,id=%lld,numsite[c=%d]=%d,trimL[c]=%d,trimR[c]=%d (Gmap[mapid]->numsite[c]=%d)\n",
		       C,m,mapid,Xcontig[m].id,c,Xcontig[m].numsite[c],Xcontig[m].trimL[c],Xcontig[m].trimR[c],pmap->numsite[c]);
		printf("  Make sure PairWise alignment, Assembler and refineA all use the same -mres -minsites -maxsites -minSNR -minlen -maxlen\n");
		printf("  Make sure input BNX files are the same for PairWise alignment and refineA\n");
		fflush(stdout);

		printf("  Gmap[mapid=%d]->site[c=%d][]:\n",mapid,c);
		for(int i = 1; i <= pmap->numsite[c] + 1; i++)
		  printf("\t i= %d: site= %0.6f\n", i, pmap->site[c][i]);
		/*		printf("  Xcontig[m=%d].site[c=%d][]:\n", m, c);
		for(int i = 1; i <= Xcontig[m].numsite[c] + 1; i++)
		printf("\t i= %d: site= %0.6f\n", i, Xcontig[m].site[c][i]);*/
		fflush(stdout);exit(1);
	      }
	    }
	  }
	}
      }


    }// (AVOID_BNX && bnxfile)
    
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
    
      int linecnt = 1,len,cnt, found = 0;
      for(;fgets(buf,LINESIZ,fp) != NULL; linecnt++){
	if((len = strlen(buf)) >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	  printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
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
    }

    if(parametersfile && ResBinsFile && ResBins[0] > 0)/* initialize rawsite[] and apply bias correction to site[] */
      BiasCorrect(Gmap,0,startmaps, 1);

    for(register int i = contig_start; i <= contig_end; i++){
      if(VERB){
	printf("Refining contig[%d]:numsite=%d\n",i,contig[i].numsite[0]);
	fflush(stdout);
      }
      gmap = Gmap;
      if(extend){
	contig[i].left = contig[i].site[0][0];
	contig[i].right = contig[i].site[0][contig[i].numsite[0]];
      }
      contig[i].totalsites = -1;

      if(SAVE_MAPSET){/* save mapset as <PREFIX>_contig<N>.bnx */
	Cmap **goodmaps = new Cmap *[contig[i].nummaps];
	for(int k = 0; k < contig[i].nummaps; k++)
	  goodmaps[k] = Gmap[contig[i].contig[k].mapid];
	char prefix[PATH_MAX];
	sprintf(prefix,"%s_contig%d",output_prefix,i+1);
	output_bnx(prefix,goodmaps,0,contig[i].nummaps-1,0,NULL,NULL,-1);
	delete [] goodmaps;
      }

      int Lfrozen = 0, Rfrozen = 0;
      contig[i].id = i+1;
      refine(&contig[i],1,Lfrozen,Rfrozen);

      /* output draft consensus for contig[i] (as .cmap file) */
      output_draft(&contig[i],i+1,output_prefix,0,0,0);
    }

    /* free up memory */
  Lfree:
    if(VERB>=2){
      printf("Freeing Gmap[] and contig[] memory:maxmaps=%d,nummaps=%d,numcontigs=%d\n",
	     maxmaps,nummaps,numcontigs);
      fflush(stdout);
    }

    if(LEAKDEBUG){
      if(contig){
	delete [] contig;
	contig = NULL;
      }
      if(id2mapid){
	delete [] id2mapid;
	id2mapid = NULL;
      }
    }
    if(VERB>=2){
      dumpmemmap();      
      fflush(stdout);
    }

    Assembler_cleanup();

    exit(0);
  }

  if(colors==1 && !usecolor)/* HERE : check why this segfaults when colors==2 */
    mapcompact(nummaps, maxmaps, Gmap);

  if(parametersfile && ResBinsFile && ResBins[0] > 0)/* initialize rawsite[] and apply bias correction to site[] */
    BiasCorrect(Gmap,0,startmaps,1);

  if(VERB>=2){
    dumpmemmap();
    fflush(stdout);
  }

  if(!num_align_files){
    printf("No .align input file specified\n");
    fflush(stdout);exit(1);
  }
  input_align(num_align_files,align_filename);

  if(VERB>=2){
    dumpmemmap();
    fflush(stdout);
  }

  if(hash_filename){/* display ROC curve and exit */
    hashROC(hash_filename);
    exit(0);
  }

  Cmap **origmap = Gmap;

  try {
    /* Now build graph from set of maps and alignments */
    /* results are in node[0..numnodes-1], see Assembler_graph.h */
#if CALIGN_SMALL   
    if(DEBUG) assert(alignmentA != NULL && alignment == NULL);
    graph_build(Gmap,nummaps,alignmentA, numaligns);
#else
    if(DEBUG) assert(alignment != NULL && alignmentA == NULL);
    graph_build(Gmap,nummaps,alignment, numaligns);
#endif
    //  graph_edgecnt(0);
    //  graph_display(1,1,(char *)"complete");
    //    exit(0);

    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    if(MaxCoverageDump){
      extern void DumpMapAlignments();
      DumpMapAlignments();
      exit(0);
    }

    /* trim map ends and detect chimeric maps based on local coverage from pairwise alignments */
    graph_nodetrim(0);
    //  graph_edgecnt(0);
    //  graph_display(1,1,(char *)"nodetrim");
    //  fflush(stdout);exit(1);

    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    // save alignment[] data to disk
    double stime = wtime();
    if(VERB){
      printf("Swapping out %lu alignments in %d blocks : wall= %0.6f\n",numaligns,num_align_blocks,stime);
      fflush(stdout);
    }
    size_t totmemsiz = 0;
    for(int i = 0; i < num_align_blocks; i++)
      totmemsiz += align_block[i]->swap(output_prefix,i);
    if(VERB){
      double etime = wtime();
      printf("Swapping out %lu alignments (%0.4f Gb, %0.2f Mb/sec, time= %0.6f) : wall= %0.6f\n",numaligns,totmemsiz * 1e-9, totmemsiz * 1e-6 / (etime - stime), etime-stime,etime);
      fflush(stdout);
    }

    // HERE HERE : save Gmap[] data to disk
    Gmap = NULL;
    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    // expand CnodeA -> Cnode, CedgeA -> Cedge (if needed)
    graph_convert();
    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    /* trim spurious edges from graph, based on mutual connectivity counts of the two nodes */
    graph_edgetrim(0);
    //  graph_edgecnt(0);
    //  graph_display(1,1,(char *)"edgetrim");
    //  fflush(stdout);exit(1);
    
    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    for(int i = 0; i < 100;i++){
      if(VERB){
	printf("i=%d:",i);
	fflush(stdout);
      }
      if(!graph_chimtrim(0))
	break;
    }

    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    /* free up memory allocated in graph_chimtrim() */
    extern int *edgecnt_old, *edgecnt_new;
    if(edgecnt_old){
      delete [] edgecnt_old;
      delete [] edgecnt_new;
      edgecnt_old = edgecnt_new = 0;
    }
    extern int **Pmark, maxpid;
    extern Cnode ***Pparent;
    if(Pmark){
      for(register int tid = 0; tid < maxpid; tid++){
	if(Pmark[tid]){
	  delete [] Pmark[tid];
	  Pmark[tid] = 0;
	  delete [] Pparent[tid];
	  Pparent[tid] = 0;
	}
      }
      delete [] Pmark;
      delete [] Pparent;
      Pmark = 0;
    }
    //    graph_display(1,1,(char *)"chimtrim");
    //  fflush(stdout);exit(1);

    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    /* remove redundant edges and update coverage(redundancy) 
       estimates of remaining edges and nodes (Algorithm 4.1.1 in Nyugen's thesis) */
    graph_redundancy(0);
    //    graph_display(1,1,(char *)"redundancy");
    //  fflush(stdout);exit(1);

    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    // compact/reallocate graph structure to save memory
    graph_compact();

    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    double origEdgeCoverageRatio = EdgeCoverageRatio;
    int origMinCoverage = MinCoverage;
    EdgeCoverageRatio  /= 4;
    MinCoverage = max(1,MinCoverage/2);  
    // WAS  MinCoverage = min(min(2,MinCoverage),max(1,MinCoverage/2));

    if(VERB>=2){
      printf("origMinCoverage=%d,MinCoverage=%d\n",origMinCoverage,MinCoverage);
      fflush(stdout);
    }

    graph_prune(0,0);/* get rid of weak edges and chimeric nodes right away */
    //  graph_display(1,1,(char *)"prune1");
    //  exit(0);

    /* Repeat Chain Collapsal and Graph Bulges until no further progress is possible */
    int progress = 3;
    double origMaxBulge = MaxBulge;
    int bulge = -1;
    int k = 0, kmax = max(999,(int)floor(MaxBulge/100.0 + 0.5));
    for(register int j=0; j < 1000;j++){  
      MaxBulge = min(origMaxBulge,(k+1)*100.0);
      int BFSdepth = k+1;
      if(VERB){
	printf("j=%d:k=%d,MaxBulge=%0.1f,BFSdepth=%d,MinCoverage=%d\n",j,k,MaxBulge,BFSdepth,MinCoverage);
	fflush(stdout);
      }

      /* Represent linear chains as single Super-Edges (with pointers to the start/end of the original chain)
	 (Chain Collapsal in Nygen's thesis) Coverage of Super-Edges is based on average coverage of original chain,
	 which must lie within the same coverage range (one of : 1-1.9, 2-3.9, 4-7.9, 8-15.9, 16-31.9, 32+) */
      if(!graph_collapse(MinCoverage,0) && k >= 5){
	if(--progress <= 0)
	  break;
      } else
	progress = 3;
      //    if(j==0) graph_display(1,1,"collapse1");

      /* Check alternate parallel paths for consistency and delete path with lower coverage */
      if(!(bulge = graph_bulge(BFSdepth, 0)) && MaxBulge >= origMaxBulge - 0.001){
	if(--progress <= 0)
	  break;
      } else
	progress = 3;
      if(bulge < 500 && k < kmax)
	k++;

      /*    if(j==4){
	    graph_display(1,1);
	    //      graph_components();
	    exit(0);
	    }*/

      /* prune weak edges, chimeric nodes and isolated nodes */
      if(!graph_prune(SideChain,0) && k >= 4){
	if(--progress <= 0)
	  break;
      } else
	progress = 3;

    }

    //  graph_display(1,1,(char *)"After_loop");
    //  graph_bulge(1000,1);
    //  fflush(stdout);exit(1);

    EdgeCoverageRatio = origEdgeCoverageRatio;
    MinCoverage = origMinCoverage;

    graph_prune(SideChain,0);
    //    graph_display(1,1,(char *)"final_prune");
    //  fflush(stdout);exit(1);

    graph_collapse(MinCoverage,0);

    // if (Gmap[0]->startloc > 0.0)
    //  graph_display(1,1,(char *)"final");

    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }
    
    if(VERB){
      printf("Calling malloc_trim(0) : wtime=%0.6f\n",wtime());
      fflush(stdout);
    }
    malloc_trim(0);
    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    // HERE HERE : restore Gmap[] and alignment[] data from disk
    Gmap = origmap;
    stime = wtime();
    if(VERB){
      printf("Swapping in %lu alignments: wall= %0.6f\n",numaligns,stime);
      fflush(stdout);
    }
    totmemsiz = 0;
    for(int i = 0; i < num_align_blocks; i++)
      totmemsiz += align_block[i]->restore(i);
    if(VERB){
      double etime = wtime();
      printf("Swapping in %lu alignments (%0.4f Gb, %0.2f Mb/sec, time= %0.6f) : wall= %0.6f\n",numaligns,totmemsiz * 1e-9, totmemsiz * 1e-6 / (etime - stime), etime-stime,etime);
      fflush(stdout);
    }

    if(VERB>=2){
      dumpmemmap();
      fflush(stdout);
    }

    /* locate connected components and the longest paths in them */
    graph_components();

    if(VERB>=2){
      dumpmemmap();      
      fflush(stdout);
    }

    Assembler_cleanup();
  } catch (exception &e){
    cout << e.what() << endl;
    printf("Assembler(): exception thrown in Assembler_graph.cpp\n");
    dumpmemmap();    
    fflush(stdout);
    assert(0);
  }

  exit(0);
}
