#include "globals.h"
#include "parameters.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/parameters.cpp 11543 2020-08-26 18:53:06Z tanantharaman $");

/* parameters shared by RefAligner and Assembler */

double MININTERVAL= 0.020; /* minimum interval size in kb (used for end intervals of input maps and split maps) */

const char *tmpfs = "/scratch";/* use this folder (if it exists) for temporary files. Should be a fast filesystem (at least 750 Mb/sec) */

int num_files = 0;/* number of -i input files */
char *vfixx_filename[MAXFILES];/* input -i filenames (has -ref filenames appended to end starting at vfixx_filename[num_files]) */
char *spots_filename[MAXFILES+1]={NULL,NULL};/* input -ref file name */
int num_reffiles = 0;/* number of input -ref files */
char *output_prefix = 0;/* output file prefix (already includes _contigN with N derived from -id, if present) */
char *draft_prefix = 0;/* draft output file prefix (_contigN to be added for each contig output if -id is specified OR Assembler contig id is known) */
char *output_filter = NULL; /* regular expression to filter output in ALLOW mode*/
regex_t *output_filter_reg = NULL; /* regular expression data structure */
char *output_veto_filter = NULL; /* regular expression to filter output in DENY mode */
regex_t *output_veto_filter_reg = NULL; /* regular expression data structure */
char *stdout_file = 0;/* If != 0 : File to which stdout is redirected */
char *stderr_file = 0;/* If != 0 : File to which stderr is redirected */

int NFSdelay = 61;/* maximum possible delay in seconds between file fclose() on one NFS client and file becoming visible (due to directory attribute cache) on another NFS client */

int BNXVersion = -1;/* BNX version :
			-1 == Unknown Version
			0 == Version 0.1
			1 == Version 1.x */
int BNXMinorVersion = -1;/* Only used with BNXVersion == 1 to distinguish 1.2 from 1.3 */
int BNXnoRunData = 0;/* If 1 : Don't output RunData header lines (if present from BNX 1.x input) to any BNX output files (downgrades format to BNX 1.0) */
int BNXheaderFXY = -1;/* Molecule headers have "Column StartFOV StartX StartY EndFOV EndX EndY" as additional fields : set on input, all BNX input files must match :
			 -1 : Unknown (no BNX files read in)
			 0 : First BNX file did NOT have these extra columns
			 1 : First BNX file had these extra columns */

int BNXheaderChim = -1;/* Molecule headers have 6 extra columns "FeatureMinPos FeatureMinScore FeatureMaxPos FeatureMaxScore FeatureMinScaled FeatureMaxScaled" : set on input, all BNX input files must match :
			  -1 : Unknown (no BNX files read in)
			  0 : One of the previous BNX files did NOT have these extra columns
			  1 : All previous BNX files had these extra columns */

int BNXRunHeaderChim = -1;/* RunData header have 2 extra columns : set on input, all BNX input files must match :
			  -1 : Unknown (no BNX files read in)
			  0 : one of the previous BNX file did NOT have these extra columns
			  1 : All previous BNX files had these extra columns */


/* NOTE during merge BNXminlen,BNXminSNR are replaced by the larger of : 
   1. The largest input value in any input BNX being merged (Not used if there is only a single Run Data line with the same value)
   2. The -minlen -minSNR command line option (if specified) 
*/
double BNXminlen = 0.0;/* For BNX version 1.x (if != 0.0) : global minimum length (kb) used to generate current BNX input (independent of per-run values in Run Data strings) */
double BNXminSNR[MAXCOLOR] = {0.0};/* for BNX version 1.x (if != 0.0) : global minimum SNR threshold value applied to labels (independent of per-run values in Run Data strings) */

int MapType = -1;/* Type of map read in from -i input files :
		    -1 : unknown (use filename suffix to infer : .vfixx or .bnx == Molecule Map, .cmap = Assembly/Contig Map)
		     0 : Molecule Map
		     1 : Assembly/Contig Map */

int MapTrueSite = 0;/* 1 : input .bnx has truesite information (one extra line per color, with two values per map site) 
		       -1 : No truesite information in .bnx, but site mapped to original input is tracked
		     */
int MapSNR = 0; /* 1 if input .bnx has SNR information (one extra line per color after map sites lines) */
int MapIntensity = 0;/* 1 if input .bnx has Intensity information (one extra line per color after map sites lines) */
int MapStitched = 0;/* 1 if input .bnx has stitch information (one extra line per color after map sites lines) */
int MapStitchLoc = 0;/* 1 if input .bnx has stitch Location information (one extra line per color after map sites lines) */
int MapPSFWidth = 0;/* 1 if input .bnx has PSFWidth site information (one extra line per color after map sites lines) */
int MapImageLoc = 0;/* 1 if input .bnx has Image Location information (one extra line per color after map sites lines) */

int MapSNROrig,MapIntensityOrig,MapStitchedOrig,MapStitchLocOrig,MapPSFWidthOrig,MapImageLocOrig;/* original command line values (1 iff corresponding information is Required for all input bnx maps) */

int Tracksites = 0;/* Use TrueSite[] values to track original map index as map is transformed by -mres (or -rres) etc */
char *QueryOrigin = NULL;/* If Tracksite != 0 : original input map name (or names separated by commas) for query maps */
char *RefOrigin = NULL;/* If Tracksite != 0 : original input map name (or names separated by commas) for reference maps */
int RefTrueSite = 0;/* > 0 : RefOrigin and corresponding truesite information read in from one or more -ref input files 
		       NOTE : MapTrueSite is set of if truesite information is read in from one or more -i input files */

int QXfilter = 0;/* If 1 : filter out QX codes other than those specified on command line (see MapSNROrig etc) */
int QXerror= 0;/* If 1 : exit with error message if different BNX files have different sets of QX codes (otherwise just discard QX codes not present in all BNX files) */
int QXmismatch = 0;/* If != 0, some files have different QX codes than 1st file */

int CmapSNR = 0;/* If 1 : output complete SNR distribution for consensus CMAPs (in addition to Gmean and logSNRsd). Requires MapSNR==1 */
int CmapChimQuality = -1;/* 1 : Some Cmap input had ChimQuality data (0 means none, -1 means unknown : no input Cmap files read, 2 means SegDupL .. OutlierQuality,ChimNorm present, 3 means FragSd,ExpSd,FragCov,FragChiSq present) */
int CmapMask = -1;/* 1 : Some Cmap input had Mask data (0 means none, -1 means unknown : no input Cmap files read) */

char *BackgroundSNR = 0;/* If != 0 : This is the name of the file with the background SNR distribution and is to be used to update Pvalue computations when both maps have SNR distributions (must be consensus maps) */
int bgSNRcnt[MAXCOLOR];
double *bgSNR[MAXCOLOR];/* If BackgroundSNR != 0 : bgSNR[0..bgSNRcnt-1] is the background SNR distribution read in from the file named BackgroundSNR */

double bpp = 0.500;/* command line value */
double BppStepsiz = 0.0;/* bpp initial value scan  step size */
int BppSteps = 0;/* bpp initial value scan number of steps (in each direction from the original value) */
int BppMincnt = 10;/* bpp initial value scan is triggered when fewer than this many molecules align */
int NoBpp = 0;/* Do not apply bpp value to rescale -i input files (see -NoBpp) */
double PixelLen = 0.0;/* bpp value read from vfixx file (If != 0) and modified subsequently by command line bpp value and -M (from updated bpp estimates) */
double origPixelLen = 0.0;/* original PixelLen value before -bpp input */
double refPixelLen = 0.0;/* PixelLen value when reference map is read and de-res'd (After -bpp input BUT before -M iterations change PixelLen) */


int bnxerror = 0;/* If map from -i input .bnx file has sites beyond map ends, exit with an error message (instead of printing a WARNING and extending map ends as with .cmap) */
int bnxmergeerror = 0; /* If map from multiple -i input BNX 1.0 files has missing or duplicate Run Data lines exist with error message (otherwise just print warning message) */
int bnxBigid = 0; /* If 1 : Use 19 digit id numbers, if possible (if fields overflow, renumber all ids sequentially starting at 1)
		     If 0 : Avoid changing id numbers (if duplicates are present, renumber all ids sequentially starting at 1) */

double EndTrimCov = 7.9;/* trim back ends of draft consensus while coverage is below this level */
double EndTrimFrac = 0.0;/* If > 0.0 : Estimate length of unlabeled ends based on this percentile of distance from last aligned label to molecule end or next label */
double EndTrimCnt = 0.0;/* If > 0.0 : Mark end as non-extendible if at least this many molecules agreed (within sizing error) on length of unlabeled ends */
double EndTrimWt = 0.0;/* If > 0.0 : Use this instead if RefineMinMapwt to decide which molecules to count for EndTrimFrac and EndTrimCnt */

double EndTrimCov2 = 0.0;/* Further trim back ends of draft consensus until no Internal regions with coverage below this level remains, but not more than EndTrimLen2 kb beyond original ends
			    (typically EndTrimCov2 <= EndTrimCov) */
double EndTrimLen2 = 0.0;

double EndTrimNoRefine = 0.0; /* If > 0.0 : Don't refine consensus more than this distance beyond initial estimate of where coverage drops below EndTrimCov */

int EndTrimCQ = 0;/* see -EndTrimCQ */

double EndTrimRatio = 0.0; /* If > 0.0 : After trimming ends based on EndTrimCov (and before applying EndTrimFrac) keep trimming back the last label (at either end) if the its label coutn is <= EndTrimRatio times the label cnt of the next label */

double res[MAXCOLOR]; /* Query resolution in pixels (PixelLen) */
double resSD[MAXCOLOR];/* Query resolution SD in pixels (PixelLen) */

double FP[MAXCOLOR];
double FN[MAXCOLOR];
double R[MAXCOLOR];
double SD[MAXCOLOR+1];
double SF[MAXCOLOR+1];
double SR[MAXCOLOR+1];
double SE[MAXCOLOR+1];
double MA_mean;/* mean misalignment value (drift bias) between 2 colors */
double MA_SF,MA_SD;/* NOTE MA_SF and MA_SD are for misalignment between 2 colors */

double origFP[MAXCOLOR], origFN[MAXCOLOR], origSF[MAXCOLOR+1], origSD[MAXCOLOR+1], origSR[MAXCOLOR+1], origSE[MAXCOLOR+1],origres[MAXCOLOR],origresSD[MAXCOLOR];


/* NOTE SF[2] and SD[2] are for misalignment between 2 colors */

double MinFP = 0.01;/* see MINFP in constants.h : can be modified with -MinFP */
double MaxFP = 4.0;/* see MAXFP in constants.h : can be modified with -MaxFP */

double MinFN = 0.001;/* see MINFN in constants.h : can be modified with -MinFN */
double MaxFN = 0.40;/* see MAXFN in constants.h : can be modified with -MaxFN */

double MinSF = 0.001;/* see MINSF in constants.h : can be modified with -MinSF */
double MaxSF = 1.000;/* see MAXSF in constants.h : can be modified with -MaxSF */

double MinSD = -0.50;/* See MINSD in constants.h : can be modified with -MinSD */
double MaxSD = 0.50;/* See MAXSD in constants.h : can be modified with -MaxSD */
double MinSD_Mult = 0.99;/* see MINSD_MULT in constants.h : can be molified with -MinSD_Mult */

double MinSR = 0.00;/* see MINSR in constants.h : can be modified with -MinSR */
double MaxSR = 0.50;/* see MAXSR in constants.h : can be modified with -MaxSR */

double MinSE = 0.00;/* see MINSR in constants.h : can be modified with -MinSR */
double MaxSE = 0.00;/* see MAXSE in constants.h : can be modified with -MaxSE */

double MinRes = 0.01;/* minimum value of -res parameter */
double MaxRes = 9.99;/* maximum value of -res parameter */

double MinResSD = 0.01;/* minimum value of -resSD parameter */
double MaxResSD = 9.99;/* maximum value of -resSD parameter */

double maxresbias = 0.0;/* If != 0, intervals below this value are corrected for -ve bias caused by limited resolution. Uses piecewise Linear regression of bias vs log(x) with RESBINS bins */ 
int ResBins[MAXCOLOR] = {0,0}, origResBins[MAXCOLOR] = {0,0};/* number of piecewise linear segments to use */
int ResBinsFile = 0;/* If != 0 : ResBins and maxresbias were read with -readparameters */
double resbias[MAXCOLOR][RESBINS+1];/* resbias[c][i=0..ResBins-1] is the bias (in KB) of raw intervals x at  x = resbiasX[c][i] for color c */
double resbiasX[MAXCOLOR][RESBINS+1];/* resbiasX[c][ResBins] == maxresbias and resbiasX[c][0] == mres*0.5 and the other knots in the piecewise Linear model are selected to be seperated by similar number of sample values */
/* NOTE : if resbias parameters were read in using -readparameters, the correction is applied after -i map input. The best (additional) correction is also estimated using -M */

int RA_MIN_MEM = 1 /* (USE_MIC ? 3 : 1) */;/* value of MIN_MEM in refalign.cpp, set by -RAmem */
int RA_MIN_TIME = 0;/* free memory if memory use exceeds 90% of -maxmem OR every RA_MIN_TIME seconds if memory use exceeds RA_MAX_MEM times -maxmem */
int QP_MIN_MEM = (USE_MIC ? 3 : 1);/* value of MIN_MEM in qprobeval.cpp, set by -QPmem */
int MP_MIN_MEM = (USE_MIC ? 3 : 1);/* value of MIN_MEM in mprobeval.cpp, set by -MPmem */
int MP_MIN_TIME = 0;
double RA_MAX_MEM = 0.60;/* free memory every RA_MIN_TIME, but only if memory use exceeds this fraction of the threshold at which memory is freed at every oppertunity */

int RAscoremem = 1;/* If 1, deallocate score table memory before each refine() call */

int ResEstimate = 0;/* If != 0 : Estimate res and resSD */

int minSNRestimate = 0;/* If != 0 : Estimate minSNR[] in range SNRmin..SNRmax */
double SNRmin,SNRmax,SNRstep = -1.0, SNRminratio = 0.0;
int SNRmaxmaps = 0;/* If != 0 AND SNRstep > 0.0 : subset number of maps to no mroe than SNRmaxmaps which performing minSNR scan */

int outlierRate = 0;/* output outlier and Endoutlier rate as 2 additional columns in .err file */
int perLabelScore = 0;/* output score per label in .err file */
int LabelDensityType = 0;/* see -LabelDensityType */

int ScanCorrection = 0;/* If > 0 : Apply per Scan scaling correction to all -i input maps. Value is number of gradient updates per -M iteration (Only applied with BNX 1.0 if UniqueScans > 1) */
CscanCorrection *ScanScale[MAXCOLOR] = {0,0};/* Scaling factor applied to different scans : typically the weighted average scaling value will be close to 1.0 since this is applied in addition to global -bpp correction
					    Indexed as ScanScale[c][0..UniqueScans-1] : NOTE these parameters are NOT saved in .errbin or read in with -readparameters. Instead the rescaled .bnx is saved as <pr

efix>_rescaled.bnx
					    (without the -bpp scaling or bias correction) */
int ScanCorrectMinMaps = 1;/* Minimum number of aligned maps in a single cohort/scan to update the scaling factor */
double ScanCorrectMinSNR = -1.0;/* If >= 0.0 : -minSNR value applied just before _rescaled.bnx output */
double ScanCorrectMinLen = -1.0;/* If >= 0.0 : -minlen value applied just before _rescaled.bnx output */

int ScanFilter = 0;/* see -ScanFilter */
double ScanFilterRm = 0.0;
double ScanFilterRb = 99.0;
double ScanFilterRe = 99.0;

double MaxCoverageAN = 0.0;/* If > 0.0 : reduce coverage to no more than this amount before _rescaled.bnx output (relative to -ref size) */
double MaxCoverageLenAN = 0.0;/* If > BNXminlen : reduce coverage by first increasing BNXminlen up to MaxCoverageLen, then (if needed) randomly subsampling the input */

/* NOTE : currently only one set of scaling factors ScanScale[0][] is computed for 2 colors and 
   applied equally to both colors */

struct CRunDataSort *RunDataSort = 0;/* information about Runs from BNX 1.0 input */

char *parametersfile = 0;/* If != 0 : filename (.errbin) from which to read error parameters include bias parameters */

int floatcov = 1;/* output Coverage and Occurrence columms in .cmap as float rather than int */

int usecolor = 0;/* If > 0 : all input data must have >= 2 colors and only the specified color will be kept and the other color ignored (except for _q.cmap or _q.bnx output) */
int hashcolor = 0;/* If > 0 : If input data has >= 2 colors AND usecolor==0, then only the specified color is used by the hashtable (Currently required if colors >= 2 && usecolor==0) */
int startcolors = 1;/* original value of colors */
int refcolors = 0;/* number of colors in -ref inputs (may be 2, even if colors=1,  in which case 2-color .err files will be output if -usecolor was used) */

double mres = 0.001;
int mres_defined = 0; /* If > 0 : -mres was specified on commandline */
int mres_override = 1;/* If > 1 : -mres -mresSD on commandline overrides value from -readparameters */
double mresSD = 0.0;
double SDrange = 3.0; /* assume likelihood of gaussian more than this many SD from mean is negligible */
double HSDrange = 2.0;/* value of SDrange used in hash.cpp : usually between 1.0 and SDrange to improve performance of hashtable */
double HResMul = 0.0;/* As alternative to HSDrange * resSD use HRESmul * res (whichever is larger) : only used in hashtable */
double rres = 0.0;/* reference maps are de-resed at rres * 0.500 */
double cres = 0.0;/* If > rres : Check refined consensus maps : if 2 or 3 sites closer than (this many multiples of 0.5 kb) have fewer than cresN molecules then merge the site pair that is closest and repeat */
double ZERO_MINKB = 0.001;/* If != 0.0 : minimum distance between labels enforced during refinement (otherwise minimum distance enforced is rres * 0.500 kb) */

int cresN = 3;
double cresF = 0.10;
int cresFix = 0;// When merging 2 labels with -cres, check if deleting either label scores better than merging them
double cresMaxLPdrop = -1.0;/* If >= 0.0 : Don't merge sites with -cres if LP drops by more than cresMaxLPdrop (see -cresMaxLPdrop) */

double cresHapL = -1.0;// see -cresHapMaxLPdrop (1st arg)
double cresHapMaxLPdrop = 0.0;// see -cresHapMaxLPdrop (2nd arg)

int rresDefer = 0;// When Haplotyping defer applying -rres until after Haplotyping (this is always done for -cres) 
int refineDefer = 0;// When Haplotyping defer refinement and proceed directly to Haplotyping

double cresRefine = 1;/* If 1 : Refine sizes after each cres merging stage (otherwise only Refine sizes once at end of last cres merging stage) */

int AlignScore = 0;

/* defaults are for Reference Alignment */
double Kbias = 2.2;
double KLbias = 2.1;
double KFbias = 2.1;
double Gbias = 1.0;
double Beta = 1.0;

double KbiasRef = 2.2;
double KLbiasRef = 2.1;
double KFbiasRef = 2.1;

/* input file based XmapCount, XmapSites and XmapLength (valid if XmapCount > 0). Read in from a file using -XmapStatRead. Is output to a file (based on current bnx mapset) by using option -XmapStatWrite */
int XmapCount[MAXCOLOR];
long long XmapSites[MAXCOLOR];
double XmapLength[MAXCOLOR],XmapGtheta[MAXCOLOR];
char *XmapStatRead = 0,*XmapStatWrite=0;/* If != 0: input/output XmapCount,XmapSits and XmapLength from a file based on command line options -XmapStatRead and -XmapStatWrite */

double Ch = 0.0;/* chimerism probability */
double ChFdr = 0.0;/* chimerism False Discovery Rate */

double Poutlier = 0.0003;/* outlier probability per true positive site */
double PoutlierS = 0.01;/* outlier probabbility per true positive site for intervals with a stitch (if known) */
double PoutlierEnd = 0.0;/* outlier probability per true positive site for alignment ends */
int PoutlierEnd_init = 0;/* If PoutlierEnd was specified on the command line */

int EndOutlierPV = 0;/* If 1, correct pvalue for number of end outliers in alignment */

int isLinear = 0;/* if genome is linear */

double ScoreThreshold = 0.0;/* Score Threshold : only alignments with score above this threshold are output */
double LogPvThreshold = 6.0;/* Pvalue Threshold, expressed on the -log10 scale */
int AlignedSiteThreshold = 2;/* minimum aligned sites */
double AlignedLengthThreshold = 0.0;/* minimum aligned length */
int AlignedEndOutlierThreshold = 2;/* maximum number of End Outliers */
double AlignedOutlierThreshold = 1000.0;/* If < 999.0 : maximum internal outlier sizing error (alignments filtered out otherwise) */
int AlignedOutlierLabels = 1000;/* internal outliers with at least this many misaligned labels will cause alignments to be filtered out */

double AlignedFractionThreshold = 0.0;/* see -F option */
double AlignedLabelDensityThreshold = 0.0;/* If > 0.0 : see -D option */

int RefineEndOutlierThreshold = 2;/* maximum number of End Outliers allowed during refinement (otherwise weight of map is reduced to RefineEndOutlierWt except during final UpdateMap) */
double RefineEndOutlierWt = 0.001;
double RefineEndOutlierLen = 0.0;/* Only count End Outliers larger than this size (in both ref & qry) */
int RefineEndOutlierN = 0;/* Only count End Outliers with this many labels (in both ref & qry) */
int RefineEndOutlierType = 0;/* If 1 (UNLESS -refine 3 -extonly) weights are NOT restored when computing -Haplotype, -Chimquality or -contigsplit : useful for debugging -extsplit */

double RefineMinMapwt = 0.30;/* minimum map weight to output in refined .xmap (based on -BestRefWT -LRbias and -TB and Haplotyping Allelic balance) */
double UnrefineMinMapwt = 0.0; /* minimum map weight to output in regular (unrefined) .xmap (based on -BestRefWT only) */

int FilterXmap = 1;/* see -FilterXmap */
int FilterXmapWT = 0;/* see -FilterXmapWT */
int FilterXmapSplit = 0;/* see -FilterXmapSplit */

int RefineOverwrite = -1;/* Overwrite pre-existing refined contig */

double ScoreThreshold2 = -1000.0;/* Score Threshold for 2nd best alignment (see -SecondBest or -MultiMatch) */
double LogPvThreshold2 = 0.0;/* Pvalue Threshold for 2nd best alignment, expresssed on the -log10 scale (see -SecondBest or -MultiMatch) */
int AlignedSiteThreshold2 = 0;
double AlignedLengthThreshold2 = 0.0;
int AlignedEndOutlierThreshold2 = 2;
double AlignedOutlierThreshold2 = 1000.0;
int AlignedOutlierLabels2 = 1000;

int MinDupSites = 2;/* Minimum number of label in (non-inverted) Duplication : must be at least 1 */

int CutFlip = 0;/* If > 0 : Try to align each internal outlier of any primary alignment with this many site intervals in reverse orientation using the same query and reference region expanded
		   by at most CutFlipExpand labels. The Flipped alignment must satisfy at least the 2nd set of thresholds (LogPvThreshold2 etc) */

double CutFlipMaxErr = 0.10;/* maximum relative sizing error of internal outlier checked for inversions */
int CutFlipExpand = 0;/* Max number of labels by which to expand outlier intervals during inversion */
int CutFlipEndExpand = 0;/* If > 0 : Also apply CutFlip to EndOutliers and expand reference in End (non-aligned) direction by this many labels (typically CutFlipEndExpand > CutFlipExpand) */
int CutFlipDupExpand = 0;/* If > 0 : Also apply CutFlip to find Inverted Duplications and expand reference in either neighborhood by this man labels larger than outlier (insertion) query interval */
double CutFlipMinErr = 0.8;/* Minimum value of (X - Y)/(X+Y) for Insertions that are checked for Inverted Duplications */
int CutFlipEndFilter = -1;/* If > 0 : Suppress CutFlip Endoutlier inversions if there are already matchgroups that extend into the endoutlier region of either qry or ref by at least "CutFlipEndFlilter" Labels */

double InversionNormalLogPv = -1.0;/* If > 0 : suppress Inversion if normal matchgroups LogPv is worse than this value */
double InversionNormalLogPv5 = 15.0;/* enhanced threshold for non-inverted MG LogPv if inverted to non-inverted MG size ratio is >= 5.0 */

double ScoreThreshold3 = -1000.0;
double LogPvThreshold3 = 0.0;
double ScoreThreshold4 = -1000.0;
double LogPvThreshold4 = 0.0;

int CutFlipBC = 0;/* Apply local Bonferroni correction to Pvalue for CutFlip inversion matchgroups : this mainly affects confidence of end outlier CutFlip inversions based on 
		     difference between labels in endoutlier vs labels in inversion matchgroup */

int CutFlipSkip = 0;/* If > 0 : skip checking for CutFlip inversions in outliers with at least CutFlipSkip labels (typically > RefSplit), assuming these inversion matchgroups have already been found
		     */
double CutFlipMerge = 0.0;/* If > 0.0 : Merge outliers seperated by aligned region with score less then -log(CutFlipMerge), before applying CutFlip to the merged region (If <= 0, use Rsplit instead) */

double CutFlipMaxOverlap = 0.5;/* Maximum number of labels overlap between matchgroups for a CutFlip Inversion (as a fraction of aligned labels in inversion matchgroup) */
double CutFlipMaxMiss = 0.0;/* If > CutFlip : Maximum number of misaligned labels in each breakpoint gap for CutFlip Inversions (as a fraction of aligned labels in inversion matchgroup) */
double CutFlipOverlapFilter = 0.0;/* If > 0 : Do partial filtering of overlapped matchgroups based on this overlap ratio & CutFlipFilterConf, before doing 2-matchgroup inversion (CutFlip) calls */
double CutFlipFilterConf = 0.0;/* see CutFlipOverlapFilter */

double InvOverlapFilterQryRatio = 0.99;/* If > 0 : suppress filtering of opposite orientation overlapped matchgroups (with overlap exceeding OverlapFilterQryRatio),
					  unless overlap in query exceeds InvOverlapFilterQryRatio (in kb) AND InvOverlapFilterQryLabRatio (in Labels) AND confidence of larger matchgroup exceeds 
					  that of smaller matchgroup by InvOverlapFilterConfDelta. The saved smaller matchgroup can only be used to call inversions */
double InvOverlapFilterConfDelta = 10.0;
double InvOverlapFilterQryLabRatio = 0.0;

double OverlapFilterQryRatio = 0.7;/* If > 0 : Filter smaller matchgroup overlapped in query by larger matchgroup by OverlapFilterQryRatio (kb ratio) AND OverlapFilterQryLabRatio (Label ratio) 
				      provided larger matchgroup Confidence exceeds smaller matchgroup confidence by at least OverlapFilterConfDelta */
double OverlapFilterQryLabRatio = 0.0;
double OverlapFilterConfDelta = 0.0;

int OverlapFilterSameRef = 0; /* If > 0 : Don't apply OverlapFilter unless ref contigs for the overlapped MGs are the same */


double MergeMG = 0.0;/* If > 0 : After running -OverlapFilter, try to merge neighboring matchgroups if they are the closest matchgroups with no overlap and no other matchgroups overlap the gap region */

double SmallDupMaxOverlap = 0.0/*0.20*/;/* Max Overlap (fraction of aligned labels) of small matchgroup with non-outlier region of query in larger matchgroup, for small Duplication calls */
double SmallDupMinOverlap = 0.0/*0.80*/;/* Min Overlap (fraction of aligned labels) of small matchgroup with non-outlier region of query in larger matchgroup, for small Duplication calls */
double SmallDupMaxMiss = 0.20;/* Max fraction of aligned labels of small matchgroup that do not match an aligned label in the reference of the larger matchgroup, for small Duplication calls */
double SmallDupConfBC = 4.0;/* Minimum confidence of small matchgroup after Bonferroni correction by duplicate region gap in labels */

double LogPvThresholdTE = 6.0;/* Pvalue Threshold corresponding to -TE option for extensions */
int AlignedSiteThresholdE = 2;/* see -AE */
double AlignedLengthThresholdE = 0.0;/* see -LE */

double LogPvThresholdTS = 6.0;/* Log10 Pvalue Threshold correspondign to -TS option for E&S contigs */
int AlignedSiteThresholdS = 2;/* see -AS */
double AlignedLengthThresholdS = 0.0;/* see -LS */

int TEinit = 0;/* If -TE option was specified */
int AEinit = 0;/* If -AE option was specified */
int LEinit = 0;/* If -LE option was specified */

int TSinit = 0;/* If -TS option was specified */
int ASinit = 0;/* If -AS option was specified */
int LSinit = 0;/* If -LS option was specified */

int TBcnt = 0;
double LogPvThresholdTB[16];/* LogPvThresholdTB[0..TBcnt-1] : Pvalue Threshold(s) corresponding to -TB option for extension maps after preliminary refinement */
double TBmult = 0.0;/* If != 0 : Replace LogPvThreadholdTB[0..TBcnt-1] by a sequence starting at LogPvThresholdTB[0] and decreasing in multiple of TBmult until LogPvThresholdTB[TBcnt-1] is reached.
			If the last value is not equal to LogPvThresholdTB[TBcnt-1], that is added to the sequence */
int TBinit = 0;/* If -TB option was specified */

double MinLen = 0.0;
double MaxLen = 0.0;/* If > 0 : Maximum -i input Map length in kb */
int PostAdjust = 1;/* Apply Minlen/Maxlen after doing -bpp adjustment */

int MinSites = 1;
int MaxSites = 0;/* If > 0 : delete maps with more than this many sites */
int MaxSites2 = 0;/* If > 0 : delete maps with more than this many sites just before output of _rescaled.bnx in refalign.cpp */
double MaxSiteDensity = 0.0;/* If > 0 : delete maps with more than this many sites / 100kb */
double MinSiteDensity = 0.0;/* If > 0 : delete maps with less than this many sites / 100kb */

double MaxEndSiteDensity = 0.0;/* If > 0 : delete maps with more than this many sites / 100kb in either Map End of at least MaxEndSiteDensityLength kb */
double MaxEndSiteDensityLength = 50.0;

double MaxEnd = -1.0;/* If >= 0 : Trim (reduce) unlabeled molecule ends to no more than MaxEnd kb */
double MaxInterval = -1.0;/* If > 0 : Split molecules with internal intervals larger than MaxInterval kb. */
int refMaxInterval = 0;/* If 1 : Also apply MaxInterval to reference CMAP input, otherwise only to query BNX input (and query CMAP input if there is NO reference)
			  If applied to ref OR query CMAP input, the maps may be split into multiple maps. With BNX input only the largest piece is retained, so number of maps remain unchanged */

double MaxIntensity = 0.0;/* If > 0 : Maximum AvgIntensity value of molecule backbone */
double minSNR[MAXCOLOR] = {0.0,0.0};/* minSNR[c] is the minimum SNR filter for color c (c = 0 or 1) */
double maxPSFWidth[MAXCOLOR] = {0.0,0.0};/* maxPSFWidth[c] is the maximum PSFWidth for color c (c = 0 or 1) */

int MaxThreads = -1;/* limit number of threads (-1 means : not specified on command line) */
int TotalThreads = 0;/* If TotalThreads > 0 and TotalSlots = 0: Total Threads of all jobs on machine : scale MaxMem and HMaxMem by MaxThreads/TotalThreads (1 means use actual total threads on machine) */
int TotalSlots = 0;/* If TotalThreads > 0 && TotalSlots > 0 : scale MaxMem and HmaxMem by MaxThreads/(TotalThreads*OMP_NUM_THREADS/TotalSlots) */
int RefineThreads = -1;/* limit number of threads during refinement (-1 means : not specified on command line and defaults to same as MaxThreads) */
int BoostThreads = 0;/* If > 0 : Increase RefineThreads to this number after elapsed time exceeds BoostThreadsTime. Do NOT re-adjust maxmem,hashmaxmem due to TotalThreads */
int BoostThreadsTime = 0;

double MaxMem = 128.0;/* If > 0 : Maximum memory in Gigabytes */
double HMaxMem = 0.0;/* If > 0 : Maximum memory in Gigabytes for hashtable (defaults to MaxMem) */
double MaxVirtMem = -1.0;/* If > 0 : Maximum Virtual memory in Gigabytes (0 means unlimited, -1 means (TotalThreads ? 3 : 1) x MaxMem)  */
double MaxMemIncrease = 0.0;/* If > 0 : set MaxMem to at MemSiz - MaxMemIncrease, where MemSiz is system memory in Gigabytes */


int colors = 1;/* number of colors in input files (typically based on -i inputs after processing -usecolor) */
int colors_set = 0;/* > 0 IFF number of colors has been specified by -colors OR first input file (-i or -readparameters) has been successfully read in confirming the value of -colors parameter.
		      All subsequent input files must agree with this value */

char *Nickase[MAXCOLOR];/* color names based on first input query (or ref) file read in : only used in input_vfixx.cpp */

int DELTA_Y = 6; /* maximum span of aligned sites in reference */
int DELTA_Y2 = 0;/* If > 0 : smaller value of DELTA_Y to use in all but the last -M iterations of refalign (to achieve greater speed) */
int DELTA_X = 4; /* maximum span of aligned sites in Nanomap */
int DELTA_X2 = 0;/* If > 0 : smaller value of DELTA_X to use in all but the last -M iterations of refalign (to achieve greater speed) */
int DELTA = 4;  /* maximum span of aligned sites for pairwise alignment */
int KMAX = 6; /* maximum number of consecutive references sites to merge due to resolution limit is KMAX+1 */

int deltaXRef = 0;/* If > 0 : modify DELTA_X during refinement to this value */
int deltaYRef = 0;/* If > 0 : modify DELTA_Y during refinement to this value */
int deltaExtXRef = 0;/* If > DELTA_X, allow X (Query) intervals to be this large as long as Y (Ref) intervals are limited to DELTA_Y */
int deltaExtYRef = 0;/* If > DELTA_Y, allow Y (Ref) intervals to be this large as long as X (Query) intervals are limited to DELTA_X */
int deltaEndRef = 0;/* If > 0 : allow non-outlier ends with up to deltaExtXRef + U - 1 qry labels (or deltaExtYRef + U - 1 ref labels), instead of just deltaXRef,deltaYRef labels */

int outlierExtend = 0;/* If > 0 : increase DELTA, DELTA_X and DELTA_Y as needed to try to enlarge internal outliers up to this many sites beyond the current outlier boundaries */
int outlierExtendLimit = 0;/* If > 0 : limit increases in DELTA, DELTA_X or DELTA_Y to never exceed this value */
int outlierBC = 0;/* If > 0 : Apply BonFerroni correction to outlier penalty */
RFLOAT outlierMax = 1e+30;/* Do NOT allow outliers with delta |X-Y| that exceeds this value in kb (only applies if maptype==0) */
double outlierLambda = 1e+30;/* If OUTLIER_LTYPE==0 : outlier likelihood is proportional to exp(-(deltaX+deltaY)/outlierLambda))
			       If OUTLIER_LTYPE==1 : outlier likelihood is proportional to exp(-|deltaX-deltaY|/outlierLambda)) */
RFLOAT OutlierType0 = 0.0;/* see OUTLIER_TYPE in RGentigRefScore.h */
RFLOAT OutlierType1 = 0.0;/* see OUTLIER_TYPE1 in RGentigRefScore.h */
int OutlierType = 1;/* see OUTLIER_TYPE in RGentigScore.h */

int outlierNorm0 = 0;/* see -outlierNorm0 (or LAMBDA_FIX0 in RGentigRefScore.h) */
int outlierNorm1 = 0;/* see -outlierNorm1 (or LAMBDA_FIX1 in RGentigRefScore.h) */
int outlierNormRef = 1;/* see -outlierNormRef (or LAMBDA_FIX in RGentigScore.h) */

double thetaScale = 1.0;/* see -thetaScale */
double thetaScaleRef = 1.0;/* see -thetaScaleRef */

double RefineStep1 = 0.10;/* Initial iteration minimum sample sizing increment before UpdateMap() is called (except during RefineA). Smaller values produce more accurate consensus maps but longer run times */
double RefineStep2 = 0.20;/* Initial iteration minimum sample sizing increment if RefineA OR after UpdateMap() is called */
double RefineStep3 = 0.10;/* Initial & Subsequent iteration sizing increment (when convergence failure is detected) Before UpdateMap() is called (except during refineA) */
double RefineStep4 = 0.20;/* Initial & Subsequent iteration sizing increment (when convergence failure is detected) if refineA OR after UpdateMap() is called */
double RefineStep5 = 0.20;/* Absolute minimum sizing increment (in kb) Before UpdateMap() is called (except during refineA) */
double RefineStep6 = 0.20;/* Absolute minimum sizing increment (in kb) if refineA OR after UpdateMap() is called */
int RefineStepCnt = 3;/* max number of steps to use with RefineStep3 or RefineStep4 (the number of steps is NOT limited with RefineStep1 or RefineStep2) */
double RefineStep7 = 0.10;/* minimum relative sizing change in a neighboring interval that may trigger RefineStepCnt x 2 extra steps */
int RefineFix = 0;/* Turns on adjustments and bug fixes in Refinment introduced between revision 4535 and 4607 */

double RefineMinRange1 = 0.01;/* see -RefineSizeRange */
double RefineMaxRange1 = 2.00;
double RefineMinRange2 = 0.70;
double RefineMaxRange2 = 1.40;

int SiteMerge = 1;/* see -SiteMerge */

int IndelScore = 0;/* If 1 : modify score parameters -outlierLambda, -outlierType, -outlier, -outlierEnd, -LRbias during Final Haplotype Indel Scoring (also used to delete them, if score is -ve) */

double outlierLambdaRef = 0.0;/* If > 0 : Replaces outlierLambda during refinement */
int outlierLambdaSwitch = 0;/* force outlierLambdaRef to outlierLambdaLabel during refinement when adding/deleting labels (or Haplotype SNPs)
                               force outlierLambdaRef to outlierLambdaSize during refinement when adjusting sizing (or Haplotype Indels) */
double outlierLambdaLabel = 1e+30;/* 0.0 means initial value of outlierLambdaRef */
double outlierLambdaSize = 0.0;/* 0.0 means initial value of outlierLambdaRef */
double IndelDelOutlierLambda = 0.0;/* 0.0 means no change during Indel Deletion */

int OutlierTypeSwitch = 0;/* force OutlierType to OutliertypeLabel during refinement when adding/deleting labels (or Haplotype SNPs) 
			     force OutlierType to OutliertypeSize during refinement when adjusting sizing (or Haplotype Indels) */
int OutlierTypeLabel = 1;/* -1 means original value of OutlierType */
int OutlierTypeSize = 0;/* -1 means original value of OutlierType */
int IndelDelOutlierType = -1;/* -1 means no change during Indel Deletion */

int PoutlierSwitch = 0;/* force Poutlier to PoutlierLabel during refinement when adding/deleting labels (or Haplotype SNPs)
			 force Poutlier to PoutlierSize during refinement when adjusting sizing (or Haplotype Indels) */
double PoutlierLabel = -1.0; 
double PoutlierSize = -1.0;
double IndelDelPoutlier = -1.0;/* < 0.0 means no change during Indel Deletion */

int PoutlierEndSwitch = 0;/* force PoutlierEnd to PoutlierEndLabel during refinement when adding/deleting labels (or Haplotype SNPs)
			     force PoutlierEnd to PoutlierEndSize during refinement when adjusting sizing (or Haplotype Indels) */
double PoutlierEndLabel = -1.0;
double PoutlierEndSize = -1.0;
double IndelDelPoutlierEnd = -1.0;/* < 0.0 means no change during Indel Deletion */

double LRbias = 0.0;/* LR bias added to all LR values during refinement (but not before adding/deleting sites in the extension region, if extend >= 2) */
int LRbiasSwitch = 0;/* force LRbias to LRbiasLabel during refinement when adding/deleting labels (or Haplotype SNPs)
			force LRbias to LRbiasSize during refinement when adjusting sizing (or Haplotype Indels) */
double LRbiasLabel = 0.0;/* 0.0 means use regular value of LRbias */
double LRbiasSize = 0.0;/* 0.0 means use regular value of LRbias */
double IndelDelLRbias = 0.0;/* 0.0 means no change during Indel Deletion */

int MinSFSwitch = 0;/* force sf to at least MinSFLabel during refinement when adding/deleting labels (or Haplotype SNPs)
		       force sf to at least MinSFSize during refinment when adjusting sizing (or Haplotype Indels) */
double MinSFLabel = 0.0;
double MinSFSize = 0.0;
int MinSFLabelIter = 0;
int MinSFSizeIter = 0;

double MinSFIndelDel = 0.0;/* minimum sf value used during Indel Deletion */

int MinSRSwitch = 0;/* force sr to at least MinSRLabel during refinement when adding/deleting labels (or Haplotype SNPs)
		       force sr to at least MinSRSize during refinment when adjusting sizing (or Haplotype Indels) */
double MinSRLabel = 0.0;
double MinSRSize = 0.0;
int MinSRLabelIter = 0;
int MinSRSizeIter = 0;

double MinSRIndelDel = 0.0;/* minimum sr value used during Indel Deletion */

int ViterbiSwitch = 0;/* see -ViterbiSwitch option */
double ViterbiLabel = 0.0;
double ViterbiSize = 0.0;
int ViterbiLabelIter = 0;
int ViterbiSizeIter = 0;

int RangeSwitch = 8;/* Larger value of RefineRange to use when updating alignments during refinement in UpdateMap() */


double RefineEps1 = 1e-8;/* log-likelihood precision required for endoutlier terms (UL,UR) during refinement (mprobeval) */
double RefineEps2 = 0.0;/* log-likelihood precision required for total alignment of one molecule (LR) during refinement (mprobeval) */
double RefineEps3 = 0.0;/* log-likelihood precision required for sub alignment terms (AL,AR) during refinement (mprobeval) */

int extend = 0;/* allow nanomaps to extend beyond reference (2 to also refine extended reference map) */
int extendonly = 0;/* refine extension region only */
double MaxExtend = 0.0;/* If > 0 : Limit extension refinement to MaxExtend (in kb) at each end (to speed up refinement) */
int UseEndoutliers = 1/* WAS 0 */;/* If > 0 : Use the endoutlier region of molecules to generate hypothesis sites in consensus during refinement */
double EndLen = 100.0;/* If extendonly > 0 : The part of the original contig ends (in kb) that will be refined along with the extension region 
		         NOTE : If EndLen <= 0.0, then extend = 0 */
double EndLen2 = 0.0;/* 2nd arg of -extonly */
double MinOverlapRatio = 0.0;/* If extend > 0 && Refine > 0 : Ignore extending nanomaps with less than the specified ratio of their length overlapping the original reference contig */
                              /* In Assembler : minimum pairwise overlap as a fraction of the smaller map */

double extendSplit = 0.0;/* If > 0 : Clone part of contig up to a target site where at least max(extendSplit, extendFrac * coverage) endoutliers in the direction away from the cloned part are present.
			    Don't count endoutlier shorter than extendSplitLen (kb) OR whose alignment ends more than extendSplitFN sites (AND extendSplitFNlen kb) before target site or 
			    more than extendSplitFP sites after target site (suspected breakpoint).
			    If multiple nearby labels satisfy the rule, use the one with the largest number of endoutliers.
			 */
double extendFrac = 0.15;
int extendSplitFN = 6;
double extendSplitFNlen = 60.0;
int extendSplitFP = 2;
double extendSplitLen = 30.0;
int extendSplitSpacing = 3;
double extendSplitSpacingX = 2.0;
int extendSplitEndSpacing = 3;
double extendSplitEndSpacingX = 2.0;

double extendWT = 0.0;
double extendLen = 20.0;/* If extendWT > 0 : Alignments that extend extendLen beyond end of contig (not counting unlabeled mol end) have their align->mapWT increased to at least extendWT (see -BestRefWT) */
int extendFN = 6;
double extendFNlen = 60.0;
double extendMinWT = 0.0;
double extendMaxOutlierKB = 1000.0;
double extendMaxOutlierKBsplit = 1000.0;
int extendMaxOutlierLabels = 1000;
int extendMaxOutlierLabelsSplit = 1000;
double extendMaxWT = 1.0;

double extendWTdup = 1;/* If 0, only raise the weight of one extending alignment per molecule */
double extendWTend = 1;/* If 0, -extendWT only applies to split contig extending alignments (not regular contig ends) */

double splitWT = 0.01;/* If splitWT > 0.0 : Increase weight of alignments (that exceed -extendMinWT) to at least splitWT when computing coverage and endoutlier counts to determine splits.
			Also if extendWT > 0 : Apply extendWT (increase align->mapWT to extendWT) for endoutlier maps moved to splits as extensions */

double extsplitMinWT = 0.01;/* skip alignments with WT below this value when computing -extsplit coverage */
double refineMinWT = 0.0;/* skip alignments with WT below this value during refinement */
double refineWT = 0.0;/* If refineMinWT > 0 AND refineWT > 0, alignments with WT above refineMinWT will have their WT increased to at least refineWT */
int refineWT_N = 999;/* see 2nd value of -refineWT */
double refineWT_NL = 0.0;/* see 3rd value of -refineWT */
double refineWT_MaxOutlierKB = 999.0;/* see 4th value of -refineWT */
int refineWT_MaxOutlierLabels = 1000;/* see 4th value of -refineWT */

double splitFiltMinWT = 0.0;/* see -splitFilt */
int splitFiltFN = -1;
double splitFiltFNlen = -1.0;
int splitFiltFP = -1;

int keepsplitmaps = 0;/* If > 0 : keep maps with endoutliers that are used in a contig split in the original contig as well, but at reduced weight based on -BestRefWT */
int extTrim = 0;/* see -extTrim */


int ExtsplitSpacingFix = 1;/* see EXTSPLIT_SPACING_FIX in refalign.cpp */
int ExtsplitVerifyPeak = 1;/* see EXTSPLIT_VERIFYPEAK in refalign.cpp */

int extsplitNE = 1;/* If 1 : include non-endoutliers that align only to the split contig region in the split contig
		      If 0 : Don't include non-endoutliers in split contig, and limit -extonly and -maxExtend value for each split contig (see -extsplitNE) */
double extsplitNEL = 0.0;/* L value of -extsplitNE */
double extsplitNEE = 0.0;/* E value of -extsplitNE */

double extsplitOutlier = 0.0;/* If > 0 : Also include outliers of least this size difference if either boundary matches a split point (along with endoutliers) */
int extsplitOutlierLabels = 1000;
double extsplitOutlierLS = 0.0;
int extsplitOutlierAS = 0;
double extsplitOutlierKB[MAX_KB_BINS+2];/* extsplitOutlierKB[1..KBcnt] correspond to L1 .. Ln in the -extsplitOutlier option, extsplitOutlierKB[0] == extsplitOutlier, extsplitOutlierKB[KBcnt+1] = 9999 */
int KBcnt = 0;

int KCcnt = 0;
double extsplitOutlierC[MAX_KB_BINS];/* extsplitOutlierC[0..KCcnt-1] are the args of -extsplitOutlierC :
				        extsplitOutlierC[0] applies to insertion Bins B = 0 .. KCcnt-2 and extsplitOutlierC[T = 1..KCcnt-1] applies to deletion bin -T */


double CloneLimitLen = 0.0;/* If > 0.0 : Limit length of -extsplit cloned part of contig to no more than this many kb or CloneLimitN labels (whichever is greater) */
int CloneLimitN = 0;
int CloneTrim = 1;/* If > 0 : Allow -extsplit split contig's non-extension end to be trimmed by -EndTrim (see -CloneTrim) */

int HasLoc = 1;/* If startloc & endloc should be read in from Version 0.1 BNX input (if present) */

int InDel = 0;/* If > 0 : output alignment outliers in .indel file along with .xmap file (if -ref was used) */
int InDelEnds = 0; /* If > 0 (AND Indel > 0) : Also output end outliers in .indel file */
int InDelNoOverlap = 0;/* If > 0 (AND Indel > 0) : If the same query region aligns to two different reference region, do not output any indels from the overlapped region of the query */
double InDelNoOverlapT = 0.0;/* minimum confidence for matchgroups to suppress indels with IndelNoOverlap > 0 */
double InDelMinConf = 0.0;/* If > 0 (AND Indel > 0) : Indels must have minimum Left or Right alignment confidence of InDelMinConf, otherwise they are ignored. Useful to filter out indels from CutFlipPaste matchgroups */

double IndelChiSq = 0.0;/* During -indel calling ignore outliers where smaller map interval has ChiSq value below this value */
double IndelFrac = 100.0;/* During -inel calling ignore outliers where smaller map interval has OutlierFrac value above this value */

double IndelMaxSize = 35.0;/* Only ignore deletions with size (delta) below this many kb */
double IndelMaxCov = 40.0;/* Only ignore outliers with MolCov value below ~MaxCov */
double IndelMinSf = 0.450 * 1000.0;/* Only ignore outliers with MolSd above ~MinSf + ~MinSr * Interval */
double IndelMinSr = 0.02 * 1000.0;
double IndelMinInterval = 2.0;/* Only ignore outliers within intervals above this many kb */
double IndelSdRatio = 1.2;
double IndelMinCov = 3.0;
double IndelMinCovScale = 1.5;
int IndelWin = 0;
double IndelFragile = 999.0;
double IndelSegDup = 999.0;
double IndelChimNorm = 0.0;
double IndelMinCovRatio = 0.0;


double BreakChiSq = 0.0;/* During -indel calling ignore outliers where smaller map interval has ChiSq value below this value */
double BreakFrac = 100.0;/* During -indel calling ignore outliers where smaller map interval has OutlierFrac value above this value */

double BreakMaxSize = 35.0;/* Only ignore outliers in Qry Intervls below this many kb */
double BreakMaxCov = 40.0;/* Only ignore outliers with MolCov value below ~MaxCov */
double BreakMinSf = 0.450 * 1000.0;/* Only ignore outliers with MolSd above ~MinSf + ~MinSr * Interval */
double BreakMinSr = 0.02 * 1000.0;
double BreakMinInterval = 2.0;/* Only ignore outliers within Qry intervals above this many kb */
double BreakSdRatio = 1.2;
double BreakMinCov = 3.0;
double BreakMinCovScale = 1.5;
double BreakMinCovRatio = 1.2;
int NoBreak = 0;
double BreakFragile = 100.0;
double BreakSegDup = 100.0;
double BreakChimNorm = 0.0;

/* Parameters only used by RefAligner */

int NoSplit = 2;
int NoSplitSet = 0;/* 1 IFF -nosplit was specified on command line */

int PairSplit = 0;
double Psplit = 0.0;
int PairsplitExtend = 0;/* If PairSplit : When breaking molecules at internal outliers extend both pieces to include the outlier region (including any other outliers nearby) */
int PairsplitMerge = 0;/* If PairSplit : At the end try to merge overlapped or adjacent alignments */
int PairsplitMinOutlier = 0;/* If PairSplit  : Do NOT break alignments at internal outliers with less than this many site intervals on both maps (intended to work with -indel) */
int PairsplitPV = 0;/* If PairSplit : Prioritize alignments by LogPV instead of Likelihood score */
int PairsplitRef = 2;/* If PairSplit : treat first -i input map set as in-silico reference :
			>=1 : don't apply -bpp to reference
			>=2 : assume zero sizing errors in reference */
int PairsplitSeparateQueries = 1;/* If PairSplit : treat multiple queries as if each query were processed in a seperate command with the same set of reference maps */

int PairwiseOutlierFix = 1;/* If > 0 : fixes several issues with the handling of outliers and endoutliers during pairwise alignment */


int RefIndex = 0;/* Reference Maps are from vfixx_filename[0..RefIndex-1] : used for .xmap output */
int QueryIndex = 0;/* Query Maps are from vfixx_filename[QueryIndex..num_files-1] : used for .xmap output */

int CmapMerge = 0;
int CmapMergeBnx = 0;/* output -merge maps as .bnx */
int CmapSort = 0; /* Order in which maps are output when using -merge :
		     0 : Output maps in any order (typically the same order as they were input) 
		     1 : Sort maps in increasing order map size
		     2 : Sort maps in decreasing order map size
		     3 : Sort maps in increasing order of number of sites
		     4 : Sort maps in decreasing order of number of sites
		     5 : Sort maps in increasing order of map id
		     6 : Randomize order of maps (helpful to avoid memory bottlenecks)
		     7 : Sort maps in increasing order of GlobalScanNumber (ie in order of RunIndex,ScanNumber) see -sort-runindexinc
		     8 : See -sort-runindexShuffle : uses CmapSortMinLen as <MinLen> value
		  */
                
double CmapSortMinLen = 0.0;

int CmapSwap = 0;/* If >= 1 : swap the two colors in the input maps */

int RandSeed = 1;
int CmapStart = -1;/* If >= 1, the subset of maps to output (based on map order after sorting) is CmapStart .. CmapEnd */
double CmapEnd = -1.0;/* see -subset : If between 0 and 1 that fraction of the total maps is output */
double CmapFrac = 0.0;/* see -subset : output at least this fraction of total maps  */
int perCohort = 0;/* see -subset : output at least this many maps per cohort bin (same UniqueScanId and RunIndex/ScanNumber combination) */

double CmapGB = 0.0;/* If > 0.0, the subset of maps to output (based on map order after sorting) is 1 .. N, where N is the value such that the molecules 1 .. N that exceed CmapMinLen kb, add up to CmapGB * 1000000 kb */
double CmapMinLen = 0.0;
int CmapGBrnd = 1;/* If 1 : round up N to include all molecules from the same cohort (RunIndex) */

int FinalSort = 0;/* If >= 1, resort both -i and -ref maps (after -i maps have been processed with -sort* or -randomize and -subset) */

int CmapBin= 0, CmapNumBins = 0;/* If >= 1, the subset of map to output (based on partition the map into CmapNumBins roughtly equal sized sections and outputing section CmapBin */
int CmapSplit = 0;/* If >=1, split input maps into single map files */

double RangeLeft = 0.0, RangeRight = 0.0;
double InversionLeft = 0.0, InversionRight = 0.0;


double PairMerge = 0;/* output pairwise merged maps in CMAP format when alignment satisfies -S -T -A thresholds and has minimum overlap in kb >= PairMerge.
			Each input map must be in a seperate input file and have distinct ID values
			Unmerged maps are written to PREFIX_contig<N>.cmap, where <N> is the ID of the map.
			Merged map pairs (starting with highest Pvalue) are written to PREFIX_contig<M>.cmap, where <M> is the smaller
			of the two map IDs (the higher map ID will be missing from the output files). The merged map transitions between the
			original maps at the midpoint of the alignment. Only applies to pairwise alignment */
double PairMergeOverlapped = -1.0;/* If >= 0.0 : Reduced minimum overlap in kb when one map completely overlaps another map */
double PairMergeOverlappedMargin = 0.001;/* see 2nd arg of -pairmergeOverlapped */
double PairMergeOverlappedSmargin = 1.0;/* see 3rd arg of -pairmergeOverlapped */
double PairMergeOverlappedNmargin = 1;/* see 4th arg of -pairmergeOverlapped */
double PairMergeOverlappedEmargin = 0.1;/* see 5th arg of -pairmergeOverlapped */

double PairMergeOutlierChiSq = 0.0;/* During -pairmerge ignore outliers where smaller map interval has ChiSq value below this value */
double PairMergeOutlierFrac = 100.0;/* During -pairmerge ignore outliers where smaller map interval has OutlierFrac value above this value */

double PairMergeOutlierMaxSize = 35.0;/* Only ignore outliers with size (delta) below this many kb */
double PairMergeOutlierMaxCov = 40.0;/* Only ignore outliers with MolCov value below ~MaxCov */
double PairMergeOutlierMinSf = 0.450 * 1000.0;/* Only ignore outliers with MolSd above ~MinSf + ~MinSr * Interval */
double PairMergeOutlierMinSr = 0.02 * 1000.0;
double PairMergeOutlierMinInterval = 2.0;/* Only ignore outliers within intervals above this many kb */
double PairMergeOutlierSdRatio = 1.2;
double PairMergeOutlierMinCov = 3.0;
double PairMergeOutlierMinCovScale = 1.5;

double PairMergeChimNorm = 0.0;
double PairMergeOutlierMinCovRatio = 0.0;

double PairMergeMaxEnd = 0.20;/* Apply PairMerge only if each unaligned end (due to -endoutlier) does not exceed "overlap region" * PairMergeMaxEnd */
double PairMergeMaxOutlier = -1.0; /* If >= 0 : Apply PairMerge only if the largest internal outlier does not exceed PairMergeMaxOutlier (in kb) on either map */
double PairMergeMaxEndKB = 0.0;/* Apply PairMerge only if each unaligned end (due to -endoutlier) does not exceed PairMergeMaxEndKB (in kb) OR does not exceed PairMergeMaxEndN labels (at either end)  */
double PairMergeMaxEndExpandKB = 0.0;/* Allow outliers near end to be removed by merging with endoutlier, as long as endoutlier does not expand to more than this value (defaults to PairMergeMaxEndKB */
int PairMergeMaxEndN = 0;
int PairMergeMaxOutlierM = 999;

int PairMergeBigMapsFirst = 0;

int PairMergeMaxOutlierBalanced = 1;/* If outlier contains this many or more misaligned labels, use inversion size instead of Indel size (if larger) as outlier size */

int PairMergeExcludeOutliers = 0;/* If 1 : Exclude internal outliers in the computation of the alignment overlap size (in kb) : 
				    Overlap size is used by various suboptions of PairMerge */

int PairMergePreferNGS = 0; /* When merging preserve as much of the NGS map (StdDev=0) as possible */
int PairMergeXmap = 0;/* output .xmap files during pairmerge */
int PairMergeRepeat = 0;/* If >= 1 : Repeat PairMerge operation until no more maps merge (avoids file IO between each operation) */

double TrimmedNoMerge = 0.0;/* If > 0 : Don't merge any pair of maps if that eliminates a map end marked as trimmed and non-extendable (if Mask[0][1] or Mask[0][N+1]), unless 
			       one map is completely overlapped by the other (except for an extension of no more than TrimmedNoMergeMaxEnd) */
double TrimmedNoMergeMaxEnd = 0.0;
size_t TrimmedNoMergeBoth = 0;/* If != 0 : require both overlapped ends to be marked with bitwise intersecting values (that also intersect with TrimmmedNoMergeBoth), before preventing merges */
size_t TrimmedNoMergeSegDupMask = 0; /* If !=0 : Disallow merges if the overlap region Labels of either map have Mask bits that intersect with this value, unless the overlap extends beyond the Masked region
				                 by at least TrimmedNoMergeSegDupFlank labels in the non-extending direction */
int TrimmedNoMergeSegDupFlank = 5;

int TrimmedNoMergeStrict = 0;// see -TrimmedNoMergeStrict

double Truncate = 0.0;/* If > 0.0 : For any map pair with one endoutlier and at least Truncate kb and TruncateN labels aligned with no outliers larger than TruncateMaxOutlier,
			 truncate the shorter map from the end without the endoutlier to no more than Trucate kb or TruncateN labels (whichever is greater) and mark the trimmed 
			 end non-extendable using Mask[0][1] or Mask[0][N+1] */
int TruncateN = 0;
double TruncateMaxOutlier = 0.0;
int TruncateFix = 1;/* see 4th arg to -Truncate */
int TruncateMaxOutlierM = 3;/* M value (max misaligned labels for non-outlier interval) and 5th arg of -Truncate */

double MaxPalindrome = 0.0;/* If > 0.0 : -MaxPalindrome <L> value */
int MaxPalindromeN = 0; /* -MaxPalindrome <N> value */
double MaxPalindromeSL = 0.0;/* -MaxPalindrome <SL> value */
size_t MaxPalindromeMask = END_NOEXT;/* -MaxPalindrome D value */

double MaxInvPalindrome = 0.0;/* If > 0.0 : -MaxInvPalindrome <L> value */
int MaxInvPalindromeN = 0;/* -MaxInvPalindrome <N> value */
double MaxInvPalindromeSL = 0.0;/* -MaxInvPalindrome <SL> value */
size_t MaxInvPalindromeMask = END_NOEXT;/* -MaxInvPalindrome D value */

double SplitSegDup = 0.0;/* If > 0.0 : -SplitSegDup <L> value */
int SplitSegDupN = 12; /* -SplitSegDup <N> value */
double SplitSegDupE = 200.0;/* -SplitSegDup <E> value */
int SplitSegDupEN = 5;/* SplitSegDup <EN> value */
double SplitSegDupC = 1.5;/* -SplitSegDup <C> value */ 
double SplitSegDupLM = 180;/* -SplitSegDup <LM> value */ 
int SplitSegDupNM = 18;/* -SplitSegDup <NM> value */
double SplitSegDupMaxOutlier = 5.0;/* -SplitSegDup <Mout> value */

size_t SegDupMask = 0x4;/* If > 0, mark overlap region of segdups (even if smaller than SplitSegDup or SplitSegDupN) with this Mask value for all labels of both maps */

int SplitSegDup_OneSided = 1;/* K value of -SplitSegDup_OneSided <K> <L> */
double SplitSegDup_OneSided_MaxLen = 300.0; /* L value of -SplitSegDup_OneSided <K> <L> */

double BreakEndOutlier = 0.0;/* If > 0.0 : -BreakEndOutlier <L> value */
int BreakEndOutlierN = 20;/* -BreakEndOutlier <N> value */
double BreakEndOutlierE = 100.0;/* -BreakEndOutlier <E> value */
int BreakEndOutlierEN = 10;/* -BreakEndOutlier <EN> value */

double BreakOutlier = 0.0;/* If > 0.0 : -BreakOutlier <L> value */
int BreakOutlierN = 20;/* -BreakOutlier <N> value */
double BreakOutlierI = 100.0;/* -BreakOutlier <I> value */
int BreakOutlierIN = 10;/* -BreakOutlier <IN> value */

int PairMergeHmap = 0;/* If >= 1 : After no more map pairs can merge, re-examing all pairwise alignments and try to combine as many as possible into Haplotype Maps and output them as HMAPs
			      First merge Bubbles (fully overlapped smaller map that could not be merged with larger map due to an internal outlier larger than PairMergeMaxOutlier
			      Next try to combine matched pairs pairs of endoutliers indicated a large incomplete insertion, and treat as a single insertion bubble (will be broken by Haplotyper if not completed).
			      Requires either -pairmergeRepeat OR that pairmerge was previously run until no more map pairs can merge. */
double TruncateHmap = -1.0;/* see -pairmergeHmap option : defaults to Truncate */
int TruncateNHmap = -1;/* see -pairmergeHmap option : defaults to TruncateN */

double PairMergeEndHmap = 0.0;/* If > PairMerge : Minimum overlap required for merging End Overlaps (with Bubbles) : useful to reduce higher risk of chimeric merges */
double MinMapLenHmap = 0.0;/* Minimum map lengths in kb for end overlaps with bubbles : useful to limit higher risk of chimeric merges to high payoff merges likely to increase N50 significantly */



int PairMergeIDrenumber = 0;/* If >= 1 : renumber output contig IDs sequentially from 1 : may have to output contigs that otherwise could be symbolic links */

int XmapInterleaved = 1;/* If >= 1 : Use interleaved label ids in 2-color XMAP alignments and with -Tracksites */


char *ForceMerge = NULL;/* Filename which contains one pairwise merge specification per line (6 integers) */

int SecondBest = 0;/* to output 2nd best alignment score and right end location in .map and .xmap */
double SecondSpacing = 1.0;/* spacing (distance of rightmost aligned sites) of 2nd best alignment from best alignment as multiple of Query Map length */
int SinglelineXmap = 0;/* output 2nd best match on same line as best match (non-standard XMAP format) */

int MultiMatches = 0;/* to output all alignments that score above ScoreThreshold2 & LogPvThreshold2 and have at least MultiMatches unique labels vs any higher scoring match */
double MultiMatchesDelta = 0.0;/* If > 0.0 : Also output overlapped alignments that do not have enough unique labels, provided the offset between maps differs by at least <D> kb */

int MultiMatchesUniqueQuery = 0;/* see -MultiMatchesUniqueQuery 1st arg */
int MultiMatchesUniqueQueryInv = 0;/* see -MultiMatchesUniqueQuery 2nd arg */
int MultiMatchesUniqueQueryTrim = 0;/* see -MultiMatchesUniqueQuery 3rd arg */

int RefSplit = 0; /* If > 0 : Split each alignment (but not maps) at internal outliers with Pvalue < Psplit/bc and at least RefSplitMinLabels site intervals in both maps
		     Note that with RefSplit outliers can consist of any alignment interval with -ve score < log(Psplit) provided no subinterval with score > log(1.0/Rsplit) is included
		     and the first and last subinterval has a -ve score. 
		     NOTE : bc == (outlierBC ? m*n : 1) and m and n are the number of site intervals involved in query and reference respectively.
		  */
double Rsplit = 0.0;/* see RefSplit : If 0, use Psplit instead */
int RefSplitMinLabels = 0;/* see RefSplit */
int MultiMatchesFilter = 0;/* If > 0 : Filter matchgroups by LogPvThreshold2 and AlignedSiteThreshold2 before simimilarity check, to avoid discarding matchgroup with low score but high confidence */
int RefSplitStitch = 0;/* If > 0 : After RefSplit, check if any remaining matchgroups can be stitched together (end to end). Requires that filtering by all thresholds is only applied either before any matchgroups are trimmed
			  (see MultiMatchesFilter) or again at the very end (after stitching any trimmed matchgroups) */
int MultiMatchesRev = 0;/* If >= 1 : Run MultiMatches twice with reference in either orientation and then eliminate duplicate matchgroups */
int MultiMatchesTotScore = 0;/* If >= 1 : Resolve overlapped matchgroups to maximize the total score of unique portions by avoiding internal outliers at the expense of losing the highest scoring single alignment
                                If >= 2 : Do this before applying -MultiMatchesDelta filter, in case alignments have multiple outliers whose sizes cancel each other so the ends have only a small offset.
                                          This takes extra time, so only do this if at least one outlier is MultiMatchesTotDelta kb or larger (should match size of smallest tandem duplication to be detected) */
double MultiMatchesTotDelta = 20.0;/* see MultiMatchesTotScore  */

double MultiMatchesDupesMinscore = 5.0;/* In MultiMatchesDupes() : Minimum absolute score of shared alignment, before breaking up lower scoring alignment */
double MultiMatchesDupesMinscoreRel = 0.5;/* In MultiMatchesDupes() : Minimum scores of shared alignment as share of total alignment score, before breaking up lower scoring alignment */
double MultiMatchesDupesMMPen = 5.0;/* In MultiMatchesDupes() : If MultiMatchesTotScore > 0 : Minimum total score improvment required before breaking up higher scoring alignment */

int bin = 1, numbins = 1, skipcnt = 0;/* partial pairwise alignment parameters */

int Kbiasinit=0,KLbiasinit=0,KFbiasinit=0,Gbiasinit=0,Betainit=0;

int RefRepeats = -1;/* number of times to repeat parameter estimation in refalign or pairalign (for pairalign -M 2+ is same as -M 1)*/
int RefRepeats2 = 1;/* number of times to repeat calls to hash_generate() followed by refalign_allpairs() */
int SaveRepeat = -1;/* If >= 0 : .errbin file will be based on iteration SaveRepeat,RefRepeat2 (instead of last iteration RefRepeat,RefRepeat2) */
double RefRepeatErr = 1e-5;/* stop parameter estimation if no parameter changed by more than this fraction */

int Refine = 0;/* apply refinement step to consensus at end of refalign() and output refined consensus
		  1 : trivial refinement (no contig map change, compute coverage and ChimQuality only)
		  2 : Full refinement, but don't keep maps for regions that will not be refined (see -extonly)
		  3 : Same as 2, but also keep maps fro regions that will not be refined (so they show up in .xmap)
	        */
int RefineRange = 8;/* Adjust RANGE parameters used during refinement */

int RefineXRange = 0;/* If != 0 : extrapolate alignments (when needed) using Xtheta * max(RefineRange , max(deltaXRef,deltaExtXRef)) instead of just Xtheta * RefineRange */
int RefineAlignExt = 0;/* If != 0 : make sure alignments are extended to the ends of each molecule (makes score less dependent on recent reposition call) */
int RefineFillMap = 1;/* see FILL_MAP in refine.h */
int RefineFillHmap = 1;/* see FILL_HMAP in refine.h */
int RefineSwitchFade = 0;/* If != 0 : turn of RefineSwitch* options after initial Haplotyping (when stringency of Indels is increased) so scoring functions becomes almost the same for Indels and SNPs */

int LocalTheta = 0;/* If != 0 : Use local estimate of consensus label interval size instead of Xlambda or Xtheta when adding labels with SCANRANGE,SPLITRANGE,EXPANDSPLIT in refine.h */
double MaxLocalTheta = 0.0;/* If > 0.0, Use this value as max local estimate of consensus label interval size */
int LocalRange = 0;/* If != 0 : Use label interval size of current molecule instead of global Xtheta in setlimit() : see RANGE_Y=0 in refine.h */

double ScanRange = 3;/* see -ScanRange and SCANRANGE in refine.h */
double HScanRange = 3;/* see -HScanRange and HSCANRANGE in refine.h */

double MaxContigSiteDensity = 0.0;/* If > 0 : skip refinment (use -refine 1) for contigs with more than this many sites / 100kb */

long long CMapID = 0;/* CMapID used in output of refined consensus map (If Refine > 0) */
long long startCMapID = 0;/* Value of -id command line option (often a group id of multiple ref ids that will replace CMapID by turn) */
int SplitRef = 1;/* 1 IFF each -ref input contained at most 1 map */
int SplitMap = 1;/* 1 IFF each -i input file contained at most 1 map */
double ContigSplitRatio = 0.0;/* If > 0, split input reference map at internal site interval if sitecntFN drops below ContigSplitRatio * "median fragcnt" (coverage) AND below ContigSplitRatio2 * "median sitecntFN"  */
double ContigSplitRatio2 = 0.0;
double ContigSplitRatio3 = 1.5;/* The split location must have sitecntFN at least this factor below the peak values on either side */
double MinSplitLen = 0.0;/* minimum length of contig split */
int SplitRev = 0;/* Revised -contigsplit code */
int SplitSite = 0;/* If > 0 : Apply -contigsplit only to coverage at site locations (and not between sites, where coverage is more noisy). Particularly useful with -TrimNorm */
long long SplitCnt = 0;/* Revised -contigsplit contig id assignment : Determine new contig ids by repeatedly adding the value in ContigCntFile to the original contig ID */
double SplitFragileD = 0.0;/* If > 0.0 : <D> value of -splitFragile <D> <R> <E> */
double SplitFragileR = 0.0;/*            <R> value of -splitFragile <D> <R> <E> */
int SplitFragileE = 0;     /*            <E> value of -splitFragile <D> <R> <E> */

double HapSitePvalue = 0.0;/* If HapSitePvalue > 0.0 : during refinement haplotype sites will be detected when haplotype Pvalue < HapSitePvalue */
double HapIndelPvalue = 0.0;/* If HapIndelPvalue > 0.0 : during refinement haplotype Indels will be detected when haplotype Pvalue < HapIndelPvalue */
double PhasePvalue = 0.01;/* If HapSitePvalue > 0.0 && PhasePvalue > 0.0: during refinement haplotype SNPs (& indels) will be Phased (with the next one) whenever this can be done with phasing Pvalue < PhasePvaluea */
double HapSiteRes = 0.0;/* If HapSiteRes > 0.0 : Suppress SNPs within this distance (in kb) of another site on same Allele. The value should be related to the imaging resolution res*PixelLen */
double HapSiteResDelay = 0.0;/* If HapSiteResDelay > 0.0 : Suppress SNPs (new or converted) within this distance (in kb) of another site on same Allele during the first HapSiteResN iterations.
			         If HapSiteResL = 2 : Also suppresss SNPs within this distance of another site on the opposite Allele during the firs tHapSiteResN iterations.
				 If HapSiteResA = 1 : Only suppress new SNPs (NOT SNPs created by converting HapSite[i] = 3 -> 1 or 2)
			     */
int HapSiteResN = 11;
int HapSiteResL = 0;
int HapSiteResA = 0;

double HapSiteUnphasedSNR = 0.0;/* Filter out Unphased HapSites with with geometric mean SNR below this value */
double HapSiteUnphasedPvalue = 0.0;/* If > 0.0 : Filter out only those Unphased HapSites with Pvalue > HapUnphased Pvalue */

double HapMinCoverage = 0.0; /* Filter out HapSite or HapIndel if local coverage is below HapMinCoverage (molecules) */
double HapAlleleCoverage = 0.0;/* Filter out HapSite or HapIndel if minor allele coverage as a fraction of local total coverage is below HapMinorCoverage (fraction) */
double HapAlleleProb = 0.7;/* Molecules are assigned to an allele when likelihood (conditional on molecule belonging to either of 2 alleles) exceeds HapAlleleProb */

double HapMinCovSize = 0.0;/* If > 0.0 : Suppress HapAlleleCoverage filter if Hap Indel size is larger than HapMinCovSize (large Hap Indels may have fewer molecules 
			      of the larger Allele, if the smaller Allele was assembled first and Haplotyping was only run during Final Refinement) */
double HapMinCovPvalue = 0.0;/* If > 0.0 : Suppress HapAlleleCoverage & HapMinCoverage filter if Hap SNP has Pvalue < HapMinCovPvalue */
double HapAlleleCovInt = 0.0;/* Increase HapAlleleCoverage for Indel by HapAlleleCovInt * "largest label interval in region of indel / 100kb" */
double HapMinCovPvalueL = 0.0;/* If > 0.0 : Suppress HapMinCoverage filter if Hap Indel has Pvalue < HapMinCovPvalueL */
double HapMinCovAbs = 0.0;/* If > 0.0 : Always filter out Hap Indels if Minor Allele absolute coverage drops below this value (except if HapMinCovAbsM,HapMinCovN disables it) */
double HapMinCovAbsM = 99.0;
int HapMinCovN = 999;/* Never filter out Hap Indels overlapped by more than this many SNPs AND more than HapMinCovAbsM SNPs/100kb */


double HapMinCovPL1 = 1.0;/* If < 1.0 : Don't filter out HapIndel if ChiSquare Pvalue for sizes of minor allele is above HapMinCovPL1 */
double HapMinCovPL2 = 0.0;/* If > 0.0 : Always filter out HapIndel if ChiSquare Pvalue for sizes of minor allele are below HapMinCovPL2 */
int HapMinCovL = 0;/* If 0 : Both versions of Pvalue must agree (both must exceed PL1 or be below PL2) before -HapMinCov action is overridden 
		      If 1 : Only Pvalue based on autoNoise sizing error estimate is used.
		      If 2 : Only Pvalue based on major allele sizing error is used */
double HapMinCovCL = 0.0;/* If > 0 : Minimum coverage of minor Allele required for -HapMinCovPL override to occur */

int checkMinorIsBigger = 1;/* If 0 : Allow HapMinCov to apply to Hap Indels even if minor Allele is bigger in size
			      If 1 : HapMinCov does NOT apply to Hap Indles IF minor Allele is bigger in size (assumes DNA knots only shrink DNA)
			      NOTE : HapMinCov filter always applies of total coverage is below HapMinCoverage */

double HapEndTrim = 0.0;/* Trim Allele A ends up to the first label with HapEndTrim < (covA + N) * (covA+covB+covC + 2N)/(covA+covB + 2N), where N = HapEndTrimMin */
double HapEndTrimMin = 1.0;/* Pseodo-Count correction for covA and covB */
double HapEndTrimMaxFrac = 0.5;/* Maximum fraction of Labels that may be trimmed (otherwise undo all per-Allele trimming and rely on -EndTrim) */
double HapEndTrimRatio = 0.20;
double HapEndTrimLen = 99.9;

int FastSnpsStartup = 0;/* If > 0 : If SNPs/Indels are already present start Haplotyping with this iteration instead of 0 (should be no more than HMIN_ITER-1) */
int FastSnps = 1;/* Don't apply fast SNPs detection iteration for the first FastSnps iterations : should be at least 1 to avoid phasing errors */
int FastIndel = 1;/* Don't apply fast HapIndel detection iteration for the first FastIndel Indel iterations : should be at least 1 to avoid phasing errors in absence of SNPs */
int HapIndelHiRes = 1;/* Number of interations to use high res sizing ladder before reverting gradually to 1.4142 x sizing ladder. Larger values improve convergence at the expense of runtime */
int MAX_DELTA_OVERSAMPLE = 2;/* how many times the number of samples can be reduced 2x after initial Indel scan has been completed (see LP_INDEL_SKIP)
				Original sizing ladder step ratio = pow(2.0, 1/(2<<MAX_DELTA_OVERSAMPLE)) */
double MaxHapDelta = 200.0;/* largest Delta (half of largest haplotype Indel size) in kb that is checked during Haplotypeing (any value above 175 may be treated same as 175 kb) */

int SplitHmap = 0;/* If 1 : Input maps that are already Haplotyped are split into two maps corresponding to the 2 Alleles (See RefineHmap on two ways to proceeed thereafter)
		    If 0 : Input maps that are already Haplotyped are treated as a regular Cmap that is the average of the 2 Allele maps.  */
int RefineHmap = 0;/* Requires SplitHap == 1:
		     If 0 : Input maps that are Haplotyped are separated into 2 Allele maps and each one Refined + Haplotyped independently of the other.
		     If 1 : Input maps that are Haplotype are preserved and Haplotyped further starting from the previous Hap map. In this case non-Haplotype refinement
		            will NOT be performed before Haplotype refinement
		  */

int HapSiteMin = 1, HapIndelMin = 1;/* Only output separate <N>1.cmaps and <N>2.cmaps and <N>0.hmap IF at least HapSiteMin different Haplotype Sites and HapIndelMin different Haplotype Indels were found 
				       (other wise output just one <N>0.cmap). Decision is made after filtering HapSites and HapIndels based on -HapUnphased or -HapMinCov */
int HapSiteWin = 0; /* Only output seperate Allele maps if at least HapSiteWin Het Sites are present in some consensus interval of size HapSitewinSize (kb) or less */
double HapSiteWinSize = 100.0;
double HapIndelMinSum = 0.0; /* Sum of absolute Hap Indel sizes must add up to at least HapIndelMinSum, otherwise ignore Hap Indels when deciding if seperate Allele maps should be output */
double HapIndelMinDensity = 0.0;
int HapIndelSkipEndSNP = 1;
double MinHapCovThresh = 0.0;/* see -HapThresh 8th arg */
double MinHapRelCovThresh = 0.0;/* see -HapThresh 9th arg */

int HapIndelMerge = 0;/* If > 0 : Treat consecutive Hap Indels seperated by Hap Sites as a single indel with a single HapIndelPvalue penalty */

int HapIndelDelskip = 0; /* If > 0 : Perform expensive Hap Indel deletion step then skip next HapIndelDelskip times if LP change was less than HapSkipEps + HapSkipEpsPermap * MD */
int HapUpdatemapSkip = 0;/* If > 0 : Perform expensive Hap Updatemap then skip next HapUpdatemapSkip times if LP change was less than HapSkipEps + HapSkipEpsPermap * MD */

double HapSkipEps = 1.0;
double HapSkipEpsPermap = 0.0;

double IntervalEps1 = 1.0;/* If LP change for all interval size changes is less than this amount before -rres & -cres, terminate */
double IntervalEps2 = 2.0;/* If LP change for all interval size changes is less than this amount after -rres & -cres, terminate */
double HapIndelEps1 = 1.0;/* If LP change for all Hap Indel size changes is less than this amount before -rres & -cres, terminate */
double HapIndelEps2 = 1.0;/* If LP change for all Hap Indel size changes is less than this amount after -rres & -cres, terminate */
double HapLabelEps1 = 0.1;/* if LP change for all Hap Label changes is less than this amount before -rres & -cres, terminate */
double HapLabelEps2 = 1.0;/* if LP change for all Hap Label changes is less than this amount after -rres & -cres, terminate */

int HapFast = 0; /*  If >= 1 : Speed up Haplotyping by recomputing likelihood only where consensus was changed (will still occassionally check full likelihood and backtrack if needed)
		     If >= 2 : Avoid global re-alignment of molecules to handle large changes on consensus  */
int HapLocalStop = 0;/* Terminate Haplotyping locally in large contigs and only work on regions that are still making progress in at least one Allele, based on -HapMinLL */
int HapLocalScan = 1;/* value of DELTA_SCAN in refine.cpp */
int HapLocalScan2 = 16;

int HapStageIter = 16;/* see -HapMaxIter 16 0.25 */
double HapStageFrac = 0.25;

int HREVERSAL = 2; /* see -HapPhaseRev 2 1 */
int FINAL_HREVERSAL = 1; 
int HAP_MAX_DERES = 2; /* see -HapCres 2 0 */
int DERES_STOP = 0;
int FILTER_STOP = 0; /* see -HapFilterStop 0 */

int SKIPMAP_ADD = 2;// HERE HERE : make command line variable

int HapMaxDeltaIter = 5;/* Maximum number of iterations in a row with only sizing changes (no HapIndel or SNP changes) before termination is forced */
int HapMaxSwitchIter = 100;/* Maximum Hap Indel iterations before *Switch options must be turned off */

int HapFilterLast = 0;/* If >= 1 : Apply Haplotype filters last (after -rres,-cres based label merging) */
int HapMergeSNPs = 0;/* If >= 1 : Merge nearby SNPs on opposite Alleles into a regular site, if score improves */


char *ContigCntFile = 0;/* If ContigSplitRatio > 0, this file will contain a single integer <N> that must be atomically incremented to obtain the 
			   name of any additional output contig as PREFIX_contig<N>_refined.cmap */

int Cfirst = 0; /* see option -first */
int CfirstPairs = 0;/* see option -firstPairs */

int EVERB =2; /* <= 1 : Enable display of each missing/false cut (in last iteration)
		>= 2 : Disable display of each missing/falst cut */

double err_dist = 2.0;/* count errors (fp and fn) that are isolated from true positives by at least this distance in kb */

char *hash_prefix = 0;/* output .hash or .hashbin file name prefix */
char *hash_filename= 0;/* input .hash or .hashbin file name */

int hash_threshold = 0;/* hash score threshold : only use pairs with hash score at least this large */
int hashsplit = 0;/* If != 0 : Program just splits hashtable into multiple pieces */
double hashdelta = 0.0;/* If != 0 : limit relative offset of aligned maps to value in hashtable +- this multiple of SD of length of smaller map */
double hashdeltaAdjust = 0.0;/* If != 0 : increase hashdelta for a particular map pair, if needed, in increments of hashdeltaAdjust */
double hashdeltaLim = 0.0; /* If != 0 : Limit increases in hashdelta to no more than hashdeltaLim */
int HASHRANGE = 0;/* If != 0 : base hashdelta range on number of sites rather than distance : faster when reference label density is higher than query label density */

int scale_output = 0;/* output scaling factor for each Molecule in .map file */

int COVERAGE_TRIM = 3; /* skip first and last few aligned sites of each Molecule when computing ChimQuality site coverage : this will expose shallow regions spanned only by chimeric molecules 
			     that tend to have few aligned sites on the "bad" side of the chimeric junction (see refine.cpp) */
double COVERAGE_TRIM_LEN = 0.0;/* Also skip at least this many kb of each aligned map end */
int mCNT = 0;/* number of alternate values for COVERAGE_TRIM_LEN */
double mCOVERAGE_TRIM_LEN[MAX_CHIM];/* If mCNT > 0 : multiple alternate values for COVERAGE_TRIM_LEN : mCOVERAGE_TRIM_LEN[0] == COVERAGE_TRIM_LEN */

double TrimNorm = -1;/* If >= 0 : Normalize ChimQuality by total number of molecules that extend beyond both sides of the target site by at least COVERAGE_TRIM_LEN kb and COVERAGE_TRIM real sites (must be > 0) :
		              Coverage will be expressed on scale 0-100 (%). The TrimNorm value is used as a pseudocount prior to avoid breaking low coverage regions (eg fragile sites) */
double TrimNormLen = 0.0;/* If > 0 : Use largest normalization count for TrimNorm within TrimNormLen kb */

int TrimNormChim = -1;/* If >= 0 : Try to suppress counting chimeric molecules during TrimNorm by not counting molecules whose nearest aligned site is more than TrimNormChim real sites (in Hcuts[]) away from the target site
			This will avoid breaking chimeric sites, but if set too small may start to fail breaking chimeric contigs as well due to false negative sites. A value of 1 or 2 should be OK */
int TrimNormFP = -1;/* If >= 0 : Try to suppress counting chimeric molecules during TrimNorm by not counting molecules whose alignment ends more than TrimNormFP sites beyond target site, but does NOT extend CovTrimLen/COVERAGE_TRIM sites beyond target site */

int ChimQualityFix = 0;/* some updates to ChimQuality computation */

double FragileQualityEnd = 0.0;/* If != 0.0 : Compute Fragile Site Quality Score, by counting molecules whose ends are within this distance (in kb) from the suspected fragile site AND have no more than FragileQualityFN sites near the end or no more than FragileQualityFP sites beyond the suspected fragile site */
int FragileQualityFN = 0;
int FragileQualityFP = 0;

int OutlierQuality = 0;/* If != 0 : Compute Outlier Quality score, by counting fraction of molecules with outliers crossing each interval
			  If >= 2 : Also compute Variance Quality score, by computing sample variance of non-outliers crossing each interval */

int HapMapWT = 0;/* If != 1 : Give full weight to maps in both Alleles, so that Quality scores of each Allele contig reflect all molecules, includes from the other Alelle */

double TrimNormMin = 0.0;/* If TrimNorm >= 0 : Force Normalized Trimmed Coverage to zero if type N1 coverage is below specified value */
double TrimNormMinCov = 0.0;/* If > 0.0, only apply TrimNormMin if raw coverage is above TrimNormMinCov */
double TrimNormMinEnd = 0.0;/* Don't compute quality scores near ends of contigs until N1 coverage is above TrimNormMinEnd */
double TrimNormMinRel = 1.0;/* If < 1.0 : only apply TrimNormMin if N1 is below TrimNormMinRel x "max value of N1 in either direction" */

int TrimNormEnd = 0;/* If TrimNorm >= 0  && TrimNormEnd > 0 : Don't include molecules in Trimmed Coverage (numerator) if either end has an endoutlier with at least TrimNormEnd unaligned sites  */
int TrimNormMed = 0;/* If TrimNorm >= 0  && TrimNormMed > 0 : Use TrimNormMed instead of median of Trimmed Coverage to compute threshold (multiplied by ContigSplitRatio2)  */
int TrimOutlier = 0; /* If > 0 : For molecules alignments with outliers, treat each aligned region between outliers as a seperate molecule (Helps break consensus with spurious indels caused by coincident outliers, but this can also be handled with -outlierLambdaRef )*/
double TrimOutlierKB = 0.0;
int TrimOutlierLabels = 1000;

int REPLACE_COV = 0;/* If > 0 : output trimmed coverage (see COVERAGE_TRIM) in place of actual coverage in output_draft() */

double TrimFactor = 0.0;/* used by Refine as initial trimfactor for estimated trimmed mean of interval sizes */

int MultiMode = 1;/* Enable Multi Modal sizing optimization */
int Mprobeval = 1;/* Enable new faster refinement code */
int MprobevalTest = 0;/* Enable slower more thorough Mprobeval refinement by adjusting a number of hardwired parameters */

char *MappedPrefix = 0;/* If !=0 : output those -i input maps that align with reference above threshold to <MappedPrefix>.cmap */
char *UnMappedPrefix = 0;/* If !=0 : output those -i input maps that did NOT align with reference above threshold to <UnMappedPrefix>.cmap */
char *PartMappedPrefix = 0;/* If !=0 : output those -i input maps that were only partly aligned (see -endoutlier) with a reference to <PartMappedPrefix>.cmap */
int MappedUnsplit = -1;  /* If -1 use default heuristic, otherwise specifies if mapped output is NOT split (1) or split (0) into seperate files per reference contig (ignored if -grouped was specified) */
char *groupfile = 0; /* If !=0 : mapped output is split into seperate files based on groups of reference contig ids specified in groupfile (one group per line, 1st line is group1 etc) */

double CenterMapped = 0.0;/* If != 0 then MappedPrefix output is limited to only those -i input maps that align inside the center region of contigs not within this distance of contig ends */

int MaxCov = 0;/* If > 0 : limit coverage at each consensus location to <C> by discarding maps with the lowest alignment Pvalue */
int MaxCovWt = 0;/* If > 0 : Use weighted coverage instead of raw coverage with MaxCov : weights are based on -BestRefWt */

int PoutlierEndCnt = 0;
double PoutlierEndRefine[16];/* PoutlierEndRefine[0..PoutlierEndCnt-1] is the sequence of outlier end probabilities used during successive phases of refinement if extend >= 2.
				 The values of PoutlierEnd will be added at the start unless the first value is smaller. The value 1.0 will be used for the final alignment and coverage computation.
			         Using a sequence from stringent to least stringent is useful to weed out mismatched extensions and (hopefully) converge to the consensus estimate of the extension */
double PoutlierEndFinal = 0.0001;/* Value of PoutlierEnd used during final alignment computation (if greater than PoutlierEndRefine[PoutlierEndCnt-1]) : does not change the consensus map but
				    can change the coverage profile and help reveal chimeric junctions by using PoutlierEndFinal = 1.0 */
int PoutlierCnt = 0;
double PoutlierRefine[16];/* matching series of Poutlier values corresponding to PoutlierEndRefine. Lengths must match */

double biasWT = 1.0;/* weight of theoretical score bias in refalign by this factor (1 = full bias) (0 = no bias, refalign score will match Viterbi version of refine score) */
double PRbiasWT = 0.0;/* weight of theoretical score bias correction for misresolved sites (1 = full bias, 0 = no bias) */
double biasWTend = 1.0;/* If 1 : apply biasWT to endoutliers (this keeps the total Bias for each molecule independent of the alignment, but makes it harder to fairly score chimeric molecules or references) */
double biasWToutlier = 1.0;/* If 1 : apply biasWT to outliers this keeps the total Bias for each molecule independent of the alignment, but makes it harder to fairly score molecules with outliers or SVs) */
double biasWTrefine = 0.0;/* biasWT used during refinement except during initial extension refinement (NOTE : biasWTend is always 0 during refinement, but biasWToutlier is same as during refalign) */
double PRbiasWTrefine = 0.0;/* weight of theoretical score bias correction for misresolved sites (1 = full bias, 0 = no bias) */
int biasWTrefineType = 0;/* 0 : Use biasWTrefine only when computing final coverage and Chimeric Quality scores
			    1 : apply biasWTrefine during refinement (except during initial refinement)
			    2 : apply biasWTrefine during refinement (including initial refinement)
			 */
int RTHETA_FIX = 1;/* Add sqrt(MD/sitecnt) term to FL() and 1/sqrt(MD/sitecnt) for FLE() etc : reduces score of endoutliers vs refalign (turn off to better match endoutliers in refalign) */

double DELTA_RES= 0.01; /* resolution in kb of MultiMode sizing optimization (2x faster if resolution is set to 0.01, but LP does not converge to global maxima as well as with 0.001) */
double DELTA_REL= 0.01;/* relative resolution of MultiMode sizing optimization */
double Site_Pen= 3.0;/* penalty per consensus map site in LP computation : requires to balance Maximum Likelihood bias in favor of too many sites (cf degrees of freedom penalty) */
double Site_Pen_Ext = 0.0;/* penalty Site_Pen used during initial refinement stage : typically lower than Site_Pen */

double FN_Refine = 0.0;/* If > 1.0 : Multiply FN value by this factor during refinement */
double FN_Ext = 0.0;/* If > max(1.0,FN_Refine) : Multiply FN value by this factor during initial refinement stage */
double FNmin = 0.0;
double FNminExt = 0.0;

double MinSF_Refine = 0.0;/* If > 0.0 : increase sf value to at least this value during refinement */
double MinSF_Ext = 0.0; /* If > MinSF_Refine : increase sf value to at least this value during initial refinement stage */

double MinSR_Refine = 0.0;/* If > 0.0 : increase sr value to at least this value during refinement */
double MinSR_Ext = 0.0; /* If > MinSR_Refine : increase sr value to at least this value during initial refinement stage */

double skip_dist= 0.0;/* speed up addition/deletion of sites during refinement by checking only every skip_dist (kb) */
double skip_dist_Ext = 0.500;/* If > skip_dist : value of skip_dist used during initial refinement stages */
double skip_score = 20.0;/* speed up addition/deletion of sites during refinement by NOT rechecking changes that drop LP by LP_SKIP as frequently */
double skip_score_Ext = 10.0;/* If < skip_score : value of skip_score used during intial refinement states */

int dosvdetect = 0; /* call output_smap */

int smap_zygosity = 1; //put zygosity flag in smap
int smapSize = 1;// add SVsize in smap
int smapFreq = 0;// add SVfreq in smap
int smapTransOrient = 0;// add Orientation column in smap (current NA unless it is a translocation)
double smap_Indel_Ref_Alignment_OverlapFrac = 0.5;/* minimum fractional overlap in ref of Alignment (without SV) with an SV (in another contig/Allele) to flag the SV as heterozygous */
double smap_Indel_Ref_SV_Overlap = -20.0;/* minimum overlap in ref (in kb) of 2 SVs for them to be considered the same. Can be -ve, which means a gap between ref regions is permitted (useful to handle repeat expansion or contraction, which may be called ambigously in different Alleles). */
int smap_Del_Ref_SV_OverlapIncrease = 1;/* For deletions increase value of sma;_Indel_Ref_SV_Overlap by the size of the deletion */
double smap_Indel_Size_MatchRatio = 0.8;/* max ratio of Indel Sizes (smaller size / larger size) for 2 indels to be considered the same */

double smap_confbylen = 0.07; /* smap: minimum confidence (in xmap) per kb length of matchgroup : matchgroups that fail are not used to call any SVs */
double smap_ConfByInterval = 0.0;/* smap: minimum confidence (in xmap) per aligned Label Interval of matchgroup : matchgroups that fail are not used to call any SVs */
double smap_side_alignkb = 50.; /* smap: minimum length (kb) of aligned regions (ie, xmap entries) on both sides of an sv */
double smap_sv_sizekb    = 5000.; /* smap: maximum abs difference in query and reference gap size for indels & inversions, otherwise they are changed into intra_trans */
double smap_inv_MaxGap = 0.;/* If > 0 : maximum reference or query gap size for inversions, otherwise they are changed into intra_trans */
double sv_indel_maxoverlap = 50.0; // max size (kb) of ref overlap allowed (& removed) for indels : Obsolete, now only used in pairInversionBP (see smap_indel_maxqryoverlap instead)
double smap_del_maxqryoverlap = 140.0;// max size (kb) of query overlap allowed for deletions : larger overlaps are assumed to be FP repeat (or segdup) compression
double smap_del_bigqryoverlap = 110.0;// If deletion query overlap is between this value (kb) and smap_del_maxqryoverlap, reduce confidence of deletion to minimum value (smap_min_conf)

double svInsMaxGapOverlap = 0.5;// max fractional overlap of qry gap allowed for insertions 
double svDelMaxGapOverlap = 0.5;// max fractional overlap of ref gap allowed for deletions

int    smap_min_end_nlabels = 5; /* smap: minimum number of un-aligned labels required for an 'end' type SV */
double smap_min_end_lenkb   = 50.; /* smap: minimum un-aligned length (kb) required for an 'end' type SV */
//double smap_inversion_ratio = 0.0; /* smap: the ratio of gap on the reference vs query must be > this */
double smap_translocation_maxsep = 500.; /* smap: for transloactions with no overlap, maximum separation on query */
double smap_translocation_ovlfr  = 0.3; /* smap: for translocations with overlap, maximum overlap length/matchgroup length */
double smap_translocation_ovlmx  = 140.; /* smap: for translocations with overlap, maximum overlap length */
double smap_InvMaxOverlapSize = 140.0;/* smap: If > 0.0 : For inversions with query overlap, maximum overlap length in kb */
double smap_indel_tiny_maxsize   = 1.; /* smap: tiny type indel is any indel whose size is <= this, in kb */
double smap_translocation_mask_tol = 30.; /* smap: when using mask file, add this as fixed tolerance in kb */
double smap_min_conf = 3.0; /* smap: do not report SVs with conf < this (unless not defined) */

int smap_DupSplit_maxoverlap = 2;/* maximum query overlap (in Labels) for duplication_split calls */
int smap_DupSplit_maxgap = 30;/* maximum query gap (in Labels) for duplication_split calls */
int smap_DupSplit_maxend = 3;/* maximum labels in query beyond the ends of the two matchgroups for duplication_split calls */
double smap_DupSplit_MinQryLen = 200;/* minimum length in query matchgroup to allow query to be extrapolated for a Duplication call, provided the unaligned portion of the end being extrapolated is less the smap_DupSplit_maxend sites */
int smap_DupSplit_MinQrySites = 15;/* minimum labels in query matchgroup to allow query to be extrapolated for a Duplication call, provided the unaligned portion of the end being extrapolated is less the smap_DupSplit_maxend sites */
int smap_DupSplit_singleEnd = 3;/* 0 : duplication_split calls with one end that cannot be extended (due to smap_DupSplit_maxend OR ~MinQryLen OR ~MinQrySites), are NOT allowed.
				   >= 1 : duplication_split calls with one end that cannot be extended, are allowed if MGs overlap on reference
				   >= 2 : duplication_split calls with one end that cannot be extended, are also allowed if nonextendable MGs is at least 50% in kb of total ref interval
				   >= 3 : duplication_split calls with one end that cannot be extended are always allowed.
				 */

int smap_Outlier_minQrySites = 4;/* For Small inversion, duplication or translocation  the minumum number of sites (1 + intervals) in Outlier of larger Matchgroup */

int smap_Dup_minRefOverlapSites = 4;/* For Duplication the reference overlap (in sites) must exceed the Qry overlap (in sites) by at least this number */
double smap_Dup_maxQryGapMult = 0.8;/* Max Qry Gap size (in sites) as a multiple of the reference overlap (in sites) to call a Duplication */
double smap_Dup_maxSpacing = 0.0;/* Maximum Gap (in kb) on Query between Duplication events (should correspond to max distance that can be phased) */

double smap_InvDup_minRefOverlapFrac = 0.5;/* For Inverted Duplication the reference overlap (in sites) must be at least this fraction (F) of the smaller matchgroup reference size (in sites)
					      Also used to suppress inversions if overlap fraction in Sites is >= F for either ref or qry (but see -svInv_maxRefOverlapFrac & -svInv_maxQryOverlapFrac) */
double smap_Inv_maxRefOverlapFrac = 0.0;/* If > 0 : suppress inversions if overlap fraction in Sites on Reference is >= F */
double smap_Inv_maxQryOverlapFrac = 0.0;/* If > 0 : suppress inversions if overlap fraction in Sites on Query is >= F */


double smap_InvDup_maxQryGapMult = 0.8;/* For Inverted Duplication the maximum qry gap (or overlap) in sites cannot exceed this fraction of the reference overlap (in sites) */
double smap_InvDup_maxSpacing = 0.0;/* Maximum Gap (in kb) on Query between Inverse Duplication events (should correspond to max distance that can be phased reliably ???) */

double smap_Inv_MaxTrimFrac = 1.0;/* In deciding between Inverted Duplication or Inversion, rule out Inversion if more than this fraction 
				     of Inversion matchgroup must be trimmed away or ignored for Inversion call */
double smap_Inv_MinTrimFrac = 0.0;/* In deciding between Inverted Duplication or Inversion, rule out Inverted Duplication if less than this fraction 
				     of Inversion matchgroup must be trimmed away or ignored for Inversion call */

int smap_Inv_minRefSites = 4;/* minimum number of non-overlapped sites in reference of inverted matchgroup needed to call inversion */
int smap_Inv_minQrySites = 4;/* minimum number of non-overlapped sites in query of inverted matchgroup needed to call inversion */

double smap_Inv_ConfMinIntervalSize = 1.0;/* minimum interval size to include in interval count to determine confidence of small paired inversions */
int smap_Inv_ConfMinNumIntervals = 2;/* minimum number of intervals (above smap_Inv_ConfMinIntervalSize) to output small paired inversions */

int smap_IndelConfBySize = 0;/* If 1, base Indel confidence on size difference only (ignoring misaligned labels) */
double smap_SmallDelConfirm = 0.0;/* For Deletions in ref intervals below this size, confirm deletion by including neighboring intevals and update center confidence if smaller */
double smap_IndelRangeExpandRes = 0.8;/* max range (in kb) of misresolved labels allowed at boundary of Indel ref range : expand range until this condition is satisfied and adjust Indel size estimate downward, if needed */
double smap_IndelRangeExpandSep = 1.7;/* min seperation (in kb) of boundary labels from next ref label (in either direction) : expand range until this condition is satisfied and adjust Indel size estimate down */
double smap_IndelRangeMinSVchange = 0.20;/* min relative change in SV size due to range expansion : If SV size does not change by this amount the range expansion is undone (but SVsize update retained) */

double smap_GapOverlapRatio = 0.0;/* minimum fractional overlap of 3rd MG by gap in both ref and qry gap of a MG pair with same qry and ref, before MG pair is not a valid SV candidate */

double smap_PrimaryOrientation = 0.0;/* A MG has a dominant orientation IFF total Qry sizes for all MGs with same qryid and refid in one orientation exceeds the total Qry sizes for the other orientation
					by this ratio. */
double smap_TransMaxOverlap = 0.7;/* If > 0 : Disallow translocation if 3rd MG overlaps one of the 2 trans MGs by at least a fraction R of the trans MG qry size and maps back to the other reference
				       NOTE : IF 3rd MG confidence exceeds minimum confidence for translocation MGs, include overlap with the gap in the total overlap fraction */
double smap_TransMaxGapOverlap = 0.7;/* If > 0 : Disallow translocation if a 3rd MG overlaps the query gap of the 2 translocation MGs by a fraction G of the 3rd MG qry size.
					NOTE : 3rd MG confidence must exceed minimum confidence for translocation MGs */

double smap_Repeat_Tolerance = 1.7;/* SV simple repeat tolerance as multiple of expected sizing error (sqrt(SF*SF+SR*SR*X*X), X = average of intervals) */
int smap_Repeat_Minele = 2;/* Minimum consecutive repeat intervals */
double smap_Repeat_ConfPenalty = 1.5;/* If > 0 : Compute pro-rata confidence of non-repeat intervals only and subtract this penalty per repeat region AND only mark with repeat if confidence drops below -T */

int smap_MATCHGROUP_TRIM = 3;/* see MATCHGROUP_TRIM in output_smap.cpp */
int smap_MATCHGROUP_PARTIAL_TRIM = 1;/* see MATCHGROUP_PARTIAL_TRIM in output_smap.cpp */

char *conf_filename = NULL;/* filename for input of confidence scaling table */
double svConfMinsize = 0.500;/* If indel size is below this output svConfNAvalue instead of value from confidence scaling table */
double svConfNAvalue = -1.00;

double Palindromic_Delta = 2.0;/* mark xmap entry as palindromic if best alignment in opposite orientation had confidence within this Delta of best orientation */

int doscaffold = 0; /* process scaffold file */

char *svcheck = 0;/* check SVs in input .indel or .smap (currently only indels) by alignment -i input maps against two versions of each SV neighbor hood and compute a Likelihood ratio for/against the SV being present.
		     The consensus maps are read in from the query map file named in svcheck[] .indel/.smap file with the corrected regions read in from the reference map file named in svcheck[].
		     Outputs corrected versions of each query map as <prefix>_contig<NNN>.cmap along with updated Pvalues for SVs checked in <prefix>.indel (otherwise same as input file svcheck[]) */
double svcheckKB = 300.0;/* SV neighborhood size in kb (on each side of SV) */
double svcheckMerge = 200.0;/* Try to merge SVs that are seperated by less than this many KB on the same alignment */
double svcheck_confidence = 0.0;/* correct SV if its updated confidence drops below this value */

char *svcompare1 = 0, *svcompare2 = 0;/* compare SVs in these two .indel or .smap files : display total SVs and number that overlap by at least SVcompareKB bases or SVcompareFrac (fraction) */
double SVcompareKB = 0.001, SVcompareFrac = 0.10, SVDeltaMatch = 0.0, SVConf = 0.0, SVXmapOverlap = 0.0;
int SVSiteMatch = 0;
char *SVfile = NULL;

int RepeatMaxShift = 0;/* Maximum alignment shift (in aligned sites) checked for repeat regions */
double RepeatPvalueRatio = 0.5;/* Ratio of sizing error (Chi-Square) P-values for shifted alignment vs original alignment required to confirm repeat region */
int RepeatRec = 0;/* Use improved RepeatMask implementation */
double RepeatLogPvalueRatio = 1.0/* TRY 0.7 */;/* minimum Ratio of Log10(Pvalue) for shifted alignment vs original alignment required to confirm repeat */
double RepeatKbShiftRatio = 0.0 /* TRY 0.5 */; /* minimum offset shift (kb) as fraction of alignment right end offset shift : this must be satisfied along entire alignment */
double RepeatKbShiftRatio2 = 3.0;/* maximum offset shift (kb) as fraction of alignment right end offset shift : this must be satisfied along entire alignment */
int RepeatMaxLabelDrop = -1;/* If >= 0 : maximum drop in aligned labels for shifted alignment for repeats is Shift + RepeatMaxLabelDrop */
double RepeatLogPvalueRatioPerLabel = 0.0;/* The minimum Ratio of "Log10(Pvalue) per label interval" for the shifted alignment vs the original alignment required to confirm repeat */

double RepeatShift = 0.0;/* With pairwise alignment if > 0 : check for repeats in each Molecule by aligning it with itself with a minimum offset of RepeatShift */

int ScaleDelta = 0;/* Scale Molecule size up/down by D = ScaleDeltaSize * bpp * (1 ... ScaleDelta) (by multiplying sizes by (1+D) and 1/(1+D)) */
double ScaleDeltaSize = 0.05;
int ScaleDeltaBPP = 0;/* If != 0 : Correctly handle bpp with ScaleDelta by treating Scaling as temporary (to find the alignment) and estimating parameters based on unscaled molecules */

int MapScale = 0;/* If > 0 : Apply per map scaling correction to all -i input maps. Applied before each -M iteration based on the previous iteration's alignment (excluding outliers).
			 The first -M iteration is run twice using the first sub-iteration to apply a scaling correction */
int MapScaleDelta = 0;/* replaces ScaleDelta during first sub-iteration of alignment if MapScale > 0 */
double MapScaleDeltaSize = 0.02;/* replaced ScaleDeltaSize during first sub-iteration of alignment if MapScale > 0 */

int ForceOverwrite = 0;/* 0 : Exit with error if output file already exists.
			  1 : Force overwrite of output files that already exist */
int Mfast = 0;/* see FAST in refalign.cpp */
double MapRate = 0.0;/* If > 0 : Limit Mapping rate to specified value by making -T more stringent and keeping only the best alignments at each iteration of -M */
double MaxLogPvThreshold = 10000.0, MaxScoreThreshold = 10000.0;/* Limit -T or -S adjustment to be no more stringent than these values */

int relaxBnds = 0;/* see -relaxBnds */
int biaswtDefer = 0;/* see -biaswtDefer */


int BestRef = 0;/* only use the alignment with the best refid for each mapid in refalign and refalign2, based on score */
int BestRefPV = 1;/* If 1 : modify BestRef to pick the best refid based on LogPV */
double BestRefPV_MaxConf = 9999.0;/* limit LogPV used to no more than this value when computing molecule weight (see -BestRefWT) */

double BestRefMargin = 0.0;/* If BestRef : prefer earlier contig ids over later contig ids by this margin of score (or logPV, if BestRefPV) */
double BestRefExt = 0.0;/* If > 0.0 : Only suppress extension alignments IF there is a center alignment (on a another ref contig NOT within BestRefExt (kb) of ends) with a log10(Pvalue) that
			   exceeds that of the extension alignment by BestRefExtMargin (can be -ve) */
double BestRefExtMargin = 0.0;
double BestRefOutlier = 1000.0; /* Only suppress alignments with internal outliers sizing error above BestRefOutlier */
int BestRefWT = 0;/* Alternative to BestRef : Only downweight molecule alignments based on relative likelihood (exp(score)) OR Pvalue */
int AddFragsWT = 0;/* Modify BestRefWT by adding alignment weights of all alignments of the same query map with the same ref map in the same orientation, treating them as fragments of one large alignment */

int XMAPmapWT = 0;/* see -XMAPmapWT : If 1 output mapWT column in .xmap right of alignment pairs */

int SkipIDcnt = 0;
long long *SkipID = 0;/* with -merge exclude specified set of map IDs from output */

int SelectIDcnt = -1, SelectIDmax = 0;
long long *SelectID = NULL;/* with -merge output a subset of maps based on the specified set of map IDs */

int SelectScanId = -1;/* If >= 0 : select a subset of maps that have a specified UniqueScanID in the range SelectScanId .. SelectScanId2 (UniqueScanID is equal to RunId field in Saphyr BNX) */
int SelectScanId2 = -1;

int BreakIDinc = 0;
int BreakIDcnt = 0;
long long *BreakID = NULL;/* IDs to be broken into fragments are BreakID[i=0..BreakIDcnt-1] */
int *BreakPointCnt = NULL;
double **BreakPoint = NULL;/* Breakpoint locations are BreakPoint[i=0..BreakIDcnt-1][0..BreakPointCnt[i]-1] */

int WinBreak = 0;/* If > 0 : break apart each input contig into overlapped fragments of size WinBreakLen (kb) each shifted by WinBreakShift (kb) from the previous one
		    OutputID = InputID * WinBreak + PieceID */
double WinBreakLen = 0.0;
double WinBreakShift = 0.0;

long setMask = -1;/* If >= 0 : set all map ends to have this mask value */
long orMask = 0;/* If > 0 : OR value into all map end mask values */

int BestAlignments = 0;/* If > 0 : output only the best BestAlignments alignments based on logPV */
int FirstAlignments = 0;/* If > 0 : output the first FirstAlignments alignments that are above all 3 thresholds */ 

int PVres = 0;/* If > 0 : Apply a correct to Pvalue to reflect limited resolution (-res & -bpp) */
int AlignRes = 0;/* If > 0 : use more thorough -res based alignment (may result in 3x as many alignments and 10x slower refinement) */
double AlignResMult = 1.0;/* If AlignRes > 0 : limit -res range to (res+SDrange*resSD)*AlignResMult */

double EndScoreThresh = 10.0;/* If pairwise score without right end score (Send(),Sbnd()) is below bestscore (so far) by this threshold (as multiple of missing site penalty), skip calling Send() and Sbnd() */
double EndScoreThresh2 = 2.0; /* If pairwise score without Sbnd() is below bestscore (so far) by this threshold (as multiple of missing site penalty), skip calling vSbnd() */

int fast_mode = 0; /* Values >0 turn on various heuristic optimizations */
float fast_mode_tolerance=4; /* Tolerance to use for pair oracle. This is in units of variance computed using sf/sd */

int XmapLen = 0;/* If > 0 : Add query and ref map length columns to .xmap output */

int NoStat = -1;/* If > 0 : Skip computation of statistics (requires ALIGN_COMPRESS) */

double XmapNoQueryOverlap = 0.0;/* If > 0 : suppress alignments that overlap in the query map by this many kb with higher Pvalue alignments of the same query map, processed in descending order of logPV */

double XmapFilterDelta = 0.0, XmapFilterThresh1 = 0.0, XmapFilterThresh2 = 0.0;/* NOTE:  XmapFilterDelta <= XmapFilterThresh2 <= XmapFilterThresh1 */
/* If XmapFilter* > 0.0 : Modify logPV threshold to XmapFilterThresh2 - XmapFilterDelta and then filter out final matchgroups as follows:
   1. Apply -XmapNoQueryOverlap (if != 0) while flagging each ref contig that had matchgroups removed with the highest confidence value for any matchgroup removed.
   2. If a ref contig has multiple (remaining) matchgroups, only keep the best one AND only if its logPV is better than 2nd best logPV by at least XmapFilterDelta, OR confidence  >= XmapFilterThresh1
   3. If a matchgroup had higher confidence matchgroups with the same ref removed (in step 1) only keep if its confidence is better than XmapFilterThresh1, 
      otherwise only keep if confidence is better than XmapFilterThresh2.
*/

int XmapUnique = 0;/* If > 0 : suppress alignments for query maps that have less than this many sites (on query map) distinct from every higher Pvalue alignment of the same query map */

int XmapChim = 0;/* If > 0 : check for chimeric query maps that align to multiple reference maps in .xmap output and have a minimum of this many disjoint sites in each region */
double XmapChimSpacing = 0.0;/* If > 0.0 : Also check for two alignments to same reference map in .xmap output separated by at least this many kb */

int CovNorm = 0;/* If > 0 : normalize coverage in _r.cmap and refined .cmap to reflect average underlying map size */
double CovLambda = 0.0;/* If CovNorm : The mean map size of the underlying exponential distribution : actual mean map size = CovLambda + MinLen */

int align_format = 0;/* 0 : use text format for .align output
			1 : use binary format (except headers) for .align output
			NOTE : .align input file will have a line "# binary format" IF the format is binary */

int CoverageIntRcmap = 0;/* 1 : extend FRAGCOV_FIX >= 1 to output of _r.cmap in output_csites() */

/* HashTable parameters */
int HashGen = 0;
int HashWin = 5;
int HashScore = 3;
double SDmax = 2.2;
double SDrms = 1.2;
double RelErr = 0.05;
double OffsetKB = 3.0;
int NumErrI= 1, NumErrQ= 1;/* maximum number of false/missing sites per alignment window for inserted and and probe map respectively */
int NumResQ= 0;/* maximum number of unresolved sites per alignment window for probe maps */

double HashSF = 0.10, HashSFscale = 0.5;
double HashSD = 0.03, HashSDscale = 0.5;
double HashSR = 0.025, HashSRscale = 0.5;


int ShiftOffset1 = 31;/* see -hashGrouped */
int GroupedScore = 3;/* see -hashGrouped */
int hashMaxCov = 0;/* see MAXHIT_COVERAGE in hash.h */
int OffsetSpread = 2;/* offset spread used when pairwise==0 : See OFFSET_SPREAD in hash.h */

int hashScaleDelta = 1;/* see -hashScaleDelta */

int hashkeys = 0;/* see -hashkeys */

int Hash_Bits = 0; /* If > 0 : overrides HASH_BITS in hash.h */
int MHash_Bits = 0; /* If > 0 : overrides MHASH_BITS in hash.h */
int HHash_Bits = 0; /* If > 0 : overrides HHASH_BITS in hash.h */
int SubsetMaps = 0; /* If > 0 && -ref : Output subset of query maps with HashTable matches in <prefix>_hash.bnx. */
int DeferMerge = 0; /* If > 0 : Output a list of filenames that can be combined into the HashTable binary by Merge Sorting the files (which can then be overlapped with alignment) */
int HashInsertThreads = 1;/* Number of threads to use for multithreaded HashTable insertion */
int HashQueryThreads = 0; /* If > 0 : Reduce number of threads from MaxThreads to this value for hashtable (except during HashTable insertion) */

int HashMultiMatch = 0;/* Keep track of multiple hashtable matches per map pair, with offset differing by at least this value */
int HashMultiMatchMax = 0;/* maximum number of hashtable matches per map pair considered in each orientation (the best scoring matches are used) (0 means no limit)*/
int HashMultiMatchMaxDelta = 0;/* If > 0 : In each orientation for each map pair keep only hashtable matches within this hashtable score of the best hashtable score */

int HashTxt = 0;/* 1 : output text version of HashTable file */
int HashMatrix = 0;/* 1 : output matrix of matching intervals (including self matches) */

int HashBest = -1;/* If >= 0 : For each query map only use hashtable matches within HashBest of the best hashscore */

int HashMaxHits = 0;/* If > 0 : For each query map only use best HashMax hashtable matches : rounds up HashMax to include all matches with lowest score included */
int HashHitsCov = 0;/* 2nd arg of -HashMaxHits */
int HashHitsGenSiz = 0;/* 3rd arg of -HashMaxHits */

int HashGC = 0;/* If > 0 : Reclaim memory of MatchTable every HashGC probe labels, for probe maps with more than this many labels. */
int HashT2 = 0;/* If > 0 : Use 2 threads per probe map (nested parallelism over alignment orientation) */

int ChimQuality = 0;/* output Chimeric junction Quality Score at each site of Consensus after refinement :
		       score ranges from 0 to 100 and low scores mean more likely to be Chimeric junction at that site */
int noQuality = 1; /* default is to _not_ do quality */

/* read in .errbin file to parameters[0] and copy the results to global parameter space */
void input_err(char *filename, Cparameters *parameters)
{
  FILE *fp;
  if((fp = fopen(filename,"rb"))==NULL){
    char *err = strerror(errno);
    printf("input_err: Failed to open input file %s:errno=%d:%s\n",filename, errno, err);
    printf("Current Time= %s\n", currentDateTime());
    sleep(1);// sleep 1 second so this error cannot happen multiple times within the same second
    fflush(stdout);exit(1);
  }

  long pos_cur, pos_end;
  size_t i;
  
  pos_cur = ftell(fp);
  fseek(fp, 0L, SEEK_END);
  pos_end = ftell(fp);
  fseek(fp, pos_cur, SEEK_SET);
  if(pos_end - pos_cur != sizeof(Cparameters)){
    printf("size of %s (%ld bytes) does not match sizeof(Cparameters)=%lu bytes (Cparameters changed since file was created ?)\n",filename, pos_end - pos_cur, sizeof(Cparameters));
    fflush(stdout);

    if(pos_end - pos_cur == sizeof(Cparameters1003)){
      printf("Trying to read as Cparameters version 1003\n");
      fflush(stdout);

      Cparameters1003 p1003;
      i = fread(&p1003,sizeof(Cparameters1003),1,fp);
      if(i != 1){
	char *err = strerror(errno);
	printf("Failed to read %s (i=%lu):errno=%s\n", filename, i, err);
	printf("sizeof(Cparameters1003)=%lu\n",sizeof(Cparameters1003));
	fflush(stdout);exit(1);
      }
      fclose(fp);
      if(p1003.version != 1003){
	printf("Binary file %s was generated with version=%d, expecting %d : Cparameters has changed since file was generated\n",filename,parameters->version,1003);
	fflush(stdout);exit(1);
      }
      if(p1003.MaxResBins != 64+7){
	printf("Binary file %s generated with RESBINS=%d, expecting 64+7\n",filename, p1003.MaxResBins);
	fflush(stdout);exit(1);
      }
      
      /* copy from p1003 to parameters[0] */
      parameters->version = Cparameters_Version;
      parameters->colors = p1003.colors;
      parameters->mapcnt = p1003.mapcnt;
      parameters->ATmapcnt = p1003.ATmapcnt;
      parameters->sitecnt = p1003.sitecnt;
      parameters->sumX = 0.0;
      parameters->totlen = 0.0;
      parameters->MA_mean = p1003.MA_mean;
      parameters->PixelLen = p1003.PixelLen;
      parameters->PixelLenSD = p1003.PixelLenSD;
      parameters->mres = p1003.mres;
      parameters->mresSD = p1003.mresSD;
      parameters->logLR = p1003.logLR;
      parameters->ATlogLR = p1003.ATlogLR;

      parameters->outlierRate = p1003.outlierRate;
      parameters->EndoutlierRate = p1003.EndoutlierRate;
      
      parameters->maxresbias = p1003.maxresbias;
      parameters->MaxResBins = p1003.MaxResBins;
      for(int c = 0; c < MAXCOLOR;c++){
	parameters->ResBins[c] = p1003.ResBins[c];
	for(int B = 0; B <= 64+7;B++){
	  parameters->resbias[c][B] = p1003.resbias[c][B];
	  parameters->resbiasX[c][B] = p1003.resbiasX[c][B];
	}
	
	parameters->res[c] = p1003.res[c];
	parameters->resSD[c] = p1003.resSD[c];
	
	parameters->FP[c] = p1003.FP[c];
	parameters->FPfrac[c] = p1003.FPfrac[c];
	parameters->FN[c] = p1003.FN[c];

	parameters->LabelDensity[c] = p1003.LabelDensity[c];
	parameters->LabelDensityX[c] = 0.0;
	parameters->minSNR[c] = p1003.minSNR[c];
      }
      for(int c = 0; c <= MAXCOLOR;c++){
	parameters->SF[c] = p1003.SF[c];
	parameters->SD[c] = p1003.SD[c];
	parameters->SR[c] = p1003.SR[c];
	parameters->SE[c] = p1003.SE[c];
      }
      goto Lconverted;
    }

    exit(1);
  }

  i = fread(parameters,sizeof(Cparameters),1,fp);
  if(i != 1){
    char *err = strerror(errno);
    printf("Failed to read %s (i=%lu):errno=%s\n", filename, i, err);
    printf("sizeof(Cparameters)=%lu\n",sizeof(Cparameters));
    fflush(stdout);exit(1);
  }
  fclose(fp);

  if(parameters->version != Cparameters_Version && !(parameters->version == 1003 && Cparameters_Version == 1004)){
    printf("Binary file %s was generated with version=%d, expecting %d : Cparameters has changed since file was generated\n",filename,parameters->version,Cparameters_Version);
    fflush(stdout);exit(1);
  }

  if(parameters->MaxResBins != RESBINS){
    printf("Binary file %s generated with RESBINS=%d cannot be read with RESBINS=%d\n",filename, parameters->MaxResBins, RESBINS);
    printf("  mapcnt=%d,aligned=%d\n",parameters->mapcnt,parameters->ATmapcnt);
    printf("  FP=%0.6f, FN=%0.6f, SF= %0.6f, SD=%0.6f, SR=%0.6f, SE=%0.6f\n", parameters->FP[0],parameters->FN[0],parameters->SF[0],parameters->SD[0],parameters->SR[0],parameters->SE[0]);
    printf("  sizeof(Cparameters)=%lu\n", sizeof(Cparameters));
    fflush(stdout);exit(1);
  }

 Lconverted:

  if(colors != parameters->colors){
    if(colors_set){
      if(parameters->colors==2){
	printf("Binary file %s generated with 2 colors data must be used with 2-color molecules and reference (colors=%d)\n",filename,colors);
	fflush(stdout);exit(1);
      } else {
	printf("Binary file %s generated with %d color data cannot be used with 2-color molecules/reference (colors=%d)\n",filename,parameters->colors,colors);
	fflush(stdout);exit(1);
      }
    } 

    if(VERB/* HERE HERE >=2 */){
      printf("input_err(%s):setting colors = %d, colors_set= 1\n",filename,parameters->colors);
      fflush(stdout);
    }

    colors = parameters->colors;
    colors_set = 1;
  }

  /* copy parameters to global space */
  if(colors > 1){
    MA_mean = parameters->MA_mean;
    MA_SF = parameters->SF[colors];
    MA_SD = parameters->SD[colors];
  }
  bpp = parameters->PixelLen;
  if(!(mres_defined && mres_override)){
    mres = parameters->mres;
    mresSD = parameters->mresSD;
  }
  for(int c = 0; c < colors; c++){
    res[c] = parameters->res[c];
    resSD[c] = parameters->resSD[c];

    FP[c] = parameters->FP[c];
    FN[c] = parameters->FN[c];
    SF[c] = parameters->SF[c];
    SD[c] = parameters->SD[c];
    if(QUADRATIC_VARIANCE)
      SR[c] = parameters->SR[c];
    if(RES_VARIANCE)
      SE[c] = parameters->SE[c];

    minSNR[c] = parameters->minSNR[c];
  }

  ResBinsFile = 0;
  for(int c = 0; c < colors; c++){
    int NumBins = parameters->ResBins[c];
    if(NumBins > 0){
      origResBins[c] = max(ResBins[c],NumBins);
      ResBins[c] = NumBins;
      maxresbias = parameters->maxresbias;
      ResBinsFile = 1;

      for(int Bin = 0; Bin <= NumBins; Bin++){
	resbias[c][Bin] = parameters->resbias[c][Bin];
	resbiasX[c][Bin] = parameters->resbiasX[c][Bin];
	if(resbiasX[c][Bin] > maxresbias){
	  printf("ERROR in parameters file %s: c=%d,resbiasX[c][%d]=%0.4f,maxresbias=%0.4f,ResBins[c]=%d/%d,ResBinsFile=%d,parameters->ResBins[c]=%d\n",
		 filename,c,Bin,resbiasX[c][Bin],maxresbias,ResBins[c],origResBins[c],ResBinsFile,parameters->ResBins[c]);
	  exit(1);
	}
	if(Bin == 0 && resbiasX[c][Bin] >= maxresbias){
	  printf("ERROR in parameters file %s:c=%d,resbiasX[c][0]=%0.4f,maxresbias=%0.4f : cannot be same\n",filename,c,resbiasX[c][0],maxresbias);
	  exit(1);
	}
	if(Bin > 0 && resbiasX[c][Bin] <= resbiasX[c][Bin-1]){
	  if(VERB){
	    printf("WARNING: In parameters file %s:c=%d,resbiasX[c][%d]=%0.4f,resbiasX[c][%d]=%0.4f (should be different) : reducing ResBins[c] from %d to %d (maxresbias=%0.4f,origResBins[c]=%d)\n",
		   filename, c, Bin-1, resbiasX[c][Bin-1], Bin, resbiasX[c][Bin], ResBins[c], (resbiasX[c][Bin-1] >= maxresbias) ? Bin-1 : Bin, maxresbias,origResBins[c]);
	    fflush(stdout);
	  }
	  if(resbiasX[c][Bin-1] >= maxresbias)
	    ResBins[c] = Bin-1;
	  else {
	    ResBins[c] = Bin;
	    resbiasX[c][Bin] = maxresbias;
	  }
	  break;
	}
      }
    }
    if(DEBUG) assert(ResBins[c] <= RESBINS);
  }

  if(VERB/* HERE >=2 */){/* display parameters read in */
    printf("Parameters read in from %s:\n",filename);
    if(mres_defined && mres_override)
      printf("\tbpp=%0.10f (unchanged mres=%0.8f mresSD=%0.8f)\n",
	     bpp*1000.0, mres, mresSD);
    else
      printf("\tbpp=%0.10f mres=%0.8f mresSD=%0.8f\n",
	     bpp*1000.0, mres, mresSD);
    for(int c = 0; c < colors; c++){
      printf("\tcolor=%d: res=%0.8f,resSD=%0.8f,FP=%0.8f FN=%0.8f sf=%0.9f sd=%0.10f sr=%0.10f se=%0.9f minSNR=%0.3f\n", c+1, res[c],resSD[c],FP[c], FN[c], SF[c], SD[c], SR[c], SE[c], minSNR[c]);
      if(ResBins[c] > 0){
	printf("\tmaxresbias=%0.10f ResBins=%d/%d:\n", maxresbias, ResBins[c], origResBins[c]);
	if(parameters->ResBins[c] > 0){
	  printf("\t\tBin\tsize\tbias\n");
	  for(int Bin = 0; Bin <= ResBins[c]; Bin++)
	    if(fabs(resbias[c][Bin]) >= 1e-7)
	      printf("\t\t%d\t%0.5f\t%0.7f\n", Bin, resbiasX[c][Bin], resbias[c][Bin]);
	}
      }
    }
    if(colors > 1)
      printf("\tMA_mean=%0.9f MA_sf=%0.9f MA_sd=%0.10f\n",MA_mean,MA_SF,MA_SD);
    fflush(stdout);
  }

  for(int c = 0; c < colors; c++){
    /* save original error parameters before applying bounds */
    origFP[c] = FP[c];
    origFN[c] = FN[c];
    origSF[c] = SF[c];
    origSD[c] = SD[c];
    origSR[c] = SR[c];
    origSE[c] = SE[c];
    origres[c] = res[c];
    origresSD[c] = resSD[c];

    /* adjust initial values based on bounds */
    if(FP[c] < MINFP){
      printf("FP[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,FP[c],MINFP,MAXFP,MINFP);
      FP[c] = MINFP;
    } else if(FP[c] > MAXFP){
      printf("FP[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,FP[c],MINFP,MAXFP,MAXFP);
      FP[c] = MAXFP;
    }
    if(FN[c] < MINFN){
      printf("FN[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,FN[c],MINFN,MAXFN,MINFN);
      FN[c] = MINFN;
    } else if(FN[c] > MAXFN){
      printf("FN[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,FN[c],MINFN,MAXFN,MAXFN);
      FN[c] = MAXFN;
    }
    if(SF[c] < MINSF){
      printf("SF[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SF[c],MINSF,MAXSF,MINSF);
      SF[c] = MINSF;
    } else if(SF[c] > MAXSF){
      printf("SF[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SF[c],MINSF,MAXSF,MAXSF);
      SF[c] = MAXSF;
    }
    if(SD[c] < MINSD){
      printf("SD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SD[c],MINSD,MAXSD,MINSD);
      SD[c] = MINSD;
    } else if(SD[c] > MAXSD){
      printf("SD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SD[c],MINSD,MAXSD,MAXSD);
      SD[c] = MAXSD;
    }
    if(SR[c] < MINSR){
      printf("SR[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SR[c],MINSR,MAXSR,MINSR);
      SR[c] = MINSR;
    } else if(SR[c] > MAXSR){
      printf("SR[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SR[c],MINSR,MAXSR,MAXSR);
      SR[c] = MAXSR;
    }

    if(SE[c] < MINSE){
      printf("SE[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SE[c],MINSE,MAXSE,MINSE);
      SE[c] = MINSE;
    } else if(SE[c] > MAXSE){
      printf("SE[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SE[c],MINSE,MAXSE,MAXSE);
      SE[c] = MAXSE;
    }
    if(res[c] < MinRes){
      printf("res[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,res[c],MinRes,MaxRes,MinRes);
      res[c] = MinRes;
    } else if(res[c] > MaxRes){
      printf("res[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,res[c],MinRes,MaxRes,MaxRes);
      res[c] = MaxRes;
    }
    if(resSD[c] < MinResSD){
      printf("resSD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,resSD[c],MinResSD,MaxResSD,MinResSD);
      resSD[c] = MinResSD;
    } else if(resSD[c] > MaxResSD){
      printf("resSD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,resSD[c],MinResSD,MaxResSD,MaxResSD);
      resSD[c] = MaxResSD;
    }

    if(SR[c] < 0.0){
      printf("ERROR:SD[%d]=%0.6f is not valid (must specify +ve sr)\n",c,SD[c]);
      fflush(stdout);exit(1);
    }
    if(SF[c] < 0.0){
      printf("ERROR:SF[%d]=%0.6f is not valid (must specify +ve sr)\n",c,SF[c]);
      fflush(stdout);exit(1);
    }

    if(SD[c] < -sqrt(2.0*SF[c]*SR[c]) * MINSD_MULT){
      printf("WARNING:SD[%d]=%0.6f is not valid for SF=%0.6f,SR=%0.6f: changing to %0.6f\n",c,SD[c],SF[c],SR[c],-sqrt(2.0*SF[c]*SR[c])*MINSD_MULT);
      fflush(stdout);
      SD[c] = -sqrt(2.0*SF[c]*SR[c]) * MINSD_MULT;
    }
    fflush(stdout);
  }
}

#ifdef WIN32
int output_filename_matched(const char *filename)
{
	return 1;// we need to implement regex matching 
}
#else

/* returns 0 IFF the filename fails to match the output_filter OR matches output_veto_filter : otherwise returns 1 */
int output_filename_matched(const char *filename)
{

  if(output_veto_filter != NULL && output_veto_filter_reg == NULL) {
    output_veto_filter_reg = (regex_t*) calloc(1, sizeof(*output_veto_filter_reg));
    if(regcomp(output_veto_filter_reg, output_veto_filter, REG_EXTENDED|REG_ICASE|REG_NOSUB)) {
      int myerrno = errno;

      fprintf(stderr,"regcomp():Failed to compile regular expression %s:\n%s\n",output_veto_filter,strerror(myerrno));
      fflush(stderr);

      regfree(output_veto_filter_reg);
      free(output_veto_filter_reg);
      output_veto_filter_reg = NULL;
      output_veto_filter = NULL;
    }
  }
  
  if(output_veto_filter_reg != NULL && !regexec(output_veto_filter_reg, filename, 0, NULL, 0)){
    if(VERB>=1+RELEASE){
      printf("filename='%s', output_veto_filter='%s' : output suppressed\n",filename,output_veto_filter);
      fflush(stdout);
    }
    return 0;
  }

  if(output_filter != NULL && output_filter_reg == NULL) {
    output_filter_reg= (regex_t*) calloc(1, sizeof(*output_filter_reg));
    if(regcomp(output_filter_reg, output_filter, REG_EXTENDED|REG_ICASE|REG_NOSUB)) {
      int myerrno = errno;

      fprintf(stderr,"regcomp():Failed to compile regular expression %s:\n%s\n",output_filter,strerror(myerrno));
      fflush(stderr);

      regfree(output_filter_reg);
      free(output_filter_reg);
      output_filter_reg = NULL;
      output_filter = NULL;
    }
  }

  if(output_filter_reg != NULL && regexec(output_filter_reg, filename, 0, NULL, 0)){
    if(VERB>=1+RELEASE && output_filter != NULL){
      printf("filename='%s', output_filter='%s' : output suppressed\n",filename,output_filter);
      fflush(stdout);
    }

    return 0;
  }

  if(VERB>=2 && output_filter != NULL){
    printf("filename='%s', output_filter='%s' : output NOT suppressed\n",filename,output_filter);
    fflush(stdout);
  }
  return 1;
}
#endif

/* swap parameters for color 1 with parameters for color "usecolor" */
void swapparam(int usecolor, Cparameters *p, int iterations)
{
  if(usecolor == 1)
    return;
  int iuse = usecolor - 1;
  double tmp = res[0]; res[0] = res[iuse]; res[iuse] = tmp;
  tmp = resSD[0]; resSD[0] = resSD[iuse]; resSD[iuse] = tmp;
  tmp = FP[0]; FP[0] = FP[iuse]; FP[iuse] = tmp;
  tmp = FN[0]; FN[0] = FN[iuse]; FN[iuse] = tmp;
  tmp = SF[0]; SF[0] = SF[iuse]; SF[iuse] = tmp;
  tmp = SD[0]; SD[0] = SD[iuse]; SD[iuse] = tmp;
  tmp = SR[0]; SR[0] = SR[iuse]; SR[iuse] = tmp;
  tmp = SE[0]; SE[0] = SE[iuse]; SE[iuse] = tmp;
  if(ResBins[iuse] > 0 || ResBins[0] > 0){
    int itmp = ResBins[0]; ResBins[0] = ResBins[iuse]; ResBins[iuse] = itmp;
    int BINmax = max(ResBins[0],ResBins[iuse]);
    for(int Bin = 0; Bin <= BINmax; Bin++){
      tmp = resbiasX[0][Bin]; resbiasX[0][Bin] = resbiasX[iuse][Bin]; resbiasX[iuse][Bin] = tmp;
      tmp = resbias[0][Bin]; resbias[0][Bin] = resbias[iuse][Bin]; resbias[iuse][Bin] = tmp;
    }
  }
  tmp = minSNR[0]; minSNR[0] = minSNR[iuse]; minSNR[iuse] = tmp;

  if(p != NULL)
    for(int i = 0; i < iterations; i++){
      tmp = p[i].res[0]; p[i].res[0] = p[i].res[iuse]; p[i].res[iuse] = tmp;

      tmp = p[i].resSD[0]; p[i].resSD[0] = p[i].resSD[iuse]; p[i].resSD[iuse] = tmp;
      tmp = p[i].FP[0]; p[i].FP[0] = p[i].FP[iuse]; p[i].FP[iuse] = tmp;
      tmp = p[i].FN[0]; p[i].FN[0] = p[i].FN[iuse]; p[i].FN[iuse] = tmp;
      tmp = p[i].SF[0]; p[i].SF[0] = p[i].SF[iuse]; p[i].SF[iuse] = tmp;
      tmp = p[i].SD[0]; p[i].SD[0] = p[i].SD[iuse]; p[i].SD[iuse] = tmp;
      tmp = p[i].SR[0]; p[i].SR[0] = p[i].SR[iuse]; p[i].SR[iuse] = tmp;
      tmp = p[i].SE[0]; p[i].SE[0] = p[i].SE[iuse]; p[i].SE[iuse] = tmp;
      if(p[i].ResBins[iuse] > 0 || p[i].ResBins[0] > 0){
	int itmp = p[i].ResBins[0]; p[i].ResBins[0] = p[i].ResBins[iuse]; p[i].ResBins[iuse] = itmp;
	int BINmax = max(p[i].ResBins[0],p[i].ResBins[iuse]);
	for(int Bin = 0; Bin <= BINmax; Bin++){
	  tmp = p[i].resbiasX[0][Bin]; p[i].resbiasX[0][Bin] = p[i].resbiasX[iuse][Bin]; p[i].resbiasX[iuse][Bin] = tmp;
	  tmp = p[i].resbias[0][Bin]; p[i].resbias[0][Bin] = p[i].resbias[iuse][Bin]; p[i].resbias[iuse][Bin] = tmp;
	}
      }
      tmp = p[i].minSNR[0]; p[i].minSNR[0] = p[i].minSNR[iuse]; p[i].minSNR[iuse] = tmp;
    }

}


