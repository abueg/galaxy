#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <sys/types.h>
#ifdef WIN32
#include <regex>
#define regex_t regex
#else
#include <regex.h>
#endif

#include "globals.h"

static Ident parameters_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/parameters.h 11543 2020-08-26 18:53:06Z tanantharaman $");

extern double MININTERVAL; /* minimum interval size in kb (used for end intervals of input maps and split maps) */

extern const char *tmpfs;/* use this folder (if it exists) for temporary files. Should be a fast filesystem (at least 750 Mb/sec) */

extern char *vfixx_filename[MAXFILES];/* input .vfixx (or .cmap) file name */
extern int num_files;/* number of input files (-i and -ref) */
extern char *spots_filename[MAXFILES+1];/* input .spots file name */
extern int num_reffiles;/* number of input -ref files */
extern char *output_prefix;/* output file prefix (already includes _contigN with N derived from -id, if present) */
extern char *draft_prefix;/* draft output file prefix (_contig<N>_refined to be added for each contig output if -id is specified OR Assembler contig id is known) */
extern char *output_filter; /* regular expression to filter output */
extern regex_t *output_filter_reg; /* regular expression data structure */
extern char *output_veto_filter; /* regular expression to filter output */
extern regex_t *output_veto_filter_reg; /* regular expression data structure */
extern char *stdout_file;/* If != 0 : File to which stdout is redirected */
extern char *stderr_file;/* If != 0 : File to which stderr is redirected */

extern int NFSdelay;/* maximum possible delay in seconds between file fclose() on one NFS client and file becoming visible (due to directory attribute cache) on another NFS client */

extern int BNXVersion;/* BNX version :
			-1 == Unknown Version
			0 == Version 0.1
			1 == Version 1.0 */
extern int BNXMinorVersion;/* Only used with BNXVersion == 1 to distinguish 1.2 from 1.3 */
extern int BNXnoRunData;/* If 1 : Don't output RunData header lines (if present from BNX 1.x input) to any BNX output files (downgrades format to BNX 1.0) */
extern int BNXheaderFXY;/* Molecule headers have "Column StartFOV StartX StartY EndFOV EndX EndY" as additional fields : set on input, all BNX input files must match :
			 -1 : Unknown (no BNX files read in)
			 0 : First BNX file did NOT have these extra columns
			 1 : First BNX file had these extra columns */


extern int BNXheaderChim;/* Molecule headers have 6 extra columns "FeatureMinPos FeatureMinScore FeatureMaxPos FeatureMaxScore FeatureMinScaled FeatureMaxScaled" : set on input, all BNX input files must match :
			  -1 : Unknown (no BNX files read in)
			  0 : One of the previous BNX files did NOT have these extra columns
			  1 : All previous BNX files had these extra columns */

extern int BNXRunHeaderChim;/* RunDdata header have 2 extra columns : set on input, all BNX input files must match :
			  -1 : Unknown (no BNX files read in)
			  0 : one of the previous BNX file did NOT have these extra columns
			  1 : All previous BNX files had these extra columns */

/* NOTE during merge BNXminlen,BNXminSNR are replaced by the larger of : 
   1. The smallest input value in any input BNX being merged.
   2. The -minlen -minSNR command line option (if specified) 
*/
extern double BNXminlen;/* For BNX version 1.x (if != 0.0) : global minimum length (kb) used to generate current BNX input (independent of per-run values in Run Data strings) */
extern double BNXminSNR[MAXCOLOR];/* for BNX version 1.x (if != 0.0) : global minimum SNR threshold value applied to labels (independent of per-run values in Run Data strings) */

extern int MapType;/* Type of map read in from -i input files :
		      -1 : unknown (use filename suffix to infer : .vfixx or .bnx == Molecule map, .cmap = Assembly/Contig map)
		      0 : Molecule map
		      1 : Assembly/Contig map */

extern int MapTrueSite;/* 1 if input .bnx has truesite information (one extra line per color, with two values per map site) */
extern int MapSNR; /* 1 if input .bnx has SNR information (one extra line per color after map sites lines) */
extern int MapIntensity;/* 1 if input .bnx has Intensity information (one extra line per color after map sites lines) */
extern int MapStitched;/* 1 if input .bnx has stitch information (one extra line per color after map sites lines) */
extern int MapStitchLoc;/* 1 if input .bnx has stitch Location information (one extra line after map sites lines) */
extern int MapPSFWidth;/* 1 if input .bnx has PSFWidth site information (one extra line per color after map sites lines) */
extern int MapImageLoc;

extern int MapSNROrig,MapIntensityOrig,MapStitchedOrig,MapStitchLocOrig,MapPSFWidthOrig,MapImageLocOrig;/* original command line values (1 iff corresponding information is Required for all input bnx maps) */

extern int Tracksites;/* Use TrueSite[] values to track original map index as map is transformed by -mres (or -rres) etc */
extern char *QueryOrigin;/* If Tracksite != 0 : original input map name (or names separated by commas) for query maps */
extern char *RefOrigin;/* If Tracksite != 0 : original input map name (or names separated by commas) for reference maps */
extern int RefTrueSite;/* > 0 : RefOrigin and corresponding truesite information read in from one or more -ref input files 
		       NOTE : MapTrueSite is set of if truesite information is read in from one or more -i input files */


extern int QXfilter;/* If 1 : filter out QX codes other than those specified on command line (see MapSNROrig etc) */
extern int QXerror;/* If 1 : exit with error message if different BNX files have different sets of QX codes (otherwise just discard QX codes not present in all BNX files) */
extern int QXmismatch;/* If 1, some files have different QX codes than 1st file */

extern int CmapSNR;/* If 1 : output complete SNR distribution for consensus CMAPs (in addition to Gmean and logSNRsd). Requires MapSNR==1 */

extern int CmapChimQuality;/* 1 : Some Cmap input had ChimQuality data (0 means none, -1 means unknown : no input Cmap files read) */
extern int CmapMask;/* 1 : Some Cmap input had Mask data (0 means none, -1 means unknown : no input Cmap files read) */

extern char *BackgroundSNR;/* If != 0 : This is the name of the file with the background SNR distribution and is to be used to update Pvalue computations when both maps have SNR distributions (must be consensus maps) */
extern int bgSNRcnt[MAXCOLOR];
extern double *bgSNR[MAXCOLOR];/* If BackgroundSNR != 0 : bgSNR[c=0..colors-1][0..bgSNRcnt-1] is the background SNR distribution for color c read in from the file named BackgroundSNR */

extern double bpp;/* command line value */
extern double BppStepsiz;/* bpp initial value scan  step size */
extern int BppSteps;/* bpp initial value scan number of steps (in each direction from the original value) */
extern int BppMincnt;/* bpp initial value scan is triggered when fewer than this many molecules align */
extern int NoBpp;/* Do not apply bpp value to rescale -i input files (see -NoBpp) */
extern double PixelLen;/* bpp value read from vfixx file (If != 0) and modified subsequently by command line bpp value and -M (from updated bpp estimates) */

extern int bnxerror;/* If map from -i input .bnx file has sites beyond map ends, exit with an error message (instead of printing a WARNING and extending map ends as with .cmap) */
extern int bnxmergeerror; /* If map from multiple -i input BNX 1.0 files has missing or duplicate Run Data lines exist with error message (otherwise just print warning message) */
extern int bnxBigid; /* If 1 : Use 19 digit id numbers, if possible (if fields overflow, renumber all ids sequentially starting at 1)
		     If 0 : Avoid changing id numbers (if duplicates are present, renumber all ids sequentially starting at 1) */


extern double EndTrimCov;/* trim ends of draft consensus with coverage below this level */
extern double EndTrimFrac;/* If > 0.0 : Estimate length of unlabeled ends based on this percentile of distance from last aligned label to molecule end or next label */
extern double EndTrimCnt;/* If > 0.0 : Mark end as non-extendible if at least this many molecules agreed (within sizing error) on length of unlabeled ends */
extern double EndTrimWt;/* If > 0.0 : Use this instead if RefineMinMapwt to decide which molecules to count for EndTrimFrac and EndTrimCnt */

extern double EndTrimCov2;/* Further trim back ends of draft consensus until no Internal regions with coverage below this level remains, but not more than EndTrimLen2 kb beyond original ends
			    (typically EndTrimCov2 <= EndTrimCov) */
extern double EndTrimLen2;

extern double EndTrimNoRefine; /* If > 0.0 : Don't refine consensus more than this distance beyond initial estimate of where coverage drops below EndTrimCov */
extern int EndTrimCQ;/* see -EndTrimCQ */

extern double EndTrimRatio; /* If > 0.0 : After trimming ends based on EndTrimCov (and before applying EndTrimFrac) keep trimming back the last label (at either end) if the its label coutn is <= EndTrimRatio times the label cnt of the next label */

extern double res[MAXCOLOR]; /* Query resolution in pixels (PixelLen) */
extern double resSD[MAXCOLOR];/* Query resolution SD in pixels (PixelLen) : reference sites at a distance of X are assumed to resolve with probabily 0.5(1+erf((X-res)/(resSD*sqrt(2)))).
		       If resSD = 0, sites closer than "res" are always merged and other sites are never merged. */

extern double FP[MAXCOLOR];/* False Positive Density per 100kb (one per color) */
extern double FN[MAXCOLOR];/* False Negative probability (one per color) */
extern double R[MAXCOLOR];/* Misaligned sites Density per 100kb (one per color) (for pairwise alignment only) */
extern double SF[MAXCOLOR+1];/* scaling error sigma in kb (one per color) */
extern double SD[MAXCOLOR+1];/* scaling error sigma in root-kb (one per color) */
extern double SR[MAXCOLOR+1];/* scaling error sigma as a ratio (one per color) */
extern double SE[MAXCOLOR+1];/* scaling error for resolution based variance */
extern double MA_mean;/* mean misalignment value (drift bias) between 2 colors */
extern double MA_SF,MA_SD;/* MA_SF and MA_SD are for misalignment between 2 colors */

extern double MinFP;/* see MINFP in constants.h : can be modified with -MinFP */
extern double MaxFP;/* see MAXFP in constants.h : can be modified with -MaxFP */

extern double MinFN;/* see MINFN in costants.h : can be modified with -MinFN */
extern double MaxFN;/* see MAXFN in constants.h : can be modified with -MaxFN */

extern double MinSF;/* see MINSF in constants.h : can be modified with -MinSF */
extern double MaxSF;/* see MAXSF in constants.h : can be modified with -MaxSF */

extern double MinSD;/* See MINSD in constants.h : can be modified with -MinSD */
extern double MaxSD;/* See MAXSD in constants.h : can be modified with -MaxSD */
extern double MinSD_Mult;/* see MINSD_MULT in constants.h : can be molified with -MinSD_Mult */

extern double MinSR;/* see MINSR in constants.h : can be modified with -MinSR */
extern double MaxSR;/* see MAXSR in constants.h : can be modified with -MaxSR */

extern double MinSE;/* see MINSR in constants.h : can be modified with -MinSR */
extern double MaxSE;/* see MAXSE in constants.h : can be modified with -MaxSE */

extern double MinRes;/* minimum value of -res parameter */
extern double MaxRes;/* maximum value of -res parameter */

extern double MinResSD;/* minimum value of -resSD parameter */
extern double MaxResSD;/* maximum value of -resSD parameter */

extern double maxresbias;/* If != 0, intervals below this value are corrected for -ve bias caused by limited resolution. Uses piecewise Linear regression of bias vs log(x) with RESBINS bins */ 
extern int ResBins[MAXCOLOR],origResBins[MAXCOLOR];/* number of piecewise linear segments to use */
extern int ResBinsFile;/* If != 0 : ResBins and maxresbias were read with -readparameters */
extern double resbias[MAXCOLOR][RESBINS+1];/* resbias[i=0..ResBins-1] is the bias (in KB) of raw intervals x at  x = resbiasX[i] */
extern double resbiasX[MAXCOLOR][RESBINS+1];/* resbiasX[ResBins] == maxresbias and resbiasX[0] == mresKB and the other knots in the piecewise Linear model are selected to be seperated by similar number of sample values */

extern int RA_MIN_MEM;/* value of MIN_MEM in refalign.cpp */
extern int RA_MIN_TIME;/* see -RAmem in RefAligner.cpp help message */
extern int QP_MIN_MEM;/* value of MIN_MEM in qprobeval.cpp, set by -QPmem */
extern int MP_MIN_MEM;/* value of MIN_MEM in mprobeval.cpp, set by -MPmem */
extern int MP_MIN_TIME;/* see -MPmem in RefAligner.cpp help message */
extern double RA_MAX_MEM;/* free memory every RA_MIN_TIME, but only if memory use exceeds this fraction of the threshold at which memory is freed at every oppertunity */

extern int RAscoremem;/* If 1, deallocate score table memory before each refine() call */

extern int ResEstimate;/* If != 0 : Estimate res and resSD */

extern int minSNRestimate;/* If != 0 : Estimate minSNR[] in range SNRmin..SNRmax */
extern double SNRmin,SNRmax,SNRstep,SNRminratio;
extern int SNRmaxmaps;/* If != 0 AND SNRstep > 0.0 : subset number of maps to no mroe than SNRmaxmaps which performing minSNR scan */

extern int outlierRate;/* output outlier and Endoutlier rate as 2 additional columns in .err file */
extern int perLabelScore;/* output score per label in .err file */
extern int LabelDensityType;/* see -LabelDensityType */

extern int ScanCorrection;/* If > 0 : Apply per Scan scaling correction to all -i input maps. Value is number of gradient updates per -M iteration (Only applied with BNX 1.0 if UniqueScans > 1) */
extern int ScanCorrectMinMaps;/* Minimum number of aligned maps in a single cohort/scan to update the scaling factor */

extern double ScanCorrectMinSNR;/* If >= 0.0 : -minSNR value applied just before _rescaled.bnx output */
extern double ScanCorrectMinLen;/* If >= 0.0 : -minlen value applied just before _rescaled.bnx output */

extern int ScanFilter;/* see -ScanFilter */
extern double ScanFilterRm;
extern double ScanFilterRb;
extern double ScanFilterRe;

extern double MaxCoverageAN;/* If > 0.0 : reduce coverage to no more than this amount before _rescaled.bnx output (relative to -ref size) */
extern double MaxCoverageLenAN;/* If > BNXminlen : reduce coverage by first increasing BNXminlen up to MaxCoverageLen, then (if needed) randomly subsampling the input */

typedef struct {
  double scale;/* scale factor of current giter */
  double cumscale;/* cumulative product from all giter */
  int RunIndex;/* Run index into RunDataList[] and RunDataSort[] */
  int ScanNumber;/* ScanNumber for all scans with the same RunIndex (in time order) */

  int totmols;/* total number of molecules in this scan (that were sampled) */
  double totlen;/* total length of molecules in this scan */

  int mols;/* molecules that mapped in this scan */
  double len;/* total length of molecules that mapped in this scan */
  double mappedlen;/* total length of the part of the molecules that mapped in this scan */

  double scaleSD;  /* scaling SD (fractional standard deviation) */

} CscanCorrection;
extern CscanCorrection *ScanScale[MAXCOLOR];/* Scaling factor applied to different scans : typically the weighted average scaling value will be close to 1.0 since this is applied in addition to global -bpp correction
					       Indexed as ScanScale[c][0..UniqueScans-1] : NOTE these parameters are NOT saved in .errbin or read in with -readparameters. Instead the rescaled .bnx is saved as <prefix>_rescaled.bnx
					       (without the -bpp scaling) */
extern struct CRunDataSort *RunDataSort;/* information about Runs from BNX 1.0 input */

extern char *parametersfile;/* filename (.errbin) from which to read error parameters include bias parameters */

extern int floatcov;/* output Coverage and Occrance columsn in .cmap as float rather than int */

extern int usecolor;/* If > 0 : all input data must have 2 colors and only the specified color will be kept and the other color discarded */
extern int hashcolor;/* If > 0 : If input data has >= 2 colors AND usecolor==0, then only the specified color is used by the hashtable (Currently required if colors >= 2 && usecolor==0) */
extern int startcolors;/* original value of colors */
extern int refcolors;/* number of colors in -ref inputs (may be 2, even if colors=1,  in which case 2-color .err files will be output if -usecolor was used) */

extern double mres;/* map resolution in pixels : sites closer than this in input are merged */
extern int mres_defined; /* If > 0 : -mres was specified on commandline (and overrides value from -readparameters) */
extern int mres_override;/* If > 1 : -mres -mresSD on commandline overrides value from -readparameters */
extern double mresSD;/* Standard deviation for mres : If mresSD > 0, sites at a distance of X are resolved with probability 0.5(1+erf((X-mres)/(mresSD*sqrt(2))).
			If mresSD = 0, sites close than "mres" are always merged and other sites are never merged. */
extern double SDrange; /* assume likelihood of gaussian more than this many SD from mean is negligible */
extern double HSDrange;/* value of SDrange used in hash.cpp : usually between 1.0 and SDrange to improve performance of hashtable */
extern double HResMul;/* As alternative to HSDrange * resSD use HRESmul * res (whichever is larger) : only used in hashtable */
extern double rres;/* If > 0.0 : reference maps are de-resed at rres * 0.500
		       If <= 0.0 : reference maps are de-resed at min(0.0, res - SDrange * resSD)*PixelLen */
extern double cres;/* If > rres : Check refined consensus maps : if 2 or 3 sites closer than this many kb have fewer than cresN molecules then merge the site pair that is closest and repeat */

extern double ZERO_MINKB;/* If != 0.0 : minimum distance between labels enforced during refinement (otherwise minimum distance enforced is rres * 0.500 kb) */

extern int cresN;
extern double cresF;
extern int cresFix;// When merging 2 labels with -cres, check if deleting either label scores better than merging them
extern double cresMaxLPdrop;/* If >= 0.0 : Don't merge sites with -cres if LP drops by more than cresMaxLPdrop (see -cresMaxLPdrop) */

extern double cresHapL;// see -cresHapMaxLPdrop (1st arg)
extern double cresHapMaxLPdrop;// see -cresHapMaxLPdrop (2nd arg)

extern int rresDefer;// When Haplotyping defer applying -rres until after Haplotyping (this is always done for -cres) 
extern int refineDefer;// When Haplotyping defer refinement and proceed directly to Haplotyping

extern double cresRefine;/* If 1 : Refine sizes after each cres merging stage (otherwise only Refine sizes once at end of last cres merging stage) */

extern int AlignScore;/* IFF 1 output detailed alignment scores in .align */

extern double Kbias;/* alignment error bias factor K */
extern double KLbias;/* alignment error bias factor for length scaled bias */
extern double KFbias;/* alignment error bias factor for per fragment (aligned interval) bias */
extern double Gbias; /* scale misaligned cuts penalty by this heuristic factor (typically Gbias >= 1) */
extern double Beta;/* Nyugen's Beta parameter */

extern double KbiasRef;
extern double KLbiasRef;
extern double KFbiasRef;

/* input file based XmapCount, XmapSites and XmapLength (valid if XmapCount > 0). Read in from a file using -XmapStatRead. Is output to a file (based on current bnx mapset) by using option -XmapStatWrite */
extern int XmapCount[MAXCOLOR];
extern long long XmapSites[MAXCOLOR];
extern double XmapLength[MAXCOLOR],XmapGtheta[MAXCOLOR];
extern char *XmapStatRead,*XmapStatWrite;/* If != 0: input/output XmapCount,XmapSits and XmapLength from a file based on command line options -XmapStatRead and -XmapStatWrite */

extern double Ch;/* chimerism probability */
extern double ChFdr;/* chimerism False Discovery Rate */

extern double Poutlier;/* outlier probability per true positive site */
extern double PoutlierS;/* outlier probabbility per true positive site for intervals with a stitch (if known) */
extern double PoutlierEnd;/* outlier probability per true positive site for alignment ends */
extern int PoutlierEnd_init;/* If PoutlierEnd was specified on the command line */

extern int EndOutlierPV;/* If 1, correct pvalue for number of end outliers in alignment */

extern int isLinear;/* if genome is linear */

extern double ScoreThreshold;/* Score Threshold : only alignments with score above this threshold are output */
extern double LogPvThreshold;/* Pvalue Threshold, expressed on the -log10 scale */
extern int AlignedSiteThreshold;/* minimum aligned sites for chimeric/local alignments  */
extern double AlignedLengthThreshold;/* minimum aligned length */
extern int AlignedEndOutlierThreshold;/* maximum number of End Outliers */
extern double AlignedOutlierThreshold;/* If >= 0.0 : maximum internal outlier size */
extern int AlignedOutlierLabels;/* internal outliers with at least this many misaligned labels will cause alignments to be filtered out */

extern double AlignedFractionThreshold;/* see -F option */
extern double AlignedLabelDensityThreshold;/* If > 0.0 : see -D option */

extern int RefineEndOutlierThreshold;/* maximum number of End Outliers allowed during refinement (otherwise weight of map is reduced to RefineEndOutlierWt except during final UpdateMap) */
extern double RefineEndOutlierWt;
extern double RefineEndOutlierLen;/* Only count End Outliers larger than this size */
extern int RefineEndOutlierN;/* Only count End Outliers with this many labels (in both ref & qry) */
extern int RefineEndOutlierType;/* If 1 (UNLESS -refine 3 -extonly) weights are NOT restored when computing -Haplotype, -Chimquality or -contigsplit : useful for debugging -extsplit */

extern double RefineMinMapwt;/* minimum map weight to output in refined .xmap */
extern double UnrefineMinMapwt; /* minimum map weight to output in regular (unrefined) .xmap (based on -BestRefWT only) */

extern int FilterXmap;/* see -FilterXmap */
extern int FilterXmapWT;/* see -FilterXmapWT */
extern int FilterXmapSplit;/* see -FilterXmapSplit */

extern int RefineOverwrite;/* Overwrite pre-existing refined contig */

extern double ScoreThreshold2;/* Score Threshold for 2nd best alignment (see -SecondBest) */
extern double LogPvThreshold2;/* Pvalue Thhresold for 2nd best alignmen t, expresssed on the -log10 scale (see -SecondBest) */
extern int AlignedSiteThreshold2;
extern double AlignedLengthThreshold2;
extern int AlignedEndOutlierThreshold2;
extern double AlignedOutlierThreshold2;
extern int AlignedOutlierLabels2;

extern int MinDupSites;/* Minimum number of label in (non-inverted) Duplication : must be at least 1 */

extern int CutFlip;/* If > 0 : Try to align each internal outlier of any primary alignment with this many site intervals in reverse orientation using the same query and reference region expanded
		   by at most CutFlipExpand labels. The Flipped alignment must satisfy at least the 2nd set of thresholds (LogPvThreshold2 etc) */
extern double CutFlipMaxErr;/* maximum relative sizing error of internal outlier checked for inversions */
extern int CutFlipExpand;/* Max number of labels by which to expand outlier intervals during inversion */
extern int CutFlipEndExpand;/* If > 0 : Also apply CutFlip to EndOutliers and expand reference in End (non-aligned) direction by this many labels (typically CutFlipEndExpand > CutFlipExpand) */
extern int CutFlipDupExpand;/* If > 0 : Also apply CutFlip to find Inverted Duplications and expand reference in either neighborhood by this man labels larger than outlier (insertion) query interval */
extern double CutFlipMinErr;/* Minimum value of (X - Y)/(X+Y) for Insertions that are checked for Inverted Duplications */
extern int CutFlipEndFilter;/* If > 0 : Suppress CutFlip Endoutlier inversions if there are already matchgroups that extend into the endoutlier region of either qry or ref by at least "CutFlipEndFlilter" Labels */

extern double InversionNormalLogPv;/* If > 0 : suppress Inversion if normal matchgroups LogPv is worse than this value */
extern double InversionNormalLogPv5;/* If > InversionNormalLogPv : suppress Inversion with 5:1 inversion/non-inversion MG size ratio if non-inv MG LogPv is worse than this value */

extern double ScoreThreshold3;
extern double LogPvThreshold3;
extern double ScoreThreshold4;
extern double LogPvThreshold4;

extern int CutFlipBC;/* Apply local Bonferroni correction to Pvalue for CutFlip inversion matchgroups : this mainly affects confidence of end outlier CutFlip inversions based on 
			difference between labels in endoutlier vs labels in inversion matchgroup */

extern int CutFlipSkip;/* If > 0 : skip checking for CutFlip inversions in outliers with at least CutFlipSkip labels (typically > RefSplit), assuming these inversion matchgroups have already been found */
extern double CutFlipMerge;/* If > 0.0 : Merge outliers seperated by aligned region with score less then -log(CutFlipMerge), before applying CutFlip to the merged region (If <= 0, use Rsplit instead) */

extern double CutFlipMaxOverlap;/* Maximum number of labels overlap between matchgroups for a CutFlip Inversion (as a fraction of aligned labels in inversion matchgroup) */
extern double CutFlipMaxMiss;/* If > 0 : Maximum number of misaligned labels in each breakpoint gap for CutFlip Inversions (as a fraction of aligned labels in inversion matchgroup) */
extern double CutFlipOverlapFilter;/* If > 0 : Do partial filtering of overlapped matchgroups based on this overlap ratio & CutFlipFilterConf, before doing 2-matchgroup inversion (CutFlip) calls */
extern double CutFlipFilterConf;/* see CutFlipOverlapFilter */

extern double InvOverlapFilterQryRatio;/* If > 0 : suppress filtering of opposite orientation overlapped matchgroups (with overlap exceeding OverlapFilterQryRatio),
					  unless overlap in query exceeds InvOverlapFilterQryRatio (in kb) AND InvOverlapFilterQryLabRatio (in Labels) AND confidence of larger matchgroup exceeds 
					  that of smaller matchgroup by InvOverlapFilterConfDelta. The saved smaller matchgroup can only be used to call inversions */
extern double InvOverlapFilterConfDelta;
extern double InvOverlapFilterQryLabRatio;

extern double OverlapFilterQryRatio;/* If > 0 : Filter smaller matchgroup overlapped in query by larger matchgroup by OverlapFilterQryRatio (kb ratio) AND OverlapFilterQryLabRatio (Label ratio) 
				      provided larger matchgroup Confidence exceeds smaller matchgroup confidence by at least OverlapFilterConfDelta */
extern double OverlapFilterQryLabRatio;
extern double OverlapFilterConfDelta;

extern int OverlapFilterSameRef; /* If > 0 : Don't apply OverlapFilter unless ref contigs for the overlapped MGs are the same */

extern double MergeMG;/* If > 0 : After running -OverlapFilter, try to merge neighboring matchgroups if they are the closest matchgroups with no overlap and no other matchgroups overlap the gap region */

extern double SmallDupMaxOverlap;/* Max Overlap (fraction of aligned labels) of small matchgroup with non-outlier region of query in larger matchgroup (same orientation), for small Duplication calls */
extern double SmallDupMinOverlap;/* Min Overlap (fraction of aligned labels) of small matchgroup with non-outlier region of query in larger matchgroup (same orientation), for small Duplication calls */
extern double SmallDupMaxMiss;/* Max fraction of aligned labels of small matchgroup that do not match an aligned label in the reference of the larger matchgroup (same orientatation), for small Duplication calls */
extern double SmallDupConfBC;/* Minimum confidence of small matchgroup after Bonferroni correction by duplicate region gap in labels */

extern double LogPvThresholdTE;/* Pvalue Threshold corresponding to -TE option for extensions */
extern int AlignedSiteThresholdE;/* see -AE */
extern double AlignedLengthThresholdE;/* see -LE */

extern double LogPvThresholdTS;/* Log10 Pvalue Threshold correspondign to -TS option for E&S contigs */
extern int AlignedSiteThresholdS;/* see -AS */
extern double AlignedLengthThresholdS;/* see -LS */

extern int TEinit;/* If -TE option was specified */
extern int AEinit;/* If -AE option was specified */
extern int LEinit;/* If -LE option was specified */

extern int TSinit;/* If -TS option was specified */
extern int ASinit;/* If -AS option was specified */
extern int LSinit;/* If -LS option was specified */

extern int TBcnt;
extern double LogPvThresholdTB[16];/* LogPvThresholdTB[0..TBcnt-1] : Pvalue Threshold(s) corresponding to -TB option for extension maps after preliminary refinement */
extern double TBmult;/* If != 0 : Replace LogPvThreadholdTB[0..TBcnt-1] by a sequence starting at LogPvThresholdTB[0] and decreasing in multiple of TBmult until LogPvThresholdTB[TBcnt-1] is reached.
			If the last value is not equal to LogPvThresholdTB[TBcnt-1], that is added to the sequence */
extern int TBinit;/* If -TB option was specified */


extern double MinLen;/* Minimum -i input map length in kb */
extern double MaxLen;/* Maximum -i input map length in kb */
extern int PostAdjust;/* Apply Minlen/Maxlen after doing -bpp adjustment */

extern int MinSites;/* Minimum -i input map sites */
extern int MaxSites;/* If > 0 : delete maps with more than this many sites */
extern int MaxSites2;/* If > 0 : delete maps with more than this many sites just before output of _rescaled.bnx in refalign.cpp */
extern double MaxSiteDensity;/* If > 0 : delete maps with more than this many sites / 100kb */
extern double MinSiteDensity;/* If > 0 : delete maps with less than this many sites / 100kb */

extern double MaxEndSiteDensity;/* If > 0 : delete maps with more than this many sites / 100kb in either Map End of at least MaxEndSiteDensityLength kb */
extern double MaxEndSiteDensityLength;

extern double MaxEnd;/* If >= 0 : Trim (reduce) unlabeled molecule ends to no more than MaxEnd kb */
extern double MaxInterval;/* If >=0 : Split molecules with internal intervals larger than MaxInterval kb */
extern int refMaxInterval;/* If > 0 : Also apply MaxInterval to reference CMAP input, otherwise only to query BNX input (and query CMAP input if there is NO reference)
			     If applied to ref OR query CMAP input, the maps may be split into multiple maps. With BNX input only the largest piece is retained, so number of maps remain unchanged */

extern double MaxIntensity;/* If > 0 : Maximum AvgIntensity value of molecule backbone */
extern double minSNR[MAXCOLOR];/* minSNR[c] is the minimum SNR filter for color c (c = 0 or 1) (see -minSNR) */
extern double maxPSFWidth[MAXCOLOR];/* maxPSFWidth[c] is the maximum PSFWidth for color c (c = 0 or 1) */

extern int TotalThreads;/* If TotalThreads > 0 and TotalSlots = 0: Total Threads of all jobs on machine : scale MaxMem and HMaxMem by MaxThreads/TotalThreads (TotalThreads = 1 means OMP_NUM_THREADS) */
extern int TotalSlots;/* If TotalThreads > 0 && TotalSlots > 0 : scale MaxMem and HmaxMem by MaxThreads/(TotalThreads*OMP_NUM_THREADS/TotalSlots) */
extern int MaxThreads;/* limit number of threads (if > 0) */
extern int RefineThreads;/* limit number of threads during refinement (-1 means : not specified on command line and defaults to same as MaxThreads) */
extern int BoostThreads;/* If > 0 : Increase MaxThreads to this number after elapsed time exceeds BoostThreadsTime. Do NOT re-adjust maxmem,hashmaxmem due to TotalThreads */
extern int BoostThreadsTime;

extern double MaxMem;/* If > 0 : Maximum memory in Gigabytes */
extern double HMaxMem;/* If > 0 : Maximum memory in Gigabytes for hashtable (defaults to MaxMem) */
extern double MaxVirtMem;/* If > 0 : Maximum Virtual memory in Gigabytes */
extern double MaxMemIncrease;/* If > 0 : Increase MaxMem to at least MemSiz - MaxMemIncrease, where MemSiz is system memory in Gigabytes */

extern int colors;/* number of colors in input files (typically based on -i inputs after processing -usecolor) */
extern int colors_set;/* > 0 IFF number of colors has been specified by the header of the first input file OR the first input file has been successfully read in confirming the value of -colors parameter.
		      All subsequent input files must agree with this value */

extern char *Nickase[MAXCOLOR];/* color names based on first input query (or ref) file read in : only used in input_vfixx.cpp */

extern int DELTA_Y; /* maximum span of aligned sites in reference */
extern int DELTA_Y2;/* If > 0 : smaller value of DELTA_Y to use in all but the last -M iterations of refalign (to achieve greater speed) */
extern int DELTA_X; /* maximum span of aligned sites in Nanomap */
extern int DELTA_X2;/* If > 0 : smaller value of DELTA_X to use in all but the last -M iterations of refalign (to achieve greater speed) */
extern int DELTA;  /* maximum span of aligned sites for pairwise alignment */
extern int KMAX; /* maximum number of consecutive references sites to merge due to resolution limit is KMAX+1 */

extern int deltaXRef;/* If > 0 : replaces DELTA_X during refinement */
extern int deltaYRef;/* If > 0 : replaces DELTA_Y during refinement */
extern int deltaExtXRef;/* If > DELTA_X, allow X (Query) intervals to be this large as long as Y (Ref) intervals are limited to DELTA_Y */
extern int deltaExtYRef;/* If > DELTA_Y, allow Y (Ref) intervals to be this large as long as X (Query) intervals are limited to DELTA_X */

extern int deltaEndRef;/* If > 0 : allow non-outlier ends with up to deltaExtXRef + U - 1 qry labels (or deltaExtYRef + U - 1 ref labels), instead of just deltaXRef,deltaYRef labels */

extern int outlierExtend;/* If > 0 : increase DELTA_X and DELTA_Y as needed to try enlarge internal outliers up to this many sites beyond the current outlier boundaries */
extern int outlierExtendLimit;/* If > 0 : limit increases in DELTA, DELTA_X or DELTA_Y to never exceed this value */
extern int outlierBC;/* If > 0 : Apply BonFerroni correction to to outlier penalty */
extern RFLOAT outlierMax;/* If > 0.0 : Do NOT allow outliers with delta |X-Y| that exceeds this value in kb */
extern double outlierLambda;/* If OUTLIER_LTYPE==0 : outlier likelihood is proportional to exp(-(deltaX+deltaY)/outlierLambda))
			       If OUTLIER_LTYPE==1 : outlier likelihood is proportional to exp(-|deltaX-deltaY|/outlierLambda)) */
extern double outlierLambdaRef;/* If > 0 : Replaces outlierLambda during refinement */

extern RFLOAT OutlierType0;/* see OUTLIER_TYPE in RGentigRefScore.h */
extern RFLOAT OutlierType1;/* see OUTLIER_TYPE1 in RGentigRefScore.h */
extern int OutlierType;/* see OUTLIER_TYPE in RGentigScore.h */

extern int outlierNorm0;/* see -outlierNorm0 (or LAMBDA_FIX0 in RGentigRefScore.h) */
extern int outlierNorm1;/* see -outlierNorm1 (or LAMBDA_FIX1 in RGentigRefScore.h) */
extern int outlierNormRef;/* see -outlierNormRef (or LAMBDA_FIX in RGentigScore.h) */

extern double thetaScale;/* see -thetaScale */
extern double thetaScaleRef;/* see -thetaScaleRef */

extern double RefineStep1;/* Initial iteration sizing increment before UpdateMap() is called (except during RefineA). Smaller values produce more accurate consensus maps but longer run times */
extern double RefineStep2;/* Initial iteration sizing increment if RefineA OR after UpdateMap() is called */
extern double RefineStep3;/* Subsequent iteration sizing increment (when convergence failure is detected) Before UpdateMap() is called (except during refineA) */
extern double RefineStep4;/* Subsequent iteration sizing increment (when convergence failure is detected) if refineA OR after UpdateMap() is called */
extern double RefineStep5;/* Absolute minimum sizing increment (in kb) Before UpdateMap() is called (except during refineA) */
extern double RefineStep6;/* Absolute minimum sizing increment (in kb) if refineA OR after UpdateMap() is called */
extern int RefineStepCnt;/* max number of steps to use with RefineStep3 or RefineStep4 (the number of steps is NOT limited with RefineStep1 or RefineStep2) */
extern double RefineStep7;/* minimum relative sizing change in a neighboring interval that may trigger RefineStepCnt x 2 extra steps */
extern int RefineFix;/* Turns on adjustments and bug fixes in Refinment introduced between revision 4535 and 4607 */

extern double RefineMinRange1;/* see -RefineSizeRange */
extern double RefineMaxRange1;
extern double RefineMinRange2;
extern double RefineMaxRange2;

extern int SiteMerge;/* see -SiteMerge */

extern int IndelScore;/* If 1 : modify score parameters -outlierLambda, -outlierType, -outlier, -outlierEnd, -LRbias during Final Haplotype Indel Scoring (also used to delete them, if score is -ve) */

extern double outlierLambdaRef;/* If > 0 : Replaces outlierLambda during refinement */
extern int outlierLambdaSwitch;/* force outlierLambdaRef to outlierLambdaLabel during refinement when adding/deleting labels (or Haplotype SNPs)
                               force outlierLambdaRef to outlierLambdaSize during refinement when adjusting sizing (or Haplotype Indels) */
extern double outlierLambdaLabel;/* 0.0 means initial value of outlierLambdaRef */
extern double outlierLambdaSize;/* 0.0 means initial value of outlierLambdaRef */
extern double IndelDelOutlierLambda;/* 0.0 means no change during Indel Deletion */

extern int OutlierTypeSwitch;/* force OutlierType to OutliertypeLabel during refinement when adding/deleting labels (or Haplotype SNPs) 
			     force OutlierType to OutliertypeSize during refinement when adjusting sizing (or Haplotype Indels) */
extern int OutlierTypeLabel; 
extern int OutlierTypeSize;
extern int IndelDelOutlierType;/* -1 means no change during Indel Deletion */

extern int PoutlierSwitch;/* force Poutlier to PoutlierLabel during refinement when adding/deleting labels (or Haplotype SNPs)
			 force Poutlier to PoutlierSize during refinement when adjusting sizing (or Haplotype Indels) */
extern double PoutlierLabel; 
extern double PoutlierSize;
extern double IndelDelPoutlier;/* < 0.0 means no change during Indel Deletion */

extern int PoutlierEndSwitch;/* force PoutlierEnd to PoutlierEndLabel during refinement when adding/deleting labels (or Haplotype SNPs)
			     force PoutlierEnd to PoutlierEndSize during refinement when adjusting sizing (or Haplotype Indels) */
extern double PoutlierEndLabel; 
extern double PoutlierEndSize;
extern double IndelDelPoutlierEnd;/* < 0.0 means no change during Indel Deletion */

extern double LRbias;/* LR bias added to all LR values during refinement (but not before adding/deleting sites in the extension region, if extend >= 2) */
extern int LRbiasSwitch;/* force LRbias to LRbiasLabel during refinement when adding/deleting labels (or Haplotype SNPs)
			force LRbias to LRbiasSize during refinement when adjusting sizing (or Haplotype Indels) */
extern double LRbiasLabel;/* 0.0 means use regular value of LRbias */
extern double LRbiasSize;/* 0.0 means use regular value of LRbias */
extern double IndelDelLRbias;/* 0.0 means no change during Indel Deletion */

extern int MinSFSwitch;/* force sf to at least MinSFLabel during refinement when adding/deleting labels (or Haplotype SNPs)
			  force sf to at least MinSFSize during refinment when adjusting sizing (or Haplotype Indels) */
extern double MinSFLabel;
extern double MinSFSize;
extern int MinSFLabelIter;
extern int MinSFSizeIter;

extern double MinSFIndelDel;/* minimum sf value used during Indel Deletion */

extern int MinSRSwitch;/* force sr to at least MinSRLabel during refinement when adding/deleting labels (or Haplotype SNPs)
			  force sr to at least MinSRSize during refinment when adjusting sizing (or Haplotype Indels) */
extern double MinSRLabel;
extern double MinSRSize;
extern int MinSRLabelIter;
extern int MinSRSizeIter;

extern double MinSRIndelDel;/* minimum sr value used during Indel Deletion */

extern int RangeSwitch;/* Larger value of RefineRange to use when updating alignments during refinement in UpdateMap() */

extern int ViterbiSwitch;/* see -ViterbiSwitch option */
extern double ViterbiLabel;
extern double ViterbiSize;
extern int ViterbiLabelIter;
extern int ViterbiSizeIter;


extern double RefineEps1;/* log-likelihood precision required for endoutlier terms (UL,UR) during refinement (mprobeval) */
extern double RefineEps2;/* log-likelihood precision required for total alignment of one molecule (LR) during refinement (mprobeval) */
extern double RefineEps3;/* log-likelihood precision required for sub alignment terms (AL,AR) during refinement (mprobeval) */

extern int extend;/* If reference map is incomplete and query map alignments can extend beyond end of reference map :
		     1 : Allow maps to extend during during alignment & refinement, but do not extend original reference map during refinement
		     2 : Also cause original reference map to refined in the extension region */

extern int extendonly;/* refine extension region only (only makes sense with extend = 2 ) */
extern double MaxExtend;/* If > 0 : Limit extension refinement to MaxExtend (in kb) at each end (to speed up refinement) */
extern int UseEndoutliers;/* If > 0 : Use the endoutlier region of molecules to generate hypothesis sites in consensus during refinement */
extern double EndLen;/* If extendonly > 0 : The part of the original contig ends (in kb) that will be refined along with the extension region */
extern double EndLen2;/* 2nd arg of -extonly */
extern double MinOverlapRatio;/* If extend > 0 && Refine > 0 : Ignore extending nanomaps with less than the specified ratio of their length overlapping the original reference contig */
                              /* In Assembler : minimum pairwise overlap as a fraction of the smaller map */

extern double extendSplit;/* If > 0 : Clone part of contig up to a target site where at least max(extendSplit, extendFrac * coverage) endoutliers in the direction away from the cloned part are present.
			     Don't count endoutliers shorter than extendSplitLen or whose alignments ends more than extendSplitFN sites before target site or 
			     more than extendSplitFP sites after target site (suspected breakpoint).
			     If multiple nearby labels satisfy the rule, use the one with the largest number of endoutliers.
			  */
extern double extendFrac;
extern int extendSplitFN;
extern double extendSplitFNlen;
extern int extendSplitFP;
extern double extendSplitLen;
extern int extendSplitSpacing;
extern double extendSplitSpacingX;
extern int extendSplitEndSpacing;
extern double extendSplitEndSpacingX;

extern double extendWT;
extern double extendLen;/* Alignments that extend this distance beyond end of contig (not counting unlabeled mol end) have their align->mapWT increased to extendWT (see -BestRefWT) */
extern int extendFN;
extern double extendFNlen;
extern double extendMinWT;
extern double extendMaxOutlierKB;
extern double extendMaxOutlierKBsplit;
extern int extendMaxOutlierLabels;
extern int extendMaxOutlierLabelsSplit;
extern double extendMaxWT;

extern double extendWTdup;/* If 0, only raise the weight of one extending alignment per molecule */
extern double extendWTend;/* If 0, -extendWT only applies to split contig extending alignments */

extern double splitWT;/* If splitWT > 0.0 : Increase weight of alignments to at least splitWT when computing coverage and endoutlier counts to determine splits.
			Also if extendWT > 0 : Apply extendWT (increase align->mapWT to extendWT) for endoutlier maps moved to splits as extensions */
extern double extsplitMinWT;/* skip alignments with WT below this value when computing -extsplit coverage */
extern double refineMinWT;/* skip alignments with WT below this value during refinement */

extern double refineWT;/* If refineMinWT > 0 AND refineWT > 0, alignments with WT above refineMinWT will have their WT increased to at least refineWT */
extern int refineWT_N;/* see 2nd value of -refineWT */
extern double refineWT_NL;/* see 3rd value of -refineWT */
extern double refineWT_MaxOutlierKB;/* see 4th value of -refineWT */
extern int refineWT_MaxOutlierLabels;/* see 5th value of -refineWT */

extern double splitFiltMinWT;/* see -splitFilt */
extern int splitFiltFN;
extern double splitFiltFNlen;
extern int splitFiltFP;

extern int keepsplitmaps;/* If > 0 : keep maps with endoutliers that are used in a contig split in the original contig as well, but at reduced weight based on -BestRefWT */
extern int extTrim;/* see -extTrim (EXTSPLIT_KEEPMAPS in refalign.cpp) */

extern int ExtsplitSpacingFix;/* see -extsplitSpacingFix & EXTSPLIT_SPACING_FIX in refalign.cpp */
extern int ExtsplitVerifyPeak;/* see -extsplitVerifyPeak & EXTSPLIT_VERIFYPEAK in refalign.cpp */

extern int extsplitNE;/* If 1 : include non-endoutliers that align only to the split contig region in the split contig
			 If 0 : Don't include non-endoutliers in split contig, and limit -extonly and -maxExtend value for each split contig (see -extsplitNE) */
extern double extsplitNEL;/* L value of -extsplitNE */
extern double extsplitNEE;/* E value of -extsplitNE */

extern double extsplitOutlier;/* If > 0 : Also include outliers of least this size difference if either boundary matches a split point (along with endoutliers) */
extern int extsplitOutlierLabels;
extern double extsplitOutlierLS;
extern int extsplitOutlierAS;
extern double extsplitOutlierKB[MAX_KB_BINS+2];/* extsplitOutlierKB[1..KBcnt] correspond to L1 .. Ln in the -extsplitOutlier option, extsplitOutlierKB[0] == extsplitOutlier */
extern int KBcnt;

extern int KCcnt;
extern double extsplitOutlierC[MAX_KB_BINS];/* extsplitOutlierC[0..KCcnt-1] are the args of -extsplitOutlierC */

extern double CloneLimitLen;/* If > 0.0 : Limit length of -extsplit cloned part of contig to no more than this many kb or CloneLimitN labels (whichever is greater) */
extern int CloneLimitN;
extern int CloneTrim;/* If > 0 : Allow -extsplit split contig's non-extension end to be trimmed by -EndTrim (see -CloneTrim) */

extern int HasLoc;/* If startloc & endloc should be read in from Version 0.1 BNX input (if present) */

extern int InDel;/* If > 0 : output .indel file along with .xmap file (if -ref was used) */
extern int InDelEnds; /* If > 0 (AND Indel > 0) : Also output end outliers in .indel file */
extern int InDelNoOverlap;/* If > 0 (AND Indel > 0) : If the same query region aligns to two different reference region, do not output any indels fromt he overlapped region of the query */
extern double InDelNoOverlapT;/* minimum confidence for matchgroups to suppress indels with IndelNoOverlap > 0 */
extern double InDelMinConf;/* If > 0 (AND Indel > 0) : Indels must have minimum Left or Right alignment confidence of InDelMinConf, otherwise they are ignored. Useful to filter out indels from CutFlipPaste matchgroups */

extern double IndelChiSq;/* During -indel calling ignore outliers where smaller map interval has ChiSq value below this value */
extern double IndelFrac;/* During -inel calling ignore outliers where smaller map interval has OutlierFrac value above this value */

extern double IndelMaxSize;/* Only ignore outliers with size (delta) below this many kb */
extern double IndelMaxCov;/* Only ignore outliers with MolCov value below ~MaxCov */
extern double IndelMinSf;/* Only ignore outliers with MolSd above ~MinSf + ~MinSr * Interval */
extern double IndelMinSr;
extern double IndelMinInterval;/* Only ignore outliers within intervals above this many kb */
extern double IndelSdRatio;
extern double IndelMinCov;
extern double IndelMinCovScale;
extern int IndelWin;
extern double IndelFragile;
extern double IndelSegDup;
extern double IndelChimNorm;
extern double IndelMinCovRatio;


extern double BreakChiSq;/* During CMAP input break maps on intervals with ChiSq value below this value */
extern double BreakFrac;/* During CMAP input break maps on intervals with OutlierFrac value above this value */

extern double BreakMaxSize;/* Only ignore outliers with size (delta) below this many kb */
extern double BreakMaxCov;/* Only ignore outliers with MolCov value below ~MaxCov */
extern double BreakMinSf;/* Only ignore outliers with MolSd above ~MinSf + ~MinSr * Interval */
extern double BreakMinSr;
extern double BreakMinInterval;/* Only ignore outliers within intervals above this many kb */
extern double BreakSdRatio;
extern double BreakMinCov;
extern double BreakMinCovScale;
extern double BreakMinCovRatio;
extern int NoBreak;
extern double BreakFragile;
extern double BreakSegDup;
extern double BreakChimNorm;

struct CRunDataSort {
  int index;/* index into RunDataList[] */
  char *F12;/* Fields 1 & 2 (tab seperated) (including the "# Run Data<tab>" before Field 1) */
  char *date;/* Field 3.1 */
  char *time;/* Field 3.2 */
  char *AM_PM;/* Field 3.3 */
  char *F456;/* Fields 4, 5 & 6 (tab seperated) */
  int NumberOfScans;/* Field 7 */
  char *F8;/* Field 8 */
  long long ChipId;/* Field 8.x (last subfield seperated by ,) */
  int FlowCell;/* Field 9 */
  char *F10;/* Field 10 (NULL if not present) */
  double MinMoleculeLength;/* Field 11 (0.0 if not present) */
  double MinLabelSNR[MAXCOLOR];/* Field 12 (or 12 and 13 if colors==2) : 0.0 if not present */

  /* optional fields FeatureMinSlope and FeatureMaxSlope */
  double FeatureMinSlope;
  double FeatureMaxSlope;

  int RunId;/* RunId based on RunIndex */
  int ScanIdStart;/* Sum of NumberOfScans value for lower RunId (same ChipId,FlowCell) : To compute ScanId of molecule, add its (ScanNumber-1) to this value */
  int UniqueScanIdStart;/* Sum of NumberOfScans for all previous entries : To compute UniqueScanId of molecule, add its (ScanNumber-1) to this value */
};

/* class to track parameter values during iterative re-estimation */
#define Cparameters_Version 1004 // NOTE : Please update this number each time the below struct is changed

typedef struct {
  int version; // value of Cparameters_Version of RefAligner binary used to generate this binary

  int colors;/* number of colors used */
  int mapcnt;/* number of map alignments */
  int ATmapcnt;/* number of map alignments above threshold */
  int sitecnt;/* number of ref site intervals in maps that aligned above threshold (including internal outliers) */
  double sumX;/* sum of aligned qry lengths (including internal outlier) */
  double totlen;/* sum of qry map lengths */
  double MA_mean; // MA_SF, MA_SD are stored in SF[colors], SD[colors] respectively
  double PixelLen;/* bpp */
  double PixelLenSD;/* standard deviation of bpp */
  double mres,mresSD;
  double logLR;/* negative log Likelihood Ratio per map that aligned */
  double ATlogLR;/* negative log Likelihood Ratio per map that aligned */
  
  /* bias parameters and observed outlier Rates*/
  double outlierRate, EndoutlierRate;

  double maxresbias;
  int MaxResBins;/* value of RESBINS of RefAligner binary used to generate this binary */
  int ResBins[MAXCOLOR];/* actual size of polynomial must be <= MaxResBins */
  double resbias[MAXCOLOR][RESBINS+1];
  double resbiasX[MAXCOLOR][RESBINS+1];
  
  double res[MAXCOLOR];
  double resSD[MAXCOLOR];

  double FP[MAXCOLOR];
  double FPfrac[MAXCOLOR];
  double FN[MAXCOLOR];
  double SF[MAXCOLOR+1];
  double SD[MAXCOLOR+1];
  double SR[MAXCOLOR+1];
  double SE[MAXCOLOR+1];

  double LabelDensity[MAXCOLOR];/* reference label density of aligned regions */
  double LabelDensityX[MAXCOLOR];/* query label density of aligned regions */
  double minSNR[MAXCOLOR];
} Cparameters; 

typedef struct {
  int version; // value of Cparameters_Version of RefAligner binary used to generate this binary

  int colors;/* number of colors used */
  int mapcnt;/* number of map (splits) that aligned */
  int ATmapcnt;/* number of maps that aligned above threshold */
  int sitecnt;/* number of sites in maps that aligned above threshold */
  double MA_mean; // MA_SF, MA_SD are stored in SF[colors], SD[colors] respectively
  double PixelLen;/* bpp */
  double PixelLenSD;/* standard deviation of bpp */
  double mres,mresSD;
  double logLR;/* negative log Likelihood Ratio per map that aligned */
  double ATlogLR;/* negative log Likelihood Ratio per map that aligned */
  
  /* bias parameters and observed outlier Rates*/
  double outlierRate, EndoutlierRate;

  double maxresbias;
  int MaxResBins;/* value of RESBINS of RefAligner binary used to generate this binary */
  int ResBins[MAXCOLOR];/* actual size of polynomial must be <= MaxResBins */
  double resbias[MAXCOLOR][64+7+1];
  double resbiasX[MAXCOLOR][64+7+1];
  
  double res[MAXCOLOR];
  double resSD[MAXCOLOR];

  double FP[MAXCOLOR];
  double FPfrac[MAXCOLOR];
  double FN[MAXCOLOR];
  double SF[MAXCOLOR+1];
  double SD[MAXCOLOR+1];
  double SR[MAXCOLOR+1];
  double SE[MAXCOLOR+1];

  double LabelDensity[MAXCOLOR];
  double minSNR[MAXCOLOR];
} Cparameters1003; 

/* Parameters only used by RefAligner */

extern int NoSplit; /* limit splitting of nanomaps due to local/chimeric alignment : 
		       1 : Do not re-split previously split map fragment 
		       2 : do not do any map splitting or chimeric map detection */
extern int NoSplitSet;/* 1 IFF -nosplit was specified on command line */

extern int PairSplit;/* != 0 :  compute multiple local alignments during pairwise alignments (requires (Ch>0 && ChFdr>0) OR PoutlierEnd > 0.0)
			 < 0 : disallow splits on 1st set of maps (reference maps) */
extern double Psplit;/* -ve scores larger than -log(1/Psplit) are used to break apart local alignments internally (only used if PairSplit > 0) */
extern int PairsplitExtend;/* If PairSplit : When breaking molecules at internal outliers extend both pieces to include the outlier region (including any other outliers nearby) */
extern int PairsplitMerge;/* If PairSplit : At the end try to merge overlapped or adjacent alignments */
extern int PairsplitMinOutlier;/* If PairSplit  : Do NOT break alignments at interval outliers with less than this many site intervals on both maps (intended to work with -indel) */
extern int PairsplitPV;/* If PairSplit : Prioritize alignments by LogPV instead of Likelihood score */
extern int PairsplitRef;/* If PairSplit : treat first -i input map set as in-silico reference with assumed zero errors */
extern int PairsplitSeparateQueries;/* If PairSplit : treat multiple queries as if each query were processed in a seperate command with the same set of reference maps */

extern int PairwiseOutlierFix;/* If > 0 : fixes several issues with the handling of outliers and endoutliers during pairwise alignment */

extern int RefIndex;/* If PairSplit : Reference Map last index into vfixx_filename[0..num_files-1] */
extern int QueryIndex;/* If PairSplit : Query Maps first index into vfixx_filename[0..num_files-1] */

extern int CmapMerge;/* -merge utility selected */
extern int CmapMergeBnx;/* output -merge maps as .bnx */
extern int CmapSort; /* Order in which maps are output when using -merge :
			0 : Output maps in any order (typically the same order as they were input) 
			1 : Sort maps in increasing order map size
			2 : Sort maps in increasing order of number of sites
			3 : Sort maps in increasing order of map id
			4 : Randomize order of maps (helpful to avoid memory bottlenecks) */
extern double CmapSortMinLen;

extern int CmapSwap;/* If >= 1 : swap the two colors in the input maps */

extern int RandSeed;
extern int CmapStart;/* If >= 1, the subset of maps to output (based on map order after sorting) */
extern double CmapEnd;
extern double CmapFrac;/* see -subset : output at least this fraction of total maps  */
extern int perCohort;/* see -subset : output at least this many maps per cohort bin (same UniqueScanId and RunIndex/ScanNumber combination) */

extern double CmapGB;/* If > 0.0, the subset of maps to output (based on map order after sorting) is 1 .. N, where N is the value such that the molecules 1 .. N that exceed CmapMinLen kb, add up to CmapGB * 1000000 kb */
extern double CmapMinLen;
extern int CmapGBrnd;/* If 1 : round up N to include all molecules from the same cohort (RunIndex) */

extern int FinalSort;/* If >= 1, resort both -i and -ref maps (after -i maps have been processed with -sort* or -randomize and -subset) */

extern int CmapBin, CmapNumBins;/* If >= 1, the subset of map to output (based on partition the map into CmapNumBins roughtly equal sized sections and outputing section CmapBin */
extern int CmapSplit;/* If >=1, split input maps into single map files */

extern double RangeLeft,RangeRight; /* If RangeRight > RangeLeft : The range (in kb) of a single input map to output as .cmap file */
extern double InversionLeft, InversionRight; /* If InversionRight > InversionLeft : The range (in kb) of the single input map to be inverted in the out */

extern double PairMerge;/* output pairwise merged maps in CMAP format when alignment satisfies -S -T -A thresholds and has minimum overlap in kb >= PairMerge.
			Each input map must be in a seperate input file and have distinct ID values
			Unmerged maps are written to PREFIX_contig<N>.cmap, where <N> is the ID of the map.
			Merged map pairs (starting with highest Pvalue) are written to PREFIX_contig<M>.cmap, where <M> is the smaller
			of the two map IDs (the higher map ID will be missing from the output files). The merged map transitions between the
			original maps at the midpoint of the alignment */
extern double PairMergeOverlapped;/* If >= 0.0 : Reduced minimum overlap in kb when one map completely overlaps another map */
extern double PairMergeOverlappedMargin;/* see 2nd arg of -pairmergeOverlapped */
extern double PairMergeOverlappedSmargin;/* see 3rd arg of -pairmergeOverlapped */
extern double PairMergeOverlappedNmargin;/* see 4th arg of -pairmergeOverlapped */
extern double PairMergeOverlappedEmargin;/* see 5th arg of -pairmergeOverlapped */

extern double PairMergeOutlierChiSq;/* During -pairmerge ignore outliers where either map interval has ChiSq value below this value */
extern double PairMergeOutlierFrac;/* During -pairmerge ignore outliers where either map interval has OutlierFrac value above this value */

extern double PairMergeOutlierMaxSize;/* Only ignore outliers with size (delta) below this many kb */
extern double PairMergeOutlierMaxCov;/* Only ignore outliers with MolCov value below ~MaxCovs */
extern double PairMergeOutlierMinSf;/* Only ignore outliers with MolSd above ~MinSf + ~MinSr * Interval */
extern double PairMergeOutlierMinSr;
extern double PairMergeOutlierMinInterval;/* Only ignore outliers within intervals above this many kb */
extern double PairMergeOutlierSdRatio;
extern double PairMergeOutlierMinCov;
extern double PairMergeOutlierMinCovScale;

extern double PairMergeChimNorm;
extern double PairMergeOutlierMinCovRatio;


extern double PairMergeMaxEnd;/* Apply PairMerge only if each unaligned end (due to -endoutlier) does not exceed "overlap region" * PairMergeMaxEnd */
extern double PairMergeMaxOutlier; /* If >= 0 : Apply PairMerge only if the largest internal outlier does not exceed PairMergeMaxOutlier (in kb) on either map */
extern double PairMergeMaxEndKB;/* Apply PairMerge only if each unaligned end (due to -endoutlier) does not exceed PairMergeMaxEndKB (in kb) OR does not exceed PairMergeMaxEndN labels (at either end)  */
extern int PairMergeMaxEndN;
extern int PairMergeMaxOutlierM;

extern int PairMergeBigMapsFirst;

extern double PairMergeMaxEndExpandKB;/* Allow outliers near end to be removed by merging with endoutlier, as long as endoutlier does not expand to more than this value (defaults to PairMergeMaxEndKB) */

extern int PairMergeMaxOutlierBalanced;/* If outlier contains this many or more misaligned labels, use inversion size instead of Indel size (if larger) as outlier size */

extern int PairMergeExcludeOutliers;/* If 1 : Exclude internal outliers in the computation of the alignment overlap size (in kb) : 
				    Overlap size is used by various suboptions of PairMerge */

extern int PairMergePreferNGS; /* When merging preserve as much of the NGS map (StdDev=0) as possible */
extern int PairMergeXmap;/* output .xmap files during pairmerge */
extern int PairMergeRepeat;/* If >= 1 : Repeat PairMerge operation until no more maps merge (avoids file IO between each operation) */

extern double TrimmedNoMerge;/* If >= 1 : Don't merge any pair of maps if that eliminates a map end marked as trimmed and non-extendable (if Mask[0][1] or Mask[0][N+1]) */
extern double TrimmedNoMergeMaxEnd;
extern size_t TrimmedNoMergeBoth;/* require both overlapped ends to be marked, before preventing merges */
extern size_t TrimmedNoMergeSegDupMask; /* If !=0 : Disallow merges if the overlap region of either map has Mask bits that intersect with this value, unless the overlap extends beyond the Masked region
					   by at least TrimmedNoMergeSegDupFlank labels in the non-extending direction */
extern int TrimmedNoMergeSegDupFlank;

extern int TrimmedNoMergeStrict;// see -TrimmedNoMergeStrict

extern double Truncate;/* If > 0.0 : For any map pair with one endoutlier and at least Truncate kb and TruncateN labels aligned with no outliers larger than TruncateMaxOutlier,
			 truncate the shorter map from the end without the endoutlier to no more than Trucate kb or TruncateN labels (whichever is greater) and mark the trimmed 
			 end non-extendable using Mask[0][1] or Mask[0][N+1] */
extern int TruncateN;
extern double TruncateMaxOutlier;
extern int TruncateFix;/* see 4th arg to -Truncate */
extern int TruncateMaxOutlierM;/* M value (max misaligned labels for non-outlier interval) of -Truncate */

extern double MaxPalindrome;/* If > 0.0 : -MaxPalindrome <L> value */
extern int MaxPalindromeN; /* -MaxPalindrom <N> value */
extern double MaxPalindromeSL;/* -MaxPalindrom <SL> value */
extern size_t MaxPalindromeMask;/* -MaxPalindrome D value */

extern double MaxInvPalindrome;/* If > 0.0 : -MaxInvPalindrome <L> value */
extern int MaxInvPalindromeN;/* -MaxInvPalindrome <N> value */
extern double MaxInvPalindromeSL;/* -MaxInvPalindrome <SL> value */
extern size_t MaxInvPalindromeMask;/* -MaxInvPalindrome D value */

extern double SplitSegDup;/* If > 0.0 : -SplitSegDup <L> value */
extern int SplitSegDupN; /* -SplitSegDup <N> value */
extern double SplitSegDupE;/* -SplitSegDup <E> value */
extern int SplitSegDupEN;/* SplitSegDup <EN> value */
extern double SplitSegDupC;/* -SplitSegDup <C> value */ 
extern double SplitSegDupLM;/* -SplitSegDup <LM> value */ 
extern int SplitSegDupNM;/* -SplitSegDup <NM> value */
extern double SplitSegDupMaxOutlier;/* -SplitSegDup <Mout> value */

extern size_t SegDupMask;/* If > 0, mark overlap region of segdups (even if smaller than SplitSegDup or SplitSegDupN) with this Mask value for all labels of both maps */

extern int SplitSegDup_OneSided;/* K value of -SplitSegDup_OneSided <K> <L> */
extern double SplitSegDup_OneSided_MaxLen; /* L value of -SplitSegDup_OneSided <K> <L> */

extern double BreakEndOutlier;/* If > 0.0 : -BreakEndOutlier <L> value */
extern int BreakEndOutlierN;/* -BreakEndOutlier <N> value */
extern double BreakEndOutlierE;/* -BreakEndOutlier <E> value */
extern int BreakEndOutlierEN;/* -BreakEndOutlier <EN> value */

extern double BreakOutlier;/* If > 0.0 : -BreakOutlier <L> value */
extern int BreakOutlierN;/* -BreakOutlier <N> value */
extern double BreakOutlierI;/* -BreakOutlier <I> value */
extern int BreakOutlierIN;/* -BreakOutlier <IN> value */

extern int PairMergeHmap;/* If >= 1 : After no more map pairs can merge, re-examing all pairwise alignments and try to combine as many as possible into Haplotype Maps and output them as HMAPs
			    First merge Bubbles (fully overlapped smaller map that could not be merged with larger map due to an internal outlier larger than PairMergeMaxOutlier
			    Next try to combine matched pairs pairs of endoutliers indicated a large incomplete insertion, and treat as a single insertion bubble (will be broken by Haplotyper if not completed).
			 Requires either -pairmergeRepeat OR that pairmerge was previously run until no more map pairs can merge */
extern double TruncateHmap;/* see -pairmergeHmap option : defaults to Truncate */
extern int TruncateNHmap;/* see -pairmergeHmap option : defaults to TruncateN */
extern double PairMergeEndHmap;/* If > PairMerge : Minimum overlap required for merging End Overlaps (with Bubbles) : useful to reduce higher risk of chimeric merges */
extern double MinMapLenHmap;/* Minimum map lengths in kb for end overlaps with bubbles : useful to limit higher risk of chimeric merges to high payoff merges likely to increase N50 significantly */


extern int PairMergeIDrenumber;/* If >= 1 : renumber output contig IDs sequentially from 1 : may have to output contigs that otherwise could be symbolic links */

extern int PairMergeIDrenumber;/* If >= 1 : renumber output contig IDs sequentially from 1 : may have to output contigs that otherwise could be symbolic links */

extern int XmapInterleaved;/* If >= 1 : Use interleaved label ids in 2-color XMAP alignments */

extern char *ForceMerge;/* Filename which contains one pairwise merge specification per line (6 integers) */

extern double origPixelLen;/* original PixelLen value before -bpp input */
extern double refPixelLen;/* PixelLen value when reference map is read and de-res'd */

extern int SecondBest;/* to output 2nd best alignment score and right end location in .map and .xmap */
extern double SecondSpacing;/* spacing (distance of rightmost aligned sites) of 2nd best alignment from best alignment as multiple of Query Map length */
extern int SinglelineXmap;/* output 2nd best match on same line as best match (non-standard XMAP format) */

extern int MultiMatches;/* to output all alignments that score above ScoreThreshold2 & LogPvThreshold2 and have at least MultiMatches unique labels vs any higher scoring match */
extern double MultiMatchesDelta;/* If > 0.0 : Also output overlapped alignments that do not have enough unique labels, provided the offset between maps differs by at least <D> kb */

extern int MultiMatchesUniqueQuery;/* see -MultiMatchesUniqueQuery 1st arg */
extern int MultiMatchesUniqueQueryInv;/* see -MultiMatchesUniqueQuery 2nd arg */
extern int MultiMatchesUniqueQueryTrim;/* see -MultiMatchesUniqueQuery 3rd arg */

extern int RefSplit; /* If > 0 : Split each alignment (but not maps) at internal outliers with Pvalue < Psplit/bc and at least RefSplitMinLabels site intervals in both maps
			Note that with RefSplit outliers can consist of any alignment interval with -ve score < log(Psplit) provided no subinterval with score > log(1.0/Rsplit) is included
			and the first and last subinterval has a -ve score. 
			NOTE : bc == (outlierBC ? m*n : 1) and m and n are the number of site intervals involved in query and reference respectively.
		     */
extern double Rsplit;/* see RefSplit : If 0, use Psplit instead */
extern int RefSplitMinLabels;/* see RefSplit */
extern int MultiMatchesFilter;/* If > 0 : Filter matchgroups by LogPvThreshold2 and AlignedSiteThreshold2 before simimilarity check, to avoid discarding matchgroup with low score but high confidence */
extern int RefSplitStitch;/* If > 0 : After RefSplit, check if any remaining matchgroups can be stitched together (end to end). Requires that filtering by all thresholds is only applied either before any matchgroups are trimmed
			  (see MultiMatchesFilter) or again at the very end (after stitching any trimmed matchgroups) */
extern int MultiMatchesRev;/* If > 0 : Run MultiMatches twice with reference in either orientation and then eliminate duplicate matchgroups */
extern int MultiMatchesTotScore;/* If >= 1 : Resolve overlapped matchgroups to maximize the total score of unique portions by avoiding internal outliers at the expense of losing the highest scoring single alignment
				   If >= 2 : Do this before applying -MultiMatchesDelta filter, in case alignments have multiple outliers whose sizes cancel each other so the ends have only a small offset.
				             This takes extra time, so only do this if at least one outlier is MultiMatchesTotDelta kb or larger (should match size of smallest tandem duplication to be detected) */
extern double MultiMatchesTotDelta;/* see MultiMatchesTotScore  */

extern double MultiMatchesDupesMinscore;/* In MultiMatchesDupes() : Minimum absolute score of shared alignment, before breaking up lower scoring alignment */
extern double MultiMatchesDupesMinscoreRel;/* In MultiMatchesDupes() : Minimum scores of shared alignment as share of total alignment score, before breaking up lower scoring alignment */
extern double MultiMatchesDupesMMPen;/* In MultiMatchesDupes() : If MultiMatchesTotScore > 0 : Minimum total score improvment required before breaking up higher scoring alignment */


extern int bin, numbins, skipcnt;/* partial pairwise alignment parameters : output specified bin out of 1..numbins */

extern int Kbiasinit,KLbiasinit,KFbiasinit,Gbiasinit,Betainit;

extern int RefRepeats;/* number of times to repeat parameter estimation in refalign */
extern int RefRepeats2;/* number of times to repeat calls to hash_generate() followed by refalign_allpairs() */
extern int SaveRepeat;/* If >= 0 : .errbin file will be based on iteration SaveRepeat,RefRepeat2 (instead of last iteration RefRepeat,RefRepeat2) */
extern double RefRepeatErr;/* stop parameter estimation if no parameter changed by more than this fraction */

extern int Refine;/* apply refinement step to consensus at end of refalign() and output refined consensus */
extern int RefineRange;/* Adjust RANGE parameters used during refinement */

extern int RefineXRange;/* If != 0 : extrapolate alignments (when needed) using both Xtheta * (RefineRange + max(deltaXRef,deltaExtXRef)) instead of just Xtheta * RefineRange */
extern int RefineAlignExt;/* If != 0 : make sure alignments are extended to the ends of each molecule (makes score less dependent on recent reposition call) */
extern int RefineFillMap;/* see FILL_MAP in refine.h */
extern int RefineFillHmap;/* see FILL_HMAP in refine.h */
extern int RefineSwitchFade;/* If != 0 : turn of RefineSwitch* options after initial Haplotyping (when stringency of Indels is increased) so scoring functions becomes almost the same for Indels and SNPs */

extern int LocalTheta;/* If != 0 : Use local estimate of consensus label interval size instead of Xlambda or Xtheta when adding labels with SCANRANGE,SPLITRANGE,EXPANDSPLIT in refine.h */
extern double MaxLocalTheta;/* If > 0.0, Use this value as max local estimate of consensus label interval size */
extern int LocalRange;/* If != 0 : Use label interval size of current molecule instead of global Xtheta in setlimit() : see RANGE_Y=0 in refine.h */

extern double ScanRange;/* see -ScanRange and SCANRANGE in refine.h */
extern double HScanRange;/* see -HScanRange and HSCANRANGE in refine.h */

extern double MaxContigSiteDensity;/* If > 0 : skip refinment (use -refine 1) for contigs with more than this many sites / 100kb */ // NEW 

extern long long CMapID;/* CMapID used in output of refined consensus map */
extern long long startCMapID;/* Value of -id command line option (often a group id of multiple ref ids that will replace CMapID by turn) */
extern int SplitRef;/* 1 IFF each -ref input contained at most 1 map */
extern int SplitMap;/* 1 IFF each -i input file contained at most 1 map */
extern double ContigSplitRatio;/* If > 0, split input reference map at internal site interval if sitecntFN drops below ContigSplitRatio * "median fragcnt" (coverage) AND below ContigSplitRatio2 * "median sitecntFN"  */
extern double ContigSplitRatio2;
extern double ContigSplitRatio3;/* The split location must have sitecntFN at least this factor below the peak values on either side */
extern double MinSplitLen;/* minimum length of contig split interval to output */
extern int SplitRev;/* If > 0 : Revised -contigsplit code */
extern int SplitSite;/* If > 0 : Apply -contigsplit only to coverage at site locations (and not between sites, where coverage is more noisy). Particularly useful with -TrimNorm */
extern long long SplitCnt;/* Revised -contigsplit contig id asignment : Determine new contig ids by repeatedly adding the value in ContigCntFile to the original contig ID */
extern double SplitFragileD;/* If > 0.0 : <D> value of -splitFragile <D> <R> <E> */
extern double SplitFragileR;/*            <R> value of -splitFragile <D> <R> <E> */
extern int SplitFragileE;   /*            <E> value of -splitFragile <D> <R> <E> */

extern double HapSitePvalue;/* If HapSitePvalue > 0.0 : during refinement haplotype sites will be detected when haplotype Pvalue < HapSitePvalue */
extern double HapIndelPvalue;/* If HapIndelPvalue > 0.0 : during refinement haplotype Indels will be detected when haplotype Pvalue < HapIndelPvalue */
extern double PhasePvalue;/* If HapSitePvalue > 0.0 && PhasePvalue > 0.0: during refinement haplotype SNPs (& indels) will be Phased (with the next one) whenever this can be done with phasing Pvalue < PhasePvaluea */
extern double HapSiteRes;/* If HapSiteRes > 0.0 : Suppress SNPs within this distance (in kb) of another site. The value should be related to the imaging resolution res*PixelLen*/
extern double HapSiteResDelay;/* If HapSiteResDelay > 0.0 : Suppress SNPs within this distance (in kb) of another site on same Allele during the first HapSiteResN iterations.
				  If HapSiteResL = 2 : Also suppresss SNPs within this distance of another site on the opposite Allele during the firs tHapSiteResN iterations. */
extern int HapSiteResN;
extern int HapSiteResL;
extern int HapSiteResA;

extern double HapSiteUnphasedSNR;/* Filter out Unphased HapSites with with geometric mean SNR below this value */
extern double HapSiteUnphasedPvalue;/* If > 0.0 : Filter out Unphased HapSites or HapIndel with Pvalue < HapUnphased Pvalue */
extern double HapMinCoverage; /* Filter out HapSite or HapIndel if local coverage is below HapMinCoverage (molecules) */
extern double HapAlleleCoverage;/* Filter out HapSite or HapIndel if minor allele coverage as a fraction of local total coverage is below HapMinorCoverage (fraction) */
extern double HapAlleleCovInt;/* Increase HapAlleleCoverage for Indel by HapAlleleCovInt * "largest label interval in region of indel / 100kb" */
extern double HapAlleleProb;/* Molecules are assigned to an allele when likelihood (conditional on molecule belonging to either of 2 alleles) exceeds HapAlleleProb */
extern double HapMinCovSize;/* If > 0.0 : Suppress HapAlleleCoverage filter if Hap Indel size is larger than HapMinCovSize (large Hap Indels may have fewer molecules 
			       of the larger Allele, if the smaller Allele was assembled first and Haplotyping was only run during Final Refinement) */
extern double HapMinCovPvalue;/* If > 0.0 : Suppress HapAlleleCoverage & HapMinCoverage filter if Hap Indel has Pvalue > HapMinCovPvalue */
extern double HapMinCovPvalueL;/* If > 0.0 : Suppress HapMinCoverage filter if Hap Indel has Pvalue < HapMinCovPvalueL */
extern double HapMinCovAbs;/* If > 0.0 : Always filter out Hap Indels if Minor Allele absolute coverage drops below this value (except if HapMinCovAbsM,HapMinCovN disables it) */
extern double HapMinCovAbsM;
extern int HapMinCovN;/* Never filter out Hap Indels overlapped by more than this many SNPs AND more than HapMinCovAbsM SNPs/100kb */

extern double HapMinCovPL1;/* If < 1.0 : Don't filter out HapIndel if ChiSquare Pvalue for sizes of minor allele is above HapMinCovPL1 */
extern double HapMinCovPL2;/* If > 0.0 : Always filter out HapIndel if ChiSquare Pvalue for sizes of minor allele are below HapMinCovPL2 */
extern int HapMinCovL;/* If 0 : Both versions of Pvalue must agree (both must exceed PL1 or be below PL2) before -HapMinCov action is overridden 
			 If 1 : Only Pvalue based on autoNoise sizing error estimate is used.
			 If 2 : Only Pvalue based on major allele sizing error is used */
extern double HapMinCovCL;/* If > 0 : Minimum coverage of minor Allele required for -HapMinCovPL override to occur */

extern int checkMinorIsBigger;/* If 0 : Allow HapMinCov to apply to Hap Indels even if minor Allele is bigger in size
				 If 1 : HapMinCov does NOT apply to Hap Indles IF minor Allele is bigger in size (assumes DNA knots only shrink DNA)
				          NOTE : HapMinCov filter always applies of total coverage is below HapMinCoverage */

extern double HapEndTrim;/* Trim Allele A ends up to the first label with HapEndTrim < (covA + N) * (covA+covB+covC + 2N)/(covA+covB + 2N), where N = HapEndTrimMin */
extern double HapEndTrimMin;/* Pseodo-Count correction for covA and covB */
extern double HapEndTrimMaxFrac;/* Maximum fraction of Labels that may be trimmed (otherwise undo all per-Allele trimming and rely on -EndTrim) */
extern double HapEndTrimRatio;
extern double HapEndTrimLen;

extern int FastSnpsStartup;/* If > 0 : If SNPs/Indels are already present start Haplotyping with this iteration instead of 0 */
extern int FastSnps;/* Don't apply fast SNPs detection iteration for the first FastSnps iterations : should be at least 1 to avoid phasing errors */
extern int FastIndel;/* Don't apply fast HapIndel detection iteration for the first FastIndel Indel iterations : should be at least 1 to avoid phasing errors in absence of SNPs */

extern int HapIndelHiRes;/* Number of interations to use high res sizing ladder before reverting gradually to 1.4142 x sizing ladder. Larger values improve convergence at the expense of runtime */
extern int MAX_DELTA_OVERSAMPLE;/* how many times the number of samples can be reduced 2x after initial Indel scan has been completed (see LP_INDEL_SKIP)
				   Original sizing ladder step ratio = pow(2.0, 1/(2<<MAX_DELTA_OVERSAMPLE)) */
extern double MaxHapDelta;// largest Delta (half of largest haplotype Indel size) in kb that is checked during Haplotypeing

extern int SplitHmap;/* If 1 : Input maps that are already Haplotyped are split into two maps corresponding to the 2 Alleles (See RefineHmap on two ways to proceeed thereafter)
		       If 0 : Input maps that are already Haplotyped are treated as a regular Cmap that is the average of the 2 Allele maps.  */
extern int RefineHmap;/* Requires SplitHap == 1:
			If 0 : Input maps that are Haplotyped are separated into 2 Allele maps and each one Refined + Haplotyped independently of the other.
			If 1 : Input maps that are Haplotype are preserved and Haplotyped further starting from the previous Hap map. In this case non-Haplotype refinement
		               will NOT be performed before Haplotype refinement
		     */

extern int HapSiteMin, HapIndelMin;/* Only output separate <N>1.cmaps and <N>2.cmaps and <N>0.hmap IF at least HapSiteMin different Haplotype Sites and HapIndelMin different Haplotype Indels were found 
				       (other wise output just one <N>0.cmap). Decision is made after filtering HapSites and HapIndels based on -HapUnphased or -HapMinCov */
extern int HapSiteWin; /* Only output seperate Allele maps if at least HapSiteWin Het Sites are present in some consensus interval of size HapSitewinSize (kb) or less */
extern double HapSiteWinSize;
extern double HapIndelMinSum; /* Sum of absolute Hap Indel sizes must add up to at least HapIndelMinSum, otherwise ignore Hap Indels when deciding if seperate Allele maps should be output */
extern double HapIndelMinDensity;
extern int HapIndelSkipEndSNP;
extern double MinHapCovThresh;/* see -HapThresh 8th arg */
extern double MinHapRelCovThresh;/* see -HapThresh 9th arg */

extern int HapIndelMerge;/* If > 0 : Treat consecutive Hap Indels seperated Hap Sites as a single indel with a single HapIndelPvalue penalty */

extern int HapIndelDelskip; /* If > 0 : Perform expensive Hap Indel deletion step then skip next HapIndelDelskip times if LP change was less than HapSkipEps + HapSkipEpsPermap * MD */
extern int HapUpdatemapSkip;/* If > 0 : Perform expensive Hap Updatemap then skip next HapUpdatemapSkip times if LP change was less than HapSkipEps + HapSkipEpsPermap * MD */
extern double HapSkipEps;
extern double HapSkipEpsPermap;

extern double IntervalEps1;/* If LP change for all interval size changes is less than this amount before -rres & -cres, terminate */
extern double IntervalEps2;/* If LP change for all interval size changes is less than this amount after -rres & -cres, terminate */
extern double HapIndelEps1;/* If LP change for all Hap Indel size changes is less than this amount before -rres & -cres, terminate */
extern double HapIndelEps2;/* If LP change for all Hap Indel size changes is less than this amount after -rres & -cres, terminate */
extern double HapLabelEps1;/* if LP change for all Hap Label changes is less than this amount before -rres & -cres, terminate */
extern double HapLabelEps2;/* if LP change for all Hap Label changes is less than this amount after -rres & -cres, terminate */

extern int HapFast; /*  If >= 1 : Speed up Haplotyping by recomputing likelihood only where consensus was changed (will still occassionally check full likelihood and backtrack if needed)
		     If >= 2 : Avoid global re-alignment of molecules to handle large changes on consensus  */
extern int HapLocalStop;/* Terminate Haplotyping locally in large contigs and only work on regions that are still making progress in at least one Allele, based on -HapMinLL */
extern int HapLocalScan;/* value of DELTA_SCAN in refine.cpp */
extern int HapLocalScan2;

extern int HapStageIter;
extern double HapStageFrac;

extern int HREVERSAL;
extern int FINAL_HREVERSAL;
extern int HAP_MAX_DERES;
extern int DERES_STOP;
extern int FILTER_STOP;

extern int SKIPMAP_ADD;

extern int HapMaxDeltaIter;/* Maximum number of iterations in a row with only sizing changes (no HapIndel or SNP changes) before termination is forced */
extern int HapMaxSwitchIter;/* Maximum Hap Indel iterations before *Switch options must be turned off */

extern int HapFilterLast;/* If >= 1 : Apply Haplotype filters last (after -rres,-cres based label merging) */
extern int HapMergeSNPs;/* If >= 1 : Merge nearby SNPs on opposite Alleles into a regular site, if score improves */

extern char *ContigCntFile;/* If ContigSplitRatio > 0, this file will contain a single integer <N> that must be atomically incremented to obtain the 
			      name of any additional output contig as PREFIX_contig<N>_refined.cmap */

extern int Cfirst; /* see option -first */
extern int CfirstPairs;/* see option -firstPairs */

extern int EVERB; /* <= 1 : Enable display of each missing/false cut (in last iteration)
		     >= 2 : Disable display of each missing/falst cut */

extern double err_dist;/* count errors (fp and fn) that are isolated from true positives by at least this distance in kb */

extern char *hash_prefix;/* output .hash or .hashbin file name prefix */
extern char *hash_filename;/* input .hash file name : list of map id pairs with hash scores */
extern int hash_threshold;/* hash score threshold : only use pairs with hash score at least this large */
extern int hashsplit;/* If != 0 : Program just splits hashtable into multiple pieces */
extern double hashdelta;/* If != 0 : limit relative offset of aligned maps to value in hashtable +- this multiple of SD of length of smaller map */
extern double hashdeltaAdjust;/* If ! = : increase hashdelta for a particular map pair, if needed, in increments of hashdeltaAdjust */
extern double hashdeltaLim; /* If != 0 : Limit increases in hashdelta to no more than hashdeltaLim */
extern int HASHRANGE;/* If != 0 : based hashdelta range on number of sites rather than distance : faster when reference label density is higher than query label density */

extern int scale_output;/* output scaling factor for each query map in .map file */

extern int COVERAGE_TRIM; /* skip first and last few aligned sites of each map when computing site coverage : this will expose shallow regions spanned by chimeric molecules 
			     that tend to have few aligned sites on the "bad" side of the chimeric junction (see refine.cpp) */
extern double COVERAGE_TRIM_LEN;/* Also skip at least this many kb of each aligned map end */
extern int mCNT;/* number of alternate values for COVERAGE_TRIM_LEN */
extern double mCOVERAGE_TRIM_LEN[MAX_CHIM];/* If mCNT > 0 : multiple alternate values for COVERAGE_TRIM_LEN : mCOVERAGE_TRIM_LEN[0] == COVERAGE_TRIM_LEN */

extern double TrimNorm;/* If != 0 : Normalize Trimmed Coverage by total number of molecules that extend beyond both sides of the target site by at least COVERAGE_TRIM_LEN kb and COVERAGE_TRIM real sites (must be > 0) :
		              Coverage will be expressed on scale 0-100 (%). The TrimNorm value is used as a pseudocount prior to avoid breaking low coverage regions (eg fragile sites) */
extern double TrimNormLen;/* If > 0 : Use largest normalization count for TrimNorm within TrimNormLen kb */

extern int TrimNormChim;/* If >= 0 : Try to suppress counting chimeric molecules during TrimNorm by not counting molecules whose nearest aligned site is more than TrimNormChim real sites (in Hcuts[]) away from the target site
			This will avoid breaking chimeric sites, but if set too small may start to fail breaking chimeric contigs as well due to false negatives. A value of 1 or 2 should be OK */
extern int TrimNormFP;/* If >= 0 : Try to suppress counting chimeric molecules during TrimNorm by not counting molecules whose alignment ends more than TrimNormFP sites beyond target site, but does NOT extend CovTrimLen sites beyond target site */

extern int ChimQualityFix;/* some updates to ChimQuality computation */

extern double FragileQualityEnd;/* If != 0.0 : Compute Fragile Site Quality Score, by counting molecules whose ends are within this distance (in kb) from the suspected fragile site AND have no more than FragileQualityFN sites near the end or no more than FragileQualityFP sites beyond the suspected fragile site */
extern int FragileQualityFN;
extern int FragileQualityFP;

extern int OutlierQuality;/* If != 0 : Compute Outlier Quality score, by counting fraction of molecules with outliers crossing each interval */

extern int HapMapWT;/* If != 1 : Give full weight to maps in both Alleles, so that Quality scores of each Allele contig reflect all molecules, includes from the other Alelle */

extern double TrimNormMin;/* If TrimNorm >= 0 : Force Normalized Trimmed Coverage to zero if type N1 coverage is below specified value */
extern double TrimNormMinCov;/* If > 0.0, only apply TrimNormMin if raw coverage is above TrimNormMinCov */
extern double TrimNormMinEnd;/* Don't compute quality scores near ends of contigs until N1 coverage is above TrimNormMinEnd */
extern double TrimNormMinRel;/* If < 1.0 : only apply TrimNormMin if N1 is below TrimNormMinRel x "max value of N1 in either direction" */

extern int TrimNormEnd;/* If TrimNorm >= 0  && TrimNormEnd > 0 : Don't include molecules in Trimmed Coverage (numerator) if either end has an endoutlier with at least TrimNormEnd unaligned sites  */
extern int TrimNormMed;/* If TrimNorm >= 0  && TrimNormMed > 0 : Use TrimNormMed instead of median of Trimmed Coverage to compute threshold (multiplied by ContigSplitRatio2)  */
extern int TrimOutlier; /* If > 0 : For molecules alignments with outliers, treat each aligned region between outliers as a seperate molecule (Helps break consensus with spurious indels caused by coincident outliers) */
extern double TrimOutlierKB;
extern int TrimOutlierLabels;


extern int REPLACE_COV;/* If > 0 : output trimmed coverage (see COVERAGE_TRIM) in place of actual coverage in output_draft() */

extern double TrimFactor;/* used by Refine as initial trimfactor */

extern int MultiMode;/* Enable Multi Modal sizing optimization */
extern int Mprobeval;/* Enable new faster refinement code */
extern int MprobevalTest;/* Enable slower more thorough Mprobeval refinement by adjusting a number of hardwired parameters */

extern char *MappedPrefix;/* If !=0 : output those -i input maps that align with reference above threshold to <MappedPrefix>.cmap */
extern char *UnMappedPrefix;/* If !=0 : output those -i input maps that did NOT align with reference above threshold to <UnMappedPrefix>.cmap */
extern char *PartMappedPrefix;/* If !=0 : output those -i input maps that were only partly aligned (see -endoutlier) with a reference to <PartMappedPrefix>.cmap */
extern int MappedUnsplit;  /* If -1 use default heuristic, otherwise specifies if mapped output is NOT split (1) or split (0) into seperate files per reference contig (ignored if -grouped was specified) */
extern char *groupfile; /* If !=0 : mapped output is split into seperate files based on groups of reference contig ids specified in groupfile (one group per line, 1st line is group1 etc) */

extern double CenterMapped;/* If != 0 then MappedPrefix output is limited to only those -i input maps that align inside the center region of contigs not within this distance of contig ends */

extern int MaxCov;/* If > 0 : limit coverage at each consensus location to <C> by discarding maps with the lowest alignment Pvalue */
extern int MaxCovWt;/* If > 0 : Use weighted coverage instead of raw coverage with MaxCov : weights are based on -BestRefWt */

extern int PoutlierEndCnt;
extern double PoutlierEndRefine[16];/* PoutlierEndRefine[0..PoutlierEndCnt-1] is the sequence of outlier end probabilities used during successive phases of refinement. 
				    The values of PoutlierEnd will be added at the start unless the first value is smaller. The value 1.0 will be used for the final alignment and coverage computation.
				    Using a sequence from stringent to least stringent is useful to weed out mismatched extensions and (hopefully) converge to the consensus estimate of the extension */
extern double PoutlierEndFinal;/* Value of PoutlierEnd used during final alignment computation (if greater than PoutlierEndRefine[PoutlierEndCnt-1]) : does not change the consensus map but
				  can change the coverage profile and help reveal chimeric junctions by using PoutlierEndFinal = 1.0 */
extern int PoutlierCnt;
extern double PoutlierRefine[16];

extern double biasWT;/* weight theoretical score bias in refalign by this factor (1 = full bias) (0 = no bias, refalign score will match refine score) */
extern double PRbiasWT;/* weight of theoretical score bias correction for misresolved sites (1 = full bias, 0 = no bias) */
extern double biasWTend;/* If 1 : apply biasWT to endoutliers (this keeps the total Bias for each molecule independent of the alignment, but makes it harder to fairly score chimeric molecules or references) */
extern double biasWToutlier;/* If 1 : apply biasWT to outliers this keeps the total Bias for each molecule independent of the alignment, but makes it harder to fairly score molecules with outliers or SVs) */
			
extern double biasWTrefine;/* biasWT used during refinement except during initial extension refinement */
extern double PRbiasWTrefine;/* weight of theoretical score bias correction for misresolved sites (1 = full bias, 0 = no bias) */
extern double PRbiasWTrefine;/* weight of theoretical score bias correction for misresolved sites (1 = full bias, 0 = no bias) */
extern int biasWTrefineType;/* 0   : Use biasWTrefine only when computing final coverage and Chimeric Quality scores
			     > 0 : apply biasWTrefine during refinement (except during initial extension refinement) */

extern int RTHETA_FIX;/* Add sqrt(R*Theta) term to FL() and 1/sqrt(R*Theta) for FLE() etc : increases likelihood of outliers vs refalign (turn off to better match endoutliers in refalign) */

extern double DELTA_RES; /* resolution of MultiMode sizing optimization (2x faster if resolution is set to 0.01, but LP does not converge to global maxima as well as with 0.001) */
extern double DELTA_REL;/* relative resolution of MultiMode sizing optimization */
extern double Site_Pen;/* penalty per consensus map site in LP computation : requires to balance Maximum Likelihood bias in favor of too many sites (cf degrees of freedom penalty) */
extern double Site_Pen_Ext;/* penalty Site_Pen used during initial extension refinement stage : typically lower than Site_Pen */

extern double FN_Refine;/* If > 1.0 : Multiply FN value by this factor during refinement */
extern double FN_Ext;/* If > max(1.0,FN_Refine) : Multiply FN value by this factor during initial refinement stage */
extern double FNmin;
extern double FNminExt;

extern double MinSF_Refine;/* If > 0.0 : increase sf value to at least this value during refinement */
extern double MinSF_Ext; /* If > MinSF_Refine : increase sf value to at least this value during initial refinement stage */

extern double MinSR_Refine;/* If > 0.0 : increase sr value to at least this value during refinement */
extern double MinSR_Ext; /* If > MinSR_Refine : increase sr value to at least this value during initial refinement stage */

extern double skip_dist;/* speed up addition/deletion of sites during refinement by checking only every skip_dist (kb) */
extern double skip_dist_Ext;/* If > skip_dist : value of skip_dist used during initial refinement stages */
extern double skip_score;/* speed up addition/deletion of sites during refinement by NOT rechecking changes that drop LP by LP_SKIP as frequently */
extern double skip_score_Ext;/* If < skip_score : value of skip_score used during intial refinement states */

extern int dosvdetect; /* call output_smap */
extern bool dosv_simpleoverlap; /* smap: apply simple overlap criteria -- on by default */

extern int smap_zygosity; //put zygosity flag in smap
extern int smapSize;// add SVsize in smap
extern int smapFreq;// add SVfreq in smap
extern int smapTransOrient;// add Orientation column in smap (current NA unless it is a translocation)

extern double smap_Indel_Ref_Alignment_OverlapFrac;/* minimum fractional overlap in ref of Alignment (without SV) with an SV (in another contig/Allele) to flag the SV as heterozygous */
extern double smap_Indel_Ref_SV_Overlap;/* minimum overlap in ref (in kb) of 2 SVs for them to be considered the same. Can be -ve, which means a gap between ref regions is permitted (useful to handle repeat expansion or contraction, which may be called ambigously in different Alleles). */
extern int smap_Del_Ref_SV_OverlapIncrease;/* For deletions increase value of sma;_Indel_Ref_SV_Overlap by the size of the deletion */
extern double smap_Indel_Size_MatchRatio;/* max ratio of Indel Sizes (smaller size / larger size) for 2 indels to be considered the same */

extern double smap_confbylen; /* smap: minimum confidence (in xmap) per kb length of matchgroup : matchgroups that fail are not used to call any SVs */
extern double smap_ConfByInterval;/* smap: minimum confidence (in xmap) per aligned Label Interval of matchgroup : matchgroups that fail are not used to call any SVs */
extern double smap_side_alignkb; /* smap: minimum length (kb) of aligned regions (ie, xmap entries) on both sides of an sv */
extern double smap_sv_sizekb; /* smap: maximum abs difference in query and reference gap size for indels & inversions, otherwise they are changed into intra_trans */
extern double smap_inv_MaxGap;/* If > 0 : maximum reference or query gap size for inversions, otherwise they are changed into intra_trans */

extern double smap_sv_sizekb; /* smap: maximum length (kb) of sv (between xmap entries) -- no longer a maximum: now changes indel to translocation*/
extern double sv_indel_maxoverlap; // Obsolete (but see pairInversionBP() and smap_del_maxqryoverlap) 
extern double smap_del_maxqryoverlap;// max size (kb) of query overlap allowed (& removed) for deletions : larger overlaps are assumed to be FP repeat (or segdup) compression
extern double smap_del_bigqryoverlap;// If deletion query overlap is between this value (kb) and smap_del_maxqryoverlap, reduce confidence of deletion to minimum value (smap_min_conf)

extern double svInsMaxGapOverlap;// max fractional overlap of qry gap allowed for insertions 
extern double svDelMaxGapOverlap;// max fractional overlap of ref gap allowed for deletions

extern int    smap_min_end_nlabels; /* smap: minimum number of un-aligned labels required for an 'end' type SV */
extern double smap_min_end_lenkb; /* smap: minimum un-aligned length (kb) required for an 'end' type SV */
//extern double smap_inversion_ratio; /* smap: the ratio of gap on the reference vs query must be > this */
extern double smap_translocation_maxsep; /* smap: for transloactions with no overlap, maximum separation on query */
extern double smap_translocation_ovlfr; /* smap: for translocations with overlap, maximum overlap length/matchgroup length */
extern double smap_translocation_ovlmx; /* smap: for translocations with overlap, maximum overlap length */
extern double smap_InvMaxOverlapSize;/* smap: If > 0.0 : For inversions with query overlap, maximum overlap length in kb */
extern double smap_indel_tiny_maxsize; /* smap: tiny type indel is any indel whose size is <= this, in kb */
extern double smap_simpleoverlap_fraction; /* smap: in simpleoverlap mode, overlap of more than this means exclude */
extern double smap_translocation_mask_tol; /* smap: when using mask file, add this as fixed tolerance in kb */
extern double smap_min_conf; /* smap: do not report SVs with conf < this (unless not defined) */

extern int smap_DupSplit_maxgap;/* maximum query gap (in Labels) for duplication_split calls */
extern int smap_DupSplit_maxoverlap;/* maximum query overlap (in Labels) for duplication_split calls */
extern int smap_DupSplit_maxend;/* maximum labels in query beyond the ends of the two matchgroups for duplication_split calls */
extern double smap_DupSplit_MinQryLen;/* minimum length in query matchgroup to allow query to be extrapolated for a Duplication call, provided the unaligned portion of the end being extrapolated is less the smap_DupSplit_maxend sites */
extern int smap_DupSplit_MinQrySites;/* minimum labels in query matchgroup to allow query to be extrapolated for a Duplication call, provided the unaligned portion of the end being extrapolated is less the smap_DupSplit_maxend sites */
extern int smap_DupSplit_singleEnd;/* 0 : duplication_split calls with one end that cannot be extended (due to smap_DupSplit_maxend OR ~MinQryLen OR ~MinQrySites), are NOT allowed.
				      >= 1 : duplication_split calls with one end that cannot be extended, are allowed if MGs overlap on reference
				      >= 2 : duplication_split calls with one end that cannot be extended, are also allowed if nonextendable MGs is at least 50% in kb of total ref interval
				      >= 3 : duplication_split calls with one end that cannot be extended are always allowed.
				   */

extern int smap_Outlier_minQrySites;/* For Small inversion, duplication or translocation  the minumum number of site intervals in Outlier of larger Matchgroup */

extern int smap_Dup_minRefOverlapSites;/* For Duplication the reference overlap (in sites) must be at least this number of sites */
extern double smap_Dup_maxQryGapMult;/* For Duplication the reference overlap (in sites) must be at least this fraction of the smaller matchgroup reference size (in sites) */
extern double smap_Dup_maxSpacing;/* Maximum Gap (in kb) on Query between Duplication events (should correspond to max distance that can be phased) */

extern double smap_InvDup_minRefOverlapFrac;/* For Inverted Duplication the reference overlap (in sites) must be at least this fraction of the smaller matchgroup reference size (in sites) 
					      Also used to suppress inversions if overlap fraction in Sites is >= F for either ref or qry (but see -svInv_maxRefOverlapFrac & -svInv_maxQryOverlapFrac) */
extern double smap_Inv_maxRefOverlapFrac;/* If > 0 : suppress inversions if overlap fraction in Sites on Reference is >= F */
extern double smap_Inv_maxQryOverlapFrac;/* If > 0 : suppress inversions if overlap fraction in Sites on Query is >= F */

extern double smap_InvDup_maxQryGapMult;/* For Inverted Duplication the maximum qry gap (or overlap) in sites cannot exceed this fraction of the reference overlap (in sites) */
extern double smap_InvDup_maxSpacing;/* Maximum Gap (in kb) on Query between Inverse Duplication events (should correspond to max distance that can be phased reliably ???) */

extern double smap_Inv_MaxTrimFrac;/* In deciding between Inverted Duplication or Inversion, rule out Inversion if more than this fraction 
				     of Inversion matchgroup must be trimmed away or ignored for Inversion call */
extern double smap_Inv_MinTrimFrac;/* In deciding between Inverted Duplication or Inversion, rule out Inverted Duplication if less than this fraction 
				     of Inversion matchgroup must be trimmed away or ignored for Inversion call */

extern int smap_Inv_minRefSites;/* minimum number of non-overlapped sites in reference of inverted matchgroup needed to call inversion */
extern int smap_Inv_minQrySites;/* minimum number of non-overlapped sites in query of inverted matchgroup needed to call inversion */

extern double smap_Inv_ConfMinIntervalSize;/* minimum interval size to include in interval count to determine confidence of small paired inversions */
extern int smap_Inv_ConfMinNumIntervals;/* minimum number of intervals (above smap_Inv_ConfMinIntervalSize) to output small paired inversions */

extern int smap_IndelConfBySize;/* If 1, base Indel confidence on size difference only (ignoring misaligned labels) */
extern double smap_SmallDelConfirm;/* For Deletions in ref intervals below this size, confirm deletion by including neighboring intevals and update center confidence if smaller */
extern double smap_IndelRangeExpandRes;/* max range (in kb) of misresolved labels allowed at boundary of Indel ref range : expand range until this condition is satisfied and adjust Indel size estimate downward, if needed */
extern double smap_IndelRangeExpandSep;/* min seperation (in kb) of boundary labels from next ref label (in either direction) : expand range until this condition is satisfied and adjust Indel size estimate down */
extern double smap_IndelRangeMinSVchange;/* min relative change in SV size due to range expansion : If SV size does not change by this amount the range expansion is undone */

extern double smap_GapOverlapRatio;/* minimum fractional overlap of 3rd MG by gap in both ref and qry gap of a MG pair with same qry and ref, before MG pair is not a valid SV candidate */

extern double smap_PrimaryOrientation;/* A MG has a dominant orientation IFF total Qry sizes for all MGs with same qryid and refid in one orientation exceeds the total Qry sizes for the other orientation
					       by this ratio. */
extern double smap_TransMaxOverlap;/* If > 0 : Disallow translocation if 3rd MG overlaps one of the 2 trans MGs by at least a fraction R of the trans MG qry size and maps back to the other reference
				       NOTE : IF 3rd MG confidence exceeds minimum confidence for translocation MGs, include overlap with the gap in the total overlap fraction */
extern double smap_TransMaxGapOverlap;/* If > 0 : Disallow translocation if a 3rd MG overlaps the query gap of the 2 translocation MGs by a fraction G of the 3rd MG qry size.
					NOTE : 3rd MG confidence must exceed minimum confidence for translocation MGs */

extern double smap_Repeat_Tolerance;/* SV simple repeat tolerance as multiple of expected sizing error (sqrt(SF*SF+SR*SR*X*X), X = average of intervals) */
extern int smap_Repeat_Minele;/* Minimum consecutive repeat intervals */
extern double smap_Repeat_ConfPenalty;/* If > 0 : Compute pro-rata confidence of non-repeat intervals only and subtract this penalty per repeat region AND only mark with repeat if confidence drops below -T */

extern int smap_MATCHGROUP_TRIM;/* see MATCHGROUP_TRIM in output_smap.cpp */
extern int smap_MATCHGROUP_PARTIAL_TRIM;/* see MATCHGROUP_PARTIAL_TRIM in output_smap.cpp */

extern char *conf_filename;/* filename for input of confidence scaling table */
extern double svConfMinsize;/* If indel size is below this output svConfNAvalue instead of value from confidence scaling table */
extern double svConfNAvalue;

extern double Palindromic_Delta;/* mark xmap entry as palindromic if best alignment in opposite orientation had confidence within this Delta of normal threshold (LogPvThreshold/log(10)) */

extern int doscaffold; /* process scaffold file */

extern char *svcheck;/* check SVs in input .indel or .smap (currently only indels) by aligning -i input maps against two versions of each SV neighborhood and compute a Likelihood ratio for/against the SV being present.
			The consensus maps are read in from the query map file named in svcheck[] .indel/.smap file with the corrected regions read in from the reference map file named in svcheck[].
			Outputs corrected versions of each query map as <prefix>_contig<NNN>.cmap along with updated Pvalues for SVs checked in <prefix>.indel (otherwise same as input file svcheck[]) */
extern double svcheckKB;/* SV neighborhood size in kb (on each side of SV) */
extern double svcheckMerge;/* Try to merge SVs that are seperated by less than this many KB on the same alignment */
extern double svcheck_confidence;/* correct SV if its updated confidence drops below this value */

extern xmapEntry *SVtocheck;/* The SVs to be checked : SVtocheck[i = 0 .. numrefmaps/2 - 1] corresponds to refmap[2*i .. 2i+1] */

extern char *svcompare1, *svcompare2;/* compare SVs in these two .indel or .smap files : display total SVs and number that overlap by at least SVcompareKB bases or SVcompareFrac (fraction) */
extern double SVcompareKB, SVcompareFrac, SVDeltaMatch, SVConf,SVXmapOverlap;
extern int SVSiteMatch;
extern char *SVfile;

extern int RepeatMaxShift;/* Maximum alignment shift (in aligned sites) checked for repeat regions */
extern double RepeatPvalueRatio;/* Ratio of sizing error (Chi-Square) P-values for shifted alignment vs original alignment required to rule out repeat region */
extern int RepeatRec;/* Use improved RepeatMask implementation */
extern double RepeatLogPvalueRatio;/* minimum Ratio of Log10(Pvalue) for shifted alignment vs original alignment required to confirm repeat region */
extern double RepeatKbShiftRatio; /* minimum offset shift (kb) as fraction of alignment right end offset shift : this must be satisfied along entire alignment */
extern double RepeatKbShiftRatio2;/* maximum offset shift (kb) as fraction of alignment right end offset shift : this must be satisfied along entire alignment */
extern int RepeatMaxLabelDrop;/* If >= 0 : maximum drop in aligned labels for shifted alignment for repeats is Shift + RepeatMaxLabelDrop */
extern double RepeatLogPvalueRatioPerLabel;/* The minimum Ratio of "Log10(Pvalue) per label interval" for the shifted alignment vs the original alignment required to confirm repeat */

extern double RepeatShift;/* With pairwise alignment if > 0 : check for repeats in each Molecule by aligning it with itself with a minimum offset of RepeatShift */

extern int ScaleDelta;/* Scale Molecule size up/down by D = ScaleDeltaSize * bpp * (1 ... ScaleDelta) (by multiplying sizes by (1+D) and 1/(1+D)) */
extern double ScaleDeltaSize;
extern int ScaleDeltaBPP;/* If != 0 : Correctly handle bpp with ScaleDelta by treating Scaling as temporary (to find the alignment) and estimating parameters based on unscaled molecules */

extern int MapScale;/* If > 0 : Apply per map scaling correction to all -i input maps. Applied before each -M iteration based on the previous iteration's alignment (excluding outliers).
			    The first -M iteration is run twice using the first sub-iteration to apply a scaling correction */
extern int MapScaleDelta;/* replaces ScaleDelta during first sub-iteration of alignment if MapScale > 0 */
extern double MapScaleDeltaSize;/* replaced ScaleDeltaSize during first sub-iteration of alignment if MapScale > 0 */

extern int ForceOverwrite;/* 0 : Exit with error if output file already exists.
			     1 : Force overwrite of output files that already exist */

extern int Mfast;/* see FAST in refalign.cpp */
extern double MapRate;/* If > 0 : Limit Mapping rate to specified value by making -T more stringent and keeping only the best alignments at each iteration of -M */
extern double MaxLogPvThreshold, MaxScoreThreshold;/* Limit -T or -S adjustment to be no more stringent than these values */

extern int relaxBnds;/* see -relaxBnds */
extern int biaswtDefer;/* see -biaswtDefer */

extern int BestRef;/* only use the alignment with the best refid for each mapid in refalign and refalign2 */
extern int BestRefPV;/* If 1 : modify BestRef to pick the best refid based on logPV */
extern double BestRefPV_MaxConf;/* limit LogPV used to no more than this value when computing molecule weight (see -BestRefWT) */
extern double BestRefMargin;/* If BestRef : prefer earlier contig ids over later contig ids by this margin of score (or logPV, if BestRefPV) */
extern double BestRefExt;/* If > 0.0 : Only suppress extension alignments IF there is a center alignment (on a another ref contig NOT within BestRefExt (kb) of ends) with a log10(Pvalue) that
			   exceeds that of the extension alignment by BestRefExtMargin (can be -ve) */
extern double BestRefExtMargin;
extern double BestRefOutlier; /* Only suppress alignments with internal outliers sizing error above BestRefOutlier */
extern int BestRefWT;/* Alternative to BestRef : Only downweight molecule alignments based on relative likelihood (exp(score)) OR Pvalue */
extern int AddFragsWT;/* Modify BestRefWT by adding alignment weights of all alignments of the same query map with the same ref map in the same orientation, treating them as fragments of one large alignment */

extern int XMAPmapWT;/* see -XMAPmapWT : If 1 output mapWT column in .xmap right of alignment pairs */

extern int SkipIDcnt;
extern long long *SkipID;/* with -merge exclude specified set of map IDs from output */

extern int SelectIDcnt, SelectIDmax;
extern long long *SelectID;/* with -merge output a subset of maps based on the specified set of map IDs */

extern int SelectScanId;
extern int SelectScanId2;

extern int BreakIDinc;
extern int BreakIDcnt;
extern long long *BreakID;/* IDs to be broken into fragments are BreakID[i=0..BreakIDcnt-1] */
extern int *BreakPointCnt;
extern double **BreakPoint;/* Breakpoint locations are BreakPoint[i=0..BreakIDcnt-1][0..BreakPointCnt[i]-1] */

extern int WinBreak;/* If > 0 : break apart each input contig into overlapped fragments of size WinBreakLen (kb) each shifted by WinBreakShift (kb) from the previous one
		       OutputID = InputID * WinBreak + PieceID. Value of WinBreak will be at least the number of fragments of largest input contig */
extern double WinBreakLen;
extern double WinBreakShift;

extern long setMask;/* If >= 0 : set all map ends to have this mask value */
extern long orMask;/* If > 0 : OR value into all map end mask values */

extern int BestAlignments;/* If > 0 : output only the best BestAlignments alignments based on logPV */
extern int FirstAlignments;/* If > 0 : output the first FirstAlignments alignments that are above all 3 thresholds */ 

extern int PVres; /* If > 0 : Apply a correct to Pvalue to reflect limited resolution (-res & -bpp) */
extern int AlignRes;/* If > 0 : use more thorough -res based alignment (may result in 3x as many alignments and 10x slower refinement) */
extern double AlignResMult;/* If AlignRes > 0 : limit -res range to (res+SDrange*resSD)*AlignResMult */

extern double EndScoreThresh;/* If pairwise score without right end score (Send(),Sbnd()) is below bestscore (so far) by this threshold, skip calling Send() and Sbnd() */
extern double EndScoreThresh2; /* If pairwise score without Sbnd() is below bestscore (so far) by this threshold, skip calling Sbnd() */

extern int fast_mode; /* values >0 turn on various heuristic optimizations */
extern float fast_mode_tolerance;  /* Tolerance to use for pair oracle. This is in units of variance computed using sf/sd */

extern int XmapLen;/* If > 0 : Add query and ref map length columns to .xmap output */

extern int NoStat;/* If > 0 : Skip computation of statistics (requires ALIGN_COMPRESS) */

extern double XmapNoQueryOverlap;/* If > 0 : suppress alignments that overlap in the query map by this many kb with higher Pvalue alignments of the same query map, processed in descending order of logPV */

extern double XmapFilterDelta, XmapFilterThresh1, XmapFilterThresh2;/* NOTE:  XmapFilterDelta <= XmapFilterThresh2 <= XmapFilterThresh1 */
/* If XmapFilter* > 0.0 : Modify logPV threshold to XmapFilterThresh2 - XmapFilterDelta and then filter out final matchgroups as follows:
   1. Apply -XmapNonOverlapped (if != 0) while flagging each ref contig that had matchgroups removed with the highest confidence value for any matchgroup removed.
   2. If a ref contig has multiple (remaining) matchgroups, only keep the best one AND only if its logPV is better than 2nd best logPV by at least XmapFilterDelta, OR confidence  >= XmapFilterThresh1
   3. If a matchgroup had higher confidence matchgroups with the same ref removed (in step 1) only keep if its confidence is better than XmapFilterThresh1, 
      otherwise only keep if confidence is better than XmapFilterThresh2.
*/
extern int XmapUnique;/* If > 0 : suppress alignments for query maps that have less than this many sites (on query map) distinct from every higher Pvalue alignment of the same query map */

extern int XmapChim;/* If > 0 : check for chimeric query maps that align to multiple reference maps in .xmap output and have a minimum of this many disjoint sites in each region */
extern double XmapChimSpacing;/* If > 0.0 : Also check for two alignments to same reference map in .xmap output separated by at least this many kb */

extern int CovNorm;/* If > 0 : normalize coverage in _r.cmap and refined .cmap to reflect average underlying map size */
extern double CovLambda;/* If CovNorm : The mean map size of the underlying exponential distribution : actual mean map size = CovLambda + MinLen */

extern int align_format;/* 0 : use text format for .align output
			   1 : use binary format (except headers) for .align output
			   NOTE : .align input file will have a line "# binary = 1" IF the format is binary */

extern int CoverageIntRcmap;/* 1 : extend FRAGCOV_FIX >= 1 to output of _r.cmap in output_csites() */

/* HashTable parameters */
extern int HashGen;/* true if hashtable file is to be generated */
extern int HashWin;/* Window size (not including false/missing sites) */
extern int HashScore;/* hash score threshold : only save hashtable matches with score at least this large */
extern double SDmax;/* maximum Sizing Error of each interval as multiple of -sf -sd based values */
extern double SDrms;/* maximum RMS Sizing Error over entire window as multiple of -sf -sd based values */
extern double RelErr;/* maximum relative sizing error of each interval (if greater than SDmax based value) */
extern double OffsetKB;/* maximum alignment offset error in Kb */
extern int NumErrI,NumErrQ;/* maximum number of false/missing sites per alignment window for inserted and and probe map respectively */
extern int NumResQ;/* maximum number of unresolved sites per alignment window for probe maps */

extern double HashSF, HashSFscale;
extern double HashSD, HashSDscale;
extern double HashSR, HashSRscale;

extern int ShiftOffset1;/* see -hashGrouped */
extern int hashMaxCov;
extern int OffsetSpread;/* offset spread used when pairwise==0 : See OFFSET_SPREAD in hash.h */

extern int hashScaleDelta;/* see -hashScaleDelta */

extern int hashkeys;/* see -hashkeys */

extern int Hash_Bits; /* If > 0 : overrides HASH_BITS in hash.h */
extern int MHash_Bits; /* If > 0 : overrides MHASH_BITS in hash.h */
extern int HHash_Bits; /* If > 0 : overrides HHASH_BITS in hash.h */
extern int SubsetMaps; /* If > 0 && -ref : Output subset of query maps with HashTable matches in <prefix>_hash.bnx. With SplitHash > 0, output file name becomes <prefix>R<R>_hash.bnx */
extern int DeferMerge; /* If > 0 : Output a list of filenames that can be combined into the HashTable binary by Merge Sorting the files (which can then be overlapped with alignment) */
extern int HashInsertThreads;/* If > 1 :  Number of threads to use for multithreaded HashTable insertion */
extern int HashQueryThreads; /* If > 0 : Reduce number of threads from MaxThreads to this value for hashtable (except during HashTable insertion) */

extern int HashMultiMatch;/* Keep track of multiple hashtable matches per map pair, with offset differing by at least this value */
extern int HashMultiMatchMax;/* maximum number of hashtable matches per map pair considered (the best scoring matches are used) */
extern int HashMultiMatchMaxDelta;/* If > 0 : In each orientation for each map pair keep only hashtable matches within this hashtable score of the best hashtable score */

extern int HashTxt;/* 1 : output text version of HashTable file */
extern int HashMatrix;/* 1 : output matrix of matching intervals (including self matches) */

extern int HashBest;/* If >= 0 : For each query map only use hashtable matches within HashBest of the best hashscore */

extern int HashMaxHits;/* If > 0 : For each query map only use best HashMax hashtable matches : rounds up HashMax to include all matches with lowest score included */
extern int HashHitsCov;/* 2nd arg of -HashMaxHits */
extern int HashHitsGenSiz;/* 3rd arg of -HashMaxHits */


extern int HashGC;/* If > 0 : Reclaim memory of MatchTable every HashGC probe labels, for probe maps with more than this many labels. */
extern int HashT2;/* If > 0 : Use 2 threads per probe map (nested parallelism over alignment orientation) */

extern int ChimQuality;/* output Chimeric junction Quality Score at each site of Consensus after refinement :
			  score ranges from 0 to 100 and low scores mean more likely to be Chimeric junction at that site */

extern int noQuality;/* disable quality score during refinement (for backward compatibility) */

extern void output_err(char * vfixx_filename, Cparameters *parameters, int iterations, double *MappingRate, double CoverageMult, int MIN_PV, int MAX_PV);
extern void input_err(char *filename, Cparameters *parameters);
extern int output_filename_matched(const char *filename);
extern void swapparam(int usecolor, Cparameters *parameters = NULL, int iterations = 0);

#endif
