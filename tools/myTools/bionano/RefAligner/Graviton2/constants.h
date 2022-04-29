#ifndef CONSTANTS_H
#define CONSTANTS_H

#define VERB 1 /* 1 : verbose output (OK for normal use)
		  2 : very verbose output (for debugging only) 
	       */
#ifndef DEBUG
#define DEBUG 1 /* 0 : no assertion checking
		   1 : assertion checking (under 10% slowdown)
		   2 : more assertion checking (400% slowdonw) */
#endif

#define PAGE 4096 // size of page in bytes
#define HUGEPAGE /* (512*PAGE) */ 0 // size of 2M huge page (0 to disable transparent huge page madvise)


#ifndef USE_MMAP
#define USE_MMAP 0 // Means use huge_pages via mmap (useful for MIC, but not tested)
#endif

#define OMP_DEBUG 0 // (USE_MIC ? (RELEASE ? 0  : 1) : 0)

#ifndef RELEASE
#define RELEASE 0 /* If 1 : turns off some risky assertions. Should be 1 for customer release binaries, 0 in all other cases */
#endif


#ifndef LEAKDEBUG
#ifdef VALGRIND
#include <valgrind/memcheck.h>
#define LEAKDEBUG 1
#else
#ifdef ASAN
#define LEAKDEBUG 1
#else
#define LEAKDEBUG 0 /* 0 : don't bother freeing memory before end of program 
					      1 : free memory at end of of program (to help catch memory leaks in valgrind or with ASAN) */
#endif
#endif
#endif

#define MSAN false
#define TSAN false

#ifdef CLANG

#if defined(__has_feature)
#  if __has_feature(memory_sanitizer)  // code that builds only under MemorySanitizer
#include <sanitizer/msan_interface.h> 
#undef MSAN
#define MSAN true
#  endif
#endif

#if defined(__has_feature)
#  if __has_feature(thread_sanitizer) // code that builds only under ThreadSanitizer
#undef TSAN
#define TSAN true
#  endif
#  endif

#endif // CLANG


#define MAX_KB_BINS 20 /* maximum KB bins used by -extsplitOutlier */

#define CALIGN_END 0 /* -1 : Lij1,Rij1,Lij2,Rij2 are adjusted for Endoutliers to reflect last aligned label.
			 0 : Add LijY,RijY,LijX,RijX to correspond to values of Lij1,Rij1,Lij2,Rij2 before adjustment
			 1 : Lij1,Rij1,Lij2,Rij2 are NOT adjusted for Endoutliers and always reflect last overlapped label : LijY etc are not defined. (NOT YET IMPLEMENTED) */

#define MULTIMATCHES_FIX 1 /* Consider 2 matchgroups distinct if they are non-overlapping at both ends of query or reference (Don't require that the matchgroups have MultiMatch non-overlapped labels, so -MultiMatch 30 -A 9 can be used) */

#define FRAGCOV_FIX 2 /* TRY 1 or 2 for Solve 3.5 */ /* >= 1 : base fragcnt and fragcntT on label interval instead of label
							>= 2 : Don't extend Lij,Rij (last labels overlapped by a molecule) to include virtual labels with Hdel[i] != 0 */
#define FRAGCOV_FIX_CSITES CoverageIntRcmap /* extend FRAGCOV_FIX >= 1 to output_csites() */

#define COVERAGE_FIX 1 /* WAS47 0 */ /* In output_csites() and output_draft() change Coverage value to include ends of query maps beyond first/last aligned label */

#define PAIRMERGE_HUGE 9999999.0 /* largest supported pairmerge overlap (in kb) */

#define REFSPLIT_FIX 2 /* >= 1 : fixed matchgroup overlap code to apply RefSplit to split each matchgroup to avoid discarding 2 good sub-submatchgroups seperated by a very -ve scoring outlier region
			  >= 2 : Also split matchgroups (see -Refsplit) and merge -ve scoring regions into a single outlier after applying -MultiMatchesFilter */

#define FRATE_FIX0 1 // fix bug relating to Frate penalty in refalign() when using OutlierType0 (for BNX query maps) : X * Frate penalty is now part of Gaussian term (NOT misaligned label penalty) to avoid including it in outlier score
#define FRATE_FIX1 1 // fix bug relating to Frate penalty in refalign() when using OutlierType1 (for CMAP query maps) : X * Frate penalty is now part of Gaussian term (NOT misaligned label penalty) to avoid including it in outlier score. (NOTE : Always excluded from Gaussian term when computing Indel score)
#define FRATE_FIX 1 // fix bug relating to Frate penalty in refine() when using OutlierTypeRef (NOTE : sometimes results in too many indels and fewer SNPs ?)


#define PVALUE_FIX 0 /* 1 :  Force LogPvalue to 0 if first or last aligned interval is an outlier
			2 :  Just penalize LogPvalue if first or last aligned interval is an outlier by treating the outlier misaligned labels as misaligned labels of
			     the non-outlier region (rather than the outlier) */

#define OUTLIER_CONF_FIX 2 /* WAS49 1 */ /* include score of unaligned end, if they are not endoutliers, in end confidence value for outliers/indels, provided there are at least
					    OUTLIER_CONF_FIX aligned labels at that end */

#define DAG_CHECK 0 /* check for cycles in coverage update DAG (for debugging : see cov_update() in Assembler_graph.cpp) */

#define END_NOEXT 0x1 // Used in Cmap::Mask[0][1,N+1] field to flag Contigs Ends that are not to be extended (during extension refinement) due to E&S split (-CloneLimitLen or -Truncate) or Palindromic trimming or chr end detection 
#define END_NOEXT2 0x2 // Used in Cmap::Mask[0][1,N+1] field to save value of END_NOEXT flag during -SplitSegDup phase of pairmerge 
#define END_CHR    0x20 // Used in Cmap::Mask[0][1,N+1] field to mark chromosome end (typically also marked with END_NOEXT) to block all merges with -TrimmedNoMerge, regardless of overlap size

#define SEQID 0 /* replace id numbers in vfixx file by sequential order number  : useful if vfixx files are concatenated and have duplicate id's 
		   This always happens if there are any duplicate ids or id values of 0 */

#define MAX_CHIM 16 /* maximum number of ChimQuality window sizes (see -ChimQuality 3 WT Win1 Win2 Win3 Win4 ... ) */

#define MAXFILES (USE_MIC ? 1024*1024 : 1024*1024*8) /* maximum number of -i or -a input files */

#define SIMPLE_BIAS 2 /* WAS 1 */ /* Type of -resbias correction:
				     0 : old log-linear fit : Cannot be debugged using REFDEBUG
				     1 : simple piece-wise constant fit
				     2 : piece-wise linear fit */
#define BIAS_TRACE 0 /* trace BiasCorrect and Scaling correction of map with id = BIAS_TRACE_ID */
#define BIAS_TRACE_ID 7585LL

#define SCORE_APPROX 1 /* >=1 : approximate Yc() so deltas are exact
		          >=2 : approximate PM() to use only a single call to Pn() (and gradients are exact) */
#define KMAX_FIX 0 /* TRY 1 (30% slower ?) */ /* limit Y[I]-Y[I-K] by AlignResMult * (res + 3 * resSD) but NOT also Y[I]-Y[I-1], .. Y[I-K+1]-Y[I-K] by (res + 3*resSD) */

#define OUTLIER_PV 1  /* enable Outlier based adjustment of Pvalue */

#define OUTLIER_MAX 1  /* enable command line option -outlierMAX (slight slowdown) */

#define ISLINEAR 0 /* enable -islinear command line arg */

#define ROC 0 /* collect statistics of all alignments for ROC curve (enable only for testing, since it wastes memory and time)
	         based on simulated data (or for refalign, assumes best score above threshold == true alignment) */

#define RESBINS (64+7) /* maximum number of piecewise-linear segments in bias estimation, make sure actual arrays are 64-byte aligned */

#define GAMMA Gbias /* see parameters.h */


#define RMS_MIN 0.05 /* flag alignments with average normalized RMS error below this, provided number of aligned sites >= RMS_SITES (to disable set to -ve value) */
#define RMS_SITES 5 /* minimum aligned sites, before applying this test */
#define RMS_MIS 2 /* maximum misaligned sites, before applying this test */
#define RMS_PVALUE 1e-8 /* maximum Pvalue of ChiSq distribution */

#define LINESIZ (512*1024) /* maximum length of input lines */
#define MAXSITES (32*1024) /* maximum sites per input line */

#define MAXCOLOR 2 /* maximum number of colors supported */

#define MINSITE 63 // (USE_MIC ? 63 : 255) /* pre-allocate memory for this many sites : should be larger than number of sites in most BNX maps */
#define DELTA_MAX 64 /* maximum value of DELTA_X or DELTA_Y */

#define K_MAX 15 /* maximum value of KMAX */

extern int DELTA_Y; /* maximum span of aligned sites in reference */
extern int DELTA_X; /* maximum span of aligned sites in Nanomap */
extern int DELTA;  /* maximum span of aligned sites for pairwise alignment */
extern int KMAX;   /* maximum number of consecutive references sites to merge due to resolution limit is KMAX+1 */

#define DELTA_YEND (10*DELTA_Y) /* maximum Y sites that can overlap DELTA_X intervals on X (in Sbnd()) */

#define QUADRATIC_VARIANCE  1 /* add quadatic term to variance */
#define RES_VARIANCE 1 /* add resolution term to variance : Does not seem to help */

#define REFSCORETYPE 4 /* 1 Valouev's score
			  2 Thomas's original score
			  3 Thomas's refactored score
			  4 Thomas's 2nd refactored score */

#define RESSD 0 /* 1 : use older simplified resSD based score (only for K>=2 in Sm() if SCORE_APPROX <= 1) */

#define PAIRSCORETYPE 5 /* 1 Thomas's original score
			   2 Nguyen's score
			   3 John's refactored score
			   4 Thomas's refactored score
			   5 thomas's 2nd refactored score
			*/

#define END_CHIMQUAL -1.0 /* Value of ChimQuality at ends of contigs where true value is unknown */

#define MINOVERLAP 0.38 /* minimum actual overlap to count as possible true positive for normalizing ROC curves */

#define ENDEXTEND MININTERVAL /* Extend ends of maps from .bnx, .vfixx or .cmap file past end of first/last site by at least this amount (With .bnx input If bnxerror > 0, site beyond the ends are not allowed, but sites
				 at the end are allowed and these ends are extended by ENDEXTEND). For various reason this value must be > 0 */

#define MINSCORE -1.0e+30 /* invalid alignment score (can be used as float OR double) */

/* fields widths for BNX 1.0 molecule id */
#define MOLID_WIDTH 1000000LL
#define SCANNUMBER_WIDTH 1000LL
#define FLOWCELL_WIDTH 10LL
#define CHIPID_MAX 922000000LL

/* field widths for vfixx molecule id (see test_generator.cpp) */
#define CHR_WIDTH 100LL /* decimal field width of map chromosome number (in modulo notation) */
#define LOC_WIDTH 1000000000LL /* decimal field width of map Location (in modulo notation) */
#define LEN_WIDTH 10000000LL   /* decimal field width of map size (in modulo notation) */
#define FLAG_WIDTH 10LL /* decimal field width of map flags (3 boolean flags per decimal) */

#define CHR_OFFSET 1000000000LL /* bp location offset as multiple of chromosome number */

/* bounds on error parameters */
#define MINFP MinFP
#define MAXFP MaxFP

#define MINFN MinFN
#define MAXFN MaxFN

#define MINSD (QUADRATIC_VARIANCE ? MinSD : 0.001)
#define MAXSD MaxSD

#define MINSF MinSF
#define MAXSF (QUADRATIC_VARIANCE ? MaxSF : 0.25)

#define MINSR (QUADRATIC_VARIANCE ? MinSR : 0.00)
#define MAXSR (QUADRATIC_VARIANCE ? MaxSR : 0.00)

#define MINSE (RES_VARIANCE ? MinSE : 0.00)
#define MAXSE (RES_VARIANCE ? MaxSE : 0.00)


#define MINSD_MULT MinSD_Mult /* WAS 0.99 */ /* SD must be at least equal to -sqrt(2*SF*SR) * MINSD_MULT */

#define ALIGN_COMPRESS 1 /* WAS 0 */ /* compress alignments[] array used by refalign[] */

#define LARGE_NEGATIVE   (-1e35)
#define LARGE_POSITIVE    (1e35)
		     
#include "ident.h" 

static Ident constants_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/constants.h 11165 2020-06-12 18:38:43Z tanantharaman $");

#endif
