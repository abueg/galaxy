#ifndef REFINE_H
#define REFINE_H

static Ident refine_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/refine.h 11414 2020-08-06 06:18:41Z tanantharaman $");

#define FLOAT_EXPLOG 0 /* try USE_MFLOAT */ /* Use float precision for log(1.0 + exp(minLP-maxLP)) */
#define HAP_VITERBI 1 /* try 1 */ /* Use Viterbi approxiation for combining likelihoods for 2 Alleles : less rebust, but more sensistive to small indels */

#if FLOAT_EXPLOG==0
static double HAP_VITERBI_WT = 1.0;
#else
static float HAP_VITERBI_WT = 1.0;
#endif

/* combine Log likelihood ratio for molecule with 2 Alleles : logH == log(0.5) */
//#pragma omp declare simd uniform(logH)
static inline double HapLP(const double LP1, const double LP2, const double logH)
{
#if HAP_VITERBI == 1
  double maxLP = max(LP1,LP2);
  double minLP = min(LP1,LP2);
#if FLOAT_EXPLOG==0
  return maxLP + HAP_VITERBI_WT * (logH + log(1.0 + exp(minLP - maxLP)));
#else
  return maxLP + HAP_VITERBI_WT * (logH + logF(1.0f + expF(minLP - maxLP)));
#endif

#else // !HAP_VITERBI

  double maxLP = max(LP1,LP2);
  double minLP = min(LP1,LP2);
#if FLOAT_EXPLOG==0
  return maxLP + logH + log(1.0 + exp(minLP - maxLP));
#else
  return maxLP + logH + logF(1.0f + expF(minLP - maxLP));
#endif

#endif // !HAP_VITERBI
}

#define ANNEAL_SITE_PEN 1 /* During extension refinement site add phase use Site_Pen == 0 */

#define REMOVE_CLOSE_SITES 1 /* After refinement remove nearby sites closer than rres * 0.500 */

#define MEDIAN_SKIPOUTLIER 1 /* When computing median sizes ignore outliers if !oldmap (old default was 1, and is the same as 0 IF WEIGHTED_MEAN < 0) */

#define LP_MINDELTA 1e-8 /* minimum LP improvement threshold : improves repeatability if this exceeds the random LP variation due to multithreading based order of maps */

#define TRIM_FIX 1 /* fix sitecntFN[] and sitecntFNnorm[] computation to base COVERAGE_TRIM_LEN distances on reference Y[] (NOT map X[]) */

#define FAST_RESIZE 1 /* If MULTIMODE, don't use non-Mprobeval SizesEstimate() (FAST_RESIZE==0 is no longer fully supported) */

#define FILL_MAP RefineFillMap /* In UpdateMap, fill mapping map[m][] from Xm[] to Hcuts[] by extrapolation (beyond last aligned sites at each end),
				  to avoid unstable site deletion, when only a single site in Xm[] remains mapped to Y[] and that site is deleted.
				  Don't do this during final UpdateMap (when coverage profiles are computed) to avoid aligning sites in Xm[] to Hcuts[i] with Hdel[i]==1 */
#define FILL_HMAP RefineFillHmap /* Value of FILL_MAP when UpdateMap() is called from HaploType() */

#define NODELETE 0 /* disable site deletion of original sites */

#define ADAPTIVE 2 /* 1 : Use adaptive expansion of Goldmean search range (as long as range is >= 0.100) when maximum appears to have drifted outside the range 
		      2 : Apply for any size range (slower) (should use a gradient based search that is faster when range is small enough to be unimodal, eg 0.100) */
#define GM_ITER ((DELTA_RES < 0.009) ? (DELTA_RES < 0.0009 ? 500 : 300) : 100) /* maximum Goldenmean search iterations */

#define QUICKUPDATE 1 /* Update fixYdel[] as soon as possible (0 no longer supported) */
#define DELTASLOW 1 /* 2x slower mprobeval() based sizing changes, but perhaps more stable (only applies if QUICKUPDATE==0) */

#define SLOWMULTIMODE MprobevalTest  /* 4x slower MultiMode code to better find global maximum (usually to diagnose if regular mode is good enough) */

#define FRAGCNT_FIX 1 /* count coverage only between aligned sites (exluding non-outlier ends) */
#define FRAGCNT_FIX2 1 /* for K to zero when nmapMD[I] < 0 so quality score convention matches normal UpdateMap() convention */

#define NEW 1 /* also consider end outliers that are not near the ends of X */

#define SITEMERGE SiteMerge /* try to merge nearby sites in consensus (slow) (see -SiteMerge) */

int probeval::RANGE = 8; /* TRY smaller value 2,3,4 after setting RANGE_UPDATE. Reducing RANGE increase the risk of local maxima when scoring changes with mprobeval()
			    since alternate alignments that differ by more than RANGE labels are not considered. Consider using a larger value initially
			    (see -RefineRange for setting this value on command line) */
#define RANGE_Y1 RANGE /* RANGE_Y value used during site interval updates */
#define RANGE_Y2 0           /* RANGE_Y value used during site add/delete */
int probeval::RANGE_Y = RANGE_Y1; /* maximum change (shift) in alignment allowed in Y & X sites (if previous alignment is known)
				  If RANGE_Y=0, instead change in X & Y sites is based on corresponding distance of RANGE sites on X (but see FIX_RANGE_Y)
				  NOTE for mprobevalp() : RANGE_Y == RANGE is most accurate with site interval changes in mprobevalwinSize(), 
				  while RANGE_Y=0 is most accurate for site additions or deletions in mprobevalwin() */

#define RANGE_UPDATE RangeSwitch /* If != 0 : larger value of RANGE to use when updating alignments for use by UpdateMap() */

/* #define OUTLIER_TYPE_SWITCH OutlierTypeSwitch */ /* switch to OutlierType = 1 before mprobeval() with site add/delete and 
						 switch to OutlierType = 0 before mprobeval() with interval size changes */
/* #define OUTLIER_LAMBDA_SWITCH outlierLambdaSwitch */ /* switch to outlierLambda = 9999 before mprobeval with site add/delete */

#define SLOWADDDELETE 0 /* TEST 1 : why is 1 worse in some cases ? */

#define LRANGE_INIT (SLOWADDDELETE ? 1000 : 10)
int probeval::Lrange = LRANGE_INIT;
#define LRANGE_HINIT (SLOWADDDELETE ? 1000 : 10) /* Lrange is set to this value during haplotype detection (to reduce change of invalid LP) */
#define LRANGE_FIX 2            /* 0 : change in Y is assumed to only affect LR of maps that overlap within LRANGE*Xtheta of the previous overlap range
				   1 : Change in Y is assumed to only affect LR of maps that overlap within LRANGE sites of Y from the previous overlap range.
				   2 : Use entire overlap range (as determined by RANGE & RANGE_Y & FIX_RANGE_Y and EXTRAPOLATE=0) before extending it by LRANGE sites on Y */

#define SCANRANGE ScanRange /* WAS98 (SLOWADDDELETE ? 1000 : 3) */ /* multiple of Xlambda to scan ahead/back for improvements in Y before confirming them (to pick locally best improvements) */
#define SPLITRANGE 100 /* split mprobevalwin() when site ranges with improvements are separated by SPLITRANGE*Xtheta or more */
#define EXPANDSPLIT2 20 /* During split, expand each split range by EXPANDSPLIT2*Xtheta from sites that improve previous value (even if left unchanged), to reflect dependencies between changes 
			   NOTE : EXPANDSPLIT2+EXPANDSPLIT should be less than (SPLITRANGE/2) */
#define EXPANDSPLIT 10 /* expand any split range by EXPANDSPLIT*Xtheta from sites that were changed, to reflect dependencies between changes */

#define LOCAL_THETA LocalTheta /* WAS 0 */ /* Use local estimate of consensus label interval size instead of Xlambda or Xtheta with SCANRANGE,SPLITRANGE,EXPANDSPLIT,EXPANDSPLIT2 */

#define HLP_MINDELTA 1e-6 /* WAS 1e-7 */ /* minimum LP improvement threshold during haplotyping : improves repeatability if this exceeds the random LP variation due to multithreading based order of maps */

#define HSCANRANGE HScanRange /* WAS98 3 */ /* multiple of Xlamda to scan ahead/back for improvements in Haplotypes before confirming them (to pick locally best improvements) */
#define HINDEL (HapIndelPvalue ? 1 : 0) /* Check for indels (NOTE: only fixed values in Initial_Delta[0..INITIAL_DELTA_RANGE-1] are checked) */
#define HSIZES (HapIndelPvalue ? 3 : 0) /* Optimize interval sizes :
					   >= 1 : Only applies until HapSites[] has no changes during a single iteration
					   >= 2 : Applies, until HapDelta[] has no suggested changes (see LP_INDEL_MINDELTA)
					   >= 3 : Applies, until also Delta[] has no suggested changes (see LP_INTERVAL_MINDELTA)
					          NOTE : 3 behaves same as 2 during 1st stage if HINDEL_PHASED != 0
						         4 applies unconditionally to HINDEL_PHASED
					  
					*/
#define HINDEL_NONSYMETRIC 1 /* WAS 0 */ /* also check for non-symetric indels (indel + size change combinations) */
#define HINDEL_UPDATE 1 /* Update non-zero indel size even during initial scan (rather than only during indel size optimization stage) */

#define HINDEL_MERGE (HapIndelMerge ? 2 : 0)  /* 1 : Enable various improvements after revision r4535 to fix bugs or support larger indels (not clear if this works without 2)
						 2 : Treat multiple Haplotype Indels between two consecutive non-haplotype labels as a single Haplotype, but just leave them
						     distributed any way that does not generate negative sub-intervals, but merge them on output into the leftmore sub-interval.
					      */
/* NOTE : the following options can work without HINDEL_MERGE, but this would break backwards compatibility, so they are controlled by -RefineFix (defaults to 0) */
#define HINDEL_SPLIT_FIX (RefineFix >= 2 ? 1 : 0) /* ignore risk of splitting HapIndel by a new site but rewind index if LP drops : Does not always help */
#define DUP_FIX RefineFix /* suppress adding Hap Sites on opposite Alleles at the same location (should instead add non-Hap site) */
#define HMAP_SWAP RefineFix /* try to swap alignments between alleles for each map before and while trying to flip phases + */
#define DELTA_EXPAND RefineFix /* expand Delta values range when switching to different Pvalues in phase 2 */
#define HMAP_TMPFIX RefineFix /*  Make the fixes to map1,map2 in hprobeval using temporary memory (otherwise the changes are permanent, which results in some Consensus changes not being reversible) */
#define HPROBEVAL_SPREAD 1 /* WAS RefineFix */ /* >= 1 : Before calling hprobeval() hsetmap() spreads out all Delta[] and HapDelta[] values between actual labels into subintervals of Hcuts[]
						  so hprobeval() sees the same Hcuts[] that mprobeval() sees */
#define OUTLIERREF_FIX RefineFix /* New Poutlier selection based on -outlierRef instead of matching -endoutlierRef then reverting to original value
						      NOTE : If -outlierRef is not specified, the only effect of turning this option off is to force the last -outlierRef value to the regular
						      outlier value and to force all the outlierRef values to be at least as stringent as the regular outlier value */
#define MPROBEVAL_FIX RefineFix /* adjustment to initial scan of sizes per interval in mprobevalwinResize() ++ */
#define HMAP_MINDISTANCE (RefineFix ? 0.001 : 0.0) /* If != 0 : In repositionH, the minimum distance between any pair of labels in Hcuts[] that are present (via HapSite) in the same Allele */
#define REPOSITIONH_FIX RefineFix /* fix some bugs in repositionH */


#define HPROBEVAL_BALANCE 0 /* If != 0 : For any given SNP change or new Indel check how many molecules increase vs decrease the score (significantly).
			       Then apply the same heuristic as HapMinCov to reject the new Indel. This is more practive than apply HapMinCov at the end but
			       also more risky since it may reject a Hap Indel prematurely (NOT YET IMPLEMENTED) */

#define HLP_MAXDROP 0.20 /* With HINDEL_MERGE : If LP before Indel iterations drops more than this much (per map)  below its high water mark (of current stage), terminate current stage of HaploTypeing
			    This is a failsafe in case something goes very wrong. Some drop in LP may happen due to the changing Scoring parameters used in alternate iterations */
#define HLP_DROP_CNT 5 /* With HINDEL_MERGE : limit on number of Indel iterations that LP was below HWM (drop count), before terminating current stage of HaploTyping (and reseting both HWM and drop count) */

#define HMIN_ITER 25 /* HERE HERE HERE make commandline variable */ /* first iteration in which we consider indel/intervals (0 based) */
#define HINDELITER(iter) ((HINDEL || HSIZES) && iter >= HMIN_ITER && ((iter)%2)==1) /* If 0 : Current iteration handles SNPs only
										       If 1 : Current iteration handles indels and interval sizes only */
#define HSNP_FAST 2 /* Enable faster way of computing predicted change in LP from changing SNPs (HapSite[]) :
		       1 : Do it both ways and compare the results (to DEBUG)
		       2 : Only do it the fast way */
#define HINDEL_FAST 2 /* Enable faster way of computing predicted change in LP from changing sizes (Delta[],HapDelta[]) : 
					   1 : Do it both ways and compare the results (to DEBUG)
					   2 : Only do it the fast way */

#define HINDEL_PHASED (HINDEL_MERGE ? 0.001 : 0.01) /* Perform Indels in 2 phases : with HINDEL_PHASED then original HapIndelPvalue (in stages)
						       to get more precise estimate of indel size before performing Pvalue check with original HapIndelPvalue (SLOW if Pvalue is larger than 0.01) */

#define HSNP_PHASED (HINDEL_MERGE ? 0.1 : 0.011) /* Use max(HSNP_PHASED,sqrt(HapSitePvalue)) during 1st phase (see HINDEL_PHASED) */

#define HINDEL_ANNEAL (HINDEL_MERGE ? 0.01 : 0.0) /* If != 0 : Maximum change factor in HapIndelPvalue every other iteration, to avoid deleting too many indels simultaneously */


#define HINDEL_REVERSE 1000 /* TRY (HMIN_ITER + 2*(LP_INDEL_SKIP + MAX_DELTA_OVERSAMPLE-1) + 16) */ /* try to reverse indels no later than this iteration, even in Phase 1 */


#define ITER_INC HapStageIter /* fixed number of additional iterations per Haplotyping phase */
#define ITER_FRAC HapStageFrac /* variable number of additional iterations per Haplotypeing phase as multiple of max(N1,N2) */

extern int HREVERSAL; /* >=1 : Try to reverse the phase of all Haplotypes to the right of Hcuts[i] : This is repeated HREVERSAL times */
extern int FINAL_HREVERSAL;/* Repeat phase reversal one more time before terminating Haplotyping loop (after applying all DERES and FILTER stages) : needed to allow filtering of unphased SNPs */
extern int HAP_MAX_DERES; /* maximum number of times to apply -rres -cres during Haplotypeing */
extern int DERES_STOP; /* don't resume Haplotyping after final DERES stage */
extern int FILTER_STOP; /* don't resume Haplotyping after Filter & Trim stage */

/* indel values (kb) during initial haplotype scan (to detect if an interval has a haplotype indel) */
//#define INITIAL_DELTA_RANGE 12 // binary ladder
//static double Max_Initial_Delta[INITIAL_DELTA_RANGE] = {0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4};

//#define INITIAL_DELTA_RANGE 24 // 41.4% ladder
//static double Max_Initial_Delta[INITIAL_DELTA_RANGE] = {0.05, 0.0707, 0.1, 0.1414, 0.2, 0.2828, 0.4, 0.5657, 0.8, 1.1314, 1.6, 2.2627, 3.2, 4.5255, 6.4, 9.0510, 12.8, 18.1019, 25.6, 36.2039, 51.2, 72.4077, 102.4, 144.8155};

//#define INITIAL_DELTA_RANGE 48 // 18.9% ladder
//static double Max_Initial_Delta[INITIAL_DELTA_RANGE] = {0.0707, 0.0837, 0.1, 0.1189, 0.14, 0.17, 0.2, 0.24, 0.28, 0.34, 0.4, 0.48, 0.56, 0.68, 0.8, 0.96, 1.12, 1.36, 1.6, 1.92, 2.24, 2.72, 3.2, 3.84, 4.48, 5.44, 6.4, 7.68, 8.96, 10.88, 12.8, 15.36, 17.92, 21.76, 25.6, 30.72, 35.84, 43.52, 51.2, 61.44, 71.68, 87.04, 102.4, 122.88, 143.36, 174.08};

// #define INITIAL_DELTA_RANGE 96 // 9.1% ladder
// static double Max_Initial_Delta[INITIAL_DELTA_RANGE] = {0.07, 0.0775, 0.085, 0.0925, 0.1, 0.11, 0.12, 0.13, 0.14, 0.155, 0.17, 0.185, 0.2, 0.22, 0.24, 0.26, 0.28, 0.31, 0.34, 0.37, 0.4, 0.44, 0.48, 0.52, 0.56, 0.62, 0.68, 0.74, 0.8, 0.88, 0.96, 1.04, 1.12, 1.24, 1.36, 1.48, 1.6, 1.76, 1.92, 2.08, 2.24, 2.48, 2.72, 2.96, 3.2, 3.52, 3.84, 4.16, 4.48, 4.96, 5.44, 5.92, 6.4, 7.04, 7.68, 8.32, 8.96, 9.92, 10.88, 11.84, 12.8, 14.08, 15.36, 16.64, 17.92, 19.84, 21.76, 23.68, 25.6, 28.16, 30.72, 33.28, 35.84, 39.68, 43.52, 47.36, 51.2, 56.32, 61.44, 66.56, 71.68, 79.36, 87.04, 94.72, 102.4, 112.64, 122.88, 133.12, 143.36, 158.72, 174.08};

//#define INITIAL_DELTA_RANGE 192 // 4.4% ladder
// static double Max_Initial_Delta[INITIAL_DELTA_RANGE] = {0.0707, 0.0738, 0.0771, 0.0805, 0.0841, 0.0878, 0.0917, 0.0958, 0.1, 0.1044, 0.1091, 0.1139, 0.1189, 0.1242, 0.1297, 0.1354, 0.1414, 0.1477, 0.1542, 0.1610, 0.1682, 0.1756, 0.1834, 0.1915, 0.2, 0.2089, 0.2181, 0.2278, 0.2378, 0.2484, 0.2594, 0.2709, 0.2828, 0.2954, 0.3084, 0.3221, 0.3364, 0.3513, 0.3668, 0.3830, 0.4, 0.4177, 0.4362, 0.4555, 0.4757, 0.4967, 0.5187, 0.5417, 0.5657, 0.5907, 0.6169, 0.6442, 0.6727, 0.7025, 0.7336, 0.7661, 0.8, 0.8354, 0.8724, 0.9110, 0.9514, 0.9935, 1.0375, 1.0834, 1.1314, 1.1815, 1.2338, 1.2884, 1.3454, 1.4050, 1.4672, 1.5322, 1.6, 1.6708, 1.7448, 1.8221, 1.9027, 1.9870, 2.0749, 2.1668, 2.2627, 2.3629, 2.4675, 2.5768, 2.6909, 2.8100, 2.9344, 3.0643, 3.2, 3.3417, 3.4896, 3.6441, 3.8055, 3.9739, 4.1499, 4.3336, 4.5255, 4.7258, 4.9351, 5.1536, 5.3817, 5.6200, 5.8688, 6.1287, 6.4, 6.6834, 6.9792, 7.2882, 7.6109, 7.9479, 8.2998, 8.6672, 9.0510, 9.4517, 9.8701, 10.3071, 10.7635, 11.24, 11.7377, 12.2573, 12.8, 13.3667, 13.9585, 14.5765, 15.2219, 15.8958, 16.5995, 17.3345, 18.1019, 18.9034, 19.7403, 20.6143, 21.5269, 22.48, 23.4753, 24.5146, 25.6, 26.7334, 27.9170, 29.1530, 30.4437, 31.7916, 33.1991, 34.6689, 36.2039, 37.8067, 39.4806, 41.2286, 43.0539, 44.96, 46.9506, 49.0293, 51.2, 53.4668, 55.8340, 58.3060, 60.8874, 63.5831, 66.3982, 69.3379, 72.4077, 75.6135, 78.9612, 82.4571, 86.1078, 89.9201, 93.9012, 98.0586, 102.4, 106.9336, 111.6680, 116.6120, 121.7748, 127.1662, 132.7964, 138.6758, 144.8155, 151.2270, 157.9224, 164.9142, 172.2156};

/* 3 if INITTIAL_DELTA_RANGE==181, 2 if INITIAL_DELTA_RANGE==91, 1 if INITIAL_DELTA_RANGE==46, 0 if INITIAL_DELTA_RANGE <= 22 */


static double *Max_Initial_Delta = NULL; 
// extern int MAX_DELTA_OVERSAMPLE;/* how many times the number of samples can be reduced 2x after initial Indel scan has been completed (see LP_INDEL_SKIP) */
#define INITIAL_DELTA_RANGE ((24 << MAX_DELTA_OVERSAMPLE) + 1)
int probeval::localINITIAL_DELTA_RANGE;/* see HaploTypeR() : <= ((24 << DELTA_OVERSAMPLE) + 1) <= INITIAL_DELTA_RANGE */
// extern double MaxDelta; // largest HapDelta value to use during Haplotyping, defaults to 175kb

#define START_DELTA (FastIndel ? 0.050 : 0.100) // WAS50 0.050 /* value of Max_Initial_Delta[0] */
#define MIN_DELTA 0.0016 // WAS50 (/*NEW41*/DELTA_OVERSAMPLE ? START_DELTA * 0.5 : 0.0012) /* smallest value of Initial_Delta[IN][0] */
#define MAX_DELTA 0.100 // WAS50 (/*NEW41*/DELTA_OVERSAMPLE ? START_DELTA : 0.100) /* largest value of Initial_Delta[IN][0] */

#define LOW_VALUES (5 << DELTA_OVERSAMPLE) /* Keep exactly this many values smaller than the last active value (subject to MIN_DELTA & MAX_DELTA)  */
#define HIGH_VALUES (5 << DELTA_OVERSAMPLE) /* If DELTA_OVERSAMPLE <= 0 : Use (at most) this many values larger than the last active value, ignoring any larger values */

#define MIN_INDEL_SIZE 1e-6 /* WAS 0.001 */ /* Count Indels (and add HapIndelPvalue based prior penalty) only if their size is greater than specified size. NOTE : MUST BE > 0 */

#define LP_INTERVAL_MINDELTA IntervalEps /* adjust INITIAL_DELTA for interval sizing only if potential LP improvement is at least this large. Also terminate if no improvements this large were noticed */

#define LP_INDEL_MINDELTA  ((HINDEL_MERGE==0 || HapIndelPvalue == origHapIndelPvalue || DELTA_STOP) ? HapIndelEps : HapIndelEps2) 
                           /* adjust INITIAL_DELTA for indels only if potential LP improvement is at least this large. Also terminate if no improvements this large were noticed */

#define LP_SITE_MINDELTA HapLabelEps /* Terminate if no improvements this large were noticed for site addition/deletion */

static double LP_SKIP = 100.0; /* Values set by -skip_score : If adding/deleting site i causes LP to drop by this amount, skip subsequent checks of site i */

#define FAST_SKIP 0 // ((SLOWADDDELETE || HapSitePvalue > 0.0 || extendonly) ? 0 : 1) /* TRY 0 during extension refinement : slower but more accurate */ /* Once a site addition is skipped (due to LP_SKIP) don't recheck that site again even after sizing changes or Pendoutlier changes. This resets after unaligned sites are re-positioned (once for each -endoutlierRef value) */

#define INIT_TRIM TrimFactor /* initial trim factor, may be set more aggressively */
#define INIT2_TRIM TrimFactor /* 2nd trim factor, default trim factor unless lower values score better */
#define EARLY_MINTRIM TrimFactor*0.5  /* Do not reduce trimfactor below this value until final optimization */
#define WEIGHTED_MEAN -0.01 /* TrimmedMean should use outlier weights when trimfactor <= WEIGHTED_MEAN (less robust if used early on, but helps convergence to a local maxima) */
#define KERNEL_MEAN 1 /* Use Kernel Density mean estimate, instead of Trimmed Mean */

#define MEAN_FIX 1 /* when a fragment is part of a larger alignment interval, downweight weight of sizing estimate by ratio of fragment size to alignment interval */
#define FORCE_RESIZE 0 /* force initial sizing estimate even if it lowers LP : the result is more robust but does not fully account for outlier scoring */
#define CENTERED_MEAN 0 /* If > 0 : Avoid using the K outermost alignment intervals at each end of the alignment for sizing, since they tend to be less reliable (K == CENTERED_MEAN) */
#define MULTIMODAL_TRIGGER 100.0 /* If MultiMode : Trigger slower multimodel Likelihood optimization when preliminary LP drops by this amount (or more)
				    NOTE : If multimode_force[], multimodel optimization is always applied, regardless of the value of the preliminary LP */

#define SPEEDFIX 1 /* experimental faster refinement settings */

#define MULTIMODE_FIX 1 /* WAS 0 */ /* Don't use sample values, just rely on MULTIMODE_MAXRANGE,MULTIMODE_MINRANGE (and largest,smallest sample value) and MULTIMODE_MP_MIN,MULTIMODE_MP_MAX */

#define MULTIMODE_DELTA ((RefineA || !OldMap) ? RefineStep2 : RefineStep1) /* minimum spacing (ratio difference) between sample values for global multimode scan : 
									      will increase speed when sample values are too closely spaced (high coverage) */
// WAS (SLOWMULTIMODE ? 0.05 : ((SPEEDFIX && (RefineA || !OldMap)) ? 0.20 : 0.10))

#define MULTIMODE_DELTA_MIN (SLOWMULTIMODE ? 0.10 : (RefineA || !OldMap) ? RefineStep6 : RefineStep5 /* 0.20 */) /* absolute value of minimum spacing in kb (see MULTIMODE_DELTA and MULTIMODE_MP), used for small intervals */

#define MULTIMODE_MINKB 1 /* WAS 0 */ /* Add MinKB to global scan values : this can avoid having to set a very small value of MULTIMODE_MINRANGE */

#define MULTIMODE_MINRANGE (SLOWMULTIMODE ? 0.001 : !SPEEDFIX ? (OldMap ? 0.5 : 0.01) : ((RefineA || !OldMap) ? RefineMinRange2 : RefineMinRange1)) 
                                                /* If MPROBEVAL : increase initial range of values sampled (if needed) to at least extend down to this ratio from the current best value 
							  at relative intervals of 1-MULTIMODE_MP_MIN (or absolute interval of MULTIMODE_DELTA_MIN, whichever is larger) */
#define MULTIMODE_MAXRANGE (SLOWMULTIMODE ? 5.0 : !SPEEDFIX ? (OldMap ? 2.0 : 3.0) : ((RefineA || !OldMap) ? RefineMaxRange2 : RefineMaxRange1)) 
                                                /* If MPROBEVAL : increase initial range of values sampled (if needed) to at least extend up to this ratio from the current best value 
							  at relative intervals of 1+MULTIMODE_MP_MAX */
#define MULTIMODE_MP RefineStepCnt /* WAS (SLOWMULTIMODE ? 14 : 3) */   /* If MPROBEVAL: During GoldenMean search use this many additional test points in each direction per interval (in addition to the main test point),
							    (at powers of 1-MULTIMODE_MP_MIN and 1+MULTIMODE_MP_MAX) during initial global jump (with no limit on number of points, but range limited by
							    MULTIMODE_MINRANGE,MULTIMODE_MAXRANGE) OR if total LP change of last iteration was above
							    MULTIMODE_MP_DELTA (AND neighboring interval was changed by at least the fraction MULTIMODE_MP_SCAN) OR 
							    if the test points fall within current search bounds Xlow[Q]..Xhigh[Q]. 
							    If MPBLOCK : During initial global jump if singlechange=1 it is only applied to intervals with LP improvements 
							    of at least 0.8 * MULTIMODE_MP_DELTA2 OR intervals next to the single changed interval OR intervals with peak LP at the end of the test point range. */
#define MULTIMODE_MP_MIN ((RefineA || !OldMap) ? RefineStep4 : RefineStep3)  /* see MULTIMODE_MP */
// WAS (SLOWMULTIMODE ? 0.05 : (RefineA || !OldMap) ? 0.20 : MULTIMODE_FIX ? 0.05 : 0.10)
#define MULTIMODE_MP_MAX ((RefineA || !OldMap) ? RefineStep4 : RefineStep3)  /* see MULTIMODE_MP */
// WAS MULTIMODE_MP_MAX (SLOWMULTIMODE ? 0.05 : (RefineA || !OldMap) ? 0.20 : MULTIMODE_FIX ? 0.05 : 0.10)

#define MULTIMODE_MP_SCAN RefineStep7 /* see MULTIMODE_MP */

// WAS MULTIMODE_MP_SCAN (SLOWMULTIMODE ? 0.05 : 0.10) /* see MULTIMODE_MP */

#define MPBLOCK 1 /* see MULTIMODE_MP */
#define MULTIMODE_MP_DELTA 50.0 /* see MULTIMODE_MP */
#define MULTIMODE_MP_DELTA1 1.0 /* If best single interval change causes LP to increase between MULTIMODE_MP_DELTA1 and MULTIMODE_MP_DELTA2 block changes to
							     the neighboring MULTIMODE_MP_CONFLICT intervals, unless those intervals have even better LP increases. This will speed up convergence */
#define MULTIMODE_MP_DELTA2 2000.0 /* WAS 60.0 */ /* while best single interval change causes LP to increase by this amount, force singlechange = 2 to reduce backtracking */
#define MULTIMODE_MP_CONFLICT 2 /* HERE HERE should be adaptive, WAS 1 */ /* see MULTIMODE_MP_DELTA1 */

#define MIN_DELTA1 0.01

#define KD_MIN 2.0 /* If MultiMode : skip checking sample values with Kernel Density below this value (ignored if MPROBEVAL) */


#define RESCAN_FIX 1 /* fix termination of "rescan" phase of site add/delete : ilast needs to be moved ahead to end of rescan range to ensure rescan range is fully checked one last time */

#define BLOCK_RESIZE 1.2 /* block increasing size of fragments based on KDmean unless current size is more than BLOCK_RESIZE * resKB during SizesEstimate 
					(-MultiMode will override and correctly handle small fragments near resolution limit where Likelihood is very sensitive to resolution model) */
#define BLOCK_RESOLVE 3 /* don't apply BLOCK_RESIZE if at least BLOCK_RESOLVE actual alignment pairs span the target interval with wt >= 0.99 : use the lesser of the original computed size or the mean of the alignment pairs */

static double SKIP_DIST=  0.0; /* Value set by -skip_dist : Don't check site additions more frequently than every SKIP_DIST kb */

//#define MERGE_WRAP (SLOWADDDELETE ? 1 : 0) /* wrap around ends of contigs when merging sites : 2x slower and rarely any better */
#define MERGE_FAST 2 /* faster merging by using previous bestLPA[] to avoid re-computing maps that are not overlapping sites being merged */

#define FAST_EXTEND 2 /* 1 : speed up 1st stage of extension refinement by not allowing site exchanges/migration (add one site, delete another nearby site) until 2nd stage with sizing adjustment
		         2 : Also speed up 2nd stage of extension refinement the same way */

#define UPDATEMAP 1 /* Update alignment location repeatedly (rather than just once at the end). Not done during initial site add/delete if extend >= 2 until last -endoutlierRef value */

#define ERRPLOT 0 /* print out x,y,wt triples from final consensus alignment (and trimmed coverage information from UpdateMap(), and more detailed Mprobeval logs) */

#define SVERB 0 /* more verbose LPdelta -Mprobeval output */

// static double ZERO_MINKB = 0.001; /* force minKB to this small value to allow interval values near zero to be checked : now in parameters.cpp */

int rverb = 0;
int tverb = 0;

#include "RGentigScore.h"

class Cydist {
public:
  double val;
  double wt;
};

extern double AddDelete(double LP, /* value of LP for current Y[] (may reflect an outdated map[],mapK[]) */
			double *TBmapWT, /* If != 0 : TBmapWT[m] is a the -TB based weight for map m to be applied when summing up the total log(LR+LRbias) value */
			double *bestLPA,/* current best per map LP values */
			double *newLPA,/* per map LP values (return value) */
			double *mapWT,
			int n, double *Hcuts, int *Hdel,/* complete consensus map */
			int &N, double *Y,/* current consensus map */
			int MD, /* number of maps */
			int *MX, double **X, /* X[m=0..MD-1][j=0..MX[m]+1] : map sites */
			int lc,/* If left end of Y is a linear chromosome end */
			int rc,/* If right end of Y is a linear chromosome end */
			int **map, /* map[m=0..MD-1][j=0..MX[m]+1] is old index in Hcuts[0..n+1] of X[m][j], 
					       map[MD][j=0..n+1] is index in Y[0..N+1] of Hcuts[j] */
			int **mapK, /* mapK[m=0..MD-1][j=0..MX[m]+1] is old K index offset in Hcuts[0..n+1] of X[m][j] */ 
			Csetlimit *limit,
			int **nmap, /* nmap[m=0..MD-1][j=0..MX[m]+1] is new index in Y[0..N+1] of X[m][j], 
						nmap[MD][j=0..N+1] is index in Hcuts[0..n+1] of Y[j] */
			int **nmapK, /* nmapK[m=0..M-1][j=0..MX[m]+1] will be new K index offset in Y[0..N+1] of X[m][j]
						 However nmapK[0] == Lij, mapK[MX[m]+1] == Rij (-ve if ends are outliers) */
			Ccontig *pcontig, /* pointer to complete contig information for debugging */
			double &trimfactor, /* The trimmed mean parameters : smallest value 0 corresponds to normal mean
					      largest value 0.5 corresponds to median */
			int resize, /* If != 0 :  resize map instervals using SizesEstimate() */
			int &MultiMode, /* Only used if resize != 0 */
			int *skip, /* skip[i=1..n]: > 0 if site Hcuts[i] need not be checked for addition/deletion */
			int &Lfrozen, int &Rfrozen, /* refine() args */
			int merge /* If != 0 : try to merge labels closer than maxKB, but only if that improves LP */
			);

extern double mprobevalwin(double LP, /* value of LP for current Y[] (may reflect an outdated map[],mapK[]) */
			   double *bestLPA,/* current best per map LP values */
			   double *newLPA,/* per map LP values (return value) */
			   double *mapWT,
			   int n, double *Hcuts, int *Hdel,/* complete consensus map */
			   int &N, double *Y,/* current consensus map */
			   int MD, /* number of maps */
			   int *MX, double **X, /* X[m=0..MD-1][j=0..MX[m]+1] : map sites */
			   int lc,/* If left end of Y is a linear chromosome end */
			   int rc,/* If right end of Y is a linear chromosome end */
			   int **map, /* map[m=0..MD-1][j=0..MX[m]+1] is old index in Hcuts[0..n+1] of X[m][j], 
						  map[MD][j=0..n+1] is index in Y[0..N+1] of Hcuts[j] */
			   int **mapK, /* mapK[m=0..MD-1][j=0..MX[m]+1] is old K index offset in Hcuts[0..n+1] of X[m][j] */ 
			   Csetlimit *limit,
			   int **nmap, /* nmap[m=0..MD-1][j=0..MX[m]+1] is new index in Y[0..N+1] of X[m][j], 
						   nmap[MD][j=0..N+1] is index in Hcuts[0..n+1] of Y[j] */
			   int **nmapK, /* nmapK[m=0..M-1][j=0..MX[m]+1] will be new K index offset in Y[0..N+1] of X[m][j]
						    However nmapK[0] == Lij, mapK[MX[m]+1] == Rij (-ve if ends are outliers) */
			   double *TBmapWT, /* If != 0 : TBmapWT[m] is a the -TB based weight for map m to be applied when summing up the total log(LR+LRbias) value(s) */
			   int imin, int imax, /* Only interested in changes in the region Hcuts[imin] .. Hcuts[imax] */
			   Ccontig *pcontig, /* pointer to complete contig information for debugging */
			   double &trimfactor, /* The trimmed mean parameters : smallest value 0 corresponds to normal mean
						  largest value 0.5 corresponds to median */
			   int resize, /* If != 0 :  resize map instervals using SizesEstimate() */
			   int &MultiMode, /* Only used if resize != 0 */
			   int *skip, /* skip[i=1..n]: > 0 if site Hcuts[i] need not be checked for addition/deletion */
			   int Lfrozen, int Rfrozen, /* refine() args */
			   int &AddDeleteCnt /* cumulative count of number of sites added/deleted */
			   );


extern double mprobevalwinResize(int n, double *Hcuts, /* Hcuts[i0..n+1] : complete consensus map */
				 double qLPstart, /* value of LP for current Y[] */
				 double *bestLPA,/* current best per map LP values */
				 double *newLPA,/* per map LP values (return value) */
				 double *mapWT,
				 int &N, double *Y,/* current consensus map */
				 int MD, /* number of maps */
				 int *MX, double **X, /* X[m=0..MD-1][j=0..MX[m]+1] : map sites */
				 int lc,/* If left end of Y is a linear chromosome end */
				 int rc,/* If right end of Y is a linear chromosome end */
				 int **map, /* map[m=0..MD-1][j=0..MX[m]+1] is old index in Hcuts[0..n+1] of X[m][j], 
							map[MD][j=0..n+1] is index in Y[0..N+1] of Hcuts[j] */
				 int **mapK, /* mapK[m=0..MD-1][j=0..MX[m]+1] is old K index offset in Hcuts[0..n+1] of X[m][j] */ 
				 Csetlimit *limit,
				 int **nmap, /* nmap[m=0..MD-1][j=0..MX[m]+1] is new index in Y[0..N+1] of X[m][j], 
							 nmap[MD][j=0..N+1] is index in Hcuts[0..n+1] of Y[j] */
				 int **nmapK, /* nmapK[m=0..M-1][j=0..MX[m]+1] will be new K index offset in Y[0..N+1] of X[m][j]
							  However nmapK[0] == Lij, mapK[MX[m]+1] == Rij (-ve if ends are outliers) */
				 double *TBmapWT, /* If != 0 : TBmapWT[m] is a the -TB based weight for map m to be applied when summing up the total log(LR+LRbias) value(s) */
				 int *Ycnt,
				 Cydist **Ydist,
				 double *newY,/* newY[0..N+1] is the quick estimate based on KDmean */
				 int Imin, int Imax, /* Only interested in changes in to intervals Y[Imin..Imax] */
				 int *MMforce,
				 double *fixY,
				 double *fixYdel,
				 Ccontig *pcontig, /* pointer to complete contig information for debugging */
				 int nochange,
				 int LfrozenI, int RfrozenI,
				 double &trimfactor, /* The trimmed mean parameters : smallest value 0 corresponds to normal mean
							largest value 0.5 corresponds to median */
				 int &Lfrozeni, int &Rfrozeni
				 );


extern double SizesEstimate(double LP, /* value of LP for current Y[] (may reflect an outdated map[],mapK[]) */
			    double *TBmapWT, /* If != 0 : TBmapWT[m] is a the -TB based weight for map m to be applied when summing up the total log(LR+LRbias) value(s) */
			    double *newLPA,/* per map LP values (return value) */
			    double *mapWT,
			    int n, double *Hcuts, int *Hdel,/* complete consensus map */
			    int N, double *Y,/* current consensus map */
			    int MD, /* number of maps */
			    int *MX, double **X, /* X[m=0..MD-1][j=0..MX[m]+1] : map sites */
			    int lc,/* If left end of Y is a linear chromosome end */
			    int rc,/* If right end of Y is a linear chromosome end */
			    int **map, /* map[m=0..MD-1][j=0..MX[m]+1] is old index in Hcuts[0..n+1] of X[m][j], 
						 map[MD][j=0..n+1] is index in Y[0..N+1] of Hcuts[j] */
			    int **mapK, /* mapK[m=0..MD-1][j=0..MX[m]+1] is old K index offset in Hcuts[0..n+1] of X[m][j] */ 
			    Csetlimit *limit,
			    int **nmap, /* nmap[m=0..MD-1][j=0..MX[m]+1] is new index in Y[0..N+1] of X[m][j], 
						    nmap[MD][j=0..N+1] is index in Hcuts[0..n+1] of Y[j] */
			    int **nmapK, /* nmapK[m=0..M-1][j=0..MX[m]+1] will be new K index offset in Y[0..N+1] of X[m][j]
						     However nmapK[0] == Lij, mapK[MX[m]+1] == Rij (-ve if ends are outliers) */
			    int imin, int imax, /* Only interested in changes in the region Hcuts[imin] .. Hcuts[imax] */
			    Ccontig *pcontig, /* pointer to complete contig information for debugging */
			    double trimfactor, /* The trimmed mean parameters : smallest value 0 corresponds to normal mean
						  largest value 0.5 corresponds to median */
			    int oldmap, /* If 1 use old mapping (map,mapK) instead of new mapping (nmap,nmapK) to compute draft consensus */
			    int nochange,/* If 1 don't make any changes to Y[], Hcuts[], but DO return the new LP value */
			    int &failed, /* Set to 1 IF Y[] was not updated, 0 otherwise : 
					     NOTE : Even if Y[] was not updated, return value may be better than LP, if LP did not
					     reflect the latest update to map[],mapK[] */
			    int &mapupdated, /* Set to 1 IFF UpdateMap() was called successfully (at least once) */
			    int &Lfrozen, int &Rfrozen
			  );

/* actual Haplotype refinement function  : input is a Haplotyped map */
extern double HaploTypeR(int n, double *Hcuts, int *&HapSite, double *&HapDelta,/* complete Haplotype map pair specification : implies Hcuts1, Hcuts2, Y1, Y2, N1, N2, see hsetmap() */
			 int lc,/* If left end of Y is a linear chromosome end */
			 int rc,/* If right end of Y is a linear chromosome end */
			 int MD, /* number of maps */
			 int *MX, double **X, /* X[m=0..MD-1][j=0..MX[m]+1] : map sites */
			 int **map1, /* map1[m=0..MD-1][j=0..MX[m]+1] is original index in Hcuts1[0..n+1] of X[m][j], 
					map1[MD][j=0..n+1] is index in Y1[0..N1+1] of Hcuts1[j] (not initialized) */
			 int **map2, /* map2[m=0..MD-1][j=0..MX[m]+1] is original index in Hcuts2[0..n+1] of X[m][j], 
					map2[MD][j=0..n+1] is index in Y2[0..N2+1] of Hcuts2[j] (not initialized) */
			 int **mapK1, /* mapK1[m=0..MD-1][j=0..MX[m]+1] is original K index offset in Hcuts1[0..n+1] of X[m][j] */ 
			 int **mapK2, /* mapK2[m=0..MD-1][j=0..MX[m]+1] is original K index offset in Hcuts2[0..n+1] of X[m][j] */ 
			 Csetlimit *limit1,
			 Csetlimit *limit2,
			 double *TBmapWT, /* If != 0 : TBmapWT[m] is a the -TB based weight for map m to be applied when summing up the total log(LR+LRbias) value */
			 Ccontig *pcontig, /* pointer to contig with input map ids in contig->contig[0..MD].mapid, and initial pcontig->HapSite,pcontig->HapDelta (aliased as HapSite,HapDelta args)
					      Final refined Haplotype Information will be in pcontig->HapSite, pcontig->HapDelta, pcontig->HapSiteScore etc
					      and final map alignments to Hcuts[] in pcontig->sitecnt,pcontig->sitecntL (based on best scoring Allele for each map) */
			 int *skip, /* skip[1..n] : skip[i] > 0 if site Hcuts[i],HapSite[i] should not be changed (typically these sites have HapSite[i]== 0, or are in a frozen region) */
			 int& Lfrozen, int& Rfrozen /* refine() args : may change when Hcuts[] changes in repositionH() */
			 );

/* wrapper function to real function HaplotypeR() : input is a non-Haploytpe map */
extern double HaploType(int n, double *Hcuts, int *Hdel,/* complete consensus map : implies Y[], N via setmap() */
			int lc,/* If left end of Y is a linear chromosome end */
			int rc,/* If right end of Y is a linear chromosome end */
			int MD, /* number of maps */
			int *MX, double **X, /* X[m=0..MD-1][j=0..MX[m]+1] : map sites */
			int **map, /* map[m=0..MD-1][j=0..MX[m]+1] is original index in Hcuts[0..n+1] of X[m][j], 
				      map[MD][j=0..n+1] is index in Y[0..N+1] of Hcuts[j] (not initialized) */
			int **mapK, /* mapK[m=0..MD-1][j=0..MX[m]+1] is old K index offset in Hcuts[0..n+1] of X[m][j] */ 
			double *TBmapWT, /* If != 0 : TBmapWT[m] is a the -TB based weight for map m to be applied when summing up the total log(LR+LRbias) value */
			Ccontig *pcontig, /* pointer to complete contig information : Haplotypes are returned in pcontig->Haplotype[] */
			int *skip,
			int &Lfrozen, int &Rfrozen /* refine() args */
			);

#endif
