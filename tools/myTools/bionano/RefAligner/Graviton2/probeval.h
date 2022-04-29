#ifndef PROBEVAL_H
#define PROBEVAL_H

static Ident probeval_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/probeval.h 11431 2020-08-07 19:52:29Z tanantharaman $");

#define MPROBEVAL_MAPORDER 30 /* In mprobeval() : order maps by location in contig (from left to right, based on limit->ileft,iright) to improve memory reference locality for newLPdelta[] and newLPa[] arrays
				 Modifies MIN_MEM to not shring Fmem[],Imem[],LRdeltaD[] memory since this memory needed will fluctuate randomly. 
				 Instead call madvice(MADV_DONTNEED) on Fmem[],Imem[],LRdeltaD[] memory at least once every MP_MIN_TIME seconds */
#define MPROBEVAL_SWAPOUT 1 /* swap out CdoubleBlock memory blocks in mprobeval() (beyond 1 block per thread) to reduce memory footprint of mprobeval() for large contigs (when called from HaploTypeR with Allele > 0) */

#define MPROBEVAL_LAZYRESTORE 2 /* WAS USE_MIC ? 2 : 0) */ /* Restore swapped out CdoubleBlock memory from mprobeval() calls in HaploTypeR, only when first needed (currently only when HINDELITER(iter))
							      If >= 2 : Also swap back out memory blocks from left end of contigs whenever percentage of swapped out memory blocks drops below 
							                LAZYRESTORE_MIN until the percentage swapped out memory blocks rises back to LAZYRESTORE_MAX
								        Only do this if VmRSS exceeds VMRSS_MIN percent of -maxmem (to avoid wasting time if there is plenty of memory) */
#define LAZYRESTORE_MIN 50
#define LAZYRESTORE_MAX 60
#define VMRSS_MIN 60

#define SETLIMIT_CACHE 1
#define SETLIMIT_CACHE_FORCE 0 /* If SETLIMIT_CACHE detects an overlap with [Ymin..Ymax], force alignment to be computed, even if setlimit() returns a changed alignment that no longer overlaps [Ymin..Ymax] */

#define SKIP_DIST_FIX 1 /* WAS 0 */ /* fix skip_sites() handling of duplicate sites */

#define VECTOR_AVX 1 // Enable AVX vector code (if USE_AVX && !USE_AVX512 && USE_SSE && USE_MFLOAT)
#define VECTOR_MIC 1 // Enable MIC/AVX512 vector code (if (USE_MIC || USE_AVX512) && USE_MFLOAT )
#define VECTOR_DEBUG 0 // Debug intrinsic vector code (1 : high-level, 2 : low-level, slower)

#define VECTOR2_ADD 1 // seperate out critical AL() and AR() calls suitable for vectorization (provides speedup on its own, but may use larger jmin..jmax ranges and generate optimistic LR predictions)

#define VECTOR2_DEBUG 0 // Debug intrinsic vector code for AL() and AR() calls (1 : high-level, 2 : low-level slower)

#define VEC_SIZ ((USE_MIC || USE_AVX512) ? 16 : 8)
#define VPAD ((USE_MIC || USE_AVX512) ? 15 : USE_AVX ? 7 : 0) /* padding of arrays by up to this many elements to allow vectorization */
#define PADDING ((USE_MIC || USE_AVX512) ? 64 : USE_AVX ? 32 : 0) /* Align and Pad Allocated memory blocks to nearest multiple of this many bytes : 
								     should correspond to size of mm512 vector with MIC or AVX512 and mm256 with AVX */		   
#define VEC_MSIZ ((USE_MIC || USE_AVX512) ? 3 : 3) /* minimum vector size for mprobeval_delta() */// HERE HERE : optimize MIC & AVX512
#define VEC_MSIZ2 ((USE_MIC || USE_AVX512) ? 3 : 2) /* minimum vector size for mprobeval_add() */// HERE HERE : optimize MIC & AVX512

#define VEC_MSIZQ ((USE_MIC || USE_AVX512) ? 1 : 1) /* minimum vector size for qprobeval3D_GTH() */// HERE HERE : optimize

#define M512_ABS 1 /* Use mm512_abs_ps() instruction */

#define FLIK_FIX 1 /* Avoid underflow in float version of FLIK() by passing FLJ value, which is often >> 1. Also applies to FR*() */
#define URANGE_FIX 1 /* Don't assume that largest FLIK() value is obtained with U value at Lij. Also applies to FR*() for Rij 
			1 : Fixes that don't affect run time signficantly
			2 : Fixes that degrade run time significantly (to be fixed later)
		      */
#define ERFC_FIX 1 /* Use erfc instead of erf when args are expected to be > 0 */
#define URANGE_QFIX 1 /* make URANGE computation in qprobeval better match mprobeval */

#define HAPSTOP_DEBUG 0 /* more verbose output to debug -HapLocalStop */
#define HAPSTOP_FIX 1 /* skip some loops with zero iterations when -HapLocalStop is used */

#define LRCUM_FIX 0 /* HERE TRY 1 */ /* new code to avoid LRcum underflow */

#define MEM_MMAP 0 /* try to MMAP, possibly using huge pages */

#define ALIGN 16 // alignment of (I,J) based 2-D arrays in elements (int or double) : size in bytes should be a multiple of cache line size (must be power of 2)
#define ALIGN2 (PAGE >> 2) /* WAS (MIC ? 16 : 1024) */ // alignment of (K,I,J) based 3-D arrays in elements (int or MFLOAT) : size in bytes should be a multiple of the page OR cache line size (must be power of 2)

#define DELTA_LIM 0 /* TRY 16 */ /* see DELTA_MIN. Also used to increase DELTA_X during extension refinement (until sizing is adjusted) */ 
#define DELTA_MIN 0 /* If != 0 : minimum sites lookahead on either X or Y that must be matched in distance by lookahead on the other map (up to DELTA_LIM). If so the lookahead
		       on the original map is reduced to DELTA_MIN. NOTE if DELTA_MIN > 0, mprobeval sees inconsistent ranges in opposite scanning directions, better to increase
		       DELTA_X during extension refinement to DELTA_LIM */
#define DELTA_STOP HapLocalStop /* stop optimizing delta values if there have been no improvements in current or DELTA_SCAN neighboring bins */
#define END_FIX 2 /* FL and FU terms are extended into the I-K..I region up to Yc(I,K) (if K > 0) : Also EXTEND is applied when END_VITERBI==0*/
#define END_VITERBI 0 /* Use Viterbi approximation for ends of alignments (as in refalign()), rather than more accurate UL,UR computation described in Refinement Document */
#define EPS RefineEps1 /* targetted likelihood accuracy : terms that do not contribute at least this fraction to the UL,UR likelihood terms may be ignored in mprobeval */
#define EPS2 RefineEps2 /* targetted likelihood accuracy : AL,AR terms that do not contribute at least this fraction to total likelihood may be ignored in mprobeval */
#define EPS3 RefineEps3 /* targetted likelihood accuracy : terms that do not contribute at least this fraction to AR terms may be ignored in mprobeval */
#define EXTEND 500 /* If > 0 : Allow X map to extend beyond consensus Y by EXTEND_MAX + DELTA_X sites : treat Y as being extended to match length but without additional sites.
		     (this causes a penalty for any extension due to misaligned sites (limited by PoutlierEnd), unlike in RefAligner)
		     EXTEND must be a constant = the largest value of EXTEND_MAX (names should be reversed!)
		     Would work better with map LR outlier detection (see LRbias) so extension is not driven by outlier maps */
#define EXTEND_MAX (extend ? Extend : 0) /* limit on alignment extension. Will cause large sizing penalty (limited by PoutlierEnd) for maps that extend end of the contig beyond this limit.
						     Even for maps that extend less than this limit, sites in the extension region are treated as misaligned sites (also limited by PoutlierEnd)]. */

#define FIX_U deltaEndRef /* Apply FR() or FL() for labels within deltaExtXRef+FIX_U-1,deltaExtYRef+FIX_U-1 of end instead of DELTA_X,DELTA_Y (deltaXRef,deltaYRef) : slightly slower but more accurate */
extern int DELTA_XU, DELTA_YU;// initialized at start of qprobeval() and mprobeval()


#define EXTRAPOLATE 0 /* Use more agressive extrapolation of missing Y alignments from all other alignments, rather than just the nearest alignment in either direction */
#define FIX_RANGE_Y 1 /* when RANGE_Y==0, avoid using RANGE site(s) as max alignment error (use RANGE*Xtheta distance as max alignment error instead), 
			 so addition/deletion of sites does not change the alignment range (as far as possible) and hence makes mprobeval more accurate */

#define LOCAL_RANGE LocalRange /* WAS 0 */ /* Use RANGE * (X[M]-X[1])/(M-1) instead of RANGE * Xtheta in setlimit() : modifies RANGE_Y=0 etc */

#define GLOBAL global /* 0 : Constrain alignments based on previous best known alignment and limit changes to RANGE labels (or (RANGE+1)*Xtheta kb)
			 1 : Use global alignment : much slower, generally only needed with repeats (or initial alignment) */
#define HDEBUG 0
#define HEAP_MINIMIZE 1 /* WAS 0 */ /* minimize real memory (at expense of more calls to malloc()/free() in non-multithreaded section) */
#define HFIX 3 /* 1 : fix mprobeval bug triggered by out of order delta values during Haplotyping (and bug in extrapolation of setlimit())
		  2 : Call UpdateMap() in HaploType()
		  3 : Call repositionH() in HaploType() */
#define LRANGE Lrange  /* see LRANGE_FIX */

#define LTRACE 0
#define LVERB 0 /* verbose for debugging (currently display coverage profile) */
#define MAP_RES 0.010 /* treat consensus sites in Hcuts[] as indistinguishable if they are within this distance : they will be replaced by a single site */


/* NOTE : MDEBUG requires use of gcc compiler that supports __float128 */
#define MDEBUG 0 /* If MPROBEVAL : check LP values return by mprobeval() using qprobeval (then exit if LP error exceeds MDEBUG_ERR in refine.cpp) */
// extern int MDEBUG;// Must initialize to 1 for arrays needed to be allocated, then set to 0 before calling mprobeval that are not being checked

#define MDEBUG_F -1 /* display LR terms factored around interval MDEBUG_F..MDEBUG_F+1 AFTER deleting site MDEBUG_S when map m==MDEBUG_M && N==MDEBUG_N-1. Typically MDEBUG_F == MDEBUG_S - 1 */
#define MDEBUG_I -1 /* For site deletion use ((N==MDEBUG_N || I < MDEBUG_S) ? I : I-1), where I is the intended value of MDEBUG_I before site deletion */
#define MDEBUG_K -1
#define MDEBUG_J -1
#define MDEBUG_M -1 /* >= 0 : display LR terms >= MDEBUG_MIN when map m==MDEBUG_M AND (see below) */
//extern int MDEBUG_M;
#define MDEBUG_MIN 1e-305
/* more detailed tracing in mprobeval of A[MDEBUG_I][MDEBUG_K][MDEBUG_J].AL,UL,AR,UR IFF MDEBUG_M >= 0 && MDEBUG_I >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 */
#define MDEBUG_N  -1 /* N value for mprobeval() being tested */
#define MDEBUG_N1 -1 /* (MDEBUG_N-1)*/
#define MDEBUG_N2 -1 /* (MDEBUG_N) */
#define MDEBUG_N3 -1 /* MDEBUG_N */ /* Also trace in qprobeval() when MDEBUG_N3 == N */

#define MDEBUG_S  -1 /* display LR terms factored around deletion of site MDEBUG_S when map m==MDEBUG_M && N==MDEBUG_N */

#define MDEBUG_SETLIMIT ((MDEBUG && MDEBUG_M >= 0/* && rverb */) ? 2 : 0) /* trace calls to setlimit(), when MDEBUG_M>=0 && MDEBUG_M==m (1 = basic trace, 2 = verbose trace) */
#define MDEBUG_ST -1 /* trace LRsite[MDEBUG_ST] (all terms >= MDEBUG_MIN) : Useful to track down LRSITE_ERR assertion failures */
// extern int MDEBUG_ST;

#define MDEBUG_T -1 /* display LR terms factored for addition at site addcnt[MDEBUG_T] (in interval Y[MDEBUG_TS..MDEBUG_TS+1]) when map m==MDEBUG_M and N==MDEBUG_N */
#define MDEBUG_TS -1 /* display LR terms factored around site MDEBUG_TS+1 when map m==MDEBUG_M and N==MDEBUG_N+1 */
#define MDEBUG_TT -1 /* LRaddT[MDEBUG_TS][MDEBUG_TT] term to be tracked along with sum(LRadd[0..MDEBUG_TS]) when m==MDEBUG_m && N==MDEBUG_N && MDEBUG_T >= 0 */
                                       /* If MDEBUG_T < 0, then instead track LRdeltaD[MDEBUG_TS][MDEBUG_TT] */
#define MTRACE 0  /* trace AL*UR sum if MDEBUG_M == m and MDEBUG_N == N (also AR*UL sum if DEBUG>=2) in qprobeval and mprobeval */

#define MFLOAT1 MFLOAT // used in qprobeval for FP arrays
#define MFLOAT2 MFLOAT // used in mprobeval for FP arrays

#define OUTLIER(x,y) (maptype || OUTLIER_DELTA(x-y)) /* only test alignments for outliers if they satisfy this condition : x = query, y = reference, maptype = 0(bnx) or 1 (Cmap) */
#define OUTLIER_FIX 1 /* fix Lij & Rij values when there are end outliers : instead of -ve values set to correct overlap end, the boundary of (I-K..I) (only in qprobeval) */
#define OUTLIER_MAXBC 0.0 /* 0.0 : Enable backward compatible scoring when outlierMax >= 1000.0
			     > 0.0 : assume outlierMax is always enabled and reduced to min(OUTLIER_MAXBC,outlierMax) : recommend 100.0 (slightly faster when outlierMax is always <= OUTLIER_MAXBC) */
#define OUTLIER_TYPE OutlierType 
/* OUTLIER type :
   0 : Exclude both misaligned sites and sizing error in outlier (requires large DELTA_X=DELTA_Y >= n to evaluate exactly, unless outlierLambda < 100)
   1 : Exclude only sizing error in outlier (still has accuracy problems due to allowing very large (in kb) outliers with same penalty as smaller size outliers, 
       Possible fixes : Use -outlierMax OR -outlierLambda (-outlierBC NOT implemented during refinement) */
#define PM_SPLIT 1 /* keep PM term seperate from FA term (to better match Refinement Document) */
#define SCALE_MONOTONIC 0 /* LogScale[] factors are required to be monotonic (requires less renormalizations, but may underflow with USE_MFLOAT==1) */
#define SITE_PEN Site_Pen /* log-likelihood prior penalty per site in consensus map :  Assuming at least one map supports each site, SITE_PEN should be at least:
			     log(sqrt(2*PI*(SF[0]*SF[0]+SD[0]*SD[0]*fragsiz))/F)/K = 2.9 for fragsiz = 50kb, F = 1.0/100kb, sf=0.2, sd=0.2, K=2.0. ???
			     Note that when using large LRbias values such that most molecules score below LRbias, SITE_PEN can result in most sites being deleted spuriously */
extern int SKIPMAP_ADD; /* Extend skip map heuristic in mprobeval.cpp to include addloc[0..addcnt-1]
			   Only applied for even iteration mprobeval calls from HaploTypeR() (with Ymin < 0.0) :
			   1 : Just skip computing addcnt > 0 code for maps that satisfy the condition.
			   2 : Skip all computation for maps that satisfy the condition
			       Once this is turned on (with 2) it is not possible to compute HapSiteScore[] and SiteScore[] for all I=1..N, since mprobeval newLPd[m][I] values will only
			       be reliable in regions that were still marked active (addcnt[t] are present in the interval Y[I-1] .. Y[I+1]). This requires the fastHapSite() etc code to
			       be updated to skip those regions when checking for HapSite[] changes as well as when updating HapSiteScore[i],SiteScore[i].
			*/

#define STRACE -1
#define VMEM_MINIMIZE (USE_MIC ? 128 : 256) /* minimize virtual memory (at expense of more calls to mmap()/madvise()/munmap() in multithreaded section) 
					       To reduce calls to mmap()/madvise()/munmap() allocates blocks of VMEM_MINIMIZE x 1024 doubles at a time (per thread) */

#define MAX_INDEL 4 /* maximum indels per interval on a single allele (if there are more, we may not be able to reverse them with HINDEL_PHASED) */


#if USE_MFLOAT
#define OUTLIER_DELTA(err) (!OUTLIER_MAX || fabsf(err) < outlierMaxBC)
#define OUTLIER_OLD(outlierMaxBC)  (!OUTLIER_MAX || (!OUTLIER_MAXBC && outlierMaxBC >= 1000.0f))    /* test if we need to support outlierMax >= 1000.0 the old way */
#else
#define OUTLIER_DELTA(err) (!OUTLIER_MAX || fabs(err) < outlierMaxBC) 
#define OUTLIER_OLD(outlierMaxBC)  (!OUTLIER_MAX || (!OUTLIER_MAXBC && outlierMaxBC >= 1000.0))    /* test if we need to support outlierMax >= 1000.0 the old way */
#endif

#if USE_EPOW
  #if USE_VEXP && USE_MFLOAT
    #include "vec_exp.h"
    #define expF vec_expf
//   #define erfcf vec_erfcf
  #elif USE_SEXP && USE_MFLOAT
    #include "simd_math.h"
    #define expF expapprox
//  #define logF logapprox
  #else
    #include "epow.h"
    #define expF epowF
//    #define exp epow
  #endif
#else
  #if USE_MFLOAT >= 1
    #define expF expf
  #else
    #define expF exp
  #endif
#endif

#if USE_MFLOAT
#define erfF erff
#define erfcF erfcf
#define sqrtF sqrtf
#define fabsF fabsf
//#define logF logf // HERE HERE
#define logF log 
#else
#define erfF erf
#define erfcF erfc
#define sqrtF sqrt
#define fabsF fabs
#define logF log
#endif


class Csetlimit {
 public:
  int ileft, iright;// Hcuts[ileft .. iright] correspond to Y[IMIN-LRANGE .. IMAX+LRANGE], the range where changes in Y[] will affect the current molecule's predicted likelihood due to changes from D[I] or addloc[] in mprobeval
  int Ileft,Iright;/* Y[Ileft .. Iright] correspond to Hcuts[ileft .. iright] : these values may not be valid since the labels may have been deleted. 
		      If so the nearest surrounding labels (from Hcuts[ileft..iright]) can be used (or set ileft = -1 and let setlimit() recompute the correct values) */
  double Xleft, Xright;// Yleft = Hcuts[ileft] - Xleft and Yright = Hcuts[iright] + Xright : Yleft..Yright are the range where changes in Y[] can affect the current molecule's likelihood in qprobeval
  // NOTE : all values are marked invalid with ileft = -1
};

int setlimit(int m,
	     int n, double *Hcuts,
	     int N, double *Y,
	     int M, double *X,
	     int *map,
	     int *mapK,
	     Csetlimit *limit,
	     int *mapMD,
	     int *nmapMD, /* mapping from Y[I] to Hcuts[i] */
	     double Ymin, double Ymax,/* If >= 0.0 : If X does not overlap this range +- RANGE*Xlambda, return 0 (without computing results) */
	     int &IMIN,
	     int &IMAX,
	     int *Jmin,
	     int *Jmax,
	     int *Imin,
	     int *Imax,
	     double * &JXmin,
	     double * &JXmax,
	     int * &Jmid,
	     int &IMIN1,
	     int &IMAX1,
	     Ccontig *pcontig, /* pointer to complete contig information for debugging */
	     int force 
	     );

#define SCALE_MONOTONIC 0 /* LogScale[] factors are required to be monotonic (requires less renormalizations, but may underflow with USE_MFLOAT==1) */

#undef DBL_MAX
#define DBL_MAX 4.494232e+307 // largest value that can be inverted without underflow!

#define MaxFloat (USE_MFLOAT ? FLT_MAX : DBL_MAX) // largest floating point number used
#define MinFloat (USE_MFLOAT ? FLT_MIN : DBL_MIN) // smallest positive floating point number (NOT denormalized)
#define MaxLR (USE_MFLOAT ? (SCALE_MONOTONIC ? 1e+18f : /* WAS 1e+15f */ 1e+19f) : 1e+100) // Largest value of A.LR(I,*,*) AFTER LogScale[I] normalization (largest value BEFORE normalization can be MaxFloat)
#define MinLR (USE_MFLOAT ? 1e-28f : 1e-120) // Largest ratio drop in A.LR(I) (OR A.AR(I)) along best alignment from HWM to right (OR left) end (while LogScale[] remains at HWM value).
#define MaxFA (USE_MFLOAT ? /* WAS 3e+20f */ 3e+15f : 1e+200) // largest value of FA : FA is reduced to this value in FA (or FAJ), which lowers scores of very large internal intervals
#define MaxUR (USE_MFLOAT ? /* WAS 1e+20f */ 1e+15f : 1e+100) // Largest value of A.UR  : UR and UL are reduced to this value in FR and FL (or FLJ), which lowers scores of very large end intervals
#define MinUR (USE_MFLOAT ? 1e-28f : 1e-280) // Smallest value of A.UR (effectively smallest value of PoutlierEnd)
#define LRCUM_MAX (USE_MFLOAT ? 1.0 : (MDEBUG ? 1e+200: 1.0)) // Rescale LRcum whenever it exceeds this value (OR to satisfy MaxScale)
#define MaxScale (USE_MFLOAT ? /* WAS (SCALE_MONOTONIC ? 3.0f : 3e+8f) */ 1e+100 : 1.0e+100) // MaxScale is largest value of Scale(I) == exp(LogScale[I] - LS(I)) 

/* NOTE : To guarantee no overlows or underflows in qprobeval() or mprobeval() :
   1. During A.LR 3D recurrance LR(I,K,J) = PM(I,K) * Sum(G) exp(LogScale[G]-LogScale[I])) * Sum(T,H) LR(G,T,H) * FA(IKJ,GTH)
         Since PM <= 1 and exp(LogScale[G]-LogScale[I]) <= 1 we need :    MaxLR * MaxFA < MaxFloat
	 NOTE : If SCALE_MONOTONIC==0, then exp(LogScale[G]-LogScale[I]) may be > 1, but in that case early re-normalization of LogScale[I] is used to recover from overflow
   2. During LRcum accumulation LRcum(I) = LRcum(I-1) + Sum(K,J) LR(I,K,J) * UR(I,K,J) * exp(LogScale[I] - LS(I)), hence we need :    MaxLR * MaxUR * MaxScale < DBL_MAX
      2b. Lowering MaxScale, will implicitly lower LRCUM_MAX (by raising LS) down to as low as MaxLR*MaxScale/MaxFloat. Hence to avoid underflow of LRcum we need :    MaxLR * MaxScale > 1.0
      2c. To avoid underflow of AR*UL or AL*UR we need :      (MaxLR * MinLR) * MinUR > MinFloat  (NOTE : Since LRcum uses doubles, may no longer apply)
              (assuming AR,AL >=  MaxLR * MinLR after normalization at best alignment ends, given that AR,AL should be MaxLR at HWM point if SCALE_MONOTONIC==1
	      If SCALE_MONOTONIC==0, then underflow is much less likely, so MaxLR can be reduced)
      2d. To avoid overflow of AR*UL or AL*UR we need: MaxLR * MaxUR < MaxFloat (NOTE : Most of the LRcum code uses doubles, so this may not apply anymore)
   3. MinLR << 1.0 is required when refining extensions with very stringent -endoutlier. Hence MinLR should be comparable to MinUR and as small as possible (does not apply if SCALE_MONOTONIC==0)
   4. With SCALE_MONOTONIC==0 to avoid NaN values in LR(G,T,H) * FA(IKJ,GTH) * exp(LogScale[G]-LogScale[I]) because intermediate values overflow, even if final product is 0, 
                  exp(LogScale[G]-LogScale[I]) cannot exceeed MaxFloat/max(MaxFA,MaxLR) (NOTE : Currently only applies if scaling in qprobeval is in inner loop, otherwise it
		  can be as large as MaxFloat (or even DBL_MAX if Scaling factor is not cast to float in intrinsic code)
        MaxLR could be reduced (down to 1.0) compared to SCALE_MONOTONIC==1, though this risks underflow for alignments that are (so far) sub-optimal, but may turn out to become the best alignment later.
	A value of 1e15 may be a better compromise with handling multiple good alignments at different locations that require smaller values of MaxLR.
   5. mprobeval requires AL(I,K,J) * AR(I,K,J) to not overflow for any I,K,J : this requires MaxLR * MaxLR < MaxFloat
*/

extern int rverb;

extern int tverb, stage;// for debugging

namespace probeval {

extern MFLOAT logMaxScale;// log(MaxScale)
extern MFLOAT logMaxUR;// log(MaxUR)
extern double Smin;// smallest possible interval sizing error = sqrt(min(SF[0]^2, SF[0]^2 + SD[0]^3|SD[0]|/(4*SR[0]^2)))
extern double MaxFAG;// MaxFA/(Smin * RefPen * sqrt(0.5))
extern MFLOAT logMaxFAG;// log(MaxFAG)
extern double MaxFAG2;// MaxFA/RefPen
extern MFLOAT logMaxFAG2;// log(MaxFAG2)

// HERE HERE : convert Cprtab into FLAT_MEMORY arrays indexed as <field>(I,K) == <field>_base[I * PRsiz + K], 
class Cprtab {// array indexed by K (or n-1) and I : NOTE that K index range is typically smaller than for n-1, but the array is allocated for the larger index range wasting memory. Fix with FLAT_MEMORY arrays
 public:
  MFLOAT PM;/* PM(0,I,K,Y) */
  MFLOAT Pr;/* Pr(Y[I]-Y[I-n]) */
};

extern double *PrBiasWT;/* PrBiasWT[I = 0..N] = -PRbiasWT * biasWT * sum(i=0..I) Pr(yi) * log(Pr(yi)) + (1-Pr(yi)) * log(1-Pr(yi)), where yi == Y[i+1]-Y[i] */

extern double AlignResMax;
extern double aNaN;
extern MFLOAT biasWToutlierF;
extern MFLOAT biasWToutlierM1;
extern MFLOAT E2var;
extern MFLOAT EBiasWT;
extern int Extend;
extern MFLOAT F2var;
extern MFLOAT FBiasWT;
extern MFLOAT FBiasWTm[];
extern MFLOAT *FBiasWTmRev;
extern double *FnPow;
extern MFLOAT *FnPowf;
extern double FpPow[];
extern MFLOAT Frate;
extern MFLOAT FrateBiasWT;
extern int global;
extern int globalfallback;
extern int gN;
extern lightweight_heap *heap;
extern double InvRtheta;
extern MFLOAT InvRthetaF;
extern MFLOAT InvRtheta2;
extern double IresSD;
extern float IresSDf;
extern int kmax;
extern int *Kmax;
extern int KMAX_APPROX;
extern MFLOAT *LogFnPow;
extern MFLOAT LogFpPow[];
extern MFLOAT *LogFpPowRev;

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT 
 extern __m256 ZeroHigh_mm256[DELTA_MAX + EXTEND + 1];
 extern __m256 v256_F2var, v256_S2var, v256_R2var, v256_E2var, v256_RefPen, v256_zero, v256_0p5, v256_1p0, v256_1p5, v256_2p0;
 extern __m256 v256_biasWToutlierF, v256_biasWToutlierM1;
 extern __m256 v256_Frate, v256_PoutlierM1, v256_FrateBiasWT;
 extern __m256 v256_OutlierLambdaInv, v256_PoutlierXthetaBYlambda, v256_PoutlierXthetaBYmax, v256_PoutlierPIbXtheta, v256_PoutlierPIBYmax;
 extern __m256 v256_outlierNegMaxBC;

 extern __m256 v256_logMaxScale, v256_logMaxUR, v256_logMaxFAG, v256_logMaxFAG2;
 extern __m256 v256_RthetaF, v256_InvRthetaF, v256_InvRtheta2;

 extern __m256d v256d_0p5, v256d_1p0;

#elif VECTOR_MIC && ( USE_MIC || USE_AVX512 ) && USE_MFLOAT
 extern __mmask16 MaskLow_mask16[DELTA_MAX + EXTEND + 1];

 extern __m512 v512_F2var, v512_S2var, v512_R2var, v512_E2var, v512_RefPen, v512_1p0, v512_0p5, v512_1p5, v512_2p0;
 extern __m512d v512d_0p5;
 extern __m512 v512_biasWToutlierF, v512_biasWToutlierM1;
 extern __m512 v512_Frate, v512_PoutlierM1, v512_FrateBiasWT;
 extern __m512 v512_OutlierLambdaInv, v512_PoutlierXthetaBYlambda, v512_PoutlierXthetaBYmax, v512_PoutlierPIbXtheta, v512_PoutlierPIBYmax;
 extern __m512 v512_outlierNegMaxBC, v512_outlierMaxBC;

 extern __m512 v512_logMaxScale, v512_logMaxUR, v512_logMaxFAG, v512_logMaxFAG2;
 extern __m512 v512_LOG2E, v512_MIN_EXP2ARG, v512_MAX_EXP2ARG;
 extern __m512 v512_RthetaF, v512_InvRthetaF;
#endif

extern int Lrange;
extern double maxdist;


extern double minKB;
extern MFLOAT OutlierLambdaInv;
extern MFLOAT outlierMaxBC;
extern MFLOAT PoutlierEndM1;
extern MFLOAT PoutlierM1;
extern MFLOAT PoutlierPIbXtheta;
extern MFLOAT PoutlierPIBYmax;
extern MFLOAT PoutlierXthetaBYlambda;
extern MFLOAT PoutlierXthetaBYmax;
extern int PRNsiz;
extern int PRsiz;
extern Cprtab **PRtab;
extern double pTP;
extern MFLOAT R2var;
extern int RANGE;
extern int RANGE_Y;
extern MFLOAT RefPen;
extern double resKB;
extern float resKBf;
extern double Rtheta;
extern MFLOAT RthetaF;
extern MFLOAT S2var;
extern MFLOAT *VARtab;
extern long long VmSize, VmRSS, VmSwap, VmPeak, VmHWM;
extern double Xtheta;
extern double Ylambda;

extern int localINITIAL_DELTA_RANGE;/* see HaploTypeR() : ((24 << DELTA_OVERSAMPLE) + 1) <= INITIAL_DELTA_RANGE */

void score_initOutlier();

 class CdoubleBlock {/* a block of up to 4G doubles */
 public:
   double *block;/* pointer to start of memory block (allocated with malloc()) */
   unsigned int siz;/* size of memory block in doubles */

   unsigned short RefLoc;/* only define if swapout != 0 : current contig location (rounded to nearest multiple of 100kb) when this block was swapped out. 
			    Useful with MPROBEVAL_MAPORDER to allow SwappedAllocations_restore() to read in blocks from right to left (Last Out First In) order */
   short swapout;// 1 IFF this block is currently swapped out to FILE fp starting at specified offset, 0 if block is in memory but not on disk, -1 if block is both in memory and on disk (swapped back in)
   int tid;// unique id for this block's fp (currently the tid that created this fp): can be used to multithread subsequent block read & write operations

   FILE *fp;/* only defined if swapout != 0 (see swapout) */
   size_t offset;/* only defined if swapout != 0 (see swapout) */

   void writeout(unsigned int tsiz, FILE *Ifp, size_t Ioffset, double IRefLoc, char *filename, int Itid){
     fp = Ifp;
     offset = Ioffset;
     tid = Itid;

     /* write out current block to disk using FILE pointer fp at offset bytes : fp and current block are assumed to not be accessible to any other thread during this call
	NOTE : tsiz is the size of the remaining unused part of the memory block : block[siz-tsiz .. siz-1] */
     siz -= tsiz;
     size_t n = fwrite_LOOP(block, sizeof(double), siz, fp);
     if(n != siz){
       #pragma omp critical
       {
	 printf("fwrite= %lu: write of %u doubles to %s failed(tid=%d)\n", n, siz, filename,Itid);
	 fflush(stdout);exit(1);
       }
     }
     swapout = 1;
     RefLoc = min(64000.0, IRefLoc * 0.01);/* as multiple of 100kb */
     if(VERB>=3 /* && tid==0 */){
       #pragma omp critical
       {
	 printf("tid=%d:Wrote out %u/%u doubles to %s (fp=%p,fd=%d) at offset=%lu bytes:siz=%u: cum wall= %0.6f secs\n",tid,siz,siz+tsiz,filename,fp,fileno(fp),offset*sizeof(double),siz,wtime());
	 fflush(stdout);
       }
     }
   };

   void readin() {/* read in the current block from disk : assumed to be called in sequential code or critical section of parallel code OR each thread reads only blocks for the same "tid" value 
		     (no locks on fp or this memory block are used) */
     if(DEBUG) assert(swapout > 0 && fp != NULL);
     if(fseek(fp, offset * sizeof(double), SEEK_SET)){
       int eno = errno;
       char *err = strerror(eno);
       printf("CdoubleBlock::readin():fseek() for fd=%d to offset=%lu bytes failed (block=%p,RefLoc=%d):errno=%d:%s\n",fileno(fp),offset * sizeof(double),block,(int)RefLoc,eno,err);
       fflush(stdout);exit(1);
     }

     size_t len = fread(block, sizeof(double), siz, fp);
     if(len != siz){
       int eno = errno;
       char *err = strerror(eno);
       printf("CdoubleBlock::readin():fread() of %u doubles from fd=%d at offset=%lu bytes failed (ret=%lu,block=%p,RefLoc=%d):errno=%d:%s\n",siz,fileno(fp),offset * sizeof(double),len,block,(int)RefLoc,eno,err);
       fflush(stdout);exit(1);
     }

     swapout = -1;
   };
 };

 double mprobeval(int n, double *Hcuts, int N, double *origY, int MD, int *MX, double **X, int lc, int rc, int **map, int **mapK, Csetlimit *limit, int *nmapMD, double *TBmapWT, double *newLP, double *oldLP, double Ymin, double Ymax, int addcnt, double *addloc, int *D, double **delta, double *LPdelG, double *LPaddG, double *LPgrad[2], double **LPdelta, double **newLPd, double *newLPd0, int **newLPdPr, int *DminM, int *DmaxM, double **newLPa, int **newLPaPr, int *TminM, int *TmaxM, double **newLPdelta, int *Dcum, int *AminM, int *AmaxM, Ccontig *pcontig, char *skipmap, int Allele, CdoubleBlock **ppBlock = 0);
		 
 double qprobeval(int n, double *Hcuts, int N, double *Y, int MD, int *MX, double **X, int lc, int rc, int **map, int **mapK, Csetlimit *limit, int **nmap, int **nmapK, double *TBmapWT, double *newLP, double *mapWT, double *oldLP, double Ymin, double Ymax, double *LogPV, int **poutlier, Ccontig *pcontig, int *forcemap = 0);

#include <sys/mman.h>

 class CdoubleAllocation {
 public:
   unsigned int mcnt;// number of blocks available (not initialized) : blocks[0..mcnt-1]
   unsigned int cnt;// number of blocks used and initialized/allocated : blocks[0..cnt-1]
   CdoubleBlock *blocks;// NOTE : If tid >= 0 : first block is part of a larger block shared by all threads that is located on global DoubleAllocations vector
   
   double *threadmem;/* threadmem[0..tsiz-1] is remainder of last block blocks[cnt-1].block */
   unsigned int tsiz;/* size of threadmem[0..tsiz-1] : same as blocks[cnt-1].block[siz-tsiz .. siz-1], where siz == block[cnt-1].siz */
   int tid;/* current OpenMP thread id */
   int Allele;/* 0, 1 or 2 : If > 0, blocks[(tid>=0 ? 1 : 0) .. cnt-2] will have been swapped out and can be restored by calling readin()*/

   FILE *fp;/* If != NULL : open FILE pointer that has been used to save a copy of blocks[0..cnt-2] (all except the last block) */
   size_t offset;/* offset to end of FILE as a multiple of sizeof(double) : sum of blocks[0..cnt-2].siz */
   char *fbuf;// binary buffer for fp
   char filename[PATH_MAX];/* name of file that has been used to save a copy of blocks[0..cnt-2] */

   CdoubleAllocation(){
     mcnt = cnt = tsiz = 0;
     blocks = NULL;
     threadmem = NULL;
     tid = -1;
     fp = NULL;
     fbuf = NULL;
   };
   CdoubleAllocation(int Itid, int IAllele, double *Ithreadmem, unsigned int Itsiz, char *prefix, long long id){
     tid = Itid;
     Allele = IAllele;
     threadmem = Ithreadmem;
     tsiz = Itsiz;
     mcnt = 1023;
     blocks = new CdoubleBlock[mcnt];
     cnt = 0;
     if(threadmem){
       if(DEBUG) assert(tid >= 0);
       cnt = 1;
       blocks[0].block = threadmem;
       blocks[0].siz = tsiz;
       blocks[0].swapout = 0;
     }

     fp = NULL;
     fbuf = NULL;
     filename[0] = '\0';

     /* construct unique binary filename : file will be created only if needed */
     if(DEBUG) assert(prefix);
     sprintf(filename,"%s_p%lld_id%lld_A%d_t%d.mpbin",prefix,(long long)pid, id,Allele,tid);
     if(DEBUG) assert(strlen(filename) < PATH_MAX-1);     

     //     (void)unlink(filename);// unlink any previous version of the file
   }
   
   void fileopen() {
     /* initialize output binary file for swapping memory blocks of current thread */
     if(DEBUG) assert(filename[0]);

     if((fp = fopen(filename,"w+")) == NULL){
       int eno = errno;
       char *err = strerror(eno);
       #pragma omp critical
       {
         printf("failed to open file %s for writing:errno=%d:%s\n",filename,eno,err);
         fflush(stdout);exit(1);
       }
     }
     offset = 0;
     
     size_t fbufsiz = VMEM_MINIMIZE * 1024 * sizeof(double);
     fbuf = new char[fbufsiz];
     setbuffer(fp,fbuf,fbufsiz);
   };

   // NOTE : the functions readin() is assumed to be called from within a critical section or non-parallel section
   void readin(){// read in all paged out blocks (NOTE : last block is never paged out and first block is typically not paged out since it is part of a larger block shared by all threads)
     if(DEBUG) assert(Allele != 0);
     if(DEBUG) assert(fp != NULL);

     size_t siz = 0;
     for(int i = cnt-1; --i >= 0; ){
       CdoubleBlock *p = &blocks[i];
       if(p->swapout > 0){
	 p->readin();
	 siz += p->siz;
       }
     }
     if(VERB>=2 && cnt > 1){
       printf("tid=%d:read in %u/%u blocks totalling %lu doubles from %s\n",tid, cnt-2,cnt,siz,filename);
       fflush(stdout);
     }
   };

   ~CdoubleAllocation(){
     if(blocks){
       if(cnt > 0){// NOTE : If tid >= 0 : first block is part of a larger memory block that must be deallocated by caller
	 for(unsigned int i = (tid >= 0) ? 1 : 0; i < cnt; i++){
	   if(munmap(blocks[i].block, blocks[i].siz * sizeof(double))){
	     int eno = errno;
	     char *err = strerror(eno);
	     printf("munmap(%p,%lu) failed: errno=%d:%s\n",blocks[i].block,blocks[i].siz * sizeof(double), eno, err);
	     dumpmemmap();
	     fflush(stdout);exit(1);
	   }
	 }
       }
       delete [] blocks;
     }

     (void)unlink(filename);
     if(fp){
       (void)fclose(fp);
       delete [] fbuf;
       if(VERB>=3){
	 #pragma omp critical
	 {
	   printf("tid=%d: unlinked swapping file %s\n",tid,filename);
	   fflush(stdout);
	 }
       }
     }
   }
 };

 extern std::vector<double *> DoubleAllocations; // see mprobeval.cpp
 extern std::vector<CdoubleAllocation *> SwappedAllocations; // see mprobeval.cpp
 extern CdoubleBlock **blockList;
 extern int blockListSwapped, blockListLen,blockListMax;

 void SwappedAllocations_ListSort();// see mprobeval.cpp
 void SwappedAllocations_restore();// see mprobeval.cpp
 void SwappedAllocations_DONTNEED(int MinSwapCnt);// see mprobeval.cpp
 int SwappedAllocations_swapcnt();// see mprobeval.cpp

 void DoubleAllocations_free();// see mprobeval.cpp

static inline MFLOAT G(MFLOAT x, MFLOAT y)
{
  MFLOAT err = x-y;

  MFLOAT var = F2var + S2var * y;
  if(QUADRATIC_VARIANCE)
    var += R2var * y * y;
  MFLOAT Ivar = ((MFLOAT)1.0)/var;

  return RefPen * sqrtF(Ivar) * expF(-err*err*Ivar);
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT

static inline float horizontal_add_mm256 (__m256 a)
{
  __m256 t1 = _mm256_hadd_ps(a,a);
  __m256 t2 = _mm256_hadd_ps(t1,t1);
  __m128 t3 = _mm256_extractf128_ps(t2,1);
  __m128 t4 = _mm_add_ss(_mm256_castps256_ps128(t2),t3);
  return _mm_cvtss_f32(t4);        
}

/* AVX256 exp function */
static inline __m256 expF_mm256(__m256 v_X)
{
  float X[8] __attribute__((aligned(32)));
  float R[8] __attribute__((aligned(32)));
  _mm256_store_ps(X, v_X);

  #pragma omp simd aligned(R, X : 32)
  for(int i = 0; i < 8; i++)
    R[i] = expF(X[i]);

  return _mm256_load_ps(R);
}

/* AVX256 erf function */
static inline __m256 erfF_mm256(__m256 v_X)
{
  float X[8];// __attribute__((aligned(32)));// NOTE : alignment seems to hinder the compiler
  float R[8];// __attribute__((aligned(32)));
  _mm256_store_ps(X, v_X);

  #pragma omp simd // aligned(R, X : 32)
  for(int i = 0; i < 8; i++)
    R[i] = erfF(X[i]);

  return _mm256_load_ps(R);
}

/* AVX256 erfc function */
static inline __m256 erfcF_mm256(__m256 v_X)
{
  float X[8];// __attribute__((aligned(32)));// NOTE : alignment seems to hinder the compiler
  float R[8];// __attribute__((aligned(32)));
  _mm256_store_ps(X, v_X);

  #pragma omp simd // aligned(R, X : 32)
  for(int i = 0; i < 8; i++)
    R[i] = erfcF(X[i]);

  return _mm256_load_ps(R);
}

static inline __m256 G_mm256(__m256 x, __m256 y)
{
  __m256 err = x - y; 

  __m256 var = v256_F2var + v256_S2var * y;
  if(QUADRATIC_VARIANCE)
    var += v256_R2var * y * y;

  /* Fast Inverse Square Root + One Newton Raphson iteration */
  __m256 Ivar1 = _mm256_rsqrt_ps(var);
  __m256 Ivar = Ivar1 * (v256_1p5 - v256_0p5 * var * Ivar1 *Ivar1);

  err *= Ivar;

  return v256_RefPen * Ivar * expF_mm256(-err*err);
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

/* convert 8+8 doubles in Ylo, Yhi to 16 float vector, with Ylo values in lower half of vector */
static inline __m512 DoubleToFloat_mm512(__m512d Ylo, __m512d Yhi)
{
#if USE_MIC
  __m512 RetLo = _mm512_cvtpd_pslo(Ylo);
  __m512 RetHi = _mm512_cvtpd_pslo(Yhi);
  __m512 Ret = _mm512_mask_permute4f128_ps(RetLo, 0xff00, RetHi, (0x1 << 6));
#elif USE_AVX512
  __m256 RetLo = _mm512_cvtpd_ps(Ylo);
  __m256 RetHi = _mm512_cvtpd_ps(Yhi);
  __m512 Ret = _mm512_insertf32x8(_mm512_castps256_ps512(RetLo),RetHi, 1);
#endif

#if VECTOR_DEBUG>=2 
  double YloV[8]; _mm512_storeu_pd(YloV, Ylo);
  double YhiV[8]; _mm512_storeu_pd(YhiV, Yhi);
  float RetV[16]; _mm512_storeu_ps(RetV, Ret);
  for(int i = 0; i < 8; i++){
    float y = YloV[i];
    if(!(fabs(y - RetV[i]) <= fabs(y) * 1e-7 + 1e-35)){
      printf("i=%d: Ylo[i]= %0.10e, float(Ylo[i])= %0.10e, RetV[i]= %0.10e\n", i, Ylo[i], y, RetV[i]);
      fflush(stdout);
      assert(fabs(y - RetV[i]) <= fabs(y) * 1e-7 + 1e-35);
    }  
    
    y = YhiV[i];
    if(!(fabs(y - RetV[i+8]) <= fabs(y) * 1e-7 + 1e-35)){
      printf("i=%d: Yhi[i]= %0.10e, float(Yhi[i])= %0.10e, RetV[i+8]= %0.10e\n", i, Yhi[i], y, RetV[i+8]);
      fflush(stdout);
      assert(fabs(y - RetV[i+8]) <= fabs(y) * 1e-7 + 1e-35);
    }  
  }
#endif

  return Ret;
}

// #define LOG2Ef    1.442695040888963407359924681001892137426645954152985934135f // log2(e)
// #define LOG2E_2P24f   24204406.323122970759869f // log2(e) * 2^24

// #define MIN_EXP2ARG -128.0f
// #define MAX_EXP2ARG  0x7f.ffffffp0f

/* AVX512 exp function */
static inline __m512 expF_mm512(__m512 v_X)
{
#if 0
  __m512 v_Xlog2 = _mm512_mul_ps(v_X, v512_LOG2E);
  __m512i v_Xlog2x2p24 = _mm512_cvtfxpnt_round_adjustps_epi32(v_Xlog2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_24);
  return _mm512_exp223_ps(v_Xlog2x2p24);
#else

  float X[16] __attribute__((aligned(64)));// HERE HERE : does the alignment help or hinder the compiler
  float R[16] __attribute__((aligned(64)));
  _mm512_store_ps(X, v_X);

  #pragma omp simd aligned(R, X : 64)
  for(int i = 0; i < 16; i++)
    R[i] = expF(X[i]);

  return _mm512_load_ps(R);
#endif
}

/* AVX512 erf function */
static inline __m512 erfF_mm512(__m512 v_X)
{
  float X[16] __attribute__((aligned(64)));// HERE HERE : check if alignment seems to hinder the compiler
  float R[16] __attribute__((aligned(64)));
  _mm512_store_ps(X, v_X);

  #pragma omp simd aligned(R, X : 64)
  for(int i = 0; i < 16; i++)
    R[i] = erfF(X[i]);

  return _mm512_load_ps(R);
}
/* AVX512 erfc function */
static inline __m512 erfcF_mm512(__m512 v_X)
{
  float X[16] __attribute__((aligned(64)));// HERE HERE : check if alignment seems to hinder the compiler
  float R[16] __attribute__((aligned(64)));
  _mm512_store_ps(X, v_X);

  #pragma omp simd aligned(R, X : 64)
  for(int i = 0; i < 16; i++)
    R[i] = erfcF(X[i]);

  return _mm512_load_ps(R);
}

static inline __m512 rsqrt_mm512(__m512 x)
{
#if USE_MIC
  return _mm512_rsqrt23_ps(x);// 23 bits precision
#elif __AVX512ER__
  return _mm512_rsqrt28_ps(x);// 28 bits precision
#else
  if(DEBUG) assert(USE_AVX512);
  __m512 xInv = _mm512_rsqrt14_ps(x);// 14 bits precision
  return xInv * (v512_1p5 - v512_0p5 * x * xInv * xInv);// Newton Raphson iteration
#endif
}

static inline __m512 G_mm512(__m512 x, __m512 y)
{
  __m512 err = _mm512_sub_ps(x, y); 

  __m512 var = _mm512_add_ps(v512_F2var, _mm512_mul_ps(v512_S2var, y));
  if(QUADRATIC_VARIANCE)
    var = _mm512_add_ps(var, _mm512_mul_ps(v512_R2var, _mm512_mul_ps(y, y)));

  __m512 Ivar = rsqrt_mm512(var);

  err = _mm512_mul_ps(err, Ivar);

  __m512 Vexp = expF_mm512(- _mm512_mul_ps(err,err));
  return _mm512_mul_ps(_mm512_mul_ps(v512_RefPen, Ivar), Vexp);
}

#endif // VECTOR_AVX || VECTOR_MIC

/* see Refinement document for function GE() */
template<int ERFC> // If ERFC==1, use 1-erfc(x) instead of erf(x)
static inline MFLOAT GE(MFLOAT x, MFLOAT y, MFLOAT b, MFLOAT IvarX, int yextend)
{
  if(DEBUG>=2) assert(END_FIX || !yextend);
  if(DEBUG>=2 && !(0.0 <= b && b <= y)){
    printf("\n GE(x=%0.6f,y=%0.6f,b=%0.6f,IvarX=%0.8f,yextend=%d)\n",x,y,b,IvarX,yextend);
    fflush(stdout);
    assert(0.0 <= b && b <= y);
  }
  MFLOAT xy = x - y;
  MFLOAT xyb = xy + b;
  MFLOAT xybIvarX = xyb * IvarX;

  if(yextend)
    return ERFC ? 1.0f - 0.5f * erfcF(xybIvarX) : 0.5f + 0.5f * erfF(xybIvarX);
		
  MFLOAT xyIvarX = xy * IvarX;
  MFLOAT ret = 0.5f * (ERFC ? erfcF(xyIvarX) - erfcF(xybIvarX) : erfF(xybIvarX) - erfF(xyIvarX));

  if((DEBUG>=2 && !(ret >= -(USE_MFLOAT ? 1e-7 : 1e-10))) || (VERB>=3 && rverb >= 2)){
    printf("\tGE(x=%0.6f,y=%0.6f,b=%0.6f,IvarX=%0.8f,yextend=%d): ret=%0.10e, xybIvarX=%0.8f,xyIvarX=%0.8f:\n\t erff(xybIvarX)=%0.10e,erff(xyIvarX)=%0.10e,erfcf(xybIvarX=%0.10e,erfcf(xyIvarX)=%0.10e,erf(xybIvarX)=%0.10e,erf(xyIvarX)=%0.10e\n",
	   x,y,b,IvarX,yextend,ret,xybIvarX,xyIvarX,erff(xybIvarX),erff(xyIvarX),erfcf(xybIvarX),erfcf(xyIvarX),erf(xybIvarX),erf(xyIvarX));
    fflush(stdout);
    assert(ret >= -(USE_MFLOAT ? 1e-7 : 1e-10));
  }
  return max(((MFLOAT)0.0), ret);
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT

/* AVX256 version of GE() */
template<int ERFC> // If ERFC==1, use 1-erfc(x) instead of erf(x) 
static inline __m256 GE_mm256(__m256 x, __m256 y, __m256 b, __m256 IvarX, int yextend)
{
  __m256 xy = x - y;
  __m256 xyb = xy + b;
  __m256 xybIvarX = xyb * IvarX;

  __m256 ret, xyIvarX;
  if(yextend)
    ret = ERFC ? v256_1p0 - v256_0p5 * erfcF_mm256(xybIvarX) : v256_0p5 + v256_0p5 * erfF_mm256(xybIvarX);
  else {
    xyIvarX = xy * IvarX;
    __m256 ret1 = v256_0p5 * ( ERFC ? erfcF_mm256(xyIvarX) - erfcF_mm256(xybIvarX) : erfF_mm256(xybIvarX) - erfF_mm256(xyIvarX) );
    ret = _mm256_max_ps( _mm256_setzero_ps(), ret1);
  }

#if VECTOR2_DEBUG >= 2
  float retV[8]; _mm256_storeu_ps(retV, ret);
  float xV[8]; _mm256_storeu_ps(xV, x);
  float yV[8]; _mm256_storeu_ps(yV, y);
  float bV[8]; _mm256_storeu_ps(bV, b);
  float IvarXV[8]; _mm256_storeu_ps(IvarXV, IvarX);
  for(int i = 0; i < 8; i++){
    float ge = GE<ERFC>(xV[i], yV[i], bV[i], IvarX[i], yextend);
    if(!(fabs(ge - retV[i]) < fabs(ge) * 1e-5 + 1e-9)){
      float xybIvarXV[8]; _mm256_storeu_ps(xybIvarXV, xybIvarX);
      float xyIvarXV[8]; _mm256_storeu_ps(xyIvarXV, xyIvarX);
      printf("i=%d:GE_mm256(x,y,b,IvarX,yext)[i]= %0.8e, GE(x[i],y[i],b[i],IvarX[i])= %0.8e\n",
        i, retV[i], ge);
      printf("\t x[i]= %0.6f, y[i]= %0.6f, b[i]= %0.6f, IvarX[i]= %0.8e, yext= %2d: (x-y+b)*IvarX= %0.8e, (x-y)*IvarX= %0.8e\n", 
        xV[i], yV[i], bV[i], IvarXV[i], yextend, (xV[i]-yV[i]+bV[i])*IvarXV[i], (xV[i]-yV[i])*IvarXV[i]);
      printf("\t erff((x-y+b)*IvarX)= %0.10e, erff((x-y)*IvarX)= %0.10e, erf((x-y+b)*IvarX)= %0.10e, erf((x-y)*IvarX)= %0.10e\n",
        erfF((xV[i]-yV[i]+bV[i])*IvarXV[i]), erfF((xV[i]-yV[i])*IvarXV[i]), erf((xV[i]-yV[i]+bV[i])*IvarXV[i]), erf((xV[i]-yV[i])*IvarXV[i]));
      printf("\t xybIvarX[i]= %0.8e, xyIvarX[i]= %0.8e, erf(xybIvarX[i])= %0.10e, erf(xyIvarX[i])= %0.10e\n",
        xybIvarXV[i], xyIvarXV[i], erfF(xybIvarXV[i]), erfF(xyIvarXV[i]));
      fflush(stdout);
      assert(fabs(ge - retV[i]) < fabs(ge) * 1e-5 + 1e-9);
    }
  }
#endif

  return ret;
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

/* AVX512 version of GE() */
template<int ERFC> // If ERFC==1, use 1-erfc(x) instead of erf(x) 
static inline __m512 GE_mm512(__m512 x, __m512 y, __m512 b, __m512 IvarX, int yextend)
{
  __m512 xy = _mm512_sub_ps(x, y);
  __m512 xyb = _mm512_add_ps(xy, b);
  __m512 xybIvarX = _mm512_mul_ps(xyb, IvarX);

  __m512 ret, xyIvarX;
  if(yextend)
    ret = ERFC ? _mm512_sub_ps(v512_1p0, _mm512_mul_ps(v512_0p5, erfcF_mm512(xybIvarX)))
      : _mm512_add_ps(v512_0p5, _mm512_mul_ps(v512_0p5, erfF_mm512(xybIvarX)));
  else {
    xyIvarX = _mm512_mul_ps(xy, IvarX);
    __m512 ret1 = _mm512_mul_ps(v512_0p5, _mm512_sub_ps(ERFC ? erfcF_mm512(xyIvarX) : erfF_mm512(xybIvarX),
							ERFC ? erfcF_mm512(xybIvarX) : erfF_mm512(xyIvarX)));
    ret = _mm512_max_ps( _mm512_setzero_ps(), ret1);
  }

#if VECTOR2_DEBUG >= 2 || DEBUG>=2
  float retV[16]; _mm512_storeu_ps(retV, ret);
  float xV[16]; _mm512_storeu_ps(xV, x);
  float yV[16]; _mm512_storeu_ps(yV, y);
  float bV[16]; _mm512_storeu_ps(bV, b);
  float IvarXV[16]; _mm512_storeu_ps(IvarXV, IvarX);
  for(int i = 0; i < 16; i++){
    float ge = GE<ERFC>(xV[i], yV[i], bV[i], IvarX[i], yextend);
    if(!(fabs(ge - retV[i]) < fabs(ge) * 1e-4 + 1e-7)){
      float xybIvarXV[16]; _mm512_storeu_ps(xybIvarXV, xybIvarX);
      float xyIvarXV[16]; _mm512_storeu_ps(xyIvarXV, xyIvarX);
      printf("i=%d:GE_mm512(x,y,b,IvarX,yext)[i]= %0.8e, GE(x[i],y[i],b[i],IvarX[i])= %0.8e\n",
        i, retV[i], ge);
      printf("\t x[i]= %0.6f, y[i]= %0.6f, b[i]= %0.6f, IvarX[i]= %0.8e, yext= %2d: (x-y+b)*IvarX= %0.8e, (x-y)*IvarX= %0.8e\n", 
        xV[i], yV[i], bV[i], IvarXV[i], yextend, (xV[i]-yV[i]+bV[i])*IvarXV[i], (xV[i]-yV[i])*IvarXV[i]);
      printf("\t erff((x-y+b)*IvarX)= %0.10e, erff((x-y)*IvarX)= %0.10e, erf((x-y+b)*IvarX)= %0.10e, erf((x-y)*IvarX)= %0.10e\n",
        erfF((xV[i]-yV[i]+bV[i])*IvarXV[i]), erfF((xV[i]-yV[i])*IvarXV[i]), erf((xV[i]-yV[i]+bV[i])*IvarXV[i]), erf((xV[i]-yV[i])*IvarXV[i]));
      printf("\t xybIvarX[i]= %0.8e, xyIvarX[i]= %0.8e, erf(xybIvarX[i])= %0.10e, erf(xyIvarX[i])= %0.10e\n",
        xybIvarXV[i], xyIvarXV[i], erfF(xybIvarXV[i]), erfF(xyIvarXV[i]));
      fflush(stdout);
      assert(fabs(ge - retV[i]) < fabs(ge) * 1e-4 + 1e-7);
    }
  }
#endif

  return ret;
}

#endif // VECTOR_AVX || VECTOR_MIC

/* version of GE without precomputed IvarX */
template<int ERFC> // If ERFC==1, use 1-erfc(x) instead of erf(x) 
static inline MFLOAT GE(MFLOAT x, MFLOAT y, MFLOAT b, int yextend)
{
  if(DEBUG>=2) assert(END_FIX || !yextend);
  if(DEBUG>=2 && !(0.0 <= b && b <= y)){
    printf("\nGE(x=%0.6f,y=%0.6f,b=%0.6f,yextend=%d)\n",x,y,b,yextend);
    fflush(stdout);
    assert(0.0 <= b && b <= y);
  }
  MFLOAT var = F2var + S2var * x;
  if(QUADRATIC_VARIANCE)
    var += R2var * x * x;
  MFLOAT IvarX = ((MFLOAT)1.0)/sqrtF(((MFLOAT)2.0)*var);
  return GE<ERFC>(x,y,b,IvarX, yextend);
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT && DEBUG < 2

/* AVX256 version of GE() without precompute IvarX */
template<int ERFC> // If ERFC==1, use 1-erfc(x) instead of erf(x) 
static inline __m256 GE_mm256(__m256 x, __m256 y, __m256 b, int yextend)
{
  __m256 var = v256_F2var + v256_S2var * x;
  if(QUADRATIC_VARIANCE)
    var += v256_R2var * x * x;
  var *= v256_2p0;

  /* Fast Inverse Square Root with One Newton Rapson iteration */
  __m256 Ivar1 = _mm256_rsqrt_ps(var);
  __m256 IvarX = Ivar1 * (v256_1p5 - v256_0p5 * var * Ivar1 *Ivar1);

  return GE_mm256<ERFC>(x,y,b,IvarX, yextend);
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

/* AVX512 version of GE() without precompute IvarX */
template<int ERFC> // If ERFC==1, use 1-erfc(x) instead of erf(x) 
static inline __m512 GE_mm512(__m512 x, __m512 y, __m512 b, int yextend)
{
  __m512 var = _mm512_add_ps(v512_F2var, _mm512_mul_ps(v512_S2var, x));
  if(QUADRATIC_VARIANCE)
    var = _mm512_add_ps(var, _mm512_mul_ps(v512_R2var, _mm512_mul_ps(x, x)));
  var = _mm512_mul_ps(var, v512_2p0);

  __m512 IvarX = rsqrt_mm512(var);

  return GE_mm512<ERFC>(x,y,b,IvarX, yextend);
}

#endif // VECTOR_AVX || VECTOR_MIC

/* Pn(y) == 1-Pr(y) */
static inline double Pn(register double y)
{
  return 0.5*erfc((y - resKB)*IresSD);
}

static inline float Pn(register float y)
{
  return 0.5*erfcf((y - resKB)*IresSD);
}

static inline double Pr(double y)
{
  return 0.5*(1.0+erf((y - resKB)*IresSD));
}

static inline float Pr(float y)
{
  return 0.5f * (1.0f + erff((y - resKBf) * IresSDf));
}

#pragma omp declare simd uniform(resKBf, IresSDf)
static inline float Pr(float y, float resKBf, float IresSDf)
{
  return 0.5f * (1.0f + erff((y - resKBf) * IresSDf));
}

/* FA() corresponsd to Sint() in RGentigRefScore */

/* compute FA() */
static inline MFLOAT FA(MFLOAT X,
			MFLOAT Y,
			int m,
			int n,
			int J,
			int I,
			int K,/* right end of interval */
			int T,/* left end of interval : only used with RES_VARIANCE */
			int N,
			MFLOAT resR2,/* VARtab[K*N + I] : only used with RES_VARIANCE */
			MFLOAT resL2,/* VARtab[T*N + I-K-n] : only used with RES_VARIANCE */
			MFLOAT FBiasWTxm,
			MFLOAT LogFpPowm1,
			MFLOAT LogFnPown1,
			MFLOAT PrBiasWT,/* PrBiasWT[I] - PrBiasWT[I-K-n] */ // NOT USED
			MFLOAT PRikn,/* PRtab[n-1][I-K].Pr */
			double *H,/* Reference array : debugging only */
			const int outlier, /* != 0 : check for outliers */
			const int OUT_TYPE, /* value of global OUTLIER_TYPE */
			const int OUT_OLD)   /* value of global OUTLIER_OLD(outlierMaxBC) */
{
  if(DEBUG>=2 && !(n>=1 && n <= max(DELTA_Y,deltaExtYRef) + (DELTA_MIN ? DELTA_LIM : 0) +1/*for deleted site*/)){
    #pragma omp critical
    {
      printf("\nFA:X=%0.6f,Y=%0.6f,m=%d,n=%d,J=%d,I=%d,K=%d,outlier=%d:gN=%d,DELTA_Y=%d,deltaExtYRef=%d,DELTA_MIN=%d,DELTA_LIM=%d\n",
	     X,Y,m,n,J,I,K,outlier,gN,DELTA_Y,deltaExtYRef,DELTA_MIN,DELTA_LIM);
      fflush(stdout);
      assert(n>=1 && n <= max(DELTA_Y,deltaExtYRef) + (DELTA_MIN ? DELTA_LIM : 0) +1/*for deleted site*/);
    }
  }
  if(DEBUG>=2) assert(0 <= K && K <= kmax);
  if(DEBUG>=2 && !(I-K-n >= 1 && I <= gN)){
    printf("FA:X=%0.3f,Y=%0.3f,m=%d,n=%d,J=%d,I=%d,K=%d,outlier=%d:gN=%d\n",
	   X,Y,m,n,J,I,K,outlier,gN);
    fflush(stdout);
    assert(I-K-n >= 1 && I <= gN);
  }
  if(DEBUG>=2 && !(PRikn >= 0.0 && fabs(PRikn - Pr(H[I-K]-H[I-K-n])) <= PRikn * (USE_MFLOAT ? 1e-6 : 1e-10) + 1e-35)){
    #pragma omp critical
    {
      printf("\nFA:X=%0.6f,Y=%0.6f,m=%d,n=%d,J=%d,I=%d,K=%d,outlier=%d,OUT_TYPE=%d,OUT_OLD=%d,gN=%d:PRikn=%0.10e,log(PRikn)=%0.10f,H[I-K]-H[I-K-n]=%0.10f,Pr(H[I-K]-H[I-K-n])=%0.10e,log(Pr())=%0.10f\n",
	     X,Y,m,n,J,I,K,outlier,OUT_TYPE,OUT_OLD,gN,PRikn,log(PRikn),H[I-K]-H[I-K-n],Pr(H[I-K]-H[I-K-n]),log(Pr(H[I-K]-H[I-K-n])));
      fflush(stdout);
      assert(PRikn >= 0.0 && fabs(PRikn - Pr(H[I-K]-H[I-K-n])) <= PRikn * (USE_MFLOAT ? 1e-6 : 1e-10) + 1e-35);
    }
  }
  if(DEBUG>=2) assert(m <= max(DELTA_X,deltaExtXRef)+EXTEND_MAX+1);

  MFLOAT var = F2var + S2var * Y;
  if(QUADRATIC_VARIANCE)
     var += R2var * Y * Y;
  if(RES_VARIANCE){
    if(DEBUG>=2 && !(I-K-n-T >= 1)){
      printf("FA:X=%0.3f,Y=%0.3f,m=%d,n=%d,J=%d,I=%d,K=%d,T=%d,N=%d,outlier=%d:gN=%d,kmax=%d\n",
	     X,Y,m,n,J,I,K,T,N,outlier,gN,kmax);
      fflush(stdout);
      assert(I-K-n-T >= 1);
    }
    if(DEBUG>=2) assert(0 <= T && T <= kmax);
    if(DEBUG>=2){
      double resR = H[I] - H[I-K];
      double resL = H[I-K-n] - H[I-K-n-T];
      MFLOAT resvar = resL*resL + resR*resR;
      if(!(fabs(resR2+resL2 - resvar) <= 1e-10 + 1e-6 * resvar)){
        printf("FAIK:I=%d,K=%d,n=%d,T=%d:H[I-K,I]= %0.6f,%0.6f,H[I-K-n,I-K-n-T]=%0.6f,%0.6f, resvar=%0.10f,resR2=%0.10f(true=%0.10f),resL2=%0.10f(true=%0.10f)\n",
          I,K,n,T,H[I-K],H[I],H[I-K-n],H[I-K-n-T],resvar,resR2,resR*resR,resL2,resL*resL);
	fflush(stdout);
	assert(fabs(resR2+resL2 - resvar) <= 1e-10 + 1e-6 * resvar);
      }
    }
    var += E2var * (resR2 + resL2);
  }

  if(DEBUG>=2) assert(var > 0.0);
  MFLOAT Ivar = sqrtF(1.0f/var);
  if(DEBUG>=2) assert(isfinite(Ivar));
  if(DEBUG>=2) assert(PRikn <= 1.0);
  MFLOAT FAik = /* FnPowf[n-1] * */ RefPen * PRikn;// NOTE : FnPow[n-1] has been moved into expF(LogFnPow[n-1])

  MFLOAT err = (X-Y) * Ivar;
  MFLOAT Bias = X * FrateBiasWT /* + PrBiasWT */ + FBiasWTxm;
  MFLOAT OutPen = -fabsF(X-Y) * OutlierLambdaInv;
  if(DEBUG>=2 && biasWT == 0.0) assert(Bias == 0.0);

  MFLOAT err2 = err*err;
  MFLOAT logerrA = Bias + LogFpPowm1;// WAS8 - X*Frate;
  if(!FRATE_FIX)
    logerrA -= X*Frate;
  MFLOAT logerr = logerrA - err2;
  if(FRATE_FIX)
    logerr -= X*Frate;

  // NOTE MaxFA cannot be applied to Fn ^^ (n-1), since mprobeval assumes FAJ changes by exactly a multiple of Fn if n changes by 1 
  MFLOAT logerrb = LogFnPown1 + min(logerr, logMaxFAG);// avoid overflow for large intervals, especially with -biasWT 0

  MFLOAT Vexp = expF(logerrb);
  if(DEBUG>=2) assert(isfinite(Vexp) && Vexp <= MaxFAG);
  MFLOAT LR = FAik * Ivar * Vexp;
  if(DEBUG>=2) assert(isfinite(LR) && LR <= MaxFA);

  if(DEBUG>=2 && rverb >= 2){
    int origrverb = rverb;
    rverb = 0;
    printf("FA:X=%0.10f,Y=%0.10f,m=%d,n=%d,J=%d,I=%d,K=%d,T=%d,N=%d,OUT_TYPE=%d,OUT_OLD=%d,outlier=%d:\n",X,Y,m,n,J,I,K,T,N,OUT_TYPE,OUT_OLD,outlier);
    printf("\t var=%0.10e,resR2=%0.10e,resL2=%0.10e,F2var=%0.10e,S2var=%0.10e,R2var=%0.10e,E2var=%0.10e\n",var,resR2,resL2,F2var,S2var,R2var,E2var);
    printf("\t Ivar=%0.8e,FAik=%0.8e,Bias=%0.8e,OutPen=%0.8e,err=%0.8e,err2=%0.8e,logerrA=%0.8e,logerr=%0.8e,logerrb=%0.8e,Vexp=%0.8e,LR=%0.8e\n",Ivar,FAik,Bias,OutPen,err,err2,logerrA,logerr,logerrb,Vexp,LR);
    fflush(stdout);
    rverb = origrverb;
  }

  if(outlier && OUTLIER_DELTA(X-Y)){
    if(!OUT_TYPE){
      MFLOAT LRout = expF(OutPen + Bias * biasWToutlierF);
      if(DEBUG>=2 && rverb >= 2){
	int origrverb = rverb;
        rverb = 0;
	printf("\t LRout=%0.8e,PoutlierXthetaBYmax=%0.8e,PoutlierM1=%0.8e\n",
            LRout, PoutlierXthetaBYmax,PoutlierM1);
	printf("\t FrateBiasWT= %0.8e, FBiasWTxm= %0.8e, Frate= %0.8e, FBiasWT= %0.8e,biasWToutlier= %0.8f\n",FrateBiasWT,FBiasWTxm,Frate,FBiasWT,biasWToutlier);
	fflush(stdout);
	rverb = origrverb;
      }
      if(DEBUG>=2) assert(isfinite(LRout) && LRout <= MaxFA);
      if(OUT_OLD)
	return LRout * PoutlierXthetaBYlambda + LR * PoutlierM1;
      else
	return LRout * PoutlierXthetaBYmax + LR * PoutlierM1;
    } else {
      MFLOAT logerr2 = logerrA - Bias * biasWToutlierM1 /* NEW */ + OutPen;
      MFLOAT logerr2b = LogFnPown1 + min(logerr2, logMaxFAG2);// avoid overflow for large intervals
      MFLOAT Vexp2 = expF(logerr2b);
      if(DEBUG>=2 && !(isfinite(Vexp2) && Vexp2 <= MaxFAG2 * 1.001)){
	#pragma omp critical
	{
	  printf("FA:X=%0.10f,Y=%0.10f,m=%d,n=%d,J=%d,I=%d,K=%d,T=%d,N=%d,OUT_TYPE=%d,OUT_OLD=%d:\n",X,Y,m,n,J,I,K,T,N,OUT_TYPE,OUT_OLD);
	  printf("\t var=%0.10e,resR2=%0.10e,resL2=%0.10e,F2var=%0.10e,S2var=%0.10e,R2var=%0.10e,E2var=%0.10e\n",var,resR2,resL2,F2var,S2var,R2var,E2var);
	  printf("\t Ivar=%0.8e,FAik=%0.8e,Bias=%0.8e,OutPen=%0.8e,err=%0.8e,err2=%0.8e,logerr=%0.8e,Vexp=%0.8e,LR=%0.8e\n",Ivar,FAik,Bias,OutPen,err,err2,logerr,Vexp,LR);
	  printf("\t logerr2= %0.8e, logerr2b= %0.8e, Vexp2= %0.8e, logMAXFAG2= %0.8f, MaxFAG2= %0.8e, LogFnPown1=%0.8e\n",logerr2,logerr2b,Vexp2,logMaxFAG2,MaxFAG2, LogFnPown1);
	  fflush(stdout);
	  assert(isfinite(Vexp2) && Vexp2 <= MaxFAG2 * 1.001);
	}
      }
      MFLOAT LRout = FAik * Vexp2;
      if(DEBUG>=2 && rverb >= 2){
	int origrverb = rverb;
        rverb = 0;
	printf("\t loggerr2=%0.8e,Vexp2=%0.8e,LRout=%0.8e,PoutlierPIBYmax=%0.8e,PoutlierM1=%0.8e\n",
            logerr + err2 - Bias * biasWToutlierM1, Vexp2, LRout, PoutlierPIBYmax,PoutlierM1);
	printf("\t FrateBiasWT= %0.8e, FBiasWTxm= %0.8e, Frate= %0.8e, FBiasWT= %0.8e,biasWToutlier= %0.8f\n",FrateBiasWT,FBiasWTxm,Frate,FBiasWT,biasWToutlier);
	fflush(stdout);
	rverb = origrverb;
      }
      if(DEBUG>=2) assert(isfinite(LRout) && LRout <= MaxFA);

      if(OUT_OLD)
	return LRout * PoutlierPIbXtheta + LR * PoutlierM1;
      else
	return LRout * PoutlierPIBYmax + LR * PoutlierM1;
    }
  } else
    return LR;
}

/* compute the part of FA() that does not depend on X, m or J */
static inline MFLOAT FAIK(MFLOAT Y,
			  MFLOAT PRikn,
			  int n,
			  int I,
			  int K,/* right end of interval */
			  int T, /* left end of interval : only used with RES_VARIANCE */
                          int N, /* only used for debugging */
			  MFLOAT resR2,/* VARtab[K*N + I] : only used with RES_VARIANCE */
			  MFLOAT resL2,/* VARtab[T*N + I-K-n] : only used with RES_VARIANCE */
			  MFLOAT &Ivar,/* return value = 1.0/(F2var + S2var*Y + R2var*Y*Y) */
			  double *H)/* Reference array : only used for debugging */
{
  if(DEBUG>=2) assert(n >= 1 && n <= max(DELTA_Y,deltaExtYRef) + (DELTA_MIN ? DELTA_LIM : 0) +1/*for deleted site*/);
  if(DEBUG>=2) assert(0 <= K && K <= kmax);
  if(DEBUG>=2) assert(I-K-n >= 1 && H[I-K-n] > 0.0 && I <= PRNsiz);
  if(DEBUG>=2) assert(Y > 0.0);
  if(DEBUG>=2 && !(PRikn >= 0.0 && fabs(PRikn - Pr(H[I-K]-H[I-K-n])) <= PRikn * (USE_MFLOAT ? 1e-6 : 1e-10) + 1e-35 )){
    printf("\nFAIK(Y=%0.3f,n=%d,I=%d,K=%d):PRikn=%0.10e,PRtab[n-1][I-K]=%0.10e,Pr(H[I-K]-H[I-K-n])=%0.10e, H[I-K]=%0.3f,H[I-K-n]=%0.3f,PRsiz=%d,PRNsiz=%d\n",
	   Y,n,I,K,PRikn,PRtab[n-1][I-K].Pr,Pr(H[I-K]-H[I-K-n]), H[I-K], H[I-K-n], PRsiz,PRNsiz);
    fflush(stdout);
    assert(PRikn >= 0.0 && fabs(PRikn - Pr(H[I-K]-H[I-K-n])) <= PRikn *(USE_MFLOAT ? 1e-6 : 1e-10) + 1e-35);
  }
  if(DEBUG>=2) assert(PRikn <= 1.0);

  MFLOAT var = F2var + S2var * Y;
  if(QUADRATIC_VARIANCE)
    var += R2var * Y * Y;
  if(RES_VARIANCE){
    if(DEBUG>=2) assert(I-K-n-T >= 1);
    if(DEBUG>=2) assert(0 <= T && T <= kmax);
    if(DEBUG>=2){
      double resR = H[I] - H[I-K];
      double resL = H[I-K-n] - H[I-K-n-T];
      MFLOAT resvar = resL*resL + resR*resR;
      if(!(fabs(resR2+resL2 - resvar) <= 1e-10 + 1e-6 * resvar)){
        printf("FAIK:I=%d,K=%d,n=%d,T=%d:H[I-K,I]= %0.6f,%0.6f,H[I-K-n,I-K-n-T]=%0.6f,%0.6f, resvar=%0.10f,resR2=%0.10f(true=%0.10f),resL2=%0.10f(true=%0.10f)\n",
          I,K,n,T,H[I-K],H[I],H[I-K-n],H[I-K-n-T],resvar,resR2,resR*resR,resL2,resL*resL);
	fflush(stdout);
	assert(fabs(resR2+resL2 - resvar) <= 1e-10 + 1e-6 * resvar);
      }
    }
    var += E2var * (resR2 + resL2);
  }

  if(DEBUG>=2) assert(var > 0.0);
  Ivar = sqrtF(1.0f / var);
  if(DEBUG>=2 && !(isfinite(Ivar) && Ivar <= 0.71/Smin)){
    printf("var=%0.10f,F2var=%0.10f,S2var=%0.10f,R2var=%0.10f,y=%0.8f,SF[0]=%0.10f,SD[0]=%0.10f,SR[0]=%0.10f,Ivar = %0.10f,Smin=%0.10f\n",var,F2var,S2var,R2var,Y,SF[0],SD[0],SR[0],Ivar,Smin);
    fflush(stdout);
    assert(isfinite(Ivar) && Ivar <= 0.71/Smin);
  }

  if(DEBUG>=2 && rverb){
    int origrverb = rverb;
    rverb = 0;
    printf("FAIK:Y=%0.10f,n=%d,I=%d,K=%d,T=%d,N=%d:\n\t var=%0.10e,resR2=%0.10e,resL2=%0.10e,F2var=%0.10e,S2var=%0.10e,R2var=%0.10e,E2var=%0.10e,Ivar=%0.8e\n", 
	   Y,n,I,K,T,N,var,resR2,resL2,F2var,S2var,R2var,E2var,Ivar);
    fflush(stdout);
    rverb = origrverb;
  }

  return /* FnPowf[n-1] * */ RefPen * PRikn;// NOTE : FnPow[n-1] has been moved into FAJ as expF(LogFnPow[n-1])
}

/* version of FAIK() with interval Y[I-K] - Y[I-K-n] increased to newY */
static inline MFLOAT FAIKD(MFLOAT Y, /* new value of Yc(Y,I,K) - Yc(I-K-n,T) */
			   int n,
			   int I,
			   int K,/* right end of interval */
			   int T,/* left end of interval : only used with RES_VARIANCE */
			   MFLOAT newY,/* new assumed value of Y[I-K] - Y[I-K-n] to compute Pr() */
			   MFLOAT &Ivar,/* return value = 1.0/(F2var + S2var*Y + R2var*Y*Y) */
			   double *H)/* Reference array */
  
{
  if(DEBUG>=2) assert(n >= 1 && n <= max(DELTA_Y,deltaExtYRef) + (DELTA_MIN ? DELTA_LIM : 0) +1/*for deleted site*/);
  if(DEBUG>=2) assert(0 <= K && K <= kmax);
  if(DEBUG>=2) assert(I-K-n >= 1 && I <= PRNsiz);
  if(DEBUG>=2) assert(Y > 0.0);
  if(DEBUG>=2) assert(newY > 0.0);

  MFLOAT var = F2var + S2var * Y;
  if(QUADRATIC_VARIANCE)
     var += R2var * Y * Y;
  if(RES_VARIANCE){
    if(DEBUG>=2) assert(I-K-n-T >= 1);
    if(DEBUG>=2) assert(0 <= T && T <= KMAX);
    MFLOAT resR = H[I] - H[I-K];
    MFLOAT resL = H[I-K-n] - H[I-K-n-T];
    MFLOAT resvar = resL*resL + resR*resR;
    var += E2var * resvar;
  }

  if(DEBUG>=2) assert(var > 0.0);
  Ivar = sqrtF(((MFLOAT)1.0)/var);
  if(DEBUG>=2) assert(isfinite(Ivar));

  return /* FnPowf[n-1] *  */ RefPen * Pr(newY);// NOTE : FnPow[n-1] has been moved into FAJ as expF(LogFnPow[n-1])
}

/* version of FAIKD without reference array */
static inline MFLOAT FAIKD(MFLOAT Y, /* new value of Yc(Y,I,K) - Yc(I-K-n,T) */
			   MFLOAT resR2,/* VARtab[T*N + P] == Y[P]-Y[P-T] : only used with RES_VARIANCE */
			   MFLOAT resL2,/* VARtab[K*N + I] == Y[I]-Y[I-K] : only used with RES_VARIANCE */
			   MFLOAT newY,/* new assumed value of Y[I-K] - Y[I-K-n] to compute Pr() */
			   MFLOAT &Ivar)/* return value = 1.0/(F2var + S2var*Y + R2var*Y*Y) */
{
  MFLOAT var = F2var + S2var * Y;
  if(QUADRATIC_VARIANCE)
     var += R2var * Y * Y;
  if(RES_VARIANCE)
    var += E2var * (resR2 + resL2);

  Ivar = sqrtF(1.0f / var);

  return RefPen * Pr(newY);// NOTE : FnPow[n-1] has been moved into FAJ as expF(LogFnPow[n-1])
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT

static inline __m256 Pr_mm256(__m256 v_Y)
{
  float Y8[8] __attribute__((aligned(32)));
  float P8[8] __attribute__((aligned(32)));

  _mm256_store_ps(Y8, v_Y);

  #pragma omp simd aligned(Y8,P8 : 32)
  for(int i = 0; i < 8; i++)
    P8[i] = Pr(Y8[i], resKBf, IresSDf);

  return _mm256_load_ps(P8);
}

static inline __m256d Pr_mm256d(__m256d v_Y)
{
  double Y4[4] __attribute__((aligned(32)));
  double P4[8] __attribute__((aligned(32)));

  _mm256_store_pd(Y4, v_Y);

  #pragma omp simd aligned(Y4,P4 : 32)
  for(int i = 0; i < 4; i++)
    P4[i] = Pr(Y4[i]);

  return _mm256_load_pd(P4);
}

static inline __m256d log_mm256d(__m256d v_Y)
{
  double Y4[4] __attribute__((aligned(32)));
  double P4[8] __attribute__((aligned(32)));

  _mm256_store_pd(Y4, v_Y);

  #pragma omp simd aligned(Y4,P4 : 32)
  for(int i = 0; i < 4; i++)
    P4[i] = log(Y4[i]);

  return _mm256_load_pd(P4);
}

/* AVX256 version of FAIKD() */
static inline __m256 FAIKD_mm256(__m256 v_Y, /* new value of Yc(Y,I,K) - Yc(I-K-n,T) */
                                  MFLOAT resR2,/* VARtab[T*N + P] : only used with RES_VARIANCE */
				  __m256 v_resL2,/* VARtab[K*N + I] : only used with RES_VARIANCE */
				  __m256 v_newY,/* new assumed value of Y[I-K] - Y[I-K-n] to compute Pr() */
				  __m256 &v_Ivar)/* return value = 1.0/(F2var + S2var*Y + R2var*Y*Y + E2var * (resR2 + resL2)) */
  
{
  __m256 v_var = v256_F2var + v256_S2var * v_Y;
#if VECTOR_DEBUG >= 2
  __m256 v_varA = v_var;
#endif

  if(QUADRATIC_VARIANCE)
    v_var += v256_R2var * v_Y * v_Y;

#if VECTOR_DEBUG >= 2
  __m256 v_varB = v_var;
#endif

  if(RES_VARIANCE){
    __m256 v_resR2 = _mm256_set1_ps(resR2);
    v_var += v256_E2var * (v_resR2 + v_resL2);
  }

  /* Fast Inverse Square Root with One Newton Rapson iteration */
  __m256 Ivar1 = _mm256_rsqrt_ps(v_var);
  v_Ivar = Ivar1 * (v256_1p5 - v256_0p5 * v_var * Ivar1 * Ivar1);

  __m256 v_R = v256_RefPen * Pr_mm256(v_newY);

#if VECTOR_DEBUG >= 2
  /* debug results using non-vector version of FAIKD() */
  float Y[8]; _mm256_storeu_ps(Y, v_Y);
  float resL2[8]; _mm256_storeu_ps(resL2, v_resL2);
  float newY[8]; _mm256_storeu_ps(newY, v_newY);
  float var8[8]; _mm256_storeu_ps(var8, v_var);
  float IvarA8[8]; _mm256_storeu_ps(IvarA8, Ivar1);
  float Ivar8[8]; _mm256_storeu_ps(Ivar8, v_Ivar);
  float R8[8]; _mm256_storeu_ps(R8, v_R);
  for(int i = 0; i < 8; i++){
    MFLOAT Ivar;
    MFLOAT R = FAIKD(Y[i], resR2, resL2[i], newY[i], Ivar);
    MFLOAT var = F2var + S2var * Y[i];
    if(QUADRATIC_VARIANCE)
      var += R2var * Y[i] * Y[i];
    if(RES_VARIANCE)
      var += E2var * (resR2 + resL2[i]);

    if(!(Ivar > 0.0 && Ivar8[i] > 0.0 && fabs(Ivar - Ivar8[i]) <= Ivar * 1e-6)){
      #pragma omp critical
      {
	printf("\nFAIKD_mm256:i=%d:Y= %0.6f,resR2= %0.6e,resL2= %0.6e,newY= %0.6f: Ivar= %0.10e(var= %0.10e), Ivar8[i]= %0.10e(var[i]= %0.10e,IvarA= %0.10e),R= %0.6e,R[i]= %0.6e\n",
	       i,Y[i], resR2, resL2[i], newY[i], Ivar, var, Ivar8[i], var8[i], IvarA8[i], R, R8[i]);
	float F2var8[8]; _mm256_storeu_ps(F2var8, v256_F2var);
	float S2var8[8]; _mm256_storeu_ps(S2var8, v256_S2var);
	float R2var8[8]; _mm256_storeu_ps(R2var8, v256_R2var);
	float E2var8[8]; _mm256_storeu_ps(E2var8, v256_E2var);
	printf("\t F2var=%0.7e(%0.7e),S2var=%0.7e(%0.7e),R2var=%0.7e(%0.7e),E2var=%0.7e(%0.7e)\n",F2var,F2var8[i],S2var,S2var8[i],R2var,R2var8[i],E2var,E2var8[i]);
	float varA8[8]; _mm256_storeu_ps(varA8, v_varA);
	float varB8[8]; _mm256_storeu_ps(varB8, v_varB);
	for(int j = 0; j < 8; j++)
	  printf("\t j=%d:Y[j]= %0.6f, F2var[j]= %0.7e, S2var[j]= %0.7e, R2var[j]= %0.7e: varA[j]= %0.7e, varB[j]= %0.7e, var[j]= %0.7e\n", j, Y[j], F2var8[j],S2var8[j],R2var8[j],varA8[j],varB8[j],var8[j]);
	fflush(stdout);
	assert(Ivar > 0.0 && Ivar8[i] > 0.0 && fabs(Ivar - Ivar8[i]) <= Ivar * 1e-6);
      }
    }
    if(!(fabs(R - R8[i]) <= R * 1e-6 + 1e-9)){
      #pragma omp critical
      {
	printf("\nFAIKD_mm256:i=%d:Y= %0.6f,resR2= %0.6e,resL2= %0.6e,newY= %0.6f:Ivar= %0.6e,Ivar8[i]= %0.6e, R= %0.7e, R8[i]= %0.7e\n",
	       i,Y[i], resR2, resL2[i], newY[i], Ivar, Ivar8[i], R, R8[i]);
	fflush(stdout);
	assert(fabs(R - R8[i]) <= R * 1e-6 + 1e-9);
      }
    }
  }
#endif

  return v_R;
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

static inline __m512 Pr_mm512(__m512 v_Y)
{
  float Y[16] __attribute__((aligned(64)));// HERE HERE : does this help or hinder the compiler ?
  float P[16] __attribute__((aligned(64)));

  _mm512_store_ps(Y, v_Y);

  #pragma omp simd aligned(Y,P : 64)
  for(int i = 0; i < 16; i++)
    P[i] = Pr(Y[i], resKBf, IresSDf);

  return _mm512_load_ps(P);
}

static inline __m512d Pr_mm512d(__m512d v_Y)
{
  double Y[8] __attribute__((aligned(64)));// HERE HERE : does this help or hinder the compiler ?
  double P[8] __attribute__((aligned(64)));

  _mm512_store_pd(Y, v_Y);

  #pragma omp simd aligned(Y,P : 64)
  for(int i = 0; i < 8; i++)
    P[i] = Pr(Y[i]);

  return _mm512_load_pd(P);
}

static inline __m512d log_mm512d(__m512d v_Y)
{
  double Y[8] __attribute__((aligned(64)));// HERE HERE : does this help or hinder the compiler ?
  double P[8] __attribute__((aligned(64)));

  _mm512_store_pd(Y, v_Y);

  #pragma omp simd aligned(Y,P : 64)
  for(int i = 0; i < 8; i++)
    P[i] = log(Y[i]);

  return _mm512_load_pd(P);
}

/* AVX512 version of FAIKD() */
static inline __m512 FAIKD_mm512(__m512 v_Y, /* new value of Yc(Y,I,K) - Yc(I-K-n,T) */
                                  MFLOAT resR2,/* VARtab[T*N + P] : only used with RES_VARIANCE */
				  __m512 v_resL2,/* VARtab[K*N + I] : only used with RES_VARIANCE */
				  __m512 v_newY,/* new assumed value of Y[I-K] - Y[I-K-n] to compute Pr() */
				  __m512 &v_Ivar)/* return value = 1.0/(F2var + S2var*Y + R2var*Y*Y + E2var * (resR2 + resL2)) */
  
{
  __m512 v_var = _mm512_add_ps(v512_F2var, _mm512_mul_ps(v512_S2var, v_Y));
#if VECTOR_DEBUG >= 2
  __m512 v_varA = v_var;
#endif

  if(QUADRATIC_VARIANCE){
    __m512 v_Ysq = _mm512_mul_ps(v_Y,v_Y);
    v_var = _mm512_add_ps(v_var, _mm512_mul_ps(v512_R2var, v_Ysq));
  }

#if VECTOR_DEBUG >= 2
  __m512 v_varB = v_var;
#endif
  if(RES_VARIANCE){
    __m512 v_resR2 = _mm512_set1_ps(resR2);
    v_var = _mm512_add_ps(v_var,_mm512_mul_ps(v512_E2var, _mm512_add_ps(v_resR2, v_resL2)));
  }

  v_Ivar = rsqrt_mm512(v_var);

  __m512 v_PR = Pr_mm512(v_newY);
  
  __m512 v_R = _mm512_mul_ps(v512_RefPen, v_PR);

#if VECTOR_DEBUG >= 2
  /* debug results using non-vector version of FAIKD() */
  float Y[16]; _mm512_storeu_ps(Y, v_Y);
  float resL2[16]; _mm512_storeu_ps(resL2, v_resL2);
  float newY[16]; _mm512_storeu_ps(newY, v_newY);
  float varV[16]; _mm512_storeu_ps(varV, v_var);
  float IvarV[16]; _mm512_storeu_ps(IvarV, v_Ivar);
  float RV[16]; _mm512_storeu_ps(RV, v_R);
  for(int i = 0; i < 16; i++){
    MFLOAT Ivar;
    MFLOAT R = FAIKD(Y[i], resR2, resL2[i], newY[i], Ivar);
    MFLOAT var = F2var + S2var * Y[i];
    if(QUADRATIC_VARIANCE)
      var += R2var * Y[i] * Y[i];
    if(RES_VARIANCE)
      var += E2var * (resR2 + resL2[i]);

    if(!(Ivar > 0.0 && IvarV[i] > 0.0 && fabs(Ivar - IvarV[i]) <= Ivar * 1e-6)){
      #pragma omp critical
      {
	printf("\nFAIKD_mm512:i=%d:Y= %0.6f,resR2= %0.6e,resL2= %0.6e,newY= %0.6f: Ivar= %0.10e(var= %0.10e), IvarV[i]= %0.10e(varV[i]= %0.10ee),R= %0.6e,RV[i]= %0.6e\n",
	       i,Y[i], resR2, resL2[i], newY[i], Ivar, var, IvarV[i], varV[i], R, RV[i]);
	float F2varV[16]; _mm512_storeu_ps(F2varV, v512_F2var);
	float S2varV[16]; _mm512_storeu_ps(S2varV, v512_S2var);
	float R2varV[16]; _mm512_storeu_ps(R2varV, v512_R2var);
	float E2varV[16]; _mm512_storeu_ps(E2varV, v512_E2var);
	printf("\t F2var=%0.7e(%0.7e),S2var=%0.7e(%0.7e),R2var=%0.7e(%0.7e),E2var=%0.7e(%0.7e)\n",F2var,F2varV[i],S2var,S2varV[i],R2var,R2varV[i],E2var,E2varV[i]);
	float varAV[16]; _mm512_storeu_ps(varAV, v_varA);
	float varBV[16]; _mm512_storeu_ps(varBV, v_varB);
	for(int j = 0; j < 16; j++)
	  printf("\t j=%d:Y[j]= %0.6f, F2var[j]= %0.7e, S2var[j]= %0.7e, R2var[j]= %0.7e: varA[j]= %0.7e, varB[j]= %0.7e, var[j]= %0.7e\n", 
		 j, Y[j], F2varV[j],S2varV[j],R2varV[j],varAV[j],varBV[j],varV[j]);
	fflush(stdout);
	assert(Ivar > 0.0 && IvarV[i] > 0.0 && fabs(Ivar - IvarV[i]) <= Ivar * 1e-6);
      }
    }
    if(!(fabs(R - RV[i]) <= R * 1e-6 + 1e-9)){
      #pragma omp critical
      {
	printf("\nFAIKD_mm512:i=%d:Y= %0.6f,resR2= %0.6e,resL2= %0.6e,newY= %0.6f:Ivar= %0.6e,IvarV[i]= %0.6e, R= %0.7e, RV[i]= %0.7e\n",
	       i,Y[i], resR2, resL2[i], newY[i], Ivar, IvarV[i], R, RV[i]);
	fflush(stdout);
	assert(fabs(R - RV[i]) <= R * 1e-6 + 1e-9);
      }
    }
  }
#endif

  return v_R;
}

#endif // VECTOR_AVX || VECTOR_MIC

/* modified FAIK to NOT use PRtab[].Pr (only used when H[] is temporarily changed and multiple threads are involved) */
/* NOTE : H[I-K-n-T .. I] includes a new added label and can be H[0] or H[N+1] */
static inline MFLOAT FAIKP(MFLOAT Y,
			   int n,
			   int I,
			   int K,/* right end of interval */
			   int T,/* left end of interval : only used with RES_VARIANCE */
			   MFLOAT &Ivar,/* return value = 1.0/(F2var + S2var*Y + R2var*Y*Y) */
			   double *H)/* Reference array */
  
{
  if(DEBUG>=2) assert(n >= 1 && n <= max(DELTA_Y,deltaExtYRef) + (DELTA_MIN ? DELTA_LIM : 0) +1/*for deleted site*/);
  if(DEBUG>=2) assert(0 <= K && K <= KMAX);
  if(DEBUG>=2 && !(I-K-n >= 0 && H[I-K-n] > 0.0 && I <= gN+1 && I <= PRNsiz)){
    printf("FAIKP(Y= %0.4f,n=%d,I=%d,K=%d,T=%d):PRNsiz=%d,gN=%d,H[I-K-n]= %0.4f\n",Y,n,I,K,T,PRNsiz,gN,H[max(0,min(PRNsiz,I-K-n))]);
    fflush(stdout);
    assert(I-K-n >= 0 && H[I-K-n] > 0.0 && I <= gN+1 && I <= PRNsiz);
  }
  if(DEBUG>=2) assert(Y>0.0);

  MFLOAT var = F2var + S2var * Y;
  if(QUADRATIC_VARIANCE)
     var += R2var * Y * Y;
  if(RES_VARIANCE){
    if(DEBUG>=2) assert(I-K-n-T >= 0);
    if(DEBUG>=2) assert(0 <= T && T <= KMAX);
    MFLOAT resR = H[I] - H[I-K];
    MFLOAT resL = H[I-K-n] - H[I-K-n-T];
    MFLOAT resvar = resL*resL + resR*resR;
    var += E2var * resvar;
  }

  if(DEBUG>=2) assert(var > 0.0);
  Ivar = sqrtF(((MFLOAT)1.0)/var);
  if(DEBUG>=2) assert(isfinite(Ivar));

  return /* FnPowf[n-1] *  */ RefPen * Pr(H[I-K]-H[I-K-n]);// NOTE : FnPow[n-1] has been moved in FAJ as expF(LogFnPow[n-1])
}

/* modified FAIKP, which does not index H[] */
static inline MFLOAT FAIKP(MFLOAT Y,
			   MFLOAT resR,/* Y[I] - Y[I-K] */
			   MFLOAT resL,/* Y[I-K-n] - Y[I-K-n-T] */
			   MFLOAT Yint,/* Y[I-K] - Y[I-K-n] */
			   int n,
			   int I,
			   int K,/* right end of interval */
			   int T,/* left end of interval : only used with RES_VARIANCE */
			   MFLOAT &Ivar,/* return value = 1.0/(F2var + S2var*Y + R2var*Y*Y) */
			   double *H)/* Reference array : only for debugging*/
  
{
  if(DEBUG>=2) assert(n >= 1 && n <= max(DELTA_Y,deltaExtYRef) + (DELTA_MIN ? DELTA_LIM : 0) +1/*for deleted site*/);
  if(DEBUG>=2) assert(0 <= K && K <= KMAX);
  if(DEBUG>=2 && !(I-K-n >= 0 && H[I-K-n] > 0.0 && I <= gN+1 && I <= PRNsiz)){
    printf("FAIKP(Y= %0.4f,n=%d,I=%d,K=%d,T=%d):PRNsiz=%d,gN=%d,H[I-K-n]= %0.4f\n",Y,n,I,K,T,PRNsiz,gN,H[max(0,min(PRNsiz,I-K-n))]);
    fflush(stdout);
    assert(I-K-n >= 0 && H[I-K-n] > 0.0 && I <= gN+1 && I <= PRNsiz);
  }
  if(DEBUG>=2) assert(Y > 0.0);
  if(DEBUG>=2) assert(isfinite(resR));
  if(DEBUG>=2) assert(isfinite(resL));
  if(DEBUG>=2) assert(isfinite(Yint));

  MFLOAT var = F2var + S2var * Y;
  if(QUADRATIC_VARIANCE)
     var += R2var * Y * Y;
  if(RES_VARIANCE){
    if(DEBUG>=2) assert(I-K-n-T >= 0);
    if(DEBUG>=2) assert(0 <= T && T <= KMAX);
    MFLOAT resvar = resL*resL + resR*resR;
    var += E2var * resvar;
  }

  if(DEBUG>=2) assert(var > 0.0);
  Ivar = sqrtF(((MFLOAT)1.0)/var);
  if(DEBUG>=2) assert(isfinite(Ivar));

  MFLOAT r = RefPen * Pr(Yint);
  if(DEBUG>=2) assert(isfinite(r));
  return r;
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT && DEBUG < 2

// AVX256 version of FAIKP (without indexing Y[]) 
static inline __m256 FAIKP_mm256(__m256 Y,
                                 __m256 resR,/* Y[I] - Y[I-K] */
				 __m256 resL,/* Y[I-K-n] - Y[I-K-n-T] */
				 __m256 Yint,/* Y[I-K] - Y[I-K-n] */
				 __m256 &Ivar)/* return value = 1.0/sqrt(F2var + S2var*Y + R2var*Y*Y) */
{
  __m256 var = v256_F2var + v256_S2var * Y;

  if(QUADRATIC_VARIANCE)
     var += v256_R2var * Y * Y;
  
  if(RES_VARIANCE){
    __m256 resvar = resL*resL + resR*resR;
    var += v256_E2var * resvar;
  }

  /* Fast Inverse Square Root with One Newton Rapson iteration */
  __m256 Ivar1 = _mm256_rsqrt_ps(var);
  Ivar = Ivar1 * (v256_1p5 - v256_0p5 * var * Ivar1 *Ivar1);

  return v256_RefPen * Pr_mm256(Yint);
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

// AVX512 version of FAIKP (without indexing Y[]) 
static inline __m512 FAIKP_mm512(__m512 Y,
                                 __m512 resR,/* Y[I] - Y[I-K] */
				 __m512 resL,/* Y[I-K-n] - Y[I-K-n-T] */
				 __m512 Yint,/* Y[I-K] - Y[I-K-n] */
				 __m512 &Ivar)/* return value = 1.0/sqrt(F2var + S2var*Y + R2var*Y*Y) */
{
  __m512 var = _mm512_add_ps(v512_F2var, _mm512_mul_ps(v512_S2var, Y));

  if(QUADRATIC_VARIANCE)
    var = _mm512_add_ps(var, _mm512_mul_ps(v512_R2var, _mm512_mul_ps(Y,  Y)));
  
  if(RES_VARIANCE){
    __m512 resvar = _mm512_add_ps(_mm512_mul_ps(resL,resL) , _mm512_mul_ps(resR,resR));
    var = _mm512_add_ps(var, _mm512_mul_ps(v512_E2var, resvar));
  }

  Ivar = rsqrt_mm512(var);

  return _mm512_mul_ps(v512_RefPen , Pr_mm512(Yint));
}

#endif // VECTOR_MIC

/* compute FA() based on previously computed FAIK() value */
static inline MFLOAT FAJ(MFLOAT X,
			 MFLOAT Y,
			 int m,/* debug only */
			 int n,/* debug only */
			 MFLOAT FBiasWTxm,
			 MFLOAT LogFpPowm1,
			 MFLOAT LogFnPown1,
			 MFLOAT PrBiasWT,/* PrBiasWT[I-K] - PrBiasWT[I-K-n] */ // NOTE : not used, see note in RGentigScore.h
			 MFLOAT FAik,/* previous return value of FAIK() */
			 MFLOAT Ivar, /* sqrt(1.0/(F2var+S2var*Y + ...)) computed in FAIK() */
			 const int outlier, /* != 0 : check for outliers (compile time constant) */  
			 const int OUT_TYPE, /* value of global OUTLIER_TYPE */
			 const int OUT_OLD   /* value of global OUTLIER_OLD(outlierMaxBC) */
			 /* , int verb = 0 */)
{
  if(DEBUG>=2) assert(m <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX + 1);
  if(DEBUG>=2) assert(isfinite(Ivar));

  MFLOAT err = (X-Y) * Ivar;
  if(DEBUG>=2 && !(PrBiasWT >= 0.0)){
    printf("FAJ:X=%0.6f,Y=%0.6f,m=%d,n=%d:PrBiasWT= %0.6f\n",X,Y,m,n,PrBiasWT);
    fflush(stdout);
    assert(PrBiasWT >= 0.0);
  }
  MFLOAT Bias = X * FrateBiasWT + /* PrBiasWT + */ FBiasWTxm;

  MFLOAT OutPen = -fabsF(X-Y) * OutlierLambdaInv;
  if(DEBUG>=2 && biasWT <= 0.0) assert(Bias <= 0.0);

  MFLOAT err2 = err*err;
  MFLOAT logerrA = Bias + LogFpPowm1; // WAS8 - X*Frate;
  if(!FRATE_FIX)
    logerrA -= X*Frate;
  MFLOAT logerr = logerrA - err2;
  if(FRATE_FIX)
    logerr -= X*Frate;

  // NOTE MaxFA cannot be applied to Fn ^^ (n-1), since mprobeval assumes FAJ changes by exactly a multiple of Fn if n changes by 1 
  logerr = LogFnPown1 + min(logerr, logMaxFAG);// avoid overflow for large intervals, especially with -biasWT 0

  MFLOAT Vexp = expF(logerr);
  if((DEBUG>=2 && !isfinite(Vexp)) /* || verb */){
    printf("FAJ:X=%0.8f,Y=%0.8f,m=%d,n=%d,outlier=%d,OUT_TYPE=%d,OUT_OLD=%d: LogFnPown1= %0.8e, LogFpPowm1= %0.8e, PrBiasWT= %0.8e\n",X,Y,m,n,outlier,OUT_TYPE,OUT_OLD, LogFnPown1, LogFpPowm1, PrBiasWT);
    printf("\t logerrA= %0.8e, logerr= %0.8e, logMaxFAG= %0.8e, Vexp= %0.8e, MaxFloat= %0.8e\n",logerrA, logerr, logMaxFAG, Vexp, MaxFloat);
    fflush(stdout);
    assert(isfinite(Vexp));
  }
  if(DEBUG>=2 && !(FAik * FN[0] <= RefPen * 1.00001)){
    printf("\nFAJ: FAik= %0.8e, RefPen= %0.8e\n",FAik, RefPen);
    fflush(stdout);
    assert(FAik * FN[0] <= RefPen * 1.000001);
  }
  MFLOAT LR = FAik * Ivar * Vexp;
  if(DEBUG>=2) assert(isfinite(LR) && LR <= MaxFA);
  if((DEBUG>=2 && rverb) /* || verb */){
    printf("FAJ:X=%0.8f,Y=%0.8f,m=%d,n=%d,outlier=%d,OUT_TYPE=%d,OUT_OLD=%d:\n",X,Y,m,n,outlier,OUT_TYPE,OUT_OLD);
    printf("\t Ivar=%0.8e,FAik=%0.8e(RefPen=%0.8e),Bias=%0.8e,OutPen=%0.8e,err=%0.8e,err2=%0.8e,logerr=%0.8e,Vexp=%0.8e,LR=%0.8e,Frate=%0.8f\n",
	   Ivar,FAik,RefPen,Bias,OutPen,err,err2,logerr,Vexp,LR,Frate);
    fflush(stdout);
  }

  if(outlier && OUTLIER_DELTA(X-Y)){
    if(!OUT_TYPE){
      MFLOAT LRout = expF(OutPen + Bias * biasWToutlierF);
      if(DEBUG>=2) assert(isfinite(LRout) && LRout <= MaxFA);
      if(OUT_OLD)
	return LRout * PoutlierXthetaBYlambda + LR * PoutlierM1;
      else
	return LRout * PoutlierXthetaBYmax + LR * PoutlierM1;
    } else {/* OUTLIER_TYPE==1 */
      MFLOAT logerr2 = logerrA  /* NEW */+ OutPen - Bias * biasWToutlierM1;
      logerr2 = LogFnPown1 + min(logerr2, logMaxFAG2);// avoid overflow for large intervals
      MFLOAT Vexp2 = expF(logerr2);
      if(DEBUG>=2) assert(isfinite(Vexp2) && Vexp2 <= MaxFAG2);

      MFLOAT LRout = FAik * Vexp2;
      if((DEBUG>=2 && (rverb || !(isfinite(LRout) && LRout <= MaxFA))) /* || verb */){
	printf("\t loggerr2=%0.8e,Vexp2=%0.8e,LRout=%0.8e(FAik=%0.8e,RefPen=%0.8e),LR= %0.8e, PoutlierPIBYmax=%0.8e,PoutlierM1=%0.8e, MaxFAG2=%0.8e,MaxFA=%0.8e\n",
	       logerr2, Vexp2, LRout, FAik,RefPen,LR, PoutlierPIBYmax,PoutlierM1,MaxFAG2,MaxFA);
	fflush(stdout);
      }
      if(DEBUG>=2) assert(isfinite(LRout) && LRout <= MaxFA);

      if(OUT_OLD)
	return LRout * PoutlierPIbXtheta + LR * PoutlierM1;
      else
	return LRout * PoutlierPIBYmax + LR * PoutlierM1;
    }
  } else
    return LR;
}

// Vector friendly version without debug args */

#pragma omp declare simd
template<int outlier, int OUT_TYPE, int OUT_OLD>
static inline MFLOAT FAJ(MFLOAT X,
			 MFLOAT Y,
			 MFLOAT FBiasWTxm,
			 MFLOAT LogFpPowm1,
			 MFLOAT LogFnPown1,
			 MFLOAT PrBiasWT,/* PrBiasWT[I-K] - PrBiasWT[I-K-n] */// NOTE : not used, see note in RGentigScore.h
			 MFLOAT FAik,/* previous return value of FAIK() */
			 MFLOAT Ivar) /* sqrt(1.0/(F2var+S2var*Y + ...)) computed in FAIK() */
{
  if(DEBUG>=2) assert(isfinite(Ivar));

  MFLOAT err = (X-Y) * Ivar;
  if(DEBUG>=2 && !(PrBiasWT >= 0.0)){
    printf("FAJ: Y=%0.4f, PrBiasWT= %0.8e\n",Y,PrBiasWT);
    fflush(stdout);
    assert(PrBiasWT >= 0.0);
  }
  MFLOAT Bias = X * FrateBiasWT /* + PrBiasWT */ + FBiasWTxm;
  MFLOAT OutPen = -fabsF(X-Y) * OutlierLambdaInv;
  if(DEBUG>=2 && biasWT <= 0.0) assert(Bias <= 0.0);

  MFLOAT err2 = err*err;
  MFLOAT logerrA =  Bias + LogFpPowm1;// WAS8 - X*Frate;
  if(!FRATE_FIX)
    logerrA -= X*Frate;
  MFLOAT logerr = logerrA - err2;
  if(FRATE_FIX)
    logerr -= X*Frate;

  // NOTE MaxFA cannot be applied to Fn ^^ (n-1), since mprobeval assumes FAJ changes by exactly a multiple of Fn if n changes by 1 
  logerr = LogFnPown1 + min(logerr, logMaxFAG);// avoid overflow for large intervals, especially with -biasWT 0

  MFLOAT Vexp = expF(logerr);
  if(DEBUG>=2) assert(isfinite(Vexp) && Vexp <= MaxFAG);
  if(DEBUG>=2) assert(FAik * FN[0] <= RefPen * 1.000001);
  MFLOAT LR = FAik * Ivar * Vexp;
  if(DEBUG>=2 && !(isfinite(LR) && LR <= MaxFA)){
    printf("FAJ(X= %0.6f, Y= %0.6f, FBiasWTxm= %0.8e,LogFpPowm1=%0.8e,LogFnPown1=%0.8e,FAik=%0.8e,Ivar=%0.8e):\n",
	   X,Y,FBiasWTxm,LogFpPowm1,LogFnPown1,FAik,Ivar);
    printf("\t Bias= %0.8e, err= %0.8e, OutPen= %0.8e, err2= %0.8e, logerrA= %0.8e, logerr= %0.8e, Vexp= %0.8e: LR= %0.8e\n",
	   Bias, err, OutPen, err2, logerrA, logerr, Vexp, LR);
    printf("\t logerr= %0.8e, exp(logerr)= %0.8e, expf(logerr)= %0.8e\n", logerr, exp(logerr), expf(logerr));
    fflush(stdout);
    assert(isfinite(LR) && LR <= MaxFA);
  }

  if(outlier && OUTLIER_DELTA(X-Y)){
    if(!OUT_TYPE){
      MFLOAT LRout = expF(OutPen + Bias * biasWToutlierF);
      if(DEBUG>=2) assert(isfinite(LRout) && LRout <= MaxFA);
      if(OUT_OLD)
	return LRout * PoutlierXthetaBYlambda + LR * PoutlierM1;
      else
	return LRout * PoutlierXthetaBYmax + LR * PoutlierM1;
    } else {/* OUTLIER_TYPE==1 */
      MFLOAT logerr2 = logerrA  /* NEW */+ OutPen - Bias * biasWToutlierM1;
      logerr2 = LogFnPown1 + min(logerr2, logMaxFAG2);// avoid overflow for large intervals
      MFLOAT Vexp2 = expF(logerr2);
      if(DEBUG>=2) assert(isfinite(Vexp2) && Vexp2 <= MaxFAG2);

      MFLOAT LRout = FAik * Vexp2;
      if(DEBUG>=2) assert(isfinite(LRout) && LRout <= MaxFA);
      if(OUT_OLD)
	return LRout * PoutlierPIbXtheta + LR * PoutlierM1;
      else
	return LRout * PoutlierPIBYmax + LR * PoutlierM1;
    }
  } else
    return LR;
}

#if 0
/* same as FAJ<>(), but a slightly different SIMD declaration */

#pragma omp declare simd uniform(X, FBiasWTxm, LogFpPowm1, LogFnPown1)
template<int outlier, int OUT_TYPE, int OUT_OLD>
  static inline MFLOAT FAJd(MFLOAT X,
			    MFLOAT Y,
			    MFLOAT FBiasWTxm,
			    MFLOAT LogFpPowm1,
			    MFLOAT LogFnPown1,
			    MFLOAT PrBiasWT,/* PrBiasWT[I-K] - PrBiasWT[I-K-n] */ // NOTE : not used
			    MFLOAT FAik,/* previous return value of FAIK() */
			    MFLOAT Ivar) /* sqrt(1.0/(F2var+S2var*Y + ...)) computed in FAIK() */
{
  if(DEBUG>=2) assert(isfinite(Ivar));

  MFLOAT err = (X-Y) * Ivar;
  //  if(DEBUG>=2) assert(PrBiasWT >= 0.0);
  MFLOAT Bias = X * FrateBiasWT /* + PrBiasWT */ + FBiasWTxm;
  MFLOAT OutPen = -fabsF(X-Y) * OutlierLambdaInv;
  if(DEBUG>=2 && biasWT <= 0.0) assert(Bias <= 0.0);

  MFLOAT err2 = err*err;
  MFLOAT logerrA = Bias + LogFpPowm1;// WAS8 - X*Frate;
  if(!FRATE_FIX)
    logerrA -= X*Frate;
  MFLOAT logerr = logerrA - err2;
  id(FRATE_FIX)
    logerr -= X*Frate;

  // NOTE MaxFA cannot be applied to Fn ^^ (n-1), since mprobeval assumes FAJ changes by exactly a multiple of Fn if n changes by 1 
  logerr = LogFnPown1 + min(logerr, logMaxFAG);// avoid overflow for large intervals, especially with -biasWT 0

  MFLOAT Vexp = expF(logerr);
  if(DEBUG>=2) assert(isfinite(Vexp) && Vexp <= MaxFAG);
  if(DEBUG>=2) assert(FAik * FN[0] <= RefPen * 1.000001);
  MFLOAT LR = FAik * Ivar * Vexp;
  if(DEBUG>=2) assert(isfinite(LR) && LR <= MaxFA);

  if(outlier && OUTLIER_DELTA(X-Y)){
    if(!OUT_TYPE){
      MFLOAT LRout = expF(OutPen + Bias * biasWToutlierF);
      if(DEBUG>=2) assert(isfinite(LRout) && LRout <= MaxFA);
      if(OUT_OLD)
	return LRout * PoutlierXthetaBYlambda + LR * PoutlierM1;
      else
	return LRout * PoutlierXthetaBYmax + LR * PoutlierM1;
    } else {/* OUTLIER_TYPE==1 */
      MFLOAT logerr2 = logerrA /* NEW */+ OutPen - Bias * biasWToutlierM1;
      logerr2 = LogFnPown1 + min(logerr2, logMaxFAG2);// avoid overflow for large intervals
      MFLOAT Vexp2 = expF(logerr2);
      if(DEBUG>=2) assert(isfinite(Vexp2) && Vexp2 <= MaxFAG2);
      MFLOAT LRout = FAik * Vexp2;
      if(DEBUG>=2) assert(isfinite(LRout) && LRout <= MaxFA);
      if(OUT_OLD)
	return LRout * PoutlierPIbXtheta + LR * PoutlierM1;
      else
	return LRout * PoutlierPIBYmax + LR * PoutlierM1;
    }
  } else
    return LR;
}
#endif // FAJd

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT

/* AVX256 version of FAJ<>() */
template<int outlier, int OUT_TYPE, int OUT_OLD>
static inline __m256 FAJ_mm256(float X,
			       __m256 v_Y,
			       int m, // only used with DEBUG>=2 or VECTOR_DEBUG
			       int n, // only used with DEBUG>=2 or VECTOR_DEBUG
			       float FBiasWTxm,
			       float LogFpPowm1,
			       __m256 v_LogFnPown1,
			       __m256 v_PrBiasWT, // NOTE : not used
			       __m256 v_FAik,/* previous return value of FAIK() */
			       __m256 v_Ivar) /* sqrt(1.0/(F2var+S2var*Y + ...)) computed in FAIK() */
{
  __m256 v_X = _mm256_set1_ps(X);
  __m256 v_Bias = v_X * v256_FrateBiasWT /* + v_PrBiasWT */ + _mm256_set1_ps(FBiasWTxm);
  __m256 v_LogFpPowm1 = _mm256_set1_ps(LogFpPowm1);

  __m256 v_err = v_Ivar * (v_X - v_Y);// Ivar * (X - Y)
  __m256 v_fabsXmY = _mm256_min_ps(v_X , v_Y) - _mm256_max_ps(v_X, v_Y);// -fabsF(X-Y) = min(X,Y) - max(X,Y)
  __m256 v_OutPen = v256_OutlierLambdaInv * v_fabsXmY;
  __m256 v_err2 = v_err * v_err;

  __m256 v_logerrA = v_Bias + v_LogFpPowm1; // WAS8 - v_X * v256_Frate;
  if(!FRATE_FIX)
    v_logerrA -= v_X * v256_Frate;
  __m256 v_logerr = v_logerrA - v_err2;
  if(FRATE_FIX)
    v_logerr -= v_X * v256_Frate;

  v_logerr = v_LogFnPown1 + _mm256_min_ps(v_logerr, v256_logMaxFAG);// avoid overflow for large intervals, especially with -biasWT 0

  __m256 v_Vexp = expF_mm256(v_logerr);
  __m256 v_LR = v_FAik * v_Ivar * v_Vexp;

  if(!outlier)
    return v_LR;

#if (VECTOR_DEBUG >= 2 || VECTOR2_DEBUG >= 2)
  __m256 v_logerr1;
#endif

  __m256 v_LR2, v_logerr2, v_LRout, v_Vexp2;
  if(!OUT_TYPE){
    v_logerr2 = v_OutPen + v_Bias * v256_biasWToutlierF;
    v_LRout = expF_mm256(v_logerr2);

    if(OUT_OLD)
      v_LR2 = v_LRout * v256_PoutlierXthetaBYlambda + v_LR * v256_PoutlierM1;
    else
      v_LR2 = v_LRout * v256_PoutlierXthetaBYmax + v_LR * v256_PoutlierM1;
  } else {/* OUTLIER_TYPE==1 */
    v_logerr2 = v_logerrA  /* NEW */+ v_OutPen - v_Bias * v256_biasWToutlierM1;

#if (VECTOR_DEBUG >= 2 || VECTOR2_DEBUG >= 2)
    v_logerr1 = v_logerr2;
#endif
    v_logerr2 = v_LogFnPown1 + _mm256_min_ps(v_logerr2, v256_logMaxFAG2);// avoid overflow for large intervals

    v_Vexp2 = expF_mm256(v_logerr2);

    v_LRout = v_FAik * v_Vexp2;
    if(OUT_OLD)
      v_LR2 = v_LRout * v256_PoutlierPIbXtheta + v_LR * v256_PoutlierM1;
    else
      v_LR2 = v_LRout * v256_PoutlierPIBYmax + v_LR * v256_PoutlierM1;
  }

  __m256 v_R = _mm256_blendv_ps(v_LR, v_LR2, 
				_mm256_cmp_ps(v_fabsXmY, v256_outlierNegMaxBC, _CMP_GT_OQ));// (-fabs(X-Y) < -outlierMaxBC) ? LR2 : LR;

#if (VECTOR_DEBUG >= 2 || VECTOR2_DEBUG >= 2)
  float Y[8]; _mm256_storeu_ps(Y, v_Y);
  float LogFnPown1[8]; _mm256_storeu_ps(LogFnPown1, v_LogFnPown1);
  float PrBiasW[8]; _mm256_storeu_ps(PrBiasW, v_PrBiasWT);
  float FAik[8]; _mm256_storeu_ps(FAik, v_FAik);
  float Ivar[8]; _mm256_storeu_ps(Ivar, v_Ivar);
  float R8[8]; _mm256_storeu_ps(R8, v_R);

  for(int i = 0; i < 8; i++){
    MFLOAT R = FAJ<outlier, OUT_TYPE, OUT_OLD>(X, Y[i], FBiasWTxm, LogFpPowm1, LogFnPown1[i], PrBiasW[i], FAik[i], Ivar[i]);
    if(!(R >= 0.0 && R8[i] >= 0.0 && fabs(R - R8[i]) < R * 1e-5 + 1e-35)){
      #pragma omp critical
      {
	printf("\nFAJ_mm256:i=%d:m=%d,n=%d,X=%0.6f,Y=%0.6f,FBiasWTxm=%0.6e,LogFpPowm1=%0.6e,LogFnPown1=%0.6e,PrBiasWT=%0.6e,FAik=%0.6e,Ivar=%0.6e:R8[i]= %0.7e, R= %0.7e\n",
	       i,m,n,X,Y[i],FBiasWTxm,LogFpPowm1,LogFnPown1[i],PrBiasW[i],FAik[i],Ivar[i],R8[i],R);

	MFLOAT err = (X - Y[i]) * Ivar[i];
	float errV[8]; _mm256_storeu_ps(errV, v_err);

	MFLOAT Bias = X * FrateBiasWT /* + PrBiasW[i] */ + FBiasWTxm;
	float BiasV[8]; _mm256_storeu_ps(BiasV, v_Bias);
	
	MFLOAT OutPen = -fabsF(X - Y[i]) * OutlierLambdaInv;
	float OutPenV[8]; _mm256_storeu_ps(OutPenV, v_OutPen);

	MFLOAT err2 = err*err;
	float err2V[8]; _mm256_storeu_ps(err2V, v_err2);

	printf("\t err= %0.7e(errV[i]= %0.7e), Bias= %0.7e(BiasV[i]=%0.7e), OutPen= %0.7e(OutPenV[i]= %0.7e), err2= %0.7e(err2V[i]= %0.7e)\n", 
	       err, errV[i], Bias, BiasV[i], OutPen, OutPenV[i], err2, err2V[i]);
	
	float OutlierLambdaInvV[8]; _mm256_storeu_ps(OutlierLambdaInvV, v256_OutlierLambdaInv);
	float fabsXmYV[8]; _mm256_storeu_ps(fabsXmYV, v_fabsXmY);
	/*	printf("\t OutlierLambdaInv= %0.6e (OutlierLambdaInvV[i]= %0.6e), -fabsF(X-Y)= %0.6e(fabsXmYV[i]= %0.6e)\n",
		OutlierLambdaInv, OutlierLambdaInvV[i], -fabsF(X-Y[i]), fabsXmYV[i]);*/

	MFLOAT logerrA = Bias + LogFpPowm1; // WAS8 - X*Frate ;
	if(!FRATE_FIX)
	  logerrA -= X*Frate;
	MFLOAT logerr = logerrA - err2;
	if(FRATE_FIX)
	  logerr -= X*Frate;

	logerr = LogFnPown1[i] + min(logerr, logMaxFAG);// avoid overflow for large intervals, especially with -biasWT 0
	float logerrAV[8]; _mm256_storeu_ps(logerrAV, v_logerrA);
	float logerrV[8]; _mm256_storeu_ps(logerrV, v_logerr);

	MFLOAT Vexp = expF(logerr);
	float VexpV[8]; _mm256_storeu_ps(VexpV, v_Vexp);

	MFLOAT LR = FAik[i] * Ivar[i] * Vexp;
	float LRV[8]; _mm256_storeu_ps(LRV, v_LR);

	printf("\t logerrA= %0.7e(logerrAV[i]= %0.7e), logerr= %0.7e(logerrV[i]= %0.7e), Vexp= %0.7e(VexpV[i]= %0.7e), LR= %0.7e(LRV[i]= %0.7e)\n", 
	       logerrA, logerrAV[i], logerr, logerrV[i], Vexp, VexpV[i], LR, LRV[i]);
	printf("\t\t X= %0.7f, Frate= %0.7e, LogFpPowm1[i]= %0.7e, logFnPown1= %0.7e\n", X, Frate, LogFpPowm1, LogFnPown1[i]);

	MFLOAT LR2;
	if(!OUT_TYPE){
	  MFLOAT logerr2 = OutPen + Bias  * biasWToutlierF;
	  float logerr2V[8]; _mm256_storeu_ps(logerr2V, v_logerr2);

	  MFLOAT LRout = expF(logerr2);
	  float LRoutV[8]; _mm256_storeu_ps(LRoutV, v_LRout);

	  if(OUT_OLD)
	    LR2 = LRout * PoutlierXthetaBYlambda + LR * PoutlierM1;
	  else
	    LR2 = LRout * PoutlierXthetaBYmax + LR * PoutlierM1;
	  float LR2V[8]; _mm256_storeu_ps(LR2V, v_LR2);

	  printf("\t logerr2= %0.7e(logerr2V[i]= %0.7e), LRout= %0.7e(LRoutV[i]= %0.7e), LR=%0.7e,LR2= %0.7e(LR2V[i]= %0.7e),OUT_OLD=%d\n",
		 logerr2,logerr2V[i],LRout, LRoutV[i], LR, LR2, LR2V[i], OUT_OLD);

	} else {/* OUTLIER_TYPE==1 */
	  MFLOAT logerr1 = logerrA - Bias * biasWToutlierM1 /* NEW */ + OutPen;
	  float logerr1V[8]; _mm256_storeu_ps(logerr1V, v_logerr1);
	  float logMaxFAG2V[8]; _mm256_storeu_ps(logMaxFAG2V, v256_logMaxFAG2);
	  float biasWToutlierM1V[8]; _mm256_storeu_ps(biasWToutlierM1V, v256_biasWToutlierM1);
	  float OutPenV[8]; _mm256_storeu_ps(OutPenV, v_OutPen);

	  printf("\t logerr1= %0.7e(logerr1V[i]= %0.7e), logMaxFAG2= %0.7e(logMaxFAG2V[i]= %0.7e), biasWToutlierM1= %0.7e, (biasWToutlierM1V[i]= %0.7e), OutPen= %0.7e(OutPenV[i]=%0.7e)\n",
		 logerr1,logerr1V[i], logMaxFAG2, logMaxFAG2V[i], biasWToutlierM1, biasWToutlierM1V[i], OutPen, OutPenV[i]);

	  MFLOAT logerr2 = LogFnPown1[i] + min(logerr1, logMaxFAG2);
	  float logerr2V[8]; _mm256_storeu_ps(logerr2V, v_logerr2);

	  MFLOAT Vexp2 = expF(logerr2);
	  float Vexp2V[8]; _mm256_storeu_ps(Vexp2V, v_Vexp2);

	  MFLOAT LRout = FAik[i] * Vexp2;
	  float LRoutV[8]; _mm256_storeu_ps(LRoutV, v_LRout);

	  if(OUT_OLD)
	    LR2 = LRout * PoutlierPIbXtheta + LR * PoutlierM1;
	  else
	    LR2 = LRout * PoutlierPIBYmax + LR * PoutlierM1;
	  float LR2V[8]; _mm256_storeu_ps(LR2V, v_LR2);

	  printf("\t logerr2= %0.7e(logerr2V[i]= %0.7e), Vexp2= %0.7f(Vexp2V[i]= %0.7e), LRout= %0.7e(LRoutV[i]= %0.7e), LR2= %0.7e(LR2V[i]= %0.7e),OUT_OLD=%d\n",
		 logerr2, logerr2V[i], Vexp2, Vexp2V[i], LRout, LRoutV[i], LR2, LR2V[i],OUT_OLD);
	}
	
	float outlierNegMaxBCV[8]; _mm256_storeu_ps(outlierNegMaxBCV, v256_outlierNegMaxBC);
	printf("\t OUTLIER_DELTA(X-Y[i])= %d,fabsXmYV[i]= %0.6f, outlierMaxBC= %0.6e(outlierNegMaxBCV[i]= %0.6e)\n", 
	       OUTLIER_DELTA(X-Y[i]) ? 1 : 0, fabsXmYV[i], outlierMaxBC, outlierNegMaxBCV[i]);

	fflush(stdout);
	assert(R >= 0.0 && R8[i] >= 0.0 && fabs(R - R8[i]) < R * 1e-5 + 1e-35);
      }
    }
  }

#endif

  return v_R;
}

/* FAJ_mm256 for use in qprobeval() : slightly different vector args */
template<int outlier, int OUT_TYPE, int OUT_OLD>
static inline __m256 FAJ_mm256(__m256 v_X,
			       __m256 v_Y,
			       __m256 v_FBiasWTxm,
			       __m256 v_LogFpPowm1,
			       __m256 v_LogFnPown1,
			       __m256 v_PrBiasWT,// NOTE : not used
			       __m256 v_FAik,/* previous return value of FAIK() */
			       __m256 v_Ivar) /* sqrt(1.0/(F2var+S2var*Y + ...)) computed in FAIK() */
{
  __m256 v_Bias = v_X * v256_FrateBiasWT /* + v_PrBiasWT */ + v_FBiasWTxm;
  __m256 v_err = v_Ivar * (v_X - v_Y);// Ivar * (X - Y)
  __m256 v_fabsXmY = _mm256_min_ps(v_X , v_Y) - _mm256_max_ps(v_X, v_Y);// -fabsF(X-Y) = min(X,Y) - max(X,Y)
  __m256 v_OutPen = v256_OutlierLambdaInv * v_fabsXmY;

  __m256 v_err2 = v_err * v_err;
  __m256 v_logerrA = v_Bias + v_LogFpPowm1; // WAS8 - v_X * v256_Frate;
  if(!FRATE_FIX)
    v_logerrA -= v_X * v256_Frate;
  __m256 v_logerr = v_logerrA - v_err2;
  if(FRATE_FIX)
    v_logerr -= v_X * v256_Frate;

  v_logerr =  v_LogFnPown1 + _mm256_min_ps(v_logerr, v256_logMaxFAG);// avoid overflow for large intervals, especially with -biasWT 0

  __m256 v_Vexp = expF_mm256(v_logerr);
  __m256 v_LR = v_FAik * v_Ivar * v_Vexp;

  if(!outlier)
    return v_LR;

  __m256 v_LR2, v_logerr2, v_LRout, v_Vexp2;
  if(!OUT_TYPE){
    v_logerr2 = v_OutPen + v_Bias * v256_biasWToutlierF;
    v_LRout = expF_mm256(v_logerr2);

    if(OUT_OLD)
      v_LR2 = v_LRout * v256_PoutlierXthetaBYlambda + v_LR * v256_PoutlierM1;
    else
      v_LR2 = v_LRout * v256_PoutlierXthetaBYmax + v_LR * v256_PoutlierM1;
  } else {/* OUTLIER_TYPE==1 */
    v_logerr2 = v_logerrA  /* NEW */ + v_OutPen - v_Bias * v256_biasWToutlierM1;
    v_logerr2 = v_LogFnPown1 + _mm256_min_ps(v_logerr2, v256_logMaxFAG2);// avoid overflow for large intervals
    v_Vexp2 = expF_mm256(v_logerr2);

    v_LRout = v_FAik * v_Vexp2;
    if(OUT_OLD)
      v_LR2 = v_LRout * v256_PoutlierPIbXtheta + v_LR * v256_PoutlierM1;
    else
      v_LR2 = v_LRout * v256_PoutlierPIBYmax + v_LR * v256_PoutlierM1;
  }

  return _mm256_blendv_ps(v_LR, v_LR2, 
			  _mm256_cmp_ps(v_fabsXmY, v256_outlierNegMaxBC, _CMP_GT_OQ));// (-fabs(X-Y) < -outlierMaxBC) ? LR2 : LR;
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT


/* AVX512 version of FAJ<>() */
template<int outlier, int OUT_TYPE, int OUT_OLD>
static inline __m512 FAJ_mm512(float X,
			       __m512 v_Y,
			       int m, // only used with DEBUG>=2 or VECTOR_DEBUG
			       int n, // only used with DEBUG>=2 or VECTOR_DEBUG
			       float FBiasWTxm,
			       float LogFpPowm1,
			       __m512 v_LogFnPown1,
			       __m512 v_PrBiasWT,// NOTE : not used
			       __m512 v_FAik,/* previous return value of FAIK() */
			       __m512 v_Ivar) /* sqrt(1.0/(F2var+S2var*Y + ...)) computed in FAIK() */
{
  __m512 v_X = _mm512_set1_ps(X);
  //  __m512 v_Bias = _mm512_add_ps(_mm512_mul_ps(v_X, v512_FrateBiasWT), _mm512_add_ps(v_PrBiasWT,_mm512_set1_ps(FBiasWTxm)));
  __m512 v_Bias = _mm512_add_ps(_mm512_mul_ps(v_X, v512_FrateBiasWT), _mm512_set1_ps(FBiasWTxm));
  __m512 v_LogFpPowm1 = _mm512_set1_ps(LogFpPowm1);

#if M512_ABS==0
  __m512 v_err = _mm512_mul_ps(v_Ivar, _mm512_sub_ps(v_X, v_Y));// Ivar * (X - Y)
  __m512 v_fabsXmY = _mm512_sub_ps( _mm512_min_ps(v_X, v_Y), _mm512_max_ps(v_X, v_Y));// -fabsF(X-Y)
#else // NEW : use mm512_abs_ps
  __m512 v_XmY = _mm512_sub_ps(v_X, v_Y);
  __m512 v_err = _mm512_mul_ps(v_Ivar, v_XmY);
  __m512 v_fabsXmY = _mm512_abs_ps(v_XmY); // fabsF(X-Y)
#endif
  __m512 v_OutPen = _mm512_mul_ps(v512_OutlierLambdaInv, v_fabsXmY);
  __m512 v_err2 = _mm512_mul_ps(v_err, v_err);

  __m512 v_logerrA, v_logerr;
  if(!FRATE_FIX){
    v_logerrA = _mm512_sub_ps(_mm512_add_ps(v_Bias, v_LogFpPowm1), _mm512_mul_ps(v_X, v512_Frate));// logerrA = Bias + LogFpPowm1 - X*Frate
    v_logerr = _mm512_sub_ps(v_logerrA, v_err2); // logerr =  logerrA - err2;
  } else { // FRATE_FIX
    v_logerrA = _mm512_add_ps(v_Bias, v_LogFpPowm1);// logerrA = Bias + logFpPowm1
    v_logerr = _mm512_sub_ps(_mm512_sub_ps(v_logerrA, v_err2), _mm512_mul_ps(v_X, v512_Frate)); // logerr =  logerrA - err2 - X*Frate;
  }
  v_logerr = _mm512_add_ps(v_LogFnPown1, _mm512_min_ps(v_logerr, v512_logMaxFAG));// avoid overflow for large intervals, especially with -biasWT 0

  __m512 v_Vexp = expF_mm512(v_logerr);
  __m512 v_LR = _mm512_mul_ps(_mm512_mul_ps(v_FAik, v_Ivar), v_Vexp);// result without outlier */

  if(!outlier)
    return v_LR;

  __m512 v_LR2, v_logerr2, v_LRout, v_Vexp2;
  if(!OUT_TYPE){
#if M512_ABS==0
    v_logerr2 = _mm512_add_ps(v_OutPen, _mm512_mul_ps(v_Bias, v512_biasWToutlierF));
#else
    v_logerr2 = _mm512_sub_ps(_mm512_mul_ps(v_Bias, v512_biasWToutlierF), v_OutPen);
#endif
    v_LRout = expF_mm512(v_logerr2);
    if(OUT_OLD)
      v_LR2 = _mm512_add_ps(_mm512_mul_ps(v_LRout, v512_PoutlierXthetaBYlambda), _mm512_mul_ps(v_LR, v512_PoutlierM1));
    else
      v_LR2 = _mm512_add_ps(_mm512_mul_ps(v_LRout, v512_PoutlierXthetaBYmax), _mm512_mul_ps(v_LR, v512_PoutlierM1));
  } else {/* OUTLIER_TYPE==1 */
#if M512_ABS==0
    v_logerr2 = _mm512_sub_ps(_mm512_add_ps(v_logerrA, /* NEW */v_OutPen), _mm512_mul_ps(v_Bias, v512_biasWToutlierM1));
#else
    v_logerr2 = _mm512_sub_ps(_mm512_sub_ps(v_logerrA, /* NEW */v_OutPen), _mm512_mul_ps(v_Bias, v512_biasWToutlierM1));
#endif
    v_logerr2 = _mm512_add_ps(v_LogFnPown1, _mm512_min_ps(v_logerr2, v512_logMaxFAG2));// avoid overflow for large intervals

    v_Vexp2 = expF_mm512(v_logerr2);

    v_LRout = _mm512_mul_ps(v_FAik, v_Vexp2);
    if(OUT_OLD)
      v_LR2 = _mm512_add_ps(_mm512_mul_ps(v_LRout, v512_PoutlierPIbXtheta), _mm512_mul_ps(v_LR, v512_PoutlierM1));
    else
      v_LR2 = _mm512_add_ps(_mm512_mul_ps(v_LRout, v512_PoutlierPIBYmax), _mm512_mul_ps(v_LR, v512_PoutlierM1));
  }

  __m512 v_R = _mm512_mask_mov_ps(v_LR, 
#if M512_ABS==0
				  _mm512_cmp_ps_mask(v_fabsXmY, v512_outlierNegMaxBC, _MM_CMPINT_GT), // (-fabs(X-Y) > -outlierMaxBC) ? LR2 : LR;
#else
				  _mm512_cmp_ps_mask(v_fabsXmY, v512_outlierMaxBC, _MM_CMPINT_LT), // (fabs(X-Y) < outlierMaxBC) ? LR2 : LR;
#endif
				  v_LR2); 

#if VECTOR_DEBUG >= 2
  float Y[16]; _mm512_storeu_ps(Y, v_Y);
  float LogFnPown1[16]; _mm512_storeu_ps(LogFnPown1, v_LogFnPown1);
  float PrBiasW[16]; _mm512_storeu_ps(PrBiasW, v_PrBiasWT);
  float FAik[16]; _mm512_storeu_ps(FAik, v_FAik);
  float Ivar[16]; _mm512_storeu_ps(Ivar, v_Ivar);
  float RV[16]; _mm512_storeu_ps(RV, v_R);

#pragma novector
  for(int i = 0; i < 16; i++){
    MFLOAT R = FAJ(X, Y[i], m, n, FBiasWTxm, LogFpPowm1, LogFnPown1[i], PrBiasW[i], FAik[i], Ivar[i], outlier, OUT_TYPE, OUT_OLD);
    if((R > 0.0 || RV[i] > 0.0) && !(R >= 0.0 && RV[i] >= 0.0 && fabs(R - RV[i]) < R * 1e-4 + 1e-35)){
      #pragma omp critical
      {
	printf("\nFAJ_mm512:i=%d:m=%d,n=%d,X=%0.6f,Y=%0.6f,FBiasWTxm=%0.6e,LogFpPowm1=%0.6e,LogFnPown1=%0.6e,PrBiasWT=%0.6e,FAik=%0.6e,Ivar=%0.6e:RV[i]= %0.7e, R= %0.7e\n",
           i,m,n,X,Y[i],FBiasWTxm,LogFpPowm1,LogFnPown1[i],PrBiasW[i],FAik[i],Ivar[i],RV[i],R);
	MFLOAT err = (X - Y[i]) * Ivar[i];
	float errV[16]; _mm512_storeu_ps(errV, v_err);

#if M512_ABS==0
	MFLOAT OutPen = -fabsF(X - Y[i]) * OutlierLambdaInv;
#else
	MFLOAT OutPen = fabsF(X - Y[i]) * OutlierLambdaInv;
#endif
	float OutPenV[16]; _mm512_storeu_ps(OutPenV, v_OutPen);

	MFLOAT err2 = err*err;
	float err2V[16]; _mm512_storeu_ps(err2V, v_err2);

	printf("\t err= %0.7e(errV[i]= %0.7e), OutPen= %0.7e(OutPenV[i]= %0.7e), err2= %0.7e(err2V[i]= %0.7e)\n", err, errV[i], OutPen, OutPenV[i], err2, err2V[i]);
	
	float OutlierLambdaInvV[16]; _mm512_storeu_ps(OutlierLambdaInvV, v512_OutlierLambdaInv);
	float fabsXmYV[16]; _mm512_storeu_ps(fabsXmYV, v_fabsXmY);
#if M512_ABS==0
	printf("\t OutlierLambdaInv= %0.6e (OutlierLambdaInvV[i]= %0.6e), -fabsF(X-Y)= %0.6e(fabsXmYV[i]= %0.6e)\n",
	       OutlierLambdaInv, OutlierLambdaInvV[i], -fabsF(X-Y[i]), fabsXmYV[i]);
#else
	printf("\t OutlierLambdaInv= %0.6e (OutlierLambdaInvV[i]= %0.6e), fabsF(X-Y)= %0.6e(fabsXmYV[i]= %0.6e)\n",
	       OutlierLambdaInv, OutlierLambdaInvV[i], fabsF(X-Y[i]), fabsXmYV[i]);
#endif

	MFLOAT Bias = X * FrateBiasWT /* + PrBiasW[i] */ + FBiasWTxm;
	MFLOAT logerrA = Bias + LogFpPowm1;// WAS8 - X*Frate;
	if(!FRATE_FIX)
	  logerrA -= X*Frate;
	MFLOAT logerr = logerrA - err2;
	if(FRATE_FIX)
	  logerr -= X*Frate;
	logerr = LogFnPown1[i] + min(logerr, logMaxFAG);// avoid overflow for large intervals, especially with -biasWT 0
	float logerrAV[16]; _mm512_storeu_ps(logerrAV, v_logerrA);
	float logerrV[16]; _mm512_storeu_ps(logerrV, v_logerr);

	MFLOAT Vexp = expF(logerr);
	float VexpV[16]; _mm512_storeu_ps(VexpV, v_Vexp);

	MFLOAT LR = FAik[i] * Ivar[i] * Vexp;
	float LRV[16]; _mm512_storeu_ps(LRV, v_LR);

	float BiasV[16]; _mm512_storeu_ps(BiasV,v_Bias);
	float FrateBiasWTV[16]; _mm512_storeu_ps(FrateBiasWTV,v512_FrateBiasWT);

	printf("\t logerrA= %0.7e(logerrAV[i]= %0.7e), logerr= %0.7e(logerrV[i]= %0.7e), Vexp= %0.7e(VexpV[i]= %0.7e), LR= %0.7e (LRV[i]= %0.7e)\n", 
           logerrA, logerrAV[i], logerr, logerrV[i], Vexp, VexpV[i], LR, LRV[i]);
	printf("\t\t Bias= %0.7e(BiasV[i]= %0.7e), X= %0.7f, Frate= %0.7e, LogFpPowm1= %0.7e, logFnPown1[i]= %0.7e\n", Bias, BiasV[i], X, Frate, LogFpPowm1, LogFnPown1[i]);
	printf("\t\t FrateBiasWT=%0.7e (FrateBiasWTV[i]= %0.7e), FBiasWTxm= %0.7e\n", FrateBiasWT,FrateBiasWTV[i],FBiasWTxm);
	MFLOAT LR2;
	if(!OUT_TYPE){
#if M512_ABS==0
	  MFLOAT logerr2 = OutPen + Bias  * biasWToutlierF;
#else
	  MFLOAT logerr2 = Bias  * biasWToutlierF - OutPen;
#endif
	  float logerr2V[16]; _mm512_storeu_ps(logerr2V, v_logerr2);

	  MFLOAT LRout = expF(logerr2);
	  float LRoutV[16]; _mm512_storeu_ps(LRoutV, v_LRout);

	  if(OUT_OLD)
	    LR2 = LRout * PoutlierXthetaBYlambda + LR * PoutlierM1;
	  else
	    LR2 = LRout * PoutlierXthetaBYmax + LR * PoutlierM1;
	  float LR2V[16]; _mm512_storeu_ps(LR2V, v_LR2);

	  printf("\t logerr2= %0.7e(logerr2V[i]= %0.7e), LRout= %0.7e(LRoutV[i]= %0.7e), LR=%0.7e,LR2= %0.7e(LR2V[i]= %0.7e),OUT_OLD=%d\n",
		 logerr2,logerr2V[i],LRout, LRoutV[i], LR, LR2, LR2V[i], OUT_OLD);

	} else {/* OUTLIER_TYPE==1 */
#if M512_ABS==0
	  MFLOAT logerr2 = (logerrA + OutPen) - Bias * biasWToutlierM1;
#else
	  MFLOAT logerr2 = (logerrA - OutPen) - Bias * biasWToutlierM1;
#endif
	  logerr2 = LogFnPown1[i] + min(logerr2, logMaxFAG2);
	  float logerr2V[16]; _mm512_storeu_ps(logerr2V, v_logerr2);

	  MFLOAT Vexp2 = expF(logerr2);
	  float Vexp2V[16]; _mm512_storeu_ps(Vexp2V, v_Vexp2);

	  MFLOAT LRout = FAik[i] * Vexp2;
	  float LRoutV[16]; _mm512_storeu_ps(LRoutV, v_LRout);

	  if(OUT_OLD)
	    LR2 = LRout * PoutlierPIbXtheta + LR * PoutlierM1;
	  else
	    LR2 = LRout * PoutlierPIBYmax + LR * PoutlierM1;
	  float LR2V[16]; _mm512_storeu_ps(LR2V, v_LR2);

	  printf("\t logerr2= %0.7e(logerr2V[i]= %0.7e), Vexp2= %0.7f(Vexp2V[i]= %0.7e), LRout= %0.7e(LRoutV[i]= %0.7e), LR2= %0.7e(LR2V[i]= %0.7e),OUT_OLD=%d\n",
		 logerr2, logerr2V[i], Vexp2, Vexp2V[i], LRout, LRoutV[i], LR2, LR2V[i],OUT_OLD);
	}
	
#if M512_ABS==0
	float outlierNegMaxBCV[16]; _mm512_storeu_ps(outlierNegMaxBCV, v512_outlierNegMaxBC);
	printf("\t OUTLIER_DELTA(X-Y[i])= %d,fabsXmYV[i]= %0.6f, outlierMaxBC= %0.6e(outlierNegMaxBCV[i]= %0.6e)\n", 
	       OUTLIER_DELTA(X-Y[i]) ? 1 : 0, fabsXmYV[i], outlierMaxBC, outlierNegMaxBCV[i]);
#else
	float outlierMaxBCV[16]; _mm512_storeu_ps(outlierMaxBCV, v512_outlierMaxBC);
	printf("\t OUTLIER_DELTA(X-Y[i])= %d,fabsXmYV[i]= %0.6f, outlierMaxBC= %0.6e (outlierMaxBC[i]= %0.6e)\n", 
	       OUTLIER_DELTA(X-Y[i]) ? 1 : 0, fabsXmYV[i], outlierMaxBC, outlierMaxBCV[i]);
#endif

	//	MFLOAT Rb = FAJ(X, Y[i], m, n, FBiasWTxm, LogFpPowm1, LogFnPown1[i], FAik[i], Ivar[i], outlier, OUT_TYPE, OUT_OLD, 1);
	//	printf("\t Rb= %0.7e\n",Rb);
	fflush(stdout);
	assert(R >= 0.0 && RV[i] >= 0.0 && fabs(R - RV[i]) < R * 1e-4 + 1e-35);
      }
    }
  }

#endif

  return v_R;
}

/* FAJ_mm512 for use in qprobeval() : slightly different vector args */
template<int outlier, int OUT_TYPE, int OUT_OLD>
static inline __m512 FAJ_mm512(__m512 v_X,
			       __m512 v_Y,
			       __m512 v_FBiasWTxm,
			       __m512 v_LogFpPowm1,
			       __m512 v_LogFnPown1,
			       __m512 v_PrBiasWT,// NOTE : not used
			       __m512 v_FAik,/* previous return value of FAIK() */
			       __m512 v_Ivar) /* sqrt(1.0/(F2var+S2var*Y + ...)) computed in FAIK() */
{
  //  __m512 v_Bias = _mm512_add_ps(_mm512_mul_ps(v_X, v512_FrateBiasWT), _mm512_add_ps(v_PrBiasWT, v_FBiasWTxm));
  __m512 v_Bias = _mm512_add_ps(_mm512_mul_ps(v_X, v512_FrateBiasWT), v_FBiasWTxm);
  __m512 v_XmY = _mm512_sub_ps(v_X, v_Y);
  __m512 v_err = _mm512_mul_ps(v_Ivar, v_XmY);
#if 0 // HERE HERE M512_ABS // use mm512_abs
  __m512 v_fabsXmY = _mm512_abs_ps(v_XmY); // fabsF(X-Y)
#else
  __m512 v_fabsXmY = _mm512_sub_ps(_mm512_min_ps(v_X , v_Y), _mm512_max_ps(v_X, v_Y));// -fabsF(X-Y) = min(X,Y) - max(X,Y)
#endif
  __m512 v_OutPen = _mm512_mul_ps(v512_OutlierLambdaInv, v_fabsXmY);

  __m512 v_err2 = _mm512_mul_ps(v_err, v_err);

  __m512 v_logerrA, v_logerr;
  if(!FRATE_FIX){
    v_logerrA = _mm512_sub_ps(_mm512_add_ps(v_Bias, v_LogFpPowm1), _mm512_mul_ps(v_X, v512_Frate));  // logerrA =  Bias + LogFpPowm1 - X*Frate;
    v_logerr = _mm512_sub_ps(v_logerrA, v_err2); // logerr =  logerrA - err2;
  } else {/* FRATE_FIX */
    v_logerrA = _mm512_add_ps(v_Bias, v_LogFpPowm1);  // logerrA =  Bias + LogFpPowm1
    v_logerr = _mm512_sub_ps(_mm512_sub_ps(v_logerrA, v_err2), _mm512_mul_ps(v_X, v512_Frate)); // logerr =  logerrA - err2 - X*Frate;
  }

  v_logerr = _mm512_add_ps(v_LogFnPown1, _mm512_min_ps(v_logerr, v512_logMaxFAG));// avoid overflow for large intervals, especially with -biasWT 0

  __m512 v_Vexp = expF_mm512(v_logerr);
  __m512 v_LR = _mm512_mul_ps(_mm512_mul_ps(v_FAik, v_Ivar), v_Vexp);// result without outlier */

  if(!outlier)
    return v_LR;

  __m512 v_LR2, v_logerr2, v_LRout, v_Vexp2;
  if(!OUT_TYPE){
#if 0 // M512_ABS // using mm512_abs
    v_logerr2 = _mm512_sub_ps(_mm512_mul_ps(v_Bias, v512_biasWToutlierF), v_OutPen);
#else
    v_logerr2 = _mm512_add_ps(v_OutPen, _mm512_mul_ps(v_Bias, v512_biasWToutlierF));
#endif
    v_LRout = expF_mm512(v_logerr2);

    if(OUT_OLD)
      v_LR2 = _mm512_add_ps(_mm512_mul_ps(v_LRout, v512_PoutlierXthetaBYlambda), _mm512_mul_ps(v_LR, v512_PoutlierM1));
    else
      v_LR2 = _mm512_add_ps(_mm512_mul_ps(v_LRout, v512_PoutlierXthetaBYmax), _mm512_mul_ps(v_LR, v512_PoutlierM1));
  } else {/* OUTLIER_TYPE==1 */
#if 0 // M512_ABS // using mm512_abs_ps
    v_logerr2 = _mm512_sub_ps(_mm512_sub_ps(v_logerrA, /* NEW */v_OutPen), _mm512_mul_ps(v_Bias, v512_biasWToutlierM1));
#else
    v_logerr2 = _mm512_sub_ps(_mm512_add_ps(v_logerrA, /* NEW */v_OutPen), _mm512_mul_ps(v_Bias, v512_biasWToutlierM1));
#endif
    v_logerr2 = _mm512_add_ps(v_LogFnPown1, _mm512_min_ps(v_logerr2, v512_logMaxFAG2));// avoid overflow for large intervals

    v_Vexp2 = expF_mm512(v_logerr2);

    v_LRout = _mm512_mul_ps(v_FAik, v_Vexp2);
    if(OUT_OLD)
      v_LR2 = _mm512_add_ps(_mm512_mul_ps(v_LRout, v512_PoutlierPIbXtheta), _mm512_mul_ps(v_LR, v512_PoutlierM1));
    else
      v_LR2 = _mm512_add_ps(_mm512_mul_ps(v_LRout, v512_PoutlierPIBYmax), _mm512_mul_ps(v_LR, v512_PoutlierM1));
  }

  __m512 v_R = _mm512_mask_mov_ps(v_LR, 
#if 0 // M512_ABS // using mm512_abs_ps
				  _mm512_cmp_ps_mask(v_fabsXmY, v512_outlierMaxBC, _MM_CMPINT_LT), // (fabs(X-Y) < outlierMaxBC) ? LR2 : LR;
#else
				  _mm512_cmp_ps_mask(v_fabsXmY, v512_outlierNegMaxBC, _MM_CMPINT_GT), // (-fabs(X-Y) > outlierNegMaxBC) ? LR2 : LR;
#endif
				  v_LR2); 

#if VECTOR_DEBUG >= 2
  float Y[16]; _mm512_storeu_ps(Y, v_Y);
  float X[16]; _mm512_storeu_ps(X, v_X);
  float FBiasWTxm[16]; _mm512_storeu_ps(FBiasWTxm, v_FBiasWTxm);
  float LogFpPowm1[16]; _mm512_storeu_ps(LogFpPowm1, v_LogFpPowm1);
  float LogFnPown1[16]; _mm512_storeu_ps(LogFnPown1, v_LogFnPown1);
  float PrBiasW[16]; _mm512_storeu_ps(PrBiasW, v_PrBiasWT);
  float FAik[16]; _mm512_storeu_ps(FAik, v_FAik);
  float Ivar[16]; _mm512_storeu_ps(Ivar, v_Ivar);
  float RV[16]; _mm512_storeu_ps(RV, v_R);

  for(int i = 0; i < 16; i++){
    MFLOAT R = FAJ<outlier, OUT_TYPE, OUT_OLD>(X[i], Y[i], FBiasWTxm[i], LogFpPowm1[i], LogFnPown1[i], PrBiasW[i], FAik[i], Ivar[i]);
    if(isfinite(R) && isfinite(RV[i]) && (R > 0.0 || RV[i] > 0.0) && !(R >= 0.0 && RV[i] >= 0.0 && fabs(R - RV[i]) < R * 1e-4 + 1e-35)){
      #pragma omp critical
      {
	printf("\nFAJ_mm512:i=%d:X[i]=%0.7f,Y[i]=%0.7f,FBiasWTxm[i]=%0.7e,LogFpPowm1[i]=%0.7e,LogFnPown1[i]=%0.7e,PrBiasW[i]=%0.7e,FAik[i]=%0.7e,Ivar[i]=%0.7e:RV[i]= %0.7e, R= %0.7e\n",
	       i,X[i],Y[i],FBiasWTxm[i],LogFpPowm1[i],LogFnPown1[i],PrBiasW[i],FAik[i],Ivar[i],RV[i],R);
	MFLOAT err = (X[i] - Y[i]) * Ivar[i];
	float errV[16]; _mm512_storeu_ps(errV, v_err);

#if 0 // M512_ABS // using mm512_abs_ps
	MFLOAT OutPen = fabsF(X[i] - Y[i]) * OutlierLambdaInv;
#else
	MFLOAT OutPen = -fabsF(X[i] - Y[i]) * OutlierLambdaInv;
#endif
	float OutPenV[16]; _mm512_storeu_ps(OutPenV, v_OutPen);

	MFLOAT err2 = err*err;
	float err2V[16]; _mm512_storeu_ps(err2V, v_err2);

	printf("\t err= %0.7e(errV[i]= %0.7e), OutPen= %0.7e(OutPenV[i]= %0.7e), err2= %0.7e(err2V[i]= %0.7e)\n", err, errV[i], OutPen, OutPenV[i], err2, err2V[i]);
	
	float OutlierLambdaInvV[16]; _mm512_storeu_ps(OutlierLambdaInvV, v512_OutlierLambdaInv);
	float fabsXmYV[16]; _mm512_storeu_ps(fabsXmYV, v_fabsXmY);
	/*	printf("\t OutlierLambdaInv= %0.6e (OutlierLambdaInvV[i]= %0.6e), -fabsF(X-Y)= %0.6e(fabsXmYV[i]= %0.6e)\n",
		OutlierLambdaInv, OutlierLambdaInvV[i], -fabsF(X-Y[i]), fabsXmYV[i]);*/

	float BiasV[16]; _mm512_storeu_ps(BiasV, v_Bias);
	float FrateBiasWTV[16]; _mm512_storeu_ps(FrateBiasWTV, v512_FrateBiasWT);

	MFLOAT Bias = X[i] * FrateBiasWTV[i] /* + PrBiasW[i] */ + FBiasWTxm[i];
	MFLOAT logerrA = Bias + LogFpPowm1[i]; // WAS8 - X[i]*Frate;
	if(!FRATE_FIX)
	  logerrA -= X[i]*Frate;
	MFLOAT logerr = logerrA - err2;
	if(FRATE_FIX)
	  logerr -= X[i]*Frate;

	logerr = LogFnPown1[i] + min(logerr, logMaxFAG);// avoid overflow for large intervals, especially with -biasWT 0

	float logerrAV[16]; _mm512_storeu_ps(logerrAV, v_logerrA);
	float logerrV[16]; _mm512_storeu_ps(logerrV, v_logerr);

	MFLOAT Vexp = expF(logerr);
	float VexpV[16]; _mm512_storeu_ps(VexpV, v_Vexp);

	MFLOAT LR = FAik[i] * Ivar[i] * Vexp;
	float LRV[16]; _mm512_storeu_ps(LRV, v_LR);

	printf("\t logerrA= %0.7e(logerrAV[i]= %0.7e), logerr= %0.7e(logerrV[i]= %0.7e), Vexp= %0.7e(VexpV[i]= %0.7e), LR= %0.7e (LRV[i]= %0.7e)\n", 
           logerrA, logerrAV[i], logerr, logerrV[i], Vexp, VexpV[i], LR, LRV[i]);
	printf("\t\t Bias= %0.7e(BiasV[i]=%0.7e), X[i]= %0.7f, Frate= %0.7e, FrateBiasWT= %0.7e(FrateBiasWTV[i]= %0.7e), LogFpPowm1[i]= %0.7e, logFnPown1[i]= %0.7e\n", 
	       Bias, BiasV[i], X[i], Frate, FrateBiasWT, FrateBiasWTV[i], LogFpPowm1[i], LogFnPown1[i]);

	MFLOAT LR2;
	if(!OUT_TYPE){
#if 0 // M512_ABS // using mm512_abs_ps
	  MFLOAT logerr2 = Bias  * biasWToutlierF - OutPen;
#else
	  MFLOAT logerr2 = OutPen + Bias  * biasWToutlierF;
#endif
	  float logerr2V[16]; _mm512_storeu_ps(logerr2V, v_logerr2);

	  MFLOAT LRout = expF(logerr2);
	  float LRoutV[16]; _mm512_storeu_ps(LRoutV, v_LRout);

	  if(OUT_OLD)
	    LR2 = LRout * PoutlierXthetaBYlambda + LR * PoutlierM1;
	  else
	    LR2 = LRout * PoutlierXthetaBYmax + LR * PoutlierM1;
	  float LR2V[16]; _mm512_storeu_ps(LR2V, v_LR2);

	  printf("\t logerr2= %0.7e(logerr2V[i]= %0.7e), LRout= %0.7e(LRoutV[i]= %0.7e), LR=%0.7e,LR2= %0.7e(LR2V[i]= %0.7e),OUT_TYPE=%d,OUT_OLD=%d\n",
		 logerr2,logerr2V[i],LRout, LRoutV[i], LR, LR2, LR2V[i], OUT_TYPE,OUT_OLD);

	} else {/* OUTLIER_TYPE==1 */

#if 0 // M512_ABS // using mm512_abs_ps
	  MFLOAT logerr2 = (logerrA - OutPen) - Bias * biasWToutlierM1;
#else
	  MFLOAT logerr2 = (logerrA + OutPen) - Bias * biasWToutlierM1;
#endif
	  logerr2 = LogFnPown1[i] + min(logerr2, logMaxFAG2);
	  float logerr2V[16]; _mm512_storeu_ps(logerr2V, v_logerr2);

	  MFLOAT Vexp2 = expf(logerr2);//expF(logerr2);
	  float Vexp2V[16]; _mm512_storeu_ps(Vexp2V, v_Vexp2);

	  MFLOAT LRout = FAik[i] * Vexp2;
	  float LRoutV[16]; _mm512_storeu_ps(LRoutV, v_LRout);

	  if(OUT_OLD)
	    LR2 = LRout * PoutlierPIbXtheta + LR * PoutlierM1;
	  else
	    LR2 = LRout * PoutlierPIBYmax + LR * PoutlierM1;
	  float LR2V[16]; _mm512_storeu_ps(LR2V, v_LR2);

	  printf("\t logerr2= %0.7e(logerr2V[i]= %0.7e), Vexp2= %0.7e(Vexp2V[i]= %0.7e), LRout= %0.7e(LRoutV[i]= %0.7e), LR2= %0.7e(LR2V[i]= %0.7e),OUT_TYPE=%d,OUT_OLD=%d\n",
		 logerr2, logerr2V[i], Vexp2, Vexp2V[i], LRout, LRoutV[i], LR2, LR2V[i],OUT_TYPE,OUT_OLD);

	  float PoutlierPIBYmaxV[16]; _mm512_storeu_ps(PoutlierPIBYmaxV, v512_PoutlierPIBYmax);
	  float PoutlierM1V[16]; _mm512_storeu_ps(PoutlierM1V, v512_PoutlierM1);
	  printf("\t PoutlierPIBYmax= %0.8e(PoutlierPIBYmaxV[i]= %0.8e), PoutlierM1= %0.8e (PoutlierM1V[i]= %0.8e)\n",
		 PoutlierPIBYmax, PoutlierPIBYmaxV[i], PoutlierM1, PoutlierM1V[i]);
	}
	
#if 0 // M512_ABS // using mm512_abs_ps
	float outlierMaxBCV[16]; _mm512_storeu_ps(outlierMaxBCV, v512_outlierMaxBC);
	printf("\t OUTLIER_DELTA(X[i]-Y[i])= %d,fabsXmYV[i]= %0.6f, outlierMaxBC= %0.6e (outlierMaxBC[i]= %0.6e)\n", 
	       OUTLIER_DELTA(X[i]-Y[i]) ? 1 : 0, fabsXmYV[i], outlierMaxBC, outlierMaxBCV[i]);
#else
	float outlierNegMaxBCV[16]; _mm512_storeu_ps(outlierNegMaxBCV, v512_outlierNegMaxBC);
	printf("\t OUTLIER_DELTA(X[i]-Y[i])= %d,fabsXmYV[i]= %0.6f, outlierMaxBC= %0.6e(outlierNegMaxBCV[i]= %0.6e)\n", 
	       OUTLIER_DELTA(X[i]-Y[i]) ? 1 : 0, fabsXmYV[i], outlierMaxBC, outlierNegMaxBCV[i]);
#endif

	//	MFLOAT Rb = FAJv<outlier, OUT_TYPE, OUT_OLD>(X[i], Y[i], FBiasWTxm[i], LogFpPowm1[i], LogFnPown1[i], FAik[i], Ivar[i]);
	//	printf("\t Rb= %0.7f\n",Rb);
	fflush(stdout);
	assert(R >= 0.0 && RV[i] >= 0.0 && fabs(R - RV[i]) < R * 1e-4 + 1e-35);
      }
    }
  }

#endif

  return v_R;
}

#endif // VECTOR_MIC

#if 0
// Vector friendly version without debug args */

#pragma omp declare simd
template<int outlier, int OUT_TYPE, int OUT_OLD>
  static inline MFLOAT FA(MFLOAT X,
			  MFLOAT Y,
			  MFLOAT resR2,/* VARtab[K*N + I] : only used with RES_VARIANCE */
			  MFLOAT resL2,/* VARtab[T*N + I-K-n] : only used with RES_VARIANCE */
			  MFLOAT FBiasWTxm,
			  MFLOAT LogFpPowm1,
			  MFLOAT LogFnPown1,
			  MFLOAT PrBiasWT,// NOTE : not used
			  MFLOAT PRikn)
{
  MFLOAT var = F2var + S2var * Y;
  if(QUADRATIC_VARIANCE)
    var += R2var * Y * Y;
  if(RES_VARIANCE)
    var += E2var * (resR2 + resL2);

  if(DEBUG>=2) assert(var > 0.0);
  MFLOAT Ivar = sqrtF(((MFLOAT)1.0)/var);
  if(DEBUG>=2) assert(isfinite(Ivar));
  if(DEBUG>=2) assert(PRikn <= 1.0);
  MFLOAT FAik = RefPen * PRikn;// NOTE : FnPow[n-1] has been moved into expF(LogFnPow[n-1])

  MFLOAT err = (X-Y) * Ivar;
  MFLOAT Bias = X * FrateBiasWT /* + PrBiasWT */ + FBiasWTxm;
  MFLOAT OutPen = -fabsF(X-Y) * OutlierLambdaInv;
  if(DEBUG>=2 && biasWT == 0.0) assert(Bias == 0.0);

  MFLOAT err2 = err*err;
  MFLOAT logerrA = Bias + LogFpPowm1;// WAS8 - X*Frate;
  if(!FRATE_FIX)
    logerrA -= X*Frate;
  MFLOAT logerr = logerrA - err2;
  if(FRATE_FIX)
    logerr -= X*Frate;

  // NOTE MaxFA cannot be applied to Fn ^^ (n-1), since mprobeval assumes FAJ changes by exactly a multiple of Fn if n changes by 1 
  MFLOAT logerrb = LogFnPown1 + min(logerr, logMaxFAG);// avoid overflow for large intervals, especially with -biasWT 0

  MFLOAT Vexp = expF(logerrb);
  if(DEBUG>=2) assert(isfinite(Vexp) && Vexp <= MaxFAG);
  MFLOAT LR = FAik * Ivar * Vexp;
  if(DEBUG>=2) assert(isfinite(LR) && LR <= MaxFA);

  if(outlier && OUTLIER_DELTA(X-Y)){
    if(!OUT_TYPE){
      MFLOAT LRout = expF(OutPen + Bias * biasWToutlierF);
      if(DEBUG>=2) assert(isfinite(LRout) && LRout <= MaxFA);
      if(OUT_OLD)
	return LRout * PoutlierXthetaBYlambda + LR * PoutlierM1;
      else
	return LRout * PoutlierXthetaBYmax + LR * PoutlierM1;
    } else {
      MFLOAT logerr2 = logerrA - Bias * biasWToutlierM1 + OutPen;
      MFLOAT logerr2b = LogFnPown1 + min(logerr2, logMaxFAG2);// avoid overflow for large intervals
      MFLOAT Vexp2 = expF(logerr2b);
      if(DEBUG>=2) assert(isfinite(Vexp2) && Vexp2 <= MaxFAG2);
      MFLOAT LRout = FAik * Vexp2;
      if(DEBUG>=2 && rverb){
        rverb = 0;
	printf("\t loggerr2=%0.8e,Vexp2=%0.8e,LRout=%0.8e,PoutlierPIBYmax=%0.8e,PoutlierM1=%0.8e\n",
            logerr + err2 - Bias * biasWToutlierM1, Vexp2, LRout, PoutlierPIBYmax,PoutlierM1);
	fflush(stdout);
	rverb = 1;
      }
      if(DEBUG>=2) assert(isfinite(LRout) && LRout <= MaxFA);

      if(OUT_OLD)
	return LRout * PoutlierPIbXtheta + LR * PoutlierM1;
      else
	return LRout * PoutlierPIBYmax + LR * PoutlierM1;
    }
  } else
    return LR;
}

#endif // 0

#if USE_MFLOAT==0
#define FLJf FLJ
#else // use MFLOAT==float instead of double
static inline MFLOAT FLJf(int J,
			  MFLOAT x, /* Xm[J] */
			  MFLOAT &IvarX)/* return value 1/sqrt(var) */
{
  if(DEBUG>=2 && !(J-1 <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX)){
    printf("FLJ(J=%d):DELTA_X=%d,deltaExtX=%d,EXTEND_MAX=%d\n",J,DELTA_X,deltaExtXRef,EXTEND_MAX);
    fflush(stdout);
    assert(J-1 <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX);
  }
  if(DEBUG>=2) assert(J > 0);
  MFLOAT Bias = x * FrateBiasWT + (J-1) * FBiasWT + EBiasWT;
  if(DEBUG>=2 && biasWT <= 0.0) assert(Bias <= 0.0);

  MFLOAT var = F2var + S2var * x;
  if(QUADRATIC_VARIANCE)
    var += R2var * x * x;
  IvarX = ((MFLOAT)1.0)/sqrtF(((MFLOAT)2.0) * var);

  MFLOAT logerr = -x*Frate + Bias + LogFpPow[J-1] ;
  logerr = min(logerr, logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0 (NOTE : slightly more stringent than for FR, to support FLJ)
  MFLOAT Vexp = expF(logerr);
  if(DEBUG>=2 && (!isfinite(Vexp) || rverb) ){
    printf("FLJ:J=%d,x=%0.6f:var=%0.8e,IvarX=%0.8e,LogFpPow[J-1]=%0.8f,Frate=%0.8e,Bias=%0.8e,expF(-x*Frate+Bias+LogFpPow[J-1])=%0.8e\n",
	   J,x,var,IvarX,LogFpPow[J-1],Frate,Bias,Vexp);
    fflush(stdout);
    assert(isfinite(Vexp));
  }

  return Vexp * Rtheta;
}
#endif

#if 0
/* FLE() corresponds to the end alignment gaussian term G(D[J],y) of FL() (U==0) */
static inline double FLE(int I,
			 int K,
			 int J,
			 int U,
			 int V,
			 double y,/* value of Yc(H,I,K) */
			 double *D,/* Data map */
			 double *H)/* Reference map */
{
  if(DEBUG>=2) assert(U==0);
  if(DEBUG>=2) assert(J-1 <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX);
  if(DEBUG>=2) assert(I-K-V <= gN);
  double x = D[J];
  double Bias = x * FrateBiasWT + (J-1) * FBiasWT + EBiasWT;
  if(DEBUG>=2 && biasWT <= 0.0) assert(Bias <= 0.0);

  return expF(-x*Frate + Bias) * FnPow[I-K-V] * FpPow[J-1] * G(D[J],y) * InvRtheta;
}
#endif

/* FL() corresponds to term FL in White Paper without the end alignment gaussian term G(D[J],y) (NOTE: ignores RES_VARIANCE) */
/* For greater speed, can be computed in parts : 
MFLOAT IvarX;
MFLOAT FLj = FLJ(J,x,IvarX);
...
MFLOAT FL =  FLj * FLIK(I,K,J,U,V,y,D[J],y-H[U],H[V]-H[U],H,IvarX,lc); */

#pragma GCC push_options
#pragma GCC optimize("no-associative-math")

template<int ERFC>
static inline MFLOAT FL(int I,
			int K,
			int J,
			int U,
			int V,
			double y,/* value Yc(Y,I,K) : debug only */
			MFLOAT x,/* value of Xm[J] */
			MFLOAT Yiu,/* Yc(Y,I,K) - Y[U] */
			MFLOAT Yvu,/* Y[V] - Y[U] */
			double *Y,/* Reference map : debug only */
			int lc)/* 1 IFF left end if Y is chromosome end */
{
  if(DEBUG>=2) assert(J-1 <= DELTA_X + EXTEND_MAX);
  if(DEBUG>=2) assert(I-K-V <= gN);
  if(DEBUG>=2) assert(END_FIX || lc);
  if(DEBUG>=2) assert(V <= I-K);

  MFLOAT Bias = x * FrateBiasWT + FBiasWTm[J-1] + EBiasWT;
  if(DEBUG>=2 && biasWT <= 0.0) assert(Bias <= 0.0);

  MFLOAT logerr = -x*Frate + Bias + LogFpPow[J-1];
  logerr = min(logerr, logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0 (NOTE : slightly more stringent than for FR, to support FLJ)
  MFLOAT Vexp = expF(logerr);
  if(DEBUG>=2) assert(isfinite(Vexp));

  if(END_FIX>=2 && V==I-K){
    if(DEBUG>=2) assert(y >= Y[U]);
    return  Vexp * Rtheta * GE<ERFC>(x, Yiu, Yiu, (EXTEND_MAX>0 && !lc && U<=0 ) ? 1 : 0);
  }

  if(DEBUG>=2) assert(U <= V);
  if(DEBUG>=2) assert(Y[U] <= Y[V]);
  if(DEBUG>=2) assert(0 <= Yvu && Yvu <= Yiu);
  MFLOAT GEval =  Rtheta * GE<ERFC>(x, Yiu, Yvu, (EXTEND_MAX>0 && !lc && U<=0 ) ? 1 : 0);

  if(FLIK_FIX) {
    return (Vexp * FnPowf[I-K-V]) * GEval;
  } else
    return Vexp * FnPowf[I-K-V] * GEval;
}
//#pragma GCC pop_options

//#pragma GCC push_options
//#pragma GCC optimize("no-associative-math")

/* FLE() / FLJ() */
static inline MFLOAT FLEIK(int I,
			   int K,
			   int J,
			   int U,
			   int V,
			   MFLOAT y,/* value of Yc(H,I,K) */
			   MFLOAT x,/* value of Xm[J] */
			   MFLOAT FLj) // only used with FLIK_FIX to avoid underflow due to FnPowf[I-K-V]
{
  if(DEBUG>=2) assert(U==0);
  if(DEBUG>=2) assert(I-K-V <= gN);
  MFLOAT Gval = G(x,y) * InvRtheta2;

  if(FLIK_FIX)
    return (FLj * FnPowf[I-K-V]) * Gval;
  else
    return FnPowf[I-K-V] * Gval;
}

// #pragma GCC pop_options

// #pragma GCC push_options
// #pragma GCC optimize("no-associative-math")

template<int ERFC>
static inline MFLOAT FLIK(int I,
			  int K,
			  int J,
			  int U,
			  int V,
			  double y,/* value Yc(Y,I,K) : debug only */
			  MFLOAT x,/* Xm[J] */
			  MFLOAT Yiu,/* Yc(Y,I,K) - Y[U] */
			  MFLOAT Yvu,/* Y[V] - Y[U] */
			  double *Y,/* Reference map : debug only */
			  MFLOAT IvarX,
			  MFLOAT FLj,// only used with FLIK_FIX to multiply result and avoid underflow
			  int lc)/* 1 IFF left end if H is chromosome end */
{
  if(DEBUG>=2) assert(isfinite(FLj) && 0.0f <= FLj && FLj <= MaxUR * 3.0);
  if(DEBUG>=2) assert(I-K-V <= gN);
  if(DEBUG>=2) assert(END_FIX || lc);
  if(DEBUG>=2 && !(V <= I-K)){
    #pragma omp critical
    {
      printf("FLIK:I=%d,K=%d,J=%d,U=%d,V=%d,y=%0.6f,x=%0.6f,Yiu=%0.6f,Yvu=%0.6f,Y[U]=%0.6f,Y[V]=%0.6f,IvarX=%0.8e,lc=%d\n",
	     I,K,J,U,V,y,x,Yiu,Yvu,Y[U],Y[V],IvarX,lc);
      fflush(stdout);
      assert(V <= I-K);
    }
  }
  if(DEBUG>=2) assert(U <= V);
  if(DEBUG>=2) assert(Y[U] <= Y[V]);
  if(END_FIX>=2 && V==I-K){
    if(DEBUG>=2) assert(Y[U] <= y);
    if(DEBUG>=2) assert(0.0f <= Yiu);
    MFLOAT GEval = GE<ERFC>(x, Yiu, Yiu, IvarX, (EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0);
    if(DEBUG>=2) assert(GEval <= (MFLOAT)1.0001);
    return FLIK_FIX ? FLj * GEval : GEval;
  }
  if(VERB>=2 && rverb >= 2){
    rverb = 0;// avoid triggering another printout from GE()
    printf("\nFLIK:I=%d,K=%d,J=%d,U=%d,V=%d,y=%0.6f,Yiu=%0.6f,Yvu=%0.6f,IvarX=%0.8e,FnPow[I-K-V]=%0.8e,GE(x,Yiu,Yvu,IvarX,%d)=%0.8e\n",
	   I,K,J,U,V,y,Yiu,Yvu,IvarX,FnPow[I-K-V],(EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0, GE<ERFC>(x, Yiu, Yvu, IvarX, (EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0));
    fflush(stdout);
    rverb = 2;
  }
  if(DEBUG>=2) assert(0 <= Yvu && Yvu <= Yiu);

  MFLOAT GEval = GE<ERFC>(x, Yiu, Yvu, IvarX, (EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0);
  if(DEBUG>=2) assert(GEval <= (MFLOAT)1.0001);

  if(FLIK_FIX) {
    if(DEBUG>=2) assert(isfinite(FnPowf[I-K-V]) && (MFLOAT)0.0 <= FnPowf[I-K-V] && FnPowf[I-K-V] <= (MFLOAT)1.0);
    return (FLj * FnPowf[I-K-V]) * GEval;
  } else
    return FnPowf[I-K-V] * GEval;
}

template<int ERFC>
static inline double FLIKd(int I,
			   int K,
			   int J,
			   int U,
			   int V,
			   double y,/* value Yc(Y,I,K) : debug only */
			   MFLOAT2 x,/* Xm[J] */
			   MFLOAT2 Yiu,/* Yc(Y,I,K) - Y[U] */
			   MFLOAT2 Yvu,/* Y[V] - Y[U] */
			   double *Y,/* Reference map : debug only */
			   MFLOAT2 IvarX,
			   double FLj,// only used with FLIK_FIX to multiply result and avoid underflow
			   int lc)/* 1 IFF left end if H is chromosome end */
{
  if(DEBUG>=2) assert(I-K-V <= gN);
  if(DEBUG>=2) assert(END_FIX || lc);
  if(DEBUG>=2 && !(V <= I-K)){
    #pragma omp critical
    {
      printf("FLIK:I=%d,K=%d,J=%d,U=%d,V=%d,y=%0.6f,x=%0.6f,Yiu=%0.6f,Yvu=%0.6f,Y[U]=%0.6f,Y[V]=%0.6f,IvarX=%0.8e,lc=%d\n",
	     I,K,J,U,V,y,x,Yiu,Yvu,Y[U],Y[V],IvarX,lc);
      fflush(stdout);
      assert(V <= I-K);
    }
  }
  if(DEBUG>=2) assert(U <= V);
  if(DEBUG>=2) assert(Y[U] <= Y[V]);
  if(END_FIX>=2 && V==I-K){
    if(DEBUG>=2) assert(Y[U] <= y);
    if(DEBUG>=2) assert(0.0f <= Yiu);
    MFLOAT GEval = GE<ERFC>(x, Yiu, Yiu, IvarX, (EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0);
    if(DEBUG>=2) assert(GEval <= (MFLOAT)1.0001);
    return FLIK_FIX ? FLj * GEval : GEval;
  }
  if(VERB>=2 && rverb >= 2){
    rverb = 0;// avoid triggering another printout from GE()
    printf("\nFLIK:I=%d,K=%d,J=%d,U=%d,V=%d,y=%0.6f,Yiu=%0.6f,Yvu=%0.6f,IvarX=%0.8e,FnPow[I-K-V]=%0.8e,GE(x,Yiu,Yvu,IvarX,%d)=%0.8e\n",
	   I,K,J,U,V,y,Yiu,Yvu,IvarX,FnPow[I-K-V],(EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0, GE<ERFC>(x, Yiu, Yvu, IvarX, (EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0));
    fflush(stdout);
    rverb = 2;
  }
  if(DEBUG>=2) assert(0 <= Yvu && Yvu <= Yiu);

  MFLOAT GEval = GE<ERFC>(x, Yiu, Yvu, IvarX, (EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0);
  if(DEBUG>=2) assert(GEval <= (MFLOAT)1.0001);

  if(FLIK_FIX) {
    if(DEBUG>=2) assert(isfinite(FnPowf[I-K-V]) && (MFLOAT)0.0 <= FnPowf[I-K-V] && FnPowf[I-K-V] <= (MFLOAT)1.0);
    return (FLj * FnPow[I-K-V]) * GEval;
  } else
    return FnPow[I-K-V] * GEval;
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT

/* AVX256 version of FLEIK */
static inline __m256 FLEIK_mm256(int I,
				 int K,
				 int U,
				 int V,
				 __m256 Yik,/* value of Yc(H,I,K) */
				 __m256 Xj,/* value of Xm[J] */
				 __m256 FLj) // only used with FLIK_FIX to avoid underflow due to FnPowf[I-K-V]
{
  if(FLIK_FIX) {
    //    __m256 Gval = InvRtheta2 * G_mm256(Xj, Yik);
    __m256 Gval = v256_InvRtheta2 * G_mm256(Xj, Yik);
    return (FLj * _mm256_set1_ps(FnPowf[I-K-V])) *  Gval;
  } else
    return (_mm256_set1_ps(FnPowf[I-K-V]) * v256_InvRtheta2) * G_mm256(Xj, Yik);
}

/* AVX256 version of FLIK */
template<int ERFC>
static inline __m256 FLIK_mm256(int I,
				int K,
				int U,
				int V,
				__m256 Xj,/* Xm[J] */
				__m256 Yiu,/* Yc(Y,I,K) - Y[U] */
				__m256 Yvu,/* Y[V] - Y[U] */
				__m256 IvarX,
				__m256 FLj,// only used with FLIK_FIX to multiply result and avoid underflow
				int lc)/* 1 IFF left end if H is chromosome end */
{
  if(END_FIX>=2 && V==I-K){
    __m256 GEval = GE_mm256<ERFC>(Xj, Yiu, Yiu, IvarX, (EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0);
    return FLIK_FIX ? FLj * GEval : GEval;
  }

  __m256 GEval = GE_mm256<ERFC>(Xj, Yiu, Yvu, IvarX, (EXTEND_MAX > 0 && U <= 0 && !lc) ? 1 : 0);
  if(FLIK_FIX)
    return (FLj * _mm256_set1_ps(FnPowf[I-K-V])) * GEval;
  else
    return _mm256_set1_ps(FnPowf[I-K-V]) * GEval;
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

/* AVX512 version of FLEIK */
static inline __m512 FLEIK_mm512(int I,
				 int K,
				 int U,
				 int V,
				 __m512 Yik,/* value of Yc(H,I,K) */
				 __m512 Xj,/* value of Xm[J] */
				 __m512 FLj) // only used with FLIK_FIX to avoid underflow due to FnPowf[I-K-V]
{
  if(FLIK_FIX) {
    __m512 Gval = _mm512_mul_ps(_mm512_set1_ps(InvRtheta2), G_mm512(Xj, Yik));
    return _mm512_mul_ps(_mm512_mul_ps(FLj, _mm512_set1_ps(FnPowf[I-K-V])), Gval);
  } else
    return _mm512_mul_ps(_mm512_set1_ps(FnPowf[I-K-V] * InvRtheta2) , G_mm512(Xj,Yik));
}

/* AVX512 version of FLIK */
template<int ERFC>
static inline __m512 FLIK_mm512(int I,
				int K,
				int U,
				int V,
				__m512 Xj,/* Xm[J] */
				__m512 Yiu,/* Yc(Y,I,K) - Y[U] */
				__m512 Yvu,/* Y[V] - Y[U] */
				__m512 IvarX,
				__m512 FLj,// only used with FLIK_FIX to multiply result and avoid underflow
				int lc)/* 1 IFF left end if H is chromosome end */
{
  if(END_FIX>=2 && V==I-K){
    __m512 GEval = GE_mm512<ERFC>(Xj, Yiu, Yiu, IvarX, (EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0);
    return FLIK_FIX ? _mm512_mul_ps(FLj, GEval) : GEval;
  }

  __m512 FnPowfIKV = _mm512_set1_ps(FnPowf[I-K-V]);
  __m512 GEval = GE_mm512<ERFC>(Xj, Yiu, Yvu, IvarX, (EXTEND_MAX > 0 && U <= 0 && !lc) ? 1 : 0);

  if(FLIK_FIX)
    return _mm512_mul_ps(_mm512_mul_ps(FLj, FnPowfIKV), GEval);
  else
    return _mm512_mul_ps(FnPowfIKV, GEval);
}

#endif // VECTOR_AVX || VECTOR_MIC

// #pragma GCC_pop_options

// #pragma GCC push_options
// #pragma GCC optimize("no-associative-math")

/* FLIKA() corresponds to FLIK() with location of H[U] or H[V] modified to an added site */
template<int ERFC>
static inline double FLIKA(int I,
			    int K,
			    int J,
			    int U,
			    int V,
			    double y,/* value Yc(Y,I,K) : debug only */
			    MFLOAT2 Xj, /* Xm[J] */
			    double *Y,/* Reference map */
			    double Yu,/* Reference map location Y[U] (as modified by new added site location) : debug only */
			    double Yv,/* Reference map location Y[V] (as modified by new added site location) : debug only */
			    MFLOAT2 Yiu,/* Yc(Y,I,K) - Yu */
			    MFLOAT2 Yvu,/* Yv - Yu */
			    MFLOAT2 IvarX,
			    double FLj,// HERE : try MFLOAT2 FLj : add assertion that FLJ <= MaxUR * 3. 
			    int lc)/* 1 IFF left end if H is chromosome end */
 {
  if(DEBUG>=2) assert(I-K-V <= gN);
  if(DEBUG>=2) assert(END_FIX || lc);
  if(DEBUG>=2) assert(Yv <= Y[I-K]);
  if(END_FIX>=2 && Yv >= Y[I-K]){
    if(DEBUG>=2) assert(Yu <= y);
    MFLOAT GEval = GE<ERFC>(Xj, Yiu, Yiu, IvarX, (EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0);
    return FLIK_FIX ? FLj * GEval : GEval;
  }
  if(DEBUG>=2) assert(Yu <= Yv);
  if(DEBUG>=2) assert(0 <= Yvu && Yvu <= Yiu);
  MFLOAT GEval = GE<ERFC>(Xj, Yiu, Yvu, IvarX, (EXTEND_MAX>0 && U <= 0 && !lc) ? 1 : 0);
  
  if(FLIK_FIX) {
    return (FLj * FnPowf[I-K-V]) * GEval;
  } else
    return FnPowf[I-K-V] * GEval;
}

#pragma GCC pop_options

static inline double FLJ(int J,
			 double x, /* Xm[J] */
			 double &IvarX) /* return value 1/sqrt(var) */
{
  if(DEBUG>=2 && !(J-1 <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX)){
    printf("FLJ(J=%d):DELTA_X=%d,deltaExtX=%d,EXTEND_MAX=%d\n",J,DELTA_X,deltaExtXRef,EXTEND_MAX);
    fflush(stdout);
    assert(J-1 <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX);
  }
  if(DEBUG>=2) assert(J > 0);
  double Bias = x * FrateBiasWT + (J-1) * FBiasWT + EBiasWT;
  if(DEBUG>=2 && biasWT <= 0.0) assert(Bias <= 0.0);

  double var = F2var + S2var * x;
  if(QUADRATIC_VARIANCE)
    var += R2var * x * x;
  IvarX = 1.0/sqrt(2.0 * var);

  double logerr = -x*Frate + Bias + LogFpPow[J-1];
  logerr = min(logerr, (double)logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0 (NOTE : slightly more stringent than for FR, to support FLJ)
  double Vexp = exp(logerr);
  if(DEBUG>=2 && (!isfinite(Vexp) || rverb) ){
    printf("FLJ:J=%d,x=%0.6f:var=%0.8e,IvarX=%0.8e,FpPow[J-1]=%0.8f,Frate=%0.8e,Bias=%0.8e,expF(-x*Frate+Bias)=%0.8e\n",
	   J,x,var,IvarX,FpPow[J-1],Frate,Bias,Vexp);
    fflush(stdout);
    assert(isfinite(Vexp));
  }

  return Vexp * Rtheta;
}

#pragma GCC push_options
#pragma GCC optimize("no-associative-math")

/* FR() corresponds to the term FR in White Paper without the end alignment gaussian term G(D[M+1]-D[J],H[U]-y) */
/* NOTE : ignores RES_VARIANCE */
template<int ERFC>
static inline MFLOAT FR(int I,
			int K,
			int J,
			int U,
			int V,
			double y,/* value Yc(H,I,K) : debug only */
			MFLOAT x,/* D[M+1] - D[J] */
			MFLOAT Yiu,/* H[U] - y */
			MFLOAT Yvu,/* H[U] - H[V] */
			double *D,/* Data map : debug only */
			double *H, /* Reference map : debug only */
			int M, /* number of sites in Data map */
			int N, /* number of sites in Reference map */
			int rc) /* If right end of H is a chromosome end */
{
  if(DEBUG>=2) assert(M-J <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX);
  if(DEBUG>=2) assert(V-I <= N);
  if(DEBUG>=2) assert(END_FIX || rc);
  if(DEBUG>=2 && !(V >= I)){
    printf("FR:I=%d,K=%d,J=%d,U=%d,V=%d,y=%0.6f,x=%0.6f,Yiu=%0.6f,Yvu=%0.6f,M=%d,N=%d,rc=%d\n",
	   I,K,J,U,V,y,x,Yiu,Yvu,M,N,rc);
    fflush(stdout);
    assert(V >= I);
  }
  if(DEBUG>=2) assert(V <= U);
  if(DEBUG>=2) assert(U <= N+1);

  //  MFLOAT x = D[M+1] - D[J];
  MFLOAT Bias = x * FrateBiasWT + FBiasWTm[M-J] + EBiasWT;
  if(DEBUG>=2 && biasWT == 0.0 && !(Bias == 0.0)){
    printf("FR: biasWT=%0.6e, FrateBiasWT=%0.6e,FBiasWT=%0.6e,X=%0.6e,M=%d,J=%d:Bias=%0.6e\n",
	   biasWT,FrateBiasWT,FBiasWT,x,M,J,Bias);
    fflush(stdout);
    assert(Bias == 0.0);
  }

  if(END_FIX>=2 && V==I){
    if(DEBUG>=2) assert(y <= H[U]);
    MFLOAT logerr = -x * Frate + Bias + LogFpPow[M-J];
    logerr = min(logerr, logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0
    MFLOAT Vexp = expF(logerr);
    if(DEBUG>=2) assert(isfinite(Vexp));
    MFLOAT GEval = GE<ERFC>(x,Yiu,Yiu, (EXTEND_MAX>0 && !rc && U > N) ? 1 : 0);
    return (Vexp * RthetaF) * GEval;
  }
  if(DEBUG>=2) assert(H[V] <= H[U]);

  if(DEBUG>=2 && rverb){
    printf("FR:I=%d,K=%d,N=%d,J=%d,M=%d,U=%d,V=%d,x=%0.6f,y=%0.6f,H[U]=%0.6f,H[V]=%0.6f:expF=%0.8e(Frate=%0.8e,Bias=%0.8e),FnPow[V-I]=%0.8e,FpPow[M-J]=%0.8e\n",
	   I,K,N,J,M,U,V,x,y,H[U],H[V],expF(-x*Frate+Bias),Frate,Bias,FnPow[V-I],FpPow[M-J]);
    printf("\t GE:H[U]-y=%0.6f,H[U]-H[V]=%0.6f,rc=%d:GE=%0.8e\n",H[U]-y,H[U]-H[V],rc,
	   GE<ERFC>(x,Yiu, Yvu, (EXTEND_MAX>0 && !rc && U > N) ? 1 : 0));
    fflush(stdout);
  }

  MFLOAT logerr = -x * Frate + Bias + LogFpPow[M-J];

  // NOTE MaxUR cannot be applied to Fn ^^ (n-1), since mprobeval assumes FR changes by exactly a multiple of Fn if V changes by 1 
  logerr = LogFnPow[V-I] + min(logerr, logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0

  MFLOAT Vexp = expF(logerr);
  if(DEBUG>=2) assert(isfinite(Vexp));
  if(DEBUG>=2) assert(0 <= Yvu && Yvu <= Yiu);
  MFLOAT GEval = GE<ERFC>(x,Yiu, Yvu, (EXTEND_MAX>0 && !rc && U > N) ? 1 : 0);

  return (Vexp * RthetaF) * GEval;
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT && DEBUG < 2

/* AVX256 version of FR() */
template<int ERFC>
static inline __m256 FR_mm256(int I,
				int K,
				int J,
				int U,
				int V,
				__m256 x,/* Xm[M+1] - Xm[J] */
				__m256 Yiu,/* H[U] - Yik */
				__m256 Yvu,/* H[U] - H[V] */
				int M, /* number of sites in Data map */
				int N, /* number of sites in Reference map */
				int rc) /* If right end of H is a chromosome end */
{
  if(DEBUG>=2) assert(M-J <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX);
  if(DEBUG>=2) assert(V-I <= N);
  if(DEBUG>=2) assert(END_FIX || rc);
  if(DEBUG>=2 && !(V >= I)){
    printf("FR_mm256:I=%d,K=%d,J=%d,U=%d,V=%d,M=%d,N=%d,rc=%d\n",
	   I,K,J,U,V,M,N,rc);
    fflush(stdout);
    assert(V >= I);
  }
  if(DEBUG>=2) assert(V <= U);
  if(DEBUG>=2) assert(U <= N+1);

  __m256 Bias = x * v256_FrateBiasWT + _mm256_set1_ps(FBiasWTm[M-J] + EBiasWT);

  if(END_FIX>=2 && V==I){
    __m256 logerr = Bias + _mm256_set1_ps(LogFpPow[M-J]) - x * v256_Frate ;
    logerr = _mm256_min_ps(logerr, v256_logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0
    __m256 Vexp = expF_mm256(logerr);
    return Vexp * v256_RthetaF * GE_mm256<ERFC>(x, Yiu, Yiu, (EXTEND_MAX > 0 && !rc && U > N) ? 1 : 0);
  }
  __m256 logerr = Bias + _mm256_set1_ps(LogFpPow[M-J]) - x * v256_Frate ;

  logerr = _mm256_set1_ps(LogFnPow[V-I]) + _mm256_min_ps(logerr, v256_logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0
  __m256 Vexp = expF_mm256(logerr);
  __m256 GEval = GE_mm256<ERFC>(x, Yiu, Yvu, (EXTEND_MAX > 0 && !rc && U > N) ? 1 : 0);
  return (Vexp * v256_RthetaF) * GEval;
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

/* AVX512 version of FR() */
template<int ERFC>
static inline __m512 FR_mm512(int I,
				int K,
				int J,
				int U,
				int V,
				__m512 x,/* Xm[M+1] - Xm[J] */
				__m512 Yiu,/* H[U] - Yik */
				__m512 Yvu,/* H[U] - H[V] */
				int M, /* number of sites in Data map */
				int N, /* number of sites in Reference map */
				int rc) /* If right end of H is a chromosome end */
{
  if(DEBUG>=2) assert(M-J <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX);
  if(DEBUG>=2) assert(V-I <= N);
  if(DEBUG>=2) assert(END_FIX || rc);
  if(DEBUG>=2 && !(V >= I)){
    printf("FR_mm256:I=%d,K=%d,J=%d,U=%d,V=%d,M=%d,N=%d,rc=%d\n",
	   I,K,J,U,V,M,N,rc);
    fflush(stdout);
    assert(V >= I);
  }
  if(DEBUG>=2) assert(V <= U);
  if(DEBUG>=2) assert(U <= N+1);

  __m512 Bias = _mm512_add_ps(_mm512_mul_ps(x , v512_FrateBiasWT), _mm512_set1_ps(FBiasWTm[M-J] + EBiasWT));

  if(END_FIX>=2 && V==I){
    __m512 logerr = _mm512_sub_ps(_mm512_add_ps(Bias, _mm512_set1_ps(LogFpPow[M-J])), _mm512_mul_ps(x, v512_Frate)) ;
    logerr = _mm512_min_ps(logerr, v512_logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0
    __m512 Vexp = expF_mm512(logerr);
    return _mm512_mul_ps(_mm512_mul_ps(Vexp, v512_RthetaF), GE_mm512<ERFC>(x, Yiu, Yiu, (EXTEND_MAX > 0 && !rc && U > N) ? 1 : 0));
  }
  __m512 logerr = _mm512_sub_ps(_mm512_add_ps(Bias, _mm512_set1_ps(LogFpPow[M-J])), _mm512_mul_ps(x, v512_Frate)) ;

  logerr = _mm512_add_ps(_mm512_set1_ps(LogFnPow[V-I]), _mm512_min_ps(logerr, v512_logMaxUR));// avoid overflow for large end intervals, especially with -biasWT 0
  __m512 Vexp = expF_mm512(logerr);

  return _mm512_mul_ps(_mm512_mul_ps(Vexp, v512_RthetaF), GE_mm512<ERFC>(x, Yiu, Yvu, (EXTEND_MAX > 0 && !rc && U > N) ? 1 : 0));
}

#endif // VECTOR_AVX || VECTOR_MIC

// HERE : try float version of FRA 
/* FRA() corresponds to FR() with location of H[U] or H[V] modified by added site to HU and HV respectively */
template<int ERFC>
static inline double FRA(int I,
			 int K,
			 int J,
			 int U,
			 int V,
			 double y,/* value Yc(H,I,K) */
			 double *D,/* Data map */
			 double *H, /* Reference map */
			 double HV, /* Reference map location H[V] (as modified by new added site location) */
			 double HU, /* Reference map location H[U] (as modified by new added site location) */
			 int M, /* number of sites in Data map */
			 int N, /* number of sites in Reference map */
			 int rc) /* If right end of H is a chromosome end */
{
  if(DEBUG>=2) assert(M-J <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX);
  if(DEBUG>=2) assert(V-I <= N);
  if(DEBUG>=2) assert(END_FIX || rc);
  if(DEBUG>=2) assert(HV >= H[I]);

  double x = D[M+1] - D[J];
  double Bias = x * FrateBiasWT + (M-J) * FBiasWT + EBiasWT;
  if(DEBUG>=2 && biasWT <= 0.0) assert(Bias <= 0.0);

  if(END_FIX>=2 && HV <= H[I]){
    double logerr = -x*Frate + Bias + LogFpPow[M-J];
    logerr = min(logerr, (double)logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0
    double Vexp = exp(logerr);
    if(DEBUG>=2) assert(isfinite(Vexp));

    if(DEBUG>=2) assert(y <= HU);
    return Vexp * Rtheta * GE<ERFC>(D[M+1]-D[J],HU - y, HU - y, (EXTEND_MAX>0 && !rc && U > N) ? 1 : 0);
  }

  double logerr = -x*Frate + Bias + LogFpPow[M-J];
  logerr = LogFnPow[V-I] + min(logerr, (double)logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0
  double Vexp = exp(logerr);
  if(DEBUG>=2) assert(isfinite(Vexp));
  if(DEBUG>=2) assert(HV <= HU);
  if(DEBUG>=2) assert(0 <= HU-HV && HU-HV <= HU-y);
  double GEval = GE<ERFC>(D[M+1]-D[J],HU - y, HU - HV, (EXTEND_MAX>0 && !rc && U > N) ? 1 : 0);

  if(DEBUG>=2){
    double val = (Vexp * Rtheta) * GEval;
    if(!(isfinite(val) && val >= 0.0)){
      printf("FRA(I=%d,K=%d,J=%d,U=%d,V=%d,y=%0.8f,HV=%0.8f,HU=%0.8f,M=%d,N=%d,rc=%d):val= %0.8e (Vexp= %0.8e, Rtheta= %0.8e, GEval= %0.8e)\n",
	     I,K,J,U,V,y,HV,HU,M,N,rc,val,Vexp,Rtheta,GEval);
      fflush(stdout);
      assert(isfinite(val) && val >= 0.0);
    }
  }

  return (Vexp * Rtheta) * GEval;
}

#pragma GCC_pop_options


/* FRE() corresponds to the end alignment gaussian term G(D[M+1]-D[J],H[U]-y) of FR() (U==N+1) */
static inline double FRE(int I,
			 int K,
			 int J,
			 int U,
			 int V,
			 double y,/* value of Yc(H,I,K) */
			 double *D,/* Data map */
			 double *H,/* reference map */
			 int M) /* sites in Data map */
{
  if(DEBUG>=2) assert(U==gN+1);
  if(DEBUG>=2) assert(M-J <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX);
  if(DEBUG>=2) assert(V-I <= gN);

  double x = D[M+1] - D[J];
  double Bias = x * FrateBiasWT + (M-J) * FBiasWT + EBiasWT;
  if(DEBUG>=2 && biasWT <= 0.0) assert(Bias <= 0.0);

  double logerr = -x*Frate + Bias + LogFpPow[M-J];
  logerr = LogFnPow[V-I] + min(logerr, (double)logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0
  double Vexp = exp(logerr);
  if(DEBUG>=2) assert(isfinite(Vexp));

  return Vexp * G(x, H[U]-y) * InvRtheta;
}

/* faster float version of FRE() */
static inline MFLOAT FRE(int I,
			 int K,
			 int J,
			 int U,
			 int V,
			 MFLOAT x,/* D[M+1] - D[J] */
			 MFLOAT Yui,/* H[U] - Yc(H,I,K) */
			 int M) /* sites in Data map */
{
  if(DEBUG>=2) assert(U==gN+1);
  if(DEBUG>=2) assert(M-J <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX);
  if(DEBUG>=2) assert(V-I <= gN);

  MFLOAT Bias = x * FrateBiasWT + FBiasWTm[M-J] + EBiasWT;
  if(DEBUG>=2 && biasWT <= 0.0) assert(Bias <= 0.0);

  MFLOAT logerr = -x*Frate + Bias + LogFpPow[M-J];
  logerr = LogFnPow[V-I] + min(logerr, logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0
  MFLOAT Vexp = expF(logerr);
  if(DEBUG>=2) assert(isfinite(Vexp));

  return Vexp * InvRthetaF * G(x, Yui);
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT

/* AVX256 float version of FRE() */
static inline __m256 FRE_mm256(int I,
				 int K,
				 int J,
				 int U,
				 int V,
				 __m256 x,/* D[M+1] - D[J] */
				 __m256 Yui,/* H[U] - Yc(H,I,K) */
				 int M) /* sites in Data map */
{
  __m256 Bias = x * v256_FrateBiasWT + _mm256_set1_ps(FBiasWTm[M-J] + EBiasWT);

  __m256 logerr = Bias + _mm256_set1_ps(LogFpPow[M-J]) - x * v256_Frate;
  logerr = _mm256_set1_ps(LogFnPow[V-I]) + _mm256_min_ps(logerr, v256_logMaxUR);// avoid overflow for large end intervals, especially with -biasWT 0
  __m256 Vexp = expF_mm256(logerr);

  return Vexp * v256_InvRthetaF * G_mm256(x, Yui);
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

/* AVX512 float version of FRE() */
static inline __m512 FRE_mm512(int I,
				 int K,
				 int J,
				 int U,
				 int V,
				 __m512 x,/* D[M+1] - D[J] */
				 __m512 Yui,/* H[U] - Yc(H,I,K) */
				 int M) /* sites in Data map */
{
  __m512 Bias = _mm512_add_ps(_mm512_mul_ps(x , v512_FrateBiasWT), _mm512_set1_ps(FBiasWTm[M-J] + EBiasWT));

  __m512 logerr = _mm512_sub_ps(_mm512_add_ps(Bias , _mm512_set1_ps(LogFpPow[M-J])), _mm512_mul_ps(x, v512_Frate));
  logerr = _mm512_add_ps(_mm512_set1_ps(LogFnPow[V-I]),_mm512_min_ps(logerr, v512_logMaxUR));// avoid overflow for large end intervals, especially with -biasWT 0
  __m512 Vexp = expF_mm512(logerr);

  return _mm512_mul_ps(_mm512_mul_ps(Vexp, v512_InvRthetaF), G_mm512(x, Yui));
}

#endif // VECTOR_AVX || VECTOR_MIC

static __attribute__ ((noinline)) double PM(int J, int I, int K, double *Y)
{
  double p = pTP;
  if(DEBUG>=2) assert(0.0 <= p && p <= 1.0);
  if(K <= 0)
    return p;

  if(DEBUG>=2) assert(1 <= K && K <= KMAX+1 && KMAX <= K_MAX);
  if(DEBUG>=2) assert(SCORE_APPROX <= 1);

  // NEW CODE : only K(K+1)/2 calls to Pn() instead of (K+1)/2 * 2^(K-1) calls AND always exact
  double q = 1.0 - p;
  // WAS93  double PMK[K+1], PowQ[K+2]; // triggers compiler bug in GCC-7.4 and GCC-8.3
  double PMK[K_MAX + 1],PowQ[K_MAX + 2];
  PowQ[0] = 1.0;

  for(int k = 0; k < K; k++){
    PowQ[k+1] = q * PowQ[k];
    /* PowQ[t = 0..k+1] = pow(q,t) */
    double Yk = Y[I-k-1];
    double sum = Pn(Y[I]-Yk) * PowQ[k];
    if(DEBUG>=2) assert(isfinite(sum) && sum >= 0.0 && sum <= 1.0);
    for(int t = 1; t <= k; t++)
      sum += Pn(Y[I-t] - Yk) * PowQ[k-t] * p * PMK[t-1];
    PMK[k] = sum;

    if(DEBUG>=2 && !(isfinite(PMK[k]) && PMK[k] >= 0.0 && PMK[k] <= 1.0)){
      //      #pragma omp critical
      {
	printf("\nI=%d,K=%d, Y[I]= %0.8f, Y[I-K]= %0.8f, p= %0.8f, q= %0.8f: k=%d, PMK[k]= %0.10e\n",I,K,Y[I],Y[I-K],p,q,k,PMK[k]);
	fflush(stdout);

	for(int t = 0; t < k; t++){
	  printf("\t t=%d:PMK[t]= %0.10e(%p), PowQ[t]= %0.10e(%p)\n", t, PMK[t],&PMK[t],PowQ[t],&PowQ[t]);
	  fflush(stdout);
	}

	double mysum = Pn(Y[I]-Yk) * PowQ[k];
	printf("\t Y[I-k-1]= %0.8f, Yk= %0.8f, Pn(Y[I]-Yk)= %0.10e, PowQ[k]= %0.10e(%p): sum= %0.10e\n",Y[I-k-1], Yk, Pn(Y[I]-Yk), PowQ[k], &PowQ[k], mysum);
	fflush(stdout);

        for(int t = 1; t <= k; t++){
	  mysum += Pn(Y[I-t] - Yk) * PowQ[k-t] * p * PMK[t-1];
	  printf("\t   t=%d: Y[I-t]= %0.8f, Pn(Y[I-t] - Yk)= %0.10e, PowQ[k-t]= %0.10e(%p), PMK[t-1]= %0.10e(%p) : sum -> %0.10e\n",
		 t,Y[I-t],Pn(Y[I-t] - Yk),PowQ[k-t],&PowQ[k-1],PMK[t-1],&PMK[t-1],mysum);
	  fflush(stdout);
	}
	assert(isfinite(PMK[k]) && PMK[k] >= 0.0 && PMK[k] <= 1.0);
      }
    }
  }

  if(DEBUG>=2 && !(isfinite(PMK[K-1]) && PMK[K-1] >= 0.0 && PMK[K-1] <= 1.0)){
    #pragma omp critical
    {
      int k = K-1;
      double Yk = Y[I-k-1];

      printf("\nI=%d,K=%d,Y[I]=%0.8f,Y[I-K]=%0.8f,p= %0.8f, q= %0.8f: k=%d, PMK[k]= %0.10e(%p)\n",I,K,Y[I],Y[I-K],p,q,k,PMK[k],&PMK[k]);
      fflush(stdout);

      for(int t = 0; t < k; t++){
	printf("\t t=%d:PMK[t]= %0.10e(%p), PowQ[t]= %0.10e(%p)\n", t, PMK[t],&PMK[t],PowQ[t],&PowQ[t]);
	fflush(stdout);
      }

      double mysum = Pn(Y[I]-Yk) * PowQ[k];
      printf("\t Y[I-k-1]= %0.8f, Yk= %0.8f, Pn(Y[I]-Yk)= %0.10e, PowQ[k]= %0.10e(%p): sum= %0.10e\n",Y[I-k-1], Yk, Pn(Y[I]-Yk), PowQ[k], &PowQ[k], mysum);
      fflush(stdout);
    
      for(int t = 1; t <= k; t++){
	mysum += Pn(Y[I-t] - Yk) * PowQ[k-t] * p * PMK[t-1];
	printf("\t   t=%d: Y[I-t]= %0.8f, Pn(Y[I-t] - Yk)= %0.10e, PowQ[k-t]= %0.10e(%p), PMK[t-1]= %0.10e(%p) : sum -> %0.10e\n",
	       t,Y[I-t],Pn(Y[I-t] - Yk),PowQ[k-t],&PowQ[k-1],PMK[t-1],&PMK[t-1],mysum);
	fflush(stdout);
      }

      assert(isfinite(PMK[K-1]) && PMK[K-1] >= 0.0 && PMK[K-1] <= 1.0);    
    }
  }

  if(PRbiasWT*biasWT > 0.0){
    double cum = 0.0;
    for(int k = I-K+1; k <= I; k++){
      double y = Y[k] - Y[k-1];
      double pr = Pr(y);
      if(DEBUG>=2) assert(isfinite(pr) && 0.0 <= pr && pr <= 1.0);
      double prbias = (pr <= 0.0 || pr >= 1.0) ? 0.0 : pr * log(pr) + (1.0 - pr) * log(1.0 - pr);
      if(DEBUG>=2) assert(isfinite(prbias) && prbias <= 0.0);
      cum -= prbias;
      if(DEBUG>=2) assert(isfinite(cum) && cum >= 0.0);
    }

    return p * p * PMK[K-1] * exp(cum * PRbiasWT * biasWT);
  }

  return p * p * PMK[K-1];
}

static __attribute__ ((noinline)) double PM2(int J, int I, int K, double *Y)// with more debug code than PM()
{
  double p = pTP;
  if(DEBUG>=2) assert(0.0 <= p && p <= 1.0);
  if(K <= 0)
    return p;

  if(DEBUG) assert(SCORE_APPROX <= 1);

  // NEW CODE : only K(K+1)/2 calls to Pn() instead of (K+1)/2 * 2^(K-1) calls AND always exact
  double q = 1.0 - p;
  double PMK[K+1], PowQ[K+2];
  if(DEBUG>=2) assert(1 <= K && K <= KMAX+1);
  PowQ[0] = 1.0;
  for(int k = 0; k < K; k++){
    PowQ[k+1] = q * PowQ[k];
    /* PowQ[t = 0..k+1] = pow(q,t) */
    double Yk = Y[I-k-1];
    double sum = Pn(Y[I]-Yk) * PowQ[k];
    if(DEBUG>=2) assert(isfinite(sum) && sum >= 0.0 && sum <= 1.0);
    for(int t = 1; t <= k; t++)
      sum += Pn(Y[I-t] - Yk) * PowQ[k-t] * p * PMK[t-1];
    PMK[k] = sum;

    if(DEBUG && !(isfinite(PMK[k]) && PMK[k] >= 0.0 && PMK[k] <= 1.0)){
      //      #pragma omp critical
      {
	printf("\nI=%d,K=%d, Y[I]= %0.8f, Y[I-K]= %0.8f, p= %0.8f, q= %0.8f: k=%d, PMK[k]= %0.10e\n",I,K,Y[I],Y[I-K],p,q,k,PMK[k]);
	fflush(stdout);

	for(int t = 0; t < k; t++){
	  printf("\t t=%d:PMK[t]= %0.10e, PowQ[t]= %0.10e\n", t, PMK[t],PowQ[t]);
	  fflush(stdout);
	}

	double mysum = Pn(Y[I]-Yk) * PowQ[k];
	printf("\t Y[I-k-1]= %0.8f, Pn(Y[I]-Yk)= %0.10e, PowQ[k]= %0.10e: sum= %0.10e\n",Y[I-k-1], Pn(Y[I]-Yk), PowQ[k], mysum);
	fflush(stdout);

        for(int t = 1; t <= k; t++){
	  mysum += Pn(Y[I-t] - Yk) * PowQ[k-t] * p * PMK[t-1];
	  printf("\t   t=%d: Y[I-t]= %0.8f, Pn(Y[I-t] - Yk)= %0.10e, PowQ[k-t]= %0.10e, PMK[t-1]= %0.10e : sum -> %0.10e\n",
		 t,Y[I-t],Pn(Y[I-t] - Yk),PowQ[k-t],PMK[t-1],mysum);
	  fflush(stdout);
	}
	assert(isfinite(PMK[k]) && PMK[k] >= 0.0 && PMK[k] <= 1.0);
      }
    }
  }
  if(PRbiasWT*biasWT > 0.0){
    double cum = 0.0;
    for(int k = I-K+1; k <= I; k++){
      double y = Y[k] - Y[k-1];
      double pr = Pr(y);
      if(DEBUG>=2) assert(isfinite(pr) && 0.0 <= pr && pr <= 1.0);
      double prbias = (pr <= 0.0 || pr >= 1.0) ? 0.0 : pr * log(pr) + (1.0 - pr) * log(1.0 - pr);
      if(DEBUG>=2) assert(isfinite(prbias) && prbias <= 0.0);
      cum -= prbias;
      if(DEBUG>=2) assert(isfinite(cum) && cum >= 0.0);
    }

    return p * p * PMK[K-1] * exp(cum * PRbiasWT * biasWT);
  }

  return p * p * PMK[K-1];
}

} // namespace probeval

#endif
