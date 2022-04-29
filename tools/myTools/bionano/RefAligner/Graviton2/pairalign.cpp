#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#ifdef _MSC_VER
#define strtoll _strtoi64
#include <malloc.h>
#define alloca _alloca 
#else
#include <alloca.h>
#endif

#ifdef VALGRIND
#include <valgrind/memcheck.h>
#endif

#ifdef WIN32
#define copysign _copysign
#endif

#ifdef CLANG
#include <sanitizer/msan_interface.h>
#endif

#include <sys/mman.h>
// #include <omp.h>

#include "globals.h"
#include "parameters.h"
#include "Calign.h"
#include "hash.h"
#include "Ccontig.h"

#ifdef __SSE3__
#include <pmmintrin.h>
#endif

#ifdef  __INTEL_COMPILER
#define CAST (const char *)
#else
#ifdef WIN32
#define CAST (const char *)
#else
#define CAST
#endif
#endif

#ifdef WIN32
static double denom = 0.0;
static PFLOAT aNaN = 0.0/denom;
#else
static PFLOAT aNaN = nan("NaN");
#endif

//#undef DEBUG
//#define DEBUG 2

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/pairalign.cpp 11549 2020-08-28 17:04:37Z tanantharaman $");

#define DIAGONAL 1 /* HERE HERE TRY 2 */ /* > 0 : with hashtable offset, diagonalize alignment array to increase speed and (if >= 2) reduce real memory usage */
#define MIN_MEM 0 /* RA_MIN_MEM */ /* minimize virtual memory usage by possibly reallocating memory in each call of pairalignXY() */
#define MIN_MEM_DEBUG 0

#define ENDOUTLIER_FIX 1 /* WAS158 0 */ /* update A.G to -1 when A.score is updated to non-endoutlier value (> C) : some slowdown, 
					   but required to avoid misclassifying alignment left ends as endoutliers */

#define IDMAX 1000000000000000LL /* 1e15 which is around MAXINT64/9223) */

#define MAX_ALLOCA 128 /* maximum doubles to allocate on stack, to avoid stack overflow */

#define THREAD_DEBUG 0 /* count number of calls to madvise(),getmem(),getstatm(),ReleaseMemory(),WaitForMemory() and OMP critical and OMP atomic */

#define TIME_SPREAD 0 /* TRY 1 */ /* spread out calls to GetStatm() and madvise() so threads are not synchronized : slightly slower ! */

#define PAIRMERGE_SORTED 0 /* sort pairmerge contigs from largest to smallest : simplifies debugging but reverses effect of -randomize */

#define OUTLIEREXTEND_FIX 1 /* with -outlierExtend : Try to MERGE current interval with a previous outlier interval : slightly more merges, performance slowdown */
#define PAIRMERGE_PARALLEL_FIX 0 /* WAS158 1 */ /* apply fixed algorithm with defered handling of outlierExtend with G==I */

#define VMEMCHECK 0 /* keep track of vmemTotsiz */
#define MEMCHECK 99  /* Call WaitForMemory() at start of each pairalignXY(), instead of at start of pairalign() (Only applied if RA_MIN_MEM <= 7)
		       If >= 2 : Also recheck after initializing arrays in pairalignXY (up to MEMCHECK-1 times), subject to VMRSS_CHECK */
#define MEMRESERVE 1 /* reserve VmRSS memory using global variable VmRSSreserved */

#define VMRSS_SLEEP 0.050 /*  WAS 1.0 */ /* Sleep Paused thread for VMRSS_SLEEP seconds before re-checking VmRSS etc */
#define VMRSS_CHECK 0.020 /* Update VmRSS no more than once every VMRSS_CHECK seconds (in each thread) by checking /proc/self/statm (but also check at thread peak memory if nMemSiz > VMRSS_PEAKCHECK) */
#define VMRSS_CHECK2 2.000 /* WAS114 0.200 */ /* Update VmRSS no more than once every VMRSS_CHECK2 seconds (in any thread) while VmRSS+VmSwap <= 0.9 * RA_MAX_MEM * maxmem (see VMRSS_CHECK otherwise) */
#define VMRSS_PEAKCHECK 0 /* 2000000000LL */ /* see VMRSS_CHECK */
#define VMSWAP_CHECK 60.000 /* Update VmSwap no more than once every VMSWAP_CHECK seconds (in each thread) by checking /proc/self/status */

#define PAIRMERGE_SEGDUP_ONESIDED SplitSegDup_OneSided /* Check for SegDup alignments with only one side having an endoutlier over E kb ( & EN labels) : see -SplitSegDup
							  If >= 2: The end without the endoutlier cannot be a truncated end (due to extend & split) : only applies if endoutlier < SplitSegDup_OneSided_MaxLen */

#define PAIRMERGE_REPEATMASK 1 /* 1 : suppress repeat masking if PairMerge and one map completely overlaps the other OR (SplitSegDup > 0 and one endoutlier is > SplitSegDupE) */

#define PAIRMERGE_FIX 1 /* ignore alignments that violate -pairmerge limits on size of internal and end outliers (in case other orientation has one but with lower score)
			   NOTE1 : When -pairmergeHmap is active, internal outlier limit is not checked.
			   NOTE2 : When -SplitSegDup is active, the endoutlier limit is not checked and internal outlier limits have been adjusted (typically to a smaller value)

			   >= 1 : NEW24 : Partly enabled : Never checks endoutlier limits  : limits this option to violations of internal outlier limits (disables completely when -pairmergeHmap is active)
			   >= 2 : NEW24 : Fully enabled.
			*/

#define PAIRMERGE_OUTLIEREXTEND_FIX 1 /* apply 2nd call to pairalignXY() with outlierExtend recursively (required if PAIRMERGE_FIX) */

#define PAIRMERGE_SEGDUP_OVERLAPPED_FIX 1 /* NEW24 : allow merge of fully overlapped maps even when -SegDupSplit is active */

#define TRUNCATE_FIX TruncateFix /* Don't truncate map at one end if it is the longer map at the other end */

#define OUTLIER_FIXES PairwiseOutlierFix /* enable various fixes to scoring of outliers (see command line parameter -PairwiseOutlierFix, on by default) */

#include "NGentigPairScore.h"

// This scoring function works poorly with very low FP/FN, eg FP=0.2, FN=0.04, sf=0, sd=0.04 and -outlier 0.0001 (see E.Coli alignments pair2.fixx) : find out why (possible error in seperating Pen and Bias, or
//    use of fixed K values independent of FP rate)

#define MAX_ALLOCA 128 /* maximum doubles to allocate on stack, to avoid stack overflow */

#define PAIRSPLIT_OVERLAP 1 /* WAS 0 */ /* include the internal outliers boundary sites in the extension region (with -pairsplitExtend) */
#define PAIRSPLIT_EXTEND PairsplitExtend /* extend center cloned regions by -ve scoring region, to handle internal outliers that have been split into pieces and may not have the right boundaries */

#define FIX_CONSTRAINT 1 /* TRY 1,2 */ /* > 0 : more conservative RepeatRec constraint */

#define HASH_DEBUG 0 /* To debug HASH_STREAM with verbose buffer update messages */

#define SEND_THRESH 1 /* enable use of EndScoreThresh  */

#define DEFER_BP 1 /* defer computation of DP backpointer to allow recurrance to be faster (no branch penalty, just a max operator) */
#define DEFER_G 2 /* >= 1 : defer computation of A.G[I][J] == -1
		     >= 2 : defer update of corresponding A.Lijx[I][J], A.Lijy[I][J] to reflect Sbnd() */

#define ENDFIX 2 /* allow for sizing error at ends (Sbnd) (1 = allow for sizing error at right end, 2 = allow for sizing errors at both ends) */
#define ENDINIT 1 /* limit end initialization to cases where X OR Y is within DELTA sites of end */

#define LOCAL 1  /* allow local alignments (ignore left/right end alignment score if it is worse than ChimScore) */


#define MAXDELTA 2.0 /* maximum overlap error as multiple of Lambda to be classified as TP */

#define MINDELTA 20.0 /* minimum overlap error to be classified as FP */

#define MISALIGN_BINS 25 /* number of points on misaligned cuts vs aligned interval size graph (fewer points will produce more reliable Y-axis values, but less X-axis resolution) */

#define MINSITES AlignedSiteThreshold /* minimum sites for a chimeric fragment : smaller fragments are ignored since it is not possible to align them reliably */

#define PVERB 0 /* Display detailed alignment information for all alignments above threshold */

#define TRACE 0 /* Display details of all alignments, even those NOT above threshold */
#define REF_ID 119LL
#define QRY_ID 100037LL
#define QRY_ID2 11LL 
#define QRY_ID3 -1LL

#define EVERB ((TRACE && ((YYmap[Yid]->id == REF_ID && XXmap[Xid]->id == QRY_ID) || (YYmap[Yid]->id == QRY_ID && XXmap[Xid]->id == REF_ID) || \
                          (YYmap[Yid]->id == REF_ID && XXmap[Xid]->id == QRY_ID2) || (YYmap[Yid]->id == QRY_ID2 && XXmap[Xid]->id == REF_ID)  || \
                          (YYmap[Yid]->id == REF_ID && XXmap[Xid]->id == QRY_ID3) || (YYmap[Yid]->id == QRY_ID3 && XXmap[Xid]->id == REF_ID) )) ? 1 : 0)


#define I_TRACE -1 // 2534
#define J_TRACE -1 // 1246
#define G_TRACE -1
#define H_TRACE -1

#define LVERB 0 /* Display details of multiple local alignment (implies PVERB==1) */

#if USE_MIC
#define memcpy _mm512_memcpy
#endif

extern double mtime();/* microsecond timer function : see refine.cpp */
extern double wtime();

Cmap **YYmap = 0,**XXmap= 0;
int numY= 0,numX= 0;

static long long IDmax = 0;/* If PairMerge : largest ID value of original ids */

static int repeatcnt = 0;/* global count of alignments that were supressed due to -RepeatMask */

/* sort alignments in order of id2, id1, orientation, sites1[0], sites1[U-1],sites2[0],sites2[U-1] (null alignments last) */
static int IdOrSiteInc(register Calign **p1, register Calign **p2)
{
  Calign *align1 = p1[0];
  Calign *align2 = p2[0];
  if(!align1) return 1;
  if(!align2) return -1;

  Cmap *Ymap1 = YYmap[align1->mapid1];
  Cmap *Ymap2 = YYmap[align2->mapid1];
  Cmap *Xmap1 = XXmap[align1->mapid2];
  Cmap *Xmap2 = XXmap[align2->mapid2];

  int ret = (Xmap1->id > Xmap2->id) ? 1 : (Xmap1->id < Xmap2->id) ? -1 :
    (Ymap1->id > Ymap2->id) ? 1 : (Ymap1->id < Ymap2->id) ? -1 :
    (align1->orientation > align2->orientation) ? 1 : (align1->orientation < align2->orientation) ? -1 : 0;
  if(ret)
    return ret;

  int left1 = align1->sites1 ? (align1->sites1[0] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int left2 = align2->sites1 ? (align2->sites1[0] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  ret = (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  int numpairs1 = align1->numpairs;
  int numpairs2 = align2->numpairs;
  int right1 = align1->sites1 ? (align1->sites1[numpairs1-1] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int right2 = align2->sites1 ? (align2->sites1[numpairs2-1] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  ret = (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
  if(ret)
    return ret;

  int M1 = Xmap1->numsite[0];
  int M2 = Xmap2->numsite[0];
  int Xshift1 = Xmap1->origmap ? max(1,Xmap1->left[0]) - 1 : 0;
  int Xshift2 = Xmap2->origmap ? max(1,Xmap2->left[0]) - 1 : 0;
  int LJ1 = align1->sites2[0], LJ2 = align2->sites2[0];
  int RJ1 = align1->sites2[numpairs1-1], RJ2 = align2->sites2[numpairs2-1];
  left1 = align1->sites2==NULL ? 0 : align1->orientation ? Xshift1 + M1 + 1 - RJ1 : Xshift1 + LJ1;
  left2 = align2->sites2==NULL ? 0 : align2->orientation ? Xshift2 + M2 + 1 - RJ2 : Xshift2 + LJ2;
  ret = (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;
  
  right1 = align1->sites2==NULL ? 0 : align1->orientation ? Xshift1 + M1 + 1 - LJ1 : Xshift1 + RJ1;
  right2 = align2->sites2==NULL ? 0 : align2->orientation ? Xshift2 + M2 + 1 - LJ2 : Xshift2 + RJ2;
  ret = (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
  return ret;

}

static double ChimScore = -1.0+100;

typedef CHashMatch Chashpair;
static Chashpair *hashpairs1 = 0, *hashpairs2 = 0;
static Chashpair *nexthash1 = 0, *nexthash2 = 0;/* nexthash points into hashpairs[0..numhashpairs-1] to the next unused entry */
static size_t maxhashpairs1 = 0, maxhashpairs2 = 0;
static size_t numhashpairs1 = 0, numhashpairs2 = 0;

static int *YidList1 = 0, *YidList2 = 0;
static int *XidList1 = 0, *XidList2 = 0;
static Chashpair **phashList1 = 0, **phashList2 = 0;
static int NumPair1 = 0, NumPair2 = 0, MaxPair = 0;

static int BinaryFormat = 0;
static FILE *HashFP = NULL;
static size_t MatchMax = 0;
static size_t totalhashpairs = 0;

static long long *Ymapid2id = 0;
static long long *Xmapid2id = 0;

static double AvgInterval = 0.0;

/** aligned interval statistics */
class Cinterval {
public:
  double x,y;
};

/** misaligned sites distribution by aligned interval size */
class Cmisalign {
public:
  double x;/* size of aligned interval (for each map seperately) */
  int m;/* number of misaligned sites */
};

static int xinc(Cmisalign *p1, Cmisalign *p2)
{
  return (p1->x > p2->x) ? 1 : (p1->x < p2->x) ? -1 : 0;
}

static int CmapIdInc(Cmap **p1, Cmap **p2)
{
  /* Cannot use p1[0]->id - p2[0]->id since that may be out of range for integer */
  return (p1[0]->id > p2[0]->id) ? 1 : (p1[0]->id < p2[0]->id) ? -1 : 0;
}

static int CmapLenDec(Cmap **p1, Cmap **p2)
{
  Cmap *pmap1 = *p1;
  Cmap *pmap2 = *p2;
  int N1 = pmap1->numsite[0];
  int N2 = pmap2->numsite[0];
  int Len1 = pmap1->site[0][N1+1];
  int Len2 = pmap2->site[0][N2+1];

  return (Len1 < Len2) ? 1 : (Len1 > Len2) ? -1 : 0;
}

/** Dynamic programming array cell : only used for per cell sizing */
class DP {
public:
  PFLOAT score;/**< alignment score for best alignment between Y[0..I] & X[0..J] with Y[I] aligned with X[J] */
  PFLOAT Uscore;/**< score with nothing except (I,J) aligned (G < 0, H < 0) */
  int G;/**< previous site in Y in best alignment (< 0 if none) */
  int H;/**< previous site in X in best alignment (undefined if G < 0) */
  int Rijy;/**< rightmost site/end in Y that overlaps X when Y[I] is aligned with X[J] */
  int Lijy;/**< leftmost site/end in Y that overlaps X when Y[I] is aligned with X[J] */
  int Rijx;/**< rightmost site/end in X that overlaps Y when Y[I] is aligned with X[J] */
  int Lijx;/**< leftmost site/end in X that overlaps Y when Y[I] is aligned with X[J] */
  /* If ENDINIT >= 1:
     Lijx,Lijy are only defined if (I <= DELTA || J <= DELTA)
     Rijx,Rijy are only defined if (I >= N-DELTA || J >= M-DELTA) */
};

static long long DPsizF = 2;/* number of PFLOAT fields in DP */
static long long DPsizI = 6;/* number of int fields in DP */

/** Array based Dynamic programming data */
class AL {
public:
  PFLOAT **score;
  PFLOAT **Uscore;
  int **G;
  int **H;
  int **Rijy;
  int **Lijy;
  int **Rijx;
  int **Lijx;

  int *JMIN, *JMAX;
  int *Imin,*Imax;

  long long M_stride;
  long long array_stride;
  //  inline int G(int I, int J) { return (&G[1][-M_stride])[M_stride*I + J]; };
};
#define FLAT_A(ptr, I, J)	(ptr[1][(I-1)*A.M_stride + J])
#define FLAT_DELX(ptr, n, J, M)  (ptr[1][J + (n-1) * COMPUTE_STRIDE4(M-1)])

#define PAIRHASH_BITS 20

static class CPairHash {
public:
  int Xid;
  int Yid;
  CPairHash *next;
} **PairHash = 0;/* PairHash[0..(1<<PAIRHASH_BITS)-1] */

static void PairHash_init()
{
  if(DEBUG) assert(PAIRHASH_BITS <= 31);
  if(!PairHash){
    PairHash = new CPairHash*[1<<PAIRHASH_BITS];
    for(int i = 1<<PAIRHASH_BITS; --i >= 0;)
      PairHash[i] = 0;
  }
}

static void PairHash_free()
{
  if(PairHash){
    for(int i = 1<<PAIRHASH_BITS; --i >= 0;){
      for(CPairHash *q,*p = PairHash[i];p;p = q){
	q = p->next;
	delete p;
      }
    }
    delete [] PairHash;
    PairHash = 0;
  }
}

static inline int PairHash_hashmap(int Yid, int Xid)
{
  return ((Yid<<8) ^ (Yid >> (PAIRHASH_BITS-8)) ^ Xid) & ((1<<PAIRHASH_BITS)-1);
}

static inline void PairHash_insert(int Yid, int Xid)
{
  CPairHash **qq = &PairHash[PairHash_hashmap(Yid,Xid)];

  for(CPairHash *p = *qq; p; p = p->next)
    if(p->Yid == Yid && p->Xid == Xid)
      return;/* already in hashtable */

  /* insert new entry at *qq  */
  CPairHash *q = new CPairHash;
  q->Yid = Yid;
  q->Xid = Xid;
  q->next = *qq;
  *qq = q;
}

static inline int PairHash_find(int Yid, int Xid)
{
  CPairHash **qq = &PairHash[PairHash_hashmap(Yid,Xid)];

  for(CPairHash *p = *qq; p; p = p->next)
    if(p->Yid == Yid && p->Xid == Xid)
      return 1;/* found in hashtable */

  return 0;
}

/** take existing alignment and return pvalue logFP() */
static double alignFP(Calign *p, FLOAT *Y, FLOAT *X, int N, int M,
		      int orientation, int scaleID, int Yid, int Xid, double lscore, int verb)
{
  int U = p->numpairs;
  
  /* allocate alignment fragment sizes */
  FLOAT *Xfrag, *Yfrag;
  if(U > MAX_ALLOCA){
    Xfrag = new FLOAT[2*(U-1)];
    Yfrag = &Xfrag[U-1];
  } else {
    Xfrag = (FLOAT *) alloca((U-1)*sizeof(FLOAT));
    Yfrag = (FLOAT *) alloca((U-1)*sizeof(FLOAT));
  }
  int misscnt = 0;
  int falsecnt = 0;
  int OutlierCnt = 0;/* number of internal outliers */
  double OutlierMN = 1.0;/* product over all outlier of m x n (the number of intervals in the internal outlier aligned segment) */

  int I= -1, J= -1, G= -1,H= -1;
  int F = 0;/* segment (aligned interval) count */

  for(int T = 0; T < U; T++, G = I, H = J){
    I = p->sites1[T];
    J = p->sites2[T];
    if(T > 0){
      if(DEBUG) assert(G > 0 && H > 0);
      if(p->iscore[T] > p->outscore[T] + Poutlier + (PFLOAT)0.01 || J <= H || I <= G){
	if(LVERB>=2 ||verb/* HERE HERE */)
	  printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:iscore[T]=%0.6f,outscore[T]=%0.6f\n",T,U,I,J,G,H,p->iscore[T],p->outscore[T]);
	if(OUTLIER_PV){
	  OutlierCnt++;
	  int m = max(1,J-H);
	  int n = max(1,I-G);
	  if(DEBUG>=2 && !outlierExtend && !(J > H && I > G)){
	    printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:iscore[T]=%0.6f,outscore[T]=%0.6f\n",T,U,I,J,G,H,p->iscore[T],p->outscore[T]);
	    fflush(stdout);
	    assert(J>H);
	    assert(I>G);
	  }
	  if(DEBUG>=2 && !(OutlierMN*m*n > 0)){
	    printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:iscore[T]=%0.6f,outscore[T]=%0.6f:m=%d,n=%d,OutlierMN=%0.1f -> %0.1f,OutlierCnt=%d\n",
		   T,U,I,J,G,H,p->iscore[T],p->outscore[T],m,n,OutlierMN,OutlierMN*m*n,OutlierCnt);
	    fflush(stdout);
	    assert(OutlierMN*m*n > 0);
	  }
	  OutlierMN *= (double)(m*n);
	  if(DEBUG>=2 && !(OutlierMN>0)){
	    printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:iscore[T]=%0.6f,outscore[T]=%0.6f:m=%d,n=%d,OutlierMN=%0.1f\n",
		   T,U,I,J,G,H,p->iscore[T],p->outscore[T],m,n,OutlierMN);
	    fflush(stdout);
	    assert(OutlierMN > 0);
	  }
	}
	continue;
      }
      if(LVERB>=2 ||verb/* HERE HERE */)
	printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:F=%d:y=%0.3f,x=%0.3f,misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	       T,U,I,J,G,H,F,Y[I]-Y[G],X[J]-X[H],misscnt,misscnt+(I-G-1),falsecnt,falsecnt + (J-H-1),OutlierCnt,OutlierMN);
      Yfrag[F] = Y[I] - Y[G];
      Xfrag[F] = X[J] - X[H];
      if(DEBUG>=1+RELEASE && !(Yfrag[F] > 0.0 && Xfrag[F] > 0.0 && I > G && J > H)){
	printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:F=%d:y=%0.3f,x=%0.3f,misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	       T,U,I,J,G,H,F,Y[I]-Y[G],X[J]-X[H],misscnt,misscnt+(I-G-1),falsecnt,falsecnt + (J-H-1),OutlierCnt,OutlierMN);
	fflush(stdout);
	assert(Yfrag[F] > 0.0);
	assert(Xfrag[F] > 0.0);
	assert(I > G);
	assert(J > H);
      }
      F++;
      misscnt += I-G-1;
      falsecnt += J-H-1;
    } else {/* left end */
      int Lijy = p->Lij1;
      int Lijx = p->Lij2;
      if(DEBUG>=1+RELEASE && !(Lijy >= 0 && Lijy <= I)){
        #pragma omp critical
        {
	  printf("alignFP:Yid=%d,Xid=%d,or=%d,sc=%d:N=%d,M=%d,score=%0.4f:U=%d,I=%d,J=%d,Lijy=%d,Lijx=%d,Lend=%d,misscnt=%d,falsecnt=%d\n",
		 Yid,Xid,orientation,scaleID,N,M,lscore,U,I,J,Lijy,Lijx,p->Lend,misscnt,falsecnt);
	  fflush(stdout);
	  assert(Lijy >= 0 && Lijy <= I);
	}
      }
      if(DEBUG>=1+RELEASE) 	assert(Lijx >= 0 && Lijx <= J);
      if(LVERB>=2 || verb/* HERE HERE */)
	printf("T=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:misscnt=%d->%d,falsecnt=%d->%d\n",T,U,I,J,Lijy,Lijx,misscnt,misscnt+(I-max(1,Lijy)),falsecnt,falsecnt + (J-max(1,Lijx)));

      if(DEBUG>=1+RELEASE && !(I-max(1,Lijy) >= 0)){
        #pragma omp critical
        {
	  printf("alignFP:Yid=%d,Xid=%d,or=%d,sc=%d:N=%d,M=%d,score=%0.4f:I=%d,Lijy=%d,misscnt=%d\n",
		   Yid,Xid,orientation,scaleID,N,M,lscore,I,Lijy,misscnt);
	  fflush(stdout);
	  assert(I - max(1,Lijy) >= 0);
	}
      }
      misscnt += I-max(1,Lijy);
      falsecnt += J-max(1,Lijx);
    }
  }
  if(DEBUG && !((Poutlier <= 0.0) ? F== U-1 : F <= U-1)){
    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:Poutlier=%0.4e,F=%d,U=%d\n",
	   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,Poutlier,F,U);
    fflush(stdout);
    assert((Poutlier <= 0.0) ? F== U-1 : F <= U-1);
  }

  /* right end */
  int Rijy = p->Rij1;
  int Rijx = p->Rij2;
  if(LVERB>=2 || verb/* HERE HERE */)
    printf("T=%d/%d:I=%d,J=%d:F=%d,N=%d,M=%d,Rijy=%d,Rijx=%d,Rend=%d:misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	   U,U,I,J,F,N,M,Rijy,Rijx,p->Rend,misscnt,misscnt+min(N,Rijy)-I,falsecnt,falsecnt + min(M,Rijx)-J,OutlierCnt,OutlierMN);
  if(DEBUG>=1+RELEASE && !(min(N,Rijy)-I >= 0)){
    #pragma omp critical
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num ();
#endif

      printf("tid=%d:Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:T=%d/%d:I=%d,J=%d:F=%d,N=%d,M=%d,Rijy=%d,Rijx=%d:misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	     tid,Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,U,U,I,J,F,N,M,Rijy,Rijx,misscnt,misscnt+min(N,Rijy)-I,falsecnt,falsecnt + min(M,Rijx)-J,OutlierCnt,OutlierMN);
      for(int T = 0; T < U; T++)
	printf("   T=%d/%d:I=%d,J=%d\n",T,U,p->sites1[T],p->sites2[T]);
      fflush(stdout);
      assert(min(N,Rijy)-I >= 0);
    }
  }
  misscnt += min(N,Rijy)-I;
  if(DEBUG>=1+RELEASE) assert(min(M,Rijx)-J >= 0);
  falsecnt += min(M,Rijx)-J;

  if(DEBUG) assert(OutlierMN > 0);

  double logOutlierMN = 0.0;
  if(OUTLIER_PV && OutlierMN > 1e+300){
    G = H = -1;
    for(int T=0; T < U; T++, G = I, H = J){
      I = p->sites1[T];
      J = p->sites2[T];
      if(T > 0 && p->iscore[T] > p->outscore[T] + Poutlier + (PFLOAT)0.01){
	int m = J-H;
	int n = I-G;
	logOutlierMN += log((double)(m*n));
      }
    }
  } else
    logOutlierMN = log(OutlierMN);

  int EndOutlierCnt = (p->Lend <= -2 ? 1 : 0) + (p->Rend <= -2 ? 1 : 0);
  int EndOutlierNL = 1;
  if(OUTLIER_PV && EndOutlierPV && EndOutlierCnt){
    if(EndOutlierCnt>=2){
      int L = min(p->sites1[U-1] - p->sites1[0], p->sites2[U-1] - p->sites2[0]);
      int T = min(N, M);
      EndOutlierNL = T - L + 1;
    } else
      EndOutlierNL = 2;
  }

  if(DEBUG>=1+RELEASE && !(misscnt>=0)){
    #pragma omp critical
    {
      printf("alignFP:Yid=%d,Xid=%d,or=%d,sc=%d:N=%d,M=%d,score=%0.4f:I=%d,J=%d,Rijy=%d,misscnt=%d\n",
	       Yid,Xid,orientation,scaleID,N,M,lscore,I,J,Rijy,misscnt);
      fflush(stdout);
      assert(misscnt >= 0);
    }
  }
  double logFP = pvalue(Xfrag,Yfrag, F, misscnt,falsecnt,OutlierCnt,logOutlierMN, EndOutlierNL, Ylambda, res[0]*PixelLen, verb, Yid, Xid,orientation);
  if(DEBUG && !(isfinite(logFP))){
    (void)pvalue(Xfrag,Yfrag,F, misscnt, falsecnt, OutlierCnt, logOutlierMN, EndOutlierNL,Ylambda,res[0]*PixelLen,1,Yid,Xid,orientation);
    fprintf(stderr,"alignFP:I=%d,J=%d:logFP=%e,misscnt=%d(Rijy=%d),falsecnt=%d(Rijx=%d),Outliers=%d,MN=%0.1f\n",
	    I,J,logFP,misscnt,Rijy,falsecnt,Rijx,OutlierCnt,OutlierMN);
    fflush(stderr);
    assert(isfinite(logFP));
  }
  if(BackgroundSNR){
    printf("BackgroundSNR not supported during alignment merger with -pairsplit\n");
    fflush(stdout);exit(1);
  }

  if(U > MAX_ALLOCA)
    delete [] Xfrag;

  return logFP;
}


/** extract alignment with right end at Y[bestI],X[bestJ] from DP array A.field[][] and return pvalue logFP() */
static double alignFP(AL A, int bestI, int bestJ, int bestR,
		      FLOAT *Y, FLOAT *X, int N, int M,
		      int orientation, int scaleID, int Yid, int Xid, PFLOAT **delY, PFLOAT **delX, double lscore, int verb, int Ascorevalid, int myoutlierExtend, int myoutlierBC, int tid,
		      int *JMIN, int *JMAX, int IMIN, int IMAX)
{
  if(VERB/* HERE HERE >=2 */ && verb){
    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,bestR=%d,lscore=%0.6f\n",Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,bestR,lscore);
    fflush(stdout);
  }

  /* allocate space to store alignment in reverse order */
  int *Ilist, *Jlist, *outlier;
  if(M > MAX_ALLOCA){
    Ilist = new int[M*3];
    Jlist = &Ilist[M];
    outlier = &Ilist[M*2];
  } else {
    Ilist = (int *)alloca(M*sizeof(int));
    Jlist = (int *)alloca(M*sizeof(int));
    outlier = (int *)alloca(M*sizeof(int));
  }
  int K = 0;

  int I = bestI;
  int J = bestJ;

  if(myoutlierExtend){
    while(I >= 0){
      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
      if(DEBUG>=2){
	assert(I >= 1 && I <= N);
	assert(J >= 1 && J <= M);
      }

      Ilist[K] = I;
      Jlist[K] = J;
      int G = FLAT_A(A.G, I, J);
      if(G > 0){
	int H = FLAT_A(A.H, I, J);
	if(DEBUG>=2) assert(H > 0);
	if(DEBUG>=2) assert(K+1 < M);
	if(DEBUG>=2) assert(IMIN <= G && G <= IMAX && JMIN[G] <= H && H <= JMAX[G]);

	PFLOAT y = Y[I] - Y[G];
	PFLOAT x = X[J] - X[H];
	if(Poutlier > 0.0){
	  PFLOAT Bias,Pen;
	  SintDetail(x,y, J-H,I-G,J,I,Bias,Pen);
	  PFLOAT OutPen = OutlierPenalty;
	  if(myoutlierBC)
	    OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-G];
	  PFLOAT iscore = OutlierBias + (OUTLIER_DELTA(x-y) ? max(Pen+Bias,OutPen + Bias*biasWToutlierF) : Pen + Bias);
	  if(DEBUG>=2 && Ascorevalid){
	    PFLOAT iscore2 = FLAT_A(A.score, I, J) - FLAT_A(A.score, G, H);
	    if(!(fabs(iscore2 - iscore) < 0.001)){
	      #pragma omp critical
	      {
		printf("alignFP:Yid=%d,Xid=%d,or=%d,sc=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,Ascorevalid=%d\n",
		       Yid,Xid,orientation,scaleID,bestI,bestJ,bestR,N,M,Ascorevalid);
		printf("K=%d:I=%d,J=%d,G=%d,H=%d,y=%0.4f,x=%0.4f,OutlierBias=%0.6f,Bias=%0.6f,Pen=%0.6f,OutPen=%0.6f,iscore=%0.6f,iscore2=%0.6f (A.score[I][J]=%0.6f,A.score[G][H]=%0.6f), OUTLIER_DELTA(x-y)=%d\n",
		       K,I,J,G,H,y,x,OutlierBias,Bias,Pen,OutPen,iscore,iscore2,FLAT_A(A.score,I,J),FLAT_A(A.score,G,H), OUTLIER_DELTA(x-y));
		printf("\t OutlierPenalty= %0.6f, biasWToutlierF= %0.6f, myoutlierExtend=%d,myoutlierBC=%d\n",OutlierPenalty,biasWToutlierF,myoutlierExtend,myoutlierBC);
		fflush(stdout);
		assert(fabs(iscore2 - iscore) < 0.001);
	      }
	    }
	  }
	  PFLOAT outscore = OutlierBias + Bias + Pen;
	  outlier[K+1] = (iscore > outscore + Poutlier + (PFLOAT)0.01 || (myoutlierExtend && (J <= H || I <= G))) ? 1 : 0;
	  if(VERB/* HERE HERE >=2 */ && verb)
	    printf("K=%d:I=%d,J=%d,G=%d,H=%d:x=%0.3f,y=%0.3f,iscore=%0.6f,outscore=%0.6f:outlier[K+1]=%d\n",K,I,J,G,H,x,y,iscore,outscore,outlier[K+1]);
	} else
	  outlier[K+1] = 0;
	J = H;
      }
      K++;
      I = G;
    }
  } else {/* General Case : myoutlierExtend==0 */
    while(I >= 0){
      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
      if(DEBUG>=2){
	assert(I >= 1 && I <= N);
	assert(J >= 1 && J <= M);
      }

      Ilist[K] = I;
      Jlist[K] = J;
      int G = FLAT_A(A.G, I, J);
      if(G > 0){
	int H = FLAT_A(A.H, I, J);
	if(DEBUG>=2) assert(IMIN <= G && G <= IMAX && JMIN[G] <= H && H <= JMAX[G]);
	if(DEBUG>=2) assert(1 <= H && K+1 < M);
	PFLOAT y = delY[I][G-I+DELTA];//Y[I] - Y[G];
	if(DEBUG>=2) assert(fabs(Y[I] - Y[G] - y) < 0.001);
	if(DEBUG>=2 && !(1 <= H && H < J && J-H <= DELTA && J <= M)){
	  #pragma omp critical
	  {
	    printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d,M=%d\n",Yid,Xid,orientation,scaleID,I,J,G,H,M);
	    fflush(stdout);
	    assert(1 <= H && H < J && J-H <= DELTA && J <= M);
	  }
	}
	if(DEBUG>=2) assert(&FLAT_DELX(delX, J-H, J, M) == &delX[J-H][J]);
	PFLOAT x = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];//X[J] - X[H];
	if(DEBUG>=2) assert(fabs(X[J] - X[H] - x) < 0.001);
	if(Poutlier > 0.0){
	  PFLOAT Bias,Pen;
	  SintDetail(x,y, J-H,I-G,J,I,Bias,Pen);
	  PFLOAT OutPen = OutlierPenalty;
	  if(myoutlierBC)
	    OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-G];
	  PFLOAT iscore = OutlierBias + (OUTLIER_DELTA(x-y) ? max(Pen + Bias, OutPen + Bias*biasWToutlierF) : Pen + Bias);
	  if(DEBUG>=2 && Ascorevalid){
	    PFLOAT iscore2 = FLAT_A(A.score, I, J) - FLAT_A(A.score, G, H);
	    if(!(fabs(iscore2 - iscore) < 0.001)){
	      #pragma omp critical
	      {
		printf("alignFP:Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,Ascorevalid=%d\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,bestR,N,M,Ascorevalid);
		printf("K=%d:I=%d,J=%d,G=%d,H=%d,y=%0.4f,x=%0.4f,OutlierBias=%0.6f,Bias=%0.6f,Pen=%0.6f,OutPen=%0.6f,iscore=%0.6f,iscore2=%0.6f (A.score[I][J]=%0.6f,A.score[G][H]=%0.6f), OUTLIER_DELTA(x-y)=%d\n",
		       K,I,J,G,H,y,x,OutlierBias,Bias,Pen,OutPen,iscore,iscore2,FLAT_A(A.score,I,J),FLAT_A(A.score,G,H), OUTLIER_DELTA(x-y));
		printf("\t OutlierPenalty= %0.6f, biasWToutlierF= %0.6f, myoutlierExtend=%d,&FLAT_A(A.score,I,J)=%p,A.score=%p,tid=%d\n",
		       OutlierPenalty,biasWToutlierF,myoutlierExtend,&FLAT_A(A.score,I,J),A.score,tid);
		fflush(stdout);
		assert(fabs(iscore2 - iscore) < 0.001);
	      }
	    }
	  }
	  PFLOAT outscore = OutlierBias + Bias + Pen;
	  outlier[K+1] = (iscore > outscore + Poutlier + (PFLOAT)0.01) ? 1 : 0;
	  if(VERB/* HERE HERE >=2 */ && verb)
	    printf("K=%d:I=%d,J=%d,G=%d,H=%d:x=%0.3f,y=%0.3f,iscore=%0.6f,outscore=%0.6f:outlier[K+1]=%d\n",K,I,J,G,H,x,y,iscore,outscore,outlier[K+1]);
	} else
	  outlier[K+1] = 0;
	J = H;
      }
      K++;
      I = G;
    }
  }/* General Case : myoutlierExtend==0 */

  int EndOutlierCnt = (bestR <= -2 ? 1 : 0) + (I <= -2 ? 1 : 0);

  if(DEBUG) assert(K>0);
  if(DEBUG) assert(K <= M);
    
  /* allocate alignment fragment sizes */
  FLOAT *Xfrag, *Yfrag;
  if(K > MAX_ALLOCA){
    Xfrag = new FLOAT[2*(K-1)];
    Yfrag = &Xfrag[K-1];
  } else {
    Xfrag = (FLOAT *)alloca((K-1)*sizeof(FLOAT));
    Yfrag = (FLOAT *)alloca((K-1)*sizeof(FLOAT));
  }

  int misscnt = 0;
  int falsecnt = 0;
  int OutlierCnt = 0;/* number of internal outliers */
  double OutlierMN = 1.0;/* product over all outlier of m x n (the number of intervals in the internal outlier aligned segment) */

  int G= -1,H= -1;
  int F = 0;/* segment (aligned interval) count */
  for(int T= 0; T < K; T++, G = I, H = J){
    I = Ilist[K-1-T];
    J = Jlist[K-1-T];
    if(T>0){
      if(DEBUG) assert(G > 0 && H > 0);
      if(outlier[K-T]){
	if(LVERB>=2)
	  printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:outlier[K-T]=%d\n",T,K,I,J,G,H,outlier[K-T]);
	if(OUTLIER_PV){
	  OutlierCnt++;
	  int m = max(1,J - H);
	  int n = max(1,I - G);
	  if(DEBUG>=2 && !myoutlierExtend && !(J > H && I > G)){
	    printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:outlier[K-T]=%d\n",T,K,I,J,G,H,outlier[K-T]);
	    fflush(stdout);
	    assert(J > H);
	    assert(I > G);
	  }
	  if(DEBUG>=2 && !(OutlierMN*m*n > 0)){
	    printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:outlier[K-T]=%d:m=%d,n=%d,OutlierMN=%0.1f -> %0.1f,OutlierCnt=%d\n",
		   T,K,I,J,G,H,outlier[K-T],m,n,OutlierMN,OutlierMN*m*n,OutlierCnt);
	    fflush(stdout);
	    assert(OutlierMN*m*n > 0);
	  }
	  OutlierMN *= (double)(m*n);
	  if(DEBUG>=2 && !(OutlierMN>0)){
	    printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:outlier[K-T]=%d:m=%d,n=%d,OutlierMN=%0.1f\n",
		   T,K,I,J,G,H,outlier[K-T],m,n,OutlierMN);
	    fflush(stdout);
	    assert(OutlierMN > 0);
	  }
	}
	continue;
      }
      if(LVERB>=2)
	printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:F=%d:y=%0.3f,x=%0.3f,misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	       T,K,I,J,G,H,F,Y[I]-Y[G],X[J]-X[H],misscnt,misscnt+(I-G-1),falsecnt,falsecnt + (J-H-1),OutlierCnt,OutlierMN);
      Yfrag[F] = Y[I] - Y[G];
      Xfrag[F] = X[J] - X[H];
      if(DEBUG>=1+RELEASE && !(Yfrag[F] > 0.0 && Xfrag[F] > 0.0 && I > G && J > H)){
	PFLOAT Bias,Pen;
	SintDetail(Xfrag[F],Yfrag[F], J-H, I-G, J, I, Bias, Pen);
	PFLOAT OutPen = OutlierPenalty;
	if(myoutlierBC)
	  OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-G];
	PFLOAT iscore = OutlierBias + (OUTLIER_DELTA(Xfrag[F]-Yfrag[F]) ? max(Pen + Bias, OutPen + Bias*biasWToutlierF) : Pen + Bias);
	PFLOAT outscore = OutlierBias + Bias + Pen;

	printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:F=%d:y=%0.3f,x=%0.3f,misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	       T,K,I,J,G,H,F,Y[I]-Y[G],X[J]-X[H],misscnt,misscnt+(I-G-1),falsecnt,falsecnt + (J-H-1),OutlierCnt,OutlierMN);
	printf("  outlier[K-T]=%d:iscore=%0.6f,outscore=%0.6f,Bias=%0.6f,Pen=%0.6f,OutPen=%0.6f,OutlierBias=%0.6f\n",outlier[K-T],iscore,outscore,Bias,Pen,OutPen,OutlierBias);
	fflush(stdout);
	assert(Yfrag[F] > 0.0);
	assert(Xfrag[F] > 0.0);
	assert(I > G);
	assert(J > H);
      }
      F++;
      misscnt += I-G-1;
      falsecnt += J-H-1;
    } else {/* left end */
      if(DEBUG>=2 && !(Ch > 0.0 || PoutlierEnd > 0.0 || (RepeatRec && RepeatMaxShift > 0 && RepeatPvalueRatio > 0.0)) && !(A.G[I][J] > -2)){
	#pragma omp critical
	{
	  printf("alignFP:Yid=%d,Xid=%d,or=%d,sc=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,score=%0.4f:U=%d,I=%d,J=%d,Ch=%0.6e,PoutlierEnd=%0.6e,A.G[I][J]=%d,A.score[I][J]=%0.6f\n",
		 Yid,Xid,orientation,scaleID,bestI,bestJ,bestR,N,M,lscore,K,I,J,Ch,PoutlierEnd,A.G[I][J],A.score[I][J]);
	  fflush(stdout);
	  assert(A.G[I][J] > -2);
	}
      }
      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
      int Lijy = (FLAT_A(A.G, I, J) <= -2) ? I : FLAT_A(A.Lijy, I, J);
      int Lijx = (FLAT_A(A.G, I, J) <= -2) ? J : FLAT_A(A.Lijx, I, J);
      if(DEBUG>=2 && !(Lijy >= 0 && Lijy <= I)){
        #pragma omp critical
        {
	  printf("alignFP:Yid=%d,Xid=%d,or=%d,sc=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,score=%0.4f:U=%d,I=%d,J=%d,Lijy=%d,Lijx=%d,misscnt=%d,falsecnt=%d\n",
		   Yid,Xid,orientation,scaleID,bestI,bestJ,bestR,N,M,lscore,K,I,J,Lijy,Lijx,misscnt,falsecnt);
	  fflush(stdout);
	  assert(Lijy >= 0 && Lijy <= I);
	}
      }
      if(DEBUG>=2) assert(Lijx >= 0 && Lijx <= J);
      if(LVERB>=2)
	printf("T=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:misscnt=%d->%d,falsecnt=%d->%d\n",T,K,I,J,Lijy,Lijx,misscnt,misscnt+(I-max(1,Lijy)),falsecnt,falsecnt + (J-max(1,Lijx)));

      misscnt += I-max(1,Lijy);
      if(DEBUG>=1+RELEASE && !(misscnt>=0 && (I >= max(1,Lijy)))){
        #pragma omp critical
        {
	  printf("alignFP:Yid=%d,Xid=%d,or=%d,sc=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,score=%0.4f:I=%d,Lijy=%d,misscnt=%d\n",
		   Yid,Xid,orientation,scaleID,bestI,bestJ,bestR,N,M,lscore,I,Lijy,misscnt);
	  fflush(stdout);
	  assert(misscnt >= 0);
	  assert(I >= max(1,Lijy));
	}
      }
      falsecnt += J-max(1,Lijx);
      if(DEBUG>=1+RELEASE) 	assert(Lijx >= 0 && max(1,Lijx) <= J);
    }
  }
  if(DEBUG) assert(I==Ilist[0] && I==bestI);
  if(DEBUG) assert(J==Jlist[0] && J==bestJ);
  if(DEBUG) assert((Poutlier <= 0.0) ? F== K-1 : F <= K-1);

  /* right end */
  if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
  if(DEBUG>=2 && !(Ch > 0.0 || PoutlierEnd > 0.0)) assert(bestR > -2);
  int Rijy = (bestR <= -2) ? I : FLAT_A(A.Rijy, I, J);
  int Rijx = (bestR <= -2) ? J : FLAT_A(A.Rijx, I, J);
  if(DEBUG>=2) assert(I <= Rijy && Rijy <= N+1);
  if(DEBUG>=2) assert(J <= Rijx && Rijx <= M+1);
  if(LVERB>=2)
    printf("T=%d/%d:I=%d,J=%d:F=%d,N=%d,M=%d,Rijy=%d,Rijx=%d:misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	   K,K,I,J,F,N,M,Rijy,Rijx,misscnt,misscnt+min(N,Rijy)-I,falsecnt,falsecnt + min(M,Rijx)-J,OutlierCnt,OutlierMN);
  if(DEBUG>=1+RELEASE && !(min(N,Rijy)-I >= 0)){
    #pragma omp critical
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num ();
#endif
      printf("tid=%d:Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:T=%d/%d:I=%d,J=%d:F=%d,N=%d,M=%d,Rijy=%d,Rijx=%d:misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	     tid,Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,K,K,I,J,F,N,M,Rijy,Rijx,misscnt,misscnt+min(N,Rijy)-I,falsecnt,falsecnt + min(M,Rijx)-J,OutlierCnt,OutlierMN);
      printf("bestI=%d,bestJ=%d,bestR=%d,FLAT_A(A.Rijy,I,J)=%d (%p)\n",bestI,bestJ,bestR, FLAT_A(A.Rijy,I,J),  &FLAT_A(A.Rijy,I,J));
      for(int T = 0; T < K; T++)
	printf("   T=%d/%d:I=%d,J=%d\n",T,K,Ilist[K-1-T],Jlist[K-1-T]);
      fflush(stdout);
      assert(min(N,Rijy)-I >= 0);
    }
  }
  misscnt += min(N,Rijy)-I;
  if(DEBUG>=1+RELEASE) assert(min(M,Rijx)-J >= 0);
  falsecnt += min(M,Rijx)-J;

  if(DEBUG) assert(OutlierMN > 0);

  double logOutlierMN = 0.0;
  if(OUTLIER_PV && OutlierMN > 1e+300){
    G = H = -1;
    for(int T=0;T < K;T++, G = I, H = J){
      I = Ilist[K-1-T];
      J = Jlist[K-1-T];
      if(T>0){
	if(outlier[K-T]){
	  int m = J-H;
	  int n = I-G;
	  logOutlierMN += log((double)(m*n));
	}
      }
    }
  } else
    logOutlierMN = log(OutlierMN);

  int EndOutlierNL = 1;
  if(OUTLIER_PV && EndOutlierPV && EndOutlierCnt){
    if(EndOutlierCnt >= 2){
      int L = min(Ilist[0] - Ilist[K-1], Jlist[0] - Jlist[K-1]);
      int T = min(N, M);
      EndOutlierNL = T - L + 1;
    } else
      EndOutlierNL = 2;
  }

  if(DEBUG>=1+RELEASE && !(misscnt>=0)){
    #pragma omp critical
    {
      printf("alignFP:Yid=%d,Xid=%d,or=%d,sc=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,score=%0.4f:I=%d,J=%d,Rijy=%d,misscnt=%d\n",
	     Yid,Xid,orientation,scaleID,bestI,bestJ,bestR,N,M,lscore,I,J,Rijy,misscnt);
      fflush(stdout);
      assert(misscnt >= 0);
    }
  }
  double logFP = pvalue(Xfrag,Yfrag, F, misscnt,falsecnt,OutlierCnt,logOutlierMN,EndOutlierNL,Ylambda, res[0]*PixelLen,verb, Yid, Xid,orientation);
  if(DEBUG && !(isfinite(logFP))){
    (void)pvalue(Xfrag,Yfrag,F, misscnt,falsecnt,OutlierCnt,logOutlierMN,EndOutlierNL,Ylambda,res[0]*PixelLen,1,Yid,Xid,orientation);
    fprintf(stderr,"alignFP:I=%d,J=%d:logFP=%e,misscnt=%d(Rijy=%d),falsecnt=%d(Rijx=%d),Outliers=%d,MN=%0.1f\n",
	    I,J,logFP,misscnt,Rijy,falsecnt,Rijx,OutlierCnt,OutlierMN);
    fflush(stderr);
    assert(isfinite(logFP));
  }
  if(BackgroundSNR){/* compute pvalue for site matches based on Kalmagorov-Smirnov distributions */
    double HalfByLog10 = 0.5/log(10.0);
    if(!YYmap[Yid]->SNRcnt[0]){
      printf("Missing SNR values for id=%lld (needed for SNR based pvalues)\n",YYmap[Yid]->id);
      fflush(stdout);exit(1);
    }
    if(!XXmap[Xid]->SNRcnt[0]){
      printf("Missing SNR values for id=%lld (needed for SNR based pvalues)\n",XXmap[Xid]->id);
      fflush(stdout);exit(1);
    }

    for(int U = 0; U < K; U++){
      int I = Ilist[U];
      int J = Jlist[U];
      int JR = orientation ? M+1-J : J;

      if(!YYmap[Yid]->SNRcnt[0][I]){
	if(VERB>=2){
	  printf("Missing SNR values for id=%lld, site=%d/%d (skipping SNR based Pvalue adjustment)\n",YYmap[Yid]->id,I,N);
	  fflush(stdout);
	}
	continue;
      }

      if(!XXmap[Xid]->SNRcnt[0][JR]){
	if(VERB>=2){
	  printf("Missing SNR values for id=%lld, site=%d/%d (skipping SNR based Pvalue adjustment)\n",XXmap[Xid]->id,JR,M);
	  fflush(stdout);
	}
	continue;
      }

      if(DEBUG>=2){
	for(int i = 1; i < YYmap[Yid]->SNRcnt[0][I]; i++)
	  if(!(YYmap[Yid]->SNRdist[0][I][i] >= YYmap[Yid]->SNRdist[0][I][i-1])){
	    printf("Yid=%d(id=%lld),I=%d:i=%d/%d,SNRdist[i-1]=%0.8e,SNRdist[i]=%0.8e\n",
		   Yid,YYmap[Yid]->id,I,i,YYmap[Yid]->SNRcnt[0][I],YYmap[Yid]->SNRdist[0][I][i-1],YYmap[Yid]->SNRdist[0][I][i]);
	    fflush(stdout);
	    assert(YYmap[Yid]->SNRdist[0][I][i] >= YYmap[Yid]->SNRdist[0][I][i-1]);
	  }
	for(int i = 1; i < XXmap[Xid]->SNRcnt[0][JR]; i++)
	  if(!(XXmap[Xid]->SNRdist[0][JR][i] >= XXmap[Xid]->SNRdist[0][JR][i-1])){
	    printf("Xid=%d,J=%d:i=%d/%d,SNRdist[i-1]=%0.8e,SNRdist[i]=%0.8e\n",
		   Xid,JR,i,XXmap[Xid]->SNRcnt[0][J],XXmap[Xid]->SNRdist[0][JR][i-1],XXmap[Xid]->SNRdist[0][JR][i]);
	    fflush(stdout);
	    assert(XXmap[Xid]->SNRdist[0][JR][i] >= XXmap[Xid]->SNRdist[0][JR][i-1]);
	  }
      }

      double Pv12 = KSPvalue(YYmap[Yid]->SNRcnt[0][I],YYmap[Yid]->SNRdist[0][I],XXmap[Xid]->SNRcnt[0][JR],XXmap[Xid]->SNRdist[0][JR],verb);
      double PvBg1 = KSPvalue(YYmap[Yid]->SNRcnt[0][I],YYmap[Yid]->SNRdist[0][I],bgSNRcnt[0],bgSNR[0],verb);
      double PvBg2 = KSPvalue(XXmap[Xid]->SNRcnt[0][JR],XXmap[Xid]->SNRdist[0][JR],bgSNRcnt[0],bgSNR[0],verb);

      if(LVERB || verb)
	printf("Yid=%d,Xid=%d,or=%d,sc=%d:SNR match for I=%d,J=%d(JR=%d):n1=%d(gmSNR=%0.2f),n2=%d(gmSNR=%0.2f),nBG=%d: LogPv12=%0.6f,LogPvBg1=%0.6f,LogPvBg2=%0.6f: -LogPV += %0.4f\n",
		 Yid,Xid,orientation,scaleID,I,J,JR,YYmap[Yid]->SNRcnt[0][I],YYmap[Yid]->SNRgmean[0][I],XXmap[Xid]->SNRcnt[0][JR], XXmap[Xid]->SNRgmean[0][JR], bgSNRcnt[0], Pv12/log(10.0),PvBg1/log(10.0),PvBg2/log(10.0), -HalfByLog10*(PvBg1+PvBg2 - 2.0*Pv12));

      logFP -= HalfByLog10*(PvBg1 + PvBg2 - 2.0*Pv12);
    }
    if(LVERB || verb)
      fflush(stdout);
  }

  if(M > MAX_ALLOCA)
    delete [] Ilist;
  if(K > MAX_ALLOCA)
    delete [] Xfrag;

  return logFP;
}

static int VecMessage = 0;
static int NoVecMessage = 0;

//#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

class CThreadMem {
 public:
  PFLOAT *Fmem;
  int *Imem;// typically Imem == &((char *)Fmem)[Fsiz * sizeof(PFLOAT)]
  long long MemSiz;/* size of mmap allocation memory in current thread : base pointer is Fmem */
  long long MemUsed;/* estimated actively used part of Fmem[0..MemSiz-1] is Fmem[0..MemUsed-1] */
  long long MemReserved;/* amount by which global VmRSSreserved was incremented */
  long long VmSize;
  long long VmRSS;
  long long VmSwap;
  double lasttime;/* The value of wtime() the last time madvise() was called on Fmem[0..MemSiz-1] : updated at least every VMRSS_CHECK sec (if memory use is above 60% of -maxmem) */
  double lastGetStatm;/* The value of wtime() the last time VmSize & VmRSS were updated from /proc/<pid>/statm : updated at least every VMRSS_CHECK seconds*/
  double lastGetMem;/* The value of wtime() the last time VmSwap was updated from /proc/<pid>/status by calling getmem() : updated at least every VMSWAP_CHECK seconds */
  int fd;/* file descriptor for /proc/<pid>/statm */

  size_t madvise_calls;
  size_t getmem_calls;
  size_t getstatm_calls;
  size_t wait_calls;
  size_t release_calls;
  size_t critical_calls;
  size_t atomic_calls;

  size_t pad1,pad2,pad3,pad4,pad5;/* round to 3 cache lines = 3 * 64 bytes */

  CThreadMem(){
    Fmem = NULL;
    Imem = NULL;
    VmSwap = VmSize = VmRSS = MemSiz = MemUsed = MemReserved = 0;
    lasttime = lastGetStatm = lastGetMem = 0.0;

    madvise_calls = getmem_calls = getstatm_calls = wait_calls = release_calls = critical_calls = atomic_calls = 0;

    char filename[PATH_MAX];
    pid_t pid = getpid();
    sprintf(filename,"/proc/%lld/statm",(long long)pid);
    fd = open(filename,O_RDONLY);
    if(fd < 0){
      char *err = strerror(errno);
      printf("CThreadMem(): Failed to open for read %s: %s\n",filename,err);
      fflush(stdout);
      char buf[BUFSIZ];
      sprintf(buf,"cat /proc/%lld/statm",(long long)pid);
      printf("$ %s\n",buf);
      fflush(stdout);
      int ret = system(buf);
      if(ret){
	printf("system(%s) failed: return code = %d\n", buf, ret);
	fflush(stdout);
      }
      fflush(stdout);  exit(1);
    }
    
    char buf[128];
    ssize_t r = read(fd,buf,(size_t)127);
    if(r <= 0){
      char *err = strerror(errno);
      printf("CThreadMem(): read of %s failed (r=%ld,fd=%d): %s\n",filename,r,fd,err);
      fflush(stdout);exit(1);
    }
    buf[r] = '\0';
    
    size_t rVmSize, rVmRSS;
    int cnt = sscanf(buf,"%lu\t%lu", &rVmSize, & rVmRSS);
    if(cnt != 2){
      printf("CThreadMem(): read+sscanf of /proc/self/statm failed (r=%ld,cnt=%d): buf=%s\n",r,cnt,buf);
      fflush(stdout);exit(1);
    }

    VmSize = rVmSize * PAGE;
    VmRSS = rVmRSS * PAGE;

    lastGetStatm = wtime();
    if(THREAD_DEBUG)
      getstatm_calls++;
  };

  void GetStatm(double wt){
    off_t s = lseek(fd, 0L, SEEK_SET);
    if(s < 0){
      char *err = strerror(errno);
      printf("CThreadMem::GetStatm(): lseek to start of /proc/self/statm failed (s=%lld): %s\n",(long long int)s,err);
      fflush(stdout);exit(1);
    }

    char buf[128];
    ssize_t r = read(fd,buf,(size_t)127);
    if(r <= 0){
      char *err = strerror(errno);
      printf("CThreadMem(): read of /proc/self/statm failed (r=%ld,fd=%d): %s\n",r,fd,err);
      fflush(stdout);exit(1);
    }
    buf[r] = '\0';
    
    size_t rVmSize, rVmRSS;
    int cnt = sscanf(buf,"%lu\t%lu", &rVmSize, & rVmRSS);
    if(cnt != 2){
      printf("CThreadMem(): read+sscanf of /proc/self/statm failed (r=%ld,cnt=%d): buf=%s\n",r,cnt,buf);
      fflush(stdout);exit(1);
    }

    VmSize = rVmSize * PAGE;
    VmRSS = rVmRSS * PAGE;
    
    lastGetStatm = wt;
    if(THREAD_DEBUG)
      getstatm_calls++;
  };

  void GetMem(double wt){
    getmem(VmSize,VmRSS,VmSwap);

    lastGetStatm = lastGetMem = wt;
    if(THREAD_DEBUG)
      getmem_calls++;
  };

  ~CThreadMem(){
    if(fd >= 0){
      if(close(fd) < 0){
	char *err = strerror(errno);
	printf("close() for /proc/self/statm failed: %s\n",err);
	fflush(stdout);
      }
    }
  }
};
static CThreadMem *ThreadMem = NULL;

static long long MaxMemSiz = 0;// maximum real memory
static long long vmemTotsiz = 0;// current virtual memory usage across all threads (only used if VMEMCHECK)
static int PausedThreads = 0, NumThreads = 0;// current number of threads paused waiting for memory : accessed using omp critical(PausedThreads)
static long long MaxMemSiz2;// MaxMemSiz * 0.9 : Max Value of VmRSS + VmSwap + nMemSiz for large jobs

static long long VmRSSreserved = 0;/* shared variable to track reserved memory by all threads */

/* this function is called after completing a pairwise alignment (if -RAmem was specified) to release memory */
static void ReleaseMemory(CThreadMem *Mem, long long nMemSiz, int tid, int Yid, int Xid, int N, int M, PFLOAT* &Fmem, int force = 0)
{
  if(THREAD_DEBUG) Mem->release_calls++;

  if(RA_MIN_TIME && RA_MIN_MEM < 8 && Mem->MemUsed > 0){/* call madvise(Fmem,Mem->MemSiz,MAV_DONTNEED) once every RA_MIN_TIME sec (in each thread) */
    if(DEBUG/* HERE >=2 */) assert(nMemSiz <= Mem->MemUsed);

    double wt = wtime();

    if(force || wt >= Mem->lasttime + VMRSS_CHECK){

      long long sVmSize,sVmRSS,sVmSwap,sVmHWM;// local copies, so global shared values are not modified
      
      if(VERB>=3)
	getmem(sVmSize,sVmRSS,sVmSwap,&sVmHWM);

      if(wt > Mem->lastGetMem + VMSWAP_CHECK){/* update both VmRSS & VmSwap */
	if(THREAD_DEBUG) Mem->getmem_calls++;
	Mem->GetMem(wt);
      } else if(wt > Mem->lastGetStatm + VMRSS_CHECK && (wt > Mem->lastGetStatm + VMRSS_CHECK2 || Mem->VmRSS + Mem->VmSwap >= MaxMemSiz2 * RA_MAX_MEM)){/* update only VmRSS */
	if(THREAD_DEBUG) Mem->getstatm_calls++;
	Mem->GetStatm(wt);
      }

      if(Mem->VmRSS + Mem->VmSwap >= MaxMemSiz2 * RA_MAX_MEM){// NEW71 : don't wast time releasing memory if we are using less than RA_MAX_MEM * maxmem * 0.9 

	if(THREAD_DEBUG) Mem->madvise_calls++;

	if(madvise(Fmem, Mem->MemSiz /* Mem->MemUsed */, MADV_DONTNEED)){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("madvise(%p,%lld,MADV_DONTNEED) failed(MemUsed=%lld):errno=%d:%s\n",Fmem,Mem->MemSiz,Mem->MemUsed,eno,err);
	  fflush(stdout);exit(1);
	}
	Mem->lasttime = wt;

	if(VERB>=3 /* || force*/){
	  double nwt = wtime();
	  long long rVmSize,rVmRSS,rVmSwap,rVmHWM;
	
	  getmem(rVmSize,rVmRSS,rVmSwap,&rVmHWM);

#pragma omp critical(PausedThreads)
	  {
	    if(PausedThreads > 0 || force){
	      if(force)
		printf("tid=%d: final ReleaseMemory: madvise on %0.4f/%0.4f Gb : VmSize= %0.4f, VmRSS= %0.4f -> %0.4f(HWM=%0.2f), VmSwap= %0.4f Gb: Paused=%d:wall=%0.6f(elapsed=%0.6f)\n",
		       tid, Mem->MemUsed*1e-9,Mem->MemSiz*1e-9,rVmSize*1e-9,sVmRSS*1e-9,rVmRSS*1e-9, rVmHWM*1e-9, rVmSwap*1e-9, PausedThreads,nwt - wt, nwt);
	      else
		printf("tid=%d: ReleaseMemory: madvise on %0.4f/%0.4f Gb : VmSize= %0.4f, VmRSS= %0.4f -> %0.4f(HWM=%0.2f), VmSwap= %0.4f Gb: (After Yid=%lld,Xid=%lld,N=%d,M=%d),Paused=%d:wall=%0.6f(elapsed=%0.6f)\n",
		       tid, Mem->MemUsed*1e-9,Mem->MemSiz*1e-9,rVmSize*1e-9,sVmRSS*1e-9,rVmRSS*1e-9, rVmHWM*1e-9, rVmSwap*1e-9, YYmap[Yid]->id,XXmap[Xid]->id,N,M,PausedThreads,nwt - wt, nwt);
	      fflush(stdout);
	    }
	  }
	}

	Mem->MemUsed = 0;
      }
    }
  }
  
  if(VMEMCHECK){
    if(DEBUG) Mem->critical_calls++;

    #pragma omp critical(vmemTotsiz)
    {
      vmemTotsiz -= nMemSiz;

      if(VERB>=3 || (DEBUG>=2 && vmemTotsiz < 0)){
	long long rVmSize,rVmRSS,rVmSwap,rVmHWM;// local copies, so global shared values are not modified
	getmem(rVmSize,rVmRSS,rVmSwap,&rVmHWM);

	printf("tid=%d: ReleaseMemory: vmemTotsiz=%lld:After releasing %0.6f Gb: VmRSS= %0.4f(HWM=%0.2f),VmSwap= %0.4fGb: wt=%0.6f\n", 
	  tid, vmemTotsiz, nMemSiz*1e-9, rVmRSS*1e-9, rVmHWM*1e-9,rVmSwap*1e-9,wtime());
	fflush(stdout);
	if(DEBUG) assert(vmemTotsiz >= 0);
      }
    }
  }
}

#include <time.h>

/* Pause current thread if not enough memory is available. */
static void WaitForMemory(CThreadMem *Mem, long long nMemSiz, int tid, int Yid, int Xid, int N, int M, PFLOAT* &Fmem, double RSSratio = 1.0)
{
  if(THREAD_DEBUG) Mem->wait_calls++;

  struct timespec tm;
  tm.tv_sec = floor(VMRSS_SLEEP);
  tm.tv_nsec = (VMRSS_SLEEP - tm.tv_sec) * 1e9;

  int waitcnt = 0;
  double waittime = 0.0;

  long long maxVmRSS2 = MaxMemSiz2;// maximum VmRSS+VmSwap allowed for large memory alignments (must fit within 90% of memory)
  long long maxVmRSS1 = MaxMemSiz;// maxmimum VmRSS+VmSwap allowed for jobs expected to use less than 10% of the per thread amount (so total will fit within 100% of memory)

  bool SmallJob = nMemSiz < (MaxMemSiz / (NumThreads * 10.0));
  long long maxVmRSS = SmallJob ? maxVmRSS1 : maxVmRSS2;// maximum VmRSS+VmSwap allowed for current pairwise alignment

  bool madvise_called = (RA_MIN_MEM >= 8) ? true : false;/* Don't call madvise if RA_MIN_MEM >= 8 */
  bool memfound = false;
  int ThreadNotPaused = 1;

  int TICKS = floor(1.0 / VMRSS_SLEEP + 0.5);
  int TICKS60 = 60*TICKS;
  int RA_TICKS = 10 * TICKS;

  while(!memfound) {
    long long RSSneeded = max(0LL, nMemSiz - Mem->MemUsed) * RSSratio;

    double wt = wtime();

    if(wt > Mem->lastGetMem + VMSWAP_CHECK){/* update both VmRSS & VmSwap */
      if(THREAD_DEBUG) Mem->getmem_calls++;
      Mem->GetMem(wt);
      if(VERB>=2 && !ThreadNotPaused){
        #pragma omp critical(PausedThreads)
	{
	  printf("tid=%d: MaxMem= %0.1f (Lim=%0.2f): GetMem: VmRSS= %0.4f, VmSwap= %0.2f Gb (Yid=%lld,Xid=%lld,N=%d,M=%d),waittime= %0.3f,Paused Threads=%d:wt=%0.6f,%0.6f,%0.6f\n",
	     tid, MaxMemSiz*1e-9, maxVmRSS*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, YYmap[Yid]->id,XXmap[Xid]->id,N,M,waittime,PausedThreads,Mem->lastGetMem,Mem->lastGetStatm,wt);
	  fflush(stdout);
	}
      }
    } else if(wt > Mem->lastGetStatm + VMRSS_CHECK && (wt > Mem->lastGetStatm + VMRSS_CHECK2 || Mem->VmRSS + Mem->VmSwap >= MaxMemSiz2 * RA_MAX_MEM)){/* update only VmRSS */
      if(THREAD_DEBUG) Mem->getstatm_calls++;
      Mem->GetStatm(wt);
      if(VERB>=2 && !ThreadNotPaused){
        #pragma omp critical(PausedThreads)
	{
	  printf("tid=%d: MaxMem= %0.1f (Lim=%0.2f): GetStatm: VmRSS= %0.4f, VmSwap= %0.2f Gb (Yid=%lld,Xid=%lld,N=%d,M=%d),waittime= %0.3f, Paused Threads=%d:wt=%0.6f,%0.6f,%0.6f\n",
	    tid, MaxMemSiz*1e-9, maxVmRSS*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, YYmap[Yid]->id,XXmap[Xid]->id,N,M,waittime,PausedThreads,Mem->lastGetMem,Mem->lastGetStatm,wt);
	  fflush(stdout);
	}
      }
    }

    long long myVmRSSreserved = 0;
    if(MEMRESERVE){
      if(THREAD_DEBUG) Mem->atomic_calls++;

      #pragma omp atomic read
      myVmRSSreserved = VmRSSreserved;
    }

    if(!madvise_called && Mem->MemUsed > 0 && (Mem->VmRSS + Mem->VmSwap + myVmRSSreserved + RSSneeded > maxVmRSS || 
					       (Mem->VmRSS + Mem->VmSwap + myVmRSSreserved + RSSneeded > maxVmRSS * RA_MAX_MEM && wt >= Mem->lasttime + RA_MIN_TIME))){
      long long sVmSize,sVmRSS,sVmSwap,sVmHWM;// local copies, so global shared values are not modified
      if(VERB>=3)
	getmem(sVmSize,sVmRSS,sVmSwap,&sVmHWM);

      if(THREAD_DEBUG) Mem->madvise_calls++;

      if(madvise(Fmem, Mem->MemSiz /* HERE HERE Mem->MemUsed */, MADV_DONTNEED)){
        int eno = errno;
        char *err = strerror(eno);
        printf("madvise(%p,%lld,MADV_DONTNEED) failed MemUsed=%lld:errno=%d:%s\n",Fmem,Mem->MemSiz,Mem->MemUsed,eno,err);
        fflush(stdout);exit(1);
      }
      madvise_called = true;
      Mem->lasttime = wt;

      double origwt = wt;
      wt = wtime();

      if(VERB>=3 && PausedThreads > 0){
        long long rVmSize,rVmRSS,rVmSwap,rVmHWM;
        getmem(rVmSize,rVmRSS,rVmSwap,&rVmHWM);

        #pragma omp critical(PausedThreads)
        {
	  if(PausedThreads > 0 || nMemSiz > 1000000000LL){
	    printf("tid=%d: After madvise on %0.4f Gb : VmRSS= %0.4f -> %0.4f (HWM=%0.2f), VmSwap= %0.4f, Lim= %0.4f, Used=%0.4f, Reserved=%0.4f Gb, : (Yid=%lld,Xid=%lld,N=%d,M=%d):Paused=%d:wall=%0.6f(elapsed=%0.6f)\n",
	      tid, Mem->MemUsed*1e-9,sVmRSS*1e-9,rVmRSS*1e-9, rVmHWM*1e-9, rVmSwap*1e-9, maxVmRSS*1e-9,Mem->MemUsed*1e-9,myVmRSSreserved*1e-9,YYmap[Yid]->id,XXmap[Xid]->id,N,M,PausedThreads,wt-origwt,wt);
	    fflush(stdout);
          }
        }
      }

      Mem->MemUsed = 0;
      RSSneeded = nMemSiz * RSSratio;

      /* update VmRSS again : this should not happen that often, since madvice is not called often unless running out of memory */
      if(THREAD_DEBUG) Mem->getstatm_calls++;
      Mem->GetStatm(wt);
      if(VERB>=2 && !ThreadNotPaused){
        #pragma omp critical(PausedThreads)
	{
	  if(PausedThreads > 0){
	    printf("tid=%d: MaxMem= %0.1f (Lim=%0.2f): GetStatm: VmRSS= %0.4f, VmSwap= %0.2f, nMem=%0.4f Gb (Yid=%lld,Xid=%lld,N=%d,M=%d),waittime= %0.3f, Paused Threads=%d:wt=%0.6f,%0.6f,%0.6f\n",
	      tid, MaxMemSiz*1e-9, maxVmRSS*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, nMemSiz*1e-9, YYmap[Yid]->id,XXmap[Xid]->id,N,M,waittime,PausedThreads,Mem->lastGetMem,Mem->lastGetStatm,wt);
	    fflush(stdout);
          }
	}
      }

      if(MEMRESERVE){
	if(THREAD_DEBUG) Mem->atomic_calls++;

        #pragma omp atomic read
        myVmRSSreserved = VmRSSreserved;
      }
    }

    int myPausedThreads;
    if(THREAD_DEBUG) Mem->atomic_calls++;

    #pragma omp atomic read
    myPausedThreads = PausedThreads;

    long long origMemUsed = Mem->MemUsed;
    if(Mem->VmRSS + Mem->VmSwap + myVmRSSreserved + RSSneeded <= maxVmRSS || myPausedThreads + ThreadNotPaused >= NumThreads || (VMEMCHECK && vmemTotsiz <= 0)){
      Mem->MemUsed = max(Mem->MemUsed, nMemSiz);
      if(MEMRESERVE){
        if(DEBUG/* HERE >= 2*/) assert(Mem->MemReserved == 0);
        Mem->MemReserved = (Mem->MemUsed - origMemUsed) * RSSratio;
	
	/* NOTE : "omp atomic" is faster than "omp critical" but risks allowing too many threads to reserve the same memory, hence must recheck VmRSS after new memory is initialized  */

	#pragma omp atomic capture
	myVmRSSreserved = VmRSSreserved += Mem->MemReserved;

	if(THREAD_DEBUG) Mem->atomic_calls++;
      }
      memfound = true;

      if(VMEMCHECK){
	if(DEBUG>=2) assert(vmemTotsiz >= 0);
	vmemTotsiz += nMemSiz;
      }
      if(!ThreadNotPaused){
	ThreadNotPaused = 1;

	if(THREAD_DEBUG) Mem->atomic_calls++;
        #pragma omp atomic capture
	myPausedThreads = --PausedThreads;
      }
    } else if(ThreadNotPaused){
      if(DEBUG>=2) assert(!memfound);
      ThreadNotPaused = 0;

      if(THREAD_DEBUG) Mem->atomic_calls++;
      #pragma omp atomic capture
      myPausedThreads = ++PausedThreads;
    }

    if(memfound){
      if(DEBUG && waittime > 0.0) assert(origMemUsed == 0);

      if((VERB && (waittime > 10.001 /* || nMemSiz > VMRSS_PEAKCHECK */)) || (VMEMCHECK && vmemTotsiz < nMemSiz)){
        if(THREAD_DEBUG) Mem->critical_calls++;

        #pragma omp critical(PausedThreads)
	{
	  printf("tid=%d: MaxMem= %0.1f (Lim=%0.2f), VmRSS= %0.4f, VmSwap= %0.2f, Used= %0.4f->%0.4f, Reserved= %0.4f->%0.4f Gb: waited %0.2f secs for %0.4f Gb (Yid=%lld,Xid=%lld,N=%d,M=%d),Paused Threads=%d:wt=%0.6f,%0.6f,%0.6f\n",
	         tid, MaxMemSiz*1e-9, maxVmRSS*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, origMemUsed*1e-9,Mem->MemUsed*1e-9,(myVmRSSreserved - Mem->MemReserved)*1e-9, myVmRSSreserved*1e-9,
		 waittime, RSSneeded*1e-9, YYmap[Yid]->id,XXmap[Xid]->id,N,M,TSAN ? myPausedThreads : PausedThreads,Mem->lastGetMem,Mem->lastGetStatm,wt);
	  if(VERB/* HERE >=3 */){
	    long long VmSize,VmRSS,VmSwap,VmHWM;// local copies, so global shared values are not modified
	    getmem(VmSize,VmRSS,VmSwap,&VmHWM);
	    printf("\t Actual VmRSS= %0.4f(HWM= %0.4f), VmSwap= %0.4f, VmSize= %0.4f Gb:wall=%0.6f\n",VmRSS*1e-9,VmHWM*1e-9,VmSwap*1e-9,VmSize*1e-9,wtime());
          }
	  fflush(stdout);
        }
      }
      break;
    }

    if(DEBUG) assert(MaxMemSiz > 0);

    if(!madvise_called && Mem->MemUsed > 0){
      long long sVmSize,sVmRSS,sVmSwap,sVmHWM;// local copies, so global shared values are not modified
      if(VERB>=3)
	getmem(sVmSize,sVmRSS,sVmSwap,&sVmHWM);

      if(THREAD_DEBUG) Mem->madvise_calls++;

      if(madvise(Fmem, Mem->MemSiz /* HERE HERE Mem->MemUsed */, MADV_DONTNEED)){
        int eno = errno;
        char *err = strerror(eno);
        printf("madvise(%p,%lld,MADV_DONTNEED) failed MemUsed=%lld:errno=%d:%s\n",Fmem,Mem->MemSiz,Mem->MemUsed,eno,err);
        fflush(stdout);exit(1);
      }
      madvise_called = true;
      Mem->lasttime = wt;

      if(VERB>=3){
        double origwt = wt;
        wt = wtime();
	
        long long rVmSize,rVmRSS,rVmSwap,rVmHWM;// local copies, so global shared values are not modified
        getmem(rVmSize,rVmRSS,rVmSwap,&rVmHWM);

        #pragma omp critical
        {
          printf("tid=%d: After madvise on %0.4f Gb : VmRSS= %0.4f -> %0.4f(HWM=%0.2f), VmSwap= %0.4f Gb: (Yid=%lld,Xid=%lld,N=%d,M=%d):wall=%0.6f(elapsed=%0.6f)\n",
		 tid, Mem->MemUsed*1e-9,sVmRSS*1e-9,rVmRSS*1e-9, rVmHWM*1e-9, rVmSwap*1e-9, YYmap[Yid]->id,XXmap[Xid]->id,N,M,wt-origwt,wt);
	  fflush(stdout);
        }
      }

      Mem->MemUsed = 0;
      RSSneeded = nMemSiz * RSSratio;
    }

    if(DEBUG && MEMRESERVE) assert(Mem->MemReserved == 0 && Mem->MemUsed == 0);

    if(VERB/* HERE >=2 */ && (waitcnt % TICKS60) == RA_TICKS){
      if(THREAD_DEBUG) Mem->critical_calls++;

      #pragma omp critical(PausedThreads)
      {
	printf("tid=%d: MaxMem= %0.1f (Lim=%0.2f), VmRSS= %0.4f, VmSwap= %0.2f, Reserved= %0.4f Gb: waiting %0.2f secs for %0.4f/%0.4f Gb (Yid=%lld,Xid=%lld,N=%d,M=%d),Paused Threads=%d:wt=%0.6f,%0.6f,%0.6f\n",
	   tid, MaxMemSiz*1e-9, maxVmRSS*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, myVmRSSreserved*1e-9, waittime, RSSneeded*1e-9, nMemSiz*1e-9,
	YYmap[Yid]->id,XXmap[Xid]->id,N,M,TSAN ? myPausedThreads : PausedThreads,Mem->lastGetMem,Mem->lastGetStatm,wt);
	if(THREAD_DEBUG)
	  printf("\t wait=%lu,release=%lu,madvise=%lu:getmem=%lu,getstatm=%lu,critical=%lu,atomic=%lu\n",Mem->wait_calls,Mem->release_calls,Mem->madvise_calls,
		 Mem->getmem_calls,Mem->getstatm_calls, Mem->critical_calls, Mem->atomic_calls);
	if(VERB/* HERE >=3 */){
	  long long VmSize,VmRSS,VmSwap,VmHWM;// local copies, so global shared values are not modified
	  getmem(VmSize,VmRSS,VmSwap,&VmHWM);
	  printf("\t Actual VmRSS= %0.4f(HWM= %0.4f), VmSwap= %0.4f Gb:wall=%0.6f\n",VmRSS*1e-9,VmHWM*1e-9,VmSwap*1e-9,wtime());
	}
	fflush(stdout);
      }
    }

    // sleep VMRSS_SLEEP seconds 
    struct timespec rem;
    nanosleep(&tm,&rem);

    waitcnt++;
    waittime += VMRSS_SLEEP;
  } // while(!memfound)
}

/** The main pairwise Dynamic Programming routine. 
Hardware accellerations done with USE_MIC and USE_AVX preprocessor directives.
*/
void pairalignXY(FLOAT *Y, int N,
		 FLOAT *X, int M,
		 AL &A,/* includes preallocated memory JMIN[1..maxN],JMAX[1..maxN], Imin[1..maxM],Imax[1..maxM], for use with hashtable offset information for diagonalized arrays */
		 Calign *align,/* Best alignment (with score) is stored here */
		 int Yid,int Xid,
		 PFLOAT **delY,/* delY[I=2..N][n=max(0,DELTA+1-I)..DELTA-1] == Y[I] - Y[I-DELTA+n] */
		 PFLOAT **delYx, /* delY but aligned like delX */
		 PFLOAT **delX,/* delX[J=2..M][m=max(0,DELTA+1-J)..DELTA-1] == X[J] - X[J-DELTA+m] 
				  Except If (USE_SSE && USE_AVX || USE_MIC) && USE_PFLOAT : delX[n=1..DELTA][J=1+n..M] = X[J] - X[J-n] */
		 Chashpair *phash,/* If != 0 pointer to hashtable entry : Not currently used */
		 int orientation, int scaleID,
                 int myoutlierExtend, int myoutlierBC,
                 int tid,
	         PFLOAT * &Fmem, // preallocated memory Fmem == Mem->Fmem[0 .. Mem->MemSiz - 1] to be used by A.score[][],A.Uscore[][] etc */
  	         CThreadMem *Mem,
	         long long &StrideMax /* If (DIAGONAL>=2 && MIN_MEM && hashdelta) : the current allocated number of elements per array */
		 )
{
  PFLOAT C = ChimScore;
  int localtype = (PoutlierEnd > 0.0) ? -3 : -2;

  int Merge = 1;/* Will be set to 0 if PairMerge AND best alignment (in current orientation) has too large internal or external outliers for merging (see PAIRMERGE_FIX) */

  int XleftClipped = 0, XrightClipped = 0, YleftClipped = 0, YrightClipped = 0;
  if(OUTLIER_FIXES && PairSplit){/* check if any map ends are Clipped (truncated) ends */
    Cmap *Ymap = YYmap[Yid];
    Cmap *Xmap = XXmap[Xid];
    if(Ymap->origmap){
      if(Ymap->left[0] > 0)
	YleftClipped = 1;
      if(Ymap->right[0] < Ymap->origmap->numsite[0])
	YrightClipped = 1;
    }
    if(Xmap->origmap){
      if(Xmap->left[0] > 0)
	(orientation ? XrightClipped : XleftClipped) = 1;
      if(Xmap->right[0] < Xmap->origmap->numsite[0])
	(orientation ? XleftClipped : XrightClipped) = 1;
    }
  }

  if(EVERB){
    //    #pragma omp critical
    {
      if(phash)
	printf("pairalignXY:Yid=%d(%lld),Xid=%d(%lld),or=%d,sc=%d,N=%d,M=%d:phash:id1=%d,id2=%d,or=%d,sc=%d,score=%d,offset=%d,myoutlierExtend=%d,myoutlierBC=%d,tid=%d:wtime=%0.6f\n",
	       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,N,M,phash->id1,phash->id2,phash->orientation,phash->scaleID,phash->hashscore,phash->offset,myoutlierExtend,myoutlierBC,tid,wtime());
      else
	printf("pairalignXY:Yid=%d(%lld),Xid=%d(%lld),or=%d,sc=%d,N=%d,M=%d:XleftClipped=%d,XrightClipped=%d,YleftClipped=%d,YrightClipped=%d,C=%0.6f, myoutlierExtend=%d,myoutlierBC=%d,tid=%d:wtime=%0.6f\n",
	       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,N,M,XleftClipped,XrightClipped,YleftClipped,YrightClipped,C, myoutlierExtend,myoutlierBC,tid,wtime());
      fflush(stdout);
    }
  }

  int IMIN = 1, IMAX = N;
  int jmin = 1, jmax = M;
  int *JMIN = A.JMIN, *JMAX = A.JMAX;
  int *Imin = A.Imin, *Imax = A.Imax;

  int Jmin1 = jmin, Jmax1 = jmax;/* only used if DIAGONAL >= 2 and is the allocated J range for index I==IMIN.
				    Allocated J range for other I are:
				    If (Jmax1-Jmin1+1 <  jmax-jmin+1) : Jmin1+(I-IMIN) ... Jmax1+(I-IMIN)
				    If (Jmax1-Jmin1+1 >= jmax-jmin+1) : Jmin1 ... Jmax1 
				 */

  if(phash && hashdelta){
    double offset = phash->offset;

    if(!(OFFSET_SCALE_FIX && hashScaleDelta) && scaleID != 0 ){    /* need to slightly correct offset if scaleID != 0 */
      FLOAT scale = ScaleFactor[scaleID];
      offset -= 0.5 * X[M+1] * (1.0 - 1.0 / scale);
    }

    double deltaoffset = Ylambda * DELTA;
    double offsetSD = max(deltaoffset, hashdelta * 2.0 * OffsetKB);
    int ioffsetSD = floor(hashdelta + 0.9999);
    ioffsetSD = max(ioffsetSD, max(DELTA_X,DELTA_Y));

    int startjmin = jmin;
    int startjmax = jmax;
    if(ShiftOffset1 <= 30){  /* update jmin and jmax based on phash->MinLoc2 and phash->MaxLoc2 */
      for(;jmin < jmax; jmin++)
	if(X[jmin+1] + offsetSD/* NEW78 */ >= phash->MinLoc2)
	  break;
      for(;jmin < jmax; jmax--)
	if(X[jmax-1] - offsetSD <= phash->MaxLoc2)
	  break;

      Jmin1 = jmin;
      Jmax1 = jmax;
    }

    if(HASHRANGE){
      int origIMIN = IMIN, origIMAX = IMAX;
      while(IMIN < IMAX && Y[IMIN+1]-X[jmin] < offset)
	IMIN++;
      while(IMIN < IMAX && Y[IMAX-1]-X[jmax] > offset)
	IMAX--;
      IMIN = max(origIMIN/* WAS 1 */, IMIN - ioffsetSD);
      IMAX = min(origIMAX/* WAS N */, IMAX + ioffsetSD);
    } else {
      while(IMIN < IMAX && Y[IMIN+1]-X[jmin] < offset - offsetSD)
	IMIN++;
      while(IMIN < IMAX && Y[IMAX-1]-X[jmax] > offset + offsetSD)
	IMAX--;
    }

    for(int I = IMIN; I <= IMAX; I++){
      JMIN[I] = jmin;
      JMAX[I] = jmax;
    }

    if(DIAGONAL){
      for(int I = IMIN; I <= IMAX; I++){
	if(HASHRANGE){
	  while(JMIN[I] < JMAX[I] && Y[I] - X[JMIN[I]+1] > offset)
	    JMIN[I]++;
	  while(JMIN[I] < JMAX[I] && Y[I] - X[JMAX[I]-1] < offset)
	    JMAX[I]--;
	  JMIN[I] = max(jmin, JMIN[I] - ioffsetSD);
	  JMAX[I] = min(jmax, JMAX[I] + ioffsetSD);
	} else {
	  while(JMIN[I] < JMAX[I] && Y[I] - X[JMIN[I]+1] > offset + offsetSD)
	    JMIN[I]++;
	  while(JMIN[I] < JMAX[I] && Y[I] - X[JMAX[I]-1] < offset - offsetSD)
	    JMAX[I]--;
	}
      }
      
      long long jcnt1 = 0, jcnt2 = 0;      
      if(MIN_MEM_DEBUG)
	for(int I = IMIN; I <= IMAX; I++)
	  jcnt1 += JMAX[I] - JMIN[I] + 1;

      /* enforce monotonic restrictions IMIN .. IMAX :
	 1. JMIN[I=IMAX..IMIN] is a non-increasing function as I decreases.
	 2. JMAX[I=IMIN..IMAX] is a non-decreasing function as I increase.
      */
      for(register int I = IMAX; --I >= IMIN;)
	JMIN[I] = min(JMIN[I],JMIN[I+1]);
      for(register int I = IMIN; ++I <= IMAX;)
	JMAX[I] = max(JMAX[I],JMAX[I-1]);

      if(MIN_MEM_DEBUG || RA_MIN_TIME)
	for(int I = IMIN; I <= IMAX; I++)
	  jcnt2 += JMAX[I] - JMIN[I] + 1;

      long long origStrideMax = StrideMax;

      if(DIAGONAL>=2){
	/* Re-allocate A->field[I][J] based on JMIN[],JMAX[] for more compact real memory usage */
	printf("DIAGONAL >= 2 : not yet implemented\n");
	fflush(stdout);exit(1);
      }
    }
  } else {// no hashtable information used */
    for(int I = IMIN; I <= IMAX; I++){
      JMIN[I] = jmin;
      JMAX[I] = jmax;
    }
  }

  if(DEBUG>=2){
    for(int I = IMIN; I <= IMAX; I++){
      assert(1 <= JMIN[I] && JMIN[I] <= JMAX[I] && JMAX[I] <= M);
    }
  }

  /* compute Imin[],Imax[] from IMIN,IMAX,JMIN[],JMAX[] : same code as at end of setlimit() */
  for(int J = 1; J <= M; J++){
    Imin[J] = N+1;
    Imax[J] = 0;
  }
  for(int I = IMIN; I <= IMAX; I++){
    for(int J = JMIN[I]; J <= JMAX[I]; J++){
      Imin[J] = min(Imin[J],I);
      Imax[J] = max(Imax[J],I);
    }
  }

  if(DEBUG>=2)
    for(register int J = 1; J <= M; J++)
      for(register int I = Imin[J]; I <= Imax[J]; I++)
	assert(JMIN[I] <= J && J <= JMAX[I]);

  for(int J = 1; J < JMIN[IMIN]; J++){
    Imin[J] = 1;
    Imax[J] = 0;
  }
  for(int J = JMAX[IMAX]; ++J <= M;){
    Imin[J] = N+1;
    Imax[J] = N;
  }

  for(int I = 1; I < IMIN; I++){
    JMIN[I] = 1;
    JMAX[I] = 0;
  }
  for(int I = IMAX; ++I <= N; ){
    JMIN[I] = M+1;
    JMAX[I] = M;
  }

  if(DEBUG>=2 && !phash){
    for(int I = IMIN; I <= IMAX; I++){
      assert(JMIN[I] == 1);
      assert(JMAX[I] == M);
    }    
    for(int J = 1; J <= M; J++){
      assert(Imin[J] == IMIN);
      assert(Imax[J] == IMAX);
    }
  }

  size_t Fsiz = A.array_stride * DPsizF + 16LL;
  size_t Isiz = A.array_stride * DPsizI + 16LL;
  //  size_t Fsiz = N * ((M + 64LL) + 16LL) * DPsizF + 16LL;
  //  size_t Isiz = N * ((M + 64LL) + 16LL) * DPsizI + 16LL;
  long long nMemSiz = ((Fsiz * sizeof(PFLOAT) + Isiz * sizeof(int)) + (PAGE-1)) & ~(PAGE-1);

  /* Compute expected ratio of Resident Memory to allocated memory (nMemSiz)*/
  double RSSratio = ENDINIT ? 0.125 + /*NEW29*/0.25*(min(1.0,2048.0/N) + min(1.0,2048.0/M)) : 0.625;/* score + (if ENDINIT==0) Rijy,Lijy,Rijx,Lijx fully used */
  if((DEFER_BP && !myoutlierExtend) || RepeatRec)
    RSSratio += 0.125;/* Uscore fully used */
  if(!(DEFER_BP && !myoutlierExtend && !myoutlierBC && DELTA==4))/* non-vectorized code */
    RSSratio += 0.25;/* G,H fully used */
  if(DEBUG) assert(RSSratio <= 1.0);

  int retries = 0;/* number of times we jumped to LwaitForMemory */

 LwaitForMemory:
  if(MEMCHECK && RA_MIN_TIME && MaxMemSiz > 0 && RA_MIN_MEM < 8)
    WaitForMemory(Mem,nMemSiz,tid,Yid,Xid,N,M,Fmem, RSSratio);    

  int *Ilist, *Jlist, *Ilist2, *Jlist2;
  FLOAT *iscore, *outscore;
  PFLOAT *Yf,*Xf,*Yr,*Xr;
  char *mem;
  if(retries == 0){
    mem = NULL;

    if(M < MAX_ALLOCA && N < MAX_ALLOCA){
      /* allocate space to store alignment in reverse order */
      Ilist = (int *)alloca(M*sizeof(int));
      Jlist = (int *)alloca(M*sizeof(int));
      iscore = (FLOAT *)alloca((M+1)*sizeof(FLOAT));
      outscore = (FLOAT *)alloca((M+1)*sizeof(FLOAT));

      /* allocate additional space for RepeatRec computation */
      Ilist2 = (int *)alloca(M*sizeof(int));
      Jlist2 = (int *)alloca(M*sizeof(int));

      /* allocate space for Maps using PFLOAT instead of double : only used near ends of maps to avoid loss of precision */
      Yf = (PFLOAT*)alloca(sizeof(PFLOAT)*(N+1));/* precompute Yf[I] = (PFLOAT)Y[I] */
      Xf = (PFLOAT*)alloca(sizeof(PFLOAT)*(M+1));/* precompute Xf[J] = (PFLOAT)[J] */
      Yr = (PFLOAT*)alloca(sizeof(PFLOAT)*(N+1));/* precompute Yr[I] = (PFLOAT)(Y[N+1]-Y[I]) */
      Xr = (PFLOAT*)alloca(sizeof(PFLOAT)*(M+1));/* precompute Xr[J] = (PFLOAT)(X[M+1]-X[J]) */
    } else {
      size_t siz = 4ul * M * sizeof(int) + 2ul * (M+1) * sizeof(FLOAT) + 2ul * (N+M+2) *sizeof(PFLOAT);
      mem = (char *)malloc(siz);  // HERE HERE : pre-allocate this memory as part of Fmem + Imem
      if(!mem){
	printf("pairalignXY:malloc(%lu) failed (M=%d,N=%d)\n",siz,M,N);
	fflush(stdout);exit(1);
      }

      size_t cnt = 0;

      Ilist = (int *)&mem[cnt]; cnt += M*sizeof(int);
      Jlist = (int *)&mem[cnt]; cnt += M*sizeof(int);
      iscore = (FLOAT *)&mem[cnt]; cnt += (M+1)*sizeof(FLOAT);
      outscore = (FLOAT *)&mem[cnt]; cnt += (M+1)*sizeof(FLOAT);

      Ilist2 = (int *)&mem[cnt]; cnt += M*sizeof(int);
      Jlist2 = (int *)&mem[cnt]; cnt += M*sizeof(int);

      Yf = (PFLOAT *)&mem[cnt]; cnt += (N+1)*sizeof(PFLOAT);
      Xf = (PFLOAT *)&mem[cnt]; cnt += (M+1)*sizeof(PFLOAT);
      Yr = (PFLOAT *)&mem[cnt]; cnt += (N+1)*sizeof(PFLOAT);
      Xr = (PFLOAT *)&mem[cnt]; cnt += (M+1)*sizeof(PFLOAT);

      if(DEBUG) assert(cnt == siz);
    }

    for(int I = 1; I <= N; I++){
      Yf[I] = Y[I];
      Yr[I] = Y[N+1] - Y[I];
    }
    for(int J = 1; J <= M; J++){
      Xf[J] = X[J];
      Xr[J] = X[M+1] - X[J];
    }
  }

  if(DEBUG>=2){/* initialize all memory with invalid values */
    for(int I = IMIN; I <= IMAX; I++){
      for(int J = JMIN[I]; J <= JMAX[I]; J++){
	FLAT_A(A.score,I,J) = aNaN;
	FLAT_A(A.Uscore,I,J) = aNaN;
	FLAT_A(A.G,I,J) = -555;
	FLAT_A(A.H,I,J) = -555;
	FLAT_A(A.Rijy,I,J) = -555;
	FLAT_A(A.Lijy,I,J) = -555;
	FLAT_A(A.Rijx,I,J) = -555;
	FLAT_A(A.Lijx,I,J) = -555;
      }
    }
  }

  /* initialize A.score[I=1..N][J=1..M] = C (typically MINSCORE) and G = localtype */
  if(RepeatShift > 0.0){
    for(int I = IMIN; I <= IMAX; I++){
      int J, jmax = min(JMAX[I], I-1);
      for(J = JMIN[I]; J <= jmax; J++){
	if(Y[I] - X[J] < RepeatShift)
	  break;
	FLAT_A(A.score, I, J) = C;
      }
      for(; J <= JMAX[I]; J++)
	FLAT_A(A.score, I, J) = MINSCORE;
    }
  } else {
    for(int I = IMIN; I <= IMAX; I++)
      memset_float(&FLAT_A(A.score,I,JMIN[I]), C, JMAX[I]-JMIN[I]+1);
    /*      for(int J= JMIN[I]; J <= JMAX[I]; J++)
	    FLAT_A(A.score, I, J) = C;*/

    if(EVERB && 0 < I_TRACE && I_TRACE <= N && 0 < J_TRACE && J_TRACE <= M){
      printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d, A.score(I,J)= %0.6f, ChimScore= %0.6f\n",Yid,Xid,orientation,scaleID,I_TRACE,J_TRACE,FLAT_A(A.score,I_TRACE,J_TRACE),ChimScore);
      fflush(stdout);
    }
  }

  if(!DEFER_BP || myoutlierExtend){
    for(int I = IMIN; I <= IMAX; I++)
      memset_int(&FLAT_A(A.G,I,JMIN[I]), localtype, JMAX[I]-JMIN[I]+1);
    /*    for(int J = 1; J <= M; J++)
      FLAT_A(A.G, 1, J) = localtype;
    for(int I = 2; I <= N; I++)
    FLAT_A(A.G, I, 1) = localtype;*/

#ifdef VALGRIND
    for(int I = IMIN; I <= IMAX; I++)
      memset_int(&FLAT_A(A.H,I,JMIN[I]), localtype, JMAX[I]-JMIN[I]+1);
    /*    for(int J = 1; J <= M; J++)
      FLAT_A(A.H, 1, J) = localtype;
    for(int I = 2; I <= N; I++)
    FLAT_A(A.H, I, 1) = localtype;*/
#endif
  }

  double wt;
  if(MEMCHECK>=3 && !MEMRESERVE && RA_MIN_TIME && MaxMemSiz > 0 && RA_MIN_MEM < 8 && retries < MEMCHECK-1 && 
     ((wt = wtime()) > Mem->lastGetStatm + VMRSS_CHECK && (wt > Mem->lastGetStatm + VMRSS_CHECK2 || Mem->VmRSS + Mem->VmSwap >= MaxMemSiz2 * RA_MAX_MEM))){/* recheck VmRSS */
    //    long long maxVmRSS2 = MaxMemSiz2;// maximum VmRSS+VmSwap allowed for large memory alignments (must fit within 90% of memory)
    //    long long maxVmRSS1 = MaxMemSiz;// maxmimum VmRSS+VmSwap allowed for jobs expected to use less than 10% of the per thread amount (so total will fit within 100% of memory)
    //    bool SmallJob = nMemSiz < (MaxMemSiz / (NumThreads * 10.0));
    //    long long maxVmRSS = SmallJob ? maxVmRSS1 : maxVmRSS2;// maximum VmRSS+VmSwap allowed for current pairwise alignment

    if(THREAD_DEBUG) Mem->getstatm_calls++;
    Mem->GetStatm(wt);

    if(Mem->VmRSS + Mem->VmSwap > MaxMemSiz /* WAS108 maxVmRSS */){
      if(VERB>=2){
        #pragma omp critical(PausedThreads)
	{
	  if(PausedThreads > 0 || (VMRSS_PEAKCHECK && Mem->MemUsed > VMRSS_PEAKCHECK)){
	    printf("tid=%d: MaxMem= %0.1f: GetStatm: VmRSS= %0.4f, VmSwap= %0.2f, MemUsed= %0.4f(nMem=%0.4f) (Yid=%lld,Xid=%lld,N=%d,M=%d),retries=%d(A), Paused Threads=%d:wt=%0.6f,%0.6f\n",
		   tid, MaxMemSiz*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, Mem->MemUsed*1e-9, nMemSiz*1e-9,YYmap[Yid]->id,XXmap[Xid]->id,N,M,retries,PausedThreads,Mem->lastGetMem,Mem->lastGetStatm);
	    fflush(stdout);
          }
	}
      }
      retries++;
      
      if(MEMRESERVE){
	if(THREAD_DEBUG) Mem->atomic_calls++;

        #pragma omp atomic update
	VmRSSreserved -= Mem->MemReserved;
      
	Mem->MemReserved = 0;
      }

      goto LwaitForMemory;
    }
  }

  /* initialize A.Lijy[I=1..N][J=1..M] as leftmost site in Y that overlaps X with Y[I] aligned with X[J] */
  /* If ENDINIT : only cases where I <= DELTA OR J <= DELTA are initialized */

  /* first estimate LijyJmax, the smallest value J = 1..M so that X[J] >= Y[imax(J)] OR M+1 if no such value J exists : use binary search to find this value between 1 and M */
  /* NOTE : imax(J) == (J <= DELTA) ? N : min(DELTA,N) */
  /* NOTE : for all J >= LijyJmax, Lijy == 0. Conversly for all J < LijyJmax, Lijx > 0 for I == imax(J) (but possibly NOT for any I in range IMIN..IMAX) */
  int LijyJmax = M+1;
  if(X[M] >= Y[!ENDINIT ? N : (M <= DELTA) ? N : min(DELTA,N)]){/* LijyJmax <= M */
    if(X[1] >= Y[!ENDINIT ? N : (1 <= DELTA) ? N : min(DELTA,N)])/* LijyJmax <= 1 */
      LijyJmax = 1;
    else {
      int lb = 1, ub = M;      
      while(lb + 1 < ub){
	int J = (lb + ub)/2;
	int imax = !ENDINIT ? N : (J <= DELTA) ? N : min(DELTA,N);
	if(X[J] >= Y[imax])
	  ub = J;
	else
	  lb = J;
      }
      if(DEBUG) assert(lb + 1 == ub);
      if(DEBUG>=2) assert(X[ub] >= Y[!ENDINIT ? N : (ub <= DELTA) ? N : min(DELTA,N)]);
      if(DEBUG>=2) assert(X[lb] < Y[!ENDINIT ? N : (lb <= DELTA) ? N : min(DELTA,N)]);
      LijyJmax = ub;
    }
  }
  // NOTE : LijyJmax is typically small and close to DELTA+1
  
  for(int J=1; J < LijyJmax; J++){
    int imax = !ENDINIT ? Imax[J] : (J <= DELTA) ? Imax[J] : min(DELTA,Imax[J]);
    
    int Lijy = 0;
    for(int I= Imin[J]; I <= imax; I++){ /* update Lijy by incrementing it until Y[Lijy] overlaps X */
      FLOAT delta = Y[I] - X[J];
      while(delta > Y[Lijy])
	Lijy++;
      if(DEBUG>=2) assert(Lijy <= I);
      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
      FLAT_A(A.Lijy, I, J) = Lijy;
    }
    if(DEBUG>=2 && !(phash && hashdelta) && !(Lijy > 0)){
      #pragma omp critical
      {
        printf("tid=%d: Yid=%lld,Xid=%lld,N=%d,M=%d: J=%d, Lijy=%d (LijyJmax=%d, DELTA=%d, imax=%d)\n",tid,YYmap[Yid]->id,XXmap[Xid]->id,N,M,J,Lijy,LijyJmax,DELTA, imax);
	printf("\t X[J]= %0.4f, Y[imax]= %0.4f\n",X[J],Y[imax]);
	fflush(stdout);
	assert(Lijy > 0);
      }
    }
  }

  int DELTA_M = min(M,DELTA);
  if(LijyJmax <= DELTA_M){
    for(int I= IMIN; I <= IMAX; I++){
      int Jmax = min(DELTA_M, JMAX[I]);
      for(int J = max(JMIN[I], LijyJmax); J <= Jmax; J++){
	if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
	FLAT_A(A.Lijy, I, J) = 0;		
      }
    }
  }

  int J = max(LijyJmax, DELTA_M + 1);
  if(J <= M){
    int imax = !ENDINIT ? IMAX : min(DELTA,IMAX);

    for(int I= IMIN; I <= imax; I++){
      int jmin = max(J, JMIN[I]);
      int jmax = JMAX[I];
      if(jmin <= jmax)
	memset_int(&FLAT_A(A.Lijy, I, jmin), 0, jmax-jmin+1);
    }
  }

  /* initialize A.Lijx.[I=1..N][J=1..M] as leftmost site in X that overlap Y with Y[I] aligned with X[J] */
  /* If ENDININT : only cases where I <= DELTA OR J <= DELTA are initialized */
  
  /* first estimate LijxImax, the smallest value I = 1..N so that X[Jmax(I)] <= Y[I] (or N+1 if no such I exists): use binary search to find this value between 1 and N */
  /* NOTE : Jmax(I) == (I <= DELTA) ? M : min(DELTA,M) */
  /* NOTE : for all I >= LijxImax, Lijx == 0. Conversely for all I < LijxImax, Lijx > 0 for some J == Jmax(I) (but possibly NOT for any J within JMIN[I] .. JMAX[I] */
  int LijxImax = N+1;
  if(Y[N] >= X[!ENDINIT ? M : (N <= DELTA) ? M : min(DELTA,M)]){/* LijxImax <= N */
    if(Y[1] >= X[!ENDINIT ? M : (1 <= DELTA) ? M : min(DELTA,M)])/* LijxImax <= 1 */
      LijxImax = 1;
    else {
      int lb = 1, ub = N;
      while(lb + 1 < ub){
        int I = (lb + ub)/2;
	int Jmax = !ENDINIT ? M : (I <= DELTA) ? M : min(DELTA,M);
	if(Y[I] >= X[Jmax])
	  ub = I;
	else
	  lb = I;
      }
      if(DEBUG) assert(lb + 1 == ub);
      if(DEBUG>=2) assert(Y[ub] >= X[!ENDINIT ? M : (ub <= DELTA) ? M : min(DELTA,M)]);
      if(DEBUG>=2) assert(Y[lb] < X[!ENDINIT ? M : (lb <= DELTA) ? M : min(DELTA,M)]);
      LijxImax = ub;
    }
  }

  // Note LijxImax is typically small and close to DELTA + 1
  for(int I= IMIN; I < LijxImax && I <= IMAX; I++){
    int jmax = !ENDINIT ? JMAX[I] : (I <= DELTA) ? JMAX[I] : min(DELTA,JMAX[I]);
    int Lijx = 0;

    for(int J= JMIN[I]; J <= jmax; J++){ /* update Lijx by incrementing it until X[Lijx] overlaps Y */
      FLOAT delta = X[J]-Y[I];
      while(delta > X[Lijx])
	Lijx++;
      if(DEBUG>=2) assert(Lijx <= J);
      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
      FLAT_A(A.Lijx, I, J) = Lijx;
    }
    if(DEBUG>=2 && !(phash && hashdelta)) assert(Lijx > 0);
  }

  int DELTA_N = min(IMAX,DELTA);
  if(LijxImax <= DELTA_N){
    for(int I = max(IMIN, LijxImax); I <= DELTA_N; I++)
      memset_int(&FLAT_A(A.Lijx, I, JMIN[I]), 0, JMAX[I]-JMIN[I]+1);
  }

  int I = max(max(IMIN,LijxImax), DELTA_N + 1);
  if(I <= IMAX){
    for(; I <= IMAX; I++){
      int jmax = !ENDINIT ? JMAX[I] : min(DELTA,JMAX[I]);
      if(jmax >= JMIN[I])
	memset_int(&FLAT_A(A.Lijx, I, JMIN[I]), 0, jmax - JMIN[I] + 1);
    }
  }

  if(DEBUG>=2){
    for(int I = IMIN; I <= IMAX; I++){
      int Jmax = !ENDINIT ? JMAX[I] : (I <= DELTA) ? JMAX[I] : min(DELTA,JMAX[I]);
      for(int J = JMIN[I]; J <= Jmax; J++){
	int Lijx = FLAT_A(A.Lijx, I, J);
	int Lijy = FLAT_A(A.Lijy, I, J);
	assert(Lijx <= 0 || Lijy <= 0);
      }
    }
  }

  if(MEMCHECK>=3 && !MEMRESERVE && RA_MIN_TIME && MaxMemSiz > 0 && RA_MIN_MEM < 8 && retries < MEMCHECK-1 && 
     ((wt = wtime()) > Mem->lastGetStatm + VMRSS_CHECK && (wt > Mem->lastGetStatm + VMRSS_CHECK2 || Mem->VmRSS + Mem->VmSwap >= MaxMemSiz2 * RA_MAX_MEM))){/* recheck VmRSS */
    //    long long maxVmRSS2 = MaxMemSiz2;// maximum VmRSS+VmSwap allowed for large memory alignments (must fit within 90% of memory)
    //    long long maxVmRSS1 = MaxMemSiz;// maxmimum VmRSS+VmSwap allowed for jobs expected to use less than 10% of the per thread amount (so total will fit within 100% of memory)
    //    bool SmallJob = nMemSiz < (MaxMemSiz / (NumThreads * 10.0));
    //    long long maxVmRSS = SmallJob ? maxVmRSS1 : maxVmRSS2;// maximum VmRSS+VmSwap allowed for current pairwise alignment

    if(THREAD_DEBUG) Mem->getstatm_calls++;
    Mem->GetStatm(wt);

    if(Mem->VmRSS + Mem->VmSwap > MaxMemSiz /* WAS108 maxVmRSS*/){
      if(VERB>=2){
        #pragma omp critical(PausedThreads)
	{
	  if(PausedThreads > 0 || (VMRSS_PEAKCHECK && Mem->MemUsed > VMRSS_PEAKCHECK)){
	    printf("tid=%d: MaxMem= %0.1f: GetStatm: VmRSS= %0.4f, VmSwap= %0.2f, MemUsed= %0.4f(nMem=%0.4f) (Yid=%lld,Xid=%lld,N=%d,M=%d),retries=%d(B), Reserved= %0.4f -> %0.4f Gb, PausedThreads=%d:wt=%0.6f,%0.6f\n",
		   tid, MaxMemSiz*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, Mem->MemUsed*1e-9, nMemSiz*1e-9,YYmap[Yid]->id,XXmap[Xid]->id,N,M,retries,
	           VmRSSreserved*1e-9, (VmRSSreserved - Mem->MemReserved)*1e-9, PausedThreads,Mem->lastGetMem,Mem->lastGetStatm);
	    fflush(stdout);
          }
	}
      }
      retries++;

      if(MEMRESERVE){
	if(THREAD_DEBUG) Mem->atomic_calls++;

        #pragma omp atomic update
	VmRSSreserved -= Mem->MemReserved;
      
	Mem->MemReserved = 0;
      }

      goto LwaitForMemory;
    }
  }

#if USE_MIC
  float *scoreSend = (float *)alloca((M+1)*sizeof(scoreSend));
  float *scoreSend1 = &scoreSend[-1];
#endif
  
  /* update A.score[I=1..N][J=1..M] with left ends unaligned (but NOT endoutlier) and Y[I] aligned with X[J] (along with G) */
  /* only cases where I <= DELTA OR J <= DELTA are initialized */
  for(int I= IMIN; I <= LijxImax && I <= IMAX; I++){

    int Jmax = (I <= DELTA) ? JMAX[I] : min(DELTA,JMAX[I]);
    if(RepeatShift > 0.0){
      Jmax = min(I-1,Jmax);
      while(Jmax >= 1 && Y[I]-X[Jmax] < RepeatShift)
	Jmax--;
    }

    int *LijxI = &FLAT_A(A.Lijx, I, 0), *LijyI = &FLAT_A(A.Lijy, I, 0);
    PFLOAT *scoreI = &FLAT_A(A.score, I, 0);

#if USE_MIC
    int score_size = Jmax;
    __mmask16 mask;
    __m512 vIdx=_mm512_add_ps(_mm512_set1_ps(I+2), v512_ascending0123);

    //prefetch_array(&Xf[1], score_size*sizeof(*Xf), _MM_HINT_T0);
    
    __m512 v_YfI=_mm512_set1_ps(Yf[I]), v1=_mm512_set1_ps(1.0f);
    __m512i v_limit=_mm512_sub_epi32(_mm512_set1_epi32(Jmax), v512i_ascending0123);
    
    /* make first pass of the loops to speedup computation */
    for(int J= JMIN[I]; J <= Jmax; J += 16){
      //      mask = (J < Jmax - 15) ? 0xffff : (0xffff>>(J+15-Jmax));
      mask= _mm512_cmp_epi32_mask(_mm512_set1_epi32(J), v_limit, _MM_CMPINT_LE);
      __m512i LijxIJ = _mm512_maskloadu_epi32(&LijxI[J], mask);
      __m512i LijyIJ = _mm512_maskloadu_epi32(&LijyI[J], mask);
      __m512 LijxIJd = _mm512_cvtfxpnt_round_adjustepi32_ps(LijxIJ, _MM_FROUND_CUR_DIRECTION, _MM_EXPADJ_NONE);
      __m512 LijyIJd = _mm512_cvtfxpnt_round_adjustepi32_ps(LijyIJ, _MM_FROUND_CUR_DIRECTION, _MM_EXPADJ_NONE);
      __m512 v_score = Send512(_mm512_gmin_ps(_mm512_maskloadu_ps(&Xf[J], mask),
	                                        v_YfI
	                                       ),
                                 _mm512_sub_ps(_mm512_add_ps(vIdx, _mm512_set1_ps(J)),
					       _mm512_add_ps(_mm512_gmax_ps(v1, LijxIJd), _mm512_gmax_ps(v1, LijyIJd))
					      )
			         );
      if(JMIN[I] == 1)
	_mm512_mask_store_ps(&scoreSend1[J], mask, v_score);
      else /* JMIN[I] > 1 : requires unaligned store for scoreSend1[J] */
	_mm512_maskstoreu_ps(&scoreSend1[J], mask, v_score);
    }
#endif
    
    int DELTA_Jmax = min(Jmax,LijyJmax);

    int J = JMIN[I];
    for(; J <= DELTA_Jmax; J++){
      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
      int Lijx = LijxI[J];
      int Lijy = LijyI[J];
      if(DEBUG>=2) assert(Lijx >= 0);
      if(DEBUG>=2) assert(Lijy >= 0);
      PFLOAT ascore;
      if(ENDFIX>=2){
	PFLOAT SmIJ = Sm(J,I);
	if(Lijx <= 0){/* left end of X fully overlapped by Y */
	  if(!XleftClipped){
#if USE_MIC
	    ascore = scoreSend1[J] + SmIJ;
#else
	    ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I+1-max(1,Lijy))) + SmIJ;
#endif
	    ascore = max(ascore,C);
	    for(int K= Lijy; K < I; K++) {
	      PFLOAT bound = Sbnd(Xf[J],Y[I]-Y[K],(PFLOAT)(J+I-K)) + SmIJ;
	      if(bound <= ascore)
		break;
	      ascore = bound;
	    }
	  } else
	    ascore = scoreI[J];
	} else {/* left end of Y fully overlapped by X */
	  if(DEBUG>=2) assert(Lijy == 0);
	  if(!YleftClipped){
#if USE_MIC
	    ascore = scoreSend1[J] + SmIJ;
#else
	    ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-max(1,Lijx)+I)) + SmIJ;
#endif
	    ascore = max(ascore,C);
	    for(int K= Lijx; K < J; K++){
	      PFLOAT bound = Sbnd(X[J]-X[K],Yf[I],(PFLOAT)(J-K+I)) + SmIJ;
	      if(bound <= ascore)
		break;
	      ascore = bound;
	    }
	  } else
	    ascore = scoreI[J];
	} 
      } else {
	if((Lijx <= 0) ? !XleftClipped : !YleftClipped){
#if USE_MIC
	  ascore = scoreSend1[J] + Sm(J,I);
#else
	  ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-max(1,Lijx)+I+1-max(1,Lijy))) + Sm(J,I);
#endif
	  ascore = max(ascore,C);
	} else
	  ascore = scoreI[J];
      }
      if(DEBUG>=2) assert(ascore >= C);
      scoreI[J] = ascore;
      if(ENDOUTLIER_FIX && ascore > C)
	FLAT_A(A.G, I, J) = -1;// NOT endoutlier
    }

    for(; J <= Jmax; J++){
      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
      int Lijx = LijxI[J];
      if(DEBUG>=2) assert(Lijx >= 0);
      if(Lijx > 0)
	break;
      if(DEBUG>=2) assert(Lijx == 0 && LijyI[J] > 0);

      /* left end of X fully overlapped by Y */
      PFLOAT ascore;
      if(!XleftClipped){
	PFLOAT SmIJ = Sm(J,I);
#if USE_MIC
	ascore = scoreSend1[J] + SmIJ;
#else
	ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I)) + SmIJ;
#endif
	ascore = max(ascore,C);
	if(ENDFIX>=2){
	  for(int K= LijyI[J]; K < I; K++){
	    PFLOAT bound = Sbnd(Xf[J],Y[I]-Y[K],(PFLOAT)(J+I-K)) + SmIJ;
	    if(bound <= ascore)
	      break;
	    ascore = bound;
	  }
	}
      } else
	ascore = scoreI[J];

      if(DEBUG>=2) assert(ascore >= C);
      scoreI[J] = ascore;
      if(ENDOUTLIER_FIX && ascore > C)
	FLAT_A(A.G, I, J) = -1;
    }
    for(; J <= Jmax; J++){
      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
      int Lijx = LijxI[J];
      if(DEBUG>=2) assert(Lijx > 0 && LijyI[J] == 0);
      /* left end of Y fully overlapped by X */
      PFLOAT ascore;
      PFLOAT SmIJ = Sm(J,I);
      if(!YleftClipped){
#if USE_MIC
	ascore = scoreSend1[J] + SmIJ;
#else
	ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-Lijx+I)) + SmIJ;
#endif
	ascore = max(ascore,C);
	if(ENDFIX>=2){
	  for(int K= Lijx; K < J; K++){
	    PFLOAT bound = Sbnd(X[J]-X[K],Yf[I],(PFLOAT)(J+I-K)) + SmIJ;
	    if(bound <= ascore)
	      break;
	    ascore = bound;
	  }
	} 
      } else
	ascore = scoreI[J];
      if(DEBUG>=2) assert(ascore >= C);
      scoreI[J] = ascore;
      if(ENDOUTLIER_FIX && ascore > C)
	FLAT_A(A.G, I, J) = -1;
    }
  }

  I = max(IMIN, LijxImax + 1);

  if(DEBUG && I <= IMAX && (I <= DELTA || min(DELTA,JMAX[I]) >= JMIN[I]) && !(A.Lijy[I][(I<=DELTA) ? JMAX[I] : min(DELTA,JMAX[I])] > 0)){
    #pragma omp critical
    {
      int Jmax = (I<=DELTA) ? JMAX[I] : min(DELTA,JMAX[I]);
      printf("pairalignXY:Yid=%d(%lld),Xid=%d(%lld),or=%d,sc=%d:I=%d(IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d),N=%d,M=%d,A.Lijy[I][Jmax=%d] = %d\n",
	     Yid,YYmap[Yid]->id, Xid, XXmap[Xid]->id,orientation,scaleID,I,IMIN,IMAX,JMIN[I],JMAX[I],N,M,Jmax, A.Lijy[I][Jmax]);
      printf("\t LijxImax= %d, LijyJmax= %d, DELTA=%d\n",LijxImax,LijyJmax,DELTA);
      fflush(stdout);
      assert(A.Lijy[I][J] > 0);
    }
  }

  for(; I <= IMAX; I++){
    int Jmax = (I <= DELTA) ? JMAX[I] : min(DELTA,JMAX[I]);
    if(RepeatShift > 0.0){
      Jmax = min(I-1,Jmax);
      while(Jmax >= 1 && Y[I]-X[Jmax] < RepeatShift)
	Jmax--;
    }

    int *LijyI = &FLAT_A(A.Lijy, I, 0);
    PFLOAT *scoreI = &FLAT_A(A.score, I, 0);

#if USE_MIC
    int score_size = Jmax;
    __mmask16 mask;
    __m512 vIdx=_mm512_add_ps(_mm512_set1_ps(I+1), v512_ascending0123);

    //prefetch_array(&Xf[1], score_size*sizeof(*Xf), _MM_HINT_T0);
    __m512 v_YfI=_mm512_set1_ps(Yf[I]), v1=_mm512_set1_ps(1.0f);
    __m512i v_limit=_mm512_sub_epi32(_mm512_set1_epi32(Jmax), v512i_ascending0123);
    
    /* make first pass of the loops to speedup computation */
    for(int J= JMIN[I]; J <= Jmax; J += 16){
      //mask= J<Jmax-7 ? 0xffff : (0xFFFF>>(J+7-Jmax));
      mask= _mm512_cmp_epi32_mask(_mm512_set1_epi32(J), v_limit, _MM_CMPINT_LE);
      __m512i LijyIJ = _mm512_maskloadu_epi32(&LijyI[J], mask);
      __m512 LijyIJd = _mm512_cvtfxpnt_round_adjustepi32_ps(LijyIJ, _MM_FROUND_CUR_DIRECTION, _MM_EXPADJ_NONE);
      __m512 v_score = Send512(_mm512_gmin_ps(_mm512_maskloadu_ps(&Xf[J], mask),
					       v_YfI
					       ),
			       _mm512_sub_ps(_mm512_add_ps(vIdx, _mm512_set1_ps(J)),
					     _mm512_gmax_ps(v1, LijyIJd)
					    )
			      );
      if(JMIN[I] == 1)
	_mm512_mask_store_ps(&scoreSend1[J], mask, v_score);
      else
	_mm512_maskstoreu_ps(&scoreSend1[J], mask, v_score);
    }
#endif
      
    for(int J= JMIN[I]; J <= Jmax; J++){
      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
      int Lijy = LijyI[J];
      if(DEBUG>=2) assert(Lijy > 0 && FLAT_A(A.Lijx,I,J) == 0);

      /* left end of X fully overlapped by Y */
      PFLOAT ascore;
      if(!XleftClipped){
	PFLOAT SmIJ = Sm(J,I);
#if USE_MIC
	ascore = scoreSend1[J] + SmIJ;
#else
	ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I+1-Lijy)) + SmIJ;
#endif
	ascore = max(ascore,C);
	if(ENDFIX>=2){
	  for(int K= Lijy; K < I; K++) {
	    PFLOAT bound = Sbnd(Xf[J],Y[I]-Y[K],(PFLOAT)(J+I-K)) + SmIJ;
	    if(bound <= ascore)
	      break;
	    ascore = bound;
	  }
	}
      } else
	ascore = scoreI[J];
      if(DEBUG>=2) assert(ascore >= C);
      scoreI[J] = ascore;
      if(ENDOUTLIER_FIX && ascore > C)
	FLAT_A(A.G, I, J) = -1;
    }
  } // for(; I <= IMAX; I++)

  /* initialize A[I=1..N][J= 1 to M].Rijy as rightmost site in Y that overlaps X when Y[I] is aligned with X[J] */
  /* If ENDINIT : only cases where I >= N+1-DELTA OR J >= M+1-DELTA are initialized */

  /* first estimate LijyJmin, the largest value J = M .. 1 so that X[M+1]-X[J] + 1e-8 >= Y[N+1] - Y[imin(J)] (0 if no such value J): use binary search to find this value between M and 1 */
  /* NOTE : imin(J) == (J >= M+1-DELTA) ? 1 : max(1, N+1-DELTA) */
  /* NOTE : For all J <= LijyJmin, Rijy == N+1. Conversly for all J > LijyJmin, Rijy <= N for I == imin(J) (but possibly NOT for any I in IMIN..IMAX) */
  int LijyJmin = 0;
  if(X[M+1] - X[1] + 1e-8 >= Y[N+1] - Y[!ENDINIT ? 1 : (1 >= M+1-DELTA) ? 1 : max(1,N+1-DELTA)]){/* LijyJmin >= 1 */
    if(X[M+1] - X[M] + 1e-8 >= Y[N+1] - Y[!ENDINIT ? 1 : (M >= M+1-DELTA) ? 1 : max(1,N+1-DELTA)])/* LijyJmin >= M */
      LijyJmin = M;
    else {
      int lb = 1, ub = M;
      while(lb + 1 < ub){
	int J = (lb + ub)/2;
	int imin = !ENDINIT ? 1 : (J >= M+1-DELTA) ? 1 : max(1,N+1-DELTA);
	if(X[M+1] - X[J] + 1e-8 >= Y[N+1] - Y[imin])
	  lb = J;
	else
	  ub = J;
      }
      if(DEBUG) assert(lb + 1 == ub);
      if(DEBUG>=2) assert(X[M+1] - X[lb] + 1e-8 >= Y[N+1] - Y[!ENDINIT ? 1 : (lb >= M+1-DELTA) ? 1 : max(1,N+1-DELTA)]);
      if(DEBUG>=2) assert(X[M+1] - X[ub] + 1e-8 < Y[N+1] - Y[!ENDINIT ? 1 : (ub >= M+1-DELTA) ? 1 : max(1,N+1-DELTA)]);
      LijyJmin = lb;
    }
  }

  // NOTE LijyJmin is typically close to M-DELTA-1
  for(int J = M; J > LijyJmin; J--){
    int imin = !ENDINIT ? Imin[J] : (J >= M+1-DELTA) ? Imin[J] : max(Imin[J],N+1-DELTA);
    
    int Rijy = N+1;
    int *RijyJ = &FLAT_A(A.Rijy, 0, J);
    long long stride = A.M_stride;

    for(int I= Imax[J]; I >= imin; I--){
      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);

      /* update Rijy by decrementing it until Y[Rijy] overlaps X */
      FLOAT delta = Y[I] + X[M+1] - X[J] + 1e-8;
      while(Y[Rijy] > delta)
	Rijy--;
      if(DEBUG>=2 && !(Rijy >= I)){
	printf("pairalignXY:Yid=%d(%lld),Xid=%d(%lld),or=%d,sc=%d:J=%d,I=%d,N=%d,M=%d:Rijy=%d,Y[I]= %0.4f,X[M+1]=%0.4f,X[J]= %0.4f, Y[Rijy]= %0.4f, Y[Rijy+1]= %0.4f, delta= %0.4f\n",
	       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,J,I,N,M,Rijy,Y[I],X[M+1],X[J],Y[Rijy],Y[Rijy+1],delta);
	fflush(stdout);
	assert(Rijy >= I);
      }
      RijyJ[I*stride] = Rijy;
    }
      
    if(DEBUG>=2 && !(phash && hashdelta)) assert(Rijy <= N);
  }

  int Rijy = N+1;
  DELTA_M = max(1, M+1-DELTA);

  if(DELTA_M <= LijyJmin){
    for(int I= IMIN; I <= IMAX; I++){
      int Jmax = min(LijyJmin, JMAX[I]);
      for(int J = max(DELTA_M, JMIN[I]); J <= Jmax; J++)
	FLAT_A(A.Rijy, I, J) = Rijy;
    }
  }

  J = min(LijyJmin, DELTA_M - 1);
  if(J >= 1){
    int imin = !ENDINIT ? IMIN : max(IMIN,N+1-DELTA);
    for(int I= IMAX; I >= imin; I--){
      int jmin = JMIN[I];
      int jmax = min(J, JMAX[I]);
      if(jmin <= jmax)
	memset_int(&FLAT_A(A.Rijy, I, jmin), Rijy, jmax-jmin+1);
    }
  }

  /* initialize A[I=1..N][J= 1 to M].Rijx as rightmost site in X that overlaps Y when Y[I] is aligned with X[J] */
  /* If ENDINIT : only cases where I >= N+1-DELTA OR J >= M+1-DELTA are initialized */

  /* first estimate LijxImin, the largest value I = N .. 1 so that X[M+1] - X[Jmin(I)] <= Y[N+1] - Y[I] + 1e-8 (or 0 if no such I exists) : use binary search to find this value between 1 and N */
  /* NOTE : Jmin(I) == (I >= N+1-DELTA) ? 1 : max(1,M+1-DELTA) */
  /* NOTE : for all I <= LijxImin, Rijx == M+1. Conversely for all I > LijxImin, Rijx <= M for J == Jmin(I) (but possible not for any J in JMIN[I]..JMAX[I]) */
  int LijxImin = 0;  
  if(X[M+1] - X[!ENDINIT ? 1 : (N <= DELTA) ? 1 : max(1, M+1-DELTA)] <= Y[N+1] - Y[1] + 1e-8){/* LijxImin >= 1 */
    if(X[M+1] - X[!ENDINIT ? 1 : (1 <= DELTA) ? 1 : max(1, M+1-DELTA)] <= Y[N+1] - Y[N] + 1e-8) /* LijxImin >= N */
      LijxImin = N;
    else {
      int lb = 1, ub = N;
      while(lb + 1 < ub){
        int I = (lb + ub)/2;
	int Jmin = !ENDINIT ? 1 : (I >= N+1-DELTA) ? 1 : max(1, M+1-DELTA);
	if(X[M+1] - X[Jmin] <= Y[N+1] - Y[I] + 1e-8)
	  lb = I;
	else
	  ub = I;
      }
      if(DEBUG) assert(lb + 1 == ub);
      if(DEBUG>=2) assert(X[M+1] - X[!ENDINIT ? 1 : (lb >= N+1-DELTA) ? 1 : max(1,M+1-DELTA)] <= Y[N+1] - Y[lb] + 1e-8);
      if(DEBUG>=2) assert(X[M+1] - X[!ENDINIT ? 1 : (ub >= N+1-DELTA) ? 1 : max(1,M+1-DELTA)] > Y[N+1] - Y[ub] + 1e-8);
      LijxImin = lb;
    }
  }

  // NOTE LijxImax is typically around N-DELTA
  for(int I = IMAX; I > LijxImin && I >= IMIN; I--){
    int Jmin = !ENDINIT ? JMIN[I] : (I >= N+1-DELTA) ? JMIN[I] : max(JMIN[I],M+1-DELTA);

    int Rijx = M+1;
    int *RijxI = &FLAT_A(A.Rijx, I, 0);
    for(int J= JMAX[I]; J >= Jmin; J--){
      /* update Rijx by decrementing it until X[Rijx] overlaps Y */
      FLOAT delta = X[J] + Y[N+1] - Y[I] + 1e-8;
      while(X[Rijx] > delta)
	Rijx--;
      if(DEBUG>=1+RELEASE)	assert(Rijx >= J);
      RijxI[J] = Rijx;
    }

    if(DEBUG>=2 && !(phash && hashdelta)) assert(Rijx <= M);
  }


  int Rijx = M+1;
  DELTA_N = max(IMIN, N+1-DELTA);

  if(DELTA_N <= LijxImin){
    int imax = min(LijxImin, IMAX);
    for(int I = DELTA_N; I <= imax; I++)
      if(JMIN[I] <= JMAX[I])
	memset_int(&FLAT_A(A.Rijx, I, JMIN[I]), Rijx, JMAX[I]-JMIN[I]+1);
  }

  I = min(LijxImin, DELTA_N - 1);
  if(I >= IMIN){
    for(int i = IMIN; i <= I; i++){
      int Jmin = !ENDINIT ? JMIN[i] : max(JMIN[i],M+1-DELTA);
      if(Jmin <= JMAX[i])
	memset_int(&FLAT_A(A.Rijx, i, Jmin), Rijx, JMAX[i]-Jmin+1);
    }
  }

  if(DEBUG>=2){
    for(int I = IMIN; I <= IMAX; I++){
      int Jmin = !ENDINIT ? JMIN[I] : max(JMIN[I],M+1-DELTA);
      for(int J = Jmin; J <= JMAX[I]; J++){
	int Rijx = FLAT_A(A.Rijx, I, J);
	int Rijy = FLAT_A(A.Rijy, I, J);
	if(!(Rijx > M || Rijy > N) || Rijx < J || Rijy < I){
	  printf("Yid=%d(%lld),Xid=%d(%lld):I=%d,J=%d:Rijx=%d,Rijy=%d,M=%d,N=%d,DELTA=%d,A.M_stride=%lld\n",
		 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,I,J,Rijx,Rijy,M,N,DELTA,A.M_stride);
	  printf("  Y[I]=%0.8f,Y[N+1]=%0.8f,X[J]=%0.8f,X[M+1]=%0.8f\n",Y[I],Y[N+1],X[J],X[M+1]);
	  fflush(stdout);
	  assert(Rijx > M || Rijy > N);
	  assert(Rijx >= J);
	  assert(Rijy >= I);
	}
      }
    }
  }

  if(VERB>=3 && Yid==69 && Xid==70 && orientation==0 && N >= 19 && M >= 48103){
    #pragma omp critical
    {
      printf("tid=%d:Yid=%d,Xid=%d,or=%d,sc=%d:N=%d,M=%d:FLAT_A(A.Rijy,19,48103)= %d (%p)\n",
	     tid,Yid,Xid,orientation,scaleID,N,M,FLAT_A(A.Rijy,19,48103), &FLAT_A(A.Rijy,19,48103));
      fflush(stdout);
    }
  }

  if(MEMCHECK>=3 && !MEMRESERVE && RA_MIN_TIME && MaxMemSiz > 0 && RA_MIN_MEM < 8 && retries < MEMCHECK-1 && 
     ((wt = wtime()) > Mem->lastGetStatm + VMRSS_CHECK && (wt > Mem->lastGetStatm + VMRSS_CHECK2 || Mem->VmRSS + Mem->VmSwap >= MaxMemSiz2 * RA_MAX_MEM))){/* recheck VmRSS */
    //    long long maxVmRSS2 = MaxMemSiz2;// maximum VmRSS+VmSwap allowed for large memory alignments (must fit within 90% of memory)
    //    long long maxVmRSS1 = MaxMemSiz;// maxmimum VmRSS+VmSwap allowed for jobs expected to use less than 10% of the per thread amount (so total will fit within 100% of memory)
    //    bool SmallJob = nMemSiz < (MaxMemSiz / (NumThreads * 10.0));
    //    long long maxVmRSS = SmallJob ? maxVmRSS1 : maxVmRSS2;// maximum VmRSS+VmSwap allowed for current pairwise alignment

    if(THREAD_DEBUG) Mem->getstatm_calls++;
    Mem->GetStatm(wt);

    if(Mem->VmRSS + Mem->VmSwap > MaxMemSiz /* WAS108 maxVmRSS */){
      if(VERB>=2){
        #pragma omp critical(PausedThreads)
	{
	  if(PausedThreads > 0 || (VMRSS_PEAKCHECK && Mem->MemUsed > VMRSS_PEAKCHECK)){
	    printf("tid=%d: MaxMem= %0.1f: GetStatm: VmRSS= %0.4f, VmSwap= %0.2f, MemUsed= %0.4f(nMem=%0.4f) (Yid=%lld,Xid=%lld,N=%d,M=%d),retries=%d(C), Paused Threads=%d:wt=%0.6f,%0.6f\n",
		   tid, MaxMemSiz*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, Mem->MemUsed*1e-9, nMemSiz*1e-9,YYmap[Yid]->id,XXmap[Xid]->id,N,M,retries,PausedThreads,Mem->lastGetMem,Mem->lastGetStatm);
	    fflush(stdout);
          }
	}
      }
      retries++;

      if(MEMRESERVE){
        if(THREAD_DEBUG) Mem->atomic_calls++;

        #pragma omp atomic update
        VmRSSreserved -= Mem->MemReserved;
      
        Mem->MemReserved = 0;
      }

      goto LwaitForMemory;
    }
  }

  if((DEFER_BP && !myoutlierExtend) || RepeatRec){ /* save score as Uscore */
    for(int I= IMIN; I <= IMAX; I++)
      memcpy(&FLAT_A(A.Uscore, I, JMIN[I]), &FLAT_A(A.score, I, JMIN[I]), (JMAX[I]-JMIN[I]+1)*sizeof(PFLOAT));

    if(EVERB && IMIN <= I_TRACE && I_TRACE <= IMAX && JMIN[I_TRACE] <= J_TRACE && J_TRACE <= JMAX[I_TRACE]){
      printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d, A.Uscore(I,J)= %0.6f, ChimScore= %0.6f\n",Yid,Xid,orientation,scaleID,I_TRACE,J_TRACE,FLAT_A(A.Uscore,I_TRACE,J_TRACE),ChimScore);
      fflush(stdout);
    }
  }

  if(VERB>=3 && Yid==69 && Xid==70 && orientation==0 && N >= 19 && M >= 48103){
    #pragma omp critical
    {
      printf("B:tid=%d:Yid=%d,Xid=%d,or=%d,sc=%d:N=%d,M=%d:FLAT_A(A.Rijy,19,48103)= %d (%p)\n",
	  tid,Yid,Xid,orientation,scaleID,N,M,FLAT_A(A.Rijy,19,48103), &FLAT_A(A.Rijy,19,48103));
      fflush(stdout);
    }
  }

  if(DEFER_BP && !(phash && hashdelta) && !myoutlierExtend && !myoutlierBC && DELTA==4){/* vectorized code will be used */
    if(VERB && (TSAN || !VecMessage)){
      #pragma omp critical
      {
        if(!VecMessage){
	  printf("pairalignXY: Using Vectorized Code: hashdelta=%0.1f,myoutlierExtend=%d, myoutlierBC=%d, DELTA=%d,tid=%d: wall time= %0.6f secs\n", 
	                      phash ? hashdelta : 0.0,myoutlierExtend, myoutlierBC, DELTA, tid, wtime());
          fflush(stdout);
	  VecMessage = 1;
        }
      }
    }
  } else { /* non-vectorized code will be used */
    if(VERB && (TSAN || !NoVecMessage)){
      #pragma omp critical
      {
	if(!NoVecMessage){
	  printf("pairalignXY: Using Non-Vectorized Code: hashdelta= %0.1f, myoutlierExtend=%d .. %d, myoutlierBC=%d, DELTA=%d,tid=%d: wall time= %0.6f secs\n", 
	            phash ? hashdelta : 0.0, myoutlierExtend, outlierExtendLimit, myoutlierBC, DELTA, tid, wtime());
	  fflush(stdout);
	  
	  NoVecMessage = 1;
	}
      }
    }
  }

  //  long long myVmRSSreserved;
  long long origMemReserved = Mem->MemReserved;

  if(MEMRESERVE){
    if(THREAD_DEBUG) Mem->atomic_calls++;

    if(VERB>=3 && MEMCHECK <= 1 && Mem->MemReserved >= 10000000LL){
      #pragma omp critical(PausedThreads)
      {
	if(PausedThreads > 0){
          printf("tid=%d: MaxMem= %0.2f: VmRSS= %0.4f, VmSwap= %0.2f, Reserved= %0.4f -> %0.4f, Used= %0.4f(nMem=%0.4f) (Yid=%lld,Xid=%lld,N=%d,M=%d),retries=%d(D), PausedThreads=%d:wt=%0.6f,%0.6f,%0.6f\n",
	    tid, MaxMemSiz*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, VmRSSreserved*1e-9, (VmRSSreserved - Mem->MemReserved)*1e-9, Mem->MemUsed*1e-9, nMemSiz*1e-9,
	       YYmap[Yid]->id,XXmap[Xid]->id,N,M,retries,PausedThreads,Mem->lastGetMem,Mem->lastGetStatm,wtime());
	  fflush(stdout);
        }
      }
    }

    #pragma omp atomic update
    VmRSSreserved -= Mem->MemReserved;
      
    Mem->MemReserved = 0;
  }

  if(MEMCHECK>=2 && RA_MIN_TIME && MaxMemSiz > 0 && RA_MIN_MEM < 8 && retries < MEMCHECK-1 && 
     (((wt = wtime()) > Mem->lastGetStatm + VMRSS_CHECK && (wt > Mem->lastGetStatm + VMRSS_CHECK2 || Mem->VmRSS + Mem->VmSwap >= MaxMemSiz2 * RA_MAX_MEM)) ||
	 (VMRSS_PEAKCHECK && Mem->MemUsed > VMRSS_PEAKCHECK))){/* recheck VmRSS */

    if(THREAD_DEBUG) Mem->getstatm_calls++;
    Mem->GetStatm(wt);

    if(Mem->VmRSS + Mem->VmSwap /* WAS29 + (MEMRESERVE ? myVmRSSreserved : 0) */ > MaxMemSiz){
      if(VERB>=2 || retries > 0){
        #pragma omp critical(PausedThreads)
	{
	  if(PausedThreads > 0 || (VMRSS_PEAKCHECK && Mem->MemUsed > VMRSS_PEAKCHECK)){
	    printf("tid=%d: MaxMem= %0.2f: GetStatm: VmRSS= %0.4f, VmSwap= %0.2f, Reserved= %0.4f -> %0.4f, Used= %0.4f(nMem=%0.4f) (Yid=%lld,Xid=%lld,N=%d,M=%d),retries=%d(D), Paused Threads=%d:wt=%0.6f,%0.6f,%0.6f\n",
		   tid, MaxMemSiz*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, (VmRSSreserved + origMemReserved)*1e-9, VmRSSreserved*1e-9, Mem->MemUsed*1e-9, nMemSiz*1e-9,
		   YYmap[Yid]->id,XXmap[Xid]->id,N,M,retries,PausedThreads,Mem->lastGetMem,Mem->lastGetStatm,wtime());
	    if(VERB/* HERE >=3 */){
	      long long VmSize,VmRSS,VmSwap,VmHWM;// local copies, so global shared values are not modified
	      getmem(VmSize,VmRSS,VmSwap,&VmHWM);
	      printf("\t Actual VmRSS= %0.4f(HWM= %0.4f), VmSwap= %0.4f, VmSize= %0.4f Gb:wall=%0.6f\n",VmRSS*1e-9,VmHWM*1e-9,VmSwap*1e-9,VmSize*1e-9,wtime());
	    }
	    fflush(stdout);
          }
	}
      }
      retries++;

      goto LwaitForMemory;
    }

    if(VERB>=3){
      #pragma omp critical(PausedThreads)
      {
      if(PausedThreads > 0 /* || nMemSiz > 1000000000LL */){
	  printf("tid=%d: MaxMem= %0.2f: GetStatm: VmRSS= %0.4f, VmSwap= %0.2f, Reserved=%0.4f -> %0.4f, Used= %0.4f(nMem=%0.4f) (Yid=%lld,Xid=%lld,N=%d,M=%d), Paused Threads=%d:wt=%0.6f,%0.6f\n",
	    tid, MaxMemSiz*1e-9, Mem->VmRSS*1e-9, Mem->VmSwap*1e-9, (VmRSSreserved + origMemReserved)*1e-9, VmRSSreserved*1e-9,Mem->MemUsed*1e-9, nMemSiz*1e-9,
	    YYmap[Yid]->id,XXmap[Xid]->id,N,M,PausedThreads,Mem->lastGetMem,Mem->lastGetStatm);
	  fflush(stdout);
	}
      }
    }
  }

  /* Dynamic programming recurance */
  for(int I = max(2,IMIN); I <= IMAX; I++){
	  
    int Jmin = max(2,JMIN[I]), Jmax = JMAX[I];
    if(RepeatShift > 0.0 && !orientation){/* decrease Jmax to below I until Y[I] - X[Jmax] >= RepeatShift */
      Jmax = I-1;
      while(Jmax >= Jmin && Y[I]-X[Jmax] < RepeatShift)
	Jmax--;
    }
    
    if(DEFER_BP && !(phash && hashdelta) && !myoutlierExtend && !myoutlierBC && DELTA==4 && I > 4 && A.M_stride <= (long long)MASK_LL(31)){/* critical case, customized to help compiler to unroll the G and H loops */
      // NOTE : Currently only (USE_SSE && USE_AVX && USE_PFLOAT) OR (USE_MIC && USE_PFLOAT) use intrinsics
      // Normally USE_SSE is defined as 1 IFF  __SSE__ and __SSE2__ are defined, but can also be forced to 0 in Makefile for debugging.
      // NOTE : USE_AVX==0 (only used in pairalign.cpp) means NEVER use intrinsics for AVX. USE_AVX==1 means use intrinsincs for AVX IF available (AND USE_SSE==1) 
      // Normally USE_AVX is defined as 1 IFF __AVX__ is defined, but can also be forced to 0 in Makefile for debugging.

      int jmax = min(4,Jmax);
      int Gmin = I-4;// NOTE : Gmin >= IMIN

 #if !(USE_SSE && USE_AVX && USE_PFLOAT) && !(USE_MIC && USE_PFLOAT)

      /* Original path */    

      int J = Jmin;
      for(; J <= jmax; J++){
	PFLOAT score = FLAT_A(A.score, I, J);
	
	for(int G= I; --G >= Gmin;){
	  PFLOAT deltaY = delY[I][G-I+4];// Y[I] - Y[G];
	  if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	  for(int H= J; --H > 0;){
	    if(DEBUG>=2) assert(fabs(X[J]-X[H] - delX[J][H-J+4]) < 0.001);
	    PFLOAT newscore = FLAT_A(A.score, G, H) + Sint(delX[J][H-J+4],deltaY,J-H,I-G);
	    score = PMAX(score,newscore);
	  }
	}
	FLAT_A(A.score, I, J) = score;
      }

      for(; J <= Jmax; J++){
	PFLOAT score = FLAT_A(A.score, I, J);

	for(int G= I; --G >= Gmin;){
	  PFLOAT deltaY = delY[I][G-I+4];
	  if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	  for(int H= J; --H >= J - 4;){
	    if(DEBUG>=2) assert(fabs(X[J]-X[H] - delX[J][H-J+4]) < 0.001);
	    PFLOAT newscore = FLAT_A(A.score, G, H) + Sint(delX[J][H-J+4],deltaY,J-H,I-G);
	    score = PMAX(score,newscore);
	  }
	}	
	FLAT_A(A.score, I, J) = score;
      }
 #else  // (USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)

      int J;

   // NOTE : This block only handles the first 4 iterations of J, the others are handled later

   #if (USE_SSE && USE_AVX && USE_PFLOAT) 
      if(DEBUG) assert(USE_SSE && USE_AVX && USE_PFLOAT);

      float fm[4]={4, 3, 2, 1}, tmp[4];
      __m128 vfm=_mm_load_ps(fm);

      if(DEBUG>=2) assert(jmax == min(4,Jmax));

      for(int J = Jmin; J <= jmax; J++){

	__m128 vns = _mm_set_ps1(FLAT_A(A.score, I, J));
	for(int H= J; --H > 0;){
	  PFLOAT deltaX = delX[J-H][J];// X[J] - X[H];
	 
//	  for(int G= I; --G >= I - 4;){
//	    PFLOAT deltaY = delY[I][G-I+4];// Y[I] - Y[G];
	  __m128 vy = _mm_loadu_ps(&delY[I][0]);// X[J] - X[J-4] ... X[J] - X[J-1] : only the last J-1 values of the vector are valid 
	  tmp[0] = A.score[I-4][H];
	  tmp[1] = A.score[I-3][H];
	  tmp[2] = A.score[I-2][H];
	  tmp[3] = A.score[I-1][H];
	  vns = _mm_max_ps(vns, _mm_add_ps(_mm_load_ps(tmp), Sint128(_mm_set_ps1(deltaX), vy, _mm_set_ps1((float)(J-H)), vfm)));
//	  }
	}
	vns=_mm_max_ps(vns, _mm_movehl_ps(vns, vns));

#ifndef __SSE3__
	vns=_mm_max_ps(vns, _mm_shuffle_ps(vns, vns, _MM_SHUFFLE(1,1,1,1)));
#else
	vns=_mm_max_ss(vns, _mm_movehdup_ps(vns));
#endif

	_mm_store_ss(&(FLAT_A(A.score, I, J)), vns);
      }

      J = jmax+1;

  #else // (USE_MIC && USE_PFLOAT)
      if(DEBUG) assert(USE_MIC && USE_PFLOAT);

      //Sint512_prefetch();
      int Hmin = 1;
      PFLOAT *delYI = delY[I];

      __m512i Hvi, Gvi, Gv, Hv, IGv, JHv, Iv, Jv, AvIdx;
      __m512 delXv, delYv, Av, score_v;
      __mmask16 mask;
      Hvi = _mm512_set_epi32(0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3);
      Gvi = _mm512_set_epi32(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3);

      Gv = _mm512_add_epi32(_mm512_set1_epi32(Gmin), Gvi);
      Hv = _mm512_add_epi32(_mm512_set1_epi32(Hmin), Hvi);
      AvIdx = _mm512_add_epi32(_mm512_mullo_epi32(_mm512_sub_epi32(Gv, _mm512_set1_epi32(1)), _mm512_set1_epi32(A.M_stride)), Hv);
//	Av = _mm512_mask_i32gather_ps(v512_PZERO, mask, AvIdx, A.score[1], 4);
      Av= _mm512_i32gather_ps(AvIdx, A.score[1], 4);
	
      J = Jmin;
      {
	PFLOAT *AscoreIJmin = &FLAT_A(A.score, I, Jmin);
	PFLOAT *delXJmin = &(delX[1][Jmin]);

	if(DEBUG>=2) assert(jmax == min(4,Jmax));

	IGv = _mm512_sub_epi32(_mm512_set1_epi32(I), Gv);
	delYv = _mm512_i32gather_ps(_mm512_sub_epi32(_mm512_set1_epi32(DELTA), IGv), delYI, 4); 
	
#define BLOCKJ(Ja) if(Ja <= (jmax - Jmin)){ \
	 int J = Ja + Jmin;					\
	 __m512 score_orig_v = _mm512_set1_ps(AscoreIJmin[Ja]);	\
	 JHv = _mm512_sub_epi32(_mm512_set1_epi32(J), Hv);		\
	 mask = _mm512_cmp_epi32_mask(JHv, _mm512_set1_epi32(0), _MM_CMPINT_GT); \
	 delXv = _mm512_mask_i32gather_ps(v512_PZERO, mask, _mm512_mullo_epi32(_mm512_set1_epi32(COMPUTE_STRIDE4(M-1)), _mm512_sub_epi32(JHv, _mm512_set1_epi32(1))), &delXJmin[Ja], 4); \
	 score_v = _mm512_mask_add_ps(score_orig_v, mask, Av, Sint512b(delXv, delYv, _mm512_cvtfxpnt_round_adjustepu32_ps(_mm512_add_epi32(JHv, IGv), _MM_FROUND_CUR_DIRECTION, _MM_EXPADJ_NONE))); \
	 AscoreIJmin[Ja] = _mm512_reduce_gmax_ps(_mm512_gmax_ps(score_orig_v, score_v)); \
	}
      
	BLOCKJ(3);
	BLOCKJ(2);
	BLOCKJ(1);
	BLOCKJ(0);
	J = jmax+1;
      }
      
  #endif // else (USE_MIC && USE_PFLOAT)
 
  //  continue with else (USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT) : Not all J iterations up to Jmax-1 have been completed


   #if (USE_SSE && USE_AVX && USE_PFLOAT) 
      if(DEBUG) assert(!USE_MIC);
 
      if(DEBUG>=2) assert(J == jmax+1);

      int iterations = (Jmax - jmax) / 8;

      for(int J = jmax+1; J <= Jmax-7; J += 8){

	__m256 v_score = _mm256_loadu_ps(&FLAT_A(A.score, I, J));
	//	__m256 v_Xj = _mm256_loadu_ps(&X[J]);
	for(int G= I; --G >= Gmin;){
	  PFLOAT deltaY = delY[I][G-I+4];// Y[I] - Y[G];
	  if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	  __m256 v_deltaY = _mm256_broadcast_ss(&deltaY);
	  float *AscoreGJ = &FLAT_A(A.score, G, J);
          #define BLOCK(H) { \
 	    __m256 v_prevscore = _mm256_loadu_ps(&AscoreGJ[-H]); \
	    __m256 v_newscore = _mm256_add_ps(v_prevscore, Sint256(_mm256_loadu_ps(&delX[H][J]), v_deltaY, H, I-G));	\
	    v_score = _mm256_max_ps(v_score, v_newscore); \
	    }
	    
	    BLOCK(1);
	    BLOCK(2);
	    BLOCK(3);
	    BLOCK(4);
          #undef BLOCK
	}	
	_mm256_storeu_ps(&FLAT_A(A.score, I, J), v_score);
      }

      if(DEBUG>=2){
        for(J = jmax+1; J <= Jmax-7; J += 8);

        if(DEBUG && !(J == jmax+1 + iterations * 8)){
	  #pragma omp critical
	  {
	    printf("jmax=%d,Jmax=%d,J=%d,iterations=%d\n",jmax,Jmax,J,iterations);
	    fflush(stdout);
	    assert(J == jmax+1 + iterations * 8);
	  }
        }
      } else
	J = jmax+1 + iterations * 8;

      if(J <= Jmax) {
	unsigned int mask[14] = {0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0, 0, 0, 0, 0, 0, 0};
	__m256i mask256 = _mm256_loadu_si256((__m256i*) &(mask[6+J-Jmax]));

	__m256 v_score = _mm256_maskload_ps(&FLAT_A(A.score, I, J), mask256);
	for(int G= I; --G >= Gmin;){
	  PFLOAT deltaY = delY[I][G-I+4];// Y[I] - Y[G];
	  if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	  __m256 v_deltaY = _mm256_broadcast_ss(&deltaY);
	  float *AscoreGJ = &FLAT_A(A.score, G, J);

          #define BLOCK(H) { \
 	    __m256 v_prevscore = _mm256_maskload_ps(&AscoreGJ[-H], mask256); \
	    __m256 v_newscore = _mm256_add_ps(v_prevscore, Sint256(_mm256_maskload_ps(&delX[H][J], mask256), v_deltaY, H, I-G));	\
	    v_score = _mm256_max_ps(v_score, v_newscore); \
	    }
	    
	    BLOCK(1);
	    BLOCK(2);
	    BLOCK(3);
	    BLOCK(4);
          #undef BLOCK
	}	
	_mm256_maskstore_ps(&FLAT_A(A.score, I, J), mask256, v_score);
	J = Jmax+1;
      }

  #endif // USE_SSE && USE_AVX along with USE_PFLOAT

  #if (USE_MIC && USE_PFLOAT)
      delYI = &delY[I][4];
      //Sint512_prefetch();
      for(/* J = jmax+1*/; J <= Jmax-15; J += 16){
	__m512 v_score = _mm512_loadu_ps(&FLAT_A(A.score, I, J));
#define DELX_BLOCK(H) __m512 delX##H = _mm512_loadu_ps(&FLAT_DELX(delX, H, J, M));
	
	DELX_BLOCK(1);
	DELX_BLOCK(2);
	DELX_BLOCK(3);
	DELX_BLOCK(4);
	
#undef DELX_BLOCK

#define BLOCKH(G, H) { \
	  __m512 v_prevscore = _mm512_loadu_ps(&AscoreGJ[-H]);		\
	  __m512 v_newscore = _mm512_add_ps(v_prevscore, Sint512b(delX##H, v_deltaY, _mm512_set1_ps((float)(H+G)))); \
	  v_score = _mm512_gmax_ps(v_score, v_newscore);		\
	}

#define BLOCKG(G) { \
	  PFLOAT deltaY = delYI[-G];		      \
	  __m512 v_deltaY = _mm512_set1_ps(deltaY);	\
	  float *AscoreGJ = &FLAT_A(A.score, I-G, J);	\
	  						\
	  BLOCKH(G, 1);					\
	  BLOCKH(G, 2);					\
	  BLOCKH(G, 3);					\
	  BLOCKH(G, 4);					\
	}	
	
	BLOCKG(1);
	BLOCKG(2);
	BLOCKG(3);
	BLOCKG(4);

#undef BLOCK
#undef BLOCKG

	_mm512_storeu_ps(&FLAT_A(A.score, I, J), v_score);
      } // for J <= Jmax-15
      
      if(J <= Jmax) {
	//Sint512_prefetch();
      __mmask16 mask16 = _mm512_int2mask((1<<(Jmax-J+1))-1); 
      __m512 v_score = _mm512_maskloadu_ps(&FLAT_A(A.score, I, J), mask16);
	
#define DELX_BLOCK(H) __m512 delX##H = _mm512_maskloadu_ps(&FLAT_DELX(delX, H, J, M), mask16);
	
	DELX_BLOCK(1);
	DELX_BLOCK(2);
	DELX_BLOCK(3);
	DELX_BLOCK(4);
	
#undef DELX_BLOCK

          #define BLOCK(G, H) { \
	  __m512 v_prevscore = _mm512_maskloadu_ps(&AscoreGJ[-H], mask16); \
	  __m512 v_newscore = _mm512_add_ps(v_prevscore, Sint512b(delX##H, v_deltaY, _mm512_set1_ps((float)(H+G)))); \
	  v_score = _mm512_gmax_ps(v_score, v_newscore);		\
	}

#define BLOCK2(G) { \
	  __m512 v_deltaY = _mm512_set1_ps(delYI[-G]);	\
	  PFLOAT *AscoreGJ = &FLAT_A(A.score, I-G, J);	\
							\
	  BLOCK(G, 1);					\
	  BLOCK(G, 2);					\
	  BLOCK(G, 3);					\
	  BLOCK(G, 4);					\
	}

	BLOCK2(1);
	BLOCK2(2);
	BLOCK2(3);
	BLOCK2(4);
#undef BLOCK
#undef BLOCK2
	_mm512_maskstoreu_ps(&FLAT_A(A.score, I, J), mask16, v_score);
	J = Jmax+1;
      } // if(J <= Jmax)

  #endif // USE_MIC

 #endif // (USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)

    } else if (USE_MIC && USE_PFLOAT && DEFER_BP && !(phash && hashdelta) && !myoutlierExtend && !myoutlierBC && DELTA==4 && A.M_stride <= (long long)MASK_LL(31)) {/* case with I <= 4 optimized with USE_MIC && USE_PFLOAT */
#if USE_MIC && USE_PFLOAT
      if(DEBUG>=2) assert(I <= 4);
      
      PFLOAT *delYI = delY[I];
      int H, G;
      PFLOAT deltaX[4], deltaY[4], score, tmp[16];
      int Gmin = max(I-DELTA,1);// max(I-DELTA,IMIN);

      //Sint512_prefetch();
      __m512i Hvi, Gvi, Gv, Hv, IGv, JHv, Iv, Jv, Gv_stride;
      __m512 delXv, delYv, Av, score_v;
      __mmask16 mask, mask_IGv;
      Hvi = _mm512_set_epi32(0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3);
      Gvi = _mm512_set_epi32(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3);
      Gv = _mm512_add_epi32(_mm512_set1_epi32(Gmin), Gvi);
      IGv = _mm512_sub_epi32(_mm512_set1_epi32(I), Gv);
      Gv_stride = _mm512_mullo_epi32(_mm512_sub_epi32(Gv, _mm512_set1_epi32(1)), _mm512_set1_epi32(A.M_stride));
      mask_IGv = _mm512_cmp_epi32_mask(IGv, _mm512_setzero_epi32(), _MM_CMPINT_GT);
      delYv = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), mask_IGv, _mm512_sub_epi32(_mm512_set1_epi32(DELTA), IGv), delYI, 4);
	
      for(int J = Jmin; J <= Jmax; J++){
	score = FLAT_A(A.score, I, J);
	int Hmin = max(J-DELTA,1);// max(J-DELTA,JMIN[I]);

	Hv = _mm512_add_epi32(_mm512_set1_epi32(Hmin), Hvi);
	JHv = _mm512_sub_epi32(_mm512_set1_epi32(J), Hv);
	mask = _mm512_kand(mask_IGv, _mm512_cmp_epi32_mask(JHv, _mm512_setzero_epi32(), _MM_CMPINT_GT));
	delXv = _mm512_mask_i32gather_ps(
					 _mm512_setzero_ps(), 
					 mask, 
					 _mm512_mullo_epi32(
							    _mm512_set1_epi32(COMPUTE_STRIDE4(M-1)), 
							    _mm512_sub_epi32(JHv, 
									     _mm512_set1_epi32(1))),
					 &(delX[1][J]), 4);
	Av = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), mask, _mm512_add_epi32(Gv_stride, Hv), A.score[1], 4);
	score_v = _mm512_mask_add_ps(_mm512_set1_ps(score), mask, Av, 
				     Sint512b(delXv, delYv, 
					      _mm512_cvtfxpnt_round_adjustepu32_ps(
										   _mm512_add_epi32(JHv, IGv), _MM_FROUND_CUR_DIRECTION, _MM_EXPADJ_NONE)));
	PFLOAT newscore = _mm512_reduce_gmax_ps(score_v);
	score = max(newscore,score);

	FLAT_A(A.score, I, J) = score;
      }// for(J=Jmin; J <= Jmax; J++)
#endif // USE_MIC && USE_PFLOAT
    } else {/* general case without (DEFER_BP && !myoutlierExtend) etc */
      PFLOAT OutPen = OutlierBias + OutlierPenalty;
      PFLOAT OutPen2 = OutPen + 1e-9;
      PFLOAT *delYI = delY[I];

      if(DEBUG>=2 && !(RepeatShift > 0.0 && !orientation)) assert(Jmin == max(2,JMIN[I]) && Jmax == JMAX[I]);

      for(int J = Jmin; J <= Jmax; J++){
	PFLOAT score = FLAT_A(A.score, I, J);
	if(!DEFER_BP || myoutlierExtend){
	  PFLOAT deltaX, deltaY, newscore;
	  int bestG = localtype, bestH = localtype;

	  int Gmin = max(I-DELTA,IMIN), Hmin = max(J-DELTA,JMIN[Gmin]);
	  
	  for(int H= J; --H >= Hmin;){
	    if(DEBUG>=2 && !(1 <= H && H < J && J-H <= DELTA && J <= JMAX[I])){
	      #pragma omp critical
	      {
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:I=%d,J=%d,H=%d,M=%d:tid=%d\n",
	           Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,H,M, tid);
		fflush(stdout);
		assert(1 <= H && H < J && J-H <= DELTA && J <= JMAX[I]);
	      }
            }

	    deltaX = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];// X[J] - X[H]

	    if(DEBUG>=2 && !(fabs(X[J] - X[H] - deltaX) < 0.001)){
	      #pragma omp critical
	      {
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:I=%d,J=%d,H=%d,M=%d:deltaX= %0.4f, X[J]= %0.4f, X[H]= %0.4f (X= %p, Y= %p)):tid=%d\n",
	           Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,H,M,deltaX,X[J],X[H], X, Y, tid);
		fflush(stdout);
		assert(fabs(X[J] - X[H] - deltaX) < 0.001);
	      }
	    }

	    for(int G= Gmin; G < I; G++){
	      if(H > JMAX[G])
		continue;
	      if(H < JMIN[G])
		break;
	      deltaY = delYI[G-I+DELTA];// Y[I] - Y[G] 
	      if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	      newscore = FLAT_A(A.score, G, H) + Sint(deltaX,deltaY,J-H,I-G,myoutlierBC);

	      if(EVERB && I==I_TRACE && J==J_TRACE /* && ((G==G_TRACE && H==H_TRACE) || newscore > score) */){
	        #pragma omp critical
		{
		  printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d:deltaX=%0.4f,deltaY=%0.4f,A.score[G][H]=%0.8f,myoutlierBC=%d,Sint()=%0.8f,newscore=%0.8f,score=%0.8f:tid=%d\n",
			 Yid,Xid,orientation,scaleID,I,J,G,H,deltaX,deltaY,FLAT_A(A.score,G,H),myoutlierBC,Sint(deltaX,deltaY,J-H,I-G,myoutlierBC),newscore,score,tid);
		  fflush(stdout);
		}
	      }

	      if(newscore > score){
		score = newscore;
		bestG = G;
		bestH = H;
	      }

	      if(OUTLIEREXTEND_FIX && myoutlierExtend){// NEW20 : try to merge with any interval that ends at G and H : just try to merge with previous aligned interval
		int G2 = FLAT_A(A.G,G,H);
		if(G2 > 0){
		  int H2 = FLAT_A(A.H,G,H);
		  if(DEBUG>=2 && !(IMIN <= G2 && G2 < G && JMIN[G2] <= H2 && H2 < H)){
                    #pragma omp critical
		    {
		      printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,Gmin=%d,Hmin=%d:G=%d,H=%d,G2=%d,H2=%d:tid=%d\n", Yid,Xid,orientation,scaleID,I,J,Gmin,Hmin,G,H,G2,H2,tid);
		      fflush(stdout);
		    }
		    assert(IMIN <= G2 && G2 < G && JMIN[G2] <= H2 && H2 < H);
		  }
		  if(!outlierExtendLimit || max(I-G2,J-H2) <= outlierExtendLimit){
		    PFLOAT deltaX = X[J]-X[H2];
		    PFLOAT deltaY = Y[I]-Y[G2];
		    // WAS107		    PFLOAT OutPen = OutPen2 + OutlierPenaltyBC[J-H2] + OutlierPenaltyBC[I-G2];
		    PFLOAT pscore = FLAT_A(A.score,G2,H2);
		    newscore = pscore + Sint(deltaX,deltaY,J-H2,I-G2,myoutlierBC);// WAS107 SintOutlier(deltaX,deltaY,OutPen);

		    if(EVERB && I==I_TRACE && J==J_TRACE /* && ((G==G_TRACE && H==H_TRACE) || newscore > score) */ && newscore > score){
                      #pragma omp critical
		      {
			printf("1:Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d,G2=%d,H2=%d:deltaX=%0.4f,deltaY=%0.4f,A.score[G2][H2]=%0.8f,Sint()=%0.8f,newscore=%0.8f,score=%0.8f:tid=%d\n",
			       Yid,Xid,orientation,scaleID,I,J,G,H,G2,H2,deltaX,deltaY,FLAT_A(A.score,G2,H2), Sint(deltaX,deltaY,J-H2,I-G2,myoutlierBC),newscore,score,tid);
			fflush(stdout);
		      }
		    }

		    if(newscore > score){
		      score = newscore;
		      bestG = G2;
		      bestH = H2;
		    }
		  }
		} // if(G2 > 0)
	      } // if(OUTLIEREXTEND_FIX && myoutlierExtend)
	    } //for(int G = Gmin; G < I; G++)
	  }// for(int H = J; --H >= Jmin;)

	  if(DEBUG>=2 && !(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J))){
	    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:I=%d,J=%d,score=%0.3f,bestG=%d,bestH=%d,localtype=%d:tid=%d\n",
		   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,score,bestG,bestH,localtype,tid);
	    fflush(stdout);
	    assert(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J));
	  }

	  if(myoutlierExtend){
	    /* try to merge with any previous outlier interval that end at I or J */	    

	    /* try outliers that end at G==I && H < J  (provided G > IMIN && H >= 2) */
	    int G = I;
	    if(G > IMIN){
	      int hmin = max(J-DELTA,max(JMIN[G],2));// max(Hmin,2)
	      for(int H = J; --H >= hmin;){
		if(DEBUG>=2) assert(H <= JMAX[G]);
		int G2 = FLAT_A(A.G,G,H);
		if(G2 > 0){
		  int H2 = FLAT_A(A.H,G,H);
		  if(DEBUG>=2) assert(IMIN <= G2 && G2 < G && JMIN[G2] <= H2 && H2 < H);

		  if(!outlierExtendLimit || max(I-G2,J-H2) <= outlierExtendLimit){
		    deltaX = X[J]-X[H2];
		    deltaY = Y[I]-Y[G2];
		    // WAS107		    PFLOAT OutPen = OutPen2 + OutlierPenaltyBC[J-H2] + OutlierPenaltyBC[I-G2];
		    PFLOAT pscore = FLAT_A(A.score,G2,H2);
		    newscore = pscore + Sint(deltaX,deltaY,J-H2,I-G2,myoutlierBC);// WAS107 SintOutlier(deltaX,deltaY,OutPen);

		    if(EVERB && I==I_TRACE && J==J_TRACE /* && ((G==G_TRACE && H==H_TRACE) || newscore > score) */){
                      #pragma omp critical
		      {
			printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d,G2=%d,H2=%d:deltaX=%0.4f,deltaY=%0.4f,A.score[G2][H2]=%0.8f,Sint()=%0.8f,newscore=%0.8f,score=%0.8f:tid=%d\n",
			       Yid,Xid,orientation,scaleID,I,J,G,H,G2,H2,deltaX,deltaY,FLAT_A(A.score,G2,H2), Sint(deltaX,deltaY,J-H2,I-G2,OutPen),newscore,score,tid);
			fflush(stdout);
		      }
		    }

		    if(newscore > score){
		      score = newscore;
		      bestG = G2;
		      bestH = H2;
		    }
		  }
		}
	      }
	    }

	    if(DEBUG>=2 && !(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J))){
	      #pragma omp critical
	      {
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:I=%d,J=%d,score=%0.3f,bestG=%d,bestH=%d,localtype=%d:tid=%d\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,score,bestG,bestH,localtype,tid);
		fflush(stdout);
		assert(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J));
	      }
	    }

	    /* try outliers that end at G < I && H==J (provided G > IMIN && H >= 2) */
	    int H = J;
	    if(H >= 2){
	      int gmin = max(I-DELTA,IMIN+1);
	      for(int G = I; --G >= gmin;){
		if(H > JMAX[G])
		  break;
		if(H < JMIN[G])
		  continue;
		int G2 = FLAT_A(A.G,G,H);
		if(G2 > 0){
		  int H2 = FLAT_A(A.H,G,H);
		  if(DEBUG>=2 && !(1 <= G2 && G2 < G && 1 <= H2 && H2 < H)){
                    #pragma omp critical
		    {
		      printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,Gmin=%d,Hmin=%d:G=%d,H=%d,G2=%d,H2=%d:tid=%d\n", Yid,Xid,orientation,scaleID,I,J,Gmin,Hmin,G,H,G2,H2,tid);
		      fflush(stdout);
		    }
		    assert(1 <= G2 && G2 < G && 1 <= H2 && H2 < H);
		  }
		  if(!outlierExtendLimit || max(I-G2,J-H2) <= outlierExtendLimit){
		    deltaX = X[J]-X[H2];
		    deltaY = Y[I]-Y[G2];
		    // WAS107		    PFLOAT OutPen = OutPen2 + OutlierPenaltyBC[J-H2] + OutlierPenaltyBC[I-G2];
		    newscore = FLAT_A(A.score,G2,H2) + Sint(deltaX,deltaY,J-H2,I-G2,myoutlierBC); // WAS107 SintOutlier(deltaX,deltaY,OutPen);

		    if(EVERB && I==I_TRACE && J==J_TRACE /* && ((G==G_TRACE && H==H_TRACE) || newscore > score) */){
                      #pragma omp critical
		      {
			printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d,G2=%d,H2=%d:deltaX=%0.4f,deltaY=%0.4f,A.score[G2][H2]=%0.8f,Sint()=%0.8f,newscore=%0.8f,score=%0.8f:tid=%d\n",
			       Yid,Xid,orientation,scaleID,I,J,G,H,G2,H2,deltaX,deltaY,FLAT_A(A.score,G2,H2), Sint(deltaX,deltaY,J-H2,I-G2,myoutlierBC),newscore,score,tid);
			fflush(stdout);
		      }
		    }

		    if(newscore > score){
		      score = newscore;
		      bestG = G2;
		      bestH = H2;
		    }
		  }
		}
	      }
	    }

	    if(DEBUG>=2 && !(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]))){
	      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:I=%d,J=%d,score=%0.3f,bestG=%d,bestH=%d,localtype=%d:tid=%d\n",
		     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,score,bestG,bestH,localtype,tid);
	      fflush(stdout);
	      assert(bestG < 0 || (1 <= bestG && bestG <= I && 1 <= bestH && bestH <= J && bestH <= JMAX[bestG]));
	    }

	    if(bestG >= 1){

	      PFLOAT pscore = FLAT_A(A.score, bestG, bestH);

	      if(VERB>=3 && Yid==1 && Xid==0 && orientation==0 && I==14938 && J==38){
		printf("A:Yid=%d,Xid=%d:I=%d,J=%d,bestG=%d,bestH=%d,score=%0.6f (pscore=%0.6f,OutPen=%0.6f)\n",
		       Yid,Xid,I,J,bestG,bestH,score,pscore,OutPen);
		fflush(stdout);
	      }

	      if(score <= pscore + SintOutlier(X[J] - X[bestH], Y[I] - Y[bestG], OutPen2 + OutlierPenaltyBC[J - bestH] + OutlierPenaltyBC[I - bestG])
		 /* && alignFP(A, bestG, bestH, localtype,  Y, X, N, M, orientation, Yid, Xid, delY, delX, pscore, 0,JMIN,JMAX,IMIN,IMAX) > LogPvThreshold */){ /* try to extend outlier interval to the left by reducing bestG,bestH by up to myoutlierExtend sites */

		int origbestG = bestG, origbestH = bestH;

		for(int cnt = 0; cnt < 100; cnt++){

		  int Gmin2 = max(IMIN, origbestG - myoutlierExtend);
		  int Hmin2 = max(JMIN[IMIN], origbestH - myoutlierExtend);
		  if(outlierExtendLimit){
		    Gmin2 = max(I - outlierExtendLimit, Gmin2);
		    Hmin2 = max(J - outlierExtendLimit, Hmin2);
		  }
		  for(int G = origbestG; G >= Gmin2; G--){
		    PFLOAT deltaY = Y[I] - Y[G];
		    PFLOAT *AscoreG = &FLAT_A(A.score,G,0);
		    int hmax = (G==origbestG) ? min(origbestH-1, JMAX[G]) : min(origbestH,JMAX[G]);
		    int hmin = max(Hmin2, JMIN[G]);
		    for(int H = hmax; H >= hmin; H--){
		      if(DEBUG>=2) assert(G < origbestG || H < origbestH);
		      if(DEBUG>=2) assert(JMIN[G] <= H && H <= JMAX[G]);
		      PFLOAT deltaX = X[J] - X[H];
		      // WAS107		      PFLOAT Pen = OutPen + OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-G];
		      PFLOAT newscore = AscoreG[H] + Sint(deltaX,deltaY,J-H,I-G,myoutlierBC);// WAS107 SintOutlier(deltaX, deltaY, Pen);

		      if(EVERB && I==I_TRACE && J==J_TRACE /* && ((G==G_TRACE && H==H_TRACE) || newscore > score)*/){
                        #pragma omp critical
			{
			  printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,cnt=%d,G=%d,H=%d:deltaX=%0.4f,deltaY=%0.4f,A.score[G][H]=%0.8f,Sint()=%0.8f,newscore=%0.8f,score=%0.8f(G=%d,H=%d):tid=%d\n",
				 Yid,Xid,orientation,scaleID,I,J,cnt,G,H,deltaX,deltaY,FLAT_A(A.score,G,H), Sint(deltaX,deltaY,J-H,I-G,myoutlierBC),newscore,score,bestG,bestH,tid);
			  fflush(stdout);
		        }
		      }

		      if(newscore > score){
			if(VERB>=2 && (Y[I]-Y[bestG] > X[J]-X[bestH] ? H < origbestH : G < origbestG)){
			  printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,cnt=%d,G=%d->%d,H=%d->%d:deltaX=%0.4f->%0.4f,deltaY=%0.4f->%0.4f,A.score[G][H]=%0.8f,Sint()=%0.8f,newscore=%0.8f,score=%0.8f(G=%d,H=%d):(Outlier not extended)\n",
				 Yid,Xid,orientation,scaleID,I,J,cnt,origbestG,G,origbestH, H,X[J]-X[origbestH],deltaX,Y[I]-Y[origbestG],deltaY,
				 FLAT_A(A.score,G,H), Sint(deltaX,deltaY,J-H,I-G,myoutlierBC),newscore,score,bestG,bestH);
			  fflush(stdout);
			}

			score = newscore;
			bestG = G;
			bestH = H;
		      }
		    }
		  }

		  if(DEBUG>=2 && !(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]))){
		    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:I=%d,J=%d,score=%0.3f,bestG=%d,bestH=%d,localtype=%d:tid=%d\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,score,bestG,bestH,localtype,tid);
		    fflush(stdout);
		    assert(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]));
		  }
		  if(DEBUG>=2) assert(bestG <= origbestG && bestH <= origbestH);

		  if(VERB>=3 && Yid==1 && Xid==0 && orientation==0 && I==14938 && J==38 && bestG >= 1){
		    printf("B:  cnt=%d:Yid=%d,Xid=%d:I=%d,J=%d,bestG=%d,bestH=%d,score=%0.6f (pscore=%0.6f,OutPen=%0.6f)\n",
			   cnt,Yid,Xid,I,J,bestG,bestH,score,FLAT_A(A.score,bestG,bestH),OutPen);
		    fflush(stdout);
		  }

		  /* also check if best interval before (G,H) can be merged with this outlier interval : include checking G==I (if H <= min(J-1,origbestH)) and H==J (with G <= min(I-1,origbestG)) */
		  int gmin2 = max(IMIN+1,Gmin2), hmin2 = max(2,Hmin2);
		  for(int G = I; G >= gmin2; G--){
		    if(G < I) G = min(origbestG,G);
		    int hmin = max(hmin2,JMIN[G]);
		    int Hmax2 = (G <= min(I-1,origbestG)) ? J : min(J-1,origbestH);

		    for(int H = min(Hmax2,JMAX[G]); H >= hmin; H--){
		      if(H < J) H = min(origbestH,H);
		      if(DEBUG>=2) assert(G <= min(I-1,origbestG) || H <= min(J-1,origbestH));
		      int G2 = FLAT_A(A.G, G, H);
		      if(G2 > 0){
			int H2 = FLAT_A(A.H, G, H);
			if(DEBUG>=2 && !(JMIN[G2] <= H2 && H2 <= JMAX[G2] && H2 <= H && IMIN <= G2 && G2 <= G && (G2 < G || H2 < H))){
			  printf("Yid=%d(%lld),Xid=%d(%lld):I=%d,J=%d:G=%d,G2=%d,H=%d,H2=%d,origbestG=%d,origbestH=%d\n",
				 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,I,J,G,G2,H,H2,origbestG,origbestH);
			  fflush(stdout);
			  assert(JMIN[G2] <= H2 && H2 <= JMAX[G2] && H2 <= H && IMIN <= G2 && G2 <= G && (G2 < G || H2 < H));
			}
			if(!outlierExtendLimit || max(I-G2,J-H2) <= outlierExtendLimit){
			  newscore = FLAT_A(A.score, G2, H2) + Sint(X[J]-X[H2],Y[I]-Y[G2],J-H2,I-G2,myoutlierBC);// WAS107 SintOutlier(X[J]-X[H2],Y[I]-Y[G2],OutPen + OutlierPenaltyBC[J - H2] + OutlierPenaltyBC[I - G2]);
			  if(newscore > score + 1e-6){
			    if(EVERB && I==I_TRACE && J==J_TRACE /* && G==G_TRACE && H==H_TRACE */){
			      printf("Yid=%d(%lld),Xid=%d(%lld):I=%d,J=%d:cnt=%d,G=%d,G2=%d,H=%d,H2=%d,bestG=%d,origbestG=%d,bestH=%d,origbestH=%d:score=%0.6f->%0.6f(Sint=%0.6f):tid=%d\n",
				     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,I,J,cnt,G,G2,H,H2,bestG,origbestG,bestH, origbestH,score,newscore,
				     Sint(X[J]-X[H2],Y[I]-Y[G2],J-H2,I-G2,myoutlierBC),tid);
			      fflush(stdout);
			    }
			    score = newscore;
			    bestG = G2;
			    bestH = H2;
			  }
			}
		      }
		    }
		  }

		  if(DEBUG>=2 && !(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]))){
		    #pragma omp critical
		    {
		      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:I=%d,J=%d,score=%0.3f,bestG=%d,bestH=%d(cnt=%d,G=%d..%d,H=%d..%d):tid=%d\n",
			     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,score,bestG,bestH,cnt,origbestG,Gmin2,origbestH,Hmin2,tid);
		      fflush(stdout);
		      assert(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]));
		    }
		  }

		  if(origbestG == bestG && origbestH == bestH)
		    break;

		  if(VERB>=3 && Yid==1 && Xid==0 && orientation==0 && I==14938 && J==38 && bestG >= 1){
		    printf("C:  cnt=%d:Yid=%d,Xid=%d:I=%d,J=%d,bestG=%d,bestH=%d,score=%0.6f (pscore=%0.6f,OutPen=%0.6f)\n",
			   cnt,Yid,Xid,I,J,bestG,bestH,score,FLAT_A(A.score,bestG,bestH),OutPen);
		    fflush(stdout);
		  }

		  if(DEBUG>=2 && !(bestG <= origbestG || bestH <= origbestH)){
		    printf("Yid=%d(%lld),Xid=%d(%lld):I=%d,J=%d:origbestG=%d,bestG=%d,origbestH=%d,bestH=%d\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,I,J,origbestG,bestG,origbestH,bestH);
		    fflush(stdout);
		    assert(bestG <= origbestG || bestH <= origbestH);
		  }

		  origbestG = bestG;
		  origbestH = bestH;
		}/* for(int cnt = 0; cnt < 100; cnt++) */
	      } else {/* try to merge non-outlier interval with previous interval (in case previous interval is an outlier) */
		int G = FLAT_A(A.G,bestG,bestH);
		if(G > 0){
		  int H = FLAT_A(A.H,bestG,bestH);
		  if(!outlierExtendLimit || max(I-G,J-H) <= outlierExtendLimit){
		    PFLOAT deltaX = X[J]-X[H], deltaY = Y[I]-Y[G];
		    // WAS107		    PFLOAT Pen = OutPen + OutlierPenaltyBC[I-G] + OutlierPenaltyBC[J-H];
		    PFLOAT newscore = FLAT_A(A.score,G,H) + Sint(deltaX,deltaY,J-H,I-G,myoutlierBC);// WAS107 SintOutlier(deltaX, deltaY, Pen);

		    if(EVERB && I==I_TRACE && J==J_TRACE /* && ((G==G_TRACE && H==H_TRACE) || newscore > score)*/){
                      #pragma omp critical
		      {
			printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d:deltaX=%0.4f,deltaY=%0.4f,A.score[G][H]=%0.8f,Sint()=%0.8f,newscore=%0.8f,score=%0.8f(G=%d,H=%d):tid=%d\n",
			       Yid,Xid,orientation,scaleID,I,J,G,H,deltaX,deltaY,FLAT_A(A.score,G,H),Sint(deltaX,deltaY,J-H,I-G,myoutlierBC),newscore,score,bestG,bestH,tid);
			fflush(stdout);
		      }
		    }

		    if(newscore > score){
		      score = newscore;
		      bestG = G;
		      bestH = H;
		      if(DEBUG>=2 && !(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]))){
			printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:I=%d,J=%d,score=%0.3f,bestG=%d,bestH=%d,localtype=%d:tid=%d\n",
			       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,score,bestG,bestH,localtype,tid);
			fflush(stdout);
			assert(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]));
		      }
		    }
	          }
		}
	      }
	    }
	  }

	  if(DEBUG>=2) assert(isfinite(score));
	  if(DEBUG>=2 && !(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]))){
	    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:I=%d,J=%d:bestG=%d,bestH=%d,score=%0.6f:tid=%d\n",
		   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,bestG,bestH,score,tid);
	    fflush(stdout);
	    assert(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]));
	  }

	  if(EVERB && I== I_TRACE && J== J_TRACE){
            printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d:score=%0.3f,bestG=%d,bestH=%d\n",Yid,Xid,orientation,scaleID,I,J,score,bestG,bestH);
	    fflush(stdout);
	  }

	  FLAT_A(A.score, I, J) = score;
	  FLAT_A(A.H, I, J) = bestH;
	  FLAT_A(A.G, I, J) = bestG;// NOTE : write bestG last

	} else {/* no need to track bestG,bestH */

	  int Gmin = max(I-DELTA,IMIN), Hmin = max(J-DELTA,JMIN[Gmin]);
	  PFLOAT deltaX, deltaY, newscore;

	  if(DEBUG>=1+RELEASE) assert(!myoutlierBC);// NOTE : uncomment if this assertion is triggered : myoutlierBC is handled in the loop below

	  for(int H= J; --H >= Hmin;){
	    if(DEBUG>=2 && !(1 <= H && H < J && J-H <= DELTA && J <= M)){
	      #pragma omp critical
	      {
	        printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,H=%d,M=%d (Jmin=%d,Jmax=%d)\n",Yid,Xid,orientation,scaleID,I,J,H,M,Jmin,Jmax);
		fflush(stdout);
		assert(1 <= H && H < J && J-H <= DELTA && J <= M);
              }
	    }
	    if(DEBUG>=2 && ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT))) assert(&FLAT_DELX(delX, J-H, J, M) == &delX[J-H][J]);
	    deltaX = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];// X[J] - X[H]
	    if(DEBUG>=2) assert(fabs(X[J] - X[H] - deltaX) < 0.001);
	    int Gmin2 = max(Gmin,Imin[H]);
	    int Gmax2 = min(I-1,Imax[H]);
	    for(int G= Gmin2; G <= Gmax2; G++){
	      deltaY = delYI[G-I+DELTA];// Y[I] - Y[G] 
	      if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	      newscore = FLAT_A(A.score, G, H) + Sint(deltaX,deltaY,J-H,I-G,myoutlierBC);
	      if(EVERB && I==I_TRACE && J==J_TRACE /* && G==I_TRACE-1 && H==J_TRACE-1 */){
	        #pragma omp critical
		{
		  printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d:deltaX=%0.4f,deltaY=%0.4f,A.score[G][H]=%0.8f,myoutlierBC=%d,Sint()=%0.8f,newscore=%0.8f,score=%0.8f\n",
			 Yid,Xid,orientation,scaleID,I,J,G,H,deltaX,deltaY,FLAT_A(A.score,G,H),myoutlierBC,Sint(deltaX,deltaY,J-H,I-G,myoutlierBC),newscore,score);
		  fflush(stdout);
		}
	      }
	      
	      score = PMAX(score,newscore);
	    }
	  }

	  if(EVERB && ((I==I_TRACE && J==J_TRACE) || (I == G_TRACE && J == H_TRACE))){
	    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:I=%d,J=%d,A.score[I][J]= %0.6f(%p) -> %0.6f,FLAT_A(A.score,I,J)= %0.6f(%p),A.score=%p,tid=%d\n",
		   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,A.score[I][J],&A.score[I][J],score,FLAT_A(A.score,I,J),&FLAT_A(A.score,I,J),A.score,tid);
	    fflush(stdout);
          }

	  FLAT_A(A.score, I, J) = score;
	}
      } // omp for J = Jmin .. Jmax
    } // general case without (DEFER_BP && !myoutlierExtend)
  } // for I = 2 .. N

  if(VERB>=3 && Yid==69 && Xid==70 && orientation==0 && N >= 19 && M >= 48103){
    #pragma omp critical
    {
      printf("C:tid=%d:Yid=%d,Xid=%d,or=%d,sc=%d:N=%d,M=%d:FLAT_A(A.Rijy,19,48103)= %d (%p)\n",
	     tid,Yid,Xid,orientation,scaleID,N,M,FLAT_A(A.Rijy,19,48103), &FLAT_A(A.Rijy,19,48103));
      fflush(stdout);
    }
  }

  /* locate best alignment result */
  /* If ENDINIT only check results where I >= N+1-DELTA OR J >= M+1-DELTA */
  PFLOAT bestscore = ScoreThreshold; // MINSCORE;
  int bestI = -1, bestJ = -1;
  int bestR = localtype;/* right end alignment is local */

  if(C > ((PFLOAT)MINSCORE)){/* consider right end local alignments */
    for(int I= IMIN; I <= IMAX; I++){
      for(int J = JMIN[I]; J <= JMAX[I]; J++){
        PFLOAT score = FLAT_A(A.score, I, J) + C;
	if(score > bestscore){
	  bestscore = score;
	  bestI = I;
	  bestJ = J;
	}
      }
    }
  }
  
  //prefetch_array(Xr, sizeof(*Xr)*4, _MM_HINT_T1);
  {
    for(int I = IMAX; I >= IMIN; I--){
      int Jmin = (!ENDINIT) ? JMIN[I] : (I >= N+1-DELTA) ? JMIN[I] : max(JMIN[I],M+1-DELTA);
      int Jmax = JMAX[I];
      PFLOAT Ydiff = Yr[I];// Y[N+1]-Y[I];
      int * RijxI = &FLAT_A(A.Rijx, I, 0), *RijyI = &FLAT_A(A.Rijy, I, 0);

      for(int J= Jmin; J <= Jmax; J++){
        PFLOAT Ascore = FLAT_A(A.score, I, J);
	PFLOAT thresh;
	if(SEND_THRESH && EndThresh < (thresh = bestscore - Ascore))
	  continue;

	int Rijx = RijxI[J], Rijy = RijyI[J];
	if(DEBUG>=2 && !(Rijx > M || Rijy > N)){
          printf("Yid=%d(%lld),Xid=%d(%lld):I=%d,J=%d:Rijx=%d,Rijy=%d,M=%d,N=%d,DELTA=%d\n",
		 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,I,J,Rijx,Rijy,M,N,DELTA);
	  printf("  Y[I]=%0.8f,Y[N+1]=%0.8f,X[J]=%0.8f,X[M+1]=%0.8f\n",Y[I],Y[N+1],X[J],X[M+1]);
	  fflush(stdout);
	  assert(Rijx > M || Rijy > N);
	}
	if((Rijy > N && YrightClipped) || (Rijx > M && XrightClipped))
	  continue;

	//      FLOAT score = Send(PMIN(X[M+1]-X[J], Ydiff), min(M, Rijx)+1-J, min(N, Rijy)+1-I);
	PFLOAT score = Send(PMIN(Xr[J], Ydiff),(PFLOAT)(min(M, Rijx)+1-J + min(N, Rijy)+1-I));
	if(ENDFIX){/* update Rijx,Rijy based on Sbnd() */
	  if(SEND_THRESH && score < thresh - EndThresh2)
	    continue;
    
	  if(Rijx > M) {
	    double a = M+1-I-J;
	    int K;
	    for(K= Rijy; K > I; K--){
	      PFLOAT bound = Sbnd(Xr[J],Y[K]-Y[I], K+a);
	      if(bound <= score)
		break;// HERE : see ENDFIX2 in refalign.cpp
	      score = bound;
	    }
	    RijyI[J] = Rijy = K;
	  }
	  if(Rijy > N) {
	    double a = N+1-I-J;
	    int K;
	    for(K= Rijx;K > J;K--){
	      PFLOAT bound = Sbnd(X[K]-X[J], Ydiff,K+a);
	      if(bound <= score)
		break;// HERE : see ENDFIX2 in refalign.cpp
	      score = bound;
            }
	    RijxI[J] = Rijx = K;
          }
        }

	Ascore += score;
	if(VERB>=2 && orientation && I==50 && J==13){
          printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f(p->score=%0.4f)\n",
            Yid,YYmap[Yid]->id, Xid,XXmap[Xid]->id, orientation, scaleID,I, J, Ascore, FLAT_A(A.score, I, J));
	  fflush(stdout);
        }
	if(Ascore > bestscore){
          bestscore = Ascore;
	  bestI = I;
	  bestJ = J;
	  bestR = -1;/* normal right end */
        }
      }
    }
  }

  if(VERB>=3 && Yid==69 && Xid==70 && orientation==0 && N >= 19 && M >= 48103){
    #pragma omp critical
    {
      printf("D:tid=%d:Yid=%d,Xid=%d,or=%d,sc=%d:N=%d,M=%d:FLAT_A(A.Rijy,19,48103)= %d (%p)\n",
	     tid,Yid,Xid,orientation,scaleID,N,M,FLAT_A(A.Rijy,19,48103), &FLAT_A(A.Rijy,19,48103));
      fflush(stdout);
    }
  }

  if(bestI < 0){/* no valid alignment found : should never happen with -S -1e30 */
    if(VERB >= (EVERB ? 1 : 2)){
      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:No alignment found (bestscore=%0.3f,bestI=%d)\n",
	     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestscore,bestI);
      fflush(stdout);
    }
    if(mem)
      free(mem);
    if(MEMCHECK && RA_MIN_TIME && RA_MIN_MEM < 8) ReleaseMemory(Mem,nMemSiz,tid, Yid, Xid, N, M, Fmem);

    return;
  }

  if(VERB >= (EVERB ? 1 : 2)){
    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,bestscore=%0.3f\n",
	   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,bestscore);
    fflush(stdout);
  }

  if(DEFER_BP){/* compute G,H values along best alignment path */
    int I = bestI;
    int J = bestJ;
    if(DEBUG>=2) assert(bestI > 0 && bestJ > 0);

    int bestG = -1, bestH = -1;

    PFLOAT score;
    while(I >= 0){
      if(DEBUG>=2) assert(I > 0 && J > 0);
      if(myoutlierExtend){/* A.G and A.H has already been computed */
	score = FLAT_A(A.score,I,J);
	bestG = FLAT_A(A.G,I,J);
	bestH = FLAT_A(A.H,I,J);

	if(DEBUG>=2 && !(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]))){
	  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,bestscore=%0.3f:I=%d,J=%d:score=%0.3f,bestG=%d,bestH=%d,localtype=%d\n",
		 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,bestscore,I,J,score,bestG,bestH,localtype);
	  fflush(stdout);
	  assert(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]));
	}
      } else {
	/* update A.G[I][J] & A.H[I][J] by repeating recurrance computation */
	score = FLAT_A(A.Uscore, I, J);
	if(DEBUG>=2) assert(isfinite(score));
	bestG = localtype;

	int Gmin = max(I-DELTA, IMIN);
	for(int G= I; --G >= Gmin;){
	  PFLOAT deltaY = delY[I][G-I+DELTA];// Y[I] - Y[G] 
	  if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	  PFLOAT *AscoreG = &FLAT_A(A.score, G, 0);
	  // WAS107	  PFLOAT a = J+I-G;
	  
	  int Hmin = max(J-DELTA, JMIN[G]);
	  int Hmax = min(J-1, JMAX[G]);

	  for(int H= Hmax; H >= Hmin; H--){
	    if(DEBUG>=2 && !(JMIN[G] <= H && H < J && H <= JMAX[G] && J-H <= DELTA && J <= M)){
	      #pragma omp critical
	      {
	        printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d,M=%d\n",Yid,Xid,orientation,scaleID,I,J,G,H,M);
		fflush(stdout);
		assert(JMIN[G] <= H && H < J && H <= JMAX[G] && J-H <= DELTA && J <= M);
              }
	    }
	    PFLOAT deltaX = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];// X[J] - X[H]
	    if(DEBUG>=2) assert(fabs(X[J] - X[H] - deltaX) < 0.001);
	    PFLOAT newscore = AscoreG[H] +  Sint(deltaX,deltaY,J-H,I-G,myoutlierBC); // WAS107 Sint(deltaX,deltaY,(PFLOAT)(a-H));
	    if(EVERB && I== I_TRACE && J== J_TRACE /* && ((G== G_TRACE && H== H_TRACE) || newscore > score) */){
	      #pragma omp critical
	      {
		printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d:deltaX=%0.4f,deltaY=%0.4f,A.score[G][H]=%0.8f,Sint()=%0.8f,newscore=%0.8f,score=%0.8f\n",
		       Yid,Xid,orientation,scaleID,I,J,G,H,deltaX,deltaY,AscoreG[H],Sint(deltaX,deltaY,J-H,I-G,myoutlierBC),newscore,score);
		fflush(stdout);
	      }
	    }
	    if(newscore > score){
	      score = newscore;
	      bestG = G;
	      bestH = H;
	    }
	  }
	}
	if(DEBUG>=2 && !(fabs(score - A.score[I][J]) < (PFLOAT)(USE_PFLOAT ? (bestscore > 0.0 ? 1e-3 : 1e-2) :1e-8))){
          #pragma omp critical
	  {
	    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),orientation=%d,sc=%d:I=%d,J=%d,G=%d,H=%d:score=%0.8f,A.score[I][J]=%0.8f,bestscore=%0.3f,myoutlierExtend=%d\n",
		   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,I,J,bestG,bestH,score,A.score[I][J],bestscore,myoutlierExtend);
	    fflush(stdout);
	    assert(fabs(score - FLAT_A(A.score, I, J)) < (PFLOAT)(USE_PFLOAT ? (bestscore > 0.0 ? 1e-3 : 1e-2):1e-8));
	  }
	}
      }

      if(bestG < 0){/* value of bestH is undefined */
	if(DEFER_G>=1 && score > C){
	  bestG = -1;// why not localtype ???
	  if(DEFER_G >= 2){/* check if actual value of Lijx,Lijy needs to be incremented */
	    if(DEBUG>=2) assert(I <= DELTA || J <= DELTA);
	    int Lijx = FLAT_A(A.Lijx,I,J);
	    int Lijy = FLAT_A(A.Lijy,I,J);
	    if(DEBUG>=2) assert(0 <= Lijx && Lijx <= J);
	    if(DEBUG>=2) assert(0 <= Lijy && Lijy <= I);
	    if(DEBUG>=2) assert(Lijx <= 0 || Lijy <= 0);
	    if((Lijx <= 0) ? !XleftClipped : !YleftClipped){
	      PFLOAT SmIJ = Sm(J,I);
	      PFLOAT ascore = C;

	      if(Lijx <= 0){
		ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I+1-max(1,Lijy))) + SmIJ;
		ascore = max(ascore,C);
		for(int K= Lijy; K < I; ) {
		  PFLOAT bound = Sbnd(X[J],Y[I]-Y[K],(PFLOAT)(J+I-K)) + SmIJ;
		  if(bound <= ascore)
		    break;
		  ascore = bound;
		  Lijy = ++K;
		}
		FLAT_A(A.Lijy, I, J) = Lijy;
	      } else {
		if(DEBUG>=2) assert(Lijy <= 0);
		ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-max(1,Lijx)+I)) + SmIJ;
		ascore = max(ascore,C);
		for(int K= Lijx; K < J; ){
		  PFLOAT bound = Sbnd(X[J]-X[K],Y[I],(PFLOAT)(J+I-K)) + SmIJ;
		  if(bound <= ascore)
		    break;
		  ascore = bound;
		  Lijx = ++K;
		}
		FLAT_A(A.Lijx, I, J) = Lijx;
	      }
	    }
	  }
	}
	FLAT_A(A.G, I, J) = bestG;
	if(DEBUG>=2 && !(Ch > 0.0 || PoutlierEnd > 0.0)) assert(A.G[I][J] > -2);
	if(DEBUG>=2 && !(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]))){
	  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,bestscore=%0.3f:I=%d,J=%d,score=%0.3f,bestG=%d,bestH=%d,localtype=%d\n",
		 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,bestscore,I,J,score,bestG,bestH,localtype);
	  fflush(stdout);
	  assert(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]));
	}
      } else {
	if(DEBUG>=2) assert(bestG > 0 && bestH > 0);
	FLAT_A(A.G, I, J) = bestG;
	FLAT_A(A.H, I, J) = bestH;
	if(DEBUG>=2 && !(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]))){
	  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,bestscore=%0.3f:I=%d,J=%d:score=%0.3f,bestG=%d,bestH=%d,localtype=%d\n",
		 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,bestscore,I,J,score,bestG,bestH,localtype);
	  fflush(stdout);
	  assert(bestG < 0 || (IMIN <= bestG && bestG <= I && JMIN[bestG] <= bestH && bestH <= J && bestH <= JMAX[bestG]));
	}
	J = bestH;
      }
      I = bestG;
    }
  }

  if(VERB>=3 && Yid==69 && Xid==70 && orientation==0 && N >= 19 && M >= 48103){
    #pragma omp critical
    {
      printf("D:tid=%d:Yid=%d,Xid=%d,or=%d,sc=%d:N=%d,M=%d:FLAT_A(A.Rijy,19,48103)= %d (%p)\n",
	tid,Yid,Xid,orientation,scaleID,N,M,FLAT_A(A.Rijy,19,48103), &FLAT_A(A.Rijy,19,48103));
      fflush(stdout);
    }
  }

  if(EVERB && !(bestscore  > align->score)){
    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestscore=%0.4f,bestI=%d,bestJ=%d:align->score=%0.4f,align->logPV=%0.2f,align->numpairs=%d\n",
	Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestscore,bestI,bestJ,align->score,align->logPV,align->numpairs);
    fflush(stdout);
  }

  double origscore = align->score;

  if(bestscore > origscore || EVERB){ /* backtrack through array to retrieve best alignment starting at right end A[bestI][bestJ] */
    if(DEBUG>=2)
      for(int U = 0; U <= M; U++)
	outscore[U] = iscore[U] = aNaN;

    int I = bestI;
    int J = bestJ;

    int K = 0;
    
    iscore[0] = bestscore - FLAT_A(A.score, I, J);
    outscore[0] = iscore[0];

    if(EVERB){
      #pragma omp critical
      {
        printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,bestR=%d,bestscore=%0.6f:iscore[0]=%0.6f\n",
	       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,bestR,bestscore,iscore[0]);
	fflush(stdout);
      }
    }

    while(I >= 0){
      if(DEBUG>=2){
	assert(I >= 1 && I <= N);
	assert(J >= 1 && J <= M);
	assert(K < M);
      }

      Ilist[K] = I;
      Jlist[K] = J;
      int G = FLAT_A(A.G, I, J);
      int H = FLAT_A(A.H, I, J);
      if(G > 0){
	if(DEBUG>=2) assert(IMIN <= G && H >= JMIN[G] && H <= JMAX[G]);
	iscore[K+1] = FLAT_A(A.score, I, J) - FLAT_A(A.score, G, H);
	PFLOAT y = Y[I] - Y[G];
	PFLOAT x = X[J] - X[H];
	if(Poutlier > 0.0){
	  PFLOAT Bias,Pen;
	  SintDetail(x, y, J - H, I - G, J, I, Bias, Pen);
	  outscore[K+1] = OutlierBias + Bias + Pen;
	  if(EVERB || DEBUG>=2){
	    PFLOAT OutPen = OutlierPenalty;
	    if(myoutlierBC)
	      OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-G];
	    PFLOAT iscoreB = OutlierBias + (OUTLIER_DELTA(x-y) ? max(Pen + Bias, OutPen + Bias*biasWToutlierF) : Pen + Bias);	    

	    if(DEBUG>=2 && !(fabs(iscoreB - iscore[K+1]) < 0.001)){
	      #pragma omp critical
	      {
	        printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d:OutlierBias=%0.6f,OutlierPenalty=%0.6f:A.score[I,J]=%0.6f,A.score[G,H]=%0.6f\n",
	               Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,OutlierBias,OutlierPenalty,FLAT_A(A.score,I,J),FLAT_A(A.score,G,H));
	        printf("  K=%d:I=%d,J=%d,G=%d,H=%d:x=%0.4f,y=%0.4f,OutlierBias=%0.6f,Bias=%0.6f,Pen=%0.6f,OutPen=%0.6f,iscore[K+1]=%0.6f,outscore[K+1]=%0.6f,iscoreB=%0.6f\n",
		       K,I,J,G,H,x,y,OutlierBias,Bias,Pen,OutPen,iscore[K+1],outscore[K+1],iscoreB);
	        fflush(stdout);
		assert(fabs(iscoreB - iscore[K+1]) < 0.001);
	      }
	    }
	  }
	  if(DEBUG>=2) assert(iscore[K+1] >= outscore[K+1] - 0.001);
	} else
	  outscore[K+1] = iscore[K+1];
      } else {/* left end segment */
	iscore[K+1] = FLAT_A(A.score, I, J);
	outscore[K+1] = iscore[K+1];
      }
      K++;
      if(EVERB>=2){
	if(G > 0){
	  PFLOAT y = Y[I] - Y[G];
	  PFLOAT x = X[J] - X[H];
	  printf("K=%d: I=%d,J=%d,G=%d,H=%d:y= %0.4f, x=%0.4f, iscore[K]= %0.6f, outscore[K]= %0.6f\n",K,I,J,G,H,y,x,iscore[K],outscore[K]);
	} else
	  printf("K=%d: I=%d,J=%d,G=%d,H=%d: iscore[K]= %0.6f, outscore[K]= %0.6f\n",K,I,J,G,H,iscore[K],outscore[K]);	  
      }
      I = G;
      J = H;
    }
    if(DEBUG) assert(K>0);

    if(VERB>=3 && Yid==69 && Xid==70 && orientation==0 && N >= 19 && M >= 48103){
      #pragma omp critical
      {
	printf("E:tid=%d:Yid=%d,Xid=%d,or=%d,sc=%d:N=%d,M=%d:FLAT_A(A.Rijy,19,48103)= %d (%p)\n",
	       tid,Yid,Xid,orientation,scaleID,N,M,FLAT_A(A.Rijy,19,48103), &FLAT_A(A.Rijy,19,48103));
	fflush(stdout);
      }
    }

    double bestLogPV = alignFP(A,bestI,bestJ,bestR,Y,X,N,M,orientation,scaleID,Yid,Xid,delY,delX,bestscore, /* (Yid==285 && Xid==22611 && orientation==0) ? 1 : 0 */ 0, 1, myoutlierExtend, myoutlierBC,tid,
			       JMIN,JMAX,IMIN,IMAX);

    if(PAIRMERGE_OUTLIEREXTEND_FIX && PairMerge && PoutlierEnd > 0.0 && outlierExtend && !myoutlierExtend && bestscore > ScoreThreshold && bestLogPV > LogPvThreshold && K >= AlignedSiteThreshold){ /* repeat alignment with outlierExtend */
      if(EVERB){
        printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestscore=%0.4f,bestI=%d,bestJ=%d,bestLogPV=%0.2f,pairs=%d:align->score=%0.4f,align->logPV=%0.2f,align->numpairs=%d\n",
		 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestscore,bestI,bestJ,bestLogPV,K,align->score,align->logPV,align->numpairs);
        fflush(stdout);
      }

      if(mem)   /* free local memory */
	free(mem);
      if(MEMCHECK && RA_MIN_TIME && RA_MIN_MEM < 8) ReleaseMemory(Mem,nMemSiz,tid, Yid, Xid, N, M, Fmem);

      return pairalignXY(Y,N,X,M,A,align,Yid,Xid,delY,delYx,delX,phash,orientation, scaleID, outlierExtend, outlierBC,tid,Fmem,Mem,StrideMax);
    }

    if(PAIRMERGE_FIX && PairMerge /* && !(SplitSegDup > 0.0 && !Cfirst && PairMergeHmap <= 0) */){

      /* compute size of largest internal outlier and endoutlier */
      double Xscale = 1.0;
      double Yscale = 1.0;
      if(!(spots_filename[0] || (svcheck && numrefmaps > 0)) && fabs(PixelLen - origPixelLen) > 1e-12 /* NEW */)
	Yscale = Xscale = origPixelLen/PixelLen;

      int LI = Ilist[K-1];
      int LJ = Jlist[K-1];
      int RI = Ilist[0];
      int RJ = Jlist[0];

      double overlap = (Y[RI]-Y[LI] + X[RJ]-X[LJ]) * 0.5;
      double leftEnd = (I >= -1) ? 0.0 : min(Y[LI], X[LJ]);
      int leftEndN = (I >= -1) ? 0 : (Y[LI] < X[LJ]) ? LI-1 : LJ-1;
      double rightEnd = (bestR >= -1) ? 0.0 : min(Y[N+1]-Y[RI], X[M+1] - X[RJ]);
      int rightEndN = (bestR >= -1) ? 0 : (Y[N+1]-Y[RI] < X[M+1] - X[RJ]) ? N-RI : M-RJ;
      double maxOutlier = 0.0;
      int maxOutlierM = 0;
      double overlapExclude = 0.0;

      if((PairMergeMaxOutlier >= 0.0 || PairMergeExcludeOutliers) && !(overlap > (PairMergeOverlapped >= 0.0 ? PairMergeOverlapped : PairMerge))){/* locate largest outlier (and subtract outliers from overlap) */
	int lastI = LI, I;
	int lastJ = LJ, J;

	for(int k = 1; k < K; lastI = I, lastJ = J, k++){
	  I = Ilist[K-1-k];
	  J = Jlist[K-1-k];
	  
	  /* outlier lies between sites lastI .. I on reference and lastJ .. J on query */
	  double qrystartpos = X[lastJ] * Xscale; /* QryStartPos */
	  double qryendpos   = X[J] * Xscale; /* QryEndPos */
	  double refstartpos = Y[lastI] * Yscale; /* RefStartPos */
	  double refendpos   = Y[I] * Yscale; /* RefEndPos */
	  double Outlier = fabs(refendpos - refstartpos - qryendpos + qrystartpos);// Indel size
	  int OutlierM = J-lastJ + I-lastI - 2;// misaligned labels

	  int outlier = (iscore[K-k] > outscore[K-k] + (RFLOAT)0.01 || (myoutlierExtend && (lastJ >= J || lastI >= I)) || Outlier > PairMergeMaxOutlier || OutlierM > PairMergeMaxOutlierM) ? 1 : 0;
	  if(DEBUG/* HERE HERE >=2 */) assert(iscore[K-k] >= outscore[K-k] - 0.001);
	  if(!outlier)/* not an outlier */
	    continue;

	  if(PairMergeExcludeOutliers)
	    overlapExclude += (refendpos - refstartpos + qryendpos - qrystartpos) * 0.5;

	  if(PairMergeMaxOutlier < 0)
	    continue;
	  
	  if(PairMergeOutlierChiSq > 0.0 && CmapChimQuality >= 3){// skip if outlier variance is too large or has too many molecule outliers
	    if(refendpos - refstartpos < qryendpos - qrystartpos){
	      Cmap *Ymap = YYmap[Yid];
	      float *OutlierFrac = Ymap->OutlierFrac[0];
	      float *FragSd = Ymap->FragSd[0];
	      float *ExpSd = Ymap->ExpSd[0];
	      float *FragCov = Ymap->FragCov[0];
	      float *FragChiSq = Ymap->FragChiSq[0];
	      if(DEBUG) assert(OutlierFrac != NULL && FragChiSq != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL);

	      float *ChimQuality = Ymap->ChimQuality[0];
	      float *ChimNorm = Ymap->ChimNorm[0];
	  
	      double LRange = (X[J] - X[lastJ]) * Xscale;
	      double LRangeRatio = LRange / AvgInterval;

	      int i = lastI;
	      for(; i < I; i++){
		double SRange = (Y[i+1] - Y[i]) * Yscale;
		double FragCovi = FragCov[i];
		if(PairMergeChimNorm > 0.0 && ChimNorm[i] > 0.0 && ChimNorm[i+1] > 0.0)
		  FragCovi = min(FragCovi, min(ChimNorm[i]*ChimQuality[i],ChimNorm[i+1]*ChimQuality[i+1]) * PairMergeChimNorm * 0.01);

		if(Outlier < PairMergeOutlierMaxSize &&
		   (OutlierFrac[i] > PairMergeOutlierFrac || 
		    (FragCovi < PairMergeOutlierMaxCov && 
		     LRange > PairMergeOutlierMinInterval &&
		     ((FragChiSq[i] < PairMergeOutlierChiSq &&
		      FragSd[i] > ExpSd[i]+ PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
		      FragSd[i] > ExpSd[i] * PairMergeOutlierSdRatio) ||		     
		      FragCovi < PairMergeOutlierMinCov + 
		      ((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != J-lastJ) ? min(PairMergeOutlierMinCovRatio,LRangeRatio) : LRangeRatio) *  PairMergeOutlierMinCovScale))))
		  break;
	      }
	      if(i < I)
		continue;// outlier will be ignored
	    } else {
	      Cmap *Xmap = XXmap[Xid];
	      float *OutlierFrac = Xmap->OutlierFrac[0];
	      float *FragSd = Xmap->FragSd[0];
	      float *ExpSd = Xmap->ExpSd[0];
	      float *FragCov = Xmap->FragCov[0];
	      float *FragChiSq = Xmap->FragChiSq[0];
	      if(DEBUG) assert(OutlierFrac != NULL && FragChiSq != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL);

	      float *ChimQuality = Xmap->ChimQuality[0];
	      float *ChimNorm = Xmap->ChimNorm[0];

	      double LRange = (Y[I] - Y[lastI]) * Yscale;
	      double LRangeRatio = LRange / AvgInterval;

	      int jmin = orientation ? M+1-J : lastJ;
	      int jmax = orientation ? M+1-lastJ : J;
	      int j = jmin;
	      for(; j < jmax; j++){
		double SRange = (X[j+1] - X[j]) * Xscale; // WAS (orientation ? X[M+1-j] - X[M-j] : X[j+1] - X[j]) * Xscale;
		double FragCovj = FragCov[j];
		if(PairMergeChimNorm > 0.0 && ChimNorm[j] > 0.0 && ChimNorm[j+1] > 0.0)
		  FragCovj = min(FragCovj, min(ChimNorm[j]*ChimQuality[j],ChimNorm[j+1]*ChimQuality[j+1]) * PairMergeChimNorm * 0.01);

		if(Outlier < PairMergeOutlierMaxSize &&
		   (OutlierFrac[j] > PairMergeOutlierFrac || 
		    (FragCovj < PairMergeOutlierMaxCov && 
		     LRange > PairMergeOutlierMinInterval &&
		     ((FragChiSq[j] < PairMergeOutlierChiSq &&
		       FragSd[j] > ExpSd[j] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
		       FragSd[j] > ExpSd[j] * PairMergeOutlierSdRatio) ||		     
		      FragCovj < PairMergeOutlierMinCov + 
		      ((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != J-lastJ)?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale))))
		break;
	      }
	      if(j < jmax)
		continue;// outlier will be ignored
	    }
	  }

	  double YLeft = (Y[lastI+1]-Y[lastI]) * Yscale;
	  double YRight = (Y[I] - Y[I-1]) * Yscale;
	  double XLeft = (X[lastJ+1] - X[lastJ]) * Xscale;
	  double XRight = (X[J] - X[J-1]) * Xscale;

	  if(I > lastI + 2){/* lower bounds for inversion/substitution size */
	    double Middle = (Y[I-1] - Y[lastI+1]) * Yscale;
	    double MinOutlier = Middle + max(0.0, max(YLeft-XLeft,YRight-XRight));
	    Outlier = max(Outlier, MinOutlier);
	  }
	  if(J > lastJ + 2){/* lower bounds for inversion/substitution size */
	    double Middle = (X[J-1] - X[lastJ+1]) * Xscale;
	    double MinOutlier = Middle + max(0.0, max(XLeft-YLeft,XRight-YRight));
	    Outlier = max(Outlier, MinOutlier);
	  }
	  
	  if(PairMergeMaxOutlier >= 0.0 && PairMergeMaxEndExpandKB > 0.0 && PairMergeMaxEndExpandKB < PAIRMERGE_HUGE && (Outlier > PairMergeMaxOutlier || OutlierM > PairMergeMaxOutlierM)){/* try to include Outlier in Left or Right endoutlier */
	    if(min(lastI - LI, lastJ - LJ) <= min(RI - I, RJ - J)){/* try to include outlier in Left endoutlier */
	      if(min(Y[I],X[J]) <= PairMergeMaxEndExpandKB || (Y[I] < X[J] ? I : J) - 1 < PairMergeMaxEndN){/* "move" start of alignment to (I,J) */
		LI = I;
		LJ = J;
		overlap = (Y[RI]-Y[LI] + X[RJ]-X[LJ]) * 0.5;
		leftEnd = min(Y[LI], X[LJ]);
		leftEndN = (Y[LI] < X[LJ]) ? LI - 1 : LJ - 1;
		continue;
	      }
	    } else {/* try to include outlier in Right endoutlier */
	      if(min(Y[N+1]-Y[lastI], X[M+1] - X[lastJ]) <= PairMergeMaxEndExpandKB || (Y[N+1]-Y[RI] < X[M+1]-X[RJ] ? N-RI : M-RJ) < PairMergeMaxEndN){/* terminate alignment at lastI,lastJ */
		RI = lastI;
		RJ = lastJ;
		overlap = (Y[RI]-Y[LI] + X[RJ]-X[LJ]) * 0.5;
		rightEnd = min(Y[N+1]-Y[RI], X[M+1] - X[RJ]);
		rightEndN = (Y[N+1]-Y[RI] < X[M+1] - X[RJ]) ? N-RI : M-RJ; 
		break;
	      }
	    }
          }


	  maxOutlier = max(maxOutlier, Outlier);
	  maxOutlierM = max(maxOutlierM, OutlierM);
	}
      }

      if(PairMergeExcludeOutliers)
	overlap -= overlapExclude;

      double MinOverlap = PairMerge;
      double OverlappedMargin = PairMergeOverlappedMargin;
      int LeftSmallerY = (Y[LI] < X[LJ] + OverlappedMargin ) ? 1 : 0;
      int LeftSmallerX = (Y[LI] > X[LJ] - OverlappedMargin ) ? 1 : 0;
      int RightSmallerY = (Y[N+1] - Y[RI] < (X[M+1] - X[RJ]) + OverlappedMargin) ? 1 : 0;
      int RightSmallerX = (Y[N+1] - Y[RI] > (X[M+1] - X[RJ]) - OverlappedMargin) ? 1 : 0;
      int Overlapped = ((LeftSmallerY && RightSmallerY) || (LeftSmallerX && RightSmallerX)) ? 1 : 0;

      if(PairMergeOverlapped >= 0.0 && Overlapped)
	MinOverlap = PairMergeOverlapped;

      int merge = (overlap > MinOverlap && (PAIRMERGE_FIX<=1 || SplitSegDup > 0.0 || leftEndN < PairMergeMaxEndN || leftEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB)) &&
		   (PAIRMERGE_FIX<=1 || SplitSegDup > 0.0 || rightEndN < PairMergeMaxEndN || rightEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB))) ? 1 : 0;/* pairmerge candidate found ? */

      // NOTE : if SplitSegDup > 0.0, all regular merges have already completed so the endoutlier restriction is suppressed

      if(merge && !(PairMergeHmap > 0) && PairMergeMaxOutlier >= 0.0 && (maxOutlier > PairMergeMaxOutlier || maxOutlierM > PairMergeMaxOutlierM))
	merge = 0;/* cannot merge because internal outlier is too large */

      if(!merge)
	Merge = 0;

      if(EVERB){
        printf("pairalignXY:Yid=%d(%lld),Xid=%d(%lld),or=%d,sc=%d:bestscore= %0.6f:Merge=%d:overlap=%0.3f(excl=%0.3f),MinOverlap=%0.3f,maxOutlier=%0.3f,%d,EndOutliers=%0.3f,%0.3f:align:or=%d,sc=%d,score=%0.6f,logPV=%0.2f,Nrange=%d,LI=%d,LJ=%d,RI=%d,RJ=%d,N=%d,M=%d\n",
	       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestscore,Merge,overlap,overlapExclude,MinOverlap,maxOutlier,maxOutlierM,leftEnd,rightEnd,
	       align->orientation,align->scaleID,align->score,align->logPV,align->Nrange,LI,LJ,RI,RJ,N,M);
	fflush(stdout);
      }
      if(DEBUG>=2){
        assert(overlap >= 0.0);
	assert(RI >= LI);
	assert(RJ >= LJ);
      }

      if(!Merge && bestscore > align->score && align->Nrange)
	bestscore = MINSCORE;
    }

    Calign *origalign = NULL;
    if(PairSplit && Psplit > 0.0 && align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= AlignedSiteThreshold){
      origalign = new Calign[1];
      copy(origalign, align, 1, 1);
    }

    if(bestscore > origscore){
      if(Palindromic_Delta > 0.0 && PairSplit && align->score > ScoreThreshold && align->orientation != orientation){/* save alignment from opposite orientation as 2nd best alignment */
	if(DEBUG) assert(orientation == 1);
	delete [] align->align2; align->align2 = NULL;
	Calign *align2 = align->align2 = new Calign[1];
	copy(align2, align, 1, 1);
	align2->score = align->score;
	align2->orientation = align->orientation;
	align2->scaleID = align->scaleID;
	align2->logPV = align->logPV;
	//	align->I2 = align->sites1[alignK-1];
	//	align->K2 = align->sitesK1[alignK-1];
	//	align->J2 = align->sites2[alignK-1];
	align2->Lend = align->Lend;
	align2->Rend = align->Rend;
	//	align->numpairs2 = alignK;
      }

      align->score = bestscore;
      align->logPV = bestLogPV;
      align->orientation = orientation;
      align->scaleID = scaleID;
      align->Lend = I;
      align->Rend = bestR;
      if(PAIRMERGE_FIX && PairMerge)
	align->Nrange = Merge;

      I = Ilist[K-1];// align->sites1[0];
      J = Jlist[K-1];// align->sites2[0];

      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);

      align->Lij1 = (FLAT_A(A.G, I, J) <= -2) ? I : FLAT_A(A.Lijy, I, J);
      align->Lij2 = (FLAT_A(A.G, I, J) <= -2) ? J : FLAT_A(A.Lijx, I, J);
      align->Rij1 = (bestR <= -2) ? bestI : FLAT_A(A.Rijy, bestI, bestJ);
      align->Rij2 = (bestR <= -2) ? bestJ : FLAT_A(A.Rijx, bestI, bestJ);

      align->repeat = -1;
      if(RepeatRec && RepeatMaxShift > 0 && RepeatLogPvalueRatio > 0.0)
	align->repeat = 0;

      /* (re)allocate best alignment */
      align->allfree();/* resets align->numpairs = 0 (if needed : typically not needed) */
      if(DEBUG) assert(align->score > MINSCORE);
      if(PairSplit || EVERB || (align->score > ScoreThreshold && align->logPV > LogPvThreshold && K >= AlignedSiteThreshold)){/* save alignment */
	if(DEBUG) assert(K > 0);
	align->numpairs = K;
	align->expand_arrays(K);

	if(DEBUG && !(min(M,align->Rij2) >= max(1,align->Lij2) + align->numpairs - 1)){
	  #pragma omp critical
	  {
	    printf("pairalignXY:Yid=%d(%lld),Xid=%d(%lld),or=%d,sc=%d:bestscore= %0.6f,bestI=%d,bestJ=%d,bestR=%d, Rijy=%d,Rijx=%d; I=%d,J=%d,G=%d,Lijy=%d,Lijx=%d,K=%d\n",
		   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestscore,bestI,bestJ,bestR,
		   FLAT_A(A.Rijy,bestI,bestJ),FLAT_A(A.Rijx,bestI,bestJ),I,J,FLAT_A(A.G,I,J),FLAT_A(A.Lijy,I,J),FLAT_A(A.Lijx,I,J),K);
	    
	    printf("M=%d,align:Rij2=%d,Lij2=%d,numpairs=%d\n",M,align->Rij2,align->Lij2,align->numpairs);
	    fflush(stdout);
	    assert(min(M,align->Rij2) >= max(1,align->Lij2) + align->numpairs - 1);
	  }
	}
	if(DEBUG) assert(min(N,align->Rij1) >= max(1,align->Lij1) + align->numpairs - 1);

	/* copy alignment from Ilist[],Jlist[] in reverse order */
        #pragma ivdep
	for(int I= 0; I < K; I++){
	  align->sites1[I] = Ilist[K-1-I];
	  align->sites2[I] = Jlist[K-1-I];
	  align->iscore[I] = iscore[K-I];
	  align->outscore[I] = outscore[K-I];
	  if(DEBUG>=2) assert(isfinite(align->outscore[I]));
	  if(DEBUG>=2) assert(align->iscore[I] >= align->outscore[I] - 0.001);
	}
	align->iscore[K] = iscore[0];
	align->outscore[K] = outscore[0];
	if(DEBUG>=2) assert(isfinite(align->outscore[K]));

	/* check if bestscore is a repeat */
	if(RepeatRec && RepeatMaxShift > 0 && RepeatPvalueRatio > 0.0){
	  // if(DEBUG) assert(!PairSplit);
	  if(VERB && (VERB>=2 || EVERB)){
            #pragma omp critical
	    {
	      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:M=%d,N=%d:checking for repeats\n",
		     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,M,N);
	      fflush(stdout);
	    }
	  }

	  double minPvalue = LogPvThreshold * RepeatLogPvalueRatio /* WAS89 LogPvThreshold + log(RepeatPvalueRatio)/log(10.0)*/;

	  if(DEBUG>=2) assert(IMIN <= bestI && bestI <= IMAX && JMIN[bestI] <= bestJ && bestJ <= JMAX[bestI]);

	  int bestI2 = bestI;
	  int bestJ2 = bestJ;

	  for(bestJ2 = max(1, align->sites2[max(0,K-1-RepeatMaxShift)]); !align->repeat && bestJ2 < bestJ; bestJ2++){/* back track best alignment from (bestI2,bestJ2) and compute its Pvalue */	
	  // WAS	  for(bestJ2 = max(1,bestJ - RepeatMaxShift); !align->repeat && bestJ2 < bestJ; bestJ2++){/* back track best alignment from (bestI2,bestJ2) and compute its Pvalue */
	    if(!(JMIN[bestI2] <= bestJ2 && bestJ2 <= JMAX[bestI2]))
	      continue;

	    double MinShiftKb = (X[bestJ]-X[bestJ2]) * RepeatKbShiftRatio; /* minimum offset shift along entire alignment */
	    double MaxShiftKb = (X[bestJ]-X[bestJ2]) * RepeatKbShiftRatio2;/* maximum offset shift along entier alignment (to avoid creating new outliers) */
	    int I = bestI2;
	    int J = bestJ2;

	    /* First compute bestscore2, bestR2 with rightmost alignment at (I,J) */
	    int bestR2 = localtype;
	    PFLOAT bestscore2 = C;
	    PFLOAT Ascore = FLAT_A(A.score, I, J);/* NOTE: Ascore will be added to bestscore2 last, unlike for bestscore at bestI,bestJ */

	    int *RijxI = &FLAT_A(A.Rijx, I, 0), *RijyI = &FLAT_A(A.Rijy, I, 0);
	    PFLOAT Ydiff = Yr[I];// Y[N+1]-Y[I];

	    if(EVERB>=2 && I==33 && J==46){
	      printf("Yid=%d(id=%lld),Xid=%d(id=%lld(),bestI=%d,bestJ=%d,bestI2=%d,bestJ2=%d:I=%d,J=%d:original Rijx[I][J]=%d,Rijy[I][J]=%d,Send(%0.8f,%d+%d)=%0.6f\n",
		     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,bestI,bestJ,bestI2,bestJ2,I,J, RijxI[J], RijyI[J],
		     PMIN(Xr[J],Ydiff),min(M,RijxI[J])+1-J, min(N,RijyI[J])+1-I, Send(PMIN(Xr[J],Ydiff),(PFLOAT)(min(M,RijxI[J])+1-J + min(N,RijyI[J])+1-I)));
	      fflush(stdout);
	    }

	    if(ENDINIT>=1 && I <= N-DELTA && J <= M-DELTA){/* must initialize Rijx,Rijy since this case was previously not initialized */
	      int Rijy = N+1;
	      int Rijx = M+1;

	      /* update Rijy by decrementing it until Y[Rijy] overlaps X */
	      FLOAT delta = Y[I] + X[M+1] - X[J];
	      while(Y[Rijy] > delta)
		Rijy--;
	      if(DEBUG>=2) assert(Rijy >= I);
	      RijyI[J] = Rijy;
	    
	      /* update Rijx by decrementing it until X[Rijx] overlaps Y */
	      delta = X[J] + Y[N+1] - Y[I];
	      while(X[Rijx] > delta)
		Rijx--;
	      if(DEBUG>=2)	assert(Rijx >= J);
	      RijxI[J] = Rijx;

	      if(EVERB>=2 && I==33 && J==46){
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),bestI=%d,bestJ=%d,bestI2=%d,bestJ2=%d:I=%d,J=%d:Updated Rijx[I][J]=%d,Rijy[I][J]=%d\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,bestI,bestJ,bestI2,bestJ2,I,J, RijxI[J], RijyI[J]);
		fflush(stdout);
	      }
	    }

	    PFLOAT score = Send(PMIN(Xr[J], Ydiff), (PFLOAT)(min(M, RijxI[J])+1-J + min(N, RijyI[J])+1-I));
	    if(EVERB>=2 && I==33 && J==46){
	      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),bestI=%d,bestJ=%d,bestI2=%d,bestJ2=%d:I=%d,J=%d:Initial Rijx[I][J]=%d,Rijy[I][J]=%d,Send(%0.8f,%d+%d)=%0.6f\n",
		     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,bestI,bestJ,bestI2,bestJ2,I,J, RijxI[J], RijyI[J],
		     PMIN(Xr[J],Ydiff),min(M,RijxI[J])+1-J, min(N,RijyI[J])+1-I, Send(PMIN(Xr[J],Ydiff),(PFLOAT)(min(M,RijxI[J])+1-J + min(N,RijyI[J])+1-I)));
	      fflush(stdout);
	    }

	    if(ENDFIX){/* update Rijx,Rijy based on Sbnd() as usual */
	      if(RijxI[J] > M) {
		double a = M+1-I-J;
		int K;
		for(K= RijyI[J]; K > I; K--){
		  PFLOAT bound = Sbnd(Xr[J],Y[K]-Y[I], K+a);
		  if(bound <= score)
		    break;// HERE : see ENDFIX2 in refalign.cpp
		  if(EVERB>=2 && I==33 && J==46){
		    printf("    I=%d,J=%d,N=%d,M=%d:K=%d,score=%0.6f,Sbnd(Xr[J]=%0.8f,Y[K]-Y[I]=%0.8f,min(M,Rijx)+1-J=%d,min(N,K-1)+1-I=%d(K+a=%0.1f))=%0.6f:Rijy[I][J]=%d->%d\n",
			   I,J,N,M,K,score,Xr[J],Y[K]-Y[I],min(M,RijxI[J])+1-J,min(N,K-1)+1-I,K+a,Sbnd(Xr[J],Y[K]-Y[I],K+a),K,K-1);
		    fflush(stdout);
		  }
		  score = bound;
		}
		RijyI[J] = K;
	      }
	      if(RijyI[J] > N) {
		double a = N+1-I-J;
		int K;
		for(K= RijxI[J];K > J;K--){
		  PFLOAT bound = Sbnd(X[K]-X[J], Ydiff,K+a);
		  if(bound <= score)
		    break;// HERE : see ENDFIX2 in refalign.cpp
		  score = bound;
		}
		RijxI[J] = K;
	      }
	      if(EVERB>=2 && I==33 && J==46){
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),bestI=%d,bestJ=%d,bestI2=%d,bestJ2=%d:I=%d,J=%d:Final Rijx[I][J]=%d,Rijy[I][J]=%d\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,bestI,bestJ,bestI2,bestJ2,I,J, RijxI[J], RijyI[J]);
		fflush(stdout);
	      }
	    }
	    if(score > bestscore2){
	      bestscore2 = score;
	      bestR2 = -1;
	    }
	    bestscore2 += Ascore;

	    if(VERB && (VERB>=2 || EVERB>=2)){
              #pragma omp critical
	      {
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,M=%d,N=%d:checking for repeats\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,M,N);
		fflush(stdout);
	      }
	    }

	    /* Now backtrack to extract best alignment starting at (I,J)==(bestI2,bestJ2) constrained by requiring J to be less than in the best alignment for the same I */
	    int U = K-1;/* keep (I2 == align->sites1[U]) <= I, and J > (J2 == align->sites2[T2]) to satisfy the repeat constraint
			   Also keep MinShiftKb <= Y[I2] - Y[I] + X[J] - X[J2] <= MaxShiftKb */
	    int K2 = 0;

	    int bestG = -1, bestH = -1;
	    if(DEBUG>=2) assert(I==bestI2 && J==bestJ2);

	    while(I >= 0){
	      Ilist2[K2] = I;
	      Jlist2[K2++] = J;

	      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);

	      /* update A.G[I][J] & A.H[I][J] by repeating recurrance computation constrained by J being less than in the best alignment for the same I */
	      /* NOTE : this needs to be done even if DEFER_BP==0, since original A.G & A.H values would have been computed without the constraint */
	      /* NOTE2 : A.score[I][J] is NOT updated, since we will rely on Pvalue to decide if the overall alignment found is good enough */

	      if(FIX_CONSTRAINT && (FIX_CONSTRAINT>=2 || bestJ-bestJ2 >= 2)){
		while(U > 0 && align->sites1[U] > I)
		  U--;
	      } else {
		while(U > 0 && align->sites1[U-1] >= I)
		  U--;
	      }

	      if(DEBUG>=2 && MinShiftKb > 0.0 && align->sites1[U] >= I && !(J < align->sites2[U])){
                #pragma omp critical 
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,M=%d,N=%d:checking for repeats\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,M,N);
		  printf("  I=%d,J=%d:U=%d,align->sites1[U]=%d,align->sites2[U]=%d\n",I,J,U,align->sites1[U],align->sites2[U]);
		  fflush(stdout);
		  assert(J < align->sites2[U]);
		}
	      }

	      PFLOAT score = FLAT_A(A.Uscore, I, J);
	      if(DEBUG>=2 && !(isfinite(score))){
                #pragma omp critical 
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,M=%d,N=%d:checking for repeats\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,M,N);
		  printf("  I=%d,J=%d:U=%d,align->sites1[U]=%d,align->sites2[U]=%d,A.Uscore[I][J]=%0.8e\n",I,J,U,align->sites1[U],align->sites2[U],score);
		  fflush(stdout);
		  assert(isfinite(score));
		}
	      }
	      bestG = localtype;

	      int Gmin = max(I-DELTA,IMIN);
	      int T3 = U;

	      for(int G= I; --G >= Gmin;){
		int Hmin = max(J-DELTA,JMIN[G]);

		FLOAT XG;/* best estimate of original alignment location on X[] corresponding to Y[G] */
		int origT3 = T3;
		while(T3 < K-1 && align->sites1[T3] < G)
		  T3++;
		while(T3 > 0 && align->sites1[T3-1] >= G)
		  T3--;
		int I3 = align->sites1[T3];
		int J3 = align->sites2[T3];

		if(DEBUG>=2 && !(T3 >= K-1 || I3 >= G)){
		  #pragma omp critical
		  {
		    printf("While Checking for repeats:RepeatRec=%d,MaxShift=%d,LogPvalueRatio=%0.4f,%0.4f,KbShiftRatio=%0.3f,%0.3f,MaxDrop=%d,extend=%d,Lij2=%d,Rij2=%d,M=%d(score=%0.6f,LogPV=%0.2f,numpairs=%d,repeat=%d)\n",
			   RepeatRec,RepeatMaxShift,RepeatLogPvalueRatio,RepeatLogPvalueRatioPerLabel,RepeatKbShiftRatio,RepeatKbShiftRatio2,RepeatMaxLabelDrop,
			   extend,align->Lij2,align->Rij2,M,align->score,align->logPV,align->numpairs,align->repeat);
		    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,maptype=%d:\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation, scaleID, bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,maptype);
		    printf("K=%d,U=%d,I=%d,J=%d,G=%d:T3=%d->%d:I3=%d,J3=%d\n",
			   K,U,I,J,G,origT3, T3,I3,J3);
		    fflush(stdout);
		    assert(T3 >= K-1 || I3 >= G);
		  }
		}
			
		int T4 = T3;
		while(T4 > 0 && align->sites1[T4] > G)
		  T4--;
		int I4 = align->sites1[T4];
		
		if(DEBUG>=2 && !(T4 <= 0 || I4 <= G)){
		  #pragma omp critical
		  {
		    printf("While Checking for repeats:RepeatRec=%d,MaxShift=%d,LogPvalueRatio=%0.4f,%0.4f,KbShiftRatio=%0.3f,%0.3f,MaxDrop=%d,extend=%d,Lij2=%d,Rij2=%d,M=%d(score=%0.6f,LogPV=%0.2f,numpairs=%d,repeat=%d)\n",
			   RepeatRec,RepeatMaxShift,RepeatLogPvalueRatio,RepeatLogPvalueRatioPerLabel,RepeatKbShiftRatio,RepeatKbShiftRatio2,RepeatMaxLabelDrop,
			   extend,align->Lij2,align->Rij2,M,align->score,align->logPV,align->numpairs,align->repeat);
		    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,maptype=%d:\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID, bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,maptype);
		    printf("K=%d,U=%d,I=%d,J=%d,G=%d:T3=%d->%d,T4=%d:I3=%d,J3=%d,I4=%d\n",
			   K,U,I,J,G,origT3, T3, T4, I3,J3, I4);
		    fflush(stdout);
		    assert(T4 <= 0 || I4 <= G);
		  }
		}
		origT3 = T3;

		if(G < I4){
		  int J4 = align->sites2[T4];
		  if(I4  <= G)
		    XG = X[J4];
		  else
		    XG = X[J4] - (Y[I4] - Y[G]);
		} else if(I3 > I4){
		  if(DEBUG>=2) assert(I4 <= G && G <= I3);
		  if(DEBUG>=2) assert(Y[I3] > Y[I4]);
		  FLOAT Gfrac = (Y[G] - Y[I4])/(Y[I3] - Y[I4]);
		  int J4 = align->sites2[T4];
		  XG = X[J4] + Gfrac * (X[J3] - X[J4]);
		} else  {
		  if(DEBUG>=2) assert(I3 == I4);
		  if(DEBUG>=2) assert(T3 == T4);
		  XG = X[J3];
		}

		if(FIX_CONSTRAINT && (FIX_CONSTRAINT>=2 || bestJ-bestJ2 >= 2)){
		  while(T3 > 0 && align->sites1[T3] > G)
		    T3--;
		} else {
		  while(T3 > 0 && align->sites1[T3-1] >= G)
		    T3--;
		}
		I3 = align->sites1[T3];
		J3 = align->sites2[T3];

		/* H must satistfy the constraint: H < J3 AND -MaxShiftKb <= Y[I3]-Y[G] + X[H] - XG <= -MinShiftKb */
		FLOAT XHmax = XG - Y[I3];
		FLOAT XHmin = XHmax;
		int Hmax = min(J,J3);/* H must be less than this value */
		if(I3 < G && T3 < K-1){
		  int I3 = align->sites1[T3+1];
		  int J3 = align->sites2[T3+1];
		  if(X[J3] - Y[I3] > XHmax)
		    XHmax = X[J3] - Y[I3];
		  if(X[J3] - Y[I3] < XHmin)
		    XHmin = X[J3] - Y[I3];
		  if(Hmax < J3)
		    Hmax = min(J,J3);
		}
		if(G < I3 && T3 > 0){
		  int I3 = align->sites1[T3-1];
		  int J3 = align->sites2[T3-1];
		  if(X[J3] - Y[I3] > XHmax)
		    XHmax = X[J3] - Y[I3];
		  if(X[J3] - Y[I3] < XHmin)
		    XHmin = X[J3] - Y[I3];
		  if(Hmax < J3)
		    Hmax = min(J,J3);
		}
		XHmax += Y[G] - MinShiftKb;/* X[H] must be <= this value */
		XHmin += Y[G] - MaxShiftKb;/* X[H] must be >= this value */

		PFLOAT deltaY = delY[I][G-I+DELTA];// Y[I] - Y[G] 
		if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
		PFLOAT *AscoreG = &FLAT_A(A.score, G, 0);
		// WAS107		PFLOAT a = J+I-G;
	
		Hmax = min(Hmax,JMAX[G]+1);

		for(int H= Hmin; H < Hmax; H++){
		  if(X[H] < XHmin && H < Hmax)
		    continue;
		  if(X[H] > XHmax)
		    break;

		  if(DEBUG>=2) assert(IMIN <= G && G <= IMAX && JMIN[G] <= H && H <= JMAX[G]);

   	          if(DEBUG>=2 && !(1 <= H && H < J && J-H <= DELTA && J <= M)){
	            #pragma omp critical
	            {
	              printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d,M=%d\n",Yid,Xid,orientation,scaleID,I,J,G,H,M);
		      fflush(stdout);
		      assert(1 <= H && H < J && J-H <= DELTA && J <= M);
	            }
	          }
		  PFLOAT deltaX = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];// X[J] - X[H]
		  if(DEBUG>=2) assert(fabs(X[J] - X[H] - deltaX) < 0.001);
		  PFLOAT newscore = AscoreG[H] +  Sint(deltaX,deltaY,J-H,I-G,myoutlierBC);// WAS107 Sint(deltaX,deltaY,(PFLOAT)(a-H));
		  if(newscore > score){
		    if(EVERB>=2)
		      printf("I=%d,J=%d,G=%d,H=%d:T3=%d,I3=%d,J3=%d:score=%0.6f->%0.6f\n",
			     I,J,G,H,T3,align->sites1[T3],align->sites2[T3],score,newscore);
		    bestG = G;
		    bestH = H;
		    score = newscore;
		  }
		}
	      }// for(int G= I; --G >= Gmin;)

	      if(bestG < 0){/* Left End : value of bestH is undefined */
		if(DEFER_G >= 1 && score > C){
		  bestG = -1;
		  if(DEFER_G >= 2){/* check if actual value of Lijx,Lijy needs to be incremented */
		    int Lijx = FLAT_A(A.Lijx,I,J);
		    int Lijy = FLAT_A(A.Lijy,I,J);
		    PFLOAT SmIJ = Sm(J,I);
		    PFLOAT ascore = C;

		    if(Lijx <= 0){
		      ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I+1-max(1,Lijy))) + SmIJ;
		      ascore = max(ascore,C);
		      for(int K= Lijy; K < I; ) {
			PFLOAT bound = Sbnd(X[J],Y[I]-Y[K],(PFLOAT)(J+I-K)) + SmIJ;
			if(bound <= ascore)
			  break;
			ascore = bound;
			Lijy = ++K;
		      }
		      FLAT_A(A.Lijy, I, J) = Lijy;
		    } else {
		      if(DEBUG>=2) assert(Lijy <= 0);
		      ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-max(1,Lijx)+I)) + SmIJ;
		      ascore = max(ascore,C);
		      for(int K= Lijx; K < J; ){
			PFLOAT bound = Sbnd(X[J]-X[K],Y[I],(PFLOAT)(J+I-K)) + SmIJ;
			if(bound <= ascore)
			  break;
			ascore = bound;
			Lijx = ++K;
		      }
		      FLAT_A(A.Lijx, I, J) = Lijx;
		    }
		  }
		}
		FLAT_A(A.G, I, J) = bestG;
		J = -1;
	      } else {
		FLAT_A(A.G, I, J) = bestG;
		FLAT_A(A.H, I, J) = bestH;
		J = bestH;
	      }
	      I = bestG;
	    }
	  
	    if(DEBUG>=2) assert(I < 0 || (IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]));

	    if(DEBUG>=2){
	      if(DEBUG>=2 && bestR2 >= -1) assert(IMIN <= bestI2 && bestI2 <= IMAX && JMIN[bestI2] <= bestJ2 && bestJ2 <= JMAX[bestI2]);

	      int Rijy = (bestR2 <= -2) ? bestI2 : FLAT_A(A.Rijy, bestI2, bestJ2);
	      int Rijx = (bestR2 <= -2) ? bestJ2 : FLAT_A(A.Rijx, bestI2, bestJ2);

	      if(!(min(N,Rijy)-bestI2 >= 0)){
                #pragma omp critical
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,Rijy=%d,Rijx=%d,M=%d,N=%d\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,Rijy,Rijx,M,N);
		  for(int T2 = 0; T2 < K2; T2++)
		    printf("   T=%d/%d:I=%d,J=%d\n",T2,K2,Ilist2[K2-1-T2],Jlist2[K2-1-T2]);
		  fflush(stdout);
		  assert(min(N,Rijy)-bestI2 >= 0);
		}
	      }
	    }

	    double logPV2 = alignFP(A,bestI2,bestJ2,bestR2,Y,X,N,M,orientation,scaleID,Yid,Xid,delY,delX,bestscore2, /* (EVERB && bestJ-bestJ2==2) ? 1 : 0*/ 0, 0, myoutlierExtend, myoutlierBC, tid,
				    JMIN,JMAX,IMIN,IMAX);
	    if(logPV2 >= max(align->logPV * RepeatLogPvalueRatio, minPvalue) &&
	       (RepeatMaxLabelDrop < 0 || K2 >= align->numpairs - abs(bestJ - bestJ2) - RepeatMaxLabelDrop) &&
	       (logPV2 / (K2-1)) >= (align->logPV / (align->numpairs-1)) * RepeatLogPvalueRatioPerLabel ){

	      align->repeat = bestJ - bestJ2;
	      if(DEBUG) assert(align->repeat > 0);

	      if(VERB>=2 || EVERB){
                #pragma omp critical
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.6f,numpairs=%d:bestI2=%d,bestJ2=%d,score2=%0.4f,logPV2=%0.6f (repeat=%d confirmed),K2=%d\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,bestscore2,logPV2,align->repeat,K2);
		  int I2 = -1, J2= -1,G2= -1,H2= -1;
		  for(int T2 = 0; T2 < K2; T2++, G2=I2, H2=J2){
		    I2 = Ilist2[K2-1-T2];
		    J2 = Jlist2[K2-1-T2];
		    if(DEBUG>=2) assert(IMIN <= I2 && I2 <= IMAX && JMIN[I2] <= J2 && J2 <= JMAX[I2]);
		    if(T2 > 0){
		      PFLOAT score = FLAT_A(A.score, I2, J2) - FLAT_A(A.score, G2, H2);
		      PFLOAT y = Y[I2] - Y[G2];
		      PFLOAT x = X[J2] - X[H2];
		      PFLOAT Bias, Pen;
		      SintDetail(x,y, J2-H2, I2-G2, J2, I2, Bias, Pen);
		      printf("   T=%d/%d:I=%d,J=%d (Y=%0.4f X=%0.4f):score=%0.4f(Pen=%0.4f,Bias=%0.4f):A.score=%0.4f\n",
			     T2,K2,I2,J2, y, x, score, Pen, Bias, FLAT_A(A.score,I2,J2));
		    } else { /* left end */
		      PFLOAT score = FLAT_A(A.score, I2, J2);
		      printf("   T2=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:A.score=%0.4f\n",T2,K2,I2, J2, A.Lijy[I2][J2],A.Lijx[I2][J2],score);
		    }
		  }
		  /* right end */
		  int Rijy = (I2 >= N-DELTA || J2 >= M-DELTA) ? A.Rijy[I2][J2] : N+1;
		  int Rijx = (I2 >= N-DELTA || J2 >= M-DELTA) ? A.Rijx[I2][J2] : M+1;
		  if(bestR2 <= -2){
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.4f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.4f,ChimScore=%0.6f,total score=%0.8e\n",
			   I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],ChimScore,FLAT_A(A.score, I2, J2)+ChimScore);
		  } else if(X[M+1]-X[J2] <= Y[N+1]-Y[I2] && Rijy <= N && Y[Rijy+1]-Y[I2] < X[M+1]-X[J2]){
		    if(DEBUG) assert(Rijx > M);
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.8f,%0.8f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],X[M+1]-X[J2],Y[Rijy+1]-Y[I2],min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			   Sbnd(X[M+1]-X[J2],Y[Rijy+1]-Y[I2],M+1-J2,Rijy+1-I2),FLAT_A(A.score, I2, J2) + Sbnd(X[M+1]-X[J2],Y[Rijy+1]-Y[I2],M+1-J2,Rijy+1-I2));
		  } else if(X[M+1]-X[J2] >= Y[N+1]-Y[I2] && Rijx <= M && X[Rijx+1]-X[J2] < Y[N+1]-Y[I2]){
		    if(DEBUG) assert(Rijy > N);
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],
			   X[Rijx+1]-X[J2],Y[N+1]-Y[I2],min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			   Sbnd(X[Rijx+1]-X[J2],Y[N+1]-Y[I2],Rijx+1-J2,N+1-I2),
			   FLAT_A(A.score, I2, J2) + Sbnd(X[Rijx+1]-X[J2],Y[N+1]-Y[I2],Rijx+1-J2,N+1-I2));
		  } else
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Send(%0.8f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],
			   min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			   Send(min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2),
			   FLAT_A(A.score, I2, J2) + Send(min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2));

		  printf("  Best alignment:\n");
		  int I = -1, J = -1, G = -1, H =  -1;
		  for(int T = 0; T < K; T++, G=I, H=J){
		    I = Ilist[K-1-T];
		    J = Jlist[K-1-T];
		    if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
		    if(T > 0){/* internal interval */
		      PFLOAT score = FLAT_A(A.score, I, J) - FLAT_A(A.score, G, H);
		      PFLOAT y = Y[I] - Y[G];
		      PFLOAT x = X[J] - X[H];
		      PFLOAT Bias, Pen;
		      SintDetail(x, y, J-H, I-G, J, I, Bias, Pen);
		      printf("   T=%d/%d:I=%d,J=%d (Y=%0.4f X=%0.4f):score=%0.4f(Pen=%0.4f,Bias=%0.4f):A.score=%0.4f\n",
			     T,K,I,J, y, x, score, Pen, Bias, FLAT_A(A.score, I, J));
		    } else { /* left end */
		      PFLOAT score = FLAT_A(A.score, I2, J2);
		      printf("   T=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:A.score=%0.4f\n", T,K,I, J, A.Lijy[I][J],A.Lijx[I][J],score);
		    }
		  }
		  /* right end */
		  Rijy = (I >= N-DELTA || J >= M-DELTA) ? A.Rijy[I][J] : N+1;
		  Rijx = (I >= N-DELTA || J >= M-DELTA) ? A.Rijx[I][J] : M+1;

		  if(bestR <= -2){
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.4f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.4f,ChimScore=%0.6f,total score=%0.8e\n",
			   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],ChimScore,FLAT_A(A.score, I, J)+ChimScore);
		  } else if(X[M+1]-X[J] <= Y[N+1]-Y[I] && Rijy <= N && Y[Rijy+1]-Y[I] < X[M+1]-X[J]){
		    if(DEBUG) assert(Rijx > M);
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],X[M+1]-X[J],Y[Rijy+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			   Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I),FLAT_A(A.score, I, J) + Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I));
		  } else if(X[M+1]-X[J] >= Y[N+1]-Y[I] && Rijx <= M && X[Rijx+1]-X[J] < Y[N+1]-Y[I]){
		    if(DEBUG) assert(Rijy > N);
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
			   X[Rijx+1]-X[J],Y[N+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			   Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I),
			   FLAT_A(A.score, I, J) + Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I));
		  } else
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Send(%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
			   min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			   Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I),
			   FLAT_A(A.score, I, J) + Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I));
		
		  fflush(stdout);
		}
	      }
	      if(DEBUG) assert(align->repeat > 0);
	      break;
	    } else if(VERB>=2 || (EVERB>=2 && bestJ-bestJ2==2)){
              #pragma omp critical
	      {
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.6f,numpairs=%d:bestI2=%d,bestJ2=%d,score2=%0.4f,logPV2=%0.6f (No repeat at shiftJ =%d)\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,bestscore2,logPV2,bestJ-bestJ2);
		int I2 = -1, J2= -1,G2= -1,H2= -1;
		for(int T2 = 0; T2 < K2; T2++, G2=I2, H2=J2){
		  I2 = Ilist2[K2-1-T2];
		  J2 = Jlist2[K2-1-T2];
		  if(DEBUG>=2) assert(IMIN <= I2 && I2 <= IMAX && JMIN[I2] <= J2 && J2 <= JMAX[I2]);
		  if(T2 > 0){
		    PFLOAT score = FLAT_A(A.score, I2, J2) -FLAT_A(A.score, G2, H2);
		    PFLOAT y = Y[I2] - Y[G2];
		    PFLOAT x = X[J2] - X[H2];
		    PFLOAT Bias, Pen;
		    SintDetail(x,y, J2-H2, I2-G2, J2, I2, Bias, Pen);
		    printf("   T=%d/%d:I=%d,J=%d (Y=%0.4f X=%0.4f):score=%0.4f(Pen=%0.4f,Bias=%0.4f):A.score=%0.4f\n",
			   T2,K2,I2,J2, y, x, score, Pen, Bias, FLAT_A(A.score,I2,J2));
		  } else { /* left end */
		    PFLOAT score = FLAT_A(A.score, I2, J2);
		    printf("   T2=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:A.score=%0.4f\n",T2,K2,I2, J2, A.Lijy[I2][J2],A.Lijx[I2][J2],score);
		  }
		}
		/* right end */
		int Rijy = (I2 >= N-DELTA || J2 >= M-DELTA) ? A.Rijy[I2][J2] : N+1;
		int Rijx = (I2 >= N-DELTA || J2 >= M-DELTA) ? A.Rijx[I2][J2] : M+1;
		if(bestR2 <= -2){
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.4f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.4f,ChimScore=%0.6f,total score=%0.8e\n",
			 I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],ChimScore,FLAT_A(A.score, I2, J2)+ChimScore);
		} else if(X[M+1]-X[J2] <= Y[N+1]-Y[I2] && Rijy <= N && Y[Rijy+1]-Y[I2] < X[M+1]-X[J2]){
		  if(DEBUG) assert(Rijx > M);
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],X[M+1]-X[J2],Y[Rijy+1]-Y[I2],min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			 Sbnd(X[M+1]-X[J2],Y[Rijy+1]-Y[I2],M+1-J2,Rijy+1-I2),FLAT_A(A.score, I2, J2) + Sbnd(X[M+1]-X[J2],Y[Rijy+1]-Y[I2],M+1-J2,Rijy+1-I2));
		} else if(X[M+1]-X[J2] >= Y[N+1]-Y[I2] && Rijx <= M && X[Rijx+1]-X[J2] < Y[N+1]-Y[I2]){
		  if(DEBUG) assert(Rijy > N);
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],
			 X[Rijx+1]-X[J2],Y[N+1]-Y[I2],min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			 Sbnd(X[Rijx+1]-X[J2],Y[N+1]-Y[I2],Rijx+1-J2,N+1-I2),
			 FLAT_A(A.score, I2, J2) + Sbnd(X[Rijx+1]-X[J2],Y[N+1]-Y[I2],Rijx+1-J2,N+1-I2));
		} else
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Send(%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],
			 min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			 Send(min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2),
			 FLAT_A(A.score, I2, J2) + Send(min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2));

		printf("  Best alignment:\n");
		int I = -1, J = -1, G = -1, H =  -1;
		for(int T = 0; T < K; T++, G=I, H=J){
		  I = Ilist[K-1-T];
		  J = Jlist[K-1-T];
		  if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
		  if(T > 0){/* internal interval */
		    PFLOAT score = FLAT_A(A.score, I, J) - FLAT_A(A.score, G, H);
		    PFLOAT y = Y[I] - Y[G];
		    PFLOAT x = X[J] - X[H];
		    PFLOAT Bias, Pen;
		    SintDetail(x, y, J-H, I-G, J, I, Bias, Pen);
		    printf("   T=%d/%d:I=%d,J=%d (Y=%0.4f X=%0.4f):score=%0.4f(Pen=%0.4f,Bias=%0.4f):A.score=%0.4f\n",
			   T,K,I,J, y, x, score, Pen, Bias, FLAT_A(A.score, I, J));
		  } else { /* left end */
		    PFLOAT score = FLAT_A(A.score, I2, J2);
		    printf("   T=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:A.score=%0.4f\n", T,K,I, J, A.Lijy[I][J],A.Lijx[I][J],score);
		  }
		}
		/* right end */
		Rijy = (I >= N-DELTA || J >= M-DELTA) ? A.Rijy[I][J] : N+1;
		Rijx = (I >= N-DELTA || J >= M-DELTA) ? A.Rijx[I][J] : M+1;

		if(bestR <= -2){
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.4f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.4f,ChimScore=%0.6f,total score=%0.8e\n",
			 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],ChimScore,FLAT_A(A.score, I, J)+ChimScore);
		} else if(X[M+1]-X[J] <= Y[N+1]-Y[I] && Rijy <= N && Y[Rijy+1]-Y[I] < X[M+1]-X[J]){
		  if(DEBUG) assert(Rijx > M);
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],X[M+1]-X[J],Y[Rijy+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			 Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I),FLAT_A(A.score, I, J) + Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I));
		} else if(X[M+1]-X[J] >= Y[N+1]-Y[I] && Rijx <= M && X[Rijx+1]-X[J] < Y[N+1]-Y[I]){
		  if(DEBUG) assert(Rijy > N);
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
			 X[Rijx+1]-X[J],Y[N+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			 Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I),
			 FLAT_A(A.score, I, J) + Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I));
		} else
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Send(%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
			 min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			 Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I),
			 FLAT_A(A.score, I, J) + Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I));
		
		fflush(stdout);
	      }
	    }
	  } // 	for(bestJ2 = max(1,bestJ - RepeatMaxShift); bestJ2 < bestJ; bestJ2++)

	  if(DEBUG && !align->repeat && !(bestJ2 == bestJ)){
	    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:M=%d,N=%d:while checking for repeats\n",
		   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,M,N);
	    printf("bestJ2=%d,bestJ=%d,RepeatMaxShift=%d\n",bestJ2,bestJ,RepeatMaxShift);
	    fflush(stdout);
	    assert(bestJ2 == bestJ);
	  }

	  for(bestI2 = max(1, align->sites1[max(0,K-1-RepeatMaxShift)]); !align->repeat && bestI2 < bestI; bestI2++){/* back track best alignment from (bestI2,bestJ2) and compute its Pvalue */
	    // WAS	  for(bestI2 = max(1,bestI - RepeatMaxShift); !align->repeat && bestI2 < bestI; bestI2++){/* back track best alignment from (bestI2,bestJ2) and compute its Pvalue */
	    if(!(IMIN <= bestI2 && bestI2 <= IMAX && JMIN[bestI2] <= bestJ2 && bestJ2 <= JMAX[bestI2]))
	      continue;

	    double MinShiftKb = (Y[bestI] - Y[bestI2]) * RepeatKbShiftRatio; /* minimum offset shift along entire alignment */
	    double MaxShiftKb = (Y[bestI] - Y[bestI2]) * RepeatKbShiftRatio2; /* maximum offset shift along entire alignment (to avid creating new outliers) */

	    int I = bestI2;
	    int J = bestJ2;

	    /* First compute bestscore2, bestR2 with rightmost alignment at (I,J) */
	    int bestR2 = localtype;
	    PFLOAT bestscore2 = C;
	    PFLOAT Ascore = FLAT_A(A.score, I, J);/* NOTE: Ascore will be added to bestscore2 last, unlike for bestscore,bestI,bestJ */

	    int * RijxI = &FLAT_A(A.Rijx, I, 0), *RijyI = &FLAT_A(A.Rijy, I, 0);
	    PFLOAT Ydiff = Yr[I];// Y[N+1]-Y[I];
	    if(ENDINIT>=1 && I <= N-DELTA && J <= M-DELTA){/* must initialize Rijx,Rijy since this case was previously not initialized */
	      // WAS	  if(ENDINIT>=1 && (I <= N-DELTA || J <= M-DELTA)){/* must initialize Rijx,Rijy since this case was previously not initialized */
	      int Rijy = N+1;
	      int Rijx = M+1;

	      /* update Rijy by decrementing it until Y[Rijy] overlaps X */
	      FLOAT delta = Y[I] + X[M+1] - X[J];
	      while(Y[Rijy] > delta)
		Rijy--;
	      if(DEBUG>=2) assert(Rijy >= I);
	      RijyI[J] = Rijy;
	    
	      /* update Rijx by decrementing it until X[Rijx] overlaps Y */
	      delta = X[J] + Y[N+1] - Y[I];
	      while(X[Rijx] > delta)
		Rijx--;
	      if(DEBUG>=2)	assert(Rijx >= J);
	      RijxI[J] = Rijx;
	    }

	    PFLOAT score = Send(PMIN(Xr[J], Ydiff), (PFLOAT)(min(M, RijxI[J])+1-J + min(N, RijyI[J])+1-I));

	    if(ENDFIX){
	      if(RijxI[J] > M) {
		double a = M+1-I-J;
		int K;
		for(K= RijyI[J]; K > I; K--){
		  PFLOAT bound = Sbnd(Xr[J],Y[K]-Y[I], K+a);
		  if(bound <= score)
		    break;// HERE : see ENDFIX2 in refalign.cpp
		  score = bound;
		}
		RijyI[J] = K;
	      }
	      if(RijyI[J] > N) {
		double a = N+1-I-J;
		int K;
		for(K= RijxI[J];K > J;K--){
		  PFLOAT bound = Sbnd(X[K]-X[J], Ydiff,K+a);
		  if(bound <= score)
		    break;// HERE : see ENDFIX2 in refalign.cpp
		  score = bound;
		}
		RijxI[J] = K;
	      }
	    }
	    if(score > bestscore2){
	      bestscore2 = score;
	      bestR2 = -1;
	    }
	    bestscore2 += Ascore;

	    if(VERB && (VERB>=2 || EVERB>=2)){
              #pragma omp critical
	      {
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,M=%d,N=%d:checking for repeats\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,M,N);
		fflush(stdout);
	      }
	    }

	    /* Now backtrack to extract best alignment starting at (I,J)==(bestI2,bestJ2) constrained by requiring I to be less than in the best alignment for the same J */
	    int U = K-1;/* keep (I2 == align->sites1[U]) <= I, and J > (J2 == align->sites2[U]) to satisfy the repeat constraint
			   Also keep MinShiftKb <= Y[I2] - Yc[I] + X[J] - X[J2] <= MaxShiftKb */
	    int K2 = 0;

	    int bestG = -1, bestH = -1;
	    while(I >= 0){
	      Ilist2[K2] = I;
	      Jlist2[K2++] = J;

	      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);

	      /* update A.G[I][J] & A.H[I][J] by repeating recurrance computation constained by requiring I to be less than in the best alignment for the same J */
	      /* NOTE : this needs to be done even if DEFER_BP==0, since original A.G & A.H values would have been computed without the constraint */	    

	      if(FIX_CONSTRAINT && (FIX_CONSTRAINT>=2 || bestI-bestI2 >= 2)){
		while(U > 0 && align->sites1[U] > I)
		  U--;
	      } else {
		while(U > 0 && align->sites1[U-1] >= I)
		  U--;
	      }

	      PFLOAT score = FLAT_A(A.Uscore, I, J);
	      bestG = localtype;

	      int Gmin = max(I-DELTA,IMIN);
	      int T3 = U;

	      for(int G= I; --G >= Gmin;){
		int Hmin = max(J-DELTA,JMIN[G]);

		FLOAT XG;/* best estimate of original alignment location on X[] corresponding to Y[G] */
		int origT3 = T3;
		while(T3 < K-1 && align->sites1[T3] < G)
		  T3++;
		while(T3 > 0 && align->sites1[T3-1] >= G)
		  T3--;
		int I3 = align->sites1[T3];
		int J3 = align->sites2[T3];

		if(DEBUG>=2 && !(T3 >= K-1 || I3 >= G)){
		  #pragma omp critical
		  {
		    printf("While Checking for repeats:RepeatRec=%d,MaxShift=%d,LogPvalueRatio=%0.4f,%0.4f,KbShiftRatio=%0.3f,%0.3f,MaxDrop=%d,extend=%d,Lij2=%d,Rij2=%d,M=%d(score=%0.6f,LogPV=%0.2f,numpairs=%d,repeat=%d)\n",
			   RepeatRec,RepeatMaxShift,RepeatLogPvalueRatio,RepeatLogPvalueRatioPerLabel,RepeatKbShiftRatio,RepeatKbShiftRatio2,RepeatMaxLabelDrop,
			   extend,align->Lij2,align->Rij2,M,align->score,align->logPV,align->numpairs,align->repeat);
		    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,maptype=%d:\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID, bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,maptype);
		    printf("K=%d,U=%d,I=%d,J=%d,G=%d:T3=%d->%d:I3=%d,J3=%d\n",
			   K,U,I,J,G,origT3, T3,I3,J3);
		    fflush(stdout);
		    assert(T3 >= K-1 || I3 >= G);
		  }
		}
			
		int T4 = T3;
		while(T4 > 0 && align->sites1[T4] > G)
		  T4--;
		int I4 = align->sites1[T4];
		
		if(DEBUG>=2 && !(T4 <= 0 || I4 <= G)){
		  #pragma omp critical
		  {
		    printf("While Checking for repeats:RepeatRec=%d,MaxShift=%d,LogPvalueRatio=%0.4f,%0.4f,KbShiftRatio=%0.3f,%0.3f,MaxDrop=%d,extend=%d,Lij2=%d,Rij2=%d,M=%d(score=%0.6f,LogPV=%0.2f,numpairs=%d,repeat=%d)\n",
			   RepeatRec,RepeatMaxShift,RepeatLogPvalueRatio,RepeatLogPvalueRatioPerLabel,RepeatKbShiftRatio,RepeatKbShiftRatio2,RepeatMaxLabelDrop,
			   extend,align->Lij2,align->Rij2,M,align->score,align->logPV,align->numpairs,align->repeat);
		    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,maptype=%d:\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID, bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,maptype);
		    printf("K=%d,U=%d,I=%d,J=%d,G=%d:T3=%d->%d,T4=%d:I3=%d,J3=%d,I4=%d\n",
			   K,U,I,J,G,origT3, T3, T4, I3,J3, I4);
		    fflush(stdout);
		    assert(T4 <= 0 || I4 <= G);
		  }
		}
		origT3 = T3;

		if(G < I4){
		  int J4 = align->sites2[T4];
		  if(I4  <= G)
		    XG = X[J4];
		  else
		    XG = X[J4] - (Y[I4] - Y[G]);
		} else if(I3 > I4){
		  if(DEBUG>=2) assert(I4 <= G && G <= I3);
		  if(DEBUG>=2) assert(Y[I3] > Y[I4]);
		  FLOAT Gfrac = (Y[G] - Y[I4])/(Y[I3] - Y[I4]);
		  int J4 = align->sites2[T4];
		  XG = X[J4] + Gfrac * (X[J3] - X[J4]);
		} else  {
		  if(DEBUG>=2) assert(I3 == I4);
		  if(DEBUG>=2) assert(T3 == T4);
		  XG = X[J3];
		}

		if(FIX_CONSTRAINT && (FIX_CONSTRAINT>=2 || bestI-bestI2 >= 2)){
		  while(T3 > 0 && align->sites1[T3] > G)
		    T3--;
		} else {
		  while(T3 > 0 && align->sites1[T3-1] >= G)
		    T3--;
		}
		I3 = align->sites1[T3];
		J3 = align->sites2[T3];

		/* H must satisfy the constraint: H > J3 AND MinShiftKb <= Y[I3] - Y[G] + X[H] - XG <= MaxShiftKb */
		FLOAT XHmin = XG - Y[I3];
		FLOAT XHmax = XHmin;
		int HminG = max(Hmin,J3+1);
		if(I3 < G && T3 < K-1){
		  int I3 = align->sites1[T3+1];
		  int J3 = align->sites2[T3+1];
		  if(X[J3] - Y[I3] < XHmin)
		    XHmin = X[J3] - Y[I3];
		  if(X[J3] - Y[I3] > XHmax)
		    XHmax = X[J3] - Y[I3];
		  if(HminG > J3+1)
		    HminG = max(Hmin,J3+1);
		}
		if(G < I3 && T3 > 0){
		  int I3 = align->sites1[T3-1];
		  int J3 = align->sites2[T3-1];
		  if(X[J3] - Y[I3] < XHmin)
		    XHmin = X[J3] - Y[I3];
		  if(X[J3] - Y[I3] > XHmax)
		    XHmax = X[J3] - Y[I3];
		  if(HminG > J3+1)
		    HminG = max(Hmin,J3+1);
		}
		XHmin += MinShiftKb + Y[G]; /* X[H] must be >= this value */
		XHmax += MaxShiftKb + Y[G]; /* X[H] must be <= this value */

		PFLOAT deltaY = delY[I][G-I+DELTA];// Y[I] - Y[G] 
		if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
		PFLOAT *AscoreG = &FLAT_A(A.score, G, 0);
		// WAS107		PFLOAT a = J+I-G;
		int HmaxG = min(J-1, JMAX[G]);

		for(int H= HmaxG; H >= HminG; H--){
		  if(DEBUG>=2) assert(IMIN <= G && G <= IMAX && JMIN[G] <= H && H <= JMAX[G]);

		  if(X[H] > XHmax && H > HminG)
		    continue;
		  if(X[H] < XHmin)
		    break;

   	          if(DEBUG>=2 && !(1 <= H && H < J && J-H <= DELTA && J <= M)){
	            #pragma omp critical
	            {
	              printf("Yid=%d,Xid=%d,or=%d,sc=%d:I=%d,J=%d,G=%d,H=%d,M=%d\n",Yid,Xid,orientation,scaleID,I,J,G,H,M);
		      fflush(stdout);
		      assert(1 <= H && H < J && J-H <= DELTA && J <= M);
	            }
	          }
		  PFLOAT deltaX = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];// X[J] - X[H]
		  if(DEBUG>=2) assert(fabs(X[J] - X[H] - deltaX) < 0.001);
		  PFLOAT newscore = AscoreG[H] +  Sint(deltaX,deltaY,J-H,I-G,myoutlierBC);// WAS107 Sint(deltaX,deltaY,(PFLOAT)(a-H));
		  if(newscore > score){
		    score = newscore;
		    bestG = G;
		    bestH = H;
		  }
		}
	      }

	      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);

	      if(bestG < 0){/* Left End : value of bestH is undefined */
		if(DEFER_G >= 1 && score > C){
		  bestG = -1;
		  if(DEFER_G >= 2){/* check if actual value of Lijx,Lijy needs to be incremented */
		    int Lijx = FLAT_A(A.Lijx,I,J);
		    int Lijy = FLAT_A(A.Lijy,I,J);
		    PFLOAT SmIJ = Sm(J,I);
		    PFLOAT ascore = C;

		    if(Lijx <= 0){
		      ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I+1-max(1,Lijy))) + SmIJ;
		      ascore = max(ascore,C);
		      for(int K= Lijy; K < I; ) {
			PFLOAT bound = Sbnd(X[J],Y[I]-Y[K],(PFLOAT)(J+I-K)) + SmIJ;
			if(bound <= ascore)
			  break;
			ascore = bound;
			Lijy = ++K;
		      }
		      FLAT_A(A.Lijy, I, J) = Lijy;
		    } else {
		      if(DEBUG>=2) assert(Lijy <= 0);
		      ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-max(1,Lijx)+I)) + SmIJ;
		      ascore = max(ascore,C);
		      for(int K= Lijx; K < J; ){
			PFLOAT bound = Sbnd(X[J]-X[K],Y[I],(PFLOAT)(J+I-K)) + SmIJ;
			if(bound <= ascore)
			  break;
			ascore = bound;
			Lijx = ++K;
		      }
		      FLAT_A(A.Lijx, I, J) = Lijx;
		    }
		  }
		}
		FLAT_A(A.G, I, J) = bestG;
		J = -1;
	      } else {
		FLAT_A(A.G, I, J) = bestG;
		FLAT_A(A.H, I, J) = bestH;
		J = bestH;
	      }
	      I = bestG;

	      if(VERB && (VERB>=2 || EVERB>=2)){
                #pragma omp critical
		{
		  printf("  K2=%d:I=Ilist2[K2]=%d,J=Jlist[K2]=%d,A.G[I][J]=%d,A.H[I][J]=%d,score=%0.6f\n",
			 K2-1,Ilist2[K2-1],Jlist[K2-1],A.G[Ilist2[K2-1]][Jlist2[K2-1]],A.H[Ilist2[K2-1]][Jlist2[K2-1]],score);
		  fflush(stdout);
		}
	      }
	    }

	    if(DEBUG>=2){
	      if(DEBUG>=2 && bestR2 >= -1) assert(IMIN <= bestI2 && bestI2 <= IMAX && JMIN[bestI2] <= bestJ2 && bestJ2 <= JMAX[bestI2]);
	      int Rijy = (bestR2 <= -2) ? bestI2 : FLAT_A(A.Rijy, bestI2, bestJ2);
	      int Rijx = (bestR2 <= -2) ? bestJ2 : FLAT_A(A.Rijx, bestI2, bestJ2);
	      if(DEBUG>=2) assert(bestI2 <= Rijy && Rijy <= N+1);
	      if(DEBUG>=2) assert(bestJ2 <= Rijx && Rijx <= M+1);
	      
	      if(!(min(N,Rijy)-bestI2 >= 0)){
                #pragma omp critical
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,Rijy=%d,Rijx=%d,M=%d,N=%d\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,Rijy,Rijx,M,N);
		  for(int T2 = 0; T2 < K2; T2++)
		    printf("   T=%d/%d:I=%d,J=%d\n",T2,K2,Ilist2[K2-1-T2],Jlist2[K2-1-T2]);
		  fflush(stdout);
		  assert(min(N,Rijy)-bestI2 >= 0);
		}
	      }
	    }

	    double logPV2;
	    logPV2 = alignFP(A,bestI2,bestJ2,bestR2,Y,X,N,M,orientation,scaleID,Yid,Xid,delY,delX,bestscore2,/* (Yid==285 && Xid==22611 && orientation==0) ? 1 : 0 */ 0, 0, myoutlierExtend,myoutlierBC,tid,
			     JMIN,JMAX,IMIN,IMAX);

	    if(logPV2 >= max(align->logPV * RepeatLogPvalueRatio, minPvalue) +
	       (RepeatMaxLabelDrop < 0 || K2 >= align->numpairs - abs(bestI - bestI2) - RepeatMaxLabelDrop) &&
	       (logPV2 / (K2-1)) >= (align->logPV / (align->numpairs-1)) * RepeatLogPvalueRatioPerLabel){

	      align->repeat = bestI - bestI2;
	      if(DEBUG) assert(align->repeat > 0);
	      if(VERB>=2 || EVERB){
                #pragma omp critical
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.6f,numpairs=%d:bestI2=%d,bestJ2=%d,score2=%0.4f,logPV2=%0.6f (repeat=%d confirmed)\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,bestscore2,logPV2,align->repeat);
		  for(int T2 = 0; T2 < K2; T2++)
		    if(T2 > 0)
		      printf("   T=%d/%d:I=%d,J=%d (Y=%0.4f X=%0.4f)\n",T2,K2,Ilist2[K2-1-T2],Jlist2[K2-1-T2],
			     Y[Ilist2[K2-1-T2]] - Y[Ilist2[K2-T2]], X[Jlist2[K2-1-T2]] - X[Jlist2[K2-T2]]);
		    else
		      printf("   T=%d/%d:I=%d,J=%d\n",T2,K2,Ilist2[K2-1-T2],Jlist2[K2-1-T2]);
		  printf("  Best alignment:\n");
		  for(int T = 0; T < K; T++)
		    if(T > 0)
		      printf("   T=%d/%d:I=%d,J=%d (Y=%0.4f X=%0.4f)\n",T,K,Ilist[K-1-T],Jlist[K-1-T],
			     Y[Ilist[K-1-T]] - Y[Ilist[K-T]], X[Jlist[K-1-T]] - X[Jlist[K-T]]);
		    else
		      printf("   T=%d/%d:I=%d,J=%d\n",T,K,Ilist[K-1-T],Jlist[K-1-T]);
		  fflush(stdout);
		}
	      }
	    } else if(VERB>=2 || EVERB>=2){
              #pragma omp critical
	      {
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,score2=%0.4f,logPV2=%0.2f (No repeat at shiftI =%d)\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,bestI2,bestJ2,bestscore2,logPV2, bestI - bestI2);
		for(int T2 = 0; T2 < K2; T2++)
		  if(T2 > 0)
		    printf("   T=%d/%d:I=%d,J=%d (Y=%0.3f X=%0.3f)\n",T2,K2,Ilist2[K2-1-T2],Jlist2[K2-1-T2],
			   Y[Ilist2[K2-1-T2]] - Y[Ilist2[K2-T2]], X[Jlist2[K2-1-T2]] - X[Jlist2[K2-T2]]);
		  else
		    printf("   T=%d/%d:I=%d,J=%d\n",T2,K2,Ilist2[K2-1-T2],Jlist2[K2-1-T2]);
		printf("  Best alignment:\n");
		for(int T = 0; T < K; T++)
		  if(T > 0)
		    printf("   T=%d/%d:I=%d,J=%d (Y=%0.3f X=%0.3f)\n",T,K,Ilist[K-1-T],Jlist[K-1-T],
			   Y[Ilist[K-1-T]] - Y[Ilist[K-T]], X[Jlist[K-1-T]] - X[Jlist[K-T]]);
		  else
		    printf("   T=%d/%d:I=%d,J=%d\n",T,K,Ilist[K-1-T],Jlist[K-1-T]);
		fflush(stdout);
	      }
	    }
	  }

	  if(VERB>=2 && align->repeat > 0){
	    #pragma omp critical
	    {
	      repeatcnt++;

	      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:repeat detected(repeat=%d), repeatcnt=%d\n",
		     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,bestI,bestJ,align->score,align->logPV,align->numpairs,align->repeat, repeatcnt);
	      fflush(stdout);
	    }
	  } else if(DEBUG>=2 && align->repeat > 0){

	    #pragma omp atomic
	    repeatcnt++;

	  }

	}// if(RepeatRec ... )

      } /* else alignment site information will not be used */

      I = Ilist[K-1];// align->sites1[0];
      J = Jlist[K-1];// align->sites2[0];

    } //  if(bestscore > origscore)

    if(EVERB && (K <= 0 || align->numpairs <= 0)){
      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,sc=%d:score=%0.4f,logPV=%0.2f,pairs=%d,numpairs=%d\n",Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,scaleID,align->score,align->logPV,K,align->numpairs);
      fflush(stdout);
    }

    if(VERB >= ((LVERB||PVERB||EVERB) ? 1 : 2) && K > 0 && (VERB>=2 || EVERB || (align->score > ScoreThreshold && align->logPV >LogPvThreshold && align->numpairs >= AlignedSiteThreshold)
							  /* || (YYmap[Yid]->id == 1769LL && XXmap[Xid]->id==1031LL) */) /* && !orientation && bestscore > 0.0*/)
      #pragma omp critical
      {
      double trueoffset = XXmap[Xid]->startloc - YYmap[Yid]->startloc;// It would be better to use endloc to compute trueoffset
      double offset = Y[bestI]-X[bestJ];
      double delta = fabs(offset-trueoffset);
      double logPV = alignFP(A,bestI,bestJ,bestR,Y,X,N,M,orientation,scaleID,Yid,Xid,delY,delX,bestscore, 1, 1, myoutlierExtend,myoutlierBC,tid,JMIN,JMAX,IMIN,IMAX);
      int repeat,Lend,Rend;
      if(!(EVERB && bestscore <= origscore)){
	repeat = align->repeat;
	Lend = align->Lend;
	Rend = align->Rend;
      } else {/* EVERB && bestscore <= origscore */
	repeat = -2;
	Lend = I;
	Rend = bestR;
      }
      if(1 || delta > MAXDELTA*Ylambda){
	printf("Yid=%d(id=%lld,left=%d,right=%d,N=%d),Xid=%d(id=%lld,left=%d,right=%d,M=%d):score=%0.4f,logPV=%0.4f,delta=%0.2f:or=%d,sc=%d,pairs=%d,repeat=%d,Lend=%d,Rend=%d:align:or=%d,sc=%d,score=%0.6f,logPV=%0.2f\n",
	       Yid,YYmap[Yid]->id,YYmap[Yid]->origmap ? YYmap[Yid]->left[0] : 0, YYmap[Yid]->origmap ? YYmap[Yid]->right[0] : N+1, YYmap[Yid]->origmap ? YYmap[Yid]->origmap->numsite[0] : N,
	       Xid,XXmap[Xid]->id,XXmap[Xid]->origmap ? XXmap[Xid]->left[0] : 0, XXmap[Xid]->origmap ? XXmap[Xid]->right[0] : M+1, XXmap[Xid]->origmap ? XXmap[Xid]->origmap->numsite[0] : M,
	       bestscore,logPV,delta,orientation,scaleID,K,repeat,Lend,Rend,align->orientation,align->scaleID,align->score,align->logPV);

	int I= -1,J= -1,G= -1, H = -1;
	for(int T=0; T < K; T++, G=I, H=J){
	  if(!(VERB >= (EVERB ? 1 : 2)) || bestscore > origscore){
	    I = align->sites1[T];
	    J = align->sites2[T];
	  } else {
	    I = Ilist[K-1-T];
	    J = Jlist[K-1-T];
	  }
	  if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);

	  printf("T=%2d:I=%2d,J=%2d:",T,I,J);
	  if(T<=0){ /* Left End: display Send(min(X[J],Y[I]),J+1-max(1,Lijx),I+1-max(1,Lijy))+Sm(J,I), A[I,J] */
	    if(Lend <= -2){
	      printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I]=%6.3f,ChimScore=%0.4f,Sm(J,I)=%0.4f,A[I,J]=%0.6f,p->G=%d,p->H=%d\n",
		     A.Lijy[I][J],A.Lijx[I][J],X[J],Y[I],ChimScore,Sm(J,I),FLAT_A(A.score, I, J),A.G[I][J],A.H[I][J]);
	    } else if(X[J] <= Y[I] && A.Lijy[I][J] > 0 && Y[I]-Y[A.Lijy[I][J]-1] < X[J]){
	      printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I]=%6.3f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.4f,Sm(J,I)=%0.4f,A[I,J]=%0.6f,p->G=%d,p->H=%d\n",
		     A.Lijy[I][J],A.Lijx[I][J],X[J],Y[I],X[J],Y[I]-Y[A.Lijy[I][J]-1],J,I+1-A.Lijy[I][J],
		     Sbnd(X[J],Y[I]-Y[A.Lijy[I][J]-1],(PFLOAT)(J + I+1-A.Lijy[I][J])),Sm(J,I),FLAT_A(A.score, I, J),A.G[I][J],A.H[I][J]);
	      if(DEBUG) assert(A.Lijx[I][J] <= 0);
	    } else if(X[J] >= Y[I] && A.Lijx[I][J] > 0 && X[J]-X[A.Lijx[I][J]-1] < Y[I]){
	      printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I]=%6.3f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.4f,Sm(J,I)=%0.4f,A[I,J]=%0.6f\n",
		     A.Lijy[I][J],A.Lijx[I][J],X[J],Y[I],X[J]-X[A.Lijx[I][J]-1],Y[I],J+1-A.Lijx[I][J],I,
		     Sbnd(X[J]-X[A.Lijx[I][J]-1],Y[I],(PFLOAT)(J+1-A.Lijx[I][J] + I)),Sm(J,I),FLAT_A(A.score, I, J));
	      if(DEBUG) assert(A.Lijy[I][J] <= 0);
	    } else
	      printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I]=%6.3f,Send(%0.3f,%d,%d)=%0.4f,Sm(J,I)=%0.4f,A[I,J]=%0.6f\n",
		     A.Lijy[I][J],A.Lijx[I][J],X[J],Y[I],min(X[J],Y[I]),J+1-max(1,A.Lijx[I][J]),I+1-max(1,A.Lijy[I][J]),
		     Send(min(X[J],Y[I]),(PFLOAT)(J+1-max(1,A.Lijx[I][J]) + I+1- max(1,A.Lijy[I][J]))),Sm(J,I),FLAT_A(A.score, I, J));
	  } else {/* internal interval : display Sint(X[J]-X[H],Y[I]-Y[G],J-H,I-G,myoutlierBC) + Sm(J,I), A[I,J] */
	    if(DEBUG) assert(G>=0 && H>=0);
	    PFLOAT x = X[J]-X[H];
	    PFLOAT y = Y[I]-Y[G];
	    PFLOAT sint = Sint(x,y,J-H,I-G,myoutlierBC);
	    PFLOAT Pen,Bias;
	    SintDetail(x,y,J-H,I-G,J,I,Bias,Pen);
	    double myvar = varF + var * (x+y);
	    if(QUADRATIC_VARIANCE)
	      myvar += varR * (x*x+y*y);
	    //	    printf("G=%2d,H=%2d,X[J]=%6.3f,Y[I]=%6.3f,X[J]-X[H]=%0.3f,Y[I]-Y[G]=%0.3f,norm=%0.4f,Sint(X,Y,%d,%d)=%0.6f(Pen=%0.6f,Bias=%0.6f),A[I,J]=%0.6f\n",
	    printf("G=%2d,H=%2d,X[J]=%6.8f,Y[I]=%6.8f,X[J]-X[H]=%0.8f,Y[I]-Y[G]=%0.8f,norm=%0.8f,Sint(X,Y,%d,%d)=%0.8f(Pen=%0.8f,Bias=%0.6f),A[I,J]=%0.6f\n",
		   G,H,X[J],Y[I],x,y,(x-y)*(x-y)/myvar,J-H,I-G,sint,Pen,Bias,FLAT_A(A.score, I, J));
	    /*	    printf("   varF2=%0.6f,var2=%0.6f, MatchScore2=%0.6f,MissPenalty=%0.6f,resKB2=%0.3f,XYlen=%0.3f,Frate=%0.6f,m=%d,n=%d\n",
		   varF2,var2,MatchScore2,MissPenalty,resKB2,max(0.0,X[J]-X[H]-resKB2)+max(0.0,Y[I]-Y[G]-resKB2),Frate,J-H,I-G);
	    printf("   varF=%0.6f,var=%0.6f, MatchScore=%0.6f,MissPenalty=%0.6f,resKB2=%0.3f,XYlen=%0.3f,Frate=%0.6f\n",
	    varF,var,MatchScore,MissPenalty,resKB2,max(0.0,X[J]-X[H]-resKB2)+max(0.0,Y[I]-Y[G]-resKB2),Frate);*/
	  }
	}
	/* Right End: display Rij,Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I) */
	int Rijy = (I>= N-DELTA || J >= M-DELTA) ? A.Rijy[I][J] : N+1;
	int Rijx = (I>= N-DELTA || J >= M-DELTA) ? A.Rijx[I][J] : M+1;
	if(Rend <= -2){
	  printf("     I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,ChimScore=%0.6f,total score=%0.8e\n",
		 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I], ChimScore,FLAT_A(A.score, I, J)+ChimScore);
	} else if(X[M+1]-X[J] <= Y[N+1]-Y[I] && Rijy <= N && Y[Rijy+1]-Y[I] < X[M+1]-X[J]){
	  if(DEBUG) assert(Rijx > M);
	  printf("     I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
		 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I], X[M+1]-X[J],Y[Rijy+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
		 Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I), FLAT_A(A.score, I, J) + Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I));
	} else if(X[M+1]-X[J] >= Y[N+1]-Y[I] && Rijx <= M && X[Rijx+1]-X[J] < Y[N+1]-Y[I]){
	  if(DEBUG) assert(Rijy > N);
	  printf("     I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
		 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
		 X[Rijx+1]-X[J],Y[N+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
		 Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I),
		 FLAT_A(A.score, I, J) + Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I));
	} else
	  printf("     I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Send(%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
		 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
		 min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,
		 Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I),
		 FLAT_A(A.score, I, J) + Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I));

	fflush(stdout);
      }
    }

    if(!EVERB || bestscore > origscore){

      if(PAIRSPLIT_EXTEND && PairSplit){
	align->refstart = align->sites1[0];
	align->start = align->sites2[0];
	align->refend = align->sites1[align->numpairs-1];
	align->end = align->sites2[align->numpairs-1];
	if(DEBUG>=2) assert(align->refstart <= align->refend);
      }

      if(Psplit > 0.0 && PairSplit){
	/* check if alignment has an internal region with outscore < -log(1/Psplit) : if so limit alignment to larger 
	   +ve scoring alignment region + end penalty of ChimScore (other region will be found after chimeric split, if large enough) */
	int U = align->numpairs;
	FLOAT *pscore = align->outscore;
	FLOAT *iscore = align->iscore;
	double score = 0.0;/* cumulative outscore */
	int iLWM = -1;/* most recent proposed alignment index break (-1 = left end, U = right end)) :
			 iLWM is either left end (-1) or at the right end of an internal region with outscore < -log(1/Psplit) (AND site interval >= PairsplitMinOutlier) 
			                     AND (if OUTLIER_FIXES) no subinterval with score >= pthreshold */
	double LWM = 0.0;/* sum(pscore[0..iLWM]) */
	int iHWM = -1;/* highest values of sum(pscore[0..iHWM]) with iLWM <= iHWM */
	double HWM = 0.0;/* sum(pscore[0..iHWM]) */
	double best = 0.0;/* best score interval = sum(pscore[bestLeft + 1 .. bestRight]) */
	int bestLeft= -1,bestRight= -1;/* range of aligned site index with best score interval */
	double threshold = log(1.0/Psplit);/* If -outlierBC, this threshold is increased by |OutlierPenaltyBC[m]+OutlierPenaltyBC[n]| for -ve scoring regions. 
					      The threshold for +ve scoring regions is unchanged */
	if(LVERB && align->logPV > LogPvThreshold && align->score > ScoreThreshold && align->numpairs >= AlignedSiteThreshold){
	  printf("Looking for splits:numpairs=%d:PairsplitMinOutlier=%d,myoutlierBC=%d\n",U,PairsplitMinOutlier,myoutlierBC);
	  fflush(stdout);
	}
	for(int T = 0; T <= U;T++){
	  /* PRE: -1 <= iLWM <= iHWM <= T-1 :
	     LWM = sum(pscore[0..iLWM])
	     HWM = sum(pscore[0..iHWM])
	     score = sum(pscore[0..T-1])
	     NOTE : If iLWM >= 0 : iLWM is at the right end of the last internal region (right end < T) with outscore < -log(1/Psplit) (AND site interval >= PairsplitMinOutlier) OR left end if none so far
	  */

	  double inc = pscore[T];
	  if(inc < 0.0 && PairsplitMinOutlier && 0 < T && T < U){
	    int IL = align->sites1[T-1], IR = align->sites1[T];
	    int JL = align->sites2[T-1], JR = align->sites2[T];
	    if(IR - IL < PairsplitMinOutlier && JR - JL < PairsplitMinOutlier)
	      inc = iscore[T] + 1e-3;
	  }
	
	  if(DEBUG>=2) assert(isfinite(pscore[T]));

	  if(VERB>=3 || EVERB || (LVERB && align->logPV > LogPvThreshold && align->score > ScoreThreshold && align->numpairs >= AlignedSiteThreshold)){
	    printf("T=%d:iLWM=%d,iHWM=%d:LWM=%0.3f,HWM=%0.3f,score=%0.3f,inc=%0.3f:best=%0.3f(%d..%d)\n",T,iLWM,iHWM,LWM,HWM,score,inc,best,bestLeft,bestRight);
	    fflush(stdout);
	  }

	  if(inc > 0.0 && iHWM >= 0){      /* check if HWM - score > threshold : If so, schedule break in alignment from iHWM to T-1 */
	    double nthreshold = threshold;
	    if(myoutlierBC) {
	      int IL = align->sites1[iHWM], IR = align->sites1[T-1];
	      int JL = align->sites2[iHWM], JR = align->sites2[T-1];
	      nthreshold -= OutlierPenaltyBC[IR-IL] + OutlierPenaltyBC[JR-JL];
	    }
	    if(HWM - score > nthreshold){
	      if(HWM - LWM > best){ /* If HWM-LWM > best : record interval from iLWM to iHWM as new best local alignment interval */
		bestLeft = iLWM;
		bestRight = iHWM;
		best = HWM - LWM;
	      }
	      iLWM = iHWM = T-1;
	      LWM = HWM = score;
	    }
	  }
	  score += inc;
	  if(score < LWM && !(OUTLIER_FIXES && HWM - LWM > threshold)){/* extends the previous (potential) split interval (bestRight to iLWM) to T */
	    if(!(OUTLIER_FIXES && HWM - LWM > threshold)){
	      LWM = HWM = score;
	      iLWM = iHWM = T;
	    }
	  } else if(score > HWM){
	    HWM = score;
	    iHWM = T;
	  }

	  /* POST: iLWM <= iHWM <= T : 
	     LWM = sum(pscore[0..iLWM])
	     HWM = sum(pscore[0..iHWM])
	     score = sum(pscore[0..T]) */
	}
	// NOTE : The following assertion can fail if Poutlier > 0 since cumulative score does not include internal outlier correction :
	if(DEBUG && !Poutlier && !(HWM-score <= fabs(ChimScore) + (USE_PFLOAT ? 1e-3 : 1e-8))){
	  printf("\nWhile localizing alignment:Yid=%d(id=%lld,left=%d,right=%d,N=%d),Xid=%d(id=%lld,left=%d,right=%d,M=%d):align->score=%0.8f->%0.8f,bestLeft=%d(I=%d,J=%d),bestRight=%d(I=%d,J=%d),U=%d\n",
		 Yid,YYmap[Yid]->id,YYmap[Yid]->origmap ? YYmap[Yid]->left[0] : 0, YYmap[Yid]->origmap ? YYmap[Yid]->right[0] : N+1,N,
		 Xid,XXmap[Xid]->id,XXmap[Xid]->origmap ? XXmap[Xid]->left[0] : 0, XXmap[Xid]->origmap ? XXmap[Xid]->right[0] : M+1,M,
		 align->score,best+(bestLeft >= 0 ? ChimScore:0.0) + (bestRight < U ? ChimScore : 0),
		 bestLeft,bestLeft>=0 ? align->sites1[bestLeft] : -1, bestLeft>=0 ? align->sites2[bestLeft] : -1, 
		 bestRight,bestRight<U ? align->sites1[bestRight] : -1, bestRight<U ? align->sites2[bestRight] : -1 , U);
	  printf("iLWM=%d,iHWM=%d,U=%d:LWM=%0.8f,HWM=%0.8f,score=%0.8f,best=%0.8f,HWM-score=%0.8f,ChimScore=%0.8f(abs=%0.8f),threshold=%0.8f\n",
		 iLWM,iHWM,U,LWM,HWM,score,best,HWM-score,ChimScore,fabs(ChimScore),threshold);
	  double cumscore = 0.0;
	  for(int t = 0; t <= U; t++){
	    cumscore += pscore[t];
	    printf("pscore[%d]=%0.8f,cum score=%0.8f\n",t,(double)pscore[t],cumscore);
	  }
	  fflush(stdout);
	
	  assert(HWM-score <= fabs(ChimScore) + (USE_PFLOAT ? 1e-3:1e-8));// otherwise alignment would have terminated at HWM (but see NOTE above)
	}

	if(OUTLIER_FIXES && iHWM >= 0){
	  if(HWM - score > threshold){
	    if(HWM - LWM > best){ /* If HWM-LWM > best : record interval from iLWM to iHWM as new best local alignment interval */
	      bestLeft = iLWM;
	      bestRight = iHWM;
	      best = HWM - LWM;
	    }
	    iLWM = iHWM = U;
	    LWM = HWM = score;
	  }
	}

	if(align->Rend <= -2){
	  if(score - ChimScore - LWM > best){
	    bestLeft = iLWM;
	    bestRight = U-1;
	    best = score - ChimScore - LWM;
	  }
	} else {
	  if(score - LWM > best){
	    bestLeft = iLWM;
	    bestRight = U;
	    best = score - LWM;
	  }
	}

	if(VERB>=3 || EVERB || (LVERB && align->logPV > LogPvThreshold && align->score > ScoreThreshold && align->numpairs >= AlignedSiteThreshold)){
	  printf("  Final: iLWM=%d,iHWM=%d:LWM=%0.3f,HWM=%0.3f,score=%0.3f:best=%0.3f(%d..%d)U=%d\n",iLWM,iHWM,LWM,HWM,score,best,bestLeft,bestRight,U);
	  fflush(stdout);
	}

	/* update alignment to only include the best local alignment interval bestLeft .. bestRight */
	if(align->numpairs >= 2 && min(U-1,bestRight) >= max(0,bestLeft)){
	  if(DEBUG && PoutlierEnd <= Psplit && !(bestLeft >= (align->Lend <= -2 ? 0 : -1))) {
	    printf("While localizing alignment:Yid=%d(id=%lld,left=%d,right=%d,N=%d),Xid=%d(id=%lld,left=%d,right=%d,M=%d):align->score=%0.6f -> %0.6f,bestLeft=%d(I=%d,J=%d),bestRight=%d(I=%d,J=%d),U=%d\n",
		   Yid,YYmap[Yid]->id,YYmap[Yid]->origmap ? YYmap[Yid]->left[0] : 0, YYmap[Yid]->origmap ? YYmap[Yid]->right[0] : N+1,N,
		   Xid,XXmap[Xid]->id,XXmap[Xid]->origmap ? XXmap[Xid]->left[0] : 0, XXmap[Xid]->origmap ? XXmap[Xid]->right[0] : M+1,M,
		   align->score,best+((bestLeft >= 0) ? ChimScore:0.0) + ((bestRight < U) ? ChimScore : 0),
		   bestLeft, (bestLeft>=0) ? align->sites1[bestLeft] : -1, (bestLeft>=0) ? align->sites2[bestLeft] : -1, 
		   bestRight, (bestRight<U) ? align->sites1[bestRight] : -1, (bestRight<U) ? align->sites2[bestRight] : -1 , U);
	    printf("iLWM=%d,iHWM=%d,LWM=%0.6f,HWM=%0.6f,score=%0.6f,Psplit=%0.4e,PoutlierEnd=%0.4e\n",iLWM,iHWM,LWM,HWM,score,Psplit,PoutlierEnd);
	    printf("align->Lend=%d,align->Rend=%d:outscore[0]=%0.6f,outscore[U=%d]=%0.6f\n",
		   align->Lend,align->Rend,align->outscore[0],U,align->outscore[U]);
	    fflush(stdout);
	    assert(bestLeft >= (align->Lend <= -2 ? 0 : -1));
	  }
	  if(DEBUG && PoutlierEnd <= Psplit) assert(bestRight <= (align->Rend <= -2 ? U-1 : U));
	}

	if((bestLeft >= 0 || bestRight < U) && min(U-1,bestRight) >= max(0,bestLeft)){
	  if(PAIRSPLIT_EXTEND){
	    align->refstart = (bestLeft >= 0) ? align->sites1[bestLeft] : align->sites1[0];
	    align->start = (bestLeft >= 0) ? align->sites2[bestLeft] : align->sites2[0];
	    align->refend = (bestRight < U) ? align->sites1[bestRight] : align->sites1[U-1];
	    align->end = (bestRight < U) ? align->sites2[bestRight] : align->sites2[U-1];
	    if(DEBUG>=2) assert(align->refstart <= align->refend);

	    if(bestLeft > 0){/* extend refstart,start into -ve scoring region to the left */
	      double score = 0.0;
	      double LWM = 0.0;
	      int iLWM = bestLeft;
	      for(int T = bestLeft; T > 0; T--){
		score += pscore[T];
		if(score >= LWM + threshold)
		  break;
		if(score < LWM){
		  LWM = score;
		  iLWM = T-1;
		} 
	      }
	      if(DEBUG) assert(align->sites1[iLWM] <= align->refstart);
	      if(DEBUG) assert(align->sites2[iLWM] <= align->start);
	      align->refstart = align->sites1[iLWM];
	      align->start = align->sites2[iLWM];
	    }
	    if(bestRight < U-1){/* extend refend,end into -ve scoring region to the right */
	      double score = 0.0;
	      double LWM = 0.0;
	      int iLWM = bestRight;
	      for(int T = bestRight; T < U-1; T++ /* WAS T-- */){
		score += pscore[T+1];
		if(score > LWM + threshold)
		  break;
		if(score < LWM){
		  LWM = score;
		  iLWM = T+1;
		}
	      }
	      if(DEBUG) assert(align->sites1[iLWM] >= align->refend);
	      if(DEBUG) assert(align->sites2[iLWM] >= align->end);
	      align->refend = align->sites1[iLWM];
	      align->end = align->sites2[iLWM];
	    }
	  }// if(PAIRSPLIT_EXTEND)

	  /* Update pscore[0...U] to pscore[0..min(U-1,bestRight)-max(0,bestLeft)+1] */
	  double oldscore = align->score;
	  double oldlogPV = align->logPV;
	  align->score = best + (bestLeft >= 0 ? ChimScore : 0) + (bestRight < U ? ChimScore : 0);

	  if(DEBUG) assert(min(U-1,bestRight)-max(0,bestLeft)+1 <= U);
	  if(bestLeft >= 0){
	    align->Lend = -2;
	    pscore[0] = align->iscore[0] = ChimScore;
	    if(bestLeft > 0){
	      align->sites1[0] = align->sites1[bestLeft];
	      align->sites2[0] = align->sites2[bestLeft];
	      for(int T = bestLeft+1; T <= bestRight; T++){
		pscore[T-bestLeft] = pscore[T];
		align->iscore[T-bestLeft] = align->iscore[T];
	  
		if(T < U){
		  align->sites1[T-bestLeft] = align->sites1[T];
		  align->sites2[T-bestLeft] = align->sites2[T];
		}
	      }
	    }
	  }
	  if(bestRight < U){
	    align->Rend = -2;
	    pscore[bestRight-max(0,bestLeft)+1] = align->iscore[bestRight-max(0,bestLeft)+1] = ChimScore;
	  }
	  align->numpairs = min(U-1,bestRight)-max(0,bestLeft)+1;
	  if(DEBUG) assert(align->numpairs >= 1);
	  if(align->Lend <= -2){
	    align->Lij1 = align->sites1[0];
	    align->Lij2 = align->sites2[0];
	  }
	  if(align->Rend <= -2){
	    align->Rij1 = align->sites1[align->numpairs-1];
	    align->Rij2 = align->sites2[align->numpairs-1];
	  }

	  /* update align->logPV */
	  bestI = align->sites1[align->numpairs-1];
	  bestJ = align->sites2[align->numpairs-1];
	  bestR = align->Rend;
	  bestscore = align->score;
	  FLAT_A(A.G, align->sites1[0], align->sites2[0]) = -2;
	  align->logPV = alignFP(A,bestI,bestJ,bestR,Y,X,N,M,orientation,scaleID,Yid,Xid,delY,delX,bestscore, 0, 1, myoutlierExtend,myoutlierBC,tid,JMIN,JMAX,IMIN,IMAX);
	  if(LVERB){
	    printf("Localizing alignment:Yid=%d(id=%lld,left=%d,right=%d,N=%d),Xid=%d(id=%lld,left=%d,right=%d,M=%d):align->score=%0.6f->%0.6f,logPV=%0.4f->%0.4f,bestLeft=%d(I=%d,J=%d),bestRight=%d(I=%d,J=%d),U=%d\n",
		   Yid,YYmap[Yid]->id,YYmap[Yid]->origmap ? YYmap[Yid]->left[0] : 0, YYmap[Yid]->origmap ? YYmap[Yid]->right[0] : N+1,N,
		   Xid,XXmap[Xid]->id,XXmap[Xid]->origmap ? XXmap[Xid]->left[0] : 0, XXmap[Xid]->origmap ? XXmap[Xid]->right[0] : M+1,M,
		   oldscore,align->score,oldlogPV,align->logPV,
		   bestLeft, (bestLeft>=0) ? align->sites1[0] : 0, (bestLeft>=0) ? align->sites2[0] : 0, 
		   bestRight, (bestRight<U) ? align->sites1[align->numpairs-1] : N+1, (bestRight<U) ? align->sites2[align->numpairs-1] : M+1 , U);
	    if(PAIRSPLIT_EXTEND)
	      printf("\t extLeft=%d,%d,extRight=%d,%d,iLWM=%d,iHWM=%d,LWM=%0.6f,HWM=%0.6f,score=%0.6f\n",align->refstart,align->start,align->refend,align->end,iLWM,iHWM,LWM,HWM,score);
	    else {
	      int IL = align->sites1[0], IR = align->sites1[align->numpairs-1];
	      FLOAT *origY = (YYmap[Yid]->origmap ? YYmap[Yid]->origmap : YYmap[Yid])->site[0];
	      int left = YYmap[Yid]->origmap ? YYmap[Yid]->left[0] : 0;
	      printf("\t iLWM=%d,iHWM=%d,LWM=%0.6f,HWM=%0.6f,score=%0.6f:alignment=(I=%d,J=%d) to (I=%d,J=%d),ref=%0.4f to %0.4f\n",iLWM,iHWM,LWM,HWM,score,
		     align->sites1[0],align->sites2[0],align->sites1[align->numpairs-1],align->sites2[align->numpairs-1], origY[IL+left],origY[IR+left]);
	    }
	    fflush(stdout);
	  }
	}
	if(PAIRSPLIT_EXTEND && DEBUG){
	  assert(0 < align->refstart);
	  assert(0 < align->start);
	  assert(align->refstart <= align->refend);
	  assert(align->start <= align->end);
	  assert(align->refend <= N);
	  assert(align->end <= M);
	}

	if(origalign && origalign->score > align->score){
	  if(EVERB){
	    printf("Restoring origalign:Yid=%d,Xid=%d,or=%d,sc=%d,numpairs=%d,score=%0.6f,logPV=%0.2f,repeat=%d\n",
		   origalign->mapid1,origalign->mapid2,origalign->orientation,origalign->scaleID,origalign->numpairs,origalign->score,origalign->logPV,origalign->repeat);
	    fflush(stdout);
	  }
	  
	  // same as align->allfree(), but preserve align->align2[]
	  if(align->allocated_size > 0 && align->block_id < 0){
	    delete [] align->sites1; align->sites1 = NULL;
	    delete [] align->sitesK1; align->sitesK1 = NULL;
	    delete [] align->sites2; align->sites2 = NULL;
	    delete [] align->iscore; align->iscore = NULL;
	    delete [] align->outscore; align->outscore = NULL;
	  }
	  align->allocated_size = 0;
	  copy(align, origalign, 1, 1);
	}
      }// if(Psplit > 0.0 && PairSplit)
    } // if(!EVERB || bestscore > origscore)
    if(origalign)
      delete [] origalign;
  } else if(bestscore <= align->score){
    if(Palindromic_Delta > 0.0 && PairSplit && orientation != align->orientation && bestscore > MINSCORE && 
       (!align->align2 || bestscore > align->align2->score + 1e-8)){ /* backtrack through array to retrieve best alignment in current orientation as 2nd best alignment */

      int I = bestI;
      int J = bestJ;

      int K = 0;
    
      iscore[0] = bestscore - FLAT_A(A.score, I, J);
      outscore[0] = iscore[0];
      while(I >= 0){
	if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && JMIN[I] <= J && J <= JMAX[I]);
	if(DEBUG>=2){
	  assert(I >= 1 && I <= N);
	  assert(J >= 1 && J <= M);
	  assert(K < M);
	}

	Ilist[K] = I;
	Jlist[K] = J;
	int G = FLAT_A(A.G, I, J);
	int H = FLAT_A(A.H, I, J);
	if(G > 0){
	  if(DEBUG>=2) assert(H > 0);
	  iscore[K+1] = FLAT_A(A.score, I, J) - FLAT_A(A.score, G, H);
	  PFLOAT y = Y[I] - Y[G];
	  PFLOAT x = X[J] - X[H];
	  if(Poutlier > 0.0){
	    PFLOAT Bias,Pen;
	    SintDetail(x, y, J - H, I - G, J, I, Bias, Pen);
	    outscore[K+1] = OutlierBias + Bias + Pen;
	  } else
	    outscore[K+1] = iscore[K+1];
	  if(DEBUG>=2) assert(iscore[K+1] >= outscore[K+1] - 0.001);
	} else {/* left end segment */
	  iscore[K+1] = FLAT_A(A.score, I, J);
	  outscore[K+1] = iscore[K+1];
	}
	K++;
	I = G;
	J = H;
      }
      if(DEBUG) assert(K>0);

      delete [] align->align2; align->align2 = NULL;
      Calign *align2 = align->align2 = new Calign[1];
      copy(align2, align, 1, 1);
      align2->score = bestscore;
      align2->logPV = alignFP(A,bestI,bestJ,bestR,Y,X,N,M,orientation,scaleID,Yid,Xid,delY,delX,bestscore, /* (Yid==285 && Xid==22611 && orientation==0) ? 1 : 0 */ 0, 1, myoutlierExtend,myoutlierBC,tid,
			      JMIN,JMAX,IMIN,IMAX);
      align2->orientation = orientation;
      align2->scaleID = scaleID;
      //      align->I2 = bestI;
      //      align->J2 = bestJ;
      align2->Lend = I;
      align2->Rend = bestR;
      //      align->numpairs2 = K;
    }
  }

  if(mem)
    free(mem);
  if(MEMCHECK && RA_MIN_TIME && RA_MIN_MEM < 8) ReleaseMemory(Mem,nMemSiz,tid, Yid,Xid,N,M,Fmem);
}


/** find the best alignment between refmap and nanomap (in either orientation for nanomap) and store the
   resulting alignment in align 
   Checks Chashpair for alignment modifications and runs pairalignXY() to 
   find alignment in each orientation.
   
   */
void pairalign(Cmap *Ymap,
	       Cmap *Xmap,
	       Calign * &galign,/**< reference to pointer located in alignment[] array : may be NULL */
	       int Yid,int Xid,
	       Chashpair *phash,/**< If != 0 pointer to hashtable entry */
	       AL &A, /**< A.<field>[I=1..N][J=1..M] : only the 1st dimension (pointers) is allocated */
	       PFLOAT * &Fmem,/**< preallocated memory to be used by A.score[][],A.Uscore[][] */
	       int * &Imem, /**< preallocated memory to be used by A.<field>[][], for <field> = G,H,Rijy,Lijy,Rijx,Lijx */
	       CThreadMem *Mem, // typically Fmem == Mem->Fmem AND Imem == Mem->Imem
	       long long &StrideMax,/* If (MIN_MEM && DIAGONAL>=2 && hashdelta) : the current allocated number of elements per array */
               int tid
	       )
{
  if(VERB>=2){
    #pragma omp critical
    {
      printf("pairalign:Yid=%d,Xid=%d,galign=%p\n",Yid,Xid,galign);
      fflush(stdout);
    }
  }

  if(colors != 1){
    printf("alignment for multi-color data not yet implemented\n");
    fflush(stdout);exit(1);
  }

  if(DEBUG && !(Xmap->deltaX[0] && Xmap->deltaR[0] && Ymap->deltaX[0])){
    #pragma omp critical
    {
      printf("Yid=%d(%lld,%lld),Xid=%d(%lld,%lld):Xmap->deltaX[0]=%p,Xmap->deltaR[0]=%p,Ymap->deltaX[0]=%p\n",
	     Yid,Ymap->id,YYmap[Yid]->id,Xid,Xmap->id,XXmap[Xid]->id,Xmap->deltaX[0],Xmap->deltaR[0],Ymap->deltaX[0]);
      fflush(stdout);
      assert(Xmap->deltaX[0] && Xmap->deltaR[0] && Ymap->deltaX[0]);
    }
  }

  int N = Ymap->numsite[0];
  int M = Xmap->numsite[0];
  FLOAT *X = Xmap->site[0];
  FLOAT *Y = Ymap->site[0];

  if(DEBUG>=2 && galign && galign->numpairs > 0){
    assert(galign->sites1[0] >= 1 && galign->sites1[galign->numpairs-1] <= N);
    assert(galign->sites2[0] >= 1 && galign->sites2[galign->numpairs-1] <= M);
  }

  Calign *align, Lalign;
  if(PairSplit){
    if(DEBUG) assert(galign != 0);
    align = galign;
  } else /* use a local temporary variable : only save it to galign if score is good */
    align = &Lalign;

  align->mapid1 = Yid;
  align->mapid2 = Xid;
  align->score = MINSCORE;
  if(PairSplit){
    delete [] align->align2;
    align->align2 = NULL;
  }
  //  align->logPV = 0.0;
  align->numpairs = 0;
	       
  if(DEBUG>=2 && !(align->mapid1 == Yid)){
    #pragma omp critical
    {
      printf("pairalign:Yid=%d(id=%lld),Xid=%d(id=%lld):N=%d,M=%d\n",
	     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,N,M);
      printf("align->mapid1=%d,align->mapid2=%d\n",align->mapid1,align->mapid2);
      fflush(stdout);
      assert(align->mapid1 == Yid);
    }
  }

  if(N <= 0 || M <= 0)
    return;/* no alignment possible */

  if(DEBUG && RA_MIN_TIME == 0) assert(Fmem != NULL && Imem != NULL);

  // WAS A.M_stride = ((M + 0x1f) & ~0x1f) + 16LL; /* This is necessary to make better use of set-associative cache ??? */
  A.M_stride = M;
  A.array_stride = N * A.M_stride + 16LL;// Arrays are offset by 4 elements to allow vectors to write up to 4 elements before start and 12 after end 

  size_t Fsiz = A.array_stride * DPsizF + 16LL;
  size_t Isiz = A.array_stride * DPsizI + 16LL;
  //  size_t Fsiz = N * ((M + 64LL) + 16LL) * DPsizF + 16LL;
  //  size_t Isiz = N * ((M + 64LL) + 16LL) * DPsizI + 16LL;
  long long nMemSiz = ((Fsiz * sizeof(PFLOAT) + Isiz * sizeof(int)) + (PAGE-1)) & ~(PAGE-1);

  if(RA_MIN_TIME && MaxMemSiz > 0 && (!MEMCHECK || RA_MIN_MEM >= 8))/* Pause this thread if not enough memory is available. */
    WaitForMemory(Mem,nMemSiz,tid,Yid,Xid,N,M,Fmem);

  if(RA_MIN_MEM >= 8){// always reallocate just the right amount of memory using mmap
    char *newmem = (char *)mmap(NULL, nMemSiz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if(newmem == MAP_FAILED){
      int eno = errno;
      char *err = strerror(eno);
	
      #pragma omp critical
      {
        printf("pairalign: mmap of %lld bytes failed:errno=%d:%s\n", nMemSiz,eno,err);
	dumpmemmap();
	fflush(stdout);exit(1);
      }
    }
    
    Fmem = (PFLOAT *) newmem;
    Imem = (int *) &newmem[Fsiz * sizeof(PFLOAT)];
    Mem->MemSiz = nMemSiz;
  }

  if(DEBUG && !(Fmem != 0 && Imem != 0)){
    #pragma omp critical
    {
      printf("tid=%d: Fmem= %p, Imem= %p, MemSiz= %lld (RA_MIN_MEM=%d,RA_MIN_TIME=%d): wall time= %0.6f\n",tid,Fmem,Imem,Mem->MemSiz,RA_MIN_MEM,RA_MIN_TIME,wtime());
      fflush(stdout);
      assert(Fmem != 0 && Imem != 0);
    }
  }

  /* re-allocated DP arrays A.<field>[I=1..N][J=1..M] */
  A.G[0] = NULL;
  A.H[0] = NULL;
  A.Rijy[0] = NULL;
  A.Lijy[0] = NULL;
  A.Rijx[0] = NULL;
  A.Lijx[0] = NULL;
  A.score[0] = NULL;
  A.Uscore[0] = NULL;
  
  if(DEBUG && N > maxN){
    #pragma omp critical
    {
      fprintf(stderr, "Yid=%d(id=%lld),Xid=%d(id=%lld): M=%d,maxM=%lld,N=%d,maxN=%lld,numY=%d\n", Yid,Ymap->id,Xid,Xmap->id,M,maxM,N, maxN,numY);
      fflush(stderr);
      assert(N <= maxN);
    }
  }
  if(DEBUG && M > maxM){
    #pragma omp critical
    {
      printf("Yid=%d(id=%lld),Xid=%d(id=%lld): M=%d,maxM=%lld,N=%d,maxN=%lld,numX=%d\n", Yid,Ymap->id,Xid,Xmap->id,M, maxM,N,maxN,numX);
      fflush(stdout);
      assert(M <= maxM);
    }
  }

  //prefetch_array(Imem, A.array_stride*4*6+64, _MM_HINT_T1);
  //prefetch_array(Fmem, A.array_stride*4*2+64, _MM_HINT_T1);
 
#pragma ivdep
  for(int I = 1; I <= N; I++){
    A.G[I] = &Imem[(I-1) * A.M_stride + 3];
    A.H[I] = &Imem[(I-1) * A.M_stride + 3 + A.array_stride];
    A.Rijy[I] = &Imem[(I-1) * A.M_stride + 3 + A.array_stride * 2];
    A.Lijy[I] = &Imem[(I-1) * A.M_stride + 3 + A.array_stride * 3];
    A.Rijx[I] = &Imem[(I-1) * A.M_stride + 3 + A.array_stride * 4];
    A.Lijx[I] = &Imem[(I-1) * A.M_stride + 3 + A.array_stride * 5];

    A.score[I] = &Fmem[(I-1)*A.M_stride + 3];
    A.Uscore[I] = &Fmem[(I-1)*A.M_stride + 3 + A.array_stride];
  }

  if(DEBUG>=2){
    for(size_t i = 0; i < Fsiz; i++)
      Fmem[i] = aNaN;
    for(size_t i = 0; i < Isiz; i++)
#ifdef VALGRIND
      Imem[i] = -555555555;
#else
      Imem[i] = -555;
#endif
  } 
  // NOTE : touching all memory defeats taking advantage of sparse/gradual memory usage

#ifdef VALGRIND
  VALGRIND_MAKE_MEM_UNDEFINED(Fmem,Fsiz*sizeof(PFLOAT));
  VALGRIND_MAKE_MEM_UNDEFINED(Imem,Isiz*sizeof(int));
#endif
#ifdef CLANG
#if defined(__has_feature)
#  if __has_feature(memory_sanitizer)  // code that builds only under MemorySanitizer
  __msan_poison(Fmem,(Fsiz*sizeof(PFLOAT)));
  __msan_poison(Imem,(Isiz*sizeof(int)));
#  endif
#endif
#endif
  
  if(PAIRMERGE_FIX && PairMerge){
    align->Nrange = 0;/* used to flag valid alignments based on -pairmerge outlier/endoutlier constraints */
    //    goutlierExtend = outlierExtend;
    //    goutlierBC = outlierBC;
  }

  FLOAT *Xrev = NULL;
  if (M > MAX_ALLOCA)
    Xrev = new FLOAT[M+2];
  else
    Xrev = (FLOAT *)alloca((M+2)*sizeof(FLOAT));

  PFLOAT **deltaX = NULL, **deltaR = NULL, *fmemblock = NULL;
  if(NumScaleFactor > 1){
    if(((USE_SSE && USE_AVX) || USE_MIC) && USE_PFLOAT){// allocate local arrays deltaX[1..DELTA][2..M] and deltaR[1..DELTA][2..M]
      deltaX = new PFLOAT*[2*(DELTA+1)];
      deltaR = &deltaX[DELTA+1];
      
      int M_stride = COMPUTE_STRIDE4(M-1);
      if(VERB>=2 && tid==0){
	printf("M=%d,M_stride=%d\n",M,M_stride);
	fflush(stdout);
      }

      fmemblock = new PFLOAT[2 * DELTA * M_stride];
      if(DEBUG>=2)
	for(int i = 0; i < 2*DELTA*M_stride; i++)
	  fmemblock[i] = aNaN;

      size_t fcnt = 0;
      for(int n = 1; n <= DELTA; n++){
	deltaX[n] = &fmemblock[fcnt-2];	fcnt += M_stride;
      }
      for(int n = 1; n <= DELTA; n++){
	deltaR[n] = &fmemblock[fcnt-2]; fcnt += M_stride;
      }
      if(DEBUG) assert(fcnt == 2ul * DELTA * M_stride);
    } else {
      deltaX = new PFLOAT*[2*(M+1)];
      deltaR = &deltaX[M+1];
      
      fmemblock = new PFLOAT[2 * DELTA * (M-1)];
      size_t fcnt = 0;
      for(int J = 2; J <= M; J++){
	deltaX[J] = &fmemblock[fcnt]; fcnt += DELTA;
	deltaR[J] = &fmemblock[fcnt]; fcnt += DELTA;
      }
      if(DEBUG) assert(fcnt == 2ul * DELTA * (M-1));
    }
  }

  /* if(NumScaleFactor > 1 && !(PAIRMERGE_OUTLIEREXTEND_FIX && PairMerge)){
    #pragma omp critical
    {
      if(NumScaleFactor > 1){
	printf("-ScaleDelta not implemented with pairwise alignment without -pairmerge: ignoring -ScaleDelta %0.3f %d\n",ScaleDeltaSize,ScaleDelta);
	fflush(stdout);
	
	ScaleDelta = 0;
	NumScaleFactor = 1;
      }
    }
    } */

  for(int scaleID = 0; scaleID < NumScaleFactor; scaleID++){

    if(Yid == Xid && scaleID)
      continue;

    FLOAT scale = ScaleFactor[scaleID];
    if(DEBUG && scaleID == 0) assert(scale == 1.0);

    if(scaleID > 0){
      for(int J = 0; J <= M+1; J++)
	Xrev[J] = X[J] * scale;

      if(((USE_SSE && USE_AVX) || USE_MIC) && USE_PFLOAT){// compute values in local arrays deltaX[1..DELTA][2..M] and deltaR[1..DELTA][2..M]
	for(int n = 1; n <= DELTA; n++){
	  PFLOAT *pdeltaXn = deltaX[n];
	  for(int J = n+1; J <= M; J++)
	    pdeltaXn[J] = (X[J] - X[J-n]) * scale;

	  PFLOAT *pdeltaRn = deltaR[n];
	  for(int J = n+1; J <= M; J++)
	    pdeltaRn[J] = (X[M+1-J+n] - X[M+1-J]) * scale;
	}
      } else {
	for(int J = 2; J <= M; J++){
	  PFLOAT *pdeltaXJ = deltaX[J];
	  FLOAT Xj = X[J];
	  FLOAT *pX = &X[J-DELTA];
	  int s = 0;
	  if(J < DELTA+1)
	    s = DELTA+1-J;
	  for(int n = s; n < DELTA; n++)
	    pdeltaXJ[n] = (Xj - pX[n]) * scale;
	  
	  PFLOAT *pdeltaRJ = deltaR[J];
	  Xj = X[M+1-J];
	  pX = &X[M+1+DELTA-J];
	  for(int n = s; s < DELTA; n++)
	    pdeltaRJ[n] = (pX[-n] - Xj) * scale;
	}
      }      
    }

    int myoutlierExtend = outlierExtend;
    int myoutlierBC = outlierBC;
    if(PoutlierEnd > 0.0){
      myoutlierExtend = 0;// Save time : First check if there is any alignment without outlierExtend (using endoutliers if needed)
      myoutlierBC = 0;
    }

    /* first check X in normal orientation and update best alignment */
    if((!phash || (!phash->orientation && (hashScaleDelta <= 1 || phash->scaleID == scaleID))) && Xmap->id != Ymap->id /* WAS430 Xid != Yid */)
      pairalignXY(Y,N,scaleID ? Xrev : X,M,A,align,Yid,Xid, Ymap->deltaY[0], Ymap->deltaX[0], scaleID ? deltaX : Xmap->deltaX[0], phash, 0, scaleID, myoutlierExtend, myoutlierBC,tid,Fmem,Mem,StrideMax);
  
    if(DEBUG>=2 && align->numpairs > 0){
      assert(align->sites1[0] >= 1 && align->sites1[align->numpairs-1] <= N);
      assert(align->sites2[0] >= 1 && align->sites2[align->numpairs-1] <= M);
    }

    if(VERB>=2){
      printf("Completed pairalignXY:Yid=%d(id=%lld),Xid=%d(id=%lld),or=0,sc=%d:N=%d,M=%d,tid=%d at wtime=%0.6f\n",  Yid,Ymap->id,Xid,Xmap->id,scaleID,N,M,tid,wtime());
      fflush(stdout);
    }

    if(DEBUG>=2 && !(align->mapid1 == Yid)){
#pragma omp critical
      {
	printf("pairalign:Yid=%d(id=%lld),Xid=%d(id=%lld):N=%d,M=%d\n", Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,N,M);
	printf("align->mapid1=%d,align->mapid2=%d\n",align->mapid1,align->mapid2);
	fflush(stdout);
	assert(align->mapid1 == Yid);
      }
    }

#ifdef VALGRIND
    VALGRIND_MAKE_MEM_UNDEFINED(Fmem,Fsiz*sizeof(PFLOAT));
    VALGRIND_MAKE_MEM_UNDEFINED(Imem,Isiz*sizeof(int));
#endif

    if(!(RepeatShift > 0.0 && align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= AlignedSiteThreshold)){
      /* next check X in reversed orientation and update best alignment */
      for(int J=0;J <= M+1;J++)
	Xrev[J] = (X[M+1]-X[M+1-J]) * scale;

      if(!phash || (phash->orientation && (hashScaleDelta <= 1 || phash->scaleID == scaleID)))
	pairalignXY(Y,N, Xrev, M, A, align, Yid, Xid, Ymap->deltaY[0], Ymap->deltaX[0], scaleID ? deltaR : Xmap->deltaR[0], phash, 1, scaleID, myoutlierExtend, myoutlierBC, tid, Fmem, Mem,StrideMax);
      else if (phash[0].id1==phash[1].id1 && phash[0].id2==phash[1].id2 && (HashBest <= 0 || phash[1].hashscore < max(besthashscore[Yid],besthashscore[Xid]) - HashBest)/* NEW111 */ 
									    && phash[1].orientation && (hashScaleDelta <= 1 || phash[1].scaleID == scaleID))
	pairalignXY(Y,N, Xrev, M, A, align, Yid, Xid, Ymap->deltaY[0], Ymap->deltaX[0], scaleID ? deltaR : Xmap->deltaR[0], &phash[1], 1, scaleID, myoutlierExtend, myoutlierBC, tid, Fmem, Mem,StrideMax);
    }

    if(!(PAIRMERGE_OUTLIEREXTEND_FIX && PairMerge) &&   // Else the following code has already been called recursively in pairalignXY before PairMerge test
       PoutlierEnd > 0.0 && outlierExtend && !myoutlierExtend && align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= AlignedSiteThreshold){// repeat alignment with outlierExtend
      if(DEBUG) assert(scaleID == 0);

      if(align->orientation)
	pairalignXY(Y,N,Xrev,M,A,align,Yid,Xid,Ymap->deltaY[0], Ymap->deltaX[0], Xmap->deltaR[0]/* WAS Ymap->deltaR[0] */, 0, 1, scaleID, outlierExtend, outlierBC,tid,Fmem,Mem,StrideMax);
      else
	pairalignXY(Y,N,X,M,A,align,Yid,Xid, Ymap->deltaY[0], Ymap->deltaX[0], Xmap->deltaX[0], 0, 0, scaleID, outlierExtend, outlierBC,tid,Fmem,Mem,StrideMax);
    }

    if(VERB>=2){
      printf("Completed pairalignXY:Yid=%d(id=%lld),Xid=%d(id=%lld),or=1,sc=%d:N=%d,M=%d,tid=%d at wtime=%0.6f\n",  Yid,Ymap->id,Xid,Xmap->id,scaleID,N,M,tid,wtime());
      fflush(stdout);
    }
  }// for scaleID = 0 .. 

  if(DEBUG>=2 && align->numpairs > 0){
    assert(align->sites1[0] >= 1 && align->sites1[align->numpairs-1] <= N);
    assert(align->sites2[0] >= 1 && align->sites2[align->numpairs-1] <= M);
  }

  if(DEBUG>=2 && !(align->mapid1 == Yid)){
    #pragma omp critical
    {
      printf("pairalign:Yid=%d(id=%lld),Xid=%d(id=%lld):N=%d,M=%d\n", Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,N,M);
      printf("align->mapid1=%d,align->mapid2=%d\n",align->mapid1,align->mapid2);
      fflush(stdout);
      assert(align->mapid1 == Yid);
    }
  }

  /* If PairSplit chimeric/local splits will be checked in caller */

  if(VERB>=2 /* || EVERB */){
    #pragma omp critical
    {
      if(align->numpairs >= 1)
	printf("pairalign:Yid=%d(id=%lld),Xid=%d(id=%lld):Score=%0.6f,logPV=%0.2f,RI=%d,RJ=%d,orientation=%d,aligned sites=%d,repeat=%d,offset=%0.3f\n",
	       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,align->score,align->logPV,align->sites1[align->numpairs-1],align->sites2[align->numpairs-1],
	       align->orientation,align->numpairs,align->repeat,
	       Y[align->sites1[align->numpairs-1]] - (align->orientation? Xrev : X)[align->sites2[align->numpairs-1]]);
      else
	printf("pairalign:Yid=%d(id=%lld),Xid=%d(id=%lld):orientation=%d,aligned sites=%d,repeat=%d\n",
	       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,align->orientation,align->numpairs,align->repeat);
      fflush(stdout);
    }
  }

  if(!PairSplit){
    if(align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= AlignedSiteThreshold){/* save local align[0] to galign[0] */
      if(DEBUG>=2 && !(align->mapid1 == Yid)){
	#pragma omp critical
	{
	  printf("pairalign:Yid=%d(id=%lld),Xid=%d(id=%lld):Score=%0.6f,logPV=%0.2f,RI=%d,RJ=%d,orientation=%d,aligned sites=%d,repeat=%d,offset=%0.3f\n",
		 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,align->score,align->logPV,align->sites1[align->numpairs-1],align->sites2[align->numpairs-1],
		 align->orientation,align->numpairs,align->repeat,
		 Y[align->sites1[align->numpairs-1]] - (align->orientation? Xrev : X)[align->sites2[align->numpairs-1]]);
	  printf("align->mapid1=%d,align->mapid2=%d\n",align->mapid1,align->mapid2);
	  fflush(stdout);
	  assert(align->mapid1 == Yid);
	}
      }
      if(DEBUG>=2) assert(align->mapid2 == Xid);

      if(!galign)
	galign = new Calign[1];
      /* swap *galign and *align : *align is on the stack and will automatically be deallocated by stack */
      Calign tmp = *galign;// NOTE : tmp ends up with copy of default initialized galign and can be safely deallocated by stack
      *galign = *align;
      *align = tmp;
      if(DEBUG>=2){
	assert(galign->numpairs > 0);
	assert(galign->sites1[0] >= 1 && galign->sites1[galign->numpairs-1] <= N);
	assert(galign->sites2[0] >= 1 && galign->sites2[galign->numpairs-1] <= M);
      }
      if(VERB>=2){
	printf("pairalign:Yid=%d,Xid=%d:numpairs=%d,I=%d..%d(N=%d),J=%d..%d(M=%d),sites1=%p,galign=%p (score=%0.3f,LogPV=%0.3f,repeat=%d)\n",
	       Yid,Xid,galign->numpairs,galign->sites1[0],galign->sites1[galign->numpairs-1],N,galign->sites2[0],galign->sites2[galign->numpairs-1],M,galign->sites1,galign,
	       galign->score,galign->logPV,galign->repeat);
        fflush(stdout);
      }
      if(DEBUG>=2) assert(galign->mapid1 == Yid);
      if(DEBUG>=2) assert(galign->mapid2 == Xid);
    } else {
      if(DEBUG && !EVERB && !(align->numpairs==0)){
	#pragma omp critical
	{
	  printf("pairalign:Yid=%d,Xid=%d:numpairs=%d,logPV=%0.2f,score=%0.6f\n",Yid,Xid,align->numpairs,align->logPV,align->score);
	  printf("  ScoreThreshold=%0.6f,LogPvThreshold=%0.2f,AlignedSiteThreshold=%d\n",ScoreThreshold,LogPvThreshold,AlignedSiteThreshold);
	  fflush(stdout);
	  assert(align->numpairs==0);
	}
      }
      //      align->allfree();
    }
  }

  if(DEBUG>=2 && galign && galign->numpairs > 0){
    if(VERB>=2){
      printf("Yid=%d,Xid=%d:numpairs=%d,I=%d..%d(N=%d),J=%d..%d(M=%d),sites1=%p,galign=%p (score=%0.3f,LogPV=%0.3f)\n",
	     Yid,Xid,galign->numpairs,galign->sites1[0],galign->sites1[galign->numpairs-1],N,galign->sites2[0],galign->sites2[galign->numpairs-1],M,galign->sites1,galign,
	     galign->score,galign->logPV);
      fflush(stdout);
    }
    assert(galign->sites1[0] >= 1 && galign->sites1[galign->numpairs-1] <= N);
    assert(galign->sites2[0] >= 1 && galign->sites2[galign->numpairs-1] <= M);
  }

  if(M > MAX_ALLOCA)
    delete [] Xrev;

  if(NumScaleFactor > 1){/* free up deltaX,deltaR,fmemblock */
    delete [] fmemblock;
    delete [] deltaX;
  }

  if(RA_MIN_TIME){/* release memory */
    if(DEBUG) assert(Fmem != NULL && Imem != NULL);
    if(DEBUG) assert(Mem->MemSiz > 0);
    
    if(!MEMCHECK || RA_MIN_MEM >= 8)
      ReleaseMemory(Mem,nMemSiz,tid,Yid,Xid,N,M,Fmem);

    if(RA_MIN_MEM >= 8){// use munmap() to free memory 
      if(Fmem && munmap(Fmem, Mem->MemSiz)){
	int eno = errno;
	char *err = strerror(eno);
	printf("munmap(%p,%lld) failed (MemUsed=%lld): errno=%d:%s\n",Fmem, Mem->MemSiz, Mem->MemUsed, eno, err);
	dumpmemmap();
	fflush(stdout);exit(1);
      }
      Fmem = NULL;
      Imem = NULL;
      Mem->MemSiz = 0;
      Mem->MemUsed = 0;
    }
  }
}

static void maxhashalloc(size_t num, size_t numhashpairs, size_t &maxhashpairs, Chashpair * &hashpairs)
{
  size_t newmax = num;
  if(newmax <= maxhashpairs)
    return;
  if(maxhashpairs > 0 && newmax < maxhashpairs * 2)
    newmax = maxhashpairs * 2;

  if(VERB && HASH_DEBUG){
    printf("maxhashalloc:&hashpairs=%p,hashpairs=%p,maxhashpairs=%llu->%llu, size=%llu bytes\n",&hashpairs, hashpairs,(unsigned long long)maxhashpairs,(unsigned long long)newmax,(unsigned long long)newmax*sizeof(Chashpair));
    fflush(stdout);
  }

  if(maxhashpairs <= 0)
    hashpairs = new Chashpair[newmax+1];/* always keep space for one extra entry in case the last entry is not part of a matched pair (both orientations) */
  else {
    Chashpair *orighashpairs = hashpairs;
    hashpairs = new Chashpair[newmax+1];
    for(size_t i = 0; i <= numhashpairs; i++)
      hashpairs[i] = orighashpairs[i];
    delete [] orighashpairs;
  }
  maxhashpairs = newmax;
  if(VERB && HASH_DEBUG){
    printf("maxhashalloc:hashpairs=%p,maxhashpairs=%llu: cum wall time=%0.6f secs\n",hashpairs,(unsigned long long)maxhashpairs, wtime());
    fflush(stdout);
  }
}

static int ChashIdInc(Chashpair *p1, Chashpair *p2)
{
  return (p1->id1 > p2->id1) ? 1 : (p1->id1 < p2->id1) ? -1 : 
    (p1->id2 > p2->id2) ? 1 : (p1->id2 < p2->id2) ? -1 : 
    (p1->orientation > p2->orientation) ? 1 : (p1->orientation < p2->orientation) ? -1 : 0;
}

static int gnumthreads = 1;

#define MAXALIGN_THRESHOLD 0 // If > 0 : report largest alignments per map if over this number */
static int maxaligns_permap = 0;

/** stream in enough hashtable entries to fill the entire table hashpairs[0..maxhashpairs-1] (before filtering by hashscore).
    Then copy hashtable entries to YidList[],XidList[],pHashList[] until we have MaxPair Entries (or we reach end of hashtable) and return number of entries.

    First check if the 2nd buffer has data : if so just copy from 2nd buffer instead of reading from disk */
static void loadhash(Chashpair* &nexthash, Chashpair* &hashpairs, size_t &numhashpairs, size_t maxhashpairs, int *YidList, int *XidList, Chashpair **phashList, int &NumPair)
{
  /* make sure any existing entries are at start of hashpairs[] : typically there will be at most one remaining entry */
  Chashpair *p = hashpairs;
  Chashpair *end = &hashpairs[numhashpairs];
  if(nexthash > p) { /* copy existing entries to start of hashpairs[] */
    if(VERB && HASH_DEBUG && end > nexthash){
      printf("copying %llu entries from hashpairs%s[%llu..%llu] to hashpairs%s[0..%llu]: cum cpu time=%0.6f, wall time=%0.6f secs\n",
	     (unsigned long long)(end-nexthash), hashpairs==hashpairs2 ? "2" : "", (unsigned long long)(nexthash-hashpairs), (unsigned long long)(end-hashpairs)-1, hashpairs==hashpairs2 ? "2" : "", 
	     (unsigned long long)(end-nexthash)-1, mtime(), wtime());
      fflush(stdout);
    }

    int len = end-nexthash;

#ifndef VALGRIND
    if((USE_MIC && USE_PFLOAT) ? (nexthash - hashpairs)*sizeof(Chashpair) >= 512 : (nexthash-hashpairs) >= len)
      memcpy(hashpairs, nexthash, len*sizeof(Chashpair));
    else
#endif
      memmove(hashpairs, nexthash, len*sizeof(Chashpair));

    p = &hashpairs[end-nexthash];
    p->id1 = 0;
  } else
    p = end;

  size_t cur = p - hashpairs;/* entries left */
  size_t num = maxhashpairs - cur;/* space left */
  nexthash = hashpairs;
  numhashpairs = cur;

  if(hashpairs != hashpairs2){/* Try to swap 2nd buffer and first buffer */
    if(DEBUG) assert(hashpairs == hashpairs1);
    Chashpair *end2 = &hashpairs2[numhashpairs2];
    if(nexthash2 < end2 || NumPair2 > 0){
      if(VERB && HASH_DEBUG){
	printf("Swapping Buffers: NumPair1=%d,Matches1=%llu, NumPair2=%d, Matches2=%llu: cum CPU time=%0.6f, wall time=%0.6f\n", 
	       NumPair1, (unsigned long long)(&hashpairs1[numhashpairs1] - nexthash1), NumPair2, (unsigned long long)(end2-nexthash2), mtime(), wtime());
	fflush(stdout);
      }
      if(DEBUG) assert(nexthash1 == &hashpairs1[numhashpairs1]);
      if(DEBUG) assert(cur == 0);
      if(DEBUG) assert(YidList == YidList1 && XidList == XidList1 && phashList == phashList1);
      if(DEBUG) assert(NumPair1 == 0);

      /* swap the two sets of buffers */
      Chashpair *hashtmp = hashpairs1; hashpairs1 = hashpairs2; hashpairs2 = hashtmp;
      hashtmp = nexthash1; nexthash1 = nexthash2; nexthash2 = hashtmp;
      if(DEBUG) assert(maxhashpairs1 == maxhashpairs2);
      size_t tmpSize = numhashpairs1; numhashpairs1 = numhashpairs2; numhashpairs2 = tmpSize;

      int *tmpList = YidList1; YidList1 = YidList2; YidList2 = tmpList;
      tmpList = XidList1; XidList1 = XidList2; XidList2 = tmpList;
      Chashpair **phashtmp = phashList1; phashList1 = phashList2; phashList2 = phashtmp;
      int tmpPair = NumPair1;NumPair1 = NumPair2; NumPair2 = tmpPair;
      
      /* copy any remaining data in hashpairs1[] to hashpairs2[] and reset hashpairs1[] */
      if(DEBUG) assert(nexthash2 == &hashpairs2[numhashpairs2]);
      nexthash2 = &hashpairs2[numhashpairs2 = 0];
      if(DEBUG) assert(NumPair2 == 0);

      Chashpair *end1 = &hashpairs1[numhashpairs1];
      while(nexthash1 < end1)
	*nexthash2++ = *nexthash1++;
      nexthash2->id1 = 0;

      numhashpairs2 = nexthash2 - hashpairs2;

      return;
    }
  }

  int HashEof = hash_eof(HashFP);
  if(HashEof && cur==0)
    return;

  if(!HashEof){
    int linecnt = 0;/* Not Used for Binary File input */
    if(VERB && HASH_DEBUG){
      printf("cur=%llu,num=%llu,maxhashpairs=%llu:Calling hash_read(HashFP,0,&hashpairs[cur],num,linecnt): cum CPU time=%0.6f, wall time=%0.6f secs\n",
	     (unsigned long long)cur,(unsigned long long)num,(unsigned long long)maxhashpairs, mtime(), wtime());
      fflush(stdout);
    }
    if(DEBUG) assert(p == &hashpairs[cur]);
    numhashpairs += hash_read(HashFP, 0, p, num, linecnt);
    hashpairs[numhashpairs].id1 = 0;
    HashEof = hash_eof(HashFP);
    if(VERB && HASH_DEBUG && numhashpairs > cur){
      printf("Read in next %llu Matches (total=%llu/%llu) from %s (EOF=%d) into hashpairs%d[]: cum cpu time=%0.6f, wall time=%0.6f secs\n", 
	     (unsigned long long)(numhashpairs-cur), (unsigned long long)(totalhashpairs + numhashpairs - cur), (unsigned long long)MatchMax, 
	     hash_filename, HashEof, hashpairs!=hashpairs2 ? 1 : 2, mtime(), wtime());
      fflush(stdout);
    }
    if(DEBUG>=2){    /* check if hash pairs are in ascending order of id1,id2,orientation */
      size_t i;
      for(i = 1; i < numhashpairs; i++)
	if(ChashIdInc(&hashpairs[i-1],&hashpairs[i]) > 0)
	  break;
      if(i < numhashpairs){
	printf("Hashtable must be sorted in ascending order of id1,id2,orientation\n");
	for(size_t j = i-1; j <= i; j++)
	  printf("hash[%llu]:id1=%lld,id2=%lld,orientation=%d\n",(unsigned long long)(totalhashpairs + j),(long long)hashpairs[j].id1,(long long)hashpairs[j].id2,hashpairs[j].orientation);
	fflush(stdout);exit(1);
      }
    }
  }
	  
  /* remove entries with score < hash_threshold */
  size_t j = cur;
  for(size_t i = cur; i < numhashpairs; i++){
    if(hashpairs[i].hashscore < hash_threshold)
      continue;
    hashpairs[j++] = hashpairs[i];
  }
  totalhashpairs += numhashpairs - cur;
  if(VERB>=2){
    printf("After filtering %lu matches for score >= %d : remaining %lu matches (+ %lu previous matches),nexthash= %lu,NumPair=%d,MaxPair=%d,maxhashpairs=%lu\n", 
	   numhashpairs-cur, hash_threshold, j-cur, cur, nexthash-hashpairs,NumPair,MaxPair,maxhashpairs);

    fflush(stdout);
  }
  numhashpairs = j;
  hashpairs[numhashpairs].id1 = 0;

  /* copy hashtable entries to YidList[], XidList[], phashList[] */
  if(DEBUG && hashpairs == hashpairs1) assert(NumPair1 == NumPair);
  if(DEBUG && hashpairs == hashpairs2) assert(NumPair2 == NumPair);

  /* just copy hashtable entries to YidList[],XidList[],phashList[] until we reach end of hashtable buffer */
  end = &hashpairs[numhashpairs];

  int lastYid = -1, lastXid = -1;
  int TNumPair = NumPair;
  int align_cnt = 0;

  for(; NumPair < MaxPair;){
    if(VERB>=3){
      printf("nexthash= %lu, end= %lu, NumPair= %d, HashEof= %d:\n",nexthash-hashpairs,end-hashpairs,NumPair1,HashEof);
      fflush(stdout);
    }

    if(nexthash >= end-1){/* next to end of current read buffer : check if this is EOF */
      if(!HashEof)
	break;/* reuse last entry in next Yid iteration since it may come in a pair of entries (both orientations) that are needed together */
      /* EOF : force break out of Yid loop, unless the last entry still needs to be processed */
      if(nexthash >= end)
	break;
    }

    if(DEBUG>=2) assert(0 <= nexthash->id1 && nexthash->id1 < numY && YYmap[nexthash->id1]->id >= 0);
    if(DEBUG>=2) assert(0 <= nexthash->id2 && nexthash->id2 < numX && XXmap[nexthash->id2]->id >= 0);

    if(MAXALIGN_THRESHOLD){
      if(nexthash->id1 != lastYid){
	if(lastYid != -1){
	  if(align_cnt > maxaligns_permap){
	    maxaligns_permap = align_cnt;
	    if(align_cnt >= MAXALIGN_THRESHOLD){
	      printf("Yid= %d: align_cnt= %d, NumPair=%d\n",lastYid,align_cnt,NumPair);
	      for(int i = 0; i < align_cnt; i++){
		Chashpair *phash = phashList[i + NumPair - align_cnt];
		if(HashBest > 0)
		  printf("\t i=%d: Yid=%d(id=%lld),Xid=%d(id=%lld), or=%d: hashscore= %d, offset= %d kb, bestscore= %d,%d\n",
			 i, phash->id1,YYmap[phash->id1]->id,phash->id2,XXmap[phash->id2]->id,phash->orientation,phash->hashscore,phash->offset,besthashscore[phash->id1],besthashscore[phash->id2]);
		else
		  printf("\t i=%d: Yid=%d(id=%lld),Xid=%d(id=%lld), or=%d: hashscore= %d, offset= %d kb\n",
			 i, phash->id1,YYmap[phash->id1]->id,phash->id2,XXmap[phash->id2]->id,phash->orientation,phash->hashscore,phash->offset);
	      }
	      fflush(stdout);
	    }
	  }
	}
	align_cnt = 0;
      }
    }

    if(HashBest > 0){ // NEW111
      if(nexthash->hashscore < max(besthashscore[nexthash->id2],besthashscore[nexthash->id1]) - HashBest){
	lastYid = nexthash->id1;
	lastXid = nexthash->id2;

	nexthash++;
	continue;
      }
    }

    if(MAXALIGN_THRESHOLD)
      align_cnt++;

    YidList[NumPair] = lastYid = nexthash->id1;
    XidList[NumPair] = lastXid = nexthash->id2;

    phashList[NumPair] = nexthash++;
    TNumPair = ++NumPair;

    /* skip next entry in hashtable if it is the same pair of maps as lastYid,lastXid :
       This happens when both orientation are in the HashTable */
    if(nexthash < end && nexthash->id1 == lastYid && nexthash->id2 == lastXid)
      nexthash++;
  }
  if(DEBUG) assert(HashBest > 0 ? NumPair >= TNumPair : NumPair == TNumPair);

  if(VERB && HASH_DEBUG && HASH_STREAM>=2){
    printf("Converted into %d map pairs (%llu Matches remaining):cum CPU time=%0.6f,wall time=%0.6f\n",
	   NumPair, (unsigned long long)(&hashpairs[numhashpairs] - nexthash), mtime(), wtime());
    fflush(stdout);
  }

  if(hashpairs == hashpairs1){
    /* copy any remaining data in hashpairs1[] to hashpairs2[] and reset hashpairs1[] */
    if(DEBUG) assert(nexthash2 == &hashpairs2[numhashpairs2]);
    nexthash2 = &hashpairs2[numhashpairs2 = 0];
    if(DEBUG) assert(NumPair2 == 0);

    Chashpair *end1 = &hashpairs1[numhashpairs1];
    while(nexthash1 < end1)
      *nexthash2++ = *nexthash1++;
    nexthash2->id1 = 0;
    numhashpairs2 = nexthash2 - hashpairs2;
  }
}

#if 0
/* scan ahead in sorted array to find specified id pair : return pointer to Chashpair if found, 0 if not */
static Chashpair *findhash(Chashpair * &nexthash1, int id1, int id2)
{
  if(DEBUG) assert(!HASH_STREAM && !HashFP);

  Chashpair *phash = nexthash1, *end = &hashpairs1[numhashpairs1];
  int found = 0;
  for(; phash < end; phash++){
    if(phash->id1 < id1)
      continue;
    if(phash->id1 > id1)
      break;
    if(phash->id2 < id2)
      continue;
    found =  (phash->id2 > id2) ? 0 : 1;
    break;
  }

  nexthash1 = phash;

  return found ? phash : 0;
}
#endif

/* sort alignments with Yid==Xid before those with Yid != Xid, and then
   IF PairMergeBigMapsFirst : sort in descending order of size of large map, and then
   sort alignments in descending order of logPV */
static int PVdec2(register Calign **pp1, register Calign **pp2)
{
  register Calign *p1 = *pp1, *p2 = *pp2;

  int same1 = (p1->mapid1 == p1->mapid2) ? 1 : 0;
  int same2 = (p2->mapid1 == p2->mapid2) ? 1 : 0;
  if (same1 < same2) 
    return 1;
  if (same1 > same2)
    return -1;
  if (PairMergeBigMapsFirst){
    Cmap *Ymap1 = YYmap[p1->mapid1];
    Cmap *Ymap2 = YYmap[p2->mapid1];
    FLOAT *Y1 = Ymap1->site[0];
    FLOAT *Y2 = Ymap2->site[0];
    int N1 = Ymap1->numsite[0];
    int N2 = Ymap2->numsite[0];

    Cmap *Xmap1 = XXmap[p1->mapid2];
    Cmap *Xmap2 = XXmap[p2->mapid2];
    FLOAT *X1 = Xmap1->site[0];
    FLOAT *X2 = Xmap2->site[0];
    int M1 = Xmap1->numsite[0];
    int M2 = Xmap2->numsite[0];

    FLOAT BigLen1 = max(Y1[N1+1],X1[M1+1]);
    FLOAT BigLen2 = max(Y2[N2+1],X2[M2+1]);

    if(BigLen1 < BigLen2) 
      return 1;
    if(BigLen1 > BigLen2)
      return -1;
  }
  

  return (p1->logPV < p2->logPV) ? 1 : (p1->logPV > p2->logPV) ? -1 : 0;
}

/* try to merge maps in Gmap[0 .. nummaps-1] based on alignment information in alignment[0..numaligns-1] and append all merged or modified maps to Gmap[] */
/* If PairMergeHmap > 0, verify that no more maps can be merged or modified. Instead try to combine maps to form .hmaps :
        1. Look for map pairs where one map completely overlaps another but there is an outlier larger than PairMergeMaxOutlier. Treat the internal outlier as a Het Indel 
	   For end overlaps, defer merges, but mark each end of the map with an end overlap (with pointer to alignment)
	2. Look for Two map pairs with endoutliers larger than PairMergeMaxEndKB which share a common larger map and where the smaller maps are fully overlapped by the larger maps
	       and lie close to each other on the larger map : Assume the two form an incomplete Het Insertion. (NOT YET IMPLEMENTED)
        3. After all merges from previous steps, repeatedly loop over remaining maps and for each map with only one end overlap :
	    A. Merge it with the other map, leaving the other map only (in case the other map has 2 end overlaps) : See step 1 above, except don't select the larger map to keep.
	    B. Mark remaining other map as no longer having the end overlap that was just completed.
   Returns number of new maps created by splitting (see -SplitSegDup)
 */
int pairmerge(int pairmergeIter)
{
  int splitcnt = 0;

  int *scaleIDcnt = (int *)alloca(NumScaleFactor * sizeof(int));
  for(int k = 0; k < NumScaleFactor; k++)
    scaleIDcnt[k] = 0;

  /* compact alignment array */
  size_t aligncnt = 0;
  for(size_t i=0; i < numaligns; i++){
    Calign *p = alignment[i];
    if(!p || p->numpairs <= 0){
      if(VERB>=3){
	printf("Alignment %lu/%lu: p= %p, numpairs= %d\n",i,numaligns,p,p ? p->numpairs : 0);
	fflush(stdout);
      }
      continue;
    }
    int Yid = p->mapid1, Xid = p->mapid2;
    if(p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold || (PAIRMERGE_REPEATMASK && p->repeat > 0)){
      if(EVERB){
	printf("Alignment %lu/%lu: Yid=%d,id=%lld, Xid=%d,id=%lld, or=%d,sc=%d:numpairs=%d,score=%0.6f,logPV=%0.2f,repeat=%d: discarding due to -S %0.3f -T %0.2f -A %d (or repeat > 0)\n",
	       i,numaligns,Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,p->orientation,(NumScaleFactor > 1) ? p->scaleID : 0, 
	       p->numpairs,p->score,p->logPV,p->repeat,ScoreThreshold,LogPvThreshold,AlignedSiteThreshold);
	fflush(stdout);
      }
      continue;
    }

    scaleIDcnt[p->scaleID]++;

    if(EVERB>=3){
      printf("Alignment %lu/%lu: Yid=%d,id=%lld, Xid=%d,id=%lld, or=%d:numpairs=%d,score=%0.6f,logPV=%0.2f,repeat=%d: moving to alignment[%lu]\n",
	     i,numaligns,Yid,YYmap[Yid]->id, Xid,XXmap[Xid]->id, p->orientation,p->numpairs,p->score,p->logPV,p->repeat,aligncnt);
      fflush(stdout);
    }

    if(aligncnt < i){
      Calign *tmp = alignment[aligncnt];
      alignment[aligncnt] = p;
      alignment[i] = tmp;
    }
    aligncnt++;
  }
  if(VERB/* HERE >=2 */){
    printf("pairmerge: Removed below threshold alignments : numaligns= %lu -> %lu\n",numaligns, aligncnt);
    if(NumScaleFactor > 1)
      for(int i = 0; i < NumScaleFactor; i++)
	printf("scaleID=%d,scale=%0.4f: aligncnt=%d\n", i, ScaleFactor[i], scaleIDcnt[i]);
    fflush(stdout);
  }
  numaligns = aligncnt;

  /* sort all alignments with mapid1==mapid2 before those with mapid1 != mapid2,
     otherwise, IF PairMergeBigMapsFirst > 0 : in descending order of size of large map,
     otherwise in descending order of logPV */
  qsort(alignment,numaligns,sizeof(Calign *),(intcmp*) PVdec2);

  double Xscale = 1.0;
  double Yscale = 1.0;
  if(!(spots_filename[0] || (svcheck && numrefmaps > 0)) && fabs(PixelLen - origPixelLen) > 1e-12 /* NEW */){
    Yscale = Xscale = origPixelLen/PixelLen;
    if(PairSplit)/* -i inputs that correspond to reference is not scaled by bpp/500 */
      Yscale = 1.0;
  }

  if(colors==1){

    if(VERB/* HERE >=2 */){
      printf("pairmerge: numaligns=%lu\n",numaligns);
      fflush(stdout);
    }

    for(size_t i=0; i < numaligns; i++){
      Calign *p = alignment[i];

      int K = p->numpairs;
      int Yid = p->mapid1, Xid = p->mapid2;
      Cmap *Ymap = YYmap[Yid];
      Cmap *Xmap = XXmap[Xid];
      int N = Ymap->numsite[0];
      int M = Xmap->numsite[0];
      FLOAT *Y = Ymap->site[0];
      FLOAT *X = Xmap->site[0];

      if(DEBUG>=2 && CmapChimQuality <= 2){
	assert(Xmap->FragSd[0] == NULL);
	assert(Xmap->ExpSd[0] == NULL);
	assert(Xmap->FragCov[0] == NULL);
	assert(Xmap->FragChiSq[0] == NULL);

	assert(Ymap->FragSd[0] == NULL);
	assert(Ymap->ExpSd[0] == NULL);
	assert(Ymap->FragCov[0] == NULL);
	assert(Ymap->FragChiSq[0] == NULL);
      }

      if(MaxContigSiteDensity > 0.0){
	if(N > 1 && (N-1)*100.0 > MaxContigSiteDensity * (Y[N] - Y[1])){
	  if(VERB){
	    printf("Skipping Alignment %llu/%llu:Ymap->id=%lld, Site Density = %0.2f / 100kb (MaxContigSiteDensity= %0.1f)\n",
		   (unsigned long long)i, (unsigned long long)numaligns, Ymap->id,(N-1)*100.0/(Y[N]-Y[1]),MaxContigSiteDensity);
	    fflush(stdout);
	  }
	  continue;
	}
	if(M > 1 && (M-1)*100.0 > MaxContigSiteDensity * (X[M] - X[1])){
	  if(VERB){
	    printf("Skipping Alignment %llu/%llu:Ymap->id=%lld, Site Density = %0.2f / 100kb (MaxContigSiteDensity= %0.1f)\n",
		   (unsigned long long)i, (unsigned long long)numaligns, Xmap->id,(M-1)*100.0/(X[M]-X[1]), MaxContigSiteDensity);
	    fflush(stdout);
	  }
	  continue;
	}
      }

      if(DEBUG) assert(Ymap->origmap == NULL);
      if(DEBUG) assert(Xmap->origmap == NULL);      

      int LI = p->sites1[0];
      int LJ = p->sites2[0];
      int RI = p->sites1[K-1];
      int RJ = p->sites2[K-1];
      if(DEBUG/* HERE >=2 */  && !(0 < LI && LI <= RI && RI <= N && 0 < LJ && LJ <= RJ && RJ <= M)){
	if(PairMergeRepeat)
	  printf("%d:",pairmergeIter);
	printf("Alignment %llu/%llu:Ymap=%d,id=%lld(paired=%d),Xmap=%d,id=%lld(paired=%d):score=%0.6f,logPV=%0.2f,or=%d,sc=%d,pairs=%d,Y[LI=%d]=%0.3f,Y[RI=%d]=%0.3f,Y[N=%d]=%0.3f,YL=0x%lx,YR=0x%lx\n\t X[LJ=%d]=%0.3f,X[RJ=%d]=%0.3f,X[M=%d]=%0.3f,XL=0x%lx,XR=0x%lx\n",
	       (unsigned long long)i,(unsigned long long)numaligns, Ymap->mapid, Ymap->id, Ymap->paired, Xmap->mapid, Xmap->id,Xmap->paired,p->score,p->logPV,p->orientation, p->scaleID, p->numpairs,
	       LI,Y[LI],RI,Y[RI],  N, Y[N],Ymap->Mask[0] ? Ymap->Mask[0][1] : 0, Ymap->Mask[0] ? Ymap->Mask[0][N+1] : 0,
	       LJ, p->orientation ? (X[M+1] - X[M+1 - LJ]) : X[LJ], RJ, p->orientation ? (X[M+1] - X[M+1 - RJ]) : X[RJ], Xmap->numsite[0], X[Xmap->numsite[0]],
	       Xmap->Mask[0] ? Xmap->Mask[0][p->orientation ? M+1 : 1] : 0, Xmap->Mask[0] ? Xmap->Mask[0][p->orientation ? 1 : M+1] : 0);
	printf("\t Y[N+1]=%0.3f, X[M+1]=%0.3f\n",Y[N+1], X[Xmap->numsite[0]+1]);
	fflush(stdout);

	assert(0 < LI && LI <= RI && RI <= N);
	assert(0 < LJ && LJ <= RJ && RJ <= M);
      }

      /* average center to center distance over all aligned sites */
      double distanceYtoX = 0;
      /* adjust location to be relative to origmap */
      double centerY = Y[N+1]*0.5;
      double centerX = X[M+1]*0.5;
      for(int index = 0; index < K; index++){
	int I = p->sites1[index];
	int J = p->sites2[index];
	if(VERB>=3 && Ymap->id==3 && Xmap->id==9){
	  printf("Yid=%lld,Xid=%lld,or=%d,N=%d(%d),M=%d(%d),centerX=%0.4f,centerY=%0.4f,index=%d/%d,I=%d,J=%d: YtoX= %0.4f\n",
		 Ymap->id,Xmap->id,p->orientation,N,N,M,M,centerX,centerY,index,K,I,J,
		 (Y[I]-centerY) + (centerX - (p->orientation ? (X[M+1]-X[M+1-J]) : X[J])));
	  fflush(stdout);
	}
	distanceYtoX += (Y[I]-centerY) + (centerX - (p->orientation ? (X[M+1]-X[M+1-J]) : X[J]));
      }
      distanceYtoX /= K;

      double overlap = (Y[RI]-Y[LI] + (p->orientation ? X[M+1-LJ] - X[M+1-RJ] : X[RJ]-X[LJ]))*0.5;
      double leftEnd = (p->Lend >= -1) ? 0.0 : min(Y[LI], (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]));
      int leftEndN = (p->Lend >= -1) ? 0 : (Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ])) ? LI-1 : LJ-1;
      double rightEnd = (p->Rend >= -1) ? 0.0 : min(Y[N+1]-Y[RI], (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]));
      int rightEndN = (p->Rend >= -1) ? 0 : (Y[N+1]-Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ])) ? N-RI : M-RJ;
      double maxOutlier = 0.0;
      int maxOutlierM = 0;
      int IleftOutlier = 0;
      int IrightOutlier = 0;
      int JleftOutlier = 0;
      int JrightOutlier = 0;

      int Ileft = N+1;/* Ref Label of left end of leftmost outlier over PairMergeMaxOutlier (WAS40 MaxInvPalindromeSL) */
      int Iright = 0;/* Ref Label of right end of rightmost outlier over PairMergeMaxOutlier (WAS40 MaxInvPalindromeSL) */

      if(PairMergeMaxOutlier >= 0.0 || PairMergeExcludeOutliers){/* locate largest outlier (and subtract outliers from overlap) */
	int lastI = LI, I;
	int lastJ = LJ, J;
	for(int k = 1; k < K; lastI = I, lastJ = J, k++){
	  I = p->sites1[k];
	  J = p->sites2[k];
	  
	  /* outlier lies between sites lastI .. I on reference and lastJ .. J on query */
	  if(DEBUG>=2) assert(J >= lastJ);
	  if(DEBUG>=2) assert(I >= lastI);
	  double qrystartpos = (p->orientation ? (X[M + 1] - X[M + 1 - lastJ]) : X[lastJ]) * Xscale; /* QryStartPos */
	  double qryendpos   = (p->orientation ? (X[M + 1] - X[M + 1 - J]) : X[J])*Xscale; /* QryEndPos */
	  double refstartpos = Y[lastI]*Yscale; /* RefStartPos */
	  double refendpos   = Y[I]*Yscale; /* RefEndPos */
	  if(DEBUG>=2) assert(qryendpos >= qrystartpos);// If this fails, fix next line
	  double Outlier = fabs(refendpos - refstartpos - qryendpos + qrystartpos);// indel size
	  int OutlierM = J-lastJ + I-lastI - 2;// misaligned labels

	  int outlier = (p->iscore[k] > p->outscore[k] + (RFLOAT)0.01 || (outlierExtend && (lastJ >= J || lastI >= I)) || Outlier > PairMergeMaxOutlier || OutlierM > PairMergeMaxOutlierM) ? 1 : 0;

	  if(EVERB || (DEBUG>=2 && !(p->iscore[k] >= p->outscore[k] - (RFLOAT)0.01))){
	    printf("    Ymap->id=%lld,Xmap->id=%lld,or=%d,N=%d,M=%d,Xscale=%0.6f,Yscale=%0.6f:I=%d..%d(%0.3f..%0.3f),J=%d..%d(%0.3f..%0.3f),Outlier=%0.3f,%d,outlier=%d:k=%d/%d,p->iscore[k]=%0.6f,outscore[k]=%0.6f\n",
		   Ymap->id,Xmap->id,p->orientation,N,M,Xscale,Yscale,lastI,I,refstartpos,refendpos,lastJ,J,qrystartpos,qryendpos,Outlier,OutlierM,outlier,k,K,p->iscore[k],p->outscore[k]);
	    if(CmapChimQuality >= 3){
	      int jmin = p->orientation ? M+1-J : lastJ;
	      int jmax = p->orientation ? M+1-lastJ : J;
	      printf("\t Ymap:MolSd=%0.1f..%0.1f,ExpSd=%0.1f..%0.1f,MolCov=%0.2f..%0.2f,ChiSq=%0.3e..%0.3e\n\t Xmap:MolSd=%0.1f..%0.1f,ExpSd=%0.1f..%0.1f,MolCov=%0.2f..%0.2f,ChiSq=%0.3e..%0.3e\n",
		     Ymap->FragSd[0][lastI],Ymap->FragSd[0][I-1],Ymap->ExpSd[0][lastI],Ymap->ExpSd[0][I-1],Ymap->FragCov[0][lastI],Ymap->FragCov[0][I-1],Ymap->FragChiSq[0][lastI],Ymap->FragChiSq[0][I-1],
		     Xmap->FragSd[0][jmin],Xmap->FragSd[0][jmax-1],Xmap->ExpSd[0][jmin],Xmap->ExpSd[0][jmin-1],
		     Xmap->FragCov[0][jmin],Xmap->FragCov[0][jmax-1],Xmap->FragChiSq[0][jmin],Xmap->FragChiSq[0][jmax-1]);
	    }
	    fflush(stdout);
	    if(DEBUG>=2) assert(p->iscore[k] >= p->outscore[k] - (RFLOAT)0.01);
	  }

	  if(!outlier)/* not an outlier */
	    continue;
	  
	  double origoverlap  = overlap;

	  if(PairMergeExcludeOutliers)
	    overlap -= (refendpos - refstartpos + fabs(qryendpos - qrystartpos)) * 0.5;

	  if(PairMergeMaxOutlier < 0)
	    continue;

	  if(PairMergeOutlierChiSq > 0.0 && CmapChimQuality >= 3){// skip if outlier variance is too large or has too many molecule outliers
	    if(refendpos - refstartpos < qryendpos - qrystartpos){
	      float *OutlierFrac = Ymap->OutlierFrac[0];
	      float *FragSd = Ymap->FragSd[0];
	      float *ExpSd = Ymap->ExpSd[0];
	      float *FragCov = Ymap->FragCov[0];
	      float *FragChiSq = Ymap->FragChiSq[0];
	      if(DEBUG>=2) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);
	      
	      float *ChimQuality = Ymap->ChimQuality[0];
	      float *ChimNorm = Ymap->ChimNorm[0];

	      double LRange = (p->orientation ? X[M+1-lastJ] - X[M+1-J] : X[J] - X[lastJ]) * Xscale;
	      double LRangeRatio = LRange / AvgInterval;

	      int i = lastI;
	      for(; i < I; i++){
		double SRange = (Y[i+1] - Y[i]) * Yscale;
		double FragCovi = FragCov[i];
		if(PairMergeChimNorm > 0.0 && ChimNorm[i] > 0.0 && ChimNorm[i+1] > 0.0)
		  FragCovi = min(FragCovi, min(ChimNorm[i]*ChimQuality[i],ChimNorm[i+1]*ChimQuality[i+1]) * PairMergeChimNorm * 0.01);

		if(Outlier < PairMergeOutlierMaxSize &&
		   (OutlierFrac[i] > PairMergeOutlierFrac || 
		    (FragCovi < PairMergeOutlierMaxCov && 
		     LRange > PairMergeOutlierMinInterval &&
		     ((FragChiSq[i] < PairMergeOutlierChiSq &&
		       FragSd[i] > ExpSd[i] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
		       FragSd[i] > ExpSd[i] * PairMergeOutlierSdRatio) ||		     
		      FragCovi < PairMergeOutlierMinCov + 
		      ((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != abs(J-lastJ))?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale)))){
		  if(EVERB){
		    printf("\ti=%d:I=%d..%d:LRange= %0.4f, Ratio= %0.4f, SRange= %0.4f, OutlierFrac= %0.3f, FragChiSq= %0.4e, Outlier= %0.4f, FragCov= %0.2f(N1=%0.2f,%0.2f), FragSd= %0.1f,ExpSd= %0.1f : skipping outlier\n",
			   i,lastI,I,LRange, LRangeRatio, SRange, OutlierFrac[i], FragChiSq[i], Outlier, FragCov[i], ChimNorm[i]*ChimQuality[i]*0.01,ChimNorm[i+1]*ChimQuality[i+1]*0.01, FragSd[i], ExpSd[i]);
		    fflush(stdout);
		  }
		  break;
		}
	      }
	      if(i < I)
		continue;// outlier will be ignored
	    } else {
	      float *OutlierFrac = Xmap->OutlierFrac[0];
	      float *FragSd = Xmap->FragSd[0];
	      float *ExpSd = Xmap->ExpSd[0];
	      float *FragCov = Xmap->FragCov[0];
	      float *FragChiSq = Xmap->FragChiSq[0];
	      if(DEBUG>=2) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);

	      float *ChimQuality = Xmap->ChimQuality[0];
	      float *ChimNorm = Xmap->ChimNorm[0];

	      double LRange = (Y[I] - Y[lastI]) * Yscale;
	      double LRangeRatio = LRange / AvgInterval;

	      int jmin = p->orientation ? M+1-J : lastJ;
	      int jmax = p->orientation ? M+1-lastJ : J;
	      int j = jmin;
	      for(; j < jmax; j++){
		double SRange = (X[j+1] - X[j]) * Xscale;
		double FragCovj = FragCov[j];
		if(PairMergeChimNorm > 0.0 && ChimNorm[j] > 0.0 && ChimNorm[j+1] > 0.0)
		  FragCovj = min(FragCovj, min(ChimNorm[j]*ChimQuality[j],ChimNorm[j+1]*ChimQuality[j+1]) * PairMergeChimNorm * 0.01);

		if(Outlier < PairMergeOutlierMaxSize &&
		   (OutlierFrac[j] > PairMergeOutlierFrac || 
		    (FragCovj < PairMergeOutlierMaxCov && 
		     LRange > PairMergeOutlierMinInterval &&
		     ((FragChiSq[j] < PairMergeOutlierChiSq &&
		       FragSd[j] > ExpSd[j] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
		       FragSd[j] > ExpSd[j] * PairMergeOutlierSdRatio) ||		     
		      FragCovj < PairMergeOutlierMinCov +
		      ((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != abs(J-lastJ))?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale)))){
		  if(EVERB){
		    printf("\t j=%d(%d..%d):J=%d..%d:LRange= %0.4f(Ratio= %0.3f), SRange= %0.4f, OutlierFrac= %0.3f, FragChiSq= %0.4e, Outlier= %0.4f, FragCov= %0.2f(N1=%0.2f,%0.2f), FragSd= %0.1f, ExpSd= %0.1f : skipping outlier\n",
			   j,jmin,jmax,lastJ,J,LRange,LRangeRatio,SRange, OutlierFrac[j], FragChiSq[j], Outlier, FragCov[j], ChimNorm[j]*ChimQuality[j]*0.01,ChimNorm[j+1]*ChimQuality[j+1]*0.01,FragSd[j],ExpSd[j]);
		    printf("\t    MaxOutlierFrac= %0.3f, MaxCov= %0.2f, ChiSq= %0.4e, MinSf= %0.4f, MinSr= %0.4f, SdRatio= %0.4f, MinInterval= %0.4f\n",
			   PairMergeOutlierFrac,PairMergeOutlierMaxCov,PairMergeOutlierChiSq,PairMergeOutlierMinSf,PairMergeOutlierMinSr,PairMergeOutlierSdRatio,PairMergeOutlierMinInterval);
		    fflush(stdout);
		  }

		  break;
		}
	      }
	      if(j < jmax)
		continue;// outlier will be ignored
	    }
	  }

	  double YLeft = (Y[lastI+1]-Y[lastI]) * Yscale;
	  double YRight = (Y[I] - Y[I-1]) * Yscale;
	  double XLeft = (p->orientation ? X[M+1-lastJ] - X[M+1-(lastJ+1)] : X[lastJ+1] - X[lastJ]) * Xscale;
	  double XRight = (p->orientation ? X[M+1-(J-1)] - X[M+1-J] : X[J] - X[J-1]) * Xscale;

	  if(I - lastI + J - lastJ - 2 >= PairMergeMaxOutlierBalanced){
	    if(I > lastI + 2){
	      double Middle = (Y[I-1] - Y[lastI+1]) * Yscale;
	      double MinOutlier = Middle + max(0.0, max(YLeft-XLeft,YRight-XRight));
	      Outlier = max(Outlier, MinOutlier);
	    }
	    if(J > lastJ + 2){
	      double Middle = (p->orientation ? X[M+1-(lastJ+1)] - X[M+1-(J-1)] : X[J-1] - X[lastJ+1]) * Xscale;
	      double MinOutlier = Middle + max(0.0, max(XLeft-YLeft,XRight-YRight));
	      Outlier = max(Outlier, MinOutlier);
	    }
	  }
	  if(MaxInvPalindrome && Outlier >= PairMergeMaxOutlier /* WAS40 MaxInvPalindromeSL */){
	    Ileft = min(Ileft, lastI);
	    Iright = max(Iright, I);
	  }
	  
	  if(EVERB/* HERE HERE >=2 */){
	    printf("    Ymap->id=%lld,Xmap->id=%lld,or=%d,M=%d,Xscale=%0.6f,Yscale=%0.6f:I=%d..%d(%0.3f..%0.3f),J=%d..%d(%0.3f..%0.3f),delta=%0.3f:overlap= %0.3f -> %0.3f, maxOutlier= %0.3f -> %0.3f,Ileft=%d,Iright=%d\n",
		   Ymap->id,Xmap->id,p->orientation,M,Xscale,Yscale,lastI,I,refstartpos,refendpos,lastJ,J,qrystartpos,qryendpos, fabs(refendpos - refstartpos - qryendpos + qrystartpos),
		   origoverlap, overlap, maxOutlier, Outlier,Ileft,Iright);
	    fflush(stdout);
	  }

	  if(Outlier > maxOutlier){
	    maxOutlier = Outlier;
	    IleftOutlier = lastI;
	    IrightOutlier = I;
	    JleftOutlier = lastJ;
	    JrightOutlier = J;
	  }
	  maxOutlierM = max(maxOutlierM,OutlierM);
	}
      }

      /* end truncation flags */
      size_t YL = Ymap->Mask[0] ? Ymap->Mask[0][1] : 0;
      size_t YR = Ymap->Mask[0] ? Ymap->Mask[0][N+1] : 0;
      size_t XL = Xmap->Mask[0] ? Xmap->Mask[0][p->orientation ? M+1 : 1] : 0;
      size_t XR = Xmap->Mask[0] ? Xmap->Mask[0][p->orientation ? 1 : M+1] : 0;

      if(VERB/* HERE >=2 */){
	if(PairMergeRepeat)
	  printf("%d:",pairmergeIter);
	printf("Alignment %llu/%llu:Ymap=%p:mapid=%d,id=%lld(paired=%d),Xmap=%p:mapid=%d,id=%lld(paired=%d):score=%0.6f,logPV=%0.2f,or=%d,sc=%d,pairs=%d,Y[LI=%d]=%0.3f,Y[RI=%d]=%0.3f,Y[N+1=%d]=%0.3f,YL=0x%lx,YR=0x%lx\n\t X[LJ=%d]=%0.3f,X[RJ=%d]=%0.3f,X[M+1=%d]=%0.3f,XL=0x%lx,XR=0x%lx,PairMerge=%0.1f,%0.2f,%0.1f,%0.1f,%d,overlap=%0.3f,Ends=%0.3e,%0.3e(%d,%d), maxOutlier=%0.3e(I=%d..%d),%d,Ileft=%d,Iright=%d\n",
	       (unsigned long long)i,(unsigned long long)numaligns, Ymap, Ymap->mapid, Ymap->id, Ymap->paired, Xmap, Xmap->mapid, Xmap->id,Xmap->paired,p->score,p->logPV,p->orientation, p->scaleID, p->numpairs,
	       LI,Y[LI],RI,Y[RI],  N+1, Y[N+1], YL, YR, LJ, p->orientation ? (X[M+1] - X[M+1 - LJ]) : X[LJ], RJ, p->orientation ? (X[M+1] - X[M+1 - RJ]) : X[RJ], M+1, X[M+1], XL, XR,
	       PairMerge, PairMergeMaxEnd, PairMergeMaxOutlier,PairMergeMaxEndKB, PairMergeMaxEndN, overlap, leftEnd, rightEnd, leftEndN,rightEndN, maxOutlier, IleftOutlier,IrightOutlier,maxOutlierM,Ileft, Iright);
	if(VERB>=3)
	  printf("\t Y[N+1]=%0.3f, X[M+1]=%0.3f\n",Y[N+1], X[Xmap->numsite[0]+1]);
	fflush(stdout);
      }

      if(Xmap->paired || Ymap->paired){
	if(DEBUG) assert(p->mapid1 != p->mapid2);// self alignments for -MaxPalindrome or -MaxInvPalindrome are checked first

	continue;/* cannot do anything if one of the maps has already been modified or removed */
      }

      double MinOverlap = PairMerge;
      double OverlappedMargin = PairMergeOverlappedMargin;
      int LeftSmallerY = (Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]) + OverlappedMargin ) ? 1 : 0;
      int LeftSmallerX = (Y[LI] > (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]) - OverlappedMargin ) ? 1 : 0;
      int RightSmallerY = (Y[N+1] - Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]) + OverlappedMargin) ? 1 : 0;
      int RightSmallerX = (Y[N+1] - Y[RI] > (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]) - OverlappedMargin) ? 1 : 0;
      int Overlapped = ((LeftSmallerY && RightSmallerY) || (LeftSmallerX && RightSmallerX)) ? 1 : 0;

      int merge = 0;
      if(!Cfirst && Ymap->mapid == Xmap->mapid){/* handle -MaxPalindrome or -MaxInvPalindrome */
	if(DEBUG) assert(Xmap == Ymap);
	if(DEBUG) assert(MaxPalindrome > 0.0 || MaxInvPalindrome > 0.0);
	if(DEBUG) assert(pairmergeIter == 1 || SplitSegDup);
	if(DEBUG && !(PairMergeHmap <= 0)){
	  printf("pairmergeIter= %d, PairMergeHmap= %d\n", pairmergeIter, PairMergeHmap);
	  fflush(stdout);
	  assert(PairMergeHmap <= 0);
	}

	if(MaxPalindrome > 0.0){
	  merge = (overlap >= MaxPalindrome && K >= MaxPalindromeN && (leftEndN < PairMergeMaxEndN || leftEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB)) &&
		   (rightEndN < PairMergeMaxEndN || rightEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB))) ? 1 : 0;/* Palindromic overlap */

	  if(PairMergeMaxOutlier >= 0.0 && (maxOutlier > PairMergeMaxOutlier || maxOutlierM > PairMergeMaxOutlierM))
	    merge = 0;/* cannot be matching Palindrome because internal outlier is too large */

	  if(merge){/* break the contig at midpoint of Palindromic (overlap) region + MaxPalindromeSL kb */
	    Ymap->paired = 1+nummaps;

	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;
	    pmap->origid = Ymap->origid;
	    pmap->id = Ymap->id;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Ymap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;

	    double PalMidpoint = (Y[LI] + Y[RI]) * 0.5;
	    if(PalMidpoint < Y[N+1] * 0.5){/* trim off left end at PalMidpoint - MaxPalindromeSL kb */
	      double BreakPoint = PalMidpoint - MaxPalindromeSL;// HERE HERE HERE : should extend left from first label left of BreakPoint

	      pmap->trim(BreakPoint, Y[N+1]);
	      if(!pmap->Mask[0]){
		pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
	      }
	      pmap->Mask[0][1] |= MaxPalindromeMask;// WAS303 END_NOEXT;
	      if(VERB/* HERE >=2 */){
		printf("    Trimming off %0.3f kb from left end of Y due to Palindrome from %0.3f to %0.3f and left end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
		       BreakPoint, Y[LI],Y[RI], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid);
		fflush(stdout);
	      }
	    } else {/* trim off right end at PalMidpoint + MaxPalindromeSL kb */
	      double BreakPoint = PalMidpoint + MaxPalindromeSL;// HERE HERE HERE : should extend right from first label right of BreakPoint
	      
	      pmap->trim(0.0,BreakPoint);
	      if(!pmap->Mask[0]){
		pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
	      }
	      pmap->Mask[0][pmap->numsite[0]+1] |= MaxPalindromeMask;// WAS303 END_NOEXT;

	      if(VERB/* HERE >=2 */){
		printf("    Trimming off right end of Y from %0.3f kb to %0.3f kb due to Palindrome from %0.3f to %0.3f and right end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
		       BreakPoint, Y[N+1], Y[LI],Y[RI], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid);
		fflush(stdout);
	      }
	    }
	    continue;
	  }
	}

	if(MaxInvPalindrome > 0.0){

	  /* first check for case with no large endoutlier and large internal outlier(s) */
	  /* range of large internal outliers (over PairMergeMaxOutlier) is Y[Ileft .. Iright] */
	  double Yendoverlap = (Y[RI]-Y[LI]) - /* NEW94 */max(0.0,Y[Iright] - Y[Ileft]);

	  merge = (Yendoverlap /* WAS40 overlap */ >= MaxInvPalindrome * 2.0 && (RI-LI - /* NEW94 */max(0, Iright - Ileft))/* WAS40 K */ >= MaxInvPalindromeN * 2
		   && (leftEndN < PairMergeMaxEndN || leftEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB))
		   && (rightEndN < PairMergeMaxEndN || rightEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB))) ? 1 : 0;

	  if(merge && Y[Iright] - Y[Ileft] /* WAS40 maxOutlier */ >= MaxInvPalindromeSL){
	    if(DEBUG) assert(0 < Ileft && Ileft < Iright && Iright <= N);

	    Ymap->paired = 1+nummaps;

	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;
	    pmap->id = Ymap->id;
	    pmap->origid = Ymap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Ymap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;

	    if(Y[Ileft] < Y[N+1]-Y[Iright]){/* trim off left end up to Y[Ileft] */
	      pmap->trim(Y[Ileft] - MININTERVAL, Y[N+1]);
	      if(!pmap->Mask[0]){
		pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
	      }
	      pmap->Mask[0][1] |= MaxInvPalindromeMask;// WAS303 END_NOEXT;

	      if(VERB/* HERE >=2 */){
		printf("    Trimming off %0.3f kb from left end of Y due to InvPalindrome outlier from %0.3f to %0.3f and left end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
		       Y[Ileft]-0.02, Y[Ileft], Y[Iright], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid);
		fflush(stdout);
	      }
	    } else {/* trim off right end beyond Y[Iright] */
	      pmap->trim(0.0,Y[Iright] + MININTERVAL);
	      if(!pmap->Mask[0]){
		pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
	      }
	      pmap->Mask[0][pmap->numsite[0]+1] |= MaxInvPalindromeMask;// WAS303 END_NOEXT;

	      if(VERB/* HERE >=2 */){
		printf("    Trimming off right end of Y from %0.3f kb to %0.3f kb due to InvPalindrome outlier from %0.3f to %0.3f and right end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
		       Y[Iright] + 0.02, Y[N+1], Y[Ileft],Y[Iright], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid);
		fflush(stdout);
	      }
	    }
	    continue;
	  }
	  
	  /* next check for case with up to one large endoutlier (if Het Inversion was too large to be a single outlier) */
	  merge = (Yendoverlap /* WAS40 overlap */ >= MaxInvPalindrome && (RI-LI - /*NEW94*/max(0, Iright - Ileft)) /* WAS40 K */ >= MaxInvPalindromeN 
		   && (leftEndN < PairMergeMaxEndN || leftEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB)
		       || rightEndN < PairMergeMaxEndN || rightEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB))) ? 1 : 0;

	  if(PairMergeMaxOutlier >= 0.0 && (maxOutlier > PairMergeMaxOutlier || maxOutlierM > PairMergeMaxOutlierM))
	    merge = 0;/* cannot be matching Inverted Palindrome because internal outlier is too large */

	  if(merge){// cut at either Y[LI] or Y[N+1-RJ /* WAS40 LJ */ ], whichever results in the largest piece
	    Ymap->paired = 1+nummaps;

	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;
	    pmap->id = Ymap->id;
	    pmap->origid = Ymap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Ymap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;

	    if(max(Y[LI], Y[N+1]-Y[LI]) > max(Y[N+1-RJ], Y[N+1] - Y[N+1-RJ])){ /* cut at Y[LI] */
	      if(Y[LI] > Y[N+1] - Y[LI]){/* keep Y[0] .. Y[LI] + 0.02 */
		pmap->trim(0.0,Y[LI] + MININTERVAL);
		if(!pmap->Mask[0]){
		  pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		  memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
		}
		pmap->Mask[0][pmap->numsite[0]+1] |= MaxInvPalindromeMask;// WAS303 END_NOEXT;

		if(VERB/* HERE >=2 */){
		  printf("    Trimming off right end of Y from %0.3f kb to %0.3f kb due to InvPalindrome from %0.3f to %0.3f and right end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
			 Y[LI] + 0.02, Y[N+1], min(Y[LI],Y[N+1-RJ]),max(Y[LI],Y[N+1-RJ]), pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid);
		  fflush(stdout);
		}
	      } else {/* keep Y[LI]-0.02 .. Y[N+1] */
		pmap->trim(Y[LI] - MININTERVAL, Y[N+1]);
		if(!pmap->Mask[0]){
		  pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		  memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
		}
		pmap->Mask[0][1] |= MaxInvPalindromeMask;// WAS303 END_NOEXT;

		if(VERB/* HERE >=2 */){
		  printf("    Trimming off %0.3f kb from left end of Y due to InvPalindrome from %0.3f to %0.3f and left end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
			 Y[LI]-0.02, min(Y[LI],Y[N+1-RJ]), max(Y[LI],Y[N+1-RJ]), pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid);
		  fflush(stdout);
		}
	      }
	    } else {/* cut at Y[N+1-RJ] */
	      if(Y[N+1-RJ] > Y[N+1] - Y[N+1-RJ]){/* keep Y[0] .. Y[N+1-RJ] + 0.02 */
		pmap->trim(0.0,Y[N+1-RJ] + MININTERVAL);
		if(!pmap->Mask[0]){
		  pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		  memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
		}
		pmap->Mask[0][pmap->numsite[0]+1] |= MaxInvPalindromeMask; // WAS303 END_NOEXT;

		if(VERB/* HERE >=2 */){
		  printf("    Trimming off right end of Y from %0.3f kb to %0.3f kb due to Palindrome from %0.3f to %0.3f and right end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
			 Y[N+1-RJ] + 0.02, Y[N+1], min(Y[LI],Y[N+1-RJ]),max(Y[LI],Y[N+1-RJ]), pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid);
		  fflush(stdout);
		}
	      } else {/* keep Y[N+1-RJ] - 0.02 .. Y[N+1] */
		pmap->trim(Y[N+1-RJ] - MININTERVAL, Y[N+1]);
		if(!pmap->Mask[0]){
		  pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		  memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
		}
		pmap->Mask[0][1] |= MaxInvPalindromeMask; // WAS303 END_NOEXT;

		if(VERB/* HERE >=2 */){
		  printf("    Trimming off %0.3f kb from left end of Y due to InvPalindrome from %0.3f to %0.3f and left end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
			 Y[N+1-RJ]-0.02, min(Y[LI],Y[N+1-RJ]), max(Y[LI],Y[N+1-RJ]), pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid);
		  fflush(stdout);
		}
	      }
	    }
	    continue;
	  }
	}

	continue;/* Self alignments are only used to handle -MaxPalindronme and -MaxInvPalindrome by breaking the contig (if needed) */
      }

      if(PairMergeOverlapped > 0.0 && Overlapped)
	MinOverlap = PairMergeOverlapped;

      merge = (overlap >= MinOverlap && (leftEndN < PairMergeMaxEndN || leftEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB)) &&
	       (rightEndN < PairMergeMaxEndN || rightEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB))) ? 1 : 0;/* pairmerge candidate found ? */

      if(PairMergeMaxOutlier >= 0.0 && (maxOutlier > PairMergeMaxOutlier || maxOutlierM > PairMergeMaxOutlierM))
	merge = 0;/* cannot merge because internal outlier is too large */

      if(EVERB){
	printf("\t merge=%d: overlap=%0.3f/%0.3f, leftEnd= %0.3f,%d, rightEnd= %0.3f,%d (%0.3f,%d), maxOut= %0.3f,%d, Overlapped=%d(LeftY=%0.3f,%d,LeftX=%0.3f,%d,RightY=%0.3f,%d,RightX=%0.3f,%d,Margin=%0.3f)\n",
	       merge,overlap,MinOverlap,leftEnd,leftEndN,rightEnd,rightEndN, max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB), PairMergeMaxEndN, maxOutlier,maxOutlierM,Overlapped,
	       Y[LI], LeftSmallerY, p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ], LeftSmallerX, Y[N+1]-Y[RI], RightSmallerY, p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ], RightSmallerX, OverlappedMargin);
	fflush(stdout);
      }

      if(merge && TrimmedNoMerge && (Xmap->Mask[0] || Ymap->Mask[0])){/* check if -TrimmedNoMerge blocks merging */

	/* Left end of Y is marked and overlaps ("left" end of) X */ 
	int Yleft = (Ymap->Mask[0] && (Ymap->Mask[0][1] & END_NOEXT) && Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ])) ? Ymap->Mask[0][1] : 0;

	/* "Left" end of X is marked and overlaps (left end of) Y */
	int Xleft = (Xmap->Mask[0] && (Xmap->Mask[0][p->orientation ? M+1 : 1] & END_NOEXT) && Y[LI] > (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ])) ? Xmap->Mask[0][p->orientation ? M+1:1] : 0;

	/* Right end of Y is marked and overlaps ("right" end of) X */
	int Yright = (Ymap->Mask[0] && (Ymap->Mask[0][N+1] & END_NOEXT) && Y[N+1]-Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ])) ? Ymap->Mask[0][N+1] : 0;
	
	/* "Right" end of X is marked and overlaps (right end of) Y */
	int Xright = (Xmap->Mask[0] && (Xmap->Mask[0][p->orientation ? 1 : M+1] & END_NOEXT) && Y[N+1]-Y[RI] > (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ])) ? Xmap->Mask[0][p->orientation ? 1 : M+1] : 0;

	int MarkedEnds = (Yleft ? 1 : 0) + (Xleft ? 1 : 0) + (Yright ? 1 : 0) + (Xright ? 1 : 0);
	int MarkedEnds2 = ((Yleft & END_CHR) ? 1 : 0) + ((Xleft & END_CHR) ? 1 : 0) + ((Yright & END_CHR) ? 1 : 0) + ((Xright & END_CHR) ? 1 : 0);/* Number of chr ends (or fragile site ends) */
	size_t MarkedMask = (Yleft ? Yleft : 0xffffffff) & (Xleft ? Xleft : 0xffffffff) & (Yright ? Yright : 0xffffffff) & (Yleft ? Yleft : 0xffffffff) & TrimmedNoMergeBoth;
	if(DEBUG) assert(MarkedEnds <= 2);/* only 2 out of 4 ends can be overlaped by the other map */

	int LeftSmallerY = (Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]) + TrimmedNoMergeMaxEnd) ? 1 : 0;
	int LeftSmallerX = (Y[LI] > (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]) - TrimmedNoMergeMaxEnd) ? 1 : 0;
	int RightSmallerY = (Y[N+1] - Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]) + TrimmedNoMergeMaxEnd) ? 1 : 0;
	int RightSmallerX = (Y[N+1] - Y[RI] > (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]) - TrimmedNoMergeMaxEnd) ? 1 : 0;
	bool FullyOverlapped = ((LeftSmallerY && RightSmallerY) || (LeftSmallerX && RightSmallerX));

	if(DEBUG && TrimmedNoMergeBoth && PairMergeHmap > 0){
	  printf("-TrimmedNoMerge non-zero 3rd arg 0x%lx not compatible with -pairmergeHmap %d ... \n",TrimmedNoMergeBoth,PairMergeHmap);
	  fflush(stdout);exit(1);
	}

	if(TrimmedNoMergeStrict && !FullyOverlapped && ((Ymap->Mask[0] && ((Ymap->Mask[0][N+1] & END_NOEXT) || (Ymap->Mask[0][1] & END_NOEXT))) ||
							(Xmap->Mask[0] && ((Xmap->Mask[0][M+1] & END_NOEXT) || (Xmap->Mask[0][1] & END_NOEXT))))){
	  if(VERB/* HERE >=2 */){
	    printf("   Skipped merging Ymap->id=%lld,Xmap->id=%lld due to TrimmedNoMergeStrict=%d,overlap=%0.3f,%d, Xleft=%x,Xright=%x,Yleft=%x,Yright=%d,SmallerY:L=%d,R=%d;SmallerX:L=%d,R=%d;Both=%lx,Ends=%d%lx,MaskY=%lx,%lx,MaskX=%lx,%lx\n",
		   Ymap->id,Xmap->id,TrimmedNoMergeStrict,overlap,p->numpairs,Xleft,Xright,Yleft,Yright,LeftSmallerY,RightSmallerY,LeftSmallerX,RightSmallerX,TrimmedNoMergeBoth,MarkedEnds,MarkedMask,
		   Ymap->Mask[0] ? Ymap->Mask[0][1] : 0, Ymap->Mask[0] ? Ymap->Mask[0][N+1] : 0, 
		   Xmap->Mask[0] ? Xmap->Mask[0][p->orientation ? M+1 : 1] : 0, Xmap->Mask[0] ? Xmap->Mask[0][p->orientation ? 1 : M+1] : 0);
	      fflush(stdout);
	    }
	    merge = 0;
	    if(PairMergeHmap > 0 && SplitSegDup > 0.0)
	      continue;// NEW3 : avoid merging end overlaps with marked ends with PairMergeHmap to avoid chimeric contigs
	}

	if(MarkedEnds && (overlap < TrimmedNoMerge || MarkedEnds2) && (!TrimmedNoMergeBoth || (MarkedEnds >= 2 && MarkedMask))){
	  if(!FullyOverlapped){/* smaller map is NOT fully overlapped by a larger map (even allowing for TrimmedNoMergeMaxEnd wiggle room) */
	    if(VERB/* HERE >=2 */){
	      printf("   Skipped merging Ymap->id=%lld,Xmap->id=%lld due to TrimmedNoMerge= %0.3f,overlap=%0.3f,%d, Xleft=%x,Xright=%x,Yleft=%x,Yright=%d,SmallerY:L=%d,R=%d;SmallerX:L=%d,R=%d;Both=%lx,Ends=%d,Mask=%lx\n",
		     Ymap->id,Xmap->id,TrimmedNoMerge,overlap,p->numpairs,Xleft,Xright,Yleft,Yright,LeftSmallerY,RightSmallerY,LeftSmallerX,RightSmallerX,TrimmedNoMergeBoth,MarkedEnds,MarkedMask);
	      fflush(stdout);
	    }
	    merge = 0;
	    if(PairMergeHmap > 0 && SplitSegDup > 0.0)
	      continue;// NEW3 : avoid merging end overlaps with marked ends with PairMergeHmap to avoid chimeric contigs
	  }
	}
	
	if(TrimmedNoMergeSegDupMask && !FullyOverlapped){/* Disallow merge if the overlap region lables of either map has Mask bits that intersect with TrimmedNoMergeSegDupMask,
							    unless the overlap extends beyond Masked region by at least TrimmedNoMergeSegDupFlank labels in the non-extending region */
	  if(DEBUG && TrimmedNoMergeSegDupMask && PairMergeHmap > 0.0){
	    printf("-TrimmedNoMerge non-zero 4th arg 0x%lx not compatible with -pairmergeHmap %d ... \n",TrimmedNoMergeSegDupMask,PairMergeHmap);
	    fflush(stdout);exit(1);
	  }

	  int mLI1 = N+1, mRI1 = -1, mLI2 = N+1, mRI2 = -1, mLJ1 = M+1, mRJ1 = -1, mLJ2 = M+1, mRJ2 = -1;/* range of labels that are masked on Y or X : leftmost and rightmost contigous ranges */
	  if(Ymap->Mask[0]){
	    for(int I = LI; I <= RI; I++){
	      if(Ymap->Mask[0][I] & TrimmedNoMergeSegDupMask){
		mLI1 = min(I,mLI1);
		mRI2 = max(I,mRI2);
	      }
	    }
	    for(int I = mLI1 + 1; I <= mRI2; I++){
	      if(!(Ymap->Mask[0][I] & TrimmedNoMergeSegDupMask)){// break if there are TrimmedNoMergeSegDupFlank unmarked labels i = I .. I + TrimmedNoMergeSegDupFlank -1 
		int cnt = 1;
		for(int i = I+1; i <= mRI2 && i < I + TrimmedNoMergeSegDupFlank; i++){
		  if((Ymap->Mask[0][i] & TrimmedNoMergeSegDupMask))
		    break;
		  cnt++;
		}
		if(cnt >= TrimmedNoMergeSegDupFlank)
		  break;
	      }
	      mRI1 = max(I,mRI1);
	    }
	    for(int I = mRI2 - 1; I >= mLI1; I--){
	      if(!(Ymap->Mask[0][I] & TrimmedNoMergeSegDupMask)){// break if there are TrimmedNoMergeSegDupFlank unmarked labels i = I .. I - TrimmedNoMergeSegDupFlank + 1
		int cnt = 1;
		for(int i = I-1; i >= mLI1 && i > I - TrimmedNoMergeSegDupFlank; i--){
		  if((Ymap->Mask[0][i] & TrimmedNoMergeSegDupMask))
		    break;
		  cnt++;
		}
		if(cnt >= TrimmedNoMergeSegDupFlank)
		  break;
	      }
	      mLI2 = min(I, mLI2);
	    }
	  }

	  if(Xmap->Mask[0]){
	    if(p->orientation){
	      for(int J = LJ; J <= RJ; J++){
		if(Xmap->Mask[0][M+1-J] & TrimmedNoMergeSegDupMask){
		  mLJ1 = min(J,mLJ1);
		  mRJ2 = max(J,mRJ2);
		}
	      }
	      for(int J = mLJ1 + 1; J <= mRJ2; J++){
		if(!(Xmap->Mask[0][M+1-J] & TrimmedNoMergeSegDupMask)){// break if there are TrimmedNoMergeSegDupFlank unmarked labels j = J .. J + TrimmedNoMergeSegDupFlank -1 
		  int cnt = 1;
		  for(int j = J+1; j <= mRJ2 && j < J + TrimmedNoMergeSegDupFlank; j++){
		    if((Xmap->Mask[0][M+1-j] & TrimmedNoMergeSegDupMask))
		      break;
		    cnt++;
		  }
		  if(cnt >= TrimmedNoMergeSegDupFlank)
		    break;
		}
		mRJ1 = max(J,mRJ1);
	      }
	      for(int J = mRJ2 - 1; J >= mLJ1; J--){
		if(!(Xmap->Mask[0][M+1-J] & TrimmedNoMergeSegDupMask)){// break if there are TrimmedNoMergeSegDupFlank unmarked labels j = J .. J - TrimmedNoMergeSegDupFlank + 1
		  int cnt = 1;
		  for(int j = J-1; j >= mLJ1 && j > J - TrimmedNoMergeSegDupFlank; j--){
		    if((Xmap->Mask[0][M+1-j] & TrimmedNoMergeSegDupMask))
		      break;
		    cnt++;
		  }
		  if(cnt >= TrimmedNoMergeSegDupFlank)
		    break;
		}
		mLJ2 = min(J, mLJ2);
	      }
	    } else {
	      for(int J = LJ; J <= RJ; J++){
		if(Xmap->Mask[0][J] & TrimmedNoMergeSegDupMask){
		  mLJ1 = min(J,mLJ1);
		  mRJ2 = max(J,mRJ2);
		}
	      }
	      for(int J = mLJ1 + 1; J <= mRJ2; J++){
		if(!(Xmap->Mask[0][J] & TrimmedNoMergeSegDupMask)){// break if there are TrimmedNoMergeSegDupFlank unmarked labels j = J .. J + TrimmedNoMergeSegDupFlank -1 
		  int cnt = 1;
		  for(int j = J+1; j <= mRJ2 && j < J + TrimmedNoMergeSegDupFlank; j++){
		    if((Xmap->Mask[0][j] & TrimmedNoMergeSegDupMask))
		      break;
		    cnt++;
		  }
		  if(cnt >= TrimmedNoMergeSegDupFlank)
		    break;
		}
		mRJ1 = max(J,mRJ1);
	      }
	      for(int J = mRJ2 - 1; J >= mLJ1; J--){
		if(!(Xmap->Mask[0][J] & TrimmedNoMergeSegDupMask)){// break if there are TrimmedNoMergeSegDupFlank unmarked labels j = J .. J - TrimmedNoMergeSegDupFlank + 1
		  int cnt = 1;
		  for(int j = J-1; j >= mLJ1 && j > J - TrimmedNoMergeSegDupFlank; j--){
		    if((Xmap->Mask[0][j] & TrimmedNoMergeSegDupMask))
		      break;
		    cnt++;
		  }
		  if(cnt >= TrimmedNoMergeSegDupFlank)
		    break;
		}
		mLJ2 = min(J, mLJ2);
	      }
	    }
	  }

	  if(VERB>=2 && i==18){
	    printf("RightSmallerY=%d,LeftSmallerY=%d,RightSmallerX=%d,LeftSmallerX=%d,FullyOverlapped=%d:mLI1=%d,mRI1=%d,mLI2=%d,mRI2=%d,mLJ1=%d,mRJ1=%d,mLJ2=%d,mRJ2=%d\n",
		   RightSmallerY,LeftSmallerY,RightSmallerX,LeftSmallerX,FullyOverlapped,mLI1,mRI1,mLI2,mRI2,mLJ1,mRJ1,mLJ2,mRJ2);
	    fflush(stdout);
	  }

	  if(RightSmallerY && LeftSmallerX){
	    if((mLI2 <= mRI2 && !(LI <= mLI2 - TrimmedNoMergeSegDupFlank)) ||
	       (mLJ1 <= mRJ1 && !(RJ >= mRJ1 + TrimmedNoMergeSegDupFlank))){
	      if(VERB/* HERE >=2 */){
		printf("   Skipped merging Ymap->id=%lld,Xmap->id=%lld due to TrimmedNoMerge SegDup= 0x%lx, Flank= %d : overlap=%0.3f,%d, Xleft=%x,Xright=%x,Yleft=%x,Yright=%d,SmallerY:L=%d,R=%d;SmallerX:L=%d,R=%d;mLI=%d,%d,mRI=%d,%d,mLJ=%d,%d,mRJ=%d,%d\n",
		       Ymap->id,Xmap->id,TrimmedNoMergeSegDupMask,TrimmedNoMergeSegDupFlank,overlap,p->numpairs,Xleft,Xright,Yleft,Yright,LeftSmallerY,RightSmallerY,LeftSmallerX,RightSmallerX,
		       mLI1,mLI2,mRI1,mRI2,mLJ1,mLJ2,mRJ1,mRJ2);
		fflush(stdout);
	      }
	      merge = 0;
	      continue;
	    }
	  } else if(RightSmallerX && LeftSmallerY){
	    if((mLI1 <= mRI1 && !(RI >= mRI1 + TrimmedNoMergeSegDupFlank)) ||
	       (mLJ2 <= mRJ2 && !(LJ <= mLJ2 - TrimmedNoMergeSegDupFlank))){
	      if(VERB/* HERE >=2 */){
		printf("   Skipped merging Ymap->id=%lld,Xmap->id=%lld due to TrimmedNoMerge SegDup= 0x%lx, Flank= %d : overlap=%0.3f,%d, Xleft=%x,Xright=%x,Yleft=%x,Yright=%d,SmallerY:L=%d,R=%d;SmallerX:L=%d,R=%d;mLI=%d,%d,mRI=%d,%d,mLJ=%d,%d,mRJ=%d,%d\n",
		       Ymap->id,Xmap->id,TrimmedNoMergeSegDupMask,TrimmedNoMergeSegDupFlank,overlap,p->numpairs,Xleft,Xright,Yleft,Yright,LeftSmallerY,RightSmallerY,LeftSmallerX,RightSmallerX,
		       mLI1,mLI2,mRI1,mRI2,mLJ1,mLJ2,mRJ1,mRJ2);
		fflush(stdout);
	      }
	      merge = 0;
	      continue;
	    }
	  }
	}
      } // if(merge ...)
      
      if(BreakEndOutlier > 0.0 && SplitSegDup <= 0.0 && PairMergeHmap <= 0){
        bool LeftEndOutlier = (leftEnd > BreakEndOutlierE && leftEndN > BreakEndOutlierEN);
        bool RightEndOutlier = (rightEnd > BreakEndOutlierE && rightEndN > BreakEndOutlierEN);

	if(overlap >= BreakEndOutlier && K >= BreakEndOutlierN && (LeftEndOutlier || RightEndOutlier)){/* Assume endoutliers ends are chimeric and break them off */
	  
	  if(DEBUG) assert(!Ymap->paired);

	  double SplitLen = 0.0;/* largest broken off piece of Ymap (or Xmap) so far */

	  if(LeftEndOutlier){  /* break off left end of Ymap */
	    Ymap->paired = 1+nummaps;
	    
	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;

	    if(IDmax < IDMAX)
	      pmap->id = Ymap->id + IDmax;
	    else {
	      #pragma omp critical
	      {
		pmap->id = ++IDmax;
	      }
	    }
	    if(DEBUG && pmap->id <= 0){
	      printf("pmap->id= %lld -> %lld, IDmax= %lld\n",Ymap->id, pmap->id, IDmax);
	      fflush(stdout);
	      assert(pmap->id > 0);
	    }

	    pmap->origid = Ymap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Ymap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;
	    
	    pmap->trim(0.0, Y[LI] + MININTERVAL);

	    SplitLen = pmap->site[0][pmap->numsite[0]+1];
	    
	    if(VERB/* HERE >=2 */){
	      printf("    Breaking off left end of Y up to %0.3f due to suspected Chimeric Alignment from %0.3f to %0.3f. New map: N=%d,Len=%0.3f,(mapid=%d,id=%lld),splits=%d\n",
		     Y[LI], Y[LI],Y[RI], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	      fflush(stdout);
	    }
	  }
	  
	  if(RightEndOutlier){/* break off right end of Ymap */
	    if(Ymap->paired)
	      splitcnt++;

	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;

	    if(IDmax < IDMAX){
	      pmap->id = Ymap->id + IDmax;
	      if(Ymap->paired)
		pmap->id += IDmax;
	    } else {
	      #pragma omp critical
	      {
		pmap->id = ++IDmax;
	      }
	    }
	    if(DEBUG && pmap->id <= 0){
	      printf("pmap->id= %lld -> %lld, IDmax= %lld\n",Ymap->id, pmap->id, IDmax);
	      fflush(stdout);
	      assert(pmap->id > 0);
	    }

	    pmap->origid = Ymap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Ymap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;

	    if(Y[N+1] - Y[RI] > SplitLen)
	      Ymap->paired = nummaps;

	    pmap->trim(Y[RI] - MININTERVAL, Y[N+1]);

	    SplitLen = max(SplitLen, pmap->site[0][pmap->numsite[0]+1]);

	    if(VERB/* HERE >=2 */){
	      printf("    Breaking off right end of Y from %0.3f kb to %0.3f due to suspected Chimeric alignment from %0.3f to %0.3f. New map: N=%d,Len=%0.3f(mapid=%d,id=%lld),splits=%d\n",
		     Y[RI], Y[N+1], Y[LI],Y[RI],pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	      fflush(stdout);
	    }
	  }

	  if(Ymap->paired){/* Save remainder of Ymap */
	    splitcnt++;

	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;
	    pmap->id = Ymap->id;
	    pmap->origid = Ymap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Ymap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;
	    
	    double leftend = LeftEndOutlier ? Y[LI] - MININTERVAL : 0.0;
	    double rightend = RightEndOutlier ? Y[RI] + MININTERVAL : Y[N+1];
	    if(rightend - leftend > SplitLen)
	      Ymap->paired = nummaps;

	    pmap->trim(leftend, rightend);

	    if(VERB/* HERE >=2 */){
	      printf("    Trimming endoutliers of Y and keeping only region from %0.3f to %0.3f. New map: N=%d,Len=%0.3f(mapid=%d,id=%lld),splits=%d\n",
		     leftend,rightend, pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	      fflush(stdout);
	    }
	  }

	  if(DEBUG) assert(!Xmap->paired);

	  SplitLen = 0.0;

	  if(LeftEndOutlier){ /* left end of Xmap */
	    Xmap->paired = 1+nummaps;
	    
	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Xmap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;

	    if(IDmax < IDMAX)
	      pmap->id = Xmap->id + IDmax;
	    else {
	      #pragma omp critical
	      {
		pmap->id = ++IDmax;
	      }
	    }
	    if(DEBUG && pmap->id <= 0){
	      printf("pmap->id= %lld -> %lld, IDmax= %lld\n",Xmap->id, pmap->id, IDmax);
	      fflush(stdout);
	      assert(pmap->id > 0);
	    }

	    pmap->origid = Xmap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Xmap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;
	    
	    double leftend = p->orientation ? max(0.0, X[M+1-LJ] - MININTERVAL) : 0.0;
	    double rightend = p->orientation ? X[M+1] : min(X[M+1], X[LJ] + MININTERVAL);

	    pmap->trim(leftend, rightend);

	    SplitLen = pmap->site[0][pmap->numsite[0]+1];

	    if(VERB/* HERE >=2 */){
	      printf("    Breaking off %s end of X from %0.3f to %0.3f due to suspected Chimeric alignment from %0.3f to %0.3f. New map: N=%d,Len=%0.3f(mapid=%d,id=%lld),splits=%d\n",
		     p->orientation ? "right" : "left", leftend, rightend, p->orientation ? X[M+1-RJ] : X[LJ], p->orientation ? X[M+1-LJ] : X[RJ], 
		     pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	      fflush(stdout);
	    }
	  }

	  if(RightEndOutlier){ /* right end of Xmap */
	    if(Xmap->paired)
	      splitcnt++;
	    
	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Xmap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;

	    if(IDmax < IDMAX){
	      pmap->id = Xmap->id + IDmax;
	      if(Xmap->paired)
		pmap->id += IDmax;
	    } else {
	      #pragma omp critical
	      {
		pmap->id = ++IDmax;
	      }
	    }
	    if(DEBUG && pmap->id <= 0){
	      printf("pmap->id= %lld -> %lld, IDmax= %lld\n",Xmap->id, pmap->id, IDmax);
	      fflush(stdout);
	      assert(pmap->id > 0);
	    }

	    pmap->origid = Xmap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Xmap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;

	    if((p->orientation ? X[M+1-RJ] : X[M+1]-X[RJ]) > SplitLen)
	      Xmap->paired = nummaps;

	    double leftend = p->orientation ? 0.0 : max(0.0,X[RJ] - MININTERVAL);
	    double rightend = p->orientation ? min(X[M+1], X[M+1-RJ] + MININTERVAL) : X[M+1];

	    pmap->trim(leftend, rightend);

	    SplitLen = max(SplitLen, pmap->site[0][pmap->numsite[0]+1]);

	    if(VERB/* HERE >=2 */){
	      printf("    Breaking off %s end of X from %0.3f kb to %0.3f due to suspected Chimeric alignment from %0.3f to %0.3f. New map: N=%d,Len=%0.3f(mapid=%d,id=%lld),splits=%d\n",
		     p->orientation ? "left" : "right", leftend, rightend, p->orientation ? X[M+1-RJ] : X[LJ], p->orientation ? X[M+1-LJ] : X[RJ],
		     pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	      fflush(stdout);
	    }
	  }

          if(Xmap->paired){/* Save remainder of Xmap */
	    splitcnt++;

	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Xmap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;
	    pmap->id = Xmap->id;
	    pmap->origid = Xmap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Ymap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;
	    
	    double leftend = p->orientation ? (RightEndOutlier ? X[M+1-RJ] - MININTERVAL : 0.0) : (LeftEndOutlier ? X[LJ] - MININTERVAL : 0.0);
	    double rightend = p->orientation ? (LeftEndOutlier ? X[M+1-LJ] + MININTERVAL : X[M+1]) : (RightEndOutlier ? X[RJ] + MININTERVAL : X[M+1]);
	    if(rightend - leftend > SplitLen)
	      Xmap->paired = nummaps;

	    pmap->trim(leftend, rightend);

	    if(VERB/* HERE >=2 */){
	      printf("    Trimming endoutliers of X and keeping only region from %0.3f to %0.3f. New map: N=%d,Len=%0.3f(mapid=%d,id=%lld),splits=%d\n",
		     leftend, rightend, pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	      fflush(stdout);
	    }
	  }

	  if(DEBUG && !RELEASE) assert(Xmap->paired && Ymap->paired);

	  continue;

	} // if(overlap >= SplitSegDup && K >= SplitSegDupN && SegDupEnds)
      } // if(BreakEndOutlier ... )

      if(BreakOutlier > 0.0 && SplitSegDup <= 0.0 && PairMergeHmap <= 0 && 
	   overlap >= BreakOutlier && K >= BreakOutlierN && maxOutlier >= BreakOutlierI){
	
	double SplitLen = 0.0;/* largest broken off piece of Ymap (or Xmap) so far */

	// break Ymap at largest internal outlier from IleftOutlier to IrightOutlier

	/* break off left end of Ymap up to IleftOutlier */
	Ymap->paired = 1+nummaps;

	maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	if(DEBUG) assert(PairSplit <= 0.0);
	if(!Cfirst || CfirstPairs==2)
	  XXmap = YYmap = Gmap;
	else {
	  YYmap = &Gmap[0];
	  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	}

	Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	pmap->mapid = nummaps++;
	pmap->fileid = -1;
	pmap->paired = 0;
	pmap->centercloned = 1;
	pmap->origmap = 0;

	if(IDmax < IDMAX)
	  pmap->id = Ymap->id + IDmax;
	else {
          #pragma omp critical
	  {
	    pmap->id = ++IDmax;
	  }
	}
	if(DEBUG && pmap->id <= 0){
	  printf("pmap->id= %lld -> %lld, IDmax= %lld\n",Ymap->id, pmap->id, IDmax);
	  fflush(stdout);
	  assert(pmap->id > 0);
	}

	pmap->origid = Ymap->origid;
	for(int c = 0; c < colors; c++)
	  pmap->Nickase[c] = Ymap->Nickase[c];
	pmap->startloc = pmap->endloc = pmap->len = 0.0;
	    
	pmap->trim(0.0, Y[IleftOutlier] + MININTERVAL);

	SplitLen = pmap->site[0][pmap->numsite[0]+1];
	    
	if(VERB/* HERE >=2 */){
	  printf("    Breaking off left end of Y up to %0.3f due to suspected Chimeric/Misscaled Internal Outlier from %0.3f to %0.3f. New map: N=%d,Len=%0.3f,(mapid=%d,id=%lld),splits=%d\n",
		 Y[IleftOutlier], Y[IleftOutlier],Y[IrightOutlier], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	  fflush(stdout);
	}

	/* break off right end of Ymap beyond Y[IrightOutlier] */
	if(Ymap->paired)
	  splitcnt++;

	maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	if(DEBUG) assert(PairSplit <= 0.0);
	if(!Cfirst || CfirstPairs==2)
	  XXmap = YYmap = Gmap;
	else {
	  YYmap = &Gmap[0];
	  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	}

	pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	pmap->mapid = nummaps++;
	pmap->fileid = -1;
	pmap->paired = 0;
	pmap->centercloned = 1;
	pmap->origmap = 0;

	if(IDmax < IDMAX){
	  pmap->id = Ymap->id + IDmax;
	  if(Ymap->paired)
	    pmap->id += IDmax;
	} else {
          #pragma omp critical
	  {
	    pmap->id = ++IDmax;
	  }
	}
	if(DEBUG && pmap->id <= 0){
	  printf("pmap->id= %lld -> %lld, IDmax= %lld\n",Ymap->id, pmap->id, IDmax);
	  fflush(stdout);
	  assert(pmap->id > 0);
	}

	pmap->origid = Ymap->origid;
	for(int c = 0; c < colors; c++)
	  pmap->Nickase[c] = Ymap->Nickase[c];
	pmap->startloc = pmap->endloc = pmap->len = 0.0;

	if(Y[N+1] - Y[IrightOutlier] > SplitLen)
	  Ymap->paired = nummaps;

	pmap->trim(Y[IrightOutlier] - MININTERVAL, Y[N+1]);

	SplitLen = max(SplitLen, pmap->site[0][pmap->numsite[0]+1]);

	if(VERB/* HERE >=2 */){
	  printf("    Breaking off right end of Y from %0.3f kb to %0.3f due to suspected Chimeric/Scaling Outlier from %0.3f to %0.3f. New map: N=%d,Len=%0.3f(mapid=%d,id=%lld),splits=%d\n",
		 Y[IrightOutlier], Y[N+1], Y[IleftOutlier],Y[IrightOutlier],pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	  fflush(stdout);
	}

	/* Save remainder of Ymap : Y[IleftOutlier .. IrightOutlier] */
	splitcnt++;

	maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	if(DEBUG) assert(PairSplit <= 0.0);
	if(!Cfirst || CfirstPairs==2)
	  XXmap = YYmap = Gmap;
	else {
	  YYmap = &Gmap[0];
	  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	}

	pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	pmap->mapid = nummaps++;
	pmap->fileid = -1;
	pmap->paired = 0;
	pmap->centercloned = 1;
	pmap->origmap = 0;
	pmap->id = Ymap->id;
	pmap->origid = Ymap->origid;
	for(int c = 0; c < colors; c++)
	  pmap->Nickase[c] = Ymap->Nickase[c];
	pmap->startloc = pmap->endloc = pmap->len = 0.0;
	    
	double leftend = Y[IleftOutlier] - MININTERVAL;
	double rightend = Y[IrightOutlier] + MININTERVAL;
	if(rightend - leftend > SplitLen)
	  Ymap->paired = nummaps;

	pmap->trim(leftend, rightend);

	if(VERB/* HERE >=2 */){
	  printf("    Cloning largest Outlier region of Y from %0.3f to %0.3f. New map: N=%d,Len=%0.3f(mapid=%d,id=%lld),splits=%d\n",
		 leftend,rightend, pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	  fflush(stdout);
	}

	// break up Xmap at largest internal outlier from JleftOutlier to JrightOutlier 

	/* break off left end of Xmap up to JleftOutlier */
	Xmap->paired = 1+nummaps;
	    
	maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	if(DEBUG) assert(PairSplit <= 0.0);
	if(!Cfirst || CfirstPairs==2)
	  XXmap = YYmap = Gmap;
	else {
	  YYmap = &Gmap[0];
	  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	}

	pmap = Gmap[nummaps] = new Cmap(Xmap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	pmap->mapid = nummaps++;
	pmap->fileid = -1;
	pmap->paired = 0;
	pmap->centercloned = 1;
	pmap->origmap = 0;

	if(IDmax < IDMAX)
	  pmap->id = Xmap->id + IDmax;
	else {
          #pragma omp critical
	  {
	    pmap->id = ++IDmax;
	  }
	}
	if(DEBUG && pmap->id <= 0){
	  printf("pmap->id= %lld -> %lld, IDmax= %lld\n",Xmap->id, pmap->id, IDmax);
	  fflush(stdout);
	  assert(pmap->id > 0);
	}

	pmap->origid = Xmap->origid;
	for(int c = 0; c < colors; c++)
	  pmap->Nickase[c] = Xmap->Nickase[c];
	pmap->startloc = pmap->endloc = pmap->len = 0.0;
	    
	leftend = p->orientation ? X[M+1- JleftOutlier] - MININTERVAL : 0.0;
	rightend = p->orientation ? X[M+1]                            : X[JleftOutlier] + MININTERVAL;

	pmap->trim(leftend, rightend);

	SplitLen = pmap->site[0][pmap->numsite[0]+1];

	if(VERB/* HERE >=2 */){
	  printf("    Breaking off %s end of X from %0.3f to %0.3f due to suspected Chimeric/Scaling Outlier from %0.3f to %0.3f. New map: N=%d,Len=%0.3f(mapid=%d,id=%lld),splits=%d\n",
		 p->orientation ? "right" : "left", leftend, rightend, p->orientation ? X[M+1-JrightOutlier] : X[JleftOutlier], p->orientation ? X[M+1-JleftOutlier] : X[JrightOutlier], 
		 pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	  fflush(stdout);
	}

	/* break off right end of Xmap beyond Outlier JleftOutlier .. JrightOutlier */
	if(Xmap->paired)
	  splitcnt++;
	    
	maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	if(DEBUG) assert(PairSplit <= 0.0);
	if(!Cfirst || CfirstPairs==2)
	  XXmap = YYmap = Gmap;
	else {
	  YYmap = &Gmap[0];
	  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	}

	pmap = Gmap[nummaps] = new Cmap(Xmap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	pmap->mapid = nummaps++;
	pmap->fileid = -1;
	pmap->paired = 0;
	pmap->centercloned = 1;
	pmap->origmap = 0;

	if(IDmax < IDMAX){
	  pmap->id = Xmap->id + IDmax;
	  if(Xmap->paired)
	    pmap->id += IDmax;
	} else {
          #pragma omp critical
	  {
	    pmap->id = ++IDmax;
	  }
	}
	if(DEBUG && pmap->id <= 0){
	  printf("pmap->id= %lld -> %lld, IDmax= %lld\n",Xmap->id, pmap->id, IDmax);
	  fflush(stdout);
	  assert(pmap->id > 0);
	}

	pmap->origid = Xmap->origid;
	for(int c = 0; c < colors; c++)
	  pmap->Nickase[c] = Xmap->Nickase[c];
	pmap->startloc = pmap->endloc = pmap->len = 0.0;

	if((p->orientation ? X[M+1- JrightOutlier] : X[M+1] - X[JrightOutlier]) > SplitLen)
	  Xmap->paired = nummaps;

	leftend = p->orientation ? 0.0 : X[JrightOutlier] - MININTERVAL;
	rightend = p->orientation ? X[M+1 - JrightOutlier] + MININTERVAL : X[M+1];

	pmap->trim(leftend, rightend);

	SplitLen = max(SplitLen, pmap->site[0][pmap->numsite[0]+1]);

	if(VERB/* HERE >=2 */){
	  printf("    Breaking off %s end of X from %0.3f kb to %0.3f due to suspected Chimeric/Scaling Outlier from %0.3f to %0.3f. New map: N=%d,Len=%0.3f(mapid=%d,id=%lld),splits=%d\n",
		 p->orientation ? "left" : "right", leftend, rightend, p->orientation ? X[M+1-JrightOutlier] : X[JleftOutlier], p->orientation ? X[M+1-JleftOutlier] : X[JrightOutlier],
		 pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	  fflush(stdout);
	}
	
	/* Save Outlier region of Xmap */
	splitcnt++;

	maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	if(DEBUG) assert(PairSplit <= 0.0);
	if(!Cfirst || CfirstPairs==2)
	  XXmap = YYmap = Gmap;
	else {
	  YYmap = &Gmap[0];
	  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	}

	pmap = Gmap[nummaps] = new Cmap(Xmap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	pmap->mapid = nummaps++;
	pmap->fileid = -1;
	pmap->paired = 0;
	pmap->centercloned = 1;
	pmap->origmap = 0;
	pmap->id = Xmap->id;
	pmap->origid = Xmap->origid;
	for(int c = 0; c < colors; c++)
	  pmap->Nickase[c] = Ymap->Nickase[c];
	pmap->startloc = pmap->endloc = pmap->len = 0.0;
	    
	leftend = p->orientation ? X[M+1- JrightOutlier] - MININTERVAL : X[JleftOutlier] - MININTERVAL;
	rightend = p->orientation ? X[M+1- JleftOutlier] + MININTERVAL : X[JrightOutlier] + MININTERVAL;
	if(rightend - leftend > SplitLen)
	  Xmap->paired = nummaps;

	pmap->trim(leftend, rightend);

	if(VERB/* HERE >=2 */){
	  printf("    Cloning largest Outlier region of X from %0.3f to %0.3f. New map: N=%d,Len=%0.3f(mapid=%d,id=%lld),splits=%d\n",
		 leftend, rightend, pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	  fflush(stdout);
	}

	if(DEBUG && !RELEASE) assert(Xmap->paired && Ymap->paired);

	continue;
      } // if(BreakOutlier > 0.0 ... )

      if(SplitSegDup > 0.0 /* WAS && !Cfirst */ && PairMergeHmap <= 0){
	// if(DEBUG && !RELEASE) assert(!merge); // NOTE : can fail if repeat masking is being used, due to a bug that causes the result to change in rare cases if X,Y are reversed

        bool LeftEndSegDup = (leftEnd > SplitSegDupE && leftEndN > SplitSegDupEN);
        bool RightEndSegDup = (rightEnd > SplitSegDupE && rightEndN > SplitSegDupEN);
	double LeftLen = -1.0, RightLen = -1.0;
	int LeftN = -1, RightN = -1;
        int LeftLongerY = -1, RightLongerY = -1;/* If Y end is longer than X end on Left or Right end of alignment */
	if(PAIRMERGE_SEGDUP_ONESIDED){
	  LeftLen = max(Y[LI], (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]));
	  LeftN = (Y[LI] > (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ])) ? LI-1 : LJ-1;
	  LeftLongerY = (Y[LI] >= LeftLen) ? 1 : 0;

	  RightLen = max(Y[N+1] - Y[RI], (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]));
	  RightN = (Y[N+1] - Y[RI] >  (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ])) ? N - RI : M - RJ;
	  RightLongerY = (Y[N+1] - Y[RI] >= RightLen) ? 1 : 0;
	}
	
	size_t YL2 = Ymap->Mask[0] ? Ymap->Mask[0][1] & END_NOEXT2 : 0;
	size_t YR2 = Ymap->Mask[0] ? Ymap->Mask[0][N+1] & END_NOEXT2 : 0;
	size_t XL2 = Xmap->Mask[0] ? Xmap->Mask[0][p->orientation ? M+1 : 1] & END_NOEXT2 : 0;
	size_t XR2 = Xmap->Mask[0] ? Xmap->Mask[0][p->orientation ? 1 : M+1] & END_NOEXT2 : 0;

	bool SegDupEnds = !PAIRMERGE_SEGDUP_ONESIDED ? (LeftEndSegDup && RightEndSegDup) :
	  ((LeftEndSegDup || RightEndSegDup) && min(LeftLen,RightLen) > SplitSegDupE && min(LeftN,RightN) > SplitSegDupEN
	   && (PAIRMERGE_SEGDUP_ONESIDED <= 1 || (!LeftEndSegDup ? (rightEnd >= SplitSegDup_OneSided_MaxLen || (LeftLongerY ? !XL2 : !YL2)) : 
						                   (leftEnd >= SplitSegDup_OneSided_MaxLen || (RightLongerY ? !XR2 : !YR2)))));

	if(VERB>=3){
	  printf("    LeftEndSegDup=%d,RightEndSegDup=%d,LeftLen=%0.3f,%d,RightLen=%0.3f,%d(SplitSegDupE=%0.3f,%d),leftEnd=%0.3f,rightEnd=%0.3f,SplitSegDup_OneSided=%d,%0.2f,LeftLongerY=%d,RightLongerY=%d\n",
		 LeftEndSegDup,RightEndSegDup,LeftLen,LeftN,RightLen,RightN,SplitSegDupE,SplitSegDupEN,leftEnd,rightEnd,PAIRMERGE_SEGDUP_ONESIDED,SplitSegDup_OneSided_MaxLen,LeftLongerY,RightLongerY);
	  fflush(stdout);
	}

	if(SegDupEnds && SegDupMask > 0){/* mark overlap region of both maps with segdup marker */
	  if(VERB/* HERE HERE >=2 */){
	    if(!p->orientation)
	      printf("    Marking SegDup region on Y from %0.3f to %0.3f, and X from %0.3f to %0.3f with Mask |= 0x%lx\n",Y[LI],Y[RI],X[LJ],X[RJ],SegDupMask);
	    else
	      printf("    Marking SegDup region on Y from %0.3f to %0.3f, and X from %0.3f to %0.3f with Mask |= 0x%lx\n",Y[LI],Y[RI],X[M+1-RJ],X[M+1-LJ],SegDupMask);
	    fflush(stdout);
	  }
	  Ymap->centercloned = 1;// mark Ymap as modified, so it is actually output instead of linked to the previous CMAP
	  if(!Ymap->Mask[0]){
	    Ymap->Mask[0] = new size_t[N + 2];
	    memset(Ymap->Mask[0], 0, (N+2) * sizeof(size_t));
	  }
	  for(int I = LI; I <= RI; I++)
	    Ymap->Mask[0][I] |= SegDupMask;

	  Xmap->centercloned = 1;// mark Xmap as modified, so it is actually output instead of linked to the previous CMAP
	  if(!Xmap->Mask[0]){
	    Xmap->Mask[0] = new size_t[M + 2];
	    memset(Xmap->Mask[0], 0, (M+2) * sizeof(size_t));	    
	  }
	  if(p->orientation){
	    for(int J = M+1-RJ; J <= M+1 -LJ; J++)
	      Xmap->Mask[0][J] |= SegDupMask;
	  } else {
	    for(int J = LJ; J <= RJ; J++)
	      Xmap->Mask[0][J] |= SegDupMask;
	  }
	}

	if(VERB/* HERE HERE >=2 */ && overlap >= SplitSegDup){
	  if(PAIRMERGE_SEGDUP_ONESIDED)
	    printf("    Found suspected SegDup on Y from %0.3f to %0.3f : overlap=%0.3f, K=%d, LeftEndSegDup=%d, RightEndSegDup=%d,LeftLen=%0.3f kb,LeftN=%d(Ybig=%d),RightLen=%0.3f kb,RightN=%d(Ybig=%d),XL2=%lu,YL2=%lu,XR2=%lu,YR2=%lu\n",
		   Y[LI],Y[RI], overlap, K, LeftEndSegDup ? 1 : 0, RightEndSegDup ? 1 : 0, LeftLen, LeftN, LeftLongerY, RightLen, RightN, RightLongerY, XL2,YL2,XR2,YR2);
	  else
	    printf("    Found suspected SegDup on Y from %0.3f to %0.3f : overlap=%0.3f, K=%d, LeftEndSegDup=%d, RightEndSegDup=%d\n",
		   Y[LI],Y[RI], overlap, K, LeftEndSegDup ? 1 : 0, RightEndSegDup ? 1 : 0);
	  fflush(stdout);
	}

	if(overlap >= SplitSegDup && K >= SplitSegDupN && SegDupEnds){  /* Assume overlap region is a SegDup and split both maps */

	  if((leftEnd > PairMergeMaxEndKB && leftEndN >= PairMergeMaxEndN) || 
	     (LeftSmallerX && (!PAIRMERGE_SEGDUP_ONESIDED || (LeftLen > PairMergeMaxEndKB && LeftN >= PairMergeMaxEndN)))){  /* left end of Ymap */
	    Ymap->paired = 1+nummaps;
	    
	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;
	    pmap->id = Ymap->id;
	    pmap->origid = Ymap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Ymap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;
	    
	    pmap->trim(0.0, Y[RI]+ MININTERVAL);
	    if(!pmap->Mask[0]){
	      pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
	      memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
	    }
	    pmap->Mask[0][pmap->numsite[0]+1] |= END_NOEXT;

	    if(VERB/* HERE >=2 */){
	      printf("    Cloning left end of Y up to %0.3f due to suspected SegDup from %0.3f to %0.3f and right end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d,id=%lld),splits=%d\n",
		     Y[RI],Y[LI],Y[RI], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	      fflush(stdout);
	    }
	  }

	  if((rightEnd > PairMergeMaxEndKB && rightEndN >= PairMergeMaxEndN) || 
	     (RightSmallerX && (!PAIRMERGE_SEGDUP_ONESIDED || (RightLen > PairMergeMaxEndKB && RightN >= PairMergeMaxEndN)))){/* right end of Ymap */
	    if(Ymap->paired)
	      splitcnt++;

	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;

	    pmap->id = Ymap->id;
	    if(Ymap->paired){
	      if(IDmax < IDMAX)
		pmap->id += IDmax;
	      else {
                #pragma omp critical
	        {
	          pmap->id = ++IDmax;
	        }
	      }

	      if(DEBUG && pmap->id <= 0){
		printf("pmap->id= %lld -> %lld, IDmax= %lld\n",Ymap->id, pmap->id, IDmax);
		fflush(stdout);
		assert(pmap->id > 0);
	      }
	    }

	    pmap->origid = Ymap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Ymap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;

	    if(!Ymap->paired || Y[N+1] - Y[LI] > Y[RI])
	      Ymap->paired = nummaps;

	    pmap->trim(Y[LI] - MININTERVAL, Y[N+1]);
	    if(!pmap->Mask[0]){
	      pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
	      memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
	    }
	    pmap->Mask[0][1] |= END_NOEXT;
	    if(VERB/* HERE >=2 */){
	      printf("    Cloning right end of Y from %0.3f kb to %0.3f due to suspected SegDup from %0.3f to %0.3f and left end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d,id=%lld),splits=%d\n",
		     Y[LI], Y[N+1], Y[LI],Y[RI],pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	      fflush(stdout);
	    }
	  }
	  
	  if((leftEnd > PairMergeMaxEndKB && leftEndN >= PairMergeMaxEndN) || 
	     (LeftSmallerY && (!PAIRMERGE_SEGDUP_ONESIDED || (LeftLen > PairMergeMaxEndKB && LeftN >= PairMergeMaxEndN)))){ /* left end of Xmap */
	    Xmap->paired = 1+nummaps;
	    
	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Xmap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;
	    pmap->id = Xmap->id;
	    pmap->origid = Xmap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Xmap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;
	    
	    double leftend = p->orientation ? max(0.0, X[M+1-RJ] - MININTERVAL) : 0.0;
	    double rightend = p->orientation ? X[M+1] : min(X[M+1], X[RJ] + MININTERVAL);

	    pmap->trim(leftend, rightend);
	    if(!pmap->Mask[0]){
	      pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
	      memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
	    }
	    if(p->orientation)
	      pmap->Mask[0][1] |= END_NOEXT;
	    else
	      pmap->Mask[0][pmap->numsite[0]+1] |= END_NOEXT;

	    if(VERB/* HERE >=2 */){
	      printf("    Cloning %s end of X from %0.3f to %0.3f due to suspected SegDup from %0.3f to %0.3f and %s end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d,id=%lld),splits=%d\n",
		     p->orientation ? "right" : "left", leftend, rightend, p->orientation ? X[M+1-RJ] : X[LJ], p->orientation ? X[M+1-LJ] : X[RJ], p->orientation ? "left" : "right",
		     pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	      fflush(stdout);
	    }
	  }

	  if((rightEnd > PairMergeMaxEndKB && rightEndN >= PairMergeMaxEndN) || 
	     (RightSmallerY && (!PAIRMERGE_SEGDUP_ONESIDED || (RightLen > PairMergeMaxEndKB && RightN >= PairMergeMaxEndN)))){ /* right end of Xmap */
	    if(Xmap->paired)
	      splitcnt++;
	    
	    maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	    /* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	    if(DEBUG) assert(PairSplit <= 0.0);
	    if(!Cfirst || CfirstPairs==2)
	      XXmap = YYmap = Gmap;
	    else {
	      YYmap = &Gmap[0];
	      XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	    }

	    Cmap *pmap = Gmap[nummaps] = new Cmap(Xmap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
	    pmap->mapid = nummaps++;
	    pmap->fileid = -1;
	    pmap->paired = 0;
	    pmap->centercloned = 1;
	    pmap->origmap = 0;

	    pmap->id = Xmap->id;
	    if(Xmap->paired){
	      if(IDmax < IDMAX)
		pmap->id += IDmax;
	      else {
                #pragma omp critical
		{
	          pmap->id = ++IDmax;
	        }
	      }
	      if(DEBUG && pmap->id <= 0){
		printf("pmap->id= %lld -> %lld, IDmax= %lld\n",Xmap->id, pmap->id, IDmax);
		fflush(stdout);
		assert(pmap->id > 0);
	      }
	    }

	    pmap->origid = Xmap->origid;
	    for(int c = 0; c < colors; c++)
	      pmap->Nickase[c] = Xmap->Nickase[c];
	    pmap->startloc = pmap->endloc = pmap->len = 0.0;

	    if(!Xmap->paired || (p->orientation ? (X[M+1-LJ] > X[M+1]-X[M+1-RJ]) : (X[M+1]-X[LJ] > X[RJ])))
	      Xmap->paired = nummaps;

	    double leftend = p->orientation ? 0.0 : max(0.0,X[LJ] - MININTERVAL);
	    double rightend = p->orientation ? min(X[M+1], X[M+1-LJ] + MININTERVAL) : X[M+1];

	    pmap->trim(leftend, rightend);
	    if(!pmap->Mask[0]){
	      pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
	      memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
	    }
	    if(p->orientation)
	      pmap->Mask[0][pmap->numsite[0]+1] |= END_NOEXT;
	    else
	      pmap->Mask[0][1] |= END_NOEXT;
	    if(VERB/* HERE >=2 */){
	      printf("    Cloning %s end of X from %0.3f kb to %0.3f due to suspected SegDup from %0.3f to %0.3f and %s end marked non-extendable. New map: N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d,id=%lld),splits=%d\n",
		     p->orientation ? "left" : "right", leftend, rightend, p->orientation ? X[M+1-RJ] : X[LJ], p->orientation ? X[M+1-LJ] : X[RJ], p->orientation ? "right" : "left",
		     pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid,pmap->id,splitcnt);
	      fflush(stdout);
	    }
	  }

	  if(DEBUG && !RELEASE) assert(Xmap->paired && Ymap->paired);

	  continue;

	} else if(overlap >= SplitSegDupLM && K >= SplitSegDupNM && max(leftEnd,rightEnd) <= SplitSegDupE && 
		  (PairMergeMaxOutlier <= 0.0 || (maxOutlier < PairMergeMaxOutlier && maxOutlierM < PairMergeMaxOutlierM))){/* merge & split */



	  // HERE HERE HERE
	}
      }// if (SplitSegDup > 0.0 && PairMergeHmap <= 0)

      if(DEBUG) assert(!Xmap->paired && !Ymap->paired);

      if(!merge && (PairMergeHmap <= 0 ? (Truncate > 0.0 && overlap > Truncate && K >= TruncateN) : /* check if we should truncate the shorter end of an overlapped map */
		                         (overlap >= TruncateHmap && K >= TruncateNHmap))){/* OR check if Haplotype merge is needed : in both cases ignore outliers smaller than TruncateMaxOutlier */
	int UL = 0;/* locate left alignment region p->sites*[0..UL] that has no outliers larger than TruncateMaxOutlier kb or more than TruncateMaxOutlierM misaligned labels */
	double LeftLen = 0.0;

	if(leftEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB) || leftEndN < PairMergeMaxEndN){
	  /* scan from left end till first outlier larger than TruncateMaxOutlier kb or TruncateMaxOutlierM misaligned labels is encountered */
	  if(DEBUG) assert(K == p->numpairs);
	  int lastI = LI, I;
	  int lastJ = LJ, J;
	  for(int k = 1; k < K;lastI = I, lastJ = J, k++){
	    I = p->sites1[k];
	    J = p->sites2[k];
	    
	    /* outlier lies between sites lastI .. I on reference and lastJ .. J on query */
	    if(DEBUG>=2) assert(J >= lastJ);
	    if(DEBUG>=2) assert(I >= lastI);
	    double qrystartpos = (p->orientation ? (X[M + 1] - X[M + 1 - lastJ]) : X[lastJ]) * Xscale; /* QryStartPos */
	    double qryendpos   = (p->orientation ? (X[M + 1] - X[M + 1 - J]) : X[J])*Xscale; /* QryEndPos */
	    double refstartpos = Y[lastI]*Yscale; /* RefStartPos */
	    double refendpos   = Y[I]*Yscale; /* RefEndPos */
	    if(DEBUG>=2) assert(refendpos >= refstartpos);
	    if(DEBUG) assert(qryendpos >= qrystartpos);// If this fails, fix next line
	    double maxOutlier = fabs(refendpos - refstartpos - qryendpos + qrystartpos);
	    int maxOutlierM = I-lastI + J-lastJ - 2;
	    
	    int outlier = (p->iscore[k] > p->outscore[k] + (RFLOAT)0.01 || (outlierExtend && (lastJ >= J || lastI >= I)) || maxOutlier > TruncateMaxOutlier || maxOutlierM > TruncateMaxOutlierM) ? 1 : 0;

	    if(outlier && PairMergeOutlierChiSq > 0.0 && CmapChimQuality >= 3){// skip if outlier variance is too large or has too many molecule outliers
	      double Outlier = maxOutlier;
	      if(refendpos - refstartpos < qryendpos - qrystartpos){
		float *OutlierFrac = Ymap->OutlierFrac[0];
		float *FragSd = Ymap->FragSd[0];
		float *ExpSd = Ymap->ExpSd[0];
		float *FragCov = Ymap->FragCov[0];
		float *FragChiSq = Ymap->FragChiSq[0];
		if(DEBUG) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);

		float *ChimQuality = Ymap->ChimQuality[0];
		float *ChimNorm = Ymap->ChimNorm[0];

		double LRange = (p->orientation ? X[M+1-lastJ] - X[M+1-J] : X[J] - X[lastJ]) * Xscale;
		double LRangeRatio = LRange / AvgInterval;

		int i = lastI;
		for(; i < I; i++){
		  double SRange = (Y[i+1] - Y[i]) * Yscale;
		  double FragCovi = FragCov[i];
		  if(PairMergeChimNorm > 0.0 && ChimNorm[i] > 0.0 && ChimNorm[i+1] > 0.0)
		    FragCovi = min(FragCovi, min(ChimNorm[i]*ChimQuality[i],ChimNorm[i+1]*ChimQuality[i+1]) * PairMergeChimNorm * 0.01);

		  if(Outlier < PairMergeOutlierMaxSize &&
		     (OutlierFrac[i] > PairMergeOutlierFrac || 
		      (FragCovi < PairMergeOutlierMaxCov && 
		       LRange > PairMergeOutlierMinInterval &&
		       ((FragChiSq[i] < PairMergeOutlierChiSq &&
			 FragSd[i] > ExpSd[i] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
			 FragSd[i] > ExpSd[i] * PairMergeOutlierSdRatio) ||		     
			FragCovi < PairMergeOutlierMinCov + 
			((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != abs(J-lastJ))?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale))))
		    break;
		}
		if(i < I)
		  outlier = 0;// outlier will be ignored
	      } else {
		float *OutlierFrac = Xmap->OutlierFrac[0];
		float *FragSd = Xmap->FragSd[0];
		float *ExpSd = Xmap->ExpSd[0];
		float *FragCov = Xmap->FragCov[0];
		float *FragChiSq = Xmap->FragChiSq[0];
		if(DEBUG) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);

		float *ChimQuality = Xmap->ChimQuality[0];
		float *ChimNorm = Xmap->ChimNorm[0];

		double LRange = (Y[I] - Y[lastI]) * Yscale;
		double LRangeRatio = LRange / AvgInterval;

		int jmin = p->orientation ? M+1-J : lastJ;
		int jmax = p->orientation ? M+1-lastJ : J;
		int j = jmin;
		for(; j < jmax; j++){
		  double SRange = (X[j+1] - X[j]) * Xscale;
		  double FragCovj = FragCov[j];
		  if(PairMergeChimNorm > 0.0 && ChimNorm[j] > 0.0 && ChimNorm[j+1] > 0.0)
		    FragCovj = min(FragCovj, min(ChimNorm[j] * ChimQuality[j], ChimNorm[j+1] * ChimQuality[j+1]) * PairMergeChimNorm * 0.01);

		  if(Outlier < PairMergeOutlierMaxSize &&
		     (OutlierFrac[j] > PairMergeOutlierFrac || 
		      (FragCovj < PairMergeOutlierMaxCov && 
		       LRange > PairMergeOutlierMinInterval &&
		       ((FragChiSq[j] < PairMergeOutlierChiSq &&
			 FragSd[j] > ExpSd[j] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
			 FragSd[j] > ExpSd[j] * PairMergeOutlierSdRatio) ||		     
			FragCovj < PairMergeOutlierMinCov + 
			((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != abs(J-lastJ))?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale))))
		    break;
		}
		if(j < jmax)
		  outlier = 0;// outlier will be ignored
	      }
	    }

	    if(outlier){
	      double YLeft = (Y[lastI+1]-Y[lastI]) * Yscale;
	      double YRight = (Y[I] - Y[I-1]) * Yscale;
	      double XLeft = (p->orientation ? X[M+1-lastJ] - X[M+1-(lastJ+1)] : X[lastJ+1] - X[lastJ]) * Xscale;
	      double XRight = (p->orientation ? X[M+1-(J-1)] - X[M+1-J] : X[J] - X[J-1]) * Xscale;

	      if(I - lastI + J - lastJ - 2 >= PairMergeMaxOutlierBalanced){
		if(I > lastI + 2){
		  double Middle = (Y[I-1] - Y[lastI+1]) * Yscale;
		  double MinOutlier = Middle + max(0.0, max(YLeft-XLeft,YRight-XRight));
		  maxOutlier = max(maxOutlier, MinOutlier);
		}
		if(J > lastJ + 2){
		  double Middle = (p->orientation ? X[M+1-(lastJ+1)] - X[M+1-(J-1)] : X[J-1] - X[lastJ+1]) * Xscale;
		  double MinOutlier = Middle + max(0.0, max(XLeft-YLeft,XRight-YRight));
		  maxOutlier = max(maxOutlier, MinOutlier);
		}
	      }
	      if(maxOutlier > TruncateMaxOutlier || maxOutlierM > TruncateMaxOutlierM){// WAS107 	      if(maxOutlier > TruncateMaxOutlier){
		if(VERB>=3){
		  printf("    Computed LeftLen=%0.3f,UL=%d:k=%d/%d,Y[%d..%d]=%0.3f..%0.3f,X[%d..%d]=%0.3f..%0.3f,iscore[k]=%0.6f,outscore[k]=%0.6f,outlierExtend=%d,outlier=%0.3f kb,%d labels,(LI=%d,LJ=%d)\n",
			 LeftLen,UL,k,K,lastI,I,refstartpos,refendpos,lastJ,J, qrystartpos,qryendpos,p->iscore[k],p->outscore[k],outlierExtend,maxOutlier,maxOutlierM,LI,LJ);
		  fflush(stdout);
		}
		break;
	      }
	    }

	    LeftLen += (refendpos - refstartpos + qryendpos - qrystartpos) * 0.5;
	    UL = k;/* update left side region with no large outlier to p->sites1[0..UL],p->sites2[0..UL] */
	    if(VERB>=4){
	      printf("    Updated LeftLen=%0.3f,UL=%d:k=%d/%d,Y[%d..%d]=%0.3f..%0.3f,X[%d..%d]=%0.3f..%0.3f,iscore[k]=%0.6f,outscore[k]=%0.6f,outlierExtend=%d,outlier=%d\n",
		     LeftLen,UL,k,K,lastI,I,refstartpos,refendpos,lastJ,J, qrystartpos,qryendpos,p->iscore[k],p->outscore[k],outlierExtend,outlier);
	      fflush(stdout);
	    }
	  }
	}
	if(DEBUG) assert(UL < K);

	int UR = K-1;/* locate right alignment region p->sites*[UR .. K-1] that has no outliers larger than TruncateMaxOutlier */ 
	double RightLen = 0.0;

	if(rightEnd <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB) || rightEndN < PairMergeMaxEndN){
	  /* scan from right end till first outlier larger than TruncateMaxOutlier is encountered */
	  if(UL < K-1){
	    UR = K - 1;
	    RightLen = 0.0;

	    int I = RI, lastI;
	    int J = RJ, lastJ;
	    for(int k = K; --k > 0;I = lastI, J = lastJ){
	      lastI = p->sites1[k-1];
	      lastJ = p->sites2[k-1];
	      if(DEBUG>=2) assert(J >= lastJ);
	      if(DEBUG>=2) assert(I >= lastI);
	      
	      /* outlier lies between sites lastI .. I on reference and lastJ .. J on query */
	      double qrystartpos = (p->orientation ? (X[M + 1] - X[M + 1 - lastJ]) : X[lastJ]) * Xscale; /* QryStartPos */
	      double qryendpos   = (p->orientation ? (X[M + 1] - X[M + 1 - J]) : X[J])*Xscale; /* QryEndPos */
	      double refstartpos = Y[lastI]*Yscale; /* RefStartPos */
	      double refendpos   = Y[I]*Yscale; /* RefEndPos */
	      if(DEBUG>=2) assert(refendpos >= refstartpos);
	      if(DEBUG) assert(qryendpos >= qrystartpos);// If this fails, fix next line
	      double maxOutlier = fabs(refendpos - refstartpos - qryendpos + qrystartpos);
	      int maxOutlierM = I-lastI + J-lastJ - 2;// NEW107
	      
	      // WAS107	      int outlier = (p->iscore[k] > p->outscore[k] + (RFLOAT)0.01 || (outlierExtend && (lastJ >= J || lastI >= I))) ? 1 : 0;
	      int outlier = (p->iscore[k] > p->outscore[k] + (RFLOAT)0.01 || (outlierExtend && (lastJ >= J || lastI >= I)) || maxOutlier > TruncateMaxOutlier || maxOutlierM > TruncateMaxOutlierM) ? 1 : 0;

	      if(outlier && PairMergeOutlierChiSq > 0.0 && CmapChimQuality >= 3){// skip if outlier variance is too large or has too many molecule outliers
		double Outlier = maxOutlier;
		if(refendpos - refstartpos < qryendpos - qrystartpos){
		  float *OutlierFrac = Ymap->OutlierFrac[0];
		  float *FragSd = Ymap->FragSd[0];
		  float *ExpSd = Ymap->ExpSd[0];
		  float *FragCov = Ymap->FragCov[0];
		  float *FragChiSq = Ymap->FragChiSq[0];
		  if(DEBUG) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);

		  float *ChimQuality = Ymap->ChimQuality[0];
		  float *ChimNorm = Ymap->ChimNorm[0];

		  double LRange = (p->orientation ? X[M+1-lastJ] - X[M+1-J] : X[J] - X[lastJ]) * Xscale;
		  double LRangeRatio = LRange / AvgInterval;

		  int i = lastI;
		  for(; i < I; i++){
		    double SRange = (Y[i+1] - Y[i]) * Yscale;
		    double FragCovi = FragCov[i];
		    if(PairMergeChimNorm > 0.0 && ChimNorm[i] > 0.0 && ChimNorm[i+1] > 0.0)
		      FragCovi = min(FragCovi, min(ChimNorm[i] * ChimQuality[i], ChimNorm[i+1] * ChimQuality[i+1]) * PairMergeChimNorm * 0.01);

		    if(Outlier < PairMergeOutlierMaxSize &&
		       (OutlierFrac[i] > PairMergeOutlierFrac || 
			(FragCovi < PairMergeOutlierMaxCov && 
			 LRange > PairMergeOutlierMinInterval &&
			 ((FragChiSq[i] < PairMergeOutlierChiSq &&
			   FragSd[i] > ExpSd[i] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
			   FragSd[i] > ExpSd[i] * PairMergeOutlierSdRatio) ||		     
			  FragCovi < PairMergeOutlierMinCov + 
			  ((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != abs(J-lastJ))?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale))))
		      break;
		  }
		  if(i < I)
		    outlier = 0;// outlier will be ignored
		} else {
		  float *OutlierFrac = Xmap->OutlierFrac[0];
		  float *FragSd = Xmap->FragSd[0];
		  float *ExpSd = Xmap->ExpSd[0];
		  float *FragCov = Xmap->FragCov[0];
		  float *FragChiSq = Xmap->FragChiSq[0];
		  if(DEBUG) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);

		  float *ChimQuality = Xmap->ChimQuality[0];
		  float *ChimNorm = Xmap->ChimNorm[0];

		  double LRange = (Y[I] - Y[lastI]) * Yscale;
		  double LRangeRatio = LRange / AvgInterval;

		  int jmin = p->orientation ? M+1-J : lastJ;
		  int jmax = p->orientation ? M+1-lastJ : J;
		  int j = jmin;
		  for(; j < jmax; j++){
		    double SRange = (X[j+1] - X[j]) * Xscale;
		    double FragCovj = FragCov[j];
		    if(PairMergeChimNorm > 0.0 && ChimNorm[j] > 0.0 && ChimNorm[j+1] > 0.0)
		      FragCovj = min(FragCovj, min(ChimNorm[j]*ChimQuality[j],ChimNorm[j+1]*ChimQuality[j+1]) * PairMergeChimNorm * 0.01);

		    if(Outlier < PairMergeOutlierMaxSize &&
		       (OutlierFrac[j] > PairMergeOutlierFrac || 
			(FragCovj < PairMergeOutlierMaxCov && 
			 LRange > PairMergeOutlierMinInterval &&
			 ((FragChiSq[j] < PairMergeOutlierChiSq &&
			   FragSd[j] > ExpSd[j] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
			   FragSd[j] > ExpSd[j] * PairMergeOutlierSdRatio) ||		     
			  FragCovj < PairMergeOutlierMinCov + 
			  ((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != abs(J-lastJ))?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale))))
		      break;
		  }
		  if(j < jmax)
		    outlier = 0;// outlier will be ignored
		}
	      }

	      if(outlier){
		double YLeft = (Y[lastI+1]-Y[lastI]) * Yscale;
		double YRight = (Y[I] - Y[I-1]) * Yscale;
		double XLeft = (p->orientation ? X[M+1-lastJ] - X[M+1-(lastJ+1)] : X[lastJ+1] - X[lastJ]) * Xscale;
		double XRight = (p->orientation ? X[M+1-(J-1)] - X[M+1-J] : X[J] - X[J-1]) * Xscale;

		if(I - lastI + J - lastJ - 2 >= PairMergeMaxOutlierBalanced){
		  if(I > lastI + 2){
		    double Middle = (Y[I-1] - Y[lastI+1]) * Yscale;
		    double MinOutlier = Middle + max(0.0, max(YLeft-XLeft,YRight-XRight));
		    maxOutlier = max(maxOutlier, MinOutlier);
		  }
		  if(J > lastJ + 2){
		    double Middle = (p->orientation ? X[M+1-(lastJ+1)] - X[M+1-(J-1)] : X[J-1] - X[lastJ+1]) * Xscale;
		    double MinOutlier = Middle + max(0.0, max(XLeft-YLeft,XRight-YRight));
		    maxOutlier = max(maxOutlier, MinOutlier);
		  }
		}
		if(maxOutlier > TruncateMaxOutlier || maxOutlierM > TruncateMaxOutlierM){// WAS107		if(maxOutlier > TruncateMaxOutlier){
		  if(VERB>=3){
		    printf("    Computed RightLen=%0.3f:k=%d/%d,Y[%d..%d]=%0.3f..%0.3f,X[%d..%d]=%0.3f..%0.3f,iscore[k]=%0.6f,outscore[k]=%0.6f,outlierExtend=%d,outlier=%0.3f kb,%d labels\n",
			   RightLen,k,K,lastI,I,refstartpos,refendpos,lastJ,J,qrystartpos,qryendpos,p->iscore[k],p->outscore[k],outlierExtend,maxOutlier,maxOutlierM);
		    fflush(stdout);
		  }
		  break;
		}
	      }
	      RightLen += (refendpos - refstartpos + qryendpos - qrystartpos) * 0.5;
	      UR = k-1;/* update right side region with no large outlier to p->sites1[UR..K-1],p->sites2[UR..K-1] */
	      if(VERB>=4){
		printf("    Updated RightLen=%0.3f:k=%d/%d,Y[%d..%d]=%0.3f..%0.3f,X[%d..%d]=%0.3f..%0.3f,iscore[k]=%0.6f,outscore[k]=%0.6f,outlierExtend=%d,outlier=%d\n",
		       RightLen,k,K,lastI,I,refstartpos,refendpos,lastJ,J,qrystartpos,qryendpos,p->iscore[k],p->outscore[k],outlierExtend,outlier);
		fflush(stdout);
	      }
	    }
	  } else {
	    UR = 0;
	    RightLen = LeftLen;
	  }
	}

	if(VERB/* HERE >=2 */){
	  printf("    -Truncate %0.1f %d %0.2f %d: leftEnd=%0.3f(kb),LeftLen=%0.3f(kb),UL=0..%d(Y[%d..%d],X[%d..%d]); rightEnd=%0.3f(kb),RightLen=%0.3f(kb),UR=%d..%d(Y[%d..%d],X[%d..%d])\n",
		 Truncate,TruncateN,TruncateMaxOutlier,TruncateFix, leftEnd,LeftLen,UL,p->sites1[0],p->sites1[UL],p->sites2[0],p->sites2[UL],
		 rightEnd,RightLen,UR,K-1,p->sites1[UR],p->sites1[K-1],p->sites2[UR],p->sites2[K-1]);
	  fflush(stdout);
	}

	if(PairMergeHmap <= 0 && Truncate > 0.0 && LeftLen > Truncate && UL+1 > TruncateN){ /* try to trim left end of alignment so overlap with no outliers larger than TruncateMaxOutlier is reduced to Truncate */
	  /* scan from left end till LeftLen has been reduced to smallest size possible */
	  int lastI = LI, I;
	  int lastJ = LJ, J;
	  int kL = 1;
	  for(; kL < UL; lastI = I, lastJ = J, kL++){
	    I = p->sites1[kL];
	    J = p->sites2[kL];
	    /* interval lies between sites lastI .. I on reference and lastJ .. J on query */
	    double qrystartpos = (p->orientation ? (X[M + 1] - X[M + 1 - lastJ]) : X[lastJ]) * Xscale; /* QryStartPos */
	    double qryendpos   = (p->orientation ? (X[M + 1] - X[M + 1 - J]) : X[J])*Xscale; /* QryEndPos */
	    double refstartpos = Y[lastI]*Yscale; /* RefStartPos */
	    double refendpos   = Y[I]*Yscale; /* RefEndPos */
	      
	    double delta = (refendpos - refstartpos + qryendpos - qrystartpos) * 0.5;
	    if(LeftLen - delta < Truncate || UL-kL < TruncateN){
	      kL--;/* The truncated aligned region is p->sites1[kL..UL] */
	      break;
	    }
	    LeftLen -= delta;
	  }

	  if(VERB/* HERE >=2 */ && kL > 0){
	    printf("    LeftLen reduced to %0.3f(kb),kl=%d..UL=%d(Y[%d..%d],X[%d..%d])\n",LeftLen,kL,UL,p->sites1[kL],p->sites1[UL],
		   p->orientation ? M+1 - p->sites2[kL] : p->sites2[kL], p->orientation ? M+1 - p->sites2[UL] : p->sites2[UL]);
	    fflush(stdout);
	  }

	  if(kL > 0){
	    if(!(TrimmedNoMerge && Ymap->Mask[0] && (Ymap->Mask[0][1] & END_NOEXT)) && Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ])){
	      if(!(TRUNCATE_FIX && Y[N+1]-Y[RI] > (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]))){  /* truncate map Y[0 .. p->sites1[kL]] and keep Y[p->sites1[kL] .. N+1] */
		Ymap->paired = 1+nummaps;
		int L = p->sites1[kL];
		if(DEBUG) assert(L > 1);/* since this is not the first aligned label of Y */
		double leftend = max(0.0,Y[L] - MININTERVAL);
		double rightend = Y[N+1];
		
		maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
		/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
		if(DEBUG) assert(PairSplit <= 0.0);
		if(!Cfirst || CfirstPairs==2)
		  XXmap = YYmap = Gmap;
		else {
		  YYmap = &Gmap[0];
		  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
		}

		Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
		pmap->mapid = nummaps++;
		pmap->fileid = -1;
		pmap->paired = 0;
		pmap->centercloned = 1;
		pmap->origmap = 0;
		pmap->id = Ymap->id;
		pmap->origid = Ymap->origid;
		for(int c = 0; c < colors; c++)
		  pmap->Nickase[c] = Ymap->Nickase[c];
		pmap->startloc = pmap->endloc = pmap->len = 0.0;

		pmap->trim(leftend, rightend);
		if(!pmap->Mask[0]){
		  pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		  memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
		}
		pmap->Mask[0][1] |= END_NOEXT;

		if(VERB/* HERE >=2 */){
		  printf("    Trimming off %d labels from left end of Y up to %0.3f and marking end as non-extendable (anchor is orig Y[%d..%d]=%0.3f..%0.3f): N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
			 L-1, leftend, L, p->sites1[UL], Y[L],Y[p->sites1[UL]], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1], pmap->mapid);
		  fflush(stdout);
		}

		if(DEBUG>=2){/* make sure no intermediate labels have non-zero Mask value */
		  for(int i = 2; i <= pmap->numsite[0];i++){
		    if(pmap->Mask[0][i]){
		      printf("WARNING: pmap->Mask[0][%d/%d] = %lx (should be 0)\n",i,pmap->numsite[0],pmap->Mask[0][i]);
		      fflush(stdout);
		      assert(pmap->Mask[0][i] == (size_t)0);
		    }
		  }
		}

		continue;
	      }
	    } else if(p->orientation ? (!(TrimmedNoMerge && Xmap->Mask[0] && (Xmap->Mask[0][M+1] & END_NOEXT)) && Y[LI] > X[M+1] - X[M+1-LJ]) 
		      : (!(TrimmedNoMerge && Xmap->Mask[0] && (Xmap->Mask[0][1] & END_NOEXT)) && Y[LI] > X[LJ])){
	      if(!(TRUNCATE_FIX && Y[N+1]-Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]))){  /* truncate map X[0 .. p->sites2[kL]]  and keep X[p->sites2[kL] .. M+1]*/
		Xmap->paired = 1+nummaps;
		int L = p->sites2[kL];
		if(DEBUG) assert(L > 1);/* since this is not the first aligned label of X */

		double leftend = p->orientation ? 0.0 : max(0.0,X[L] - MININTERVAL);
		double rightend = p->orientation ? min(X[M+1], X[M+1-L] + MININTERVAL) : X[M+1];

		maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
		/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
		if(DEBUG) assert(PairSplit <= 0.0);
		if(!Cfirst || CfirstPairs==2)
		  XXmap = YYmap = Gmap;
		else {
		  YYmap = &Gmap[0];
		  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
		}
		if(DEBUG) assert(Ymap == YYmap[p->mapid1]);
		if(DEBUG) assert(Xmap == XXmap[p->mapid2]);		

		Cmap *pmap = Gmap[nummaps] = new Cmap(Xmap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
		pmap->mapid = nummaps++;
		pmap->fileid = -1;
		pmap->paired = 0;
		pmap->centercloned = 1;
		pmap->origmap = 0;
		pmap->id = Xmap->id;
		pmap->origid = Xmap->origid;
		for(int c = 0; c < colors; c++)
		  pmap->Nickase[c] = Ymap->Nickase[c];
		pmap->startloc = pmap->endloc = pmap->len = 0.0;

		pmap->trim(leftend, rightend);
		if(!pmap->Mask[0]){
		  pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		  memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
		}
		if(p->orientation)
		  pmap->Mask[0][pmap->numsite[0]+1] |= END_NOEXT;
		else
		  pmap->Mask[0][1] |= END_NOEXT;

		if(VERB/* HERE >=2 */){
		  if(p->orientation)
		    printf("    Trimming off %d labels up to %0.3f from right end of X and marking end as non-extendable (anchor is orig X[%d..%d]=%0.3f..%0.3f): M=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
			   L-1, rightend, M+1-p->sites2[UL], M+1-L, X[M+1-p->sites2[UL]], X[M+1-L], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],
			   pmap->Mask[0][pmap->numsite[0]+1],pmap->mapid);
		  else
		    printf("    Trimming off %d labels up to %0.3f from left end of X and marking end as non-extendable (anchor is orig X[%d..%d]=%0.3f..%0.3f): M=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
			   L-1, leftend, L, p->sites2[UL], X[L], X[p->sites2[UL]], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1],pmap->mapid);
		  fflush(stdout);
		}
		
		if(DEBUG>=2){/* make sure no intermediate labels have non-zero Mask value */
		  for(int i = 2; i <= pmap->numsite[0];i++){
		    if(pmap->Mask[0][i]){
		      printf("WARNING: pmap->Mask[0][%d/%d] = %lx (should be 0)\n",i,pmap->numsite[0],pmap->Mask[0][i]);
		      fflush(stdout);
		      assert(pmap->Mask[0][i] == (size_t)0);
		    }
		  }
		}

		continue;
	      }
	    }
	  }
	}

	if(PairMergeHmap <= 0 && Truncate > 0.0/* NEW440 */ && RightLen > Truncate && K-UR > TruncateN){  /* try to trim right end of alignment so overlap with no outliers larger than TruncateMaxOutlier is reduced to Truncate */
	  /* scan from right end till RightLen has been reduced to smallest size possible */
	  int I = RI, lastI;
	  int J = RJ, lastJ;
	  int kR = K;
	  for(; --kR > UR;I = lastI, J = lastJ){
	    lastI = p->sites1[kR-1];
	    lastJ = p->sites2[kR-1];

	    /* interval lies between sites lastI .. I on reference and lastJ .. J on query */
	    double qrystartpos = (p->orientation ? (X[M + 1] - X[M + 1 - lastJ]) : X[lastJ]) * Xscale; /* QryStartPos */
	    double qryendpos   = (p->orientation ? (X[M + 1] - X[M + 1 - J]) : X[J])*Xscale; /* QryEndPos */
	    double refstartpos = Y[lastI]*Yscale; /* RefStartPos */
	    double refendpos   = Y[I]*Yscale; /* RefEndPos */
		
	    double delta = (refendpos - refstartpos + qryendpos - qrystartpos) * 0.5;		
	    if(RightLen - delta < Truncate || kR-UR < TruncateN)
	      break;/* The truncated aligned region is p->sites1[UR..kR] */
	    RightLen -= delta;
	  }
	  if(VERB/* HERE >=2 */ && kR < K-1){
	    printf("    RightLen reduced to %0.3f(kb),UR=%d..kR=%d(Y[%d..%d],X[%d..%d])\n",RightLen,UR,kR,p->sites1[kR],p->sites1[UR],
		   p->orientation ? M+1-p->sites2[kR]:p->sites2[kR],p->orientation ? M+1-p->sites2[UR] : p->sites2[UR]);
	    fflush(stdout);
	  }

	  if(kR < K-1){
	    if(!(TrimmedNoMerge && Ymap->Mask[0] && (Ymap->Mask[0][N+1] & END_NOEXT)) && Y[N+1] - Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ])){
	      if(!(TRUNCATE_FIX && Y[LI] > (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]))){ /* truncate Ymap and keep Y[0 .. p->sites1[kR]] */
		Ymap->paired = 1+nummaps;
		int R = p->sites1[kR];
		if(DEBUG) assert(R < N);/* since this is not the last aligned label of Y */
		double leftend = 0.0;
		double rightend = min(Y[N+1], Y[R] + MININTERVAL);
		
		maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
		/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
		if(DEBUG) assert(PairSplit <= 0.0);
		if(!Cfirst || CfirstPairs==2)
		  XXmap = YYmap = Gmap;
		else {
		  YYmap = &Gmap[0];
		  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
		}

		Cmap *pmap = Gmap[nummaps] = new Cmap(Ymap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this memory
		pmap->mapid = nummaps++;
		pmap->fileid = -1;
		pmap->paired = 0;
		pmap->centercloned = 1;
		pmap->origmap = 0;
		pmap->id = Ymap->id;
		pmap->origid = Ymap->origid;
		for(int c = 0; c < colors; c++)
		  pmap->Nickase[c] = Ymap->Nickase[c];
		pmap->startloc = pmap->endloc = pmap->len = 0.0;

		pmap->trim(leftend, rightend);
		if(!pmap->Mask[0]){
		  pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		  memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
		}
		pmap->Mask[0][pmap->numsite[0] + 1] |= END_NOEXT;

		if(VERB/* HERE >=2 */){
		  printf("  Trimming off %d labels from right end of Y up to %0.3f and marking end as non-extendable (anchor is orig Y[%d..%d]=%0.3f..%0.3f): N=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
			 N-R, rightend, p->sites1[UR], R, Y[p->sites1[UR]], Y[R], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1],pmap->mapid);
		  fflush(stdout);
		}
		
		if(DEBUG>=2){/* make sure no intermediate labels have non-zero Mask value */
		  for(int i = 2; i <= pmap->numsite[0];i++){
		    if(pmap->Mask[0][i]){
		      printf("WARNING: pmap->Mask[0][%d/%d] = %lx (should be 0)\n",i,pmap->numsite[0],pmap->Mask[0][i]);
		      fflush(stdout);
		      assert(pmap->Mask[0][i] == (size_t)0);
		    }
		  }
		}

		continue;
	      }
	    } else if(p->orientation ? (!(TrimmedNoMerge && Xmap->Mask[0] && (Xmap->Mask[0][1] & END_NOEXT)) && Y[N+1] - Y[RI] > X[M+1-RJ])
		      : (!(TrimmedNoMerge && Xmap->Mask[0] && (Xmap->Mask[0][M+1] & END_NOEXT)) && Y[N+1] - Y[RI] > X[M+1] - X[RJ])){
	      if(!(TRUNCATE_FIX && Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]))){ /* truncate Xmap and keep X[0 .. p->sites2[kR]] */
		Xmap->paired = 1+nummaps;
		int R = p->sites2[kR];
		if(DEBUG) assert(R < M);/* since this is not the last aligned label of X */
		double leftend = p->orientation ? max(0.0, X[M+1-R] - MININTERVAL) : 0.0;
		double rightend = p->orientation ? X[M+1] : min(X[M+1], X[R] + MININTERVAL);
		
		maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
		/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
		if(DEBUG) assert(PairSplit <= 0.0);
		if(!Cfirst || CfirstPairs==2)
		  XXmap = YYmap = Gmap;
		else {
		  YYmap = &Gmap[0];
		  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
		}

		Cmap *pmap = Gmap[nummaps] = new Cmap(Xmap);// NOTE this will waste the space allocated by maxmapalloc() for the Cmaps and leak this newly allocated memory
		pmap->mapid = nummaps++;
		pmap->fileid = -1;
		pmap->paired = 0;
		pmap->centercloned = 1;
		pmap->origmap = 0;
		pmap->id = Xmap->id;
		pmap->origid = Xmap->origid;
		for(int c = 0; c < colors; c++)
		  pmap->Nickase[c] = Ymap->Nickase[c];
		pmap->startloc = pmap->endloc = pmap->len = 0.0;

		pmap->trim(leftend, rightend);
		if(!pmap->Mask[0]){
		  pmap->Mask[0] = new size_t[pmap->numsite[0] + 2];
		  memset(pmap->Mask[0], 0, (pmap->numsite[0] + 2) * sizeof(size_t));
		}
		if(p->orientation)
		  pmap->Mask[0][1] |= END_NOEXT;
		else
		  pmap->Mask[0][pmap->numsite[0]+1] |= END_NOEXT;

		if(VERB/* HERE >=2 */){
		  if(p->orientation)
		    printf("   Trimming off %d labels up to %0.3f from left end of X and marking end as non-extendable (anchor is orig X[%d..%d]=%0.3f..%0.3f): M=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
			   M-R, leftend, M+1-R, M+1-p->sites2[UR], X[M+1-R], X[M+1-p->sites2[UR]], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1],pmap->mapid);
		  else
		    printf("   Trimming off %d labels starting at %0.3f from right end of X and marking end as non-extendable (anchor region= X[%d..%d]=%0.3f..%0.3f): M=%d,Len=%0.3f,Mask=%lx,%lx(mapid=%d)\n",
			   M-R, rightend, p->sites2[UR], R, X[p->sites2[UR]], X[R], pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1], pmap->Mask[0][1],pmap->Mask[0][pmap->numsite[0]+1],pmap->mapid);
		  fflush(stdout);
		}

		if(DEBUG>=2){/* make sure no intermediate labels have non-zero Mask value */
		  for(int i = 2; i <= pmap->numsite[0];i++){
		    if(pmap->Mask[0][i]){
		      printf("WARNING: pmap->Mask[0][%d/%d] = %lx (should be 0)\n",i,pmap->numsite[0],pmap->Mask[0][i]);
		      fflush(stdout);
		      assert(pmap->Mask[0][i] == (size_t)0);
		    }
		  }
		}

		continue;
	      }
	    } // else check "right" end of Xmap
	  } // if(kR < K-1) : right end of alignment was trimmed
	} // if(PairMergeHmap <= 0 && RightLen > Truncate && K-UR > TruncateN)

	if(PairMergeHmap > 0 && LeftLen >= TruncateHmap && UL+1 > TruncateNHmap && RightLen >= TruncateHmap && K-UR > TruncateNHmap /* WAS && maxOutlier > PairMergeMaxOutlier */){
	  if(DEBUG) assert(maxOutlier > PairMergeMaxOutlier || maxOutlierM > PairMergeMaxOutlierM || overlap <= PairMerge || (TrimmedNoMerge && overlap < TrimmedNoMerge));// otherwise it should have merged

	  /* check if one map is NOT completely overlapped by the other */
	  if(!Overlapped){
	    if(SplitSegDup > 0.0 && max(leftEnd,rightEnd) > max(SplitSegDupE, PairMergeMaxEndKB)){// NEW3
	      if(VERB/* HERE >=2 */){
		printf("\t Skipping due to end overlap with endoutliers= %0.3f, %0.3f kb (must not exceed %0.3f)\n",
		       leftEnd, rightEnd, PairMergeMaxEndKB);
		fflush(stdout);
	      }
	      continue;
	    }

	    if(TrimmedNoMerge && SplitSegDup > 0.0){/* check if -TrimmedNoMerge blocks merging */
	      /* Left end of Y is marked and overlaps "left" end of X */ 
	      int Yleft = (Ymap->Mask[0] && (Ymap->Mask[0][1] & END_NOEXT) && Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ])) ? 1 : 0;

	      /* "Left" end of X is marked and overlaps left end of Y */
	      int Xleft = (Xmap->Mask[0] && (Xmap->Mask[0][p->orientation ? M+1 : 1] & END_NOEXT) && Y[LI] > (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ])) ? 1 : 0;

	      /* Right end of Y is marked and overlaps "right" end of X */
	      int Yright = (Ymap->Mask[0] && (Ymap->Mask[0][N+1] & END_NOEXT) && Y[N+1]-Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ])) ? 1 : 0;
	
	      /* "Right" end of X is marked and overlaps right end of Y */
	      int Xright = (Xmap->Mask[0] && (Xmap->Mask[0][p->orientation ? 1 : M+1] & END_NOEXT) && Y[N+1]-Y[RI] > (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ])) ? 1 : 0;

	      if((Yleft || Xleft || Yright || Xright) && overlap < TrimmedNoMerge ){
		int LeftSmallerY = (Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]) + TrimmedNoMergeMaxEnd) ? 1 : 0;
		int LeftSmallerX = (Y[LI] > (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]) - TrimmedNoMergeMaxEnd) ? 1 : 0;
		int RightSmallerY = (Y[N+1] - Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]) + TrimmedNoMergeMaxEnd) ? 1 : 0;
		int RightSmallerX = (Y[N+1] - Y[RI] > (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]) - TrimmedNoMergeMaxEnd) ? 1 : 0;
		if(!((LeftSmallerY && RightSmallerY)||(LeftSmallerX && RightSmallerX))){/* smaller map is NOT fully overlapped by a larger map (even allowing for TrimmedNoMergeMaxEnd wiggle room) */
		  if(VERB/* HERE >=2 */){
		    printf("   Skipped merging Ymap->id=%lld,Xmap->id=%lld due to TrimmedNoMerge = %0.3f, overlap=%0.3f,%d, Xleft=%d,Xright=%d,Yleft=%d,Yright=%d,SmallerY:L=%d,R=%d\n",
			   Ymap->id,Xmap->id,TrimmedNoMerge,overlap,p->numpairs,Xleft,Xright,Yleft,Yright,LeftSmallerY,RightSmallerY);
		    fflush(stdout);
		  }
		  continue;// avoid merging end overlaps with ends marked by SplitSegDup to avoid undoing effect of SplitSegDup and re-forming chimeric contigs
		}
	      }
	    }

	    if(LeftSmallerX && !LeftSmallerY && RightSmallerY && !RightSmallerX){
	      if(VERB/* HERE >=2 */){
		printf("End overlap of Ymap->id=%lld with Xmap->id=%lld : right end of Ymap overlaps %s end of Xmap\n",Ymap->id,Xmap->id,p->orientation ? "right":"left");
		fflush(stdout);
	      }
	      if(overlap < PairMergeEndHmap || min(Y[N+1],X[M+1]) <  MinMapLenHmap){
		if(VERB/* HERE >=2 */){
		  printf("\t Skipping due to overlap= %0.3f (must be at least %0.3f) and Ylen=%0.3f,Xlen=%0.3f (must be at least %0.3f)\n",
		       overlap, PairMergeEndHmap, Y[N+1],X[M+1],MinMapLenHmap);
		  fflush(stdout);
		}
		continue;
	      }

	      Calign *q = Ymap->rightoverlap;
	      if(q){/* check if existing overlap is better (larger than Xmap) */
		Cmap *Ymap2 = YYmap[q->mapid1];
		Cmap *Xmap2 = XXmap[q->mapid2];
		int N2 = Ymap2->numsite[0];
		int M2 = Xmap2->numsite[0];
		if(DEBUG) assert(Ymap == Ymap2 || Ymap == Xmap2);
		if((Ymap == Ymap2 ? Xmap2->site[0][M2+1] : Ymap2->site[0][N2+1]) > X[M+1]){
		  if(VERB/* HERE >=2 */){
		    if(Ymap==Ymap2)
		      printf("\t right end of Ymap already overlapped with id=%lld,len=%0.3f (q= %p)\n",Xmap2->id,Xmap2->site[0][M2+1],q);
		    else
		      printf("\t right end of Ymap already overlapped with id=%lld,len=%0.3f (q= %p)\n",Ymap2->id,Ymap2->site[0][N2+1],q);
		    fflush(stdout);
		  }
		  continue;
		}
	      }
	      Calign *r = p->orientation ? Xmap->rightoverlap : Xmap->leftoverlap;
	      if(r){/* check if existing overlap is better (larger than Ymap)  */
		Cmap *Ymap2 = YYmap[r->mapid1];
		Cmap *Xmap2 = XXmap[r->mapid2];
		int N2 = Ymap2->numsite[0];
		int M2 = Xmap2->numsite[0];
		if(DEBUG) assert(Xmap == Ymap2 || Xmap == Xmap2);
		if((Xmap == Ymap2 ? Xmap2->site[0][M2+1] : Ymap2->site[0][N2+1]) > Y[N+1]){
		  if(VERB/* HERE >=2 */){
		    if(Xmap==Ymap2)
		      printf("\t %s end of Xmap already overlapped with id=%lld,len=%0.3f (r= %p)\n",p->orientation ? "right" : "left", Xmap2->id,Xmap2->site[0][M2+1], r);
		    else
		      printf("\t %s end of Xmap already overlapped with id=%lld,len=%0.3f (r= %p)\n",p->orientation ? "right" : "left", Ymap2->id,Ymap2->site[0][N2+1], r);
		    fflush(stdout);
		  }
		  continue;
		}

		/* clean up previous overlap alignment pointers in Xmap */
		if(Xmap == Ymap2){
		  if(Xmap2->leftoverlap == r){
		    if(VERB/* HERE >=2 */){
		      printf("\t Xmap2->id= %lld,len=%0.3f: r=leftoverlap= %p -> 0\n",Xmap2->id,Xmap2->site[0][M2+1], Xmap2->leftoverlap);
		      fflush(stdout);
		    }
		    Xmap2->leftoverlap = NULL;
		  } else {
		    if(VERB/* HERE >=2 */){
		      printf("\t Xmap2->id= %lld,len=%0.3f: r=rightoverlap= %p -> 0\n",Xmap2->id,Xmap2->site[0][M2+1], Xmap2->rightoverlap);
		      fflush(stdout);
		    }
		    if(DEBUG) assert(Xmap2->rightoverlap == r);
		    Xmap2->rightoverlap = NULL;
		  }
		} else {
		  if(Ymap2->leftoverlap == r){
		    Ymap2->leftoverlap = NULL;
		    if(VERB/* HERE >=2 */){
		      printf("\t Ymap2->id= %lld,len=%0.3f: r=leftoverlap= %p -> 0\n",Ymap2->id,Ymap2->site[0][N2+1], Ymap2->leftoverlap);
		      fflush(stdout);
		    }
		  } else {
		    if(VERB/* HERE >=2 */){
		      printf("\t Ymap2->id= %lld,len=%0.3f: r=rightoverlap= %p -> 0\n",Ymap2->id,Ymap2->site[0][N2+1], Ymap2->rightoverlap);
		      fflush(stdout);
		    }
		    if(DEBUG) assert(Ymap2->rightoverlap == r);
		    Ymap2->rightoverlap = NULL;
		  }
		}
	      }

	      if(q){ /* clean up previous overlap alignment pointers in Ymap */
		Cmap *Ymap2 = YYmap[q->mapid1];
		Cmap *Xmap2 = XXmap[q->mapid2];
		int N2 = Ymap2->numsite[0];
		int M2 = Xmap2->numsite[0];
		if(DEBUG) assert(Ymap == Ymap2 || Ymap == Xmap2);
		if(Ymap == Ymap2){
		  if(Xmap2->leftoverlap == q){
		    if(VERB/* HERE >=2 */){
		      printf("\t Xmap2->id= %lld,len=%0.3f: q=leftoverlap= %p -> 0\n", Xmap2->id,Xmap2->site[0][M2+1], Xmap2->leftoverlap);
		      fflush(stdout);
		    }
		    Xmap2->leftoverlap = NULL;
		  } else {
		    if(VERB/* HERE >=2 */){
		      printf("\t Xmap2->id= %lld,len=%0.3f: q=rightoverlap= %p -> 0\n",Xmap2->id,Xmap2->site[0][M2+1], Xmap2->rightoverlap);
		      fflush(stdout);
		    }
		    if(DEBUG) assert(Xmap2->rightoverlap == q);
		    Xmap2->rightoverlap = NULL;
		  }
		} else {
		  if(Ymap2->leftoverlap == q){
		    if(VERB/* HERE >=2 */){
		      printf("\t Ymap2->id= %lld,len=%0.3f: q=leftoverlap= %p -> 0\n",Ymap2->id,Ymap2->site[0][N2+1], Ymap2->leftoverlap);
		      fflush(stdout);
		    }
		    Ymap2->leftoverlap = NULL;
		  } else {
		    if(VERB/* HERE >=2 */){
		      printf("\t Ymap2->id= %lld,len=%0.3f: q=rightoverlap= %p -> 0\n",Ymap2->id,Ymap2->site[0][N2+1], Ymap2->rightoverlap);
		      fflush(stdout);
		    }
		    if(DEBUG) assert(Ymap2->rightoverlap == q);
		    Ymap2->rightoverlap = NULL;
		  }
		}
	      }
	      
	      Ymap->rightoverlap = p;
	      if(p->orientation)
		Xmap->rightoverlap = p;
	      else
		Xmap->leftoverlap = p;
	      if(VERB/* HERE >=2*/){
		printf("Postponing merging Ymap->id=%lld with Xmap->id=%lld : right end of Ymap overlaps %s end of Xmap (q= %p -> %p, r= %p -> %p),Xmap->contig=%p,Ymap->contig=%p\n",
		       Ymap->id,Xmap->id,p->orientation ? "right":"left",q,p,r,p,Xmap->contig,Ymap->contig);
		fflush(stdout);
	      }
	    }

	    if(LeftSmallerY && !LeftSmallerX && RightSmallerX && !RightSmallerY){
	      if(VERB/* HERE >=2*/){
		printf("EndOverlap of Ymap->id=%lld with Xmap->id=%lld : left end of Ymap overlaps %s end of Xmap\n",Ymap->id,Xmap->id,p->orientation ? "left":"right");
		fflush(stdout);
	      }

	      Calign *q = Ymap->leftoverlap;
	      if(q){/* check if existing overlap is better (larger than Xmap) */
		Cmap *Ymap2 = YYmap[q->mapid1];
		Cmap *Xmap2 = XXmap[q->mapid2];
		int N2 = Ymap2->numsite[0];
		int M2 = Xmap2->numsite[0];
		if(DEBUG) assert(Ymap == Ymap2 || Ymap == Xmap2);
		if((Ymap == Ymap2 ? Xmap2->site[0][M2+1] : Ymap2->site[0][N2+1]) > X[M+1]){
		  if(VERB/* HERE >=2 */){
		    if(Ymap==Ymap2)
		      printf("\t left end of Ymap already overlapped with id=%lld,len=%0.3f (q= %p)\n",Xmap2->id,Xmap2->site[0][M2+1],q);
		    else
		      printf("\t left end of Ymap already overlapped with id=%lld,len=%0.3f (q= %p)\n",Ymap2->id,Ymap2->site[0][N2+1],q);
		    fflush(stdout);
		  }
		  continue;
		}
	      }
	      Calign *r = p->orientation ? Xmap->leftoverlap : Xmap->rightoverlap;
	      if(r){/* check if existing overlap is better (larger than Ymap)  */
		Cmap *Ymap2 = YYmap[r->mapid1];
		Cmap *Xmap2 = XXmap[r->mapid2];
		int N2 = Ymap2->numsite[0];
		int M2 = Xmap2->numsite[0];
		if(DEBUG) assert(Xmap == Ymap2 || Xmap == Xmap2);
		if((Xmap == Ymap2 ? Xmap2->site[0][M2+1] : Ymap2->site[0][N2+1]) > Y[N+1]){
		  if(VERB/* HERE >=2 */){
		    if(Xmap==Ymap2)
		      printf("\t %s end of Xmap already overlapped with id=%lld,len=%0.3f (r= %p)\n",p->orientation ? "left" : "right", Xmap2->id,Xmap2->site[0][M2+1], r);
		    else
		      printf("\t %s end of Xmap already overlapped with id=%lld,len=%0.3f (r= %p)\n",p->orientation ? "left" : "right", Ymap2->id,Ymap2->site[0][N2+1], r);
		    fflush(stdout);
		  }
		  continue;
		}

		/* clean up previous overlap alignment pointers in Xmap */
		if(Xmap == Ymap2){
		  if(Xmap2->leftoverlap == r){
		    if(VERB/* HERE >=2 */){
		      printf("\t Xmap2->id= %lld,len=%0.3f: r=leftoverlap= %p -> 0\n",Xmap2->id,Xmap2->site[0][M2+1], Xmap2->leftoverlap);
		      fflush(stdout);
		    }
		    Xmap2->leftoverlap = NULL;
		  } else {
		    if(VERB/* HERE >=2 */){
		      printf("\t Xmap2->id= %lld,len=%0.3f: r=rightoverlap= %p -> 0\n",Xmap2->id,Xmap2->site[0][M2+1], Xmap2->rightoverlap);
		      fflush(stdout);
		    }
		    if(DEBUG) assert(Xmap2->rightoverlap == r);
		    Xmap2->rightoverlap = NULL;
		  }
		} else {
		  if(Ymap2->leftoverlap == r){
		    if(VERB/* HERE >=2 */){
		      printf("\t Ymap2->id= %lld,len=%0.3f: r=leftoverlap= %p -> 0\n",Ymap2->id,Ymap2->site[0][N2+1], Ymap2->leftoverlap);
		      fflush(stdout);
		    }
		    Ymap2->leftoverlap = NULL;
		  } else {
		    if(VERB/* HERE >=2 */){
		      printf("\t Ymap2->id= %lld,len=%0.3f: r=rightoverlap= %p -> 0\n",Ymap2->id,Ymap2->site[0][N2+1], Ymap2->rightoverlap);
		      fflush(stdout);
		    }
		    if(DEBUG) assert(Ymap2->rightoverlap == r);
		    Ymap2->rightoverlap = NULL;
		  }
		}
	      }

	      if(q){		/* clean up previous overlap alignment pointers in Xmap */
		Cmap *Ymap2 = YYmap[q->mapid1];
		Cmap *Xmap2 = XXmap[q->mapid2];
		int N2 = Ymap2->numsite[0];
		int M2 = Xmap2->numsite[0];
		if(DEBUG) assert(Ymap == Ymap2 || Ymap == Xmap2);
		if(Ymap == Ymap2){
		  if(Xmap2->leftoverlap == q){
		    if(VERB/* HERE >=2 */){
		      printf("\t Xmap2->id= %lld,len=%0.3f: q=leftoverlap= %p -> 0\n", Xmap2->id,Xmap2->site[0][M2+1], Xmap2->leftoverlap);
		      fflush(stdout);
		    }
		    Xmap2->leftoverlap = NULL;
		  } else {
		    if(VERB/* HERE >=2 */){
		      printf("\t Xmap2->id= %lld,len=%0.3f: q=rightoverlap= %p -> 0\n",Xmap2->id,Xmap2->site[0][M2+1], Xmap2->rightoverlap);
		      fflush(stdout);
		    }
		    if(DEBUG) assert(Xmap2->rightoverlap == q);
		    Xmap2->rightoverlap = NULL;
		  }
		} else {
		  if(Ymap2->leftoverlap == q){
		    if(VERB/* HERE >=2 */){
		      printf("\t Ymap2->id= %lld,len=%0.3f: q=leftoverlap= %p -> 0\n",Ymap2->id,Ymap2->site[0][N2+1], Ymap2->leftoverlap);
		      fflush(stdout);
		    }
		    Ymap2->leftoverlap = NULL;
		  } else {
		    if(VERB/* HERE >=2 */){
		      printf("\t Ymap2->id= %lld,len=%0.3f: q=rightoverlap= %p -> 0\n",Ymap2->id,Ymap2->site[0][N2+1], Ymap2->rightoverlap);
		      fflush(stdout);
		    }
		    if(DEBUG) assert(Ymap2->rightoverlap == q);
		    Ymap2->rightoverlap = NULL;
		  }
		}
	      }

	      Ymap->leftoverlap = p;
	      if(p->orientation)
		Xmap->leftoverlap = p;
	      else
		Xmap->rightoverlap = p;
	      if(VERB/* HERE >=2*/){
		printf("Postponing merging Ymap->id=%lld with Xmap->id=%lld : left end of Ymap overlaps %s end of Xmap (q= %p -> %p, r= %p -> %p),Xmap->contig=%p,Ymap->contig=%p\n",
		       Ymap->id,Xmap->id,p->orientation ? "left":"right",q,p,r,p,Xmap->contig,Ymap->contig);
		fflush(stdout);
	      }
	    }
	    continue;
	  }

	  if(!Ymap->contig){
	    if(DEBUG) assert(!(Ymap->centercloned & 2));
	    Ymap->contig = new Ccontig(Ymap);/* convert Ymap into a trivial Hmap with no Het events */
	  }
	  if(!Xmap->contig){
	    if(DEBUG) assert(!(Xmap->centercloned & 2));
	    Xmap->contig = new Ccontig(Xmap);/* convert Xmap into a trivial Hmap with no Het events */
	  }

	  /* Always merge the smaller map into the larger map : If the smaller map is already an Hmap, both Alleles need to be merged into the larger map's 2 Alleles */
	  if(Y[N+1] > X[M+1]){/* merge Xmap into Ymap */
	    if(VERB/* HERE >=2*/){
	      printf("Merging Ymap->id=%lld with Xmap->id=%lld into Hmap\n",Ymap->id,Xmap->id);
	      fflush(stdout);
	    }
	    if(merge_contigs(Ymap->contig,Xmap->contig,p, outlierExtend)){/* failed to merge */
	      if(VERB){
		printf("Failed to merge Ymap->id=%lld with Xmap->id=%lld into Hmap\n",Ymap->id,Xmap->id);
		fflush(stdout);
	      }
	      if(!(Xmap->centercloned & 2)){
		delete Xmap->contig;
		Xmap->contig = NULL;
	      }
	      if(!(Ymap->centercloned & 2)){
		delete Ymap->contig;
		Ymap->contig = NULL;
	      }
	    } else {
	      Ymap->contig->id = Ymap->id;// WAS min(Ymap->id, Xmap->id);
	      if(VERB){
		int N = Ymap->contig->numsite[0];
		printf("Merged Ymap->mapid=%d(id=%lld,paired=%d) with Xmap->mapid=%d(id=%lld,paired=%d->%d) into Hmap(Ymap->contig->id=%lld,sites=%d,len=%0.3f)\n",
		       Ymap->mapid,Ymap->id,Ymap->paired,Xmap->mapid,Xmap->id,Xmap->paired,1+p->mapid1, Ymap->contig->id,N,Ymap->contig->site[0][N+1]);
		fflush(stdout);
	      }
	      Xmap->paired = 1 + p->mapid1;/* Mark as used up */
	      Ymap->centercloned |= 2;/* Mark as modified/Hmap */
	      // WAS	      Ymap->id = min(Ymap->id, Xmap->id);
	    }
	  } else {/* merge Ymap into Xmap */
	    if(VERB/* HERE >=2 */){
	      printf("Merging Xmap->id=%lld with Ymap->id=%lld into Hmap\n",Xmap->id,Ymap->id);
	      fflush(stdout);
	    }

	    Calign *q = new Calign;
	    q->mapid1 = p->mapid2;
	    q->mapid2 = p->mapid1;
	    int U = p->numpairs;
	    q->numpairs = U;
	    q->orientation = p->orientation;

	    if(q->orientation){/* need to create alignment of Xmap with flipped Ymap */
	      if(VERB/* HERE >=2 */){
		printf("Merging Xmap->id=%lld with flipped Ymap->id=%lld into Hmap\n",Xmap->id,Ymap->id);
		if(VERB>=3 && Xmap->id==1790 && Xmap->contig->contig[0].numsite[0]==177)
		  printf("\t Xmap->contig:contig[0].numsite[0]=%d, sitemap[0][0][70]=%d,sitemap[0][0][177]=%d\n",
			 Xmap->contig->contig[0].numsite[0],Xmap->contig->sitemap[0][0][70],Xmap->contig->sitemap[0][0][177]);
		fflush(stdout);
	      }

	      q->allocated_size = U;
	      q->sites1 = new int[U];
	      q->sites2 = new int[U];
	      q->iscore = new FLOAT[U+1];
	      q->outscore = new FLOAT[U+1];

	      for(int t = 0; t < U; t++){
		q->sites1[U-1-t] = M + 1 - p->sites2[t];
		q->sites2[U-1-t] = N + 1 - p->sites1[t];
	      }
	      for(int t = 0; t <= U; t++){
		q->iscore[U-t] = p->iscore[t];
		q->outscore[U-t] = p->outscore[t];
	      }
	      if(merge_contigs(Xmap->contig,Ymap->contig,q,outlierExtend)){
		if(VERB){
		  printf("Failed to merge Xmap->id=%lld with flipped Ymap->id=%lld (or=%d) into Hmap\n",Xmap->id,Ymap->id, p->orientation);
		  fflush(stdout);
		}
		if(!(Xmap->centercloned & 2)){
		  delete Xmap->contig;
		  Xmap->contig = NULL;
		}
		if(!(Ymap->centercloned & 2)){
		  delete Ymap->contig;
		  Ymap->contig = NULL;
		}
	      } else {
		Xmap->contig->id = Xmap->id; // WAS min(Xmap->id, Ymap->id);
		if(VERB){
		  int N = Xmap->contig->numsite[0];
		  printf("Merged Xmap->mapid=%d(id=%lld,paired=%d) with flipped Ymap->mapid=%d(id=%lld,paired=%d->%d) (or=%d) as Hmap(Xmap->contig:id=%lld,sites=%d,len=%0.3f)\n",
			 Xmap->mapid,Xmap->id,Xmap->paired,Ymap->mapid,Ymap->id,Ymap->paired,1+p->mapid2, p->orientation, Xmap->contig->id,N,Xmap->contig->site[0][M+1]);
		  fflush(stdout);
		}
		Ymap->paired = 1 + p->mapid2;
		Xmap->centercloned |= 2;
		// WAS		Xmap->id = min(Xmap->id, Ymap->id);
	      }
	    } else {
	      if(VERB/* HERE >=2 */){
		printf("Merging Xmap->id=%lld with unflipped Ymap->id=%lld into Hmap\n",Xmap->id,Ymap->id);
		fflush(stdout);
	      }
	      /* first create a copy of p == alignment[i] in which sites1[] and sites2[] are swapped */
	      q->sites1 = p->sites2;
	      q->sites2 = p->sites1;
	      q->iscore = p->iscore;
	      q->outscore = p->outscore;
	      
	      if(merge_contigs(Xmap->contig,Ymap->contig,q,outlierExtend)){
		if(VERB){
		  printf("Failed to merge Xmap->id=%lld with unflipped Ymap->id=%lld (or=%d) as Hmap\n",Xmap->id,Ymap->id, p->orientation);
		  fflush(stdout);
		}
		if(!(Xmap->centercloned & 2)){
		  delete Xmap->contig;
		  Xmap->contig = NULL;
		}
		if(!(Ymap->centercloned & 2)){
		  delete Ymap->contig;
		  Ymap->contig = NULL;
		}
	      } else {
		Xmap->contig->id = Xmap->id; // WAS min(Xmap->id, Ymap->id);
		if(VERB){
		  int N = Xmap->contig->numsite[0];
		  printf("Merged Xmap->mapid=%d(id=%lld,paired=%d) with unflipped Ymap->mapid=%d(id=%lld,paired=%d->%d) (or=%d) as Hmap(Xmap->contig:id=%lld,sites=%d,len=%0.3f)\n",
			 Xmap->mapid,Xmap->id,Xmap->paired,Ymap->mapid,Ymap->id, Ymap->paired,1+p->mapid2,p->orientation, Xmap->contig->id, N, Xmap->contig->site[0][N+1]);
		  fflush(stdout);
		}
		Ymap->paired = 1 + p->mapid2;
		Xmap->centercloned |= 2;
		// WAS		Xmap->id = min(Xmap->id, Ymap->id);
	      }
	      q->sites1 = q->sites2 = NULL;// should not be needed since q->allocated_size = 0;
	    }
	    delete q;
	  }
	  continue;
	}// if(PairMergeHmap > 0 && LeftLen >= TruncateHmap && RightLen >= TruncateHmap ...)
      } // if(!merge && Truncate > 0.0 && overlap >= TruncateHmap && K >= TruncateNHmap)

      if(!merge && PairMergeHmap > 0 && overlap >= TruncateHmap && (maxOutlier > PairMergeMaxOutlier || maxOutlierM > PairMergeMaxOutlierM)){
	/* check if we can merge shorter map that is fully overlapped by larger map into an Hmap */
        /* next check if there is only a single large endoutlier and we can find another map that has a matching single large endoutlier with the same large map */
	if(min(leftEnd,rightEnd) <= max(overlap * PairMergeMaxEnd, PairMergeMaxEndKB)){
	  // HERE HERE

	  //	  printf("WARNING: Generation of merged Hmap from endoutlier alignment not yet implemented\n");
	  //	  fflush(stdout);
	}
      }

      if(VERB>=3){
        printf("merge=%d: SplitSegDup= %0.3f, Overlapped=%d\n",merge,SplitSegDup,Overlapped);
	fflush(stdout);
      }

      if(merge && (SplitSegDup <= 0.0 || Overlapped)){ /* create a new map based on the merger of X and Y */
	if(PairMergeHmap > 0){
	  printf("WARNING: with PairMergeHmap= %d contigs should have been merged previously (TrimmedNoMerge= %0.3e %0.3f): skipping\n", PairMergeHmap, TrimmedNoMerge, TrimmedNoMergeMaxEnd);
	  fflush(stdout);

	  continue;
	}

	/*	if(SplitSegDup > 0.0 && !((XL | XR | YL | YR) & END_NOEXT)){
	  printf("WARNING: with SplitSegDup= %0.3f contigs should have been merged previously (TrimmedNoMerge= %0.3e %0.3f, XL=%lu,XR=%lu,YL=%lu,YR=%lu): skipping\n", 
		 SplitSegDup, TrimmedNoMerge, TrimmedNoMergeMaxEnd,XL,XR,YL,YR);
	  fflush(stdout);
	  
	  assert(!(merge && SplitSegDup > 0.0));
	  }*/

	Xmap->paired = 1 + nummaps;
	Ymap->paired = 1 + nummaps;

	maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);

	/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	if(DEBUG) assert(PairSplit <= 0.0);
	if(!Cfirst || CfirstPairs==2)
	  XXmap = YYmap = Gmap;
	else {
	  YYmap = &Gmap[0];
	  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	}

	Cmap *pmap = Gmap[nummaps];
	if(DEBUG && !PairMergeIDrenumber) assert(pmap != NULL);

	if(!pmap)
	  pmap = Gmap[nummaps] = new Cmap;// NOTE : this memory will be leaked.

	pmap->mapid = nummaps++;
	pmap->fileid = -1;
	pmap->paired = 0;
	pmap->centercloned = 1;
	pmap->origmap = 0;
	if(Ymap->id < Xmap->id){
	  pmap->id = Ymap->id;
	  pmap->origid = Ymap->origid;
	} else {
	  pmap->id = Xmap->id;
	  pmap->origid = Xmap->origid;
	}

	for(int c = 0; c < colors; c++){
	  pmap->Nickase[c] = Xmap->Nickase[c];
	  pmap->numsite[c] = 0;
	}
	pmap->startloc = pmap->endloc = pmap->len = 0.0;

	int I = p->sites1[K/2], J = p->sites2[K/2];

	if(VERB>=2 && Ymap->id==1776 && Xmap->id==11367){
	  printf("Yid=%lld,Xid=%lld,or=%d:centerX=%0.6f,centerY=%0.6f,distanceYtoX=%0.6f,I=%d,J=%d: merging (nummaps=%d)\n",
		 Ymap->id,Xmap->id,p->orientation,centerX,centerY,distanceYtoX,I,J,nummaps);
	  fflush(stdout);
	}

	/* If ends match, prefer to interpret smaller map overlapped by larger map, to maximize chance of just keeping the larger map */
	if(Y[N+1] - X[M+1] >= -PairMergeOverlappedEmargin && (Y[N+1] - X[M+1] >= PairMergeOverlappedSmargin || N >= M + PairMergeOverlappedNmargin) &&
	   distanceYtoX >= centerX - centerY - PairMergeOverlappedMargin && distanceYtoX <= centerY - centerX + PairMergeOverlappedMargin)
	  goto L_X_inside_Y;

	if(X[M+1] - Y[N+1] >= -PairMergeOverlappedEmargin && (X[M+1] - Y[N+1] >= PairMergeOverlappedSmargin || N >= M + PairMergeOverlappedNmargin) &&
	   distanceYtoX <= centerX - centerY + PairMergeOverlappedMargin && distanceYtoX >= centerY - centerX - PairMergeOverlappedMargin)
	  goto L_Y_inside_X;

	if(SplitSegDup > 0.0){
	  if(DEBUG>=1+RELEASE) assert(Overlapped);
	  if(LeftSmallerY && RightSmallerY)
	    goto L_Y_inside_X;
	  if(DEBUG>=1+RELEASE) assert(LeftSmallerX && RightSmallerX);
	  goto L_X_inside_Y;
	}

	if(distanceYtoX > centerX - centerY){/* left end of Y comes first */

	  if(distanceYtoX <= centerY - centerX){ /* right end of Y comes last */
	  L_X_inside_Y:
	    int xstart = 0, xend = 0, ystart = 0, yend = 0;
	    if(PairMergePreferNGS){  /*  check if X has lower SD in some region of the alignment, if so copy X in that aligned region */
	      double *VARdiff = new double[K];
	      for(int k = 0; k < K-1; k++){
		int i = p->sites1[k];
		int j = p->sites2[k];
		int I = p->sites1[k+1];
		int J = p->sites2[k+1];
		double varY = 0, varX = 0.0;
		for(int t = i; t < I; t++){
		  double SDy = Ymap->siteSD[0][t];
		  varY += SDy*SDy;
		}
		for(int t = j; t < J; t++){
		  double SDx = p->orientation ? Xmap->siteSD[0][M-t] : Xmap->siteSD[0][t];
		  varX += SDx*SDx;
		}
		VARdiff[k] = varY - varX;
	      }
	      /* Now locate interval of X (xstart .. xend) with greatest reduction in variance over Y (largest sum of Vardiff[k = 0 .. K-2])*/
	      double sum = 0.0;/* cumulative sum of VARdiff[0..k] */
	      double sumLWM = 0.0;/*  Low Water Mark of cumulative sum of VARdiff[0..kLWM-1] */
	      int kLWM = 0;/* k value of Low Water Mark */
	      double maxdiff = 0.0;/* largest value of sum-sumLWM */
	      for(int k = 0; k < K-1; k++){
		sum += VARdiff[k];
		if(sum < sumLWM){
		  sumLWM = sum;
		  kLWM = k+1;
		}
		if(sum - sumLWM > maxdiff){
		  maxdiff = sum - sumLWM;
		  ystart = p->sites1[kLWM];
		  xstart = p->sites2[kLWM];
		  yend = p->sites1[k+1];
		  xend = p->sites2[k+1];
		}
	      }
	      if(VERB/* HERE >=2 */ && xstart < xend){
		printf("    Retaining part of overlapped X from X[%d]=%0.3f to X[%d]=%0.3f (replacing Y[%d]=%0.3f to Y[%d]=%0.3f\n",
		       xstart, Xmap->site[0][p->orientation ? (M+1-xstart) : xstart],
		       xend, Xmap->site[0][p->orientation ? (M+1-xend) : xend],
		       ystart, Ymap->site[0][ystart], yend, Ymap->site[0][yend]);
		fflush(stdout);
	      }
	    }

	    if(xstart < xend){/* copy Y from left end to ystart, then X from xstart to xend, then Y from yend to right end */
	      /* compute number of sites and allocate memory */
	      int NumSites = N - (yend - ystart) + (xend - xstart);
	      for(int c = 0; c < colors; c++){
		if(pmap->site[c]){
		  if(!pmap->blockmem) delete [] pmap->site[c];
		  if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		  if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		  if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
		  if(CmapChimQuality >= 1){
		    if(pmap->ChimQuality[c]) delete [] pmap->ChimQuality[c];
		    if(CmapChimQuality >= 2){
		      if(pmap->ChimNorm[c]) delete [] pmap->ChimNorm[c];
		      if(pmap->SegDupL[c]) delete [] pmap->SegDupL[c];
		      if(pmap->SegDupR[c]) delete [] pmap->SegDupR[c];
		      if(pmap->FragileEndL[c]) delete [] pmap->FragileEndL[c];
		      if(pmap->FragileEndR[c]) delete [] pmap->FragileEndR[c];
		      if(pmap->OutlierFrac[c]) delete [] pmap->OutlierFrac[c];
		      if(CmapChimQuality >= 3){
			if(pmap->FragSd[c]) delete [] pmap->FragSd[c];
			if(pmap->ExpSd[c]) delete [] pmap->ExpSd[c];
			if(pmap->FragCov[c]) delete [] pmap->FragCov[c];
			if(pmap->FragChiSq[c]) delete [] pmap->FragChiSq[c];
		      }
		    }
		  }
		}
		if(pmap->Mask[c]){
		  delete [] pmap->Mask[c];
		  pmap->Mask[c] = NULL;
		}
		if(pmap->SNRcnt[c]){
		  if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		  if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		  if(pmap->SNRdist[c]){
		    for(int i = 0; i <= pmap->numsite[0]+1; i++)
		      if(pmap->SNRdist[c][i]) {
			if(VERB>=3){
			  printf("Deleting Gmap[%d](%p)->SNRdist[%d][%d]=%p(SNRcnt=%d)\n",pmap->mapid,pmap,c,i,pmap->SNRdist[c][i],pmap->SNRcnt[c][i]);
			  fflush(stdout);
			}
			delete [] pmap->SNRdist[c][i];
		      }
		    delete [] pmap->SNRdist[c];
		  }
		  delete [] pmap->SNRcnt[c];
		}

		pmap->numsite[c] = NumSites;
		pmap->site[c] = new FLOAT[NumSites+2];
		pmap->siteSD[c] = new double[NumSites+2];
		pmap->sitecov[c] = new float[NumSites+2];
		pmap->sitecnt[c] = new float[NumSites+2];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c] = new float[NumSites+2];
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c] = new float[NumSites+2];
		    pmap->SegDupL[c] = new float[NumSites+2];
		    pmap->SegDupR[c] = new float[NumSites+2];
		    pmap->FragileEndL[c] = new float[NumSites+2];
		    pmap->FragileEndR[c] = new float[NumSites+2];
		    pmap->OutlierFrac[c] = new float[NumSites+2];
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c] = new float[NumSites+2];
		      pmap->ExpSd[c] = new float[NumSites+2];
		      pmap->FragCov[c] = new float[NumSites+2];
		      pmap->FragChiSq[c] = new float[NumSites+2];
		    }
		  }
		}
		if(Ymap->Mask[c] || Xmap->Mask[c]){
		  pmap->Mask[c] = new size_t[NumSites+2];
		  memset(pmap->Mask[c], 0, (NumSites+2)*sizeof(size_t));
		}
		pmap->SNRcnt[c] = 0;
		if(Ymap->SNRcnt[c] || Xmap->SNRcnt[c]){
		  pmap->SNRcnt[c] = new int[NumSites+2];
		  pmap->SNRgmean[c] = new double[NumSites+2];
		  pmap->lnSNRsd[c] = new double[NumSites+2];
		  pmap->SNRdist[c] = new double*[NumSites+2];
		}

		/* left end */
		pmap->site[c][0] = 0.0;
		pmap->siteSD[c][0] = 0.0;
		pmap->sitecov[c][0] = 1;
		pmap->sitecnt[c][0] = 1;
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c][0] = END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c][0] = END_CHIMQUAL;
		    pmap->SegDupL[c][0] = END_CHIMQUAL;
		    pmap->SegDupR[c][0] = END_CHIMQUAL;
		    pmap->FragileEndL[c][0] = 0.0;
		    pmap->FragileEndR[c][0] = 0.0;
		    pmap->OutlierFrac[c][0] = 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c][0] = 0.0;
		      pmap->ExpSd[c][0] = 0.0;
		      pmap->FragCov[c][0] = 0.0;
		      pmap->FragChiSq[c][0] = 1.0;
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][0] = 0;
		  pmap->SNRgmean[c][0] = 0.0;
		  pmap->lnSNRsd[c][0] = 0.0;
		  pmap->SNRdist[c][0] = NULL;
		}

		/* right end */
		pmap->siteSD[c][NumSites+1] = 0.0;
		pmap->sitecov[c][NumSites+1] = 1;
		pmap->sitecnt[c][NumSites+1] = 1;
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c][NumSites+1] = END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c][NumSites+1] = END_CHIMQUAL;
		    pmap->SegDupL[c][NumSites+1] = END_CHIMQUAL;
		    pmap->SegDupR[c][NumSites+1] = END_CHIMQUAL;
		    pmap->FragileEndL[c][NumSites+1] = 0.0;
		    pmap->FragileEndR[c][NumSites+1] = 0.0;
		    pmap->OutlierFrac[c][NumSites+1] = 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c][NumSites+1] = 0.0;
		      pmap->ExpSd[c][NumSites+1] = 0.0;
		      pmap->FragCov[c][NumSites+1] = 0.0;
		      pmap->FragChiSq[c][NumSites+1] = 1.0;
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][NumSites+1] = 0;
		  pmap->SNRgmean[c][NumSites+1] = 0.0;
		  pmap->lnSNRsd[c][NumSites+1] = 0.0;
		  pmap->SNRdist[c][NumSites+1] = NULL;
		}

		if(pmap->Mask[c] && Ymap->Mask[c]){
		  pmap->Mask[c][1] |= Ymap->Mask[c][1];
		  pmap->Mask[c][NumSites+1] |= Ymap->Mask[c][N+1];
		}
	      }
	      pmap->blockmem = 0;
	      FLOAT *Z = pmap->site[0];

	      /* copy left end of Ymap from site 0 to site ystart */
	      register int i = 0;
	      for(; i <= ystart; i++){
		Z[i] = Y[i];
		if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		pmap->siteSD[0][i] = Ymap->siteSD[0][i];
		pmap->sitecov[0][i] = Ymap->sitecov[0][i];
		pmap->sitecnt[0][i] = Ymap->sitecnt[0][i];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[0][i] = Ymap->ChimQuality[0] ? Ymap->ChimQuality[0][i] : END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[0][i] = Ymap->ChimNorm[0] ? Ymap->ChimNorm[0][i] : END_CHIMQUAL;
		    pmap->SegDupL[0][i] = Ymap->SegDupL[0] ? Ymap->SegDupL[0][i] : END_CHIMQUAL;
		    pmap->SegDupR[0][i] = Ymap->SegDupR[0] ? Ymap->SegDupR[0][i] : END_CHIMQUAL;
		    pmap->FragileEndL[0][i] = Ymap->FragileEndL[0] ? Ymap->FragileEndL[0][i] : 0.0;
		    pmap->FragileEndR[0][i] = Ymap->FragileEndR[0] ? Ymap->FragileEndR[0][i] : 0.0;
		    pmap->OutlierFrac[0][i] = Ymap->OutlierFrac[0] ? Ymap->OutlierFrac[0][i] : 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[0][i] = Ymap->FragSd[0] ? Ymap->FragSd[0][i] : 0.0;
		      pmap->ExpSd[0][i] = Ymap->ExpSd[0] ? Ymap->ExpSd[0][i] : 0.0;
		      pmap->FragCov[0][i] = Ymap->FragCov[0] ? Ymap->FragCov[0][i] : 0.0;
		      pmap->FragChiSq[0][i] = Ymap->FragChiSq[0] ? Ymap->FragChiSq[0][i] : 1.0;
		    }
		  }
		}
		if(Ymap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = Ymap->SNRcnt[0][i];
		  pmap->SNRgmean[0][i] = Ymap->SNRgmean[0][i];
		  pmap->lnSNRsd[0][i] = Ymap->lnSNRsd[0][i];
		  pmap->SNRdist[0][i] = 0;
		  if(pmap->SNRcnt[0][i]){
		    pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		    if(VERB>=3){
		      printf("Allocating Gmap[%d](%p)->SNRdist[%d][%d] = %p (SNRcnt=%d)\n",pmap->mapid,pmap,0,i,pmap->SNRdist[0][i],pmap->SNRcnt[0][i]);
		      fflush(stdout);
		    }

		    for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
		      pmap->SNRdist[0][i][t] = Ymap->SNRdist[0][i][t];
		  }
		} else if(pmap->SNRcnt[0]){/* default to empty SNR distribution */
		  pmap->SNRcnt[0][i] = 0;
		  pmap->SNRgmean[0][i] = 0.0;
		  pmap->lnSNRsd[0][i] = 0.0;
		  pmap->SNRdist[0][i] = NULL;
		}
		if(VERB>=3 && Xmap->id==13905){
		  printf("Copying site Y[i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f,siteSD=%0.6f\n",
			 i,N, Y[i],Ymap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i], pmap->siteSD[0][i]);
		  fflush(stdout);
		}
	      }
	      if(DEBUG) assert(i-1 == ystart);

	      if(p->orientation){/* copy Xmap from M+1-xstart down to M+1-xend */
		int kstart = M+1-xstart;
		int kend = M+1-xend;
		double Xoffset = Z[i-1] + X[kstart];
		for(int j = kstart; --j >= kend; i++){
		  Z[i] = Xoffset - X[j];
		  if(DEBUG) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[0][i-1] = Xmap->siteSD[0][j];/* see definition of siteSD[] */
		  pmap->sitecov[0][i] = Xmap->sitecov[0][j];
		  pmap->sitecnt[0][i] = Xmap->sitecnt[0][j];
		  if(CmapChimQuality >= 1){
		    pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][j] : END_CHIMQUAL;
		    if(CmapChimQuality >= 2){
		      pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][j] : END_CHIMQUAL;
		      pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][j] : END_CHIMQUAL;
		      pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][j] : END_CHIMQUAL;
		      pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][j] : 0.0;
		      pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][j] : 0.0;
		      pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][j] : 0.0;
		      if(CmapChimQuality >= 3){
			pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][j] : 0.0;
			pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][j] : 0.0;
			pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][j] : 0.0;
			pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][j] : 1.0;
		      }
		    }
		  }
		  if(Xmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][j];
		    pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][j];
		    pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][j];
		    pmap->SNRdist[0][i] = 0;
		    if(pmap->SNRcnt[0][i]){
		      pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		      if(VERB>=3){
			printf("Allocating Gmap[%d](%p)->SNRdist[%d][%d] = %p (SNRcnt=%d)\n",pmap->mapid,pmap,0,i,pmap->SNRdist[0][i],pmap->SNRcnt[0][i]);
			fflush(stdout);
		      }

		      for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
			pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][j][t];
		    }
		  } else if(pmap->SNRcnt[0]){/* use default empty SNR distribution */
		    pmap->SNRcnt[0][i] = 0;
		    pmap->SNRgmean[0][i] = 0.0;
		    pmap->lnSNRsd[0][i] = 0.0;
		    pmap->SNRdist[0][i] = NULL;
		  }
		  if(VERB>=3 && Xmap->id==13905){
		    printf("Copying site X[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f,siteSD=%0.6f\n",
			   j, M, X[j],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i],pmap->siteSD[0][i-1]);
		    fflush(stdout);
		  }
		}
		if(DEBUG && !(i-1 == ystart + (xend-xstart))){
		  printf("N=%d,M=%d,N=%d,M=%d,I=%d,J=%d,i=%d,NumSites=%d,orientation=%d\n",
			 N,M,N,M,I,J,i,NumSites,p->orientation);
		  fflush(stdout);
		  assert(i-1 == ystart + (xend-xstart));
		}

	      } else {/* copy Xmap from xstart to xend */

		int kstart = xstart;
		int kend = xend;
		double Xoffset = Z[i-1] - X[kstart];
		for(int j = kstart; ++j <= kend; i++){
		  Z[i] = Xoffset + X[j];
		  if(DEBUG) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[0][i-1] = Xmap->siteSD[0][j-1];/* see definition of siteSD[] */
		  pmap->sitecov[0][i] = Xmap->sitecov[0][j];
		  pmap->sitecnt[0][i] = Xmap->sitecnt[0][j];
		  if(CmapChimQuality >= 1){
		    pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][j] : END_CHIMQUAL;
		    // WAS		    pmap->ChimQuality[0][i] = Xmap->sitecnt[0] ? Xmap->sitecnt[0][j] : -1.0;
		    if(CmapChimQuality >= 2){
		      pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][j] : END_CHIMQUAL;
		      pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][j] : END_CHIMQUAL;
		      pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][j] : END_CHIMQUAL;
		      pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][j] : 0.0;
		      pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][j] : 0.0;
		      pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][j] : 0.0;
		      if(CmapChimQuality >= 3){
			pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][j] : 0.0;
			pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][j] : 0.0;
			pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][j] : 0.0;
			pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][j] : 1.0;
		      }
		    }
		  }
		  if(Xmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][j];
		    pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][j];
		    pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][j];
		    pmap->SNRdist[0][i] = NULL;
		    if(pmap->SNRcnt[0][i]){
		      pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		      if(VERB>=3){
			printf("Allocating Gmap[%d](%p)->SNRdist[%d][%d] = %p (SNRcnt=%d)\n",pmap->mapid,pmap,0,i,pmap->SNRdist[0][i],pmap->SNRcnt[0][i]);
			fflush(stdout);
		      }
		      for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
			pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][j][t];
		    }
		  } else if(pmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = 0;
		    pmap->SNRgmean[0][i] = 0.0;
		    pmap->lnSNRsd[0][i] = 0.0;
		    pmap->SNRdist[0][i] = NULL;
		  }
		  if(VERB>=3 && Xmap->id==13905){
		    printf("Copying site X[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f,siteSD=%0.6f\n",
			   j, M, X[j],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i],pmap->siteSD[0][i-1]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == ystart + (xend-xstart));
	      }

	      /* copy Ymap from yend to right end site N+1 */
	      int kstart = yend;
	      int kend = N+1;
	      double Yoffset = Z[i-1] - Y[kstart];
	      for(int j = kstart; ++j <= kend; i++){
		Z[i] = Yoffset + Y[j];
		if(DEBUG && i > 0 && !(Z[i] >= Z[i-1])){
		  printf("yend=%d,N=%d,kstart=%d\n",yend,N,kstart);
		  printf("j=%d,i=%d:Z[i]=%0.6f,Z[i-1]=%0.6f:Y[kstart]=%0.6f,Y[j]=%0.6f\n",
			 j,i,Z[i],Z[i-1],Y[kstart],Y[j]);
		  fflush(stdout);
		  assert(Z[i] >= Z[i-1]);
		}
		pmap->siteSD[0][i-1] = Ymap->siteSD[0][j-1];/* see definition of siteSD[] */
		pmap->sitecov[0][i] = Ymap->sitecov[0][j];
		pmap->sitecnt[0][i] = Ymap->sitecnt[0][j];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[0][i] = Ymap->ChimQuality[0] ? Ymap->ChimQuality[0][j] : END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[0][i] = Ymap->ChimNorm[0] ? Ymap->ChimNorm[0][j] : END_CHIMQUAL;
		    pmap->SegDupL[0][i] = Ymap->SegDupL[0] ? Ymap->SegDupL[0][j] : END_CHIMQUAL;
		    pmap->SegDupR[0][i] = Ymap->SegDupR[0] ? Ymap->SegDupR[0][j] : END_CHIMQUAL;
		    pmap->FragileEndL[0][i] = Ymap->FragileEndL[0] ? Ymap->FragileEndL[0][j] : 0.0;
		    pmap->FragileEndR[0][i] = Ymap->FragileEndR[0] ? Ymap->FragileEndR[0][j] : 0.0;
		    pmap->OutlierFrac[0][i] = Ymap->OutlierFrac[0] ? Ymap->OutlierFrac[0][j] : 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[0][i] = Ymap->FragSd[0] ? Ymap->FragSd[0][j] : 0.0;
		      pmap->ExpSd[0][i] = Ymap->ExpSd[0] ? Ymap->ExpSd[0][j] : 0.0;
		      pmap->FragCov[0][i] = Ymap->FragCov[0] ? Ymap->FragCov[0][j] : 0.0;
		      pmap->FragChiSq[0][i] = Ymap->FragChiSq[0] ? Ymap->FragChiSq[0][j] : 1.0;
		    }		    
		  }
		}
		if(Ymap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = Ymap->SNRcnt[0][j];
		  pmap->SNRgmean[0][i] = Ymap->SNRgmean[0][j];
		  pmap->lnSNRsd[0][i] = Ymap->lnSNRsd[0][j];
		  pmap->SNRdist[0][i] = NULL;
		  if(pmap->SNRcnt[0][i]){
		    pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		    if(VERB>=3){
		      printf("Allocating Gmap[%d](%p)->SNRdist[%d][%d] = %p (SNRcnt=%d)\n",pmap->mapid,pmap,0,i,pmap->SNRdist[0][i],pmap->SNRcnt[0][i]);
		      fflush(stdout);
		    }

		    for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
		      pmap->SNRdist[0][i][t] = Ymap->SNRdist[0][j][t];
		  }
		} else if(pmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = 0;
		  pmap->SNRgmean[0][i] = 0.0;
		  pmap->lnSNRsd[0][i] = 0.0;
		  pmap->SNRdist[0][i] = NULL;		
		}
		if(VERB>=3 && Xmap->id==13905){
		  printf("Copying site Y[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f,siteSD=%0.6f\n",
			 j,N, Y[j],Ymap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i],pmap->siteSD[0][i-1]);
		  fflush(stdout);
		}
	      }
	      if(DEBUG) assert(i-1 == NumSites+1);
	    } else { /* just copy Y */

	      if(PairMergeRepeat){ // new code : instead of copying Y, just restore original Y
		if(VERB){
		  printf("%d:Created merged map (id=%lld,len=%0.3f) based on alignment between Y mapid=%d(id=%lld,len=%0.3f) and X mapid=%d(id=%lld,len=%0.3f):mapid=%d,id=%lld,len=%0.3f,sites=%d (Just keeping Ymap=%p)\n",
			 pairmergeIter,Ymap->id, Ymap->site[0][N+1], Ymap->mapid,Ymap->id, Ymap->site[0][N+1], 
			 Xmap->mapid,Xmap->id, Xmap->site[0][Xmap->numsite[0]+1],Ymap->mapid,Ymap->id,Ymap->len,Ymap->numsite[0],Ymap);
		  //		printf("%d:Restoring Ymap(id=%lld,mapid=%d) discarding Xmap(id=%lld,mapid=%d)\n",pairmergeIter,Ymap->id,Ymap->mapid,Xmap->id,Xmap->mapid);
		  fflush(stdout);
		}
		Ymap->paired = 0;
		if(DEBUG) assert(Ymap->mapid == p->mapid1);
		Xmap->paired = 1 + Ymap->mapid;
		nummaps--; 

		continue;/* next Gmap[i] */
	      }

	      /* compute number of sites and allocate memory : NOTE Cmap->orig* fields are not changed and become invalid/undefined, flagged with centercloned=1 */
	      register int NumSites = N;
	      for(register int c = 0; c < colors; c++){
		if(pmap->site[c]){
		  if(!pmap->blockmem) delete [] pmap->site[c];
		  if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		  if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		  if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
		  if(CmapChimQuality >= 1){
		    if(pmap->ChimQuality[c]) delete [] pmap->ChimQuality[c];
		    if(CmapChimQuality >= 2){
		      if(pmap->ChimNorm[c]) delete [] pmap->ChimNorm[c];
		      if(pmap->SegDupL[c]) delete [] pmap->SegDupL[c];
		      if(pmap->SegDupR[c]) delete [] pmap->SegDupR[c];
		      if(pmap->FragileEndL[c]) delete [] pmap->FragileEndL[c];
		      if(pmap->FragileEndR[c]) delete [] pmap->FragileEndR[c];
		      if(pmap->OutlierFrac[c]) delete [] pmap->OutlierFrac[c];
		      if(CmapChimQuality >= 3){
			if(pmap->FragSd[c]) delete [] pmap->FragSd[c];
			if(pmap->ExpSd[c]) delete [] pmap->ExpSd[c];
			if(pmap->FragCov[c]) delete [] pmap->FragCov[c];
			if(pmap->FragChiSq[c]) delete [] pmap->FragChiSq[c];
		      }
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		  if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		  if(pmap->SNRdist[c]){
		    for(int i = 0; i <= pmap->numsite[0]+1; i++)
		      if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		    delete [] pmap->SNRdist[c];
		  }
		  delete [] pmap->SNRcnt[c];
		}
		if(pmap->Mask[c]){
		  delete [] pmap->Mask[c];
		  pmap->Mask[c] = NULL;
		}

		pmap->numsite[c] = NumSites;
		pmap->site[c] = new FLOAT[NumSites+2];
		pmap->siteSD[c] = new double[NumSites+2];
		pmap->sitecov[c] = new float[NumSites+2];
		pmap->sitecnt[c] = new float[NumSites+2];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c] = new float[NumSites+2];
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c] = new float[NumSites+2];
		    pmap->SegDupL[c] = new float[NumSites+2];
		    pmap->SegDupR[c] = new float[NumSites+2];
		    pmap->FragileEndL[c] = new float[NumSites+2];
		    pmap->FragileEndR[c] = new float[NumSites+2];
		    pmap->OutlierFrac[c] = new float[NumSites+2];
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c] = new float[NumSites+2];
		      pmap->ExpSd[c] = new float[NumSites+2];
		      pmap->FragCov[c] = new float[NumSites+2];
		      pmap->FragChiSq[c] = new float[NumSites+2];
		    }
		  }
		}
		pmap->SNRcnt[c] = 0;
		if(Ymap->SNRcnt[c]){
		  pmap->SNRcnt[c] = new int[NumSites+2];
		  pmap->SNRgmean[c] = new double[NumSites+2];
		  pmap->lnSNRsd[c] = new double[NumSites+2];
		  pmap->SNRdist[c] = new double*[NumSites+2];
		}
		if(Ymap->Mask[c]){
		  pmap->Mask[c] = new size_t[NumSites+2];
		  memcpy(pmap->Mask[c],Ymap->Mask[c],(NumSites+2)*sizeof(size_t));
		}

		/* left end */
		pmap->site[c][0] = 0.0;
		pmap->siteSD[c][0] = 0.0;
		pmap->sitecov[c][0] = 1;
		pmap->sitecnt[c][0] = 1;
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c][0] = END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c][0] = END_CHIMQUAL;
		    pmap->SegDupL[c][0] = END_CHIMQUAL;
		    pmap->SegDupR[c][0] = END_CHIMQUAL;
		    pmap->FragileEndL[c][0] = 0.0;
		    pmap->FragileEndR[c][0] = 0.0;
		    pmap->OutlierFrac[c][0] = 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c][0] = 0.0;
		      pmap->ExpSd[c][0] = 0.0;
		      pmap->FragCov[c][0] = 0.0;
		      pmap->FragChiSq[c][0] = 1.0;
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][0] = 0;
		  pmap->SNRgmean[c][0] = 0.0;
		  pmap->lnSNRsd[c][0] = 0.0;
		  pmap->SNRdist[c][0] = NULL;
		}

		/* right end */
		pmap->siteSD[c][NumSites+1] = 0.0;
		pmap->sitecov[c][NumSites+1] = 1;
		pmap->sitecnt[c][NumSites+1] = 1;
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c][NumSites+1] = END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c][NumSites+1] = END_CHIMQUAL;
		    pmap->SegDupL[c][NumSites+1] = END_CHIMQUAL;
		    pmap->SegDupR[c][NumSites+1] = END_CHIMQUAL;
		    pmap->FragileEndL[c][NumSites+1] = 0.0;
		    pmap->FragileEndR[c][NumSites+1] = 0.0;
		    pmap->OutlierFrac[c][NumSites+1] = 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c][NumSites+1] = 0.0;
		      pmap->ExpSd[c][NumSites+1] = 0.0;
		      pmap->FragCov[c][NumSites+1] = 0.0;
		      pmap->FragChiSq[c][NumSites+1] = 1.0;
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][NumSites+1] = 0;
		  pmap->SNRgmean[c][NumSites+1] = 0.0;
		  pmap->lnSNRsd[c][NumSites+1] = 0.0;
		  pmap->SNRdist[c][NumSites+1] = NULL;
		}
	      }
	      pmap->blockmem = 0;

	      register FLOAT *Z = pmap->site[0];

	      /* copy Ymap from site 0 to site N+1 */
	      register int i = 0;
	      for(; i <= N+1; i++){
		Z[i] = Y[i];
		if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		pmap->siteSD[0][i] = Ymap->siteSD[0][i];
		pmap->sitecov[0][i] = Ymap->sitecov[0][i];
		pmap->sitecnt[0][i] = Ymap->sitecnt[0][i];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[0][i] = Ymap->ChimQuality[0] ? Ymap->ChimQuality[0][i] : END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[0][i] = Ymap->ChimNorm[0] ? Ymap->ChimNorm[0][i] : END_CHIMQUAL;
		    pmap->SegDupL[0][i] = Ymap->SegDupL[0] ? Ymap->SegDupL[0][i] : END_CHIMQUAL;
		    pmap->SegDupR[0][i] = Ymap->SegDupR[0] ? Ymap->SegDupR[0][i] : END_CHIMQUAL;
		    pmap->FragileEndL[0][i] = Ymap->FragileEndL[0] ? Ymap->FragileEndL[0][i] : 0.0;
		    pmap->FragileEndR[0][i] = Ymap->FragileEndR[0] ? Ymap->FragileEndR[0][i] : 0.0;
		    pmap->OutlierFrac[0][i] = Ymap->OutlierFrac[0] ? Ymap->OutlierFrac[0][i] : 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[0][i] = Ymap->FragSd[0] ? Ymap->FragSd[0][i] : 0.0;
		      pmap->ExpSd[0][i] = Ymap->ExpSd[0] ? Ymap->ExpSd[0][i] : 0.0;
		      pmap->FragCov[0][i] = Ymap->FragCov[0] ? Ymap->FragCov[0][i] : 0.0;
		      pmap->FragChiSq[0][i] = Ymap->FragChiSq[0] ? Ymap->FragChiSq[0][i] : 1.0;
		    }
		  }
		}
		if(Ymap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = Ymap->SNRcnt[0][i];
		  pmap->SNRgmean[0][i] = Ymap->SNRgmean[0][i];
		  pmap->lnSNRsd[0][i] = Ymap->lnSNRsd[0][i];
		  pmap->SNRdist[0][i] = NULL;
		  if(pmap->SNRcnt[0][i]){
		    pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		    for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
		      pmap->SNRdist[0][i][t] = Ymap->SNRdist[0][i][t];
		  }
		} else if(pmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = 0;
		  pmap->SNRgmean[0][i] = 0.0;
		  pmap->lnSNRsd[0][i] = 0.0;
		  pmap->SNRdist[0][i] = NULL;
		}
		if(VERB>=3){
		  printf("Copying site Y[i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			 i,N, Y[i],Ymap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		  fflush(stdout);
		}
	      }
	      if(DEBUG) assert(i-1 == NumSites+1);	  
            }
	  } else { /* right end of X comes last : copy left end of Y and right end of X */

	    if(PairMergePreferNGS){/* location transition point that retains the most NGS sites (with StdDev==0) */
	      double *VARdiff = new double[K];
	      double VARmin = 0.0, VARsum = 0.0;
	      int kmin = 0;
	      /* NOTE : siteSD[0][i] is the sizing error for interval site[0][i+1] - site[0][i] */
	      for(int k = 0; k < K-1; k++){
		int i = p->sites1[k];
		int j = p->sites2[k];
		int I = p->sites1[k+1];
		int J = p->sites2[k+1];
		double varY = 0, varX = 0.0;
		for(int t = i; t < I; t++){
		  double SDy = Ymap->siteSD[0][t];
		  varY += SDy*SDy;
		  if(VERB>=4)
		    printf("       t=%d:SDy=%0.6f,varY=%0.6e (or=%d, Ymap->siteSD[0][%d]=%0.6f)\n",
			   t,SDy,varY,p->orientation,t,Ymap->siteSD[0][t]);
		}
		for(int t = j; t < J; t++){
		  double SDx = p->orientation ? Xmap->siteSD[0][M-t] : Xmap->siteSD[0][t];
		  varX += SDx*SDx;
		  if(VERB>=4){
		    int index = (p->orientation ? M-t : t);
		    printf("       t=%d:SDx=%0.6f,varX=%0.6e (or=%d, Xmap->siteSD[0][%d]=%0.6f)\n",
			   t,SDx,varX,p->orientation,index,Xmap->siteSD[0][index]);
		  }
		}

		VARsum += VARdiff[k] = varY - varX;
		if(VARsum < VARmin){
		  VARmin = VARsum;
		  kmin = k+1;
		}
		if(VERB>=3)
		  printf("   k=%d/%d:i=%d,j=%d,I=%d,J=%d:varX=%0.6e,varY=%0.6e,VARsum=%0.6e(VARmin=%0.6e,kmin=%d)\n",k,K,i,j,I,J,varX,varY,VARsum,VARmin,kmin);
	      }
	      delete [] VARdiff;
	      if(VERB/* HERE >=2*/ ){
		printf("    Changing crossover point from (Y[%d]=%0.3f,X[%d]=%0.3f) to (Y[%d]=%0.3f,X[%d]=%0.3f)\n",
		       I,Ymap->site[0][I], p->orientation ? (M+1-J) : J, Xmap->site[0][p->orientation ? (M+1-J) : J], 
		       p->sites1[kmin], Ymap->site[0][p->sites1[kmin]],
		       p->orientation ? (M+1-p->sites2[kmin]) : p->sites2[kmin],    Xmap->site[0][p->orientation ? (M+1-p->sites2[kmin]) : p->sites2[kmin]]);
		fflush(stdout);
	      }
	      I = p->sites1[kmin];
	      J = p->sites2[kmin];
	    }

	    /* compute number of sites and allocate memory */
	    register int NumSites = I + (p->orientation ? (M+1-J)-1 : M - J);
	    for(register int c = 0; c < colors; c++){
	      if(pmap->site[c]){
		if(!pmap->blockmem) delete [] pmap->site[c];
		if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
		if(CmapChimQuality >= 1){
		  if(pmap->ChimQuality[c]) delete [] pmap->ChimQuality[c];
		  if(CmapChimQuality >= 2){
		    if(pmap->ChimNorm[c]) delete [] pmap->ChimNorm[c];
		    if(pmap->SegDupL[c]) delete [] pmap->SegDupL[c];
		    if(pmap->SegDupR[c]) delete [] pmap->SegDupR[c];
		    if(pmap->FragileEndL[c]) delete [] pmap->FragileEndL[c];
		    if(pmap->FragileEndR[c]) delete [] pmap->FragileEndR[c];
		    if(pmap->OutlierFrac[c]) delete [] pmap->OutlierFrac[c];
		    if(CmapChimQuality >= 3){
		      if(pmap->FragSd[c]) delete [] pmap->FragSd[c];
		      if(pmap->ExpSd[c]) delete [] pmap->ExpSd[c];
		      if(pmap->FragCov[c]) delete [] pmap->FragCov[c];
		      if(pmap->FragChiSq[c]) delete [] pmap->FragChiSq[c];
		    }
		  }
		}
	      }
	      if(pmap->SNRcnt[c]){
		if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		if(pmap->SNRdist[c]){
		  for(int i = 0; i <= pmap->numsite[0]+1; i++)
		    if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		  delete [] pmap->SNRdist[c];
		}
		delete [] pmap->SNRcnt[c];
	      }
	      if(pmap->Mask[c]){
		delete [] pmap->Mask[c];
		pmap->Mask[c] = NULL;
	      }

	      pmap->numsite[c] = NumSites;
	      pmap->site[c] = new FLOAT[NumSites+2];
	      pmap->siteSD[c] = new double[NumSites+2];
	      pmap->sitecov[c] = new float[NumSites+2];
	      pmap->sitecnt[c] = new float[NumSites+2];
	      if(CmapChimQuality >= 1){
		pmap->ChimQuality[c] = new float[NumSites+2];
		if(CmapChimQuality >= 2){
		  pmap->ChimNorm[c] = new float[NumSites+2];
		  pmap->SegDupL[c] = new float[NumSites+2];
		  pmap->SegDupR[c] = new float[NumSites+2];
		  pmap->FragileEndL[c] = new float[NumSites+2];
		  pmap->FragileEndR[c] = new float[NumSites+2];
		  pmap->OutlierFrac[c] = new float[NumSites+2];
		  if(CmapChimQuality >= 3){
		    pmap->FragSd[c] = new float[NumSites+2];
		    pmap->ExpSd[c] = new float[NumSites+2];
		    pmap->FragCov[c] = new float[NumSites+2];
		    pmap->FragChiSq[c] = new float[NumSites+2];
		  }
		}
	      }
	      pmap->SNRcnt[c] = 0;
	      if(Ymap->SNRcnt[c] || Xmap->SNRcnt[c]){
		pmap->SNRcnt[c] = new int[NumSites+2];
		pmap->SNRgmean[c] = new double[NumSites+2];
		pmap->lnSNRsd[c] = new double[NumSites+2];
		pmap->SNRdist[c] = new double*[NumSites+2];
	      }
	      if(Ymap->Mask[c] || Xmap->Mask[c]){
		pmap->Mask[c] = new size_t[NumSites+2];
		memset(pmap->Mask[c],0,(NumSites+2)*sizeof(size_t));
	      }

	      /* left end */
	      pmap->site[c][0] = 0.0;
	      pmap->siteSD[c][0] = 0.0;
	      pmap->sitecov[c][0] = 1;
	      pmap->sitecnt[c][0] = 1;
	      if(CmapChimQuality >= 1){
		pmap->ChimQuality[c][0] = END_CHIMQUAL;
		if(CmapChimQuality >= 2){
		  pmap->ChimNorm[c][0] = END_CHIMQUAL;
		  pmap->SegDupL[c][0] = END_CHIMQUAL;
		  pmap->SegDupR[c][0] = END_CHIMQUAL;
		  pmap->FragileEndL[c][0] = 0.0;
		  pmap->FragileEndR[c][0] = 0.0;
		  pmap->OutlierFrac[c][0] = 0.0;
		  if(CmapChimQuality >= 3){
		    pmap->FragSd[c][0] = 0.0;
		    pmap->ExpSd[c][0] = 0.0;
		    pmap->FragCov[c][0] = 0.0;
		    pmap->FragChiSq[c][0] = 1.0;
		  }
		}
	      }
	      if(pmap->SNRcnt[c]){
		pmap->SNRcnt[c][0] = 0;
		pmap->SNRgmean[c][0] = 0.0;
		pmap->lnSNRsd[c][0] = 0.0;
		pmap->SNRdist[c][0] = NULL;
	      }
	      if(pmap->Mask[c] && Ymap->Mask[c])
		pmap->Mask[c][1] |= Ymap->Mask[c][1];

	      /* right end */
	      pmap->siteSD[c][NumSites+1] = 0.0;
	      pmap->sitecov[c][NumSites+1] = 1;
	      pmap->sitecnt[c][NumSites+1] = 1;
	      if(CmapChimQuality >= 1){
		pmap->ChimQuality[c][NumSites+1] = END_CHIMQUAL;
		if(CmapChimQuality >= 2){
		  pmap->ChimNorm[c][NumSites+1] = END_CHIMQUAL;
		  pmap->SegDupL[c][NumSites+1] = END_CHIMQUAL;
		  pmap->SegDupR[c][NumSites+1] = END_CHIMQUAL;
		  pmap->FragileEndL[c][NumSites+1] = 0.0;
		  pmap->FragileEndR[c][NumSites+1] = 0.0;
		  pmap->OutlierFrac[c][NumSites+1] = 0.0;
		  if(CmapChimQuality >= 3){
		    pmap->FragSd[c][NumSites+1] = 0.0;
		    pmap->ExpSd[c][NumSites+1] = 0.0;
		    pmap->FragCov[c][NumSites+1] = 0.0;
		    pmap->FragChiSq[c][NumSites+1] = 1.0;
		  }
		}
	      }
	      if(pmap->SNRcnt[c]){
		pmap->SNRcnt[c][NumSites+1] = 0;
		pmap->SNRgmean[c][NumSites+1] = 0.0;
		pmap->lnSNRsd[c][NumSites+1] = 0.0;
		pmap->SNRdist[c][NumSites+1] = NULL;
	      }
	    }
	    pmap->blockmem = 0;
	    register FLOAT *Z = pmap->site[0];

	    /* copy left end of Ymap from site 0 to site I */
	    register int i = 0;
	    for(; i <= I; i++){
	      Z[i] = Y[i];
	      if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
	      pmap->siteSD[0][i] = Ymap->siteSD[0][i];
	      pmap->sitecov[0][i] = Ymap->sitecov[0][i];
	      pmap->sitecnt[0][i] = Ymap->sitecnt[0][i];
	      if(CmapChimQuality >= 1){
		pmap->ChimQuality[0][i] = Ymap->ChimQuality[0] ? Ymap->ChimQuality[0][i] : END_CHIMQUAL;
		if(CmapChimQuality >= 2){
		  pmap->ChimNorm[0][i] = Ymap->ChimNorm[0] ? Ymap->ChimNorm[0][i] : END_CHIMQUAL;
		  pmap->SegDupL[0][i] = Ymap->SegDupL[0] ? Ymap->SegDupL[0][i] : END_CHIMQUAL;
		  pmap->SegDupR[0][i] = Ymap->SegDupR[0] ? Ymap->SegDupR[0][i] : END_CHIMQUAL;
		  pmap->FragileEndL[0][i] = Ymap->FragileEndL[0] ? Ymap->FragileEndL[0][i] : 0.0;
		  pmap->FragileEndR[0][i] = Ymap->FragileEndR[0] ? Ymap->FragileEndR[0][i] : 0.0;
		  pmap->OutlierFrac[0][i] = Ymap->OutlierFrac[0] ? Ymap->OutlierFrac[0][i] : 0.0;
		  if(CmapChimQuality >= 3){
		    pmap->FragSd[0][i] = Ymap->FragSd[0] ? Ymap->FragSd[0][i] : 0.0;
		    pmap->ExpSd[0][i] = Ymap->ExpSd[0] ? Ymap->ExpSd[0][i] : 0.0;
		    pmap->FragCov[0][i] = Ymap->FragCov[0] ? Ymap->FragCov[0][i] : 0.0;
		    pmap->FragChiSq[0][i] = Ymap->FragChiSq[0] ? Ymap->FragChiSq[0][i] : 1.0;
		  }
		}
	      }
	      if(Ymap->SNRcnt[0]){
		pmap->SNRcnt[0][i] = Ymap->SNRcnt[0][i];
		pmap->SNRgmean[0][i] = Ymap->SNRgmean[0][i];
		pmap->lnSNRsd[0][i] = Ymap->lnSNRsd[0][i];
		pmap->SNRdist[0][i] = NULL;
		if(pmap->SNRcnt[0][i]){
		  pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		  for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
		    pmap->SNRdist[0][i][t] = Ymap->SNRdist[0][i][t];
		}
	      } else if(pmap->SNRcnt[0]){
		pmap->SNRcnt[0][i] = 0;
		pmap->SNRgmean[0][i] = 0.0;
		pmap->lnSNRsd[0][i] = 0.0;
		pmap->SNRdist[0][i] = NULL;
	      }
	      if(VERB>=3){
		printf("Copying site Y[i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
		       i,N, Y[i],Ymap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		fflush(stdout);
	      }
	    }

	    if(p->orientation){/* copy right end of Xmap from M+1-J downto 0 */
	      register int k = M+1-J;
	      register double Xstart = Z[i-1] + X[k];
	      for(register int j = k; --j >= 0; i++){
		Z[i] = Xstart - X[j];
		if(DEBUG) assert(Z[i] >= Z[i-1]);
		pmap->siteSD[0][i-1] = Xmap->siteSD[0][j];/* see definition of siteSD[] */
		pmap->sitecov[0][i] = Xmap->sitecov[0][j];
		pmap->sitecnt[0][i] = Xmap->sitecnt[0][j];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][j] : END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][j] : END_CHIMQUAL;
		    pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][j] : END_CHIMQUAL;
		    pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][j] : END_CHIMQUAL;
		    pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][j] : 0.0;
		    pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][j] : 0.0;
		    pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][j] : 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][j] : 0.0;
		      pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][j] : 0.0;
		      pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][j] : 0.0;
		      pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][j] : 1.0;
		    }
		  }
		}
		if(Xmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][j];
		  pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][j];
		  pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][j];
		  pmap->SNRdist[0][i] = NULL;
		  if(pmap->SNRcnt[0][i]){
		    pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		    for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
		      pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][j][t];
		  }
		} else if(pmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = 0;
		  pmap->SNRgmean[0][i] = 0.0;
		  pmap->lnSNRsd[0][i] = 0.0;
		  pmap->SNRdist[0][i] = NULL;
		}
		if(VERB>=3){
		  printf("Copying site X[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			 j, M, X[j],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		  fflush(stdout);
		}
	      }
	      if(DEBUG && !(i-1 == NumSites+1)){
		printf("N=%d,M=%d,N=%d,M=%d,I=%d,J=%d,i=%d,NumSites=%d,orientation=%d\n",
		       N,M,N,M,I,J,i,NumSites,p->orientation);
		fflush(stdout);
		assert(i-1 == NumSites+1);
	      }
	      if(pmap->Mask[0] && Xmap->Mask[0])
		pmap->Mask[0][NumSites+1] |= Xmap->Mask[0][1];
	    } else {/* copy right end of Xmap from  J to M+1 */

	      register int k = J;
	      register double Xstart = Z[i-1] - X[k];
	      for(register int j = k; ++j <= M+1; i++){
		Z[i] = Xstart + X[j];
		if(DEBUG) assert(Z[i] >= Z[i-1]);
		pmap->siteSD[0][i-1] = Xmap->siteSD[0][j-1];/* see definition of siteSD[] */
		pmap->sitecov[0][i] = Xmap->sitecov[0][j];
		pmap->sitecnt[0][i] = Xmap->sitecnt[0][j];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][j] : END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][j] : END_CHIMQUAL;
		    pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][j] : END_CHIMQUAL;
		    pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][j] : END_CHIMQUAL;
		    pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][j] : 0.0;
		    pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][j] : 0.0;
		    pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][j] : 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][j] : 0.0;
		      pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][j] : 0.0;
		      pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][j] : 0.0;
		      pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][j] : 1.0;
		    }
		  }
		}
		if(Xmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][j];
		  pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][j];
		  pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][j];
		  pmap->SNRdist[0][i] = NULL;
		  if(pmap->SNRcnt[0][i]){
		    pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		    for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
		      pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][j][t];
		  }
		} else if(pmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = 0;
		  pmap->SNRgmean[0][i] = 0.0;
		  pmap->lnSNRsd[0][i] = 0.0;
		  pmap->SNRdist[0][i] = NULL;
		}
		if(VERB >=3){
		  printf("Copying site X[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			 j, M, X[j],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		  fflush(stdout);
		}
	      }
	      if(DEBUG) assert(i-1 == NumSites+1);
	      if(pmap->Mask[0] && Xmap->Mask[0])
		pmap->Mask[0][NumSites+1] |= Xmap->Mask[0][M+1];
	    }
	  }
	} else { /* left end of X comes first */

	  if(distanceYtoX >= centerY - centerX){/* right end of X comes Last */
	  L_Y_inside_X:
	    int xstart = 0, xend = 0, ystart = 0, yend = 0;
	    if(PairMergePreferNGS){ /* check if Y has lower SD in some region of the alignment, if so copy Y in that aligned region */
	      double *VARdiff = new double[K];
	      for(int k = 0; k < K-1; k++){
		int i = p->sites1[k];
		int j = p->sites2[k];
		int I = p->sites1[k+1];
		int J = p->sites2[k+1];
		double varY = 0, varX = 0.0;
		for(int t = i; t < I; t++){
		  double SDy = Ymap->siteSD[0][t];
		  varY += SDy*SDy;
		}
		for(int t = j; t < J; t++){
		  double SDx = p->orientation ? Xmap->siteSD[0][M-t] : Xmap->siteSD[0][t];
		  varX += SDx*SDx;
		}
		VARdiff[k] = varX - varY;
	      }
	      /* Now locate interval of Y (ystart .. yend) with greatest reduction in variance over X (largest sum of Vardiff[k = 0 .. K-2])*/
	      double sum = 0.0;/* cumulative sum of VARdiff[0..k] */
	      double sumLWM = 0.0;/*  Low Water Mark of cumulative sum of VARdiff[0..k] */
	      int kLWM = 0;/* k value of Low Water Mark */
	      double maxdiff = 0.0;/* largest value of sum-sumLWM */
	      for(int k = 0; k < K-1; k++){
		sum += VARdiff[k];
		if(sum < sumLWM){
		  sumLWM = sum;
		  kLWM = k;
		}
		if(sum - sumLWM > maxdiff){
		  maxdiff = sum - sumLWM;
		  ystart = p->sites1[kLWM];
		  xstart = p->sites2[kLWM];
		  yend = p->sites1[k+1];
		  xend = p->sites2[k+1];
		}
	      }
	      if(VERB/* HERE >=2 */ && ystart < yend){
		printf("Retaining part of overlapped Y from Y[%d]=%0.3f to Y[%d]=%0.3f (replacing X[%d]=%0.3f to X[%d]=%0.3f\n",
		       ystart, Ymap->site[0][ystart], yend, Ymap->site[0][yend],
		       xstart, Xmap->site[0][p->orientation ? (M+1-xstart) : xstart],
		       xend, Xmap->site[0][p->orientation ? (M+1-xend) : xend]);
		fflush(stdout);
	      }
	      delete [] VARdiff;
	    }

	    if(ystart < yend){/* copy X from left end to xstart, then Y from ystart to yend, then X from xend to right end */
	      /* compute number of sites and allocate memory */
	      int NumSites = M - (xend - xstart) + (yend - ystart);
	      for(register int c = 0; c < colors; c++){
		if(pmap->site[c]){
		  if(!pmap->blockmem) delete [] pmap->site[c];
		  if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		  if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		  if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
		  if(CmapChimQuality >= 1){
		    if(pmap->ChimQuality[c]) delete [] pmap->ChimQuality[c];
		    if(CmapChimQuality >= 2){
		      if(pmap->ChimNorm[c]) delete [] pmap->ChimNorm[c];
		      if(pmap->SegDupL[c]) delete [] pmap->SegDupL[c];
		      if(pmap->SegDupR[c]) delete [] pmap->SegDupR[c];
		      if(pmap->FragileEndL[c]) delete [] pmap->FragileEndL[c];
		      if(pmap->FragileEndR[c]) delete [] pmap->FragileEndR[c];
		      if(pmap->OutlierFrac[c]) delete [] pmap->OutlierFrac[c];
		      if(CmapChimQuality >= 3){
			if(pmap->FragSd[c]) delete [] pmap->FragSd[c];
			if(pmap->ExpSd[c]) delete [] pmap->ExpSd[c];
			if(pmap->FragCov[c]) delete [] pmap->FragCov[c];
			if(pmap->FragChiSq[c]) delete [] pmap->FragChiSq[c];
		      }
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		  if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		  if(pmap->SNRdist[c]){
		    for(int i = 0; i <= pmap->numsite[0]+1; i++)
		      if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		    delete [] pmap->SNRdist[c];
		  }
		  delete [] pmap->SNRcnt[c];
		}
		if(pmap->Mask[c]){
		  delete [] pmap->Mask[c];
		  pmap->Mask[c] = NULL;
		}

		pmap->numsite[c] = NumSites;
		pmap->site[c] = new FLOAT[NumSites+2];
		pmap->siteSD[c] = new double[NumSites+2];
		pmap->sitecov[c] = new float[NumSites+2];
		pmap->sitecnt[c] = new float[NumSites+2];
		if(CmapChimQuality>=1){
		  pmap->ChimQuality[c] = new float[NumSites+2];
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c] = new float[NumSites+2];
		    pmap->SegDupL[c] = new float[NumSites+2];
		    pmap->SegDupR[c] = new float[NumSites+2];
		    pmap->FragileEndL[c] = new float[NumSites+2];
		    pmap->FragileEndR[c] = new float[NumSites+2];
		    pmap->OutlierFrac[c] = new float[NumSites+2];
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c] = new float[NumSites+2];
		      pmap->ExpSd[c] = new float[NumSites+2];
		      pmap->FragCov[c] = new float[NumSites+2];
		      pmap->FragChiSq[c] = new float[NumSites+2];
		    }
		  }
		}
		pmap->SNRcnt[c] = 0;
		if(Ymap->SNRcnt[c] || Xmap->SNRcnt[c]){
		  pmap->SNRcnt[c] = new int[NumSites+2];
		  pmap->SNRgmean[c] = new double[NumSites+2];
		  pmap->lnSNRsd[c] = new double[NumSites+2];
		  pmap->SNRdist[c] = new double*[NumSites+2];
		}
		if(Ymap->Mask[c] || Xmap->Mask[c]){
		  pmap->Mask[c] = new size_t[NumSites+2];
		  memset(pmap->Mask[c],0,(NumSites+2)*sizeof(size_t));
		}

		/* left end */
		pmap->site[c][0] = 0.0;
		pmap->siteSD[c][0] = 0.0;
		pmap->sitecov[c][0] = 1;
		pmap->sitecnt[c][0] = 1;
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c][0] = END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c][0] = END_CHIMQUAL;
		    pmap->SegDupL[c][0] = END_CHIMQUAL;
		    pmap->SegDupR[c][0] = END_CHIMQUAL;
		    pmap->FragileEndL[c][0] = 0.0;
		    pmap->FragileEndR[c][0] = 0.0;
		    pmap->OutlierFrac[c][0] = 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c][0] = 0.0;
		      pmap->ExpSd[c][0] = 0.0;
		      pmap->FragCov[c][0] = 0.0;
		      pmap->FragChiSq[c][0] = 1.0;
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][0] = 0;
		  pmap->SNRgmean[c][0] = 0.0;
		  pmap->lnSNRsd[c][0] = 0.0;
		  pmap->SNRdist[c][0] = NULL;
		}

		/* right end */
		pmap->siteSD[c][NumSites+1] = 0.0;
		pmap->sitecov[c][NumSites+1] = 1;
		pmap->sitecnt[c][NumSites+1] = 1;
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c][NumSites+1] = END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c][NumSites+1] = END_CHIMQUAL;
		    pmap->SegDupL[c][NumSites+1] = END_CHIMQUAL;
		    pmap->SegDupR[c][NumSites+1] = END_CHIMQUAL;
		    pmap->FragileEndL[c][NumSites+1] = 0.0;
		    pmap->FragileEndR[c][NumSites+1] = 0.0;
		    pmap->OutlierFrac[c][NumSites+1] = 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c][NumSites+1] = 0.0;
		      pmap->ExpSd[c][NumSites+1] = 0.0;
		      pmap->FragCov[c][NumSites+1] = 0.0;
		      pmap->FragChiSq[c][NumSites+1] = 1.0;
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][NumSites+1] = 0;
		  pmap->SNRgmean[c][NumSites+1] = 0.0;
		  pmap->lnSNRsd[c][NumSites+1] = 0.0;
		  pmap->SNRdist[c][NumSites+1] = NULL;
		}
	      }
	      pmap->blockmem = 0;
	      FLOAT *Z = pmap->site[0];
	    
	      int i = 0;
	      if(p->orientation){/* reversed X orientation */

		/* copy left end of Xmap from site M+1 downto M+1-xstart */
		double Xoffset = X[M+1];
		for(register int j = M+1; j >= M+1-xstart; j--, i++){
		  Z[i] = Xoffset - X[j];
		  if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[0][i] = Xmap->siteSD[0][j-1]; /* see definition of siteSD[] */
		  pmap->sitecov[0][i] = Xmap->sitecov[0][j];
		  pmap->sitecnt[0][i] = Xmap->sitecnt[0][j];
		  if(CmapChimQuality >= 1){
		    pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][j] : END_CHIMQUAL;
		    if(CmapChimQuality >= 2){
		      pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][j] : END_CHIMQUAL;
		      pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][j] : END_CHIMQUAL;
		      pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][j] : END_CHIMQUAL;
		      pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][j] : 0.0;
		      pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][j] : 0.0;
		      pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][j] : 0.0;
		      if(CmapChimQuality >= 3){
			pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][j] : 0.0;
			pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][j] : 0.0;
			pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][j] : 0.0;
			pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][j] : 1.0;
		      }
		    }
		  }
		  if(Xmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][j];
		    pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][j];
		    pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][j];
		    pmap->SNRdist[0][i] = NULL;
		    if(pmap->SNRcnt[0][i]){
		      pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		      for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
			pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][j][t];
		    }
		  } else if(pmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = 0;
		    pmap->SNRgmean[0][i] = 0.0;
		    pmap->lnSNRsd[0][i] = 0.0;
		    pmap->SNRdist[0][i] = NULL;
		  }
		  if(VERB>=3){
		    printf("Copying site X[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   j, M, X[j],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == (M+1)-(M+1-xstart));
		if(pmap->Mask[0] && Xmap->Mask[0])
		  pmap->Mask[0][1] |= Xmap->Mask[0][M+1];

		/* copy Ymap from site ystart to yend */
		int j = ystart;
		double Yoffset = Z[i-1] - Y[j];
		for(; ++j <= yend; i++){
		  Z[i] = Yoffset + Y[j];
		  if(DEBUG) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[0][i-1] = Ymap->siteSD[0][j-1];/* see definition of siteSD[] */
		  pmap->sitecov[0][i] = Ymap->sitecov[0][j];
		  pmap->sitecnt[0][i] = Ymap->sitecnt[0][j];
		  if(CmapChimQuality >= 1){
		    pmap->ChimQuality[0][i] = Ymap->ChimQuality[0] ? Ymap->ChimQuality[0][j] : END_CHIMQUAL;
		    if(CmapChimQuality >= 2){
		      pmap->ChimNorm[0][i] = Ymap->ChimNorm[0] ? Ymap->ChimNorm[0][j] : END_CHIMQUAL;
		      pmap->SegDupL[0][i] = Ymap->SegDupL[0] ? Ymap->SegDupL[0][j] : END_CHIMQUAL;
		      pmap->SegDupR[0][i] = Ymap->SegDupR[0] ? Ymap->SegDupR[0][j] : END_CHIMQUAL;
		      pmap->FragileEndL[0][i] = Ymap->FragileEndL[0] ? Ymap->FragileEndL[0][j] : 0.0;
		      pmap->FragileEndR[0][i] = Ymap->FragileEndR[0] ? Ymap->FragileEndR[0][j] : 0.0;
		      pmap->OutlierFrac[0][i] = Ymap->OutlierFrac[0] ? Ymap->OutlierFrac[0][j] : 0.0;
		      if(CmapChimQuality >= 3){
			pmap->FragSd[0][i] = Ymap->FragSd[0] ? Ymap->FragSd[0][j] : 0.0;
			pmap->ExpSd[0][i] = Ymap->ExpSd[0] ? Ymap->ExpSd[0][j] : 0.0;
			pmap->FragCov[0][i] = Ymap->FragCov[0] ? Ymap->FragCov[0][j] : 0.0;
			pmap->FragChiSq[0][i] = Ymap->FragChiSq[0] ? Ymap->FragChiSq[0][j] : 1.0;
		      }
		    }
		  }
		  if(Ymap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = Ymap->SNRcnt[0][j];
		    pmap->SNRgmean[0][i] = Ymap->SNRgmean[0][j];
		    pmap->lnSNRsd[0][i] = Ymap->lnSNRsd[0][j];
		    pmap->SNRdist[0][i] = NULL;
		    if(pmap->SNRcnt[0][i]){
		      pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		      for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
			pmap->SNRdist[0][i][t] = Ymap->SNRdist[0][j][t];
		    }
		  } else if(pmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = 0;
		    pmap->SNRgmean[0][i] = 0.0;
		    pmap->lnSNRsd[0][i] = 0.0;
		    pmap->SNRdist[0][i] = NULL;
		  }
		  if(VERB>=3){
		    printf("Copying site Y[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   j,N, Y[j],Ymap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == (M+1)-(M+1-xstart) + (yend-ystart));

		/* copy right end of Xmap from site M+1-xend downto 0 */
		int k = M+1-xend;
		Xoffset = Z[i-1] + X[k];
		for(register int j = k; --j >= 0; i++){
		  Z[i] = Xoffset - X[j];
		  if(DEBUG) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[0][i-1] = Xmap->siteSD[0][j];/* see definition of siteSD[] */
		  pmap->sitecov[0][i] = Xmap->sitecov[0][j];
		  pmap->sitecnt[0][i] = Xmap->sitecnt[0][j];
		  if(CmapChimQuality >= 1){
		    pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][j] : END_CHIMQUAL;
		    if(CmapChimQuality >= 2){
		      pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][j] : END_CHIMQUAL;
		      pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][j] : END_CHIMQUAL;
		      pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][j] : END_CHIMQUAL;
		      pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][j] : 0.0;
		      pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][j] : 0.0;
		      pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][j] : 0.0;
		      if(CmapChimQuality >= 3){
			pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][j] : 0.0;
			pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][j] : 0.0;
			pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][j] : 0.0;
			pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][j] : 1.0;
		      }
		    }
		  }
		  if(Xmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][j];
		    pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][j];
		    pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][j];
		    pmap->SNRdist[0][i] = NULL;
		    if(pmap->SNRcnt[0][i]){
		      pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		      for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
			pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][j][t];
		    }
		  } else if(pmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = 0;
		    pmap->SNRgmean[0][i] = 0.0;
		    pmap->lnSNRsd[0][i] = 0.0;
		    pmap->SNRdist[0][i] = NULL;
		  }
		  if(VERB>=3){
		    printf("Copying site X[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   j, M, X[j],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG && !(i-1 == NumSites+1)){
		  printf("i=%d,NumSites=%d,M=%d,xstart=%d,xend=%d,ystart=%d,yend=%d\n",i,NumSites,M,xstart,xend,ystart,yend);
		  fflush(stdout);
		  assert(i-1 == NumSites+1);
		}
		if(pmap->Mask[0] && Xmap->Mask[0])
		  pmap->Mask[0][NumSites+1] |= Xmap->Mask[0][1];

	      } else {/* normal X orientation */

		/* copy left end of Xmap from site 0 to xstart */
		if(pmap->Mask[0] && Xmap->Mask[0])
		  pmap->Mask[0][1] |= Xmap->Mask[0][1];

		for(; i <= xstart; i++){
		  Z[i] = X[i];
		  if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[0][i] = Xmap->siteSD[0][i];
		  pmap->sitecov[0][i] = Xmap->sitecov[0][i];
		  pmap->sitecnt[0][i] = Xmap->sitecnt[0][i];
		  if(CmapChimQuality >= 1){
		    pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][i] : END_CHIMQUAL;
		    if(CmapChimQuality >= 2){
		      pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][i] : END_CHIMQUAL;
		      pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][i] : END_CHIMQUAL;
		      pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][i] : END_CHIMQUAL;
		      pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][i] : 0.0;
		      pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][i] : 0.0;
		      pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][i] : 0.0;
		      if(CmapChimQuality >= 3){
			pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][i] : 0.0;
			pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][i] : 0.0;
			pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][i] : 0.0;
			pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][i] : 1.0;
		      }
		    }
		  }
		  if(Xmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][i];
		    pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][i];
		    pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][i];
		    pmap->SNRdist[0][i] = NULL;
		    if(pmap->SNRcnt[0][i]){
		      pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		      for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
			pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][i][t];
		    }
		  } else if(pmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = 0;
		    pmap->SNRgmean[0][i] = 0.0;
		    pmap->lnSNRsd[0][i] = 0.0;
		    pmap->SNRdist[0][i] = NULL;
		  }
		  if(VERB>=3){
		    printf("Copying site X[i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   i, M, X[i],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == xstart);
	      
		/* copy Ymap from site ystart to yend */
		int j = ystart;
		double Yoffset = Z[i-1] - Y[j];
		for(; ++j <= yend; i++){
		  Z[i] = Yoffset + Y[j];
		  if(DEBUG) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[0][i-1] = Ymap->siteSD[0][j-1];/* see definition of siteSD[] */
		  pmap->sitecov[0][i] = Ymap->sitecov[0][j];
		  pmap->sitecnt[0][i] = Ymap->sitecnt[0][j];
		  if(CmapChimQuality >= 1){
		    pmap->ChimQuality[0][i] = Ymap->ChimQuality[0] ? Ymap->ChimQuality[0][j] : END_CHIMQUAL;
		    if(CmapChimQuality >= 2){
		      pmap->ChimNorm[0][i] = Ymap->ChimNorm[0] ? Ymap->ChimNorm[0][j] : END_CHIMQUAL;
		      pmap->SegDupL[0][i] = Ymap->SegDupL[0] ? Ymap->SegDupL[0][j] : END_CHIMQUAL;
		      pmap->SegDupR[0][i] = Ymap->SegDupR[0] ? Ymap->SegDupR[0][j] : END_CHIMQUAL;
		      pmap->FragileEndL[0][i] = Ymap->FragileEndL[0] ? Ymap->FragileEndL[0][j] : 0.0;
		      pmap->FragileEndR[0][i] = Ymap->FragileEndR[0] ? Ymap->FragileEndR[0][j] : 0.0;
		      pmap->OutlierFrac[0][i] = Ymap->OutlierFrac[0] ? Ymap->OutlierFrac[0][j] : 0.0;
		      if(CmapChimQuality >= 3){
			pmap->FragSd[0][i] = Ymap->FragSd[0] ? Ymap->FragSd[0][j] : 0.0;
			pmap->ExpSd[0][i] = Ymap->ExpSd[0] ? Ymap->ExpSd[0][j] : 0.0;
			pmap->FragCov[0][i] = Ymap->FragCov[0] ? Ymap->FragCov[0][j] : 0.0;
			pmap->FragChiSq[0][i] = Ymap->FragChiSq[0] ? Ymap->FragChiSq[0][j] : 1.0;
		      }
		    }
		  }
		  if(Ymap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = Ymap->SNRcnt[0][j];
		    pmap->SNRgmean[0][i] = Ymap->SNRgmean[0][j];
		    pmap->lnSNRsd[0][i] = Ymap->lnSNRsd[0][j];
		    pmap->SNRdist[0][i] = NULL;
		    if(pmap->SNRcnt[0][i]){
		      pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		      for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
			pmap->SNRdist[0][i][t] = Ymap->SNRdist[0][j][t];
		    }
		  } else if(pmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = 0;
		    pmap->SNRgmean[0][i] = 0.0;
		    pmap->lnSNRsd[0][i] = 0.0;
		    pmap->SNRdist[0][i] = NULL;
		  }
		  if(VERB>=3){
		    printf("Copying site Y[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   j,N, Y[j],Ymap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == xstart + (yend-ystart));
	      
		/* copy left end of Xmap from site xend to M+1 */
		int k = xend;
		double Xoffset = Z[i-1] - X[k];
		for(register int j = k; ++j <= M+1; i++){
		  Z[i] = Xoffset + X[j];
		  if(DEBUG) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[0][i-1] = Xmap->siteSD[0][j-1];/* see definition of siteSD[] */
		  pmap->sitecov[0][i] = Xmap->sitecov[0][j];
		  pmap->sitecnt[0][i] = Xmap->sitecnt[0][j];
		  if(CmapChimQuality >= 1){
		    pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][j] : END_CHIMQUAL;
		    if(CmapChimQuality >= 2){
		      pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][j] : END_CHIMQUAL;
		      pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][j] : END_CHIMQUAL;
		      pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][j] : END_CHIMQUAL;
		      pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][j] : 0.0;
		      pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][j] : 0.0;
		      pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][j] : 0.0;
		      if(CmapChimQuality >= 3){
			pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][j] : 0.0;
			pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][j] : 0.0;
			pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][j] : 0.0;
			pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][j] : 1.0;
		      }
		    }
		  }
		  if(Xmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][j];
		    pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][j];
		    pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][j];
		    pmap->SNRdist[0][i] = NULL;
		    if(pmap->SNRcnt[0][i]){
		      pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		      for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
			pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][j][t];
		    }
		  } else if(pmap->SNRcnt[0]){
		    pmap->SNRcnt[0][i] = 0;
		    pmap->SNRgmean[0][i] = 0.0;
		    pmap->lnSNRsd[0][i] = 0.0;
		    pmap->SNRdist[0][i] = NULL;
		  }
		  if(VERB >=3){
		    printf("Copying site X[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   j, M, X[j],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == NumSites+1);
		if(pmap->Mask[0] && Xmap->Mask[0])
		  pmap->Mask[0][NumSites+1] |= Xmap->Mask[0][M+1];
	      } /* normal X orientation */
	    } else {/* just copy X */

	      if(PairMergeRepeat){ // new code : instead of copying X, just restore original X
		if(VERB){
		  printf("%d:Created merged map (id=%lld,len=%0.3f) based on alignment between Y mapid=%d(id=%lld,len=%0.3f) and X mapid=%d(id=%lld,len=%0.3f):mapid=%d,id=%lld,len=%0.3f,sites=%d (Just keeping Xmap=%p)\n",
			 pairmergeIter,Xmap->id, Xmap->site[0][Xmap->numsite[0]+1], Ymap->mapid,Ymap->id, Ymap->site[0][Ymap->numsite[0]+1], 
			 Xmap->mapid,Xmap->id, Xmap->site[0][Xmap->numsite[0]+1],Xmap->mapid,Xmap->id,Xmap->len,Xmap->numsite[0], Xmap);
		  /*		printf("%d:Restoring Xmap(id=%lld,mapid=%d) discarding Ymap(id=%lld,mapid=%d)\n",pairmergeIter,Xmap->id,Xmap->mapid,Ymap->id,Ymap->mapid);*/
		  fflush(stdout);
		}
		Xmap->paired = 0;
		if(DEBUG) assert(Xmap->mapid == p->mapid2);
		Ymap->paired = 1 + Xmap->mapid;
		nummaps--;

		continue;/* next Gmap[i] */
	      }

	      /* compute number of sites and allocate memory : NOTE Cmap->orig* fields are not changed and become invalid/undefined, flagged with centercloned=1 */
	      register int NumSites = M;
	      for(register int c = 0; c < colors; c++){
		if(pmap->site[c]){
		  if(!pmap->blockmem) delete [] pmap->site[c];
		  if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		  if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		  if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
		  if(CmapChimQuality >= 1){
		    if(pmap->ChimQuality[c]) delete [] pmap->ChimQuality[c];
		    if(CmapChimQuality >= 2){
		      if(pmap->ChimNorm[c]) delete [] pmap->ChimNorm[c];
		      if(pmap->SegDupL[c]) delete [] pmap->SegDupL[c];
		      if(pmap->SegDupR[c]) delete [] pmap->SegDupR[c];
		      if(pmap->FragileEndL[c]) delete [] pmap->FragileEndL[c];
		      if(pmap->FragileEndR[c]) delete [] pmap->FragileEndR[c];
		      if(pmap->OutlierFrac[c]) delete [] pmap->OutlierFrac[c];
		      if(CmapChimQuality >= 3){
			if(pmap->FragSd[c]) delete [] pmap->FragSd[c];
			if(pmap->ExpSd[c]) delete [] pmap->ExpSd[c];
			if(pmap->FragCov[c]) delete [] pmap->FragCov[c];
			if(pmap->FragChiSq[c]) delete [] pmap->FragChiSq[c];
		      }
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		  if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		  if(pmap->SNRdist[c]){
		    for(int i = 0; i <= pmap->numsite[0]+1; i++)
		      if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		    delete [] pmap->SNRdist[c];
		  }
		  delete [] pmap->SNRcnt[c];
		}
		if(pmap->Mask[c]){
		  delete [] pmap->Mask[c];
		  pmap->Mask[c] = NULL;
		}

		pmap->numsite[c] = NumSites;
		pmap->site[c] = new FLOAT[NumSites+2];
		pmap->siteSD[c] = new double[NumSites+2];
		pmap->sitecov[c] = new float[NumSites+2];
		pmap->sitecnt[c] = new float[NumSites+2];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c] = new float[NumSites+2];
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c] = new float[NumSites+2];
		    pmap->SegDupL[c] = new float[NumSites+2];
		    pmap->SegDupR[c] = new float[NumSites+2];
		    pmap->FragileEndL[c] = new float[NumSites+2];
		    pmap->FragileEndR[c] = new float[NumSites+2];
		    pmap->OutlierFrac[c] = new float[NumSites+2];
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c] = new float[NumSites+2];
		      pmap->ExpSd[c] = new float[NumSites+2];
		      pmap->FragCov[c] = new float[NumSites+2];
		      pmap->FragChiSq[c] = new float[NumSites+2];
		    }
		  }
		}
		pmap->SNRcnt[c] = 0;
		if(Xmap->SNRcnt[c]){
		  pmap->SNRcnt[c] = new int[NumSites+2];
		  pmap->SNRgmean[c] = new double[NumSites+2];
		  pmap->lnSNRsd[c] = new double[NumSites+2];
		  pmap->SNRdist[c] = new double*[NumSites+2];
		}
		if(Xmap->Mask[c]){
		  pmap->Mask[c] = new size_t[NumSites+2];
		  memcpy(pmap->Mask[c],Xmap->Mask[c],(NumSites+2)*sizeof(size_t));
		}

		/* left end */
		pmap->site[c][0] = 0.0;
		pmap->siteSD[c][0] = 0.0;
		pmap->sitecov[c][0] = 1;
		pmap->sitecnt[c][0] = 1;
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c][0] = END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c][0] = END_CHIMQUAL;
		    pmap->SegDupL[c][0] = END_CHIMQUAL;
		    pmap->SegDupR[c][0] = END_CHIMQUAL;
		    pmap->FragileEndL[c][0] = 0.0;
		    pmap->FragileEndR[c][0] = 0.0;
		    pmap->OutlierFrac[c][0] = 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c][0] = 0.0;
		      pmap->ExpSd[c][0] = 0.0;
		      pmap->FragCov[c][0] = 0.0;
		      pmap->FragChiSq[c][0] = 1.0;
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][0] = 0;
		  pmap->SNRgmean[c][0] = 0.0;
		  pmap->lnSNRsd[c][0] = 0.0;
		  pmap->SNRdist[c][0] = NULL;
		}

		/* right end */
		pmap->siteSD[c][NumSites+1] = 0.0;
		pmap->sitecov[c][NumSites+1] = 1;
		pmap->sitecnt[c][NumSites+1] = 1;
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[c][NumSites+1] = END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[c][NumSites+1] = END_CHIMQUAL;
		    pmap->SegDupL[c][NumSites+1] = END_CHIMQUAL;
		    pmap->SegDupR[c][NumSites+1] = END_CHIMQUAL;
		    pmap->FragileEndL[c][NumSites+1] = 0.0;
		    pmap->FragileEndR[c][NumSites+1] = 0.0;
		    pmap->OutlierFrac[c][NumSites+1] = 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[c][NumSites+1] = 0.0;
		      pmap->ExpSd[c][NumSites+1] = 0.0;
		      pmap->FragCov[c][NumSites+1] = 0.0;
		      pmap->FragChiSq[c][NumSites+1] = 1.0;
		    }
		  }
		}
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][NumSites+1] = 0;
		  pmap->SNRgmean[c][NumSites+1] = 0.0;
		  pmap->lnSNRsd[c][NumSites+1] = 0.0;
		  pmap->SNRdist[c][NumSites+1] = NULL;
		}
	      }
	      pmap->blockmem = 0;
	      register FLOAT *Z = pmap->site[0];

	      register int i = 0;
	      for(; i <= M+1; i++){
		Z[i] = X[i];
		if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		pmap->siteSD[0][i] = Xmap->siteSD[0][i];
		pmap->sitecov[0][i] = Xmap->sitecov[0][i];
		pmap->sitecnt[0][i] = Xmap->sitecnt[0][i];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][i] : END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][i] : END_CHIMQUAL;
		    pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][i] : END_CHIMQUAL;
		    pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][i] : END_CHIMQUAL;
		    pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][i] : 0.0;
		    pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][i] : 0.0;
		    pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][i] : 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][i] : 0.0;
		      pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][i] : 0.0;
		      pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][i] : 0.0;
		      pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][i] : 1.0;
		    }
		  }
		}
		if(pmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][i];
		  pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][i];
		  pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][i];
		  pmap->SNRdist[0][i] = 0;
		  if(pmap->SNRcnt[0][i]){
		    pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		    for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
		      pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][i][t];
		  }
		}
		if(VERB>=3){
		  printf("Copying site X[i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			 i, M, X[i],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		  fflush(stdout);
		}
	      }
	      if(DEBUG) assert(i-1 == NumSites+1);
	    }
	  } else { /* copy left end of X and right end of Y */

	    if(PairMergePreferNGS){/* location transition point that retains the most NGS sites (with StdDev==0) */
	      double *VARdiff = new double[K];
	      double VARmin = 0.0, VARsum = 0.0;
	      int kmin = 0;
	      for(int k = 0; k < K-1; k++){
		int i = p->sites1[k];
		int j = p->sites2[k];
		int I = p->sites1[k+1];
		int J = p->sites2[k+1];
		double varY = 0, varX = 0.0;
		for(int t = i; t < I; t++){
		  double SDy = Ymap->siteSD[0][t];
		  varY += SDy*SDy;
		  if(VERB>=4)
		    printf("       t=%d:SDy=%0.6f,varY=%0.6e (or=%d, Ymap->siteSD[0][%d]=%0.6f)\n",
			   t,SDy,varY,p->orientation,t,Ymap->siteSD[0][t]);
		}
		for(int t = j; t < J; t++){
		  double SDx = p->orientation ? Xmap->siteSD[0][M-t] : Xmap->siteSD[0][t];
		  varX += SDx*SDx;
		  if(VERB>=4){
		    int index = (p->orientation ? M-t : t);
		    printf("       t=%d:SDx=%0.6f,varX=%0.6e (or=%d, Xmap->siteSD[0][%d]=%0.6f)\n",
			   t,SDx,varX,p->orientation,index,Xmap->siteSD[0][index]);
		  }
		}

		VARsum += VARdiff[k] = varX - varY;
		if(VARsum < VARmin){
		  VARmin = VARsum;
		  kmin = k+1;
		}
		if(VERB>=3)
		  printf("   k=%d/%d:i=%d,j=%d,I=%d,J=%d:varX=%0.6e,varY=%0.6e,VARsum=%0.6e(VARmin=%0.6e,kmin=%d)\n",k,K,i,j,I,J,varX,varY,VARsum,VARmin,kmin);
	      }
	      delete [] VARdiff;
	      if(VERB/* HERE >=2 */){
		printf("Changing crossover point from (Y[%d]=%0.3f,X[%d]=%0.3f) to (Y[%d]=%0.3f,X[%d]=%0.3f)\n",
		       I,Ymap->site[0][I], p->orientation ? (M+1-J) : J, Xmap->site[0][p->orientation ? (M+1-J) : J], 
		       p->sites1[kmin], Ymap->site[0][p->sites1[kmin]],
		       p->orientation ? (M+1-p->sites2[kmin]) : p->sites2[kmin], Xmap->site[0][p->orientation ? (M+1-p->sites2[kmin]) : p->sites2[kmin]]);
		fflush(stdout);
	      }
	      I = p->sites1[kmin];
	      J = p->sites2[kmin];
	    }

	    /* compute number of sites and allocate memory */
	    int NumSites = J + N - I;
	    for(int c = 0; c < colors; c++){
	      if(pmap->site[c]){
		if(!pmap->blockmem) delete [] pmap->site[c];
		if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
		if(CmapChimQuality >= 1){
		  if(pmap->ChimQuality[c]) delete [] pmap->ChimQuality[c];
		  if(CmapChimQuality >= 2){
		    if(pmap->ChimNorm[c]) delete [] pmap->ChimNorm[c];
		    if(pmap->SegDupL[c]) delete [] pmap->SegDupL[c];
		    if(pmap->SegDupR[c]) delete [] pmap->SegDupR[c];
		    if(pmap->FragileEndL[c]) delete [] pmap->FragileEndL[c];
		    if(pmap->FragileEndR[c]) delete [] pmap->FragileEndR[c];
		    if(pmap->OutlierFrac[c]) delete [] pmap->OutlierFrac[c];
		    if(CmapChimQuality >= 3){
		      if(pmap->FragSd[c]) delete [] pmap->FragSd[c];
		      if(pmap->ExpSd[c]) delete [] pmap->ExpSd[c];
		      if(pmap->FragCov[c]) delete [] pmap->FragCov[c];
		      if(pmap->FragChiSq[c]) delete [] pmap->FragChiSq[c];
		    }
		  }
		}
	      }
	      if(pmap->SNRcnt[c]){
		if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		if(pmap->SNRdist[c]){
		  for(int i = 0; i <= pmap->numsite[0]+1; i++)
		    if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		  delete [] pmap->SNRdist[c];
		}
		delete [] pmap->SNRcnt[c];
	      }
	      if(pmap->Mask[c]){
		delete [] pmap->Mask[c];
		pmap->Mask[c] = NULL;
	      }

	      pmap->numsite[c] = NumSites;
	      pmap->site[c] = new FLOAT[NumSites+2];
	      pmap->siteSD[c] = new double[NumSites+2];
	      pmap->sitecov[c] = new float[NumSites+2];
	      pmap->sitecnt[c] = new float[NumSites+2];
	      if(CmapChimQuality >= 1){
		pmap->ChimQuality[c] = new float[NumSites+2];
		if(CmapChimQuality >= 2){
		  pmap->ChimNorm[c] = new float[NumSites+2];
		  pmap->SegDupL[c] = new float[NumSites+2];
		  pmap->SegDupR[c] = new float[NumSites+2];
		  pmap->FragileEndL[c] = new float[NumSites+2];
		  pmap->FragileEndR[c] = new float[NumSites+2];
		  pmap->OutlierFrac[c] = new float[NumSites+2];
		  if(CmapChimQuality >= 3){
		    pmap->FragSd[c] = new float[NumSites+2];
		    pmap->ExpSd[c] = new float[NumSites+2];
		    pmap->FragCov[c] = new float[NumSites+2];
		    pmap->FragChiSq[c] = new float[NumSites+2];
		  }
		}
	      }
	      pmap->SNRcnt[c] = 0;
	      if(Xmap->SNRcnt[c] || Ymap->SNRcnt[c]){
		pmap->SNRcnt[c] = new int[NumSites+2];
		pmap->SNRgmean[c] = new double[NumSites+2];
		pmap->lnSNRsd[c] = new double[NumSites+2];
		pmap->SNRdist[c] = new double*[NumSites+2];
	      }
	      if(Ymap->Mask[c] || Xmap->Mask[c]){
		pmap->Mask[c] = new size_t[NumSites+2];
		memset(pmap->Mask[c],0,(NumSites+2)*sizeof(size_t));
	      }

	      /* left end */
	      pmap->site[c][0] = 0.0;
	      pmap->siteSD[c][0] = 0.0;
	      pmap->sitecov[c][0] = 1;
	      pmap->sitecnt[c][0] = 1;
	      if(CmapChimQuality >= 1){
		pmap->ChimQuality[c][0] = END_CHIMQUAL;
		if(CmapChimQuality >= 2){
		  pmap->ChimNorm[c][0] = END_CHIMQUAL;
		  pmap->SegDupL[c][0] = END_CHIMQUAL;
		  pmap->SegDupR[c][0] = END_CHIMQUAL;
		  pmap->FragileEndL[c][0] = 0.0;
		  pmap->FragileEndR[c][0] = 0.0;
		  pmap->OutlierFrac[c][0] = 0.0;
		  if(CmapChimQuality >= 3){
		    pmap->FragSd[c][0] = 0.0;
		    pmap->ExpSd[c][0] = 0.0;
		    pmap->FragCov[c][0] = 0.0;
		    pmap->FragChiSq[c][0] = 1.0;
		  }
		}
	      }
	      if(pmap->SNRcnt[c]){
		pmap->SNRcnt[c][0] = 0;
		pmap->SNRgmean[c][0] = 0.0;
		pmap->lnSNRsd[c][0] = 0.0;
		pmap->SNRdist[c][0] = NULL;
	      }

	      /* right end */
	      pmap->siteSD[c][NumSites+1] = 0.0;
	      pmap->sitecov[c][NumSites+1] = 1;
	      pmap->sitecnt[c][NumSites+1] = 1;
	      if(CmapChimQuality >= 1){
		pmap->ChimQuality[c][NumSites+1] = END_CHIMQUAL;
		if(CmapChimQuality >= 2){
		  pmap->ChimNorm[c][NumSites+1] = END_CHIMQUAL;
		  pmap->SegDupL[c][NumSites+1] = END_CHIMQUAL;
		  pmap->SegDupR[c][NumSites+1] = END_CHIMQUAL;
		  pmap->FragileEndL[c][NumSites+1] = 0.0;
		  pmap->FragileEndR[c][NumSites+1] = 0.0;
		  pmap->OutlierFrac[c][NumSites+1] = 0.0;
		  if(CmapChimQuality >= 3){
		    pmap->FragSd[c][NumSites+1] = 0.0;
		    pmap->ExpSd[c][NumSites+1] = 0.0;
		    pmap->FragCov[c][NumSites+1] = 0.0;
		    pmap->FragChiSq[c][NumSites+1] = 1.0;
		  }
		}
	      }
	      if(pmap->SNRcnt[c]){
		pmap->SNRcnt[c][NumSites+1] = 0;
		pmap->SNRgmean[c][NumSites+1] = 0.0;
		pmap->lnSNRsd[c][NumSites+1] = 0.0;
		pmap->SNRdist[c][NumSites+1] = NULL;
	      }
	    }

	    register FLOAT *Z = pmap->site[0];
	    register int i = 0;
	    if(p->orientation){ /* copy left end of Xmap from site M+1 downto M+1-J */
	      if(Xmap->Mask[0])
		pmap->Mask[0][1] |= Xmap->Mask[0][M+1];

	      register double Xstart = X[M+1];
	      for(register int j = M+1; j >= M+1-J; j--, i++){
		Z[i] = Xstart - X[j];
		if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		pmap->siteSD[0][i] = Xmap->siteSD[0][j-1]; /* see definition of siteSD[] */
		pmap->sitecov[0][i] = Xmap->sitecov[0][j];
		pmap->sitecnt[0][i] = Xmap->sitecnt[0][j];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][j] : END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][j] : END_CHIMQUAL;
		    pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][j] : END_CHIMQUAL;
		    pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][j] : END_CHIMQUAL;
		    pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][j] : 0.0;
		    pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][j] : 0.0;
		    pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][j] : 0.0;
		    if(CmapChimQuality >= 3){
 		      pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][j] : 0.0;
 		      pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][j] : 0.0;
 		      pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][j] : 0.0;
 		      pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][j] : 1.0;
		    }
		  }
		}
		if(Xmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][j];
		  pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][j];
		  pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][j];
		  pmap->SNRdist[0][i] = NULL;
		  if(pmap->SNRcnt[0][i]){
		    pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		    for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
		      pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][j][t];
		  }
		} else if(pmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = 0;
		  pmap->SNRgmean[0][i] = 0.0;
		  pmap->lnSNRsd[0][i] = 0.0;
		  pmap->SNRdist[0][i] = NULL;
		}
		if(VERB>=3){
		  printf("Copying site X[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			 j, M, X[j],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		  fflush(stdout);
		}
	      }
	      if(DEBUG) assert(i-1 == (M+1)-(M+1-J));
	    } else { /* copy left end of Xmap from site 0 to J */

	      if(Xmap->Mask[0])
		pmap->Mask[0][1] |= Xmap->Mask[0][1];
	      for(; i <= J; i++){
		Z[i] = X[i];
		if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		pmap->siteSD[0][i] = Xmap->siteSD[0][i];
		pmap->sitecov[0][i] = Xmap->sitecov[0][i];
		pmap->sitecnt[0][i] = Xmap->sitecnt[0][i];
		if(CmapChimQuality >= 1){
		  pmap->ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][i] : END_CHIMQUAL;
		  if(CmapChimQuality >= 2){
		    pmap->ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][i] : END_CHIMQUAL;
		    pmap->SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][i] : END_CHIMQUAL;
		    pmap->SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][i] : END_CHIMQUAL;
		    pmap->FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][i] : 0.0;
		    pmap->FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][i] : 0.0;
		    pmap->OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][i] : 0.0;
		    if(CmapChimQuality >= 3){
		      pmap->FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][i] : 0.0;
		      pmap->ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][i] : 0.0;
		      pmap->FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][i] : 0.0;
		      pmap->FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][i] : 1.0;
		    }
		  }
		}
		if(Xmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = Xmap->SNRcnt[0][i];
		  pmap->SNRgmean[0][i] = Xmap->SNRgmean[0][i];
		  pmap->lnSNRsd[0][i] = Xmap->lnSNRsd[0][i];
		  pmap->SNRdist[0][i] = NULL;
		  if(pmap->SNRcnt[0][i]){
		    pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		    for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
		      pmap->SNRdist[0][i][t] = Xmap->SNRdist[0][i][t];
		  }
		} else if(pmap->SNRcnt[0]){
		  pmap->SNRcnt[0][i] = 0;
		  pmap->SNRgmean[0][i] = 0.0;
		  pmap->lnSNRsd[0][i] = 0.0;
		  pmap->SNRdist[0][i] = NULL;
		}
		if(VERB>=3){
		  printf("Copying site X[i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			 i, M, X[i],Xmap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		  fflush(stdout);
		}
	      }
	      if(DEBUG) assert(i-1 == J);
	    }

	    /* copy right end of Ymap from site I to N+1 */
	    register int j = I;
	    register double Ystart = Z[i-1] - Y[j];
	    for(; ++j <= N+1; i++){
	      Z[i] = Ystart + Y[j];
	      if(DEBUG) assert(Z[i] >= Z[i-1]);
	      pmap->siteSD[0][i-1] = Ymap->siteSD[0][j-1];/* see definition of siteSD[] */
	      pmap->sitecov[0][i] = Ymap->sitecov[0][j];
	      pmap->sitecnt[0][i] = Ymap->sitecnt[0][j];
	      if(CmapChimQuality >= 1){
		pmap->ChimQuality[0][i] = Ymap->ChimQuality[0] ? Ymap->ChimQuality[0][j] : END_CHIMQUAL;
		if(CmapChimQuality >= 2){
		  pmap->ChimNorm[0][i] = Ymap->ChimNorm[0] ? Ymap->ChimNorm[0][j] : END_CHIMQUAL;
		  pmap->SegDupL[0][i] = Ymap->SegDupL[0] ? Ymap->SegDupL[0][j] : END_CHIMQUAL;
		  pmap->SegDupR[0][i] = Ymap->SegDupR[0] ? Ymap->SegDupR[0][j] : END_CHIMQUAL;
		  pmap->FragileEndL[0][i] = Ymap->FragileEndL[0] ? Ymap->FragileEndL[0][j] : 0.0;
		  pmap->FragileEndR[0][i] = Ymap->FragileEndR[0] ? Ymap->FragileEndR[0][j] : 0.0;
		  pmap->OutlierFrac[0][i] = Ymap->OutlierFrac[0] ? Ymap->OutlierFrac[0][j] : 0.0;
		  if(CmapChimQuality >= 3){
		    pmap->FragSd[0][i] = Ymap->FragSd[0] ? Ymap->FragSd[0][j] : 0.0;
		    pmap->ExpSd[0][i] = Ymap->ExpSd[0] ? Ymap->ExpSd[0][j] : 0.0;
		    pmap->FragCov[0][i] = Ymap->FragCov[0] ? Ymap->FragCov[0][j] : 0.0;
		    pmap->FragChiSq[0][i] = Ymap->FragChiSq[0] ? Ymap->FragChiSq[0][j] : 1.0;
		  }
		}
	      }
	      if(Ymap->SNRcnt[0]){
		pmap->SNRcnt[0][i] = Ymap->SNRcnt[0][j];
		pmap->SNRgmean[0][i] = Ymap->SNRgmean[0][j];
		pmap->lnSNRsd[0][i] = Ymap->lnSNRsd[0][j];
		pmap->SNRdist[0][i] = NULL;
		if(pmap->SNRcnt[0][i]){
		  pmap->SNRdist[0][i] = new double[pmap->SNRcnt[0][i]];
		  for(int t = 0; t < pmap->SNRcnt[0][i]; t++)
		    pmap->SNRdist[0][i][t] = Ymap->SNRdist[0][j][t];
		}
	      } else if(pmap->SNRcnt[0]){
		pmap->SNRcnt[0][i] = 0;
		pmap->SNRgmean[0][i] = 0.0;
		pmap->lnSNRsd[0][i] = 0.0;
		pmap->SNRdist[0][i] = NULL;
	      }
	      if(VERB>=3){
		printf("Copying site Y[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
		       j,N, Y[j],Ymap->id, i, pmap->numsite[0], Z[i], pmap->sitecov[0][i], pmap->sitecnt[0][i]);
		fflush(stdout);
	      }
	    }
	    if(DEBUG) assert(i-1 == NumSites+1);
	    if(Ymap->Mask[0])
	      pmap->Mask[0][NumSites + 1] |= Ymap->Mask[0][N + 1];
	  }
	}
	pmap->len = pmap->site[0][pmap->numsite[0]+1];
	if(VERB){
	  if(PairMergeRepeat)
	    printf("%d:",pairmergeIter);
	  printf("Created merged map (id=%lld,len=%0.3f) based on alignment between Y mapid=%d(id=%lld,len=%0.3f) and X mapid=%d(id=%lld,len=%0.3f):mapid=%d,id=%lld,len=%0.3f,sites=%d\n",
		 pmap->id, pmap->len, Ymap->mapid,Ymap->id, Ymap->site[0][Ymap->numsite[0]+1], 
		 Xmap->mapid,Xmap->id, Xmap->site[0][Xmap->numsite[0]+1],pmap->mapid,pmap->id,pmap->len,pmap->numsite[0]);
	  //	  printf(" Map Merge: Yid= %lld Xid= %lld or= %d : id= %lld\n",Ymap->id,Xmap->id,p->orientation, pmap->id);
#ifdef ASAN
	  fflush(stdout);
#endif
	}
      }
    }
  } else { // colors==2
    //    printf("colors=%d\n",colors);fflush(stdout);    exit(1);

    if(CmapChimQuality >= 1){
      printf("CMAP ChimQuality not yet supported for 2 color maps\n");
      fflush(stdout);exit(1);
    }
    if(TrimmedNoMerge){
      printf("-TrimmedNoMerge not yet supported for 2 color maps\n");
      fflush(stdout);exit(1);
    }
    for(size_t i=0; i < numaligns; i++){
      Calign *p = &alignment[i][1];

      Cmap *Ymap = YYmap[p[-1].mapid1];
      Cmap *Xmap = XXmap[p[-1].mapid2];
      int *NN = Ymap->numsite;
      int *MM = Xmap->numsite;
      FLOAT **YY = Ymap->site;
      FLOAT **XX = Xmap->site;

      int KK[2], LI[2],LJ[2],RI[2],RJ[2];
      for(int c = 0; c < colors; c++){
	KK[c] = p[c].numpairs;
	LI[c] = p[c].sites1[0];
	LJ[c] = p[c].sites2[0];
	RI[c] = p[c].sites1[KK[c]-1];
	RJ[c] = p[c].sites2[KK[c]-1];
      }

      /* average center to center distance over all aligned sites */
      double distanceYtoX = 0;
      double centerY = YY[0][NN[0]+1]*0.5;
      double centerX = XX[0][MM[0]+1]*0.5;
      for(int c = 0; c < colors; c++){
	for(int index = KK[c]; --index >= 0;){
	  int I = p[c].sites1[index];
	  int J = p[c].sites2[index];
	  distanceYtoX += (YY[c][I] - centerY) + (centerX - (p[c].orientation ? (XX[c][MM[c]+1] - XX[c][MM[c]+1-J]) : XX[c][J]));
	}
      }
      distanceYtoX /= max(1,KK[0]+KK[1]);

      double overlap[2], leftEnd[2],rightEnd[2], maxOutlier = 0.0; int maxOutlierM = 0;// HERE HERE (2-color)
      for(int c = 0; c < colors; c++){
	overlap[c] = (YY[c][RI[c]] - YY[c][LI[c]] + (p[c].orientation ? XX[c][MM[c]+1-LJ[c]] - XX[c][MM[c]+1-RJ[c]] : XX[c][RJ[c]] - XX[c][LJ[c]]))*0.5;
	leftEnd[c] = min(YY[c][LI[c]], (p[c].orientation ? XX[c][MM[c]+1] - XX[c][MM[c]+1-LJ[c]] : XX[c][LJ[c]]));
	rightEnd[c] = min(YY[c][NN[c]+1]-YY[c][RI[c]], (p->orientation ? XX[c][MM[c]+1-RJ[c]] : XX[c][MM[c]+1] - XX[c][RJ[c]]));
      }
      if(PairMergeMaxOutlier >= 0.0){/* locate largest outlier */
	for(int c = 0; c < colors; c++){
	  int lastI = LI[c];
	  int lastJ = LJ[c];
	  int I,J;
	  for(int k = 1; k < KK[c]; lastI = I, lastJ = J, k++){
	    I = p[c].sites1[k];
	    J = p[c].sites2[k];
	    if(DEBUG/* HERE HERE >=2 */) assert(J >= lastJ);
	    if(DEBUG/* HERE HERE >=2 */) assert(I >= lastI);
	  
	    /* outlier lies between sites lastI .. I on reference and lastJ .. J on query */
	    int qrystartidx = p->orientation ? MM[c] + 1 - lastJ : lastJ;
	    int qryendidx   = p->orientation ? MM[c] + 1 - J : J;
	    double qrystartpos = XX[c][qrystartidx]*1000.0*Xscale; /* QryStartPos */
	    double qryendpos   = XX[c][qryendidx  ]*1000.0*Xscale; /* QryEndPos */
	    double refstartpos = YY[c][lastI]*1000.0*Yscale; /* RefStartPos */
	    double refendpos   = YY[c][I]*1000.0*Yscale; /* RefEndPos */
	    if(DEBUG/* HERE HERE >=2 */) assert(refendpos >= refstartpos);
	    if(DEBUG) assert(qryendpos >= qrystartpos);// If this fails, fix next line
	    double Outlier = fabs(refendpos - refstartpos - fabs(qryendpos - qrystartpos));// indel size
	    int OutlierM = J-lastJ + I-lastI - 2;// misaligned labels

	    int outlier = (p[c].iscore[k] > p[c].outscore[k] + (RFLOAT)0.01 || (outlierExtend && (lastJ >= J || lastI >= I))) ? 1 : 0;
	    if(!outlier)/* not an outlier */
	      continue;
	    
#if 0 // Update for 2 color
	    if(PairMergeOutlierChiSq > 0.0 && CmapChimQuality >= 3){// skip if outlier variance is too large or has too many molecule outliers
	      if(refendpos - refstartpos < fabs(qryendpos - qrystartpos)){
		float *OutlierFrac = Ymap->OutlierFrac[c];
		float *FragSd = Ymap->FragSd[c];
		float *ExpSd = Ymap->ExpSd[c];
		float *FragCov = Ymap->FragCov[c];
		float *FragChiSq = Ymap->FragChiSq[c];
		if(DEBUG) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);
		int i = lastI;
		for(; i < I; i++){
		  double Yrange = (Y[i+1] - Y[i]) * Yscale;		  
		  if((FragChiSq[i] < PairMergeOutlierChiSq || OutlierFrac[i] > PairMergeOutlierFrac)
		     Outlier < PairMergeOutlierMaxSize &&
		     FragCov[i] < PairMergeOutlierMaxCov && 
		     FragSd[i] > PairMergeOutlierMinSf + PairMergeOutlierMinSr * Yrange &&
		     Yrange > PairMergeOutlierMinInterval)
		    break;
		}
		if(i < I)
		  continue;// outlier will be ignored
	      } else {
		float *OutlierFrac = Xmap->OutlierFrac[c];
		float *FragSd = Xmap->FragSd[c];
		float *FragCov = Xmap->FragCov[c];
		float *FragChiSq = Xmap->FragChiSq[c];
		if(DEBUG) assert(OutlierFrac != NULL && FragSd != NULL && FragCov != NULL && FragChiSq != NULL);
		int jmin = p->orientation ? M+1-J : lastJ;
		int jmax = p->orientation ? M+1-lastJ : J;
		int j = jmin;
		for(; j < jmax; j++){
		  double Xrange = (X[j+1] - X[j]) * Xscale;		  
		  if((FragChiSq[j] < PairMergeOutlierChiSq || OutlierFrac[j] > PairMergeOutlierFrac) &&
		     Outlier < PairMergeOutlierMaxSize &&
		     FragCov[j] < PairMergeOutlierMaxCov && 
		     FragSd[j] > PairMergeOutlierMinSf + PairMergeOutlierMinSr * Xrange &&
		     Xrange > PairMergeOutlierMinInterval)
		    break;
		}
		if(j < jmax)
		  continue;// outlier will be ignored
	      }
	    }

	    if(I - lastI + J - lastJ - 2 >= PairMergeMaxOutlierBalanced){
	      double YLeft = (Y[lastI+1]-Y[lastI]) * Yscale;
	      double YRight = (Y[I] - Y[I-1]) * Yscale;
	      double XLeft = (p->orientation ? X[M+1-lastJ] - X[M+1-(lastJ+1)] : X[lastJ+1] - X[lastJ]) * Xscale;
	      double XRight = (p->orientation ? X[M+1-(J-1)] - X[M+1-J] : X[J] - X[J-1]) * Xscale;

	      if(I > lastI + 2){
		double Middle = (Y[I-1] - Y[lastI+1]) * Yscale;
		double MinOutlier = Middle + max(0.0, max(YLeft-XLeft,YRight-XRight));
		Outlier = max(Outlier, MinOutlier);
	      }
	      if(J > lastJ + 2){
		double Middle = (p->orientation ? X[M+1-(lastJ+1)] - X[M+1-(J-1)] : X[J-1] - X[lastJ+1]) * Xscale;
		double MinOutlier = Middle + max(0.0, max(XLeft-YLeft,XRight-YRight));
		Outlier = max(Outlier, MinOutlier);
	      }
	    }

#endif
	    maxOutlier = max(maxOutlier, Outlier);
	    maxOutlierM = max(maxOutlierM, OutlierM);
	  }
	}
      }

      if(VERB){
	if(PairMergeRepeat)
	  printf("%d:",pairmergeIter);
	printf("Alignment %llu/%llu:Ymap->id=%lld(paired=%d),Xmap->id=%lld(paired=%d):score=%0.6f,logPV=%0.2f,or=%d,pairs=%d+%d,Y[LI=%d,%d]=%0.3f,%0.3f,Y[RI=%d,%d]=%0.3f,%0.3f,Y[N=%d,%d]=%0.3f,%0.3f\n",
	       (unsigned long long)i,(unsigned long long)numaligns, Ymap->id, Ymap->paired, Xmap->id,Xmap->paired,p[-1].score,p[-1].logPV,p[-1].orientation, p[0].numpairs, p[1].numpairs,
	       LI[0],LI[1],YY[0][LI[0]],YY[1][LI[1]],RI[0],RI[1],YY[0][RI[0]],YY[1][RI[1]],
	       Ymap->numsite[0], Ymap->numsite[1], YY[0][Ymap->numsite[0]], YY[1][Ymap->numsite[1]]);
	printf("\t X[LI=%d,%d]=%0.3f,%0.3f,X[RI=%d,%d]=%0.3f,%0.3f,X[M=%d,%d]=%0.3f,%0.3f,PairMerge=%0.1f,%0.4f,%0.1f,%0.1f,overlap=%0.3f,%0.3f,leftEnd=%0.3f,%0.3f,rightEnd=%0.3f,%0.3f,maxOutlier=%0.3f\n",
	       LJ[0], LJ[1], p[0].orientation ? (XX[0][MM[0]+1] - XX[0][MM[0]+1 - LJ[0]]) : XX[0][LJ[0]],
	       p[1].orientation ? (XX[1][MM[1]+1] - XX[1][MM[1]+1 - LJ[1]]) : XX[1][LJ[1]], 
	       RJ[0], RJ[1], p[0].orientation ? (XX[0][MM[0]+1] - XX[0][MM[0]+1 - RJ[0]]) : XX[0][RJ[0]], 
	       p[1].orientation ? (XX[1][MM[1]+1] - XX[1][MM[1]+1 - RJ[1]]) : XX[1][RJ[1]], 
	       Xmap->numsite[0], Xmap->numsite[1], XX[0][Xmap->numsite[0]], XX[1][Xmap->numsite[1]], 
	       PairMerge, PairMergeMaxEnd, PairMergeMaxOutlier, PairMergeMaxEndKB, overlap[0], overlap[1], leftEnd[0], leftEnd[1], rightEnd[0], rightEnd[1], maxOutlier);
	if(VERB>=3)
	  printf("\t Y[N+1]=%0.3f, X[M+1]=%0.3f\n",YY[0][Ymap->numsite[0]+1], XX[0][Xmap->numsite[0]+1]);
#ifdef ASAN
	fflush(stdout);
#endif
      }

      if(!Xmap->paired && !Ymap->paired && (overlap[0]+overlap[1])*0.5 > PairMerge && 
	 max(leftEnd[0]+leftEnd[1],rightEnd[0]+rightEnd[1]) <= max((overlap[0]+overlap[1]) * PairMergeMaxEnd,PairMergeMaxEndKB)
	 && (PairMergeMaxOutlier < 0.0 || !(maxOutlier > PairMergeMaxOutlier /* HERE HERE || maxOutlierM > PairMergeMaxOutlierM*/))){
	/* create a new map based on the merger of X and Y */

	// convert each 2-color map into an interleaved format : array of (c,index), corresponding to Ymap->site[c][index] and in ascending order of this value
	Cinterleaved *Ysite = new Cinterleaved[NN[0]+NN[1]+2];
	int *Yindex[2];/* Yindex[c][I=1..NN[c]] is the index into Ysite[] corresponding to Ymap->site[c][I] */
	Yindex[0] = new int[NN[0]+NN[1]+2];
	Yindex[1] = &Yindex[0][NN[0]+1];

	Ysite[0].c = 0;
	Ysite[0].index = 0;
	int Ysites = 1;
	int I[2] = {1,1};
	for(;I[0] <= NN[0] || I[1] <= NN[1];){
	  if(I[0] > NN[0] || (I[1] <= NN[1] && YY[0][I[0]] > YY[1][I[1]])){/* next site is YY[1][I[1]] */
	    Yindex[1][I[1]] = Ysites;
	    Ysite[Ysites].c = 1;
	    Ysite[Ysites++].index = I[1]++;
	  } else {/* next site is YY[0][I[0]] */
	    if(DEBUG) assert(I[1] > NN[1] || (I[0] <= NN[0] && YY[1][I[1]] >= YY[0][I[0]]));
	    Yindex[0][I[0]] = Ysites;
	    Ysite[Ysites].c = 0;
	    Ysite[Ysites++].index = I[0]++;
	  }
	}
	Ysite[Ysites].c = 0;
	Ysite[Ysites--].index = NN[0]+1;

	Cinterleaved *Xsite = new Cinterleaved[MM[0]+MM[1]+2];
	int *Xindex[2];/* Xindex[c][I=1..MM[c]] is the index into Xsite[] corresponding to Xmap->site[c][I] */
	Xindex[0] = new int[MM[0]+MM[1]+2];
	Xindex[1] = &Xindex[0][MM[0]+1];

	Xsite[0].c = 0;
	Xsite[0].index = 0;
	int Xsites = 1;
	I[0] = I[1] = 1;
	for(;I[0] <= MM[0] || I[1] <= MM[1];){
	  if(I[0] > MM[0] || (I[1] <= MM[1] && XX[0][I[0]] > XX[1][I[1]])){/* next site is XX[1][I[1]] */
	    Xindex[1][I[1]] = Xsites;
	    Xsite[Xsites].c = 1;
	    Xsite[Xsites++].index = I[1]++;
	  } else {/* next site is XX[0][I[0]] */
	    if(DEBUG) assert(I[1] > MM[1] || (I[0] <= MM[0] && XX[1][I[1]] >= XX[0][I[0]]));
	    Xindex[0][I[0]] = Xsites;
	    Xsite[Xsites].c = 0;
	    Xsite[Xsites++].index = I[0]++;
	  }
	}
	Xsite[Xsites].c = 0;
	Xsite[Xsites--].index = MM[0]+1;

	// convert the 2-color alignment into an interleaved format : array of (c,index1,index2)
	Cinterleaved *alignIL = new Cinterleaved[KK[0]+KK[1]];
	int alignILcnt = 0;
	int K[2] = {0,0};
	for(; K[0] < KK[0] || K[1] < KK[1];){
	  if(K[0] >= KK[0]){/* last alignment K[1] */
	    alignIL[alignILcnt].c = 1;
	    alignIL[alignILcnt++].index = K[1]++;
	  } else if(K[1] >= KK[1]){/* last alignment K[0] */
	    alignIL[alignILcnt].c = 0;
	    alignIL[alignILcnt++].index = K[0]++;
	  } else {/* find which aligned pair is closer to left end based on average distance to left ends */
	    int i0 = p[0].sites1[K[0]];
	    int j0 = p[0].sites2[K[0]];
	    int i1 = p[1].sites1[K[1]];
	    int j1 = p[1].sites2[K[1]];
	    FLOAT y0 = YY[0][i0];
	    FLOAT y1 = YY[1][i1];
	    FLOAT x0 = p[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-j0] : XX[0][j0];
	    FLOAT x1 = p[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-j1] : XX[1][j1];
	    if(x0 + y0 > x1 + y1){
	      alignIL[alignILcnt].c = 1;
	      alignIL[alignILcnt++].index = K[1]++;
	    } else {
	      alignIL[alignILcnt].c = 0;
	      alignIL[alignILcnt++].index = K[0]++;
	    }
	  }
	}
	if(DEBUG) assert(alignILcnt == KK[0]+KK[1]);

	Xmap->paired = 1 + nummaps;
	Ymap->paired = 1 + nummaps;

	maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	/* Update XXmap[], YYmap[] in case Gmap[] was reallocated */
	if(DEBUG) assert(PairSplit <= 0.0);
	if(!Cfirst || CfirstPairs==2)
	  XXmap = YYmap = Gmap;
	else {
	  YYmap = &Gmap[0];
	  XXmap = &Gmap[(CfirstPairs==1) ? 0 : Cfirst];
	}

	Cmap *pmap = Gmap[nummaps];
	pmap->mapid = nummaps++;
	pmap->fileid = -1;
	pmap->paired = 0;
	pmap->origmap = 0;
	if(Ymap->id < Xmap->id){
	  pmap->id = Ymap->id;
	  pmap->origid = Ymap->origid;
	} else {
	  pmap->id = Xmap->id;
	  pmap->origid = Xmap->origid;
	}

	for(int c = 0; c < colors; c++){
	  pmap->Nickase[c] = Xmap->Nickase[c];
	  pmap->numsite[c] = 0;
	}
	pmap->startloc = pmap->endloc = pmap->len = 0.0;

	//	int I = p->sites1[K/2], J= p->sites2[K/2];

	if(distanceYtoX > centerX - centerY){/* left end of Y comes first */
	  if(distanceYtoX < centerY - centerX){ /* right end of Y comes last */
	    int kLWM = 0, kHWM = 0;/* If kHWM > kLWM : region of X to copy is given by alignIL[kLWM..kHWM] :  XX[cstart][xstart[cstart]] .. XX[cend][xend[cend]] */

	    int xstart[2], xend[2], ystart[2], yend[2], cstart= -1, cend = -1;
	    FLOAT Xstart,Ystart,Xend,Yend;
	    
	    if(PairMergePreferNGS){  /*  check if X has lower SD in some region of the alignment, if so copy X in that aligned region */
	      double *VARdiff = new double[KK[0]+KK[1]];

	      for(int k = 0; k < KK[0]+KK[1] - 1; k++){
		int u = alignIL[k+1].index;
		if(u <= 0){
		  VARdiff[k] = 0.0;
		  continue;
		}
		int c = alignIL[k+1].c;
		int i = p[c].sites1[u-1];
		int j = p[c].sites2[u-1];
		int I = p[c].sites1[u];
		int J = p[c].sites2[u];
		double varY = 0, varX = 0.0;
		for(int t = i; t < I; t++){
		  double SDy = Ymap->siteSD[c][t];
		  varY += SDy*SDy;
		}
		for(int t = j; t < J; t++){
		  double SDx = p->orientation ? Xmap->siteSD[c][MM[c] - t] : Xmap->siteSD[c][t];
		  varX += SDx*SDx;
		}
		VARdiff[k] = varY - varX;
	      }

	      /* Now locate interval of X (xstart .. xend) with greatest reduction in variance over Y (largest sum of VARdiff[k = 0 .. KK[0]+KK[1]-2])*/
	      double sum = 0.0;/* cumulative sum of VARdiff[0..k] */
	      double sumLWM = 0.0;/*  Low Water Mark of cumulative sum of VARdiff[0..kLWM-1] */
	      double maxdiff = 0.0;/* largest value of sum-sumLWM */
	      // int kLWM = 0;/* k value of Low Water Mark */
	      // int kHWM = 0;/* k value of subsequent High Water Mark */
	      for(int k = 0; k < KK[0]+KK[1]-1; k++){
		sum += VARdiff[k];
		if(sum < sumLWM){
		  sumLWM = sum;
		  kHWM = kLWM = k+1;
		}
		if(sum - sumLWM > maxdiff){
		  maxdiff = sum - sumLWM;
		  kHWM = k+1;
		}
	      }
	      if(kLWM < kHWM){/* compute ystart[2],xstart[2],xend[2],yend[2] */
		if(DEBUG) assert(maxdiff > 0.0);
		cstart = alignIL[kLWM].c;
		ystart[cstart] = p[cstart].sites1[alignIL[kLWM].index];
		xstart[cstart] = p[cstart].sites2[alignIL[kLWM].index];
		Ystart = YY[cstart][ystart[cstart]];
		Xstart = p->orientation ? -XX[cstart][MM[cstart]+1-xstart[cstart]] : XX[cstart][xstart[cstart]];

		/* compute nearest previous site for other color c == 1-cstart */
		int c = 1-cstart;
		if(DEBUG) ystart[c] = xstart[c] = -1;

		for(int I = 1; I <= NN[c]+1; I++)
		  if(YY[c][I] > Ystart){
		    ystart[c] = I-1;
		    break;
		  }
		if(DEBUG) assert(ystart[c] >= 0);

		for(int J = 1; J <= MM[c]+1; J++)
		  if((p->orientation ? -XX[c][MM[c]+1-J] : XX[c][J]) > Xstart){
		    xstart[c] = J-1;
		    break;
		  }
		if(DEBUG) assert(xstart[c] >= 0);

		cend = alignIL[kHWM].c;
		yend[cend] = p[cend].sites1[alignIL[kHWM].index];
		xend[cend] = p[cend].sites2[alignIL[kHWM].index];
		Yend = YY[cend][yend[cend]];
		Xend = p->orientation ? -XX[cend][MM[cend]+1-xend[cend]] : XX[cend][xend[cend]];

		/* compute nearest previous site for other color c == 1-cstart */
		c = 1-cend;
		if(DEBUG) yend[c] = xend[c] = -1;
		
		for(int I = 1; I <= NN[c]+1; I++)
		  if(YY[c][I] > Yend){
		    yend[c] = I-1;
		    break;
		  }
		if(DEBUG) assert(yend[c] >= 0);

		for(int J = 1; J <= MM[c]+1; J++)
		  if((p->orientation ? -XX[c][MM[c]+1-J] : XX[c][J]) > Xend){
		    xend[c] = J-1;
		    break;
		  }
		if(DEBUG) assert(xend[c] >= 0);
	      }

	      if(VERB/* HERE >=2 */ && kLWM < kHWM){
		printf("Retaining part of overlapped X from XX[%d][%d]=%0.3f to XX[%d][%d]=%0.3f (replacing YY[%d][%d]=%0.3f to YY[%d][%d]=%0.3f\n",
		       cstart, xstart[cstart], Xmap->site[cstart][p->orientation ? (MM[cstart]+1-xstart[cstart]) : xstart[cstart]],
		       cend, xend[cend], Xmap->site[cend][p->orientation ? (MM[cend]+1-xend[cend]) : xend[cend]],
		       cstart, ystart[cstart], Ymap->site[cstart][ystart[cstart]], cend, yend[cend], Ymap->site[cend][yend[cend]]);
		fflush(stdout);
	      }
	    }
	    if(kLWM < kHWM){/* copy Y from left end to Ystart == YY[cstart][ystart[cstart], then X from Xstart == XX[cstart][xstart[cstart]] to 
			       Xend == XX[cend][xend[cend]], then Y from Yend == YY[cend][yend[cend]] to right end */

	      FLOAT Yoffset = YY[cstart][ystart[cstart]];
	      FLOAT Xoffset = Yoffset + (p->orientation ? XX[cstart][MM[cstart]+1-xstart[cstart]] : -XX[cstart][xstart[cstart]]);
	      FLOAT Yoffset2 = Yoffset + (p->orientation ? XX[cstart][MM[cstart]+1-xstart[cstart]] - XX[cend][MM[cend]+1-xend[cend]] : XX[cend][xend[cend]] - XX[cstart][xstart[cstart]]) - YY[cend][yend[cend]];

	      /* compute number of sites and allocate memory */
	      int NumSites[2];
	      for(int c = 0; c < colors; c++){
		NumSites[c] = NN[c] - (yend[c] - ystart[c]) + (xend[c] - xstart[c]);
		if(pmap->site[c]){
		  if(!pmap->blockmem) delete [] pmap->site[c];
		  if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		  if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		  if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
		}
		if(pmap->SNRcnt[c]){
		  if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		  if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		  if(pmap->SNRdist[c]){
		    for(int i = 0; i <= pmap->numsite[0]+1; i++)
		      if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		    delete [] pmap->SNRdist[c];
		  }
		  delete [] pmap->SNRcnt[c];
		}
		pmap->numsite[c] = NumSites[c];
		pmap->site[c] = new FLOAT[NumSites[c]+2];
		pmap->siteSD[c] = new double[NumSites[c]+2];
		pmap->sitecov[c] = new float[NumSites[c]+2];
		pmap->sitecnt[c] = new float[NumSites[c]+2];
		pmap->SNRcnt[c] = 0;
		if(Ymap->SNRcnt[c] || Xmap->SNRcnt[c]){
		  pmap->SNRcnt[c] = new int[NumSites[c]+2];
		  pmap->SNRgmean[c] = new double[NumSites[c]+2];
		  pmap->lnSNRsd[c] = new double[NumSites[c]+2];
		  pmap->SNRdist[c] = new double*[NumSites[c]+2];
		}

		/* left end */
		pmap->site[c][0] = 0.0;
		pmap->siteSD[c][0] = 0.0;
		pmap->sitecov[c][0] = 1;
		pmap->sitecnt[c][0] = 1;
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][0] = 0;
		  pmap->SNRgmean[c][0] = 0.0;
		  pmap->lnSNRsd[c][0] = 0.0;
		  pmap->SNRdist[c][0] = NULL;
		}

		/* right end */
		pmap->siteSD[c][NumSites[c]+1] = 0.0;
		pmap->sitecov[c][NumSites[c]+1] = 1;
		pmap->sitecnt[c][NumSites[c]+1] = 1;
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][NumSites[c]+1] = 0;
		  pmap->SNRgmean[c][NumSites[c]+1] = 0.0;
		  pmap->lnSNRsd[c][NumSites[c]+1] = 0.0;
		  pmap->SNRdist[c][NumSites[c]+1] = NULL;
		}
		pmap->blockmem = 0;
		FLOAT *Z = pmap->site[c];

		/* copy left end of Ymap from site 0 to site ystart[c] */
		int i = 0;
		for(; i <= ystart[c]; i++){
		  Z[i] = YY[c][i];
		  if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[c][i] = Ymap->siteSD[c][i];
		  pmap->sitecov[c][i] = Ymap->sitecov[c][i];
		  pmap->sitecnt[c][i] = Ymap->sitecnt[c][i];
		  if(Ymap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = Ymap->SNRcnt[c][i];
		    pmap->SNRgmean[c][i] = Ymap->SNRgmean[c][i];
		    pmap->lnSNRsd[c][i] = Ymap->lnSNRsd[c][i];
		    pmap->SNRdist[c][i] = 0;
		    if(pmap->SNRcnt[c][i]){
		      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			pmap->SNRdist[c][i][t] = Ymap->SNRdist[c][i][t];
		    }
		  } else if(pmap->SNRcnt[c]){/* default to empty SNR distribution */
		    pmap->SNRcnt[c][i] = 0;
		    pmap->SNRgmean[c][i] = 0.0;
		    pmap->lnSNRsd[c][i] = 0.0;
		    pmap->SNRdist[c][i] = NULL;
		  }
		  if(VERB>=3){
		    printf("Copying site YY[%d][i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   c, i, NN[c], YY[c][i],Ymap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == ystart[c]);
		if(DEBUG && c==cstart) assert(fabs(Z[i-1] - Yoffset) < 1e-8);

		if(p->orientation){/* copy Xmap from MM[c]+1-xstart[c] down to MM[c]+1-xend[c] */
		  int kstart = MM[c]+1-xstart[c];
		  int kend = MM[c]+1-xend[c];
		  //		double Xoffset = Z[i-1] + X[kstart];
		  for(int j = kstart; --j >= kend; i++){
		    Z[i] = Xoffset - XX[c][j];
		    if(DEBUG) assert(Z[i] >= Z[i-1]);
		    pmap->siteSD[c][i-1] = Xmap->siteSD[c][j];/* see definition of siteSD[] */
		    pmap->sitecov[c][i] = Xmap->sitecov[c][j];
		    pmap->sitecnt[c][i] = Xmap->sitecnt[c][j];
		    if(Xmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][j];
		      pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][j];
		      pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][j];
		      pmap->SNRdist[c][i] = 0;
		      if(pmap->SNRcnt[c][i]){
			pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
			for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			  pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][j][t];
		      }
		    } else if(pmap->SNRcnt[c]){/* use default empty SNR distribution */
		      pmap->SNRcnt[c][i] = 0;
		      pmap->SNRgmean[c][i] = 0.0;
		      pmap->lnSNRsd[c][i] = 0.0;
		      pmap->SNRdist[c][i] = NULL;
		    }
		    if(VERB>=3){
		      printf("Copying site X[%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			     c,j, MM[c], XX[c][j],Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		      fflush(stdout);
		    }
		  }
		} else {/* copy Xmap from xstart to xend */
		  int kstart = xstart[c];
		  int kend = xend[c];
		  //		double Xoffset = Z[i-1] - X[kstart];
		  for(int j = kstart; ++j <= kend; i++){
		    Z[i] = Xoffset + XX[c][j];
		    if(DEBUG) assert(Z[i] >= Z[i-1]);
		    pmap->siteSD[c][i-1] = Xmap->siteSD[c][j-1];/* see definition of siteSD[] */
		    pmap->sitecov[c][i] = Xmap->sitecov[c][j];
		    pmap->sitecnt[c][i] = Xmap->sitecnt[c][j];
		    if(Xmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][j];
		      pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][j];
		      pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][j];
		      pmap->SNRdist[c][i] = NULL;
		      if(pmap->SNRcnt[c][i]){
			pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
			for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			  pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][j][t];
		      }
		    } else if(pmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = 0;
		      pmap->SNRgmean[c][i] = 0.0;
		      pmap->lnSNRsd[c][i] = 0.0;
		      pmap->SNRdist[c][i] = NULL;
		    }
		    if(VERB>=3){
		      printf("Copying site X[%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			     c, j, MM[c], XX[c][j],Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		      fflush(stdout);
		    }
		  }
		} // else p->orientation == 0

		/* copy Ymap from yend[c] + 1 to right end site N+1 */
		int kstart = yend[c];
		int kend = NN[c]+1;
		for(int j = kstart; ++j <= kend; i++){
		  Z[i] = Yoffset2 + YY[c][j];
		  if(DEBUG) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[c][i-1] = Ymap->siteSD[c][j-1];/* see definition of siteSD[] */
		  pmap->sitecov[c][i] = Ymap->sitecov[c][j];
		  pmap->sitecnt[c][i] = Ymap->sitecnt[c][j];
		  if(Ymap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = Ymap->SNRcnt[c][j];
		    pmap->SNRgmean[c][i] = Ymap->SNRgmean[c][j];
		    pmap->lnSNRsd[c][i] = Ymap->lnSNRsd[c][j];
		    pmap->SNRdist[c][i] = NULL;
		    if(pmap->SNRcnt[c][i]){
		      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			pmap->SNRdist[c][i][t] = Ymap->SNRdist[c][j][t];
		    }
		  } else if(pmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = 0;
		    pmap->SNRgmean[c][i] = 0.0;
		    pmap->lnSNRsd[c][i] = 0.0;
		    pmap->SNRdist[c][i] = NULL;		
		  }
		  if(VERB>=3){
		    printf("Copying site Y[c=%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   c, j, NN[c], YY[c][j],Ymap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == NumSites[c]+1);
	      }// c = 0..colors-1
	    } else { /* just copy Y */
	      if(PairMergeRepeat){ // new code : instead of copying Y, just restore original Y
		if(VERB){
		  printf("%d:Created merged map (id=%lld,len=%0.3f) based on alignment between Y mapid=%d(id=%lld,len=%0.3f) and X mapid=%d(id=%lld,len=%0.3f):mapid=%d,id=%lld,len=%0.3f,sites=%d (Just keeping Ymap=%p)\n",
			 pairmergeIter,Ymap->id, Ymap->site[0][Ymap->numsite[0]+1], Ymap->mapid,Ymap->id, Ymap->site[0][Ymap->numsite[0]+1], 
			 Xmap->mapid,Xmap->id, Xmap->site[0][Xmap->numsite[0]+1], Ymap->mapid,Ymap->id,Ymap->len,Ymap->numsite[0],Ymap);
		  //		printf("%d:Restoring Ymap(id=%lld,mapid=%d) discarding Xmap(id=%lld,mapid=%d)\n",pairmergeIter,Ymap->id,Ymap->mapid,Xmap->id,Xmap->mapid);
		  fflush(stdout);
		}
		Ymap->paired = 0;
		if(DEBUG) assert(Ymap->mapid == p->mapid1);
		Xmap->paired = 1 + Ymap->mapid;
		nummaps--;
		continue;/* next Gmap[i] */
	      }

	      /* compute number of sites and allocate memory : NOTE orig* fields are not changed and become invalid/undefined, flagged with centercloned=1 */
	      int NumSites[2];
	      for(int c = 0; c < colors; c++){
		NumSites[c] = NN[c];
		if(pmap->site[c]){
		  if(!pmap->blockmem) delete [] pmap->site[c];
		  if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		  if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		  if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
		}
		if(pmap->SNRcnt[c]){
		  if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		  if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		  if(pmap->SNRdist[c]){
		    for(int i = 0; i <= pmap->numsite[0]+1; i++)
		      if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		    delete [] pmap->SNRdist[c];
		  }
		  delete [] pmap->SNRcnt[c];
		}
		pmap->numsite[c] = NumSites[c];
		pmap->site[c] = new FLOAT[NumSites[c]+2];
		pmap->siteSD[c] = new double[NumSites[c]+2];
		pmap->sitecov[c] = new float[NumSites[c]+2];
		pmap->sitecnt[c] = new float[NumSites[c]+2];
		pmap->SNRcnt[c] = 0;
		if(Ymap->SNRcnt[c]){
		  pmap->SNRcnt[c] = new int[NumSites[c]+2];
		  pmap->SNRgmean[c] = new double[NumSites[c]+2];
		  pmap->lnSNRsd[c] = new double[NumSites[c]+2];
		  pmap->SNRdist[c] = new double*[NumSites[c]+2];
		}

		/* left end */
		pmap->site[c][0] = 0.0;
		pmap->siteSD[c][0] = 0.0;
		pmap->sitecov[c][0] = 1;
		pmap->sitecnt[c][0] = 1;
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][0] = 0;
		  pmap->SNRgmean[c][0] = 0.0;
		  pmap->lnSNRsd[c][0] = 0.0;
		  pmap->SNRdist[c][0] = NULL;
		}

		/* right end */
		pmap->siteSD[c][NumSites[c]+1] = 0.0;
		pmap->sitecov[c][NumSites[c]+1] = 1;
		pmap->sitecnt[c][NumSites[c]+1] = 1;
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][NumSites[c]+1] = 0;
		  pmap->SNRgmean[c][NumSites[c]+1] = 0.0;
		  pmap->lnSNRsd[c][NumSites[c]+1] = 0.0;
		  pmap->SNRdist[c][NumSites[c]+1] = NULL;
		}

		pmap->blockmem = 0;

		FLOAT *Z = pmap->site[c];

		/* copy Ymap from site 0 to site NN[c]+1 */
		int i = 0;
		for(; i <= NN[c]+1; i++){
		  Z[i] = YY[c][i];
		  if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[c][i] = Ymap->siteSD[c][i];
		  pmap->sitecov[c][i] = Ymap->sitecov[c][i];
		  pmap->sitecnt[c][i] = Ymap->sitecnt[c][i];
		  if(Ymap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = Ymap->SNRcnt[c][i];
		    pmap->SNRgmean[c][i] = Ymap->SNRgmean[c][i];
		    pmap->lnSNRsd[c][i] = Ymap->lnSNRsd[c][i];
		    pmap->SNRdist[c][i] = NULL;
		    if(pmap->SNRcnt[c][i]){
		      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			pmap->SNRdist[c][i][t] = Ymap->SNRdist[c][i][t];
		    }
		  } else if(pmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = 0;
		    pmap->SNRgmean[c][i] = 0.0;
		    pmap->lnSNRsd[c][i] = 0.0;
		    pmap->SNRdist[c][i] = NULL;
		  }
		  if(VERB>=3){
		    printf("Copying site Y[c=%d][i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   c, i, NN[c], YY[c][i],Ymap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == NumSites[c]+1);	  
	      }
	    }
	  } else { /* copy left end of Y and right end of X */
	    int kmin = (KK[0]+KK[1])/2;/* default transition is at alignIL[kmin] */
	    if(PairMergePreferNGS){/* location transition point that retains the most NGS sites (with StdDev==0) */
	      double *VARdiff = new double[KK[0]+KK[1]];
	      /* NOTE : siteSD[c][i] is the sizing error for interval site[c][i+1] - site[c][i] */
	      for(int k = 0; k < KK[0]+KK[1] - 1; k++){
		int u = alignIL[k+1].index;
		if(u <= 0){
		  VARdiff[k] = 0.0;
		  continue;
		}
		int c = alignIL[k+1].c;

		int i = p[c].sites1[u-1];
		int j = p[c].sites2[u-1];
		int I = p[c].sites1[u];
		int J = p[c].sites2[u];
		double varY = 0, varX = 0.0;

		for(int t = i; t < I; t++){
		  double SDy = Ymap->siteSD[c][t];
		  varY += SDy*SDy;
		}
		for(int t = j; t < J; t++){
		  double SDx = p->orientation ? Xmap->siteSD[c][MM[c]-t] : Xmap->siteSD[c][t];
		  varX += SDx*SDx;
		}
		VARdiff[k] = varY - varX;
	      }
	      double VARmin = 0.0, VARsum = 0.0;
	      int origkmin = kmin;
	      kmin = 0;

	      for(int k = 0; k < KK[0]+KK[1]-1; k++){
		VARsum += VARdiff[k];
		if(VARsum < VARmin){
		  VARmin = VARsum;
		  kmin = k+1;
		}
	      }
	      delete [] VARdiff;
	      if(VERB/* HERE >=2 */){
		int c1 = alignIL[origkmin].c;
		int u1 = alignIL[origkmin].index;
		int I1 = p[c1].sites1[u1];
		int J1 = p[c1].sites1[u1];
		int cmin = alignIL[kmin].c;
		int umin = alignIL[kmin].index;
		int Imin = p[cmin].sites1[umin];
		int Jmin = p[cmin].sites2[umin];

		printf("Changing crossover point from (YY[%d][%d]=%0.3f,X[%d]=%0.3f) to (YY[%d][%d]=%0.3f,X[%d]=%0.3f)\n",
		       c1, I1, Ymap->site[c1][I1], p->orientation ? (MM[c1]+1-J1) : J1, Xmap->site[c1][p->orientation ? (MM[c1]+1-J1) : J1], 
		       cmin, Imin, Ymap->site[cmin][Imin],
		       p->orientation ? (MM[cmin]+1-Jmin) : Jmin, Xmap->site[cmin][p->orientation ? (MM[cmin]+1-Jmin) : Jmin]);
		fflush(stdout);
	      }
	    }

	    int xmin[2],ymin[2];
	    int cmin = alignIL[kmin].c;
	    int umin = alignIL[kmin].index;
	    ymin[cmin] = p[cmin].sites1[umin];
	    xmin[cmin] = p[cmin].sites2[umin];
	    FLOAT Ymin = YY[cmin][ymin[cmin]];
	    FLOAT Xmin = p->orientation ? -XX[cmin][MM[cmin]+1-xmin[cmin]] : XX[cmin][xmin[cmin]];

	    /* compute nearest previous site for other colors c == 1-cmin */
	    int c = 1-cmin;
	    if(DEBUG) ymin[c] = xmin[c] = -1;	    

	    for(int I = 1; I <= NN[c]+1; I++)
	      if(YY[c][I] > Ymin){
		ymin[c] = I-1;
		break;
	      }
	    if(DEBUG) assert(ymin[c] >= 0);

	    for(int J = 1; J <= MM[c]+1; J++)
	      if((p->orientation ? -XX[c][MM[c]+1-J] : XX[c][J]) > Xmin){
		xmin[c] = J-1;
		break;
	      }
	    if(DEBUG) assert(xmin[c] >= 0);

	    /* copy Y from left end to YY[cmin][ymin[cmin]], then X from XX[cmin][xmin[cmin]] to right end */
	    FLOAT Xstart = Ymin + (p->orientation ? XX[cmin][MM[cmin]+1-xmin[cmin]] : -XX[cmin][xmin[cmin]]);

	    /* compute number of sites and allocate memory */
	    int NumSites[2];
	    for(int c = 0; c < colors; c++){
	      NumSites[c] = ymin[c] + (p->orientation ? (MM[c] + 1 - xmin[c]) - 1 : MM[c] - xmin[c]);
	      if(pmap->site[c]){
		if(!pmap->blockmem) delete [] pmap->site[c];
		if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
	      }
	      if(pmap->SNRcnt[c]){
		if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		if(pmap->SNRdist[c]){
		  for(int i = 0; i <= pmap->numsite[0]+1; i++)
		    if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		  delete [] pmap->SNRdist[c];
		}
		delete [] pmap->SNRcnt[c];
	      }
	      pmap->numsite[c] = NumSites[c];
	      pmap->site[c] = new FLOAT[NumSites[c]+2];
	      pmap->siteSD[c] = new double[NumSites[c]+2];
	      pmap->sitecov[c] = new float[NumSites[c]+2];
	      pmap->sitecnt[c] = new float[NumSites[c]+2];
	      pmap->SNRcnt[c] = 0;
	      if(Ymap->SNRcnt[c] || Xmap->SNRcnt[c]){
		pmap->SNRcnt[c] = new int[NumSites[c]+2];
		pmap->SNRgmean[c] = new double[NumSites[c]+2];
		pmap->lnSNRsd[c] = new double[NumSites[c]+2];
		pmap->SNRdist[c] = new double*[NumSites[c]+2];
	      }

	      /* left end */
	      pmap->site[c][0] = 0.0;
	      pmap->siteSD[c][0] = 0.0;
	      pmap->sitecov[c][0] = 1;
	      pmap->sitecnt[c][0] = 1;
	      if(pmap->SNRcnt[c]){
		pmap->SNRcnt[c][0] = 0;
		pmap->SNRgmean[c][0] = 0.0;
		pmap->lnSNRsd[c][0] = 0.0;
		pmap->SNRdist[c][0] = NULL;
	      }

	      /* right end */
	      pmap->siteSD[c][NumSites[c]+1] = 0.0;
	      pmap->sitecov[c][NumSites[c]+1] = 1;
	      pmap->sitecnt[c][NumSites[c]+1] = 1;
	      if(pmap->SNRcnt[c]){
		pmap->SNRcnt[c][NumSites[c]+1] = 0;
		pmap->SNRgmean[c][NumSites[c]+1] = 0.0;
		pmap->lnSNRsd[c][NumSites[c]+1] = 0.0;
		pmap->SNRdist[c][NumSites[c]+1] = NULL;
	      }

	      pmap->blockmem = 0;
	      FLOAT *Z = pmap->site[c];

	      /* copy left end of Ymap from site 0 to site ymin[c] */
	      int i = 0;
	      for(; i <= ymin[c]; i++){
		Z[i] = YY[c][i];
		if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		pmap->siteSD[c][i] = Ymap->siteSD[c][i];
		pmap->sitecov[c][i] = Ymap->sitecov[c][i];
		pmap->sitecnt[c][i] = Ymap->sitecnt[c][i];
		if(Ymap->SNRcnt[c]){
		  pmap->SNRcnt[c][i] = Ymap->SNRcnt[c][i];
		  pmap->SNRgmean[c][i] = Ymap->SNRgmean[c][i];
		  pmap->lnSNRsd[c][i] = Ymap->lnSNRsd[c][i];
		  pmap->SNRdist[c][i] = NULL;
		  if(pmap->SNRcnt[c][i]){
		  pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		  for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
		    pmap->SNRdist[c][i][t] = Ymap->SNRdist[c][i][t];
		  }
		} else if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][i] = 0;
		  pmap->SNRgmean[c][i] = 0.0;
		  pmap->lnSNRsd[c][i] = 0.0;
		  pmap->SNRdist[c][i] = NULL;
		}
		if(VERB>=3){
		  printf("Copying site Y[c=%d][i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			 c, i, NN[c], YY[c][i], Ymap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		  fflush(stdout);
		}
	      }

	      if(p->orientation){/* copy right end of Xmap from MM[c]+1-xmin[c]-1 downto 0 */
		int k = MM[c] + 1 - xmin[c];
		// double Xstart = Z[i-1] + X[k];
		//		double Xstart = Ymin + XX[cmin][MM[cmin]+1-xmin[cmin]];
		for(int j = k; --j >= 0; i++){
		  Z[i] = Xstart - XX[c][j];
		  if(DEBUG) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[c][i-1] = Xmap->siteSD[c][j];/* see definition of siteSD[] */
		  pmap->sitecov[c][i] = Xmap->sitecov[c][j];
		  pmap->sitecnt[c][i] = Xmap->sitecnt[c][j];
		  if(Xmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][j];
		    pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][j];
		    pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][j];
		    pmap->SNRdist[c][i] = NULL;
		    if(pmap->SNRcnt[c][i]){
		      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][j][t];
		    }
		  } else if(pmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = 0;
		    pmap->SNRgmean[c][i] = 0.0;
		    pmap->lnSNRsd[c][i] = 0.0;
		    pmap->SNRdist[c][i] = NULL;
		  }
		  if(VERB>=3){
		    printf("Copying site X[c=%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   c, j, MM[c], XX[c][j],Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == NumSites[c]+1);
	      } else {/* copy right end of Xmap from  xmin[c]+1 to MM[c]+1 */
		int k = xmin[c];
		// double Xstart = Z[i-1] - X[k];
		// double Xstart = Ymin - XX[cmin][xmin[cmin]];
		for(int j = k; ++j <= MM[c] + 1; i++){
		  Z[i] = Xstart + XX[c][j];
		  if(DEBUG) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[c][i-1] = Xmap->siteSD[c][j-1];/* see definition of siteSD[] */
		  pmap->sitecov[c][i] = Xmap->sitecov[c][j];
		  pmap->sitecnt[c][i] = Xmap->sitecnt[c][j];
		  if(Xmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][j];
		    pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][j];
		    pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][j];
		    pmap->SNRdist[c][i] = NULL;
		    if(pmap->SNRcnt[c][i]){
		      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][j][t];
		    }
		  } else if(pmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = 0;
		    pmap->SNRgmean[c][i] = 0.0;
		    pmap->lnSNRsd[c][i] = 0.0;
		    pmap->SNRdist[c][i] = NULL;
		  }
		  if(VERB >=3){
		    printf("Copying site X[c=%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   c, j, MM[c], XX[c][j], Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == NumSites[c]+1);
	      }// else p->alignment==0
	    }// c = 0..colors-1
	  } /* copy left end of Y and right end of X */
	} else { /* left end of X comes first */
	  if(distanceYtoX > centerY - centerX){/* right end of X comes Last */
	    int kLWM = 0, kHWM = 0;/* If kHWM > kLWM : region of Y to copy is given by alignIL[kLWM..kHWM] : YY[cstart][ystart[cstart]] .. YY[cend][yend[cend]] */

	    int xstart[2], xend[2], ystart[2], yend[2], cstart = -1, cend = -1;
	    FLOAT Xstart,Ystart,Xend,Yend;

	    if(PairMergePreferNGS){ /* check if Y has lower SD in some region of the alignment, if so copy Y in that aligned region */
	      double *VARdiff = new double[KK[0]+KK[1]];

	      for(int k = 0; k < KK[0]+KK[1] - 1; k++){
		int u = alignIL[k+1].index;
		if(u <= 0){
		  VARdiff[k] = 0.0;
		  continue;
		}
		int c = alignIL[k+1].c;

		int i = p[c].sites1[u-1];
		int j = p[c].sites2[u-1];
		int I = p[c].sites1[u];
		int J = p[c].sites2[u];
		double varY = 0, varX = 0.0;
		for(int t = i; t < I; t++){
		  double SDy = Ymap->siteSD[c][t];
		  varY += SDy*SDy;
		}
		for(int t = j; t < J; t++){
		  double SDx = p->orientation ? Xmap->siteSD[c][MM[c]-t] : Xmap->siteSD[c][t];
		  varX += SDx*SDx;
		}
		VARdiff[k] = varX - varY;
	      }

	      /* Now locate interval of Y (ystart .. yend) with greatest reduction in variance over X (largest sum of Vardiff[k = 0 .. KK[0]+KK[1]-2])*/
	      double sum = 0.0;/* cumulative sum of VARdiff[0..k] */
	      double sumLWM = 0.0;/*  Low Water Mark of cumulative sum of VARdiff[0..kLWM-1] */
	      double maxdiff = 0.0;/* largest value of sum-sumLWM */
	      for(int k = 0; k < KK[0]+KK[1]-1; k++){
		sum += VARdiff[k];
		if(sum < sumLWM){
		  sumLWM = sum;
		  kHWM = kLWM = k+1;
		}
		if(sum - sumLWM > maxdiff){
		  maxdiff = sum - sumLWM;
		  kHWM = k+1;
		}
	      }
	      delete [] VARdiff;

	      if(kLWM < kHWM){/* compute ystart[2],xstart[2],xend[2],yend[2] */
		if(DEBUG) assert(maxdiff > 0.0);
		cstart = alignIL[kLWM].c;
		ystart[cstart] = p[cstart].sites1[alignIL[kLWM].index];
		xstart[cstart] = p[cstart].sites2[alignIL[kLWM].index];
		Ystart = YY[cstart][ystart[cstart]];
		Xstart = p->orientation ? -XX[cstart][MM[cstart]+1-xstart[cstart]] : XX[cstart][xstart[cstart]];

		/* compute nearest previous site for other color c == 1-cstart */
		int c = 1-cstart;
		if(DEBUG) ystart[c] = xstart[c] = -1;

		for(int I = 1; I <= NN[c]+1; I++)
		  if(YY[c][I] > Ystart){
		    ystart[c] = I-1;
		    break;
		  }
		if(DEBUG) assert(ystart[c] >= 0);

		for(int J = 1; J <= MM[c]+1; J++)
		  if((p->orientation ? -XX[c][MM[c]+1-J] : XX[c][J]) > Xstart){
		    xstart[c] = J-1;
		    break;
		  }
		if(DEBUG) assert(xstart[c] >= 0);

		cend = alignIL[kHWM].c;
		yend[cend] = p[cend].sites1[alignIL[kHWM].index];
		xend[cend] = p[cend].sites2[alignIL[kHWM].index];
		Yend = YY[cend][yend[cend]];
		Xend = p->orientation ? -XX[cend][MM[cend]+1-xend[cend]] : XX[cend][xend[cend]];

		/* compute nearest previous site for other color c == 1-cstart */
		c = 1-cend;
		if(DEBUG) yend[c] = xend[c] = -1;
		
		for(int I = 1; I <= NN[c]+1; I++)
		  if(YY[c][I] > Yend){
		    yend[c] = I-1;
		    break;
		  }
		if(DEBUG) assert(yend[c] >= 0);

		for(int J = 1; J <= MM[c]+1; J++)
		  if((p->orientation ? -XX[c][MM[c]+1-J] : XX[c][J]) > Xend){
		    xend[c] = J-1;
		    break;
		  }
		if(DEBUG) assert(xend[c] >= 0);
	      }

	      if(VERB/* HERE >=2 */ && kLWM < kHWM){
		printf("Retaining part of overlapped Y from YY[%d][%d]=%0.3f to YY[%d][%d]=%0.3f (replacing XX[%d][%d]=%0.3f to XX[%d][%d]=%0.3f\n",
		       cstart, ystart[cstart], Ymap->site[cstart][ystart[cstart]], cend, yend[cend], Ymap->site[cend][yend[cend]],
		       cstart, xstart[cstart], Xmap->site[cstart][p->orientation ? (MM[cstart]+1-xstart[cstart]) : xstart[cstart]],
		       cend, xend[cend], Xmap->site[cend][p->orientation ? (MM[cend]+1-xend[cend]) : xend[cend]]);
		fflush(stdout);
	      }
	    }
	    if(kLWM < kHWM){/* copy X from left end to XX[cstart][xstart[cstart]] , then Y from YY[cstart][ystart[cstart]] to
			       YY[cend][yend[cend]], then X from XX[cend][xend[cend]] to right end */
	      FLOAT Zoffset = p->orientation ? XX[cstart][MM[cstart]+1] - XX[cstart][MM[cstart]+1-xstart[cstart]] : XX[cstart][xstart[cstart]];
	      FLOAT Yoffset = Zoffset - YY[cstart][ystart[cstart]];
	      FLOAT Xoffset2 = Yoffset + YY[cend][yend[cend]] + (p->orientation ? XX[cend][MM[cend]+1-xend[cend]] : -XX[cend][xend[cend]]);

	      /* compute number of sites and allocate memory */
	      int NumSites[2];
	      for(int c = 0; c < colors; c++){
		NumSites[c] = MM[c] - (xend[c] - xstart[c]) + (yend[c] - ystart[c]);
		if(pmap->site[c]){
		  if(!pmap->blockmem) delete [] pmap->site[c];
		  if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		  if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		  if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
		}
		if(pmap->SNRcnt[c]){
		  if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		  if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		  if(pmap->SNRdist[c]){
		    for(int i = 0; i <= pmap->numsite[0]+1; i++)
		      if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		    delete [] pmap->SNRdist[c];
		  }
		  delete [] pmap->SNRcnt[c];
		}
		pmap->numsite[c] = NumSites[c];
		pmap->site[c] = new FLOAT[NumSites[c]+2];
		pmap->siteSD[c] = new double[NumSites[c]+2];
		pmap->sitecov[c] = new float[NumSites[c]+2];
		pmap->sitecnt[c] = new float[NumSites[c]+2];
		pmap->SNRcnt[c] = 0;
		if(Ymap->SNRcnt[c] || Xmap->SNRcnt[c]){
		  pmap->SNRcnt[c] = new int[NumSites[c]+2];
		  pmap->SNRgmean[c] = new double[NumSites[c]+2];
		  pmap->lnSNRsd[c] = new double[NumSites[c]+2];
		  pmap->SNRdist[c] = new double*[NumSites[c]+2];
		}

		/* left end */
		pmap->site[c][0] = 0.0;
		pmap->siteSD[c][0] = 0.0;
		pmap->sitecov[c][0] = 1;
		pmap->sitecnt[c][0] = 1;
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][0] = 0;
		  pmap->SNRgmean[c][0] = 0.0;
		  pmap->lnSNRsd[c][0] = 0.0;
		  pmap->SNRdist[c][0] = NULL;
		}

		/* right end */
		pmap->siteSD[c][NumSites[c]+1] = 0.0;
		pmap->sitecov[c][NumSites[c]+1] = 1;
		pmap->sitecnt[c][NumSites[c]+1] = 1;
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][NumSites[c]+1] = 0;
		  pmap->SNRgmean[c][NumSites[c]+1] = 0.0;
		  pmap->lnSNRsd[c][NumSites[c]+1] = 0.0;
		  pmap->SNRdist[c][NumSites[c]+1] = NULL;
		}

		pmap->blockmem = 0;
		FLOAT *Z = pmap->site[c];
	    
		int i = 0;
		if(p->orientation){/* reversed X orientation */
		  /* copy left end of Xmap from site M+1 downto M+1-xstart */
		  FLOAT Xoffset = XX[cstart][MM[cstart]+1];
		  for(int j = MM[c]+1; j >= MM[c]+1-xstart[c]; j--, i++){
		    Z[i] = Xoffset - XX[c][j];
		    if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		    pmap->siteSD[c][i] = Xmap->siteSD[c][j-1]; /* see definition of siteSD[] */
		    pmap->sitecov[c][i] = Xmap->sitecov[c][j];
		    pmap->sitecnt[c][i] = Xmap->sitecnt[c][j];
		    if(Xmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][j];
		      pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][j];
		      pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][j];
		      pmap->SNRdist[c][i] = NULL;
		      if(pmap->SNRcnt[c][i]){
			pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
			for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			  pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][j][t];
		      }
		    } else if(pmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = 0;
		      pmap->SNRgmean[c][i] = 0.0;
		      pmap->lnSNRsd[c][i] = 0.0;
		      pmap->SNRdist[c][i] = NULL;
		    }
		    if(VERB>=3){
		      printf("Copying site XX[c=%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			     c, j, MM[c], XX[c][j], Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		      fflush(stdout);
		    }
		  }
		  if(DEBUG) assert(i-1 == (MM[c]+1)-(MM[c]+1-xstart[c]));
		  if(DEBUG && c==cstart) assert(fabs(Z[i-1] - Zoffset) < 1e-8);

		  /* copy Ymap from site ystart[c]+1 to yend[c] */
		  int j = ystart[c];
		  //		  double Yoffset = Z[i-1] - Y[j];
		  for(; ++j <= yend[c]; i++){
		    Z[i] = Yoffset + YY[c][j];
		    if(DEBUG) assert(Z[i] >= Z[i-1]);
		    pmap->siteSD[c][i-1] = Ymap->siteSD[c][j-1];/* see definition of siteSD[] */
		    pmap->sitecov[c][i] = Ymap->sitecov[c][j];
		    pmap->sitecnt[c][i] = Ymap->sitecnt[c][j];
		    if(Ymap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = Ymap->SNRcnt[c][j];
		      pmap->SNRgmean[c][i] = Ymap->SNRgmean[c][j];
		      pmap->lnSNRsd[c][i] = Ymap->lnSNRsd[c][j];
		      pmap->SNRdist[c][i] = NULL;
		      if(pmap->SNRcnt[c][i]){
			pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
			for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			  pmap->SNRdist[c][i][t] = Ymap->SNRdist[c][j][t];
		      }
		    } else if(pmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = 0;
		      pmap->SNRgmean[c][i] = 0.0;
		      pmap->lnSNRsd[c][i] = 0.0;
		      pmap->SNRdist[c][i] = NULL;
		    }
		    if(VERB>=3){
		      printf("Copying site Y[c=%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			     c, j, NN[c], YY[c][j], Ymap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		      fflush(stdout);
		    }
		  }
		  if(DEBUG) assert(i-1 == (MM[c]+1)-(MM[c]+1-xstart[c]) + (yend[c]-ystart[c]));

		  /* copy right end of Xmap from site MM[c]+1-xend[c]-1 downto 0 */
		  int k = MM[c]+1-xend[c];
		  //		  Xoffset = Z[i-1] + X[k];
		  for(int j = k; --j >= 0; i++){
		    Z[i] = Xoffset2 - XX[c][j];
		    if(DEBUG) assert(Z[i] >= Z[i-1]);
		    pmap->siteSD[c][i-1] = Xmap->siteSD[c][j];/* see definition of siteSD[] */
		    pmap->sitecov[c][i] = Xmap->sitecov[c][j];
		    pmap->sitecnt[c][i] = Xmap->sitecnt[c][j];
		    if(Xmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][j];
		      pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][j];
		      pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][j];
		      pmap->SNRdist[c][i] = NULL;
		      if(pmap->SNRcnt[c][i]){
			pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
			for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			  pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][j][t];
		      }
		    } else if(pmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = 0;
		      pmap->SNRgmean[c][i] = 0.0;
		      pmap->lnSNRsd[c][i] = 0.0;
		      pmap->SNRdist[c][i] = NULL;
		    }
		    if(VERB>=3){
		      printf("Copying site X[c=%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			     c, j, MM[c], XX[c][j], Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		      fflush(stdout);
		    }
		  }
		  if(DEBUG) assert(i-1 == NumSites[c]+1);
		} else {/* normal X orientation */
		  /* copy left end of Xmap from site 0 to xstart[c] */
		  for(; i <= xstart[c]; i++){
		    Z[i] = XX[c][i];
		    if(DEBUG) assert(Z[i] >= Z[i-1]);
		    pmap->siteSD[c][i] = Xmap->siteSD[c][i];
		    pmap->sitecov[c][i] = Xmap->sitecov[c][i];
		    pmap->sitecnt[c][i] = Xmap->sitecnt[c][i];
		    if(Xmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][i];
		      pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][i];
		      pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][i];
		      pmap->SNRdist[c][i] = NULL;
		      if(pmap->SNRcnt[c][i]){
			pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
			for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			  pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][i][t];
		      }
		    } else if(pmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = 0;
		      pmap->SNRgmean[c][i] = 0.0;
		      pmap->lnSNRsd[c][i] = 0.0;
		      pmap->SNRdist[c][i] = NULL;
		    }
		    if(VERB>=3){
		      printf("Copying site X[c=%d][i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			     c, i, MM[c], XX[c][i],Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		      fflush(stdout);
		    }
		  }
		  if(DEBUG) assert(i-1 == xstart[c]);
	      
		  /* copy Ymap from site ystart[c]+1 to yend[c] */
		  int j = ystart[c];
		  // double Yoffset = Z[i-1] - Y[j];
		  for(; ++j <= yend[c]; i++){
		    Z[i] = Yoffset + YY[c][j];
		    if(DEBUG) assert(Z[i] >= Z[i-1]);
		    pmap->siteSD[c][i-1] = Ymap->siteSD[c][j-1];/* see definition of siteSD[] */
		    pmap->sitecov[c][i] = Ymap->sitecov[c][j];
		    pmap->sitecnt[c][i] = Ymap->sitecnt[c][j];
		    if(Ymap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = Ymap->SNRcnt[c][j];
		      pmap->SNRgmean[c][i] = Ymap->SNRgmean[c][j];
		      pmap->lnSNRsd[c][i] = Ymap->lnSNRsd[c][j];
		      pmap->SNRdist[c][i] = NULL;
		      if(pmap->SNRcnt[c][i]){
			pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
			for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			  pmap->SNRdist[c][i][t] = Ymap->SNRdist[c][j][t];
		      }
		    } else if(pmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = 0;
		      pmap->SNRgmean[c][i] = 0.0;
		      pmap->lnSNRsd[c][i] = 0.0;
		      pmap->SNRdist[c][i] = NULL;
		    }
		    if(VERB>=3){
		      printf("Copying site Y[c=%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			     c, j, NN[c], YY[c][j], Ymap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		      fflush(stdout);
		    }
		  }
		  if(DEBUG) assert(i-1 == xstart[c] + (yend[c]-ystart[c]));
	      
		  /* copy left end of Xmap from site xend[c]+1 to MM[c]+1 */
		  int k = xend[c];
		  //		  double Xoffset = Z[i-1] - X[k];
		  for(int j = k; ++j <= MM[c]+1; i++){
		    Z[i] = Xoffset2 + XX[c][j];
		    if(DEBUG) assert(Z[i] >= Z[i-1]);
		    pmap->siteSD[c][i-1] = Xmap->siteSD[c][j-1];/* see definition of siteSD[] */
		    pmap->sitecov[c][i] = Xmap->sitecov[c][j];
		    pmap->sitecnt[c][i] = Xmap->sitecnt[c][j];
		    if(Xmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][j];
		      pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][j];
		      pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][j];
		      pmap->SNRdist[c][i] = NULL;
		      if(pmap->SNRcnt[c][i]){
			pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
			for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			  pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][j][t];
		      }
		    } else if(pmap->SNRcnt[c]){
		      pmap->SNRcnt[c][i] = 0;
		      pmap->SNRgmean[c][i] = 0.0;
		      pmap->lnSNRsd[c][i] = 0.0;
		      pmap->SNRdist[c][i] = NULL;
		    }
		    if(VERB >=3){
		      printf("Copying site X[c=%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			     c, j, MM[c], XX[c][j], Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		      fflush(stdout);
		    }
		  }
		  if(DEBUG) assert(i-1 == NumSites[c]+1);
		} /* normal X orientation */
	      }// c = 0 .. colors-1
	    } else {/* just copy X */
	      if(PairMergeRepeat){ // new code : instead of copying X, just restore inal X
		if(VERB){
		  printf("%d:Created merged map (id=%lld,len=%0.3f) based on alignment between Y mapid=%d(id=%lld,len=%0.3f) and X mapid=%d(id=%lld,len=%0.3f):mapid=%d,id=%lld,len=%0.3f,sites=%d (Just keeping Xmap=%p)\n",
			 pairmergeIter,Xmap->id, Xmap->site[0][Xmap->numsite[0]+1], Ymap->mapid,Ymap->id, Ymap->site[0][Ymap->numsite[0]+1], 
			 Xmap->mapid,Xmap->id, Xmap->site[0][Xmap->numsite[0]+1],Xmap->mapid,Xmap->id,Xmap->len,Xmap->numsite[0],Xmap);
		  /*		printf("%d:Restoring Xmap(id=%lld,mapid=%d) discarding Ymap(id=%lld,mapid=%d)\n",pairmergeIter,Xmap->id,Xmap->mapid,Ymap->id,Ymap->mapid);*/
		  fflush(stdout);
		}
		Xmap->paired = 0;
		if(DEBUG) assert(Xmap->mapid == p->mapid2);
		Ymap->paired = 1 + Xmap->mapid;
		nummaps--;
		continue;/* next Gmap[i] */
	      }

	      /* compute number of sites and allocate memory : NOTE Cmap->orig* fields are not changed and become invalid/undefined, flagged with centercloned=1 */
	      int NumSites[2];
	      for(int c = 0; c < colors; c++){
		NumSites[c] = MM[c];
		if(pmap->site[c]){
		  if(!pmap->blockmem) delete [] pmap->site[c];
		  if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		  if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		  if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
		}
		if(pmap->SNRcnt[c]){
		  if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		  if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		  if(pmap->SNRdist[c]){
		    for(int i = 0; i <= pmap->numsite[0]+1; i++)
		      if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		    delete [] pmap->SNRdist[c];
		  }
		  delete [] pmap->SNRcnt[c];
		}
		pmap->numsite[c] = NumSites[c];
		pmap->site[c] = new FLOAT[NumSites[c]+2];
		pmap->siteSD[c] = new double[NumSites[c]+2];
		pmap->sitecov[c] = new float[NumSites[c]+2];
		pmap->sitecnt[c] = new float[NumSites[c]+2];
		pmap->SNRcnt[c] = 0;
		if(Xmap->SNRcnt[c]){
		  pmap->SNRcnt[c] = new int[NumSites[c]+2];
		  pmap->SNRgmean[c] = new double[NumSites[c]+2];
		  pmap->lnSNRsd[c] = new double[NumSites[c]+2];
		  pmap->SNRdist[c] = new double*[NumSites[c]+2];
		}

		/* left end */
		pmap->site[c][0] = 0.0;
		pmap->siteSD[c][0] = 0.0;
		pmap->sitecov[c][0] = 1;
		pmap->sitecnt[c][0] = 1;
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][0] = 0;
		  pmap->SNRgmean[c][0] = 0.0;
		  pmap->lnSNRsd[c][0] = 0.0;
		  pmap->SNRdist[c][0] = NULL;
		}

		/* right end */
		pmap->siteSD[c][NumSites[c]+1] = 0.0;
		pmap->sitecov[c][NumSites[c]+1] = 1;
		pmap->sitecnt[c][NumSites[c]+1] = 1;
		if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][NumSites[c]+1] = 0;
		  pmap->SNRgmean[c][NumSites[c]+1] = 0.0;
		  pmap->lnSNRsd[c][NumSites[c]+1] = 0.0;
		  pmap->SNRdist[c][NumSites[c]+1] = NULL;
		}

		pmap->blockmem = 0;
		FLOAT *Z = pmap->site[c];

		/* copy Xmap from site 0 to site MM[c]+1 */
		int i = 0;
		for(; i <= MM[c]+1; i++){
		  Z[i] = XX[c][i];
		  if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[c][i] = Xmap->siteSD[c][i];
		  pmap->sitecov[c][i] = Xmap->sitecov[c][i];
		  pmap->sitecnt[c][i] = Xmap->sitecnt[c][i];
		  if(pmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][i];
		    pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][i];
		    pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][i];
		    pmap->SNRdist[c][i] = 0;
		    if(pmap->SNRcnt[c][i]){
		      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][i][t];
		    }
		  }
		  if(VERB>=3){
		    printf("Copying site X[c=%d][i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   c, i, MM[c], XX[c][i], Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == NumSites[c]+1);
	      } /* c = 0..colors-1 */
	    }
	  } else { /* copy left end of X and right end of Y */
	    int kmin = (KK[0]+KK[1])/2;/* default transition is at alignIL[kmin] */
	    if(PairMergePreferNGS){/* location transition point that retains the most NGS sites (with StdDev==0) */
	      double *VARdiff = new double[KK[0]+KK[1]];
	      for(int k = 0; k < KK[0]+KK[1] - 1; k++){
		int u = alignIL[k+1].index;
		if(u <= 0){
		  VARdiff[k] = 0.0;
		  continue;
		}
		int c = alignIL[k+1].c;

		int i = p[c].sites1[u-1];
		int j = p[c].sites2[u-1];
		int I = p[c].sites1[u];
		int J = p[c].sites2[u];
		double varY = 0, varX = 0.0;
		for(int t = i; t < I; t++){
		  double SDy = Ymap->siteSD[c][t];
		  varY += SDy*SDy;
		}
		for(int t = j; t < J; t++){
		  double SDx = p->orientation ? Xmap->siteSD[c][MM[c]-t] : Xmap->siteSD[c][t];
		  varX += SDx*SDx;
		}
		VARdiff[k] = varX - varY;
	      }

	      double VARmin = 0.0, VARsum = 0.0;
	      int origkmin = kmin;
	      kmin = 0;
	      for(int k = 0; k < KK[0]+KK[1] - 1; k++){
		VARsum += VARdiff[k];
		if(VARsum < VARmin){
		  VARmin = VARsum;
		  kmin = k+1;
		}
	      }
	      delete [] VARdiff;
	      if(VERB/* HERE >=2 */){
		int c1 = alignIL[origkmin].c;
		int u1 = alignIL[origkmin].index;
		int I1 = p[c1].sites1[u1];
		int J1 = p[c1].sites1[u1];
		int cmin = alignIL[kmin].c;
		int umin = alignIL[kmin].index;
		int Imin = p[cmin].sites1[umin];
		int Jmin = p[cmin].sites2[umin];

		printf("Changing crossover point from (YY[%d][%d]=%0.3f,X[%d]=%0.3f) to (YY[%d][%d]=%0.3f,X[%d]=%0.3f)\n",
		       c1, I1, Ymap->site[c1][I1], p->orientation ? (MM[c1]+1-J1) : J1, Xmap->site[c1][p->orientation ? (MM[c1]+1-J1) : J1], 
		       cmin, Imin, Ymap->site[cmin][Imin], p->orientation ? (MM[cmin]+1-Jmin) : Jmin, Xmap->site[cmin][p->orientation ? (MM[cmin]+1-Jmin) : Jmin]);
		fflush(stdout);
	      }
	    }

	    int xmin[2],ymin[2];
	    int cmin = alignIL[kmin].c;
	    int umin = alignIL[kmin].index;
	    ymin[cmin] = p[cmin].sites1[umin];
	    xmin[cmin] = p[cmin].sites2[umin];
	    FLOAT Ymin = YY[cmin][ymin[cmin]];
	    FLOAT Xmin = p->orientation ? -XX[cmin][MM[cmin]+1-xmin[cmin]] : XX[cmin][xmin[cmin]];

	    /* compute nearest previous site for other colors c == 1-cmin */
	    int c = 1-cmin;
	    if(DEBUG) ymin[c] = xmin[c] = -1;	    

	    for(int I = 1; I <= NN[c]+1; I++)
	      if(YY[c][I] > Ymin){
		ymin[c] = I-1;
		break;
	      }
	    if(DEBUG) assert(ymin[c] >= 0);

	    for(int J = 1; J <= MM[c]+1; J++)
	      if((p->orientation ? -XX[c][MM[c]+1-J] : XX[c][J]) > Xmin){
		xmin[c] = J-1;
		break;
	      }
	    if(DEBUG) assert(xmin[c] >= 0);
	    
	    /* copy X from left end to XX[cmin][xmin[cmin]], then Y from YY[cmin][ymin[cmin]] to right end */
	    FLOAT Ystart = (p->orientation ? XX[cmin][MM[cmin]+1] - XX[cmin][MM[cmin]+1-xmin[cmin]] : XX[cmin][xmin[cmin]]) - Ymin;

	    /* compute number of sites and allocate memory */
	    int NumSites[2];
	    for(int c = 0; c < colors; c++){
	      NumSites[c] = (p->orientation ? (MM[c]+1)-(MM[c] + 1 - xmin[c]) : xmin[c]) + NN[c] - (ymin[c]) ;
	      if(pmap->site[c]){
		if(!pmap->blockmem) delete [] pmap->site[c];
		if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
		if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
		if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
	      }
	      if(pmap->SNRcnt[c]){
		if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		if(pmap->SNRdist[c]){
		  for(int i = 0; i <= pmap->numsite[0]+1; i++)
		    if(pmap->SNRdist[c][i]) delete [] pmap->SNRdist[c][i];
		  delete [] pmap->SNRdist[c];
		}
		delete [] pmap->SNRcnt[c];
	      }
	      pmap->numsite[c] = NumSites[c];
	      pmap->site[c] = new FLOAT[NumSites[c]+2];
	      pmap->siteSD[c] = new double[NumSites[c]+2];
	      pmap->sitecov[c] = new float[NumSites[c]+2];
	      pmap->sitecnt[c] = new float[NumSites[c]+2];
	      pmap->SNRcnt[c] = 0;
	      if(Xmap->SNRcnt[c] || Ymap->SNRcnt[c]){
		pmap->SNRcnt[c] = new int[NumSites[c]+2];
		pmap->SNRgmean[c] = new double[NumSites[c]+2];
		pmap->lnSNRsd[c] = new double[NumSites[c]+2];
		pmap->SNRdist[c] = new double*[NumSites[c]+2];
	      }

	      /* left end */
	      pmap->site[c][0] = 0.0;
	      pmap->siteSD[c][0] = 0.0;
	      pmap->sitecov[c][0] = 1;
	      pmap->sitecnt[c][0] = 1;
	      if(pmap->SNRcnt[c]){
		pmap->SNRcnt[c][0] = 0;
		pmap->SNRgmean[c][0] = 0.0;
		pmap->lnSNRsd[c][0] = 0.0;
		pmap->SNRdist[c][0] = NULL;
	      }

	      /* right end */
	      pmap->siteSD[c][NumSites[c]+1] = 0.0;
	      pmap->sitecov[c][NumSites[c]+1] = 1;
	      pmap->sitecnt[c][NumSites[c]+1] = 1;
	      if(pmap->SNRcnt[c]){
		pmap->SNRcnt[c][NumSites[c]+1] = 0;
		pmap->SNRgmean[c][NumSites[c]+1] = 0.0;
		pmap->lnSNRsd[c][NumSites[c]+1] = 0.0;
		pmap->SNRdist[c][NumSites[c]+1] = NULL;
	      }

	      pmap->blockmem = 0;
	      FLOAT *Z = pmap->site[c];

	      int i = 0;
	      if(p->orientation){ /* copy left end of Xmap from site MM[c]+1 downto MM[c]+1-xmin[c] */
		double Xstart = XX[c][MM[c]+1];
		for(int j = MM[c]+1; j >= MM[c]+1-xmin[c]; j--, i++){
		  Z[i] = Xstart - XX[c][j];
		  if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[c][i] = Xmap->siteSD[c][j-1]; /* see definition of siteSD[] */
		  pmap->sitecov[c][i] = Xmap->sitecov[c][j];
		  pmap->sitecnt[c][i] = Xmap->sitecnt[c][j];
		  if(Xmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][j];
		    pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][j];
		    pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][j];
		    pmap->SNRdist[c][i] = NULL;
		    if(pmap->SNRcnt[c][i]){
		      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][j][t];
		    }
		  } else if(pmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = 0;
		    pmap->SNRgmean[c][i] = 0.0;
		    pmap->lnSNRsd[c][i] = 0.0;
		    pmap->SNRdist[c][i] = NULL;
		  }
		  if(VERB>=3){
		    printf("Copying site X[c=%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   c, j, MM[c], XX[c][j], Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == (MM[c]+1)-(MM[c]+1-xmin[c]));
	      } else { /* copy left end of Xmap from site 0 to J */
		for(; i <= xmin[c]; i++){
		  Z[i] = XX[c][i];
		  if(DEBUG) assert(Z[i] >= Z[i-1]);
		  pmap->siteSD[c][i] = Xmap->siteSD[c][i];
		  pmap->sitecov[c][i] = Xmap->sitecov[c][i];
		  pmap->sitecnt[c][i] = Xmap->sitecnt[c][i];
		  if(Xmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = Xmap->SNRcnt[c][i];
		    pmap->SNRgmean[c][i] = Xmap->SNRgmean[c][i];
		    pmap->lnSNRsd[c][i] = Xmap->lnSNRsd[c][i];
		    pmap->SNRdist[c][i] = NULL;
		    if(pmap->SNRcnt[c][i]){
		      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
			pmap->SNRdist[c][i][t] = Xmap->SNRdist[c][i][t];
		    }
		  } else if(pmap->SNRcnt[c]){
		    pmap->SNRcnt[c][i] = 0;
		    pmap->SNRgmean[c][i] = 0.0;
		    pmap->lnSNRsd[c][i] = 0.0;
		    pmap->SNRdist[c][i] = NULL;
		  }
		  if(VERB>=3){
		    printf("Copying site X[c=%d][i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			   c, i, MM[c], XX[c][i], Xmap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		    fflush(stdout);
		  }
		}
		if(DEBUG) assert(i-1 == xmin[c]);
	      } // else p->orientation == 0

	      /* copy right end of Ymap from site ymin[c] to NN[c] + 1 */
	      int j = ymin[c];
	      // double Ystart = Z[i-1] - Y[j];
	      for(; ++j <= NN[c] + 1; i++){
		Z[i] = Ystart + YY[c][j];
		if(DEBUG) assert(Z[i] >= Z[i-1]);
		pmap->siteSD[c][i-1] = Ymap->siteSD[c][j-1];/* see definition of siteSD[] */
		pmap->sitecov[c][i] = Ymap->sitecov[c][j];
		pmap->sitecnt[c][i] = Ymap->sitecnt[c][j];
		if(Ymap->SNRcnt[c]){
		  pmap->SNRcnt[c][i] = Ymap->SNRcnt[c][j];
		  pmap->SNRgmean[c][i] = Ymap->SNRgmean[c][j];
		  pmap->lnSNRsd[c][i] = Ymap->lnSNRsd[c][j];
		  pmap->SNRdist[c][i] = NULL;
		  if(pmap->SNRcnt[c][i]){
		    pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		    for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
		      pmap->SNRdist[c][i][t] = Ymap->SNRdist[c][j][t];
		  }
		} else if(pmap->SNRcnt[c]){
		  pmap->SNRcnt[c][i] = 0;
		  pmap->SNRgmean[c][i] = 0.0;
		  pmap->lnSNRsd[c][i] = 0.0;
		  pmap->SNRdist[c][i] = NULL;
		}
		if(VERB>=3){
		  printf("Copying site Y[c=%d][j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
			 c, j, NN[c], YY[c][j], Ymap->id, i, pmap->numsite[c], Z[i], pmap->sitecov[c][i], pmap->sitecnt[c][i]);
		  fflush(stdout);
		}
	      }
	      if(DEBUG) assert(i-1 == NumSites[c]+1);
	    }// c = 0 .. colors-1
	  } // copy left end of X and right end of Y

	  pmap->len = pmap->site[0][pmap->numsite[0]+1];
	  if(VERB){
	    if(PairMergeRepeat)
	      printf("%d:",pairmergeIter);
	    printf("Created merged map (id=%lld,len=%0.3f) based on alignment between Y mapid=%d(id=%lld,len=%0.3f) and X mapid=%d(id=%lld,len=%0.3f):mapid=%d,id=%lld,len=%0.3f,sites=%d\n",
		   pmap->id, pmap->len, Ymap->mapid,Ymap->id, Ymap->site[0][Ymap->numsite[0]+1], 
		   Xmap->mapid,Xmap->id, Xmap->site[0][Xmap->numsite[0]+1],pmap->mapid,pmap->id,pmap->len,pmap->numsite[0]);
	    //	    printf(" Map Merge: Yid= %lld Xid= %lld or= %d : id= %lld\n",Ymap->id,Xmap->id,p->orientation, pmap->id);
#ifdef ASAN
	    fflush(stdout);
#endif
	  }
	} // left end of X comes first
	delete [] Ysite;
	delete [] Xsite;
	delete [] Yindex[0];
	delete [] Xindex[0];
      } // if(PairMerge ... )
    }// for i = 0 .. numaligns-1
  } // colors==2

  if(PairMergeHmap <= 0)
    return splitcnt;

  if(VERB){
    printf("Checking for Haplotype end overlaps\n");
    fflush(stdout);
  }

  int iter= 0, progress = 1;

  while(progress){
    progress = 0;
    iter++;
    int overlapcnt = 0;

    for(int i = 0; i < nummaps; i++){
      Cmap *pmap = Gmap[i];
      if(pmap->paired)
	continue;
      if(!pmap->leftoverlap && !pmap->rightoverlap)
	continue;
      if(pmap->leftoverlap && pmap->rightoverlap){
	overlapcnt++;
	continue;
      }

      if(pmap->leftoverlap){
	if(DEBUG) assert(!pmap->rightoverlap);
	Calign *p = pmap->leftoverlap;
	
	//	int K = p->numpairs;
	Cmap *Ymap = YYmap[p->mapid1];
	Cmap *Xmap = XXmap[p->mapid2];
	int N = Ymap->numsite[0];
	int M = Xmap->numsite[0];
	//	FLOAT *Y = Ymap->site[0];
	//	FLOAT *X = Xmap->site[0];

	if(DEBUG>=1+RELEASE) assert(pmap == Xmap || pmap == Ymap);
	if(pmap == Xmap && Ymap->paired){
	  if(VERB>=2){
	    printf("NOT merging Ymap->id=%lld with Xmap->id=%lld into Hmap based on %s end of Ymap overlapping left end of Xmap:Ymap->paired=%d(Xmap->contig=%p,Ymap->contig=%p)\n",
		   Ymap->id,Xmap->id, p->orientation ? "left" : "right", Ymap->paired,Xmap->contig,Ymap->contig);
	    fflush(stdout);
	  }
	  continue;// NEW79
	}
	if(pmap == Ymap && Xmap->paired){
	  if(VERB>=2){
	    printf("NOT merging Xmap->id=%lld with Ymap->id=%lld (or=%d) into Hmap based on left end of Ymap overlapping %s end of Xmap:Xmap->paired=%d(Xmap->contig=%p,Ymap->contig=%p)\n",
		   Xmap->id,Ymap->id,p->orientation, p->orientation ? "left" : "right",Xmap->paired,Xmap->contig,Ymap->contig);
	    fflush(stdout);
	  }
	  continue;// NEW79
	}

	if(!Ymap->contig){
	  if(DEBUG) assert(!(Ymap->centercloned & 2));
	  Ymap->contig = new Ccontig(Ymap);/* convert Ymap into a trivial Hmap with no Het events */
	}
	if(!Xmap->contig){
	  if(DEBUG) assert(!(Xmap->centercloned & 2));
	  Xmap->contig = new Ccontig(Xmap);/* convert Xmap into a trivial Hmap with no Het events */
	}

	if(pmap == Xmap){/* merge Xmap/pmap into Ymap */
	  if(DEBUG) assert((p->orientation ? Ymap->leftoverlap : Ymap->rightoverlap) == p);
	  if(VERB/* HERE >=2 */){
	    printf("Merging Ymap->id=%lld with Xmap->id=%lld into Hmap based on %s end of Ymap overlapping left end of Xmap\n",Ymap->id,Xmap->id, p->orientation ? "left" : "right");
	    fflush(stdout);
	  }
	  if(merge_contigs(Ymap->contig,Xmap->contig,p, outlierExtend)){/* failed to merge */
	    if(VERB){
	      printf("Failed to merge Ymap->id=%lld with Xmap->id=%lld into Hmap\n",Ymap->id,Xmap->id);
	      fflush(stdout);
	    }
	    if(!(Xmap->centercloned & 2)){
	      delete Xmap->contig;
	      Xmap->contig = NULL;
	    }
	    if(!(Ymap->centercloned & 2)){
	      delete Ymap->contig;
	      Ymap->contig = NULL;
	    }
	  } else {
	    Ymap->contig->id = Ymap->id;// WAS min(Ymap->id, Xmap->id);
	    if(VERB){
	      int N = Ymap->contig->numsite[0];
	      printf("Merged Ymap->mapid=%d(id=%lld,paired=%d) with Xmap->mapid=%d(id=%lld,paired=%d->%d) into Hmap(Ymap->contig->id=%lld,sites=%d,len=%0.3f)\n",
		     Ymap->mapid,Ymap->id,Ymap->paired,Xmap->mapid,Xmap->id,Xmap->paired,1+p->mapid1, Ymap->contig->id,N,Ymap->contig->site[0][N+1]);
	      fflush(stdout);
	    }
	    Xmap->paired = 1 + p->mapid1;/* Mark as used up */
	    Ymap->centercloned |= 2;/* Mark as modified/Hmap */
	    progress++;
	  }

	  pmap->leftoverlap = NULL;
	  if(p->orientation)
	    Ymap->leftoverlap = NULL;
	  else
	    Ymap->rightoverlap = NULL;

	} else /* pmap == Ymap */{/* merge Ymap/pmap into Xmap */
	  if(VERB/* HERE >=2 */){
	    printf("Merging Xmap->id=%lld with Ymap->id=%lld (or=%d) into Hmap based on left end of Ymap overlapping %s end of Xmap\n",Xmap->id,Ymap->id,p->orientation, p->orientation ? "left" : "right");
	    fflush(stdout);
	  }	  
	  if(DEBUG) assert(pmap == Ymap);
	  if(DEBUG) assert((p->orientation ? Xmap->leftoverlap : Xmap->rightoverlap) == p);

	  Calign *q = new Calign;
	  q->mapid1 = p->mapid2;
	  q->mapid2 = p->mapid1;
	  int U = p->numpairs;
	  q->numpairs = U;
	  q->orientation = p->orientation;

	  if(q->orientation){/* need to create alignment of Xmap with flipped Ymap */
	    if(VERB/* HERE >=2 */){
	      printf("Merging Xmap->id=%lld with flipped Ymap->id=%lld into Hmap\n",Xmap->id,Ymap->id);
	      fflush(stdout);
	    }

	    q->allocated_size = U;
	    q->sites1 = new int[U];
	    q->sites2 = new int[U];
	    q->iscore = new FLOAT[U+1];
	    q->outscore = new FLOAT[U+1];

	    for(int t = 0; t < U; t++){
	      q->sites1[U-1-t] = M + 1 - p->sites2[t];
	      q->sites2[U-1-t] = N + 1 - p->sites1[t];
	    }
	    for(int t = 0; t <= U; t++){
	      q->iscore[U-t] = p->iscore[t];
	      q->outscore[U-t] = p->outscore[t];
	    }
	    if(merge_contigs(Xmap->contig,Ymap->contig,q,outlierExtend)){
	      if(VERB){
		printf("Failed to merge Xmap->id=%lld with flipped Ymap->id=%lld (or=%d) into Hmap\n",Xmap->id,Ymap->id, p->orientation);
		fflush(stdout);
	      }
	      if(!(Xmap->centercloned & 2)){
		delete Xmap->contig;
		Xmap->contig = NULL;
	      }
	      if(!(Ymap->centercloned & 2)){
		delete Ymap->contig;
		Ymap->contig = NULL;
	      }
	    } else {
	      Xmap->contig->id = Xmap->id; // WAS min(Xmap->id, Ymap->id);
	      if(VERB){
		int M = Xmap->contig->numsite[0];
		printf("Merged Xmap->mapid=%d(id=%lld,paired=%d) with flipped Ymap->mapid=%d(id=%lld,paired=%d->%d) (or=%d) as Hmap(Xmap->contig->id=%lld,sites=%d,len=%0.3f)\n",
		       Xmap->mapid,Xmap->id,Xmap->paired,Ymap->mapid,Ymap->id,Ymap->paired,1+p->mapid2, p->orientation, Xmap->contig->id,M,Xmap->contig->site[0][M+1]);
		fflush(stdout);
	      }
	      Ymap->paired = 1 + p->mapid2;
	      Xmap->centercloned |= 2;
	      // WAS	      Xmap->id = min(Xmap->id, Ymap->id);
	      progress++;
	    }
	  } else /* !q->orientation */{
	    if(VERB/* HERE >=2 */){
	      printf("Merging Xmap->id=%lld with unflipped Ymap->id=%lld into Hmap\n",Xmap->id,Ymap->id);
	      fflush(stdout);
	    }
	    /* first create a copy of p == alignment[i] in which sites1[] and sites2[] are swapped */
	    q->sites1 = p->sites2;
	    q->sites2 = p->sites1;
	    q->iscore = p->iscore;
	    q->outscore = p->outscore;
	    
	    if(merge_contigs(Xmap->contig,Ymap->contig,q,outlierExtend)){
	      if(VERB){
		printf("Failed to merge Xmap->id=%lld with unflipped Ymap->id=%lld (or=%d) as Hmap\n",Xmap->id,Ymap->id, p->orientation);
		fflush(stdout);
	      }
	      if(!(Xmap->centercloned & 2)){
		delete Xmap->contig;
		Xmap->contig = NULL;
	      }
	      if(!(Ymap->centercloned & 2)){
		delete Ymap->contig;
		Ymap->contig = NULL;
	      }
	    } else {
	      Xmap->contig->id = Xmap->id; // WAS min(Xmap->id, Ymap->id);
	      if(VERB){
		int M = Xmap->contig->numsite[0];
		printf("Merged Xmap->mapid=%d(id=%lld,paired=%d) with unflipped Ymap->mapid=%d(id=%lld,paired=%d->%d) (or=%d) as Hmap(Xmap->contig->id=%lld,sites=%d,len=%0.3f)\n",
		       Xmap->mapid,Xmap->id,Xmap->paired,Ymap->mapid,Ymap->id,Ymap->paired,1+p->mapid2, p->orientation, Xmap->contig->id,M,Xmap->contig->site[0][M+1]);
		fflush(stdout);
	      }
	      Ymap->paired = 1 + p->mapid2;
	      Xmap->centercloned |= 2;
	      // WAS	      Xmap->id = min(Xmap->id, Ymap->id);
	      progress++;
	    }
	    q->sites1 = q->sites2 = NULL;// should not be needed since q->allocated_size = 0;
	  }
	  delete q;

	  pmap->leftoverlap = NULL;
	  if(p->orientation)
	    Xmap->leftoverlap = NULL;
	  else
	    Xmap->rightoverlap = NULL;

	} // pmap == Ymap

      } else { /* pmap->rightoverlap */
	if(DEBUG) assert(pmap->rightoverlap);
	Calign* p = pmap->rightoverlap;
	
	//	int K = p->numpairs;
	Cmap *Ymap = YYmap[p->mapid1];
	Cmap *Xmap = XXmap[p->mapid2];
	int N = Ymap->numsite[0];
	int M = Xmap->numsite[0];
	//	FLOAT *Y = Ymap->site[0];
	//	FLOAT *X = Xmap->site[0];

	if(DEBUG>=1+RELEASE) assert(pmap == Xmap || pmap == Ymap);	
	if(pmap == Xmap && Ymap->paired){
	  if(VERB>=2){
	    printf("NOT Merging Ymap->id=%lld with Xmap->id=%lld (or=%d) into Hmap based on %s end of Ymap (%p) overlapping right end of Xmap (p= %p):Ymap->paired=%d(Xmap->contig=%p,Ymap->contig=%p)\n", 
		   Ymap->id,Xmap->id, p->orientation, p->orientation ? "right" : "left", p->orientation ? Ymap->rightoverlap : Ymap->leftoverlap, p, Ymap->paired,Xmap->contig,Ymap->contig);
	    fflush(stdout);
	  }
	  continue;// NEW79
	}
	if(pmap == Ymap && Xmap->paired){
	  if(VERB>=2){
	    printf("NOT Merging Xmap->id=%lld with Ymap->id=%lld into Hmap based on right end of Ymap overlapping %s end of Xmap:Xmap->paired=%d(Xmap->contig=%p,Ymap->contig=%p)\n",
		   Xmap->id,Ymap->id,p->orientation ? "right" : "left",Xmap->paired,Xmap->contig,Ymap->contig);
	    fflush(stdout);
	  }	  
	  continue;// NEW79
	}

	if(!Ymap->contig){
	  if(DEBUG) assert(!(Ymap->centercloned & 2));
	  Ymap->contig = new Ccontig(Ymap);/* convert Ymap into a trivial Hmap with no Het events */
	}
	if(!Xmap->contig){
	  if(DEBUG) assert(!(Xmap->centercloned & 2));
	  Xmap->contig = new Ccontig(Xmap);/* convert Xmap into a trivial Hmap with no Het events */
	}

	if(pmap == Xmap){/* merge Xmap/pmap into Ymap */
	  if(VERB/* HERE >=2 */){
	    printf("Merging Ymap->id=%lld with Xmap->id=%lld (or=%d) into Hmap based on %s end of Ymap (%p) overlapping right end of Xmap (p= %p)\n", 
		   Ymap->id,Xmap->id, p->orientation, p->orientation ? "right" : "left", p->orientation ? Ymap->rightoverlap : Ymap->leftoverlap, p);
	    fflush(stdout);
	  }
	  if(DEBUG) assert((p->orientation ? Ymap->rightoverlap : Ymap->leftoverlap) == p);
	  if(merge_contigs(Ymap->contig,Xmap->contig,p, outlierExtend)){/* failed to merge */
	    if(VERB){
	      printf("Failed to merge Ymap->id=%lld with Xmap->id=%lld into Hmap\n",Ymap->id,Xmap->id);
	      fflush(stdout);
	    }
	    if(!(Xmap->centercloned & 2)){
	      delete Xmap->contig;
	      Xmap->contig = NULL;
	    }
	    if(!(Ymap->centercloned & 2)){
	      delete Ymap->contig;
	      Ymap->contig = NULL;
	    }
	  } else {
	    Ymap->contig->id = Ymap->id; // WAS min(Ymap->id, Xmap->id);
	    if(VERB){
	      int N = Ymap->contig->numsite[0];
	      printf("Merged Ymap->mapid=%d(id=%lld,paired=%d) with Xmap->mapid=%d(id=%lld,paired=%d->%d) into Hmap(Ymap->contig->id=%lld,sites=%d,len=%0.2f)\n",
		     Ymap->mapid,Ymap->id,Ymap->paired,Xmap->mapid,Xmap->id,Xmap->paired,1+p->mapid1, Ymap->contig->id,N,Ymap->contig->site[0][N+1]);
	      fflush(stdout);
	    }
	    Xmap->paired = 1 + p->mapid1;/* Mark as used up */
	    Ymap->centercloned |= 2;/* Mark as modified/Hmap */
	    // WAS	    Ymap->id = min(Ymap->id, Xmap->id);
	    progress++;
	  }

	  pmap->rightoverlap = NULL;
	  if(p->orientation)
	    Ymap->rightoverlap = NULL;
	  else
	    Ymap->leftoverlap = NULL;

	} else /* pmap == Ymap */{/* merge Ymap/pmap into Xmap */
	  if(DEBUG) assert(pmap == Ymap);

	  if(DEBUG) assert((p->orientation ? Xmap->rightoverlap : Xmap->leftoverlap) == p);
	  if(VERB/* HERE >=2 */){
	    printf("Merging Xmap->id=%lld with Ymap->id=%lld into Hmap based on right end of Ymap overlapping %s end of Xmap\n",Xmap->id,Ymap->id,p->orientation ? "right" : "left");
	    fflush(stdout);
	  }	  

	  Calign *q = new Calign;
	  q->mapid1 = p->mapid2;
	  q->mapid2 = p->mapid1;
	  int U = p->numpairs;
	  q->numpairs = U;
	  q->orientation = p->orientation;

	  if(q->orientation){/* need to create alignment of Xmap with flipped Ymap */
	    if(VERB/* HERE >=2 */){
	      printf("Merging Xmap->id=%lld with flipped Ymap->id=%lld into Hmap\n",Xmap->id,Ymap->id);
	      fflush(stdout);
	    }

	    q->allocated_size = U;
	    q->sites1 = new int[U];
	    q->sites2 = new int[U];
	    q->iscore = new FLOAT[U+1];
	    q->outscore = new FLOAT[U+1];

	    for(int t = 0; t < U; t++){
	      q->sites1[U-1-t] = M + 1 - p->sites2[t];
	      q->sites2[U-1-t] = N + 1 - p->sites1[t];
	    }
	    for(int t = 0; t <= U; t++){
	      q->iscore[U-t] = p->iscore[t];
	      q->outscore[U-t] = p->outscore[t];
	    }
	    if(merge_contigs(Xmap->contig,Ymap->contig,q,outlierExtend)){
	      if(VERB){
		printf("Failed to merge Xmap->id=%lld with flipped Ymap->id=%lld (or=%d) into Hmap\n",Xmap->id,Ymap->id, p->orientation);
		fflush(stdout);
	      }
	      if(!(Xmap->centercloned & 2)){
		delete Xmap->contig;
		Xmap->contig = NULL;
	      }
	      if(!(Ymap->centercloned & 2)){
		delete Ymap->contig;
		Ymap->contig = NULL;
	      }
	    } else {
	      Xmap->contig->id = Xmap->id; // WAS min(Xmap->id, Ymap->id);
	      if(VERB){
		int M = Xmap->contig->numsite[0];
		printf("Merged Xmap->mapid=%d(id=%lld,paired=%d) with flipped Ymap->mapid=%d(id=%lld,paired=%d->%d) (or=%d) as Hmap(Xmap->contig->id=%lld,sites=%d,len=%0.3f)\n",
		       Xmap->mapid,Xmap->id,Xmap->paired,Ymap->mapid,Ymap->id,Ymap->paired,1+p->mapid2, p->orientation, Xmap->contig->id,M,Xmap->contig->site[0][M+1]);
		fflush(stdout);
	      }
	      Ymap->paired = 1 + p->mapid2;
	      Xmap->centercloned |= 2;
	      // WAS	      Xmap->id = min(Xmap->id, Ymap->id);
	      progress++;
	    }
	  } else {
	    if(VERB/* HERE >=2 */){
	      printf("Merging Xmap->id=%lld with unflipped Ymap->id=%lld into Hmap\n",Xmap->id,Ymap->id);
	      fflush(stdout);
	    }
	    /* first create a copy of p == alignment[i] in which sites1[] and sites2[] are swapped */
	    q->sites1 = p->sites2;
	    q->sites2 = p->sites1;
	    q->iscore = p->iscore;
	    q->outscore = p->outscore;
	    
	    if(merge_contigs(Xmap->contig,Ymap->contig,q,outlierExtend)){
	      if(VERB){
		printf("Failed to merge Xmap->id=%lld with unflipped Ymap->id=%lld (or=%d) as Hmap\n",Xmap->id,Ymap->id, p->orientation);
		fflush(stdout);
	      }
	      if(!(Xmap->centercloned & 2)){
		delete Xmap->contig;
		Xmap->contig = NULL;
	      }
	      if(!(Ymap->centercloned & 2)){
		delete Ymap->contig;
		Ymap->contig = NULL;
	      }
	    } else {
	      Xmap->contig->id = Xmap->id; // WAS min(Xmap->id, Ymap->id);
	      if(VERB){
		int M = Xmap->contig->numsite[0];
		printf("Merged Xmap->mapid=%d(id=%lld,paired=%d) with unflipped Ymap->mapid=%d(id=%lld,paired=%d->%d) (or=%d) as Hmap(Xmap->contig->id=%lld,sites=%d,len=%0.3f)\n",
		       Xmap->mapid,Xmap->id,Xmap->paired,Ymap->mapid,Ymap->id,Ymap->paired,1+p->mapid2, p->orientation, Xmap->contig->id,M,Xmap->contig->site[0][M+1]);
		fflush(stdout);
	      }
	      Ymap->paired = 1 + p->mapid2;
	      Xmap->centercloned |= 2;
	      // WAS	      Xmap->id = min(Xmap->id, Ymap->id);
	      progress++;
	    }
	    q->sites1 = q->sites2 = NULL;// should not be needed since q->allocated_size = 0;
	  }
	  delete q;

	  pmap->rightoverlap = NULL;
	  if(p->orientation)
	    Xmap->rightoverlap = NULL;
	  else
	    Xmap->leftoverlap = NULL;

	} // pmap == Ymap

      } // pmap->rightoverlap
    } // i = 0 .. nummaps-1

    if(VERB){
      printf("iter=%d: found %d maps with double overlap, merged %d maps\n", iter, overlapcnt, progress);
      fflush(stdout);
    }
  } // while(progress)

  return splitcnt;
}

static int threads_displayed = 0;

static Cmap **MergedMaps = NULL;/* used with PairMergeIDrenumber to save pointers to merged maps : they have links (origmap) to output contig that includes them (if any) */
static int MergedNummaps = 0, MergedMaxmaps = 0;
class CID {
public:
  long long newID;
  long long oldID;
};

static int CIDinc(CID *p1, CID *p2)
{
  return (p1->newID > p2->newID) ? 1 : (p1->newID < p2->newID) ? -1 : 
    (p1->oldID > p2->oldID) ? 1 : (p1->oldID < p2->oldID) ? -1 : 0;
}

//static char buf[LINESIZ];

static int A_numthreads= -1;
static long long A_maxN= -1;
static AL *A_threads = NULL;

static int mem_numthreads= -1;
static long long mem_maxN= -1, mem_maxM= -1;

/** Find best alignment between each nanomap pair (in either orientation for the map). For each
   pair the best alignment is appended to the alignment list alignment[0..numaligns]. 
   Long function > 2000 lines.
   
   - Condition input, injest hash results
   - Allocate memory for maps, multithreaded dynamic programming.
   - Perform all-to-all alignment, call pairalign() for each nanomap pair.
      + Pairsplit, and chimeric functionality for sub-map alignments
   - Calculate statistics
   - Output align file.
    */

int pairalign_allpairs(int mapstart,int mapend, int pairmergeIter, int pairmergeIterLast)
{
  double wtstart = wtime();

  if(RA_MIN_MEM >= 8 && !RA_MIN_TIME){
    printf("-RAmem %d %d not supported : using -RAmem %d 1 instead\n",RA_MIN_MEM,RA_MIN_TIME,RA_MIN_MEM);
    fflush(stdout);

    RA_MIN_TIME = 1;
  }

  if(numbins > 1 || skipcnt > 0 || bin != 1){
    printf("-partial is no longer supported\n");
    fflush(stdout);exit(1);
  }

  if(colors > 1){
    if(colors > 2){
      printf("pairalign_allpairs: colors(=%d) > 2 not supported\n",colors);
      fflush(stdout);exit(1);
    }
    if(DEBUG) assert(!usecolor);
    return pairalign_allpairs2(mapstart, mapend, pairmergeIter);
  }

  if(AlignedEndOutlierThreshold < 2){
    printf("-E not yet implemented for pairwise alignment\n");
    fflush(stdout);exit(1);
  }
  if(AlignedOutlierThreshold < 999.0){
    printf("-I not yet implemented for pairwise alignment\n");
    fflush(stdout);exit(1);
  }

  if(DEBUG) assert(mapstart == 0 && mapend == nummaps - 1);
  repeatcnt = 0;

  if(DEBUG>=2 && !hash_filename && !(PairMergeRepeat && pairmergeIter > 1)){
    for(int mapid = 0; mapid < nummaps; mapid++){
      if(!(Gmap[mapid]->mapid == mapid)){
	printf("mapid=%d/%d:Gmap[mapid]=%p,Gmap[mapid]->mapid=%d,Gmap[mapid]->id=%lld\n",
	       mapid,nummaps,Gmap[mapid],Gmap[mapid]->mapid,Gmap[mapid]->id);
	for(int m = 0; m < min(nummaps,10); m++)
	  printf("m=%d:Gmap[m]=%p,Gmap[m]->mapid=%d,Gmap[m]->id=%lld\n",
		 m,Gmap[m],Gmap[m]->mapid,Gmap[m]->id);
	fflush(stdout);
	assert(Gmap[mapid]->mapid == mapid);
      }
    }
  }

  if(!PairSplit)
    PairsplitSeparateQueries = 0;

  if(usecolor && PairSplit){
    printf("2 color -i input with -usecolor %d and -pairsplit (or -sv) not yet supported\n", usecolor);
    fflush(stdout);exit(1);
  }

  if(VERB>=2){
    printf("pairalign:PVERB=%d,LVERB=%d\n",PVERB,LVERB);
    fflush(stdout);
  }

  if(ROC){
    printf("ROC=%d no longer supported\n", ROC);
    fflush(stdout);exit(1);
  }
  if(DEBUG && DEFER_G != 2){
    printf("DEFER_G=%d not supported\n",DEFER_G);
    fflush(stdout);exit(1);
  }
  if(DEBUG && !DEFER_BP){
    printf("DEFER_BP=%d not supported\n",DEFER_BP);
    fflush(stdout);exit(1);
  }
  if(DEBUG && !LOCAL){
    printf("LOCAL=%d not supported\n",LOCAL);
    fflush(stdout);exit(1);
  }
  /*  if(DEBUG && !ENDFIX){
    printf("ENDFIX=%d not supported\n",ENDFIX);
    fflush(stdout);exit(1);
    }*/
  if(DEBUG && PAIRSCORETYPE != 5){
    printf("PAISCORETYPE=%d not supported(must be 5)\n",PAIRSCORETYPE);
    fflush(stdout);exit(1);
  }


  if(!Kbiasinit)
    Kbias = 1.5;
  if(!KLbiasinit)
    KLbias = 1.5;
  if(!KFbiasinit)
    KFbias = 1.5;
  if(!Gbiasinit)
    Gbias = 1.0;
  if(!Betainit)
    Beta = 0.0;

  if(resSD[0] > 0.0){
    printf("resSD = %0.6f : Not implemented for pairwise alignment: using resSD=0\n", resSD[0]);
    resSD[0] = 0.0;
  }

  /* adjust initial values based on bounds */
  if(FP[0] < MINFP){
    printf("FP=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",FP[0],MINFP,MAXFP,MINFP);
    FP[0] = MINFP;
  } else if(FP[0] > MAXFP){
    printf("FP=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",FP[0],MINFP,MAXFP,MAXFP);
    FP[0] = MAXFP;
  }
  if(FN[0] < MINFN){
    printf("FN=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",FN[0],MINFN,MAXFN,MINFN);
    FN[0] = MINFN;
  } else if(FN[0] > MAXFN){
    printf("FN=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",FN[0],MINFN,MAXFN,MAXFN);
    FN[0] = MAXFN;
  }
  if(SD[0] < MINSD){
    printf("SD=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",SD[0],MINSD,MAXSD,MINSD);
    SD[0] = MINSD;
  } else if(SD[0] > MAXSD){
    printf("SD=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",SD[0],MINSD,MAXSD,MAXSD);
    SD[0] = MAXSD;
  }
  if(SF[0] < MINSF){
    printf("SF=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",SF[0],MINSF,MAXSF,MINSF);
    SF[0] = MINSF;
  } else if(SF[0] > MAXSF){
    printf("SF=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",SF[0],MINSF,MAXSF,MAXSF);
    SF[0] = MAXSF;
  }

  if(PairSplit){    /* initialize PairHash[] */
    PairHash_init();
    if(DEBUG) assert(PairHash);
  }

  if(PoutlierEnd > 0.9999)
    PoutlierEnd = 0.9999;

  /* initialize statistics counter */
  long long mapcnt=0;/* map pairs aligned above threshold */
  double len=0.0;/* sum of aligned lengths (based on overlap region) */
  double lenB = 0.0;/* sum of aligned lengths (based on reference) - resKB2 */
  long long sitecnt=0, misscnt=0;/* count of total and misaligned sites in aligned regions of all pairs */
  long long segcnt=0;/* count of aligned intervals */
  double errnorm=0.0;/* sum of (x-y)*(x-y)/(2*SF[0]*SF[0] + fabs(SD[0])*SD[0]*(x+y) + SR[0]*SR[0]*(x*x+y*y))  */

  long long outliertot = 0;/* all intervals that are candidates for outliers (eg x < y) */
  long long outliercnt = 0;
  double outliersum = 0.0;

  long long poscnt=0;/* count of map pairs with +ve score */

  long long tcnt =0;/* total pairwise alignment attempts */
  long long paircnt = 0;/* number of pairs with minimum overlap of MINOVERLAP * min(len1,len2) */
    
  size_t numerrors = 0;
  size_t maxerrors = 16;
  Cinterval *errors = new Cinterval[maxerrors];

  size_t maxmisaligns = 16;
  Cmisalign *misaligns = new Cmisalign[maxmisaligns];
  size_t nummisaligns = 0;

  gnumthreads = 1;
  int maxthreads = 1;
  #ifdef _OPENMP
  maxthreads = MaxThreads; // WAS omp_get_max_threads();
  //  maxthreads = min(maxthreads,MaxThreads);
  gnumthreads = maxthreads;
  #endif

  score_init(0,nummaps-1);// also computes maxN,maxM over all maps
  
  //  double ExpOutlierPenalty = exp(OutlierPenalty);
  ChimScore = MINSCORE;

  if(Ch > 0.0 && ChFdr > 0.0){
    ChimScore = -log(/* max(0.001,1.0-ChFdr) */ 1.0/(Ch*ChFdr));
    if(ChimScore > 0.0)
      ChimScore = 0.0;
  }
  if(PoutlierEnd > 0.0){
    double OutlierEndPen = OutlierEndBias + OutlierEndPenalty;
    if(OutlierEndPen > ChimScore)
      ChimScore = OutlierEndPen;
  }
  if(PairSplit && ChimScore > -0.001)
    ChimScore = -0.001;// PairSplit assumes ChimScore < 0 
  if(VERB){
    printf("ChimScore=%0.5e:Cfirst=%d(Pairs=%d),PairSplit=%d,Psplit=%0.6e,Poutlier=%0.6f,Ch=%0.6f,ChFdr=%0.4f,PoutlierEnd=%0.6f\n",
	   ChimScore,Cfirst,CfirstPairs,PairSplit,Psplit,Poutlier,Ch,ChFdr,PoutlierEnd);
    fflush(stdout);
  }

  if(PairMerge && PairMergeOutlierChiSq > 0.0 && CmapChimQuality >= 3){/* compute average interval size */
    int cnt = 0;
    double sum = 0.0;

    for(int i = 0; i < nummaps; i++){
      Cmap *pmap = Gmap[i];
      int  M = pmap->numsite[0];
      cnt += M;
      sum += pmap->site[0][M+1];
    }

    AvgInterval = max(PairMergeOutlierMinInterval, sum / max(1,cnt));
    if(VERB){
      printf("AvgInterval = %0.4f kb computed from %d Xmaps\n",AvgInterval,nummaps);
      fflush(stdout);
    }
  }

  if(PairMerge/* WAS17 && pairmergeIter==1*/){/* compute IDmax : may increase with pairmergeIter */
    long long prevIDmax = IDmax;

    IDmax = 0;
    for(int i = 0; i < nummaps; i++)
      IDmax = max(IDmax, Gmap[i]->id);
    long long origIDmax = IDmax;
    if(IDmax < IDMAX){/* round IDmax up to nearest power of 10 (but stop if IDmax exceeds IDMAX) */
      long long val = 10LL;
      while(val < IDmax)
	val *= 10LL;
      IDmax = val;
    }
    if(VERB/* HERE >=2 */){
      if(pairmergeIter==1)
	printf("%d:IDmax = %lld -> %lld\n",pairmergeIter,origIDmax, IDmax);
      else
	printf("%d:IDmax = %lld -> %lld -> %lld\n",pairmergeIter,prevIDmax, origIDmax, IDmax);
      fflush(stdout);
    }
    if(DEBUG) assert(IDmax > 0);
  }

  if(!Cfirst || CfirstPairs == 2){
    if(DEBUG && !Cfirst && nummaps > 1 && !((PairMergeHmap <= 0 && SplitSegDup <= 0.0 && (HashGen || pairmergeIter==1)) || ((PairMergeHmap > 0 || SplitSegDup > 0.0) && pairmergeIter > 1))){
      printf("Cfirst=%d,CfirstPairs=%d,pairmergeIter=%d,PairMergeHmap=%d\n",  Cfirst,CfirstPairs,pairmergeIter,PairMergeHmap);
      fflush(stdout);
      assert((PairMergeHmap <= 0 && SplitSegDup <= 0.0 && (HashGen || pairmergeIter==1)) || ((PairMergeHmap > 0 || SplitSegDup > 0.0) && pairmergeIter > 1));
    }
    if((hash_filename && !HashGen) || (PAIRMERGE_SORTED && PairMerge && !Cfirst))    /* sort maps[0..nummaps-1] by id in ascending order */
      qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapIdInc);
    if(DEBUG) assert(!PairSplit);
    YYmap = XXmap = &Gmap[0];
    numX = numY = nummaps;
    if(Cfirst && CfirstPairs == 2){/* compute separate maxM for XXmap[Cfirst .. numX-1] */
      maxM = 0;
      for(int Xid = Cfirst; Xid < numX; Xid++)
	if(XXmap[Xid]->numsite[0] > maxM)
	  maxM = XXmap[Xid]->numsite[0];
    }
    if(VERB){
      if(Cfirst && CfirstPairs == 2){
	printf("map[0..%d] duplicated as XXmap[0..%d] and YYmap[0..%d] (only XXmap[%d..%d] will be used):maxN=%lld,maxM=%lld (Cfirst=%d,CfirstPairs=%d)\n",
	       nummaps-1,numX-1,numY-1, Cfirst, numX-1,maxN,maxM,Cfirst,CfirstPairs);
	if(VERB>=2 && pairmergeIter > 1){
	  for(int Yid = 0; Yid < numY;Yid++){
	    Cmap *Ymap = YYmap[Yid];
	    int N = Ymap->numsite[0];
	    printf("YYmap[%d]=%p:mapid=%d,id=%lld,len=%0.4f,N=%d\n",Yid,Ymap,Ymap->mapid,Ymap->id,Ymap->site[0][N+1],N);
	  }
	  for(int Xid = Cfirst; Xid < numX;Xid++){
	    Cmap *Xmap = XXmap[Xid];
	    int M = Xmap->numsite[0];
	    printf("XXmap[%d]=%p:mapid=%d,id=%lld,len=%0.4f,M=%d\n",Xid,Xmap,Xmap->mapid,Xmap->id,Xmap->site[0][M+1],M);
	  }
	}
      } else
	printf("map[%d..%d] duplicated as XXmap[0..%d] and YYmap[0..%d]:maxN=maxM=%lld\n",
	       0,nummaps-1,numX-1,numY-1,maxN);
      fflush(stdout);
    }
  } else { /* split maps into separate arrays XXmap[0..numX-1] and YYmap[0..numY-1] */
    if(DEBUG) assert(Cfirst > 0);

    if((hash_filename && !HashGen) || (PAIRMERGE_SORTED && PairMerge && pairmergeIter <= 1)){    /* sort maps in each mapset by id in ascending order */
      qsort(Gmap,Cfirst,sizeof(Cmap *),(intcmp*)CmapIdInc);
      qsort(&Gmap[Cfirst],nummaps-Cfirst,sizeof(Cmap *),(intcmp*)CmapIdInc);

      if(PairMerge && Gmap[Cfirst-1]->id >= Gmap[Cfirst]->id){/* in absense of hash_filename, this can happen, but check for duplicate id values in two map sets */
	Cid2mapid *id2mapid = new Cid2mapid[Cfirst];
	for(int i = Cfirst; --i >= 0;){
	  Cid2mapid *p = &id2mapid[i];
	  p->id = Gmap[i]->id;
	  p->mapid = i;
        }

	for(int i = Cfirst; i < nummaps; i++){
	  long long id2 = Gmap[i]->id;
	  int index1 = findid(id2,id2mapid,Cfirst);
	  if(index1 >= 0){
	    printf("-pairmerge with -first requires all ids in 1st set to be distinct from all ids in 2nd set:\n");
	    printf("  1st set includes id=%lld from file %s\n",Gmap[index1]->id, vfixx_filename[Gmap[index1]->fileid]);
	    printf("  2nd set includes id=%lld from file %s\n",Gmap[i]->id, vfixx_filename[Gmap[i]->fileid]);
	    fflush(stdout);exit(1);
	  }
	}
	delete [] id2mapid;
      }
    }

    numY = numrefmaps = Cfirst;
    //    if(CfirstPairs == 2)
    //      numY = numrefmaps = nummaps;
    if(PairSplit){/* must use separate array, so both arrays can grow */
      if(VERB>=2){
	printf("Calling maxmapalloc(%d,%d,refmap=%p,0)\n",Cfirst,maxrefmaps,refmap);
	fflush(stdout);
      }
      maxmapalloc(numY,maxrefmaps,refmap,0,1);
      YYmap = refmap;
      for(int i = 0; i < 0 + numY; i++){
	refmap[i-0] = Gmap[i];
	Gmap[i] = 0;
      }
    } else
      YYmap= &Gmap[0];

    if(CfirstPairs == 1){
      XXmap = &Gmap[0];
      numX = nummaps;
    } else {
      XXmap = &Gmap[Cfirst];
      numX = nummaps - Cfirst;
    }

    /* compute seperate maxM, maxN */
    maxM = maxN = 0;
    for(int Xid = 0; Xid < numX; Xid++)
      if(XXmap[Xid]->numsite[0] > maxM)
	maxM = XXmap[Xid]->numsite[0];
    for(int Yid = 0; Yid < numY; Yid++)
      if(YYmap[Yid]->numsite[0] > maxN)
	maxN = YYmap[Yid]->numsite[0];
    if(VERB){
      printf("Gmap[0..%d] renamed as YYmap[0..%d] and Gmap[%d..%d] renamed as XXmap[0..%d](maxN=%lld,maxM=%lld)\n",
	     numY-1,numY-1,(CfirstPairs == 1) ? 0 : Cfirst,nummaps-1,numX-1,maxN,maxM);
      printf("\t Gmap[0]=%p,YYmap[0]=%p,Gmap[%d]=%p,XXmap[0]=%p, YYmap=%p,XXmap=%p\n",
	     Gmap[0],YYmap[0],(CfirstPairs==1) ? 0 : Cfirst,Gmap[(CfirstPairs==1) ? 0 : Cfirst],XXmap[0],YYmap,XXmap);
      if(VERB>=2 && pairmergeIter > 1){
	for(int Yid = 0; Yid < numY;Yid++)
	  printf("YYmap[%d]:mapid=%d,id=%lld\n",Yid,YYmap[Yid]->mapid,YYmap[Yid]->id);
	for(int Xid = 0; Xid < numX;Xid++)
	  printf("XXmap[%d]:mapid=%d,id=%lld\n",Xid,XXmap[Xid]->mapid,XXmap[Xid]->id);
      }
      fflush(stdout);
    }
  }

  for(int Xid = 0; Xid < numX; Xid++){
    XXmap[Xid]->mapid = Xid;
    if(!PairMerge){
      XXmap[Xid]->centercloned = 0;
      XXmap[Xid]->paired = 0;
    }
    XXmap[Xid]->origmap = 0;
  }
  if(Cfirst && CfirstPairs <= 0){
    for(int Yid = 0; Yid < numY; Yid++){
      YYmap[Yid]->mapid = Yid;
      if(!PairMerge){
	YYmap[Yid]->centercloned = 0;
	YYmap[Yid]->paired = 0;
      }
      YYmap[Yid]->origmap = 0;
    }
  } else if(DEBUG/* HERE >=2 */){
    for(int Yid = 0; Yid < numY; Yid++){
      assert(YYmap[Yid]->mapid == Yid);
      if(!PairMerge){
	assert(YYmap[Yid]->centercloned == 0);
	assert(YYmap[Yid]->paired == 0);
      }
      assert(YYmap[Yid]->origmap == 0);
    }
  }

  if(VERB>=2 && pairmergeIter){/* count number of maps with paired==1 */
    int cnt = 0, merged= 0;
    for(int i = 0; i <= nummaps-1; i++){
      if(Gmap[i]->centercloned)
	merged++;
      if(Gmap[i]->paired)
	cnt++;
    }
    printf("pairmerge Iteration %d: %d maps, including %d merged maps and %d newly merged maps\n", pairmergeIter, nummaps-1-0+1, merged, cnt);
    fflush(stdout);
  }

  int orignumY = numY;
  int orignumX = numX;

  if(PairsplitSeparateQueries){
    /* NOTE : With multiple Query maps, if we want to handle each Query map separately, we need to :
       1. Initialize each Query map's Xid field to equal its mapid.
       2. Make one copy of each original Yid for each original Xid, and initialize the Xid field with the Xid's mapid.
       3. Propagate Xid of Query or Reference map to map splits.
       4. Skip alignments between any map pairs that do not have the same Xid field.
    */
    if(DEBUG) assert(PairSplit);
    if(DEBUG) assert(YYmap == refmap);

    /* Initialize each Query map's Xid field to equal its mapid */
    for(int Xid = 0; Xid < numX; Xid++)
      XXmap[Xid]->Xid = XXmap[Xid]->mapid;

    /* Make one copy of each original Yid for each original Xid AND initialize the Xid field */
    maxmapalloc(numY*numX,maxrefmaps,refmap,0,1);
    YYmap = refmap;
    for(int Yid = 0; Yid < numY; Yid++)
      YYmap[Yid]->Xid = 0;
    for(int Xid = 1; Xid < numX; Xid++){
      for(int Yid = 0; Yid < numY;Yid++){
	int nYid = numY*Xid + Yid;
	if(nYid < numY){
	  printf("numX=%d,numY=%d,Yid=%d: numX*numY+Yid=%d overflowed 32 bit signed int\n",numX,numY,Yid,nYid);
	  fflush(stdout);exit(1);
	}
	if(DEBUG) assert(YYmap[Yid]->origmap == NULL);
	YYmap[nYid][0] = YYmap[Yid][0];/* just copy the struct since only the scalar fields will differ */
	YYmap[nYid]->Xid = Xid;
	YYmap[nYid]->mapid = nYid;
	YYmap[nYid]->origmap = YYmap[Yid];
	if(DEBUG) assert(YYmap[nYid]->origmap->origmap == NULL);
	YYmap[nYid]->left[0] = 0;
	YYmap[nYid]->right[0] = YYmap[Yid]->numsite[0]+1;
      }
    }
    numrefmaps = numY *= numX;
  }
  int orignumYnumX = numY;

  int LocalCnt = 0;/* count of local alignments (if PairSplit) */

  //  int binM1 = 0;
  //  int binM1N = -1;

  if(MaxMem > 0){
    MaxMemSiz = (double)MaxMem * 1024.0 * 1024.0 * 1024.0;
    if(RA_MIN_TIME)
      MaxMemSiz2 = MaxMemSiz * 0.9;
  }
  vmemTotsiz = 0;

  #ifdef _OPENMP
  if(!RA_MIN_TIME && MaxMem > 0 && maxM > 0){/* NOTE : maxM,maxN is the largest Xmap,Ymap size */
    double memperthread = (double)(maxN+1) * (double)(maxM  * sizeof(DP)) + sizeof(AL);
    double maxmem = (double)MaxMem * 1024.0 * 1024.0 * 1024.0;
    int limitthreads = max(1,(int)floor(maxmem/memperthread + 0.5));
    if(limitthreads < maxthreads){
      if(VERB){
	printf("max number of threads reduced from %d to %d due to MaxMem=%0.1f Gbytes, per thread memory = %0.3f Gbytes:maxN=%lld,maxM=%lld\n",
	       maxthreads,limitthreads, MaxMem, memperthread/(1024.0*1024.0*1024.0),maxN,maxM);
	fflush(stdout);
      }
      if(DEBUG) assert(limitthreads >= 1);
      maxthreads = limitthreads;
    } else if(VERB){
      printf("max number of threads not reduced from %d due to MaxMem=%0.1f Gbytes, per thread memory = %0.3f Gbytes:maxN=%lld,maxM=%lld\n",
	     maxthreads, MaxMem, memperthread/(1024.0*1024.0*1024.0),maxN,maxM);
      fflush(stdout);
    }
  }
  #endif // _OPENMP

  if((MaxPalindrome > 0.0 || MaxInvPalindrome > 0.0) && PairMerge <= 0){
    printf("ERROR: -MaxPalindrome or -MaxInvPalindrome only supported with -PairMerge\n");
    fflush(stdout);exit(1);
  }

  if(DEBUG && !HashGen && (MaxPalindrome > 0.0 || MaxInvPalindrome > 0.0) && PairMerge && !Cfirst && PairMergeHmap <= 0) 
    assert(!(PairMergeRepeat && pairmergeIter > 1 && !SplitSegDup)); 

  int numPal = ((MaxPalindrome > 0.0 || MaxInvPalindrome > 0.0) && PairMerge && !Cfirst && (pairmergeIter <= 1 || SplitSegDup) && PairMergeHmap <= 0) ? numY : 0;

  double dMaxPair = min(max((numY + 0.01)*(numX + 0.01),nummaps*(nummaps-1.0)*0.5) + numPal,
			max(nummaps * 2.0 + numPal, maxthreads*((FirstAlignments > 0) ? 1024.0 : 1024.0*64.0)));
  if(dMaxPair > MASK(31)){
    printf("ERROR: MaxPair 31-bit overflow: dMaxPair = %0.1f cannot be converted to signed int (31-bit): FirstAlignments=%d,nummaps=%d,numY=%d,numX=%d,numPal=%d\n", 
	   dMaxPair,FirstAlignments,nummaps,numY,numX,numPal);
    fflush(stdout);
    assert(dMaxPair <= MASK(31));
    exit(1);
  }
  assert(dMaxPair >= 0.0);
  MaxPair = floor(dMaxPair);
  MaxPair = max(1, MaxPair);

  //  MaxPair /= 16;
  if(VERB/* HERE >=2 */){
    printf("maxthreads=%d,MaxPair=%d,FirstAlignments=%d,nummaps=%d,numY=%d,numX=%d,numPal=%d\n",maxthreads,MaxPair,FirstAlignments,nummaps,numY,numX,numPal);
    fflush(stdout);
  }
  if(DEBUG) assert(MaxPair > 0);

  BinaryFormat = 0;

  if(hash_filename){/* read in hash pairs */
#if HASH_STREAM < 2
  printf("HASH_STREAM=%d no longer supported (see hash.h)\n",HASH_STREAM);
  fflush(stdout);exit(1);
#endif
    if(!HashGen){
      printf("Previously computed hash table no longer supported : use -hashgen\n");
      fflush(stdout);exit(1);
    }

    if(DEBUG) assert(hashpairs1 == 0 && numhashpairs1 == 0);

    size_t L = strlen(hash_filename);
    if(L >= strlen(".hashbin") && !strcmp(&hash_filename[L-strlen(".hashbin")],".hashbin")){/* binary file input */
      BinaryFormat = 1;/* used later to enable streaming file input */

      /* HASH_STREAM : set up streamed IO */
      if(RepeatShift > 0.0){
	printf("-hash not supported with -Repeat\n");
	fflush(stdout);exit(1);
      }
      if(PairSplit){
	printf("-hash not yet supported with -pairsplit\n");
	fflush(stdout);exit(1);
      }

      /* set up input buffers */
      size_t num = MaxPair;
      numhashpairs1 = maxhashpairs1 = 0;
      maxhashalloc(num,numhashpairs1,maxhashpairs1,hashpairs1);

      numhashpairs2 = maxhashpairs2 = 0;
      maxhashalloc(num,numhashpairs2,maxhashpairs2,hashpairs2);

      MatchMax = hash_open(HashFP, hash_filename);

    } else {/* Old text format no longer supported */
      printf("Non-binary Hashtable format no longer supported for %s\n",hash_filename);
      fflush(stdout);exit(1);
    }

    hashpairs1[numhashpairs1].id1 = 0;
  }
  nexthash1 = hashpairs1;
  nexthash2 = hashpairs2;
  if(hashsplit){
    printf("-hashsplit no longer supported\n");
    fflush(stdout);exit(1);
  }
  if(VERB>=2){
    printf("hashpairs1=%p,RepeatShift= %0.3f, PairSplit=%d\n",hashpairs1,RepeatShift, PairSplit);
    fflush(stdout);
  }

  size_t numaligns_compressed = 0;/* last alignment[] index already checked for NULL entries (if PairSplit) */
  size_t prev_aligncnt = 0;/* accumulated count of alignments above threshold and chimpair == 0 in alignment[0..numaligns_compress-1] (if PairSplit) */
  size_t numaligns_hashed = 0;/* last alignment[] index already added to PairHash table (if PairSplit) */

  int wcnt = 0;
  while(1){/* With -pairsplit : repeatedly locate best local alignment amongst maps (excluding pairs already confirmed with chimpair = -1, or below threshold (or discarded) with chimpair = 1) */
    wcnt++;

    //    LVERB = (wcnt >= 2 ) ? 1 : 0;

    if(PairSplit){/* update MaxPair, incase numX,numY changed */
      MaxPair = min(max(numY*numX*1.0,nummaps*(nummaps-1.0)*0.5), max(nummaps*2.0, maxthreads*((FirstAlignments > 0) ? 1024.0 : 1024.0*64.0)));
      if(VERB>=2){
	printf("maxthreads=%d,MaxPair=%d,FirstAlignments=%d,nummaps=%d\n",maxthreads,MaxPair,FirstAlignments,nummaps);
	fflush(stdout);
      }
      if(DEBUG) assert(MaxPair > 0);
    }

    if(LVERB){
      int badcnt = 0;
      for(size_t i = 0; i < numaligns; i++)
	if(alignment[i] && alignment[i]->chimpair == 1)
	  badcnt++;
      printf("Found %llu local alignments (%d confirmed,%d ruled out): checking for additional alignments amongst numY=%d,numX=%d maps\n",
	     (unsigned long long)numaligns,LocalCnt,badcnt,numY,numX);
      fflush(stdout);
    }

    YidList1 = new int[MaxPair];
    XidList1 = new int[MaxPair];
    phashList1 = new Chashpair*[MaxPair];
    NumPair1 = 0;
  
    YidList2 = new int[MaxPair];
    XidList2 = new int[MaxPair];
    phashList2 = new Chashpair*[MaxPair];
    NumPair2 = 0;

    size_t PairCnt = 0, LastPairCnt = 0;
    int LastRepeat = 0;

    size_t numaligns_start = 0;

    if(numPal > 0){/* add Yid,Yid to start of list */    
      if(DEBUG) assert(MaxPalindrome > 0.0 || MaxInvPalindrome > 0.0);
      for(int Yid = 0; Yid < numY; Yid++){
	if(DEBUG) assert(NumPair1 < MaxPair);
	phashList1[NumPair1] = 0;
	YidList1[NumPair1] = Yid;
	XidList1[NumPair1++] = Yid;
      }
    }

    for(int Yid = 0;Yid < numY; Yid++){
      if(LVERB && ((Cfirst && CfirstPairs != 1) ? (CfirstPairs==2 ? max(Cfirst,Yid+1) : 0) : Yid+1) <= numX-1){
	printf("aligning YYmap[%d/%d](id=%lld,centercloned=%d) with XXmap[%d..%d]:Cfirst=%d,%d\n", Yid, numY, YYmap[Yid]->id,YYmap[Yid]->centercloned, Cfirst ? 0 : Yid+1, numX-1,Cfirst,CfirstPairs);
	fflush(stdout);
      }
      if(DEBUG>=2) assert(YYmap[Yid]->id >= 0);
      if(VERB>=2){
	printf("Yid=%d(id=%lld),numY=%d\n",Yid,YYmap[Yid]->id,numY);
	fflush(stdout);
      }

      int Yunpaired,YYid;
      if(DEBUG>=2){
	Yunpaired = (Cfirst && pairmergeIter > 1 && !YYmap[Yid]->paired) ? 1 : 0;
	YYid = YYmap[Yid]->id;
      }

      if(FirstAlignments > 0 && mapcnt >= FirstAlignments)
	break;

      if(RepeatShift <= 0.0 && !nexthash1 && !PairSplit){ /* fast special case (without hashtable) */
	int Xid = (Cfirst && CfirstPairs != 1) ? (CfirstPairs==2 ? max(Cfirst,Yid+1) : 0) : Yid+1;
	
	// Tried multi-threading the following loop but it slows down
	for(; Xid < numX; Xid++){
	  if(LVERB>=2){
	    printf("aligning YYmap[%d/%d](centercloned=%d) with XXmap[%d/%d](centercloned=%d)\n", Yid, numY, YYmap[Yid]->centercloned, Xid, numX, XXmap[Xid]->centercloned);
	    fflush(stdout);
	  }

	  if(DEBUG>=2 && (Yunpaired && !XXmap[Xid]->paired)){
	    printf("aligning YYmap[%d/%d](paired=%d) with XXmap[%d/%d](paired=%d),Yunpaired=%d\n",Yid,numY,YYmap[Yid]->paired,Xid,numX,XXmap[Xid]->paired,Yunpaired);
	    fflush(stdout);
	    assert(!(Yunpaired && !XXmap[Xid]->paired));
	  }
	  if(DEBUG>=2) assert(YYid != XXmap[Xid]->id);

	  phashList1[NumPair1] = 0;
	  YidList1[NumPair1] = Yid;
	  XidList1[NumPair1++] = Xid;

	  if(VERB>=2){
	    printf("Yid=%d,Xid=%d:NumPair=%d,MaxPair=%d,numX=%d,numY=%d\n",Yid,Xid,NumPair1,MaxPair,numX,numY);
	    fflush(stdout);
	  }
	}
      } else if(HASH_STREAM && HashFP /* && RepeatShift <= 0.0 && !PairSplit*/){/* fast special case (with streamed hashtable) */
	if(DEBUG) assert(nexthash1);
	if(DEBUG) assert(BinaryFormat==1);
	if(DEBUG && maxhashpairs1 != (size_t)MaxPair){
	  printf("maxhashpairs1=%llu, MaxPair=%d\n",(unsigned long long)maxhashpairs1,MaxPair);
	  fflush(stdout);
	  assert(maxhashpairs1 == (size_t)MaxPair);
	}

	if(VERB>=2){
	  printf("NumPair=%d,MaxPair=%d:nexthash = &hashpairs[%llu],numhashpairs=%llu,maxhashpairs=%llu,totalhashpairs=%llu\n",
		 NumPair1,MaxPair,(unsigned long long)(nexthash1-hashpairs1),(unsigned long long)numhashpairs1,(unsigned long long)maxhashpairs1,(unsigned long long)totalhashpairs);
	  fflush(stdout);
	}

	/* first check if we need to read in more hashtable entries */
	if(MaxPair - NumPair1 >= nummaps)
	  loadhash(nexthash1, hashpairs1, numhashpairs1, maxhashpairs1, YidList1, XidList1, phashList1, NumPair1);

	if(hash_eof(HashFP))
	  Yid = numY - 1;
      } else {/* general case without hashtable */
	Cmap *Ymap = YYmap[Yid];
	    
	if(PairSplit){
	  
	  if(Ymap->centercloned)
	    continue;

	  /* insert previously computed alignments into alignment hashtable */
	  /* only insert alignments not previously inserted : alignment[numaligns_hashed..numaligns-1] */
	  for(size_t i = numaligns_hashed; i < numaligns; i++){
	    Calign *palign = alignment[i];
	    if(palign)
	      PairHash_insert(palign->mapid1,palign->mapid2);
	  }
	  numaligns_hashed = numaligns;
	}

	int Xid = (RepeatShift > 0.0) ? Yid : (Cfirst && CfirstPairs !=1) ? (CfirstPairs==2 ? max(Cfirst,Yid+1) : 0) : Yid+1;
	
	for(; (RepeatShift > 0.0) ? (Xid <= Yid) : (Xid < numX); Xid++){
	  if(LVERB>=2){
	    printf("aligning YYmap[%d/%d](id=%lld,centercloned=%d) with XXmap[%d/%d](id=%lld,centercloned=%d)\n", 
		   Yid, numY, YYmap[Yid]->id,YYmap[Yid]->centercloned, Xid, numX, XXmap[Xid]->id,XXmap[Xid]->centercloned);
	    fflush(stdout);
	  }

	  Cmap *Xmap = XXmap[Xid];
      
	  if(PairSplit){
	    if(Xmap->centercloned)
	      continue;

	    if(PairsplitSeparateQueries && Xmap->Xid != Ymap->Xid)
	      continue;

	    if(PairHash_find(Yid,Xid))
	      continue;
	  }

	  phashList1[NumPair1] = 0;
	  if(DEBUG>=2 && RepeatShift > 0.0) assert(Xid==Yid);
	  YidList1[NumPair1] = Yid;
	  XidList1[NumPair1++] = Xid;
	  if(VERB>=2){
	    printf("Yid=%d,Xid=%d:NumPair=%d,MaxPair=%d\n",Yid,Xid,NumPair1,MaxPair);
	    fflush(stdout);
	  }
	}
      } /* General case generating YidList[],XidList[],phashList[] */

      if(!(HASH_STREAM>=2 && HashFP) && NumPair1 < MaxPair && (RepeatShift > 0.0 || MaxPair - NumPair1 >= nummaps) && Yid+1 < numY)
      	continue;/* accumulate more Yid,Xid pairs with next Yid */

      if(NumPair1 <= 0)
	continue;

      if(VERB /* && !PairSplit */ && ((HASH_STREAM && HashFP) || PairCnt + NumPair1 > LastPairCnt + maxthreads)){
	if(RepeatShift <= 0.0){
	  int Yid_start = YidList1[0];
	  int Yid_end = YidList1[NumPair1-1];
	  int Xid_start = (Cfirst && CfirstPairs !=1) ? (CfirstPairs==2 ? max(Cfirst,Yid_start+1) : 0) : Yid_start+1;

	  double completion = Cfirst ? ((double)min(Yid_start,Yid_end))/numY : 1.0 - (double)(numY-min(Yid_start,Yid_end))*(numY-min(Yid_start,Yid_end)-1.0)/(numY*(numY-1.0));
	  printf("%d:aligning YYmap[%d..%d] with XXmap[%d..%d], NumPair=%d/%d, Cfirst=%d,%d (previous alignments found=%lld/%llu),Yid=%d/%d cum CPU time=%0.6f, wall time=%0.6f %4.1f%% complete\n", 
		 PairMergeRepeat ? pairmergeIter : wcnt,min(Yid_start,Yid_end),max(Yid_start, Yid_end), Xid_start, numX-1, NumPair1, MaxPair, Cfirst, CfirstPairs, mapcnt, (unsigned long long)PairCnt, 
		 Yid,numY,mtime(), wtime(), completion*100.0);
	}
	if(RepeatShift > 0.0)
	  printf("Checking repeats in XXmap[%d..%d] (total repeats found= %lld)\n", LastRepeat, XidList1[NumPair1-1], mapcnt);
	fflush(stdout);


	LastPairCnt = NumPair1 + PairCnt;
	LastRepeat = XidList1[NumPair1-1]+1;
      }
      PairCnt += NumPair1;

      if(!PairSplit){  /* compress alignment[numaligns_start .. numaligns-1] */
	size_t j = numaligns_start;
	for(size_t k = j; k < numaligns;k++){
	  if(!alignment[k])
	    continue;
	  if(j < k){
	    if(DEBUG>=2) assert(alignment[j]==0);
	    alignment[j] = alignment[k];
	    alignment[k] = 0;
	  }
	  j++;
	}
	if(VERB>=2){
	  printf("Yid=%d/%d:reduced numaligns from %llu to %llu, NumPair=%d\n",
		 Yid,numY,(unsigned long long)numaligns,(unsigned long long)j,NumPair1);
	  fflush(stdout);
	}
	numaligns = j;
      }
      numaligns_start = numaligns;/* no need to compress alignment[0..numalign_start-1] */

      /* preallocate space for alignment pointers, to minimize allocation in multi-threaded region */
      maxalignallocNULL(numaligns+NumPair1,alignment,numaligns,maxaligns);

      int numthreads = 1;
      #ifdef _OPENMP
      numthreads = maxthreads;
      if(numthreads > NumPair1)
	numthreads = max(1,NumPair1);
      #endif

      int block = (PairSplit || PairMerge || RA_MIN_TIME) ? 1 : min(16,max(1,NumPair1/(64*numthreads)));

      if(VERB /* && numthreads > 1 */ && !threads_displayed) {
	printf("Yid=%d:nummaps=%d:Using %u threads,block size=%d(NumPair=%d,MaxPair=%d)\n",Yid,nummaps,numthreads,block,NumPair1,MaxPair);
	fflush(stdout);
	threads_displayed = 1;
      }

      //      int origYid = Yid;
      
      try {/* following code may throw out-of-memory exception */
	if(A_numthreads < numthreads || A_maxN < maxN || A_maxN < maxM) {
	  if(VERB>=2){
	    printf("Allocating A_numthreads:maxN=%lld,maxM=%lld,A_maxN=%lldnumthreads=%d\n",maxN,maxM,A_maxN,numthreads);
	    fflush(stdout);
	  }

	  for(int i=0; i < A_numthreads; i++) {
	    delete [] A_threads[i].score;
	    delete [] A_threads[i].Uscore;
	    delete [] A_threads[i].G;
	    delete [] A_threads[i].H;
	    delete [] A_threads[i].Rijy;
	    delete [] A_threads[i].Lijy;
	    delete [] A_threads[i].Rijx;
	    delete [] A_threads[i].Lijx;

	    delete [] A_threads[i].JMIN;
	    delete [] A_threads[i].JMAX;
	    delete [] A_threads[i].Imin;
	    delete [] A_threads[i].Imax;
	  }
	  if(A_threads) delete [] A_threads;
		
	  if(numthreads > A_numthreads) A_numthreads = numthreads;
	  if(maxN > A_maxN) A_maxN = maxN;
	  if(maxM > A_maxN) A_maxN = maxM;
		
	  A_threads = new AL[A_numthreads];
	  for(int i= 0; i < A_numthreads; i++) {
	    A_threads[i].score = new PFLOAT*[A_maxN+1];
	    A_threads[i].Uscore = new PFLOAT*[A_maxN+1];
	    A_threads[i].G = new int*[A_maxN+1];
	    A_threads[i].H = new int*[A_maxN+1];
	    A_threads[i].Rijy = new int*[A_maxN+1];
	    A_threads[i].Lijy = new int*[A_maxN+1];
	    A_threads[i].Rijx = new int*[A_maxN+1];
	    A_threads[i].Lijx = new int*[A_maxN+1];

	    A_threads[i].JMIN = new int[A_maxN+1];
	    A_threads[i].JMAX = new int[A_maxN+1];
	    A_threads[i].Imin = new int[A_maxN+1];
	    A_threads[i].Imax = new int[A_maxN+1];
	  }
	}
      
	if(mem_numthreads < numthreads || mem_maxN < maxN || mem_maxM < maxM) {
	  if(VERB>=2){
	    printf("Allocating mem_thread_arrays: mem_maxN=%lld mem_maxM=%lld, mem_numthreads=%d\n", 
		   mem_maxN, mem_maxM, mem_numthreads);
	    fflush(stdout);
	  }
	  for(int i=0; i < mem_numthreads; i++) {
	    if(RA_MIN_TIME){/* use munmap to free Fmem/Imem */
	      if(ThreadMem[i].Fmem && munmap(ThreadMem[i].Fmem, ThreadMem[i].MemSiz)){
		int eno = errno;
		char *err = strerror(eno);
		printf("munmap(%p,%lld) failed: errno=%d:%s\n",ThreadMem[i].Fmem, ThreadMem[i].MemSiz, eno, err);
		dumpmemmap();
		fflush(stdout);exit(1);
	      }
	      ThreadMem[i].Fmem = NULL;
	      ThreadMem[i].Imem = NULL;
	      ThreadMem[i].MemSiz = 0;
	    } else {
	      delete [] ThreadMem[i].Fmem;
	      delete [] ThreadMem[i].Imem;
	    }
	  }

	  if(numthreads > mem_numthreads) mem_numthreads = numthreads;
	  if(maxN > mem_maxN) mem_maxN = maxN;
	  if(maxM > mem_maxM) mem_maxM = maxM;
	     
	  if(DEBUG) assert(sizeof(CThreadMem) == 64 * 3);

	  ThreadMem = new CThreadMem[A_numthreads];

	  double wt = wtime();

	  for(int i=0; i < A_numthreads; i++) {
	    long long Fsiz = mem_maxN * ((mem_maxM + 64LL) + 16) * DPsizF + 16;
	    if(DEBUG) assert(Fsiz >= 0);
	    long long Isiz = mem_maxN * ((mem_maxM + 64LL) + 16) * DPsizI + 16;
	    if(DEBUG) assert(Isiz >= 0);

	    if(RA_MIN_TIME){/* on demand mmap or madvise based memory allocation */
	      if(RA_MIN_MEM >= 8){/* on demand mmap based memory allocation */
		//		ThreadMem[i].Fmem = NULL;
		//		ThreadMem[i].Imem = NULL;
		//		ThreadMem[i].MemSiz = 0;
	      } else {/* one time mmap and on demand madvise based memory allocation */
		long long MemSiz  = ((Fsiz * sizeof(PFLOAT) + Isiz * sizeof(int)) + (PAGE-1)) & ~(PAGE-1);
		
		char *newmem = (char *)mmap(NULL, MemSiz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if(newmem == MAP_FAILED || newmem == NULL){
		  int eno = errno;
		  char *err = strerror(eno);
		  
                  #pragma omp critical
		  {
		    printf("pairalign: mmap of %lld bytes failed:errno=%d:%s\n", MemSiz,eno,err);
		    dumpmemmap();
		    fflush(stdout);exit(1);
		  }
		}
		
		ThreadMem[i].Fmem = (PFLOAT *) newmem;
		ThreadMem[i].Imem = (int *) &newmem[Fsiz * sizeof(PFLOAT)];
		ThreadMem[i].MemSiz = MemSiz;
		ThreadMem[i].MemUsed = 0;
		ThreadMem[i].lasttime = wt;
		if(TIME_SPREAD){
		  double Inv_numthreads = (double)i/A_numthreads;
                  ThreadMem[i].lasttime -= RA_MIN_TIME * Inv_numthreads;
		  ThreadMem[i].lastGetStatm -= VMRSS_CHECK * Inv_numthreads;
		  ThreadMem[i].lastGetMem -= VMSWAP_CHECK * Inv_numthreads;
                }
		if(VERB>=2){
		  printf("Allocated using mmap for thread %d: Fmem=%p,Imem=%p,MemSiz=%lld (RA_MIN_MEM=%d,RA_MIN_TIME=%d)\n",
			 i,ThreadMem[i].Fmem,ThreadMem[i].Imem,ThreadMem[i].MemSiz,RA_MIN_MEM,RA_MIN_TIME);
		  fflush(stdout);
		}
	      }
	    } else {// static (one time) memory allocation
	      ThreadMem[i].Fmem = new PFLOAT[Fsiz];
	      ThreadMem[i].Imem = new int[Isiz];
	      ThreadMem[i].MemUsed = ThreadMem[i].MemSiz = Fsiz * sizeof(PFLOAT) + Isiz * sizeof(int);
	      if(VERB>=2){
		printf("tid=%d:Fmem=%p (Fsiz= %lld), Imem=%p (Isiz= %lld)\n",i, ThreadMem[i].Fmem, Fsiz, ThreadMem[i].Imem, Isiz);
		fflush(stdout);
	      }
	    }
	  } /* for i = 0 .. A_numthreads -1 */

	  long long VmSize,VmRSS,VmSwap,VmHWM;
	  getmem(VmSize,VmRSS,VmSwap,&VmHWM);	  
	  vmemTotsiz = 0;
	  //	  vmemTotsiz = VmRSS + VmSwap;
	  if(VERB/* HERE HERE >=2 */){
            printf("Allocation complete: maxN=%lld maxM=%lld:vmemTotsiz= %0.4f, VmSize= %0.4f, VmRSS= %0.4f(VmHWM= %0.4f), VmSwap= %0.4f\n", 
		   mem_maxN, mem_maxM, vmemTotsiz*1e-9, VmSize*1e-9, VmRSS*1e-9,VmHWM*1e-9,VmSwap*1e-9);
	    fflush(stdout);
	  }
	}
      } catch (exception& e){
	cout << e.what() << endl;
	printf("Memory Allocation failed:mem_maxN=%lld mem_maxM=%lld,A_numthreads=%d\n",
	       mem_maxN,mem_maxM,A_numthreads);
	fflush(stdout);exit(1);
      }

      if(DEBUG) assert(NumPair1 <= MaxPair);

      if(VERB>=2 && numthreads > 1){
	printf("Entering parallel section\n");
	fflush(stdout);
      }
  
      FLOAT varSF2 = 2.0*SF[0]*SF[0];
      FLOAT varSD = fabs(SD[0]) * SD[0];
      FLOAT varSR = SR[0]*SR[0];

      if(VERB/* HERE >=2 */){
	printf("NumPair=%d,dynamic block=%d,numthreads= %d: alignments so far=%lld,maxN=%lld,maxM=%lld\n",NumPair1,block,numthreads,mapcnt,maxN,maxM);
	fflush(stdout);
      }

      NumThreads = numthreads;
      PausedThreads = 0;
      VmRSSreserved = 0;

      int LoopCnt = 0;

      #pragma omp parallel num_threads(numthreads) if (numthreads > 1)
      //      #pragma omp parallel num_threads(numthreads) proc_bind(spread) if (numthreads > 1)
      {
	int myLoopCnt = 0;

	int tid = 0;
        #ifdef _OPENMP
	tid = omp_get_thread_num ();
        #endif

	/* allocate thread local DP array A[I=1..maxN][J=1..maxM] */
	AL A = A_threads[tid];
	PFLOAT * &Fmem = ThreadMem[tid].Fmem;
	int * &Imem = ThreadMem[tid].Imem;
	CThreadMem *Mem = &ThreadMem[tid];

	long long StrideMax = 0;/* initial allocation (elements per array) in current thread */

	/* thread local variables */
	long long Tmapcnt = 0, Tpaircnt=0;
	long long Ttcnt = 0;// thread local versions of mapcnt,paircnt and tcnt
	long long Tposcnt = 0, Tsitecnt = 0, Tmisscnt = 0, Toutliertot = 0, Toutliercnt = 0, Tsegcnt = 0;
	double Tlen = 0.0, TlenB = 0.0, Toutliersum = 0.0, Terrnorm = 0.0;

	// NOTE : the following loop uses nowait, so memory is deallocated by threads that have completed, so memory starvation cannot occur with RA_MIN_TIME

	#pragma omp for nowait schedule(dynamic,block)
	for(int PairId = 0; PairId < NumPair1; PairId++){
	  if(OMP_DEBUG) myLoopCnt++;
	  if(HASH_STREAM>=2 && !tid && HashFP /* && RepeatShift <= 0.0 && !PairSplit */ && NumPair2 == 0){
	    /* Read ahead from disk into a 2nd buffer in the master thread. In Non-parallel section swap the 2 buffers instead of reading from disk */
	    loadhash(nexthash2,hashpairs2,numhashpairs2,maxhashpairs2,YidList2,XidList2,phashList2,NumPair2);
	  }

	  Chashpair *phash = phashList1[PairId];
	  int Yid = YidList1[PairId];
	  int Xid = XidList1[PairId];
	  if(DEBUG && (Yid<0 || Xid<0)) {
	    printf("INTERNAL ERROR: negative id encountered: Xid=%d Yid=%d PairId=%d,NumPair=%d\n", Yid, Xid, PairId,NumPair1);
	    fflush(stdout);
	    assert(Yid >= 0 && Xid >= 0);
	  }
	  //	  if(!EVERB) continue;
	  Cmap *Ymap = YYmap[Yid];
	  Cmap *Xmap = XXmap[Xid];
	  int alignid = numaligns_start + PairId; 
	  Calign * &align = alignment[alignid];

	  if(VERB>=2){
	    printf("PairId=%d/%d: Yid=%d(id=%lld),Xid=%d(id=%lld),N=%d,M=%d,tid=%d:Starting at wtime=%0.6f\n",PairId,NumPair1,Yid,Ymap->id,Xid,Xmap->id,Ymap->numsite[0],Xmap->numsite[0],tid,wtime());
	    fflush(stdout);
	  }

	  if(LVERB || EVERB>=2){
	    printf("aligning YYmap[%d/%d](centercloned=%d) with XXmap[%d/%d](centercloned=%d)\n", Yid, numY, YYmap[Yid]->centercloned, Xid, numX, XXmap[Xid]->centercloned);
	    fflush(stdout);
	  }

	  if(PairSplit){
	    if(DEBUG>=2) assert(PairHash);
	    if(!align)
	      align = new Calign[colors==1 ? 1 : colors+1];
	  } else {
	    if(DEBUG)assert(align == 0);// NOTE : if this fails, need to initialize align->numpairs = 0, in case pairalign() return without making any changes to align
	  }

	  if(VERB>=2){
            #pragma omp critical
	    {
	      tcnt++;
	      printf("Yid=%d,Xid=%d:startloc=%0.3f,%0.3f,endloc=%0.3f,%0.3f,pairs cnt=%lld\n",
		     Yid,Xid,Ymap->startloc,Xmap->startloc,Ymap->endloc,Xmap->endloc,tcnt);
	      fflush(stdout);
	    }
	  } else 
	    Ttcnt++;

	  if(DEBUG>=2 && (Xmap->id < 0 || Ymap->id < 0)){
	    #pragma omp critical
	    {
	      printf("Calling pairalign:Yid=%d/%d(id=%lld),Xid=%d/%d(id=%lld), PairId=%d/%d, alignid=%d\n",
	         Yid,numY,Ymap->id,Xid,numX,Xmap->id,PairId,NumPair1,alignid);
	      fflush(stdout);
	      assert(Xmap->id >= 0);
	      assert(Ymap->id >= 0);
	    }
	  }

	  if(PAIRSPLIT_EXTEND && PairSplit){
	    align->refstart = 0;
	    align->start = 0;
	    align->refend = Ymap->numsite[0]+1;
	    align->end = Xmap->numsite[0]+1;
	  }

	  if(DEBUG>=2) assert(0 <= Xid && Xid < numX);
	  if(DEBUG>=2) assert(0 <= Yid && Yid < numY);

	  if(VERB && !(PairId % 100000)){
	    printf("PairId=%d/%d:Yid=%d,Xid=%d, tid=%d, wall time= %0.6f secs\n",PairId,NumPair1,Yid,Xid,tid,wtime());
	    fflush(stdout);
	  }

	  pairalign(Ymap, Xmap, align, Yid, Xid, phash, A, Fmem,Imem, Mem,StrideMax,tid);

	  if(DEBUG>=2) assert(align==alignment[alignid]);
	  if(DEBUG>=2 && PairSplit) assert(align && align->mapid1 == Yid && align->mapid2 == Xid);

	  if((DEBUG>=2 || VERB>=2) && align && align->numpairs > 0){
	    int N = Ymap->numsite[0];
	    int M = Xmap->numsite[0];
	    if(!(align->sites1[0] >= 1 && align->sites1[align->numpairs-1] <= N)){
	      printf("pairalign_allpairs:Yid=%d,Xid=%d:align->numpairs=%d,I=%d..%d(N=%d),align->sites1=%p,align=%p,align->logPV=%0.2f,align->repeat=%d,align->IsRepeatRegion()=%d\n",
		     Yid,Xid,align->numpairs,align->sites1[0],align->sites1[align->numpairs-1],N,align->sites1,align,align->logPV,align->repeat,align->IsRepeatRegion());
	      fflush(stdout);
	      assert(align->sites1[0] >= 1 && align->sites1[align->numpairs-1] <= N);
	    }
	    assert(align->sites2[0] >= 1 && align->sites2[align->numpairs-1] <= M);
	  }


	  if(VERB>=2){
	    printf("PairId=%d/%d: Yid=%d(id=%lld),Xid=%d(id=%lld),N=%d,M=%d,tid=%d:Completed at wtime=%0.6f\n",PairId,NumPair1,Yid,Ymap->id,Xid,Xmap->id,Ymap->numsite[0],Xmap->numsite[0],tid,wtime());
	    fflush(stdout);
	  }

	  if(align && align->numpairs > 0 && align->score > MINSCORE){
	    if(DEBUG) assert(align->mapid1 == Yid);
	    if(DEBUG) assert(align->mapid2 == Xid);

	    /* accumulate true positive count (paircnt), if possible */
	    FLOAT *Y = Ymap->site[0];
	    FLOAT *X = Xmap->site[0];
	    int N = Ymap->numsite[0];
	    int M = Xmap->numsite[0];

	    if(DEBUG>=2){
	      int Lijy = max(1,align->Lij1);
	      int Lijx = max(1,align->Lij2);
	      int Rijy= min(N,align->Rij1);
	      int Rijx= min(M,align->Rij2);
	      if(DEBUG) assert(Rijx >= Lijx + align->numpairs - 1);
	      if(DEBUG && !(Rijy >= Lijy + align->numpairs - 1)){
		#pragma omp critical
		{
		  printf("Yid=%d(%lld),Xid=%d(%lld):N=%d,M=%d:alignid=%d:align->numpairs=%d,align->score=%0.4f,align->logPV=%0.2f:align->Rij1=%d,align->align->Lij1=%d,Lijy=%d,Rijy=%d,align->repeat=%d,align->IsRepeatRegion()=%d\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,N,M,alignid,align->numpairs,align->score,align->logPV,align->Rij1,align->Lij1,Lijy,Rijy,align->repeat,align->IsRepeatRegion());
		  fflush(stdout);
		  assert(Rijy >= Lijy + align->numpairs - 1);
		}
	      }
	    }

	    double trueoffset = Xmap->startloc - Ymap->startloc;// should use endloc to compute trueoffset
	    double trueoverlap = (trueoffset>=0.0) ? min(Y[N+1]-trueoffset,X[M+1]) : min(X[M+1]+trueoffset,Y[N+1]);

	    align->trueoffset = (XXmap[0]->startloc == 0.0) ? aNaN : trueoffset;
	    align->trueoverlap = (XXmap[0]->startloc == 0.0) ? aNaN : trueoverlap;

	    align->trueTP = (XXmap[0]->startloc==0.0) ? -1 : 0;
	    if(!(XXmap[0]->startloc==0.0) && trueoverlap >= (X[M+1]+Y[N+1])*(0.5*MINOVERLAP)){
	      align->trueTP = 1;

	      Tpaircnt++;
	    }

	    if(align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= AlignedSiteThreshold && 
 	                  ((PAIRMERGE_REPEATMASK && PairMerge) || !align->IsRepeatRegion())){/* good alignment */
	      if(DEBUG>=2){
		assert(align->numpairs > 0);
		assert(align->sites1[0] >= 1 && align->sites1[align->numpairs-1] <= N);
		assert(align->sites2[0] >= 1 && align->sites2[align->numpairs-1] <= M);
	      }

	      if(VERB>=2){
		#pragma omp atomic
		mapcnt++;
	      } else
		Tmapcnt++;

	      if(VERB>=2){
		#pragma omp critical
		{
		  printf("Yid=%lld,Xid=%lld:number of alignments = %lld\n",YYmap[Yid]->id,XXmap[Xid]->id,Tmapcnt+mapcnt);
		  fflush(stdout);
		}
	      }
	    } else { /* poor quality (below threshold) alignment */
	      if(!PairSplit){/* minimize memory use : free up align since alignment will not be output or used subsequently */
		/* free up memory of align : this will only happen when IsRepeatRegion() discarded the alignment */
		if(DEBUG) assert(align);
		align->allfree();
		delete [] align;
		align = 0; 
		if(DEBUG) assert(alignment[alignid]==0);
	      } else /* mark for later deletion : need to keep until while(1) loop is exited so this pair (Xid,Yid) is not repeately aligned */
		align->chimpair = 1;
	    }

	    if(!align)
	      continue;

	    /* collect some more alignment statistics */
	    int I = align->sites1[0];
	    int J = align->sites2[0];
	    int K = align->numpairs;
	    int Imax = align->sites1[K-1];
	    int Jmax = align->sites2[K-1];
	    int Lijy = max(1,align->Lij1);
	    int Lijx = max(1,align->Lij2);
	    int Rijy= min(N,align->Rij1);
	    int Rijx= min(M,align->Rij2);
	    if(DEBUG) assert(Rijx >= Lijx + align->numpairs - 1);
	    if(DEBUG) assert(Rijy >= Lijy + align->numpairs - 1);

	    if(align->score > 0.0)
	      Tposcnt++;

	    int XYsites = 0;
	    if(align->Lend > -2){
	      double xy = min(Y[I],align->orientation ? X[M+1]-X[M+1-J] : X[J]);
	      Tlen += 2.0 * xy;
	      TlenB += 2.0 * max(0.0, xy - resKB2);
	      XYsites += I-Lijy + J-Lijx;
	    }
	    if(align->Rend > -2){
	      double xy = min(Y[N+1]-Y[Imax], align->orientation ? X[M+1-Jmax] : X[M+1]-X[Jmax]);
	      Tlen += 2.0 * xy;
	      TlenB += 2.0 * max(0.0, xy - resKB2);
	      XYsites += Rijy-Imax + Rijx-Jmax;
	    }
		
	    Tsitecnt += XYsites + 1;
	    Tmisscnt += XYsites;


	    int G = I;
	    int H = J;
	    for(int T=1;T < align->numpairs;T++, G = I, H = J){

	      I = align->sites1[T];
	      J = align->sites2[T];

	      FLOAT y = Y[I] - Y[G];
	      FLOAT x = (align->orientation==0) ? X[J] - X[H] : X[M+1-H] - X[M+1-J];
	      int m = J - G;
	      int n = I - H;

	      if(Poutlier > 0.0){
		Toutliertot++;
		PFLOAT Bias,Pen;
		SintDetail(x,y,m,n,J,I,Bias,Pen);
		double OutPenBias = OutlierPenalty + Bias*biasWToutlier;
		double ExpOutlierPenalty = exp(OutPenBias);
		double val = ExpOutlierPenalty/(ExpOutlierPenalty+exp(Pen+Bias));
		if(DEBUG) assert(isfinite(val));
		Toutliersum += val;
		if(DEBUG) assert(isfinite(Toutliersum));
		if(!OUTLIER_DELTA(x-y) || Pen + Bias < OutPenBias){
		  Toutliercnt++;
		  continue;// don't include this interval in error parameter estimation
		}
	      }

	      Tlen += x+y;
	      TlenB += max(ZERO,x - resKB2) + max(ZERO,y - resKB2);
	      Tsitecnt += n+m;
	      Tmisscnt += n+m-2;

	      FLOAT err = x-y;
	      FLOAT var = QUADRATIC_VARIANCE ? varSF2 + varSD*(x+y) + varSR*(x*x+y*y) : varSF2 + varSD*(x+y);
	      Terrnorm += err*err/var;
	      Tsegcnt++;

	      if(RefRepeats > 0){  /* append error information of current interval to errors[0..numerrors-1] */
                #pragma omp critical
		{
		  if(numerrors >= maxerrors){
		    int newmax = maxerrors*2;
		    Cinterval *newerrors = new Cinterval[newmax];
		    for(size_t i= 0;i < numerrors;i++)
		      newerrors[i] = errors[i];
		    delete [] errors;
		    errors = newerrors;
		    maxerrors = newmax;
		  }

		  Cinterval *perr = &errors[numerrors++];
		  perr->x = x;
		  perr->y = y;
		  if(DEBUG) assert(x > 0.0 && y > 0.0);

		  /* append misaligned sites information of current interval to misaligns[0..nummisaligns-1] */
		  if(nummisaligns >= maxmisaligns-1){
		    int newmax = maxmisaligns*2;
		    Cmisalign *newmisaligns = new Cmisalign[newmax];
		    for(size_t i= 0; i < nummisaligns; i++)
		      newmisaligns[i] = misaligns[i];
		    delete [] misaligns;
		    misaligns = newmisaligns;
		    maxmisaligns = newmax;
		  }

		  Cmisalign *pmis = &misaligns[nummisaligns];
		  pmis[0].x = y;
		  pmis[0].m = (I-G)-1;
		  pmis[1].x = x;
		  pmis[1].m = (J-H)-1;
		  nummisaligns += 2;
		}
	      }
	    }  
	  }
	} /* omp for nowait PairId = 0 .. NumPair -1 */

	if(RA_MIN_TIME) {/* fully release memory  : previous ReleaseMemory() calls only released memory every RA_MIN_TIME msec (and only if VmRSS exceeds RA_MAX_MEM * maxmem * 0.9) */
	  /* Make sure last thread is not stuck in WaitForMemory() loop, even if it requires more than -maxmem memory : 
	     Incrementing PausedThreads for completed threads allows WaitForMemory() to recognize when all threads are either paused (in WaitForMemory()) or finished */
          #pragma omp atomic update
	  PausedThreads++;

	  ReleaseMemory(Mem, Mem->MemUsed, tid, -1, -1, 0,0,Fmem, 1/* force*/);
        }

	#pragma omp critical
	{
	  mapcnt += Tmapcnt;
	  paircnt += Tpaircnt;
	  tcnt += Ttcnt;
	  poscnt += Tposcnt;
	  sitecnt += Tsitecnt;
	  misscnt += Tmisscnt;
	  outliertot += Toutliertot;
	  outliercnt += Toutliercnt;
	  len += Tlen;
	  lenB += TlenB;
	  outliersum += Toutliersum;
	  segcnt += Tsegcnt;
	  errnorm += Terrnorm;
	  if(OMP_DEBUG) LoopCnt += myLoopCnt;
	}
      } // end of parallel section

      PausedThreads = 0;
      
      if(DEBUG && RA_MIN_TIME && RA_MIN_MEM < 8){
	if(VMEMCHECK && DEBUG>=1+RELEASE && vmemTotsiz){
	  printf("vmemTotsiz= %lld\n",vmemTotsiz);
	  fflush(stdout);
	  assert(vmemTotsiz == 0);
        }
	if(MEMRESERVE && DEBUG>=1+RELEASE && VmRSSreserved > 0){
	  printf("VmRSSreserved= %lld\n",VmRSSreserved);
	  fflush(stdout);
	  assert(VmRSSreserved==0);
	}
	
	if(VERB/* HERE HERE >=2 */ && PairMerge){
	  long long VmSize,VmRSS,VmSwap,VmHWM;// local copies, so global shared values are not modified
	  getmem(VmSize,VmRSS,VmSwap,&VmHWM);

	  printf("After parallel alignment loop: VmRSS= %0.4f(HWM= %0.4f), VmSwap= %0.4f, VmSize= %0.4f Gb: CPU time= %0.6f, wall time= %0.6f secs\n",
		 VmRSS*1e-9,VmHWM*1e-9,VmSwap*1e-9,VmSize*1e-9,mtime(),wtime());
	  fflush(stdout);
	}
      }

      if(OMP_DEBUG) assert(LoopCnt == NumPair1);

      if(VERB & (VERB>=2 || HASH_DEBUG) && numthreads > 1){
	printf("%d:Completed parallel section, numaligns=%lld,mapcnt=%lld,NumPair1=%d,Yid=%d: cum CPU time=%0.6f, wall time=%0.6f\n", 
	       PairMergeRepeat ? pairmergeIter : wcnt,(long long)numaligns,mapcnt,NumPair1,Yid,mtime(), wtime());
	fflush(stdout);
      }

      if(DEBUG>=2 && PairSplit){
	for(size_t A = 0; A < numaligns; A++){
	  Calign *align = alignment[A];
	  assert(align != NULL);
	  if(!(0 <= align->mapid1 && align->mapid1 < numY)){
	    printf("alignment[%lld]=0x%p: mapid1=%d,mapid2=%d,score=%0.6f,logPV=%0.2f,numpairs=%d\n",
		   (long long)A,align,align->mapid1,align->mapid2,align->score,align->logPV,align->numpairs);
	    fflush(stdout);
	    assert(0 <= align->mapid1 && align->mapid1 < numY);
	  }
	  assert(0 <= align->mapid2 && align->mapid2 < numX);
	}
	for(size_t A = numaligns; A < numaligns + NumPair1; A++){
	  Calign *align = alignment[A];
	  if(align){
	    assert(0 <= align->mapid1 && align->mapid1 < numY);
	    assert(0 <= align->mapid2 && align->mapid2 < numX);
	  }
	}
      }

      numaligns += NumPair1;
      NumPair1 = 0;
    }/* Yid = 0 .. numY-1 */

    if(MAXALIGN_THRESHOLD){
      printf("max alignments per map = %d\n", maxaligns_permap);
      fflush(stdout);
    }

    delete [] YidList1;
    delete [] XidList1;
    delete [] phashList1;
    XidList1 = YidList1 = NULL;
    phashList1 = NULL;
    if(HASH_STREAM>=2){
      delete [] YidList2;
      delete [] XidList2;
      delete [] phashList2;
      YidList2 = XidList2 = NULL;
      phashList2 = NULL;
    }
    if(HashFP){
      hash_close(HashFP,1);
      HashFP = NULL;
      MatchMax = 0;
      totalhashpairs = 0;
    }

    if(VERB>=2 && !PairSplit){
      if(FirstAlignments > 0)
	printf("alignments found=%lld (terminated after %d alignments or more)\n", mapcnt,FirstAlignments);
      else
	printf("%d:Total good alignments found=%lld (numaligns=%lu): cum wall time=%0.6f\n", PairMergeRepeat ? pairmergeIter : wcnt, mapcnt, numaligns, wtime());
     fflush(stdout);
    }

    if(!PairSplit)
      break;
    
    if(FirstAlignments > 0)
      break;

    /* scan all alignments for highest scoring local alignment : exclude any that are below threshold OR (pairs already confirmed with chimpair = -1, or below threshold/discarded with chimpair = 1) */
    /* also compress the alignment[] array to remove all NULL entries to save memory */
    size_t B = (LVERB ? 0 : numaligns_compressed), aligncnt = (LVERB ? 0 : prev_aligncnt);// numaligns_compressed in number of alignments previously checked for NULL entries or below threshold etc */
    size_t orignumaligns_hashed = numaligns_hashed;
    for(size_t A = B; A < numaligns; A++){
      Calign *align = alignment[A];
      if(!align)
	continue;

      if(DEBUG>=2) assert(0 <= align->mapid1 && align->mapid1 < numY);
      if(DEBUG>=2) assert(0 <= align->mapid2 && align->mapid2 < numX);
      if(B < A){/* swap alignment[A] and alignment[B] */
	if(DEBUG>=2) assert(alignment[B]==NULL);
	alignment[B] = align;
	alignment[A] = NULL;
      }
      if(orignumaligns_hashed == A+1)
	numaligns_hashed = B+1;
      B++;

      if(align->chimpair)
	continue;
      if(!(align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= AlignedSiteThreshold)){
	align->chimpair = 1;
	continue;
      }

      int Yid = align->mapid1;
      int Xid = align->mapid2;
      int orientation = align->orientation;
      if(VERB>=2 || EVERB){
	aligncnt++;
	if(EVERB){
	  printf("Alignment %lu (alignment[%lu]):Yid=%d,Xid=%d,or=%d:numpairs=%d,score=%0.6f,logPV=%0.2f\n",
		 aligncnt,A,Yid,Xid,orientation,align->numpairs,align->score,align->logPV);
	  fflush(stdout);
	}
      }
    }

    if(VERB>=2){
      printf("%d:Total alignments found=%lld:numX=%d,numY=%d,numaligns=%lld->%lld(good=%lld)\n", wcnt, mapcnt,numX,numY,(long long int)numaligns,(long long int)B,(long long int)aligncnt);
      fflush(stdout);
    }
    numaligns_compressed = numaligns = B;
    prev_aligncnt = aligncnt;

    if(DEBUG>=2){
      for(size_t A = 0; A < numaligns; A++){
	Calign *align = alignment[A];
	assert(align != NULL);
	if(!(0 <= align->mapid1 && align->mapid1 < numY)){
	  printf("alignment[%lld]=0x%p: mapid1=%d,mapid2=%d,score=%0.6f,logPV=%0.2f,numpairs=%d\n",
		 (long long)A,align,align->mapid1,align->mapid2,align->score,align->logPV,align->numpairs);
	  fflush(stdout);
	  assert(0 <= align->mapid1 && align->mapid1 < numY);
	}
	assert(0 <= align->mapid2 && align->mapid2 < numX);
      }
    }

    size_t splitcnt = 0;/* number of best alignments found and used to split maps */

    int numthreads = 1;
    #ifdef _OPENMP
    numthreads = maxthreads;
    if(numthreads > NumPair1)
      numthreads = max(1,NumPair1);
    #endif

    if(VERB>=2){
      printf("numthreads=%d,PairsplitSeparateQueries=%d,orignumX=%d\n",numthreads,PairsplitSeparateQueries,orignumX);
      fflush(stdout);
    }

    int rootXidInc = (PairsplitSeparateQueries ? 1 : orignumX);

    int LoopCnt = 0;

    #pragma omp parallel num_threads(numthreads) if(numthreads > 1 && PairsplitSeparateQueries && !LVERB)
    {
      int myLoopCnt = 0;

      #pragma omp for schedule(dynamic,1)
      for(int rootXid = 0; rootXid < orignumX; rootXid += rootXidInc){
	if(OMP_DEBUG) myLoopCnt++;
	/* find best alignment with query rootXid == YYmap[Yid]->Xid */

      Lnextbest:
	Calign *bestalign = 0;
	FLOAT bestscore = LARGE_NEGATIVE;
	for(size_t A = 0; A < numaligns; A++){
	  Calign *align = alignment[A];
	  if(DEBUG>=2) assert(align);
	  if(DEBUG>=2) assert(0 <= align->mapid1 && align->mapid1 < numY);
	  if(DEBUG>=2) assert(0 <= align->mapid2 && align->mapid2 < numX);

	  if(align->chimpair)
	    continue;
	  if(PairsplitSeparateQueries && XXmap[align->mapid2]->Xid != rootXid)
	    continue;

	  if(DEBUG>=2) assert(align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= AlignedSiteThreshold);

	  if(PairsplitPV){
	    if(align->logPV > bestscore){
	      bestscore = align->logPV;
	      bestalign = align;
	    }
	  } else {
	    if(align->score > bestscore){
	      bestscore = align->score;
	      bestalign = align;
	    }
	  }
	}

	if(!bestalign)/* no more alignments above threshold for query map XXmap[rootXid] */
	  continue;/* next rootXid */

        #pragma omp atomic
	splitcnt++;

	Calign *align = bestalign;
	int Yid = align->mapid1;
	if(DEBUG && PairsplitSeparateQueries) assert(YYmap[Yid]->Xid == rootXid);
	int Xid = align->mapid2;
	Cmap *Ymap = YYmap[Yid];
	Cmap *Xmap = XXmap[Xid];
	int K = align->numpairs;
	int N = Ymap->numsite[0];
	int M = Xmap->numsite[0];
	FLOAT *Y = Ymap->site[0];
	FLOAT *X = Xmap->site[0];

	LocalCnt++;
	align->chimpair = -1;   /*  marks this alignment as confirmed (for output) and will be ignored in subsequent iterations of this loop
				    NOTE : If either Xid or Yid are split, this alignment will be changed to chimpair=1 (discarded alignments, will not be output) */
	if(LVERB){
	  printf("Query=%d:Next best local alignment %d:Yid=%d(id=%lld,left=%d,right=%d,origN=%d,N=%d),Xid=%d(id=%lld,left=%d,right=%d,origM=%d,M=%d),OR=%d:score=%0.4f,logPV=%0.4f,sites=%d:(%d,%d) to (%d,%d),Lend=%d,Rend=%d\n",
		 rootXid,LocalCnt,Yid,Ymap->id,Ymap->origmap ? Ymap->left[0] : 0, Ymap->origmap ? Ymap->right[0] : N+1, Ymap->origmap ? Ymap->origmap->numsite[0] : N, N,
		 Xid,Xmap->id,Xmap->origmap ? Xmap->left[0] : 0, Xmap->origmap ? Xmap->right[0] : M+1, Xmap->origmap ? Xmap->origmap->numsite[0] : M, M,
		 align->orientation, align->score,align->logPV,K,align->sites1[0],align->sites2[0],align->sites1[K-1],align->sites2[K-1],align->Lend,align->Rend);
	  fflush(stdout);
	}

	/* track center regions, in case they need to be cloned (to detect repeats) */
	int Yleft = 0, Yright = N+1;
	int Xleft = 0, Xright = M+1;

	if(PAIRSPLIT_EXTEND && DEBUG){
	  assert(0 < align->refstart);
	  assert(0 < align->start);
	  assert(align->refstart <= align->refend);
	  assert(align->start <= align->end);
	  assert(align->refend <= N);
	  assert(align->end <= M);
	}

	if(align->Rend <= -1 || align->Lend <= -1){/* check if Ymap,Xmap need to be broken up */
	  int I = align->sites1[0];
	  int J = align->sites2[0];
	  int Imax = align->sites1[K-1];
	  int Jmax = align->sites2[K-1];

	  if(PAIRSPLIT_EXTEND){
	    if(align->Lend <= -1 && (align->refstart < I || align->start < J)){/* flag the left end as if it was deleted (even if AlignedSiteThreshold or map->paired prevents it) to make
										  sure the center is cloned and the alignment rechecked/extended */
	      if(!Ymap->paired){
		if(LVERB){
		  printf("Flagged left end of Yid=%d(id=%lld) with Yleft=%d\n",Yid,YYmap[Yid]->id,I);
		  fflush(stdout);
		}
		Yleft = I;
	      }
	      if(!Xmap->paired){
		if(LVERB){
		  printf("Flagged left end of Xid=%d(id=%lld) with Xleft=%d\n",Xid,XXmap[Xid]->id,J);
		  fflush(stdout);
		}
		Xleft = J;
	      }
	    }
	    if(align->Rend <= -1 && (align->refend > Imax || align->end > Jmax)){/* flag the right end as if it was deleted (even if AlignedSiteThreshold or map->paired prevents it) to make
										    sure the center is cloned and the alignment rechecked/extended */
	      if(!Ymap->paired){
		if(LVERB){
		  printf("Flagged right end of Yid=%d(id=%lld) with Yright=%d\n",Yid,YYmap[Yid]->id,Imax);
		  fflush(stdout);
		}
		Yright = Imax;
	      }
	      if(!Xmap->paired){
		if(LVERB){
		  printf("Flagged right end if Xid=%d(id=%lld) with Xright=%d\n",Xid,XXmap[Xid]->id,Jmax);
		  fflush(stdout);
		}
		Xright = Jmax;
	      }
	    }
	  }

	  if(align->Lend <= -1 && I-1 >= AlignedSiteThreshold && !Ymap->paired){/* save copy of left end of Ymap up to site I-1 and append it to end of YYmaps[]/refmap[] */
	    Yleft = I;
	    int Iend = (PAIRSPLIT_OVERLAP && PAIRSPLIT_EXTEND && align->Lend <= -1 && (align->refstart < I || align->start < J)) ? min(N,I+1) : I;/* extend map to site Iend-1 */

	    Cmap *newmap = NULL;

#pragma omp critical(numrefmaps)
	    {
	      maxmapalloc(numrefmaps+1,maxrefmaps,refmap,0,1);
	      YYmap = refmap;/* in case refmap[] was reallocated */
	      newmap = refmap[numrefmaps++];
	      if(!newmap)
		refmap[numrefmaps-1] = newmap = new Cmap;
	      newmap->mapid = numY;
	      //	    YYmap[numY++] = newmap;
	      if(DEBUG) assert(YYmap[numY] == newmap);
	      numY++;
	      if(DEBUG) assert(numY == numrefmaps);
	    }
	    newmap->allfree();
	    if(DEBUG) assert(newmap->blockmem==0);
	    newmap->centercloned = 0;
	    newmap->paired = 0;
	    newmap->origmap = Ymap;
	    if(PairsplitSeparateQueries)
	      newmap->Xid = Ymap->Xid;
	    newmap->id = Ymap->id;
	    newmap->fileid = Ymap->fileid;
	    newmap->numsite[0] = Iend - 1;
	    newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
	    if(LEAKDEBUG>=2){
	      printf("refmap[%d]=%p:id=%lld,site[0]->%p,numsite[0]=%d\n",
		     numrefmaps-1,newmap,newmap->id,newmap->site[0],newmap->numsite[0]);
	      fflush(stdout);
	    }
	    newmap->site[0][0] = 0.0;
	    newmap->left[0] = 0;
	    newmap->right[0] = Iend - 1;
	    for(int i = 1; i < Iend; i++)
	      newmap->site[0][i] = Y[i];
	    newmap->site[0][Iend] = newmap->site[0][Iend - 1] + MININTERVAL;

	    newmap->startloc = Ymap->startloc;
	    newmap->endloc = Ymap->endloc;
	    newmap->len = Ymap->len * (newmap->site[0][Iend]/Y[N+1]);
	    newmap->fpcnt = -1;

	    newmap->deltaXinit();

	    if(Ymap->origmap){
	      if(Ymap->left[0] > 0 && newmap->left[0] == 0)
		newmap->left[0] = 1;
	      if(Ymap->right[0] <= Ymap->origmap->numsite[0] && newmap->right[0] > N)
		newmap->right[0] = N;
	      newmap->origmap = Ymap->origmap;
	      if(DEBUG) assert(newmap->left[0] >= 0);
	      newmap->left[0] += max(0,Ymap->left[0] - 1);
	      newmap->right[0] += max(0,Ymap->left[0] - 1);
	      if(DEBUG && !(newmap->origmap->origmap == 0)){
		printf("Split off left end of Yid=%d left of I=%d(%0.3f) to YYmap[%d]:id=%lld(file=%d),left=%d,right=%d,origN=%d,N=%d\n",
		       Yid,Iend,Y[Iend],numY-1,Ymap->id, Ymap->fileid,newmap->left[0],newmap->right[0],newmap->origmap->numsite[0],newmap->numsite[0]);
		if(DEBUG) assert(newmap->origmap!=0);
		printf("YYmap[Yid]->origmap=%p,newmap->origmap=%p,newmap->origmap->origmap=%p\n",YYmap[Yid]->origmap,newmap->origmap,newmap->origmap->origmap);
		fflush(stdout);
		assert(newmap->origmap->origmap == 0);
	      }
	    }

	    if(LVERB){
	      printf("Split off left end of Yid=%d left of I=%d(%0.3f) to YYmap[%d]:id=%lld(file=%d),left=%d,right=%d,origN=%d,N=%d\n",
		     Yid,Iend,Y[Iend],numY-1,Ymap->id, Ymap->fileid,newmap->left[0],newmap->right[0],newmap->origmap->numsite[0],newmap->numsite[0]);
	      fflush(stdout);
	    }
	    if(DEBUG>=2 && newmap->origmap && newmap->numsite[0] >= 2){/* check that mapping is consistent by checking the interval sizes */
	      int shift = max(1,newmap->left[0])-1;
	      assert(newmap->site[0][0] == 0.0);
	      assert(newmap->site[0][1] > 0.0);
	      assert(newmap->site[0][newmap->numsite[0]+1] > newmap->site[0][newmap->numsite[0]]);
	      for(int T = 1; T < newmap->numsite[0];T++){
		double newsize = newmap->site[0][T+1] - newmap->site[0][T];
		double oldsize = newmap->origmap->site[0][shift + T+1] - newmap->origmap->site[0][shift + T]; 
		if(fabs(newsize - oldsize) >= 0.001 || newsize < 0.0 || oldsize < 0.0){
		  for(int K = 1; K < newmap->numsite[0]; K++)
		    printf("YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f); YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f)\n",
			   numY-1,K,newmap->site[0][K],K+1,newmap->site[0][K+1],newmap->site[0][K+1]-newmap->site[0][K],
			   Yid,shift+K,newmap->origmap->site[0][shift+K],shift+K+1,newmap->origmap->site[0][shift+K+1],
			   newmap->origmap->site[0][shift+K+1]-newmap->origmap->site[0][shift+K]);
		  printf("YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f); YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f),shift=%d,newmap->left=%d\n",
			 numY-1,T,newmap->site[0][T],T+1,newmap->site[0][T+1],newsize,
			 Yid,shift+T,newmap->origmap->site[0][shift+T],shift+T+1,newmap->origmap->site[0][shift+T+1],oldsize,shift,newmap->left[0]);
		  fflush(stdout);
		  assert(!(fabs(newsize-oldsize) >= 0.001));
		}
	      }
	    }
	  }

	  if(align->Lend <= -1 && J-1 >= AlignedSiteThreshold && !Xmap->paired){/* save copy of left end of Xmap up to site J-1 and append it to end of XXmaps[]/map[] */
	    Xleft = J;
	    int Jend = (PAIRSPLIT_OVERLAP && PAIRSPLIT_EXTEND && align->Lend <= -1 && (align->refstart < I || align->start < J)) ? min(M,J+1) : J;/* extend map to site Jend-1 */

	    Cmap *newmap = NULL;

#pragma omp critical(nummaps)
	    {
	      maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	      XXmap = &Gmap[Cfirst];/* in case Gmap[] was reallocated */
	      newmap = Gmap[nummaps++];
	      if(!newmap)
		Gmap[nummaps-1] = newmap = new Cmap;
	      newmap->mapid = numX;
	      //	    XXmap[numX++] = newmap;
	      if(DEBUG) assert(XXmap[numX] == newmap);
	      numX++;
	      if(DEBUG) assert(numX == nummaps - Cfirst);
	    }

	    newmap->allfree();
	    if(DEBUG) assert(newmap->blockmem==0);
	    newmap->centercloned = 0;
	    newmap->paired = 0;
	    newmap->origmap = Xmap;
	    if(PairsplitSeparateQueries)
	      newmap->Xid = Xmap->Xid;
	    newmap->id = Xmap->id;
	    newmap->fileid = Xmap->fileid;
	    newmap->numsite[0] = Jend - 1;
	    newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
	    if(LEAKDEBUG>=2){
	      printf("Gmap[%d]=%p:id=%lld,site[0]->%p,numsite[0]=%d\n",
		     nummaps-1,newmap,newmap->id,newmap->site[0],newmap->numsite[0]);
	      fflush(stdout);
	    }
	    newmap->site[0][0] = 0.0;
	    if(!align->orientation){
	      newmap->left[0] = 0;
	      newmap->right[0] = Jend - 1;
	      for(int i= 1; i < Jend; i++)
		newmap->site[0][i] = X[i];
	      newmap->site[0][Jend] = newmap->site[0][Jend - 1] + MININTERVAL;
	    } else {
	      newmap->right[0] = M+1;
	      newmap->left[0] = M+1-(Jend - 1);
	      if(DEBUG) assert(newmap->left[0] >= 0);
	      for(int i = 1; i <= Jend; i++)
		newmap->site[0][i] = MININTERVAL + X[newmap->left[0] + i-1] - X[newmap->left[0]];
	    }
	
	    newmap->startloc = Xmap->startloc;
	    newmap->endloc = Xmap->endloc;
	    newmap->len = Xmap->len * (newmap->site[0][Jend]/X[M+1]);
	    newmap->fpcnt = -1;

	    newmap->deltaXinit();

	    if(Xmap->origmap){
	      if(Xmap->left[0] > 0 && newmap->left[0] == 0)
		newmap->left[0] = 1;
	      if(Xmap->right[0] <= Xmap->origmap->numsite[0] && newmap->right[0] > M)
		newmap->right[0] = M;
	      newmap->origmap = Xmap->origmap;
	      if(DEBUG) assert(newmap->left[0] >= 0);
	      newmap->left[0] += max(0,Xmap->left[0] - 1);
	      newmap->right[0] += max(0,Xmap->left[0] - 1);
	      if(DEBUG) assert(newmap->origmap->origmap == 0);
	    }
	    if(LVERB){
	      printf("Split off left end of Xid=%d(flip=%d) left of J=%d(%0.3f) to XXmap[%d]:id=%lld(file=%d),left=%d,right=%d,origM=%d,M=%d\n",
		     Xid,align->orientation,Jend,align->orientation ? X[M+1]-X[M+1-Jend] : X[Jend], numX-1,Xmap->id,Xmap->fileid,newmap->left[0],newmap->right[0],newmap->origmap->numsite[0],newmap->numsite[0]);
	      fflush(stdout);
	    }
	    if(DEBUG>=2 && newmap->origmap && newmap->numsite[0] >= 2){/* check that mapping is consistent by checking the interval sizes */
	      assert(newmap->site[0][0] == 0.0);
	      assert(newmap->site[0][1] > 0.0);
	      assert(newmap->site[0][newmap->numsite[0]+1] > newmap->site[0][newmap->numsite[0]]);
	      for(int T = 1; T < newmap->numsite[0];T++){
		int shift = max(1,newmap->left[0])-1;
		double newsize = newmap->site[0][T+1] - newmap->site[0][T];
		double oldsize = newmap->origmap->site[0][shift + T+1] - newmap->origmap->site[0][shift + T]; 
		if(fabs(newsize - oldsize) >= 0.001 || newsize < 0.0 || oldsize < 0.0){
		  printf("XXmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f); YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f)\n",
			 numX-1,T,newmap->site[0][T],T+1,newmap->site[0][T+1],newsize,
			 Xid,shift+T,newmap->origmap->site[0][shift+T],shift+T+1,newmap->origmap->site[0][shift+T+1],oldsize);
		  fflush(stdout);
		  assert(!(fabs(newsize-oldsize) >= 0.001));
		}
	      }
	    }
	  }

	  if(align->Rend <= -1 && N-Imax >= AlignedSiteThreshold && !Ymap->paired && PairSplit > 0){/* save copy of right end of Ymap starting at site Imax+1 and append it to end of YYmaps[]/refmap[] */
	    Yright = Imax;
	    int Istart = (PAIRSPLIT_OVERLAP && PAIRSPLIT_EXTEND && align->Rend <= -1 && (align->refend > Imax || align->end > Jmax)) ? max(1,Imax-1) : Imax;/* start at site Istart+1 */
	    Cmap *newmap = NULL;

#pragma omp critical(numrefmaps)
	    {
	      maxmapalloc(numrefmaps+1,maxrefmaps,refmap,0,1);
	      YYmap = refmap;/* in case refmap[] was reallocated */
	      newmap = refmap[numrefmaps++];
	      if(!newmap)
		refmap[numrefmaps-1] = newmap = new Cmap;
	      newmap->mapid = numY;
	      //	    YYmap[numY++] = newmap;
	      if(DEBUG) assert(YYmap[numY] == newmap);
	      numY++;
	      if(DEBUG) assert(numY == numrefmaps);
	    }
	    newmap->allfree();
	    if(DEBUG) assert(newmap->blockmem==0);
	    newmap->centercloned = 0;
	    newmap->paired = 0;
	    newmap->origmap = Ymap;
	    if(PairsplitSeparateQueries)
	      newmap->Xid = Ymap->Xid;
	    newmap->id = Ymap->id;
	    newmap->fileid = Ymap->fileid;
	    newmap->numsite[0] = N - Istart;
	    newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
	    if(LEAKDEBUG>=2){
	      printf("refmap[%d]=%p:id=%lld,site[0]->%p,numsite[0]=%d\n",
		     numrefmaps-1,newmap,newmap->id,newmap->site[0],newmap->numsite[0]);
	      fflush(stdout);
	    }
	    newmap->site[0][0] = 0.0;
	    newmap->left[0] = Istart + 1;
	    if(DEBUG) assert(newmap->left[0] >= 0);
	    newmap->right[0] = N+1;
	    for(int i = 1; i <= N+1 - Istart; i++)
	      newmap->site[0][i] = MININTERVAL + Y[newmap->left[0] + i-1] - Y[newmap->left[0]];
	
	    newmap->startloc = Ymap->startloc;
	    newmap->endloc = Ymap->endloc;
	    newmap->len = Ymap->len * (newmap->site[0][N-Istart+1]/Y[N+1]);
	    newmap->fpcnt = -1;

	    newmap->deltaXinit();

	    if(Ymap->origmap){
	      if(Ymap->left[0] > 0 && newmap->left[0] == 0)
		newmap->left[0] = 1;
	      if(Ymap->right[0] <= Ymap->origmap->numsite[0] && newmap->right[0] > N)
		newmap->right[0] = N;
	      newmap->origmap = Ymap->origmap;
	      if(DEBUG) assert(newmap->left[0] >= 0);
	      newmap->left[0] += max(0,Ymap->left[0] - 1);
	      newmap->right[0] += max(0,Ymap->left[0] - 1);
	      if(DEBUG) assert(newmap->origmap->origmap == 0);
	    }
	    if(LVERB){
	      printf("Split off right end of Yid=%d right of I=%d(%0.3f) as YYmap[%d]:id=%lld(file=%d),left=%d,right=%d,origN=%d,N=%d\n",
		     Yid,Istart,Y[Istart],numY-1,Ymap->id,Ymap->fileid,newmap->left[0],newmap->right[0],newmap->origmap->numsite[0],newmap->numsite[0]);
	      fflush(stdout);
	    }
	    if(DEBUG>=2 && newmap->origmap && newmap->numsite[0] >= 2){/* check that mapping is consistent by checking the interval sizes */
	      assert(newmap->site[0][0] == 0.0);
	      assert(newmap->site[0][1] > 0.0);
	      assert(newmap->site[0][newmap->numsite[0]+1] > newmap->site[0][newmap->numsite[0]]);
	      for(int T = 1; T < newmap->numsite[0];T++){
		int shift = max(1,newmap->left[0])-1;
		double newsize = newmap->site[0][T+1] - newmap->site[0][T];
		double oldsize = newmap->origmap->site[0][shift + T+1] - newmap->origmap->site[0][shift + T]; 
		if(fabs(newsize - oldsize) >= 0.001 || newsize < 0.0 || oldsize < 0.0){
		  printf("YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f); YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f)\n",
			 numY-1,T,newmap->site[0][T],T+1,newmap->site[0][T+1],newsize,
			 Yid,shift+T,newmap->origmap->site[0][shift+T],shift+T+1,newmap->origmap->site[0][shift+T+1],oldsize);
		  fflush(stdout);
		  assert(!(fabs(newsize-oldsize) >= 0.001));
		}
	      }
	    }
	  }

	  if(align->Rend <= -1 && M-Jmax >= AlignedSiteThreshold && !Xmap->paired){/* save copy of right end of Xmap starting at site Jmax+1 and append it to end of XXmaps[]/map[] */
	    Xright = Jmax;
	    int Jstart = (PAIRSPLIT_OVERLAP && PAIRSPLIT_EXTEND && align->Rend <= -1 && (align->refend > Imax || align->end > Jmax)) ? max(1,Jmax-1) : Jmax;/* start at site Jstart+1 */
	    Cmap *newmap = NULL;
        

#pragma omp critical(nummaps)
	    {
	      maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	      XXmap = &Gmap[Cfirst];/* in case Gmap[] was reallocated */
	      newmap = Gmap[nummaps++];
	      if(!newmap)
		Gmap[nummaps-1] = newmap = new Cmap;
	      newmap->mapid = numX;
	      //	    XXmap[numX++] = newmap;
	      if(DEBUG) assert(XXmap[numX] == newmap);
	      numX++;
	      if(DEBUG) assert(numX == nummaps - Cfirst);
	    }
	    newmap->allfree();
	    if(DEBUG) assert(newmap->blockmem==0);
	    newmap->centercloned = 0;
	    newmap->paired = 0;
	    newmap->origmap = Xmap;
	    if(PairsplitSeparateQueries)
	      newmap->Xid = Xmap->Xid;
	    newmap->id = Xmap->id;
	    newmap->fileid = Xmap->fileid;
	    newmap->numsite[0] = M - Jstart;
	    newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
	    if(LEAKDEBUG>=2){
	      printf("Gmap[%d]=%p:id=%lld,site[0]->%p,numsite[0]=%d\n",
		     nummaps-1,newmap,newmap->id,newmap->site[0],newmap->numsite[0]);
	      fflush(stdout);
	    }
	    newmap->site[0][0] = 0.0;
	    if(!align->orientation){
	      newmap->left[0] = Jstart + 1;
	      if(DEBUG) assert(newmap->left[0] >= 0);
	      newmap->right[0] = M+1;
	      for(int i = 1; i <= M+1 - Jstart; i++)
		newmap->site[0][i] = MININTERVAL + X[newmap->left[0] + i-1] - X[newmap->left[0]];
	    } else {
	      newmap->left[0] = 0;
	      newmap->right[0] = M - Jstart;
	      for(int i = 1; i <= M - Jstart; i++)
		newmap->site[0][i] = X[i];
	      newmap->site[0][M - Jstart + 1] = newmap->site[0][M - Jstart] + MININTERVAL;
	    }

	    newmap->startloc = Xmap->startloc;
	    newmap->endloc = Xmap->endloc;
	    newmap->len = Xmap->len * (newmap->site[0][M - Jstart + 1]/X[M+1]);
	    newmap->fpcnt = -1;

	    newmap->deltaXinit();

	    if(Xmap->origmap){
	      if(Xmap->left[0] > 0 && newmap->left[0] == 0)
		newmap->left[0] = 1;
	      if(Xmap->right[0] <= Xmap->origmap->numsite[0] && newmap->right[0] > M)
		newmap->right[0] = M;

	      newmap->origmap = Xmap->origmap;
	      if(DEBUG) assert(newmap->left[0] >= 0);
	      newmap->left[0] += max(0, Xmap->left[0] - 1);
	      newmap->right[0] += max(0, Xmap->left[0] - 1);
	      if(DEBUG) assert(newmap->origmap->origmap == 0);
	    }

	    if(LVERB){
	      printf("Split off right end of Xid=%d(flip=%d) right of J=%d(%0.3f) to XXmap[%d]:id=%lld(file=%d),left=%d,right=%d,origM=%d,M=%d\n",
		     Xid,align->orientation,Jstart,align->orientation ? X[M+1]-X[M+1-Jstart] : X[Jstart],numX-1,Xmap->id,Xmap->fileid,newmap->left[0],newmap->right[0],newmap->origmap->numsite[0],newmap->numsite[0]);
	      fflush(stdout);
	    }

	    if(DEBUG>=2 && newmap->origmap && newmap->numsite[0] >= 2){/* check that mapping is consistent by checking the interval sizes */
	      assert(newmap->site[0][0] == 0.0);
	      assert(newmap->site[0][1] > 0.0);
	      assert(newmap->site[0][newmap->numsite[0]+1] > newmap->site[0][newmap->numsite[0]]);
	      for(int T = 1; T < newmap->numsite[0];T++){
		int shift = max(1,newmap->left[0])-1;
		double newsize = newmap->site[0][T+1] - newmap->site[0][T];
		double oldsize = newmap->origmap->site[0][shift + T+1] - newmap->origmap->site[0][shift + T]; 
		if(fabs(newsize - oldsize) >= 0.001 || newsize < 0.0 || oldsize < 0.0){
		  printf("XXmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f); YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f)\n",
			 numX-1,T,newmap->site[0][T],T+1,newmap->site[0][T+1],newsize,
			 Xid,shift+T,newmap->origmap->site[0][shift+T],shift+T+1,newmap->origmap->site[0][shift+T+1],oldsize);
		  fflush(stdout);
		  assert(!(fabs(newsize-oldsize) >= 0.001));
		}
	      }
	    }
	  }
	}

	if(!((Yleft > 0 || Yright <= N) || (Xleft > 0 || Xright <= M))){/* no splitting hence no center cloning : lock Xmap and Ymap to prevent further splitting */
	  Xmap->paired = Yid + 1;
	  Ymap->paired = Xid + 1;
	  if(LVERB){
	    printf("Locked Xid=%d(id=%lld) and Yid=%d(id=%lld) from further splitting\n",Xid,XXmap[Xid]->id,Yid,YYmap[Yid]->id);
	    fflush(stdout);
	  }
	  goto Lnextbest;/* No need to redo alignments: just locate next best alignment above threshold(s) */

	} else {
	  /* filter all alignments with either Yid or Xid if they were split (their center will be cloned) : use chimpair=1 to flag invalid alignments */
	  /* NOTE : Perhaps only alignments involving the center region to be cloned (or a part thereof) should be filtered out ??? */
	  for(size_t A = 0; A < numaligns; A++){
	    Calign *palign = alignment[A];
	    if(!palign || palign->chimpair)
	      continue;
	    if(palign->mapid1 == Yid && (Yleft > 0 || Yright <= N)){
	      if(LVERB){
		int IL = palign->sites1[0];
		int IR = palign->sites1[palign->numpairs-1];
		printf("supressing alignment[%lld](score=%0.6f,logPV=%0.2f,numpairs=%d) between Yid=%d(id=%lld) and Xid=%d(id=%lld) (Yrange = %d to %d) due to cloning of Yid from Yleft=%d to Yright=%d, N=%d\n",
		       (long long)A,palign->score,palign->logPV,palign->numpairs,Yid,YYmap[Yid]->id,palign->mapid2,XXmap[palign->mapid2]->id, IL, IR, Yleft,Yright, N);
		fflush(stdout);
	      }
	      palign->chimpair = 1;
	      continue;
	    }
	    if(palign->mapid2 == Xid && (Xleft > 0 || Xright <= M)){
	      if(LVERB){
		int JL = palign->sites2[0];
		int JR = palign->sites2[palign->numpairs-1];
		printf("supressing alignment[%lld](score=%0.6f,logPV=%0.2f,numpairs=%d) between Yid=%d(id=%lld) and Xid=%d(id=%lld),or=%d, (Xrange = %d to %d) due to cloning of Xid (orientation=%d) from Xleft=%d to Xright=%d, M=%d\n",
		       (long long)A,palign->score,palign->logPV,palign->numpairs,palign->mapid1,YYmap[palign->mapid1]->id,Xid,XXmap[Xid]->id, palign->orientation, JL, JR, align->orientation, Xleft, Xright, M);
		fflush(stdout);
	      }
	      palign->chimpair = 1;
	      continue;
	    }
	  }

	  /* Clone center regions */
	  align->chimpair = 1;/* this alignment will be found again with the Cloned center region(s) */
	  if(LVERB){
	    printf("supressing alignment (score=%0.6f,logPV=%0.2f,numpairs=%d) between Yid=%d(id=%lld) and Xid=%d(id=%lld),or=%d due to center cloning of Yid and/or Xid\n",
		   align->score,align->logPV,align->numpairs,align->mapid1,YYmap[align->mapid1]->id,align->mapid2,XXmap[align->mapid2]->id,align->orientation);
	    fflush(stdout);
	  }

	  LocalCnt--;
	  if(Yleft > 0 || Yright <= N){
	    if(DEBUG) assert(!Ymap->paired);
	    if(DEBUG>=2) assert(min(N,Yright)-max(1,Yleft) < N-1);/* make sure center cloned region is smaller than original */

	    if(PAIRSPLIT_EXTEND){/* limit extensions to no more than AlignedSiteThreshold-1 sites */
	      if(DEBUG && !(align->refstart >= 0)){
		printf("Yleft=%d,Yright=%d,N=%d:bestalign->refstart=%d,bestalign->refend=%d\n",Yleft,Yright,N,align->refstart,align->refend);
		fflush(stdout);
		assert(align->refstart >= 0);
	      }
	      if(DEBUG && !(align->refend <= N+1)){
		printf("Yleft=%d,Yright=%d,N=%d:align->refstart=%d,align->refend=%d\n",Yleft,Yright,N,align->refstart,align->refend);
		fflush(stdout);
		assert(align->refend <= N+1);
	      }
	      int origYleft = Yleft;
	      int origYright = Yright;

	      Yleft = min(Yleft,max(Yleft-AlignedSiteThreshold+1,align->refstart));
	      Yright = max(Yright,min(Yright+AlignedSiteThreshold-1,align->refend));

	      /* make sure the center cloned region is smaller than the original, to avoid infinite loops */
	      if(Yleft < origYleft)
		Yleft = max(2, Yleft);
	      if(Yright > origYright && N > 1)
		Yright = min(N-1,Yright);

	      if(LVERB){
		printf("Yleft=%d->%d,Yright=%d->%d,N=%d\n",origYleft,Yleft,origYright,Yright,N);
		fflush(stdout);
	      }
	    }
	    Cmap *newmap = NULL;

#pragma omp critical(numrefmaps)
	    {
	      Ymap->centercloned = numY;
	      maxmapalloc(numrefmaps+1,maxrefmaps,refmap,0,1);
	      YYmap = refmap;/* in case refmap[] was reallocated */
	      newmap = refmap[numrefmaps++];
	      if(!newmap)
		refmap[numrefmaps-1] = newmap = new Cmap;
	      newmap->mapid = numY;
	      //	    YYmap[numY++] = newmap;
	      if(DEBUG) assert(YYmap[numY] == newmap);
	      numY++;
	      if(DEBUG) assert(numY == numrefmaps);
	    }
	    newmap->allfree();
	    if(DEBUG) assert(newmap->blockmem==0);
	    newmap->centercloned = 0;
	    newmap->origmap = Ymap;
	    if(PairsplitSeparateQueries)
	      newmap->Xid = Ymap->Xid;
	    newmap->id = Ymap->id;
	    newmap->fileid = Ymap->fileid;
	    newmap->numsite[0] = min(N,Yright)-max(1,Yleft)+1;
	    newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
	    if(LEAKDEBUG>=2){
	      printf("refmap[%d]=%p:id=%lld,site[0]->%p,numsite[0]=%d\n",
		     numrefmaps-1,newmap,newmap->id,newmap->site[0],newmap->numsite[0]);
	      fflush(stdout);
	    }
	    newmap->site[0][0] = 0.0;
	    newmap->left[0] = Yleft;
	    newmap->right[0] = Yright;

	    double YLextend = MININTERVAL;
	    double YRextend = MININTERVAL;
	    if(OUTLIER_FIXES && align->Lend == -1 && Yleft > 0){/* make sure left end of new map is extended to at least overlap the left end of X */
	      if(DEBUG && !PAIRSPLIT_EXTEND && AlignedSiteThreshold >= 8) assert(Xleft <= 0);/* If Xleft > 0 && Yleft > 0, then both X and Y were trimmed on the left, hence there should be an endoutlier with align->Lend <= -2 */
	      int IL = align->sites1[0];
	      int JL = align->sites2[0];
	      double Yend = Y[IL] - (align->orientation ? X[M+1] - X[M+1-JL] : X[JL]) - 0.001;
	      YLextend = max(YLextend, Y[Yleft] - Yend);
	    }
	    if(OUTLIER_FIXES && align->Rend == -1 && Yright <= N){/* make sure right end is extended to at least overlap the right end of X */
	      if(DEBUG && !PAIRSPLIT_EXTEND && AlignedSiteThreshold >= 8 && !(Xright > M)){/* If Xright <= M && Yright <= N, then both X and Y were trimmed on the right, hence there should be an endoutlier with align->Rend  <= -2 */
		printf("Splitting off center of Yid=%d from Yleft=%d(%0.3f) to Yright=%d(%0.3f),N=%d(%0.3f): Xid=%d, Xleft=%d(%0.3f), Xright=%d(%0.3f,M=%d(%0.3f)\n",
		       Yid,Yleft,Y[Yleft],Yright,Y[Yright],N,Y[N],Xid,
		       Xleft, align->orientation ? X[M+1] - X[M+1-Xleft] : X[Xleft],
		       Xright, align->orientation ? X[M+1] - X[M+1-Xright] : X[Xright], 
		       M, align->orientation ? X[M+1] - X[1] : X[M]);
		printf("  align->Rend=%d, align->numpairs=K=%d,align->sites1[0,K-1]=%d,%d,align->sites2[0,K-1]=%d,%d,align->orientation=%d\n",
		       align->Rend,K,align->sites1[0],align->sites1[K-1],align->sites2[0],align->sites2[K-1],align->orientation);
		fflush(stdout);

		assert(Xright > M);
	      }
	      int IR = align->sites1[K - 1];
	      int JR = align->sites2[K - 1];
	      double Yend = Y[IR] + (align->orientation ? X[M+1-JR] : X[M+1] - X[JR]) + 0.001;
	      YRextend = max(YRextend, Yend - Y[Yright]);
	    }

	    if(Yleft == 0){/* copy left end up to Yright */
	      for(int i = 1; i <= Yright; i++)
		newmap->site[0][i] = Y[i];
	      newmap->site[0][Yright+1] = Y[Yright] + YRextend;
	    } else if(Yright == N+1){/* copy right end from Yleft (see above with Imax == Yleft-1)*/
	      for(int i = 1; i <= N+2 - Yleft; i++)
		newmap->site[0][i] = YLextend + Y[Yleft + i-1] - Y[Yleft];
	    } else {
	      for(int i = Yleft; i <= Yright; i++)
		newmap->site[0][i-Yleft+1] = YLextend + Y[i] - Y[Yleft];
	      newmap->site[0][newmap->numsite[0]+1] = YLextend + Y[Yright] - Y[Yleft] + YRextend;
	    }
	    newmap->startloc = Ymap->startloc;
	    newmap->endloc = Ymap->endloc;
	    newmap->len = Ymap->len * (newmap->site[0][newmap->numsite[0]+1]/Y[N+1]);
	    newmap->fpcnt = -1;

	    newmap->deltaXinit();

	    if(Ymap->origmap){
	      if(Ymap->left[0] > 0 && newmap->left[0] == 0)
		newmap->left[0] = 1;
	      if(Ymap->right[0] <= Ymap->origmap->numsite[0] && newmap->right[0] > N)
		newmap->right[0] = N;
	      newmap->origmap = Ymap->origmap;
	      if(DEBUG) assert(newmap->left[0] >= 0);
	      newmap->left[0] += max(0,Ymap->left[0] - 1);
	      newmap->right[0] += max(0,Ymap->left[0] - 1);
	      if(DEBUG) assert(newmap->origmap->origmap == 0);
	    }
	    if(LVERB){
	      printf("Split off center of Yid=%d from I=%d(%0.3f) to I=%d(%0.3f) to YYmap[%d]:id=%lld(file=%d),left=%d,right=%d,origN=%d,N=%d,YLextend=%0.3f,YRextend=%0.3f,\n",
		     Yid,Yleft,Y[Yleft],Yright,Y[Yright],numY-1,Ymap->id, Ymap->fileid,newmap->left[0],newmap->right[0],newmap->origmap->numsite[0],newmap->numsite[0],YLextend,YRextend);
	      fflush(stdout);
	    }
	    if(DEBUG>=2 && newmap->origmap && newmap->numsite[0] >= 2){/* check that mapping is consistent by checking the interval sizes */
	      int shift = max(1,newmap->left[0])-1;
	      assert(newmap->site[0][0] == 0.0);
	      assert(newmap->site[0][1] > 0.0);
	      assert(newmap->site[0][newmap->numsite[0]+1] > newmap->site[0][newmap->numsite[0]]);
	      for(int T = 1; T < newmap->numsite[0];T++){
		double newsize = newmap->site[0][T+1] - newmap->site[0][T];
		double oldsize = newmap->origmap->site[0][shift + T+1] - newmap->origmap->site[0][shift + T]; 
		if(fabs(newsize - oldsize) >= 0.001 || newsize < 0.0 || oldsize < 0.0){
		  for(int K = 1; K < newmap->numsite[0]; K++)
		    printf("YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f); YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f)\n",
			   numY-1,K,newmap->site[0][K],K+1,newmap->site[0][K+1],newmap->site[0][K+1]-newmap->site[0][K],
			   Yid,shift+K,newmap->origmap->site[0][shift+K],shift+K+1,newmap->origmap->site[0][shift+K+1],
			   newmap->origmap->site[0][shift+K+1]-newmap->origmap->site[0][shift+K]);
		  printf("YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f); YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f),shift=%d,newmap->left[0]=%d\n",
			 numY-1,T,newmap->site[0][T],T+1,newmap->site[0][T+1],newsize,
			 Yid,shift+T,newmap->origmap->site[0][shift+T],shift+T+1,newmap->origmap->site[0][shift+T+1],oldsize,shift,newmap->left[0]);
		  fflush(stdout);
		  assert(!(fabs(newsize-oldsize) >= 0.001));
		}
	      }
	    }
	  }
	  if(Xleft > 0 || Xright <= M){
	    if(DEBUG) assert(!Xmap->paired);

	    /* make sure center cloned region is smaller than the original */
	    if(DEBUG>=2 && !(min(M,Xright)-max(1,Xleft) < M-1)){
	      printf("Yid=%d,Xid=%d:Xleft=%d,Xright=%d,M=%d\n",Yid,Xid,Xleft,Xright,M);
	      fflush(stdout);
	      assert(Xright-Xleft < M-1);
	    }

	    if(PAIRSPLIT_EXTEND){/* limit extensions to no more than AlignedSiteThreshold-1 sites */
	      if(DEBUG) assert(align->start >= 0);
	      if(DEBUG) assert(align->end <= M+1);
	      int origXleft = Xleft;
	      int origXright = Xright;
	      Xleft = min(Xleft,max(Xleft-AlignedSiteThreshold+1,align->start));
	      Xright = max(Xright,min(Xright+AlignedSiteThreshold-1,align->end));

	      /* make sure the extended center cloned region is smaller than the original, to avoid infinite loops */
	      if(Xleft < origXleft)
		Xleft = max(2, Xleft);
	      if(Xright > origXright && M > 1)
		Xright = min(M-1,Xright);

	      if(LVERB){
		printf("Xleft=%d->%d,Xright=%d->%d,M=%d\n",origXleft,Xleft,origXright,Xright,M);
		fflush(stdout);
	      }
	    }
	    Cmap *newmap = NULL;

#pragma omp critical(nummaps)
	    {
	      Xmap->centercloned = numX;
	      maxmapalloc(nummaps+1,maxmaps,Gmap,0,1);
	      XXmap = &Gmap[Cfirst];/* in case Gmap[] was reallocated */
	      newmap = Gmap[nummaps++];
	      if(!newmap)
		Gmap[nummaps-1] = newmap = new Cmap;
	      newmap->mapid = numX;
	      //	    XXmap[numX++] = newmap;
	      if(DEBUG) assert(XXmap[numX] == newmap);
	      numX++;
	      if(DEBUG && !(numX == nummaps - Cfirst)){
		printf("numX=%d,nummaps=%d,mapstart=%d,Cfirst=%d,%d\n",numX,nummaps,0,Cfirst,CfirstPairs);
		fflush(stdout);
		assert(numX == nummaps - Cfirst);
	      }
	    }
	    newmap->allfree();
	    if(DEBUG) assert(newmap->blockmem==0);
	    newmap->centercloned = 0;
	    newmap->origmap = Xmap;
	    if(PairsplitSeparateQueries)
	      newmap->Xid = Xmap->Xid;
	    newmap->id = Xmap->id;
	    newmap->fileid = Xmap->fileid;
	    newmap->numsite[0] = min(M,Xright)-max(1,Xleft)+1;
	    newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
	    if(LEAKDEBUG>=2){
	      printf("map[%d]=%p:id=%lld,site[0]->%p,numsite[0]=%d\n",
		     nummaps-1,newmap,newmap->id,newmap->site[0],newmap->numsite[0]);
	      fflush(stdout);
	    }
	    newmap->site[0][0] = 0.0;
	    double XLextend = MININTERVAL;
	    double XRextend = MININTERVAL;
	    if(!align->orientation){
	      newmap->left[0] = Xleft;
	      newmap->right[0] = Xright;

	      if(OUTLIER_FIXES && align->Lend == -1 && Xleft > 0){/* make sure left end of new map is extended to at least overlap the left end of Y */
		if(DEBUG && !PAIRSPLIT_EXTEND && AlignedSiteThreshold >= 8) assert(Yleft <= 0);/* If Yleft > 0 && Xleft > 0, then both X and Y were trimmed on the left, hence there should be an endoutlier with align->Lend <= -2 */
		int IL = align->sites1[0];
		int JL = align->sites2[0];
		double Xend = X[JL] - Y[IL] - 0.001;
		XLextend = max(XLextend, X[Xleft] - Xend);
	      }
	      if(OUTLIER_FIXES && align->Rend == -1 && Xright <= M){/* make sure right end of new map is extended to at least overlap the right end of Y */
		if(DEBUG && !PAIRSPLIT_EXTEND && AlignedSiteThreshold >= 8) assert(Yright > N);/* If Yright <= N && Xright <= M, then both X and Y were trimmed on the right, hence there should be an endoutlier with align->Rend <= -2 */
		int IR = align->sites1[K - 1];
		int JR = align->sites2[K - 1];
		double Xend = X[JR] + Y[N+1] - Y[IR] + 0.001;
		XRextend = max(XRextend, Xend - X[Xright]);
	      }
	      if(Xleft == 0){/* copy left end only (see above with J-1 == Xright) */
		for(int i = 1; i <= Xright; i++)
		  newmap->site[0][i] = X[i];
		newmap->site[0][Xright+1] = newmap->site[0][Xright] + XRextend;
	      } else if(Xright == M+1){/* copy right end only (see above with Jmax == Xleft - 1) */
		for(int i = 1; i <= M+2 - Xleft; i++)
		  newmap->site[0][i] = XLextend + X[newmap->left[0] + i-1] - X[newmap->left[0]];
	      } else {
		for(int i = Xleft; i <= Xright; i++)
		  newmap->site[0][i-Xleft+1] = XLextend + X[i] - X[Xleft];
		newmap->site[0][newmap->numsite[0]+1] = XLextend + X[Xright] - X[Xleft] + XRextend;
	      }
	    } else {
	      newmap->right[0] = M+1 - Xleft;
	      newmap->left[0] = M+1 - Xright;
	      if(DEBUG) assert(newmap->right[0] >= 0);
	      if(DEBUG) assert(newmap->left[0] >= 0);
	  
	      if(OUTLIER_FIXES && align->Lend == -1 && Xleft > 0){/* make sure left end of new map (reversed) is extended to at least overlap the left end of Y */
		if(DEBUG && !PAIRSPLIT_EXTEND && AlignedSiteThreshold >= 8) assert(Yleft <= 0);/* If Yleft > 0 && Xleft > 0, then both X(reversed) and Y were trimmed on the left, hence there should be an endoutlier with align->Lend <= -2 */
		int IL = align->sites1[0];
		int JL = align->sites2[0];
		double Xend = (X[M+1] - X[M+1 - JL]) - Y[IL] - 0.001;
		XLextend = max(XLextend, X[M+1] - X[M+1 - Xleft] - Xend);
	      }
	      if(OUTLIER_FIXES && align->Rend == -1 && Xright <= M){/* make sure right end of new map (reversed) is extended to at least overlap the right end of Y */
		if(DEBUG && !PAIRSPLIT_EXTEND && AlignedSiteThreshold >= 8) assert(Yright > N);/* If Yright <= N && Xright <= M, then both X(reversed) and Y were trimmed on the right, hence there should be an endoutlier with align->Rend <= -2 */
		int IR = align->sites1[K - 1];
		int JR = align->sites2[K - 1];
		double Xend = X[M+1] - X[M+1 - JR] + Y[N+1] - Y[IR] + 0.001;
		XRextend = max(XRextend, Xend - (X[M+1] - X[M+1 - Xright]));
	      }

	      if(Xleft == 0){/* copy left end only of reversed X[] (OR right end of X[]) */
		for(int i = 1; i <= Xright+1; i++)
		  newmap->site[0][i] = XRextend + X[newmap->left[0] + i-1] - X[newmap->left[0]];
	      } else if(Xright == M+1){/* copy right end only of reversed X[] (or left end of X[], see above with Jmax == Xleft - 1) */
		for(int i = 1; i <= M+1-Xleft; i++)
		  newmap->site[0][i] = X[i];
		newmap->site[0][M+2-Xleft] = XLextend + newmap->site[0][M+1-Xleft];
	      } else {
		for(int i = Xleft; i <= Xright; i++)
		  newmap->site[0][Xright+1-i] = XRextend + X[M+1-i] - X[newmap->left[0]];
		newmap->site[0][newmap->numsite[0]+1] = XRextend + X[newmap->right[0]] - X[newmap->left[0]] + XLextend;
	      }
	    }

	    //	printf("newmap->left[0]=%d,newmap->right[0]=%d,newmap->numsite=%d\n",newmap->left[0],newmap->right[0],newmap->numsite[0]);

	    newmap->startloc = Xmap->startloc;
	    newmap->endloc = Xmap->endloc;
	    newmap->len = Xmap->len * (newmap->site[0][newmap->numsite[0]+1]/X[M+1]);
	    newmap->fpcnt = -1;

	    newmap->deltaXinit();

	    if(Xmap->origmap){
	      if(Xmap->left[0] > 0 && newmap->left[0] == 0)
		newmap->left[0] = 1;
	      if(Xmap->right[0] <= Xmap->origmap->numsite[0] && newmap->right[0] > M)
		newmap->right[0] = M;

	      newmap->origmap = Xmap->origmap;
	      if(DEBUG) assert(newmap->left[0] >= 0);
	      newmap->left[0] += max(0, Xmap->left[0] - 1);
	      newmap->right[0] += max(0, Xmap->left[0] - 1);
	      if(DEBUG) assert(newmap->origmap->origmap == 0);
	    }

	    if(LVERB){
	      printf("Split off center of Xid=%d(flip=%d) from J=%d(%0.3f) to %d(%0.3f) to XXmap[%d]:id=%lld(file=%d),left=%d,right=%d,origM=%d,M=%d,XLextend=%0.3f,XRextend=%0.3f\n",
		     Xid,align->orientation,Xleft,align->orientation ? X[M+1]-X[M+1-Xleft] : X[Xleft],
		     Xright,align->orientation ? X[M+1]-X[M+1-Xright] : X[Xright],numX-1,Xmap->id,Xmap->fileid,newmap->left[0],newmap->right[0],newmap->origmap->numsite[0],newmap->numsite[0],XLextend,XRextend);
	      fflush(stdout);
	    }

	    if(DEBUG>=2 && newmap->origmap && newmap->numsite[0] >= 2){/* check that mapping is consistent by checking the interval sizes */
	      assert(newmap->site[0][0] == 0.0);
	      assert(newmap->site[0][1] > 0.0);
	      assert(newmap->site[0][newmap->numsite[0]+1] > newmap->site[0][newmap->numsite[0]]);
	      for(int T = 1; T < newmap->numsite[0];T++){
		int shift = max(1,newmap->left[0])-1;
		double newsize = newmap->site[0][T+1] - newmap->site[0][T];
		double oldsize = newmap->origmap->site[0][shift + T+1] - newmap->origmap->site[0][shift + T]; 
		if(fabs(newsize - oldsize) >= 0.001 || newsize < 0.0 || oldsize < 0.0){
		  printf("XXmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f); YYmap[%d]->site[%d]=%0.4f,site[%d]=%0.4f (delta=%0.4f)\n",
			 numX-1,T,newmap->site[0][T],T+1,newmap->site[0][T+1],newsize,
			 Xid,shift+T,newmap->origmap->site[0][shift+T],shift+T+1,newmap->origmap->site[0][shift+T+1],oldsize);
		  fflush(stdout);
		  assert(!(fabs(newsize-oldsize) >= 0.001));
		}
	      }
	    }
	  }
	}
      }// for(rootXid = 0; rootXid < orignumX; rootXid++)
      if(OMP_DEBUG && myLoopCnt > 0){
        #pragma omp critical
	{
	  LoopCnt += myLoopCnt;
        }
      }
    }
    if(OMP_DEBUG && rootXidInc == 1) assert(LoopCnt == orignumX);// verify OMP loop works

    if(VERB>=2 && PairSplit){
      printf("After splitting maps of %lu best alignment for %d original Query id's: numX=%d,numY=%d : cum wall time=%0.6f\n", splitcnt, orignumX,numX,numY,wtime());
      fflush(stdout);
    }
    if(!splitcnt)
      break;
  }// while(1)

  if(PairSplit){
    /* remove alignments with chimpair=1 (these were previously filtered out) */
    size_t U = 0;
    for(size_t A = 0; A < numaligns; A++){
      Calign *align = alignment[A];
      if(!align)
	continue;
      alignment[A] = 0;
      if(DEBUG) assert(align->chimpair != 0);
      if(align->chimpair==1){
	align->allfree();
	delete [] align;
	continue;
      }
      if(DEBUG) assert(alignment[U] == 0);
      alignment[U++] = align;
    }
    numaligns = U;
    if(VERB){
      printf("Final number of (above threshold) local alignments = %llu\n",(unsigned long long)numaligns);
      fflush(stdout);
    }

    if(PairsplitMerge){/* try to merge neighboring alignments */

      FLOAT *origXrevMem = new FLOAT[maxM + 2];
      int localtype = (PoutlierEnd > 0.0) ? -3 : -2;

      if(VERB >= (LVERB?1:2)){
	printf("-pairsplitMerge: examining %lu alignments\n",numaligns);
	for(size_t i = 0; i < numaligns; i++){
	  Calign *p = alignment[i];
	  printf("alignment[%lu]=%p:", i,p);
	  if(!p){
	    printf("\n");
	    continue;
	  }
	  int rid = p->mapid1;
	  int mid = p->mapid2;
	  printf("rid=%d(%lld),mid=%d(%lld),or=%d\n",rid,YYmap[rid]->id,mid,XXmap[mid]->id,p->orientation);
	}
	fflush(stdout);
      }

      /* sort alignments in order of id2, id1, orientation, sites1[0],sites1[U-1],sites2[0],sites2[U-1] (null alignments last) : this way neighboring alignments will be next to each other */
      qsort(&alignment[0],numaligns,sizeof(Calign *), (intcmp*)IdOrSiteInc);
      
      double threshold = -log(1.0/Psplit);

      for(size_t i = 1; i < numaligns; i++){
       Calign *p = alignment[i-1];
       for(size_t i2 = i; i2 <= numaligns; i2++){
	Calign *q = alignment[i2];
	if(VERB >= 2 /* (LVERB?1:2)*/){
	  printf("i=%lu,i2=%lu:p=alignment[i-1]=%p, q=alignment[i2]=%p\n",i,i2,p,q);
	  fflush(stdout);
	}
	if(DEBUG < 2 && i2 > i) break; // NOTE : comment out this line (or set DEBUG>=2) to check all pairs of alignments (to debug sorting order of alignments)

	if(!p || !q)
	  continue;

	int rid = p->mapid1;
	int mid = p->mapid2;
	int rid2 = q->mapid1;
	int mid2 = q->mapid2;

	if(VERB >= (LVERB?1:2)){
	  printf("i=%lu,i2=%lu: rid=%d(%lld),rid2=%d(%lld),mid=%d(%lld),mid2=%d(%lld),or=%d,or2=%d: checking pairsplitMerge\n",
            i, i2,rid, YYmap[rid]->id, rid2, YYmap[rid2]->id, mid, XXmap[mid]->id, mid2, XXmap[mid2]->id, p->orientation, q->orientation);
	  fflush(stdout);
	}

	Cmap *Ymap = YYmap[rid];
	Cmap *Xmap = XXmap[mid];
	Cmap *Ymap2 = YYmap[rid2];
	Cmap *Xmap2 = XXmap[mid2];

	if(Ymap->id != Ymap2->id || Xmap->id != Xmap2->id || p->orientation != q->orientation)
	  continue;
	
	int U = p->numpairs;
	/*	Cmap *origmap = Xmap;
	while(origmap->origmap){
	  Xmap->origmap = origmap = origmap->origmap;
	  Xmap->left[0] += max(0,origmap->left[0]-1);
	  Xmap->right[0] += max(0,origmap->left[0]-1);
	  } */
	int M = Xmap->numsite[0];
	int N = Ymap->numsite[0];
	//	FLOAT *Y = Ymap->site[0];
	//	FLOAT *X = Xmap->site[0];
	if(DEBUG) assert(Xmap->origmap==0 || Xmap->origmap->origmap==0);/* If this fails fix origmap link in pairalign() to point to root map */
	if(DEBUG) assert(Ymap->origmap==0 || Ymap->origmap->origmap==0);/* If this fails fix origmap link (same as for Xmap above) to point to root map */
	Cmap *origXmap = Xmap->origmap ? Xmap->origmap : Xmap;
	Cmap *origYmap = Ymap->origmap ? Ymap->origmap : Ymap;
	FLOAT *origX = origXmap->site[0];
	FLOAT *origY = origYmap->site[0];
	int Yshift = Ymap->origmap ? max(1,Ymap->left[0]) - 1 : 0;
	int Xshift = Xmap->origmap ? max(1,Xmap->left[0]) - 1 : 0;
	int LI = p->sites1[0];
	int LJ = p->sites2[0];
	int RI = p->sites1[U-1];
	int RJ = p->sites2[U-1];
	int refstart = Yshift + LI;
	int refend = Yshift + RI;
	int qrystart = p->orientation ? Xshift + M + 1 - LJ : Xshift + LJ;
	int qryend   = p->orientation ? Xshift + M + 1 - RJ : Xshift + RJ;

	int M2 = Xmap2->numsite[0];
	int N2 = Ymap2->numsite[0];
	//	FLOAT *Y2 = Ymap2->site[0];
	//	FLOAT *X2 = Xmap2->site[0];
	if(DEBUG) assert(Xmap2->origmap==0 || Xmap2->origmap->origmap==0);/* If this fails fix origmap link in pairalign() to point to root map */
	if(DEBUG) assert(Ymap2->origmap==0 || Ymap2->origmap->origmap==0);/* If this fails fix origmap link (same as for Xmap above) to point to root map */
	Cmap *origXmap2 = Xmap2->origmap ? Xmap2->origmap : Xmap2;
	Cmap *origYmap2 = Ymap2->origmap ? Ymap2->origmap : Ymap2;
	FLOAT *origX2 = origXmap2->site[0];
	FLOAT *origY2 = origYmap2->site[0];
	int Yshift2 = Ymap2->origmap ? max(1,Ymap2->left[0]) - 1 : 0;
	int Xshift2 = Xmap2->origmap ? max(1,Xmap2->left[0]) - 1 : 0;

	int U2 = q->numpairs;
	int LI2 = q->sites1[0];
	int LJ2 = q->sites2[0];
	int RI2 = q->sites1[U2-1];
	int RJ2 = q->sites2[U2-1];
	int refstart2 = Yshift2 + LI2;
	int refend2 = Yshift2 + RI2;
	int qrystart2 = q->orientation ? Xshift2 + M2 + 1 - LJ2 : Xshift2 + LJ2;
	int qryend2   = q->orientation ? Xshift2 + M2 + 1 - RJ2 : Xshift2 + RJ2;

	if(DEBUG) assert(origXmap2 == origXmap);
	if(DEBUG) assert(origYmap2 == origYmap);
	int origM = origXmap->numsite[0];
	int origN = origYmap->numsite[0];

	double *origXrev = origX;
	if(p->orientation){
	  origXrev = origXrevMem;
	  for(int J = 0; J <= origM + 1; J++)
	    origXrev[J] = origX[origM + 1] - origX[origM + 1 -J];
	}

	int mapid1 = origYmap->mapid;
	int mapid2 = origXmap->mapid;
	if(DEBUG && !(0 <= mapid1 && mapid1 < numrefmaps && YYmap[mapid1] == origYmap)){
	  if(0 <= mapid1 && mapid1 < numY)
	    printf("origYmap=%p,origYmap->mapid=mapid1=%d,YYmap[mapid1]=%p,YYmap[mapid1]->mapid=%d(id=%lld)\n",
		   origYmap,origYmap->mapid,YYmap[mapid1],YYmap[mapid1]->mapid,YYmap[mapid1]->id);
	  else
	    printf("origYmap=%p,origYmap->mapid=mapid1=%d,numY=%d\n",
		   origYmap,origYmap->mapid,numY);
	  assert(0 <= mapid1 && mapid1 < numY && YYmap[mapid1] == origYmap);
	}
	if(DEBUG && !(0 <= mapid2 && mapid2 < numX && XXmap[mapid2] == origXmap)){
	  if(0 <= mapid2 && mapid2 < numX)
	    printf("origXmap=%p,origXmap->mapid=mapid2=%d,XXmap[mapid2]=%p,XXmap[mapid2]->mapid=%d(id=%lld)\n",
		   origXmap,origXmap->mapid,XXmap[mapid2],XXmap[mapid2]->mapid,XXmap[mapid2]->id);
	  else
	    printf("origXmap=%p,origXmap->mapid=mapid2=%d,numX=%d\n",
		   origXmap,origXmap->mapid,numX);
	  fflush(stdout);
	  assert(0 <= mapid2 && mapid2 < numX && XXmap[mapid2] == origXmap);
	}

	PFLOAT Bias,Pen;/* For use as return values by SintDetail() */

	if(VERB >= (LVERB?1:2)){
	  printf("Yid=%lld,Xid=%lld,or=%d:N=%d,M=%d,N2=%d,M2=%d:LI=%d(%d,%0.3f),RI=%d(%d,%0.3f),LJ=%d(%d,%0.3f),RJ=%d(%d,%0.3f):\n\t LI2=%d(%d,%0.3f),RI2=%d(%d,%0.3f),LJ2=%d(%d,%0.3f),RJ2=%d(%d,%0.3f)\n",
		 Ymap->id,Xmap->id,p->orientation,N,M,N2,M2,LI,refstart,origY[refstart],RI,refend,origY[refend],LJ,qrystart,origX[qrystart],RJ,qryend,origX[qryend],
                 LI2,refstart2,origY2[refstart2],RI2,refend2,origY2[refend2],LJ2,qrystart2,origX2[qrystart2],RJ2,qryend2,origX2[qryend2]);
	  fflush(stdout);
	}

	if(refend < refstart2){/* p is left of q */
	  if(VERB >= (LVERB ? 1 : 2)){
	    printf("  Neighboring alignments with gap\n");
	    fflush(stdout);
	  }
	  if(q->orientation ? (qrystart2 >= qryend) : (qrystart2 <= qryend)){
	    if(VERB >= (LVERB ? 1 : 2)){            
              printf("    Crossed alignment (query has -ve gap size)\n");
	      fflush(stdout);
            }
	    continue;/* crossed alignment */
          }

	  /* check if the gap is a valid alignment interval (NOT an outlier) */
	  PFLOAT x = fabs(origX2[qrystart2]-origX[qryend]);
	  PFLOAT y = origY2[refstart2]-origY[refend];
	  SintDetail(x, y, abs(qrystart2-qryend), refstart2-refend, qrystart2, refstart2, Bias, Pen);
	  PFLOAT outscore = OutlierBias + Bias + Pen;
	  PFLOAT nthreshold = threshold;
	  PFLOAT OutPen = OutlierPenalty;
	  PFLOAT BCscore = OutlierPenaltyBC[refstart2-refend] + OutlierPenaltyBC[abs(qrystart2-qryend)]; 
	  if(outlierBC) {
	    nthreshold += BCscore;
	    OutPen += BCscore;
	  }
	  PFLOAT iscore = OutlierBias + (OUTLIER_DELTA(x-y) ? max(Pen + Bias, OutPen + Bias*biasWToutlierF) : Pen + Bias);
	  if(outscore < nthreshold){ /* An outlier at -pairsplit P-value */
	    if(VERB >= (LVERB ? 1 : 2)){            
              printf("    Psplit gap: nthreshold=%0.6f,iscore=%0.6f,outscore=%0.6f, x=%0.4f,y=%0.4f,m=%d,n=%d,Pen=%0.6f,Bias=%0.6f\n", nthreshold, iscore, outscore, x, y, abs(qrystart2-qryend), refstart2-refend, Pen, Bias);
              fflush(stdout);
            }

	    continue;
          }

	  /* merge alignment q and p */
	  Calign *r = new Calign[colors==1 ? 1 : colors+1];
	  r->expand_arrays(U+U2);
	  r->numpairs = U + U2;
	  r->mapid1 = mapid1;
	  r->mapid2 = mapid2;
	  r->orientation = p->orientation;
	  r->scaleID = p->scaleID;
	  if(DEBUG) assert(YYmap[mapid1]->numsite[0] == origN);
	  if(DEBUG) assert(XXmap[mapid2]->numsite[0] == origM);

	  if(VERB >= (LVERB ? 1 : 2)){
	    printf("p is left of q and gap is a valid alignment interval (iscore=%0.6f,outscore=%0.6f) : combining into one alignment (p->numpairs=%d,q->numpairs=%d,r->numpairs=%d)\n",
		   iscore,outscore,p->numpairs,q->numpairs,r->numpairs);
	    fflush(stdout);
	  }
          if(DEBUG)assert(i2 == i);

	  /* copy p */
	  for(int t = 0; t < U; t++){
	    r->sites1[t] = Yshift + p->sites1[t];
	    r->sites2[t] = (p->orientation ? origM - Xshift - M : Xshift) + p->sites2[t];
	    r->iscore[t] = p->iscore[t];
	    r->outscore[t] = p->outscore[t];
	  }

	  /* copy iscore,outscore of gap region */
	  r->iscore[U] = iscore;
	  r->outscore[U] = outscore;

	  /* copy q */
	  for(int t = 0; t < U2; t++){
	    r->sites1[t+U] = Yshift2 + q->sites1[t];
	    r->sites2[t+U] = (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->sites2[t];
	    r->iscore[t+U+1] = q->iscore[t+1];
	    r->outscore[t+U+1] = q->outscore[t+1];
	  }
	  
	  if(DEBUG/* HERE >=2 */){
	    for(int t = 1; t < r->numpairs; t++){
	      assert(r->sites1[t] > r->sites1[t-1]);
	      assert(r->sites2[t] > r->sites2[t-1]);
	    }
	  }

	  if(Xshift || Yshift){
	    r->Lend = localtype;
	    r->Lij1 = r->sites1[0];
	    r->Lij2 = r->sites2[0];
	  } else {
	    r->Lend = (Xshift || Yshift) ? localtype : p->Lend;
	    r->Lij1 = Yshift + p->Lij1;
	    r->Lij2 = Xshift + p->Lij2; // WAS (p->orientation ? origM - Xshift - M : Xshift) + p->Lij2;
	  }
	  if(Xshift2 + M2 < origM || Yshift2 + N2 < origN){
	    r->Rend = localtype;
	    r->Rij1 = r->sites1[r->numpairs - 1];
	    r->Rij2 = r->sites2[r->numpairs - 1];
	  } else {
	    r->Rend = q->Rend;
	    r->Rij1 = Yshift2 + q->Rij1;
	    r->Rij2 = Xshift2 + q->Rij2; // WAS (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->Rij2;
	  }
	  r->score = 0;
	  for(int t = 0; t <= U + U2; t++)
	    r->score += r->iscore[t];
	  r->logPV = alignFP(r, origY, origXrev, origN, origM, p->orientation, p->scaleID, mapid1, mapid2, r->score, 0);
	  //	  r->logPV = p->logPV + q->logPV;

	  delete [] p;
	  delete [] q;
	  p = alignment[i-1] = r;
	  q = alignment[i2] = 0;
	} else if(refend2 < refstart){/* p is right of q */
	  if(VERB >= (LVERB ? 1 : 2)){
	    printf("  Neighboring alignments with gap\n");
	    fflush(stdout);
	  }
	  if(q->orientation ? (qrystart >= qryend2) : (qrystart <= qryend2)){
	    if(VERB >= (LVERB ? 1 : 2)){            
              printf("    Crossed alignment (query has -ve gap size)\n");
	      fflush(stdout);
            }
	    continue;/* crossed alignment */
	  }

	  /* check if the gap is a valid alignment interval (NOT an outlier) */
	  PFLOAT x = fabs(origX[qrystart]-origX2[qryend2]);
	  PFLOAT y =  origY[refstart]-origY2[refend2];
	  SintDetail(x, y, abs(qrystart - qryend2), refstart - refend2, qrystart, refstart, Bias, Pen);
	  FLOAT outscore = OutlierBias + Bias + Pen;
	  PFLOAT nthreshold = threshold;
	  PFLOAT OutPen = OutlierPenalty;
	  PFLOAT BCscore = OutlierPenaltyBC[refstart - refend2] + OutlierPenaltyBC[abs(qrystart - qryend2)];
	  if(outlierBC){
	    nthreshold += BCscore;
	    OutPen += BCscore;
	  }
	  FLOAT iscore = OutlierBias + (OUTLIER_DELTA(x-y) ? max(Pen + Bias, OutPen + Bias*biasWToutlierF) : Pen + Bias);
	  if(outscore < nthreshold){ /* An outlier at -pairsplit P-value */
	    if(VERB >= (LVERB ? 1 : 2)){            
              printf("    Psplit gap: nthreshold=%0.6f,iscore=%0.6f,outscore=%0.6f, x=%0.4f,y=%0.4f,m=%d,n=%d,Pen=%0.6f,Bias=%0.6f\n", 
		     nthreshold, iscore, outscore, x, y, abs(qrystart-qryend2), refstart-refend2, Pen, Bias);
              fflush(stdout);
            }

	    continue;
	  }

	  /* merge alignment q and p */	  
	  Calign *r = new Calign[colors==1 ? 1 : colors+1];
	  r->expand_arrays(U+U2);
	  r->numpairs = U + U2;
	  r->mapid1 = mapid1;
	  r->mapid2 = mapid2;
	  r->orientation = p->orientation;

	  if(VERB >= (LVERB ? 1 : 2)){
	    printf("p is right of q and gap is a valid alignment interval (iscore=%0.6f,outscore=%0.6f) : combining into one alignment (p->numpairs=%d,q->numpairs=%d,r->numpairs=%d)\n",
		   iscore,outscore,p->numpairs,q->numpairs,r->numpairs);
	    fflush(stdout);
	  }
          if(DEBUG) assert(i2 == i);

	  /* copy q */
	  for(int t = 0; t < U2; t++){
	    r->sites1[t] = Yshift2 + q->sites1[t];
	    r->sites2[t] = (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->sites2[t];
	    r->iscore[t] = q->iscore[t];
	    r->outscore[t] = q->outscore[t];
	  }

	  /* copy iscore,outscore of gap region */
	  r->iscore[U2] = iscore;
	  r->outscore[U2] = outscore;

	  /* copy p */
	  for(int t = 0; t < U; t++){
	    r->sites1[t+U2] = Yshift + p->sites1[t];
	    r->sites2[t+U2] = (p->orientation ? origM - Xshift - M : Xshift) + p->sites2[t];
	    r->iscore[t+U2+1] = p->iscore[t+1];
	    r->outscore[t+U2+1] = p->outscore[t+1];
	  }

	  if(DEBUG/* HERE >=2 */){
	    for(int t = 1; t < r->numpairs; t++){
	      assert(r->sites1[t] > r->sites1[t-1]);
	      assert(r->sites2[t] > r->sites2[t-1]);
	    }
	  }

	  if(Xshift2 || Yshift2){
	    r->Lend = localtype;
	    r->Lij1 = r->sites1[0];
	    r->Lij2 = r->sites2[0];
	  } else {
	    r->Lend = q->Lend;
	    r->Lij1 = Yshift2 + q->Lij1;
	    r->Lij2 = Xshift2 + q->Lij2; // WAS (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->Lij2;
	  }
	  if(Xshift + M < origM || Yshift + N < origN){
	    r->Rend = localtype;
	    r->Rij1 = r->sites1[r->numpairs - 1];
	    r->Rij2 = r->sites2[r->numpairs - 1];
	  } else{
	    r->Rend = p->Rend;
	    r->Rij1 = Yshift + p->Rij1;
	    r->Rij2 = Xshift + p->Rij2; // WAS (p->orientation ? origM - Xshift - M : Xshift) + p->Rij2;
	  }
	  r->score = 0;
	  for(int t = 0; t <= U + U2; t++)
	    r->score += r->iscore[t];
	  r->logPV = alignFP(r, origY, origXrev, origN, origM, p->orientation, p->scaleID, mapid1, mapid2, r->score, 0);
	  //	  r->logPV = p->logPV + q->logPV;

	  delete [] p;
	  delete [] q;
	  p = alignment[i-1] = r;
	  q = alignment[i2] = 0;

	} else {/* overlapped or adjoining with shared site in Y */
	  if(VERB >= (LVERB ? 1 : 2)){
	    printf("  Overlapped or Adjoining alignments\n");
	    fflush(stdout);
	  }
	  /* first check if one alignment completely overlaps the other */
	  if(refstart <= refstart2 && refend2 <= refend  &&
	     (p->orientation ? (qryend <= qryend2 && qrystart2 <= qrystart) : (qrystart <= qrystart2 && qryend2 <= qryend))){/* p completely overlaps q */
	    if(VERB >= (LVERB ? 1 : 2)){
	      printf("Complete Overlap of q (deleting)\n");
	      fflush(stdout);
	    }
            if(DEBUG)assert(i2 == i);
	    delete [] q;
	    q = alignment[i2] = 0;
	  } else if(refstart2 <= refstart && refend <= refend2 &&
		    (p->orientation ? (qryend2 <= qryend && qrystart <= qrystart2) : (qrystart2 <= qrystart && qryend <= qryend2))){/* q completely overlaps by p */
	    if(VERB >= (LVERB ? 1 : 2)){
	      printf("Complete Overlap of p (deleting)\n");
	      fflush(stdout);
	    }
            if(DEBUG)assert(i2 <= i);
	    delete [] p;
	    p = alignment[i-1] = 0;
	  } else if(refstart < refstart2 && refend < refend2 && 
		    (p->orientation ? (qrystart > qrystart2 && qryend > qryend2) : (qrystart < qrystart2 && qryend < qryend2))){/* p left of q */
	    if(VERB >= (LVERB ? 1 : 2)){
	      printf("p overlaps q from left\n");
	      fflush(stdout);
	    }
	    /* check if overlapped alignment region matches */
	    int t1 = 0, t2 = 0;
	    /* scan ahead in Yshift + p->sites1[t1] until it matches refstart2 */
	    while(Yshift + p->sites1[t1] != refstart2 && t1 < p->numpairs)
	      t1++;
	    if(t1 >= p->numpairs)
	      continue;/* no match found */
	    if(DEBUG) assert(Yshift + p->sites1[t1] == Yshift2 + q->sites1[t2]);
	    if(p->orientation ? (Xshift + M - p->sites2[t1] != Xshift2 + M2 - q->sites2[t2]) : (Xshift + p->sites2[t1] != Xshift2 + q->sites2[t2]))
	      continue;/* not a match */
	    int origt1 = t1;
	    for(;++t1 < p->numpairs && ++t2 < q->numpairs;)
	      if(Yshift + p->sites1[t1] != Yshift2 + q->sites1[t2] || 
		 (p->orientation ? (Xshift + M - p->sites2[t1] != Xshift2 + M2 - q->sites2[t2]) : (Xshift + p->sites2[t1] != Xshift2 + q->sites2[t2])))
		break;
	    if(t1 < p->numpairs && t2 < q->numpairs)
	      continue;/* not a complete match */

	    t1 = origt1;
	    if(DEBUG) assert(t1 > 0);

	    /* merge p (up to sites t1-1) with q */
	    Calign *r = new Calign[colors==1 ? 1 : colors+1];
	    r->numpairs = t1 + q->numpairs;
	    r->expand_arrays(r->numpairs);
	    r->mapid1 = mapid1;
	    r->mapid2 = mapid2;
	    r->orientation = p->orientation;

	    if(VERB >= (LVERB ? 1 : 2)){
	      printf("p overlaps q from left and overlapped region matches : combining into one alignment (t1=%d,p->numpairs=%d,q->numpairs=%d,r->numpairs=%d)\n",
		     t1, p->numpairs,q->numpairs,r->numpairs);
	      printf("p->Lij=%d,%d,p->Lend=%d,p->Rij=%d,%d,p->Rend=%d,q->Lij=%d,%d,q->Lend=%d,q->Rij=%d,%d,q->Rend=%d\n",
		     p->Lij1,p->Lij2,p->Lend,p->Rij1,p->Rij2,p->Rend,q->Lij1,q->Lij2,q->Lend,q->Rij1,q->Rij2,q->Rend);
	      fflush(stdout);
	    }
            if(DEBUG)assert(i2 <= i);

	    /* copy p (sites 0..t1-1) */
	    for(int t = 0; t < t1; t++){
	      r->sites1[t] = Yshift + p->sites1[t];
	      r->sites2[t] = (p->orientation ? origM - Xshift - M : Xshift) + p->sites2[t];
	      r->iscore[t] = p->iscore[t];
	      r->outscore[t] = p->outscore[t];
	    }
	    /* copy q */
	    for(int t = 0; t < q->numpairs; t++){
	      r->sites1[t1 + t] = Yshift2 + q->sites1[t];
	      r->sites2[t1 + t] = (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->sites2[t];
	      r->iscore[t1 + t] = q->iscore[t];
	      r->outscore[t1 + t] = q->outscore[t];
	    }
	    r->iscore[r->numpairs] = q->iscore[q->numpairs];
	    r->outscore[r->numpairs] = q->outscore[q->numpairs];

	    if(DEBUG/* HERE >=2 */){
	      for(int t = 1; t < r->numpairs; t++){
		assert(r->sites1[t] > r->sites1[t-1]);
		assert(r->sites2[t] > r->sites2[t-1]);
	      }
	    }

	    if(Xshift || Yshift){
	      r->Lend = localtype;
	      r->Lij1 = r->sites1[0];
	      r->Lij2 = r->sites2[0];
	    } else {
	      r->Lend = p->Lend;
	      r->Lij1 = Yshift + p->Lij1;
	      r->Lij2 = Xshift + p->Lij2; // WAS (p->orientation ? origM - Xshift - M : Xshift) + p->Lij2;
	    }
	    if(Xshift2 + M2 < origM || Yshift2 + N2 < origN){
	      r->Rend = localtype;
	      r->Rij1 = r->sites1[r->numpairs - 1];
	      r->Rij2 = r->sites2[r->numpairs - 1];
	    } else {
	      r->Rend = q->Rend;
	      r->Rij1 = Yshift2 + q->Rij1;
	      r->Rij2 = Xshift2 + q->Rij2; // WAS (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->Rij2;
	    }

	    r->score = 0;
	    for(int t = 0; t <= U + U2; t++)
	      r->score += r->iscore[t];
	    r->logPV = alignFP(r, origY, origXrev, origN, origM, p->orientation, p->scaleID, mapid1, mapid2, r->score, 0);
	    
	    delete [] p;
	    delete [] q;
	    p = alignment[i-1] = r;
	    q = alignment[i2] = 0;
	  } else if(refstart > refstart2 && refend > refend2 && 
		    (p->orientation ? (qrystart < qrystart2 && qryend < qryend2) : (qrystart > qrystart2 && qryend > qryend2))){ /* q left of p */
	    if(VERB >= (LVERB?1:2)){
	      printf("q overlaps p from left\n");
	      fflush(stdout);
	    }

	    /* check if overlapped alignment region matches */
	    int t1 = 0, t2 = 0;
	    /* scan ahead in Yshift2 + q->sites1[t1] until it matches refstart */
	    while(Yshift2 + q->sites1[t1] != refstart && t1 < q->numpairs)
	      t1++;
	    if(t1 >= q->numpairs)
	      continue;/* no match found */
	    int origt1 = t1;
	    for(;++t1 < q->numpairs && ++t2 < p->numpairs;)
	      if(Yshift2 + q->sites1[t1] != Yshift + p->sites1[t2] || 
		 (p->orientation ? (Xshift2 + M2 - q->sites2[t1] != Xshift + M - p->sites2[t2]) : (Xshift2 + q->sites2[t1] != Xshift + p->sites2[t2])))
		break;
	    if(t1 < q->numpairs && t2 < p->numpairs)
	      continue;/* not a complete match */
	    
	    t1 = origt1;
	    if(DEBUG) assert(t1 > 0);
	    
	    /* merge q (up to sites t1-1) with p */
	    Calign *r = new Calign[colors==1 ? 1 : colors+1];
	    r->numpairs = t1 + p->numpairs;
	    r->expand_arrays(r->numpairs);
	    r->mapid1 = mapid1;
	    r->mapid2 = mapid2;
	    r->orientation = q->orientation;

	    if(VERB >= (LVERB?1:2)){
	      printf("q overlaps p from left and overlapped region matches : combining into one alignment (t1=%d,p->numpairs=%d,q->numpairs=%d,r->numpairs=%d)\n",
		     t1, p->numpairs,q->numpairs,r->numpairs);
	      fflush(stdout);
	    }
            if(DEBUG)assert(i2 <= i);

	    /* copy q (sites 0..t1-1) */
	    for(int t = 0; t < t1; t++){
	      r->sites1[t] = Yshift2 + q->sites1[t];
	      r->sites2[t] = (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->sites2[t];
	      r->iscore[t] = q->iscore[t];
	      r->outscore[t] = q->outscore[t];
	    }
	    
	    /* copy p */
	    for(int t = 0; t < p->numpairs; t++){
	      r->sites1[t1 + t] = Yshift + p->sites1[t];
	      r->sites2[t1 + t] = (p->orientation ? origM - Xshift - M : Xshift) + p->sites2[t];
	      r->iscore[t1 + t] = p->iscore[t];
	      r->outscore[t1 + t] = p->outscore[t];
	    }
	    r->iscore[r->numpairs] = p->iscore[p->numpairs];
	    r->outscore[r->numpairs] = p->outscore[p->numpairs];

	    if(DEBUG/* HERE >=2 */){
	      for(int t = 1; t < r->numpairs; t++){
		assert(r->sites1[t] > r->sites1[t-1]);
		assert(r->sites2[t] > r->sites2[t-1]);
	      }
	    }

	    if(Xshift2 || Yshift2){
	      r->Lend = localtype;
	      r->Lij1 = r->sites1[0];
	      r->Lij2 = r->sites2[0];
	    } else {
	      r->Lend = q->Lend;
	      r->Lij1 = Yshift2 + q->Lij1;
	      r->Lij2 = Xshift2 + q->Lij2; // WAS (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->Lij2;
	    }
	    if(Xshift + M < origM || Yshift + N < origN){
	      r->Rend = localtype;
	      r->Rij1 = r->sites1[r->numpairs - 1];
	      r->Rij2 = r->sites2[r->numpairs - 1];
	    } else {
	      r->Rend = p->Rend;
	      r->Rij1 = Yshift + p->Rij1;
	      r->Rij2 = Xshift + p->Rij2; // WAS (p->orientation ? origM - Xshift - M : Xshift) + p->Rij2;
	    }

	    r->score = 0.0;
	    for(int t = 0; t <= U + U2; t++)
	      r->score += r->iscore[t];
	    r->logPV = alignFP(r, origY, origXrev, origN, origM, p->orientation, p->scaleID, mapid1, mapid2, r->score, 0);
	    
	    delete [] p;
	    delete [] q;
	    p = alignment[i-1] = r;
	    q = alignment[i2] = 0;
	  }
	}
       }
      }
      delete [] origXrevMem;
    }
  }

  if(VERB){/* display statistics */
    printf("maps=%d:TP pairs=%lld/%lld, total pairs=%lld, positive scores=%lld:\n",
	   nummaps,paircnt,mapcnt,tcnt,poscnt);
    printf("  overlaps: len=%0.3f, sites=%lld, misscnt=%lld (expected %0.1f)\n", 
	   len, sitecnt, misscnt, 2.0*len*(FP[0]*0.01 + FN[0]/Ytheta));
    if(segcnt > 0)
      printf("  segcnt=%lld,norm=%0.4f\n",
	     segcnt,errnorm/segcnt);
    if(Poutlier>0.0)
      printf("  outliers=%lld/%lld (rate=%0.6f)\n",outliercnt,outliercnt+segcnt,outliersum/(outliertot + 0.001));
    fflush(stdout);


    if(RefRepeats > 0){ // WARNING : does not support QUADRATIC_VARIANCE correctly
      double resKB = mres*PixelLen*2.0;
      if(numerrors>0){
	/* generate updated estimate of sizing error parameters Var(x+y) = A + B(x+y) + R(x^2 + y^2) */
	double A = SF[0]*SF[0]*2.0;
	double B = fabs(SD[0])*SD[0];
	double R = SR[0]*SR[0];
	for(int iter=0; iter < 30; iter++){
	  double Lsum=0.0;/* sum of log(A+B(x+y)+R(x^2+y^2)) */
	  double errsum=0.0;/* sum of (x-y)^2/(A+B(x+y)+R(x^2+y^2)) */
	  Cinterval *perr = errors;
	  for(size_t i = 0;i < numerrors;i++,perr++){      
	    double err = perr->y - perr->x;
	    double sum = perr->y + perr->x;
	    err *= err;
	    double var = A + B*sum;
	    if(QUADRATIC_VARIANCE)
	      var += R*(perr->x * perr->x + perr->y * perr->y);
	    Lsum += log(var);
	    errsum += err/var;
	  }

	  double delta = errsum/numerrors;

	  if(VERB>=2){
	    printf("iter=%d:cnt=%lu,sf=%0.6f,sd=%0.6f:LL=%0.6f,norm=%0.6f:sf->%0.6f,sd->%0.6f,sr->%0.6f\n",
		   iter,numerrors,sqrt(A*0.5),sqrt(B),-(Lsum+errsum)*0.5,delta,sqrt(A*delta*0.5),sqrt(B*delta),sqrt(R*delta));
	    fflush(stdout);
	  }

	  B *= delta;
	  A *= delta;
	  R *= delta;

	  /* update estimate of A */
	  double bestLL = MINSCORE;
	  double bestA = A;

	  for(int Aiter=0; Aiter < 3; Aiter++){
	    double Lsum= 0.0;/* sum of log(A+B(x+y)+R(x^2+y^2)) */
	    double Isum= 0.0;/* sum of 1/(A+B(x+y)+R(x^2+y^2)) */
	    double Isum2=0.0;/* sum of 1/(A+B(x+y)+R(x^2+y^2))^2 */
	    double errsum= 0.0; /* sum of (x-y)^2/(A+B(x+y)+R(x^2+y^2)) */
	    double errsum2= 0.0;/* sum of (x-y)^2/(A+B(x+y)+R(x^2+y^2))^2 */
	    double errsum3= 0.0;/* sum of (x-y)^2/(A+B(x+y)+R(x^2+y^2))^3 */
	    Cinterval *perr = errors;
	    for(size_t i = 0;i < numerrors;i++,perr++){
	      double err = perr->y - perr->x;
	      double sum = perr->y + perr->x;
	      err *= err;
	      double var = A + B*sum;
	      if(QUADRATIC_VARIANCE)
		var += R*(perr->x * perr->x + perr->y * perr->y);
	      Lsum += log(var);
	  
	      double Ivar = var = 1.0/var;
	  
	      errsum += err * Ivar;
	      Isum += Ivar;
	  
	      Ivar *= var;

	      errsum2 += err*Ivar;
	      Isum2 += Ivar;
	  
	      errsum3 += err*Ivar*var;
	    }
	    
	    double LL = -(Lsum+errsum)*0.5;
	    double LL1 = 0.5*(errsum2-Isum);
	    double LL2 = 0.5*Isum2-errsum3;

	    if(LL < bestLL){
	      if(VERB>=2){
		printf("iter=%d,%d:sf=%0.8f,sd=%0.8f,sr=%0.8f:LL=%0.8f,LL1=%0.8f,LL2=%0.8f:bestLL=%0.8f,backtracking (sf->%0.8f)\n",
		       Aiter,iter,sqrt(A),copysign(sqrt(fabs(B)),B),sqrt(R),LL, LL1, LL2, bestLL,sqrt(0.5*(A+bestA)));
		fflush(stdout);
	      }
	      A = (A+bestA)*0.5;
	      if(LL < bestLL - 1e-6 && fabs(A-bestA) > 1e-6)
		Aiter--;/* add an iteration if drop is LL was significant */
	      continue;
	    }
	    bestA = A;
	    bestLL = LL;

	    double maxdelta = A * 0.5 + 0.01;
	    double delta = (LL2 < 0.0) ? min(maxdelta,fabs(LL1/LL2)) : maxdelta;
	    if(A + delta < 0.0)
	      delta = -A;
	
	    if(VERB>=2){
	      printf("iter=%d,%d:sf=%0.6f,sd=%0.4f:LL=%0.6f,LL1=%0.3f,LL2=%0.3f:delta=%0.6f(sf->%0.6f)\n",
		     Aiter,iter,sqrt(A*0.5),sqrt(B),-(Lsum+errsum)*0.5,0.5*(errsum2-Isum),0.5*Isum2-errsum3,
		     delta,sqrt((A+delta)*0.5));
	      fflush(stdout);
	    }
	    A += delta;
	    if(QUADRATIC_VARIANCE){
	      if(DEBUG) assert(R >= 0.0);
	      if(B < 0.0 && R >= 0.0)
		A = max(A, B*B/(4.0*R));
	    }
	  }// Aiter = 0..2 */
	  A = bestA;
	}// iter = 0..29
	if(VERB) {
	  printf("  var(y) = sf^2 + sd^2 y : sf = %0.6f -> %0.6f (kb), sd = %0.6f -> %0.6f (kb^1/2), sr=%0.6f -> %0.6f\n", SF[0], sqrt(A*0.5), SD[0], sqrt(B), SR[0], sqrt(R));
	  fflush(stdout);
	}
      }

      if(nummisaligns > 0){/* plot misaligned cuts vs aligned interval size AND estimate scaling for FP & FN */

	/* sort aligned intervals by size */
	qsort(misaligns,nummisaligns,sizeof(Cmisalign), (intcmp*)xinc);
	if(VERB>=2)
	  printf("IntervalSize     MisalignedSites    Density    Density(above %0.2fkb) Samples\n", resKB);
	double msumtot = 0.0,xpsumtot=0.0;
	for(int bin=0;bin < MISALIGN_BINS;bin++){
	  int jmin = (bin*nummisaligns)/MISALIGN_BINS;
	  int jmax = (bin == MISALIGN_BINS-1) ? nummisaligns : ((bin+1)*nummisaligns)/MISALIGN_BINS;
	  if(jmax > jmin){
	    double xsum=0.0,xpsum=0.0;
	    int msum = 0;
	    Cmisalign *pmis = &misaligns[jmin];
	    for(int j = jmin; j < jmax; j++,pmis++){
	      xsum += pmis->x;
	      msum += pmis->m;
	      if(pmis->x > resKB)
		xpsum += pmis->x - resKB;
	    }
	    if(DEBUG) assert(xsum>0.0);
	    if(VERB>=2 && xsum > 0.0 && xpsum > 0.0)
	      printf("%0.3f\t%0.4f\t%0.4f\t%0.4f\t%d\n",
		     xsum/(jmax-jmin), ((double)msum)/(jmax-jmin), msum/xsum, msum/xpsum, jmax-jmin);
	    msumtot += msum;
	    xpsumtot += xpsum;
	  }
	}

	if(xpsumtot>0.0){
	  double delta = mres*PixelLen;
	  double pFN = FN[0];
	  double pTP = 1.0 - pFN;
	  double F = FP[0]*0.01;
	  double lambda = Ylambda;
	  double theta = Ytheta;
	  double Itheta = 1.0/(theta - delta - var*0.5);
	  double R = F+pTP*pFN/lambda;
	  double newR = msumtot/xpsumtot;
	  double newF = F*newR/(F+(Itheta-F)*pFN);
	  double newpTP = (Itheta - newR)/(Itheta - newF);

	  printf("  R= %0.6f -> %0.6f, FP= %0.6f -> %0.6f, FN= %0.6f -> %0.6f \n", R, newR, F*100.0, newF*100.0, FN[0], 1.0 - newpTP);
	  fflush(stdout);

	  if(DEBUG) assert(fabs(newF - F*newR/R) <= newF*0.001);
	}
      }
    }
  }

  double Xscale = 1.0;
  double Yscale = 1.0;
  if(fabs(PixelLen - origPixelLen) > 1e-12){
    Yscale = Xscale = origPixelLen/PixelLen;
    if(PairSplit)/* -i inputs that correspond to reference is not scaled by bpp/500 */
      Yscale = 1.0;
  }

  if(PAIRMERGE_REPEATMASK && PairMerge && RepeatRec && RepeatMaxShift > 0 && RepeatLogPvalueRatio > 0.0){
      /* suppress repeat masking if no extension-merging is possible:
         1. One map completely overlaps the other, OR
	 2. SplitSegDup > 0.0 && overlap > SplitSegDup && largest endoutlier > SplitSegDupE
      */

    if(VERB/* HERE >=2 */){
      printf("Checking to suppress repeats when one map completely overlaps the other\n");
      fflush(stdout);
    }
    
    for(size_t i=0; i < numaligns; i++){
      Calign *p = alignment[i];
      if(!p || p->numpairs <= 0)
	continue;
      int Yid = p->mapid1;
      int Xid = p->mapid2;

      if(p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold || p->repeat <= 0){
	if(EVERB){
	  printf("Alignment %lu/%lu: Yid=%d,id=%lld, Xid=%d,id=%lld, or=%d,sc=%d:numpairs=%d,score=%0.6f,logPV=%0.2f,repeat=%d: Skipping due to -S %0.3f -T %0.2f -A %d (or repeat = 0)\n",
		 i,numaligns,Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,p->orientation,p->scaleID,p->numpairs,p->score,p->logPV,p->repeat,ScoreThreshold,LogPvThreshold,AlignedSiteThreshold);
	  fflush(stdout);
	}
	continue;
      }
      
      int K = p->numpairs;
      Cmap *Ymap = YYmap[Yid];
      Cmap *Xmap = XXmap[Xid];
      if(DEBUG) assert(!Ymap->origmap);
      if(DEBUG) assert(!Xmap->origmap);
      int N = Ymap->numsite[0];
      int M = Xmap->numsite[0];
      FLOAT *Y = Ymap->site[0];
      FLOAT *X = Xmap->site[0];
      
      int LI = p->sites1[0];
      int LJ = p->sites2[0];
      int RI = p->sites1[K-1];
      int RJ = p->sites2[K-1];
      if(DEBUG/* HERE >=2 */){
	assert(0 < LI && LI <= RI && RI <= N);
	assert(0 < LJ && LJ <= RJ && RJ <= M);
      }

      double overlap = (Y[RI]-Y[LI] + (p->orientation ? X[M+1-LJ] - X[M+1-RJ] : X[RJ]-X[LJ]))*0.5;
      double leftEnd = (p->Lend >= -1) ? 0.0 : min(Y[LI], (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]));
      //      int leftEndN = (p->Lend >= -1) ? 0 : (Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ])) ? LI-1 : LJ-1;
      double rightEnd = (p->Rend >= -1) ? 0.0 : min(Y[N+1]-Y[RI], (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]));
      //      int rightEndN = (p->Rend >= -1) ? 0 : (Y[N+1]-Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ])) ? N-RI : M-RJ;

      double OverlappedMargin = 0.0;
      int LeftSmallerY = (Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]) + OverlappedMargin ) ? 1 : 0;
      int LeftSmallerX = (Y[LI] > (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]) - OverlappedMargin ) ? 1 : 0;
      int RightSmallerY = (Y[N+1] - Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]) + OverlappedMargin) ? 1 : 0;
      int RightSmallerX = (Y[N+1] - Y[RI] > (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]) - OverlappedMargin) ? 1 : 0;
      int Overlapped = ((LeftSmallerY && RightSmallerY) || (LeftSmallerX && RightSmallerX)) ? 1 : 0;

      if(Overlapped || (SplitSegDup > 0.0 && !Cfirst && PairMergeHmap <= 0 && overlap > SplitSegDup && max(leftEnd, rightEnd) > SplitSegDupE)){
	if(EVERB){
          if(Overlapped)
	     printf("Alignment %lu/%lu: Yid=%d,id=%lld, Xid=%d,id=%lld, or=%d:numpairs=%d,score=%0.6f,logPV=%0.2f,repeat=%d -> 0: Repeat suppressed since one map completely overlaps the other\n",
		    i,numaligns,Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,p->orientation,p->numpairs,p->score,p->logPV,p->repeat);
	  else
	     printf("Alignment %lu/%lu: Yid=%d,id=%lld, Xid=%d,id=%lld, or=%d:numpairs=%d,score=%0.6f,logPV=%0.2f,repeat=%d -> 0: Repeat suppressed since SplitSegDup and no merger is possible\n",
	  	 i,numaligns,Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,p->orientation,p->numpairs,p->score,p->logPV,p->repeat);
	  fflush(stdout);
	}
	p->repeat = 0;
      }
    }
  }

  if(PairMerge && PairMergeMaxOutlier > 0.0 && PairMergeMaxEndExpandKB > 0.0){/* try to optimize alignment by expanding endoutlier up to PairMergeMaxEndExpandKB to eliminate 
										 internal outliers > PairMergeMaxOutlier located near the ends of the original alignment */
    if(VERB/* HERE >=2 */){
      printf("Checking to expand endoutlier up to %0.3f kb to eliminate internal outliers > %0.3f kb\n",PairMergeMaxEndExpandKB,PairMergeMaxOutlier);
      fflush(stdout);
    }

    FLOAT *XrevMem = new FLOAT[maxM + 2];

    for(size_t i=0; i < numaligns; i++){
      Calign *p = alignment[i];
      if(!p || p->numpairs <= 0)
	continue;
      int Yid = p->mapid1;
      int Xid = p->mapid2;
      if(p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold || p->repeat > 0){
	if(EVERB){
	  printf("Alignment %lu/%lu: Yid=%d,id=%lld, Xid=%d,id=%lld, or=%d:numpairs=%d,score=%0.6f,logPV=%0.2f,repeat=%d: Skipping due to -S %0.3f -T %0.2f -A %d\n",
		 i,numaligns,Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,p->orientation,p->numpairs,p->score,p->logPV,p->repeat,ScoreThreshold,LogPvThreshold,AlignedSiteThreshold);
	  fflush(stdout);
	}
	continue;
      }

      int K = p->numpairs;
      Cmap *Ymap = YYmap[Yid];
      Cmap *Xmap = XXmap[Xid];
      int N = Ymap->numsite[0];
      int M = Xmap->numsite[0];
      FLOAT *Y = Ymap->site[0];
      FLOAT *X = Xmap->site[0];

      FLOAT *Xrev = X;
      if(p->orientation){
	Xrev = XrevMem;
	for(int J = 0; J <= M + 1; J++)
	  Xrev[J] = X[M+1] - X[M+1 - J];
      }

      int LI = p->sites1[0];
      int LJ = p->sites2[0];
      int RI = p->sites1[K-1];
      int RJ = p->sites2[K-1];
      if(DEBUG/* HERE >=2 */){
	assert(0 < LI && LI <= RI && RI <= N);
	assert(0 < LJ && LJ <= RJ && RJ <= M);
      }

      if(PairMergeMaxEndExpandKB >= PAIRMERGE_HUGE){/* avoid unlimited endoutlier expansion if original alignment does not satisfy conditions for -SplitSegDup */
	double leftEnd = (p->Lend >= -1) ? 0.0 : min(Y[LI], (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]));
	int leftEndN = (p->Lend >= -1) ? 0 : (Y[LI] < (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ])) ? LI-1 : LJ-1;
	double rightEnd = (p->Rend >= -1) ? 0.0 : min(Y[N+1]-Y[RI], (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]));
	int rightEndN = (p->Rend >= -1) ? 0 : (Y[N+1]-Y[RI] < (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ])) ? N-RI : M-RJ;

        bool LeftEndSegDup = (leftEnd > SplitSegDupE && leftEndN > SplitSegDupEN);
        bool RightEndSegDup = (rightEnd > SplitSegDupE && rightEndN > SplitSegDupEN);
	double LeftLen = -1.0, RightLen = -1.0;
	int LeftN = -1, RightN = -1;
        int LeftLongerY = -1, RightLongerY = -1;/* If Y end is longer than X end on Left or Right end of alignment */
	if(PAIRMERGE_SEGDUP_ONESIDED){
	  LeftLen = max(Y[LI], (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]));
	  LeftN = (Y[LI] > (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ])) ? LI-1 : LJ-1;
	  LeftLongerY = (Y[LI] >= LeftLen) ? 1 : 0;

	  RightLen = max(Y[N+1] - Y[RI], (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]));
	  RightN = (Y[N+1] - Y[RI] >  (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ])) ? N - RI : M - RJ;
	  RightLongerY = (Y[N+1] - Y[RI] >= RightLen) ? 1 : 0;
	}
	
	size_t YL2 = Ymap->Mask[0] ? Ymap->Mask[0][1] & END_NOEXT2 : 0;
	size_t YR2 = Ymap->Mask[0] ? Ymap->Mask[0][N+1] & END_NOEXT2 : 0;
	size_t XL2 = Xmap->Mask[0] ? Xmap->Mask[0][p->orientation ? M+1 : 1] & END_NOEXT2 : 0;
	size_t XR2 = Xmap->Mask[0] ? Xmap->Mask[0][p->orientation ? 1 : M+1] & END_NOEXT2 : 0;

	bool SegDupEnds = !PAIRMERGE_SEGDUP_ONESIDED ? (LeftEndSegDup && RightEndSegDup) :
	  ((LeftEndSegDup || RightEndSegDup) && min(LeftLen,RightLen) > SplitSegDupE && min(LeftN,RightN) > SplitSegDupEN
	   && (PAIRMERGE_SEGDUP_ONESIDED <= 1 || (!LeftEndSegDup ? (rightEnd >= SplitSegDup_OneSided_MaxLen || (LeftLongerY ? !XL2 : !YL2)) : 
						                   (leftEnd >= SplitSegDup_OneSided_MaxLen || (RightLongerY ? !XR2 : !YR2)))));
	if(VERB>=3){
	  printf("Ymap->id=%lld,Xmap->id=%lld,or=%d:numpairs=%d: SegDupEnds=%d\n", Ymap->id,Xmap->id,p->orientation,K,SegDupEnds);
	  printf("\t LeftEndSegDup=%d,RightEndSegDup=%d,LeftLen=%0.3f,%d,RightLen=%0.3f,%d(SplitSegDupE=%0.3f,%d),leftEnd=%0.3f,rightEnd=%0.3f,SplitSegDup_OneSided=%d,%0.2f,LeftLongerY=%d,RightLongerY=%d\n",
		 LeftEndSegDup,RightEndSegDup,LeftLen,LeftN,RightLen,RightN,SplitSegDupE,SplitSegDupEN,leftEnd,rightEnd,PAIRMERGE_SEGDUP_ONESIDED,SplitSegDup_OneSided_MaxLen,LeftLongerY,RightLongerY);
	  fflush(stdout);
	}

	if(!SegDupEnds)
	  continue;
      }

      /* locate outliers near left end and remove from alignment, provided leftEnd does not increase beyond PairMergeMaxEndExpandKB */
      int lastI = LI, I;
      int lastJ = LJ, J;
      int UL = 0;
      for(int k = 1; k < K; lastI = I, lastJ = J, k++){
	I = p->sites1[k];
	J = p->sites2[k];
	if(DEBUG>=2) assert(J >= lastJ);
	if(DEBUG>=2) assert(I >= lastI);
	  
	if(min(Y[I], (p->orientation ? X[M+1] - X[M+1-J] : X[J])) > PairMergeMaxEndExpandKB &&
	   (Y[I] < (p->orientation ? X[M+1] - X[M+1-J] : X[J]) ? I-1 : J-1) >= PairMergeMaxEndN)/* Cannot "move" start of alignment to (I,J) */
	  break;

	/* outlier lies between sites lastI .. I on reference and lastJ .. J on query */
	double qrystartpos = (p->orientation ? (X[M + 1] - X[M + 1 - lastJ]) : X[lastJ]) * Xscale; /* QryStartPos */
	double qryendpos   = (p->orientation ? (X[M + 1] - X[M + 1 - J]) : X[J])*Xscale; /* QryEndPos */
	double refstartpos = Y[lastI]*Yscale; /* RefStartPos */
	double refendpos   = Y[I]*Yscale; /* RefEndPos */
	if(DEBUG>=2) assert(refendpos >= refstartpos);
	if(DEBUG) assert(qryendpos >= qrystartpos);// If this fails, fix next line
	double maxOutlier = fabs(refendpos - refstartpos - qryendpos + qrystartpos);
	int maxOutlierM = I-lastI + J-lastJ - 2;// misaligned labels

	int outlier = (p->iscore[k] > p->outscore[k] + (RFLOAT)0.01 || (outlierExtend && (lastJ >= J || lastI >= I)) || maxOutlier > PairMergeMaxOutlier || maxOutlierM > PairMergeMaxOutlierM) ? 1 : 0;
	if(!outlier)/* not an outlier */
	  continue;

	if(PairMergeOutlierChiSq > 0.0 && CmapChimQuality >= 3){// skip if outlier variance is too large or has too many molecule outliers
	  double Outlier = maxOutlier;
	  if(refendpos - refstartpos < qryendpos - qrystartpos){
	    float *OutlierFrac = Ymap->OutlierFrac[0];
	    float *FragSd = Ymap->FragSd[0];
	    float *ExpSd = Ymap->ExpSd[0];
	    float *FragCov = Ymap->FragCov[0];
	    float *FragChiSq = Ymap->FragChiSq[0];
	    if(DEBUG>=2) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);

	    float *ChimQuality = Ymap->ChimQuality[0];
	    float *ChimNorm = Ymap->ChimNorm[0];

	    double LRange = (p->orientation ? X[M+1-lastJ] - X[M+1-J] : X[J] - X[lastJ]) * Xscale;
	    double LRangeRatio = LRange / AvgInterval;

	    int i = lastI;
	    for(; i < I; i++){
	      double SRange = (Y[i+1] - Y[i]) * Yscale;
	      double FragCovi = FragCov[i];
	      if(PairMergeChimNorm > 0.0 && ChimNorm[i] > 0.0 && ChimNorm[i+1] > 0.0)
		FragCovi = min(FragCovi, min(ChimNorm[i] * ChimQuality[i], ChimNorm[i+1] * ChimQuality[i+1]) * PairMergeChimNorm * 0.01);

	      if(Outlier < PairMergeOutlierMaxSize &&
		 (OutlierFrac[i] > PairMergeOutlierFrac || 
		  (FragCovi < PairMergeOutlierMaxCov && 
		   LRange > PairMergeOutlierMinInterval &&
		   ((FragChiSq[i] < PairMergeOutlierChiSq &&
		     FragSd[i] > ExpSd[i] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
		     FragSd[i] > ExpSd[i] * PairMergeOutlierSdRatio) ||		     
		    FragCovi < PairMergeOutlierMinCov + 
		    ((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != abs(J-lastJ))?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale))))
		break;
	    }
	    if(i < I)
	      continue;// outlier will be ignored
	  } else {
	    float *OutlierFrac = Xmap->OutlierFrac[0];
	    float *FragSd = Xmap->FragSd[0];
	    float *ExpSd = Xmap->ExpSd[0];
	    float *FragCov = Xmap->FragCov[0];
	    float *FragChiSq = Xmap->FragChiSq[0];
	    if(DEBUG>=2) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);

	    float *ChimQuality = Xmap->ChimQuality[0];
	    float *ChimNorm = Xmap->ChimNorm[0];

	    double LRange = (Y[I] - Y[lastI]) * Yscale;
	    double LRangeRatio = LRange / AvgInterval;

	    int jmin = p->orientation ? M+1-J : lastJ;
	    int jmax = p->orientation ? M+1-lastJ : J;
	    int j = jmin;
	    for(; j < jmax; j++){
	      double SRange = (X[j+1] - X[j]) * Xscale;
	      double FragCovj = FragCov[j];
	      if(PairMergeChimNorm > 0.0 && ChimQuality[j] > 0.0 && ChimNorm[j+1] > 0.0)
		FragCovj = min(FragCovj, min(ChimNorm[j] * ChimQuality[j], ChimNorm[j+1] * ChimQuality[j+1]) * PairMergeChimNorm * 0.01);

	      if(Outlier < PairMergeOutlierMaxSize &&
		 (OutlierFrac[j] > PairMergeOutlierFrac || 
		  (FragCovj < PairMergeOutlierMaxCov && 
		   LRange > PairMergeOutlierMinInterval &&
		   ((FragChiSq[j] < PairMergeOutlierChiSq &&
		     FragSd[j] > ExpSd[j] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
		     FragSd[j] > ExpSd[j] * PairMergeOutlierSdRatio) ||		     
		    FragCovj < PairMergeOutlierMinCov +
		    ((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != abs(J-lastJ))?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale))))
		break;
	    }
	    if(j < jmax)
	      continue;// outlier will be ignored
	  }
	}

	double YLeft = (Y[lastI+1]-Y[lastI]) * Yscale;
	double YRight = (Y[I] - Y[I-1]) * Yscale;
	double XLeft = (p->orientation ? X[M+1-lastJ] - X[M+1-(lastJ+1)] : X[lastJ+1] - X[lastJ]) * Xscale;
	double XRight = (p->orientation ? X[M+1-(J-1)] - X[M+1-J] : X[J] - X[J-1]) * Xscale;

	if(I - lastI + J - lastJ - 2 >= PairMergeMaxOutlierBalanced){
	  if(I > lastI + 2){
	    double Middle = (Y[I-1] - Y[lastI+1]) * Yscale;
	    double MinOutlier = Middle + max(0.0, max(YLeft-XLeft,YRight-XRight));
	    maxOutlier = max(maxOutlier, MinOutlier);
	  }
	  if(J > lastJ + 2){
	    double Middle = (p->orientation ? X[M+1-(lastJ+1)] - X[M+1-(J-1)] : X[J-1] - X[lastJ+1]) * Xscale;
	    double MinOutlier = Middle + max(0.0, max(XLeft-YLeft,XRight-YRight));
	    maxOutlier = max(maxOutlier, MinOutlier);
	  }
	}

	if(maxOutlier > PairMergeMaxOutlier || maxOutlierM > PairMergeMaxOutlier){/* move start of alignment to (I,J) */
          if(EVERB){
            printf("\t k=%d:left end outlier at Y[%d..%d]= %0.3f..%0.3f,Y[N+1]=%0.3f,X[%d..%d]= %0.3f..%0.3f,X[M+1]=%0.3f:iscore=%0.6f,outscore=%0.6f:Outlier=%0.3f kb:UL= %d -> %d\n",
		   k, lastI,I,Y[lastI],Y[I],Y[N+1],lastJ,J,(p->orientation ? X[M+1]-X[M+1-lastJ] : X[lastJ]),(p->orientation ? X[M+1]-X[M+1-J] : X[J]),X[M+1],
		   p->iscore[k],p->outscore[k],maxOutlier,UL,k);
	    fflush(stdout);
	  }
	  UL = k;
	}
      }

      int UR = K-1;
      I = RI;
      J = RJ;
      for(int k = K; --k > 0; I = lastI, J = lastJ){
	lastI = p->sites1[k-1];
	lastJ = p->sites2[k-1];
	if(DEBUG>=2) assert(J >= lastJ);
	if(DEBUG>=2) assert(I >= lastI);

	if(min(Y[N+1]-Y[lastI], (p->orientation ? X[M+1-lastJ] : X[M+1] - X[lastJ])) > PairMergeMaxEndExpandKB && 
	   (Y[N+1]-Y[lastI] < (p->orientation ? X[M+1-lastJ] : X[M+1] - X[lastJ]) ? N-lastI : M-lastJ) >= PairMergeMaxEndN)/* Cannot "move" end of alignment to lastI,lastJ */
	  break;

	/* outlier lies between sites lastI .. I on reference and lastJ .. J on query */
	double qrystartpos = (p->orientation ? (X[M + 1] - X[M + 1 - lastJ]) : X[lastJ]) * Xscale; /* QryStartPos */
	double qryendpos   = (p->orientation ? (X[M + 1] - X[M + 1 - J]) : X[J])*Xscale; /* QryEndPos */
	double refstartpos = Y[lastI]*Yscale; /* RefStartPos */
	double refendpos   = Y[I]*Yscale; /* RefEndPos */
	if(DEBUG>=2) assert(refendpos >= refstartpos);
	if(DEBUG) assert(qryendpos >= qrystartpos);// If this fails, fix next line
	      
	double maxOutlier = fabs(refendpos - refstartpos - qryendpos + qrystartpos);// WAS440 fabs(refendpos - refstartpos -qryendpos - qrystartpos)
	int maxOutlierM = I-lastI + J-lastJ - 2;// misaligned labels

	int outlier = (p->iscore[k] > p->outscore[k] + (RFLOAT)0.01 || (outlierExtend && (lastJ >= J || lastI >= I)) || maxOutlier > PairMergeMaxOutlier || maxOutlierM > PairMergeMaxOutlierM) ? 1 : 0;
	if(!outlier)/* not an outlier */
	  continue;
	      
	if(PairMergeOutlierChiSq > 0.0 && CmapChimQuality >= 3){// skip if outlier variance is too large or has too many molecule outliers
	  double Outlier = maxOutlier;
	  if(refendpos - refstartpos < qryendpos - qrystartpos){
	    float *OutlierFrac = Ymap->OutlierFrac[0];
	    float *FragSd = Ymap->FragSd[0];
	    float *ExpSd = Ymap->ExpSd[0];
	    float *FragCov = Ymap->FragCov[0];
	    float *FragChiSq = Ymap->FragChiSq[0];
	    if(DEBUG>=2) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);
	    
	    float *ChimQuality = Ymap->ChimQuality[0];
	    float *ChimNorm = Ymap->ChimNorm[0];

	    double LRange = (p->orientation ? X[M+1-lastJ] - X[M+1-J] : X[J] - X[lastJ]) * Xscale;
	    double LRangeRatio = LRange / AvgInterval;

	    int i = lastI;
	    for(; i < I; i++){
	      double SRange = (Y[i+1] - Y[i]) * Yscale;
	      double FragCovi = FragCov[i];
	      if(PairMergeChimNorm > 0.0 && ChimNorm[i] > 0.0 && ChimNorm[i+1] > 0.0)
		FragCovi = min(FragCovi, min(ChimNorm[i] * ChimQuality[i], ChimNorm[i+1] * ChimQuality[i+1]) * PairMergeChimNorm * 0.01);

	      if(Outlier < PairMergeOutlierMaxSize &&
		 (OutlierFrac[i] > PairMergeOutlierFrac || 
		  (FragCovi < PairMergeOutlierMaxCov && 
		   LRange > PairMergeOutlierMinInterval &&
		   ((FragChiSq[i] < PairMergeOutlierChiSq &&
		     FragSd[i] > ExpSd[i] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
		     FragSd[i] > ExpSd[i] * PairMergeOutlierSdRatio) ||		     
		    FragCovi < PairMergeOutlierMinCov + 
		    ((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != abs(J-lastJ))?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale))))
		break;
	    }
	    if(i < I)
	      continue;// outlier will be ignored
	  } else {
	    float *OutlierFrac = Xmap->OutlierFrac[0];
	    float *FragSd = Xmap->FragSd[0];
	    float *ExpSd = Xmap->ExpSd[0];
	    float *FragCov = Xmap->FragCov[0];
	    float *FragChiSq = Xmap->FragChiSq[0];
	    if(DEBUG>=2) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);

	    float *ChimQuality = Xmap->ChimQuality[0];
	    float *ChimNorm = Xmap->ChimNorm[0];

	    double LRange = (Y[I] - Y[lastI]) * Yscale;
	    double LRangeRatio = LRange / AvgInterval;

	    int jmin = p->orientation ? M+1-J : lastJ;
	    int jmax = p->orientation ? M+1-lastJ : J;
	    int j = jmin;
	    for(; j < jmax; j++){
	      double SRange = (X[j+1] - X[j]) * Xscale;
	      double FragCovj = FragCov[j];
	      if(PairMergeChimNorm > 0.0 && ChimNorm[j] > 0.0 && ChimNorm[j+1] > 0.0)
		FragCovj = min(FragCovj, min(ChimNorm[j]*ChimQuality[j],ChimNorm[j+1]*ChimQuality[j+1]) * PairMergeChimNorm * 0.01);

	      if(Outlier < PairMergeOutlierMaxSize &&
		 (OutlierFrac[j] > PairMergeOutlierFrac || 
		  (FragCovj < PairMergeOutlierMaxCov && 
		   LRange > PairMergeOutlierMinInterval &&
		   ((FragChiSq[j] < PairMergeOutlierChiSq &&
		     FragSd[j] > ExpSd[j] + PairMergeOutlierMinSf + PairMergeOutlierMinSr * SRange &&
		     FragSd[j] > ExpSd[j] * PairMergeOutlierSdRatio) ||		     
		    FragCovj < PairMergeOutlierMinCov + 
		    ((PairMergeOutlierMinCovRatio > 0.0 && I-lastI != abs(J-lastJ))?min(PairMergeOutlierMinCovRatio,LRangeRatio):LRangeRatio) *  PairMergeOutlierMinCovScale))))
		break;
	    }
	    if(j < jmax)
	      continue;// outlier will be ignored
	  }
	}

	double YLeft = (Y[lastI+1]-Y[lastI]) * Yscale;
	double YRight = (Y[I] - Y[I-1]) * Yscale;
	double XLeft = (p->orientation ? X[M+1-lastJ] - X[M+1-(lastJ+1)] : X[lastJ+1] - X[lastJ]) * Xscale;
	double XRight = (p->orientation ? X[M+1-(J-1)] - X[M+1-J] : X[J] - X[J-1]) * Xscale;

	if(I - lastI + J - lastJ - 2 >= PairMergeMaxOutlierBalanced){
	  if(I > lastI + 2){
	    double Middle = (Y[I-1] - Y[lastI+1]) * Yscale;
	    double MinOutlier = Middle + max(0.0, max(YLeft-XLeft,YRight-XRight));
	    maxOutlier = max(maxOutlier, MinOutlier);
	  }
	  if(J > lastJ + 2){
	    double Middle = (p->orientation ? X[M+1-(lastJ+1)] - X[M+1-(J-1)] : X[J-1] - X[lastJ+1]) * Xscale;
	    double MinOutlier = Middle + max(0.0, max(XLeft-YLeft,XRight-YRight));
	    maxOutlier = max(maxOutlier, MinOutlier);
	  }
	}
	if(maxOutlier > PairMergeMaxOutlier || maxOutlierM > PairMergeMaxOutlierM){/* move end of alignment to (lastI,lastJ) */
          if(EVERB){
            printf("\t k=%d:right end outlier at Y[%d..%d]= %0.3f..%0.3f,Y[N+1]=%0.3f,X[%d..%d]= %0.3f..%0.3f,X[M+1]=%0.3f:iscore=%0.6f,outscore=%0.6f:Outlier=%0.3f kb:UR= %d -> %d\n",
		   k, lastI,I,Y[lastI],Y[I],Y[N+1],lastJ,J,(p->orientation ? X[M+1]-X[M+1-lastJ] : X[lastJ]),(p->orientation ? X[M+1]-X[M+1-J] : X[J]),X[M+1],
		   p->iscore[k],p->outscore[k],maxOutlier,UR,k-1);
	    fflush(stdout);
	  }
	  UR = k-1;
        }
      }

      if(PairMergeMaxEndExpandKB >= PAIRMERGE_HUGE && UR < UL){/* use different logic since UR < UL is possible */
	if(VERB/* HERE >=2 */ && EVERB){
	  if(max(UR + 1, K - UL) < AlignedSiteThreshold)
	    printf("Ymap->id=%lld,Xmap->id=%lld,or=%d:numpairs=%d,UL=%d(I=%d,J=%d),UR=%d(I=%d,J=%d): cannot fix alignment\n",
		   Ymap->id,Xmap->id,p->orientation,K,UL,p->sites1[UL],p->sites2[UL],UR,p->sites1[UR],p->sites2[UR]);
	  else if(UR + 1 > K-UL)
	    printf("Ymap->id=%lld,Xmap->id=%lld,or=%d:numpairs=%d,UL=%d(I=%d,J=%d),UR=%d(I=%d,J=%d): Truncating alignment to 0..UR\n",
		   Ymap->id,Xmap->id,p->orientation,K,UL,p->sites1[UL],p->sites2[UL],UR,p->sites1[UR],p->sites2[UR]);
	  else
	    printf("Ymap->id=%lld,Xmap->id=%lld,or=%d:numpairs=%d,UL=%d(I=%d,J=%d),UR=%d(I=%d,J=%d): Truncating alignment to UL .. U-1\n",
		   Ymap->id,Xmap->id,p->orientation,K,UL,p->sites1[UL],p->sites2[UL],UR,p->sites1[UR],p->sites2[UR]);
	  fflush(stdout);
	}

	if(max(UR + 1, K -UL) < AlignedSiteThreshold)
	  continue;/* cannot fix alignment to remove outliers > PairMergeMaxOutlier */

	if(UR + 1 > K - UL)
	  UL = 0;
	else
	  UR = K-1;

      } else {/* normal finite value of PairMergeMaxEndExpandKB */

	if(VERB/* HERE >=2 */ && EVERB){
	  if(UR-UL <= 0 || UR-UL < AlignedSiteThreshold)
	    printf("Ymap->id=%lld,Xmap->id=%lld,or=%d:numpairs=%d,UL=%d(I=%d,J=%d),UR=%d(I=%d,J=%d) : cannot fix alignment\n",
		   Ymap->id,Xmap->id,p->orientation,K,UL,p->sites1[UL],p->sites2[UL],UR,p->sites1[UR],p->sites2[UR]);
	  else if(!(UL > 0 || UR < K-1)){
	    if(EVERB)
	      printf("Ymap->id=%lld,Xmap->id=%lld,or=%d:numpairs=%d,UL=%d(I=%d,J=%d),UR=%d(I=%d,J=%d) : no change\n",
		     Ymap->id,Xmap->id,p->orientation,K,UL,p->sites1[UL],p->sites2[UL],UR,p->sites1[UR],p->sites2[UR]);
	  } else
	    printf("Ymap->id=%lld,Xmap->id=%lld,or=%d:numpairs=%d,UL=%d(I=%d,J=%d),UR=%d(I=%d,J=%d) : Truncating alignment to UL..UR\n",
		   Ymap->id,Xmap->id,p->orientation,K,UL,p->sites1[UL],p->sites2[UL],UR,p->sites1[UR],p->sites2[UR]);
	  fflush(stdout);
	}

	if(UR - UL <= 0 || UR - UL < AlignedSiteThreshold)
	  continue;/* cannot fix alignment to remove outliers > PairMergeMaxOutlier */
      }
      if(UL > 0 || UR < K-1){/* truncate alignment and recompute logPV  */
	double origscore = p->score;
	double origLogPV = p->logPV;

	//	p->logPV = alignFP(p, Y, Xrev, N, M, p->orientation, p->scaleID,Yid, Xid, p->score, 1);
	//	printf("\t original logPV= %0.6f (%0.6f)\n",origLogPV, p->logPV);
	//	fflush(stdout);

	int k = 0;
	if(UL > 0){
	  p->Lend = -2;
	  p->Lij1 = p->sites1[UL];
	  p->Lij2 = p->sites2[UL];
	  p->iscore[0] = p->outscore[0] = ChimScore;
	}

	for(int U = UL; U <= UR; U++, k++){
	  p->sites1[k] = p->sites1[U];
	  p->sites2[k] = p->sites2[U];
	  p->iscore[k+1] = p->iscore[U+1];
	  p->outscore[k+1] = p->outscore[U+1];
	}
	if(UR < K-1){
	  p->Rend = -2;
	  p->Rij1 = p->sites1[k-1];
	  p->Rij2 = p->sites2[k-1];
	  p->iscore[k] = p->outscore[k] = ChimScore;
	}
	p->numpairs = k;

	p->score = 0.0;
	for(int t = 0; t <= k; t++)
	  p->score += p->iscore[t];

	if(p->numpairs <= 1){
	  p->logPV = 0.0;
	  continue;
	}

	p->logPV = alignFP(p, Y, Xrev, N, M, p->orientation, p->scaleID,Yid, Xid, p->score, EVERB ? 1 : 0);

	if(EVERB){
	  printf("Ymap->id=%lld,Xmap->id=%lld,or=%d:numpairs=%d,UL=%d(I=%d,J=%d),UR=%d(I=%d,J=%d) : score= %0.6f -> %0.6f, logPV= %0.6f -> %0.6f\n",
		 Ymap->id,Xmap->id,p->orientation,K,UL,p->sites1[UL],p->sites2[UL],UR,p->sites1[UR],p->sites2[UR], origscore,p->score,origLogPV, p->logPV);
	  fflush(stdout);
	}
      } // if(UL > 0 || UR < K-1)
    } // for i = 0 .. numaligns - 1

    delete [] XrevMem;
  } // if (PairMerge && PairMergeMaxOutlier > 0.0 && PairMergeMaxEndExpandKB > 0.0)

  int orignummaps = nummaps;
  
  if(PairMerge){/* reset paired flag for all input maps */
    for(int i=0; i < numX; i++)
      XXmap[i]->paired = 0;
    for(int i=0; i < numY; i++)
      YYmap[i]->paired = 0;
  }

  char prefix[PATH_MAX];
  if(PairMerge && PairMergeRepeat)
    sprintf(prefix,"%sM%d",output_prefix,pairmergeIter);
  else
    strcpy(prefix,output_prefix);
  if(PairSplit || (/* PairMerge && */ PairMergeXmap) || RepeatShift > 0.0) {
    if(mres > 0.001 || mresSD > 0.0 || fabs(PixelLen - origPixelLen) > 1e-12){
      /* create reduced resolution Reference and Query map files : <prefix>_r.cmap and <prefix>_q.cmap */
      char cmapprefix[PATH_MAX];
      sprintf(cmapprefix,"%s_r",prefix);
      output_cmap(cmapprefix,YYmap,0,orignumY-1);

      sprintf(cmapprefix,"%s_q",prefix);
      output_cmap(cmapprefix,XXmap,0,orignumX-1);
    } else if(RepeatShift > 0.0){
      RefIndex = num_files;
      QueryIndex = 0;
    }
    output_xmap(prefix, vfixx_filename[0], 0, numaligns,false);/* output .xmap file */
    if(PairSplit && dosvdetect)
      output_smap(output_prefix, vfixx_filename[0]);/* output .smap file */
  }

  if(VERB/* >=2 */){
    printf("pairalign:numX=%d,numY=%d,nummaps=%d,numrefmaps=%d,numaligns=%llu\n",
	    numX,numY,nummaps,numrefmaps,(unsigned long long)numaligns);
    fflush(stdout);
  }

  int fmaps = 0;
  int splitcnt = 0;
  int modified = 0;/* number of maps that have been modified (segdup mask) or changed : forces map to be output, instead of linked to previous file */
  int changed = 0;/* number of maps that have been changed : will be a new map in Gmap[orignummaps .. nummaps-1] */

  if(PairMerge){
    if(VERB>=3){
      printf("Before pairmergeIter=%d: nummaps=%d,orignummaps=%d\n",pairmergeIter,nummaps,orignummaps);
      for(int i = 0; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	int M = pmap->numsite[0];
	printf("map[%d]=%p : mapid=%d,id=%lld,paired=%d,len=%0.4f,M=%d\n",
	       i,pmap,pmap->mapid,pmap->id,pmap->paired,pmap->site[0][M+1],M);
      }

      fflush(stdout);
    }

    splitcnt = pairmerge(pairmergeIter);/* generates merged maps at Gmap[nummaps...] and tags merged maps with paired=1 and new/modified maps with centercloned=1) */
    // NOTE : If PairMergeHmap > 0, instead verify that no maps can merge and instead create .hmaps by combining maps that align with large outliers  (or matched endoutliers)

    Cmap **nmap = new Cmap*[nummaps];
    double N50Sum= 0.0, N50Len = 0.0, totlen = 0.0;

    fmaps = modified = changed = 0;
    for(int i = 0; i < nummaps;i++){
      Cmap *pmap = Gmap[i];
      if(DEBUG) assert(!pmap->origmap);
      if(pmap->paired)
	continue;
      
      nmap[fmaps++] = pmap;
      int M = pmap->numsite[0];
      totlen += pmap->site[0][M+1];

      if(DEBUG && !(pmap->id > 0)){
	printf("i=%d/%d:map[i]->id = %lld, mapid= %d, paired=%d, centercloned=%d, orignummaps=%d\n",i,nummaps,pmap->id,pmap->mapid,pmap->paired,pmap->centercloned,orignummaps);
	fflush(stdout);
	assert(pmap->id > 0);
      }

      if(pmap->centercloned){
	if(DEBUG/* HERE >=2 */ && pairmergeIter==1 && !SegDupMask && !(i >= orignummaps)){
	  printf("Gmap[%d]->centercloned=%d:orignummaps=%d,nummaps=%d\n",i,Gmap[i]->centercloned,orignummaps,nummaps);
	  fflush(stdout);
	  assert(i >= orignummaps);
	}
	modified++;
	if(i >= orignummaps)
	  changed++;
      }
    }

    extern int PmapLenDec(Cmap **p1, Cmap **p2);

    qsort(nmap, fmaps, sizeof(Cmap*), (intcmp*)PmapLenDec);/* sort map pointers in descending order of length */
    for(int i = 0; i < fmaps; i++){
      Cmap *pmap = nmap[i];
      int M = pmap->numsite[0];
      if((N50Sum += pmap->site[0][M+1]) >= totlen * 0.5){
	N50Len = pmap->site[0][M+1];
	break;
      }
    }
    if(DEBUG && totlen > 0.0) assert(N50Len > 0.0);
    delete [] nmap;

    if(VERB){
      int merged = orignummaps - fmaps + splitcnt;
      printf("pairmerge: iter=%d: merged %d map pairs, split %d maps (%d maps remain, including %d modified and %d changed maps) totlen= %0.3f Mb (avg= %0.3f kb, N50= %0.4f Mb): cum wall time=%0.6f\n",
	     pairmergeIter, merged, splitcnt, fmaps, modified, changed, totlen * 0.001, totlen / max(fmaps, 1), N50Len * 0.001, wtime());
      //      printf("    nummaps=%d, PairMergeRepeat=%d,PairMergeHmap=%d\n",nummaps,PairMergeRepeat,PairMergeHmap);
      fflush(stdout);
    }
    if(DEBUG) assert(changed == nummaps - orignummaps);
  }

  if(!PairSplit && !PairMerge)
    output_align(prefix,0);/* output .align file */

  int PairMergeProgress = 0;
  if(PairMerge){
    if(VERB>=3){
      printf("After pairmergeIter=%d: nummaps=%d,orignummaps=%d\n",pairmergeIter,nummaps,orignummaps);
      for(int i = 0; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	int M = pmap->numsite[0];
	printf("Gmap[%d]=%p : mapid=%d,id=%lld,paired=%d,len=%0.4f,M=%d\n",
	       i,pmap,pmap->mapid,pmap->id,pmap->paired,pmap->site[0][M+1],M);
      }
      printf("\n");
      fflush(stdout);
    }

    if(DEBUG/* HERE >=2 */){/* make sure there are no duplicate IDs (excluding maps with paired != 0) */
      Cmap **smap = new Cmap*[nummaps];
      int cnt = 0;

      for(int i = 0; i < nummaps; i++){
	Gmap[i]->mapid = i;
	if(Gmap[i]->paired)
	  continue;
	smap[cnt++] = Gmap[i];
      }
      
      qsort(smap,cnt,sizeof(Cmap *),(intcmp*)CmapIdInc);
      
      for(int i = 1; i < cnt; i++){
	if(smap[i]->id == smap[i-1]->id){
	  printf("%d:Duplicate map ids:\n",pairmergeIter);
	  printf("\tmap[%d] : id=%lld, M=%d, Len= %0.3f kb\n",smap[i-1]->mapid,smap[i-1]->id,smap[i-1]->numsite[0],smap[i-1]->site[0][smap[i-1]->numsite[0]+1]);
	  for(int j = i; j < cnt; j++)
	    if(smap[j]->id == smap[i-1]->id)
	      printf("\tmap[%d] : id=%lld, M=%d, Len= %0.3f kb\n",smap[j]->mapid,smap[j]->id,smap[j]->numsite[0],smap[j]->site[0][smap[j]->numsite[0]+1]);	      
	  fflush(stdout);
	  assert(smap[i]->id  > smap[i-1]->id);
	}
      }

      delete [] smap;
    }

    extern double origSplitSegDup;
    extern int PairMergeIDrenumbered;

    PairMergeProgress = max(max(orignummaps - fmaps + splitcnt,splitcnt), nummaps - orignummaps); // WAS167 nummaps - orignummaps;
    if(VERB/* HERE HERE >=2 */){
      printf("iter=%d/%d: PairMergeProgress= %d: PairMergeRepeat=%d,PairMergeHmap=%d,PairMergePreferNGS=%d,SplitSegDup=%0.3f,origSplitSegDup=%0.3f\n",
	     pairmergeIter,pairmergeIterLast,PairMergeProgress,PairMergeRepeat,PairMergeHmap,PairMergePreferNGS,SplitSegDup,origSplitSegDup);
      fflush(stdout);
    }

    if((!(PairMerge && PairMergeRepeat && PairMergeProgress) || PairMergeHmap > 0 || 
	(pairmergeIter >= pairmergeIterLast && (!PairMergePreferNGS || SplitSegDup)))

       && !(PairMerge && PairMergeRepeat && !SplitSegDup && origSplitSegDup)
       && PairMergeHmap >= 0){/* output merged/modified and unchanged maps */

      extern int PairMergeCmapOutput;
      PairMergeCmapOutput = 1;

      int linkcnt = 0;/* number of files linked to previous copy instead of being written */
      double origwt = wtime();

      if(PairMergeIDrenumber){/* renumber all remaining (output) contig ids sequentially */
	PairMergeIDrenumbered = 1;
	long long nextid = 1;
	char IDmap[PATH_MAX];
	sprintf(IDmap,"%s.cidmap",output_prefix);
	FILE *fpID;
	if((fpID = fopen(IDmap,"w")) == NULL){
	  printf("Cannot write file %s\n",IDmap);
	  fflush(stdout);exit(1);
	}
	printversion(fpID);
	fprintf(fpID,"# newID oldID (NOTE: there may be multiple lines with the same newID: multiple input maps can be merged into one output map)\n");
	if(SplitSegDup)
	  fprintf(fpID,"#            (NOTE: there may be multiple lines with the same oldID if input maps were split at suspected segdup regions)\n");

	long long *newid = new long long[nummaps];
	CID *CIDmap = new CID[MergedNummaps + nummaps];
	int CIDcnt = 0;

	/* sort remaining maps Gmap[0..nummaps-1] in descending order of size */
	Cmap** nmap = new Cmap*[nummaps];
	for(int i = 0; i < nummaps; i++)
	  nmap[i] = Gmap[i];
	qsort(nmap,/* WAS Cfirst*/ nummaps,sizeof(Cmap *),(intcmp*)CmapLenDec);

	if(VERB>=2){
	  printf("Sorted %d/%d remaining maps in descending order of size\n", /* WAS Cfirst */ nummaps, nummaps);
	  fflush(stdout);
	}

	/* First compute new id number of remaining maps to be output : nmap[0..nummaps-1] with paired==0 */
	for(int i = 0; i < nummaps; i++){
	  Cmap *pmap = nmap[i];
	  if(DEBUG) assert(pmap == Gmap[pmap->mapid]);
	  newid[pmap->mapid] = -1;
	  if(VERB>=2){
	    int M = pmap->numsite[0];
	    printf("nmap[%d]=%p : mapid=%d,id=%lld,paired=%d,len=%0.4f,M=%d, origmap=%p, nextid=%lld\n",
		   i, pmap, pmap->mapid, pmap->id, pmap->paired, pmap->site[0][M+1], M, pmap->origmap, nextid);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(!pmap->origmap);

	  if(pmap->paired)
	    continue;
	  newid[pmap->mapid] = nextid++;
	}
	delete [] nmap;

	/* output mapping from previously merged maps MergedMaps[0..MergedNummaps-1]->id to new output map ids */
	for(int i = 0; i < MergedNummaps; i++){
	  Cmap *pmap = MergedMaps[i];
	  if(DEBUG) assert(pmap->paired);

	  /* first follow origmap links to one of the remaining maps (map[0..nummaps-1]) */
	  Cmap *qmap = pmap->origmap;
	  if(DEBUG) assert(qmap);
	  int cnt = 0;
	  if(VERB>=2){
	    int M = pmap->numsite[0];
	    printf("MergedMaps[%d]=%p : mapid=%d,id=%lld,paired=%d,len=%0.4f,M=%d, origmap=%p\n",
		   i, pmap, pmap->mapid, pmap->id, pmap->paired, pmap->site[0][M+1], M, pmap->origmap);
	    fflush(stdout);
	  }

	  while(qmap->origmap && cnt < MergedNummaps){
	    qmap = qmap->origmap;
	    if(VERB >=2){
	      int N = qmap->numsite[0];
	      printf("\t qmap=%p : mapid=%d,id=%lld,paired=%d,len=%0.4f,N=%d,origmap=%p\n",
		       qmap,qmap->mapid,qmap->id,qmap->paired,qmap->site[0][N+1],N,qmap->origmap);
	      fflush(stdout);
	    }
	    cnt++;
	  }
	  if(DEBUG) assert(qmap->origmap == NULL);

	  if(DEBUG && !(0 <= qmap->mapid && qmap->mapid < nummaps)){
	    int M = pmap->numsite[0];
	    printf("MergedMaps[%d]=%p : mapid=%d,id=%lld,paired=%d,len=%0.4f,M=%d, origmap=%p\n",
		   i, pmap, pmap->mapid, pmap->id, pmap->paired, pmap->site[0][M+1], M, pmap->origmap);

	    int N = qmap->numsite[0];
	    printf("\t qmap=%p : mapid=%d,id=%lld,paired=%d,len=%0.4f,N=%d,origmap=%p\n",
		   qmap,qmap->mapid,qmap->id,qmap->paired,qmap->site[0][N+1],N,qmap->origmap);
	    fflush(stdout);

	    assert(0 <= qmap->mapid && qmap->mapid < nummaps);
	  }
	  
	  /* next traverse paired links across merges in current iteration to (largest) output map */
	  int outid = qmap->mapid;
	  if(DEBUG) assert(qmap == Gmap[outid]);
	  cnt = 0;
	  while(qmap->paired && cnt < nummaps){
	    qmap = Gmap[outid = qmap->paired - 1];
	    if(VERB>=2){
	      int M = qmap->numsite[0];
	      printf("\t qmap=%p : mapid=%d,id=%lld,paired=%d,outid=%d,len=%0.4f,M=%d,origmap=%p\n",
		     qmap,qmap->mapid,qmap->id,qmap->paired,outid,qmap->site[0][M+1],M,qmap->origmap);
	      fflush(stdout);
	    }
	    if(DEBUG) assert(qmap->origmap == NULL);
	    cnt++;
	  }
	  if(DEBUG) assert(qmap->paired == 0);

	  if(VERB>=2 || !(qmap == Gmap[outid])){
	    int M = qmap->numsite[0];
	    printf("\t final qmap=%p : mapid=%d,id=%lld,paired=%d,len=%0.4f,M=%d: newid[outid=%d]=%lld, pmap->id=%lld\n",
		   qmap, qmap->mapid,qmap->id,qmap->paired,qmap->site[0][M+1],M,outid,newid[outid],pmap->id);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(qmap == Gmap[outid]);
	  if(DEBUG) assert(newid[outid] > 0);
	  if(DEBUG) assert(CIDcnt < MergedNummaps + nummaps);
	  CIDmap[CIDcnt].newID = newid[outid];
	  CIDmap[CIDcnt++].oldID = pmap->origid;
	  //	  fprintf(fpID,"%lld %lld\n",newid[outid], pmap->id);
	}

	/* output mapping from remaining maps Gmap[0..nummaps-1] to new output map ids */
	for(int i = 0; i < nummaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(DEBUG) assert(!pmap->origmap);

	  /* traverse paired links across merges in current iteration to (largest) output map */
	  Cmap *qmap = pmap;
	  int outid = i, cnt = 0;
	  while(qmap->paired && cnt < nummaps){
	    qmap = Gmap[outid = qmap->paired - 1];
	    if(VERB>=2){
	      int M = qmap->numsite[0];
	      printf("\t qmap=%p : mapid=%d,id=%lld,paired=%d,outid=%d,len=%0.4f,M=%d,origmap=%p\n",
		     qmap,qmap->mapid,qmap->id,qmap->paired,outid,qmap->site[0][M+1],M,qmap->origmap);
	      fflush(stdout);
	    }
	    if(DEBUG) assert(qmap->origmap == NULL);
	    cnt++;
	  }
	  if(VERB>=2){
	    int M = pmap->numsite[0];
	    printf("map[%d]=%p : mapid=%d,id=%lld,paired=%d,outid=%d,len=%0.4f,M=%d: newid[%d]=%lld, pmap->id=%lld\n",
		   i, pmap, pmap->mapid,pmap->id,pmap->paired,outid,pmap->site[0][M+1],M,outid,newid[outid],pmap->id);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(qmap->paired == 0);
	  if(DEBUG) assert(qmap == Gmap[outid]);
	  if(DEBUG) assert(newid[outid] > 0);
	  if(DEBUG) assert(CIDcnt < MergedNummaps + nummaps);
	  CIDmap[CIDcnt].newID = newid[outid];
	  CIDmap[CIDcnt++].oldID = pmap->origid;
	  //	  fprintf(fpID,"%lld %lld\n", newid[outid], pmap->id);
	}

	if(DEBUG) assert(CIDcnt <= MergedNummaps + nummaps);

	qsort(CIDmap,CIDcnt,sizeof(CID),(intcmp*)CIDinc);

	/* remove duplicates caused by intermediate nodes in the merge graph */
	int j = 0;
	for(int i = 1; i < CIDcnt; i++){
	  if(VERB>=3){
	    printf("CIDmap[j=%d]:newID=%lld,oldID=%lld; CIDmap[i=%d]:newID=%lld,oldID=%lld\n",
		   j,CIDmap[j].newID,CIDmap[j].oldID,i,CIDmap[i].newID,CIDmap[i].oldID);
	    fflush(stdout);
	  }
	  if(CIDmap[i].newID == CIDmap[j].newID && CIDmap[i].oldID == CIDmap[j].oldID)
	    continue;
	  CIDmap[++j] = CIDmap[i];
	}
	if(VERB/* HERE >=2 */){
	  printf("Reduced number of ID mappings from %d to %d\n",CIDcnt,min(CIDcnt,j+1));
	  fflush(stdout);
	}
	CIDcnt = min(CIDcnt,j+1);

	for(int i = 0; i < CIDcnt; i++)
	  fprintf(fpID,"%lld %lld\n",CIDmap[i].newID, CIDmap[i].oldID);

	FILEclose(fpID);

	if(VERB){
	  printf("Renumbered output contig ids from 1 to %lld (mapping from original to new IDs in %s)\n",nextid-1,IDmap);
	  fflush(stdout);
	}
	
	/* update output map and hmap IDs to new IDs */
	for(int i = 0; i < nummaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(DEBUG) assert(pmap->origmap == NULL);
	  if(pmap->paired){
	    if(VERB>=2){
	      printf("map[%d]: id = %lld, mapid=%d, paired=%d (skipping)\n",i, pmap->id, pmap->mapid, pmap->paired);
	      fflush(stdout);
	    }
	    continue;
	  }
	  if(DEBUG) assert(newid[i] > 0);
	  if(VERB>=2){
	    printf("map[%d]: id = %lld -> %lld, mapid=%d\n",i, pmap->id, newid[i], pmap->mapid);
	    fflush(stdout);
	  }
	  pmap->id = newid[i];
	  if(pmap->contig && pmap->contig->nummaps)
	    pmap->contig->id = pmap->id;
	}

	delete [] newid;
	delete [] CIDmap;
      }

      int writecnt = 0;

      // HERE HERE  #pragma omp parallel for schedule(dynamic,1) num_threads(gnumthreads)
      for(int i = 0; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	if(DEBUG>=2 && i >= orignummaps) assert(!pmap->paired);
	if(DEBUG) assert(pmap->origmap == NULL);
	if(pmap->paired){
	  if(VERB>=2){
	    printf("Skipping output of Gmap[%d]: paired=%d\n",i,pmap->paired);
	    fflush(stdout);
	  }
	  continue;
	}

	char prefix[PATH_MAX];
	sprintf(prefix,"%s_contig%lld",output_prefix, pmap->id);

	if(DEBUG && PairMergeHmap > 0 && !(pmap->centercloned & 2) && pmap->contig){
	  int N = pmap->numsite[0];
	  printf("i=%d/%d: mapid=%d, id=%lld, len= %0.3f, N=%d, paired=%d, centercloned=%d\n",i,nummaps,pmap->mapid,pmap->id, 
		 pmap->site[0][N+1],N, pmap->paired,pmap->centercloned);
	  fflush(stdout);

	  assert(!(PairMergeHmap > 0 && !(pmap->centercloned & 2) && pmap->contig));
	}

	/* allow unmodified Mask values for unmodified contig to be linked to original file (if !PairMergeIDrenumber) OR have undo mres (otherwise). 
	   NOTE : Any modification to Mask should have set pmap->centercloned */

	if(i < orignummaps && !pmap->centercloned){/* unchanged contig map without Mask */
          if(DEBUG && PairMergeHmap > 0) assert(!pmap->contig);

#ifndef WIN32
	  if(SplitMap && !PairMergeIDrenumber && !strstr(prefix,"/dev/null")){ /* create a link to the original file */
            if(DEBUG) assert(pmap->fileid < num_files);
	    char filename[PATH_MAX];
	    sprintf(filename,"%s.cmap",prefix);	
	    
	    FILE *fp = fopen(filename,"r");
	    if(ForceOverwrite && fp != NULL){
              #pragma omp critical
	      {
	        printf("WARNING:Output file %s already exists\n",filename);
		fflush(stdout);
		fclose(fp);

		/* unlink the filename */
		errno = 0;
		if(unlink(filename)){
	          int eno = errno;
		  char *err = strerror(eno);
		  printf("unlink(%s) failed:errno=%d:%s\n",filename,eno,err);
		  fflush(stdout);exit(1);
                }
	      }
            }
	      
            errno = 0;
	    if(!link(vfixx_filename[pmap->fileid],filename)){
              if(VERB/* >=2 */){
                printf("created new link %s to %s (errno=%d)\n",filename,vfixx_filename[pmap->fileid],errno);
		fflush(stdout);
              }
	      linkcnt++;
	      continue;
            }

	    int eno = errno;
	    char *err = strerror(eno);
            fprintf(stderr,"link(%s,%s) failed: %s\n",vfixx_filename[pmap->fileid],filename,err);
	    fflush(stderr);
	  }
#endif

	  if(mres > 0.0){/* restore original map resolution */
	    for(int c=0; c < colors; c++){
	      if(pmap->origsite[c]==NULL)
		continue;/* no original map resolution : must have already been restored */

	      /* free up regular arrays */
	      if(!pmap->blockmem)
		delete [] pmap->site[c]; 
	      pmap->site[c] = NULL;
	      if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
	      if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
	      if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
	      if(CmapChimQuality >= 1){
		if(pmap->ChimQuality[c]) delete [] pmap->ChimQuality[c];
		if(CmapChimQuality >= 2){
		  if(pmap->ChimNorm[c]) delete [] pmap->ChimNorm[c];
		  if(pmap->SegDupL[c]) delete [] pmap->SegDupL[c];
		  if(pmap->SegDupR[c]) delete [] pmap->SegDupR[c];
		  if(pmap->FragileEndL[c]) delete [] pmap->FragileEndL[c];
		  if(pmap->FragileEndR[c]) delete [] pmap->FragileEndR[c];
		  if(pmap->OutlierFrac[c]) delete [] pmap->OutlierFrac[c];
		  if(CmapChimQuality >= 3){
		    if(pmap->FragSd[c]) delete [] pmap->FragSd[c];
		    if(pmap->ExpSd[c]) delete [] pmap->ExpSd[c];
		    if(pmap->FragCov[c]) delete [] pmap->FragCov[c];
		    if(pmap->FragChiSq[c]) delete [] pmap->FragChiSq[c];
		  }
		}
	      }
	      if(pmap->Mask[c]){ delete [] pmap->Mask[c]; pmap->Mask[c] = NULL; }
	      if(pmap->SNRcnt[c]){
		if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		if(pmap->SNRdist[c]){
		  for(int t = 0; t <= pmap->numsite[c]+1; t++)
		    if(pmap->SNRdist[c][t]){
		      if(VERB>=3){
			printf("Calling delete[] Gmap[%d](%p)->SNRdist[%d][t=%d]=%p -> 0,blockmem=%d\n",i,pmap,c,t,pmap->SNRdist[c][t],pmap->blockmem);
			fflush(stdout);
		      }
		  
		      delete [] pmap->SNRdist[c][t];
		      pmap->SNRdist[c][t] = NULL;
		    }
		  delete [] pmap->SNRdist[c];
		  pmap->SNRdist[c] = NULL;
		}
		delete [] pmap->SNRcnt[c];
		pmap->SNRcnt[c] = NULL;
	      }

	      pmap->numsite[c] = pmap->orignumsite[c];
	      if(DEBUG) assert(pmap->origsite[c]);

	      /* just copy pointers */
	      pmap->blockmem = -1;
	      pmap->site[c] = pmap->origsite[c];
	      pmap->siteSD[c] = pmap->origsiteSD[c];
	      pmap->sitecov[c] = pmap->origsitecov[c];
	      if(CmapChimQuality >= 1){
		if(pmap->origChimQuality[c]) pmap->ChimQuality[c] = pmap->origChimQuality[c];
		if(CmapChimQuality >= 2){
		  if(pmap->origChimNorm[c]) pmap->ChimNorm[c] = pmap->origChimNorm[c];
		  if(pmap->origSegDupL[c]) pmap->SegDupL[c] = pmap->origSegDupL[c];
		  if(pmap->origSegDupR[c]) pmap->SegDupR[c] = pmap->origSegDupR[c];
		  if(pmap->origFragileEndL[c]) pmap->FragileEndL[c] = pmap->origFragileEndL[c];
		  if(pmap->origFragileEndR[c]) pmap->FragileEndR[c] = pmap->origFragileEndR[c];
		  if(pmap->origOutlierFrac[c]) pmap->OutlierFrac[c] = pmap->origOutlierFrac[c];
		  if(CmapChimQuality >= 3){
		    if(pmap->origFragSd[c]) pmap->FragSd[c] = pmap->origFragSd[c];
		    if(pmap->origExpSd[c]) pmap->ExpSd[c] = pmap->origExpSd[c];
		    if(pmap->origFragCov[c]) pmap->FragCov[c] = pmap->origFragCov[c];
		    if(pmap->origFragChiSq[c]) pmap->FragChiSq[c] = pmap->origFragChiSq[c];
		  }
		}
	      }
	      if(pmap->origMask[c]) pmap->Mask[c] = pmap->origMask[c];

	      pmap->sitecnt[c] = pmap->origsitecnt[c];
	      pmap->SNRcnt[c] = pmap->origSNRcnt[c];
	      pmap->SNRdist[c] = pmap->origSNRdist[c];

	      pmap->SNRgmean[c] = pmap->origSNRgmean[c];
	      pmap->lnSNRsd[c] = pmap->origlnSNRsd[c];
	      if(pmap->SNRcnt[c])
		for(int t = 0; t <= pmap->orignumsite[c]+1; t++){
		  if(VERB>=3){
		    printf("map[%d]->SNRdist[%c][%d]= %p -> %p\n",i,c,t,pmap->SNRdist[c][t],pmap->origSNRdist[c][t]);
		    fflush(stdout);
		  }

		  pmap->SNRdist[c][t] = pmap->origSNRdist[c][t];
		}

	      /* reset backup pointers */
	      pmap->origsite[c] = 0;
	      pmap->origsiteSD[c] = 0;
	      pmap->origsitecov[c] = 0;
	      pmap->origsitecnt[c] = 0;
	      if(CmapChimQuality >= 1){
		pmap->origChimQuality[c] = NULL;
		if(CmapChimQuality >= 2){
		  pmap->origChimNorm[c] = NULL;
		  pmap->origSegDupL[c] = NULL;
		  pmap->origSegDupR[c] = NULL;
		  pmap->origFragileEndL[c] = NULL;
		  pmap->origFragileEndR[c] = NULL;
		  pmap->origOutlierFrac[c] = NULL;
		  if(CmapChimQuality >= 3){
		    pmap->origFragSd[c] = NULL;
		    pmap->origExpSd[c] = NULL;
		    pmap->origFragCov[c] = NULL;
		    pmap->origFragChiSq[c] = NULL;
		  }
		}
	      }
	      pmap->origMask[c] = NULL;
	      if(pmap->origSNRcnt[c]){
		pmap->origSNRcnt[c] = 0;
		pmap->origSNRdist[c] = 0;
		pmap->origSNRgmean[c] = 0;
		pmap->origlnSNRsd[c] = 0;
	      }
	    }
	  }
	}

	if(pmap->contig && pmap->contig->nummaps <= 1){
	  if(VERB>=2){
	    printf("map[%d] : contig->nummaps= %d (mapid=%d,id=%lld, contig->id=%lld) : converting Hmap to regular map\n",
		   i,pmap->contig->nummaps,pmap->mapid,pmap->id,pmap->contig->id);
	    fflush(stdout);
	  }
	  delete pmap->contig;
	  pmap->contig = 0;
	}

	if(VERB/* HERE >=2 */){
          #pragma omp critical
          {
	    int N = pmap->numsite[0];
	    if(i < orignummaps && !pmap->centercloned)
	      printf("writing out original Gmap[%d]=%p: mapid=%d,paired=%d,id=%lld,len=%0.3f,N=%d as %s.cmap\n",i,pmap,pmap->mapid,pmap->paired,pmap->id,pmap->site[0][N+1],N,prefix);
	    else if(!(pmap->contig && pmap->contig->nummaps > 1))
	      printf("writing out merged Gmap[%d]=%p: mapid=%d,paired=%d,id=%lld,len=%0.3f,N=%d as %s.cmap\n",i,pmap,pmap->mapid,pmap->paired,pmap->id, pmap->site[0][N+1],N,prefix);
	    else {
	      int N = pmap->contig->numsite[0];
	      printf("writing out merged Gmap[%d]=%p: mapid=%d,paired=%d,id=%lld,len=%0.3f,N=%d (nummaps=%d) as HaploTyped %s.cmap\n",
		     i,pmap,pmap->mapid,pmap->paired,pmap->id,pmap->contig->site[0][N+1],N,pmap->contig->nummaps, prefix);
	    }
	    fflush(stdout);
          }
	}

	#pragma omp atomic
	writecnt++;

	if(pmap->contig && pmap->contig->nummaps > 1){
          if(DEBUG) assert(PairMergeHmap > 0);
	  output_hmap(prefix, Gmap, i, i, 0);
        } else
	  output_cmap(prefix, Gmap, i, i);
      }

      if(VERB){
	double wt = wtime();
        printf("pairmerge: iter=%d: %d Files written (and %d linked): wall time=%0.6f (cum=%0.6f)\n",pairmergeIter, writecnt, linkcnt, wt - origwt, wt);
	fflush(stdout);
      }
    } else {/* prepare for pairalign_allpairs() to be called again */
      /* NOTE : score_free must be called here before the arrays are rearranged */
      score_free(XXmap,numX);XXmap = 0;numX= 0;
      score_free(YYmap,numY);YYmap = 0;numY= 0;

      if(PairMergeIDrenumber){ /* (re)allocate array MergedMaps[0..MergedNummaps + orignummaps - 1] */
        if(DEBUG>=2){
          for(int i = 0; i < MergedNummaps; i++){
	    if(!MergedMaps[i]->paired || !MergedMaps[i]->origmap){
              printf("pairmergeIter=%d,i=%d,MergedNummaps=%d/%d:MergedMaps[i]:paired=%d,id=%lld\n",
	        pairmergeIter,i, MergedNummaps,MergedMaxmaps,MergedMaps[i]->paired,MergedMaps[i]->id);
	      if(MergedMaps[i]->origmap)
		printf("\t MergedMaps[i]->origmap:mapid=%d,id=%lld,paired=%d\n",MergedMaps[i]->origmap->mapid,MergedMaps[i]->origmap->id,MergedMaps[i]->origmap->paired);
	      fflush(stdout);
              assert(MergedMaps[i]->paired);
              assert(MergedMaps[i]->origmap);
            }
          }
        }
        if(VERB>=2 && MergedNummaps > 45){
	  Cmap *pmap = MergedMaps[45];
	  printf("MergedNummaps=%d/%d:MergedMaps[45]=%p:paired=%d,mapid=%d,id=%lld\n",
	    MergedNummaps,MergedMaxmaps,pmap,pmap->paired,pmap->mapid,pmap->id);
	  fflush(stdout);
        }

        maxmapalloc(MergedNummaps + orignummaps, MergedMaxmaps, MergedMaps, 0, 0);

        if(DEBUG>=2){
          for(int i = 0; i < MergedNummaps; i++){
	    if(!MergedMaps[i]->paired || !MergedMaps[i]->origmap){
              printf("pairmergeIter=%d,i=%d,MergedNummaps=%d/%d:MergedMaps[i]:paired=%d,id=%lld\n",
	        pairmergeIter,i, MergedNummaps,MergedMaxmaps,MergedMaps[i]->paired,MergedMaps[i]->id);
	      if(MergedMaps[i]->origmap)
		printf("\t MergedMaps[i]->origmap:mapid=%d,id=%lld,paired=%d\n",MergedMaps[i]->origmap->mapid,MergedMaps[i]->origmap->id,MergedMaps[i]->origmap->paired);
	      fflush(stdout);
              assert(MergedMaps[i]->paired);
              assert(MergedMaps[i]->origmap);
            }
          }
        }
        if(VERB>=2 && MergedNummaps > 45){
	  Cmap *pmap = MergedMaps[45];
	  printf("MergedNummaps=%d/%d:MergedMaps[45]=%p:paired=%d,mapid=%d,id=%lld\n",
	    MergedNummaps,MergedMaxmaps,pmap,pmap->paired,pmap->mapid,pmap->id);
	  fflush(stdout);
        }
      }

      int j = 0;
      //      int newCfirst = -1;
      for(int i = 0; i < orignummaps; i++){
	Cmap *pmap = Gmap[i];
	if(pmap->paired){
          if(PairMergeIDrenumber){/* append pmap to MergedMaps[0..MergeNummaps-1] */
	    pmap->origmap = Gmap[pmap->paired - 1];
	    if(VERB>=2){
/*	      printf("i=%d/%d:MergedNummaps=%d/%d:MergedMaps[%d]= %p -> (paired=%d,mapid=%d,id=%lld,origmap=%p)\n",
		      i,orignummaps,MergedNummaps,MergedMaxmaps,MergedNummaps,pmap,pmap->paired,pmap->mapid,pmap->id,pmap->origmap);*/
	      printf("map[%d] = MergedMaps[%d] -> (paired=%d,mapid=%d,id=%lld)\n",
		     i, MergedNummaps,pmap->paired,pmap->mapid,pmap->id);
	      fflush(stdout);
	    }
	    if(DEBUG) assert(0 < pmap->paired && pmap->paired <= nummaps);
	    if(DEBUG) assert(MergedNummaps < MergedMaxmaps);
	    MergedMaps[MergedNummaps++] = pmap;
          }
        }
      }
      for(int i = 0; i < orignummaps; i++){
	Cmap *pmap = Gmap[i];
	if(pmap->paired){
          if(PairMergeIDrenumber)
	    Gmap[i] = NULL;
	  continue;
        }
	if(j < i){
          Cmap *tmp = Gmap[j];
          Gmap[j] = Gmap[i];
	  Gmap[i] = tmp;
	}
	if(PairMergeIDrenumber){
	  if(VERB>=2){
	    int M = Gmap[j]->numsite[0];
	    if(j < i)
	      printf("map[%d]=%p: copied from Gmap[%d]: mapid=%d->%d,id=%lld,paired=%d,len=%0.4f,M=%d\n",
	        j, Gmap[j], i, Gmap[j]->mapid, j, Gmap[j]->id, Gmap[j]->paired, Gmap[j]->site[0][M+1], M);
	    else
	      printf("map[%d]=%p : mapid=%d,id=%lld,paired=%d,len=%0.4f,M=%d\n",
	         j, Gmap[j], Gmap[j]->mapid, Gmap[j]->id, Gmap[j]->paired, Gmap[j]->site[0][M+1], M);
	    fflush(stdout);
	  }
	  Gmap[j]->mapid = j;
	}
	j++;
      }
      Cfirst = j;
      for(int i = orignummaps; i < nummaps; i++){
        if(DEBUG) assert(Gmap[i]->paired==0);
	Gmap[i]->paired = 1;// marks map as newly merged : in next round all pairs must include at least one such map
	if(j < i){
	  Cmap *tmp = Gmap[j];
	  Gmap[j] = Gmap[i];
	  Gmap[i] = tmp;
	}
	if(PairMergeIDrenumber){
	  if(VERB>=2){
	    int M = Gmap[j]->numsite[0];
	    if(j < i)
	      printf("Gmap[%d]=%p: copied from Gmap[%d]: mapid=%d->%d,id=%lld,paired=0->%d,len=%0.4f,M=%d\n",
	        j, Gmap[j], i, Gmap[j]->mapid, j, Gmap[j]->id, Gmap[j]->paired, Gmap[j]->site[0][M+1], M);
	    else
	      printf("map[%d]=%p : mapid=%d,id=%lld,paired=0->%d,len=%0.4f,M=%d\n",
	         j, Gmap[j], Gmap[j]->mapid, Gmap[j]->id, Gmap[j]->paired, Gmap[j]->site[0][M+1], M);
	    fflush(stdout);
	  }
	  Gmap[j]->mapid = j;
        }
	j++;
      }
      nummaps = j;
      if(DEBUG>=1+RELEASE) assert(nummaps == fmaps);
      int unmodified = Cfirst;
      if(!Cfirst && PairMergeProgress){// NEW : rare case when there are no un-modified maps left
	if(DEBUG) assert(nummaps > 0);
	Cfirst = 1;
	//	map[0]->paired = 0;// is this needed ?
      }

      CfirstPairs = 2;
      if(VERB/* HERE >=2 */){
	printf("PairMergeRepeat:Cfirst -> %d (original unmodified maps=%d), nummaps= %d -> %d, numaligns=%lu\n",Cfirst,unmodified,orignummaps, nummaps,numaligns);
	if(PairMergeIDrenumber)
	  printf("\t Saved MergedNummaps=%d/%d\n",MergedNummaps,MergedMaxmaps);
	fflush(stdout);
      }
      if(DEBUG>=2){
        for(int i = 0; i < MergedNummaps; i++){
          if(!MergedMaps[i]->paired || !MergedMaps[i]->origmap){
            printf("pairmergeIter=%d,i=%d,MergedNummaps=%d/%d:MergedMaps[i]=%p:paired=%d,id=%lld\n",
	      pairmergeIter,i, MergedNummaps,MergedMaxmaps,MergedMaps[i],MergedMaps[i]->paired,MergedMaps[i]->id);
	    if(MergedMaps[i]->origmap)
	      printf("\t MergedMaps[i]->origmap:mapid=%d,id=%lld,paired=%d\n",MergedMaps[i]->origmap->mapid,MergedMaps[i]->origmap->id,MergedMaps[i]->origmap->paired);
	    fflush(stdout);
            assert(MergedMaps[i]->paired);
            assert(MergedMaps[i]->origmap);
          }
        }
      }
    } // if((PairMergeRepeat && nummaps > orignummaps) || (!SplitSegDup && origSplitSegDup) || PairMergeHmap < 0)
  }

  if(PairsplitSeparateQueries)/* reset copied array pointers to avoid call free() on them */
    for(int Yid = orignumY; Yid < orignumYnumX; Yid++)
      YYmap[Yid]->init();

  if(LEAKDEBUG || PairMerge || RA_MIN_TIME){
    if(!PairMergeProgress){
      score_free(XXmap,numX);XXmap = 0;numX= 0;
      score_free(YYmap,numY);YYmap = 0;numY= 0;
    }

    delete [] misaligns;
    delete [] errors;

    if(hashpairs1){
      delete [] hashpairs1;
      hashpairs1 = NULL;
      numhashpairs1 = maxhashpairs1 = 0;
    }
    if(hashpairs2){
      delete [] hashpairs2;
      hashpairs2 = NULL;
      numhashpairs2 = maxhashpairs2 = 0;
    }
    if(Ymapid2id){
      delete [] Ymapid2id;
      Ymapid2id = 0;
    }
    if(Cfirst && Xmapid2id){
      delete [] Xmapid2id;
      Xmapid2id = 0;
    }

    if(PairSplit)
      PairHash_free();

    if(A_threads) {
      /* free memory of *thread_arrays* */
      for(int i=0;i<A_numthreads;i++) {
	delete [] A_threads[i].score;
	delete [] A_threads[i].Uscore;
	delete [] A_threads[i].G;
	delete [] A_threads[i].H;
	delete [] A_threads[i].Rijy;
	delete [] A_threads[i].Lijy;
	delete [] A_threads[i].Rijx;
	delete [] A_threads[i].Lijx;
      }
      delete [] A_threads;
    }
    A_threads = NULL;

    if(ThreadMem){
      if(RA_MIN_TIME){/* use munmap to free Fmem/Imem for all threads */
	for(int i=0; i < mem_numthreads; i++){
	  if(VERB>=2){
#pragma omp critical
	    {
	      printf("tid=%d: Freeing Fmem=%p,Imem=%p, MemSiz= %lld (RA_MIN_MEM=%d,RA_MIN_TIME=%d)\n",i,ThreadMem[i].Fmem,ThreadMem[i].Imem,ThreadMem[i].MemSiz,RA_MIN_MEM,RA_MIN_TIME);
	      fflush(stdout);
	    }
	  }
          if(ThreadMem[i].Fmem && munmap(ThreadMem[i].Fmem, ThreadMem[i].MemSiz)){
	    int eno = errno;
	    char *err = strerror(eno);
	    printf("munmap(%p,%lld) failed (i=%d): errno=%d:%s\n",ThreadMem[i].Fmem, ThreadMem[i].MemSiz, i, eno, err);
	    dumpmemmap();
	    fflush(stdout);exit(1);
	  }
	}
	if(VERB>=2){
	  long long VmSize,VmRSS,VmSwap,VmHWM;// local copies, so global shared values are not modified
	  getmem(VmSize,VmRSS,VmSwap,&VmHWM);

	  printf("After freeing Fmem[] memory: VmRSS= %0.4f(HWM= %0.4f), VmSwap= %0.4f, VmSize= %0.4f Gb: wall time= %0.6f secs\n",VmRSS*1e-9,VmHWM*1e-9,VmSwap*1e-9,VmSize*1e-9,wtime());
	  fflush(stdout);
	}
      } else {// if(RA_MIN_TIME == 0)
	for(int i=0; i < mem_numthreads; i++){
	  delete [] ThreadMem[i].Fmem;
	  delete [] ThreadMem[i].Imem;
	}
      }
      for(int i = 0; i < mem_numthreads; i++){
	ThreadMem[i].Fmem = NULL;
	ThreadMem[i].Imem = NULL;
      }
    }

    if(THREAD_DEBUG){
      double wt = wtime();
      if(VERB && wt > wtstart){
	size_t madvise_totcalls = 0, getmem_totcalls = 0, getstatm_totcalls = 0, critical_totcalls = 0, atomic_totcalls=0, wait_totcalls=0,release_totcalls = 0;
	for(int i = 0; i < A_numthreads; i++){
	  madvise_totcalls += ThreadMem[i].madvise_calls;
	  getmem_totcalls += ThreadMem[i].getmem_calls;
	  getstatm_totcalls += ThreadMem[i].getstatm_calls;
	  critical_totcalls += ThreadMem[i].critical_calls;
	  atomic_totcalls += ThreadMem[i].atomic_calls;
	  wait_totcalls += ThreadMem[i].wait_calls;
	  release_totcalls += ThreadMem[i].release_calls;
	}
	printf("calls:madvise= %lu(%0.1f/sec), getmem= %lu(%0.1f/sec), getstatm= %lu(%0.1f/sec) critical= %lu(%0.1f/sec), atomic= %lu(%0.f/sec), wait=%lu,release=%lu: elapsed time= %0.6f (wt=%0.6f)\n",
	       madvise_totcalls, madvise_totcalls/(wt-wtstart), getmem_totcalls, getmem_totcalls/(wt-wtstart), getstatm_totcalls, getstatm_totcalls/(wt-wtstart), 
	       critical_totcalls, critical_totcalls/(wt-wtstart),  atomic_totcalls, atomic_totcalls/(wt-wtstart), wait_totcalls, release_totcalls, wt-wtstart,wt);
	fflush(stdout);
      }
    }

    if(ThreadMem){
      delete [] ThreadMem;
      ThreadMem = 0;
    }

    A_numthreads = mem_numthreads = 0;
    A_maxN = mem_maxN = mem_maxM = 0;
  }
  return PairMergeProgress;
}
