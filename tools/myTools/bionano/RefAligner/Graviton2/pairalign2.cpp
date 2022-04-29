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

// #include <omp.h>

#include "globals.h"
#include "parameters.h"
#include "Calign.h"
#include "hash.h"

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
static PFLOAT aNaN = 0/denom;
#else
static PFLOAT aNaN = nan("NaN");
#endif

//#undef DEBUG
//#define DEBUG 2

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/pairalign2.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

#define OUTLIER_FIXES PairwiseOutlierFix /* enable various fixes to scoring of outliers (see command line parameter -PairwiseOutlierFix, on by default) */

#include "NGentigPairScore2.h"

// This scoring function works poorly with very low FP/FN, eg FP=0.2, FN=0.04, sf=0, sd=0.04 and -outlier 0.0001 (see E.Coli alignments pair2.fixx) : find out why (possible error in seperating Pen and Bias, or
//    use of fixed K values independent of FP rate)

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

#define EVERB 0 /* (((YYmap[Yid]->id == 3766LL && XXmap[Xid]->id == 19LL)||(YYmap[Yid]->id == 19LL && XXmap[Xid]->id == 3766LL)) ? 1 : 0) */ /* Display details of all alignments, even those NOT above threshold */

#define LVERB 0 /* Display details of multiple local alignment (implies PVERB==1) */

#define FLAT(A, field, I, J)	((A) -> field [1][(I-1) * (A)->M_stride + J])
#define FLAT_DELX(ptr, n, J, M)  (ptr[1][J+(n-1)*COMPUTE_STRIDE4(M-1)])

#if USE_MIC
#define memcpy _mm512_memcpy
#endif

extern double mtime();/* microsecond timer function : see refine.cpp */
extern double wtime();

extern Cmap **YYmap,**XXmap;
extern int numY,numX;

static int repeatcnt = 0;/* global count of alignments that were supressed due to -RepeatMask */

#if 0
/* sort alignments in order of id2, id1, orientation, sites1[0], sites1[U-1] (null alignments last) */
static int IdOrSiteInc(register Calign **p1, register Calign **p2)
{
  register Calign *align1 = p1[0];
  register Calign *align2 = p2[0];
  if(!align1) return 1;
  if(!align2) return -1;

  register Cmap *Ymap1 = YYmap[align1->mapid1];
  register Cmap *Ymap2 = YYmap[align2->mapid1];
  register Cmap *Xmap1 = XXmap[align1->mapid2];
  register Cmap *Xmap2 = XXmap[align2->mapid2];

  int ret = (Xmap1->id > Xmap2->id) ? 1 : (Xmap1->id < Xmap2->id) ? -1 :
    (Ymap1->id > Ymap2->id) ? 1 : (Ymap1->id < Ymap2->id) ? -1 :
    (align1->orientation > align2->orientation) ? 1 : (align1->orientation < align2->orientation) ? -1 : 0;
  if(ret)
    return ret;

  register int left1 = align1->sites1 ? (align1->sites1[0] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  register int left2 = align2->sites2 ? (align2->sites1[0] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  ret = (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  int numpairs1 = align1->numpairs;
  int numpairs2 = align2->numpairs;
  int right1 = align1->sites1 ? (align1->sites1[numpairs1-1] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int right2 = align2->sites2 ? (align2->sites1[numpairs2-1] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  return (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
}
#endif

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

static int maxerrors[MAXCOLOR+1], numerrors[MAXCOLOR+1];
static Cinterval *errors[MAXCOLOR+1];

static int maxmisaligns[MAXCOLOR], nummisaligns[MAXCOLOR];
static Cmisalign *misaligns[MAXCOLOR];

static int xinc(Cmisalign *p1, Cmisalign *p2)
{
  return (p1->x > p2->x) ? 1 : (p1->x < p2->x) ? -1 : 0;
}

static int CmapIdInc(Cmap **p1, Cmap **p2)
{
  /* Cannot use p1[0]->id - p2[0]->id since that may be out of range for integer */
  return (p1[0]->id > p2[0]->id) ? 1 : (p1[0]->id < p2[0]->id) ? -1 : 0;
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

/** Array based Dynamic programming data */
class AL {
public:
  PFLOAT **score;
  PFLOAT **Uscore;
  PFLOAT **Rscore;
  int **G;
  int **H;
  int **Rijy;
  int **Lijy;
  int **Rijx;
  int **Lijx;
  long long M_stride;
  long long array_stride;
};

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

#if 0
/** take existing alignment and return pvalue logFP() */
static double alignFP(int c, Calign *p, FLOAT *Y, FLOAT *X, int N, int M,
		      int orientation, int Yid, int Xid, double lscore, int verb)
{
  int U = p->numpairs;
  
  /* allocate alignment fragment sizes */
  FLOAT *Xfrag = (FLOAT *) alloca((U-1)*sizeof(FLOAT));
  FLOAT *Yfrag = (FLOAT *) alloca((U-1)*sizeof(FLOAT));
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
	if(LVERB>=2)
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
      if(LVERB>=2)
	printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:F=%d:y=%0.3f,x=%0.3f,misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	       T,U,I,J,G,H,F,Y[I]-Y[G],X[J]-X[H],misscnt,misscnt+(I-G-1),falsecnt,falsecnt + (J-H-1),OutlierCnt,OutlierMN);
      Yfrag[F] = Y[I] - Y[G];
      Xfrag[F] = X[J] - X[H];
      if(DEBUG>=2 && !(Yfrag[F] > 0.0 && Xfrag[F] > 0.0)){
	printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:F=%d:y=%0.3f,x=%0.3f,misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	       T,U,I,J,G,H,F,Y[I]-Y[G],X[J]-X[H],misscnt,misscnt+(I-G-1),falsecnt,falsecnt + (J-H-1),OutlierCnt,OutlierMN);
	fflush(stdout);
	assert(Yfrag[F] > 0.0);
	assert(Xfrag[F] > 0.0);
      }
      F++;
      misscnt += I-G-1;
      falsecnt += J-H-1;
    } else {/* left end */
      int Lijy = p->Lij1;
      int Lijx = p->Lij2;
      if(DEBUG>=2 && !(Lijy >= 0 && Lijy <= I)){
        #pragma omp critical
        {
	  printf("alignFP:Yid=%d,Xid=%d,or=%d:N=%d,M=%d,score=%0.4f:U=%d,I=%d,J=%d,Lijy=%d,Lijx=%d,misscnt=%d,falsecnt=%d\n",
		 Yid,Xid,orientation,N,M,lscore,U,I,J,Lijy,Lijx,misscnt,falsecnt);
	  fflush(stdout);
	  assert(Lijy >= 0 && Lijy <= I);
	}
      }
      if(DEBUG>=2) 	assert(Lijx >= 0 && Lijx <= J);
      if(LVERB>=2)
	printf("T=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:misscnt=%d->%d,falsecnt=%d->%d\n",T,U,I,J,Lijy,Lijx,misscnt,misscnt+(I-max(1,Lijy)),falsecnt,falsecnt + (J-max(1,Lijx)));

      misscnt += I-max(1,Lijy);
      falsecnt += J-max(1,Lijx);
      if(DEBUG>=2 && !(misscnt>=0)){
        #pragma omp critical
        {
	  printf("alignFP:Yid=%d,Xid=%d,or=%d:N=%d,M=%d,score=%0.4f:I=%d,Lijy=%d,misscnt=%d\n",
		 Yid,Xid,orientation,N,M,lscore,I,Lijy,misscnt);
	  fflush(stdout);
	  assert(misscnt >= 0);
	}
      }
    }
  }
  if(DEBUG) assert((Poutlier <= 0.0) ? F== U-1 : F <= U-1);

  /* right end */
  int Rijy = p->Rij1;
  int Rijx = p->Rij2;

  if(LVERB>=2)
    printf("T=%d/%d:I=%d,J=%d:F=%d,N=%d,M=%d,Rijy=%d,Rijx=%d:misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	   U,U,I,J,F,N,M,Rijy,Rijx,misscnt,misscnt+min(N,Rijy)-I,falsecnt,falsecnt + min(M,Rijx)-J,OutlierCnt,OutlierMN);
  if(DEBUG>=2 && !(min(N,Rijy)-I >= 0)){
    #pragma omp critical
    {
      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d:T=%d/%d:I=%d,J=%d:F=%d,N=%d,M=%d,Rijy=%d,Rijx=%d:misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,U,U,I,J,F,N,M,Rijy,Rijx,misscnt,misscnt+min(N,Rijy)-I,falsecnt,falsecnt + min(M,Rijx)-J,OutlierCnt,OutlierMN);
      for(int T = 0; T < U; T++)
	printf("   T=%d/%d:I=%d,J=%d\n",T,U,p->sites1[T],p->sites2[T]);
      fflush(stdout);
      assert(min(N,Rijy)-I >= 0);
    }
  }
  if(DEBUG>=2) assert(min(M,Rijx)-J >= 0);
  misscnt += min(N,Rijy)-I;
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

  if(DEBUG>=2 && !(misscnt>=0)){
    #pragma omp critical
    {
      printf("alignFP:Yid=%d,Xid=%d,or=%d:N=%d,M=%d,score=%0.4f:I=%d,J=%d,Rijy=%d,misscnt=%d\n",
	     Yid,Xid,orientation,N,M,lscore,I,J,Rijy,misscnt);
      fflush(stdout);
      assert(misscnt >= 0);
    }
  }
  double logFP = pvalue(Xfrag,Yfrag, F, misscnt,falsecnt,OutlierCnt,logOutlierMN, EndOutlierNL, Ylambda[c], res[c]*PixelLen, verb, Yid, Xid,orientation);
  if(DEBUG && !(isfinite(logFP))){
    (void)pvalue(Xfrag,Yfrag,F, misscnt, falsecnt, OutlierCnt, logOutlierMN, EndOutlierNL,Ylambda[c],res[c]*PixelLen,1,Yid,Xid,orientation);
    fprintf(stderr,"alignFP:I=%d,J=%d:logFP=%e,misscnt=%d(Rijy=%d),falsecnt=%d(Rijx=%d),Outliers=%d,MN=%0.1f\n",
	    I,J,logFP,misscnt,Rijy,falsecnt,Rijx,OutlierCnt,OutlierMN);
    fflush(stderr);
    assert(isfinite(logFP));
  }
  if(BackgroundSNR){
    printf("BackgroundSNR not supported during alignment merger with -pairsplit\n");
    exit(1);
  }
  return logFP;
}
#endif


/** extract alignment with right end at Y[I],X[J] from DP array A.field[][] and return pvalue logFP() */
static double alignFP(int c, AL *A, int bestI, int bestJ, int bestR,
		      FLOAT *Y, FLOAT *X, int N, int M,
		      int orientation, int Yid, int Xid, PFLOAT **delY, PFLOAT **delX, double lscore, int verb, int Ascorevalid)
{
  /* allocate space to store alignment in reverse order */
  int *Ilist = (int *)alloca(M*sizeof(int));
  int *Jlist = (int *)alloca(M*sizeof(int));
  int *outlier = (int *)alloca(M*sizeof(int));
  int K = 0;

  int I = bestI;
  int J = bestJ;

  if(outlierExtend){
    while(I >= 0){
      if(DEBUG){
	assert(I >= 1 && I <= N);
	assert(J >= 1 && J <= M);
      }

      Ilist[K] = I;
      Jlist[K] = J;
      int G = FLAT(A, G, I, J);
      if(G > 0){
	int H = FLAT(A,H, I, J);
	if(DEBUG>=2) assert(H > 0);
	if(DEBUG>=2) assert(K+1 < M);
	PFLOAT y = Y[I] - Y[G];
	PFLOAT x = X[J] - X[H];
	if(Poutlier > 0.0){
	  PFLOAT Bias,Pen;
	  SintDetail(x,y, J-H,I-G,J,I,Bias,Pen,c);
	  PFLOAT OutPen = OutlierPenalty;
	  if(outlierBC)
	    OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-G];
	  PFLOAT iscore = Bias + OutlierBias + (OUTLIER_DELTA(x-y) ? max(Pen,OutPen) : Pen);
	  if(DEBUG>=2 && Ascorevalid){
	    PFLOAT iscore2 = FLAT(A,score, I, J) - FLAT(A,score, G, H);
	    if(!(fabs(iscore2 - iscore) < 0.001)){
	      #pragma omp critical
	      {
		printf("alignFP:Yid=%d,Xid=%d,or=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,Ascorevalid=%d\n",
		       Yid,Xid,orientation,bestI,bestJ,bestR,N,M,Ascorevalid);
		printf("K=%d:I=%d,J=%d,G=%d,H=%d,y=%0.4f,x=%0.4f,Bias=%0.6f,Pen=%0.6f,iscore=%0.6f,iscore2=%0.6f (A.score[I][J]=%0.6f,A.score[G][H]=%0.6f), OUTLIER_DELTA(x-y)=%d\n",
		       K,I,J,G,H,y,x,Bias,Pen,iscore,iscore2,FLAT(A,score,I,J),FLAT(A,score,G,H), OUTLIER_DELTA(x-y));
		fflush(stdout);
		assert(fabs(iscore2 - iscore) < 0.001);
	      }
	    }
	  }
	  PFLOAT outscore = Bias + Pen;
	  outlier[K+1] = (iscore > outscore + Poutlier + (PFLOAT)0.01 || (outlierExtend && (J <= H || I <= G))) ? 1 : 0;
	  if(VERB>=2 && verb)
	    printf("K=%d:x=%0.3f,y=%0.3f,iscore=%0.6f,outscore=%0.6f:outlier[K+1]=%d\n",K,x,y,iscore,outscore,outlier[K+1]);
	} else
	  outlier[K+1] = 0;
	J = H;
      }
      K++;
      I = G;
    }
  } else {/* General Case : outlierExtend==0 */
    while(I >= 0){
      if(DEBUG){
	assert(I >= 1 && I <= N);
	assert(J >= 1 && J <= M);
      }

      Ilist[K] = I;
      Jlist[K] = J;
      int G = FLAT(A,G, I, J);
      if(G > 0){
	int H = FLAT(A,H, I, J);
	if(DEBUG>=2) assert(H > 0);
	if(DEBUG>=2) assert(K+1 < M);
	PFLOAT y = delY[I][G-I+DELTA];//Y[I] - Y[G];
	if(DEBUG>=2) assert(fabs(Y[I] - Y[G] - y) < 0.001);
	PFLOAT x = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];//X[J] - X[H];
	if(DEBUG>=2) assert(fabs(X[J] - X[H] - x) < 0.001);
	if(Poutlier > 0.0){
	  PFLOAT Bias,Pen;
	  SintDetail(x,y, J-H,I-G,J,I,Bias,Pen,c);
	  PFLOAT OutPen = OutlierPenalty;
	  if(outlierBC)
	    OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-G];
	  PFLOAT iscore = Bias + OutlierBias + (OUTLIER_DELTA(x-y) ? max(Pen,OutPen) : Pen);
	  if(DEBUG>=2 && Ascorevalid){
	    PFLOAT iscore2 = FLAT(A,score, I, J) - FLAT(A,score, G, H);
	    assert(fabs(iscore2 - iscore) < 0.001);
	  }
	  PFLOAT outscore = Bias + Pen;
	  outlier[K+1] = (iscore > outscore + Poutlier + (PFLOAT)0.01) ? 1 : 0;
	  if(VERB>=2 && verb)
	    printf("K=%d:x=%0.3f,y=%0.3f,iscore=%0.6f,outscore=%0.6f:outlier[K+1]=%d\n",K,x,y,iscore,outscore,outlier[K+1]);
	} else
	  outlier[K+1] = 0;
	J = H;
      }
      K++;
      I = G;
    }
  }/* General Case : outlierExtend==0 */

  int EndOutlierCnt = (bestR <= -2 ? 1 : 0) + (I <= -2 ? 1 : 0);

  if(DEBUG) assert(K>0);
  if(DEBUG) assert(K <= M);
    
  /* allocate alignment fragment sizes */
  FLOAT *Xfrag = (FLOAT *) alloca((K-1)*sizeof(FLOAT));
  FLOAT *Yfrag = (FLOAT *) alloca((K-1)*sizeof(FLOAT));
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
	  if(DEBUG>=2 && !outlierExtend && !(J > H && I > G)){
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
      if(DEBUG>=2 && !(Yfrag[F] > 0.0 && Xfrag[F] > 0.0)){
	PFLOAT Bias,Pen;
	SintDetail(Xfrag[F],Yfrag[F], J-H, I-G, J, I, Bias, Pen,c);
	PFLOAT OutPen = OutlierPenalty;
	if(outlierBC)
	  OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-G];
	PFLOAT iscore = Bias + OutlierBias + (OUTLIER_DELTA(Xfrag[F]-Yfrag[F]) ? max(Pen,OutPen) : Pen);
	PFLOAT outscore = Bias + Pen;

	printf("T=%d/%d:I=%d,J=%d:G=%d,H=%d:F=%d:y=%0.3f,x=%0.3f,misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	       T,K,I,J,G,H,F,Y[I]-Y[G],X[J]-X[H],misscnt,misscnt+(I-G-1),falsecnt,falsecnt + (J-H-1),OutlierCnt,OutlierMN);
	printf("  outlier[K-T]=%d:iscore=%0.6f,outscore=%0.6f,Bias=%0.6f,Pen=%0.6f\n",outlier[K-T],iscore,outscore,Bias,Pen);
	fflush(stdout);
	assert(Yfrag[F] > 0.0);
	assert(Xfrag[F] > 0.0);
      }
      F++;
      misscnt += I-G-1;
      falsecnt += J-H-1;
    } else {/* left end */
      if(DEBUG>=2 && !(Ch > 0.0 || PoutlierEnd > 0.0 || (RepeatRec && RepeatMaxShift > 0 && RepeatPvalueRatio > 0.0)) && !(A->G[I][J] > -2)){
	#pragma omp critical
	{
	  printf("alignFP:Yid=%d,Xid=%d,or=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,score=%0.4f:U=%d,I=%d,J=%d,Ch=%0.6e,PoutlierEnd=%0.6e,A.G[I][J]=%d,A.score[I][J]=%0.6f\n",
		 Yid,Xid,orientation,bestI,bestJ,bestR,N,M,lscore,K,I,J,Ch,PoutlierEnd,A->G[I][J],A->score[I][J]);
	  fflush(stdout);
	  assert(A->G[I][J] > -2);
	}
      }
      int Lijy = (FLAT(A,G, I, J) <= -2) ? I : FLAT(A,Lijy, I, J);
      int Lijx = (FLAT(A,G, I, J) <= -2) ? J : FLAT(A,Lijx, I, J);
      if(DEBUG>=2 && !(Lijy >= 0 && Lijy <= I)){
        #pragma omp critical
        {
	  printf("alignFP:Yid=%d,Xid=%d,or=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,score=%0.4f:U=%d,I=%d,J=%d,Lijy=%d,Lijx=%d,misscnt=%d,falsecnt=%d\n",
		 Yid,Xid,orientation,bestI,bestJ,bestR,N,M,lscore,K,I,J,Lijy,Lijx,misscnt,falsecnt);
	  fflush(stdout);
	  assert(Lijy >= 0 && Lijy <= I);
	}
      }
      if(DEBUG>=2) 	assert(Lijx >= 0 && Lijx <= J);
      if(LVERB>=2)
	printf("T=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:misscnt=%d->%d,falsecnt=%d->%d\n",T,K,I,J,Lijy,Lijx,misscnt,misscnt+(I-max(1,Lijy)),falsecnt,falsecnt + (J-max(1,Lijx)));

      misscnt += I-max(1,Lijy);
      falsecnt += J-max(1,Lijx);
      if(DEBUG>=2 && !(misscnt>=0)){
        #pragma omp critical
        {
	  printf("alignFP:Yid=%d,Xid=%d,or=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,score=%0.4f:I=%d,Lijy=%d,misscnt=%d\n",
		 Yid,Xid,orientation,bestI,bestJ,bestR,N,M,lscore,I,Lijy,misscnt);
	  fflush(stdout);
	  assert(misscnt >= 0);
	}
      }
    }
  }
  if(DEBUG) assert(I==Ilist[0] && I==bestI);
  if(DEBUG) assert(J==Jlist[0] && J==bestJ);
  if(DEBUG) assert((Poutlier <= 0.0) ? F== K-1 : F <= K-1);

  /* right end */
  if(DEBUG>=2 && !(Ch > 0.0 || PoutlierEnd > 0.0)) assert(bestR > -2);
  int Rijy = (bestR <= -2) ? I : FLAT(A,Rijy, I, J);
  int Rijx = (bestR <= -2) ? J : FLAT(A,Rijx, I, J);
  if(LVERB>=2)
    printf("T=%d/%d:I=%d,J=%d:F=%d,N=%d,M=%d,Rijy=%d,Rijx=%d:misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	   K,K,I,J,F,N,M,Rijy,Rijx,misscnt,misscnt+min(N,Rijy)-I,falsecnt,falsecnt + min(M,Rijx)-J,OutlierCnt,OutlierMN);
  if(DEBUG>=2 && !(min(N,Rijy)-I >= 0)){
    #pragma omp critical
    {
      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d:T=%d/%d:I=%d,J=%d:F=%d,N=%d,M=%d,Rijy=%d,Rijx=%d:misscnt=%d->%d,falsecnt=%d->%d,Outliers=%d,MN=%0.1e\n",
	     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,K,K,I,J,F,N,M,Rijy,Rijx,misscnt,misscnt+min(N,Rijy)-I,falsecnt,falsecnt + min(M,Rijx)-J,OutlierCnt,OutlierMN);
      printf("bestI=%d,bestJ=%d\n",bestI,bestJ);
      for(int T = 0; T < K; T++)
	printf("   T=%d/%d:I=%d,J=%d\n",T,K,Ilist[K-1-T],Jlist[K-1-T]);
      fflush(stdout);
      assert(min(N,Rijy)-I >= 0);
    }
  }
  if(DEBUG>=2) assert(min(M,Rijx)-J >= 0);
  misscnt += min(N,Rijy)-I;
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

  if(DEBUG>=2 && !(misscnt>=0)){
    #pragma omp critical
    {
      printf("alignFP:Yid=%d,Xid=%d,or=%d:bestI=%d,bestJ=%d,bestR=%d,N=%d,M=%d,score=%0.4f:I=%d,J=%d,Rijy=%d,misscnt=%d\n",
	     Yid,Xid,orientation,bestI,bestJ,bestR,N,M,lscore,I,J,Rijy,misscnt);
      fflush(stdout);
      assert(misscnt >= 0);
    }
  }
  double logFP = pvalue(Xfrag,Yfrag, F, misscnt,falsecnt,OutlierCnt,logOutlierMN,EndOutlierNL,Ylambda[c], res[c]*PixelLen,verb, Yid, Xid,orientation);
  if(DEBUG && !(isfinite(logFP))){
    (void)pvalue(Xfrag,Yfrag,F, misscnt,falsecnt,OutlierCnt,logOutlierMN,EndOutlierNL,Ylambda[c],res[c]*PixelLen,1,Yid,Xid,orientation);
    fprintf(stderr,"alignFP:I=%d,J=%d:logFP=%e,misscnt=%d(Rijy=%d),falsecnt=%d(Rijx=%d),Outliers=%d,MN=%0.1f\n",
	    I,J,logFP,misscnt,Rijy,falsecnt,Rijx,OutlierCnt,OutlierMN);
    fflush(stderr);
    assert(isfinite(logFP));
  }
  if(BackgroundSNR){/* compute pvalue for site matches based on Kalmagorov-Smirnov distributions */
    double HalfByLog10 = 0.5/log(10.0);
    if(!YYmap[Yid]->SNRcnt[0]){
      printf("Missing SNR values for id=%lld (needed for SNR based pvalues)\n",YYmap[Yid]->id);
      exit(1);
    }
    if(!XXmap[Xid]->SNRcnt[0]){
      printf("Missing SNR values for id=%lld (needed for SNR based pvalues)\n",XXmap[Xid]->id);
      exit(1);
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
	printf("Yid=%d,Xid=%d,or=%d:SNR match for I=%d,J=%d(JR=%d):n1=%d(gmSNR=%0.2f),n2=%d(gmSNR=%0.2f),nBG=%d: LogPv12=%0.6f,LogPvBg1=%0.6f,LogPvBg2=%0.6f: -LogPV += %0.4f\n",
	       Yid,Xid,orientation,I,J,JR,YYmap[Yid]->SNRcnt[0][I],YYmap[Yid]->SNRgmean[0][I],XXmap[Xid]->SNRcnt[0][JR], XXmap[Xid]->SNRgmean[0][JR], bgSNRcnt[0], Pv12/log(10.0),PvBg1/log(10.0),PvBg2/log(10.0), -HalfByLog10*(PvBg1+PvBg2 - 2.0*Pv12));

      logFP -= HalfByLog10*(PvBg1 + PvBg2 - 2.0*Pv12);
    }
    if(LVERB || verb)
      fflush(stdout);
  }

  return logFP;
}

/** The main pairwise Dynamic Programming routine. 
Hardware accellerations done with USE_MIC and USE_AVX preprocessor directives.
*/
void pairalignXY(FLOAT **YY, int *NN,
		 FLOAT **XX, int *MM,
		 AL *AA,
		 Calign *align,/* Best alignment (with score) is stored at align[-1] and align[0 .. colors-1] */
		 int Yid,int Xid,
		 PFLOAT ***delYY,/* delYY[c][I=2..NN[c]][n=max(0,DELTA+1-I)..DELTA-1] == Y[I] - Y[I-DELTA+n] */
		 PFLOAT ***delXX,/* delXX[c][J=2..MM[c]][m=max(0,DELTA+1-J)..DELTA-1] == X[J] - X[J-DELTA+m] 
				  Except If (USE_SSE && USE_AVX || USE_MIC) && USE_PFLOAT : delXX[c][n=1..DELTA][J=1+n..M] = X[J] - X[J-n] */
		 Chashpair *phash,/* If != 0 pointer to hashtable entry : Not currently used */
		 int orientation
		 )
{
  PFLOAT C = ChimScore;
  int localtype = (PoutlierEnd > 0.0) ? -3 : -2;

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

  if(VERB >= (EVERB?1:2)){
    if(phash)
      printf("pairalignXY:Yid=%d(%lld),Xid=%d(%lld),orientation=%d:phash:id1=%d,id2=%d,or=%d,score=%d,offset=%d\n",
	     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,phash->id1,phash->id2,phash->orientation,phash->hashscore,phash->offset);
    else
      printf("pairalignXY:Yid=%d(%lld),Xid=%d(%lld),orientation=%d:XleftClipped=%d,XrightClipped=%d,YleftClipped=%d,YrightClipped=%d,C=%0.6f\n",
	     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,XleftClipped,XrightClipped,YleftClipped,YrightClipped,C);
    fflush(stdout);
  }
  if(phash && hashdelta){
    #pragma omp critical
    {
      printf("-hashdelta Not implemented for pairwise alignment\n");
      exit(1);
    }
  }

  PFLOAT *Yfc[2], *Xfc[2], *Yrc[2],*Xrc[2];

  /* first align each color seperately (without performing back trace) */
  for(int c = 0; c < colors; c++){
    FLOAT *Y = YY[c];
    FLOAT *X = XX[c];
    int N = NN[c];
    int M = MM[c];
    AL *A = &AA[c];
    PFLOAT **delY = delYY[c];
    PFLOAT **delX = delXX[c];

    PFLOAT *Yf = Yfc[c] = (PFLOAT*)alloca(sizeof(PFLOAT)*(N+1));/* precompute Yf[I] = (PFLOAT)Y[I] */
    PFLOAT *Xf = Xfc[c] = (PFLOAT*)alloca(sizeof(PFLOAT)*(M+1));/* precompute Xf[J] = (PFLOAT)[J] */
    PFLOAT *Yr = Yrc[c] = (PFLOAT*)alloca(sizeof(PFLOAT)*(N+1));/* precompute Yr[I] = (PFLOAT)(Y[N+1]-Y[I]) */
    PFLOAT *Xr = Xrc[c] = (PFLOAT*)alloca(sizeof(PFLOAT)*(M+1));/* precompute Xr[J] = (PFLOAT)(X[M+1]-X[J]) */
    for(int I = 1; I <= N; I++){
      Yf[I] = Y[I];
      Yr[I] = Y[N+1] - Y[I];
    }
    for(int J = 1; J <= M; J++){
      Xf[J] = X[J];
      Xr[J] = X[M+1] - X[J];
    }

    /* initialize A.score[I=1..N][J=1..M] = C (typically MINSCORE) and G = localtype */
    if(RepeatShift > 0.0){
      for(int I = 1; I <= N; I++){
	int J;
	for(J = 1; J < I; J++){
	  if(Y[I] - X[J] < RepeatShift)
	    break;
	  FLAT(A,score, I, J) = C;
	}
	for(; J <= M; J++)
	  FLAT(A,score, I, J) = MINSCORE;
      }
    } else {
      for(int J=1; J <= M; J++)
	FLAT(A,score, 1, J) = C;

      for(int I = 2; I <= N; I++)
	memcpy(&FLAT(A,score, I, 1), &FLAT(A,score, 1, 1), M*sizeof(PFLOAT));
    }

    if(!DEFER_BP || outlierExtend){
      for(int J = 1; J <= M; J++)
	FLAT(A,G, 1, J) = localtype;
      for(int I = 2; I <= N; I++)
	FLAT(A,G, I, 1) = localtype;
#ifdef VALGRIND
      for(int J = 1; J <= M; J++)
	FLAT(A,H, 1, J) = localtype;
      for(int I = 2; I <= N; I++)
	FLAT(A,H, I, 1) = localtype;
#endif
    }

    /* initialize A->Lijy[I=1..N][J=1..M] as leftmost site in Y that overlaps X with Y[I] aligned with X[J] */
    /* only cases where I <= DELTA OR J <= DELTA are initialized */

    int J;
    for(J=1; J <= M; J++){
      int Imax = !ENDINIT ? N : (J <= DELTA) ? N : min(DELTA,N);
    
      int Lijy = 0;
      for(int I=1; I <= Imax; I++){ /* update Lij by incrementing it until Y[Lijy] overlaps X */
	FLOAT delta = Y[I] - X[J];
	while(delta > Y[Lijy])
	  Lijy++;
	if(DEBUG>=2)	assert(Lijy <= I);
	FLAT(A,Lijy, I, J) = Lijy;
      }
      if(!Lijy)
	break;/* all remaining Lijy values will be 0 */
    }
    int LijyJmax = min(M,J);
    int DELTA_M = min(M,DELTA);// DELTA_J = min(M,J);
    for(J++; J <= DELTA_M; J++)
      for(int I=1; I <= N; I++)
	FLAT(A,Lijy, I, J) = 0;		
    if(J <= M){
      int Imax = ENDINIT ? N : min(DELTA,N);
      for(int I=1; I <= Imax; I++)
	memset_int(&FLAT(A,Lijy, I, J), 0, M-J+1);		
    }

    /* initialize A->Lijx.[I=1..N][J=1..M] as leftmost site in X that overlap Y with Y[I] aligned with X[J] */
    /* only cases where I <= DELTA OR J <= DELTA are initialized */
    int I;
    for(I=1; I <= N; I++){
      int Jmax = !ENDINIT ? M : (I <= DELTA) ? M : min(DELTA,M);
      int Lijx = 0;
      for(int J=1; J <= Jmax; J++){ /* update Lijx by incrementing it until X[Lijx] overlaps Y */
	FLOAT delta = X[J]-Y[I];
	while(delta > X[Lijx])
	  Lijx++;
	if(DEBUG>=2)	assert(Lijx <= J);
	FLAT(A,Lijx, I, J) = Lijx;
      }
      if(!Lijx)
	break;/* all remaining Lijx values will be 0 */
    }

    int LijxImax = min(N,I);/* Lijx == 0 for I > LijxImax  */
    int DELTA_N = min(N,DELTA); // DELTA_I = min(N,I);
    for(I++; I <= DELTA_N; I++)
      memset_int(&FLAT(A,Lijx, I, 1), 0, M);
    if(I <= N){
      int Jmax = !ENDINIT ? M : min(DELTA,M);
      for(; I <= N; I++)
	memset_int(&FLAT(A,Lijx, I, 1), 0, Jmax);
    }

    if(DEBUG>=2){
      for(int I = 1; I <= N && I <= DELTA; I++)
	for(int J = 1; J <= M && J <= DELTA; J++){
	  int Lijx = FLAT(A,Lijx, I, J);
	  int Lijy = FLAT(A,Lijy, I, J);
	  assert(Lijx <= 0 || Lijy <= 0);
	}
    }

#if USE_MIC
    float *scoreSend = (float *)alloca((M+1)*sizeof(scoreSend));
    float *scoreSend1 = &scoreSend[-1];
#endif
  
    /* initialize A->score[I=1..N][J=1..M] with left ends unaligned and Y[I] aligned with X[J] */
    /* only cases where I <= DELTA OR J <= DELTA are initialized */
    for(I=1; I <= LijxImax; I++){
      int Jmax = (I <= DELTA) ? M : min(DELTA,M);
      if(RepeatShift > 0.0){
	Jmax = min(I-1,Jmax);
	while(Jmax >= 1 && Y[I]-X[Jmax] < RepeatShift)
	  Jmax--;
      }

      int *LijxI = &FLAT(A,Lijx, I, 0), *LijyI = &FLAT(A,Lijy, I, 0);
      PFLOAT *scoreI = &FLAT(A,score, I, 0);

#if USE_MIC
      int score_size = Jmax;
      __mmask16 mask;
      __m512 vIdx = _mm512_add_ps(_mm512_set1_ps(I+2), v512_ascending0123);

      //prefetch_array(&Xf[1], score_size*sizeof(*Xf), _MM_HINT_T0);
    
      __m512 v_YfI = _mm512_set1_ps(Yf[I]); 
      __m512 v1 = _mm512_set1_ps(1.0f);
      __m512i v_limit = _mm512_sub_epi32(_mm512_set1_epi32(Jmax), v512i_ascending0123);
    
      for(int J= 1; J <= Jmax; J += 16){
	/* make first pass of the loops to speedup computation */
	//      mask = (J < Jmax - 15) ? 0xffff : (0xffff>>(J+15-Jmax));
      mask = _mm512_cmp_epi32_mask(_mm512_set1_epi32(J), v_limit, _MM_CMPINT_LE);
      __m512i LijxIJ = _mm512_maskloadu_epi32(&LijxI[J], mask);
      __m512i LijyIJ = _mm512_maskloadu_epi32(&LijyI[J], mask);
      __m512 LijxIJd = _mm512_cvtfxpnt_round_adjustepi32_ps(LijxIJ, _MM_FROUND_CUR_DIRECTION, _MM_EXPADJ_NONE);
      __m512 LijyIJd = _mm512_cvtfxpnt_round_adjustepi32_ps(LijyIJ, _MM_FROUND_CUR_DIRECTION, _MM_EXPADJ_NONE);
      __m512 v_score = Send512(
			   _mm512_gmin_ps(
					  _mm512_maskloadu_ps(&Xf[J], mask),
					  v_YfI
					  ),
			   _mm512_sub_ps(
					 _mm512_add_ps(
						       vIdx, 
						       _mm512_set1_ps(J)),
					 _mm512_add_ps(
						       _mm512_gmax_ps(v1, LijxIJd),
						       _mm512_gmax_ps(v1, LijyIJd))
					 ),
			   c);
      _mm512_mask_store_ps(&scoreSend1[J], mask, v_score);
      }
#endif // USE_MIC
      
      int DELTA_Jmax = min(Jmax,LijyJmax);

      int J;
      for(J=1; J <= DELTA_Jmax; J++){
	int Lijx = LijxI[J];
	int Lijy = LijyI[J];
	PFLOAT ascore;
	if(ENDFIX>=2){
	  PFLOAT SmIJ = Sm(J,I);
	  if(Lijx <= 0){/* left end of X fully overlapped by Y */
	    if(!XleftClipped){
#if USE_MIC
	      ascore = scoreSend1[J] + SmIJ;
#else
	      ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I+1-max(1,Lijy)),c) + SmIJ;
#endif
	      ascore = max(ascore,C);
	      for(int K= Lijy; K < I; K++) {
		PFLOAT bound = Sbnd(Xf[J],Y[I]-Y[K],(PFLOAT)(J+I-K),c) + SmIJ;
		if(bound <= ascore)
		  break;// HERE : see ENDFIX2 in refalign.cpp
		ascore = bound;
	      }
	    } else
	      ascore = scoreI[J];
	  } else {/* left end of Y fully overlapped by X */
	    if(DEBUG>=2) assert(Lijy <= 0);
	    if(!YleftClipped){
#if USE_MIC
	      ascore = scoreSend1[J] + SmIJ;
#else
	      ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-max(1,Lijx)+I),c) + SmIJ;
#endif
	      ascore = max(ascore,C);
	      for(int K= Lijx; K < J; K++){
		PFLOAT bound = Sbnd(X[J]-X[K],Yf[I],(PFLOAT)(J-K+I),c) + SmIJ;
		if(bound <= ascore)
		  break;// HERE : see ENDFIX2 in refalign.cpp
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
	    ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-max(1,Lijx)+I+1-max(1,Lijy)),c) + Sm(J,I);
#endif
	    ascore = max(ascore,C);
	  } else
	    ascore = scoreI[J];
	}
	if(DEBUG>=2) assert(ascore >= C);
	scoreI[J] = ascore;
      }

      for(; J <= Jmax; J++){
	int Lijx = LijxI[J];
	if(Lijx > 0)
	  break;
	/* left end of X fully overlapped by Y */
	PFLOAT ascore;
	if(!XleftClipped){
	  PFLOAT SmIJ = Sm(J,I);
#if USE_MIC
	  ascore = scoreSend1[J] + SmIJ;
#else
	  ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I),c) + SmIJ;
#endif
	  ascore = max(ascore,C);
	  if(ENDFIX>=2){
	    for(int K= Lijx; K < J; K++){
	      PFLOAT bound = Sbnd(X[J]-X[K],Yf[I],(PFLOAT)(J+I-K),c) + SmIJ;
	      if(bound <= ascore)
		break;// HERE : see ENDFIX2 in refalign.cpp
	      ascore = bound;
	    }
	  }
	} else
	  ascore = scoreI[J];

	if(DEBUG>=2) assert(ascore >= C);
	scoreI[J] = ascore;
      }
      for(; J <= Jmax; J++){
	int Lijx = LijxI[J];
	if(DEBUG>=2) assert(Lijx > 0);
	/* left end of Y fully overlapped by X */
	PFLOAT ascore;
	PFLOAT SmIJ = Sm(J,I);
	if(!YleftClipped){
#if USE_MIC
	  ascore = scoreSend1[J] + SmIJ;
#else
	  ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-Lijx+I),c) + SmIJ;
#endif
	  ascore = max(ascore,C);
	  if(ENDFIX>=2){
	    for(int K= Lijx; K < J; K++){
	      PFLOAT bound = Sbnd(X[J]-X[K],Yf[I],(PFLOAT)(J+I-K),c) + SmIJ;
	      if(bound <= ascore)
		break;// HERE : see ENDFIX2 in refalign.cpp
	      ascore = bound;
	    }
	  } 
	} else
	  ascore = scoreI[J];
	if(DEBUG>=2) assert(ascore >= C);
	scoreI[J] = ascore;
      }
    }

    if(DEBUG && I <= N) assert(A->Lijy[I][(I<=DELTA) ? M : min(DELTA,M)] > 0);

    for(; I <= N; I++){
      int Jmax = (I <= DELTA) ? M : min(DELTA,M);
      if(RepeatShift > 0.0){
	Jmax = min(I-1,Jmax);
	while(Jmax >= 1 && Y[I]-X[Jmax] < RepeatShift)
	  Jmax--;
      }

      int *LijyI = &FLAT(A,Lijy, I, 0);
      PFLOAT *scoreI = &FLAT(A,score, I, 0);

#if USE_MIC
      int score_size = Jmax;
      __mmask16 mask;
      __m512 vIdx = _mm512_add_ps(_mm512_set1_ps(I+1), v512_ascending0123);

      //prefetch_array(&Xf[1], score_size*sizeof(*Xf), _MM_HINT_T0);
      __m512 v_YfI = _mm512_set1_ps(Yf[I]);
      __m512 v1 = _mm512_set1_ps(1.0f);
      __m512i v_limit = _mm512_sub_epi32(_mm512_set1_epi32(Jmax), v512i_ascending0123);
    
      for(int J= 1; J <= Jmax; J+=16){
	/* make first pass of the loops to speedup computation */
	//mask= J<Jmax-7 ? 0xff : (0xFF>>(J+7-Jmax));
	mask = _mm512_cmp_epi32_mask(_mm512_set1_epi32(J), v_limit, _MM_CMPINT_LE);
	__m512i LijyIJ = _mm512_maskloadu_epi32(&LijyI[J], mask);
	__m512 LijyIJd = _mm512_cvtfxpnt_round_adjustepi32_ps(LijyIJ, _MM_FROUND_CUR_DIRECTION, _MM_EXPADJ_NONE);
	__m512 v_score = Send512(
				 _mm512_gmin_ps(
						_mm512_maskloadu_ps(&Xf[J], mask),
						v_YfI
						),
				 _mm512_sub_ps(
					       _mm512_add_ps(
							     vIdx, 
							     _mm512_set1_ps(J)),
					       _mm512_gmax_ps(v1, LijyIJd)
					       ),
				 c);
	_mm512_mask_store_ps(&scoreSend1[J], mask, v_score);
      }
#endif
      
      for(int J=1; J <= Jmax; J++){
	int Lijy = LijyI[J];
	if(DEBUG>=2) assert(Lijy > 0);
	/* left end of X fully overlapped by Y */
	PFLOAT ascore;
	if(!XleftClipped){
	  PFLOAT SmIJ = Sm(J,I);
#if USE_MIC
	  ascore = scoreSend1[J] + SmIJ;
#else
	  ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I+1-Lijy),c) + SmIJ;
#endif
	  ascore = max(ascore,C);
	  if(ENDFIX>=2){
	    for(int K= Lijy; K < I; K++) {
	      PFLOAT bound = Sbnd(Xf[J],Y[I]-Y[K],(PFLOAT)(J+I-K),c) + SmIJ;
	      if(bound <= ascore)
		break;// HERE : see ENDFIX2 in refalign.cpp
	      ascore = bound;
	    }
	  }
	} else
	  ascore = scoreI[J];
	if(DEBUG>=2) assert(ascore >= C);
	scoreI[J] = ascore;
      }
    }

    /* initialize A[I=1..N][J= 1 to M].Rijy as rightmost site in Y that overlaps X when
       Y[I] is aligned with X[J] */
    /* If ENDINIT : only cases where I >= N+1-DELTA OR J >= M+1-DELTA are initialized */
    int Rijy = -1;
    for(int J= M; J >= 1; J--){
      int Imin = !ENDINIT ? 1 : (J >= M+1-DELTA) ? 1 : max(1,N+1-DELTA);
    
      if(Rijy==N+1) {
	if(J>=M+1-DELTA) {
	  for(int I= N; I >= Imin; I--)
	    FLAT(A,Rijy, I, J) = Rijy;
	} else {
	  for(int I= N; I >= Imin; I--)
	    memset_int(&FLAT(A,Rijy, I, 1), Rijy, J);
	  break;
	}
      } else {
	Rijy = N+1;
	int *RijyJ = &FLAT(A,Rijy, 0, J);
	int stride = A->M_stride;
	for(int I= N; I >= Imin; I--){
	  /* update Rijy by decrementing it until Y[Rijy] overlaps X */
	  FLOAT delta = Y[I] + X[M+1] - X[J] + 1e-8;
	  while(Y[Rijy] > delta)
	    Rijy--;
	  if(DEBUG>=2)	assert(Rijy >= I);
	  RijyJ[I*stride] = Rijy;
	}
      }
    }

    /* initialize A[I=1..N][J= 1 to M].Rijx as rightmost site in X that overlaps Y when
       Y[I] is aligned with X[J] */
    /* If ENDINIT: only cases where I >= N+1-DELTA OR J >= M+1-DELTA are initialized */
    int Rijx = -1;
    for(int I= N; I >= 1; I--){
      int Jmin = !ENDINIT ? 1 : (I >= N+1-DELTA) ? 1 : max(1,M+1-DELTA);
      if(Rijx==M+1) {
	if(I>=N+1-DELTA) {
	  memset_int(&FLAT(A,Rijx, I, Jmin), Rijx, M-Jmin+1);
	} else {
	  for(;I >= 1; I--)
	    memset_int(&FLAT(A,Rijx, I, Jmin), Rijx, M-Jmin+1);
	  break;
	}
      } else {
	Rijx = M+1;
	int *RijxI = &FLAT(A,Rijx, I, 0);
	for(int J= M; J >= Jmin; J--){
	  /* update Rijx by decrementing it until X[Rijx] overlaps Y */
	  FLOAT delta = X[J] + Y[N+1] - Y[I] + 1e-8;
	  while(X[Rijx] > delta)
	    Rijx--;
	  if(DEBUG>=2)	assert(Rijx >= J);
	  RijxI[J] = Rijx;
	}
      }
    }

    if(DEBUG>=2 && ENDINIT){
      for(int I = max(1,N+1-DELTA); I <= N; I++)
	for(int J = max(1,M+1-DELTA); J <= M; J++){
	  int Rijx = FLAT(A,Rijx, I, J);
	  int Rijy = FLAT(A,Rijy, I, J);
	  if(!(Rijx > M || Rijy > N)){
	    printf("Yid=%d(%lld),Xid=%d(%lld):I=%d,J=%d:Rijx=%d,Rijy=%d,M=%d,N=%d,DELTA=%d\n",
		   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,I,J,Rijx,Rijy,M,N,DELTA);
	    printf("  Y[I]=%0.8f,Y[N+1]=%0.8f,X[J]=%0.8f,X[M+1]=%0.8f\n",Y[I],Y[N+1],X[J],X[M+1]);
	    fflush(stdout);
	    assert(Rijx > M || Rijy > N);
	  }
	}
    }

    if((DEFER_BP && !outlierExtend) || RepeatRec) /* save score as Uscore */
      for(int I=1; I <= N; I++)
	memcpy(&FLAT(A,Uscore, I, 1), &FLAT(A,score, I, 1), M*sizeof(FLAT(A,score, I, 1)));

    /* Dynamic programming recurance */
    for(int I = 2; I <= N; I++){
      //    prefetch_array(&A, sizeof(A), _MM_HINT_T0);
      //    prefetch_array(&FLAT(A,score, I, 0), M*sizeof(FLAT(A,score, I, 1)), _MM_HINT_T0);
	  
      int Jmin = 2, Jmax = M;
      if(RepeatShift > 0.0 && !orientation){/* decrease Jmax to below I until Y[I] - X[Jmax] >= RepeatShift */
	Jmax = I-1;
	while(Jmax >= Jmin && Y[I]-X[Jmax] < RepeatShift)
	  Jmax--;
      }
    
      if(DEFER_BP && !outlierExtend && !outlierBC && DELTA==4 && I > 4){/* critical case, customized to help compiler to unroll the G and H loops */
	// NOTE : Currently only (USE_SSE && USE_AVX && USE_PFLOAT) OR (USE_MIC && USE_PFLOAT) use intrinsics

 #if !(USE_SSE && USE_AVX && USE_PFLOAT) && !(USE_MIC && USE_PFLOAT)

	/* Original path */    
	int J = Jmin;
	for(int jmax = min(4,Jmax); J <= jmax; J++){
	  PFLOAT score = FLAT(A,score, I, J);
	  for(int G= I; --G >= I - 4;){
	    PFLOAT deltaY = delY[I][G-I+4];// Y[I] - Y[G];
	    if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	    for(int H= J; --H > 0;){
	      if(DEBUG>=2) assert(fabs(X[J]-X[H] - delX[J][H-J+4]) < 0.001);
	      PFLOAT newscore = FLAT(A,score, G, H) + Sint(delX[J][H-J+4],deltaY,J-H,I-G,c);
	      score = PMAX(score,newscore);
	    }
	  }
	  FLAT(A,score, I, J) = score;
	}

	for(; J <= Jmax; J++){
	  PFLOAT score = FLAT(A,score, I, J);
	  for(int G= I; --G >= I - 4;){
	    PFLOAT deltaY = delY[I][G-I+4];
	    if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	    for(int H= J; --H >= J - 4;){
	      if(DEBUG>=2) assert(fabs(X[J]-X[H] - delX[J][H-J+4]) < 0.001);
	      PFLOAT newscore = FLAT(A,score, G, H) + Sint(delX[J][H-J+4],deltaY,J-H,I-G,c);
	      score = PMAX(score,newscore);
	    }
	  }	
	  FLAT(A,score, I, J) = score;
	}
#else  // (USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)

	// NOTE : This block only handles the first 4 iterations of J, the others are handled later

  #if (USE_SSE && USE_AVX && USE_PFLOAT) 
	if(DEBUG) assert(USE_SSE && USE_AVX && USE_PFLOAT);

	float fm[4]={4, 3, 2, 1}, tmp[4];
	__m128 vfm=_mm_load_ps(fm);

	int J = Jmin;
	for(int jmax = min(4,Jmax); J <= jmax; J++){
	  __m128 vns = _mm_set_ps1(FLAT(A,score, I, J));
	  for(int H= J; --H > 0;){
	    PFLOAT deltaX = delX[J-H][J];// X[J] - X[H];
	 
	    //	  for(int G= I; --G >= I - 4;){
	    //	    PFLOAT deltaY = delY[I][G-I+4];// Y[I] - Y[G];
	    __m128 vy = _mm_loadu_ps(&delY[I][0]);// X[J] - X[J-4] ... X[J] - X[J-1] : only the last J-1 values of the vector are valid 
	    tmp[0] = A->score[I-4][H];
	    tmp[1] = A->score[I-3][H];
	    tmp[2] = A->score[I-2][H];
	    tmp[3] = A->score[I-1][H];
	    vns = _mm_max_ps(vns, _mm_add_ps(_mm_load_ps(tmp), Sint128(_mm_set_ps1(deltaX), vy, _mm_set_ps1((float)(J-H)), vfm,c)));
	    //	  }
	  }
	  vns=_mm_max_ps(vns, _mm_movehl_ps(vns, vns));

#ifndef __SSE3__
	  vns=_mm_max_ps(vns, _mm_shuffle_ps(vns, vns, _MM_SHUFFLE(1,1,1,1)));
#else
	  vns=_mm_max_ss(vns, _mm_movehdup_ps(vns));
#endif

	  _mm_store_ss(&(FLAT(A,score, I, J)), vns);
	}

#else // (USE_MIC && USE_PFLOAT)
	if(DEBUG) assert(USE_MIC && USE_PFLOAT);

	//Sint512_prefetch();
	int Gmin = I-4, Hmin = 1;
	PFLOAT *delYI = delY[I];

	__m512i Hvi, Gvi, Gv, Hv, IGv, JHv, Iv, Jv, AvIdx;
	__m512 delXv, delYv, Av, score_v;
	__mmask16 mask;
	Hvi = _mm512_set_epi32(0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3);
	Gvi = _mm512_set_epi32(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3);

	Gv = _mm512_add_epi32(_mm512_set1_epi32(Gmin), Gvi);
	Hv = _mm512_add_epi32(_mm512_set1_epi32(Hmin), Hvi);
	AvIdx = _mm512_add_epi32(_mm512_mullo_epi32(_mm512_sub_epi32(Gv, _mm512_set1_epi32(1)), _mm512_set1_epi32(A->M_stride)), Hv);
	//	Av = _mm512_mask_i32gather_ps(v512_PZERO, mask, AvIdx, A->score[1], 4);
	Av= _mm512_i32gather_ps(AvIdx, A->score[1], 4);
	
	int J = Jmin;
	{
	  PFLOAT *AscoreIJmin = &FLAT(A,score, I, Jmin);
	  PFLOAT *delXJmin = &(delX[1][Jmin]);
	  int jmax = min(4,Jmax);
	  IGv = _mm512_sub_epi32(_mm512_set1_epi32(I), Gv);
	  delYv = _mm512_i32gather_ps(_mm512_sub_epi32(_mm512_set1_epi32(DELTA), IGv), delYI, 4); 
	
#define BLOCKJ(Ja) if(Ja <= (jmax - Jmin)){				\
	    int J = Ja + Jmin;						\
	    __m512 score_orig_v = _mm512_set1_ps(AscoreIJmin[Ja]);	\
	    JHv = _mm512_sub_epi32(_mm512_set1_epi32(J), Hv);		\
	    mask = _mm512_cmp_epi32_mask(JHv, _mm512_set1_epi32(0), _MM_CMPINT_GT); \
	    delXv = _mm512_mask_i32gather_ps(v512_PZERO, mask, _mm512_mullo_epi32(_mm512_set1_epi32(COMPUTE_STRIDE4(M-1)), _mm512_sub_epi32(JHv, _mm512_set1_epi32(1))), &delXJmin[Ja], 4); \
	    score_v = _mm512_mask_add_ps(score_orig_v, mask, Av, Sint512b(delXv, delYv, _mm512_cvtfxpnt_round_adjustepu32_ps(_mm512_add_epi32(JHv, IGv), _MM_FROUND_CUR_DIRECTION, _MM_EXPADJ_NONE), c)); \
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
 
	for(; J <= Jmax-7; J+=8){
	  __m256 v_score = _mm256_loadu_ps(&FLAT(A,score, I, J));
	  //	__m256 v_Xj = _mm256_loadu_ps(&X[J]);
	  for(int G= I; --G >= I - 4;){
	    PFLOAT deltaY = delY[I][G-I+4];// Y[I] - Y[G];
	    if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	    __m256 v_deltaY = _mm256_broadcast_ss(&deltaY);
	    float *AscoreGJ = &FLAT(A,score, G, J);
#define BLOCK(H) {							\
	      __m256 v_prevscore = _mm256_loadu_ps(&AscoreGJ[-H]);	\
	      __m256 v_newscore = _mm256_add_ps(v_prevscore, Sint256(_mm256_loadu_ps(&delX[H][J]), v_deltaY, H, I-G,c)); \
	      v_score = _mm256_max_ps(v_score, v_newscore);		\
	    }
	    
	    BLOCK(1);
	    BLOCK(2);
	    BLOCK(3);
	    BLOCK(4);
#undef BLOCK
	  }	
	  _mm256_storeu_ps(&FLAT(A,score, I, J), v_score);
	}

	if(J <= Jmax) {
	  unsigned int mask[14] = {0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0, 0, 0, 0, 0, 0, 0};
	  __m256i mask256 = _mm256_loadu_si256((__m256i*) &(mask[6+J-Jmax]));

	  __m256 v_score = _mm256_maskload_ps(&FLAT(A,score, I, J), mask256);
	  for(int G= I; --G >= I - 4;){
	    PFLOAT deltaY = delY[I][G-I+4];// Y[I] - Y[G];
	    if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	    __m256 v_deltaY = _mm256_broadcast_ss(&deltaY);
	    float *AscoreGJ = &FLAT(A,score, G, J);

#define BLOCK(H) {							\
	      __m256 v_prevscore = _mm256_maskload_ps(&AscoreGJ[-H], mask256); \
	      __m256 v_newscore = _mm256_add_ps(v_prevscore, Sint256(_mm256_maskload_ps(&delX[H][J], mask256), v_deltaY, H, I-G,c)); \
	      v_score = _mm256_max_ps(v_score, v_newscore);		\
	    }
	    
	    BLOCK(1);
	    BLOCK(2);
	    BLOCK(3);
	    BLOCK(4);
#undef BLOCK
	  }	
	  _mm256_maskstore_ps(&FLAT(A,score, I, J), mask256, v_score);
	  J = Jmax+1;
	}

  #endif // USE_SSE && USE_AVX along with USE_PFLOAT

  #if (USE_MIC && USE_PFLOAT)
	delYI = &delY[I][4];
	//Sint512_prefetch();

	for(; J <= Jmax-15; J+=16){
	  __m512 v_score = _mm512_loadu_ps(&FLAT(A,score, I, J));
#define DELX_BLOCK(H) __m512 delX##H = _mm512_loadu_ps(&FLAT_DELX(delX, H, J, M));
	
	  DELX_BLOCK(1);
	  DELX_BLOCK(2);
	  DELX_BLOCK(3);
	  DELX_BLOCK(4);
	
#undef DELX_BLOCK

#define BLOCKH(G, H) {							\
	    __m512 v_prevscore = _mm512_loadu_ps(&AscoreGJ[-H]);	\
	    __m512 v_newscore = _mm512_add_ps(v_prevscore, Sint512b(delX##H, v_deltaY, _mm512_set1_ps((float)(H+G)), c)); \
	    v_score = _mm512_gmax_ps(v_score, v_newscore);		\
	  }

#define BLOCKG(G) {					\
	    PFLOAT deltaY = delYI[-G];			\
	    __m512 v_deltaY = _mm512_set1_ps(deltaY);	\
	    float *AscoreGJ = &FLAT(A,score, I-G, J);	\
	  						\
	    BLOCKH(G, 1);				\
	    BLOCKH(G, 2);				\
	    BLOCKH(G, 3);				\
	    BLOCKH(G, 4);				\
	  }	
	
	  BLOCKG(1);
	  BLOCKG(2);
	  BLOCKG(3);
	  BLOCKG(4);

#undef BLOCK
#undef BLOCKG

	  _mm512_storeu_ps(&FLAT(A,score, I, J), v_score);
	} // for J <= Jmax-15
      
	if(J <= Jmax) {
	  //Sint512_prefetch();
	  __mmask16 mask16 = _mm512_int2mask((1<<(Jmax-J+1))-1); 
	  __m512 v_score = _mm512_maskloadu_ps(&FLAT(A,score, I, J), mask16);
	
#define DELX_BLOCK(H) __m512 delX##H = _mm512_maskloadu_ps(&FLAT_DELX(delX, H, J, M), mask16);
	
	  DELX_BLOCK(1);
	  DELX_BLOCK(2);
	  DELX_BLOCK(3);
	  DELX_BLOCK(4);
	
#undef DELX_BLOCK

#define BLOCK(G, H) {							\
	    __m512 v_prevscore = _mm512_maskloadu_ps(&AscoreGJ[-H], mask16); \
	    __m512 v_newscore = _mm512_add_ps(v_prevscore, Sint512b(delX##H, v_deltaY, _mm512_set1_ps((float)(H+G)), c)); \
	    v_score = _mm512_gmax_ps(v_score, v_newscore);		\
	  }

#define BLOCK2(G) {						\
	    __m512 v_deltaY = _mm512_set1_ps(delYI[-G]);	\
	    PFLOAT *AscoreGJ = &FLAT(A,score, I-G, J);	\
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
	  _mm512_maskstoreu_ps(&FLAT(A,score, I, J), mask16, v_score);
	  J = Jmax+1;
	} // if(J <= Jmax)

  #endif // USE_MIC && USE_PFLOAT

#endif // (USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)

      } else if (USE_MIC && USE_PFLOAT && DEFER_BP && !outlierExtend && !outlierBC && DELTA==4) {/* case with I <= 4 optimized with USE_MIC && USE_PFLOAT */
#if USE_MIC && USE_PFLOAT
	if(DEBUG>=2) assert(I <= 4);
      
	PFLOAT *delYI = delY[I];
	int H, G;
	PFLOAT deltaX[4], deltaY[4], score, tmp[16];
	int Gmin = max(I-DELTA,1);

	//Sint512_prefetch();
	__m512i Hvi, Gvi, Gv, Hv, IGv, JHv, Iv, Jv, Gv_stride;
	__m512 delXv, delYv, Av, score_v;
	__mmask16 mask, mask_IGv;
	Hvi = _mm512_set_epi32(0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3);
	Gvi = _mm512_set_epi32(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3);
	Gv = _mm512_add_epi32(_mm512_set1_epi32(Gmin), Gvi);
	IGv = _mm512_sub_epi32(_mm512_set1_epi32(I), Gv);
	Gv_stride = _mm512_mullo_epi32(_mm512_sub_epi32(Gv, _mm512_set1_epi32(1)), _mm512_set1_epi32(A->M_stride));
	mask_IGv = _mm512_cmp_epi32_mask(IGv, _mm512_setzero_epi32(), _MM_CMPINT_GT);
	delYv = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), mask_IGv, _mm512_sub_epi32(_mm512_set1_epi32(DELTA), IGv), delYI, 4);
	
	for(int J = Jmin; J <= Jmax; J++){
	  score = FLAT(A,score, I, J);
	  int Hmin = max(J-DELTA,1);

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
	  Av = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), mask, _mm512_add_epi32(Gv_stride, Hv), A->score[1], 4);
	  score_v = _mm512_mask_add_ps(_mm512_set1_ps(score), mask, Av, 
				       Sint512b(delXv, delYv, 
						_mm512_cvtfxpnt_round_adjustepu32_ps(
										     _mm512_add_epi32(JHv, IGv), _MM_FROUND_CUR_DIRECTION, _MM_EXPADJ_NONE), c));
	  PFLOAT newscore = _mm512_reduce_gmax_ps(score_v);
	  score = max(newscore,score);

	  FLAT(A,score, I, J) = score;
	}// for(J=Jmin; J <= Jmax; J++)
#endif // USE_MIC && USE_PFLOAT
      } else {/* general case without DEFER_BP etc */
	PFLOAT OutPen = MatchBiasOutlier[c] + OutlierPenalty;
	PFLOAT OutPen2 = OutPen + 1e-9;
	PFLOAT *delYI = delY[I];
	for(int J = Jmin; J <= Jmax; J++){
	  PFLOAT score = FLAT(A,score, I, J);
	  int bestG, bestH;
	  if(!DEFER_BP || outlierExtend){
	    PFLOAT deltaX, deltaY, newscore;
	    bestG = bestH = localtype;

	    int Gmin = max(I-DELTA,1), Hmin = max(J-DELTA,1);
	    for(int H= J; --H >= Hmin;){
	      deltaX = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];// X[J] - X[H]
	      if(DEBUG>=2) assert(fabs(X[J] - X[H] - deltaX) < 0.001);
	      for(int G= Gmin; G < I; G++){
		deltaY = delYI[G-I+DELTA];// Y[I] - Y[G] 
		if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
		newscore = FLAT(A,score, G, H) + Sint(deltaX,deltaY,J-H,I-G,outlierBC,c);
		if(EVERB && I==6495 && J==150 && G==6492 && H==149){
#pragma omp critical
		  {
		    printf("Yid=%d,Xid=%d,or=%d:I=%d,J=%d,G=%d,H=%d:deltaX=%0.4f,deltaY=%0.4f,A->score[G][H]=%0.8f,outlierBC=%d,Sint()=%0.8f,newscore=%0.8f,score=%0.8f\n",
			   Yid,Xid,orientation,I,J,G,H,deltaX,deltaY,FLAT(A,score,G,H),outlierBC,Sint(deltaX,deltaY,J-H,I-G,outlierBC),newscore,score);
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

	    if(DEBUG>=2) assert(bestG == localtype || bestG >= 1);

	    if(outlierExtend){
	      /* try to merge with any previous outlier interval that end at I or J */	    

	      /* try outliers that end at G==I && H < J */
	      int G = I;
	      if(G >= 2){
		for(int H = J; --H >= Hmin;){
		  int G2 = FLAT(A,G,G,H);
		  if(G2 > 0){
		    int H2 = FLAT(A,H,G,H);
		    if(DEBUG>=2) assert(G2 >= 1);
		    if(!outlierExtendLimit || max(I-G2,J-H2) <= outlierExtendLimit){
		      newscore = FLAT(A,score,G2,H2) + SintOutlier(X[J]-X[H2],Y[I]-Y[G2],OutPen2 + OutlierPenaltyBC[J-H2] + OutlierPenaltyBC[I-G2],c);
		      if(newscore > score){
			score = newscore;
			bestG = G2;
			bestH = H2;
		      }
		    }
		  }
		}
	      }

	      if(DEBUG>=2) assert(bestG == localtype || bestG >= 1);

	      /* try outliers that end at G < I && H==J */
	      int H = J;
	      if(H >= 2){
		for(int G = I; --G >= Gmin;){
		  int G2 = FLAT(A,G,G,H);
		  if(G2 > 0){
		    int H2 = FLAT(A,H,G,H);
		    if(DEBUG>=2) assert(G2 >= 1);
		    if(!outlierExtendLimit || max(I-G2,J-H2) <= outlierExtendLimit){
		      newscore = FLAT(A,score,G2,H2) + SintOutlier(X[J]-X[H2],Y[I]-Y[G2],OutPen2 + OutlierPenaltyBC[J-H2] + OutlierPenaltyBC[I-G2],c);
		      if(newscore > score){
			score = newscore;
			bestG = G2;
			bestH = H2;
		      }
		    }
		  }
		}
	      }

	      if(DEBUG>=2) assert(bestG == localtype || bestG >= 1);

	      if(bestG >= 1){

		PFLOAT pscore = FLAT(A,score, bestG, bestH);

		if(VERB>=3 && Yid==1 && Xid==0 && orientation==0 && I==14938 && J==38){
		  printf("A:Yid=%d,Xid=%d:I=%d,J=%d,bestG=%d,bestH=%d,score=%0.6f (pscore=%0.6f,OutPen=%0.6f)\n",
			 Yid,Xid,I,J,bestG,bestH,score,pscore,OutPen);
		  fflush(stdout);
		}

		if(score <= pscore + SintOutlier(X[J] - X[bestH], Y[I] - Y[bestG], OutPen2 + OutlierPenaltyBC[J - bestH] + OutlierPenaltyBC[I - bestG],c)
		   /* && alignFP(c,A, bestG, bestH, localtype,  Y, X, N, M, orientation, Yid, Xid, delY, delX, pscore, 0) > LogPvThreshold */){ /* try to extend outlier interval to the left by reducing bestG,bestH by up to outlierExtend sites */

		  int origbestG = bestG, origbestH = bestH;

		  for(int cnt = 0; cnt < 100; cnt++){

		    int Gmin2 = max(1, origbestG - outlierExtend);
		    int Hmin2 = max(1, origbestH - outlierExtend);
		    if(outlierExtendLimit){
		      Gmin2 = max(I - outlierExtendLimit, Gmin2);
		      Hmin2 = max(J - outlierExtendLimit, Hmin2);
		    }
		    for(int G = origbestG; G >= Gmin2; G--){
		      PFLOAT deltaY = Y[I] - Y[G];
		      PFLOAT *AscoreG = &FLAT(A,score,G,0);
		      for(int H = origbestH; H >= Hmin2; H--){
			PFLOAT newscore = AscoreG[H] + SintOutlier(X[J] - X[H], deltaY, OutPen + OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-G], c);
			if(newscore > score){
			  score = newscore;
			  bestG = G;
			  bestH = H;
			}
		      }
		    }

		    if(DEBUG>=2) assert(bestG <= origbestG && bestH <= origbestH);
		    if(DEBUG>=2) assert(bestG == localtype || bestG >= 1);

		    if(VERB>=3 && Yid==1 && Xid==0 && orientation==0 && I==14938 && J==38 && bestG >= 1){
		      printf("B:  cnt=%d:Yid=%d,Xid=%d:I=%d,J=%d,bestG=%d,bestH=%d,score=%0.6f (pscore=%0.6f,OutPen=%0.6f)\n",
			     cnt,Yid,Xid,I,J,bestG,bestH,score,FLAT(A,score,bestG,bestH),OutPen);
		      fflush(stdout);
		    }

		    /* also check if best interval before (G,H) can be merged with this outlier interval : include checking G==I (if H <= min(J-1,origbestH)) and H==J (with G <= min(I-1,origbestG)) */
		    for(int G = I; G >= Gmin2; G--){
		      if(G < I) G = min(origbestG,G);
		      int Hmax2 = (G <= min(I-1,origbestG)) ? J : min(J-1,origbestH);
		      for(int H = Hmax2; H >= Hmin2; H--){
			if(H < J) H = min(origbestH,H);
			if(DEBUG>=2) assert(G <= min(I-1,origbestG) || H <= min(J-1,origbestH));
			int G2 = FLAT(A,G, G, H);
			if(G2 > 0){
			  int H2 = FLAT(A,H, G, H);
			  if(DEBUG>=2 && !(G2 < G || H2 < H)){
			    printf("Yid=%d(%lld),Xid=%d(%lld):I=%d,J=%d:G=%d,G2=%d,H=%d,H2=%d,origbestG=%d,origbestH=%d\n",
				   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,I,J,G,G2,H,H2,origbestG,origbestH);
			    fflush(stdout);
			    assert(G2 < G || H2 < H);
			  }
			  if(!outlierExtendLimit || max(I-G2,J-H2) <= outlierExtendLimit){
			    newscore = FLAT(A,score, G2, H2) + SintOutlier(X[J]-X[H2],Y[I]-Y[G2],OutPen + OutlierPenaltyBC[J - H2] + OutlierPenaltyBC[I - G2],c);
			    if(newscore > score + 1e-6){
			      /*			  if(G2==I && H2 == J-1){
							  printf("Yid=%d(%lld),Xid=%d(%lld):I=%d,J=%d:G=%d,G2=%d,H=%d,H2=%d,bestG=%d,origbestG=%d,bestH=%d,origbestH=%d:score=%0.6f->%0.6f(Sint=%0.6f,OutPen=%0.6f)\n",
							  Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,I,J,G,G2,H,H2,bestG,origbestG,bestH, origbestH,score,newscore,
							  SintOutlier(X[J]-X[H2],Y[I]-Y[G2],OutPen+OutlierPenaltyBC[J-H2]+PoutlierPenaltyBC[I-G2]),OutPen);
							  fflush(stdout);
							  }*/
			      score = newscore;
			      bestG = G2;
			      bestH = H2;
			    }
			  }
			}
		      }
		    }

		    if(DEBUG>=2) assert(bestG == localtype || bestG >= 1);

		    if(origbestG == bestG && origbestH == bestH)
		      break;

		    if(VERB>=3 && Yid==1 && Xid==0 && orientation==0 && I==14938 && J==38 && bestG >= 1){
		      printf("C:  cnt=%d:Yid=%d,Xid=%d:I=%d,J=%d,bestG=%d,bestH=%d,score=%0.6f (pscore=%0.6f,OutPen=%0.6f)\n",
			     cnt,Yid,Xid,I,J,bestG,bestH,score,FLAT(A,score,bestG,bestH),OutPen);
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
		  }/* while(1) */
		} else {/* try to merge non-outlier interval with previous interval (in case previous interval is an outlier) */
		  int G = FLAT(A,G,bestG,bestH);
		  if(G > 0){
		    int H = FLAT(A,H,bestG,bestH);
		    if(!outlierExtendLimit || max(I-G,J-H) <= outlierExtendLimit){
		      PFLOAT newscore = FLAT(A,score,G,H) + SintOutlier(X[J]-X[H],Y[I]-Y[G],OutPen + OutlierPenaltyBC[I-G] + OutlierPenaltyBC[J-H],c);
		      if(newscore > score){
			score = newscore;
			bestG = G;
			bestH = H;
			if(DEBUG>=2) assert(bestG == localtype || bestG >= 1);

			if(VERB>=3 && Yid==1 && Xid==0 && orientation==0 && I==14938 && J==38){
			  printf("D:  Yid=%d,Xid=%d:I=%d,J=%d,bestG=%d,bestH=%d,score=%0.6f (pscore=%0.6f,OutPen=%0.6f)\n",
				 Yid,Xid,I,J,bestG,bestH,score,FLAT(A,score,bestG,bestH),OutPen);
			  fflush(stdout);
			}
		      }
		    }
		  }
		}
	      }
	    }

	    if(DEBUG>=2) assert(isfinite(score));
	    if(DEBUG>=2) assert(bestG == localtype || bestG >= 1);

	    FLAT(A,score, I, J) = score;
	    FLAT(A,G, I, J) = bestG;
	    FLAT(A,H, I, J) = bestH;
	  } else {/* no need to track bestG,bestH */
	    int Gmin = max(I-DELTA,1), Hmin = max(J-DELTA,1);
	    PFLOAT deltaX, deltaY, newscore;

	    for(int H= J; --H >= Hmin;){
	      deltaX = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];// X[J] - X[H]
	      if(DEBUG>=2 && !(fabs(X[J] - X[H] - deltaX) < 0.001)){
		#pragma omp critical
		{
		  if((USE_SSE && USE_AVX) || USE_MIC)
		    printf("Yid=%d,Xid=%d(id=%lld),or=%d,c=%d:J=%d,H=%d,M=%d:deltaX= %0.6f, X[J]=%0.6f,X[H]=%0.6f:X[J]-X[H]=%0.6f, FLAT_DELX(delX, J-H, J, M)=%0.6f\n",
			   Yid,Xid,XXmap[Xid]->id,orientation,c,J,H,M,deltaX,X[J],X[H],X[J]-X[H],FLAT_DELX(delX,H-H,J,M));
		  else
		    printf("Yid=%d,Xid=%d(id=%lld),or=%d,c=%d:J=%d,H=%d,M=%d:deltaX= %0.6f, X[J]=%0.6f,X[H]=%0.6f:X[J]-X[H]=%0.6f, delX[J][H-J+DELTA]=%0.6f\n",
			   Yid,Xid,XXmap[Xid]->id,orientation,c,J,H,M,deltaX,X[J],X[H],X[J]-X[H],delX[J][H-J+DELTA]);
		  fflush(stdout);
		  assert(fabs(X[J] - X[H] - deltaX) < 0.001);
		}
	      }
	      for(int G= Gmin; G < I; G++){
		deltaY = delYI[G-I+DELTA];// Y[I] - Y[G] 
		if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
		newscore = FLAT(A,score, G, H) + Sint(deltaX,deltaY,J-H,I-G,c);
		score = PMAX(score,newscore);
	      }
	    }
	    FLAT(A,score, I, J) = score;
	  }
	} // J = Jmin .. Jmax
      } // general case without DEFER_BP etc
    } // I = 2 .. N
  }// c = 0..colors-1

  /* locate best alignment right ends for both colors simultaneously */
  /* If ENDINIT only check results where I >= N+1-DELTA OR J >= M+1-DELTA */
  if(DEBUG>=2) assert(colors==2);

  int bestI[2]= {-1,-1}, bestJ[2] = {-1,-1};/* index refers to color */
  int bestR[2] = {localtype,localtype};/* bestR[c] = right end alignment type (defaults to local) */
  //  double bestYi = 0.0;
  PFLOAT bestscore = ScoreThreshold; // MINSCORE;

  /* compute range of sites Qmin[J]..Qmax[J] for 2nd color of XX that is near site J for 1st color of XX */
  int *Qmin = (int *)alloca((MM[0]+1)*sizeof(int));
  int *Qmax = (int *)alloca((MM[0]+1)*sizeof(int));
  double Yrange = (DELTA+1)*Ytheta12;
  double Xrange = Yrange;
  int qmin = 1, qmax = 0;

  for(int J = 1; J <= MM[0]; J++){
    FLOAT Xj = XX[0][J];
    FLOAT *X = XX[1];
    int M = MM[1];

    FLOAT Xmin = Xj - Xrange;
    while(qmin <= M && X[qmin] < Xmin)
      qmin++;
    qmax = max(qmin-1,qmax);

    FLOAT Xmax = Xj + Xrange;
    while(qmax < M && X[qmax+1] <= Xmax)
      qmax++;

    Qmin[J] = qmin;
    Qmax[J] = qmax;
  }

  int pmin = 1, pmax = 0;
  for(int I = 1; I <= NN[0]; I++){
    double bscore = MINSCORE;
    int bJ = -1, bR = -1, bP = -1, bQ = -1, bS = -1;

    int *RijxI = &FLAT(&AA[0], Rijx, I, 0), *RijyI = &FLAT(&AA[0], Rijy, I, 0);
    PFLOAT YdiffI = Yrc[0][I];

    /* update pmin,pmax */
    double Ymin = YY[0][I] - Yrange;
    double Ymax = YY[0][I] + Yrange;
    while(pmin <= NN[1] && YY[1][pmin] < Ymin)
      pmin++;
    pmax = max(pmin-1,pmax);
    while(pmax < NN[1] && YY[1][pmax+1] <= Ymax)
      pmax++;

    FLOAT Yi = YY[0][I];
    AL *p = &AA[0];
    for(int J = 1; J <= MM[0]; J++){
      FLOAT Xj = XX[0][J];
      for(int P = pmin; P <= pmax; P++){
	if(DEBUG>=2) assert(YY[1][P] >= Ymin && YY[1][P] <= Ymax);

	int *RijxP = &FLAT(&AA[1], Rijx, P, 0);
	int *RijyP = &FLAT(&AA[1], Rijy, P, 0);
	PFLOAT YdiffP = Yrc[1][P];

	FLOAT y1 = YY[1][P];
	PFLOAT y = y1 - Yi;

	AL *q = &AA[1];
	for(int Q = Qmin[J]; Q <= Qmax[J]; Q++){
	  FLOAT x1 = XX[1][Q];
	  PFLOAT x = x1 - Xj;
	  PFLOAT len;
	  if(COLOR_SHIFT)
	    len = max(x1, Xj) + max(y1,Yi);
	  double Lscore = FLAT(p,score,I,J) + FLAT(q,score,P,Q) + SMA(x,y, len, orientation);

	  if(EVERB && !orientation && I==62 && J==47 && P==48 && Q==44){
	    #pragma omp critical
	    {
	      printf("Yid=%d,Xid=%d,or=%d:I=%d,J=%d,P=%d,Q=%d,N=%d,%d,M=%d,%d:AA[0].score[I][J]= %0.6f, AA[1].score[P][Q]= %0.6f, x=%0.4f,y=%0.4f,len=%0.4f:SMA(x,y,len,or)=%0.6f: Lscore=%0.6f,bestscore=%0.6f(I=%d,%d,J=%d,%d,R=%d,%d)\n",
		     Yid,Xid,orientation,I,J,P,Q,NN[0],NN[1],MM[0],MM[1],AA[0].score[I][J],AA[1].score[P][Q], x, y, len, SMA(x,y,len,orientation), Lscore, bestscore, bestI[0],bestI[1],bestJ[0],bestJ[1],bestR[0],bestR[1]);
	      if(I==62 && J==47 && P==48 && Q==48)
		printf("    EndThresh= %0.6f + %0.6f\n",EndThresh[0],EndThresh[1]);
	      fflush(stdout);
	    }
	  }

	  PFLOAT thresh;
	  if(SEND_THRESH && EndThresh[0]+EndThresh[1] < (thresh = bestscore - Lscore))
	    continue;

	  PFLOAT Rscore = C, Sscore = C;/* default to unaligned ends */
	  int R = localtype, S = localtype;

	  /* compute right end score for aligned end for color 0 */
	  int N = NN[0], M = MM[0];
	  FLOAT *Y = YY[0], *X = XX[0];
	  PFLOAT *Xr = Xrc[0];
	  if(I >= N+1-DELTA || J >= M+1-DELTA){
	    int Rijx = RijxI[J], Rijy = RijyI[J];
	    if(DEBUG>=2) assert(Rijx > M || Rijy > N);
	    if(EVERB && !orientation && I==62 && J==47 && P==48 && Q==48){
	      #pragma omp critical
	      {
		printf("    AA[0].Rijx[I][J]=%d,AA[0].Rijy[I][J]=%d:M=%d,N=%d,YrightClipped=%d,XrightClipped=%d\n",AA[0].Rijx[I][J],AA[0].Rijy[I][J],M,N,YrightClipped,XrightClipped);
		fflush(stdout);
	      }
	    }

	    if(!((Rijy > N && YrightClipped) || (Rijx > M && XrightClipped))){
	      PFLOAT score = Send(PMIN(Xr[J], YdiffI), (PFLOAT)(min(M, Rijx)+1-J + min(N, Rijy)+1-I), 0);
	      if(ENDFIX){/* update Rijx,Rijy based on Sbnd() */
		if(!(SEND_THRESH && score < thresh - (EndThresh2[0]+EndThresh2[1]))){
		  if(Rijx > M) {
		    double a = M+1-I-J;
		    int K;
		    for(K= Rijy; K > I; K--){
		      PFLOAT bound = Sbnd(Xr[J],Y[K]-Y[I], K+a, 0);
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
		      PFLOAT bound = Sbnd(X[K]-X[J], YdiffI,K+a, 0);
		      if(bound <= score)
			break;// HERE : see ENDFIX2 in refalign.cpp
		      score = bound;
		    }
		    RijxI[J] = Rijx = K;
		  }
		}
	      }
	      if(EVERB && !orientation && I==62 && J==47 && P==48 && Q==48){
	        #pragma omp critical
		{
		  printf("    score=%0.6f,AA[0].Rijx[I][J] -> %d, AA[0].Rijy[I][J] -> %d\n",score, AA[0].Rijx[I][J],AA[0].Rijy[I][J]);
		  fflush(stdout);
		}
	      }

	      if(score > Rscore){
		Rscore = score;
		R = -1;/* end type : -1 = normal, 0 = chr end, -2 = local, -3 = outlier */
	      }
	    }
	  }
	  
	  /* compute right end score for aligned end for color 1 */
	  N = NN[1]; M = MM[1];
	  Y = YY[1]; X = XX[1];
	  Xr = Xrc[1];
	  if(P >= N+1-DELTA || Q >= M+1-DELTA){
	    int Rijx = RijxP[Q], Rijy = RijyP[Q];
	    if(DEBUG>=2) assert(Rijx > M || Rijy > N);
	    if(EVERB && !orientation && I==62 && J==47 && P==48 && Q==48){
	      #pragma omp critical
	      {
		printf("    AA[1].Rijx[P][Q]=%d,AA[1].Rijy[P][Q]=%d:M=%d,N=%d,YrightClipped=%d,XrightClipped=%d\n",AA[1].Rijx[P][Q],AA[1].Rijy[P][Q],M,N,YrightClipped,XrightClipped);
		fflush(stdout);
	      }
	    }
	    if(!((Rijy > N && YrightClipped) || (Rijx > M && XrightClipped))){
	      PFLOAT score = Send(PMIN(Xr[Q], YdiffP), (PFLOAT)(min(M, Rijx)+1-Q + min(N, Rijy)+1-P), 1);
	      if(ENDFIX){/* update Rijx,Rijy based on Sbnd() */
		if(!(SEND_THRESH && Rscore + score < thresh - EndThresh2[1])){
		  if(Rijx > M) {
		    double a = M+1-P-Q;
		    int K;
		    for(K= Rijy; K > P; K--){
		      PFLOAT bound = Sbnd(Xr[Q],Y[K]-Y[P], K+a, 1);
		      if(bound <= score)
			break;// HERE : see ENDFIX2 in refalign.cpp
		      score = bound;
		    }
		    RijyP[Q] = Rijy = K;
		  }
		  if(Rijy > N) {
		    double a = N+1-P-Q;
		    int K;
		    for(K= Rijx;K > Q;K--){
		      PFLOAT bound = Sbnd(X[K]-X[Q], YdiffP,K+a, 1);
		      if(bound <= score)
			break;// HERE : see ENDFIX2 in refalign.cpp
		      score = bound;
		    }
		    RijxP[Q] = Rijx = K;
		  }
		}
	      }
	      if(EVERB && !orientation && I==62 && J==47 && P==48 && Q==48){
	        #pragma omp critical
		{
		  printf("    score=%0.6f,AA[1].Rijx[P][Q] -> %d, AA[1].Rijy[P][Q] -> %d\n",score, AA[1].Rijx[P][Q],AA[1].Rijy[P][Q]);
		  fflush(stdout);
		}
	      }
	      if(score > Sscore){
		Sscore = score;
		S = -1;/* end type : -1 = normal, 0 = chr end, -2 = local, -3 = outlier */
	      }
	    }
	  }
	  if(DEBUG>=2) assert(isfinite(Rscore));
	  if(DEBUG>=2) assert(isfinite(Sscore));
	  FLAT(&AA[0],Rscore,I,J) = Rscore;
	  FLAT(&AA[1],Rscore,P,Q) = Sscore;

	  double nscore = Lscore + Rscore + Sscore;
	  if((EVERB && !orientation && I==62 && J==47 && P==48 && Q==44) || (VERB>=3 && nscore > bscore && Yid==0 && Xid==204 && nscore > 0.0)){
	    printf("Yid=%d,Xid=%d,or=%d:nscore=%0.6f,I=%d,J=%d,P=%d,Q=%d:Rscore=%0.6e,Sscore=%0.6e,R=%d,S=%d:bscore=%0.6f(bJ=%d,bP=%d,bQ=%d,bR=%d,bS=%d)\n",
		   Yid,Xid,orientation,nscore,I,J,P,Q,Rscore,Sscore,R,S,bscore,bJ,bP,bQ,bR,bS);
	    fflush(stdout);
	  }
	  if(nscore > bscore){
	    bscore = nscore;
	    bJ = J;
	    bP = P;
	    bQ = Q;
	    bR = R;
	    bS = S;
	  }
	}
      }
    }

    /* HERE HERE : need to allow for cases where one color has no aligned sites (just use Send() for color without aligned sites + trivial SMA()) : this will only be better if there is no sites at all to align for one color */

    if((VERB>=3 || EVERB) && bscore > bestscore + 1e-8){
      #pragma omp critical
      {
	printf("result for Yid=%d,Xid=%d,or=%d:I=%d,YY[0][I]=%0.3f:bscore=%0.6f,bJ=%d,bR=%d,bP=%d,bQ=%d,bS=%d(bestscore=%0.6f)\n",
	       Yid,Xid,orientation,I,YY[0][I],bscore,bJ,bR,bP,bQ,bS,bestscore);
	fflush(stdout);
      }
    }

    //    double bYi = (bJ >= 1 ) ? YY[0][I] : -1e30;

    if(bscore > bestscore + 1e-8){
      bestscore = bscore;
      bestI[0] = I;
      bestJ[0] = bJ;
      bestR[0] = bR;
      //      bestYi = bYi;
      bestI[1] = bP;
      bestJ[1] = bQ;
      bestR[1] = bS;
      if(DEBUG){
	for(int c = 0; c < colors; c++){
	  assert(bestI[c] >= 1 && bestI[c] <= NN[c]);
	  assert(bestJ[c] >= 1 && bestJ[c] <= MM[c]);
	  assert(bestR[c] <= -2 || bestJ[c] >= MM[c]+1 - DELTA_X || bestI[c] >= NN[c]+1-DELTA_Y);
	}
      }
    }
  }

  if(bestI[0] < 0 && bestI[1] < 0){/* no valid alignment found : should never happen with -S -1e30 */
    if(VERB >= (EVERB ? 1 : 2)){
      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d:No alignment found (bestscore=%0.3f,bestI=%d,%d)\n",
	     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,bestscore,bestI[0],bestI[1]);
      fflush(stdout);
    }
    return;
  }

  if(DEFER_BP){/* compute G,H values along best alignment path */
    for(int c = 0; c < colors; c++){
      FLOAT *Y = YY[c];
      FLOAT *X = XX[c];
      //      int N = NN[c];
      int M = MM[c];
      AL *A = &AA[c];
      PFLOAT **delY = delYY[c];
      PFLOAT **delX = delXX[c];
      PFLOAT *Xf = Xfc[c];
      PFLOAT *Yf = Yfc[c];

      int I = bestI[c];
      int J = bestJ[c];
      int bestG = -1, bestH = -1;

      PFLOAT score;
      while(I >= 0){
	if(DEBUG>=2) assert(bestI[c] > 0 && bestJ[c] > 0);
	if(outlierExtend){/* A->G and A->H has already been computed */
	  score = FLAT(A,score,I,J);
	  bestG = FLAT(A,G,I,J);
	  bestH = FLAT(A,H,I,J);
	  if(DEBUG>=2) assert(bestG == localtype || bestG >= 1);
	} else {
	  /* update A->G[I][J] & A->H[I][J] by repeating recurrance computation */
	  score = FLAT(A,Uscore, I, J);
	  if(DEBUG>=2) assert(isfinite(score));
	  bestG = localtype;

	  int Gmin = max(I-DELTA,1), Hmin = max(J-DELTA,1);
	  for(int G= I; --G >= Gmin;){
	    PFLOAT deltaY = delY[I][G-I+DELTA];// Y[I] - Y[G] 
	    if(DEBUG>=2) assert(fabs(Y[I]-Y[G] - deltaY) < 0.001);
	    PFLOAT *AscoreG = &FLAT(A,score, G, 0);
	    PFLOAT a = J+I-G;
	  
	    for(int H= J; --H >= Hmin;){
	      PFLOAT deltaX = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];// X[J] - X[H]
	      if(DEBUG>=2) assert(fabs(X[J] - X[H] - deltaX) < 0.001);
	      PFLOAT newscore = AscoreG[H] +  Sint(deltaX,deltaY,(PFLOAT)(a-H),c);
	      if(EVERB && I==6495 && J==150 && G==6492 && H==149){
                #pragma omp critical
		{
		  printf("Yid=%d,Xid=%d,or=%d,c=%d:I=%d,J=%d,G=%d,H=%d:deltaX=%0.4f,deltaY=%0.4f,A->score[G][H]=%0.8f,Sint()=%0.8f,newscore=%0.8f,score=%0.8f\n",
			 Yid,Xid,orientation,c,I,J,G,H,deltaX,deltaY,AscoreG[H],Sint(deltaX,deltaY,(PFLOAT)(a-H),c),newscore,score);
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
	  if(DEBUG>=2 && !(fabs(score - A->score[I][J]) < (PFLOAT)(USE_PFLOAT ? (bestscore > 0.0 ? 1e-3 : 1e-2) :1e-8))){
            #pragma omp critical
	    {
	      printf("Yid=%d,Xid=%d,orientation=%d,c=%d:I=%d,J=%d,G=%d,H=%d:score=%0.8f,A->score[I][J]=%0.8f,bestscore=%0.3f\n",
		     Yid,Xid,orientation,c,I,J,bestG,bestH,score,A->score[I][J],bestscore);
	      fflush(stdout);
	      assert(fabs(score - FLAT(A,score, I, J)) < (PFLOAT)(USE_PFLOAT ? (bestscore > 0.0 ? 1e-3 : 1e-2):1e-8));
	    }
	  }
	}

	if(bestG < 0){/* value of bestH is undefined */
	  if(DEFER_G>=1 && score > C){
	    bestG = -1;
	    if(DEFER_G >= 2){/* check if actual value of Lijx,Lijy needs to be incremented */
	      if(DEBUG>=2) assert(I <= DELTA || J <= DELTA);
	      int Lijx = FLAT(A,Lijx,I,J);
	      int Lijy = FLAT(A,Lijy,I,J);
	      if(DEBUG>=2) assert(0 <= Lijx && Lijx <= J);
	      if(DEBUG>=2) assert(0 <= Lijy && Lijy <= I);
	      if(DEBUG>=2) assert(Lijx <= 0 || Lijy <= 0);
	      if((Lijx <= 0) ? !XleftClipped : !YleftClipped){
		PFLOAT SmIJ = Sm(J,I);
		PFLOAT ascore = C;
		
		if(Lijx <= 0){
		  ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I+1-max(1,Lijy)),c) + SmIJ;
		  ascore = max(ascore,C);
		  for(int K= Lijy; K < I; ) {
		    PFLOAT bound = Sbnd(X[J],Y[I]-Y[K],(PFLOAT)(J+I-K),c) + SmIJ;
		    if(bound <= ascore)
		      break;// HERE HERE : see ENDFIX2 in refalign.cpp
		    ascore = bound;
		    Lijy = ++K;
		  }
		  FLAT(A,Lijy, I, J) = Lijy;
		} else {
		  if(DEBUG>=2) assert(Lijy <= 0);
		  ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-max(1,Lijx)+I),c) + SmIJ;
		  ascore = max(ascore,C);
		  for(int K= Lijx; K < J; ){
		    PFLOAT bound = Sbnd(X[J]-X[K],Y[I],(PFLOAT)(J+I-K),c) + SmIJ;
		    if(bound <= ascore)
		      break;// HERE HERE : see ENDFIX2 in refalign.cpp
		    ascore = bound;
		    Lijx = ++K;
		  }
		  FLAT(A,Lijx, I, J) = Lijx;
		}
	      }
	    }
	  }
	  FLAT(A,G, I, J) = bestG;
	  if(DEBUG>=2 && !(Ch > 0.0 || PoutlierEnd > 0.0)) assert(A->G[I][J] > -2);
	} else {
	  if(DEBUG>=2) assert(bestG > 0 && bestH > 0);
	  FLAT(A,G, I, J) = bestG;
	  FLAT(A,H, I, J) = bestH;
	  J = bestH;
	}
	I = bestG;
      } // while(I >= 0)
    } // c = 0 .. colors -1
  } // if(DEFER_BP)

  if(VERB >= (EVERB ? 1 : 2) && !(bestscore  > align->score)){
    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d:bestscore=%0.4f,bestI=%d,%d,bestJ=%d,%d:align[-1].score=%0.4f,align[-1].logPV=%0.2f,align[-1].numpairs=%d\n",
	   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,bestscore,bestI[0],bestI[1],bestJ[0],bestJ[1],align[-1].score,align[-1].logPV,align[-1].numpairs);
    fflush(stdout);
  }

  double origscore = align[-1].score;
  int KK[2] = {-1,-1}, bestIL[2] = {-1,-1}, bestJL[2] = {-1,-1};

  if(bestscore > origscore || VERB >= (EVERB ? 1 : 2)){ /* backtrack through array to retrieve best alignment starting at right end A[bestI][bestJ] */
    /* allocate space to store alignment in reverse order */
    int M = max(MM[0],MM[1]);
    int *Ilist = (int *)alloca(M*sizeof(int));
    int *Jlist = (int *)alloca(M*sizeof(int));
    FLOAT *iscore = (FLOAT *)alloca((M+1)*sizeof(FLOAT));
    FLOAT *outscore = (FLOAT *)alloca((M+1)*sizeof(FLOAT));

    if(DEBUG>=2)
      for(int U = 0; U <= M; U++)
	outscore[U] = iscore[U] = aNaN;

    for(int c = 0; c < colors; c++){
      FLOAT *Y = YY[c];
      FLOAT *X = XX[c];
      int N = NN[c];
      int M = MM[c];
      AL *A = &AA[c];

      PFLOAT **delY = delYY[c];
      PFLOAT **delX = delXX[c];

      PFLOAT *Xf = Xfc[c];
      PFLOAT *Yf = Yfc[c];
      PFLOAT *Xr = Xrc[c];
      PFLOAT *Yr = Yrc[c];

      int I = bestI[c];
      int J = bestJ[c];

      int K = 0;
    
      /* right end segment */
      outscore[0] = iscore[0] = FLAT(A,Rscore,I,J);
      if(DEBUG>=2 && !(isfinite(outscore[0]))){
	#pragma omp critical
	{
	  printf("Yid=%d,Xid=%d,or=%d,c=%d:bestI=%d,bestJ=%d,A->Rscore[I][J]=%0.6e,outscore[0]=%0.6e, bestscore=%0.6f\n",
		 Yid,Xid,orientation,c,bestI[c],bestJ[c],A->Rscore[I][J],outscore[0],bestscore);
	  fflush(stdout);
	  assert(isfinite(outscore[0]));
	}
      }

      if(EVERB){
        #pragma omp critical
        {
	  printf("Yid=%d,Xid=%d,or=%d,c=%d:bestI=%d,bestJ=%d:iscore[0]=%0.6f\n",Yid,Xid,orientation,c,bestI[c],bestJ[c],iscore[0]);
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
	int G = FLAT(A,G, I, J);
	int H = FLAT(A,H, I, J);
	if(G > 0){
	  if(DEBUG>=2) assert(H > 0);
	  iscore[K+1] = FLAT(A,score, I, J) - FLAT(A,score, G, H);
	  PFLOAT y = Y[I] - Y[G];
	  PFLOAT x = X[J] - X[H];
	  if(Poutlier > 0.0){
	    PFLOAT Bias,Pen;
	    SintDetail(x, y, J - H, I - G, J, I, Bias, Pen,c);
	    outscore[K+1] = Bias + Pen;
	    if(EVERB || DEBUG>=2){
	      PFLOAT OutPen = OutlierPenalty;
	      if(outlierBC)
		OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-G];
	      PFLOAT iscoreB = Bias + OutlierBias + (OUTLIER_DELTA(x-y) ? max(Pen,OutPen) : Pen);	    

	      if(DEBUG>=2 && !(fabs(iscoreB - iscore[K+1]) < 0.001)){
	        #pragma omp critical
	        {
		  printf("Yid=%d,Xid=%d,or=%d,c=%d:bestI=%d,bestJ=%d:OutlierBias=%0.6f,OutlierPenalty=%0.6f\n",Yid,Xid,orientation,c,bestI[c],bestJ[c],OutlierBias,OutlierPenalty);
		  printf("  K=%d:I=%d,J=%d,G=%d,H=%d:x=%0.4f,y=%0.4f,Bias=%0.6f,Pen=%0.6f,iscore[K+1]=%0.6f,outscore[K+1]=%0.6f,iscoreB=%0.6f\n",
			 K,I,J,G,H,x,y,Bias,Pen,iscore[K+1],outscore[K+1],iscoreB);
		  fflush(stdout);
		  assert(fabs(iscoreB - iscore[K+1]) < 0.001);
		}
	      }
	    }
	  } else
	    outscore[K+1] = iscore[K+1];
	} else {/* left end segment */
	  iscore[K+1] = FLAT(A,score, I, J);
	  outscore[K+1] = iscore[K+1];
	}
	K++;
	I = G;
	J = H;
      }
      KK[c] = K;
      if(DEBUG && bestscore > origscore && !(K > 0)){
	#pragma omp critical
	{
	  printf("Yid=%d,Xid=%d,or=%d,c=%d:K=%d,bestI[c]=%d,bestJ[c]=%d,bestscore=%0.6f,origscore=%0.8e\n",Yid,Xid,orientation,c,K,bestI[c],bestJ[c],bestscore,origscore);
	  fflush(stdout);
	  assert(K>0);
	}
      }

      if(bestscore > origscore){
	if(DEBUG>=2 && !(Ch > 0.0 || PoutlierEnd > 0.0) && !(bestR[c] > -2)){
	  #pragma omp critical
	  {
	    printf("Yid=%d,Xid=%d,or=%d,c=%d:K=%d,bestI[c]=%d,bestJ[c]=%d,bestR[c]=%d,bestscore=%0.6f,origscore=%0.8e,ChimScore=%0.8e\n",Yid,Xid,orientation,c,K,bestI[c],bestJ[c],bestR[c],bestscore,origscore,ChimScore);
	    fflush(stdout);
	    assert(bestR[c] > -2);
	  }
	}

	align[c].score = bestscore;
	align[c].logPV = alignFP(c,A,bestI[c],bestJ[c],bestR[c],Y,X,N,M,orientation,Yid,Xid,delYY[c],delXX[c],bestscore,0, 1);
	align[c].orientation = orientation;
	align[c].Lend = I;
	align[c].Rend = bestR[c];

	I = Ilist[K-1];// align->sites1[0];
	J = Jlist[K-1];// align->sites2[0];

	align[c].Lij1 = (FLAT(A,G, I, J) <= -2) ? I : FLAT(A,Lijy, I, J);
	align[c].Lij2 = (FLAT(A,G, I, J) <= -2) ? J : FLAT(A,Lijx, I, J);
	align[c].Rij1 = (bestR[c] <= -2) ? bestI[c] : FLAT(A,Rijy, bestI[c], bestJ[c]);
	align[c].Rij2 = (bestR[c] <= -2) ? bestJ[c] : FLAT(A,Rijx, bestI[c], bestJ[c]);

	/* (re)allocate best alignment */
	align[c].allfree();/* resets align->numpairs = 0 (if needed : typically not needed) */
	if(DEBUG>=2) assert(align[c].score > MINSCORE);
	align[c].numpairs = K;
	align[c].expand_arrays(K);

	if(DEBUG) assert(min(M,align[c].Rij2) >= max(1,align[c].Lij2) + align[c].numpairs - 1);
	if(DEBUG) assert(min(N,align[c].Rij1) >= max(1,align[c].Lij1) + align[c].numpairs - 1);

	/* copy alignment from Ilist[],Jlist[] in reverse order */
        #pragma ivdep
	for(int I= 0; I < K; I++){
	  align[c].sites1[I] = Ilist[K-1-I];
	  align[c].sites2[I] = Jlist[K-1-I];
	  align[c].iscore[I] = iscore[K-I];
	  align[c].outscore[I] = outscore[K-I];
	  if(DEBUG>=2) assert(isfinite(align[c].outscore[I]));
	}
	align[c].iscore[K] = iscore[0];
	align[c].outscore[K] = outscore[0];
	if(DEBUG>=2) assert(isfinite(align[c].outscore[K]));

	align[c].repeat = -1;

	/* check if bestscore is a repeat (for color c) */
	if(RepeatRec && RepeatMaxShift > 0 && RepeatPvalueRatio > 0.0 && align[c].numpairs >= AlignedSiteThreshold){
	  // if(DEBUG) assert(!PairSplit);
	  if(VERB && (VERB>=2 || EVERB)){
            #pragma omp critical
	    {
	      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:M=%d,N=%d:checking for repeats\n",
		     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,M,N);
	      fflush(stdout);
	    }
	  }

	  double minPvalue = LogPvThreshold + log(RepeatPvalueRatio)/log(10.0);

	  int *Ilist2 = (int *)alloca(M*sizeof(int));
	  int *Jlist2 = (int *)alloca(M*sizeof(int));

	  align[c].repeat = 0;

	  int bestI2 = bestI[c];
	  int bestJ2 = bestJ[c];

	  for(bestJ2 = max(1, align[c].sites2[max(0,K-1-RepeatMaxShift)]); !align[c].repeat && bestJ2 < bestJ[c]; bestJ2++){/* back track best alignment from (bestI2,bestJ2) and compute its Pvalue */	
	    // WAS  for(bestJ2 = max(1,bestJ[c] - RepeatMaxShift); !align[c].repeat && bestJ2 < bestJ[c]; bestJ2++){/* back track best alignment from (bestI2,bestJ2) and compute its Pvalue */
	    double MinShiftKb = (X[bestJ[c]]-X[bestJ2]) * RepeatKbShiftRatio; /* minimum offset shift along entire alignment */
	    double MaxShiftKb = (X[bestJ[c]]-X[bestJ2]) * RepeatKbShiftRatio2;/* maximum offset shift along entier alignment (to avoid creating new outliers) */
	    int I = bestI2;
	    int J = bestJ2;

	    /* First compute bestscore2, bestR2 with rightmost alignment at (I,J) */
	    int bestR2 = localtype;
	    PFLOAT bestscore2 = C;
	    PFLOAT Ascore = FLAT(A,score, I, J);/* NOTE: Ascore will be added to bestscore2 last, unlike for bestscore at bestI,bestJ */

	    int *RijxI = &FLAT(A,Rijx, I, 0), *RijyI = &FLAT(A,Rijy, I, 0);
	    PFLOAT Ydiff = Yr[I];// Y[N+1]-Y[I];

	    if(EVERB>=2 && I==33 && J==46){
	      printf("Yid=%d(id=%lld),Xid=%d(id=%lld(),or=%d,c=%d:bestI=%d,bestJ=%d,bestI2=%d,bestJ2=%d:I=%d,J=%d:original Rijx[I][J]=%d,Rijy[I][J]=%d,Send(%0.8f,%d+%d)=%0.6f\n",
		     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation, c, bestI[c],bestJ[c],bestI2,bestJ2,I,J, RijxI[J], RijyI[J],
		     PMIN(Xr[J],Ydiff),min(M,RijxI[J])+1-J, min(N,RijyI[J])+1-I, Send(PMIN(Xr[J],Ydiff),(PFLOAT)(min(M,RijxI[J])+1-J + min(N,RijyI[J])+1-I),c));
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
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,bestI2=%d,bestJ2=%d:I=%d,J=%d:Updated Rijx[I][J]=%d,Rijy[I][J]=%d\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],bestI2,bestJ2,I,J, RijxI[J], RijyI[J]);
		fflush(stdout);
	      }
	    }

	    PFLOAT score = Send(PMIN(Xr[J], Ydiff), (PFLOAT)(min(M, RijxI[J])+1-J + min(N, RijyI[J])+1-I),c);
	    if(EVERB>=2 && I==33 && J==46){
	      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,bestI2=%d,bestJ2=%d:I=%d,J=%d:Initial Rijx[I][J]=%d,Rijy[I][J]=%d,Send(%0.8f,%d+%d)=%0.6f\n",
		     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],bestI2,bestJ2,I,J, RijxI[J], RijyI[J],
		     PMIN(Xr[J],Ydiff),min(M,RijxI[J])+1-J, min(N,RijyI[J])+1-I, Send(PMIN(Xr[J],Ydiff),(PFLOAT)(min(M,RijxI[J])+1-J + min(N,RijyI[J])+1-I),c));
	      fflush(stdout);
	    }

	    if(ENDFIX){/* update Rijx,Rijy based on Sbnd() as usual */
	      if(RijxI[J] > M) {
		double a = M+1-I-J;
		int K;
		for(K= RijyI[J]; K > I; K--){
		  PFLOAT bound = Sbnd(Xr[J],Y[K]-Y[I], K+a,c);
		  if(bound <= score)
		    break;// HERE HERE : see ENDFIX2 in refalign.cpp
		  if(EVERB>=2 && I==33 && J==46){
		    printf("    I=%d,J=%d,N=%d,M=%d:K=%d,score=%0.6f,Sbnd(Xr[J]=%0.8f,Y[K]-Y[I]=%0.8f,K+a=%0.1f,c=%d)=%0.6f:Rijy[I][J]=%d->%d\n",
			   I,J,N,M,K,score,Xr[J],Y[K]-Y[I],K+a,c,Sbnd(Xr[J],Y[K]-Y[I],K+a,c),K,K-1);
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
		  PFLOAT bound = Sbnd(X[K]-X[J], Ydiff,K+a,c);
		  if(bound <= score)
		    break;// HERE HERE : see ENDFIX2 in refalign.cpp
		  score = bound;
		}
		RijxI[J] = K;
	      }
	      if(EVERB>=2 && I==33 && J==46){
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,bestI2=%d,bestJ2=%d:I=%d,J=%d:Final Rijx[I][J]=%d,Rijy[I][J]=%d\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],bestI2,bestJ2,I,J, RijxI[J], RijyI[J]);
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
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,M=%d,N=%d:checking for repeats\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,M,N);
		fflush(stdout);
	      }
	    }

	    /* Now backtrack to extract best alignment starting at (I,J)==(bestI2,bestJ2) constrained by requiring J to be less than in the best alignment for the same I */
	    int U = K-1;/* keep (I2 == align[c].sites1[U]) <= I, and J > (J2 == align[c].sites2[T2]) to satisfy the repeat constraint
			   Also keep MinShiftKb <= Y[I2] - Y[I] + X[J] - X[J2] <= MaxShiftKb */
	    int K2 = 0;

	    int bestG = -1, bestH = -1;
	    if(DEBUG>=2) assert(I==bestI2 && J==bestJ2);

	    while(I >= 0){
	      Ilist2[K2] = I;
	      Jlist2[K2++] = J;

	      /* update A->G[I][J] & A->H[I][J] by repeating recurrance computation constrained by J being less than in the best alignment for the same I */
	      /* NOTE : this needs to be done even if DEFER_BP==0, since original A->G & A->H values would have been computed without the constraint */
	      /* NOTE2 : A->score[I][J] is NOT updated, since we will rely on Pvalue to decide if the overall alignment found is good enough */

	      if(FIX_CONSTRAINT && (FIX_CONSTRAINT>=2 || bestJ[c] - bestJ2 >= 2)){
		while(U > 0 && align[c].sites1[U] > I)
		  U--;
	      } else {
		while(U > 0 && align[c].sites1[U-1] >= I)
		  U--;
	      }

	      if(DEBUG>=2 && MinShiftKb > 0.0 && align[c].sites1[U] >= I && !(J < align[c].sites2[U])){
                #pragma omp critical 
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,M=%d,N=%d:checking for repeats\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,M,N);
		  printf("  I=%d,J=%d:U=%d,align[c].sites1[U]=%d,align[c].sites2[U]=%d\n",I,J,U,align[c].sites1[U],align[c].sites2[U]);
		  fflush(stdout);
		  assert(J < align[c].sites2[U]);
		}
	      }

	      PFLOAT score = FLAT(A,Uscore, I, J);
	      if(DEBUG>=2 && !(isfinite(score))){
                #pragma omp critical 
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,M=%d,N=%d:checking for repeats\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align->score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,M,N);
		  printf("  I=%d,J=%d:U=%d,align->sites1[U]=%d,align[c].sites2[U]=%d,A[c].Uscore[I][J]=%0.8e\n",I,J,U,align[c].sites1[U],align[c].sites2[U],score);
		  fflush(stdout);
		  assert(isfinite(score));
		}
	      }
	      bestG = localtype;

	      int Gmin = max(I-DELTA,1), Hmin = max(J-DELTA,1);
	      int T3 = U;

	      for(int G= I; --G >= Gmin;){
		FLOAT XG;/* best estimate of original alignment location on X[] corresponding to Y[G] */
		int origT3 = T3;
		while(T3 < K-1 && align[c].sites1[T3] < G)
		  T3++;
		while(T3 > 0 && align[c].sites1[T3-1] >= G)
		  T3--;
		int I3 = align[c].sites1[T3];
		int J3 = align[c].sites2[T3];

		if(DEBUG>=2 && !(T3 >= K-1 || I3 >= G)){
		  #pragma omp critical
		  {
		    printf("While Checking for repeats:RepeatRec=%d,RepeatMaxShift=%d,RepeatLogPvalueRatio=%0.6f,RepeatKbShiftRatio=%0.2f,extend=%d,Lij2=%d,Rij2=%d,M=%d(score=%0.6f,LogPV=%0.2f,numpairs=%d,repeat=%d)\n",
			   RepeatRec,RepeatMaxShift,RepeatLogPvalueRatio,RepeatKbShiftRatio,extend,align[c].Lij2,align[c].Rij2,M,align[c].score,align[c].logPV,align[c].numpairs,align[c].repeat);
		    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,maptype=%d:\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c, bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,maptype);
		    printf("K=%d,U=%d,I=%d,J=%d,G=%d:T3=%d->%d:I3=%d,J3=%d\n",
			   K,U,I,J,G,origT3, T3,I3,J3);
		    fflush(stdout);
		    assert(T3 >= K-1 || I3 >= G);
		  }
		}
			
		int T4 = T3;
		while(T4 > 0 && align[c].sites1[T4] > G)
		  T4--;
		int I4 = align[c].sites1[T4];
		
		if(DEBUG>=2 && !(T4 <= 0 || I4 <= G)){
		  #pragma omp critical
		  {
		    printf("While Checking for repeats:c=%d,RepeatRec=%d,RepeatMaxShift=%d,RepeatLogPvalueRatio=%0.6f,RepeatKbShiftRatio=%0.2f,extend=%d,Lij2=%d,Rij2=%d,M=%d,repeat=%d\n",
			   c,RepeatRec,RepeatMaxShift,RepeatLogPvalueRatio,RepeatKbShiftRatio,extend,align[c].Lij2,align[c].Rij2,M,align[c].repeat);
		    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,maptype=%d:\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c, bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,maptype);
		    printf("K=%d,U=%d,I=%d,J=%d,G=%d:T3=%d->%d,T4=%d:I3=%d,J3=%d,I4=%d\n",
			   K,U,I,J,G,origT3, T3, T4, I3,J3, I4);
		    fflush(stdout);
		    assert(T4 <= 0 || I4 <= G);
		  }
		}
		origT3 = T3;

		if(G < I4){
		  int J4 = align[c].sites2[T4];
		  if(I4  <= G)
		    XG = X[J4];
		  else
		    XG = X[J4] - (Y[I4] - Y[G]);
		} else if(I3 > I4){
		  if(DEBUG>=2) assert(I4 <= G && G <= I3);
		  if(DEBUG>=2) assert(Y[I3] > Y[I4]);
		  FLOAT Gfrac = (Y[G] - Y[I4])/(Y[I3] - Y[I4]);
		  int J4 = align[c].sites2[T4];
		  XG = X[J4] + Gfrac * (X[J3] - X[J4]);
		} else  {
		  if(DEBUG>=2) assert(I3 == I4);
		  if(DEBUG>=2) assert(T3 == T4);
		  XG = X[J3];
		}

		if(FIX_CONSTRAINT && (FIX_CONSTRAINT>=2 || bestJ[c] - bestJ2 >= 2)){
		  while(T3 > 0 && align[c].sites1[T3] > G)
		    T3--;
		} else {
		  while(T3 > 0 && align[c].sites1[T3-1] >= G)
		    T3--;
		}
		I3 = align[c].sites1[T3];
		J3 = align[c].sites2[T3];

		/* H must satistfy the constraint: H < J3 AND -MaxShiftKb <= Y[I3]-Y[G] + X[H] - XG <= -MinShiftKb */
		FLOAT XHmax = XG - Y[I3];
		FLOAT XHmin = XHmax;
		int Hmax = min(J,J3);/* H must be less than this value */
		if(I3 < G && T3 < K-1){
		  int I3 = align[c].sites1[T3+1];
		  int J3 = align[c].sites2[T3+1];
		  if(X[J3] - Y[I3] > XHmax)
		    XHmax = X[J3] - Y[I3];
		  if(X[J3] - Y[I3] < XHmin)
		    XHmin = X[J3] - Y[I3];
		  if(Hmax < J3)
		    Hmax = min(J,J3);
		}
		if(G < I3 && T3 > 0){
		  int I3 = align[c].sites1[T3-1];
		  int J3 = align[c].sites2[T3-1];
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
		PFLOAT *AscoreG = &FLAT(A, score, G, 0);
		PFLOAT a = J+I-G;
	
		for(int H= Hmin; H < Hmax; H++){
		  if(X[H] < XHmin && H < Hmax)
		    continue;
		  if(X[H] > XHmax)
		    break;
		  PFLOAT deltaX = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];// X[J] - X[H]
		  if(DEBUG>=2) assert(fabs(X[J] - X[H] - deltaX) < 0.001);
		  PFLOAT newscore = AscoreG[H] +  Sint(deltaX,deltaY,(PFLOAT)(a-H),c);
		  if(newscore > score){
		    if(EVERB)
		      printf("c=%d:I=%d,J=%d,G=%d,H=%d:T3=%d,I3=%d,J3=%d:score=%0.6f->%0.6f\n",
			     c,I,J,G,H,T3,align[c].sites1[T3],align[c].sites2[T3],score,newscore);
		    bestG = G;
		    bestH = H;
		    score = newscore;
		  }
		}
	      }
	      if(bestG < 0){/* Left End : value of bestH is undefined */
		if(DEFER_G >= 1 && score > C){
		  bestG = -1;
		  if(DEFER_G >= 2){/* check if actual value of Lijx,Lijy needs to be incremented */
		    int Lijx = FLAT(A,Lijx,I,J);
		    int Lijy = FLAT(A,Lijy,I,J);
		    PFLOAT SmIJ = Sm(J,I);
		    PFLOAT ascore = C;

		    if(Lijx <= 0){
		      ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I+1-max(1,Lijy)),c) + SmIJ;
		      ascore = max(ascore,C);
		      for(int K= Lijy; K < I; ) {
			PFLOAT bound = Sbnd(X[J],Y[I]-Y[K],(PFLOAT)(J+I-K),c) + SmIJ;
			if(bound <= ascore)
			  break;// HERE : see ENDFIX2 in refalign.cpp
			ascore = bound;
			Lijy = ++K;
		      }
		      FLAT(A,Lijy, I, J) = Lijy;
		    } else {
		      if(DEBUG>=2) assert(Lijy <= 0);
		      ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-max(1,Lijx)+I),c) + SmIJ;
		      ascore = max(ascore,C);
		      for(int K= Lijx; K < J; ){
			PFLOAT bound = Sbnd(X[J]-X[K],Y[I],(PFLOAT)(J+I-K),c) + SmIJ;
			if(bound <= ascore)
			  break;// HERE : see ENDFIX2 in refalign.cpp
			ascore = bound;
			Lijx = ++K;
		      }
		      FLAT(A,Lijx, I, J) = Lijx;
		    }
		  }
		}
		FLAT(A,G, I, J) = bestG;
		J = -1;
	      } else {
		FLAT(A,G, I, J) = bestG;
		FLAT(A,H, I, J) = bestH;
		J = bestH;
	      }
	      I = bestG;
	    }
	  
	    if(DEBUG>=2){
	      int Rijy = (bestR2 <= -2) ? bestI2 : FLAT(A,Rijy, bestI2, bestJ2);
	      int Rijx = (bestR2 <= -2) ? bestJ2 : FLAT(A,Rijx, bestI2, bestJ2);

	      if(!(min(N,Rijy)-bestI2 >= 0)){
                #pragma omp critical
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,Rijy=%d,Rijx=%d,M=%d,N=%d\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,Rijy,Rijx,M,N);
		  for(int T2 = 0; T2 < K2; T2++)
		    printf("   T=%d/%d:I=%d,J=%d\n",T2,K2,Ilist2[K2-1-T2],Jlist2[K2-1-T2]);
		  fflush(stdout);
		  assert(min(N,Rijy)-bestI2 >= 0);
		}
	      }
	    }

	    double logPV2 = alignFP(c,A,bestI2,bestJ2,bestR2,Y,X,N,M,orientation,Yid,Xid,delY,delX,bestscore2, /* (EVERB && bestJ-bestJ2==2) ? 1 : 0*/ 0, 0);
	    if(logPV2 >= max(align[c].logPV * RepeatLogPvalueRatio, minPvalue)){
	      align[c].repeat = bestJ[c] - bestJ2;
	      if(DEBUG) assert(align->repeat > 0);

	      if(VERB>=2 || EVERB){
                #pragma omp critical
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.6f,numpairs=%d:bestI2=%d,bestJ2=%d,score2=%0.4f,logPV2=%0.6f (repeat=%d confirmed)\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,bestscore2,logPV2,align[c].repeat);
		  int I2 = -1, J2= -1,G2= -1,H2= -1;
		  for(int T2 = 0; T2 < K2; T2++, G2=I2, H2=J2){
		    I2 = Ilist2[K2-1-T2];
		    J2 = Jlist2[K2-1-T2];
		    if(T2 > 0){
		      PFLOAT score = FLAT(A,score, I2, J2) -FLAT(A,score, G2, H2);
		      PFLOAT y = Y[I2] - Y[G2];
		      PFLOAT x = X[J2] - X[H2];
		      PFLOAT Bias, Pen;
		      SintDetail(x,y, J2-H2, I2-G2, J2, I2, Bias, Pen,c);
		      printf("   T=%d/%d:I=%d,J=%d (Y=%0.4f X=%0.4f):score=%0.4f(Pen=%0.4f,Bias=%0.4f):A->score=%0.4f\n",
			     T2,K2,I2,J2, y, x, score, Pen, Bias, FLAT(A,score,I2,J2));
		    } else { /* left end */
		      PFLOAT score = FLAT(A,score, I2, J2);
		      printf("   T2=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:A->score=%0.4f\n",T2,K2,I2, J2, A->Lijy[I2][J2],A->Lijx[I2][J2],score);
		    }
		  }
		  /* right end */
		  int Rijy = (I2 >= N-DELTA || J2 >= M-DELTA) ? A->Rijy[I2][J2] : N+1;
		  int Rijx = (I2 >= N-DELTA || J2 >= M-DELTA) ? A->Rijx[I2][J2] : M+1;
		  if(bestR2 <= -2){
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.4f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.4f,ChimScore=%0.6f,total score=%0.8e\n",
			   I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],ChimScore,FLAT(A,score, I2, J2)+ChimScore);
		  } else if(X[M+1]-X[J2] <= Y[N+1]-Y[I2] && Rijy <= N && Y[Rijy+1]-Y[I2] < X[M+1]-X[J2]){
		    if(DEBUG) assert(Rijx > M);
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.8f,%0.8f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],X[M+1]-X[J2],Y[Rijy+1]-Y[I2],min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			   Sbnd(X[M+1]-X[J2],Y[Rijy+1]-Y[I2],M+1-J2,Rijy+1-I2,c),FLAT(A,score, I2, J2) + Sbnd(X[M+1]-X[J2],Y[Rijy+1]-Y[I2],M+1-J2,Rijy+1-I2,c));
		  } else if(X[M+1]-X[J2] >= Y[N+1]-Y[I2] && Rijx <= M && X[Rijx+1]-X[J2] < Y[N+1]-Y[I2]){
		    if(DEBUG) assert(Rijy > N);
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],
			   X[Rijx+1]-X[J2],Y[N+1]-Y[I2],min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			   Sbnd(X[Rijx+1]-X[J2],Y[N+1]-Y[I2],Rijx+1-J2,N+1-I2,c),
			   FLAT(A,score, I2, J2) + Sbnd(X[Rijx+1]-X[J2],Y[N+1]-Y[I2],Rijx+1-J2,N+1-I2,c));
		  } else
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Send(%0.8f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],
			   min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			   Send(min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,c),
			   FLAT(A,score, I2, J2) + Send(min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,c));

		  printf("  Best alignment:\n");
		  int I = -1, J = -1, G = -1, H =  -1;
		  for(int T = 0; T < K; T++, G=I, H=J){
		    I = Ilist[K-1-T];
		    J = Jlist[K-1-T];
		    if(T > 0){/* internal interval */
		      PFLOAT score = FLAT(A,score, I, J) - FLAT(A,score, G, H);
		      PFLOAT y = Y[I] - Y[G];
		      PFLOAT x = X[J] - X[H];
		      PFLOAT Bias, Pen;
		      SintDetail(x, y, J-H, I-G, J, I, Bias, Pen,c);
		      printf("   T=%d/%d:I=%d,J=%d (Y=%0.4f X=%0.4f):score=%0.4f(Pen=%0.4f,Bias=%0.4f):A->score=%0.4f\n",
			     T,K,I,J, y, x, score, Pen, Bias, FLAT(A,score, I, J));
		    } else { /* left end */
		      PFLOAT score = FLAT(A,score, I2, J2);
		      printf("   T=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:A->score=%0.4f\n", T,K,I, J, A->Lijy[I][J],A->Lijx[I][J],score);
		    }
		  }
		  /* right end */
		  Rijy = (I >= N-DELTA || J >= M-DELTA) ? A->Rijy[I][J] : N+1;
		  Rijx = (I >= N-DELTA || J >= M-DELTA) ? A->Rijx[I][J] : M+1;

		  if(bestR[c] <= -2){
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.4f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.4f,ChimScore=%0.6f,total score=%0.8e\n",
			   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],ChimScore,FLAT(A,score, I, J)+ChimScore);
		  } else if(X[M+1]-X[J] <= Y[N+1]-Y[I] && Rijy <= N && Y[Rijy+1]-Y[I] < X[M+1]-X[J]){
		    if(DEBUG) assert(Rijx > M);
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],X[M+1]-X[J],Y[Rijy+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			   Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I,c),FLAT(A,score, I, J) + Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I,c));
		  } else if(X[M+1]-X[J] >= Y[N+1]-Y[I] && Rijx <= M && X[Rijx+1]-X[J] < Y[N+1]-Y[I]){
		    if(DEBUG) assert(Rijy > N);
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
			   X[Rijx+1]-X[J],Y[N+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			   Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I,c),
			   FLAT(A,score, I, J) + Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I,c));
		  } else
		    printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Send(%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
			   min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			   Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,c),
			   FLAT(A,score, I, J) + Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,c));
		
		  fflush(stdout);
		}
	      }
	      if(DEBUG) assert(align[c].repeat > 0);
	      break;
	    } else if(VERB>=2 || (EVERB/* HERE >=2 && bestJ[c] - bestJ2==2 */ )){
              #pragma omp critical
	      {
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.6f,numpairs=%d:bestI2=%d,bestJ2=%d,score2=%0.4f,logPV2=%0.6f (No repeat at shiftJ =%d)\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,bestscore2,logPV2,bestJ[c] - bestJ2);
		int I2 = -1, J2= -1,G2= -1,H2= -1;
		for(int T2 = 0; T2 < K2; T2++, G2=I2, H2=J2){
		  I2 = Ilist2[K2-1-T2];
		  J2 = Jlist2[K2-1-T2];
		  if(T2 > 0){
		    PFLOAT score = FLAT(A,score, I2, J2) -FLAT(A,score, G2, H2);
		    PFLOAT y = Y[I2] - Y[G2];
		    PFLOAT x = X[J2] - X[H2];
		    PFLOAT Bias, Pen;
		    SintDetail(x,y, J2-H2, I2-G2, J2, I2, Bias, Pen,c);
		    printf("   T=%d/%d:I=%d,J=%d (Y=%0.4f X=%0.4f):score=%0.4f(Pen=%0.4f,Bias=%0.4f):A.score=%0.4f\n",
			   T2,K2,I2,J2, y, x, score, Pen, Bias, FLAT(A,score,I2,J2));
		  } else { /* left end */
		    PFLOAT score = FLAT(A,score, I2, J2);
		    printf("   T2=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:A.score=%0.4f\n",T2,K2,I2, J2, A->Lijy[I2][J2],A->Lijx[I2][J2],score);
		  }
		}
		/* right end */
		int Rijy = (I2 >= N-DELTA || J2 >= M-DELTA) ? A->Rijy[I2][J2] : N+1;
		int Rijx = (I2 >= N-DELTA || J2 >= M-DELTA) ? A->Rijx[I2][J2] : M+1;
		if(bestR2 <= -2){
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.4f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.4f,ChimScore=%0.6f,total score=%0.8e\n",
			 I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],ChimScore,FLAT(A,score, I2, J2)+ChimScore);
		} else if(X[M+1]-X[J2] <= Y[N+1]-Y[I2] && Rijy <= N && Y[Rijy+1]-Y[I2] < X[M+1]-X[J2]){
		  if(DEBUG) assert(Rijx > M);
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],X[M+1]-X[J2],Y[Rijy+1]-Y[I2],min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			 Sbnd(X[M+1]-X[J2],Y[Rijy+1]-Y[I2],M+1-J2,Rijy+1-I2,c),FLAT(A,score, I2, J2) + Sbnd(X[M+1]-X[J2],Y[Rijy+1]-Y[I2],M+1-J2,Rijy+1-I2,c));
		} else if(X[M+1]-X[J2] >= Y[N+1]-Y[I2] && Rijx <= M && X[Rijx+1]-X[J2] < Y[N+1]-Y[I2]){
		  if(DEBUG) assert(Rijy > N);
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],
			 X[Rijx+1]-X[J2],Y[N+1]-Y[I2],min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			 Sbnd(X[Rijx+1]-X[J2],Y[N+1]-Y[I2],Rijx+1-J2,N+1-I2,c),
			 FLAT(A,score, I2, J2) + Sbnd(X[Rijx+1]-X[J2],Y[N+1]-Y[I2],Rijx+1-J2,N+1-I2,c));
		} else
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Send(%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I2,J2,N,M,Rijy,Rijx,X[M+1]-X[J2],Y[N+1]-Y[I2],Y[I2+1]-Y[I2],
			 min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,
			 Send(min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,c),
			 FLAT(A,score, I2, J2) + Send(min(X[M+1]-X[J2],Y[N+1]-Y[I2]),min(M,Rijx)+1-J2,min(N,Rijy)+1-I2,c));

		printf("  Best alignment:\n");
		int I = -1, J = -1, G = -1, H =  -1;
		for(int T = 0; T < K; T++, G=I, H=J){
		  I = Ilist[K-1-T];
		  J = Jlist[K-1-T];
		  if(T > 0){/* internal interval */
		    PFLOAT score = FLAT(A,score, I, J) - FLAT(A,score, G, H);
		    PFLOAT y = Y[I] - Y[G];
		    PFLOAT x = X[J] - X[H];
		    PFLOAT Bias, Pen;
		    SintDetail(x, y, J-H, I-G, J, I, Bias, Pen,c);
		    printf("   T=%d/%d:I=%d,J=%d (Y=%0.4f X=%0.4f):score=%0.4f(Pen=%0.4f,Bias=%0.4f):A.score=%0.4f\n",
			   T,K,I,J, y, x, score, Pen, Bias, FLAT(A,score, I, J));
		  } else { /* left end */
		    PFLOAT score = FLAT(A,score, I2, J2);
		    printf("   T=%d/%d:I=%d,J=%d:Lijy=%d,Lijx=%d:A.score=%0.4f\n", T,K,I, J, A->Lijy[I][J],A->Lijx[I][J],score);
		  }
		}
		/* right end */
		Rijy = (I >= N-DELTA || J >= M-DELTA) ? A->Rijy[I][J] : N+1;
		Rijx = (I >= N-DELTA || J >= M-DELTA) ? A->Rijx[I][J] : M+1;

		if(bestR[c] <= -2){
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.4f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.4f,ChimScore=%0.6f,total score=%0.8e\n",
			 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],ChimScore,FLAT(A,score, I, J)+ChimScore);
		} else if(X[M+1]-X[J] <= Y[N+1]-Y[I] && Rijy <= N && Y[Rijy+1]-Y[I] < X[M+1]-X[J]){
		  if(DEBUG) assert(Rijx > M);
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],X[M+1]-X[J],Y[Rijy+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			 Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I,c), FLAT(A,score, I, J) + Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I,c));
		} else if(X[M+1]-X[J] >= Y[N+1]-Y[I] && Rijx <= M && X[Rijx+1]-X[J] < Y[N+1]-Y[I]){
		  if(DEBUG) assert(Rijy > N);
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
			 X[Rijx+1]-X[J],Y[N+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			 Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I,c),
			 FLAT(A,score, I, J) + Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I,c));
		} else
		  printf("    Right End:I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Send(%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
			 I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
			 min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,
			 Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,c),
			 FLAT(A,score, I, J) + Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,c));
		
		fflush(stdout);
	      }
	    }
	  } // 	for(bestJ2 = max(1,bestJ - RepeatMaxShift); bestJ2 < bestJ; bestJ2++)

	  if(DEBUG && !align[c].repeat) assert(bestJ2 == bestJ[c]);

	  for(bestI2 = max(1, align[c].sites1[max(0,K-1-RepeatMaxShift)]); !align[c].repeat && bestI2 < bestI[c]; bestI2++){/* back track best alignment from (bestI2,bestJ2) and compute its Pvalue */
	    // WAS	  for(bestI2 = max(1,bestI[c] - RepeatMaxShift); !align[c].repeat && bestI2 < bestI[c]; bestI2++){/* back track best alignment from (bestI2,bestJ2) and compute its Pvalue */
	    double MinShiftKb = (Y[bestI[c]] - Y[bestI2]) * RepeatKbShiftRatio; /* minimum offset shift along entire alignment */
	    double MaxShiftKb = (Y[bestI[c]] - Y[bestI2]) * RepeatKbShiftRatio2; /* maximum offset shift along entire alignment (to avid creating new outliers) */

	    int I = bestI2;
	    int J = bestJ2;

	    /* First compute bestscore2, bestR2 with rightmost alignment at (I,J) */
	    int bestR2 = localtype;
	    PFLOAT bestscore2 = C;
	    PFLOAT Ascore = FLAT(A,score, I, J);/* NOTE: Ascore will be added to bestscore2 last, unlike for bestscore,bestI,bestJ */

	    int * RijxI = &FLAT(A,Rijx, I, 0), *RijyI = &FLAT(A,Rijy, I, 0);
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

	    PFLOAT score = Send(PMIN(Xr[J], Ydiff), (PFLOAT)(min(M, RijxI[J])+1-J + min(N, RijyI[J])+1-I),c);

	    if(ENDFIX){
	      if(RijxI[J] > M) {
		double a = M+1-I-J;
		int K;
		for(K= RijyI[J]; K > I; K--){
		  PFLOAT bound = Sbnd(Xr[J],Y[K]-Y[I], K+a,c);
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
		  PFLOAT bound = Sbnd(X[K]-X[J], Ydiff,K+a,c);
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
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,M=%d,N=%d:checking for repeats\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,M,N);
		fflush(stdout);
	      }
	    }

	    /* Now backtrack to extract best alignment starting at (I,J)==(bestI2,bestJ2) constrained by requiring I to be less than in the best alignment for the same J */
	    int U = K-1;/* keep (I2 == align[c].sites1[U]) <= I, and J > (J2 == align[c].sites2[U]) to satisfy the repeat constraint
			   Also keep MinShiftKb <= Y[I2] - Yc[I] + X[J] - X[J2] <= MaxShiftKb */
	    int K2 = 0;

	    int bestG = -1, bestH = -1;
	    while(I >= 0){
	      Ilist2[K2] = I;
	      Jlist2[K2++] = J;

	      /* update A->G[I][J] & A->H[I][J] by repeating recurrance computation constained by requiring I to be less than in the best alignment for the same J */
	      /* NOTE : this needs to be done even if DEFER_BP==0, since original A->G & A->H values would have been computed without the constraint */	    

	      if(FIX_CONSTRAINT && (FIX_CONSTRAINT>=2 || bestI[c] - bestI2 >= 2)){
		while(U > 0 && align[c].sites1[U] > I)
		  U--;
	      } else {
		while(U > 0 && align[c].sites1[U-1] >= I)
		  U--;
	      }

	      PFLOAT score = FLAT(A,Uscore, I, J);
	      bestG = localtype;

	      int Gmin = max(I-DELTA,1), Hmin = max(J-DELTA,1);
	      int T3 = U;

	      for(int G= I; --G >= Gmin;){
		FLOAT XG;/* best estimate of original alignment location on X[] corresponding to Y[G] */
		int origT3 = T3;
		while(T3 < K-1 && align[c].sites1[T3] < G)
		  T3++;
		while(T3 > 0 && align[c].sites1[T3-1] >= G)
		  T3--;
		int I3 = align[c].sites1[T3];
		int J3 = align[c].sites2[T3];

		if(DEBUG>=2 && !(T3 >= K-1 || I3 >= G)){
		  #pragma omp critical
		  {
		    printf("While Checking for repeats:RepeatRec=%d,RepeatMaxShift=%d,RepeatLogPvalueRatio=%0.6f,RepeatKbShiftRatio=%0.2f,extend=%d,Lij2=%d,Rij2=%d,M=%d(repeat=%d)\n",
			   RepeatRec,RepeatMaxShift,RepeatLogPvalueRatio,RepeatKbShiftRatio,extend,align[c].Lij2,align[c].Rij2,M,align[c].repeat);
		    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,maptype=%d:\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c, bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,maptype);
		    printf("K=%d,U=%d,I=%d,J=%d,G=%d:T3=%d->%d:I3=%d,J3=%d\n",
			   K,U,I,J,G,origT3, T3,I3,J3);
		    fflush(stdout);
		    assert(T3 >= K-1 || I3 >= G);
		  }
		}
			
		int T4 = T3;
		while(T4 > 0 && align[c].sites1[T4] > G)
		  T4--;
		int I4 = align[c].sites1[T4];
		
		if(DEBUG>=2 && !(T4 <= 0 || I4 <= G)){
		  #pragma omp critical
		  {
		    printf("While Checking for repeats:RepeatRec=%d,RepeatMaxShift=%d,RepeatLogPvalueRatio=%0.6f,RepeatKbShiftRatio=%0.2f,extend=%d,Lij2=%d,Rij2=%d,M=%d(repeat=%d)\n",
			   RepeatRec,RepeatMaxShift,RepeatLogPvalueRatio,RepeatKbShiftRatio,extend,align[c].Lij2,align[c].Rij2,M,align[c].repeat);
		    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,maptype=%d:\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c, bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,maptype);
		    printf("K=%d,U=%d,I=%d,J=%d,G=%d:T3=%d->%d,T4=%d:I3=%d,J3=%d,I4=%d\n",
			   K,U,I,J,G,origT3, T3, T4, I3,J3, I4);
		    fflush(stdout);
		    assert(T4 <= 0 || I4 <= G);
		  }
		}
		origT3 = T3;

		if(G < I4){
		  int J4 = align[c].sites2[T4];
		  if(I4  <= G)
		    XG = X[J4];
		  else
		    XG = X[J4] - (Y[I4] - Y[G]);
		} else if(I3 > I4){
		  if(DEBUG>=2) assert(I4 <= G && G <= I3);
		  if(DEBUG>=2) assert(Y[I3] > Y[I4]);
		  FLOAT Gfrac = (Y[G] - Y[I4])/(Y[I3] - Y[I4]);
		  int J4 = align[c].sites2[T4];
		  XG = X[J4] + Gfrac * (X[J3] - X[J4]);
		} else  {
		  if(DEBUG>=2) assert(I3 == I4);
		  if(DEBUG>=2) assert(T3 == T4);
		  XG = X[J3];
		}

		if(FIX_CONSTRAINT && (FIX_CONSTRAINT>=2 || bestI[c] - bestI2 >= 2)){
		  while(T3 > 0 && align[c].sites1[T3] > G)
		    T3--;
		} else {
		  while(T3 > 0 && align[c].sites1[T3-1] >= G)
		    T3--;
		}
		I3 = align[c].sites1[T3];
		J3 = align[c].sites2[T3];

		/* H must satisfy the constraint: H > J3 AND MinShiftKb <= Y[I3] - Y[G] + X[H] - XG <= MaxShiftKb */
		FLOAT XHmin = XG - Y[I3];
		FLOAT XHmax = XHmin;
		int HminG = max(Hmin,J3+1);
		if(I3 < G && T3 < K-1){
		  int I3 = align[c].sites1[T3+1];
		  int J3 = align[c].sites2[T3+1];
		  if(X[J3] - Y[I3] < XHmin)
		    XHmin = X[J3] - Y[I3];
		  if(X[J3] - Y[I3] > XHmax)
		    XHmax = X[J3] - Y[I3];
		  if(HminG > J3+1)
		    HminG = max(Hmin,J3+1);
		}
		if(G < I3 && T3 > 0){
		  int I3 = align[c].sites1[T3-1];
		  int J3 = align[c].sites2[T3-1];
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
		PFLOAT *AscoreG = &FLAT(A,score, G, 0);
		PFLOAT a = J+I-G;

		for(int H= J; --H >= HminG;){
		  if(X[H] > XHmax && H > HminG)
		    continue;
		  if(X[H] < XHmin)
		    break;
		  PFLOAT deltaX = ((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)) ? FLAT_DELX(delX, J-H, J, M) : delX[J][H-J+DELTA];// X[J] - X[H]
		  if(DEBUG>=2) assert(fabs(X[J] - X[H] - deltaX) < 0.001);
		  PFLOAT newscore = AscoreG[H] +  Sint(deltaX,deltaY,(PFLOAT)(a-H),c);
		  if(newscore > score){
		    score = newscore;
		    bestG = G;
		    bestH = H;
		  }
		}
	      }
	      if(bestG < 0){/* Left End : value of bestH is undefined */
		if(DEFER_G >= 1 && score > C){
		  bestG = -1;
		  if(DEFER_G >= 2){/* check if actual value of Lijx,Lijy needs to be incremented */
		    int Lijx = FLAT(A,Lijx,I,J);
		    int Lijy = FLAT(A,Lijy,I,J);
		    PFLOAT SmIJ = Sm(J,I);
		    PFLOAT ascore = C;

		    if(Lijx <= 0){
		      ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+I+1-max(1,Lijy)),c) + SmIJ;
		      ascore = max(ascore,C);
		      for(int K= Lijy; K < I; ) {
			PFLOAT bound = Sbnd(X[J],Y[I]-Y[K],(PFLOAT)(J+I-K),c) + SmIJ;
			if(bound <= ascore)
			  break;// HERE : see ENDFIX2 in refalign.cpp
			ascore = bound;
			Lijy = ++K;
		      }
		      FLAT(A,Lijy, I, J) = Lijy;
		    } else {
		      if(DEBUG>=2) assert(Lijy <= 0);
		      ascore = Send(min(Xf[J],Yf[I]),(PFLOAT)(J+1-max(1,Lijx)+I),c) + SmIJ;
		      ascore = max(ascore,C);
		      for(int K= Lijx; K < J; ){
			PFLOAT bound = Sbnd(X[J]-X[K],Y[I],(PFLOAT)(J+I-K),c) + SmIJ;
			if(bound <= ascore)
			  break;// HERE : see ENDFIX2 in refalign.cpp
			ascore = bound;
			Lijx = ++K;
		      }
		      FLAT(A,Lijx, I, J) = Lijx;
		    }
		  }
		}
		FLAT(A,G, I, J) = bestG;
		J = -1;
	      } else {
		FLAT(A,G, I, J) = bestG;
		FLAT(A,H, I, J) = bestH;
		J = bestH;
	      }
	      I = bestG;

	      if(VERB && (VERB>=2 || EVERB>=2)){
                #pragma omp critical
		{
		  printf("  K2=%d:I=Ilist2[K2]=%d,J=Jlist[K2]=%d,A->G[I][J]=%d,A->H[I][J]=%d,score=%0.6f\n",
			 K2-1,Ilist2[K2-1],Jlist[K2-1],A->G[Ilist2[K2-1]][Jlist2[K2-1]],A->H[Ilist2[K2-1]][Jlist2[K2-1]],score);
		  fflush(stdout);
		}
	      }
	    }

	    if(DEBUG>=2){
	      int Rijy = (bestR2 <= -2) ? bestI2 : FLAT(A,Rijy, bestI2, bestJ2);
	      int Rijx = (bestR2 <= -2) ? bestJ2 : FLAT(A,Rijx, bestI2, bestJ2);

	      if(!(min(N,Rijy)-bestI2 >= 0)){
                #pragma omp critical
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,Rijy=%d,Rijx=%d,M=%d,N=%d\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,Rijy,Rijx,M,N);
		  for(int T2 = 0; T2 < K2; T2++)
		    printf("   T=%d/%d:I=%d,J=%d\n",T2,K2,Ilist2[K2-1-T2],Jlist2[K2-1-T2]);
		  fflush(stdout);
		  assert(min(N,Rijy)-bestI2 >= 0);
		}
	      }
	    }

	    double logPV2;
	    logPV2 = alignFP(c,A,bestI2,bestJ2,bestR2,Y,X,N,M,orientation,Yid,Xid,delY,delX,bestscore2,/* (Yid==285 && Xid==22611 && orientation==0) ? 1 : 0 */ 0, 0);

	    if(logPV2 >= max(align[c].logPV * RepeatLogPvalueRatio, minPvalue)){
	      align[c].repeat = bestI[c] - bestI2;
	      if(DEBUG) assert(align[c].repeat > 0);
	      if(VERB>=2 || EVERB){
                #pragma omp critical
		{
		  printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.6f,numpairs=%d:bestI2=%d,bestJ2=%d,score2=%0.4f,logPV2=%0.6f (repeat=%d confirmed)\n",
			 Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,bestscore2,logPV2,align[c].repeat);
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
	    } else if(VERB>=2 || EVERB/* HERE >=2 */){
              #pragma omp critical
	      {
		printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:bestI2=%d,bestJ2=%d,score2=%0.4f,logPV2=%0.2f (No repeat at shiftI =%d)\n",
		       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,bestI2,bestJ2,bestscore2,logPV2, bestI[c] - bestI2);
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

	  if(VERB>=2 && align[c].repeat > 0){
	    #pragma omp critical
	    {
	      repeatcnt++;

	      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d,c=%d:bestI=%d,bestJ=%d,score=%0.4f,logPV=%0.2f,numpairs=%d:repeat detected(repeat=%d), repeatcnt=%d\n",
		     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,c,bestI[c],bestJ[c],align[c].score,align[c].logPV,align[c].numpairs,align[c].repeat, repeatcnt);
	      fflush(stdout);
	    }
	  } else if(DEBUG/* HERE HERE >=2 */ && align[c].repeat > 0){

	    #pragma omp atomic
	    repeatcnt++;

	  }

	}// if(RepeatRec ... )

      } // if(bestscore > origscore)

      bestIL[c] = I = Ilist[K-1];// align[c].sites1[0];
      bestJL[c] = J = Jlist[K-1];// align[c].sites2[0];
    } // c = 0 .. colors-1

    /* compute combined score(s) : only score,orientation,logPV,numpairs are valid for align[-1] */
    align[-1].score = bestscore;
    align[-1].orientation = orientation;
    align[-1].logPV = align[0].logPV + align[1].logPV;/* excludes misalignment information */
    align[-1].numpairs = align[0].numpairs + align[1].numpairs;
    align[-1].repeat = ((align[0].numpairs < MINSITES || align[0].repeat) && (align[1].numpairs < MINSITES || align[1].repeat)) ? 1 : 0;

    if(VERB >= (EVERB?1:2) && align->numpairs <= 0){
      printf("Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d:score=%0.4f,logPV=%0.2f,pairs=%d+%d,numpairs=%d,repeat=%d\n",
	     Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,orientation,align[-1].score,align[-1].logPV,KK[0],KK[1],align[-1].numpairs,align[-1].repeat);
      fflush(stdout);
    }

    if(VERB >= ((LVERB||PVERB||EVERB) ? 1 : 2) && KK[0]+KK[1] > 0 && (VERB>=2  || (align[-1].score > ScoreThreshold && align[-1].logPV > LogPvThreshold && align[-1].numpairs >= MINSITES) || EVERB)){
     #pragma omp critical
     {
      printf("Yid=%d(id=%lld,left=%d,%d,right=%d,%d,N=%d,%d),Xid=%d(id=%lld,left=%d,%d,right=%d,%d,M=%d,%d),or=%d:score=%0.4f,logPV=%0.4f,aligned sites=%d,repeat=%d:\n",
	     Yid,YYmap[Yid]->id,YYmap[Yid]->origmap ? YYmap[Yid]->left[0] : 0, YYmap[Yid]->origmap ? YYmap[Yid]->left[1] : 0, 
	     YYmap[Yid]->origmap ? YYmap[Yid]->right[0] : NN[0]+1, YYmap[Yid]->origmap ? YYmap[Yid]->right[1] : NN[1]+1,
	     YYmap[Yid]->origmap ? YYmap[Yid]->origmap->numsite[0] : NN[0], YYmap[Yid]->origmap ? YYmap[Yid]->origmap->numsite[1] : NN[1],
	     Xid,XXmap[Xid]->id,XXmap[Xid]->origmap ? XXmap[Xid]->left[0] : 0, XXmap[Xid]->origmap ? XXmap[Xid]->left[1] : 0,
	     XXmap[Xid]->origmap ? XXmap[Xid]->right[0] : MM[0]+1,	     XXmap[Xid]->origmap ? XXmap[Xid]->right[1] : MM[1]+1,
	     XXmap[Xid]->origmap ? XXmap[Xid]->origmap->numsite[0] : MM[0],	     XXmap[Xid]->origmap ? XXmap[Xid]->origmap->numsite[1] : MM[1],
	     orientation, bestscore, align[-1].logPV, KK[0]+KK[1], align[-1].repeat);

      if(!(VERB >= (EVERB ? 1 : 2)) || bestscore > origscore){
	FLOAT y0 = YY[0][bestI[0]];
	FLOAT y1 = YY[1][bestI[1]];
	FLOAT x0 = XX[0][bestJ[0]];
	FLOAT x1 = XX[1][bestJ[1]];
	FLOAT len;
	if(COLOR_SHIFT)
	  len = max(x1,x0) + max(y1,y0);
	
	printf(" SMA(%0.4f,%0.4f,%0.4f,%d)=%0.6f:MA_mean=%0.4f,SF[2]=%0.8e,SD[2]=%0.8e\n",
	       x1-x0,y1-y0, COLOR_SHIFT ? len : 0.0, orientation,SMA(x1-x0,y1-y0,COLOR_SHIFT ? len : 0.0, orientation), MA_mean,SF[2],SD[2]);

	for(int c = 0; c < colors; c++){
	  FLOAT *Y = YY[c];
	  FLOAT *X = XX[c];
	  int N = NN[c];
	  int M = MM[c];
	  AL *A = &AA[c];
	  
	  PFLOAT **delY = delYY[c];
	  PFLOAT **delX = delXX[c];
	  
	  double offset = Y[bestI[c]] - X[bestJ[c]];

	  int repeat,Lend,Rend;
	  repeat = align[c].repeat;
	  Lend = align[c].Lend;
	  Rend = align[c].Rend;

	  double logPV = alignFP(c,A,bestI[c],bestJ[c],bestR[c],Y,X,N,M,orientation,Yid,Xid,delY,delX,bestscore, 1, 1);

	  printf("  c=%d:score=%0.4f,logPV=%0.4f,pairs=%d,repeat=%d,offset=%0.4f,Lend=%d,Rend=%d:\n",
		 c, align[c].score,logPV,KK[c],repeat,offset,Lend,Rend);

	  int I= -1,J= -1,G= -1, H = -1;
	  for(int T=0; T < KK[c]; T++, G=I, H=J){
	    I = align[c].sites1[T];
	    J = align[c].sites2[T];

	    printf("T=%2d:I=%2d,J=%2d:",T,I,J);
	    if(T <= 0){ /* Left End: display Send(min(X[J],Y[I]),J+1-max(1,Lijx),I+1-max(1,Lijy),c)+Sm(J,I), A[I,J] */
	      if(Lend <= -2){
		printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I]=%6.3f,ChimScore=%0.4f,Sm(J,I)=%0.4f,A[I,J]=%0.6f,p->G=%d,p->H=%d\n",
		       A->Lijy[I][J],A->Lijx[I][J],X[J],Y[I],ChimScore,Sm(J,I),FLAT(A,score, I, J),A->G[I][J],A->H[I][J]);
	      } else if(X[J] <= Y[I] && A->Lijy[I][J] > 0 && Y[I]-Y[A->Lijy[I][J]-1] < X[J]){
		printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I]=%6.3f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.4f,Sm(J,I)=%0.4f,A[I,J]=%0.6f,p->G=%d,p->H=%d\n",
		       A->Lijy[I][J],A->Lijx[I][J],X[J],Y[I],X[J],Y[I]-Y[A->Lijy[I][J]-1],J,I+1-A->Lijy[I][J],
		       Sbnd(X[J],Y[I]-Y[A->Lijy[I][J]-1],(PFLOAT)(J + I+1-A->Lijy[I][J]),c),Sm(J,I),FLAT(A,score, I, J),A->G[I][J],A->H[I][J]);
		if(DEBUG) assert(A->Lijx[I][J] <= 0);
	      } else if(X[J] >= Y[I] && A->Lijx[I][J] > 0 && X[J]-X[A->Lijx[I][J]-1] < Y[I]){
		printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I]=%6.3f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.4f,Sm(J,I)=%0.4f,A[I,J]=%0.6f\n",
		       A->Lijy[I][J],A->Lijx[I][J],X[J],Y[I],X[J]-X[A->Lijx[I][J]-1],Y[I],J+1-A->Lijx[I][J],I,
		       Sbnd(X[J]-X[A->Lijx[I][J]-1],Y[I],(PFLOAT)(J+1-A->Lijx[I][J] + I),c),Sm(J,I),FLAT(A,score, I, J));
		if(DEBUG) assert(A->Lijy[I][J] <= 0);
	      } else
		printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I]=%6.3f,Send(%0.3f,%d,%d)=%0.4f,Sm(J,I)=%0.4f,A[I,J]=%0.6f\n",
		       A->Lijy[I][J],A->Lijx[I][J],X[J],Y[I],min(X[J],Y[I]),J+1-max(1,A->Lijx[I][J]),I+1-max(1,A->Lijy[I][J]),
		       Send(min(X[J],Y[I]),(PFLOAT)(J+1-max(1,A->Lijx[I][J]) + I+1- max(1,A->Lijy[I][J])),c),Sm(J,I),FLAT(A,score, I, J));
	    } else {/* internal interval : display Sint(X[J]-X[H],Y[I]-Y[G],J-H,I-G,c) + Sm(J,I), A[I,J] */
	      if(DEBUG) assert(G>=0 && H>=0);
	      PFLOAT x = X[J]-X[H];
	      PFLOAT y = Y[I]-Y[G];
	      PFLOAT sint = Sint(x,y,J-H,I-G,c);
	      PFLOAT Pen,Bias;
	      SintDetail(x,y,J-H,I-G,J,I,Bias,Pen,c);
	      double myvar = varF[c] + var[c] * (x+y);
	      if(QUADRATIC_VARIANCE)
		myvar += varR[c] * (x*x+y*y);
	      printf("G=%2d,H=%2d,X[J]=%6.3f,Y[I]=%6.3f,X[J]-X[H]=%0.3f,Y[I]-Y[G]=%0.3f,norm=%0.4f,Sint(X,Y,%d,%d)=%0.6f(Pen=%0.6f,Bias=%0.6f),A[I,J]=%0.6f\n",
		     G,H,X[J],Y[I],x,y,(x-y)*(x-y)/myvar,J-H,I-G,sint,Pen,Bias,FLAT(A,score, I, J));
	      /*	    printf("   varF2=%0.6f,var2=%0.6f, MatchScore2=%0.6f,MissPenalty=%0.6f,resKB2=%0.3f,XYlen=%0.3f,Frate=%0.6f,m=%d,n=%d\n",
			    varF2,var2,MatchScore2,MissPenalty,resKB2,max(0.0,X[J]-X[H]-resKB2)+max(0.0,Y[I]-Y[G]-resKB2),Frate,J-H,I-G);
			    printf("   varF=%0.6f,var=%0.6f, MatchScore=%0.6f,MissPenalty=%0.6f,resKB2=%0.3f,XYlen=%0.3f,Frate=%0.6f\n",
			    varF,var,MatchScore,MissPenalty,resKB2,max(0.0,X[J]-X[H]-resKB2)+max(0.0,Y[I]-Y[G]-resKB2),Frate);*/
	    }
	  }
	  /* Right End: display Rij,Send(min(X[M+1]-X[J],Y[N+1]-Y[I],c),min(M,Rijx)+1-J,min(N,Rijy)+1-I) */
	  int Rijy = (I>= N-DELTA || J >= M-DELTA) ? A->Rijy[I][J] : N+1;
	  int Rijx = (I>= N-DELTA || J >= M-DELTA) ? A->Rijx[I][J] : M+1;
	  if(Rend <= -2){
	    printf("     I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,ChimScore=%0.6f,total score=%0.8e\n",
		   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I], ChimScore,FLAT(A,score, I, J)+ChimScore);
	  } else if(X[M+1]-X[J] <= Y[N+1]-Y[I] && Rijy <= N && Y[Rijy+1]-Y[I] < X[M+1]-X[J]){
	    if(DEBUG) assert(Rijx > M);
	    printf("     I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
		   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I], X[M+1]-X[J],Y[Rijy+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
		   Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I,c), FLAT(A,score, I, J) + Sbnd(X[M+1]-X[J],Y[Rijy+1]-Y[I],M+1-J,Rijy+1-I,c));
	  } else if(X[M+1]-X[J] >= Y[N+1]-Y[I] && Rijx <= M && X[Rijx+1]-X[J] < Y[N+1]-Y[I]){
	    if(DEBUG) assert(Rijy > N);
	    printf("     I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
		   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
		   X[Rijx+1]-X[J],Y[N+1]-Y[I],min(M,Rijx)+1-J,min(N,Rijy)+1-I,
		   Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I,c),
		   FLAT(A,score, I, J) + Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I],Rijx+1-J,N+1-I,c));
	  } else
	    printf("     I=%d,J=%d,N=%d,M=%d:Rijy=%d,Rijx=%d,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I]=%0.3f,Y[I+1]-Y[I]=%0.6f,Send(%0.3f,%d,%d)=%0.6f,total score=%0.8e\n",
		   I,J,N,M,Rijy,Rijx,X[M+1]-X[J],Y[N+1]-Y[I],Y[I+1]-Y[I],
		   min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,
		   Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,c),
		   FLAT(A,score, I, J) + Send(min(X[M+1]-X[J],Y[N+1]-Y[I]),min(M,Rijx)+1-J,min(N,Rijy)+1-I,c));

	  fflush(stdout);
	}
      } // bestscore > origscore
     } // omp critical
    } // if(VERB ... ) 

    if(!(VERB >= (EVERB ? 1 : 2)) || bestscore > origscore){

      if(PairSplit){
	printf("-pairsplit not yet implemented for 2 colors\n");
	exit(1);
      }

#if 0  // HERE HERE : PairSplit conversion to 2-colors : Require both colors to agree on split location ??? Use intersection of split interval ???

      if(PAIRSPLIT_EXTEND && PairSplit){
	align->refstart = align->sites1[0];
	align->start = align->sites2[0];
	align->refend = align->sites1[align->numpairs-1];
	align->end = align->sites2[align->numpairs-1];
	if(DEBUG>=2) assert(align->refstart <= align->refend);
      }

      if(Psplit > 0.0 && PairSplit){/* NOTE : it would be more efficient for small indels to do this splitting just before output of the alignments into the XMAP (see -indel -pairsplitMinOutlier) */
	/* check if alignment has an internal region with outscore < -log(1/Psplit) : if so limit alignment to larger 
	   +ve scoring alignment region + end penalty of ChimScore (other region will be found after chimeric split, if large enough) .
	   Note : With PairSplit an internal Outlier Penalty will not work as well as a seperate parameter to specify a maximum internal -ve score (ie score < - log(1/Psplit) */
	int U = align->numpairs;
	FLOAT *pscore = align->outscore;
	FLOAT *iscore = align->iscore;
	double score = 0.0;/* cumulative outscore */
	int iLWM = -1;/* most recent proposed alignment index break (-1 = left end, U = right end)) :
			 iLWM is either left end (-1) or at the right end of an internal region with outscore < -log(1/Psplit) (AND site interval >= PairsplitMinOutlier) */
	double LWM = 0.0;/* sum(pscore[0..iLWM]) */
	int iHWM = -1;/* highest values of sum(pscore[0..iHWM]) with iLWM <= iHWM */
	double HWM = 0.0;/* sum(pscore[0..iHWM]) */
	double best = 0.0;/* best score interval = sum(pscore[bestLeft + 1 .. bestRight]) */
	int bestLeft= -1,bestRight= -1;/* range of aligned site index with best score interval */
	double threshold = log(1.0/Psplit);/* If -outlierBC, this threshold is increased by |OutlierPenaltyBC[m]+OutlierPenaltyBC[n]| */
	if(LVERB && align->logPV > LogPvThreshold && align->score > ScoreThreshold && align->numpairs >= MINSITES){
	  printf("Looking for splits:numpairs=%d:PairsplitMinOutlier=%d,outlierBC=%d\n",U,PairsplitMinOutlier,outlierBC);
	  fflush(stdout);
	}
	for(int T = 0; T <= U;T++){
	  /* PRE: -1 <= iLWM <= iHWM <= T-1 :
	     LWM = sum(pscore[0..iLWM])
	     HWM = sum(pscore[0..iHWM])
	     score = sum(pscore[0..T-1])
	     NOTE : If iLWM >= 0 : iLWM is at the right end of the last (rightmost < T) internal region with outscore < -log(1/Psplit) (AND site interval >= PairsplitMinOutlier) OR left end if none so far
	  */

	  double inc = pscore[T];
	  if(inc < 0.0 && PairsplitMinOutlier && 0 < T && T < U){
	    int IL = align->sites1[T-1], IR = align->sites1[T];
	    int JL = align->sites2[T-1], JR = align->sites2[T];
	    if(IR - IL < PairsplitMinOutlier && JR - JL < PairsplitMinOutlier)
	      inc = iscore[T] + 1e-3;
	  }
	
	  if(DEBUG>=2) assert(isfinite(pscore[T]));

	  if(VERB>=3 || (LVERB && align->logPV > LogPvThreshold && align->score > ScoreThreshold && align->numpairs >= MINSITES)){
	    printf("T=%d:iLWM=%d,iHWM=%d:LWM=%0.3f,HWM=%0.3f,score=%0.3f,inc=%0.3f:best=%0.3f(%d..%d)\n",T,iLWM,iHWM,LWM,HWM,score,inc,best,bestLeft,bestRight);
	    fflush(stdout);
	  }

	  if(inc > 0.0 && iHWM >= 0){      /* check if HWM - score > threshold : If so, schedule break in alignment from iHWM to T-1 */
	    double nthreshold = threshold;
	    if(outlierBC) {
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

	if(VERB>=3 || (LVERB && align->logPV > LogPvThreshold && align->score > ScoreThreshold && align->numpairs >= MINSITES)){
	  printf("  Final: iLWM=%d,iHWM=%d:LWM=%0.3f,HWM=%0.3f,score=%0.3f:best=%0.3f(%d..%d)U=%d\n",iLWM,iHWM,LWM,HWM,score,best,bestLeft,bestRight,U);
	  fflush(stdout);
	}

	/* update alignment to only include the best local alignment interval bestLeft .. bestRight */
	if(align->numpairs >= 2){
	  if(DEBUG && !(bestLeft >= (align->Lend <= -2 ? 0 : -1))){
	    printf("While localizing alignment:Yid=%d(id=%lld,left=%d,right=%d,N=%d),Xid=%d(id=%lld,left=%d,right=%d,M=%d):align->score=%0.6f -> %0.6f,bestLeft=%d(I=%d,J=%d),bestRight=%d(I=%d,J=%d),U=%d\n",
		   Yid,YYmap[Yid]->id,YYmap[Yid]->origmap ? YYmap[Yid]->left[0] : 0, YYmap[Yid]->origmap ? YYmap[Yid]->right[0] : N+1,N,
		   Xid,XXmap[Xid]->id,XXmap[Xid]->origmap ? XXmap[Xid]->left[0] : 0, XXmap[Xid]->origmap ? XXmap[Xid]->right[0] : M+1,M,
		   align->score,best+((bestLeft >= 0) ? ChimScore:0.0) + ((bestRight < U) ? ChimScore : 0),
		   bestLeft, (bestLeft>=0) ? align->sites1[bestLeft] : -1, (bestLeft>=0) ? align->sites2[bestLeft] : -1, 
		   bestRight, (bestRight<U) ? align->sites1[bestRight] : -1, (bestRight<U) ? align->sites2[bestRight] : -1 , U);
	    printf("iLWM=%d,iHWM=%d,LWM=%0.6f,HWM=%0.6f,score=%0.6f\n",iLWM,iHWM,LWM,HWM,score);
	    printf("align->Lend=%d,align->Rend=%d:outscore[0]=%0.6f,outscore[U=%d]=%0.6f\n",
		   align->Lend,align->Rend,align->outscore[0],U,align->outscore[U]);
	    fflush(stdout);
	    assert(bestLeft >= (align->Lend <= -2 ? 0 : -1));
	  }
	  if(DEBUG) assert(bestRight <= (align->Rend <= -2 ? U-1 : U));
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
	  FLAT(A,G, align->sites1[0], align->sites2[0]) = -2;
	  align->logPV = alignFP(A,bestI,bestJ,bestR,Y,X,N,M,orientation,Yid,Xid,delY,delX,bestscore, 0, 1);
	  if(LVERB){
	    printf("Localizing alignment:Yid=%d(id=%lld,left=%d,right=%d,N=%d),Xid=%d(id=%lld,left=%d,right=%d,M=%d):align->score=%0.6f->%0.6f,logPV=%0.4f->%0.4f,bestLeft=%d(I=%d,J=%d),bestRight=%d(I=%d,J=%d),U=%d\n",
		   Yid,YYmap[Yid]->id,YYmap[Yid]->origmap ? YYmap[Yid]->left[0] : 0, YYmap[Yid]->origmap ? YYmap[Yid]->right[0] : N+1,N,
		   Xid,XXmap[Xid]->id,XXmap[Xid]->origmap ? XXmap[Xid]->left[0] : 0, XXmap[Xid]->origmap ? XXmap[Xid]->right[0] : M+1,M,
		   oldscore,align->score,oldlogPV,align->logPV,
		   bestLeft, (bestLeft>=0) ? align->sites1[0] : 0, (bestLeft>=0) ? align->sites2[0] : 0, 
		   bestRight, (bestRight<U) ? align->sites1[align->numpairs-1] : N+1, (bestRight<U) ? align->sites2[align->numpairs-1] : M+1 , U);
	    if(PAIRSPLIT_EXTEND)
	      printf("\t extLeft=%d,%d,extRight=%d,%d,iLWM=%d,iHWM=%d,LWM=%0.6f,HWM=%0.6f,score=%0.6f\n",align->refstart,align->start,align->refend,align->end,iLWM,iHWM,LWM,HWM,score);
	    else
	      printf("\t iLWM=%d,iHWM=%d,LWM=%0.6f,HWM=%0.6f,score=%0.6f:alignment=(I=%d,J=%d) to (I=%d,J=%d)\n",iLWM,iHWM,LWM,HWM,score,
		     align->sites1[0],align->sites2[0],align->sites1[align->numpairs-1],align->sites2[align->numpairs-1]);
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
      }
#endif // PairSplit code
    }
  } // if(bestscore > origscore || VERB ...)
}


static int threads_displayed = 0;

static int A_numthreads= 0;
static long long A_maxN[MAXCOLOR] = {-1, -1};
static AL **A_threads = NULL;/* A_threads[0..A_numthreads-1][0..colors-1] */
static AL *ALmem = NULL;/* ALmem[0 .. (A_numthreads*colors)-1]; */

static int mem_numthreads= 0;
static long long Fmem_siz[MAXCOLOR] = {0,0}, Imem_siz[MAXCOLOR] = {0,0};
static PFLOAT **Fmem_threads = NULL;/* Fmem_threads[0..mem_numthreads-1][0..Fmem_siz*mem_numthreads-1] */
static int **Imem_threads = NULL;/* Imem_threads[0..mem_numthreads-1][0..Imem_siz*mem_numthreads-1] */
static long long mem_maxN[MAXCOLOR]= {-1,-1}, mem_maxM[MAXCOLOR]= {-1,-1};

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
	       AL *A, /**< A[c=0..colors-1].<field>[I=1..N][J=1..M] : Only the c and I demensions are allocated */
	       PFLOAT *Fmem,/**< preallocated memory to be used by A[c].score[][],A[c].Uscore[][], A[c].Rscore[][] */
	       int *Imem /**< preallocated memory to be used by A[c].<field>[][], for <field> = G,H,Rijy,Lijy,Rijx,Lijx */
	       )
{
  if(VERB>=2){
    #pragma omp critical
    {
      printf("pairalign:Yid=%d,Xid=%d,galign=%p\n",Yid,Xid,galign);
      fflush(stdout);
    }
  }

  if(DEBUG>=2){
    for(int c = 0; c < colors; c++){
      if(DEBUG && !(Xmap->deltaX[c] && Xmap->deltaR[c] && Ymap->deltaX[c])){
        #pragma omp critical
        {
          printf("Yid=%d(%lld,%lld),Xid=%d(%lld,%lld),c=%d:Xmap->deltaX[c]=%p,Xmap->deltaR[c]=%p,Ymap->deltaX[c]=%p\n",
		 Yid,Ymap->id,YYmap[Yid]->id,Xid,Xmap->id,XXmap[Xid]->id,c,Xmap->deltaX[c],Xmap->deltaR[c],Ymap->deltaX[c]);
	  fflush(stdout);
	  assert(Xmap->deltaX[c] && Xmap->deltaR[c] && Ymap->deltaX[c]);
	}
      }
    }
  }

  int *NN = Ymap->numsite;
  int *MM = Xmap->numsite;
  FLOAT **XX = Xmap->site;
  FLOAT **YY = Ymap->site;

  if(DEBUG>=2 && galign && galign->numpairs > 0){
    for(int c = 0; c < colors; c++){
      assert(galign[c+1].sites1[0] >= 1 && galign[c+1].sites1[galign[c+1].numpairs-1] <= NN[c]);
      assert(galign[c+1].sites2[0] >= 1 && galign[c+1].sites2[galign[c+1].numpairs-1] <= MM[c]);
    }
  }

  Calign *align, Lalign[MAXCOLOR+1];
  if(PairSplit){
    if(DEBUG) assert(galign != 0);
    align = galign;
  } else /* use a local temporary variable : only save it to galign if score is good */
    align = &Lalign[0];
  align++;

  for(int c = -1; c < colors; c++){
    align[c].mapid1 = Yid;
    align[c].mapid2 = Xid;
    align[c].orientation = -1;
    align[c].score = MINSCORE;
  }
  if(PairSplit){
    delete [] align->align2;
    align->align2 = NULL;
  }
  align[-1].numpairs = 0;

  if(NN[0]+NN[1] <= 0 || MM[0]+MM[1] <= 0)
    return;/* no alignment possible */

  /* re-allocated DP arrays A[c=0..colors-1].<field>[I=1..N][J=1..M] */
  for(int c = 0; c < colors; c++){
    A[c].G[0] = A[c].H[0] = A[c].Rijy[0] = A[c].Lijy[0] = A[c].Rijx[0] = A[c].Lijx[0] = NULL;
    A[c].score[0] = A[c].Uscore[0] = A[c].Rscore[0] = NULL;

    A[c].M_stride= ((MM[c] + 0x1f) & ~0x1f) + 16; /* This is necessary to make better use of set-associative cache */
    A[c].array_stride = A[c].M_stride * NN[c] + 16;
    if(DEBUG) assert(A[c].array_stride * 3LL <= Fmem_siz[c]);
    if(DEBUG && !(A[c].array_stride * 6LL <= Imem_siz[c])){
      #pragma omp critical
      {
	printf("Yid=%d,Xid=%d:c=%d,MM[c]=%d,NN[c]=%d,A[c].M_stride=%lld, A[c].array_stride=%lld,Imem_siz[c]=%lld\n",Yid,Xid,c,MM[c],NN[c],A[c].M_stride,A[c].array_stride,Imem_siz[c]);
	fflush(stdout);
	assert(A[c].array_stride * 6LL <= Imem_siz[c]);
      }
    }

    if(DEBUG && NN[c] > maxN[c]){
      fprintf(stderr, "**** Aieee ! INTERNAL ERROR: c=%d:N[c] > maxN[c] (%d > %lld)\n", c, NN[c], maxN[c]);
      fflush(stderr);
      assert(NN[c] <= maxN[c]);
    }
    if(DEBUG && MM[c] > maxM[c]){
      fprintf(stderr, "**** Aieee ! INTERNAL ERROR: c=%d:M[c] > maxM[c] (%d > %lld)\n", c, MM[c], maxM[c]);
      fflush(stderr);
      assert(MM[c] <= maxM[c]);
    }

    long long Fmem_offset = c ? Fmem_siz[0] : 0;
    long long Imem_offset = c ? Imem_siz[0] : 0;

#pragma ivdep
    for(int I = 1; I <= NN[c]; I++){
      A[c].G[I] = &Imem[Imem_offset + (I-1)*A[c].M_stride + 3];
      A[c].H[I] = &Imem[Imem_offset + (I-1)*A[c].M_stride + 3 + A[c].array_stride];
      A[c].Rijy[I] = &Imem[Imem_offset + (I-1)*A[c].M_stride + 3 + A[c].array_stride * 2];
      A[c].Lijy[I] = &Imem[Imem_offset + (I-1)*A[c].M_stride + 3 + A[c].array_stride * 3];
      A[c].Rijx[I] = &Imem[Imem_offset + (I-1)*A[c].M_stride + 3 + A[c].array_stride * 4];
      A[c].Lijx[I] = &Imem[Imem_offset + (I-1)*A[c].M_stride + 3 + A[c].array_stride * 5];

      A[c].score[I] = &Fmem[Fmem_offset + (I-1)*A[c].M_stride + 3];
      A[c].Uscore[I] = &Fmem[Fmem_offset + (I-1)*A[c].M_stride + 3 + A[c].array_stride];
      A[c].Rscore[I] = &Fmem[Fmem_offset + (I-1)*A[c].M_stride + 3 + A[c].array_stride * 2];
    }
  } // c = 0..colors-1

  if(DEBUG>=2){
    for(long long i = Fmem_siz[0]+Fmem_siz[1]; --i >= 0;)
      Fmem[i] = aNaN;
    for(long long i = Imem_siz[0]+Imem_siz[1]; --i >= 0;){
#ifdef VALGRIND
      Imem[i] = -555555555;
#else
      Imem[i] = -5;
#endif
    }
  }

#ifdef VALGRIND
  VALGRIND_MAKE_MEM_UNDEFINED(Fmem, (Fmem_siz[0]+Fmem_siz[1]) * sizeof(PFLOAT));
  VALGRIND_MAKE_MEM_UNDEFINED(Imem, (Imem_siz[0]+Imem_siz[1]) * sizeof(int));
#endif
  
  /* first check X in normal orientation and update best alignment */
  if(!phash || !phash->orientation)
    pairalignXY(YY,NN,XX,MM,A,align,Yid,Xid, Ymap->deltaY, Xmap->deltaX, phash, 0);
  
  if(DEBUG>=2 && align[-1].numpairs > 0){
    for(int c = 0; c < colors; c++){
      assert(align[c].numpairs > 0);
      assert(align[c].sites1[0] >= 1 && align[c].sites1[align[c].numpairs-1] <= NN[c]);
      assert(align[c].sites2[0] >= 1 && align[c].sites2[align[c].numpairs-1] <= MM[c]);
    }
  }

#ifdef VALGRIND
  VALGRIND_MAKE_MEM_UNDEFINED(Fmem, (Fmem_siz[0]+Fmem_siz[1]) * sizeof(PFLOAT));
  VALGRIND_MAKE_MEM_UNDEFINED(Imem, (Imem_siz[0]+Imem_siz[1]) * sizeof(int));
#endif

  FLOAT *Xrev0 = (FLOAT *)alloca((MM[0]+2)*sizeof(FLOAT));
  FLOAT *Xrev1 = (FLOAT *)alloca((MM[1]+2)*sizeof(FLOAT));
  FLOAT *XXrev[2] = {Xrev0,Xrev1};

  if(!(RepeatShift > 0.0 && align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= MINSITES)){

    /* next check X in reversed orientation and update best alignment */
    for(int c = 0; c < colors; c++)
      for(int J = 0; J <= MM[c] + 1; J++)
	XXrev[c][J] = XX[c][MM[c] + 1] - XX[c][MM[c] + 1 - J];

    if(!phash || phash->orientation)
      pairalignXY(YY, NN, XXrev, MM, A, align, Yid, Xid, Ymap->deltaY, Xmap->deltaR, phash, 1);
    else if (phash[0].id1==phash[1].id1 && phash[0].id2==phash[1].id2 && phash[1].orientation)
      pairalignXY(YY, NN, XXrev, MM, A, align, Yid, Xid, Ymap->deltaY, Xmap->deltaR, &phash[1], 1);
  }

  if(DEBUG>=2 && align[-1].numpairs > 0){
    for(int c = 0; c < colors; c++){
      assert(align[c].sites1[0] >= 1 && align[c].sites1[align[c].numpairs-1] <= NN[c]);
      assert(align[c].sites2[0] >= 1 && align[c].sites2[align[c].numpairs-1] <= MM[c]);
    }
  }

  /* If PairSplit chimeric/local splits will be checked in caller */

  if(VERB>=2 /* || EVERB */){
    #pragma omp critical
    {
      if(align[-1].numpairs >= 1)
	printf("pairalign:Yid=%d(id=%lld),Xid=%d(id=%lld):Score=%0.6f,logPV=%0.2f,orientation=%d,aligned sites=%d,repeat=%d\n",
	       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,align[-1].score,align[-1].logPV,
	       align[-1].orientation,align[-1].numpairs,align[-1].repeat);
      else
	printf("pairalign:Yid=%d(id=%lld),Xid=%d(id=%lld):orientation=%d,aligned sites=%d,repeat=%d\n",
	       Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,align[-1].orientation,align[-1].numpairs,align[-1].repeat);
      fflush(stdout);
    }
  }

  if(!PairSplit){
    if(align[-1].score > ScoreThreshold && align[-1].logPV > LogPvThreshold && align[-1].numpairs >= MINSITES){/* save local align to galign */
      align--;
      if(DEBUG) assert(align == &Lalign[0]);

      if(!galign)
	galign = new Calign[1+colors];

      for(int c = 0; c <= colors; c++){
	/* swap galign[c] and align[c] : align[c] is on the stack and will automatically be deallocated */
	Calign tmp = galign[c];
	galign[c] = align[c];
	align[c] = tmp;
	if(DEBUG>=2 && c > 0 && galign[c].numpairs > 0){
	  assert(galign[c].sites1[0] >= 1 && galign[c].sites1[galign[c].numpairs-1] <= NN[c-1]);
	  assert(galign[c].sites2[0] >= 1 && galign[c].sites2[galign[c].numpairs-1] <= MM[c-1]);
	}
      }
      align = galign;
      align++;
	
      if(VERB>=2){
	printf("pairalign:Yid=%d,Xid=%d:numpairs=%d,N=%d,%d,M=%d,%d,(score=%0.3f,LogPV=%0.3f,repeat=%d)\n",
	       Yid,Xid,align[-1].numpairs,NN[0],NN[1],MM[0],MM[1], align[-1].score,align[-1].logPV,align[-1].repeat);
        fflush(stdout);
      }
    } else {
      //      if(DEBUG && !EVERB) assert(align[-1].numpairs==0);
      //      align[-1].allfree();
    }
  }

  if(DEBUG>=2 && align[-1].numpairs > 0){
    for(int c = 0; c < colors; c++){
      assert(align[c].sites1[0] >= 1 && align[c].sites1[align[c].numpairs-1] <= NN[c]);
      assert(align[c].sites2[0] >= 1 && align[c].sites2[align[c].numpairs-1] <= MM[c]);
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
    printf("maxhashalloc:hashpairs=%p,maxhashpairs=%llu: wall time=%0.6f secs\n",hashpairs,(unsigned long long)maxhashpairs, wtime());
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

/** stream in enough hashtable entries to fill the entire table hashpairs[0..maxhashpairs-1] (before filtering by hashscore).
    Then copy hashtable entries to YidList[],XidList[],pHashList[] until we have MaxPair Entries (or we reach end of hashtable) and return number of entries.
    First check if the 2nd buffer has data : if so just copy from 2nd buffer instead of reading from disk */
static int loadhash(Chashpair* &nexthash, Chashpair* &hashpairs, size_t &numhashpairs, size_t maxhashpairs, int *YidList, int *XidList, Chashpair **phashList)
{
  /* make sure any existing entries are at start of hashpairs[] : typically there will be at most one remaining entry */
  Chashpair *p = hashpairs;
  Chashpair *end = &hashpairs[numhashpairs];
  if(nexthash > p) { /* copy existing entries to start of hashpairs[] */
    if(VERB && HASH_DEBUG && end > nexthash){
      printf("copying %llu entries from hashpairs%s[%llu..%llu] to hashpairs%s[0..%llu]: cpu time=%0.6f, wall time=%0.6f secs\n",
	     (unsigned long long)(end-nexthash), hashpairs==hashpairs2 ? "2" : "", (unsigned long long)(nexthash-hashpairs), (unsigned long long)(end-hashpairs)-1, hashpairs==hashpairs2 ? "2" : "", 
	     (unsigned long long)(end-nexthash)-1, mtime(), wtime());
      fflush(stdout);
    }

    long int len = end - nexthash;

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
	printf("Swapping Buffers: NumPair1=%d,Matches1=%llu, NumPair2=%d, Matches2=%llu: CPU time=%0.6f, wall time=%0.6f\n", 
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

      return NumPair1;
    }
  }

  int HashEof = hash_eof(HashFP);
  if(HashEof && cur==0)
    return 0;

  if(!HashEof){
    int linecnt = 0;/* Not Used for Binary File input */
    if(VERB && HASH_DEBUG){
      printf("cur=%llu,num=%llu,maxhashpairs=%llu:Calling hash_read(HashFP,0,&hashpairs[cur],num,linecnt): CPU time=%0.6f, wall time=%0.6f secs\n",
	     (unsigned long long)cur,(unsigned long long)num,(unsigned long long)maxhashpairs, mtime(), wtime());
      fflush(stdout);
    }
    if(DEBUG) assert(p == &hashpairs[cur]);
    numhashpairs += hash_read(HashFP, 0, p, num, linecnt);
    hashpairs[numhashpairs].id1 = 0;
    HashEof = hash_eof(HashFP);
    if(VERB && HASH_DEBUG && numhashpairs > cur){
      printf("Read in next %llu Matches (total=%llu/%llu) from %s (EOF=%d) into hashpairs%d[]: cpu time=%0.6f, wall time=%0.6f secs\n", 
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
	exit(1);
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
  numhashpairs = j;
  hashpairs[numhashpairs].id1 = 0;

  /* copy hashtable entries to YidList[], XidList[], phashList[] */
  if(DEBUG && hashpairs == hashpairs1) assert(NumPair1 == 0);
  if(DEBUG && hashpairs == hashpairs2) assert(NumPair2 == 0);

  /* just copy hashtable entries to YidList[],XidList[],phashList[] until we reach end of hashtable buffer */
  end = &hashpairs[numhashpairs];

  long long lastYid = -1, lastXid = -1;
  int NumPair, TNumPair = 0;
  for(NumPair= 0; NumPair < MaxPair; NumPair++){
    if(nexthash >= end-1){/* next to end of current read buffer : check if this is EOF */
      if(!HashEof)
	break;/* reuse last entry in next Yid iteration since it may come in a pair of entries (both orientations) that are needed together */
      /* EOF : force break out of Yid loop, unless the last entry still needs to be processed */
      if(nexthash >= end)
	break;
    }

    if(DEBUG>=2) assert(0 <= nexthash->id1 && nexthash->id1 < numY && YYmap[nexthash->id1]->id >= 0);
    if(DEBUG>=2) assert(0 <= nexthash->id2 && nexthash->id2 < numX && XXmap[nexthash->id2]->id >= 0);

    YidList[NumPair] = lastYid = nexthash->id1;
    XidList[NumPair] = lastXid = nexthash->id2;

    phashList[NumPair] = nexthash++;
    TNumPair = NumPair+1;

    /* skip next entry in hashtable if it is the same pair of maps as lastYid,lastXid :
       This happens when both orientation are in the HashTable */
    if(nexthash < end && nexthash->id1 == lastYid && nexthash->id2 == lastXid)
      nexthash++;
  }
  if(DEBUG) assert(NumPair == TNumPair);

  if(VERB && HASH_DEBUG && HASH_STREAM>=2){
    printf("Converted into %d map pairs (%llu Matches remaining):cpu time=%0.6f,wall time=%0.6f\n",
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

  return NumPair;
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

extern void pairmerge(int pairmergeIter);/* see pairalign.cpp : NOTE : 2-color code should be put in a separate function that is local here */

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
int pairalign_allpairs2(int mapstart,int mapend, int pairmergeIter)
{
  if(DEBUG) assert(mapstart == 0 && mapend == nummaps - 1);
  repeatcnt = 0;

  if(DEBUG>=2 && !hash_filename){
    for(int mapid = 0; mapid < nummaps; mapid++){
      if(!(Gmap[mapid]->mapid == mapid)){
	printf("mapid=%d:Gmap[mapid]=%p,Gmap[mapid]->mapid=%d,Gmap[mapid]->id=%lld,nummaps=%d\n",
	       mapid,Gmap[mapid],Gmap[mapid]->mapid,Gmap[mapid]->id,nummaps);
	for(int m = 0; m < 10; m++)
	  printf("m=%d:Gmap[m]=%p,Gmap[m]->mapid=%d,Gmap[m]->id=%lld\n",
		 m,Gmap[m],Gmap[m]->mapid,Gmap[m]->id);
	fflush(stdout);
	assert(Gmap[mapid]->mapid == mapid);
      }
    }
  }

  if(!PairSplit)
    PairsplitSeparateQueries = 0;

  if(VERB>=2){
    printf("pairalign:PVERB=%d,LVERB=%d\n",PVERB,LVERB);
    fflush(stdout);
  }
  if(ROC){
    printf("ROC=%d no longer supported\n", ROC);
    exit(1);
  }
  if(DEBUG && DEFER_G != 2){
    printf("DEFER_G=%d not supported\n",DEFER_G);
    exit(1);
  }
  if(DEBUG && !DEFER_BP){
    printf("DEFER_BP=%d not supported\n",DEFER_BP);
    exit(1);
  }
  if(DEBUG && !LOCAL){
    printf("LOCAL=%d not supported\n",LOCAL);
    exit(1);
  }
  /*  if(DEBUG && !ENDFIX){
    printf("ENDFIX=%d not supported\n",ENDFIX);
    exit(1);
    }*/
  if(DEBUG && PAIRSCORETYPE != 5){
    printf("PAISCORETYPE=%d not supported(must be 5)\n",PAIRSCORETYPE);
    exit(1);
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

  for(int c = 0; c < colors; c++){
    if(resSD[c] > 0.0){
      printf("c=%d:resSD[c] = %0.6f : Not implemented for pairwise alignment: using resSD[c]=0\n", c,resSD[c]);
      resSD[c] = 0.0;
    }

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
    if(SD[c] < MINSD){
      printf("SD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SD[c],MINSD,MAXSD,MINSD);
      SD[c] = MINSD;
    } else if(SD[c] > MAXSD){
      printf("SD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SD[c],MINSD,MAXSD,MAXSD);
      SD[c] = MAXSD;
    }
    if(SF[c] < MINSF){
      printf("SF[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SF[c],MINSF,MAXSF,MINSF);
      SF[c] = MINSF;
    } else if(SF[c] > MAXSF){
      printf("SF[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c,SF[c],MINSF,MAXSF,MAXSF);
      SF[c] = MAXSF;
    }
  }

  if(PairSplit){    /* initialize PairHash[] */
    PairHash_init();
    if(DEBUG) assert(PairHash);
  }

  if(PoutlierEnd > 0.9999)
    PoutlierEnd = 0.9999;

  /* initialize statistics counter */
  int mapcnt=0;/* map pairs aligned above threshold */
  double len[2] = {0.0,0.0};/* sum of aligned lengths (based on overlap region) */
  double lenB[2] = {0.0,0.0};/* sum of aligned lengths (based on reference) - resKB2 */
  int sitecnt[2] = {0,0}, misscnt[2] = {0,0};/* count of total and misaligned sites in aligned regions of all pairs */
  int segcnt[2] = {0,0};/* count of aligned intervals */
  double errnorm[2] = {0.0,0.0};/* sum of (x-y)*(x-y)/(2*SF[c]*SF[c] + fabs(SD[c])*SD[c]*(x+y) + SR[c]*SR[c]*(x*x+y*y))  */

  int outliertot[2] = {0,0};/* all intervals that are candidates for outliers (eg x < y) */
  int outliercnt[2] = {0,0};
  double outliersum[2] = {0.0, 0.0};

  int poscnt = 0;/* count of map pairs with +ve score */

  long long tcnt = 0;/* total pairwise alignment attempts */
  int paircnt = 0;/* number of pairs with minimum overlap of MINOVERLAP * min(len1,len2) */
    
  for(int c = 0; c < colors;c++){
    numerrors[c] = 0;
    maxerrors[c] = 1024;
    errors[c] = new Cinterval[maxerrors[c]];

    nummisaligns[c] = 0;
    maxmisaligns[c] = 1024;
    misaligns[c] = new Cmisalign[maxmisaligns[c]];
  }

  gnumthreads = 1;
  int maxthreads = 1;
  #ifdef _OPENMP
  maxthreads = omp_get_max_threads();
  maxthreads = min(maxthreads,MaxThreads);
  gnumthreads = maxthreads;
  #endif

  score_init(0,nummaps-1);
  
  double ExpOutlierPenalty = exp(OutlierPenalty);
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

  if(!Cfirst || CfirstPairs == 2){
    if(hash_filename || PairMerge)    /* sort maps[0..nummaps-1] by id in ascending order */
      qsort(YYmap,numY,sizeof(Cmap *),(intcmp*)CmapIdInc);
    if(DEBUG) assert(!PairSplit);
    YYmap = XXmap = &Gmap[0];
    if(CfirstPairs == 2){/* compute separte maxM for XXmap[Cfirst .. numX-1] */
      for(int c = 0; c < colors; c++)
	maxM[c] = 0;
      for(int Xid = 0; Xid < numX; Xid++)
	for(int c = 0; c < colors; c++)
	  if(XXmap[Xid]->numsite[c] > maxM[c])
	    maxM[c] = XXmap[Xid]->numsite[c];
    }
    numX = numY = nummaps;
    if(VERB){
      if(CfirstPairs == 2)
	printf("Gmap[%d..%d] duplicated as XXmap[0..%d] and YYmap[0..%d] (only XXmap[%d..%d] will be used):maxN=%lld,%lld,maxM=%lld,%lld\n",
	       0,nummaps-1,numX-1,numY-1, Cfirst, numX-1,maxN[0],maxN[1],maxM[0],maxM[1]);
      else
	printf("Gmap[%d..%d] duplicated as XXmap[0..%d] and YYmap[0..%d]:maxN=maxM=%lld,%lld\n",
	       0,nummaps-1,numX-1,numY-1,maxN[0],maxN[1]);
      fflush(stdout);
    }
  } else { /* copy Gmap[0..Cfirst-1] to refmap[0..Cfirst-1] */
    if(DEBUG) assert(Cfirst > 0);

    if(hash_filename || PairMerge){    /* sort maps by id in ascending order */
      qsort(Gmap,Cfirst,sizeof(Cmap *),(intcmp*)CmapIdInc);
      qsort(&Gmap[Cfirst],nummaps-Cfirst,sizeof(Cmap *),(intcmp*)CmapIdInc);
    }

    numY = numrefmaps = Cfirst;
    if(PairSplit){
      if(VERB>=2){
	printf("Calling maxmapalloc(%d,%d,refmap=%p,0)\n",Cfirst,maxrefmaps,refmap);
	fflush(stdout);
      }
      maxmapalloc(numY,maxrefmaps,refmap,0,1);
      YYmap = refmap;
      for(int i = 0; i < Cfirst; i++){
	refmap[i] = Gmap[i];
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
    for(int c = 0; c < colors; c++)
      maxM[c] = maxN[c] = 0;
    for(int Xid = 0; Xid < numX; Xid++)
      for(int c = 0; c < colors; c++)
	if(XXmap[Xid]->numsite[c] > maxM[c])
	  maxM[c] = XXmap[Xid]->numsite[c];
    for(int Yid = 0; Yid < numY; Yid++)
      for(int c = 0; c < colors; c++)
	if(YYmap[Yid]->numsite[c] > maxN[c])
	  maxN[c] = YYmap[Yid]->numsite[c];
    if(VERB){
      printf("Gmap[0..%d] renamed as YYmap[0..%d] and Gmap[%d..%d] renamed as XXmap[0..%d](maxN=%lld,%lld,maxM=%lld,%lld)\n",
	     Cfirst-1,numY-1,Cfirst,nummaps-1,numX-1,maxN[0],maxN[1],maxM[0],maxM[1]);
      if(VERB>=2){
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
	  exit(1);
	}
	if(DEBUG) assert(YYmap[Yid]->origmap == NULL);
	*YYmap[nYid] = *YYmap[Yid];/* just copy the struct since only the scalar fields will differ */
	YYmap[nYid]->Xid = Xid;
	YYmap[nYid]->mapid = nYid;
	YYmap[nYid]->origmap = YYmap[Yid];
	if(DEBUG) assert(YYmap[nYid]->origmap->origmap == NULL);

	for(int c = 0; c < colors; c++){
	  YYmap[nYid]->left[c] = 0;
	  YYmap[nYid]->right[c] = YYmap[Yid]->numsite[c]+1;
	}
      }
    }
    numrefmaps = numY *= numX;
  }

  int LocalCnt = 0;/* count of local alignments (if PairSplit) */

  //  int binM1 = bin - 1;
  //  int binM1N = bin - 1 - numbins;

  #ifdef _OPENMP
  if(MaxMem > 0 && max(maxM[0],maxM[1]) > 0){/* NOTE : maxM,maxN is the largest Xmap,Ymap size */
    double memperthread = 0.0;
    for(int c = 0; c < colors; c++)
      memperthread += (double)(maxN[c]+1) * (double)(maxM[c]  * sizeof(DP)) + sizeof(AL);
    double maxmem = (double)MaxMem * 1024.0 * 1024.0 * 1024.0;
    int limitthreads = max(1,(int)floor(maxmem/memperthread + 0.5));
    if(limitthreads < maxthreads){
      if(VERB){
	printf("max number of threads reduced from %d to %d due to MaxMem=%0.1f Gbytes, per thread memory = %0.3f Gbytes:maxN=%lld,%lld,maxM=%lld,%lld\n",
	       maxthreads,limitthreads, MaxMem, memperthread/(1024.0*1024.0*1024.0),maxN[0],maxN[1],maxM[0],maxM[1]);
	fflush(stdout);
      }
      if(DEBUG) assert(limitthreads >= 1);
      maxthreads = limitthreads;
    } else if(VERB){
      printf("max number of threads not reduced from %d due to MaxMem=%0.1f Gbytes, per thread memory = %0.3f Gbytes:maxN=%lld,%lld,maxM=%lld,%lld\n",
	     maxthreads, MaxMem, memperthread/(1024.0*1024.0*1024.0),maxN[0],maxN[1],maxM[0],maxM[1]);
      fflush(stdout);
    }
  }
  #endif // _OPENMP

  MaxPair = min(max(numY*numX*1.0,nummaps*(nummaps-1.0)*0.5), max(nummaps*2.0, maxthreads*((FirstAlignments > 0) ? 1024.0 : 1024.0*64.0)));
  if(VERB>=2){
    printf("maxthreads=%d,MaxPair=%d,FirstAlignments=%d,nummaps=%d,numY=%d,numX=%d\n",maxthreads,MaxPair,FirstAlignments,nummaps,numY,numX);
    fflush(stdout);
  }
  if(DEBUG) assert(MaxPair > 0);

  BinaryFormat = 0;

  if(hash_filename){/* read in hash pairs */
#if HASH_STREAM < 2
    printf("HASH_STREAM=%d no longer supported (see hash.h)\n",HASH_STREAM);
    exit(1);
#endif
    if(!HashGen){
      printf("Previously computed hash table no longer supported : use -hashgen\n");
      exit(1);
    }

    if(DEBUG) assert(hashpairs1 == 0 && numhashpairs1 == 0);

    size_t L = strlen(hash_filename);
    if(L >= strlen(".hashbin") && !strcmp(&hash_filename[L-strlen(".hashbin")],".hashbin")){/* binary file input */
      BinaryFormat = 1;/* used later to enable streaming file input */

      /* HASH_STREAM : set up streamed IO */
      if(RepeatShift > 0.0){
	printf("-hash not supported with -Repeat\n");
	exit(1);
      }
      if(PairSplit){
	printf("-hash not yet supported with -pairsplit\n");
	exit(1);
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
      exit(1);
    }

    hashpairs1[numhashpairs1].id1 = 0;
  }
  nexthash1 = hashpairs1;
  nexthash2 = hashpairs2;
  if(hashsplit){
    printf("-hashsplit no longer supported\n");
    exit(1);
  }

  //  size_t numaligns_compressed = 0;/* last alignment[] index already checked for NULL entries (if PairSplit) */
  //  size_t prev_aligncnt = 0;/* accumulated count of alignments above threshold and chimpair == 0 in alignment[0..numaligns_compress-1] (if PairSplit) */
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

    for(int Yid = 0;Yid < numY; Yid++){
      if(LVERB && ((Cfirst && CfirstPairs != 1) ? (CfirstPairs==2 ? max(Cfirst,Yid+1) : 0) : Yid+1) <= numX-1){
	printf("aligning YYmap[%d/%d](id=%lld,centercloned=%d) with XXmap[%d..%d]\n", Yid, numY, YYmap[Yid]->id,YYmap[Yid]->centercloned, Cfirst ? 0 : Yid+1, numX-1);
	fflush(stdout);
      }
      if(DEBUG>=2) assert(YYmap[Yid]->id >= 0);
      if(VERB>=2){
	printf("Yid=%d(id=%lld)\n",Yid,YYmap[Yid]->id);
	fflush(stdout);
      }

      int Yunpaired,YYid;
      if(DEBUG>=2){
	Yunpaired = (pairmergeIter > 1 && !YYmap[Yid]->paired) ? 1 : 0;
	YYid = YYmap[Yid]->id;
      }

      if(FirstAlignments > 0 && mapcnt >= FirstAlignments)
	break;

      if(RepeatShift <= 0.0 && !nexthash1 && !PairSplit){ /* fast special case (without hashtable) */
	int Xid = (Cfirst && CfirstPairs != 1) ? (CfirstPairs==2 ? max(Cfirst,Yid+1) : 0) : Yid+1;
	
	// Tried multi-threading the following loop but it slows down
	for(; Xid < numX; Xid++, NumPair1++){
	  if(LVERB>=2){
	    printf("aligning YYmap[%d/%d](centercloned=%d) with XXmap[%d/%d](centercloned=%d)\n", Yid, numY, YYmap[Yid]->centercloned, Xid, numX, XXmap[Xid]->centercloned);
	    fflush(stdout);
	  }

#if 1 // goes with new code at end of pairalign : currently does not work
	  if(DEBUG>=2 && (Yunpaired && !XXmap[Xid]->paired)){
	    printf("aligning YYmap[%d/%d](paired=%d) with XXmap[%d/%d](paired=%d),Yunpaired=%d\n",Yid,numY,YYmap[Yid]->paired,Xid,numX,XXmap[Xid]->paired,Yunpaired);
	    fflush(stdout);
	    assert(!(Yunpaired && !XXmap[Xid]->paired));
	  }
#else
	  if(Yunpaired && !XXmap[Xid]->paired)
	    continue;
#endif
	  if(DEBUG>=2) assert(YYid != XXmap[Xid]->id);

	  phashList1[NumPair1] = 0;
	  YidList1[NumPair1] = Yid;
	  XidList1[NumPair1] = Xid;

	  if(VERB>=2){
	    printf("Yid=%d,Xid=%d:NumPair=%d,MaxPair=%d\n",Yid,Xid,NumPair1+1,MaxPair);
	    fflush(stdout);
	  }
	}
      } else if(HASH_STREAM && HashFP /* && RepeatShift <= 0.0 && !PairSplit*/){/* fast special case (with streamed hashtable) : used for pairwise molecule alignment */
	if(DEBUG>=2) assert(nexthash1);
	if(DEBUG>=2) assert(BinaryFormat==1);
	if(DEBUG>=2) assert(numbins==1);
	if(DEBUG>=2) assert(RepeatShift <= 0.0);
	if(DEBUG>=2) assert(!PairSplit);

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
	if(!NumPair1)
	  NumPair1 = loadhash(nexthash1, hashpairs1, numhashpairs1, maxhashpairs1, YidList1,XidList1,phashList1);
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

      if(VERB /* && !PairSplit */ && ((HASH_STREAM && HashFP) || PairCnt + NumPair1 > LastPairCnt + maxthreads)){
	if(RepeatShift <= 0.0){
	  double completion = Cfirst ? ((double)YidList1[0])/numY : 1.0 - (double)(numY-YidList1[0])*(numY-YidList1[0]-1.0)/(numY*(numY-1.0));
	  if(DEBUG && NumPair1 > 0) assert(YidList1[0] <= YidList1[NumPair1-1]);
	  printf("%d:aligning YYmap[%d..%d] with XXmap[%d..%d], NumPair=%d/%d (previous alignments found=%d/%llu) cpu time=%0.6f, wall time=%0.6f %4.1f%% complete\n", 
		 wcnt,YidList1[0], YidList1[NumPair1-1], Cfirst ? 0 : max(skipcnt,YidList1[0]+1), numX-1, NumPair1, MaxPair, mapcnt, (unsigned long long)PairCnt, mtime(), wtime(), completion*100.0);
	}
	if(RepeatShift > 0.0)
	  printf("Checking repeats in XXmap[%d..%d] (total repeats found=%d)\n", LastRepeat, XidList1[NumPair1-1], mapcnt);
	fflush(stdout);
	LastPairCnt = NumPair1 + PairCnt;
	LastRepeat = XidList1[NumPair1-1]+1;
      }
      if(NumPair1 <= 0)
      	continue;
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
	  printf("reduced numaligns from %llu to %llu, NumPair=%d\n",
		 (unsigned long long)numaligns,(unsigned long long)j,NumPair1);
	  fflush(stdout);
	}
	numaligns = j;
      }
      numaligns_start = numaligns;/* no need to compress alignment[0..numalign_start-1] */

      /* preallocate space for alignment pointers, to minimize allocation in multi-threaded region */
      maxalignallocNULL(numaligns+NumPair1, alignment, numaligns, maxaligns);

      int numthreads = 1;
      #ifdef _OPENMP
      numthreads = maxthreads;
      if(numthreads > NumPair1)
	numthreads = max(1,NumPair1);
      #endif

      int block = (PairSplit || PairMerge) ? 1 : min(16,max(1,NumPair1/(64*numthreads)));

      if(VERB && numthreads > 1 && !threads_displayed) {
	printf("Yid=%d:nummaps=%d:Using %d threads,block size=%d(NumPair=%d,MaxPair=%d)\n",Yid,nummaps,numthreads,block,NumPair1,MaxPair);
	fflush(stdout);
	threads_displayed = 1;
      }

      //      int origYid = Yid;
      
      try {/* following code may throw out-of-memory exception */
	if(A_numthreads < numthreads || A_maxN[0] < maxN[0] || A_maxN[1] < maxN[1]) {
	  if(VERB>=2){
	    printf("Allocating A_threads:maxN=%lld,%lld,numthreads=%d\n",maxN[0],maxN[1],numthreads);
	    fflush(stdout);
	  }

	  for(int i=0; i < A_numthreads; i++) {
	    for(int c = 0; c < colors; c++){
	      delete [] A_threads[i][c].score;
	      delete [] A_threads[i][c].Uscore;
	      delete [] A_threads[i][c].Rscore;
	      delete [] A_threads[i][c].G;
	      delete [] A_threads[i][c].H;
	      delete [] A_threads[i][c].Rijy;
	      delete [] A_threads[i][c].Lijy;
	      delete [] A_threads[i][c].Rijx;
	      delete [] A_threads[i][c].Lijx;
	    }
	  }
	  delete [] ALmem; ALmem = NULL;
	  delete [] A_threads; A_threads = NULL;
		
	  if(numthreads > A_numthreads) A_numthreads = numthreads;
	  for(int c = 0; c < colors; c++)
	    if(maxN[c] > A_maxN[c]) A_maxN[c] = maxN[c];
		
	  A_threads = new AL*[A_numthreads];
	  ALmem = new AL[A_numthreads*colors];
	  for(int i= 0; i < A_numthreads; i++) {
	    A_threads[i] = &ALmem[i*colors];
	    for(int c = 0; c < colors; c++){
	      A_threads[i][c].score = new PFLOAT*[A_maxN[c] + 1];
	      A_threads[i][c].Uscore = new PFLOAT*[A_maxN[c] + 1];
	      A_threads[i][c].Rscore = new PFLOAT*[A_maxN[c] + 1];
	      A_threads[i][c].G = new int*[A_maxN[c] + 1];
	      A_threads[i][c].H = new int*[A_maxN[c] + 1];
	      A_threads[i][c].Rijy = new int*[A_maxN[c] + 1];
	      A_threads[i][c].Lijy = new int*[A_maxN[c] + 1];
	      A_threads[i][c].Rijx = new int*[A_maxN[c] + 1];
	      A_threads[i][c].Lijx = new int*[A_maxN[c] + 1];
	    }
	  }
	}
      
	if(mem_numthreads < numthreads || mem_maxN[0] < maxN[0] || mem_maxN[1] < maxN[1] || mem_maxM[0] < maxM[0] || mem_maxM[1] < maxM[1]) {
	  if(VERB>=2){
	    printf("Allocating mem_thread_arrays: mem_maxN=%lld,%lld mem_maxM=%lld,%lld mem_numthreads=%d\n", 
		   mem_maxN[0],mem_maxN[1],mem_maxM[0],mem_maxM[1], mem_numthreads);
	    fflush(stdout);
	  }
	  for(int i=0; i < mem_numthreads; i++) {
	    delete [] Fmem_threads[i];
	    delete [] Imem_threads[i];
	  }
	  if(Fmem_threads) free(Fmem_threads);
	  if(Imem_threads) free(Imem_threads);

	  if(numthreads > mem_numthreads) mem_numthreads = numthreads;
	  for(int c = 0; c < colors; c++){
	    if(maxN[c] > mem_maxN[c]) mem_maxN[c] = maxN[c];
	    if(maxM[c] > mem_maxM[c]) mem_maxM[c] = maxM[c];
	  }
	     
	  Fmem_threads = new PFLOAT*[mem_numthreads];
	  Imem_threads = new int*[mem_numthreads];

	  for(int c = 0; c < colors; c++){
	    Fmem_siz[c] = mem_maxN[c] * ((mem_maxM[c] + 64) + 16) * 3LL + 16;
	    if(DEBUG) assert(Fmem_siz[c] >= 0);
	  }

	  for(int c = 0; c < colors; c++){
	    Imem_siz[c] = mem_maxN[c] * ((mem_maxM[c] + 64) + 16) * 6LL + 16;
	    if(DEBUG) assert(Imem_siz[c] >= 0);
	  }

	  Fmem_threads[0] = new PFLOAT[(Fmem_siz[0]+Fmem_siz[1]) * mem_numthreads];
	  Imem_threads[0] = new int[(Imem_siz[0]+Imem_siz[1]) * mem_numthreads];
	  for(int t= 1; t < mem_numthreads; t++) {
	    Fmem_threads[t] = &Fmem_threads[0][t * (Fmem_siz[0] + Fmem_siz[1])];
	    Imem_threads[t] = &Imem_threads[0][t * (Imem_siz[0] + Imem_siz[1])];
	  }
	  if(VERB>=2){
	    printf("Allocation complete: mem_maxN=%lld,%lld mem_maxM=%lld,%lld\n", mem_maxN[0], mem_maxN[1], mem_maxM[0], mem_maxM[1]);
	    fflush(stdout);
	  }
	}
      } catch (exception& e){
	cout << e.what() << endl;
	printf("Memory Allocation failed:mem_maxN=%lld,%lld mem_maxM=%lld,%lld,A_numthreads=%d\n",
	       mem_maxN[0],mem_maxN[1],mem_maxM[0],mem_maxM[1],A_numthreads);
	exit(1);
      }

      if(DEBUG) assert(NumPair1 <= MaxPair);

      if(VERB>=2 && numthreads > 1){
	printf("Entering parallel section\n");
	fflush(stdout);
      }
  
      FLOAT varSF2[MAXCOLOR], varSD[MAXCOLOR],varSR[MAXCOLOR];
      for(int c = 0; c < colors; c++){
	varSF2[c] = 2.0*SF[c]*SF[c];
	varSD[c] = fabs(SD[c]) * SD[c];
	varSR[c] = SR[c]*SR[c];
      }

      #pragma omp parallel num_threads(numthreads) if (numthreads > 1)
      {
	int tid = 0;
        #ifdef _OPENMP
	tid = omp_get_thread_num ();
        #endif

	/* allocate thread local DP array A[c=0..colors-1][I=1..maxN[c]][J=1..maxM[c]] */
	AL *A = A_threads[tid];
	PFLOAT *Fmem = Fmem_threads[tid];
	int *Imem = Imem_threads[tid];

	/* thread local variables */
	int Tmapcnt = 0, Tpaircnt=0;
	long long Ttcnt = 0;// thread local versions of mapcnt,paircnt and tcnt
	int Tposcnt = 0;
	int Tsitecnt[2] = {0,0};
	int Tmisscnt[2] = {0,0};
	int Toutliertot[2] = {0,0};
	int Toutliercnt[2] = {0,0};
	int Tsegcnt[2] = {0,0};
	double Tlen[2] = {0.0,0.0};
	double TlenB[2] = {0.0,0.0};
	double Toutliersum[2] = {0.0,0.0};
	double Terrnorm[2] = {0.0,0.0};

	if(VERB>=2){
	  printf("NumPair=%d,dynamic block=%d,numthreads=%d:alignments so far=%d\n",NumPair1,block,numthreads,mapcnt);
	  fflush(stdout);
	}

	#pragma omp for schedule(dynamic,block)
	for(int PairId = 0; PairId < NumPair1; PairId++){

	  if(HASH_STREAM>=2 && !tid && HashFP /* && RepeatShift <= 0.0 && !PairSplit */ && NumPair2 == 0){
	    /* Read ahead from disk into a 2nd buffer in the master thread. In Non-parallel section swap the 2 buffers instead of reading from disk */
	    NumPair2 = loadhash(nexthash2,hashpairs2,numhashpairs2,maxhashpairs2,YidList2,XidList2,phashList2);
	  }

	  Chashpair *phash = phashList1[PairId];
	  int Yid = YidList1[PairId];
	  int Xid = XidList1[PairId];
	  if(DEBUG && (Yid<0 || Xid<0)) {
	    printf("INTERNAL ERROR: negative id encountered: Xid=%d Yid=%d PairId=%d,NumPair=%d\n", Yid, Xid, PairId,NumPair1);
	    fflush(stdout);
	    assert(Yid >= 0 && Xid >= 0);
	  }
	  Cmap *Ymap = YYmap[Yid];
	  Cmap *Xmap = XXmap[Xid];
	  size_t alignid = numaligns_start + PairId; 
	  Calign * &align = alignment[alignid];

	  if(LVERB){
	    printf("aligning YYmap[%d/%d](centercloned=%d) with XXmap[%d/%d](centercloned=%d)\n", Yid, numY, YYmap[Yid]->centercloned, Xid, numX, XXmap[Xid]->centercloned);
	    fflush(stdout);
	  }

	  if(PairSplit){
	    if(DEBUG>=2) assert(PairHash);
	    if(!align)
	      align = new Calign[1+colors];
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
	      printf("Calling pairalign:Yid=%d/%d(id=%lld),Xid=%d/%d(id=%lld)\n",
		     Yid,numY,Ymap->id,Xid,numX,Xmap->id);
	      fflush(stdout);
	      assert(Xmap->id >= 0);
	      assert(Ymap->id >= 0);
	    }
	  }

#if 0 // HERE HERE : update PairSplit code 
	  if(PAIRSPLIT_EXTEND && PairSplit){
	    align->refstart = 0;
	    align->start = 0;
	    align->refend = Ymap->numsite[0]+1;
	    align->end = Xmap->numsite[0]+1;
	  }
#endif

	  if(DEBUG>=2) assert(0 <= Xid && Xid < numX);
	  if(DEBUG>=2) assert(0 <= Yid && Yid < numY);

	  pairalign(Ymap, Xmap, align, Yid, Xid, phash, A, Fmem, Imem);

	  if(DEBUG>=2) assert(align==alignment[alignid]);
	  if(DEBUG>=2 && PairSplit) assert(align && align[-1].mapid1 == Yid && align[-1].mapid2 == Xid);

	  if((DEBUG>=2 || VERB>=2) && align && align->numpairs > 0){
	    for(int c = 0; c < colors; c++){
	      int N = Ymap->numsite[c];
	      int M = Xmap->numsite[c];
	      if(!(align[1+c].sites1[0] >= 1 && align[1+c].sites1[align[1+c].numpairs-1] <= N)){
		printf("pairalign_allpairs:Yid=%d,Xid=%d,c=%d:align[c].numpairs=%d,I=%d..%d(N=%d),align[c].sites1=%p,align[c]=%p,align[c].logPV=%0.2f,align[c].repeat=%d,IsRepeatRegion()=%d\n",
		       Yid,Xid,c,align[1+c].numpairs,align[1+c].sites1[0],align[1+c].sites1[align[1+c].numpairs-1],N,align[1+c].sites1,&align[c+1],align[1+c].logPV,align[1+c].repeat,align->IsRepeatRegion());
		fflush(stdout);
		assert(align[1+c].sites1[0] >= 1 && align[1+c].sites1[align[1+c].numpairs-1] <= N);
	      }
	      assert(align[1+c].sites2[0] >= 1 && align[1+c].sites2[align[1+c].numpairs-1] <= M);
	    }
	  }

	  if(align && align->numpairs > 0 && align->score > MINSCORE){
	    /* accumulate true positive count (paircnt), if possible */
	    FLOAT **YY = Ymap->site;
	    FLOAT **XX = Xmap->site;
	    int *NN = Ymap->numsite;
	    int *MM = Xmap->numsite;

	    if(DEBUG>=2){
	      for(int c = 0; c < colors; c++){
		int Lijy = max(1,align[c+1].Lij1);
		int Lijx = max(1,align[c+1].Lij2);
		int Rijy= min(NN[c],align[c+1].Rij1);
		int Rijx= min(MM[c],align[c+1].Rij2);
		if(DEBUG) assert(Rijx >= Lijx + align[c+1].numpairs - 1);
		if(DEBUG && !(Rijy >= Lijy + align[c+1].numpairs - 1)){
		  #pragma omp critical
		  {
		    printf("Yid=%d(%lld),Xid=%d(%lld),or=%d,c=%d:N=%d,M=%d:alignid=%lu:align[c].numpairs=%d,align[c].score=%0.4f,align[c].logPV=%0.2f:align[c].Rij1=%d,align[c].Lij1=%d,Lijy=%d,Rijy=%d,align[c].repeat=%d,IsRepeatRegion()=%d\n",
			   Yid,YYmap[Yid]->id,Xid,XXmap[Xid]->id,align->orientation,c,NN[c],MM[c],alignid,align[c+1].numpairs,align[c+1].score,align[c+1].logPV,align[c+1].Rij1,align[c+1].Lij1,Lijy,Rijy,align[c].repeat,align->IsRepeatRegion());
		    fflush(stdout);
		    assert(Rijy >= Lijy + align[c+1].numpairs - 1);
		  }
		}
	      }
	    }

	    double trueoffset = Xmap->startloc - Ymap->startloc;// should use endloc to compute trueoffset
	    double trueoverlap = (trueoffset > 0.0) ? min(YY[0][NN[0]+1]-trueoffset, XX[0][MM[0]+1]) : min(XX[0][MM[0]+1]+trueoffset, YY[0][NN[0]+1]);

	    align->trueoffset = (XXmap[0]->startloc == 0.0) ? aNaN : trueoffset;
	    align->trueoverlap = (XXmap[0]->startloc == 0.0) ? aNaN : trueoverlap;

	    align->trueTP = (XXmap[0]->startloc==0.0) ? -1 : 0;
	    if(!(XXmap[0]->startloc==0.0) && trueoverlap >= (XX[0][MM[0]+1] + YY[0][NN[0]+1])*(0.5*MINOVERLAP)){
	      align->trueTP = 1;

	      Tpaircnt++;
	    }

	    if(align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= MINSITES && !align->IsRepeatRegion()){/* good alignment */
	      if(DEBUG>=2){
		for(int c = 0; c < colors; c++){
		  assert(align[c+1].numpairs > 0);
		  assert(align[c+1].sites1[0] >= 1 && align[c+1].sites1[align[c+1].numpairs-1] <= NN[c]);
		  assert(align[c+1].sites2[0] >= 1 && align[c+1].sites2[align[c+1].numpairs-1] <= MM[c]);
		}
	      }

	      if(VERB>=2){
		#pragma omp atomic
		mapcnt++;
	      } else
		Tmapcnt++;

	      if(VERB>=2){
		#pragma omp critical
		{
		  printf("Yid=%lld,Xid=%lld:number of alignments = %d\n",YYmap[Yid]->id,XXmap[Xid]->id,Tmapcnt+mapcnt);
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

	    if(!align || PairSplit)
	      continue;

	    if(DEBUG>=2) assert(align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= MINSITES && !align->IsRepeatRegion());

	    /* collect some more alignment statistics */

	    if(align->score > 0.0)
	      Tposcnt++;

	    for(int c = 0; c < colors; c++){
	      double *Y = YY[c];
	      double *X = XX[c];
	      int N = NN[c];
	      int M = MM[c];
	      int I = align[c+1].sites1[0];
	      int J = align[c+1].sites2[0];
	      int K = align[c+1].numpairs;
	      int Imax = align[c+1].sites1[K-1];
	      int Jmax = align[c+1].sites2[K-1];
	      int Lijy = max(1,align[c+1].Lij1);
	      int Lijx = max(1,align[c+1].Lij2);
	      int Rijy= min(N,align[c+1].Rij1);
	      int Rijx= min(M,align[c+1].Rij2);
	      if(DEBUG) assert(Rijx >= Lijx + align[c+1].numpairs - 1);
	      if(DEBUG) assert(Rijy >= Lijy + align[c+1].numpairs - 1);


	      int XYsites = 0;
	      if(align[c+1].Lend > -2){
		double xy = min(Y[I],align[c+1].orientation ? X[M+1]-X[M+1-J] : X[J]);
		Tlen[c] += 2.0 * xy;
		TlenB[c] += 2.0 * max(0.0, xy - resKB2[c]);
		XYsites += I-Lijy + J-Lijx;
	      }
	      if(align[c+1].Rend > -2){
		double xy = min(Y[N+1]-Y[Imax], align[c+1].orientation ? X[M+1-Jmax] : X[M+1]-X[Jmax]);
		Tlen[c] += 2.0 * xy;
		TlenB[c] += 2.0 * max(0.0, xy - resKB2[c]);
		XYsites += Rijy - Imax + Rijx-Jmax;
	      }
		
	      Tsitecnt[c] += XYsites + 1;
	      Tmisscnt[c] += XYsites;

	      int G = I;
	      int H = J;
	      for(int T = 1; T < K; T++, G= I, H= J){

		I = align[c+1].sites1[T];
		J = align[c+1].sites2[T];

		FLOAT y = Y[I] - Y[G];
		FLOAT x = (align[c+1].orientation==0) ? X[J] - X[H] : X[M+1-H] - X[M+1-J];
		int m = J - G;
		int n = I - H;
		
		if(Poutlier > 0.0){
		  Toutliertot[c]++;
		  PFLOAT Bias,Pen;
		  SintDetail(x,y,m,n,J,I,Bias,Pen,c);
		  double val = ExpOutlierPenalty/(ExpOutlierPenalty+exp(Pen));
		  if(DEBUG) assert(isfinite(val));
		  Toutliersum[c] += val;
		  if(DEBUG) assert(isfinite(Toutliersum[c]));
		  if(!OUTLIER_DELTA(x-y) || Pen < OutlierPenalty){
		    Toutliercnt[c]++;
		    continue;// don't include this interval in error parameter estimation
		  }
		}
		
		Tlen[c] += x+y;
		TlenB[c] += max(ZERO,x - resKB2[c]) + max(ZERO,y - resKB2[c]);
		Tsitecnt[c] += n+m;
		Tmisscnt[c] += n+m-2;
		
		FLOAT err = x-y;
		FLOAT var = varSF2[c] + varSD[c]*(x+y);
		if(QUADRATIC_VARIANCE)
		  var += varSR[c]*(x*x+y*y);
		Terrnorm[c] += err*err/var;
		Tsegcnt[c]++;
		
		if(RefRepeats > 0){  /* append error information of current interval to errors[0..numerrors-1] */
                  #pragma omp critical
		  {
		    if(numerrors[c] >= maxerrors[c]){
		      int newmax = maxerrors[c]*2;
		      Cinterval *newerrors = new Cinterval[newmax];
		      for(int i= 0; i < numerrors[c]; i++)
			newerrors[i] = errors[c][i];
		      delete [] errors[c];
		      errors[c] = newerrors;
		      maxerrors[c] = newmax;
		    }
		    
		    Cinterval *perr = &errors[c][numerrors[c]++];
		    perr->x = x;
		    perr->y = y;
		    if(DEBUG) assert(x > 0.0 && y > 0.0);
		    
		    /* append misaligned sites information of current interval to misaligns[c][0..nummisaligns[c]-1] */
		    if(nummisaligns[c] >= maxmisaligns[c]-1){
		      int newmax = maxmisaligns[c]*2;
		      Cmisalign *newmisaligns = new Cmisalign[newmax];
		      for(int i=0; i < nummisaligns[c];i++)
			newmisaligns[i] = misaligns[c][i];
		      delete [] misaligns[c];
		      misaligns[c] = newmisaligns;
		      maxmisaligns[c] = newmax;
		    }
		    
		    Cmisalign *pmis = &misaligns[c][nummisaligns[c]];
		    pmis[0].x = y;
		    pmis[0].m = (I-G)-1;
		    pmis[1].x = x;
		    pmis[1].m = (J-H)-1;
		    nummisaligns[c] += 2;
		  }// omp critical
		} // if(RefRepeats > 0)
	      } //  T = 0 .. K-1
	    } // c = 0 .. colors-1
	  } // if (align && align->numpairs > 0 ... )
	} /* omp for PairId = 0 .. NumPair -1 */

	#pragma omp critical
	{
	  mapcnt += Tmapcnt;
	  paircnt += Tpaircnt;
	  tcnt += Ttcnt;
	  poscnt += Tposcnt;
	  for(int c = 0; c < colors; c++){
	    sitecnt[c] += Tsitecnt[c];
	    misscnt[c] += Tmisscnt[c];
	    outliertot[c] += Toutliertot[c];
	    outliercnt[c] += Toutliercnt[c];
	    len[c] += Tlen[c];
	    lenB[c] += TlenB[c];
	    outliersum[c] += Toutliersum[c];
	    segcnt[c] += Tsegcnt[c];
	    errnorm[c] += Terrnorm[c];
	  }
	}
      } // end of parallel section

      if(VERB & (VERB>=2 || HASH_DEBUG) && numthreads > 1){
	printf("%d:Completed parallel section, numaligns=%lld,NumPair1=%d,Yid=%d: cpu time=%0.6f, wall time=%0.6f\n", wcnt,(long long int)numaligns,NumPair1,Yid,mtime(), wtime());
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
    
    delete [] YidList1;
    delete [] XidList1;
    delete [] phashList1;
    if(HASH_STREAM>=2){
      delete [] YidList2;
      delete [] XidList2;
      delete [] phashList2;
    }
    if(HashFP)
      hash_close(HashFP,1);

    if(VERB/* HERE HERE >=2 && !PairSplit*/){
      if(FirstAlignments > 0)
	printf("alignments found=%d (terminated after %d alignments or more)\n", mapcnt,FirstAlignments);
      else
	printf("%d:Total alignments found=%d: wall time=%0.6f\n", wcnt, mapcnt, wtime());
      if(VERB/* HERE HERE >=2 */ && RepeatRec && RepeatMaxShift > 0 && RepeatPvalueRatio > 0.0)
	printf("Alignments suppressed due to -RepeatMask %d %0.2e -RepeatRec %0.2f %0.2f %0.2f : %d\n", RepeatMaxShift, RepeatPvalueRatio, RepeatLogPvalueRatio, RepeatKbShiftRatio,RepeatKbShiftRatio2,repeatcnt);
     fflush(stdout);
    }

    if(!PairSplit)
      break;
    
    if(FirstAlignments > 0)
      break;

    printf("-pairsplit not yet implemented for 2-color data\n");
    exit(1);
#if 0 // HERE HERE : convert PairSplit code to 2-color

    /* scan all alignments for highest scoring local alignment : exclude any that are below threshold OR (pairs already confirmed with chimpair = -1, or below threshold/discarded with chimpair = 1) */
    /* also compress the alignment[] array to remove all NULL entries to save memory */
    size_t B = numaligns_compressed, aligncnt = prev_aligncnt;// numaligns_compressed in number of alignments previously checked for NULL entries or below threshold etc */
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
      if(!(align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= MINSITES)){
	align->chimpair = 1;
	continue;
      }
      if(VERB/* HERE HERE >=2 */) aligncnt++;
    }

    if(VERB/* HERE HERE >=2 */){
      printf("%d:Total alignments found=%d:numX=%d,numY=%d,numaligns=%lld->%lld(good=%lld)\n", wcnt, mapcnt,numX,numY,(long long int)numaligns,(long long int)B,(long long int)aligncnt);
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

    #pragma omp parallel for schedule(dynamic,1) num_threads(numthreads) if(numthreads > 1 && PairsplitSeparateQueries && !LVERB)
    for(int rootXid = 0; rootXid < orignumX; rootXid += rootXidInc){

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

	if(DEBUG>=2) assert(align->score > ScoreThreshold && align->logPV > LogPvThreshold && align->numpairs >= MINSITES);

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
	  if(align->Lend <= -1 && (align->refstart < I || align->start < J)){/* flag the left end as if it was deleted (even if MINSITES or map->paired prevents it) to make
										sure the center is cloned and the alignment rechecked/extended */
	    Yleft = I;
	    Xleft = J;
	  }
	  if(align->Rend <= -1 && (align->refend > Imax || align->end > Jmax)){/* flag the right end as if it was deleted (even if MINSITES or map->paired prevents it) to make
										  sure the center is cloned and the alignment rechecked/extended */
	    Yright = Imax;
	    Xright = Jmax;
	  }
	}

	if(align->Lend <= -1 && I-1 >= MINSITES && !Ymap->paired){/* save copy of left end of Ymap up to site I-1 and append it to end of YYmaps[]/refmap[] */
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
	    YYmap[numY++] = newmap;
	    if(DEBUG) assert(numY == numrefmaps);
	  }
	  newmap->allfree();
	  newmap->centercloned = 0;
	  newmap->paired = 0;
	  newmap->origmap = Ymap;
	  if(PairsplitSeparateQueries)
	    newmap->Xid = Ymap->Xid;
	  newmap->id = Ymap->id;
	  newmap->fileid = Ymap->fileid;
	  newmap->numsite[0] = Iend - 1;
	  newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
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

	if(align->Lend <= -1 && J-1 >= MINSITES && !Xmap->paired){/* save copy of left end of Xmap up to site J-1 and append it to end of XXmaps[]/map[] */
	  Xleft = J;
	  int Jend = (PAIRSPLIT_OVERLAP && PAIRSPLIT_EXTEND && align->Lend <= -1 && (align->refstart < I || align->start < J)) ? min(M,J+1) : J;/* extend map to site Jend-1 */

	  Cmap *newmap = NULL;

          #pragma omp critical(nummaps)
	  {
	    maxmapalloc(nummaps+1,maxmaps,map,0,1);
	    XXmap = &map[Cfirst];/* in case map[] was reallocated */
	    newmap = map[nummaps++];
	    if(!newmap)
	      map[nummaps-1] = newmap = new Cmap;
	    newmap->mapid = numX;
	    XXmap[numX++] = newmap;
	    if(DEBUG) assert(numX == nummaps - Cfirst);
	  }

	  newmap->allfree();
	  newmap->centercloned = 0;
	  newmap->paired = 0;
	  newmap->origmap = Xmap;
	  if(PairsplitSeparateQueries)
	    newmap->Xid = Xmap->Xid;
	  newmap->id = Xmap->id;
	  newmap->fileid = Xmap->fileid;
	  newmap->numsite[0] = Jend - 1;
	  newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
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

	if(align->Rend <= -1 && N-Imax >= MINSITES && !Ymap->paired && PairSplit > 0){/* save copy of right end of Ymap starting at site Imax+1 and append it to end of YYmaps[]/refmap[] */
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
	    YYmap[numY++] = newmap;
	    if(DEBUG) assert(numY == numrefmaps);
	  }
	  newmap->allfree();
	  newmap->centercloned = 0;
	  newmap->paired = 0;
	  newmap->origmap = Ymap;
	  if(PairsplitSeparateQueries)
	    newmap->Xid = Ymap->Xid;
	  newmap->id = Ymap->id;
	  newmap->fileid = Ymap->fileid;
	  newmap->numsite[0] = N - Istart;
	  newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
	  newmap->site[0][0] = 0.0;
	  newmap->left[0] = Istart + 1;
	  if(DEBUG) assert(newmap->left >= 0);
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

	if(align->Rend <= -1 && M-Jmax >= MINSITES && !Xmap->paired){/* save copy of right end of Xmap starting at site Jmax+1 and append it to end of XXmaps[]/map[] */
	  Xright = Jmax;
	  int Jstart = (PAIRSPLIT_OVERLAP && PAIRSPLIT_EXTEND && align->Rend <= -1 && (align->refend > Imax || align->end > Jmax)) ? max(1,Jmax-1) : Jmax;/* start at site Jstart+1 */
	  Cmap *newmap = NULL;
        

          #pragma omp critical(nummaps)
	  {
	    maxmapalloc(nummaps+1,maxmaps,map,0,1);
	    XXmap = &map[Cfirst];/* in case map[] was reallocated */
	    newmap = map[nummaps++];
	    if(!newmap)
	      map[nummaps-1] = newmap = new Cmap;
	    newmap->mapid = numX;
	    XXmap[numX++] = newmap;
	    if(DEBUG) assert(numX == nummaps - Cfirst);
	  }
	  newmap->allfree();
	  newmap->centercloned = 0;
	  newmap->paired = 0;
	  newmap->origmap = Xmap;
	  if(PairsplitSeparateQueries)
	    newmap->Xid = Xmap->Xid;
	  newmap->id = Xmap->id;
	  newmap->fileid = Xmap->fileid;
	  newmap->numsite[0] = M - Jstart;
	  newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
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
	Xmap->paired = Yid+1;
	Ymap->paired = Xid+1;
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
	      int JR = palign->sites2[align->numpairs-1];
	      printf("supressing alignment[%lld](score=%0.6f,logPV=%0.2f,numpairs=%d) between Yid=%d(id=%lld) and Xid=%d(id=%lld,orientation=%d) (Xrange = %d to %d) due to cloning of Xid (orientation=%d) from Xleft=%d to Xright=%d, M=%d\n",
		     (long long)A,palign->score,palign->logPV,palign->numpairs,palign->mapid1,YYmap[palign->mapid1]->id,Xid,XXmap[Xid]->id, palign->orientation, JL, JR, align->orientation, Xleft, Xright, M);
	      fflush(stdout);
	    }
	    palign->chimpair = 1;
	    continue;
	  }
	}

	/* Clone center regions */
	align->chimpair = 1;/* this alignment will be found again with the Cloned center region(s) */
	LocalCnt--;
	if(Yleft > 0 || Yright <= N){
	  if(DEBUG>=2) assert(min(N,Yright)-max(1,Yleft) < N-1);/* make sure center cloned region is smaller than original */

	  if(PAIRSPLIT_EXTEND){/* limit extensions to no more than MINSITES-1 sites */
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

	    Yleft = min(Yleft,max(Yleft-MINSITES+1,align->refstart));
	    Yright = max(Yright,min(Yright+MINSITES-1,align->refend));

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
	    YYmap[numY++] = newmap;
	    if(DEBUG) assert(numY == numrefmaps);
	  }
	  newmap->allfree();
	  newmap->centercloned = 0;
	  newmap->origmap = Ymap;
	  if(PairsplitSeparateQueries)
	    newmap->Xid = Ymap->Xid;
	  newmap->id = Ymap->id;
	  newmap->fileid = Ymap->fileid;
	  newmap->numsite[0] = min(N,Yright)-max(1,Yleft)+1;
	  newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
	  newmap->site[0][0] = 0.0;
	  newmap->left[0] = Yleft;
	  newmap->right[0] = Yright;

	  double YLextend = MININTERVAL;
	  double YRextend = MININTERVAL;
	  if(OUTLIER_FIXES && align->Lend == -1 && Yleft > 0){/* make sure left end of new map is extended to at least overlap the left end of X */
	    if(DEBUG && !PAIRSPLIT_EXTEND && MINSITES >= 8) assert(Xleft <= 0);/* If Xleft > 0 && Yleft > 0, then both X and Y were trimmed on the left, hence there should be an endoutlier with align->Lend <= -2 */
	    int IL = align->sites1[0];
	    int JL = align->sites2[0];
	    double Yend = Y[IL] - (align->orientation ? X[M+1] - X[M+1-JL] : X[JL]) - 0.001;
	    YLextend = max(YLextend, Y[Yleft] - Yend);
	  }
	  if(OUTLIER_FIXES && align->Rend == -1 && Yright <= N){/* make sure right end is extended to at least overlap the right end of X */
	    if(DEBUG && !PAIRSPLIT_EXTEND && MINSITES >= 8 && !(Xright > M)){/* If Xright <= M && Yright <= N, then both X and Y were trimmed on the right, hence there should be an endoutlier with align->Rend  <= -2 */
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
	  /* make sure center cloned region is smaller than the original */
	  if(DEBUG>=2 && !(min(M,Xright)-max(1,Xleft) < M-1)){
	    printf("Yid=%d,Xid=%d:Xleft=%d,Xright=%d,M=%d\n",Yid,Xid,Xleft,Xright,M);
	    fflush(stdout);
	    assert(Xright-Xleft < M-1);
	  }

	  if(PAIRSPLIT_EXTEND){/* limit extensions to no more than MINSITES-1 sites */
	    if(DEBUG) assert(align->start >= 0);
	    if(DEBUG) assert(align->end <= M+1);
	    int origXleft = Xleft;
	    int origXright = Xright;
	    Xleft = min(Xleft,max(Xleft-MINSITES+1,align->start));
	    Xright = max(Xright,min(Xright+MINSITES-1,align->end));

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
	    maxmapalloc(nummaps+1,maxmaps,map,0,1);
	    XXmap = &map[Cfirst];/* in case map[] was reallocated */
	    newmap = map[nummaps++];
	    if(!newmap)
	      map[nummaps-1] = newmap = new Cmap;
	    newmap->mapid = numX;
	    XXmap[numX++] = newmap;
	    if(DEBUG && !(numX == nummaps - Cfirst)){
	      printf("numX=%d,nummaps=%d,Cfirst=%d\n",numX,nummaps,Cfirst);
	      fflush(stdout);
	      assert(numX == nummaps - Cfirst);
	    }
	  }
	  newmap->allfree();
	  newmap->centercloned = 0;
	  newmap->origmap = Xmap;
	  if(PairsplitSeparateQueries)
	    newmap->Xid = Xmap->Xid;
	  newmap->id = Xmap->id;
	  newmap->fileid = Xmap->fileid;
	  newmap->numsite[0] = min(M,Xright)-max(1,Xleft)+1;
	  newmap->site[0] = new FLOAT[newmap->numsite[0]+2];
	  newmap->site[0][0] = 0.0;
	  double XLextend = MININTERVAL;
	  double XRextend = MININTERVAL;
	  if(!align->orientation){
	    newmap->left[0] = Xleft;
	    newmap->right[0] = Xright;

	    if(OUTLIER_FIXES && align->Lend == -1 && Xleft > 0){/* make sure left end of new map is extended to at least overlap the left end of Y */
	      if(DEBUG && !PAIRSPLIT_EXTEND && MINSITES >= 8) assert(Yleft <= 0);/* If Yleft > 0 && Xleft > 0, then both X and Y were trimmed on the left, hence there should be an endoutlier with align->Lend <= -2 */
	      int IL = align->sites1[0];
	      int JL = align->sites2[0];
	      double Xend = X[JL] - Y[IL] - 0.001;
	      XLextend = max(XLextend, X[Xleft] - Xend);
	    }
	    if(OUTLIER_FIXES && align->Rend == -1 && Xright <= M){/* make sure right end of new map is extended to at least overlap the right end of Y */
	      if(DEBUG && !PAIRSPLIT_EXTEND && MINSITES >= 8) assert(Yright > N);/* If Yright <= N && Xright <= M, then both X and Y were trimmed on the right, hence there should be an endoutlier with align->Rend <= -2 */
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
	      if(DEBUG && !PAIRSPLIT_EXTEND && MINSITES >= 8) assert(Yleft <= 0);/* If Yleft > 0 && Xleft > 0, then both X(reversed) and Y were trimmed on the left, hence there should be an endoutlier with align->Lend <= -2 */
	      int IL = align->sites1[0];
	      int JL = align->sites2[0];
	      double Xend = (X[M+1] - X[M+1 - JL]) - Y[IL] - 0.001;
	      XLextend = max(XLextend, X[M+1] - X[M+1 - Xleft] - Xend);
	    }
	    if(OUTLIER_FIXES && align->Rend == -1 && Xright <= M){/* make sure right end of new map (reversed) is extended to at least overlap the right end of Y */
	      if(DEBUG && !PAIRSPLIT_EXTEND && MINSITES >= 8) assert(Yright > N);/* If Yright <= N && Xright <= M, then both X(reversed) and Y were trimmed on the right, hence there should be an endoutlier with align->Rend <= -2 */
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
    if(VERB>=2 && PairSplit){
      printf("After splitting maps of %lu best alignment for %d original Query id's: numX=%d,numY=%d : wall time=%0.6f\n", splitcnt, orignumX,numX,numY,wtime());
      fflush(stdout);
    }
    if(!splitcnt)
      break;
#endif
  }// while(1)

#if 0 // HERE HERE : convert -pairsplit code to 2-color 
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

      /* sort alignments in order of id2, id1, orientation, sites1[0],sites1[U-1] (null alignments last) : this way neighboring alignments will be next to each other */
      qsort(&alignment[0],numaligns,sizeof(Calign *), (intcmp*)IdOrSiteInc);
      
      double threshold = -log(1.0/Psplit);

      for(size_t i = 1; i < numaligns; i++){
	Calign *p = alignment[i-1];
	Calign *q = alignment[i];
	if(VERB >= (LVERB?1:2)){
	  printf("i=%lu:p=alignment[i-1]=%p, q=alignment[i]=%p\n",i,p,q);
	  fflush(stdout);
	}

	if(!p || !q)
	  continue;

	int rid = p->mapid1;
	int mid = p->mapid2;
	int rid2 = q->mapid1;
	int mid2 = q->mapid2;

	if(VERB >= (LVERB?1:2)){
	  printf("i=%lu: rid=%d(%lld),rid2=%d(%lld),mid=%d(%lld),mid2=%d(%lld),or=%d,or2=%d\n",i, rid, YYmap[rid]->id, rid2, YYmap[rid2]->id, mid, XXmap[mid]->id, mid2, XXmap[mid2]->id, p->orientation, q->orientation);
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
	  printf("Yid=%lld,Xid=%lld,or=%d:N=%d,M=%d,N2=%d,M2=%d:LI=%d(%d),RI=%d(%d),LJ=%d(%d),RJ=%d(%d):LI2=%d(%d),RI2=%d(%d),LJ2=%d(%d),RJ2=%d(%d)\n",
		 Ymap->id,Xmap->id,p->orientation,N,M,N2,M2,LI,refstart,RI,refend,LJ,qrystart,RJ,qryend,LI2,refstart2,RI2,refend2,LJ2,qrystart2,RJ2,qryend2);
	  fflush(stdout);
	}

	if(refend < refstart2){/* p is left of q */
	  if(VERB >= (LVERB ? 1 : 2)){
	    printf("  Neighboring alignments with gap\n");
	    fflush(stdout);
	  }
	  /* check if the gap is a valid alignment interval (NOT an outlier) */
	  PFLOAT x = fabs(origX2[qrystart2]-origX[qryend]);
	  PFLOAT y = origY2[refstart2]-origY[refend];
	  SintDetail(x, y, abs(qrystart2-qryend), refstart2-refend, qrystart2, refstart2, Bias, Pen);
	  PFLOAT outscore = Bias + Pen;
	  PFLOAT nthreshold = threshold;
	  if(outlierBC)
	    nthreshold -= OutlierPenaltyBC[refstart2-refend] + OutlierPenaltyBC[qrystart2-qryend];
	  if(outscore < nthreshold) /* An outlier at -pairsplit P-value */
	    continue;
	  PFLOAT iscore = Bias + OutlierBias + (OUTLIER_DELTA(x-y) ? max(Pen,OutlierPenalty) : Pen);

	  /* merge alignment q and p */
	  Calign *r = new Calign[colors==1 ? 1 : colors+1];
	  r->expand_arrays(U+U2);
	  r->numpairs = U + U2;
	  r->mapid1 = mapid1;
	  r->mapid2 = mapid2;
	  r->orientation = p->orientation;
	  if(DEBUG) assert(YYmap[mapid1]->numsite[0] == origN);
	  if(DEBUG) assert(XXmap[mapid2]->numsite[0] == origM);

	  if(VERB >= (LVERB ? 1 : 2)){
	    printf("p is left of q and gap is a valid alignment interval : combining into one alignment (p->numpairs=%d,q->numpairs=%d,r->numpairs=%d)\n",
		   p->numpairs,q->numpairs,r->numpairs);
	    fflush(stdout);
	  }

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
	  if(Xshift || Yshift){
	    r->Lend = localtype;
	    r->Lij1 = r->sites1[0];
	    r->Lij2 = r->sites2[0];
	  } else {
	    r->Lend = (Xshift || Yshift) ? localtype : p->Lend;
	    r->Lij1 = Yshift + p->Lij1;
	    r->Lij2 = (p->orientation ? origM - Xshift - M : Xshift) + p->Lij2;
	  }
	  if(Xshift2 + M2 < origM || Yshift2 + N2 < origN){
	    r->Rend = localtype;
	    r->Rij1 = r->sites1[r->numpairs - 1];
	    r->Rij2 = r->sites2[r->numpairs - 1];
	  } else {
	    r->Rend = q->Rend;
	    r->Rij1 = Yshift2 + q->Rij1;
	    r->Rij2 = (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->Rij2;
	  }
	  r->score = 0;
	  for(int t = 0; t <= U + U2; t++)
	    r->score += r->iscore[t];
	  r->logPV = alignFP(r, origY, origXrev, origN, origM, p->orientation, mapid1, mapid2, r->score, 0);
	  //	  r->logPV = p->logPV + q->logPV;

	  delete [] p;
	  delete [] q;
	  alignment[i-1] = r;
	  alignment[i] = 0;
	} else if(refend2 < refstart){/* p is right of q */
	  if(VERB >= (LVERB ? 1 : 2)){
	    printf("  Neighboring alignments with gap\n");
	    fflush(stdout);
	  }

	  /* check if the gap is a valid alignment interval (NOT an outlier) */
	  PFLOAT x = fabs(origX[qrystart]-origX2[qryend2]);
	  PFLOAT y =  origY[refstart]-origY2[refend2];
	  SintDetail(x, y, abs(qrystart - qryend2), refstart - refend2, qrystart, refstart, Bias, Pen);
	  FLOAT outscore = Bias + Pen;
	  PFLOAT nthreshold = threshold;
	  if(outlierBC)
	    nthreshold -= OutlierPenaltyBC[refstart - refend2] + OutlierPenaltyBC[qrystart - qryend2];
	  if(outscore < nthreshold) /* An outlier at -pairsplit P-value */
	    continue;
	  FLOAT iscore = Bias + OutlierBias + (OUTLIER_DELTA(x-y) ? max(Pen,OutlierPenalty) : Pen);

	  /* merge alignment q and p */	  
	  Calign *r = new Calign[colors==1 ? 1 : colors+1];
	  r->expand_arrays(U+U2);
	  r->numpairs = U + U2;
	  r->mapid1 = mapid1;
	  r->mapid2 = mapid2;
	  r->orientation = p->orientation;

	  if(VERB >= (LVERB ? 1 : 2)){
	    printf("p is right of q and gap is a valid alignment interval : combining into one alignment (p->numpairs=%d,q->numpairs=%d,r->numpairs=%d)\n",
		   p->numpairs,q->numpairs,r->numpairs);
	    fflush(stdout);
	  }

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
	  if(Xshift2 || Yshift2){
	    r->Lend = localtype;
	    r->Lij1 = r->sites1[0];
	    r->Lij2 = r->sites2[0];
	  } else {
	    r->Lend = q->Lend;
	    r->Lij1 = Yshift2 + q->Lij1;
	    r->Lij2 = (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->Lij2;
	  }
	  if(Xshift + M < origM || Yshift + N < origN){
	    r->Rend = localtype;
	    r->Rij1 = r->sites1[r->numpairs - 1];
	    r->Rij2 = r->sites2[r->numpairs - 1];
	  } else{
	    r->Rend = p->Rend;
	    r->Rij1 = Yshift + p->Rij1;
	    r->Rij2 = (p->orientation ? origM - Xshift - M : Xshift) + p->Rij2;
	  }
	  r->score = 0;
	  for(int t = 0; t <= U + U2; t++)
	    r->score += r->iscore[t];
	  r->logPV = alignFP(r, origY, origXrev, origN, origM, p->orientation, mapid1, mapid2, r->score, 0);
	  //	  r->logPV = p->logPV + q->logPV;

	  delete [] p;
	  delete [] q;
	  alignment[i-1] = r;
	  alignment[i] = 0;

	} else {/* overlapped or adjoining with shared site in Y */
	  if(VERB >= (LVERB ? 1 : 2)){
	    printf("  Overlapped or Adjoining alignments\n");
	    fflush(stdout);
	  }
	  /* first check if one alignment completely overlaps the other */
	  if(refstart <= refstart2 && refend2 <= refend  &&
	     (p->orientation ? (qryend <= qryend2 && qrystart2 <= qrystart) : (qrystart <= qrystart2 && qryend2 <= qryend))){/* p overlapped by q */
	    if(VERB >= (LVERB ? 1 : 2)){
	      printf("Complete Overlap of q (deleting)\n");
	      fflush(stdout);
	    }
	    delete [] q;
	    alignment[i] = 0;
	  } else if(refstart2 <= refstart && refend <= refend2 &&
		    (p->orientation ? (qryend2 <= qryend && qrystart <= qrystart2) : (qrystart2 <= qrystart && qryend <= qryend2))){/* q overlapped by p */
	    if(VERB >= (LVERB ? 1 : 2)){
	      printf("Complete Overlap of p (deleting)\n");
	      fflush(stdout);
	    }
	    delete [] p;
	    alignment[i-1] = 0;
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

	    if(Xshift || Yshift){
	      r->Lend = localtype;
	      r->Lij1 = r->sites1[0];
	      r->Lij2 = r->sites2[0];
	    } else {
	      r->Lend = p->Lend;
	      r->Lij1 = Yshift + p->Lij1;
	      r->Lij2 = (p->orientation ? origM - Xshift - M : Xshift) + p->Lij2;
	    }
	    if(Xshift2 + M2 < origM || Yshift2 + N2 < origN){
	      r->Rend = localtype;
	      r->Rij1 = r->sites1[r->numpairs - 1];
	      r->Rij2 = r->sites2[r->numpairs - 1];
	    } else {
	      r->Rend = q->Rend;
	      r->Rij1 = Yshift2 + q->Rij1;
	      r->Rij2 = (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->Rij2;
	    }

	    r->score = 0;
	    for(int t = 0; t <= U + U2; t++)
	      r->score += r->iscore[t];
	    r->logPV = alignFP(r, origY, origXrev, origN, origM, p->orientation, mapid1, mapid2, r->score, 0);
	    
	    delete [] p;
	    delete [] q;
	    alignment[i-1] = r;
	    alignment[i] = 0;
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
	    r->orientation = p->orientation;

	    if(VERB>=2){
	      printf("q overlaps p from left and overlapped region matches : combining into one alignment (t1=%d,p->numpairs=%d,q->numpairs=%d,r->numpairs=%d)\n",
		     t1, p->numpairs,q->numpairs,r->numpairs);
	      fflush(stdout);
	    }

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

	    if(Xshift2 || Yshift2){
	      r->Lend = localtype;
	      r->Lij1 = r->sites1[0];
	      r->Lij2 = r->sites2[0];
	    } else {
	      r->Lend = q->Lend;
	      r->Lij1 = Yshift2 + q->Lij1;
	      r->Lij2 = (q->orientation ? origM - Xshift2 - M2 : Xshift2) + q->Lij2;
	    }
	    if(Xshift + M < origM || Yshift + N < origN){
	      r->Rend = localtype;
	      r->Rij1 = r->sites1[r->numpairs - 1];
	      r->Rij2 = r->sites2[r->numpairs - 1];
	    } else {
	      r->Rend = p->Rend;
	      r->Rij1 = Yshift + p->Rij1;
	      r->Rij2 = (p->orientation ? origM - Xshift - M : Xshift) + p->Rij2;
	    }

	    r->score = 0;
	    for(int t = 0; t <= U + U2; t++)
	      r->score += r->iscore[t];
	    r->logPV = alignFP(r, origY, origXrev, origN, origM, p->orientation, mapid1, mapid2, r->score, 0);
	    
	    delete [] p;
	    delete [] q;
	    alignment[i-1] = r;
	    alignment[i] = 0;
	  }
	}
      }
      delete [] origXrevMem;
    }
  }
#endif // PairSplit code

  if(VERB){/* display statistics */
    printf("maps=%d:TP pairs=%d/%d, total pairs=%lld, positive scores=%d:\n",
	   nummaps,paircnt,mapcnt,tcnt,poscnt);
    for(int c = 0; c < colors; c++){
      printf("  c=%d:overlaps:len=%0.3f, sites=%d, misscnt=%d (expected %0.1f)\n", 
	     c,len[c], sitecnt[c], misscnt[c], 2.0*len[c]*(FP[c]*0.01 + FN[c]/Ytheta[c]));
      if(segcnt[c] > 0)
	printf("  segcnt=%d,norm=%0.4f\n",
	     segcnt[c],errnorm[c]/segcnt[c]);
      if(Poutlier > 0.0)
	printf("  outliers=%d/%d (rate=%0.6f)\n",outliercnt[c],outliercnt[c]+segcnt[c],outliersum[c]*1.0/outliertot[c]);
    }
    fflush(stdout);

    if(RefRepeats > 0){ // WARNING : does not support QUADRATIC_VARIANCE correctly
      double resKB = mres*PixelLen*2.0;
      for(int c = 0; c < colors; c++){
	if(numerrors[c] > 0){
	  /* generate updated estimate of sizing error parameters Var(x+y) = A + B(x+y) + R(x^2 + y^2) */
	  double A = SF[c]*SF[c]*2.0;
	  double B = fabs(SD[c])*SD[c];
	  double R = SR[c]*SR[c];
	  for(int iter=0; iter < 30; iter++){
	    double Lsum= 0.0;/* sum of log(A+B(x+y)+R(x^2+y^2)) */
	    double errsum= 0.0;/* sum of (x-y)^2/(A+B(x+y)+R(x^2+y^2)) */
	    Cinterval *perr = errors[c];
	    for(int i = 0;i < numerrors[c];i++,perr++){      
	      double err = perr->y - perr->x;
	      double sum = perr->y + perr->x;
	      err *= err;
	      double var = A+B*sum;
	      if(QUADRATIC_VARIANCE)
		var += R*(perr->x * perr->x + perr->y * perr->y);
	      Lsum += log(var);
	      errsum += err/var;
	    }
	    double delta = errsum/numerrors[c];

	    if(VERB>=2){
	      printf("iter=%d,c=%d:cnt=%d,sf=%0.6f,sd=%0.6f:LL=%0.6f,norm=%0.6f:sf->%0.6f,sd->%0.6f,sr->%0.6f\n",
		     iter,c,numerrors[c],sqrt(A*0.5),sqrt(B),-(Lsum+errsum)*0.5,delta,sqrt(A*delta*0.5),sqrt(B*delta),sqrt(R*delta));
	      fflush(stdout);
	    }

	    B *= delta;
	    A *= delta;
	    R *= delta;

	    /* update estimate of A */
	    double bestLL = MINSCORE;
	    double bestA = A;

	    for(int Aiter=0; Aiter < 3; Aiter++){
	      double Lsum= 0.0;/* sum of log(A+B(x+y)) */
	      double errsum= 0.0;/* sum of (x-y)^2/(A+B(x+y)) */
	      double Isum= 0.0;/* sum of 1/(A+B(x+y)) */
	      double errsum2= 0.0;/* sum of (x-y)^2/(A+B(x+y))^2 */
	      double Isum2= 0.0;/* sum of 1/(A+B(x+y))^2 */
	      double errsum3= 0.0;/* sum of (x-y)^2/(A+B(x+y))^3 */
	      Cinterval *perr = errors[c];
	      for(int i = 0;i < numerrors[c]; i++,perr++){
		double err = perr->y - perr->x;
		double sum = perr->y + perr->x;
		err *= err;
		double var = A+B*sum;
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
		  printf("iter=%d,%d,c=%d:sf=%0.8f,sd=%0.8f,sr=%0.8f:LL=%0.8f,LL1=%0.8f,LL2=%0.8f:bestLL=%0.8f,backtracking (sf->%0.8f)\n",
			 Aiter,iter,c,sqrt(A),copysign(sqrt(fabs(B)),B),sqrt(R),LL, LL1, LL2, bestLL,sqrt(0.5*(A+bestA)));
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
		printf("iter=%d,%d,c=%d:sf=%0.6f,sd=%0.4f:LL=%0.6f,LL1=%0.3f,LL2=%0.3f:delta=%0.6f(sf->%0.6f)\n",
		       Aiter,iter,c,sqrt(A*0.5),sqrt(B),-(Lsum+errsum)*0.5,0.5*(errsum2-Isum),0.5*Isum2-errsum3,
		       delta,sqrt((A+delta)*0.5));
		fflush(stdout);
	      }

	      A += delta;
	      if(QUADRATIC_VARIANCE){
		if(DEBUG) assert(R >= 0.0);
		if(B < 0.0 && R >= 0.0)
		  A = max(A, B*B/(4.0*R));
	      }
	    } // Aiter = 0..2 */
	    A = bestA;
	  } // iter = 0..29
	  if(VERB) {
	    printf("  c=%d:var(y) = sf^2 + sd^2 y + sr^2 y^2 : sf = %0.6f -> %0.6f (kb), sd = %0.6f -> %0.6f (kb^1/2), sr= %0.6f -> %0.6f\n", c, SF[c], sqrt(A*0.5), SD[c], sqrt(B), SR[c], sqrt(R));
	    fflush(stdout);
	  }
	} // if(numerrors[c] > 0)

	if(nummisaligns[c] > 0){ /* plot misaligned cuts vs aligned interval size AND estimate scaling of FP[c] & FN[c] */
	  /* sort aligned intervals by size */
	  qsort(misaligns[c],nummisaligns[c],sizeof(Cmisalign), (intcmp*)xinc);
	  if(VERB>=2)
	    printf("IntervalSize     MisalignedSites    Density    Density(above %0.2fkb) Samples\n", resKB);
	  double msumtot = 0.0,xpsumtot=0.0;
	  for(int bin=0;bin < MISALIGN_BINS;bin++){
	    int jmin = (bin*nummisaligns[c])/MISALIGN_BINS;
	    int jmax = (bin == MISALIGN_BINS-1) ? nummisaligns[c] : ((bin+1)*nummisaligns[c])/MISALIGN_BINS;
	    if(jmax > jmin){
	      double xsum=0.0,xpsum=0.0;
	      int msum = 0;
	      Cmisalign *pmis = &misaligns[c][jmin];
	      for(int j = jmin; j < jmax; j++, pmis++){
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
	    double pFN = FN[c];
	    double pTP = 1.0 - pFN;
	    double F = FP[c]*0.01;
	    double lambda = Ylambda[c];
	    double theta = Ytheta[c];
	    double Itheta = 1.0/(theta - delta - var[c]*0.5);
	    double R = F+pTP*pFN/lambda;
	    double newR = msumtot/xpsumtot;
	    double newF = F*newR/(F+(Itheta-F)*pFN);
	    double newpTP = (Itheta - newR)/(Itheta - newF);

	    printf("  c=%d:R= %0.6f -> %0.6f, FP= %0.6f -> %0.6f, FN= %0.6f -> %0.6f \n", c, R, newR, F*100.0, newF*100.0, FN[c], 1.0 - newpTP);
	    fflush(stdout);

	    if(DEBUG) assert(fabs(newF - F*newR/R) <= newF*0.001);
	  }
	}// if(nummisaligns[c] > 0)
      } // c = 0 .. colors-1
    } // if(RefRepeats > 0)
  } // if(VERB)

  int orignummaps = nummaps;
  
  if(PairMerge){/* reset paired flag for all input maps */
    for(int i=0; i < numX; i++)
      XXmap[i]->paired = 0;
    for(int i=0; i < numY; i++)
      YYmap[i]->paired = 0;
  }

  if(RepeatShift > 0.0){
    RefIndex = num_files;
    QueryIndex = 0;
    output_xmap2(output_prefix, 0, numaligns);/* output .xmap file */
  }

  char prefix[PATH_MAX];
  if(PairMerge && PairMergeRepeat)
    sprintf(prefix,"%sM%d",output_prefix,pairmergeIter);
  else
    strcpy(prefix,output_prefix);
  if(PairSplit || (PairMerge && PairMergeXmap)) {
    if(mres > 0.001 || mresSD > 0.0 || fabs(PixelLen - 0.500) > 1e-12){
      /* create reduced resolution Reference and Query map files : <prefix>_r.cmap and <prefix>_q.cmap */
      char cmapprefix[PATH_MAX];
      sprintf(cmapprefix,"%s_r",prefix);
      output_cmap(cmapprefix,YYmap,0,orignumY-1);

      sprintf(cmapprefix,"%s_q",prefix);
      output_cmap(cmapprefix,XXmap,0,orignumX-1);
    }
    output_xmap2(prefix, 0 , numaligns);/* output .xmap file */
    if(PairSplit && dosvdetect)
      output_smap(output_prefix, vfixx_filename[0]);/* output .smap file */
  }

  if(VERB/* >=2 */){
    printf("pairalign:numX=%d,numY=%d,nummaps=%d,numaligns=%llu\n",
	   numX,numY,nummaps,(unsigned long long)numaligns);
    fflush(stdout);
  }

  if(PairMerge)
    pairmerge(pairmergeIter);/* generates merged maps at map[nummaps...] and tags merged maps with paired=1) */

  if(!PairSplit && !PairMerge)
    output_align(output_prefix,0);/* output .align file */

  int PairMergeProgress = 0;
  if(PairMerge){
    if(!(PairMergeRepeat && nummaps > orignummaps)){/* output merged and unmerged maps */
      int linkcnt = 0;/* number of files linked to previous copy instead of being written */
      double origwt = wtime();

      #pragma omp parallel for schedule(dynamic,1) num_threads(gnumthreads)
      for(int i = 0; i < nummaps; i++){
        Cmap *pmap = Gmap[i];
	if(DEBUG>=2 && i >= orignummaps) assert(!pmap->paired);
	if(pmap->origmap || pmap->paired)
	  continue;

	char prefix[PATH_MAX];
	sprintf(prefix,"%s_contig%lld",output_prefix, pmap->id);

	if(i < orignummaps && !pmap->centercloned){/* unchanged contig map */

#ifndef WIN32
	  if(SplitMap && !strstr(prefix,"/dev/null")){ /* create a link to the original file */
            if(DEBUG) assert(pmap->fileid < num_files);
	    char filename[PATH_MAX];
	    sprintf(filename,"%s.cmap",prefix);

	    FILE *fp = fopen(filename,"r");
	    if(ForceOverwrite && fp != NULL){
	      printf("WARNING:Output file %s already exists\n",filename);
	      fflush(stdout);
	      fclose(fp);

	      /* unlink the filename */
	      errno = 0;
	      if(unlink(filename)){
                int eno = errno;
		char *err = strerror(eno);
		printf("unlink(%s) failed:errno=%d:%s\n",filename,eno,err);
		printf("Unable to write new version of %s\n",filename);
		exit(1);
              }
            }

	    errno = 0;
	    if(!link(vfixx_filename[pmap->fileid],filename)){
              if(VERB>=2){
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
	      if(pmap->origsite[c] == NULL)
		continue;/* no original map resolution : must have already been restored */

	      /* free up regular arrays */
	      if(!pmap->blockmem)
		delete [] pmap->site[c]; 
	      pmap->site[c] = NULL;
	      if(pmap->siteSD[c]) delete [] pmap->siteSD[c];
	      if(pmap->sitecov[c]) delete [] pmap->sitecov[c];
	      if(pmap->sitecnt[c]) delete [] pmap->sitecnt[c];
	      if(pmap->SNRcnt[c]){
		if(pmap->SNRgmean[c]) delete [] pmap->SNRgmean[c];
		if(pmap->lnSNRsd[c]) delete [] pmap->lnSNRsd[c];
		if(pmap->SNRdist[c]){
		  for(int t = 0; t <= pmap->numsite[c]+1; t++)
		    if(pmap->SNRdist[c][t]){
		      if(VERB>=3){
			printf("Calling delete[] Gmap[%d]->SNRdist[%d][t=%d]=%p -> 0\n",i,c,t,pmap->SNRdist[c][t]);
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
	      pmap->blockmem = 0;
	      pmap->site[c] = pmap->origsite[c];
	      pmap->siteSD[c] = pmap->origsiteSD[c];
	      pmap->sitecov[c] = pmap->origsitecov[c];
	      pmap->sitecnt[c] = pmap->origsitecnt[c];
	      pmap->SNRcnt[c] = pmap->origSNRcnt[c];
	      pmap->SNRdist[c] = pmap->origSNRdist[c];

	      pmap->SNRgmean[c] = pmap->origSNRgmean[c];
	      pmap->lnSNRsd[c] = pmap->origlnSNRsd[c];
	      if(pmap->SNRcnt[c]){
		for(int t = 0; t <= pmap->orignumsite[c]+1; t++){
		  if(VERB>=3){
		    printf("map[%d]->SNRdist[%c][%d]= %p -> %p\n",i,c,t,pmap->SNRdist[c][t],pmap->origSNRdist[c][t]);
		    fflush(stdout);
		  }

		  pmap->SNRdist[c][t] = pmap->origSNRdist[c][t];
		}
	      }

	      /* reset backup pointers */
	      pmap->origsite[c] = NULL;
	      pmap->origsiteSD[c] = NULL;
	      pmap->origsitecov[c] = NULL;
	      pmap->origsitecnt[c] = NULL;
	      if(pmap->origSNRcnt[c]){
		pmap->origSNRcnt[c] = NULL;
		pmap->origSNRdist[c] = NULL;
		pmap->origSNRgmean[c] = NULL;
		pmap->origlnSNRsd[c] = NULL;
	      }
	    }
	  }
	}

        if(VERB>=2){
          #pragma omp critical
	  {
	    if(i < orignummaps && !pmap->centercloned)
	      printf("writing out original mapid=%d,id=%lld,len=%0.3f as %s.cmap\n",i,pmap->id,pmap->site[0][pmap->numsite[0]+1],prefix);
	    else
	      printf("writing out merged mapid=%d,id=%lld as %s.cmap\n",i,pmap->id,prefix);
	    fflush(stdout);
          }
	}

	output_cmap(prefix, Gmap, i, i);
      }    

      if(VERB){
        double wt = wtime();
	printf("pairmerge: iter=%d: %d Files written (including %d linked): wall time=%0.6f (cum=%0.6f)\n", pairmergeIter, nummaps, linkcnt, wt - origwt, wt);
	fflush(stdout);
      }

#if 0
      for(int i = orignummaps; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	if(pmap->origmap)
	  continue;
	if(DEBUG) assert(!pmap->paired);
	sprintf(filename,"%s_contig%lld",output_prefix, pmap->id);

	if(VERB/* >= 2*/){
	  printf("writing out merged mapid=%d,id=%lld as %s.cmap\n",i,pmap->id,filename);
	  fflush(stdout);
	}

	output_cmap(filename, Gmap, i, i);
      }
#endif

    } else {/* prepare for pairalign_allpairs2() to be called again */
      PairMergeProgress = nummaps - orignummaps;

      /* NOTE : score_free must be called here before the arrays are rearranged */
      score_free(XXmap,numX);XXmap = 0;numX= 0;
      score_free(YYmap,numY);YYmap = 0;numY= 0;

#if 1 // new code : Always append merged maps at end of Gmap[] array and set Cfirst so Gmap[Cfirst..nummaps-1] are the merged maps, then set CfirstPairs=2
      int j = 0;
      //      int newCfirst = -1;
      for(int i = 0; i < orignummaps; i++){
	Cmap *pmap = Gmap[i];
	if(pmap->paired)
	  continue;
	if(j < i){
	  Cmap *tmp = Gmap[j];
	  Gmap[j] = Gmap[i];
	  Gmap[i] = tmp;
	}
	j++;
      }
      Cfirst = j;
      for(int i = orignummaps; i < nummaps; i++){
	Gmap[i]->paired = 1;// marks map as newly merged : in next round all pairs must include at least one such map
	if(j < i){
	  Cmap *tmp = Gmap[j];
	  Gmap[j] = Gmap[i];
	  Gmap[i] = tmp;
	}
	j++;
      }
      nummaps = j;
      CfirstPairs = 2;
      
#else // oldcode

      /* remove maps with paired = 1 from Gmap[0..orignummaps-1] */
      if(!Cfirst || CfirstPairs == 2){
	/* If Cfirst && CfirstPairs == 2 : include newly merged maps in second map set Gmap[Cfirst .. nummaps-1] */
	int j = 0;
	int newCfirst = -1;
	for(int i = 0; i < orignummaps; i++){
	  if(i == Cfirst)
	    newCfirst = j;
	  Cmap *pmap = Gmap[i];
	  if(pmap->paired)
	    continue;
	  if(j < i){
	    Cmap *tmp = Gmap[j];
	    Gmap[j] = Gmap[i];
	    Gmap[i] = tmp;
	  }
	  j++;
	}
	if(DEBUG) assert(newCfirst >= 0 && newCfirst <= Cfirst);
	for(int i = orignummaps; i < nummaps; i++){
	  Gmap[i]->paired = 1;// marks map as newly merged : in next round all pairs must include at least one such map
	  if(j < i){
	    Cmap *tmp = Gmap[j];
	    Gmap[j] = Gmap[i];
	    Gmap[i] = tmp;
	  }
	  j++;
	}
	nummaps = j;
	Cfirst = newCfirst;
      } else if(CfirstPairs == 1){/* must include newly merged maps in first map set Gmap[0..Cfirst-1] */
	/* save 2nd map set Gmap[Cfirst .. orignummaps-1] in mapset2[0..orignummaps-1-Cfirst] */
	Cmap **mapset2 = new Cmap*[orignummaps - Cfirst];
	int j = 0;
	for(int i = Cfirst; i < orignummaps; i++){
	  mapset2[j++] = Gmap[i];
	  Gmap[i] = NULL;
	}
	if(VERB>=2){
	  printf("2nd map set is Gmap[%d .. %d] : moved to mapset2[0..%d]\n",
		 Cfirst,orignummaps-1,j-1);
	  fflush(stdout);
	}
	
	j = 0;
	/* remove maps with paired==1 from 1st mapset (map[0..Cfirst-1]) */
	for(int i = 0; i < Cfirst; i++){
	  Cmap *pmap = Gmap[i];
	  if(pmap->paired)
	    continue;
	  if(j < i){
	    Cmap *tmp = Gmap[j];
	    Gmap[j] = Gmap[i];
	    Gmap[i] = tmp;
	  }
	  j++;
	}
	if(VERB>=2){
	  printf("1st map set Gmap[0..%d] : reduced to Gmap[0..%d] after removed merged maps\n",Cfirst-1,j-1);
	  fflush(stdout);
	}

	/* append newly merged maps to 1st mapset (now in Gmap[0..j-1]) */
	for(int i = orignummaps; i < nummaps; i++){
	  Gmap[i]->paired = 1;// marks map as newly merged : in next round all pairs must include at least one such map
	  if(j < i){
	    Cmap *tmp = Gmap[j];
	    Gmap[j] = Gmap[i];
	    Gmap[i] = tmp;
	  }
	  j++;
	}
	int newCfirst = j;
	
	if(VERB>=2){
	  printf("After appending merged maps Gmap[%d..%d], extended 1st map set is Gmap[0..%d] : Cfirst = %d -> %d\n",
		 orignummaps,nummaps-1,j-1,Cfirst,newCfirst);
	  fflush(stdout);
	}

	/* move all null values in Gmap[j .. nummaps-1] to end of this range */
	int k = j;
	for(int i = k; i < nummaps; i++){
	  if(Gmap[i] == NULL)
	    continue;
	  if(k < i){
	    if(DEBUG) assert(Gmap[k] == NULL);
	    Gmap[k] = Gmap[i];
	    Gmap[i] = NULL;
	  }
	  k++;
	}
	/* Gmap[k .. nummaps-1] will now be all nulls : move the 2nd map set pointers into this space */
	if(DEBUG) assert(nummaps - k == orignummaps - Cfirst);
	for(int i = 0; i < orignummaps - Cfirst; i++)
	  Gmap[k+i] = mapset2[i];
	
	/* filter 2nd mapset for maps with paired==1 and append it to new 1st mapset Gmap[0..newCfirst-1] == Gmap[0..j-1] */
	for(int i = k; i < nummaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(pmap->paired)
	    continue;
	  if(j < i){
	    Cmap *tmp = Gmap[j];
	    Gmap[j] = Gmap[i];
	    Gmap[i] = tmp;
	  }
	  j++;
	}
	if(VERB>=2){
	  printf("After filtering and appending 2nd map set Gmap[%d..%d], combined sets are Gmap[0..%d]\n",
		 k,nummaps-1,j-1);
	  fflush(stdout);
	}

	nummaps = j;
	Cfirst = newCfirst;
      } else if(CfirstPairs == 0){/* place newly merged maps in first map set and all other maps in 2nd map set and set CfirstPairs to 1 */
	/* save both map set Gmap[0 .. orignummaps-1] in mapset2[0..orignummaps-1] */
	Cmap **mapset2 = new Cmap*[orignummaps];
	int j = 0;
	for(int i = 0; i < orignummaps; i++){
	  mapset2[j++] = Gmap[i];
	  Gmap[i] = NULL;
	}
	if(VERB>=2){
	  printf("original map sets Gmap[0 .. %d] : moved to mapset2[0..%d]\n",
		 orignummaps-1,j-1);
	  fflush(stdout);
	}
	
	j = 0;

	/* append newly merged maps to Gmap[0..j-1] */
	for(int i = orignummaps; i < nummaps; i++){
	  Gmap[i]->paired = 1;// marks map as newly merged : in next round all pairs must include at least one such map
	  if(j < i){
	    Cmap *tmp = Gmap[j];
	    Gmap[j] = Gmap[i];
	    Gmap[i] = tmp;
	  }
	  j++;
	}
	int newCfirst = j;
	
	if(VERB>=2){
	  printf("After appending merged maps Gmap[%d..%d], 1st map set is Gmap[0..%d] : Cfirst = %d -> %d\n",
		 orignummaps,orignummaps-1,j-1,Cfirst,newCfirst);
	  fflush(stdout);
	}

	/* move all null values in Gmap[j .. nummaps-1] to end of this range */
	int k = j;
	for(int i = k; i < nummaps; i++){
	  if(Gmap[i] == NULL)
	    continue;
	  if(k < i){
	    if(DEBUG) assert(Gmap[k] == NULL);
	    Gmap[k] = Gmap[i];
	    Gmap[i] = NULL;
	  }
	  k++;
	}
	/* Gmap[k .. nummaps-1] will now be all nulls : move the two map sets mapset2[0..orignummaps-1] into this space */
	if(DEBUG) assert(nummaps - k == orignummaps);
	for(int i = 0; i < orignummaps; i++){
	  if(DEBUG) assert(Gmap[k+i] == NULL);
	  Gmap[k+i] = mapset2[i];
	}
	
	/* filter both original mapets at Gmap[k..nummaps-1] and append to new 1st mapset Gmap[0..newCfirst-1] == Gmap[0..j-1] */
	for(int i = k; i < nummaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(pmap->paired)
	    continue;
	  if(j < i){
	    Cmap *tmp = Gmap[j];
	    Gmap[j] = Gmap[i];
	    Gmap[i] = tmp;
	  }
	  j++;
	}
	if(VERB>=2){
	  printf("After filtering and appending original two map sets Gmap[%d..%d], combined sets are Gmap[0..%d]\n",
		 k,nummaps-1,j-1);
	  fflush(stdout);
	}

	nummaps = j;
	Cfirst = newCfirst;
	CfirstPairs = 1;
      }
#endif // old code
    } // if(PairMergeRepeat && nummaps > orignummaps 
  }

  if(PairsplitSeparateQueries)/* reset copied array pointers to avoid call free() on them */
    for(int Yid = orignumY; Yid < numY; Yid++)
      YYmap[Yid]->init();

  if(LEAKDEBUG || PairMergeProgress){
    if(!PairMergeProgress){
      score_free(XXmap,numX);XXmap = 0;numX= 0;
      score_free(YYmap,numY);YYmap = 0;numY= 0;
    }

    for(int c = 0; c < colors; c++){
      delete [] misaligns[c];
      delete [] errors[c];
    }

    if(hashpairs1){
      delete [] hashpairs1;
      hashpairs1 = 0;
    }
    if(hashpairs2){
      delete [] hashpairs2;
      hashpairs2 = 0;
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

    if(A_threads) {      /* free memory of *thread_arrays* */
      for(int i= 0; i < A_numthreads; i++) {
	for(int c = 0; c < colors; c++){
	  delete [] A_threads[i][c].score;
	  delete [] A_threads[i][c].Uscore;
	  delete [] A_threads[i][c].Rscore;
	  delete [] A_threads[i][c].G;
	  delete [] A_threads[i][c].H;
	  delete [] A_threads[i][c].Rijy;
	  delete [] A_threads[i][c].Lijy;
	  delete [] A_threads[i][c].Rijx;
	  delete [] A_threads[i][c].Lijx;
	}
      }
      delete [] ALmem; ALmem = NULL;
      delete [] A_threads; A_threads = NULL;
    }

    if(Fmem_threads){
      delete [] Fmem_threads[0];
      delete [] Fmem_threads; Fmem_threads = NULL;
    }
    if(Imem_threads) {
      delete [] Imem_threads[0];
      delete [] Imem_threads; Imem_threads = NULL;
    }
  
    A_numthreads = mem_numthreads = 0;
    for(int c = 0; c < colors; c++)
      A_maxN[c] = mem_maxN[c] = mem_maxM[c] = 0;
  }
  return PairMergeProgress;
}
	      
