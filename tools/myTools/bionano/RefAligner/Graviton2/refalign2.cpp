#include <sys/types.h>
#include <sys/stat.h>
#ifndef WIN32
#include <unistd.h>
#else
#define copysign _copysign
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <omp.h>

#include "constants.h"
#include "parameters.h"
#include "globals.h"
#include "Calign.h"
#include "Ccontig.h"
#include "hash.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/trunk/refalign.cpp 1141 2013-01-28 18:18:43Z tanantharaman $");

// #undef DEBUG
// #define DEBUG 2

#define FENCE (DEBUG>=2 ? 1024 : 0) /* size of fence before and after each row of data : currently only used with DIAGNONAL >=2 */

#define OUTLIER_MARGIN (USE_RFLOAT ? 0.001 : 1e-8) /* margin used to classify an interval as outlier during parameter estimation */
#define SCORE_MARGIN (USE_RFLOAT ? 1e-4 : 1e-10) /* relative score accuracy margin (to trigger inconsistency assertions) */

#define RESDATA ResEstimate /* optimize res and resSD values */

#define STITCH 0 /* support stitching signal (still being developed) */

#define MIS_VITERBI 1 /* WAS 0 */ /* use Viterbi scoring of misaligned sites */

#define LFLOAT FLOAT /* WAS RFLOAT */ /* floating point precision used by logLR() and resLL() for parameter estimation */
#define LZERO ZERO /* WAS RZERO */

#define FN_GM (SCORE_APPROX <= 1 && !MIS_VITERBI) /* Use Golden Mean search to optimize FN (instead of just counting misaligned sites) : slower but more accurate for high density data with many misresolved sites */
#define MINMEM 1 /* WAS 0 */ /* minimize memory usage by reallocating lightweight heap each iteration of giter (and if >=2 after iteration iteration of refid) */

#define FIX_OUTLIER_MIS 1 /* fix scoring of outliers to not including misalignment scores Sm() (from between alignment intervals) */

#define DIAGONAL 2 /* > 0 : with hashtable offset, diagonalize alignment array to increase speed and (if >= 2) reduce real memory usage */
#define MIN_MEM (!USE_MIC) /* minimize memory usage by reallocating memory in each call of refalignXYsd() : only used with DIAGONAL>=2 */

#undef TRACE
#define TRACE 0 /* (giter==RefRepeats-1) */ /* Trace alignment scoring for specified rmap->id, nanomap->id, orientation, I,K,J */

#define REF_TRACE 1LL
#define MAP_TRACE 28316LL
#define OR_TRACE 0

#define C_TRACE -1
#define I_TRACE -1
#define K_TRACE -1
#define J_TRACE -1

#define PVERB ((TRACE && YYmap[refid]->id==REF_TRACE && XXmap[mapid]->id==MAP_TRACE && orientation==OR_TRACE) ? 1 : 0) /* Display detailed alignment information for all alignments above threshold */

#define BESTREF_PV BestRefPV /* with -BestRef chose best reference based on LogPV instead of score */

#define RESBIAS_MINSAMPLE 128 /* minimum number of samples per bin */
#define RESBIAS_MAXSLOPE 0.99 /* maximum linear slope allowed */

#define FIX_CONSTRAINT 1 /* TRY 1,2 */ /* > 0 : more conservative RepeatRec constraint */

#define HASH_DEBUG 0 /* To debug HASH_STREAM with verbose buffer update messages (2 = very verbose, each match is displayed) */

#define FLAT_MEMORY 1 /* Use flat memory for AL.field(I,K,J) : faster (if vectorized) but uses 4x more virtual memory */

#define DEFER_BP 0 /* defer computation of DP backpointer to allow recurrance to be faster (no branch penalty, just a max operator) : HERE HERE Not yet implemented */

#define SCORE_ALLREF 1 /* average scoring function over all reference maps */

#define QQPLOT 0 /* output a QQ-plot of error norms to stdout (first iteration only) */
#define PROFILE 0 /* If > 0, print coverage profile for above threshold hits to stdout (exclude last (PROFILE-1) sites for each alignment to better reveal gaps) */

#define REFSCORE 1 /* new scoring method with Bias a function of X length and sites M only */

#define FN_VITERBI 0 /* try to compensate false negative estimate for Viterbi greediness (if resSD > 0 && K > 0) Does not seem to help */

#define NEW 2  /* 1 : handle end outliers using local alignment if PoutlierEnd > 0.0 (works just like -local PoutlierEnd 1 without Bias for outlier ends) \
		  2 : Add Bias for outlier ends if PoutlierEnd > 0.0 (makes alignment independent of Bias factor) */

#define RESDATA_FIX 1 /* Use updated x interval estimates to reflect changes due to -resbias and scaling (bpp) */

#define SCAN_FRATE_FIX 1 /* include -Frate*Xlengths in LL when doing per-scan scaling */

#define WITH_RESBIAS 1 /* Specify REFDEBUG debugging test case (0 = assertion at line 8823, 1 = asssertion failure with -resbias) */

#define SDDEBUG 0 /* more detailed output of SD,SF based logLR() updates */

#define REFDEBUG ((DEBUG>=3 && REFSCORE && SIMPLE_BIAS && FAST_OUTLIER && COLOR_SHIFT==2) ? (SDDEBUG ? 2 : 1) : 0)

#define REFDEBUG_STRICT 0 /* Guarantee that logLR() always improves by :
			     >=1 : Estimate Kmax[] only once before first -M iteration
				   Disable lookahead in -resEstimate (see resLL())
                             >=2 : Restrict C scaling factor to values <= 1.0 (see NOTE at Cscale declaration in refalign_allpairs2)
			           Disable bounding of error parameters
			           Disable -resbias
			  */

#include "RGentigRefScore2.h"

static double origScoreThreshold = -10000.0;/* ScoreThreshold before applying -MapRate */
static double origLogPvThreshold = 10.0;/* LogPvThreshold before applying -MapRate */

#define ENDFIX 3 /* WAS 2 */  /* 1 : allow for sizing error at ends
				 2 : Also use ends with Sbnd() in sizing error estimation 
				 3 : Also include ends with Send() in sizing error estimation, when y < C*x, for greater agreement with logLR() */

#define ENDFIX2 1 /* WAS 0 */ /* more thorough, bug slightly slower, checking for sizing error at ends of alignments */

#define FAST Mfast /* 1 : limit alignment to best refid,orientation and scaling from previous iteration, if score was > FAST_SCORE \
		      2 : Also limit alignment to best alignment location from previous iteration +- DELTA_Y */
#define FAST_SCORE 20.0  /* WAS 10.0 for 1 color */

#define FP_DIST err_dist /* seperately count fp and fn sites that are at least this distance from the nearest aligned site : NOT YET IMPLEMENTED */

#if USE_MIC
#define memcpy _mm512_memcpy
#endif

extern Cmap **YYmap,**XXmap;

extern double zscore(double pvalue);

/** consensus site linked list */
class Csite {
public:
  Csite *next;/**< next site in contig */
  FLOAT loc;/**< location of site relative to left end of original reference (can extend beyond ends of reference map) */
  int I;/**< index into reference map (or -1 if site is not present on reference map, -2 if this is a map end) */
  int site;/**< index into contig */
};

/** aligned interval statistics */
class Cinterval {
public:
  double x,y;/**< query and Reference interval sizes that aligned with each other */
  double len;/**< either same as y (c < 2) or length (c==2) */
  double resvar;/**< resolution error (errL*errL + errR*errR) : Only used if RES_VARIANCE >= 1*/
  double err;  /**< NOTE normalized error = (x-y - MA[c]*end)/sqrt(SF[c]*SF[c]+len*SD[c]*SD[c]+SR[c]*SR[c]*y*y+SE[c]*SE[c]*resvar) */
  double rawx;/* Query map interval size based on rawsite[c][0] instead of set[0] and its log() value */
  int mapid;/**< map that aligned with reference is Gmap[mapid] */
  int L,R;/**< left and right end of aligned query interval Gmap[mapid]->site[c][L..R] (if c < color)
	     right ends of query aligns for colors 0 and 1 Gmap[mapid]->site[0][L] or Gmap[mapid]->site[1][R] (c==colors) 
	     NOTE : x == end*(Gmap[mapid]->site[1][R] - Gmap[mapid]->site[0][L]) */
  short end;/**< For color < colors :  1 IFF this is an end interval, with sizing error term due to Sbnd(), hence no log(SF[c]*SF[c]+SD[c]*SD[c]*y+SR[c]*SR[c]*y*y+SE[c]*SE[c]*resvar) term
		For color == colors : -1 or 1, depending on orientatation of Molecule */
  short stitch;/**< 0 : regular interval without stitch point
		!= 0 : this interval has an image stitch point, -1 or 1 depending on orientation of Query map */
#if REFDEBUG >= 1
  double origx;
  double origlen;
#endif
};

static double aNaN = nan("NaN");

class Cresdata {
public:
  Calign *align;/* Y == refmap[align->mapid1]->site[c] */
  int U;/* score segment in alignment (0 == left end, 1 == first internal alignment interval etc) */

  int I,K;/* interval that did NOT resolve is Y[I-K .. I] */
  int n,T;/* If n != 0 : interval that resolved is Y[I-K-n .. I-K] and its size is Yc(Y,I,K) - Yc(Y,I-K-n,T) */
  int J,m;/* If n != 0 : interval X[J-m..J] that aligned with Yc(Y,I,K) - Yc(Y,I-K-n,T) (NOTE: J is always valid to support REFDEBUG) */
  int outlier;/* If n != 0 : Force interval to be scored as outlier, or NOT (required with REFDEBUG_STRICT) */
  double x;/* If n != 0 size of interval X[J-m..J] */

  /* Following fields only valid if Ip != 0 and describe an alternate alignment with unresolved interval that is one site larger than Y[I-K..I] */
  int Ip;/* Alternate larger interval that did NOT resolve is Y[Ip-K-1 .. Ip] */
  int I1,K1,n1,T1,J1,m1;/* If n1 != 0 : original interval that resolved if Ip,Kp is NOT used (same as I,K,n,T,J,m if I==Ip) */
  int I2,K2,n2,T2;/* If n1 != 0 : alternate interval that resolved if IP,Kp is used with same J1,m1 */
  double x1;/* If n1 != 0 : corresponding size of interval X[J1-m1..J1] */

  /* NOTE 
     Let Sint  == (n ? Sint(x, Yc(Y,I,K) - Yc(Y,I-K-n,T), m,n,J,I,K,T,Y,c) : Sm(0,I,K,Y,c))
         Sintp == (n ? Sint(x, Yc(Y,Ip,K+1) - Yc(Y,I-K-n,T), m, n-1 + (Ip-I),J,Ip,K+1,T,Y,c) : Sm(0,Ip,K+1,Y,c))
         Sint1 == Sint(x1, Yc(Y,I1,K1) - Yc(Y,I1-K1-n1,T1), m1, n1, J1, I1, K1,T1,Y,c)
         Sint2 == Sint(x1, Yc(Y,I2,K2) - Yc(Y,I2-K2-n2,T2), m1, n2, J1, I2, K2,T2,Y,c)

     score is larger of the following two expressions of res,resKB (If Ip == 0, only the first expression is used) :
     
     Sint

     n1 ? Sintp + Sint2 - Sint1 : Sm(0,Ip,K+1,Y) + FnPenalty

  */
  Cresdata()
  {
    if(DEBUG/* HERE >=2 */){
      align = NULL;
      U = I = K = n = T = J = m = -999;
      outlier = -1;
      Ip = -999;
      I1 = K1 = n1 = T1 = J1 = m1 = -999;
      I2 = K2 = n2 = T2 = -999;
      
      x = x1 = aNaN;
    }
  }
};

//static Ccontig *gcontig = 0;

static lightweight_heap *rs_heap = NULL;

static CHashMatch *hashpairs1 = 0, *hashpairs2 = 0;
static CHashMatch *nexthash1 = 0, *nexthash2 = 0;/* nexthash1 points into hashpairs1[0..numhashpairs1-1] to the next unused entry */
static size_t maxhashpairs1 = 0, maxhashpairs2 = 0;
static size_t numhashpairs1 = 0, numhashpairs2 = 0;

static int refid1 = -1, refid2 = -1;/* current and next scheduled first index into refmap[] (last index used is returned by loadhash2()) */
static CHashMatch **hashmatch1 = 0, **hashmatch2 = 0;/* hashmatch1[mapid=0..nummaps-1] points to the first location in hashpairs1[] for HashTable matches for (refid1,mapid) */

static int *YidListMem = NULL;
static int *YidList1 = 0, *YidList2 = 0;
static int *XidList1 = 0, *XidList2 = 0;
static CHashMatch **phashListMem = NULL;
static CHashMatch **phashList1 = 0, **phashList2 = 0;
static int NumPair1 = 0, NumPair2 = 0, MaxPair = 0;/* HERE : try changing to size_t to reduce warnings */
static int refidend1 = -1, refidend2 = -1;/* current and next scheduled index range into refmap[] are refid1..refidend1 and refid2..refidend2 respectively */

static int BinaryFormat = 0;    
static FILE *HashFP = NULL;
static size_t MatchMax = 0;
static size_t totalhashpairs = 0;

static void maxhashalloc(size_t num, size_t numhashpairs, size_t &maxhashpairs, CHashMatch * &hashpairs)
{
  size_t newmax = num;
  if(newmax <= maxhashpairs)
    return;
  if(maxhashpairs > 0 && newmax < maxhashpairs * 2)
    newmax = maxhashpairs * 2;

  if(VERB && HASH_DEBUG){
    printf("maxhashalloc:&hashpairs=%p,hashpairs=%p,maxhashpairs=%llu->%llu, size=%llu bytes\n",&hashpairs, hashpairs,(unsigned long long)maxhashpairs,(unsigned long long)newmax,(unsigned long long)newmax*sizeof(CHashMatch));
    fflush(stdout);
  }

  if(maxhashpairs <= 0)
    hashpairs = new CHashMatch[newmax+1];
  else {
    CHashMatch *orighashpairs = hashpairs;
    hashpairs = new CHashMatch[newmax+1];
    for(size_t i = 0; i < numhashpairs; i++)
      hashpairs[i] = orighashpairs[i];
    hashpairs[numhashpairs].id1 = 0;
    delete [] orighashpairs;
  }
  maxhashpairs = newmax;
  if(VERB && HASH_DEBUG){
    printf("maxhashalloc:hashpairs=%p,maxhashpairs=%llu: wall time=%0.6f secs\n",hashpairs,(unsigned long long)maxhashpairs, wtime());
    fflush(stdout);
  }
}

#if 0
static int SVQrystartposInc(xmapEntry **pp1, xmapEntry **pp2)
{
  double start1 = pp1[0]->qrystartpos;
  double start2 = pp2[0]->qrystartpos;
  double end1 = pp1[0]->qryendpos;
  double end2 = pp2[0]->qryendpos;
  return (start1 > start2) ? 1 : (start1 < start2) ? -1 : 
    (end1 > end2) ? 1 : (end1 < end2) ? -1 : 0;
}
#endif

/* qsort() intcmp function to sort CHashMatch array in increasing order of id2,orientation,(-hashscore),offset */
static inline int ChashIdInc(CHashMatch *p1, CHashMatch *p2)
{
  return CHashMatchIncId1(p1,p2);
}

/** load data into specified hash buffer (hashpairs[0..numhashpairs-1], plus a invalid entry at hashpairs[numhashpairs]) 
    Then copy hashtable entries to YidList[],XidList[],pHashList[] until we have MaxPair Entries (or we reach end of hashtable).
    First expected refid value (on YidList[]) is specified.
    Make sure not to start new refid value (on YidList[]) unless all mapid values are available and will fit in XidList[] : Return last refid value processed (or -1 if none) */
static int loadhash2(int refid, int nummaps, CHashMatch* &nexthash, CHashMatch *hashpairs, size_t &numhashpairs, size_t maxhashpairs, int *YidList, int *XidList, CHashMatch **phashList, int &NumPair)
{
  if(DEBUG) assert(BinaryFormat && HashFP);

  /* check if the 2nd buffer already has the data we need */
  if(hashpairs == hashpairs1 && refid == refid2 && refidend2 >= refid2 ){/* just swap the two buffers */
    if(VERB && HASH_DEBUG){
      #pragma omp critical
      {
	if(nexthash2 < &hashpairs2[numhashpairs2])
	  printf("Swapping Buffers: Matches1=%llu/%llu,NumPair1=%d,refid1=%d..%d, Matches2=%llu/%llu(id1=%lld..%lld,id2=%lld..%lld),NumPair2=%d,refid2=%d..%d, refid=%d: CPU time=%0.6f, wall time=%0.6f\n",
		 (unsigned long long)(&hashpairs1[numhashpairs1]-nexthash1), (unsigned long long)numhashpairs1, NumPair1,refid1, refidend1, 
		 (unsigned long long)(&hashpairs2[numhashpairs2]-nexthash2), (unsigned long long)numhashpairs2, 
		 (long long)nexthash2->id1,(long long)hashpairs2[numhashpairs2-1].id1,
		 (long long)nexthash2->id2,(long long)hashpairs2[numhashpairs2-1].id2, 
		 NumPair2, refid2, refidend2, refid, mtime(), wtime());
	else
	  printf("Swapping Buffers: Matches1=%llu/%llu,NumPair1=%d,refid1=%d..%d, Matches2=%llu/%llu,NumPair2=%d,refid2=%d..%d, refid=%d: CPU time=%0.6f, wall time=%0.6f\n",
		 (unsigned long long)(&hashpairs1[numhashpairs1]-nexthash1), (unsigned long long)numhashpairs1, NumPair1,refid1, refidend1, 
		 (unsigned long long)(&hashpairs2[numhashpairs2]-nexthash2), (unsigned long long)numhashpairs2, NumPair2, refid2, refidend2, refid, mtime(), wtime());
	fflush(stdout);
      }
    }

    /* there should be no data left in 1st buffer */
    if(DEBUG) assert(nexthash1 == &hashpairs1[numhashpairs1]);
    if(DEBUG) assert(YidList == YidList1 && XidList == XidList1 && phashList == phashList1);
    if(DEBUG) assert(NumPair1 == 0);
    
    /* swap the two sets of buffers */
    CHashMatch *hashtmp = hashpairs1; hashpairs1 = hashpairs2; hashpairs2 = hashtmp;
    hashtmp = nexthash1; nexthash1 = nexthash2; nexthash2 = hashtmp;
    if(DEBUG) assert(maxhashpairs1 == maxhashpairs2);
    size_t tmpSize = numhashpairs1; numhashpairs1 = numhashpairs2; numhashpairs2 = tmpSize;

    int *tmpList = YidList1; YidList1 = YidList2; YidList2 = tmpList;
    tmpList = XidList1; XidList1 = XidList2; XidList2 = tmpList;
    CHashMatch **phashtmp = phashList1; phashList1 = phashList2; phashList2 = phashtmp;
    int tmpPair = NumPair1;NumPair1 = NumPair2; NumPair2 = tmpPair;

    refid1 = refid2; refid2 = -1;
    refidend1 = refidend2; refidend2 = -1;

    /* NOTE : there will sometimes be surplus data in hashpairs1[] starting at nexthash1 that is not in YidList1[],XidList1[],phashList1[] */

    return refidend1;
  }
  
  if(DEBUG) assert(nexthash2 == &hashpairs2[numhashpairs2]);/* there should be no data left in 2nd buffer */

  if(hashpairs==hashpairs1){
    CHashMatch *end1 = &hashpairs1[numhashpairs1];
    if(nexthash1 < end1 && nexthash1 > hashpairs1){/* copy surplus data in 1st buffer to start of buffer : should happen rarely */
      if(VERB && HASH_DEBUG){
	#pragma omp critical
	{
	  printf("copying %llu entries from hashpairs1[%llu..%llu] to hashpairs1[0..%llu] (id1=%lld..%lld,id2=%lld..%lld): cpu time=%0.6f, wall time=%0.6f secs\n",
		 (unsigned long long)(end1-nexthash1), (unsigned long long)(nexthash1-hashpairs1), (unsigned long long)(end1-hashpairs1)-1, (unsigned long long)(end1-nexthash1)-1, 
		 (long long)nexthash1->id1, (long long)end1[-1].id1, (long long)nexthash1->id2, (long long)end1[-1].id2, mtime(), wtime());
	  fflush(stdout);
	}
      }

      int len = end1-nexthash1;

#ifndef VALGRIND
      if(USE_MIC ? (nexthash1 - hashpairs1)*sizeof(CHashMatch) >= 512 : (nexthash1-hashpairs1) >= len)
	memcpy(hashpairs1, nexthash1, len * sizeof(CHashMatch));
      else
#endif
	memmove(hashpairs1, nexthash1, len * sizeof(CHashMatch));

      numhashpairs1 = len;
      nexthash1 = hashpairs1;

      hashpairs1[numhashpairs1].id1 = -1;    /* initialize invalid entry at hashpairs[numhashpairs] */
    }
  }

  if(hashpairs==hashpairs2){
    CHashMatch *end1 = &hashpairs1[numhashpairs1];
    if(nexthash1 < end1){/* copy surplus data in 1st buffer to start of 2nd buffer : This happens only during multi-threaded section */
      if(VERB && HASH_DEBUG){
	#pragma omp critical
	{
	  printf("copying %llu entries from hashpairs1[%llu..%llu] to hashpairs2[0..%llu] (id1=%lld..%lld,id2=%lld..%lld): cpu time=%0.6f, wall time=%0.6f secs\n",
		 (unsigned long long)(end1-nexthash1), (unsigned long long)(nexthash1-hashpairs1), (unsigned long long)(end1-hashpairs1)-1, (unsigned long long)(end1-nexthash1)-1, 
		 (long long)nexthash1->id1, (long long)end1[-1].id1, (long long)nexthash1->id2, (long long)end1[-1].id2, mtime(), wtime());
	  fflush(stdout);
	}
      }
      size_t len = end1-nexthash1;
      memcpy(hashpairs2, nexthash1, len * sizeof(CHashMatch));
      numhashpairs2 = len;
      nexthash2 = hashpairs2;
      numhashpairs1 = 0;
      nexthash1 = hashpairs1;

      hashpairs2[numhashpairs2].id1 = -1;    /* initialize invalid entry at hashpairs[numhashpairs] */
    }
  }
  
  /* update refid1 or refid2 */
  if(hashpairs==hashpairs1){
    refid1 = refid;
    refidend1 = -1;
  } else {
    if(DEBUG) assert(hashpairs==hashpairs2);
    refid2 = refid;
    refidend2 = -1;
  }

  /* first fill up current buffer with data from HashFP */
  int HashEof = hash_eof(HashFP);
  if(HashEof && numhashpairs <= 0)
    return -1;

  size_t cur = numhashpairs;/* entries left */
  size_t num = maxhashpairs - numhashpairs;/* space left */

  if(!HashEof){
    int linecnt = 0;/* Not Used for Binary File input */
    if(VERB && HASH_DEBUG){
      #pragma omp critical
      {
	printf("loadhash2(refid=%d):numhashpairs=%llu,maxhashpairs=%llu,num=%llu:Calling hash_read(Hash,0,&hashpairs[numhashpairs],num,linecnt):CPU time=%0.6f, wall time=%0.6f secs\n",
	       refid,(unsigned long long)numhashpairs,(unsigned long long)maxhashpairs,(unsigned long long)(maxhashpairs - numhashpairs),mtime(),wtime());
	fflush(stdout);
      }
    }
    numhashpairs += hash_read(HashFP, 0, &hashpairs[numhashpairs], num, linecnt);
    HashEof = hash_eof(HashFP);
    if(VERB && HASH_DEBUG && numhashpairs > cur){
      #pragma omp critical
      {
	printf("loadhash2(refid=%d):Read in next %llu Matches (total=%llu/%llu) (EOF=%d) into hashpairs%d[]: cpu time=%0.6f, wall time=%0.6f secs (totalhashpairs=%llu,cur=%llu,nexthash=%llu,numhashpairs=%llu,id1=%lld..%lld,id2=%lld..%lld)\n", 
	       refid,(unsigned long long)(numhashpairs-cur), (unsigned long long)(totalhashpairs + numhashpairs - cur), (unsigned long long)MatchMax, HashEof, hashpairs==hashpairs1 ? 1 : 2, mtime(), wtime(), 
	       (unsigned long long)totalhashpairs,(unsigned long long)cur,(unsigned long long)(nexthash-hashpairs),(unsigned long long)numhashpairs,
	       (long long)nexthash->id1, (long long)hashpairs[numhashpairs-1].id1, (long long)nexthash->id2, (long long)hashpairs[numhashpairs-1].id2);
	for(size_t i = cur; i < numhashpairs; i++)
	  if(XXmap[hashpairs[i].id2]->id == MAP_TRACE){
	    printf("hashpairs[%llu]:id1=%d(%lld),id2=%d(%lld),orientation=%d,offset=%d,hashscore=%d\n",
		   (unsigned long long)i,hashpairs[i].id1,YYmap[hashpairs[i].id1]->id,hashpairs[i].id2,XXmap[hashpairs[i].id2]->id,
		   hashpairs[i].orientation,hashpairs[i].offset,hashpairs[i].hashscore);
	    fflush(stdout);
	  }
	fflush(stdout);
      }
    }
    if(DEBUG>=2){    /* check if hash pairs are in ascending order of id1,id2,orientation,offset */
      size_t i;
      for(i = 1; i < numhashpairs; i++)
	if(ChashIdInc(&hashpairs[i-1],&hashpairs[i]) > 0)
	  break;
      if(i < numhashpairs){
	printf("Hashtable must be sorted in ascending order of id1,id2,orientation,-hashscore,offset\n");
	printf("totalhashpairs=%llu,numhashpairs=%llu,cur=%llu\n",(unsigned long long)totalhashpairs,(unsigned long long)numhashpairs,(unsigned long long)cur);
	for(size_t j = i-1; j <= i; j++)
	  printf("hashpairs[%llu]:id1=%lld,id2=%lld,orientation=%d,offset=%d\n",(unsigned long long)j,(long long)hashpairs[j].id1,(long long)hashpairs[j].id2,hashpairs[j].orientation,hashpairs[j].offset);
	fflush(stdout);
	assert(0);
	exit(1);
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
    if(VERB && HASH_DEBUG && j < numhashpairs){
      #pragma omp critical
      {
	printf("Reduced matches from %llu to %llu due the hashscore < %d\n",(unsigned long long)numhashpairs, (unsigned long long)j, hash_threshold);
	fflush(stdout);
      }
    }
    numhashpairs = j;

    hashpairs[numhashpairs].id1 = -1;    /* initialize invalid entry at hashpairs[numhashpairs] */
  }

  CHashMatch *orignext = nexthash;
  if(VERB && HASH_DEBUG){
    #pragma omp critical
    {
      printf("refid=%d:nexthash= &hashpairs%d[%llu]\n",refid,hashpairs==hashpairs1 ? 1 : 2, (unsigned long long)(nexthash-hashpairs));
      fflush(stdout);
    }
  }

  /* use data in current buffer to create YidList[], XidList[], phashList[] */
  if(DEBUG && hashpairs == hashpairs1) assert(NumPair1 == 0);
  if(DEBUG && hashpairs == hashpairs2) assert(NumPair2 == 0);

  CHashMatch *end = &hashpairs[numhashpairs];
  if(DEBUG && nexthash < end && !(nexthash->id1 >= refid)){
    printf("refid=%d:nexthash= &hashpairs%d[%llu],end-nexthash=%llu,nexthash->id1=%d,refid1=%d,refid2=%d\n",
	   refid,hashpairs==hashpairs1 ? 1 : 2, (unsigned long long)(nexthash-hashpairs), (unsigned long long)(end-nexthash),nexthash->id1,refid1,refid2);
    fflush(stdout);
    assert(nexthash->id1 >= refid);
  }
  int origrefid = refid;
  int mapid = -1;
  int id1max = hashpairs[numhashpairs-1].id1;
  if(id1max <= nexthash->id1 && !HashEof){/* need to enlarge buffers to hold complete set of all hashmatches for nexthash->id1 */
    printf("refid=%d:nexthash->id1=%d,numhashpairs=%lu,hashpairs[numhashpairs-1].id1=%d,nummaps=%d,MaxPair=%d\n",
	   refid,nexthash->id1,numhashpairs,hashpairs[numhashpairs-1].id1,nummaps,MaxPair);
    printf("loadhash2(): Increase size of hashtable buffer maxhashpairs=%lu to hold all hashtable matches from any single refid\n", maxhashpairs);
    exit(1);
  }

  //  int id2max = hashpairs[numhashpairs-1].id2;
  for(NumPair= 0; NumPair < MaxPair; ){
    if(nexthash >= end-1){/* next to end of current read buffer : check if this is EOF */
      if(!HashEof){
	if(DEBUG/* HERE >=2 */) assert(nexthash >= end);/* NOTE : this should be the case since the last id1 at &end[-1] == hashpairs[numhashpairs-1] should not be copied unless HashEof */
	break;/* reuse last entry in next Yid iteration since it may come in a pair of entries (both orientations) that are needed together */
      }
      /* EOF : force break out of Yid loop, unless the last entry still needs to be processed */
      if(nexthash >= end)
	break;
    }

    if(DEBUG>=2) assert(origrefid <= nexthash->id1 && nexthash->id1 < numrefmaps && refmap[nexthash->id1]->id >= 0);
    if(DEBUG>=2) assert(0 <= nexthash->id2 && nexthash->id2 < nummaps && gmap[nexthash->id2]->id >= 0);

    /* don't start a new refid (id1) unless we can complete all mapid (id2) values for this new refid  : this is typically handled at end of loop by backtracking, to avoid terminating prematurely */
    if(!HashEof && nexthash->id1 >= id1max) {/* last id will be incomplete unless HashEof==1 */
      if(VERB && (VERB>=2 || HASH_DEBUG)){
	#pragma omp critical
	{
	  printf("refid=%d..%d:nexthash->id1=%d,end-nexthash=%llu,nummaps=%d,MaxPair=%d,NumPair=%d\n",
		 origrefid,refid,nexthash->id1,(unsigned long long)(end-nexthash),nummaps,MaxPair,NumPair);
	  fflush(stdout);
	}
      }
      if(DEBUG && !(nexthash->id1 > refid)){
	#pragma omp critical
	{
	  printf("refid=%d..%d:nexthash->id1=%d,id1max=%d,end-nexthash=%llu,numhashpairs=%llu/%llu,nexthash-hashpairs=%llu,orignexthash-hashpairs=%llu,nummaps=%d,MaxPair=%d,NumPair=%d\n",
		 origrefid,refid,nexthash->id1,id1max,(unsigned long long)(end-nexthash),(unsigned long long)numhashpairs, (unsigned long long)((hashpairs==hashpairs1) ? maxhashpairs1 : maxhashpairs2),
		 (unsigned long long)(nexthash-hashpairs),(unsigned long long)(orignext-hashpairs),nummaps,MaxPair,NumPair);
	  fflush(stdout);
	  assert(nexthash->id1 > refid);
	}
      }
      break;
    }

    YidList[NumPair] = refid = nexthash->id1;
    XidList[NumPair] = mapid = nexthash->id2;

    phashList[NumPair++] = nexthash++;

    /* skip all subsequent entries in hashtable with the same pair of maps as lastYid,lastXid :
       This happens when both orientation or different offsets are in the HashTable */
    while(nexthash < end && nexthash->id1 == refid && nexthash->id2 == mapid)
      nexthash++;
  }

  /* make sure last refid is complete (can be incomplete if NumPair >= MaxPair) */
  if(DEBUG && !HashEof) assert(nexthash < end);
			       
  if((!HashEof || nexthash < end) && nexthash > orignext && nexthash->id1 == nexthash[-1].id1){/* need to undo last incomplete refid */
    if(DEBUG) assert(NumPair >= MaxPair);
    if(DEBUG) assert(refid > origrefid);
    int id1 = nexthash->id1;
    if(DEBUG) assert(id1 == refid);
    while(nexthash[-1].id1 >= id1)
      nexthash--;
    while(NumPair > 0 && YidList[NumPair-1] >= refid)
      NumPair--;
    refid--;
  }
  if(VERB && HASH_DEBUG){
    #pragma omp critical
    {
      printf("Used %llu/%llu matches from hashpairs%d[%llu..%llu] to create YidList[],XidList[] with nummaps=%d: NumPair=%d/%d alignments for refid=%d..%d, nexthash=%llu(id1=%d,id2=%d): cpu time=%0.6f, wall time=%0.6f\n", 
	     (unsigned long long)(nexthash-orignext), (unsigned long long)(end-orignext), hashpairs==hashpairs1 ? 1 : 2, (unsigned long long)(orignext-hashpairs), (unsigned long long)(end-hashpairs) -1, 
	     nummaps, NumPair, MaxPair, origrefid, refid, (unsigned long long)(nexthash-hashpairs), nexthash->id1, nexthash->id2, mtime(), wtime());
      fflush(stdout);
    }
  }
  if(DEBUG && nexthash < end && nexthash > orignext) assert(nexthash->id1 > nexthash[-1].id1);

  if(DEBUG>=2){
    for(int i = 0; i < NumPair; i++){
      int rid = YidList[i];
      assert(origrefid <= rid && rid <= refid);
    }
  }

  /* update refidend1 or refidend2 */
  if(hashpairs==hashpairs1){
    refidend1 = refid;
  } else {
    if(DEBUG) assert(hashpairs==hashpairs2);
    refidend2 = refid;
  }

  return refid;
}

static inline int CalignLogPVDec(Calign **pp1, Calign **pp2)
{
  double logPV1 = pp1[0]->logPV;
  double logPV2 = pp2[0]->logPV;

  return (logPV1 < logPV2) ? 1 : (logPV1 > logPV2) ? -1 : 0;
}

static inline int CalignMapidInc(Calign **pp1, Calign **pp2)
{
  int id1 = pp1[0]->mapid2;
  int id2 = pp2[0]->mapid2;

  return (id1 > id2) ? 1 : (id1 < id2) ? -1 : 0;
}

static inline int CalignScoreDec(Calign **pp1, Calign **pp2)
{
  double score1 = pp1[0]->score;
  double score2 = pp2[0]->score;

  return (score1 < score2) ? 1 : (score1 > score2) ? -1 : 0;
}

static int ErrorRawxInc(Cinterval *p1, Cinterval *p2)
{
  return (p1->rawx > p2->rawx) ? 1 : (p1->rawx < p2->rawx) ? -1 : 0;
}

/* Sort in ascending order of mapid, L */
static int ErrorIdInc(Cinterval *p1, Cinterval *p2)
{
  return (p1->mapid > p2->mapid) ? 1 : (p1->mapid < p2->mapid) ? -1 :
    (p1->L > p2->L) ? 1 : (p1->L < p2->L) ? -1 : 0;
}

static int DoubleInc(double *p1, double *p2)
{
  double val1 = *p1;
  double val2 = *p2;
  return (val1 > val2) ? 1 : (val1 < val2) ? -1 : 0;
}

static int CmapIdInc(Cmap **p1, Cmap **p2)
{
  /* Cannot return p1[0]->id - p2[0]->id since that may be out of range for integer */
  return (p1[0]->id > p2[0]->id) ? 1 : (p1[0]->id < p2[0]->id) ? -1 : 0;
}

#if 0
static int logPVdec(Calign **p1, Calign **p2)
{
  Calign *align1 = p1[0];
  Calign *align2 = p2[0];
  return (align1->logPV < align2->logPV) ? 1 : (align1->logPV > align2->logPV) ? -1 : 0;
}
#endif

class CmapPv {
public:
  int contigid;/**< index into gcontig->contig[] */
  double LogPv;
  int discard;
};

#if 0
static int mapPVinc(register CmapPv *p1, register CmapPv *p2)
{
  if(p1->LogPv > p2->LogPv)
    return 1;
  if(p1->LogPv < p2->LogPv)
    return -1;
  register int mapid1 = Gmap[gcontig->contig[p1->contigid].mapid]->mapid;
  register int mapid2 = Gmap[gcontig->contig[p2->contigid].mapid]->mapid;
  if(DEBUG>=2) assert(mapid1 != mapid2);
  return (mapid1 > mapid2) ? 1 : (mapid1 < mapid2) ? -1 : 0;
}
#endif

#if 0
static int errinc(Cinterval *p1, Cinterval *p2)
{
  return (p1->err > p2->err) ? 1 : (p1->err < p2->err) ? -1 : 0;
}
#endif

static TFLOAT ChimScore = MINSCORE;

/* compute log Likelihood for current res & resSD */
static double resLL(Cresdata *resdata, int numresdata, int numthreads, double *Tarray1, int c, int verb)
{
  assert(FIX_OUTLIER_MIS);
  if(VERB && verb){
    printf("resLL:numresdata=%d,c=%d,res=%0.6f,resSD[c]=%0.6f\n",numresdata,c,res[c],resSD[c]);
    fflush(stdout);
  }
  int lastmapid = -1, lastrefid = -1, lastor = -1;

  /* update score_init() parameters resKB and IresSD */
  resKB[c] = res[c] * PixelLen;
  IresSD[c] = 1.0/(resSD[c] * PixelLen * sqrt(2.0));

  double Lsum = 0.0;
  for(int tid = 0; tid < numthreads; tid++)
    Tarray1[tid] = 0.0;
      
  #pragma omp parallel num_threads(numthreads) if(numthreads > 1 && !verb)
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num ();
#endif

    double myLsum = 0.0;

    #pragma omp for schedule(static,64)
    for(int i = 0; i < numresdata; i++){
      Cresdata *pres = &resdata[i];
      Calign *align = pres->align;
      int refid = align->mapid1;
      int mapid = align->mapid2;
      if(DEBUG/* HERE >=2 */) assert(0 <= refid && refid < numrefmaps);
      if(DEBUG/* HERE >=2 */ && !(0 <= mapid && (mapid < nummaps || (startmaps <= mapid && mapid < totalmaps)))){
	printf("resdata[%d]:align->mapid1=%d,numrefmaps=%d\n",i,resdata[i].align->mapid1,numrefmaps);
	fflush(stdout);
	assert(0 <= mapid && (mapid < nummaps || (startmaps <= mapid && mapid < totalmaps)));
      }
      int orientation = align->orientation;
      if(VERB && verb && (mapid != lastmapid || refid != lastrefid)){
	if(lastmapid >= 0){
  	  printf("mapid=%d,refid=%d,or=%d:cum Lsum=%0.6f\n",
            lastmapid,lastrefid,lastor,myLsum);
	  fflush(stdout);
        }

	lastmapid = mapid;
	lastrefid = refid;
	lastor = orientation;
      }
      FLOAT *Y = refmap[refid]->site[c];
      double LL1 = pres->n ? 
	(REFDEBUG_STRICT ? 
	 Sint(pres->x, Yc(Y,pres->I,pres->K) - Yc(Y,pres->I - pres->K - pres->n, pres->T), pres->m, pres->n, pres->J, pres->I, pres->K, pres->T, Y, pres->outlier, c) : 
	 Sint(pres->x, Yc(Y,pres->I,pres->K) - Yc(Y,pres->I - pres->K - pres->n, pres->T), pres->m, pres->n, pres->J, pres->I, pres->K, pres->T, Y, c)) : 
	Sm(0,pres->I,pres->K,Y,c);
      double LL2 = -1e+10;
      LL1 = max(LL2,LL1);
      if(!REFDEBUG_STRICT && pres->Ip > 0){
	if(DEBUG) assert(pres->Ip == pres->I || pres->Ip == pres->I + 1);
	if(pres->n1){
	  double Sintp = pres->n ? 
	    Sint(pres->x, Yc(Y,pres->Ip,pres->K + 1) - Yc(Y,pres->I - pres->K - pres->n, pres->T), pres->m, pres->n - 1 + pres->Ip - pres->I, pres->J, pres->Ip, pres->K + 1, pres->T, Y,c) : 
	    Sm(0, pres->Ip, pres->K + 1, Y,c);
	  double Sint1 = Sint(pres->x1, Yc(Y,pres->I2,pres->K2) - Yc(Y,pres->I2 - pres->K2 - pres->n2, pres->T2), pres->m1, pres->n2, pres->J1, pres->I2, pres->K2, pres->T2, Y,c);
	  double Sint2 = Sint(pres->x1, Yc(Y,pres->I1,pres->K1) - Yc(Y,pres->I1 - pres->K1 - pres->n1, pres->T1), pres->m1, pres->n1, pres->J1, pres->I1, pres->K1, pres->T1, Y,c);
	  LL2 = max(LL2,Sintp) + max(LL2,Sint1) - max(LL2,Sint2);
	} else {
	  double Sm2 = Sm(0, pres->Ip, pres->K + 1, Y,c) + FnPenalty[c];
	  LL2 = max(LL2, Sm2);
	}
      }
      myLsum += max(LL1,LL2);

      LFLOAT Pen = 0.0;
      if(REFDEBUG_STRICT && pres->U == 0){/* duplicate code from logLR() for left end score */
	FLOAT *X = Gmap[mapid]->site[c];
	int M = Gmap[mapid]->numsite[c];
	int I = pres->I;
	int K = pres->K;
	int J = pres->J;
	int Lij = align->Lij1;
	int Lijx = align->Lij2;
	LFLOAT scale = ScaleDeltaBPP ? 1.0 : (align->scaleID > 1) ? ScaleFactor[align->scaleID] : 1.0;

	LFLOAT x = scale * ((orientation==0) ? X[J] : X[M+1] - X[M+1 - J]);
	LFLOAT xLijx_1 = (Lijx <= 0) ? 0.0 : scale * ((align->orientation==0) ? X[Lijx-1] : X[M+1] - X[M+1-(Lijx-1)]);

	int n = I - K - Lij + 1;
	int m = J;
	if(DEBUG>=2) assert(n >= 1);
	if(DEBUG>=2) assert(m >= 1);

	LFLOAT xB = max(LZERO, x - resKB2[c]);
	if(VERB>=2 && verb)
	  printf("\t M=%d,I=%d,K=%d,J=%d,Lij=%d,scale=%0.4f,x=%0.3f,n=%d,m=%d,xB=%0.3f:Yc(I,K)=%0.3f,",M,I,K,J,Lij,scale,x,n,m,xB,Yc(Y,I,K));

	if(align->Lend <= -2){/* end outlier */
	  Pen = ChimScore;
	} else if(ENDFIX && x <= Yc(Y,I,K) && Lij > 0 && Yc(Y,I,K)-Y[Lij-1] < x){ /* Sbnd(x,Y[I,K]-Y[Lij-1],n,m) : NOTE RES_VARIANCE is ignored */
	  LFLOAT y = Yc(Y,I,K)-Y[Lij-1];
	  LFLOAT var = QUADRATIC_VARIANCE ? (SF[c]*SF[c] + fabs(SD[c])*SD[c]*y + SR[c]*SR[c]*y*y) : (SF[c]*SF[c] + fabs(SD[c])*SD[c]*y);
	  LFLOAT err = x-y;
	  Pen = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[c]) */ - err*err*0.5/var;
	  if(VERB>=3 && mapid==6 && c==1 && verb){
            printf("Pen=%0.6f:xB=%0.4f,n=%d,m=%d,M=%d,I=%d,K=%d,J=%d,Lij=%d,Y[Lij-1]=%0.4f,Yc=%0.4f:err=%0.8f,var=%0.8f,x=%0.4f,y=%0.4f,Frate[c]=%0.8f,FnPenalty[c]=%0.8f,FpPenalty[c]=%0.8f\n",
                Pen,xB,n,m,M,I,K,J,Lij,Y[Lij-1],Yc(Y,I,K),err,var,x,y,Frate[c],FnPenalty[c],FpPenalty[c]);
	    if(orientation)
	      printf("    X[M+1]=%0.4f,X[M+1-J]=%0.4f\n",X[M+1],X[M+1-J]);
	    for(int i = max(1,Lij-10); i <= I; i++)
	      printf("Y[%d]= %0.4f\n",i,Y[i]);
	    fflush(stdout);
          }
	} else if(ENDFIX && extend && x >= Yc(Y,I,K) && Lijx > 0 && x - xLijx_1 < Yc(Y,I,K)){/* Sbnd(x-xLijx_1, Yc(Y,I,K), J+1-Lijx, I-K, c) */
	  LFLOAT y = Yc(Y,I,K);
	  LFLOAT var = SF[c]*SF[c] + fabs(SD[c])*SD[c]*y;
	  if(QUADRATIC_VARIANCE)
	    var += SR[c]*SR[c]*y*y;
	  x = x - xLijx_1;
	  xB = max(LZERO,x-resKB2[c]);
	  LFLOAT err = x - y;
	  m = J+1-Lijx;
	  n = I-K;
	  if(DEBUG>=2) assert(m >= 1);
	  if(DEBUG>=2) assert(n >= 1);
	  Pen = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[0]) */ - err*err*0.5/var;
	} else { /* Send(min(x,Yc(Y,I,K)),J-1-max(1,Lijx),I-K+1 - max(1,Lij),c) */
	  if(extend){
	    x = min(x,Yc(Y,I,K));
	    xB = max(LZERO,x-resKB2[c]);
	    m = J+1-max(1,Lijx);
	    n = I-K+1 - max(1,Lij);
	    if(DEBUG>=2) assert(m >= 1);
	    if(DEBUG>=2) assert(n >= 1);
	  }
	  Pen = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[c]) */;
	  if(VERB>=3 && mapid==6 && c==1 && verb){
            printf("Pen=%0.6f:xB=%0.4f,n=%d,m=%d,M=%d,I=%d,K=%d,J=%d,Lij=%d,Y[Lij-1]=%0.4f,Yc=%0.4f:Frate[c]=%0.8f,FnPenalty[c]=%0.8f,FpPenalty[c]=%0.8f\n",
               Pen,xB,n,m,M,I,K,J,Lij,Y[Lij-1],Yc(Y,I,K),Frate[c],FnPenalty[c],FpPenalty[c]);
	    fflush(stdout);
          }
	}
      }
      if(REFDEBUG_STRICT && pres->U == align->numpairs - 1){/* duplicate code from logLR() for right end score */
	FLOAT *X = Gmap[mapid]->site[c];
	int M = Gmap[mapid]->numsite[c];
	int N = refmap[refid]->numsite[c];
	int H = pres->I;
	int D = pres->K;
	int G = pres->J;
	int Rij = align->Rij1;
	int Rijx = align->Rij2;
	LFLOAT scale = ScaleDeltaBPP ? 1.0 : (align->scaleID > 1) ? ScaleFactor[align->scaleID] : 1.0;
	
	LFLOAT x = scale * ((align->orientation==0) ? X[M+1]-X[G] : X[M + 1 - G]);
	LFLOAT xRijx_1 = (Rijx > M) ? 0.0 : scale * ((align->orientation==0) ? X[M+1] - X[Rijx+1] : X[M+1-(Rijx+1)]);

	int n = Rij - H + 1;
	int m = M + 1 - G;
	if(DEBUG>=2) assert(n >= 1);
	if(DEBUG>=2) assert(m >= 1);

	LFLOAT xB = max(LZERO,x-resKB2[c]);
	if(align->Rend <= -2){/* end outier */
	  Pen = ChimScore;
	} else if(ENDFIX && x <= Y[N+1]-Yc(Y,H,D) && Rij <= N && Y[Rij+1]-Yc(Y,H,D) < x){  /* Sbnd(x,Y[Rij+1]-Yc(Y,H,D),M+1-G,Rij+1-I,c) */
	  LFLOAT y = Y[Rij+1]-Yc(Y,H,D);
	  LFLOAT var = SF[c]*SF[c] + fabs(SD[c])*SD[c]*y + SR[c]*SR[c]*y*y;
	  LFLOAT err = x-y;
	  Pen = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[c]) */ - err*err*0.5/var;
	  if(VERB>=3 && mapid==6 && c==1 && verb){
            printf("Pen=%0.6f:xB=%0.4f,n=%d,m=%d,N=%d,M=%d,H=%d,D=%d,G=%d,Rij=%d:err=%0.8f,var=%0.8f,x=%0.4f,y=%0.4f,Frate[c]=%0.8f,FnPenalty[c]=%0.8f,FpPenalty[c]=%0.8f\n",
		   Pen,xB,n,m,N,M,H,D,G,Rij,err,var,x,y,Frate[c],FnPenalty[c],FpPenalty[c]);
	    fflush(stdout);
          }
	} else if(ENDFIX && extend && x >= Y[N+1] - Yc(Y,H,D) && Rijx <= M && x - xRijx_1 < Y[N+1]-Yc(Y,H,D)){/* Sbnd(x-xRijx_1, Y[N+1]-Yc(Y,H,D), Rijx+1-G,min(N,Rij)+1-H,c) */
	  LFLOAT y = Y[N+1] - Yc(Y,H,D);
	  LFLOAT var = SF[c]*SF[c] + fabs(SD[c])*SD[c]*y;
	  if(QUADRATIC_VARIANCE)
	    var += SR[c]*SR[c]*y*y;
	  x = x - xRijx_1;
	  xB = max(LZERO,x-resKB2[c]);
	  LFLOAT err = x - y;
	  m = Rijx + 1 - G;
	  n = min(N,Rij) + 1 - H;
	  if(DEBUG>=2) assert(m >= 1);
	  if(DEBUG>=2) assert(n >= 1);
	  Pen = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[0]) */ - err*err*0.5/var;
        } else {
	  if(extend){
	    x = min(x, Y[N+1]-Yc(Y,H,D));
	    xB = max(LZERO,x-resKB2[c]);
	    m = min(M,Rijx)+1-G;
	    n = min(N,Rij)+1-H;
	    if(DEBUG>=2) assert(m >= 1);
	    if(DEBUG>=2) assert(n >= 1);
	  }
	  Pen = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[c]) */;
	  if(VERB>=3 && mapid==6 && c==1 && verb){
            printf("Pen=%0.6f:xB=%0.4f,n=%d,m=%d,N=%d,M=%d,H=%d,D=%d,G=%d,Rij=%d:Frate[c]=%0.8f,FnPenalty[c]=%0.8f,FpPenalty[c]=%0.8f\n",
                 Pen,xB,n,m,N,M,H,D,G,Rij,Frate[c],FnPenalty[c],FpPenalty[c]);
	    fflush(stdout);
          }
        }
      }

      if(REFDEBUG_STRICT)
	myLsum += Pen;

      if((VERB>=2 && verb && mapid==6 && c==1) || (DEBUG>=2 && (!(isfinite(LL2) || (LL2 < 0.0 && isinf(LL2))) || (!(isfinite(LL1) || (LL1 < 0.0 && isinf(LL1)))) || !isfinite(Pen) || !(isfinite(myLsum) || (myLsum < 0.0 && isinf(myLsum)))))){
	#pragma omp critical
	{
	  printf("mapid=%d,refid=%d,or=%d,U=%d:I=%d,K=%d,n=%d",
		 mapid,refid,orientation,pres->U,pres->I,pres->K,pres->n);
	  double sm =  Sm(0,pres->I,pres->K,Y,c);
	  if(pres->n){
	    double sint = (REFDEBUG_STRICT ? 
			   Sint(pres->x, Yc(Y,pres->I,pres->K) - Yc(Y,pres->I - pres->K - pres->n, pres->T), pres->m, pres->n, pres->J, pres->I, pres->K, pres->T, Y, pres->outlier, c) : 
			   Sint(pres->x, Yc(Y,pres->I,pres->K) - Yc(Y,pres->I - pres->K - pres->n, pres->T), pres->m, pres->n, pres->J, pres->I, pres->K, pres->T, Y, c));

	    printf("(T=%d,J=%d,m=%d,x=%0.3f,y=%0.3f(%0.3f ... %0.3f,%0.3f),out=%d,Sint=%0.6f, Sm=%0.6f)",
		   pres->T,pres->J,pres->m,pres->x,Yc(Y,pres->I,pres->K) - Yc(Y,pres->I - pres->K - pres->n, pres->T),
		   Yc(Y,pres->I - pres->K - pres->n, pres->T), Y[pres->I - pres->K], Y[pres->I],pres->outlier,sint,sm);
	    fflush(stdout);
	    if(DEBUG) assert(isfinite(sint) || (!isnan(sint) && sint >= -1e+10));
	  } else
	    printf("(Sm=%0.6f)", sm);
	  printf("\n\tIp=%d",pres->Ip);
	  if(!REFDEBUG_STRICT && pres->Ip){
	    if(!pres->n1)
	      printf("(Sm=%0.6f,FnPen=%0.6f)",Sm(0, pres->Ip, pres->K + 1, Y,c), FnPenalty[c]);
	    else
	      printf(",I1=%d,K1=%d,n1=%d,T1=%d,J1=%d,m1=%d,x1=%0.3f,y1=%0.3f,I2=%d,K2=%d,n2=%d,T2=%d,y2=%0.3f(Sintp=%0.6f,Sint1=%0.6f,Sint2=%0.6f):LL2=%0.8f,",
		     pres->I1,pres->K1,pres->n1,pres->T1,pres->J1,pres->m1,pres->x1,Yc(Y,pres->I1,pres->K1) - Yc(Y,pres->I1 - pres->K1 - pres->n1, pres->T1),
		     pres->I2,pres->K2,pres->n2,pres->T2,Yc(Y,pres->I2,pres->K2) - Yc(Y,pres->I2 - pres->K2 - pres->n2, pres->T2),
                     (pres->n ? 
		      Sint(pres->x, Yc(Y,pres->Ip,pres->K + 1) - Yc(Y,pres->I - pres->K - pres->n, pres->T), pres->m, pres->n - (pres->I - pres->Ip + 1), pres->J, pres->I, pres->K + 1, pres->T, Y,c) : 
		      Sm(0, pres->Ip, pres->K + 1, Y,c)),
		     Sint(pres->x1, Yc(Y,pres->I1,pres->K1) - Yc(Y,pres->I1 - pres->K1 - pres->n1, pres->T1), pres->m1, pres->n1, pres->J1, pres->I1, pres->K1, pres->T1, Y,c),
		     Sint(pres->x1, Yc(Y,pres->I2,pres->K2) - Yc(Y,pres->I2 - pres->K2 - pres->n2, pres->T2), pres->m1, pres->n2, pres->J1, pres->I2, pres->K2, pres->T2, Y,c), LL2); 
	  }
	  printf(":LL1=%0.6f:LL=%0.6f,Pen=%0.6f(delta=%0.6f):LLsum=%0.6f (score=%0.6f,logPV=%0.2f,numpairs=%d)\n",LL1,max(LL1,LL2),Pen,max(LL1,LL2)+Pen,myLsum, align->score, align->logPV, align->numpairs);
	  fflush(stdout);
	  if(DEBUG/* HERE >=2 */) assert(isfinite(LL1) || (LL1 < 0.0 && isinf(LL1)));
	  if(DEBUG/* HERE >=2 */) assert(isfinite(LL2) || (LL2 < 0.0 && isinf(LL2)));
	  if(DEBUG/* HERE >=2 */) assert(isfinite(Pen));
	  if(DEBUG/* HERE >=2 */) assert(isfinite(myLsum) || (myLsum < 0.0 || isinf(myLsum)));
	}
      }
    }
    
    if(VERB && verb && lastmapid >= 0){/* output last molecule summary */
      printf("mapid=%d,refid=%d,or=%d:cum Lsum=%0.6f\n",
        lastmapid,lastrefid,lastor,myLsum);
      fflush(stdout);
    }

    if(DEBUG) assert(isfinite(myLsum) || (myLsum < 0.0 || isinf(myLsum)));
    Tarray1[tid] = myLsum;
  }
  qsort(Tarray1, numthreads, sizeof(double), (intcmp *)DoubleInc);
  for(int tid = 0; tid < numthreads; tid++)
    Lsum += Tarray1[tid];
  
  if(VERB>=2 || verb){
    printf("resLL:numresdata=%d,res=%0.6f,resSD[c]=%0.6f:LLsum=%0.6f\n",numresdata,res[0],resSD[0],Lsum);
    fflush(stdout);
  }
  return Lsum;
}

/* return true if number of EndOutliers <= AlignedEndoutlierThreshold */
static inline bool AlignedEndOutlier(Calign *align)
{
  int EndOutlierCnt = ((align[0].Lend <= -2 || align[1].Lend <= -2) ? 1 : 0) + ((align[0].Rend <= -2 || align[1].Rend <= -2) ? 1 : 0);
  return (EndOutlierCnt <= AlignedEndOutlierThreshold);
}

static inline FLOAT AlignedLength(Calign *align, FLOAT **YY)
{
  //  if(AlignedLengthThreshold <= 0.0)
  //    return 1.0;
  int U1 = align[1].numpairs;
  int U2 = align[2].numpairs;
  if(DEBUG) assert(U1 >= 0);
  if(DEBUG) assert(U2 >= 0);
  FLOAT len1 = YY[0][align[1].sites1[U1 - 1]] - YY[0][align[1].sites1[0]];
  FLOAT len2 = YY[1][align[2].sites1[U2 - 1]] - YY[1][align[2].sites1[0]];
  if(DEBUG) assert((U1 > 1) ? (len1 > 0.0) : (len1 >= 0.0));
  if(DEBUG) assert((U2 > 1) ? (len2 > 0.0) : (len2 >= 0.0));
  return max(len1,len2);
}

/* return true of all alignment thresholds are satisfied */
static int AlignedThreshold(Calign *align, FLOAT **YY, double ScoreThresh = ScoreThreshold, double LogPvThresh = LogPvThreshold)
{
  return (align->numpairs >= AlignedSiteThreshold && align->score > ScoreThresh && align->logPV > LogPvThresh && 
	  max(align[1].maxoutlier,align[2].maxoutlier) <= AlignedOutlierThreshold && AlignedEndOutlier(&align[1]) && AlignedLength(align,YY) >= AlignedLengthThreshold) ? 1 : 0;
}

static int maxerrors[MAXCOLOR+1],numerrors[MAXCOLOR+1];
static Cinterval *errors[MAXCOLOR+1];

static int maxresdata[MAXCOLOR], numresdata[MAXCOLOR];
static Cresdata *resdata[MAXCOLOR] = {NULL,NULL};

static void appenderror(int c,/**< color (or colors+1 for misalignment term) */
			double x,
			double y,
			double len,
			double resvar,
			int mapid,
			int L, int R,
			int end,
			int stitch)
{
  if(numerrors[c] >= maxerrors[c]){
    int newmax = max(1024,maxerrors[c]*2);
    Cinterval *newerrors = new Cinterval[newmax];
    for(int i=0;i<numerrors[c];i++)
      newerrors[i] = errors[c][i];
    delete [] errors[c];
    errors[c] = newerrors;
    maxerrors[c] = newmax;
  }
  Cinterval *perr = &errors[c][numerrors[c]++];
  perr->x = x;
  perr->y = y;
  perr->len = len;
  if(RES_VARIANCE)
    perr->resvar = resvar;
  perr->mapid = mapid;
  perr->L = L;
  perr->R = R;
  perr->end = end;
  if(STITCH)
    perr->stitch = stitch;
}
			
static void appendres(int c,
		      Calign *align,
		      FLOAT *X, FLOAT *Y,
		      int M, int N, 
		      int I, int K, int J, 
		      int H, int D, int G, int T,
		      FLOAT x, double escale, int outlier
		      )
{
  if(numresdata[c] >= maxresdata[c]){
    int newmax = max(1024,numresdata[c]*2);
    Cresdata *newresdata = new Cresdata[newmax];
    memcpy(newresdata,resdata[c],numresdata[c]*sizeof(Cresdata));
    delete [] resdata[c];
    resdata[c] = newresdata;
    maxresdata[c] = newmax;
  }
  Cresdata *pres = &resdata[c][numresdata[c]++];

  pres->align = &align[c];
  if(DEBUG>=2) assert(0 <= align[-1].mapid1 && align[-1].mapid1 < numrefmaps);
  if(DEBUG>=2) assert(0 <= align[-1].mapid2 && (align[-1].mapid2 < nummaps || (NoSplit <= 1 && startmaps <= align[-1].mapid2 && align[-1].mapid2 < totalmaps)));
  if(VERB>=2){
    printf("c=%d:resdata[c][%d]:&align[c]=%p:align[-1].mapid1=%d,align[-1].mapid2=%d\n",c,numresdata[c]-1,pres->align,align[-1].mapid1,align[-1].mapid2);
    fflush(stdout);
  }

  pres->I = I;
  pres->K = K;
  pres->J = J;// to support REFDEBUG
  pres->U = T;

  FLOAT delY1 = (I-K > 1) ? Y[I-K] - Y[I-K-1] : 9999.0;
  FLOAT delY2 = (I < N) ? Y[I+1]-Y[I] : 9999.0;

  if(T <= 0){/* left end  align site */

    pres->n = 0;

    if(delY1 < delY2){/* expand deres interval Y[I-K..I] to the left to Y[I-K-1..I] */
      pres->Ip = I;
      pres->n1 = 0;
    } else { /* expand deres interval Y[I-K..I] to the right to Y[I-K .. I+1] */
      if(align[c].numpairs <= T+1){/* no further aligned sites */
	if(I >= N)/* cannot expand deres interval to right, since there are no more sites */
	  pres->Ip = 0;
	else {
	  pres->Ip = I+1;
	  pres->n1 = 0;
	}
      } else {/* next aligned site is present */
	int In = align[c].sites1[1];
	int Kn = align[c].sitesK1[1];
	int Jn = align[c].sites2[1];
	if(I+1 >= In-Kn) /* cannot expand deres interval to right, since there is no unaligned site before next aligned site */
	  pres->Ip = 0;
	else {
	  pres->Ip = I+1;
	  
	  /* original resolved interval */
	  pres->I1 = In;
	  pres->K1 = Kn;
	  pres->n1 = In - Kn - I;
	  pres->T1 = K;
	  pres->J1 = Jn;
	  pres->m1 = Jn - J;
	  pres->x1 = escale * (align[c].orientation ? X[M+1-J] - X[M+1-Jn] : X[Jn] - X[J]);
	  
	  /* alternate resolved interval */
	  pres->I2 = In;
	  pres->K2 = Kn;
	  pres->n2 = In - Kn - I - 1;
	  pres->T2 = K + 1;
	}
      }
    }
  } else {/* right aligned site of an aligned interval */

    pres->m = J-G;
    pres->n = I-K-H;
    pres->T = D;
    pres->x = x;
    if(REFDEBUG_STRICT)
      pres->outlier = outlier;

    if(delY1 < delY2){/* try to expand deres interval Y[I-K..I] to the left to Y[I-K-1..I] */
      if(I-K-1 >= H) /* cannot expand deres interval to left, since there is no unaligned site before previous aligned site */
	pres->Ip = 0;
      else {
	pres->Ip = I;
	
	/* original resolved interval */
	pres->I1 = I;
	pres->K1 = K;
	pres->n1 = I - K - H;
	pres->T1 = D;
	pres->J1 = J;
	pres->m1 = J - G;
	pres->x1 = x;
	
	/* alternate resolved interval */
	pres->I2 = I;
	pres->K2 = K +1;
	pres->n2 = I - K - 1 - H;
	pres->T2 = D;
      }
    } else {
      if(align[c].numpairs <= T + 1){/* no further aligned sites */
	if(I >= N)/* cannot expand deres interval to right, since there are no more sites */
	  pres->Ip = 0;
	else {
	  pres->Ip = I + 1;
	  pres->n1 = 0;
	}
      } else {/* next aligned site is present */
	int In = align[c].sites1[T+1];
	int Kn = align[c].sitesK1[T+1];
	int Jn = align[c].sites2[T+1];
	if(I + 1 >= In - Kn) /* cannot expand deres interval to right, since there is no unaligned site before next aligned site */
	  pres->Ip = 0;
	else {
	  pres->Ip = I + 1;
	  
	  /* original resolved interval Y[I .. In - Kn] */
	  pres->I1 = In;
	  pres->K1 = Kn;
	  pres->n1 = In - Kn - I;
	  pres->T1 = K;
	  pres->J1 = Jn;
	  pres->m1 = Jn - J;
	  pres->x1 = escale * (align[c].orientation ? X[M+1-J] - X[M+1-Jn] : X[Jn] - X[J]);
	  
	  /* alternate resolved interval Y[I+1 .. In-Kn] */
	  pres->I2 = In;
	  pres->K2 = Kn;
	  pres->n2 = In - Kn - I - 1;
	  pres->T2 = K + 1;
	}
      }
    }
  }
}

class DP {
public:
  TFLOAT scoreN;/**< best alignment score that does NOT have a left unaligned end (L at left end must be >= -1 cf PoutlierEnd == 0) */
  TFLOAT scoreEnd;/**< alignment score for best alignment between Y[0..I] & X[0..J] with Y[I-T..I] aligned with X[J] (LEnd at left end must be <= -2) */
  TFLOAT Lscore;/**< best left hand alignment end score with Y[I-T..I] aligned with X[J] with type (L) >= -1 */
  TFLOAT LscoreEnd;/**< best Lscore with type (LEnd) <= -2 */
  TFLOAT Rscore;/**< best right hand alignment end score with Y[I-T..I] aligned with X[J] with type (R) >= -1 */
  TFLOAT RscoreEnd;/**< best Rscore with type (REnd) <= -2 */

  int G;/**< previous site in Y in best alignment (G <= 0 if none : value is type of left end score. In this case scoreN == Lscore && G==L) */
  int T;/**< sites Y[G-T .. G] were merged to align with X[H] (If resSD>0) */
  int H;/**< previous site in X in best alignment (if G > 0) */

  int GEnd,TEnd,HEnd;/**< previous site(s) for scoreEnd (GEnd <= -2 if none : value is type of left end score. In this case scoreEnd == LscoreEnd && GEnd == LEnd) */

  int Lij;/**< leftmost site/end in Y that overlaps X when Y[I-K..I] is aligned with X[J] (for debugging) */
  int Rij;/**< rightmost site/end in Y that overlaps X when Y[I-K..I] is aligned with X[J] */
  int Lijx;/**< same as Lij for site/end in Y : only used if extend==1 */
  int Rijx;/**< same as Rij for site/end in Y : only used if extend==1 */

  int L;/**< type of left hand alignment end score Lscore (>= -1) */
  int LEnd;/**< type of left hand alignment end score LscoreEnd ( <= -2) */
  int R;/**< type of right hand alignment end score Rscore (>= -1)*/
  int REnd;/**< type of right hand alignment end score RscoreEnd ( <= -2) */
};

/* Flat Array based Dynamic programming data */
class AL2 {
#if FLAT_MEMORY >= 1

public:

  long long Kstride;/**< distance between K index increments == maxN*maxM */
  long long Istride;/**< distance between I index increments == M */
  int kmax;/**< largest value of Kmax[1..N] */
  /* NOTE : the J index stride is always 1 */
  
#define ARRAY(type,field) type *field##_base; inline type& field(int I, int K, int J) { return field##_base[I*Istride + K*Kstride + J];}

  ARRAY(TFLOAT, scoreN);/**< best alignment score that does NOT have a left unaligned end (L at left end must be >= -1 cf PoutlierEnd == 0) */
  ARRAY(int,G);/**< previous site in Y in best alignment (G <= 0 if none : value is type of left end score. In this case scoreN == Lscore && G==L) */
  ARRAY(int,T);/**< sites Y[G-T .. G] were merged to align with X[H] (If resSD>0) */
  ARRAY(int,H);/**< previous site in X in best alignment (if G > 0) */

  ARRAY(TFLOAT, scoreEnd);/**< alignment score for best alignment between Y[0..I] & X[0..J] with Y[I-T..I] aligned with X[J] (LEnd at left end must be <= -2) */
  ARRAY(int, GEnd); ARRAY(int, TEnd); ARRAY(int, HEnd);/**< previous site(s) for scoreEnd (GEnd <= -2 if none : value is type of left end score. In this case scoreEnd == LscoreEnd && GEnd == LEnd) */

  ARRAY(int, Lij);/**< leftmost site/end in Y that overlaps X when Y[I-K..I] is aligned with X[J] (for debugging) */
  ARRAY(int, Rij);/**< rightmost site/end in Y that overlaps X when Y[I-K..I] is aligned with X[J] */
  ARRAY(int, Lijx);/**< same as Lij for site/end in Y : only used if extend==1 */
  ARRAY(int, Rijx);/**< same as Rij for site/end in Y : only used if extend==1 */
  ARRAY(TFLOAT, Lscore);/**< best left hand alignment end score with Y[I-T..I] aligned with X[J] with type (L) >= -1 */
  ARRAY(int, L);/**< type of left hand alignment end score Lscore */
  ARRAY(TFLOAT, LscoreEnd);/**< best Lscore with type (LEnd) <= -2 */
  ARRAY(int, LEnd);
  ARRAY(TFLOAT, Rscore);/**< best right hand alignment end score with Y[I-T..I] aligned with X[J] with type (R) >= -1 */
  ARRAY(int, R);/**< type of right hand alignment end score Rscore */
  ARRAY(TFLOAT, RscoreEnd);/**< best Rscore with type (REnd) <= -2 */
  ARRAY(int, REnd);/**< type of right hand alignment end score RscoreEnd */

#undef ARRAY

  AL2(){
    scoreN_base = scoreEnd_base = Lscore_base = LscoreEnd_base = Rscore_base = RscoreEnd_base = NULL;
    G_base = T_base = H_base = GEnd_base = TEnd_base = HEnd_base = Lij_base = Rij_base = Lijx_base = Rijx_base = L_base = LEnd_base = R_base = REnd_base = NULL;
  }

#else // FLAT_MEMORY==0

public:
  DP ***array;
  AL2(){
    array = NULL;
  }
#define ARRAY(type, field) inline type& field(int I, int K, int J) { return array[I][K][J].field;}
  
  ARRAY(TFLOAT, scoreN);/**< best alignment score that does NOT have a left unaligned end (L at left end must be >= -1 cf PoutlierEnd == 0) */
  ARRAY(int,G);/**< previous site in Y in best alignment (G <= 0 if none : value is type of left end score. In this case scoreN == Lscore && G==L) */
  ARRAY(int,T);/**< sites Y[G-T .. G] were merged to align with X[H] (If resSD>0) */
  ARRAY(int,H);/**< previous site in X in best alignment (if G > 0) */

  ARRAY(TFLOAT, scoreEnd);/**< alignment score for best alignment between Y[0..I] & X[0..J] with Y[I-T..I] aligned with X[J] (LEnd at left end must be <= -2) */
  ARRAY(int, GEnd); ARRAY(int, TEnd); ARRAY(int, HEnd);/**< previous site(s) for scoreEnd (GEnd <= -2 if none : value is type of left end score. In this case scoreEnd == LscoreEnd && GEnd == LEnd) */

  ARRAY(int, Lij);/**< leftmost site/end in Y that overlaps X when Y[I-K..I] is aligned with X[J] (for debugging) */
  ARRAY(int, Rij);/**< rightmost site/end in Y that overlaps X when Y[I-K..I] is aligned with X[J] */
  ARRAY(int, Lijx);/**< same as Lij for site/end in Y : only used if extend==1 */
  ARRAY(int, Rijx);/**< same as Rij for site/end in Y : only used if extend==1 */
  ARRAY(TFLOAT, Lscore);/**< best left hand alignment end score with Y[I-T..I] aligned with X[J] with type (L) >= -1 */
  ARRAY(int, L);/**< type of left hand alignment end score Lscore */
  ARRAY(TFLOAT, LscoreEnd);/**< best Lscore with type (LEnd) <= -2 */
  ARRAY(int, LEnd);
  ARRAY(TFLOAT, Rscore);/**< best right hand alignment end score with Y[I-T..I] aligned with X[J] with type (R) >= -1 */
  ARRAY(int, R);/**< type of right hand alignment end score Rscore */
  ARRAY(TFLOAT, RscoreEnd);/**< best Rscore with type (REnd) <= -2 */
  ARRAY(int, REnd);/**< type of right hand alignment end score RscoreEnd */

#undef ARRAY  

#endif // FLAT_MEMORY==0
};

static long long DPsizF = 6;/* number of TFLOAT fields in DP */
static long long DPsizI = 14;/* number of int fields in DP */
/* NOTE :  DPsizF*sizeof(TFLOAT) + DPsizI*sizeof(int) should match sizeof(DP) */

static int chimpaircnt;/* count how many chimeric pairs are next to each other (within 20 * Ylambda) */

static int giter;

static long long pcnt[MAXCOLOR] = {0,0};/**< sum over I=1..N[c] of KKmax[c][I]+1  (for debugging memory allocation) */

/** extract alignment with right end at Y[I],X[J] from A[][] and call logFP() : version with resSD > 0.0 */
static double alignFPsd(AL2 *A, int I, int J,int K,
			FLOAT *Y, FLOAT *X,
			int N, int M, int R,int L,
			int orientation, int refid, int mapid, int c, int verb)
{
  /* allocate space to store alignment in reverse order */
  int *Ilist = new int[M];
  int *Klist = new int[M];
  int *Jlist = new int[M];
  int *outlier = new int[M];
  int U = 0;
  //  int origI = I, origJ = J, origK = K;

  while(I > 0){
    if(DEBUG){
      assert(I-K >= 1 && I <= N);
      if(!(J >= 1 && J <= M)){
        #pragma omp critical
	{
	  printf("alignFPsd,mapid=%d,or=%d,c=%d:U=%d,I=%d,K=%d,J=%d,M=%d,N=%d\n",
		 mapid,orientation,c,U,I,K,J,M,N);
	  fflush(stdout);
	}
	assert(J >= 1 && J <= M);
      }
    }
    //    DP *p = &A[I][K][J];

    Ilist[U] = I;
    Klist[U] = K;
    Jlist[U] = J;
    if(L >= -1){
      int G = A->G(I,K,J);
      int T = A->T(I,K,J);
      int H = A->H(I,K,J);
      if(G > 0){
	if(DEBUG) assert(T >= 0);
	if(DEBUG) assert(H > 0);
	if(DEBUG) assert(U+1 < M);
	double x = X[J] - X[H];
	double y = Yc(Y,I,K) - Yc(Y,G,T);
	if(maptype){
	  TFLOAT Bias,Pen,PenSm;
	  SintDetail(x,y, J - H, I - K - G,J,I,K,T,PRtab[c][refid],Y,c,Bias,Pen,PenSm,0);
	  //	  double iscore = p->scoreN - A[p->G][p->T][p->H].scoreN;
	  RFLOAT OutPen = OutlierPenalty[c];
	  if(outlierBC)
	    OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-K-G];
	  if(OUTLIER_LTYPE==0)
	    OutPen -= (x+y) * OutlierLambdaInv;
	  else
	    OutPen -= fabs(x-y) * OutlierLambdaInv;
	  TFLOAT iscore = Bias + OutlierBias + (FIX_OUTLIER_MIS ? max(Pen,OutPen) + PenSm : max(Pen+PenSm, OutPen));
	  TFLOAT outscore = Bias + Pen + PenSm;
	  outlier[U+1] = (iscore > outscore + (TFLOAT)0.01) ? 1 : 0;
	} else if(OUTLIER_DELTA(x-y)){
	  TFLOAT Bias,Pen,PenSm;
	  SintDetail(x,y, J - H, I - K - G,J,I,K,T,PRtab[c][refid],Y,c,Bias,Pen,PenSm,0);
	  //	  double iscore = p->scoreN - A[p->G][p->T][p->H].scoreN;
	  TFLOAT OutPen = OutlierPenalty[c];
	  if(outlierBC)
	    OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-K-G];
	  OutPen -= fabs(x-y) * OutlierLambdaInv;
	  TFLOAT iscore = Bias + OutlierBias + (FIX_OUTLIER_MIS ? max(Pen,OutPen) + PenSm : max(Pen+PenSm, OutPen));
	  TFLOAT outscore = Bias + Pen + PenSm;
	  outlier[U+1] = (iscore > outscore + (TFLOAT)0.01) ? 1 : 0;
	} else
	  outlier[U+1] = 0;
      }
      I = G;
      K = T;
      J = H;
    } else {/* L <= -2 */
      int GEnd = A->GEnd(I,K,J);
      int TEnd = A->TEnd(I,K,J);
      int HEnd = A->HEnd(I,K,J);
      if(GEnd > 0){
	if(DEBUG) assert(TEnd >= 0);
	if(DEBUG) assert(HEnd > 0);
	if(DEBUG) assert(U+1 < M);
	FLOAT x = X[J] - X[HEnd];
	FLOAT y = Yc(Y,I,K) - Yc(Y,GEnd,TEnd);
	if(maptype){
	  TFLOAT Bias,Pen,PenSm;
	  SintDetail(x,y, J- HEnd,I-K- GEnd,J,I,K,TEnd,PRtab[c][refid],Y,c,Bias,Pen,PenSm,0);
	  //	  double iscore = A->scoreEnd(I,K,J) - A->scoreend(GEnd,TEnd,HEnd);
	  TFLOAT OutPen = OutlierPenalty[c];
	  if(outlierBC)
	    OutPen += OutlierPenaltyBC[J- HEnd] + OutlierPenaltyBC[I-K- GEnd];
	  if(OUTLIER_LTYPE==0)
	    OutPen -= (x+y) * OutlierLambdaInv;
	  else
	    OutPen -= fabs(x-y) * OutlierLambdaInv;

	  TFLOAT iscore = Bias + OutlierBias + (FIX_OUTLIER_MIS ? max(Pen,OutPen) + PenSm : max(Pen+PenSm,OutPen));
	  TFLOAT outscore = Bias + Pen + PenSm;
	  outlier[U+1] = (iscore > outscore + (TFLOAT)0.01) ? 1 : 0;
	} else if(OUTLIER_DELTA(x-y)){
	  TFLOAT Bias,Pen,PenSm;
	  SintDetail(x,y, J- HEnd,I-K- GEnd,J,I,K,TEnd,PRtab[c][refid],Y,c,Bias,Pen,PenSm,0);
	  //	  double iscore = A->scoreEnd(I,K,J) - A->scoreend(GEnd,TEnd,HEnd);
	  TFLOAT OutPen = OutlierPenalty[c];
	  if(outlierBC)
	    OutPen += OutlierPenaltyBC[J- HEnd] + OutlierPenaltyBC[I-K- GEnd];
	  OutPen -= fabs(x-y) * OutlierLambdaInv;

	  TFLOAT iscore = Bias + OutlierBias + (FIX_OUTLIER_MIS ? max(Pen,OutPen) + PenSm : max(Pen+PenSm,OutPen));
	  TFLOAT outscore = Bias + Pen + PenSm;
	  outlier[U+1] = (iscore > outscore + (TFLOAT)0.01) ? 1 : 0;
	} else
	  outlier[U+1] = 0;
      }
      I = GEnd;
      K = TEnd;
      J = HEnd;
    }
    U++;
  }
  if(DEBUG) assert(U>0);
    
  /* allocate alignment fragment sizes */
  FLOAT *Xfrag = new FLOAT[U-1];
  FLOAT *Yfrag = new FLOAT[U-1];
  int misscnt = 0;
  double falsecnt = 0.0;
  int OutlierCnt = 0;/* number of internal outliers */
  double OutlierMN = 1.0;/* product over all outlier of m x n (the number of intervals in the internal outlier aligned segment) */

  int G= I, H= J, D= K;
  int F = 0;/* segment (aligned interval) count */
  int EndOutlierCnt = 0;
  for(int T=0;T < U;T++, G = I, H = J, D = K){
    I = Ilist[U-1-T];
    K = Klist[U-1-T];
    J = Jlist[U-1-T];
    if(PVres >= 2 && K > 0)
      for(int k = 0; k < K; k++)
	falsecnt += Pr(Y[I-k] - Y[I-k-1],c);
    if(T>0){
      if(DEBUG)
	assert(G > 0 && H > 0);
      if(outlier[U-T]){
	if(OUTLIER_PV){
	  OutlierCnt++;
	  int m = J-H;
	  int n = I-K-G;
	  OutlierMN *= (double)(m*n);
	  if(DEBUG>=2) assert(OutlierMN > 0);
	}	
	continue;
      }
      Yfrag[F] = Yc(Y,I,K) - Yc(Y,G,D);
      Xfrag[F] = X[J] - X[H];
      if(DEBUG>=2) assert(Yfrag[F] > 0.0);
      if(DEBUG>=2) assert(Xfrag[F] > 0.0);
      F++;
      misscnt += I-K-G-1;
      falsecnt += J-H-1;
    } else {/* left end */
      int Lij = (G <= -2) ? I-K : max(1,A->Lij(I,K,J));
      int Lijx = (G <= -2) ? J : extend ? max(1,A->Lijx(I,K,J)) : 1;
      if(VERB>=2 && verb){
	printf("I=%d,K=%d,J=%d:G=%d,Lij=%d,Lijx=%d,A[I][K][J].Lij=%d,A[I][K][J].Lijx=%d\n",
	       I,K,J,G,Lij,Lijx,A->Lij(I,K,J),A->Lijx(I,K,J));
	fflush(stdout);
      }
      if(I-K > Lij)
	misscnt += I-K-Lij;
      if(J > Lijx)
	falsecnt += J-Lijx;
      if(G <= -2)
	EndOutlierCnt++;
    }
  }
  if(DEBUG) assert(I==Ilist[0]);
  if(DEBUG) assert(K==Klist[0]);
  if(DEBUG) assert(J==Jlist[0]);
  /* right end */
  //  int Rijx = extend ? min(M,A[I][K][J].Rijx) : (R <= -2 ? J : M);
  int Rij = (R <= -2) ? I : min(N,A->Rij(I,K,J));
  int Rijx = (R <= -2) ? J : extend ? min(M,A->Rijx(I,K,J)) : M;
  if(VERB>=2 && verb){
    printf("I=%d,K=%d,J=%d:R=%d,Rij=%d,Rijx=%d,A[I][K][J].Rij=%d,A[I][K][J].Rijx=%d\n",
	   I,K,J,R,Rij,Rijx,A->Rij(I,K,J),A->Rijx(I,K,J));
    fflush(stdout);
  }
  if(I < Rij)
    misscnt += Rij-I;
  if(J < Rijx)
    falsecnt += Rijx-J;
  if(DEBUG) assert((Poutlier <= 0.0) ? F== U-1 : F <= U-1);
  if(R <= -2)
    EndOutlierCnt++;
  
  double logOutlierMN = 0.0;
  if(OUTLIER_PV && OutlierMN > 1e+300){
    G = H = -1; D = 0;
    for(int T=0;T < U;T++, G = I, H = J, D = K){
      I = Ilist[U-1-T];
      K = Klist[U-1-T];
      J = Jlist[U-1-T];
      if(T>0){
	if(outlier[U-T]){
	  int m = J-H;
	  int n = I-K-G;
	  logOutlierMN += log((double)(m*n));
	}
      }
    }
  } else
    logOutlierMN = log(OutlierMN);

  int EndOutlierNL = 1;
  if(OUTLIER_PV && EndOutlierPV && EndOutlierCnt >= 1){
    if(EndOutlierCnt >= 2){
      int L = min(Ilist[U-1] - Ilist[0], Jlist[U-1] - Jlist[0]);
      int T = min(N, M);
      EndOutlierNL = T - L + 1;
    } else
      EndOutlierNL = 2;
  }

  double logFP = pvalue(Xfrag,Yfrag,F, misscnt,falsecnt,OutlierCnt,logOutlierMN,EndOutlierNL,Ylambda[c],res[c]*PixelLen,verb,refid,mapid,orientation);
  if(DEBUG && !(isfinite(logFP))){
    (void)pvalue(Xfrag,Yfrag,F, misscnt,falsecnt,OutlierCnt,logOutlierMN,EndOutlierNL,Ylambda[c],res[c]*PixelLen,1,refid,mapid,orientation);
    fprintf(stderr,"alignFP:I=%d,J=%d:logFP=%e,misscnt=%d(Rij=%d),falsecnt=%0.1f(Rijx=%d),Outliers=%d,MN=%0.1f\n",
	    I,J,logFP,misscnt,Rij,falsecnt,Rijx,OutlierCnt,OutlierMN);
    fflush(stderr);
    assert(isfinite(logFP));
  }

  delete [] Ilist;
  delete [] Klist;
  delete [] Jlist;
  delete [] outlier;
  delete [] Xfrag;
  delete [] Yfrag;
  return logFP;
}

/* NOTE : outliers can cause left ends of alignments to be color-misaligned by over 100kb, even though right ends might have small color-misalignment, since only right end color-misalignment is scored */

/** 2-color version of refalignXY for resSD > 0 with 3-D recurrance based alignment */
static void refalignXYsd(FLOAT **YY, int *NN,
			 FLOAT **XX, int *MM,
			 int orientation, 
			 int scaleID,
			 int **KKmax,
			 AL2 *AA,
#if FLAT_MEMORY
			 TFLOAT **Fmem, /**< preallocated memory Fmem[c=0..1][0 ... pcnt[c]*maxM[c]*DPsizF - 1] to be used by A[c].field[I][K][] of type TFLOAT */
			 int **Imem, /**< preallocated memory Imem[c=0..1][0 ... pcnt[c]*maxM[c]*DPsizI - 1] to be used by A->field[I][K][] of type int */
			 long long *StrideMax,/* If (MIN_MEM && DIAGONAL>=2 && hashdelta) : the current allocated number of elements per array */
			 int **Jmin, int **Jmax, /**< preallocated memory Jmin[c=0..1][1..maxN[c]], Jmax[c=0..1][1..maxN[c]] used with phash offset to create diagonalized arrays */
			 int **Imin, int **Imax,  /**< preallocated memory Imin[c][1..maxM[c]], Imax[c][1..maxM[c]] used with phash offset to create diagonalized arrays */
			 int tid,
			 long long *limN,
#endif
			 Calign *align,
			 Cmap *rmap, Cmap *nanomap,
			 int refid, int mapid,/**< NOTE : due to multithreading, Gmap[] or refmap[] may be outdated cache values, but nanomap and rmap will correctly represent Gmap[mapid] & refmap[refid] */
			 CXPen *XPen,
			 CHashMatch *phash, /**< If != 0, limit map offset to within hashdelta of phash->offset (in kb) */
			 int *Qmin, int *Qmax)
{
  XPen_init(XX,MM,XPen);

  int isLinear = 0;
  
  int localtype = (NEW==1 && PoutlierEnd > 0.0) ? -3 : -2;
  TFLOAT OutlierEndPen =  OutlierEndBias + OutlierEndPenalty;

  int IImin[MAXCOLOR],IImax[MAXCOLOR];
  for(int c = 0; c < colors; c++){
    IImin[c] = 1;
    IImax[c] = NN[c];
  }
  int Jmin1[2] = {1,1},Jmax1[2] = {MM[0],MM[1]};/* Jmin1[c] and Jmax1[c] are used if DIAGONAL >= 2 and are the allocated J range for index I==IImin[c] :
						   allocated ranges for other I are Jmin1[c]+(I-IImin[c]) ... Jmax1[c]+(I-IImin[c]) */

  if(FAST && Gmap[mapid]->align[0].mapid1 >= 0){
    for(int c = 0; c < colors; c++){
      Calign_brief *palign = &Gmap[mapid]->align[c];
      if(palign->mapid1>=0 && orientation == palign->orientation && scaleID == palign->scaleID && (!BestRef || palign->mapid1 == refid)){
	if(FAST>=2 && Gmap[mapid]->align->score > FAST_SCORE){
	  int N = NN[c];
	  IImin[c] = max(1,palign->sites1_0 - DELTA_Y);
	  IImax[c] = min(N,palign->sites1_M + DELTA_Y);
	}
      } else {
	IImin[c] = 1;
	IImax[c] = 0;
      }
    }
  }

  /* Note: I - Kmax[c][I] must be monotonically non-decreasing functions of I (unless Kmax[I] < 0)*/
  if(DEBUG>=2){
    for(int c = 0; c < colors; c++){
      int *Kmax = KKmax[c];
      int Tmax = 1;
      for(int I = IImin[c]; I <= IImax[c]; I++){
	if(Kmax[I] < 0)
	  continue;
	assert(I-Kmax[I] >= 1);
	if(!(I-Kmax[I] >= Tmax)){
	  printf("refalignXYsd:mapid=%d,orientation=%d\n",mapid,orientation);
	  printf("c=%d,I=%d,Kmax[c][I]=%d,Tmax=%d\n",c,I,Kmax[I],Tmax);
	  assert(I-Kmax[I] >= Tmax);
	}
	if(I-Kmax[I] > Tmax)
	  Tmax = I-Kmax[I];
      }
    }
  }

#ifdef VALGRIND
  int ival = -555555555;
#else
  int ival = -5555;
#endif

  if(DEBUG>=2 && DIAGONAL <= 1){/* set all A->*score(I,K,J)  = nan("NaN") */
    for(int c = 0; c < colors; c++){
      for(int I = 1; I <= NN[c]; I++){
	for(int K = 0; K <= KKmax[c][I]; K++){
	  for(int J = 1; J <= MM[c]; J++){
	    AA[c].scoreN(I,K,J) = AA[c].scoreEnd(I,K,J) = AA[c].Lscore(I,K,J) =
	      AA[c].LscoreEnd(I,K,J) = AA[c].Rscore(I,K,J) = AA[c].RscoreEnd(I,K,J) = aNaN;
	    AA[c].REnd(I,K,J) = AA[c].LEnd(I,K,J) = AA[c].R(I,K,J) = AA[c].L(I,K,J) = AA[c].Rijx(I,K,J) = AA[c].Lijx(I,K,J) = AA[c].Rij(I,K,J) =  AA[c].Lij(I,K,J) = 
	      AA[c].GEnd(I,K,J) = AA[c].TEnd(I,K,J) = AA[c].HEnd(I,K,J) = AA[c].H(I,K,J) = AA[c].T(I,K,J) = AA[c].G(I,K,J) = ival;
	  }
	}
      }
    }
  }

  if(phash && hashdelta){/* use hashtable information to limit range IImin[c] .. IImax[c] (If DIAGONAL : also adjust Jmin[c][I]..Jmax[c][I] and reallocate AA[c].field[][][]) */
    double offset = phash->offset;
    /* HERE HERE : need to correct offset if scaleID != 0 */
    if(nanomap->origmap){/* correct phash->offset since it is based on original map of which nanomap is a sub-fragment */
      if(DEBUG) assert(HashMultiMatch);
      Cmap *origmap = nanomap->origmap;
      FLOAT *origX = origmap->site[0];
      if(DEBUG) assert(!origmap->origmap);
      if(!orientation){
	int L = nanomap->left[0];
	offset = phash->offset + origX[L] - XX[0][1];
      } else {
	int R = nanomap->right[0];
	int origM = origmap->numsite[0];
	offset = phash->offset + (origX[origM+1] - origX[R]) - (XX[0][MM[0]+1] - XX[0][MM[0]]);
      }
    }


    for(int c = 0; c < colors; c++){
      FLOAT *Y = YY[c];
      FLOAT *X = XX[c];
      int N = NN[c];
      int M = MM[c];
      int *Kmax = KKmax[c];
      AL2 *A = &AA[c];

      int &IMIN = IImin[c];
      int &IMAX = IImax[c];
      int *JMIN = Jmin[c];
      int *JMAX = Jmax[c];

      double len = min(Y[N+1],X[M+1]);
      double deltaoffset =  max(Ylambda[c]*DELTA_Y,Xlambda[c]*DELTA_X);
      double offsetSD = max(deltaoffset, hashdelta * max(Ylambda[c],Xlambda[c])) + 5.0*sqrt(SF[c]*SF[c] + fabs(SD[c])*SD[c]*len + SR[c]*SR[c]*len*len);
      if(PVERB){
	printf("refalignXYsd:refid=%d(id=%lld),mapid=%d(id=%lld),or=%d,c=%d:phash=%p(offset=%d,offset=%0.1f,hashscore=%d,orientation=%d):offsetSD=%0.3f,Y[N+1]=%0.3f,X[M+1]=%0.3f,Ylambda=%0.3f,Xlambda=%0.3f,SF=%0.6f,SD=%0.6f,SR=%0.6f\n",
	       refid,rmap->id,mapid,nanomap->id,orientation,c,phash,phash->offset,offset,phash->hashscore,phash->orientation, offsetSD, Y[N+1],X[M+1],Ylambda[c],Xlambda[c],SF[c],SD[c],SR[c]);
	fflush(stdout);
      }

      while(IMIN < IMAX && Y[IMIN+1]-X[1] < offset - offsetSD)
	IMIN++;
      while(IMIN < IMAX && Y[IMAX-1]-X[M] > offset + offsetSD)
	IMAX--;

      for(int I = IMIN; I <= IMAX; I++){
	JMIN[I] = 1;
	JMAX[I] = M;
      }
      if(DIAGONAL){
	for(int I = IMIN; I <= IMAX; I++){
	  while(JMIN[I] < JMAX[I] && Y[I] - X[JMIN[I]+1] > offset + offsetSD)
	    JMIN[I]++;
	  while(JMIN[I] < JMAX[I] && Y[I] - X[JMAX[I]-1] < offset - offsetSD)
	    JMAX[I]--;
	}

	/* enforce monotonic restrictions IMIN .. IMAX :
	   1. JMIN[I=IMAX..IMIN] is a non-increasing function as I decreases.
	   2. JMAX[I=IMIN..IMAX] is a non-decreasing function as I increase.
	*/
	for(int I = IMAX; --I >= IMIN;)
	  JMIN[I] = min(JMIN[I],JMIN[I+1]);
	for(int I = IMIN; ++I <= IMAX;)
	  JMAX[I] = max(JMAX[I],JMAX[I-1]);

	if(DIAGONAL>=2){
	  /* Re-allocate A->field(I,K,0)[J] based on JMIN[],JMAX[] for more compact real memory usage */
	  /* expand band so boundaries are exactly diagnonal (slope = 1:1) */
	  Jmin1[c] = 1;
	  Jmax1[c] = 0;
	  /*	  Jmin1[c] = JMIN[IMIN];
		  Jmax1[c] = JMAX[IMIN];*/
	  for(int I = IMIN; I <= IMAX; I++){
	    Jmin1[c] = min(Jmin1[c],JMIN[I]-(I-IMIN));
	    Jmax1[c] = max(Jmax1[c],JMAX[I]-(I-IMIN));
	  }

	  int kmax = Kmax[IMIN];
	  for(int I = IMIN+1; I <= IMAX; I++)
	    kmax = max(kmax,Kmax[I]);

	  int deltaJ = Jmax1[c] - Jmin1[c] + 1;
	  int origdeltaJ = deltaJ;
	  if(deltaJ >= M){/* revert to non-diagonal method base(I=IMIN..IMAX,K=0..kmax,J=1..M) = base[I*Istride + K*Kstride + J - firstindex] */
	    Jmin1[c] = 1;
	    Jmax1[c] = M;
	    deltaJ = M;
	    A->Istride = M;
	  } else
	    A->Istride = deltaJ - 1;

	  if(FENCE > 0){
	    Jmin1[c] -= FENCE;
	    Jmax1[c] += FENCE;
	    deltaJ += 2*FENCE;
	    origdeltaJ += 2*FENCE;
	    A->Istride += 2*FENCE;
	  }

	  A->kmax = kmax;
	  A->Kstride = (IMAX - IMIN + 1LL) * deltaJ;
	  long long stride = (1LL + A->kmax) * A->Kstride;
	  long long firstindex = IMIN * A->Istride + Jmin1[c];

	  long long index1 = IMIN * A->Istride + (IMIN <= IMAX ? JMIN[IMIN] : 1) ;
	  long long index2 = IMAX * A->Istride + kmax * A->Kstride + (IMIN <= IMAX ? JMAX[IMAX] : 0);
	  if((VERB>=2 && origdeltaJ != deltaJ) || (DEBUG && IMIN <= IMAX && !(firstindex <= index1 && index2 < firstindex + stride))){
            #pragma omp critical
	    {
	      printf("tid=%d,rid=%lld,mid=%lld,or=%d,c=%d:deltaJ=%d->%d,M=%d,N=%d,IMIN=%d,IMAX=%d,firstindex=%lld(%lld..%lld),stride=%lld:Jmin1[c]=%d,JMIN[IMIN]=%d,JMAX[IMAX]=%d\n",
		     tid,rmap->id,nanomap->id,orientation,c,origdeltaJ,deltaJ,M,N,IMIN,IMAX, firstindex, index1, index2, stride,Jmin1[c],(IMIN<=IMAX ? JMIN[IMIN] : 1),(IMIN<=IMAX ? JMAX[IMAX] : 0));
	      fflush(stdout);

	      if(DEBUG && IMIN <= IMAX) assert(firstindex <= index1 && index2 < firstindex + stride);
	    }
	  }

	  if(MIN_MEM && (stride > StrideMax[c] || stride*2 < StrideMax[c])){
	    free(Fmem[c]);
	    long long origStrideMax = StrideMax[c];
	    StrideMax[c] = stride*3/2;
	    Fmem[c] = (TFLOAT *)malloc(StrideMax[c] * (DPsizF*sizeof(TFLOAT) + DPsizI*sizeof(int)));
	    if(!Fmem[c]){
	      printf("malloc(%llu) failed : c=%d,stride=%lld,StrideMax=%lld\n",(unsigned long long)(StrideMax[c] * (DPsizF*sizeof(TFLOAT) + DPsizI*sizeof(int))), c, stride, StrideMax[c]);
	      exit(1);
	    }
	    Imem[c] = (int *) &Fmem[c][DPsizF * StrideMax[c]];
	    if(VERB>=2){
              #pragma omp critical
	      {
		printf("tid=%d,rid=%lld,mid=%lld,or=%d,c=%d:Allocated %llu bytes for Fmem[] & Imem[] : StrideMax=%lld->%lld,stride=%lld,Fmem=%p,Imem=%p,firstindex=%lld(%lld..%lld),IMIN=%d,IMAX=%d,Istride=%lld,Jmin1=%d,Jmax1=%d\n",
		       tid, rmap->id,nanomap->id,orientation,c,(unsigned long long)(StrideMax[c] * (DPsizF*sizeof(TFLOAT) + DPsizI*sizeof(int))), origStrideMax, StrideMax[c], stride, Fmem[c], Imem[c], 
		       firstindex, index1, index2, IMIN,IMAX,A->Istride,Jmin1[c],Jmax1[c]);
		fflush(stdout);
	      }
	    }
	  }

	  long long mcnt = 0;

          #define FBLOCK(field) {                             \
	    A->field##_base = &Fmem[c][mcnt - firstindex];    \
	    mcnt += stride;                                   \
	  }

	  FBLOCK(scoreN);
	  FBLOCK(scoreEnd);
	  FBLOCK(Lscore);
	  FBLOCK(LscoreEnd);
	  FBLOCK(Rscore);
	  FBLOCK(RscoreEnd);

          #undef FBLOCK

	  if(DEBUG && !(mcnt <= stride * DPsizF)){
	    #pragma omp critical
	    {
	      printf("refid=%d,mapid=%d,or=%d,c=%d:N=%d,M=%d,maxM=%lld,maxN=%lld:deltaJ=%d,kmax=%d,Istride=%lld,Kstride=%lld,pcnt=%lld,stride=%lld,DPsizF=%lld\n",
		     refid,mapid,orientation,c,N,M,maxM[c],maxN[c],deltaJ,kmax,A->Istride,A->Kstride,pcnt[c],stride,DPsizF);
	      fflush(stdout);
	      assert(mcnt <= stride * DPsizF);
	    }
	  }

	  mcnt = 0;
	  if(VERB>=2){
	    printf("refid=%lld,mapid=%lld,or=%d,c=%d:re-allocating A:IMIN=%d,IMAX=%d,Jmin1=%d,Jmax1=%d,deltaJ=%d,kmax=%d,Istride=%lld,Kstride=%lld,stride=%lld,firstindex=%lld\n",
		   rmap->id,nanomap->id,orientation,c,IMIN,IMAX,Jmin1[c],Jmax1[c],deltaJ,kmax,A->Istride,A->Kstride,stride,firstindex);
	    fflush(stdout);
	  }

          #define BLOCK(field) {				\
	    A->field##_base = &Imem[c][mcnt - firstindex];	\
	    mcnt += stride;					\
	  }
	  BLOCK(G);
	  BLOCK(T);
	  BLOCK(H);
	  BLOCK(GEnd);
	  BLOCK(TEnd);
	  BLOCK(HEnd);
	  BLOCK(Lij);
	  BLOCK(Rij);
	  BLOCK(Lijx);
	  BLOCK(Rijx);
	  BLOCK(L);
	  BLOCK(LEnd);
	  BLOCK(R);
	  BLOCK(REnd);

          #undef BLOCK

	  if(DEBUG) assert(mcnt <= stride * DPsizI);

	  if(DEBUG>=2){/* set all TFLOAT fields to nan("NaN") and all int fields to ival */
	    assert(stride <= StrideMax[c]);
	    for(long long t = DPsizF * StrideMax[c]; --t >= 0;)
	      Fmem[c][t] = aNaN;
	    for(long long t = DPsizI * StrideMax[c]; --t >= 0;)
	      Imem[c][t] = ival;
	    //	  memset_int(Imem[c], ival, stride * DPsizI);
	  }
	}
      }

      if(PVERB>=2){
        #pragma omp critical
	{
	  size_t cnt = 0;
	  if(DEBUG) assert(1 <= IMIN && IMAX <= N);
	  for(int I = IMIN; I <= IMAX; I++)
	    cnt += JMAX[I]-JMIN[I]+1;
	  printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d,c=%d:offset=%0.3f,offsetSD=%0.3f,IMIN=%d,IMAX=%d,N=%d,M=%d,cnt=%llu/%llu/%llu,Kstride=%lld,kmax=%d\n",
		 refid,rmap->id,mapid,nanomap->id,orientation,c,offset,offsetSD,IMIN,IMAX,N,M,(unsigned long long)cnt,
		 (unsigned long long)(IMAX-IMIN+1)*M,(unsigned long long)N*M,A->Kstride,A->kmax);
	  if(PVERB>=3)
	    for(int I = IMIN; I <= IMAX; I++)
	      printf("I=%d:JMIN[I]= %d, JMAX[I]= %d\n",I,JMIN[I],JMAX[I]);
	  fflush(stdout);
	}
      }
    }// c = 0..colors-1
  } else { // !(phash && hashdelta)
    for(int c = 0; c < colors; c++){
      int M = MM[c];
      int *JMIN = Jmin[c];
      int *JMAX = Jmax[c];
      for(int I = IImin[c]; I <= IImax[c]; I++){
	JMIN[I] = 1;
	JMAX[I] = M;
      }
    }
  }

  if(DEBUG>=2)
    for(int c = 0; c < colors; c++)
      for(int I = IImin[c]; I <= IImax[c]; I++){
	if(!(1 <= Jmin[c][I] && Jmin[c][I] <= Jmax[c][I]+1 && Jmax[c][I] <= MM[c])){
	  #pragma omp critical
	  {
	    printf("refid=%d,mapid=%d,or=%d,c=%d:IMIN=%d,IMAX=%d:I=%d,JMIN[I]=%d,JMAX[I]=%d,N=%d,M=%d\n",
		   refid,mapid,orientation,c,IImin[c],IImax[c],I,Jmin[c][I],Jmax[c][I],NN[c],MM[c]);
	    fflush(stdout);
	    assert(1 <= Jmin[c][I] && Jmin[c][I] <= Jmax[c][I] && Jmax[c][I] <= MM[c]);
	  }
	}
      }

  /* compute Imin[],Imax[] from IImin,IImax,JMIN[],JMAX[] : same code as at end of setlimit() */
  for(int c = 0; c < colors; c++){
    int M = MM[c];
    int N = NN[c];
    for(int J = 1; J <= M; J++){
      Imin[c][J] = N+1;
      Imax[c][J] = 0;
    }
    for(int I = IImin[c]; I <= IImax[c]; I++){
      for(int J = Jmin[c][I]; J <= Jmax[c][I]; J++){
	Imin[c][J] = min(Imin[c][J],I);
	Imax[c][J] = max(Imax[c][J],I);
      }
    }

    if(DEBUG>=2)
      for(int J = 1; J <= M; J++)
	for(int I = Imin[c][J]; I <= Imax[c][J]; I++)
	  assert(Jmin[c][I] <= J && J <= Jmax[c][I]);

    if(DEBUG>=2 && !phash){
      for(int I = IImin[c]; I <= IImax[c]; I++){
	assert(Jmin[c][I] == 1);
	assert(Jmax[c][I] == M);
      }    
      for(int J = 1; J <= M; J++){
	assert(Imin[c][J] == IImin[c]);
	assert(Imax[c][J] == IImax[c]);
      }
    }
    if(PVERB>=2){
      #pragma omp critical
      {
	printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d,c=%d:IMIN=%d,IMAX=%d,N=%d,M=%d\n",
	       refid,rmap->id,mapid,nanomap->id,orientation,c,IImin[c],IImax[c],N,M);
	fflush(stdout);
      }
    }
  }// c = 0..colors-1

  /* first align each color seperately (without performing back trace) */
  for(int c = 0; c < colors; c++){
    FLOAT *Y = YY[c];
    FLOAT *X = XX[c];
    int N = NN[c];
    int M = MM[c];
    int *Kmax = KKmax[c];
    AL2 *A = &AA[c];
    int IMIN = IImin[c];
    int IMAX = IImax[c];
    int *JMIN = Jmin[c];
    int *JMAX = Jmax[c];
    Cprtab **PRtabY = PRtab[c][refid];
    CYPen ****YPenR = YPen[c][refid];    

    /* initialize A->scoreEnd(I=1..N, K=0..Kmax[I], J=1..M) (== LscoreEnd) as endoutlier alignment scores (type = -3) */
    /* Also initialize scoreN & Lscore to MINSCORE (type = -1) */
    for(int I = IMIN; I <= IMAX; I++){
      if(DEBUG>=2) assert(Kmax[I] >= 0 && Kmax[I] <= I-1 && Kmax[I] <= KMAX);
      for(int K = 0; K <= Kmax[I]; K++){
	TFLOAT sm = Sm(0,I,K,Y,c);
	TFLOAT score = ChimScore + sm;
	int jmin = JMIN[I], jmax = JMAX[I];

	for(int J= jmin; J <= jmax; J++){
	  A->LscoreEnd(I,K,J) = A->scoreEnd(I,K,J) = score;/* local alignment score does not include bias for excluded end of X */
	  A->LEnd(I,K,J) = A->GEnd(I,K,J) = localtype;/* local alignment (-2) OR outlier (-3) */
	  if(NEW>=2 && PoutlierEnd > 0.0){
	    TFLOAT score = BiasEnd2(X[J],J,c) + sm + OutlierEndPen;
	    if(score > A->scoreEnd(I,K,J)){
	      A->LscoreEnd(I,K,J) = A->scoreEnd(I,K,J) = score;
	      A->LEnd(I,K,J) = A->GEnd(I,K,J) = -3;/* outlier */
	    }
	  }

	  /* initialize Lscore,L,LscoreEnd,LEnd */
	  A->scoreN(I,K,J) = A->Lscore(I,K,J) = MINSCORE;
	  A->G(I,K,J) = A->L(I,K,J) = -1;
	  A->T(I,K,J) = 0;

	  /* initialize Lij,Rij,Lijx,Rijx */
	  A->Lij(I,K,J) = I-K;
	  if(DEBUG>=2) assert(A->Lij(I,K,J) >= 0);
	  A->Rij(I,K,J) = I;
	  if(DEBUG>=2) assert(A->Rij(I,K,J) >= I && A->Rij(I,K,J) <= N+1);
	  A->Lijx(I,K,J) = A->Rijx(I,K,J) = J;
	  if(DEBUG>=2) assert(A->Rijx(I,K,J) >= J && A->Rijx(I,K,J) <= M+1);
	}
      }
    }

    if(extend){
      /* initialize A->Lijx(I=IMIN..IMAX, K=0..Kmax[I], J= JMIN[I]..JMAX[I]) as leftmost site in X that overlaps Y with Y[I-K .. I] aligned with X[J] */
      /* only cases where I-K <= DELTA_Y OR J <= DELTA_X are initialized (all other cases can be assumed to use Lijx = 0) */
      int I = IMIN;
      for(; I <= IMAX; I++){
	int Lijx = 0;
	for(int K= 0; K <= Kmax[I]; K++){
	  Lijx = 0;
	  int Jmax = (I-K <= DELTA_Y) ? JMAX[I] : min(DELTA_X, JMAX[I]);
	  FLOAT y = Yc(Y,I,K);
	  for(int J= JMIN[I]; J <= Jmax; J++){ /* update Lijx by incrementing it until X[Lijx] overlaps Y[I,K] */
	    while(X[J]-X[Lijx] > y)
	      Lijx++;
	    if(DEBUG>=2) assert(Lijx <= J && !(X[J]-X[Lijx] > y) && (Lijx==0 || X[J]-X[Lijx-1] > y));
	    A->Lijx(I,K,J) = Lijx;
	  }
	}
	if(Lijx==0 && (!DIAGONAL || JMAX[I] >= (I-Kmax[I] <= DELTA_Y ? M : DELTA_X)))
	  break;
      }
      I++;
      for(; I <= IMAX; I++){
	for(int K= 0; K <= Kmax[I]; K++){
	  int Jmax = (I-K <= DELTA_Y) ? JMAX[I] : min(DELTA_X, JMAX[I]);
	  if(Jmax >= JMIN[I])
	    memset_int(&A->Lijx(I,K,JMIN[I]),0,Jmax-JMIN[I]+1);
	}
      }
    }

    /* initialize A->Lij(I= Imin[J]..Imax[J], K=0..Kmax[I], J=1..M) as leftmost site in Y that overlaps X with Y[I-K .. I] aligned with X[J] */
    /* update A[I=1..N][K=0..Kmax[I]][J=1..M].scoreN (== Lscore) with left ends unaligned and Y[I-K .. I] aligned with X[J] (type = -1) */
    /* only cases where J <= DELTA_X (OR extend && I-K <= DELTA_Y) are initialized */
    if(extend){/* this is the commonly used case */
      int J = 1;
      for(; J <= M; J++){
	int Lij = 0;
	int ImKmax = (J <= DELTA_X) ? Imax[c][J] : min(DELTA_Y, Imax[c][J]);
	int maxLij = 0;
	for(int ImK = max(1, Imin[c][J] - KMAX); ImK <= ImKmax; ImK++){
	  int Lijk = Lij;
	  int I;
	  for(int K = max(0, Imin[c][J] - ImK); (I = ImK + K) <= Imax[c][J] && K <= Kmax[I]; K++){
	    if(DEBUG>=2) assert(I >= Imin[c][J]);
	    FLOAT y = Yc(Y,I,K);
	    TFLOAT score;
	    TFLOAT *AscoreNIKJ = &A->scoreN(I,K,J);
	    TFLOAT *ALscoreIKJ = &A->Lscore(I,K,J);

	    /* update Lijk by incrementing it until Y[Lijk] overlaps X */
	    while(y - X[J] > Y[Lijk])
	      Lijk++;
	    if(DEBUG>=2) assert(Lijk <= I && !(y - X[J] > Y[Lijk]) && (Lijk <= 0 || y - X[J] > Y[Lijk-1]));
	    int Lijky = A->Lij(I,K,J) = min(Lijk, I - K);
	    if(K==0)
	      Lij = Lijk;
	    maxLij = max(maxLij, Lijk);

	    int Lijkx = A->Lijx(I,K,J);

	    TFLOAT Smjik = PRtabY[K][I].Sm /* Sm(J,I,K,Y,c)*/;
	    if((score = Send(min(X[J],y),J+1-max(1,Lijkx),I-K+1-max(1,Lijky),c) + Smjik) > AscoreNIKJ[0]){
	      AscoreNIKJ[0] = ALscoreIKJ[0] = score;
	      A->G(I,K,J) = A->L(I,K,J) = -1;/* no furuther aligned sites, but no local alignment */
	    }
	    if(ENDFIX && Lijkx <= 0){
	      for(int U= Lijky; U < I-K; U++){
		TFLOAT bound = Sbnd(X[J],y-Y[U],J,I-K-U,c) + Smjik;
		if(bound > AscoreNIKJ[0]){
		  AscoreNIKJ[0] = ALscoreIKJ[0] = bound;
		  if(DEBUG>=2) assert(isfinite(A->scoreN(I,K,J)));
		  A->G(I,K,J) = A->L(I,K,J) = -1;/* no further aligned sites, but no local alignment */
		  A->Lij(I,K,J) = U+1;/* NOTE : dont update Lijk, since it is being used by next p's Lij estimate */
		  if(!ENDFIX2)
		    continue;
		}
		if(!ENDFIX2)
		  break;
	      }
	    }
	    if(ENDFIX && Lijk <= 0){
	      for(int U= Lijkx;U < J; U++){
		double bound = Sbnd(((double)(X[J]-X[U])),y,J-U,I-K,c) + Smjik;
		if(bound > AscoreNIKJ[0]){
		  AscoreNIKJ[0] = ALscoreIKJ[0] = bound;
		  A->G(I,K,J) = A->L(I,K,J) = -1;/* no further aligned sites, but no local alignment */
		  Lijkx = U+1;
		  if(!ENDFIX2)
		    continue;
		}
		if(!ENDFIX2)
		  break;
	      }
	      A->Lijx(I,K,J) = Lijkx;
	    }

	    if(ISLINEAR && isLinear){/* try to align left ends of X and Y */
	      double Lend = Sint((TFLOAT)X[J],y,J,I-K,J,I,K,0,PRtabY,Y,c)/* includes Sm(J,I,K,Y,c)*/;
	      if(Lend > AscoreNIKJ[0]){
		AscoreNIKJ[0] = ALscoreIKJ[0] = Lend;
		if(DEBUG>=2) assert(isfinite(A->scoreN(I,K,J)));
		A->G(I,K,J) = A->L(I,K,J) = 0;
	      }
	    }
	  } /* K loop */
	} /* ImK loop */
	if(maxLij <= 0 && (!DIAGONAL || Imax[c][J] >= (J <= DELTA_X ? N : DELTA_Y)))
	  break;/* all subsequent Lij values will be 0 */
      }/* J loop */

      for(J++; J <= M; J++){/* all Lij values will be 0  */
	int ImKmax = (J <= DELTA_X) ? Imax[c][J] : min(DELTA_Y, Imax[c][J]);
	for(int ImK = max(1, Imin[c][J] - KMAX); ImK <= ImKmax; ImK++){
	  int I;
	  for(int K = max(0, Imin[c][J] - ImK); (I = ImK + K) <= Imax[c][J] && K <= Kmax[I]; K++){
	    if(DEBUG>=2) assert(I >= Imin[c][J]);
	    FLOAT y = Yc(Y,I,K);

	    A->Lij(I,K,J) = 0;

	    int Lijkx = A->Lijx(I,K,J);
	    if(DEBUG>=2) assert(Lijkx > 0);
	    TFLOAT Smjik = PRtabY[K][I].Sm /* Sm(J,I,K,Y)*/;
	    TFLOAT score = Send(min(X[J],y),J+1-Lijkx,I-K,c) + Smjik;
	    if(ENDFIX){
	      for(int U= Lijkx; U < J; U++){
		RFLOAT bound = Sbnd(X[J]-X[U],y,J-U,I-K,c) + Smjik;
		if(bound > score){
		  score = bound;
		  Lijkx = U+1;
		  if(!ENDFIX2)
		    continue;
		}
		if(!ENDFIX2)
		  break;
	      }
	      A->Lijx(I,K,J) = Lijkx;
	    }
	    if(score > A->scoreN(I,K,J)){
	      A->scoreN(I,K,J) = A->Lscore(I,K,J) = score;
	      if(DEBUG>=2) assert(isfinite(A->scoreN(I,K,J)));
	      A->G(I,K,J) = A->L(I,K,J) = -1;/* no further aligned sites, but no local alignment */
	    }

	    if(ISLINEAR && isLinear){/* try to align left ends of X and Y */
	      RFLOAT Lend = Sint((RFLOAT)X[J],y,J,I-K,J,I,K,0,PRtabY,Y,c)/* includes Sm(J,I,K,Y)*/;
	      if(Lend > A->scoreN(I,K,J)){
		A->scoreN(I,K,J) = A->Lscore(I,K,J) = Lend;
		if(DEBUG>=2) assert(isfinite(A->scoreN(I,K,J)));
		A->G(I,K,J) = A->L(I,K,J) = 0;
	      }
	    }
	  } /* K loop */
	} /* ImK loop */
      } /* J loop */
    } else {/* extend == 0 : no need to update Lijx */
      int Jmax = min(M, DELTA_X);
      for(int J=1; J <= Jmax; J++){
	int Lij = 0;
	int ImKmax = Imax[c][J];
	for(int ImK = max(1, Imin[c][J] - KMAX); ImK <= ImKmax; ImK ++){
	  int Lijk = Lij;
	  int I;
	  for(int K = max(0,Imin[c][J] - ImK); (I = ImK + K) <= Imax[c][J] && K <= Kmax[I]; K++){
	    if(DEBUG>=2) assert(I >= Imin[c][J]);
	    FLOAT y = Yc(Y,I,K);
	    TFLOAT score;

	    /* update Lijk by incrementing it until Y[Lijk] overlaps X[J] */
	    while(y - Y[Lijk] > X[J])
	      Lijk++;
	    if(DEBUG>=2) assert(Lijk <= I && !(y-Y[Lijk] > X[J]) && (Lijk==0 || y-Y[Lijk-1] > X[J]));
	    A->Lij(I,K,J) = min(Lijk,I-K);

	    if(K==0)
	      Lij = Lijk;

	    int Lijkx = 0;

	    if(Lijk > 0){
	      TFLOAT Smjik = PRtabY[K][I].Sm /* Sm(J,I,K,Y,c)*/;
	      if((score = Send(min(X[J],y),J+1-max(1,Lijkx),I-K+1-max(1,A->Lij(I,K,J)),c) + Smjik) > A->scoreN(I,K,J)){
		A->scoreN(I,K,J) = A->Lscore(I,K,J) = score;
		if(DEBUG>=2) assert(isfinite(A->scoreN(I,K,J)));
		A->G(I,K,J) = A->L(I,K,J) = -1;/* no further aligned sites, but no local alignment */
	      }
	      if(ENDFIX && Lijkx <= 0){
		for(int U= A->Lij(I,K,J);U < I-K; U++){
		  TFLOAT bound = Sbnd(X[J],y-Y[U],J,I-K-U,c) + Smjik;
		  if(bound > A->scoreN(I,K,J)){
		    A->scoreN(I,K,J) = A->Lscore(I,K,J) = bound;
		    if(DEBUG>=2) assert(isfinite(A->scoreN(I,K,J)));
		    A->G(I,K,J) = A->L(I,K,J) = -1;/* no further aligned sites, but no local alignment */
		    A->Lij(I,K,J) = U+1;/* NOTE : dont update Lijk, since it is being used by next p's Lij estimate */
		    if(!ENDFIX2)
		      continue;
		  }
		  if(!ENDFIX2)
		    break;
		}
	      }
	      if(DEBUG>=2) assert(!(ENDFIX && Lijk <= 0));
	    } /* Lijk > 0 */
	    if(ISLINEAR && isLinear){/* try to align left ends of X and Y */
	      double Lend = Sint((TFLOAT)X[J],y,J,I-K,J,I,K,0,PRtabY,Y,c)/* includes Sm(J,I,K,Y,c)*/;
	      if(Lend > A->scoreN(I,K,J)){
		A->scoreN(I,K,J) = A->Lscore(I,K,J) = Lend;
		if(DEBUG>=2) assert(isfinite(A->scoreN(I,K,J)));
		A->G(I,K,J) = A->L(I,K,J) = 0;
	      }
	    }
	  }
	}
      }
    }/* extend == 0 */

    if(DEBUG>=2){
      for(int I = IMIN; I <= IMAX; I++){
	for(int K = 0; K <= Kmax[I]; K++){
	  for(int J= JMIN[I]; J <= JMAX[I]; J++){
	    if(VERB>=2 && mapid==0  && orientation && c==1 && I==55 && K==1 && J==2){
	      printf("mapid=%d,or=%d,c=%d:I=%d,K=%d,J=%d:score=%0.6f(End=%0.6f),Lscore=%0.6f(End=%0.6f),G=%d,GEnd=%d(L=%d,LEnd=%d)\n",
		     mapid,orientation,c,I,K,J,AA[c].scoreN(I,K,J),AA[c].scoreEnd(I,K,J),AA[c].Lscore(I,K,J),AA[c].LscoreEnd(I,K,J),AA[c].G(I,K,J),AA[c].GEnd(I,K,J),AA[c].L(I,K,J),AA[c].LEnd(I,K,J));
	      fflush(stdout);
	    }
	    if(!(AA[c].Lscore(I,K,J) <= AA[c].scoreN(I,K,J))){
	      #pragma omp critical
	      {
		printf("mapid=%d,or=%d,c=%d:I=%d,K=%d,J=%d:scoreN=%0.6f(End=%0.6f),Lscore=%0.6f(End=%0.6f),G=%d,GEnd=%d(L=%d,LEnd=%d)\n",
		       mapid,orientation,c,I,K,J,AA[c].scoreN(I,K,J),AA[c].scoreEnd(I,K,J),AA[c].Lscore(I,K,J),AA[c].LscoreEnd(I,K,J),AA[c].G(I,K,J),AA[c].GEnd(I,K,J),AA[c].L(I,K,J),AA[c].LEnd(I,K,J));
		fflush(stdout);
		assert(AA[c].Lscore(I,K,J) <= AA[c].scoreN(I,K,J));
	      }
	    }
	    assert(AA[c].LscoreEnd(I,K,J) <= AA[c].scoreEnd(I,K,J));
	  }
	}
      }
    }

    /* initialize A->Rij(I= Imin[J]..Imax[J], K=0..Kmax[I], J= M+1-DELTA_X to M) as rightmost site/end in Y that overlaps X when
       Y[I-K..I] is aligned with X[J] */
    /* only cases where J >= M+1-DELTA_X (OR extend && I >= N+1-DELTA_Y) are initialized */
    int Jmin = (!extend && M >= DELTA_X) ? M+1-DELTA_X : 1;
    if(PVERB>=3 && c== C_TRACE){
      #pragma omp critical
      {
	printf("refid=%d,mapid=%d,or=%d:c=%d:Jmin=%d,M=%d\n",  refid,mapid,orientation,c,Jmin,M);
	fflush(stdout);
      }
    }
    
    for(int J= Jmin; J <= M; J++){
      int Rij = N+1;
      int imin = (J >= M+1-DELTA_X) ? Imin[c][J] : max(Imin[c][J], N+1-DELTA_Y);
      if(PVERB>=3 && c== C_TRACE){
        #pragma omp critical
	{
	  printf("refid=%d,mapid=%d,or=%d:c=%d:J=%d,I=%d .. %d,M=%d,N=%d\n",  refid,mapid,orientation,c,J,Imax[c][J],imin,M,N);
	  fflush(stdout);
	}
      }

      for(int I= Imax[c][J]; I >= imin; I--){
	int Rijk = Rij;
	int *ARijI0J = &A->Rij(I,0,J);
	long long Kstride = A->Kstride;

	for(int K= 0; K <= Kmax[I]; K++){
	  double y = Yc(Y,I,K);
	  /* update Rijk by decrementing it until Y[Rijk] overlaps X[J] */
	  while(Y[Rijk] > X[M+1]-X[J]+y)
	    Rijk--;
	  if(PVERB>=3 && c== C_TRACE && I == I_TRACE && K== K_TRACE && J == J_TRACE){
	    #pragma omp critical
	    {
	      printf("refid=%d,mapid=%d,or=%d:c=%d,I=%d,K=%d,J=%d,N=%d,M=%d:Rij=%d,Rijk=%d,A->Rij(I,K,J)->%d\n",
		     refid,mapid,orientation,c,I,K,J,N,M,Rij,Rijk,max(Rijk,I));
	      fflush(stdout);
	    }
	  }
	  if(DEBUG>=2) assert(Rijk >= I-K && !(Y[Rijk] > X[M+1] - X[J] + y) && (Rijk > N || Y[Rijk+1] > X[M+1] - X[J] + y));

	  *ARijI0J = max(Rijk,I);
	  ARijI0J += Kstride;

	  if(K==0)
	    Rij = Rijk;
	}
      }
    }

    if(extend){
      /* initialize A->Rijx(I=1..N, K=0..Kmax[I], J= JMIN[I] to JMAX[I]) as rightmost site in X that overlaps Y when
	 Y[I-K .. I] is aligned with X[J] */
      /* only cases where I >= N+1-DELTA_Y OR J >= M+1-DELTA_X are initialized */
      for(int I= IMIN; I <= IMAX; I++){
	for(int K= 0; K <= Kmax[I]; K++){
	  double y = Y[N+1] - Yc(Y,I,K);
	  int Rijkx = M+1;
	  int Jmin = (I >= N+1 - DELTA_Y) ? JMIN[I] : max(JMIN[I],M+1 - DELTA_X);
	  for(int J= JMAX[I]; J >= Jmin; J--){
	    /* update Rijkx by decrementing it until X[Rijkx] overlaps Y */
	    while(X[Rijkx]-X[J] > y)
	      Rijkx--;
	    if(DEBUG>=2) assert(Rijkx >= J && !(X[Rijkx]-X[J] > y) && (Rijkx > M || X[Rijkx+1]-X[J] > y));
	    A->Rijx(I,K,J) = Rijkx;
	  }
	}
      }
    }

    /* Dynamic programming 3-D recurrance */
    if(VERB >=3 && !orientation)
      printf("refalignXYsd:IMIN=%d,IMAX=%d\n",IMIN,IMAX);

    for(int I=IMIN; I <= IMAX; I++){
      int jmin = JMIN[I];
      int jmax = JMAX[I];
      for(int K= 0; K <= Kmax[I]; K++){
	double Yik = Yc(Y,I,K);

	for(int J= jmin; J <= jmax; J++){
	  if(extendonly && Refine==2)/* skip if map is not within (2*Ylambda + EndLen) kb of extending beyond the end of the reference */
	    if(Yik - X[J] > 2.0 * Ylambda12 + EndLen && Y[N+1]-Yik - (X[M+1]-X[J]) > 2.0 * Ylambda12 + EndLen)
	      continue;

	  double Xj = X[J];
	  TFLOAT score = A->scoreN(I,K,J);
	  TFLOAT scoreEnd = A->scoreEnd(I,K,J);
	  int bestG = A->G(I,K,J);
	  int bestT = A->T(I,K,J);
	  int bestH = A->H(I,K,J);
	  int bestGEnd = A->GEnd(I,K,J);
	  int bestTEnd = A->TEnd(I,K,J);
	  int bestHEnd = A->HEnd(I,K,J);
	  int Gmin = max(IMIN,I-K-DELTA_Y);
	  if(VERB>=3 && giter==1 && mapid==1875 && !orientation && c==0 && ((I==5 && K==0 && J==7))){
	    printf("mapid=%d,orientation=%d,c=%d:I=%d,K=%d,J=%d:score=%0.6f,Gmin=%d\n",mapid,orientation,c,I,K,J,score,Gmin);
	    fflush(stdout);
	  }
	  TFLOAT OutlierPenalty3 = MINSCORE;
	  if(outlierExtend){
	    OutlierPenalty3 = OutlierPenalty2[c];
	    if(FIX_OUTLIER_MIS){
	      if(DEBUG>=2) assert(K <= KMAX && I-K >= 1);
	      OutlierPenalty3 += PRtabY[K][I].Sm;// Sm(0,I,K,Y);
	    }
	  }

	  for(int G= I - K; --G >= Gmin;){
	    int Hmin = max(JMIN[G], J - DELTA_X);
	    int Hmax = min(JMAX[G], J - 1);
	    for(int T = Kmax[G]; T >= 0; T--){
	      if(DEBUG>=2) assert(G-T >= 1 && G <= N);
	      TFLOAT deltaY,Ivar,Sm;
	      TFLOAT penY = SintY(I-K-G,I,K,T,deltaY,Ivar,Sm,YPenR);
	      if(DEBUG>=2) assert(fabs(deltaY - (Yc(YY[c],I,K) - Yc(YY[c],G,T))) < ((sizeof(TFLOAT)==8) ? 1e-10 : 1e-7) * YY[c][I]);
	      TFLOAT OutPenBC = OutlierPenaltyBC[I-K-G];
	      int H = Hmax;
	      for(; H >= Hmin; H--){
		if(DEBUG>=2) assert(H >= max(1,JMIN[G]) && H <= min(J,JMAX[G]));
		int m = J-H;
		TFLOAT deltaX = XPen->deltaX[c][m][J];
		if(DEBUG>=2) assert(fabs(deltaX - (Xj - X[H])) < max((TFLOAT)1.0,deltaX) * ((sizeof(TFLOAT)==8) ? 1e-10 : 1e-7) * Xj);
		TFLOAT prevscore = A->scoreN(G,T,H);
		if(DEBUG>=2){
		  TFLOAT Pen = (J-H) * FpPenalty[c] + deltaX * Frate[c];
		  if(!(fabs(Pen - XPen->Pen[c][m][J]) < 1e-5)){
		    #pragma omp critical
		    {
		      printf("c=%d,I=%d,K=%d,J=%d,G=%d,T=%d,H=%d,m=%d:deltaX=%0.6f,FpPenalty[c]=%0.6f,Frate[c]=%0.8f:Pen=%0.6f,XPen->Pen[c][m][J]=%0.6f\n",
			     c,I,K,J,G,T,H,m,deltaX,FpPenalty[c],Frate[c],Pen,XPen->Pen[c][m][J]);
		      fflush(stdout);
		      assert(fabs(Pen - XPen->Pen[c][m][J]) < 1e-5);
		    }
		  }
		}
		TFLOAT newscore = SintX(deltaX,deltaY,penY,Sm,Ivar,XPen->Pen[c][m][J],XPen->OutlierBias[c][m][J],XPen->PenBias[c][m][J],OutPenBC + OutlierPenaltyBC[m],c);
		TFLOAT nscore;
		if(DEBUG>=2) nscore = Sint(deltaX,deltaY,J-H,I-K-G,J,I,K,T,PRtabY,Y,c);
		if((VERB>=3 && TRACE && rmap->id == REF_TRACE && nanomap->id == MAP_TRACE && orientation == OR_TRACE && c==0 && I==27 && K==0 && J==4) ||
		   (DEBUG>=2 && !(fabs(nscore - newscore) < 1e-5 * max(((TFLOAT)1.0),fabs(nscore))))){
		  #pragma omp critical
		  {
		    TFLOAT Pen,PenSm,Bias;
		    SintDetail(deltaX,deltaY,J-H,I-K-G,J,I,K,T,PRtabY,Y,c,Bias,Pen,PenSm,1);
		    TFLOAT OutPen = OutlierPenalty[c];
		    if(outlierBC)
		      OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-K-G];
		    TFLOAT iscore = Bias + OutlierBias + (FIX_OUTLIER_MIS ? max(Pen,OutPen) + PenSm : max(Pen+PenSm, OutPen));
		    
		    if(DEBUG>=2){
		      printf("I=%d,K=%d,J=%d,G=%d,T=%d,H=%d:deltaX=%0.3f,deltaY=%0.3f,A[G][T][H].score=%0.6f,SintX=%0.6f,Sint=%0.6f(Bias=%0.6f,Pen=%0.6f,Sm=%0.6f,OutPen=%0.6f:iscore=%0.6f),total=%0.6f(best=%0.6f,G=%d,T=%d,H=%d)\n",
			     I,K,J,G,T,H,deltaX,deltaY,prevscore,newscore,nscore,Bias,Pen,PenSm,OutPen,iscore,prevscore+newscore,score,bestG,bestT,bestH);
		      printf("     penY=%0.6f,Ivar=%0.6f,Sm=%0.6f;XPen:Pen=%0.6f,OutlierBias=%0.6f,PenBias=%0.6f,outlier=%d\n",
			     penY,Ivar,Sm,XPen->Pen[c][m][J],XPen->OutlierBias[c][m][J],XPen->PenBias[c][m][J],OUTLIER(maptype,deltaX,deltaY));
		    } else
		      printf("I=%d,K=%d,J=%d,G=%d,T=%d,H=%d:deltaX=%0.3f,deltaY=%0.3f,A[G][T][H].score=%0.6f,SintX=%0.6f(Bias=%0.6f,Pen=%0.6f,Sm=%0.6f,OutPen=%0.6f:iscore=%0.6f),total=%0.6f(best=%0.6f,G=%d,T=%d,H=%d)\n",
			     I,K,J,G,T,H,deltaX,deltaY,prevscore,newscore,Bias,Pen,PenSm,OutPen,iscore,prevscore+newscore,score,bestG,bestT,bestH);

		    fflush(stdout);
		    if(DEBUG>=2) assert(fabs(nscore - newscore) < 1e-5 * max(((TFLOAT)1.0), fabs(nscore)));
		  }
		}
		if((prevscore += newscore) > score){
		  score = prevscore;
		  bestG = G;
		  bestT = T;
		  bestH = H;
		}
		double prevscoreEnd = A->scoreEnd(G,T,H);
		if((prevscoreEnd += newscore) > scoreEnd){
		  scoreEnd = prevscoreEnd;
		  bestGEnd = G;
		  bestTEnd = T;
		  bestHEnd = H;
		}
	      }// H = Hmax .. Hmin
	    } // T = 0 .. Kmax[G]
	  } // G = I-K-1 .. Gmin
	  if(outlierExtend){// HERE HERE : adapt from refalign.cpp
	    printf("-outlierExtend not yet implemented for 2 colors\n");
	    exit(1);
	  }

	  A->scoreN(I,K,J) = score;
	  A->scoreEnd(I,K,J) = scoreEnd;
	  A->G(I,K,J) = bestG;
	  A->T(I,K,J) = bestT;
	  A->H(I,K,J) = bestH;
	  A->GEnd(I,K,J) = bestGEnd;
	  A->TEnd(I,K,J) = bestTEnd;
	  A->HEnd(I,K,J) = bestHEnd;
	  if(VERB>=3 && giter==1 && refid == 0 && mapid==0 && !orientation && c==0 && ((I==219 && K==0 && J==1) /* ||
										       (I==6 && K==0 && J==8) ||
										       (I==7 && K==0 && J==9) ||
										       (I==8 && K==0 && J==10) ||
										       (I==9 && K==0 && J==11) ||
										       (I==10 && K==0 && J==12) ||
										       (I==11 && K==0 && J==13) ||
										       (I==12 && K==0 && J==15) ||
										       (I==14 && K==1 && J==16) ||
										       (I==15 && K==0 && J==17) ||
										       (I==16 && K==0 && J==18)*/ )){
	    printf("mapid=%d,orientation=%d,c=%d:I=%d,K=%d,J=%d:score=%0.6f,G=%d,T=%d,H=%d\n",mapid,orientation,c,I,K,J,score,A->G(I,K,J),A->T(I,K,J),A->H(I,K,J));
	    fflush(stdout);
	  }
	}
      }
    }

    if(DEBUG>=2){
      for(int I = IMIN; I <= IMAX; I++){
	for(int K = 0; K <= Kmax[I]; K++){
	  for(int J= JMIN[I]; J <= JMAX[I]; J++){
	    if(VERB>=2 && mapid==0  && orientation && c==1 && I==55 && K==1 && J==2){
	      printf("After Rec:mapid=%d,or=%d,c=%d:I=%d,K=%d,J=%d:score=%0.6f(End=%0.6f),Lscore=%0.6f(End=%0.6f),G=%d,GEnd=%d(L=%d,LEnd=%d)\n",
		     mapid,orientation,c,I,K,J,A->scoreN(I,K,J),A->scoreEnd(I,K,J),A->Lscore(I,K,J),A->LscoreEnd(I,K,J),A->G(I,K,J),A->GEnd(I,K,J),A->L(I,K,J),A->LEnd(I,K,J));
	      fflush(stdout);
	    }
	    assert(A->Lscore(I,K,J) <= A->scoreN(I,K,J));
	    assert(A->LscoreEnd(I,K,J) <= A->scoreEnd(I,K,J));
	  }
	}
      }
    }

    /* compute right hand alignment end score Rscore,RscoreEnd */
    for(int I= IMIN; I <= IMAX; I++){
      for(int K = 0; K <= Kmax[I]; K++){
	double Yik = Yc(Y,I,K);

	/* Initialize A->Rscore(I,K,J=JMIN[I]..JMAX[I]) & consider right end local alignments */
	for(int J = JMIN[I]; J <= JMAX[I]; J++){
	  A->Rscore(I,K,J) = MINSCORE;
	  A->R(I,K,J) = -1;
	  A->RscoreEnd(I,K,J) = ChimScore;
	  A->REnd(I,K,J) = localtype;
	  if(NEW>=2 && PoutlierEnd > 0.0){
	    double Rscore = BiasEnd(X[M+1]-X[J],M+1-J,c) + OutlierEndPen;
	    if(Rscore > A->RscoreEnd(I,K,J)){
	      A->RscoreEnd(I,K,J) = Rscore;
	      A->REnd(I,K,J) = -3;/* outlier */
	    }
	  }
	}

	double y = Y[N+1] - Yik;
	int Jmin = (extend && I >= N+1-DELTA_Y) ? JMIN[I] : max(JMIN[I],M+1-DELTA_X);
	int Jmax = JMAX[I];
	for(int J= Jmin ; J <= Jmax; J++){
	  if(PVERB>=3 && c== C_TRACE && I == I_TRACE && K== K_TRACE && J == J_TRACE){
	    #pragma omp critical
	    {
	      printf("refid=%d,mapid=%d,or=%d:c=%d,I=%d,K=%d,J=%d,N=%d,M=%d:A->Rij(I,K,J)=%d,A->Rijx(I,K,J)=%d,A->Rscore(I,K,J)=%0.8e,A->R(I,K,J)=%d,A->RscoreEnd(I,K,J)=%0.8e,A->REnd(I,K,J)=%d\n",
		     refid,mapid,orientation,c,I,K,J,N,M,A->Rij(I,K,J),A->Rijx(I,K,J),A->Rscore(I,K,J),A->R(I,K,J),A->RscoreEnd(I,K,J),A->REnd(I,K,J));
	      fflush(stdout);
	    }
	  }
	  if((extend || A->Rij(I,K,J) <= N)){
	    int Rijx = (extend ? A->Rijx(I,K,J) : M+1);
	    if(DEBUG>=2) assert(Rijx >= J && Rijx <= M+1);
	    double score = Send(min((double)(X[M+1]-X[J]),y),min(M,Rijx)+1-J,min(N,A->Rij(I,K,J))+1-I,c);
	    int linear = -1;/* end type : -1 = normal, 0 = chr end, -2 = local, -3 = outlier */
	    if(ENDFIX && Rijx > M)
	      for(int U= A->Rij(I,K,J); U > I; U--){
		double bound = Sbnd(X[M+1]-X[J],Y[U]-Yik,M+1-J,U-I,c);
		if(bound > score){
		  if(PVERB>=3 && c== C_TRACE && I == I_TRACE && K== K_TRACE && J == J_TRACE){
		    #pragma omp critical
		    {
		      printf("refid=%d,mapid=%d,or=%d:c=%d,I=%d,K=%d,J=%d,N=%d,M=%d:U=%d,A->Rij(I,K,J)=%d->%d,A->Rijx(I,K,J)=%d:score=%0.8f,bound=%0.8f\n",
			     refid,mapid,orientation,c,I,K,J,N,M,U,A->Rij(I,K,J),U-1,A->Rijx(I,K,J),score,bound);
		      fflush(stdout);
		    }
		  }
		  score = bound;
		  A->Rij(I,K,J) = U-1;
		  if(!ENDFIX2)
		    continue;
		}
		if(!ENDFIX2)
		  break;
	      }
	    if(ENDFIX && A->Rij(I,K,J) > N && N+1-I <= DELTA_YEND){/* only possible if extend==1 */
	      if(DEBUG>=2) assert(extend);
	      for(int U= Rijx;U > J;U--){
		double bound = Sbnd(X[U]-X[J],y,U-J,N+1-I,c);
		if(bound > score){
		  score = bound;
		  Rijx = U-1;
		  if(!ENDFIX2)
		    continue;
		}
		if(!ENDFIX2)
		  break;
	      }
	      if(DEBUG>=2) assert(Rijx >= J && Rijx <= M+1);
	      A->Rijx(I,K,J) = Rijx;
	    }
	    if(isLinear){/* try to align right ends of X and Y */
	      double Rend = Sint(X[M+1]-X[J],y,M+1-J,N+1-I,J,I,K,0,PRtabY,Y,c)-PRtabY[K][I].Sm /* Sm(J,I,K,Y,c)*/;
	      if(Rend > score){
		score = Rend;
		linear = 0;
	      }
	    }
	    if(score > A->Rscore(I,K,J)){
	      if(DEBUG){
		assert(I-K >= 1 && I <= N);
		assert(J >= 1 && J <= M);
	      }
	      A->Rscore(I,K,J) = score;
	      A->R(I,K,J) = linear;
	    }
	  }
	}
      }
    }
  } // c = 0 .. color-1

  /* locate best alignment right ends for both colors simultaneously */
  if(DEBUG>=2) assert(colors==2);

  int bestI[2]= {-1,-1}, bestK[2], bestJ[2], bestR[2], bestL[2] = {1,1};/* index refers to color, bestR = right end alignment type, bestL = left end alignment type */
  double bestYik= 0.0;
  double bestscore = MINSCORE;

  /* compute range of sites Qmin[J]..Qmax[J] for 2nd color of XX that is near site J for 1st color of XX */
  //  int *Qmin = new int[MM[0]+1];
  //  int *Qmax = new int[MM[0]+1];
  double Yrange = (DELTA_Y+1)*Ylambda12;
  double Xrange = (DELTA_X+1)*Xtheta12;
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

  int *Ilist = 0,*Klist = 0,*Jlist = 0;
  double *iscore = 0,*outscore = 0;

  if(VERB >= (PVERB?1:2)){
    printf("mapid=%d,or=%d:IImin[0]=%d,IImax[0]=%d,IImin[1]=%d,IImax[1]=%d\n",
	   mapid,orientation,IImin[0],IImax[0],IImin[1],IImax[1]);
    fflush(stdout);
  }

  int pmin = IImin[1], pmax = max(0,IImin[1]-1);

  for(int I = IImin[0]; I <= IImax[0]; I++){
    double bscore = MINSCORE;
    int bK = 0, bJ = -1, bR = -1, bL = -1,   bP = -1, bT = 0, bQ = -1, bS = -1;

    if(VERB>=3){
      printf("result for mapid=%d:I=%d,IImax=%d,%d,M=%d,%d:YY[0][I]=%0.3f,bestscore=%0.6f,bestI=%d,%d,bestK=%d,%d,bestJ=%d,%d,bestR=%d,%d,bestL=%d,%d\n",
	     mapid,I,IImax[0],IImax[1],MM[0],MM[1],YY[0][I],bestscore,bestI[0],bestI[1],bestK[0],bestK[1],bestJ[0],bestJ[1],bestR[0],bestR[1],bestL[0],bestL[1]);
      fflush(stdout);
    }

    /* update pmin,pmax */
    double Ymin = YY[0][I] - Yrange;
    double Ymax = YY[0][I] + Yrange;
    while(pmin <= IImax[1] && YY[1][pmin] < Ymin)
      pmin++;
    pmax = max(pmin-1,pmax);
    while(pmax < IImax[1] && YY[1][pmax+1] <= Ymax)
      pmax++;

    for(int K = 0; K <= KKmax[0][I]; K++){
      double Yik = Yc(YY[0],I,K);
      AL2 *p = &AA[0];
      int jmax = Jmax[0][I];
      for(int J = Jmin[0][I]; J <= jmax; J++){
	double Xj = XX[0][J];
	for(int P = pmin; P <= pmax; P++){
	  /* constrain P so YY[1][P] is within Yrange of YY[0][I] : HERE the P ranges should be precomputed as a function of I, since they are the same for all mapid */
	  if(DEBUG>=2) assert(YY[1][P] >= Ymin && YY[1][P] <= Ymax);

	  for(int T = 0; T <= KKmax[1][P]; T++){
	    double y1 = Yc(YY[1],P,T);
	    double y = y1 - Yik;
	    double len;
	    if(COLOR_SHIFT >= 2)
	      len = max(y1,Yik);
	    AL2 *q = &AA[1];
	    int qmax = min(Qmax[J], Jmax[1][P]);
	    for(int Q = max(Qmin[J], Jmin[1][P]); Q <= qmax; Q++){
	      if(COLOR_SHIFT == 1)
		len = max(XX[1][Q], Xj);
	      double score = SMA(XX[1][Q]-Xj,y, (COLOR_SHIFT > 0 ? len : 0.0), orientation);
	      double Lscore = p->scoreN(I,K,J) + q->scoreN(P,T,Q) + score;
	      double LscoreEnd = p->scoreEnd(I,K,J) + q->scoreEnd(P,T,Q) + score;
	      double Rscore = p->Rscore(I,K,J) + q->Rscore(P,T,Q);
	      double RscoreEnd = p->RscoreEnd(I,K,J) + q->RscoreEnd(P,T,Q);
	      double nscore = max(Lscore,LscoreEnd) + max(Rscore,RscoreEnd);
	      if((VERB>=3 && nscore > 13.0 && nscore > bscore && nscore > bestscore && refid==0 && mapid==0 && orientation /* && I==258 && K==0 && J==38 && P==115 && T==0 && Q==21*/) 
		 ||(VERB >=2 && nscore > bscore && ((p->Rscore(I,K,J) < p->RscoreEnd(I,K,J) && q->Rscore(P,T,Q) > q->RscoreEnd(P,T,Q))||
						    (p->Rscore(I,K,J) > p->RscoreEnd(I,K,J) && q->Rscore(P,T,Q) < q->RscoreEnd(P,T,Q))||
						    (p->Lscore(I,K,J) < p->LscoreEnd(I,K,J) && q->Lscore(P,T,Q) > q->LscoreEnd(P,T,Q))||
						    (p->Lscore(I,K,J) > p->LscoreEnd(I,K,J) && q->Lscore(P,T,Q) < q->LscoreEnd(P,T,Q))))){
		printf("mapid=%d,or=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:p->score=%0.8e(End=%0.8e),p->Rscore=%0.8e(End=%0.8e),q->score=%0.8e(End=%0.8e),q->Rscore=%0.8e(End=%0.8e),x=%0.3f,%0.3f(%0.3f),y=%0.3f,%0.3f(%0.3f),SMA=%0.6f:Lscore=%0.8e(End=%0.8e),Rscore=%0.8e(End=%0.8e):nscore=%0.6f,bscore=%0.6f,bestscore=%0.6f\n",
		       mapid,orientation,I,K,J,P,T,Q,p->scoreN(I,K,J),p->scoreEnd(I,K,J),p->Rscore(I,K,J),p->RscoreEnd(I,K,J),
		       q->scoreN(P,T,Q),q->scoreEnd(P,T,Q),q->Rscore(P,T,Q),q->RscoreEnd(P,T,Q),Xj,XX[1][Q],XX[1][Q]-Xj,Yik,Yc(YY[1],P,T), y, score,
		       Lscore-score,LscoreEnd-score,Rscore,RscoreEnd,nscore,bscore,bestscore);
		fflush(stdout);
	      }
	      if(nscore > bscore){
		bscore = nscore;
	       	bK = K;
		bJ = J;
		if(DEBUG>=2) assert(1 <= bJ && bJ <= MM[0]);
		if(Rscore > RscoreEnd){
		  bR = p->R(I,K,J);
		  if(DEBUG>=2) assert(isLinear ? (bR >= -1) : (bR == -1));
		  bS = q->R(P,T,Q);
		  if(DEBUG>=2) assert(isLinear ? (bS >= -1) : (bS == -1));
		} else {
		  bR = p->REnd(I,K,J);
		  if(DEBUG>=2) assert(bR <= -2);
		  bS = q->REnd(P,T,Q);
		  if(DEBUG>=2) assert(bS <= -2);
		}
		bL = (Lscore > LscoreEnd) ? -1 : -2;
		bP = P;
		bT = T;
		bQ = Q;
	      }
	    }
	  }
	}
      }
    }

    /* HERE : need to allow for cases where one color has no aligned sites (just use Send() for color without aligned sites + trivial SMA()) : this will only be better if there is no sites at all to align for one color */

    if(VERB>=3 && bscore > bestscore + 1e-8){
      printf("result for mapid=%d,or=%d:I=%d,YY[0][I]=%0.3f:bscore=%0.6f,bK=%d,bJ=%d,bR=%d,bL=%d:bP=%d,bT=%d,bQ=%d,bS=%d(bestscore=%0.6f)\n",
	     mapid,orientation,I,YY[0][I],bscore,bK,bJ,bR,bL,bP,bT,bQ,bS,bestscore);
      fflush(stdout);
    }

    double bYik = (bJ >= 1 ) ? Yc(YY[0],I,bK) : -1e100;
    //    double bestExtend = (bestJ[0] < 0) ? 1e100 : max(XX[0][bestJ[0]] - bestYik,XX[0][MM[0]+1]-XX[0][bestJ[0]] - (YY[0][NN[0]+1]-bestYik));
    //    if(bscore >= bestscore - 1e-8 && bscore > MINSCORE && (bscore > bestscore + 1e-8 || (bestExtend > 0.0 && bestExtend > max(XX[0][bJ] - bYik,XX[0][MM[0]+1]-XX[0][bJ] - (YY[0][NN[0]+1]-bYik))))){
    if(bscore > bestscore + 1e-8){
      if((DEBUG>=2 || VERB>=2) /* && bL <= -1 */ ){
	/* backtrace the alignment to make sure the leftmost (local or global) end type matches for both colors */
	/* HERE : outliers can cause left ends of alignments to be color-misaligned by over 100kb, even though right ends might have small color-misalignment, since only right end color-misalignment is scored
	   Also the left endoutlier may be much larger for one color than the other */
	if(!Ilist){    /* allocate space to store alignment in reverse order */
	  int M = max(MM[0],MM[1]);
	  Ilist = new int[M];
	  Klist = new int[M];
	  Jlist = new int[M];
	  iscore = new double[M+1];
	  outscore = new double[M+1];
	}
	double Score[MAXCOLOR], ScoreEnd[MAXCOLOR], LScore[MAXCOLOR],LScoreEnd[MAXCOLOR];
	int LI[MAXCOLOR], L[MAXCOLOR], LEnd[MAXCOLOR], II[MAXCOLOR],KK[MAXCOLOR],JJ[MAXCOLOR];
	for(int c = 0; c < colors; c++){
	  int N = NN[c];
	  int M = MM[c];
	  AL2 *A = &AA[c];

	  int rI = (c==0) ? I : bP;
	  int K = (c==0) ? bK : bT;
	  int J = (c==0) ? bJ : bQ;
	  int U = 0;
	
	  while(rI > 0){
	    //	    DP *p = &A[rI][K][J];
	    if(DEBUG>=2){
	      assert(rI-K >= 1 && rI <= N);
	      if(!(J >= 1 && J <= M)){
		#pragma omp critical
		{
		  printf("refid=%d,mapid=%d,or=%d:I=%d,c=%d:rI=%d,K=%d,J=%d,M=%d,N=%d,bscore=%0.8e,bestscore=%0.8e,bK=%d,bJ=%d;bP=%d,bT=%d,bQ=%d\n",
			 refid,mapid,orientation,I,c,rI,K,J,M,N,bscore,bestscore,bK,bJ,bP,bT,bQ);
		  fflush(stdout);
		  assert(J >= 1 && J <= M);
		}
	      }
	    }
	    Ilist[U] = rI;
	    Klist[U] = K;
	    Jlist[U] = J;
	    U++;
	    int G,T,H;
	    if(bL >= -1){
	      G = A->G(rI,K,J);
	      T = A->T(rI,K,J);
	      H = A->H(rI,K,J);
	    } else {
	      G = A->GEnd(rI,K,J);
	      T = A->TEnd(rI,K,J);
	      H = A->HEnd(rI,K,J);
	    }
	    rI = G;
	    K = T;
	    J = H;
	  }

	  if(DEBUG>=2) assert(U>0);
	  if(DEBUG>=2) assert(bL >= -1 ? (rI >= -1 && rI <= (isLinear ? -1 : 0)) : (rI <= -2));

	  int i = II[c] = Ilist[U-1];
	  int k = KK[c] = Klist[U-1];
	  int j = JJ[c] = Jlist[U-1];
	  //	  DP *p = &A[II[c]][KK[c]][JJ[c]];
	  Score[c] = A->scoreN(i,k,j);
	  ScoreEnd[c] = A->scoreEnd(i,k,j);
	  LI[c] = rI;
	  LScore[c] = A->Lscore(i,k,j);
	  L[c] = A->L(i,k,j);
	  LScoreEnd[c] = A->LscoreEnd(i,k,j);
	  LEnd[c] = A->LEnd(i,k,j);

	  if(DEBUG>=2) assert(A->scoreN(i,k,j) >= A->Lscore(i,k,j) - 1e-8);
	  if(DEBUG>=2) assert(A->scoreEnd(i,k,j) >= A->LscoreEnd(i,k,j) - 1e-8);
	}

	if(VERB>=2 && bscore > ScoreThreshold && (fabs(YY[0][II[0]] - YY[1][II[1]]) >= max(Ylambda[0],Ylambda[1]) * 10.0 ||
						  !((LI[0] <= -2 && LI[1] <= -2)||(LI[0] >= -1 && LI[1] >= -1)))){
          #pragma omp critical
	  {
	    printf("result for mapid=%d,or=%d:N=%d,%d,M=%d,%d:I=%d:bscore=%0.6f,bK=%d,bJ=%d,bR=%d,bL=%d,bP=%d,bT=%d,bQ=%d,bS=%d,Score=%0.6f,%0.6f,ScoreEnd=%0.6e,%0.6e,LScore=%0.6f,%0.6f(L=%d,%d)\n\t LEndScore=%0.6e,%0.6e(LEnd=%d,%d)(bestscore=%0.6f),II=%d,%d,KK=%d,%d,JJ=%d,%d,YY[II]=%0.3f,%0.3f(err=%0.3f),Ylambda=%0.3f,%0.3f,LI=%d,%d!\n",
		   mapid,orientation,NN[0],NN[1],MM[0],MM[1],I,bscore,bK,bJ,bR,bL,bP,bT,bQ,bS,Score[0],Score[1],ScoreEnd[0],ScoreEnd[1],LScore[0],LScore[1],L[0],L[1],LScoreEnd[0],LScoreEnd[1],LEnd[0],LEnd[1],bestscore,
		   II[0],II[1],KK[0],KK[1],JJ[0],JJ[1],YY[0][II[0]],YY[1][II[1]],fabs(YY[0][II[0]] - YY[1][II[1]]), Ylambda[0],Ylambda[1],LI[0],LI[1]);
	    fflush(stdout);
	  }
	}
	if(DEBUG>=2) assert((LI[0] <= -2 && LI[1] <= -2) || (LI[0] >= -1 && LI[1] >= -1));/* left end type must match even if the locations do not match */
      }
      if(bscore > bestscore + 1e-8){
	bestscore = bscore;
	bestI[0] = I;
	bestK[0] = bK;
	bestJ[0] = bJ;
	bestR[0] = bR;
	bestL[0] = bL;
	bestYik = bYik;
	bestI[1] = bP;
	bestK[1] = bT;
	bestJ[1] = bQ;
	bestR[1] = bS;
	bestL[1] = bL;
	if(DEBUG){
	  for(int c = 0; c < colors; c++){
	    assert(bestI[c] >= IImin[c] && bestI[c]-bestK[c] >= 1 && bestI[c] <= IImax[c]);
	    assert(bestJ[c] >= 1 && bestJ[c] <= MM[c]);
	    assert(bestR[c] <= -2 || bestJ[c] >= MM[c]+1 - DELTA_X || (extend && bestI[c] >= NN[c]+1-DELTA_Y));
	  }
	}
      }
    }
  }

  if(bestI[0] < 0 && bestI[1] < 0){/* no valid alignment found */
    if(Ilist){
      delete [] Ilist;
      delete [] Klist;
      delete [] Jlist;
      delete [] iscore;
      delete [] outscore;
    }
    return;
  }

  double MinSpacing = XX[0][MM[0]+1] * SecondSpacing;
  int bestI2[2]= {-1,-1},bestK2[2],bestJ2[2], bestR2[2], bestL2[2]= {1,1};/* index refers to color, bestR2 = right end alignment type, bestL2 = left end alignment type */
  double bestYik2 = 0.0;
  double bestscore2 = MINSCORE;

  if(SecondBest){/* locate 2nd best alignment that is at least X[M+1] * SecondSpacing from best alignment (right end location) */
    int pmin = IImin[1], pmax = max(0,IImin[1]-1);

    for(int I= IImin[0]; I <= IImax[0]; I++){
      double bscore = MINSCORE;
      int bK = 0, bJ = -1, bR = -1, bL = -1, bP = -1, bT = 0, bQ = -1, bS = -1;

      /* update pmin,pmax */
      double Ymin = YY[0][I] - Yrange;
      double Ymax = YY[0][I] + Yrange;
      while(pmin <= IImax[1] && YY[1][pmin] < Ymin)
	pmin++;
      pmax = max(pmin-1,pmax);
      while(pmax < IImax[1] && YY[1][pmax+1] <= Ymax)
	pmax++;

      for(int K = 0; K <= KKmax[0][I]; K++){
	double Yik = Yc(YY[0],I,K);
	if(fabs(Yik - bestYik) < MinSpacing)
	  continue;
	AL2 *p = &AA[0];
	int jmax = Jmax[0][I];
	for(int J= Jmin[0][I] ; J <= jmax; J++){
	  double Xj = XX[0][J];

	  for(int P = pmin; P <= pmax; P++){
	    /* constrain P so YY[1][P] is within Yrange == (DELTA_Y+1)*Ylambda12 of  YY[0][I] */
	    if(DEBUG>=2) assert(YY[1][P] >= Ymin && YY[1][P] <= Ymax);
	    /*	    if(YY[1][P] < YY[0][I] - (DELTA_Y+1)*Ylambda12)
	      continue;
	    if(YY[1][P] > YY[0][I] + (DELTA_Y+1)*Ylambda12)
	    break;*/
	    
	    for(int T = 0; T <= KKmax[1][P]; T++){
	      double y1 = Yc(YY[1],P,T);
	      double y = y1 - Yik;
	      double len;
	      if(COLOR_SHIFT==2)
		len = max(y1,Yik);
	      AL2 *q = &AA[1];
	      int qmax = min(Qmax[J], Jmax[1][P]);
	      for(int Q = max(Qmin[J], Jmin[1][P]); Q <= qmax; Q++){
		if(COLOR_SHIFT==1)
		  len = max(XX[1][Q],Xj);
		double score = SMA(XX[1][Q]-Xj, y, (COLOR_SHIFT > 0 ? len : 0.0), orientation);
		double Lscore = p->scoreN(I,K,J) + q->scoreN(P,T,Q) + score;
		double LscoreEnd = p->scoreEnd(I,K,J) + q->scoreEnd(P,T,Q) + score;
		double Rscore = p->Rscore(I,K,J) + q->Rscore(P,T,Q);
		double RscoreEnd = p->RscoreEnd(I,K,J) + q->RscoreEnd(P,T,Q);
		double nscore = max(Lscore,LscoreEnd) + max(Rscore,RscoreEnd);
		if(nscore > bscore + 1e-8){
		  bK = K;
		  bJ = J;
		  if(Rscore > RscoreEnd){
		    bR = p->R(I,K,J);
		    bS = q->R(P,T,Q);
		    if(DEBUG/* HERE >=2 */) assert(isLinear ? (bR >= -1) : (bR == -1));
		    if(DEBUG/* HERE >=2 */) assert(isLinear ? (bS >= -1) : (bS == -1));
		  } else {
		    bR = p->REnd(I,K,J);
		    bS = q->REnd(P,T,Q);
		    if(DEBUG/* HERE >=2 */) assert(bR <= -2);
		    if(DEBUG/* HERE >=2 */) assert(bS <= -2);
		  }
		  bL = (Lscore > LscoreEnd) ? -1 : -2;
		  bP = P;
		  bT = T;
		  bQ = Q;
		  if(VERB>=2 && refid == 0 && mapid == 0 && orientation == 0 && nscore > bestscore2 + 1e-8){
		    #pragma omp critical
		    {
		      printf("refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:score=%0.6f,scoreN=%0.6f,%0.6f,scoreEnd=%0.6e,%0.6e,Rscore=%0.6f,%0.6f,RscoreEnd=%0.6e,%0.6e\n",
			     refid,mapid,orientation,I,K,J,P,T,Q,score,p->scoreN(I,K,J),q->scoreN(P,T,Q),p->scoreEnd(I,K,J),q->scoreEnd(P,T,Q),p->Rscore(I,K,J),q->Rscore(P,T,Q),p->RscoreEnd(I,K,J),q->RscoreEnd(P,T,Q));
		      printf("    bestscore2=%0.6f -> %0.6f:bestI2=%d,%d,bestK2=%d,%d,bestJ2=%d,%d,bestR2=%d,%d,bestL2=%d,%d,bestYik2=%0.4f(bestYik=%0.4f)\n",
			     bestscore2,nscore,I,bP,bK,bT,bJ,bQ,bR,bS,bL,bL,(bJ >= 1) ? Yc(YY[0],I,bK) : -1e5,bestYik);
		      fflush(stdout);
		    }
		  }
		  bscore = nscore;
		}
	      }
	    }
	  }
	}
      }

      double bYik = (bJ >= 1 ) ? Yc(YY[0],I,bK) : -1e100;
      //      double bestExtend = (bestJ2[0] < 0) ? 1e100 : max(XX[0][bestJ2[0]] - bestYik,XX[0][MM[0]+1]-XX[0][bestJ2[0]] - (YY[0][NN[0]+1]-bestYik2));
      //      if(bscore >= bestscore2 - 1e-8 && bscore > MINSCORE && (bscore > bestscore2 + 1e-8 || (bestExtend > 0.0 && bestExtend > max(XX[0][bJ] - bYik,XX[0][MM[0]+1]-XX[0][bJ] - (YY[0][NN[0]+1]-bYik))))){
      if(bscore > bestscore2 + 1e-8){
	if(VERB>=2 && refid == 0 && mapid == 0 && orientation == 0){
	  printf("refid=%d,mapid=%d,or=%d:bestscore2=%0.6f -> %0.6f\n",refid,mapid,orientation,bestscore2,bscore);
	  fflush(stdout);
	}
	bestscore2 = bscore;
	bestI2[0] = I;
	bestK2[0] = bK;
	bestJ2[0] = bJ;
	bestR2[0] = bR;
	bestL2[0] = bL;
	bestYik2 = bYik;
	bestI2[1] = bP;
	bestK2[1] = bT;
	bestJ2[1] = bQ;
	bestR2[1] = bS;
	bestL2[1] = bL;
	if(DEBUG/* HERE >=2 */){
	  if(!((bestR2[0] >= -1) ? (bestR2[1] >= -1) : (bestR2[1] < -1))){
	    #pragma omp critical
	    {
	      printf("refid=%d,mapid=%d,or=%d:bestR2[0] -> %d, bestR2[1] -> %d, I=%d,bK=%d,bJ=%d,bL=%d,bestscore2 -> %0.6f\n",refid,mapid,orientation,bestR2[0],bestR2[1],I,bK,bJ,bL,bestscore2);
	      fflush(stdout);
	      assert((bestR2[0] >= -1) ? (bestR2[1] >= -1) : (bestR2[1] < -1));
	    }
	  }
	  assert((bestR2[0] >= -1) ? (bestR2[0] == AA[0].R(bestI2[0],bestK2[0],bestJ2[0])) : (bestR2[0] == AA[0].REnd(bestI2[0],bestK2[0],bestJ2[0])));
	  assert((bestR2[1] >= -1) ? (bestR2[1] == AA[1].R(bestI2[1],bestK2[1],bestJ2[1])) : (bestR2[1] == AA[1].REnd(bestI2[1],bestK2[1],bestJ2[1])));
	  assert(bestL2[0] == bestL2[1]);
	  for(int c = 0; c < colors; c++){
	    assert(bestI2[c] >= IImin[c] && bestI2[c]-bestK2[c] >= 1 && bestI2[c] <= IImax[c]);
	    assert(bestJ2[c] >= 1 && bestJ2[c] <= MM[c]);
	    if(!(bestR2[c] <= -2 || bestJ2[c] >= MM[c]+1 - DELTA_X || (extend && bestI2[c] >= NN[c]+1-DELTA_Y))){
	      #pragma omp critical
	      {
		printf("refid=%d,mapid=%d,or=%d,c=%d:bestR2[c]=%d,bestI2[c]=%d,bestK2[c]=%d,bestJ2[c]=%d,MM[c]=%d,NN[c]=%d,DELTA_X=%d,DELTA_Y=%d\n",
		       refid,mapid,orientation,c,bestR2[c],bestI2[c],bestK2[c], bestJ2[c], MM[c],NN[c],DELTA_X,DELTA_Y);
		fflush(stdout);
		assert(bestR2[c] <= -2 || bestJ2[c] >= MM[c]+1 - DELTA_X || (extend && bestI2[c] >= NN[c]+1-DELTA_Y));
	      }
	    }
	  }
	}
      }
    }
  }

  if(bestscore > align[-1].score + 1e-8){ /* backtrack through array for each color to retrieve best alignment starting at right end AA[c][bestI[c]][bestK[c]][bestJ[c]] */
    if(SecondBest && align[-1].score > MINSCORE){  /* check if we need to save previous best alignment as 2nd best alignment */
      int numpairs1 = align->numpairs;
      double alignYik = (numpairs1 <= 0) ? -1e100 : Yc(YY[0],align->sites1[numpairs1-1],align->sitesK1[numpairs1-1]);
      if(DEBUG) assert(bestYik >= 0.0);
      Calign *align2;
      if(align[-1].orientation != orientation || fabs(alignYik - bestYik) >= MinSpacing){
	delete [] align[-1].align2; align[-1].align2 = new Calign[3];
	align2 = &align[-1].align2[1]; 
	align2[-1].mapid1 = refid;
	align2[-1].mapid2 = mapid;

	align2[-1].score = align[-1].score;
	align2[-1].orientation = align[-1].orientation;
	align2[-1].scaleID = align[-1].scaleID;
	align2[-1].logPV = align[-1].logPV;
	for(int c = 0; c < colors; c++){
	  copy(&align2[c],&align[c], 1, 1);

	  //	  align2[c].score = align[c].score;
	  //	  align2[c].orientation = align[c].orientation;
	  //	  align2[c].logPV = align[c].logPV;
	  //	  align2[c].Lend = align[c].Lend;
	  //	  align2[c].Rend = align[c].Rend;

	  //	  align2[c].Lij1 = align[c].Lij1;
	  //	  align2[c].Lij2 = align[c].Lij2;
	  //	  align2[c].Rij1 = align[c].Rij1;
	  //	  align2[c].Rij2 = align[c].Rij2;
	}
	align2[-1].numpairs = align[0].numpairs + align[1].numpairs;
      } else if((align2 = align[-1].align2) && align2->orientation == orientation){/* check if existing 2nd best alignment is now invalid (too close to new best alignment) */
	int numpairs1 = align2[1].numpairs;
	double alignYik2 = (numpairs1 <= 0) ? -1e100 : Yc(YY[0],align2[1].sites1[numpairs1 - 1],align2[1].sitesK1[numpairs1 - 1]);
	if(fabs(alignYik2 - bestYik) < MinSpacing){
	  delete [] align[-1].align2; align[-1].align2 = NULL;
	}
      }
    }

    if(DEBUG){
      //      DP *p[2];
      //      for(int c = 0; c < colors; c++)
      //	p[c] = &AA[c][bestI[c]][bestK[c]][bestJ[c]];
      double y0 = Yc(YY[0],bestI[0],bestK[0]);
      double y1 = Yc(YY[1],bestI[1],bestK[1]);
      double x0 = XX[0][bestJ[0]];
      double x1 = XX[1][bestJ[1]];
      double score = SMA(x1-x0, y1-y0, (COLOR_SHIFT==2) ? max(y0,y1) : max(x0,x1), orientation);
      double Lscore = AA[0].scoreN(bestI[0],bestK[0],bestJ[0]) + AA[1].scoreN(bestI[1],bestK[1],bestJ[1]) + score;
      double LscoreEnd = AA[0].scoreEnd(bestI[0],bestK[0],bestJ[0]) + AA[1].scoreEnd(bestI[1],bestK[1],bestJ[1]) + score;
      double Rscore = AA[0].Rscore(bestI[0],bestK[0],bestJ[0]) + AA[1].Rscore(bestI[1],bestK[1],bestJ[1]);
      double RscoreEnd = AA[0].RscoreEnd(bestI[0],bestK[0],bestJ[0]) + AA[1].RscoreEnd(bestI[1],bestK[1],bestJ[1]);
      if(VERB>=2 && ((AA[0].Rscore(bestI[0],bestK[0],bestJ[0]) > AA[0].RscoreEnd(bestI[0],bestK[0],bestJ[0]) && AA[1].Rscore(bestI[1],bestK[1],bestJ[1]) < AA[1].RscoreEnd(bestI[1],bestK[1],bestJ[1])) || 
		     (AA[0].Rscore(bestI[0],bestK[0],bestJ[0]) < AA[0].RscoreEnd(bestI[0],bestK[0],bestJ[0]) && AA[1].Rscore(bestI[1],bestK[1],bestJ[1]) > AA[1].RscoreEnd(bestI[1],bestK[1],bestJ[1])) ||
		     (AA[0].Lscore(bestI[0],bestK[0],bestJ[0]) > AA[0].LscoreEnd(bestI[0],bestK[0],bestJ[0]) && AA[1].Lscore(bestI[1],bestK[1],bestJ[1]) < AA[1].LscoreEnd(bestI[1],bestK[1],bestJ[1])) ||
		     (AA[0].Lscore(bestI[0],bestK[0],bestJ[0]) < AA[0].LscoreEnd(bestI[0],bestK[0],bestJ[0]) && AA[1].Lscore(bestI[1],bestK[1],bestJ[1]) > AA[1].LscoreEnd(bestI[1],bestK[1],bestJ[1])))){
        #pragma omp critical
	{
	  printf("mapid=%d,or=%d:bestscore=%0.6f:SMA=%0.6f,Lscore=%0.6f(End=%0.6f),Rscore=%0.6f(End=%0.6f),best I=%d,%d,K=%d,%d,J=%d,%d,R=%d,%d,L=%d,N=%d,%d,M=%d,%d,p[0]:Rscore=%0.6f(End=%0.6f),score=%0.6f(End=%0.6f),p[1]:Rscore=%0.6f(End=%0.6f),score=%0.6f(End=%0.6f)\n",
		 mapid,orientation,bestscore,score,Lscore-score,LscoreEnd-score,Rscore,RscoreEnd,bestI[0],bestI[1],bestK[0],bestK[1],bestJ[0],bestJ[1],bestR[0],bestR[1],bestL[0],NN[0],NN[1],MM[0],MM[1],
		 AA[0].Rscore(bestI[0],bestK[0],bestJ[0]),AA[0].RscoreEnd(bestI[0],bestK[0],bestJ[0]),AA[0].scoreN(bestI[0],bestK[0],bestJ[0]),AA[0].scoreEnd(bestI[0],bestK[0],bestJ[0]),
		 AA[1].Rscore(bestI[1],bestK[1],bestJ[1]),AA[1].RscoreEnd(bestI[1],bestK[1],bestJ[1]),AA[1].scoreN(bestI[1],bestK[1],bestJ[1]),AA[1].scoreEnd(bestI[1],bestK[1],bestJ[1]));
	  fflush(stdout);
	}
      }

      assert((bestR[0] >= -1) ? (bestR[1] >= -1) : (bestR[1] < -1));
      assert((bestR[0] >= -1) ? (bestR[0] == AA[0].R(bestI[0],bestK[0],bestJ[0])) : (bestR[0] == AA[0].REnd(bestI[0],bestK[0],bestJ[0])));
      assert((bestR[1] >= -1) ? (bestR[1] == AA[1].R(bestI[1],bestK[1],bestJ[1])) : (bestR[1] == AA[1].REnd(bestI[1],bestK[1],bestJ[1])));
      assert(bestL[0] == bestL[1]);

      double bscore = ((bestR[0] >= -1) ? Rscore : RscoreEnd) + ((bestL[0] >= -1) ? Lscore : LscoreEnd);
      assert(fabs(bscore - bestscore) < 1e-8);
    }

    if(!Ilist){    /* allocate space to store alignment in reverse order */
      int M = max(MM[0],MM[1]);
      Ilist = new int[M];
      Klist = new int[M];
      Jlist = new int[M];
      iscore = new double[M+1];
      outscore = new double[M+1];
    }

    for(int c = 0; c < colors; c++){
      FLOAT *Y = YY[c];
      FLOAT *X = XX[c];
      int N = NN[c];
      int M = MM[c];
      int *Kmax = KKmax[c];
      AL2 *A = &AA[c];
      Cprtab **PRtabY = PRtab[c][refid];

      int I = bestI[c];
      int K = bestK[c];
      int J = bestJ[c];
      if(DEBUG>=2) assert(J >= 1 && J <= M);

      int U = 0;

      /* right end segment */
      // NOTE : if PoutlierEnd > 0.0, bestR[c] <= -2 may be caused either by local alignment or PoutlierEnd
      outscore[0] = iscore[0] = (bestR[c] >= -1 ? A->Rscore(I,K,J) : A->RscoreEnd(I,K,J));

      while(I > 0){
	if(DEBUG){
	  assert(I-K >= 1 && I <= N);
	  if(!(J >= 1 && J <= M)){
	    #pragma omp critical
	    {
	      printf("refid=%d,mapid=%d,or=%d,c=%d:N=%d,M=%d:U=%d,I=%d,K=%d,J=%d\n",refid,mapid,orientation,c,N,M,U,I,K,J);
	      fflush(stdout);
	      assert(J >= 1 && J <= M);
	    }
	  }
	}
	Ilist[U] = I;
	Klist[U] = K;
	Jlist[U] = J;
	if(bestL[c] >= -1){
	  int G = A->G(I,K,J);
	  int T = A->T(I,K,J);
	  int H = A->H(I,K,J);
	  if(G > 0){
	    if(DEBUG>=2) assert(H >= 1 && H <= M);
	    if(DEBUG>=2) assert(T >= 0);
	    iscore[U+1] = A->scoreN(I,K,J) - A->scoreN(G, T, H);
	    double x = X[J] - X[H];
	    double y = Yc(Y,I,K) - Yc(Y,G,T);
	    if(Poutlier > 0.0 && OUTLIER(maptype,x,y)){
	      double Bias,Pen,PenSm;
	      SintDetail(x,y, J-H,I-K-G,J,I,K,T,PRtabY,Y,c,Bias,Pen,PenSm,0);
	      outscore[U+1] = Bias + Pen + PenSm;
	    } else
	      outscore[U+1] = iscore[U+1];
	  } else {/* left end segment */
	    // NOTE : if PoutlierEnd > 0.0, p->G <= -2 may be caused either by local alignment or PoutlierEnd
	    outscore[U+1] = iscore[U+1] = A->scoreN(I,K,J);
	  }
	  I = G;
	  K = T;
	  J = H;
	  if(DEBUG>=2 && I > 0) assert(J >= 1 && J <= M);
	} else {/* bestL[c] <= -2 */
	  int GEnd = A->GEnd(I,K,J);
	  int TEnd = A->TEnd(I,K,J);
	  int HEnd = A->HEnd(I,K,J);
	  if(GEnd > 0){
	    if(DEBUG>=2) assert(HEnd >= 1 && HEnd <= M);
	    if(DEBUG) assert(TEnd >= 0 && HEnd > 0);
	    iscore[U+1] = A->scoreEnd(I,K,J) - A->scoreEnd(GEnd,TEnd,HEnd);
	    double x = X[J]-X[HEnd];
	    double y = Yc(Y,I,K) - Yc(Y,GEnd,TEnd);
	    if(Poutlier > 0.0 && OUTLIER(maptype,x,y)){
	      double Bias,Pen,PenSm;
	      SintDetail(x,y, J-HEnd,I-K-GEnd,J,I,K,TEnd,PRtabY,Y,c,Bias,Pen,PenSm, 0);
	      outscore[U+1] = Bias + Pen + PenSm;
	    } else
	      outscore[U+1] = iscore[U+1];
	  } else {/* left end segment */
	    // NOTE : if PoutlierEnd > 0.0, A->G(I,K,J) <= -2 may be caused either by local alignment or PoutlierEnd
	    outscore[U+1] = iscore[U+1] = A->scoreEnd(I,K,J);
	  }
	  I = GEnd;
	  K = TEnd;
	  J = HEnd;
	  if(DEBUG>=2 && I > 0) assert(J >= 1 && J <= M);
	}
	U++;
      }

      if(DEBUG) assert(U > 0);
      if(DEBUG) assert(I <= (isLinear ? -1 : 0));
      if(DEBUG) assert((bestL[c] >= -1) ? I >= -1 : I <= -2);

      //      DP *p = &A[bestI[c]][bestK[c]][bestJ[c]];
      int bI = bestI[c], bK = bestK[c], bJ = bestJ[c];
      align[c].score = ((bestL[c] >= -1) ? A->scoreN(bI,bK,bJ) : A->scoreEnd(bI,bK,bJ)) + ((bestR[c] >= -1) ? A->Rscore(bI,bK,bJ) : A->RscoreEnd(bI,bK,bJ));
      if(VERB>=2 && mapid==0 && refid==0 && orientation){
	printf("refid=%d,mapid=%d,or=%d,c=%d:score=%0.6f,Lscore=%0.6f,Rscore=%0.6f:bI=%d,bK=%d,bJ=%d\n",
	       refid,mapid,orientation,c,align[c].score, ((bestL[c] >= -1) ? A->scoreN(bI,bK,bJ) : A->scoreEnd(bI,bK,bJ)), ((bestR[c] >= -1) ? A->Rscore(bI,bK,bJ) : A->RscoreEnd(bI,bK,bJ)),bI,bK,bJ);
	fflush(stdout);
      }
      align[c].orientation = orientation;
      align[c].scaleID = scaleID;
      align[c].logPV = alignFPsd(A,bI,bJ,bK,Y,X,N,M,bestR[c],bestL[c],orientation,refid,mapid,c,0);
      align[c].Lend = I;
      align[c].Rend = bestR[c];

      I = Ilist[U-1];
      K = Klist[U-1];
      J = Jlist[U-1];

      if(DEBUG>=2) assert(align[c].Lend <= -2 || J <= DELTA_X || (extend && I-K <= DELTA_Y));
      align[c].Lij1 = (align[c].Lend <= -2) ? I-K : A->Lij(I,K,J);
      if(DEBUG) assert(align[c].Lij1 <= I-K && align[c].Lij1 >= 0);
      align[c].Lij2 = (align[c].Lend <= -2) ? J : extend ? A->Lijx(I,K,J) : 0;
      if(DEBUG) assert(align[c].Lij2 <= J && align[c].Lij2 >= 0);
      if(DEBUG>=2) assert(align[c].Rend <= -2 || bJ >= M+1 - DELTA_X || (extend && bI >= N+1-DELTA_Y));
      align[c].Rij1 = (align[c].Rend <= -2) ? bI : A->Rij(bI,bK,bJ);
      if(DEBUG && !(align[c].Rij1 >= bI && align[c].Rij1 <= N+1)){
	printf("refalignXYsd:mapid=%d,orientation=%d,c=%d:bestI=%d,bestJ=%d,bestK=%d,bestR=%d,bestL=%d,bestscore=%0.6f,Lscore=%0.6f(End=%0.6f),Rscore=%0.6f(End=%0.6f):align->Rij1=%d,N=%d,M=%d,align->Rend=%d\n",
	       mapid,orientation,c,bestI[c],bestJ[c],bestK[c],bestR[c],bestL[c],bestscore,A->scoreN(bI,bK,bJ),A->scoreEnd(bI,bK,bJ),
	       A->Rscore(bI,bK,bJ),A->RscoreEnd(bI,bK,bJ),align[c].Rij1,N,M,align[c].Rend);
	printf("A[bestI][bestK][bestJ].Rij=%d\n",A->Rij(bI,bK,bJ));
	fflush(stdout);
	assert(align[c].Rij1 >= bestI[c] && align[c].Rij1 <= N+1);
      }
      align[c].Rij2 = (align[c].Rend <= -2) ? bestJ[c] : extend ? A->Rijx(bI,bK,bJ) : M+1;
      if(DEBUG && !(align[c].Rij2 >= bestJ[c] && align[c].Rij2 <= M+1)){
	printf("refalignXYsd:mapid=%d,orientation=%d,c=%d:bestI=%d,bestJ=%d,bestK=%d,bestR=%d,bestL=%d,bestscore=%0.6f,Lscore=%0.6f(End=%0.6f),Rscore=%0.6f(End=%0.6f):align->Rij2=%d,N=%d,M=%d,align->Rend=%d\n",
	       mapid,orientation,c,bestI[c],bestJ[c],bestK[c],bestR[c],bestL[c],bestscore,A->scoreN(bI,bK,bJ),A->scoreEnd(bI,bK,bJ),A->Rscore(bI,bK,bJ),A->RscoreEnd(bI,bK,bJ),align[c].Rij2,N,M,align[c].Rend);
	printf("A[bestI][bestK][bestJ].Rijx=%d,Kmax[bestI]=%d\n",A->Rijx(bI,bK,bJ),Kmax[bestI[c]]);
	fflush(stdout);
	assert(align[c].Rij2 >= bestJ[c] && align[c].Rij2 <= M+1);
      }

      /* reallocate best alignment */
      align[c].allfree();
      if(DEBUG) assert(U > 0);
      align[c].expand_arrays(U);
      align[c].numpairs = U;

      /* copy alignment from Ilist[],Jlist[] in reverse order */
      for(int I=0;I<U;I++){
	align[c].sites1[I] = Ilist[U-1-I];
	align[c].sitesK1[I] = Klist[U-1-I];
	align[c].sites2[I] = Jlist[U-1-I];
	align[c].iscore[I] = iscore[U-I];
	align[c].outscore[I] = outscore[U-I];
      }
      align[c].iscore[U] = iscore[0];
      align[c].outscore[U] = outscore[0];

      align[c].maxoutlier = 0.0;
      align[c].noutliers = 0;
      align[c].maxoutlierLabels = 0;
      int lastI = align[c].sites1[0];
      int lastK = align[c].sitesK1[0];
      int lastJ = align[c].sites2[0];
      for(int t = 1; t < U; lastI = I, lastK = K, lastJ = J, t++){
	int I = align[c].sites1[t];
	int K = align[c].sitesK1[t];
	int J = align[c].sites2[t];
	if(align[c].outscore[t] + (RFLOAT)0.01 < align[c].iscore[t]){
	  align[c].noutliers++;
	  double deltaY = Yc(Y,I,K) - Yc(Y,lastI,lastK);
	  double deltaX = X[J] - X[lastJ];
	  if(DEBUG) assert(deltaX > 0.0);
	  double delta = fabs(deltaY-deltaX);
	  align[c].maxoutlier = max(delta,align[c].maxoutlier);
	  align[c].maxoutlierLabels = max(I-K-lastI + J-lastJ -2, align[c].maxoutlierLabels);
	}
      }
    } // c = 0 .. colors-1

    /* compute combined score(s) : only score,orientation,logPV,numpairs are valid for align[-1] */
    align[-1].score = bestscore;
    align[-1].orientation = orientation;
    align[-1].scaleID = scaleID;
    align[-1].logPV = align[0].logPV + align[1].logPV;/* excludes misalignment information */
    if(VERB>=2 && refmap[refid]->id == 1 && Gmap[mapid]->id == 101000333LL && !orientation && align->align2){
      Calign *align2 = &align->align2[1];
      printf("A:refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:align=%p(logPV=%0.2f:%0.2f,%0.2f),align2=%p(logPV=%0.2f:%0.2f,%0.2f)\n",
	     refid,refmap[refid]->id,mapid,Gmap[mapid]->id,orientation,
	     align, align[-1].logPV, align[0].logPV, align[1].logPV,
	     align2,align2[-1].logPV,align2[0].logPV,align2[1].logPV);
      fflush(stdout);
    }

    align[-1].numpairs = align[0].numpairs + align[1].numpairs;

    /* HERE    if(RepeatMaxShift > 0 && giter == RefRepeats-1 && align[-1].numpairs >= AlignedSiteThreshold && AlignedEndOutlier(align) && align[-1].score > ScoreThreshold && align[-1].logPV > LogPvThreshold)
	  align[1-].repeat = alignRepeatSD(align,A,bestI,bestJ,bestK,Y,X,N,M,bestR,orientation,refid,mapid);*/

    if(VERB>= (PVERB?1:2) && bestscore > MINSCORE && (VERB/* HERE >=2 */ ||  (bestscore > ScoreThreshold && align[-1].logPV > LogPvThreshold && align[-1].numpairs >= AlignedSiteThreshold)) /* fabs(offset-trueoffset) > 2*Ylambda || align->orientation*/){
      #pragma omp critical 
      { 
	//    double bestYik = Yc(YY[0],bestI[0],bestK[0]);
	double offset = bestYik - XX[0][bestJ[0]];
	Cmap *pmap = Gmap[mapid];
	while(pmap->origmap)
	  pmap = pmap->origmap;

        printf("refid=%d,mapid=%d(id=%lld),or=%d:score=%0.6f(c1=%0.6f,I=%d,K=%d,J=%d,L=%d,R=%d,c2=%0.6f,P=%d,T=%d,Q=%d,L=%d,R=%d,ma=%0.6f,y=%0.3f,x=%0.3f),logPV=%0.4f(%0.4f,%0.4f):orientation=%d,scaleID=%d(%0.4f),offset=%0.3f,aligned sites=%d+%d, Lend=%d,%d,Rend=%d,%d,Cscore=%0.4e,fpcnt=%d(%d):%s\n",
	       refid,mapid,Gmap[mapid]->id,orientation,bestscore,align[0].score,bestI[0],bestK[0],bestJ[0],bestL[0],bestR[0],align[1].score,bestI[1],bestK[1],bestJ[1],bestL[1],bestR[1],bestscore-align[0].score-align[1].score,
	       Yc(YY[1],bestI[1],bestK[1])-Yc(YY[0],bestI[0],bestK[1]),XX[1][bestJ[1]]-XX[0][bestJ[0]],
	       align[-1].logPV,align[0].logPV,align[1].logPV,orientation,scaleID,scaleID ? ScaleFactor[scaleID] : 1.0,
	       offset,align[0].numpairs,align[1].numpairs,align[0].Lend,align[1].Lend,align[0].Rend,align[1].Rend,ChimScore,
	       MM[0]+MM[1]-align[-1].numpairs, Gmap[mapid]->fpcnt, pmap->name ? pmap->name : "");
	double y0 = Yc(YY[0],bestI[0],bestK[0]);
	double y1 = Yc(YY[1],bestI[1],bestK[1]);
	double x0 = XX[0][bestJ[0]];
	double x1 = XX[1][bestJ[1]];

	printf("SMA(%0.4f,%0.4f,%d)=%0.6f:MA_mean=%0.4f,SF[2]=%0.8e,FM2var=%0.8e,SD[2]=%0.8e,SM2var=%0.8e\n",
	       XX[1][bestJ[1]]-XX[0][bestJ[0]],Yc(YY[1],bestI[1],bestK[1])-Yc(YY[0],bestI[0],bestK[1]),orientation,
	       SMA(x1-x0,y1-y0,(COLOR_SHIFT==2) ? max(y0,y1) : max(x0,x1),orientation), MA_mean,SF[2],FM2var,SD[2],SM2var);

	for(int c = 0; c < colors; c++){
	  FLOAT *Y = YY[c];
	  FLOAT *X = XX[c];
	  int N = NN[c];
	  int M = MM[c];
	  AL2 *A = &AA[c];
	  Cprtab **PRtabY = PRtab[c][refid];

	  printf("color=%d:",c+1);
	  (void)alignFPsd(A,bestI[c],bestJ[c],bestK[c],Y,X,N,M,bestR[c],bestL[c],orientation,refid,mapid,c,1);

	  int I= -1,J= -1 ,K=0, G= -1,H= -1, D = 0;
	  int U = align[c].numpairs;
	  for(int T=0; T < U; T++, G=I, H=J, D=K){
	    I = align[c].sites1[T];/* Ilist[U-1-T] */
	    K = align[c].sitesK1[T];/* Klist[U-1-T] */
	    J = align[c].sites2[T];/* Jlist[U-1-T] */
	    //	    DP *p = &A[I][K][J];
	    printf("T=%2d:I=%2d,K=%2d,J=%2d:",T,I,K,J);
	    if(T<=0){ /* Left End: display Send(min(X[J],Y[I-K..I]),J+1-max(1,Lijx),I+1-max(1,Lij),c)+Sm(J,I,K), A[I,K,J].score */
	      if(DEBUG && align[c].Lend > -2) assert(J <= DELTA_X || (extend && I-K <= DELTA_Y));
	      int Lij = (align[c].Lend <= -2) ? I-K : A->Lij(I,K,J);
	      if(DEBUG && align[c].Lend > -2 && extend) assert(I-K <= DELTA_Y || J <= DELTA_X);
	      int Lijx = (align[c].Lend <= -2) ? J : extend ? A->Lijx(I,K,J) : 0;
	      if(align->Lend <= -2 && (NEW==1 || align->Lend == -2))
		printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I,K]=%6.3f,C=%0.6f(%d),Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
		       Lij,Lijx,X[J],Yc(Y,I,K),ChimScore,localtype,Sm(J,I,K,Y,c),ChimScore+Sm(J,I,K,Y,c),(bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
	      else if(align->Lend <= -3)
		printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I,K]=%6.3f,EndPen=%0.6f,BiasEnd=%0.6f,Sm(J,I,K)=%0.6f,score=%0.6f,A[I,K,J]=%0.6f\n",
		       Lij,Lijx,X[J],Yc(Y,I,K),OutlierEndPen,BiasEnd2(X[J],J,c),Sm(J,I,K,Y,c),BiasEnd2(X[J],J,c)+OutlierEndPen+Sm(J,I,K,Y,c),(bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
	      else if(ENDFIX && X[J] <= Yc(Y,I,K) && Lij > 0 && Yc(Y,I,K) - Y[Lij-1] < X[J]){
		TFLOAT score = Sbnd(X[J],Yc(Y,I,K)-Y[Lij-1],J,I-K+1-Lij,c);
		printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I]=%6.3f,Sbnd(X[J],Y[I,K]-Y[Lij-1],J,I-K+1-Lij)=%0.6f,Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
		       Lij,Lijx,X[J],Yc(Y,I,K),score,Sm(J,I,K,Y,c), max(ChimScore,score+(TFLOAT)Sm(J,I,K,Y,c)), (bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		if(DEBUG) assert(Lijx <= 0);
	      } else if(ENDFIX && extend && X[J] >= Yc(Y,I,K) && Lijx > 0 && X[J]-X[Lijx-1] < Yc(Y,I,K)){
		TFLOAT score = Sbnd(X[J]-X[Lijx-1],Yc(Y,I,K),J+1-Lijx,I-K,c);
		printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I,K]=%6.3f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
		       Lij,Lijx,X[J],Y[I],X[J]-X[Lijx-1],Yc(Y,I,K),J+1-Lijx,I-K, score,Sm(J,I,K,Y,c), max(ChimScore,score+(TFLOAT)Sm(J,I,K,Y,c)), (bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		if(DEBUG) assert(Lij <= 0);
	      } else {
		TFLOAT score = Send(X[J],J+1-max(1,Lijx) ,I-K+1 - max(1,Lij),c);
		printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I,K]=%6.3f,Send(min(X[J],Y[I,K]),J+1-max(1,Lijx),I-K+1-max(1,Lij))=%0.6f,Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
		       Lij,Lijx,X[J],Yc(Y,I,K),score,Sm(J,I,K,Y,c), max(ChimScore,Lij > 0 ? score+(TFLOAT)Sm(J,I,K,Y,c) : (TFLOAT)MINSCORE), (bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
	      }
	    } else {/* internal interval : display Sint(X[J]-X[H],Y[I,K]-Y[G,D],J-H,I-K-G) + Sm(J,I,K), A[I,K,J].score */
	      if(DEBUG) assert(G>=0 && H>=0 && D>=0);
	      TFLOAT x = X[J]-X[H];
	      TFLOAT y = Yc(Y,I,K)-Yc(Y,G,D);
	      TFLOAT Pen,Bias,PenSm;
	      SintDetail(x,y,J-H,I-K-G,J,I,K,D,PRtabY,Y,c,Bias,Pen,PenSm,0);
	      TFLOAT OutPen = OutlierPenalty[c];
	      if(outlierBC)
		OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-K-G];
	      //	      double sint = Sint(x,y,J-H,I-K-G,J,I,K,D,PRtabY,Y,c);
	      TFLOAT sint = Bias + OutlierBias + (FIX_OUTLIER_MIS ? max(Pen,OutPen) + PenSm : max(Pen+PenSm, OutPen));
	      double var = SF[c]*SF[c] + fabs(SD[c])*SD[c] * y;
	      if(QUADRATIC_VARIANCE)
		var += SR[c]*SR[c]*y*y;
	      if(RES_VARIANCE){
		FLOAT resR = Y[I] - Y[I-K];
		FLOAT resL = Y[G] - Y[G-D];
		FLOAT resvar = resL*resL + resR*resR;
		var += SE[c]*SE[c]*resvar;
	      }
	      if(align[c].orientation)
		printf("G=%2d,H=%2d,D=%2d,X[J]=%0.3f,Y[I,K]=%0.3f,X[J]-X[H]=%0.3f,Y[I,K]-Y[G,D]=%0.3f,norm=%0.4f,Sint(X,Y,%d,%d)=%0.6f(Pen=%0.6f,Bias=%0.6f,Sm=%0.6f),A[I,K,J]=%0.6f,(X[M+1]-X[J])/C=%0.3f\n",
		       G,H,D,X[J],Yc(Y,I,K),x,y,(x-y)*(x-y)/var,J-H,I-K-G,sint,Pen,Bias,PenSm,(bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J), (X[M+1]-X[J])*origPixelLen/PixelLen);
	      else
		printf("G=%2d,H=%2d,D=%2d,X[J]=%0.3f,Y[I,K]=%0.3f,X[J]-X[H]=%0.3f,Y[I,K]-Y[G,D]=%0.3f,norm=%0.4f,Sint(X,Y,%d,%d)=%0.6f(Pen=%0.6f,Bias=%0.6f,Sm=%0.6f),A[I,K,J]=%0.6f,X[J]/C=%0.3f\n",
		       G,H,D,X[J],Yc(Y,I,K),x,y,(x-y)*(x-y)/var,J-H,I-K-G,sint,Pen,Bias,PenSm,(bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J), X[J]*origPixelLen/PixelLen);
	    }
	  }
	  //	  DP *p = &A[I][K][J];      
	  if(DEBUG && align[c].Rend > -2) assert(J >= ((extend && I >= N+1-DELTA_Y) ? 1 : max(1,M+1-DELTA_X)));
	  int Rij = (align[c].Rend <= -2) ? I : A->Rij(I,K,J);
	  if(DEBUG && extend && align[c].Rend > -2) assert(I>=N+1-DELTA_Y || J >= M+1-DELTA_X);
	  int Rijx = (align[c].Rend <= -2) ? J : extend ? A->Rijx(I,K,J) : M+1;

	  /* Right End: display Rij,Send(X[M+1]-X[J],min(M,Rijx)+1-J,min(N,Rij)+1-I,c) */
	  if(align->Rend <= -2 && (NEW==1 || align->Rend == -2))
	    printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rijk=%d,Rijkx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,C=%0.6f(%d),total score=%0.6f, X[M+1]/C=%0.3f\n",
		   I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), ChimScore, localtype, ((bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + ChimScore, X[M+1]*origPixelLen/PixelLen);
	  else if(align->Rend <= -3)
	    printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rijk=%d,Rijkx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,EndPen=%0.6f,BiasEnd=%0.6f,total score=%0.6f, X[M+1]/C=%0.3f\n",
		   I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), OutlierEndPen, BiasEnd2(X[M+1]-X[J],M+1-J,c), 
		   ((bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + OutlierEndPen + BiasEnd2(X[M+1]-X[J],M+1-J,c), X[M+1]*origPixelLen/PixelLen);
	  else if(ENDFIX && X[M+1]-X[J] <= Y[N+1]-Yc(Y,I,K) && Rij <= N && Y[Rij+1]-Yc(Y,I,K) < X[M+1]-X[J]){
	    TFLOAT score = Sbnd(X[M+1]-X[J],Y[Rij+1]-Yc(Y,I,K),min(M,Rijx)+1-J,Rij+1-I,c);
	    printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rij=%d,Rijx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Yc(Y,I,K)=%0.3f,Sbnd(X[M+1]-X[J],Y[Rij+1]-Y[I,K],min(M,Rijx)+1-J,Rij+1-I)=%0.6f,Rscore=%0.6f,total score=%0.6f,X[M+1]/C=%0.3f\n",
		   I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), score, ((bestR[c] >= -1) ? A->Rscore(I,K,J) : A->RscoreEnd(I,K,J)),
		   ((bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + ((bestR[c] >= -1) ? A->Rscore(I,K,J) : A->RscoreEnd(I,K,J)), X[M+1]*origPixelLen/PixelLen);
	  } else if(ENDFIX && X[M+1]-X[J] >= Y[N+1]-Yc(Y,I,K) && Rijx <= M && X[A->Rijx(I,K,J)+1]-X[J] < Y[N+1]-Yc(Y,I,K)){
	    TFLOAT score = Sbnd(X[Rijx+1]-X[J],Y[N+1]-Yc(Y,I,K),min(M,Rijx)+1-J,min(N,Rij)+1-I,c);
	    printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rij=%d,Rijx=%d,X[M+1]=%0.4f,Y[N+1]=%0.4f,X[Rijx+1]-X[J]=%0.4f,Y[N+1]-Y[I,K]=%0.4f,Sbnd(X[Rijx+1]-X[J],Y[N+1]-Y[I,K],Rijx+1-J,min(N,Rij)+1-I)=%0.6f,Rscore=%0.6f,total score=%0.6f,X[M+1]/C=%0.3f\n",
		   I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[Rijx+1]-X[J],Y[N+1]-Yc(Y,I,K), score, ((bestR[c] >= -1) ? A->Rscore(I,K,J) : A->RscoreEnd(I,K,J)),
		   ((bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + ((bestR[c] >= -1) ? A->Rscore(I,K,J) : A->RscoreEnd(I,K,J)) , X[M+1]*origPixelLen/PixelLen);
	  } else {
	    TFLOAT score = Send(min(X[M+1]-X[J],Y[N+1]-Yc(Y,I,K)),min(M,Rijx)+1-J,min(N,Rij)+1-I,c);
	    printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rij=%d,Rijx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,Send(min(X[M+1]-X[J],Y[N+1]-Y[I,K]),min(M,Rijx)+1-J,min(N,Rij)+1-I)=%0.6f,Rscore=%0.6f,total score=%0.6f,X[M+1]/C=%0.3f\n",
		   I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), score, ((bestR[c] >= -1) ? A->Rscore(I,K,J) : A->RscoreEnd(I,K,J)),
		   ((bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + ((bestR[c] >= -1) ? A->Rscore(I,K,J) : A->RscoreEnd(I,K,J))/* max(ChimScore, (!extend && Rij > N) ? (TFLOAT)MINSCORE : score)*/,
		   X[M+1]*origPixelLen/PixelLen);
	  }
	  fflush(stdout);
	  if(DEBUG) assert(Rij >= I && Rij <= N+1);
	}
      } // pragma omp critical
    } // if(VERB...) 
  }

  if(SecondBest && bestscore2 > MINSCORE && (!align[-1].align2 || bestscore2 > align[-1].align2->score + 1e-8)){/* check if 2nd best alignment should be saved as overall 2nd best alignment */
    int U1 = align->numpairs;
    int I2 = align->sites1[U1 - 1];
    int K2 = align->sitesK1[U1 - 1];
    FLOAT alignYik = Yc(YY[0],I2,K2);
    if(align[-1].orientation != orientation || fabs(alignYik - bestYik2) >= MinSpacing){   /* backtrack through array for each color to retrieve 2nd best alignment starting at right end AA[c][bestI2[c]][bestK2[c]][bestJ2[c]] */

      delete [] align[-1].align2; align[-1].align2 = new Calign[3];
      Calign *align2 = &align[-1].align2[1];
      align2[-1].mapid1 = refid;
      align2[-1].mapid2 = mapid;

      if(DEBUG) assert(bestI2[0] >= 0 && bestI2[1] >= 0);
      if(DEBUG){
	//	DP *p[2];
	//	for(int c = 0; c < colors; c++)
	//	  p[c] = &AA[c][bestI2[c]][bestK2[c]][bestJ2[c]];
	double y0 = Yc(YY[0],bestI2[0],bestK2[0]);
	double y1 = Yc(YY[1],bestI2[1],bestK2[1]);
	double x0 = XX[0][bestJ2[0]];
	double x1 = XX[1][bestJ2[1]];
	double score = SMA(x1-x0, y1-y0, (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0), orientation);
	double Lscore = AA[0].scoreN(bestI2[0],bestK2[0],bestJ2[0]) + AA[1].scoreN(bestI2[1],bestK2[1],bestJ2[1]) + score;
	double LscoreEnd = AA[0].scoreEnd(bestI2[0],bestK2[0],bestJ2[0]) + AA[1].scoreEnd(bestI2[1],bestK2[1],bestJ2[1]) + score;
	double Rscore = AA[0].Rscore(bestI2[0],bestK2[0],bestJ2[0]) + AA[1].Rscore(bestI2[1],bestK2[1],bestJ2[1]);
	double RscoreEnd = AA[0].RscoreEnd(bestI2[0],bestK2[0],bestJ2[0]) + AA[1].RscoreEnd(bestI2[1],bestK2[1],bestJ2[1]);

	assert((bestR2[0] >= -1) ? (bestR2[1] >= -1) : (bestR2[1] < -1));
	assert((bestR2[0] >= -1) ? (bestR2[0] == AA[0].R(bestI2[0],bestK2[0],bestJ2[0])) : (bestR2[0] == AA[0].REnd(bestI2[0],bestK2[0],bestJ2[0])));
	assert((bestR2[1] >= -1) ? (bestR2[1] == AA[1].R(bestI2[1],bestK2[1],bestJ2[1])) : (bestR2[1] == AA[1].REnd(bestI2[1],bestK2[1],bestJ2[1])));
	assert(bestL2[0] == bestL2[1]);

	double bscore = ((bestR2[0] >= -1) ? Rscore : RscoreEnd) + ((bestL2[0] >= -1) ? Lscore : LscoreEnd);
	if(!(fabs(bscore - bestscore2) < 1e-8)){
	  printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:bestI2=%d,%d,bestK2=%d,%d,bestJ2=%d,%d,bestR2=%d,%d,bestL2=%d,%d:\n",
		 refid,refmap[refid]->id,mapid,Gmap[mapid]->id,orientation,bestI2[0],bestI2[1],bestK2[0],bestK2[1],bestJ2[0],bestJ2[1],bestR2[0],bestR2[1],bestL2[0],bestL2[1]);
	  printf("    scoreN=%0.6f+%0.6f,score=%0.6f,scoreEnd=%0.6e,%0.6e,Rscore=%0.6f,%0.6f,RscoreEnd=%0.6e,%0.6e,bscore=%0.6f,bestscore2=%0.6f\n",
		 AA[0].scoreN(bestI2[0],bestK2[0],bestJ2[0]),AA[1].scoreN(bestI2[1],bestK2[1],bestJ2[1]),score,AA[0].scoreEnd(bestI2[0],bestK2[0],bestJ2[0]),
		 AA[1].scoreEnd(bestI2[1],bestK2[1],bestJ2[1]),AA[0].Rscore(bestI2[0],bestK2[0],bestJ2[0]),AA[1].Rscore(bestI2[1],bestK2[1],bestJ2[1]),
		 AA[0].RscoreEnd(bestI2[0],bestK2[0],bestJ2[0]),AA[1].RscoreEnd(bestI2[1],bestK2[1],bestJ2[1]),bscore,bestscore2);
	  printf("    XX[0][bestJ2[0]]=%0.4f,XX[1][bestJ2[1]]=%0.4f,y0=%0.4f,y1=%0.4f\n",
		 XX[0][bestJ2[0]],XX[1][bestJ2[1]],y0,y1);
	  fflush(stdout);
	  assert(fabs(bscore - bestscore2) < 1e-8);
	}
      }

      if(!Ilist){    /* allocate space to store alignment in reverse order */
	int M = max(MM[0],MM[1]);
	Ilist = new int[M];
	Klist = new int[M];
	Jlist = new int[M];
	iscore = new double[M+1];
	outscore = new double[M+1];
      }

      for(int c = 0; c < colors; c++){
	FLOAT *Y = YY[c];
	FLOAT *X = XX[c];
	int N = NN[c];
	int M = MM[c];
	AL2 *A = &AA[c];
	Cprtab **PRtabY = PRtab[c][refid];

	int I = bestI2[c];
	int K = bestK2[c];
	int J = bestJ2[c];

	int U = 0;

	/* right end segment */
	// NOTE : if PoutlierEnd > 0.0, bestR2[c] <= -2 may be caused either by local alignment or PoutlierEnd
	outscore[0] = iscore[0] = (bestR2[c] >= -1 ? A->Rscore(I,K,J) : A->RscoreEnd(I,K,J));
	while(I > 0){
	  //	  DP *p = &A[I][K][J];
	  if(DEBUG){
	    assert(I-K >= 1 && I <= N);
	    assert(J >= 1 && J <= M);
	  }
	  Ilist[U] = I;
	  Klist[U] = K;
	  Jlist[U] = J;
	  if(bestL2[c] >= -1){
	    int G = A->G(I,K,J);
	    int T = A->T(I,K,J);
	    int H = A->H(I,K,J);
	    if(G > 0){
	      if(DEBUG) assert(T >= 0 && H > 0);
	      iscore[U+1] = A->scoreN(I,K,J) - A->scoreN(G,T,H);
	      double x = X[J] - X[H];
	      double y = Yc(Y,I,K) - Yc(Y,G,T);
	      if(Poutlier > 0.0 && OUTLIER(maptype,x,y)){
		double Bias,Pen,PenSm;
		SintDetail(x,y, J-H,I-K-G,J,I,K,T,PRtabY,Y,c,Bias,Pen,PenSm,0);
		outscore[U+1] = Bias + Pen + PenSm;
	      } else
		outscore[U+1] = iscore[U+1];
	    } else {/* left end segment */
	      // NOTE : if PoutlierEnd > 0.0, G <= -2 may be caused either by local alignment or PoutlierEnd
	      outscore[U+1] = iscore[U+1] = A->scoreN(I,K,J);
	    }
	    I = G;
	    K = T;
	    J = H;
	  } else {/* bestL2[c] <= -2 */
	    int GEnd = A->GEnd(I,K,J);
	    int TEnd = A->TEnd(I,K,J);
	    int HEnd = A->HEnd(I,K,J);
	    if(GEnd > 0){
	      if(DEBUG) assert(TEnd >= 0 && HEnd > 0);
	      iscore[U+1] = A->scoreEnd(I,K,J) - A->scoreEnd(GEnd,TEnd,HEnd);
	      double x = X[J]-X[HEnd];
	      double y = Yc(Y,I,K) - Yc(Y,GEnd,TEnd);
	      if(Poutlier > 0.0 && OUTLIER(maptype,x,y)){
		double Bias,Pen, PenSm;
		SintDetail(x,y, J-HEnd,I-K-GEnd,J,I,K,TEnd,PRtabY,Y,c,Bias,Pen,PenSm,0);
		outscore[U+1] = Bias + Pen + PenSm;
	      } else
		outscore[U+1] = iscore[U+1];
	    } else {/* left end segment */
	      // NOTE : if PoutlierEnd > 0.0, GEnd <= -2 may be caused either by local alignment or PoutlierEnd
	      outscore[U+1] = iscore[U+1] = A->scoreEnd(I,K,J);
	    }
	    I = GEnd;
	    K = TEnd;
	    J = HEnd;
	  }
	  U++;
	}

	if(DEBUG) assert(U > 0);
	if(DEBUG && !(I <= (isLinear ? -1 : 0))){
	  printf("refid=%d,mapid=%d,or=%d:c=%d:I=%d,K=%d,J=%d,isLinear=%d,bestL2=%d,%d\n",
		 refid,mapid,orientation,c,I,K,J,isLinear,bestL2[0],bestL2[1]);
	  printf("   bestI2=%d,%d,bestK2=%d,%d,bestJ2=%d,%d:U=%d:\n",bestI2[0],bestI2[1],bestK2[0],bestK2[1],bestJ2[0],bestJ2[1],U);
	  for(int t = 0; t < U; t++){
	    int I = Ilist[U-1-t];
	    int K = Klist[U-1-t];
	    int J = Jlist[U-1-t];
	    //	    DP *p = &A[I][K][J];
	    printf("t=%d:I=%d,K=%d,J=%d:G=%d,T=%d,H=%d,GEnd=%d,TEnd=%d,HEnd=%d,iscore=%0.6f,outscore=%0.6f\n",
		   t,I,K,J, A->G(I,K,J),A->T(I,K,J),A->H(I,K,J), A->GEnd(I,K,J),A->TEnd(I,K,J),A->HEnd(I,K,J), iscore[U-t], outscore[U-t]);
	  }
	  printf("t=%d: iscore=%0.6f, outscore=%0.6f\n", U, iscore[0], outscore[0]);
	  fflush(stdout);
	  assert(I <= (isLinear ? -1 : 0));
	}
	if(DEBUG) assert((bestL2[c] >= -1) ? I >= -1 : I <= -2);

	//	DP *p = &A[bestI2[c]][bestK2[c]][bestJ2[c]];
	int bI2 = bestI2[c], bK2 = bestK2[c], bJ2 = bestJ2[c];
	align2[c].score = ((bestL2[c] >= -1) ? A->scoreN(bI2,bK2,bJ2) : A->scoreEnd(bI2,bK2,bJ2)) + (bestR2[c] >= -1 ? A->Rscore(bI2,bK2,bJ2) : A->RscoreEnd(bI2,bK2,bJ2));
	align2[c].orientation = orientation;
	align2[c].scaleID = scaleID;
	align2[c].logPV = alignFPsd(A,bI2,bJ2,bK2,Y,X,N,M,bestR2[c],bestL2[c],orientation,refid,mapid,c,0);
	align2[c].Lend = I;
	align2[c].Rend = bestR2[c];

	I = Ilist[U-1];
	K = Klist[U-1];
	J = Jlist[U-1];

	if(DEBUG>=2) assert(align2[c].Lend <= -2 || J <= DELTA_X || (extend && I-K <= DELTA_Y));
	align2[c].Lij1 = (align2[c].Lend <= -2) ? I-K : A->Lij(I,K,J);
	if(DEBUG) assert(align2[c].Lij1 <= I-K && align2[c].Lij1 >= 0);
	align2[c].Lij2 = (align2[c].Lend <= -2) ? J : extend ? A->Lijx(I,K,J) : 0;
	if(DEBUG) assert(align2[c].Lij2 <= J && align2[c].Lij2 >= 0);
	if(DEBUG>=2) assert(align2[c].Rend <= -2 || bestJ2[c] >= M+1 - DELTA_X || (extend && bestI2[c] >= N+1-DELTA_Y));
	align2[c].Rij1 = (align2[c].Rend <= -2) ? bestI2[c] : A->Rij(bI2,bK2,bJ2);
	if(DEBUG) assert(align2[c].Rij1 >= bestI2[c] && align2[c].Rij1 <= N+1);
	align2[c].Rij2 = (align2[c].Rend <= -2) ? bestJ2[c] : extend ? A->Rijx(bI2,bK2,bJ2) : M+1;
	if(DEBUG)	assert(align2[c].Rij2 >= bestJ2[c] && align2[c].Rij2 <= M+1);

	/* reallocate best alignment */
	align2[c].allfree();
	if(DEBUG) assert(U > 0);
	align2[c].expand_arrays(U);
	align2[c].numpairs = U;

	/* copy alignment from Ilist[],Jlist[] in reverse order */
	for(int I=0;I<U;I++){
	  align2[c].sites1[I] = Ilist[U-1-I];
	  align2[c].sitesK1[I] = Klist[U-1-I];
	  align2[c].sites2[I] = Jlist[U-1-I];
	  align2[c].iscore[I] = iscore[U-I];
	  align2[c].outscore[I] = outscore[U-I];
	}
	align2[c].iscore[U] = iscore[0];
	align2[c].outscore[U] = outscore[0];
	
	align2[c].maxoutlier = 0.0;
	align2[c].noutliers = 0;
	align2[c].maxoutlierLabels = 0;
	int lastI = align2[c].sites1[0];
	int lastK = align2[c].sitesK1[0];
	int lastJ = align2[c].sites2[0];
	for(int t = 1; t < U; lastI = I, lastK = K, lastJ = J, t++){
	  int I = align2[c].sites1[t];
	  int K = align2[c].sitesK1[t];
	  int J = align2[c].sites2[t];
	  if(align2[c].outscore[t] + (RFLOAT)0.01 < align2[c].iscore[t]){
	    align2[c].noutliers++;
	    double deltaY = Yc(Y,I,K) - Yc(Y,lastI,lastK);
	    double deltaX = X[J] - X[lastJ];
	    if(DEBUG) assert(deltaX > 0.0);
	    double delta = fabs(deltaY-deltaX);
	    align2[c].maxoutlier = max(delta,align2[c].maxoutlier);
	    // HERE HERE HERE
	  }
	}

      }// c = 0 .. colors - 1

      /* compute combined score(s) : only score2,orientation2,logPV2,numpairs2 are valid for align[-1] */
      align2[-1].score = bestscore2;
      align2[-1].orientation = orientation;
      align2[-1].scaleID = scaleID;
      align2[-1].logPV = align2[0].logPV + align2[1].logPV;/* excludes misalignment information */
      if(VERB>=2 && refmap[refid]->id == 1 && Gmap[mapid]->id == 101000333LL && !orientation){
	printf("B:refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:align=%p(score=%0.6f,logPV=%0.2f:%0.2f,%0.2f),align2=%p(score=%0.6f,logPV=%0.2f:%0.2f,%0.2f)\n",
	       refid,refmap[refid]->id,mapid,Gmap[mapid]->id,orientation,
	       align, align[-1].score, align[-1].logPV, align[0].logPV, align[1].logPV,
	       align2, align2[-1].score, align2[-1].logPV,align2[0].logPV,align2[1].logPV);
	printf("    bestI=%d,%d,bestK=%d,%d,bestJ=%d,%d,bestR=%d,%d,bestscore=%0.6f; bestI2=%d,%d,bestK2=%d,%d,bestJ2=%d,%d,bestR2=%d,%d,bestscore2=%0.6f\n",
	       bestI[0],bestI[1],bestK[0],bestK[1],bestJ[0],bestJ[1],bestR[0],bestR[1],bestscore,
	       bestI2[0],bestI2[1],bestK2[0],bestK2[1],bestJ2[0],bestJ2[1],bestR2[0],bestR2[1],bestscore2);
	fflush(stdout);
      }
      align2[-1].numpairs = align2[0].numpairs + align2[1].numpairs;

      if(VERB>=2 && bestscore2 > MINSCORE && (VERB>=2 || (/* giter == RefRepeats-1 && */ mapid == 90 && refid == 0 && !orientation /* && bestscore > ScoreThreshold && align[-1].logPV > LogPvThreshold && align[-1].numpairs >= AlignedSiteThreshold */)/* ||  (mapid==1254 && align->orientation) */ /* || fabs(offset-trueoffset) > 2*Ylambda || align->orientation*/)){
        #pragma omp critical 
	{ 
	  double offset = bestYik2 - XX[0][bestJ2[0]];
	  double trueoffset = Gmap[mapid]->endloc - XX[0][MM[0]+1];
	  double delta = offset-trueoffset;
	  Cmap *pmap = Gmap[mapid];
	  while(pmap->origmap)
	    pmap = pmap->origmap;

	  printf("2nd Best:refid=%d,mapid=%d(id=%lld),score=%0.4f(c1=%0.4f,I=%d,K=%d,J=%d,R=%d,c2=%0.4f,P=%d,T=%d,Q=%d,S=%d,ma=%0.4f),logPV=%0.4f(%0.4f,%0.4f):orientation=%d,scaleid=%d(%0.4f),trueflip=%d,trueoffset=%0.3f,delta=%0.3f,aligned sites=%d+%d, Lend=%d,%d,Rend=%d,%d,fpcnt=%d(%d):%s\n",
		 refid,mapid,Gmap[mapid]->id,bestscore2,align2[0].score,bestI2[0],bestK2[0],bestJ2[0],bestR2[0],align2[1].score,bestI2[1],bestK2[1],bestJ2[1],bestR2[1],bestscore2 - align2[0].score - align2[1].score,
		 align2[-1].logPV,align2[0].logPV,align2[1].logPV,orientation,scaleID, scaleID ? ScaleFactor[scaleID] : 1.0,
		 Gmap[mapid]->flipped,trueoffset,delta,align2[0].numpairs,align2[1].numpairs,align2[0].Lend,align2[1].Lend,align2[0].Rend,align2[1].Rend,
		 MM[0]+MM[1]-align2[-1].numpairs, Gmap[mapid]->fpcnt,pmap->name ? pmap->name : "");

	  for(int c = 0; c < colors; c++){
	    FLOAT *Y = YY[c];
	    FLOAT *X = XX[c];
	    int N = NN[c];
	    int M = MM[c];
	    AL2 *A = &AA[c];
	    Cprtab **PRtabY = PRtab[c][refid];

	    printf("color=%d:\n",c+1);
	    fflush(stdout);

	    (void)alignFPsd(A,bestI2[c],bestJ2[c],bestK2[c],Y,X,N,M,bestR2[c],bestL2[c],orientation,refid,mapid,c,1);
	    int I= -1,J= -1 ,K=0, G= -1,H= -1, D = 0;
	    int U = align2[c].numpairs;

	    for(int T=0; T < U; T++, G=I, H=J, D=K){
	      I = align2[c].sites1[T];
	      K = align2[c].sitesK1[T];
	      J = align2[c].sites2[T];
	      //	      DP *p = &A[I][K][J];
	      printf("T=%2d:I=%2d,K=%2d,J=%2d:",T,I,K,J);
	      if(T<=0){ /* Left End: display Send(min(X[J],Y[I-K..I]),J+1-max(1,Lijx),I+1-max(1,Lij))+Sm(J,I,K), A[I,K,J].score */
		if(DEBUG && align2[c].Lend > -2 && !(J <= DELTA_X || (extend && I-K <= DELTA_Y))){
		  printf("c=%d,align2[c].Lend=%d,I=%d,K=%d,J=%d,extend=%d,DELTA_X=%d,DELTA_Y=%d\n",
			 c,align2[c].Lend,I,K,J,extend,DELTA_X,DELTA_Y);
		  fflush(stdout);
		  assert(J <= DELTA_X || (extend && I-K <= DELTA_Y));
		}
		int Lij = (align2[c].Lend <= -2) ? I-K : A->Lij(I,K,J);
		if(DEBUG && align2[c].Lend > -2 && extend) assert(I-K <= DELTA_Y || J <= DELTA_X);
		int Lijx = (align2[c].Lend <= -2) ? J : extend ? A->Lijx(I,K,J) : 0;

		if(align2[c].Lend <= -2 && (NEW==1 || align2[c].Lend == -2))
		  printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I,K]=%6.3f,C=%0.6f(%d),Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
			 Lij,Lijx,X[J],Yc(Y,I,K),ChimScore,localtype,Sm(J,I,K,Y,c),ChimScore+Sm(J,I,K,Y,c),(bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		else if(align2[c].Lend <= -3)
		  printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I,K]=%6.3f,EndPen=%0.6f,BiasEnd=%0.6f,Sm(J,I,K)=%0.6f,score=%0.6f,A[I,K,J]=%0.6f\n",
			 Lij,Lijx,X[J],Yc(Y,I,K),OutlierEndPen,BiasEnd2(X[J],J,c),Sm(J,I,K,Y,c),BiasEnd2(X[J],J,c)+OutlierEndPen+Sm(J,I,K,Y,c),(bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		else if(ENDFIX && X[J] <= Yc(Y,I,K) && Lij > 0 && Yc(Y,I,K) - Y[Lij-1] < X[J]){
		  TFLOAT score = Sbnd(X[J],Yc(Y,I,K)-Y[Lij-1],J,I-K+1-Lij,c);
		  printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I]=%6.3f,Sbnd(X[J],Y[I,K]-Y[Lij-1],J,I-K+1-Lij)=%0.6f,Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
			 Lij,Lijx,X[J],Yc(Y,I,K),score,Sm(J,I,K,Y,c), max(ChimScore,score+(TFLOAT)Sm(J,I,K,Y,c)), (bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		  if(DEBUG) assert(Lijx <= 0);
		} else if(ENDFIX && extend && X[J] >= Yc(Y,I,K) && Lijx > 0 && X[J]-X[Lijx-1] < Yc(Y,I,K)){
		  TFLOAT score = Sbnd(X[J]-X[Lijx-1],Yc(Y,I,K),J+1-Lijx,I-K,c);
		  printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I,K]=%6.3f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
			 Lij,Lijx,X[J],Y[I],X[J]-X[Lijx-1],Yc(Y,I,K),J+1-Lijx,I-K, score,Sm(J,I,K,Y,c), max(ChimScore,score+(TFLOAT)Sm(J,I,K,Y,c)), (bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		  if(DEBUG) assert(Lij <= 0);
		} else {
		  TFLOAT score = Send(X[J],J+1-max(1,Lijx) ,I-K+1 - max(1,Lij),c);
		  printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I,K]=%6.3f,Send(min(X[J],Y[I,K]),J+1-max(1,Lijx),I-K+1-max(1,Lij))=%0.6f,Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
			 Lij,Lijx,X[J],Yc(Y,I,K),score,Sm(J,I,K,Y,c), max(ChimScore,Lij > 0 ? score+(TFLOAT)Sm(J,I,K,Y,c) : (TFLOAT)MINSCORE), (bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		}
	      } else {/* internal interval : display Sint(X[J]-X[H],Y[I,K]-Y[G,D],J-H,I-K-G) + Sm(J,I,K), A[I,K,J].score */
  		if(DEBUG) assert(G>=0 && H>=0 && D>=0);
		TFLOAT x = X[J]-X[H];
		TFLOAT y = Yc(Y,I,K)-Yc(Y,G,D);
		TFLOAT Pen,Bias,PenSm;
		SintDetail(x,y,J-H,I-K-G,J,I,K,D,PRtabY,Y,c,Bias,Pen,PenSm,0);
		TFLOAT OutPen = OutlierPenalty[c];
		if(outlierBC)
		  OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-K-G];
		//		double sint = Sint(x,y,J-H,I-K-G,J,I,K,D,PRtabY,Y,c);
		TFLOAT sint = Bias + OutlierBias + (FIX_OUTLIER_MIS ? max(Pen,OutPen) + PenSm : max(Pen+PenSm,OutPen));
		double var = SF[c]*SF[c] + fabs(SD[c])*SD[c] * y;
		if(QUADRATIC_VARIANCE)
		  var += SR[c]*SR[c]*y*y;
		if(RES_VARIANCE){
		  FLOAT resR = Y[I] - Y[I-K];
		  FLOAT resL = Y[G] - Y[G-D];
		  FLOAT resvar = resL*resL + resR*resR;
		  var += SE[c]*SE[c]*resvar;
		}
		if(align2[c].orientation)
		  printf("G=%2d,H=%2d,D=%2d,X[J]=%0.3f,Y[I,K]=%0.3f,X[J]-X[H]=%0.3f,Y[I,K]-Y[G,D]=%0.3f,norm=%0.4f,Sint(X,Y,%d,%d)=%0.6f(Pen=%0.6f,Bias=%0.6f),A[I,K,J]=%0.6f,(X[M+1]-X[J])/C=%0.3f\n",
			 G,H,D,X[J],Yc(Y,I,K),x,y,(x-y)*(x-y)/var,J-H,I-K-G,sint,Pen,Bias,(bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J), (X[M+1]-X[J])*origPixelLen/PixelLen);
		else
		  printf("G=%2d,H=%2d,D=%2d,X[J]=%0.3f,Y[I,K]=%0.3f,X[J]-X[H]=%0.3f,Y[I,K]-Y[G,D]=%0.3f,norm=%0.4f,Sint(X,Y,%d,%d)=%0.6f(Pen=%0.6f,Bias=%0.6f),A[I,K,J]=%0.6f,X[J]/C=%0.3f\n",
			 G,H,D,X[J],Yc(Y,I,K),x,y,(x-y)*(x-y)/var,J-H,I-K-G,sint,Pen,Bias,(bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J), X[J]*origPixelLen/PixelLen);
		//		printf("G=%2d,H=%2d,D=%2d,X[J]=%0.3f,Y[I,K]=%0.3f,X[J]-X[H]=%0.3f,Y[I]-Y[G]=%0.3f,norm=%0.4f,Sint(X,Y,%d,%d)=%0.6f,A[I,K,J]=%0.6f\n",
		//		       G,H,D,X[J],Yc(Y,I,K),x,y,(x-y)*(x-y)/var,J-H,I-K-G,sint, (bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
	      }
	    }
	    //	    DP *p = &A[I][K][J];      
	    if(DEBUG && align2[c].Rend > -2) assert(J >= ((extend && I >= N+1-DELTA_Y) ? 1 : max(1,M+1-DELTA_X)));
	    int Rij = (align2[c].Rend <= -2) ? I : A->Rij(I,K,J);
	    if(DEBUG && extend && align2[c].Rend > -2) assert(I>=N+1-DELTA_Y || J >= M+1-DELTA_X);
	    int Rijx = (align2[c].Rend <= -2) ? J : extend ? A->Rijx(I,K,J) : M+1;

	    /* Right End: display Rij,Send(X[M+1]-X[J],min(M,Rijx)+1-J,min(N,Rij)+1-I) */
	    if(align2[c].Rend <= -2 && (NEW==1 || align2[c].Rend == -2))
	      printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rijk=%d,Rijkx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,C=%0.6f(%d),total score=%0.6f, X[M+1]/C=%0.3f\n",
		     I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), ChimScore, localtype, ((bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + ChimScore, X[M+1]*origPixelLen/PixelLen);
	    else if(align2[c].Rend <= -3)
	      printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rijk=%d,Rijkx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,EndPen=%0.6f,BiasEnd=%0.6f,total score=%0.6f, X[M+1]/C=%0.3f\n",
		     I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), OutlierEndPen, BiasEnd2(X[M+1]-X[J],M+1-J,c), ((bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + OutlierEndPen+BiasEnd2(X[M+1]-X[J],M+1-J,c), X[M+1]*origPixelLen/PixelLen);
	    else if(ENDFIX && X[M+1]-X[J] <= Y[N+1]-Yc(Y,I,K) && Rij <= N && Y[Rij+1]-Yc(Y,I,K) < X[M+1]-X[J]){
	      TFLOAT score = Sbnd(X[M+1]-X[J],Y[Rij+1]-Yc(Y,I,K),min(M,Rijx)+1-J,min(N,Rij)+1-I,c);
	      printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rij=%d,Rijx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Yc(Y,I,K)=%0.3f,Sbnd(X[M+1]-X[J],Y[Rij+1]-Yc(Y,I,K),M+1-J,Rij+1-I)=%0.6f,total score=%0.6f,X[M+1]/C=%0.3f\n",
		     I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), score, ((bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + max(ChimScore, score), X[M+1]*origPixelLen/PixelLen);
	    } else if(ENDFIX && X[M+1]-X[J] >= Y[N+1]-Yc(Y,I,K) && Rijx <= M && X[A->Rijx(I,K,J)+1]-X[J] < Y[N+1]-Yc(Y,I,K)){
	      TFLOAT score = Sbnd(X[Rijx+1]-X[J],Y[N+1]-Yc(Y,I,K),min(M,Rijx)+1-J,min(N,Rij)+1-I,c);
	      printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rij=%d,Rijx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,Sbnd(X[M+1]-X[J],Y[Rij+1]-Y[I,K],M+1-J,Rij+1-I)=%0.6f,total score=%0.6f,X[M+1]/C=%0.3f\n",
		     I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), score,  ((bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + max(ChimScore, score), X[M+1]*origPixelLen/PixelLen);
	    } else {
	      TFLOAT score = Send(min(X[M+1]-X[J],Y[N+1]-Yc(Y,I,K)),min(M,Rijx)+1-J,min(N,Rij)+1-I,c);
	      printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rij=%d,Rijx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,Send(min(X[M+1]-X[J],Y[N+1]-Y[I,K]),min(M,Rijx)+1-J,min(N,Rij)+1-I)=%0.6f,total score=%0.6f,X[M+1]/C=%0.3f\n",
		     I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), score, ((bestL2[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + max(ChimScore, (!extend && Rij > N) ? (TFLOAT)MINSCORE : score), X[M+1]*origPixelLen/PixelLen);
	    }
	    fflush(stdout);
	    if(DEBUG) assert(Rij >= I && Rij <= N+1);
	    /*      if(DEBUG>=2)
		    exit(1);*/
	  }
	} // pragma omp critical 
      } // if(VERB...)
    } // backtrack to retrive 2nd best alignment
  }

  if(SecondBest && bestscore > MINSCORE && (!align[-1].align2 || bestscore > align[-1].align2->score + 1e-8)){/* check if current best alignment should be saved as overall 2nd best alignment */
    int numpairs1 = align->numpairs;
    FLOAT alignYik = (numpairs1 <= 0) ? -1e100 : Yc(YY[0],align->sites1[numpairs1-1],align->sitesK1[numpairs1-1]);
    if(align[-1].orientation != orientation || fabs(alignYik - bestYik) >= MinSpacing){/* backtrack through array to retrieve best alignment in current orientation as overall 2nd best alignment */
      if(DEBUG) assert(bestscore <= align[-1].score + 1e-8);

      delete [] align[-1].align2; align[-1].align2 = new Calign[3];
      Calign *align2 = &align[-1].align2[1];
      align2[-1].mapid1 = refid;
      align2[-1].mapid2 = mapid;

      if(DEBUG) assert(bestI[0] >= 0 && bestI[1] >= 0);
      if(DEBUG){
	//	DP *p[2];
	//	for(int c = 0; c < colors; c++)
	//	  p[c] = &AA[c][bestI[c]][bestK[c]][bestJ[c]];
	double y0 = Yc(YY[0],bestI[0],bestK[0]);
	double y1 = Yc(YY[1],bestI[1],bestK[1]);
	double x0 = XX[0][bestJ[0]];
	double x1 = XX[1][bestJ[1]];
	double score = SMA(x1-x0, y1-y0, (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0), orientation);
	double Lscore = AA[0].scoreN(bestI[0],bestK[0],bestJ[0]) + AA[1].scoreN(bestI[1],bestK[1],bestJ[1]) + score;
	double LscoreEnd = AA[0].scoreEnd(bestI[0],bestK[0],bestJ[0]) + AA[1].scoreEnd(bestI[1],bestK[1],bestJ[1]) + score;
	double Rscore = AA[0].Rscore(bestI[0],bestK[0],bestJ[0]) + AA[1].Rscore(bestI[1],bestK[1],bestJ[1]);
	double RscoreEnd = AA[0].RscoreEnd(bestI[0],bestK[0],bestJ[0]) + AA[1].RscoreEnd(bestI[1],bestK[1],bestJ[1]);
	double bscore = ((bestR[0] >= -1) ? Rscore : RscoreEnd) + ((bestL[0] >= -1) ? Lscore : LscoreEnd);
	assert(fabs(bscore - bestscore) < 1e-8);
      }

      if(!Ilist){    /* allocate space to store alignment in reverse order */
	int M = max(MM[0],MM[1]);
	Ilist = new int[M];
	Klist = new int[M];
	Jlist = new int[M];
	iscore = new double[M+1];
	outscore = new double[M+1];
      }

      for(int c = 0; c < colors; c++){
	FLOAT *Y = YY[c];
	FLOAT *X = XX[c];
	int N = NN[c];
	int M = MM[c];
	AL2 *A = &AA[c];
	Cprtab **PRtabY = PRtab[c][refid];

	int I = bestI[c];
	int K = bestK[c];
	int J = bestJ[c];

	int U = 0;

	/* right end segment */
	// NOTE : if PoutlierEnd > 0.0, bestR[c] <= -2 may be caused either by local alignment or PoutlierEnd
	outscore[0] = iscore[0] = (bestR[c] >= -1 ? A->Rscore(I,K,J) : A->RscoreEnd(I,K,J));

	while(I > 0){
	  //	  DP *p = &A[I][K][J];
	  if(DEBUG){
	    assert(I-K >= 1 && I <= N);
	    assert(J >= 1 && J <= M);
	  }
	  Ilist[U] = I;
	  Klist[U] = K;
	  Jlist[U] = J;
	  if(bestL[c] >= -1){
	    int G = A->G(I,K,J);
	    int T = A->T(I,K,J);
	    int H = A->H(I,K,J);
	    if(G > 0){
	      if(DEBUG) assert(T >= 0 && H > 0);
	      iscore[U+1] = A->scoreN(I,K,J) - A->scoreN(G,T,H);
	      TFLOAT x = X[J] - X[H];
	      TFLOAT y = Yc(Y,I,K) - Yc(Y,G,T);
	      if(Poutlier > 0.0 && OUTLIER(maptype,x,y)){
		TFLOAT Bias,Pen,PenSm;
		SintDetail(x,y, J-H,I-K-G,J,I,K,T,PRtabY,Y,c,Bias,Pen,PenSm,0);
		outscore[U+1] = Bias + Pen + PenSm;
	      } else
		outscore[U+1] = iscore[U+1];
	    } else {/* left end segment */
	      // NOTE : if PoutlierEnd > 0.0, G <= -2 may be caused either by local alignment or PoutlierEnd
	      outscore[U+1] = iscore[U+1] = A->scoreN(I,K,J);
	    }
	    I = G;
	    K = T;
	    J = H;
	  } else {
	    int GEnd = A->GEnd(I,K,J);
	    int TEnd = A->TEnd(I,K,J);
	    int HEnd = A->HEnd(I,K,J);
	    if(GEnd > 0){
	      if(DEBUG) assert(TEnd >= 0 && HEnd > 0);
	      iscore[U+1] = A->scoreEnd(I,K,J) - A->scoreEnd(GEnd,TEnd,HEnd);
	      TFLOAT x = X[J] - X[HEnd];
	      TFLOAT y = Yc(Y,I,K) - Yc(Y,GEnd,TEnd);
	      if(Poutlier > 0.0 && OUTLIER(maptype,x,y)){
		TFLOAT Bias,Pen,PenSm;
		SintDetail(x,y, J-HEnd,I-K-GEnd,J,I,K,TEnd,PRtabY,Y,c,Bias,Pen,PenSm,0);
		outscore[U+1] = Bias + Pen + PenSm;
	      } else
		outscore[U+1] = iscore[U+1];
	    } else {/* left end segment */
	      // NOTE : if PoutlierEnd > 0.0, GEnd <= -2 may be caused either by local alignment or PoutlierEnd
	      outscore[U+1] = iscore[U+1] = A->scoreEnd(I,K,J);
	    }
	    I = GEnd;
	    K = TEnd;
	    J = HEnd;
	  }
	  U++;
	}

	if(DEBUG) assert(U > 0);
	if(DEBUG) assert(I <= (isLinear ? -1 : 0));
	if(DEBUG) assert((bestL[c] >= -1) ? I >= -1 : I <= -2);

	//	DP *p = &A[bestI[c]][bestK[c]][bestJ[c]];
	int bI = bestI[c], bK = bestK[c], bJ = bestJ[c];
	align2[c].score = ((bestL[c] >= -1) ? A->scoreN(bI,bK,bJ) : A->scoreEnd(bI,bK,bJ)) + ((bestR[c] >= -1) ? A->Rscore(bI,bK,bJ) : A->RscoreEnd(bI,bK,bJ));
	align2[c].orientation = orientation;
	align2[c].scaleID = scaleID;
	align2[c].logPV = alignFPsd(A,bI,bJ,bK,Y,X,N,M,bestR[c],bestL[c],orientation,refid,mapid,c,0);
	align2[c].Lend = I;
	align2[c].Rend = bestR[c];

	I = Ilist[U-1];
	K = Klist[U-1];
	J = Jlist[U-1];

	if(DEBUG>=2) assert(align2[c].Lend <= -2 || J <= DELTA_X || (extend && I-K <= DELTA_Y));
	align2[c].Lij1 = (align2[c].Lend <= -2) ? I-K : A->Lij(I,K,J);
	if(DEBUG) assert(align2[c].Lij1 <= I-K && align2[c].Lij1 >= 0);
	align2[c].Lij2 = (align2[c].Lend <= -2) ? J : extend ? A->Lijx(I,K,J) : 0;
	if(DEBUG) assert(align2[c].Lij2 <= J && align2[c].Lij2 >= 0);
	if(DEBUG>=2) assert(align2[c].Rend <= -2 || bestJ[c] >= M+1 - DELTA_X || (extend && bestI[c] >= N+1-DELTA_Y));
	align2[c].Rij1 = (align2[c].Rend <= -2) ? bestI[c] : A->Rij(bI,bK,bJ);
	if(DEBUG) assert(align2[c].Rij1 >= bestI[c] && align2[c].Rij1 <= N+1);
	align2[c].Rij2 = (align2[c].Rend <= -2) ? bestJ[c] : extend ? A->Rijx(bI,bK,bJ) : M+1;
	if(DEBUG)	assert(align2[c].Rij2 >= bestJ[c] && align2[c].Rij2 <= M+1);

	/* reallocate best alignment */
	align2[c].allfree();
	if(DEBUG) assert(U > 0);
	align2[c].expand_arrays(U);
	align2[c].numpairs = U;

	/* copy alignment from Ilist[],Jlist[] in reverse order */
	for(int I=0;I<U;I++){
	  align2[c].sites1[I] = Ilist[U-1-I];
	  align2[c].sitesK1[I] = Klist[U-1-I];
	  align2[c].sites2[I] = Jlist[U-1-I];
	  align2[c].iscore[I] = iscore[U-I];
	  align2[c].outscore[I] = outscore[U-I];
	}
	align2[c].iscore[U] = iscore[0];
	align2[c].outscore[U] = outscore[0];

	align2[c].maxoutlier = 0.0;
	align2[c].noutliers = 0;
	int lastI = align2[c].sites1[0];
	int lastK = align2[c].sitesK1[0];
	int lastJ = align2[c].sites2[0];
	for(int t = 1; t < U; lastI = I, lastK = K, lastJ = J, t++){
	  int I = align2[c].sites1[t];
	  int K = align2[c].sitesK1[t];
	  int J = align2[c].sites2[t];
	  if(align2[c].outscore[t] + (RFLOAT)0.01 < align2[c].iscore[t]){
	    align2[c].noutliers++;
	    double deltaY = Yc(Y,I,K) - Yc(Y,lastI,lastK);
	    double deltaX = X[J] - X[lastJ];
	    double delta = fabs(deltaY-deltaX);
	    align2[c].maxoutlier = max(delta,align2[c].maxoutlier);
	  }
	}
      }// c = 0 .. colors -1

      /* compute combined score(s) : only score2,orientation2,logPV2,numpairs2 are valid for align[-1] */
      align2[-1].score = bestscore2;
      align2[-1].orientation = orientation;
      align2[-1].scaleID = scaleID;
      align2[-1].logPV = align2[0].logPV + align2[1].logPV;/* excludes misalignment information */
      if(VERB>=2 && refmap[refid]->id == 1 && Gmap[mapid]->id == 101000333LL && !orientation){
	printf("C:refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:align=%p(logPV=%0.2f:%0.2f,%0.2f),align2=%p(logPV=%0.2f:%0.2f,%0.2f)\n",
	       refid,refmap[refid]->id,mapid,Gmap[mapid]->id,orientation,
	       align, align[-1].logPV, align[0].logPV, align[1].logPV,
	       align2,align2[-1].logPV,align2[0].logPV,align2[1].logPV);
	fflush(stdout);
      }
      align2[-1].numpairs = align2[0].numpairs + align2[1].numpairs;

      if(VERB>=2 && bestscore > MINSCORE && (VERB>=2 || (giter == RefRepeats-1/* && Gmap[mapid]->id == 35008162*/ /* && bestscore > ScoreThreshold && align[-1].logPV > LogPvThreshold && align[-1].numpairs >= AlignedSiteThreshold && AlignedEndOutlier(align) */)/* ||  (mapid==1254 && align->orientation) */ /* || fabs(offset-trueoffset) > 2*Ylambda || align->orientation*/)){
        #pragma omp critical 
	{ 
	  double offset = bestYik2 - XX[0][bestJ[0]];
	  double trueoffset = Gmap[mapid]->endloc - XX[0][MM[0]+1];
	  double delta = offset-trueoffset;
	  Cmap *pmap = Gmap[mapid];
	  while(pmap->origmap)
	    pmap = pmap->origmap;

	  printf("2nd Best(2):refid=%d,mapid=%d(id=%lld),score=%0.4f(c1=%0.4f,I=%d,K=%d,J=%d,R=%d,c2=%0.4f,P=%d,T=%d,Q=%d,S=%d,ma=%0.4f),logPV=%0.4f(%0.4f,%0.4f):orientation=%d,scaleid=%d(%0.4f),trueflip=%d,trueoffset=%0.3f,delta=%0.3f,aligned sites=%d+%d, Lend=%d,%d,Rend=%d,%d,fpcnt=%d(%d):%s\n",
		 refid,mapid,Gmap[mapid]->id,bestscore,align2[0].score,bestI[0],bestK[0],bestJ[0],bestR[0],align2[1].score,bestI[1],bestK[1],bestJ[1],bestR[1],bestscore2 - align2[0].score - align2[1].score,
		 align2[-1].logPV,align2[0].logPV,align2[1].logPV,orientation,scaleID, scaleID ? ScaleFactor[scaleID] : 1.0,
		 Gmap[mapid]->flipped,trueoffset,delta,align2[0].numpairs,align2[1].numpairs,align2[0].Lend,align2[1].Lend,align2[0].Rend,align2[1].Rend,
		 MM[0]+MM[1] - align2[-1].numpairs, Gmap[mapid]->fpcnt,pmap->name ? pmap->name : "");

	  for(int c = 0; c < colors; c++){
	    FLOAT *Y = YY[c];
	    FLOAT *X = XX[c];
	    int N = NN[c];
	    int M = MM[c];
	    AL2 *A = &AA[c];
	    Cprtab **PRtabY = PRtab[c][refid];

	    printf("color=%d:\n",c+1);
	    fflush(stdout);

	    (void)alignFPsd(A,bestI2[c],bestJ2[c],bestK2[c],Y,X,N,M,bestR2[c],bestL2[c],orientation,refid,mapid,c,1);
	    int I= -1,J= -1 ,K=0, G= -1,H= -1, D = 0;
	    int U = align2[c].numpairs;

	    for(int T=0; T < U; T++, G=I, H=J, D=K){
	      I = align2[c].sites1[T];
	      K = align2[c].sitesK1[T];
	      J = align2[c].sites2[T];
	      printf("T=%2d:I=%2d,K=%2d,J=%2d:",T,I,K,J);
	      if(T<=0){ /* Left End: display Send(min(X[J],Y[I-K..I]),J+1-max(1,Lijx),I+1-max(1,Lij))+Sm(J,I,K), A[I,K,J].score */
		if(DEBUG && align2[c].Lend > -2) assert(J <= DELTA_X || (extend && I-K <= DELTA_Y));
		int Lij = (align2[c].Lend <= -2) ? I-K : A->Lij(I,K,J);
		if(DEBUG && align2[c].Lend > -2 && extend) assert(I-K <= DELTA_Y || J <= DELTA_X);
		int Lijx = (align2[c].Lend <= -2) ? J : extend ? A->Lijx(I,K,J) : 0;

		if(align2[c].Lend <= -2 && (NEW==1 || align2[c].Lend == -2))
		  printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I,K]=%6.3f,C=%0.6f(%d),Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
			 Lij,Lijx,X[J],Yc(Y,I,K),ChimScore,localtype,Sm(J,I,K,Y,c),ChimScore+Sm(J,I,K,Y,c), (bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		else if(align2[c].Lend <= -3)
		  printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I,K]=%6.3f,EndPen=%0.6f,BiasEnd=%0.6f,Sm(J,I,K)=%0.6f,score=%0.6f,A[I,K,J]=%0.6f\n",
			 Lij,Lijx,X[J],Yc(Y,I,K),OutlierEndPen,BiasEnd2(X[J],J,c),Sm(J,I,K,Y,c),BiasEnd2(X[J],J,c)+OutlierEndPen+Sm(J,I,K,Y,c),(bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		else if(ENDFIX && X[J] <= Yc(Y,I,K) && Lij > 0 && Yc(Y,I,K) - Y[Lij-1] < X[J]){
		  TFLOAT score = Sbnd(X[J],Yc(Y,I,K)-Y[Lij-1],J,I-K+1-Lij,c);
		  printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I]=%6.3f,Sbnd(X[J],Y[I,K]-Y[Lij-1],J,I-K+1-Lij)=%0.6f,Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
			 Lij,Lijx,X[J],Yc(Y,I,K),score,Sm(J,I,K,Y,c), max(ChimScore,score+(TFLOAT)Sm(J,I,K,Y,c)), (bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		  if(DEBUG) assert(Lijx <= 0);
		} else if(ENDFIX && extend && X[J] >= Yc(Y,I,K) && Lijx > 0 && X[J]-X[Lijx-1] < Yc(Y,I,K)){
		  TFLOAT score = Sbnd(X[J]-X[Lijx-1],Yc(Y,I,K),J+1-Lijx,I-K,c);
		  printf("Lijy=%d,Lijx=%d,X[J]=%0.6f,Y[I,K]=%6.3f,Sbnd(%0.3f,%0.3f,%d,%d)=%0.6f,Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
			 Lij,Lijx,X[J],Y[I],X[J]-X[Lijx-1],Yc(Y,I,K),J+1-Lijx,I-K, score,Sm(J,I,K,Y,c), max(ChimScore,score+(TFLOAT)Sm(J,I,K,Y,c)), (bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		  if(DEBUG) assert(Lij <= 0);
		} else {
		  TFLOAT score = Send(X[J],J+1-max(1,Lijx) ,I-K+1 - max(1,Lij),c);
		  printf("Lij=%d,Lijx=%d,X[J]=%6.3f,Y[I,K]=%6.3f,Send(min(X[J],Y[I,K]),J+1-max(1,Lijx),I-K+1-max(1,Lij))=%0.6f,Sm(J,I,K)=%0.4f,score=%0.6f,A[I,K,J]=%0.6f\n",
			 Lij,Lijx,X[J],Yc(Y,I,K),score,Sm(J,I,K,Y,c), max(ChimScore,Lij > 0 ? score+(TFLOAT)Sm(J,I,K,Y,c) : (TFLOAT)MINSCORE), (bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
		}
	      } else {/* internal interval : display Sint(X[J]-X[H],Y[I,K]-Y[G,D],J-H,I-K-G) + Sm(J,I,K), A[I,K,J].score */
		if(DEBUG) assert(G>=0 && H>=0 && D>=0);
		TFLOAT x = X[J]-X[H];
		TFLOAT y = Yc(Y,I,K)-Yc(Y,G,D);
		TFLOAT Pen,Bias,PenSm;
		SintDetail(x,y,J-H,I-K-G,J,I,K,D,PRtabY,Y,c,Bias,Pen,PenSm,0);
		TFLOAT OutPen = OutlierPenalty[c];
		if(outlierBC)
		  OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-K-G];
		//		double sint = Sint(x,y,J-H,I-K-G,J,I,K,D,PRtabY,Y,c);
		TFLOAT sint = Bias + OutlierBias + (FIX_OUTLIER_MIS ? max(Pen,OutPen) + PenSm : max(Pen+PenSm,OutPen));
		double var = SF[c]*SF[c] + fabs(SD[c])*SD[c] * y;
		if(QUADRATIC_VARIANCE)
		  var += SR[c]*SR[c]*y*y;
		if(RES_VARIANCE){
		  FLOAT resR = Y[I] - Y[I-K];
		  FLOAT resL = Y[G] - Y[G-D];
		  FLOAT resvar = resL*resL + resR*resR;
		  var += SE[c]*SE[c]*resvar;
		}
		if(align2[c].orientation)
		  printf("G=%2d,H=%2d,D=%2d,X[J]=%0.3f,Y[I,K]=%0.3f,X[J]-X[H]=%0.3f,Y[I,K]-Y[G,D]=%0.3f,norm=%0.4f,Sint(X,Y,%d,%d)=%0.6f(Pen=%0.6f,Bias=%0.6f,Sm=%0.6f),A[I,K,J]=%0.6f,(X[M+1]-X[J])/C=%0.3f\n",
			 G,H,D,X[J],Yc(Y,I,K),x,y,(x-y)*(x-y)/var,J-H,I-K-G,sint,Pen,Bias,PenSm, (bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J), (X[M+1]-X[J])*origPixelLen/PixelLen);
		else
		  printf("G=%2d,H=%2d,D=%2d,X[J]=%0.3f,Y[I,K]=%0.3f,X[J]-X[H]=%0.3f,Y[I,K]-Y[G,D]=%0.3f,norm=%0.4f,Sint(X,Y,%d,%d)=%0.6f(Pen=%0.6f,Bias=%0.6f,Sm=%0.6f),A[I,K,J]=%0.6f,X[J]/C=%0.3f\n",
			 G,H,D,X[J],Yc(Y,I,K),x,y,(x-y)*(x-y)/var,J-H,I-K-G,sint,Pen,Bias,PenSm, (bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J), X[J]*origPixelLen/PixelLen);
		//		printf("G=%2d,H=%2d,D=%2d,X[J]=%0.3f,Y[I,K]=%0.3f,X[J]-X[H]=%0.3f,Y[I]-Y[G]=%0.3f,norm=%0.4f,Sint(X,Y,%d,%d)=%0.6f,A[I,K,J]=%0.6f\n",
		//		       G,H,D,X[J],Yc(Y,I,K),x,y,(x-y)*(x-y)/var,J-H,I-K-G,sint, (bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J));
	      }
	    }
	    //	    DP *p = &A[I][K][J];      
	    if(DEBUG && align2[c].Rend > -2) assert(J >= ((extend && I >= N+1-DELTA_Y) ? 1 : max(1,M+1-DELTA_X)));
	    int Rij = (align2[c].Rend <= -2) ? I : A->Rij(I,K,J);
	    if(DEBUG && extend && align2[c].Rend > -2) assert(I>=N+1-DELTA_Y || J >= M+1-DELTA_X);
	    int Rijx = (align2[c].Rend <= -2) ? J : extend ? A->Rijx(I,K,J) : M+1;

	    /* Right End: display Rij,Send(X[M+1]-X[J],min(M,Rijx)+1-J,min(N,Rij)+1-I) */
	    if(align2[c].Rend <= -2 && (NEW==1 || align2[c].Rend == -2))
	      printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rijk=%d,Rijkx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,C=%0.6f(%d),total score=%0.6f, X[M+1]/C=%0.3f\n",
		     I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), ChimScore, localtype, ((bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + ChimScore, X[M+1]*origPixelLen/PixelLen);
	    else if(align2[c].Rend <= -3)
	      printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rijk=%d,Rijkx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,EndPen=%0.6f,BiasEnd=%0.6f,total score=%0.6f, X[M+1]/C=%0.3f\n",
		     I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), OutlierEndPen, BiasEnd2(X[M+1]-X[J],M+1-J,c), ((bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + OutlierEndPen+BiasEnd2(X[M+1]-X[J],M+1-J,c), X[M+1]*origPixelLen/PixelLen);
	    else if(ENDFIX && X[M+1]-X[J] <= Y[N+1]-Yc(Y,I,K) && Rij <= N && Y[Rij+1]-Yc(Y,I,K) < X[M+1]-X[J]){
	      TFLOAT score = Sbnd(X[M+1]-X[J],Y[Rij+1]-Yc(Y,I,K),min(M,Rijx)+1-J,min(N,Rij)+1-I,c);
	      printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rij=%d,Rijx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Yc(Y,I,K)=%0.3f,Sbnd(X[M+1]-X[J],Y[Rij+1]-Yc(Y,I,K),M+1-J,Rij+1-I)=%0.6f,total score=%0.6f,X[M+1]/C=%0.3f\n",
		     I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), score, ((bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + max(ChimScore, score), X[M+1]*origPixelLen/PixelLen);
	    } else if(ENDFIX && X[M+1]-X[J] >= Y[N+1]-Yc(Y,I,K) && Rijx <= M && X[A->Rijx(I,K,J)+1]-X[J] < Y[N+1]-Yc(Y,I,K)){
	      TFLOAT score = Sbnd(X[Rijx+1]-X[J],Y[N+1]-Yc(Y,I,K),min(M,Rijx)+1-J,min(N,Rij)+1-I,c);
	      printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rij=%d,Rijx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,Sbnd(X[M+1]-X[J],Y[Rij+1]-Y[I,K],M+1-J,Rij+1-I)=%0.6f,total score=%0.6f,X[M+1]/C=%0.3f\n",
		     I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), score,  ((bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + max(ChimScore, score), X[M+1]*origPixelLen/PixelLen);
	    } else {
	      TFLOAT score = Send(min(X[M+1]-X[J],Y[N+1]-Yc(Y,I,K)),min(M,Rijx)+1-J,min(N,Rij)+1-I,c);
	      printf("I=%d,K=%d,J=%d,N=%d,M=%d:Rij=%d,Rijx=%d,X[M+1]=%0.3f,Y[N+1]=%0.3f,X[M+1]-X[J]=%0.3f,Y[N+1]-Y[I,K]=%0.3f,Send(min(X[M+1]-X[J],Y[N+1]-Y[I,K]),min(M,Rijx)+1-J,min(N,Rij)+1-I)=%0.6f,total score=%0.6f,X[M+1]/C=%0.3f\n",
		     I,K,J,N,M,Rij,Rijx,X[M+1],Y[N+1],X[M+1]-X[J],Y[N+1]-Yc(Y,I,K), score, ((bestL[c] >= -1) ? A->scoreN(I,K,J) : A->scoreEnd(I,K,J)) + max(ChimScore, (!extend && Rij > N) ? (TFLOAT)MINSCORE : score), X[M+1]*origPixelLen/PixelLen);
	    }
	    fflush(stdout);
	    if(DEBUG) assert(Rij >= I && Rij <= N+1);
	    /*      if(DEBUG>=2)
		    exit(1);*/
	  }
	} // pragma omp critical 
      } // if(VERB...)
    } // backtrack through array to retrieve current best alignment as overall 2nd best alignment 
  }

  if(Ilist){
    delete [] Ilist;
    delete [] Klist;
    delete [] Jlist;
    delete [] iscore;
    delete [] outscore;
  }
}

/* version of refalign() that can handle resSD > 0.0 using 3-D recurrance based alignment */
static void refalignSD(Cmap *rmap,
		       Cmap *nanomap,
		       size_t &alignid, /**< reference to index into alignment[] array. NOTE : alignment[alignid] may be NULL.
					If phash && NoSplit==1 && HashMultiMatch up to 3 alignment results may be returned in alignment[alignid..alignid+2] */
		       //		       Calign *align,/* each color's best alignment in align[0..colors-1] : combined score in align[-1] */
		       int **Kmax, /* Kmax[0..colors-1] */
		       int refid,int mapid, /* id's used to identify rmap & nanomap in result (align) */
		       CHashMatch *phash, /**< If !=0 pointer to hashtable entry */
		       CXPen *XPen, /* per thread score table */
		       int *Qmin, int *Qmax,/* preallocated memory Qmin[1..maxM],Qmax[1..maxM] */
		       AL2 *A, /* A[c=0..color-1]].array[I=1..N][K=0..Kmax[I]][(J=1..maxM)] : only 1st 3 dimensions (c,I,K) are allocated */
#if FLAT_MEMORY
		       TFLOAT **Fmem, /**< preallocated memory Fmem[c=0..1][0 ... pcnt[c]*maxM[c]*DPsizF - 1] to be used by A[c].field[I][K][] of type TFLOAT */
		       int **Imem, /**< preallocated memory Imem[c=0..1][0 ... pcnt[c]*maxM[c]*DPsizI - 1] to be used by A->field[I][K][] of type int */
		       long long *StrideMax,/* If (MIN_MEM && DIAGONAL>=2 && hashdelta) : the current allocated number of elements per array */
		       int **JMIN, int **JMAX, /**< preallocated memory Jmin[c=0..1][1..maxN[c]], Jmax[c=0..1][1..maxN[c]] used with phash offset to create diagonalized arrays */
		       int **Imin, int **Imax,  /**< preallocated memory Imin[c][1..maxM[c]], Imax[c][1..maxM[c]] used with phash offset to create diagonalized arrays */
		       int tid,
		       long long *limN,/* use instead of maxN : limN[c] == (MINMEM ? NN[c] : maxN[c]) */
#else
		       DP **Amem, /* Amem[c=0..color-1] is preallocated memory to be used by A[c].array[I][K][] */
#endif
		       int nosplit /* local value used in place of global NoSplit */
		       )
{
  Calign* &galign = alignment[alignid++];

  if(DEBUG)assert(resSD[0] > 0.0 && resSD[1] > 0.0);

  // if(DEBUG && galign) assert(galign->numpairs == 0);  
  if(galign){
    for(int c = 0; c <= colors; c++)
      galign[c].allfree();
    if(ALIGN_COMPRESS){
      delete [] galign;
      galign = NULL;
    }
  }

  Calign *align, Lalign[MAXCOLOR+1];
  if(galign != 0)
    align = galign;
  else
    align = &Lalign[0];/* a local temporary variable : alignment will be saved to galign IFF score is good */
  align++;

  for(int c = -1; c < colors; c++){
    align[c].mapid1 = refid;
    align[c].mapid2 = mapid;
    align[c].Lend = align[c].Rend = -1;
    align[c].orientation = -1;
    align[c].sites1 = align[c].sitesK1 = align[c].sites2 = NULL;
    align[c].score = MINSCORE;
    if(RepeatMaxShift)
      align[c].repeat = 0;
  }
  if(SecondBest){
    delete [] align[-1].align2;
    align[-1].align2 = NULL;
  }

  int *NN = rmap->numsite;
  FLOAT **YY = rmap->site;
  int *MM = nanomap->numsite;
  FLOAT **XX = nanomap->site;

  /* Do memory allocation for dynamic programming array A */
#if FLAT_MEMORY
  /* NOTE : if (DIAGONAL>=2 && phash && hashdelta) A->field[][][] will be re-allocated in refalignXYsd() using malloc()/free() */
  if(!(DIAGONAL >= 2 && phash && hashdelta)){  /* re-allocate memory for A->field[][][] */
    for(int c = 0; c < colors; c++){
      long long mcnt = 0;
      A[c].Kstride = limN[c] * maxM[c];
      A[c].Istride = MM[c];
      long long stride = (1LL + A[c].kmax) * A[c].Kstride;

#define FALLOC(field) { A[c].field##_base = &Fmem[c][mcnt - A[c].Istride - 1]; mcnt += stride; }
      FALLOC(scoreN);
      FALLOC(scoreEnd);
      FALLOC(Lscore);
      FALLOC(LscoreEnd);
      FALLOC(Rscore);
      FALLOC(RscoreEnd);

#undef FALLOC
      if(DEBUG) assert(mcnt <= (1LL + A[c].kmax) * limN[c] * maxM[c] * DPsizF);

      mcnt = 0;
#define IALLOC(field)  A[c].field##_base = &Imem[c][mcnt - A[c].Istride - 1]; mcnt += stride;
      IALLOC(G);
      IALLOC(T);
      IALLOC(H);
      IALLOC(GEnd);
      IALLOC(TEnd);
      IALLOC(HEnd);
      IALLOC(Lij);
      IALLOC(Rij);
      IALLOC(Lijx);
      IALLOC(Rijx);
      IALLOC(L);
      IALLOC(LEnd);
      IALLOC(R);
      IALLOC(REnd);
#undef IALLOC
      if(DEBUG) assert(mcnt <= (1LL + A[c].kmax) * limN[c] * maxM[c] * DPsizI);
    }
  }
#else
  for(int c = 0; c < colors; c++){
    long long mcnt = 0;
    for(int I=1; I <= NN[c]; I++){
      for(int K=0; K <= Kmax[c][I]; K++){
	A[c].array[I][K] = &Amem[c][mcnt-1];/* A[c].array[I][K][J=1..MM[c]] */
	mcnt += MM[c];
      }
    }
  }
#endif

  FLOAT *Xrev[MAXCOLOR];
  Xrev[0] = (FLOAT *)alloca((MM[0]+MM[1]+4)*sizeof(FLOAT));
  Xrev[1] = &Xrev[0][MM[0]+2];

  if(DEBUG) assert(nanomap->mapid == mapid);
  Cmap *origmap = nanomap;
  while(origmap->origmap)
    origmap = origmap->origmap;
  int origmapid = origmap->mapid;

  if(VERB>=2 && phash && hashdelta && nanomap->id == MAP_TRACE){
    CHashMatch *hash = phash;
    printf("refalignSD:refid=%d,mapid=%d(id=%lld),origmapid=%d,hashdelta=%0.1f:hashtable entries:\n",refid,mapid,Gmap[mapid]->id,origmapid,hashdelta);
    while(hash->id1 == refid && hash->id2 == origmapid){
      printf("phash=%p:mapid1=%d,mapid2=%d,hashscore=%d,offset=%d,orientation=%d\n",
	     hash,hash->id1,hash->id2,hash->hashscore,hash->offset,hash->orientation);
      hash++;
    }
    printf("[phash=%p:mapid1=%d,mapid2=%d,hashscore=%d,offset=%d,orientation=%d]\n",
	   hash,hash->id1,hash->id2,hash->hashscore,hash->offset,hash->orientation);
    fflush(stdout);
  }

  for(int scaleID = 0; scaleID < NumScaleFactor; scaleID++){
    FLOAT scale = ScaleFactor[scaleID];
    if(DEBUG && scaleID == 0 && NumScaleFactor > 0) assert(scale == 1.0);

    CHashMatch *origphash = phash;

    /* first check X in normal orientation and update best alignment */
    for(int cnt = 0; HashMultiMatchMax <= 0 || cnt < HashMultiMatchMax; cnt++){ /* if(phash && hashdelta) loop over multiple offsets */

      if(cnt && (cnt%1000)==0){
	printf("WARNING:refid=%d,mapid=%d(id=%lld),or=0:over %d hashtable offsets:phash=%p:id1=%d,id2=%d,hashscore=%d,offset=%d,orientation=%d\n",
	       refid,mapid,nanomap->id,cnt,phash,phash->id1,phash->id2,phash->hashscore,phash->offset,phash->orientation);
	fflush(stdout);
      }

      if(VERB>=2 && phash){
	printf("cnt1=%d,scaleID=%d:phash=%p:id1=%d,id2=%d,hashscore=%d,offset=%d,orientation=%d\n",
	       cnt,scaleID,phash,phash->id1,phash->id2,phash->hashscore,phash->offset,phash->orientation);
	fflush(stdout);
      }

      if(phash && (phash->id1 != refid || phash->id2 != origmapid || phash->orientation))
	break;

      if(scaleID > 0){
	for(int c = 0; c < colors; c++)
	  for(register int J=0; J <= MM[c]+1;J++)
	    Xrev[c][J] = XX[c][J]*scale;
	refalignXYsd(YY, NN, Xrev, MM, 0, scaleID, Kmax, A, 
#if FLAT_MEMORY
		     Fmem, Imem, StrideMax, JMIN, JMAX, Imin, Imax,tid, limN,
#endif
		     align, rmap, nanomap, refid, mapid, XPen, phash, Qmin, Qmax);
      } else
	refalignXYsd(YY, NN, XX, MM, 0, scaleID, Kmax, A, 
#if FLAT_MEMORY
		     Fmem, Imem, StrideMax, JMIN, JMAX, Imin, Imax,tid, limN,
#endif
		     align, rmap, nanomap, refid, mapid, XPen, phash, Qmin, Qmax);
  
      if(!(phash && hashdelta))
	break;
      phash++;

      if(!hashdelta)
	while(phash->id1 == refid && phash->id2 == origmapid && !phash->orientation)
	  phash++;
      if(DEBUG && !(phash->id1 != refid || phash->id2 != origmapid || phash->orientation)) assert(HashMultiMatch);
    }

    /* next check X in reversed orientation and update best alignment */
    for(int cnt = 0; HashMultiMatchMax <= 0 || cnt < HashMultiMatchMax ; cnt++){ /* if(phash && hashdelta) loop over multiple offsets */
      if(cnt && (cnt%1000)==0){
	printf("WARNING:refid=%d,mapid=%d(id=%lld),or=1:over %d hashtable offsets:phash=%p:id1=%d,id2=%d,hashscore=%d,offset=%d,orientation=%d\n",
	       refid,mapid,nanomap->id,cnt,phash,phash->id1,phash->id2,phash->hashscore,phash->offset,phash->orientation);
	fflush(stdout);
      }

      if(VERB>=2 && phash){
	printf("cnt2=%d,scaleID=%d:phash=%p:id1=%d,id2=%d,hashscore=%d,offset=%d,orientation=%d\n",
	       cnt,scaleID,phash,phash->id1,phash->id2,phash->hashscore,phash->offset,phash->orientation);
	fflush(stdout);
      }

      if(phash){
	if(phash->id1 != refid || phash->id2 != origmapid)
	  break;
	if(!phash->orientation){
	  if(DEBUG && HashMultiMatchMax <= 0) assert(phash->orientation);
	  continue;/* skip excess matches with orientation == 0 */
	}
      }

      for(int c = 0; c < colors; c++)
	for(int J=0;J <= MM[c]+1;J++)
	  Xrev[c][J] = (XX[c][MM[c] + 1] - XX[c][MM[c] + 1 - J]) * scale;
      refalignXYsd(YY, NN, Xrev, MM, 1, scaleID, Kmax, A, 
#if FLAT_MEMORY
		   Fmem, Imem, StrideMax, JMIN, JMAX, Imin, Imax,tid,limN,
#endif
		   align, rmap, nanomap, refid, mapid, XPen, phash, Qmin, Qmax);

      if(!(phash && hashdelta))
	break;
      phash++;

      if(!hashdelta)
	while(phash->id1 == refid && phash->id2 == origmapid && !phash->orientation)
	  phash++;
      if(DEBUG && !(phash->id1 != refid || phash->id2 != origmapid)) assert(HashMultiMatch);
    }

    phash = origphash;

    if(0 && AlignedThreshold(&align[-1], YY, ScoreThreshold, LogPvThreshold))
      break;
  }

  align[-1].chimpair = 0;

  if(VERB>=2 && nanomap->id == MAP_TRACE){
    #pragma omp critical
    {
      printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d,,score=%0.3f/%0.3f,logPV=%0.3f/%0.3f,numpairs=%d/%d,Length=%0.3f/%0.3f,repeat=%d,IsRepeatRegion=%d\n",
	     refid,rmap->id, mapid,nanomap->id, align[-1].orientation,align[-1].score,ScoreThreshold,align[-1].logPV,LogPvThreshold,align[-1].numpairs,AlignedSiteThreshold, 
	     align[-1].numpairs ? AlignedLength(align,YY) : 0.0, AlignedLengthThreshold, align->repeat, 0 /* align->IsRepeatRegion() ? 1 : 0*/);
      fflush(stdout);
    }
  }

  if(AlignedThreshold(&align[-1], YY, ScoreThreshold, LogPvThreshold)){
    
    if(VERB>=2 && nanomap->id == MAP_TRACE){
      #pragma omp critical
      {
	printf("alignment[%lu]:refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:score=%0.4f,logPV=%0.3f,pairs=%d(%d,%d),nosplit=%d,repeat=%d,IsRepeatRegion=%d,Lend=%d,%d,Rend=%d,%d\n",
	       alignid-1,align->mapid1,rmap->id,align->mapid2,Gmap[align->mapid2]->id,align[-1].orientation,align[-1].score,align[-1].logPV,align[-1].numpairs,align[0].numpairs,align[1].numpairs,
	       nosplit,align[-1].repeat,0 /*align->IsRepeatRegion()*/, align[0].Lend,align[1].Lend,align[0].Rend,align[1].Rend);
	if(SecondBest && align[-1].align2){
	  Calign *align2 = align[-1].align2;
	  printf("    2nd Best:refid=%d,mapid=%d(id=%lld),or=%d:score=%0.4f,logPV=%0.3f,pairs=%d\n",
		 align2->mapid1,align2->mapid2,Gmap[align2->mapid2]->id,align2->orientation,align2->score,align2->logPV,align2->numpairs);
	}
	fflush(stdout);
      }
    }

    if(!galign){/* save local Lalign to galign */
      align--;
      if(DEBUG) assert(align == &Lalign[0]);

      galign = new Calign[1+colors];

      for(int c = 0; c <= colors; c++){
	/* swap galign[c] and align[c] : align[c] is on the stack and will automatically be deallocated */
	Calign tmp = galign[c];
	galign[c] = align[c];
	align[c] = tmp;
	if(DEBUG && c > 0 && galign[c].numpairs > 0){
	  if(!(galign[c].sites1[0] >= 1 && galign[c].sites1[galign[c].numpairs-1] <= NN[c-1])){
	    printf("refid=%d,mapid=%d:c=%d,galign[c].numpairs=%d,NN[c]=%d\n",refid,mapid,c,galign[c].numpairs,NN[c-1]);
	    printf("  galign[c].sites1[0]=%d\n",galign[c].sites1[0]);
	    printf("  galign[c].sites1[%d]=%d\n",galign[c].numpairs-1,galign[c].sites1[galign[c].numpairs-1]);
	    fflush(stdout);
	    assert(galign[c].sites1[0] >= 1 && galign[c].sites1[galign[c].numpairs-1] <= NN[c-1]);
	  }
	  if(!(galign[c].sites2[0] >= 1 && galign[c].sites2[galign[c].numpairs-1] <= MM[c-1])){
	    printf("refid=%d,mapid=%d:c=%d,galign[c].numpairs=%d,MM[c]=%d\n",refid,mapid,c,galign[c].numpairs,MM[c-1]);
	    printf("  galign[c].sites2[0]=%d\n",galign[c].sites2[0]);
	    printf("  galign[c].sites2[%d]=%d\n",galign[c].numpairs-1,galign[c].sites2[galign[c].numpairs-1]);
	    fflush(stdout);
	    assert(galign[c].sites2[0] >= 1 && galign[c].sites2[galign[c].numpairs-1] <= MM[c-1]);
	  }
	}
      }
      if(DEBUG) assert(galign->numpairs > 0);
      align = galign;

      align++;
    }

    align--;

    /* check if mapid is a chimeric fragment : if so check how close the original fragment aligned (tandem alignment) */
    if(nanomap->origmap){
      /* locate original alignment in alignment[0..numaligns-2] */
      Cmap *origmap = Gmap[mapid]->origmap;
      Calign *origalign = 0;
      int i;
      for(i = numaligns-1;--i >= 0;)
	if((origalign = alignment[i]) && Gmap[origalign->mapid2] ==  origmap && origalign->mapid1 == align->mapid1)
	  break;
      if(i < 0){
	fprintf(stderr,"refid=%d,mapid=%d:Could not find alignment for original fragment\n",refid,mapid);
	exit(1);
      }
      if(nosplit && DEBUG) assert(!Gmap[origalign->mapid2]->origmap);

      int tandem = 0;
      if(align->orientation == origalign->orientation){/* if orientations don't match, the fragments cannot be tandem alignments */
	if(!align->orientation){
	  int left[2],right[2],origM[2], LI[2],LK[2],LJ[2],RI[2],RK[2],RJ[2],origLI[2],origLK[2],origLJ[2],origRI[2],origRK[2],origRJ[2];
	  FLOAT **origX = origmap->site;

	  for(int c = 0; c < colors; c++){
	    left[c] = nanomap->left[c];
	    right[c] = nanomap->right[c];
	    origM[c] = Gmap[origalign->mapid2]->numsite[c];

	    if(DEBUG) assert(align[1+c].numpairs >= 1);// HERE HERE : fix when alignments without labels in one color are possible

	    LI[c] = align[1+c].sites1[0];
	    LK[c] = align[1+c].sitesK1[0];
	    LJ[c] = align[1+c].sites2[0] + max(0, left[c] - 1);
	    RI[c] = align[1+c].sites1[align[1+c].numpairs-1];
	    RK[c] = align[1+c].sitesK1[align[1+c].numpairs-1];
	    RJ[c] = align[1+c].sites2[align[1+c].numpairs-1] + max(0, left[c] -1);
	  
	    if(DEBUG) assert(origalign[1+c].numpairs >= 1);// HERE HERE : fix when alignments without labels in one color are possible

	    origLI[c] = origalign[1+c].sites1[0];
	    origLK[c] = origalign[1+c].sitesK1[0];
	    origLJ[c] = origalign[1+c].sites2[0];
	    origRI[c] = origalign[1+c].sites1[origalign[1+c].numpairs-1];
	    origRK[c] = origalign[1+c].sitesK1[origalign[1+c].numpairs-1];
	    origRJ[c] = origalign[1+c].sites2[origalign[1+c].numpairs-1];
	  }
	  double Ygap = 0.0;
	  for(int c = 0; c < colors; c++)
	    Ygap += Yc(YY[c],LI[c],LK[c]) - Yc(YY[c],origRI[c],origRK[c]);
	  Ygap *= 0.5;

	  if(Ygap >= 0.0){
	    double Xgap = 0;
	    for(int c = 0; c < colors; c++){
	      //	    if(DEBUG) assert(LI[c] >= origRI[c]);
	      Xgap += origX[c][LJ[c]] - origX[c][origRJ[c]];
	    }
	    Xgap *= 0.5;

	    /*	    if(VERB>=2 && Gmap[mapid]->id==96){
	      printf("refid=%d,mapid=%d:origid=%d(M=%d,%d),OR=%d,left=%d,%d,right=%d,%d:Ygap=%0.3f,Xgap=%0.3f,Ylambda=%0.3f\n",
		     refid,mapid,origalign->mapid2,origM[0],origM[1],
		     align->orientation,left[0],left[1],right[0],right[1],Ygap,Xgap,Ylambda);
	      fflush(stdout);
	      }*/

	    if(Xgap > 0.0 && Xgap < Ygap + Ylambda[0]+Ylambda[1] && Ygap-Xgap <= 10.0*(Ylambda[0]+Ylambda[1])){
	      if(VERB>=2){
		printf("refid=%d,mapid=%d:origid=%d(M=%d,%d),OR=%d,left=%d,%d,right=%d,%d:Ygap=%0.3f,Xgap=%0.3f\n",
		       refid,mapid,origalign->mapid2,origM[0],origM[1],
		       align->orientation,left[0],left[1],right[0],right[1],Ygap,Xgap);
		fflush(stdout);
	      }
	      tandem++;

              #pragma omp atomic
	      chimpaircnt++;

	      if(DEBUG) assert(!align->chimpair);
	      align->chimpair = 1;
	      //	      if(Xgap > 0.0)
	      //		if(DEBUG) assert(LJ > origRJ);
	    }
	  } else {
	    Ygap = 0.0;
	    for(int c = 0; c < colors; c++)
	      Ygap += Yc(YY[c],origLI[c],origLK[c]) - Yc(YY[c],RI[c],RK[c]);
	    Ygap *= 0.5;

	    if(Ygap >= 0.0){
	      double Xgap = 0.0;
	      for(int c = 0; c < colors; c++){
		//	      if(DEBUG) assert(origLI[c] >= RI[c]);
		Xgap += origX[c][origLJ[c]] - origX[c][RJ[c]];
	      }
	      Xgap *= 0.5;

	      /*	      if(VERB>=2 && Gmap[mapid]->id==96){
		printf("refid=%d,mapid=%d:origid=%d(M=%d),OR=%d,left=%d,right=%d:Yngap=%0.3f(%d,%d..%d,%d),Xngap=%0.3f(%d..%d),Ylambda=%0.3f\n",
		       refid,mapid,origalign->mapid2,origM,
		       align->orientation,left,right,Ygap,RI,RK,origLI,origLK,Xgap,RJ,origLJ,Ylambda);
		fflush(stdout);
		}*/
	      if(Xgap > 0.0 && Xgap < Ygap + (Ylambda[0]+Ylambda[1]) && Ygap - Xgap <= 10.0*(Ylambda[0]+Ylambda[1])){
		if(VERB>=2){
		  printf("refid=%d,mapid=%d:origid=%d(M=%d,%d),OR=%d,left=%d,%d,right=%d,%d:Yngap=%0.3f,Xngap=%0.3f\n",
			 refid,mapid,origalign->mapid2,origM[0],origM[1],
			 align->orientation,left[0],left[1],right[0],right[1],Ygap,Xgap);
		  fflush(stdout);
		}
		tandem++;

                #pragma omp atomic
		chimpaircnt++;

		if(DEBUG) assert(!align->chimpair);
		align->chimpair = 1;
		//		if(Xgap > 0.0)
		//		  if(DEBUG) assert(origLJ > RJ);
	      }
	    }
	  }
	} else {/* align->orientation == 1 */
	  int left[2],right[2],origM[2], origN[2], LI[2],LK[2],LJ[2],RI[2],RK[2],RJ[2],origLI[2],origLK[2],origLJ[2],origRI[2],origRK[2],origRJ[2];
	  FLOAT **origX = origmap->site;

	  for(int c = 0; c < colors; c++){
	    left[c] = nanomap->left[c];
	    right[c] = nanomap->right[c];
	    origN[c] = rmap->numsite[c];
	    origM[c] = Gmap[origalign->mapid2]->numsite[c];

	    if(DEBUG) assert(align[1+c].numpairs >= 1);// HERE HERE : fix when alignments without labels in one color are possible

	    LI[c] = align[1+c].sites1[0];
	    LK[c] = align[1+c].sitesK1[0];
	    LJ[c] = (MM[c] + 1 - align[1+c].sites2[0]) + max(0,left[c]-1);
	    RI[c] = align[1+c].sites1[align[1+c].numpairs-1];
	    RK[c] = align[1+c].sitesK1[align[1+c].numpairs-1];
	    RJ[c] = (MM[c] + 1 - align[1+c].sites2[align[1+c].numpairs-1]) + max(0,left[c] -1);

	    if(DEBUG) assert(origalign[1+c].numpairs >= 1);// HERE HERE : fix when alignments without labels in one color are possible

	    origLI[c] = origalign[1+c].sites1[0];
	    origLK[c] = origalign[1+c].sitesK1[0];
	    origRI[c] = origalign[1+c].sites1[origalign[1+c].numpairs-1];
	    origRK[c] = origalign[1+c].sitesK1[origalign[1+c].numpairs-1];
	    origLJ[c] = origM[c] + 1 - origalign[1+c].sites2[0];
	    origRJ[c] = origM[c] + 1 - origalign[1+c].sites2[origalign[1+c].numpairs-1];

	    if(DEBUG) assert(origRI[c]-origRK[c] >= 1 && origRI[c] <= origN[c]);
	    if(DEBUG) assert(origRJ[c] >= 1 && origRJ[c] <= origM[c]);
	  }
	  double Ygap = 0.0;
	  for(int c = 0; c < colors; c++)
	    Ygap += Yc(YY[c],LI[c],LK[c]) - Yc(YY[c],origRI[c],origRK[c]);
	  Ygap *= 0.5;

	  if(Ygap >= 0.0){
	    double Xgap = 0.0;
	    for(int c = 0; c < colors; c++){
	      //	    if(DEBUG) assert(LI[c] >= origRI[c]);
	      Xgap += origX[c][origRJ[c]] - origX[c][LJ[c]];
	    }
	    Xgap *= 0.5;

	    /*	    if(VERB>=2 && Gmap[mapid]->id==96){
	      printf("refid=%d,mapid=%d:origid=%d(M=%d),OR=%d,left=%d,right=%d:Ygap=%0.3f(origRI=%d,%d..LI=%d,%d),Xgap=%0.3f(LJ=%d..origRJ=%d),Ylambda=%0.3f\n",
		     refid,mapid,origalign->mapid2,origM,
		     align->orientation,left,right,Ygap,origRI,origRK,LI,LK,Xgap,LJ,origRJ,Ylambda);
	      fflush(stdout);
	      }*/

	    if(Xgap > 0.0 && Xgap < Ygap + (Ylambda[0]+Ylambda[1]) && Ygap-Xgap <= 10.0*(Ylambda[0]+Ylambda[1])){
	      if(VERB>=2){
		printf("refid=%d,mapid=%d:origid=%d(M=%d,%d),OR=%d,left=%d,%d,right=%d,%d:Ygap=%0.3f,Xgap=%0.3f\n",
		       refid,mapid,origalign->mapid2,origM[0],origM[1],
		       align->orientation,left[0],left[1],right[0],right[1],Ygap,Xgap);
		fflush(stdout);
	      }
	      tandem++;

              #pragma omp atomic
  	      chimpaircnt++;

	      if(DEBUG) assert(!align->chimpair);
	      align->chimpair = 1;
	      //	      if(Xgap > 0.0)
	      //		if(DEBUG) assert(origRJ > LJ);
	    }
	  } else {
	    Ygap = 0.0;
	    for(int c = 0; c < colors; c++)
	      Ygap += Yc(YY[c],origLI[c],origLK[c]) - Yc(YY[c],RI[c],RK[c]);
	    Ygap *= 0.5;

	    if(Ygap >= 0.0){
	      double Xgap = 0.0;
	      for(int c = 0; c < colors; c++){
		//	      if(DEBUG) assert(origLI[c] >= RI[c]);
		Xgap += origX[c][RJ[c]] - origX[c][origLJ[c]];
	      }
	      Xgap *= 0.5;

	      /*	      if(VERB && Gmap[mapid]->id==96){
		printf("refid=%d,mapid=%d:origid=%d(M=%d),OR=%d,left=%d,right=%d:Yngap=%0.3f(RI=%d,%d..origLI=%d,%d),Xngap=%0.3f(origLJ=%d..RJ=%d),Ylambda=%0.3f\n",
		       refid,mapid,origalign->mapid2,origM,
		       align->orientation,left,right,Ygap,RI,RK,origLI,origLK,Xgap,origLJ,RJ,Ylambda);
		fflush(stdout);
		}*/
	      if(Xgap > 0.0 && Xgap < Ygap + (Ylambda[0]+Ylambda[1]) && Ygap - Xgap <= 10.0*(Ylambda[0]+Ylambda[1])){
		if(VERB>=2){
		  printf("refid=%d,mapid=%d:origid=%d(M=%d,%d),OR=%d,left=%d,%d,right=%d,%d:Yngap=%0.3f,Xngap=%0.3f\n",
			 refid,mapid,origalign->mapid2,origM[0],origM[1],
			 align->orientation,left[0],left[1],right[0],right[1],Ygap,Xgap);
		  fflush(stdout);
		}
		tandem++;

                #pragma omp atomic
		chimpaircnt++;

		if(DEBUG) assert(!align->chimpair);
		align->chimpair = 1;
		//		if(Xgap > 0.0)
		//		  if(DEBUG) assert(RJ > origLJ);
	      }
	    }
	  }
	}
      } 
      if(VERB >= EVERB && !tandem){/* not a tandem alignment : a true chimerism */
	printf("chimeric alignment: refid=%d,mapid=%d(id=%lld):orientation=%d,origorientation=%d\n",
	       refid,origalign->mapid2,Gmap[origalign->mapid2]->id,align->orientation,origalign->orientation);
	fflush(stdout);
      }
    }

    /* check if chimeric split is required */
    if((align[1].Rend <= -2 || align[1].Lend <= -2) && (!nosplit || !nanomap->origmap) && nosplit <= 1){
      int *MM = nanomap->numsite;
      int J[2] = {align[1].sites2[0],align[2].sites2[0]};
      if(align[1].Lend <= -2 && J[0] + J[1] - 2 >= AlignedSiteThreshold){ /* save copy of left remainder of Molecule and append it to the end of maps[] */
	/* left end is cut off 0.001 kb left of sites J[0],J[1] (whichever is leftmost) */
	if(DEBUG) assert(align[2].Lend <= -2);

	Cmap *origmap = Gmap[mapid];
	Cmap *newmap;

	FLOAT **XX = nanomap->site;

	/* locate rightmost label J1[c]-1 in each color c = 0..1 that is left of aligned region */
	int J1[2];
	for(int c = 0;c < colors; c++){
	  int c2 = 1-c;
	  if(!align->orientation){
	    for(J1[c] = 1; J1[c] <= MM[c]; J1[c]++)
	      if(XX[c][J1[c]] >= XX[c2][J[c2]])
		break;
	  } else {
	    for(J1[c] = 1; J1[c] <= MM[c]; J1[c]++)
	      if(XX[c][MM[c]+1] - XX[c][MM[c]+1 - J1[c]] >= XX[c2][MM[c2]+1] - XX[c2][MM[c2]+1 - J[c2]])
		break;
	  }
	  J1[c] = min(J1[c], J[c]);
	}

        #pragma omp critical 	
	{
	  maxmapalloc(totalmaps+1,maxmaps,Gmap,0,1);
	  XXmap = Gmap;
	  newmap = Gmap[totalmaps];
	  if(!newmap)
	    Gmap[totalmaps] = newmap = new Cmap;
	  newmap->mapid = totalmaps++;
	}

	newmap->allfree();
	newmap->origmap = origmap;
	newmap->id = origmap->id;

	for(int c = 0; c < colors; c++){
	  newmap->numsite[c] = J1[c] - 1;
	  newmap->site[c] = new FLOAT[newmap->numsite[c] + 2];
	  newmap->site[c][0] = 0.0;

	  if(!align->orientation){
	    FLOAT len = min(XX[0][J1[0]],XX[1][J1[1]]) - 0.001;
	    newmap->left[c] = 0;
	    newmap->right[c] = J1[c] - 1;
	    for(int i = 1; i < J1[c]; i++)
	      newmap->site[c][i] = XX[c][i];
	    newmap->site[c][J1[c]] = len;
	    if(DEBUG>=2)
	      for(int i = 1; i <= newmap->numsite[c]+1; i++)
		assert(newmap->site[c][i] > newmap->site[c][i-1]);
	  } else {
	    FLOAT len = min(XX[0][MM[0]+1] - XX[0][MM[0]+1 - J1[0]], XX[1][MM[1]+1] - XX[1][MM[1]+1 - J1[1]]) - 0.001;
	    newmap->right[c] = MM[c] + 1;
	    newmap->left[c] = MM[c]+1 - (J1[c]-1);
	    for(int i = 1; i <= J1[c]; i++)
	      newmap->site[c][i] = XX[c][newmap->left[c] + i-1] - (XX[c][MM[c]+1] - len);
	    if(DEBUG>=2)
	      for(int i = 1; i <= newmap->numsite[c]+1; i++)
		assert(newmap->site[c][i] > newmap->site[c][i-1]);
	  }
	}

	newmap->len = nanomap->len * (newmap->site[0][newmap->numsite[0]+1]/XX[0][MM[0]+1]);
	newmap->fpcnt = -1;

	if(VERB>=2){
	  if(origmap->name)
	    printf("refid=%d,mapid=%d(id=%lld):Lend=%d,%d,Rend=%d,%d,score=%0.6f(%0.6f,%0.6f),or=%d:left fragment at maps[%d],sites=%d,%d(left=%d,%d,right=%d,%d,M=%d,%d,len=%0.3f kb),name=%s\n",
		   refid,mapid,origmap->id,align[1].Lend,align[2].Lend,align[1].Rend,align[2].Rend,align->score,align[1].score,align[2].score,align->orientation, 
		   totalmaps-1,newmap->numsite[0],newmap->numsite[1],newmap->left[0],newmap->left[1],newmap->right[0],newmap->right[1],MM[0],MM[1],newmap->len, origmap->name);
	  else
	    printf("refid=%d,mapid=%d(id=%lld):Lend=%d,%d,Rend=%d,%d,score=%0.6f(%0.6f,%0.6f),or=%d:left fragment at maps[%d],sites=%d,%d(left=%d,%d,right=%d,%d,M=%d,%d,len=%0.3f kb)\n",
		   refid,mapid,origmap->id,align[1].Lend,align[2].Lend,align[1].Rend,align[2].Rend,align->score,align[1].score,align[2].score,align->orientation, 
		   totalmaps-1,newmap->numsite[0],newmap->numsite[1],newmap->left[0],newmap->left[1],newmap->right[0],newmap->right[1],MM[0],MM[1],newmap->len);
	}

	if(phash && hashdelta && HashMultiMatch){ /* use same hashtable entries to try to alignment newmap */
	  size_t origalignid = alignid;

	  refalignSD(rmap, newmap, alignid, Kmax, refid, newmap->mapid, phash, XPen, Qmin, Qmax, A, Fmem, Imem, StrideMax, JMIN, JMAX, Imin, Imax, tid, limN, 2);

	  if(VERB>=2){
	    Calign *align = alignment[origalignid];
	    if(align && AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold)){
              #pragma omp critical
	      {
		printf("split alignment[%lu]:refid=%d,mapid=%d(id=%lld),or=%d:score=%0.4f,logPV=%0.3f,pairs=%d\n",
		       origalignid,align->mapid1,align->mapid2,Gmap[align->mapid2]->id,align->orientation,align->score,align->logPV,align->numpairs);
		fflush(stdout);
	      }
	    }
	  }
	}
      }

      if(align[1].Rend <= -2 && (MM[0]+MM[1] - J[0] - J[1]) >= AlignedSiteThreshold){/* save copy of right remainder of Molecule and append it to the end of maps[] */
	/* right end is cut off 0.001 kb right of sites J[0],J[1] (whichever is righmost) */

	if(DEBUG) assert(align[2].Rend <= -2);
	Cmap *origmap = nanomap;
	Cmap *newmap;

	/* locate leftmost label J1[c] + 1 in each color c = 0..1 that is right of the aligned region */
	int J1[2];
	FLOAT **XX = nanomap->site;
	for(int c = 0; c < colors; c++){
	  int c2 = 1-c;
	  if(!align->orientation){
	    for(J1[c] = MM[c]; J1[c] >= 1; J1[c]--)
	      if(XX[c][J1[c]] <= XX[c2][J[c2]])
		break;
	  } else {
	    for(J1[c] = MM[c]; J1[c] >= 1; J1[c]--)
	      if(XX[c][MM[c]+1] - XX[c][MM[c]+1 - J1[c]] <= XX[c2][MM[c2]+1] - XX[c2][MM[c2]+1-J[c2]])
		break;
	  }
	  J1[c] = max(J1[c],J[c]);
	}

        #pragma omp critical 
	{ 
	  maxmapalloc(totalmaps+1,maxmaps,Gmap,0,1);
	  XXmap = Gmap;
	  newmap = Gmap[totalmaps];
	  if(!newmap)
	    Gmap[nummaps] = newmap = new Cmap;
	  newmap->mapid = totalmaps++;
	}

	newmap->allfree();
	newmap->origmap = origmap;
	newmap->id = origmap->id;

	for(int c = 0; c < colors; c++){
	  newmap->numsite[c] = MM[c] - J1[c];
	  newmap->site[c] = new FLOAT[newmap->numsite[c] + 2];
	  newmap->site[c][0] = 0.0;

	  if(!align->orientation){
	    FLOAT len = min(XX[0][MM[0]+1] - XX[0][J1[0]], XX[1][MM[1]+1] - XX[1][J1[1]]) - 0.001;
	    newmap->left[c] = J1[c] + 1;
	    newmap->right[c] = MM[c] + 1;
	    for(int i = 1; i <= MM[c]+1-J1[c]; i++)
	      newmap->site[c][i] = len - (XX[c][MM[c]+1] - XX[c][J1[c] + i]);
	    if(DEBUG>=2)
	      for(int i = 1; i <= newmap->numsite[c]+1; i++)
		assert(newmap->site[c][i] > newmap->site[c][i-1]);
	  } else {
	    FLOAT len = min(XX[0][MM[0]+1-J1[0]],XX[1][MM[1]+1-J1[1]]) - 0.001;
	    newmap->left[c] = 0;
	    newmap->right[c] = MM[c] - J1[c];
	    for(int i = 1; i <= MM[c] - J1[c]; i++)
	      newmap->site[c][i] = origmap->site[c][i];
	    newmap->site[c][MM[c]-J1[c]+1] = len;
	    if(DEBUG>=2)
	      for(int i = 1; i <= newmap->numsite[c]+1; i++)
		assert(newmap->site[c][i] > newmap->site[c][i-1]);
	  }
	}

	newmap->len = nanomap->len * (newmap->site[0][newmap->numsite[0]+1]/XX[0][MM[0]+1]);
	newmap->fpcnt = -1;

	if(VERB>=2){
	  if(origmap->name)
	    printf("refid=%d,mapid=%d(id=%lld):Lend=%d,%d,Rend=%d,%d,score=%0.6f(%0.6f,%0.6f),or=%d:right fragment at maps[%d],sites=%d,%d(left=%d,%d,right=%d,%d,M=%d,%d,len=%0.3f kb),name=%s\n",
		   refid,mapid,origmap->id,align[1].Lend,align[2].Lend,align[1].Rend,align[2].Rend,align->score,align[1].score,align[2].score,align->orientation,
		   totalmaps-1,newmap->numsite[0],newmap->numsite[1],newmap->left[0],newmap->left[1],newmap->right[0], newmap->right[1], MM[0], MM[1], newmap->len, origmap->name);
	  else
	    printf("refid=%d,mapid=%d(id=%lld):Lend=%d,%d,Rend=%d,%d,score=%0.6f(%0.6f,%0.6f),or=%d:right fragment at maps[%d],sites=%d,%d(left=%d,%d,right=%d,%d,M=%d,%d,len=%0.3f kb)\n",
		   refid,mapid,origmap->id,align[1].Lend,align[2].Lend,align[1].Rend,align[2].Rend,align->score,align[1].score,align[2].score,align->orientation,
		   totalmaps-1,newmap->numsite[0],newmap->numsite[1],newmap->left[0],newmap->left[1],newmap->right[0], newmap->right[1], MM[0], MM[1], newmap->len);
	}
	if(phash && hashdelta && HashMultiMatch){ /* use same hashtable entries to try to alignment newmap */
	  size_t origalignid = alignid;

	  refalignSD(rmap, newmap, alignid, Kmax, refid, newmap->mapid, phash, XPen, Qmin, Qmax, A, Fmem, Imem, StrideMax, JMIN, JMAX, Imin, Imax, tid, limN, 2);

	  if(VERB>=2){
	    Calign *align = alignment[origalignid];
	    if(align && AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold)){
              #pragma omp critical
	      {
		printf("split alignment[%lu]:refid=%d,mapid=%d(id=%lld),or=%d:score=%0.4f,logPV=%0.3f,pairs=%d\n",
		       origalignid,align->mapid1,align->mapid2,Gmap[align->mapid2]->id,align->orientation,align->score,align->logPV,align->numpairs);
		fflush(stdout);
	      }
	    }
	  }
	}
      }
    }
  } /* good alignment */
}

static double logLR(int mapid, int refid, int c,
		    FLOAT *X, FLOAT *Y, 
		    int M,int N,
		    Calign *align,
		    int I, int J, int K,/* K==0 unless resSD > 0 */
		    int RJ,/* right end of alignment */
		    int Lij, int Rij, int Lijx, int Rijx,
		    int update,/**< 0 if called to verify original logLR, 1 if an error parameter update is being computed : 
				  If update==1 : freeze the original outlier calls, but also check LR with updated outlier calls to make sure it is equal or better than with original outlier calls
				*/
		    int verb,
		    double *LLsd,/* If != 0 : set to logLR component that is a function of SD and SF */
		    double *Xlen)/* If != 0 : set to sum of X interval lengths that contributed to Frate terms */
{
  if(DEBUG) assert(I == align[c].sites1[0]);
  if(DEBUG) assert(K == align[c].sitesK1[0]);
  if(DEBUG) assert(J == align[c].sites2[0]);

  LFLOAT F = FP[c]*0.01;
  LFLOAT scale = ScaleDeltaBPP ? 1.0 : (align->scaleID > 1) ? ScaleFactor[align->scaleID] : 1.0;
  LFLOAT XlenB = 0.0;
  int Jmin = ((NEW >= 2) ? (align[c].Lend != -2) : (align[c].Lend > -2)) ? 0 : J;
  int Jmax = ((NEW >= 2) ? (align[c].Rend != -2) : (align[c].Rend > -2)) ? M+1 : RJ;
  for(int j= Jmin;j < Jmax;j++)
    XlenB += max((FLOAT)0.0,X[j+1]-X[j] - resKB2[c]);
  XlenB *= scale;
  LFLOAT TotalBias = (Jmax-Jmin-2)*FBias[c] + EBias[c]*((Jmin==0?1:0)+(Jmax==M+1?1:0));
  TotalBias += XlenB*(Frate[c] + LBias[c]);
  TotalBias *= biasWT;

  LFLOAT logLR1 = align[c].score - TotalBias;
  if(0 && !REFDEBUG && !LLsd)
    return logLR1;

  if(DEBUG) assert(NEW >= 2);
  if(DEBUG) assert(resKB2[c] == 0.0);

  LFLOAT Bscale = 1.0-BETA;

  /* The following variables are used for full LR evaluation */
  LFLOAT logLR2 = 0.0;/* 2nd way to compute -log(LR) */
  LFLOAT score2 = 0.0;/* 2nd way to compute S(X) = -log(LR(X)) -E(-log(LR(X))|X) */
  LFLOAT TotBias2 = 0.0;/* 2nd way to compute TotalBias = -E(-log(LR(X))|X)*/
  
  /* The following variables are used with update==1 to compute LR with outlier calls frozen as specified in align[] */
  LFLOAT logLR3 = 0.0;
  LFLOAT score3 = 0.0;
  LFLOAT score3sd = 0.0;/* score3sd is the part of score3 that is a function of SD or SF */
  LFLOAT Xlen3 = 0.0;/* Xlen3 is the sum of X interval lengths that contribute to Frate based terms in score3 */

  /* left end fragment */
  LFLOAT x = scale * ((align[c].orientation==0) ? X[J] : X[M+1] - X[M+1-J]);
  LFLOAT xLijx_1 = (Lijx <= 0) ? 0.0 : scale * ((align[c].orientation==0) ? X[Lijx-1] : X[M+1] - X[M+1-(Lijx-1)]);
  int n = I - K - Lij + 1;
  int m = J;
  if(DEBUG>=2) assert(n >= 1);
  if(DEBUG>=2) assert(m >= 1);

  LFLOAT term;
  LFLOAT xB = max(0.0,x - resKB2[c]);
  if(VERB && verb){
    printf("logLR:refid=%d(id=%lld),mapid=%d(id=%lld),orientation=%d,c=%d:res[c]=%0.6f,resSD[c]=%0.6f,Xtheta=%0.8f,F=%0.6f,FN=%0.6f,scale=%0.8f,X[]=%p\n",
	   refid,refmap[refid]->id,mapid,Gmap[mapid]->id,align[c].orientation,c,res[c],resSD[c],Xtheta[c],F,FN[c],scale,X);
    fflush(stdout);
  }

  LFLOAT sm = Sm(J,I,K,Y,c);
  double sdterm = 0.0;
  LFLOAT Bias = biasWT * (xB*(Frate[c] + LBias[c]) + (m-1)*FBias[c] + EBias[c]);

  if(align->Lend <= -2){/* end outlier */
    Bias = BiasEnd2(x,J,c);
    term = ChimScore;
  } else if(ENDFIX && x <= Yc(Y,I,K) && Lij > 0 && Yc(Y,I,K)-Y[Lij-1] < x){ /* Sbnd(x,Y[I,K]-Y[Lij-1],n,m) ( NOTE: RES_VARIANCE is ignored) */
    LFLOAT y = Yc(Y,I,K)-Y[Lij-1];
    LFLOAT var = SF[c]*SF[c] + fabs(SD[c])*SD[c]*y;
    if(QUADRATIC_VARIANCE)
      var += SR[c]*SR[c]*y*y;
    LFLOAT err = x-y;
    term = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[c]) */ - err*err*0.5/var;
    if(SDDEBUG){
      score3sd += sdterm = -0.5*err*err/var;
      Xlen3 += xB;
    }
  } else if(ENDFIX && extend && x >= Yc(Y,I,K) && Lijx > 0 && x - xLijx_1 < Yc(Y,I,K)){/* Sbnd(x-xLijx_1, Yc(Y,I,K), J+1-Lijx, I-K) */
    LFLOAT y = Yc(Y,I,K);
    LFLOAT var = SF[c]*SF[c] + fabs(SD[c])*SD[c]*y;
    if(QUADRATIC_VARIANCE)
      var += SR[c]*SR[c]*y*y;
    x = x - xLijx_1;
    xB = max(LZERO,x - resKB2[c]);
    LFLOAT err = x - y;
    m = J+1-Lijx;
    n = I-K;
    if(DEBUG>=2) assert(m >= 1);
    if(DEBUG>=2) assert(n >= 1);
    Bias = biasWT * (xB*(Frate[c] + LBias[c]) + (m-1)*FBias[c] + EBias[c]);
    term = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[c]) */ - err*err*0.5/var;
    if(SDDEBUG){
      score3sd += sdterm = -0.5*err*err/var;
      Xlen3 += xB;
    }
  } else { /* Send(min(x,Yc(Y,I,K),J-1-max(1,Lijx),I-K+1 - max(1,Lij)) */
    if(extend){
      x = min(x,Yc(Y,I,K));
      xB = max(LZERO,x - resKB2[c]);
      m = J+1-max(1,Lijx);
      n = I-K+1 - max(1,Lij);
      if(DEBUG>=2) assert(m >= 1);
      if(DEBUG>=2) assert(n >= 1);
      Bias = biasWT * (xB*(Frate[c]+LBias[c]) + (m-1)*FBias[c] + EBias[c]);
    }
    term = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[c]) */;
    Xlen3 += xB;
  }

  logLR2 += term + sm;
  score2 += term + sm + Bias;
  logLR3 += term + sm;
  score3 += term + sm + Bias;

  if((DEBUG && !isfinite(term)) || (VERB && verb)){
    if(SCORE_APPROX>=2)
      printf("Lend:I=%d,K=%d,J=%d,Lij=%d,n=%d,m=%d,xB=%0.4f,Pen term=%0.6f(PenErr=%0.8f),Bias=%0.6f,Sm=%0.6f:\n\tlogLR2=%0.6f,score2=%0.6f\n",
	     I,K,J,Lij,n,m,xB,term,sdterm,Bias,sm,logLR2,score2);
    else
      printf("Lend:I=%d,K=%d,J=%d,Lij=%d,n=%d,m=%d,xB=%0.6f,Pen term=%0.6f(PenErr=%0.8f),Bias=%0.6f,Sm=%0.6f:\n\tlogLR2=%0.6f,score2=%0.6f,LLsd=%0.6f,Xlen=%0.6f\n",
	     I,K,J,Lij,n,m,xB,term,sdterm,Bias,sm,logLR2,score2,score3sd,Xlen3);
    for(int U = 0; U <= K; U++)
      printf("Y[I-%d]=%0.3f,",U,Y[I-U]);
    printf("\n");
    fflush(stdout);
    if(DEBUG) assert(isfinite(term));
  }
  TotBias2 += Bias;

  int H = I;
  int D = K;
  int G = J;
  for(int T=1;T < align[c].numpairs;T++, H = I, D = K, G = J){
    I = align[c].sites1[T];
    K = align[c].sitesK1[T];
    J = align[c].sites2[T];
    LFLOAT y = Yc(Y,I,K) - Yc(Y,H,D);
    x = scale * ((align[c].orientation==0) ? X[J] - X[G] : X[M+1-G] - X[M+1-J]);

    n = I - K - H;
    m = J - G;
    LFLOAT var = SF[c]*SF[c] + fabs(SD[c])*SD[c]*y;
    if(QUADRATIC_VARIANCE)
      var += SR[c]*SR[c]*y*y;
    if(RES_VARIANCE){
      FLOAT resR = Y[I] - Y[I-K];
      FLOAT resL = Y[I-K-n] - Y[I-K-n-D];
      FLOAT resvar = resL*resL + resR*resR;
      var += SE[c]*SE[c] * resvar;
    }
    LFLOAT err = x-y;
    LFLOAT Gtermsd = 0.5*log(2.0*M_PI*var);
    LFLOAT Gterm = log(Xtheta[c]) - Gtermsd;
    logLR1 += Bscale * Gterm;// Adjustment for BETA != 1.0

    xB = max(LZERO,x - resKB2[c]);
    term = Gterm - err*err*0.5/var - xB*Frate[c] - (n-1)*FnPenalty[c] - 
      (m-1)*FpPenalty[c] /* + log(1.0-FN[c])*/ + log(Pr(Y[I-K]-Y[H],c));
    LFLOAT sm = Sm(J,I,K,Y,c);
    if(SDDEBUG){
      score3sd += sdterm = Gterm - err*err*0.5/var;
      Xlen3 += xB;
    }
    Bias = biasWT * (xB*(Frate[c] + LBias[c]) + m * FBias[c]);
    LFLOAT rterm = term + sm;
    LFLOAT rterm3 = rterm;

    if(OUTLIER(maptype,x,y)){/* adjust rterm for outlier */
      LFLOAT OutPen = OutlierPenalty[c];
      if(outlierBC)
	OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
      if(FIX_OUTLIER_MIS)
	rterm = OutlierBias + max(term,OutPen) + sm;
      else
	rterm = OutlierBias + max(term + sm,OutPen);
    }
    logLR2 += rterm;
    score2 += rterm + Bias;

    /* adjust rterm3 & score3sd for previously classified outlier */
    LFLOAT OutPen = 0.0;
    if(align[c].outscore[T] + OUTLIER_MARGIN < align[c].iscore[T]){/* previously classified as outlier : evaluate as such */
      OutPen = OutlierPenalty[c];      
      if(outlierBC)
	OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
      if(FIX_OUTLIER_MIS)
	rterm3 = OutlierBias + OutPen + sm;
      else
	rterm3 = OutlierBias + OutPen;
      if(SDDEBUG){/* don't count sd,Xlen term for previously classified outlier, so value is comparable with value computed from errors[] */
	score3sd -= sdterm;
	Xlen3 -= xB;
      }
    }
    logLR3 += rterm3;
    score3 += rterm3 + Bias;

    if((DEBUG && !isfinite(term)) || verb){
      if(SCORE_APPROX>=2){
	if(update)
	  printf("T=%d:I=%d,K=%d,J=%d,H=%d,D=%d,G=%d:n=%d,m=%d,xB=%0.4f,Pen=%0.6f(Gterm=%0.6f,ErrPen=%0.6f,x=%0.4f,y=%0.4f,X*Frate=%0.6f,FnPen=%0.6f,FpPen=%0.6f,PenPr=%0.6f),Bias=%0.6f,Sm=%0.6f:\n\tlogLR2=%0.6f(delta=%0.6f),score2=%0.6f(delta=%0.6f),logLR3=%0.6f(delta=%0.6f)\n",
		 T, I, K, J, H, D, G, n, m, xB, term, Gterm, -err*err*0.5/var, x, y, -xB*Frate[c], -(n-1)*FnPenalty[c], -(m-1)*FpPenalty[c], log(Pr(Y[I-K]-Y[H],c)), Bias, sm, logLR2, rterm, score2, rterm + Bias, logLR3, rterm3);
	else
	  printf("T=%d:I=%d,K=%d,J=%d,H=%d,D=%d,G=%d:n=%d,m=%d,xB=%0.4f,Pen=%0.6f(Gterm=%0.6f(var=%0.6f,Xtheta=%0.6f),PenSD=%0.6f,ErrPen=%0.6f,x=%0.4f,y=%0.4f,X*Frate=%0.6f,FnPen=%0.6f,FpPen=%0.6f,PenPr=%0.6f),Bias=%0.6f,Sm=%0.6f:\n\tlogLR2=%0.6f(delta=%0.6f),score2=%0.6f(delta=%0.6f)\n",
		 T, I, K, J, H, D, G, n, m, xB, term, Gterm,var,Xtheta[c], sdterm, -err*err*0.5/var, x, y, -xB*Frate[c], -(n-1)*FnPenalty[c], -(m-1)*FpPenalty[c], log(Pr(Y[I-K]-Y[H],c)), Bias, sm, logLR2, rterm, score2, rterm + Bias);

      } else {
	if(update)
	  printf("T=%d:I=%d,K=%d,J=%d,H=%d,D=%d,G=%d:n=%d,m=%d,xB=%0.3f,Pen=%0.6f(Gterm=%0.6f,PenErr=%0.8f,var=%0.8f,x=%0.6f,y=%0.6f,Pr=%0.6f),Bias=%0.6f,Sm=%0.6f:\n\tlogLR2=%0.6f(delta=%0.6f),score2=%0.6f(delta=%0.6f),logLR3=%0.6f(delta=%0.6f,OutPen=%0.6f,delY=%0.6f,outscore=%0.6f,iscore=%0.6f),LLsd=%0.6f,X[%d]=%0.8f,X[%d]=%0.8f,Xlen=%0.6f\n",
		 T, I, K, J, H, D, G, n, m, xB, term, Gterm, -err*err*0.5/var, var, x, y, log(Pr(Y[I-K]-Y[H],c)), Bias, sm, logLR2, rterm, score2, rterm + Bias, logLR3, rterm3, OutPen,
		 Y[I]-Y[I-K],align[c].outscore[T],align[c].iscore[T],score3sd, align[c].orientation ? M+1-G : G, align[c].orientation ? X[M+1-G] : X[G] , align[c].orientation ? M+1-J : J, align[c].orientation ? X[M+1-J] : X[J], Xlen3);
	else
	  printf("T=%d:I=%d,K=%d,J=%d,H=%d,D=%d,G=%d:n=%d,m=%d,xB=%0.3f,Pen=%0.6f(Gterm=%0.6f,PenErr=%0.8f,var=%0.7f,x=%0.6f,y=%0.6f,Pr=%0.6f),Bias=%0.6f,Sm=%0.6f:\n\tlogLR2=%0.6f(delta=%0.6f),score2=%0.6f(delta=%0.6f)\n",
		 T, I, K, J, H, D, G, n, m, xB, term, Gterm, -err*err*0.5/var, var, x, y, log(Pr(Y[I-K]-Y[H],c)), Bias, sm, logLR2, rterm, score2, rterm + Bias);
      }
      //      printf("    FP[c]=%0.6f,F=%0.8f,LBias=%0.8f,Frate[c]=%0.8f,1/Xtheta[c]=%0.8f,FBias[c]=%0.8f\n",FP[c],F,LBias[c],Frate[c],1.0/Xtheta[c],FBias[c]);
      fflush(stdout);
      if(DEBUG) assert(isfinite(term));
    }
    TotBias2 += Bias;
  }
  /* right end fragment */
  if(DEBUG)assert(H==I && D==K && G==J);
  x = scale * ((align[c].orientation==0) ? X[M+1]-X[G] : X[M + 1 - G]);
  LFLOAT xRijx_1 = (Rijx > M) ? 0.0 : scale * ((align[c].orientation==0) ? X[M+1] - X[Rijx+1] : X[M+1-(Rijx+1)]);
  n = Rij - H + 1;
  m = M + 1 - G;

  xB = max(LZERO,x - resKB2[c]);
  sdterm = 0.0;
  Bias = biasWT * (xB*(Frate[c] + LBias[c]) + (m-1)*FBias[c] + EBias[c]);

  if(align->Rend <= -2){/* end outlier */
    Bias = BiasEnd2(x,M+1-J,c);
    term = ChimScore;
  } else if(ENDFIX && x <= Y[N+1]-Yc(Y,H,D) && Rij <= N && Y[Rij+1]-Yc(Y,H,D) < x){  /* Sbnd(x,Y[Rij+1]-Yc(Y,H,D),M+1-G,Rij+1-I) : RES_VARIANCE is ignored */
    double y = Y[Rij+1]-Yc(Y,H,D);
    double var = SF[c]*SF[c] + fabs(SD[c])*SD[c]*y;
    if(QUADRATIC_VARIANCE)
      var += SR[c]*SR[c]*y*y;
    double err = x-y;
    if(VERB && verb){
      printf("  Sbnd():H=%d,D=%d,G=%d,Rij=%d,N=%d,n=%d,m=%d,Y[N+1]=%0.6f,Y[Rij+1]=%0.6f,Yc(Y,H,D)=%0.6f,y=%0.6f,x=%0.6f,var=%0.6f\n",
	     H,D,G,Rij,N,n,m,Y[N+1],Y[Rij+1],Yc(Y,H,D),y,x,var);
      fflush(stdout);
    }
    term = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[c]) */ - err*err*0.5/var;
    if(SDDEBUG){
      score3sd += sdterm = -0.5 * err*err/var;
      Xlen3 += xB;
    }
  } else if(ENDFIX && extend && x >= Y[N+1] - Yc(Y,I,K) && Rijx <= M && x - xRijx_1 < Y[N+1]-Yc(Y,I,K)){/* Sbnd(x-xRijx_1, Y[N+1]-Yc(Y,I,K), Rijx+1-J,min(N,Rij)+1-I) */
    LFLOAT y = Y[N+1] - Yc(Y,I,K);
    LFLOAT var = SF[c]*SF[c] + fabs(SD[c])*SD[c]*y;
    if(QUADRATIC_VARIANCE)
      var += SR[c]*SR[c]*y*y;
    x = x - xRijx_1;
    xB = max(LZERO,x - resKB2[c]);
    LFLOAT err = x - y;
    m = Rijx + 1 - J;
    n = min(N,Rij) + 1 - I;
    if(DEBUG>=2) assert(m >= 1);
    if(DEBUG>=2) assert(n >= 1);
    if(VERB && verb){
      printf("  Sbnd(x,y,m,n,c)=%0.6f:H=%d,D=%d,G=%d,Rij=%d,Rijx=%d,n=%d,m=%d,X[M+1]-X[Rijx+1]=%0.6f,Y[N+1]=%0.6f,Yc(Y,H,D)=%0.6f,Y[N+1]-Yc(Y,I,K)=%0.6f,X[Rijx+1]-X[G]=%0.6f,var=%0.6f\n",
	     Sbnd(x,y,m,n,c), H,D,G,Rij,Rijx,n,m,xRijx_1,Y[N+1],Yc(Y,H,D),y,x,var);
      fflush(stdout);
    }
    Bias = biasWT * (xB*(Frate[c]+LBias[c]) + (m-1)*FBias[c] + EBias[c]);
    term = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[c]) */ - err*err*0.5/var;
    if(SDDEBUG){
      score3sd += sdterm = -0.5 * err*err/var;
      Xlen3 += xB;
    }
  } else {
    if(extend){
      x = min(x, Y[N+1]-Yc(Y,I,K));
      xB = max(LZERO,x - resKB2[c]);
      m = min(M,Rijx)+1-J;
      n = min(N,Rij)+1-I;
      if(DEBUG>=2) assert(m >= 1);
      if(DEBUG>=2) assert(n >= 1);
      Bias = biasWT * (xB*(Frate[c]+LBias[c]) + (m-1)*FBias[c] + EBias[c]);
    }
    term = -Frate[c]*xB - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] /* + 0.5*log(1.0-FN[c]) */;
    if(SDDEBUG)
      Xlen3 += xB;
  }
  logLR2 += term;
  score2 += term + Bias;
  logLR3 += term;
  score3 += term + Bias;
  if((DEBUG && !isfinite(term)) || (VERB && verb)){
    if(SCORE_APPROX>=2)
      printf("Rend=%d:H=%d,D=%d,G=%d,Rij=%d:n=%d,m=%d,xB=%0.3f,Pen=%0.6f(PenErr=%0.8f),Bias=%0.6f:\n\tlogLR1=%0.6f,logLR2=%0.6f,score2=%0.6f\n",
	     align->Rend,H,D,G,Rij,n,m,xB,term,sdterm,Bias,logLR1,logLR2,score2);
    else
      printf("Rend=%d:H=%d,D=%d,G=%d,Rij=%d:n=%d,m=%d,xB=%0.6f,Pen=%0.6f(PenErr=%0.8f),Bias=%0.6f:\n\tlogLR1=%0.6f,logLR2=%0.6f,score2=%0.6f,LLsd=%0.6f,Xlen=%0.6f,logLR3=%0.6f(delta=%0.6f)\n",
	     align->Rend,H,D,G,Rij,n,m,xB,term,sdterm,Bias,logLR1,logLR2,score2,score3sd,Xlen3,logLR3,term);
    printf("    Y[Rij]-Y[H,D]=%0.3f,Y[Rij-1]-Y[H,D]=%0.3f,Y[Rij+1]-Y[H,D]=%0.3f,x=%0.8f\n",
	   Y[Rij]-Yc(Y,H,D),Y[Rij-1]-Yc(Y,H,D),Y[Rij+1]-Yc(Y,H,D),x);
    fflush(stdout);
    if(DEBUG) assert(isfinite(term));
  }
  TotBias2 += Bias;

  if((REFDEBUG || DEBUG/* HERE HERE >=2 */) && !align->scaleID && !update && 
     fabs(align[c].score - score2) >= SCORE_MARGIN * max(1.0,fabs(score2)) * sqrt((double)align[c].numpairs)){
    printf("refid=%d(id=%lld),mapid=%d(id=%lld),orientation=%d,c=%d:logLR1=%0.6f,logLR2=%0.6f,logLR3=%0.6f,score=%0.6f,logPV=%0.2f,numpairs=%d,score2=%0.6f,TotBias=%0.6f,TotBias2=%0.6f\nresKB2=%0.4f,Xlen=%0.3f,XlenB=%0.3f\n",
	   refid,refmap[refid]->id,mapid,Gmap[mapid]->id,align[c].orientation,c,logLR1,logLR2,logLR3,align[c].score,align[c].logPV,align[c].numpairs,score2,TotalBias,TotBias2,resKB2[c],scale*X[M+1],XlenB);
    if(Poutlier > 0.0)
      printf("Poutlier=%0.6f,OutlierPenalty[c]=%0.4e,OutlierBias=%0.4f\n",
	     Poutlier,OutlierPenalty[c],OutlierBias);
    if(VERB && !verb)
      (void)logLR(mapid,refid,c,X,Y,M,N,align,align[c].sites1[0],align[c].sites2[0],align[c].sitesK1[0], RJ,Lij,Rij,Lijx,Rijx,update,1,LLsd,0);
    fflush(stdout);
    assert(fabs(align[c].score - score2) < SCORE_MARGIN * max(1.0,fabs(score2)) * sqrt((double)align[c].numpairs));
    exit(1);
  }
  if(REFDEBUG && update && !(logLR2 > logLR3 - SCORE_MARGIN * max(1.0,fabs(logLR2)) * sqrt((double)align->numpairs))){
    printf("refid=%d(id=%lld),mapid=%d(id=%lld),orientation=%d,c=%d:logLR1=%0.6f,logLR2=%0.6f,logLR3=%0.6f,score=%0.6f,numpairs=%d,score2=%0.6f,TotBias=%0.6f,TotBias2=%0.6f\nresKB2=%0.4f,Xlen=%0.3f,XlenB=%0.3f\n",
	   refid,refmap[refid]->id,mapid,Gmap[mapid]->id,align[c].orientation,c,logLR1,logLR2,logLR3,align[c].score,align[c].numpairs,score2,TotalBias,TotBias2,resKB2[c],scale*X[M+1],XlenB);
    if(VERB && !verb)
      (void)logLR(mapid,refid,c,X,Y,M,N,align,align[c].sites1[0],align[c].sites2[0], align[c].sitesK1[0], RJ,Lij,Rij,Lijx,Rijx,update,1,LLsd,0);
    fflush(stdout);
    assert(logLR2 > logLR3 - SCORE_MARGIN * max(1.0,fabs(logLR2)) * sqrt((double)align[c].numpairs));
    exit(1);
  }

  if(SDDEBUG && LLsd)
    *LLsd = score3sd;
  if(SDDEBUG && Xlen)
    *Xlen = Xlen3;

  return logLR3; // (REFDEBUG || update) ? logLR3 : logLR2;
}

/* compute log Likelihood for alternate value of FN */
static double FNtotLL(double newFN, int c, int numthreads, size_t numaligns, Calign **alignment, double *Tarray1)
{
  /* update score_init() terms that depend on FN : pFN,pTP,pTPd,log_pTP,log_qTPd, Xlambda, Ylambda, FBias, LBias, EBias, FnPenalty : assumes SCORE_Y==0 */
  double pFN = FN[c] = newFN;
  pTP[c] = pTPd[c] = 1.0 - pFN;
  log_pTP[c] = log_pTPd[c] = log(pTPd[c]);
  log_qTPd[c] = log(pFN);
  double F = FP[c]*0.01;
  double var = fabs(SD[c])*SD[c];
  double varF = SF[c]*SF[c];
  double varR = SR[c]*SR[c];
  Xlambda[c] = pTP[c]/(1.0/(Xtheta[c] - var*0.5) - F);/* HERE : does not include correction for SR[],SE[] */
  if(XmapCount[0] && RefRepeats==1){
    Ylambda[c] = Xlambda[c];
    Ygm[c] = XmapGtheta[0] * Xlambda[c]/Xtheta[c];
  }

  double YgmBpTP = (Ygm[c]/pTP[c]) * (Xlambda[c]/Ylambda[c]);
  double Sbias = log(Xtheta[c]) - 0.5*log(2.0*M_PI*(varF+var*YgmBpTP+varR*YgmBpTP*YgmBpTP));
  FBias[c] = 0.5*KBIAS - KFBIAS*pFN*log(pFN)/pTP[c] - log_pTP[c] - Sbias; /* Bias per site interval in X (overlapping Y)  = Abias - Sbias in document */
  LBias[c] = -KLBIAS * F*log(F*Xtheta[c]) - F* FBias[c];/* Bias per kb in X (overlapping Y) : In document LBias is larger by Frate = F - 1/Xtheta */
  EBias[c] = 0.5*(Sbias + FBias[c]);/* Bias per aligned end , corresponding to Abias/2 in the document. This is larger by KBIAS/4 compared to the document, so the final
				       score will be KBIAS/2 larger than in the document */
  if(resSD[c] > 0.0) /* PRbias =  Bias per site interval in X for log(Pr(y)) and Sm()-log(pTP) : will be estimated based on average of actual alignments in refalign */
    FBias[c] += PRbias[c];

  FnPenalty[c] = -log(pFN);/* missing cut penalty*/
  if(GAMMA != 1.0)
    FnPenalty[c] *= GAMMA;

  for(int tid = 0; tid < numthreads; tid++)
    Tarray1[tid] = 0.0;

  #pragma omp parallel num_threads(numthreads) if (numthreads > 1)
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num ();
#endif
    double myLRsum = 0.0;

    #pragma omp for schedule(dynamic,1)
    for(size_t i = 0; i < numaligns; i++){
      Calign *align = alignment[i];
      if(!align || align->numpairs <= 1)
	continue;
      int rid = align->mapid1;
      int mid = align->mapid2;
      if(BestRef){
	Cmap *origmap = Gmap[mid];
	while(origmap->origmap)
	  origmap = origmap->origmap;
	if(origmap->align->mapid1 != align->mapid1)
	  continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
      }
      FLOAT *X = Gmap[mid]->site[c];
      int M = Gmap[mid]->numsite[c];
      FLOAT *Y = refmap[rid]->site[c];
      int N = refmap[rid]->numsite[c];
      double LRsd = 0.0;
      double LR = logLR(mid,rid,c,X,Y,M,N,align,align->sites1[0], align->sites2[0], align->sitesK1[0], align->sites2[align->numpairs-1], align->Lij1,align->Rij1,align->Lij2,align->Rij2,1,0,&LRsd,NULL);
      if(MapRate > 0.0 && !(align->score > ScoreThreshold))
	continue;
      myLRsum += LR;
    }
    Tarray1[tid] = myLRsum;
  }

  qsort(Tarray1, numthreads, sizeof(double), (intcmp *)DoubleInc);
  double LRsum = 0.0;
  for(int tid = 0; tid < numthreads; tid++)
    LRsum += Tarray1[tid];

  return LRsum;
}

static double scoresum;
static double logPVsum;
static int NumpairsSum;
static double logLRsum, logLRsdsum, totXlen;

static double ATscoresum;
static double ATlogPVsum;
static int ATNumpairsSum;
static double ATlogLRsum;

static double *logLRarray = NULL, *logLRarrayC[MAXCOLOR], *logSDarray[MAXCOLOR], *Xlen[MAXCOLOR];

/* The alignments for each refmap i = 0..numrefmaps-1 are located in alignment[numalign_start[i]..numalign_end[i]-1] : useful when there are multiple reference maps */
extern int *orignummaps;
extern int warning2;

static int thread_printed = 0;

extern double mtime(),wtime();/* see refine.cpp */

#define MAX_PV 50
static double RawMappingRatePV[MAX_PV+1];
static double MappingRatePV[MAX_PV+1];

#if SIMPLE_BIAS <= 1
/* compute

   biassum = AVG((rawx - y)/var(y))
   Ivarsum = AVG(1/var(y))

   AVG(z) means average of z over errors[BinStart .. BinEnd-1] */

static void biascompute(double &biassum, double &Ivarsum, Cinterval *errors, int Bin, int BinStart, int BinEnd, double *Tarray1, double *Tarray2, int numthreads,
			double A, double B, double R, double E)
{
  for(int tid = 0; tid < numthreads; tid++)
    Tarray1[tid] = Tarray2[tid] = 0.0;

  // non-deterministic with gcc 4.6 unless -fno-unsafe-math-optimizations is used
  #pragma omp parallel num_threads(numthreads)
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num ();
#endif

    double mybiassum = 0.0, myIvarsum = 0.0;
    
    #pragma omp for schedule(static,64)
    for(int i = BinStart; i < BinEnd; i++){
      Cinterval *perr = &errors[i];
      double y = perr->y;
      if(ENDFIX>=3 && perr->end && y >= perr->x){
	if(VERB>=2)
	  printf("Bin=%d,i=%d:end=%d,x=%0.4f,rawx=%0.4f,y=%0.4f:skipping due to ENDFIX>=3\n",
		 Bin,i,perr->end,perr->x,perr->rawx,perr->y);
	continue;
      }
      double var = A + B * y;
      if(QUADRATIC_VARIANCE)
	var += R*y*y;
      if(RES_VARIANCE)
	var += E * perr->resvar;
      double Ivar = 1.0/var;
      double err = perr->rawx - y;
      mybiassum += err * Ivar;
      myIvarsum += Ivar;
      if(VERB>=2)
	printf("Bin=%d,i=%d:end=%d,x=%0.4f,rawx=%0.4f,y=%0.4f,err=%0.4f,var=%0.8f,Ivar=%0.8f:biassum=%0.8f,Ivarsum=%0.8f\n",
	       Bin,i,perr->end,perr->x,perr->rawx,perr->y,err,var,Ivar,biassum,Ivarsum);
    }

    Tarray1[tid] = mybiassum;
    Tarray2[tid] = myIvarsum;
  }
  qsort(Tarray1, numthreads, sizeof(double), (intcmp *)DoubleInc);
  qsort(Tarray2, numthreads, sizeof(double), (intcmp *)DoubleInc);

  if(VERB>=2){
    double tbiassum = 0.0, tIvarsum = 0.0;
    printf("biascompute:Bin=%d,BinStart=%d,BinEnd=%d,numthreads=%d\n",Bin,BinStart,BinEnd,numthreads);
    for(int tid = 0; tid < numthreads; tid++){
      tbiassum += Tarray1[tid];
      tIvarsum += Tarray2[tid];
      printf("Bin=%d:tid=%d:biassum=%0.18f(delta=%0.18f),Ivarsum=%0.18f(delta=%0.18f)\n",
	     Bin,tid,tbiassum, Tarray1[tid], tIvarsum, Tarray2[tid]);
    }
    fflush(stdout);
  }

  /* NOTE : above summation with printf is deterministic , the below identical summation is NOT ! */
  biassum = Ivarsum = 0.0;
  for(int tid = 0; tid < numthreads; tid++){
    biassum += Tarray1[tid];
    Ivarsum += Tarray2[tid];
  }
  if(DEBUG>=2) assert(isfinite(biassum));
  if(DEBUG>=2) assert(isfinite(Ivarsum));
}
#endif // SIMPLE_BIAS <= 1

/** Find best alignment between each nanomap and reference pair (in either orientation for the map). For each
   pair the best alignment is appended to the alignment list alignment[0..numaligns]
   If maps are chimeric, pieces of them may be appended to the end of maps[
*/

void refalign_allpairs2(int& nummaps,
			int& numrefmaps,
			double *pSNRtotLL)
{
  if(DEBUG) assert(colors == 2);

  if(AlignedOutlierLabels < 1000 || AlignedOutlierLabels2 < 1000){
    printf("refalign_allpairs2: -I 2nd arg and -MultiMatches <M2> not yet implemented\n");
    fflush(stdout);exit(1);
  }

  rs_heap = new lightweight_heap(0,0);

  if(pSNRtotLL != NULL){
    printf("-minSNRestimate : global scan not implemented for 2 colors\n");
    fflush(stdout);
  }
  if(BestRefExt){
    printf("-BestRefExt %0.3f not yet implemented for 2 colors\n",BestRefExt);
    exit(1);
  }
  if(DEBUG && !(DPsizI*sizeof(int) + DPsizF*sizeof(TFLOAT) == sizeof(DP))){
    printf("DPsizI= %lld, DPsizF= %lld, DPsizI*sizeof(int)+DPsizF*sizeof(TFLOAT)==%lld, sizeof(DP)=%lu\n",
	   DPsizI,DPsizF, DPsizI*sizeof(int)+DPsizF*sizeof(TFLOAT),sizeof(DP));
    fflush(stdout);
    assert(DPsizI*sizeof(int) + DPsizF*sizeof(TFLOAT) == sizeof(DP));
  }

  if(VERB>=2){
    printf("outlierExtend=%d,outlierBC=%d\n",outlierExtend,outlierBC);
    fflush(stdout);
  }
  if(RepeatMaxShift > 0){
    printf("-RepeatMask not yet implemtned for 2 colors\n");
    exit(1);
  }
  if(outlierExtend){// HERE HERE
    printf("-outlierExtend not yet implemented for 2 colors\n");
    exit(1);
  }

  if(ScanCorrection > 0 && UniqueScans > 1 /* HERE HERE && nummaps < totalmaps && !(hash_filename || NoSplit == 2) */){
    //    printf("-ScanCorrection with -subset is not support without EITHER -hash OR -nosplit 2\n");
    printf("-ScanCorrection is not yet supported with colors=%d\n",colors);
    exit(1);
  }
  if(!(ScanCorrection > 0 && UniqueScans > 1))
    startmaps = totalmaps = nummaps;

  if(Refine && NoSplit != 2){
    printf("-refine %d not supported without -nosplit 2\n",Refine);
    exit(1);
  }
  if(Refine)
    startmaps = totalmaps = nummaps;

  if(!FLAT_MEMORY){
    printf("FLAT_MEMORY==0 no longer supported\n");
    exit(1);
  }

  if(!ALIGN_COMPRESS){
    printf("ALIGN_COMPRESS==0 no longer supported\n");
    exit(1);
  }

  if(!SCORE_ALLREF){
    printf("SCORE_ALLREF=0 no longer supported\n");
    exit(1);
  }

  if(!GOODMAPS){
    printf("GOODMAPS=0 no longer supported\n");
    exit(1);
  }

  if(PoutlierS != Poutlier && MapStitched){
    printf("WARNING:Seperate -outlier pvalue for stitched intervals not yet implmented: using same pvalue\n");
    fflush(stdout);
  }

  if(DEBUG) assert(startmaps == totalmaps);

  origScoreThreshold = ScoreThreshold;/* global to support REFDEBUG */
  origLogPvThreshold = LogPvThreshold;

  if(REFDEBUG){
    if(!(origScoreThreshold <= -9999.0) || !(origLogPvThreshold <= -0.99) || !(AlignedSiteThreshold <= 2) || !(AlignedEndOutlierThreshold >= 2) || AlignedOutlierThreshold < 9999.0 || AlignedLengthThreshold > 0.0 || extend || PoutlierEnd > 0 || (Ch > 0.0 && ChFdr > 0.0) || NoSplit != 2){
      printf("REFDEBUG requires -S -1e4 -T 10 -A 2 -L 0 -E 2 -extend 0 -endoutlier 0.0 -nosplit 2\n");
      exit(1);
    }
    if(REFDEBUG && MININTERVAL > 1e-5){
      printf("REFDEBUG requires MININTERLVAL <= 0.00001 (in constants.h)\n");
      exit(1);
    }
  }

  if(VERB>=2){
    printf("refalign:nummaps=%d,totalmaps=%d,DIAGONAL=%d:start time=%0.6f(Elapsed=%0.6f),REFDEBUG=%d\n",nummaps,totalmaps,DIAGONAL,mtime(),wtime(),REFDEBUG);
    fflush(stdout);
  }

  //  int origCMapID = CMapID;

  if(!ISLINEAR && isLinear){
    printf("isLinear=%d not supported (set ISLINEAR in constants.h)\n",isLinear);
    exit(1);
  }

  if(NEW <= 0){
    printf("NEW=0 no longer supported\n");
    exit(1);
  }

  if(NEW && Ch > PoutlierEnd){
    if(VERB){
      printf("Using -PoutlierEnd %0.6e due to -Ch\n", Ch);
      fflush(stdout);
    }
    PoutlierEnd = Ch;
  }

  if(DEBUG>=2 && (usecolor || colors >= 2)){
    for(int i= 0; i < nummaps; i++){
      Cmap *pmap = Gmap[i];
      FLOAT len1 = pmap->site[0][pmap->numsite[0]+1];
      FLOAT len2 = pmap->site[1][pmap->numsite[1]+1];
      if(!(fabs(len1 - len2) < 1e-3)){
	printf("refalign: Gmap[%d]: id=%lld, M1=%d,len1=%0.4f, M2=%d,len2=%0.4f\n",i, pmap->id, pmap->numsite[0], len1, pmap->numsite[1], len2);
	fflush(stdout);
	assert(fabs(len1 - len2) < 1e-3);
      }
    }
  }
  if(DEBUG>=2){
    for(int mapid = 0; mapid < nummaps; mapid++){
      if(!(Gmap[mapid]->mapid == mapid)){
	printf("mapid=%d:Gmap[mapid]=%p,Gmap[mapid]->mapid=%d,Gmap[mapid]->id=%lld\n",
	       mapid,Gmap[mapid],Gmap[mapid]->mapid,Gmap[mapid]->id);
	for(int m = 0; m < 10; m++)
	  printf("m=%d:Gmap[m]=%p,Gmap[m]->mapid=%d,Gmap[m]->id=%lld\n",
		 m,Gmap[m],Gmap[m]->mapid,Gmap[m]->id);
	fflush(stdout);
	assert(Gmap[mapid]->mapid == mapid);
      }
    }
  }

  if(DEBUG && REFSCORETYPE != 4){
    printf("REFSCORETYPE=%d not supported(must be 4)\n",REFSCORETYPE);
    exit(1);
  }

  if(numrefmaps <= 0){
    printf("refalign_allpairs2: Insufficient data: nummaps=%d, numrefmaps=%d\n",nummaps,numrefmaps);
    exit(1);
  }

  if(NoStat && RefRepeats > 1){
    printf("Cannot use -nostat or -refine with -M %d\n", RefRepeats);
    exit(1);
  }

  double totlen = 0.0;
  for(int m = 0; m < nummaps; m++){
    Cmap *pmap = Gmap[m];
    int M = pmap->numsite[0];
    totlen += pmap->site[0][M+1];
  }
  double totreflen = 0.0;
  for(int m = 0; m < numrefmaps; m++){
    Cmap *pmap = refmap[m];
    int N = pmap->numsite[0];
    totreflen += pmap->site[0][N+1];
  }
  double CoverageMult = totlen/max(0.001,totreflen);
  if(totalmaps > nummaps){/* adjust for use of -subset */
    if(VERB){
      printf("Adjusting Coverage estimate for use of -subset : totalmaps=%d,nummaps=%d\n",totalmaps,nummaps);
      fflush(stdout);
    }
    CoverageMult *= totalmaps;
    CoverageMult /= nummaps;
  }

  int numthreads = 1;/* NOTE : inside the giter & refid loops local numthreads variables are used
			to customize numthreads to the current number iterations */
#ifdef _OPENMP
  numthreads = omp_get_max_threads();
  numthreads = min(numthreads,MaxThreads);
#endif

  double *Tarray1 = new double[numthreads];
  double *Tarray2 = new double[numthreads];
  double *Tarray3 = new double[numthreads];
  double *Tarray4 = new double[numthreads];
  double *Tarray5 = new double[numthreads];
  double *Tarray6 = new double[numthreads];

  extern Cmap **YYmap,**XXmap;
  extern int numY,numX;
  numY = numrefmaps;
  YYmap = refmap;
  numX = nummaps;
  XXmap = Gmap;

  /* -ScanScaling variables */
  int ScanAligned = ((UniqueScans*3*sizeof(double) + THREAD_PADDING - 1) & (~(THREAD_PADDING-1)))/sizeof(double);
  if(DEBUG) assert(ScanAligned >= UniqueScans * 3);

  //  double *Sxsum[MAXCOLOR],*Sxysum[MAXCOLOR],*Sxxsum[MAXCOLOR];
  double *ThetaDelta[MAXCOLOR],**Txy[MAXCOLOR],**Txx[MAXCOLOR],**Tx[MAXCOLOR];
  if(ScanCorrection){
    for(int c = 0; c < colors; c++){
      ThetaDelta[c] = new double[UniqueScans*4];
      //      Sxsum[c] = &ThetaDelta[c][UniqueScans];
      //      Sxysum[c] = &ThetaDelta[c][UniqueScans*2];
      //      Sxxsum[c] = &ThetaDelta[c][UniqueScans*3];
      Txy[c] = new double*[numthreads*3];
      Txx[c] = &Txy[c][numthreads];
      Tx[c] = &Txy[c][numthreads*2];

      for(int t = 0; t < numthreads; t++){
	double *mem = new double[ScanAligned];
	Txy[c][t] = mem;
	Txx[c][t] = &mem[UniqueScans];
	Tx[c][t] = &mem[UniqueScans*2];
      }
    }
  }

  if(ScanCorrection && UniqueScans > 1 && NoSplit != 2){
    printf("-ScanCorrection only supported with -nosplit 2\n");
    exit(1);
  }

  int Imaxresbias = 0, *SizeToBin[MAXCOLOR] = {0,0};
  int rawsitemaps = 0;
  if(VERB>=2){
    printf("NoStat=%d,maxresbias=%0.3f,mres=%0.3f,parametersfile=%p\n",NoStat,maxresbias,mres,parametersfile);
    fflush(stdout);
  }

  if(VERB>=2){
    printf("PixelLen=%0.4f,origPixelLen=%0.4f\n",PixelLen,origPixelLen);
    fflush(stdout);
  }

  if(!NoStat && maxresbias > mres * 0.5){/* initialize bias parameters */
    if(ScaleDelta && !ScaleDeltaBPP){
      printf("WARNING: -resbias with -ScaleDelta requires -ScaleDeltaBPP\n");
      if(REFDEBUG)
	exit(1);
    }

    if(0 && !giter2 && (!parametersfile || !ResBinsFile)){
      /* default initialization of ResBins[] and resbiasX[] resbias[] is not done in RefAligner */
    }
    if(rawsitemaps < totalmaps){
      /* allocate rawsite[c] and copy site[c] values to it */
      /* NOTE : applied to all maps including -subset maps (if -ScanScaling) */
      rawsitealloc(Gmap,rawsitemaps,totalmaps);
      rawsitemaps = totalmaps;
    }

    if(DEBUG>=2 && (usecolor || colors >= 2)){
      for(int i= 0; i < nummaps; i++){
	Cmap *pmap = Gmap[i];
	assert(fabs(pmap->rawsite[0][pmap->numsite[0]+1] - pmap->rawsite[1][pmap->numsite[1]+1]) < 1e-3);
      }
    }
  }

  if(DEBUG) assert(sizeof(RFLOAT) <= sizeof(FLOAT));

  if(VERB>=2){
    printf("PixelLen=%0.4f,origPixelLen=%0.4f\n",PixelLen,origPixelLen);
    fflush(stdout);
  }

  if(extend && !(ScoreThreshold >= 0.0 || LogPvThreshold >= 7.0 )){
    printf("WARNING:-ve score threshold=%e without -T 1e-7 (or better) not recommended with -extend (Since any nanomap can achieve a 0 score by extending the entire map beyond the reference)\n",ScoreThreshold);
    fflush(stdout);
  }

  /* adjust initial values based on bounds */
  for(int c = 0; c < colors; c++){
    if(FP[c] < MINFP){
      printf("FP[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,FP[c],MINFP,MAXFP,MINFP);
      FP[c] = MINFP;
    } else if(FP[c] > MAXFP){
      printf("FP[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,FP[c],MINFP,MAXFP,MAXFP);
      FP[c] = MAXFP;
    }
    if(FN[c] < MINFN){
      printf("FN[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,FN[c],MINFN,MAXFN,MINFN);
      FN[c] = MINFN;
    } else if(FN[c] > MAXFN){
      printf("FN[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,FN[c],MINFN,MAXFN,MAXFN);
      FN[c] = MAXFN;
    }
    if(SF[c] < MINSF){
      printf("SF[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,SF[c],MINSF,MAXSF,MINSF);
      SF[c] = MINSF;
    } else if(SF[c] > MAXSF){
      printf("SF[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,SF[c],MINSF,MAXSF,MAXSF);
      SF[c] = MAXSF;
    }
    if(SR[c] < max(0.0,MINSR)){
      printf("SR[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,SR[c],MINSR,MAXSR,MINSR);
      SR[c] = max(0.0,MINSR);
    } else if(SR[c] > MAXSR){
      printf("SR[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,SR[c],MINSR,MAXSR,MAXSR);
      SR[c] = MAXSR;
    }
    if(SD[c] < MINSD){
      printf("SD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,SD[c],MINSD,MAXSD,MINSD);
      SD[c] = MINSD;
    } else if(SD[c] > MAXSD){
      printf("SD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,SD[c],MINSD,MAXSD,MAXSD);
      SD[c] = MAXSD;
    }

    if((SR[c] < 0.0 || SD[c] < -sqrt(2.0*SF[c]*SR[c])) * MINSD_MULT > 0.0){
      if(SR[c] < 0.0)
	printf("ERROR:SD[%d]=%0.6f is not valid for SF=%0.6f,SR=%0.6f (must specify +ve sr with -ve sd)\n",c+1,SD[c],SF[c],SR[c]);
      else
	printf("WARNING:SD[%d]=%0.6f is not valid for SF=%0.6f,SR=%0.6f: changing to %0.6f\n",c+1,SD[c],SF[c],SR[c],-sqrt(2.0*SF[c]*SR[c])*MINSD_MULT);
      fflush(stdout);
      assert(SR[c] >= 0.0);
      SD[c] = -sqrt(2.0*SF[c]*SR[c]) * MINSD_MULT;
    }
    if(SE[c] < MINSE){
      printf("SE[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,SE[c],MINSE,MAXSE,MINSE);
      SE[c] = MINSE;
    } 
    if(SE[c] > MAXSE){
      printf("SE[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",c+1,SE[c],MINSE,MAXSE,MAXSE);
      SE[c] = MAXSE;
    }
    PRbias[c] = 0.0;
  }

  /* initialize pmap->nalign,align */
  for(int i= 0; i < startmaps; i++){
    Cmap *pmap = Gmap[i];
    for(int c = 0; c < colors; c++){
      pmap->align[c].reset();
      pmap->nalign[c].reset();
    }
  }

  /* first check the average query and reference interval size to see if FP and FN are in the right ballpark */
  for(int c = 0; c < colors; c++){
    double rawtheta = 0.0;
    int rawXsites = 0;
    for(int i= 0; i < startmaps; i++){
      Cmap *pmap = Gmap[i];
      FLOAT *X = pmap->site[c];
      int N = pmap->numsite[c];
      rawtheta += X[N+1];
      rawXsites += N;
    }
    if(rawXsites > 0)
      rawtheta /= rawXsites;
    double rawlambda = 0.0;
    int rawYsites = 0;
    for(int refid = 0;refid < numrefmaps;refid++){
      Cmap *rmap = refmap[refid];
      FLOAT *Y = rmap->site[c];
      int N = rmap->numsite[c];
      rawlambda += Y[N+1];
      rawYsites += N;
    }
    if(rawYsites > 0)
      rawlambda /= rawYsites;
    if(rawXsites > 0 && rawYsites > 0){/* check if FP and FN need to be adjusted upwards */
      if(rawtheta < rawlambda*0.5){
	printf("WARNING:color=%d:raw theta=%0.3f,lambda=%0.3f : theta is too small for a reasonable false positive rate\n",
	       c+1,rawtheta,rawlambda);
	fflush(stdout);
	//      exit(1);
      }
    }
  }

  int maxthreads = 1;
  #ifdef _OPENMP
  maxthreads = omp_get_max_threads();
  maxthreads = min(maxthreads,MaxThreads);
  #endif

  MaxPair = maxthreads * 1024 * 8 + nummaps;
  if(VERB && (VERB>=2 || HASH_DEBUG)){
    printf("MaxPair=%d,maxthreads=%d,nummaps=%d\n",MaxPair,maxthreads,nummaps);
    fflush(stdout);
  }
  if(DEBUG) assert(MaxPair > 0);/* make sure it hasn't wrapped around 32 bit signed number */

  if(hash_filename){
#if HASH_STREAM < 2
    printf("HASH_STREAM=%d no longer supported (see hash.h)\n",HASH_STREAM);
    fflush(stdout);
#endif

    if(idrenumbered && !HashGen){
      printf("hash table requires unique +ve id numbers for input Maps\n");
      exit(1);
    }
    if(NoSplit < 0){
      printf("hash table does not support -nosplit 0\n");
      exit(1);
    }

    if(!HashGen){ /* sort maps by id in ascending order and renumber Ids */
      qsort(refmap,numrefmaps,sizeof(Cmap *),(intcmp*)CmapIdInc);
      qsort(Gmap,nummaps,sizeof(Cmap *),(intcmp*)CmapIdInc);
    
      for(int m = 0; m < nummaps; m++)
	Gmap[m]->mapid = m;
      for(int m = 0; m < numrefmaps; m++)
	refmap[m]->mapid = m;
    }

    if(DEBUG) assert(hashpairs1 == 0 && numhashpairs1 == 0);

    size_t L = strlen(hash_filename);
    if(L >= strlen(".hashbin") && !strcmp(&hash_filename[L-strlen(".hashbin")],".hashbin")){/* binary file input */
      BinaryFormat = 1;
      
      /* set up streaming hashtable IO */
      /* set up input buffers */
      size_t num = maxMatchesPerRefid + 1;
      if(DEBUG) assert(maxMapsPerRefid <= MASK(31));
      if((size_t)MaxPair < maxMapsPerRefid)
	MaxPair = maxMapsPerRefid;
      if((size_t)MaxPair > num)
	num = MaxPair;
      if(VERB){
	printf("maxMatchesPerRefid=%lu, maxMapsPerRefid=%lu: num=%lu,MaxPair=%d\n",maxMatchesPerRefid,maxMapsPerRefid,num,MaxPair);
	fflush(stdout);
      }
      
      YidListMem = new int[MaxPair*4];
      YidList1 = &YidListMem[0];
      XidList1 = &YidListMem[MaxPair];
      phashListMem = new CHashMatch*[MaxPair*2];
      phashList1 = &phashListMem[0];
      NumPair1 = 0;

      YidList2 = &YidListMem[MaxPair*2];
      XidList2 = &YidListMem[MaxPair*3];
      phashList2 = &phashListMem[MaxPair];
      NumPair2 = 0;

      numhashpairs1 = maxhashpairs1 = 0;
      maxhashalloc(num+1,numhashpairs1,maxhashpairs1,hashpairs1);// Allocate one extra space so phash[1] for last phash == &hashpairs1[maxhashpairs1-1] is valid
      nexthash1 = hashpairs1;

      numhashpairs2 = maxhashpairs2 = 0;
      maxhashalloc(num+1,numhashpairs2,maxhashpairs2,hashpairs2);
      nexthash2 = hashpairs2;

      if(VERB>=2){
	printf("Allocated Hash Buffers:MaxPair=%d (%llu bytes for 2 buffers)\n",MaxPair,(unsigned long long)(MaxPair*2*(sizeof(CHashMatch) + sizeof(int)*2+sizeof(CHashMatch*))));
	fflush(stdout);
      }

      // The Following hash_open() will be called before every -M iteration
      // MatchMax = hash_open(HashFP, hash_filename);

    } else {/* text format no longer supported */
      printf("refalign_allpairs:Text format Hashtable in %s not supported\n",hash_filename);
      exit(1);
    }
  }

  logLRsum = LARGE_NEGATIVE;
  logLRsdsum = LARGE_NEGATIVE;
  totXlen = LARGE_NEGATIVE;
  logLRarray = NULL;
  for(int c = 0; c < colors; c++)
    logLRarrayC[c] = logSDarray[c] = Xlen[c] = NULL;

  if(RefRepeats < 1)
    RefRepeats = 1;

  Cparameters *parameters = new Cparameters[RefRepeats+1];

  numalign_start = new size_t[numrefmaps+1];
  numalign_end = new size_t[numrefmaps+1];
  orignummaps = new int[numrefmaps];

  double PixelLenSD = 0.0;

  if(Refine && RefRepeats > 1){
    printf("-M %d not compatible with -refine %d\n",RefRepeats, Refine);
    exit(1);
  }

  int kmax[MAXCOLOR];
  int sumN = 0;/* sum of N[c]+1 */
  for(int c = 0; c < colors; c++){
    maxN[c] = 0;
    for(int refid = 0; refid < numrefmaps; refid++){
      Cmap *rmap = refmap[refid];
      if(DEBUG) assert(rmap != 0);
      int N = rmap->numsite[c];
      sumN += N+1;
      if(N > maxN[c])
	maxN[c] = N;
    }
  }
  int ***RKKmax = new int**[numrefmaps];/* RKKmax[refid 0..numrefmaps-1][c=0..colors-1][1..N[c]] */
  int **RKmem = new int*[numrefmaps*colors];
  int *RKKmem = new int[sumN];
  for(int refid = 0; refid < numrefmaps; refid++)
    RKKmax[refid] = &RKmem[refid * colors];
  
  /* NOTE: since res & resSD can change, the following computation of maxdist,AlignResMax and Kmax[] should be repeated before each score_init() */
  double maxdist[MAXCOLOR],AlignResMax[MAXCOLOR];

  if(!(ResEstimate && !REFDEBUG_STRICT)){
    int cntN = 0;
    for(int c = 0; c < colors; c++){
      maxdist[c] = (res[c] + SDrange * resSD[c])*PixelLen;
      AlignResMax[c] = maxdist[c] * AlignResMult;
      kmax[c] = 0;
      
      long long ksum = 0, Nsum = 0;
      for(int refid = 0; refid < numrefmaps; refid++){
	Cmap *rmap = refmap[refid];
	FLOAT *Y = rmap->site[c];
	int N = rmap->numsite[c];

	if(DEBUG>=2){/* check that site Y[0..N+1] are monotonic */
	  for(int I = 0; I <= N; I++)
	    if(DEBUG && !(Y[I+1] >= Y[I])){
	      printf("refid=%d:c=%d,I=%d,N=%d:Y[I]=%0.8f,Y[I+1]=%0.8f\n",refid,c,I,N,Y[I],Y[I+1]);
	      fflush(stdout);
	      assert(Y[I+1] >= Y[I]);
	    }
	}

	int *Kmax = RKKmax[refid][c] = &RKKmem[cntN];   cntN += N+1;

	for(int I = 1; I <= N; I++){
	  int K = 0;
	  if(SCORE_APPROX >= 2 || !AlignRes){/* determine Kmax[I] based on total distance Y[I] - Y[I-K] < maxdist[c] */
	    while(I-K > 1){
	      if(Y[I] - Y[I-K-1] >= maxdist[c])
		break;
	      K++;
	    }
	  } else {
	    while(I-K > 1){
	      if(Y[I-K] - Y[I-K-1] >= maxdist[c] || Y[I]-Y[I-K-1] >= AlignResMax[c])
		break;
	      K++;
	    }
	  }
	  if(K > KMAX){
	    if(VERB >= 2){
	      printf("c=%d,refid=%d:Kmax[%d]=%d(reducing to %d):res[c]=%0.3f,resSD[c]=%0.3f:",c,refid,I,Kmax[I],KMAX, res[0],resSD[0]);
	      for(int U = 0; U <= K; U++)
		printf("Y[%d]=%0.3f,",I-U,Y[I-U]);
	      printf("\n");
	    }
	    K = KMAX;
	  }
	  Kmax[I] = K;
	  if(DEBUG>=2 && resSD[c] <= 0.001 && Kmax[I]>=2){
	    printf("Kmax[%d]=%d:c=%d,res[c]=%0.3f,resSD[c]=%0.3f:",I,Kmax[I],c,res[c],resSD[c]);
	    for(int U = 0; U <= K; U++)
	      printf("Y[%d]=%0.3f,",I-U,Y[I-U]);
	    printf("\n");
	  }
	}

	for(int I = 1; I <= N; I++){
	  ksum += Kmax[I];
	  if(Kmax[I] > kmax[c])
	    kmax[c] = Kmax[I];
	}
	Nsum += N;
      }
      if(DEBUG) assert(cntN <= sumN);
      if(VERB){
	printf("c=%d:kmax=%d,average Kmax[I]=%0.3f\n",
	       c, kmax[c], ((double) ksum)/Nsum);
	fflush(stdout);
      }
    }
  }

  if(Refine && LogPvThresholdTE < LogPvThreshold)
    LogPvThreshold = LogPvThresholdTE;

  for(int c = 0; c < colors;c++)
    parameters[0].FPfrac[c] = 0.0;

  double origScoreThreshold2 = ScoreThreshold;
  double origLogPvThreshold2 = LogPvThreshold;

  int *scaleIDcnt = new int[max(1,NumScaleFactor)];

  for(giter = 0; giter < RefRepeats; giter++){
    /* NOTE : if NoSplit < 2, each iterations will create a seperate set of map splits that accumulate at the end of Gmap[startmaps ... totalmaps-1] */

    if(DEBUG>=2){
      for(int mapid = 0; mapid < nummaps; mapid++){
        if(!(Gmap[mapid]->mapid == mapid)){
 	  printf("mapid=%d:Gmap[mapid]=%p,Gmap[mapid]->mapid=%d,Gmap[mapid]->id=%lld\n",
		mapid,Gmap[mapid],Gmap[mapid]->mapid,Gmap[mapid]->id);
	  for(int m = 0; m < 10; m++)
	    printf("m=%d:Gmap[m]=%p,Gmap[m]->mapid=%d,Gmap[m]->id=%lld\n",
		  m,Gmap[m],Gmap[m]->mapid,Gmap[m]->id);
	  fflush(stdout);
	  assert(Gmap[mapid]->mapid == mapid);
        }
      }
    }

    ScoreThreshold = origScoreThreshold2;
    LogPvThreshold = origLogPvThreshold2;/* restore original LogPvThreshold, in case it was change by -MapRate */

    /* save current parameter values */
#pragma novector
    for(int c = 0; c < colors;c++){
      parameters[giter].res[c] = res[c];
      parameters[giter].resSD[c] = resSD[c];

      parameters[giter].FP[c] = FP[c];
      parameters[giter].FN[c] = FN[c];
      parameters[giter].SF[c] = SF[c];
      parameters[giter].SD[c] = SD[c];
      if(QUADRATIC_VARIANCE)
	parameters[giter].SR[c] = SR[c];
      if(RES_VARIANCE)
	parameters[giter].SE[c] = SE[c];
      parameters[giter].minSNR[c] = minSNR[c];

      parameters[giter].ResBins[c] = ResBins[c];
      for(int Bin = 0; Bin <= ResBins[c]; Bin++){
	parameters[giter].resbias[c][Bin] = resbias[c][Bin];
	parameters[giter].resbiasX[c][Bin] = resbiasX[c][Bin];
	if(DEBUG/* HERE >=2 */) assert(maxresbias >= resbiasX[c][Bin] - 1e-6);	
      }
      if(DEBUG/* HERE >=2 */ && ResBins[c] > 0) assert(maxresbias >= resbiasX[c][ResBins[c]] - 1e-6);
    }
    parameters[giter].SF[colors] = SF[colors];
    parameters[giter].SD[colors] = SD[colors];
    parameters[giter].MA_mean = MA_mean;
    parameters[giter].PixelLen = PixelLen;
    parameters[giter].PixelLenSD = PixelLenSD;

    parameters[giter].mres = mres; // NOTE : mres does not change during autonoize
    parameters[giter].mresSD = mresSD; // NOTE : mresSD does not change during autonoise

    parameters[giter].maxresbias = maxresbias;

    printf("%d:", giter+1);
    fflush(stdout);

    numaligns = 0;
    nummaps = startmaps;
    for(int mapid = 0; mapid < nummaps; mapid++)
      Gmap[mapid]->origmap = 0;/* mark map as not a chimeric fragment */

    for(int c = 0; c <= colors; c++)
      numerrors[c] = 0;
    for(int c = 0; c < colors; c++)
      numresdata[c] = 0;

    /* initialize statistics counter */
    scoresum = 0.0;
    logPVsum = 0.0;
    NumpairsSum = 0.0;
    double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};

    ATscoresum = 0.0;
    ATlogPVsum = 0.0;
    ATNumpairsSum = 0.0;
    ATlogLRsum = 0.0;
    if(REFDEBUG && SDDEBUG)
      logLRsdsum = 0.0;

    int mapcnt = 0;/* Total number of maps that were aligned (over all reference maps : same map can align with multiple reference maps) */
    double resK2[MAXCOLOR];
    double lenY[MAXCOLOR];/* sum of aligned lengths including internal outlier regions */
    double len[MAXCOLOR];/* sum of aligned lengths (excluding internal outlier regions) */
    double lenB[MAXCOLOR];/* sum of aligned lengths - resK2 (excluding internal outlier regions) */
    int sitecnt[MAXCOLOR], fp[MAXCOLOR]; double fn[MAXCOLOR];/* count of total, false-negative and false-positive sites in 
								aligned region of reference */
    int fpI[MAXCOLOR], fnI[MAXCOLOR];/* false-positives and false-negative sites that are NOT within FP_DIST of true-positive site */

    int segcnt[MAXCOLOR];/* count of aligned intervals */
    double PrPen[MAXCOLOR];/* sum of penalties for log(Pr(y)) and Sm()-log(p) */

    for(int color = 0; color < colors; color++){
      resK2[color] = 2.0*res[color]*PixelLen;
      lenY[color] = len[color] = lenB[color] = 0.0;
      sitecnt[color] = fp[color] = 0; fn[color] = 0.0;
      fpI[color] = fnI[color] = 0;
      segcnt[color] = 0;
      PrPen[color] = 0.0;
    }

    int outliertot = 0;/* all intervals that are candidates for outliers (eg x < y) */
    int outliercnt = 0;
    int outlierEndcnt = 0;/* number of ends that were outliers (Lend or Rend <= -3) (per color) */
    double outliersum = 0.0;

    int chimcnt = 0;/* count of chimeric sites (x colors) at least AlignedSiteThreshold from either end (for best alignments) */
    int chimconf = 0;/* count of confirmed chimeric fragments (with good alignment to seperate location) */
    int chimconfFP = 0;/* count of confirmed chimeric fragments that are FP (if ground truth is known) */
    chimpaircnt = 0;/* count how many chimeric pairs (see chimconf) are next to each other in same orientation (within 20 * Ylambda) */

    size_t errcnt = 0, totcnt=0, tcnt =0;

    //    int stitchcnt = 0, stitchOut = 0;/* count of internal intervals with stitch and intervals with both a stitch and outlier (if stitch information is present) */
    //    int stitchEnd = 0, stitchEndOut = 0;/* count of end intervals with stitch and intervals with both a stitch and outlier (if stitch information is present) */
#if 0
    int sitecntS[MAXCOLOR], fpS[MAXCOLOR]; double fnS[MAXCOLOR],lenS[MAXCOLOR];/* FP and FN estimates for stitched intervals */
    for(int c = 0; c < colors; c++){
      sitecntS[c] = fpS[c] = 0;
      fnS[c] = lenS[c] = 0.0;
    }
#endif

    for(int i = 0; i < max(1,NumScaleFactor); i++)
      scaleIDcnt[i] = 0;

    if(ResEstimate && !REFDEBUG_STRICT){/* recompute maxdist,AlignResMax & Kmax[] since res,resSD may have changed */
      int cntN = 0;
      for(int c = 0; c < colors; c++){
	maxdist[c] = (res[c] + SDrange * resSD[c])*PixelLen;
	AlignResMax[c] = maxdist[c] * AlignResMult;
	kmax[c] = 0;

	long long ksum = 0, Nsum = 0;
	for(int refid = 0; refid < numrefmaps; refid++){
	  Cmap *rmap = refmap[refid];
	  FLOAT *Y = rmap->site[c];
	  int N = rmap->numsite[c];

	  if(DEBUG>=2){/* check that site Y[0..N+1] are monotonic */
	    for(int I = 0; I <= N; I++)
	      if(DEBUG && !(Y[I+1] >= Y[I])){
		printf("refid=%d:I=%d,N=%d:Y[I]=%0.8f,Y[I+1]=%0.8f\n",refid,I,N,Y[I],Y[I+1]);
		fflush(stdout);
		assert(Y[I+1] >= Y[I]);
	      }
	  }

	  int *Kmax = RKKmax[refid][c] = &RKKmem[cntN];   cntN += N+1;

	  for(int I = 1; I <= N; I++){
	    int K = 0;
	    if(SCORE_APPROX >= 2 || !AlignRes){/* determine Kmax[I] based on total distance Y[I] - Y[I-K] < maxdist[c] */
	      while(I-K > 1){
		if(Y[I] - Y[I-K-1] >= maxdist[c])
		  break;
		K++;
	      }
	    } else {
	      while(I-K > 1){
		if(Y[I-K] - Y[I-K-1] >= maxdist[c] || Y[I]-Y[I-K-1] >= AlignResMax[c])
		  break;
		K++;
	      }
	    }
	    if(K > KMAX){
	      if(VERB >= 2){
		printf("c=%d:Kmax[%d]=%d(reducing to %d):res[c]=%0.3f,resSD[c]=%0.3f:",c,I,Kmax[I],KMAX,res[c],resSD[c]);
		for(int U = 0; U <= K; U++)
		  printf("Y[%d]=%0.3f,",I-U,Y[I-U]);
		printf("\n");
	      }
	      K = KMAX;
	    }
	    Kmax[I] = K;
	    if(DEBUG>=2 && resSD[c] <= 0.001 && Kmax[I]>=2){
	      printf("c=%d:Kmax[%d]=%d:res[c]=%0.3f,resSD[c]=%0.3f:",c,I,Kmax[I],res[c],resSD[c]);
	      for(int U = 0; U <= K; U++)
		printf("Y[%d]=%0.3f,",I-U,Y[I-U]);
	      printf("\n");
	    }
	  }

	  for(int I = 1; I <= N; I++){
	    ksum += Kmax[I];
	    if(Kmax[I] > kmax[c])
	    kmax[c] = Kmax[I];
	  }
	  Nsum += N;
	}
	if(DEBUG) assert(cntN <= sumN);
	if(VERB){
	  printf("c=%d:kmax=%d,average Kmax[I]=%0.3f\n",
		 c,kmax[c], ((double) ksum)/Nsum);
	  fflush(stdout);
	}
      }
    }

    if(VERB>=2){
      printf("refalign:nummaps=%d,startmaps=%d,totalmaps=%d:Before score_init: cpu time=%0.6f, wall time=%0.6f\n",nummaps,startmaps,totalmaps,mtime(),wtime());
      fflush(stdout);
    }

    (void)score_init(refmap,0,numrefmaps,nummaps, RKKmax);

    if(VERB>=2){
      printf("refalign:nummaps=%d,startmaps=%d,totalmaps=%d:After score_init (mem=%0.1f Mbytes): cpu time=%0.6f, wall time=%0.6f\n",nummaps,startmaps, totalmaps, score_init_mem/(1024.0*1024.0), mtime(),wtime());
      fflush(stdout);
    }
    ChimScore = OutlierEndBias + OutlierEndPenalty;

    if(nexthash1){/* go back to start of HashTable data */
      if(VERB>=2){
	printf("Calling hash_open:cpu time=%0.6f, wall time=%0.6f\n",mtime(),wtime());
	fflush(stdout);
      }
      if(HashFP){
	(void)hash_close(HashFP,0);
	numhashpairs1 = numhashpairs2 = 0;
	NumPair1 = NumPair2 = 0;
	nexthash1 = hashpairs1;
	refid1 = refidend1 = -1;
	if(HASH_STREAM>=2){
	  nexthash2 = hashpairs2;
	  refid2 = refidend2 = -1;
	}
	HashFP = NULL;
      }
      MatchMax = hash_open(HashFP, hash_filename);
      if(DEBUG) assert(HashFP != NULL);
      if(VERB>=2){
	printf("   Finished hash_open:MatchMax=%llu:cpu time=%0.6f, wall time=%0.6f\n",(unsigned long long)MatchMax,mtime(),wtime());
	printf("nexthash1=%p,FirstAlignments=%d,NoSplit=%d\n",nexthash1,FirstAlignments,NoSplit);
	fflush(stdout);
      }
      totalhashpairs = 0;
    }

    if(DEBUG>=2){
      for(int mapid = 0; mapid < nummaps; mapid++){
        if(!(Gmap[mapid]->mapid == mapid)){
 	  printf("mapid=%d:Gmap[mapid]=%p,Gmap[mapid]->mapid=%d,Gmap[mapid]->id=%lld\n",
		mapid,Gmap[mapid],Gmap[mapid]->mapid,Gmap[mapid]->id);
	  for(int m = 0; m < 10; m++)
	    printf("m=%d:Gmap[m]=%p,Gmap[m]->mapid=%d,Gmap[m]->id=%lld\n",
		  m,Gmap[m],Gmap[m]->mapid,Gmap[m]->id);
	  fflush(stdout);
	  assert(Gmap[mapid]->mapid == mapid);
        }
      }
    }

    if(DEBUG>=2) assert(0 <= giter && giter <= RefRepeats-1);

    if(nexthash1){/* with hashtable */
      if(!FLAT_MEMORY){
	printf("refalign2.cpp: No support for hashtable without FLAT_MEMORY\n");
	exit(1);
      }

      for(int mapid = 0;mapid < nummaps;mapid++){
	Cmap *nanomap = Gmap[mapid];
	nanomap->incscale = 1.0;
	nanomap->incwt = 0.0;
      }

      int refidend = -1;
   
      for(int refid = 0; refid < numrefmaps; refid = refidend + 1){

	refidend = refidend1;
	if(refid != refid1){ /* Compute YidList1[],XidList1[],phashList1[] (may just copy/swap from 2nd buffer set) */
	  refidend = loadhash2(refid, nummaps, nexthash1, hashpairs1, numhashpairs1, maxhashpairs1, YidList1, XidList1, phashList1, NumPair1);
	  if(DEBUG) assert(refid==refid1);
	}
	if(VERB>=2){
	  printf("refid=%d(refid1=%d),refidend=%d(refidend1=%d), NumPair1=%d\n",
		 refid,refid1,refidend,refidend1,NumPair1);
	  fflush(stdout);
	}
	if(refidend < refid){
	  for(int rid = refid; rid <= numrefmaps; rid++)
	    numalign_start[rid] = numalign_end[rid] = numaligns;	    
	  break;
	}
	if(NumPair1 <= 0) {
	  for(int rid = refid; rid <= refidend; rid++)
	    numalign_start[rid] = numalign_end[rid] = numaligns;	    
	  continue;
	}

	if(DEBUG) assert(YidList1[0] >= refid && YidList1[NumPair1-1] <= refidend);

	/* determine pcnt[c] = sum(I=1..NN[c]) KKmax[c][I]+1 for largest reference map in range refid .. refidend */	
	long long limN[MAXCOLOR];
	int alignblock = (NoSplit <= 1 && HashMultiMatch) ? 3 : 1;
	size_t orignumaligns = numaligns;
	numaligns += NumPair1 * alignblock;

	for(int c = 0; c < colors; c++){
	  pcnt[c] = 0;
	  limN[c] = 0;
	  for(int rid = refid; rid <= refidend; rid++){
	    Cmap *rmap = refmap[rid];
	    int N = rmap->numsite[c];
	    int *Kmax = RKKmax[rid][c];
	    if(N > limN[c])
	       limN[c] = N;
	    
	    long long pcnt_rid = 0;
	    for(int I=1; I <= N; I++)
	      pcnt_rid += Kmax[I] +1;
	    pcnt[c] = max(pcnt[c], pcnt_rid);
	  }
	  if(DEBUG) assert(pcnt[c] >= 0 && pcnt[c] <= MAXINT);
	}

	if(VERB>=2){
	  printf("refid=%d..%d,nummaps=%d:numaligns=%llu->%llu,NumPair1=%d,refid1=%d..%d,refid2=%d..%d\n",refid,refidend,nummaps,
		 (unsigned long long)orignumaligns,(unsigned long long)numaligns,NumPair1,refid1,refidend1,refid2,refidend2);
	  fflush(stdout);
	}

	maxalignallocNULL(numaligns,alignment,numaligns,maxaligns);

	/* determine number of threads (numthreads) to use */
	int numthreads = 1;
        #ifdef _OPENMP
	numthreads = omp_get_max_threads();
	numthreads = min(numthreads,MaxThreads);
	if(MaxMem > 0 && pcnt[0]+pcnt[1] > 0 && maxM[0]+maxM[1] > 0){/* NOTE : maxM, the largest map size, is computed in score_init() */
	  size_t DPsiz = DPsizF * sizeof(TFLOAT) + DPsizI * sizeof(int);
	  double virtperthread = 0.0, memperthread = 0.0;
	  for(int c = 0; c < colors; c++){
	    virtperthread += (1.0 + kmax[c]) * limN[c] * ( maxM[c] * DPsiz) + sizeof(AL2);// virtual memory used per thread
	    memperthread += (double)pcnt[c] * (double)(maxM[c]  * DPsiz) + sizeof(AL2);// Estimate of real memory used per thread
	  }
	  if(MIN_MEM && DIAGONAL>=2 && hashdelta){ /* reduce estimated virtual and real memory by 4x */
	    virtperthread *= 0.25;
	    memperthread *= 0.25;
	  }
	  extern int numbufs;
	  double HashMem = MaxPair*2*(sizeof(CHashMatch) + sizeof(int)*2+sizeof(CHashMatch*)) + numbufs * HASHBUF_SIZ * 2 * sizeof(CHashMatch);
	  double maxmem = (double)MaxMem * 1024.0 * 1024.0 * 1024.0 - score_init_mem - HashMem;
	  int limitthreads = max(1,(int)floor(maxmem/virtperthread + 0.5));
	  //	  int limitthreads = max(1,(int)floor(maxmem/memperthread + 0.5));
	  if(limitthreads < numthreads){
	    if(VERB){
	      printf("refid=%d..%d: number of threads reduced from %d to %d due to MaxMem=%0.3f Gbytes, Score=%0.3f Gbytes, HashBuf=%0.3f Gbytes, per thread virtual memory = %0.3f Gbytes, real memory ~ %0.3f Gbytes:N<=%lld,%lld,pcnt=%lld,%lld(maxN=%lld,%lld,kmax=%d,%d),maxM=%lld,%lld\n",
		     refid,refidend, numthreads,limitthreads, MaxMem, score_init_mem/(1024.0*1024.0*1024.0), HashMem/(1024.0*1024.0*1024.0), 
		     virtperthread/(1024.0*1024.0*1024.0),memperthread/(1024.0*1024.0*1024.0),limN[0],limN[1],pcnt[0],pcnt[1],maxN[0],maxN[1],kmax[0],kmax[1],maxM[0],maxM[1]);
	      fflush(stdout);
	    }
	    if(DEBUG) assert(limitthreads >= 1);
	    numthreads = limitthreads;
	  } else if(VERB>=2){
	    printf("refid=%d..%d: MaxMem=%0.3f Gbytes, Score=%0.3f Gbytes, HashBuf=%0.3f Gbytes, per thread virtual = %0.3f Gbytes, real memory ~ %0.3f Gbytes:N<=%lld,%lld,pcnt=%lld,%lld(maxN=%lld,%lld,kmax=%d,%d),maxM=%lld,%lld\n",
		   refid,refidend, MaxMem, score_init_mem/(1024.0*1024.0*1024.0), HashMem/(1024.0*1024.0*1024.0),
		   virtperthread/(1024.0*1024.0*1024.0),memperthread/(1024.0*1024.0*1024.0),limN[0],limN[1],pcnt[0],pcnt[1],maxN[0],maxN[1],kmax[0],kmax[1],maxM[0],maxM[1]);
	    fflush(stdout);
	  }
	}
	if(numthreads > NumPair1)
	  numthreads = max(1,NumPair1);

        #endif

	int block = min(16,max(1,NumPair1/(64*numthreads)));
	if(VERB && numthreads > 1 && (VERB>=2 || numthreads != thread_printed)){
	  printf("refid=%d..%d:Using %d threads(NumPair1=%d,block=%d)\n",refid,refidend,numthreads,NumPair1,block);
	  fflush(stdout);
	  thread_printed = numthreads;
	}

	if(DEBUG) assert(numthreads >= 1);
	if(VERB>=2){
	  printf("refid=%d..%d/%d:kmax=%d,%d,maxN=%lld,%lld,maxM=%lld,%lld\n",refid,refidend,numrefmaps,kmax[0],kmax[1],maxN[0],maxN[1],maxM[0],maxM[1]);
	  fflush(stdout);
	}
	size_t totmem = 0;
	size_t DPsiz = DPsizF * sizeof(TFLOAT) + DPsizI * sizeof(int);
	for(int c = 0; c < colors; c++){
	  pcnt[c] = (1LL + kmax[c]) * maxN[c];
	  if(MIN_MEM && DIAGONAL>=2 && hashdelta)/* (re)allocate most memory in refalignXYsd() using malloc()/free() */
	    totmem += (maxN[c] + maxM[c]) * (2LL*sizeof(int *));
	  else
	    totmem += (maxN[c] + maxM[c]) * (2LL*sizeof(int *)) + pcnt[c] * maxM[c] * DPsiz;
	}
	if(MIN_MEM && DIAGONAL>=2 && hashdelta)/* (re)allocate most memory in refalignXYsd() using malloc()/free() */
	  rs_heap->reset_lightweight_heap(totmem +/*alignment fudge*/ 2*(4+1)*64, numthreads);
	else
	  rs_heap->reset_lightweight_heap(totmem +/*alignment fudge*/ 2*(6+1)*64, numthreads);

	if(VERB>=2){
	  printf("refid=%d..%d/%d: start of parallel section. CPU time=%0.6f wall time=%0.6f\n",refid,refidend,numrefmaps,mtime(),wtime());
	  fflush(stdout);
	}

	int lastPairId = 0;

	/* Multi-threaded work section : alignments of YidList1[0..NumPair1-1] with XidList1[0..NumPair1-1] */
	#pragma omp parallel num_threads(numthreads) if(numthreads > 1)
	{
	  int tid = 0;
          #ifdef _OPENMP
	  tid = omp_get_thread_num ();
          #endif
	  
	  /* each thread needs its own score table memory and Qmin[],Qmax[] */
	  CXPen XPenMem;
	  CXPen *XPen = &XPenMem;
	  int *Qmin = 0, *Qmax = 0;
	  XPen_alloc(XPen,Qmin,Qmax);

	  /* allocate DP array A[c=0..1].<field>(I=1..maxN[c], K=0..kmax[c], J=1..maxM[c]) */	  
	  AL2 A[2];
	  TFLOAT *Fmem[2];
	  int *Imem[2], *JMIN[2],*JMAX[2],*Imin[2],*Imax[2];
	  long long StrideMax[2] = {0,0};

	  for(int c = 0; c < colors; c++){
	    A[c].kmax = kmax[c];
	    if(sizeof(A[c].Kstride) < sizeof(long long) && (long long)maxN[c] * (long long)maxM[c] * (long long) kmax[c] > MASK(31)){
	      printf("c=%d:maxN=%lld,maxM=%lld,kmax=%d: maxN * maxM * kmax exceeds signed integer range : change Kstride to long long in refalign.cpp\n",c,maxN[c],maxM[c],A[c].kmax);
	      exit(1);
	    }
	    if(sizeof(A[c].Istride) < sizeof(long long) && (long long)maxM[c] > MASK(31)){
	      printf("c=%d:maxM=%lld: maxM exceeds signed integer range : change Istride to long long in refalign.cpp\n",c,maxM[c]);
	      exit(1);
	    }
	    
	    JMIN[c] = rs_heap->alloc_int(maxN[c],tid) - 1;
	    JMAX[c] = rs_heap->alloc_int(maxN[c],tid) - 1;
	    Imin[c] = rs_heap->alloc_int(maxM[c],tid) - 1;
	    Imax[c] = rs_heap->alloc_int(maxM[c],tid) - 1;
	    
	    StrideMax[c] = 0;
	    Fmem[c] = NULL;
	    Imem[c] = NULL;
	    if(!(MIN_MEM && DIAGONAL>=2 && hashdelta)){
	      StrideMax[c] = (1LL + kmax[c]) * maxN[c] * maxM[c];
	      Fmem[c] = rs_heap->alloc_TFLOAT(maxN[c] * (1LL + kmax[c]) * maxM[c] * DPsizF, tid);
	      Imem[c] = rs_heap->alloc_int(maxN[c] * (1LL + kmax[c]) * maxM[c] * DPsizI, tid);
	    }

	    if(VERB>=2){
              #pragma omp critical
	      {
		printf("thread=%d/%d,c=%d:memory allocated:block size=%d maps, pcnt=%lld,N<=%lld,kmax=%d,maxM=%lld,maxN=%lld, NumPair1=%d\n",tid,numthreads,c,block,pcnt[c],limN[c],kmax[c],maxM[c],maxN[c],NumPair1);
		fflush(stdout);
	      }
	    }
	  } // c = 0.. colors-1

	  int Tmapcnt = 0;

	  #pragma omp for schedule(dynamic,block)
	  for(int PairId = 0; PairId < NumPair1; PairId++){
	    CHashMatch *phash = phashList1[PairId];
	    int rid = YidList1[PairId];
	    int mapid = XidList1[PairId];
	    if(DEBUG>=2 && !(rid >= refid && rid <= refidend)){
	      #pragma omp critical
	      {
		printf("PaidId=%d/%d:rid=%d/%d,mapid=%d/%d,refid=%d..%d,refid1=%d..%d\n",
		       PairId,NumPair1,rid,numrefmaps,mapid,nummaps,refid,refidend,refid1,refidend1);
		fflush(stdout);
		assert(rid >= refid && rid <= refidend);
	      }
	    }

	    if(!tid && refidend+1 < numrefmaps && refid2 < refidend+1){
	      /* Read ahead into 2nd buffer in the master thread. In Non-parallel section swap the 2 buffers instead of reading data */
	      int end2 = loadhash2(refidend+1, nummaps, nexthash2, hashpairs2, numhashpairs2, maxhashpairs2, YidList2, XidList2, phashList2, NumPair2);
	      if(DEBUG) assert(refid2 >= refidend+1);
	      if(DEBUG) assert(refidend2 == end2);
	    }

	    size_t alignid = orignumaligns + PairId * alignblock;
	    size_t origalignid = alignid;
	    if(DEBUG/* HERE >=2 */ && !(alignid + alignblock <= numaligns)){
	      #pragma omp critical
	      {
		printf("PairId=%d/%d:rid=%d/%d,mapid=%d/%d:alignid=%llu,alignblock=%d,numaligns=%llu,numalign_start[rid]=%llu\n",
		       PairId,NumPair1,rid,numrefmaps,mapid,nummaps,(unsigned long long)alignid,alignblock,(unsigned long long)numaligns,(unsigned long long)numalign_start[rid]);
		fflush(stdout);
		assert(alignid < numaligns);
	      }
	    }
	    
	    for(int t = 0; t < alignblock; t++){
	      Calign *p = alignment[alignid+t];
	      if(p){
		p->allfree();
		for(int c = 1; c <= colors; c++)
		  p[c].allfree();
		if(ALIGN_COMPRESS){
		  delete [] p;
		  alignment[alignid+t] = NULL;
		}
	      }
	    }

	    if(FirstAlignments && mapcnt+Tmapcnt >= FirstAlignments)
	      continue;	  // NOTE : multithreaded loop is not determininistic since mapcnt is not synchronized

	    Cmap *rmap = refmap[rid];
	    int *NN = rmap->numsite;

	    Cmap *nanomap = Gmap[mapid];
	    if(DEBUG>=2) assert(nanomap != 0);
	    if(DEBUG>=2) assert(nanomap->mapid == mapid);
	    int *MM = nanomap->numsite;
	    //	    FLOAT **XX = nanomap->site;

	    if(VERB>=2 && PairId > lastPairId + 1000){
	      #pragma omp critical
	      {
		if(PairId > lastPairId + 1000){
		  printf("tid=%d:PairId=%d/%d:rid=%d,mapid=%d,N=%d,%d,M=%d,%d:Calling refalignSD\n",tid,PairId,NumPair1,rid,mapid,NN[0],NN[1],MM[0],MM[1]);
		  fflush(stdout);
		  lastPairId = PairId;
		}
	      }
	    }

	    int **KKmax = RKKmax[rid];
	    if(DEBUG>=2) assert(KKmax != 0);

	    refalignSD(rmap, nanomap, alignid, KKmax, rid, mapid, phash, XPen, Qmin, Qmax, A, Fmem, Imem, StrideMax, JMIN, JMAX, Imin, Imax, tid, maxN, NoSplit);

	    if(DEBUG>=2) assert(alignid <= origalignid + alignblock);

	    for(int t = 0; t < alignblock; t++){
	      Calign * &align = alignment[origalignid + t];

	      if(align){
		if(DEBUG) assert(align->mapid1 == rid);
		if(DEBUG && t==0) assert(align->mapid2 == mapid);
		if(DEBUG) assert(Gmap[align->mapid2]->id == Gmap[mapid]->id);
		
		if(AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold)){
		  if(t == 0){
		    Cmap *origmap = nanomap;
		    while(origmap->origmap)
		      origmap = origmap->origmap;
		    if(origmap->nalign->mapid1 < 0 || (BESTREF_PV ? (align->logPV > origmap->nalign->logPV + BestRefMargin * (rid > origmap->nalign->mapid1 ? 1.0 : -1.0)) :
						       (align->score > origmap->nalign->score + BestRefMargin * (rid > origmap->nalign->mapid1 ? 1.0 : -1.0)))){
		      /* save best alignment for next iteration */
		      for(int c = 0; c < colors; c++)
			origmap->nalign[c].update(align[c+1], NN[c]);
		      /* replace score & logPV of color 0 by total score & logPV */
		      origmap->nalign->score = align->score;
		      origmap->nalign->logPV = align->logPV;
		    }
		  }

		  Tmapcnt++;

		  if(FirstAlignments && Tmapcnt >= 16){
                    #pragma omp critical
		    {
		      mapcnt += Tmapcnt;
		      Tmapcnt = 0;
		    }
		  }
		}
	      }

	      if(VERB>=2 && nanomap->id == MAP_TRACE /* && giter == RefRepeats-1 */){
                #pragma omp critical
		{
		  if(!align || align->numpairs <= 0){
		    if(VERB>=2 && t==0)
		      printf("refid=%d,mapid=%d(id=%lld),M=%d,%d:no alignment found(thread=%d),t=%d\n",
			     rid,mapid,Gmap[mapid]->id,MM[0],MM[1],tid,t);
		  } else if (!AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold)){/* no good alignment */
		    if(VERB>=2 && t==0)
		      printf("refid=%d,mapid=%d(id=%lld),M=%d,%d:no good alignment:score=%0.4f,logPV=%0.3f,orientation=%d,scale=%d(%0.4f),pairs=%d(thread=%d),t=%d\n",
			     rid,mapid,Gmap[mapid]->id, MM[0],MM[1],align->score,align->logPV, align->orientation, align->scaleID, align->scaleID ? ScaleFactor[align->scaleID] : 1.0,
			     align->numpairs,tid,t);
		  } else {
		    printf("refid=%d,mapid=%d,M=%d,%d:score=%0.4f,logPV=%0.3f,orientation=%d,scale=%d(%0.4f),pairs=%d,(thread=%d),mapcnt=%d,t=%d,repeat=%d\n",
			   rid,mapid,MM[0],MM[1],align->score,align->logPV, align->orientation, align->scaleID, align->scaleID ? ScaleFactor[align->scaleID] : 1.0, align->numpairs,tid,mapcnt,t,align->repeat);
		    if(SecondBest && align->align2){
		      Calign *align2 = align->align2;
		      printf("    2nd Best:refid=%d,mapid=%d(id=%lld),or=%d:score=%0.4f,logPV=%0.3f,pairs=%d\n",
			     align2->mapid1,align2->mapid2,Gmap[align2->mapid2]->id,align2->orientation,align2->score,align2->logPV,align2->numpairs);
		    }
		  }
		  fflush(stdout);
		}
	      } // pragma omp critical
	    }
	  } /* parallel for mapid = 0 .. nummaps -1 */

          #pragma omp critical
	  {
	    mapcnt += Tmapcnt;
	    Tmapcnt=0;
	  }

	  /* free per-thread memory */
	  if(MIN_MEM && DIAGONAL>=2 && hashdelta){
	    for(int c = 0; c < colors; c++){/* NOTE : Fmem[c][] includes Imem[c] */
	      free(Fmem[c]); 
	      Fmem[c] = NULL;
	      Imem[c] = NULL;
	    }
	  }
	  XPen_free(XPen,Qmin,Qmax);
	} // parallel

	if(VERB){
	  printf("[refid=%d..%d/%d:total alignments=%d (alignids=%llu..%llu,nummaps=%d): cum CPU time=%0.6f wall time=%0.6f secs\n",
		 refid,refidend,numrefmaps,mapcnt,(unsigned long long)orignumaligns,(unsigned long long)numaligns,nummaps,mtime(),wtime());
	  fflush(stdout);
	}

	if(DEBUG) NumPair1 = 0;

	if(ALIGN_COMPRESS){/* compress alignment[numalign_start[refid]..numalign_end[refidend]-1] to remove NULL pointers */
	  size_t j = numalign_start[refid] = orignumaligns;
	  int rid = refid;
	  size_t k = j;
	  for(; k < numaligns; k++){
	    if(!alignment[k]){
	      k++;
	      break;
	    }
	    for(; alignment[j]->mapid1 > rid; rid++)
	      numalign_start[rid+1] = numalign_end[rid] = j;
	    if(VERB>=3 && giter == RefRepeats - 1){
	      Calign *p = alignment[j];
	      printf("alignment[%lu]=%p:refid=%d,mapid=%d,or=%d:numpairs=%d,score=%0.6f,logPV=%0.2f\n",
		     j,p,p->mapid1,p->mapid2,p->orientation,p->numpairs,p->score,p->logPV);
	      fflush(stdout);
	    }
	    j++;
	  }
	  for(; k < numaligns; k++){
	    if(!alignment[k])
	      continue;
	    if(DEBUG>=2) assert(j < k);
	    if(DEBUG>=2) assert(alignment[j]==0);
	    alignment[j] = alignment[k];
	    alignment[k] = 0;
	    for(; alignment[j]->mapid1 > rid; rid++)
	      numalign_start[rid+1] = numalign_end[rid] = j;
	    if(VERB>=3 && giter == RefRepeats - 1){
	      Calign *p = alignment[j];
	      printf("alignment[%lu]=%p:refid=%d,mapid=%d,or=%d:numpairs=%d,score=%0.6f,logPV=%0.2f\n",
		     j,p,p->mapid1,p->mapid2,p->orientation,p->numpairs,p->score,p->logPV);
	      fflush(stdout);
	    }
	    j++;
	  }
	  if(VERB>=2 && j < numaligns){
	    printf("refid=%d..%d:reduced numaligns from %llu to %llu(rid=%d): CPU time=%0.6f wall time=%0.6f\n",
		   refid,refidend,(unsigned long long)numaligns,(unsigned long long)j,rid,mtime(),wtime());
	    fflush(stdout);
	  }
	  numalign_end[rid] = numaligns = j;
	  while(++rid < numrefmaps)
	    numalign_start[rid] = numalign_end[rid] = numaligns;
	  if(VERB>=2){
	    printf("After removing null alignments from refid=%d..%d (orignumaligns=%llu)\n",refid,refidend,(unsigned long long)orignumaligns);
	    for(int r = refid; r <= refidend; r++)
	      printf("   refid=%d:alignid=%llu..%llu\n",r, (unsigned long long)numalign_start[r],(unsigned long long)numalign_end[r]-1);
	    fflush(stdout);
	  }
	}
      }
    } else { /* No hashtable */

      if(DEBUG) assert(!HashMultiMatch);
      if(DEBUG) assert(!nexthash1);
      if(DEBUG) assert(!HashFP);
      if(DEBUG) assert(nummaps == startmaps);

      for(int refid = 0;refid < numrefmaps;refid++){
	Cmap *rmap = refmap[refid];
	//	FLOAT **YY = rmap->site;
	int *NN = rmap->numsite;
	int **KKmax = RKKmax[refid];
	if(DEBUG>=2) assert(KKmax != 0);

	orignummaps[refid] = totalmaps;/* skip maps from nummaps..orignummaps[refid]-1 for current refid (Not needed if NoSplit >= 2)*/
	numalign_start[refid] = numalign_end[refid] = numaligns;

	if(VERB && numrefmaps > 1){
	  printf("refid=%d(id=%lld),sites=%d+%d:\n",refid,refmap[refid]->id,NN[0],NN[1]);
	  fflush(stdout);
	}
	if(VERB>=2){
	  printf("iter=%d:refid=%d:KKmax=x%p\n",giter,refid,KKmax);
	  fflush(stdout);
	}

	numaligns += nummaps;
	maxalignallocNULL(numaligns,alignment,numaligns,maxaligns);
	numalign_end[refid] = numaligns;

	if(DEBUG/* HERE >=2 */)
	  for(int c = 0; c < colors; c++)
	    assert(resSD[c] > 0.0 && KKmax[c] != 0);

	/* determine pcnt[c] = sum(I=1..NN[c]) KKmax[c][I]+1 */
	/* NOTE: With FLAT_MEMORY pcnt is only used to estimate real memory and adjusted the number of threads */
	long long limN[MAXCOLOR];
	for(int c = 0; c < colors; c++){
	  pcnt[c] = 0;

	  int *Kmax = KKmax[c];
	  for(int I=1; I <= NN[c]; I++)
	    pcnt[c] += Kmax[I]+1;
	  if(DEBUG) assert(pcnt[c] >= 0 && pcnt[c] <= MAXINT);

	  limN[c] = MINMEM ? NN[c] : maxN[c];
	}

	/* determine number of threads (numthreads) to use */
	int numthreads = 1;
        #ifdef _OPENMP
	numthreads = omp_get_max_threads();
	numthreads = min(numthreads,MaxThreads);
#if FLAT_MEMORY
	if(MaxMem > 0 && pcnt[0]+pcnt[1] > 0 && maxM[0]+maxM[1] > 0){ /* NOTE : maxM[], the largest map size, is computed in score_init() */
	  size_t DPsiz = DPsizF * sizeof(TFLOAT) + DPsizI * sizeof(int);
	  double virtperthread = 0.0, memperthread = 0.0;
	  for(int c = 0; c < colors; c++){
	    virtperthread += (1.0 + kmax[c]) * limN[c] * ( maxM[c] * DPsiz) + sizeof(AL2);// virtual memory used per thread
	    memperthread += (double)pcnt[c] * (double)(maxM[c]  * DPsiz) + sizeof(AL2);// Estimate of real memory used per thread
	  }
	  double maxmem = (double)MaxMem * 1024.0 * 1024.0 * 1024.0 - score_init_mem;
	  int limitthreads = max(1,(int)floor(maxmem/virtperthread + 0.5));
	  //	  int limitthreads = max(1,(int)floor(maxmem/memperthread + 0.5));
	  if(limitthreads < numthreads){
	    if(VERB){
	      printf("refid=%d(id=%lld): number of threads reduced from %d to %d due to MaxMem=%0.1f Gbytes, Score=%0.3f Gbytes, per threads virtual memory= %0.3f Gbytes, real memory = %0.3f Gbytes\n",
		     refid,refmap[refid]->id, numthreads,limitthreads, MaxMem, score_init_mem/(1024.0*1024.0*1024.0),
		     virtperthread/(1024.0*1024.0*1024.0), memperthread/(1024.0*1024.0*1024.0));
	      fflush(stdout);
	    }
	    numthreads = limitthreads;
	  }
	}
#else // FLAT_MEMORY == 0
	if(MaxMem > 0 && pcnt[0]+pcnt[1] > 0 && maxM[0]+maxM[1] > 0){ /* NOTE : maxM[], the largest map size, is computed in score_init() */
	  double memperthread = 0;
	  for(int c = 0; c < colors; c++)
	    memperthread += (double)pcnt[c] * (double)(maxM[c]  * sizeof(DP) + sizeof(DP*));
	  double maxmem = (double)MaxMem * 1024.0 * 1024.0 * 1024.0;
	  int limitthreads = (int)max(1.0,floor(maxmem/memperthread + 0.5));
	  if(limitthreads < numthreads){
	    if(VERB){
	      printf("refid=%d: number of threads reduced from %d to %d due to MaxMem=%0.1f Gbytes, per thread memory = %0.3f Gbytes\n",
		     refid,numthreads,limitthreads, MaxMem, memperthread/(1024.0*1024.0*1024.0));
	      fflush(stdout);
	    }
	    numthreads = limitthreads;
	  }
	}
#endif
	if(numthreads > nummaps)
	  numthreads = max(1,nummaps);
        #endif // OPEN_MP

	int block = min(3,max(1,nummaps/(64*numthreads)));
	if(VERB && numthreads > 1 && (VERB>=2 || numthreads != thread_printed)){
	  printf("refid=%d:Using %d threads(nummaps=%d,block=%d)\n",refid,numthreads,nummaps,block);
	  fflush(stdout);
	  thread_printed = numthreads;
	}

#if FLAT_MEMORY
	size_t totmem = 0;
	size_t DPsiz = DPsizF * sizeof(TFLOAT) + DPsizI * sizeof(int);
	for(int c = 0; c < colors; c++){
	  pcnt[c] = (1LL + kmax[c]) * limN[c];
	  totmem += (limN[c] + maxM[c]) * (2LL*sizeof(int *)) + pcnt[c] * maxM[c] * DPsiz;
	}
	int num_allocations = 2 * (6+1);
	rs_heap->reset_lightweight_heap(totmem + /* alignment fudge */ num_allocations * CACHE_LINE_GUARD*sizeof(float), numthreads);
#endif	

 	/* Multi-threaded work section : alignments of refid with mapid = 0..nummaps-1 */
	#pragma omp parallel num_threads(numthreads) 
	{
	  int tid = 0;
          #ifdef _OPENMP
	  tid = omp_get_thread_num ();
          #endif

	  /* each thread needs its own score table memory */
	  CXPen XPenMem;
	  CXPen *XPen = &XPenMem;
	  int *Qmin = 0, *Qmax = 0;
	  XPen_alloc(XPen,Qmin,Qmax);

#if FLAT_MEMORY == 0
	  /* allocate DP array AmemIKJ[c][I=1..NN[c]][K=0..RKKmax[refid][c][I]][J=1..maxM[c]] */
	  DP ***AmemIKJ[MAXCOLOR], **AmemIK[MAXCOLOR], *Amem[MAXCOLOR];
	  for(int c = 0; c < colors; c++){
	    int *Kmax = RKKmax[refid][c];
	    if(DEBUG) assert(maxM[c] <= MAXINT);
	    AmemIKJ[c] = new DP**[NN[c]+1];
	    AmemIK[c] = new DP*[pcnt[c]];
	    Amem[c] = new DP[pcnt[c] * maxM[c]];
	    int cnt = 0;
	    for(int I=1; I <= NN[c]; I++){
	      AmemIKJ[c][I] = &AmemIK[c][cnt];
	      cnt += Kmax[I]+1;
	    }
	    if(DEBUG) assert(cnt == pcnt[c]);
	  }
	  /* package AmemIKJ[c=0..colors-1] into A[c].array */
	  AL2 A[MAXCOLOR];
	  for(int c = 0; c < colors; c++)
	    A[c].array = AmemIKJ[c];

	  if(VERB>=2){
            #pragma omp critical
	    {
	      printf("thread=%d/%d:memory allocated:block size=%d maps, pcnt=%lld,%lld, maxM=%lld,%lld\n",
		     tid,numthreads,block,pcnt[0],pcnt[1],maxM[0],maxM[1]);
	      fflush(stdout);
	    }
	  }
#else // FLAT_MEMORY == 1
	  /* allocate DP array A[c=0..1].<field>(I=1..NN[c], K=0..RKKmax[refid][c][I], J=1..maxM[c]) */
	  AL2 A[2];
	  TFLOAT *Fmem[2];
	  int *Imem[2], *JMIN[2],*JMAX[2],*Imin[2],*Imax[2];
	  long long StrideMax[2] = {0,0};/* will not be used without hashtable */
	  for(int c = 0; c < colors; c++){
	    A[c].kmax = kmax[c];
	    if(sizeof(A[c].Kstride) < sizeof(long long) && (long long)limN[c] * (long long)maxM[c] * (long long) kmax[c] > MASK(31)){
	      printf("c=%d:limN=%lld,maxM=%lld,kmax=%d: limN * maxM * kmax exceeds signed integer range : change Kstride to long long in refalign.cpp\n",c,limN[c],maxM[c],kmax[c]);
	      exit(1);
	    }
	    if(sizeof(A[c].Istride) < sizeof(long long) && (long long)maxM[c] > MASK(31)){
	      printf("maxM=%lld: maxM exceeds signed integer range : change Istride to long long in refalign.cpp\n",maxM[c]);
	      exit(1);
	    }
	    Fmem[c] = rs_heap->alloc_TFLOAT(pcnt[c] * maxM[c] * DPsizF, tid);
	    Imem[c] = rs_heap->alloc_int(pcnt[c] * maxM[c] * DPsizI, tid);

	    JMIN[c] = rs_heap->alloc_int(limN[c],tid) - 1;
	    JMAX[c] = rs_heap->alloc_int(limN[c],tid) - 1;
	    Imin[c] = rs_heap->alloc_int(maxM[c],tid) - 1;
	    Imax[c] = rs_heap->alloc_int(maxM[c],tid) - 1;

	    if(VERB>=2){
              #pragma omp critical
	      {
		printf("thread=%d/%d,c=%d:memory allocated:block size=%d maps, pcnt=%lld,N=%d,kmax=%d,maxM=%lld,limN=%lld,Fmem=%p,Imem=%p,Imax=%p,end=%p\n",
		       tid,numthreads,c,block,pcnt[c],NN[c],kmax[c],maxM[c],limN[c],Fmem[c],Imem[c],Imax[c],Imax[c]+maxM[c]*sizeof(int));
		fflush(stdout);
	      }
	    }
	  }// c = 0 .. colors-1

#endif

	  int Tmapcnt = 0;

	  /* NOTE: If NoSplit <= 1 : loop must be broken into two parts since loop bounds cannot change inside parallel loop */
	  if(DEBUG) assert(startmaps <= orignummaps[refid]);

          #pragma omp for schedule(dynamic,block)
	  for(int mapid = 0;mapid < startmaps;mapid++){
	    if(FirstAlignments && mapcnt + Tmapcnt >= FirstAlignments)
	      continue;	  // NOTE : multithreaded loop is not determininistic since mapcnt is not synchronized

	    size_t alignid = numalign_start[refid] + mapid;
	    if(DEBUG) assert(alignid < numaligns);
	    if(DEBUG) assert(mapid < orignummaps[refid]);
	    Calign * &align = alignment[alignid];// NOTE : pointer value at alignment[alignid] may change in refalignSD !

	    if(VERB && nummaps > 1000*numthreads && !(mapid%(1000*numthreads))){
              #pragma omp critical
	      {
		mapcnt += Tmapcnt;
		Tmapcnt = 0;
      		if(HASH_DEBUG)
		  printf("[refid=%d/%d,mapid=%d/%d,alignments=%d+]: CPU time=%0.6f, wall time=%0.6f secs\n",refid,numrefmaps,mapid,nummaps,mapcnt, mtime(),wtime());
		else
		  printf("[refid=%d/%d,mapid=%d/%d]\n",refid,numrefmaps,mapid,nummaps);
		fflush(stdout);
	      } // omp critical
	    }

	    Cmap *nanomap = Gmap[mapid];
	    if(DEBUG>=2) assert(nanomap != 0);
	    if(DEBUG>=2 && !(nanomap->mapid == mapid)){
	      #pragma omp critical
	      {
		printf("mapid=%d:Gmap[mapid]->mapid=%d\n",mapid,Gmap[mapid]->mapid);
		fflush(stdout);
		assert(nanomap->mapid == mapid);
	      }
	    }
	    //	    int *MM = nanomap->numsite;
	    //	    FLOAT **XX = nanomap->site;
	    nanomap->incscale = 1.0;
	    nanomap->incwt = 0.0;

	    if(DEBUG) assert(KKmax != 0);

	    size_t origalignid;
	    if(DEBUG>=2) origalignid = alignid;

#if FLAT_MEMORY
	    refalignSD(rmap,nanomap,alignid,KKmax,refid,mapid,0,XPen,Qmin,Qmax,A,Fmem,Imem,StrideMax,JMIN,JMAX,Imin,Imax,tid,limN,NoSplit);
#else
	    refalignSD(rmap,nanomap,alignid,KKmax,refid,mapid,0,XPen,Qmin,Qmax,A,Amem,NoSplit);
#endif

	    if(DEBUG>=2) assert(align == alignment[origalignid]);
	    if(DEBUG && align){
	      for(int lc = 0; lc < colors; lc++){
		if(DEBUG && !(align[1 + lc].mapid1 == refid)){
                  #pragma omp critical
		  {
		    printf("refid=%d,mapid=%d:lc=%d,align[1+lc].mapid1=%d,align[1+lc].mapid2=%d\n",
			   refid,mapid,lc,align[1+lc].mapid1,align[1+lc].mapid2);
		    fflush(stdout);
		    assert(align[1 + lc].mapid1 == refid);
		  }
		}
		if(DEBUG) assert(align[1 + lc].mapid2 == mapid);
	      }
	    }

	    if(align && AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold)){/* good alignment */
	      Cmap *origmap = nanomap;
	      while(origmap->origmap)
		origmap = origmap->origmap;

	      if(origmap->nalign->mapid1 < 0 || (BESTREF_PV ? (align->logPV > origmap->nalign->logPV + BestRefMargin * (refid > origmap->nalign->mapid1 ? 1.0 : -1.0)) :
						 (align->score > origmap->nalign->score + BestRefMargin * (refid > origmap->nalign->mapid1 ? 1.0 : -1.0)))){
		/* save best alignment for next iteration */
		for(int c = 0; c < colors; c++)
		  origmap->nalign[c].update(align[c+1], NN[c]);
		/* replace score & logPV of color 0 by total score & logPV */
		origmap->nalign->score = align->score;
		origmap->nalign->logPV = align->logPV;
	      }

	      Tmapcnt++;

	      if(FirstAlignments && Tmapcnt >= 16){
                #pragma omp critical
		{
		  mapcnt += Tmapcnt;
		  Tmapcnt = 0;
		}
	      }
	    }

	    if(VERB>=2 && nanomap->id == MAP_TRACE && giter == RefRepeats-1){
              #pragma omp critical
	      {
		if(!align || align->numpairs <= 0)
		  printf("refid=%d,mapid=%d(id=%lld):no alignment found(thread=%d)\n",
			 refid,mapid,Gmap[mapid]->id,tid);
		else if (!AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold))/* no good alignment */
		  printf("refid=%d,mapid=%d(id=%lld):no good alignment:score=%0.4f,logPV=%0.3f,orientation=%d,scale=%d(%0.4f),pairs=%d(thread=%d)\n",
			 refid,mapid,Gmap[mapid]->id,align->score,align->logPV, align->orientation, align->scaleID, align->scaleID ? ScaleFactor[align->scaleID] : 1.0, align->numpairs,tid);
		else
		  printf("refid=%d,mapid=%d:score=%0.4f,logPV=%0.3f,orientation=%d,scale=%d(%0.4f)pairs=%d(thread=%d),mapcnt=%d\n",
			 refid,mapid,align->score,align->logPV, align->orientation, align->scaleID, align->scaleID ? ScaleFactor[align->scaleID] : 1.0, align->numpairs, tid, mapcnt+Tmapcnt);
		fflush(stdout);
	      } // pragma omp critical
	    }
	  } /* parallel for mapid = 0 .. nummaps -1 */

	  int start = orignummaps[refid];
	  int lim = totalmaps;
	  if(DEBUG && NoSplit >= 2) assert(start == lim);

	  if(VERB>=2){
	    #pragma omp critical
	    {
	      mapcnt += Tmapcnt;
	      Tmapcnt = 0;
	    }

	    if(numthreads > 1){
              #pragma omp barrier
	    }

	    #pragma omp master
	    {
	      printf("[refid=%d:total alignments=%d/%llu (alignids=%llu..%llu,nummaps=%d,totalmaps=%d)],start=%d,lim=%d: cum CPU time=%0.6f wall time=%0.6f secs\n",
		     refid,Tmapcnt + mapcnt,(unsigned long long)numaligns,(unsigned long long)numalign_start[refid],(unsigned long long)numalign_end[refid]-1,nummaps,totalmaps,start,lim,mtime(),wtime());
	      fflush(stdout);
	    }
	  }

	  while(start < lim){/* loop if more map splits are produced */
	    if(FirstAlignments && mapcnt + Tmapcnt > FirstAlignments)
	      break;

	    if(VERB>=2){
              #pragma omp master
	      {
		printf("refid=%d/%d(id=%lld),sites=%d+%d:processing map splits (%d..%d):alignids so far=%llu..%llu\n",
		       refid,numrefmaps,refmap[refid]->id,NN[0],NN[1],start,lim-1, (unsigned long long)numalign_start[refid],(unsigned long long)numalign_end[refid]-1);
		fflush(stdout);
	      }
	    }

	    /* allocate enough memory for alignment[] array before entering multithreaded loop : HERE HERE use defered allocation as in refalign */
	    size_t alignidstart = numalign_start[refid] + start - (orignummaps[refid] - nummaps);
	    #pragma omp master
	    {
	      size_t alignid = numalign_start[refid] + lim - (orignummaps[refid] - nummaps);
	      if(VERB>=2){
		printf("maxalignallocNULL(%llu):refid=%d,numalign_start[refid]=%llu,lim=%d,orignummaps[refid]=%d,nummaps=%d,start=%d,numaligns=%llu\n",
		       (unsigned long long)alignid,refid,(unsigned long long)numalign_start[refid],lim,orignummaps[refid],nummaps,start,(unsigned long long)numaligns);
		fflush(stdout);
	      }
	      maxalignallocNULL(alignid,alignment,numaligns,maxaligns);
	      numalign_end[refid] = numaligns = max(numaligns,alignid);
	      
#if 0 // old pre-allocation code
	      for(int mapid = start;mapid < lim;mapid++){
		size_t alignid = numalign_start[refid] + mapid - (mapid >= startmaps ? orignummaps[refid] - startmaps : 0);
		if(DEBUG) assert(alignid < numalign_end[refid]);
		Calign *align = alignment[alignid] = new Calign[colors+1];
		align->allfree();
		for(int c = 1; c <= colors; c++)
		  align[c].allfree();/* sets align[c].numpairs = 0 */
	      }
#endif
	    }

	    if(numthreads > 1){
              #pragma omp barrier
	    }

	    int block = min(3,max(1,(lim-start)/(64*numthreads)));

            #pragma omp for schedule(dynamic,block)
	    for(int mapid = start;mapid < lim;mapid++){

	      size_t alignid = numalign_start[refid] + mapid - (mapid >= startmaps ? orignummaps[refid] - startmaps : 0);
	      if(DEBUG>=2) assert(0 <= alignid && alignid < numaligns);

	      Calign * &align = alignment[alignid];// NOTE : pointer value may change in refalignSD !
	      Cmap *nanomap = Gmap[mapid];
	      if(DEBUG>=2) assert(nanomap != 0);
	      if(DEBUG>=2) assert(nanomap->mapid == mapid);
	      //	      int *MM = nanomap->numsite;
	      //	      FLOAT **XX = nanomap->site;

	      if(VERB && nummaps > 1000 && !(mapid%1000)){
                #pragma omp critical
		{
		  mapcnt += Tmapcnt;
		  Tmapcnt = 0;
		  printf("[refid=%d,mapid=%d/%d,alignments=%d]\n",refid,mapid,nummaps,mapcnt);
		  fflush(stdout);
		} // omp critical
	      }

	      nanomap->incscale = 0.0;
	      nanomap->incwt = 0.0;

	      if(DEBUG) assert(KKmax != 0);

	      //	      size_t origalignid = alignid;

#if FLAT_MEMORY
	      refalignSD(rmap,nanomap,alignid,KKmax,refid,mapid,0,XPen,Qmin,Qmax,A,Fmem,Imem,StrideMax,JMIN,JMAX,Imin,Imax,tid,limN,NoSplit);
#else
	      refalignSD(rmap,nanomap,alignid,KKmax,refid,mapid,0,XPen,Qmin,Qmax,A,Amem,NoSplit);
#endif

	      if(DEBUG && align)
		for(int c = 0; c < colors; c++){
		  assert(align[1+c].mapid1 == refid);
		  assert(align[1+c].mapid2 == mapid);
		}

	      if(align && AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold))/* good alignment */
		Tmapcnt++;
	      if(VERB>=2 && giter == RefRepeats-1){
                #pragma omp critical
		{
		  if(!align || align->numpairs <= 0)
		    printf("refid=%d,mapid=%d(id=%lld):no alignment found(thread=%d)\n",
			   refid,mapid,Gmap[mapid]->id,tid);
		  else if (!AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold)) /* no good alignment */
		    printf("refid=%d,mapid=%d(id=%lld):no good alignment:score=%0.4f,logPV=%0.3f,orientation=%d,scale=%d(%0.4f),pairs=%d(thread=%d)\n",
			   refid,mapid,Gmap[mapid]->id,align->score,align->logPV, align->orientation, align->scaleID, align->scaleID ? ScaleFactor[align->scaleID] : 1.0, align->numpairs,tid);
		  else
		    printf("refid=%d,mapid=%d:score=%0.4f,logPV=%0.3f,orientation=%d,scale=%d(%0.4f)pairs=%d(thread=%d)\n",
			   refid,mapid,align->score,align->logPV, align->orientation, align->scaleID, align->scaleID ? ScaleFactor[align->scaleID] : 1.0, align->numpairs, tid);
		  fflush(stdout);
		} // pragma omp critical
	      } // VERB>=2
	    } /* 2nd parallel for mapid = 0 .. nummaps -1 */

	    if(VERB>=2){
	      if(numthreads > 1){
                #pragma omp barrier
	      }
	      #pragma omp master
	      {
		printf("Split Alignments from %d map splits:\n",lim-start);
		for(size_t i = alignidstart; i < numaligns; i++){
		  Calign *align = alignment[i];
		  if(!align || align->numpairs <= 0)
		    continue;
		  if(AlignedThreshold(align, Gmap[align->mapid1]->site, ScoreThreshold, LogPvThreshold))
		    printf("%lu:refid=%d,mapid=%d(id=%lld),or=%d:score=%0.4f,logPV=%0.3f,pairs=%d\n",
			   i,align->mapid1,align->mapid2,Gmap[align->mapid2]->id,align->orientation,align->score,align->logPV,align->numpairs);
		}
		fflush(stdout);
	      }
	    }

	    start = lim;
	    lim = totalmaps;
	    if(DEBUG && NoSplit > 0) assert(start >= lim);
	  } /* while (start < lim) */

          #pragma omp critical
	  {
	    mapcnt += Tmapcnt;
	    Tmapcnt=0;
	  }

	  if(totalmaps > lim && !warning2){
            #pragma omp critical 
	    {
	      warning2 = 1;
	      printf("WARNING: parallel loop bound exceeded : may have missed some alignments:refid=%d,totalmaps=%d,lim=%d,start=%d,orignummaps=%d\n",refid,totalmaps,lim,start,orignummaps[refid]);
	      fflush(stdout);
	    }
	  }

	  /* free per-thread memory */
#if FLAT_MEMORY==0
	  for(int c = 0; c < colors; c++){
	    delete [] Amem[c]; Amem[c] = NULL;
	    delete [] AmemIK[c]; AmemIK[c] = NULL;
	    delete [] AmemIKJ[c]; AmemIKJ[c] = NULL;
	  }
#endif
	  XPen_free(XPen,Qmin,Qmax);
	} /* end of parallel section */

	if(VERB){
	  printf("[refid=%d:total alignments=%d/%lu (alignids=%lu..%lu)]\n",refid,mapcnt,numaligns,numalign_start[refid],numalign_end[refid]-1);
	  fflush(stdout);
	}

	if(MINMEM>=2){// does not help since score_init() memory prevents release of virtual memory
          delete rs_heap; 
	  rs_heap = new lightweight_heap(0,0);
	  if(VERB>=2){
	    printf("Reset lightweight heap\n");
	    fflush(stdout);
	  }
	}

	if(ALIGN_COMPRESS){/* compress alignment[numalign_start[refid]..numalign_end[refid]-1] to remove NULL pointers */
	  size_t j = numalign_start[refid];

	  if(REFDEBUG){ /* First remove below threshold alignments */
	    for(size_t k = j; k < numaligns; k++){
	      Calign *align = alignment[k];
	      if(!align)
		continue;
	      align++;
	      int origgoodalign = ((!REFDEBUG || (align[0].numpairs >= 2 && align[1].numpairs >=2)) && AlignedThreshold(&align[-1], rmap->site, origScoreThreshold, origLogPvThreshold)) ? 1 : 0;
	      if(DEBUG && !REFDEBUG && !origgoodalign){
		printf("refid=%d,mapid=%d,or=%d:score=%0.6f,logPV=%0.6f,numpairs=%d(%d,%d)\n",
		       align[-1].mapid1,align[-1].mapid2,align[-1].orientation,align[-1].score,align[-1].logPV,align[-1].numpairs,align[0].numpairs,align[1].numpairs);
		fflush(stdout);
		assert(origgoodalign);
	      }
	      if(!origgoodalign){
		delete [] alignment[k];
		alignment[k] = NULL;
	      }
	    }
	  }

	  size_t k = j;
	  for(; k < numaligns; k++){
	    if(!alignment[k]){
	      k++;
	      break;
	    }
	    j++;
	  }
	  for(; k < numaligns; k++){
	    if(!alignment[k])
	      continue;
	    if(DEBUG>=2) assert(j < k);
	    if(DEBUG>=2) assert(alignment[j]==0);
	    alignment[j] = alignment[k];
	    alignment[k] = 0;
	    j++;
	  }
	  if(VERB>=2 && j < numaligns){
	    printf("reduced numalign_end[%d] from %llu to %llu (numalign_start[%d]=%llu)\n",refid,(unsigned long long)numalign_end[refid],(unsigned long long)j,refid,(unsigned long long)numalign_start[refid]);

	    for(size_t t = 0; t < j; t++){
	      Calign *align = alignment[t];
	      if(!align || align->numpairs <= 0)
		continue;
	      if(AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold))
		printf("%lu:refid=%d,mapid=%d(id=%lld),or=%d:score=%0.4f,logPV=%0.3f,pairs=%d\n",
		       t,align->mapid1,align->mapid2,Gmap[align->mapid2]->id,align->orientation,align->score,align->logPV,align->numpairs);
	    }
	    fflush(stdout);
	  }
	  numalign_end[refid] = numaligns = j;
	}
      } /* refid 0..numrefmaps-1 */
    } /* General case (no hashtable) */

    if(MINMEM){
      delete rs_heap; 
      rs_heap = new lightweight_heap(0,0);
      if(VERB>=2){
	printf("Reset lightweight heap\n");
	fflush(stdout);
      }
    }

    //    if(VERB>=2)/* sort alignment by mapid2 (for debugging) */
    //      qsort(alignment, numaligns, sizeof(Calign*), (intcmp*)CalignMapidInc);

    if(MapRate > 0.0 || giter == RefRepeats-1 || BestRef){/* adjust LogPvThreshold to reduce effective mapping rate (if needed) to no more than MapRate */
      Calign **sortalign = new Calign*[numaligns];
      for(size_t i = 0; i < numaligns; i++)
	sortalign[i] = alignment[i];
      if(BestRef){/* replace logPV by 0 for alignments that will be filtered out by BestRef */
	if(giter==RefRepeats-1){/* Compute Raw Mapping rate */
	  qsort(sortalign, numaligns, sizeof(Calign*), (intcmp*)CalignLogPVDec);
	  size_t realaligns = numaligns;

	  // new code to provide Log10PV values that are integers 
	  double LogPv = floor(max(0.0,origLogPvThreshold) + 0.5);
	  double maxLogPv = MAX_PV + 0.01;
	  for(int k = 0; k < LogPv - 0.01 && k <= MAX_PV; k++)
	    RawMappingRatePV[k] = ((double)realaligns)/nummaps;
	  int t = realaligns-1;
	  for(; LogPv <= maxLogPv; LogPv += 1.0){
	    while(t >= 0 && sortalign[t]->logPV < LogPv)
	      t--;
	    if(t < 0)
	      break;
	    int k = (int)floor(LogPv+0.5);
	    RawMappingRatePV[k] = ((double)t)/nummaps;
	  }
	}

	for(size_t i = 0; i < numaligns; i++){
	  Calign *align = alignment[i];
	  if(DEBUG) assert(align != 0);
	  int refid = align->mapid1;
	  int mapid = align->mapid2;
	  if(DEBUG) assert(mapid >= 0);
	  Cmap *nanomap = Gmap[mapid];
	  Cmap *origmap = nanomap;
	  while(origmap->origmap)
	    origmap = origmap->origmap;
	  if(origmap->nalign->mapid1 != refid){	  
	    if(VERB>=3){
	      printf("alignment[%lu]=%p:refid=%d,mapid=%d,or=%d:numpairs=%d,nanomap->mapid=%d,origmap->mapid=%d,origmap->nalign->mapid1=%d,score=%0.6f->MINSCORE,logPV=%0.2f->-100.00\n",
		     i,align,align->mapid1,align->mapid2,align->orientation,align->numpairs,nanomap->mapid,origmap->mapid,origmap->nalign->mapid1,align->score,align->logPV);
	      fflush(stdout);
	    }
	    align->logPV = -100.0;
	    align->score = MINSCORE;
	  }
	}
      }
      if(MapRate > 0.0){
	if(BestRefPV){/* adjust logPV threshold */
	  qsort(sortalign, numaligns, sizeof(Calign*), (intcmp*)CalignLogPVDec);
	  size_t goalcnt = floor(nummaps * MapRate + 0.5);
	  if(goalcnt < numaligns && sortalign[goalcnt]->logPV > 0.0){
	    LogPvThreshold = sortalign[goalcnt]->logPV;
	    if(VERB){
	      printf("Changing LogPV threshold from %0.2f to %0.2f\n", origLogPvThreshold, LogPvThreshold);
	      fflush(stdout);
	    }
	  }
	} else {/* adjust score threshold */
	  qsort(sortalign, numaligns, sizeof(Calign*), (intcmp*)CalignScoreDec);
	  size_t goalcnt = floor(nummaps * MapRate + 0.5);
	  if(VERB>=2){
	    printf("MapRate=%0.4f,numaligns=%llu,nummaps=%d,goalcnt=%llu\n",MapRate,(unsigned long long)numaligns,nummaps,(unsigned long long)goalcnt);
	    fflush(stdout);
	  }
	  if(goalcnt < numaligns && sortalign[goalcnt]->score >= min(0.0,origScoreThreshold)){
	    ScoreThreshold = sortalign[goalcnt]->score;
	    if(VERB){
	      printf("Changing score threshold from %0.2f to %0.2f\n", origScoreThreshold, ScoreThreshold);
	      fflush(stdout);
	    }
	  }
	}
      }
      if(giter==RefRepeats-1){/* display effect of Pvalue cutoff on mapping rate */
	if(!(MapRate > 0.0 && BestRefPV))
	  qsort(sortalign, numaligns, sizeof(Calign*), (intcmp*)CalignLogPVDec);
	size_t realaligns = numaligns;
	if(BestRef){
	  for(realaligns=0;realaligns < numaligns; realaligns++)
	    if(sortalign[realaligns]->logPV <= 0.0)
	      break;
	  if(DEBUG && NoSplit==2 && !(realaligns <= (size_t)nummaps)){
	    printf("After BestRef filtering of alignments: realaligns=%llu,nummaps=%d\n",
		   (unsigned long long)realaligns,nummaps);
	    fflush(stdout);
	    assert(realaligns <= (size_t)nummaps);
	  }
	}

	// new code to provide Log10PV values that are integers 
	double LogPv = floor(max(0.0,origLogPvThreshold) + 0.5);
	double maxLogPv = MAX_PV + 0.01;
	for(int k = 0; k < LogPv - 0.01 && k <= MAX_PV; k++)
	  MappingRatePV[k] = ((double)realaligns)/nummaps;
	if(VERB && RefRepeats > 1)
	  printf("Log10PV threshold vs Mapping Rate (numaligns=%llu,LogPV range=%0.1f .. %0.1f,MapRate=%0.4f):\n",(unsigned long long)realaligns,LogPv,maxLogPv,MapRate);
	int t = realaligns-1;
	for(; LogPv <= maxLogPv; LogPv += 1.0){
	  while(t >= 0 && sortalign[t]->logPV < LogPv)
	    t--;
	  if(t < 0)
	    break;
	  int k = (int)floor(LogPv+0.5);
	  MappingRatePV[k] = ((double)t)/nummaps;
	  if(VERB && RefRepeats > 1){
	    if(BestRef)
	      printf("Log10PV=%0.2f MappingRate=%0.3f (coverage = %0.2f), MappingRate without BestRef=%0.3f\n", LogPv, MappingRatePV[k], MappingRatePV[k] * CoverageMult,RawMappingRatePV[k]);
	    else
	      printf("Log10PV=%0.2f MappingRate=%0.3f (coverage = %0.2f)\n", LogPv, MappingRatePV[k], MappingRatePV[k] * CoverageMult);
	  }
	}

	fflush(stdout);
      }
      delete [] sortalign;
    }

    if(DEBUG>=2) assert(0 <= giter && giter <= RefRepeats-1);

    if(VERB && BestRef){
      printf("Removing duplicate alignments of same map to different references:\n");
      fflush(stdout);
    }

    if(Refine>=2){
      if(DEBUG) assert(RefRepeats == 1);
      delete rs_heap; 
      rs_heap = NULL;
      if(VERB>=2){
	printf("Freed lightweight heap\n");
	fflush(stdout);
      }
    }

    mapcnt = 0;
    for(int refid = 0;refid < numrefmaps;refid++){
      Cmap *rmap = refmap[refid];
      FLOAT **YY = rmap->site;
      int *NN = rmap->numsite;

      if(NN[0]+NN[1] <= max(1,AlignedSiteThreshold) || YY[0][NN[0]] <= YY[0][1] + max(max(res[0],res[1])*PixelLen, AlignedLengthThreshold)){
	if(Refine){
	  if(VERB){
	    printf("Cannot Refine refmap %d: sites=%d,%d(AlignedSiteThreshold=%d),Y[1]=%0.3f,Y[N]=%0.3f (Empty refined map will be output)\n",
		   refid, NN[0], NN[1], AlignedSiteThreshold,YY[0][1],YY[0][NN[0]]);
	    fflush(stdout);
	  }

	  /* HERE : perhaps it would be better to copy the unrefined CMAP to preserve short regions between fragile sites that might merge due to SNR information */
	  int origCMapID = CMapID;
	  CMapID = refmap[refid]->id;
	  
	  char filename[PATH_MAX];
	  strcpy(filename,draft_prefix);
	  int i = strlen(filename);
	  sprintf(&filename[i],"_contig%lld.cmap",CMapID);
	  
	  if(!checkFile(filename)){
	    if(VERB){
	      printf("Generating empty %s (contig%lld draft consensus map in .cmap format)\n",filename, CMapID);
	      fflush(stdout);
	    }
	    
	    FILE *fp;
	    if((fp = fopen(filename,"w"))==NULL){
	      int eno = errno;
	      char *err = strerror(eno);
	      printf("failed to open file %s for writing contig%lld draft consensus(as .cmap):errno=%d:%s\n",filename,CMapID,eno,err);
	      exit(1);
	    }
	    FILEclose(fp);
	  }
	  CMapID = origCMapID;
	}

	continue;
      }

      if(VERB>=2){
	printf("refid=%d/%d:PoutlierEnd=%0.6f:OutlierEndBias=%0.6f,OutlierEndPenalty=%0.6f,OutlierPenalty=%0.6f,%0.6f:ChimScore=%0.6f\n",
	       refid,numrefmaps,PoutlierEnd,OutlierEndBias,OutlierEndPenalty,OutlierPenalty[0],OutlierPenalty[1],ChimScore);
	fflush(stdout);
      }

      //      double ExpOutlierPenalty = exp(OutlierPenalty);

      //      int origmapcnt = mapcnt;

      if(!NoStat){      /*  The following autonoise code is not needed during Refinement or if -nostat */

	if(REFDEBUG && !logLRarray){
	  logLRarray = new double[nummaps];
	  for(int c = 0; c < colors; c++){
	    logLRarrayC[c] = new double[nummaps];
	    logSDarray[c] = new double[nummaps];
	    Xlen[c] = new double[nummaps];
	  }
	  for(int i = 0; i < nummaps; i++){
	    logLRarray[i] = -1e5;
	    for(int c = 0; c < colors; c++){
	      logLRarrayC[c][i] = -1e5;
	      logSDarray[c][i] = 1e5;
	      Xlen[c][i] = -1.0;
	    }
	  }
	}

	Cmap *rmap = refmap[refid];// refmap[align->mapid1]

	for(size_t alignid = numalign_start[refid]; alignid < numalign_end[refid]; alignid++){
	  Calign *align = alignment[alignid];
	  if(DEBUG) assert(align != 0);
	  if(DEBUG && !(align->mapid1 == refid)){
	    printf("giter=%d,giter2=%d:refid=%d,numalign_start[refid]=%lu,numalign_end[refid]=%lu,alignid=%lu,align->mapid1=%d,align->mapid2=%d\n",
		   giter,giter2,refid,numalign_start[refid],numalign_end[refid],alignid,align->mapid1,align->mapid2);
	    fflush(stdout);
	    assert(align->mapid1 == refid);
	  }
	  if(DEBUG/* HERE >=2 */) assert(0 <= refid && refid < numrefmaps);
	  int mapid = align->mapid2;
	  if(DEBUG) assert(0 <= mapid);
	  if(DEBUG && !hash_filename) assert(mapid < startmaps || (NoSplit <= 1 && orignummaps[refid] <= mapid && (refid >= numrefmaps-1 ? mapid < totalmaps : mapid < orignummaps[refid+1])));
	  if(DEBUG && hash_filename && !(mapid < nummaps || (NoSplit <= 1 && startmaps <= mapid && mapid < totalmaps))){
	    printf("alignid=%lu:refid=%d,mapid1=%d,mapid2=%d(id=%lld),nummaps=%d,startmaps=%d,totalmaps=%d,NoSplit=%d\n",
		   alignid,refid,align->mapid1,align->mapid2,Gmap[align->mapid2]->id,nummaps,startmaps,totalmaps,NoSplit);
	    fflush(stdout);
	    assert(mapid < nummaps || (NoSplit <= 1 && startmaps <= mapid && mapid < totalmaps));
	  }

	  Cmap *nanomap = Gmap[mapid];
	  Cmap *origmap = nanomap;
	  while(origmap->origmap)
	    origmap = origmap->origmap;
	  if(VERB>=2 && align && AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold)){
	    printf("refid=%d/%d,mapid=%d/%d:align score=%0.6f,logPV=%0.2f,numpairs=%d,or=%d: best refid=%d:%s\n",
		   refid,numrefmaps,mapid,nummaps,align->score,align->logPV,align->numpairs,align->orientation, origmap->nalign->mapid1,
		   (BestRef && origmap->nalign->mapid1 != refid) ? "Skipping suboptimal refid" : "");
	    fflush(stdout);
	  }
	  if(BestRef && origmap->nalign->mapid1 != refid){
	    if(VERB>=2){
	      printf("BestRef Skipping: refid=%d/%d,mapid=%d/%d:align score=%0.6f,logPV=%0.2f,numpairs=%d,or=%d: best refid=%d\n",
		   refid,numrefmaps,mapid,nummaps,align->score,align->logPV,align->numpairs,align->orientation, origmap->nalign->mapid1);
	      fflush(stdout);
	    }
	    continue;/* skip alignment unless it is with the reference with the best alignment score with mapid */
	  }

	  if(DEBUG && align && !(align->mapid1 == refid)){
	    if(!hash_filename)
	      printf("refid=%d,mapid=%d,alignid=%llu:align->mapid1=%d,align->mapid2=%d,align->logPV=%0.2f,align->score=%0.3f,align->numpairs=%d:nummaps=%d,startmaps=%d,totalmaps=%d,orignummaps[refid]=%d,numalign_start[refid]=%llu,numalign_end[refid]=%llu\n",
		     refid,mapid,(unsigned long long)alignid,align->mapid1,align->mapid2,align->logPV,align->score,align->numpairs,nummaps,startmaps,totalmaps,orignummaps[refid],(unsigned long long)numalign_start[refid],(unsigned long long)numalign_end[refid]);
	    else
	      printf("refid=%d,mapid=%d,alignid=%llu:align->mapid1=%d,align->mapid2=%d,align->logPV=%0.2f,align->score=%0.3f,align->numpairs=%d:nummaps=%d,startmaps=%d,totalmaps=%d,numalign_start[refid]=%llu,numalign_end[refid]=%llu\n",
		     refid,mapid,(unsigned long long)alignid,align->mapid1,align->mapid2,align->logPV,align->score,align->numpairs,nummaps,startmaps,totalmaps,(unsigned long long)numalign_start[refid],(unsigned long long)numalign_end[refid]);
	    fflush(stdout);
	    assert(align->mapid1 == refid);
	  }
	  if(DEBUG) assert(align->mapid2 == mapid);

	  double scale = (align && align->numpairs > 0 && align->scaleID) ? ScaleFactor[align->scaleID] : 1.0;	    
	  double escale = ScaleDeltaBPP ? 1.0 : scale;
	  int *MM = nanomap->numsite;
	  FLOAT **XX = nanomap->site;

	  //	  if(AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold))/* good alignment */
	  //	    mapcnt++;

	  if(VERB>=2 && giter == RefRepeats-1){
	    if(align->numpairs <= 0)
	      printf("refid=%d,mapid=%d(id=%lld):no alignment found\n",
		     refid,mapid,Gmap[mapid]->id);
	    else if (!AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold)) /* no good alignment */
	      printf("refid=%d,mapid=%d(id=%lld):no good alignment:score=%0.4f,logPV=%0.3f,orientation=%d,scale=%d(%0.4f),pairs=%d\n",
		     refid,mapid,Gmap[mapid]->id,align->score,align->logPV, align->orientation, align->scaleID, align->scaleID ? ScaleFactor[align->scaleID] : 1.0, align->numpairs);
	    else
	      printf("refid=%d,mapid=%d:score=%0.4f,logPV=%0.3f,orientation=%d,scale=%d(%0.4f)pairs=%d,mapcnt=%d\n",
		     refid,mapid,align->score,align->logPV, align->orientation, align->scaleID, align->scaleID ? ScaleFactor[align->scaleID] : 1.0, align->numpairs, mapcnt);
	    fflush(stdout);
	  }

	  if(VERB >= EVERB && giter == RefRepeats-1 && align->numpairs <= 0){
	    printf("refid=%d,mapid=%d(id=%lld):no alignment found\n",
		   refid,mapid,Gmap[mapid]->id);
	    fflush(stdout);
	  }

	  align++;

	  if(1){
	    /* accumulate basic statistics for all alignments 
	       Also re-estimate error parameters from alignments scoring above the 3 thresholds :
	       NOTE : We only look at alignments scoring above the 4 thresholds (without adjustment with -MapRate) (origgoodalign==1).
	              However Autonoise only uses the subset of those alignments after adjustment with -MapRate (goodalign==1).
		      Unlike in refalign(), currently alignments below 4 thresholds (without adjustment with -MapRate) may be present and must be skipped.
	    */
	    int origgoodalign = ((!REFDEBUG || (align[0].numpairs >= 2 && align[1].numpairs >=2)) && AlignedThreshold(&align[-1], rmap->site, origScoreThreshold, origLogPvThreshold)) ? 1 : 0;
	    int goodalign = (origgoodalign && align[-1].score > ScoreThreshold && align[-1].logPV > LogPvThreshold) ? 1 : 0;

	    if(DEBUG && ALIGN_COMPRESS)
	      assert(origgoodalign);

	    if(!origgoodalign){
	      align[-1].score = MINSCORE;
	      continue;
	    }

	    tcnt++;

	    if(goodalign)
	      mapcnt++;

	    if(DEBUG && !BestRefPV && !((MapRate > 0.0 && !(align[-1].score > ScoreThreshold)) == (!goodalign))){
	      printf("refid=%d,mapid=%d,or=%d:MapRate=%0.6f,align[-1].score=%0.6f,ScoreThreshold=%0.6f,goodalign=%d\n",
		     refid,mapid,align[-1].orientation,MapRate,align[-1].score,ScoreThreshold,goodalign);
	      fflush(stdout);
	      assert((MapRate > 0.0 && !(align[-1].score > ScoreThreshold)) == (!goodalign));
	    }

	    //	    double trueoffset = nanomap->endloc - X[M+1];
	    double LR[3] = {0.0,0.0,0.0}, LRsd[2] = {0.0,0.0};

	    double y2sumC = 0.0;
	    double xysumC = 0.0;
	    double ysumC = 0.0;

	    /* first handle alignment for each color seperately */
	    for(int c = 0; c < colors; c++){
	      if(DEBUG) assert(align[c].orientation == align[-1].orientation);
	      if(DEBUG) assert(align[c].mapid1 == align[-1].mapid1);
	      if(DEBUG) assert(align[c].mapid2 == align[-1].mapid2);

	      int N = NN[c];
	      int M = MM[c];
	      FLOAT *Y = YY[c];
	      FLOAT *X = XX[c];

	      int U = align[c].numpairs;

	      int Lij = align[c].Lij1;
	      int Rij = align[c].Rij1;
	      int Lijx = align[c].Lij2;
	      int Rijx = align[c].Rij2;

	      Cprtab **PRtabY = PRtab[c][refid];

	      //	      double offset = 0.0;
	      int I=0,K=0,J=0,RI=0,RK=0,RJ=0;

	      if(U >= 1){
		I = align[c].sites1[0];
		K = align[c].sitesK1[0];
		J = align[c].sites2[0];
		RI = align[c].sites1[U-1];
		RK = align[c].sitesK1[U-1];
		RJ = align[c].sites2[U-1];

		//		offset = Yc(Y,I,K) - (align[c].orientation==0 ? X[J] : X[M+1]-X[M+1-J]);
	      }

	      scoresum += align[c].score;
	      logPVsum += align[c].logPV;
	      NumpairsSum += U;

	      
	      LR[c] = logLR(mapid,refid,c,X,Y,M,N,align,I,J,K,RJ,Lij,Rij,Lijx,Rijx,0,
			    /* (mapid==0 && c==0) ? 1 : 0 */ 0, SDDEBUG ? &LRsd[c] : NULL,NULL);
	      if(DEBUG) assert(isfinite(LR[c]));

	      /* NOTE : LR[0]+LR[1]+LR[2] will be compared against logLRarray[mapid] after exiting the c = 0 .. color-1 loop */

	      if(VERB && (VERB>=2 /* || (giter==0 && mapid==0 && c==0)*/)){
		if(logLRarray)
		  printf("iter=%d:refid=%d,mapid=%d(id=%lld),c=%d:or=%d,logLR()=%0.6f(previous alignment logLR=%0.6f),score=%0.6f,LRsum=%0.6f,N=%d,M=%d:best I=%d,K=%d,J=%d,Y[I]=%0.3f\n",
			 giter,refid,mapid,Gmap[mapid]->id,c,align[c].orientation,LR[c],logLRarrayC[c][mapid],align[c].score,LRsum,N,M, RI, RK, RJ, Y[RI]);
		else
		  printf("iter=%d:refid=%d,mapid=%d(id=%lld),c=%d:or=%d,logLR()=%0.6f,score=%0.6f,LRsum=%0.6f,N=%d,M=%d:best I=%d,K=%d,J=%d,Y[I]=%0.3f\n",
			 giter,refid,mapid,Gmap[mapid]->id,c,align[c].orientation,LR[c],align[c].score,LRsum,N,M,  RI, RK, RJ, Y[RI]);
		if(VERB>=2){
		  printf("  X[]=");
		  for(int j = 1; j <= M+1; j++)
		    printf("%0.3f,", X[j]);
		  printf("\n");
		  printf("N=%d:Y[1]=%0.3f,Y[N+1]=%0.3f\n",N,Y[1],Y[N+1]);
		}
		fflush(stdout);
	      }

	      if(DEBUG && !(isfinite(align[c].score) && isfinite(align[c].logPV))){
		printf("refid=%d,mapid=%d,c=%d,M=%d:invalid score or logPV:score=%0.4f,logPV=%0.3f,orientation=%d,pairs=%d,Lend=%d,Rend=%d\n",
		       refid,mapid,c,M,align[c].score,align[c].logPV, align[c].orientation,U,align[c].Lend,align[c].Rend);
		assert(isfinite(align[c].score));
		assert(isfinite(align[c].logPV));
	      }

	      if(VERB >= EVERB && giter == RefRepeats-1 && c==0 && !goodalign){/* no good alignment */
		printf("refid=%d,mapid=%d(id=%lld),M=%d+%d:no good alignment:score=%0.4f,logPV=%0.3f,orientation=%d,pairs=%d,Lend=%d,%d,Rend=%d,%d,len=%0.3f\n",
		       refid,mapid,Gmap[mapid]->id,MM[0],MM[1],align[-1].score,align[-1].logPV, align[-1].orientation,align[-1].numpairs,
		       align[0].Lend, align[1].Lend, align[0].Rend, align[1].Rend, X[M+1]*origPixelLen/PixelLen);
		fflush(stdout);
	      }

	      if(goodalign){/* good alignment */
		if(VERB >= EVERB && giter == RefRepeats-1 && c==0){
		  printf("refid=%d,mapid=%d,M=%d+%d:score=%0.4f,logPV=%0.3f,orientation=%d,pairs=%d,Lend=%d,%d,Rend=%d,%d (cnt=%d)\n",
			 refid,mapid,MM[0],MM[1],align[-1].score,align[-1].logPV, align[-1].orientation,align[-1].numpairs,
			 align[0].Lend,align[1].Lend,align[0].Rend,align[1].Rend, mapcnt);
		  fflush(stdout);
		}
		if(DEBUG/* >=2 */) assert(align->scaleID < max(1,NumScaleFactor));
		scaleIDcnt[align->scaleID]++;

		ATscoresum += align[c].score;
		ATlogPVsum += align[c].logPV;
		ATNumpairsSum += U;

		if(align[c].Lend <= -2 && align[c].sites2[0]-1 >= AlignedSiteThreshold)
		  chimcnt++;
		if(align[c].Rend <= -2 &&  M-align[c].sites2[U-1] >= AlignedSiteThreshold)
		  chimcnt++;
		if(align[c].Lend <= -2)
		  outlierEndcnt++;
		if(align[c].Rend <= -2)
		  outlierEndcnt++;
		if(c==0 && Gmap[align->mapid2]->origmap){
		  chimconf++;
		  Cmap *pmap = Gmap[align->mapid2]->origmap;
		  if(!NoSplit)
		    while(pmap->origmap)
		      pmap = pmap->origmap;
		  if(pmap->name && !(strstr(pmap->name,"R_") || strstr(pmap->name,"C_")))
		    chimconfFP++;
		}

		if(VERB >= EVERB && giter == RefRepeats-1){
		  printf("mapid=%d(id=%lld),c=%d:alignment range Y = %0.3f to %0.3f (X = %0.3f to %0.3f,len=%0.3f)\n", mapid, Gmap[mapid]->id,c,
			 Y[I-K]/* - (align[c].Lend <= -2 ? 0.0 : (align[c].orientation ? X[M+1]-X[M+1-J] : X[J]))*/,
			 Y[RI]/* + (align[c].Rend <= -2 ? 0.0 : (align[c].orientation ? X[M+1-RJ] : X[M+1]-X[RJ]))*/,
			 scale * (/*align[c].Lend > -2 ? 0.0 :*/ (align[c].orientation ? X[M+1]-X[M+1-J] : X[J])) * origPixelLen/PixelLen,
			 scale * (/*align[c].Rend > -2 ? X[M+1] :*/ (align[c].orientation ? X[M+1]-X[M+1-RJ] : X[RJ])) * origPixelLen/PixelLen, X[M+1]*origPixelLen/PixelLen);
		  fflush(stdout);
		}

		int Ysites = 0;


		/* left end interval */
		if(align[c].Lend > -2){
		  double x = escale * (align[c].orientation ? X[M+1]-X[M+1-J] : X[J]);
		  lenY[c] += x;
		  len[c] += x;
		  lenB[c] += max(0.0,x - resK2[c]);

		  for(int t = max(1,Lij); t < I-K; t++){
		    if(Y[I-K] - Y[t] >= FP_DIST){
		      fnI[c]++;
		      if(VERB>= EVERB && giter==RefRepeats-1){
			if(K>0)
			  printf("mapid=%d,c=%d:missing cut at Y[%d]=%0.3f(first aligned site X[%d]=%0.3f,Y[%d..%d]=%0.3f..%0.3f)\n",
				 mapid,c,t,Y[t],J,x * origPixelLen/PixelLen, I-K,I,Y[I-K],Y[I]);
			else
			  printf("mapid=%d,c=%d:missing cut at Y[%d]=%0.3f(first aligned site X[%d]=%0.3f,Y[%d]=%0.3f)\n",
				 mapid,c,t,Y[t],J,x * origPixelLen/PixelLen, I,Y[I]);
		      }
		    }
		  }
		  for(int t = max(1,Lijx); t < J; t++){
		    double Xt = escale * (align[c].orientation ? X[M+1]-X[M+1-t] : X[t]);
		    double XJ = escale * (align[c].orientation ? X[M+1]-X[M+1-J] : X[J]);
		    fp[c]++;
		    if(XJ - Xt >= FP_DIST){
		      fpI[c]++;
		      Xt *= origPixelLen/PixelLen;
		      XJ *= origPixelLen/PixelLen;
		      if(VERB>= EVERB && giter==RefRepeats-1){
			if(K>0)
			  printf("mapid=%d,c=%d:false cut at X[%d]=%0.3f(first aligned site X[%d]=%0.3f,Y[%d..%d]=%0.3f..%0.3f)\n",
				 mapid,c,t,Xt,J, XJ, I-K,I,Y[I-K],Y[I]);
			else
			  printf("mapid=%d,c=%d:false cut at X[%d]=%0.3f(first aligned site X[%d]=%0.3f,Y[%d]=%0.3f)\n",
				 mapid,c,t,Xt, J, XJ,I,Y[I]);
		      }
		    }
		  }
		  //		  fflush(stdout);

		  Ysites += I-K-max(1,Lij);
		  if(RESDATA)/* save left end aligned site information */
		    appendres(c,align,X,Y,M,N,I,K,J,0,0,0,0,0.0,escale,0);

		  double y = -1.0;
		  if(DEBUG>=2) assert(0 < Lij && Lij <= I-K);
		  if(DEBUG>=2) assert(Yc(Y,I,K) - Y[Lij-1] > 0.0);
		  if(ENDFIX>=2 && x <= Yc(Y,I,K) && Lij > 0 && ((y = Yc(Y,I,K)-Y[Lij-1]) < x || ENDFIX>=3))/* Sbnd(x,Y[I,K]-Y[Lij-1],n,m) */
		    appenderror(c,x,y,y,0.0,mapid, align[c].orientation ? M+1-J : 0, align[c].orientation  ? M+1 : J, 1, 0);
		}

		/* right end interval */
		if(align[c].Rend > -2){
		  double x = escale * (align[c].orientation ? X[M+1-RJ] : X[M+1]-X[RJ]);
		  lenY[c] += x;
		  len[c] += x;
		  lenB[c] += max(0.0,x - resK2[c]);
		  for(int t = RI+1; t <= min(N,Rij); t++){
		    if(Y[t]-Y[RI] >= FP_DIST){
		      fnI[c]++;
		      if(VERB>= EVERB && giter==RefRepeats-1){
			if(RK>0)
			  printf("mapid=%d,c=%d:missing cut at Y[%d]=%0.3f(last aligned site X[%d]=%0.3f,Y[%d..%d]=%0.3f..%0.3f)\n",
				 mapid,c,t,Y[t],RJ,x*origPixelLen/PixelLen,RI-RK,RI,Y[RI-RK],Y[RI]);
			else
			  printf("mapid=%d,c=%d:missing cut at Y[%d]=%0.3f(last aligned site X[%d]=%0.3f,Y[%d]=%0.3f)\n",
				 mapid,c,t,Y[t],RJ,x*origPixelLen/PixelLen,RI,Y[RI]);
		      }
		    }
		  }

		  for(int t = RJ+1; t <= min(M,Rijx); t++){
		    double Xt = escale * (align[c].orientation ? X[M+1]-X[M+1-t] : X[t]);
		    double XJ = escale * (align[c].orientation ? X[M+1]-X[M+1-RJ] : X[RJ]);
		    fp[c]++;
		    if(Xt - XJ >= FP_DIST){
		      fpI[c]++;
		      Xt *= origPixelLen/PixelLen;
		      XJ *= origPixelLen/PixelLen;
		      if(VERB>= EVERB && giter==RefRepeats-1){
			if(RK > 0)
			  printf("mapid=%d,c=%d:false cut at X[%d]=%0.3f(last aligned site X[%d]=%0.3f,Y[%d..%d]=%0.3f..%0.3f)\n",
				 mapid, c, t, Xt, RJ, XJ, RI-RK, RI, Y[RI-RK],Y[RI]);
			else
			  printf("mapid=%d,c=%d:false cut at X[%d]=%0.3f(last aligned site X[%d]=%0.3f,Y[%d]=%0.3f)\n",
				 mapid, c, t, Xt, RJ, XJ,  RI, Y[RI]);
		      }
		    }
		  }
		  //		  fflush(stdout);

		  Ysites += min(N,Rij)-RI;
		  double y = -1.0;
		  if(ENDFIX>=2 && x <= Y[N+1]-Yc(Y,RI,RK) && Rij <= N && ((y = Y[Rij+1]-Yc(Y,RI,RK)) < x || ENDFIX>=3))/* Sbnd(x,Y[Rij+1]-Yc(Y,RI,RK),M+1-RJ,Rij+1-RI) */
		    appenderror(c,x,y,y,0.0,mapid,align[c].orientation ? 0 : RJ, align[c].orientation ? M+1-RJ : M+1, 1, 0);
		}
		sitecnt[c] += Ysites + 1 + (RESSD ? K : ((SCORE_APPROX && K) ? ((SCORE_APPROX>=2 || MIS_VITERBI) ? K : 1) : 0));
		fn[c] += Ysites + ((RESSD==0 && SCORE_APPROX>=2 && K >= 2) ? K-1 : 0);
		if(FN_VITERBI && K){
		  double p = pTP[c];
		  double A = p*p/((2-p)*(2-p));
		  double B = 2*(1-p)*(1-p)/((2-p)*(2-p));
		  if(DEBUG) assert( A <= 1.0 && B <= 1.0);
		  sitecnt[c] += B-A;
		  fn[c] += B;
		}

		if(K > 0)/* add penalty for Sm()-log(pTP) to PrPen */
		  PrPen[c] += PRtabY[K][I].Sm/*Sm(J,I,K,Y,c)*/ - log_pTP[c];
		int H = I;
		int D = K;
		int G = J;
		double ysum = 0.0;/* sum of y for this mapid,color */
		double xysum = 0.0;/* sum of xy/var for this mapid,color */
		double x2sum = 0.0;/* sum of xx/var for this mapid,color */
		double y2sum = 0.0;/* sum of yy/var for this mapid,color */
		for(int T=1; T < U; T++, H=I, D=K, G=J){
		  I = align[c].sites1[T];
		  K = align[c].sitesK1[T];
		  J = align[c].sites2[T];
		  double y = Yc(Y,I,K) - Yc(Y,H,D);
		  double x = escale * ((align[c].orientation==0) ? X[J] - X[G] : X[M+1-G] - X[M+1-J]);
		  int m = J - G;
		  int n = I - K - H;

		  lenY[c] += y;

		  if(Poutlier > 0.0 && OUTLIER(maptype,x,y)){
		    outliertot++;
		    double Bias,Pen,PenSm;
		    SintDetail(x,y,m,n,J,I,K,D,PRtabY,Y,c,Bias,Pen,PenSm,0);
		    if(!FIX_OUTLIER_MIS)
		      Pen += PenSm;
		    double OutPen = OutlierPenalty[c];
		    if(outlierBC)
		      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];

		    double ExpOutPen = exp(OutPen);
		    double val = ExpOutPen/(ExpOutPen + exp(Pen));
		    if(DEBUG) assert(isfinite(val));
		    outliersum += val;
		    if(DEBUG) assert(isfinite(outliersum));

		    if(align[c].outscore[T] + OUTLIER_MARGIN < align[c].iscore[T]){
		      if(DEBUG && REFDEBUG) assert(Pen < OutlierPenalty[c]);
		      outliercnt++;
		      if(VERB >= EVERB && giter==RefRepeats-1){
			printf("mapid=%d,c=%d:outlier alignment between Y[%d..%d]=%0.3f..%0.3f and X[%d..%d]=%0.3f..%0.3f\n",
			       mapid,c,H,I,Yc(Y,H,D),Yc(Y,I,K),G,J,
			       escale*(align[c].orientation ? X[M+1]-X[M+1-G] : X[G])*origPixelLen/PixelLen,
			       escale*(align[c].orientation ? X[M+1]-X[M+1-J] : X[J])*origPixelLen/PixelLen);
			fflush(stdout);
		      }
		      if(FIX_OUTLIER_MIS){/* include Sm() term at right end, since it is not part of outlier */
			Cprtab *pr = &PRtabY[K][I];
			PrPen[c] += pr->Sm - log_pTP[c];
			sitecnt[c] += 1 + (RESSD ? K : ((SCORE_APPROX && K) ? ((SCORE_APPROX>=2 || MIS_VITERBI) ? K : 1) : 0));
			fn[c] += ((RESSD==0 && SCORE_APPROX>=2 && K >= 2) ? K-1 : 0);
			if(FN_VITERBI && K){
			  double p = pTP[c];
			  double A = p*p/((2-p)*(2-p));
			  double B = 2*(1-p)*(1-p)/((2-p)*(2-p));
			  if(DEBUG) assert( A <= 1.0 && B <= 1.0);
			  sitecnt[c] += B-A;
			  fn[c] += B;
			}
		      }
		      if(RESDATA)/* save res/resKB based information for current interval (including outliers) */
			appendres(c,align,X,Y,M,N,I,K,J,H,D,G,T,x,escale, 1);

		      continue;// don't include this interval in error parameter estimation
		    }
		  }
		  if(RESDATA)/* save res/resKB based information for current interval (including outliers) */
		    appendres(c,align,X,Y,M,N,I,K,J,H,D,G,T,x,escale, 0);

		  Cprtab *pr = &PRtabY[K][I];
		  PrPen[c] += pr->LogPr + pr->Sm - log_pTP[c];
		  len[c] += x;
		  lenB[c] += max(0.0,x - resK2[c]);

		  for(int t = H+1; t < I-K; t++){
		    if(Y[t]-Y[H] >= FP_DIST && Y[I-K]-Y[t] >= FP_DIST){
		      fnI[c]++;
		      if(VERB >= EVERB && giter==RefRepeats-1){
			if(D > 0 || K > 0)
			  printf("mapid=%d,c=%d:missing cut at Y[%d]=%0.3f (between aligned sites X[%d]=%0.3f,Y[%d..%d]=%0.3f..%0.3f and X[%d]=%0.3f,Y[%d..%d]=%0.3f..%0.3f)\n",
				 mapid,c,t,Y[t],G,escale*(align[c].orientation ? X[M+1]-X[M+1-G] : X[G])*origPixelLen/PixelLen,H-D,H,Y[H-D],Y[H], J,(align[c].orientation ? X[M+1]-X[M+1-J]: X[J])*origPixelLen/PixelLen,I-K,I,Y[I-K],Y[I]);
			else
			  printf("mapid=%d,c=%d:missing cut at Y[%d]=%0.3f (between aligned sites X[%d]=%0.3f,Y[%d]=%0.3f and X[%d]=%0.3f,Y[%d]=%0.3f)\n",
				 mapid,c,t,Y[t],G,escale*(align[c].orientation ? X[M+1]-X[M+1-G] : X[G])*origPixelLen/PixelLen,H,Y[H], J,(align[c].orientation ? X[M+1]-X[M+1-J]: X[J])*origPixelLen/PixelLen,I,Y[I]);
		      }
		    }
		  }
		  for(int t = G+1; t < J; t++){
		    double Xt = escale*(align[c].orientation ? X[M+1]-X[M+1-t] : X[t]);
		    double XG = escale*(align[c].orientation ? X[M+1]-X[M+1-G] : X[G]);
		    double XJ = escale*(align[c].orientation ? X[M+1]-X[M+1-J] : X[J]);
		    if(Xt-XG >= FP_DIST && XJ - Xt >= FP_DIST){
		      fpI[c]++;
		      Xt *= origPixelLen/PixelLen;
		      XG *= origPixelLen/PixelLen;
		      XJ *= origPixelLen/PixelLen;
		      if(VERB >= EVERB && giter==RefRepeats-1){
			if(D > 0 || K > 0)
			  printf("mapid=%d,c=%d:false cut at X[%d]=%0.3f (between aligned sites X[%d]=%0.3f,Y[%d..%d]=%0.3f..%0.3f and X[%d]=%0.3f,Y[%d..%d]=%0.3f..%0.3f)\n",
				 mapid, c, t, Xt, G, XG ,H-D,H,Y[H-D],Y[H], J, XJ, I-K,I,Y[I-K],Y[I]);
			else
			  printf("mapid=%d,c=%d:false cut at X[%d]=%0.3f (between aligned sites X[%d]=%0.3f,Y[%d]=%0.3f and X[%d]=%0.3f,Y[%d]=%0.3f)\n",
				 mapid,c, t, Xt, G, XG, H, Y[H], J, XJ, I, Y[I]);
		      }
		    }
		  }
		  fflush(stdout);

		  sitecnt[c] += n + (RESSD ? K : ((SCORE_APPROX && K) ? ((SCORE_APPROX>=2 || MIS_VITERBI) ? K : 1) : 0));
		  fn[c] += n-1 + ((RESSD==0 && SCORE_APPROX>=2 && K >= 2) ? K-1 : 0);
		  if(FN_VITERBI && K){
		    double p = pTP[c];
		    double A = p*p/((2-p)*(2-p));
		    double B = 2*(1-p)*(1-p)/((2-p)*(2-p));
		    if(DEBUG) assert( A <= 1.0 && B <= 1.0);
		    sitecnt[c] += B-A;
		    fn[c] += B;
		  }
		  fp[c] += m-1;

		  double var = SF[c]*SF[c] + fabs(SD[c])*SD[c]*y;
		  if(QUADRATIC_VARIANCE)
		    var += SR[c]*SR[c]*y*y;
		  double resvar = 0.0;
		  if(RES_VARIANCE){
		    double errL = Y[H] - Y[H-D];
		    double errR = Y[I] - Y[I-K];
		    resvar = errL*errL + errR*errR;
		  }
		  var += resvar * SE[c]*SE[c];

		  /* Update bppSD estimate */
		  double Ivar = 1.0/var;
		  ysum += y;
		  y2sum += y * y * Ivar;
		  xysum += y * x * Ivar;
		  x2sum += x * x * Ivar;
		  if(VERB>=3 && Gmap[mapid]->id == 2 && giter2==RefRepeats2-1 && giter==RefRepeats-1){
		    printf("    c=%d,T=%d/%d:x=%0.4f,y=%0.4f,sd=%0.6f,Ivar=%0.8f:y2sum=%0.8f,xysum=%0.8f,x2sum=%0.8f\n",
			   c,T,U,x,y,sqrt(var),Ivar,y2sum,xysum,x2sum);
		    fflush(stdout);
		  }

		  segcnt[c]++;

		  /* append error information of current interval to errors[c][0..numerrors[c]-1] */
		  appenderror(c,x,y,y,resvar,mapid,align[c].orientation ? M+1-J : G, align[c].orientation ? M+1-G : J, 0,0);
		}// for(T=1; T < U; T++)

		if(DEBUG) assert(H==RI && G==RJ && D==RK);
		if(DEBUG>=2) assert(align[-1].numpairs >= AlignedSiteThreshold);

		y2sumC += y2sum;
		xysumC += xysum;
		ysumC += ysum;

		if(VERB>=2 && Gmap[mapid]->id == 2 && giter2==RefRepeats2-1 && giter==RefRepeats-1){
		  printf("mapid=%d,c=%d:score=%0.2f,logPV=%0.3f,numpairs=%d: local scaling = %0.6f(sum=%0.6f), wt=%0.6f(sum=%0.6f)\n",
			 mapid,c,align[c].score,align[c].logPV,align[c].numpairs, y2sum/xysum, nanomap->incscale, ysum, nanomap->incwt);
		  fflush(stdout);
		}
	      } // if(goodalign)
	    } /* c = 0 .. colors-1 */

	    if(/* align[-1].numpairs >= AlignedSiteThreshold && */ xysumC > 0.0){/* update local Molecule scaling factor */
	      nanomap->incscale += y2sumC/xysumC;
	      nanomap->incwt += ysumC;
	    }

	    /* next handle color misalignment penalty term */
	    if(colors == 2){
	      int U = align[0].numpairs;
	      int I = align[0].sites1[U-1];
	      int K = align[0].sitesK1[U-1];
	      int J = align[0].sites2[U-1];
	      U = align[1].numpairs;
	      int P = align[1].sites1[U-1];
	      int T = align[1].sitesK1[U-1];
	      int Q = align[1].sites2[U-1];
	      double y0 = Yc(YY[0],I,K);
	      double y1 = Yc(YY[1],P,T);
	      double x0 = escale * (align[-1].orientation ? (XX[0][MM[0]+1] - XX[0][MM[0]+1-J]) : XX[0][J]);
	      double x1 = escale * (align[-1].orientation ? (XX[1][MM[1]+1] - XX[1][MM[1]+1-Q]) : XX[1][Q]);
	      if(DEBUG>=2) assert(fabs(XX[0][MM[0]+1]-XX[1][MM[1]+1]) < 1e-8);
	      if(DEBUG>=2) assert(fabs(YY[0][NN[0]+1]-YY[1][NN[1]+1]) < 1e-8);
	      double y = y1 - y0;
	      double x = x1 - x0;
	      double len = (COLOR_SHIFT==2) ? max(y0,y1) : max(x0,x1);

	      if(goodalign){
		if(COLOR_SHIFT==0)
		  appenderror(colors,x,y,fabs(y),0.0,mapid,J,Q, (align[-1].orientation) ? -1 : 1,0);
		else
		  appenderror(colors,x,y,len,0.0,mapid,J,Q, (align[-1].orientation) ? -1 : 1,0);
	      }
	      LR[2] = SMA(x,y,len,align[-1].orientation);

	      if((DEBUG && !isfinite(LR[2])) || (VERB && (VERB>=2 /* || (giter<=1 && mapid==16)*/))){
		printf("refid=%d,mapid=%d,or=%d,c=2:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,y=%0.4f,len=%0.4f,SMA=%0.6f:MA_mean=%0.4f,SF[2]=%0.8e,FM2var=%0.8e,SD[2]=%0.8e,SM2var=%0.8e,RefPen12=%0.8e,Xtheta12=%0.4f\n",
		       refid,mapid,align[-1].orientation,I,K,J,P,T,Q,x0,x1,y0,y1,y,len,LR[2],MA_mean,SF[2],FM2var,SD[2],SM2var,RefPen12,Xtheta12);
		fflush(stdout);
		assert(isfinite(LR[2]));
	      }

	      double LRtot = LR[0]+LR[1]+LR[2];
	      if(DEBUG) assert(isfinite(LRtot));
	      if(REFDEBUG && !HashGen && (numrefmaps <= 1 || (BestRef && !BestRefPV)) && logLRarray && logLRarray[mapid] > LARGE_NEGATIVE && (!(LRtot > logLRarray[mapid] - (USE_RFLOAT ? 1e-4 : 1e-6)))){
		printf("WARNING: refid=%d(id=%lld),mapid=%d(id=%lld),orientation=%d:logLR[0,1,2]=%0.6f,%0.6f,%0.6f(LRtot=%0.6f, previous alignment=%0.6f(logLRC[0,1]=%0.6f,%0.6f)\n",
		       refid, refmap[refid]->id,mapid,Gmap[mapid]->id,align[-1].orientation,LR[0],LR[1],LR[2],LRtot,logLRarray[mapid], logLRarrayC[0][mapid], logLRarrayC[1][mapid]);
		for(int c = 0; c < colors; c++){
		  int N = NN[c];
		  int M = MM[c];
		  FLOAT *Y = YY[c];
		  FLOAT *X = XX[c];
		  int U = align[c].numpairs;
		  int Lij = align[c].Lij1;
		  int Rij = align[c].Rij1;
		  int Lijx = align[c].Lij2;
		  int Rijx = align[c].Rij2;
		  int I = align[c].sites1[0];
		  int K = align[c].sitesK1[0];
		  int J = align[c].sites2[0];
		  int RJ = align[c].sites2[U-1];

		  (void)logLR(mapid,refid,c,X,Y,M,N,align,I,J,K,RJ,Lij,Rij,Lijx,Rijx,0,1,NULL,NULL);
		}
		printf("MA:orientation=%d,x=%0.3f,y=%0.3f,len=%0.3f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,SMA=%0.6f\n",align[-1].orientation,x,y,len,I,K,J,P,T,Q,LR[2]);
		printf("SMA(x,y,len,orientation)=%0.6f:MA_mean=%0.4f,SF[2]=%0.8e,FM2var=%0.8e,SD[2]=%0.8e,SM2var=%0.8e\n",
		       SMA(x,y,len,align[-1].orientation),MA_mean,SF[2],FM2var,SD[2],SM2var);
		fflush(stdout);

		if(!HashGen) assert(LRtot > logLRarray[mapid] - (USE_RFLOAT ? 1e-4 : 1e-6));
	      }
	      if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		LRsum += LRtot;
		for(int lc = 0; lc <= colors; lc++)
		  LRsumC[lc] += LR[lc];
	      }
	      if(REFDEBUG){
		logLRarray[mapid] = LRtot;
		for(int c = 0; c < colors; c++){
		  logLRarrayC[c][mapid] = LR[c];
		  logSDarray[c][mapid] = LRsd[c];
		}
	      }
	      if(DEBUG) assert(isfinite(LRsum));
	      if(goodalign){/* good alignment */
		ATlogLRsum += LRtot;
		logLRsdsum += LRsd[0]+LRsd[1];
	      }
	    }

	    if(VERB>=2){
	      if(align[-1].numpairs <= 1) {
		printf("refid=%d,mapid=%d:no valid alignment:score=%0.4f,orientation=%d,pairs=%d,Lend=%d,%d,Rend=%d,%d\n",
		       refid,mapid,align[-1].score,align[-1].orientation,align[-1].numpairs,align[0].Lend,align[1].Lend,align[0].Rend,align[1].Rend);
		exit(1);
	      } else {
		int I = align[0].sites1[0];
		int K = align[0].sitesK1[0];
		int J = align[0].sites2[0];
		int M = MM[0];
		FLOAT *X = XX[0];
		FLOAT *Y = YY[0];
		double offset = Yc(Y,I,K) - scale * (align->orientation==0 ? X[J] : X[M+1]-X[M+1-J]);
		double trueoffset = nanomap->endloc - scale * X[M+1];
		if(fabs(offset- trueoffset) > 2*Ylambda12)
		  errcnt++;
		totcnt++;
		printf("refid=%d,mapid=%d,N=%d,M=%d:best alignment:score=%0.4f,logPV=%0.4f,orientation=%d,pairs=%d,loc=%0.3f,offset=%0.3f%c(%lu/%lu):Lend=%d,Rend=%d,chimcnt=%d,confirmed=%d(FP=%d),pairs=%d,id=%lld,origid=%lld\n",
		       refid,mapid,NN[0],MM[0],align[-1].score,align[-1].logPV,align[-1].orientation,align[-1].numpairs,offset,trueoffset,
		       fabs(offset - trueoffset) > 2*Ylambda12 ? '!' : ' ',errcnt,totcnt,
		       align[0].Lend,align[0].Rend,chimcnt,chimconf,chimconfFP,chimpaircnt,Gmap[mapid]->id,Gmap[mapid]->origmap ? Gmap[mapid]->origmap->id : -1);
	      }
	      fflush(stdout);
	    }
	  }/* if (align[-1].numpairs >= 1) */
	} /* loop over alignid/mapid */
      } /* if ( !NoStat) */ else {/* Update mapcnt : cumulative number of maps with good alignment to refid */

	for(size_t alignid = numalign_start[refid]; alignid < numalign_end[refid]; alignid++){
	  Calign *align = alignment[alignid];
	  if(DEBUG) assert(align != 0);
	  int mapid = align->mapid2;
	  if(DEBUG) assert(mapid >= 0);
	  if(DEBUG && !hash_filename) assert(mapid < startmaps || (NoSplit <= 1 && orignummaps[refid] <= mapid && (refid >= numrefmaps-1 ? mapid < totalmaps : mapid < orignummaps[refid+1])));
	  if(DEBUG && hash_filename) assert(mapid < nummaps || (NoSplit <= 1 && startmaps <= mapid && mapid < totalmaps));
	  Cmap *nanomap = Gmap[mapid];
	  Cmap *origmap = nanomap;
	  while(origmap->origmap)
	    origmap = origmap->origmap;
	  if(BestRef && origmap->nalign->mapid1 != refid)
	    continue;/* skip alignment unless it is with the reference with the best alignment score with mapid */
	  if(DEBUG && align && !(align->mapid1 == refid)){
	    if(!hash_filename)
	      printf("refid=%d,mapid=%d,alignid=%llu:align->mapid1=%d,align->mapid2=%d,align->logPV=%0.2f,align->score=%0.3f,align->numpairs=%d:nummaps=%d,startmaps=%d,totalmaps=%d,orignummaps[refid]=%d,numalign_start[refid]=%llu,numalign_end[refid]=%llu\n",
		     refid,mapid,(unsigned long long)alignid,align->mapid1,align->mapid2,align->logPV,align->score,align->numpairs,nummaps,startmaps,totalmaps,orignummaps[refid],
		     (unsigned long long)numalign_start[refid],(unsigned long long)numalign_end[refid]);
	    else
	      printf("refid=%d,mapid=%d,alignid=%llu:align->mapid1=%d,align->mapid2=%d,align->logPV=%0.2f,align->score=%0.3f,align->numpairs=%d:nummaps=%d,startmaps=%d,totalmaps=%d,numalign_start[refid]=%llu,numalign_end[refid]=%llu\n",
		     refid,mapid,(unsigned long long)alignid,align->mapid1,align->mapid2,align->logPV,align->score,align->numpairs,nummaps,startmaps,totalmaps,
		     (unsigned long long)numalign_start[refid],(unsigned long long)numalign_end[refid]);
	    fflush(stdout);
	    assert(align->mapid1 == refid);
	  }
	  if(DEBUG && align) assert(align->mapid2 == mapid);

	  Cmap *rmap = refmap[refid];
	  if(align && AlignedThreshold(align, rmap->site, ScoreThreshold, LogPvThreshold))/* good alignment */
	    mapcnt++;
	}
      } // if(NoStat)

      if(VERB>=2 && BestRef){
	printf("[refid=%d:alignments=%d/%lu (alignids=%lu..%lu)]\n",refid,mapcnt,numaligns,numalign_start[refid],numalign_end[refid]-1);
	fflush(stdout);
      }

      if(Refine && giter == RefRepeats-1){/* create a refined reference map and output it */
	/* HERE HERE : Copy code from refalign.cpp and update to 2 colors */
	printf("-refine not implemented with %d colors\n",colors);
	exit(1);
      }/* if(Refine && ... ) */
    } /* for ( refid = 0; refid < numrefmaps; refid++) */
    if(DEBUG) assert(numalign_end[numrefmaps-1] >= numaligns);

    /* save nanomap->nalign to nanomap->align for use during next giter with FAST (or during output with BestRef) */
    for(int mapid = 0;mapid < startmaps;mapid++){
      Cmap *nanomap = Gmap[mapid];      
      if(nanomap->nalign->mapid1>=0){
	for(int c = 0; c < colors; c++) {
	  nanomap->align[c].update(nanomap->nalign[c]);
	  nanomap->nalign[c].reset();
	}
      }
    }

    if(DEBUG>=2) assert(0 <= giter && giter <= RefRepeats-1);

    if(giter == RefRepeats-1){
      if(VERB/* HERE >=2 */){
	printf("refalign:nummaps=%d,startmaps=%d,totalmaps=%d,refmaps=%d:end time=%0.6f(Elapsed=%0.6f)\n",nummaps,startmaps,totalmaps,numrefmaps,mtime(),wtime());
	fflush(stdout);
      }
      int orignummaps = nummaps;
      nummaps = totalmaps;/* most output routines expect numaps to be the total size of Gmap[] array */

      if(svcheck){
	printf("-svcheck not yet implemented for 2 colors\n");
	exit(1);
	/* HERE HERE : copy code from refalign.cpp and update to 2 colors */
      }

      if(DEBUG>=2) assert(giter == RefRepeats-1);

      int split_maps = (MappedUnsplit >= 0) ? !MappedUnsplit : (SplitRef && CMapID > 0) ? 1 : (svcheck && numrefmaps > 0) ? 1 : 0;
      if(VERB>=2){
	printf("split_maps=%d\n",split_maps);
	fflush(stdout);
      }
      qsort(refmap,numrefmaps,sizeof(Cmap *), (intcmp*)CmapIdInc);
      output_csites(0, numrefmaps-1, output_prefix, split_maps, alignment, numalign_start, numalign_end, false); /* output condensed reference map as *_q.cmap */

      if(BestAlignments > 0){
	printf("WARNING:-bestalignments not implemtned for 2 colors\n");
	fflush(stdout);
      }

      if(Ch > 0.0 && ChFdr > 0.0){
	printf("WARNING: .chim outout file not yet implemented for 2 colors\n");
	fflush(stdout);
	//	output_chimmaps(output_prefix);/* output .chim file */
      }

      if(DEBUG>=2) assert(giter == RefRepeats-1);

      if(MappedPrefix){/* output maps that scored above threshold to <MappedPrefix>.cmap */
	if(VERB){
	  printf("Starting output of mapped bnx files: CPU time=%0.6f, wall time=%0.6f\n", mtime(),wtime());
	  fflush(stdout);
	}

	if(groupfile) {
	  int numgroups = input_groups(groupfile,refmap,numrefmaps);/* refmap[0..numrefmaps-1]->paired are set to group ids ranging from 1 .. numgroups */
	  int numLL = (numgroups + 63)/64;/* number of 64bit words required to represent a bitmap of size = numgroups (bits) : NOTE : bits runs 0..numgroups-1*/
	  unsigned long long *MapSave = new unsigned long long[nummaps * numLL];/* bitmap used to keep track if which maps have been saved to each group output (to avoid duplicates) */
	  memset(MapSave,0,nummaps*numLL*sizeof(unsigned long long));

	  FILE **groupFP  =  new FILE*[numgroups+1];
	  char **groupbuffer = new char*[numgroups+1];
	  omp_lock_t *lock = new omp_lock_t[numgroups+1];/* omp file locks */
	  for(int i = 1; i <= numgroups; i++){
	    groupFP[i] = NULL;
	    groupbuffer[i] = NULL;
	    #ifdef _OPENMP
	    omp_init_lock(&lock[i]);
	    #endif
	  }

	  int nthreads = max(1,min(numthreads,numrefmaps));

	  #pragma omp parallel num_threads(nthreads) if(nthreads > 1)
	  {
	    int tid = -1;
            #ifdef _OPENMP
	    tid = omp_get_thread_num ();
            #endif

	    /* assemble all good maps for one refid into array goodmap[0..goodmapcnt-1] */
	    char MappedIdPrefix[PATH_MAX];
	    Cmap **goodmaps = new Cmap *[mapcnt];
	    
	    #pragma omp for schedule(dynamic,1)
	    for(int refid = 0; refid < numrefmaps; refid++){
	      Cmap *rmap = refmap[refid];
	      int group = rmap->paired;
	      if(VERB>=2){
		printf("refid=%d:group=%d (-S %0.3f -T %0.6e -A %d -L %0.2f, -E %d)\n",refid,group,ScoreThreshold,LogPvThreshold,AlignedSiteThreshold, AlignedLengthThreshold, AlignedEndOutlierThreshold);
		fflush(stdout);
	      }
	      if(DEBUG) assert(0 < group && group <= numgroups);
	      int groupWord = (group-1)/64;
	      int groupBit = (group-1)%64;
	      unsigned long long groupMask = 1LL << groupBit;
	      size_t align_start = numalign_start[refid], align_end = numalign_end[refid];
	      int numgoodmaps = 0;
	      for(size_t i = align_start; i < align_end; i++){
		Calign *p = alignment[i];
		if(!p || p->numpairs <= 0)
		  continue;
		if(DEBUG) assert(p->mapid1 == refid);
		Cmap *rmap = refmap[refid];
		if(!AlignedThreshold(p, rmap->site, ScoreThreshold, LogPvThreshold))
		  continue;
		Cmap *Xmap = Gmap[p->mapid2];
		if(Xmap->origmap)
		  continue;/* skip map fragments : the best alignment of a map fragment should always be worse than that of the full molecules */
		if(BestRef){
		  Cmap *origmap = Xmap;
		  while(origmap->origmap)
		    origmap = origmap->origmap;
		  if(origmap->align->mapid1 != p->mapid1)
		    continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
		}
		//		if(BestRef && Xmap->align->mapid1 != p->mapid1)
		//		  continue;/* skip alignment unless it is with the reference with the best alignment score with this map Xmap */

		if(CenterMapped > 0.0){/* skip alignment if it overlaps anywhere within CenterMapped kb of reference contig ends */
		  int c = 0;
		  for(;c < colors; c++){
		    double scale = p->scaleID ? ScaleFactor[p->scaleID] : 1.0;	    
		    FLOAT *X = Xmap->site[c];
		    FLOAT *Y = rmap->site[c];
		    int M = Xmap->numsite[c];
		    int N = rmap->numsite[c];
		    int I = p[c+1].sites1[0];
		    int K = p[c+1].sitesK1[0];
		    int J = p[c+1].sites2[0];
		    int RI = p->sites1[p[c+1].numpairs-1];
		    int RK = p->sitesK1[p[c+1].numpairs-1];
		    int RJ = p->sites2[p[c+1].numpairs-1];
		    double leftend = (p[c+1].Lend > -2) ? Yc(Y,I,K) - scale*(p[c+1].orientation ? X[M+1]-X[M+1-J] : X[J]) : Yc(Y,I,K);
		    double rightend = (p[c+1].Rend > -2) ?  Yc(Y,RI,RK) + scale*(p[c+1].orientation ? X[M+1-RJ] : X[M+1] - X[RJ]) : Yc(Y,RI,RK);
		    if((leftend < CenterMapped) || (rightend > Y[N+1] - CenterMapped))
		      break;
		  }
		  if(c < colors)
		    continue;
		}

		if(DEBUG>=2) assert(Xmap->mapid == p->mapid2);
		goodmaps[numgoodmaps++] = Xmap;
	      }

	      if((DEBUG && !(numgoodmaps <= mapcnt)) || VERB>=2){
		printf("refid=%d:numgoodmaps=%d,mapcnt=%d,nummaps=%d,startmaps=%d,totalmaps=%d\n",refid,numgoodmaps,mapcnt,orignummaps,startmaps,totalmaps);
		fflush(stdout);
		assert(numgoodmaps <= mapcnt);
	      }

	      /* check for duplicate maps (requires read/write of shared bitmap MapSave[]) */
	      #pragma omp critical(MapSavelock)
	      {
		for(int i = 0; i < numgoodmaps; i++){
		  int Mindex = goodmaps[i]->mapid * numLL + groupWord;
		  if(VERB>=2){
		    printf("refid=%d:group=%d,goodmaps[%d]:mapid=%d(%lld),Mindex=%d,groupMask=%llu,MapSave[Mindex]=%llu\n",refid,group,i,goodmaps[i]->mapid,goodmaps[i]->id,Mindex,groupMask,MapSave[Mindex]);
		    fflush(stdout);
		  }
		  if(MapSave[Mindex] & groupMask){
		    if(VERB>=2){
		      printf("refid=%d:deleted mapid=%d(%lld) from group=%d\n",refid,goodmaps[i]->mapid,goodmaps[i]->id,group);
		      fflush(stdout);
		    }
		    goodmaps[i] = 0;
		  } else
		    MapSave[Mindex] |= groupMask;
		}
	      }

	      int j = 0;
	      for(int i = 0; i < numgoodmaps; i++){
		if(!goodmaps[i])
		  continue;
		goodmaps[j++] = goodmaps[i];
	      }
	      if(j < numgoodmaps){
		if(VERB>=2){
		  printf("refid=%d:numgoodmaps=%d->%d due to duplicate maps in same group (groupWord=%d,groupBit=%d,groupMask=%llu)\n",refid,numgoodmaps,j,groupWord,groupBit,groupMask);
		  fflush(stdout);
		}
		numgoodmaps = j;
	      }

	      sprintf(MappedIdPrefix,"%s_group%d",MappedPrefix, group);

	      #ifdef _OPENMP
	      omp_set_lock(&lock[group]);/* lock file */
	      #endif

	      if(maptype==0 || MapType >= 0 || CmapMergeBnx){
		output_bnx(MappedIdPrefix,goodmaps,0,numgoodmaps-1, 1, &groupFP[group], &groupbuffer[group], tid);
	      } else {
		printf("CMAP output of grouped _mapped files not implemented: Use -maptype or -bnx to output BNX file\n");
		exit(1);
		//	output_cmap(MappedIdPrefix,goodmaps,0,numgoodmaps-1, &groupFP[group]);
	      }

	      #ifdef _OPENMP
	      omp_unset_lock(&lock[group]);/* unlock file */
	      #endif

	    }// refid = 0 .. numrefmaps-1
	  } // parallel

          delete [] MapSave;
	  for(int i = 1; i <= numgroups; i++){
	    if(groupFP[i])
	      FILEclose(groupFP[i]);
	    if(groupbuffer[i])
	      free(groupbuffer[i]);
	    #ifdef _OPENMP
	    omp_destroy_lock(&lock[i]);
	    #endif
	  }
	  delete [] groupFP;
	  delete [] groupbuffer;

	} else {// No groupfile

	  /* assemble all good maps for one refid into array goodmaps[0..numgoodmaps-1] */
	  char MappedIdPrefix[PATH_MAX];
	  Cmap **goodmaps = new Cmap *[mapcnt];

	  int refmax = split_maps ? numrefmaps : 1;

	  for(int refid = 0; refid < refmax; refid++){
	    size_t align_start = 0, align_end = numaligns;
	    if(!split_maps)
	      for(int i = 0; i < nummaps; i++)
		Gmap[i]->paired = 0;
	    else {
	      align_start = numalign_start[refid];
	      align_end = numalign_end[refid];
	    }

	    int numgoodmaps = 0;

	    for(size_t i = align_start; i < align_end; i++){
	      Calign *p = alignment[i];
	      if(!p || p->numpairs <= 0)
		continue;
	      if(DEBUG && split_maps) assert(p->mapid1 == refid);
	      Cmap *rmap = refmap[p->mapid1];
	      if(!AlignedThreshold(p, rmap->site, ScoreThreshold, LogPvThreshold))
		continue;
	      Cmap *Xmap = Gmap[p->mapid2];
	      if(Xmap->origmap)
		continue;/* skip map fragments : the best alignment of a map fragment should always be worse than that of the full molecules */
	      if(BestRef){
		Cmap *origmap = Xmap;
	        while(origmap->origmap)
		  origmap = origmap->origmap;
		if(origmap->align->mapid1 != p->mapid1)
		  continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	      }
	      //	      if(BestRef && Xmap->align->mapid1 != p->mapid1)
	      //		continue;/* skip alignment unless it is with the reference with the best alignment score with this map Xmap */

	      if(CenterMapped > 0.0){/* skip alignment if it overlaps anywhere within CenterMapped kb of reference contig ends */
	        double scale = p->scaleID ? ScaleFactor[p->scaleID] : 1.0;	    
		int c = 0;
		for(;c < colors; c++){
		  FLOAT *X = Xmap->site[c];
		  FLOAT *Y = rmap->site[c];
		  int M = Xmap->numsite[c];
		  int N = rmap->numsite[c];
		  int I = p[c+1].sites1[0];
		  int K = p[c+1].sitesK1[0];
		  int J = p[c+1].sites2[0];
		  int RI = p[c+1].sites1[p[c+1].numpairs-1];
		  int RK = p[c+1].sitesK1[p[c+1].numpairs-1];
		  int RJ = p[c+1].sites2[p->numpairs-1];
		  double leftend = (p[c+1].Lend > -2) ? Yc(Y,I,K) - scale*(p[c+1].orientation ? X[M+1]-X[M+1-J] : X[J]) : Yc(Y,I,K);
		  double rightend = (p[c+1].Rend > -2) ?  Yc(Y,RI,RK) + scale*(p[c+1].orientation ? X[M+1-RJ] : X[M+1] - X[RJ]) : Yc(Y,RI,RK);
		  if((leftend < CenterMapped) || (rightend > Y[N+1] - CenterMapped))
		    break;
		}
		if(c < colors)
		  continue;
	      }

	      if(!split_maps){
		if(DEBUG && BestRef) assert(!Xmap->paired);
		if(Xmap->paired)
		  continue;
	      }
	      Xmap->paired = 1;
	      goodmaps[numgoodmaps++] = Xmap;
	    }
	    if(DEBUG && !(numgoodmaps <= mapcnt)){
	      printf("numgoodmaps=%d,mapcnt=%d,nummaps=%d\n",numgoodmaps,mapcnt,nummaps);
	      fflush(stdout);
	      assert(numgoodmaps <= mapcnt);
	    }
	    if(!split_maps)
	      strcpy(MappedIdPrefix,MappedPrefix);
	    else
	      sprintf(MappedIdPrefix,"%s_contig%lld",MappedPrefix, refmap[refid]->id);

	    if(maptype==0 || MapType >= 0 || CmapMergeBnx)
	      output_bnx(MappedIdPrefix,goodmaps,0,numgoodmaps-1, 1, NULL, NULL,-1);
	    else
	      output_cmap(MappedIdPrefix,goodmaps,0,numgoodmaps-1);

	  }// refid = 0 .. numrefmaps-1
	  delete [] goodmaps;
	} // no groupfile

	if(VERB){
	  printf("Finished output of mapped bnx files: CPU time=%0.6f, wall time=%0.6f\n", mtime(),wtime());
	  fflush(stdout);
	}
      }

      if(UnMappedPrefix){/* output maps that scored below threshold to <UnMappedPrefix>.cmap */
	char UnMappedIdPrefix[PATH_MAX];

	/* assemble all good maps into array badmap[0..badmapcnt-1] */
	Cmap **badmaps = new Cmap *[nummaps];
	for(int i = 0; i < nummaps; i++)
	  Gmap[i]->paired = 0;

	for(int refid = 0; refid < numrefmaps; refid++){
	  size_t align_start = 0, align_end = numaligns;
	  for(int i = 0; i < nummaps; i++)
	    Gmap[i]->paired = 0;
	  if(split_maps){
	    align_start = numalign_start[refid];
	    align_end = numalign_end[refid];
	  }

	  int numbadmaps = 0;

	  for(size_t i = align_start; i < align_end; i++){
	    Calign *p = alignment[i];
	    if(DEBUG && p && split_maps) assert(p->mapid1 == refid);
	    if(!p || p->numpairs <= 0)
	      continue;
	    Cmap *rmap = refmap[p->mapid1];
	    if(!AlignedThreshold(p, rmap->site, ScoreThreshold, LogPvThreshold))
	      continue;
	    Cmap *Xmap = Gmap[p->mapid2];
	    if(Xmap->origmap)
	      continue;/* skip map fragments : the best alignment of a map fragment should always be worse than that of the full molecules */
	    if(BestRef){
	      Cmap *origmap = Xmap;
	      while(origmap->origmap)
		origmap = origmap->origmap;
	      if(origmap->align->mapid1 != p->mapid1)
		continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	    }
	    //	    if(BestRef && Xmap->align->mapid1 != p->mapid1)
	    //	      continue;/* skip alignment unless it is with the reference with the best alignment score with this map Xmap */

	    if(!split_maps){
	      if(DEBUG && BestRef) assert(!Xmap->paired);
	      if(Xmap->paired)
		continue;
	    }
	    Xmap->paired = 1;
	  }
	  for(int i = 0; i < nummaps; i++){
	    Cmap *pmap = Gmap[i];
	    if(!pmap->paired)
	      badmaps[numbadmaps++] = pmap;
	  }
	  if(DEBUG) assert(numbadmaps <= nummaps);
	  if(!split_maps)
	    strcpy(UnMappedIdPrefix,UnMappedPrefix);
	  else
	    sprintf(UnMappedIdPrefix,"%s_contig%lld",UnMappedPrefix, refmap[refid]->id);

	  if(maptype==0 || MapType >= 0 || CmapMergeBnx)
	    output_bnx(UnMappedIdPrefix,badmaps,0,numbadmaps-1, 1, NULL, NULL,-1);
	  else
	    output_cmap(UnMappedIdPrefix,badmaps,0,numbadmaps-1);

	  if(!split_maps)
	    break;
	}
	delete [] badmaps;
      }

      if(PartMappedPrefix){/* output maps that scored above threshold but aligned only partly due to -endoutlier to <PartMappedPrefix>.cmap */
	char PartMappedIdPrefix[PATH_MAX];
	//	int split_maps = (MappedUnsplit >= 0) ? !MappedUnsplit : split_maps;/* If -mapped outputs are to be split */

	/* assemble all good maps with endoutliers into array partmaps[0..numpartmaps-1] */
	Cmap **partmaps = new Cmap *[nummaps];
	for(int refid = 0; refid < numrefmaps; refid++){
	  size_t align_start = 0, align_end = numaligns;
	  if(!split_maps)
	    for(int i = 0; i < nummaps; i++)
	      Gmap[i]->paired = 0;
	  else {
	    align_start = numalign_start[refid];
	    align_end = numalign_end[refid];
	  }

	  int numpartmaps = 0;

	  for(size_t i = align_start; i < align_end; i++){
	    Calign *p = alignment[i];
	    if(DEBUG && p && split_maps) assert(p->mapid1 == refid);
	    if(!p || p->numpairs <= 0)
	      continue;

	    Cmap *rmap = refmap[p->mapid1];
	    if(!AlignedThreshold(p, rmap->site, ScoreThreshold, LogPvThreshold))
	      continue;
	    if(p[1].Lend > -2 && p[1].Rend > -2 && p[2].Lend > -2 && p[2].Rend > -2)
	      continue;
	    
	    Cmap *Xmap = Gmap[p->mapid2];
	    if(Xmap->origmap)
	      continue;/* skip map fragments : the best alignment of a map fragment should always be worse than that of the full molecules */

	    if(BestRef){
	      Cmap *origmap = Xmap;
	      while(origmap->origmap)
		origmap = origmap->origmap;
	      if(origmap->align->mapid1 != p->mapid1)
		continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	    }
	    //	    if(BestRef && Xmap->align->mapid1 != p->mapid1)
	    //	      continue;/* skip alignment unless it is with the reference with the best alignment score with this map Xmap */

	    if(!split_maps){
	      if(DEBUG && BestRef) assert(!Xmap->paired);
	      if(Xmap->paired)
		continue;
	    }
	    Xmap->paired = 1;
	    partmaps[numpartmaps++] = Xmap;
	  }
	  if(DEBUG) assert(numpartmaps <= nummaps);
	  if(!split_maps)
	    strcpy(PartMappedIdPrefix,PartMappedPrefix);
	  else
	    sprintf(PartMappedIdPrefix,"%s_contig%lld",PartMappedPrefix, refmap[refid]->id);

	  if(maptype==0 || MapType >= 0 || CmapMergeBnx)
	    output_bnx(PartMappedIdPrefix,partmaps,0,numpartmaps-1, 1, NULL, NULL,-1);
	  else
	    output_cmap(PartMappedIdPrefix,partmaps,0,numpartmaps-1);

	  if(!split_maps)
	    break;
	}
	delete [] partmaps;
      }

      if(DEBUG>=2) assert(giter == RefRepeats-1);

      if(((Ch > 0.0 && ChFdr > 0.0) || (NEW && PoutlierEnd >= 0.0)) && giter2 == RefRepeats2 - 1){/* output XMAP and _r.cmap and _q.cmap */
	for(int refid = 0; refid < numrefmaps; refid++)
	  refmap[refid]->origmap = 0;/* reference maps are never chimeric fragments in refalign.cpp */

	extern Cmap **YYmap,**XXmap;
	extern int numY,numX;
	numY = numrefmaps;
	YYmap = refmap;
	numX = nummaps;
	XXmap = Gmap;

	/* need to create single CMAP file of input vfixx file (or multiple CMAP files) from possibly rescaled & filtered maps */
	char refid_prefix[PATH_MAX],cmapprefix[PATH_MAX],vfixx[PATH_MAX];
	
	/* assemble array of pointers to maps that have an alignment above all 3 thresholds into array goodmap[0..goodmapcnt-1] */
	Cmap **goodmaps = new Cmap *[mapcnt];

	for(int refid = 0; refid < numrefmaps; refid++){
	  size_t align_start = 0, align_end = numaligns;
	  for(int i = 0; i < nummaps; i++)
	    Gmap[i]->paired = 0;
	  if(split_maps) {
	    align_start = numalign_start[refid];
	    align_end = numalign_end[refid];
	  }
	  int numgoodmaps = 0;
	  for(size_t i = align_start; i < align_end; i++){
	    Calign *p = alignment[i];
	    if(DEBUG && p && split_maps) assert(p->mapid1 == refid);
	    if(!p || p->numpairs <= 0)
	      continue;
	    Cmap *rmap = refmap[p->mapid1];
	    if(!AlignedThreshold(p, rmap->site, ScoreThreshold, LogPvThreshold))
	      continue;
	    Cmap *Xmap = Gmap[p->mapid2];
	    //	    if(BestRef && Xmap->align->mapid1 != p->mapid1)
	    //	      continue;/* skip alignment unless it is with the reference with the best alignment score with this map Xmap */

	    Cmap *origmap = Xmap;
	    while(origmap->origmap)
	      origmap = origmap->origmap;
	    if(BestRef && origmap->align->mapid1 != p->mapid1)
	      continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */

	    if(!split_maps){
	      if(DEBUG && BestRef && NoSplit==2 && Xmap->paired){
		printf("refid=%d,align_start=%lu,align_end=%lu,i=%lu:p->mapid1=%d,p->mapid2=%d:Xmap->paired = %d (expected 0),origmap->align->mapid1=%d\n",
		       refid,align_start,align_end,i,p->mapid1,p->mapid2,Xmap->paired,origmap->align->mapid1);
		fflush(stdout);
		assert(!Xmap->paired);
	      }
	    }
	    if(origmap->paired)
	      continue;
	    origmap->paired = 1 + p->mapid1;
	    if(DEBUG>=2) assert(origmap->paired);
	    goodmaps[numgoodmaps++] = origmap;
	  }
	  if(DEBUG && !(numgoodmaps <= mapcnt)){
	    printf("refid=%d(%lld):numgoodmaps=%d,mapcnt=%d,nummaps=%d\n",refid,YYmap[refid]->id,numgoodmaps,mapcnt,orignummaps);
	    fflush(stdout);
	    assert(numgoodmaps <= mapcnt);
	  }

	  if(!split_maps)
	    strcpy(refid_prefix,output_prefix);
	  else
	    sprintf(refid_prefix,"%s_contig%lld",output_prefix,YYmap[refid]->id);

	  sprintf(cmapprefix,"%s_q",refid_prefix);

	  qsort(goodmaps,numgoodmaps,sizeof(Cmap *), (intcmp*)CmapIdInc);

	  if(CmapMergeBnx){
	    output_bnx(cmapprefix,goodmaps,0,numgoodmaps-1, 0,  NULL, NULL, -1);
	    sprintf(vfixx,"%s.bnx",cmapprefix);	  /* output_xmap() expects Query Map file name in vfixx_filename[0] */
	  } else {
	    output_cmap(cmapprefix,goodmaps,0,numgoodmaps-1);
	    sprintf(vfixx,"%s.cmap",cmapprefix);	  /* output_xmap() expects Query Map file name in vfixx_filename[0] */
	  }

	  int orig_num_files = num_files;
	  char *orig_vfixx_filename0 = vfixx_filename[0];
	  vfixx_filename[0] = strdup(vfixx);
	  num_files = 1;

	  if(VERB>=2){
	    printf("refid=%d:calling output_xmap:split_maps=%d,refid_prefix=%s\n",
		refid,split_maps,refid_prefix);
	    fflush(stdout);
          }

	  output_xmap2(refid_prefix,align_start,align_end);/* output alignments (above 3 thresholds) as .xmap file */

	  free(vfixx_filename[0]);
	  vfixx_filename[0] = orig_vfixx_filename0;
	  num_files = orig_num_files;

	  if(!split_maps)
	    break;
	} // refid = 0 .. numrefmaps-1

	delete [] goodmaps;
      }

      nummaps = orignummaps;
    } // if(giter == RefRepeats-1)

    if(DEBUG>=2) assert(0 <= giter && giter <= RefRepeats-1);

    if(NoStat){
      score_free(0,numrefmaps);
      continue;/* skip parameter display/estimation/output */
    }

    double startPixelLen = PixelLen;
    double FPorig[MAXCOLOR], FNorig[MAXCOLOR], SForig[MAXCOLOR], SDorig[MAXCOLOR], SRorig[MAXCOLOR], SEorig[MAXCOLOR];
    int fragcnt = 0;/* total sitecnt[] over all colors */
    for(int c = 0; c < colors; c++){
      fragcnt += sitecnt[c];
      FPorig[c] = FP[c];
      FNorig[c] = FN[c];
      SForig[c] = SF[c];
      SDorig[c] = SD[c];
      SRorig[c] = SR[c];
      SEorig[c] = SE[c];
    }

    size_t lcnt = mapcnt;

    if(VERB){/* display statistics */
      printf("refmaps=%d,maps=%d,totalmaps=%d,%d:\n",
	     numrefmaps,nummaps,startmaps,totalmaps);

      if(biasWT==0.0)
	printf("  %lu alignments:score=%0.4f,logPV=%0.4f,pairs=%0.2f,log10LR=%0.6f(LR=%0.6f= %0.6f,%0.6f,%0.6f)\n",
	       tcnt, (tcnt <= 0) ? 0.0 : scoresum/tcnt, (tcnt <= 0) ? 0.0 : logPVsum/tcnt,
	       (tcnt <= 0) ? 0.0 : ((double)NumpairsSum)/tcnt, (tcnt <= 0) ? 0.0 : LRsum/(tcnt*log(10.0)), (tcnt <= 0) ? 0.0 : LRsum/tcnt, 
	       (tcnt <= 0) ? 0.0 : LRsumC[0]/tcnt, (tcnt <= 0) ? 0.0 : LRsumC[1]/tcnt, (tcnt <= 0) ? 0.0 : LRsumC[2]/tcnt);  
      else
	printf("  %lu alignments:score=%0.4f,logPV=%0.4f\n",
	       tcnt, (tcnt <= 0) ? 0.0 : scoresum/tcnt, (tcnt <= 0) ? 0.0 : logPVsum/tcnt);

      printf("  %d above threshold:score=%0.4f,logPV=%0.4f,pairs=%0.2f,logLR=%0.6f(%0.6f/site)\n",
	     mapcnt, (mapcnt <= 0) ? 0.0 : ATscoresum/mapcnt, (mapcnt <= 0) ? 0.0 : ATlogPVsum/mapcnt,
	     (mapcnt <= 0) ? 0.0 : ((double)ATNumpairsSum)/mapcnt, (mapcnt <= 0) ? 0.0 : ATlogLRsum/mapcnt,
	     (fragcnt-mapcnt <= 0) ? 0.0 : ATlogLRsum/((fragcnt-mapcnt)));
      if(Poutlier>0.0)
	printf("  outliers=%d/%d (rate=%0.6f), endoutliers=%d/%d (rate=%0.6f)\n",
	       outliercnt,outliercnt+segcnt[0]+segcnt[1],outliersum/max(1,outliertot),outlierEndcnt,outliercnt+outlierEndcnt+segcnt[0]+segcnt[1], outlierEndcnt/((double)(outliercnt+outlierEndcnt+segcnt[0]+segcnt[1])));
      if((Ch > 0.0 && ChFdr > 0.0) || (NEW && PoutlierEnd > 0.0))
	printf("  chimeric sites=%d/%d(%0.6f),chimeric frags=%d(FP=%d,close pairs=%d)):\n",
	       chimcnt/colors,fragcnt/colors,((double)chimcnt)/fragcnt,chimconf,chimconfFP,chimpaircnt);
      fflush(stdout);
    }

    if(DEBUG>=2) assert(0 <= giter && giter <= RefRepeats-1);

    parameters[giter].mapcnt = tcnt;
    parameters[giter].logLR = (tcnt <= 0 ? 0.0 : LRsum/tcnt);
    parameters[giter].ATmapcnt = mapcnt;
    parameters[giter].ATlogLR = (mapcnt <= 0 ? 0.0 : ATlogLRsum/mapcnt);
    parameters[giter].sitecnt = sitecnt[0]+sitecnt[1];
    parameters[giter].outlierRate = outliersum/max(1,outliertot);
    parameters[giter].EndoutlierRate = outlierEndcnt/((double)(outliercnt+outlierEndcnt + segcnt[0]+segcnt[1]));
    for(int c = 0; c < colors; c++)
      parameters[giter].LabelDensity[c] = sitecnt[c]*100.0/len[c];

    if(REFDEBUG && logLRsum > LARGE_NEGATIVE && (!(ATlogLRsum > logLRsum -(USE_RFLOAT ? 1e-4 : 1e-6) * mapcnt) || VERB)){
      if(SDDEBUG)
	printf("giter=%d:After alignment update:score=%0.6f -> %0.6f,mapcnt=%llu(logLRsd= %0.6f)\n",
	       giter,logLRsum/max(1,mapcnt),ATlogLRsum/max(1,mapcnt),(unsigned long long)mapcnt, logLRsdsum/max(1,mapcnt));
      else
	printf("giter=%d:After alignment update:score=%0.6f -> %0.6f,mapcnt=%llu\n",giter,logLRsum/max(1,mapcnt),ATlogLRsum/max(1,mapcnt),(unsigned long long)mapcnt);
      if(HashGen && !(ATlogLRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6) * mapcnt))
	printf("WARNING: Can happen with -hashgen\n");
      fflush(stdout);

      if(!HashGen) assert(ATlogLRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6) * mapcnt);
    }
    logLRsum = ATlogLRsum;
    if(DEBUG) assert(isfinite(logLRsum));

    for(int c = 0; c < colors; c++){
      /* update estimates of PRbias */
      if(sitecnt[c]-fn[c]+fp[c] > 0){
	if(VERB){
	  printf("  c=%d:PrPen=%0.6f,Xsites=%0.1f: PRbias = %0.6f -> %0.6f\n",
		 c, PrPen[c], sitecnt[c]-fn[c]+fp[c], PRbias[c], PrPen[c]/(sitecnt[c]-fn[c]+fp[c]));
	  fflush(stdout);
	}
	if(0)// TODO : enable PRbias updates
	  PRbias[c] = PrPen[c]/(sitecnt[c]-fn[c]+fp[c]);
      }
      /* update estimates of FP[] */
      if(len[c]*FP[c] > 0.0){
	if(VERB){
	  if(sitecnt[c]<=0)
	    printf("  c=%d:%s=%0.3f(%0.3f),%0.3f w/outliers,FPcnt=%d(isolated=%d) : FP = %0.6f -> %0.6f(unbiased=%0.6f) (per 100kb), FPrate=%0.6f\n",
		   c,REFSCORE?"sumX":"sumY",len[c],lenB[c],lenY[c],fp[c],fpI[c], FP[c], fp[c]*100.0/len[c], fp[c]*100.0/lenB[c],0.0);
	  else
	    printf("  c=%d:%s=%0.3f(%0.3f),%0.3f w/outliers,FPcnt=%d(isolated=%d) : FP = %0.6f -> %0.6f(unbiased=%0.6f) (per 100kb), FPrate=%0.6f\n",
		   c,REFSCORE?"sumX":"sumY",len[c],lenB[c],lenY[c],fp[c],fpI[c], FP[c], fp[c]*100.0/len[c], fp[c]*100.0/lenB[c], ((double)fp[c])/sitecnt[c]);
	  fflush(stdout);
	}
	FP[c] = fp[c]*100.0/len[c];
      }
      parameters[giter+1].FPfrac[c] = 0.0;
      if(sitecnt[c] > 0)
	parameters[giter+1].FPfrac[c] = ((double)fp[c])/sitecnt[c];
    }
    if(VERB>=2){
      printf("giter=%d:REFDEBUG=%d,numrefmaps=%d,BestRef=%d,BestRefPV=%d\n",giter, REFDEBUG,numrefmaps,BestRef,BestRefPV);
      fflush(stdout);
    }

    if(REFDEBUG && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
      if(DEBUG) assert(logLRarray);
      score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
      double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
      double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
      lcnt= 0;
      double *origLRarray = new double[nummaps];
      double *origSDarray[MAXCOLOR];
      for(int mid = 0; mid < nummaps; mid++)
	origLRarray[mid] = logLRarray[mid];
      for(int c = 0; c < colors; c++){
	origSDarray[c] = new double[nummaps];
	for(int mid = 0; mid < nummaps; mid++){
	  origSDarray[c][mid] = logSDarray[c][mid];
	}
      }

      double origLRsum = 0.0;
      double origSDsum[MAXCOLOR] = {0.0,0.0};

      for(size_t i = 0; i < numaligns; i++){
	double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	Calign *align = alignment[i];
	if(!align || align->numpairs <= 1)
	  continue;
	int mid = align->mapid2;
	if(DEBUG) assert(!Gmap[mid]->origmap);
	if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
	  if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
	  continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	}
	if(align->score <= origScoreThreshold)
	  continue;
	int rid = align->mapid1;
	if(DEBUG>=2){
	  int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
	  assert(origgoodalign);
	  int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
	  assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	}
	FLOAT **XX = Gmap[mid]->site;
	int *MM = Gmap[mid]->numsite;
	FLOAT **YY = refmap[rid]->site;
	//	int *NN = Gmap[rid]->numsite;
	align++;
	logLRarray[mid] = 0.0;
	for(int c = 0; c < colors; c++){
	  FLOAT *X = XX[c];
	  double M = MM[c];
	  FLOAT *Y = YY[c];
	  int N = refmap[rid]->numsite[c];
	  int U = align[c].numpairs;

	  LR[c] = logLR(mid,rid,c,X,Y,M,N,align,align[c].sites1[0], align[c].sites2[0], align[c].sitesK1[0], align[c].sites2[U-1], align[c].Lij1,align[c].Rij1,align[c].Lij2,align[c].Rij2,1,0,SDDEBUG ? &LRsd[c] : 0,NULL);
	  logLRarrayC[c][mid] = LR[c];
	  logSDarray[c][mid] = LRsd[c];
	  if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
	    continue;
	  LRsum += LR[c];
	  if(SDDEBUG){
	    LRsdsum += LRsd[c];
	    LRsdsumC[c] += LRsd[c];
	  }
	}
	int U = align[0].numpairs;
	int I = align[0].sites1[U-1];
	int K = align[0].sitesK1[U-1];
	int J = align[0].sites2[U-1];
	U = align[1].numpairs;
	int P = align[1].sites1[U-1];
	int T = align[1].sitesK1[U-1];
	int Q = align[1].sites2[U-1];
	double y0 = Yc(YY[0],I,K);
	double y1 = Yc(YY[1],P,T);
	double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);
	
	LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
	if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	  LRsum += LR[2];
	  for(int lc = 0; lc <= colors; lc++)
	    LRsumC[lc] += LR[lc];
	}

	logLRarray[mid] = LR[0] + LR[1] + LR[2];
	lcnt++;
	if(VERB>=2 && giter==0){
	  if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	    origLRsum += origLRarray[mid];
	    for(int lc = 0; lc < colors; lc++)
	      origSDsum[lc] += origSDarray[lc][mid];
	  }
	  if(mid==265){
	    printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)(delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f)\n",
		   i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
		   logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
	    printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
		   logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
		   LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
		   LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
	    printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);
	  }
	}
      }
      if(SDDEBUG)
	printf("After FP update:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f= %0.6f,%0.6f)\n",
	       (unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt,LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
      else
	printf("After FP update:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
      fflush(stdout);

      if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	if(VERB){
	  double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
	  for(size_t i = 0; i < numaligns; i++){
	    Calign *align = alignment[i];
	    if(!align || align->numpairs <= 1)
	      continue;
	    int mid = align->mapid2;
	    if(DEBUG>=2) assert(!Gmap[mid]->origmap);
	    if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
	      continue;
	    if(MapRate > 0.0 && !(align->score > ScoreThreshold))
	      continue;
	    sum1 += origLRarray[mid];
	    sum2 += logLRarray[mid];
	    sum3 += origSDarray[0][mid] + origSDarray[1][mid];
	    sum4 += logSDarray[0][mid] + logSDarray[1][mid];
	    printf("i=%lu/%lu:rid=%d,mid=%d,id=%lld,or=%d:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
		   i, numaligns, align->mapid1, mid, Gmap[mid]->id, align->orientation, origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
		   origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
		   logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
	  }
	  fflush(stdout);
        }
	      
	assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
      }

      if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	assert(lcnt == tcnt);
      }
      logLRsum = LRsum;
      if(SDDEBUG)
	logLRsdsum = LRsdsum;
      score_free(0,numrefmaps);

      delete [] origLRarray;
      for(int c = 0; c < colors; c++){
	delete [] origSDarray[c];
      }
    }

    /* update estimates of FN[] */
    for(int c = 0; c < colors; c++){
      if(sitecnt[c]*FN[c] > 0.0){
	double origFN = FN[c];
	FN[c] = ((double)fn[c])/((double)sitecnt[c]);
	if(VERB){
	  printf("  c=%d:sitecnt=%d,FNcnt=%d(isolated=%d) : FN =%0.6f -> %0.6f\n",c,sitecnt[c],(int)floor(fn[c]+0.5),fnI[c],origFN, FN[c]);
	  fflush(stdout);
	}

	FN[c] = min(FN[c],MAXFN);
	FN[c] = max(FN[c],MINFN);

	/* if SCORE_APPROX <= 1 && MIS_VITERBI==0 this simple estimate of FN may underestimate : use Golden Mean search and calls to logLR() to correct it */
	if(FN_GM){/* use goldenmean search and calls to totLL() to estimate FN */
	  origFN = FN[c];

	  #ifdef _OPENMP
	  numthreads = omp_get_max_threads();
	  numthreads = min(numthreads,MaxThreads);
	  if(numthreads > mapcnt)
	    numthreads = max(1,mapcnt);
	  #endif

	  double bestLL = FNtotLL(FN[c], c, numthreads, numaligns, alignment, Tarray1);
	  double origLL = bestLL;

	  double phi = (1.0 + sqrt(5.0))*0.5;
	  double resphi = 2.0 - phi;
	
	  double Low = MINFN;
	  double High = MAXFN;
	  double Mid = FN[c];
	  if(VERB>=2){
	    printf("FN[%d]=%0.6f:LL=%0.6f, Low=%0.6f,High=%0.6f\n", c,FN[c], bestLL, Low, High);
	    fflush(stdout);
	  }

	  while(max(High - Mid, Mid - Low) > 0.000001){
	    double Lsum;
	    if(High - Mid > Mid - Low){
	      FN[c] = Mid + resphi * (High - Low);
	      Lsum = FNtotLL(FN[c],c,numthreads,numaligns,alignment,Tarray1);

	      if(Lsum > bestLL){
		bestLL = Lsum;
		Low = Mid;
		Mid = FN[c];
	      } else
		High = FN[c];
	    } else {
	      FN[c] = Mid - resphi * (High - Low);
	      Lsum = FNtotLL(FN[c],c,numthreads,numaligns,alignment,Tarray1);

	      if(Lsum > bestLL){
		bestLL = Lsum;
		High = Mid;
		Mid = FN[c];
	      } else 
		Low = FN[c];
	    }
	    if(VERB>=2){
	      printf("FN[%d]=%0.6f:LL=%0.6f (best:FN=%0.6f,LL=%0.6f), Low=%0.6f,High=%0.6f\n",c, FN[c], Lsum, Mid, bestLL, Low, High);
	      fflush(stdout);
	    }
	  }
	  FN[c] = Mid;
	  if(VERB){
	    printf("\tGolden Mean Search: FN[%d]= %0.6f -> %0.6f (score= %0.6f -> %0.6f, lcnt= %lu)\n", c, origFN, FN[c], origLL/max((size_t)1,lcnt), bestLL/max((size_t)1,lcnt), lcnt );
	    fflush(stdout);
	  }
	}
      }
    }
    if(REFDEBUG && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
      if(DEBUG) assert(logLRarray);
      score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
      double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
      double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
      lcnt= 0;
      double *origLRarray = new double[nummaps];
      double *origSDarray[MAXCOLOR];
      for(int mid = 0; mid < nummaps; mid++)
	origLRarray[mid] = logLRarray[mid];
      for(int c = 0; c < colors; c++){
	origSDarray[c] = new double[nummaps];
	for(int mid = 0; mid < nummaps; mid++){
	  origSDarray[c][mid] = logSDarray[c][mid];
	}
      }

      double origLRsum = 0.0;
      double origSDsum[MAXCOLOR] = {0.0,0.0};

      for(size_t i = 0; i < numaligns; i++){
	double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	Calign *align = alignment[i];
	if(!align || align->numpairs <= 1)
	  continue;
	int mid = align->mapid2;
	if(DEBUG) assert(!Gmap[mid]->origmap);
	if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
	  if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
	  continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	}
	if(align->score <= origScoreThreshold)
	  continue;
	int rid = align->mapid1;
	if(DEBUG>=2){
	  int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
	  assert(origgoodalign);
	  int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
	  assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	}
	FLOAT **XX = Gmap[mid]->site;
	int *MM = Gmap[mid]->numsite;
	FLOAT **YY = refmap[rid]->site;
	//	int *NN = Gmap[rid]->numsite;
	align++;
	logLRarray[mid] = 0.0;
	for(int c = 0; c < colors; c++){
	  FLOAT *X = XX[c];
	  double M = MM[c];
	  FLOAT *Y = YY[c];
	  int N = refmap[rid]->numsite[c];
	  int U = align[c].numpairs;

	  LR[c] = logLR(mid,rid,c,X,Y,M,N,align,align[c].sites1[0], align[c].sites2[0], align[c].sitesK1[0], align[c].sites2[U-1], align[c].Lij1,align[c].Rij1,align[c].Lij2,align[c].Rij2,1,0, SDDEBUG ? &LRsd[c] : 0,NULL);
	  logLRarrayC[c][mid] = LR[c];
	  logSDarray[c][mid] = LRsd[c];
	  if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
	    continue;
	  LRsum += LR[c];
	  if(SDDEBUG){
	    LRsdsum += LRsd[c];
	    LRsdsumC[c] += LRsd[c];
	  }
	}

	int U = align[0].numpairs;
	int I = align[0].sites1[U-1];
	int K = align[0].sitesK1[U-1];
	int J = align[0].sites2[U-1];
	U = align[1].numpairs;
	int P = align[1].sites1[U-1];
	int T = align[1].sitesK1[U-1];
	int Q = align[1].sites2[U-1];
	double y0 = Yc(YY[0],I,K);
	double y1 = Yc(YY[1],P,T);
	double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);
	
	LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
	if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	  LRsum += LR[2];
	  for(int lc = 0; lc <= colors; lc++)
	    LRsumC[lc] += LR[lc];
	}

	logLRarray[mid] = LR[0] + LR[1] + LR[2];
	lcnt++;
	if(VERB>=2 && giter==0){
	  if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	    origLRsum += origLRarray[mid];
	    for(int lc = 0; lc < colors; lc++)
	      origSDsum[lc] += origSDarray[lc][mid];
	  }
	  if(mid==265){
	    printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)(delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f)\n",
		   i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
		   logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
	    printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
		   logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
		   LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
		   LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
	    printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);
	  }
	}
      }
      if(SDDEBUG)
	printf("After FN update:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f= %0.6f,%0.6f)\n",
	       (unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt, LRsdsumC[0]/lcnt, LRsdsumC[1]/lcnt);
      else
	printf("After FN update:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
      fflush(stdout);

      if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	if(VERB){
	  double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
	  for(size_t i = 0; i < numaligns; i++){
	    Calign *align = alignment[i];
	    if(!align || align->numpairs <= 1)
	      continue;
	    int mid = align->mapid2;
	    if(DEBUG>=2) assert(!Gmap[mid]->origmap);
	    if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
	      continue;
	    if(MapRate > 0.0 && !(align->score > ScoreThreshold))
	      continue;
	    sum1 += origLRarray[mid];
	    sum2 += logLRarray[mid];
	    sum3 += origSDarray[0][mid] + origSDarray[1][mid];
	    sum4 += logSDarray[0][mid] + logSDarray[1][mid];
	    printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
		   i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
		   origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
		   logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
	  }
	  fflush(stdout);
        }
	      
	assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
      }

      if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	assert(lcnt == tcnt);
      }
      logLRsum = LRsum;
      if(SDDEBUG)
	logLRsdsum = LRsdsum;
      score_free(0,numrefmaps);

      delete [] origLRarray;
      for(int c = 0; c < colors; c++){
	delete [] origSDarray[c];
      }
    }

    if(numerrors[0]+numerrors[1]+numerrors[2] > 0){

      double A[MAXCOLOR+1],B[MAXCOLOR+1],R[MAXCOLOR+1],E[MAXCOLOR+1],MA[MAXCOLOR+1];
      for(int c = 0; c <= colors; c++){
	A[c] = SF[c]*SF[c];
	B[c] = fabs(SD[c])*SD[c];
	MA[c] = 0.0;
	if(c < colors){
	  R[c] = SR[c]*SR[c];
	  E[c] = SE[c]*SE[c];
	} else
	  R[c] = E[c] = 0.0;
      }
      MA[colors] = MA_mean;

      double Amin = MINSF*MINSF, Amax = MAXSF*MAXSF;
      double Bmin = fabs(MINSD)*MINSD, Bmax = MAXSD*MAXSD;
      double Rmin = MINSR*MINSR, Rmax = MAXSR*MAXSR;
      double Emin = MINSE*MINSE, Emax = MAXSE*MAXSE;

      if(REFDEBUG_STRICT < 2 && maxresbias > mres * 0.5){

	if(VERB/* HERE >=2 */){
	  printf("Estimating -resbias parameters:wall time=%0.6f\n",wtime());
	  fflush(stdout);
	}

	/* save original bias values, incase we need to backtrack */
	double startresbiasX[MAXCOLOR][RESBINS+1],startresbias[MAXCOLOR][RESBINS+1];
	int startResBins[MAXCOLOR];

	Imaxresbias = (int)floor(maxresbias*1000.0+0.5);

	for(int c = 0; c < colors; c++){
	  startResBins[c] = ResBins[c];
	  for(int Bin = 0; Bin <= ResBins[c]; Bin++){
	    startresbiasX[c][Bin] = resbiasX[c][Bin];
	    startresbias[c][Bin] = resbias[c][Bin];
	  }

	  ResBins[c] = origResBins[c];/* in case last iteration reduced ResBins */
	  if(VERB>=2){
	    printf("giter=%d,c=%d:origResBins[c]=%d,startResBins[c]=%d\n",giter,c,origResBins[c],startResBins[c]);
	    fflush(stdout);
	  }

	  SizeToBin[c] = new int[Imaxresbias + 1];
	}

	if(DEBUG) assert(!NoStat);

	if(rawsitemaps < totalmaps){/* in case maps were split */
	  rawsitealloc(Gmap,rawsitemaps,totalmaps);
	  rawsitemaps = totalmaps;
	}

#ifdef _OPENMP
	numthreads = omp_get_max_threads();
	numthreads = min(numthreads,MaxThreads);
	if(numthreads > max(numerrors[0],numerrors[1])/64)
	  numthreads = max(1,max(numerrors[0],numerrors[1])/64);
#endif

	int maxindex[MAXCOLOR];
        int Start[MAXCOLOR][RESBINS], End[MAXCOLOR][RESBINS];
	double meanbias[MAXCOLOR][RESBINS+1];

	/* precompute rawx, lograwx */
	for(int c = 0; c < colors; c++){

	  if(numerrors[c] <= 0){      /* set reasonable default values corresponding to no bias */
	    ResBins[c] = 1;
	    resbiasX[c][0] = mres * 0.5;
	    resbiasX[c][1] = maxresbias;
	    for(int Bin = 0; Bin <= ResBins[c]; Bin++)
	      resbias[c][Bin] = 0.0;

	    for(int i = 0; i <= Imaxresbias; i++)
	      SizeToBin[c][i] = 0;

	    continue;
	  }

  	  #pragma omp parallel for num_threads(numthreads) schedule(static,64)
	  for(int i = 0; i < numerrors[c]; i++){
	    Cinterval *perr = &errors[c][i];
	    if(ENDFIX>=2 && perr->end){
	      perr->rawx = maxresbias + 10.0;
	      continue;
	    }
	    double *rawsite = Gmap[perr->mapid]->rawsite[c];
	    if(DEBUG && !rawsite){
	      #pragma omp critical
	      {
	        printf("errors[%d]:mapid=%d,Gmap[mapid]->rawsite=%p\n",i,perr->mapid,rawsite);
		fflush(stdout);
		assert(rawsite != 0);
	      }
	    }
	    perr->rawx = rawsite[perr->R] - rawsite[perr->L];

	    if(DEBUG>=2 && !(perr->rawx > 0.0)){
	      #pragma omp critical
	      {
		printf("perror[%d]:mapid=%d(id=%lld,M=%d),L=%d,R=%d:rawsite[L]=%0.4f,rawsite[R]=%0.4f,rawx=%0.4f\n",
		       i,perr->mapid,Gmap[perr->mapid]->id,Gmap[perr->mapid]->numsite[c],perr->L,perr->R,rawsite[perr->L],rawsite[perr->R],perr->rawx);
		fflush(stdout);
		assert(perr->rawx > 0.0);
	      }
	    }
	  }

	  /* estimate bias parameters */
	
	  /* first sort intervals in ascending order of rawx */
	  qsort(errors[c],numerrors[c],sizeof(Cinterval),(intcmp*)ErrorRawxInc);

	  /* locate largest size <= maxresbias */
	  maxindex[c] = 0;
	  for(; maxindex[c] < numerrors[c]; maxindex[c]++)
	    if(errors[c][maxindex[c]].rawx > maxresbias)
	      break;
	  maxindex[c]--;
	  if(DEBUG) assert(maxindex[c] < numerrors[c]);
	
	  /* initialize resbiasX[c][0..ResBins[c]] : assume equal interval size per Bin */
	  if(maxindex[c] < ResBins[c] * RESBIAS_MINSAMPLE){
	    if(VERB){
	      printf("Reduced ResBins[%d] from %d to %d due to limited total number of samples = %d\n",
		     c,ResBins[c], max(1,maxindex[c]/RESBIAS_MINSAMPLE),maxindex[c]);
	      fflush(stdout);
	    }
	    ResBins[c] = max(1,maxindex[c]/RESBIAS_MINSAMPLE);
	  }
	  resbiasX[c][ResBins[c]] = maxresbias;
	  if(numerrors[c] > 0){
	    resbiasX[c][0] = errors[c][0].rawx;
	    if(DEBUG) assert(resbiasX[c][0] >= 0.0);
	  } else
	    resbiasX[c][0] = 0.0;
	  Start[c][0] = 0;
	  double BinWidth = (resbiasX[c][ResBins[c]]-resbiasX[c][0])/((double)ResBins[c]);
	  if(VERB/* HERE >=2 */ && numerrors[c] > 0){
	    printf("c=%d:maxresbias=%0.4f,ResBins=%d,errors[0].rawx=%0.4f(mapid=%d,id=%lld),numerrors=%d,maxindex=%d\n",c,maxresbias,ResBins[c],errors[c][0].rawx,errors[c][0].mapid,Gmap[errors[c][0].mapid]->id,numerrors[c],maxindex[c]);
	    printf("   BinWidth = %0.4f\n", BinWidth);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(isfinite(BinWidth) && BinWidth > 0.0);
	  for(int Bin = 1; Bin < ResBins[c]; Bin++){
	    double LB = resbiasX[c][Bin] = min(maxresbias,resbiasX[c][Bin-1] + BinWidth);
	    if(DEBUG>=2) assert(resbiasX[c][Bin] <= maxresbias);

	    /* locate the start of this Bin */
	    int BinStart = Start[c][Bin-1];
	    for(;BinStart <= maxindex[c]; BinStart++)
	      if(errors[c][BinStart].rawx >= LB)
		break;
	    if(BinStart < Start[c][Bin-1] + RESBIAS_MINSAMPLE)/* make sure each bin has at least RESBIAS_MINSAMPLE data points */
	      BinStart = Start[c][Bin-1] + RESBIAS_MINSAMPLE;

	    /* avoid placing boundary point between two almost identical sample values */
	    if(BinStart > 0)
	      for(;BinStart < maxindex[c]; BinStart++)
		if(errors[c][BinStart-1].rawx + 1e-8 < errors[c][BinStart].rawx)
		  break;
	    if(BinStart >= maxindex[c]+1 || (errors[c][BinStart].rawx >= maxresbias)){
	      if(VERB){
		printf("c=%d:Reduced ResBins from %d to %d due to reaching maxresbias=%0.4f in Bin %d(possibly limited number of quantized values)\n",
 	               c,ResBins[c], Bin-1, maxresbias, Bin-1);
		if(VERB>=2)
		  printf("  Bin=%d,Start[Bin-1]=%d,LB=%0.4f,BinStart=%d,errors[BinStart].rawx=%0.4f,maxresbias=%0.4f,maxindex=%d\n",
			 Bin,Start[c][Bin-1],LB,BinStart,errors[c][BinStart].rawx,maxresbias,maxindex[c]);
		fflush(stdout);
              }
	      resbiasX[c][Bin-1] = resbiasX[c][ResBins[c]];
	      ResBins[c] = Bin-1;
	    } else {
	      if(DEBUG) assert(errors[c][BinStart].rawx >= resbiasX[c][Bin]);
	      resbiasX[c][Bin] = min(maxresbias,0.5*(errors[c][BinStart-1].rawx + errors[c][BinStart].rawx));/* avoid placing Bin boundary close to actual sample value */
	      if(VERB>=2 || (DEBUG && !(resbiasX[c][Bin] > resbiasX[c][Bin-1]))){
		printf("c=%d:resbiasX[c][%d]=%0.8f : BinStart=%d, errors[c][BinStart-1].rawx=%0.8f, errors[c][BinStart].rawx=%0.8f\n",c,Bin,resbiasX[c][Bin],BinStart,errors[c][BinStart-1].rawx,errors[c][BinStart].rawx);
		fflush(stdout);
		if(DEBUG) assert(resbiasX[c][Bin] > resbiasX[c][Bin-1]);
	      }
	      Start[c][Bin] = BinStart;
	      if(DEBUG) assert(BinStart < maxindex[c]+1);
	    }
	    if(VERB>=2){
	      printf("c=%d:Bin=%d:resBiasX[c][Bin]=%0.4f,Start[c][Bin]=%d(maxresbias=%0.4f,maxindex[c]=%d,ResBins[c]=%d)\n",c,Bin,resbiasX[c][Bin],Start[c][Bin],maxresbias,maxindex[c],ResBins[c]);
	      fflush(stdout);
	    }
	    if(DEBUG && Bin < ResBins[c]) assert(resbiasX[c][Bin] > resbiasX[c][Bin-1]);
	  }
	  End[c][ResBins[c]-1] = maxindex[c]+1;
	  for(int Bin= ResBins[c]-1; --Bin >= 0;)
	    End[c][Bin] = Start[c][Bin+1];

	  if(VERB>=2 && nummaps > 0){
	    printf("giter=%d,c=%d:While estimating bias parameters(2):Gmap[0]->numsite[c]=M=%d,Gmap[0]->site[c][M]=%0.3f,Gmap[0]->rawsite[c][1]=%0.3f,Gmap[0]->rawsite[c][M]=%0.3f\n",
		   giter,c,Gmap[0]->numsite[c],Gmap[0]->site[c][Gmap[0]->numsite[c]],Gmap[0]->rawsite[c][1],Gmap[0]->rawsite[c][Gmap[0]->numsite[c]]);
	    printf("    resbiasX[c][ResBins[c]=%d]=%0.4f,maxresbias=%0.4f\n",ResBins[c],resbiasX[c][ResBins[c]],maxresbias);
	    fflush(stdout);
	  }
	  if(DEBUG>=2) assert(resbiasX[c][ResBins[c]] <= maxresbias);
	
#if SIMPLE_BIAS==2 // weighted piecewise linear fit 
	  /* For k = 0..ResBins[c], compute tri-diagonal matrix coeficients:
	     AA[k] = Sum(i=Start[c][k-1]..End[c][k-1]-1) (1-T(i,k))*T(i,k)/var(i) (With AA[0] = 0)
	     BB[k] = Sum(i=Start[c][k-1]..End[c][k-1]-1) T(i,k)^2/var(i)   +   Sum(i=Start[c][k]..End[c][k]-1) U(i,k)^2/var(i)
	     CC[k] = Sum(i=Start[c][k]..End[c][k]-1) (1-U(i,k))*U(i,k)/var(i)     (With C[ResBins[c]] = 0)
	     DD[k] = Sum(i=Start[c][k-1]..End[c][k-1]-1) (x(i)-y(i))*T(i,k)/var(i)   +   Sum(i=Start[c][k]..End[c][k]-1) (x(i)-y(i))U(i,k)/var(i)

	     Where:
	     T(i,k) = (x(i) - resbiasX[c][k-1])/(resbiasX[c][k] - resbiasX[c][k-1]), for k = 1..ResBins[c] (T(i,0) == 0)
	     U(i,k) = (resbiasX[c][k+1] - x(i))/(resbiasX[c][k+1] - resbiasX[c][k]), for k = 0..ResBins[c]-1 (U(i,ResBins[c]) == 0)
	     x(i) == errors[c][i].rawx
	     y(i) == errors[c][i].y
	     var(i) == A[c] + B[c] * y(i) + R[c]*y(i)^2 + E[c] * errors[i]->resvar

	     Then solve the tri-diagnonal set of equations for resbias[c][k = 0..ResBins[c]] as follows:
	     
	     nC[0] = CC[0]/BB[0]
	     nC[k] = CC[k]/(BB[k] - AA[k] * nC[k-1]), k = 1..ResBins[c]
	     nD[0] = DD[0]/BB[0]
	     nD[k] = (DD[k] - AA[k] * nD[k-1])/(BB[k] - AA[k] * nC[k-1]), k = 1 ... ResBins[c]

	     resbias[c][ResBins[c]] = nD[ResBins[c]]
	     resbias[c][k] = nD[k] - nC[k] * resbias[c][k+1], k = ResBins[c]-1 ... 0

	     Then compute corrected errors[i].x as follows :

	     For k = 0...ResBins[c]-1, 
	       For i = Start[c][k]..End[c][k] -1
	         errors[c][i].x = x(i) - (resbias[c][k] + (resbias[c][k+1]-resbias[c][k])(x(i)-resbiasX[c][k])/(resbiasX[c][k+1]-resbiasX[c][k]))
	   */
	  if(1){	  /* compute coeficients of tridiagonal linear equations */
	    double AA[RESBINS+1], BB[RESBINS+1], CC[RESBINS+1], DD[RESBINS+1];
	    for(int k = 0; k <= ResBins[c]; k++){
	      AA[k] = BB[k] = CC[k] = DD[k] = 0.0;
	      meanbias[c][k] = 0.0;
	      if(k > 0){
		if(DEBUG) assert(resbiasX[c][k] > resbiasX[c][k-1]);
		double Tinv = 1.0/(resbiasX[c][k] - resbiasX[c][k-1]);
		for(int i = Start[c][k-1]; i < End[c][k-1]; i++){// HERE : multithread
		  Cinterval *perr = &errors[c][i];
		  double x = perr->rawx;
		  double y = perr->y;
		  double var = A[c] + B[c] * y;
		  if(QUADRATIC_VARIANCE)
		    var += R[c]*y*y;
		  if(RES_VARIANCE)
		    var += E[c] * perr->resvar;
		  if(DEBUG>=2) assert(isfinite(var) && var > 0.0);
		  double Tik = (x - resbiasX[c][k-1])*Tinv;
		  double TIvar = Tik/var;
		
		  AA[k] += (1.0 - Tik) * TIvar;
		  BB[k] += Tik * TIvar;
		  DD[k] += (x - y) * TIvar;
		}
	      }
	      if(k < ResBins[c]){
		double Wsum = 0.0;
		if(DEBUG) assert(resbiasX[c][k+1] > resbiasX[c][k]);
		double Uinv = 1.0/(resbiasX[c][k+1] - resbiasX[c][k]);
		for(int i = Start[c][k]; i < End[c][k]; i++){// HERE : multithread
		  Cinterval *perr = &errors[c][i];
		  double x = perr->rawx;
		  double y = perr->y;
		  double err = x - y;
		  double var = A[c] + B[c] * y;
		  if(QUADRATIC_VARIANCE)
		    var += R[c]*y*y;
		  if(RES_VARIANCE)
		    var += E[c] * perr->resvar;
		  if(DEBUG>=2) assert(isfinite(var) && var > 0.0);
		  double Uik = (resbiasX[c][k+1] - x) * Uinv;
		  double Ivar = 1.0/var;
		  double UIvar = Uik * Ivar;

		  BB[k] += Uik * UIvar;
		  CC[k] += (1.0 - Uik) * UIvar;
		  DD[k] += err * UIvar;

		  meanbias[c][k] += err * Ivar;
		  Wsum += Ivar;
		}	      
		meanbias[c][k] /= Wsum;
	      }
	    }
	  
	    /* solve tri-diagonoal linear equations */
	    double nC[RESBINS+1], nD[RESBINS+1];

	    if(DEBUG) assert(BB[0] > 0.0);
	    double Binv = 1.0/BB[0];
	    nC[0] = CC[0] * Binv;
	    nD[0] = DD[0] * Binv;
	    for(int k = 1; k <= ResBins[c]; k++){
	      double Bexp = BB[k] - AA[k] * nC[k-1];
	      if(Bexp == 0.0){
		printf("-resbias equations are singular (too few molecules ?) : please rerun without -resbias (k=%d,BB[k]=%0.8e,AA[k]=%0.8e,nC[k-1]=%0.8e)\n", k,BB[k],AA[k],nC[k-1]);
		fflush(stdout);
		if(DEBUG) assert(Bexp != 0.0);
		exit(1);
	      }
	      Binv = 1.0/Bexp;
	      nC[k] = CC[k] * Binv;
	      nD[k] = (DD[k] - AA[k] * nD[k-1]) * Binv;
	    }
	  
	    /* constrain solution so highest few resbias[c][K..ResBins[c]] are <= 0.0 where K is lowest value with meanbias[c][K] >= 0 */
	    int K = 0;
	    for(;K <= ResBins[c]; K++)
	      if(meanbias[c][K] >= 0.0)
		break;
	    for(int k = K; k <= ResBins[c]; k++)
	      resbias[c][k] = 0.0;

	    if(K > ResBins[c])
	      resbias[c][K = ResBins[c]] = min(0.0,nD[ResBins[c]]);
	    for(int k = K; --k >= 0;){
	      resbias[c][k] = nD[k] - nC[k] * resbias[c][k+1];
	      if(resbias[c][k+1] == 0.0)
		resbias[c][k] = min(0.0,resbias[c][k]);
	    }

	    if(VERB>=2){
	      for(int k = ResBins[c]; k >= 0; k--)
		printf("c=%d,k=%d:resbias[c][k]=%0.8f at resbiasX[c][k]=%0.8f\n",c,k,resbias[c][k], resbiasX[c][k]);
	      fflush(stdout);
	    }
	  }

#else // SIMPLE_BIAS <= 1
	  
	  /* now estimate resbias[c][0..ResBins[c]] using piece-wise contant fit */
	  resbias[c][ResBins[c]] = 0.0;/* assume bias is 0 at or above size cutoff : will be handled (approximately) by -bpp adjustment if not correct */

	  for(int Bin = ResBins[c]; --Bin >= 0; ){
	    int BinStart = Start[c][Bin];
	    int BinEnd = End[c][Bin];
	    if(DEBUG) assert(0 <= BinStart && BinStart < BinEnd && BinEnd <= maxindex[c]+1);


	    /* resbias[c][Bin] = AVG((rawx - y)/var(y))/AVG(1/var(y)) */

	    double biassum = 0.0, Ivarsum = 0.0;

	    biascompute(biassum, Ivarsum, errors[c], Bin, BinStart, BinEnd, Tarray1, Tarray2, numthreads, A[c], B[c], R[c], E[c]);

	    if(VERB>=2){
	      printf("c=%d,Bin=%d:biassum=%0.18f,Ivarsum=%0.18f,BinStart=%d,BinEnd=%d,maxindex[c]=%d,numerrors[c]=%d\n",
		     c,Bin,biassum,Ivarsum,BinStart,BinEnd,maxindex[c],numerrors[c]);
	      fflush(stdout);
	    }

	    meanbias[c][Bin] = biassum / Ivarsum;
	    if(Bin+1 == ResBins[c])
	      resbias[c][Bin+1] = 0.0;
	    resbias[c][Bin] = meanbias[c][Bin];
	  }
#endif // SIMPLE_BIAS <= 1

	  if(VERB>=2 && nummaps > 0){
	    printf("giter=%d,c=%d:While estimating bias parameters(3):Gmap[0]->numsite[c]=M=%d,Gmap[0]->site[c][M]=%0.3f,Gmap[0]->rawsite[c][1]=%0.3f,Gmap[0]->rawsite[c][M]=%0.3f\n",
		   giter,c,Gmap[0]->numsite[c],Gmap[0]->site[c][Gmap[0]->numsite[c]],Gmap[0]->rawsite[c][1],Gmap[0]->rawsite[c][Gmap[0]->numsite[c]]);
	    fflush(stdout);
	  }
	  
	  if(VERB>=3){
	    printf("\nRaw bias values:\n");
	    for(int Bin = ResBins[c];--Bin >= 0; )
	      printf("c=%d,Bin=%d:size=%0.4f kb, bias=%7.4f (mean bias between %0.3f .. %0.3f : %0.4f, samples=%d)\n",
		     c, Bin, resbiasX[c][Bin], resbias[c][Bin], resbiasX[c][Bin],resbiasX[c][Bin+1], meanbias[c][Bin], End[c][Bin] - Start[c][Bin]);
	  }
	  if(DEBUG>=2) assert(resbiasX[c][ResBins[0]] <= maxresbias);

	  /* suppress +ve bias values : let these be handled by bpp estimation */

	  // Also force it to be monotonic :
	  resbias[c][0] = min(0.0,resbias[c][0]);
	  for(int Bin = 1; Bin <= ResBins[c]; Bin++){
	    resbias[c][Bin] = min(0.0,resbias[c][Bin]);
	    resbias[c][Bin] = max(resbias[c][Bin-1],resbias[c][Bin]);
	  }
	
	  /* create a table to quickly determine the Bin of any interval size up to maxresbias by the size expressed in bases : should match code in BiasCorrect() in Cmap.cpp*/
	  if(DEBUG) assert(resbiasX[c][ResBins[c]] == maxresbias);

	  int Ibot = (int)floor(resbiasX[c][0]*1000.0 + 0.5), Itop = 0;
	  for(int i = 0; i < Ibot; i++)
	    SizeToBin[c][i] = 0;
	  for(int Bin = 0; Bin < ResBins[c]; Bin++, Ibot = Itop){
	    Itop = (int)floor(resbiasX[c][Bin+1]*1000.0 + 0.5);
	    if(DEBUG/* HERE >=2 */ && Bin && !(Ibot >= 0)){
	      for(int B = 0; B <= ResBins[c]; B++)
		printf("c=%d,Bin=%d:resbiasX[c][Bin]=%0.4f,resbias[0][Bin]=%0.4f\n",c,B,resbiasX[c][B],resbias[c][B]);
	      printf("Bin=%d,Ibot=%d\n",Bin,Ibot);
	      fflush(stdout);
	      assert(Ibot >= 0);
	    }
	    if(DEBUG/* HERE >=2 */) assert(Ibot <= Itop && Itop <= Imaxresbias);
	    for(int i = Ibot; i < Itop; i++)
	      SizeToBin[c][i] = Bin;
	    //	  SizeToBin[c][Itop] = Bin + 1;
	  }
	  if(DEBUG) assert(Itop >= Imaxresbias);
	  for(int i = Ibot /* Itop+1 */; i <= Imaxresbias; i++)
	    SizeToBin[c][i] = ResBins[c];
	  
	  if(VERB>=2 && nummaps > 0){
	    printf("giter=%d,c=%d:While estimating bias parameters(6):Gmap[0]->numsite[c]=M=%d,Gmap[0]->site[c][M]=%0.3f,Gmap[0]->rawsite[c][1]=%0.3f,Gmap[0]->rawsite[c][M]=%0.3f\n",
		   giter,c,Gmap[0]->numsite[c],Gmap[0]->site[c][Gmap[0]->numsite[c]],Gmap[0]->rawsite[c][1],Gmap[0]->rawsite[c][Gmap[0]->numsite[c]]);
	    fflush(stdout);
	  }
	}// c = 0 .. colors-1

	/* compute Offset[c][-1..ResBins[c]], Slope[c][-1..ResBins[c]] : should match code in BiasCorrect() in Cmap.cpp */
	double OffsetMem[MAXCOLOR][RESBINS+2], SlopeMem[MAXCOLOR][RESBINS+2];
	double *Offset[MAXCOLOR], *Slope[MAXCOLOR];
	for(int c = 0; c < colors; c++){
	  Offset[c] = &OffsetMem[c][1];
	  Slope[c] = &SlopeMem[c][1];
	}
	BiasCorrectSlope(Offset,Slope);

	if(VERB/* HERE >=2 */){
	  for(int c = 0; c < colors; c++){
	    for(int Bin = ResBins[c]; --Bin >= 0;){
	      if(resbiasX[c][Bin] < 1.6)
		printf("c=%d,Bin=%d:size=%6.3f kb, bias=%0.18f (mean bias between %0.3f .. %0.3f : %0.18f, samples=%d, Offset=%0.4f,Slope=%0.6f)\n",
		       c,Bin, resbiasX[c][Bin], resbias[c][Bin], resbiasX[c][Bin],resbiasX[c][Bin+1], meanbias[c][Bin], End[c][Bin] - Start[c][Bin], Offset[c][Bin], Slope[c][Bin]);
	    }
	  }
	  fflush(stdout);
	}

	/* estimate current LL and LL with new bias estimate */
	double Lsum1[MAXCOLOR+1] = {0.0,0.0,0.0}; /* sum of log(A[c]+B[c]y+R[c]yy) excluding cases with perr->end (NOTE does not depend on x or bias parameters, except for Lsum1[2]) */
	double errsum1[MAXCOLOR+1] = {0.0,0.0,0.0};/* sum of (x - y)^2/(A[c]+B[c]y+R[c]yy) for current x (previous bias parameters) */
	double Lsum3[MAXCOLOR+1] = {0.0,0.0,0.0}, errsum3[MAXCOLOR+1] = {0.0,0.0,0.0};/* for updated x using updated bias parameters (NOTE Lsum3[0,1] should match Lsum1[0,1]) */

	double logXtheta[MAXCOLOR+1];
	for(int c = 0; c < colors;c++)
	  logXtheta[c] = 2.0*log(Xtheta[c]) - log(2.0*M_PI);/* to match terms LLsd in logLR() */
	logXtheta[colors] = 2.0*log(Xtheta12) - log(2.0*M_PI);

	for(int c = 0; c < colors; c++){
	  if(VERB>=2 && WITH_RESBIAS && giter== -1){
	    numthreads = 1;
	    qsort(errors[c],numerrors[c],sizeof(Cinterval),(intcmp*)ErrorIdInc);

	    printf("giter=%d,c=%d: Original Resbias:\n",giter,c);
	    for(int Bin = 0; Bin <= ResBins[c]; Bin++)
	      printf("Bin=%d:resbiasX[c][Bin]=%0.3f,resbias[c][Bin]=%0.6f\n",Bin,resbiasX[c][Bin],resbias[c][Bin]);
	    if(nummaps > 43){
	      Cmap *Xmap = Gmap[43];
	      FLOAT *X = Xmap->site[c];
	      FLOAT *rawX = Xmap->rawsite[c];
	      int M = Xmap->numsite[c];
	      for(int I = 0; I <= M+1; I++)
		printf("I=%d: X[I]=%0.6f, rawX[I]=%0.6f\n",
		       I, X[I], rawX[I]);
	    }
	    fflush(stdout);
	  }

	  for(int tid = 0; tid < numthreads; tid++)
	    Tarray1[tid] = Tarray2[tid] = Tarray3[tid] = 0.0;

          #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
	  {
	    int tid = 0;
            #ifdef _OPENMP
	    tid = omp_get_thread_num ();
            #endif

	    double myLsum1 = 0.0;
	    double myerrsum1 = 0.0;

            #pragma omp for schedule(static,64)
	    for(int i = 0; i < numerrors[c]; i++){
	      Cinterval *perr = &errors[c][i];
	      double y = perr->y;
	      double x = perr->x;
#if REFDEBUG >= 1
	      perr->origx = x;// for debugging : see below
#endif
	      if(ENDFIX>=3 && perr->end && y >= x){
		if(VERB>=2 && WITH_RESBIAS && giter== -1 && perr->mapid == 43){
		  printf("c=%d,i=%d:    mapid=%d,L=%d,R=%d,end=%d:x=%0.6f,y=%0.6f : skipping since y >= x\n",
			 c,i,perr->mapid,perr->L,perr->R, perr->end, x, y);
		  fflush(stdout);
		}
	      
		continue;
	      }
	      double var = A[c] + B[c]*y;
	      if(QUADRATIC_VARIANCE)
		var += R[c]*y*y;
	      if(RES_VARIANCE)
		var += E[c] * perr->resvar;
	      if(DEBUG>=2) assert(isfinite(var) && var > 0.0);
	      double err1 = y - x;
	      if(!(ENDFIX>=2 && perr->end)){
		myLsum1 += log(var);
		if(VERB>=2 || SDDEBUG)
		  myLsum1 -= logXtheta[c];
	      }
	      double Ivar = 1.0/var;
	      if(VERB>=2 && WITH_RESBIAS && giter== -1 && perr->mapid == 43){
		FLOAT *X = Gmap[perr->mapid]->site[c];
		int M = Gmap[perr->mapid]->numsite[c];
		if(!(ENDFIX>=2 && perr->end))
		  printf("c=%d,i=%d,tid=%d:    mapid=%d,L=%d,R=%d,end=%d:Lsum=%0.6f,errterm=%0.6f,LLsd=%0.6f,x=%0.6f,y=%0.6f\n",
			 c,i,tid,perr->mapid,perr->L,perr->R,perr->end,-0.5*(log(var) - logXtheta[c]), -0.5*err1*err1*Ivar,
			 -0.5*(log(var)-logXtheta[c] + err1*err1*Ivar),x,y);
		else
		  printf("c=%d,i=%d,tid=%d:    mapid=%d,L=%d,R=%d,end=%d:LLsd=%0.6f,x=%0.6f,y=%0.6f, X[L]=%0.6f,X[R]=%0.6f,X[%d]=%0.6f\n",
			 c,i,tid,perr->mapid,perr->L,perr->R,perr->end, -0.5*err1*err1*Ivar, x,y, X[perr->L], X[perr->R], M+1, X[M+1]);
		fflush(stdout);
	      }
	      myerrsum1 += err1*err1*Ivar;
	    
	      if(DEBUG>=2 && !isfinite(myerrsum1)){
		printf("c=%d,i=%d/%d:x=%0.4f,y=%0.4f,var=%0.6f,err1=%0.4f,myerrsum1=%0.8f\n",
		       c,i,numerrors[c], x, y, var, err1, myerrsum1);
		fflush(stdout);
		assert(isfinite(myerrsum1));
	      }
	    }

	    Tarray1[tid] = myLsum1;
	    Tarray3[tid] = myerrsum1;
	  }
	  for(int tid = 0; tid < numthreads; tid++){
	    Lsum1[c] += Tarray1[tid];
	    errsum1[c] += Tarray3[tid];
	  }
	  if(DEBUG>=2) assert(isfinite(Lsum1[c]));
	  if(DEBUG>=2) assert(isfinite(errsum1[c]));

	  if(VERB>=2){
	    printf("c=%d,BiasCorrect1:maxresbias=%0.10f,ResBins[c]=%d:\n",c,maxresbias,ResBins[c]);
	    for(int Bin = 0; Bin <= ResBins[c]; Bin++)
	      printf("Bin=%d:resbiasX[c][Bin]=%0.10f,resbias[c][Bin]=%0.10f\n",Bin,resbiasX[c][Bin],resbias[c][Bin]);
	    fflush(stdout);
	  }

	}// c = 0 .. colors-1
	
	/* compute Lsum1[2] and errsum1[2] based on SMA() terms and current bias parameters */
	Lsum1[2] = errsum1[2] = 0.0;
	for(int c = colors; c <= colors; c++){
	  if(VERB>=2 && giter==0){
	    printf("  giter=%d,c=%d:Before bias correction:MA[2]=%0.8f,MA_mean=%0.8f,SF[2]*SF[2]*2 = %0.8f,FM2var = %0.8f, fabs(SD[2])*SD[2]*2= %0.8f, SM2var= %0.8f, 0.5*logXtheta[2]=%0.8f, RefPen12 - 0.5*log(2)=%0.8f, log(Xtheta12) - 0.5*log(2*PI)=%0.8f\n",
		   giter,c,MA[2],MA_mean,SF[2]*SF[2]*2.0, FM2var, fabs(SD[2])*SD[2]*2.0, SM2var, 0.5*logXtheta[2], RefPen12 - 0.5*log(2.0), log(Xtheta12) - 0.5*log(2.0*M_PI));
	    fflush(stdout);
	  }

	  Cinterval *perr = errors[c];
	  for(int i = 0; i < numerrors[c]; i++){
	    double y = perr[i].y + MA[c] * perr[i].end;/* WAS perr->y - MA[c] */
	    double x = perr[i].x;
	    double len = perr[i].len;
	    double var = A[c]+B[c]*len;
	    double err = x - y;
	    Lsum1[2] += log(var);
	    if(VERB>=2 || SDDEBUG)
	      Lsum1[2] -= logXtheta[c];
	    errsum1[2] += err*err/var;
	    if(VERB>=2 && SDDEBUG && giter==0)
	      printf("  i=%d/%d:mapid=%d,end=%d:x=%0.4f,y=%0.4f,len=%0.4f(%0.4f):SMA(x,y,len,or=%d)=%0.6f,LR2=%0.6f(cum=%0.6f)\n",i,numerrors[2],perr[i].mapid,perr[i].end, x, perr[i].y,perr[i].len,len,
		     (perr[i].end < 0) ? 1 : 0, SMA(x, perr[i].y,len,(perr[i].end < 0) ? 1 : 0), -0.5*(log(var) - logXtheta[c] + err*err/var), -0.5*(Lsum1[2] + errsum1[2]));
	  }
	} // c == colors

	/* Apply full bias removal to all sites */
	/* Same as NEW CODE as in BiasCorrect in Cmap.cpp */

	int nthreads = 1;
#ifdef _OPENMP
	nthreads = omp_get_max_threads();
	nthreads = min(nthreads,MaxThreads);
	nthreads = max(1,min(nthreads, nummaps/2048));
#endif  

	/* NOTE : no need to handle -subset maps or split maps since they will not be re-aligned */
	BiasCorrect2(Gmap, 0, nummaps, nthreads, SizeToBin, Imaxresbias, Offset, Slope);

#ifdef _OPENMP
	numthreads = omp_get_max_threads();
	numthreads = min(numthreads,MaxThreads);
	if(numthreads > max(numerrors[0],numerrors[1])/64)
	  numthreads = max(1,max(numerrors[0],numerrors[1])/64);
#endif

	if(VERB>=2 && WITH_RESBIAS && giter== -1)
	  numthreads = 1;

	double LL1[MAXCOLOR+1] = {0.0,0.0,0.0}, LL3[MAXCOLOR+1] = {0.0,0.0,0.0};

	for(int c = 0; c < colors; c++){
	  /* now recompute x values in errors[] based on new site[] values and new LL value */
	  Lsum3[c] = errsum3[c] = 0.0;
	  for(int tid = 0; tid < numthreads; tid++)
	    Tarray2[tid] = Tarray3[tid] = 0.0;

	  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
	  {
	    int tid = 0;
            #ifdef _OPENMP
	    tid = omp_get_thread_num ();
            #endif

	    int lastid  = -1;
	    double Errsum1 = 0.0, Errcum1 = 0.0;/* for debugging */
	    double LogSum3 = 0.0, Errsum3 = 0.0, LogCum3 = 0.0, Errcum3 = 0.0;/* for debugging */

	    double myLsum3 = 0.0, myerrsum3 = 0.0;

	    #pragma omp for schedule(static,64)
  	    for(int i = 0; i < numerrors[c]; i++){
	      Cinterval *perr = &errors[c][i];
	      FLOAT *X = Gmap[perr->mapid]->site[c], *rawX;
	      FLOAT newx = X[perr->R] - X[perr->L], oldx;
	      if(0 && DEBUG && (rawX = Gmap[perr->mapid]->rawsite[c], oldx  = rawX[perr->R] - rawX[perr->L]) <= maxresbias && 
		 resbias[c][0] < 0.0 && !(newx >= oldx + resbias[c][0] - 1e-6 && newx <= oldx - resbias[c][0] + 1e-6)){
	        #pragma omp critical
	        {
		  printf("c=%d,error[%d]:mapid=%d:L=%d,R=%d,M=%d:newx=%0.4f,oldx=%0.4f,resbias[c][0]=%0.4f:X[L-1..R+1]=%0.4f,%0.4f,%0.4f,%0.4f rawX[L-1..R+1]=%0.4f,%0.4f,%0.4f,%0.4f\n",
			 c,i,perr->mapid,perr->L,perr->R,Gmap[perr->mapid]->numsite[c],newx,oldx,resbias[c][0],X[perr->L-1],X[perr->L],X[perr->R],X[perr->R+1],rawX[perr->L-1],rawX[perr->L],rawX[perr->R],rawX[perr->R+1]);
		  assert(newx >= oldx + resbias[c][0] - 1e-6 && newx <= oldx - resbias[c][0] + 1e-6);
		}
	      }
	      double x = perr->x = newx;
	      if(DEBUG>=2 && !(perr->x > 0.0)){
		#pragma omp critical
		{
		  printf("c=%d,i=%d/%d:mapid=%d,L=%d,R=%d,M=%d:X[L]=%0.4f,X[R]=%0.4f,newx=%0.4f\n",
			 c,i,numerrors[c],perr->mapid,perr->L,perr->R,Gmap[perr->mapid]->numsite[c],X[perr->L],X[perr->R],newx);
		  fflush(stdout);
		  assert(perr->x > 0.0);
		}
	      }
	      double y = perr->y;
	      if(ENDFIX>=3 && perr->end && y >= x){
		if(VERB>=2 && WITH_RESBIAS && giter== -1 && perr->mapid == 43){
		  printf("    c=%d,mapid=%d,L=%d,R=%d,end=%d:x=%0.6f,y=%0.6f : skipping since y >= x\n",
			 c,perr->mapid,perr->L,perr->R, perr->end, x, y);
		  fflush(stdout);
		}
		continue;
	      }
	      double var = A[c] + B[c]*y;
	      if(QUADRATIC_VARIANCE)
		var += R[c]*y*y;
	      if(RES_VARIANCE)
		var += E[c] * perr->resvar;
	      if(DEBUG>=2) assert(var > 0.0);
	      double Ivar = 1.0/var;
	      double err = y - x;

	      if(VERB>=2 && WITH_RESBIAS && giter== -1){/* display per molecule statistics ( NOTE: must call qsort() before this parallel loop and reduced threads to 1) */
		int mapid = perr->mapid;
		if(mapid != lastid){
		  Errcum1 += Errsum1;
		  LogCum3 += LogSum3;
		  Errcum3 += Errsum3;
		  if(lastid >= 0){/* display values for mol lastid */
		    printf("  c=%d,mapid=%d(id=%lld):Lsum=%0.6f,errsum=%0.6f->%0.6f,LLsd=%0.6f -> %0.6f(delta=%0.6f,cum=%0.6f -> %0.6f)\n",
			   c,lastid,Gmap[lastid]->id, -0.5*LogSum3, -0.5*Errsum3, -0.5*Errsum3, -0.5*(LogSum3+Errsum1), -0.5*(LogSum3+Errsum3),-0.5*(Errsum3-Errsum1),
			   -0.5*(LogCum3 + Errcum1), -0.5*(LogCum3 + Errcum3));
		    fflush(stdout);
		  }
		  lastid = mapid;
		  LogSum3 = Errsum1 = Errsum3 = 0.0;
		}
		if(!(ENDFIX>=2 && perr->end)){
		  LogSum3 += log(var)- logXtheta[c];
		  if(DEBUG){
		    rawX = Gmap[perr->mapid]->rawsite[c];
		    oldx  = rawX[perr->R] - rawX[perr->L];
		    assert(fabs(oldx - perr->rawx) < 1e-6);
		  }
		}
		Errsum3 += err*err*Ivar;
#if REFDEBUG>=1
		double err1 = y - perr->origx;
		Errsum1 += err1*err1*Ivar;
		if(mapid == -1){
		  int M = Gmap[perr->mapid]->numsite[c];
		  if(!(ENDFIX>=2 && perr->end))
		    printf("    c=%d,mapid=%d,L=%d,R=%d,end=%d:Lsum=%0.6f,errterm=%0.6f->%0.6f,LLsd=%0.6f->%0.6f(delta=%0.6f,cum=%0.6f -> %0.6f),x=%0.6f,y=%0.6f\n",
			   c,mapid,perr->L,perr->R,perr->end,-0.5*(log(var) - logXtheta[c]), -0.5*err1*err1*Ivar, -0.5*err*err*Ivar,
			   -0.5*(log(var)-logXtheta[c] + err1*err1*Ivar), -0.5*(log(var)-logXtheta[c] + err*err*Ivar), -0.5*(err*err*Ivar - err1*err1*Ivar),
			   -0.5*(LogSum3+Errsum1),-0.5*(LogSum3+Errsum3), x, y);
		  else
		    printf("    c=%d,mapid=%d,L=%d,R=%d,end=%d:LLsd=%0.6f->%0.6f(delta=%0.6f,cum=%0.6f -> %0.6f),x=%0.6f,y=%0.6f, X[L]=%0.6f,X[R]=%0.6f,X[%d]=%0.6f\n",
			   c,mapid,perr->L,perr->R,perr->end, -0.5*err1*err1*Ivar, -0.5*err*err*Ivar, -0.5*(err*err*Ivar - err1*err1*Ivar),
			   -0.5*(LogSum3+Errsum1),-0.5*(LogSum3+Errsum3), x, y, X[perr->L], X[perr->R], M+1, X[M+1]);
		  fflush(stdout);
		}
#endif
	      }
	      if(DEBUG>=2 && !(isfinite(err))){
	        #pragma omp critical
	        {
		  int M = Gmap[perr->mapid]->numsite[c];
		  rawX = Gmap[perr->mapid]->rawsite[c];
		  printf("c=%d,R=%d,L=%d,M=%d:X[R]=%0.3f,X[L]=%0.3f,x=%0.3f,y=%0.3f,rawX[R]=%0.3f,rawX[L]=%0.3f,err=%0.3f\n",
			 c,perr->R,perr->L,M,X[perr->R],X[perr->L],x,y,rawX[perr->R],rawX[perr->L],err);
		  printf("  mapid=%d,id=%lld\n",perr->mapid,Gmap[perr->mapid]->id);
		  for(int J = 1; J <= M+1; J++)
		    printf("J=%d:rawX[J]=%0.3f,X[J]=%0.3f\n",J,rawX[J],X[J]);
		  fflush(stdout);
		  assert(isfinite(err));
		}
	      }
	      myerrsum3 += err*err*Ivar;
	      if(!(ENDFIX>=2 && perr->end)){
		myLsum3 += log(var);
		if(VERB>=2 || SDDEBUG)
		  myLsum3 -= logXtheta[c];
	      }
	      if(DEBUG>=2) assert(isfinite(myerrsum3));
	      if(DEBUG>=2) assert(isfinite(myLsum3));
	    }
	    if(VERB>=2 && WITH_RESBIAS && giter== -1 && lastid >= 0){/* display values for lastid */
	      Errcum1 += Errsum1;
	      LogCum3 += LogSum3;
	      Errcum3 += Errsum3;
	      printf("  c=%d,mapid=%d(id=%lld):Lsum=%0.6f,errsum=%0.6f->%0.6f,LLsd=%0.6f -> %0.6f(delta=%0.6f,cum=%0.6f -> %0.6f)\n",
		     c,lastid,Gmap[lastid]->id, -0.5*LogSum3, -0.5*Errsum1, -0.5*Errsum3, -0.5*(LogSum3+Errsum1), -0.5*(LogSum3+Errsum3),-0.5*(Errsum3-Errsum1),
		     -0.5*(LogCum3 + Errcum1), -0.5*(LogCum3 + Errcum3));
	      fflush(stdout);
	    }

	    Tarray2[tid] = myLsum3;
	    Tarray3[tid] = myerrsum3;
	  }
	  for(int tid = 0; tid < numthreads; tid++){
	    Lsum3[c] += Tarray2[tid];
	    errsum3[c] += Tarray3[tid];
	  }

	  LL1[c] = -(Lsum1[c]+errsum1[c])*0.5;
	  LL3[c] = -(Lsum3[c]+errsum3[c])*0.5;
	  if(DEBUG) assert(isfinite(LL1[c]));
	  if(DEBUG) assert(isfinite(LL3[c]));
	} // c = 0 .. colors -1

	/* compute Lsum3[2] and errsum3[2] based on SMA() terms and current bias parameters */
	Lsum3[2] = errsum3[2] = 0.0;
	for(int c = colors; c <= colors; c++){
	  if(VERB>=2 && giter==0){
	    printf("  giter=%d,c=%d:After bias correction:MA[2]=%0.8f,MA_mean=%0.8f,SF[2]*SF[2]*2 = %0.8f,FM2var = %0.8f, fabs(SD[2])*SD[2]*2= %0.8f, SM2var= %0.8f, 0.5*logXtheta[2]=%0.8f, RefPen12 - 0.5*log(2)=%0.8f, log(Xtheta12) - 0.5*log(2*PI)=%0.8f\n",
		   giter,c, MA[2],MA_mean,SF[2]*SF[2]*2.0, FM2var, fabs(SD[2])*SD[2]*2.0, SM2var, 0.5*logXtheta[2], RefPen12 - 0.5*log(2.0), log(Xtheta12) - 0.5*log(2.0*M_PI));
	    fflush(stdout);
	  }

	  Cinterval *perr = errors[c];
	  for(int i = 0; i < numerrors[c]; i++){
#if REFDEBUG>=1
	    perr[i].origlen = perr[i].len;
	    perr[i].origx = perr[i].x;
#endif
	    /* first update perr[i].x and perr[i].len */
	    int L = perr[i].L;
	    int R = perr[i].R;
	    Cmap *Xmap = Gmap[perr[i].mapid];
	    FLOAT **XX = Xmap->site;
	    int *MM = Xmap->numsite;
	    double x0 = (perr[i].end < 0) ? (XX[0][MM[0]+1] - XX[0][MM[0]+1-L]) : XX[0][L];
	    double x1 = (perr[i].end < 0) ? (XX[1][MM[1]+1] - XX[1][MM[1]+1-R]) : XX[1][R];
	    perr[i].x = x1 - x0;
	    if(COLOR_SHIFT==1)
	      perr[i].len = max(x0,x1);

	    double y = perr[i].y + MA[c] * perr[i].end;
	    double x = perr[i].x;
	    double len = perr[i].len;
	    double var = A[c]+B[c]*len;
	    double err = x - y;
	    Lsum3[2] += log(var);
	    if(VERB>=2 || SDDEBUG)
	      Lsum3[2] -= logXtheta[c];
	    errsum3[2] += err*err/var;
	    if(VERB>=2 && SDDEBUG && giter==0)
	      printf("  i=%d/%d:mapid=%d,end=%d:x=%0.4f(L=%d,R=%d,M=%d,%d,x0=%0.4f,x1=%0.4f),y=%0.4f,len=%0.4f(%0.4f):SMA(x,y,len,or=%d)=%0.6f,LR2=%0.6f(cum=%0.6f)\n",
		     i,numerrors[2],perr[i].mapid,perr[i].end, x, L,R,MM[0],MM[1],x0,x1,perr[i].y,perr[i].len,len,
		     (perr[i].end < 0) ? 1 : 0, SMA(x, perr[i].y,len,(perr[i].end < 0) ? 1 : 0), -0.5*(log(var) - logXtheta[c] + err*err/var), -0.5*(Lsum3[2] + errsum3[2]));
	  }

	  LL1[c] = -(Lsum1[c] + errsum1[c])*0.5;
	  LL3[c] = -(Lsum3[c] + errsum3[c])*0.5;
	  if(DEBUG) assert(isfinite(LL1[c]));
	  if(DEBUG) assert(isfinite(LL3[c]));
	} // c == colors

	if(VERB){
	  printf("Previous LL=%0.6f(%0.6f,%0.6f,%0.6f),errsum1=%0.6f, New Bias Correction LL=%0.6f(%0.6f,%0.6f,%0.6f),Lsum3=%0.6f,errsum3=%0.6f:wall time=%0.6f\n",
		 LL1[0]+LL1[1]+LL1[2], LL1[0],LL1[1],LL1[2], errsum1[0]+errsum1[1]+errsum1[2], LL3[0]+LL3[1]+LL3[2], LL3[0],LL3[1],LL3[2], Lsum3[0]+Lsum3[1]+Lsum3[2], errsum3[0]+errsum3[1]+errsum3[2], wtime());
	  fflush(stdout);
	}

	if(LL3[0]+LL3[1]+LL3[2] < LL1[0]+LL1[1]+LL1[2]){/* backtrack to original bias values */
	  if(VERB){
	    printf("Undoing resbias changes\n");
	    fflush(stdout);
  	  }

	  for(int c = 0; c < colors; c++){
	    ResBins[c] = startResBins[c];
	    for(int Bin = 0; Bin <= ResBins[c]; Bin++){
	      resbiasX[c][Bin] = startresbiasX[c][Bin];
	      resbias[c][Bin] = startresbias[c][Bin];
	    }
	  }
	  /* NOTE : no need to handle -subset maps or split maps since they will not be re-aligned */
	  BiasCorrect(Gmap, 0, nummaps, 0);

	  for(int c = 0; c < colors; c++){
	    /* recompute errors[i].x */
	    if(REFDEBUG){
	      Lsum3[c] = errsum3[c] = 0.0;
	      for(int tid = 0; tid < numthreads; tid++)
		Tarray2[tid] = Tarray3[tid] = 0.0;
	    }

	    if(VERB>=2 && giter==10)
	      numthreads = 1;

  	    #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
	    {
	      int tid = 0;
              #ifdef _OPENMP
	      tid = omp_get_thread_num ();
              #endif

	      double myLsum3 = 0.0, myerrsum3 = 0.0;

              #pragma omp for schedule(static,64)
	      for(int i = 0; i < numerrors[c]; i++){
		Cinterval *perr = &errors[c][i];
		FLOAT *X = Gmap[perr->mapid]->site[c];
		double x = perr->x = X[perr->R] - X[perr->L];
		if(DEBUG>=2 && !(perr->x > 0.0)){
		  #pragma omp critical
		  {
		    int M = Gmap[perr->mapid]->numsite[c];
		    printf("c=%d,i=%d/%d:mapid=%d,L=%d,R=%d,M=%d:X[L]=%0.4f,X[R]=%0.4f,X[M+1]=%0.4f\n",
			   c,i,numerrors[c],perr->mapid,perr->L,perr->R,M,X[perr->L],X[perr->R],X[M+1]);
		    fflush(stdout);
		    assert(perr->x > 0.0);
		  }
		}

		if(REFDEBUG){
		  double y = perr->y;
		  if(ENDFIX>=3 && perr->end && y >= x){
		    if(VERB>=2 && WITH_RESBIAS && giter== -1 && perr->mapid == 43){
#if REFDEBUG>=1
		      printf("c=%d,i=%d,tid=%d:    mapid=%d,L=%d,R=%d:x=%0.6f,y=%0.6f,end=%d : skipping since y >= x (origx=%0.6f)\n",
			     c,i,tid,perr->mapid,perr->L,perr->R, x, y, perr->end, perr->origx);
#else
		      printf("c=%d,i=%d,tid=%d:    mapid=%d,L=%d,R=%d:x=%0.6f,y=%0.6f,end=%d : skipping since y >= x\n",
			     c,i,tid,perr->mapid,perr->L,perr->R, x, y, perr->end);
#endif
		      fflush(stdout);
		    }
		    continue;
		  }
		  double var = A[c] + B[c]*y;
		  if(QUADRATIC_VARIANCE)
		    var += R[c]*y*y;
		  if(RES_VARIANCE)
		    var += E[c] * perr->resvar;
		  if(DEBUG>=2) assert(var > 0.0);
		  double Ivar = 1.0/var;
		  double err = y - x;
		  if(VERB>=2 && WITH_RESBIAS && giter== -1 && perr->mapid == 43){
		    int M = Gmap[perr->mapid]->numsite[c];
#if REFDEBUG>=1
		    if(!(ENDFIX>=2 && perr->end))
		      printf("c=%d,i=%d,tid=%d:    mapid=%d,L=%d,R=%d,end=%d:Lsum=%0.6f,errterm=%0.6f,LLsd=%0.6f,x=%0.6f,y=%0.6f (origx=%0.6f)\n",
			     c,i,tid,perr->mapid,perr->L,perr->R,perr->end,-0.5*(log(var) - logXtheta[c]), -0.5*err*err*Ivar,
			     -0.5*(log(var)-logXtheta[c] + err*err*Ivar),x,y, perr->origx);
		    else
		      printf("c=%d,i=%d,tid=%d:    mapid=%d,L=%d,R=%d,end=%d:LLsd=%0.6f,x=%0.6f,y=%0.6f (origx=%0.6f), X[L]=%0.6f,X[R]=%0.6f,X[%d]=%0.6f\n",
			     c,i,tid, perr->mapid,perr->L,perr->R,perr->end, -0.5*err*err*Ivar, x,y, perr->origx, X[perr->L], X[perr->R], M+1, X[M+1]);
#else
		    if(!(ENDFIX>=2 && perr->end))
		      printf("c=%d,i=%d,tid=%d:    mapid=%d,L=%d,R=%d,end=%d:Lsum=%0.6f,errterm=%0.6f,LLsd=%0.6f,x=%0.6f,y=%0.6f\n",
			     c,i,tid,perr->mapid,perr->L,perr->R,perr->end,-0.5*(log(var) - logXtheta[c]), -0.5*err*err*Ivar,
			     -0.5*(log(var)-logXtheta[c] + err*err*Ivar),x,y);
		    else
		      printf("c=%d,i=%d,tid=%d:    mapid=%d,L=%d,R=%d,end=%d:LLsd=%0.6f,x=%0.6f,y=%0.6f, X[L]=%0.6f,X[R]=%0.6f,X[%d]=%0.6f\n",
			     c,i,tid, perr->mapid,perr->L,perr->R,perr->end, -0.5*err*err*Ivar, x,y, X[perr->L], X[perr->R], M+1, X[M+1]);
#endif
  
		    fflush(stdout);
		  }
		  myerrsum3 += err*err*Ivar;
		  if(!(ENDFIX>=2 && perr->end)){
		    myLsum3 += log(var);
		    if(VERB>=2 || SDDEBUG)
		      myLsum3 -= logXtheta[c];
		  }
		}// REFDEBUG
	      } // for(i=0;i<numerrors;i++)
	      if(REFDEBUG){
		Tarray2[tid] = myLsum3;
		Tarray3[tid] = myerrsum3;
	      }
	    } // parallel

	    if(REFDEBUG){
	      Lsum3[c] = errsum3[c] = 0.0;
	      for(int tid = 0; tid < numthreads; tid++){
		Lsum3[c] += Tarray2[tid];
		errsum3[c] += Tarray3[tid];
	      }
	      LL3[c] = -(Lsum3[c] + errsum3[c])*0.5;
	    }
	  } // c = 0 .. colors - 1

	  /* compute Lsum3[2] and errsum3[2] based on SMA() terms and restored bias parameters */
	  Lsum3[2] = errsum3[2] = 0.0;
	  for(int c = colors; c <= colors; c++){
	    if(VERB>=2){
	      printf("  MA[2]=%0.8f,MA_mean=%0.8f,SF[2]*SF[2]*2 = %0.8f,FM2var = %0.8f, fabs(SD[2])*SD[2]*2= %0.8f, SM2var= %0.8f, 0.5*logXtheta[2]=%0.8f, RefPen12 - 0.5*log(2)=%0.8f, log(Xtheta12) - 0.5*log(2*PI)=%0.8f\n",
		     MA[2],MA_mean,SF[2]*SF[2]*2.0, FM2var, fabs(SD[2])*SD[2]*2.0, SM2var, 0.5*logXtheta[2], RefPen12 - 0.5*log(2.0), log(Xtheta12) - 0.5*log(2.0*M_PI));
	      fflush(stdout);
	    }

	    Cinterval *perr = errors[c];

	    for(int i = 0; i < numerrors[c]; i++){
	      /* first update perr[i].x and perr[i].len */
	      int L = perr[i].L;
	      int R = perr[i].R;
	      Cmap *Xmap = Gmap[perr[i].mapid];
	      FLOAT **XX = Xmap->site;
	      int *MM = Xmap->numsite;
	      double x0 = (perr[i].end < 0) ? (XX[0][MM[0]+1] - XX[0][MM[0]+1-L]) : XX[0][L];
	      double x1 = (perr[i].end < 0) ? (XX[1][MM[1]+1] - XX[1][MM[1]+1-R]) : XX[1][R];
	      perr[i].x = x1 - x0;
	      if(COLOR_SHIFT==1)
		perr[i].len = max(x0,x1);
#if REFDEBUG>=1
	      if(fabs(perr[i].len - perr[i].origlen) > 1e-6 || fabs(perr[i].x - perr[i].origx) > 1e-6){
		printf("c=%d,i=%d/%d:mapid=%d,L=%d,R=%d,end=%d,x0=%0.6f,x1=%0.6f:x=%0.6f(orig=%0.6f),len=%0.6f(orig=%0.6f)\n",
		       c,i,numerrors[c],perr[i].mapid,L,R,perr[i].end,x0,x1,perr[i].x,perr[i].origx,perr[i].len,perr[i].origlen);
		if(perr[i].end < 0)
		  printf("   MM[0]=%d,MM[1]=%d,XX[0][MM[0]+1]=%0.6f,XX[0][MM[0]+1-L]=%0.6f,XX[1][MM[1]+1]=%0.6f,XX[1][MM[1]+1-R]=%0.6f\n",
			 MM[0],MM[1],XX[0][MM[0]+1],XX[0][MM[0]+1-L],XX[1][MM[1]+1],XX[1][MM[1]+1-R]);
		else
		  printf("   MM[0]=%d,MM[1]=%d,XX[0][L]=%0.6f,XX[1][R]=%0.6f\n",MM[0],MM[1],XX[0][L],XX[1][R]);
		fflush(stdout);
		assert(fabs(perr[i].len - perr[i].origlen) <= 1e-6);
		assert(fabs(perr[i].x - perr[i].origx) <= 1e-6);
	      }
#endif

	      if(REFDEBUG){
		double y = perr[i].y + MA[c] * perr[i].end;
		double x = perr[i].x;
		double len = perr[i].len;
		double var = A[c]+B[c]*len;
		double err = x - y;
		Lsum3[2] += log(var);
		if(VERB>=2 || SDDEBUG)
		  Lsum3[2] -= logXtheta[c];
		errsum3[2] += err*err/var;
		if(VERB>=2 && SDDEBUG && giter==0)
		  printf("  i=%d/%d:mapid=%d,end=%d:x=%0.4f,y=%0.4f,len=%0.4f(%0.4f):SMA(x,y,len,or=%d)=%0.6f,LR2=%0.6f(cum=%0.6f)\n",i,numerrors[2],perr[i].mapid,perr[i].end, x, perr[i].y,perr[i].len,len,
		     (perr[i].end < 0) ? 1 : 0, SMA(x, perr[i].y,len,(perr[i].end < 0) ? 1 : 0), -0.5*(log(var) - logXtheta[c] + err*err/var), -0.5*(Lsum3[2] + errsum3[2]));
	      }
	    }

	    LL3[c] = -(Lsum3[c] + errsum3[c])*0.5;
	    if(DEBUG) assert(isfinite(LL3[c]));
	  } // c == colors

	  if(REFDEBUG){
	    double LL4 = LL3[0]+LL3[1]+LL3[2];
	    if(VERB){
	      printf("  After Undoing resbias changes: LL=%0.6f(%0.6f,%0.6f,%0.6f)(errsum=%0.6f) :original value was LL=%0.6f(%0.6f,%0.6f,%0.6f)\n",
		     LL4,LL3[0],LL3[1],LL3[2],errsum3[0]+errsum3[1]+errsum3[2],LL1[0]+LL1[1]+LL1[2],LL1[0],LL1[1],LL1[2]);
	      for(int c = 0; c < colors; c++){
		for(int Bin = 0; Bin <= ResBins[c]; Bin++)
		  printf("c=%d,Bin=%d:resbiasX[c][Bin]=%0.3f,resbias[c][Bin]=%0.6f\n",c,Bin,resbiasX[c][Bin],resbias[c][Bin]);
		if(REFDEBUG>=3 && WITH_RESBIAS && giter == -1 && nummaps > 43){
		  Cmap *Xmap = Gmap[43];
		  FLOAT *X = Xmap->site[c];
		  FLOAT *rawX = Xmap->rawsite[c];
		  int M = Xmap->numsite[c];
		  for(int I = 0; I <= M+1; I++)
		    printf("c=%d,I=%d: X[I]=%0.6f, rawX[I]=%0.6f\n",
			   c,I, X[I], rawX[I]);
		}
		fflush(stdout);
	      }
	      if(DEBUG) assert(fabs(LL4 - (LL1[0]+LL1[1]+LL1[2])) < tcnt * 1e-6);
	    }
	  }
	}

	if(REFDEBUG && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	  if(DEBUG) assert(logLRarray);
	  (void)score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	  double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	  double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	  lcnt= 0;
	  double *origLRarray = new double[nummaps];
	  double *origSDarray[MAXCOLOR];
	  for(int mid = 0; mid < nummaps; mid++)
	    origLRarray[mid] = logLRarray[mid];
	  for(int c = 0; c < colors; c++){
	    origSDarray[c] = new double[nummaps];
	    for(int mid = 0; mid < nummaps; mid++){
	      origSDarray[c][mid] = logSDarray[c][mid];
	    }
	  }

	  double origLRsum = 0.0;
	  double origSDsum[MAXCOLOR] = {0.0,0.0};

	  if(VERB>=2 && WITH_RESBIAS && giter== -1)/* sort alignments in ascending order of mapid2 */
	    qsort(alignment, numaligns, sizeof(Calign*), (intcmp*)CalignMapidInc);

	  for(size_t i = 0; i < numaligns; i++){
	    double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	    Calign *align = alignment[i];
	    if(!align || align->numpairs <= 1)
	      continue;
	    int mid = align->mapid2;
	    if(DEBUG)assert(!Gmap[mid]->origmap);
	    if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
	      if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
	      continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	    }
	    int rid = align->mapid1;
	    if(DEBUG>=2){
	      int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
	      assert(origgoodalign);
	      int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
	      assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	    }
	    FLOAT **XX = Gmap[mid]->site;
	    int *MM = Gmap[mid]->numsite;
	    FLOAT **YY = refmap[rid]->site;
	    //	    int *NN = refmap[rid]->numsite;
	    align++;
	    logLRarray[mid] = 0.0;
	    for(int c = 0; c < colors; c++){
	      FLOAT *X = XX[c];
	      double M = MM[c];
	      FLOAT *Y = YY[c];
	      int N = refmap[rid]->numsite[c];
	      int U = align[c].numpairs;

	      LR[c] = logLR(mid,rid,c,X,Y,M,N,align,align[c].sites1[0], align[c].sites2[0], align[c].sitesK1[0], align[c].sites2[U-1], align[c].Lij1,align[c].Rij1,align[c].Lij2,align[c].Rij2,1,
			    0 /* (giter==0 && c==0 && mid == 674) ? 1 : 0 */, SDDEBUG ? &LRsd[c] : 0,NULL);
	      logLRarrayC[c][mid] = LR[c];
	      logSDarray[c][mid] = LRsd[c];
	      if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
		continue;
	      LRsum += LR[c];
	      if(SDDEBUG){
		LRsdsum += LRsd[c];
		LRsdsumC[c] += LRsd[c];
	      }
	    } // c = 0 .. colors -1

	    int U = align[0].numpairs;
	    int I = align[0].sites1[U-1];
	    int K = align[0].sitesK1[U-1];
	    int J = align[0].sites2[U-1];
	    U = align[1].numpairs;
	    int P = align[1].sites1[U-1];
	    int T = align[1].sitesK1[U-1];
	    int Q = align[1].sites2[U-1];
	    double y0 = Yc(YY[0],I,K);
	    double y1 = Yc(YY[1],P,T);
	    double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	    double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	    double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);

	    LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);

	    if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	      LRsum += LR[2];
	      for(int lc = 0; lc <= colors; lc++)
		LRsumC[lc] += LR[lc];
	    }
	    logLRarray[mid] = LR[0] + LR[1] + LR[2];
	    lcnt++;
	    if(VERB>=2 && giter==0){
	      if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		origLRsum += origLRarray[mid];
		for(int lc = 0; lc < colors; lc++)
		  origSDsum[lc] += origSDarray[lc][mid];
	      }
	      if(1){
		printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)\n\t delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f\n",
		       i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
		       logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
		/*	      printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
			      logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
			      LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
			      LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
			      printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);*/
	      }
	    }
	  }
	  if(SDDEBUG)
	    printf("After resbias updates:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		   (unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt,LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	  else
	    printf("After resbias updates:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
	  fflush(stdout);

	  if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	    if(VERB){
	      double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
	      for(size_t i = 0; i < numaligns; i++){
		Calign *align = alignment[i];
		if(!align || align->numpairs <= 1)
		  continue;
		int mid = align->mapid2;
		if(DEBUG>=2) assert(!Gmap[mid]->origmap);
		if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		  continue;
		if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		  continue;
		sum1 += origLRarray[mid];
		sum2 += logLRarray[mid];
		sum3 += origSDarray[0][mid] + origSDarray[1][mid];
		sum4 += logSDarray[0][mid] + logSDarray[1][mid];
		printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
		       i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
		       origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
		       logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
	      }
	      fflush(stdout);
	    }
	      
	    assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	  }

	  if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	    printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	    assert(lcnt == tcnt);
	  }
	  logLRsum = LRsum;
	  if(SDDEBUG)
	    logLRsdsum = LRsdsum;
	  score_free(0,numrefmaps);

	  delete [] origLRarray;
	  for(int c = 0; c < colors; c++){
	    delete [] origSDarray[c];
	  }
	}

	for(int c = 0; c < colors; c++){
	  delete [] SizeToBin[c]; SizeToBin[c] = NULL;
	}
      } // if(REFDEBUG_STRICT < 2 && maxresbias > mres * 0.5) */

      /* Print out interval data for the very last iteration */
      if(giter == RefRepeats-1 && giter2 == RefRepeats2 - 1) {
	/* HERE HERE : copy & modify code from refalign.cpp */
      }

      double C = 1.0;

      /* generate optimized scaling factor C for x : minimize MSE(x,y/C) */
      /* Note: There is a small possibility that the current alignment will become invalid if x is scaled up to overlap another site on Y.
	 The error should be small, especially if ENDFIX is used. 
	 REFDEBUG >= 2 requires special care when scaling up of maps X */

      double logXtheta[MAXCOLOR+1];
      for(int c = 0; c < colors;c++)
	logXtheta[c] = 2.0*log(Xtheta[c]) - log(2.0*M_PI);/* to match terms LLsd in logLR() */
      logXtheta[colors] = 2.0*log(Xtheta12) - log(2.0*M_PI);

      if(VERB>=2 && SDDEBUG){
	printf("  giter=%d:MA[2]=%0.8f,MA_mean=%0.8f,SF[2]*SF[2]*2 = %0.8f,FM2var = %0.8f, fabs(SD[2])*SD[2]*2= %0.8f, SM2var= %0.8f, 0.5*logXtheta[2]=%0.8f, RefPen12 - 0.5*log(2)=%0.8f, log(Xtheta12) - 0.5*log(2*PI)=%0.8f\n",
	       giter,MA[2],MA_mean,SF[2]*SF[2]*2.0, FM2var, fabs(SD[2])*SD[2]*2.0, SM2var, 0.5*logXtheta[2], RefPen12 - 0.5*log(2.0), log(Xtheta12) - 0.5*log(2.0*M_PI));
	printf("  RefPen12=%0.8f,Xtheta12=%0.8f,log(Xtheta12) - 0.5*log(PI)=%0.8f\n",RefPen12,Xtheta12,log(Xtheta12) - 0.5*log(M_PI));
	fflush(stdout);
	fflush(stdout);
      }

      for(int iter=0; iter < 30; iter++){
	/* Try to scale x (or C) by factor Cscale : 
	   Note that A & B & R & E will also be scaled by Cscale^2, while log(A+By+Ryy+E*resvar)/Xtheta is unchanged since Xtheta is also scaled by Cscale.
	   PixelLen will be scaled by Cscale (hence res,resSD are scaled by 1/Cscale, so that all Pn() terms remain the same)
	*/
	double Lsum = 0.0;/* sum of log(A+By+Ryy+E*resvar)-2*log(Xtheta) excluding cases with perr->end */
	double errsum = 0.0;/* sum of (Cx - y - MA * s)^2/(A+BL+RLL+E*resvar) */
	double xysum[3] = {0.0,0.0,0.0};/* sum of Cx(y+MA*s)/(A+BL+RLL+E*resvar) */
	double y2sum[3] = {0.0,0.0,0.0};/* sum of (y+MA*s)^2/(A+BL+RLL+E*resvar) */

	/* NOTE : L is either y (c < 2) or len (c == 2), s == the orientation of X, expressed as -1 or +1 */
	if(AlignedSiteThreshold >= 2){/* NOTE : If REFDEBUG logLR may occasionally get worse if C > 1 due to additional FN sites overlapping with unaligned ends
					 Outliers that get reclassified due to changes in C will always provide unexpected improvement : Previous outliers are excluded and can only improve over the
					 fixed outlier score.
					 Current non-outliers that become outliers are scored as non-outliers and can only improve by being treated as outliers.
					 This will only show up after new alignments are computed, so to guarantee improvement in logLR, limit each Cscale so cumulative C value is never greater than 1 */

	  if(VERB>=2 && SDDEBUG){
	    printf("giter=%d,iter=%d:Estimating rescaling factor for Xmaps:\n",giter,iter);
	    printf("  MA[2]=%0.8f,MA_mean=%0.8f,SF[2]*SF[2]*2 = %0.8f,FM2var = %0.8f, fabs(SD[2])*SD[2]*2= %0.8f, SM2var= %0.8f, 0.5*logXtheta[2]=%0.8f, RefPen12 - 0.5*log(2)=%0.8f, log(Xtheta12) - 0.5*log(2*PI)=%0.8f\n",
		   MA[2],MA_mean,SF[2]*SF[2]*2.0, FM2var, fabs(SD[2])*SD[2]*2.0, SM2var, 0.5*logXtheta[2], RefPen12 -0.5*log(2.0), log(Xtheta12) - 0.5*log(2.0*M_PI));
	    fflush(stdout);
	  }

	  for(int c = 0; c < colors; c++){/* MA[c] == 0 */
	    
	    Cinterval *perr = errors[c];

	    for(int tid = 0; tid < numthreads; tid++)
	      Tarray1[tid] = Tarray2[tid] = Tarray3[tid] = Tarray4[tid] = 0.0;

	    // NOTE : non-deterministic with gcc 4.6, even if Tarray*[] are sorted by value, unless -fno-unsafe-math-optimizations is used
	    #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
	    {
	      int tid = 0;
              #ifdef _OPENMP
	      tid = omp_get_thread_num ();
              #endif

	      double myLsum = 0.0;
	      double myerrsum = 0.0;
	      double myxysum = 0.0;
	      double myy2sum = 0.0;

              #pragma omp for schedule(static,64)
	      for(int i = 0; i < numerrors[c]; i++){
		double y = perr[i].y;
		double x = C*perr[i].x;
		if(DEBUG>=2) assert(perr[i].y == perr[i].len);
		if(ENDFIX>=3 && perr[i].end && y >= x)
		  continue;
		double var = A[c] + B[c] * y;
		if(QUADRATIC_VARIANCE)
		  var += R[c] * y * y;
		if(RES_VARIANCE)
		  var += E[c] * perr[i].resvar;
		double err = y-x;
		if((VERB>=2 || SDDEBUG) && !(ENDFIX>=2 && perr[i].end))
		  myLsum += log(var) - logXtheta[c];
		double Ivar = 1.0/var;
		myerrsum += err*err*Ivar;
		Ivar *= y;
		myy2sum += y*Ivar;
		myxysum += x*Ivar;
	      }

	      Tarray1[tid] = myLsum;
	      Tarray2[tid] = myerrsum;
	      Tarray3[tid] = myy2sum;
	      Tarray4[tid] = myxysum;
	    }
	    if(VERB>=2 || SDDEBUG)
	      qsort(Tarray1, numthreads, sizeof(double), (intcmp *)DoubleInc);
	    qsort(Tarray2, numthreads, sizeof(double), (intcmp *)DoubleInc);
	    qsort(Tarray3, numthreads, sizeof(double), (intcmp *)DoubleInc);
	    qsort(Tarray4, numthreads, sizeof(double), (intcmp *)DoubleInc);

	    for(int tid = 0; tid < numthreads; tid++){
	      if(VERB>=2 || SDDEBUG)
		Lsum += Tarray1[tid];
	      errsum += Tarray2[tid];
	      y2sum[c] += Tarray3[tid];
	      xysum[c] += Tarray4[tid];
	      if(VERB>=2){
		printf("c=%d,tid=%d:errsum=%0.18f(delta=%0.18f),y2sum= %0.18f(delta=%0.18f),xysum=%0.18f(delta=%0.18f)\n",
		       c,tid,errsum,Tarray2[tid],y2sum[c],Tarray3[tid],xysum[c],Tarray4[tid]);
		fflush(stdout);
	      }
	    }
	  }
	  if(colors==2){
	    Cinterval *perr = errors[2];
	    for(int i = numerrors[2]; --i >= 0;){
	      double y = perr[i].y + MA[2] * perr[i].end;
	      double len = perr[i].len;
	      if(COLOR_SHIFT==1)
		len *= C;
	      double x = C*perr[i].x;
	      double var = (A[2]+B[2]*len);
	      double err = x - y;
	      double Ivar = 1.0/var;
	      if(VERB>=2 || SDDEBUG)
		Lsum += log(var);
	      errsum += err*err*Ivar;
	      Ivar *= y;
	      y2sum[2] += y*Ivar;
	      xysum[2] += x*Ivar;
	      if(VERB>=2 && giter==0 && iter==0 && perr[i].mapid==0)
		printf("  i=%d/%d:mapid=%d,end=%d:x=%0.4f,y=%0.4f,MA[2]=%0.4f,len=%0.4f(len=%0.4f)\n",i,numerrors[2],perr[i].mapid,perr[i].end, x, perr[i].y,MA[2],perr[i].len,len);
	    }
	  }
	  if(VERB>=2){
	    for(int c = 0; c <= 2; c++){
	      if(xysum[c] > 0.0)
		printf("iter=%d,c=%d:Cscale=%0.6f\n",iter,c,y2sum[c]/xysum[c]);
	    }
	    fflush(stdout);
	  }

	  double Cscale = (y2sum[0]+y2sum[1]+y2sum[2])/(xysum[0]+xysum[1]+xysum[2]);
	  if(VERB>=2){
	    printf("iter=%d:cnt=%d+%d:LL=%0.6f,Cscale=%0.6f: C = %0.6f -> %0.6f\n",
		   iter,numerrors[0],numerrors[1],-(Lsum+errsum)*0.5, Cscale,C, C*Cscale);
	    fflush(stdout);
	  }
	  if(DEBUG>=2 && (VERB>=2 || SDDEBUG) && !isfinite(-(Lsum+errsum)*0.5)){
	    printf("y2sum=%0.6f,xysum=%0.6f,Lsum=%0.6f,errsum=%0.6f\n",
		   y2sum[0]+y2sum[1]+y2sum[2],xysum[0]+xysum[1]+xysum[2],Lsum,errsum);
	    for(int c = 0; c < colors; c++)
	      printf("c=%d:A=%0.6f,B=%0.6f,numerrors=%d\n",c,A[c],B[c],numerrors[c]);
	    fflush(stdout);
	    exit(1);
	  }
	  if(VERB && (VERB>=2 || SDDEBUG)){
	    double delta = (y2sum[0]+y2sum[1]+y2sum[2])*(1.0-1.0/(Cscale*Cscale)) - 2.0*(xysum[0]+xysum[1]+xysum[2])*(1.0-1.0/Cscale);
	    for(int c = 0; c < colors; c++)
	      printf("iter=%d,c=%d:sf=%0.8f -> %0.8f,sd=%0.8f -> %0.8f,sr=%0.8f -> %0.8f,se=%0.8f -> %0.8f:y2sum=%0.8f,xysum=%0.8f,LLsd=%0.6f -> %0.6f,Cscale=%0.8f:C->%0.8f\n",
		     iter,c,sqrt(A[c]),sqrt(A[c])*Cscale,copysign(sqrt(fabs(B[c])),B[c]),copysign(sqrt(fabs(B[c])),B[c])*Cscale,sqrt(R[c]),sqrt(R[c])*Cscale,sqrt(E[c]),sqrt(E[c])*Cscale,
		     y2sum[0]+y2sum[1]+y2sum[2],xysum[0]+xysum[1]+xysum[2],-(Lsum+errsum)*0.5,-0.5*(Lsum+errsum-delta),Cscale,C*Cscale);
	    fflush(stdout);
	  }

	  if(DEBUG && REFDEBUG && REFDEBUG_STRICT >= 2) assert(C <= 1.0);
	  C *= Cscale;
	  if(REFDEBUG && REFDEBUG_STRICT >= 2 && C > 1.0){
	    double correction = 1.0/C;
	    C *= correction;
	    if(DEBUG) assert(Cscale > 1.0);
	    Cscale *= correction;
	    if(DEBUG) assert(Cscale >= 1.0);
	  }
	  for(int c = 0; c <= colors; c++){
	    A[c] *= Cscale*Cscale;
	    B[c] *= Cscale*Cscale;
	    if(QUADRATIC_VARIANCE)
	      R[c] *= Cscale*Cscale;
	    if(RES_VARIANCE)
	      E[c] *= Cscale*Cscale;
	    if(REFDEBUG_STRICT < 2){
	      if(c < colors)
		A[c] = min(A[c],Amax);
	      A[c] = max(A[c],Amin);
	      if(c < colors)
		B[c] = min(B[c],Bmax);
	      B[c] = max(B[c],Bmin);
	      R[c] = min(R[c],Rmax);
	      R[c] = max(R[c],Rmin);
	      E[c] = max(E[c],Emin);
	      E[c] = min(E[c],Emax);
	    }
	    if(DEBUG) assert(E[c] >= 0.0);
	    if(DEBUG) assert(R[c] >= 0.0);
	    if(DEBUG) assert(A[c] >= 0.0);
	    if(R[c] >= 0.0 && A[c] >= 0.0)
	      B[c] = max(B[c], -2.0*sqrt(A[c]*R[c]) * MINSD_MULT * MINSD_MULT);
	  }// c = 0..colors 

	  if(REFDEBUG>=2 && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	    for(int c = 0; c < colors; c++){
	      SF[c] = sqrt(A[c]);
	      SD[c] = copysign(sqrt(fabs(B[c])),B[c]);
	      SR[c] = sqrt(R[c]);
	      SE[c] = sqrt(E[c]);
	    }
	    SF[2] = sqrt(A[2]);
	    SD[2] = sqrt(B[2]);
	    if(SD[2] < 0.001)
	      SD[2] = 0.001;
	    MA_mean = MA[2];
	    
	    if(Cscale != 1.0){ /* rescale all query map sizes by Cscale */
	      double C = Cscale;/* local value of C ! */
	      if(VERB/* HERE >=2 */){
		printf("    Rescaling all map site[] and rawsite[] values  by C=%0.12f (reflecting change in bpp)\n",C);
		fflush(stdout);
	      }

              #pragma omp parallel for num_threads(numthreads) schedule(dynamic,256) if(numthreads > 1)
	      for(int i= 0; i < nummaps; i++){
		Cmap *pmap = Gmap[i];
		if(DEBUG/* HERE >=2 */) assert(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3);
		for(int c = 0; c < colors; c++){
		  FLOAT *X = pmap->site[c];
		  int M = pmap->numsite[c];
		  for(int j= M+1; j > 0; j--)
		    X[j] *= C;
		}
		if(DEBUG/* HERE >=2 */) assert(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3);	    
	      }

	      if(maxresbias > mres * 0.5){/* also rescale rawsite[]. NOTE : applied to all maps including -subset maps (if -ScanScaling) and split maps */

		if(DEBUG) assert(rawsitemaps >= totalmaps);

                #pragma omp parallel for num_threads(numthreads) schedule(dynamic,256) if(numthreads > 1)
		for(int i= 0; i < totalmaps; i++){
		  Cmap *pmap = Gmap[i];
		  if(DEBUG/* HERE >=2 */) assert(fabs(pmap->rawsite[0][pmap->numsite[0]+1] - pmap->rawsite[1][pmap->numsite[1]+1]) < 1e-3);
		  for(int c = 0; c < colors; c++){
		    FLOAT *X = pmap->rawsite[c];
		    int M = pmap->numsite[c];
		    for(int j= M+1; j > 0; j--)
		      X[j] *= C;
		  }
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
	      }

	      double InvC = 1.0/C;
	      for(int i = 0; i < nummaps; i++)
		Gmap[i]->incscale *= InvC;
	      PixelLen *= C;
	      for(int c = 0; c < colors; c++){
		FP[c] *= InvC;
		Xtheta[c] *= C;
		res[c] *= InvC;
		resSD[c] *= InvC;

		if(REFDEBUG_STRICT < 2){
		  FP[c] = max(FP[c],MINFP);
		  FP[c] = min(FP[c],MAXFP);
		}
	      }
	      Xtheta12 *= C;

	      for(int c = 0; c < colors;c++)
		logXtheta[c] = 2.0*log(Xtheta[c]) - log(2.0*M_PI);/* to match terms LLsd in logLR() */
	      logXtheta[colors] = 2.0*log(Xtheta12) - log(2.0*M_PI);
	    }

	    if(DEBUG) assert(logLRarray);
	    (void)score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	    double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	    double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	    lcnt= 0;
	    double *origLRarray = new double[nummaps];
	    double *origSDarray[MAXCOLOR];
	    for(int mid = 0; mid < nummaps; mid++)
	      origLRarray[mid] = logLRarray[mid];
	    for(int c = 0; c < colors; c++){
	      origSDarray[c] = new double[nummaps];
	      for(int mid = 0; mid < nummaps; mid++){
		origSDarray[c][mid] = logSDarray[c][mid];
	      }
	    }

	    double origLRsum = 0.0;
	    double origSDsum[MAXCOLOR] = {0.0,0.0};

	    for(size_t i = 0; i < numaligns; i++){
	      double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	      Calign *align = alignment[i];
	      if(!align || align->numpairs <= 1)
		continue;
	      int mid = align->mapid2;
	      if(DEBUG) assert(!Gmap[mid]->origmap);
	      if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
		if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
		continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	      }
	      if(align->score <= origScoreThreshold)
		continue;
	      int rid = align->mapid1;
	      if(DEBUG>=2){
		int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
		assert(origgoodalign);
		int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
		assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	      }
	      FLOAT **XX = Gmap[mid]->site;
	      int *MM = Gmap[mid]->numsite;
	      FLOAT **YY = refmap[rid]->site;
	      //	int *NN = Gmap[rid]->numsite;
	      align++;
	      logLRarray[mid] = 0.0;
	      for(int c = 0; c < colors; c++){
		FLOAT *X = XX[c];
		double M = MM[c];
		FLOAT *Y = YY[c];
		int N = refmap[rid]->numsite[c];
		int U = align[c].numpairs;

		LR[c] = logLR(mid,rid,c,X,Y,M,N,align,align[c].sites1[0], align[c].sites2[0], align[c].sitesK1[0], align[c].sites2[U-1], align[c].Lij1,align[c].Rij1,align[c].Lij2,align[c].Rij2,1,0,SDDEBUG ? &LRsd[c] : 0, NULL);
		logLRarrayC[c][mid] = LR[c];
		logSDarray[c][mid] = LRsd[c];
		if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
		  continue;
		LRsum += LR[c];
		if(SDDEBUG){
		  LRsdsum += LRsd[c];
		  LRsdsumC[c] += LRsd[c];
		}
	      }

	      int U = align[0].numpairs;
	      int I = align[0].sites1[U-1];
	      int K = align[0].sitesK1[U-1];
	      int J = align[0].sites2[U-1];
	      U = align[1].numpairs;
	      int P = align[1].sites1[U-1];
	      int T = align[1].sitesK1[U-1];
	      int Q = align[1].sites2[U-1];
	      double y0 = Yc(YY[0],I,K);
	      double y1 = Yc(YY[1],P,T);
	      double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	      double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	      double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);
	      
	      LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
	      if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		LRsum += LR[2];
		for(int lc = 0; lc <= colors; lc++)
		  LRsumC[lc] += LR[lc];
	      }

	      logLRarray[mid] = LR[0] + LR[1] + LR[2];
	      lcnt++;
	      if(VERB>=2 && giter==0){
		if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		  origLRsum += origLRarray[mid];
		  for(int lc = 0; lc < colors; lc++)
		    origSDsum[lc] += origSDarray[lc][mid];
		}
		if(mid==265){
		  printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)(delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f)\n",
			 i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
			 logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
		  printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
			 logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
			 LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
			 LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
		  printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);
		}
	      }
	    }
	    if(SDDEBUG)
	      printf("iter=%d:After rescaling Xmaps by %0.12f:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		     iter,Cscale,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt,LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	    else
	      printf("iter=%d:After rescaling Xmaps by %0.12f:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",iter,Cscale,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
	    fflush(stdout);

	    if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	      if(VERB){
		double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
		for(size_t i = 0; i < numaligns; i++){
		  Calign *align = alignment[i];
		  if(!align || align->numpairs <= 1)
		    continue;
		  int mid = align->mapid2;
		  if(DEBUG>=2) assert(!Gmap[mid]->origmap);
		  if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		    continue;
		  if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		    continue;
		  sum1 += origLRarray[mid];
		  sum2 += logLRarray[mid];
		  sum3 += origSDarray[0][mid] + origSDarray[1][mid];
		  sum4 += logSDarray[0][mid] + logSDarray[1][mid];
		  printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
			 i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
			 origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
			 logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
		}
		fflush(stdout);
	      }
	      
	      assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	    }

	    if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	      printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	      assert(lcnt == tcnt);
	    }
	    logLRsum = LRsum;
	    if(SDDEBUG)
	      logLRsdsum = LRsdsum;
	    score_free(0,numrefmaps);

	    delete [] origLRarray;
	    for(int c = 0; c < colors; c++){
	      delete [] origSDarray[c];
	    }
	  }// if(REFDEBUG>=2 ... )
	} // if(AlignedSiteThreshold >= 2)

	if(VERB>=2 && SDDEBUG){
	  printf("giter=%d,iter=%d:Estimating MA[2]:\n",giter,iter);
	  printf("  MA[2]=%0.8f,MA_mean=%0.8f,SF[2]*SF[2]*2 = %0.8f,FM2var = %0.8f, fabs(SD[2])*SD[2]*2= %0.8f, SM2var= %0.8f, 0.5*logXtheta[2]=%0.8f, RefPen12 - 0.5*log(2)=%0.8f, log(Xtheta12) - 0.5*log(2*PI)=%0.8f\n",
		 MA[2],MA_mean,SF[2]*SF[2]*2.0, FM2var, fabs(SD[2])*SD[2]*2.0, SM2var, 0.5*logXtheta[2], RefPen12 - 0.5*log(2.0), log(Xtheta12) - 0.5*log(2.0*M_PI));
	  fflush(stdout);
	}

	/* update MA[2] */
	for(int c = 2; c <= 2; c++){
	  Cinterval *perr = errors[c];
	  double Ivarsum = 0.0;
	  double Serrsum = 0.0;

	  double Lcnt = 0.0;
	  double Lsum = 0.0;/* sum of log(A+BL) */
	  double errsum = 0.0;/* sum of (Cx-y-MS[2]*s)^2/(A+BL) */

	  for(int i = numerrors[c]; --i >= 0;){
	    double y = perr[i].y;
	    double len = perr[i].len;
	    if(COLOR_SHIFT==1)
	      len *= C;
	    double x = C * perr[i].x;
	    double var = A[c]+B[c]*len;
	    double err = x - y;
	    double Ivar = 1.0/var;
	    Serrsum += perr[i].end * err * Ivar;
	    Ivarsum += Ivar;
	    if(VERB && (VERB >= 2 || SDDEBUG)){
	      Lcnt += 1.0;
	      Lsum += log(var);
	      double terr = (err - MA[2]*perr[i].end);
	      errsum += terr*terr*Ivar;
	    }
	    if(VERB>=2 && giter==0 && iter==0 && perr[i].mapid==0)
	      printf("  i=%d/%d:mapid=%d,end=%d:x=%0.4f,y=%0.4f,MA[2]=%0.4f,len=%0.4f(%0.4f)\n",i,numerrors[c],perr[i].mapid,perr[i].end, x, y,MA[2],perr[i].len,len);
	  }
	  if(Ivarsum > 0.0){
	    double oldMA = MA[2];
	    MA[2] = Serrsum/Ivarsum;
	    if(VERB && (VERB >= 2 || SDDEBUG)){
	      Cinterval *perr = errors[c];

	      double Lcnt2 = 0.0;
	      double Lsum2 = 0.0;/* sum of log(A+BL) */
	      double errsum2 = 0.0;/* sum of (Cx-y-MS[2]*s)^2/(A+BL) */
	      
	      for(int i = numerrors[c]; --i >= 0;){
		double x = C * perr[i].x;
		double y = perr[i].y;
		double len = perr[i].len;
		if(COLOR_SHIFT)
		  len *= C;
		double var = A[c]+B[c]*len;
		double err = x - y;
		double Ivar = 1.0/var;
		Lcnt2 += 1.0;
		Lsum2 += log(var);
		double terr = (err - MA[2]*perr[i].end);
		errsum2 += terr*terr*Ivar;
		if(VERB>=2 && giter==0 && iter==0  && perr[i].mapid==0){
		  printf("  A:i=%d/%d:mapid=%d,end=%d:x=%0.4f,y=%0.4f,MA[2]=%0.4f,len=%0.4f(%0.4f)\n",i,numerrors[c],perr[i].mapid,perr[i].end, x, perr[i].y, MA[2], perr[i].len,len);
		  fflush(stdout);
		}
	      }

	      printf("iter=%d,c=2:Serrsum=%0.6f,Ivarsum=%0.6f,cnt=%0.1f,LL=%0.6f->%0.6f:MA_mean=%0.4f -> %0.4f\n",
		     iter,Serrsum,Ivarsum,Lcnt,-(Lsum+errsum)*0.5,-(Lsum2+errsum2)*0.5,oldMA,MA[2]);
	      fflush(stdout);
	    }
	  }
	}

	if(REFDEBUG>=2 && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	  for(int c = 0; c < colors; c++){
	    SF[c] = sqrt(A[c]);
	    SD[c] = copysign(sqrt(fabs(B[c])),B[c]);
	    SR[c] = sqrt(R[c]);
	    SE[c] = sqrt(E[c]);
	  }
	  SF[2] = sqrt(A[2]);
	  SD[2] = sqrt(B[2]);
	  if(SD[2] < 0.001)
	    SD[2] = 0.001;
	  MA_mean = MA[2];

	  if(DEBUG) assert(logLRarray);
	  score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	  double LRsum = 0.0,LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	  double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	  lcnt= 0;
	  double *origLRarray = new double[nummaps];
	  double *origSDarray[MAXCOLOR];
	  for(int mid = 0; mid < nummaps; mid++)
	    origLRarray[mid] = logLRarray[mid];
	  for(int c = 0; c < colors; c++){
	    origSDarray[c] = new double[nummaps];
	    for(int mid = 0; mid < nummaps; mid++){
	      origSDarray[c][mid] = logSDarray[c][mid];
	    }
	  }

	  double origLRsum = 0.0;
	  double origSDsum[MAXCOLOR] = {0.0,0.0};

	  for(size_t i = 0; i < numaligns; i++){
	    double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	    Calign *align = alignment[i];
	    if(!align || align->numpairs <= 1)
	      continue;
	    int mid = align->mapid2;
	    if(DEBUG) assert(!Gmap[mid]->origmap);
	    if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
	      if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
	      continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	    }
	    if(align->score <= origScoreThreshold)
	      continue;
	    int rid = align->mapid1;
	    if(DEBUG>=2){
	      int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
	      assert(origgoodalign);
	      int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
	      assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	    }
	    FLOAT **XX = Gmap[mid]->site;
	    int *MM = Gmap[mid]->numsite;
	    FLOAT **YY = refmap[rid]->site;
	    //	int *NN = Gmap[rid]->numsite;
	    align++;
	    logLRarray[mid] = 0.0;
	    for(int c = 0; c < colors; c++){
	      FLOAT *X = XX[c];
	      double M = MM[c];
	      FLOAT *Y = YY[c];
	      int N = refmap[rid]->numsite[c];
	      int U = align[c].numpairs;

	      LR[c] = logLR(mid,rid,c,X,Y,M,N,align,align[c].sites1[0], align[c].sites2[0], align[c].sitesK1[0], align[c].sites2[U-1], align[c].Lij1,align[c].Rij1,align[c].Lij2,align[c].Rij2,1,
			    0 /* (giter==0 && iter==1 && c==0 && mid==265) ? 1 : 0 */, SDDEBUG ? &LRsd[c] : 0,NULL);
	      logLRarrayC[c][mid] = LR[c];
	      logSDarray[c][mid] = LRsd[c];
	      if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
		continue;
	      LRsum += LR[c];
	      if(SDDEBUG){
		LRsdsum += LRsd[c];
		LRsdsumC[c] += LRsd[c];
	      }
	    }

	    int U = align[0].numpairs;
	    int I = align[0].sites1[U-1];
	    int K = align[0].sitesK1[U-1];
	    int J = align[0].sites2[U-1];
	    U = align[1].numpairs;
	    int P = align[1].sites1[U-1];
	    int T = align[1].sitesK1[U-1];
	    int Q = align[1].sites2[U-1];
	    double y0 = Yc(YY[0],I,K);
	    double y1 = Yc(YY[1],P,T);
	    double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	    double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	    double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);

	    LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
	    if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	      LRsum += LR[2];
	      for(int lc = 0; lc <= colors; lc++)
		LRsumC[lc] += LR[lc];
	    }

	    logLRarray[mid] = LR[0] + LR[1] + LR[2];
	    lcnt++;
	    if(VERB>=2 && giter==0){
	      if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		origLRsum += origLRarray[mid];
		for(int lc = 0; lc < colors; lc++)
		  origSDsum[lc] += origSDarray[lc][mid];
	      }
	      if(iter==1 || mid==265){
		printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)(delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f)\n",
		       i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
		       logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
		printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
		       logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
		       LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
		       LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
		printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);
	      }
	    }
	  }

	  if(SDDEBUG)
	    printf("iter=%d:After MA[2] update:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		   iter,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt,LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	  else
	    printf("iter=%d:After MA[2] update:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",iter,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
	  fflush(stdout);

	  if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	    if(VERB){
	      double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
	      for(size_t i = 0; i < numaligns; i++){
		Calign *align = alignment[i];
		if(!align || align->numpairs <= 1)
		  continue;
		int mid = align->mapid2;
		if(DEBUG>=2) assert(!Gmap[mid]->origmap);
		if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		  continue;
		if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		  continue;
		sum1 += origLRarray[mid];
		sum2 += logLRarray[mid];
		sum3 += origSDarray[0][mid] + origSDarray[1][mid];
		sum4 += logSDarray[0][mid] + logSDarray[1][mid];
		printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
		       i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
		       origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
		       logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
	      }
	      fflush(stdout);
	    }
	      
	    assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	  }

	  if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	    printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	    assert(lcnt == tcnt);
	  }
	  logLRsum = LRsum;
	  if(SDDEBUG)
	    logLRsdsum = LRsdsum;
	  score_free(0,numrefmaps);

	  delete [] origLRarray;
	  for(int c = 0; c < colors; c++){
	    delete [] origSDarray[c];
	  }
	}// if(REFDEBUG>=2 .. )

	for(int c = 0; c <= colors; c++){
	  /* Try to scale A[c], B[C], R[c] & E[c] by the same factor */
	  Lsum = 0.0;
	  errsum = 0.0;
	  Cinterval *perr = errors[c];
	  double Lcnt = 0.0;/* count of terms in errsum */

	  if(c < colors){
	    for(int tid = 0; tid < numthreads; tid++)
	      Tarray1[tid] = Tarray2[tid] = Tarray3[tid] = 0.0;

	    if(VERB>=2 && SDDEBUG){
	      printf("giter=%d,iter=%d,c=%d:Estimating shared scaling factor for SF[c]=%0.8f,SD[c]=%0.8f,SR[c]=%0.8f,SE[c]=%0.8f\n",giter,iter,c,SF[c],SD[c],SR[c],SE[c]);
	      printf("  MA[2]=%0.8f,MA_mean=%0.8f,SF[2]*SF[2]*2 = %0.8f,FM2var = %0.8f, fabs(SD[2])*SD[2]*2= %0.8f, SM2var= %0.8f, 0.5*logXtheta[2]=%0.8f, RefPen12 - 0.5*log(2)=%0.8f, log(Xtheta12) - 0.5*log(2*PI)=%0.8f\n",
		     MA[2],MA_mean,SF[2]*SF[2]*2.0, FM2var, fabs(SD[2])*SD[2]*2.0, SM2var, 0.5*logXtheta[2], RefPen12 - 0.5*log(2.0), log(Xtheta12) - 0.5*log(2.0*M_PI));
	      fflush(stdout);
	    }

	    // NOTE : non-deterministic with gcc 4.6, even if Tarray*[] are sorted by value, unless -fno-unsafe-math-optimizations is used
            #pragma omp parallel num_threads(numthreads) if(numthreads > 1 && !SDDEBUG)
	    {
	      int tid = 0;
              #ifdef _OPENMP
	      tid = omp_get_thread_num ();
              #endif

	      double myLsum = 0.0;
	      double myLcnt = 0.0;
	      double myerrsum = 0.0;
	  
              #pragma omp for schedule(static,64)
	      for(int i = 0; i < numerrors[c]; i++){
		double y = perr[i].y;
		double x = C * perr[i].x;
		if(ENDFIX>=3 && perr[i].end && y >= x)
		  continue;
		double err = y - x;
		double var = A[c] + B[c]*y;
		if(QUADRATIC_VARIANCE)
		  var += R[c] * y*y;
		if(RES_VARIANCE)
		  var += E[c] * perr[i].resvar;
		if(!(ENDFIX>=2 && perr[i].end)){
		  if(VERB >= 2 || SDDEBUG){
		    myLsum += log(var);
		    if(SDDEBUG)
		      myLsum -= logXtheta[c];
		  }
		  myLcnt += 1.0;
		}
		myerrsum += err*err/var;
		if(VERB>=2 && SDDEBUG && giter==0 && iter==1 && c==0){
		  printf("iter=%d,c=%d:i=%d/%d:mid=%d,L=%d,R=%d,end=%d,x=%0.4f,y=%0.4f,var=%0.8f,Lcnt=%0.1f,Lsum=%0.8f(delta=%0.8f),errsum=%0.8f(delta=%0.8f):LL=%0.8f(delta=%0.8f)\n",
			 iter,c,i,numerrors[c],perr[i].mapid,perr[i].L,perr[i].R,perr[i].end,x,y,var,Lcnt,myLsum,log(var) - logXtheta[c],myerrsum,err*err/var,
			 -(myLsum + myerrsum)*0.5, -(log(var) - logXtheta[c] + err*err/var)*0.5);
		  fflush(stdout);
		}
	      }
	      
	      if(VERB>=2 || SDDEBUG)
		Tarray1[tid] = myLsum;
	      Tarray2[tid] = myLcnt;
	      Tarray3[tid] = myerrsum;
	    }

	    if(VERB>=2 || SDDEBUG)
	      qsort(Tarray1, numthreads, sizeof(double), (intcmp *)DoubleInc);
	    qsort(Tarray2, numthreads, sizeof(double), (intcmp *)DoubleInc);
	    qsort(Tarray3, numthreads, sizeof(double), (intcmp *)DoubleInc);
	    for(int tid = 0; tid < numthreads; tid++){
	      if(VERB>=2 || SDDEBUG)
		Lsum += Tarray1[tid];
	      Lcnt += Tarray2[tid];
	      errsum += Tarray3[tid];
	    }
	  } else {// c == 2 
	    if(VERB>=2){
	      printf("  MA[2]=%0.8f,MA_mean=%0.8f,SF[2]*SF[2]*2 = %0.8f,FM2var = %0.8f, fabs(SD[2])*SD[2]*2= %0.8f, SM2var= %0.8f, 0.5*logXtheta[2]=%0.8f, RefPen12 - 0.5*log(2)=%0.8f, log(Xtheta12) - 0.5*log(2*PI)=%0.8f\n",
		     MA[2],MA_mean,SF[2]*SF[2]*2.0, FM2var, fabs(SD[2])*SD[2]*2.0, SM2var, 0.5*logXtheta[2], RefPen12 - 0.5*log(2.0), log(Xtheta12) - 0.5*log(2.0*M_PI));
	      fflush(stdout);
	    }

	    for(int i = 0; i < numerrors[c]; i++){
	      double y = perr[i].y + MA[c] * perr[i].end;/* WAS perr->y - MA[c] */
	      double x = C * perr[i].x;
	      double len = perr[i].len;
	      if(COLOR_SHIFT==1)
		len *= C;
	      double var = A[c]+B[c]*len;
	      double err = x - y;
	      if(VERB >= 2 || SDDEBUG)
		Lsum += log(var) - logXtheta[c];
	      Lcnt += 1.0;
	      errsum += err*err/var;
	      if(VERB>=2 && SDDEBUG && giter==0 && iter==1)
		printf("  i=%d/%d:mapid=%d,end=%d:x=%0.4f,y=%0.4f,len=%0.4f(%0.4f):SMA(x,y,len,or=%d)=%0.6f,LR2=%0.6f(cum=%0.6f)\n",i,numerrors[2],perr[i].mapid,perr[i].end, x, perr[i].y,perr[i].len,len,
		       (perr[i].end < 0) ? 1 : 0, SMA(x, perr[i].y,len,(perr[i].end < 0) ? 1 : 0), -0.5*(log(var) - logXtheta[c] + err*err/var), -0.5*(Lsum + errsum));
	    }
	  }

	  double delta = errsum/Lcnt;
	  if(VERB && (VERB >= 2 || SDDEBUG)){
	    if(c < colors)
	      printf("iter=%d,c=%d:tcnt=%lu,Lcnt=%0.1f/%d:LL=%0.6f(LR[c]=%0.6f -> %0.6f),norm=%0.6f:sf=%0.6f->%0.6f,sd=%0.6f->%0.6f,sr=%0.6f->%0.6f,se=%0.6f->%0.6f\n",
		     iter,c,tcnt,Lcnt,numerrors[c],-(Lsum+errsum)*0.5,-(Lsum+errsum)*0.5/tcnt,-(Lsum + Lcnt*log(delta) + errsum/delta)*0.5/tcnt,
		     delta,sqrt(A[c]),sqrt(A[c]*delta),copysign(sqrt(fabs(B[c])),B[c]),copysign(sqrt(fabs(B[c]*delta)),B[c]),
		     sqrt(R[c]),sqrt(R[c]*delta),sqrt(E[c]),sqrt(E[c]*delta));
	    else
	      printf("iter=%d,c=%d:tcnt=%lu,Lcnt=%0.1f/%d:LL=%0.6f(LR[2]=%0.6f -> %0.6f),norm=%0.6f:sf=%0.6f->%0.6f,sd=%0.6f->%0.6f\n",
		     iter,c,tcnt,Lcnt,numerrors[c],-(Lsum+errsum)*0.5,-(Lsum+errsum)*0.5/tcnt,-((Lsum + Lcnt * log(delta)) + errsum/delta)*0.5/tcnt,
		     delta,sqrt(A[c]),sqrt(A[c]*delta),copysign(sqrt(fabs(B[c])),B[c]),copysign(sqrt(fabs(B[c]*delta)),B[c]));
	    fflush(stdout);
	  }

	  if(VERB>=2 && SDDEBUG && giter==0 && iter==1 && c==2){/* print out SMA() values before AND after delta correction */
	    printf("FM2var=%0.8f -> %0.8f,SM2var=%0.8f -> %0.8f,RefPen12=%0.8f\n",FM2var,FM2var*delta,SM2var,SM2var*delta,RefPen12);
	    fflush(stdout);

	    FM2var *= delta;
	    SM2var *= delta;
	    double errsum = 0.0, errsum2 = 0.0, Lsum = 0.0, Lsum2 = 0.0;
	    for(int i = 0; i < numerrors[c]; i++){
	      double y = perr[i].y + MA[c] * perr[i].end;/* WAS perr->y - MA[c] */
	      double x = C * perr[i].x;
	      double len = perr[i].len;
	      if(COLOR_SHIFT==1)
		len *= C;
	      double var = A[c]+B[c]*len;
	      double err = x - y;
	      if(VERB >= 2 || SDDEBUG){
		Lsum += log(var) - logXtheta[c];
		Lsum2 += log(var)+log(delta) - logXtheta[c];
	      }
	      Lcnt += 1.0;
	      errsum += err*err/var;
	      errsum2 += err*err/(var * delta);

	      printf("  i=%d/%d:mapid=%d,end=%d:x=%0.4f,y=%0.4f,len=%0.4f(%0.4f):SMA(x,y,len,or=%d)=%0.6f,LR2=%0.6f->%0.6f (cum= %0.6f -> %0.6f),change=%0.6f(cum=%0.6f)\n",
		     i,numerrors[2],perr[i].mapid,perr[i].end, x, perr[i].y,perr[i].len,len,
		     (perr[i].end < 0) ? 1 : 0, SMA(x,perr[i].y,len,(perr[i].end < 0) ? 1 : 0),
		     -0.5*(log(var) - logXtheta[c] + err*err/var), -0.5*(log(var)+log(delta) - logXtheta[c] + err*err/(var*delta)), 
		     -0.5*(Lsum + errsum), -0.5*(Lsum2 + errsum2),
		     -0.5*(log(var)+log(delta) - logXtheta[c] + err*err/(var*delta)) + 0.5*(log(var) - logXtheta[c] + err*err/var),
		     -0.5*(Lsum2 + errsum2) + 0.5*(Lsum + errsum));
	    }
	    fflush(stdout);
	  }

	  A[c] *= delta;
	  B[c] *= delta;
	  if(c < colors){
	    if(QUADRATIC_VARIANCE)
	      R[c] *= delta;
	    if(RES_VARIANCE)
	      E[c] *= delta;
	  }
	  if(REFDEBUG_STRICT < 2){
	    if(c < colors)
	      A[c] = min(A[c],Amax);
	    A[c] = max(A[c],Amin);
	    if(c < colors)
	      B[c] = min(B[c],Bmax);
	    B[c] = max(B[c],Bmin);
	    if(c < colors){
	      R[c] = min(R[c],Rmax);
	      R[c] = max(R[c],Rmin);
	      E[c] = min(E[c],Emax);
	      E[c] = max(E[c],Emin);
	    }
	  }
	  if(DEBUG) assert(E[c] >= 0.0);
	  if(DEBUG) assert(R[c] >= 0.0);
	  if(DEBUG) assert(A[c] >= 0.0);
	  if(QUADRATIC_VARIANCE && R[c] >= 0.0 && A[c] >= 0.0)
	    B[c] = max(B[c], -2.0*sqrt(A[c]*R[c]) * MINSD_MULT * MINSD_MULT);

	  if(REFDEBUG>=2 && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	    SF[c] = sqrt(A[c]);
	    SD[c] = copysign(sqrt(fabs(B[c])),B[c]);
	    if(c < colors){
	      SR[c] = sqrt(R[c]);
	      SE[c] = sqrt(E[c]);
	    }

	    if(DEBUG) assert(logLRarray);
	    score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	    double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	    double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	    lcnt= 0;
	    double *origLRarray = new double[nummaps];
	    double *origSDarray[MAXCOLOR];
	    for(int mid = 0; mid < nummaps; mid++)
	      origLRarray[mid] = logLRarray[mid];
	    for(int lc = 0; lc < colors; lc++){
	      origSDarray[lc] = new double[nummaps];
	      for(int mid = 0; mid < nummaps; mid++){
		origSDarray[lc][mid] = logSDarray[lc][mid];
	      }
	    }

	    double origLRsum = 0.0;
	    double origSDsum[MAXCOLOR] = {0.0,0.0};

	    if(VERB>=2 && giter==0 && iter==1 && c==2){
	      printf("FM2var=%0.8f,SM2var=%0.8f,RefPen12=%0.8f\n",FM2var,SM2var,RefPen12);
	      fflush(stdout);
	    }

	    for(size_t i = 0; i < numaligns; i++){
	      double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	      Calign *align = alignment[i];
	      if(!align || align->numpairs <= 1)
		continue;
	      int mid = align->mapid2;
	      if(DEBUG) assert(!Gmap[mid]->origmap);
	      if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
		if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
		continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	      }
	      if(align->score <= origScoreThreshold)
		continue;
	      int rid = align->mapid1;
	      if(DEBUG>=2){
		int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
		assert(origgoodalign);
		int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
		assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	      }
	      FLOAT **XX = Gmap[mid]->site;
	      int *MM = Gmap[mid]->numsite;
	      FLOAT **YY = refmap[rid]->site;
	      //	int *NN = Gmap[rid]->numsite;
	      align++;
	      logLRarray[mid] = 0.0;
	      for(int lc = 0; lc < colors; lc++){
		FLOAT *X = XX[lc];
		double M = MM[lc];
		FLOAT *Y = YY[lc];
		int N = refmap[rid]->numsite[lc];
		int U = align[lc].numpairs;

		LR[lc] = logLR(mid,rid,lc,X,Y,M,N,align,align[lc].sites1[0], align[lc].sites2[0], align[lc].sitesK1[0], align[lc].sites2[U-1], align[lc].Lij1,align[lc].Rij1,align[lc].Lij2,align[lc].Rij2,1,
			       0 /* (SDDEBUG && giter==0 && iter==1 && lc==0 && mid==265) ? 1 : 0 */, SDDEBUG ? &LRsd[lc] : 0,NULL);
		logLRarrayC[lc][mid] = LR[lc];
		logSDarray[lc][mid] = LRsd[lc];
		if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
		  continue;
		LRsum += LR[lc];
		if(SDDEBUG){
		  LRsdsum += LRsd[lc];
		  LRsdsumC[lc] += LRsd[lc];
		}
	      }

	      int U = align[0].numpairs;
	      int I = align[0].sites1[U-1];
	      int K = align[0].sitesK1[U-1];
	      int J = align[0].sites2[U-1];
	      U = align[1].numpairs;
	      int P = align[1].sites1[U-1];
	      int T = align[1].sitesK1[U-1];
	      int Q = align[1].sites2[U-1];
	      double y0 = Yc(YY[0],I,K);
	      double y1 = Yc(YY[1],P,T);
	      double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	      double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	      double len = (COLOR_SHIFT==2) ? max(y0,y1) : max(x0,x1);
	      
	      LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
	      if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		LRsum += LR[2];
		for(int lc = 0; lc <= colors; lc++)
		  LRsumC[lc] += LR[lc];
	      }

	      logLRarray[mid] = LR[0] + LR[1] + LR[2];
	      lcnt++;
	      if(VERB>=2 && giter==0 && iter==1 && c==2){
		if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		  origLRsum += origLRarray[mid];
		  for(int lc = 0; lc < colors; lc++)
		    origSDsum[lc] += origSDarray[lc][mid];
		}
		if(1){
		  printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)(delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f(%0.6f,%0.6f))\n",
			 i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
			 logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
		  printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
			 logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
			 LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
			 LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
		  printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,SMA(x,y,len,or)=%0.6f,lcnt=%llu\n",
			 x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,SMA(x1-x0,y1-y0,len,align[-1].orientation),(unsigned long long)lcnt);
		}
	      }
	    }

	    if(c < colors){
	      if(SDDEBUG)
		printf("iter=%d,c=%d:After rescaling SF[c],SD[c],SR[c],SE[c]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		       iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt, LRsdsumC[0]/lcnt, LRsdsumC[1]/lcnt);
	      else
		printf("iter=%d,c=%d:After rescaling SF[c],SD[c],SR[c],SE[c]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",
		       iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
	    } else {
	      if(SDDEBUG)
		printf("iter=%d,c=%d:After rescaling SF[c],SD[c]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		       iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt, LRsdsumC[0]/lcnt, LRsdsumC[1]/lcnt);
	      else
		printf("iter=%d,c=%d:After rescaling SF[c],SD[c]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",
		       iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
	    }
	    fflush(stdout);

	    if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	      if(VERB){
		double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
		for(size_t i = 0; i < numaligns; i++){
		  Calign *align = alignment[i];
		  if(!align || align->numpairs <= 1)
		    continue;
		  int mid = align->mapid2;
		  if(DEBUG>=2) assert(!Gmap[mid]->origmap);
		  if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		    continue;
		  if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		    continue;
		  sum1 += origLRarray[mid];
		  sum2 += logLRarray[mid];
		  sum3 += origSDarray[0][mid] + origSDarray[1][mid];
		  sum4 += logSDarray[0][mid] + logSDarray[1][mid];
		  printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
			 i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
			 origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
			 logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
		}
		fflush(stdout);
	      }
	      
	      assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	    }

	    if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	      printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	      assert(lcnt == tcnt);
	    }
	    logLRsum = LRsum;
	    if(SDDEBUG)
	      logLRsdsum = LRsdsum;
	    score_free(0,numrefmaps);

	    delete [] origLRarray;
	    for(int c = 0; c < colors; c++){
	      delete [] origSDarray[c];
	    }
	  }// if(REFDEBUG>=2 ... )

	  /* update estimate of A[c] (seperately from B[c],R[c] or E[c]) (3 iterations with delta = -f'/f'') */
	  double bestLL = MINSCORE;
	  double bestA = A[c];

	  if(VERB>=2 && SDDEBUG){
	    printf("giter=%d,iter=%d,c=%d:Estimating scaling factor for SF[c]=%0.8f\n",giter,iter,c,SF[c]);
	    printf("  MA[2]=%0.8f,MA_mean=%0.8f,SF[2]*SF[2]*2 = %0.8f,FM2var = %0.8f, fabs(SD[2])*SD[2]*2= %0.8f, SM2var= %0.8f, 0.5*logXtheta[2]=%0.8f, RefPen12 - 0.5*log(2)=%0.8f, log(Xtheta12) - 0.5*log(2*PI)=%0.8f\n",
		   MA[2],MA_mean,SF[2]*SF[2]*2.0, FM2var, fabs(SD[2])*SD[2]*2.0, SM2var, 0.5*logXtheta[2], RefPen12 - 0.5*log(2.0), log(Xtheta12) - 0.5*log(2.0*M_PI));
	    fflush(stdout);
	  }

	  int AiterMax = 3;
	  for(int Aiter=0; Aiter < AiterMax; Aiter++){
	    double Lsum = 0.0;/* sum of log(A+By+Ryy+E*resvar) */
	    double Isum = 0.0;/* sum of 1/(A+By+Ryy+E*resvar) */
	    double Isum2 = 0.0;/* sum of 1/(A+By+Ryy+E*resvar)^2 */
	    double errsum = 0.0;/* sum of (Cx-y)^2/(A+By+Ryy+E*resvar) */
	    double errsum2 = 0.0;/* sum of (Cx-y)^2/(A+By+Ryy+E*resvar)^2 */
	    double errsum3 = 0.0;/* sum of (Cx-y)^2/(A+By+Ryy+E*resvar)^3 */
	    Cinterval *perr = errors[c];

	    if(c < colors) {
	      for(int tid = 0; tid < numthreads; tid++)
		Tarray1[tid] = Tarray2[tid] = Tarray3[tid] = Tarray4[tid] = Tarray5[tid] = Tarray6[tid] = 0.0;

	      // NOTE : non-deterministic with gcc 4.6, even if Tarray*[] are sorted by value, unless -fno-unsafe-math-optimizations is used
	      #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
	      {
		int tid = 0;
                #ifdef _OPENMP
		tid = omp_get_thread_num ();
                #endif

		double myLsum = 0.0;
		double myIsum = 0.0;
		double myIsum2 = 0.0;
		double myerrsum = 0.0;
		double myerrsum2 = 0.0;
		double myerrsum3 = 0.0;

                #pragma omp for schedule(static,64)
		for(int i = 0; i < numerrors[c]; i++){
		  double y = perr[i].y;
		  double x = C*perr[i].x;
		  if(ENDFIX>=3 && perr[i].end && y >= x)
		    continue;
		  double err = y - x;
		  err *= err;

		  double var = A[c] + B[c] * y;
		  if(QUADRATIC_VARIANCE)
		    var += R[c] * y*y;
		  if(RES_VARIANCE)
		    var += E[c] * perr[i].resvar;

		  double Ivar = 1.0/var;
		  double Ivar2 = Ivar*Ivar;

		  if(!(ENDFIX>=2 && perr[i].end)){
		    myLsum += log(var);
		    myIsum += Ivar;
		    myIsum2 += Ivar2;
		  }

		  myerrsum += err * Ivar;
		  myerrsum2 += err*Ivar2;
		  myerrsum3 += err*Ivar2*Ivar;
		}

		Tarray1[tid] = myLsum;
		Tarray2[tid] = myIsum;
		Tarray3[tid] = myIsum2;
		Tarray4[tid] = myerrsum;
		Tarray5[tid] = myerrsum2;
		Tarray6[tid] = myerrsum3;
	      }

	      qsort(Tarray1, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray2, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray3, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray4, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray5, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray6, numthreads, sizeof(double), (intcmp *)DoubleInc);

	      for(int tid = 0; tid < numthreads; tid++){
		Lsum += Tarray1[tid];
		Isum += Tarray2[tid];
		Isum2 += Tarray3[tid];
		errsum += Tarray4[tid];
		errsum2 += Tarray5[tid];
		errsum3 += Tarray6[tid];
	      }
	      
	    } else {/* c == colors */
	      for(int i = 0; i < numerrors[c]; i++){
		double y = perr[i].y + MA[c] * perr[i].end;
		double x = C*perr[i].x;
		double len = perr[i].len;
		if(COLOR_SHIFT==1)
		  len *= C;
		double var = A[c]+B[c]* len;
		double err = x - y;

		double Ivar = 1.0/var;
		double Ivar2 = Ivar*Ivar;

		Lsum += log(var);
		Isum += Ivar;
		Isum2 += Ivar2;

		err *= err;
		errsum += err * Ivar;
		errsum2 += err * Ivar2;
		errsum3 += err * Ivar2*Ivar;

		if(VERB>=2 && giter==0 && iter==0)
		  printf("  i=%d/%d:mapid=%d,end=%d:x=%0.4f,y=%0.4f,len=%0.4f(len=%0.4f)\n",i,numerrors[2],perr[i].mapid,perr[i].end, x, perr[i].y,perr[i].len,len);

	      }
	    } /* else c==2 */

	    double LL = -(Lsum+errsum)*0.5;
	    double LL1 = 0.5*(errsum2-Isum);
	    double LL2 = 0.5*Isum2-errsum3;
	    if(LL < bestLL){
	      if(VERB && (VERB >= 2 || SDDEBUG)){
		if(c < colors)
		  printf("iter=%d,%d,c=%d,sf=%0.8f,sd=%0.8f,sr=%0.8f,se=%0.8f:LL=%0.6f,LL1=%0.6f,LL2=%0.3f:backtracking (sf=%0.6f->%0.6f)\n",
			 iter,Aiter,c,sqrt(A[c]),copysign(sqrt(fabs(B[c])),B[c]),sqrt(R[c]),sqrt(E[c]),LL,LL1,LL2, sqrt(A[c]),sqrt(0.5*(A[c]+bestA)));
		else
		  printf("iter=%d,%d,c=%d,sf=%0.8f,sd=%0.8f:LL=%0.6f,LL1=%0.6f,LL2=%0.3f:backtracking (sf=%0.6f->%0.6f)\n",
			 iter,Aiter,c,sqrt(A[c]),copysign(sqrt(fabs(B[c])),B[c]),LL,LL1,LL2, sqrt(A[c]),sqrt(0.5*(A[c]+bestA)));
		fflush(stdout);
	      }
	      A[c] = (A[c] + bestA)*0.5;
	      if(LL < bestLL - 1e-6 && fabs(A[c] - bestA) > 1e-6)
		Aiter--;/* add an iteration if drop is LL was significant */
	      continue;
	    }
	    bestA = A[c];
	    bestLL = LL;

	    double maxdelta = Amax * 0.5;
	    double delta = (LL2 < 0.0) ? min(maxdelta,fabs(LL1/LL2)) : maxdelta;
	    delta = copysign(delta,LL1);
	    if(A[c] + delta < 0.0)
	      delta = -A[c];

	    if(VERB && (VERB >= 2 || SDDEBUG)){
	      printf("iter=%d,%d,c=%d,sf=%0.6f,s=%0.4f:LL=%0.6f,LL1=%0.6f,LL2=%0.3f:delta=%0.10f(maxdelta=%0.10f),A[c]=%0.10f->%0.10f(sf=%0.6f->%0.6f)\n",
		     iter,Aiter,c,sqrt(A[c]),copysign(sqrt(fabs(B[c])),B[c]),LL,LL1,LL2, delta, maxdelta,A[c],A[c]+delta,sqrt(A[c]),sqrt(A[c]+delta));
	      fflush(stdout);
	    }

	    if(fabs(delta) < 1e-10)
	      break;

	    A[c] += delta;
	    if(REFDEBUG_STRICT < 2){
	      if(c < colors)
		A[c] = min(A[c],Amax);
	      A[c] = max(A[c],Amin);
	    }
	    if(QUADRATIC_VARIANCE){
	      if(DEBUG) assert(R[c] >= 0.0);
	      if(B[c] < 0.0 && R[c] >= 0.0)
		A[c] = max(A[c], B[c]*B[c]/(4.0*R[c]));
	    }
	  }/* Aiter = 0..2 */

	  A[c] = bestA;

	  if(REFDEBUG>=2 && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	    SF[c] = sqrt(A[c]);
	    SD[c] = copysign(sqrt(fabs(B[c])),B[c]);
	    if(c < colors){
	      SR[c] = sqrt(R[c]);
	      SE[c] = sqrt(E[c]);
	    }

	    if(DEBUG) assert(logLRarray);
	    score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	    double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	    double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	    lcnt= 0;
	    double *origLRarray = new double[nummaps];
	    double *origSDarray[MAXCOLOR];
	    for(int mid = 0; mid < nummaps; mid++)
	      origLRarray[mid] = logLRarray[mid];
	    for(int lc = 0; lc < colors; lc++){
	      origSDarray[lc] = new double[nummaps];
	      for(int mid = 0; mid < nummaps; mid++){
		origSDarray[lc][mid] = logSDarray[lc][mid];
	      }
	    }

	    double origLRsum = 0.0;
	    double origSDsum[MAXCOLOR] = {0.0,0.0};

	    for(size_t i = 0; i < numaligns; i++){
	      double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	      Calign *align = alignment[i];
	      if(!align || align->numpairs <= 1)
		continue;
	      int mid = align->mapid2;
	      if(DEBUG) assert(!Gmap[mid]->origmap);
	      if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
		if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
		continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	      }
	      if(align->score <= origScoreThreshold)
		continue;
	      int rid = align->mapid1;
	      if(DEBUG>=2){
		int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
		assert(origgoodalign);
		int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
		assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	      }
	      FLOAT **XX = Gmap[mid]->site;
	      int *MM = Gmap[mid]->numsite;
	      FLOAT **YY = refmap[rid]->site;
	      //	int *NN = Gmap[rid]->numsite;
	      align++;
	      logLRarray[mid] = 0.0;
	      for(int c = 0; c < colors; c++){
		FLOAT *X = XX[c];
		double M = MM[c];
		FLOAT *Y = YY[c];
		int N = refmap[rid]->numsite[c];
		int U = align[c].numpairs;

		LR[c] = logLR(mid,rid,c,X,Y,M,N,align,align[c].sites1[0], align[c].sites2[0], align[c].sitesK1[0], align[c].sites2[U-1], align[c].Lij1,align[c].Rij1,align[c].Lij2,align[c].Rij2,1,0,SDDEBUG ? &LRsd[c] : 0,NULL);
		logLRarrayC[c][mid] = LR[c];
		logSDarray[c][mid] = LRsd[c];
		if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
		  continue;
		LRsum += LR[c];
		if(SDDEBUG){
		  LRsdsum += LRsd[c];
		  LRsdsumC[c] += LRsd[c];
		}
	      }

	      int U = align[0].numpairs;
	      int I = align[0].sites1[U-1];
	      int K = align[0].sitesK1[U-1];
	      int J = align[0].sites2[U-1];
	      U = align[1].numpairs;
	      int P = align[1].sites1[U-1];
	      int T = align[1].sitesK1[U-1];
	      int Q = align[1].sites2[U-1];
	      double y0 = Yc(YY[0],I,K);
	      double y1 = Yc(YY[1],P,T);
	      double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	      double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	      double len = (COLOR_SHIFT==2) ? max(y0,y1) : max(x0,x1);
	      
	      LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
	      if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		LRsum += LR[2];
		for(int lc = 0; lc <= colors; lc++)
		  LRsumC[lc] += LR[lc];
	      }

	      logLRarray[mid] = LR[0] + LR[1] + LR[2];
	      lcnt++;
	      if(VERB>=2 && giter==0){
		if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		  origLRsum += origLRarray[mid];
		  for(int lc = 0; lc < colors; lc++)
		    origSDsum[lc] += origSDarray[lc][mid];
		}
		if(mid==265){
		  printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)(delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f)\n",
			 i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
			 logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
		  printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
			 logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
			 LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
			 LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
		  printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);
		}
	      }
	    }

	    if(SDDEBUG)
	      printf("iter=%d:After updating SF[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		     iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt, LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	    else
	      printf("iter=%d:After updating SF[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",
		     iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt, LRsumC[1]/lcnt, LRsumC[2]/lcnt);
	    fflush(stdout);

	    if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	      if(VERB){
		double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
		for(size_t i = 0; i < numaligns; i++){
		  Calign *align = alignment[i];
		  if(!align || align->numpairs <= 1)
		    continue;
		  int mid = align->mapid2;
		  if(DEBUG>=2) assert(!Gmap[mid]->origmap);
		  if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		    continue;
		  if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		    continue;
		  sum1 += origLRarray[mid];
		  sum2 += logLRarray[mid];
		  sum3 += origSDarray[0][mid] + origSDarray[1][mid];
		  sum4 += logSDarray[0][mid] + logSDarray[1][mid];
		  printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
			 i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
			 origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
			 logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
		}
		fflush(stdout);
	      }
	      
	      assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	    }

	    if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	      printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	      assert(lcnt == tcnt);
	    }
	    logLRsum = LRsum;
	    if(SDDEBUG)
	      logLRsdsum = LRsdsum;
	    score_free(0,numrefmaps);

	    delete [] origLRarray;
	    for(int c = 0; c < colors; c++){
	      delete [] origSDarray[c];
	    }
	  } // if(REFDEBUG>=2 ... )

	  if(QUADRATIC_VARIANCE && MAXSR > MINSR && c < colors){/* update estimate of R[c] (separately from A[c], B[c] & E[c] ( 3 iterations with delta = -f'/f'') */
	    bestLL = MINSCORE;
	    double bestR = R[c];
	    for(int Riter = 0; Riter < 3; Riter++){
	      double Lsum = 0.0;/* sum of log(A+By+Ryy+E*resvar) */
	      double Isum = 0.0;/* sum of yy/(A+By+Ryy+E*resvar) */
	      double Isum2 = 0.0;/* sum of yyyy/(A+By+Ryy+E*resvar)^2 */
	      double errsum = 0.0;/* sum of (Cx-y)^2/(A+By+Ryy+E*resvar) */
	      double errsum2 = 0.0;/* sum of yy (Cx-y)^2/(A+By+Ryy+E*resvar)^2 */
	      double errsum3 = 0.0;/* sum of yyyy (Cx-y)^2/(A+By+Ryy+E*resvar)^3 */
	      Cinterval *perr = errors[c];

	      for(int tid = 0; tid < numthreads; tid++)
		Tarray1[tid] = Tarray2[tid] = Tarray3[tid] = Tarray4[tid] = Tarray5[tid] = Tarray6[tid] = 0.0;

	      // NOTE : non-deterministic with gcc 4.6, even if Tarray*[] are sorted by value, unless -fno-unsafe-math-optimizations is used
              #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
	      {
		int tid = 0;
                #ifdef _OPENMP
		tid = omp_get_thread_num ();
                #endif

		double myLsum = 0.0;
		double myIsum = 0.0;
		double myIsum2 = 0.0;
		double myerrsum = 0.0;
		double myerrsum2 = 0.0;
		double myerrsum3 = 0.0;

                #pragma omp for schedule(static,64)
		for(int i = 0; i < numerrors[c]; i++){
		  double y = perr[i].y;
		  double x = C*perr[i].x;
		  if(ENDFIX>=3 && perr[i].end && y >= x)
		    continue;
		  double err = y - x;
		  err *= err;

		  double var = A[c] + B[c] * y;
		  if(QUADRATIC_VARIANCE)
		    var += R[c] * y*y;
		  if(RES_VARIANCE)
		    var += E[c] * perr[i].resvar;
		  double Ivar = 1.0/var;
		  double yIvar = y*y*Ivar;
		  double yIvar2 = yIvar*yIvar;

		  if(!(ENDFIX>=2 && perr[i].end)){
		    myLsum += log(var);
		    if(SDDEBUG) myLsum -= logXtheta[c];
		    myIsum += yIvar;
		    myIsum2 += yIvar2;
		  }

		  myerrsum += err *= Ivar;
		  myerrsum2 += err * yIvar;
		  myerrsum3 += err * yIvar2;
		}
	    
		Tarray1[tid] = myLsum;
		Tarray2[tid] = myIsum;
		Tarray3[tid] = myIsum2;
		Tarray4[tid] = myerrsum;
		Tarray5[tid] = myerrsum2;
		Tarray6[tid] = myerrsum3;
	      }
	      qsort(Tarray1, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray2, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray3, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray4, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray5, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray6, numthreads, sizeof(double), (intcmp *)DoubleInc);

	      for(int tid = 0; tid < numthreads; tid++){
		Lsum += Tarray1[tid];
		Isum += Tarray2[tid];
		Isum2 += Tarray3[tid];
		errsum += Tarray4[tid];
		errsum2 += Tarray5[tid];
		errsum3 += Tarray6[tid];
	      }

	      double LL = -(Lsum+errsum)*0.5;
	      double LL1 = 0.5*(errsum2-Isum);
	      double LL2 = 0.5*Isum2-errsum3;
	      if(LL < bestLL){
		if(VERB && (VERB >= 2 || SDDEBUG)){
		  printf("iter=%d,%d,c=%d:sf=%0.8f,sd=%0.8f,sr=%0.8f,se=%0.8f:LL=%0.8f,LL1=%0.8f,LL2=%0.8f:backtracking (sr->%0.8f)\n",
			 iter,Riter,c,sqrt(A[c]),copysign(sqrt(fabs(B[c])),B[c]),sqrt(R[c]), sqrt(E[c]),LL, LL1, LL2, sqrt(0.5*(R[c]+bestR)));
		  fflush(stdout);
		}
		R[c] = (R[c] + bestR)*0.5;
		if(LL < bestLL - 1e-6 && fabs(R[c] - bestR) > 1e-6)
		  Riter--;/* add an iteration of drop is LL was significant */
		continue;
	      }
	      bestR = R[c];
	      bestLL = LL;

	      double maxdelta = Rmax * 0.5;
	      double delta = (LL2 < 0.0) ? min(maxdelta,fabs(LL1/LL2)) : maxdelta;
	      delta = copysign(delta,LL1);
	      if(R[c] + delta < Rmin)
		delta = Rmin - R[c];
	      if(VERB && (VERB >= 2 || SDDEBUG)){
		printf("iter=%d,%d,c=%d:sf=%0.8f,sd=%0.8f,sr=%0.8f,se=%0.8f:LL=%0.8f,LL1=%0.8f,LL2=%0.8f:delta=%0.8f(sr->%0.8f)\n",
		       iter,Riter,c,sqrt(A[c]),copysign(sqrt(fabs(B[c])),B[c]),sqrt(R[c]),sqrt(E[c]),LL, LL1, LL2, delta,sqrt(R[c]+delta));
		fflush(stdout);
	      }

	      R[c] += delta;
	      if(REFDEBUG_STRICT < 2){
		R[c] = min(R[c],Rmax);
		R[c] = max(R[c],Rmin);
	      }
	      if(DEBUG) assert(A[c] >= 0.0);
	      if(B[c] < 0.0 && A[c] >= 0.0)
		R[c] = max(R[c], B[c]*B[c]/(4.0*A[c]));
	    } /* Riter = 0..2 */
	    R[c] = bestR;

	    if(REFDEBUG>=2 && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	      SF[c] = sqrt(A[c]);
	      SD[c] = copysign(sqrt(fabs(B[c])),B[c]);
	      if(c < colors){
		SR[c] = sqrt(R[c]);
		SE[c] = sqrt(E[c]);
	      }

	      if(DEBUG) assert(logLRarray);
	      score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	      double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	      double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	      lcnt= 0;
	      double *origLRarray = new double[nummaps];
	      double *origSDarray[MAXCOLOR];
	      for(int mid = 0; mid < nummaps; mid++)
		origLRarray[mid] = logLRarray[mid];
	      for(int lc = 0; lc < colors; lc++){
		origSDarray[lc] = new double[nummaps];
		for(int mid = 0; mid < nummaps; mid++){
		  origSDarray[lc][mid] = logSDarray[lc][mid];
		}
	      }

	      double origLRsum = 0.0;
	      double origSDsum[MAXCOLOR] = {0.0,0.0};

	      for(size_t i = 0; i < numaligns; i++){
		double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
		Calign *align = alignment[i];
		if(!align || align->numpairs <= 1)
		  continue;
		int mid = align->mapid2;
		if(DEBUG) assert(!Gmap[mid]->origmap);
		if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
		  if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
		  continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
		}
		if(align->score <= origScoreThreshold)
		  continue;
		int rid = align->mapid1;
		if(DEBUG>=2){
		  int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
		  assert(origgoodalign);
		  int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
		  assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
		}
		FLOAT **XX = Gmap[mid]->site;
		int *MM = Gmap[mid]->numsite;
		FLOAT **YY = refmap[rid]->site;
		//	int *NN = Gmap[rid]->numsite;
		align++;
		logLRarray[mid] = 0.0;
		for(int c = 0; c < colors; c++){
		  FLOAT *X = XX[c];
		  double M = MM[c];
		  FLOAT *Y = YY[c];
		  int N = refmap[rid]->numsite[c];
		  int U = align[c].numpairs;

		  LR[c] = logLR(mid,rid,c,X,Y,M,N,align,align[c].sites1[0], align[c].sites2[0], align[c].sitesK1[0], align[c].sites2[U-1], align[c].Lij1,align[c].Rij1,align[c].Lij2,align[c].Rij2,1,0,SDDEBUG ? &LRsd[c] : 0,NULL);
		  logLRarrayC[c][mid] = LR[c];
		  logSDarray[c][mid] = LRsd[c];
		  if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
		    continue;
		  LRsum += LR[c];
		  if(SDDEBUG){
		    LRsdsum += LRsd[c];
		    LRsdsumC[c] += LRsd[c];
		  }
		}

		int U = align[0].numpairs;
		int I = align[0].sites1[U-1];
		int K = align[0].sitesK1[U-1];
		int J = align[0].sites2[U-1];
		U = align[1].numpairs;
		int P = align[1].sites1[U-1];
		int T = align[1].sitesK1[U-1];
		int Q = align[1].sites2[U-1];
		double y0 = Yc(YY[0],I,K);
		double y1 = Yc(YY[1],P,T);
		double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
		double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
		double len = (COLOR_SHIFT==2) ? max(y0,y1) : max(x0,x1);

		LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
		if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		  LRsum += LR[2];
		  for(int lc = 0; lc <= colors; lc++)
		    LRsumC[lc] += LR[lc];
		}

		logLRarray[mid] = LR[0] + LR[1] + LR[2];
		lcnt++;
		if(VERB>=2 && giter==0){
		  if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		    origLRsum += origLRarray[mid];
		    for(int lc = 0; lc < colors; lc++)
		      origSDsum[lc] += origSDarray[lc][mid];
		  }
		  if(mid==265){
		    printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)(delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f)\n",
			   i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
			   logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
		    printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
			   logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
			   LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
			   LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
		    printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);
		  }
		}
	      }
	      if(SDDEBUG)
		printf("iter=%d:After updating SR[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		       iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt, LRsumC[2]/lcnt, logLRsdsum/lcnt,LRsdsum/lcnt, LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	      else
		printf("iter=%d:After updating SR[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",
		       iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt, LRsumC[2]/lcnt);
	      fflush(stdout);

	      if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
		if(VERB){
		  double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
		  for(size_t i = 0; i < numaligns; i++){
		    Calign *align = alignment[i];
		    if(!align || align->numpairs <= 1)
		      continue;
		    int mid = align->mapid2;
		    if(DEBUG>=2) assert(!Gmap[mid]->origmap);
		    if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		      continue;
		    if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		      continue;
		    sum1 += origLRarray[mid];
		    sum2 += logLRarray[mid];
		    sum3 += origSDarray[0][mid] + origSDarray[1][mid];
		    sum4 += logSDarray[0][mid] + logSDarray[1][mid];
		    printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
			   i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
			   origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
			   logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
		  }
		  fflush(stdout);
		}
	      
		assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	      }

	      if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
		printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
		assert(lcnt == tcnt);
	      }
	      logLRsum = LRsum;
	      if(SDDEBUG)
		logLRsdsum = LRsdsum;
	      score_free(0,numrefmaps);

	      delete [] origLRarray;
	      for(int c = 0; c < colors; c++){
		delete [] origSDarray[c];
	      }
	    } // if(REFDEBUG>=2 ... )
	  }// if(QUADRATIC_VARIANCE && MAXSR > MINSR && c < colors)

	  if(c < colors){/* update estimate of B[c] (separately from A[c], R[c] & E[c]) (3 iterations with delta = -f'/f'') */

	    if(VERB>=2 && SDDEBUG){
	      printf("giter=%d,iter=%d,c=%d:Estimating scaling factor for SD[c]=%0.8f\n",giter,iter,c,SD[c]);
	      printf("  MA[2]=%0.8f,MA_mean=%0.8f,SF[2]*SF[2]*2 = %0.8f,FM2var = %0.8f, fabs(SD[2])*SD[2]*2= %0.8f, SM2var= %0.8f, 0.5*logXtheta[2]=%0.8f, RefPen12 - 0.5*log(2)=%0.8f, log(Xtheta12) - 0.5*log(2*PI)=%0.8f\n",
		     MA[2],MA_mean,SF[2]*SF[2]*2.0, FM2var, fabs(SD[2])*SD[2]*2.0, SM2var, 0.5*logXtheta[2], RefPen12 - 0.5*log(2.0), log(Xtheta12) - 0.5*log(2.0*M_PI));
	      fflush(stdout);
	    }

	    bestLL = MINSCORE;
	    double bestB = B[c];
	    for(int Biter = 0; Biter < 3; Biter++){
	      double Lsum = 0.0;/* sum of log(A+By+Ryy+E*resvar) */
	      double Isum = 0.0;/* sum of y/(A+By+Ryy+E*resvar) */
	      double Isum2 = 0.0;/* sum of yy/(A+By+Ryy+E*resvar)^2 */
	      double errsum = 0.0;/* sum of (Cx-y)^2/(A+By+Ryy+E*resvar) */
	      double errsum2 = 0.0;/* sum of y (Cx-y)^2/(A+By+Ryy+E*resvar)^2 */
	      double errsum3 = 0.0;/* sum of yy (Cx-y)^2/(A+By+Ryy+E*resvar)^3 */
	      Cinterval *perr = errors[c];

	      for(int tid = 0; tid < numthreads; tid++)
		Tarray1[tid] = Tarray2[tid] = Tarray3[tid] = Tarray4[tid] = Tarray5[tid] = Tarray6[tid] = 0.0;

	      // NOTE : non-deterministic with gcc 4.6, even if Tarray*[] are sorted by value, unless -fno-unsafe-math-optimizations is used
	      #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
	      {
		int tid = 0;
                #ifdef _OPENMP
		tid = omp_get_thread_num ();
                #endif

		double myLsum = 0.0;
		double myIsum = 0.0;
		double myIsum2 = 0.0;
		double myerrsum = 0.0;
		double myerrsum2 = 0.0;
		double myerrsum3 = 0.0;
		
                #pragma omp for schedule(static,64)
		for(int i = 0; i < numerrors[c]; i++){
		  double y = perr[i].y;
		  double x = C*perr[i].x;
		  if(ENDFIX>=3 && perr[i].end && y >= x)
		    continue;
		  double err = y - x;
		  err *= err;

		  double var = A[c] + B[c] * y;
		  if(QUADRATIC_VARIANCE)
		    var += R[c] * y*y;
		  if(RES_VARIANCE)
		    var += E[c] * perr[i].resvar;
		  double Ivar = 1.0/var;
		  double yIvar = y*Ivar;
		  double yIvar2 = yIvar*yIvar;

		  if(!(ENDFIX>=2 && perr[i].end)){
		    myLsum += log(var);
		    myIsum += yIvar;
		    myIsum2 += yIvar2;
		  }

		  myerrsum += err *= Ivar;
		  myerrsum2 += err * yIvar;
		  myerrsum3 += err * yIvar2;
		}
	    
		Tarray1[tid] = myLsum;
		Tarray2[tid] = myIsum;
		Tarray3[tid] = myIsum2;
		Tarray4[tid] = myerrsum;
		Tarray5[tid] = myerrsum2;
		Tarray6[tid] = myerrsum3;
	      }// parallel
	      qsort(Tarray1, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray2, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray3, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray4, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray5, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray6, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      
	      for(int tid = 0; tid < numthreads; tid++){
		Lsum += Tarray1[tid];
		Isum += Tarray2[tid];
		Isum2 += Tarray3[tid];
		errsum += Tarray4[tid];
		errsum2 += Tarray5[tid];
		errsum3 += Tarray6[tid];
	      }

	      double LL = -(Lsum+errsum)*0.5;
	      double LL1 = 0.5*(errsum2-Isum);
	      double LL2 = 0.5*Isum2-errsum3;
	      if(LL < bestLL){
		if(VERB && (VERB >= 2 || SDDEBUG)){
		  printf("iter=%d,%d,c=%d:sf=%0.8f,sd=%0.8f,sr=%0.8f,se=%0.8f:LL=%0.8f,LL1=%0.8f,LL2=%0.8f:backtracking (sd->%0.8f)\n",
			 iter,Biter,c,sqrt(A[c]),copysign(sqrt(fabs(B[c])),B[c]),sqrt(R[c]), sqrt(E[c]),LL, LL1, LL2, copysign(sqrt(0.5*fabs(B[c]+bestB)),B[c]+bestB));
		  fflush(stdout);
		}
		B[c] = (B[c] + bestB)*0.5;
		if(LL < bestLL - 1e-6 && fabs(B[c] - bestB) > 1e-6)
		  Biter--;/* add an iteration of drop is LL was significant */
		continue;
	      }
	      bestB = B[c];
	      bestLL = LL;

	      double maxdelta = Bmax * 0.5;
	      double delta = (LL2 < 0.0) ? min(maxdelta,fabs(LL1/LL2)) : maxdelta;
	      delta = copysign(delta,LL1);
	      
	      if(DEBUG) assert(R[c] >= 0.0 && A[c] >= 0.0);
	      if(B[c] + delta < -2.0*sqrt(A[c]*R[c]) * MINSD_MULT * MINSD_MULT)
		delta = -2.0*sqrt(A[c]*R[c])*MINSD_MULT*MINSD_MULT - B[c];
	      if(VERB && (VERB >= 2 || SDDEBUG)){
		printf("iter=%d,%d,c=%d:sf=%0.8f,sd=%0.8f,sr=%0.8f,se=%0.8f:LL=%0.8f,LL1=%0.8f,LL2=%0.8f:delta=%0.8f(sd->%0.8f)\n",
		       iter,Biter,c,sqrt(A[c]),copysign(sqrt(fabs(B[c])),B[c]),sqrt(R[c]),sqrt(E[c]),LL, LL1, LL2, delta,copysign(sqrt(fabs(B[c]+delta)),B[c]+delta));
		fflush(stdout);
	      }
	      B[c] += delta;
	      if(REFDEBUG_STRICT < 2){
		B[c] = min(B[c],Bmax);
		B[c] = max(B[c],Bmin);
	      }
	    } /* Biter = 0..2 */

	    B[c] = bestB;

	    if(REFDEBUG>=2 && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	      SF[c] = sqrt(A[c]);
	      SD[c] = copysign(sqrt(fabs(B[c])),B[c]);
	      if(c < colors){
		SR[c] = sqrt(R[c]);
		SE[c] = sqrt(E[c]);
	      }

	      if(DEBUG) assert(logLRarray);
	      score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	      double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	      double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	      lcnt= 0;
	      double *origLRarray = new double[nummaps];
	      double *origSDarray[MAXCOLOR];
	      for(int mid = 0; mid < nummaps; mid++)
		origLRarray[mid] = logLRarray[mid];
	      for(int lc = 0; lc < colors; lc++){
		origSDarray[lc] = new double[nummaps];
		for(int mid = 0; mid < nummaps; mid++){
		  origSDarray[lc][mid] = logSDarray[lc][mid];
		}
	      }

	      double origLRsum = 0.0;
	      double origSDsum[MAXCOLOR] = {0.0,0.0};

	      for(size_t i = 0; i < numaligns; i++){
		double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
		Calign *align = alignment[i];
		if(!align || align->numpairs <= 1)
		  continue;
		int mid = align->mapid2;
		if(DEBUG) assert(!Gmap[mid]->origmap);
		if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
		  if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
		  continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
		}
		if(align->score <= origScoreThreshold)
		  continue;
		int rid = align->mapid1;
		if(DEBUG>=2){
		  int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
		  assert(origgoodalign);
		  int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
		  assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
		}
		FLOAT **XX = Gmap[mid]->site;
		int *MM = Gmap[mid]->numsite;
		FLOAT **YY = refmap[rid]->site;
		//	int *NN = Gmap[rid]->numsite;
		align++;
		logLRarray[mid] = 0.0;
		for(int lc = 0; lc < colors; lc++){
		  FLOAT *X = XX[lc];
		  double M = MM[lc];
		  FLOAT *Y = YY[lc];
		  int N = refmap[rid]->numsite[lc];
		  int U = align[lc].numpairs;

		  LR[lc] = logLR(mid,rid,lc,X,Y,M,N,align,align[lc].sites1[0], align[lc].sites2[0], align[lc].sitesK1[0], align[lc].sites2[U-1], align[lc].Lij1,align[lc].Rij1,align[lc].Lij2,align[lc].Rij2,1,0,SDDEBUG ? &LRsd[lc] : 0,NULL);
		  logLRarrayC[lc][mid] = LR[lc];
		  logSDarray[lc][mid] = LRsd[lc];
		  if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
		    continue;
		  LRsum += LR[lc];
		  if(SDDEBUG){
		    LRsdsum += LRsd[lc];
		    LRsdsumC[lc] = LRsd[lc];
		  }
		}
		int U = align[0].numpairs;
		int I = align[0].sites1[U-1];
		int K = align[0].sitesK1[U-1];
		int J = align[0].sites2[U-1];
		U = align[1].numpairs;
		int P = align[1].sites1[U-1];
		int T = align[1].sitesK1[U-1];
		int Q = align[1].sites2[U-1];
		double y0 = Yc(YY[0],I,K);
		double y1 = Yc(YY[1],P,T);
		double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
		double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
		double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);

		LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
		if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		  LRsum += LR[2];
		  for(int lc = 0; lc <= colors; lc++)
		    LRsumC[lc] += LR[lc];
		}

		logLRarray[mid] = LR[0] + LR[1] + LR[2];
		lcnt++;
		if(VERB>=2 && giter==0 && iter==1 && c==1){
		  if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		    origLRsum += origLRarray[mid];
		    for(int lc = 0; lc < colors; lc++)
		      origSDsum[lc] += origSDarray[lc][mid];
		  }
		  if(1){
		    printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)(delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f)\n",
			   i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
			   logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
		    printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
			   logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
			   LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
			   LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
		    printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);
		  }
		}
	      }

	      if(SDDEBUG)
		printf("iter=%d:After updating SD[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		       iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt,LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	      else
		printf("iter=%d:After updating SD[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);

	      fflush(stdout);

	      if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
		if(VERB){
		  double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
		  for(size_t i = 0; i < numaligns; i++){
		    Calign *align = alignment[i];
		    if(!align || align->numpairs <= 1)
		      continue;
		    int mid = align->mapid2;
		    if(DEBUG>=2) assert(!Gmap[mid]->origmap);
		    if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		      continue;
		    if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		      continue;
		    sum1 += origLRarray[mid];
		    sum2 += logLRarray[mid];
		    sum3 += origSDarray[0][mid] + origSDarray[1][mid];
		    sum4 += logSDarray[0][mid] + logSDarray[1][mid];
		    printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
			   i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
			   origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
			   logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
		  }
		  fflush(stdout);
		}
	      
		assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	      }

	      if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
		printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
		assert(lcnt == tcnt);
	      }
	      logLRsum = LRsum;
	      if(SDDEBUG)
		logLRsdsum = LRsdsum;
	      score_free(0,numrefmaps);

	      delete [] origLRarray;
	      for(int c = 0; c < colors; c++){
		delete [] origSDarray[c];
	      }
	    } // if(REFDEBUG>=2 ... )

	  }/* if c < colors */

	  if(RES_VARIANCE && MAXSE > MINSE && c < colors){
	    bestLL = MINSCORE;
	    double bestE = E[c];

	    /* update estimate E[c] (seperately from A[c], B[c], R[c]) (3 iterations with delta = -f'/f'') */
	    for(int Eiter = 0; Eiter < 3; Eiter++){
	      double Lsum = 0.0;/* sum of log(A+By+Ryy+E*resvar) */
	      double Isum = 0.0;/* sum of resvar/(A+By+Ryy+E*resvar) */
	      double Isum2 = 0.0;/* sum of resvar^2/(A+By+Ryy+E*resvar)^2 */
	      double errsum = 0.0;/* sum of (Cx-y)^2/(A+By+Ryy+E*resvar) */
	      double errsum2 = 0.0;/* sum of resvar (Cx-y)^2/(A+By+Ryy+E*resvar)^2 */
	      double errsum3 = 0.0;/* sum of resvar^2 (Cx-y)^2/(A+By+Ryy+E*resvar)^3 */
	      Cinterval *perr = errors[c];

	      for(int tid = 0; tid < numthreads; tid++)
		Tarray1[tid] = Tarray2[tid] = Tarray3[tid] = Tarray4[tid] = Tarray5[tid] = Tarray6[tid] = 0.0;

	      // NOTE : non-deterministic with gcc 4.6, even if Tarray*[] are sorted by value, unless -fno-unsafe-math-optimizations is used
              #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
	      {
		int tid = 0;
                #ifdef _OPENMP
		tid = omp_get_thread_num ();
                #endif

		double myLsum = 0.0;
		double myIsum = 0.0;
		double myIsum2 = 0.0;
		double myerrsum = 0.0;
		double myerrsum2 = 0.0;
		double myerrsum3 = 0.0;
		
                #pragma omp for schedule(static,64)
		for(int i = 0; i < numerrors[c]; i++){
		  double y = perr[i].y;
		  double x = C*perr[i].x;
		  if(ENDFIX>=3 && perr[i].end && y >= x)
		    continue;
		  double err = y - x;
		  err *= err;
		  
		  double var = A[c] + B[c] * y;
		  if(QUADRATIC_VARIANCE)
		    var += R[c] * y*y;
		  double resvar = perr[i].resvar;
		  var += E[c] * resvar;
		  double Ivar = 1.0/var;
		  double yIvar = resvar*Ivar;
		  double yIvar2 = yIvar*yIvar;
		
		  if(!(ENDFIX>=2 && perr[i].end)){
		    myLsum += log(var);
		    myIsum += yIvar;
		    myIsum2 += yIvar2;
		  }

		  myerrsum += err *= Ivar;
		  myerrsum2 += err * yIvar;
		  myerrsum3 += err * yIvar2;
		}
	    
		Tarray1[tid] = myLsum;
		Tarray2[tid] = myIsum;
		Tarray3[tid] = myIsum2;
		Tarray4[tid] = myerrsum;
		Tarray5[tid] = myerrsum2;
		Tarray6[tid] = myerrsum3;
	      }// parallel
	      qsort(Tarray1, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray2, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray3, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray4, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray5, numthreads, sizeof(double), (intcmp *)DoubleInc);
	      qsort(Tarray6, numthreads, sizeof(double), (intcmp *)DoubleInc);

	      for(int tid = 0; tid < numthreads; tid++){
		Lsum += Tarray1[tid];
		Isum += Tarray2[tid];
		Isum2 += Tarray3[tid];
		errsum += Tarray4[tid];
		errsum2 += Tarray5[tid];
		errsum3 += Tarray6[tid];
	      }

	      double LL = -(Lsum+errsum)*0.5;
	      double LL1 = 0.5*(errsum2-Isum);
	      double LL2 = 0.5*Isum2-errsum3;
	      if(LL < bestLL){
		if(VERB && (VERB >= 2 || SDDEBUG)){
		  printf("iter=%d,%d,c=%d:sf=%0.8f,sd=%0.8f,sr=%0.8f,se=%0.8f:LL=%0.8f,LL1=%0.8f,LL2=%0.8f:backtracking (se->%0.8f)\n",
			 iter,Eiter,c,sqrt(A[c]),copysign(sqrt(fabs(B[c])),B[c]),sqrt(R[c]), sqrt(E[c]),LL, LL1, LL2, sqrt((E[c] + bestE)*0.5));
		  fflush(stdout);
		}
		E[c] = (E[c]+bestE)*0.5;
		if(LL < bestLL - 1e-6 && fabs(E[c]-bestE) > 1e-6)
		  Eiter--;/* add an iteration if drop in LL was significant */
		continue;
	      }
	      bestE = E[c];
	      bestLL = LL;

	      double maxdelta = Emax * 0.5;
	      double delta = (LL2 < 0.0) ? min(maxdelta,fabs(LL1/LL2)) : maxdelta;
	      delta = copysign(delta,LL1);
	      if(E[c] + delta < Emin)
		delta = Emin - E[c];

	      if(VERB && (VERB>=2 || SDDEBUG)){
		printf("iter=%d,%d,c=%d:sf=%0.8f,sd=%0.8f,sr=%0.8f,se=%0.8f:LL=%0.8f,LL1=%0.8f,LL2=%0.8f:delta=%0.8f(se->%0.8f)\n",
		       iter,Eiter,c,sqrt(A[c]),copysign(sqrt(fabs(B[c])),B[c]),sqrt(R[c]),sqrt(E[c]),LL, LL1, LL2, delta, sqrt(E[c] + delta));
		fflush(stdout);
	      }
	      
	      E[c] += delta;
	      if(DEBUG) assert(E[c] >= 0.0 && isfinite(E[c]));
	      if(REFDEBUG_STRICT < 2){
		E[c] = min(E[c],Emax);
		E[c] = max(E[c],Emin);
	      }
	    } // Eiter = 0 ... 2
	    E[c] = bestE;

	    if(REFDEBUG>=2 && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	      SF[c] = sqrt(A[c]);
	      SD[c] = copysign(sqrt(fabs(B[c])),B[c]);
	      if(c < colors){
		SR[c] = sqrt(R[c]);
		SE[c] = sqrt(E[c]);
	      }

	      if(DEBUG) assert(logLRarray);
	      score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	      double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	      double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	      lcnt= 0;
	      double *origLRarray = new double[nummaps];
	      double *origSDarray[MAXCOLOR];
	      for(int mid = 0; mid < nummaps; mid++)
		origLRarray[mid] = logLRarray[mid];
	      for(int lc = 0; lc < colors; lc++){
		origSDarray[lc] = new double[nummaps];
		for(int mid = 0; mid < nummaps; mid++){
		  origSDarray[lc][mid] = logSDarray[lc][mid];
		}
	      }

	      double origLRsum = 0.0;
	      double origSDsum[MAXCOLOR] = {0.0,0.0};

	      for(size_t i = 0; i < numaligns; i++){
		double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
		Calign *align = alignment[i];
		if(!align || align->numpairs <= 1)
		  continue;
		int mid = align->mapid2;
		if(DEBUG) assert(!Gmap[mid]->origmap);
		if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
		  if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
		  continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
		}
		if(align->score <= origScoreThreshold)
		  continue;
		int rid = align->mapid1;
		if(DEBUG>=2){
		  int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
		  assert(origgoodalign);
		  int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
		  assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
		}
		FLOAT **XX = Gmap[mid]->site;
		int *MM = Gmap[mid]->numsite;
		FLOAT **YY = refmap[rid]->site;
		//	int *NN = Gmap[rid]->numsite;
		align++;
		logLRarray[mid] = 0.0;
		for(int c = 0; c < colors; c++){
		  FLOAT *X = XX[c];
		  double M = MM[c];
		  FLOAT *Y = YY[c];
		  int N = refmap[rid]->numsite[c];
		  int U = align[c].numpairs;

		  LR[c] = logLR(mid,rid,c,X,Y,M,N,align,align[c].sites1[0], align[c].sites2[0], align[c].sitesK1[0], align[c].sites2[U-1], align[c].Lij1,align[c].Rij1,align[c].Lij2,align[c].Rij2,1,0,SDDEBUG ? &LRsd[c] : 0,NULL);
		  logLRarrayC[c][mid] = LR[c];
		  logSDarray[c][mid] = LRsd[c];
		  if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
		    continue;
		  LRsum += LR[c];
		  if(SDDEBUG){
		    LRsdsum += LRsd[c];
		    LRsdsumC[c] += LRsd[c];
		  }
		}

		int U = align[0].numpairs;
		int I = align[0].sites1[U-1];
		int K = align[0].sitesK1[U-1];
		int J = align[0].sites2[U-1];
		U = align[1].numpairs;
		int P = align[1].sites1[U-1];
		int T = align[1].sitesK1[U-1];
		int Q = align[1].sites2[U-1];
		double y0 = Yc(YY[0],I,K);
		double y1 = Yc(YY[1],P,T);
		double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
		double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
		double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);
		
		LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
		if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		  LRsum += LR[2];
		  for(int lc = 0; lc <= colors; lc++)
		    LRsumC[lc] += LR[lc];
		}

		logLRarray[mid] = LR[0] + LR[1] + LR[2];
		lcnt++;
		if(VERB>=2 && giter==0){
		  if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		    origLRsum += origLRarray[mid];
		    for(int lc = 0; lc < colors; lc++)
		      origSDsum[lc] += origSDarray[lc][mid];
		  }
		  if(mid==265){
		    printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)(delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f)\n",
			   i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
			   logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
		    printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
			   logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
			   LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
			   LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
		    printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);
		  }
		}
	      }
	      if(SDDEBUG)
		printf("iter=%d:After updating SE[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		       iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt,LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	      else
		printf("iter=%d:After updating SE[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",iter,c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
	      fflush(stdout);

	      if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
		if(VERB){
		  double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
		  for(size_t i = 0; i < numaligns; i++){
		    Calign *align = alignment[i];
		    if(!align || align->numpairs <= 1)
		      continue;
		    int mid = align->mapid2;
		    if(DEBUG>=2) assert(!Gmap[mid]->origmap);
		    if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		      continue;
		    if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		      continue;
		    sum1 += origLRarray[mid];
		    sum2 += logLRarray[mid];
		    sum3 += origSDarray[0][mid] + origSDarray[1][mid];
		    sum4 += logSDarray[0][mid] + logSDarray[1][mid];
		    printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
			   i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
			   origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
			   logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
		  }
		  fflush(stdout);
		}
	      
		assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	      }

	      if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
		printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
		assert(lcnt == tcnt);
	      }
	      logLRsum = LRsum;
	      if(SDDEBUG)
		logLRsdsum = LRsdsum;
	      score_free(0,numrefmaps);

	      delete [] origLRarray;
	      for(int c = 0; c < colors; c++){
		delete [] origSDarray[c];
	      }
	    } // if(REFDEBUG>=2 ... )
	  }// if(RES_VARIANCE ... && c < colors)
	}// c = 0..colors 
      }// iter = 0 .. 29

      if(VERB) {
	for(int c = 0; c < colors; c++)
	  printf("  c=%d:var(y) = sf^2 + sd^2 y +sr^2 y^2 + se^2 r^2: sf = %0.6f -> %0.6f (kb), sd = %0.6f -> %0.6f (kb^1/2), sr= %0.6f -> %0.6f, se= %0.6f -> %0.6f\n", 
		 c, SF[c], sqrt(A[c]), SD[c], copysign(sqrt(fabs(B[c])),B[c]), SR[c], sqrt(R[c]), SE[c], sqrt(E[c]) );
	printf("  c=%d:var(y) = sf^2 + sd^2 y : sf = %0.6f -> %0.6f (kb), sd = %0.6f -> %0.6f (kb^1/2), bias= %0.6f -> %0.6f (kb)\n", 
	       colors, SF[colors], sqrt(A[colors]), SD[colors], sqrt(B[colors]), MA_mean, MA[colors]);

	/* NOTE : bppSD computation ignores map splits (treating each split as a seperate map) */
	double scalesum=0.0,scalesumsq=0.0, scalewtsum = 0.0;
	for(int i = 0; i < nummaps; i++){
	  double incscale = Gmap[i]->incscale;
	  double wt = Gmap[i]->incwt;
	  if(DEBUG>=2) assert(isfinite(wt));
	  if(DEBUG>=2) assert(isfinite(incscale));
	  scalesum += incscale * wt ;
	  scalesumsq += incscale*incscale * wt;
	  scalewtsum += wt;
	  if(VERB>=2 && giter2==RefRepeats2-1 && giter == RefRepeats-1){
	    printf("m=%d:mapid=%d,id=%lld:incscale=%0.8f,wt=%0.8f:scalesum=%0.8f,scalesumsq=%0.8f,scalewtsum=%0.8f\n",
		   i,Gmap[i]->mapid,Gmap[i]->id,incscale,wt,scalesum,scalesumsq,scalewtsum);
	    fflush(stdout);
	  }
	  if(DEBUG>=2) assert(isfinite(scalewtsum));
	}	  

	double scaleSD = 0.0;
	if(scalewtsum > 0.0){
	  scalesum /= scalewtsum;
	  scalesumsq /= scalewtsum;
	  scaleSD = sqrt(scalesumsq - scalesum*scalesum);
	}
	PixelLenSD = scaleSD * PixelLen * C;
	if(C != 1.0)
	  printf("  Bases per pixel = %0.2f -> %0.2f (S.D. =%0.2f) (FP = %0.6f,%0.6f -> %0.6f,%0.6f)\n", 
		 PixelLen*1000.0, PixelLen*C*1000.0, PixelLenSD * 1000.0, FP[0], FP[1], FP[0]/C, FP[1]/C);
	else
	  printf("  Bases per pixel = %0.2f (S.D. =%0.2f)\n", 
		 PixelLen*1000.0, PixelLenSD * 1000.0);
	if(NumScaleFactor > 1)
	  for(register int i = 0; i < NumScaleFactor; i++)
	    printf("scaleID=%d,scale=%0.4f:mapcnt=%d\n",i,ScaleFactor[i],scaleIDcnt[i]);
	fflush(stdout);
      }

      for(int c = 0; c < colors; c++){
	SF[c] = sqrt(A[c]);
	if(REFDEBUG_STRICT < 2){
	  SF[c] = min(SF[c],MAXSF);
	  SF[c] = max(SF[c],MINSF);
	}

	SD[c] = copysign(sqrt(fabs(B[c])), B[c]);
	if(REFDEBUG_STRICT < 2){
	  SD[c] = min(SD[c],MAXSD);
	  SD[c] = max(SD[c],MINSD);
	}

	SR[c] = sqrt(R[c]);
	if(REFDEBUG_STRICT < 2){
	  SR[c] = min(SR[c],MAXSR);
	  SR[c] = max(SR[c],MINSR);
	}

	if(DEBUG) assert(E[c] >= 0.0);
	SE[c] = sqrt(E[c]);
	if(DEBUG) assert(isfinite(SE[c]));
	if(REFDEBUG_STRICT < 2){
	  SE[c] = min(SE[c],MAXSE);
	  SE[c] = max(SE[c],MINSE);
	}
      } 
      SF[2] = sqrt(A[2]);
      SD[2] = sqrt(B[2]);
      if(SD[2] < 0.001)
	SD[2] = 0.001;
      MA_mean = MA[2];

      if(C != 1.0 && !(REFDEBUG>=2 && (numrefmaps==1 || (BestRef && !BestRefPV)))){ /* rescale all query maps */
	if(VERB/* HERE >=2 */){
	  printf("    Rescaling all map site[] and rawsite[] values  by C=%0.12f (reflecting change in bpp)\n",C);
	  fflush(stdout);
	}

	double totX[2] = {0.0,0.0};
	if(DEBUG>=2)
	  for(int i=0; i < nummaps; i++)
	    for(int c = 0; c < colors; c++)
	      totX[c] += Gmap[i]->site[c][Gmap[i]->numsite[c]+1];

        #pragma omp parallel for num_threads(numthreads) schedule(dynamic,256) if(numthreads > 1)
	for(int i= 0; i < nummaps; i++){
	  Cmap *pmap = Gmap[i];
	  if(DEBUG && colors>=2 && !(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3)){
            #pragma omp critical
	    {
	      printf("mapid=%d(id=%lld):numsite[0]=%d,numsite[1]=%d,site[0][numsite[0]+1]= %0.6f, site[1][numsite[1]+1]= %0.6f\n",
		     i,pmap->id,pmap->numsite[0],pmap->numsite[1],pmap->site[0][pmap->numsite[0]+1],pmap->site[1][pmap->numsite[1]+1]);
	      fflush(stdout);
	      assert(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3);
	    }
	  }
	  for(int c = 0; c < colors; c++){
	    FLOAT *X = pmap->site[c];
	    int M = pmap->numsite[c];
	    for(int j= M+1; j > 0; j--)
	      X[j] *= C;
	  }
	  if(VERB>=2 && pmap->id == 8400138831001000109LL){
	    #pragma omp critical
	    {
	      printf("Rescaled sites for map id=%lld by C=%0.12f\n",pmap->id,C);
	      for(int c = 0; c < colors; c++){
		FLOAT *X = pmap->site[c];
		int M = pmap->numsite[c];
		for(int j= 1; j <= M+1; j++)
		  printf("site[%d][%d]= %0.8f\n",c,j,X[j]);
	      }
	      fflush(stdout);
	    }
	  }
	  if(DEBUG/* HERE >=2 */) assert(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3);	    
	}

	if(DEBUG>=2){
	  double newtotX[2] = {0.0,0.0};
	  for(int c = 0; c < colors; c++){
	    for(int i=0; i < nummaps; i++)
	      newtotX[c] += Gmap[i]->site[c][Gmap[i]->numsite[c]+1];
	    if(fabs(newtotX[c] - totX[c]*C) >= newtotX[c] * 0.000001){
	      printf("C=%0.8f,c=%d,totX[c]=%0.4f,newtotX[c]=%0.4f(totX[c]*C=%0.4f)\n",
		     C,c,totX[c],newtotX[c],totX[c]*C);
	      exit(1);
	    }
	  }
	}

	if(maxresbias > mres * 0.5){/* also rescale rawsite[]. NOTE : applied to all maps including -subset maps (if -ScanScaling) and split maps */

	  if(DEBUG) assert(rawsitemaps >= totalmaps);

          #pragma omp parallel for num_threads(numthreads) schedule(dynamic,256) if(numthreads > 1)
	  for(int i= 0; i < totalmaps; i++){
	    Cmap *pmap = Gmap[i];
	    if(DEBUG/* HERE >=2 */) assert(fabs(pmap->rawsite[0][pmap->numsite[0]+1] - pmap->rawsite[1][pmap->numsite[1]+1]) < 1e-3);
	    for(int c = 0; c < colors; c++){
	      FLOAT *X = pmap->rawsite[c];
	      int M = pmap->numsite[c];
	      for(int j= M+1; j > 0; j--)
		X[j] *= C;
	    }

	    if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
              #pragma omp critical
	      {
		for(int c = 0; c < colors; c++){
		  int M = pmap->numsite[c];
		  FLOAT *RawX = pmap->rawsite[c];
		  FLOAT *X = pmap->site[c];
		  printf("mapid=%d(id=%lld),c=%d,M=%d:\n",pmap->mapid,pmap->id,c,M);
		  for(int j = 1; j <= M+1; j++)
		    printf("  j=%d:X[j]=%0.10f,X[j]=%0.10f(shift=%0.10f)\n",j,RawX[j],X[j],X[j]-RawX[j]);
		  if(pmap->origsite[c]){
		    int origM = pmap->orignumsite[c];
		    double *origX = pmap->origsite[c];
		    double scale = PixelLen / 0.500;
		    printf("origM=%d:bpp/500 = %0.12f\n",origM,scale);
		    for(int j = 1; j <= origM+1; j++)
		      printf("  j=%d:origX[j]=%0.10f,origX[j]*bpp/500=%0.10f\n",j,origX[j],origX[j]*scale);
		  }
		}
		fflush(stdout);
	      }
	    }

	    if(DEBUG>=2) assert(fabs(pmap->rawsite[0][pmap->numsite[0]+1] - pmap->rawsite[1][pmap->numsite[1]+1]) < 1e-3);
	  }

	  /* also scale resbias[0][] values */
	  for(int c = 0; c < colors; c++){
	    for(int Bin = 0; Bin <= ResBins[c]; Bin++){
	      resbias[c][Bin] *= C;
	      resbiasX[c][Bin] *= C;
	    }
	  }
	  maxresbias *= C;
	  if(maxresbias <= mres * 0.5)
	    maxresbias = mres * 0.5 + 1e-6;/* Fudge : should never happen */
	}

	double InvC = 1.0/C;
	for(int i = 0; i < nummaps; i++)
	  Gmap[i]->incscale *= InvC;
	PixelLen *= C;
	for(int c = 0; c < colors; c++){
	  FP[c] *= InvC;
	  Xtheta[c] *= C;
	  res[c] *= InvC;
	  resSD[c] *= InvC;

	  if(REFDEBUG_STRICT < 2){
	    FP[c] = min(FP[c],MAXFP);
	    FP[c] = max(FP[c],MINFP);
	  }
	}
	Xtheta12 *= C;

	for(int c = 0; c < colors;c++)
	  logXtheta[c] = 2.0*log(Xtheta[c]) - log(2.0*M_PI);/* to match terms LLsd in logLR() */
	logXtheta[colors] = 2.0*log(Xtheta12) - log(2.0*M_PI);
      }

      if(REFDEBUG && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	if(DEBUG) assert(logLRarray);
	score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	lcnt= 0;
	double *origLRarray = new double[nummaps];
	double *origSDarray[MAXCOLOR];
	for(int mid = 0; mid < nummaps; mid++)
	  origLRarray[mid] = logLRarray[mid];
	for(int c = 0; c < colors; c++){
	  origSDarray[c] = new double[nummaps];
	  for(int mid = 0; mid < nummaps; mid++){
	    origSDarray[c][mid] = logSDarray[c][mid];
	  }
	}

	double origLRsum = 0.0;
	double origSDsum[MAXCOLOR] = {0.0,0.0};

	for(size_t i = 0; i < numaligns; i++){
	  double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	  Calign *align = alignment[i];
	  if(!align || align->numpairs <= 1)
	    continue;
	  int mid = align->mapid2;
	  if(DEBUG) assert(!Gmap[mid]->origmap);
	  if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
	    if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
	    continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	  }
	  if(align->score <= origScoreThreshold)
	    continue;
	  int rid = align->mapid1;
	  if(DEBUG>=2){
	    int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
	    assert(origgoodalign);
	    int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
	    assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	  }
	  FLOAT **XX = Gmap[mid]->site;
	  int *MM = Gmap[mid]->numsite;
	  FLOAT **YY = refmap[rid]->site;
	  //	int *NN = Gmap[rid]->numsite;
	  align++;
	  logLRarray[mid] = 0.0;
	  for(int c = 0; c < colors; c++){
	    FLOAT *X = XX[c];
	    double M = MM[c];
	    FLOAT *Y = YY[c];
	    int N = refmap[rid]->numsite[c];
	    int U = align[c].numpairs;

	    LR[c] = logLR(mid,rid,c,X,Y,M,N,align,align[c].sites1[0], align[c].sites2[0], align[c].sitesK1[0], align[c].sites2[U-1], align[c].Lij1,align[c].Rij1,align[c].Lij2,align[c].Rij2,1,
                                 0 /* (giter==0 && c==0 && mid == 674) ? 1 : 0 */, SDDEBUG ? &LRsd[c] : 0,NULL);
	    logLRarrayC[c][mid] = LR[c];
	    logSDarray[c][mid] = LRsd[c];
	    if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
	      continue;
	    LRsum += LR[c];
	    if(SDDEBUG){
	      LRsdsum += LRsd[c];
	      LRsdsumC[c] += LRsd[c];
	    }
	  } // c = 0 .. colors -1

	  int U = align[0].numpairs;
	  int I = align[0].sites1[U-1];
	  int K = align[0].sitesK1[U-1];
	  int J = align[0].sites2[U-1];
	  U = align[1].numpairs;
	  int P = align[1].sites1[U-1];
	  int T = align[1].sitesK1[U-1];
	  int Q = align[1].sites2[U-1];
	  double y0 = Yc(YY[0],I,K);
	  double y1 = Yc(YY[1],P,T);
	  double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	  double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	  double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);

	  LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
	  if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	    LRsum += LR[2];
	    for(int lc = 0; lc <= colors; lc++)
	      LRsumC[lc] += LR[lc];
	  }

	  logLRarray[mid] = LR[0] + LR[1] + LR[2];
	  lcnt++;
	  if(VERB>=2 && giter==0){
	    if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	      origLRsum += origLRarray[mid];
	      for(int lc = 0; lc < colors; lc++)
		origSDsum[lc] += origSDarray[lc][mid];
	    }
	    if(1){
	      printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)\n\t delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f\n",
		     i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
		     logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
	      /*	      printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
		     logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
		     LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
		     LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
		     printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);*/
	    }
	  }
	}
	if(SDDEBUG)
	  printf("After SF,SD,SR,SE,bpp updates:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		 (unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt,LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	else
	  printf("After SF,SD,SR,SE,bpp updates:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
	fflush(stdout);

	if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	  if(VERB){
	    double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
	    for(size_t i = 0; i < numaligns; i++){
	      Calign *align = alignment[i];
	      if(!align || align->numpairs <= 1)
		continue;
	      int mid = align->mapid2;
	      if(DEBUG>=2) assert(!Gmap[mid]->origmap);
	      if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		continue;
	      if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		continue;
	      sum1 += origLRarray[mid];
	      sum2 += logLRarray[mid];
	      sum3 += origSDarray[0][mid] + origSDarray[1][mid];
	      sum4 += logSDarray[0][mid] + logSDarray[1][mid];
	      printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
		     i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
		     origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
		     logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
	    }
	    fflush(stdout);
	  }
	      
	  assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	}

	if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	  printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	  assert(lcnt == tcnt);
	}
	logLRsum = LRsum;
	if(SDDEBUG)
	  logLRsdsum = LRsdsum;
	score_free(0,numrefmaps);

	delete [] origLRarray;
	for(int c = 0; c < colors; c++){
	  delete [] origSDarray[c];
	}
      } // if(REFDEBUG ... )

      if(REFDEBUG && (numrefmaps==1 || (BestRef && !BestRefPV)) && maxresbias > mres * 0.5){/* check if logLRsum has changed after calling BiasCompute() */
	BiasCorrect(Gmap, 0, nummaps, 0);

	if(DEBUG) assert(logLRarray);
	score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	lcnt= 0;
	double *origLRarray = new double[nummaps];
	double *origSDarray[MAXCOLOR];
	for(int mid = 0; mid < nummaps; mid++)
	  origLRarray[mid] = logLRarray[mid];
	for(int c = 0; c < colors; c++){
	  origSDarray[c] = new double[nummaps];
	  for(int mid = 0; mid < nummaps; mid++){
	    origSDarray[c][mid] = logSDarray[c][mid];
	  }
	}

	double origLRsum = 0.0;
	double origSDsum[MAXCOLOR] = {0.0,0.0};

	for(size_t i = 0; i < numaligns; i++){
	  double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	  Calign *align = alignment[i];
	  if(!align || align->numpairs <= 1)
	    continue;
	  int mid = align->mapid2;
	  if(DEBUG) assert(!Gmap[mid]->origmap);
	  if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
	    if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
	    continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	  }
	  if(align->score <= origScoreThreshold)
	    continue;
	  int rid = align->mapid1;
	  if(DEBUG>=2){
	    int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
	    assert(origgoodalign);
	    int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
	    assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	  }
	  FLOAT **XX = Gmap[mid]->site;
	  int *MM = Gmap[mid]->numsite;
	  FLOAT **YY = refmap[rid]->site;
	  //	int *NN = Gmap[rid]->numsite;
	  align++;
	  logLRarray[mid] = 0.0;
	  for(int c = 0; c < colors; c++){
	    FLOAT *X = XX[c];
	    double M = MM[c];
	    FLOAT *Y = YY[c];
	    int N = refmap[rid]->numsite[c];
	    int U = align[c].numpairs;

	    LR[c] = logLR(mid,rid,c,X,Y,M,N,align,align[c].sites1[0], align[c].sites2[0], align[c].sitesK1[0], align[c].sites2[U-1], align[c].Lij1,align[c].Rij1,align[c].Lij2,align[c].Rij2,1,0,SDDEBUG ? &LRsd[c] : 0,NULL);
	    logLRarrayC[c][mid] = LR[c];
	    logSDarray[c][mid] = LRsd[c];
	    if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
	      continue;
	    LRsum += LR[c];
	    if(SDDEBUG){
	      LRsdsum += LRsd[c];
	      LRsdsumC[c] += LRsd[c];
	    }
	  }

	  int U = align[0].numpairs;
	  int I = align[0].sites1[U-1];
	  int K = align[0].sitesK1[U-1];
	  int J = align[0].sites2[U-1];
	  U = align[1].numpairs;
	  int P = align[1].sites1[U-1];
	  int T = align[1].sitesK1[U-1];
	  int Q = align[1].sites2[U-1];
	  double y0 = Yc(YY[0],I,K);
	  double y1 = Yc(YY[1],P,T);
	  double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	  double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	  double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);
	  
	  LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
	  if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	    LRsum += LR[2];
	    for(int lc = 0; lc <= colors; lc++)
	      LRsumC[lc] += LR[lc];
	  }

	  logLRarray[mid] = LR[0] + LR[1] + LR[2];
	  lcnt++;
	  if(VERB>=2 && giter==0){
	    if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	      origLRsum += origLRarray[mid];
	      for(int lc = 0; lc < colors; lc++)
		origSDsum[lc] += origSDarray[lc][mid];
	    }
	    if(mid==265){
	      printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)(delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f)\n",
		     i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
		     logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
	      printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
		     logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
		     LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
		     LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
	      printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);
	    }
	  }
	}
	if(SDDEBUG)
	  printf("After extra BiasCorrect:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		 (unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt,LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	else
	  printf("After extra BiasCorrect:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
	fflush(stdout);

	if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	  if(VERB){
	    double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
	    for(size_t i = 0; i < numaligns; i++){
	      Calign *align = alignment[i];
	      if(!align || align->numpairs <= 1)
		continue;
	      int mid = align->mapid2;
	      if(DEBUG>=2) assert(!Gmap[mid]->origmap);
	      if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		continue;
	      if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		continue;
	      sum1 += origLRarray[mid];
	      sum2 += logLRarray[mid];
	      sum3 += origSDarray[0][mid] + origSDarray[1][mid];
	      sum4 += logSDarray[0][mid] + logSDarray[1][mid];
	      printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
		     i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
		     origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
		     logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
	    }
	    fflush(stdout);
	  }
	      
	  assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	}

	if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	  printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	  assert(lcnt == tcnt);
	}
	logLRsum = LRsum;
	if(SDDEBUG)
	  logLRsdsum = LRsdsum;
	score_free(0,numrefmaps);

	delete [] origLRarray;
	for(int c = 0; c < colors; c++){
	  delete [] origSDarray[c];
	}
      } // if(REFDEBUG ... )

      if(ScanCorrection && UniqueScans > 1){/* per scan scaling of molecules */
	/* NOTE : per scan scaling of molecules is same for all colors */
	/* HERE HERE : update ScanScale code for 2 colors from refalign.cpp */
	printf("-ScanScaling not yet implemented for 2 colors\n");
	exit(1);
      }/* if(ScanCorrection && UniqueScans > 1) */
    } // if(numerrors[0] + numerrors[1] + numerrors[2] > 0)

    if(RESDATA && numresdata[0]+numresdata[1] > 0){/* optimize res and resSD */
      int numthreads = 1;
      #ifdef _OPENMP
      numthreads = omp_get_max_threads();
      numthreads = min(numthreads,MaxThreads);
      if(numthreads > max(numresdata[0],numresdata[1])/256)
	numthreads = max(1,max(numresdata[0],numresdata[1])/256);
      #endif

      for(int c = 0; c < colors; c++){
     	if(RESDATA_FIX){  /* resdata[] should be updated to reflect -resbias and X scaling C : simply update values from Gmap[]->site[][] */
          #pragma omp parallel for schedule(static,256) num_threads(numthreads) if(numthreads > 1)
	  for(int i = 0; i < numresdata[c]; i++){
	    Cresdata *pres = &resdata[c][i];
	    Calign *align = pres->align;
	    int mapid = align->mapid2;
	    Cmap *Xmap = Gmap[mapid];
	    FLOAT *X = Xmap->site[c];	  
	    int M = Xmap->numsite[c];
	    double escale = ScaleDeltaBPP ? 1.0 : (align->numpairs > 0 && align->scaleID) ? ScaleFactor[align->scaleID] : 1.0;
	    if(pres->n)
	      pres->x = escale * (align->orientation ? X[M+1-(pres->J-pres->m)] - X[M+1-pres->J] : X[pres->J] - X[pres->J - pres->m]);
	    if(pres->Ip && pres->n1)
	      pres->x1 = escale * (align->orientation ? X[M+1-(pres->J1-pres->m1)] - X[M+1-pres->J1] : X[pres->J1] - X[pres->J1 - pres->m1]);
	  }
	}

	double bestRes = res[c];
	double bestResSD = resSD[c];
	double bestLL = resLL(resdata[c], numresdata[c], numthreads, Tarray1, c, 0);
	double origLL = bestLL;

	if(VERB>=2){
	  printf("c=%d:res[c]=%0.4f,resSD[c]=%0.4f:LL=%0.6f: numresdata = %d\n",c,res[c],resSD[c], bestLL, numresdata[c]);
	  fflush(stdout);
	}

	/* perform goldenmean search on res[c] in range 0.50 .. 10.00 */
	double phi = (1.0 + sqrt(5.0))*0.5;
	double resphi = 2.0 - phi;

	double Low = max(0.001,bestRes - 0.200);
	double High = bestRes + 0.200;
	double Mid = bestRes;
	while(max(High - Mid, Mid - Low) > 0.0001){
	  double Lsum;
	  if(High - Mid > Mid - Low){
	    res[c] = Mid + resphi * (High - Low);
	    Lsum = resLL(resdata[c],numresdata[c],numthreads,Tarray1, c, 0);

	    if(Lsum > bestLL){
	      bestLL = Lsum;
	      Low = Mid;
	      Mid = res[c];
	    } else
	      High = res[c];
	  } else {
	    res[c] = Mid - resphi * (High - Low);
	    Lsum = resLL(resdata[c],numresdata[c],numthreads,Tarray1, c, 0);

	    if(Lsum > bestLL){
	      bestLL = Lsum;
	      High = Mid;
	      Mid = res[c];
	    } else 
	      Low = res[c];
	  }
	  if(VERB>=2){
	    printf("c=%d:res[c]=%0.4f,resSD[c]=%0.4f:LL=%0.6f (best:res[c]=%0.4f,resSD[c]=%0.4f,LL=%0.6f), res:Low=%0.4f,Mid=%0.4f,High=%0.4f\n",c,res[c],resSD[c],Lsum, Mid, bestResSD, bestLL,Low,Mid,High);
	    fflush(stdout);
	  }
	}
	if(VERB){
	  printf("  res[%d]= %0.6f -> %0.6f (score= %0.6f -> %0.6f,lcnt=%lu)\n",c,bestRes,Mid,origLL/max((size_t)1,lcnt), bestLL/max((size_t)1,lcnt),lcnt);
	  fflush(stdout);
	}
	bestRes = res[c] = Mid;
      
	resKB[c] = res[c] * PixelLen;/* in case res changed */

	if(REFDEBUG && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	  if(DEBUG) assert(logLRarray);
	  score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	  double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	  double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	  lcnt= 0;
	  double *origLRarray = new double[nummaps];
	  double *origSDarray[MAXCOLOR];
	  for(int mid = 0; mid < nummaps; mid++)
	    origLRarray[mid] = logLRarray[mid];
	  for(int lc = 0; lc < colors; lc++){
	    origSDarray[lc] = new double[nummaps];
	    for(int mid = 0; mid < nummaps; mid++){
	      origSDarray[lc][mid] = logSDarray[lc][mid];
	    }
	  }

	  double origLRsum = 0.0;
	  double origSDsum[MAXCOLOR] = {0.0,0.0};

	  for(size_t i = 0; i < numaligns; i++){
	    double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	    Calign *align = alignment[i];
	    if(!align || align->numpairs <= 1)
	      continue;
	    int mid = align->mapid2;
	    if(DEBUG) assert(!Gmap[mid]->origmap);
	    if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
	      if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
	      continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	    }
	    if(align->score <= origScoreThreshold)
	      continue;
	    int rid = align->mapid1;
	    if(DEBUG>=2){
	      int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
	      assert(origgoodalign);
	      int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
	      assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	    }
	    FLOAT **XX = Gmap[mid]->site;
	    int *MM = Gmap[mid]->numsite;
	    FLOAT **YY = refmap[rid]->site;
	    //	int *NN = Gmap[rid]->numsite;
	    align++;
	    logLRarray[mid] = 0.0;
	    for(int lc = 0; lc < colors; lc++){
	      FLOAT *X = XX[lc];
	      double M = MM[lc];
	      FLOAT *Y = YY[lc];
	      int N = refmap[rid]->numsite[lc];
	      int U = align[lc].numpairs;

	      LR[lc] = logLR(mid,rid,lc,X,Y,M,N,align,align[lc].sites1[0], align[lc].sites2[0], align[lc].sitesK1[0], align[lc].sites2[U-1], align[lc].Lij1,align[lc].Rij1,align[lc].Lij2,align[lc].Rij2,1,
			     0 /* (giter==0 && c==1 && lc==0 && mid <= 2) ? 1 : 0*/, SDDEBUG ? &LRsd[lc] : 0,NULL);
	      logLRarrayC[lc][mid] = LR[lc];
	      logSDarray[lc][mid] = LRsd[lc];
	      if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
		continue;
	      LRsum += LR[lc];
	      if(SDDEBUG){
		LRsdsum += LRsd[lc];
		LRsdsumC[lc] += LRsd[lc];
	      }
	    }

	    int U = align[0].numpairs;
	    int I = align[0].sites1[U-1];
	    int K = align[0].sitesK1[U-1];
	    int J = align[0].sites2[U-1];
	    U = align[1].numpairs;
	    int P = align[1].sites1[U-1];
	    int T = align[1].sitesK1[U-1];
	    int Q = align[1].sites2[U-1];
	    double y0 = Yc(YY[0],I,K);
	    double y1 = Yc(YY[1],P,T);
	    double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	    double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	    double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);
	    
	    LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
	    if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	      LRsum += LR[2];
	      for(int lc = 0; lc <= colors; lc++)
		LRsumC[lc] += LR[lc];
	    }

	    logLRarray[mid] = LR[0] + LR[1] + LR[2];
	    lcnt++;
	    if(VERB>=2 && giter==0 && c==1){
	      if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		origLRsum += origLRarray[mid];
		for(int lc = 0; lc < colors; lc++)
		  origSDsum[lc] += origSDarray[lc][mid];
	      }
	      if(1){
		printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)\n\t delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f\n",
		       i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
		       logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
		/*		printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
		       logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
		       LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
		       LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
		       printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);*/
	      }
	    }
	  }
	  if(SDDEBUG)
	    printf("After updating res[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		   c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt,LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	  else
	    printf("After updating res[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
	  fflush(stdout);

	  if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	    if(VERB){
	      double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
	      for(size_t i = 0; i < numaligns; i++){
		Calign *align = alignment[i];
		if(!align || align->numpairs <= 1)
		  continue;
		int mid = align->mapid2;
		if(DEBUG>=2) assert(!Gmap[mid]->origmap);
		if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		  continue;
		if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		  continue;
		sum1 += origLRarray[mid];
		sum2 += logLRarray[mid];
		sum3 += origSDarray[0][mid] + origSDarray[1][mid];
		sum4 += logSDarray[0][mid] + logSDarray[1][mid];
		printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
		       i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
		       origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
		       logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
	      }
	      fflush(stdout);
	    }
	      
	    assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	  }

	  if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	    printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	    assert(lcnt == tcnt);
	  }
	  logLRsum = LRsum;
	  if(SDDEBUG)
	    logLRsdsum = LRsdsum;
	  score_free(0,numrefmaps);

	  delete [] origLRarray;
	  for(int c = 0; c < colors; c++){
	    delete [] origSDarray[c];
	  }
	} // if(REFDEBUG ... )

	/* perform goldenmean search on resSD */
	if(VERB>=2 && giter==0 && c==1){
	  bestLL = resLL(resdata[c], numresdata[c], numthreads, Tarray1, c, 1);

	  printf("c=%d:res[c]=%0.4f,resSD[c]=%0.4f:LL=%0.6f: numresdata = %d\n",c,res[c],resSD[c], bestLL, numresdata[c]);
	  fflush(stdout);
	}

	origLL = bestLL;
	Low = max(0.001,bestResSD - 0.100);
	High = bestResSD + 0.100;
	Mid = bestResSD;
	while(max(High - Mid, Mid - Low) > 0.0001){
	  double Lsum;
	  if(High - Mid > Mid - Low){
	    resSD[c] = Mid + resphi * (High - Low);
	    Lsum = resLL(resdata[c],numresdata[c],numthreads,Tarray1,c,0);

	    if(Lsum > bestLL){
	      bestLL = Lsum;
	      Low = Mid;
	      Mid = resSD[c];
	    } else
	      High = resSD[c];
	  } else {
	    resSD[c] = Mid - resphi * (High - Low);
	    Lsum = resLL(resdata[c],numresdata[c],numthreads,Tarray1,c,0);

	    if(Lsum > bestLL){
	      bestLL = Lsum;
	      High = Mid;
	      Mid = resSD[c];
	    } else 
	      Low = resSD[c];
	  }
	  if(VERB>=2){
	    printf("c=%d:res[c]=%0.4f,resSD[c]=%0.4f:LL=%0.6f (best:res[c]=%0.4f,resSD[c]=%0.4f,LL=%0.6f), resSD:Low=%0.4f,Mid=%0.4f,High=%0.4f\n",c,res[c],resSD[c],Lsum, bestRes, Mid,bestLL,Low,Mid,High);
	    fflush(stdout);
	  }
	}
	if(VERB){
	  printf("  resSD[%d]= %0.8f -> %0.8f (score= %0.8f -> %0.8f, lcnt=%lu)\n",c,bestResSD,Mid,origLL/max((size_t)1,lcnt), bestLL/max((size_t)1,lcnt),lcnt);
	  fflush(stdout);
	}
	bestResSD = resSD[c] = Mid;

	IresSD[c] = 1.0/(resSD[c] * PixelLen * sqrt(2.0));/* in case resSD changed */

	if(REFDEBUG && (numrefmaps==1 || (BestRef && !BestRefPV))){/* check if logLRsum has improved */
	  if(DEBUG) assert(logLRarray);
	  score_init(refmap,0,numrefmaps,nummaps,RKKmax);      
	  double LRsum = 0.0, LRsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	  double LRsdsum = 0.0, LRsdsumC[MAXCOLOR+1] = {0.0,0.0,0.0};
	  lcnt= 0;
	  double *origLRarray = new double[nummaps];
	  double *origSDarray[MAXCOLOR];
	  for(int mid = 0; mid < nummaps; mid++)
	    origLRarray[mid] = logLRarray[mid];
	  for(int lc = 0; lc < colors; lc++){
	    origSDarray[lc] = new double[nummaps];
	    for(int mid = 0; mid < nummaps; mid++){
	      origSDarray[lc][mid] = logSDarray[lc][mid];
	    }
	  }

	  double origLRsum = 0.0;
	  double origSDsum[MAXCOLOR] = {0.0,0.0};

	  for(size_t i = 0; i < numaligns; i++){
	    double LR[MAXCOLOR+1], LRsd[MAXCOLOR];
	    Calign *align = alignment[i];
	    if(!align || align->numpairs <= 1)
	      continue;
	    int mid = align->mapid2;
	    if(DEBUG) assert(!Gmap[mid]->origmap);
	    if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1){
	      if(DEBUG) assert(numrefmaps > 1 && !BestRefPV);
	      continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
	    }
	    if(align->score <= origScoreThreshold)
	      continue;
	    int rid = align->mapid1;
	    if(DEBUG>=2){
	      int origgoodalign = AlignedThreshold(align, refmap[rid]->site, origScoreThreshold, origLogPvThreshold);
	      assert(origgoodalign);
	      int goodalign = (origgoodalign && align->score > ScoreThreshold && align->logPV > LogPvThreshold) ? 1 : 0;
	      assert((MapRate > 0.0 && !(align->score > ScoreThreshold)) == (!goodalign));
	    }
	    FLOAT **XX = Gmap[mid]->site;
	    int *MM = Gmap[mid]->numsite;
	    FLOAT **YY = refmap[rid]->site;
	    //	int *NN = Gmap[rid]->numsite;
	    align++;
	    logLRarray[mid] = 0.0;
	    for(int lc = 0; lc < colors; lc++){
	      FLOAT *X = XX[lc];
	      double M = MM[lc];
	      FLOAT *Y = YY[lc];
	      int N = refmap[rid]->numsite[lc];
	      int U = align[lc].numpairs;

	      LR[lc] = logLR(mid,rid,lc,X,Y,M,N,align,align[lc].sites1[0], align[lc].sites2[0], align[lc].sitesK1[0], align[lc].sites2[U-1], align[lc].Lij1,align[lc].Rij1,align[lc].Lij2,align[lc].Rij2,1,
			     /* (mid==16) ? 1 : 0*/ 0, SDDEBUG ? &LRsd[lc] : 0,NULL);
	      logLRarrayC[lc][mid] = LR[lc];
	      logSDarray[lc][mid] = LRsd[lc];
	      if(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))
		continue;
	      LRsum += LR[lc];
	      if(SDDEBUG){
		LRsdsum += LRsd[lc];
		LRsdsumC[lc] += LRsd[lc];
	      }
	    }

	    int U = align[0].numpairs;
	    int I = align[0].sites1[U-1];
	    int K = align[0].sitesK1[U-1];
	    int J = align[0].sites2[U-1];
	    U = align[1].numpairs;
	    int P = align[1].sites1[U-1];
	    int T = align[1].sitesK1[U-1];
	    int Q = align[1].sites2[U-1];
	    double y0 = Yc(YY[0],I,K);
	    double y1 = Yc(YY[1],P,T);
	    double x0 = align[-1].orientation ? XX[0][MM[0]+1] - XX[0][MM[0]+1-J] : XX[0][J];
	    double x1 = align[-1].orientation ? XX[1][MM[1]+1] - XX[1][MM[1]+1-Q] : XX[1][Q];
	    double len = (COLOR_SHIFT==2) ? max(y1,y0) : max(x1,x0);
	    
	    LR[2] = SMA(x1-x0, y1-y0, len, align[-1].orientation);
	    if(VERB>=2 && giter==0 && c==1 && mid==16){
	      printf("rid=%d,mid=%d,or=%d,c=2:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,y=%0.4f,len=%0.4f,SMA=%0.6f:MA_mean=%0.4f,SF[2]=%0.8e,FM2var=%0.8e,SD[2]=%0.8e,SM2var=%0.8e,RefPen12=%0.8e,Xtheta12=%0.4f\n",
		     rid,mid,align[-1].orientation,I,K,J,P,T,Q,x0,x1,y0,y1,y1-y0,len,LR[2],MA_mean,SF[2],FM2var,SD[2],SM2var,RefPen12,Xtheta12);
	      fflush(stdout);
	    }
	    if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
	      LRsum += LR[2];
	      for(int lc = 0; lc <= colors; lc++)
		LRsumC[lc] += LR[lc];
	    }

	    logLRarray[mid] = LR[0] + LR[1] + LR[2];
	    lcnt++;

	    if(VERB>=2 && giter==0 && c==1 && mid==16){
	      if(!(MapRate > 0.0 && !(align[-1].score > ScoreThreshold))){
		origLRsum += origLRarray[mid];
		for(int lc = 0; lc < colors; lc++)
		  origSDsum[lc] += origSDarray[lc][mid];
	      }
	      if(1){
		printf("  i=%lu/%lu:refid=%d,mapid=%d,or=%d,pairs=%d,logLR=%0.6f(%0.6f,%0.6f,%0.6f),LRsd=%0.6f(%0.6f,%0.6f)\n\t delta=%0.6f= %0.6f,%0.6f,%0.6f,LRsd=%0.6f= %0.6f,%0.6f\n",
		       i,numaligns,rid,mid,align[-1].orientation,align[-1].numpairs,LRsum,LRsumC[0],LRsumC[1],LRsumC[2],LRsdsum,LRsdsumC[0],LRsdsumC[1],
		       logLRarray[mid],LR[0],LR[1],LR[2],LRsd[0]+LRsd[1],LRsd[0],LRsd[1]);
		/*		printf("\t changes:logLR=%0.6f(prev=%0.6f), LRsd=%0.6f(prev=%0.6f)(%0.6f,%0.6f),logLR-LRsd=%0.6f(cum LR=%0.6f,LRsd=%0.6f,LR-LRsd=%0.6f)\n",
		       logLRarray[mid] - origLRarray[mid], origLRarray[mid], LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid], origSDarray[0][mid] + origSDarray[1][mid],
		       LRsd[0] - origSDarray[0][mid], LRsd[1] - origSDarray[1][mid], logLRarray[mid] - origLRarray[mid] - (LRsd[0]+LRsd[1] - origSDarray[0][mid] - origSDarray[1][mid]),
		       LRsum - origLRsum, LRsdsum - origSDsum[0] - origSDsum[1], LRsum - origLRsum - (LRsdsum - origSDsum[0] - origSDsum[1]));
		       printf("\t x0=%0.4f,x1=%0.4f,y0=%0.4f,y1=%0.4f,x=%0.4f,y=%0.4f,len=%0.4f,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d,lcnt=%llu\n",x0,x1,y0,y1,x1-x0,y1-y0,len,I,K,J,P,T,Q,(unsigned long long)lcnt);*/
	      }
	    }
	  }
	  if(SDDEBUG)
	    printf("After updating resSD[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)(logLRsd=%0.6f -> %0.6f = %0.6f,%0.6f)\n",
		   c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt,logLRsdsum/lcnt,LRsdsum/lcnt,LRsdsumC[0]/lcnt,LRsdsumC[1]/lcnt);
	  else
	    printf("After updating resSD[%d]:lcnt=%llu:logLR=%0.6f -> %0.6f(%0.6f,%0.6f,%0.6f)\n",
		   c,(unsigned long long)lcnt,logLRsum/lcnt,LRsum/lcnt,LRsumC[0]/lcnt,LRsumC[1]/lcnt,LRsumC[2]/lcnt);
	  fflush(stdout);

	  if(DEBUG && !(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt)){
	    if(VERB){
	      double sum1 = 0.0,sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
	      for(size_t i = 0; i < numaligns; i++){
		Calign *align = alignment[i];
		if(!align || align->numpairs <= 1)
		  continue;
		int mid = align->mapid2;
		if(DEBUG>=2) assert(!Gmap[mid]->origmap);
		if(BestRef && Gmap[mid]->align->mapid1 != align->mapid1)
		  continue;
		if(MapRate > 0.0 && !(align->score > ScoreThreshold))
		  continue;
		sum1 += origLRarray[mid];
		sum2 += logLRarray[mid];
		sum3 += origSDarray[0][mid] + origSDarray[1][mid];
		sum4 += logSDarray[0][mid] + logSDarray[1][mid];
		printf("i=%lu/%lu:mid=%d,id=%lld:logLRarray[mid]=%0.8f -> %0.8f (delta=%0.8f,cum=%0.8f -> %0.8f), LLsd=%0.8f -> %0.8f(delta=%0.8f,cum=%0.8f -> %0.8f)\n",
		       i,numaligns,mid,Gmap[mid]->id,origLRarray[mid],logLRarray[mid], logLRarray[mid]-origLRarray[mid],sum1, sum2,
		       origSDarray[0][mid]+origSDarray[1][mid],logSDarray[0][mid]+logSDarray[1][mid],
		       logSDarray[0][mid]+logSDarray[1][mid]-origSDarray[0][mid]-origSDarray[1][mid],sum3,sum4);
	      }
	      fflush(stdout);
	    }
	      
	    assert(LRsum > logLRsum - (USE_RFLOAT ? 1e-4 : 1e-6)*lcnt);
	  }

	  if(DEBUG && MapRate <= 0.0 && !(lcnt==tcnt)){
	    printf("lcnt=%lu,tcnt=%lu\n",lcnt,tcnt);
	    assert(lcnt == tcnt);
	  }
	  logLRsum = LRsum;
	  if(SDDEBUG)
	    logLRsdsum = LRsdsum;
	  score_free(0,numrefmaps);

	  delete [] origLRarray;
	  for(int c = 0; c < colors; c++){
	    delete [] origSDarray[c];
	  }
	} // if(REFDEBUG ... )
      }// c = 0 .. colors-1
    } // if(RESDATA ...)

    /* save updated parameter values */
#pragma novector
    for(int c = 0; c < colors;c++){
      parameters[giter+1].res[c] = res[c];
      parameters[giter+1].resSD[c] = resSD[c];

      parameters[giter+1].FP[c] = FP[c];
      parameters[giter+1].FN[c] = FN[c];
      parameters[giter+1].SF[c] = SF[c];
      parameters[giter+1].SD[c] = SD[c];
      if(QUADRATIC_VARIANCE)
	parameters[giter+1].SR[c] = SR[c];
      if(RES_VARIANCE)      
	parameters[giter+1].SE[c] = SE[c];
      parameters[giter+1].minSNR[c] = minSNR[c];

      parameters[giter+1].ResBins[c] = ResBins[c];
      for(int Bin = 0; Bin <= ResBins[c]; Bin++){
	parameters[giter+1].resbias[c][Bin] = resbias[c][Bin];
	parameters[giter+1].resbiasX[c][Bin] = resbiasX[c][Bin];
	if(DEBUG/* HERE >=2 */) assert(maxresbias >= resbiasX[c][Bin] - 1e-6);	
      }
      if(DEBUG/* HERE >=2 */ && ResBins[c] > 0) assert(maxresbias >= resbiasX[c][ResBins[c]] - 1e-6);
      for(int Bin = ResBins[c] + 1; Bin <= origResBins[c]; Bin++){
	parameters[giter+1].resbias[c][Bin] = resbias[c][ResBins[c]];
	parameters[giter+1].resbiasX[c][Bin] = resbiasX[c][ResBins[c]];
      }
    }
    parameters[giter+1].SF[colors] = SF[colors];
    parameters[giter+1].SD[colors] = SD[colors];
    parameters[giter+1].MA_mean = MA_mean;

    parameters[giter+1].PixelLen = PixelLen;
    parameters[giter+1].PixelLenSD = PixelLenSD;
    parameters[giter+1].mres = mres; // NOTE : mres does not change during autonoize
    parameters[giter+1].mresSD = mresSD; // NOTE : mresSD does not change during autonoise

    parameters[giter+1].mapcnt = 0;
    parameters[giter+1].logLR = 0.0;
    parameters[giter+1].ATmapcnt = 0;
    parameters[giter+1].ATlogLR = 0.0;
    parameters[giter+1].sitecnt = 0;
    parameters[giter+1].outlierRate = 0.0;
    parameters[giter+1].EndoutlierRate = 0.0;
    for(int c = 0; c < colors; c++)
      parameters[giter+1].LabelDensity[c] = 0.0;

    parameters[giter+1].maxresbias = maxresbias;

    if(giter < RefRepeats-1){ /* check if parameters changed significantly */
      double maxerr = 0.0;
      for(int c = 0; c < colors; c++){
	if(fabs(FP[c] - FPorig[c])/FPorig[c] > maxerr)
	  maxerr = fabs(FP[c]-FPorig[c])/FPorig[c];
	if(fabs(FN[c] - FNorig[c])/FNorig[c] > maxerr)
	  maxerr = fabs(FN[c]-FNorig[c])/FNorig[c];
	if(fabs(SF[c] - SForig[c])/SForig[c] > maxerr)
	  maxerr = fabs(SF[c] - SForig[c])/SForig[c];
	if(fabs(SD[c] - SDorig[c])/fabs(SDorig[c]) > maxerr)
	  maxerr = fabs(SD[c] - SDorig[c])/fabs(SDorig[c]);
	if(SRorig[c] > 0.0 && fabs(SR[c] - SRorig[c])/SRorig[c] > maxerr)
	  maxerr = fabs(SR[c] - SRorig[c])/SRorig[c];
	if(SEorig[c] > 0.0 && fabs(SE[c] - SEorig[c])/SEorig[c] > maxerr)
	  maxerr = fabs(SE[c] - SEorig[c])/SEorig[c];
      }
      if(fabs(PixelLen - startPixelLen)/startPixelLen > maxerr)
	maxerr = fabs(PixelLen - startPixelLen)/startPixelLen;
      
      if(maxerr < 1.0e-5){
	if(VERB){
	  printf("max relative parameter change = %0.7f : skipping remaining iterations\n", maxerr);
	  fflush(stdout);
	}
	RefRepeats = giter+2;
      }
    }// if(giter < RefRepeats-1)
    score_free(0,numrefmaps);
  }// for giter = 0 .. RefRepeats -1

  delete [] scaleIDcnt;    scaleIDcnt = 0;

  delete [] RKKmax;
  delete [] RKmem;
  delete [] RKKmem;

  delete [] Tarray1;
  delete [] Tarray2;
  delete [] Tarray3;
  delete [] Tarray4;
  delete [] Tarray5;
  delete [] Tarray6;

  if(HashFP){
    hash_close(HashFP,1);

    numhashpairs1 = numhashpairs2 = 0;
    NumPair1 = NumPair2 = 0;
    nexthash1 = hashpairs1;
    refid1 = refidend1 = -1;
    nexthash2 = hashpairs2;
    refid2 = refidend2 = -1;
    HashFP = NULL;
  }

  if(!Refine && !NoStat){
    /* output parameter values after each iteration in .err file and .errbin files */
    output_err(output_prefix,parameters,RefRepeats,MappingRatePV,CoverageMult,(int)floor(max(0.0,origLogPvThreshold) + 0.5), MAX_PV);

    if(ScanCorrection && UniqueScans > 1 && giter2 == RefRepeats2 - 1){/* save rawsite[] as <prefix>_rescaled.bnx (without -bpp or -resbias correction) */
      char basename[PATH_MAX];
      sprintf(basename,"%s_rescaled",output_prefix);

      output_bnx(basename,Gmap,0,startmaps-1,1/*raw*/, NULL, NULL, -1);
    }
  }

  delete rs_heap; rs_heap = NULL;

  delete [] hashmatch1; hashmatch1 = 0;
  refid1 = refid2 = -1;

  delete [] hashmatch2; hashmatch2 = 0;
  delete [] hashpairs1; hashpairs1 = 0;
  delete [] hashpairs2; hashpairs2 = 0;
  numhashpairs1 = numhashpairs2 = maxhashpairs1 = maxhashpairs2 = 0;

  delete [] YidListMem; YidListMem = NULL;
  delete [] phashListMem; phashListMem = NULL;
  YidList1 = YidList2 = XidList1 = XidList2 = 0;
  phashList1 = phashList2 = 0;
  NumPair1 = NumPair2 = MaxPair = 0;
  refidend1 = refidend2 = -1;

  for(int c = 0; c < colors; c++){
    delete [] SizeToBin[c]; SizeToBin[c] = NULL;
  }

  //  delete [] origmapid; origmapid = NULL;
  delete [] numalign_start;
  delete [] numalign_end;
  delete [] orignummaps;
  numalign_start = numalign_end = 0;
  orignummaps = 0;

  delete [] parameters;

  delete [] logLRarray; logLRarray = NULL;
  for(int c = 0; c < colors; c++){
    delete [] logLRarrayC[c]; logLRarrayC[c] = NULL;
    delete [] logSDarray[c]; logSDarray[c] = NULL;
    delete [] Xlen[c]; Xlen[c] = NULL;
    delete [] resdata[c]; maxresdata[c] = 0; resdata[c] = NULL;
  }
  for(int c = 0; c <= colors; c++){
    delete [] errors[c]; maxerrors[c] = 0; errors[c] = NULL;
  }

  if(ScanCorrection){
    for(int c = 0; c < colors; c++){
      delete [] ThetaDelta[c];
      for(int t = 0; t < numthreads; t++)
	delete [] Txy[c][t];
      delete [] Txy[c];
    }
  }

  /* restore original thresholds in case they were changed by -MapRate and RefRepeats2 > 1 */
  ScoreThreshold = origScoreThreshold;
  LogPvThreshold = origLogPvThreshold;
}
	      
