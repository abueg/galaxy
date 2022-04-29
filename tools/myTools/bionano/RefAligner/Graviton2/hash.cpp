#include <stdlib.h>
#include <stdio.h>
#ifndef WIN32
#include <unistd.h>
#include <sched.h>
#endif
#include <ctype.h>
#include <omp.h>
#include <float.h>

#ifndef __ARM_ARCH
#include <immintrin.h>
#include <emmintrin.h>
#ifdef __SSE3__
#include <pmmintrin.h>
#ifdef __SSSE3__
#include <tmmintrin.h>
#endif
#endif
#endif

//#include <linux/getcpu.h>

#include <sys/mman.h>

#include "parameters.h"
#include "hash.h"
#include "Ccontig.h"

#ifdef  __INTEL_COMPILER
#define CAST (const char *)
#else
#ifdef WIN32
#define CAST (const char *)
#else
#define CAST
#endif
#endif

#ifndef USE_HSSE
#define USE_HSSE 1 // Use SSE2 & SSE3 instruction intrinsics in hashtable
#endif

static RFLOAT aNaN = nan("NaN");

//#undef DEBUG
//#define DEBUG 2

#if 0 // The following constants have been converted to template parameters
#define OFFSET_T int

#define OFFSET1_T unsigned short /* good for query maps up to 320Mb (if OffsetKB = 5kb) */

#define HASHWIN 5 /* compile time value of HashWin (must match commandline value) hardwired for speed */

#define RANGE 1 /* keep track of approximate offset1 range for hashtable matches (if -hashRange is specified) */

#define HASH_SCALE 1 /* >= 1 : If -hashScaleDelta AND -ScaleDelta are specified, try different scalings of probe map returning the best scaling for each offset */
#endif

#define SCALE_SMALL 1 /* NEW173 : prefer smaller scaleID values over larger ones */

#define QUICKHASH 1 /* 1 : use two or three 1-bit tables indexed by key to filter out keys with no entries in hashtable
		       2 : use a single large 1-bit table, indexed by shuffled key (so 243 keys differ mostly in lower 10 bits), to filter out all keys with no entries in hashtable */

#define QUICKSTAT 0 /* 1 : compute average number of cache lines of indexR[] required per call to hashfindV() */
static size_t hashfindV_calls = 0, hashfindV_lines = 0;
static int hashfindV_max = 0;

extern int HashT2;/* see parameters.h */
#define MULTIMATCH_TWOTHREADS (!MSAN && HashT2) /* WAS 0 */ /* Use 2 threads for each hash_query */

extern int HashGC;/* see parameters.h */

#define STREAM_DEBUG 0 /* More verbose output of hash matches being streamed in (2 = Very verbose, each match is displayed) */

#define DETERMINISTIC 0  /* for consistency of results with different mapsets for the same rmap[mapid1]->id values, sort MatchBuf[0..MatchCnt-1] by offset, so identical scores are resolved in favor of smallests offset */

#define MAXCOLLIDE 2048 /* largest number of hashkeys per hash index */

#define NEWHASH 2  /* 1 : Use Cacheline friendly compacted hashtable (read-only) 
		      2 : Use Software Pipelining & prefetch for HashTable query */

#define SIZE_BITS 6 /* number of bits to encode quantized size of one site interval */

#define RHASH_BITS 15 /* 15 is sufficient for HashWin=5,SIZE_BITS=6 */ /* number of bits used in tables to randomize the mapping of key to hashtable index (not more than 29) */

/* The following 3 values are set in hash_init() based on numrmaps and based on Commandline overrides : The defaults shown are for numrmaps < 128000 and no commandline overrides */
static int HASH_BITS; /* number of bits to use in hashtable primary index (not more then 31) (Commandline override -Hash_Bits) */
static int MHASH_BITS; /* number of bits to use in Match Table hash index (not more than 31) (Commandline override -MHash_Bits) */
static int HHASH_BITS; /* number of bits to use in Hit Table hash index (not more than 30) (Commandline override -HHash_Bits) */

#define EXACT_SIZING 1 /* 1 : check all hashtable hits using exact sizing errors (requires larger value for -SDerr for similar TP matching rate TPprob : SDerr = sqrt(ChiSq(Hashwin,TPprob)),
			  where ChiSq(n,p) is the value at which a Chi-Square distribution of n degrees of freedom has a cumulative probability p. Typically TPprob should be >= 0.9 */
#define SIZING_BND 1 /* 1 : check all hashtable hits to ensure sizing errors are within bounds given by SDmax (requires EXACT_SIZING=1 : To not use SDerr, just set SDerr= 1000) */

#define SORT_KEYS 0  /* 1 : sort collision list by hash key to speed up searching the collision arrays
			(only useful if HASH_BITS is reduced due to limited memory, resulting in longer collision arrays) */

#define HASHSCORE_SHIFT 2 /* scale weighted hashscore by (1 << HASHSCORE_SHIFT) */

/* NOTE : hash_init() changes the weights of OFFSET_WT[2..OffsetSpread] to min(OFFSET_WT[0],OFFSET_WT[1]) if pairwise==0.
          Total weight is divided by OFFSET_WT[1] == (1<<HASHSCORE_SHIFT) */

#define OFFSET_FILL 1 /* If OFFSET_SPREAD : average around offset values that had no hits (zero scores) in case multiple neighboring offset values are non-zero */
#define OFFSET_ZEROFIX 1 /* avoid duplicate zero scores in MatchBuf, when adding a zero score entry in MatchTable : this is needed to make RANGE work */

#define MAXHIT_COVERAGE 1 /* If hashMaxCov > 0 : A quick way to suppress repeat data which slows down the HashTable (If the number of entries for a key exceeds HashMaxCov it is supressed) */

/* NOTE : redefining MAX_WIN(OFFSET1_T) to larger values then below will disabled optimized code */
#define MAX_WIN(OFFSET1_T) ((sizeof(OFFSET1_T)==2) ? 6 : 5) /* largest size of window supported (cannot be larger than 14, smaller values use less memory) */
#define MAXR_WIN 8 /* MAX_WIN rounded up to the nearest power of 2 : used for aligned arrays on the stack */

#define MAP_INDEX_TXT 0 /* MAP_INDEX_TXT==0 : -hashtxt output will contain map id instead of map index (useful for debugging) */

#define BNX_OUTPUT SubsetMaps /* If !pairwise, also output subset of input rmap[0..numrmaps-1] with hashtable matches to <prefix>_hash.bnx */

#define STAT 0 // (RELEASE ? 0 : 1) /* compute hashtable performance statistics (slight slowdown) */

#define QUICK 0 /* For debugging : terminate after querying this many maps */

static const unsigned int MINMEM= 2 /* WAS7 4 */; /* initial size of each CHashCollide allocation block in hashtable (must be power of 2 AND >= 2) */
static const unsigned int MINMEM2= 1 /*  WAS8 (USE_MIC ? 1 : 2) */ /* WAS7 8 */;/* initial size of each CHashEntry allocation block in hashtable (must be power of 2 AND >= 1) */

#define __STDC_LIMIT_MACROS
#include <stdint.h>

#define TRACE 0 /* Turn on tracing of hashscore computation for specific map pair : trace HashTable hits (that match both TraceID1 & TraceID2) */
#define TRACE2 0 /* If TRACE : Also trace HashTable inserts (that match TraceID1) and (if TRACE2 >= 2) queries (that match TraceID2) even if no matches are found */

static const long long TraceID1 = 2060LL, TraceID2 = 12LL /* 23LL */;/* Note: In the output  ID1 ID2 are reversed : Typically TraceID2 < TraceID1 if maps were sorted in ascending order of id */
static const int TraceOR2 = 0;/* Will trace hashscore computation for rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 and flip2 == TraceOR2 */
static const int TraceOffset = INT32_MIN;/* If != INT32_MIN, limit tracing to cases with offset == TraceOFFset * OffsetKB */
static const int TraceQKey[MAXR_WIN] = {63,63,63,63,63};/* trace this query key by listing all sizing variations keys used to probe the hashtable */

const double QLen = 0.01; /* Length Quantum in KB */
const double QLenInv = 1.0/QLen;

//static FLOAT maxdist = 2.0;

extern int DisJoint;/* 0 if rmap[],qmap[] are identical OR have shared maps (only combinations with rmap[]->id > qmap[]->id are aligned), 
                       1 if rmap[],qmap[] are Disjoint sets of maps with rmap[] having larger ids than qmap[] */

extern int hnumthreads;/* number of threads to use in parallel sections of hashtable */

extern double MaxFragLen;/* largest interval spanning max(NumErrI,NumErrQ) intermediate sites */
extern unsigned int *IntLen;/* IntLen[0..MaxIntLen-1] : Table to convert any size X between resKB < X <= MaxFragLen to a small integer (quantized representation of length) */

extern int MaxIntLen;
extern unsigned int SizeLen;
extern unsigned int HitHWM;
extern int MemHWM, HWMid,BestHWM;
extern long long BufHWM;
extern size_t OutHWM;
extern long long TotalMatchCnt;/* progress counter */
extern long long gHTindcnt, gHTqindcnt, gHTfindcnt, gHTcolcnt, gHTqcolcnt, gHTmatchcnt, gHTentrycnt;/* HashTable performance statistics */
extern long long gMTcnt, gMTcnt2, gMTcnt3, gMTcol, gMTqcol,gMTicnt;/* MatchTable performance statistics */
extern long long gHTcnt, gHTcol, gHTdup;/* HitTable performance statistics */
extern FLOAT **YY;/* YY[mapid1] == map[mapid1]->site[hcolor] */
extern FLOAT *Ylen;/* Ylen[mapid1] == YY[mapid1][N1 == map[mapid1]->numsite[hcolor]+1] */

//extern OFFSET_T *maxOffset1; /* Used with flip2 == 0 : maxOffset1[mapid1]  == floor(YY[mapid1][max(0,N1 - HashWin)] / OffsetKB + 0.5) */
//extern OFFSET_T *minOffset1R;/* Used with flip2 == 1 : minOffset1R[mapid1] == floor((Ylen[mapid1] - YY[mapid1][N1]) / OffsetKB + 0.5) */

extern int ShiftOffset1;/* If RANGE : Offset1 values are shifted right by this many bits, to get binned Offset1 values (see -hashGrouped) */

extern int MaxSite;
extern int MinOffset,MaxOffset;
extern long long TotWin;

extern int splitcnt;
extern int roffset;

#include <stdlib.h>
// #include <mcheck.h>

/** Hash of single window of rmap. */
template<class OFFSET1_T>
class CHashEntry {/* 32 bytes == half cacheline (if HASHWIN > 0) */
 public:
  PFLOAT win1[MAX_WIN(OFFSET1_T)];
  int mapid1;/* insert map rmap[mapid1] */
  OFFSET1_T offset1; /**< floor(Y[mapid1][site1]/OffsetKB + 0.5) */
  OFFSET1_T Roffset1; /**< floor((Ylen[mapid1] - Y[mapid1][site1])/OffsetKB + 0.5) */

  inline void init(PFLOAT *w1, int id1, OFFSET1_T off1, OFFSET1_T Roff1, int s, int e) {
    memcpy(win1,w1,HashWin*sizeof(PFLOAT));
    mapid1 = id1;
    offset1 = off1;
    Roffset1 = Roff1;
  };
};

/** Container for coincident hash entries. Used during 1st stage (map insertion), along with CHashCollide *index[] */
class CHashCollide {
 public:
  long long key;/**< 64 bit hashkey */
  unsigned int hash;/**< start index into CHashEntry EntryBlock[]. EntryBlock[hash+(0..n-1)] are in use. The smallest power of 2 >= max(n,MINMEM2) is allocated */
  unsigned int n;/**< number of CHashEntry entries (or number of collision entries for first entry pointed to by index[]) */

  void init(int N, long long k, unsigned int h) {n = N; key = k; hash = h;};
};

#define COLLIDELINE 4 /* number of HashCollide entries per CHashIndexR cache line */

/** Used during 2nd stage as index and first 4 coincident hash entries (condensed), along with CHashCollideR CollideMemR[] (for overflow) */
class CHashIndexR {
 public:
  long long key[COLLIDELINE];
  unsigned int hash[COLLIDELINE];/**< hash[i] is index into CHashEntry EntryMem[] for key[i] */
  unsigned short n[COLLIDELINE];/**< n[i] is number of HashTable entries for key[i] : EntryMem[hash[i] + (0 ..n[i]-1])] */
  unsigned int cnt;/**< number Collision entries : If cnt > COLLIDELINE, the overflow is stored in CollideMemR[next..next+cnt-COLLIDELINE-1]  */
  unsigned int next;
};

#define CHashCollideR CHashCollide

/** Single hash result. A vector of CHashResult is returned holding all hits for a single query (mapid2,site2) */
class CHashResult {
 public:
  unsigned int hash;/* index into EntryMem[] of first hashtable entry found */
  unsigned int n;/* number of hashtable entries found */
};

/* number of MatchTable entries per MatchTable cache line */
#define MATCHLINE(OFFSET_T, RANGE, HASH_SCALE) (HASH_SCALE ? (sizeof(OFFSET_T)==2 ? (RANGE ? 5 : 7) : (RANGE ? 4 : 5)) : (sizeof(OFFSET_T)==2 ? (RANGE ? 6 : 8) : (RANGE ? 5 : 6)))

/* bytes of padding required to make CMatchLine equal to 64 bytes (1 cache line) */
#define MPADDING(OFFSET_T, RANGE, HASH_SCALE) (56 - MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) * (sizeof(int) + sizeof(OFFSET_T) + (RANGE ? sizeof(unsigned short) : 0) + (HASH_SCALE ? 2 : 1) * sizeof(unsigned char)))

/** 64-byte Cache line friendly MatchTable Line with MATCHLINE entries + links to additional Lines */
template<class OFFSET_T, int RANGE, int HASH_SCALE>
class CMatchLine { 
 public:
  int mapid1[MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)];/* mapid1[0..cnt-1] */
  OFFSET_T offset[MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)];/* offset[0..cnt-1] */
  unsigned short Offset1[RANGE ? MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) : 0];/* If RANGE : Offset1[0..cnt-1] : binned Offset1 values 1 + (offset1 >> ShiftOffset1) */
  unsigned char score[MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)];/* score[0..cnt-1] : score counts saturate at 255 */
  unsigned char scaleID[HASH_SCALE ? MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) : 0];/* scaleID[0..cnt-1] : scaleID value ranging from 0 to 255 */
  char padding[MPADDING(OFFSET_T,RANGE,HASH_SCALE)];
  unsigned int cnt;/* Number of entries used in this Cache line (up to MATCHLINE) : For MatchIndex[] if this value is > MATCHLINE, then this value holds the first index in MatchOverflow[] */
  int next;/* next CMatchEntry in MatchOverflow[] ; For MatchIndex[] this value is the next index in MatchIndex[] with cnt > 0 (or -1 if none) */
  void init(int id1, OFFSET_T off, unsigned char s, unsigned short Soffset1, unsigned char ScaleID) { 
    mapid1[0] = id1; offset[0] = off; score[0] = s; cnt = 1; 
    if(RANGE) Offset1[0] = Soffset1; 
    if(HASH_SCALE) scaleID[0] = ScaleID;
  };
};

template <class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
class CMatchTable;/* advance declaration */

#include "resource.h"

/** Primary Data and methods for hash table generation of rmaps.
*
* NOTE : One copy of the CHashTable can be shared (read-only) by multiple threads each with their own copy of MatchTable & HitTable.
          HashTable is populated with maps in the master thread then never modified during the query phase */

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
class CHashTable {
 public:
  Cmap **rmap;/**< map array : rmap[mapid1] points to global map array for mapid1 */
  Cmap **qmap;/**< query maps : only needed for debugging */
  int numrmaps;/**< mapid1 = roffset..roffset+numrmaps-1 */
  int numqmaps;/**< mapid2 = 0..numqmaps-1 */

  inline long long KEYHASH1(long long key) { 
    if(QUICKHASH >= 2){// shuffle bits in key so 243 (or 729) keys differing in size mostly differ in the lower order bits, so only 2(or 8) cache lines are needed in most cases 
      if(HASHWIN==5)/* shuffle 30-bit key so low order 2-bits of each 6-bit size are the 10 lowest order bits */
	return (key & (MASK(2) | (MASK(4) << 26))) | ((key >> 4) & (MASK(2)<<2)) | ((key >> 8) & (MASK(2)<<4)) | ((key >> 12) & (MASK(2)<<6)) | ((key >> 16) & (MASK(2)<<8))
	  | ((key << 8) & (MASK(4) << 10)) | ((key << 6) & (MASK(4) << 14)) | ((key << 4) & (MASK(4) << 18)) | ((key << 2) & (MASK(4) << 22));
      else // HASHWIN==6 : shuffle 36-bit key so lower order 2-bits of each 6-bit size are the 12 lowest order bits */
	return (key & (MASK_LL(2) | (MASK_LL(4) << 32))) | ((key >> 4) & (MASK_LL(2)<<2)) | ((key >> 8) & (MASK_LL(2)<<4)) | ((key >> 12) & (MASK_LL(2)<<6)) | ((key >> 16) & (MASK_LL(2)<<8)) | ((key >> 20) & (MASK_LL(2)<<10))
	  | ((key << 10) & (MASK_LL(4) << 12)) | ((key << 8) & (MASK_LL(4) << 16)) | ((key << 6) & (MASK_LL(4) << 20)) | ((key << 4) & (MASK_LL(4) << 24)) | ((key << 2) & (MASK_LL(4) << 28));
    } else {// QUICKHASH==1 
      if(HASHWIN==5) // low order 5 bits of each 6-bit size value
	return (key & MASK(5)) | ((key >> 1) & (MASK(5) << 5)) | ((key >> 2) & (MASK(5) << 10)) | ((key >> 3) & (MASK(5) << 15)) | ((key >> 4) & (MASK(5) << 20));
      else // HASHWIN==6 : low order 4 bits of each 6-bit size value 
	return (key & MASK(4)) | ((key >> 2) & (MASK(4) << 4)) | ((key >> 4) & (MASK(4) << 8)) | ((key >> 6) & (MASK(4) << 12)) | ((key >> 8) & (MASK(4) << 16)) | ((key >> 10) & (MASK(4) << 20));
    }
  };
  inline long long KEYHASH2(long long key) {
    if(HASHWIN==5) // high order 5 bits of each 6-bit size value
      return ((key >> 1) & MASK(5)) | ((key >> 2) & (MASK(5) << 5)) | ((key >> 3) & (MASK(5) << 10)) | ((key >> 4) & (MASK(5) << 15)) | ((key >> 5) & (MASK(5) << 20));
    else // middle order 4 bits of each 6-bit size value
      return ((key >> 1) & MASK(4)) | ((key >> 3) & (MASK(4) << 4)) | ((key >> 5) & (MASK(4) << 8)) | ((key >> 7) & (MASK(4) << 12)) | ((key >> 9) & (MASK(4) << 16)) | ((key >> 11) & (MASK(4) << 20));
  };
  inline long long KEYHASH3(long long key){// HASHWIN==6 : high order 4 bits of each 6-bit size value
    return ((key >> 2) & MASK(4)) | ((key >> 4) & (MASK(4) << 4)) | ((key >> 6) & (MASK(4) << 8)) | ((key >> 8) & (MASK(4) << 12)) | ((key >> 10) & (MASK(4) << 16)) | ((key >> 12) & (MASK(4) << 20));    
  };

  unsigned long long *keytable1;/* keytable1[HASHWIN==6 ? (1<<18) : (1<<19)] : a boolean (bit) table indexed by [KEYHASH1(sizekey) & MASK(HASHWIN==6 ? 24 : 25)] to detect presence in hashtable */
  unsigned long long *keytable2;/* keytable2[HASHWIN==6 ? (1<<18) : (1<<19)] : a boolean (bit) table indexed by [KEYHASH2(sizekey) & MASK(HASHWIN==6 ? 24 : 25)] to detect presence in hashtable */
  unsigned long long *keytable3;/* keytable3[1<<18] : a boolean (bit) table indexed by [KEYHASH3(sizekey) & MASK(24)] : only used if HASHWIN==6 */

  unsigned int *rand1;/**< rand1[0..MASK(RHASH_BITS)] is a table of random numbers used by hashindex() to improve the hash mapping */
  unsigned int *rand2;/**< rand2[0..MASK(RHASH_BITS)] is a table of random numbers used by hashindex() to improve the hash mapping */

  unsigned int *index;/**< index[hashindex(key)] is 0 OR is the start index into CollideBlock[]  :
					* The first entry gives the length of the array (fields key and hash are not valid), 
					* The actual collision arrays start with the next entry */
  CHashCollide *CollideBlock;/**< block allocated CHashCollide memory array */
  CHashEntry<OFFSET1_T> *EntryBlock;/**< block allocated CHashEntry memory array for all index[] entries combined */
  unsigned int CollideBlockMax,EntryBlockMax;/* Allocated size of CollideBlock[],EntryBlock[] respectively */
  unsigned int CollideBlockCnt,EntryBlockCnt;/* index into next free entry in CollideBlock[],EntryBlock[] respectively */

  CHashIndexR *indexR;
  CHashCollideR *CollideMemR;/**< compacted CHashCollideR memory array for all indexR[] entries combined */
  CHashEntry<OFFSET1_T> *EntryMem;/**< compacted CHashEntry memory array for all indexR[] entries combined */
  unsigned int CollideCnt, EntryCnt, CollideCntMax;
  unsigned int keytable_msk;// MASK(HASHWIN==6 ? 18 : 19)
  unsigned int rhash_msk; // MASK(RHASH_BITS)
  int sorted;/**< If hashtable collission lists have been compacted & sorted */
  int cloned;/**< If tables rand1[],rand2[],MapHash1[],MapHash2[],SiteHash[],OffsetHash[] point to another hashtable from which this hashtable was initialized */

  unsigned int *MapHash1, *MapHash2, *SiteHash, *OffsetHash, *Offset1Hash, *scaleIDhash;/**< random number Tables used by MatchTable and HitTable */

  Resource<OFFSET_T, OFFSET1_T, HASHWIN, RANGE, HASH_SCALE> *MatchMemReserve;

  /* member functions */
  inline unsigned int CollideBlockAlloc(unsigned int size);/**< Fast allocation of CHashCollide array of specified size as sub-array of CollideBlock[0..CollideBlockMax-1] */
  inline unsigned int EntryBlockAlloc(unsigned int size);/**< Fast allocation of CHashEntry array of specified size as sub-array of EntryBlock[0..EntryBlockMax-1] */

  inline unsigned int hashindex(long long key);/**< compute a hash index of more than HASH_BITS bits from the larger number of bits in key */
  inline void keyinsert(long long key, PFLOAT *win1, int mapid1, OFFSET1_T offset1, OFFSET1_T Roffset1, int site, int err);/**< insert new hashtable entry into hashtable */
  void hash_insert(int mapid1, Cmap *pmap);/**< insert single map into hashtable */
  void hashsort();/**< compact and sort hashtable entries for faster retrieval */
  void hashmerge(CHashTable **hashtables, int numhashtables);/**< merge multiple hashtables[0..numhashtables-1] into capacted hashtable. Note that "this" hashtable can be one of the hashtables[] */
  void hashinsert(Cmap **rmap, int numrmaps);/**< insert array of maps into hashtable */

  inline int hashfind(long long key, CHashEntry<OFFSET1_T>* &hash, CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *matchtable);/* retrieve pointer to hashtable entires for specified key : number of entries is returned (0 if none found) */
#if NEWHASH>=2
  int hashfindV(long long key, CHashResult *hash, CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *matchtable);/* retrieve index values for all hashtable entries matching any keys == key + SizeDelta[i=0..SizeLen-1] : length of hash[] array populated is returned */
#endif
  
  CHashTable(Cmap **Rmap, int numRmaps, Cmap **Qmap, int numQmaps, CHashTable *clone) {
    sorted = 0;

    /* The following initial allocation size will be reduced if MaxMem > 0 AND/OR HashInsertThreads > 1 */
    CollideBlockMax = (1U << min(31,HASH_BITS+3)) - 1;
    EntryBlockMax = (1U << min(31,HASH_BITS+5)) - 1;
    if(VERB>=2){
      printf("hashtable=%p:HASH_BITS=%d,CollideBlockMax=%u,EntryBlockMax=%u,(1<<min(31,HASH_BITS+3))=%u, index= %p\n",this,HASH_BITS,CollideBlockMax,EntryBlockMax,(1U << min(31,HASH_BITS+3)),index);
      fflush(stdout);
    }
    size_t finalmem1 = (1U << HASH_BITS)*sizeof(CHashIndexR) + TotWin * sizeof(CHashEntry<OFFSET1_T>);// estimate of compacted hashtable size (estimated VmRSS) : ignores CHashCollideR[] overflow
    size_t finalmem2 = (1U << HASH_BITS)*sizeof(unsigned int);// estimate of per thread hashtable size : ignores CHashCollideR[] overflow
    size_t finalmem = finalmem1 + finalmem2 * HashInsertThreads;
    if(QUICKHASH && hashkeys){
      if(QUICKHASH==1)
	finalmem += (HASHWIN==5 ? 2 : 3) * (1<<18) * sizeof(unsigned long long);// estimate of memory for keytable1[],keytable2[] and keytable3[]
      if(QUICKHASH==2)
	finalmem += (1 << (HASHWIN==5 ? 24 : 30)) * sizeof(unsigned long long);// estimate of memory for keytable1[] with full indexing
    }
    size_t maxmem = (size_t)((HMaxMem > 0 ? HMaxMem : 256.0) * 1e9);
    size_t totmem = (CollideBlockMax*sizeof(CHashCollide) + EntryBlockMax*sizeof(CHashEntry<OFFSET1_T>)) * HashInsertThreads;
    if(VERB>=2){
      size_t entrymem = TotWin * sizeof(CHashEntry<OFFSET1_T>);
      printf("hashtable=%p:maxmem= %0.4f, totmem= %0.4f, finalmem= %0.4f GB (%0.4f GB + %0.4f GB /thread), HashInsertThreads=%d,clone=%p,TotWin=%lld, entrymem= %0.4f Gb\n",
	     this,maxmem*1e-9,totmem*1e-9,finalmem*1e-9,finalmem1*1e-9,finalmem2*1e-9,HashInsertThreads,clone,(long long)TotWin, entrymem * 1e-9);
      fflush(stdout);
    }
    if(totmem + finalmem > maxmem){
      if(!clone && finalmem > maxmem * 0.6){
	if(finalmem > maxmem * 0.6 && HashInsertThreads > 4){/* reduce HashInsertThreads down to 4 threads, if needed */
	  int origHashInsertThreads = HashInsertThreads;

	  HashInsertThreads = 4;
	  if(maxmem * 0.6 > finalmem1 + 4*finalmem2){
	    HashInsertThreads = (maxmem * 0.6 - finalmem1)/finalmem2;
	    if(DEBUG) assert(HashInsertThreads >= 4);
	  }

	  if(VERB){
	    printf("Reducing HashInsertThreads from %d to %d to reduce hashtable Index VmRSS from %0.4f Gb to %0.4f Gb (maxmem = %0.4f Gb, hashindex = %0.4f GB + %0.4f GB/thread)\n",
		   origHashInsertThreads,HashInsertThreads, finalmem * 1e-9, (finalmem1 + finalmem2 * HashInsertThreads) * 1e-9, maxmem * 1e-9, finalmem1 * 1e-9, finalmem2 * 1e-9);
	    fflush(stdout);
	  }
	  
	  totmem = (CollideBlockMax*sizeof(CHashCollide) + EntryBlockMax*sizeof(CHashEntry<OFFSET1_T>)) * HashInsertThreads;
	  finalmem = finalmem1 + finalmem2 * HashInsertThreads;
	}

	/* reduce HASH_BITS (by 1 or 2), if finalmem is more than 60% of maxmem */        
	int origHASH_BITS = HASH_BITS;
	size_t origfinalmem = finalmem;

	while(finalmem > maxmem * 0.6){
	  HASH_BITS--;
	  finalmem1 = (1U << HASH_BITS)*sizeof(CHashIndexR) + TotWin * sizeof(CHashEntry<OFFSET1_T>);// estimate of compacted hashtable size (estimated VmRSS) : ignores CHashCollideR[] overflow
	  finalmem2 = (1U << HASH_BITS)*sizeof(unsigned int);// estimate of per thread hashtable size : ignores CHashCollideR[] overflow
	  finalmem = finalmem1 + finalmem2 * HashInsertThreads;

	  CollideBlockMax = (1U << min(31,HASH_BITS+3)) - 1;
	  EntryBlockMax = (1U << min(31,HASH_BITS+5)) - 1;
	  totmem = (CollideBlockMax*sizeof(CHashCollide) + EntryBlockMax*sizeof(CHashEntry<OFFSET1_T>)) * HashInsertThreads;	  

	  if(HASH_BITS <= origHASH_BITS - 2)
	    break;
	}

	if(HASH_BITS < origHASH_BITS){
	  if(VERB){
	    size_t entrymem = TotWin * sizeof(CHashEntry<OFFSET1_T>);
	    printf("Reducing HASH_BITS from %d to %d to reduce hashtable Index VmRSS from %0.4f Gb to %0.4f Gb (maxmem= %0.4f Gb, entrymem= %0.4f Gb, TotWin=%lld)\n",
		   origHASH_BITS, HASH_BITS, (origfinalmem - entrymem) * 1e-9, (finalmem - entrymem) *1e-9, maxmem * 1e-9, entrymem * 1e-9, (long long)TotWin);
	    fflush(stdout);
	  }
	}

	if(finalmem > maxmem * 0.6 && HashInsertThreads > 1){/* reduce HashInsertThreads to 1 */
	  int origHashInsertThreads = HashInsertThreads;
	  HashInsertThreads = 1;
	  if(maxmem * 0.6 > finalmem1 + finalmem2){
	    HashInsertThreads = (maxmem * 0.6 - finalmem1)/finalmem2;
	    if(DEBUG) assert(HashInsertThreads >= 1);
	  }

	  if(VERB){
	    size_t entrymem = TotWin * sizeof(CHashEntry<OFFSET1_T>);
	    printf("Reducing HashInsertThreads from %d to %d to reduce hashtable Index VmRSS from %0.4f Gb to %0.4f Gb ( + entrymem= %0.4f Gb)\n",
		   origHashInsertThreads, HashInsertThreads, (finalmem - entrymem) * 1e-9, (finalmem1 + finalmem2 * HashInsertThreads - entrymem) * 1e-9, entrymem * 1e-9);
	    fflush(stdout);
	  }
	  
	  totmem = (CollideBlockMax*sizeof(CHashCollide) + EntryBlockMax*sizeof(CHashEntry<OFFSET1_T>)) * HashInsertThreads;
	  finalmem = finalmem1 + finalmem2 * HashInsertThreads;
	}
      }


      //      double ratio = ((double)(maxmem - finalmem))/((double)totmem);
      size_t origtotmem = totmem;
      unsigned int origCollideBlockMax = CollideBlockMax;
      unsigned int origEntryBlockMax = EntryBlockMax;
      while(totmem + finalmem > maxmem && totmem * 8 > maxmem && CollideBlockMax  > 255){
	CollideBlockMax /= 2;
	EntryBlockMax /= 2;
	totmem = (CollideBlockMax*sizeof(CHashCollide) + EntryBlockMax*sizeof(CHashEntry<OFFSET1_T>)) * HashInsertThreads;	
      }

      if(VERB/* HERE >=2 */ && !clone){
	printf("hashtable=%p:CollideBlockMax=%u->%u,EntryBlockMax=%u->%u (per thread), HashInsertThreads=%d, TotWin=%lld, maxmem=%0.4f Gb, InitialMem=%0.4f -> %0.4f Gb\n",
	       this,origCollideBlockMax,CollideBlockMax,origEntryBlockMax,EntryBlockMax, HashInsertThreads, (long long)TotWin, maxmem * 1e-9, 
	       (origtotmem + finalmem) * 1e-9, (totmem + finalmem) * 1e-9);
	fflush(stdout);
      }
    }

    index = new unsigned int[1U << HASH_BITS];
    if(clone){
      memset(index, 0, (1U << HASH_BITS)*sizeof(unsigned int));
    } else {
      #pragma omp parallel for num_threads(hnumthreads) schedule(static,1)
      for(int hindex = MASK(HASH_BITS-9); hindex >= 0; hindex--)
        memset(&index[hindex << 9], 0, (1 << 9)*sizeof(unsigned int));
    }

    POSIX_MEMALIGN(CollideBlock,64,sizeof(CHashCollide)*CollideBlockMax);
    POSIX_MEMALIGN(EntryBlock,64,sizeof(CHashEntry<OFFSET1_T>)*EntryBlockMax);
    CollideBlockCnt = EntryBlockCnt = 1;/* Avoid 0 which since this value is used to denote no entries in index[] */

    indexR = 0;
    CollideMemR = 0;
    EntryMem = 0;

    rmap = Rmap;
    qmap = Qmap;
    numrmaps = numRmaps;
    numqmaps = numQmaps;

    if(QUICKHASH && hashkeys){
      if(QUICKHASH==1)
	keytable_msk = MASK(HASHWIN==6 ? 18 : 19);
      if(QUICKHASH==2)
	keytable_msk = MASK(HASHWIN==6 ? 30 : 24);
    }
    rhash_msk = MASK(RHASH_BITS);

    if(!clone){
      cloned = 0;

      srandom(1);

      if(QUICKHASH && hashkeys){
        if(QUICKHASH==1){
          //          keytable1 = new unsigned long long[HASHWIN==6 ? 3 * (1 << 18) : 2 * (1 << 19)];
          POSIX_MEMALIGN(keytable1, 64, sizeof(unsigned long long) * (HASHWIN==6 ? 3 * (1 << 18) : 2 * (1 << 19)));
	  keytable2 = &keytable1[HASHWIN==6 ? (1 << 18) : (1 << 19)];
	  if(HASHWIN==6)
	    keytable3 = &keytable1[2*(1 << 18)];
        }
        if(QUICKHASH==2)
          POSIX_MEMALIGN(keytable1, 64, sizeof(unsigned long long) * (1 << (HASHWIN==6 ? 30 : 24)));
      }

      rand1 = new unsigned int[1 << (RHASH_BITS+1)];
      rand2 = &rand1[1 << RHASH_BITS];

      register unsigned int HashMask = MASK(HASH_BITS);
      for(int i = (1 << RHASH_BITS); --i >= 0;){
	rand1[i] = HashMask & random();
	rand2[i] = HashMask & random();
      }
      if(VERB>=2){
	for(int i = 0; i < 16; i++)
	  printf("i=%d:rand1[i]=%u,rand2[i]=%u\n",i,rand1[i],rand2[i]);
	fflush(stdout);
      }
    
      /* allocate/initialize MapHash1[],MapHash2[],SiteHash[],OffsetHash[] */
      MapHash1 = new unsigned int[numrmaps*2];
      MapHash2 = &MapHash1[numrmaps];

      SiteHash = new unsigned int[(MaxOffset-MinOffset+1)*2 + NumScaleFactor];
      OffsetHash = &SiteHash[MaxOffset-MinOffset+1];
      scaleIDhash = &SiteHash[(MaxOffset-MinOffset+1)*2];

      SiteHash -= MinOffset; /* SiteHash[MinOffset .. MaxOffset] */
      OffsetHash -= MinOffset; /* OffsetHash[MinOffset .. MaxOffset] */

      Offset1Hash = NULL;
      if(RANGE >= 1) // Also allocate Offset1Hash[]
	Offset1Hash = new unsigned int[1<<16];

      unsigned int *thash1 = new unsigned int[(1 << 16)*4];
      unsigned int *thash2 = &thash1[1 << 16];
      unsigned int *thash3 = &thash1[2*(1 << 16)];
      unsigned int *thash4 = &thash1[3*(1 << 16)];

      srandom(0);
      unsigned int HitMask = MASK(HHASH_BITS);
      unsigned int MatchMask = MASK(MHASH_BITS);
      for(int i = (1 << 16); --i >= 0;){
	thash1[i] = random();
	thash2[i] = random();
	thash3[i] = random();
	thash4[i] = random();
	if(RANGE >= 1)
	  Offset1Hash[i] = random() & MatchMask;
      }
      for(int i = numrmaps; --i >= 0;){
	long long rid = rmap[roffset + i]->id;
	if(DEBUG && !(rid >= 0)){
          printf("CHashTable: i=%d/%d:roffset=%d,rmap[roffset+i]->id= %lld, mapid= %d\n",i,numrmaps,roffset,rmap[roffset+i]->id,rmap[roffset+i]->mapid);
	  fflush(stdout);
	  assert(rid >= 0);
	}
	unsigned int msk = thash1[rid & MASK(16)] ^ thash2[(rid>>16) & MASK(16)] ^ thash3[(rid>>32) & MASK(16)] ^ thash4[(rid>>48) & MASK(16)];
	MapHash1[i] = msk & HitMask;
	MapHash2[i] = msk & MatchMask;
      }
      delete [] thash1;
      if(VERB>=2 || TRACE){
	for(int i = 0; i < numrmaps; i++)
	  if(VERB>=2 || (TRACE && rmap[roffset+i]->id == TraceID1))
	    printf("i=%d,id=%lld:MapHash1[i]=%u,MapHash2[i]=%u\n",i, rmap[roffset+i]->id, MapHash1[i],MapHash2[i]);
	fflush(stdout);
      }
      for(int i = 0; i <= -MinOffset || i <= MaxOffset; i++){
	unsigned int msk = random();
	if(i <= MaxOffset){
	  SiteHash[i] = msk & HitMask;
	  OffsetHash[i] = msk & MatchMask;
	}
	msk = random();
	if(i && i <= -MinOffset){
	  SiteHash[-i] = msk & HitMask;
	  OffsetHash[-i] = msk & MatchMask;
	}
      }
      if(HASH_SCALE){
	scaleIDhash[0] = 0;
        for(int i = 1; i < NumScaleFactor; i++)
	  scaleIDhash[i] = random() & MatchMask;
      }

      if(VERB>=2){
	for(int i = MinOffset; i <= MaxOffset; i++)
	  printf("i=%d:OffsetHash[i]=%u,SiteHash[i]=%u\n",
		 i,OffsetHash[i],SiteHash[i]);
	fflush(stdout);
      }
    } else {/* clone != 0 */
      cloned = 1;
      if(QUICKHASH && hashkeys){
        keytable1 = clone->keytable1;
	if(QUICKHASH==1){
          keytable2 = clone->keytable2;
          if(HASHWIN==6)
	    keytable3 = clone->keytable3;
        }
      }
      rand1 = clone->rand1;
      rand2 = clone->rand2;
      MapHash1 = clone->MapHash1;
      MapHash2 = clone->MapHash2;
      SiteHash = clone->SiteHash;
      OffsetHash = clone->OffsetHash;
      if(RANGE >= 1)
	Offset1Hash = clone->Offset1Hash;
    }
    if(VERB>=2){
      if(!cloned)
	printf("Allocated hashtable=%p:index=%p,indexR=%p,CollideMemR=%p,EntryMem=%p,rand1=%p,MapHash1=%p\n",this,index,indexR,CollideMemR,EntryMem,rand1,MapHash1);
      else
	printf("Allocated hashtable=%p:index=%p,indexR=%p,CollideMemR=%p,EntryMem=%p\n",this,index,indexR,CollideMemR,EntryMem);
      printf("\t index[0]=%u,index[1]=%u,index[%u]=%u\n",index[0],index[1],(1U << HASH_BITS)-1U,index[(1U << HASH_BITS)-1U]);
      fflush(stdout);
    }
  };

  ~CHashTable() {
    if(VERB>=2){
      if(!cloned)
	printf("deleting hashtable=%p:index=%p,indexR=%p,CollideMemR=%p,EntryMem=%p,rand1=%p,MapHash1=%p\n",this,index,indexR,CollideMemR,EntryMem,rand1,MapHash1);
      else
	printf("deleting hashtable=%p:index=%p,indexR=%p,CollideMemR=%p,EntryMem=%p\n",this,index,indexR,CollideMemR,EntryMem);
      fflush(stdout);
      printf("\t index[0]=%u,index[1]=%u,index[%u]=%u\n",index[0],index[1],(1U << HASH_BITS)-1U,index[(1U << HASH_BITS)-1U]);
      fflush(stdout);
    }
    //      mcheck_check_all();

    if(!sorted){
      if(CollideBlock) free(CollideBlock);
      if(EntryBlock) free(EntryBlock);
    }

    delete [] index;

    if(indexR) free(indexR);
    if(CollideMemR) free(CollideMemR);
    if(EntryMem) free(EntryMem);

    if(!cloned){
      if(QUICKHASH && hashkeys && keytable1)
	free(keytable1);
      delete [] rand1;
      delete [] MapHash1;
      if(SiteHash) delete [] &SiteHash[MinOffset];
      if(RANGE >= 1){
	if(Offset1Hash) delete [] Offset1Hash;
      }
    }
    //      mcheck_check_all();
  };

  size_t HashSizeRtab() { /* return total size of randum number tables in bytes */
    size_t cnt = 0;
    cnt += 2LL * (1 << RHASH_BITS)*sizeof(unsigned int);
    cnt += 2LL * numrmaps * sizeof(unsigned int);
    cnt += 2LL * (MaxOffset-MinOffset+1)*sizeof(unsigned int);
    return cnt;
  }

  size_t HashSize(){/* return total size of all tables in bytes */
    if(!sorted)
      return 0;/* not implemented */
    size_t cnt = 0;
    cnt += (1U << HASH_BITS)*sizeof(CHashIndexR);
    cnt += CollideCnt * sizeof(CHashCollideR);
    cnt += EntryCnt*sizeof(CHashEntry<OFFSET1_T>);
    cnt += HashSizeRtab();
    return cnt;
  };
};

/** Describes matching window (mapid1,offset1) and links to neighboring CHitEntry in query map. Initialized and updated by HitInsert() */
template<class OFFSET1_T>
class CHitEntry {
 public:
  int next;/* index into HitMem[] of next CHitEntry with same (Maphash1[mapid1] ^ Sitehash[offset1]) */
  int mapid1;/* IF TRACE : the high order 3 bits encode the number of misaligned + misresolved labels (site2-nsite2) (NOT YET IMPLEMENTED) */
  OFFSET1_T offset1; /* flip2 ? floor((Ylen[mapid1]-Y[mapid1][site1+HashWin+CNT(err1)])/OffsetKB + 0.5) 
		              : floor(Y[mapid1][site]/OffsetKB + 0.5) */
  void init(int id1, OFFSET1_T off1, int s, int n) {
    mapid1 = id1;
    offset1 = off1;
    next = n;
  };
};

// #define CNT(err) (err ? 1 : 0)

/** Highest scoring hash result for a given ref map. Updated by BestInsert() */
/* Also used as temporary scratch space MatchBuf[] by MatchFind() */
template<class OFFSET_T, int RANGE, int HASH_SCALE>
class CBestMatch {
 public:
  int mapid1;
  OFFSET_T offset;
  short score;
  unsigned short MinOffset1[RANGE];
  unsigned short MaxOffset1[RANGE];// Also used in MatchFind() to store original score value (before -hashdelta adjustment), so it can be restored in MatchTable for entries that are not deleted
  unsigned char scaleID[HASH_SCALE];
  /* NOTE: for more accurate offset estimation with OFFSET_SPREAD > 0 : compute weighted average of offset*OffsetKB over 1+2*OFFSET_SPREAD values with wt= OFFSET_WT[].
     Requires two additional floats in the struct : offsetSum, wtSum. This is only needed if offset will be used subsequently */
};

/* Note: There will be once instance of CMatchTable for each main thread (each main thread works on a seperate query (probe) map) */

/** Data for hash hits, storage and validation methods
* Note: There will be once instance of CMatchTable for each thread
* Note: The query map (mapid2, flip2, X[]) is fixed during the use of the MatchTable until MatchFind() is used to reset the MatchTable for the next query */

template <class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
class CMatchTable {
  unsigned int *Maphash1, *Maphash2, *Sitehash, *Offsethash, *Offset1hash, *scaleIDhash; /* pointers to correct copy (will be in local memory if -localmem was used, otherwise all threads point to same copy) */

 public:
  int HitFirst;/* index to first non-zero entry in HitIndex[] */
  int *HitIndex;/* HitIndex[MapHash1[mapid1] ^ SiteHash[offset1]] is 0 OR an index into HitMem[] of type CHitEntry */
  CHitEntry<OFFSET1_T> *HitMem;/* Memory pool HitMem[0..HitMax-1] : may grow over time */
  unsigned int HitMax;
  unsigned int HitFree;/* Index into next free location in HitMem[] */
  inline int CHitEntryAlloc();/* grab an entry from HitMem[HitFree..HitMax-1] */
  int *HitNext;/* If HitIndex[i] is nonzero, then HitNext[i] is the index of the next nonzero entry in HitIndex (or -1 if None) */

  long long *sizekeyV;/* sizekeyV[S = 0..SizeLen-1] is used to prefetch (and filter, if QUICKHASH) key + SizeDelta[S] */
  unsigned int *rand1V;/* rand1V[S = 0..SizeLen-1] is used to prefetch rand1[(key + SizeDelta[S]) & rhash_msk] */
  //  CHashIndexR *HashIndexV;/* HashIndexV[S = 0..SizeLen-1] is used to prefetch indexR[rand1V[S]] */
  CHashResult *hashV;/**< hashV[0..SizeLen-1] is used for Software pipelining and the results of hashtable lookups by hashfindV() */
  
  int MatchFirst;/* index to the first entry in MatchIndex[] with cnt > 0 (or -1 if none) */
  CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *MatchIndex;/* MatchIndex[MapHash2[mapid1]^OffsetHash[offset]] has cnt==0 or includes the first MATCHLINE entries and the index to the Overflow linked list */
  CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *MatchOverflow;/* MatchOverflow[0..MatchOverflowMax-1] : memory for Overflow linked list : content is referenced by index, so can be reallocated */
  int MatchOverflowMax;/* Cannot exceed MASK(31) , due to field next */
  int MatchOverflowNext;/* Unused area MatchOverflow[MatchOverflowNext .. MatchOverflowMax-1] */
  int MatchOverflowFree;/* Start index of free memory linked list (using field "next") within MatchOverflow[MATCHLINE+1..MatchOverflowNext-1] (or -1 if empty) */
  long long MatchCnt;/**< total matches in MatchTable */
  int minSite2;/* the smallest site2 value whose matches need to be transfered to MatchBuf[] & BestMatch[]. Normally 0, unless HashGC > 0 : In that case matches for smaller site2 values are already in BestMatch[] */

  CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *MatchBuf;/**< MatchBuf[0..MatchBufMax-1] : linear array to buffer entries in MatchTable */
  long long MatchBufMax;

  CBestMatch<OFFSET_T,RANGE,HASH_SCALE> **BestMatch;/**< BestMatch[mapid1][0..BestMatchCnt[mapid1]-1] contains information about the best score(s) for mapid1 in Match Table */
  CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *BestMatchMem;/* single block of memory used by BestMatch[][] */
  int *BestMapid;/**< BestMapid[0..BestCnt-1] points to the mapid1 values with valid entries in BestMatch[] */
  int BestCnt;
  int *BestMatchCnt;/**< BestMatchCnt[mapid1] is the number of best scores for mapid1 in Match Table : Only used if HashMultiMatch != 0*/
  int *BestMatchMax;/**< BestMatchMax[mapid1] is the amount of memory currently allocated for BestMatch[mapid1][0..BestMatchMax-1] : Only used if HashMultiMatch != 0*/

  int tid;/* which thread is using this MatchTable */
  int pid;/* process id : used to generate unique output tmp file names */
  CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *hashtable;/* pointer to shared hashtable */
  CMatchTable *matchtable;/* If != NULL : pointer to Matchtable for flip2==0. Only used if 2 threads are used for each hash_query() call */
  
  int MTHashGC;/* max(HashGC, N2 * ScaleDelta * ScaleDeltaSiz) */

  /* output buffer : Only defined if matchtable == NULL */
  FILE *fp;
  char *fbuf;
  size_t matchcnt, matchmax;
  CHashMatch *matches;
  void OutputFlush();/* flush output buffer */

  unsigned int HitMemHWM;/**< High Water Mark for HitMax (size of HitMem[] used by Hit Table) for any single query */
  int MatchMemHWM;/**< High Water Mark for MatchOverflowMax (size of MatchOverflow[] used by Match Table) for any single query (HWMmapid) */
  int HWMmapid;
  long long MatchBufHWM;/**< High Water Mark for MatchBufMax (size of MatchBuf[]) for any single query */
  int BestCntHWM;/**< High Water Mark for BestCnt (size of BestMapid[]) for any single query */
  size_t OutputHWM;/**< High Water Mark for matchmax (sizeo of matches[]) for any single query */

  inline int MatchOverflowAlloc(int tid);/* return next index into MatchOverflow[] */
  inline short MatchScoreFind(int mapid1, OFFSET_T offset, unsigned short Offset1, unsigned char scaleID);/* return score (or 0 if not found) */
  inline void MatchScoreUpdate(int mapid1, OFFSET_T offset, unsigned short Offset1, unsigned char scaleID, unsigned char score, int mapid2, int flip2);/* Update score (or Error if not found) */
  inline void BestInsert(int mapid1, OFFSET_T offset, short score, CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *pBuf, int mapid2,int flip2, int tid);/* update Best Table (last 3 args for debugging) */

  void hash_queryF(Cmap *pmap, int mapid2, int flip2, FLOAT *X, FLOAT *Xscale, int *Xflag, int N, int tid);/* generate matches for one map in orientation flip2 and append them to output buffer (or FILE MatrixFP if HashMatrix) */
  void keyquery(long long key, PFLOAT *win2, int mapid2, int flip2, FLOAT *X, int site2, unsigned int err2, OFFSET_T offset2, int scaleID);/* check all sizing variations of key and append hashtable matches to Match Table after removing duplicates */

  /* The following two functions are only called from keyquery(long long, ... ) : */
  inline void keyqueryB(const int DISJOINT, int hcnt, int mapid2, int flip2, OFFSET_T offset2,
#if USE_HSSE==1 && USE_PFLOAT==1 && USE_SSE==1 && EXACT_SIZING && MAXR_WIN==8 && (HASHWIN == 5 || HASHWIN == 6)
				    __m128 vL_LB , __m128 vL_UB , __m128 vL_varInv , __m128 vL_win2,
				    __m128 vH_LB , __m128 vH_UB , __m128 vH_varInv , __m128 vH_win2,  __m128 v_SDnormMax,
#endif
				    CHashEntry<OFFSET1_T> *hashBase, PFLOAT *win2, PFLOAT *LB, PFLOAT *UB, PFLOAT *varInv);

  inline void keyentry_hashloop(const int DISJOINT,const int PREFETCH,int i, int mapid2, int flip2, OFFSET_T offset2,
#if USE_HSSE==1 && USE_PFLOAT==1 && USE_SSE==1 && EXACT_SIZING && MAXR_WIN==8 && (HASHWIN == 5 || HASHWIN == 6)
				__m128 vL_LB , __m128 vL_UB , __m128 vL_varInv , __m128 vL_win2,
				__m128 vH_LB , __m128 vH_UB , __m128 vH_varInv , __m128 vH_win2, __m128 v_SDnormMax,
#endif
				CHashEntry<OFFSET1_T> *hashBase, PFLOAT *win2, PFLOAT *LB, PFLOAT *UB, PFLOAT *varInv);

  inline void HitInsert(int mapid1, OFFSET1_T offset1, OFFSET_T offset2, int site1,int TraceFlag);/* insert hit into Hit Table, removing duplicates : 
					     All calls share the same (mapid2,flip2,site2) until HitToMatch() is called to reset the Hit Table (and free all memory used) */
  void HitToMatch(int mapid2, int flip2, FLOAT *X, int N, int site2, OFFSET_T offset2, int scaleID, int tid);/* Copy all entries in Hit Table to Match Table, and reset Hit Table
									NOTE : If HashMatrix just output all remaining matches in Hit Table to FILE MatrixFP */

  inline unsigned char MatchInsert(int mapid1, OFFSET_T offset, unsigned char score, int mapid2, int flip2, OFFSET_T offset2, unsigned char scaleID, unsigned short Offset1, OFFSET1_T offset1);
                                                   /* add hit to Match Table (mapid2,flip2,offset2,offset1 for debugging only) : Returns updated hashscore. score is always 1, unless OFFSET_ZEROFIX */
  void MatchFind(int mapid2, int flip2, int maxSite2,int tid);/* Locate the best remaining matches (for each mapid1) in Match Table, apply OFFSET_SPREAD then append those with score > HashThresh to output buffer */ 

  long long totalMcnt, totalQcnt;/* progress counters */
  long long HTindcnt, HTqindcnt, HTfindcnt, HTcolcnt, HTqcolcnt, HTmatchcnt, HTentrycnt;/* HashTable performance statistics */
  long long MTcnt, MTcnt2, MTcnt3, MTcol, MTqcol, MTicnt;/* MatchTable performance statistics */
  long long HTcnt, HTcol, HTdup;/* HitTable performance statistics */

  /* input map information */
  int numrmaps;
  Cmap **rmap; /* rmap[mapid1 = roffset .. roffset+numrmaps-1] */  
  int numqmaps;
  Cmap **qmap; /* qmap[mapid2 = 0 .. numqmaps-1] */

  char *output_prefix;

  char filename[PATH_MAX];/* output buffer filename */

  CMatchTable(int thread_id, CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *Hashtable, int numRmaps, Cmap **Rmap, int numQmaps, Cmap **Qmap, char *prefix, CMatchTable *Matchtable){
    if(STAT) HTindcnt = HTqindcnt = HTfindcnt = HTcolcnt = HTqcolcnt = HTmatchcnt = HTentrycnt = 0;
    if(STAT) MTcnt = MTcnt2 = MTcnt3 = MTcol = MTqcol = MTicnt = HTcnt = HTcol = HTdup = 0;

    tid = thread_id;
    pid = getpid();
    hashtable = Hashtable;
    matchtable = Matchtable;
    
    /* save input map information */
    numrmaps = numRmaps;
    rmap = Rmap;
    numqmaps = numQmaps;
    qmap = Qmap;

    output_prefix = prefix;

    /* initialize local copy of pointer to random number tables MapHash1[],MapHash2[],SiteHash[],OffSetHash[] */
    Maphash1 = Hashtable->MapHash1 - roffset;
    Maphash2 = Hashtable->MapHash2 - roffset;
    Sitehash = Hashtable->SiteHash;
    Offsethash = Hashtable->OffsetHash;
    Offset1hash = Hashtable->Offset1Hash;
    scaleIDhash = Hashtable->scaleIDhash;

    /* allocate rand1V,HashIndexV,hashV buffers for Hashtable lookups */
    int SizeLenRnd = (SizeLen  + 15) & ~15;/* round SizeLen to nearest multiple of 16 */
    sizekeyV = new long long[SizeLenRnd];
    rand1V = new unsigned int[SizeLenRnd];
    //    HashIndexV = new CHashIndexR[SizeLenRnd];
    hashV = new CHashResult[SizeLenRnd];

    /* initialize output buffer */
    matchcnt = 0;
    matchmax = 4*1024;/* 64 kbytes output buffer */
    matches = new CHashMatch[matchmax];

    fbuf = 0;

    MTHashGC = HashGC;/* will be updated in hash_query() */

    if(!matchtable){
      sprintf(filename,"%s_%d_%d.bin",prefix,pid,tid);
      if(DEBUG && strlen(filename) >= PATH_MAX){
	printf("filename = %s is too long (%lu chars), must not exceed %d chars\n",filename, strlen(filename),PATH_MAX-1);
	fflush(stdout);exit(1);
      }

      unlink(filename);/* in case the file already exists : could be from failed run of same command or (unlikely) same pid on different machines */
      errno = 0;
      if((fp = fopen(filename,"w"))==NULL){
	int eno = errno;
#pragma omp critical
	{
	  char *err = strerror(eno);
	  printf("failed to open tmp file %s for writing:errno=%d:%s\n", filename,eno,err);
	  fflush(stdout); exit(1);
	}
      }
      if(VERB>=2){
#pragma omp critical
	{
	  printf("pid=%d:Opened tmp file %s for writing\n",pid,filename);
	  fflush(stdout);
	}
      }
      fbuf = new char[sizeof(CHashMatch)*matchmax];
      setbuffer(fp,fbuf, matchmax * sizeof(CHashMatch));
    }
    totalMcnt = totalQcnt = 0;

    HitFirst = -1;
    HitIndex = new int[1 << HHASH_BITS];
    memset(HitIndex, 0, (1 << HHASH_BITS)*sizeof(int));
    HitMax = (1 << (HHASH_BITS-2));
    POSIX_MEMALIGN(HitMem,64,sizeof(CHitEntry<OFFSET1_T>)*HitMax);
    HitFree = 1;/* Index 0 is reserved */
    HitNext = new int[1 << HHASH_BITS];
    memset(HitNext, 0, (1 << HHASH_BITS)*sizeof(int));

    minSite2 = 0;
    MatchFirst = -1;

    // WAS    POSIX_MEMALIGN(MatchIndex,64,sizeof(CMatchLine)*(1U << MHASH_BITS));
    size_t memsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * (1UL << MHASH_BITS) + PAGE-1) & ~(PAGE-1);
    char *newmem = (char *) mmap(NULL, memsiz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if(newmem == MAP_FAILED){
      int eno = errno;
      char *err = strerror(eno);
	    
      #pragma omp critical
      {
        printf("CMatchTable: mmap of %lu bytes failed:errno=%d:%s\n", memsiz,eno,err);
	dumpmemmap();
	fflush(stdout);exit(1);
      }
    }
    MatchIndex = (CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *)newmem;

    for(int i = MASK(MHASH_BITS); i >= 0;i--)
      MatchIndex[i].cnt = 0;
    MatchOverflowMax = (1 << (MHASH_BITS-2));
    MatchOverflowNext = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) + 1;/* Entries MatchOverflow[0..MATCHLINE] are never used */
    MatchOverflowFree = -1;

    // WAS    POSIX_MEMALIGN(MatchOverflow,64,sizeof(CMatchLine)*MatchOverflowMax);
    memsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * MatchOverflowMax + PAGE-1) & ~(PAGE-1);
    newmem = (char *) mmap(NULL, memsiz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if(newmem == MAP_FAILED){
      int eno = errno;
      char *err = strerror(eno);
	    
      #pragma omp critical
      {
        printf("CMatchTable: mmap of %lu bytes failed:errno=%d:%s\n", memsiz,eno,err);
	dumpmemmap();
	fflush(stdout);exit(1);
      }
    }
    MatchOverflow = (CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *)newmem;

    MatchCnt = 0;

    MatchBufMax = 0;
    MatchBuf = 0;

    BestMatch = new CBestMatch<OFFSET_T,RANGE,HASH_SCALE>*[numrmaps];
    BestMapid = new int[numrmaps];
    BestCnt = 0;
    if(!HashMultiMatch){
      BestMatchCnt = BestMatchMax = NULL;
      BestMatchMem = new CBestMatch<OFFSET_T,RANGE,HASH_SCALE>[numrmaps];
      for(int i = numrmaps; --i >= 0;){
	BestMatch[i] = &BestMatchMem[i];
	BestMatch[i]->mapid1 = -1;
      }
    } else {
      BestMatchCnt = new int[numrmaps];
      BestMatchMax = new int[numrmaps];
      BestMatchMem = NULL;
      for(int i = numrmaps; --i >= 0;){
	BestMatchCnt[i] = BestMatchMax[i] = 0;
	BestMatch[i] = NULL;
      }
    }
    if(splitcnt > 1){
      BestMatch -= roffset;
      if(HashMultiMatch){
        BestMatchCnt -= roffset;
	BestMatchMax -= roffset;
      }
    }

    if(VERB>=3){
      #pragma omp critical
      {
	printf("tid=%d:Allocated CMatchTable=%p:BestMatchCnt=%p,BestMatchMax=%p,BestMatchMem=%p,BestMatch=%p,BestMapid=%p,MatchOverflow=%p,MatchIndex=%p\n",
	       tid,this,BestMatchCnt,BestMatchMax,BestMatchMem,BestMatch,BestMapid,MatchOverflow,MatchIndex);
	fflush(stdout);
      }
    }

    HitMemHWM = 0;
    MatchMemHWM = 0;
    HWMmapid = -1;
    MatchBufHWM = 0;
    BestCntHWM = 0;
    OutputHWM = 0;

  };

  ~CMatchTable(){

    if(matchtable == NULL){
      OutputFlush();

      if(VERB>=3){
        #pragma omp critical
	{
  	  printf("tid=%d:Deallocating CMatchTable=%p:BestMatchCnt=%p,BestMatchMax=%p,BestMatchMem=%p,BestMatch=%p,BestMapid=%p,MatchOverflow=%p,MatchIndex=%p,fp=%p,filename=%s\n",
	       tid,this,BestMatchCnt,BestMatchMax,BestMatchMem,BestMatch,BestMapid,MatchOverflow,MatchIndex,fp,filename);
	  fflush(stdout);
        }
      }

      FILEclose(fp);

      if(VERB>=2){
        #pragma omp critical
        {
  	  printf("tid=%d:Closed tmp file %s\n",tid,filename);
	  fflush(stdout);
        }
      }
    }

    #pragma omp critical
    {
      if(HitMemHWM > HitHWM)
	HitHWM = HitMemHWM;
      if(MatchMemHWM > MemHWM){
	MemHWM = MatchMemHWM;
	HWMid = HWMmapid;
      }
      if(MatchBufHWM > BufHWM)
	BufHWM = MatchBufHWM;
      if(BestCntHWM > BestHWM)
	BestHWM = BestCntHWM;
      if(OutputHWM > OutHWM)
	OutHWM = OutputHWM;

      TotalMatchCnt += totalMcnt;
      if(STAT){
	gHTindcnt += HTindcnt;
	gHTqindcnt += HTqindcnt;
	gHTfindcnt += HTfindcnt;
	gHTcolcnt += HTcolcnt;
	gHTqcolcnt += HTqcolcnt;
	gHTmatchcnt += HTmatchcnt;
	gHTentrycnt += HTentrycnt;
	gHTcnt += HTcnt;
	gHTcol += HTcol;
	gHTdup += HTdup;
	gMTcnt += MTcnt;
	gMTcnt2 += MTcnt2;
	gMTcnt3 += MTcnt3;
	gMTcol += MTcol;
	gMTqcol += MTqcol;
	gMTicnt += MTicnt;
      }
      if(VERB && (VERB>=2 /* || STAT */)){/* display amount of block memory used */
        int tid = 0;
#ifdef _OPENMP
	tid = omp_get_thread_num();
#endif
	printf("tid=%d:CMatchTable() destructor:",tid);
	Hit_memPrint(1);
	printf("tid=%d:CMatchTable() destructor:",tid);
	Match_memPrint(1);
      }
    }
    //      mcheck_check_all();
    
    if(!matchtable && fbuf) delete [] fbuf;

    if(sizekeyV) delete [] sizekeyV;
    if(rand1V) delete [] rand1V;
    //    if(HashIndexV) delete [] HashIndexV;
    if(hashV) delete [] hashV;

    delete [] matches;

    free(HitMem);
    delete [] HitIndex;
    delete [] HitNext;

    //      mcheck_check_all();

    // WAS    free(MatchIndex);
    size_t memsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * (1UL << MHASH_BITS) + PAGE-1) & ~(PAGE-1);
    if(munmap(MatchIndex, memsiz)){
      int eno = errno;
      char *err = strerror(eno);
	
      #pragma omp critical
      {
        printf("munmap(%p,%lu) failed: errno=%d:%s\n",MatchOverflow, memsiz, eno, err);
        dumpmemmap();
        fflush(stdout);exit(1);
      }
    }
    MatchIndex = NULL;

    //      mcheck_check_all();

    if(MatchOverflow){
      // WAS    free(MatchOverflow);
      size_t memsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * MatchOverflowMax + PAGE-1) & ~(PAGE-1);
      if(munmap(MatchOverflow, memsiz)){
	int eno = errno;
	char *err = strerror(eno);
	
        #pragma omp critical
	{
          printf("munmap(%p,%lu) failed: errno=%d:%s\n",MatchOverflow, memsiz, eno, err);
	  dumpmemmap();
	  fflush(stdout);exit(1);
	}
      }
      MatchOverflow = NULL;
    }

    //      mcheck_check_all();

    if(MatchBufMax > 0){
      if(DEBUG) assert(MatchBuf != NULL);

      // WAS delete [] MatchBuf;
      size_t memsiz = (sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>) * MatchBufMax + PAGE-1) & ~(PAGE-1);
      if(munmap(MatchBuf, memsiz)){
	int eno = errno;
	char *err = strerror(eno);

        #pragma omp critical
        {
	  printf("munmap(%p,%lu) failed: errno=%d:%s\n",MatchBuf, memsiz, eno, err);
	  dumpmemmap();
	  fflush(stdout);exit(1);
	}
      }
    } else if(DEBUG) assert(MatchBuf == NULL);

    if(splitcnt > 1){
      BestMatch += roffset;
      if(HashMultiMatch){
        BestMatchCnt += roffset;
        BestMatchMax += roffset;
      }
    }

    delete [] BestMatchMem;
    delete [] BestMatch;
    delete [] BestMapid;
    delete [] BestMatchCnt;
    delete [] BestMatchMax;

    //      mcheck_check_all();
  };
  

  void hash_query(Cmap *pmap, int tid, int *Xflag, FLOAT *Xrev, FLOAT *Xscale, CMatchTable *matchtable2); /* compute all matches for one map (both orientations) in hashtable and append them to the output buffer in sorted order */

  void Hit_memPrint(int verb);/* Check mount of CHitEntry memory on freelist/used/leaked (for debugging) : exit if memory has leaked */
  void Match_memPrint(int verb);/* Check amount of CMatchEntry memory on freelist/used/leaked (for debugging) : exit if memory has leaked */
};

static inline int key1Inc(long long *p1, long long *p2)
{
  long long key1 = *p1;
  long long key2 = *p2;

  return (key1 > key2) ? 1 : (key1 < key2) ? -1 : 0;
}

#if 0
/* qsort() intcmp function to sort CHashMatch array in increasing order of id2,orientation,(-hashscore),offset */
static inline int ChashIdInc(CHashMatch *p1, CHashMatch *p2)
{
  return CHashMatchIncId1(p1,p2);
}
#endif

#define SIZES (1 << SIZE_BITS) /* number of quantized interval sizes */

#define MAPID2 -1 // 7 // -1 // 594 // 87
#define FLIP2 -1 // 1 // -1 // 1 // 1 // 0 
#define MAXSITE2 -1 // 5400 // 4800 // 600

extern double mtime();/* microsecond timer function : see refine.cpp */
extern double wtime();

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/hash.cpp 11578 2020-08-31 22:42:04Z tanantharaman $");

static unsigned char OFFSET_WT[1+OFFSET_SPREAD] = {5,4,0} /* WAS {1,1,0} */;/* weight window for different offset mismatches. WT[0] >= 1 && WT[1..OFFSET_SPREAD] >= 0 and WT[0]+2*sum(WT[1..OFFSET_SPREAD]) <= 255 */

static int hcolor = 0;/* which color to use : if colors >= 2 and hcolor = -1, use all colors */

static double varSF, varSD, varSR;
static PFLOAT fvarSF, fvarSD, fvarSR;
static PFLOAT SDnormMax;

static PFLOAT SFm,SDm,SRm,SD1,VAR1,SRInv;

static int Gpairwise = 0;/* global copy of arg "pairwise" */

static bool *MapSave = NULL;/* Used to mark maps that need to be saved */

FILE *MatrixFP = NULL;/* pointer to <prefix>.matrix : only used if HashMatrix > 0 */

int splitcnt = 1;
int roffset = 0;/* starting index in rmap[roffset .. roffset + numrmaps - 1] */
static int orig_numrmaps = 0;
static long long cumTotalMatchCnt = 0;

int hnumthreads = 1; /* number of threads to use in parallel sections of hash.cpp, hash.h */
int numbufs = 0;/* number of output buffers : typically the same as numthreads */
static char **bufnames = NULL;/* output buffers bufnames[0..numbufs-1] */
short *besthashscore = NULL;/* If -hashbest, besthashscore[i=0..numrmaps-1] is the best hashscore for any match with query (probe) mapid i (If pairalign, also best hashscore for any ref (insert) mapid i) */


double MaxFragLen = 0.0; 
unsigned int *IntLen = NULL;/* Table to convert any size X between resKB < X <= MaxFragLen to a small integer (quantized representation of length) */
int MaxIntLen = 0;

unsigned int MaxHit = 0;/* limit of HashTable entries with the same key == MAXHIT_COVERAGE/(1+HashWin) */

/* global statistics */
unsigned int HitHWM=0; int MemHWM = 0, HWMid = -1;/* High Water Mark for HitMax & MatchOverflowMax used by any Match Table for any single query (HWMid) */
long long BufHWM = 0;/* High Water Mark for MatchBufMax used by any Match Table for any single query */
int BestHWM = 0;/* High Water Mark for BestCnt used by any Match Table for any single query */
size_t OutHWM = 0; /* High Water Mark for matchmax used by any Match Table for any single query */
long long TotalMatchCnt;/* progress counter */
 long long gHTindcnt, gHTqindcnt, gHTfindcnt, gHTcolcnt, gHTqcolcnt, gHTmatchcnt, gHTentrycnt;/* HashTable performance statistics */
long long gMTcnt, gMTcnt2, gMTcnt3, gMTcol, gMTqcol,gMTicnt;/* MatchTable performance statistics */
long long gHTcnt, gHTcol, gHTdup;/* HitTable performance statistics */
FLOAT **YY = NULL;/* YY[mapid1 = 0..numrmaps-1] == rmap[mapid1]->site[hcolor] */
FLOAT *Ylen = NULL;/* Ylen[mapid1 = 0..numrmaps-1] == YY[mapid1][rmap[mapid1]->numsite[hcolor]+1] */
 /*OFFSET_T*/ void *maxOffset1; /* Used with flip2 == 0 : maxOffset1[mapid1]  == floor(YY[mapid1][max(0,N1 - HashWin)] / OffsetKB + 0.5) */
 /*OFFSET_T*/ void *minOffset1R;/* Used with flip2 == 1 : minOffset1R[mapid1] == floor((Ylen[mapid1] - YY[mapid1][N1]) / OffsetKB + 0.5) */

//int ShiftOffset1 = 15;

int DisJoint = 0;/* 0 if rmap[],qmap[] are identical OR have shared maps (only combinations with rmap[]->id > qmap[]->id are aligned),
		    1 if rmap[],qmap[] are Disjoint sets of maps with rmap[] having larger ids than qmap[]
		    NOTE: If HashMatrix==1, DisJoint is set to 1 */

int Nmax = 0;/* largest number of sites for any Probe maps (qmap[i]->numsite[hcolor]) */

unsigned int SizeLen = 0;/* 3^HashWin */
static long long *SizeDelta = NULL;/* SizeDelta[3^HashWin] : all possible hash key offsets to account for sizing error of +- 1 in every size of the window */

static double OffsetKBinv;/* 1.0/OffsetKB */
int MinOffset=0,MaxOffset=0;/* offset values are in the range MinOffset ... MaxOffset */
int MaxSite = 0;/* maximum number of sites in any rmap[] */
long long TotWin = 0;/* total number of windows over all rmap[]'s */

size_t maxMatchesPerRefid = 0, maxMapsPerRefid = 0;
int Gnumqmaps = 0;
size_t *MatchesPerRefid = NULL, *MapsPerRefid = NULL;// MatchesPerRefid[0..numqmaps-1], MapsPerRefid[0..numqmaps-1] 

static Cid2mapid *id2mapid = NULL;

static long long maxMatchMem = 0;
//static Resource *MatchMemReserve = NULL; // moved to CHashTable

static int IntDec(int *p1, int *p2)
{
  return *p2 - *p1;
}

/* sort by increasing order of number of sites */
static int CmapMapidInc(Cmap **p1, Cmap **p2)
{
  return p1[0]->mapid - p2[0]->mapid;
}

/** qsort() intcmp function to sort CHitEntry array in increasing order of mapid1,offset1,scaleID[0] */
template<class OFFSET1_T>
static int CHitMapidOffsetInc(CHitEntry<OFFSET1_T> *p1, CHitEntry<OFFSET1_T> *p2)
{
  return (p1->mapid1 > p2->mapid1) ? 1 : (p1->mapid1 < p2->mapid1) ? -1 : (p1->offset1 - p2->offset1);
}

/** qsort() intcmp function to sort CHashEntry<OFFSET1_T> array in increasing order of mapid1 */
template<class OFFSET1_T>
static int CHashEntryMapid1Inc(CHashEntry<OFFSET1_T> *p1, CHashEntry<OFFSET1_T> *p2)
{
  return (p1->mapid1 > p2->mapid1) ? 1 : (p1->mapid1 < p2->mapid1) ? -1 : 0;
}

/** qsort() intcmp function to sort CHashCollideR array in increasing order of key */
static int CHashCollideRKeyInc(CHashCollideR *p1, CHashCollideR *p2)
{
  return (p1->key > p2->key) ? 1 : (p1->key < p2->key) ? -1 : 0;
}

/** qsort() intcmp function to sort CHashMatch array in increasing order of id2,orientation,(-hashscore),offset (assume all id1 values are the same) */
static inline int CHashMatchIncId2(CHashMatch *p1, CHashMatch *p2)
{
  return (p1->id2 > p2->id2) ? 1 : (p1->id2 < p2->id2) ? -1 : 
    (p1->orientation > p2->orientation) ? 1 : (p1->orientation < p2->orientation) ? -1 :
    (p1->hashscore < p2->hashscore) ? 1 : (p1->hashscore > p2->hashscore) ? -1 :
    (p1->offset - p2->offset);
}
/** qsort() intcmp function to sort CHashMatch array in decreasing order of hashscore (assume all id1 values are the same) */
static inline int CHashMatchDecHscore(CHashMatch *p1, CHashMatch *p2)
{
  return (p1->hashscore < p2->hashscore) ? 1 : (p1->hashscore > p2->hashscore) ? -1 : 0;
}


/** qsort() intcmp function to sort Cmap* array in increasing order of id */
static int CmapIdInc(Cmap **p1, Cmap **p2)
{
  /* Cannot use p1[0]->id - p2[0]->id since that may be out of range for integer */
  return (p1[0]->id > p2[0]->id) ? 1 : (p1[0]->id < p2[0]->id) ? -1 : 0;
}

/** qsort() intcmp function to sort CBestMatch array in increasing order of offset, scaleID then (if RANGE>=1) MinOffset1,MaxOffset1 */
template <class OFFSET_T, int RANGE, int HASH_SCALE>
static int CBestMatchOffsetScaleInc(CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *p1, CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *p2)
{
  return (p1->offset != p2->offset) ? p1->offset - p2->offset :
    (HASH_SCALE && p1->scaleID[0] > p2->scaleID[0]) ? 1 : (HASH_SCALE && p1->scaleID[0] < p2->scaleID[0]) ? -1 : 
    (RANGE>=1 && p1->MinOffset1[0] > p2->MinOffset1[0]) ? 1 : (RANGE>=1 && p1->MinOffset1[0] < p2->MinOffset1[0]) ? -1 :
    (RANGE>=1 && p1->MaxOffset1[0] > p2->MaxOffset1[0]) ? 1 : (RANGE>=1 && p1->MaxOffset1[0] < p2->MaxOffset1[0]) ? -1 : 0;
}

/** qsort() intcmp function to sort CBestMatch array in increasing order of offset, then decreasing order of score, then (if SCALE_SMALL) increasing order of scaleID, 
    then (if RANGE>=1) increasing order of MinOffset1,MaxOffset1 */
template <class OFFSET_T, int RANGE, int HASH_SCALE>
static int CBestMatchOffsetInc(CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *p1, CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *p2)
{
  return (p1->offset != p2->offset) ? p1->offset - p2->offset :
    (p1->score != p2->score) ? p2->score - p1->score : 
    (HASH_SCALE && SCALE_SMALL && p1->scaleID[0] > p2->scaleID[0]) ? 1 : (HASH_SCALE && SCALE_SMALL && p1->scaleID[0] < p2->scaleID[0]) ? -1 : 
    (RANGE>=1 && p1->MinOffset1[0] > p2->MinOffset1[0]) ? 1 : (RANGE>=1 && p1->MinOffset1[0] < p2->MinOffset1[0]) ? -1 :
    (RANGE>=1 && p1->MaxOffset1[0] > p2->MaxOffset1[0]) ? 1 : (RANGE>=1 && p1->MaxOffset1[0] < p2->MaxOffset1[0]) ? -1 : 0;
}

/** qsort() intcmp function to sort CBestMatch array in decreasing order of hashscore */
template <class OFFSET_T, int RANGE, int HASH_SCALE>
static int CBestMatchScoreDec(CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *p1, CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *p2)
{
  return p2->score - p1->score;
}

/** Increase size of CHashMatch array (if needed) */
static inline void maxmatchalloc(size_t cnt, size_t newcnt, size_t &maxcnt, CHashMatch* &matches)
{
  if(newcnt <= maxcnt)
    return;
  if(newcnt < maxcnt*2)
    newcnt = maxcnt*2;
  if(newcnt < 64*1024)
    newcnt = 64*1024;
  try {
    if(maxcnt <= 0){
      matches = new CHashMatch[newcnt];
    } else {
      CHashMatch *origMatches = matches;
      matches = new CHashMatch[newcnt];
      memcpy(matches, origMatches, cnt * sizeof(CHashMatch));
      delete [] origMatches;
    }
  } catch (exception& e){
    cout << e.what() << endl;
    printf("maxmatchalloc(cnt=%lu,newcnt=%lu,maxcnt=%lu):failed to allocate newcnt CHashMatch entries (%lu bytes)\n",
	   cnt, newcnt, maxcnt, newcnt * sizeof(CHashMatch));
    fflush(stdout);
    assert(0);
  }
  maxcnt = newcnt;
}

/** allocate memory from HitMem[HitFree..HitMax-1] : returns index into HitMem[] */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline int CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::CHitEntryAlloc()
{
  if(DEBUG>=3) Hit_memPrint(0);

  if(HitFree >= HitMax){/* increase size of HitMem[] */
    unsigned int newHitMax = (HitMax * 3)/2;
    assert(newHitMax > HitFree);/* in case 32 bit HitMax wrapped around */
    CHitEntry<OFFSET1_T> *newHitMem = 0; POSIX_MEMALIGN(newHitMem,64,sizeof(CHitEntry<OFFSET1_T>)*newHitMax);
    memcpy(newHitMem,HitMem,HitMax*sizeof(CHitEntry<OFFSET1_T>));
    free(HitMem);
    HitMem = newHitMem;
    HitMax = newHitMax;
  }
  return HitFree++;
}

/** print amount of CHitEntry memory on Hit_freelist and used in Hit Table (for debugging) */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::Hit_memPrint(int verb)
{
  int usecnt = 0;
  for(int hindex = HitFirst;hindex >= 0; hindex = HitNext[hindex]){  /* Traverse the Hit Table */
    if(HitIndex[hindex]==0){
      printf("Hit_memPrint:HitFirst=%d,hindex=%d,HitIndex[hindex]=%d\n",HitFirst,hindex,HitIndex[hindex]);
      fflush(stdout);
      assert(HitIndex[hindex] != 0);
    }
    for(int qindex = HitIndex[hindex];qindex; qindex = HitMem[qindex].next)
      usecnt++;
  }
  if(verb){
    printf("tid=%d,HitMax=%d structs,Free=%d structs, used=%d structs\n",tid,HitMax,HitMax-HitFree, usecnt);
    fflush(stdout);
  }
}

/** print amount of CMatchEntry memory on Match_freelist and used in Match Table (for debugging) */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::Match_memPrint(int verb)
{
  if(verb){
    if(DEBUG && HWMmapid >= 0){
      printf("tid=%d,MatchOverflow=%d/%d,HWM=%d+%lld(%llu bytes per thread) HWMmapid=%d, MatchIndex=%u(%lu bytes per thread), BestHWM=%d,OutputHWM=%lu(%lu bytes per thread)\n",
	     tid,MatchOverflowNext-MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)-1, MatchOverflowMax-MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)-1, MatchMemHWM,MatchBufHWM,
	     MatchMemHWM*sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) + MatchBufHWM*sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), HWMmapid, (1U << MHASH_BITS), (1U << MHASH_BITS)*sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>),
	     BestCntHWM, OutputHWM, BestCntHWM*sizeof(int) + numrmaps*sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>) + OutputHWM*sizeof(CHashMatch));
    } else
      printf("tid=%d,MatchOverflow=%d/%d, MatchIndex=%u(%lu bytes per thread)\n",
	     tid,MatchOverflowNext-MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)-1, MatchOverflowMax-MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)-1, (1U << MHASH_BITS), (1U << MHASH_BITS)*sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>));
    fflush(stdout);
  }
}

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline int CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::MatchOverflowAlloc(int tid)
{
  if(MatchOverflowFree > 0){
    unsigned int ret = MatchOverflowFree;
    if(DEBUG>=2) assert(ret > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
    MatchOverflowFree = MatchOverflow[ret].next;
    if(DEBUG>=2) assert(MatchOverflowFree < 0 || MatchOverflowFree > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
    return ret;
  }

  if(MatchOverflowNext >= MatchOverflowMax){
    unsigned int origMatchOverflowMax = MatchOverflowMax;
    MatchOverflowMax = min(MASK_LL(31), max(MatchOverflowNext, MatchOverflowMax) * 3LL / 2LL);
    if(MatchOverflowNext >= MatchOverflowMax){
      printf("MatchOverflowAlloc:: MatchOverflowMax= %u cannot be increased using unsigned 32 bit int: Memory table overflow\n",MatchOverflowMax);
      fflush(stdout);exit(1);
    }

    if(hashtable->MatchMemReserve){
      if(VERB>=2){
        #pragma omp critical
	{
	  printf("tid=%d: MatchOverflowAlloc:Max=%u -> %u (%lu bytes), MatchOverflowNext=%u\n",
		 tid,origMatchOverflowMax,MatchOverflowMax,(size_t)MatchOverflowMax * sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>), MatchOverflowNext);
	  fflush(stdout);
	}
      }

      size_t memsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * origMatchOverflowMax + PAGE-1) & ~(PAGE-1);

      hashtable->MatchMemReserve->allocate((size_t)(MatchOverflowMax-origMatchOverflowMax) * sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>), 0, memsiz, this);
    }
    assert(MatchOverflowMax > MatchOverflowNext);/* in case 32 bit value wrapped around */

    // WAS     CMatchLine *NewOverflow = 0; POSIX_MEMALIGN(NewOverflow,64,sizeof(CMatchLine)*MatchOverflowMax);
    size_t memsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * MatchOverflowMax + PAGE-1) & ~(PAGE-1);
    char *newmem = (char *) mmap(NULL, memsiz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if(newmem == MAP_FAILED){
      int eno = errno;
      char *err = strerror(eno);
	    
      #pragma omp critical
      {
        printf("CMatchTable: mmap of %lu bytes failed:errno=%d:%s\n", memsiz,eno,err);
	dumpmemmap();
	fflush(stdout);exit(1);
      }
    }
    CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *NewOverflow = (CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *)newmem;

    memcpy(NewOverflow,MatchOverflow,MatchOverflowNext * sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>));

    // WAS  free(MatchOverflow);
    memsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * origMatchOverflowMax + PAGE-1) & ~(PAGE-1);
    if(munmap(MatchOverflow, memsiz)){
      int eno = errno;
      char *err = strerror(eno);

      #pragma omp critical
      {
        printf("munmap(%p,%lu) failed: errno=%d:%s\n",MatchOverflow, memsiz, eno, err);
	dumpmemmap();
	fflush(stdout);exit(1);
      }
    }

    MatchOverflow = NewOverflow;
  }

  return MatchOverflowNext++;
}

/** Insert hit into CMatchTable::HitIndex, removing duplicates (same mapid1,offset1). This table is cleared 
    after each query (mapid2,flip2, scaleID, site2) so there is only data from a single query map (mapid2), orientation (flip2), scaling (scaleID) and site window (site2) in the table.\n

    Hit table is indexed by the ref map:
    > hindex = Maphash1[mapid1] ^ Sitehash[offset1];\n
    > CMatchTable::HitIndex[hindex];
    Link between reference map sites by CMatchTable::HitNext[].\n
    Link between query map sites by ->next.
    */

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::HitInsert(int mapid1, OFFSET1_T offset1, OFFSET_T offset2, int site1, int TraceFlag)
{
  if(DEBUG>=2) assert(roffset <= mapid1 && mapid1 < roffset+numrmaps);
  if(DEBUG>=2) assert(MinOffset <= offset1 && offset1 <= MaxOffset);
  register unsigned int hindex = Maphash1[mapid1] ^ Sitehash[offset1];
  if(DEBUG>=2) assert(hindex < (1U << HHASH_BITS));
  register int index = HitIndex[hindex];
  if(STAT) HTcnt++;
  if(!index){/* start new list */
    HitIndex[hindex] = index = CHitEntryAlloc();
    if(DEBUG>=2) assert(index > 0);
    register CHitEntry<OFFSET1_T> *p = &HitMem[index];
    p->init(mapid1, offset1,/*  offset2,*/ site1, 0);
    HitNext[hindex] = HitFirst;
    HitFirst = hindex;
    if(TRACE && TRACE2 >= 2 && TraceFlag)
      printf("\tHitInsert of mapid1=%d(%lld),offset1=%d,offset2=%d at index=%d hindex=%u\n",mapid1,rmap[mapid1]->id,offset1,offset2,index,hindex);
    if(DEBUG>=3) Hit_memPrint(0);
    return;
  }
  
  /* see if matching entry is already in list p */
  while(1){
    register CHitEntry<OFFSET1_T> *p = &HitMem[index];
    if(p->mapid1==mapid1 && p->offset1 == offset1){
      if(TRACE && TRACE2 >= 2 && TraceFlag)
	printf("\tHitInsert of mapid1=%d(%lld),offset1=%d,%d: duplicate at index=%d hindex=%u\n",mapid1,rmap[mapid1]->id,offset1,offset2,index,hindex);
      if(STAT) HTdup++;
      return;
    }
    if(STAT) HTcol++;
    if(!(index = p->next))
      break;
  }

  /* append new entry to list HitIndex[hindex] */
  index = CHitEntryAlloc();
  register CHitEntry<OFFSET1_T> *q = &HitMem[index];
  q->init(mapid1,offset1,/*offset2,*/site1,HitIndex[hindex]);
  HitIndex[hindex] = index;
  if(TRACE && TRACE2 >= 2 && TraceFlag)
    printf("\tHitInsert of mapid1=%d(%lld),offset1=%d,offset2=%d at index=%d hindex=%u\n",mapid1,rmap[mapid1]->id,offset1,offset2,index,hindex);
  if(DEBUG>=3) Hit_memPrint(0);
}

/* Copy all (non-duplicate) entries in Hit Table to Match Table, and reset Hit Table */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::HitToMatch(int mapid2, int flip2, FLOAT *X, int M, int site2, OFFSET_T offset2, int scaleID, int tid)
{
  if(VERB>=2)
    printf("HitToMatch:mapid2=%d,flip2=%d,site2=%d,offset2=%d,scaleID=%d:\n",mapid2,flip2,site2,(int)offset2,scaleID);

  if(DEBUG>=2) assert(M == qmap[mapid2]->numsite[hcolor]);
  
  if(VERB && HitMax > HitMemHWM)
    HitMemHWM = HitMax;

  if(HashMatrix){/* just append sorted HitTable entries to FILE *MatrixFP */
    CHitEntry<OFFSET1_T> *HitBuf = new CHitEntry<OFFSET1_T>[HitFree];
    unsigned int HitCnt = 0;

    int nindex = HitFirst;
    for(int hindex = nindex; hindex >= 0; hindex = nindex){  /* Traverse the Hit Table */
      if(DEBUG>=2) assert(HitIndex[hindex] != 0);
      CHitEntry<OFFSET1_T> *p = &HitMem[HitIndex[hindex]],*q,*r;

      /* loop over linked list starting at p and copy all entries to MatchTable */
      for(q = p;; q = r){
	int rindex = q->next;
	r = &HitMem[rindex];
	
	/* append *q to buffer HitBuf[] */
	HitBuf[HitCnt++] = *q;

	if(!rindex)
	  break;
      }

      nindex = HitNext[hindex];

      /* reset HitIndex[hindex] */
      HitIndex[hindex] = 0;
    }
    if(DEBUG) assert(HitCnt <= HitFree);

    /* sort HitBuf[0..HitCnt-1] by mapid1,offset1 */
    qsort(HitBuf,HitCnt,sizeof(CHitEntry<OFFSET1_T>),(intcmp*)CHitMapidOffsetInc<OFFSET1_T>);

    /* append HitBuf[0..HitCnt-1] to FILE *MatrixFP */
    if(DEBUG) assert(MatrixFP != NULL);
    #pragma omp critical 
    {
      for(unsigned int i = 0; i < HitCnt; i++){
	CHitEntry<OFFSET1_T> *p = &HitBuf[i];
	fprintf(MatrixFP,"%lld %d %lld %d %d\n", qmap[mapid2]->id, (int)floor(offset2 * OffsetKB + 0.5), rmap[p->mapid1]->id, (int)floor(p->offset1 * OffsetKB + 0.5), flip2);
      }
      fflush(MatrixFP);
    }

  } else {/* HashMatrix == 0 : normal case optimized for speed */
    int nindex = HitFirst;
    _mm_prefetch(CAST &HitIndex[nindex], _MM_HINT_T0);

    for(int hindex = nindex; hindex >= 0; hindex = nindex){  /* Traverse the Hit Table */
      if(DEBUG>=2 && HitIndex[hindex]==0){
	printf("hindex=%d:HitIndex[hindex]=%d\n",hindex,HitIndex[hindex]);
	fflush(stdout);
	assert(HitIndex[hindex] != 0);
      }
      CHitEntry<OFFSET1_T> *p = &HitMem[HitIndex[hindex]],*q,*r;
      _mm_prefetch(CAST p, _MM_HINT_T0);
      _mm_prefetch(CAST &HitNext[hindex], _MM_HINT_T0);

      /* loop over linked list starting at p and copy all entries to MatchTable */

      if(RANGE && ShiftOffset1 <= 30){
	for(q = p;; q = r){
	  int rindex = q->next;
	  r = &HitMem[rindex];
	  _mm_prefetch(CAST r, _MM_HINT_T0);

	  /* copy *q to MatchTable */
	  //	  if(DEBUG) assert(offset2 == q->offset2);
	  OFFSET_T offset = offset2 - q->offset1;

	  if(DEBUG>=2) assert(offset >= MinOffset && offset <= MaxOffset);

	  unsigned char score = MatchInsert(q->mapid1,offset,1,mapid2,flip2, offset2, scaleID, 1 + (q->offset1 >> ShiftOffset1),q->offset1);

	  if(VERB>=2 || (TRACE && rmap[q->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	    printf("HitToMatch:q=%p,hindex=%u:mapid1=%d(id=%lld,len=%0.3f),mapid2=%d(id=%lld,len=%0.3f),flip2=%d,site2=%d,offset2=%u,scaleID=%u:offset=%d,offset1=%d(%u),score=%d:",
		   q,hindex,q->mapid1,rmap[q->mapid1]->id,Ylen[q->mapid1],mapid2,qmap[mapid2]->id,X[M+1],flip2,site2,offset2,scaleID,offset,q->offset1,1 + (q->offset1 >> ShiftOffset1), score);
	    printf(":X[%d..%d]=[%0.3f..%0.3f]",site2,min(M,site2+HashWin+1),X[site2],X[min(M,site2+HashWin+1)]);
	    printf("\n");
	  }

	  if(!rindex)
	    break;
	}
      } else {/* simpler code without ShiftOffset1 */
	for(q = p;; q = r){
	  int rindex = q->next;
	  r = &HitMem[rindex];
	  _mm_prefetch(CAST r, _MM_HINT_T0);

	  /* copy *q to MatchTable */
	  //	  if(DEBUG) assert(offset2 == q->offset2);
	  OFFSET_T offset = offset2 - q->offset1;

	  if(DEBUG>=2) assert(offset >= MinOffset && offset <= MaxOffset);

	  unsigned char score = MatchInsert(q->mapid1,offset,1,mapid2,flip2, offset2, scaleID, 0, q->offset1);

	  if(VERB>=2 || (TRACE && rmap[q->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	    printf("HitToMatch:q=%p,hindex=%u:mapid1=%d(id=%lld,len=%0.3f),mapid2=%d(id=%lld,len=%0.3f),flip2=%d,site2=%d,offset2=%u,scaleID=%u,offset=%d,offset1=%d,score=%d:",
		   q,hindex,q->mapid1,rmap[q->mapid1]->id,Ylen[q->mapid1],mapid2,qmap[mapid2]->id,X[M+1],flip2,site2,offset2,scaleID,offset,q->offset1,score);
	    printf(":X[%d..%d]=[%0.3f..%0.3f]",site2,min(M,site2+HashWin+1),X[site2],X[min(M,site2+HashWin+1)]);
	    printf("\n");
	  }

	  if(!rindex)
	    break;
	}
      }

      nindex = HitNext[hindex];
      _mm_prefetch(CAST &HitIndex[nindex], _MM_HINT_T0);

      /* reset HitIndex[hindex] */
      HitIndex[hindex] = 0;
    }
  }/* if(HashMatrix==0 */

  /* reset HitTable */
  HitFree = 1;
  HitFirst = -1;

  if(DEBUG>=3){
    /* verify that HitTable is empty */
    for(int hindex = (1 << HHASH_BITS);--hindex >= 0;)
      assert(HitIndex[hindex]==0);
    Hit_memPrint(0);
  }

  /* free up memory of HitTable (if needed) */
  if(HitMax > (1U << (HHASH_BITS-2))){
    free(HitMem);
    HitMax = (1U << (HHASH_BITS-2));
    POSIX_MEMALIGN(HitMem,64,sizeof(CHitEntry<OFFSET1_T>)*HitMax);
  }

  if(MTHashGC > 0 && !Gpairwise && !(site2 % MTHashGC))
    MatchFind(mapid2,flip2,site2,tid);
}

/* 
Locate the best remaining matches in MatchTable that are above hashscore threshold and append them to output buffer CMatchTable::matches using CMatchTable::BestInsert().
This completes hashtable query, so reset MatchIndex[]/MatchOverflow[], MatchBuf[] and BestMatch[]

If HashGC && (minSite2 > 0 OR maxSite2 > 0): Only offset2 values from offset2(flip2,minSite2) - 2*max_spread .. offset2(flip2,maxSite2) are in MatchTable (If maxSite2 == 0, all remaining matches are in 
MatchTable, so replace offset2(flip2,maxSite2) by Infinity). Apply the following modified steps:

    1a. If flip2 == 0 : For each mapid1 only matches with (offset2(0,minSite2) - maxOffset1[mapid1] - max_spread <= offset <= offset2(0,maxSite2) - maxOffset1[mapid1] - max_spread - 1) go to MatchBuf[]
    1b. If flip2 == 1 : For each mapid1 only matches with (offset2(1,MaxSite2) - minOffset1R[mapid1] + max_spread + 1 <= offset <= offset2(1,minSite2) - minOffset1R[mapid1] + max_spread) go to MatchBuf[]
    2. Apply OffsetSpread to MatchBuf[] (this may use matches in MatchIndex/MatchOverflow with offsets slightly outside the range in step 1)
    3. If HashMultiMatch && minSite2 > 0 : (Re)allocate BestMatch[] and add contents of MatchBuf[] with offset in target range (see 1. above) to BestMatch.
    4. Reset MatchBuf[]
    5. If HashMultiMatch && maxSite2 == 0 : Keep only locally best scoring matches (maxSite==0 means all matches are now in BestMatch).
    6a. If maxSite2 > 0 && flip2==0 : Remove all entries from MatchIndex/MatchOverflow with offset < offset2(0,maxSite2) - maxOffset1[mapid1] - 2 * max_spread
    6b. If maxSite2 > 0 && flip2==1 : Remove all entries from MatchIndex/MatchOverflow with offset > offset2(1,MaxSite2) - minOffset1R[mapid1] + 2 * max_spread
    7. If maxSite2 <= 0 : Free up MatchOverflow[] and MatchIndex[]
    8. If maxSite2 == 0 : output BestMatch[] and reset BestMatch[]
    9. If HashGC && maxSite2 > 0 : Update minSite2 = maxSite2

  Where:
       offset2(0,site2) == floor(X[site2] * minScaleFactor / OffsetKB + 0.5)
       offset2(1,site2) == floor((X[N2 + 1 - site2 - HashWin] * maxScaleFactor / OffsetKB + 0.5)
       maxScaleFactor == 1 + HASH_SCALE * ScaleDelta * ScaleDeltaSize (NOT YET IMPLEMENTED)
       minScaleFactor == 1 / maxScaleFactor (NOT YET IMPLEMENTED)
       X[] == qmap[mapid2]->site[hcolor][]
       N2 = qmap[mapid2]->numsite[hcolor]
       max_spread is the largest value delta <= OffsetSpread with OFFSET_WT[delta] > 0
 */

const int PREFETCH_D = 5;/* prefetch distance (iterations) : the number of outstanding prefetch instruction is (2 + RANGE + HASH_SCALE) * PREFETCH_D */

const int PREFETCH_D2 = 2 * PREFETCH_D;

#define PREFETCH_BM 4 // number of outstanding prefetch instructions is 3 * PREFETCH_BM 

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::MatchFind(int mapid2, int flip2, int maxSite2, int tid)
{
  register OFFSET_T *Offset1flip2 = flip2 ? (OFFSET_T*)minOffset1R : (OFFSET_T*)maxOffset1;

  if(DEBUG>=3 && !(minSite2 > 0)){
    assert(BestCnt==0);
    if(HashMultiMatch){
      for(int mapid1 = numrmaps; --mapid1 >= 0; )
	assert(BestMatchCnt[mapid1] == 0);
    } else {
      for(int mapid1 = numrmaps; --mapid1 >= 0; )
	assert(BestMatch[mapid1]->mapid1 < 0);
    }
  }
  
  if(VERB>=2 || (TRACE && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2)){
    printf("MatchFind:mapid2=%d(id=%lld),flip2=%d,minSite2=%d,maxSite2=%d,MGHashGC=%d,N2=%d:tid=%d:MatchCnt=%lld\n",
	   mapid2,qmap[mapid2]->id,flip2,minSite2,maxSite2,MTHashGC,qmap[mapid2]->numsite[hcolor],tid,MatchCnt);
    fflush(stdout);
  }

  if(VERB && MatchOverflowMax > MatchMemHWM){
    if(VERB>=2){
      printf("MatchFind:mapid2=%d(id=%lld),flip2=%d:MatchCnt=%lld:MatchMemHWM=%u -> %u\n",mapid2,qmap[mapid2]->id,flip2,MatchCnt, MatchMemHWM, MatchOverflowMax);
      fflush(stdout);
    }
    MatchMemHWM = MatchOverflowMax;
  }

  if(DEBUG>=3){/* check if MatchCnt is valid */
    long long TrueCnt1 = 0, TrueCnt2 = 0, TrueCnt3 = 0;
    int IndexCnt = 0, FreeCnt = 0;    
    register int nindex = MatchFirst;
    register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *q = &MatchIndex[nindex], *p;
    for(register int mindex = nindex; mindex >= 0; mindex = nindex){
      p = q;
      q = &MatchIndex[nindex = p->next];
      
      unsigned int cnt = p->cnt;
      if(DEBUG>=2) assert(cnt > 0);
      if(cnt > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){
	for(int i = cnt; i > 0; i = MatchOverflow[i].next){
	  if(DEBUG>=2 && !(MatchOverflow[i].cnt > 0)){
	    #pragma omp critical
	    {
	      printf("tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d:mindex=%d,nindex=%d,p->cnt=%d:i=%d,MatchOverflow[i].next=%d,cnt=%d\n",
		     tid,mapid2,flip2,minSite2,maxSite2,mindex,nindex,p->cnt,i,MatchOverflow[i].next,MatchOverflow[i].cnt);
	      fflush(stdout);
	      
	      assert(MatchOverflow[i].cnt > 0);
	    }
	  }
	}
      }

      IndexCnt++;
      TrueCnt1 += min(MATCHLINE(OFFSET_T,RANGE,HASH_SCALE),p->cnt);
    }
    p = MatchOverflow;

    for(int i = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)+1; i < MatchOverflowNext; i++){
      if(DEBUG>=2 && !HashGC && !(0 < p[i].cnt && p[i].cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE))){
        #pragma omp critical
	{
	  printf("tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d:MatchCnt=%lld,IndexCnt=%d,Overflow Lines=%d,MatchOverflow[%d]:cnt=%d,next=%d\n",
	    tid,mapid2,flip2,minSite2,maxSite2,MatchCnt,IndexCnt, MatchOverflowNext-MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)-1,i,p[i].cnt,p[i].next);
	  fflush(stdout);

	  assert(0 < p[i].cnt && p[i].cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
        }
      }
      TrueCnt2 += p[i].cnt;
    }
    if(HashGC && minSite2 > 0){/* subtract p[i].cnt from freelist starting at MatchOverflowFree */
      for(int i = MatchOverflowFree; i > 0; i = MatchOverflow[i].next){
	if(DEBUG>=2) assert(i > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	if(DEBUG>=2) assert(MatchOverflow[i].cnt==0);
	TrueCnt3 += MatchOverflow[i].cnt;
	FreeCnt++;
	if((VERB>=3 && mapid2==MAPID2 && flip2==FLIP2) || (DEBUG && FreeCnt >= MatchOverflowNext)){
          printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d: i=%d,next=%d,p[i].cnt=%u:IndexCnt=%d,FreeCnt=%d(Free=%d),TrueCnt=%lld,%lld,%lld,MatchOverflowNext=%d,MatchOverflowMax=%d\n",
		 tid,mapid2,flip2,minSite2,maxSite2,i,MatchOverflow[i].next,MatchOverflow[i].cnt,IndexCnt,FreeCnt,MatchOverflowFree,TrueCnt1,TrueCnt2,TrueCnt3,MatchOverflowNext,MatchOverflowMax);
	  fflush(stdout);
	  if(DEBUG) assert(!(FreeCnt >= MatchOverflowNext + 100));
	}
      }
      if(DEBUG) assert(!(FreeCnt >= MatchOverflowNext));
    }

    if(VERB>=2 || TrueCnt1+TrueCnt2-TrueCnt3 != MatchCnt){
      #pragma omp critical
      {
	printf("tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d:MatchCnt=%lld,TrueCnt=%lld(%lld+%lld-%lld),IndexCnt=%d,FreeCnt=%d, Overflow Lines=%d\n",
	       tid,mapid2,flip2,minSite2,maxSite2,MatchCnt, TrueCnt1+TrueCnt2-TrueCnt3, TrueCnt1, TrueCnt2, TrueCnt3,IndexCnt,FreeCnt, MatchOverflowNext-MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)-1);
	fflush(stdout);
	assert(TrueCnt1+TrueCnt2-TrueCnt3 == MatchCnt);
      }
    }
  }

  /* compute largest offset spread value with non-zero OFFSET_WT[] */
  int max_spread = 0;
  if(OFFSET_SPREAD){
    for(int delta = 1;  delta <= OffsetSpread; delta++)
      if(OFFSET_WT[delta])
	max_spread = delta;
    if(DEBUG) assert(OFFSET_WT[0] >= 1);
    int wt_sum = OFFSET_WT[0];
    for(int delta = 1; delta <= max_spread; delta++){
      if(DEBUG) assert(OFFSET_WT[delta] > 0);
      wt_sum += 2*OFFSET_WT[delta];
    }
    if(DEBUG) assert(wt_sum <= 255);
  }

  if(DEBUG && !(MatchCnt >= 0)){
    printf("MatchCnt = %lld\n",MatchCnt);
    fflush(stdout);
    assert(MatchCnt >= 0);
  }

  if(MatchCnt*(1+2*max_spread) > MatchBufMax){/* If minSite2 || maxSite2 only a subset of the Match Table will be copied, but the difference is small */
    if(MatchBufMax > 0){
      // WAS delete [] MatchBuf;
      size_t memsiz = (sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>) * MatchBufMax + PAGE-1) & ~(PAGE-1);
      if(munmap(MatchBuf, memsiz)){
	int eno = errno;
	char *err = strerror(eno);

        #pragma omp critical
        {
	  printf("munmap(%p,%lu) failed: errno=%d:%s\n",MatchBuf, memsiz, eno, err);
	  dumpmemmap();
	  fflush(stdout);exit(1);
	}
      }

      if(hashtable->MatchMemReserve)
	hashtable->MatchMemReserve->ResFree(MatchBufMax * sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), 1, this);
    }

    MatchBufMax = MatchCnt * (1LL + 2LL * max_spread);
    if(STAT && MatchBufMax > MatchBufHWM){
      if(VERB>=2){
	printf("MatchFind:mapid2=%d(id=%lld),flip2=%d:MatchCnt=%lld:MatchBufHWM=%lld -> %lld, HWMmapid2= %d -> %d\n",mapid2,qmap[mapid2]->id,flip2,MatchCnt, MatchBufHWM, MatchBufMax, HWMmapid, mapid2);
	fflush(stdout);
      }
      MatchBufHWM = MatchBufMax;
      HWMmapid = mapid2;
    }
    if(DEBUG) assert(MatchBufMax >= 0);
    if(hashtable->MatchMemReserve && MatchBufMax > 0){
      size_t memsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * MatchOverflowMax + PAGE-1) & ~(PAGE-1);

      if(VERB>=2){
        #pragma omp critical
	{
	  printf("MatchFind(mapid2=%d,flip2=%d,M=%d):tid=%d,MatchCnt=%lld,MatchBufMax=%lld,max_spread=%d,memsiz=%lu,MatchOverflowNext=%u\n",
		 mapid2,flip2,qmap[mapid2]->numsite[hcolor], tid,MatchCnt,MatchBufMax,max_spread,memsiz,MatchOverflowNext);
	  fflush(stdout);
	}
      }

      hashtable->MatchMemReserve->allocate(MatchBufMax * sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), 1, memsiz, this);
    }

    //WAS MatchBuf = new CBestMatch[MatchBufMax];
    size_t memsiz = (sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>) * MatchBufMax + PAGE-1) & ~(PAGE-1);
    char *newmem = (char *) mmap(NULL, memsiz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if(newmem == MAP_FAILED){
      int eno = errno;
      char *err = strerror(eno);
	    
      #pragma omp critical
      {
        printf("CMatchTable: mmap of %lu bytes failed:errno=%d:%s\n", memsiz,eno,err);
	dumpmemmap();
	fflush(stdout);exit(1);
      }
    }
    MatchBuf = (CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *)newmem;
  }

  FLOAT *X = qmap[mapid2]->site[hcolor];
  int N2 = qmap[mapid2]->numsite[hcolor];
  double OffsetKBinv = 1.0 / OffsetKB;
  if(DEBUG && minSite2 > 0 && !(minSite2 + HashWin <= N2)){
    #pragma omp critical
    {
      printf("MatchFind(mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d):HashWin=%u,N2=%d\n",
	     mapid2,flip2, minSite2,maxSite2,HashWin,N2);
      fflush(stdout);
      assert(minSite2 + HashWin <= N2);
    }
  }
  if(DEBUG && maxSite2 > 0) assert(maxSite2 + HashWin <= N2);
  
  OFFSET_T minOffset2, maxOffset2;
  double maxScaleFactor = HASH_SCALE ? 1.0 + ScaleDelta * ScaleDeltaSize : 1.0;
  double minScaleFactor = 1.0/maxScaleFactor;
  if(!flip2){
    minOffset2 = (minSite2 <= 0) ? -MASK((sizeof(OFFSET_T) >= 4) ? 31 : 15) : (int)floor(X[minSite2] * minScaleFactor * OffsetKBinv + 0.5) - max_spread;
    maxOffset2 = (maxSite2 <= 0) ? MASK((sizeof(OFFSET_T) >= 4) ? 31 : 15) : (int)floor(X[maxSite2] * minScaleFactor * OffsetKBinv + 0.5) - max_spread - 1;
  } else {
    minOffset2 = (maxSite2 <= 0) ? -MASK(sizeof(OFFSET_T) >= 4 ? 31 : 15) : (int)floor(X[N2 + 1 - maxSite2 - HashWin] * maxScaleFactor * OffsetKBinv + 0.5) + max_spread + 1;
    maxOffset2 = (minSite2 <= 0) ? MASK(sizeof(OFFSET_T) >= 4 ? 31 : 15) : (int)floor(X[N2 + 1 - minSite2 - HashWin] * maxScaleFactor * OffsetKBinv + 0.5) + max_spread;
  }

  int minOffset, maxOffset, minOffsetP,maxOffsetP;
  if(VERB>=2){
    minOffsetP = minOffset = MASK(31);
    maxOffsetP = maxOffset = -MASK(31);
  }

  /* start up sequential prefetch of MatchBuf[] */
  register CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *pBuf = MatchBuf;
  _mm_prefetch(CAST pBuf, _MM_HINT_T0);
  _mm_prefetch(CAST &pBuf[8], _MM_HINT_T0);

  if(!HashGC || (minSite2 <= 0 && maxSite2 <= 0)){  /* traverse the Match Table and copy all entries to MatchBuf[0..MatchCnt-1] */
    /* HERE HERE : If MatchCnt is a significant fraction (25% ?) of (1<<MHASH_BITS), sequential access of MatchIndex[] (skipping null entries) may be faster
       Also sequential access of MatchOverflow[0..MatchOverflowMax-1] is always faster (for read only access) */

    register int nindex = MatchFirst;
    register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *q = &MatchIndex[nindex], *p;
    _mm_prefetch(CAST q, _MM_HINT_T0);
    for(register int mindex = nindex; mindex >= 0; mindex = nindex){
      p = q;
      q = &MatchIndex[nindex = p->next];
      _mm_prefetch(CAST q, _MM_HINT_T0);    

      if(DEBUG>=2 && !(p->cnt > 0)){
	#pragma omp critical
	{
	  printf("MatchFind(mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d:MatchFirst=%d,mindex=%d,nindex=%d:p->cnt=%u\n",
		 mapid2,flip2,minSite2,maxSite2,MatchFirst,mindex,nindex,p->cnt);
	  fflush(stdout);
	  assert(p->cnt > 0);
	}
      }
      register int cnt = p->cnt;
      if(cnt > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){/* handle overflow linked list */
	register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *s,*r = &MatchOverflow[cnt];
	_mm_prefetch(CAST r, _MM_HINT_T0);    

	for(register int i = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE); --i >= 0;pBuf++){
	  pBuf->mapid1 = p->mapid1[i];
	  pBuf->offset = p->offset[i];
	  pBuf->score = p->score[i];
	  if(RANGE >= 1)
	    pBuf->MinOffset1[0] = p->Offset1[i];
	  if(HASH_SCALE >= 1)
	    pBuf->scaleID[0] = p->scaleID[i];
	  if(VERB>=2){
	    minOffset = min(minOffset, pBuf->offset);
	    maxOffset = max(maxOffset, pBuf->offset);
	  }
	  if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf->offset - TraceOffset) <= OffsetSpread))){
	    if(RANGE>=1)
	      printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,offset1=%u,score=%u,scaleID=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u)\n",
		     p,i,mindex,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,pBuf->MinOffset1[0], pBuf->score,
		     (HASH_SCALE>=1 && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[pBuf->mapid1],Offsethash[pBuf->offset]);
	    else
	      printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,score=%u,scaleID=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u)\n",
		     p,i,mindex,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,pBuf->score,
		     (HASH_SCALE>=1 && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[pBuf->mapid1],Offsethash[pBuf->offset]);
	  }
	}

	for(;r->next >= 0; r = s){
	  s = &MatchOverflow[r->next];
	  _mm_prefetch(CAST s, _MM_HINT_T0);

	  if(DEBUG>=2 && !(r->cnt > 0 && r->cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE))){
#pragma omp critical
	    {
	      printf("tid=%d:mapid2=%d,flip2=%d:MatchIndex[%d].cnt=%d,r= &MatchOverflow[%ld]:r->cnt=%d\n",
		     tid,mapid2,flip2,mindex,MatchIndex[mindex].cnt,r-MatchOverflow,r->cnt);
	      fflush(stdout);
	      assert(r->cnt > 0 && r->cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	    }
	  }
	  for(register int i = r->cnt; --i >= 0;pBuf++){
	    pBuf->mapid1 = r->mapid1[i];
	    pBuf->offset = r->offset[i];
	    pBuf->score = r->score[i];
	    if(RANGE >= 1)
	      pBuf->MinOffset1[0] = r->Offset1[i];
	    if(HASH_SCALE >= 1)
	      pBuf->scaleID[0] = r->scaleID[i];

	    if(VERB>=2){
	      minOffset = min(minOffset, pBuf->offset);
	      maxOffset = max(maxOffset, pBuf->offset);
	    }
	    if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf->offset - TraceOffset) <= OffsetSpread))){
	      if(RANGE >= 1)
		printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,score=%u,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u)\n",
		       r,i,mindex,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,pBuf->MinOffset1[0],pBuf->score,
		       (HASH_SCALE>=1 && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[pBuf->mapid1],Offsethash[pBuf->offset]);
	      else
		printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,score=%u,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u)\n",
		       r,i,mindex,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,pBuf->score,
		       (HASH_SCALE>=1 && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[pBuf->mapid1],Offsethash[pBuf->offset]);
	    }
	  }
	}

	if(DEBUG>=2) assert(r->cnt > 0 && r->cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	for(register int i = r->cnt; --i >= 0;pBuf++){
	  pBuf->mapid1 = r->mapid1[i];
	  pBuf->offset = r->offset[i];
	  pBuf->score = r->score[i];
	  if(RANGE >= 1)
	    pBuf->MinOffset1[0] = r->Offset1[i];
	  if(HASH_SCALE >= 1)
	    pBuf->scaleID[0] = r->scaleID[i];
	  if(VERB>=2){
	    minOffset = min(minOffset, pBuf->offset);
	    maxOffset = max(maxOffset, pBuf->offset);
	  }
	  if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf->offset - TraceOffset) <= OffsetSpread))){
	    if(RANGE>=1)
	      printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u)\n",
		     r,i,mindex,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,pBuf->MinOffset1[0], pBuf->score,
		     (HASH_SCALE>=1 && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[pBuf->mapid1],Offsethash[pBuf->offset]);
	    else
	      printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u)\n",
		     r,i,mindex,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,pBuf->score,
		     (HASH_SCALE>=1 && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[pBuf->mapid1],Offsethash[pBuf->offset]);
	  }
	}

	continue;
      }

      /* just copy p->cnt entries to MatchBuf[] */
      for(; --cnt >= 0;pBuf++){
	pBuf->mapid1 = p->mapid1[cnt];
	pBuf->offset = p->offset[cnt];
	pBuf->score = p->score[cnt];
	if(RANGE >= 1)
	  pBuf->MinOffset1[0] = p->Offset1[cnt];
	if(HASH_SCALE >= 1)
	  pBuf->scaleID[0] = p->scaleID[cnt];
	if(VERB>=2){
	  minOffset = min(minOffset, pBuf->offset);
	  maxOffset = max(maxOffset, pBuf->offset);
	}
	if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf->offset - TraceOffset) <= OffsetSpread))){
	  if(RANGE >= 1)
	    printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u:score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u)\n",
		   p,cnt,mindex,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset, pBuf->MinOffset1[0], pBuf->score,
		   (HASH_SCALE>=1 && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[pBuf->mapid1],Offsethash[pBuf->offset]);
	  else
	    printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d:score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u)\n",
		   p,cnt,mindex,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,pBuf->score,
		   (HASH_SCALE>=1 && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[pBuf->mapid1],Offsethash[pBuf->offset]);
	}
      }
    }
    if(VERB>=2){
      minOffsetP = minOffset;
      maxOffsetP = maxOffset;
    }
  } else { /* traverse the Match Table and copy only entries with minOffset2 - Offset1flip2[mapid1] <= offset <= maxOffset2 - Offset1flip2[mapid1] to MatchBuf[0..MatchCnt-1] */
    if(DEBUG) assert(HashGC);

    register int nindex = MatchFirst;
    register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *q = &MatchIndex[nindex], *p;
    _mm_prefetch(CAST q, _MM_HINT_T0);

    for(register int mindex = nindex; mindex >= 0; mindex = nindex){
      p = q;
      q = &MatchIndex[nindex = p->next];
      _mm_prefetch(CAST q, _MM_HINT_T0);    

      register int cnt = p->cnt;
      if(DEBUG>=2) assert(cnt > 0);

      if(cnt > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){/* handle overflow linked list */
	register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *s,*r = &MatchOverflow[cnt];
	_mm_prefetch(CAST r, _MM_HINT_T0);    

	for(register int i = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE); --i >= 0;){
	  int mapid1 = p->mapid1[i];
	  OFFSET_T offset = p->offset[i];
	  OFFSET_T Offset2 = offset + Offset1flip2[mapid1];
	  if(VERB>=2){
	    minOffset = min(minOffset, offset);
	    maxOffset = max(maxOffset, offset);
	  }
	  if(Offset2 < minOffset2 || Offset2 > maxOffset2){
	    if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	      if(RANGE >= 1)
		printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,score=%u,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Skipping Offset2=%d(min=%d,max=%d)\n",
		       p,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,p->Offset1[i],p->score[i],
		       (HASH_SCALE && hashScaleDelta) ? p->scaleID[i] : 0,Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	      else
		printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Skipping Offset2=%d(min=%d,max=%d)\n",
		       p,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,p->score[i],
		       (HASH_SCALE && hashScaleDelta) ? p->scaleID[i] : 0,Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	    }
	    continue;
	  }

	  pBuf->mapid1 = mapid1;
	  pBuf->offset = offset;
	  pBuf->score = p->score[i];
	  if(RANGE >= 1){
	    pBuf->MinOffset1[0] = p->Offset1[i];
	    pBuf->MaxOffset1[0] = pBuf->score;
	  }
	  if(HASH_SCALE >= 1)
	    pBuf->scaleID[0] = p->scaleID[i];
	  if(VERB>=2){
	    if(flip2==1 && offset < minOffsetP && offset < minOffset2){
	      int M = rmap[mapid1]->numsite[hcolor];
	      FLOAT *Y = rmap[mapid1]->site[hcolor];
	      printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d:score=%d,sc=%u:Offset2=%d..%d,Offset1=%d,M=%d,Y[M]=%0.4f,Y[M+1]=%0.4f\n",
		     p,cnt,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,pBuf->score,
		     (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0,
		     minOffset2,maxOffset2,Offset1flip2[mapid1],M,Y[M],Y[M+1]);
	      fflush(stdout);
	    }
	    minOffsetP = min(minOffsetP, offset);
	    maxOffsetP = max(maxOffsetP, offset);
	  }
	  if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	    if(RANGE >= 1 )
	      printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Offset2=%d(min=%d,max=%d)\n",
		     p,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,pBuf->MinOffset1[0],pBuf->score,
		     (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	    else
	      printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Offset2=%d(min=%d,max=%d)\n",
		     p,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,pBuf->score,
		     (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	  }
	  pBuf++;
	}

	for(;r->next >= 0; r = s){
	  s = &MatchOverflow[r->next];
	  _mm_prefetch(CAST s, _MM_HINT_T0);

	  if(DEBUG>=2 && !(r->cnt > 0 && r->cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE))){
            #pragma omp critical
	    {
	      printf("tid=%d:mapid2=%d,flip2=%d:MatchIndex[%d].cnt=%d,r= &MatchOverflow[%ld]:r->cnt=%d\n",
		     tid,mapid2,flip2,mindex,MatchIndex[mindex].cnt,r-MatchOverflow,r->cnt);
	      fflush(stdout);
	      assert(r->cnt > 0 && r->cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	    }
	  }
	  for(register int i = r->cnt; --i >= 0;){
	    int mapid1 = r->mapid1[i];
	    OFFSET_T offset = r->offset[i];
	    OFFSET_T Offset2 = offset + Offset1flip2[mapid1];
	    if(VERB>=2){
	      minOffset = min(minOffset, offset);
	      maxOffset = max(maxOffset, offset);
	    }
	    if(Offset2 < minOffset2 || Offset2 > maxOffset2){
	      if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	        if(RANGE >= 1)
		  printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,offset1=%u,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Skipping Offset2=%d(min=%d,max=%d)\n",
			 r,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,r->Offset1[i],r->score[i],
			 (HASH_SCALE && hashScaleDelta) ? r->scaleID[i] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
		else
		  printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Skipping Offset2=%d(min=%d,max=%d)\n",
			 r,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,r->score[i],
			 (HASH_SCALE && hashScaleDelta) ? r->scaleID[i] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
		fflush(stdout);
	      }
	      continue;
	    }

	    pBuf->mapid1 = mapid1;
	    pBuf->offset = offset;
	    pBuf->score = r->score[i];
	    if(RANGE >= 1){
	      pBuf->MinOffset1[0] = r->Offset1[i];
	      pBuf->MaxOffset1[0] = pBuf->score;
	    }
	    if(HASH_SCALE >= 1)
	      pBuf->scaleID[0] = r->scaleID[i];
	    if(VERB>=2){
	      if(flip2==1 && offset < minOffsetP && offset < minOffset2){
		int M = rmap[mapid1]->numsite[hcolor];
		FLOAT *Y = rmap[mapid1]->site[hcolor];
		printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d:score=%d,sc=%u:Offset2=%d..%d,Offset1=%d,M=%d,Y[M]=%0.4f,Y[M+1]=%0.4f\n",
		       p,cnt,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,pBuf->score,
		       (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, minOffset2,maxOffset2,Offset1flip2[mapid1],M,Y[M],Y[M+1]);
		fflush(stdout);
	      }
	      minOffsetP = min(minOffsetP, offset);
	      maxOffsetP = max(maxOffsetP, offset);
	    }
	    if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	      if(RANGE >= 1)
		printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Offset2=%d(min=%d,max=%d)\n",
		       r,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,pBuf->MinOffset1[0],pBuf->score,
		       (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	      else
		printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Offset2=%d(min=%d,max=%d)\n",
		       r,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,pBuf->score,
		       (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	    }
	    pBuf++;
	  }
	}

	if(DEBUG>=2) assert(r->cnt > 0 && r->cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));

	for(register int i = r->cnt; --i >= 0;){
	  int mapid1 = r->mapid1[i];
	  OFFSET_T offset = r->offset[i];
	  OFFSET_T Offset2 = offset + Offset1flip2[mapid1];
	  if(VERB>=2){
	    minOffset = min(minOffset, offset);
	    maxOffset = max(maxOffset, offset);
	  }
	  if(Offset2 < minOffset2 || Offset2 > maxOffset2){
	    if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	      if(RANGE >= 1)
		printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,score=%d (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Skipping Offset2=%d(min=%d,max=%d)\n",
		       r,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,r->Offset1[i],r->score[i],Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	      else
		printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,score=%d (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Skipping Offset2=%d(min=%d,max=%d)\n",
		       r,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,r->score[i],Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	      fflush(stdout);
	    }
	    continue;
	  }

	  pBuf->mapid1 = mapid1;
	  pBuf->offset = offset;
	  pBuf->score = r->score[i];
	  if(RANGE >= 1){
	    pBuf->MinOffset1[0] = r->Offset1[i];
	    pBuf->MaxOffset1[0] = pBuf->score;
	  }
	  if(HASH_SCALE >= 1)
	    pBuf->scaleID[0] = r->scaleID[i];
	  if(VERB>=2){
	    if(flip2==1 && offset < minOffsetP && offset < minOffset2){
	      int M = rmap[mapid1]->numsite[hcolor];
	      FLOAT *Y = rmap[mapid1]->site[hcolor];
	      printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d:score=%d,sc=%u:Offset2=%d..%d,Offset1=%d,M=%d,Y[M]=%0.4f,Y[M+1]=%0.4f\n",
		     p,cnt,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,pBuf->score,
		       (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, minOffset2,maxOffset2,Offset1flip2[mapid1],M,Y[M],Y[M+1]);
	      fflush(stdout);
	    }
	    minOffsetP = min(minOffsetP, offset);
	    maxOffsetP = max(maxOffsetP, offset);
	  }
	  if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	    if(RANGE >= 1)
	      printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Offset2=%d(min=%d,max=%d)\n",
		     r,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,pBuf->MinOffset1[0],pBuf->score,
		     (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	    else
	      printf("\t r=%p,i=%d,mindex=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Offset2=%d(min=%d,max=%d)\n",
		     r,i,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,pBuf->score,
		     (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	    fflush(stdout);
	  }
	  pBuf++;
	}

	continue;
      }

      /* just copy p->cnt entries to MatchBuf[] */
      for(; --cnt >= 0;){
	int mapid1 = p->mapid1[cnt];
	OFFSET_T offset = p->offset[cnt];
	OFFSET_T Offset2 = offset + Offset1flip2[mapid1];
	if(VERB>=2){
	  minOffset = min(minOffset, offset);
	  maxOffset = max(maxOffset, offset);
	}
	if(Offset2 < minOffset2 || Offset2 > maxOffset2){
	  if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	    if(RANGE >= 1)
	      printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u:score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Skipping Offset2=%d(min=%d,max=%d)\n",
		     p,cnt,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,p->Offset1[cnt],p->score[cnt],
		     (HASH_SCALE && hashScaleDelta) ? p->scaleID[cnt] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	    else
	      printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d:score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Skipping Offset2=%d(min=%d,max=%d)\n",
		     p,cnt,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,p->score[cnt],
		     (HASH_SCALE && hashScaleDelta) ? p->scaleID[cnt] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	  }
	  continue;
	}

	pBuf->mapid1 = p->mapid1[cnt];
	pBuf->offset = p->offset[cnt];
	pBuf->score = p->score[cnt];
	if(RANGE>=1){
	  pBuf->MinOffset1[0] = p->Offset1[cnt];
	  pBuf->MaxOffset1[0] = pBuf->score;
	}
	if(HASH_SCALE >= 1)
	  pBuf->scaleID[0] = p->scaleID[cnt];
	if(VERB>=2){
	  if(flip2==1 && offset < minOffsetP && offset < minOffset2){
	    int M = rmap[mapid1]->numsite[hcolor];
	    FLOAT *Y = rmap[mapid1]->site[hcolor];
	    printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d:score=%d,sc=%u:Offset2=%d..%d,Offset1=%d,M=%d,Y[M]=%0.4f,Y[M+1]=%0.4f\n",
		   p,cnt,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,pBuf->score,
		   (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, minOffset2,maxOffset2,Offset1flip2[mapid1],M,Y[M],Y[M+1]);
	    fflush(stdout);
	  }
	  minOffsetP = min(minOffsetP, offset);
	  maxOffsetP = max(maxOffsetP, offset);
	}
	if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	  if(RANGE >= 1)
	    printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u:score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Offset2=%d(min=%d,max=%d)\n",
		   p,cnt,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,pBuf->MinOffset1[0],
		   (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, p->score[cnt], Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	  else
	    printf("\t p=%p,i=%d,mindex=%d:mapid1=%d(%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d:score=%d,sc=%u (Maphash2[mapid1]=%u,Offsethash[offset]=%u):Offset2=%d(min=%d,max=%d)\n",
		   p,cnt,mindex,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,p->score[cnt],
		   (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, Maphash2[mapid1],Offsethash[offset],Offset2,minOffset2,maxOffset2);
	}
	pBuf++;
      }
    }
  }

  if(DEBUG>=1+RELEASE && HashGC && (minSite2 > 0 || maxSite2 > 0)) assert(pBuf <= &MatchBuf[MatchCnt]);

  if(DEBUG>=1+RELEASE && !(HashGC && (minSite2 > 0 || maxSite2 > 0))) assert(pBuf == &MatchBuf[MatchCnt]);

  register long long i;

  long long nMatchCnt = pBuf - MatchBuf;/* May be less than MatchCnt if (minSite2 > 0 || maxSite2 > 0) */
  long long orignMatchCnt;
  if(RANGE>=1)
    orignMatchCnt = nMatchCnt;

  if(VERB>=2){
    #pragma omp critical
    {
      printf("MatchFind(mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d):tid=%d: Finished initial copy of MatchTable to MatchBuf:MatchCnt=%lld,nMatchCnt=%lld,offset=%d..%d(%d..%d),Offset2=%d..%d\n",
	     mapid2,flip2,minSite2,maxSite2,N2,tid,MatchCnt,nMatchCnt,minOffset,maxOffset,minOffsetP,maxOffsetP,minOffset2,maxOffset2);
      fflush(stdout);
    }
  }

  if(OFFSET_SPREAD){/* update entries in MatchBuf[] Table with nearby offset values from MatchTable */
    int wt_sum = OFFSET_WT[0];
    for(int delta = 1; delta <= max_spread; delta++){
      if(DEBUG) assert(OFFSET_WT[delta] > 0);
      wt_sum += 2*OFFSET_WT[delta];
    }
    if(DEBUG) assert(wt_sum <= 255);
    if(OFFSET_WT[0] > 1){/* multiply scores in MatchBuf[] by OFFSET_WT[0] */
      if(VERB>=2 || (TRACE && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2))
	printf("OFFSET_WT[0]=%d:mapid2=%d(id=%lld),flip2=%d: multiplying all scores by %d\n",OFFSET_WT[0],mapid2,qmap[mapid2]->id,flip2,OFFSET_WT[0]);

      register short wt0 = OFFSET_WT[0];
      pBuf = MatchBuf;
      _mm_prefetch(CAST pBuf, _MM_HINT_T0);
      _mm_prefetch(CAST &pBuf[8], _MM_HINT_T0);
      for(i = 0; i < nMatchCnt; i++)
	pBuf[i].score *= wt0;
    }

    /* In order to handle scores with a pattern of ...,2,0,2,... with OFFSET_WT[1]=1, it is necessarily to add 0 score entries within max_spread of non-zero scores to MatchBuf[] */
    register CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *qBuf = &MatchBuf[nMatchCnt];/* pointer to beyond last entry appended to MatchBuf[0 .. nMatchCnt-1] with zero scores (see OFFSET_FILL) */
    for(register OFFSET_T delta = -max_spread; delta <= max_spread; delta++){ /* Look for MatchTable entries at pBuf[i].offset + delta */
      if(!delta)
	continue;
      register short wt = OFFSET_WT[abs((int)delta)];
      if(wt <= 0)
	continue;

      if(VERB>=2 || (TRACE && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2))
	printf("OFFSET_WT[%d]=%d:mapid2=%d(id=%lld),flip2=%d:nMatchCnt=%lld\n",delta,wt,mapid2,qmap[mapid2]->id,flip2,nMatchCnt);

      pBuf = MatchBuf;
      // HERE : try using small circular buffers for mapid1, offset+delta, nindex to avoid have to recompute them more than once
      _mm_prefetch(CAST pBuf, _MM_HINT_T0);
      for(i = 0; i < PREFETCH_D && i < nMatchCnt; i++){/* 1st prefetch loop */
	_mm_prefetch(CAST &Maphash2[pBuf[i].mapid1], _MM_HINT_T0);
	_mm_prefetch(CAST &Offsethash[pBuf[i].offset+delta], _MM_HINT_T0);
	if(RANGE>=1)
	  _mm_prefetch(CAST &Offset1hash[pBuf[i].MinOffset1[0]], _MM_HINT_T0);
	if(HASH_SCALE>=1)
	  _mm_prefetch(CAST &scaleIDhash[pBuf[i].scaleID[0]], _MM_HINT_T0);
      }

      for(i = 0; i < PREFETCH_D && i < nMatchCnt - PREFETCH_D; i++){/* 2nd prefetch loop */
	_mm_prefetch(CAST &Maphash2[pBuf[i+PREFETCH_D].mapid1], _MM_HINT_T0);
	_mm_prefetch(CAST &Offsethash[pBuf[i+PREFETCH_D].offset+delta], _MM_HINT_T0);
	if(RANGE >= 1)
	  _mm_prefetch(CAST &Offset1hash[pBuf[i+PREFETCH_D].MinOffset1[0]], _MM_HINT_T0);
	if(HASH_SCALE >= 1)
	  _mm_prefetch(CAST &scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]], _MM_HINT_T0);

	register unsigned int nindex = Maphash2[pBuf[i].mapid1] ^ Offsethash[pBuf[i].offset+delta];
	if(RANGE >= 1)
	  nindex ^= Offset1hash[pBuf[i].MinOffset1[0]];
	if(HASH_SCALE >= 1)
	  nindex ^= scaleIDhash[pBuf[i].scaleID[0]];
	_mm_prefetch(CAST &MatchIndex[nindex], _MM_HINT_T0);
      }

      if(minSite2 <= 0 && maxSite2 <= 0){/* original code */
	for(i = 0; i < nMatchCnt - PREFETCH_D2; i++){/* main loop */
	  if(PREFETCH_D){
	    _mm_prefetch(CAST &Maphash2[pBuf[i+PREFETCH_D2].mapid1], _MM_HINT_T0);
	    _mm_prefetch(CAST &Offsethash[pBuf[i+PREFETCH_D2].offset+delta], _MM_HINT_T0);
	    if(RANGE >= 1)
	      _mm_prefetch(CAST &Offset1hash[pBuf[i+PREFETCH_D2].MinOffset1[0]], _MM_HINT_T0);
	    if(HASH_SCALE >= 1)
	      _mm_prefetch(CAST &scaleIDhash[pBuf[i+PREFETCH_D2].scaleID[0]], _MM_HINT_T0);

	    unsigned int nindex = Maphash2[pBuf[i+PREFETCH_D].mapid1] ^ Offsethash[pBuf[i+PREFETCH_D].offset+delta];
	    if(RANGE >= 1)
	      nindex ^= Offset1hash[pBuf[i+PREFETCH_D].MinOffset1[0]];
	    if(HASH_SCALE >= 1)
	      nindex ^= scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]];
	    _mm_prefetch(CAST &MatchIndex[nindex], _MM_HINT_T0);
	  }

	  register int mapid1 = pBuf[i].mapid1;
	  register OFFSET_T offset = pBuf[i].offset + delta;
	  register short score = MatchScoreFind(mapid1,offset,(RANGE>=1) ? pBuf[i].MinOffset1[0] : 0, (HASH_SCALE >=1) ? pBuf[i].scaleID[0] : 0);

	  if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf[i].offset - TraceOffset) <= OffsetSpread))){
	    if(RANGE >= 1){
	      if(score <= (OFFSET_ZEROFIX ? -1 : 0))
		printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%d,sc=%u,delta=%d:score %d added at offset=%d, MatchBuf[%ld]\n",
		       mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,pBuf[i].MinOffset1[0],
		       (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,score,offset,qBuf-MatchBuf);
	      else
		printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%d,sc=%u,delta=%d:score=%d -> %d from offset=%d\n",
		       mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,pBuf[i].MinOffset1[0],
		       (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,pBuf[i].score,pBuf[i].score + score * wt, offset);
	    } else {
	      if(score <= (OFFSET_ZEROFIX ? -1 : 0))
		printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,sc=%u,delta=%d:score %d added at offset=%d, MatchBuf[%ld]\n",
		       mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,
		       (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,score,offset,qBuf-MatchBuf);
	      else
		printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,sc=%u,delta=%d:score=%d -> %d from offset=%d\n",
		       mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,
		       (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,pBuf[i].score,pBuf[i].score + score * wt, offset);
	    }
	  }
	  if(score > 0) {
	    pBuf[i].score += score * wt;
	  } else if(OFFSET_FILL && score <= (OFFSET_ZEROFIX ? -1 : 0)){/* NOTE : this can generate duplicate entries with score == 0, unless OFFSET_ZEROFIX */
	    qBuf->mapid1 = mapid1;
	    qBuf->offset = offset;
	    qBuf->score = 0;
	    if(HASH_SCALE>=1)
	      qBuf->scaleID[0] = pBuf[i].scaleID[0];
	    if(RANGE >= 1)
	      qBuf->MinOffset1[0] = pBuf[i].MinOffset1[0];
	    if(OFFSET_ZEROFIX)
	      (void)MatchInsert(mapid1,offset,0,mapid2,flip2, 0, HASH_SCALE>=1 ? qBuf->scaleID[0] : 0, RANGE>=1 ? qBuf->MinOffset1[0] : 0, 0);
	    qBuf++;
	  }
	}
      } else {/* restrict zero score entries added to qBuf[] to those with minOffset2 - Offset1flip2[mapid1] <= offset <= maxOffset2 - Offset1flip2[mapid1] */
	for(i = 0; i < nMatchCnt - PREFETCH_D2; i++){/* main loop */
	  if(PREFETCH_D){
	    _mm_prefetch(CAST &Maphash2[pBuf[i+PREFETCH_D2].mapid1], _MM_HINT_T0);
	    _mm_prefetch(CAST &Offsethash[pBuf[i+PREFETCH_D2].offset+delta], _MM_HINT_T0);
	    if(RANGE >= 1)
	      _mm_prefetch(CAST &Offset1hash[pBuf[i+PREFETCH_D2].MinOffset1[0]], _MM_HINT_T0);
	    if(HASH_SCALE >= 1)
	      _mm_prefetch(CAST &scaleIDhash[pBuf[i+PREFETCH_D2].scaleID[0]], _MM_HINT_T0);	    

	    unsigned int nindex = Maphash2[pBuf[i+PREFETCH_D].mapid1] ^ Offsethash[pBuf[i+PREFETCH_D].offset+delta];
	    if(RANGE >= 1)
	      nindex ^= Offset1hash[pBuf[i+PREFETCH_D].MinOffset1[0]];
	    if(HASH_SCALE >= 1)
	      nindex ^= scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]];
	    _mm_prefetch(CAST &MatchIndex[nindex], _MM_HINT_T0);
	  }

	  register int mapid1 = pBuf[i].mapid1;
	  register OFFSET_T offset = pBuf[i].offset + delta;
	  register short score = MatchScoreFind(mapid1,offset,RANGE>=1 ? pBuf[i].MinOffset1[0] : 0, HASH_SCALE>=1 ? pBuf[i].scaleID[0] : 0);

	  if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf[i].offset - TraceOffset) <= OffsetSpread))){
	    if(score <= 0)
	      printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,sc=%u,delta=%d:score %d added at offset=%d, MatchBuf[%ld]\n",
		     mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset, RANGE>=1 ? pBuf[i].MinOffset1[0]:0,
		     (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0,  delta,score,offset,qBuf-MatchBuf);
	    else
	      printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1 -> %u,sc=%u,delta=%d:score=%d to %d from offset=%d\n",
		     mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,RANGE>=1 ? pBuf[i].MinOffset1[0]:0,
		     (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,pBuf[i].score,pBuf[i].score + score * wt, offset);
	  }

	  if(score > 0){
	    pBuf[i].score += score * wt;
	  } else if(OFFSET_FILL && score <= (OFFSET_ZEROFIX ? -1 : 0)){/* NOTE : this can generate duplicate entries with score == 0, unless OFFSET_ZEROFIX */
	    OFFSET_T Offset2 = offset + Offset1flip2[mapid1];
	    if(Offset2 < minOffset2 || Offset2 > maxOffset2)
	      continue;

	    qBuf->mapid1 = mapid1;
	    qBuf->offset = offset;
	    qBuf->score = 0;
	    if(HASH_SCALE >= 1)
	      qBuf->scaleID[0] = pBuf[i].scaleID[0];
	    if(RANGE >= 1)
	      qBuf->MinOffset1[0] = pBuf[i].MinOffset1[0];
	    if(OFFSET_ZEROFIX)
	      (void)MatchInsert(mapid1,offset,0,mapid2,flip2, 0, HASH_SCALE>=1 ? qBuf->scaleID[0] : 0, RANGE>=1 ? qBuf->MinOffset1[0] : 0, 0);

	    qBuf++;
	  }
	}
      }

      for(; i < nMatchCnt - PREFETCH_D; i++){/* 1st post loop */
	unsigned int nindex = Maphash2[pBuf[i+PREFETCH_D].mapid1] ^ Offsethash[pBuf[i+PREFETCH_D].offset+delta];
	if(RANGE >= 1)
	  nindex ^= Offset1hash[pBuf[i+PREFETCH_D].MinOffset1[0]];
	if(HASH_SCALE >= 1)
	  nindex ^= scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]];
	_mm_prefetch(CAST &MatchIndex[nindex], _MM_HINT_T0);

	register int mapid1 = pBuf[i].mapid1;
	register OFFSET_T offset = pBuf[i].offset + delta;
	register short score = MatchScoreFind(mapid1,offset,RANGE>=1 ? pBuf[i].MinOffset1[0]:0, HASH_SCALE>=1 ? pBuf[i].scaleID[0] : 0);

	if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf[i].offset - TraceOffset) <= OffsetSpread))){
	  if(score <= (OFFSET_ZEROFIX ? -1 : 0))
	    printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,sc=%u,delta=%d:score %d added at offset=%d,MatchBuf[%ld]\n",
		   mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,RANGE>=1 ? pBuf[i].MinOffset1[0]:0,
		   (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,score,offset,qBuf-MatchBuf);
	  else
	    printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,sc=%u,delta=%d:score=%d to %d from offset=%d\n",
		   mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,RANGE>=1 ? pBuf[i].MinOffset1[0]:0,
		   (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,pBuf[i].score,pBuf[i].score + score * wt, offset);
	}

	if(score > 0){
	  pBuf[i].score += score * wt;
	} else if(OFFSET_FILL && score <= (OFFSET_ZEROFIX ? -1 : 0)){/* NOTE : this can generate duplicate entries with score == 0, unless OFFSET_ZEROFIX */
	  OFFSET_T Offset2 = offset + Offset1flip2[mapid1];
	  if(Offset2 < minOffset2 || Offset2 > maxOffset2)
	    continue;

	  qBuf->mapid1 = mapid1;
	  qBuf->offset = offset;
	  qBuf->score = 0;
	  if(HASH_SCALE >= 1)
	    qBuf->scaleID[0] = pBuf[i].scaleID[0];
	  if(RANGE >= 1)
	    qBuf->MinOffset1[0] = pBuf[i].MinOffset1[0];
	  if(OFFSET_ZEROFIX)
	    (void)MatchInsert(mapid1,offset,0,mapid2,flip2, 0, HASH_SCALE>=1 ? qBuf->scaleID[0]:0, RANGE>=1 ? qBuf->MinOffset1[0] : 0,0);

	  qBuf++;
	}
      }

      for(; i < nMatchCnt; i++){/* 2nd post loop */
	register int mapid1 = pBuf[i].mapid1;
	register OFFSET_T offset = pBuf[i].offset + delta;
	register short score = MatchScoreFind(mapid1,offset,RANGE>=1 ? pBuf[i].MinOffset1[0]:0,HASH_SCALE>=1?pBuf[i].scaleID[0]:0);

	if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf[i].offset - TraceOffset) <= OffsetSpread))){
	  if(score <= (OFFSET_ZEROFIX ? -1 : 0))
 	    printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,sc=%d,delta=%d:score %d added at offset=%d,MatchBuf[%ld]\n",
		   mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,RANGE>=1 ? pBuf[i].MinOffset1[0]:0,
		   (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,score,offset,qBuf-MatchBuf);
	  else
	    printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,sc=%d,delta=%d:score=%d to %d from offset=%d\n",
		   mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,RANGE>=1 ? pBuf[i].MinOffset1[0]:0,
		   (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,pBuf[i].score,pBuf[i].score + score * wt,offset);
	}

	if(score > 0){
	  pBuf[i].score += score * wt;
	} else if(OFFSET_FILL && score <= (OFFSET_ZEROFIX ? -1 : 0)){/* NOTE : this can generate duplicate entries with score == 0, unless OFFSET_ZEROFIX */
	  OFFSET_T Offset2 = offset + Offset1flip2[mapid1];
	  if(Offset2 < minOffset2 || Offset2 > maxOffset2)
	    continue;

	  qBuf->mapid1 = mapid1;
	  qBuf->offset = offset;
	  qBuf->score = 0;
	  if(HASH_SCALE >= 1)
	    qBuf->scaleID[0] = pBuf[i].scaleID[0];
	  if(RANGE >= 1)
	    qBuf->MinOffset1[0] = pBuf[i].MinOffset1[0];
	  if(OFFSET_ZEROFIX)
	    (void)MatchInsert(mapid1,offset,0,mapid2,flip2, 0, HASH_SCALE >= 1 ? qBuf->scaleID[0] : 0, RANGE >= 1 ? qBuf->MinOffset1[0] : 0, 0);

	  qBuf++;
	}
      }
    }

    long long MatchCntNull = qBuf - MatchBuf;

    if(VERB>=2 || (TRACE && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2)){
      #pragma omp critical
      {
	printf("MatchFind(mapid2=%d):tid=%d: Updated MatchBuf[] with empty slots:MatchCntNull=%lld: wall time=%0.6f\n",mapid2,tid,MatchCntNull,wtime());
	fflush(stdout);
      }
    }

    if(OFFSET_FILL && MatchCntNull > nMatchCnt){
      if(DEBUG) assert(MatchCntNull <= MatchBufMax);
      for(register OFFSET_T delta = -max_spread; delta <= max_spread; delta++){ /* Look for MatchTable entries at pBuf[nMatchCnt..MatchCntNull-1].offset + delta */
	if(!delta)
	  continue;
	register short wt = OFFSET_WT[abs((int)delta)];
	if(wt <= 0)
	  continue;

	if(VERB>=2 || (TRACE && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2))
	  printf("OFFSET_WT[%d]=%d:mapid2=%d(id=%lld),flip2=%d:MatchCnt=%lld\n",delta,wt,mapid2,qmap[mapid2]->id,flip2,MatchCnt);

	// HERE : try using small circular buffers for mapid1, offset+delta, nindex to avoid have to recompute them more than once
	_mm_prefetch(CAST pBuf, _MM_HINT_T0);
	for(i = nMatchCnt; i < nMatchCnt + PREFETCH_D && i < MatchCntNull; i++){/* 1st prefetch loop */
	  _mm_prefetch(CAST &Maphash2[pBuf[i].mapid1], _MM_HINT_T0);
	  _mm_prefetch(CAST &Offsethash[pBuf[i].offset+delta], _MM_HINT_T0);
	  if(RANGE >= 1)
	    _mm_prefetch(CAST &Offset1hash[pBuf[i].MinOffset1[0]], _MM_HINT_T0);	  
	  if(HASH_SCALE >= 1)
	    _mm_prefetch(CAST &scaleIDhash[pBuf[i].scaleID[0]], _MM_HINT_T0);	  
	}
	for(i = nMatchCnt; i < nMatchCnt + PREFETCH_D && i < MatchCntNull-PREFETCH_D; i++){/* 2nd prefetch loop */
	  _mm_prefetch(CAST &Maphash2[pBuf[i+PREFETCH_D].mapid1], _MM_HINT_T0);
	  _mm_prefetch(CAST &Offsethash[pBuf[i+PREFETCH_D].offset+delta], _MM_HINT_T0);
	  if(RANGE >= 1)
	    _mm_prefetch(CAST &Offset1hash[pBuf[i+PREFETCH_D].MinOffset1[0]], _MM_HINT_T0);
	  if(HASH_SCALE >= 1)
	    _mm_prefetch(CAST &scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]], _MM_HINT_T0);

	  register unsigned int nindex = Maphash2[pBuf[i].mapid1] ^ Offsethash[pBuf[i].offset+delta];
	  if(RANGE >= 1)
	    nindex ^= Offset1hash[pBuf[i].MinOffset1[0]];
	  if(HASH_SCALE >= 1)
	    nindex ^= scaleIDhash[pBuf[i].scaleID[0]];
	  _mm_prefetch(CAST &MatchIndex[nindex], _MM_HINT_T0);
	}

	for(i = nMatchCnt; i < MatchCntNull-PREFETCH_D2; i++){/* main loop */
	  if(PREFETCH_D){
	    _mm_prefetch(CAST &Maphash2[pBuf[i+PREFETCH_D2].mapid1], _MM_HINT_T0);
	    _mm_prefetch(CAST &Offsethash[pBuf[i+PREFETCH_D2].offset+delta], _MM_HINT_T0);
	    if(RANGE >= 1)
	      _mm_prefetch(CAST &Offset1hash[pBuf[i+PREFETCH_D2].MinOffset1[0]], _MM_HINT_T0);
	    if(HASH_SCALE >= 1)
	      _mm_prefetch(CAST &scaleIDhash[pBuf[i+PREFETCH_D2].scaleID[0]], _MM_HINT_T0);

	    unsigned int nindex = Maphash2[pBuf[i+PREFETCH_D].mapid1] ^ Offsethash[pBuf[i+PREFETCH_D].offset+delta];
	    if(RANGE >= 1)
	      nindex ^= Offset1hash[pBuf[i+PREFETCH_D].MinOffset1[0]];
	    if(HASH_SCALE >= 1)
	      nindex ^= scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]];
	    _mm_prefetch(CAST &MatchIndex[nindex], _MM_HINT_T0);
	  }

	  register int mapid1 = pBuf[i].mapid1;
	  register OFFSET_T offset = pBuf[i].offset + delta;
	  register short score = MatchScoreFind(mapid1,offset,RANGE>=1 ? pBuf[i].MinOffset1[0]:0,HASH_SCALE>=1 ? pBuf[i].scaleID[0]:0);

	  if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf[i].offset - TraceOffset) <= OffsetSpread)
	    && (!OFFSET_ZEROFIX || score > 0)))
	    printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,sc=%d,delta=%d:score=%d to %d from offset=%d\n",
		   mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,RANGE>=1 ? pBuf[i].MinOffset1[0] : 0,
		     (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,pBuf[i].score,pBuf[i].score + score * wt, offset);

	  if(!OFFSET_ZEROFIX || score > 0)
	    pBuf[i].score += score * wt;
	}
	for(; i < MatchCntNull-PREFETCH_D; i++){/* 1st post loop */
	  unsigned int nindex = Maphash2[pBuf[i+PREFETCH_D].mapid1] ^ Offsethash[pBuf[i+PREFETCH_D].offset+delta];
	  if(RANGE >= 1)
	    nindex ^= Offset1hash[pBuf[i+PREFETCH_D].MinOffset1[0]];
	  if(HASH_SCALE >= 1)
	    nindex ^= scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]];
	  _mm_prefetch(CAST &MatchIndex[nindex], _MM_HINT_T0);

	  register int mapid1 = pBuf[i].mapid1;
	  register OFFSET_T offset = pBuf[i].offset + delta;
	  register short score = MatchScoreFind(mapid1,offset, RANGE>=1 ? pBuf[i].MinOffset1[0]:0, HASH_SCALE>=1 ? pBuf[i].scaleID[0]:0);

	  if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf[i].offset - TraceOffset) <= OffsetSpread)
	     && (!OFFSET_ZEROFIX || score > 0)))
	    printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,sc=%d,delta=%d:score=%d to %d from offset=%d\n",
		   mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,RANGE>=1 ? pBuf[i].MinOffset1[0]:0,
		   (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,pBuf[i].score,pBuf[i].score+score, offset);

	  if(!OFFSET_ZEROFIX || score > 0)
	    pBuf[i].score += score * wt;
	}

	for(; i < MatchCntNull; i++){/* 2nd post loop */
	  register int mapid1 = pBuf[i].mapid1;
	  register OFFSET_T offset = pBuf[i].offset + delta;
	  register short score = MatchScoreFind(mapid1,offset, RANGE>=1 ? pBuf[i].MinOffset1[0]:0, HASH_SCALE>=1 ? pBuf[i].scaleID[0]:0);

	  if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf[i].offset - TraceOffset) <= OffsetSpread)
	     && (!OFFSET_ZEROFIX || score > 0)))
	    printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,sc=%d,delta=%d:score=%d to %d from offset=%d\n",
		   mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset, RANGE>=1 ? pBuf[i].MinOffset1[0]:0,
		   (HASH_SCALE && hashScaleDelta) ? pBuf[i].scaleID[0] : 0, delta,pBuf[i].score,pBuf[i].score + score * wt, offset);

	  if(!OFFSET_ZEROFIX || score > 0)
	    pBuf[i].score += score * wt;
	}
      }
    }
    nMatchCnt = MatchCntNull;
  }

  if(VERB>=2 || (TRACE && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2)){
    #pragma omp critical
    {
      printf("MatchFind(mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d):tid=%d: Updated MatchBuf[] with shifted copies of MatchTable : max_spread=%d: wall time=%0.6f\n",
	     mapid2,flip2,minSite2,maxSite2,tid,max_spread,wtime());
      fflush(stdout);
    }
  }

  /* for consistency of results with different mapsets for the same rmap[mapid1]->id values, sort MatchBuf[0..MatchCnt-1] by offset, so identical scores are resolved in favor of smallests offset */
  if(DETERMINISTIC){
    qsort(MatchBuf,(size_t)nMatchCnt,sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), (intcmp*)CBestMatchOffsetInc<OFFSET_T,RANGE,HASH_SCALE>);// NOTE : not needed except for deterministic results when the inserted mapsets change

    if(VERB>=2){
#pragma omp critical
      {
	printf("MatchFind(mapid2=%d):tid=%d: sorted MatchBuf[]: wall time=%0.6f\n",mapid2,tid,wtime());
	fflush(stdout);
      }
    }
  }

  short scoreMult = (1 << HASHSCORE_SHIFT);
  short hashscore;
  if(RANGE <= 0)
    hashscore = HashScore * scoreMult;
  else { // RANGE >= 1
    hashscore = GroupedScore * scoreMult;

    if(ShiftOffset1 < 31){  /* copy back min(255, score >> HASHSCORE_SHIFT) from MatchBuf[0..orignMatchCnt-1] to MatchTable */
      if(VERB>=2 || (TRACE && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2)){
	printf("MatchFind:mapid2=%d(id=%lld),flip2=%d: copying back score / %d from MatchBuf[0..%lld] to MatchTable\n",mapid2,qmap[mapid2]->id,flip2,(1<<HASHSCORE_SHIFT),orignMatchCnt-1);
	fflush(stdout);
      }
      
      pBuf = MatchBuf;
      // HERE : try using small circular buffers for mapid1, offset+delta, nindex to avoid have to recompute them more than once
      _mm_prefetch(CAST pBuf, _MM_HINT_T0);
      for(i = 0; i < PREFETCH_D && i < orignMatchCnt; i++){/* 1st prefetch loop */
	_mm_prefetch(CAST &Maphash2[pBuf[i].mapid1], _MM_HINT_T0);
	_mm_prefetch(CAST &Offsethash[pBuf[i].offset], _MM_HINT_T0);
	_mm_prefetch(CAST &Offset1hash[pBuf[i].MinOffset1[0]], _MM_HINT_T0);
	if(HASH_SCALE >= 1)
	  _mm_prefetch(CAST &scaleIDhash[pBuf[i].scaleID[0]], _MM_HINT_T0);
      }

      for(i = 0; i < PREFETCH_D && i < orignMatchCnt - PREFETCH_D; i++){/* 2nd prefetch loop */
	_mm_prefetch(CAST &Maphash2[pBuf[i+PREFETCH_D].mapid1], _MM_HINT_T0);
	_mm_prefetch(CAST &Offsethash[pBuf[i+PREFETCH_D].offset], _MM_HINT_T0);
	_mm_prefetch(CAST &Offset1hash[pBuf[i+PREFETCH_D].MinOffset1[0]], _MM_HINT_T0);
	if(HASH_SCALE >= 1)
	  _mm_prefetch(CAST &scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]], _MM_HINT_T0);

	register unsigned int nindex = Maphash2[pBuf[i].mapid1] ^ Offsethash[pBuf[i].offset] ^ Offset1hash[pBuf[i].MinOffset1[0]];
	if(HASH_SCALE >= 1)
	  nindex ^= scaleIDhash[pBuf[i].scaleID[0]];
	_mm_prefetch(CAST &MatchIndex[nindex], _MM_HINT_T0);
      }

      for(i = 0; i < orignMatchCnt - PREFETCH_D2; i++){/* main loop */
	if(PREFETCH_D){
	  _mm_prefetch(CAST &Maphash2[pBuf[i+PREFETCH_D2].mapid1], _MM_HINT_T0);
	  _mm_prefetch(CAST &Offsethash[pBuf[i+PREFETCH_D2].offset], _MM_HINT_T0);
	  _mm_prefetch(CAST &Offset1hash[pBuf[i+PREFETCH_D2].MinOffset1[0]], _MM_HINT_T0);
	  if(HASH_SCALE >= 1)
	    _mm_prefetch(CAST &scaleIDhash[pBuf[i+PREFETCH_D2].scaleID[0]], _MM_HINT_T0);

	  unsigned int nindex = Maphash2[pBuf[i+PREFETCH_D].mapid1] ^ Offsethash[pBuf[i+PREFETCH_D].offset] ^ Offset1hash[pBuf[i+PREFETCH_D].MinOffset1[0]];
	  if(HASH_SCALE >= 1)
	    nindex ^= scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]];
	  _mm_prefetch(CAST &MatchIndex[nindex], _MM_HINT_T0);
	}

	register int mapid1 = pBuf[i].mapid1;
	register OFFSET_T offset = pBuf[i].offset;
	unsigned int nscore = pBuf[i].score >> HASHSCORE_SHIFT;

	MatchScoreUpdate(mapid1, offset, pBuf[i].MinOffset1[0], pBuf[i].scaleID[0], min(255,nscore),mapid2,flip2);
      }

      for(; i < orignMatchCnt - PREFETCH_D; i++){/* 1st post loop */
	unsigned int nindex = Maphash2[pBuf[i+PREFETCH_D].mapid1] ^ Offsethash[pBuf[i+PREFETCH_D].offset] ^ Offset1hash[pBuf[i+PREFETCH_D].MinOffset1[0]];
	if(HASH_SCALE >= 1)
	  nindex ^= scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]];
	_mm_prefetch(CAST &MatchIndex[nindex], _MM_HINT_T0);

	register int mapid1 = pBuf[i].mapid1;
	register OFFSET_T offset = pBuf[i].offset;
	unsigned int nscore = pBuf[i].score >> HASHSCORE_SHIFT;
	MatchScoreUpdate(mapid1,offset,pBuf[i].MinOffset1[0], pBuf[i].scaleID[0], min(255,nscore),mapid2,flip2);
      }

      for(; i < orignMatchCnt; i++){/* 2nd post loop */
	register int mapid1 = pBuf[i].mapid1;
	register OFFSET_T offset = pBuf[i].offset;
	unsigned int nscore = pBuf[i].score >> HASHSCORE_SHIFT;
	MatchScoreUpdate(mapid1,offset,pBuf[i].MinOffset1[0], pBuf[i].scaleID[0], min(255,nscore),mapid2,flip2);
      }

      /* check neighboring entries (Offset1+1 or Offset1-1) and discard entries if scores do not add up to GroupedScore threshold in combination with either neighboring entry */
      _mm_prefetch(CAST pBuf, _MM_HINT_T0);
      for(i = 0; i < PREFETCH_D && i < nMatchCnt; i++){/* 1st prefetch loop */
	_mm_prefetch(CAST &Maphash2[pBuf[i].mapid1], _MM_HINT_T0);
	_mm_prefetch(CAST &Offsethash[pBuf[i].offset], _MM_HINT_T0);
	unsigned short Offset1 = pBuf[i].MinOffset1[0];
	if(DEBUG>=2) assert(0 < Offset1 && Offset1 < MASK(16));
	_mm_prefetch(CAST &Offset1hash[Offset1 - 1u], _MM_HINT_T0);
	_mm_prefetch(CAST &Offset1hash[Offset1 + 1u], _MM_HINT_T0);
	if(HASH_SCALE >= 1)
	  _mm_prefetch(CAST &scaleIDhash[pBuf[i].scaleID[0]], _MM_HINT_T0);
      }

      for(i = 0; i < PREFETCH_D && i < nMatchCnt - PREFETCH_D; i++){/* 2nd prefetch loop */
	_mm_prefetch(CAST &Maphash2[pBuf[i+PREFETCH_D].mapid1], _MM_HINT_T0);
	_mm_prefetch(CAST &Offsethash[pBuf[i+PREFETCH_D].offset], _MM_HINT_T0);
	unsigned short Offset1 = pBuf[i+PREFETCH_D].MinOffset1[0];
	if(DEBUG>=2) assert(0 < Offset1 && Offset1 < MASK(16));
	_mm_prefetch(CAST &Offset1hash[Offset1 - 1u], _MM_HINT_T0);
	_mm_prefetch(CAST &Offset1hash[Offset1 + 1u], _MM_HINT_T0);
	if(HASH_SCALE >= 1)
	  _mm_prefetch(CAST &scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]], _MM_HINT_T0);

	register unsigned int nindex = Maphash2[pBuf[i].mapid1] ^ Offsethash[pBuf[i].offset];
	if(HASH_SCALE >= 1)
	  nindex ^= scaleIDhash[pBuf[i].scaleID[0]];
	Offset1 = pBuf[i].MinOffset1[0];
	if(DEBUG>=2) assert(0 < Offset1 && Offset1 < MASK(16));
	_mm_prefetch(CAST &MatchIndex[nindex ^ Offset1hash[Offset1 - 1u]], _MM_HINT_T0);
	_mm_prefetch(CAST &MatchIndex[nindex ^ Offset1hash[Offset1 + 1u]], _MM_HINT_T0);
      }

      for(i = 0; i < nMatchCnt - PREFETCH_D2; i++){/* main loop */
	if(PREFETCH_D){
	  _mm_prefetch(CAST &Maphash2[pBuf[i+PREFETCH_D2].mapid1], _MM_HINT_T0);
	  _mm_prefetch(CAST &Offsethash[pBuf[i+PREFETCH_D2].offset], _MM_HINT_T0);
	  unsigned short Offset1 = pBuf[i+PREFETCH_D2].MinOffset1[0];
	  if(DEBUG>=2) assert(0 < Offset1 && Offset1 < MASK(16));
	  _mm_prefetch(CAST &Offset1hash[Offset1 - 1u], _MM_HINT_T0);
	  _mm_prefetch(CAST &Offset1hash[Offset1 + 1u], _MM_HINT_T0);
	  if(HASH_SCALE >= 1)
	    _mm_prefetch(CAST &scaleIDhash[pBuf[i+PREFETCH_D2].scaleID[0]], _MM_HINT_T0);

	  unsigned int nindex = Maphash2[pBuf[i+PREFETCH_D].mapid1] ^ Offsethash[pBuf[i+PREFETCH_D].offset];
	  if(HASH_SCALE >= 1)
	    nindex ^= scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]];
	  Offset1 = pBuf[i+PREFETCH_D].MinOffset1[0];
	  if(DEBUG>=2) assert(0 < Offset1 && Offset1 < MASK(16));
	  _mm_prefetch(CAST &MatchIndex[nindex ^ Offset1hash[Offset1 - 1u]], _MM_HINT_T0);
	  _mm_prefetch(CAST &MatchIndex[nindex ^ Offset1hash[Offset1 + 1u]], _MM_HINT_T0);
	}

	register short score = pBuf[i].score;
	if(score >= hashscore)
	  continue;
	register int mapid1 = pBuf[i].mapid1;
	register OFFSET_T offset = pBuf[i].offset;
	unsigned short Offset1 = pBuf[i].MinOffset1[0];
	short scoreA = 0;
	if(Offset1){
	  scoreA = MatchScoreFind(mapid1,offset,Offset1 - 1u, pBuf[i].scaleID[0]);
	  if(score + (scoreA << HASHSCORE_SHIFT) >= hashscore)
	    continue;
	}
	short scoreB = MatchScoreFind(mapid1,offset,Offset1 + 1u, pBuf[i].scaleID[0]);
	if(score + (scoreB << HASHSCORE_SHIFT) >= hashscore)
	  continue;

	if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf[i].offset - TraceOffset) <= OffsetSpread) 
		       && pBuf[i].score > 0)){
	  printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1= %u,sc=%d,score=%d to 0 due to GroupedScore=%d (hashscore=%d,scoreA=%d,scoreB=%d)\n",
		 mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,pBuf[i].MinOffset1[0],(HASH_SCALE>=1)?pBuf[i].scaleID[0]:0, pBuf[i].score, GroupedScore, hashscore,scoreA,scoreB);
	  fflush(stdout);
	}

	pBuf[i].score = 0;    /* Invalidate pBuf[i] */
      }

      for(; i < nMatchCnt - PREFETCH_D; i++){/* 1st post loop */
	unsigned int nindex = Maphash2[pBuf[i+PREFETCH_D].mapid1] ^ Offsethash[pBuf[i+PREFETCH_D].offset];
	if(HASH_SCALE >= 1)
	  nindex ^= scaleIDhash[pBuf[i+PREFETCH_D].scaleID[0]];
	unsigned short Offset1 = pBuf[i+PREFETCH_D].MinOffset1[0];
	if(DEBUG>=2) assert(0 < Offset1 && Offset1 < MASK(16));
	_mm_prefetch(CAST &MatchIndex[nindex ^ Offset1hash[Offset1 - 1u]], _MM_HINT_T0);
	_mm_prefetch(CAST &MatchIndex[nindex ^ Offset1hash[Offset1 + 1u]], _MM_HINT_T0);

	register short score = pBuf[i].score;
	if(score >= hashscore)
	  continue;
	register int mapid1 = pBuf[i].mapid1;
	register OFFSET_T offset = pBuf[i].offset;
	Offset1 = pBuf[i].MinOffset1[0];
	short scoreA = 0;
	if(Offset1){
	  scoreA = MatchScoreFind(mapid1,offset,Offset1 - 1u, pBuf[i].scaleID[0]);
	  if(score + (scoreA << HASHSCORE_SHIFT) >= hashscore)
	    continue;
	}
	short scoreB = MatchScoreFind(mapid1,offset,Offset1 + 1u, pBuf[i].scaleID[0]);
	if(score + (scoreB << HASHSCORE_SHIFT) >= hashscore)
	  continue;

	if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf[i].offset - TraceOffset) <= OffsetSpread)
		       && pBuf[i].score > 0)){
	  printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1= %u,sc=%d,score=%d to 0 due to GroupedScore=%d (hashscore=%d,scoreA=%d,scoreB=%d)\n",
		 mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,pBuf[i].MinOffset1[0],HASH_SCALE>=1 ?pBuf[i].scaleID[0]:0,pBuf[i].score, GroupedScore, hashscore,scoreA,scoreB);
	  fflush(stdout);
	}

	pBuf[i].score = 0;    /* Invalidate pBuf[i] */
      }
      for(; i < nMatchCnt; i++){/* 2nd post loop */
	register short score = pBuf[i].score;
	if(score >= hashscore)
	  continue;
	register int mapid1 = pBuf[i].mapid1;
	register OFFSET_T offset = pBuf[i].offset;
	unsigned short Offset1 = pBuf[i].MinOffset1[0];
	short scoreA = 0;
	if(Offset1){
	  scoreA = MatchScoreFind(mapid1,offset,Offset1 - 1u,pBuf[i].scaleID[0]);
	  if(score + (scoreA << HASHSCORE_SHIFT) >= hashscore)
	    continue;
	}
	short scoreB = MatchScoreFind(mapid1,offset,Offset1 + 1u,pBuf[i].scaleID[0]);
	if(score + (scoreB << HASHSCORE_SHIFT) >= hashscore)
	  continue;

	if(VERB>=2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(pBuf[i].offset - TraceOffset) <= OffsetSpread) 
		       && pBuf[i].score > 0)){
	  printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1= %u,sc=%d,score=%d to 0 due to GroupedScore=%d (hashscore=%d,scoreA=%d,scoreB=%d)\n",
		 mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf[i].offset,pBuf[i].MinOffset1[0],HASH_SCALE>=1?pBuf[i].scaleID[0]:0,pBuf[i].score, GroupedScore, hashscore,scoreA,scoreB);
	  fflush(stdout);
	}
	pBuf[i].score = 0;    /* Invalidate pBuf[i] */
      }

      // any remaining +ve score will be combined with other scores with same offset : the sum of at least 2 of them must be >= GroupedScore * scoreMult (see above)
      hashscore = 1;
    }// if(ShiftOffset1 < 31)
  } // RANGE >= 1

  if(VERB>=2 || (TRACE && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2)){
    printf("BestMatch:mapid2=%d(id=%lld),flip2=%d:nMatchCnt=%lld\n",mapid2,qmap[mapid2]->id,flip2,nMatchCnt);
    fflush(stdout);
  }

  /* Copy MatchBuf[0..nMatchCnt-1] to BestMatch[] */
  if(HashMultiMatch){/* (Re)allocate BestMatch[mapid][] and populate BestMapid[0..BestCnt-1] */
    if(DEBUG>=2 && minSite2 <= 0){
      for(int i = 0; i < numrmaps; i++)
	assert(BestMatchMax[roffset+i] == 0);
      assert(BestMatchMem == NULL);
    }

    pBuf = MatchBuf;
    _mm_prefetch(CAST pBuf,_MM_HINT_T0);
    _mm_prefetch(CAST &pBuf[8], _MM_HINT_T0);
    for(i = 0; i < PREFETCH_D && i < nMatchCnt; i++)  /* prefetch first PREFETCH_D BestMatch Entries */
      _mm_prefetch(CAST &BestMatch[pBuf[i].mapid1], _MM_HINT_T0);
    
    for(i = 0; i < nMatchCnt - PREFETCH_D; i++, pBuf++){
      if(PREFETCH_D) _mm_prefetch(CAST &BestMatch[pBuf[PREFETCH_D].mapid1], _MM_HINT_T0);
      if(pBuf->score >= hashscore)
	if(!BestMatchMax[pBuf->mapid1]++)
	  BestMapid[BestCnt++] = pBuf->mapid1;
    }
    for(; i < nMatchCnt; i++, pBuf++)
      if(pBuf->score >= hashscore)
	if(!BestMatchMax[pBuf->mapid1]++)
	  BestMapid[BestCnt++] = pBuf->mapid1;
    int totcnt = 0;
    for(int i = 0; i < BestCnt; i++)
      totcnt += BestMatchMax[BestMapid[i]];
    if(DEBUG) assert(totcnt >= 0);/* check for overflow of totcnt */

    int next = 0;
    if(!HashGC || minSite2 <= 0){/* original code : no re-allocation required */
      BestMatchMem = new CBestMatch<OFFSET_T,RANGE,HASH_SCALE>[totcnt];
      for(int i = 0; i < BestCnt; i++){
	int mapid1 = BestMapid[i];
	BestMatch[mapid1] = &BestMatchMem[next];
	next += BestMatchMax[mapid1];
      }
      if(VERB>=2){
        #pragma omp critical
	{
	  printf("MatchFind(mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d):tid=%d: Allocated %lu bytes for BestMatchMem[], totcnt= %d, next= %d, nMatchCnt=%lld\n",
		 mapid2,flip2,minSite2,maxSite2,tid,(size_t)totcnt * sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), totcnt,next,nMatchCnt);
	  fflush(stdout);
	}
      }
    } else {/* reallocate BestMatch[mapid1][] and copy MatchBuf[] to BestMatch[] && BestMapid[0..BestCnt-1] */
      CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *nBestMatchMem = new CBestMatch<OFFSET_T,RANGE,HASH_SCALE>[totcnt];
      int origtotcnt = 0;

      if(PREFETCH_BM){
	_mm_prefetch(CAST BestMapid, _MM_HINT_T0);// sequential access
	for(int i = 0; i < PREFETCH_BM && i < BestCnt; i++){
	  int mapid1 = BestMapid[i];
	  _mm_prefetch(CAST &BestMatchCnt[mapid1], _MM_HINT_T0);
	  _mm_prefetch(CAST &BestMatchMax[mapid1], _MM_HINT_T0);
	  _mm_prefetch(CAST &BestMatch[mapid1], _MM_HINT_T0);
	}
      }

      for(int i = 0; i < BestCnt - PREFETCH_BM; i++){
	if(PREFETCH_BM){// prefetch (BestMatchCnt,BestMatch,BestMatchMax)[BestMapid[i+1 .. i+PREFETCH_MB]]
	  int mapid1 = BestMapid[i + PREFETCH_BM];
	  _mm_prefetch(CAST &BestMatchCnt[mapid1], _MM_HINT_T0);
	  _mm_prefetch(CAST &BestMatchMax[mapid1], _MM_HINT_T0);
	  _mm_prefetch(CAST &BestMatch[mapid1], _MM_HINT_T0);
	}

	int mapid1 = BestMapid[i];
	CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *pBestMatch = &nBestMatchMem[next];
	if(BestMatchCnt[mapid1] > 0){
	  origtotcnt += BestMatchCnt[mapid1];
	  if(DEBUG>=2) assert(BestMatchCnt[mapid1] <= BestMatchMax[mapid1]);
	  memcpy(pBestMatch, BestMatch[mapid1], BestMatchCnt[mapid1] * sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>));
	}
	BestMatch[mapid1] = pBestMatch;
	next += BestMatchMax[mapid1];
      }

      for(int i = max(0, BestCnt - PREFETCH_BM); i < BestCnt; i++){
	int mapid1 = BestMapid[i];
	CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *pBestMatch = &nBestMatchMem[next];
	if(BestMatchCnt[mapid1] > 0){
	  origtotcnt += BestMatchCnt[mapid1];
	  if(DEBUG>=2) assert(BestMatchCnt[mapid1] <= BestMatchMax[mapid1]);
	  memcpy(pBestMatch, BestMatch[mapid1], BestMatchCnt[mapid1] * sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>));
	}
	BestMatch[mapid1] = pBestMatch;
	next += BestMatchMax[mapid1];
      }

      delete [] BestMatchMem;
      BestMatchMem = nBestMatchMem;

      if(VERB>=2){
        #pragma omp critical
	{
	  printf("MatchFind(mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d):tid=%d: (Re)Allocated %lu bytes for BestMatchMem[], totcnt= %d -> %d, next= %d, nMatchCnt=%lld\n",
		 mapid2,flip2,minSite2,maxSite2,tid,(size_t)totcnt * sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), origtotcnt, totcnt,next,nMatchCnt);
	  fflush(stdout);
	}
      }
    }
    if(DEBUG) assert(next == totcnt);
  } // else BestMatch[mapid][0..0] has already been allocated 

  pBuf = MatchBuf;
  _mm_prefetch(CAST pBuf,_MM_HINT_T0);
  _mm_prefetch(CAST &pBuf[8], _MM_HINT_T0);
  for(i = 0; i < PREFETCH_D && i < nMatchCnt; i++)  /* prefetch first PREFETCH_D BestMatch Entries */
    _mm_prefetch(CAST &BestMatch[pBuf[i].mapid1], _MM_HINT_T0);

  if(!HashGC || maxSite2 <= 0){/* copy all entries in MatchBuf[] to BestMatch[] */
    for(i = 0; i < nMatchCnt - PREFETCH_D; i++, pBuf++){
      if(PREFETCH_D) _mm_prefetch(CAST &BestMatch[pBuf[PREFETCH_D].mapid1], _MM_HINT_T0);
      if(pBuf->score >= hashscore){
	BestInsert(pBuf->mapid1,pBuf->offset,pBuf->score,pBuf,mapid2,flip2,tid);
      } else if((VERB >= 2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && 
			       (TraceOffset == INT32_MIN || abs(pBuf->offset - TraceOffset) <= OffsetSpread))) && (!RANGE || pBuf->score > 0)){
	printf("\t Skipping BestInsert:tid=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,sc=%d,score=%d: hashscore=%d\n", 
	       tid,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,
	       (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, pBuf->score,hashscore);
	fflush(stdout);
      }
    }
    for(; i < nMatchCnt; i++, pBuf++){
      if(pBuf->score >= hashscore){
	BestInsert(pBuf->mapid1,pBuf->offset,pBuf->score,pBuf,mapid2,flip2,tid);
      } else if((VERB >= 2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && 
			       (TraceOffset == INT32_MIN || abs(pBuf->offset - TraceOffset) <= OffsetSpread))) && (!RANGE || pBuf->score > 0)){
	printf("\t Skipping BestInsert:tid=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,sc=%d,score=%d: hashscore=%d\n", 
	       tid,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,
	       (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, pBuf->score,hashscore);
	fflush(stdout);
      }
    }
  } else {/* copy only those entries in MatchBuf[] with (minOffset2 - Offset1flip2[mapid1] <= offset <= maxOffset2 - Offset1flip2[mapid1]) to BestMatch[] */
    if(VERB>=2 && mapid2==MAPID2 && flip2==FLIP2){
      printf("tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d:nMatchCnt=%lld: copying matches with %d <= offset + Offset1(mapid) <= %d to BestMatch[]\n",
	     tid,mapid2,flip2,minSite2,maxSite2,nMatchCnt,minOffset2,maxOffset2);
      fflush(stdout);
    }

    for(i = 0; i < PREFETCH_D && i < nMatchCnt; i++)  /* prefetch first PREFETCH_D Offset1flip2[mapid] Entries */
      _mm_prefetch(CAST &Offset1flip2[pBuf[i].mapid1], _MM_HINT_T0);

    for(i = 0; i < nMatchCnt - PREFETCH_D; i++, pBuf++){/* main loop */
      if(PREFETCH_D){
        int mapid1 = pBuf[PREFETCH_D].mapid1;
	_mm_prefetch(CAST &BestMatch[mapid1], _MM_HINT_T0);
        _mm_prefetch(CAST &Offset1flip2[mapid1], _MM_HINT_T0);
      }
      if(pBuf->score >= hashscore){
	int mapid1 = pBuf->mapid1;
	OFFSET_T offset = pBuf->offset;
	OFFSET_T Offset2 = offset + Offset1flip2[mapid1];
	if(minOffset2 <= Offset2 && Offset2 <= maxOffset2){
	  BestInsert(mapid1,offset,pBuf->score,pBuf,mapid2,flip2,tid);
        } else if((VERB >= 2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && 
				 (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))) && (!RANGE || pBuf->score > 0)){
	  printf("\t Skipping BestInsert:tid=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,sc=%d,score=%d: minOffset2=%d,maxOffset2=%d,Offset2=%d,Offset1flip2[mapid1]=%d\n", 
		 tid,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,
		 (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, pBuf->score,minOffset2,maxOffset2,Offset2,Offset1flip2[mapid1]);
	  fflush(stdout);
	}
      } else if((VERB >= 2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2)) && (!RANGE || pBuf->score > 0)){
	printf("\t Skipping BestInsert:tid=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,sc=%d,score=%d: hashscore=%d\n", 
	       tid,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,
	       (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, pBuf->score,hashscore);
	fflush(stdout);
      }
    }
    for(; i < nMatchCnt; i++, pBuf++){/* post loop */
      if(pBuf->score >= hashscore){
	int mapid1 = pBuf->mapid1;
	OFFSET_T offset = pBuf->offset;
	OFFSET_T Offset2 = offset + Offset1flip2[mapid1];
	if(minOffset2 <= Offset2 && Offset2 <= maxOffset2)
	  BestInsert(mapid1,offset,pBuf->score,pBuf,mapid2,flip2,tid);
	else if(VERB >= 2 || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2)){
	  printf("\t Skipping BestInsert:tid=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,sc=%d,score=%d: minOffset2=%d,maxOffset2=%d,Offset2=%d,Offset1flip2[mapid1]=%d\n", 
		 tid,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,
		 (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, pBuf->score,minOffset2,maxOffset2,Offset2,Offset1flip2[mapid1]);
	  fflush(stdout);
	}
      } else if((VERB >= 2 || (TRACE && rmap[pBuf->mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && 
			       (TraceOffset == INT32_MIN || abs(pBuf->offset - TraceOffset) <= OffsetSpread))) && (!RANGE || pBuf->score > 0)){
	printf("\t Skipping BestInsert:tid=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,sc=%d,score=%d: hashscore=%d\n",
	       tid,pBuf->mapid1,rmap[pBuf->mapid1]->id,mapid2,qmap[mapid2]->id,flip2,pBuf->offset,
	       (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, pBuf->score,hashscore);
	fflush(stdout);
      }
    }
  }

  if(!(RANGE && ShiftOffset1 <= 30 && HashGC && maxSite2 > 0)){ /* reduce memory used in MatchBuf[] array */
    if(DEBUG && MatchCnt > 0 && !(MatchBufMax > 0)){
      printf("MatchBufMax=%lld, MatchCnt=%lld\n",MatchBufMax,MatchCnt);
      fflush(stdout);
      assert(MatchBufMax > 0);
    }
    if(MatchBufMax > 0){
      // WAS delete [] MatchBuf;
      size_t memsiz = (sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>) * MatchBufMax + PAGE-1) & ~(PAGE-1);
      if(munmap(MatchBuf, memsiz)){
	int eno = errno;
	char *err = strerror(eno);

        #pragma omp critical
	{
	  printf("munmap(%p,%lu) failed: errno=%d:%s\n",MatchBuf, memsiz, eno, err);
	  dumpmemmap();
	  fflush(stdout);exit(1);
	}
      }

      if(hashtable->MatchMemReserve)
	hashtable->MatchMemReserve->ResFree(MatchBufMax * sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), 1, this);
    }
    MatchBufMax = 0;
    MatchBuf = 0;
  }

  if(HashGC && maxSite2 > 0){/* compress memory for MatchOverflow[] & MatchIndex[] */
    OFFSET_T Offset2 = (flip2==0) ? maxOffset2 - max_spread + 1 : minOffset2 + max_spread - 1;
    /* Remove all entries from MatchIndex[] & MatchOverflow[] with (flip ? offset > Offset2 - Offset1flip2[mapid1] : offset < Offset2 - Offset1flip2[mapid1])
       If RANGE : restore score values for remaining entries from MatchBuf : (minOffset2 - Offset1flip2[mapid1] <= offset <= maxOffset2 - Offset1flip2[mapid1]) */

    if(VERB>=2 && mapid2==MAPID2 && flip2==FLIP2){
      printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d: MatchFirst=%d, MatchOverflowNext=%d,MatchOverflowFree=%d, Offset2= %d\n",
	     tid,mapid2,flip2,minSite2,maxSite2,MatchFirst,MatchOverflowNext,MatchOverflowFree, Offset2);
      fflush(stdout);
    }

    int *qnext = &MatchFirst, *pnext;
    register int qindex = *qnext;
    register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *q = &MatchIndex[qindex], *p;
    _mm_prefetch(CAST q, _MM_HINT_T0);

    register long long totcnt = 0;

    for(register int pindex = qindex; pindex >= 0; pindex = qindex){
      pnext = qnext;
      p = q;
      qnext = &p->next;
      q = &MatchIndex[qindex = *qnext];
      _mm_prefetch(CAST q, _MM_HINT_T0);    

      long long origtotcnt;
      unsigned int origkeepcnt,keepcnt, deletecnt;
      if(DEBUG>=2){
	origtotcnt = totcnt;
	origkeepcnt = keepcnt = deletecnt = 0;
      }

      register unsigned int cnt = p->cnt;
      if(DEBUG>=2) assert(cnt > 0);

      if(cnt > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){/* handle overflow linked list */
	int origFreeCnt, OverflowCnt;
	if(DEBUG>=2){/* compute Overflow List Length and original FreeList length */
	  origkeepcnt += MATCHLINE(OFFSET_T,RANGE,HASH_SCALE);
	  OverflowCnt = 0;
	  origFreeCnt = 0;
          #pragma omp critical
	  {
	    if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2 && (MAXSITE2== -1 || maxSite2==MAXSITE2)){
	      printf("tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,Free List starting at MatchOverflowFree=%d:\n",
	        tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,MatchOverflowFree);
	      fflush(stdout);
	    }
	    for(int i = MatchOverflowFree; i > 0; i = MatchOverflow[i].next){
	      origFreeCnt++;
	      if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2 && (MAXSITE2== -1 || maxSite2==MAXSITE2)){
	        printf("\t tid=%d: i=%d, MatchOverflow[i].next=%d,cnt=%d: origFreeCnt=%d\n", tid,i,MatchOverflow[i].next,MatchOverflow[i].cnt,origFreeCnt);
	        fflush(stdout);
	      }
	      if(DEBUG) assert(i > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) && i <= MatchOverflowNext);
	      if(DEBUG) assert(MatchOverflow[i].cnt == 0);
	    }
	    int Cnt = 0;
	    if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2 && (MAXSITE2== -1 || maxSite2==MAXSITE2))
	      printf("tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,Overflow linked listed starting at p->cnt= %d:\n",
	         tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->cnt);
	    for(int i = cnt; i > 0; i = MatchOverflow[i].next){
	      Cnt++;
	      if((VERB>=3 && mapid2==MAPID2 && flip2==FLIP2 && (MAXSITE2== -1 || maxSite2==MAXSITE2)) || (DEBUG && !(MatchOverflow[i].cnt > 0))){
	        printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d:i=%d, MatchOverflow[i].next=%d,cnt=%d: Cnt= %d\n", 
		       tid, mapid2,flip2,minSite2,maxSite2,i,MatchOverflow[i].next,MatchOverflow[i].cnt,Cnt);
		fflush(stdout);
	      }
	      if(DEBUG) assert(MatchOverflow[i].cnt > 0);
	      origkeepcnt += MatchOverflow[i].cnt;
	    }
	    fflush(stdout);
	  }
	}

	register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *s,*r = &MatchOverflow[cnt];
	_mm_prefetch(CAST r, _MM_HINT_T0);    

        register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *np = p, *nq = r;
        register unsigned int ncnt = 0;

	if(flip2 == 0){
	  for(register int i = 0; i < MATCHLINE(OFFSET_T,RANGE,HASH_SCALE); i++){
	    int mapid1 = p->mapid1[i];
	    OFFSET_T offset = p->offset[i];
	    if(offset + Offset1flip2[mapid1] < Offset2){
	      if(DEBUG>=2)	deletecnt++;
	      if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
		printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->cnt=%d,i=%d:Deleting mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d\n",
		       tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->cnt,i,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next);
		fflush(stdout);
	      }	  

	      continue;
	    }
	  
	    np->mapid1[ncnt] = mapid1;
	    np->offset[ncnt] = offset;
	    np->score[ncnt] = p->score[i];
	    if(RANGE >= 1)
	      np->Offset1[ncnt] = p->Offset1[i];
	    if(HASH_SCALE >= 1)
	      np->scaleID[ncnt] = p->scaleID[i];
	    ncnt++;

	    if(DEBUG>=3) keepcnt++;
	    if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	      printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->cnt=%d,i=%d:Keeping mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d\n",
		     tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->cnt,i,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next);
	      fflush(stdout);
	    }
	  }
	} else {// flip==1
	  for(register int i = 0; i < MATCHLINE(OFFSET_T,RANGE,HASH_SCALE); i++){
	    int mapid1 = p->mapid1[i];
	    OFFSET_T offset = p->offset[i];
	    if(offset + Offset1flip2[mapid1] > Offset2){
	      if(DEBUG>=2)	deletecnt++;
	      if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
		printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->cnt=%d,i=%d:Deleting mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d\n",
		       tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->cnt,i,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next);
		fflush(stdout);
	      }	  

	      continue;
	    }
	  
	    np->mapid1[ncnt] = mapid1;
	    np->offset[ncnt] = offset;
	    np->score[ncnt] = p->score[i];
	    if(RANGE >= 1)
	      np->Offset1[ncnt] = p->Offset1[i];
	    if(HASH_SCALE >= 1)
	      np->scaleID[ncnt] = p->scaleID[i];

	    ncnt++;
	    if(DEBUG>=3) keepcnt++;
	    if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	      printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->cnt=%d,i=%d:Keeping mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d\n",
		     tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->cnt,i,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next);
	      fflush(stdout);
	    }
	  }
	}

	for(; r->next >= 0; r = s){
	  if(DEBUG>=2) OverflowCnt++;
          if(DEBUG>=2) assert(r->next > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	  s = &MatchOverflow[r->next];
	  _mm_prefetch(CAST s, _MM_HINT_T0);

	  register unsigned int rcnt = r->cnt;
	  if(DEBUG>=2) assert(rcnt > 0 && rcnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));

	  if(flip2 == 0){
	    for(register unsigned int i = 0; i < rcnt; i++){
	      int mapid1 = r->mapid1[i];
	      OFFSET_T offset = r->offset[i];
	      if(offset + Offset1flip2[mapid1] < Offset2){
		if(DEBUG>=2) deletecnt++;
		if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
		  printf("  tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,r->next=%d,i=%d/%d:Deleting mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d,OverCnt=%d\n",
			 tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,r->next,i,rcnt,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next,OverflowCnt);
		  fflush(stdout);
		}
		continue;
	      }
	    
	      if(ncnt >= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){
		totcnt += MATCHLINE(OFFSET_T,RANGE,HASH_SCALE);
		if(DEBUG>=2) assert(nq != NULL);
		if(np != p)
		  np->cnt = ncnt;
		np = nq;
		if(DEBUG>=2) assert(nq->next > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
		nq = &MatchOverflow[nq->next];
		ncnt = 0;
		if(DEBUG>=2) OverflowCnt--;
	      }
	      np->mapid1[ncnt] = mapid1;
	      np->offset[ncnt] = offset;
	      np->score[ncnt] = r->score[i];  
	      if(RANGE >= 1)
		np->Offset1[ncnt] = r->Offset1[i];
	      if(HASH_SCALE >= 1)
		np->scaleID[ncnt] = r->scaleID[i];

	      ncnt++;
	      if(DEBUG>=3) keepcnt++;
	      if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
		printf("  tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,r->next=%d,i=%d/%d:Keeping mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d,OverCnt=%d\n",
		       tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,r->next,i,rcnt,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next,OverflowCnt);
		fflush(stdout);
	      }
	    }
	  } else {// flip2 == 1
	    for(register unsigned int i = 0; i < rcnt; i++){
	      int mapid1 = r->mapid1[i];
	      OFFSET_T offset = r->offset[i];
	      if(offset + Offset1flip2[mapid1] > Offset2){
		if(DEBUG>=2) deletecnt++;
		if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
		  printf("  tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,r->next=%d,i=%d/%d:Deleting mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d,OverCnt=%d\n",
			 tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,r->next,i,rcnt,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next,OverflowCnt);
		  fflush(stdout);
		}
		continue;
	      }
	    
	      if(ncnt >= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){
		totcnt += MATCHLINE(OFFSET_T,RANGE,HASH_SCALE);
		if(DEBUG>=2) assert(nq != NULL);
		if(np != p)
		  np->cnt = ncnt;
		np = nq;
		if(DEBUG>=2) assert(nq->next > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
		nq = &MatchOverflow[nq->next];
		ncnt = 0;
		if(DEBUG>=2) OverflowCnt--;
	      }
	      np->mapid1[ncnt] = mapid1;
	      np->offset[ncnt] = offset;
	      np->score[ncnt] = r->score[i];  
	      if(RANGE >= 1)
		np->Offset1[ncnt] = r->Offset1[i];
	      if(HASH_SCALE >= 1)
		np->scaleID[ncnt] = r->scaleID[i];
	      ncnt++;
	      if(DEBUG>=3) keepcnt++;
	      if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
		printf("  tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,r->next=%d,i=%d/%d:Keeping mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d,OverCnt=%d\n",
		       tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,r->next,i,rcnt,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next,OverflowCnt);
		fflush(stdout);
	      }
	    }
	  }
        }
	if(DEBUG>=2) assert(r->next < 0 && r->cnt > 0 && r->cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	if(DEBUG>=2) OverflowCnt++;

	if(flip2==0){
	  for(register unsigned int i = 0; i < r->cnt; i++){
	    int mapid1 = r->mapid1[i];
	    OFFSET_T offset = r->offset[i];
	    if(offset + Offset1flip2[mapid1] < Offset2){
	      if(DEBUG>=2) deletecnt++;
	      if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
		printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,r->next=%d,i=%d/%d:Deleting mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d\n",
		       tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,r->next,i,r->cnt,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next);
		fflush(stdout);
	      }
	      continue;
	    }

	    if(ncnt >= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){
	      totcnt += MATCHLINE(OFFSET_T,RANGE,HASH_SCALE);
	      if(DEBUG>=2) assert(nq != NULL);
	      if(np != p)
		np->cnt = ncnt;
	      np = nq;
	      nq = (nq->next < 0) ? NULL : &MatchOverflow[nq->next];
	      ncnt = 0;
	      if(DEBUG>=2) OverflowCnt--;
	    }
	    np->mapid1[ncnt] = mapid1;
	    np->offset[ncnt] = offset;
	    np->score[ncnt] = r->score[i];  
	    if(RANGE >= 1)
	      np->Offset1[ncnt] = r->Offset1[i];
	    if(HASH_SCALE >=1 )
	      np->scaleID[ncnt] = r->scaleID[i];

	    ncnt++;
	    if(DEBUG>=3) keepcnt++;
	    if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	      printf("  tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,r->next=%d,i=%d/%d:Keeping mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d,OverCnt=%d\n",
		     tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,r->next,i,r->cnt,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next,OverflowCnt);
	      fflush(stdout);
	    }
	  }
	} else {// flip2 == 1
	  for(register unsigned int i = 0; i < r->cnt; i++){
	    int mapid1 = r->mapid1[i];
	    OFFSET_T offset = r->offset[i];
	    if(offset + Offset1flip2[mapid1] > Offset2){
	      if(DEBUG>=2) deletecnt++;
	      if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
		printf("  tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,r->next=%d,i=%u/%u:Deleting mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%u,totcnt=%lld,np->next=%d,OverCnt=%d\n",
		       tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,r->next,i,r->cnt,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next,OverflowCnt);
		fflush(stdout);
	      }
	      continue;
	    }

	    if(ncnt >= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){
	      totcnt += MATCHLINE(OFFSET_T,RANGE,HASH_SCALE);
	      if(DEBUG>=2) assert(nq != NULL);
	      if(np != p)
		np->cnt = ncnt;
	      np = nq;
	      nq = (nq->next < 0) ? NULL : &MatchOverflow[nq->next];
	      ncnt = 0;
	      if(DEBUG>=2) OverflowCnt--;
	    }
	    np->mapid1[ncnt] = mapid1;
	    np->offset[ncnt] = offset;
	    np->score[ncnt] = r->score[i];  
	    if(RANGE >= 1)
	      np->Offset1[ncnt] = r->Offset1[i];
	    if(HASH_SCALE >= 1)
	      np->scaleID[ncnt] = r->scaleID[i];

	    ncnt++;
	    if(DEBUG>=3) keepcnt++;
	    if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	      printf("  tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,r->next=%d,i=%d/%d:Keeping mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld,np->next=%d,OverCnt=%d\n",
		     tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,r->next,i,r->cnt,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt,np->next,OverflowCnt);
	      fflush(stdout);
	    }
	  }
	}

	/* terminate overflow linked list (np,ncnt) and add remainder of linked list to Freelist */
	if(DEBUG>=2) assert(ncnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	if(DEBUG>=2) assert(p != r);
	if(np == r){
          r->cnt = ncnt;
	  totcnt += ncnt;
	  if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	    printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d:np==r:r->cnt -> %d, totcnt=%lld\n", 
		   tid,mapid2,flip2,minSite2,maxSite2,N2,r->cnt,totcnt);
	    fflush(stdout);
	  }
	  if(DEBUG>=2) assert(ncnt > 0);
	  if(DEBUG>=2) assert(totcnt-origtotcnt == keepcnt && origkeepcnt == keepcnt + deletecnt);
	  if(DEBUG>=2){/* check keepcnt */
	    unsigned int nkeepcnt = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE);
	    for(int i = p->cnt; i > 0; i = MatchOverflow[i].next){
	      if(DEBUG) assert(MatchOverflow[i].cnt > 0);
	      nkeepcnt += MatchOverflow[i].cnt;
	    }
	    if(DEBUG && !(nkeepcnt == keepcnt)){
	      #pragma omp critical
	      {
		printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d:np==r:r=[%ld],r->cnt -> %d, pindex=%d,p->cnt=%d,totcnt=%lld->%lld,keepcnt=%d->%d,deletecnt=%d,nkeepcnt=%d\n", 
		       tid,mapid2,flip2,minSite2,maxSite2,N2,r-MatchOverflow,r->cnt,pindex, p->cnt,origtotcnt,totcnt,origkeepcnt,keepcnt,deletecnt,nkeepcnt);
		fflush(stdout);
		nkeepcnt = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE);
		for(int i = p->cnt; i > 0; i = MatchOverflow[i].next){
		  nkeepcnt += MatchOverflow[i].cnt;
		  printf("\t\t tid=%d: i=%d, MatchOverflow[i].next=%d, MatchOverflow[i].cnt= %d, nkeepcnt=%d\n",tid,i,MatchOverflow[i].next,MatchOverflow[i].cnt,nkeepcnt);
		}
		fflush(stdout);
		assert(nkeepcnt == keepcnt);
	      }
	    }
	  }
          continue;
        }

	r->next = MatchOverflowFree;
	if(DEBUG>=2) assert(r->next < 0 || r->next > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	if(np == p){
	  if(DEBUG>=2){/* set cnt field of Free list to 0, from MatchOverflow[np->cnt] to r */
	    for(int u = np->cnt;; u = MatchOverflow[u].next){
	      MatchOverflow[u].cnt = 0;
	      if(r == &MatchOverflow[u])
		break;
	    }
	  }

	  MatchOverflowFree = np->cnt;
	  if(DEBUG>=2) assert(MatchOverflowFree > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	  int FreeCnt;
	  if(DEBUG>=2){/* make sure MatchOverflowFree linked list length is <= MatchOverflowNext - MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) */
	    FreeCnt = 0;
	    for(int i = MatchOverflowFree; i > 0; i = MatchOverflow[i].next){
	      if(DEBUG>=2) assert(i > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	      FreeCnt++;
	      if((DEBUG>=3 && ! (FreeCnt <= MatchOverflowNext - MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) && i <= MatchOverflowNext && FreeCnt <= origFreeCnt + OverflowCnt)) || (VERB>=3 && mapid2==MAPID2 && flip2==FLIP2)){
		printf(" tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d:p==np:np->cnt= %d,totcnt= %lld,OverflowFree= %d:i=%d,cnt=%d,FreeCnt=%d->%d,OverCnt=%d,OverflowNext=%d,np=[%ld],np->cnt=%d,r=[%ld],r->next=%d\n", 
	           tid,mapid2,flip2,minSite2,maxSite2,N2,ncnt,totcnt + ncnt,MatchOverflowFree,i,MatchOverflow[i].cnt,origFreeCnt,FreeCnt,OverflowCnt,MatchOverflowNext,np - MatchOverflow, np->cnt,r - MatchOverflow, r->next);
		fflush(stdout);
		assert(FreeCnt <= MatchOverflowNext - MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) + 100 && i <= MatchOverflowNext);
		assert(FreeCnt <= origFreeCnt + OverflowCnt);
	      }
	      if(DEBUG) assert(MatchOverflow[i].cnt == 0);
	    }
	    if(DEBUG) assert(FreeCnt <= MatchOverflowNext - MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	  }
	  if((VERB>=3 && mapid2==MAPID2 && flip2==FLIP2) || (DEBUG>=2 && !(FreeCnt == origFreeCnt + OverflowCnt))){
	    printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d:p==np:np->cnt -> %d, totcnt -> %lld, MatchOverflowFree -> %d,FreeCnt=%d->%d,OverCnt=%d,np=[%ld],np->next=%d,r=[%ld],r->next=%d\n", 
		   tid,mapid2,flip2,minSite2,maxSite2,N2,ncnt,totcnt + ncnt,MatchOverflowFree,origFreeCnt,FreeCnt,OverflowCnt,np - MatchOverflow, np->next,r - MatchOverflow, r->next);
	    fflush(stdout);
	    if(DEBUG) assert(FreeCnt == origFreeCnt + OverflowCnt);
	  }
	  totcnt += ncnt;
	  if(DEBUG>=2 && !(totcnt - origtotcnt == keepcnt && origkeepcnt == keepcnt + deletecnt)){
	    printf("totcnt=%lld, origtotcnt=%lld, keepcnt=%u, ncnt=%u, origkeepcnt=%u, deletecnt=%u\n",totcnt,origtotcnt,keepcnt,ncnt,origkeepcnt,deletecnt);
	    printf("\t roffset=%d,HashGC=%d,minSite2=%d,maxSite2=%d,minOffset2=%d,maxOffset2=%d,max_spread=%d\n", roffset,HashGC,minSite2,maxSite2,minOffset2,maxOffset2,max_spread);
	    fflush(stdout);
	    assert(totcnt - origtotcnt == keepcnt && origkeepcnt == keepcnt + deletecnt);
	  }
	  if(DEBUG>=2) assert(ncnt == keepcnt);
	  if(!(np->cnt = ncnt)){
	    *pnext = qindex;
	    qnext = pnext;
	    continue;
	  }
	} else { /* np !=p && np !=r */
	  if(DEBUG>=2){/* set cnt field of Free list to 0, from MatchOverflow[np->next] to r */
	    for(int u = np->next;; u = MatchOverflow[u].next){
	      MatchOverflow[u].cnt = 0;
	      if(r == &MatchOverflow[u])
		break;
	    }
	  }
	  MatchOverflowFree = np->next;
	  if(DEBUG>=2) assert(MatchOverflowFree > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	  np->cnt = ncnt;
	  np->next = -1;

	  int FreeCnt;
	  if(DEBUG>=2){/* make sure MatchOverflowFree linked list length is <= MatchOverflowNext - MATCHLINE */
	    FreeCnt = 0;
	    for(int i = MatchOverflowFree; i > 0; i = MatchOverflow[i].next){
	      if(DEBUG>=2) assert(i > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	      FreeCnt++;
	      if((DEBUG>=3 && ! (FreeCnt <= MatchOverflowNext - MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) && i <= MatchOverflowNext && FreeCnt <= origFreeCnt + OverflowCnt)) || 
		 (VERB>=3 && mapid2==MAPID2 && flip2==FLIP2 && (MAXSITE2== -1 || maxSite2==MAXSITE2))){
		printf(" tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d:p!=np!=r:np->cnt= %d,totcnt= %lld,OverflowFree= %d:i=%d,cnt=%d,FreeCnt=%d->%d,OverCnt=%d,OverflowNext=%d,np=[%ld],np->next=%d,r=[%ld],r->next=%d\n", 
   	           tid,mapid2,flip2,minSite2,maxSite2,N2,ncnt,totcnt + ncnt,MatchOverflowFree,i,MatchOverflow[i].cnt,origFreeCnt,FreeCnt,OverflowCnt,MatchOverflowNext,np - MatchOverflow, np->next,r - MatchOverflow, r->next);
		fflush(stdout);
		assert(FreeCnt <= MatchOverflowNext - MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) + 100 && i <= MatchOverflowNext);
		//		assert(FreeCnt <= origFreeCnt + OverflowCnt);
	      }
	      if(DEBUG) assert(MatchOverflow[i].cnt == 0);
	    }
	    if(DEBUG) assert(FreeCnt <= MatchOverflowNext - MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	  }
	  if((VERB>=3 && mapid2==MAPID2 && flip2==FLIP2) || (DEBUG>=2 && !(FreeCnt == origFreeCnt + OverflowCnt))){
	    printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d:p!=np!=r:np->cnt -> %d, totcnt -> %lld, MatchOverflowFree -> %d,FreeCnt=%d->%d,OverCnt=%d,np=[%ld],np->next=%d,r=[%ld],r->next=%d\n", 
		   tid,mapid2,flip2,minSite2,maxSite2,N2,ncnt,totcnt + ncnt,MatchOverflowFree,origFreeCnt,FreeCnt,OverflowCnt,np - MatchOverflow, np->next,r - MatchOverflow, r->next);
	    fflush(stdout);
	    if(DEBUG) assert(FreeCnt == origFreeCnt + OverflowCnt);
	  }
	  if(DEBUG>=2) assert(ncnt > 0);
	  totcnt += ncnt;
	  if(DEBUG>=2) assert(totcnt - origtotcnt == keepcnt && origkeepcnt == keepcnt + deletecnt);
	  if(DEBUG>=2){/* check keepcnt */
	    unsigned int nkeepcnt = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE);
	    for(int i = p->cnt; i > 0; i = MatchOverflow[i].next){
	      if(DEBUG) assert(MatchOverflow[i].cnt > 0);
	      nkeepcnt += MatchOverflow[i].cnt;
	    }
	    if(DEBUG) assert(nkeepcnt == keepcnt);
	  }
	}

	continue;
      }

      /* just check cnt entries in p */
      register unsigned int ncnt = 0;
      if(DEBUG>=2) origkeepcnt = cnt;
      if(flip2==0){
	for(register unsigned int i = 0; i < cnt; i++){
	  int mapid1 = p->mapid1[i];
	  OFFSET_T offset = p->offset[i];
	  if(offset + Offset1flip2[mapid1] < Offset2){
	    if(DEBUG>=2) deletecnt++;
	    if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	      printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,i=%d:Deleting mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld\n",
		     tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,i,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt);
	      fflush(stdout);
	    }
	    continue;
	  }

	  p->mapid1[ncnt] = mapid1;
	  p->offset[ncnt] = offset;
	  p->score[ncnt] = p->score[i];
	  if(RANGE >= 1)
	    p->Offset1[ncnt] = p->Offset1[i];
	  if(HASH_SCALE >= 1)
	    p->scaleID[ncnt] = p->scaleID[i];
	  ncnt++;
	  if(DEBUG>=2) keepcnt++;
	  if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	    printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,i=%d:Keeping mapid1=%d,offset=%d,score=%d(Offset2=%d,maxOffset1=%d):ncnt=%d,totcnt=%lld\n",
		   tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,i,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt);
	    fflush(stdout);
	  }
	}
      } else {// flip2==1
	for(register unsigned int i = 0; i < cnt; i++){
	  int mapid1 = p->mapid1[i];
	  OFFSET_T offset = p->offset[i];
	  if(offset + Offset1flip2[mapid1] > Offset2){
	    if(DEBUG>=2) deletecnt++;
	    if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	      printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,i=%d:Deleting mapid1=%d,offset=%d,score=%d(Offset2=%d,minOffset1=%d):ncnt=%d,totcnt=%lld\n",
		     tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,i,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt);
	      fflush(stdout);
	    }
	    continue;
	  }

	  p->mapid1[ncnt] = mapid1;
	  p->offset[ncnt] = offset;
	  p->score[ncnt] = p->score[i];
	  if(RANGE >= 1)
	    p->Offset1[ncnt] = p->Offset1[i];
	  if(HASH_SCALE >= 1)
	    p->scaleID[ncnt] = p->scaleID[i];
	  ncnt++;
	  if(DEBUG>=2) keepcnt++;
	  if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	    printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: pindex=%d,%d,qindex=%d,p->next=%d,p->cnt=%d,i=%d:Keeping mapid1=%d,offset=%d,score=%d(Offset2=%d,minOffset1=%d):ncnt=%d,totcnt=%lld\n",
		   tid,mapid2,flip2,minSite2,maxSite2,N2,pindex,*pnext,qindex,p->next,p->cnt,i,mapid1,offset,p->score[i],Offset2,Offset1flip2[mapid1],ncnt,totcnt);
	    fflush(stdout);
	  }
	}
      }
      totcnt += ncnt;
      if(DEBUG>=2) assert(totcnt - origtotcnt == keepcnt && origkeepcnt == keepcnt + deletecnt);
      if(DEBUG>=2) assert(ncnt == keepcnt);
      if((p->cnt = ncnt) <= 0){/* remove MatchIndex[*pnext] from MatchFirst linked list */
        *pnext = qindex;
	qnext = pnext;/* pnext stays same in next iteration */
	continue;
      }
    }

    if(VERB>=2){
      #pragma omp critical
      {
        printf("MatchFind(mapid2=%d(id=%lld),flip2=%d,minSite2=%d,maxSite2=%d,N2=%d:tid=%d:Offset2=%d,N2=%d,X[N2+1]= %0.3f,minSite2=%d -> %d, MatchCnt= %lld -> %lld\n",
	       mapid2,qmap[mapid2]->id,flip2,minSite2,maxSite2,N2,tid,Offset2,N2,X[N2+1], minSite2, maxSite2 - max_spread + 1, MatchCnt, totcnt);
	fflush(stdout);
      }
    }

    MatchCnt = totcnt;

    if(DEBUG>=2){/* verify accuracy of MatchCnt */
      long long TrueCnt1 = 0, TrueCnt2 = 0, TrueCnt3 = 0;
      int IndexCnt = 0, FreeCnt = 0;    
      register int nindex = MatchFirst;
      register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *q = &MatchIndex[nindex], *p;
      for(register int mindex = nindex; mindex >= 0; mindex = nindex){
	p = q;
	q = &MatchIndex[nindex = p->next];

	unsigned int cnt = p->cnt;
	if(DEBUG>=2) assert(cnt > 0);
	if(cnt > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){
	  for(int i = cnt; i > 0; i = MatchOverflow[i].next){
	    if(DEBUG>=2 && !(MatchOverflow[i].cnt > 0)){
	      #pragma omp critical
	      {
		printf("tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d:mindex=%d,nindex=%d,p->cnt=%d:i=%d,MatchOverflow[i].next=%d,cnt=%d\n",
		     tid,mapid2,flip2,minSite2,maxSite2,mindex,nindex,p->cnt,i,MatchOverflow[i].next,MatchOverflow[i].cnt);
		fflush(stdout);
	      
		assert(MatchOverflow[i].cnt > 0);
	      }
	    }
	  }
	}

	IndexCnt++;
	TrueCnt1 += min(MATCHLINE(OFFSET_T,RANGE,HASH_SCALE),p->cnt);

	if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	  printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: mindex=%d,nindex=%d,p->next=%d,p->cnt=%u(mapid1=%d,offset=%d,score=%u):IndexCnt=%d,FreeCnt=%d(Free=%d),TrueCnt=%lld,%lld,%lld\n",
		 tid,mapid2,flip2,minSite2,maxSite2,N2,mindex,nindex,p->next,p->cnt,p->mapid1[0],p->offset[0],p->score[0],IndexCnt,FreeCnt,MatchOverflowFree,TrueCnt1,TrueCnt2,TrueCnt3);
	  fflush(stdout);
	}
      }

      p = MatchOverflow;
      for(int i = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)+1; i < MatchOverflowNext; i++){
	if(DEBUG>=2 && !HashGC && !(0 < p[i].cnt && p[i].cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE))){
          #pragma omp critical
	  {
	    printf("tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d:MatchCnt=%lld,IndexCnt=%d,Overflow Lines=%d:MatchOverflow[i=%d]:cnt=%d,next=%d\n",
		   tid,mapid2,flip2,minSite2,maxSite2,MatchCnt,IndexCnt, MatchOverflowNext-MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)-1,i,p[i].cnt,p[i].next);
	    fflush(stdout);

	    assert(0 < p[i].cnt && p[i].cnt <= MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	  }  
        }
	TrueCnt2 += p[i].cnt;

	if(VERB>=3 && mapid2==MAPID2 && flip2==FLIP2){
	  printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: i=%d/%d,p[i].cnt=%u(mapid1=%d,offset=%d,score=%u):IndexCnt=%d,FreeCnt=%d(Free=%d),TrueCnt=%lld,%lld,%lld\n",
	    tid,mapid2,flip2,minSite2,maxSite2,N2,i,MatchOverflowNext,p[i].cnt,p[i].mapid1[0],p[i].offset[0],p[i].score[0],IndexCnt,FreeCnt,MatchOverflowFree,TrueCnt1,TrueCnt2,TrueCnt3);
	  fflush(stdout);
	}
      }

      if(HashGC) { /* subtract p[i].cnt from freelist starting at MatchOverflowFree */
        for(int i = MatchOverflowFree; i > 0; i = MatchOverflow[i].next){
	  if(DEBUG>=2) assert(i > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
	  TrueCnt3 += MatchOverflow[i].cnt;
	  FreeCnt++;
	  if((VERB>=3 && mapid2==MAPID2 && flip2==FLIP2) || (DEBUG && FreeCnt >= MatchOverflowNext)){
	    printf("\t tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d: i=%d,next=%d,p[i].cnt=%u:IndexCnt=%d,FreeCnt=%d(Free=%d),TrueCnt=%lld,%lld,%lld, MatchOverflowNext=%d,MatchOverflowMax=%d\n",
	    tid,mapid2,flip2,minSite2,maxSite2,N2,i,MatchOverflow[i].next,MatchOverflow[i].cnt,IndexCnt,FreeCnt,MatchOverflowFree,TrueCnt1,TrueCnt2,TrueCnt3,MatchOverflowNext,MatchOverflowMax);
	    fflush(stdout);
	    if(DEBUG) assert(!(FreeCnt >= MatchOverflowNext + 100));
	  }
	}
	if(DEBUG) assert(!(FreeCnt >= MatchOverflowNext));
      }

      if(VERB>=2 || TrueCnt1+TrueCnt2-TrueCnt3 != MatchCnt){
        #pragma omp critical
        {
	  printf("tid=%d:mapid2=%d,flip2=%d,minSite2=%d,maxSite2=%d,N2=%d,:MatchCnt=%lld,TrueCnt=%lld(%lld+%lld-%lld),IndexCnt=%d,FreeCnt=%d Overflow Lines=%d\n",
	    tid,mapid2,flip2,minSite2,maxSite2,N2,MatchCnt, TrueCnt1+TrueCnt2-TrueCnt3, TrueCnt1, TrueCnt2, TrueCnt3,IndexCnt, FreeCnt,MatchOverflowNext-MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)-1);
	  fflush(stdout);
	  assert(TrueCnt1+TrueCnt2-TrueCnt3 == MatchCnt);
        }
      }
    } // if(DEBUG>=3)

    if(RANGE && ShiftOffset1 <= 30){
      /* restore score values for remaining entries from MatchBuf :  (flip2 ? minOffset2 <= offset + Offset1flip2[mapid1] <= Offset2 : Offset2 <= offset + Offset1flip2[mapid1] <= maxOffset2) */ 
      int MinOffset2 = (flip2 ? minOffset2 : Offset2);
      int MaxOffset2 = (flip2 ? Offset2 : maxOffset2);

      pBuf = MatchBuf;
      _mm_prefetch(CAST pBuf, _MM_HINT_T0);
      _mm_prefetch(CAST &pBuf[8], _MM_HINT_T0);
      for(i = 0; i < PREFETCH_D && i < nMatchCnt; i++)  /* prefetch first PREFETCH_D Offset1flip2[mapid] Entries */
	_mm_prefetch(CAST &Offset1flip2[pBuf[i].mapid1], _MM_HINT_T0);

      for(i = 0; i < orignMatchCnt - PREFETCH_D; i++){/* main loop */
	if(PREFETCH_D) _mm_prefetch(CAST &Offset1flip2[pBuf[PREFETCH_D].mapid1], _MM_HINT_T0);
	register int mapid1 = pBuf[i].mapid1;
	register OFFSET_T offset = pBuf[i].offset;
	register OFFSET_T offset2 = offset + Offset1flip2[mapid1];
	if(!(MinOffset2 <= offset2 && offset2 <= MaxOffset2))
	  continue;
	MatchScoreUpdate(mapid1,offset,pBuf[i].MinOffset1[0], pBuf[i].scaleID[0], pBuf[i].MaxOffset1[0],mapid2,flip2);
      }
      for(; i < orignMatchCnt; i++){/* post loop */
	register int mapid1 = pBuf[i].mapid1;
	register OFFSET_T offset = pBuf[i].offset;
	register OFFSET_T offset2 = offset + Offset1flip2[mapid1];
	if(!(MinOffset2 <= offset2 && offset2 <= MaxOffset2))
	  continue;
	MatchScoreUpdate(mapid1,offset,pBuf[i].MinOffset1[0], pBuf[i].scaleID[0], pBuf[i].MaxOffset1[0],mapid2,flip2);
      }

      /* reduce memory used in MatchBuf[] array */

      if(DEBUG && MatchCnt > 0) assert(MatchBufMax > 0);
      if(MatchBufMax > 0){
	// WAS delete [] MatchBuf;
	size_t memsiz = (sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>) * MatchBufMax + PAGE-1) & ~(PAGE-1);
	if(munmap(MatchBuf, memsiz)){
	  int eno = errno;
	  char *err = strerror(eno);

          #pragma omp critical
	  {
	    printf("munmap(%p,%lu) failed: errno=%d:%s\n",MatchBuf, memsiz, eno, err);
	    dumpmemmap();
	    fflush(stdout);exit(1);
          }
	}

	if(hashtable->MatchMemReserve)
	  hashtable->MatchMemReserve->ResFree(MatchBufMax * sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), 1, this);
      }
      MatchBufMax = 0;
      MatchBuf = 0;
    }// if(RANGE && ShiftOffset1 <= 30)

    /* Update minSite2 */
    minSite2 = maxSite2;
    
    return;
  }

  if(DEBUG) assert(!HashGC || maxSite2 <= 0);
  if(HashGC)
    minSite2 = 0;// reset value

  long long totBest = BestCnt;
  
  if(!HashGC || maxSite2 <= 0){
    long long origtotBest = BestCnt;

    long long dupcnt = 0;

    if(HashMultiMatch){/* sort all matches for the same mapid1 by offset and filter out matches with offsets within HashMultiMatch, keeping only the locally best scoring ones.
			  HERE HERE : Scaling differences are ignored, even though the same offset values don't mean the same thing for different scalings : the scaling with the highest score is retained
			  If RANGE>=1 this can be improved : check if true (unscaled) offset values differ at either end of the matchgroup by HashMultiMatch, if so do not filter out the lower scoring one. 
		       */
                     
      origtotBest = totBest = 0;

      for(int i = 0; i < BestCnt; i++){
	int mapid1 = BestMapid[i];
	CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *q = BestMatch[mapid1];
	int cnt = BestMatchCnt[mapid1];
	if(VERB>=2)
	  origtotBest += cnt;
	if(DEBUG>=2) assert(BestMatchMax[mapid1] == cnt);
	if(HASH_SCALE>=1 && OFFSET_SCALE_FIX>=2){// probe map (ref) was scaled by ScaleFactor[p[k].scaleID[0]] : correct offset by undoing this scaling, since output assumes ref map is unscaled and query is scaled
          if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 /* && flip2 == TraceOR2 && cnt > 0 */){
	    printf("All %d matches for mapid1=%d(id=%lld), mapid2=%d(id=%lld), flip2=%d\n", cnt, mapid1, rmap[mapid1]->id, mapid2, qmap[mapid2]->id, flip2);
	    fflush(stdout);
	  }

	  for(int k = 0; k < cnt; k++){
	    int scaleID = q[k].scaleID[0];
	    int bscale = (scaleID + 1) / 2;
	    int iscale = (scaleID + 1) % 2;
	    if(DEBUG>=2 && !hashScaleDelta){
	      if(DEBUG) assert(bscale == 0);
	      if(DEBUG) assert(iscale == 1);
	    }
	    int rscaleID = (bscale==0) ? 0 : (bscale * 2 + (iscale ? 0 : 1)) - 1;
	    if(DEBUG>=2 && !(rscaleID >= 0 && (rscaleID + 1) / 2 == bscale && ((scaleID==0 && rscaleID==0) || ((rscaleID + 1) % 2 != iscale)))){
	      printf("q[k=%d].scaleID= %d, bscale=%d, iscale=%d, rscaleID=%d\n", k,scaleID,bscale,iscale,rscaleID);
	      fflush(stdout);
	      assert(rscaleID >= 0 && (rscaleID + 1) / 2 == bscale && ((bscale==0 && rscaleID==0) || ((rscaleID + 1) % 2 != iscale)));
	    }
	    if(DEBUG>=2) assert(fabs((1.0/ScaleFactor[scaleID]) - ScaleFactor[rscaleID]) < 0.001);

	    if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 /*  && flip2 == TraceOR2 */){
	      printf("  k=%d:score=%d, scaleID= %d, bscale=%d, iscale=%d, rscaleID=%d: offset= %d -> %d\n",
		     k, q[k].score, q[k].scaleID[0], bscale,iscale,rscaleID,q[k].offset, (int)floor(q[k].offset * ScaleFactor[rscaleID] + 0.5));
	      fflush(stdout);
	    }
	    q[k].offset = floor(q[k].offset * ScaleFactor[rscaleID] + 0.5);
	  }
	}

	if(HASH_SCALE>=1 || RANGE>=1){/* NEW401: first combine entries (add scores) with same offset and scaling : may be due to re-scaling OR different Offset1 values*/
	  qsort(q,cnt,sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), (intcmp*)CBestMatchOffsetScaleInc<OFFSET_T,RANGE,HASH_SCALE>);	  
	  
	  if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 /* && flip2 == TraceOR2*/){
	    if(RANGE>=1)
	      printf("After sorting by offset,scaleID, Offset1:  %d matches for mapid1=%d(id=%lld), mapid2=%d(id=%lld), flip2=%d\n", cnt, mapid1, rmap[mapid1]->id, mapid2, qmap[mapid2]->id, flip2);
	    else
	      printf("After sorting by offset, scaleID:  %d matches for mapid1=%d(id=%lld), mapid2=%d(id=%lld), flip2=%d\n", cnt, mapid1, rmap[mapid1]->id, mapid2, qmap[mapid2]->id, flip2);
	    for(int k = 0; k < cnt; k++){
	      if(RANGE>=1)
		printf("  k=%d:offset= %d, scaleID= %d, Offset1= %u..%u, score= %d\n",
		       k, q[k].offset,q[k].scaleID[0], q[k].MinOffset1[0], q[k].MaxOffset1[0],q[k].score);
	      else
		printf("  k=%d:offset= %d, scaleID= %d, score= %d\n",
		       k, q[k].offset, q[k].scaleID[0], q[k].score);
	    }
	    fflush(stdout);
	  }

	  OFFSET_T offset = q[0].offset;
	  unsigned char scaleID = q[0].scaleID[0];
	  int j = 0;

	  for(int k = j+1;k < cnt; k++){
	    if(q[k].offset == offset && q[k].scaleID[0] == scaleID){
	      q[j].score += q[k].score;
	      if(RANGE>=1){
		q[j].MinOffset1[0] = min(q[j].MinOffset1[0],q[k].MinOffset1[0]);
		q[j].MaxOffset1[0] = max(q[j].MaxOffset1[0],q[k].MaxOffset1[0]);
	      }
	      continue;
	    }

	    offset = q[k].offset;
	    scaleID = q[k].scaleID[0];
	    q[++j] = q[k];
	  }

	  if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 /* && flip2 == TraceOR2 && cnt > 0*/ ){
	    printf("Removed %d entries with duplicate offset values for mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d: cnt= %d -> %d\n",cnt-j,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,cnt,j);
	    fflush(stdout);
	  }
	  cnt = j;
	}

	qsort(q,cnt,sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), (intcmp*)CBestMatchOffsetInc<OFFSET_T,RANGE,HASH_SCALE>);

	if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 /* && flip2 == TraceOR2*/){
	  printf("After merging same offset,scaleID matches:  %d matches for mapid1=%d(id=%lld), mapid2=%d(id=%lld), flip2=%d\n", cnt, mapid1, rmap[mapid1]->id, mapid2, qmap[mapid2]->id, flip2);
	  for(int k = 0; k < cnt; k++){
	    if(RANGE>=1)
	      printf("  k=%d:offset= %d, score= %d, scaleID= %d, Offset1= %u..%u\n",
		     k, q[k].offset,q[k].score,q[k].scaleID[0], q[k].MinOffset1[0], q[k].MaxOffset1[0]);
	    else
	      printf("  k=%d:offset= %d, score= %d, scaleID= %d\n",
		     k, q[k].offset, q[k].score, q[k].scaleID[0]);
	  }
	  fflush(stdout);
	}

	if(HASH_SCALE>=1){/* merge neighboring entries with same offset. These are typically for different scalings :
			     Pick the first one (which will be the highest scoring one), breaking ties in favor of smallest scaleID value if SCALE_SMALL */
	  OFFSET_T offset = q[0].offset;
	  int j = 1;
	  for(; j < cnt; j++){	  
	    if(q[j].offset <= offset){
	      if(DEBUG>=2) assert(offset == q[j-1].offset);
	      if(DEBUG>=2) assert(q[j].score <= q[j-1].score);
	      if(DEBUG>=2 && SCALE_SMALL) assert(q[j].scaleID[0] >= q[j-1].scaleID[0]);
	      
	      break;
	    }
	    offset = q[j].offset;
	  }

	  for(int k = j+1;k < cnt; k++){
	    if(q[k].offset <= offset){
	      if(DEBUG>=2) assert(offset == q[j-1].offset);
	      if(DEBUG>=2) assert(q[k].score <= q[j-1].score);
	      if(DEBUG>=2 && SCALE_SMALL) assert(q[k].scaleID[0] >= q[j-1].scaleID[0]);
	      continue;
	    }
	    offset = q[k].offset;

	    q[j++] = q[k];
	  }

	  if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && cnt > j){
	    printf("Removed %d entries with duplicate offset values for mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d: cnt= %d -> %d\n",cnt-j,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,cnt,j);
	    fflush(stdout);
	  }
	  cnt = j;
	}

	if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 /* && flip2 == TraceOR2 && cnt > 0 */){
          #pragma omp critical
	  {
	    printf("All %d matches for mapid1=%d(id=%lld), mapid2=%d(id=%lld), flip2=%d\n", cnt, mapid1, rmap[mapid1]->id, mapid2, qmap[mapid2]->id, flip2);
	    for(int t = 0; t < cnt; t++){
	      if(RANGE>=1)
		printf("offset= %d, Offset1= %u..%u, sc=%d, score= %d\n", 
		       q[t].offset, q[t].MinOffset1[0], q[t].MaxOffset1[0], (HASH_SCALE && hashScaleDelta) ? q[t].scaleID[0] : 0, q[t].score);
	      else
		printf("offset= %d, sc=%d, score= %d\n", 
		     q[t].offset, (HASH_SCALE && hashScaleDelta) ? q[t].scaleID[0] : 0, q[t].score);
	    }
	    fflush(stdout);
	  }
	}

	if(OFFSET_FILL && !(HASH_SCALE>=1 /* && OFFSET_SCALE_FIX>=2*/) && (!OFFSET_ZEROFIX || DEBUG>=2)){ /* NEW3 : remove duplicate scores with the same offset (and scaleID) */
	  int j = 0;
	  for(int t = 0; t < cnt - 1; t++){
	    if(q[t].offset == q[t+1].offset && (HASH_SCALE <= 0 || q[t].scaleID[0] == q[t+1].scaleID[0])){
	      if(DEBUG) assert((!RANGE && !OFFSET_ZEROFIX) && q[t].score == q[t+1].score);/* duplicate entry */;
	      continue;
	    }
	    q[j++] = q[t];
	  }
	  q[j++] = q[cnt-1];
	  if(VERB>=2)
	    dupcnt += cnt-j;
	  cnt = j;
	}

	// Confirm highest score, then remove any lower scores within HashMultiMatch (same scaleID), then repeat with next highest remaining score, until all scores are either confirmed or deleted
	// Matches to be removed are marked with mapid = -1, confirmed scores are marked with mapid = -2
	if(DEBUG>=2)
	  for(int t = 0; t < cnt; t++)
	    assert(q[t].mapid1 == mapid1);

	while(1){
	  /* locate q[k] with highest score that is not yet confirmed or marked for deletion */
	  /* NOTE : score ties are broken in favor of the smaller offset (and for same offset, smallest scaleID, provided SCALE_SMALL!=0) */
	  short score = 0;
	  int k = -1;
	  for(int t = 0; t < cnt; t++){
	    if(q[t].mapid1 < 0)
	      continue;
	    if(q[t].score > score){
	      score = q[t].score;
	      k = t;
	    }
	  }
	
	  if(k < 0)
	    break;

	  q[k].mapid1 = -2;

	  int offset = q[k].offset;
	  int min_offset = offset - HashMultiMatch;
	  int max_offset = offset + HashMultiMatch;

	  if(TRACE>=2 && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2){
            #pragma omp critical
	    {
	      printf("mapid1=%d(id=%lld(,mapid2=%d(id=%lld),flip2=%d:next best score at k=%d:score= %d, offset=%d: suppressing offset = %d .. %d\n",
		     mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,k,score,offset,min_offset,max_offset);
	      fflush(stdout);
	    }
	  }

	  for(int L = k; --L >= 0;){
	    if(q[L].mapid1 < 0)
	      continue;
	    if(q[L].offset < min_offset)
	      break;

	    if(DEBUG>=2) assert(q[L].score < score);
	    if(RANGE>=1){
	      q[k].MinOffset1[0] = min(q[k].MinOffset1[0], q[L].MinOffset1[0]);
	      q[k].MaxOffset1[0] = max(q[k].MaxOffset1[0], q[L].MaxOffset1[0]);
	    }
	    q[L].mapid1 = -1;/* will be removed later */
	  }
	  for(int R = k; ++R < cnt; ){
	    if(q[R].mapid1 < 0)
	      continue;
	    if(q[R].offset > max_offset)
	      break;

	    if(DEBUG>=2) assert(q[R].score <= score);

	    if(RANGE >= 1){
	      q[k].MinOffset1[0] = min(q[k].MinOffset1[0], q[R].MinOffset1[0]);
	      q[k].MaxOffset1[0] = max(q[k].MaxOffset1[0], q[R].MaxOffset1[0]);
	    }
	    q[R].mapid1 = -1;/* will be removed later */
	  }
	}

	/* now remove matches marked for deletion and restore mapid1 for confirmed scores */
	int tot = 0;
	for(int k = 0; k < cnt; k++){
	  if(q[k].mapid1 == -1)
	    continue;
	  q[k].mapid1 = mapid1;
	  q[tot++] = q[k];
	}

	if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 /* && cnt > 0 */){
          #pragma omp critical
	  {
	    printf("mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d:tot=%d matches remain after filtering with HashMultiMatch=%d:\n",
		   mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,tot,HashMultiMatch);
	    for(int t = 0; t < tot; t++){
	      if(RANGE>=1)
		printf("\t q[%d]: offset= %d, Offset1=%u..%u, sc= %d, score= %d\n", 
		       t, q[t].offset, q[t].MinOffset1[0], q[t].MaxOffset1[0], (HASH_SCALE && hashScaleDelta) ? q[t].scaleID[0] : 0, q[t].score);
	      else
		printf("\t q[%d]: offset= %d, sc= %d, score= %d\n", 
		       t, q[t].offset, (HASH_SCALE && hashScaleDelta) ? q[t].scaleID[0] : 0, q[t].score);
	    }
	    fflush(stdout);
	  }
	}
	if(DEBUG>=2){
	  for(int u = 1; u < tot; u++){
	    if(!(q[u-1].offset + HashMultiMatch <= q[u].offset || (HASH_SCALE && q[u-1].scaleID[0] != q[u].scaleID[0]))){
              #pragma omp critical
	      {
		printf("mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,tot=%d: After filtering with HashMultiMatch=%d:\n",
		       mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,tot,HashMultiMatch);

		for(int t = 0; t < tot; t++)
		  printf("q[%d]: offset= %d, sc= %d, score= %d\n", t, q[t].offset, (HASH_SCALE && hashScaleDelta) ? q[t].scaleID[0] : 0, q[t].score);
		printf("    u=%d,q[u-1].offset=%d,q[u-1].scaleID=%d; q[u].offset=%d,q[u].scaleID=%d\n", 
	          u,q[u-1].offset,(HASH_SCALE && hashScaleDelta) ? q[u-1].scaleID[0] : 0, q[u].offset, (HASH_SCALE && hashScaleDelta) ? q[u].scaleID[0]:0);
		fflush(stdout);

		assert(q[u-1].offset + HashMultiMatch <= q[u].offset || (HASH_SCALE && q[u-1].scaleID[0] != q[u].scaleID[0]));
	      }
	    }
	  }
	}

	totBest += BestMatchCnt[mapid1] = tot;
      }
    } // if(HashMultiMatch)

    if(TRACE && qmap[mapid2]->id == TraceID2){
#pragma omp critical
      {
        printf("MatchFind(mapid2=%d(id=%lld),flip2=%d,M=%d):tid=%d:MatchCnt=%lld:mapid1's matched=%d, total matches=%lld->%lld, duplicate nulls= %lld\n",
  	  mapid2,qmap[mapid2]->id,flip2,qmap[mapid2]->numsite[hcolor], tid,MatchCnt,BestCnt,origtotBest,totBest,dupcnt);
        fflush(stdout);
      }
    }
  } // if(!MULTIMATCH_GC || maxSite2 <= 0)

  if(STAT && BestCnt > BestCntHWM)
    BestCntHWM = BestCnt;

  /* reset MatchOverflow[] & MatchIndex[] */
  if(1){ /* reduce memory used in MatchOverflow[] array */
    if(MatchOverflowMax > (1 << (MHASH_BITS-2))){
      // WAS free(MatchOverflow);
      size_t memsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * MatchOverflowMax + PAGE-1) & ~(PAGE-1);
      if(munmap(MatchOverflow, memsiz)){
	int eno = errno;
	char *err = strerror(eno);
	
        #pragma omp critical
	{
	  printf("munmap(%p,%lu) failed: errno=%d:%s\n",MatchOverflow, memsiz, eno, err);
	  dumpmemmap();
	  fflush(stdout);exit(1);
        }
      }
      
      if(VERB>=2){
        #pragma omp critical
	{
	  printf("tid=%d: MatchOverflowMax=%d -> %d\n",tid,MatchOverflowMax,(1 << (MHASH_BITS-2)));
	  fflush(stdout);
	}
      }

      if(hashtable->MatchMemReserve)
	hashtable->MatchMemReserve->ResFree((MatchOverflowMax - (1 << (MHASH_BITS-2))) * sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>), 0, this);

      MatchOverflowMax = (1 << (MHASH_BITS-2));

      // WAS POSIX_MEMALIGN(MatchOverflow,64,sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>)*MatchOverflowMax);
      memsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * MatchOverflowMax + PAGE-1) & ~(PAGE-1);
      char *newmem = (char *) mmap(NULL, memsiz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
      if(newmem == MAP_FAILED){
	int eno = errno;
	char *err = strerror(eno);
	    
#pragma omp critical
	{
	  printf("CMatchTable: mmap of %lu bytes failed:errno=%d:%s\n", memsiz,eno,err);
	  dumpmemmap();
	  fflush(stdout);exit(1);
        }
      }
      MatchOverflow = (CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *)newmem;
    }
  } // Else : No need to deallocate Match Table Overflow memory since it is accessed at the low end only

  /* Reset MatchOverflow[] */	
  MatchOverflowNext = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)+1;
  MatchOverflowFree = -1;

  /* Reset MatchIndex[] */
  register int nindex = MatchFirst;
  register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *q = &MatchIndex[nindex], *p;
  _mm_prefetch(CAST q, _MM_HINT_T0);
  for(register int mindex = nindex; mindex >= 0; mindex = nindex){
    p = q;
    q = &MatchIndex[nindex = p->next];
    _mm_prefetch(CAST q, _MM_HINT_T0);
    p->cnt = 0;
  }
  MatchFirst = -1;
  MatchCnt = 0;

  /* append BestMatch[] Table entries to output buffer matches[0..matchcnt-1], filtering out scores below hashscore threshold */
  short maxscore = MASK(15);/* largest possible hashscore we can save */
  if(matchcnt + totBest >= matchmax){
    if(VERB>=3 || (DEBUG && !(totBest > 0))){
      #pragma omp critical
      {
	printf("MatchFind(mapid2=%d,flip2=%d):tid=%d: Allocating %lu bytes for %lu -> %lld matches[], matchmax = %lu -> %lld\n",
	       mapid2,flip2,tid,(size_t)(matchcnt + totBest) * sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), matchcnt, matchcnt + totBest, matchmax, matchcnt + totBest);
	fflush(stdout);
	if(DEBUG) assert(totBest > 0);
      }
    }
    maxmatchalloc(matchcnt, matchcnt + totBest, matchmax, matches);
    if(STAT && matchmax > OutputHWM)
      OutputHWM = matchmax;
  }

  _mm_prefetch(CAST &matches[matchcnt], _MM_HINT_T0);
  if(DEBUG && OffsetSpread > 0 && !(OFFSET_WT[1] == (1<<HASHSCORE_SHIFT))){
    printf("OFFSET_WT[0,1,2]= %d, %d, %d, HASHSCORE_SHIFT= %d, OffsetSpread=%d\n",OFFSET_WT[0], OFFSET_WT[1], OFFSET_WT[1], HASHSCORE_SHIFT, OffsetSpread);
    fflush(stdout);
    assert(OFFSET_WT[1] == (1<<HASHSCORE_SHIFT));
  }

  if(RANGE)
    hashscore = HashScore * scoreMult;
  if(DEBUG) assert(hashscore == HashScore * scoreMult);

  if(HashMultiMatch){
    if(DEBUG>=1+RELEASE) assert(!Gpairwise);
    for(register int i = 0; i < BestCnt; i++){
      int mapid1 = BestMapid[i];
      int cnt = BestMatchCnt[mapid1];
      if(DEBUG>=2) assert(cnt > 0);
      CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *p = BestMatch[mapid1];

      if(HashMultiMatchMax || HashBest >= 0){
	qsort(p,cnt,sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>), (intcmp*)CBestMatchScoreDec<OFFSET_T,RANGE,HASH_SCALE>);
	if(TRACE>=2 && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2){	
	  printf("mapid2=%d(id1=%lld),mapid1=%d(id2=%lld),or=%d:%d matches in descending order of score:\n",
		 mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,flip2,cnt);
	  for(int t = 0; t < cnt; t++)
	    printf("\t t=%d:score= %d, offset= %d, scaleID= %d\n",t,p[t].score,p[t].offset,p[t].scaleID[0]);

	  fflush(stdout);
	}

	short maxscore = p->score >> HASHSCORE_SHIFT;
	if(HashBest >= 0 && maxscore > besthashscore[mapid1]){
          #pragma omp critical(besthashscore)
	  {
	    if(TRACE && rmap[mapid1]->id == TraceID1){
	      printf("  mapid2=%d(id1=%lld),mapid1=%d(id2=%lld),or=%d:maxscore=%d (at offset=%d) besthashscore[mapid1]=%d -> %d (tid=%d)\n",
	          mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,flip2,maxscore,p->offset,besthashscore[mapid1],maxscore, tid);
	      fflush(stdout);
	    }
	    besthashscore[mapid1] = max(maxscore,besthashscore[mapid1]);
	  }
	}

	if(VERB>=2){
	  #pragma omp critical
	  {
	    if(HashBest >= 0)
	      printf("  mapid2=%d(id1=%lld),mapid1=%d(id2=%lld),or=%d:maxscore=%d,besthashscore[mapid1]=%d,cnt= %d -> %d (tid=%d)\n",
	        mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,flip2,maxscore,besthashscore[mapid1], cnt, HashMultiMatchMax ? min(cnt,HashMultiMatchMax) : cnt, tid);
	    else
	      printf("  mapid2=%d(id1=%lld),mapid1=%d(id2=%lld),or=%d:maxscore=%d,best score=%d,cnt= %d -> %d (tid=%d)\n",
	        mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,flip2,maxscore,maxscore, cnt, HashMultiMatchMax ? min(cnt,HashMultiMatchMax) : cnt, tid);
	    fflush(stdout);
	  }
	}

	if(HashMultiMatchMax) {
	  cnt = min(cnt,HashMultiMatchMax);

	  if(HashMultiMatchMaxDelta){/* filter out values not within Delta of best hashscore */
	    short minscore = maxscore - HashMultiMatchMaxDelta;
	    if(HashBest >= 0)/* filter out values not within HashBest of best known score for mapid1 (so far) */
	      minscore = max(minscore, besthashscore[mapid1] - HashBest);
	    while(cnt >= 1 && (p[cnt-1].score >> HASHSCORE_SHIFT) < minscore){
	      if(TRACE>=2 && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 /* && flip2 == TraceOR2*/){
	        if(HashBest >= 0)
		  printf("filtered matches[%d]:mapid2=%d(id1=%lld),mapid1=%d(id2=%lld),or=%d,hashscore=%d(%d),offset=%d,tid=%d (bestscore=%d,maxscore=%d,minscore=%d)\n",
	             cnt-1,mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,flip2,p[cnt-1].score >> HASHSCORE_SHIFT, p[cnt-1].score,p[cnt-1].offset,tid, besthashscore[mapid1],maxscore,minscore);
		else
		  printf("filtered matches[%d]:mapid2=%d(id1=%lld),mapid1=%d(id2=%lld),or=%d,hashscore=%d(%d),offset=%d,tid=%d (maxscore=%d,minscore=%d)\n",
              	    cnt-1,mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,flip2,p[cnt-1].score >> HASHSCORE_SHIFT, p[cnt-1].score,p[cnt-1].offset,tid, maxscore,minscore);
		fflush(stdout);
	      }

	      cnt--;
	    }
	  }
	}
      }

      _mm_prefetch(CAST p, _MM_HINT_T0);
      _mm_prefetch(CAST &rmap[mapid1], _MM_HINT_T0);

      if(DEBUG>=2) assert(p->mapid1 == mapid1);
      BestMatchMax[mapid1] = BestMatchCnt[mapid1] = 0;/* reset BestMatch[mapid1] */

      for(int k = 0; k < cnt; k++){
	if(p[k].score < hashscore)
	  continue;// NEW400

	/* mapids are output in swapped order so that each output section for a single mapid2 can be sorted seperately (since mapid2 is queried in ascending order) */
	if(DEBUG>=2) assert(matchcnt < matchmax);
	CHashMatch *r = &matches[matchcnt];
	r->orientation = flip2;
	r->id1 = mapid2;
	r->id2 = mapid1;

	if(DEBUG>=2) assert(p[k].score <= maxscore);
	r->hashscore = p[k].score >> HASHSCORE_SHIFT;
	r->offset = floor(p[k].offset * OffsetKB + 0.5); /* offset already reflects reversal of mapid1,mapid2 AND flip2 */
	
	if(RANGE>=1){
	  if(DEBUG>=2) assert(p[k].MinOffset1[0] <= p[k].MaxOffset1[0]);
	  int scale = (1<<ShiftOffset1);
	  r->MinLoc2 = max(0, max(0, p[k].MinOffset1[0] - 2) * scale - 1) * OffsetKB;/* already reflects flip2 */
	  r->MaxLoc2 = min(Ylen[mapid1], (p[k].MaxOffset1[0] + 1) * scale * OffsetKB);/* already reflects flip2 */
	} else
	  r->MinLoc2 = r->MaxLoc2 = 0;

	if(HASH_SCALE>=1){
	  int scaleID = p[k].scaleID[0];
	  int bscale = (scaleID + 1) / 2;
	  int iscale = (scaleID + 1) % 2;
	  if(DEBUG>=2 && !hashScaleDelta){
	    if(DEBUG) assert(bscale == 0);
	    if(DEBUG) assert(iscale == 1);
	  }
	  int rscaleID = (bscale==0) ? 0 : (bscale * 2 + (iscale ? 0 : 1)) - 1;
	  if(DEBUG>=2 && !(rscaleID >= 0 && (rscaleID + 1) / 2 == bscale && ((scaleID==0 && rscaleID==0) || ((rscaleID + 1) % 2 != iscale)))){
	    printf("p[k=%d].scaleID= %d, bscale=%d, iscale=%d, rscaleID=%d\n", k,scaleID,bscale,iscale,rscaleID);
	    fflush(stdout);
	    assert(rscaleID >= 0 && (rscaleID + 1) / 2 == bscale && ((bscale==0 && rscaleID==0) || ((rscaleID + 1) % 2 != iscale)));
	  }
	  r->scaleID = rscaleID;

	  if(OFFSET_SCALE_FIX == 1){// probe map (ref) was scaled by ScaleFactor[p[k].scaleID[0]] : correct offset by undoing this scaling, since output assumes ref map is unscaled and query is scaled
	    // NOTE : this correction should be done earlier to correctly handle -hashMultiMatch  (see OFFSET_SCALE_FIX = 2)
	    if(DEBUG>=2) assert(fabs((1.0/ScaleFactor[scaleID]) - ScaleFactor[rscaleID]) < 0.001);

	    if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && RANGE>=1){
	      printf("  mapid2=%d(id1=%lld,mapid1=%d(id2=%lld),or=%d,hashscore=%d,offset= %0.1f -> %0.1f, sc=%d->%d, Loc2= %0.1f..%0.1f -> %0.3f..%0.1f,tid=%d\n",
		     mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,r->orientation,r->hashscore,p[k].offset * OffsetKB, p[k].offset * OffsetKB * ScaleFactor[rscaleID], 
		     (HASH_SCALE>=1) ? p[k].scaleID[0] : 0, r->scaleID, r->MinLoc2,r->MaxLoc2, r->MinLoc2 * ScaleFactor[rscaleID], r->MaxLoc2 * ScaleFactor[rscaleID],tid);
	      fflush(stdout);
	    }

	    r->offset = floor(p[k].offset * OffsetKB * ScaleFactor[rscaleID] + 0.5);

	    if(RANGE>=1){/* also scale Locations on qry map, which will be scaled */
	      r->MinLoc2 *= ScaleFactor[rscaleID];
	      r->MaxLoc2 *= ScaleFactor[rscaleID];
	    }
	  }
	} else
	  r->scaleID = 0;

	if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 /* && flip2 == TraceOR2 */){
	  if(RANGE >= 1)
	    printf("Appended matches[%lu]:mapid2=%d(id1=%lld),mapid1=%d(id2=%lld),or=%d,hashscore=%d(%d),offset=%d,sc=%d->%d,Loc2=%0.1f..%0.1f,tid=%d\n",
		   matchcnt,mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,r->orientation,r->hashscore,p[k].score,r->offset,HASH_SCALE>=1 ? p[k].scaleID[0]:0,r->scaleID,r->MinLoc2,r->MaxLoc2,tid);
	  else
	    printf("Appended matches[%lu]:mapid2=%d(id1=%lld),mapid1=%d(id2=%lld),or=%d,hashscore=%d(%d),offset=%d,sc=%d->%d,tid=%d\n",
		   matchcnt,mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,r->orientation,r->hashscore,p[k].score,r->offset,HASH_SCALE>=1 ? p[k].scaleID[0]:0,r->scaleID,tid);
	  fflush(stdout);
	}
	matchcnt++;
      }
    }
    if(1){
      delete [] BestMatchMem;
      BestMatchMem = NULL;
    }
  } else {// !HashMultiMatch
    for(register int i = 0; i < BestCnt; i++){
      register int mapid1 = BestMapid[i];
      register CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *p = BestMatch[mapid1];

      _mm_prefetch(CAST p, _MM_HINT_T0);
      //      _mm_prefetch(CAST &rmap[mapid1], _MM_HINT_T0);

      if(DEBUG>=2) assert(p->mapid1 == mapid1);
      p->mapid1 = -1;    /* reset BestMatch[mapid1] */

      if(p->score < hashscore)
	continue;// NEW400

      if(DEBUG>=2) assert(p->score <= maxscore);

      short pscore = p->score >> HASHSCORE_SHIFT;

      if(TRACE && rmap[mapid1]->id == TraceID1 && (HashBest >= 0 || (qmap[mapid2]->id == TraceID2 /* && flip2 == TraceOR2 */))){
	if(HashBest >= 0)
	  printf("BestMapid[%d/%d]:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d:offset=%d, score= %d, pscore= %d, besthashscore[mapid1]= %d, HashBest=%d\n",
		 i,BestCnt,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id, flip2, p->offset, p->score, pscore, besthashscore[mapid1], HashBest);
	else
	  printf("BestMapid[%d/%d]:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d:offset=%d,score= %d, pscore= %d\n",
		 i,BestCnt,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id, flip2, p->offset, p->score, pscore);
	fflush(stdout);
      }

      if(HashBest >= 0){// NOTE : making entire loop one critical section is slower
	short bscore;

	#pragma omp critical(besthashscore)
	{
	  besthashscore[mapid1] = bscore = max(pscore,besthashscore[mapid1]);
        }

	if(Gpairwise){// NEW111
	  short bscore2 = besthashscore[mapid2] = max(pscore,besthashscore[mapid2]);
	  bscore = min(bscore,bscore2);
	}

	if(pscore < bscore - HashBest)
	  continue;
      }

      /* mapids are output in reverse order so that each output section for a single mapid2 can be sorted seperately (since mapid2 is queried in ascending order) */
      register CHashMatch *r = &matches[matchcnt];
      r->orientation = flip2;
      r->id1 = mapid2;
      r->id2 = mapid1;
	
      r->hashscore = pscore;
      r->offset = floor(p->offset * OffsetKB + 0.5); /* offset already reflects reversal of mapid1,mapid2 AND flip2 */
      if(RANGE>=1){
	if(DEBUG>=2) assert(p->MinOffset1[0] <= p->MaxOffset1[0]);
	int scale = (1<<ShiftOffset1);
	r->MinLoc2 = max(0, max(0, p->MinOffset1[0] - 2) * scale - 1) * OffsetKB;/* already reflects flip2 */
	r->MaxLoc2 = min(Ylen[mapid1], (p->MaxOffset1[0] + 1) * scale * OffsetKB);/* already reflects flip2 */
      } else
	r->MinLoc2 = r->MaxLoc2 = 0;

      if(HASH_SCALE>=1){
        int scaleID = p->scaleID[0];
	int bscale = (scaleID + 1) / 2;
	int iscale = (scaleID + 1) % 2;
	if(!hashScaleDelta){
	  if(DEBUG) assert(bscale == 0);
	  if(DEBUG) assert(iscale == 1);
	  /*	  bscale = 0;
		  iscale = 1;*/
	}
	int rscaleID = (bscale==0) ? 0 : (bscale * 2 + (iscale ? 0 : 1)) - 1;
	if(DEBUG) assert(rscaleID >= 0 && (rscaleID + 1) / 2 == bscale && ((scaleID==0 && rscaleID==0) || ((rscaleID + 1) % 2 != iscale)));
	r->scaleID = rscaleID;

	if(OFFSET_SCALE_FIX){// probe map (ref) was scaled by ScaleFactor[p->scaleID[0]] : correct offset by undoing this scaling, since output assumes ref map is unscaled and query is scaled
          if(DEBUG>=2) assert(fabs((1.0/ScaleFactor[scaleID]) - ScaleFactor[rscaleID]) < 0.001);

	  if(VERB>=3){
	    if(RANGE>=1)
	      printf("  mapid2=%d(id1=%lld,mapid1=%d(id2=%lld),or=%d,hashscore=%d,offset= %0.1f -> %0.1f, sc=%d->%d, Loc2= %0.1f..%0.1f -> %0.3f..%0.1f,tid=%d\n",
		     mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,r->orientation,r->hashscore,p->offset * OffsetKB, p->offset * OffsetKB * ScaleFactor[rscaleID], 
		     (HASH_SCALE>=1) ? p->scaleID[0] : 0, r->scaleID, r->MinLoc2,r->MaxLoc2, r->MinLoc2 * ScaleFactor[rscaleID], r->MaxLoc2 * ScaleFactor[rscaleID],tid);
	    else
	      printf("  mapid2=%d(id1=%lld,mapid1=%d(id2=%lld),or=%d,hashscore=%d,offset= %0.1f -> %0.1f, sc=%d->%d, tid=%d\n",
		     mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,r->orientation,r->hashscore,p->offset * OffsetKB, p->offset * OffsetKB * ScaleFactor[rscaleID], 
		     (HASH_SCALE>=1) ? p->scaleID[0] : 0, r->scaleID,tid);
	    fflush(stdout);
	  }

          r->offset = floor(p->offset * OffsetKB * ScaleFactor[rscaleID] + 0.5);

	  if(RANGE>=1){/* also scale Locations on qry map, which will be scaled */
	    r->MinLoc2 *= ScaleFactor[rscaleID];
	    r->MaxLoc2 *= ScaleFactor[rscaleID];
	  }
        }
      } else
	r->scaleID = 0;
 
      if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 /* && flip2 == TraceOR2*/){
	if(RANGE>=1)
	  printf("Appended matches[%lu]:mapid2=%d(id1=%lld),mapid1=%d(id2=%lld),or=%d,hashscore=%d,offset=%d,sc=%d->%d,Loc2=%0.1f,%0.1f,tid=%d\n",
	      matchcnt,mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,r->orientation,r->hashscore,r->offset, HASH_SCALE>=1 ? p->scaleID[0]:0,r->scaleID,r->MinLoc2,r->MaxLoc2,tid);
	else
	  printf("Appended matches[%lu]:mapid2=%d(id1=%lld),mapid1=%d(id2=%lld),or=%d,hashscore=%d,offset=%d,sc=%d->%d,tid=%d\n",
	      matchcnt,mapid2,qmap[mapid2]->id,mapid1,rmap[mapid1]->id,r->orientation,r->hashscore,r->offset,HASH_SCALE>=1 ? p->scaleID[0]:0,r->scaleID,tid);
      }
      matchcnt++;
    }
  }// !HashMultiMatch

  if(TRACE && qmap[mapid2]->id == TraceID2 /* && flip2 == TraceOR2 */ ){
    #pragma omp critical
    {
      printf("MatchFind(mapid2=%d,flip2=%d,M=%d):tid=%d: copied BestMatch[] to output buffer: matchcnt=%lu,matchmax=%lu: wall time=%0.6f\n",mapid2,flip2,qmap[mapid2]->numsite[hcolor],tid,matchcnt,matchmax,wtime());
      fflush(stdout);
    }
  }

  if(DEBUG>=2){/* check that Hit Table and Match Table  and BestMatch Table have been reset */
    for(int hindex = MASK(HHASH_BITS); hindex >= 0; hindex--)
      assert(HitIndex[hindex]==0);
    assert(HitFirst == -1);
    
    for(int mindex = MASK(MHASH_BITS); mindex >= 0; mindex--)
      assert(MatchIndex[mindex].cnt==0);
    assert(MatchFirst == -1);
    assert(MatchCnt == 0);
    assert(MatchOverflowNext == MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)+1);
    
    if(HashMultiMatch){
      for(int mapid1 = numrmaps + roffset; --mapid1 >= roffset; ){
	assert(BestMatchCnt[mapid1] == 0);
	assert(BestMatchMax[mapid1] == 0);
      }
    } else {
      for(int mapid1 = numrmaps + roffset; --mapid1 >= roffset; )
	assert(BestMatch[mapid1]->mapid1 < 0);
    }
  }

  BestCnt = 0;
}

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::BestInsert(int mapid1, OFFSET_T offset, short score, CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *pBuf, int mapid2, int flip2, int tid)
{
  CBestMatch<OFFSET_T,RANGE,HASH_SCALE> *q = BestMatch[mapid1];
  if(!HashMultiMatch){// accumulate only one score (best offset) per mapid1 */
    if(DEBUG>=2) assert(q->mapid1 == -1 || q->mapid1 == mapid1);

    if(q->mapid1 < 0){/* first match with this mapid1 */
      BestMapid[BestCnt++] = mapid1;/* update list of distinct mapid1 values in BestMatch Table */
    
      if(VERB >= 2 || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	printf("\t BestInsert:tid=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,sc=%d,score=%d:BestCnt=%d\n", 
	  tid, mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,RANGE>=1 ? pBuf->MinOffset1[0] : 0, (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, score,BestCnt);
	fflush(stdout);
      }

      q->mapid1 = mapid1;
      q->offset = offset;
      q->score = score;
      if(RANGE >= 1)
	q->MaxOffset1[0] = q->MinOffset1[0] = pBuf->MinOffset1[0];
      if(HASH_SCALE >= 1)
	q->scaleID[0] = pBuf->scaleID[0];

      return;
    }

    if(RANGE >= 1){ // check if offset (and scaleID) is the same as before : add up the scores
      if(offset == q->offset && (HASH_SCALE <= 0 || pBuf->scaleID == q->scaleID)){
	if(VERB >= 2 || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread)))
	  printf("\t BestInsert:tid=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,sc=%d,Offset1=%u..%u to %u..%u,score=%d to %d\n", 
		 tid,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset, (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0,
		 q->MinOffset1[0], q->MaxOffset1[0], min(q->MinOffset1[0],pBuf->MinOffset1[0]),max(q->MaxOffset1[0],pBuf->MinOffset1[0]), q->score, q->score + score);

	q->score += score;
	q->MinOffset1[0] = min(q->MinOffset1[0],pBuf->MinOffset1[0]);
	q->MaxOffset1[0] = max(q->MaxOffset1[0],pBuf->MinOffset1[0]);

	return;
      }
    } else {
      if(DEBUG>=2) assert(!(offset == q->offset && (HASH_SCALE <= 0 || pBuf->scaleID == q->scaleID)));
    }

    if(score >= q->score && (score > q->score || 
			     offset < q->offset || /* NOTE : breaking ties in favor of smallest offset produces deterministic results when the inserted mapsets change */
			     (SCALE_SMALL && HASH_SCALE>=1 && pBuf->scaleID[0] < q->scaleID[0]) /* NOTE : breaking ties in favor of smallest scaleID */
	  )){
      if(VERB >= 2 || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	printf("\t BestInsert:tid=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,sc=%d->%d,offset=%d -> %d,Offset1=%u..%u -> %u,score=%d to %d\n", 
	       tid,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,(HASH_SCALE && hashScaleDelta) ? q->scaleID[0] : 0, (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0,
	       q->offset,offset,RANGE>=1 ? q->MinOffset1[0]:0, RANGE>=1 ? q->MaxOffset1[0]:0, RANGE>=1 ? pBuf->MinOffset1[0]:0, q->score,score);
	fflush(stdout);
      }

      q->offset = offset;
      q->score = score;
      if(HASH_SCALE >= 1)
	q->scaleID[0] = pBuf->scaleID[0];
      if(RANGE >= 1)
	q->MaxOffset1[0] = q->MinOffset1[0] = pBuf->MinOffset1[0];
    }
    return;
  }

  /* Accumulate multiple matches per mapid1 */
  if(VERB >= 2 || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
    printf("\t BestInsert:tid=%d:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u,sc=%d,score=%d:BestMatchCnt[mapid1]= %d\n", 
	   tid,mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,RANGE>=1 ? pBuf->MinOffset1[0]:0, (HASH_SCALE && hashScaleDelta) ? pBuf->scaleID[0] : 0, score,BestMatchCnt[mapid1] + 1);
    fflush(stdout);
  }

  int cnt = BestMatchCnt[mapid1];
  if(DEBUG>=2) assert(cnt == 0 || q[cnt-1].mapid1 == mapid1);

  if(RANGE>=1){ /* check if any existing entry has the same offset (and scaleID) value */
    for(int i = 0; i < cnt; i++){
      if(q[i].offset == offset && (HASH_SCALE <= 0 || pBuf->scaleID == q->scaleID)){
	if(DEBUG>=2) assert(q[i].mapid1 == mapid1);
	q[i].score += score;
	q[i].MinOffset1[0] = min(q[i].MinOffset1[0], pBuf->MinOffset1[0]);
	q[i].MaxOffset1[0] = max(q[i].MaxOffset1[0], pBuf->MinOffset1[0]);
	if(VERB >= 2 || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread)))
	  printf("\t\t tid=%d:Updating BestMatch[mapid1=%d][i=%d]:offset=%d,sc=%d,Offset1=%u..%u,score=%d\n",
		 tid,mapid1,i,offset,(HASH_SCALE && hashScaleDelta) ? q->scaleID[0] : 0, q[i].MinOffset1[0],q[i].MaxOffset1[0],q[i].score);
	return;
      }
    }
  }

  q[cnt].mapid1 = mapid1;
  q[cnt].offset = offset;
  q[cnt].score = score;
  if(RANGE >= 1)
    q[cnt].MaxOffset1[0] = q[cnt].MinOffset1[0] = pBuf->MinOffset1[0];
  if(HASH_SCALE >=1)
    q[cnt].scaleID[0] = pBuf->scaleID[0];

  BestMatchCnt[mapid1] = ++cnt;
  if(DEBUG>=2) assert(BestMatchMax[mapid1] >= cnt);
}

/* return a score in MatchTable. If no score is founds, return (OFFSET_ZEROFIX ? -1 : 0). 
   Note that return value of -1 will be treated same as a return value 0, but can be used to avoid duplicate 0 scores, with OFFSET_FILL) */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline short CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::MatchScoreFind(int mapid1, OFFSET_T offset, unsigned short Offset1, unsigned char scaleID)
{
  if(DEBUG>=2) assert(roffset <= mapid1 && mapid1 < roffset + numrmaps);
  if(DEBUG>=2) assert(MinOffset <= offset && offset <= MaxOffset);
  if(DEBUG>=2) assert(0 <= scaleID && scaleID < NumScaleFactor);

  register unsigned int mindex = Maphash2[mapid1] ^ Offsethash[offset];
  if(RANGE >= 1)
    mindex ^= Offset1hash[Offset1];
  if(HASH_SCALE >= 1)
    mindex ^= scaleIDhash[scaleID];

  if(DEBUG>=2) assert(mindex < (1U << MHASH_BITS));
  register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *p = &MatchIndex[mindex];
  if(STAT) MTcnt2++;
  register int cnt = p->cnt;
  if(!cnt)
    return (OFFSET_ZEROFIX ? -1 : 0);
  
  if(cnt > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){/* Also check Overflow linked list */
    register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *q = &MatchOverflow[cnt], *r, *s;
    _mm_prefetch(CAST q, _MM_HINT_T0);

    for(register int i = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE); --i >= 0;){
      if(p->mapid1[i] == mapid1 && p->offset[i] == offset && (RANGE<=0 || p->Offset1[i] == Offset1))
	return p->score[i];
      if(STAT) MTqcol++;/* quick collision in current cache line */
    }
    
    for(r = q; r->next >= 0; r = s){
      s = &MatchOverflow[r->next];
      _mm_prefetch(CAST s, _MM_HINT_T0);
      
      if(STAT) MTcol++; /* collision that required a new cache line */
      for(register int i = r->cnt; --i >= 0;){
	if(r->mapid1[i] == mapid1 && r->offset[i] == offset && (RANGE <= 0 || r->Offset1[i] == Offset1))
	  return r->score[i];
	if(STAT) MTqcol++;/* quick collision in current cache line */	
      }
    }    

    if(STAT) MTcol++; /* collision that required a new cache line */
    for(register int i = r->cnt; --i >= 0;){
      if(r->mapid1[i] == mapid1 && r->offset[i] == offset && (RANGE <= 0 || r->Offset1[i] == Offset1))
	return r->score[i];
      if(STAT) MTqcol++;/* quick collision in current cache line */	
    }

    return (OFFSET_ZEROFIX ? -1 : 0);
  }
  for(;--cnt >= 0;){
    if(p->mapid1[cnt] == mapid1 && p->offset[cnt] == offset && (RANGE <= 0 || p->Offset1[cnt] == Offset1))
      return p->score[cnt];
    if(STAT) MTqcol++;/* quick collision in current cache line */
  }
  return (OFFSET_ZEROFIX ? -1 : 0);
}

/* update (change) an existing score in MatchTable (if score is not in table it is an error : do nothing) */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::MatchScoreUpdate(register int mapid1, register OFFSET_T offset, unsigned short Offset1, unsigned char scaleID, unsigned char score, int mapid2, int flip2)
{
  if(DEBUG>=2) assert(roffset <= mapid1 && mapid1 < roffset + numrmaps);
  if(DEBUG>=2) assert(MinOffset <= offset && offset <= MaxOffset);
  if(DEBUG>=2) assert(0 <= scaleID && scaleID < NumScaleFactor);

  register unsigned int mindex = Maphash2[mapid1] ^ Offsethash[offset];
  if(RANGE >= 1)
    mindex ^= Offset1hash[Offset1];
  if(HASH_SCALE >= 1)
    mindex ^= scaleIDhash[scaleID];
  if(DEBUG>=2) assert(mindex < (1U << MHASH_BITS));
  register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *p = &MatchIndex[mindex];
  if(STAT) MTcnt3++;
  register int cnt = p->cnt;
  if(!cnt){/* this is an error */
    if(DEBUG>=1+RELEASE){
      printf("MatchScoreUpdate(mapid1=%d,offset=%d,Offset1=%u,score=%u): Entry not found in MatchTable\n",mapid1,offset,Offset1,score);
      fflush(stdout);
      assert(cnt > 0);
    }
    return;
  }
  if(cnt > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){/* Also check Overflow linked list */
    register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *q = &MatchOverflow[cnt], *r, *s;
    _mm_prefetch(CAST q, _MM_HINT_T0);

    for(register int i = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE); --i >= 0;){
      if(p->mapid1[i] == mapid1 && p->offset[i] == offset && (RANGE <= 0 || p->Offset1[i] == Offset1)){
	if(VERB>=2 || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	  printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u:score=%u replaced with %u\n",
		 mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,RANGE>=1 ? Offset1 : 0,p->score[i],score);
	  fflush(stdout);
	}
	p->score[i] = score;
	return;
      }
      if(STAT) MTqcol++;/* quick collision in current cache line */
    }
    
    for(r = q; r->next >= 0; r = s){
      s = &MatchOverflow[r->next];
      _mm_prefetch(CAST s, _MM_HINT_T0);
      
      if(STAT) MTcol++; /* collision that required a new cache line */
      for(register int i = r->cnt; --i >= 0;){
	if(r->mapid1[i] == mapid1 && r->offset[i] == offset && (RANGE <= 0 || r->Offset1[i] == Offset1)){
	  if(VERB>=2 || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	    printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u:score=%u replaced with %u\n",
		   mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,RANGE>=1 ? Offset1 : 0u,r->score[i],score);
	    fflush(stdout);
	  }
	  r->score[i] = score;
	  return;
	}
	if(STAT) MTqcol++;/* quick collision in current cache line */	
      }
    }    

    if(STAT) MTcol++; /* collision that required a new cache line */
    for(register int i = r->cnt; --i >= 0;){
      if(r->mapid1[i] == mapid1 && r->offset[i] == offset && (RANGE<=0 || r->Offset1[i] == Offset1)){
	if(VERB>=2 || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	  printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%u:score=%u replaced with %u\n",
		 mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,RANGE>=1 ? Offset1 : 0,r->score[i],score);
	  fflush(stdout);
	}
	r->score[i] = score;
	return;
      }

      if(STAT) MTqcol++;/* quick collision in current cache line */	
    }

    if(DEBUG>=1+RELEASE){
      printf("MatchScoreUpdate(mapid1=%d,offset=%d,Offset1=%u,score=%u): Entry not found in MatchTable\n",mapid1,offset,Offset1,score);
      fflush(stdout);
      assert(cnt > 0);
    }
    return;
  }

  for(;--cnt >= 0;){
    if(p->mapid1[cnt] == mapid1 && p->offset[cnt] == offset && (RANGE <= 0 || p->Offset1[cnt] == Offset1)){
      if(VERB>=2 || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	printf("\t mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,offset=%d,Offset1=%d:score=%u replaced with %u\n",
	       mapid1,rmap[mapid1]->id,mapid2,qmap[mapid2]->id,flip2,offset,RANGE >= 1 ? Offset1 : 0,p->score[cnt],score);
	fflush(stdout);
      }
      p->score[cnt] = score;
      return;
    }
    if(STAT) MTqcol++;/* quick collision in current cache line */
  }

  if(DEBUG>=1+RELEASE){
    printf("MatchScoreUpdate(mapid1=%d,offset=%d,Offset1=%u,score=%u): Entry not found in MatchTable\n",mapid1,offset,Offset1,score);
    fflush(stdout);
    assert(cnt > 0);
  }
  return;
}


/**
Updates the Match Table : CMatchTable::MatchIndex and accumulates the hash score for map pairs which support the same alignment offset. 

Indexing is done by reference (insertion) mapid (mapid1) and alignment offset.
> mindex = Maphash2[mapid1] ^ Offsethash[offset] ^ Offset1hash[Offset1] ^ scaleIDhash[scaleID]
> MatchIndex[mindex]

The MatchTable table will only contain entries for a single query (probe) map (mapid2) and orientation (flip2). Unlike HitTable, MatchTable can contain entries based on multiple offset2 and scaleID values.  

Allocate only necessary memory taken using the CMatchTable::MatchOverflowAlloc() method and related objects.
*/
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline unsigned char CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::MatchInsert(register int mapid1, register OFFSET_T offset, register unsigned char score /* ALWAYS 1, unless OFFSET_ZEROFIX */, int mapid2, int flip2, OFFSET_T offset2, unsigned char scaleID, unsigned short Offset1, OFFSET1_T offset1)
{
  if(STAT && score) MTicnt++;

  if(DEBUG>=2) assert(roffset <= mapid1 && mapid1 < roffset + numrmaps);
  if(DEBUG>=2) assert(MinOffset <= offset && offset <= MaxOffset);
  if(DEBUG>=2) assert(0 <= scaleID && scaleID < NumScaleFactor);
  register unsigned int mindex = Maphash2[mapid1] ^ Offsethash[offset];
  if(RANGE >= 1)
    mindex ^= Offset1hash[Offset1];
  if(HASH_SCALE >= 1)
    mindex ^= scaleIDhash[scaleID];
  register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *p = &MatchIndex[mindex];
  if(STAT) MTcnt++;

  if(!p->cnt){/* start new list */
    p->init(mapid1, offset, score, Offset1, scaleID);
    p->next = MatchFirst;
    MatchFirst = mindex;
    MatchCnt++;

    if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))
      printf("\t MatchInsert of mapid1=%d(id=%lld),offset=%d,offset1=%u,Offset1=%d,score=%d at p=%p,i=0,mindex=%u (mapid2=%d(id=%lld),flip2=%d,offset2=%u,scaleID=%u)\n", 
	     mapid1,rmap[mapid1]->id,offset,offset1,Offset1,p->score[0],p,mindex,mapid2,qmap[mapid2]->id,flip2,(unsigned int)offset2,(unsigned int)scaleID);
    return score;
  }

  long long origMTcol;
  if(STAT) origMTcol= MTcol;

  /* see if matching entry is already in MatchTable starting at p */
  const unsigned char lim = 254;
  register unsigned int cnt = p->cnt;
  if(cnt > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){/* Also check Overflow linked list */
    register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *q = &MatchOverflow[cnt],*r,*s;
    _mm_prefetch(CAST q, _MM_HINT_T0);
    
    for(register int i = MATCHLINE(OFFSET_T,RANGE,HASH_SCALE); --i >= 0;){
      if(p->mapid1[i] == mapid1 && p->offset[i] == offset && (RANGE <= 0 || p->Offset1[i] == Offset1) && (HASH_SCALE <= 0 || p->scaleID[i] == scaleID)){
        if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))
	  printf("\t MatchInsert of mapid1=%d(id=%lld),offset=%d,offset1=%u,Offset1=%u,score=%d->%d at p=%p,i=%d,mindex=%u (mapid2=%d(id=%lld),flip2=%d,offset2=%u,scaleID=%u)\n",
		 mapid1,rmap[mapid1]->id,offset,offset1,Offset1,p->score[i],p->score[i]+score,p,i,mindex,mapid2,qmap[mapid2]->id,flip2,(unsigned int)offset2, (unsigned int)scaleID);
	register unsigned char nscore = p->score[i] = min(lim, p->score[i]) + score;
	return nscore;
      }
      if(STAT) MTqcol++;/* quick collision in current cache line */
    }
    
    for(r = q; r->next >= 0; r = s){
      if(DEBUG>=2) assert(r->next == -1 || r->next > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
      s = &MatchOverflow[r->next];
      _mm_prefetch(CAST s, _MM_HINT_T0);
      
      if(STAT) MTcol++; /* collision that required a new cache line */
      if(STAT && DEBUG>=2 && (MTcol < origMTcol || MTcol > origMTcol + 10000)){
	printf("\nMatchInsert of mapid1=%d(id=%lld),offset=%d,offset1=%u,Offset1=%u,mapid2=%d,flip2=%d,offset2=%u,scaleID=%u : MTcol= %lld -> %lld, r->cnt= %u, r->next= %d : suspected cycle in MatchOverflow\n",
	       mapid1,rmap[mapid1]->id,offset,offset1,Offset1,mapid2,flip2, (unsigned int)offset2, (unsigned int)scaleID, origMTcol, MTcol, r->cnt, r->next);
	fflush(stdout);
	assert(!(MTcol < origMTcol || MTcol > origMTcol + 10100));
      }
      for(register int i = r->cnt; --i >= 0;){
	if(r->mapid1[i] == mapid1 && r->offset[i] == offset && (RANGE <= 0 || r->Offset1[i] == Offset1) && (HASH_SCALE <= 0 || r->scaleID[i] == scaleID)){
	  if((STAT && DEBUG>=2 && MTcol > origMTcol + 10000) || (TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread)))
	    printf("\t MatchInsert of mapid1=%d(id=%lld),offset=%d,offset1=%u,Offset1=%u,score=%d->%d at r=%p,i=%d,mindex=%u (mapid2=%d(id=%lld),flip2=%d,offset2=%u,scaleID=%u)\n",
		   mapid1,rmap[mapid1]->id,offset,offset1,Offset1,r->score[i],r->score[i]+score,r,i,mindex,mapid2,qmap[mapid2]->id,flip2,(unsigned int)offset2,(unsigned int)scaleID);

	  register unsigned char nscore = r->score[i] = min(lim, r->score[i]) + score;
	  return nscore;
	}
	if(STAT) MTqcol++;/* quick collision in current cache line */	
      }
    }

    if(STAT) MTcol++; /* collision that required a new cache line */
    for(register int i = r->cnt; --i >= 0;){
      if(r->mapid1[i] == mapid1 && r->offset[i] == offset && (RANGE <= 0 || r->Offset1[i] == Offset1) && (HASH_SCALE <= 0 || r->scaleID[i] == scaleID)){
	if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))
	  printf("\t MatchInsert of mapid1=%d(id=%lld),offset=%d,offset1=%u,Offset1=%u,score=%d->%d at r=%p,i=%d,mindex=%u (mapid2=%d(id=%lld),flip2=%d,offset2=%u,scaleID=%u)\n",
		 mapid1,rmap[mapid1]->id,offset,offset1,Offset1,r->score[i],r->score[i]+score,r,i,mindex,mapid2,qmap[mapid2]->id,flip2,(unsigned int)offset2,(unsigned int)scaleID);
	register unsigned char nscore = r->score[i] = min(lim, r->score[i]) + score;
	return nscore;
      }
      if(STAT) MTqcol++;/* quick collision in current cache line */	
    }

    /* append new entry at q */
    register unsigned int qcnt = q->cnt;
    if(qcnt == MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){/* extend Overflow list */
      register int next = MatchOverflowAlloc(tid);
      q = &MatchOverflow[cnt];/* in case MatchOverflow[] was reallocated */
      r = &MatchOverflow[next];
      _mm_prefetch(CAST r, _MM_HINT_T0);
      if(DEBUG>=2) assert(next > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));

      if(STAT) MTcol++;
      MatchCnt++;

      if(DEBUG>=2) assert(cnt == p->cnt);
      p->cnt = next;
      if(DEBUG>=2) assert(cnt > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));
      r->next = cnt;
      r->init(mapid1, offset, score, Offset1, scaleID);

      if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))
	printf("\t MatchInsert of mapid1=%d(id=%lld),offset=%d,offset1=%u,Offset1=%u,score=%d at r=%p,i=0,mindex=%u (mapid2=%d(id=%lld),flip2=%d,offset2=%u,scaleID=%u)\n",
	       mapid1,rmap[mapid1]->id,offset,offset1,Offset1,r->score[0],r,mindex,mapid2,qmap[mapid2]->id,flip2,(unsigned int)offset2,(unsigned int)scaleID);

      return score;
    }

    MatchCnt++;
    q->mapid1[qcnt] = mapid1;
    q->offset[qcnt] = offset;
    q->score[qcnt] = score;
    if(RANGE >= 1)
      q->Offset1[qcnt] = Offset1;
    if(HASH_SCALE >= 1)
      q->scaleID[qcnt] = scaleID;
    if(DEBUG>=2) assert(qcnt == q->cnt);
    q->cnt = qcnt+1;

    if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))
      printf("\t MatchInsert of mapid1=%d(id=%lld),offset=%d,offset1=%u,Offset1=%u,score=%d at q=%p,i=%d,mindex=%u (mapid2=%d(id=%lld),flip2=%d,offset2=%u,scaleID=%u)\n",
	     mapid1,rmap[mapid1]->id, offset,offset1,Offset1,q->score[qcnt],q,qcnt,mindex, mapid2, qmap[mapid2]->id, flip2, (unsigned int)offset2, (unsigned int)scaleID);

    return score;
  }

  for(register int i = cnt; --i >= 0;){
    if(p->mapid1[i] == mapid1 && p->offset[i] == offset && (RANGE <= 0 || p->Offset1[i] == Offset1) && (HASH_SCALE <= 0 || p->scaleID[i] == scaleID)){
      if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))
	printf("\t MatchInsert of mapid1=%d(id=%lld),offset=%d,offset1=%u,Offset1=%d,score=%d->%d at p=%p,i=%d,mindex=%u (mapid2=%d(id=%lld),flip2=%d,offset2=%u,scaleID=%u)\n",
	       mapid1,rmap[mapid1]->id,offset,offset1,Offset1,p->score[i],p->score[i]+score,p,i,mindex,mapid2,qmap[mapid2]->id,flip2,(unsigned int)offset2,(unsigned int)scaleID);

      register unsigned char nscore = p->score[i] = min(lim,p->score[i]) + score;
      return nscore;
    }
    if(STAT) MTqcol++;/* quick collision in current cache line */
  }

  /* append new entry at p */
  if(cnt == MATCHLINE(OFFSET_T,RANGE,HASH_SCALE)){/* initialize Overflow list */
    register int next = MatchOverflowAlloc(tid);
    register CMatchLine<OFFSET_T,RANGE,HASH_SCALE> *q = &MatchOverflow[next];
    _mm_prefetch(CAST q, _MM_HINT_T0);
    if(DEBUG>=2) assert(next > MATCHLINE(OFFSET_T,RANGE,HASH_SCALE));

    if(STAT) MTcol++;
    MatchCnt++;

    p->cnt = next;
    q->init(mapid1, offset, score, Offset1, scaleID);
    q->next = -1;

    if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))
      printf("\t MatchInsert of mapid1=%d(id=%lld),offset=%d,offset1=%u,Offset1=%u,score=%d at q=%p,i=0,mindex=%u (mapid2=%d(id=%lld),flip2=%d,offset2=%u,scaleID=%u)\n",
	     mapid1,rmap[mapid1]->id,offset,offset1,Offset1,q->score[0],q,mindex,mapid2,qmap[mapid2]->id,flip2,(unsigned int)offset2,(unsigned int)scaleID);

    return score;
  }

  MatchCnt++;
  p->mapid1[cnt] = mapid1;
  p->offset[cnt] = offset;
  if(RANGE >= 1)
    p->Offset1[cnt] = Offset1;
  if(HASH_SCALE >= 1)
    p->scaleID[cnt] = scaleID;
  p->score[cnt] = score;
  p->cnt = cnt+1;

  if(TRACE && rmap[mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2 && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))
    printf("\t MatchInsert of mapid1=%d(id=%lld),offset=%d,offset1=%u,Offset1=%u,score=%d at p=%p,i=%d,mindex=%u (mapid2=%d(id=%lld),flip2=%d,offset2=%u,scaleID=%u)\n",
	   mapid1,rmap[mapid1]->id,offset,offset1,Offset1,p->score[cnt],p,cnt,mindex,mapid2,qmap[mapid2]->id,flip2,(unsigned int)offset2,(unsigned int)scaleID);

  return score;
}

/** compute a hash index of more than HASH_BITS bits from the larger number of bits in key */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline unsigned int CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::hashindex(long long key)
{
  if(DEBUG>=2) assert(rand1[key & rhash_msk] <= (unsigned int)MASK(HASH_BITS));
  if(DEBUG>=2) assert(rand2[(key >> RHASH_BITS) & rhash_msk] <= (unsigned int)MASK(HASH_BITS));
  return rand1[key & rhash_msk] ^ rand2[(key >> RHASH_BITS) & rhash_msk];// sufficient for RHASH_BITS=15,HashWin=5,SIZE_BITS=6
  //    return rand1[key & rhash_msk] ^ rand2[(key >> RHASH_BITS) & rhash_msk] ^ (key >> (RHASH_BITS*2));
};


/* locate next block of specified size in CollideBlock[] : fast implementation does no Garbage collection */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline unsigned int CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::CollideBlockAlloc(unsigned int size)
{
  unsigned int origCollideBlockCnt = CollideBlockCnt;
  CollideBlockCnt += size;
  if(CollideBlockCnt > CollideBlockMax){/* double size of CollideBlock[] */
    unsigned int origCollideBlockMax = CollideBlockMax;
    do {
      if(DEBUG && CollideBlockMax > MASK(31)){	
	printf("Cannot double size of CollideBlock from %u using 32 bit unsigned int\n", CollideBlockMax);
	fflush(stdout);exit(1);
	/* NOTE : If this happens, CollideBlock[],EntryBlock[] can be turned in 2-D arrays, with first index being the low order T (eg 8) bits of the hashindex(key) value, 
	   expanding the possible size by up to (1 << T) (eg 256) times. CollideMemR[] and EntryMem[] would similarly be 2-D arrays */
      }
      CollideBlockMax = 2 * CollideBlockMax + 1;/* sizes are 2^^n - 1, to make full use of 32 bit unsigned int */
    } while (CollideBlockCnt > CollideBlockMax);

    if(VERB>=2){
      if(HashInsertThreads > 1)
	printf("CollideBlockAlloc:size=%u,CollideBlockCnt=%u,CollideBlockMax=%u->%u,hashtable=%p\n",
	       size,CollideBlockCnt,origCollideBlockMax,CollideBlockMax,this);
      else
	printf("CollideBlockAlloc:size=%u,CollideBlockCnt=%u,CollideBlockMax=%u->%u\n",
	       size,CollideBlockCnt,origCollideBlockMax,CollideBlockMax);
      fflush(stdout);
    }

    assert(CollideBlockMax >= CollideBlockCnt);
    CHashCollide *origCollideBlock = CollideBlock;
    POSIX_MEMALIGN(CollideBlock,64,sizeof(CHashCollide)*CollideBlockMax);
    memcpy(CollideBlock,origCollideBlock,sizeof(CHashCollide)*origCollideBlockCnt);
    free(origCollideBlock);
  }
  return origCollideBlockCnt;
}

/* locate next block of specified size in EntryBlock[] */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline unsigned int CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::EntryBlockAlloc(unsigned int size)
{
  unsigned int origEntryBlockCnt = EntryBlockCnt;
  EntryBlockCnt += size;
  if(EntryBlockCnt > EntryBlockMax){/* double size of EntryBlock[] */
    unsigned int origEntryBlockMax = EntryBlockMax;
    do {
      if(DEBUG && EntryBlockMax >= MASK_LL(32)){
	printf("Cannot double size of EntryBlock from %u using 32 bit unsigned int\n",EntryBlockMax);
	fflush(stdout);exit(1);
	/* NOTE : If this happens, CollideBlock[],EntryBlock[] can be turned in 2-D arrays, with first index being the low order T (eg 8) bits of the hashindex(key) value, 
	   expanding the possible size by up to (1 << T) (eg 256) times. CollideMemR[] and EntryMem[] would similarly be 2-D arrays, and EntryBlockCnt, EntryBlockMax, CollideBlockCnt, CollideBlockMax would be 1-D array */
      }
      EntryBlockMax = EntryBlockMax * 2 + 1;/* sizes are 2^^n - 1, to make full use of 32 bit unsigned int */
    } while(EntryBlockCnt > EntryBlockMax);

    if(VERB>=2){
      if(HashInsertThreads > 1)
	printf("EntryBlockAlloc:size=%u,EntryBlockCnt=%u,EntryBlockMax=%u->%u,hashtable=%p\n",
	       size, EntryBlockCnt,origEntryBlockMax,EntryBlockMax,this);
      else
	printf("EntryBlockAlloc:EntryBlockCnt=%u->%u,EntryBlockMax=%u->%u\n",
	       size, EntryBlockCnt,origEntryBlockMax,EntryBlockMax);
      fflush(stdout);
    }

    assert(EntryBlockMax >= EntryBlockCnt);
    CHashEntry<OFFSET1_T> *origEntryBlock = EntryBlock;
    POSIX_MEMALIGN(EntryBlock,64,sizeof(CHashEntry<OFFSET1_T>) * EntryBlockMax);
    if(DEBUG>=2){
      register PFLOAT z = nan("NaN");
      for(register int i = EntryBlockMax; --i >= 0;){
	register CHashEntry<OFFSET1_T> *phash = &EntryBlock[i];
	for(register int w = HashWin; --w >= 0;)
	  phash->win1[w] = z;
      }
    }
    memcpy(EntryBlock,origEntryBlock,sizeof(CHashEntry<OFFSET1_T>)*origEntryBlockCnt);
    free(origEntryBlock);
  }
  return origEntryBlockCnt;
}

#if USE_PFLOAT==0
/** compute hash key for window X[0..HashWin] */
static inline long long hashkey(FLOAT *X)
{
  int i = HashWin-1;
  int ilen = floor(X[i]*QLenInv+0.5);
  if(DEBUG>=2) assert(ilen <= MaxIntLen);
  long long key = IntLen[ilen];
  for(; --i >= 0;){
    key <<= SIZE_BITS;
    ilen = floor(X[i]*QLenInv+0.5);
    if(DEBUG>=2 && !(ilen <= MaxIntLen)){
      printf("hashkey:i=%d,X[i]=%0.3f,X[i+1]=%0.3f,ilen=%d,MaxIntLen=%d,MaxFragLen=%0.3f\n",
	     i,X[i],X[i+1],ilen,MaxIntLen,MaxFragLen);
      fflush(stdout);
      assert(ilen <= MaxIntLen);
    }
    key |= IntLen[ilen];
  }
  return key;
}

#else // USE_PFLOAT==1

static inline long long hashkey(PFLOAT *X)
{
  int i = HashWin-1;
  int ilen = floor(X[i]*QLenInv+0.5);
  if(DEBUG>=2) assert(0 <= ilen && ilen <= MaxIntLen);
  long long key = IntLen[ilen];
  for(; --i >= 0;){
    key <<= SIZE_BITS;
    ilen = floor(X[i]*QLenInv+0.5);
    if(DEBUG>=2 && !(X[i] > 0.0 && 0 <= ilen && ilen <= MaxIntLen)){
      printf("hashkey:i=%d,X[i]=%0.3f,ilen=%d,MaxIntLen=%d,MaxFragLen=%0.3f\n",
	     i,X[i],ilen,MaxIntLen,MaxFragLen);
      fflush(stdout);
      assert(X[i] > 0.0);
      assert(ilen <= MaxIntLen);
    }
    key |= IntLen[ilen];
  }
  return key;
}

#endif // USE_PFLOAT==1

/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,err2) */
//# pragma optimization_level 2
template<int HASH_SCALE>
inline void window(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int err2)
{
  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++)
      win2[w] = -1.0;

  if(!HASH_SCALE && ((USE_SSE && USE_AVX) || USE_MIC)){
    PFLOAT **delX = (flip2 ? qmap->deltaR[hcolor] : qmap->deltaX[hcolor]);
    PFLOAT *p = &delX[1][site2+1];
    if(err2){
      if(DEBUG>=2) assert(err2 <= HashWin);
      PFLOAT xint;
      for(int w = 0; w < err2-1; w++){
	xint = p[w];
	if(DEBUG>=2) assert(fabs(xint - (X[site2 + w + 1]-X[site2 + w])) <= 1e-3);
	win2[w] = xint;
      }
      xint = delX[2][site2 + err2 + 1];
      if(DEBUG>=2) assert(fabs(xint - (X[site2 + err2 + 1] - X[site2 + err2 - 1])) <= 1e-3);
      win2[err2-1] = xint;
      for(int w = err2+1; w <= HashWin; w++){
	xint = p[w];
	if(DEBUG>=2) assert(fabs(xint - (X[site2 + w + 1] - X[site2 + w])) <= 1e-3);
	win2[w - 1] = xint;
      }
    } else {
      for(int w = 0; w < HashWin; w++){
	PFLOAT xint = p[w];
	if(DEBUG>=2) assert(fabs(xint - (X[site2+w+1]-X[site2+w])) <= 1e-3);
	win2[w] = xint;
      }
    }
  } else {
    if(err2){
      if(DEBUG>=2) assert(err2 <= HashWin);
#pragma vector unaligned
      for(int w = 0; w < err2-1; w++)
	win2[w] = X[site2 + w + 1] - X[site2 + w];

      win2[err2 - 1] = X[site2 + err2 + 1] - X[site2 + err2 - 1];

#pragma vector unaligned
      for(int w = err2+1; w <= HashWin; w++)
	win2[w - 1] = X[site2 + w + 1] - X[site2 + w];

    } else {
#pragma vector unaligned
      for(int w = 0; w < HashWin; w++)
	win2[w] = X[site2+w+1]-X[site2+w];
    }
  }

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,mis2), where X[mis2..mis2+1] are sites that failed to resolve */
void windowM1(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2)
{
  if(DEBUG>=2) assert(0 <= mis2 && mis2 <= HashWin);
  if(DEBUG>=2) assert(HashWin+1 < MAXR_WIN+8);

  //  FLOAT *Z = new FLOAT[MAXR_WIN + 8];
  FLOAT Z[MAXR_WIN+8];
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2;  w++)
    Z[cnt++] = X[site2 + w];

  Z[cnt++] = Yc(X, site2 + mis2 + 1, 1);

#pragma vector unaligned
  for(int w = mis2 + 2; w <= HashWin+1; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2) assert(cnt == HashWin + 1);

  for(int w = 0; w < HashWin; w++)
    win2[w] = Z[w+1] - Z[w];

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      if(!isfinite(win2[w]) || !(win2[w] > 0.0)){
	#pragma omp critical
	{
	  printf("windowM1:mapid2=%d,flip2=%d,site2=%d,mis2=%d:w=%d/%d:Z[w]= %0.6f, Z[w+1]= %0.6f, win2[w]= %0.6f\n",
		 mapid2,flip2,site2,mis2,w,HashWin,Z[w],Z[w+1],win2[w]);
	  for(int t = site2; t <= site2 + HashWin + 1; t++)
	    printf("X[%d] = %0.6f\n",t,X[t]);
	  fflush(stdout);
	  assert(isfinite(win2[w]));
	  assert(win2[w] > 0.0);
	}
      }
    }

  //  delete [] Z;
}

/** compute a vector of interval sizes coresponding to hash windo w(mapid2,flip2,site2,mis2A,mis2B), where X[mis2A..mis2A+1] and X[mis2B..mis2B+1] are sites that failed to resolve */
void windowM2(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2A, int mis2B)
{
  if(DEBUG>=2) assert(0 <= mis2A && mis2A+1 <= mis2B && mis2B <= HashWin+1);
  if(DEBUG>=2) assert(HashWin+2 < MAXR_WIN+8);

  FLOAT Z[MAXR_WIN+8];
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2A;  w++)
    Z[cnt++] = X[site2 + w];

  if(mis2B == mis2A + 1){/* consecutive misresolved sites mis2A,mis2B */
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 2);
    goto LafterB;
  }
  /* mis2B >= mis2A + 2 */
  Z[cnt++] = Yc(X, site2 + mis2A + 1, 1);
    
  //LafterA:
#pragma vector unaligned
  for(int w = mis2A + 2; w < mis2B; w++)
    Z[cnt++] = X[site2 + w];
    
  Z[cnt++] = Yc(X, site2 + mis2B + 1, 1);

 LafterB:
#pragma vector unaligned
  for(int w = mis2B + 2; w <= HashWin+2; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2) assert(cnt == HashWin + 1);

#pragma vector unaligned
  for(int w = 0; w < HashWin; w++)
    win2[w] = Z[w+1] - Z[w];

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,mis2A,mis2B,err2), where X[mis2A..mis2A+1], X[mis2B..mis2B+1] are sites that failed to resolve and 
   site err2 (of the deres'd sites) is missing */
void windowM2E1(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2A, int mis2B, int err2)
{
  if(DEBUG>=2) assert(0 <= mis2A && mis2A+1 <= mis2B && mis2B <= HashWin+2);
  if(DEBUG>=2) assert(0 < err2 && err2 <= HashWin);

  FLOAT Z[MAXR_WIN+8];

  /* first combine misresolved sites (mis2A,mis2A+1) abd (mis2B,mis2B+1) */
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2A;  w++)
    Z[cnt++] = X[site2 + w];

  if(mis2B == mis2A + 1)/* consecutive misresolved sites */
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 2);
  else {/* mis2B >= mis2A + 2 */
    Z[cnt++] = Yc(X, site2 + mis2A + 1, 1);
    
#pragma vector unaligned
    for(int w = mis2A + 2; w < mis2B; w++)
      Z[cnt++] = X[site2 + w];
    
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 1);
  }

#pragma vector unaligned
  for(int w = mis2B + 2; w <= HashWin+3; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2) assert(cnt == HashWin + 2);

  /* next remove site err2 */
  cnt = err2;
  for(int w = err2+1; w <= HashWin+1; w++)
    Z[cnt++] = Z[w];

  if(DEBUG>=2) assert(cnt == HashWin + 1);

  for(int w = 0; w < HashWin; w++)
    win2[w] = Z[w+1] - Z[w];

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,mis2A,mis2B,err2A,err2B), where X[mis2A..mis2A+1], X[mis2B..mis2B+1] are sites that failed to resolve and 
   sites err2A and err2B (of the deres'd sites) are missing */
void windowM2E2(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2A, int mis2B, int err2A, int err2B)
{
  if(DEBUG>=2) assert(0 <= mis2A && mis2A+1 <= mis2B && mis2B <= HashWin+3);
  if(DEBUG>=2) assert(0 < err2A && err2A < err2B && err2B <= HashWin+1);

  FLOAT Z[MAXR_WIN+8];

  /* first combine misresolved sites (mis2A,mis2A+1) abd (mis2B,mis2B+1) */
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2A;  w++)
    Z[cnt++] = X[site2 + w];

  if(mis2B == mis2A + 1)/* consecutive misresolved sites */
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 2);
  else {/* mis2B >= mis2A + 2 */
    Z[cnt++] = Yc(X, site2 + mis2A + 1, 1);
    
#pragma vector unaligned
    for(int w = mis2A + 2; w < mis2B; w++)
      Z[cnt++] = X[site2 + w];
    
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 1);
  }

#pragma vector unaligned
  for(int w = mis2B + 2; w <= HashWin+4; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2) assert(cnt == HashWin + 3);

  /* next remove site's err2A & err2B */
  cnt = err2A;
  for(int w = err2A+1; w < err2B; w++)
    Z[cnt++] = Z[w];
  for(int w = err2B+1; w <= HashWin+2; w++)
    Z[cnt++] = Z[w];

  if(DEBUG>=2) assert(cnt == HashWin + 1);

  for(int w = 0; w < HashWin; w++)
    win2[w] = Z[w+1] - Z[w];

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,mis2,err2), where X[mis2..mis2+1] are sites that failed to resolve and err2 != mis2 is the missing site */
void windowM1E1(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2, int err2)
{
  if(DEBUG>=2) assert(0 <= mis2 && mis2 <= HashWin+1);
  if(DEBUG>=2) assert(0 < err2 && err2 <= HashWin);

  FLOAT Z[MAXR_WIN+8];

  /* first combine misresolved sites mis2,mis2+1 */
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2;  w++)
    Z[cnt++] = X[site2 + w];

  Z[cnt++] = Yc(X, site2 + mis2 + 1, 1);

#pragma vector unaligned
  for(int w = mis2 + 2; w <= HashWin+2; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2) assert(cnt == HashWin + 2);
  
  /* next remove site err2 */
  cnt = err2;
  for(int w = err2+1; w <= HashWin+1; w++)
    Z[cnt++] = Z[w];

  if(DEBUG>=2) assert(cnt == HashWin + 1);

  for(int w = 0; w < HashWin; w++){
    win2[w] = Z[w+1] - Z[w];
    if(DEBUG>=2) assert(win2[w] > 0.0);
  }

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}


/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,mis2A,mis2B,mis2C), where X[mis2A..mis2A+1], X[mis2B..mis2B+1] and X[mis2C..misc2C+1] are sites that failed to resolve */
void windowM3(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2A, int mis2B, int mis2C)
{
  if(DEBUG>=2) assert(0 <= mis2A && mis2A+1 <= mis2B && mis2B+1 <= mis2C && mis2C <= HashWin+2);

  FLOAT Z[MAXR_WIN+8];
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2A;  w++)
    Z[cnt++] = X[site2 + w];

  if(mis2B == mis2A + 1){/* 2 consecutive misresolved sites mis2A, mis2B */
    if(mis2C == mis2B + 1){ /* 3 consecutive misresolved sites mis2A, mis2B, mis2C */
      Z[cnt++] = Yc(X, site2 + mis2C + 1, 3);
      goto LafterC;
    }
    /* mis2C >= mis2B + 2 */
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 2);
    goto LafterB;
  } 
  /* mis2B >= mis2A + 2 */
  Z[cnt++] = Yc(X, site2 + mis2A + 1, 1);
    
  //LafterA:
#pragma vector unaligned
  for(int w = mis2A + 2; w < mis2B; w++)
    Z[cnt++] = X[site2 + w];
    
  if(mis2C == mis2B + 1){/* 2 consecutive misresolved sites mis2B, mis2C */
    Z[cnt++] = Yc(X, site2 + mis2C + 1, 2);
    goto LafterC;
  }
  /* mis2C >= mis2B + 2 */
  Z[cnt++] = Yc(X, site2 + mis2B + 1, 1);

 LafterB:
#pragma vector unaligned
  for(int w = mis2B + 2; w < mis2C; w++)
    Z[cnt++] = X[site2 + w];

  Z[cnt++] = Yc(X, site2 + mis2C + 1, 1);

 LafterC:
#pragma vector unaligned
  for(int w = mis2C + 2; w <= HashWin+3; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2) assert(cnt == HashWin + 1);
  for(int w = 0; w < HashWin; w++)
    win2[w] = Z[w+1] - Z[w];

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,mis2A,mis2B,mis2C,mis2D,mis2E), 
    where X[mis2A..mis2A+1], X[mis2B..mis2B+1],X[mis2C..mis2C+1], X[mis2D..mis2D+1] and X[mis2E..mis2E+1] are sites that failed to resolve */
void windowM5(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2A, int mis2B, int mis2C, int mis2D, int mis2E)
{
  if(DEBUG>=2) assert(0 <= mis2A && mis2A+1 <= mis2B && mis2B+1 <= mis2C && mis2C+1 <= mis2D && mis2D+1 <= mis2E && mis2E <= HashWin+4);

  FLOAT Z[MAXR_WIN+8];
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2A;  w++)
    Z[cnt++] = X[site2 + w];

  if(mis2B == mis2A + 1){/* 2 consecutive misresolved sites mis2A, mis2B */
    if(mis2C == mis2B + 1){ /* 3 consecutive misresolved sites mis2A, mis2B, mis2C */
      if(mis2D == mis2C + 1){/* 4 consecutive misresolved sites mis2A,mis2B,mis2C,mis2D */
	if(mis2E == mis2D + 1){/* 5 consecutive misresolved sites mis2A,mis2B,mis2C,mis2D,mis2E */
	  Z[cnt++] = Yc(X,site2 + mis2E + 1, 5);
	  goto LafterE;
	}
	/* mis2E >= mis2D + 2 */
	Z[cnt++] = Yc(X,site2 + mis2D + 1, 4);
	goto LafterD;
      }
      /* mis2D >= mis2C + 2 */
      Z[cnt++] = Yc(X, site2 + mis2C + 1, 3);
      goto LafterC;
    }
    /* mis2C >= mis2B + 2 */
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 2);
    goto LafterB;
  }
  /* mis2B >= mis2A + 2 */
  Z[cnt++] = Yc(X, site2 + mis2A + 1, 1);
    
  //  LafterA:
#pragma vector unaligned
  for(int w = mis2A + 2; w < mis2B; w++)
    Z[cnt++] = X[site2 + w];
    
  if(mis2C == mis2B + 1){/* 2 consecutive misresolved sites mis2B, mis2C */
    if(mis2D == mis2C + 1){/* 3 consecutive misresolved sites mis2B,mis2C,mis2D */
      if(mis2E == mis2D + 1){/* 4 consecutive misresolved sites mis2B,mis2C,mis2D,mis2E */
	Z[cnt++] = Yc(X, site2 + mis2E + 1, 4);
	goto LafterE;
      }
      /* mis2E >= mis2D + 2 */
      Z[cnt++] = Yc(X, site2 + mis2D + 1, 3);
      goto LafterD;
    }

    /* mis2D >= mis2C + 2 */
    Z[cnt++] = Yc(X, site2 + mis2C + 1, 2);
    goto LafterC;
  }
  /* mis2C >= mis2B + 2 */
  Z[cnt++] = Yc(X, site2 + mis2B + 1, 1);

 LafterB:
#pragma vector unaligned
  for(int w = mis2B + 2; w < mis2C; w++)
    Z[cnt++] = X[site2 + w];

  if(mis2D == mis2C + 1){/* 2 consecutive misresolved sites mis2C,mis2D */
    if(mis2E == mis2D + 1){/* 3 consecutive misresolved sites mis2C,mis2D,mis2E */
      Z[cnt++] = Yc(X, site2 + mis2E + 1, 3);
      goto LafterE;
    }
    /* mis2E >= mis2D + 2 */
    Z[cnt++] = Yc(X, site2 + mis2D + 1, 2);
    goto LafterD;
  }
  /* mis2D >= mis2C + 2 */
  Z[cnt++] = Yc(X, site2 + mis2C + 1, 1);

 LafterC:
#pragma vector unaligned
  for(int w = mis2C + 2; w < mis2D; w++)
    Z[cnt++] = X[site2 + w];
  if(mis2E == mis2D + 1){/* 2 consecutive misresolved sites mis2D,mis2E */
    Z[cnt++] = Yc(X, site2 + mis2E + 1, 2);
    goto LafterE;
  }
  /* mis2E >= mis2D + 2 */
  Z[cnt++] = Yc(X, site2 + mis2D + 1, 1);  

 LafterD:
#pragma vector unaligned
  for(int w = mis2D + 2; w < mis2E; w++)
    Z[cnt++] = X[site2 + w];
  
  Z[cnt++] = Yc(X, site2 + mis2E + 1, 1);

 LafterE:
  for(int w = mis2E + 2; w <= HashWin+5; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2){
    for(int w = 1; w < cnt; w++){
       if(DEBUG && !(Z[w] > Z[w-1])){
         #pragma omp critical
         {
	   printf("windowM5:mapid2=%d,flip2=%d,site2=%d,mis2=%d,%d,%d,%d,%d,HashWin=%u:cnt=%d\n",mapid2,flip2,site2,mis2A,mis2B,mis2C,mis2D,mis2E,HashWin,cnt);
	   for(int k = 0; k <= HashWin+5; k++)
	     printf("w=%d:X[site2+w]=%0.6f\n",k,X[site2+k]);
	   for(int c = 0; c < cnt; c++)
	     printf("c=%d:Z[c]=%0.6f\n",c,Z[c]);
	   fflush(stdout);
	   assert(Z[w] > Z[w-1]);
        }
      }
    }
  }

  if(DEBUG>=2 && !(cnt == HashWin+1)){
    #pragma omp critical
    {
      printf("windowM5:mapid2=%d,flip2=%d,site2=%d,mis2=%d,%d,%d,%d,%d,HashWin=%u:cnt=%d\n",mapid2,flip2,site2,mis2A,mis2B,mis2C,mis2D,mis2E,HashWin,cnt);
      for(int k = 0; k <= HashWin+5; k++)
	printf("w=%d:X[site2+w]=%0.6f\n",k,X[site2+k]);
      for(int c = 0; c < cnt; c++)
	printf("c=%d:Z[c]=%0.6f\n",c,Z[c]);
      fflush(stdout);
      assert(cnt == HashWin + 1);
    }
  }
  for(int w = 0; w < HashWin; w++){
    win2[w] = Z[w+1] - Z[w];
    if(DEBUG>=2) assert(0.0 < win2[w]);
  }

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,mis2A,mis2B,mis2C,mis2D), 
    where X[mis2A..mis2A+1], X[mis2B..mis2B+1],X[mis2C..mis2C+1] and X[mis2D..misc2D+1] are sites that failed to resolve */
void windowM4(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2A, int mis2B, int mis2C, int mis2D)
{
  if(DEBUG>=2) assert(0 <= mis2A && mis2A+1 <= mis2B && mis2B+1 <= mis2C && mis2C+1 <= mis2D && mis2D <= HashWin+3);

  FLOAT Z[MAXR_WIN+8];
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2A;  w++)
    Z[cnt++] = X[site2 + w];

  if(mis2B == mis2A + 1){/* 2 consecutive misresolved sites mis2A, mis2B */
    if(mis2C == mis2B + 1){ /* 3 consecutive misresolved sites mis2A, mis2B, mis2C */
      if(mis2D == mis2C + 1){/* 4 consecutive misresolved sites mis2A,mis2B,mis2C,mis2D */
	Z[cnt++] = Yc(X,site2 + mis2D + 1, 4);
	goto LafterD;
      }
      /* mis2D >= mis2C + 2 */
      Z[cnt++] = Yc(X, site2 + mis2C + 1, 3);
      goto LafterC;
    }
    /* mis2C >= mis2B + 2 */
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 2);
    goto LafterB;
  }
  /* mis2B >= mis2A + 2 */
  Z[cnt++] = Yc(X, site2 + mis2A + 1, 1);
    
  //  LafterA:
#pragma vector unaligned
  for(int w = mis2A + 2; w < mis2B; w++)
    Z[cnt++] = X[site2 + w];
    
  if(mis2C == mis2B + 1){/* 2 consecutive misresolved sites mis2B, mis2C */
    if(mis2D == mis2C + 1){/* 3 consecutive misresolved sites mis2B,mis2C,mis2D */
      Z[cnt++] = Yc(X, site2 + mis2D + 1, 3);
      goto LafterD;
    }

    /* mis2D >= mis2C + 1 */
    Z[cnt++] = Yc(X, site2 + mis2C + 1, 2);
    goto LafterC;
  }
  /* mis2C >= mis2B + 2 */
  Z[cnt++] = Yc(X, site2 + mis2B + 1, 1);

 LafterB:
#pragma vector unaligned
  for(int w = mis2B + 2; w < mis2C; w++)
    Z[cnt++] = X[site2 + w];

  if(mis2D == mis2C + 1){/* 2 consecutive misresolved sites mis2C,mis2D */
    Z[cnt++] = Yc(X, site2 + mis2D + 1, 2);
    goto LafterD;
  }
  /* mis2D >= mis2C + 2 */
  Z[cnt++] = Yc(X, site2 + mis2C + 1, 1);

 LafterC:
#pragma vector unaligned
  for(int w = mis2C + 2; w < mis2D; w++)
    Z[cnt++] = X[site2 + w];

  Z[cnt++] = Yc(X, site2 + mis2D + 1, 1);  

 LafterD:
#pragma vector unaligned
  for(int w = mis2D + 2; w <= HashWin+4; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2){
    for(int w = 1; w < cnt; w++){
       if(DEBUG && !(Z[w] > Z[w-1])){
         #pragma omp critical
         {
	   printf("windowM4:mapid2=%d,flip2=%d,site2=%d,mis2=%d,%d,%d,%d,HashWin=%u:cnt=%d\n",mapid2,flip2,site2,mis2A,mis2B,mis2C,mis2D,HashWin,cnt);
	   for(int k = 0; k <= HashWin+4; k++)
	     printf("w=%d:X[site2+w]=%0.6f\n",k,X[site2+k]);
	   for(int c = 0; c < cnt; c++)
	     printf("c=%d:Z[c]=%0.6f\n",c,Z[c]);
	   fflush(stdout);
	   assert(Z[w] > Z[w-1]);
        }
      }
    }
  }

  if(DEBUG>=2) assert(cnt == HashWin + 1);
  for(int w = 0; w < HashWin; w++){
    win2[w] = Z[w+1] - Z[w];
    if(DEBUG>=2) assert(0.0 < win2[w]);
  }

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,mis2A,mis2B,mis2C,mis2D,err2), 
    where X[mis2A..mis2A+1], X[mis2B..mis2B+1],X[mis2C..mis2C+1] and X[mis2D..misc2D+1] are sites that failed to resolve and site err2 (of the deres'd sites) is missing */
 void windowM4E1(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2A, int mis2B, int mis2C, int mis2D, int err2)
{
  if(DEBUG>=2) assert(0 <= mis2A && mis2A+1 <= mis2B && mis2B+1 <= mis2C && mis2C+1 <= mis2D && mis2D <= HashWin+4);
  if(DEBUG>=2) assert(0 < err2 && err2 <= HashWin);

  FLOAT Z[MAXR_WIN+8];

  /* first combine misresolved sites mis2A, mis2B, mis2C and mis2D (by combining them with the next site) */
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2A;  w++)
    Z[cnt++] = X[site2 + w];

  if(mis2B == mis2A + 1){/* 2 consecutive misresolved sites mis2A, mis2B */
    if(mis2C == mis2B + 1){ /* 3 consecutive misresolved sites mis2A, mis2B, mis2C */
      if(mis2D == mis2C + 1){/* 4 consecutive misresolved sites mis2A,mis2B,mis2C,mis2D */
	Z[cnt++] = Yc(X,site2 + mis2D + 1, 4);
	goto LafterD;
      }
      /* mis2D >= mis2C + 2 */
      Z[cnt++] = Yc(X, site2 + mis2C + 1, 3);
      goto LafterC;
    }
    /* mis2C >= mis2B + 2 */
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 2);
    goto LafterB;
  }
  /* mis2B >= mis2A + 2 */
  Z[cnt++] = Yc(X, site2 + mis2A + 1, 1);
    
  //  LafterA:
#pragma vector unaligned
  for(int w = mis2A + 2; w < mis2B; w++)
    Z[cnt++] = X[site2 + w];
    
  if(mis2C == mis2B + 1){/* 2 consecutive misresolved sites mis2B, mis2C */
    if(mis2D == mis2C + 1){/* 3 consecutive misresolved sites mis2B,mis2C,mis2D */
      Z[cnt++] = Yc(X, site2 + mis2D + 1, 3);
      goto LafterD;
    }

    /* mis2D >= mis2C + 1 */
    Z[cnt++] = Yc(X, site2 + mis2C + 1, 2);
    goto LafterC;
  }
  /* mis2C >= mis2B + 2 */
  Z[cnt++] = Yc(X, site2 + mis2B + 1, 1);

 LafterB:
#pragma vector unaligned
  for(int w = mis2B + 2; w < mis2C; w++)
    Z[cnt++] = X[site2 + w];

  if(mis2D == mis2C + 1){/* 2 consecutive misresolved sites mis2C,mis2D */
    Z[cnt++] = Yc(X, site2 + mis2D + 1, 2);
    goto LafterD;
  }
  /* mis2D >= mis2C + 2 */
  Z[cnt++] = Yc(X, site2 + mis2C + 1, 1);

 LafterC:
#pragma vector unaligned
  for(int w = mis2C + 2; w < mis2D; w++)
    Z[cnt++] = X[site2 + w];

  Z[cnt++] = Yc(X, site2 + mis2D + 1, 1);  

 LafterD:
#pragma vector unaligned
  for(int w = mis2D + 2; w <= HashWin+5; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2 && !(cnt == HashWin + 2)){
    #pragma omp critical
    {
      printf("windowM4E1:mapid2=%d,flip2=%d,site2=%d,mis2=%u,%u,%u,%u,err2=%u,HashWin=%u:cnt=%d\n",mapid2,flip2,site2,mis2A,mis2B,mis2C,mis2D,err2,HashWin,cnt);
      for(int w = 0; w <= HashWin+5; w++)
	printf("w=%d:X[site2+w]=%0.6f\n",w,X[site2+w]);
      for(int c = 0; c < cnt; c++)
	printf("c=%d:Z[c]=%0.6f\n",c,Z[c]);
      fflush(stdout);
      assert(cnt == HashWin + 2);
    }
  }
  if(DEBUG>=2){
    for(int w = 1; w < cnt; w++){
       if(DEBUG && !(Z[w] > Z[w-1])){
         #pragma omp critical
         {
	   printf("windowM4E1:mapid2=%d,flip2=%d,site2=%d,mis2=%u,%u,%u,%u,err2=%u,HashWin=%u:cnt=%d\n",mapid2,flip2,site2,mis2A,mis2B,mis2C,mis2D,err2,HashWin,cnt);
	   for(int t = 0; t <= HashWin+5; t++)
	     printf("w=%d:X[site2+w]=%0.6f\n",t,X[site2+t]);
	   for(int c = 0; c < cnt; c++)
	     printf("c=%d:Z[c]=%0.6f\n",c,Z[c]);
	   fflush(stdout);
	   assert(Z[w] > Z[w-1]);
        }
      }
    }
  }

  /* next remove site err2 */
  cnt = err2;
  for(int w = err2+1; w <= HashWin+1; w++)
    Z[cnt++] = Z[w];

  if(DEBUG>=2) assert(cnt == HashWin + 1);
  for(int w = 0; w < HashWin; w++){
    win2[w] = Z[w+1] - Z[w];
    if(DEBUG>=2 && !(win2[w] > 0.0)){
      #pragma omp critical
      {
	printf("windowM4E1:mapid2=%d,flip2=%d,site2=%d,mis2=%u,%u,%u,%u,err2=%u,HashWin=%u:cnt=%d\n",mapid2,flip2,site2,mis2A,mis2B,mis2C,mis2D,err2,HashWin,cnt);
	for(int t = 0; t <= HashWin+5; t++)
	  printf("w=%d:X[site2+w]=%0.6f\n",t,X[site2+t]);
	for(int c = 0; c < cnt; c++)
	  printf("c=%d:Z[c]=%0.6f\n",c,Z[c]);
	for(int t = 0; t <= w; t++)
	  printf("t=%d:win2[t]=%0.6f\n",t,win2[t]);
	fflush(stdout);
	assert(win2[w] > 0.0);
      }
    }
  }

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window with 1 misresolved site AND 2 missing site errors */
inline void windowM1E2(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2, int err2A, int err2B)
{
  if(DEBUG>=2) assert(0 <= mis2 && mis2 <= HashWin+2);
  if(DEBUG>=2) assert(0 < err2A  && err2A < err2B && err2B <= HashWin+1);

  FLOAT Z[MAXR_WIN+8];
  //  FLOAT *Z = new FLOAT[MAXR_WIN+8];

  /* first combine misresolved sites mis2,mis2+1 */
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2;  w++)
    Z[cnt++] = X[site2 + w];

  Z[cnt++] = Yc(X, site2 + mis2 + 1, 1);

#pragma vector unaligned
  for(int w = mis2 + 2; w <= HashWin+3; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2) assert(cnt == HashWin + 3);
  
  /* next remove sites err2A & err2B */
  cnt = err2A;
  for(int w = err2A+1; w < err2B; w++)
    Z[cnt++] = Z[w];
  for(int w = err2B+1; w <= HashWin+2; w++)
    Z[cnt++] = Z[w];

  if(DEBUG>=2) assert(cnt == HashWin + 1);

  for(int w = 0; w < HashWin; w++)
    win2[w] = Z[w+1] - Z[w];

  //  delete [] Z;
  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window with 2 errors */
template<int HASH_SCALE>
inline void window2(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int err2A, int err2B)
{
  if(DEBUG>=2) assert(0 < err2A  && err2A < err2B && err2B < HashWin+2);
  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++)
      win2[w] = -1.0;

  if(!HASH_SCALE && ((USE_SSE && USE_AVX) || USE_MIC)){
    PFLOAT **delX = (flip2 ? qmap->deltaR[hcolor] : qmap->deltaX[hcolor]);
    PFLOAT *p = &delX[1][site2+1];
    PFLOAT xint;

    for(int w = 0; w < err2A - 1; w++){
      xint = p[w];
      if(DEBUG>=2) assert(fabs(xint - (X[site2 + w + 1]-X[site2 + w])) <= 1e-3);
      win2[w] = xint;
    }
    if(err2B == err2A+1){
      xint = delX[3][site2 + err2B + 1];
      if(DEBUG>=2) assert(fabs(xint - (X[site2 + err2B + 1] - X[site2 + err2A - 1])) <= 1e-3);
      win2[err2A - 1] = xint;
      for(int w = err2B + 1; w <= HashWin+1; w++){
	xint = p[w];
	if(DEBUG>=2) assert(fabs(xint - (X[site2 + w + 1] - X[site2 + w])) <= 1e-3);
	win2[w - 2] = xint;
      }
    } else {
      xint = delX[2][site2 + err2A + 1];
      if(DEBUG>=2) assert(fabs(xint - (X[site2 + err2A + 1] - X[site2 + err2A - 1])) <= 1e-3);
      win2[err2A - 1] = xint;
      for(int w = err2A + 1; w < err2B - 1; w++){
	xint = p[w];
	if(DEBUG>=2) assert(fabs(xint - (X[site2 + w + 1] - X[site2 + w])) <= 1e-3);
	win2[w - 1] = xint;
      }
      xint = delX[2][site2 + err2B + 1];
      if(DEBUG>=2 && !(fabs(xint - (X[site2 + err2B + 1] - X[site2 + err2B - 1])) <= 1e-3)){
	printf("site2=%d,err2B=%d:xint = %0.10f, X[site2+err2B+1]-X[site2 + err2B - 1]=%0.10f, err=%0.10f\n",
	       site2,err2B,xint, X[site2 + err2B + 1] - X[site2+err2B - 1],xint - (X[site2 + err2B + 1] - X[site2 + err2B - 1]));
	fflush(stdout);
	assert(fabs(xint - (X[site2 + err2B + 1] - X[site2 + err2B - 1])) <= 1e-3);
      }
      win2[err2B - 2] = xint;
      for(int w = err2B + 1; w <= HashWin+1; w++){
	xint = p[w];
	if(DEBUG>=2) assert(fabs(xint - (X[site2 + w + 1] - X[site2 + w])) <= 1e-3);
	win2[w - 2] = xint;
      }
    }
  } else {
    for(int w = 0; w < err2A - 1; w++)
      win2[w] = X[site2 + w + 1] - X[site2 + w];
    if(err2B == err2A+1){
      win2[err2A - 1] = X[site2 + err2B + 1] - X[site2 + err2A - 1];
      for(int w = err2B + 1; w <= HashWin+1; w++)
	win2[w - 2] = X[site2 + w + 1] - X[site2 + w];
    } else {
      win2[err2A - 1] = X[site2 + err2A + 1] - X[site2 + err2A - 1];
      for(int w = err2A + 1; w < err2B - 1; w++)
	win2[w - 1] = X[site2 + w + 1] - X[site2 + w];
      win2[err2B - 2] = X[site2 + err2B + 1] - X[site2 + err2B - 1];	
      for(int w = err2B + 1; w <= HashWin+1; w++)
	win2[w - 2] = X[site2 + w + 1] - X[site2 + w];
    }
  }

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,mis2A,mis2B,mis2C,err2), 
    where X[mis2A..mis2A+1], X[mis2B..mis2B+1] and X[mis2C..misc2C+1] are sites that failed to resolve and X[mis2] is missing */
void windowM3E1(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2A, int mis2B, int mis2C, int err2)
{
  if(DEBUG>=2) assert(0 <= mis2A && mis2A+1 <= mis2B && mis2B+1 <= mis2C && mis2C <= HashWin+3);
  if(DEBUG>=2) assert(0 < err2 && err2 <= HashWin);

  FLOAT Z[MAXR_WIN+8];

  /* first combine misresolved sites mis2A, mis2B and mis2C (by combining them with the next site) */
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2A;  w++)
    Z[cnt++] = X[site2 + w];

  if(mis2B == mis2A + 1){/* 2 consecutive misresolved sites mis2A, mis2B */
    if(mis2C == mis2B + 1){ /* 3 consecutive misresolved sites mis2A, mis2B, mis2C */
      Z[cnt++] = Yc(X, site2 + mis2C + 1, 3);
      goto LafterC;
    }

    /* mis2C >= mis2B + 2 */
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 2);
    goto LafterB;
  }
  /* mis2B >= mis2A + 2 */
  Z[cnt++] = Yc(X, site2 + mis2A + 1, 1);
    
  //LafterA:
#pragma vector unaligned
  for(int w = mis2A + 2; w < mis2B; w++)
    Z[cnt++] = X[site2 + w];
    
  if(mis2C == mis2B + 1){/* 2 consecutive misresolved sites mis2B, mis2C */
    Z[cnt++] = Yc(X, site2 + mis2C + 1, 2);
    goto LafterC;
  }
  /* mis2C >= mis2B + 2 */
  Z[cnt++] = Yc(X, site2 + mis2B + 1, 1);

 LafterB:
#pragma vector unaligned
  for(int w = mis2B + 2; w < mis2C; w++)
    Z[cnt++] = X[site2 + w];

  Z[cnt++] = Yc(X, site2 + mis2C + 1, 1);

 LafterC:
#pragma vector unaligned
  for(int w = mis2C + 2; w <= HashWin+4; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2) assert(cnt == HashWin + 2);

  /* next remove site err2 */
  cnt = err2;
  for(int w = err2+1; w <= HashWin+1; w++)
    Z[cnt++] = Z[w];

  if(DEBUG>=2) assert(cnt == HashWin + 1);
  for(int w = 0; w < HashWin; w++)
    win2[w] = Z[w+1] - Z[w];

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

/** compute a vector of interval sizes coresponding to hash window (mapid2,flip2,site2,mis2A,mis2B,mis2C,err2A,err2B), 
    where X[mis2A..mis2A+1], X[mis2B..mis2B+1] and X[mis2C..misc2C+1] are sites that failed to resolve and sites err2A and err2B (of the deres'd sites) are missing */
void windowM3E2(PFLOAT *win2, FLOAT *X, Cmap *qmap, int mapid2, int flip2, int site2, int mis2A, int mis2B, int mis2C, int err2A, int err2B)
{
  if(DEBUG>=2) assert(0 <= mis2A && mis2A+1 <= mis2B && mis2B+1 <= mis2C && mis2C <= HashWin+4);
  if(DEBUG>=2) assert(0 < err2A && err2A < err2B && err2B <= HashWin+1);

  FLOAT Z[MAXR_WIN+8];

  /* first combine misresolved sites mis2A, mis2B and mis2C (by combining them with the next site) */
  int cnt = 0;

#pragma vector unaligned
  for(int w = 0; w < mis2A;  w++)
    Z[cnt++] = X[site2 + w];

  if(mis2B == mis2A + 1){/* 2 consecutive misresolved sites mis2A, mis2B */
    if(mis2C == mis2B + 1){ /* 3 consecutive misresolved sites mis2A, mis2B, mis2C */
      Z[cnt++] = Yc(X, site2 + mis2C + 1, 3);
      goto LafterC;
    }
    /* mis2C >= mis2B + 2 */
    Z[cnt++] = Yc(X, site2 + mis2B + 1, 2);
    goto LafterB;
  }
  /* mis2B >= mis2A + 2 */
  Z[cnt++] = Yc(X, site2 + mis2A + 1, 1);
    
  //LafterA:
#pragma vector unaligned
  for(int w = mis2A + 2; w < mis2B; w++)
    Z[cnt++] = X[site2 + w];
    
  if(mis2C == mis2B + 1){/* 2 consecutive misresolved sites mis2B, mis2C */
    Z[cnt++] = Yc(X, site2 + mis2C + 1, 2);
    goto LafterC;
  }
  /* mis2C >= mis2B + 2 */
  Z[cnt++] = Yc(X, site2 + mis2B + 1, 1);

 LafterB:
#pragma vector unaligned
  for(int w = mis2B + 2; w < mis2C; w++)
    Z[cnt++] = X[site2 + w];

  Z[cnt++] = Yc(X, site2 + mis2C + 1, 1);

 LafterC:
#pragma vector unaligned
  for(int w = mis2C + 2; w <= HashWin+5; w++)
    Z[cnt++] = X[site2 + w];

  if(DEBUG>=2) assert(cnt == HashWin + 3);

  /* next remove site's err2A & err2B */
  cnt = err2A;
  for(int w = err2A+1; w < err2B; w++)
    Z[cnt++] = Z[w];
  for(int w = err2B+1; w <= HashWin+2; w++)
    Z[cnt++] = Z[w];

  if(DEBUG>=2) assert(cnt == HashWin + 1);
  for(int w = 0; w < HashWin; w++)
    win2[w] = Z[w+1] - Z[w];

  if(DEBUG>=2)
    for(int w = 0; w < HashWin; w++){
      assert(isfinite(win2[w]));
      assert(win2[w] > 0.0);
    }
}

static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;

/** Top level call to populate hash table for all maps.
*/
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::hashinsert(Cmap **rmap, int numrmaps)
{
  /* insert maps into HashTable */

  if(DEBUG) assert(HashInsertThreads >= 1);
  CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> **hashtables = 0;
  try {
    hashtables = new CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *[HashInsertThreads];
  } catch (exception& e){
    cout << e.what() << endl;
    printf("hashinsert(numrmaps=%d): new CHashTable*[HashInsertThreads=%d] failed (%lu bytes)\n",
	   numrmaps, HashInsertThreads, HashInsertThreads * sizeof(CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *));
    fflush(stdout);
    assert(0);
  }
  
  hashtables[0] = this;

  #pragma omp parallel num_threads(HashInsertThreads) if(HashInsertThreads > 1)
  {
    int tid = 0;
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #endif

    if(DEBUG) assert(tid >= 0 && tid < HashInsertThreads);
    
    if(tid > 0) {
      try {
        hashtables[tid] = new CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>(rmap,numrmaps,qmap,numqmaps,this);/* One Hashtable per thread */
      } catch (exception& e){
        cout << e.what() << endl;
	printf("hashinsert(numrmaps=%d): new CHashTable() failed\n", numrmaps);
	fflush(stdout);
	assert(0);
      }   
    }
    if(VERB>=2 && HashInsertThreads > 1){
      #pragma omp barrier

      if(tid==0){
	printf("%d additional HashTables initialized: CPU time=%0.6f, wall time=%0.6f secs\n",HashInsertThreads-1, mtime(), wtime());
	fflush(stdout);
      }
    }

    #pragma omp for schedule(dynamic,8)
    for(int i = 0; i < numrmaps; i++){
      if(VERB && i > 0 && !(i % 100000)){
        #pragma omp critical
        {
          printf("Inserted %d/%d maps into hashtable: wall time= %0.6f secs\n", i, numrmaps, wtime());
	  fflush(stdout);
        }
      }
      hashtables[tid]->hash_insert(roffset + i,rmap[roffset + i]);
    }
  }
  
  if(VERB){
    printf("Hash insertion of %d maps complete: CPU time=%0.6f, wall time=%0.6f secs\n",numrmaps,mtime(),wtime());
    //    printf("Sleeping for 60 seconds\n");
    fflush(stdout);
    //    sleep(60);
  }

  //  if(HashInsertThreads <= 1)
  //    hashsort();
  //  else
  hashmerge(hashtables,HashInsertThreads); /* merge, compact and sort hashtable(s) */

  if(VERB>=2){
    printf("Hash merge/compaction/sort done : deleting %d hashtables: CPU time=%0.6f, wall time=%0.6f secs\n", HashInsertThreads-1, mtime(),wtime());
    fflush(stdout);
  }

  /* free up memory from additional hashtables */
  for(int i = 1; i < HashInsertThreads; i++)
    delete hashtables[i];
  delete [] hashtables;

  if(VERB){
    getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
    printf("Hash merge/compaction/sort/free complete, size=%lu bytes (VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb): CPU time=%0.6f, wall time=%0.6f secs\n",
	   HashSize(), VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(),wtime());
    fflush(stdout);
  }

  if(VERB>=3 && HashWin==5){/* print out hashtable contents */
    //    int Tcnt = 0;
    unsigned int Scnt = 0, Lcnt = 0;
    for(int hindex = 0; hindex <= (int)MASK(HASH_BITS); hindex++){
      register CHashIndexR *p = &indexR[hindex];
      register unsigned int cnt = p->cnt;
      if(!cnt)
	continue;
      Scnt += cnt;
      Lcnt += max(COLLIDELINE,cnt)-COLLIDELINE;
      if((hindex%1024)) // output only every 1024'th entry
	continue;
      //      if(++Tcnt > 100 )/* limit to first 100 non-empty entries */
      //	break;
      printf("hindex=%d:cnt=%u, Overflow=%u: Cumulative: cnt=%u, Overflow=%u\n",hindex,cnt,max(COLLIDELINE,cnt)-COLLIDELINE,Scnt,Lcnt);
      //      continue;

      unsigned int imax = min(cnt,COLLIDELINE);
      for(unsigned int i = 0; i < imax; i++){
	printf("  key=%lld,hash=%u,n=%u:\n", p->key[i], p->hash[i], p->n[i]);
	CHashEntry<OFFSET1_T> *phash = &EntryMem[p->hash[i]];
	for(unsigned int k = 0; k < p->n[i]; k++)
	  printf("    mapid1=%d,offset1=%d,Roffset1=%d,win1[0]=%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n",
		 phash[k].mapid1,phash[k].offset1,phash[k].Roffset1,phash[k].win1[0],phash[k].win1[1],phash[k].win1[2],phash[k].win1[3],phash[k].win1[4]);
      }
      if(cnt > COLLIDELINE){
	cnt -= COLLIDELINE;
	CHashCollideR *pr = &CollideMemR[p->next];
	for(unsigned int i = 0; i < cnt; i++){
	  printf("  key=%lld,hash=%u,n=%u:\n", pr[i].key, pr[i].hash, pr[i].n);
	  CHashEntry<OFFSET1_T> *phash = &EntryMem[pr[i].hash];
	  for(unsigned int k=0; k < pr[i].n; k++)
	    printf("    mapid1=%d,offset1=%d,Roffset1=%d,win1[0]=%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n",
		   phash[k].mapid1,phash[k].offset1,phash[k].Roffset1,phash[k].win1[0],phash[k].win1[1],phash[k].win1[2],phash[k].win1[3],phash[k].win1[4]);
	  
	}
      }
    }
    fflush(stdout);
  } // if(VERB>=3)
}

/** Top level call to populate hash table for single map.
*/
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::hash_insert(int mapid1, Cmap *pmap /* rmap[mapid1] */)
{
  if(DEBUG>=2) assert(pmap == rmap[mapid1]);
  if(DEBUG>=2) assert(roffset <= mapid1 && mapid1 < roffset + numrmaps);

  int N = pmap->numsite[hcolor];
  if(N <= HashWin)
    return;
  FLOAT *Y = YY[mapid1];// WAS13 pmap->site[hcolor];
  FLOAT len = Ylen[mapid1];
  PFLOAT win1[MAXR_WIN];

  if(TRACE2 && pmap->id == TraceID1){
    printf("hash_insert(mapid1=%d,pmap):(id=%lld),N=%d:\n",mapid1,pmap->id,N);
    for(int I = 1; I <= N+1; I++)
      printf("  Y[%d]= %0.3f (delta= %0.3f)\n",I,Y[I], Y[I] - Y[I-1]);
    fflush(stdout);
  }

  for(int site1 = 1; site1 <= N-HashWin; site1++){

    OFFSET1_T offset1 = floor(Y[site1]*OffsetKBinv + 0.5);
    OFFSET1_T Roffset1 = floor((len - Y[site1+HashWin])*OffsetKBinv + 0.5);

    if(VERB && (TRACE2 && rmap[mapid1]->id == TraceID1)){
      printf("keyinsert:mapid1=%d(id=%lld),site1=%d(Y[site1]=%0.3f,Y[site1+HashWin]=%0.3f,Ylen=%0.3f,OffsetKB=%0.3f),offset1=%d,Roffset1=%d\n",
	     mapid1,rmap[mapid1]->id,site1,Y[site1],Y[site1+HashWin],len, OffsetKB,offset1,Roffset1);
      fflush(stdout);
    }

    /* use window without site errors */
    window<0>(win1, Y, pmap, mapid1, 0, site1, 0);
    long long key = hashkey(win1);
    keyinsert(key, win1, mapid1, offset1, Roffset1, site1, 0);
    
    if(NumErrI <= 0 || site1 >= N-HashWin)
      continue;

    /* consider single missing site in window */
    Roffset1 = floor((len - Y[site1+HashWin+1])*OffsetKBinv + 0.5);

    if(VERB && (TRACE2 && rmap[mapid1]->id == TraceID1)){
      printf("keyinsert:mapid1=%d(id=%lld),site1=%d(Y[site1+HashWin+1]=%0.3f,Ylen=%0.3f,OffsetKB=%0.3f),Roffset1 -> %d\n",
	     mapid1,rmap[mapid1]->id,site1,Y[site1+HashWin+1],len, OffsetKB,Roffset1);
      fflush(stdout);
    }

    for(int err = 1; err <= HashWin; err++){
      window<0>(win1, Y, pmap, mapid1, 0, site1, err);
      key = hashkey(win1);
      keyinsert(key, win1, mapid1, offset1, Roffset1, site1, err);
    }
    
    if(NumErrI <= 1 || site1 >= N-HashWin-1)
      continue;
    
    /* consider 2 missing sites in window */
    Roffset1 = floor((len - Y[site1+HashWin+2])*OffsetKBinv + 0.5);
    for(int err1 = 1; err1 <= HashWin; err1++)
      for(int err2 = err1+1; err2 <= HashWin+1; err2++){
	window2<0>(win1, Y, pmap, mapid1, 0, site1, err1, err2);
	key = hashkey(win1);
	keyinsert(key,win1,mapid1,offset1,Roffset1,site1, err1 | (err2 << 4));
      }

    if(NumErrI <= 2 || site1 >= N-HashWin-2)
      continue;

    /* Handle 3 missing sites in window */
    printf("CHashTable::hash_insert() not yet implemented for NumErrI = %d\n", NumErrI);    
    fflush(stdout);    exit(1);
  }
}

class CKeyList {
public:
  CHashCollide *p;
  int table;/* index into hashtables[0..numhashtables-1] */
};

static int LonglongInc(register long long *p1, register long long *p2)
{
  register long long n1 = *p1;
  register long long n2 = *p2;
  /* Cannot return n1-n2 since that may be out of range for integer */
  return (n1 > n2) ? 1 : (n1 < n2) ? -1 : 0;
}

static int KeyInc(register CKeyList *p1, register CKeyList *p2)
{
  register long long n1 = p1->p->key;
  register long long n2 = p2->p->key;
  /* Cannot return n1-n2 since that may be out of range for integer */
  return (n1 > n2) ? 1 : (n1 < n2) ? -1 : 0;
}

/** merge multiple hashtables[0..numhashtables-1] into capacted hashtable. Note that "this" hashtable can be one of the hashtables[]. 
    Intended to be used with multiple hashtables that were generated by multithreaded version of hashinert() 

    NOTE : Also reduces memory used if there is only 1 hashtable using hashtable->hashmerge(&hashtable,1)

    if QUICKHASH==1 : Computes tables keytable1[],keytable2[],keytable3[]
    if QUICKHASH==2 : Computes table keytable1[]
 */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::hashmerge(CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> **hashtables, int numhashtables)
{
  /* first divide up the hashtable by index into "numthreads" approximately equal size sections : each section of all hashtables will be handled by a seperate thread */
  unsigned int *hindex_start = new unsigned int[hnumthreads+1];
  unsigned int inc = (1U << HASH_BITS)/hnumthreads, sum = 0;
  for(int t = 0; t < hnumthreads; t++){
    hindex_start[t] = sum;
    sum += inc;
  }
  hindex_start[hnumthreads] = (1U << HASH_BITS);

  /* add up the required length of HashCollide & HashEntry arrays for each section of the hashtable */
  CollideCntMax = CollideCnt = EntryCnt = 0;
  unsigned int *pCollideCnt = (unsigned int*)calloc((hnumthreads+1)*2,sizeof(unsigned int));
  unsigned int *pEntryCnt = &pCollideCnt[hnumthreads+1];
  unsigned int hcnt = 0, CollideCntTot = 0;

  #pragma omp parallel num_threads(hnumthreads) if(hnumthreads > 1)
  {
    unsigned int myCollideCnt = 0, myEntryCnt = 0, myhcnt = 0, myCollideCntTot = 0, myCollideCntMax = 0;
    int tid = 0;
    #ifdef _OPENMP
    tid = omp_get_thread_num();  /* NOTE : each thread will only work on its own section of the hashtable : hindex_start[tid] .. hindex_start[tid+1]-1 */
    #endif
    
    long long keylist[MAXCOLLIDE];

    for(unsigned int hindex = hindex_start[tid]; hindex < hindex_start[tid+1]; hindex++){
      /* Count number of CHashCollide,CHashEntry entries needed for hindex over all hashtables[0..numhashtables-1] */
      /* NOTE: MAXHIT_COVERAGE, STAT and maxEntry updates are done later */
      unsigned int Ccnt = 0, Ecnt = 0;
      for(int h = numhashtables; --h >= 0;){
	register CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *hashtable = hashtables[h];
	unsigned int pindex = hashtable->index[hindex];
	if(!pindex)
	  continue;
	if(DEBUG>=2 && !(pindex < hashtable->CollideBlockCnt)){
	  #pragma omp critical
	  {
	    printf("hindex=%d:h=%d,pindex=%u,CollideBlockCnt=%u\n",hindex,h,pindex,hashtable->CollideBlockCnt);
	    fflush(stdout);
	    assert(pindex < CollideBlockCnt);
	  }
	}
	CHashCollide *p = &hashtable->CollideBlock[pindex];
	for(int i = p->n; i > 0; i--){
	  Ecnt += p[i].n;
	  keylist[Ccnt++] = p[i].key;
	}
	if(Ccnt > MAXCOLLIDE){
	  #pragma omp critical
	  {
	    printf("Ccnt = %u : increase MAXCOLLIDE in hash.cpp (current value= %u)\n", Ccnt, MAXCOLLIDE);
	    fflush(stdout); exit(1);
	  }
	}
      }
      if(Ccnt <= 0)
	continue;

      /* Count unique keys in keylist[0..Ccnt-1] */
      qsort(keylist,Ccnt,sizeof(long long),(intcmp*)LonglongInc);
      long long key = keylist[0];
      int cnt = 1;
      for(unsigned int i = cnt; i < Ccnt; i++){
	if(keylist[i] == key)
	  continue;
	key = keylist[i];
	cnt++;
      }
      Ccnt = cnt;
      if(DEBUG>=2) assert(Ccnt>0);

      myhcnt++;
      myEntryCnt += Ecnt;
      myCollideCnt += max(COLLIDELINE,Ccnt)-COLLIDELINE;
      myCollideCntTot += Ccnt;
      myCollideCntMax = max(myCollideCntMax,Ccnt);
    } // hindex = hindex_start[tid] .. hindex_start[tid+1]-1
    pCollideCnt[tid+1] = myCollideCnt;
    pEntryCnt[tid+1] = myEntryCnt;

    #pragma omp critical
    {
      CollideCnt += myCollideCnt;
      CollideCntTot += myCollideCntTot;
      CollideCntMax = max(CollideCntMax,myCollideCntMax);
      EntryCnt += myEntryCnt;
      hcnt += myhcnt;
    }
  } // parallel

  if(DEBUG>=2){
    unsigned int sum = 0;
    for(int i = 0; i <= hnumthreads; i++){
      sum += pCollideCnt[i];
      //      printf("pCollideCnt[%d]=%u\n",i,pCollideCnt[i]);
    }
    //    printf("sum=%u,CollideCnt=%u\n",sum,CollideCnt);
    //    fflush(stdout);
    assert(sum == CollideCnt);
  }

  /* accumulate pCollideCnt[0..numthreads],pEntryCnt[0..numthreads] */
  for(int i = 1; i <= hnumthreads; i++){
    pCollideCnt[i] += pCollideCnt[i-1];
    //    printf("Cumulative pCollideCnt[%d]=%u\n",i,pCollideCnt[i]);
    pEntryCnt[i] += pEntryCnt[i-1];
  }
  //  fflush(stdout);
  if(DEBUG) assert(pCollideCnt[hnumthreads] == CollideCnt);
  if(DEBUG) assert(pEntryCnt[hnumthreads] == EntryCnt);
  
  /* allocate all memory of Compacted hashtable */
  POSIX_MEMALIGN(indexR,64,sizeof(CHashIndexR)*(1U << HASH_BITS));
  POSIX_MEMALIGN(CollideMemR,64,sizeof(CHashCollideR)*CollideCnt);
  POSIX_MEMALIGN(EntryMem,64,sizeof(CHashEntry<OFFSET1_T>)*EntryCnt);

  if(VERB && hcnt > 0){
    unsigned int CollideBlockCntTot = 0, CollideBlockMaxTot = 0, EntryBlockCntTot = 0, EntryBlockMaxTot = 0;
    for(int i = numhashtables; --i >= 0;){
      CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *hashtable = hashtables[i];
      CollideBlockCntTot += hashtable->CollideBlockCnt;
      CollideBlockMaxTot += hashtable->CollideBlockMax;
      EntryBlockCntTot += hashtable->EntryBlockCnt;
      EntryBlockMaxTot += hashtable->EntryBlockMax;
    }
    if(STAT || (MAXHIT_COVERAGE && hashMaxCov > 0))
      printf("Hash index %u/%u, Average keys per index = %0.3f (max=%d), OverFlow keys per index = %0.3f, Average Number of entries per key= %0.3f(lim=%d):CPU time=%0.6f wall time=%0.6f\n",
	     hcnt, (1U << HASH_BITS), ((double)(CollideCntTot))/hcnt, CollideCntMax, ((double)(CollideCnt))/hcnt, EntryCnt / (double)(CollideCntTot),MaxHit, mtime(), wtime());
    else
      printf("Hash index %u/%u, Average keys per index = %0.3f (max=%d), OverFlow keys per index = %0.3f, Average Number of entries per key= %0.3f:CPU time=%0.6f wall time=%0.6f\n",
	     hcnt,(1U << HASH_BITS), ((double)(CollideCntTot))/hcnt, CollideCntMax, ((double)(CollideCnt))/hcnt, EntryCnt / ((double)(CollideCntTot)), mtime(),wtime());
    if(VERB/* HERE >=2 */){
      printf("    Original:CollideBlockCnt=%u/%u,EntryBlockCnt=%u/%u (Total=%0.4f Gb). Compacted: CollideCnt=%u,EntryCnt=%u (Total=%0.4f Gb), Rand Tables=%lu bytes\n",
	     CollideBlockCntTot, CollideBlockMaxTot, EntryBlockCntTot, EntryBlockMaxTot, 
             (CollideBlockMaxTot*sizeof(CHashCollide) + EntryBlockMaxTot*sizeof(CHashEntry<OFFSET1_T>) + sizeof(unsigned int)*(1U << HASH_BITS)*numhashtables) * 1e-9,
             CollideCnt,EntryCnt, (CollideCnt*sizeof(CHashCollideR) + EntryCnt*sizeof(CHashEntry<OFFSET1_T>) + sizeof(CHashIndexR)*(1U << HASH_BITS))*1e-9, HashSizeRtab());
      //      printf("    TotWin=%lld\n",TotWin);
      //      printf("Sleeping for 60 seconds at HashInsert memory high water mark\n");
      fflush(stdout);
      //      sleep(60);
    }
    fflush(stdout);
  }

  /* now copy data from original hashtable(s) to Compacted hashtable */

  #pragma omp parallel num_threads(hnumthreads) if(hnumthreads > 1)
  {
    int tid = 0;
    #ifdef _OPENMP
    tid = omp_get_thread_num();  /* NOTE : each thread will only work on its own section of the hashtable : hindex_start[tid] .. hindex_start[tid+1]-1 */
    #endif
    unsigned int Ccnt = pCollideCnt[tid], Cmax = pCollideCnt[tid+1];/* Each thread will only write to CollideMem[Ccnt..Cmax-1] */
    unsigned int Ecnt = pEntryCnt[tid], Emax = pEntryCnt[tid+1];/* Each thread will only write to EntryMem[Ecnt..Emax-1] */
    
    CKeyList keylist[MAXCOLLIDE];/* array used to merge CHashCollide arrays from multiple hashtables for a single hindex */
    CHashCollide CollideMerge[MAXCOLLIDE];

    for(unsigned int hindex = hindex_start[tid]; hindex < hindex_start[tid+1]; hindex++){
      /* initialize indexR[hindex] */
      register CHashIndexR *pR = &indexR[hindex];
      pR->cnt = 0;
      pR->next = 0;

      /* assemble keylist[] for hindex */
      unsigned int cnt = 0;
      for(int table = numhashtables; --table >= 0;){
	register CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *hashtable = hashtables[table];
	unsigned int pindex = hashtable->index[hindex];
	if(!pindex)
	  continue;
	CHashCollide *p = &hashtable->CollideBlock[pindex];
	for(int i = p->n; i > 0; i--){
	  CKeyList *pk = &keylist[cnt++];
	  pk->p = &p[i];
	  pk->table = table;
	}
      } // hashtable = hashtables[0 .. numhashtables-1]

      if(cnt<=0)
	continue;

      /* sort keylist[] by key field */
      qsort(keylist,cnt,sizeof(CKeyList),(intcmp*)KeyInc);
      unsigned int end = 0,mcnt = 0;
      do {
	unsigned int start = end;
	long long key = keylist[start].p->key;
	for(end = start+1; end < cnt && key == keylist[end].p->key;)
	  end++;
	/* keylist[start .. end-1] share the same key and will be merged into a single CHashCollide at CollideMerge[mcnt] & EntryMem[Ecnt] */
	CHashCollide *p = &CollideMerge[mcnt++];
	p->key = key;
	p->hash = Ecnt;
	p->n = 0;
	for(unsigned int i = start; i < end; i++){
	  CKeyList *pk = &keylist[i];
	  CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *hashtable = hashtables[pk->table];
	  CHashCollide *pc = pk->p;/* points into Old hashtable */
	  int n = pc->n;
	  if(DEBUG>=2) assert(n > 0);
	  if(MAXHIT_COVERAGE && p->n + n > MaxHit)
	    break;
	  memcpy(&EntryMem[Ecnt],&hashtable->EntryBlock[pc->hash],sizeof(CHashEntry<OFFSET1_T>)*n);
	  p->n += n;
	  Ecnt += n;
	}
	if(DEBUG>=2) assert(p->n <= MASK(16));
	if(DEBUG>=2) assert(end > start);
	if(DEBUG>=2) assert(Ecnt <= Emax);
	if(VERB>=3 && hindex==0)
	  printf("hindex=%d,cnt=%u:start=%u(key=%lld,hash=%u,n=%u),end=%u,mcnt=%u:Ecnt=%u\n",hindex,cnt,start,p->key,p->hash,p->n,end,mcnt,Ecnt);
	
	if(!DisJoint){/* need to sort EntryMem[p->hash ... Ecnt-1] in ascending order of mapid1 */
	  qsort(&EntryMem[p->hash],Ecnt - p->hash,sizeof(CHashEntry<OFFSET1_T>),(intcmp*)CHashEntryMapid1Inc<OFFSET1_T>);

	  if(DEBUG>=2){
	    for(unsigned int i = p->hash+1; i < Ecnt; i++){
	      if(!(EntryMem[i].mapid1 >= EntryMem[i-1].mapid1)){
                #pragma omp critical
		{
		  printf("hindex=%d,cnt=%u:start=%u(key=%lld,hash=%u,n=%u),end=%u,mcnt=%u:Ecnt=%u\n",hindex,cnt,start,p->key,p->hash,p->n,end,mcnt,Ecnt);
		  printf("i=%d:EntryMem[i].mapid1=%d,EntryMem[i-1].mapid1=%d\n",i,EntryMem[i].mapid1,EntryMem[i-1].mapid1);
		  fflush(stdout);
		  assert(EntryMem[i].mapid1 >= EntryMem[i-1].mapid1);
		}
	      }
	    }
	  }
	}
      } while ( end < cnt );

      if(VERB>=3 && hindex==0){
	printf("hindex=%d:cnt=%u,mcnt=%u,mcnt-COLLIDELINE=%u,Ccnt=%u->%u,Cmax=%u\n",hindex,cnt,mcnt,mcnt-COLLIDELINE,Ccnt,Ccnt+mcnt-COLLIDELINE,Cmax);
	for(unsigned int i = 0; i < cnt; i++)
	  printf("  keylist[%u]:table=%d,key=%lld,hash=%u,n=%u\n",i,keylist[i].table,keylist[i].p->key,keylist[i].p->hash,keylist[i].p->n);
	for(unsigned int i = 0; i < mcnt; i++){
	  printf("  CollideMerge[%u]:key=%lld,hash=%u,n=%u\n",i,CollideMerge[i].key,CollideMerge[i].hash,CollideMerge[i].n);
	  CHashEntry<OFFSET1_T> *phash = &EntryMem[CollideMerge[i].hash];
	  if(HashWin==5)
	    for(unsigned int k = 0; k < CollideMerge[i].n; k++)
	      printf("    k=%u:mapid1=%d,offset1=%d,Roffset1=%d,win1=%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n",
		     k,phash[k].mapid1,phash[k].offset1,phash[k].Roffset1,phash[k].win1[0],phash[k].win1[1],phash[k].win1[2],phash[k].win1[3],phash[k].win1[4]);
	}
	fflush(stdout);
      }

      if(VERB>=2 && mcnt >= 2500){/* suspiciously long collision list */
	printf("Suspicously long hashtable collision list (length = %u):\n",mcnt);
	for(unsigned int i = 0; i < mcnt; i++){
	  printf("  CollideMerge[%u]:key=%lld,hash=%u,n=%u:",i,CollideMerge[i].key,CollideMerge[i].hash,CollideMerge[i].n);
	  if(HashWin==5){
	    printf("(%lld,%lld,%lld,%lld,%lld)\n",
		   (CollideMerge[i].key >> (0*SIZE_BITS)) & MASK(SIZE_BITS),
		   (CollideMerge[i].key >> (1*SIZE_BITS)) & MASK(SIZE_BITS),
		   (CollideMerge[i].key >> (2*SIZE_BITS)) & MASK(SIZE_BITS),
		   (CollideMerge[i].key >> (3*SIZE_BITS)) & MASK(SIZE_BITS),
		   (CollideMerge[i].key >> (4*SIZE_BITS)) & MASK(SIZE_BITS));
	    CHashEntry<OFFSET1_T> *phash = &EntryMem[CollideMerge[i].hash];
	    for(unsigned int k = 0; k < CollideMerge[i].n; k++)
	      printf("    k=%u:mapid1=%d,offset1=%d,Roffset1=%d,win1=(%0.3f,%0.3f,%0.3f,%0.3f,%0.3f)\n",
		     k,phash[k].mapid1,phash[k].offset1,phash[k].Roffset1,phash[k].win1[0],phash[k].win1[1],phash[k].win1[2],phash[k].win1[3],phash[k].win1[4]);
	  } else
	    printf("\n");
	}
	fflush(stdout);	exit(1);
      }

      /* copy CollideMerge[0..mcnt-1] to indexR[hindex] AND CollideMemR[Ccnt+(0..mcnt-COLLIDELINE-1)] */
      if(DEBUG>=2) assert(mcnt > 0);
      register CHashCollide *p = CollideMerge;

      register CHashIndexR *pr = &indexR[hindex];
      /* First copy up to to COLLIDELINE entries */
      int imax = min(mcnt,COLLIDELINE);
      int i = 0;
      for(; i < imax; i++){
	pr->key[i] = p[i].key;
	pr->hash[i] = p[i].hash;
	pr->n[i] = p[i].n;
      }
      if(mcnt > COLLIDELINE){
	if(DEBUG>=2) assert(i == imax && i == COLLIDELINE);
	pr->next = Ccnt;
	int mem = mcnt - COLLIDELINE;/* OverFlow[] memory required */
	if(DEBUG>=2) assert(Ccnt + mem <= Cmax);
	//	memcpy(&CollideMemR[Ccnt],&p[COLLIDELINE],mem*sizeof(CHashCollide));
	/* HERE : when CHashCollideR and CHashCollide are made the same, replace following loop with above memcpy() */
	register CHashCollideR *pC = &CollideMemR[Ccnt];
	for(unsigned int i = COLLIDELINE; i < mcnt; i++,pC++){
	  pC->key = p[i].key;
	  pC->hash = p[i].hash;
	  pC->n = p[i].n;
	}
	if(SORT_KEYS) /* sort in order of key */
	  qsort(&CollideMemR[Ccnt],mem,sizeof(CHashCollideR),(intcmp*)CHashCollideRKeyInc);
	Ccnt += mem;
      }
      pr->cnt = mcnt;

    } /* hindex = hindex_start[tid .. hindex_start[tid+1]-1 */
    if(DEBUG) assert(Ccnt == Cmax);
    if(DEBUG) assert(Ecnt <= Emax);/* some HashEntries may have been discarded due to MAXHIT_COVERAGE */
  } // parallel

  if(VERB/* HERE >=2 */){
    getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
    printf("Hash merge/compaction done: freeing memory next (VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb): CPU time=%0.6f, wall time=%0.6f secs\n",
	   VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(),wtime());
    fflush(stdout);
  }

  /* free up original hashtable(s) memory */
  for(int table = numhashtables; --table >= 0;){
    CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *hashtable = hashtables[table];
    if(hashtable->CollideBlock) { free(hashtable->CollideBlock); hashtable->CollideBlock = 0;}
    if(hashtable->EntryBlock) { free(hashtable->EntryBlock); hashtable->EntryBlock = 0;}
  }

  /* free up local memory */
  delete [] hindex_start;
  free(pCollideCnt);

  if(VERB/* HERE >=2 */){
    printf("Hashtable : completed freeing memory: CPU time=%0.6f, wall time=%0.6f secs\n", mtime(),wtime());
    fflush(stdout);
  }

  if(QUICKHASH && hashkeys){/* populate keytable1[] (AND keytable2[],keytable3[] if QUICKHASH==6 */
    if(QUICKHASH==1)
      memset(keytable1,0,(HASHWIN==6 ? 3 * (1U << 18) : 2 * (1U << 19)) * sizeof(unsigned long long));
    if(QUICKHASH==2){
      memset(keytable1,0,(1U << (HASHWIN==6 ? 30 : 24)) * sizeof(unsigned long long));
      
      if(VERB>=2){
        long long key = 0x43a63eb;
	long long key1 = KEYHASH1(key);
	printf("After zero'ing out keytable1[] : key= 0x%llx, key1= 0x%llx, keytable_msk= 0x%llx, keytable1[(key1 >> 6) & keytable_msk]= 0x%llx, bit= %lld, msk= 0x%llx, val= 0x%llx\n",
               (unsigned long long)key, (unsigned long long)key1, (unsigned long long)keytable_msk,(unsigned long long) keytable1[(key1 >> 6) & keytable_msk], 
	       (unsigned long long)(key1 & MASK(6)), (1ull << (key1 & MASK(6))),  (unsigned long long) keytable1[(key1 >> 6) & keytable_msk] & (1ull << (key1 & MASK(6))));
	fflush(stdout);
      }
    }

    for(unsigned int hindex = 0; hindex <  (1U << HASH_BITS); hindex++){
      register CHashIndexR *p = &indexR[hindex];
      register unsigned int cnt = p->cnt;
      register unsigned int imax = min(cnt, COLLIDELINE);
      for(unsigned int i = 0; i < imax; i++){
        register unsigned long long key = p->key[i];
	register unsigned long long key1 = KEYHASH1(key);

	if(QUICKHASH >= 2){
	  // #pragma omp atomic update
          keytable1[(key1 >> 6) & keytable_msk] |= (1ull << (key1 & MASK(6)));

	  if(VERB>=2 && (key == 0x43a63eb || key1== 0x4e4e8af)){
	    printf("hashmerge: hindex= 0x%x:i=%d,p->cnt=%d,p->key[i]= 0x%llx, key1= 0x%llx, keytable1[(key1 >> 6) & keytable_msk]= 0x%llx, bit= %lld, msk= 0x%llx, val= 0x%llx, p->n[i]=%d\n",
		   hindex,i,cnt,(unsigned long long)p->key[i],(unsigned long long)key1,(unsigned long long)keytable1[(key1 >> 6) & keytable_msk], 
	           key1 & MASK(6), (1ull << (key1 & MASK(6))), (unsigned long long) keytable1[(key1 >> 6) & keytable_msk] & (1ull << (key1 & MASK(6))), p->n[i]);
	    fflush(stdout);
	  }
	} else {

	  // #pragma omp atomic update
	  keytable1[key1 & keytable_msk] |= (1ull << ((key1 >> (HASHWIN==6 ? 18 : 19)) & MASK(6)));

	  register long long key2 = KEYHASH2(key);
	  // #pragma omp atomic update
	  keytable2[key2 & keytable_msk] |= (1ull << (HASHWIN==5 ? (key2 >> 19) : (key2 >> 18) & MASK(6)));

	  if(HASHWIN==6){
	    register long long key3 = KEYHASH3(key);
	    // #pragma omp atomic update
	    keytable3[key3 & keytable_msk] |= (1ull << (key3 >> 18));
	  }
	}
      }

      if(cnt <= COLLIDELINE)
	continue;

      /* find keys in overflow Collision list */
      register CHashCollideR *pr = &CollideMemR[p->next];
      cnt -= COLLIDELINE;
      
      for(register unsigned int i = 0; i < cnt;i++){
	register long long key = pr[i].key;
	register long long key1 = KEYHASH1(key);

	if(QUICKHASH >= 2){
	  // #pragma omp atomic update
          keytable1[(key1 >> 6) & keytable_msk] |= (1ull << (key1 & MASK(6)));

	  if(VERB>=2 && (key == 0x43a63eb || key1== 0x4e4e8af)){
	    printf("hashmerge: hindex= 0x%x, p->next=0x%x:i=%d,p->cnt=%d,cnt=%d,pr[i].key= 0x%llx, key1= 0x%llx, keytable[(key1 >> 6) & keytable_msk]= 0x%llx, bit= %lld, msk= 0x%llx, val= 0x%llx(pr[i].n=%u)\n",
		   hindex,p->next,i,p->cnt,cnt,(unsigned long long)pr[i].key,(unsigned long long)key1,(unsigned long long)keytable1[(key1 >> 6) & keytable_msk], 
	           key1 & MASK(6), 1ull << (key1 & MASK(6)), (unsigned long long) keytable1[(key1 >> 6) & keytable_msk] & (1ull << (key1 & MASK(6))),(unsigned int) pr[i].n);
	    fflush(stdout);
	  }
	} else {
	  
	  // #pragma omp atomic update
	  keytable1[key1 & keytable_msk] |= (1ull << ((key1 >> (HASHWIN==6 ? 18 : 19)) & MASK(6)));

	  register long long key2 = KEYHASH2(key);
	  // #pragma omp atomic update
	  keytable2[key2 & keytable_msk] |= (1ull << (HASHWIN==5 ? (key2 >> 19) : (key2 >> 18) & MASK(6)));

	  if(HASHWIN==6){
	    register long long key3 = KEYHASH3(key);
	    // #pragma omp atomic update
	    keytable3[key3 & keytable_msk] |= (1ull << (key3 >> 18));
	  }
	}
      }
    }

    if(VERB/* HERE >=2 */){
      printf("Hashtable : created 1-bit key tables: CPU time=%0.6f, wall time=%0.6f secs\n", mtime(),wtime());
      fflush(stdout);
    }
    if(VERB>=2){
      long long key = 0x43a63eb;
      long long key1 = KEYHASH1(key);
      printf("After creating keytable1[] : key= 0x%llx, key1= 0x%llx, keytable_msk= 0x%llx, keytable1[(key1 >> 6) & keytable_msk]= 0x%llx, bit= %lld, msk= 0x%llx, val= 0x%llx\n",
	     (unsigned long long)key, (unsigned long long)key1, (unsigned long long)keytable_msk,(unsigned long long) keytable1[(key1 >> 6) & keytable_msk], 
             (unsigned long long) (key1 & MASK(6)), (1ull << (key1 & MASK(6))), (unsigned long long) keytable1[(key1 >> 6) & keytable_msk] & (1ull << (key1 & MASK(6))));
      fflush(stdout);
    }
  } // if(QUICKHASH && hashkeys)

  sorted = 1;/* Hashtable should not be modified after this */
}

/** compact and sort hashtable entries for faster retrieval : also collect statistics on hashtable performance (to tune HASH_BITS) */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::hashsort()
{
  /* first add up the total length of HashCollide arrays and HashEntry arrays */
  CollideCntMax = CollideCnt = EntryCnt = 0;
  unsigned int hcnt = 0, CollideCntTot = 0;
  int maxEntry = 0;

  #pragma omp parallel num_threads(hnumthreads) if(hnumthreads > 1)
  {
    unsigned int myCollideCnt = 0, myEntryCnt = 0, myhcnt = 0, myCollideCntTot = 0, myCollideCntMax = 0, mymaxEntry = 0;

    #pragma omp for schedule(static,256)
    for(int hindex = MASK(HASH_BITS); hindex >= 0; hindex--){
      unsigned int pindex = index[hindex];
      if(!pindex)
	continue;
      CHashCollide *p = &CollideBlock[pindex];

      myhcnt++;
      unsigned int Ccnt = p->n;
      for(int i = p->n; i > 0;i--){
	unsigned int n = p[i].n;
	if(STAT && n > mymaxEntry)
	  mymaxEntry = n;
	if(MAXHIT_COVERAGE && n > MaxHit){
	  if(!STAT && n > mymaxEntry)
	    mymaxEntry = n;
	  Ccnt--;
	  continue;
	}
	myEntryCnt += n;
      }
      myCollideCnt += max(0,Ccnt-COLLIDELINE);
      myCollideCntTot += Ccnt;
      myCollideCntMax = max(myCollideCntMax,Ccnt);
    }
    #pragma omp critical
    {
      CollideCnt += myCollideCnt;
      CollideCntTot += myCollideCntTot;
      CollideCntMax = max(CollideCntMax,myCollideCntMax);
      EntryCnt += myEntryCnt;
      hcnt += myhcnt;
      maxEntry = max(maxEntry,mymaxEntry);
    }
  }

  if(VERB && hcnt > 0){
    if(VERB/* HERE >=2 */){
      printf("    Original:CollideBlockCnt=%u/%u,EntryBlockCnt=%u/%u (Total=%lu bytes). Compacted: CollideCnt=%u,EntryCnt=%u (Total=%lu bytes)\n",
	     CollideBlockCnt, CollideBlockMax, EntryBlockCnt, EntryBlockMax, 
	     CollideBlockMax*sizeof(CHashCollide) + EntryBlockMax*sizeof(CHashEntry<OFFSET1_T>) + sizeof(unsigned int)*(1U << HASH_BITS) + HashSizeRtab(), 
	     CollideCnt,EntryCnt, CollideCnt*sizeof(CHashCollideR) + EntryCnt*sizeof(CHashEntry<OFFSET1_T>) + sizeof(CHashIndexR)*(1U << HASH_BITS) + HashSizeRtab());
    }
    if(STAT || (MAXHIT_COVERAGE && hashMaxCov > 0))
      printf("Hash index %u/%u, Average keys per index = %0.3f (max=%d), OverFlow keys per index = %0.3f, Average Number of entries per key= %0.3f(max=%d,lim=%d):CPU time=%0.6f wall time=%0.6f\n",
	     hcnt, (1U << HASH_BITS), ((double)(CollideCntTot))/hcnt, CollideCntMax, ((double)(CollideCnt))/hcnt, EntryCnt / (double)(CollideCntTot),maxEntry,MaxHit,mtime(),wtime());
    else
      printf("Hash index %u/%u, Average keys per index = %0.3f (max=%d), OverFlow keys per index = %0.3f, Average Number of entries per key= %0.3f:CPU time=%0.6f wall time=%0.6f\n",
	     hcnt,(1U << HASH_BITS), ((double)(CollideCntTot))/hcnt, CollideCntMax, ((double)(CollideCnt))/hcnt, EntryCnt / ((double)(CollideCntTot)), mtime(),wtime());
    fflush(stdout);
  }

  /* allocate all memory of Compacted hashtable */
  POSIX_MEMALIGN(indexR,64,sizeof(CHashIndexR)*(1U << HASH_BITS));
  POSIX_MEMALIGN(CollideMemR,64,sizeof(CHashCollideR)*CollideCnt);
  POSIX_MEMALIGN(EntryMem,64,sizeof(CHashEntry<OFFSET1_T>)*EntryCnt);
  
  /* now copy data from original hashtable to Compacted hashtable */

  #pragma omp parallel for num_threads(hnumthreads) schedule(static,256) if(hnumthreads > 1)
  for(int hindex = MASK(HASH_BITS); hindex >= 0; hindex--){
    register CHashIndexR *p = &indexR[hindex];
    p->cnt = 0;
    p->next = 0;
  }

  unsigned int Ccnt = 0, Ecnt = 0; 

  if(!MAXHIT_COVERAGE || MaxHit > MASK(16)){
    printf("NEWHASH requires MAXHIT_COVERAGE and MaxHit <= %d\n", MASK(16));
    fflush(stdout); exit(1);
  }

  for(int hindex = MASK(HASH_BITS); hindex >= 0; hindex--){
    unsigned int pindex = index[hindex];
    if(!pindex)
      continue;
    CHashCollide *p = &CollideBlock[pindex];
    int n = p->n;

    int maxhits = 0;
    unsigned int pnsum = 0;
    for(int i = 1; i <= n; i++){
      unsigned int pn = p[i].n;
      if(MAXHIT_COVERAGE && pn > MaxHit){
	maxhits++;
	continue;
      }
      pnsum += pn;
    }
    if(DEBUG>=2) assert(n+1 > maxhits);
    int mem = max(0, n - maxhits - COLLIDELINE);/* OverFlow[] memory required */

    /* reallocate memory pointed to by index[hindex] */
    unsigned int Ecntstart,EcntEnd, Ccntstart;
    {
      Ccntstart = Ccnt;
      Ccnt += mem;
      Ecntstart = Ecnt;
      Ecnt += pnsum;
      if(DEBUG>=2) EcntEnd = Ecnt;
    }

    /* reallocate p[1..n].hash[] */
    register CHashCollideR *newCollide = &CollideMemR[Ccntstart];
    register CHashIndexR *pr = &indexR[hindex];
    register int cnt = 0;
    register int i = 1;
    for(; i <= n; i++){
      unsigned int pn = p[i].n;
      if(MAXHIT_COVERAGE && pn > MaxHit)
	continue;
      pr->key[cnt] = p[i].key;
      pr->hash[cnt] = Ecntstart;
      pr->n[cnt] = pn;
      CHashEntry<OFFSET1_T> *newhash = &EntryMem[Ecntstart];
      memcpy(newhash,&EntryBlock[p[i].hash],pn*sizeof(CHashEntry<OFFSET1_T>));
      if(DEBUG>=2){/* check that newhash[0..pn-1].win1[0..HashWin-1] is valid (NOT NaN) */
	for(register int k = pn;--k >= 0;){
	  register CHashEntry<OFFSET1_T> *phash = &newhash[k];
	  for(register int w = HashWin; --w >= 0;){
	    if(!(isfinite(phash->win1[w]))){
	      printf("hindex=%d,pindex=%u:p->n=%u,i=%d,p[i].n=%u,p[i].hash=%u,k=%d,w=%d,EntryBlock[p[i].hash+k].win1[w]=%e,phash->win1[w]=%e\n",
		     hindex,pindex,p->n,i,p[i].n,p[i].hash,k,w,EntryBlock[p[i].hash+k].win1[w],phash->win1[w]);
	      fflush(stdout);
	      assert(isfinite(phash->win1[w]));
	    }
	  }
	}
      }

      Ecntstart += pn;
      if(++cnt >= COLLIDELINE){
	pr->next = Ccntstart;
	break;
      }
    }
    for(i++; i <= n; i++){/* process OverFlow collisions */
      unsigned int pn = p[i].n;
      if(MAXHIT_COVERAGE && pn > MaxHit)
	continue;
      newCollide->key = p[i].key;
      newCollide->hash = Ecntstart;
      newCollide->n = pn;
      newCollide++;
      cnt++;

      CHashEntry<OFFSET1_T> *newhash = &EntryMem[Ecntstart];
      memcpy(newhash,&EntryBlock[p[i].hash],pn*sizeof(CHashEntry<OFFSET1_T>));
      if(DEBUG>=2){/* check that newhash[0..pn-1].win1[0..HashWin-1] is valid (NOT NaN) */
	for(register int k = pn;--k >= 0;){
	  register CHashEntry<OFFSET1_T> *phash = &newhash[k];
	  for(register int w = HashWin; --w >= 0;)
	    assert(isfinite(phash->win1[w]));
	}
      }
      Ecntstart += pn;
    }
    pr->cnt = cnt;

    if(DEBUG>=2) assert(cnt == n-maxhits);
    if(DEBUG>=2) assert(mem == max(0, cnt-COLLIDELINE));
    if(DEBUG>=2) assert(Ecntstart == EcntEnd);

    if(SORT_KEYS) /* sort OverFlow[Ccntstart .. Ccntstart + mem-1] in order of key */
      qsort(&CollideMemR[Ccntstart],mem,sizeof(CHashCollideR),(intcmp*)CHashCollideRKeyInc);
  }
  if(DEBUG) assert(Ccnt == CollideCnt);
  if(DEBUG) assert(Ecnt == EntryCnt);
  
  delete [] index;  index = 0;
  if(CollideBlock) { free(CollideBlock); CollideBlock = 0;}
  if(EntryBlock) { free(EntryBlock); EntryBlock = 0;}

  sorted = 1;/* Hashtable should not be modified after this */
}

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline int CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::hashfind(long long key, CHashEntry<OFFSET1_T>* &hash, CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>* matchtable)
{
  if(STAT) matchtable->HTindcnt++;

  unsigned int hindex = hashindex(key);
  if(DEBUG>=2) assert(hindex <= (unsigned int)MASK(HASH_BITS));
  register CHashIndexR *p = &indexR[hindex];
  register unsigned int cnt = p->cnt;
  if(!cnt)
    return 0;
  
  if(STAT) matchtable->HTfindcnt++;

  if(cnt <= COLLIDELINE){/*just look for key in indexR[hindex] */
    for(register unsigned int i = 0; i < cnt; i++){
      if(STAT) matchtable->HTqcolcnt++;
      if(p->key[i] == key){
	hash = &EntryMem[p->hash[i]];
	if(STAT) matchtable->HTmatchcnt++;
	if(STAT) matchtable->HTentrycnt += p->n[i];
	return p->n[i];
      }
    }
    return 0;
  }

  /* first look for key in indexR[hindex] */
  for(register unsigned int i = 0; i < COLLIDELINE; i++){
    if(STAT) matchtable->HTqcolcnt++;
    if(p->key[i] == key){
      hash = &EntryMem[p->hash[i]];
      if(STAT) matchtable->HTmatchcnt++;
      if(STAT) matchtable->HTentrycnt += p->n[i];
      return p->n[i];
    }
  }

  /* find key in overflow Collision list */
  register CHashCollideR *pr = &CollideMemR[p->next];
  cnt -= COLLIDELINE;

  if(!SORT_KEYS){
    for(register unsigned int i = 0; i < cnt; i++){
      if(STAT) matchtable->HTcolcnt++;
      if(pr[i].key == key){
	hash = &EntryMem[pr[i].hash];
	if(STAT) matchtable->HTmatchcnt++;
	if(STAT) matchtable->HTcolcnt--;
	if(STAT) matchtable->HTentrycnt += pr[i].n;
	return pr[i].n;
      }
    }
    return 0;
  } else { /* SORT_KEYS == 1 : use binary search on sorted key values */
    register int low = 0;
    register int high = cnt - 1;
    while(high > low){
      if(STAT) matchtable->HTcolcnt++;
      register int mid = (low+high)/2;
      register long long midKey = pr[mid].key;
      if(midKey < key)
	low = mid+1;
      else if(midKey > key)
	high = mid-1;
      else {/* found */
        if(DEBUG>=2) assert(pr[mid].key == key);
	hash = &EntryMem[pr[mid].hash];
	if(STAT) matchtable->HTmatchcnt++;
	if(STAT) matchtable->HTentrycnt += pr[mid].n;
	return pr[mid].n;
      }
    }
    if(high == low){
      if(STAT) matchtable->HTcolcnt++;
      if(pr[low].key == key){
	hash = &EntryMem[pr[low].hash];
	if(STAT) matchtable->HTmatchcnt++;
	if(STAT) matchtable->HTcolcnt--;
	if(STAT) matchtable->HTentrycnt += pr[low].n;
	return pr[low].n;
      }
    }

    if(DEBUG>=2){/* check to make sure we didn't miss a match (due to bug in binary search) */
      for(register unsigned int i = 0; i < cnt; i++){
	if(pr[i].key == key){
          #pragma omp critical 
	  {
	    printf("binary search failed to find key=%lld in pr[0..%d].key (low=%u,high=%u):\n",key,cnt-1,low,high);
	    for(unsigned int j = 0; j < cnt; j++)
	      printf("pr[%d].key=%lld\n",j,pr[j].key);
	    fflush(stdout);
	    assert(pr[i].key != key);
	  }
	}
      }
    }
  } // SORT_KEYS==1

  return 0;
}

#define PD (USE_MIC ? 4 : 8) /* Prefetch distance : Powers of 2 will be slightly faster: 1,2,4,8 */
#define PD2 (PD*2)

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
__attribute__ ((noinline)) int CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::hashfindV(register long long key, register CHashResult *hash, CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *matchtable)
{
  register int rcnt = 0;/* count how many key values had results */

  if(QUICKHASH && hashkeys){ // NEW code using matchtable->sizekeyV[],rand1V[] : slower than old code, unless QUICKHASH
    _mm_prefetch(CAST SizeDelta, _MM_HINT_T0); /* Sequential access for SizeDelta[] of under 2K, so only single prefetch needed */
    register long long *sizekeyV = matchtable->sizekeyV;
    register unsigned int *rand1V = matchtable->rand1V;

    register unsigned int Smax, Slen;

    if(QUICKHASH){
      register int numsizekey = 0;

      if(QUICKHASH==1){
	/* Use two or three 1-bit tables (24 bit index, 2 Mbytes for each table) to filter output key + SizeDelta[S] values that have no entries in hashtable.
	   1. For HASHWIN==5, index tables by [sizekey & MASK(24)] and [sizekey >> 6] : if either table is 0, skip this sizekey
	   2. For HASHWIN==6, index tables by [sizekey & MASK(24)] and [(sizekey >> 6) & MASK(24)] and [sizekey >> 12] : if either table is 0, skip this sizekey

	   sizekey values that pass are appended to matchtable->sizekeyV[0..numsizekey-1]
	*/

	// HERE HERE : use circular buffer key1V[2*PD] to store key1 = KEYHASH1(sizekey), so it is not compute more than once. key1V[] should be in L1 cache ???

	for(register unsigned int S = 0; S < PD; S++){// prefetch loop for keytable1[],keytable2[] (each are 2Mb and should be in L3 cache)
	  register long long sizekey = key + SizeDelta[S];      
	  _mm_prefetch(CAST &keytable1[KEYHASH1(sizekey) & keytable_msk], _MM_HINT_T0);
	  _mm_prefetch(CAST &keytable2[KEYHASH2(sizekey) & keytable_msk], _MM_HINT_T0);
	  if(HASHWIN==6)
	    _mm_prefetch(CAST &keytable3[KEYHASH3(sizekey) & keytable_msk], _MM_HINT_T0);
	}
	_mm_prefetch(CAST sizekeyV, _MM_HINT_T0); /* Sequential write access for sizekeyV[] of under 2K, so only single prefetch needed */
	Smax = max(0, SizeLen - PD/2);
	for(register unsigned int S = 0; S < Smax; S++){
	  register unsigned int SD = S + PD/2;
	  register long long sizekeyPD = key + SizeDelta[SD];
	  _mm_prefetch(CAST &keytable1[KEYHASH1(sizekeyPD) & keytable_msk], _MM_HINT_T0);
	  _mm_prefetch(CAST &keytable2[KEYHASH2(sizekeyPD) & keytable_msk], _MM_HINT_T0);
	  if(HASHWIN==6)
	    _mm_prefetch(CAST &keytable3[KEYHASH3(sizekeyPD) & keytable_msk], _MM_HINT_T0);
      
	  register long long sizekey = key + SizeDelta[S];
	  register long long key1 = KEYHASH1(sizekey);
	  if(!(keytable1[key1 & keytable_msk] & (1ull << ((key1 >> (HASHWIN==6 ? 18 : 19)) & MASK(6)))))
	    continue;
	  register long long key2 = KEYHASH2(sizekey);
	  if(!(keytable2[key2 & keytable_msk] & (1ull << (HASHWIN==5 ? (key2 >> 19) : (key2 >> 18) & MASK(6)))))
	    continue;
	  if(HASHWIN==6){
	    register long long key3 = KEYHASH3(sizekey);
	    if(!(keytable3[key3 & keytable_msk] & (1ull << (key3 >> 18))))
	      continue;
	  }
	  sizekeyV[numsizekey++] = sizekey;
	}

	for(register unsigned int S = Smax; S < SizeLen; S++){// postloop
	  register long long sizekey = key + SizeDelta[S];
	  register long long key1 = KEYHASH1(sizekey);
	  if(!(keytable1[key1 & keytable_msk] & (1ull << ((key1 >> (HASHWIN==6 ? 18 : 19)) & MASK(6)))))
	    continue;
	  register long long key2 = KEYHASH2(sizekey);
	  if(!(keytable2[key2 & keytable_msk] & (1ull << (HASHWIN==5 ? (key2 >> 19) : (key2 >> 18) & MASK(6)))))
	    continue;
	  if(HASHWIN==6){
	    register long long key3 = KEYHASH3(sizekey);
	    if(!(keytable3[key3 & keytable_msk] & (1ull << (key3 >> 18))))
	      continue;
	  }
	  sizekeyV[numsizekey++] = sizekey;
	}

      } else if(QUICKHASH==2){
	/* Use large 1-bit table with 30 (or 36) bit index, 128 Mbytes (or 1 Gbyte) to filter all output key + SizeDelta[S] values that have no entries in hashtable.
	   1. For HASHWIN==5, index tables by [sizekey & MASK(30)]
	   2. For HASHWIN==6, index tables by [sizekey & MASK(36)]

	   sizekey values that pass are appended to matchtable->sizekeyV[0..numsizekey-1]
	*/

	if(QUICKSTAT){ /* check effectiveness of KEYHASH1() by checking now many values of KEYHASH1(key + SizeDelta[0..SizeLen-1]) & (~MASK(9)) differ from each other
			  This is the number of cache lines that need to read from memory and will be at least 2. Compute average value and examine cases with values >> 2 */
	  #pragma omp atomic update
	  hashfindV_calls++;
	  
	  long long key1V[SizeLen];
	  
	  for(unsigned int S = 0; S < SizeLen; S++){
	    long long sizekey = key + SizeDelta[S];
	    key1V[S] = KEYHASH1(sizekey) >> 9;
	  }
	  qsort(key1V,SizeLen,sizeof(long long), (intcmp*)key1Inc);

	  int cnt = 1;
	  for(unsigned int S = 1; S < SizeLen; S++)
	    if(key1V[S] != key1V[S-1])
	      cnt++;

	  #pragma omp atomic update
	  hashfindV_lines += cnt;

	  if(cnt > hashfindV_max){
	    #pragma omp critical
	    {
	      hashfindV_max = max(cnt, hashfindV_max);

	      if(cnt >= 32){
		printf("hashfindV() required %d cache lines with QUICKHASH=%d:\n",cnt,QUICKHASH);
		for(unsigned int S = 0; S < SizeLen; S++){
		  long long sizekey = key + SizeDelta[S];
		  long long key1 = KEYHASH1(sizekey);
		  printf("S=%d/%d: sizekey= 0x%llx, key1= 0x%llx, (key1 >> 9)= 0x%llx, (key1 >> 10)= 0x%llx\n", S,SizeLen, sizekey, key1, key1 >> 9, key1 >> 10);
		}
		fflush(stdout);
	      }
	    }
	  }
	}

	// HERE HERE : use circular buffer key1V[2*PD] to store key1 = KEYHASH1(sizekey), so it is not compute more than once. key1V[] should be in L1 cache ???

	for(register unsigned int S = 0; S < PD; S++){// prefetch loop for keytable1[],keytable2[] (each are 2Mb and should be in L3 cache)
	  register long long sizekey = key + SizeDelta[S];      
	  _mm_prefetch(CAST &keytable1[(KEYHASH1(sizekey) >> 6) & keytable_msk], _MM_HINT_T0);
	}
	_mm_prefetch(CAST sizekeyV, _MM_HINT_T0); /* Sequential write access for sizekeyV[] of under 2K, so only single prefetch needed */
	Smax = max(0, SizeLen - PD/2);
	for(register unsigned int S = 0; S < Smax; S++){
	  register unsigned int SD = S + PD/2;
	  register long long sizekeyPD = key + SizeDelta[SD];
	  _mm_prefetch(CAST &keytable1[(KEYHASH1(sizekeyPD) >> 6) & keytable_msk], _MM_HINT_T0);
      
	  register long long sizekey = key + SizeDelta[S];
	  register long long key1 = KEYHASH1(sizekey);
	  if(VERB>=2 && sizekey == 0x43a63eb){
	    printf("S=%d,Smax=%d,key= 0x%llx, sizekey= 0x%llx, numsizekey=%d, keytable_msk= 0x%llx, key1= 0x%llx, keytable1[(key1 >> 6) & keytable_msk]= 0x%llx, bit=%lld, msk= 0x%llx, val=0x%llx\n",
		   S,Smax,(unsigned long long) key, (unsigned long long)sizekey,numsizekey,(unsigned long long)keytable_msk,(unsigned long long)key1, 
		   (unsigned long long)keytable1[(key1 >> 6) & keytable_msk],
		   (unsigned long long)(key1 & MASK(6)), (1ull << (key1 & MASK(6))), (unsigned long long) keytable1[(key1 >> 6) & keytable_msk] & (1ull << (key1 & MASK(6))));
	    fflush(stdout);
	  }
	  if(!(keytable1[(key1 >> 6) & keytable_msk] & (1ull << (key1 & MASK(6)))))
	    continue;

	  sizekeyV[numsizekey++] = sizekey;
	}

	for(register unsigned int S = Smax; S < SizeLen; S++){// postloop
	  register long long sizekey = key + SizeDelta[S];
	  register long long key1 = KEYHASH1(sizekey);
	  if(VERB>=2 && sizekey == 0x43a63eb){
	    printf("S=%d,Smax=%d,key= 0x%llx, sizekey= 0x%llx, numsizekey=%d, keytable_msk= 0x%llx, key1= 0x%llx, keytable1[(key1 >> 6) & keytable_msk]= 0x%llx, bit= %lld, msk= 0x%llx, val= 0x%llx\n",
		   S,Smax,(unsigned long long) key, (unsigned long long)sizekey,numsizekey,(unsigned long long)keytable_msk,(unsigned long long)key1,
		   (unsigned long long)keytable1[(key1 >> 6) & keytable_msk], 
		   (unsigned long long)(key1 & MASK(6)), (1ull << (key1 & MASK(6))), (unsigned long long)keytable1[(key1 >> 6) & keytable_msk] & (1ull << (key1 & MASK(6))));
	    fflush(stdout);
	  }
	  if(!(keytable1[(key1 >> 6) & keytable_msk] & (1ull << (key1 & MASK(6)))))
	    continue;

	  sizekeyV[numsizekey++] = sizekey;
	}

      } // QUICKHASH==2

      _mm_prefetch(CAST sizekeyV, _MM_HINT_T0); /* Sequential access for sizekeyV[] of under 2K, so only single prefetch needed */

      Smax = min(numsizekey, PD/2);
      for(register unsigned int S = 0; S < Smax; S++){ // prefetch loop for rand1[] & rand2[] (each are 64K and should be in L2 cache)
	register long long sizekey = sizekeyV[S];
	_mm_prefetch(CAST &rand1[ sizekey & rhash_msk], _MM_HINT_T0);
	_mm_prefetch(CAST &rand2[ (sizekey >> RHASH_BITS) & rhash_msk], _MM_HINT_T0);
      }
      Smax = max(0, numsizekey - PD/2);
      for(register unsigned int S = 0; S < Smax;S++){
	register unsigned int SD = S + PD/2;
	register long long sizekeyPD = sizekeyV[SD];
	_mm_prefetch(CAST &rand1[ sizekeyPD & rhash_msk], _MM_HINT_T0);
	_mm_prefetch(CAST &rand2[ (sizekeyPD >> RHASH_BITS) & rhash_msk], _MM_HINT_T0);
      
	register long long sizekey = sizekeyV[S];
	rand1V[S] = rand1[sizekey & rhash_msk] ^ rand2[(sizekey >> RHASH_BITS) & rhash_msk];
      }
      Slen = numsizekey;
      for(register unsigned int S = Smax; S < Slen; S++){// postloop
	register long long sizekey = sizekeyV[S];
	rand1V[S] = rand1[sizekey & rhash_msk] ^ rand2[(sizekey >> RHASH_BITS) & rhash_msk];
      }

    } else {// !QUICKHASH && hashkeys : this branch will never execute
      _mm_prefetch(CAST sizekeyV, _MM_HINT_T0); /* Sequential write access for sizekeyV[] of under 2K, so only single prefetch needed */

      for(register unsigned int S = 0; S < PD/2; S++){ // prefetch loop for rand1[] & rand2[] (each are 64K and should be in L2 cache)
	register long long sizekey = key + SizeDelta[S];
	_mm_prefetch(CAST &rand1[ sizekey & rhash_msk], _MM_HINT_T0);
	_mm_prefetch(CAST &rand2[ (sizekey >> RHASH_BITS) & rhash_msk], _MM_HINT_T0);
      }
      Smax = max(0, SizeLen - PD/2);
      for(register unsigned int S = 0; S < Smax;S++){
	register unsigned int SD = S + PD/2;
	register long long sizekeyPD = key + SizeDelta[SD];
	_mm_prefetch(CAST &rand1[ sizekeyPD & rhash_msk], _MM_HINT_T0);
	_mm_prefetch(CAST &rand2[ (sizekeyPD >> RHASH_BITS) & rhash_msk], _MM_HINT_T0);
      
	register long long sizekey = key + SizeDelta[S];
	rand1V[S] = rand1[sizekey & rhash_msk] ^ rand2[(sizekey >> RHASH_BITS) & rhash_msk];
	sizekeyV[S] = sizekey;
      }
      Slen = SizeLen;
      for(register unsigned int S = Smax; S < Slen; S++){// postloop
	register long long sizekey = key + SizeDelta[S];
	rand1V[S] = rand1[sizekey & rhash_msk] ^ rand2[(sizekey >> RHASH_BITS) & rhash_msk];
	sizekeyV[S] = sizekey;
      }
    } // !QUICKHASH && hashkeys

    Smax = min(Slen,PD);
    for(register unsigned int S = 0; S < Smax; S++) // prefetch loop for indexR[]
      _mm_prefetch(CAST &indexR[ rand1V[S] ], _MM_HINT_T0);// rand1V is less than 1 kb, so should still be in L1 cache

    if(STAT) matchtable->HTindcnt += SizeLen;
    if(STAT && QUICKHASH) matchtable->HTqindcnt += Slen;

    //    register CHashIndexR *HashIndexV = matchtable->HashIndexV;
    _mm_prefetch(CAST sizekeyV, _MM_HINT_T0); /* Sequential access for sizekeyV[] of under 2K, so only single prefetch needed */

    Smax = max(0, Slen - PD);
    for(register unsigned int S = 0; S < Smax; S++){
      register unsigned int SD = S + PD;/* look ahead iteration */
      _mm_prefetch(CAST &indexR[ rand1V[SD] ], _MM_HINT_T0);

      register CHashIndexR *p = &indexR[ rand1V[S] ];/* rand1V[] is less than 1 kb, so should still be in L1 cache */
      register long long sizekey = sizekeyV[S];

      register unsigned int cnt = p->cnt;// NOTE : typically a cache fault (and TLB fault as well)
      if(DEBUG>=2 && QUICKHASH==2) assert(cnt > 0);
      if(QUICKHASH <= 1 && !cnt)
	continue;

      if(STAT) matchtable->HTfindcnt++;

      if(DEBUG>=2 && !(cnt <= CollideCntMax)){
	unsigned int hindex = rand1V[S];
	printf("hindex=%u:indexR[hindex].cnt=%u,cnt=%u,CollideCntMax=%u\n",hindex,indexR[hindex].cnt,cnt,CollideCntMax);
	fflush(stdout);
	assert(cnt <= CollideCntMax);
      }

      register unsigned int i;    

      if(cnt <= COLLIDELINE){/* just look for key in indexR[hindex] */
	for(i = 0; i < cnt; i++){
	  if(STAT) matchtable->HTqcolcnt++;
	  if(p->key[i] == sizekey){
	    hash[rcnt].hash = p->hash[i];
	    if(STAT) matchtable->HTmatchcnt++;
	    if(STAT) matchtable->HTentrycnt += p->n[i];
	    hash[rcnt++].n = p->n[i];
	    break;
	  }
	}

	if(DEBUG>=2 && QUICKHASH==2 && !(i < cnt)){
#pragma omp critical
	  {
	    printf("S=%d/%d: i=%d,cnt=%d, key= 0x%llx(key1=0x%llx), sizekeyV[S]= 0x%llx(key1=0x%llx):\n",
		   S,Slen,i,cnt,(unsigned long long)key, (unsigned long long)KEYHASH1(key), (unsigned long long)sizekeyV[S], (unsigned long long)KEYHASH1(sizekeyV[S]));
	    for(unsigned int t = 0; t < cnt; t++)
	      printf("\t indexR[rand1V[S]].key[%u] = 0x%llx\n",t, (unsigned long long)indexR[rand1V[S]].key[t]);
	    fflush(stdout);
	    assert(i < cnt);
	  }
	}
	continue;
      }

      register CHashCollideR *pr = &CollideMemR[p->next];
      _mm_prefetch(CAST pr, _MM_HINT_NTA);    

      /* first look for key in indexR[hindex] */
      for(i = 0; i < COLLIDELINE; i++){
	if(STAT) matchtable->HTqcolcnt++;
	if(p->key[i] == sizekey){
	  hash[rcnt].hash = p->hash[i];
	  if(STAT) matchtable->HTmatchcnt++;
	  if(STAT) matchtable->HTentrycnt += p->n[i];
	  hash[rcnt++].n = p->n[i];
	  break;
	}
      }
      if(i < COLLIDELINE)
	continue;

      /* find key in overflow Collision list */
      cnt -= COLLIDELINE;

      for(i = 0; i < cnt; i++){
	if(STAT) matchtable->HTcolcnt++;
	if(pr[i].key == sizekey){
	  hash[rcnt].hash = pr[i].hash;
	  if(STAT) matchtable->HTmatchcnt++;
	  if(STAT) matchtable->HTentrycnt += pr[i].n;
	  hash[rcnt++].n = pr[i].n;
	  break;
	}
      }
      if(DEBUG>=2 && QUICKHASH==2) assert(i < cnt);
    }

    for(register unsigned int S = Smax; S < Slen; S++){// post loop
      register CHashIndexR *p = &indexR[ rand1V[S] ];
      register unsigned int cnt = p->cnt;      
      if(DEBUG>=2 && QUICKHASH==2) assert(cnt > 0);
      if(QUICKHASH <= 1 && !cnt)
	continue;

      register long long sizekey = sizekeyV[S];
      if(STAT) matchtable->HTfindcnt++;

      register unsigned int i;    

      if(cnt <= COLLIDELINE){/* just look for key in indexR[hindex] */
	for(i = 0; i < cnt; i++){
	  if(STAT) matchtable->HTqcolcnt++;
	  if(p->key[i] == sizekey){
	    hash[rcnt].hash = p->hash[i];
	    if(STAT) matchtable->HTmatchcnt++;
	    if(STAT) matchtable->HTentrycnt += p->n[i];
	    hash[rcnt++].n = p->n[i];
	    break;
	  }
	}

	if(DEBUG>=2 && QUICKHASH==2 && !(i < cnt)){
#pragma omp critical
	  {
	    printf("S=%d/%d: i=%d,cnt=%d, key= 0x%llx(key1=0x%llx), sizekeyV[S]= 0x%llx(key1=0x%llx):\n",
		   S,Slen,i,cnt,(unsigned long long)key, (unsigned long long)KEYHASH1(key), (unsigned long long)sizekeyV[S], (unsigned long long)KEYHASH1(sizekeyV[S]));
	    for(unsigned int t = 0; t < cnt; t++)
	      printf("\t indexR[rand1V[S]].key[%u] = 0x%llx\n",t, (unsigned long long)indexR[rand1V[S]].key[t]);
	    fflush(stdout);
	    assert(i < cnt);
	  }
	}
	continue;
      }

      register CHashCollideR *pr = &CollideMemR[p->next];
      _mm_prefetch(CAST pr, _MM_HINT_NTA);    

      /* first look for key in indexR[hindex] */
      for(i = 0; i < COLLIDELINE; i++){
	if(STAT) matchtable->HTqcolcnt++;
	if(p->key[i] == sizekey){
	  hash[rcnt].hash = p->hash[i];
	  if(STAT) matchtable->HTmatchcnt++;
	  if(STAT) matchtable->HTentrycnt += p->n[i];
	  hash[rcnt++].n = p->n[i];
	  break;
	}
      }
      if(i < COLLIDELINE)
	continue;

      /* find key in overflow Collision list */
      cnt -= COLLIDELINE;

      for(i = 0; i < cnt; i++){
	if(STAT) matchtable->HTcolcnt++;
	if(pr[i].key == sizekey){
	  hash[rcnt].hash = pr[i].hash;
	  if(STAT) matchtable->HTmatchcnt++;
	  if(STAT) matchtable->HTentrycnt += pr[i].n;
	  hash[rcnt++].n = pr[i].n;
	  break;
	}
      }
      if(DEBUG>=2 && QUICKHASH==2) assert(i < cnt);
    }

#if 0 // obsolete code using HashIndexV
    _mm_prefetch(CAST HashIndexV, _MM_HINT_T0);/* Sequential access fo HashIndexV[] of 16kb, so only single prefetch needed */
    _mm_prefetch(CAST sizekeyV, _MM_HINT_T0); /* Sequential access for sizekeyV[] of under 2K, so only single prefetch needed */

    for(register unsigned int S = 0; S < Slen; S++){
      register CHashIndexR *p = &HashIndexV[S];
      register unsigned int cnt = p->cnt;
      if(DEBUG>=2 && QUICKHASH==2) assert(cnt > 0);
      if(!cnt)
	continue;

      if(DEBUG>=2 && !(cnt <= CollideCntMax)){
	unsigned int hindex = rand1V[S];
	printf("hindex=%u:indexR[hindex].cnt=%u,cnt=%u,CollideCntMax=%u\n",hindex,indexR[hindex].cnt,cnt,CollideCntMax);
	fflush(stdout);
	assert(cnt <= CollideCntMax);
      }
    
      register long long sizekey = sizekeyV[S];

      if(STAT) matchtable->HTfindcnt++;
      register unsigned int i;    

      if(cnt <= COLLIDELINE){/* just look for key in indexR[hindex] */
	for(i = 0; i < cnt; i++){
	  if(STAT) matchtable->HTqcolcnt++;
	  if(p->key[i] == sizekey){
	    hash[rcnt].hash = p->hash[i];
	    if(STAT) matchtable->HTmatchcnt++;
	    if(STAT) matchtable->HTentrycnt += p->n[i];
	    hash[rcnt++].n = p->n[i];
	    break;
	  }
	}
	if(DEBUG>=2 && QUICKHASH==2 && !(i < cnt)){
#pragma omp critical
	  {
	    printf("S=%d/%d: i=%d,cnt=%d, key= 0x%llx(key1=0x%llx), sizekeyV[S]= 0x%llx(key1=0x%llx):\n",
		   S,Slen,i,cnt,(unsigned long long)key, (unsigned long long)KEYHASH1(key), (unsigned long long)sizekeyV[S], (unsigned long long)KEYHASH1(sizekeyV[S]));
	    for(unsigned int t = 0; t < cnt; t++)
	      printf("\t HashIndexV[S].key[%u] = 0x%llx\n",t, (unsigned long long)HashIndexV[S].key[t]);
	    fflush(stdout);
	    assert(i < cnt);
	  }
	}
	continue;
      }

      register CHashCollideR *pr = &CollideMemR[p->next];
      _mm_prefetch(CAST pr, _MM_HINT_NTA);    

      /* first look for key in indexR[hindex] */
      for(i = 0; i < COLLIDELINE; i++){
	if(STAT) matchtable->HTqcolcnt++;
	if(p->key[i] == sizekey){
	  hash[rcnt].hash = p->hash[i];
	  if(STAT) matchtable->HTmatchcnt++;
	  if(STAT) matchtable->HTentrycnt += p->n[i];
	  hash[rcnt++].n = p->n[i];
	  break;
	}
      }
      if(i < COLLIDELINE)
	continue;

      cnt -= COLLIDELINE;

      for(i = 0; i < cnt; i++){
	if(STAT) matchtable->HTcolcnt++;
	if(pr[i].key == sizekey){
	  hash[rcnt].hash = pr[i].hash;
	  if(STAT) matchtable->HTmatchcnt++;
	  if(STAT) matchtable->HTentrycnt += pr[i].n;
	  hash[rcnt++].n = pr[i].n;
	  break;
	}
      }
      if(DEBUG && QUICKHASH==2) assert(i < cnt);
    }
#endif // obsolete code using HashIndexV
  } else { // OLD CODE
    register unsigned int S, Smax;

    if(DEBUG && PD < 0) assert(PD >= 0);
#if PD > 0  
    unsigned int hindexV[PD2];/* Circular Lookahead buffer for hindex values of next PD iterations */
    _mm_prefetch(CAST hash, _MM_HINT_T0); /* Sequential access for hash[], so only single prefetch needed */
    Smax = min(PD,SizeLen);
    for(S = 0; S < Smax; S++){/* Prefetch loop */
      register long long sizekey = key + SizeDelta[S];
      register unsigned int hindex = rand1[sizekey & rhash_msk] ^ rand2[(sizekey >> RHASH_BITS) & rhash_msk];
      _mm_prefetch(CAST &indexR[hindex], _MM_HINT_NTA);
      hindexV[S] = hindex;
    }
#endif // PD > 0
    Smax = max(0, SizeLen - PD);
    for(S = 0; S < Smax; S++){/* main loop */
#if PD > 0
      register int SD = S + PD;/* look ahead iteration */
      register long long sizekeyD = key + SizeDelta[SD];
      register unsigned int hindexD = rand1[sizekeyD & rhash_msk] ^ rand2[(sizekeyD >> RHASH_BITS) & rhash_msk];
      _mm_prefetch(CAST &indexR[hindexD], _MM_HINT_NTA);
      hindexV[SD % PD2 ] = hindexD;

      register long long sizekey = key + SizeDelta[S];
      register unsigned int hindex = hindexV[S % PD2];
#else // PD <= 0
      register long long sizekey = key + SizeDelta[S];
      register unsigned int hindex = rand1[sizekey & rhash_msk] ^ rand2[(sizekey >> RHASH_BITS) & rhash_msk];
#endif
      if(STAT) matchtable->HTindcnt++;
      register CHashIndexR *p = &indexR[hindex];
      register unsigned int cnt = p->cnt;// Critical memory IO wait :  uses over 50% of total runtime for many hashtable usescase
      if(!cnt)
	continue;
      if(DEBUG>=2 && !(cnt <= CollideCntMax)){
	printf("hindex=%u:indexR[hindex].cnt=%u,cnt=%u,CollideCntMax=%u\n",hindex,indexR[hindex].cnt,cnt,CollideCntMax);
	fflush(stdout);
	assert(cnt <= CollideCntMax);
      }
    
      if(STAT) matchtable->HTfindcnt++;
      register unsigned int i;

      if(cnt <= COLLIDELINE){/* just look for key in indexR[hindex] */
	for(i = 0; i < cnt; i++){
	  if(STAT) matchtable->HTqcolcnt++;
	  if(p->key[i] == sizekey){
	    hash[rcnt].hash = p->hash[i];
	    if(STAT) matchtable->HTmatchcnt++;
	    if(STAT) matchtable->HTentrycnt += p->n[i];
	    hash[rcnt++].n = p->n[i];
	    break;
	  }
	}
	continue;
      }
      _mm_prefetch(CAST &CollideMemR[p->next], _MM_HINT_NTA);

      /* first look for key in indexR[hindex] */
      for(i = 0; i < COLLIDELINE; i++){
	if(STAT) matchtable->HTqcolcnt++;
	if(p->key[i] == sizekey){
	  hash[rcnt].hash = p->hash[i];
	  if(STAT) matchtable->HTmatchcnt++;
	  if(STAT) matchtable->HTentrycnt += p->n[i];
	  hash[rcnt++].n = p->n[i];
	  break;
	}
      }
      if(i < COLLIDELINE)
	continue;

      /* find key in overflow Collision list */
      register CHashCollideR *pr = &CollideMemR[p->next];
      cnt -= COLLIDELINE;

      for(i = 0; i < cnt; i++){
	if(STAT) matchtable->HTcolcnt++;
	if(pr[i].key == sizekey){
	  hash[rcnt].hash = pr[i].hash;
	  if(STAT) matchtable->HTmatchcnt++;
	  if(STAT) matchtable->HTentrycnt += pr[i].n;
	  hash[rcnt++].n = pr[i].n;
	  break;
	}
      }
    }

#if PD > 0
    for(; S < SizeLen; S++){/* Post-loop */
      if(STAT) matchtable->HTindcnt++;
      register long long sizekey = key + SizeDelta[S];
      register unsigned int hindex = hindexV[S % PD2];
      register CHashIndexR *p = &indexR[hindex];
      register unsigned int cnt = p->cnt;
      if(!cnt)
	continue;
    
      if(STAT) matchtable->HTfindcnt++;
      register unsigned int i;
    
      if(cnt <= COLLIDELINE){/* just look for key in indexR[hindex] */
	for(i = 0; i < cnt; i++){
	  if(STAT) matchtable->HTqcolcnt++;
	  if(p->key[i] == sizekey){
	    hash[rcnt].hash = p->hash[i];
	    if(STAT) matchtable->HTmatchcnt++;
	    if(STAT) matchtable->HTentrycnt += p->n[i];
	    hash[rcnt++].n = p->n[i];
	    break;
	  }
	}
	continue;
      }

      /* first look for key in indexR[hindex] */
      for(i = 0; i < COLLIDELINE; i++){
	if(STAT) matchtable->HTqcolcnt++;
	if(p->key[i] == sizekey){
	  hash[rcnt].hash = p->hash[i];
	  if(STAT) matchtable->HTmatchcnt++;
	  if(STAT) matchtable->HTentrycnt += p->n[i];
	  hash[rcnt++].n = p->n[i];
	  break;
	}
      }
      if(i < COLLIDELINE)
	continue;

      /* find key in overflow Collision list */
      register CHashCollideR *pr = &CollideMemR[p->next];
      cnt -= COLLIDELINE;

      for(i = 0; i < cnt; i++){
	if(STAT) matchtable->HTcolcnt++;
	if(pr[i].key == sizekey){
	  hash[rcnt].hash = pr[i].hash;
	  if(STAT) matchtable->HTmatchcnt++;
	  if(STAT) matchtable->HTentrycnt += pr[i].n;
	  hash[rcnt++].n = pr[i].n;
	  break;
	}
      }
    }
#endif // PD > 0 post-loop
  } // OLD CODE

  return rcnt;
}


#undef PD
#undef PD2

#define PD 4 /* TRY 1,2,4,8 */ /* Prefetch distance : Powers of 2 will be slightly faster: 1,2,4,8 */
#define PD2 (PD*2)

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::keyentry_hashloop(const int DISJOINT,const int PREFETCH, int i, int mapid2, int flip2, OFFSET_T offset2,
#if USE_HSSE==1 && USE_PFLOAT==1 && USE_SSE==1 && EXACT_SIZING && MAXR_WIN==8 && (HASHWIN == 5 || HASHWIN == 6)
				     __m128 vL_LB , __m128 vL_UB , __m128 vL_varInv , __m128 vL_win2,
				     __m128 vH_LB , __m128 vH_UB , __m128 vH_varInv , __m128 vH_win2, __m128 v_SDnormMax,
#endif
				     CHashEntry<OFFSET1_T> *hashBase, PFLOAT *win2, PFLOAT *LB, PFLOAT *UB, PFLOAT *varInv)
{
  register CHashEntry<OFFSET1_T> *phash;
#if PD > 0
  if(PREFETCH){
    phash = &hashBase[hashV[i+PD].hash];
    _mm_prefetch(CAST phash, _MM_HINT_T0);
    if(PD >= 4) _mm_prefetch(CAST &phash[2], _MM_HINT_T0);/* 2 entries per cache line */
    if(PD >= 6) _mm_prefetch(CAST &phash[4], _MM_HINT_NTA);/* 2 entries per cache line */
    if(PD >= 8) _mm_prefetch(CAST &phash[6], _MM_HINT_NTA);/* 2 entries per cache line */
  }
#endif // PD > 0

  register int n = hashV[i].n;
  if(DEBUG>=2) assert(n > 0);
  phash = &hashBase[hashV[i].hash];

  for(register int k = 0; k < n; k++){
#if PD > 0
    _mm_prefetch(CAST &phash[k+PD], _MM_HINT_T0);
#endif

    if(!DISJOINT)/* For single map sets (DisJoint==0) avoid self matches and duplicates of all matching pairs : we need ID1 < ID2, hence mapid1 > mapid2 */
      if(phash[k].mapid1 <= mapid2)
	continue; // break since phash[k].mapid1 is ordered in ascending order of mapid1

    if(EXACT_SIZING){ /* check if phash[k].(mapid1,site,err) is a valid sizing match by checking each interval AND computing the Chi-Square statistic (Gaussian Norm) */
      register PFLOAT *win1 = phash[k].win1; 
#if USE_HSSE==1 && USE_PFLOAT==1 && USE_SSE==1 && MAXR_WIN==8 && (HASHWIN == 5 || HASHWIN == 6)
      __m128 vL_win1 = _mm_load_ps(&win1[0]);
      __m128 vH_win1 = _mm_set1_ps(win1[4]);
      int bnderr, normerr;
      PFLOAT norm;
      if(DEBUG>=2){
	bnderr = normerr = 0;
	norm = PZERO;
	for(int w = HASHWIN; --w >= 0; ){
	  PFLOAT x1 = win1[w];
	  bnderr |= (x1 < LB[w]) || (UB[w] < x1);
	  PFLOAT err = x1 - win2[w];
	  norm += err * err * varInv[w];
	}
	normerr = (norm > SDnormMax) ? 1 : 0;
      }
      /*	  if(SIZING_BND && ((_mm_movemask_ps(_mm_cmplt_ps(vL_win1, vL_LB)) |
		  _mm_movemask_ps(_mm_cmplt_ps(vL_UB, vL_win1))) ||
		  win1[4] < LB[4] || UB[4] < win1[4])){*/
      if(SIZING_BND && (_mm_movemask_ps(_mm_cmplt_ps(vL_win1, vL_LB)) |
			_mm_movemask_ps(_mm_cmplt_ps(vL_UB, vL_win1)) |
			((_mm_movemask_ps(_mm_cmplt_ps(vH_win1, vH_LB)) |
			  _mm_movemask_ps(_mm_cmplt_ps(vH_UB, vH_win1))) & (HASHWIN == 5 ? 0x1 : 0x3)))){
	if(DEBUG>=2)  assert(bnderr);
	continue;
      }
      if(DEBUG>=2) assert(!bnderr);

      __m128 vL_err = _mm_sub_ps(vL_win1, vL_win2);
      __m128 vH_err = _mm_sub_ps(vH_win1, vH_win2);
      __m128 vL_norm = _mm_mul_ps(vL_varInv, _mm_mul_ps(vL_err,vL_err));
      __m128 vH_norm = _mm_mul_ps(vH_varInv, _mm_mul_ps(vH_err,vH_err));
      vL_norm = _mm_hadd_ps(vL_norm,vL_norm);/* horizontal add 4 floats -> 2 float */
      if(HASHWIN==6) vL_norm = _mm_add_ps(vL_norm,vH_norm);
      vL_norm = _mm_hadd_ps(vL_norm,vL_norm);/* horizontal add 2 floats -> 1 float */
      if(HASHWIN==5) vL_norm = _mm_add_ps(vL_norm,vH_norm);
      if(DEBUG>=2){
	PFLOAT v_norm;
	_mm_store_ss(&v_norm,vL_norm);
	if(fabs(v_norm-norm) > max(((PFLOAT)1.0),fabs(norm))*1e-6){
          #pragma omp critical
	  {
	    printf("Correct norm=%0.6f, v_norm=%0.6f,SDnormMax=%0.6f:\n",norm,v_norm,SDnormMax);
	    for(int w = 0; w < HASHWIN; w++)
	      printf("w=%d:win1[w]=%0.3f,win2[w]=%0.3f,varInv[w]=%0.8e\n",w,win1[w],win2[w],varInv[w]);
	    fflush(stdout);
	    assert(fabs(v_norm-norm) <= max(((PFLOAT)1.0),fabs(norm))*1e-6);
	  }
	}
      }
      if(_mm_ucomigt_ss(vL_norm, v_SDnormMax)) /* norm > SDnormMax */
	continue;
#else // Non-Vectorized case
      register PFLOAT norm = PZERO, bnderr;
      if(SIZING_BND) bnderr = PZERO;
      register int w = HASHWIN;//HashWin;
      for(; --w >= 0;){
	register PFLOAT x1 = win1[w];
	register PFLOAT bnd = max(LB[w]-x1,x1-UB[w]);
	bnderr = max(bnderr,bnd);
	//	    if(SIZING_BND) bnderr += max(PZERO, LB[w]-x1) + max(PZERO, x1 - UB[w]);
	register PFLOAT err = x1 - win2[w];
	norm += err * err * varInv[w];
      }
      if((SIZING_BND && bnderr > 0.0) || norm > SDnormMax)
	continue;
#endif // Non-Vectorized case
    } // if(EXACT_SIZING)

    /* will use offset to detect collisions (instead of site) */
    OFFSET1_T offset1 = flip2 ? phash[k].Roffset1 : phash[k].offset1;
    HitInsert(phash[k].mapid1, offset1, offset2, 0, 0);
  } // for(int k=0; k < n; k++)
}

/* NOTE : this is not the main CMatchTable::keyquery() : this one is called from the main function below */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::keyqueryB(const int DISJOINT, int hcnt, int mapid2, int flip2, OFFSET_T offset2,
#if USE_HSSE==1 && USE_PFLOAT==1 && USE_SSE==1 && EXACT_SIZING && MAXR_WIN==8 && (HASHWIN == 5 || HASHWIN == 6)
				     __m128 vL_LB , __m128 vL_UB , __m128 vL_varInv , __m128 vL_win2,
				     __m128 vH_LB , __m128 vH_UB , __m128 vH_varInv , __m128 vH_win2, __m128 v_SDnormMax,
#endif
				     CHashEntry<OFFSET1_T> *hashBase, PFLOAT *win2, PFLOAT *LB, PFLOAT *UB, PFLOAT *varInv)
{
  if(DEBUG && PD < 0) assert(PD >= 0);


  int imax;
#if PD > 0
  _mm_prefetch(CAST hashV, _MM_HINT_NTA);/* Sequential access for hashV[], so only single prefetch needed */

  imax = min(PD,hcnt);
  for(int i = 0; i < imax; i++){/* Prefetch loop */
    register CHashEntry<OFFSET1_T> *phash = &hashBase[hashV[i].hash];
    _mm_prefetch(CAST phash, _MM_HINT_NTA);
    if(PD >= 4) _mm_prefetch(CAST &phash[2], _MM_HINT_NTA);/* 2 entries per cache line */
    if(PD >= 6) _mm_prefetch(CAST &phash[4], _MM_HINT_NTA);/* 2 entries per cache line */
    if(PD >= 8) _mm_prefetch(CAST &phash[6], _MM_HINT_NTA);/* 2 entries per cache line */
  }
#endif // PD > 0    

  imax = hcnt - max(0,PD);
  int i = 0;
  for(; i < imax; i++)
    keyentry_hashloop(DISJOINT,1,i,mapid2,flip2,offset2,
#if USE_HSSE==1 && USE_PFLOAT==1 && USE_SSE==1 && EXACT_SIZING && MAXR_WIN==8 && (HASHWIN == 5 || HASHWIN == 6)
		      vL_LB,vL_UB,vL_varInv,vL_win2,vH_LB,vH_UB,vH_varInv,vH_win2,v_SDnormMax,
#endif
		      hashBase,win2,LB,UB,varInv);
  for(; i < hcnt; i++)/* post-prefetch loop */
    keyentry_hashloop(DISJOINT,0,i,mapid2,flip2,offset2,
#if USE_HSSE==1 && USE_PFLOAT==1 && USE_SSE==1 && EXACT_SIZING && MAXR_WIN==8 && (HASHWIN == 5 || HASHWIN == 6)
		      vL_LB,vL_UB,vL_varInv,vL_win2,vH_LB,vH_UB,vH_varInv,vH_win2,v_SDnormMax,
#endif
		      hashBase,win2,LB,UB,varInv);
}

/** Check all sizing variations of key and append hash matches to HitTable, removing duplicates
  * Optionally, check each hash match to make sure sizing errors are valid */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::keyquery(long long key, PFLOAT *win2, int mapid2, int flip2, FLOAT *X, int site2, unsigned int err2, OFFSET_T offset2, int scaleID)
{
  if(TRACE || TRACE2 || DEBUG>=2){/* old version of code without memory pipelining, but includes support for TRACE, TRACE2 and DEBUG>=2 */
    int TraceFlag;
    if(TRACE2)
      TraceFlag = 0;
    if(VERB>=2 || (TRACE2 >= 2 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2)){
      TraceFlag = 1;/* enable tracing of all size variations of key if it matches TraceQKey[] */
      long long ukey = key;
      printf("keyquery:mapid2=%d(id=%lld),flip2=%d,sc=%d,site2=%d,offset2=%d,err2=%u,mis2=%u,key=%lld",mapid2,qmap[mapid2]->id,flip2,scaleID,site2,offset2,err2 & MASK(8), err2 >> 8, ukey & MASK(SIZE_BITS));
      if(TraceQKey[0] != (ukey & MASK(SIZE_BITS)))
	TraceFlag = 0;
      for(unsigned int w = 1; w < HashWin; w++){
	ukey >>= SIZE_BITS;
	printf(",%lld",ukey & MASK(SIZE_BITS));
	if(TraceQKey[w] != (ukey & MASK(SIZE_BITS)))
	  TraceFlag = 0;
      }
      unsigned int E;
      if(!(E = err2)){
	printf(":X[%d..%d]=[%0.3f..%0.3f]=",site2,site2+HashWin,X[site2],X[site2+HashWin]);
	for(unsigned int w = 1; w <= HashWin; w++)
	  printf("%0.3f(%d),", X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
      } else if(E <= (int)MASK(4)){
	assert(NumErrQ >= 1 && 0 < E && E < HashWin+1);
	printf(":X[%d..%d]=[%0.3f..%0.3f]=",site2,site2+HashWin+1,X[site2],X[site2+HashWin+1]);
	for(unsigned int w = 1; w < E; w++)
	  printf("%0.3f(%d),", X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	printf("%0.3f+%0.3f(%d),", X[site2+E]-X[site2+E-1], X[site2+E+1]-X[site2+E], IntLen[(int)floor((X[site2+E+1]-X[site2+E-1])*QLenInv+0.5)]);
	for(unsigned int w = E+2; w <= HashWin+1; w++)
	  printf("%0.3f(%d),", X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
      } else if(E <= (int)MASK(8)){/* 2 errors */
	unsigned int E1 = E & MASK(4);
	unsigned int E2 = E >> 4;
	assert(NumErrQ >= 2 && 0 < E1 && E1 < E2 && E2 < HashWin+2);
	printf(":X[%d..%d]=[%0.3f..%0.3f]=",site2,site2+HashWin+2,X[site2],X[site2+HashWin+2]);
	for(unsigned int w = 1; w < E1; w++)
	  printf("%0.3f(%d),", X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	if(E2 == E1+1){
	  printf("%0.3f+%0.3f+%0.3f(%d),",X[site2+E1]-X[site2+E1-1],X[site2+E2]-X[site2+E1],X[site2+E2+1]-X[site2+E2],
		 IntLen[(int)floor((X[site2+E2+1]-X[site2+E1-1])*QLenInv+0.5)]);
	  for(unsigned int w = E2+2;w <= HashWin+2; w++)
	    printf("%0.3f(%d),",X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	} else {
	  printf("%0.3f+%0.3f(%d),",X[site2+E1] - X[site2+E1-1], X[site2+E1+1] - X[site2+E1], IntLen[(int)floor((X[site2+E1+1]-X[site2+E1-1])*QLenInv+0.5)]);
	  for(unsigned int w = E1+2; w < E2; w++)
	    printf("%0.3f(%d),", X[site2+w] - X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	  printf("%0.3f+%0.3f(%d),",X[site2+E2] - X[site2+E2-1], X[site2+E2+1] - X[site2+E2], IntLen[(int)floor((X[site2+E2+1]-X[site2+E2-1])*QLenInv+0.5)]);
	  for(unsigned int w = E2+2; w <= HashWin+2; w++)
	    printf("%0.3f(%d),",X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	}
      } else {/* handles cases with up to 2 site errors and 5 misresolved sites as a single case */
	int err = (((E)&MASK(4))?1:0) + (((E>>4)&MASK(4))?1:0) + (((E>>8)&MASK(4))?1:0) + (((E>>12)&MASK(4))?1:0) + (((E>>16)&MASK(4))?1:0) + (((E>>20)&MASK(4))?1:0) + (((E>>24)&MASK(4))?1:0);
	printf(":X[%d..%d]=[%0.3f..%0.3f]=",site2,site2+HashWin+err,X[site2],X[site2+HashWin+err]);
	for(unsigned int w = 0; w < HashWin; w++)
	  printf("%0.3f(%d),", win2[w], IntLen[(int)floor((win2[w]*QLenInv+0.5))]);
      }

      printf("\n");
      fflush(stdout);
    }

    /* precompute query window interval sizes and their variances for mapid2 (variance is approximate assuming the errors are small, to avoid divides in the inner loop) */
    PFLOAT LB[MAXR_WIN],UB[MAXR_WIN], varInv[MAXR_WIN];
    if(EXACT_SIZING){
      for(int w = HashWin; --w >= 0;){
	PFLOAT x = win2[w];
	varInv[w] = 1.0f/(QUADRATIC_VARIANCE ? fvarSF + fvarSD * x + fvarSR * x * x : fvarSF + fvarSD * x);
	if(SIZING_BND){
	  PFLOAT var = QUADRATIC_VARIANCE ? SFm + x * SDm + x * x * SRm : SFm + x * SDm;
	  UB[w] = x + sqrt(var);
	  if(QUADRATIC_VARIANCE)
	    LB[w] = (x - sqrt(var + VAR1) + SD1)*SRInv; /* Obtained by solving : x = LB[w] + sqrt(SFm + SDm * LB[w] + SRm * LB[w]*LB[w]) */
	  else
	    LB[w] = x - sqrt(var + VAR1) + SD1; /* Obtained by solving : x = LB[w] + sqrt(SFm + SDm * LB[w]) */
	  if(DEBUG>=2) assert(LB[w] < x);
	}
      }
      if(SIZING_BND && (VERB>=2 || (TRACE2 >= 2 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2))){
	printf("mapid2=%d(id=%lld),flip2=%d,sc=%d,site2=%d,err2=%u,mis2=%u:offset2=%d,window:\n",mapid2,qmap[mapid2]->id,flip2,scaleID,site2,err2 & MASK(8), err2 >> 8, offset2);
	for(unsigned int w = 0; w < HashWin; w++)
	  printf("win2[%d]=%0.3f:LB=%0.3f,UB=%0.3f\n", w, win2[w], LB[w], UB[w]);
	fflush(stdout);
      }
    }

    for(unsigned int i = 0; i < SizeLen; i++){
      long long sizekey = key + SizeDelta[i];
      CHashEntry<OFFSET1_T>* hash = 0;
      register int k = hashtable->hashfind(sizekey, hash, this);
      register CHashEntry<OFFSET1_T> *phash = hash;

      if(TRACE2 >= 2 && (TraceFlag || (k > 0 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2))){
	long long ukey = sizekey;
	printf("\tsize=%d/%d:skey=%lld", i, SizeLen, ukey & MASK(SIZE_BITS));
	for(unsigned int w = 1; w < HashWin; w++){
	  ukey >>= SIZE_BITS;
	  printf(",%lld", ukey & MASK(SIZE_BITS));
	}
	ukey = SizeDelta[i];
	for(unsigned int w = 0; w < HashWin; w++)
	  ukey += 1 << (w*SIZE_BITS);
	printf(",serr=%lld", (ukey  & MASK(SIZE_BITS)) - 1);
	for(unsigned int w = 1; w < HashWin; w++){
	  ukey >>= SIZE_BITS;
	  printf(",%lld", (ukey & MASK(SIZE_BITS)) - 1);
	}
	printf(":Hashtable hits=%d:\n",k);
	fflush(stdout);
      }

      if(VERB>=2){
	printf("key=%lld,sizekey=%lld:k=%d,phash=&EntryMem[%lu]\n",key,sizekey,k,phash-hashtable->EntryMem);
	fflush(stdout);
      }
    
      //   register int n = k;
      //       for(k = 0; k < n; k++){
      for(; --k >= 0;){
	if(DEBUG>=2 && k > 0 && (!DisJoint || HashInsertThreads <= 1)) assert(phash[k].mapid1 >= phash[k-1].mapid1);
	if(DEBUG>=2 && DisJoint && !HashMatrix && Gpairwise) assert(rmap[phash[k].mapid1]->id > qmap[mapid2]->id);

	/* For single map sets (DisJoint==0) avoid self matches and duplicates of all matching pairs : we need ID1 < ID2, hence mapid1 > mapid2 */
	if(!DisJoint && phash[k].mapid1 <= mapid2)
	  break; // break since phash[k].mapid1 is ordered in ascending order of mapid1

	register PFLOAT norm, bnderr;
	if(EXACT_SIZING){ /* check if phash[k].(mapid1,site,err) is a valid sizing match given (X, site2, err2) by checking each interval AND computing the Chi-Square statistic (Gaussian Norm) */
	  /* NOTE : this block is responsible for 45% of runtime */
	  register PFLOAT *win1 = phash[k].win1; 
	
	  norm = PZERO;
	  if(SIZING_BND)
	    bnderr = PZERO;
	  register int w = HashWin;
	  for(; --w >= 0;){
	    register PFLOAT x1 = win1[w];
	    if(SIZING_BND)
	      bnderr += max(PZERO, LB[w]-x1) + max(PZERO, x1 - UB[w]);
	    register PFLOAT err = x1 - win2[w];
	    norm += err * err * varInv[w];
	  }
	}

	int TraceFlag;
	if(VERB>=2 || TRACE || DEBUG >= 2)
	  TraceFlag = 0;

	if(!TRACE2 && EXACT_SIZING && ((SIZING_BND && bnderr > 0.0) || norm > SDnormMax))
	  continue;

	if(VERB>=2 || (TRACE && rmap[phash[k].mapid1]->id == TraceID1 && qmap[mapid2]->id == TraceID2 && flip2 == TraceOR2)){
	  unsigned int M = qmap[mapid2]->numsite[hcolor];
	  int offset = offset2 - (flip2 ? phash[k].Roffset1 : phash[k].offset1);
	  if(VERB>=2 || (TRACE && (TraceOffset == INT32_MIN || abs(offset - TraceOffset) <= OffsetSpread))){
	    TraceFlag = 1;
	    long long ukey = sizekey;

	    printf("HitInsert:mapid1=%d(id=%lld),mapid2=%d(id=%lld),flip2=%d,sc=%d,site2=%d,err2=%u,mis2=%u,size2=%d:norm=%0.4f,bnderr=%0.4f,offset1=%d,offset2=%d,offset=%d,skey=%lld",
		   phash[k].mapid1,rmap[phash[k].mapid1]->id,mapid2,qmap[mapid2]->id,flip2,scaleID, site2,err2 & MASK(8), err2 >> 8 ,
		   i,norm,bnderr,flip2 ? phash[k].Roffset1 : phash[k].offset1, offset2, offset,ukey & MASK(SIZE_BITS));
	    for(unsigned int w = 1; w < HashWin; w++){
	      ukey >>= SIZE_BITS;
	      printf(",%lld", ukey & MASK(SIZE_BITS));
	    }

	    ukey = SizeDelta[i];
	    for(unsigned int w = 0; w < HashWin; w++)
	      ukey += 1 << (w*SIZE_BITS);// make sure all ukey digits are 0,1 or 2 (instead of -1,0, or 1)
	    printf(",serr=%lld", (ukey  & MASK(SIZE_BITS)) - 1);
	    for(unsigned int w = 1; w < HashWin; w++){
	      ukey >>= SIZE_BITS;
	      printf(",%lld", (ukey & MASK(SIZE_BITS)) - 1);
	    }
	    unsigned int E;
	    printf("\n\twin1[]=%0.3f(%u)",phash[k].win1[0], IntLen[(int)floor(phash[k].win1[0]*QLenInv+0.5)]);
	    for(unsigned int w = 1; w < HashWin; w++)
	      printf(",%0.3f(%u)",phash[k].win1[w], IntLen[(int)floor(phash[k].win1[w]*QLenInv+0.5)]);
	    if(!(E = err2)){
	      printf(":X[%d..%d]=[%0.3f..%0.3f]=",site2,site2+HashWin,X[site2],X[site2+HashWin]);
	      for(unsigned int w = 1; w <= HashWin; w++)
		printf("%0.3f(%d),", X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	    } else if(E <= (int)MASK(4)){
	      assert(M >= site2+HashWin+1);
	      printf(":X[%d..%d]=[%0.3f..%0.3f]=",site2,site2+HashWin+1,X[site2],X[site2+HashWin+1]);
	      if(DEBUG && !(NumErrQ >= 1 && 0 < E && E < HashWin+1)){
		#pragma omp critical
		{
		  printf("\n\t NumErrQ=%d, mapid2=%d,flip2=%d,site2=%d,err2=%u,offset2=%d,M=%u: E=%u, HashWin=%u\n",NumErrQ,mapid2,flip2,site2,err2,offset2,M,E,HashWin);
		  fflush(stdout);
		  assert(NumErrQ >= 1 && 0 < E && E < HashWin+1);
		}
	      }
	      for(unsigned int w = 1; w < E; w++)
		printf("%0.3f(%d),", X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	      printf("%0.3f+%0.3f(%d),", X[site2+E]-X[site2+E-1], X[site2+E+1]-X[site2+E], IntLen[(int)floor((X[site2+E+1]-X[site2+E-1])*QLenInv+0.5)]);
	      for(unsigned int w = E+2; w <= HashWin+1; w++)
		printf("%0.3f(%d),", X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	    } else if(E <= (int)MASK(8)){ /* 2 errors */
	      unsigned int E1 = E & MASK(4);
	      unsigned int E2 = E >> 4;
	      assert(NumErrQ >= 2 && 0 < E1 && E1 < E2 && E2 < HashWin+2);
	      printf(":X[%d..%d]=[%0.3f..%0.3f]=",site2,site2+HashWin+2,X[site2],X[site2+HashWin+2]);
	      for(unsigned int w = 1; w < E1; w++)
		printf("%0.3f(%d),", X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	      if(E2 == E1+1){
		printf("%0.3f+%0.3f+%0.3f(%d),",X[site2+E1]-X[site2+E1-1],X[site2+E2]-X[site2+E1],X[site2+E2+1]-X[site2+E2],
		       IntLen[(int)floor((X[site2+E2+1]-X[site2+E1-1])*QLenInv+0.5)]);
		for(unsigned int w = E2+2;w <= HashWin+2; w++)
		  printf("%0.3f(%d),",X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	      } else {
		printf("%0.3f+%0.3f(%d),",X[site2+E1] - X[site2+E1-1], X[site2+E1+1] - X[site2+E1], IntLen[(int)floor((X[site2+E1+1]-X[site2+E1-1])*QLenInv+0.5)]);
		for(unsigned int w = E1+2; w < E2; w++)
		  printf("%0.3f(%d),", X[site2+w] - X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
		printf("%0.3f+%0.3f(%d),",X[site2+E2] - X[site2+E2-1], X[site2+E2+1] - X[site2+E2], IntLen[(int)floor((X[site2+E2+1]-X[site2+E2-1])*QLenInv+0.5)]);
		for(unsigned int w = E2+2; w <= HashWin+2; w++)
		  printf("%0.3f(%d),",X[site2+w]-X[site2+w-1], IntLen[(int)floor((X[site2+w]-X[site2+w-1])*QLenInv+0.5)]);
	      }
	    } else {/* handles case with up to 2 site errors and 5 misresolved sites as a single case */
	      int err = (((E)&MASK(4))?1:0) + (((E>>4)&MASK(4))?1:0) + (((E>>8)&MASK(4))?1:0) + (((E>>12)&MASK(4))?1:0) + (((E>>16)&MASK(4))?1:0) + (((E>>20)&MASK(4))?1:0) + (((E>>24)&MASK(4))?1:0);
	      printf(":X[%d..%d]=[%0.3f..%0.3f]=",site2,site2+HashWin+err,X[site2],X[site2+HashWin+err]);
	      for(unsigned int w = 0; w < HashWin; w++)
		printf("%0.3f(%d),", win2[w], IntLen[(int)floor((win2[w]*QLenInv+0.5))]);
	    }
	    printf("\n");
	    fflush(stdout);
	  }
	}

	if(TRACE2 && EXACT_SIZING && ((SIZING_BND && bnderr > 0.0) || norm > SDnormMax))
	  continue;

	/* Use offset to detect collisions (instead of site) : 
	   two sites in mapid1 could be close enough to result
	   in the same offset and producing two hits for the same window (mapid2,fip2,site2) */
	OFFSET1_T offset1 = flip2 ? phash[k].Roffset1 : phash[k].offset1;
	HitInsert(phash[k].mapid1, offset1, offset2, 0, TraceFlag);
      }
    }
    
  } else {/* !(TRACE || TRACE2 || DEBUG>=2) : new code with memory pipelining */
    /* precompute query window interval sizes and their variances for mapid2 (variance is approximate assuming the errors are small, to avoid divides in the inner loop) */
    PFLOAT *LB = (PFLOAT*) alloca(MAXR_WIN*sizeof(*LB));
    PFLOAT *UB = (PFLOAT*) alloca(MAXR_WIN*sizeof(*UB));
    PFLOAT *varInv = (PFLOAT*) alloca(MAXR_WIN*sizeof(*varInv));
    if(EXACT_SIZING){
      for(int w = HashWin; --w >= 0;){
	PFLOAT x = win2[w];
	varInv[w] = 1.0f/(QUADRATIC_VARIANCE ? fvarSF + fvarSD * x + fvarSR * x * x : fvarSF + fvarSD * x);
	if(SIZING_BND){
	  PFLOAT var = QUADRATIC_VARIANCE ? SFm + x * SDm + x * x * SRm : SFm + x * SDm;
	  UB[w] = x + sqrt(var);
	  if(QUADRATIC_VARIANCE)
	    LB[w] = (x - sqrt(var + VAR1) + SD1)*SRInv; /* Obtained by solving : x = LB[w] + sqrt(SFm + SDm * LB[w] + SRm * LB[w]*LB[w]) */
	  else
	    LB[w] = x - sqrt(var + VAR1) + SD1; /* Obtained by solving : x = LB[w] + sqrt(SFm + SDm * LB[w]) */
	}
      }
      if(SIZING_BND && VERB>=2){
	printf("mapid2=%d,flip2=%d,site2=%d,err2=%u,mis2=%u:offset2=%d,window:\n",mapid2,flip2,site2,err2 & MASK(8), err2 >> 8, offset2);
	for(unsigned int w = 0; w < HashWin; w++)
	  printf("win2[%d]=%0.3f:LB=%0.3f,UB=%0.3f\n", w, win2[w], LB[w], UB[w]);
	fflush(stdout);
      }
    }

#if USE_HSSE==1 && USE_PFLOAT==1 && USE_SSE==1 && EXACT_SIZING && MAXR_WIN==8 && (HASHWIN == 5 || HASHWIN == 6)
    __m128 vL_LB = _mm_load_ps(&LB[0]);
    __m128 vL_UB = _mm_load_ps(&UB[0]);
    __m128 vL_varInv = _mm_load_ps(&varInv[0]);
    __m128 vL_win2 = _mm_load_ps(&win2[0]);

    __m128 vH_LB = _mm_set1_ps(LB[4]);
    __m128 vH_UB = _mm_set1_ps(UB[4]);
    __m128 vH_varInv = _mm_set1_ps(varInv[4]);
    __m128 vH_win2 = _mm_set1_ps(win2[4]);

    __m128 v_SDnormMax = _mm_set1_ps(SDnormMax);
#endif

    register int hcnt = hashtable->hashfindV(key, hashV, this);
    register CHashEntry<OFFSET1_T> *hashBase = &hashtable->EntryMem[0];

    if(DisJoint)
      keyqueryB(1,hcnt,mapid2,flip2,offset2,
#if USE_HSSE==1 && USE_PFLOAT==1 && USE_SSE==1 && EXACT_SIZING && MAXR_WIN==8 && (HASHWIN == 5 || HASHWIN == 6)
		vL_LB,vL_UB,vL_varInv,vL_win2,vH_LB,vH_UB,vH_varInv,vH_win2,v_SDnormMax,
#endif
		hashBase,win2,LB,UB,varInv);
    else // if(!DisJoint)
      keyqueryB(0,hcnt,mapid2,flip2,offset2,
#if USE_HSSE==1 && USE_PFLOAT==1 && USE_SSE==1 && EXACT_SIZING && MAXR_WIN==8 && (HASHWIN == 5 || HASHWIN == 6)
		vL_LB,vL_UB,vL_varInv,vL_win2,vH_LB,vH_UB,vH_varInv,vH_win2,v_SDnormMax,
#endif
		hashBase,win2,LB,UB,varInv);

  }// !(TRACE || TRACE2 || DEBUG>=2) */
}

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
inline void CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::keyinsert(long long key, PFLOAT *win1, int mapid1, OFFSET1_T offset1, OFFSET1_T Roffset1, int site1, int err) 
{
  if(DEBUG>=2) assert(!sorted);
  if(DEBUG>=2)
    for(register int w = HashWin; --w >= 0;)
      assert(isfinite(win1[w]));

  if(VERB>=2 || (TRACE2 && rmap[mapid1]->id == TraceID1)){
    long long ukey = key;
    printf("keyinsert:mapid1=%d(id=%lld),site1=%d,offset1=%d,%d,err1=%d:key=%lld",
	   mapid1,rmap[mapid1]->id,site1,offset1,Roffset1,err, ukey & MASK(SIZE_BITS));
    for(unsigned int w = 1; w < HashWin; w++){
      ukey >>= SIZE_BITS;
      printf(",%lld", ukey & MASK(SIZE_BITS));
    }
    unsigned int E;
    FLOAT *Y = YY[mapid1];
    if(!(E = err)){
      printf(":Y[%d..%d]=[%0.3f..%0.3f]",site1,site1+HashWin,Y[site1],Y[site1+HashWin]);
      for(unsigned int w = 1; w <= HashWin; w++)
	printf(",%0.3f(%d)", Y[site1+w] - Y[site1+w-1], IntLen[(int)floor((Y[site1+w]-Y[site1+w-1])*QLenInv+0.5)]);
    } else {
      unsigned int N = rmap[mapid1]->numsite[hcolor];
      if(E <= (int)MASK(4)){
	assert(N >= site1+HashWin+1);
	assert(NumErrI >= 1 && 0 < E && E < HashWin+1);
	printf(":Y[%d..%d]=[%0.3f..%0.3f]=",site1,site1+HashWin+1,Y[site1],Y[site1+HashWin+1]);
	for(unsigned int w = 1; w < E; w++)
	  printf("%0.3f(%d),", Y[site1+w] - Y[site1+w-1], IntLen[(int)floor((Y[site1+w]-Y[site1+w-1])*QLenInv+0.5)]);
	printf("%0.3f+%0.3f(%d),", Y[site1+E] - Y[site1+E-1], Y[site1+E+1] - Y[site1+E], IntLen[(int)floor((Y[site1+E+1]-Y[site1+E-1])*QLenInv+0.5)]);	    
	for(unsigned int w = E+2; w <= HashWin+1; w++)
	  printf("%0.3f(%d),", Y[site1+w] - Y[site1+w-1], IntLen[(int)floor((Y[site1+w]-Y[site1+w-1])*QLenInv+0.5)]);
      } else {/* E > MASK(4) : 2 Errors */
	unsigned int E1 = E & MASK(4);
	unsigned int E2 = E >> 4;
	assert(N >= site1+HashWin+2);
	assert(NumErrI >= 2 && 0 < E1 && E1 < E2 && E2 < HashWin+2);
	printf(":Y[%d..%d]=[%0.3f..%0.3f]=",site1,site1+HashWin+2,Y[site1],Y[site1+HashWin+2]);
	for(unsigned int w = 1; w < E1; w++)
	  printf("%0.3f(%d),", Y[site1+w] - Y[site1+w-1], IntLen[(int)floor((Y[site1+w]-Y[site1+w-1])*QLenInv+0.5)]);
	if(E2 == E1+1){
	  printf("%0.3f+%0.3f+%0.3f(%d),",Y[site1+E1]-Y[site1+E1-1],Y[site1+E2]-Y[site1+E1],Y[site1+E2+1]-Y[site1+E2],
		 IntLen[(int)floor((Y[site1+E2+1]-Y[site1+E1-1])*QLenInv+0.5)]);
	  for(unsigned int w = E2+2;w <= HashWin+2; w++)
	    printf("%0.3f(%d),",Y[site1+w]-Y[site1+w-1], IntLen[(int)floor((Y[site1+w]-Y[site1+w-1])*QLenInv+0.5)]);
	} else {
	  printf("%0.3f+%0.3f(%d),",Y[site1+E1] - Y[site1+E1-1], Y[site1+E1+1] - Y[site1+E1], IntLen[(int)floor((Y[site1+E1+1]-Y[site1+E1-1])*QLenInv+0.5)]);
	  for(unsigned int w = E1+2; w < E2; w++)
	    printf("%0.3f(%d),", Y[site1+w] - Y[site1+w-1], IntLen[(int)floor((Y[site1+w]-Y[site1+w-1])*QLenInv+0.5)]);
	  printf("%0.3f+%0.3f(%d),",Y[site1+E2] - Y[site1+E2-1], Y[site1+E2+1] - Y[site1+E2], IntLen[(int)floor((Y[site1+E2+1]-Y[site1+E2-1])*QLenInv+0.5)]);
	  for(unsigned int w = E2+2; w <= HashWin+2; w++)
	    printf("%0.3f(%d),",Y[site1+w]-Y[site1+w-1], IntLen[(int)floor((Y[site1+w]-Y[site1+w-1])*QLenInv+0.5)]);
	}
      }
    }

    printf("\n");
    fflush(stdout);
  }

  long long hindex = hashindex(key);
  unsigned int *ppindex = &index[hindex];
  unsigned int pindex = *ppindex;
  if(!pindex){
    *ppindex = pindex = CollideBlockAlloc(MINMEM);
    CHashCollide *p = &CollideBlock[pindex];
    p->n = 1;
    unsigned int hash = EntryBlockAlloc(MINMEM2);
    CHashEntry<OFFSET1_T> *phash = &EntryBlock[hash];
    phash->init(win1,mapid1,offset1,Roffset1,site1,err);
    p[1].init(1,key,hash);
    if(DEBUG>=2){
      for(int k = p[1].n;--k >= 0;){
	register CHashEntry<OFFSET1_T> *phash = &EntryBlock[p[1].hash + k];
	for(register int w = HashWin; --w >= 0;)
	  assert(isfinite(phash->win1[w]));
      }
    }
    return;
  }
  CHashCollide *p = &CollideBlock[pindex];

  /* check if key is already present in Collision list */
  unsigned int n = p->n;
  unsigned int i = 1;
  for(; i <= n; i++)
    if(p[i].key == key)
      break;
  if(i <= n){/* found another hashentry with same key */
    if(DEBUG>=2) assert(p[i].key == key);
    if(DEBUG>=2) assert(p[i].n >= 1);
    if(DEBUG>=2) assert(p[i].hash != 0);

    p += i;

    unsigned int pn = p->n;
    if(DEBUG>=2) assert(pn >= 1);
    if(pn >= MINMEM2 && !(pn & (pn-1))){/* double size of p->hash[] array */
      if(pn > MASK(31)){
	printf("hashinsert: hashentry array length=%d cannot be doubled with 32 bit unsigned index\n",pn);
	fflush(stdout);	exit(1);
      }
      unsigned int orighash = p->hash;
      unsigned int newhash = EntryBlockAlloc(pn*2);
      memcpy(&EntryBlock[newhash],&EntryBlock[orighash],sizeof(CHashEntry<OFFSET1_T>)*pn);
      /* NOTE : no garbage collection so no need to "delete" original block of memory */
      p->hash = newhash;
    }
    CHashEntry<OFFSET1_T> *phash = &EntryBlock[p->hash + p->n++];
    phash->init(win1,mapid1,offset1,Roffset1,site1,err);
    if(DEBUG>=2){
      assert(pindex == index[hindex]);
      assert(p->n == CollideBlock[pindex + i].n);
      for(int k = p->n; --k >= 0;){
	register CHashEntry<OFFSET1_T> *phash = &EntryBlock[p->hash + k];
	for(register int w = HashWin; --w >= 0;){
	  assert(isfinite(phash->win1[w]));
	}
      }
    }
    return;
  }

  /* Create new hashentry array and add it to Collision list */
  if(++n >= MINMEM && !(n & (n-1))){/* double size of array p to 2*n */
    if(n > MASK(31)){
      printf("hashinsert: hash Collision array length=%d cannot be doubled with 32 bit unsigned index\n", n);
      fflush(stdout); exit(1);
    }
    unsigned int newpindex = CollideBlockAlloc(n*2);/* Note p is now invalid, since CollideBlock[] may have been reallocated */
    CHashCollide *newp = &CollideBlock[newpindex];
    memcpy(newp, &CollideBlock[pindex], n * sizeof(CHashCollide));
    *ppindex = pindex = newpindex;
    p = newp;
  }
  p->n = n;

  /* Initialize Collision list entry p[n] */
  unsigned int hash = EntryBlockAlloc(MINMEM2);
  CHashEntry<OFFSET1_T> *phash = &EntryBlock[hash];
  phash->init(win1,mapid1,offset1,Roffset1,site1,err);
  p[n].init(1,key,hash);
  if(DEBUG>=2){
    for(int k = p[n].n; --k >= 0;){
      register CHashEntry<OFFSET1_T> *phash = &EntryBlock[p[n].hash + k];
      for(register int w = HashWin; --w >= 0;)
	assert(isfinite(phash->win1[w]));
    }
  }
}

static void PrintStatistics(){/* display HashTable query performance statistics */
  if(STAT && gHTfindcnt > 0){/* display performance statistics */
    if(QUICKHASH && hashkeys)
      printf("\t HashTable:probes=%lld,%lld,hits=%lld,inline collisions=%lld(%0.3f/hit),collisions=%lld(%0.3f/hit),matches=%lld(%0.3f/hit),total entries=%lld(%.3f/match)\n",
	     gHTindcnt, gHTqindcnt, gHTfindcnt, gHTqcolcnt-gHTmatchcnt, ((double)(gHTqcolcnt-gHTmatchcnt))/gHTfindcnt,
	     gHTcolcnt, (double)(gHTcolcnt)/gHTfindcnt, gHTmatchcnt, ((double)gHTmatchcnt)/gHTfindcnt, gHTentrycnt, ((double)gHTentrycnt)/gHTmatchcnt);
    else
      printf("\t HashTable:probes=%lld,hits=%lld,inline collisions=%lld(%0.3f/hit),collisions=%lld(%0.3f/hit),matches=%lld(%0.3f/hit),total entries=%lld(%.3f/match)\n",
	     gHTindcnt, gHTfindcnt, gHTqcolcnt-gHTmatchcnt, ((double)(gHTqcolcnt-gHTmatchcnt))/gHTfindcnt,
	     gHTcolcnt, (double)(gHTcolcnt)/gHTfindcnt, gHTmatchcnt, ((double)gHTmatchcnt)/gHTfindcnt, gHTentrycnt, ((double)gHTentrycnt)/gHTmatchcnt);
    printf("\t HitTable:probes=%lld,collisions=%lld(%0.3f/probe),duplicates=%lld\n\t MatchTable:probes=%lld+%lld+%lld InLine collisions=%lld(%0.3f/probe), Other collisions=%lld(%0.3f/probe), inserts=%lld\n",
	   gHTcnt, gHTcol, ((double)gHTcol)/gHTcnt, gHTdup, gMTcnt, gMTcnt2, gMTcnt3, gMTqcol, ((double)gMTqcol)/(gMTcnt+gMTcnt2+gMTcnt3), gMTcol, ((double)gMTcol)/(gMTcnt+gMTcnt2+gMTcnt3), gMTicnt);
    fflush(stdout);
  }
}

/* compute hash key for quantized window delta[deltawin[0..HashWin-1]] */
static inline long long hashkey(int *deltawin,int *delta)
{
  long long key = 0;
  /* NOTE : delta[] values can be -ve hence we need to used * and + instead of << and | */
  for(int i = HashWin; --i >= 0;){
    key *= SIZES;
    key += delta[deltawin[i]];
  }
  return key;
}

/* NOTE : rmap[] is the set of maps inserted into the hashtable, qmap[] the set of maps used to probe the hashtable : 
   Typically rmap[] are the query maps (molecules) and qmap[] the reference maps (contigs)
   
   If splitcnt > 1 : only rmap[roffset .. roffset + numrmaps - 1] is used
 */

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void hash_init(Cmap *rmap[], int numrmaps, Cmap *qmap[], int numqmaps, int pairwise, char *prefix)
{
  if(DEBUG && splitcnt <= 1) assert(roffset == 0);
  if(DEBUG && splitcnt > 1) assert(rmap != qmap);

#if HASH_STREAM < 2
  printf("HASH_STREAM=%d no longer supported (see hash.h)\n",HASH_STREAM);
  fflush(stdout);exit(1);
#endif

#if NEWHASH < 2
  printf("NEWHASH=%d no longer supported\n", NEWHASH);
  fflush(stdout); exit(1);
#endif

  if(NumScaleFactor >= 256){
    printf("-ScaleDelta %0.3f %d : number of scaling factors is limited to 127\n",ScaleDeltaSize,ScaleDelta);
    fflush(stdout);exit(1);
  }

  if(TRACE>=2){
    printf("TRACE=%d no longer supported\n", TRACE);
    fflush(stdout); exit(1);
  }

  if(HashTxt && DeferMerge){
    printf("Cannot combine options -defermerge and -hashtxt with -hashgen\n");
    fflush(stdout); exit(1);
  }
  if(hash_filename && DeferMerge){
    printf("Cannot combine options -defermerge and -hash with -hashgen\n");
    fflush(stdout); exit(1);
  }

  if(NumErrQ > 2){
    printf("NumErrP=%d: Values over 2 not implemented\n", NumErrQ);
    fflush(stdout); exit(1);
  }
  if(NumErrI > 2){
    printf("ERROR:NumErrI=%d: Values over 2 not implemented\n", NumErrI);
    fflush(stdout); exit(1);
  }
  if(NumResQ > 5){
    printf("WARNING:NumResP=%d: NumResQ over 5 not implemented\n", NumResQ);
    fflush(stdout);
  }
  if(NumErrQ + NumResQ > 5){
    printf("WARNING:NumErrP=%d,NumResP=%d: Total errors over 5 per window not implemented\n",NumErrQ,NumResQ);
    fflush(stdout);
  }
  if(HashWin + max(NumResQ,NumErrQ-1) > 15){
    printf("ERROR:NumErrP=%d,NumResP=%d,HashWin=%u: HashWin cannot be larger than 15-max(NumResP,NumErrP-1)\n",NumErrQ,NumResQ,HashWin);
    fflush(stdout); exit(1);
  }
  if(sizeof(OFFSET_T) > sizeof(int)){
    printf("sizeof(OFFSET_T)=%lu must be less than or equal to sizeof(int)=%lu\n",sizeof(OFFSET_T),sizeof(int));
    fflush(stdout); exit(1);
  }

  if(QUICKSTAT)
    hashfindV_calls = hashfindV_lines = hashfindV_max = 0;

  SDnormMax = SDrms * SDrms * HashWin;
  assert(sizeof(int)==4);
  assert(sizeof(long long)==8);
  if(!(sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) == 64)){
    printf("sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>)= %lu, MATCHLINE=%d,PADDING=%lu\n",sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>), MATCHLINE(OFFSET_T,RANGE,HASH_SCALE), MPADDING(OFFSET_T,RANGE,HASH_SCALE));
    fflush(stdout);
    assert(sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>)==64);
  }
  assert(sizeof(CHashIndexR)==64);
  
  if(VERB/* HERE >=2 */){
    printf("sizeof(CMatchLine)=%lu, sizeof(CBestMatch)=%lu, sizeof(CHitEntry)=%lu, sizeof(CHashMatch)=%lu,sizeof(CHashCollide)=%lu,sizeof(CHashIndexR)=%lu\n",
	   sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>), sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>),sizeof(CHitEntry<OFFSET1_T>),sizeof(CHashMatch),sizeof(CHashCollide),sizeof(CHashIndexR));
    fflush(stdout);
  }

  if(VERB){
    printf("HashTable initialization:HashWin=%u,HashScore=%d,SDmax=%0.2f,SDrms=%0.2f(norm=%0.4f),RelErr=%0.3f,OffsetKB=%0.3f,NumErrI=%d,NumErrQ=%d,NumResQ=%d\n",
	   HashWin,HashScore,SDmax,SDrms,SDnormMax,RelErr,OffsetKB,NumErrI,NumErrQ,NumResQ);
    //    printf("sleeping for 60 seconds\n");
    printf("\t sf= %0.6f, sd= %0.6f, sr= %0.6f\n", SF[hcolor], SD[hcolor], SR[hcolor]);
    fflush(stdout);
    //    sleep(60);
  }

  if(HASHWIN && HashWin != HASHWIN){
    printf("-hashgen window size must be %d (change HASHWIN in hash.h to change)\n",HASHWIN);
    fflush(stdout); exit(1);
  }
  if(HashWin > (64/SIZE_BITS)){
    printf("hashtable cannot use a window greater than 64/SIZE_BITS = %d (HashWin=%u, SIZE_BITS=%d)\n",64/SIZE_BITS,HashWin, SIZE_BITS);
    fflush(stdout); exit(1);
  }
  if(HashWin > MAX_WIN(OFFSET1_T)){
    printf("HashWin=%u: 1st arg of -hashgen cannot exceed %d (sizeof(OFFSET1_T)=%lu)\n", HashWin, MAX_WIN(OFFSET1_T),sizeof(OFFSET1_T));
    fflush(stdout); exit(1);
  }
  if(RHASH_BITS*3 < SIZE_BITS * HashWin){
    printf("hashtable cannot use a window greater than %d (HashWin=%u)\n", (RHASH_BITS*3)/SIZE_BITS,HashWin);
    fflush(stdout); exit(1);
  }
  /*  if(NumErrQ > 1){
    printf("EXACT_SIZING not implemented for numerr > 1 (NumErrQ=%d)\n",NumErrQ);
    fflush(stdout); exit(1);
    }*/

  Gpairwise = pairwise;
  if(!pairwise){
    if(OffsetSpread >= 2)
      for(int i = 2; i <= OffsetSpread; i++)
	OFFSET_WT[i] = min(OFFSET_WT[0],OFFSET_WT[1]);
    else
      for(int i = 2; i > OffsetSpread; i--)
	OFFSET_WT[i] = 0;
  }

  /* initialize global statistics counters */
  TotalMatchCnt = 0;
  if(STAT) gHTindcnt = gHTfindcnt = gHTcolcnt = gHTmatchcnt = gHTentrycnt = 0;
  if(STAT) gMTcnt = gMTcnt2 = gMTcnt3 = gMTqcol = gMTcol = gMTicnt = gHTcnt = gHTcol = 0;

  MaxHit = (MAXHIT_COVERAGE && hashMaxCov > 0) ? hashMaxCov : MASK(16);

  if(DEBUG>=1+RELEASE && (FinalSort==4 || numrefmaps > Hapnumrefmaps)){  // make sure original qmap[] and rmap[] are in mapid order
    for(int i = 0; i < numrmaps; i++)
      assert(rmap[i + roffset]->mapid == i + roffset);

    for(int i = 0; i < numqmaps; i++)
      assert(qmap[i]->mapid == i);
  }

  if(VERB>=2){
    printf("rmap[%d]->mapid=%d, rmap[%d]->mapid=%d\n",roffset,rmap[roffset]->mapid, roffset+numrmaps-1,rmap[roffset+numrmaps-1]->mapid);
    printf("qmap[0]->mapid=%d, qmap[%d]->mapid=%d\n",qmap[0]->mapid, numqmaps-1,qmap[numqmaps-1]->mapid);
    fflush(stdout);
  }

  if(VERB>=2)
    for(int i = 0; i < numrmaps; i++)
      rmap[i + roffset]->mapid = i + roffset;

  /* sort maps by id in ascending order and renumber mapids */
  qsort(&rmap[roffset],numrmaps,sizeof(Cmap *),(intcmp*)CmapIdInc);

  /* check for duplicate or out-of-order map ids and count total number of sites in rmaps */
  long long rsites = 0;
  for(int i = 0; i < numrmaps; i++){
    rsites += rmap[roffset + i]->numsite[hcolor];
    if(DEBUG && i > 0 && !(rmap[roffset + i-1]->id < rmap[roffset + i]->id)){
      printf("Duplicate Map ids found in Map set (roffset=%d,numrmaps=%d):\n", roffset,numrmaps);
      for(int j = i-1; j <= i; j++)
	printf("rmap[%d]->id = %lld, mapid= %d\n", roffset + j, rmap[roffset + j]->id, rmap[roffset + j]->mapid);
      fflush(stdout); exit(1);
    }
  }

  if(DEBUG && !pairwise) assert(qmap != rmap);

  if(Cfirst && CfirstPairs == 1){
    printf("Cfirst=%d, CfirstPairs= 1 : Not supported\n",Cfirst);
    fflush(stdout);exit(1);
  }

  DisJoint = 0;

  if(qmap == rmap){
    if(numqmaps != numrmaps){
      printf("hash_init: rmap and qmap must be identical or non-overlapping sets of maps (numqmaps=%d,numrmaps=%d)\n",numqmaps,numrmaps);
      fflush(stdout); exit(1);
    }
    if(DEBUG) assert(roffset==0);
    for(int i = 0; i < numrmaps; i++)
      rmap[i]->mapid = i;
    if(VERB && TRACE){
      printf("hash_init: rmap == qmap : numrmaps=%d:\n",numrmaps);
      fflush(stdout);
    }
  } else if(numqmaps > 0){ /* qmap != rmap */
    if(VERB){
      printf("hash_init: Looking for matches between non-overlapping sets of maps (numqmaps=%d,numrmaps=%d,roffset=%d,numrefmaps=%d,Hapnumrefmaps=%d,FinalSort=%d,SplitHmap=%d)\n",
	     numqmaps,numrmaps,roffset,numrefmaps,Hapnumrefmaps,FinalSort,SplitHmap);
      fflush(stdout);
    }
    qsort(qmap,numqmaps,sizeof(Cmap *),(intcmp*)CmapIdInc);
    if(DEBUG){/* check for duplicate or out-of-order map ids */
      for(int i = 1; i < numqmaps; i++){
	if(!(qmap[i-1]->id < qmap[i]->id)){
	  printf("Duplicate Map ids found in Query Map set:\n");
	  for(int j = i-1; j <= i; j++)
	    printf("qmap[%d]->id = %lld\n", j, qmap[j]->id);
	  fflush(stdout);exit(1);
	}
      }
    }

    if(VERB>=2){
      printf("After sorting by id:\nqmap[0,numqmaps-1]->id= %lld, %lld, rmap[roffset,roffset+numrmaps-1]->id= %lld, %lld\n",
	     qmap[0]->id,qmap[numqmaps-1]->id,rmap[roffset]->id,rmap[roffset+numrmaps-1]->id);
      fflush(stdout);
    }

    DisJoint = 1;/* default if qmap != rmap */

    if(Cfirst && CfirstPairs && pairwise)/* PairMerge may use overlapped (or non-overlapped) subsets of the same map set */
      DisJoint = (qmap[numqmaps-1]->id < rmap[0]->id) ? 1 : 0;

    if(FinalSort==4 || numrefmaps > Hapnumrefmaps){  // re-sort qmap[] and rmap[] in original mapid order
      if(VERB/* HERE >=2 */){
        printf("re-sorting qmap[0..%d] and rmap[0..%d] in original mapid order (DisJoint=%d)\n",numqmaps-1,numrmaps-1,DisJoint);
        fflush(stdout);
      }

      qsort(&rmap[roffset],numrmaps,sizeof(Cmap *),(intcmp*)CmapMapidInc);
      if(rmap != qmap)
	qsort(qmap,numqmaps,sizeof(Cmap *),(intcmp*)CmapMapidInc);	

      if(VERB>=2){
	printf("rmap[%d]->mapid=%d, rmap[%d]->mapid=%d\n",roffset,rmap[roffset]->mapid, roffset+numrmaps-1,rmap[roffset+numrmaps-1]->mapid);
	printf("qmap[0]->mapid=%d, qmap[%d]->mapid=%d\n",qmap[0]->mapid, numqmaps-1,qmap[numqmaps-1]->mapid);
	fflush(stdout);
      }

      if(VERB>=2){
	for(int t = 0; t < numrefmaps; t++){
	  printf("refmap[%d]= %p: id=%lld, Allele=%d, contig=%p\n", t, refmap[t], refmap[t]->id, refmap[t]->Allele, refmap[t]->contig);
	  if(refmap[t]->contig && refmap[t]->contig->contig)
	    printf("\t contig->contig[0].id= %lld, contig->contig[1].id= %lld\n", refmap[t]->contig->contig[0].id, refmap[t]->contig->contig[1].id);
	}
	fflush(stdout);
      }
    }

    /* renumber mapids in qmap[] and rmap[] */
    if(qmap != rmap){
      if(DEBUG && Cfirst && CfirstPairs == 2) assert(numrmaps <= numqmaps);/* since rmap[] should be a subset of qmap[] consisting of the last numrmaps entries of qmap[] */

      if(!DisJoint){/* overlapping mapsets */
	if(DEBUG) assert(roffset == 0);
	for(int i = 0; i < numrmaps; i++)
	  rmap[i]->mapid = -1;
      }

      for(int i = 0; i < numqmaps; i++)
	qmap[i]->mapid = i;

      if(DisJoint){
	if(DEBUG>=2){
	  for(int i = 0; i < numrmaps; i++){
	    for(int j = i+1; j < numrmaps; j++){
	      if(&rmap[roffset+i]->mapid == &rmap[roffset+j]->mapid || rmap[roffset+i] == rmap[roffset+j]){
		printf("i=%d,j=%d,numrmaps=%d:rmap[i]=%p,rmap[j]=%p,&rmap[i]->mapid=%p,&rmap[j]->mapid=%p\n",
		       i,j,numrmaps,rmap[roffset+i],rmap[roffset+j],&rmap[roffset+i]->mapid,&rmap[roffset+j]->mapid);
		fflush(stdout);
		assert(&rmap[roffset+i]->mapid != &rmap[roffset+j]->mapid);
	      }
	    }
	  }
	}

	if(VERB>=2){
	  printf("Before reseting mapid values rmap[%d..%d]:\n",roffset,roffset+numrmaps-1);
	  for(int i = 0; i < 10; i++)
	    printf("i=%d:rmap[roffset+i]->mapid=%d at %p\n",i,rmap[roffset+i]->mapid,&rmap[roffset+i]->mapid);
	  fflush(stdout);
	}

	for(int i = 0; i < numrmaps; i++)
	  rmap[roffset + i]->mapid = i + roffset;

	if(VERB/* HERE >=2 */ && numrmaps > 0){
	  printf("Reset mapid values of rmap[%d..%d]: rmap[%d]->mapid=%d, rmap[%d]->mapid=%d\n",
		 roffset,roffset+numrmaps-1,roffset,rmap[roffset]->mapid,roffset+numrmaps-1,rmap[roffset+numrmaps-1]->mapid);
	  fflush(stdout);
	}
	if(DEBUG>=2){
	  int cnt = 0;
	  for(int i = 0; i < numrmaps; i++){
	    if(rmap[roffset+i]->mapid != roffset+i){
	      printf("i=%d,roffset=%d:rmap[roffset+i]->mapid=%d at %p\n",i,roffset,rmap[roffset+i]->mapid, &rmap[roffset+i]->mapid);
	      fflush(stdout);
	      //	      assert(rmap[roffset+i]->mapid == roffset+i);
	      cnt++;
	    }
	  }
	  if(cnt > 0){
	    printf("cnt=%d errors\n",cnt);
	    fflush(stdout);
	    assert(cnt <= 0);
	  }
	}
      }
      else {// overlapping mapsets
        if(DEBUG)assert(roffset == 0);

	int nextid = numqmaps, i = 0;
	for(; i < numrmaps; i++)
	  if(rmap[i]->mapid < 0)
	    break;
	if(DEBUG && Cfirst && CfirstPairs == 2) assert(i >= numrmaps);/* since rmap[] should be a subset of qmap[] consisting of the last numrmaps entries of qmap[] */

	for(;i < numrmaps; i++){
	  if(DEBUG) assert(rmap[i]->mapid < 0);
	  rmap[i]->mapid = nextid++;
	}
      }
    } else {// qmap == rmap
      if(DEBUG) assert(numrmaps == numqmaps);
      for(int i = 0; i < numrmaps; i++)
	rmap[i]->mapid = i;
    }

    if(SplitHmap){/* also renumber the Allele fields (index or mapid of other Allele) */
      Cid2mapid *id2mapid = new Cid2mapid[numrefmaps];
      for(int k = 0; k < numrefmaps; k++){
	Cid2mapid *p = &id2mapid[k];
	p->id = refmap[k]->id;
	p->mapid = k;
      }
      qsort(id2mapid,numrefmaps,sizeof(Cid2mapid),(intcmp *)idinc);

      long long idoffset = MaxContigID;// WAS (HapSitePvalue > 0.0 || HapIndelPvalue > 0.0) ? MaxContigID : 0;
      for(int k = 0; k < Hapnumrefmaps; k++){
	if(refmap[k]->contig){
	  long long id = refmap[k]->id;
	  if(DEBUG) assert(refmap[k]->contig->contig[0].id == id);

	  /* first locate other allele */
	  long long id2 = id + idoffset;
	  if(DEBUG && !(refmap[k]->contig->contig[1].id == id2)){
	    printf("k=%d,Hapnumrefmaps=%d,numrefmaps=%d:id=%lld,idoffset=%lld,id2=%lld,refmap[k]->contig->contig[0,1].id=%lld,%lld\n",
		   k,Hapnumrefmaps,numrefmaps,id,idoffset,id2,refmap[k]->contig->contig[0].id,refmap[k]->contig->contig[1].id);
	    fflush(stdout);
	    assert(refmap[k]->contig->contig[1].id == id2);
	  }

	  int k2 = findid(id2,id2mapid,numrefmaps);
	  if(k2 < 0){
	    printf("ERROR:Cannot find id2=%lld in refmap[0..%d] as matching Allele2 of refmap[%d]->id= %lld (idoffset=%lld)\n",id2,numrefmaps-1,k,id,idoffset);
	    fflush(stdout); exit(1);
	  }
	  refmap[k2]->Allele = k;
	  refmap[k]->Allele = k2;
	  refmap[k]->contig->contig[0].mapid = k;
	  refmap[k]->contig->contig[1].mapid = k2;
	}
      }
      delete [] id2mapid;
    }
  }

  if(VERB && TRACE){
    for(int i = 0; i < numrmaps; i++){
      if(rmap[roffset+i]->id == TraceID1)
	printf("rmap[%d]:id=%lld,mapid=%d\n",roffset+i,rmap[roffset+i]->id,rmap[roffset+i]->mapid);
    }
    for(int i = 0; i < numqmaps; i++){
      if(qmap[i]->id == TraceID2)
	printf("qmap[%d]:id=%lld,mapid=%d\n",i,qmap[i]->id,qmap[i]->mapid);
    }
    fflush(stdout);
  }

  if(VERB>=3 && giter2 == RefRepeats-1){/* save copy of maps as seen by hashtable */
    char basename[PATH_MAX];
    sprintf(basename,"%s_hash%d", prefix,giter2);
    output_bnx(basename,rmap,0,numrmaps-1,0,NULL,NULL,-1);

    sprintf(basename,"%s_phash%d", prefix,giter2);
    output_cmap(basename,qmap,0,numqmaps-1);
  }

  int numthreads = 1;
  #ifdef _OPENMP
  numthreads = omp_get_max_threads();
  if(HashQueryThreads > 0)
    numthreads = min(numthreads,HashQueryThreads);
  else
    numthreads = min(numthreads,MaxThreads);
  numthreads = max(1,min(numthreads,numqmaps));
  if(DEBUG && !(numthreads >= 1)){
    printf("max_threads()= %d, MaxThreads=%d, HashQueryThreads=%d, numqmaps=%d : numthreads=%d\n",
	   omp_get_max_threads(), MaxThreads, HashQueryThreads, numqmaps, numthreads);
    fflush(stdout);
    assert(numthreads >= 1);
  }
  #endif

  /* compute average number of sites in largest numthreads maps in qmap[] */
  int *numsite = new int[numqmaps];
  for(int i = 0; i < numqmaps; i++)
    numsite[i] = qmap[i]->numsite[hcolor];
  qsort(numsite,numqmaps,sizeof(int),(intcmp*)IntDec);

  int qsites = 0;
  for(int i = 0; i < numthreads; i++)
    qsites += numsite[i];
  if(DEBUG) assert(qsites >= 0);
  qsites /= numthreads;

  int qsitesT = qsites;
  if(numqmaps > numthreads){
    qsitesT = 0;
    for(int i = 0; i < numqmaps; i++)
      qsitesT += numsite[i];
    qsitesT /= numqmaps;
  }

  delete [] numsite;

  /* locate largest number of sites in Probe maps qmap[] */
  Nmax = 0;
  for(int i = 0; i < numqmaps; i++)
    Nmax = max(Nmax,qmap[i]->numsite[hcolor]);

  if(DEBUG) assert(HASHWIN==5 || HASHWIN==6);

  HASH_BITS = (HASHWIN==5) ? 21 : 22;
  MHASH_BITS = 14;
  HHASH_BITS = 14;
  int origHASH_BITS = HASH_BITS;
  int origMHASH_BITS = MHASH_BITS;
  int origHHASH_BITS = HHASH_BITS;

  if(1 /* WAS14 !(Hash_Bits > 0 || MHash_Bits > 0 || HHash_Bits > 0) */){

    /* scale HASH_BITS if number of inserted map sites (rsites) differs from 32K maps x 32 sites by at least a power of 2  */
    long long i = 32*1024*32;
    while(HASH_BITS < 31 && rsites > 2*i){
      HASH_BITS++;
      i <<= 1;
    }
    while(HASH_BITS > 14 && rsites*2 < i){
      HASH_BITS--;
      i >>= 1;
    }
    if(VERB /* && HASH_BITS != origHASH_BITS*/){
      printf("Changed HASH_BITS = %d -> %d due to %lld sites in inserted maps\n",origHASH_BITS,HASH_BITS,rsites);
      fflush(stdout);
    }

    /* scale MHASH_BITS if average probe map sites (qsites) differs from 16 sites by a power of 2 */
    long long sites = /* WAS USE_MIC */ pairwise ? max(qsitesT, qsites/4) : HashGC ? min(HashGC/* WAS22 HASHGC*2 */, qsites) : qsites;
    if(hashScaleDelta && NumScaleFactor > 1)
      sites = max(sites, (long long)((USE_MIC ? qsites : Nmax /* WAS100 min(qsites*2LL, (long long)Nmax)*/) * min(1.0, ScaleDelta * ScaleDeltaSize)));

    i = 16;
    while(MHASH_BITS < 31 && sites > 2*i){
      MHASH_BITS++;
      i <<= 1;
    }
    while(MHASH_BITS > 10 && sites*2 < i){
      MHASH_BITS--;
      i <<= 1;
    }

    if(VERB /* && MHASH_BITS != origMHASH_BITS*/){
      printf("Changed MHASH_BITS = %d -> %d due to %d average number of sites in largest %d probe maps (overall average sites = %d, Nmax=%d),sites=%lld,pairwise=%d,HashGC=%d,NumScale=%d\n",
	     origMHASH_BITS,MHASH_BITS,qsites,numthreads,qsitesT,Nmax,sites,pairwise,HashGC,NumScaleFactor);
      fflush(stdout);
    }

    if(!pairwise)
      MHASH_BITS = max(MHASH_BITS, HASH_BITS - 4);

    HASH_BITS = max(HASH_BITS,MHASH_BITS-1);

    HHASH_BITS = min(HHASH_BITS,MHASH_BITS);
    HHASH_BITS = max(HHASH_BITS,HASH_BITS-7);
    HHASH_BITS = max(HHASH_BITS,MHASH_BITS - 3);
  }

  if (Hash_Bits > 0 || MHash_Bits > 0 || HHash_Bits > 0){
    if(Hash_Bits > 0){
      //      MHASH_BITS += Hash_Bits-HASH_BITS;
      //      HHASH_BITS += Hash_Bits-HASH_BITS;
      HASH_BITS = Hash_Bits;
    }
    if (MHash_Bits > 0) MHASH_BITS = MHash_Bits;
    if (HHash_Bits > 0) HHASH_BITS = HHash_Bits;
  }

  if(VERB){
    printf("After final adjustment:HASH_BITS=%d,MHASH_BITS=%d,HHASH_BITS=%d->%d\n",HASH_BITS,MHASH_BITS,origHHASH_BITS,HHASH_BITS);
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  assert(HASH_BITS <= 31);
  assert(MHASH_BITS <= 31);
  assert(HHASH_BITS <= 31);

  if(HashMatrix) {
    char filename[PATH_MAX];
    sprintf(filename,"%s.matrix",prefix);
    checkFile(filename);

    if(VERB){
      printf("Generating %s (Alignment Matrix)\n",filename);
      fflush(stdout);
    }
    if((MatrixFP = fopen(filename,"w"))==NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("failed to open %s for writing:errno=%d:%s\n", filename,eno,err);
      fflush(stdout); exit(1);
    }
    printversion(MatrixFP);
    fprintf(MatrixFP,"# Loc1 and Loc2 below are left ends of window of %d sites if orientation2=0 (otherwise right ends of window).\n",HashWin);
    fprintf(MatrixFP,"# ID1 Loc1(kb) ID2 Loc2(kb) orientation2\n");
    fflush(MatrixFP);

    DisJoint = 1;
  }

  if(!pairwise){
    if(BNX_OUTPUT){ /* set bit=0 for all rmap[] entries using a bitmap : Those maps that have good hash hits will be set to bit=1 */
      if(DEBUG) assert(roffset == 0);
      MapSave = new bool[numrmaps];
      for(int i = numrmaps; --i >= 0;)
	MapSave[i] = false;
    }
  }

  if(HashBest>=0){
    if(roffset > 0){// NEW102
      if(DEBUG) assert(splitcnt > 1);
      if(DEBUG) assert(besthashscore);
      short *nbesthashscore = new short[roffset + numrmaps];
      memcpy(nbesthashscore, besthashscore, roffset * sizeof(short));
      memset(&nbesthashscore[roffset],0,numrmaps*sizeof(short));      
      delete [] besthashscore;
      besthashscore = nbesthashscore;
    } else {
      if(besthashscore) delete [] besthashscore;// NEW102
      besthashscore = new short[numrmaps];
      memset(besthashscore,0,numrmaps*sizeof(short));
    }
  }

  if(colors > 1){
    printf("WARNING: hashtable based on one color (%d) only\n", hcolor+1);
    fflush(stdout);
  }

  /* initialize YY[] & YLen[] */
  YY = new FLOAT*[numrmaps];
  Ylen = new FLOAT[numrmaps];
  for(int i = 0; i < numrmaps; i++){
    Cmap *pmap = rmap[roffset + i];
    if(DEBUG && !(Cfirst && CfirstPairs) && !(splitcnt > 1) && !(roffset + i == pmap->mapid)){
      printf("i=%d,roffset=%d,numrmaps=%d,rmap[roffset+i]->mapid=%d (DisJoint=%d,Cfirst=%d,CfirstPairs=%d,splitcnt=%d)\n",i,roffset,numrmaps,rmap[roffset+i]->mapid,DisJoint,Cfirst,CfirstPairs,splitcnt);
      fflush(stdout);
      assert(roffset + i == pmap->mapid);
    }
    YY[i] = pmap->site[hcolor];
    Ylen[i] = YY[i][pmap->numsite[hcolor]+1];
  }
  if(splitcnt > 1){
    YY -= roffset;
    Ylen -= roffset;
  }

  if(HashGC || OFFSET_FILL){/* Initialize maxOffset1[0..1] */
    OFFSET_T* mymaxOffset1 = new OFFSET_T[2*numrmaps];
    OFFSET_T* myminOffset1R = &mymaxOffset1[numrmaps];

    double OffsetKBinv = 1.0/OffsetKB;

    for(int i = 0; i < numrmaps; i++){
      int N1 = rmap[roffset+i]->numsite[hcolor];
      int HashWin1 = min(N1,HashWin);
      mymaxOffset1[i] = floor(YY[roffset + i][N1 - HashWin1] * OffsetKBinv + 0.5);
      myminOffset1R[i] = floor((Ylen[roffset + i] - YY[roffset + i][N1]) * OffsetKBinv + 0.5);
    }

    maxOffset1 = mymaxOffset1  - roffset;
    minOffset1R = myminOffset1R - roffset;
  }

  if(VERB/* HERE >=2 */){
    FLOAT maxdist = (res[hcolor] + max(HResMul * res[hcolor], HSDrange * resSD[hcolor])) * PixelLen;  // WAS40 (res[hcolor] + HSDrange * resSD[hcolor])*PixelLen;
    
    printf("res[0]= %0.8f, resSD[0]= %0.8f:  hcolor= %d, res[hcolor]= %0.6f, maxdist= %0.3f kb (bpp=%0.2f, -HSDrange %0.2f %0.2f)\n",
	   res[0],resSD[0], hcolor,res[hcolor], maxdist,PixelLen * 1000.0, HSDrange,HResMul);
    fflush(stdout);
  }

  double resKB = res[hcolor] * PixelLen * 0.5;

  OffsetKBinv = 1.0/OffsetKB;

  int MaxErr = 1 + max(NumErrI,NumErrQ) + NumResQ;
  int origDELTA = DELTA;
  DELTA = MaxErr;

  /* compute MaxFragLen, MinOffset,MaxOffset */
  int MaxRoffset = 0, MaxQoffset = 0;
  MaxFragLen = 0.0;
  TotWin = 0;

  /* allocate and initialize deltaX[][][] (for DELTA=2) */
  if(usecolor && hcolor > 0)
    colors = 2;
  deltaXinitBlock(rmap,roffset,roffset+numrmaps-1);
  if(usecolor && hcolor > 0)
    colors = 1;

  double qMaxScale = (HASH_SCALE && ScaleDelta) ? 1.0 + ScaleDelta * ScaleDeltaSize : 1.0;
  double rMaxScale = (HASH_SCALE && ScaleDelta && qmap == rmap) ? qMaxScale : 1.0;

  int TotWinIncMax = 0;
  int MaxId = 0;

  if(DEBUG) assert(HashWin >= 2);

  for(int i = 0; i < numrmaps; i++){
    Cmap *pmap = rmap[roffset + i];
    int N = pmap->numsite[hcolor];
    if(N <= HashWin)
      continue;
    if(sizeof(OFFSET_T)==2 && N > (int)MASK(15)){
      printf("hashtable cannot support maps with more than %d sites (%d sites for map id %lld)\n",MASK(15),N,pmap->id);
      fflush(stdout); exit(1);
    }
    long long origTotWin = TotWin;
    int TotWinInc = (N-HashWin) + (NumErrI>=1 ? max(0,N-HashWin-1)*HashWin : 0) + (NumErrI>=2 ? max(0,N-HashWin-2)*HashWin*(HashWin-1)/2 : 0);// NEW171
    if(DEBUG) assert(TotWinInc >= 0);
    if(TotWinInc > TotWinIncMax){
      TotWinIncMax = TotWinInc;
      MaxId = i;
    }
    TotWin += TotWinInc;
    if(DEBUG && TotWin < 0){
      Cmap *qmap = rmap[roffset + MaxId];
      unsigned int M = qmap->numsite[hcolor];
      double len = qmap->site[hcolor][M+1];
      printf("TotWin= %lld -> %lld : value overflowed :i=%d,numrmaps=%d,HashWin=%d, N=%u, NumErrI=%d,inc=%d,max=%d(imax=%d,id=%lld,M=%u,len=%0.3f kb)\n",
	     origTotWin,(long long)TotWin, i,numrmaps,HashWin,N,NumErrI,TotWinInc,TotWinIncMax,MaxId, qmap->id,M,len);
      fflush(stdout);
      assert(TotWin >= 0);
    }
    if(N > MaxSite)
      MaxSite = N;
    FLOAT *X = pmap->site[hcolor];
    for(int k = N-MaxErr; k > 0; k--)
      if((X[k+MaxErr] - X[k]) * rMaxScale > MaxFragLen){
	MaxFragLen = (X[k+MaxErr] - X[k]) * rMaxScale;
	if(VERB>=2){
	  printf("mapid=%d,k=%d,MaxErr=%d,N=%d:X[k]=%0.4f,X[k+MaxErr]=%0.4f,rMaxScale= %0.8f, MaxFragLen=%0.4f\n",
		 i,k,N,MaxErr,X[k],X[k+MaxErr],rMaxScale,MaxFragLen);
	  fflush(stdout);
	}
      }
    long long maxoffset = floor(X[N+1] * OffsetKBinv * rMaxScale + 0.5);
    if(maxoffset > MASK(31)){
      printf("hashtable cannot support maps larger than %0.3f kb : inserted map %d(id=%lld) has size %0.3f kb)\n",
	     MASK(31) * OffsetKB, i, pmap->id, X[N+1]);
      fflush(stdout); exit(1);
    }

    if(maxoffset > MaxRoffset)
      MaxRoffset = maxoffset;
  }
  if(VERB && numrmaps > 0){
    Cmap *pmap = rmap[roffset + MaxId];
    unsigned int N = pmap->numsite[hcolor];
    double len = pmap->site[hcolor][N+1];
    printf("TotWin=%lld, TotWinIncMax=%d(i=%d,id=%lld,N=%u,len= %0.3f kb), numrmaps=%d\n",(long long)TotWin,TotWinIncMax,MaxId,pmap->id,N,len,numrmaps);
    fflush(stdout);
  }

  if(qmap != rmap){
    if(usecolor && hcolor > 0)
      colors = 2;
    deltaXinitBlock(qmap,0,numqmaps-1);
    if(usecolor && hcolor > 0)
      colors = 1;

    for(int i = 0; i < numqmaps; i++){
      Cmap *pmap = qmap[i];
      int N = pmap->numsite[hcolor];
      if(N <= HashWin)
	continue;
      if(sizeof(OFFSET_T)==2 && N > (int)MASK(15)){
	printf("hashtable cannot support maps with more than %d sites (%d sites for map id %lld)\n",MASK(15),N,pmap->id);
	fflush(stdout); exit(1);
      }
      FLOAT *X = pmap->site[hcolor];
      for(int k = N-MaxErr; k > 0; k--)
	if((X[k+MaxErr] - X[k]) * qMaxScale > MaxFragLen)
	  MaxFragLen = (X[k+MaxErr] - X[k]) * qMaxScale;
      long long maxoffset = floor(X[N+1] * OffsetKBinv * qMaxScale + 0.5);
      if(maxoffset > MASK(31)){
	printf("hashtable cannot support maps larger than %0.3f kb : probe map %d(id=%lld) has size %0.3f kb)\n",
	       MASK(31) * OffsetKB, i, pmap->id, X[N+1]);
	fflush(stdout); exit(1);
      }
      if(maxoffset > MaxQoffset)
	MaxQoffset = maxoffset;
    }
  } else
    MaxQoffset = MaxRoffset;

  DELTA = origDELTA;

  /* MatchTable uses offset = offset2 - offset1 which is in the range (0-MaxRoffset)-OFFSET_SPREAD ... (MaxQoffset-0)+OFFSET_SPREAD.
     HitTable uses offset1 which is in the range 0 ...  MaxRoffset+OFFSET_SPREAD */
  MinOffset = -MaxRoffset - OFFSET_SPREAD;
  MaxOffset = MaxQoffset + OFFSET_SPREAD;
  MaxOffset = max(MaxOffset,MaxRoffset) + OFFSET_SPREAD;

  if(sizeof(OFFSET_T)==2 && (MaxOffset > (int)MASK(15) || MinOffset < -(1 << 15))){
    printf("Offset range of %d .. %d (as multiples of %0.1f Kb) cannot be respresented as short int (define OFFSET_T as \"int\" in hash.h)\n",MinOffset,MaxOffset, OffsetKB);
    printf("Nmax=%d, MinOffset=%d,MaxOffset=%d,MaxRoffset=%d,MaxQoffset=%d\n",Nmax,MinOffset,MaxOffset,MaxRoffset,MaxQoffset);
    printf("OffsetKBinv= %0.8f rMaxScale= %0.8f, qMaxScale= %0.8f\n",OffsetKBinv,rMaxScale,qMaxScale);
    fflush(stdout); exit(1);
  }

  if(sizeof(OFFSET1_T)==2 && (MaxRoffset + OFFSET_SPREAD) > (int)MASK(16)){
    printf("Offset1 range of 0 .. %u (as multiples of %0.1f Kb) cannot be respresented as unsigned short int (define OFFSET1_T as \"int\" in hash.h)\n",MaxRoffset+OFFSET_SPREAD,OffsetKB);
    fflush(stdout); exit(1);
  }
  if(RANGE >= 1 && sizeof(OFFSET1_T)==4 && ((MaxRoffset + OFFSET_SPREAD) >> ShiftOffset1) > (int)MASK(16)){
    printf("Offset1 range of 0 .. %u (as multiples of %0.1f Kb) cannot be respresented with 16 + %d (-hashrange) bits\n",MaxRoffset+OFFSET_SPREAD,OffsetKB, ShiftOffset1);
    fflush(stdout); exit(1);
  }

  /* sizing error parameters */
  if(pairwise){/* see NGentigPairScore.h : both maps have sizing error */
    fvarSF = varSF = 2.0 * SF[hcolor] * SF[hcolor];
    fvarSD = varSD = 2.0 * SD[hcolor] * fabs(SD[hcolor]);
    fvarSR = varSR = 2.0 * SR[hcolor] * SR[hcolor];
  } else {/* see GentigRefScore.h : only rmap has sizing error */
    fvarSF = varSF = SF[hcolor] * SF[hcolor];
    fvarSD = varSD = SD[hcolor] * fabs(SD[hcolor]);
    fvarSR = varSR = SR[hcolor] * SR[hcolor];
    if(VERB>=2){
      printf("c=%d:SF[c]=%0.10f,SD[c]=%0.10f,SR[c]=%0.10f:\n\tvarSF=%0.10e,varSD=%0.10e,varSR=%0.10e\n",hcolor,SF[hcolor],SD[hcolor],SR[hcolor],varSF,varSD,varSR);
      fflush(stdout);
    }
  }
  if(SIZING_BND){
    SDm = SDmax*SDmax*varSD;
    SFm = SDmax*SDmax*varSF;
    SRm = SDmax*SDmax*varSR;
    if(SRm >= 1.0){
      printf("SR[%d]=%0.6f,SDmax=%0.6f:SRm=%0.6f is >= 1.0 : reduce SR (with -MaxSR) or SDmax (3rd value of -hashgen)\n",hcolor,SR[hcolor],SDmax,SRm);
      fflush(stdout); exit(1);
    }
    if(DEBUG) assert(SRm < 1.0);
    SD1 = 0.5 * SDm;
    VAR1 = SD1*SD1 - SRm*SFm;
    SRInv = 1.0/(1.0-SRm);
  }

  /* compute size quantization boundaries */
  int sizecnt = 0;
  double sizequant[1 << SIZE_BITS];
  double sizeUB = sizequant[0] = resKB;
  while(++sizecnt < SIZES){
    double SE = sqrt(varSF + varSD*sizeUB + varSR*sizeUB*sizeUB) * SDmax;
    double origSE = SE;
    if(SE < sizeUB * RelErr)
      SE = sizeUB * RelErr;
    if(VERB>=2){
      printf("sizecnt=%d:LB=%0.6f:SE=%0.6f -> %0.6f (SDmax=%0.6f), UB=%0.6f\n",sizecnt,sizeUB,origSE,SE,SDmax,sizeUB + SE);
      fflush(stdout);
    }
    sizequant[sizecnt] = sizeUB += SE;
  }
  if(VERB){
    printf("Largest representable interval with %d bits = %0.3f kb, largest actual interval=%0.3f kb, qMaxScale= %0.8f\n", SIZE_BITS, sizequant[SIZES-2], MaxFragLen, qMaxScale);
    if(VERB && (VERB>=2 || TRACE))
      for(int i = 0; i < SIZES; i++)
	printf("size interval %d: %0.3f .. %0.3f\n", i, (i ? sizequant[i-1] : 0.0), sizequant[i]);
    fflush(stdout);
  }
  MaxIntLen = (int)floor(MaxFragLen * QLenInv + 0.5);
  IntLen = new unsigned int[MaxIntLen+1];

  sizecnt = 1;/* forces all IntLen[] values to between 1 and SIZES-1 : eg values below resKB will be treated same as resKB */
  for(int i = 0; i <= MaxIntLen; i++){
    double size = i * QLen;
    while(sizecnt < SIZES-1 && size >= sizequant[sizecnt])
      sizecnt++;
    IntLen[i] = sizecnt;
  }
  if(VERB>=2){
    for(int i = 0; i <= MaxIntLen; i++)
      printf("IntLen[i=%d]=%u : i*QLen=%0.3f\n",i,IntLen[i],i*QLen);
    fflush(stdout);
  }

  /* compute SizeDelta[0..SizeLen-1] : hash key offsets to account for sizing error of +- 1 in every size in the window */
  SizeLen = 3;
  for(int i = 1; i < HashWin; i++)
    SizeLen *= 3;
  assert(SizeLen < (int)MASK(31) && SizeLen > 0);
  SizeDelta = new long long[SizeLen];

  int delta[3] = {0,1,-1};/* quantized size delta values to consider per interval */
  int deltawin[MAXR_WIN];/* deltawin[0..HashWin-1] is the index into delta[] for each interval in the window */
  for(int i=0; i < HashWin; i++)
    deltawin[i] = 0;
  for(unsigned int i = 0; i < SizeLen; i++){
    if(TRACE2){
      printf("SizeDelta[%d/%d]=%d(%d)",i,SizeLen,deltawin[0],delta[deltawin[0]]);
      for(int w = 1; w < HashWin; w++)
	printf(",%d(%d)", deltawin[w],delta[deltawin[w]]);
      printf("\n");
      fflush(stdout);
    }
    SizeDelta[i] = hashkey(deltawin,delta);

    /* increment deltawin[0..HashWin-1] as if these were Base-3 digits of a number */
    for(int k = 0; k < HashWin; k++){
      if(++deltawin[k] < 3)
	break;
      deltawin[k] = 0;
    }
  }

  if(VERB>=2){
    printf("hash_init : constants & tables initialized:CPU time=%0.6f, wall time=%0.6f secs\n",mtime(),wtime());
    dumpmemmap();
    fflush(stdout);
  }

  return;

}


/** look up all matches for one map (single orientation):
   generate all possible windows and call keyquery(key,mapid2,flip2,X,site2,err,offset2)
 */
/* If probe map is flipped (flip2==1), site2 is reduced each time win2 is expended so the right end of the window stays the same :
   This way in unflipped orientation of probe map the left end of window stays the same, improving likelihood of catching duplicate Hits with the same Roffset1.
   Otherwise each perfect match case will produce an additional match for each allowed missing OR misresolved site per window : these additional
   matches are really duplicates but will not be caught by HitTable since Roffset1 is not the same */


template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::hash_queryF(Cmap *pmap, int mapid2, int flip2, FLOAT *Xrev, FLOAT *Xscale, int *Xres, int N, int tid)
{

  if(VERB>=2){
    int tid2 = 0;// Level 2 thread id
    #ifdef _OPENMP
    tid2 = omp_get_thread_num();
    #endif
    
    #pragma omp critical
    {
      printf("tid=%d,%d:mapid2=%d(id=%lld),flip2=%d,N=%d: start of hash_queryF: wtime= %0.6f\n",tid,tid2,mapid2,qmap[mapid2]->id,flip2,N,wtime());
      if(!MULTIMATCH_TWOTHREADS)
	printf("\t HashTable:probes=%lld,%lld,hits=%lld,inline col=%lld(%0.3f/hit),col=%lld(%0.3f/hit),matches=%lld(%0.3f/hit),total entries=%lld(%0.3f/match)\n",
	       HTindcnt, HTqindcnt, HTfindcnt,HTqcolcnt-HTmatchcnt,((double)(HTqcolcnt-HTmatchcnt))/max(1,HTfindcnt),
	       HTcolcnt,(double)(HTcolcnt)/max(1,HTfindcnt), HTmatchcnt, ((double)HTmatchcnt)/max(1,HTfindcnt), HTentrycnt, ((double)HTentrycnt)/max(1,HTmatchcnt));
      fflush(stdout);
    }
  }

  /* At this point Match Table is reset/empty */
  //  PFLOAT *win2 = (PFLOAT*) alloca(MAXR_WIN*sizeof(PFLOAT));
  PFLOAT win2[MAXR_WIN];

  FLOAT maxdist = (res[hcolor] + max(HResMul * res[hcolor], HSDrange * resSD[hcolor])) * PixelLen;  // WAS40 (res[hcolor] + HSDrange * resSD[hcolor])*PixelLen;

  //  double origres = res[hcolor], origresSD = resSD[hcolor];
  if(DEBUG) assert(ScaleFactor[0] == 1.0);

  if(DEBUG){
   assert(minSite2 == 0);
   assert(MatchFirst == -1);
   assert(MatchOverflowNext == MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) + 1);
   assert(MatchOverflowFree == -1);
  }

  FLOAT len = Xrev[N+1];// NOTE : Unless OFFSET_SCALE_FIX, len, offset2 are based on unscaled X. Xres is always based on unscaled X
  if(NumResQ) 
    for(int site2 = 1; site2 < N; site2++)
      Xres[site2] = (Xrev[site2+1] - Xrev[site2] < maxdist) ? 1 : 0;

  OFFSET_T offset2;

  if(DEBUG>=2)
    for(int J = 0; J <= N+1; J++)
      Xscale[J] = aNaN;

  for(int site2 = 1; site2 <= N-HashWin; site2++){
    if(VERB && (VERB>=2 || (!(site2 % 100000) /* && numrmaps >= 50000*/ ))){
      #pragma omp critical
      {
        printf("Accumulating Hits for mapid2=%d,flip2=%d,site2=%d/%d:tid=%d,MatchOverflow=%d/%d(HWM=%d+%lld),BestHWM=%d,OutHWM=%lu,HitMemHWM=%d,MatchCnt=%lld:CPU time=%0.6f, wall time=%0.6f secs\n",
	  mapid2,flip2,site2,N,tid,MatchOverflowNext,MatchOverflowMax,MatchMemHWM,MatchBufHWM,BestHWM,OutHWM,HitMemHWM,MatchCnt,mtime(),wtime());
	fflush(stdout);
      }
    }

    for(int scaleID = (HASH_SCALE && hashScaleDelta) ? NumScaleFactor : 1; --scaleID >= 0; HitToMatch(mapid2,flip2,(scaleID ? Xscale : Xrev),N,site2,offset2,scaleID,tid) ){
      FLOAT scale = ScaleFactor[scaleID];
      FLOAT *X = scaleID ? Xscale : Xrev;
      if(scaleID > 0){
	int Jmin = flip2 ? max(0, site2 - NumErrQ - NumResQ) : site2;
	int Jmax = flip2 ? min(site2 + HashWin, N+1) : min(site2 + HashWin + NumErrQ + NumResQ, N+1);

	// WAS21	for(int J=0;J <= N+1;J++)
	for(int J= Jmin;J <= Jmax;J++)
	  Xscale[J] = Xrev[J] * scale;

	if(DEBUG>=2){
	  if(Jmin > 0)
	    Xscale[Jmin - 1] = aNaN;
	  if(Jmax <= N)
	    Xscale[Jmax+1] = aNaN;
	}
      }
      FLOAT Xlen = OFFSET_SCALE_FIX ? Xrev[N+1] * scale : Xrev[N+1];

      //      res[hcolor] = origres * scale;
      //      resSD[hcolor] = origresSD * scale;

      /* use window without site errors */

      window<HASH_SCALE>(win2, X, pmap, mapid2, flip2, site2, 0);
      long long key = hashkey(win2);
      if(OFFSET_SCALE_FIX)
	offset2 = floor((flip2 ? Xlen - X[site2+HashWin] : X[site2])*OffsetKBinv + 0.5);// based on rescaled X[]
      else
	offset2 = floor((flip2 ? len - Xrev[site2+HashWin] : Xrev[site2])*OffsetKBinv + 0.5); // based on original Xrev[]
      
      if(DEBUG>=2) assert(0 <= offset2 && offset2 <= MaxOffset);

      if(VERB>=3 && flip2 && site2 == 8 && scaleID == 14){
	printf("hash_queryF: mapid2=%d,flip2=%d,site2=%d,scaleID=%d(scale= %0.3f): offset2= %d (Xlen= %0.3f->%0.3f, X[site2]=%0.3f,X[site2+HashWin]=%0.3f,OffsetKBinv=%0.3f)\n",
	       mapid2,flip2,site2,scaleID,scale,offset2,len,Xlen,X[site2],X[site2+HashWin],OffsetKBinv);
	fflush(stdout);
      }

      keyquery(key,win2,mapid2,flip2,X,site2,0,offset2,scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
    
      /* PROBFLIP_FIX : If flip2==1, site2 should be reduced each time win2 is expanded so the right end of the window stays the same : this way in unflipped orientation of probe map the left end
	 of window stays same improving likelihood of catching duplicates Hits with same site1. Otherwise each perfect match case will produce an additional match for each allowed misresolved 
	 site per window : these matches are really duplicates but will not be caught by HitTable since site1 is not the same */

      int nsite2 = site2;

      if(NumErrQ + NumResQ <= 0 || (flip2==1 ? (--nsite2 <= 0) : (nsite2 >= N-HashWin)))
	continue;

      if(NumResQ){ /* consider single misresolved site per window */
	for(int mis = 0; mis <= HashWin; mis++){
	  if(Xres[nsite2+mis]){
	    windowM1(win2, X, pmap, mapid2, flip2, nsite2, mis);
	    key = hashkey(win2);
	    keyquery(key, win2, mapid2, flip2, X, nsite2, (mis+1) << 8, offset2,scaleID);
	  }
	}
      }

      if(NumErrQ){   /* consider single missing site in window */
	for(int err = 1; err <= HashWin; err++){
	  window<HASH_SCALE>(win2, X, pmap, mapid2, flip2, nsite2, err);
	  key = hashkey(win2);
	  keyquery(key, win2, mapid2, flip2, X, nsite2, err, offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
	}
      }
    
      if(NumErrQ+NumResQ <= 1 || (flip2==1 ? (--nsite2 <= 0) : (nsite2 >= N-HashWin-1)))
	continue;
    
      if(NumResQ >= 2){/* consider 2 misresolved sites per window */
	for(int mis1 = 0; mis1 <= HashWin; mis1++){
	  if(Xres[nsite2 + mis1]){
	    for(int mis2 = mis1+1; mis2 <= HashWin+1; mis2++){
	      if(Xres[nsite2 + mis2]){
		windowM2(win2, X, pmap, mapid2, flip2, nsite2, mis1, mis2);
		key = hashkey(win2);
		keyquery(key, win2, mapid2, flip2, X, nsite2, (((mis2+1) << 4) | (mis1+1)) << 8, offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
	      }
	    }
	  }
	}
      }

      if(NumErrQ && NumResQ){/* consider 1 missing sites + 1 misresolved site per window */
	for(int mis = 0; mis <= HashWin+1; mis++){
	  if(Xres[nsite2 + mis]){
	    for(int err = 1; err <= HashWin; err++){
	      windowM1E1(win2, X, pmap, mapid2, flip2, nsite2, mis, err);
	      key = hashkey(win2);
	      keyquery(key, win2, mapid2, flip2, X, nsite2, ((mis+1) << 8) | err, offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
	    }
	  }
	}
      }

      if(NumErrQ >= 2){  /* consider 2 missing sites per window */
	for(int err1 = 1; err1 <= HashWin; err1++){
	  for(int err2 = err1 + 1; err2 <= HashWin + 1; err2++){
	    window2<HASH_SCALE>(win2, X, pmap, mapid2, flip2, nsite2, err1, err2);
	    key = hashkey(win2);
	    keyquery(key, win2, mapid2, flip2, X, nsite2, err1 | (err2 << 4), offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
	  }
	}
      }
    
      if(NumErrQ+NumResQ <= 2 || (flip2==1 ? (--nsite2 <= 0) : (nsite2 >= N-HashWin-2)))
	continue;

      if(NumResQ >= 3){/* consider 3 misresolved sites per window */
	for(int mis1 = 0; mis1 <= HashWin; mis1++){
	  if(Xres[nsite2 + mis1]){
	    for(int mis2 = mis1+1; mis2 <= HashWin + 1; mis2++){
	      if(Xres[nsite2 + mis2]){
		for(int mis3 = mis2+1; mis3 <= HashWin + 2; mis3++){
		  if(Xres[nsite2 + mis3]){
		    windowM3(win2, X, pmap, mapid2, flip2, nsite2, mis1, mis2, mis3);
		    key = hashkey(win2);
		    keyquery(key, win2, mapid2, flip2, X, nsite2, (((mis3+1) << 8) | ((mis2+1) << 4) | (mis1+1)) << 8, offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
		  }
		}
	      }
	    }
	  }
	}
      }

      if(NumErrQ && NumResQ >= 2){/* consider single missing site + 2 misresolved sites per window */
	for(int mis1 = 0; mis1 <= HashWin + 1; mis1++){
	  if(Xres[nsite2 + mis1]){
	    for(int mis2 = mis1 + 1; mis2 <= HashWin + 2; mis2++){
	      if(Xres[nsite2 + mis2]){
		for(int err = 1; err <= HashWin; err++){
		  windowM2E1(win2, X, pmap, mapid2, flip2, nsite2, mis1, mis2, err);
		  key = hashkey(win2);
		  keyquery(key, win2, mapid2, flip2, X, nsite2, ((((mis2+1) << 4)|(mis1+1)) << 8) | err, offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
		}
	      }
	    }
	  }
	}
      }

      if(NumErrQ >= 2 && NumResQ){/* consider 2 missing sites + 1 misresolved site per window */
	for(int mis = 0; mis <= HashWin + 2; mis++){
	  if(Xres[nsite2 + mis]){
	    for(int err1 = 1; err1 <= HashWin; err1++){
	      for(int err2 = err1 + 1; err2 <= HashWin + 1; err2++){
		windowM1E2(win2, X, pmap, mapid2, flip2, nsite2, mis, err1, err2);
		key = hashkey(win2);
		keyquery(key, win2, mapid2, flip2, X, nsite2, ((mis+1) << 8) | err1 | (err2 << 4), offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
	      }
	    }
	  }
	}
      }

      if(NumErrQ >= 3){   /* Handle 3 missing sites per window */
	printf("hash_queryF() not implemented for NumErrP=%d, NumResP=%d\n", NumErrQ, NumResQ);    
	fflush(stdout); exit(1);
      }

      if(NumErrQ+NumResQ <= 3 || (flip2==1 ? (--nsite2 <= 0) : (nsite2 >= N-HashWin-3)))
	continue;

      /* handle any combination of 4 total miss + misresolved sites per window */

      if(NumResQ >= 4){/* consider 4 misresolved sites per window */
	for(int mis1 = 0; mis1 <= HashWin; mis1++){
	  if(Xres[nsite2 + mis1]){
	    for(int mis2 = mis1+1; mis2 <= HashWin + 1; mis2++){
	      if(Xres[nsite2 + mis2]){
		for(int mis3 = mis2+1; mis3 <= HashWin + 2; mis3++){
		  if(Xres[nsite2 + mis3]){
		    for(int mis4 = mis3+1; mis4 <= HashWin + 3; mis4++){
		      if(Xres[nsite2 + mis4]){
			windowM4(win2,X,pmap,mapid2,flip2,nsite2,mis1,mis2,mis3,mis4);
			key = hashkey(win2);
			keyquery(key, win2, mapid2, flip2, X, nsite2, (((int) (((mis4+1)<<12)| ((mis3+1) << 8) | ((mis2+1) << 4) | (mis1+1))) << 8), 
	                         offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      if(NumErrQ && NumResQ >= 3){/* consider 1 missing site + 3 misresolved sites per window */
	for(int mis1 = 0; mis1 <= HashWin + 1; mis1++){
	  if(Xres[nsite2 + mis1]){
	    for(int mis2 = mis1+1; mis2 <= HashWin + 2; mis2++){
	      if(Xres[nsite2 + mis2]){
		for(int mis3 = mis2+1; mis3 <= HashWin + 3; mis3++){
		  if(Xres[nsite2 + mis3]){
		    for(int err = 1; err <= HashWin; err++){
		      windowM3E1(win2, X, pmap, mapid2, flip2, nsite2, mis1, mis2, mis3, err);
		      key = hashkey(win2);
		      keyquery(key, win2, mapid2, flip2, X, nsite2, ((((mis3+1) << 8) | ((mis2+1) << 4) | (mis1+1)) << 8) | err, 
	                    offset2,scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    
      if(NumErrQ >= 2 && NumResQ >= 2){/* consider 2 missing sites + 2 misresolved sites per window */
	for(int mis1 = 0; mis1 <= HashWin + 2; mis1++){
	  if(Xres[nsite2 + mis1]){
	    for(int mis2 = mis1 + 1; mis2 <= HashWin + 3; mis2++){
	      if(Xres[nsite2 + mis2]){
		for(int err1 = 1; err1 <= HashWin; err1++){
		  for(int err2 = err1 + 1; err2 <= HashWin + 1; err2++){
		    windowM2E2(win2, X, pmap, mapid2, flip2, nsite2, mis1, mis2, err1, err2);
		    key = hashkey(win2);
		    keyquery(key, win2, mapid2, flip2, X, nsite2, ((((mis2+1) << 4)|(mis1+1)) << 8) | err1 | (err2 << 4), 
	                     offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
		  }
		}
	      }
	    }
	  }
	}
      }

      if(NumErrQ+NumResQ <= 4 || (flip2==1 ? (--nsite2 <= 0) : (nsite2 >= N-HashWin-4)))
	continue;

      /* handle any combination of 5 total missing + misresolved sites per window */

      if(NumResQ >= 5){/* consider 5 misresolved sites per window */
	for(int mis1 = 0; mis1 <= HashWin; mis1++){
	  if(Xres[nsite2 + mis1]){
	    for(int mis2 = mis1+1; mis2 <= HashWin + 1; mis2++){
	      if(Xres[nsite2 + mis2]){
		for(int mis3 = mis2+1; mis3 <= HashWin + 2; mis3++){
		  if(Xres[nsite2 + mis3]){
		    for(int mis4 = mis3+1; mis4 <= HashWin + 3; mis4++){
		      if(Xres[nsite2 + mis4]){
			for(int mis5 = mis4+1; mis5 <= HashWin + 4; mis5++){
			  if(Xres[nsite2 + mis5]){
			    windowM5(win2,X,pmap,mapid2,flip2,nsite2,mis1,mis2,mis3,mis4,mis5);
			    key = hashkey(win2);
			    keyquery(key, win2, mapid2, flip2, X, nsite2, (((int) (((mis5+1)<<16)|((mis4+1)<<12)| ((mis3+1) << 8) | ((mis2+1) << 4) | (mis1+1))) << 8),
	                             offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      if(NumErrQ && NumResQ >= 4){/* consider 1 missing site + 4 misresolved sites per window */
	for(int mis1 = 0; mis1 <= HashWin + 1; mis1++){
	  if(Xres[nsite2 + mis1]){
	    for(int mis2 = mis1+1; mis2 <= HashWin + 2; mis2++){
	      if(Xres[nsite2 + mis2]){
		for(int mis3 = mis2+1; mis3 <= HashWin + 3; mis3++){
		  if(Xres[nsite2 + mis3]){
		    for(int mis4 = mis3+1; mis4 <= HashWin + 4; mis4++){
		      if(Xres[nsite2 + mis4]){
			for(int err = 1; err <= HashWin; err++){
			  windowM4E1(win2,X,pmap,mapid2,flip2,nsite2,mis1,mis2,mis3,mis4,err);
			  key = hashkey(win2);
			  keyquery(key, win2, mapid2, flip2, X, nsite2, (((int) (((mis4+1)<<12)| ((mis3+1) << 8) | ((mis2+1) << 4) | (mis1+1))) << 16) | err,
	                           offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      if(NumErrQ >= 2 && NumResQ >= 3){/* consider 2 missing sites + 3 misresolved sites per window */
	for(int mis1 = 0; mis1 <= HashWin + 2; mis1++){
	  if(Xres[nsite2 + mis1]){
	    for(int mis2 = mis1+1; mis2 <= HashWin + 3; mis2++){
	      if(Xres[nsite2 + mis2]){
		for(int mis3 = mis2+1; mis3 <= HashWin + 4; mis3++){
		  if(Xres[nsite2 + mis3]){
		    for(int err1 = 1; err1 <= HashWin; err1++){
		      for(int err2 = err1 + 1; err2 <= HashWin + 1; err2++){
			windowM3E2(win2, X, pmap, mapid2, flip2, nsite2, mis1, mis2, mis3, err1,err2);
			key = hashkey(win2);
			keyquery(key, win2, mapid2, flip2, X, nsite2, ((((mis3+1) << 8) | ((mis2+1) << 4) | (mis1+1)) << 16) | err1 | (err2 << 4), 
	                         offset2, scaleID);/* appends hash matches to Hit Table while removing duplicate (mapid1,offset1) */
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      //    printf("hash_queryF() not implemented for NumErrP=%d, NumResP=%d\n", NumErrQ, NumResQ);    
      //    fflush(stdout); exit(1);
    } // scaleID = NumScaleFactor-1 .. 0 

  }// site2 = 1 .. N-HashWin

  if(VERB>=2 && !tid){
    printf("tid=%d:Calling MatchFind(%d,%d)\n",tid,mapid2,flip2);
    fflush(stdout);
  }

  /* save best offset in match table for each mapid1 (if score >= HashThresh) */
  MatchFind(mapid2, flip2, 0, tid);

  if(VERB>=2 && !tid){
    printf("tid=%d:Finished MatchFind(%d,%d)\n",tid,mapid2,flip2);
    fflush(stdout);
  }

  /* if(NumScaleFactor > 1){
    res[hcolor] = origres;
    resSD[hcolor] = origresSD;
    } */
}

/** top level hash query function : look up all matches for one map (both orientations) */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::hash_query(Cmap *pmap, int tid, int *Xflag, FLOAT *Xrev, FLOAT *Xscale, CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *matchtable2)
{
  size_t origmatchcnt = matchcnt;/* start of buffer region */

  int M = pmap->numsite[hcolor];
  if(DEBUG) assert(M <= Nmax);

  if(M <= HashWin)
    return;

  if(HashGC && hashScaleDelta && NumScaleFactor > 1){
    MTHashGC = max(HashGC, (int)(M * min(1.0,ScaleDelta * ScaleDeltaSize)));
    if(matchtable2)
      matchtable2->MTHashGC = MTHashGC;
  }

  int mapid2 = pmap->mapid;
  if(DEBUG) assert(mapid2 < numqmaps);
  FLOAT *X = pmap->site[hcolor];

  for(int J=0; J <= M+1; J++)
    Xrev[J] = X[M+1] - X[M+1-J];

  if(!matchtable2){
    hash_queryF(pmap, mapid2, 0, X, Xscale, Xflag, M, tid);
    if(VERB>=2){
      #pragma omp critical
      {
	printf("tid=%d:mapid2=%d(id=%lld),flip=0,M=%d: matchcnt=%lu: wtime= %0.6f\n",tid,mapid2,qmap[mapid2]->id,M,matchcnt-origmatchcnt,wtime());
	if(mapid2 == MAPID2 && 0==FLIP2)
	  for(size_t i = origmatchcnt; i < matchcnt; i++)
	    printf("match[%lu/%lu]:id1=%d(id=%lld),id2=%d(id=%lld),or=%d,offset=%d,hashscore=%u\n",i-origmatchcnt, matchcnt-origmatchcnt,
		   matches[i].id1,qmap[matches[i].id1]->id,matches[i].id2,rmap[matches[i].id2]->id,matches[i].orientation,matches[i].offset,matches[i].hashscore);
	fflush(stdout);
      }
    }
    size_t origmatchcnt2 = matchcnt;

    hash_queryF(pmap, mapid2, 1, Xrev, Xscale, Xflag, M, tid);

    if(VERB>=2){
      #pragma omp critical
      {
	printf("tid=%d:mapid2=%d(id=%lld),flip=1,M=%d: matchcnt=%lu: wtime= %0.6f\n",tid,mapid2,qmap[mapid2]->id,M,matchcnt-origmatchcnt2,wtime());
	if(mapid2 == MAPID2 && 1==FLIP2)
	  for(size_t i = origmatchcnt2; i < matchcnt; i++)
	    printf("match[%lu/%lu]:id1=%d(id=%lld),id2=%d(id=%lld),or=%d,offset=%d,hashscore=%u\n",i-origmatchcnt2, matchcnt-origmatchcnt2,
		   matches[i].id1,qmap[matches[i].id1]->id,matches[i].id2,rmap[matches[i].id2]->id,matches[i].orientation,matches[i].offset,matches[i].hashscore);
	fflush(stdout);
      }
    }
  } else {
    if(DEBUG) assert(MULTIMATCH_TWOTHREADS);
    CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *Matchtable[2] = {this, matchtable2};

    // Avoid 2 threads using the same cache line
    int Nmax1 = ((Nmax + 15)/16) * 16; // NEW164
    int Nmax2 = ((Nmax + 2 + 15)/16) * 16;// NEW164

    #pragma omp parallel num_threads(2)
    {
      #pragma omp for schedule(static,1)
      for(int flip2 = 0; flip2 <= 1; flip2++){
	Matchtable[flip2] -> hash_queryF(pmap, mapid2, flip2, flip2 ? Xrev : X, &Xscale[flip2 ? Nmax2 : 0], &Xflag[flip2 ? Nmax1 : 0], M, tid);

	if(mapid2 == MAPID2 && flip2 == FLIP2){
	  #pragma omp critical
	  {
	    printf("tid=%d:mapid2=%d(id=%lld),flip2=%d,M=%u: matchcnt=%lu: wtime= %0.6f\n",tid,mapid2,qmap[mapid2]->id,flip2,M,Matchtable[flip2]->matchcnt,wtime());
	    CHashMatch *matchesF = Matchtable[flip2]->matches;
	    for(size_t i = 0; i < Matchtable[flip2]->matchcnt; i++)
	      printf("match[%lu/%lu]:id1=%d(id=%lld),id2=%d(id=%lld),or=%d,offset=%d,hashscore=%d\n",i, Matchtable[flip2]->matchcnt,
		     matchesF[i].id1,qmap[matchesF[i].id1]->id,matchesF[i].id2,rmap[matchesF[i].id2]->id,matchesF[i].orientation,matchesF[i].offset,matchesF[i].hashscore);
	    fflush(stdout);
	  }
	}
      }
    }

    if(VERB>=2){
      #pragma omp critical
      {
	printf("tid=%d:mapid2=%d,M=%d: matchcnt=%lu, matchtable2->matchcnt=%lu: wtime= %0.6f\n",tid,mapid2,M,matchcnt,matchtable2->matchcnt,wtime());
	fflush(stdout);
      }
    }

    /* append output buffer from matchtable2 to current MatchTable */
    if(matchcnt + matchtable2->matchcnt > matchmax){
      maxmatchalloc(matchcnt, matchcnt + matchtable2->matchcnt, matchmax, matches);
      if(STAT && matchmax > OutputHWM)
	OutputHWM = matchmax;
    }
    memcpy(&matches[matchcnt], matchtable2->matches, sizeof(CHashMatch) * matchtable2->matchcnt);
    matchcnt += matchtable2->matchcnt;
    matchtable2->matchcnt = 0;

    /* reduce memory used in output buffer for matchtable2 */
    if(matchtable2->matchmax > 4*1024){
      delete [] matchtable2->matches;
      matchtable2->matchmax = 4*1024;
      matchtable2->matches = new CHashMatch[matchtable2->matchmax];    
    }
  }

  if(Gpairwise && HashMaxHits > 0 && matchcnt > (size_t)HashMaxHits){/* limit matches to highest scoring HashMaxHits results */
    if(DEBUG) assert(origmatchcnt == 0);

    /* sort by (-hashscore) */
    qsort(&matches[origmatchcnt],matchcnt-origmatchcnt,sizeof(CHashMatch),(intcmp*)CHashMatchDecHscore);

    size_t i = origmatchcnt + HashMaxHits;
    int maxscore = matches[i - 1].hashscore;
    for(;i < matchcnt; i++)
      if(matches[i].hashscore < maxscore)
	break;
    if(VERB>=2){
      #pragma omp critical
      {
	printf("tid=%d:mapid2=%d(id=%lld),M=%u: HashMaxHits= %d, maxscore= %d, matchcnt= %lu -> %lu\n",tid,mapid2,qmap[mapid2]->id,M,HashMaxHits,maxscore,matchcnt,i);
	fflush(stdout);
      }
    }
    matchcnt = i;
  }

  /* now sort the output buffer region matches[origmatchcnt .. matchcnt-1] by id2, orientation,(-hashscore),offset (All id1 == qmap[mapid2]->id) */
  qsort(&matches[origmatchcnt],matchcnt-origmatchcnt,sizeof(CHashMatch),(intcmp*)CHashMatchIncId2);

  totalMcnt += matchcnt-origmatchcnt;/* thread-local counter */

  if(!((++totalQcnt) % 16)){/* avoid updating shared global counter too often */
    #pragma omp critical
    {// NOTE : cannot use atomic update, since it is also updated in another critical section
      TotalMatchCnt += totalMcnt;
    }

    totalMcnt = 0;/* reset thread-local counter */
  }

  /* count number of distinct maps in output buffer region */
  size_t MapCnt = 1;
  for(size_t i = origmatchcnt+1; i < matchcnt; i++){
    CHashMatch *r = &matches[i];
    if(r[-1].id2 != r->id2)
      MapCnt++;
  }

  if(VERB>=2){
    #pragma omp critical
    {
      printf("mapid2=%d:MapCnt=%lu,matchcnt=%lu\n",mapid2,MapCnt,matchcnt-origmatchcnt);
      fflush(stdout);
    }
  }

  if(TSAN || splitcnt > 1 || matchcnt-origmatchcnt > maxMatchesPerRefid || MapCnt > maxMapsPerRefid){
    #pragma omp critical
    {
      if(splitcnt > 1) {
	if(DEBUG) assert(0 <= mapid2 && mapid2 < Gnumqmaps);

	if(VERB>=2){
	  printf("mapid2=%d:MapCnt=%lu(cur=%lu,max=%lu),matchcnt=%lu,origmatchcnt=%lu:MatchesPerRefid[mapid2]=%lu -> %lu, maxMatchesPerRefid=%lu -> %lu\n",
		 mapid2,MapCnt,MapsPerRefid[mapid2],maxMapsPerRefid,matchcnt,origmatchcnt,MatchesPerRefid[mapid2],MatchesPerRefid[mapid2] + matchcnt - origmatchcnt,
		 maxMatchesPerRefid,max(maxMatchesPerRefid,MatchesPerRefid[mapid2]+matchcnt-origmatchcnt));
	  fflush(stdout);
	}

	if((MatchesPerRefid[mapid2] += matchcnt - origmatchcnt) > maxMatchesPerRefid)
	  maxMatchesPerRefid = MatchesPerRefid[mapid2];

	if((MapsPerRefid[mapid2] += MapCnt) > maxMapsPerRefid)
	  maxMapsPerRefid = MapsPerRefid[mapid2];
	
      } else {
	if(matchcnt - origmatchcnt > maxMatchesPerRefid){
	  if(VERB>=2){
	    printf("mapid2=%d:MapCnt=%lu(max=%lu),matchcnt=%lu,origmatchcnt=%lu:maxMatchesPerRefid=%lu -> %lu\n",
		   mapid2,MapCnt,maxMapsPerRefid,matchcnt,origmatchcnt,maxMatchesPerRefid,matchcnt-origmatchcnt);
	    fflush(stdout);
	  }
	  maxMatchesPerRefid = matchcnt - origmatchcnt;
	}
	if(MapCnt > maxMapsPerRefid)
	  maxMapsPerRefid = MapCnt;
      }
    }
  }

  OutputFlush();/* flush output buffer */

  if(VERB && (VERB>=2 || M > 100000)){
    #pragma omp critical
    {
      printf("Completed hash_query for mapid=%d(%d sites):MatchOverflowMax=%d,MatchBufMax=%lld:CPU time=%0.6f, wall time=%0.6f secs\n",mapid2,M,MatchOverflowMax,MatchBufMax,mtime(),wtime());
      fflush(stdout);
    }
  }

  // memory was freed at end of MatchFind() 
}

/* flush output buffer matches[0..matchcnt-1] (if needed) */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>::OutputFlush()
{
  if(TRACE){
    for(unsigned int i = 0; i < matchcnt; i++){
      CHashMatch *r = &matches[i];
      if(qmap[r->id1]->id == TraceID2 && rmap[r->id2]->id == TraceID1){
	if(RANGE>=1)
	  printf("Flushing matches[%u/%lu]:mapid2=id1=%d(id=%lld),mapid1=id2=%d(id=%lld),or=%d,sc=%d,hashscore=%d,offset=%d,Loc2=%0.1f..%0.1f kb, tid=%d\n",
		 i,matchcnt,r->id1,qmap[r->id1]->id,r->id2,rmap[r->id2]->id,r->orientation,r->scaleID,r->hashscore,r->offset,r->MinLoc2,r->MaxLoc2,tid);
	else
	  printf("Flushing matches[%u/%lu]:mapid2=id1=%d(id=%lld),mapid1=id2=%d(id=%lld),or=%d,sc=%d,hashscore=%d,offset=%d, tid=%d\n",
		 i,matchcnt,r->id1,qmap[r->id1]->id,r->id2,rmap[r->id2]->id,r->orientation,r->scaleID,r->hashscore,r->offset,tid);
      }
    }
  }

  if(BNX_OUTPUT && !Gpairwise){/* mark bit corresponding to (rmap[p->id2]) */
    #pragma omp critical(MapSave)
    {/* MapSave[] is a shared array : Typically only one thread will be flushing output at a time */
      for(unsigned int i = 0; i < matchcnt; i++){
	CHashMatch *p = &matches[i];
	int index = p->id2;
	MapSave[index] = true;
      }
    }
  }

  if(VERB>=3 && giter2 == RefRepeats-1 && matchcnt > 0){
    char filename[PATH_MAX];
    sprintf(filename, "%s_I%d_matches%d.txt",output_prefix, giter2, matches[0].id1);

    #pragma omp critical
    {
      printf("Writing matches for probe id=%d to %s\n", matches[0].id1, filename);
      fflush(stdout);
    }

    FILE *fp = fopen(filename,"w");
    if(fp == NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("Failed to open file %s for write:errno=%d:%s\n",filename,eno,err);
      fflush(stdout); exit(1);
    }
    fprintf(fp,"# id1 id2 offset hashscore orientation\n");
    for(size_t i = 0; i < matchcnt; i++)
      fprintf(fp,"%d %d %d %d %d\n",matches[i].id1,matches[i].id2,matches[i].offset, matches[i].hashscore,matches[i].orientation);
    fclose(fp);
  }

  if(matchcnt > 0){
    if(STREAM_DEBUG){
      printf("Appending %lu matches to stream buffer %s\n",matchcnt, filename);
      if( STREAM_DEBUG >= 2 ){
        for(unsigned int i = 0; i < matchcnt; i++){
          CHashMatch *r = &matches[i];
	  printf("Flushing matches[%u/%lu]:id1=%lld,id2=%lld,or=%d,hashscore=%d,offset=%d,tid=%d\n",
		 i,matchcnt,(long long)r->id1,(long long)r->id2,r->orientation,r->hashscore,r->offset,tid);
        }
      }
    }

    size_t n = fwrite_LOOP(matches,sizeof(CHashMatch),matchcnt,fp);
    if(n != matchcnt){
      #pragma omp critical
      {
        printf("fwrite=%lu: write of %lu Matches (%lu bytes each) to %s failed\n", 
          n, matchcnt, sizeof(CHashMatch), filename);
	fflush(stdout);exit(1);
      }
    }
    matchcnt = 0;
  }

  /* reduce memory used in output buffer */
  if(matchmax > 4*1024){
    delete [] matches;
    matchmax = 4*1024;
    matches = new CHashMatch[matchmax];    
  }
}

/* comparison function to sort CHashMatch array in increasing order of id1,id2,orientation,(-hashscore),offset.
   Special case : NULL pointers (EOF) are treated as having the largest possible values */
static inline int CHashMatchInc(CHashMatch *p1, CHashMatch *p2)
{
  if(!p1)
    return p2 ? 1 : 0;
  if(!p2)
    return -1;

  return CHashMatchIncId1(p1,p2);
}

class CReadBuffer {
  FILE *fp;
  char *fbuf;
  CHashMatch *matches;
#if STREAM_DEBUG>=2
public:
#endif
  char *filename;
  size_t matchcnt, matchmax, matchnext;/* match[matchnext] is the next unread element. If matchnext >= matchcnt, there are no more elements to read (EOF) */
  size_t totcnt;/* count of number of matches read so far (not including the current matchcnt matches in matches[] */
#if STREAM_DEBUG<2
public:
#endif
  CReadBuffer(char *bufname) {
    matchmax = HASHBUF_SIZ;
    size_t fbufsiz = sizeof(CHashMatch) * matchmax;
    fbuf = new char[fbufsiz];
    matches = new CHashMatch[matchmax];
    filename = bufname;
    if((fp = fopen(bufname,"r"))==NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("failed to open file '%s' for reading:errno=%d:%s\n",filename,eno,err);
      fflush(stdout); exit(1);
    }
    setbuffer(fp,fbuf,fbufsiz);

    /* trigger initial read */
    matchcnt = matchmax;
    matchnext = matchcnt - 1;
    totcnt = -matchcnt;
    next();
  }
  ~CReadBuffer(){
    if(fp) fclose(fp);
    delete [] matches;
    if(fbuf) delete [] fbuf;
  };
  inline CHashMatch *Match() {return (matchnext < matchcnt) ? &matches[matchnext] : NULL;}; /** returns pointer to next Match in buffer (or NULL if none)*/
  inline void next();/** advance buffer pointer by 1 (and trigger next file read, if needed) */
};

inline void CReadBuffer::next()
{
  if(matchnext >= matchcnt)/* already at EOF */
    return;
  if(++matchnext >= matchcnt){/* buffer is empty : try to read in more data */
    if(matchcnt < matchmax)/* implies the last read reached EOF */
      return;
    totcnt += matchcnt;
    matchcnt = fread(matches,sizeof(CHashMatch), matchmax, fp);
    if(matchcnt < matchmax){/* EOF or error */
      if(!feof(fp) && ferror(fp)){
	printf("read of file %s failed\n",filename);
	fflush(stdout); exit(1);
      }
      delete [] fbuf;
      fbuf = 0;
      fclose(fp);
      fp = NULL;
    }
    matchnext = 0;

    if(DEBUG>=2) /* check if Matches are in sorted order */
      for(size_t i = matchnext+1; i < matchcnt; i++)
	assert(CHashMatchInc(&matches[i-1],&matches[i]) <= 0);
  }
}

/* Below heap[] is a Min-Heap wrt the ordering of the next Match() in each buffer : 
   ordering is wrt (id1,id2,orientation) in ascending order except that a NULL Match() value (EOF) is treated as the largest possible  */

/* Pre : heap[1..heapcnt] is a Min-Heap exept for heap[1] being (possibly) too large
   Post : heap[1.heapcnt] is a Min-Heap.
   Children of heap[i] are heap[2i] and heap[2i+1]
*/
static void sift_down(CReadBuffer* heap[], int heapcnt)
{
  int i = 1, j;
  CReadBuffer* heapI = heap[i];
  CHashMatch *MatchI = heapI->Match();
  while((j = 2*i) <= heapcnt){
    if(j < heapcnt && CHashMatchInc(heap[j]->Match(),heap[j+1]->Match()) > 0){
      if(CHashMatchInc(MatchI, heap[j+1]->Match()) <= 0)
	break;
      heap[i] = heap[j+1];
      i = j+1;
    } else {
      if(CHashMatchInc(MatchI, heap[j]->Match()) <= 0)
	break;
      heap[i] = heap[j];
      i = j;
    }
  }
  heap[i] = heapI;
}

/*
   Pre : heap[1..heapcnt] is a Min-Heap except for heap[heapcnt] being (possibly) too small
   Post : heap[1..heapcnt] is a Min-Heap.

   Parent of heap[j] is heap[j/2] 
*/
static void sift_up(CReadBuffer* heap[], int heapcnt)
{
  int i = heapcnt;
  CReadBuffer* heapI = heap[i];
  CHashMatch* MatchI = heapI->Match();
  for(int j; i > 1; i = j){
    j = i/2;
    if(CHashMatchInc(MatchI , heap[j]->Match()) >= 0)
      break;
    heap[i] = heap[j];
  }
  heap[i] = heapI;
}

static CReadBuffer **bufferheap = NULL;
static CHashMatch *matches = NULL;
static size_t matchcnt = 0, matchmax = 0;
static CHashMatch lastmatch;/* for debugging */

/* unlink per-thread file buffers and free merge buffer memory */
static void hash_free()
{
  if(!(HashGen && DeferMerge))/* delete hashtable buffer files */
    for(int i=0; i < numbufs; i++)
      if(bufnames[i]){
	if(VERB>=2){
	  printf("hash_free:unlink(%s)\n",bufnames[i]);
	  fflush(stdout);
	}
	unlink(bufnames[i]);
      }

  delete [] matches;matches = NULL;
  matchcnt = matchmax = 0;
		      
  for(int i=0; i < numbufs; i++)
    if(bufnames[i])
      free(bufnames[i]);
  delete [] bufnames; bufnames = NULL;
  numbufs = hnumthreads = 0;

  delete [] besthashscore; besthashscore = NULL;
}

static void hash_outputF(char *filename, char **bufnames, int numbufs, long long TotalMatchCnt, int pairwise, Cmap *rmap[], int numrmaps, Cmap *qmap[], int numqmaps)
{
  checkFile(filename);

  FILE *fp;
  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("failed to open file %s for writing:errno=%d:%s\n",filename,eno,err);
    fflush(stdout); exit(1);
  }

  if(DeferMerge){/* just output the list of thread buffer names */

    printversion(fp);
    fprintf(fp,"# List of HashTable files (To be Merge Sorted):\n");
    for(int i = 0; i < numbufs; i++)
      fprintf(fp,"%s\n",bufnames[i]);
    (void)FILEclose(fp);    
    
    return;
  }

  if(VERB){
    printf("Generating %s (%lld Matches) by merge sorting %d buffers: CPU time=%0.6f, wall time=%0.6f secs\n",
	     filename, TotalMatchCnt, numbufs, mtime(), wtime());
    fflush(stdout);
  }
  
  if(HashTxt){
    printversion(fp);    /* write out commandline & header line */
    fprintf(fp,"# ID1 ID2 Orientation2 Score Offset(Left end of ID2 relative to Left end of ID1 in kb)\n");
  }

  /* use a heap of buffers to perform a merge sort */
  bufferheap = new CReadBuffer*[numbufs+1];
  for(int i = 1; i <= numbufs; i++)
    bufferheap[i] = new CReadBuffer(bufnames[i-1]);/* triggers the first buffer read */

  /* arrange bufferheap[1..numbufs] as a Min-Heap */
  for(int i = 2; i <= numbufs; i++)
    sift_up(bufferheap, i);

  /* Create output buffer for binary file */
  matchcnt = 0;
  matchmax = 512*1024;/* 8.0 Mbytes output buffer */
  if(!HashTxt) matches = new CHashMatch[matchmax];
  size_t fbufsiz = sizeof(CHashMatch) * matchmax;
  char *fbuf = new char[fbufsiz];
  setbuffer(fp,fbuf,fbufsiz);

  if(DEBUG>=2){
    lastmatch.id1 = -1;
    lastmatch.id2 = -1;
    lastmatch.orientation = 0;
  }

  /* Merge Sort */
  long long cnt = 0;
  for(CHashMatch *p; (p = bufferheap[1]->Match()) != NULL;){
    if(HashTxt) {
      if(MAP_INDEX_TXT)
	fprintf(fp,"%lld %lld %d %d %d\n", (long long)p->id1, (long long)p->id2, p->orientation ? -1 : 1, p->hashscore, p->offset);
      else
	fprintf(fp,"%lld %lld %d %d %d\n", qmap[p->id1]->id, rmap[p->id2]->id, p->orientation ? -1 : 1, p->hashscore, p->offset);
    } else {
      if(DEBUG>=2){
	CHashMatch *q = (matchcnt==0) ? &lastmatch : &matches[matchcnt-1];
	assert(q->id1 <= p->id1);
	if(q->id1 == p->id1){
	  assert(q->id2 <= p->id2);
	  if(q->id2 == p->id2)
	    assert(q->orientation < p->orientation);
	}
      }
      matches[matchcnt++] = *p;

      if(matchcnt >= matchmax){/* write output buffer */
	if(DEBUG>=2) lastmatch = matches[matchcnt-1];
	size_t n = fwrite_LOOP(matches,sizeof(CHashMatch),matchcnt,fp);
	if(n != matchcnt){
	  #pragma omp critical
	  {
	    printf("fwrite=%lu: write of %lu Matches (%lu bytes each) to %s failed\n",
		 n, matchcnt, sizeof(CHashMatch), filename);
	    fflush(stdout);exit(1);
	  }
	}
	matchcnt = 0;
      }
    }

    /* update bufferheap[] : may involve an actual file read, if needed */
    bufferheap[1]->next();
    sift_down(bufferheap,numbufs);
    
    if(VERB && cnt && !(cnt % 1000000000)){
      printf("Merge-Sorted %lld Billion Matches: wall time=%0.6f secs\n",cnt/1000000000, wtime());
      fflush(stdout);
    }
    cnt++;
  }
  if(matchcnt > 0){/* write last output buffer */
    size_t n = fwrite_LOOP(matches,sizeof(CHashMatch),matchcnt,fp);
    if(n != matchcnt){
      #pragma omp critical
      {
	printf("fwrite=%lu: write of %lu Matches (%lu bytes each) to %s failed\n",
	       n, matchcnt, sizeof(CHashMatch), filename);
	fflush(stdout);exit(1);
      }
    }
  }

  (void)FILEclose(fp);
  delete [] fbuf;

  if(VERB){
    printf("Wrote %lld matches to %s: CPU time=%0.6f, wall time=%0.6f secs\n", cnt, filename, mtime(),wtime());
    fflush(stdout);
  }

  /* free up input buffers */
  for(int i = 1; i <= numbufs; i++)
    delete bufferheap[i];
  delete [] bufferheap; bufferheap = NULL;
  if(!HashTxt){
    delete [] matches; matches = NULL;
    matchcnt = matchmax = 0;
  }
}

/** merge sort temporary output buffers into a single hashtable file.  If SubsetMap && !pairwise : Output subset of rmap[0..numrmaps-1] with hashtable matches to <prefix>_hash.bnx */
static void hash_output(char *prefix, char **bufnames, int numbufs, long long TotalMatchCnt, int pairwise, Cmap *rmap[], int numrmaps, Cmap *qmap[], int numqmaps)
{
  char filename[PATH_MAX];

  if(HashTxt)
    sprintf(filename,"%s.hash",prefix);
  else
    sprintf(filename,"%s.hashbin",prefix);

  if(VERB>=2){
    printf("Calling hash_outputF(%s,numbufs=%d)\n",filename,numbufs);
    fflush(stdout);
  }

  hash_outputF(filename, bufnames, numbufs, TotalMatchCnt, pairwise, rmap, numrmaps, qmap, numqmaps);

  if(BNX_OUTPUT && !pairwise){/* output <prefix>_hash.bnx of rmap[] subset with hashtable matches */
    Cmap **goodmaps = new Cmap *[numrmaps];
    int numgoodmaps = 0;
    for(int rid = 0; rid < numrmaps; rid++)
      if(MapSave[rid])
	goodmaps[numgoodmaps++] = rmap[rid];
    if(VERB>=2){
      printf("numgoodmaps=%d,numrmaps=%d\n",numgoodmaps,numrmaps);
      fflush(stdout);
    }
    sprintf(filename,"%s_hash",prefix);
    output_bnx(filename,goodmaps,0,numgoodmaps-1,0,NULL,NULL,-1);

    /* output mapid mapping */
    sprintf(filename,"%s.hashmap",prefix);
    checkFile(filename);
    if(VERB){
      printf("Generating %s (mapid mapping of %d maps)\n", filename, numgoodmaps);
      fflush(stdout);
    }
    FILE *fp;
    if((fp = fopen(filename,"w")) == NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("failed to open file %s for writing:errno=%d:%s\n",filename,eno,err);
      fflush(stdout); exit(1);
    }
    printversion(fp);
    fprintf(fp,"# List of hashtable mapids of maps in %s_hash.bnx:\n",prefix);
    fprintf(fp,"# Index Mapid\n");
    for(int i = 0; i < numgoodmaps; i++)
      fprintf(fp,"%d %d\n",i,goodmaps[i]->mapid);
    (void)FILEclose(fp);

    delete [] goodmaps;
  }
  return;
}

static char buf[LINESIZ];

static int MergeSort = 0;

/** Open a hashtable file for streaming input with hash_read() : If hashtable file is encoded as a list of files, Merge Sort is also performed as the data streams in */
/** Returns the total number of hash matches present in the hashtable file(s) */
size_t hash_open(FILE* &HashFP, char *hash_filename)
{
  HashFP = NULL;

  if(HASH_STREAM>=3 && hash_filename && !HashTxt){/* stream merge sort output */

    MergeSort = HashGen;
    if(!HashGen){
      /* check if hash_filename is encoded as a list of filenames : If so set bufnames[] to those filenames and set MergeSort=1 */
      if(!HashFP && (HashFP = fopen(hash_filename, "r"))==NULL){
	int eno = errno;
	char *err = strerror(errno);
	printf("failed to open file %s for reading:errno=%d:%s\n", hash_filename, eno, err);
	fflush(stdout); exit(1);
      }
      if(1 != fread(buf,1,1,HashFP)){/* EOF or error */
	int myerrno = errno;
	if(feof(HashFP)){
	  printf("Unable to read first character from file %s(EOF)\n",hash_filename);
	  fflush(stdout); exit(1);
	}
	if(ferror(HashFP)){
	  printf("read error trying to read first character from file %s:%s\n",hash_filename, strerror(myerrno));
	  fflush(stdout); exit(1);
	}
      }
      if(buf[0] == '#'){
	if(VERB && STREAM_DEBUG){
	  printf("%s is encoded as a list of files: wall time=%0.6f\n",hash_filename,wtime());
	  fflush(stdout);
	}
	fseek(HashFP, 0L, SEEK_SET);
	if(numbufs > 0)
	  for(int i = 0; i < numbufs; i++)
	    free(bufnames[i]);
	else {
	  if(bufnames) delete [] bufnames;
	  numbufs = hnumthreads;
	  bufnames = new char*[numbufs];
	}
	int bufcnt = 0;
	int linecnt = 1;
	FILE *fp;
	TotalMatchCnt = 0;
	for(;fgets(buf,BUFSIZ,HashFP) != NULL; linecnt++){
	  long len = strlen(buf);
	  if(len<1)continue;
	  if(buf[0]=='#' || buf[0]=='\n' || buf[0]=='\r')continue;
	  if(len >= BUFSIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	    printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",
		   buf[len-1],hash_filename,linecnt,buf);
	    fflush(stdout); exit(1);
	  }
	  buf[len-1] = '\0';
	  if(bufcnt >= numbufs){
	    char **oldbufnames = bufnames;
	    numbufs *= 2;
	    bufnames = new char*[numbufs];
	    memcpy(bufnames,oldbufnames,bufcnt*sizeof(char*));
	    delete [] oldbufnames;
	  }
	  bufnames[bufcnt++] = strdup(buf);

	  /* determine file length to update TotalMatchCnt (and make sure file exists) */
	  if((fp = fopen(buf,"r")) == NULL){
	    int eno = errno;
	    char *err = strerror(errno);
	    printf("failed to open file %s for reading:errno=%d:%s\n",buf,eno,err);
	    fflush(stdout); exit(1);
	  }
	  fseek(fp, 0L, SEEK_END);
	  len = ftell(fp);
	  if((len % sizeof(CHashMatch)) != 0){
	    printf("file %s length=%ld bytes is NOT a multiple of CHashMatch = %lu bytes\n",
		   buf, len, sizeof(CHashMatch));
	    fflush(stdout); exit(1);
	  }
	  TotalMatchCnt += len / sizeof(CHashMatch);
	  fclose(fp);
	}
	numbufs = bufcnt;
	if(VERB){
	  printf("Read in (from %s) %d buffer filenames (which contain %lld Total Matches):wall time=%0.6f\n",hash_filename,bufcnt,TotalMatchCnt,wtime());
	  fflush(stdout);
	}
	fclose(HashFP);

	MergeSort = 1;
      }
    }

    if(MergeSort){
      /* use a heap of buffers to perform a merge sort */
      if(DEBUG) assert(bufnames != NULL);
      bufferheap = new CReadBuffer*[numbufs+1];
      for(int i = 1; i <= numbufs; i++)
	bufferheap[i] = new CReadBuffer(bufnames[i-1]);/* triggers the first buffer read */
      
      /* arrange bufferheap[1..numbufs] as a Min-Heap */
      for(int i = 2; i <= numbufs; i++)
	sift_up(bufferheap, i);
      
      if(DEBUG>=2){
	lastmatch.id1 = -1;
	lastmatch.id2 = -1;
	lastmatch.orientation = 0;
      }
      
      if(VERB){
	printf("Streaming %lld Matches by merge sorting %d buffers: CPU time=%0.6f, wall time=%0.6f secs\n",
	       TotalMatchCnt, numbufs, mtime(), wtime());
	fflush(stdout);
      }
      
      HashFP = (FILE*)(-1);
      return (size_t)TotalMatchCnt;
    }
  } 

  /* stream file input */
  if(!HashFP && (HashFP = fopen(hash_filename, "r"))==NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("failed to open file %s for reading:errno=%d:%s\n", hash_filename, eno,err);
    fflush(stdout); exit(1);
  }

  /* determine file length : may not work under Windows with 32-bit long */
  fseek(HashFP, 0L, SEEK_END);
  long len = ftell(HashFP);
  fseek(HashFP, 0L, SEEK_SET);
  
  if((len % sizeof(CHashMatch)) != 0){
    printf("file %s length=%ld bytes is not a multiple of CHashMatch = %lu bytes\n",
	   hash_filename,len,sizeof(CHashMatch));
    fflush(stdout); exit(1);
  }
  size_t MatchMax = len/sizeof(CHashMatch);
  if(VERB){
    printf("Streaming in %lu Matches from %s ...\n",MatchMax,hash_filename);
    fflush(stdout);
  }

  return MatchMax;
}

/** Check if hashtable input stream has reached the end */
int hash_eof(FILE *HashFP)
{ 
  if(HASH_STREAM>=3 && MergeSort && hash_filename && !HashTxt)
    return (bufferheap[1]->Match() == NULL) ? 1 : 0;

  return feof(HashFP);
}

/** Close a hashtable input stream. NOTE : hash_close() discards(unlinks) the input data files IFF final==1 : If the HashTable data will be required again use final=0 */
void hash_close(FILE *HashFP, int final)
{
  if(HASH_STREAM>=3 && MergeSort && hash_filename && !HashTxt){
    if(VERB>=2){
      printf("hash_close:final=%d\n",final);
      fflush(stdout);
    }

    /* free up input buffers */
    for(int i = 1; i <= numbufs; i++)
      delete bufferheap[i];
    delete [] bufferheap;

    if(final)
      hash_free();/* delete buffer names and files and matches[] */
  } else
    (void) fclose(HashFP);
}

/** read hashtable input stream : trys to get the specified amount of hashtable matches. Returns the actual number obtained (if less than requested, EOF has been reached and hash_eof() will return 1) */
size_t hash_read(FILE *fp, int txt, CHashMatch* matches, size_t matchmax, int &linecnt)
{
  if(!txt){/* binary file */
    if(HASH_STREAM>=3 && MergeSort && hash_filename){/* Continue Merge Sort */
      if(VERB>=2){
	static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
	getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
	printf("hash_read: starting read of up to %lu matches, hashbest=%d: VmSize= %0.4f VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb\n",
	       matchmax,HashBest,VmSize*1e-9,VmRSS*1e-9,VmHWM*1e-9,VmSwap*1e-9);
	fflush(stdout);
      }

      size_t matchcnt = 0, totcnt = 0;
      for(CHashMatch *p; (p = bufferheap[1]->Match()) != NULL;){
#if STREAM_DEBUG>=2 // requires most fields in CReadBuffer to be public
	if(VERB){
	  if(HashBest >= 0)
	    printf("%s,cnt=%lu/%lu:id1=%d,id2=%d,orientation=%d,hashscore=%d,offset=%d: besthashscore[id2]= %d, HashBest=%d, matchcnt=%lu%s\n",
		   bufferheap[1]->filename, bufferheap[1]->totcnt + bufferheap[1]->matchnext, bufferheap[1]->totcnt + bufferheap[1]->matchcnt,
		   p->id1,p->id2,p->orientation,p->hashscore,p->offset, besthashscore[p->id2], HashBest, matchcnt, (HashBest >= 0 && p->hashscore < besthashscore[p->id2] - HashBest) ? ": Discarded" : "");
	  else
	    printf("%s,cnt=%lu/%lu:id1=%d,id2=%d,orientation=%d,hashscore=%d,offset=%d: HashBest=%d, matchcnt=%lu%s\n",
		   bufferheap[1]->filename, bufferheap[1]->totcnt + bufferheap[1]->matchnext, bufferheap[1]->totcnt + bufferheap[1]->matchcnt,
		   p->id1,p->id2,p->orientation,p->hashscore,p->offset, HashBest, matchcnt, (HashBest >= 0 && p->hashscore < besthashscore[p->id2] - HashBest) ? ": Discarded" : "");
	  fflush(stdout);
	}
#endif
	if(DEBUG>=2 && Gpairwise){
	  CHashMatch *q = (matchcnt==0) ? &lastmatch : &matches[matchcnt-1];
	  assert(q->id1 <= p->id1);
	  if(q->id1 == p->id1){
	    assert(q->id2 <= p->id2);
	    if(q->id2 == p->id2){
	      if(HashMultiMatch)
		assert(q->orientation <= p->orientation);
	      else
		assert(q->orientation < p->orientation);
	    }
	  }
	}
	totcnt++;
	if(!(HashBest >= 0 && p->hashscore < besthashscore[p->id2] - HashBest))
	  matches[matchcnt++] = *p;

	/* update bufferheap[] : may involve an actual file read, if needed */
	bufferheap[1]->next();
	sift_down(bufferheap,numbufs);

	if(matchcnt >= matchmax){
	  if(DEBUG>=2 && Gpairwise) lastmatch = matches[matchcnt-1];
	  if(VERB /* HERE HERE && STREAM_DEBUG */){
	    static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
	    getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
	    printf("hash_read: read %lu/%lu matches, hashbest=%d: VmSize= %0.4f VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb\n",
		   matchcnt,totcnt,HashBest,VmSize*1e-9,VmRSS*1e-9,VmHWM*1e-9,VmSwap*1e-9);
	    fflush(stdout);
	  }
	  return matchcnt;
	}
      }
      if(VERB /* HERE HERE && STREAM_DEBUG */){
	static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
	getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
	printf("hash_read: read %lu/%lu matches, hashbest=%d: VmSize= %0.4f VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb\n",
	       matchcnt,totcnt,HashBest,VmSize*1e-9,VmRSS*1e-9,VmHWM*1e-9,VmSwap*1e-9);
	fflush(stdout);
      }
      return matchcnt;
    }

    size_t matchcnt = fread(matches, sizeof(CHashMatch), matchmax, fp);
    int myerrno = errno;
    if(matchcnt != matchmax){/* EOF or error */
      if(!feof(fp) && ferror(fp)){
	printf("hash_read: read error while trying to read %lu bytes:%s\n",matchmax * sizeof(CHashMatch),strerror(myerrno));
	fflush(stdout); exit(1);
      }
    }
    return matchcnt;
  }

  /* read text file */
  size_t matchcnt = 0;

  for(;matchcnt < matchmax && fgets(buf,LINESIZ,fp) != NULL;linecnt++){
    int len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],hash_filename,linecnt,buf);
      fflush(stdout); exit(1);
    }
    if(buf[0] == '#')
      continue;
    CHashMatch *phash = &matches[matchcnt++];
    
    char *pt = buf,*qt;
    phash->id1 = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || phash->id1 < 0){
      printf("Invalid integer value %lld for Mol ID1 on line %d of %s:\n%s\n",(long long)phash->id1, linecnt,hash_filename,buf);
      fflush(stdout); exit(1);
    }
    pt = qt;
    phash->id2 = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || phash->id2 < 0){
      printf("Invalid integer value %lld for Mol ID2 on line %d of %s:\n%s\n",(long long)phash->id2, linecnt,hash_filename,buf);
      fflush(stdout); exit(1);
    }
    pt = qt;
    long orientation = strtol(pt,&qt,10);
    if(pt==qt || !(*qt == 0 || isspace(*qt)) || labs(orientation) != 1){
      printf("Invalid integer value %ld for orientation on line %d of %s:\n%s\n",orientation,linecnt,hash_filename,buf);
      fflush(stdout); exit(1);
    }
    phash->orientation = (orientation < 0) ? 1 : 0;

    pt = qt;
    long hashscore = strtol(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || hashscore > MASK(16)){
      printf("Invalid integer value for hash score on line %d of %s (must be <= %d):\n%s\n",
	     linecnt,hash_filename,MASK(16),buf);
      fflush(stdout); exit(1);
    }
    phash->hashscore = hashscore;


    pt = qt;
    phash->offset = strtol(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt))){
      printf("Invalid integer value for offset (in Kb) between left end of maps on line %d of %s:\n%s\n",linecnt,hash_filename,buf);
      fflush(stdout); exit(1);
    }

    if(HashBest >= 0 && hashscore < besthashscore[phash->id2] - HashBest)
      matchcnt--;
  }

  return matchcnt;
}

/** file input of entire .hash or .hashbin file (not recommended for files larger than 1 Gbyte : use hash_open() & hash_read() instead to stream the input) */
void hash_input(char *hash_filename, size_t &matchcnt, CHashMatch* &matches)
{
  FILE *fp;
  if((fp = fopen(hash_filename,"r"))==NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("failed to open file %s for reading:errno=%d:%s\n", hash_filename,eno,err);
    fflush(stdout); exit(1);
  }

  if(VERB){
    printf("Reading in Matches from %s\n", hash_filename);
    fflush(stdout);
  }

  size_t L = strlen(hash_filename);
  if(L >= strlen(".hashbin") && !strcmp(&hash_filename[L-strlen(".hashbin")],".hashbin")){/* binary file input */

    /* check if file is encoded as a list of filenames : If so exit with error message (Not supported) */
    if(1 != fread(buf,1,1,fp)){/* EOF or error */
      int myerrno = errno;
      if(feof(fp)){
	printf("Unable to read first character from file %s(EOF)\n",hash_filename);
	fflush(stdout); exit(1);
      }
      if(ferror(fp)){
	printf("read error trying to read first character from file %s:%s\n",hash_filename, strerror(myerrno));
	fflush(stdout); exit(1);
      }
    }
    if(buf[0] == '#'){
      printf("hash_input: file %s begins with '#' : List of buffer names not supported\n",hash_filename);
      fflush(stdout); exit(1);
    }

    /* determine file length : may not work under Windows with 32-bit long */
    fseek(fp, 0L, SEEK_END);
    long len = ftell(fp);
    fseek(fp, 0L, SEEK_SET);
  
    if((len % sizeof(CHashMatch)) != 0){
      printf("file %s length=%ld bytes is NOT a multiple of CHashMatch = %lu bytes\n",
	     hash_filename,len,sizeof(CHashMatch));
      fflush(stdout); exit(1);
    }
    size_t matchmax = len/sizeof(CHashMatch);
    matches = new CHashMatch[matchmax+1];
    matchcnt = fread(matches, sizeof(CHashMatch), matchmax, fp);
    if(matchcnt != matchmax){
      printf("fread()=%lu : Failed to read %lu Matches (%lu bytes each) from %s\n",
	     matchcnt, matchmax, sizeof(CHashMatch), hash_filename);
      fflush(stdout); exit(1);
    }
  } else { /* text file */
    size_t matchmax = 0;
    int linecnt = 0;
    for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++){
      int len = strlen(buf);
      if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],hash_filename,linecnt,buf);
	fflush(stdout); exit(1);
      }
      if(buf[0] == '#')
	continue;
      maxmatchalloc(matchcnt,matchcnt+1,matchmax,matches);
      CHashMatch *phash = &matches[matchcnt++];
    
      char *pt = buf,*qt;
      phash->id1 = strtoll(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt)) || phash->id1 < 0){
	printf("Invalid integer value %lld for Mol ID1 on line %d of %s:\n%s\n",(long long)phash->id1, linecnt,hash_filename,buf);
	fflush(stdout); exit(1);
      }
      pt = qt;
      phash->id2 = strtoll(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt)) || phash->id2 < 0){
	printf("Invalid integer value %lld for Mol ID2 on line %d of %s:\n%s\n",(long long)phash->id2, linecnt,hash_filename,buf);
	fflush(stdout); exit(1);
      }
      pt = qt;
      long orientation = strtol(pt,&qt,10);
      if(pt==qt || !(*qt == 0 || isspace(*qt)) || labs(orientation) != 1){
	printf("Invalid integer value %ld for orientation on line %d of %s:\n%s\n",orientation,linecnt,hash_filename,buf);
	fflush(stdout); exit(1);
      }
      phash->orientation = (orientation < 0) ? 1 : 0;

      pt = qt;
      long hashscore = strtol(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt)) || hashscore > MASK(16)){
	printf("Invalid integer value for hash score on line %d of %s (must be <= %d):\n%s\n",
	       linecnt,hash_filename,MASK(16),buf);
	fflush(stdout); exit(1);
      }
      phash->hashscore = hashscore;

      pt = qt;
      phash->offset = strtol(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt))){
	printf("Invalid integer value for offset (in Kb) between left end of maps on line %d of %s:\n%s\n",linecnt,hash_filename,buf);
	fflush(stdout); exit(1);
      }
    }
  }
  
  if(VERB){
    printf("Read in %lu Matches from %s : CPU time=%0.6f, wall time=%0.6f secs\n", 
	   matchcnt, hash_filename, mtime(), wtime());
    fflush(stdout);
  }
  fclose(fp);

  if(DEBUG>=2) /* check if Matches are in sorted order */
    for(size_t i = 1; i < matchcnt; i++)
      assert(CHashMatchInc(&matches[i-1],&matches[i]) <= 0);
}

/** Top Level function to generate a hashtable File.
hash_generate() - first generate the hash table using ref maps, then query it using query maps (NOTE: 'ref' and 'query' maps have their meaning reversed compared to alignment)
A summary of the routine:
#### Hash Phase
     - hash_init() - initialize CHashTable, key functions, constants, tables.
     - CHashTable::hashinsert()
        + for each map : CHashTable::hash_insert()
             + for each site on map : CHashTable::keyinsert() - apply hash function to window and index in CHashTable::index.
        - CHashTable::hashsort() - Compress CHashTable::index into CHashTable::indexR, sort indices.
____
#### Query Phase
     - for each query map : CMatchTable::hash_query() - find the best hit offset/score for each ref map
        + for each orientation : CMatchTable::hash_queryF() - compute sizes & variances for window
             - for each site on map: CHashTable::hashfind() - query CHashTable::indexR and return hits 
             - check sizing errors (computationally expensive)
             - (if pass) CMatchTable::HitInsert() - register hits in CMatchTable::HitIndex indexed by ref map and ref site.
             - CMatchTable::HitToMatch() - traverse CMatchTable::HitIndex, remove duplicates and accumulate scores in CMatchTable::MatchIndex
        - CMatchTable::MatchFind() - traverse CMatchTable::MatchIndex average match scores over nearby offsets, and pass them to CMatchTable::BestMatch and output highest scoring match per ref map, if above HashScore.
	
*/

template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
void hash_generateT(Cmap *rmap[], int numrmaps, Cmap *qmap[], int numqmaps,  int pairwise, char *prefix)
{
  

  if(VERB){
    printf("Start of hash_generate(): CPU time=%0.6f, wall time=%0.6f secs\n",mtime(),wtime());
    fflush(stdout);
  }

  if(DEBUG && !pairwise) assert(rmap != qmap);

  splitcnt = 1;
  
  if(HMaxMem > 0.0 && (HASH_STREAM>=3 && HashGen && hash_filename && !HashTxt) && !pairwise){ /* check if rmap[] needs to be split into sections to be processed sequentially : up to 100k maps or 3m labels per 8 Gb of -maxmem */
    long long totsites = 0;
    for(int i = numrmaps; --i >= 0;)
      totsites += rmap[i]->numsite[hcolor];

    splitcnt = max(ceil((totsites * 6.0) / (3e6 * HMaxMem)), ceil((numrmaps * 6.0) / (100e3 * HMaxMem)));
    splitcnt = max(1,splitcnt);

    if(VERB && splitcnt > 1){
      printf("processing %d maps with %lld sites in %d stages to stay within %0.1f Gb of memory\n", numrmaps,totsites, splitcnt,HMaxMem);
      fflush(stdout);
    }
  }

  orig_numrmaps = numrmaps;

  int split_numbufs = 0;
  char **split_bufnames = NULL;
  char *orig_prefix = prefix;
  maxMatchesPerRefid = maxMapsPerRefid = 0;
  if(splitcnt > 1){
    prefix = new char[strlen(orig_prefix) + 256];
    Gnumqmaps = numqmaps;
    delete [] MatchesPerRefid;
    MatchesPerRefid = new size_t[numqmaps];
    delete [] MapsPerRefid;
    MapsPerRefid = new size_t[numqmaps];
    for(int i = numqmaps; --i >= 0;)
      MatchesPerRefid[i] = MapsPerRefid[i] = 0;
  }

  cumTotalMatchCnt = 0;/* sum of TotalMatchCnt over all splits */

  if(colors >= 2){
    if(!hashcolor){
      printf("hashtable with colors=%d requires use of -hashcolor\n",colors);
      fflush(stdout); exit(1);
    }
  }
  if(hashcolor)
    hcolor = hashcolor - 1;

  if(SF[hcolor] < max(0.0,MINSF)){
    printf("SF[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,SF[hcolor],MINSF,MAXSF,MINSF);
    SF[hcolor] = max(0.0,MINSF);
  } else if(SF[hcolor] > MAXSF){
    printf("SF[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,SF[hcolor],MINSF,MAXSF,MAXSF);
    SF[hcolor] = MAXSF;
  }
  if(SR[hcolor] < max(0.0,MINSR)){
    printf("SR[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,SR[hcolor],MINSR,MAXSR,MINSR);
    SR[hcolor] = max(0.0,MINSR);
  } else if(SR[hcolor] > MAXSR){
    printf("SR[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,SR[hcolor],MINSR,MAXSR,MAXSR);
    SR[hcolor] = MAXSR;
  }
  if(SD[hcolor] < MINSD){
    printf("SD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,SD[hcolor],MINSD,MAXSD,MINSD);
    SD[hcolor] = MINSD;
  } 
  if(SD[hcolor] > MAXSD){
    printf("SD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,SD[hcolor],MINSD,MAXSD,MAXSD);
    SD[hcolor] = MAXSD;
  }
  if(SD[hcolor] < -sqrt(2.0*SF[hcolor]*SR[hcolor])){
    if(SR[hcolor] <= 0.0)
      printf("ERROR:c=%d:SD[c]=%0.6f is not valid for SF[c]=%0.6f,SR[c]=%0.6f: changing SD[c] to %0.6f (must specify +ve sr with -ve sd)\n",hcolor,SD[hcolor],SF[hcolor],SR[hcolor],-sqrt(2.0*SF[hcolor]*SR[hcolor]));
    else
      printf("WARNING:c=%d:SD[c]=%0.6f is not valid for SF[c]=%0.6f,SR[c]=%0.6f: changing SD[c] to %0.6f\n",hcolor,SD[hcolor],SF[hcolor],SR[hcolor],-sqrt(2.0*SF[hcolor]*SR[hcolor]));
    fflush(stdout);
    assert(SR[hcolor] >= 0.0);
    SD[hcolor] = -sqrt(2.0*SF[hcolor]*SR[hcolor]);
  }
#if 0 // SE[] not yet used by hashtable
  if(SE[hcolor] < MINSE){
    printf("SE[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,SE[hcolor],MINSE,MAXSE,MINSE);
    SE[hcolor] = MINSE;
  } 
  if(SE[hcolor] > MAXSE){
    printf("SE[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,SE[hcolor],MINSE,MAXSE,MAXSE);
    SE[hcolor] = MAXSE;
  }
#endif
  if(res[hcolor] < MinRes){
    printf("res[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,res[hcolor],MinRes,MaxRes,MinRes);
    res[hcolor] = MinRes;
  } else if(res[hcolor] > MaxRes){
    printf("res[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,res[hcolor],MinRes,MaxRes,MaxRes);
    res[hcolor] = MaxRes;
  }
  if(resSD[hcolor] < MinResSD){
    printf("resSD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,resSD[hcolor],MinResSD,MaxResSD,MinResSD);
    resSD[hcolor] = MinResSD;
  } else if(resSD[hcolor] > MaxResSD){
    printf("resSD[%d]=%0.6f is not in range [%0.6f .. %0.6f] : changing to %0.6f\n",hcolor,resSD[hcolor],MinResSD,MaxResSD,MaxResSD);
    resSD[hcolor] = MaxResSD;
  }

  double origSF = SF[hcolor];
  double origSD = SD[hcolor];
  double origSR = SR[hcolor];
  double origSDmax = SDmax;

  if(!HASH_SCALE && (ScaleDelta > 0 || (MapScale && MapScaleDelta > 0))){  // a tempory hack, until HASH_SCALE is implemented in the hashtable
    SF[hcolor] = max(SF[hcolor], (giter2 < RefRepeats2 - 3) ? 0.10 : 0.00);
    SD[hcolor] = max(SD[hcolor], 0.0);
    SR[hcolor] = max(SR[hcolor], (giter2 < RefRepeats2 - 3) ? 0.03 : 0.02);
  }
  
  // Use hashtable specific upper bounds to SF[] and SR[], to help reduce hashtable hits, without affecting alignments
  SF[hcolor] = max(origSF* HashSFscale, HashSF);
  SD[hcolor] = max(origSD* HashSDscale, HashSD);
  SR[hcolor] = max(origSR* HashSRscale, HashSR);

  for(int split = 0; split < splitcnt; split++){
    int splitmaps = orig_numrmaps / splitcnt;
    int start = split * splitmaps;
    int end = (split == splitcnt-1) ? orig_numrmaps : splitmaps * (split + 1);

    numrmaps = end - start;
    roffset = start;

    if(splitcnt > 1)
      sprintf(prefix, "%s_S%d", orig_prefix, split+1);
    
    if(VERB && splitcnt > 1){
      printf("split= %d/%d : processing rmap[%d..%d], output prefix= %s\n",split,splitcnt,start,end-1, prefix);
      fflush(stdout);
    }

    /* START OF ORIGINAL CODE WITHOUT splitcnt */

    hnumthreads = 1;/* global value used in hash.cpp, hash.h only */

#ifdef _OPENMP
    hnumthreads = omp_get_max_threads();
  
    if(HashQueryThreads > 0)
      hnumthreads = min(hnumthreads,HashQueryThreads);
    else
      hnumthreads = min(hnumthreads,MaxThreads);

    if(HashMatrix)
      hnumthreads = 1;

    HashInsertThreads = min(hnumthreads, HashInsertThreads);

    if(VERB){
      printf("Using %d threads (%d for map insertion)\n",hnumthreads, HashInsertThreads);
      fflush(stdout);
    }
#else
    HashInsertThreads = 1;
#endif // _OPENMP


    /* initialize hashtable and related tables and constants */
    hash_init<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>(rmap,numrmaps,qmap,numqmaps,pairwise,prefix);

    CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *hashtable = new CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>(rmap,numrmaps,qmap,numqmaps,0);

    if(VERB){
      printf("HashTable initialized:CPU time=%0.6f, wall time=%0.6f secs\n",mtime(),wtime());
      /*    printf("HASH_BITS=%d,MHASH_BITS=%d,HHASH_BITS=%d,CollideBlockMax=%u,EntryBlockMax=%u\n",
	    HASH_BITS,MHASH_BITS,HHASH_BITS,hashtable->CollideBlockMax,hashtable->EntryBlockMax);*/
      fflush(stdout);
    }

    hashtable->hashinsert(rmap,numrmaps);/* insert all maps rmap[0..numrmaps-1] into hashtable */

    if(VERB>=2){
      printf("HashTable insertion complete: sleeping for 60 seconds\n");
      fflush(stdout);
      sleep(60);
    }

    if(MULTIMATCH_TWOTHREADS && (hnumthreads % 2))
      hnumthreads++;

    if(hnumthreads > (MULTIMATCH_TWOTHREADS ? 2 : 1) * max(1,numqmaps)){
      if(VERB){
	printf("Reducing -queryThreads from %d to %d due to only %d Query Maps\n",hnumthreads, (MULTIMATCH_TWOTHREADS ? 2 : 1) * numqmaps, numqmaps);
	fflush(stdout);
      }
      hnumthreads = (MULTIMATCH_TWOTHREADS ? 2 : 1) * max(1,numqmaps);
    }

    if(HMaxMem > 0){  /* check if we need to reduce numthreads to satisfy -maxmem */
      long long hashmem = hashtable->CollideCnt*sizeof(CHashCollideR) + hashtable->EntryCnt*sizeof(CHashEntry<OFFSET1_T>) + sizeof(CHashIndexR)*(1LL << HASH_BITS) + hashtable->HashSizeRtab();
 
      long long matchmem = (1LL << HHASH_BITS)*(sizeof(int)+sizeof(CHitEntry<OFFSET1_T>)) + (1LL << MHASH_BITS)*(sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * 5LL)/4LL;

      /* The following memory estimate for MatchBuf[] is not used by all threads at the same time, so is discounted 5x (except for pairwise) vs the worst case scenario : 
	 If worst cases happens, resource.h will pause some threads and swap out their MatchBuf[] memory */
      long long matchmem2 = ((1LL << MHASH_BITS) * MATCHLINE(OFFSET_T,RANGE,HASH_SCALE) * sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>) * (pairwise ? 5 : (1+2*OffsetSpread)) * 5LL)/(4LL * 5);

      int maxthreads = (HMaxMem * (1000LL * 1000LL * 1000LL) - hashmem) / (matchmem + matchmem2);
      maxthreads = max(1, maxthreads);// NEW7
      if(MULTIMATCH_TWOTHREADS && (maxthreads % 2))
	maxthreads++;
      
      if(maxthreads < hnumthreads){
	if(VERB){
	  printf("Reducing -queryThreads from %d to %d due to -maxmem %0.1f: HashTable= %lld bytes, HitTable+MatchIndex= %lld, MatchBuf= %lld (bytes per thread)\n",
		 hnumthreads, max(1,maxthreads), HMaxMem, hashmem, matchmem, matchmem2);
	  fflush(stdout);
	}
	hnumthreads = min(hnumthreads,max(1,maxthreads));
      }

      /* set amount of global memory available for MatchFind() allocations of MatchBuf[] and MatchOverflow[] reallocations */
      maxMatchMem = HMaxMem *(1000LL * 1000LL * 1000LL) - hashmem - matchmem * hnumthreads;

      if(VERB/* HERE >=2 */){
	printf("split=%d/%d:Preallocated:HashTable= %lld bytes, HitTable+MatchIndex= %lld bytes(for %d threads)\n",split,splitcnt,hashmem,matchmem*hnumthreads,hnumthreads);
	printf("\t Sharing %lld bytes of memory for MatchBuf[] and MatchOverflow[] between %d threads\n",maxMatchMem,hnumthreads);
	fflush(stdout);
      }

      maxMatchMem = max(1024LL * 1024LL * 128LL, maxMatchMem);
      hashtable->MatchMemReserve = new Resource<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>(maxMatchMem,hnumthreads);
    }

    if(MULTIMATCH_TWOTHREADS) {/* use 2 threads per hash_query() call */
      hnumthreads = (hnumthreads + 1)/2;
#ifdef _OPENMP    
      omp_set_nested(1);
      omp_set_dynamic(0);
#endif
    }

    numbufs = hnumthreads;
    if(splitcnt <= 1)
      bufnames  = new char* [numbufs];
    else {
      int prev_split_numbufs = split_numbufs;
      char **prev_split_bufnames = split_bufnames;

      split_numbufs += hnumthreads;
      split_bufnames = new char* [split_numbufs];
      for(int i = 0; i < prev_split_numbufs; i++)
	split_bufnames[i] = prev_split_bufnames[i];

      delete [] prev_split_bufnames;
	  
      bufnames = &split_bufnames[prev_split_numbufs];
    }

    CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> **HashTables = new CHashTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>*[hnumthreads];
    for(int i = 0; i < hnumthreads; i++)
      HashTables[i] = hashtable;

    if(VERB>=2){
      printf("split=%d/%d: Starting multithreaded hashtable probe with %d threads (pairwise=%d): cum wall time= %0.6f \n",split,splitcnt,hnumthreads, pairwise, wtime());
      fflush(stdout);
    }

    double mtime1 = mtime(), wtime1 = wtime();

#pragma omp parallel num_threads(hnumthreads) if(hnumthreads > 1)
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif

      CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *matchtable = new CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>(tid, HashTables[tid], numrmaps, rmap, numqmaps, qmap, prefix, NULL);
      bufnames[tid] = strdup(matchtable->filename);/* save name of output buffer used by this thread */
      if(VERB>=2){
	#pragma omp critical
	{
	  printf("split=%d/%d, tid=%d: output bufname= %s\n",split,splitcnt,tid,bufnames[tid]);
	  fflush(stdout);
	}
      }
      CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *matchtable2 = NULL;
      if(MULTIMATCH_TWOTHREADS)
	matchtable2 = new CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE>(tid, HashTables[tid], numrmaps, rmap, numqmaps, qmap, prefix, matchtable);

      /* pre-allocate memory used by matchtable->hash_query() */
      int Nmax1 = MULTIMATCH_TWOTHREADS ? ((Nmax + 15)/16) * 16 : Nmax;// Nmax rounded to nearest multiple of 16, to avoid false sharing of cache lines if MULTIMATCH_TWOTHREADS
      int *Xflag;
      POSIX_MEMALIGN(Xflag, 64, sizeof(int) * Nmax1 * (MULTIMATCH_TWOTHREADS ? 2 : 1));

      int Nmax2 = MULTIMATCH_TWOTHREADS ? ((Nmax + 2 + 15)/16) * 16 : Nmax + 2;// Nmax+2 rounded to nearest multiple of 16, to avoid false sharing of cache lines if MULTIMATCH_TWOTHREADS
      FLOAT *Xrev;
      POSIX_MEMALIGN(Xrev, 64, sizeof(FLOAT) * Nmax2 * (MULTIMATCH_TWOTHREADS ? 3 : 2));
      FLOAT *Xscale = &Xrev[Nmax2];

      if(QUICK){
	if(VERB && numqmaps > QUICK){
	  printf("numqmaps=%d->%d\n",numqmaps,QUICK);
	  fflush(stdout);
	}
	numqmaps = min(QUICK,numqmaps);
      }

      int block = (pairwise && !PairMerge) ? max(1,min(numqmaps/(4*hnumthreads),16)) : 1;// WAS33 (pairwise==0) ? 1 : 16;
    
#pragma omp for nowait schedule(dynamic, block)
      for(int i = 0; i < numqmaps; i++){
	if(DEBUG/* HERE >=2 */) assert(qmap[i]->mapid == i);
	if(0 /* TRACE */){/* skip irrelevant maps */
	  if (qmap[i]->id != TraceID2)
	    continue;
	}
#if 0 
	if (qmap[i]->id != 3LL)
	  continue;
#endif
	if(TRACE && qmap[i]->id == TraceID2){
#pragma omp critical
	  {
	    int M = qmap[i]->numsite[hcolor];
	    FLOAT *X = qmap[i]->site[hcolor];

	    printf("i=%d:mapid2=%d(id=%lld),M=%d: Tracing Hashtable queries\n",i,qmap[i]->mapid,qmap[i]->id,M);
	    for(int I = 1; I <= M+1; I++)
	      printf("  X[%d]= %0.3f (delta= %0.3f)\n", I, X[I],X[I] - X[I-1]);
	    fflush(stdout);
	  }
	}
	if(VERB && !(i%(QUICK ? 10 : 10000)) && i > 0){
#pragma omp critical
	  {
	    printf("split=%d/%d: Queried %d/%d maps: Found %lld matches with score >= %d : CPU time=%0.6f, wall time=%0.6f secs\n",
		   split,splitcnt, i, numqmaps, TotalMatchCnt + matchtable->totalMcnt, HashScore, mtime(),wtime());
	    if(VERB && ((VERB>=2 && numqmaps > 500000) /* || STAT*/)){
	      printf("\t");
	      matchtable->Hit_memPrint(1);
	      printf("\t");
	      matchtable->Match_memPrint(1);
	    }
	    fflush(stdout);
	  }
	}
	if(VERB>=2){
#pragma omp critical
	  {
	    printf("tid=%d:Calling hash_query:i=%d,id=%lld,sites=%d,len=%0.3f\n",
		   tid,i,qmap[i]->id,qmap[i]->numsite[hcolor],qmap[i]->site[hcolor][qmap[i]->numsite[hcolor]+1]);
	    //	  printf("Sleeping for 60 seconds\n");
	    fflush(stdout);
	    //	  sleep(60);
	  }
	}
	if(DEBUG) assert(qmap[i]->mapid == i);
	matchtable->hash_query(qmap[i],tid, Xflag, Xrev, Xscale, matchtable2);

	if(VERB>=1+RELEASE && numqmaps <= 30){
#pragma omp critical
	  {
	    printf("tid=%d:Completed hash_query:split=%d/%d:i=%d,id=%lld,sites=%d,len=%0.3f,MGHashGC=%d, Found %lld matches: wtime=%0.6f\n",
		   tid,split,splitcnt,i,qmap[i]->id,qmap[i]->numsite[hcolor],qmap[i]->site[hcolor][qmap[i]->numsite[hcolor]+1],matchtable->MTHashGC,TotalMatchCnt + matchtable->totalMcnt, wtime());
	    if(STAT){
#define MT matchtable
#define MT2 matchtable2
	      if(!MT2)
		printf("\t HashTable:probes=%lld,%lld,hits=%lld,inline col=%lld(%0.3f/hit),col=%lld(%0.3f/hit),matches=%lld(%0.3f/hit),total entries=%lld(%0.3f/match)\n",
		       MT->HTindcnt, MT->HTqindcnt, MT->HTfindcnt,MT->HTqcolcnt-MT->HTmatchcnt,((double)(MT->HTqcolcnt-MT->HTmatchcnt))/MT->HTfindcnt,
		       MT->HTcolcnt,(double)(MT->HTcolcnt)/MT->HTfindcnt, MT->HTmatchcnt, ((double)MT->HTmatchcnt)/MT->HTfindcnt, MT->HTentrycnt, ((double)MT->HTentrycnt)/MT->HTmatchcnt);
	      else
		printf("\t HashTable:probes=%lld,%lld,hits=%lld,inline col=%lld(%0.3f/hit),col=%lld(%0.3f/hit),matches=%lld(%0.3f/hit),total entries=%lld(%0.3f/match)\n",
		       MT->HTindcnt+MT2->HTindcnt, MT->HTqindcnt+MT2->HTqindcnt, MT->HTfindcnt + MT2->HTfindcnt, MT->HTqcolcnt+MT2->HTqcolcnt-MT->HTmatchcnt-MT2->HTmatchcnt,
		       ((double)(MT->HTqcolcnt+MT2->HTqcolcnt-MT->HTmatchcnt-MT2->HTmatchcnt))/max(1,MT->HTfindcnt+MT2->HTfindcnt),
		       MT->HTcolcnt+MT2->HTcolcnt,(double)(MT->HTcolcnt+MT2->HTcolcnt)/max(1,MT->HTfindcnt+MT2->HTfindcnt), MT->HTmatchcnt+MT2->HTmatchcnt, 
		       ((double)(MT->HTmatchcnt+MT2->HTmatchcnt))/max(1,MT->HTfindcnt+MT2->HTfindcnt), MT->HTentrycnt+MT2->HTentrycnt, ((double)(MT->HTentrycnt+MT2->HTentrycnt))/max(1,MT->HTmatchcnt+MT2->HTmatchcnt));
#undef MT
#undef MT2
	    }
	    fflush(stdout);
	  }
	}

      }// omp for nowait
    
      free(Xflag);
      free(Xrev);
      delete matchtable;
      if(MULTIMATCH_TWOTHREADS)
	delete matchtable2;
    } // parallel

#ifdef _OPENMP
    omp_set_nested(0);
#endif

    if(VERB){
      double mtime2 = mtime();
      double wtime2 = wtime();
      printf("split=%d/%d: Queried %d maps (against %d in hashtable): Found %lld matches with score >= %d : CPU time=%0.6f(delta=%0.6f),wall time=%0.6f secs(delta=%0.6f)\n", 
	     split,splitcnt, numqmaps, numrmaps, TotalMatchCnt, HashScore, mtime2, (mtime2-mtime1), wtime2, wtime2-wtime1);
      if(MemHWM > 0){
	printf("\t MatchTable memory HWM=%d+%lld(%llu Kbytes per thread) for mapid=%d, BestHWM=%d,OutHWM=%lu(%lu Kbytes per thread), MatchIndex=%d(%lu Kbytes per thread)\n",
	       MemHWM, BufHWM, (MemHWM*sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) + BufHWM*sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>))/1024, HWMid, 
	       BestHWM,OutHWM,(BestHWM*sizeof(int) + numrmaps*sizeof(CBestMatch<OFFSET_T,RANGE,HASH_SCALE>) + OutHWM*sizeof(CHashMatch))/1024,
	       (1U << MHASH_BITS),(1U << MHASH_BITS)*sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>)/1024);
	if(MHASH_BITS < 31 && MemHWM >= MASK(MHASH_BITS+1)){
	  int mhash_bits = MHASH_BITS;
	  while(mhash_bits < 31 && MemHWM >= MASK(mhash_bits + 1))
	    mhash_bits++;
	  printf("Use -MHash_Bits %d for better performance and reduced memory usage\n",mhash_bits);
	  fflush(stdout);
	}
      }
      if(HitHWM > 0){
	printf("\t HitTable memory HWM=%d(%lu Kbytes per thread), HitIndex=%d(%lu Kbytes per thread)\n",HitHWM,HitHWM*sizeof(CHitEntry<OFFSET1_T>)/1024,(1<<HHASH_BITS),(1<<HHASH_BITS)*sizeof(int)/1024);
	if(HitHWM >= (1U << (HHASH_BITS+1))){
	  int hhash_bits = HHASH_BITS+1;
	  while(hhash_bits < 30 && HitHWM >= (1U << hhash_bits))
	    hhash_bits++;
	  printf("Use -HHash_Bits %d for better performance and reduced memory usage\n",hhash_bits-1);
	  fflush(stdout);
	}
      }
      fflush(stdout);
    }

    cumTotalMatchCnt += TotalMatchCnt;

    if(hashtable->MatchMemReserve){
      delete hashtable->MatchMemReserve;
      hashtable->MatchMemReserve = 0;
    }

    PrintStatistics();

    /* Free up remaining hashtable memory */
    if(VERB){
      printf("\t Freeing hashtable memory:wall time=%0.6f secs\n",wtime());
      if(VERB/* HERE HERE >=2 */){
	static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
	getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
	printf("\t VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb\n",  VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9);
	fflush(stdout);
      }
      fflush(stdout);
    }
    delete hashtable;
    delete [] HashTables;

    if(usecolor && hcolor > 0)
      colors = 2;

    deltaXfreeBlock(rmap,roffset,roffset+numrmaps-1);
    if(qmap != rmap)
      deltaXfreeBlock(qmap,0,numqmaps-1);

    if(usecolor && hcolor > 0)
      colors = 1;

    if(MatrixFP)
      (void)FILEclose(MatrixFP);

    if(!(HASH_STREAM>=3 && HashGen && hash_filename && !HashTxt)){
      if(splitcnt > 1){
	printf("hashtable output not supported with splitcnt=%d\n",splitcnt);
	fflush(stdout); exit(1);
      }
      /* combine results from all threads in tmp files bufnames[0..numthreads-1] using merge sort */
      /* If !pairwise : Also output subset of input rmap[] with hashtable entries in <prefix>_hash.bnx */
      if(!HashMatrix){
	if(VERB>=2){
	  printf("Calling hash_output(%s)\n",prefix);
	  fflush(stdout);
	}
	hash_output(prefix, bufnames, numbufs,TotalMatchCnt, pairwise, rmap, numrmaps, qmap, numqmaps);
      }

      hash_free();
    } /* otherwise merge sort will be performed as part of hash_open(),hash_read(),hash_close() during pairwise alignment */

    /* free memory (other than merge sort memory freed by hash_free()) */
    MaxFragLen = 0.0;
    delete [] IntLen; IntLen = NULL; 
    MaxIntLen = 0; 

    SizeLen = 0;
    delete [] SizeDelta; SizeDelta = NULL;
    delete [] &YY[roffset]; YY = NULL;
    delete [] &Ylen[roffset]; Ylen = NULL;
    if(HashGC || OFFSET_FILL){
      OFFSET_T  *mymaxOffset1 = (OFFSET_T *)maxOffset1;
      delete [] &mymaxOffset1[roffset];
      maxOffset1 = minOffset1R = NULL;
    }
    delete [] id2mapid; id2mapid = NULL;
    delete [] MapSave; MapSave = NULL;

    SF[hcolor] = origSF;
    SD[hcolor] = origSD;
    SR[hcolor] = origSR;
    SDmax = origSDmax;

    MergeSort = 0;
    if(VERB){
      printf("split=%d/%d:Completed Freeing hashtable memory:wall time=%0.6f secs\n",split,splitcnt,wtime());
      if(VERB/* HERE HERE >=2 */){
	static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
	getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
	printf("\t VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb\n",  VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9);
	fflush(stdout);
      }
      if(VERB>=2)
	dumpmemmap();
      fflush(stdout);
    }
  }// for(int split = 0; split < splitcnt; split++)

  if(splitcnt > 1){
    delete [] prefix;
    prefix = orig_prefix;
    numbufs = split_numbufs;
    bufnames = split_bufnames;

    TotalMatchCnt = cumTotalMatchCnt;
  }
}

void hash_generate(Cmap *rmap[], int numrmaps, Cmap *qmap[], int numqmaps,  int pairwise, char *prefix)
{
  hash_free();/* NEW102 : in case it was not called by refalign() */

  /* OFFSET_T depends on max number of labels in rmap[] or qmap[] */
  int Nmax = 0;
  for(int i = 0; i < numrmaps; i++){
    Cmap *pmap = rmap[i];
    Nmax = max(pmap->numsite[hcolor],Nmax);
  }
  for(int i = 0; i < numqmaps; i++){
    Cmap *pmap = qmap[i];
    Nmax = max(pmap->numsite[hcolor],Nmax);
  }
  
  /* OFFSET_T && OFFSET1_T depend on MaxRoffset and MaxQoffset */

  OffsetKBinv = 1.0/OffsetKB;

  int MaxErr = 1 + max(NumErrI,NumErrQ) + NumResQ;

  /* compute MaxFragLen, MinOffset,MaxOffset */
  int MaxRoffset = 0, MaxQoffset = 0;
  MaxFragLen = 0.0;
  TotWin = 0;

  double qMaxScale = (hashScaleDelta && ScaleDelta) ? 1.0 + ScaleDelta * ScaleDeltaSize : 1.0;
  double rMaxScale = (hashScaleDelta && ScaleDelta && qmap == rmap) ? qMaxScale : 1.0;

  double totrlen = 0.0;

  for(int i = 0; i < numrmaps; i++){
    Cmap *pmap = rmap[i];
    int N = pmap->numsite[hcolor];
    if(N <= HashWin)
      continue;
    TotWin += (N-HashWin) + (NumErrI>=1 ? (N-HashWin-1)*HashWin : 0) + (NumErrI>=2 ? (N-HashWin-2)*HashWin*(HashWin-1)/2 : 0);
    if(N > MaxSite)
      MaxSite = N;

    FLOAT *X = pmap->site[hcolor];
    if(DEBUG>=2){
      assert(X[0] == 0.0);
      for(int I = 1; I <= N+1; I++){
	if(!(X[I] >= X[I-1])){
	  printf("rmap[%d]:I=%u,N=%u,X[I-1,I]=%0.4f,%0.4f,X[N]=%0.4f,X[N+1]=%0.4f\n",
		 i,I,N,X[I-1],X[I],X[N],X[N+1]);
	  fflush(stdout);
	  assert(X[I] >= X[I-1]);
	}
      }
    }
    totrlen += X[N+1];

    for(int k = N-MaxErr; k > 0; k--)
      if((X[k+MaxErr] - X[k]) * rMaxScale > MaxFragLen)
	MaxFragLen = (X[k+MaxErr] - X[k]) * rMaxScale;

    long long maxoffset = floor(X[N+1] * OffsetKBinv * rMaxScale + 0.5);
    if(maxoffset > MASK(31)){
      printf("hashtable cannot support maps larger than %0.3f kb : inserted map %d(id=%lld) has size %0.3f kb)\n",
	     MASK(31) * OffsetKB, i, pmap->id, X[N+1]);
      fflush(stdout); exit(1);
    }

    if(maxoffset > MaxRoffset)
      MaxRoffset = maxoffset;
  }
  if(HashMaxHits > 0 && HashHitsCov > 0 && HashHitsGenSiz > 0){
    int origHashMaxHits = HashMaxHits;
    HashMaxHits = max(1.0, floor(HashHitsCov * totrlen / (HashHitsGenSiz * 1000.0) + 0.5));
    if(VERB){
      printf("Adjusting -HashMaxHits from %d to %d (HashHitsCov= %d, HashHitsGenSiz= %d Mb, totrlen= %0.3f Mb)\n",origHashMaxHits,HashMaxHits,HashHitsCov,HashHitsGenSiz,totrlen*0.001);
      fflush(stdout);
    }
  }

  if(qmap != rmap){
    for(int i = 0; i < numqmaps; i++){
      Cmap *pmap = qmap[i];
      int N = pmap->numsite[hcolor];
      FLOAT *X = pmap->site[hcolor];
      if(DEBUG>=2){
	assert(X[0] == 0.0);
	for(int I = 1; I <= N+1; I++){
	  if(!(X[I] >= X[I-1])){
	    printf("qmap[%d]:I=%u,N=%u,X[I-1,I]=%0.4f,%0.4f,X[N]=%0.4f,X[N+1]=%0.4f\n",
		   i,I,N,X[I-1],X[I],X[N],X[N+1]);
	    fflush(stdout);
	    assert(X[I] >= X[I-1]);
	  }
	}
      }
      if(N <= HashWin)
	continue;
      for(int k = N-MaxErr; k > 0; k--)
	if((X[k+MaxErr] - X[k]) * qMaxScale > MaxFragLen)
	  MaxFragLen = (X[k+MaxErr] - X[k]) * qMaxScale;
      long long maxoffset = floor(X[N+1] * OffsetKBinv * qMaxScale + 0.5);
      if(maxoffset > MASK(31)){
	printf("hashtable cannot support maps larger than %0.3f kb : probe map %d(id=%lld) has size %0.3f kb)\n",
	       MASK(31) * OffsetKB, i, pmap->id, X[N+1]);
	fflush(stdout); exit(1);
      }
      if(maxoffset > MaxQoffset)
	MaxQoffset = maxoffset;
    }
  } else
    MaxQoffset = MaxRoffset;

  /* MatchTable uses offset = offset2 - offset1 which is in the range (0-MaxRoffset)-OFFSET_SPREAD ... (MaxQoffset-0)+OFFSET_SPREAD.
     HitTable uses offset1 which is in the range 0 ...  MaxRoffset+OFFSET_SPREAD */
  MinOffset = -MaxRoffset - OFFSET_SPREAD;
  MaxOffset = MaxQoffset + OFFSET_SPREAD;
  MaxOffset = max(MaxOffset,MaxRoffset) + OFFSET_SPREAD;

  int OFFSET_T_size = (Nmax > (int)MASK(15) || MaxOffset > (int)MASK(15) || MinOffset < -(1<<15)) ? 4 : 2;
  int OFFSET1_T_size = (MaxRoffset + OFFSET_SPREAD > (int)MASK(16)) ? 4 : 2;
  int HASHWIN_val = HashWin;
  int RANGE_val = (ShiftOffset1 <= 30) ? 1 : 0;
  int HASH_SCALE_val = (hashScaleDelta && NumScaleFactor > 1) ? 1 : 0;
  if(DEBUG && HASH_SCALE_val) assert(ScaleDelta);

  if(VERB/* HERE >=2 */){
    printf("Using hashtable template: sizeof(OFFSET_T)=%d,sizeof(OFFSET1_T)=%d,HASHWIN=%d,RANGE=%d,HASH_SCALE=%d (Nmax=%d,MaxOffset=%d,MinOffset=%d,MaxRoffset=%d,MaxQoffset=%d,Scale=%d,%0.3f)\n",
	   OFFSET_T_size,OFFSET1_T_size,HASHWIN_val,RANGE_val,HASH_SCALE_val,Nmax,MaxOffset,MinOffset,MaxRoffset,MaxQoffset,NumScaleFactor,NumScaleFactor > 1 ? ScaleFactor[1]-1.0 : 0.0);
    printf("OffsetKBinv= %0.8f rMaxScale= %0.8f, qMaxScale= %0.8f\n",OffsetKBinv,rMaxScale,qMaxScale);
    fflush(stdout);
  }
  if(HashWin > ((OFFSET1_T_size == 2) ? 6 : 5)){
    printf("HashWin = %d is too large : largest possible value is MAX_WIN(OFFSET1_T)= %d (sizeof(OFFSET1_T) = %d)\n", HashWin, (OFFSET1_T_size == 2) ? 6 : 5, OFFSET1_T_size);
    fflush(stdout);exit(1);
  }
  if(HashWin < 5){
    printf("HashWin = %d : values smaller than 5 are currently not supported with hashtable templates\n", HashWin);
    fflush(stdout);exit(1);
  }


  if(OFFSET_T_size == 2 && OFFSET1_T_size == 2 && HASHWIN_val==5 && RANGE_val == 0 && HASH_SCALE_val == 0)
    hash_generateT<short, unsigned short, 5, 0, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 2 && HASHWIN_val==5 && RANGE_val == 0 && HASH_SCALE_val == 1)
    hash_generateT<short, unsigned short, 5, 0, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
#if USE_MIC==0
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 2 && HASHWIN_val==5 && RANGE_val == 1 && HASH_SCALE_val == 0)
    hash_generateT<short, unsigned short, 5, 1, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 2 && HASHWIN_val==5 && RANGE_val == 1 && HASH_SCALE_val == 1)
    hash_generateT<short, unsigned short, 5, 1, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 2 && HASHWIN_val==6 && RANGE_val == 0 && HASH_SCALE_val == 0)
    hash_generateT<short, unsigned short, 6, 0, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 2 && HASHWIN_val==6 && RANGE_val == 0 && HASH_SCALE_val == 1)
    hash_generateT<short, unsigned short, 6, 0, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 2 && HASHWIN_val==6 && RANGE_val == 1 && HASH_SCALE_val == 0)
    hash_generateT<short, unsigned short, 6, 1, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 2 && HASHWIN_val==6 && RANGE_val == 1 && HASH_SCALE_val == 1)
    hash_generateT<short, unsigned short, 6, 1, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
#if 0 // OFFSET_T_size < OFFSET1_T_size should never be needed
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 4 && HASHWIN_val==5 && RANGE_val == 0 && HASH_SCALE_val == 0)
    hash_generateT<short, int, 5, 0, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 4 && HASHWIN_val==5 && RANGE_val == 0 && HASH_SCALE_val == 1)
    hash_generateT<short, int, 5, 0, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 4 && HASHWIN_val==5 && RANGE_val == 1 && HASH_SCALE_val == 0)
    hash_generateT<short, int, 5, 1, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 4 && HASHWIN_val==5 && RANGE_val == 1 && HASH_SCALE_val == 1)
    hash_generateT<short, int, 5, 1, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 4 && HASHWIN_val==6 && RANGE_val == 0 && HASH_SCALE_val == 0)
    hash_generateT<short, int, 6, 0, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 4 && HASHWIN_val==6 && RANGE_val == 0 && HASH_SCALE_val == 1)
    hash_generateT<short, int, 6, 0, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 4 && HASHWIN_val==6 && RANGE_val == 1 && HASH_SCALE_val == 0)
    hash_generateT<short, int, 6, 1, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 2 && OFFSET1_T_size == 4 && HASHWIN_val==6 && RANGE_val == 1 && HASH_SCALE_val == 1)
    hash_generateT<short, int, 6, 1, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
#endif
#endif
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 2 && HASHWIN_val==5 && RANGE_val == 0 && HASH_SCALE_val == 0)
    hash_generateT<int, unsigned short, 5, 0, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 2 && HASHWIN_val==5 && RANGE_val == 0 && HASH_SCALE_val == 1)
    hash_generateT<int, unsigned short, 5, 0, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
#if USE_MIC==0
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 2 && HASHWIN_val==5 && RANGE_val == 1 && HASH_SCALE_val == 0)
    hash_generateT<int, unsigned short, 5, 1, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 2 && HASHWIN_val==5 && RANGE_val == 1 && HASH_SCALE_val == 1)
    hash_generateT<int, unsigned short, 5, 1, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 2 && HASHWIN_val==6 && RANGE_val == 0 && HASH_SCALE_val == 0)
    hash_generateT<int, unsigned short, 6, 0, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 2 && HASHWIN_val==6 && RANGE_val == 0 && HASH_SCALE_val == 1)
    hash_generateT<int, unsigned short, 6, 0, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 2 && HASHWIN_val==6 && RANGE_val == 1 && HASH_SCALE_val == 0)
    hash_generateT<int, unsigned short, 6, 1, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 2 && HASHWIN_val==6 && RANGE_val == 1 && HASH_SCALE_val == 1)
    hash_generateT<int, unsigned short, 6, 1, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 4 && HASHWIN_val==5 && RANGE_val == 0 && HASH_SCALE_val == 0)
    hash_generateT<int, int, 5, 0, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 4 && HASHWIN_val==5 && RANGE_val == 0 && HASH_SCALE_val == 1)
    hash_generateT<int, int, 5, 0, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 4 && HASHWIN_val==5 && RANGE_val == 1 && HASH_SCALE_val == 0)
    hash_generateT<int, int, 5, 1, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 4 && HASHWIN_val==5 && RANGE_val == 1 && HASH_SCALE_val == 1)
    hash_generateT<int, int, 5, 1, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
#if 0 // OFFSET1_T == int and HASHWIN==6 not supported
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 4 && HASHWIN_val==6 && RANGE_val == 0 && HASH_SCALE_val == 0)
    hash_generateT<int, int, 6, 0, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 4 && HASHWIN_val==6 && RANGE_val == 0 && HASH_SCALE_val == 1)
    hash_generateT<int, int, 6, 0, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 4 && HASHWIN_val==6 && RANGE_val == 1 && HASH_SCALE_val == 0)
    hash_generateT<int, int, 6, 1, 0> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
  else if(OFFSET_T_size == 4 && OFFSET1_T_size == 4 && HASHWIN_val==6 && RANGE_val == 1 && HASH_SCALE_val == 1)
    hash_generateT<int, int, 6, 1, 1> (rmap,numrmaps,qmap,numqmaps,pairwise,prefix);    
#endif
#endif
  else {
    printf("Unhandled template case: OFFSET_T_size=%d,OFFSET1_T+size=%d,HASHWIN=%d,RANGE=%d,HASH_SCALE=%d\n",OFFSET_T_size,OFFSET1_T_size,HASHWIN_val,RANGE_val,HASH_SCALE_val);
    fflush(stdout);exit(1);
  }

  if(QUICKSTAT){
    printf("hashfindV(): calls = %lu, cache lines = %lu (%0.3f lines/call), max lines= %d\n",hashfindV_calls,hashfindV_lines,((double)hashfindV_lines)/hashfindV_calls,hashfindV_max);
    fflush(stdout);
  }

  //  exit(1);
}

