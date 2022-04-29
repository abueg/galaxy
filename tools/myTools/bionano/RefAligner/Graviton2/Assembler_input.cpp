#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include <math.h>
#include <omp.h>

#include "globals.h"
#include "Assembler_parameters.h"

#if CALIGN_SMALL
#include "CalignA.h"
#define Calign CalignA
#define Calign_block CalignA_block

#define alignment alignmentA
#define align_block alignA_block
#define num_align_blocks num_alignA_blocks

#define align_blockalloc alignA_blockalloc
#else
#include "Calign.h"
#endif

#if defined(_MSC_VER)
#define strtoll _strtoi64
#endif

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Assembler_input.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

int AlignmentID_ordered = 1;

double origMinLen = 0.0;
double origLogPvThreshold = 0.0;

static int mapsdeleted = 0;

static  size_t mapfiltered = 0;
static  size_t filtered = 0;/* number of filtered alignments */
static  size_t threshcnt = 0;
static  size_t repeatcnt = 0;/* number of alignments entirely within repeat regions */
static  size_t chimcnt = 0;/* total alignments removed as chimeric */
static  size_t chimcntTP = 0;/* true positives mistaken for chimeric alignment */
static  size_t chimcntFP = 0;/* false positive mistaken for chimeric alignment */
static size_t cntTP = 0;/* confirmed true positive alignment (simulated molecules that overlap and have correct orientation : amount of overlap is not necessarily correct) */
static double PChimThreshold = 1000.0;  

static char buf[LINESIZ];

//#undef DEBUG
//#define DEBUG 2

double AlignedLength(Calign *p, Cmap *map[])
{
  //  if(AlignedLengthThreshold <= 0.0)
  //    return 1.0;
  int U = p->numpairs;
  Cmap *Ymap = map[p->mapid1];
  FLOAT *Y = Ymap->site[0];
  int N = Ymap->numsite[0];
  Cmap *Xmap = map[p->mapid2];
  FLOAT *X = Xmap->site[0];
  int M = Xmap->numsite[0];
  FLOAT len = Y[p->sites1[U-1]] - Y[p->sites1[0]];
  if(p->orientation)
    len += X[M+1 - p->sites2[0]] - X[M+1 - p->sites2[U-1]];
  else
    len += X[p->sites2[U-1]] - X[p->sites2[0]];
  if(DEBUG/* HERE >=2 */ && !(U>=0 && 0 < p->sites1[0] && p->sites1[0] <= p->sites1[U-1] && p->sites1[U-1] <= N && 
			      0 < p->sites2[0] && p->sites2[0] <= p->sites2[U-1] && p->sites2[U-1] <= M && ((U > 1) ? (len > 0.0) : (len >= 0.0)))){
    #pragma omp critical
    {
      printf("mapid1=%d(id=%lld),mapid2=%d(id=%lld),or=%d:U=%d,M=%d,N=%d,sites1[0,U-1]=%d,%d(Y=%0.3f,%0.3f),sites2[0,U-1]=%d,%d(X=%0.3f,%0.3f),len=%0.3f\n",
	     p->mapid1,map[p->mapid1]->id,p->mapid2,map[p->mapid2]->id,p->orientation,U,M,Ymap->numsite[0],p->sites1[0],p->sites1[U-1],Y[p->sites1[0]],Y[p->sites1[U-1]],
	     p->sites2[0],p->sites2[U-1],p->orientation ? X[M+1] - X[p->sites2[0]] : X[p->sites2[0]], p->orientation ? X[M+1] - X[p->sites2[U-1]] : X[p->sites2[U-1]], len);
      fflush(stdout);
      assert(U>=0);
      assert(0 < p->sites1[0] && p->sites1[0] <= p->sites1[U-1] && p->sites1[U-1] <= N);
      assert(0 < p->sites2[0] && p->sites2[0] <= p->sites2[U-1] && p->sites2[U-1] <= M);
      assert((U > 1) ? (len > 0.0) : (len >= 0.0));
    }
  }
  return  len * 0.5;
}


Cid2mapid *id2mapid = 0;

class Cid2mapidHash {/* NOTE : struct is 64 bytes long to fit in 1 or 2 cache lines */
public:
  long long id[5];
  int mapid[5];
  int cnt;

  inline void init(int C,long long Id, int Mapid){
    if(DEBUG>=2) assert(0 <= C && C < 5);
    id[C] = Id;
    mapid[C] = Mapid;
    cnt = C+1;
  };
};

static Cid2mapidHash *HashTable = 0;
static int *HashIndex = 0;
static int Next = 5, HashMax = 0;/* Next is last used entry in HashTable[6..HashMax-1] */

#define H_BITS 21

static inline int HashKey(long long id)
{
  return (id ^ (id >> H_BITS) /* ^ (id >> (2*H_BITS)) */) & MASK(H_BITS);
}

static inline int HashNext()
{
  if(++Next >= HashMax){
    int origHashMax = HashMax;
    HashMax = (HashMax * 3) / 2;
    Cid2mapidHash *newHashTable = (Cid2mapidHash *)malloc(HashMax * sizeof(Cid2mapidHash));
    if(!newHashTable){
      printf("HashNext:malloc(%llu) failed\n", (unsigned long long)(HashMax * sizeof(Cid2mapidHash)));
      exit(1);
    }
    memcpy(&newHashTable[6],&HashTable[6], (origHashMax-6) * sizeof(Cid2mapidHash));
    free(HashTable);
    HashTable = newHashTable;
  }
  return Next;
}

static void HashInit()
{
  if(HashIndex)
    free(HashIndex);
  if(!(HashIndex = (int *)calloc(1 << H_BITS, sizeof(int)))){
    printf("HashInit:calloc(%llu) failed\n",(unsigned long long)((1 << H_BITS)*sizeof(int)));
    exit(1);
  }

  if(HashTable)
    free(HashTable);
  HashMax = (1 << H_BITS);
  if(!(HashTable = (Cid2mapidHash *)malloc(HashMax*sizeof(Cid2mapidHash)))){
    printf("HashInit:malloc(%llu) failed\n", (unsigned long long)(HashMax * sizeof(Cid2mapidHash)));
    exit(1);
  }
  Next = 5;
}

/* NOTE : it is illegal to insert the same id twice */
static void HashInsert(long long id, int mapid)
{
  int key = HashKey(id);
  int index = HashIndex[key];
  if(!index){
    HashIndex[key] = index = HashNext();
    if(DEBUG>=2) assert(6 <= index && index <= Next);
    Cid2mapidHash *p = &HashTable[index];
    p->init(0,id,mapid);
  } else {
    if(DEBUG>=2) assert(6 <= index && index <= Next);
    Cid2mapidHash *p = &HashTable[index];
    int cnt = p->cnt;
    if(DEBUG>=2)/* check that id is not already in HashTable */
      for(int i = min(5,cnt); --i >= 0;)
	assert(p->id[i] != id);
    while(cnt > 5){/* hash line overflowed */
      if(DEBUG>=2) assert(cnt <= Next);
      p = &HashTable[index = cnt];
      cnt = p->cnt;
      if(DEBUG>=2)/* check that id is not already in HashTable */
	for(int i = min(5,cnt); --i >= 0;)
	  assert(p->id[i] != id);
    }
    if(DEBUG>=2) assert(p == &HashTable[index] && cnt == p->cnt);

    if(cnt >= 5){/* full */
      int nindex = HashNext();
      HashTable[index].cnt = nindex;// HashTable may have been reallocated by HashNext(), invalidating previous value of p
      p = &HashTable[nindex];
      p->init(0,id,mapid);
    } else
      p->init(cnt,id,mapid);
  }
}

/* returns mapid or -1 if not found */
int HashFindID(long long id)
{
  int key = HashKey(id);
  int index = HashIndex[key];
  if(!index)
    return -1;
  if(DEBUG>=2) assert(6 <= index && index <= Next);
  Cid2mapidHash *p = &HashTable[index];
  int cnt = p->cnt;
  for(int i = min(5,cnt); --i >= 0;)
    if(p->id[i] == id)
      return p->mapid[i];
  while(cnt > 5){/* overflowed */
    if(DEBUG>=2) assert(cnt <= Next);
    p = &HashTable[cnt];
    cnt = p->cnt;
    for(int i = min(5,cnt); --i >= 0;)
      if(p->id[i] == id)
	return p->mapid[i];
  }
  return -1;
}

/* read in binary alignment data from a single file that has already been read in to bufferbin[0..lenbin-1].

   The file contained "num" alignments with byte offset locations of offset_vec[0..num-1] into bufferbin[].

   Uses global array id2mapid[0..nummaps-1].
 */
static void input_align_binary(size_t &numaligns, size_t num, size_t *offset_vec, char *bufferbin, size_t lenbin, char *filenamebin, int fileno)
{
  int numthreads = 1;
#ifdef _OPENMP
  numthreads = omp_get_max_threads();
  numthreads = min(numthreads,MaxThreads);
  int nthreads = num/512;
  numthreads = max(1,min(numthreads, nthreads));
#endif

  long long numLL = num;

  size_t offset = 0;/* current offset in bufferbin[] */

  int *Yid = new int[num*2];
  int *Xid = &Yid[num];

  /* get the 1st id array and translate it into mapids in Yid[] */
  long long *id_vecY = (long long*)(&bufferbin[offset]);
  offset += num * sizeof(long long);
  if(DEBUG) assert(offset <= lenbin);

  #pragma omp parallel for num_threads(numthreads) schedule(static,512)
  for(long long i = 0; i < numLL; i++){
    Yid[i] = HashFindID(id_vecY[i]);
    if(VERB>=2){
      printf("i=%lld/%lld: id_vecY[i]= %lld, HashFindID= Yid[i]= %d\n", i, numLL, id_vecY[i], Yid[i]);
      fflush(stdout);
    }
  }

  if(VERB>=2){
    printf("Translated Yid[]:time=%0.6f(wall time=%0.6f)\n",mtime(),wtime());
    fflush(stdout);
  }

  /* get the 2nd id array and translate it into mapids in Xid[] */
  long long *id_vecX = (long long*)(&bufferbin[offset]);
  offset += num * sizeof(long long);
  if(DEBUG) assert(offset <= lenbin);

  #pragma omp parallel for num_threads(numthreads) schedule(static,512)
  for(long long i = 0; i < numLL; i++){
    Xid[i] = HashFindID(id_vecX[i]);
    if(VERB>=2){
      printf("i=%lld/%lld: id_vecX[i]= %lld, HashFindID= Xid[i]= %d (id_vecY[i]= %lld, Yid[i]= %d)\n", i, numLL, id_vecX[i], Xid[i], id_vecY[i], Yid[i]);
      fflush(stdout);
    }
  }

  if(VERB>=2){
    printf("Translated Xid[]:time=%0.6f(wall time=%0.6f)\n",mtime(),wtime());
    fflush(stdout);
  }

  /* get the numpairs array */
  int *numpairs = (int *)(&bufferbin[offset]);
  offset += num * sizeof(int);
  if(DEBUG) assert(offset <= lenbin);
  
  /* expand alignment array */
  maxalignallocNULL(numaligns + num, alignment, numaligns, maxaligns);

  if(VERB>=2){
    printf("called maxalignallocNULL(%lu):time=%0.6f(wall time=%0.6f)\n",numaligns + num, mtime(),wtime());
    fflush(stdout);
  }

  /* compute sum of numpairs[0..num-1] */
  size_t *numpairs_cum = new size_t[num+1];
  size_t numpairs_sum = 0;

  for(size_t i = 0; i < num; i++){
    numpairs_cum[i] = numpairs_sum;
    numpairs_sum += numpairs[i];
  }
  numpairs_cum[num] = numpairs_sum;

  if(VERB>=2){
    printf("Computed prefix sum of numpairs[]:time=%0.6f(wall time=%0.6f)\n",mtime(),wtime());
    fflush(stdout);
  }

  int block_id;
  Calign_block *block = align_blockalloc(num, numpairs_sum, 0ul, AlignScore, 0, block_id);

  size_t clen = (colors==1) ? 1 : colors + 1;

  if(DEBUG) assert(block->numaligns == num);

  #pragma omp parallel for num_threads(numthreads) schedule(static,512)
  for(long long i = 0; i < numLL; i++){
    if(DEBUG>=2) assert((i+1)*clen <= block->alignsiz);
    Calign *p = alignment[numaligns+i] = &block->block[i*clen];
    size_t nump = numpairs[i];
    size_t offset = 2 * numpairs_cum[i];
    p->block_id = block_id;
    p->allocated_size = nump;
    p->sites1 = &block->intblock[offset];
    p->sites2 = &block->intblock[offset+nump];
    if(DEBUG>=2) assert((char *)(&p->sites2[nump]) <= ((char *)block->block) + block->memsiz);
#if !CALIGN_SMALL
    p->sitesK1 = NULL;
#endif
    if(AlignScore){
      p->iscore = &block->scoreblock[offset];
      p->outscore = &block->scoreblock[offset+nump];
      if(DEBUG>=2 && !((char *)(&p->outscore[nump]) <= ((char *)block->block) + block->memsiz)){
	printf("i=%lld,num=%lu,nump=%lu,offset=%lu(numpairs_cum[i]=%lu):block->scoreblock=%p,block->block=%p,block->intblock=%p,block->memsiz=%lu\n",
	       i,num,nump,offset,numpairs_cum[i],(void *)block->scoreblock,(void *)block->block,(void *)block->intblock,block->memsiz);
	fflush(stdout);
	assert((char *)(&p->outscore[nump]) <= ((char *)block->block) + block->memsiz);
      }
    } else
      p->iscore = p->outscore = NULL;
  }

  block->nextint = 2 * numpairs_cum[num];
  if(DEBUG && !(block->nextint <= block->intsiz)){
    printf("num=%lu,numpairs_cum[num]=%lu, block->nextint=%lu, block->intsiz=%lu\n",num,numpairs_cum[num], block->nextint, block->intsiz);
    fflush(stdout);
    assert(block->nextint <= block->intsiz);
  }
  if(AlignScore){
    block->nextscore = 2 * numpairs_cum[num];
    if(DEBUG && !(block->nextscore <= block->scoresiz)){
      printf("num=%lu,numpairs_cum[num]=%lu, block->nextscore=%lu, block->scoresiz=%lu\n", num, numpairs_cum[num], block->nextscore, block->scoresiz);
      printf("    block->nextint=%lu, block->intsiz=%lu\n", block->nextint, block->intsiz);
      fflush(stdout);
      assert(block->nextscore <= block->scoresiz);
    }
  }

  if(VERB>=2){
    printf("Pre-allocated %lu alignments:time=%0.6f(wall time=%0.6f)\n",num,mtime(),wtime());
    fflush(stdout);
  }

  #pragma omp parallel num_threads(numthreads)
  {
  long long mymapfiltered = 0;
  long long mythreshcnt = 0;
  long long myrepeatcnt = 0;
  long long mycntTP = 0;
  long long mychimcnt = 0;
  long long mychimcntTP = 0;
  long long mychimcntFP = 0;
  long long myfiltered = 0;

  #pragma omp for schedule(static,512)
  for(long long i = 0; i < numLL; i++){  /* read in the alignment information */
    Calign *palign = alignment[numaligns + i];
    palign->Lend = palign->Rend = -2;
    palign->fileid = fileno;// NEW10

    palign->mapid1 = Yid[i];
    palign->mapid2 = Xid[i];

    if(palign->mapid1 < 0 || palign->mapid2 < 0){/* maps have been filtered */
      if(VERB>=3 || (DEBUG>=2 && MinLen == origMinLen && LogPvThreshold == origLogPvThreshold)){
	#pragma omp critical
	{
	  long long idY = id_vecY[i], idX = id_vecX[i];
	  int indY = HashFindID(idY);
	  int indX = HashFindID(idX);
	  printf("alignment[%llu]:i=%lld/%lld:mapid1=%d(id=%lld,Yid[i]=%d,ind=%d),mapid2=%d(id_vecX[i]=%lld,Xid[i]=%d,ind=%d) : one of the maps has been filtered\n",
		 (long long unsigned)(numaligns+i),i,numLL,palign->mapid1,id_vecY[i],Yid[i],indY,palign->mapid2,id_vecX[i],Xid[i],indX);
	  fflush(stdout);
	  //	  if(DEBUG) assert(!(MinLen == origMinLen && LogPvThreshold == origLogPvThreshold));
	}
      }
      mymapfiltered++;
      alignment[numaligns + i] = NULL;      // defer removal of filtered alignment
      continue;
    }
    
    size_t noffset = offset_vec[i];
    if(DEBUG>=2 && !(noffset >= offset && noffset < lenbin)){
      #pragma omp critical
      {
        printf("i=%lld/%lu:noffset=%lu,offset=%lu,lenbin=%lu\n",
	     i,num,noffset,offset,lenbin);
	fflush(stdout);
	assert(noffset >= offset && noffset < lenbin);/* not valid for multithreading */
      }
    }
    
    memcpy(&palign->score, &bufferbin[noffset], sizeof(double));
    noffset += sizeof(double);
    if(DEBUG>=2) assert(noffset < lenbin);
    
    memcpy(&palign->logPV, &bufferbin[noffset], sizeof(double));
    noffset += sizeof(double);
    if(DEBUG>=2) assert(noffset < lenbin);    
    
    memcpy(&palign->trueoffset, &bufferbin[noffset], sizeof(double));
    noffset += sizeof(double);
    if(DEBUG>=2) assert(noffset < lenbin);

    memcpy(&palign->trueoverlap, &bufferbin[noffset], sizeof(double));
    noffset += sizeof(double);
    if(DEBUG>=2) assert(noffset < lenbin);

    memcpy(&palign->orientation, &bufferbin[noffset], sizeof(int));
    noffset += sizeof(int);
    if(DEBUG>=2) assert(noffset < lenbin);

    memcpy(&palign->scaleID, &bufferbin[noffset], sizeof(int));
    noffset += sizeof(int);
    if(DEBUG>=2) assert(noffset <= lenbin);

    palign->trueTP = (palign->trueoffset==0.0) ? -1 : (palign->trueoverlap > 0.0) ? 1 : 0;

    palign->numpairs = numpairs[i];

    /* allocate arrays palign->sites1[] etc */
    if(!MaxCoverageDump && hash_filename)/* skip reading actual alignment information */
      continue;

    if(DEBUG>=2) assert(&palign->sites2[0] == &palign->sites1[palign->numpairs]);
    memcpy(&palign->sites1[0], &bufferbin[noffset], sizeof(int) * palign->numpairs * 2);
    noffset += sizeof(int) * palign->numpairs * 2;
    if(DEBUG>=2) assert(noffset <= lenbin);

    /*    memcpy(&palign->sites2[0], &bufferbin[noffset], sizeof(int) * palign->numpairs);
    noffset += sizeof(int) * palign->numpairs;
    if(DEBUG>=2) assert(noffset <= lenbin);*/
    
    if(AlignScore){
      if(DEBUG>=2) assert(&palign->outscore[0] == &palign->iscore[palign->numpairs]);
      memcpy(&palign->iscore[0], &bufferbin[noffset], sizeof(FLOAT) * palign->numpairs * 2);
      noffset += sizeof(FLOAT) * palign->numpairs * 2;
      if(DEBUG>=2) assert(noffset <= lenbin);

      /*      memcpy(&palign->outscore[0], &bufferbin[noffset], sizeof(FLOAT) * palign->numpairs);
      noffset += sizeof(FLOAT) * palign->numpairs;
      if(DEBUG>=2) assert(noffset <= lenbin); */
    }

    if(DEBUG && i+1 < numLL) assert(noffset == offset_vec[i+1]);
    else if(DEBUG) assert(noffset <= lenbin);

    if(DEBUG>=2)
      offset = noffset;/* only needed for non-multithreaded case */

    /* check if alignment needs to be filtered out */

    /* filter early to reduce memory footprint */
    if(!(palign->score > ScoreThreshold && 
	 palign->logPV > LogPvThreshold &&
	 palign->numpairs >= AlignedSiteThreshold &&
	 AlignedLength(palign,Gmap) >= AlignedLengthThreshold)){
      mythreshcnt++;
      alignment[numaligns + i] = NULL; // defer removal of filtered alignment
      continue;
    }

    if(AlignScore && RepeatMaxShift > 0){/* filter output alignments entirely within a repeat region */
      bool repeat = palign->IsRepeatRegion();
      if(repeat){
	myrepeatcnt++;
	alignment[numaligns + i] = NULL;// defer removal of filtered alignment
	continue;
      }
    }

    if(DEBUG>=2) assert(palign->mapid1 < nummaps);
    if(DEBUG>=2) assert(palign->mapid2 < nummaps);

    /* filter out poor and chimeric alignments */
    //    char *name1 = Gmap[palign->mapid1]->name;
    //    char *name2 = Gmap[palign->mapid2]->name;

    register double cum = 0.0;
    register double hwm = 0.0;
    register int Khwm = -1;
    register double maxdrop = 0.0;
    register int start = -1, end = -1;
    if(AlignScore){
      for(register int k = 0; k <= palign->numpairs;k++){
	register double score = palign->iscore[k];
	cum += score;
	if(cum > hwm){
	  hwm = cum;
	  Khwm = k;
	} else if(hwm - cum > maxdrop){
	  maxdrop = hwm - cum;
	  start = Khwm;
	  end = k;
	}
      }
    }

    if(VERB>=3){
      printf("alignment[%llu]:i=%lld/%lld:mapid1=%d(id=%lld),mapid2=%d(id=%lld): trueoffset= %0.3f, trueoverlap= %0.3f, trueTP=%d (cntTP=%lld)\n",
	     (long long unsigned)(numaligns+i),i,numLL,palign->mapid1,id_vecY[i],palign->mapid2,id_vecX[i],palign->trueoffset, palign->trueoverlap, palign->trueTP, mycntTP + cntTP);
      fflush(stdout);
    }

    if(palign->trueTP >= 0){/* simulated molecules : check if this alignment has correct orientation and neither molecule is chimeric */
      long long id1 = Gmap[palign->mapid1]->id;
      long long id2 = Gmap[palign->mapid2]->id;
      if(id1 > (LOC_WIDTH*LEN_WIDTH*FLAG_WIDTH) && id2 > (LOC_WIDTH*LEN_WIDTH*FLAG_WIDTH)){// see test_generator.cpp
	if((id1 & 0x4) || (id2 & 0x4))/* chimeric molecules are assumed to never align correctly */
	  palign->trueTP = 0;
	else {
	  bool TrueFlip = (id1 & 0x1) != (id2 & 0x1);
	  if(TrueFlip != (palign->orientation))
	    palign->trueTP = 0;
	}
      } else
	palign->trueTP = -1;

      mycntTP += (palign->trueTP >= 1) ? 1 : 0;
    }

    /* check for chimeric alignments with leftmost/rightmost misaligned sites exceed ENDSITES */
    /* postpone this till later ? graph_trim() can do a better job analyzing each maps local coverage based on all its pairwise alignments
       and make a better recognition of both chimeric maps as well as maps with regions that do not align : such maps need to be trimmed off,
       rather than have their alignments deleted. */
    if((min(palign->sites1[0],palign->sites2[0]) > ENDSITES ||
	min(Gmap[palign->mapid1]->numsite[0] - palign->sites1[palign->numpairs-1],
	    Gmap[palign->mapid2]->numsite[0] - palign->sites2[palign->numpairs-1]) > ENDSITES) ||
       (AlignScore && maxdrop > PChimThreshold && end - start >= Pchim_minlen)){

      mychimcnt++;

      if(VERB>=3){
	#pragma omp critical
	{
	  printf("alignment[%llu]:i=%llu/%llu:mapid1=%d(id=%lld),mapid2=%d(id=%lld):maxdrop= %0.3f, start=%d, end=%d : skipping due to -PVchim %0.2e %d\n",
		 (long long unsigned)(numaligns+i),i,numLL,palign->mapid1,id_vecY[i],palign->mapid2,id_vecX[i],maxdrop,start,end,Pchim, Pchim_minlen);
	  fflush(stdout);
	}
      }

      /* check if we deleted a true positive */
      if(palign->trueTP >= 0){
	if(palign->trueTP > 0){
	  mychimcntTP++;
	  mycntTP--;
	} else
	  mychimcntFP++;
      }
      myfiltered++;
      alignment[numaligns + i] = NULL;// defer removal of filtered alignment
      continue;
    }

    if(MinOverlapRatio > 0.0){
      Cmap *Ymap = Gmap[palign->mapid1];
      Cmap *Xmap = Gmap[palign->mapid2];
      FLOAT *Y = Ymap->site[0];
      FLOAT *X = Xmap->site[0];
      int N = Ymap->numsite[0];
      int M = Xmap->numsite[0];
      int I = palign->sites1[0];
      int J = palign->sites2[0];
      int U = palign->numpairs;
      int RI = palign->sites1[U-1];
      int RJ = palign->sites2[U-1];
      double overlap = 0.5*(Y[RI]-Y[I]) + (palign->orientation ? (X[M+1-J] - X[M+1-RJ])*0.5 + min(X[M+1]-X[M+1-J],Y[I]) + min(X[M+1-RJ],Y[N+1]-Y[RI])
					   : (X[RJ]-X[J])*0.5 + min(X[J],Y[I]) + min(X[M+1]-X[RJ],Y[N+1]-Y[RI]));
      if(overlap < min(Y[N+1],X[M+1]) * MinOverlapRatio){
	myfiltered++;
	if(palign->trueTP > 0)
	  mycntTP--;
	alignment[numaligns + i] = NULL;// defer removal of filtered alignment
	continue;
      }
    }
  } // omp for(i = 0; i < numLL; i++)

  #pragma omp critical
  {
    mapfiltered += mymapfiltered;
    threshcnt += mythreshcnt;
    repeatcnt += myrepeatcnt;
    cntTP += mycntTP;
    chimcnt += mychimcnt;
    chimcntTP += mychimcntTP;
    chimcntFP += mychimcntFP;
    filtered += myfiltered;
  }
  } // omp parallel

  if(VERB>=2){
    printf("%lu Alignments parsed:time=%0.6f(wall time=%0.6f)\n",num,mtime(),wtime());
    fflush(stdout);
  }

  size_t orignumaligns = numaligns;
  numaligns += num;

  block->start = orignumaligns;
  block->end = numaligns;

  if(DEBUG){
    extern Calign_block *align_block[MAX_ALIGNBLOCKS];
    assert(block == align_block[num_align_blocks-1]);
  }

  size_t n = align_blockcompact(num_align_blocks-1, orignumaligns, numpairs, numpairs_cum, 0 /* WAS AlignScore */, numthreads, alignment);
  numaligns = orignumaligns + n;

  delete [] numpairs_cum;
  delete [] Yid;
}

void input_align(int num_files, char *filename[])
{
  register char *pt;
  char *qt;
  FILE *fp;
  
  int numthreads = 1;
#ifdef _OPENMP
  numthreads = omp_get_max_threads();
  numthreads = min(numthreads,MaxThreads);
  numthreads = max(1,min(numthreads, nummaps/512));
#endif

  if(VERB>=2){
    printf("starting input_align:num_files=%d:wall time=%0.6f\n",num_files,wtime());
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  long long maxmapid = 0, minmapid = -1;

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {
    long long mymaxmapid = 0, myminmapid = -1;

    #pragma omp for schedule(static,512)
    for(int i= 0; i < nummaps; i++){
      if(Gmap[i]->id > mymaxmapid)
	mymaxmapid = Gmap[i]->id;
      if(myminmapid < 0 || Gmap[i]->id < myminmapid)
	myminmapid = Gmap[i]->id;
      Gmap[i]->mapid = i;    /* make sure Gmap[i]->mapid == i */
    }

    #pragma omp critical
    {
      maxmapid = max(mymaxmapid,maxmapid);
      if(minmapid < 0 || myminmapid < minmapid)
	minmapid = myminmapid;
    }
  }

  if(VERB>=2){
    printf("minmapid=%lld,maxmapid=%lld,nummaps=%d:wall time=%0.6f\n",minmapid,maxmapid,nummaps,wtime());
    /*    int mapid = 111663;
	  printf("mapid=%d:id=%lld,M=%d\n", mapid,Gmap[mapid]->id,Gmap[mapid]->numsite[0]);*/
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  /* set up translation from Gmap[i]->id to Gmap[i]->mapid */
  HashInit();
  for(int i = nummaps; --i >= 0;){
    Cmap *pmap = Gmap[i];
    if(VERB>=2 && (pmap->id == 1506108811301619991LL || pmap->id == 5603437221903379500LL)){
      printf("Calling HashInsert(id=%lld,mapid=%d): file=%s\n",pmap->id,pmap->mapid,vfixx_filename[pmap->fileid]);
      fflush(stdout);
    }
    HashInsert(pmap->id, pmap->mapid);
  }

  if(VERB>=2){
    printf("Initialized HashTable with %d map ids:wall time=%0.6f\n",nummaps,wtime());
    fflush(stdout);
  }

  Calign *palign = 0;
  size_t orignumaligns = numaligns;

  mapsdeleted = 0;/* number of maps deleted (due to -minlen -T or to reduce memory) */

  mapfiltered = 0;/* number alignments filtered due to map being filtered */
  filtered = 0;/* number of filtered alignments */
  threshcnt = 0;
  repeatcnt = 0;/* number of alignments entirely within repeat regions */
  chimcnt = 0;/* total alignments removed as chimeric */
  chimcntTP = 0;/* true positives mistaken for chimeric alignment */
  chimcntFP = 0;/* false positive mistaken for chimeric alignment */
  cntTP = 0;/* true positive alignments : only valid for simulated data */
  
  PChimThreshold = (Pchim > 0.0) ? -log(min(1.0,Pchim)) : 1000.0;

  size_t blen = 1024*1024;
  char *fbuf = (char *)malloc(blen);
  if(!fbuf){
    printf("input_align:malloc(%llu) failed\n",(unsigned long long)blen);
    exit(1);
  }

  origMinLen = MinLen;
  origLogPvThreshold = LogPvThreshold;

  double InvGB = 1.0/(1024.0*1024.0*1024.0);

  FILE *fpbin = NULL;
  char filenamebin[PATH_MAX];
  char *fbufbin = 0;
  char *bufferbin = NULL;
  size_t lenbin = 0;
  size_t *offset_vec = NULL;

  for(int fileno = 0; fileno < num_files; fileno++){
    int binary_format = 0;
    size_t num = 0;/* number of alignments : read in from text file comment (if binary_format)*/

    if((fp = fopen(filename[fileno],"r"))==NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("Failed to read %d'th input file %s:errno=%d:%s\n",fileno+1, filename[fileno], eno,err);
      exit(1);
    }
    setbuffer(fp, fbuf, blen);
  
    if(VERB){
      printf("Reading input file %d of %d (%s) : previous alignments = %llu:wall time=%0.6f\n",fileno+1,num_files,filename[fileno], (unsigned long long)numaligns,wtime());
      fflush(stdout);
    }

    int linecnt = 1;

    for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++){
   
      for(int iter=0; AlignmentFilter > 0 && (AlignmentsPerGb <= 1.0 || numaligns > ((size_t)nummaps) * AlignmentFilter || numaligns > MaxMem * AlignmentsPerGb) && iter < 10; iter++){
	double MemUse = 0.0;
	long long VmSize = 0, VmRss = 0, VmSwap = 0;
	if(AlignmentsPerGb <= 1.0){/* this value is actually the fraction of maxmem that can be used during Assembler_input() */
	  if(linecnt > 1)/* avoid calling getmem() too often */
	    break;
	  /* determine actual memory used */
	  getmem(VmSize,VmRss,VmSwap);
	  MemUse = (double)(VmRss + VmSwap) * InvGB;
	  if(MemUse < MaxMem * AlignmentsPerGb)
	    break;
	}

	if(VERB){
	  if(AlignmentsPerGb <= 1.0){
	    printf("Too many alignments for -maxmem %0.1f x %0.3f (maps=%d,numaligns=%lld,VmRSS=%0.4f,VmSwap=%0.4f Gb): Changing -minlen from %0.2f to %0.2f and log10PV from %0.2f to %0.2f\n",
		   MaxMem, AlignmentsPerGb, nummaps, (long long)numaligns, VmRss * InvGB, VmSwap * InvGB, MinLen, MinLen + minlenDelta, LogPvThreshold, LogPvThreshold + LogPvDelta);
	  } else if(numaligns > ((size_t)nummaps) * AlignmentFilter)
	    printf("Too many alignments per map (maps=%d,numaligns=%lld,MaxMem=%0.1f): Changing -minlen from %0.2f to %0.2f and log10PV from %0.2f to %0.2f\n",
		   nummaps, (long long)numaligns, MaxMem, MinLen, MinLen + minlenDelta, LogPvThreshold, LogPvThreshold + LogPvDelta);
	  else
	    printf("Too many alignments for -maxmem %0.1f (maps=%d,numaligns=%lld): Changing -minlen from %0.2f to %0.2f and log10PV from %0.2f to %0.2f\n",
		   MaxMem, nummaps, (long long)numaligns, MinLen, MinLen + minlenDelta, LogPvThreshold, LogPvThreshold + LogPvDelta);
	  if(VERB>=2)
	    dumpmemmap();
	  fflush(stdout);
	}

	if(DEBUG>=2)
	  for(size_t i = 0; i < numaligns; i++)
	    if(DEBUG && alignment[i]==0){
	      printf("i=%lu,numaligns=%lu,alignment[i]=0 !\n",i,numaligns);
	      fflush(stdout);
	      assert(alignment[i] != 0);
	    }

	MinLen += minlenDelta;
	LogPvThreshold += LogPvDelta;
	
	/* filter molecules */
	int *newmapid = new int[nummaps];/* mapping from old to new locations in Gmap[] (-1 if mapid was deleted) */

	int j = 0;
	for(int i = 0; i < nummaps; i++){
	  Cmap *p = Gmap[i];
	  if(p->site[0][p->numsite[0]+1] < MinLen){
	    newmapid[i] = -1;
	    p->allfree();
	    p->init();
	    continue;
	  }
	  newmapid[i] = j;
	  if(j < i){
	    Cmap *tmp = Gmap[j];
	    p->mapid = j;
	    Gmap[j] = p;
	    tmp->mapid = i;
	    Gmap[i] = tmp;
	  }
	  j++;
	}
	if(VERB){
	  printf("  Reduced number of maps from %d to %d due to MinLen=%0.1f\n",nummaps,j,MinLen);
	  fflush(stdout);
	}
	mapsdeleted += nummaps - j;
	nummaps = j;

	/* filter out alignments based on missing maps due to updated MinLen OR LogPvThreshold */
	/* Also remap mapids in alignments */
	size_t n = 0;
	for(size_t i = 0; i < numaligns; i++){
	  Calign *p = alignment[i];
	  if(p->logPV < LogPvThreshold){
	    threshcnt++;
	    if(p->trueTP > 0)
	      cntTP--;
	    alignment[i] = NULL;
	    continue;
	  }
	  if(newmapid[p->mapid1] < 0 || newmapid[p->mapid2] < 0){
	    mapfiltered++;
	    if(p->trueTP > 0)
	      cntTP--;
	    alignment[i] = NULL;
	    continue;
	  }
	  p->mapid1 = newmapid[p->mapid1];
	  p->mapid2 = newmapid[p->mapid2];
	  n++;
	}
	if(VERB){
	  printf("  Reduced number of alignments so far from %lld to %lld due to MinLen=%0.1f,LogPvThreshold=%0.2f\n",(long long)numaligns,(long long)n,MinLen,LogPvThreshold);
	  fflush(stdout);
	}

	/* compact all align_block[0.. num_align_blocks-1] */
	double tstart = wtime();
	if(DEBUG) assert(align_block[0]->start == 0);
	if(DEBUG) assert(align_block[num_align_blocks - 1]->end == numaligns);
	size_t m = 0;
	size_t maxaligns = 0;
	for(int i = 0; i < num_align_blocks; i++){
	  maxaligns = max(maxaligns, align_block[i]->end - align_block[i]->start);
	  if(DEBUG && i > 0) assert(align_block[i]->start == align_block[i-1]->end);
	}
	int *numpairs = new int[maxaligns];
	size_t *numpairs_cum = new size_t[maxaligns];
	for(int i = 0; i < num_align_blocks; i++){
          size_t num = align_blockcompact(i, m, numpairs, numpairs_cum, 0 /* WAS30 AlignScore */, numthreads, alignment);
	  m += num;
	}
	delete [] numpairs;
	delete [] numpairs_cum;
	if(DEBUG) assert(m == n);
	if(VERB /* HERE >=2 */){
	  double tend = wtime();
	  long long VmSize = 0, VmRss = 0, VmSwap = 0;	  
	  getmem(VmSize,VmRss,VmSwap);
	  printf("\t Compacted alignment memory: VmRSS=%0.4f, VmSwap=%0.4f Gb: wall time= %0.6f (cum = %0.6f secs)\n", 
		 VmRss*InvGB,VmSwap*InvGB, tend - tstart, tend);
	  if(VERB>=2)
	    dumpmemmap();
	  fflush(stdout);
	}

	numaligns = n;

	if(DEBUG>=2)
	  for(size_t i = 0; i < numaligns; i++)
	    if(DEBUG && alignment[i]==0){
	      printf("i=%lu,numaligns=%lu,alignment[i]=0 !\n",i,numaligns);
	      fflush(stdout);
	      assert(alignment[i] != 0);
	    }

	delete [] newmapid;

	/* recompute translation from Gmap[i]->id to Gmap[i]->mapid */
	HashInit();
	for(int i = nummaps; --i >= 0;){
	  Cmap *pmap = Gmap[i];
	  HashInsert(pmap->id, pmap->mapid);
	}
      }/* memory overflow filter loop */

      int len = strlen(buf);
      if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],filename[fileno],linecnt,buf);
	exit(1);
      }

      if(buf[0] == '#'){/* comment lines are ignored, except for "Fixx File Input:"  or "binary format" "alignments=" */
	char *key = (char *)"Fixx File Input:";
	if((pt = strstr(&buf[1],key))){
	  pt += strlen(key);
	  while(*pt && isspace(*pt))
	    pt++;
	  if(!*pt || isspace(*pt)){
	    fprintf(stderr,"Syntax error on Fixx File Input: line in %s\n%s\n", filename[fileno],buf);
	    exit(1);
	  }
	  for(qt=pt; *qt && !isspace(*qt);)
	    qt++;
	  *qt = 0;
	  /*	  if(!strstr(pt,vfixx_filename[0])){
	    printf("WARNING : Fixx File Input=%s specified in %s does not match command line map input (%s ...)\n",
		   pt,align_filename,vfixx_filename[0]);
	    fflush(stdout);
	    }*/
	  continue;
	}

	key = (char *)"binary format";
	if((pt = strstr(&buf[1],key))){
	  binary_format = 1;
	  sprintf(filenamebin,"%s.bin",filename[fileno]);
	  if((fpbin = fopen(filenamebin,"rb"))==NULL){
	    int eno = errno;
	    char *err = strerror(errno);
	    printf("failed to open file %s for reading alignment binary data:errno=%d:%s\n",filenamebin,eno,err);
	    exit(1);
	  }
	  if(!fbufbin){
	    fbufbin = (char *)malloc(blen);
	    if(!fbufbin){
	      printf("input_align:malloc(%llu) failed\n",(unsigned long long)blen);
	      exit(1);
	    }
	  }
	  setbuffer(fpbin, fbufbin, blen);

	  size_t pos_cur = ftell(fpbin);
	  fseek(fpbin, 0L, SEEK_END);
	  size_t pos_end = ftell(fpbin);
	  fseek(fpbin, pos_cur, SEEK_SET);
	  
	  if(!(bufferbin = (char *)malloc(pos_end-pos_cur+1))){
	    printf("input_align: malloc(%llu) failed\n",(unsigned long long)(pos_end-pos_cur+1));
	    fflush(stdout);exit(1);
	  }
	  lenbin = fread(bufferbin, 1, pos_end-pos_cur, fpbin);
	  if(lenbin != pos_end - pos_cur){
	    printf("input_align:fread of %llu bytes from %s failed (return value= %llu)\n", (unsigned long long)(pos_end - pos_cur), filenamebin, (unsigned long long)lenbin);
	    fflush(stdout);exit(1);
	  }
	  bufferbin[lenbin] = 0;
	  fclose(fpbin);
	  fpbin = NULL;
	  if(VERB>=2){
	    printf("Read binary file %s into memory buffer:wall time=%0.6f\n",filenamebin,wtime());
	    fflush(stdout);
	  }

	  continue;
	}
	
	key = (char *)"alignments=";	
	if((pt = strstr(&buf[1],key))){
	  if(!binary_format)
	    continue;
	  pt += strlen(key);
	  num = strtoll(pt,&qt,10);
	  if(pt==qt || !(*qt==0 || isspace(*qt))){
	    printf("Syntax error on line %d of %s : expected number of alignments\n%s\n",linecnt,filename[fileno],pt);
	    exit(1);
	  }
	  
	  if(num > 0){
	    /* read in byte offset array from fp */
	    if(offset_vec)
	      free(offset_vec);
	    offset_vec = (size_t *)malloc(num * sizeof(size_t));
	    size_t r;
	    if((r = fread(offset_vec, sizeof(size_t), num, fp)) != num){
	      printf("input_align: Failed to read %llu byte offsets after line %d in %s (return value = %llu)\n",
		     (unsigned long long)num, linecnt, filename[fileno], (unsigned long long)r);
	      exit(1);
	    }
	    if(VERB>=2){
	      printf("Read in byte offset array from %s into memory buffer(num=%lu):wall time=%0.6f\n",filename[fileno],num,wtime());
	      /*	    for(size_t i = 0; i < num; i++)
			    printf("offset[%lu] = %lu\n",i, offset_vec[i]); */
	      fflush(stdout);
	    }
	  
	    /* complete binary file input */
	    input_align_binary(numaligns, num, offset_vec, bufferbin, lenbin, filenamebin, fileno);
	  
	    free(offset_vec);
	    offset_vec = 0;
	  }

	  free(bufferbin);
	  bufferbin = 0;

	  break;/* done with this file */
	}

	continue;
      }

      long long Mol0ID, Mol1ID;

      if(buf[0] == '>'){
	if(buf[1] != '0'){
	  printf("1: Syntax error on line %d of %s(expected >0):\n%s\n",linecnt,filename[fileno],buf);
	  exit(1);
	}
	if(VERB>=2){
	  printf("Reading Header line %d of %s:\n%s\n",linecnt,filename[fileno],buf);
	  fflush(stdout);
	}

	/* read in the next alignment in alignment[numaligns] */
	maxalignallocNULL(numaligns+1,alignment,numaligns,maxaligns);
	if(!alignment[numaligns])
	  alignment[numaligns] = new Calign[colors==1 ? 1 : colors+1];
	palign = alignment[numaligns++];
	palign->allfree();

	palign->Lend = palign->Rend = -2;
	palign->fileid = fileno;

	pt = &buf[2];
	long long AlignmentID = strtoll(pt,&qt,10);
	if(pt == qt || AlignmentID <= 0){
	  printf("Invalid integer value for AlignmentID on line %d of %s:\n%s\n",linecnt,filename[fileno],buf);
	  exit(1);
	}
	if(AlignmentID != (long long)(numaligns+mapfiltered)){
	  if(AlignmentID_ordered && !threshcnt){
	    printf("WARNING: Out of order AlignmentID=%lld on line %d of %s (numaligns=%llu,mapfiltered=%llu)\n",AlignmentID,linecnt,filename[fileno],(unsigned long long)numaligns,(unsigned long long)mapfiltered);
	    AlignmentID_ordered = 0;
	  }
	}

	pt = qt;
	Mol0ID =  strtoll(pt,&qt,10);
	if(pt == qt || !(*qt==0 || isspace(*qt))){
	  printf("Invalid Mol0ID syntax on line %d of %s:\n%s\n",linecnt,filename[fileno],buf);
	  exit(1);
	}
	int filter = 0;
	if(Mol0ID > maxmapid || Mol0ID < minmapid){
	  if(MinLen <= 0 && MinSites <= 1){
	    printf("Mol0ID value %lld on line %d of %s outside known ID range of %lld .. %lld:\n%s\n",Mol0ID, linecnt,filename[fileno],minmapid, maxmapid,buf);
	    exit(1);
	  }
	  filter = 1;
	} else {
	  int index = HashFindID(Mol0ID);
	  if(index < 0){
	    if(MinLen <= 0 && MinSites <= 1){
	      printf("Mol0ID value=%lld on line %d of %s not found in %d nanomaps from %s ... \n",Mol0ID,linecnt,filename[fileno],nummaps, vfixx_filename[0]);
	      exit(1);
	    }
	    filter = 1;
	  }
	  if(VERB>=2){
	    printf("line%d of %s:Mol0ID=%lld,mapid=%d\n",linecnt,filename[fileno],Mol0ID,index);
	    fflush(stdout);
	  }
	  palign->mapid1 = index;
	}

	pt = qt;
	Mol1ID = strtoll(pt,&qt,10);
	if(pt == qt || !(*qt==0 || isspace(*qt))){
	  printf("Invalid Mol1ID syntax on line %d of %s:\n%s\n",linecnt,filename[fileno],buf);
	  exit(1);
	}
	if(Mol1ID > maxmapid || Mol1ID < minmapid){
	  if(MinLen <= 0 && MinSites <= 1){
	    printf("Mol1ID value %lld on line %d of %s outside known ID range %lld .. %lld:\n%s\n", Mol1ID, linecnt,filename[fileno],minmapid,maxmapid,buf);
	    exit(1);
	  }
	  filter = 1;
	} else {
	  int index = HashFindID(Mol1ID);
	  if(index < 0){
	    if(MinLen <= 0 && MinSites <= 1){
	      printf("Mol1ID value=%lld on line %d of %s not found in %d nanomaps from %s ...\n",Mol1ID,linecnt,filename[fileno],nummaps, vfixx_filename[0]);
	      exit(1);
	    }
	    filter = 1;
	  }
	  if(VERB>=2){
	    printf("line%d of %s:Mol1ID=%lld,mapid=%d\n",linecnt,filename[fileno],Mol1ID,index);
	    fflush(stdout);
	  }
	  palign->mapid2 = index;
	}

	if(filter){/* either map was filtered out : so skip the remaining alignment information */
	  if(VERB>=2){
            printf("Skipping alignment on lines %d..%d of %s since map was filtered out\n",linecnt, linecnt + (AlignScore ? 4 : 2), filename[fileno]);
	    fflush(stdout);
	  }

	  /* go to next line and read in the aligned sites of the first map */
	  linecnt++;
	  if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	    printf("Mol0 alignment line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	    exit(1);
	  }

	  /* go to nextline and read in the aligned sites of the 2nd map */
	  linecnt++;
	  if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	    printf("Mol1 alignment line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	    exit(1);
	  }

	  if(AlignScore){
	    linecnt++;
	    if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	      printf("Alignment score line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	      exit(1);
	    }
	    linecnt++;
	    if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	      printf("Raw Alignment score line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	      exit(1);
	    }
	  }

	  mapfiltered++;
	  numaligns--;
	  continue;
	}

	pt = qt;
	palign->score = strtod(pt,&qt);
	if(pt == qt){
	  printf("Invalid Score value on line %d of %s:\n%s\n",linecnt,filename[fileno],buf);
	  exit(1);
	}
      
	pt = qt;
	double Offset = strtod(pt,&qt);
	if(pt == qt){
	  printf("Invalid Offset value=%e on line %d of %s:\n%s\n",Offset,linecnt,filename[fileno],buf);
	  exit(1);
	}

	pt = qt;
	double Overlap = strtod(pt,&qt);
	if(pt == qt){
	  printf("Invalid Overlap value=%e on line %d of %s:\n%s\n",Overlap,linecnt,filename[fileno],buf);
	  exit(1);
	}

	pt = qt;
	int orientation = strtol(pt,&qt,10);
	if(pt == qt || abs(orientation) != 1){
	  printf("Invalid Orientation value on line %d of %s:\n%s\n",linecnt,filename[fileno],buf);
	  exit(1);
	}
	palign->orientation = (orientation < 0) ? 1 : 0;

	pt = qt;
	palign->logPV = strtod(pt,&qt);
	if(pt == qt){
	  printf("Invalid PvalueLog10 value on line %d of %s:\n%s\n",linecnt,filename[fileno],buf);
	  exit(1);
	}

	pt = qt;
	double TrueOffset = strtod(pt,&qt);
	if(pt == qt){
	  printf("Invalid TrueOffset value on line %d of %s:\n%s\n",linecnt,filename[fileno],buf);
	  exit(1);
	}

	pt = qt;
	double TrueOverlapFraction = strtod(pt,&qt);
	if(pt == qt){
	  printf("Invalid TrueOverlapFraction value on line %d of %s:\n%s\n",linecnt,filename[fileno],buf);
	  exit(1);
	}
	palign->trueTP = (TrueOffset==0.0) ? -1 : (TrueOverlapFraction > 0.0) ? 1 : 0;

	/* go to next line and read in the aligned sites of the first map */
	linecnt++;
	if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	  printf("Mol0 alignment line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	  exit(1);
	}
	if(buf[0] == '>'){
	  printf("syntax error on Mol0 alignment line on line %d of file %s (line cannot start with >):\n%s\n",
		 linecnt,filename[fileno],buf);
	  exit(1);
	}
	if(VERB>=2){
	  printf("Reading line %d of %s (Map0 aligned sites):\n%s\n",linecnt,filename[fileno],buf);
	  fflush(stdout);
	}

	/* scan over input line to count how many aligned sites are present */
	palign->numpairs = 0;
	int sum = 0;
	for(pt= buf; *pt; pt = qt){
	  sum +=  strtol(pt,&qt,10);
	  if(pt == qt)
	    break;
	  palign->numpairs++;
	}
	while(*qt && isspace(*qt))
	  qt++;
	if(*qt){
	  printf("Invalid last char(%c) at position %ld on line %d in %s:\n%s\n",*qt,(long)(qt-buf),linecnt,filename[fileno],buf);
	  exit(1);
	}
      
	palign->expand_arrays(palign->numpairs, AlignScore);

	if(!MaxCoverageDump && hash_filename){/* skip the actual alignment information */
	  /* go to nextline and read in the aligned sites of the 2nd map */
	  linecnt++;
	  if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	    printf("Mol1 alignment line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	    exit(1);
	  }
	  if(AlignScore){
	    linecnt++;
	    if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	      printf("Alignment score line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	      exit(1);
	    }
	    linecnt++;
	    if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	      printf("Raw Alignment score line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	      exit(1);
	    }
	  }
	  continue;
	}

	/* now read in the Map0 sites */
	int site = 0;
	for(pt= buf; *pt; pt = qt){
	  int siteindex = strtol(pt,&qt,10);
	  if(pt == qt)
	    break;
	  if(siteindex < 0 || siteindex >= Gmap[palign->mapid1]->numsite[0]){
	    printf("Invalid site index %d for mapid=%d(id=%lld, M=%d, file=%s) on line %d of %s (Mol0ID=%lld,Mol1ID=%lld):\n%s\n",
		   siteindex,palign->mapid1,Gmap[palign->mapid1]->id,Gmap[palign->mapid1]->numsite[0], vfixx_filename[Gmap[palign->mapid1]->fileid], linecnt, filename[fileno], Mol0ID, Mol1ID, buf);
	    exit(1);
	  }
	  if(DEBUG && !(site < palign->numpairs)){
	    printf("site=%d,palign->mapid1=%d(id=%lld),palign->mapid2=%d(id=%lld),palign->numpairs=%d (line %d of %s):\n%s\n",
		   site,palign->mapid1,Gmap[palign->mapid1]->id,palign->mapid2,Gmap[palign->mapid2]->id,palign->numpairs,linecnt,filename[fileno],buf);
	    assert(site < palign->numpairs);
	  }
	  palign->sites1[site++] = siteindex+1;
	}

	/*	if(Mol0ID == 2245000294L && Mol1ID == 8606000556LL){
		int mapid1 = 111663;
		printf("Mol0ID=%lld,Mol1ID=%lld:aligned sites1=%d: palign->sites1[0]=%d, palign->sites1[%d]=%d\n",
		Mol0ID,Mol1ID,site, palign->sites1[0]=%d,site-1,palign->sites1[site-1]);
		fflush(stdout);
		}*/

	/* go to nextline and read in the aligned sites of the 2nd map */
	linecnt++;
	if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	  printf("Mol1 alignment line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	  exit(1);
	}
	if(buf[0] == '>'){
	  printf("syntax error on Mol1 alignment line on line %d of file %s (line cannot start with >):\n%s\n",
		 linecnt,filename[fileno],buf);
	  exit(1);
	}
	if(VERB>=2){
	  printf("Reading line %d of %s (Map1 aligned sites):\n%s\n",linecnt,filename[fileno],buf);
	  fflush(stdout);
	}
      
	/* now read in the Map1 sites */
	site = 0;
	register int M = Gmap[palign->mapid2]->numsite[0];
	for(pt= buf; *pt; pt = qt){
	  int siteindex = strtol(pt,&qt,10);
	  if(pt == qt)
	    break;
	  if(siteindex < 0 || siteindex >= Gmap[palign->mapid2]->numsite[0]){
	    printf("Invalid site index %d for mapid=%d(id=%lld, M=%d, file=%s) on line %d of %s (Mol0ID=%lld,Mol1ID=%lld):\n%s\n",
		   siteindex,palign->mapid2,Gmap[palign->mapid2]->id,Gmap[palign->mapid2]->numsite[0], vfixx_filename[Gmap[palign->mapid2]->fileid], linecnt, filename[fileno], Mol0ID, Mol1ID, buf);
	    exit(1);
	  }
	  assert(site < palign->numpairs);
	  palign->sites2[site++] = palign->orientation ? M-siteindex : siteindex+1;
	}

	if(AlignScore){

	  /* go to nextline and read in the segment scores  */
	  linecnt++;
	  if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	    printf("Alignment score line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	    exit(1);
	  }
	  if(buf[0] == '>'){
	    printf("syntax error on Alignment Score line on line %d of file %s (line cannot start with >):\n%s\n",
		   linecnt,filename[fileno],buf);
	    exit(1);
	  }
	  if(VERB>=2){
	    printf("Reading line %d of %s (iscore[] values):\n%s\n",linecnt,filename[fileno],buf);
	    fflush(stdout);
	  }

	  /* now read in the iscore[] values */
	  site = 0;
	  for(pt= buf; *pt; pt = qt){
	    double score = strtod(pt,&qt);
	    if(pt == qt)
	      break;
	    assert(site <= palign->numpairs);
	    palign->iscore[site++] = score;
	  }
      
	  /* go to nextline and read in the raw segment scores  */
	  linecnt++;
	  if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	    printf("Raw Alignment score line missing, too long or not terminated in input file %s on line %d\n",filename[fileno],linecnt);
	    exit(1);
	  }
	  if(buf[0] == '>'){
	    printf("syntax error on Raw Alignment Score line on line %d of file %s (line cannot start with >):\n%s\n",
		   linecnt,filename[fileno],buf);
	    exit(1);
	  }
	  if(VERB>=2){
	    printf("Reading line %d of %s (outscore[] values):\n%s\n",linecnt,filename[fileno],buf);
	    fflush(stdout);
	  }

	  /* now read in the outscore[] values */
	  site = 0;
	  for(pt= buf; *pt; pt = qt){
	    double score = strtod(pt,&qt);
	    if(pt == qt)
	      break;
	    assert(site <= palign->numpairs);
	    palign->outscore[site++] = score;
	  }
	}
      } else {/* buf[0] != '>' */
	/* Syntax error : Perhaps we can ignore lines until next line starting with '>' */
	printf("2: Syntax error on line %d of %s(expected >0):\n%s\n",linecnt,filename[fileno],buf);
	exit(1);
      }

      /* filter early to reduce memory footprint */
      if(!(palign->score > ScoreThreshold && 
	   palign->logPV > LogPvThreshold &&
	   palign->numpairs >= AlignedSiteThreshold &&
	   AlignedLength(palign,Gmap) >= AlignedLengthThreshold)){

	threshcnt++;

	/* recycle palign */
	palign->allfree();
	numaligns--;
	continue;
      }

      if(AlignScore && RepeatMaxShift > 0){/* filter output alignments entirely within a repeat region */
	bool repeat = palign->IsRepeatRegion();
	if(repeat){
	  repeatcnt++;
	    
	  /* recycle palign */
	  numaligns--;
	  if(DEBUG) assert(palign == alignment[numaligns]);
	  palign->allfree();
	  continue;
	}
      }

      if(DEBUG>=2) assert(palign->mapid1 < nummaps);
      if(DEBUG>=2) assert(palign->mapid2 < nummaps);

      /* filter out poor and chimeric alignments */
      //      char *name1 = Gmap[palign->mapid1]->name;
      //      char *name2 = Gmap[palign->mapid2]->name;

      register double cum = 0.0;
      register double hwm = 0.0;
      register int Khwm = -1;
      register double maxdrop = 0.0;
      register int start = -1, end = -1;
      if(AlignScore){
	for(register int k = 0; k <= palign->numpairs;k++){
	  register double score = palign->iscore[k];
	  cum += score;
	  if(cum > hwm){
	    hwm = cum;
	    Khwm = k;
	  } else if(hwm - cum > maxdrop){
	    maxdrop = hwm - cum;
	    start = Khwm;
	    end = k;
	  }
	}
      }

      if(palign->trueTP >= 0){
	long long id1 = Gmap[palign->mapid1]->id;
	long long id2 = Gmap[palign->mapid2]->id;
	if(id1 > (LOC_WIDTH*LEN_WIDTH*FLAG_WIDTH) && id2 > (LOC_WIDTH*LEN_WIDTH*FLAG_WIDTH)){// see test_generator.cpp
	  if((id1 & 0x4) || (id2 & 0x4))/* chimeric molecules are assumed to never align correctly */
	    palign->trueTP = 0;
	  else {
	    bool TrueFlip = (id1 & 0x1) != (id2 & 0x1);
	    if(TrueFlip != (palign->orientation))
	      palign->trueTP = 0;
	  }
	} else
	  palign->trueTP = -1;
	cntTP += (palign->trueTP >= 1) ? 1 : 0;
      }

      /* check for chimeric alignments with leftmost/rightmost misaligned sites exceed ENDSITES */
      /* postpone this till later ? graph_trim() can do a better job analyzing each maps local coverage based on all its pairwise alignments
	 and make a better recognition of both chimeric maps as well as maps with regions that do not align : such maps need to be trimmed off,
	 rather than have their alignments deleted. */
      if((min(palign->sites1[0],palign->sites2[0]) > ENDSITES ||
	  min(Gmap[palign->mapid1]->numsite[0] - palign->sites1[palign->numpairs-1],
	      Gmap[palign->mapid2]->numsite[0] - palign->sites2[palign->numpairs-1]) > ENDSITES) ||
	 (AlignScore && maxdrop > PChimThreshold && end - start >= Pchim_minlen)){

	chimcnt++;

	/* check if we deleted a true positive */
	if(palign->trueTP >= 0){
	  if(palign->trueTP > 0){
	    chimcntTP++;
	    cntTP--;
	  } else
	    chimcntFP++;
	}
	filtered++;
	  
	/* recyle palign */
	palign->allfree();
	numaligns--;
	continue;
      }

      if(MinOverlapRatio > 0.0){
	Cmap *Ymap = Gmap[palign->mapid1];
	Cmap *Xmap = Gmap[palign->mapid2];
	FLOAT *Y = Ymap->site[0];
	FLOAT *X = Xmap->site[0];
	int N = Ymap->numsite[0];
	int M = Xmap->numsite[0];
	int I = palign->sites1[0];
	int J = palign->sites2[0];
	int U = palign->numpairs;
	int RI = palign->sites1[U-1];
	int RJ = palign->sites2[U-1];
	double overlap = 0.5*(Y[RI]-Y[I]) + (palign->orientation ? (X[M+1-J] - X[M+1-RJ])*0.5 + min(X[M+1]-X[M+1-J],Y[I]) + min(X[M+1-RJ],Y[N+1]-Y[RI])
					     : (X[RJ]-X[J])*0.5 + min(X[J],Y[I]) + min(X[M+1]-X[RJ],Y[N+1]-Y[RI]));
	if(overlap < min(Y[N+1],X[M+1]) * MinOverlapRatio){
	  filtered++;
	  if(palign->trueTP > 0)
	    cntTP--;
	  /* recyle palign */
	  palign->allfree();
	  numaligns--;
	  continue;
	}
      }
    }// for(; fgets(buf,LINESIZ,fp) != NULL; linecnt++)

    fclose(fp);

    if(DEBUG>=2)
      for(size_t i = 0; i < numaligns; i++)
	if(DEBUG && alignment[i]==0){
	  printf("i=%lu,numaligns=%lu,alignment[i]=0 !\n",i,numaligns);
	  fflush(stdout);
	  assert(alignment[i] != 0);
	}
  }

  if(fbuf){
    free(fbuf);
    fbuf = 0;
  }
  if(fbufbin){
    free(fbufbin);
    fbufbin = 0;
  }
  
  if(VERB){
    printf("Read in %llu alignments from %d files %s ... (%llu mapfiltered, %llu thresholded, %llu RepeatMasked, %llu chimeric(FP=%llu,TP=%llu)) %llu remaining (TP=%llu):wall time=%0.6f\n",
	   (unsigned long long)(numaligns - orignumaligns + filtered+threshcnt+repeatcnt+mapfiltered), num_files,filename[0], 
	   (unsigned long long)mapfiltered, (unsigned long long)threshcnt, (unsigned long long)repeatcnt, (unsigned long long)chimcnt, (unsigned long long)chimcntFP, 
	   (unsigned long long)chimcntTP, (unsigned long long)(numaligns - orignumaligns), (unsigned long long)cntTP, wtime());
    printf("Alignment size = %lu bytes, previous alignments = %lu, number of maps deleted= %d\n",sizeof(Calign),orignumaligns,mapsdeleted);
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }
  if(DEBUG>=2){
    for(size_t i = orignumaligns; i < numaligns; i++){
      Calign *p = alignment[i];
      if(!(p->Lend == -2 && p->Rend == -2 && 0 <= p->fileid && p->fileid < num_files)){
	printf("alignment[%lu]: Lend= %d, Rend=%d, fileid= %d : mapid1=%d,mapid2=%d,or=%d:numpairs=%d,score= %0.6f, logPV= %0.2f (numaligns= %lu .. %lu,num_files=%d)\n",
	       i,p->Lend,p->Rend,p->fileid,p->mapid1,p->mapid2,p->orientation,p->numpairs,p->score,p->logPV,orignumaligns,numaligns,num_files);
	fflush(stdout);
	assert(p->Lend == -2);
	assert(p->Rend == -2);
	assert(0 <= p->fileid && p->fileid < num_files);
      }
    }
  }

  // compress alignments[] and Gmap[] data to save memory
  double tstart = wtime();

  int maps_compacted = 0;
  if(mapsdeleted > 1000){
    mapcompact(nummaps, maxmaps, Gmap);
    maps_compacted = mapsdeleted;
  }

  /* compact all align_block[0.. num_align_blocks-1] */  
  size_t aligndeleted = 0;
  for(size_t i = 0; i < numaligns; i++)
    if(!alignment[i])
      aligndeleted++;

  size_t m = numaligns;

  if(aligndeleted > 1000){
    if(DEBUG) assert(align_block[0]->start == 0);
    if(DEBUG) assert(align_block[num_align_blocks - 1]->end == numaligns);
    m = 0;
    size_t maxaligns = 0;
    for(int i = 0; i < num_align_blocks; i++){
      maxaligns = max(maxaligns, align_block[i]->end - align_block[i]->start);
      if(DEBUG && i > 0) assert(align_block[i]->start == align_block[i-1]->end);
    }  
    int *numpairs = new int[maxaligns];
    size_t *numpairs_cum = new size_t[maxaligns];
    for(int i = 0; i < num_align_blocks; i++){
      size_t num = align_blockcompact(i, m, numpairs, numpairs_cum, 0 /* AlignScore */, numthreads,alignment);
      m += num;
    }
    delete [] numpairs;
    delete [] numpairs_cum;
    if(DEBUG) assert(m <= numaligns);
  }

  if(VERB /* HERE >=2 */ && (aligndeleted > 1000 || maps_compacted > 0)){
    double tend = wtime();
    long long VmSize = 0, VmRss = 0, VmSwap = 0;	  
    getmem(VmSize,VmRss,VmSwap);

    printf("\t Compacted Gmap[] and alignment[] memory: VmRSS=%0.4f, VmSwap=%0.4f Gb, maps compacted=%d/%d, nummaps= %d, numaligns=%lu -> %lu : wall time= %0.6f (cum = %0.6f secs)\n", 
      VmRss*InvGB,VmSwap*InvGB, maps_compacted, mapsdeleted, nummaps, numaligns, m,tend - tstart, tend);
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  numaligns = m;

  AlignScore = 0;// alignment scores (iscore,outscore) have been deleted after using them to filter alignments
}
