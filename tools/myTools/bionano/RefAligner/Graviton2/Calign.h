#ifndef CALIGN_H
#define CALIGN_H

#define GOODMAPS 1 /* save detailed alignment information only for good alignments (scoring above 3 thresholds) */

class Calign_block;

  /** Data container for alignment between two maps
  * 
  * NOTE : for multiple colors, there will be one alignment structure per color, plus one structure with combined score 
  * (only mapid1,mapid2,score,logPV, orientation & numpairs will be valid,
  * plus score2,logPV2,orientation2,numpairs2 if SecondBest) 
  */
class Calign {

 public:
  
  int mapid1; ///< index into map[] of 1st map (or index into refmap[] of reference map)
  int mapid2; ///< index into map[] of 2nd map (Nanomap) 
  int fileid; ///< index into align_filename[0..num_align_files-1] (or -1 if alignment was not read in from a file)
  double score; ///< score of best alignment 
  double logPV; ///< -log10(Pvalue), Pvalue without Bonferroni correction 
  int orientation; ///<  orientation of 2nd map in best alignment (0 is normal, 1 = reversed) 
  int scaleID; ///< 2nd map was scaled by ScaleFactor[scaleID] (ScaleFactor[0] == 1.0) 
  //  int color; ///< If colors > 1 : the color of this alignment (the complement alignment will include one alignment for each color) 

  int Lij1,Rij1; ///<first and last misaligned cut in 1st map (If CALIGN_ENDS==0 : for End Outliers, this is the last aligned cut)
  int Lij2,Rij2; ///<first and last misaligned cut in 2nd map (If CALIGN_ENDS==0 : for End Outliers, this is the last aligned cut)

#if CALIGN_END==0
  int LijY,RijY; /* original value of Lij1,Rij1 without End Outlier adjustment : only valid in refalign.cpp */
  int LijX,RijX; /* original value of Lij2,Rij2 without End Outlier adjustment : only valid in refalign.cpp */
#endif

  int Lend,Rend; /**< Specifies end of alignment type:
		   * - 0 : Linear chromosome end (with ends aligned)
		   * - -1 : Normal alignment end (chromosome ends not aligned but NOT a local alignment)
		   * - -2 or -3 : Local alignment end (ends right after last aligned site) : End Outlier */ 
  int numpairs; ///< number of aligned site pairs (>= 0 : 0 means no valid alignment found) 
  int allocated_size;
  int *sites1; ///< sites1[i=0..numpairs-1] are the sites in map[mapid1] that aligned 
  int *sitesK1; ///< sitesK1[i=0..numpairs-1] is the number of additional sites left of sites1[i] that merged with sites1[i] (only defined if spots_filename!=0)
  int *sites2; ///< sites2[i=0..numpairs-1] are the sites in map[mapid2] that aligned : map[mapid2] may be in reversed orientation 

  int Nrange,Mrange;/* number of contiguous labels in Ref and Qry that were considered when performing this alignment (affects Bonferroni correction for logPV) */

  double norm; ///< normalized RMS sizing error per aligned interval : used to flag suspected duplicate maps (with norm < RMS_MIN && numpairs >= RMS_SITES) 

  /* The following fields are used with -Secondbest to point to 2nd best alignment */
  /* Also used with -pairsplit in which case it refers to the best alignment with the opposite orientation */
  Calign *align2;

  /* The following fields are used with -MultiMatch, to hold multiple match groups of which the current alignment is the best one */
  int multicnt;/* number of match groups (only valid of Malign != NULL) */
  int multimax;/* size of allocated array Malign[0..multimax-1] */
  Calign **Malign;/* Malign[m=0..multicnt-1] points to the m'th alignment (for 2 colors each m has a 3-element array Malign[m=0..multicnt-1][0..2], otherwise it is a 1-element array) */
  double logPV2;/* If alignment was split with -RefSplit or CutFlip was used, the original logPV of the full alignment (split alignments are not discarded unless this value is below threshold) */

  FLOAT *iscore; /**< iscore[i=1..numpairs-1] are the scores for aligned interval sites*[i-1..i]
		 * - iscore[0] is the score for the left end of the alignment (includes endoutlier adjustment)
		 * - iscore[numpairs] is the score of the right end of the alignment (includes endoutlier edjustment) */
  FLOAT *outscore; /**< outscore[i=1..numpairs-1] is the value of iscore[i] before outlier adjustment
		      * - if(outscore[i] < iscore[i]) then the aligned interval [i-1..i] is an outlier.
		      * - outscore[0] and outscore[numpairs] correspond to end intervals and currently always match iscore[], see Lend,Rend to check if there is an EndOutlier */
  int noutliers;/* number of internal outliers : sum(i=1..numpairs-1) of (outscore[i] < iscore[i] ? 1 : 0) */
  FLOAT maxoutlier;/* largest sizing error of any internal outlier */
  int maxoutlierLabels;/* largest number of misalinged labels in any internal outlier : not necessilty the same outlier as the one with sizing error = maxoutlier */

  int repeat; ///< If > 0 : A repeat was detected  at this many times the median aligned interval size 

  double trueoffset; ///< true offset between maps (NaN if ground truth is not known)
  double trueoverlap; ///< true overlap between maps (NaN if ground truth is not known) 
  int trueTP; /**< 1 IFF trueoverlap is >= MINOVERLAP * average size of 2 maps, 0 otherwise
		* (-1 if ground truth unknown) */
  int chimpair; /**< 1 IFF map[mapid2]->origmap AND This alignment is tandem with the alignment of map[mapid2]->origmap
		  * with a gap of no more than 20 * Ylambda  */

  int *suppress; /* suppress[i=1..M] : Region of mapid2 for which SVs are suppressed : there is another alignment with mapid2 for the same region and either the other alignment has better logPV (or -indelNoOverlap was used)*/

  int start,end;/* Used as extended ends for center cloned region with PairSplit for query (along with refstart,refend for reference) */
  int refstart,refend;/* Used with PairSplit : see start,end */

  int stitch;/* Used with -MultiMatch && -RefSplitStitch : & 1 IFF left end was truncated due to overlap with higher scoring matchgroup (but no further) and could still be stitched back if original matchgroup is split by RefSplit
	                                                   & 2 IFF right end was truncated due to overlap with higher scoring matchgroup (but no further) and could still be stitched back if original matchgroup is split by RefSplit */
  int rev;/* Use with -MultiMatchRev : 1 IFF ref Y is reversed (orientation = 0 means qry X has same orientation as Y) */

  int block_id;// -1 : each array is allocated seperately, >= 0 : index into align_block[0..num_align_blocks-1] 
  //  Calign_block *blockmem;/* If 0, each array is allocated seperately,
  //			    If != 0, all memory (except align2[], Malign[] and suppress[]) is allocated from a large memory block of type Calign_block */
  double mapWT;/* Used by BestRefWT */

  /* member functions */
  Calign(){
    init();
  }
  Calign(Calign *pinit, int myscore){/* custom copy constructor (if score==1 also copy iscore[] & outscore[]) */
    extern void copy(Calign *dest, Calign *pinit, int myscore, int deepcopy);
    copy(this,pinit,myscore,1);
  }

  ~Calign(){
    allfree();
  }

  inline void init(){
    mapid1 = mapid2 = -1;/* not initialized */
    fileid = -1;
    numpairs = 0;/* no arrays allocated */
    sites1 = NULL;
    sites2 = NULL;
    sitesK1 = NULL;/* not allocated */
    iscore = NULL;/* not allocated */
    outscore = NULL;/* not allocated */
    chimpair = 0;
    scaleID = 0;
    allocated_size = 0;
    repeat = -1;
    block_id = -1;
    align2 = NULL;
    Malign = NULL;
    suppress = NULL;
    mapWT = 1.0;
    logPV = logPV2 = -1.0;
    noutliers = -1;
    maxoutlier = -1.0;
    maxoutlierLabels = -1;
    if(VERB>=3){
      printf("Calign():this=%p initialized\n", this);
      fflush(stdout);
    }
    if(DEBUG) Nrange = Mrange = -1;
#if CALIGN_END == 0
    if(DEBUG) Lij1 = Rij1 = Lij2 = Rij2 = LijY = RijY = LijX = RijX = -1;
#endif
  }

  inline void allfree(){
    if(VERB>=3){
      printf("Calign->allfree():this=%p, allocated_size=%d, mapid1=%d,mapid2=%d,or=%d,Malign=%p\n", 
	     this, allocated_size,mapid1,mapid2,orientation,Malign);
      if(Malign != NULL)
	printf(":multicnt=%d,multimax=%d,Malign[0]=%p\n",multicnt,multimax,Malign[0]);
      printf("\n");
      fflush(stdout);
    }
    if(allocated_size > 0 && block_id < 0){
      delete [] sites1; sites1 = NULL;
      delete [] sitesK1; sitesK1 = NULL;
      delete [] sites2; sites2 = NULL;
      delete [] iscore; iscore = NULL;
      delete [] outscore; outscore = NULL;
    }
    allocated_size = 0;
    if(align2){
      delete [] align2;
      align2 = NULL;
    }
    if(Malign && block_id < 0){
      if(DEBUG) assert(multimax > 0);
      for(int i = 0; i < multimax; i++)
	delete [] Malign[i];
      delete [] Malign;
      Malign = NULL;
    }
    multicnt = 0;
    if(suppress){
      delete [] suppress;
      suppress = NULL;
    }
    numpairs = 0;
    if(DEBUG>=2 && !(numpairs==0 && sitesK1==0 && iscore==0 && outscore==0)){
      printf("numpairs=%d,sitesK1=x%p,iscore=x%p,outscore=x%p\n",
	     numpairs, sitesK1, iscore, outscore);
      assert(numpairs==0 && sitesK1==0 && iscore==0 && outscore==0);
    }
    //    block_id = -1;
    /* NOTE : Do NOT initialize mapid1 etc, since allfree() is used to reallocate alignment memory in pairalign(),refalign(), after mapid1,mapid2 & orientation are recorded */
  }
  
  inline void expand_arrays(int num) {
    if(num <= allocated_size) return;
    if(allocated_size > 0 && block_id < 0){
      delete [] sites1;
      delete [] sitesK1;
      delete [] sites2;
      delete [] iscore;
      delete [] outscore;
    }
    allocated_size = num;
    sites1 = new int[allocated_size];
    sitesK1 = new int[allocated_size];
    sites2 = new int[allocated_size];
    iscore = new FLOAT[allocated_size+1];
    outscore = new FLOAT[allocated_size+1];
    block_id = -1;
    if(VERB>=3){
      printf("Calign->expand_arrays(%d):this=%p,sites1=%p,sites2=%p\n",num,this,sites1,sites2);
      fflush(stdout);
    }
  }

  inline void expand_arrays(int num, int myAlignScore) {// used in input_align() : does NOT allocate sitesK1[]
    if(num <= allocated_size) return;
    if(allocated_size > 0 && block_id < 0){
      delete [] sites1;
      delete [] sitesK1;
      delete [] sites2;
      delete [] iscore;
      delete [] outscore;
    }
    block_id = -1;
    allocated_size = num;
    sites1 = new int[allocated_size];
    sites2 = new int[allocated_size];
    sitesK1 = NULL;
    if(myAlignScore){
      iscore = new FLOAT[allocated_size+1];
      outscore = new FLOAT[allocated_size+1];
    } else
      iscore = outscore = NULL;
    if(VERB>=3){
      printf("Calign->expand_arrays(%d,%d):this=%p,sites1=%p,sites2=%p\n",num,myAlignScore,this,sites1,sites2);
      fflush(stdout);
    }
  }

  int IsRepeatRegion();/* return K > 0 iff the alignment only covers a single repeat region : shifting the alignment by K sites results in similar average sizing errors */
};


extern void maxalignalloc(size_t num, Calign **&alignment, size_t &numaligns, size_t &maxaligns);
extern void maxalignallocNULL(size_t num, Calign **&alignment, size_t &numaligns, size_t &maxaligns, int *alignment_blockid = NULL);
extern void alignfree();/* free up all memory of alignment[] array (plus block memory) */
extern void copy(Calign *dest, Calign *pinit, int score, int deepcopy);
extern void copyrev(Calign *dest, Calign *pinit, int score);/* copy with both map orientation reversed (rev is changed, orientation stays same) */
extern void align_blockfree();/* free up all block memory of alignment[] array */
extern Calign_block *align_blockalloc(size_t num,size_t sum, size_t Msiz, int AlignScore, int sitesK, int &block_id);/* allocate a new Calign_block */
extern size_t align_blockcompact(int block_index, size_t start, int *numpairs, size_t *numpairs_cum, int AlignScore, int numthreads, Calign **alignment);/* compact a Calign_block */
extern void expand_arrays(Calign *p, int num, int AlignScore, Calign_block *block, int block_index);

#define MAX_ALIGNBLOCKS MAXFILES // This should exceed the maximum number of input alignment files */

extern int num_align_blocks;
extern Calign_block *align_block[MAX_ALIGNBLOCKS];

static inline Calign *MultCalign(Calign *pinit, int score, int ncolors){/* copy alignment for multiple colors */
  if(DEBUG) assert(ncolors > 1);
  Calign *p = new Calign[ncolors];
  for(register int c = 0; c < ncolors; c++)
    copy(&p[c],&pinit[c],score,1);
  return p;
}

#include <sys/mman.h>

class Calign_block {
 public:
  size_t numaligns, alignsiz;
  size_t numpairs_sum;/* upper bound on sum of numpairs of block[].numpairs */
  size_t nextint, nextscore;// nextalign;
  size_t intsiz, scoresiz, Malignsiz, memsiz, mmapsiz;

  size_t start,end;/* index range in alignment[start .. end-1] that uses this memory block */

  Calign *block;/* block[0..numaligns*(colors==1 ? 1 : colors+1)-1] */
  int *intblock;/* enough memory for sites1[],sites2[] (and possibly sitesK1[]) */
  FLOAT *scoreblock;/* If != 0 : enough memory for iscore[],outscore[] */
  Calign **Malignblock;/* If != 0 : enough memory for Malign[] pointer array */
  char *swapfile;/* If != 0 : file that contains copy of ((char *)block)[0..memsiz-1] : If so all block[] etc mmap memory will be in MADV_DONTNEED state */
  FILE *FPswap;
  size_t blen;
  char *write_buffer;

  Calign_block(){
    numaligns = numpairs_sum = 0;
    swapfile = NULL;
    FPswap = NULL;
    write_buffer = NULL;
  };

  Calign_block(size_t num, size_t pairs_sum, size_t Msiz, int myAlignScore, int sitesK, int block_id){
    extern double mtime(),wtime();

    numaligns = num;
    numpairs_sum = pairs_sum;
    size_t clen = (colors==1) ? 1 : colors + 1;
    alignsiz = num * clen;
    intsiz = numpairs_sum * (2 + (sitesK ? 1 : 0));
    intsiz = (intsiz + 3ul) & ~3ul;/* round to nearest multiple of 4 ints == 16 bytes */
    scoresiz = (numpairs_sum + num) * (myAlignScore ? 1 : 0) * 2;
    Malignsiz = Msiz;

    size_t siz = alignsiz * sizeof(Calign);
    siz = (siz + 15ul) & ~(15ul);/* round to nearest multiple of 16 bytes */
    if(DEBUG) assert(siz >= alignsiz * sizeof(Calign));

    memsiz = siz + intsiz*sizeof(int) + scoresiz * sizeof(FLOAT) + Malignsiz * sizeof(Calign *);

    if(HUGEPAGE)
      mmapsiz = (memsiz + HUGEPAGE - 1) & ~(HUGEPAGE-1);
    else
      mmapsiz = (memsiz + PAGE - 1) & ~(PAGE-1);

    blen = USE_MIC ? 10*1024*1024 : 2*1024*1024; // If USE_MIC : must match scif_io.h buffer size to avoid triggering bug (failed fwrite())

    if(num <= 0){
      block = NULL;
      intblock = NULL;
      scoreblock = NULL;
      Malignblock = NULL;
      return;
    }

    //    char *p = (char *)malloc(memsiz);
    char *p = (char *) mmap(NULL,mmapsiz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if(p == MAP_FAILED){
      int eno = errno;
      char *err = strerror(eno);

      #pragma omp critical
      {
	printf("Calign_block: mmap of %lu bytes failed: errno=%d: %s\n",mmapsiz,eno,err);
	dumpmemmap();
	fflush(stdout);exit(1);
      }
    }
    if(HUGEPAGE){
      if(madvise(p, mmapsiz, MADV_HUGEPAGE)){
	int eno = errno;
	char *err = strerror(eno);
	printf("madvise(%p,%lu,MADV_HUGEPAGE) failed: errno=%d: %s\n",block,mmapsiz,eno,err);
	fflush(stdout);
      }
    }

    if(DEBUG>=2) assert((((size_t)p) & 0xf) == 0);
 
    char *origp = p;

    block = (Calign *)p;
    p += siz;
    if(DEBUG>=2) assert((((size_t)p) & 0xf) == 0);

    intblock = (int *)p;
    p += intsiz * sizeof(int);
    if(DEBUG>=2 && !((((size_t)p) & 0xf) == 0)){
      printf("Calign_block: num=%lu, numpairs_sum=%lu: alignsiz=%lu,sizeof(Calign)=%lu (siz=%lu), intsiz=%lu, scoresiz=%lu, memsiz=%lu\n",
	     num,numpairs_sum,alignsiz,sizeof(Calign),siz,intsiz,scoresiz,memsiz);
      printf("\t block=%p(origp=%p),intblock=%p(%lu),p=%p(%ld)\n",(void *)block,(void *)origp,(void *)intblock,((char *)intblock) - ((char *)block),
	     (void *)p, ((char *)p) - ((char *)block));
      fflush(stdout);
      assert((((size_t)p) & 0xf) == 0);
    }

    if(myAlignScore){
      scoreblock = (FLOAT *)p;
      p += scoresiz * sizeof(FLOAT);
    } else
      scoreblock = 0;

    if(Msiz){
      Malignblock = (Calign **)p;
      p += Malignsiz * sizeof(Calign *);
    } else
      Malignblock = NULL;
    if(DEBUG) assert(p == ((char *)block) + memsiz);

    if(VERB/* HERE HERE HERE >=2 */){
      printf("Calign_block: num=%lu, numpairs_sum=%lu: alignsiz=%lu,sizeof(Calign)=%lu (siz=%lu), intsiz=%lu, scoresiz=%lu, Malignsiz=%lu, memsiz=%lu\n",
	     num,numpairs_sum,alignsiz,sizeof(Calign),siz,intsiz,scoresiz,Malignsiz, memsiz);
      printf("\t block=%p(origp=%p),intblock=%p(%lu),scoreblock=%p(%lu),Malignblock=%p\n",(void *)block,(void *)origp,(void *)intblock,((char *)intblock) - ((char *)block),
	     (void *)scoreblock, ((char *)scoreblock) - ((char *)block),Malignblock);
      fflush(stdout);
    }

    int numthreads = 1;
#ifdef _OPENMP
    extern int MaxThreads;
    numthreads = omp_get_max_threads();
    numthreads = min(numthreads,MaxThreads);
    numthreads = max(1,min(numthreads, (int)(alignsiz/512)));
#endif

    #pragma omp parallel for num_threads(numthreads) schedule(static,512) if(numthreads > 1)
    for(long long i = 0; i < (long long)alignsiz; i++)
      block[i].block_id = block_id;

    if(VERB>=2){
      printf("Calign_block: allocated and initialized Calign[%lu]:time=%0.6f(wall time=%0.6f)\n", alignsiz, mtime(),wtime());
      fflush(stdout);
    }

    nextint = nextscore = 0; //    nextalign = 0;
    if(VERB>=2){
      printf("Calign_block: Allocated %lu bytes of block memory:time=%0.6f(wall time=%0.6f)\n", memsiz, mtime(),wtime());
      fflush(stdout);
    }
  };
  
  ~Calign_block(){
    if(numaligns > 0){
      //      free((char *)block);
      if(munmap(block, mmapsiz)){
        int eno = errno;
	char *err = strerror(eno);

#pragma omp critical
	{
	  printf("munmap(%p,%lu) failed: errno=%d:%s\n",block, mmapsiz, eno, err);
	  dumpmemmap();
	  fflush(stdout);exit(1);
	}
      }
      numaligns = 0;
    }
  };

  size_t swap(char *output_prefix, int i){
    double stime;
    if(VERB>=2)
      stime = wtime();

    if(numaligns <= 0){
      printf("Calign_block::swap: WARNING: numaligns= %lu, mmapsiz= %lu\n", numaligns, mmapsiz);
      fflush(stdout);
    }
    char filename[PATH_MAX];
    sprintf(filename,"%s_Calign_block%p",output_prefix,this);
    if(DEBUG && strlen(filename) >= PATH_MAX){
      printf("filename %s_Calign_block%p is too long (%lu chars): must not exceed %d chars\n",output_prefix,this,strlen(filename),PATH_MAX-1);
      fflush(stdout);exit(1);
    }
    swapfile = strdup(filename);

    unlink(swapfile);/* in case it is present from previous run */
    if((FPswap = fopen(filename,"w+")) == NULL){
      int eno = errno;
      char *err = strerror(eno);
      #pragma omp critical
      {
	printf("failed to open file %s for writing:error=%d:%s\n",filename,eno,err);
	if(eno==24){
	  printf("List of all open files (lsof)\n");
	  fflush(stdout);

	  char buf[1024];
	  pid_t pid = getpid();
#if USE_MIC
	  sprintf(buf,"ssh -t -t -oStrictHostKeyChecking=no -oBatchMode=yes host /usr/sbin/lsof -p %lld",(long long)pid);
#else
	  sprintf(buf,"LD_LIBRARY_PATH=/usr/lib64 /usr/sbin/lsof -p %lld", (long long)pid);
#endif
	  int err = system(buf);
	  if(err)
	    printf("system(%s) failed: return code = %d\n", buf, err);
	  if(0){
	    printf("Sleeping 600 seconds\n");
	    fflush(stdout);
	    sleep(600);
	  }
	}

	fflush(stdout);exit(1);
      }
    }

    if(!(write_buffer = (char *)malloc(blen))){
      printf("Calign_block::swap: malloc(%lu) failed\n",blen);
      fflush(stdout);exit(1);
    }
    setbuffer(FPswap, write_buffer, blen);

    size_t n = fwrite_LOOP(block, sizeof(char), mmapsiz, FPswap);
    if(n != mmapsiz){
      #pragma omp critical
      {
	printf("Calign_block::swap: fwrite= %lu: write of %lu bytes to %s failed\n", n, mmapsiz, swapfile);
	fflush(stdout);exit(1);
      }
    }
    fflush(FPswap);

#if 1     // enable to avoid too many open file descriptors
    fclose(FPswap);
    FPswap = NULL;
    free(write_buffer);
    write_buffer = NULL;
#endif

    if(madvise(block,mmapsiz, MADV_DONTNEED)){
      int eno = errno;
      char *err = strerror(eno);
      printf("madvise(%p,%lu,MADV_DONTNEED) failed: errno=%d: %s\n",block,mmapsiz,eno,err);
      fflush(stdout);exit(1);
    }
    
    if(VERB>=2){
      double nwt = wtime();
      printf("Calign_block %d/%d: swapped out & freed %0.4f G (%0.2f Mb/sec) from %p to %s : elapsed= %0.6f (wall time= %0.6f)\n",
	     i,num_align_blocks,mmapsiz * 1e-9, mmapsiz * 1e-6 / (nwt - stime), block, swapfile, nwt - stime, nwt);
      fflush(stdout);
    }

    return mmapsiz;
  };

  size_t restore(int i){
    double stime;
    if(VERB>=2)
      stime = wtime();
    if(!swapfile){
      printf("Calign_block: called restore() before calling swap()\n");
      fflush(stdout);exit(1);
    }
    if(FPswap == NULL){/* reopen file */
      if((FPswap = fopen(swapfile,"r+")) == NULL){
        int eno = errno;
	char *err = strerror(eno);
#pragma omp critical
	{
          printf("failed to open file %s for writing:error=%d:%s\n",swapfile,eno,err);
	  if(eno==24){
	    printf("List of all open files (lsof)\n");
	    fflush(stdout);

	    char buf[1024];
	    pid_t pid = getpid();
#if USE_MIC
	    sprintf(buf,"ssh -t -t -oStrictHostKeyChecking=no -oBatchMode=yes host /usr/sbin/lsof -p %lld",(long long)pid);
#else
	    sprintf(buf,"LD_LIBRARY_PATH=/usr/lib64 /usr/sbin/lsof -p %lld", (long long)pid);
#endif
	    int err = system(buf);
	    if(err)
	      printf("system(%s) failed: return code = %d\n", buf, err);
	  }

	  fflush(stdout);exit(1);
        }
      } /* fopen == NULL */
      if(!(write_buffer = (char *)malloc(blen))){
        printf("Calign_block::restorep: malloc(%lu) failed\n",blen);
        fflush(stdout);exit(1);
      }
    } /* FPswap == NULL */

    fseek(FPswap, 0L, SEEK_SET);
    size_t r;
    if((r = fread(block, sizeof(char), mmapsiz, FPswap)) != mmapsiz){
      printf("Calign_block: Failed to read %lu bytes from %s (r=%lu)\n",mmapsiz,swapfile,r);
      fflush(stdout);exit(1);
    }
    fclose(FPswap);
    FPswap = NULL;
    free(write_buffer);
    unlink(swapfile);

    if(VERB>=2){
      double nwt = wtime();
      printf("Calign_block %d/%d: restored %0.4f G (%0.2f Mb/sec) from %s to %p : elapsed= %0.6f (wall time= %0.6f)\n",
	     i,num_align_blocks,mmapsiz * 1e-9, mmapsiz * 1e-6 / (nwt - stime), swapfile, block, nwt - stime, nwt);
      fflush(stdout);
    }

    free(swapfile);
    swapfile = NULL;

    return mmapsiz;
  };

  inline int *int_alloc(size_t siz){
    size_t origNextint = nextint;
    if((nextint += siz) > intsiz){
      printf("Calign_block.int_alloc(%llu) failed: nextint=%llu, intsiz=%llu\n",(unsigned long long)siz,(unsigned long long)origNextint,(unsigned long long)intsiz);
      exit(1);
    }
    return &intblock[origNextint];
  }

  inline FLOAT *FLOAT_alloc(size_t siz){
    size_t origNextscore = nextscore;
    if((nextscore += siz) > scoresiz){
      printf("Calign_block.FLOAT_alloc(%llu) failed: nextscore=%llu, scoresiz=%llu\n",(unsigned long long)siz,(unsigned long long)origNextscore,(unsigned long long)scoresiz);
      exit(1);
    }
    return &scoreblock[origNextscore];
  }
};

#include "ident.h" // globals.h"
static Ident Calign_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Calign.h 10942 2020-05-02 04:56:21Z tanantharaman $");

#endif
