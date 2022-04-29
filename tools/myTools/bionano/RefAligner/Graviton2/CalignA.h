#ifndef CALIGNA_H
#define CALIGNA_H

#define GOODMAPS 1 /* save detailed alignment information only for good alignments (scoring above 3 thresholds) */

#include <stdio.h>
#include <sys/stat.h>

#include "parameters.h"

class CalignA_block;

class CalignA {/* subset of fields in Calign : only those used by Assembler */
 public:

  int mapid1; ///< index into map[] of 1st map (or index into refmap[] of reference map)
  int mapid2; ///< index into map[] of 2nd map (Nanomap) 
  int fileid; ///< index into align_filename[0..num_align_files-1] (or -1 if alignment was not read in from a file)
  double score; ///< score of best alignment 
  double logPV; ///< -log10(Pvalue), Pvalue without Bonferroni correction 

  int numpairs; ///< number of aligned site pairs (>= 0 : 0 means no valid alignment found) 

  int trueTP; // 1 IFF trueoverlap is >= MINOVERLAP * average size of 2 maps, 0 otherwise

  int *sites1; // sites1[i=0..numpairs-1] are the sites in map[mapid1] that aligned 
  int *sites2; // sites2[i=0..numpairs-1] are the sites in map[mapid2] that aligned : map[mapid2] may be in reversed orientation 
  FLOAT *iscore; /**< iscore[i=1..numpairs-1] are the scores for aligned interval sites*[i-1..i] (includes outlier adjustment)
		 * - iscore[0] is the score for the left end of the alignment (includes endoutlier adjustment)
		 * - iscore[numpairs] is the score of the right end of the alignment (includes endoutlier edjustment) */
  FLOAT *outscore; /**< outscore[i=1..numpairs-1] is the value of iscore[i] before outlier adjustment
		      * - if(outscore[i] < iscore[i]) then the aligned interval [i-1..i] is an outlier.
		      * - outscore[0] and outscore[numpairs] correspond to end intervals and currently always match iscore[] */

  double trueoffset; // true offset between maps (NaN if ground truth is not known)
  double trueoverlap; // true overlap between maps (NaN if ground truth is not known) 

  unsigned short allocated_size;
  int block_id;// -1 : each array is allocated seperately, >= 0 : index into alignA_block[0..num_alignA_blocks-1] 

  int orientation; ///<  orientation of 2nd map in best alignment (0 is normal, 1 = reversed) 
  int scaleID; ///< 2nd map was scaled by ScaleFactor[scaleID] (ScaleFactor[0] == 1.0) 
  char Lend,Rend; /**< Specifies end of alignment type:
		   * - 0 : Linear chromosome end (with ends aligned)
		   * - -1 : Normal alignment end (chromosome ends not aligned but NOT a local alignment)
		   * - -2 or -3 : Local alignment end (ends right after last aligned site) : End Outlier */ 

  /* member functions */
  CalignA(){
    init();
  }
  CalignA(CalignA *pinit, int myscore){/* custom copy constructor (if score==1 also copy iscore[] & outscore[]) */
    extern void copy(CalignA *dest, CalignA *pinit, int myscore, int deepcopy);
    copy(this,pinit,myscore,1);
  }

  ~CalignA(){
    allfree();
  }

  inline void init(){
    mapid1 = mapid2 = -1;/* not initialized */
    fileid = -1;
    numpairs = 0;/* no arrays allocated */
    sites1 = NULL;
    sites2 = NULL;
    iscore = NULL;/* not allocated */
    outscore = NULL;/* not allocated */
    scaleID = 0;
    allocated_size = 0;
    logPV = -1.0;
  }

  inline void allfree(){
    if(allocated_size > 0){
      delete [] sites1; sites1 = NULL;
      delete [] sites2; sites2 = NULL;
      delete [] iscore; iscore = NULL;
      delete [] outscore; outscore = NULL;
    }
    numpairs = 0;
    allocated_size = 0;
    /* NOTE : Do NOT initialize mapid1 etc, since allfree() is used to reallocate alignment memory in pairalign(),refalign(), after mapid1,mapid2 & orientation are recorded */
  }
  
  inline void expand_arrays(int num, int myAlignScore) {// used in input_align()
    if(num > MASK(15)){
      printf("CalignA::expand_arrays:num=%d is too large for 15 bit value : must be <= %d\n",num,MASK(15));
      fflush(stdout);exit(1);
    }
    if(num <= allocated_size)
      return;
    if(allocated_size > 0 && block_id < 0){
      delete [] sites1;
      delete [] sites2;
      delete [] iscore;
      delete [] outscore;
    }
    block_id = -1;
    allocated_size = num;
    sites1 = new int[allocated_size];
    sites2 = new int[allocated_size];
    if(myAlignScore){
      iscore = new FLOAT[allocated_size+1];
      outscore = new FLOAT[allocated_size+1];
    } else
      iscore = outscore = NULL;
  }

  int IsRepeatRegion();/* return K > 0 iff the alignment only covers a single repeat region : shifting the alignment by K sites results in similar average sizing errors */
};

extern void maxalignallocNULL(size_t num, CalignA **&alignmentA, size_t &numaligns, size_t &maxaligns);
extern void copy(CalignA *dest, Calign *pinit, int score, int deepcopy);
extern void alignA_blockfree();/* free up all block memory of alignmentA[] array */
extern CalignA_block *alignA_blockalloc(size_t num,size_t sum, size_t Msiz, int AlignScore, int sitesK, int& block_id);/* allocate a new CalignA_block */
extern size_t align_blockcompact(int block_index, size_t start, int *numpairs, size_t *numpairs_cum, int AlignScore, int numthreads, CalignA **alignmentA);/* compact a CalignA_block */
extern void expand_arrays(CalignA *p, int num, int AlignScore, CalignA_block *block, int block_index);

#define MAX_ALIGNBLOCKS MAXFILES // This should exceed the maximum number of input alignment files */

extern int num_alignA_blocks;
extern CalignA_block *alignA_block[MAX_ALIGNBLOCKS];

#include <sys/mman.h>

#include <stdlib.h>

class CalignA_block {
 public:
  size_t numaligns, alignsiz;
  size_t numpairs_sum;/* upper bound on sum of numpairs of block[].numpairs */
  size_t nextint, nextscore;// nextalign;
  size_t intsiz, scoresiz,  Malignsiz, memsiz, mmapsiz;

  size_t start,end;/* index range in alignment[start .. end-1] that uses this memory block */

  CalignA *block;/* block[0..numaligns*(colors==1 ? 1 : colors+1)-1] */
  int *intblock;/* enough memory for sites1[],sites2[] (and possibly sitesK1[]) */
  FLOAT *scoreblock;/* If != 0 : enough memory for iscore[],outscore[] */
  CalignA **Malignblock;
  char *swapfile;/* If != 0 : file that contains copy of ((char *)block)[0..memsiz-1] : If so all block[] etc mmap memory will be in MADV_DONTNEED state */
  FILE *FPswap;
  size_t blen;
  char *write_buffer;

  CalignA_block(){
    numaligns = numpairs_sum = 0;
    swapfile = NULL;
    FPswap = NULL;
    write_buffer = NULL;
  };

  CalignA_block(size_t num, size_t pairs_sum, size_t Msiz, int myAlignScore, int sitesK, short block_id){
    extern double mtime(),wtime();

    numaligns = num;
    numpairs_sum = pairs_sum;
    size_t clen = (colors==1) ? 1 : colors + 1;
    alignsiz = num * clen;
    intsiz = numpairs_sum * (2 + (sitesK ? 1 : 0));
    intsiz = (intsiz + 3ul) & ~3ul;/* round to nearest multiple of 4 ints == 16 bytes */
    scoresiz = (numpairs_sum + num) * (myAlignScore ? 1 : 0) * 2;
    Malignsiz = Msiz;

    size_t siz = alignsiz * sizeof(CalignA);
    siz = (siz + 15ul) & ~(15ul);/* round to nearest multiple of 16 bytes */
    if(DEBUG) assert(siz >= alignsiz * sizeof(CalignA));
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
	printf("CalignA_block: mmap of %lu bytes failed: errno=%d: %s\n",mmapsiz,eno,err);
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

    block = (CalignA *)p;
    p += siz;
    if(DEBUG>=2) assert((((size_t)p) & 0xf) == 0);

    intblock = (int *)p;
    p += intsiz * sizeof(int);
    if(DEBUG>=2 && !((((size_t)p) & 0xf) == 0)){
      printf("CalignA_block: num=%lu, numpairs_sum=%lu: alignsiz=%lu,sizeof(CalignA)=%lu (siz=%lu), intsiz=%lu, scoresiz=%lu, memsiz=%lu\n",
	     num,numpairs_sum,alignsiz,sizeof(CalignA),siz,intsiz,scoresiz,memsiz);
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
      Malignblock = (CalignA **)p;
      p += Malignsiz * sizeof(Calign *);
    } else
      Malignblock = NULL;

    if(DEBUG) assert(p == ((char *)block) + memsiz);

    if(VERB>=2){
      printf("CalignA_block: num=%lu, numpairs_sum=%lu: alignsiz=%lu,sizeof(CalignA)=%lu (siz=%lu), intsiz=%lu, scoresiz=%lu, memsiz=%lu\n",
	     num,numpairs_sum,alignsiz,sizeof(CalignA),siz,intsiz,scoresiz,memsiz);
      printf("\t block=%p(origp=%p),intblock=%p(%lu),scoreblock=%p(%lu)\n",(void *)block,(void *)origp,(void *)intblock,((char *)intblock) - ((char *)block),
	     (void *)scoreblock, ((char *)scoreblock) - ((char *)block));
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
      printf("CalignA_block: allocated and initialized CalignA[%lu]:time=%0.6f(wall time=%0.6f)\n", alignsiz, mtime(),wtime());
      fflush(stdout);
    }

    nextint = nextscore = 0; //    nextalign = 0;
    if(VERB>=2){
      printf("CalignA_block: Allocated %lu bytes of block memory:time=%0.6f(wall time=%0.6f)\n", memsiz, mtime(),wtime());
      fflush(stdout);
    }
  };
  
  ~CalignA_block(){
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

  ssize_t swap(char *output_prefix, int i){
    double stime;
    if(VERB>=2)
      stime = wtime();

    if(numaligns <= 0){
      printf("CalignA_block::swap: WARNING: numaligns= %lu, mmapsiz= %lu\n", numaligns, mmapsiz);
      fflush(stdout);
    }
    char filename[PATH_MAX];
    filename[0] = '\0';
    if(0 && tmpfs){ /* check if tmpfs exists */
      struct stat sb;
      if(stat(tmpfs, &sb) == 0 && S_ISDIR(sb.st_mode)){
	char cwd[BUFSIZ],dir[BUFSIZ];	
	sprintf(dir,"%s%s",tmpfs,getcwd(cwd,BUFSIZ));
	const int err = mkdir(dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(err){
	  printf("mkdir(%s) failed\n",dir);
	  fflush(stdout);
	} else {

	}
      }
    }
    if(!filename[0])
      sprintf(filename,"%s_CalignA_block%p",output_prefix,this);
    if(DEBUG && strlen(filename) >= PATH_MAX){
      printf("filename %s_CalignA_block%p is too long (%lu chars): must not exceed %d chars\n",output_prefix,this,strlen(filename),PATH_MAX-1);
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
      printf("CalignA_block::swap: malloc(%lu) failed\n",blen);
      fflush(stdout);exit(1);
    }
    setbuffer(FPswap, write_buffer, blen);

    size_t n = fwrite_LOOP(block, sizeof(char), mmapsiz, FPswap);
    if(n != mmapsiz){
      #pragma omp critical
      {
	printf("CalignA_block::swap: fwrite= %lu: write of %lu bytes to %s failed\n", n, mmapsiz, swapfile);
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
      printf("CalignA_block %d/%d: swapped out & freed %0.4f G (%0.2f Mb/sec) from %p to %s : elapsed= %0.6f (wall time= %0.6f)\n",
	     i,num_alignA_blocks,mmapsiz * 1e-9, mmapsiz * 1e-6 / (nwt - stime), block, swapfile, nwt - stime, nwt);
      fflush(stdout);
    }

    return mmapsiz;
  };

  size_t restore(int i){
    double stime;
    if(VERB>=2)
      stime = wtime();
    if(!swapfile){
      printf("CalignA_block: called restore() before calling swap()\n");
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
      printf("CalignA_block: Failed to read %lu bytes from %s (r=%lu)\n",mmapsiz,swapfile,r);
      fflush(stdout);exit(1);
    }
    fclose(FPswap);
    free(write_buffer);
    unlink(swapfile);

    if(VERB>=2){
      double nwt = wtime();
      printf("CalignA_block %d/%d: restored %0.4f G (%0.2f Mb/sec) from %s to %p : elapsed= %0.6f (wall time= %0.6f)\n",
	     i,num_alignA_blocks,mmapsiz * 1e-9, mmapsiz * 1e-6 / (nwt - stime), swapfile, block, nwt - stime, nwt);
      fflush(stdout);
    }

    free(swapfile);
    swapfile = NULL;

    return mmapsiz;
  };

  inline int *int_alloc(size_t siz){
    size_t origNextint = nextint;
    if((nextint += siz) > intsiz){
      printf("CalignA_block.int_alloc(%llu) failed: nextint=%llu, intsiz=%llu\n",(unsigned long long)siz,(unsigned long long)origNextint,(unsigned long long)intsiz);
      exit(1);
    }
    return &intblock[origNextint];
  }

  inline FLOAT *FLOAT_alloc(size_t siz){
    size_t origNextscore = nextscore;
    if((nextscore += siz) > scoresiz){
      printf("CalignA_block.FLOAT_alloc(%llu) failed: nextscore=%llu, scoresiz=%llu\n",(unsigned long long)siz,(unsigned long long)origNextscore,(unsigned long long)scoresiz);
      exit(1);
    }
    return &scoreblock[origNextscore];
  }
};

#include "ident.h" // globals.h"
static Ident CalignA_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/CalignA.h 10703 2020-03-04 18:51:07Z tanantharaman $");

#endif
