#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <omp.h>
#include <set>

using namespace std;

#include "globals.h"
#include "parameters.h"
#include "Ccontig.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/input_vfixx.cpp 10982 2020-05-07 19:18:58Z tanantharaman $");

#define MAPS_PER_THREAD 4096 /* how many maps to parse per thread in each multithreaded block */
#define FILE_BUFFER_SIZE (512*1024*1024) /* number of bytes to read into buffer for multithreaded parsing of BNX file */

#define BNX_FASTFILTER (USE_MIC ? 0 : 1) /* filter BNX maps based on -minlen -maxlen -maxsites -minsites and -selectid(f) in parse_lines */

#ifdef WIN32
typedef long long ssize_t;
#endif

int MaxImageFOV = 0;

#ifdef _MSC_VER

#define snprintf c99_snprintf

inline int c99_vsnprintf(char* str, size_t size, const char* format, va_list ap)
{
	int count = -1;

	if (size != 0)
		count = _vsnprintf_s(str, size, _TRUNCATE, format, ap);
	if (count == -1)
		count = _vscprintf(format, ap);

	return count;
}

inline int c99_snprintf(char* str, size_t size, const char* format, ...)
{
	int count;
	va_list ap;

	va_start(ap, format);
	count = c99_vsnprintf(str, size, format, ap);
	va_end(ap);

	return count;
}

#endif // _MSC_VER


extern void deres(Cmap *pmap, double resKB);/* see input.cpp */

/* A uniformly distributed random number on [0,1) */
/* 0.0 <= urandom() < 1.0 */
static double urandom(void)
{
#ifdef WIN32
	register double r = rand();
	return r / RAND_MAX; 
#else
  register double r = random();
  return r / 2147483648.0; /* 2^31 */
#endif
}

#if 0 // now in globals.h
static double Pr(register double y, register double resKB, register double IresSD)
{
  return 0.5*(1.0+erf((y - resKB)*IresSD));
}
#endif

static int doubleInc(double *p1, double *p2)
{
  return (p1[0] > p2[0]) ? 1 : (p1[0] < p2[0]) ? -1 : 0;
}

int PmapLenDec(Cmap **p1, Cmap **p2)
{
  Cmap *pmap1 = *p1;
  Cmap *pmap2 = *p2;
  int M1 = pmap1->numsite[0];
  int M2 = pmap2->numsite[0];
  double len1 = pmap1->site[0][M1 + 1];
  double len2 = pmap2->site[0][M2 + 1];
  
  return (len1 < len2) ? 1 : (len1 > len2) ? -1 : 0;
}

/* input vfixx/cmap/bnx file and return the additional number of Maps appended to Gmap[0..nummaps-1] */

//int BNXheaderChim_warning = 0;

int ChimQuality_warning = 0;
int CmapMask_warning = 0;

int scannumber_warning = 0;
int chimqual_warning = 0;
int neglen_warning = 0;
int stddev_warning = 0;
int RunIndexWarning = 0;/* warn if RunIndex is missing from BNX 1.0 files with multiple Run Data lines (RunIndex defaults to 1) */
int RunDataWarning = 0;/* warn if there are no Run Data lines in BNX 1.0 files */
int MapWtWarning = 0;
int float_warning = 0;

static char buf[LINESIZ];

struct Cline {
  char *buf;
  int linecnt;
  int RunIndex;
};

extern double mtime();
extern double wtime();

/* First backup the original uncorrected map.
   Then apply bpp (PixelLen/origPixelLen), ScanScaling, minSNR, mres, MaxEnd and MaxInterval corrections to all input maps of current input file */

static void mapcorrect(Cmap **Map, int orignummaps, int &nummaps, int &maxmaps, int fileid, int BNX)
{
  int numthreads = 1;
#ifdef _OPENMP
  numthreads = omp_get_max_threads();
  numthreads = min(numthreads,MaxThreads);
  numthreads = max(1,min(numthreads, (nummaps-orignummaps)/1024));
#endif

  int origcolors = colors;
  if(colors==1 && usecolor){
    if(usecolor == 1){
      printf("WARNING: Cannot use -usecolor %d unless all query inputs have 2 colors: ignoring -usecolor\n",usecolor);
      fflush(stdout);
    } else if(usecolor == 2){
      printf("ERROR: Cannot use -usecolor %d unless all query inputs have 2 colors\n",usecolor);
      fflush(stdout);exit(1);
    }

    usecolor = 0;

    //    if(DEBUG) assert(colors > 1 || usecolor==0);
    //    colors = 2;/* make sure orig* are allocated and saved for both colors */
  }

  if(VERB>=2){
    printf("mapcorrect:colors=%d,usecolor=%d,numthreads=%d,nummaps=%d,orignummaps=%d,mres=%0.3f:time=%0.6f(wall time=%0.6f)\n",colors,usecolor,numthreads,nummaps,orignummaps,mres,mtime(),wtime());
    printf("  MapTrueSite=%d,Tracksites=%d,MapSNR=%d,MapIntensity=%d,MapStitched=%d,MapPSFWidth=%d,MapImageLoc=%d\n",
	   MapTrueSite,Tracksites,MapSNR,MapIntensity,MapStitched,MapPSFWidth,MapImageLoc);
    fflush(stdout);
  }

  /* count total memory required for orig*[] arrays */
  ssize_t memsiz = 0;
  size_t *mapoffset = new size_t[nummaps+1];

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {
    size_t mymemsiz = 0;
    
    #pragma omp for schedule(static,512)
    for(int T = orignummaps; T < nummaps; T++){
      size_t siz = 0;
      Cmap *p = Map[T];// WAS12 Gmap[T]
      for(int c = 0; c < colors; c++){
	int numsite = p->numsite[c];
	siz += (numsite + 2) * sizeof(double);

	if(p->siteSD[c]){
	  siz += (numsite + 2)*(sizeof(double) + 2 * sizeof(float));
	  if(p->SNRcnt[c]){
	    siz += (((numsite + 2) + 1) & ~1) * sizeof(int);
	    siz += (numsite + 2) * (2 * sizeof(double) + sizeof(double *));
	    for(int J = 0; J <= numsite + 1; J++)
	      if(p->SNRcnt[c][J] > 0)
		siz += p->SNRcnt[c][J] * sizeof(double);
	  }
	}

	if(p->ChimQuality[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);
	if(p->ChimNorm[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);
	if(p->SegDupL[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);
	if(p->SegDupR[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);
	if(p->FragileEndL[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);
	if(p->FragileEndR[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);

	if(p->OutlierFrac[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);
	if(p->FragSd[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);
	if(p->ExpSd[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);
	if(p->FragCov[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);
	if(p->FragChiSq[c])
	  siz += (((numsite + 2) + 1) & ~1) *sizeof(float);
	
	if(p->Mask[c])
	  siz += (numsite+2) * sizeof(size_t);

	if(BNX){
          if(MapTrueSite || Tracksites)
	    siz += (numsite + 2) * 2 * sizeof(long long);
	  if(MapSNR)
	    siz += (numsite + 2) * sizeof(double);
	  if(MapIntensity)
	    siz += (numsite + 2) * sizeof(double);
	  if(MapStitched)
	    siz += (((numsite + 2) + 1) & ~1) * sizeof(int);
	  if(MapPSFWidth)
	    siz += (numsite + 2) * sizeof(double);
	  if(MapImageLoc){
            siz += (((numsite + 2) + 1) & ~1) * sizeof(int);
	    siz += (numsite + 2) * 2 * sizeof(double);
          }
        }
      }
      mymemsiz += siz;
      mapoffset[T] = siz;

      if(VERB>=2 && T <= 30){
        #pragma omp critical
        { 
	  printf("T=%d:colors=%d,numsite[0]=%d,siz=%lu,mymemsiz=%lu\n",T,colors,p->numsite[0],siz,mymemsiz);
	  for(int c = 0; c < colors; c++){
            if(c > 0)
	      printf("  numsite[%d]=%d\n",c,p->numsite[c]);
	    if(p->siteSD[c] && p->SNRcnt[c]){
	      for(int J = 0; J <= p->numsite[c] + 1; J++)
	        if(p->SNRcnt[c][J])
		  printf("   SNRcnt[%d][%d]=%d\n",c,J,p->SNRcnt[c][J]);
	    }
	  }
	  fflush(stdout);
        }
      }
      if(DEBUG) assert(mymemsiz < 1000LL*1000LL*1000LL*1000LL);
    }

    #pragma omp atomic
    memsiz += mymemsiz;
  }

  if(DEBUG && !(memsiz >= 0)){
    printf("mapcorrect:memsiz=%ld\n",memsiz);
    fflush(stdout);
    assert(memsiz >= 0);
  }
  if(VERB>=2){
    printf("mapcorrect:memsiz=%ld bytes\n",memsiz);
    fflush(stdout);
  }

  /* allocate a new block */
  maxblockalloc(num_blocks);
  Cmap_blocks[num_blocks] = NULL;
  if(VERB>=2 || LEAKDEBUG>=2){
    printf("mapcorrect:Cmap_blocks[%d] -> 0\n", num_blocks);
    fflush(stdout);
  }
  id_blocks[num_blocks] = NULL;
  mapcnt_blocks[num_blocks] = 0;
  char *mem = site_blocks[num_blocks] = (char *)malloc(memsiz);
  if(!mem){
    printf("mapcorrect:malloc(%llu) failed\n", (unsigned long long)memsiz);
    fflush(stdout);exit(1);
  }
  if(VERB>=2 || LEAKDEBUG>=2){
    printf("mapcorrect: allocated origsite[] map block %d with %llu bytes(site_blocks[%d]=%p):time=%0.6f(wall=%0.6f)\n",
	   num_blocks,(unsigned long long)memsiz, num_blocks, mem, mtime(), wtime());
    if(num_blocks > 1)
      printf("site_blocks[1]=%p\n",site_blocks[1]);
    fflush(stdout);
  }
  num_blocks++;

  /* compute start of each map's memory as offset into mem[] : prefix sum of mapoffset[orignummaps ... nummaps-1] */
  /* NOTE : multi-threading a prefix sum is tricky */
  ssize_t prefix_sum = 0;
  for(int T = orignummaps; T < nummaps; T++){
    ssize_t siz = mapoffset[T];
    mapoffset[T] = prefix_sum;
    prefix_sum += siz;
  }    
  mapoffset[nummaps] = prefix_sum;
  if(DEBUG) assert(prefix_sum == memsiz);

  if(VERB>=2){
    printf("mapcorrect: Computed %d map offsets:time=%0.6f(wall=%0.6f)\n",nummaps-orignummaps, mtime(), wtime());
    if(minSNRestimate)
      printf("minSNR[0]=%0.3f,minSNR[1]=%0.3f\n",minSNR[0],minSNR[1]);
    fflush(stdout);
  }

  #pragma omp parallel for num_threads(numthreads) schedule(dynamic,512)
  for(int T = orignummaps; T < nummaps; T++){
    Cmap *pmap = Map[T];
    char *pmem = &mem[mapoffset[T]];
    char *origpmem = pmem;

    /* initialize and copy original site locations (before resolution based condensation) */
    for(int c=0; c < colors; c++){
      int numsite = pmap->orignumsite[c] = pmap->numsite[c];

      pmap->origsite[c] = (double *)pmem;
      pmem += sizeof(double) * (numsite + 2);
      for(int J=0; J <= numsite+1; J++)
	pmap->origsite[c][J] = pmap->site[c][J];

      if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
	printf("mapcorrect:saving original site locations for id=%lld,c=%d,orignumsite[c] -> %d\n", pmap->id,c,numsite);
	for(int J = 1; J <= numsite+1; J++)
	  if(MapSNR)
	    printf(" origsite[c][%d] -> %0.4f, SNR = %0.3f\n",J,pmap->origsite[c][J], pmap->SNR[c][J]);
	  else
	    printf(" origsite[c][%d] -> %0.4f\n",J,pmap->origsite[c][J]);
	fflush(stdout);
      }

      if(pmap->siteSD[c]){
	pmap->origsiteSD[c] = (double *)pmem;
	pmem += sizeof(double) * (numsite + 2);

	pmap->origsitecov[c] = (float *)pmem;
	pmem += sizeof(float) * (numsite + 2);
	  
	pmap->origsitecnt[c] = (float *)pmem;
	pmem += sizeof(float) * (numsite + 2);

	for(int J=1; J <= numsite; J++){
	  pmap->origsiteSD[c][J] = pmap->siteSD[c][J];
	  pmap->origsitecov[c][J] = pmap->sitecov[c][J];
	  pmap->origsitecnt[c][J] = pmap->sitecnt[c][J];
	}

	if(pmap->SNRcnt[c]){
	  pmap->origSNRcnt[c] = (int *)pmem;
	  pmem += sizeof(int) * (((numsite + 2) + 1) & ~1);
	  pmap->origSNRgmean[c] = (double *)pmem;
	  pmem += sizeof(double) * (numsite + 2);
	  pmap->origlnSNRsd[c] = (double *)pmem;
	  pmem += sizeof(double) * (numsite + 2);
	  pmap->origSNRdist[c] = (double **)pmem;
	  pmem += sizeof(double*) * (numsite + 2);

	  for(int J = 0; J <= pmap->numsite[c]+1; J++){
	    pmap->origSNRcnt[c][J] = pmap->SNRcnt[c][J];
	    pmap->origSNRgmean[c][J] = pmap->SNRgmean[c][J];
	    pmap->origlnSNRsd[c][J] = pmap->lnSNRsd[c][J];
	    pmap->origSNRdist[c][J] = NULL;
	    if(pmap->SNRcnt[c][J] > 0){
	      pmap->origSNRdist[c][J] = (double *)pmem;
	      pmem += sizeof(double) * pmap->SNRcnt[c][J];
	      if(DEBUG>=2) assert(pmem <= &mem[memsiz]);
	      for(int t = 0; t < pmap->SNRcnt[c][J]; t++)
		pmap->origSNRdist[c][J][t] = pmap->SNRdist[c][J][t];
	    }
          }
	}
      }
      if(pmap->ChimQuality[c]){
	pmap->origChimQuality[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++){
	  pmap->origChimQuality[c][J] = pmap->ChimQuality[c][J];
	  if(DEBUG/* HERE >=2 */) assert(pmap->origChimQuality[c][J] >= -1.0);
	}
      }
      if(pmap->ChimNorm[c]){
	pmap->origChimNorm[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origChimNorm[c][J] = pmap->ChimNorm[c][J];
      }
      if(pmap->SegDupL[c]){
	pmap->origSegDupL[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origSegDupL[c][J] = pmap->SegDupL[c][J];
      }
      if(pmap->SegDupR[c]){
	pmap->origSegDupR[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origSegDupR[c][J] = pmap->SegDupR[c][J];
      }
      if(pmap->FragileEndL[c]){
	pmap->origFragileEndL[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origFragileEndL[c][J] = pmap->FragileEndL[c][J];
      }
      if(pmap->FragileEndR[c]){
	pmap->origFragileEndR[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origFragileEndR[c][J] = pmap->FragileEndR[c][J];
      }
      if(pmap->OutlierFrac[c]){
	pmap->origOutlierFrac[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origOutlierFrac[c][J] = pmap->OutlierFrac[c][J];
      }

      if(pmap->FragSd[c]){
	pmap->origFragSd[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origFragSd[c][J] = pmap->FragSd[c][J];
      }
      if(pmap->ExpSd[c]){
	pmap->origExpSd[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origExpSd[c][J] = pmap->ExpSd[c][J];
      }
      if(pmap->FragCov[c]){
	pmap->origFragCov[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origFragCov[c][J] = pmap->FragCov[c][J];
      }
      if(pmap->FragChiSq[c]){
	pmap->origFragChiSq[c] = (float *)pmem;
	pmem += sizeof(float) * (((numsite + 2) + 1) & ~1);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origFragChiSq[c][J] = pmap->FragChiSq[c][J];
      }

      if(pmap->Mask[c]){
	pmap->origMask[c] = (size_t *)pmem;
	pmem += sizeof(size_t) * (numsite + 2);
	for(int J = 0; J <= numsite+1; J++)
	  pmap->origMask[c][J] = pmap->Mask[c][J];
      }

      if(BNX){
	if(MapTrueSite || Tracksites){
	  pmap->origtruesiteL[c] = (long long *)pmem;
	  pmem += sizeof(long long) * (numsite+2);
	  pmap->origtruesiteR[c] = (long long *)pmem;
	  pmem += sizeof(long long) * (numsite+2);
	  for(int J = 1; J <= numsite; J++){
	    pmap->origtruesiteL[c][J] = pmap->truesiteL[c][J];
	    pmap->origtruesiteR[c][J] = pmap->truesiteR[c][J];
	  }
	}
	if(MapSNR){
	  pmap->origSNR[c] = (double *)pmem;
	  pmem += sizeof(double) * (numsite+2);
	  for(int J = 1; J <= numsite; J++)
	    pmap->origSNR[c][J] = pmap->SNR[c][J];
	}
	if(MapIntensity){
	  pmap->origIntensity[c] = (double *)pmem;
	  pmem += sizeof(double) * (numsite+2);
	  for(int J = 1; J <= numsite; J++)
	    pmap->origIntensity[c][J] = pmap->Intensity[c][J];
	}
	if(MapStitched){
	  pmap->origstitch[c] = (int *)pmem;
	  pmem += sizeof(int) * (((numsite+2) + 1) & ~1);
	  for(int J = 1; J <= numsite; J++)
	    pmap->origstitch[c][J] = pmap->stitch[c][J];
	}
	if(MapPSFWidth){
	  pmap->origPSFWidth[c] = (double *)pmem;
	  pmem += sizeof(double) * (numsite+2);
	  for(int J = 1; J <= numsite; J++)
	    pmap->origPSFWidth[c][J] = pmap->PSFWidth[c][J];
	}
	if(MapImageLoc){
	  pmap->origImageFOV[c] = (int *)pmem;
	  pmem += sizeof(int) * (((numsite+2) + 1) & ~1);
	  for(int J = 1; J <= numsite; J++){
	    if(DEBUG>=2) assert(pmap->ImageFOV[c] && 0 <= pmap->ImageFOV[c][J] && pmap->ImageFOV[c][J] <= MaxImageFOV);
	    pmap->origImageFOV[c][J] = pmap->ImageFOV[c][J];
	  }
	  
	  pmap->origImageX[c] = (double *)pmem;
	  pmem += sizeof(double) * (numsite + 2);
	  for(int J = 1; J <= numsite; J++)
	    pmap->origImageX[c][J] = pmap->ImageX[c][J];
	  
	  pmap->origImageY[c] = (double *)pmem;
	  pmem += sizeof(double) * (numsite + 2);
	  for(int J = 1; J <= numsite; J++)
	    pmap->origImageY[c][J] = pmap->ImageY[c][J];
	}
      }
    }
    if(DEBUG) assert(pmem <= &mem[memsiz]);
    if(DEBUG) assert(pmem == &mem[mapoffset[T+1]]);
     if(VERB>=2 && T <= 30 && (colors >= 2 || usecolor)){
       printf("T=%d:pmem=%p .. %p (new site memory interval),pmap->orignumsite[]=%d,%d,pmap->origsite[0]=%p..%p,pmap->origsite[1]=%p..%p,pmap->origImageY[0]=%p..%p,pmap->origImageY[1]=%p..%p\n",
	      T,origpmem,pmem,pmap->orignumsite[0],pmap->orignumsite[1],pmap->origsite[0],&pmap->origsite[0][pmap->orignumsite[0]+1],pmap->origsite[1],&pmap->origsite[1][pmap->orignumsite[1]+1],
	      pmap->origImageY[0],&pmap->origImageY[0][pmap->orignumsite[0]+1],pmap->origImageY[1],&pmap->origImageY[1][pmap->orignumsite[1]+1]);
      fflush(stdout);
    }
  }

  if(VERB>=2){
    printf("mapcorrect: Backed up sites and siteSD (nummaps=%d):time=%0.6f(wall=%0.6f)\n", nummaps, mtime(), wtime());
    fflush(stdout);
  }

  if(DEBUG>=3 && MapSNR){/* print map id BIAS_TRACE_ID */
    for(int m = orignummaps; m < nummaps; m++){
      Cmap *pmap = Map[m];
      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) ||
	 !((pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3) ||
	 !((pmap->origsite[0][pmap->orignumsite[0]+1] - pmap->origsite[1][pmap->orignumsite[1]+1]) < 1e-3)){
	for(int c = 0; c < colors; c++){
	  int M = pmap->numsite[c];
	  printf("mapcorrect:mapid=%d,id=%lld,c=%d,M=%d,origM=%d:\n",m,pmap->id,c,M,pmap->orignumsite[c]);
	  double *X = pmap->site[c];
	  for(int J = 0; J <= M+1; J++)
	    printf(" J=%d:site[c][J]=%0.4f,origsite[c][J]=%0.4f,SNR[c][J]=%0.3f\n",J,X[J],pmap->origsite[c][J],(1 <= J && J <= M) ? pmap->SNR[c][J] : -1.0);
	  fflush(stdout);
	}
	fflush(stdout);
	assert((pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3);
	assert((pmap->origsite[0][pmap->orignumsite[0]+1] - pmap->origsite[1][pmap->orignumsite[1]+1]) < 1e-3);
      }
    }
  }

  delete [] mapoffset;

  colors = origcolors;

  mapcorrectApply(Map, orignummaps, nummaps, maxmaps, BNX, 0, numthreads);
}

static int siteblock = -1;/* id of block memory used for restored site[] etc */

/* apply -breakChiSq to a set of CMAP input maps */
void mapbreakApply(Cmap **&Map, int &nummaps, int &maxmaps)
{
  int numthreads = 1;

#ifdef _OPENMP
  numthreads = omp_get_max_threads();
  numthreads = min(numthreads,MaxThreads);
  numthreads = max(1,min(numthreads, nummaps/64));
#endif

  int startnummaps = nummaps;

  double TotLen = 0.0;
  size_t TotSites = 0;
  long long maxID = 0;

  for(int T = 0; T < startnummaps; T++){
    Cmap *p = Map[T];
    int M = p->numsite[0];
    FLOAT *X = p->site[0];

    maxID = max(maxID, p->id);
    TotLen += X[M+1];
    TotSites += M;
  }
  double AvgInterval = TotLen / max(1.0,(double)TotSites);

  int breakcnt = 0;

  #pragma omp parallel num_threads(numthreads)
  {
    int my_breakcnt = 0;

    #pragma omp for schedule(dynamic,64)
    for(int T = 0; T < startnummaps; T++){
      Cmap *p = Map[T];
      FLOAT *X = p->site[0];
      int M = p->numsite[0], J = 1;

      float *OutlierFrac = p->OutlierFrac[0];
      float *FragSd = p->FragSd[0];
      float *ExpSd = p->ExpSd[0];
      float *FragCov = p->FragCov[0];
      float *FragChiSq = p->FragChiSq[0];
      if(DEBUG/* HERE HERE >=2 */) assert(OutlierFrac != NULL && FragChiSq != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL);

      float *ChimQuality = p->ChimQuality[0];
      float *ChimNorm = p->ChimNorm[0];
      float *SegDupL = p->SegDupL[0];
      float *SegDupR = p->SegDupR[0];
      float *FragileEndL = p->FragileEndL[0];
      float *FragileEndR = p->FragileEndR[0];

      int origbreakcnt = my_breakcnt;

      if(DEBUG && NoBreak)
	for(int J = 1; J < M; J++)
	  assert(0.0 <= FragChiSq[J] && FragChiSq[J] <= 1.0);

      for(; J < M; J++){
        double Range = X[J+1] - X[J];
	double RangeRatio = Range / AvgInterval;
	float PeakLeft = FragCov[J];
	float PeakRight = FragCov[J];
	if(BreakMinCovRatio > 1.0){
          for(int K = J; --K >= 1; )// HERE HERE HERE : terminate early
	    PeakLeft = max(FragCov[K], PeakLeft);
	  for(int K = J; ++K < M; )// HERE HERE HERE : terminate early
	    PeakRight = max(FragCov[K], PeakRight);
	}

	if(DEBUG && !NoBreak) assert(0.0 <= FragChiSq[J] && FragChiSq[J] <= 1.0);

	//	if(FragChiSq[J] < 0.0)
	//	  continue;// already suppressed from neighboring interval

	double FragCovJ = FragCov[J];
	if(BreakChimNorm > 0.0 && ChimNorm[J] > 0.0 && ChimNorm[J+1] > 0.0)
	  FragCovJ = min(FragCovJ, min(ChimNorm[J]*ChimQuality[J],ChimNorm[J+1]*ChimQuality[J+1]) * BreakChimNorm * 0.01);

	if(Range < BreakMaxSize && Range > BreakMinInterval && 
	   (OutlierFrac[J] > BreakFrac ||
	    max(FragileEndL[J]+FragileEndR[J],FragileEndL[J+1]+FragileEndR[J+1]) > BreakFragile ||
	    max(SegDupL[J]+SegDupR[J],SegDupL[J+1]+SegDupR[J+1]) > BreakSegDup ||
	    (FragCovJ < BreakMaxCov && 
	     ((0 <= FragChiSq[J] && FragChiSq[J] < BreakChiSq &&
	       FragSd[J] > ExpSd[J] + BreakMinSf + BreakMinSr * Range &&
	       FragSd[J] > ExpSd[J] * BreakSdRatio) ||		     
	      (FragCovJ < BreakMinCov + RangeRatio *  BreakMinCovScale &&
	       (NoBreak || BreakMinCovRatio <= 1.0 || min(PeakLeft,PeakRight) >= FragCovJ * BreakMinCovRatio)))))){

	  if(!NoBreak)
	    break;

	  int Jmax = min(M-1, J + (NoBreak-1));
	  int Jmin = max(1, J - (NoBreak-1));

	  /*	  if(BreakMinCovRatio > 10 && min(PeakLeft,PeakRight) < FragCovJ * BreakMinCovRatio)
		  Jmin = Jmax = J;*/

	  if(VERB/* HERE HERE >=2 */){
	    #pragma omp critical
	    {
	      double Range = X[J+1] - X[J];
	      double RangeRatio = Range / AvgInterval;
	      printf("Suppressing Deletions in Map[%d] (id=%lld) at Label intervals %d..%d\n", T,Map[T]->id, Jmin,Jmax+1);
	      printf("\t Range= %0.4f(Ratio= %0.3f), OutlierFrac= %0.1f, Fragile= %0.1f,%0.1f, SegDup=%0.1f,%0.1f,FragChiSq= %0.4e, FragCov= %0.2f, ChimNorm=%0.2f,%0.2f, FragSd= %0.1f, ExpSd= %0.1f, Peak=%0.2f,%0.2f\n",
		     Range,RangeRatio,OutlierFrac[J], FragileEndL[J]+FragileEndR[J], FragileEndL[J+1]+FragileEndR[J+1],SegDupL[J]+SegDupR[J],SegDupL[J+1]+SegDupR[J+1],FragChiSq[J],
		     FragCov[J], ChimNorm[J]*ChimQuality[J]*0.01,ChimNorm[J+1]*ChimQuality[J+1]*0.01, FragSd[J],ExpSd[J],PeakLeft,PeakRight);
	      printf("\t    MaxSize= %0.1f, MaxOutlierFrac= %0.1f, MaxCov= %0.2f, ChiSq= %0.4e, MinSf= %0.4f, MinSr= %0.4f, SdRatio= %0.4f, MinInterval= %0.4f, MinCov= %0.2f, MinCovR= %0.2f, Ctrim=%0.2f\n",
                   BreakMaxSize, BreakFrac,BreakMaxCov,BreakChiSq,BreakMinSf,BreakMinSr,BreakSdRatio,BreakMinInterval, BreakMinCov, BreakMinCovScale,BreakMinCovRatio);
	    }
	  }

	  for(int t = Jmin; t <= Jmax; t++){
            if(FragChiSq[t] >= 0.0){
              if(DEBUG>=2) assert(0.0 <= FragChiSq[t] && FragChiSq[t] <= 1.0);
	      FragChiSq[t] = -1.0 - FragChiSq[t];
            }
	    if(DEBUG) assert(FragChiSq[t] < 0.0);
          }
	} else if(VERB>=2 && Map[T]->id == 9347){
	  int Jmax = min(M-1, J + (NoBreak-1));
	  int Jmin = max(1, J - (NoBreak-1));

	    #pragma omp critical
	    {
	      double Range = X[J+1] - X[J];
	      double RangeRatio = Range / AvgInterval;
	      printf("NOT Suppressing Deletions in Map[%d] (id=%lld) at interval J=%d (intervals %d..%d)\n", T,Map[T]->id, J, Jmin,Jmax+1);
	      printf("\t Range= %0.4f(Ratio= %0.3f), OutlierFrac= %0.1f, FragChiSq= %0.4e, FragCov= %0.2f, FragSd= %0.1f, ExpSd= %0.1f, Peak=%0.2f,%0.2f\n",
		     Range,RangeRatio,OutlierFrac[J], FragChiSq[J], FragCov[J], FragSd[J],ExpSd[J],PeakLeft,PeakRight);
	      printf("\t    MaxSize= %0.1f, MaxOutlierFrac= %0.1f, MaxCov= %0.2f, ChiSq= %0.4e, MinSf= %0.4f, MinSr= %0.4f, SdRatio= %0.4f, MinInterval= %0.4f, MinCov= %0.2f, MinCovR= %0.2f, Ctrim=%0.2f\n",
                      BreakMaxSize, BreakFrac,BreakMaxCov,BreakChiSq,BreakMinSf,BreakMinSr,BreakSdRatio,BreakMinInterval, BreakMinCov, BreakMinCovScale,BreakMinCovRatio);
	    }
	}
      }
      
      if(J < M){/* Cmap needs to be broken into multiple pieces */
        int next;

        #pragma omp critical(nummaps)
	{
          maxmapalloc(nummaps+1,maxmaps,Map,0,1);// NOTE : uses "#pragma omp critical(MAX_BLOCKS)"
	  next = nummaps++;
        }
	my_breakcnt++;

	// right part will be saved in Map[next], which will be rechecked later in case it needs to be broken again

	if(VERB>=2){
	  #pragma omp critical
	  {
            double Range = X[J+1] - X[J];
	    double RangeRatio = Range / AvgInterval;

	    printf("splitting Map[%d] (id=%lld) into 2 pieces at Label interval %d..%d, with right piece replacing Map[%d],id=%lld (maxID=%lld)\n",
	         T,Map[T]->id, J,J+1, next,Map[T]->id + maxID * (origbreakcnt - my_breakcnt), maxID);
	    printf("\t Range= %0.4f(Ratio= %0.3f), OutlierFrac= %0.1f, FragChiSq= %0.4e, FragCov= %0.2f, FragSd= %0.1f, ExpSd= %0.1f\n",
		       Range,RangeRatio,OutlierFrac[J], FragChiSq[J], FragCov[J], FragSd[J],ExpSd[J]);
	    printf("\t    MaxSize= %0.1f, MaxOutlierFrac= %0.1f, MaxCov= %0.2f, ChiSq= %0.4e, MinSf= %0.4f, MinSr= %0.4f, SdRatio= %0.4f, MinInterval= %0.4f, MinCov= %0.2f, MinCovR= %0.2f\n",
		       BreakMaxSize, BreakFrac,BreakMaxCov,BreakChiSq,BreakMinSf,BreakMinSr,BreakSdRatio,BreakMinInterval, BreakMinCov, BreakMinCovScale);
	    fflush(stdout);
	  }
	}

	long long pid = p->id;
	int fileid = p->fileid;

	if(1){
	  /* compute right part */
	  Cmap *nrmap = new Cmap(p);
	  nrmap->trim(X[J] + 0.001, X[M+1]);
	  Cmap *newrmap = new Cmap(nrmap);// compacts the memory for the trimmed map
	  newrmap->id = pid + maxID * (origbreakcnt - my_breakcnt);
	  newrmap->mapid = next;
	  newrmap->fileid = fileid;
	  delete nrmap;// frees memory of original untrimmed map
	  
	  /* swap *Map[next] with *newrmap  and free sub-array memory of original *Map[next] */
	  Cmap tmp;
	  tmp = *Map[next];
	  *Map[next] = *newrmap;
	  tmp.allfree();// NOTE : this frees sub-array memory of original *Map[next] 
	  *newrmap = tmp;
	  delete newrmap;// just frees the Cmap struct
	}

	if(1){
	  /* compute left part of p == Map[T] */
	  if(DEBUG) assert(p == Map[T]);
	  p->trim(0.0, X[J+1] - 0.001);
	  Cmap *newrmapL = new Cmap(p);// compacts the mmeory for the trimmed map
	  newrmapL->id = pid;
	  newrmapL->fileid = fileid;
	  newrmapL->mapid = T;
	  
	  /* swap *Map[T] with *newrmap and free sub-array memory of original *Map[T] */
	  Cmap tmpL = *Map[T];
	  *Map[T] = *newrmapL;
	  tmpL.allfree();// NOTE : this frees sub-array memory of original *Map[T]
	  *newrmapL = tmpL;
	  delete newrmapL;// just frees the Cmap struct
	}
      } // if (J < M)
    }// omp for(int T = 0; T < numstartmaps; T++)

    #pragma omp atomic
    breakcnt += my_breakcnt;
  } // #pragma omp parallel num_threads(numthreads)

  if(nummaps > startnummaps){// check if new fragments need to be broken again 
    for(int T = startnummaps; T < nummaps; T++){
      Cmap *p = Map[T];
      FLOAT *X = p->site[0];
      int M = p->numsite[0], J = 1;

      float *OutlierFrac = p->OutlierFrac[0];
      float *FragSd = p->FragSd[0];
      float *ExpSd = p->ExpSd[0];
      float *FragCov = p->FragCov[0];
      float *FragChiSq = p->FragChiSq[0];
      if(DEBUG) assert(OutlierFrac != NULL && FragChiSq != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL);

      for(; J < M; J++){
        double Range = X[J+1] - X[J];
	double RangeRatio = Range / AvgInterval;

	if(Range < BreakMaxSize &&
	   (OutlierFrac[J] > BreakFrac ||
	   (FragCov[J] < BreakMaxCov &&
	    Range > BreakMinInterval &&
	    ((FragChiSq[J] < BreakChiSq &&
	      FragSd[J] > ExpSd[J] + BreakMinSf + BreakMinSr * Range &&
	      FragSd[J] > ExpSd[J] * BreakSdRatio) ||		     
	     FragCov[J] < BreakMinCov + RangeRatio *  BreakMinCovScale))))
	  break;
      }
      
      if(J < M){/* Cmap needs to be broken into multiple pieces */
        int next;

	//        #pragma omp critical(nummaps)
	{
          maxmapalloc(nummaps+1,maxmaps,Map,0,1);// NOTE : uses "#pragma omp critical(MAX_BLOCKS)"
	  next = nummaps++;
        }

	// right part will be saved in Map[next], which will be rechecked later in case it needs to be broken again

	if(VERB>=2){
	  //	  #pragma omp critical
	  {
            double Range = X[J+1] - X[J];
	    double RangeRatio = Range / AvgInterval;

	    printf("splitting Map[%d] (id=%lld) into 2 pieces at Label interval %d..%d, with right piece replacing Map[%d], maxID=%lld\n",
		   T,Map[T]->id, J,J+1, next,maxID);
	    printf("\t Range= %0.4f(Ratio= %0.3f), OutlierFrac= %0.1f, FragChiSq= %0.4e, FragCov= %0.2f, FragSd= %0.1f, ExpSd= %0.1f\n",
		       Range,RangeRatio,OutlierFrac[J], FragChiSq[J], FragCov[J], FragSd[J],ExpSd[J]);
	    printf("\t    MaxSize= %0.1f, MaxOutlierFrac= %0.1f, MaxCov= %0.2f, ChiSq= %0.4e, MinSf= %0.4f, MinSr= %0.4f, SdRatio= %0.4f, MinInterval= %0.4f, MinCov= %0.2f, MinCovR= %0.2f\n",
		       BreakMaxSize, BreakFrac,BreakMaxCov,BreakChiSq,BreakMinSf,BreakMinSr,BreakSdRatio,BreakMinInterval, BreakMinCov, BreakMinCovScale);
	    fflush(stdout);
	  }
	}
	breakcnt++;

	long long pid = p->id;
	int fileid = p->fileid;

	if(1){
	  /* compute right part */
	  Cmap *nrmap = new Cmap(p);
	  nrmap->trim(X[J] + 0.001, X[M+1]);
	  Cmap *newrmap = new Cmap(nrmap);// compacts the memory for the trimmed map
	  newrmap->id = pid + maxID;
	  newrmap->mapid = next;
	  newrmap->fileid = fileid;
	  delete nrmap;// frees memory of original untrimmed map
	  
	  /* swap *Map[next] with *newrmap  and free sub-array memory of original *Map[next] */
	  Cmap tmp;
	  tmp = *Map[next];
	  *Map[next] = *newrmap;
	  tmp.allfree();// NOTE : this frees sub-array memory of original *Map[next] 
	  *newrmap = tmp;
	  delete newrmap;// just frees the Cmap struct
	}

	if(1){
	  /* compute left part of p == Map[T] */
	  if(DEBUG) assert(p == Map[T]);
	  p->trim(0.0, X[J+1] - 0.001);
	  Cmap *newrmapL = new Cmap(p);// compacts the mmeory for the trimmed map
	  newrmapL->id = pid;
	  newrmapL->fileid = fileid;
	  newrmapL->mapid = T;
	  
	  /* swap *Map[T] with *newrmap and free sub-array memory of original *Map[T] */
	  Cmap tmpL = *Map[T];
	  *Map[T] = *newrmapL;
	  tmpL.allfree();// NOTE : this frees sub-array memory of original *Map[T]
	  *newrmapL = tmpL;
	  delete newrmapL;// just frees the Cmap struct
	}

	p = Map[T];
      } // if (J < M)
    }
  }


}

/* If restore==1, first restore all maps from orig*

   Then apply bpp, ScanScaling[c], minSNR, mres, MaxEnd and MaxInterval corrections to all maps Map[orignummaps .. nummaps-1]

   If ScanScale[c] || MapScale : apply scan scaling factors ScanScale[c][scan].cumscale and per Map scaling factors Map[i]->cumscale to site[] values  */

void mapcorrectApply(Cmap **Map, int orignummaps, int &nummaps, int &maxmaps, int BNX, int restore, int numthreads)
{
  if(VERB>=2){
    printf("mapcorrectApply:orignummaps=%d,nummaps=%d,restore=%d,numthreads=%d\n",orignummaps,nummaps,restore,numthreads);
    printf("  MapTrueSite=%d,Tracksites=%d,MapSNR=%d,MapIntensity=%d,MapStitched=%d,MapPSFWidth=%d,MapImageLoc=%d\n",
	   MapTrueSite,Tracksites,MapSNR,MapIntensity,MapStitched,MapPSFWidth,MapImageLoc);
    fflush(stdout);
  }
  //  if(DEBUG && restore) assert(minSNRestimate);

  if(DEBUG>=3 && (colors >= 2 || usecolor)){
    for(int i = orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];// WAS12 Gmap[i]
      double len1 = pmap->site[0][pmap->numsite[0]+1];
      double len2 = pmap->site[1][pmap->numsite[1]+1];
      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
	printf("mapcorrectApply start : i=%d:colors=%d:pmap->id=%lld:\n",i,colors,pmap->id);
	for(int c = 0; c < max(2,colors); c++)
	  printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.3f\n",
		 c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	fflush(stdout);
      }
      assert(fabs(len1-len2) < 0.001);
      
      len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
      len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
	printf("mapcorrectApply start: i=%d:colors=%d:pmap->id=%lld:\n",i,colors,pmap->id);
	for(int c = 0; c < max(2,colors); c++)
	  printf("  c=%d:pmap->orignumsite[c]=%d,pmap->origsite[c][pmap->orignumsite[c]+1] = %0.3f,pmap->origsite[c]=%p\n",
		 c,pmap->orignumsite[c],pmap->origsite[c][pmap->orignumsite[c]+1],pmap->origsite[c]);
	fflush(stdout);
      }
      assert(fabs(len1-len2) < 0.001);
    }
  }

  int origcolors = colors;
  if(colors==1 && usecolor)
    colors = 2;/* make sure orig* are restored for both colors */

  if(restore){/* restore all maps */

    /* count total memory required for orig*[] arrays */
    size_t memsiz = 0;
    size_t *mapoffset = new size_t[nummaps+1];
    
    #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
    {
      size_t mymemsiz = 0;
    
      #pragma omp for schedule(static,512)
      for(int T = orignummaps; T < nummaps; T++){
        size_t siz = 0;
        Cmap *p = Map[T];// WAS12 Gmap[T]
        for(int c = 0; c < colors; c++){
	  int numsite = p->orignumsite[c];
  	  siz += (numsite + 2) * sizeof(double);

	  if(p->origsiteSD[c]){
	    siz += (numsite + 2)*(sizeof(double) + 2 * sizeof(float));
	    if(p->origSNRcnt[c]){
	      siz += (numsite + 2) * (2 * sizeof(double) + sizeof(double *));
	      siz += sizeof(int) * (((numsite+2) + 1) & ~1);
	      for(int J = 0; J <= numsite + 1; J++)
	        if(p->origSNRcnt[c][J] > 0)
		  siz += p->origSNRcnt[c][J] * sizeof(double);
	    }
	  }
	  if(p->origChimQuality[c])
	    siz += (numsite + 2)*sizeof(float);
	  if(p->ChimNorm[c])
	    siz += (numsite + 2)*sizeof(float);
	  if(p->SegDupL[c])
	    siz += (numsite + 2)*sizeof(float);
	  if(p->SegDupR[c])
	    siz += (numsite + 2)*sizeof(float);
	  if(p->FragileEndL[c])
	    siz += (numsite + 2)*sizeof(float);
	  if(p->FragileEndR[c])
	    siz += (numsite + 2)*sizeof(float);
	  if(p->OutlierFrac[c])
	    siz += (numsite + 2)*sizeof(float);
	
	  if(p->FragSd[c])
	    siz += (numsite + 2)*sizeof(float);
	  if(p->ExpSd[c])
	    siz += (numsite + 2)*sizeof(float);
	  if(p->FragCov[c])
	    siz += (numsite + 2)*sizeof(float);
	  if(p->FragChiSq[c])
	    siz += (numsite + 2)*sizeof(float);

	  if(p->Mask[c])
	    siz += (numsite+2) * sizeof(size_t);
	  
	  if(BNX){
            if(MapTrueSite || Tracksites)
	      siz += (numsite + 2) * 2 * sizeof(long long);
	    if(MapSNR)
	      siz += (numsite + 2) * sizeof(double);
	    if(MapIntensity)
	      siz += (numsite + 2) * sizeof(double);
	    if(MapStitched)
	      siz += sizeof(int) * (((numsite+2) + 1) & ~1);
	    if(MapPSFWidth)
	      siz += (numsite + 2) * sizeof(double);
	    if(MapImageLoc){
              siz += sizeof(int) * (((numsite+2) + 1) & ~1);
	      siz += (numsite + 2) * 2 * sizeof(double);
            }
          }
        }
	mymemsiz += siz;
	mapoffset[T] = siz;
      }

      #pragma omp atomic
      memsiz += mymemsiz;
    }

    /* allocate a new block */
    maxblockalloc(num_blocks);
    Cmap_blocks[num_blocks] = NULL;
    if(VERB>=2 || LEAKDEBUG>=2){
      printf("mapcorrectApply:Cmap_blocks[%d] -> 0\n", num_blocks);
      fflush(stdout);
    }
    id_blocks[num_blocks] = NULL;
    mapcnt_blocks[num_blocks] = 0;
    char *mem = site_blocks[num_blocks] = (char *)malloc(memsiz);
    if(!mem){
      printf("mapcorrectApply:malloc(%llu) failed\n", (unsigned long long)memsiz);
      fflush(stdout);exit(1);
    }
    if(VERB>=2 || LEAKDEBUG>=2){
      printf("mapcorrectApply: allocated map block %d with %llu bytes(site_blocks[%d]=%p):time=%0.6f(wall=%0.6f)\n",
         num_blocks,(unsigned long long)memsiz, num_blocks, mem, mtime(), wtime());
      if(num_blocks > 1)
        printf("site_blocks[1]=%p\n",site_blocks[1]);
      fflush(stdout);
    }

    if(siteblock >= 0){
      if(VERB>=2 || LEAKDEBUG>=2){
        printf("mapcorrectApply: freeing previous map block %d (site_blocks[%d]=%p)\n", siteblock,siteblock,site_blocks[siteblock]);
	fflush(stdout);
      }
      free(site_blocks[siteblock]);
      site_blocks[siteblock] = NULL;
    }
    siteblock = num_blocks++;

    if(DEBUG>=3 && (colors >= 2 || usecolor)){
      for(int i = orignummaps; i < nummaps; i++){
        Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
	double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
	  printf("mapcorrectApply after block allocation : i=%d:colors=%d:pmap->id=%lld:\n",i,colors,pmap->id);
	  for(int c = 0; c < max(2,colors); c++)
	    printf("  c=%d:pmap->orignumsite[c]=%d,pmap->origsite[c][pmap->orignumsite[c]+1] = %0.3f, pmap->origsite[c]=%p\n",
		   c,pmap->orignumsite[c],pmap->origsite[c][pmap->orignumsite[c]+1],pmap->origsite[c]);
	  fflush(stdout);
	}
	assert(fabs(len1-len2) < 0.001);
      }
    }

    /* compute start of each map's memory as offset into mem[] : prefix sum of mapoffset[orignummaps ... nummaps-1] */
    /* NOTE : multi-threading a prefix sum is tricky */
    size_t prefix_sum = 0;
    for(int T = orignummaps; T < nummaps; T++){
      size_t siz = mapoffset[T];
      mapoffset[T] = prefix_sum;
      prefix_sum += siz;
    }    
    mapoffset[nummaps] = prefix_sum;
    if(DEBUG) assert(prefix_sum == memsiz);

    if(VERB>=2){
      printf("mapcorrectApply: Computed %d map offsets:time=%0.6f(wall=%0.6f)\n",nummaps-orignummaps, mtime(), wtime());
      fflush(stdout);
    }

    if(DEBUG>=3 && (colors >= 2 || usecolor)){
      for(int i = orignummaps; i < nummaps; i++){
        Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
	double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
	  printf("mapcorrectApply before site restore : i=%d:colors=%d:pmap->id=%lld:\n",i,colors,pmap->id);
	  for(int c = 0; c < max(2,colors); c++)
	    printf("  c=%d:pmap->orignumsite[c]=%d,pmap->origsite[c][pmap->orignumsite[c]+1] = %0.3f, pmap->origsite[c]=%p\n",
		   c,pmap->orignumsite[c],pmap->origsite[c][pmap->orignumsite[c]+1],pmap->origsite[c]);
	  fflush(stdout);
	}
	assert(fabs(len1-len2) < 0.001);
      }
    }

    #pragma omp parallel for num_threads(numthreads) schedule(dynamic,512)
    for(int T = orignummaps; T < nummaps; T++){
      Cmap *pmap = Map[T];
      char *pmem = &mem[mapoffset[T]];
      if(VERB>=2 && T <= 10 && (colors >= 2 || usecolor)){
        printf("T=%d:pmem=%p .. %p (new site memory interval),pmap->orignumsite[]=%d,%d,pmap->origsite[0]=%p..%p,pmap->origsite[1]=%p..%p,pmap->origImageY[0]=%p..%p,pmap->origImageY[1]=%p..%p\n",
	       T,pmem,&mem[mapoffset[T+1]],pmap->orignumsite[0],pmap->orignumsite[1],pmap->origsite[0],&pmap->origsite[0][pmap->orignumsite[0]+1],pmap->origsite[1],&pmap->origsite[1][pmap->orignumsite[1]+1],
	       pmap->origImageY[0],&pmap->origImageY[0][pmap->orignumsite[0]+1],pmap->origImageY[1],&pmap->origImageY[1][pmap->orignumsite[1]+1]);
        fflush(stdout);
      }

      if(DEBUG>=4 && (colors >= 2 || usecolor)){
        for(int i = orignummaps; i < nummaps; i++){
          Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
	  double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	  double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	  if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
            printf("mapcorrectApply before T=%d site restore (pmem=%p) : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,i,colors,pmap->id);
	    for(int lc = 0; lc < max(2,colors); lc++)
	      printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
                lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
	    fflush(stdout);
          }
          assert(fabs(len1-len2) < 0.001);
        }
      }

      /* initialize and copy original site locations (before resolution based condensation) */
      for(int c=0; c < colors; c++){
	int numsite = pmap->numsite[c] = pmap->orignumsite[c];

	if(!pmap->blockmem)
	  delete [] pmap->site[c];
	pmap->site[c] = (double *)pmem;
	pmem += sizeof(double) * (numsite + 2);

        if(DEBUG>=4 && (colors >= 2 || usecolor)){
          for(int i = orignummaps; i < nummaps; i++){
            Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
  	    double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	    double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	    if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
              printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: A : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
	      for(int lc = 0; lc < max(2,colors); lc++)
	        printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
                  lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
	      fflush(stdout);
            }
            assert(fabs(len1-len2) < 0.001);
          }
        }

	if(VERB>=2 && BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
	  printf("mapcorrectApply:restoring id=%lld,c=%d,numsite[c] -> %d\n", pmap->id,c,numsite);
	  for(int J = 1; J <= numsite+1; J++)
	    if(MapSNR)
	      printf(" site[c][%d] -> %0.4f, SNR= %0.3f\n",J,pmap->origsite[c][J], (J <= numsite) ? pmap->origSNR[c][J] : -1.0 );
	    else
	      printf(" site[c][%d] -> %0.4f\n",J,pmap->origsite[c][J]);
	  fflush(stdout);
        }
	for(int J=0; J <= numsite+1; J++)
	  pmap->site[c][J] = pmap->origsite[c][J];

        if(DEBUG>=4 && (colors >= 2 || usecolor)){
          for(int i = orignummaps; i < nummaps; i++){
            Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
  	    double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	    double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	    if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
              printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: B : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
	      for(int lc = 0; lc < max(2,colors); lc++)
	        printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
                  lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
	      fflush(stdout);
            }
            assert(fabs(len1-len2) < 0.001);
          }
        }

	if(pmap->origsiteSD[c]){
	  if(pmap->blockmem >= 0){
	    delete [] pmap->siteSD[c];
	    delete [] pmap->sitecov[c];
	    delete [] pmap->sitecnt[c];
          }
          if(DEBUG>=4 && (colors >= 2 || usecolor)){
            for(int i = orignummaps; i < nummaps; i++){
              Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
    	      double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	      double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
                printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: C : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
	        for(int lc = 0; lc < max(2,colors); lc++)
	          printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
                    lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
	        fflush(stdout);
              }
              assert(fabs(len1-len2) < 0.001);
            }
          }

    	  pmap->siteSD[c] = (double *)pmem;
	  pmem += sizeof(double) * (numsite + 2);

	  pmap->sitecov[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  
	  pmap->sitecnt[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);

	  for(int J=1; J <= pmap->numsite[c]; J++){
	    pmap->siteSD[c][J] = pmap->origsiteSD[c][J];
	    pmap->sitecov[c][J] = pmap->origsitecov[c][J];
	    pmap->sitecnt[c][J] = pmap->origsitecnt[c][J];
	  }

          if(DEBUG>=4 && (colors >= 2 || usecolor)){
            for(int i = orignummaps; i < nummaps; i++){
              Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
    	      double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	      double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
                printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: D : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
	        for(int lc = 0; lc < max(2,colors); lc++)
	          printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
                    lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
	        fflush(stdout);
              }
              assert(fabs(len1-len2) < 0.001);
            }
          }

	  if(pmap->origSNRcnt[c]){
	    if(pmap->blockmem >= 0){
	      for(int J = 0; J <= numsite + 1;  J++)
		if(pmap->SNRcnt[c][J] > 0)
		  delete [] pmap->SNRdist[c][J];
	      delete [] pmap->SNRcnt[c];
	      delete [] pmap->SNRgmean[c];
	      delete [] pmap->lnSNRsd[c];
	      delete [] pmap->SNRdist[c];
	    }
	    if(DEBUG>=4 && (colors >= 2 || usecolor)){
	      for(int i = orignummaps; i < nummaps; i++){
		Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
		double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
		double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
		if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
		  printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: E : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
		  for(int lc = 0; lc < max(2,colors); lc++)
		    printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
			   lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
		  fflush(stdout);
		}
		assert(fabs(len1-len2) < 0.001);
	      }
	    }

	    pmap->SNRcnt[c] = (int *)pmem;
	    pmem += sizeof(int) * (((numsite + 2) + 1) & ~1);
	    pmap->SNRgmean[c] = (double *)pmem;
	    pmem += sizeof(double) * (numsite + 2);
	    pmap->lnSNRsd[c] = (double *)pmem;
	    pmem += sizeof(double) * (numsite + 2);
	    pmap->SNRdist[c] = (double **)pmem;
	    pmem += sizeof(double*) * (numsite + 2);

	    for(int J = 0; J <= numsite + 1; J++){
	      pmap->SNRcnt[c][J] = pmap->origSNRcnt[c][J];
	      pmap->SNRgmean[c][J] = pmap->origSNRgmean[c][J];
	      pmap->lnSNRsd[c][J] = pmap->origlnSNRsd[c][J];
	      pmap->SNRdist[c][J] = NULL;
	      if(pmap->origSNRcnt[c][J] > 0){
		pmap->SNRdist[c][J] = (double *)pmem;
		pmem += sizeof(double) * pmap->origSNRcnt[c][J];
		if(DEBUG>=2) assert(pmem <= &mem[memsiz]);
		for(int t = 0; t < pmap->origSNRcnt[c][J]; t++)
		  pmap->SNRdist[c][J][t] = pmap->origSNRdist[c][J][t];
	      }
	    }
	    if(DEBUG>=4 && (colors >= 2 || usecolor)){
	      for(int i = orignummaps; i < nummaps; i++){
		Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
		double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
		double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
		if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
		  printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: F : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
		  for(int lc = 0; lc < max(2,colors); lc++)
		    printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
			   lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
		  fflush(stdout);
		}
		assert(fabs(len1-len2) < 0.001);
	      }
	    }
	  }
        }
	if(pmap->origChimQuality[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->ChimQuality[c];
    	  pmap->ChimQuality[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++){
            pmap->ChimQuality[c][J] = pmap->origChimQuality[c][J];
	    if(DEBUG/* HERE >=2 */) assert(pmap->ChimQuality[c][J] >= -1.0);
	  }
	}
	if(pmap->origChimNorm[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->ChimNorm[c];
    	  pmap->ChimNorm[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->ChimNorm[c][J] = pmap->origChimNorm[c][J];
	}
	if(pmap->origSegDupL[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->SegDupL[c];
    	  pmap->SegDupL[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->SegDupL[c][J] = pmap->origSegDupL[c][J];
	}
	if(pmap->origSegDupR[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->SegDupR[c];
    	  pmap->SegDupR[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->SegDupR[c][J] = pmap->origSegDupR[c][J];
	}
	if(pmap->origFragileEndL[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->FragileEndL[c];
    	  pmap->FragileEndL[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->FragileEndL[c][J] = pmap->origFragileEndL[c][J];
	}
	if(pmap->origFragileEndR[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->FragileEndR[c];
    	  pmap->FragileEndR[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->FragileEndR[c][J] = pmap->origFragileEndR[c][J];
	}
	if(pmap->origOutlierFrac[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->OutlierFrac[c];
    	  pmap->OutlierFrac[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->OutlierFrac[c][J] = pmap->origOutlierFrac[c][J];
	}

	if(pmap->origFragSd[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->FragSd[c];
    	  pmap->FragSd[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->FragSd[c][J] = pmap->origFragSd[c][J];
	}
	if(pmap->origExpSd[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->ExpSd[c];
    	  pmap->ExpSd[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->ExpSd[c][J] = pmap->origExpSd[c][J];
	}
	if(pmap->origFragCov[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->FragCov[c];
    	  pmap->FragCov[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->FragCov[c][J] = pmap->origFragCov[c][J];
	}
	if(pmap->origFragChiSq[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->FragChiSq[c];
    	  pmap->FragChiSq[c] = (float *)pmem;
	  pmem += sizeof(float) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->FragChiSq[c][J] = pmap->origFragChiSq[c][J];
	}

	if(pmap->origMask[c]){
	  if(pmap->blockmem >= 0)
	    delete [] pmap->Mask[c];
    	  pmap->Mask[c] = (size_t *)pmem;
	  pmem += sizeof(size_t) * (numsite + 2);
	  for(int J = 0; J <= pmap->numsite[c]+1; J++)
            pmap->Mask[c][J] = pmap->origMask[c][J];
	}

        if(DEBUG>=4 && (colors >= 2 || usecolor)){
          for(int i = orignummaps; i < nummaps; i++){
            Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
  	    double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	    double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	    if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
              printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: G : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
	      for(int lc = 0; lc < max(2,colors); lc++)
	        printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
                  lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
	      fflush(stdout);
            }
            assert(fabs(len1-len2) < 0.001);
          }
        }

	if(BNX){
	  if(MapTrueSite || Tracksites){
	    if(!pmap->blockmem){
	      delete [] pmap->truesiteL[c];
	      delete [] pmap->truesiteR[c];
	    }
	    pmap->truesiteL[c] = (long long *)pmem;
	    pmem += sizeof(long long) * (numsite+2);
	    pmap->truesiteR[c] = (long long *)pmem;
	    pmem += sizeof(long long) * (numsite+2);
	    for(int J = 1; J <= numsite; J++){
	      pmap->truesiteL[c][J] = pmap->origtruesiteL[c][J];
	      pmap->truesiteR[c][J] = pmap->origtruesiteR[c][J];
	    }
	  }
	  if(DEBUG>=4 && (colors >= 2 || usecolor)){
	    for(int i = orignummaps; i < nummaps; i++){
	      Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
	      double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	      double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
		printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: H : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
		for(int lc = 0; lc < max(2,colors); lc++)
		  printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
			 lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
		fflush(stdout);
	      }
	      assert(fabs(len1-len2) < 0.001);
	    }
	  }

	  if(MapSNR){
	    if(!pmap->blockmem)
	      delete [] pmap->SNR[c];
	    pmap->SNR[c] = (double *)pmem;
	    pmem += sizeof(double) * (numsite+2);
	    for(int J = 1; J <= numsite; J++)
	      pmap->SNR[c][J] = pmap->origSNR[c][J];
	  }
	  if(DEBUG>=4 && (colors >= 2 || usecolor)){
	    for(int i = orignummaps; i < nummaps; i++){
	      Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
	      double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	      double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
		printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: I : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
		for(int lc = 0; lc < max(2,colors); lc++)
		  printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
			 lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
		fflush(stdout);
	      }
	      assert(fabs(len1-len2) < 0.001);
	    }
	  }

	  if(MapIntensity){
	    if(!pmap->blockmem)
	      delete [] pmap->Intensity[c];
	    pmap->Intensity[c] = (double *)pmem;
	    pmem += sizeof(double) * (numsite+2);
	    for(int J = 1; J <= numsite; J++)
	      pmap->Intensity[c][J] = pmap->origIntensity[c][J];
	  }
	  if(DEBUG>=4 && (colors >= 2 || usecolor)){
	    for(int i = orignummaps; i < nummaps; i++){
	      Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
	      double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	      double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
		printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: J : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
		for(int lc = 0; lc < max(2,colors); lc++)
		  printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
			 lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
		fflush(stdout);
	      }
	      assert(fabs(len1-len2) < 0.001);
	    }
	  }

	  if(MapStitched){
	    if(!pmap->blockmem)
	      delete [] pmap->stitch[c];
	    pmap->stitch[c] = (int *)pmem;
	    pmem += sizeof(int) * (((numsite + 2) + 1) & ~1);
	    for(int J = 1; J <= numsite; J++)
	      pmap->stitch[c][J] = pmap->origstitch[c][J];
	  }
	  if(DEBUG>=4 && (colors >= 2 || usecolor)){
	    for(int i = orignummaps; i < nummaps; i++){
	      Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
	      double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	      double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
		printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: K : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
		for(int lc = 0; lc < max(2,colors); lc++)
		  printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
			 lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
		fflush(stdout);
	      }
	      assert(fabs(len1-len2) < 0.001);
	    }
	  }

	  if(MapPSFWidth){
	    if(!pmap->blockmem)
	      delete [] pmap->PSFWidth[c];
	    pmap->PSFWidth[c] = (double *)pmem;
	    pmem += sizeof(double) * (numsite+2);
	    for(int J = 1; J <= numsite; J++)
	      pmap->PSFWidth[c][J] = pmap->origPSFWidth[c][J];
	  }
	  if(DEBUG>=4 && (colors >= 2 || usecolor)){
	    for(int i = orignummaps; i < nummaps; i++){
	      Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
	      double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
	      double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
	      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
		printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: L : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
		for(int lc = 0; lc < max(2,colors); lc++)
		  printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
			 lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
		fflush(stdout);
	      }
	      assert(fabs(len1-len2) < 0.001);
	    }
	  }

	  if(MapImageLoc){
	    if(!pmap->blockmem){
	      delete [] pmap->ImageFOV[c];
	      delete [] pmap->ImageX[c];
	      delete [] pmap->ImageY[c];
	      if(DEBUG>=2){
		pmap->ImageFOV[c] = NULL;
		pmap->ImageX[c] = NULL;
		pmap->ImageY[c] = NULL;
	      }
	    }
	    if(DEBUG>=4 && (colors >= 2 || usecolor)){
	      for(int i = orignummaps; i < nummaps; i++){
		Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
		double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
		double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
		if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
		  printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: M1 : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
		  for(int lc = 0; lc < max(2,colors); lc++)
		    printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d\n",
			   lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
		  fflush(stdout);
		}
		assert(fabs(len1-len2) < 0.001);
	      }
	    }

	    pmap->ImageFOV[c] = (int *)pmem;
	    pmem += sizeof(int) * (((numsite + 2) + 1) & ~1);
	    if(DEBUG>=2) pmap->ImageFOV[c][0] = pmap->ImageFOV[c][numsite+1] = -1;
	    for(int J = 1; J <= numsite; J++){
	      pmap->ImageFOV[c][J] = pmap->origImageFOV[c][J];
	      if(DEBUG>=2) assert(pmap->ImageFOV[c] && 0 <= pmap->ImageFOV[c][J] &&  pmap->ImageFOV[c][J] <= MaxImageFOV);
	    }

	    if(DEBUG>=4 && (colors >= 2 || usecolor)){
	      for(int i = orignummaps; i < nummaps; i++){
		Cmap *pmap = Map[i];// WAS12 Gmap[i]
      
		double len1 = pmap->origsite[0][pmap->orignumsite[0]+1];
		double len2 = pmap->origsite[1][pmap->orignumsite[1]+1];
		if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
		  printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: M2 : i=%d:colors=%d:pmap->id=%lld:\n",T,pmem,c,i,colors,pmap->id);
		  for(int lc = 0; lc < max(2,colors); lc++)
		    printf("  lc=%d:pmap->orignumsite[lc]=%d,pmap->origsite[lc][pmap->orignumsite[lc]+1] = %0.3f, pmap->site[lc]=%p,pmap->origsite[lc]=%p,pmap->blockmem=%d,\n",
			   lc,pmap->orignumsite[lc],pmap->origsite[lc][pmap->orignumsite[lc]+1],pmap->site[lc],pmap->origsite[lc],pmap->blockmem);
		  fflush(stdout);
		}
		assert(fabs(len1-len2) < 0.001);
	      }
	    }

	    pmap->ImageX[c] = (double *)pmem;
	    pmem += sizeof(double) * (numsite + 2);
	    for(int J = 1; J <= numsite; J++)
	      pmap->ImageX[c][J] = pmap->origImageX[c][J];

	    if(DEBUG>=4 && (colors >= 2 || usecolor)){
	      for(int i = orignummaps; i < nummaps; i++){
		Cmap *p = Map[i];// WAS12 Gmap[i]
      
		double len1 = p->origsite[0][p->orignumsite[0]+1];
		double len2 = p->origsite[1][p->orignumsite[1]+1];
		if((BIAS_TRACE && p->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
		  printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: M3 : i=%d:colors=%d:p->id=%lld:\n",T,pmem,c,i,colors,p->id);
		  for(int lc = 0; lc < max(2,colors); lc++)
		    printf("  lc=%d:p->orignumsite[lc]=%d,p->origsite[lc][p->orignumsite[lc]+1] = %0.3f, p->site[lc]=%p,p->origsite[lc]=%p..%p,p->blockmem=%d,pmap->origImageY[c]=%p..%p,numsite=%d\n",
			   lc,p->orignumsite[lc],p->origsite[lc][p->orignumsite[lc]+1],p->site[lc],p->origsite[lc],&p->origsite[lc][p->orignumsite[lc]+1],p->blockmem,pmap->origImageY[c],&pmap->origImageY[c][numsite+1],numsite);
		  fflush(stdout);
		}
		assert(fabs(len1-len2) < 0.001);
	      }
	    }

	    pmap->ImageY[c] = (double *)pmem;
	    pmem += sizeof(double) * (numsite + 2);
	    for(int J = 1; J <= numsite; J++)
	      pmap->ImageY[c][J] = pmap->origImageY[c][J];
	  }
	}// if(BNX)

	if(DEBUG>=4 && (colors >= 2 || usecolor)){
	  for(int i = orignummaps; i < nummaps; i++){
	    Cmap *p = Map[i];// WAS12 Gmap[i]
      
	    double len1 = p->origsite[0][p->orignumsite[0]+1];
	    double len2 = p->origsite[1][p->orignumsite[1]+1];
	    if((BIAS_TRACE && p->id == BIAS_TRACE_ID) || (DEBUG && !(fabs(len1-len2) < 0.001))){
	      printf("mapcorrectApply during T=%d site restore (pmem=%p),c=%d: M4 : i=%d:colors=%d:p->id=%lld:\n",T,pmem,c,i,colors,p->id);
	      for(int lc = 0; lc < max(2,colors); lc++)
		printf("  lc=%d:p->orignumsite[lc]=%d,p->origsite[lc][p->orignumsite[lc]+1] = %0.3f, p->site[lc]=%p,p->origsite[lc]=%p,p->blockmem=%d\n",
		       lc,p->orignumsite[lc],p->origsite[lc][p->orignumsite[lc]+1],p->site[lc],p->origsite[lc],p->blockmem);
	      fflush(stdout);
	    }
	    assert(fabs(len1-len2) < 0.001);
	  }
	}
      }

      if(DEBUG) assert(pmem <= &mem[memsiz]);
      pmap->blockmem = -1;/* marks both site[],SNR[] etc as well as siteSD[],sitecov[],sitecnt[]...lnSNRsd[],ChimQuality[] as block allocated */
    }

    if(VERB>=2){
      printf("mapcorrectApply: restored sites (nummaps=%d):time=%0.6f(wall=%0.6f)\n", nummaps, mtime(), wtime());
      if(minSNRestimate)
	printf("   minSNR[0]=%0.3f,minSNR[1]=%0.3f\n",minSNR[0],minSNR[1]);
      fflush(stdout);
    }
    delete [] mapoffset;
  } // if(restore)

  if(DEBUG>=3 && MapSNR){/* print map id BIAS_TRACE_ID */
    for(int m = orignummaps; m < nummaps; m++){
      Cmap *pmap = Map[m];
      if((BIAS_TRACE && pmap->id == BIAS_TRACE_ID) ||
	 !((pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3) ||
	 !((pmap->origsite[0][pmap->orignumsite[0]+1] - pmap->origsite[1][pmap->orignumsite[1]+1]) < 1e-3)){
	for(int c = 0; c < colors; c++){
	  int M = pmap->numsite[c];
	  double *X = pmap->site[c];
	  printf("mapcorrectApply after restore:mapid=%d,id=%lld,c=%d,M=%d:\n",m,pmap->id,c,M);
	  for(int J = 0; J <= M+1; J++)
	    printf(" J=%d:site[c][J]=%0.4f,SNR[c][J]=%0.3f\n",J,X[J],(1 <= J && J <= M) ? pmap->SNR[c][J] : -1.0);

	  M = pmap->orignumsite[c];
	  X = pmap->origsite[c];
	  printf("mapcorrectApply after restore:mapid=%d,id=%lld,c=%d,origM=%d:\n",m,pmap->id,c,M);
	  for(int J = 0; J <= M+1; J++)
	    printf(" J=%d:origsite[c][J]=%0.4f,origSNR[c][J]=%0.3f\n",J,X[J],(1 <= J && J <= M) ? pmap->origSNR[c][J] : -1.0);
	}
	fflush(stdout);
	assert((pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3);
	assert((pmap->origsite[0][pmap->orignumsite[0]+1] - pmap->origsite[1][pmap->orignumsite[1]+1]) < 1e-3);
      }
    }
  }

  double resKB = mres * 0.500;
  double IresSD = (mresSD > 0) ? 1.0/(mresSD* 0.500 * sqrt(2.0)) : 1.0;

  if(DEBUG) assert(origPixelLen > 0.0);
  if(DEBUG) assert(PixelLen > 0.0);
  double C = PixelLen/origPixelLen;
  int start = orignummaps;
  if(PairSplit && PairsplitRef){/* Do NOT apply bpp to first set of -i maps, since these are treated as the reference */
    if(Cfirst > 0)
      start = Cfirst;
    else
      start = nummaps;
  }
  if(VERB>=2){
    printf("mapcorrectApply:PixelLen=%0.9f,origPixelLen=%0.9f,C=%0.12f,start=%d\n",PixelLen,origPixelLen, C,start);
    if(minSNR[0] > 0.0 || (colors >= 2 && minSNR[1] > 0.0)){
      if(colors >= 2)
        printf("mapcorrectApply: minSNR[0]= %0.6f, minSNR[1]= %0.6f\n",minSNR[0],minSNR[1]);
      else
        printf("mapcorrectApply: minSNR[0]= %0.6f\n",minSNR[0]);
    }
    fflush(stdout);
  }

  int mergecnt=0;
    
  size_t orignumsites=0, numsites=0, trimends = 0, breakcnt = 0;
  FLOAT trimlen = 0.0;

  if(MaxEnd >= 0.0 && MaxEnd < MININTERVAL){
    printf("Increasing -maxEnd to %0.5f\n", MININTERVAL);
    fflush(stdout);
    MaxEnd = MININTERVAL;
  }

  int startnummaps = nummaps;

  long long maxID = 0;
  for(int T = orignummaps; T < startnummaps; T++)
    maxID = max(maxID, Map[T]->id);

  #pragma omp parallel num_threads(numthreads)
  {
    size_t my_orignumsites=0, my_numsites=0, my_trimends = 0, my_breakcnt = 0;
    FLOAT my_trimlen = 0.0;

    #pragma omp for schedule(dynamic,512)
    for(int T = orignummaps; T < startnummaps; T++){
      Cmap *pmap = Map[T];
      Ccontig *pcontig = pmap->contig;
      if(DEBUG) assert(pmap->id > 0);

      if(DEBUG && NoBpp) assert(fabs(PixelLen - origPixelLen) <= 1e-12);

      // Apply bpp scaling factor PixelLen/origPixelLen */
      if(PixelLen > 0.0 && fabs(PixelLen - origPixelLen) > 1e-12 && !NoBpp && T >= start){
	for(int c=0; c < colors; c++){
	  FLOAT *X = pmap->site[c];
	  int N = pmap->numsite[c];
	  if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID && Map == Gmap)
	    printf("mapcorrectApply:mapid=%d,id=%lld,c=%d,M=%d:After adjusting for PixelLen=%0.9f:\n",T,pmap->id,c,N,PixelLen);
	  for(int j= 0; j <= N+1; j++){
	    if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID && Map == Gmap){
	      if(minSNRestimate)
		printf("site[%d][%d]= %0.10f -> %0.10f (C= %0.12f) SNR=%0.3f\n", c, j, X[j], X[j] * C, C, (1 <= j && j <= N) ? pmap->SNR[c][j] : -1.0);
	      else
		printf("site[%d][%d]= %0.10f -> %0.10f (C= %0.12f)\n", c, j, X[j], X[j] * C, C);
	    }
	    X[j] *= C;
	  }
	}	  
	pmap->len *= C;// NEW
      }

      // Filter out sites below minSNR[c]
      for(int c= 0; c < colors; c++){
	FLOAT *Y =  pmap->site[c];
	int N = pmap->numsite[c];
	my_orignumsites += N;

	if(BNX && MapSNR && minSNR[c] > 0.0){
	  if(DEBUG>=2)  assert(pmap->sitecnt[c]== NULL);/* since map is BNX */
	  if(DEBUG>=2)  assert(pmap->SNRcnt[c]== NULL);/* since map is BNX */
	  int J = 0;
	  double *SNR = pmap->SNR[c];
	  for(int I = 1; I <= N; I++){
	    if(SNR[I] < minSNR[c]){
	      if(MapStitched)
		pmap->stitch[c][I+1] |= pmap->stitch[c][I];
	      continue;
	    }
	    if(++J < I){
	      Y[J] = Y[I];
	      SNR[J] = SNR[I];
	      if(MapTrueSite || Tracksites){
                pmap->truesiteL[c][J] = pmap->truesiteL[c][I];
		pmap->truesiteR[c][J] = pmap->truesiteR[c][I];
              }
	      if(MapIntensity)
		pmap->Intensity[c][J] = pmap->Intensity[c][I];
	      if(MapStitched)
		pmap->stitch[c][J] = pmap->stitch[c][I];
	      if(MapPSFWidth)
		pmap->PSFWidth[c][J] = pmap->PSFWidth[c][I];
	      if(MapImageLoc){
		pmap->ImageFOV[c][J] = pmap->ImageFOV[c][I];
		if(DEBUG>=2) assert(pmap->ImageFOV[c] && 0 <= pmap->ImageFOV[c][J] && pmap->ImageFOV[c][J] <= MaxImageFOV);
		pmap->ImageX[c][J] = pmap->ImageX[c][I];
		pmap->ImageY[c][J] = pmap->ImageY[c][I];
	      }
	    }
	  }
	  Y[J+1] = Y[N+1];/* length of molecule remains unchanged */
	  if(BIAS_TRACE && J < N && pmap->id== BIAS_TRACE_ID && Map == Gmap){
	    printf("mapcorrectApply: map %d (id=%lld), c=%d: %d of %d sites deleted due to -minSNR[%d] %0.4f:\n",T,pmap->id,c,N-J,N,c,minSNR[c]);
	    for(int k = 0; k <= J+1; k++)
	      printf("  k=%d/%d:site=%0.6f,SNR=%0.3f\n",k,J,Y[k],(1 <= k && k <= J) ? pmap->SNR[c][k] : -1.0);
	    fflush(stdout);
	  }
	  N = pmap->numsite[c] = J;

	  if(DEBUG>=2){/* check that map is monotonic */
	    assert(Y[0] == 0.0);
	    for(int k = 0; k <= N; k++)
	      if(Y[k+1] < Y[k]){
                #pragma omp critical
		{
		  printf("mapcorrectApply:map %d (id=%lld), color %d:site%d = %0.6f, site%d = %0.6f (out of order) in %s (len=%0.6f,N=%d)\n",
			 T, pmap->id, c+1, k, Y[k], k+1, Y[k+1], vfixx_filename[pmap->fileid], pmap->len,N);
		  fflush(stdout);
		  assert(Y[k+1] >= Y[k]);
		  fflush(stdout);exit(1);
		}
	      }
	  }
	}
      }

      // Apply mres,mresSD
      for(int c= 0; c < colors; c++){
	FLOAT *Y =  pmap->site[c];
	int N = pmap->numsite[c];

	for(int I=1;I < N;I++){
	  int J,K;

	  //#pragma unroll(0)
          #pragma nounroll
	  for(J=I+1; J <= N; J++){
	    double r,p;
	    if(mresSD){
	      r = urandom();
	      p = Pr(Y[J]-Y[J-1],resKB,IresSD);
	      if(VERB>=2 && pmap->id == 7585)
		printf("Y[J]-Y[J-1]=%0.6f,mres=%0.6f kb,mresSD=%0.6f kb, IresSD=%0.6f:r=%0.6f,Pr=%0.6f\n",
		       Y[J]-Y[J-1],resKB,mresSD * 0.500, IresSD, r, p);
	      if(r < p)
		break;
	    } else if(Y[J]-Y[J-1] > resKB)
	      break;
	  }
	  if(--J > I){/* merge sites Y[I] .. Y[J] */
	    if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID && Map == Gmap){
              #pragma omp critical
	      {
		mergecnt += J-I;
		printf("mapcorrectApply:map id=%lld:merging sites %d .. %d (at %0.10f ... %0.10f) into site %0.10f sites merged=%d(sum=%d), file=%s\n",
		       pmap->id,I,J,Y[I],Y[J],Yc(Y,J,J-I),J-I,mergecnt, vfixx_filename[pmap->fileid]);
		fflush(stdout);
	      }
	    }

	    Y[I] = Yc(Y,J,J-I);
	    if(pcontig){/* must also deres the contig->site,contig->HapSite, and contig->HapDelta */
	      double *conY = pcontig->site[c];
	      conY[I] = Y[I];
	      if(pcontig->HapSite[c]){
		int *HapSite = pcontig->HapSite[c];
		for(int k = I+1; k <= J; k++)
		  HapSite[I] |= HapSite[k];
	      }
	      if(pcontig->HapDelta[c]){// NOTE : this is NOT exactly correct if HapSite[I] was changed from 1 or 2 to 3, but should be good enough if map is subsequently refined
		double *HapDelta = pcontig->HapDelta[c];
		for(int k = I+1; k < J; k++)
		  HapDelta[I] += HapDelta[k];
		HapDelta[I-1] += HapDelta[I] *= 0.5;
		HapDelta[I] += HapDelta[J];
	      }
	    }
	    double inv = 1.0/(J-I+1);
	    if(pmap->truesiteL[c]){/* take min of truesiteL[c][I..J] values (and max for truesiteR[c][I..J])*/
	      long long siteL = pmap->truesiteL[c][I];
	      long long siteR = pmap->truesiteR[c][I];
	      for(int k = I+1; k <= J; k++){
		siteL = min(siteL, pmap->truesiteL[c][k]);
		siteR = max(siteR, pmap->truesiteR[c][k]);
	      }
	      pmap->truesiteL[c][I] = siteL;
	      pmap->truesiteR[c][I] = siteR;
	    }
	    if(BNX && MapSNR){ /* average SNR[c][I..J] values */
	      for(int k = I+1; k <= J; k++)
		pmap->SNR[c][I] += pmap->SNR[c][k];
	      pmap->SNR[c][I] *= inv;
	    }
	    if(BNX && MapIntensity){/* Average Intensity[c][] values */
	      for(int k = I+1; k <= J; k++)
		pmap->Intensity[c][I] += pmap->Intensity[c][k];
	      pmap->Intensity[c][I] *= inv;
	    }
	    if(BNX && MapStitched)/* OR stitch[c][] values */
	      for(int k = I+1; k <= J; k++)
		pmap->stitch[c][I] |= pmap->stitch[c][k];
	    int IJ = (I+J)/2;
	    if(BNX && MapPSFWidth && IJ > I)/* copy values for midpoint IJ */
	      pmap->PSFWidth[c][I] = pmap->PSFWidth[c][IJ];
	    if(BNX && MapImageLoc && IJ > I){/* copy values for midpoint IJ */
	      pmap->ImageFOV[c][I] = pmap->ImageFOV[c][IJ];
	      if(DEBUG>=2) assert(pmap->ImageFOV[c] && 0 <= pmap->ImageFOV[c][I] && pmap->ImageFOV[c][I] <= MaxImageFOV);
	      pmap->ImageX[c][I] = pmap->ImageX[c][IJ];
	      pmap->ImageY[c][I] = pmap->ImageY[c][IJ];
	    }
	    if(pmap->siteSD[c]){ /* average siteSD, sitecov, sitecnt */
	      for(int k = I+1; k <= J; k++){
		pmap->siteSD[c][I] += pmap->siteSD[c][k];
		pmap->sitecov[c][I] += pmap->sitecov[c][k];
		pmap->sitecnt[c][I] += pmap->sitecnt[c][k];
	      }
	      pmap->siteSD[c][I] *= inv;
	      pmap->sitecov[c][I] *= inv;
	      pmap->sitecnt[c][I] *= inv;
	    }

	    if(pmap->ChimQuality[c]){/* average ChimQuality[c] : skip -ve values */
	      float sum = 0.0;
	      int  cnt = 0;
	      for(int k = I; k <= J; k++){
		if(pmap->ChimQuality[c][k] >= 0.0){
		  sum += pmap->ChimQuality[c][k];
		  cnt++;
		}
	      }
	      if(cnt > 0){
		pmap->ChimQuality[c][I] = sum / cnt;
		if(DEBUG/* HERE >=2 */) assert(pmap->ChimQuality[c][I] >= 0.0);
	      } else
		pmap->ChimQuality[c][I] = END_CHIMQUAL;
	    }
	    if(pmap->ChimNorm[c]){/* average ChimNorm[c] : skip -ve values */
	      float sum = 0.0;
	      int  cnt = 0;
	      for(int k = I; k <= J; k++){
		if(pmap->ChimNorm[c][k] >= 0.0){
		  sum += pmap->ChimNorm[c][k];
		  cnt++;
		}
	      }
	      if(cnt > 0)
		pmap->ChimNorm[c][I] = sum / cnt;
	      else
		pmap->ChimNorm[c][I] = END_CHIMQUAL;
	    }
	    if(pmap->SegDupL[c]){/* average SegDupL[c] : skip -ve values */
	      float sum = 0.0;
	      int  cnt = 0;
	      for(int k = I; k <= J; k++){
		if(pmap->SegDupL[c][k] >= 0.0){
		  sum += pmap->SegDupL[c][k];
		  cnt++;
		}
	      }
	      if(cnt > 0)
		pmap->SegDupL[c][I] = sum / cnt;
	      else
		pmap->SegDupL[c][I] = END_CHIMQUAL;
	    }
	    if(pmap->SegDupR[c]){/* average SegDupR[c] : skip -ve values */
	      float sum = 0.0;
	      int  cnt = 0;
	      for(int k = I; k <= J; k++){
		if(pmap->SegDupR[c][k] >= 0.0){
		  sum += pmap->SegDupR[c][k];
		  cnt++;
		}
	      }
	      if(cnt > 0)
		pmap->SegDupR[c][I] = sum / cnt;
	      else
		pmap->SegDupR[c][I] = END_CHIMQUAL;
	    }
	    if(pmap->FragileEndL[c]){/* average FragileEndL[c] : skip -ve values */
	      float sum = 0.0;
	      int  cnt = 0;
	      for(int k = I; k <= J; k++){
		if(pmap->FragileEndL[c][k] >= 0.0){
		  sum += pmap->FragileEndL[c][k];
		  cnt++;
		}
	      }
	      if(cnt > 0)
		pmap->FragileEndL[c][I] = sum / cnt;
	      else
		pmap->FragileEndL[c][I] = END_CHIMQUAL;
	    }
	    if(pmap->FragileEndR[c]){/* average FragileEndR[c] : skip -ve values */
	      float sum = 0.0;
	      int  cnt = 0;
	      for(int k = I; k <= J; k++){
		if(pmap->FragileEndR[c][k] >= 0.0){
		  sum += pmap->FragileEndR[c][k];
		  cnt++;
		}
	      }
	      if(cnt > 0)
		pmap->FragileEndR[c][I] = sum / cnt;
	      else
		pmap->FragileEndR[c][I] = END_CHIMQUAL;
	    }

	    if(pmap->OutlierFrac[c]){/* use max OutlierFrac[c] */
	      float val = END_CHIMQUAL;
	      for(int k = I; k <= J; k++)
		val = max(val,pmap->OutlierFrac[c][k]);
	      pmap->OutlierFrac[c][I] = val;
	    }
	    if(pmap->FragSd[c]){/* use max FragSd[c] */
	      float val = END_CHIMQUAL;
	      for(int k = I; k <= J; k++)
		val = max(val,pmap->FragSd[c][k]);
	      pmap->FragSd[c][I] = val;
	    }
	    if(pmap->ExpSd[c]){/* use max ExpSd[c] */
	      float val = END_CHIMQUAL;
	      for(int k = I; k <= J; k++)
		val = max(val,pmap->ExpSd[c][k]);
	      pmap->ExpSd[c][I] = val;
	    }
	    if(pmap->FragCov[c]){/* use max FragCov[c] */
	      float val = END_CHIMQUAL;
	      for(int k = I; k <= J; k++)
		val = max(val,pmap->FragCov[c][k]);
	      pmap->FragCov[c][I] = val;
	    }
	    if(pmap->FragChiSq[c]){/* use min FragChiSq[c] (highest significance) */
	      float val = 1.0;
	      for(int k = I; k <= J; k++)
		val = min(val,pmap->FragChiSq[c][k]);
	      pmap->FragChiSq[c][I] = val;
	    }

	    if(pmap->Mask[c])/* OR Mask[c] */
	      for(int k = I+1; k <= J; k++)
		pmap->Mask[c][I] |= pmap->Mask[c][k];

	    if(pmap->SNRcnt[c]){
	      for(int k = I+1; k <= J; k++){
		pmap->SNRgmean[c][I] += pmap->SNRgmean[c][k];
		pmap->lnSNRsd[c][I] += pmap->lnSNRsd[c][k];
	      }
	      pmap->SNRgmean[c][I] *= inv;
	      pmap->lnSNRsd[c][I] *= inv;

	      /* combine SNRcnt[] and SNRdist[] */
	      int cnt = 0;
	      for(int k = I; k <= J; k++){
		if(DEBUG && pmap->SNRcnt[c][k]) assert(pmap->SNRdist[c][k] != NULL);
		cnt += pmap->SNRcnt[c][k];
	      }
	      double *newSNRdist = new double[cnt];/* NOTE: this memory will be leaked, but this doesn't waste much memory, since SNRcnt[] is only used for consensus maps */
	      int ncnt = 0;
	      for(int k = I; k <= J; k++)
		for(int t = 0; t < pmap->SNRcnt[c][k]; t++)
		  newSNRdist[ncnt++] = pmap->SNRdist[c][k][t];
	      qsort(newSNRdist,cnt,sizeof(double),(intcmp*)doubleInc);	      /* sort the SNR values in ascending order */

	      if(VERB>=3){
		printf("Calling delete[] Map[%d]->SNRdist[%d][I=%d]=%p -> %p,SNRcnt[c][I]=%d->%d,J=%d\n",T,c,I,pmap->SNRdist[c][I],newSNRdist,pmap->SNRcnt[c][I],cnt,J);
		fflush(stdout);
	      }
	      if(0 /* HERE blockmem >= 0 */)
		for(int k = I; k <= J; k++)
		  delete pmap->SNRdist[c][k];
	      pmap->SNRcnt[c][I] = cnt;
	      pmap->SNRdist[c][I] = newSNRdist;
	      for(int k = I+1; k <= J; k++){
		pmap->SNRcnt[c][k] = 0;
		pmap->SNRdist[c][k] = NULL;
	      }
	    }

	    /* shift remaining sites */
	    for(K=I;J < N;){
	      Y[++K] = Y[++J];
	      if(pcontig){
		double *conY = pcontig->site[c];
		conY[K] = conY[J];
		if(pcontig->HapSite[c])
		  pcontig->HapSite[c][K] = pcontig->HapSite[c][J];
		if(pcontig->HapDelta[c])
		  pcontig->HapDelta[c][K] = pcontig->HapDelta[c][J];
	      }
	      if(pmap->truesiteL[c]){
		pmap->truesiteL[c][K] = pmap->truesiteL[c][J];
		pmap->truesiteR[c][K] = pmap->truesiteR[c][J];
	      }
	      if(BNX && MapSNR)
		pmap->SNR[c][K] = pmap->SNR[c][J];
	      if(BNX && MapIntensity)
		pmap->Intensity[c][K] = pmap->Intensity[c][J];
	      if(BNX && MapStitched)
		pmap->stitch[c][K] = pmap->stitch[c][J];
	      if(BNX && MapPSFWidth)
		pmap->PSFWidth[c][K] = pmap->PSFWidth[c][J];
	      if(BNX && MapImageLoc){
		pmap->ImageFOV[c][K] = pmap->ImageFOV[c][J];
		if(DEBUG>=2) assert(pmap->ImageFOV[c] && 0 <= pmap->ImageFOV[c][K] && pmap->ImageFOV[c][K] <= MaxImageFOV);
		pmap->ImageX[c][K] = pmap->ImageX[c][J];
		pmap->ImageY[c][K] = pmap->ImageY[c][J];
	      }
	      
	      if(pmap->siteSD[c]){
		pmap->sitecnt[c][K] = pmap->sitecnt[c][J];
		pmap->sitecov[c][K] = pmap->sitecov[c][J];
		pmap->siteSD[c][K] = pmap->siteSD[c][J];
	      }

	      if(pmap->ChimQuality[c])
		pmap->ChimQuality[c][K] = pmap->ChimQuality[c][J];
	      if(pmap->ChimNorm[c])
		pmap->ChimNorm[c][K] = pmap->ChimNorm[c][J];
	      if(pmap->SegDupL[c])
		pmap->SegDupL[c][K] = pmap->SegDupL[c][J];
	      if(pmap->SegDupR[c])
		pmap->SegDupR[c][K] = pmap->SegDupR[c][J];
	      if(pmap->FragileEndL[c])
		pmap->FragileEndL[c][K] = pmap->FragileEndL[c][J];
	      if(pmap->FragileEndR[c])
		pmap->FragileEndR[c][K] = pmap->FragileEndR[c][J];
	      if(pmap->OutlierFrac[c])
		pmap->OutlierFrac[c][K] = pmap->OutlierFrac[c][J];

	      if(pmap->FragSd[c])
		pmap->FragSd[c][K] = pmap->FragSd[c][J];
	      if(pmap->ExpSd[c])
		pmap->ExpSd[c][K] = pmap->ExpSd[c][J];
	      if(pmap->FragCov[c])
		pmap->FragCov[c][K] = pmap->FragCov[c][J];
	      if(pmap->FragChiSq[c])
		pmap->FragChiSq[c][K] = pmap->FragChiSq[c][J];

	      if(pmap->Mask[c])
		pmap->Mask[c][K] = pmap->Mask[c][J];

	      if(pmap->SNRcnt[c]){
		pmap->SNRgmean[c][K] = pmap->SNRgmean[c][J];
		pmap->lnSNRsd[c][K] = pmap->lnSNRsd[c][J];

		pmap->SNRcnt[c][K] = pmap->SNRcnt[c][J];
		if(DEBUG) assert(pmap->SNRdist[c][K] == NULL);
		if(VERB>=3){
		  printf("Updating Map[%d]->SNRdist[%d][K=%d]=%p -> %p,J=%d\n",T,c,K,pmap->SNRdist[c][K],pmap->SNRdist[c][J],J);
		  fflush(stdout);
		}
		pmap->SNRdist[c][K] = pmap->SNRdist[c][J];

		pmap->SNRcnt[c][J] = 0;
		pmap->SNRdist[c][J] = NULL;
	      }
	    }

	    /* copy last site at N+1 to K+1 */
	    Y[K+1] = Y[N+1];
	    if(pcontig)
	      pcontig->site[c][K+1] = pcontig->site[c][N+1];
	    if(pmap->Mask[c])
	      pmap->Mask[c][K+1] = pmap->Mask[c][N+1]; // NEW51
	    if(pmap->SNRcnt[c]){
	      pmap->SNRgmean[c][K+1] = pmap->SNRgmean[c][N+1];
	      pmap->lnSNRsd[c][K+1] = pmap->lnSNRsd[c][N+1];
	      pmap->SNRcnt[c][K+1] = pmap->SNRcnt[c][N+1];
	      if(VERB>=3){
		printf("Updating Map[%d]->SNRdist[%d][K=%d]=%p -> %p,N+1=%d\n",T,c,K+1,pmap->SNRdist[c][K+1],pmap->SNRdist[c][N+1],N+1);
		fflush(stdout);
	      }
	      if(DEBUG) assert(pmap->SNRdist[c][K+1] == NULL);
	      pmap->SNRdist[c][K+1] = pmap->SNRdist[c][N+1];
	      pmap->SNRcnt[c][N+1] = 0;
	      pmap->SNRdist[c][N+1] = NULL;
	    }
	    N = K;
	  }
	}

	my_numsites += N;

	if(DEBUG) assert(Y[N+1] == Y[pmap->numsite[c] + 1]);
	Y[N+1] = Y[pmap->numsite[c]+1];

	pmap->numsite[c] = N;
	if(pcontig){
	  if(DEBUG) assert(pcontig->site[c][pcontig->numsite[c]+1] == pcontig->site[c][N+1]);
	  pcontig->numsite[c] = N;
	  if(pcontig->HapDelta[0]){/* check that fabs(HapDelta[I]) does not exceed Y[I+1]-Y[I] */
	    double *HapDelta = pcontig->HapDelta[0];
	    for(int I = 0; I <= N; I++){
	      double delta = fabs(HapDelta[I]);
	      if(delta > (Y[I+1]-Y[I]) * 0.9999)
		HapDelta[I] = copysign((Y[I+1]-Y[I]) * 0.9999, HapDelta[I]);
	    }
	  }
	}

	if(DEBUG>=3){/* check that map is monotonic */
	  assert(Y[0] == 0.0);
	  for(int k = 0; k <= N; k++){
	    if(Y[k+1] < Y[k]){
              #pragma omp critical
	      {
	        printf("mapcorrectApply:map %d (id=%lld), color %d:site%d = %0.6f, site%d = %0.6f (out of order) in %s (len=%0.6f,N=%d)\n",
		   T, pmap->id, c+1, k, Y[k], k+1, Y[k+1], vfixx_filename[pmap->fileid], pmap->len,N);
		fflush(stdout);
		assert(Y[k+1] >= Y[k]);
		fflush(stdout);exit(1);
	      }
	    }
	  }
	}
      }

      if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID && Map == Gmap){/* print map id BIAS_TRACE_ID */
	for(int c = 0; c < colors; c++){
	  int M = pmap->numsite[c];
	  double *X = pmap->site[c];
	  printf("mapcorrectApply:mapid=%d,id=%lld,c=%d,M=%d:After minSNR[c]=%0.3f,mres=%0.6f kb,mresSD=%0.6fkb\n",T, pmap->id, c, M, minSNR[c],resKB, mresSD*0.500);
	  for(int J = 0; J <= M+1; J++)
	    if(MapSNR)
	      printf(" J=%d:site[c][J]=%0.6f,SNR[c][J]=%0.3f\n",J,X[J],(1 <= J && J <= M) ? pmap->SNR[c][J] : -1.0);	      
	    else
	      printf(" J=%d:site[c][J]=%0.6f\n",J,X[J]);	      
	}
	fflush(stdout);
      }

      if(ScanScale[0] || MapScale){
	FLOAT C = pmap->cumscale;
	int scan = 1;
	if(ScanScale[0]){
	  scan = pmap->UniqueScanId;
	  C *= ScanScale[0][scan].cumscale;
	}
	for(int c = 0; c < colors; c++){
	  FLOAT *X = pmap->site[c];
	  int M = pmap->numsite[c];
	  for(int j= M+1; j > 0; j--)
	    X[j] *= C;
	}
	if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID && Map == Gmap){/* print map id BIAS_TRACE_ID */
	  for(int c = 0; c < colors; c++){
	    int M = pmap->numsite[c];
	    double *X = pmap->site[c];
	    printf("mapcorrectApply:mapid=%d,id=%lld,c=%d,M=%d:After applying ScanScale[0][%d].cumscale=%0.12f\n",T,pmap->id, c, M, scan, C);
	    for(int J = 0; J <= M+1; J++)
	      printf(" J=%d:site[c][J]=%0.4f,SNR[c][J]=%0.3f\n",J,X[J],(1 <= J && J <= M) ? pmap->SNR[c][J] : -1.0);
	  }
	  fflush(stdout);
        }
      }

      // Apply MaxInterval to break molecules with larger internal intervals
      if(VERB>=2){
	printf("mapcorrectApply:mapid=%d,id=%lld:BNX=%d,num_reffiles=%d,MaxInterval= %0.4f\n",T,pmap->id,BNX,num_reffiles,MaxInterval);
	fflush(stdout);
      }
      if(MaxInterval > 0.0 && colors > 1){
        #pragma omp critical
	{
	  if(MaxInterval > 0.0){
	    printf("WARNING: -maxInterval not yet implemented for colors=%d: ignoring -maxInterval\n",colors);
	    fflush(stdout);

	    MaxInterval = 0.0;
	  }
        }
      }

      if((BNX || num_reffiles <= 0 || (refMaxInterval && pmap->fileid >= num_files)) && MaxInterval > 0.0){
	if(colors > 1){
          #pragma omp critical
	  {
	    printf("ERROR: -maxInterval not yet implemented for colors=%d\n",colors);
	    fflush(stdout);
	    //	  sleep(30);
	    exit(1);
	  }
	}
	
	Cmap *p = pmap;
	FLOAT *X = p->site[0];
	int M = p->numsite[0], J = 1;
	for(; J < M; J++){
	  if(X[J+1]-X[J] > MaxInterval){
	    if(VERB>=2){
	      #pragma omp critical
	      {
		printf("T=%d,id=%lld:M=%d,J=%d:X[J]=%0.4f,X[J+1]=%0.4f : interval=%0.4f exceeds MaxInterval=%0.3f\n",
		       T,p->id,M,J,X[J],X[J+1],X[J+1]-X[J],MaxInterval);
		fflush(stdout);
	      }
	    }
	    break;
	  }
	}

	if(J < M){
	  if(!BNX){/* Cmap needs to be broken into multiple pieces */
	    int next;

            #pragma omp critical(nummaps)
	    {
	      maxmapalloc(nummaps+1,maxmaps,Map,0,1);// NOTE : uses "#pragma omp critical(MAX_BLOCKS)"
	      next = nummaps++;
	    }

	    // right part will be saved in Map[next], which will be rechecked later in case it needs to be broken again

	    if(VERB>=2){
	      #pragma omp critical
	      {
		printf("splitting Map[%d] = %p (site[0]= %p, id=%lld) into 2 pieces, with right piece replacing Map[%d] = %p (site[0]=%p), maxID=%lld\n",
		       T,Map[T],Map[T]->site[0],Map[T]->id, next,Map[next],Map[next]->site[0],maxID);
		fflush(stdout);
	      }
	    }

	    long long pid = p->id;
	    int fileid = p->fileid;

	    if(1){
	      /* compute right part */
	      Cmap *nrmap = new Cmap(p);
	      nrmap->trim(X[J+1] - min((MaxEnd >= 0.0 ? MaxEnd : 100.0), X[J+1]-X[J]), X[M+1]);
	      Cmap *newrmap = new Cmap(nrmap);// compacts the memory for the trimmed map
	      newrmap->id = pid + maxID;
	      newrmap->mapid = next;
	      newrmap->fileid = fileid;
	      delete nrmap;// frees memory of original untrimmed map
	  
	      /* swap *Map[next] with *newrmap  and free sub-array memory of original *Map[next] */
	      Cmap tmp;
	      tmp = *Map[next];
	      *Map[next] = *newrmap;
	      tmp.allfree();// NOTE : this frees sub-array memory of original *Map[next] 
	      *newrmap = tmp;
	      delete newrmap;// just frees the Cmap struct
	    }

	    if(1){
	      /* compute left part of p == Map[T] */
	      if(DEBUG) assert(p == Map[T]);
	      p->trim(0.0, X[J] + min((MaxEnd >= 0.0 ? MaxEnd : 100.0), X[J+1]-X[J]));
	      Cmap *newrmapL = new Cmap(p);// compacts the mmeory for the trimmed map
	      newrmapL->id = pid;
	      newrmapL->fileid = fileid;
	      newrmapL->mapid = T;
	  
	      /* swap *Map[T] with *newrmap and free sub-array memory of original *Map[T] */
	      Cmap tmpL = *Map[T];
	      *Map[T] = *newrmapL;
	      tmpL.allfree();// NOTE : this frees sub-array memory of original *Map[T]
	      *newrmapL = tmpL;
	      delete newrmapL;// just frees the Cmap struct
	    }

	    p = pmap = Map[T];

	  } else {/* BNX map : keep just the larger piece : delete either labels <= J or >= J+1 */
	    if(J >= M-J){
	      p->numsite[0] = J;
	      p->len = X[J+1];
	    } else {
	      FLOAT shift = X[J];
	      X[0] = 0.0;
	      int t = 1;
	      for(J = J+1; J <= M; J++, t++){
		X[t] = X[J] - shift;
		if(MapTrueSite || Tracksites){
		  p->truesiteL[0][t] = p->truesiteL[0][J];
		  p->truesiteR[0][t] = p->truesiteR[0][J];
		}
		if(MapSNR)
		  p->SNR[0][t] = p->SNR[0][J];
		if(MapIntensity)
		  p->Intensity[0][t] = p->Intensity[0][J];
		if(MapStitched)
		  p->stitch[0][t] = p->stitch[0][J];
		if(MapPSFWidth)
		  p->PSFWidth[0][t] = p->PSFWidth[0][J];
		if(MapImageLoc){
		  p->ImageFOV[0][t] = p->ImageFOV[0][J];
		  p->ImageY[0][t] = p->ImageY[0][J];
		  p->ImageX[0][t] = p->ImageX[0][J];
		}
	      }
	      p->len = X[t] = X[M+1] - shift;
	      p->numsite[0] = t-1;
	    }
	  }

	  my_breakcnt++;
	  if(BIAS_TRACE && pmap->id== BIAS_TRACE_ID && Map == Gmap){
            int M = p->numsite[0];
            printf("mapcorrectApply: map %d (id=%lld): After truncating due to MaxInterval=%0.3f:\n",T,pmap->id,MaxInterval);
	    for(int k = 0; k <= M+1; k++)
	      if(MapSNR)
		printf("  k=%d/%d:site=%0.6f,SNR=%0.3f\n",k,M,X[k],(1 <= k && k <= M) ? p->SNR[0][k] : -1.0);
	      else
		printf("  k=%d/%d:site=%0.6f\n",k,M,X[k]);
	    fflush(stdout);
	  }
	}
      }

      // Trim ends longer than MaxEnd
      if(MaxEnd >= 0.0){
	FLOAT end = pmap->site[0][1];
        for(int c = 1; c < colors; c++)
	  end = min(end,pmap->site[c][1]);
	if(end > MaxEnd){/* trim excess left ends */
	  pmap->len = 0.0;
	  FLOAT trim = end - MaxEnd;
	  my_trimlen += trim;
	  my_trimends++;
	  for(int c = 0; c < colors; c++){
	    FLOAT *Y =  pmap->site[c];
	    int N = pmap->numsite[c];
	    for(int I = 1; I <= N+1; I++)
	      Y[I] -= trim;
	    pmap->len = max(pmap->len,Y[N+1]);// NEW
	  }
	}
	
	int N = pmap->numsite[0];
	end = pmap->site[0][N+1] - pmap->site[0][N];
	for(int c = 1; c < colors; c++){
	  N = pmap->numsite[c];
	  end = min(end,pmap->site[c][N+1] - pmap->site[c][N]);
	}
	if(end > MaxEnd){/* trim excess right ends */
	  pmap->len = 0.0;
	  FLOAT trim = end - MaxEnd;
	  my_trimlen += trim;
	  my_trimends++;
	  for(int c = 0; c < colors; c++){
	    FLOAT *Y =  pmap->site[c];
	    int N = pmap->numsite[c];
	    Y[N+1] -= trim;
	    pmap->len = max(pmap->len,Y[N+1]);
	  }
	}

	if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID && Map == Gmap){
	  for(int c = 0; c < colors; c++){
	    FLOAT *X =  pmap->site[c];
	    int N = pmap->numsite[c];
	    printf("mapcorrectApply:mapid=%d,id=%lld,c=%d,M=%d:After adjusting for maxEnd=%0.6f:\n",T,pmap->id,c,N,MaxEnd);
	    for(int j= 0; j <= N+1; j++){
	      if(MapSNR)
		printf("site[%d][%d]= %0.6f, SNR=%0.3f\n", c, j, X[j], (1 <= j && j <= N) ? pmap->SNR[c][j] : -1.0);
	      else
		printf("site[%d][%d]= %0.6f\n", c, j, X[j]);
	    }
	  }
	}
      }
    } // for T = orignummaps .. startnummaps - 1 
    
    #pragma omp critical
    {
      orignumsites += my_orignumsites;
      numsites += my_numsites;
      trimlen += my_trimlen;
      trimends += my_trimends;
      breakcnt += my_breakcnt;
    }
  }

  if((BNX || num_reffiles <= 0 || (refMaxInterval && num_reffiles > 0)) && MaxInterval > 0.0 && nummaps > startnummaps){/* check newly created maps in case they still have intervals over MaxInterval */
    for(int T = startnummaps; T < nummaps; T++){
      Cmap *p = Map[T];
      if(!BNX && num_reffiles > 0 && refMaxInterval && p->fileid < num_files)
	continue;
      if(DEBUG) assert(p->id > 0);
      FLOAT *X = p->site[0];
      int M = p->numsite[0], J = 1;

      for(; J < M; J++){
	if(X[J+1]-X[J] > MaxInterval){
	  if(VERB>=2){
	    printf("T=%d,id=%lld:M=%d,J=%d:X[J]=%0.4f,X[J+1]=%0.4f : interval=%0.4f exceeds MaxInterval=%0.3f\n",
		   T,p->id,M,J,X[J],X[J+1],X[J+1]-X[J],MaxInterval);
	    fflush(stdout);
	  }
	  break;
	}
      }
      
      if(J < M){
        int next;

        #pragma omp critical(nummaps)
	{
	  maxmapalloc(nummaps+1,maxmaps,Map,0,1);// NOTE : uses "#pragma omp critical(MAX_BLOCKS)"
	  next = nummaps++;
	}

	long long pid = p->id;
	int fileid = p->fileid;

	if(1){
	  /* compute right part */
	  Cmap *nrmap = new Cmap(p);
	  nrmap->trim(X[J+1] - min((MaxEnd >= 0.0 ? MaxEnd : 100.0), X[J+1]-X[J]), X[M+1]);
	  Cmap *newrmap = new Cmap(nrmap);// compacts the memory for the trimmed map
	  newrmap->id = pid + maxID;
	  newrmap->fileid = fileid;
	  newrmap->mapid = next;
	  delete nrmap;// frees memroy of original untrimmed map
	  

	  /* swap *Map[next] with *newrmap and free sub-array memory of original *Map[next] */
	  Cmap tmp = *Map[next];
	  *Map[next] = *newrmap;
	  tmp.allfree();
	  *newrmap = tmp;
	  delete newrmap;// just frees the Cmap struct
	}

	if(1){
	  /* compute left part of p == Map[T] */
	  if(DEBUG) assert(p == Map[T]);
	  p->trim(0.0, X[J] + min((MaxEnd >= 0.0 ? MaxEnd : 100.0), X[J+1]-X[J]));
	  Cmap *newrmapL = new Cmap(p);// compacts the mmeory for the trimmed map
	  newrmapL->id = pid;
	  newrmapL->fileid = fileid;
	  newrmapL->mapid = T;

	  /* swap *Map[T] with *newrmap AND free sub-array memory of original *Map[t] */
	  Cmap tmpL = *Map[T];
	  *Map[T] = *newrmapL;
	  tmpL.allfree();
	  *newrmapL = tmpL;
	  delete newrmapL;// just frees the Cmap struct
	}

	//	p = Map[T];
	
	breakcnt++;
      }
    }
  }

  if(VERB>=2 && trimends > 0 && trimlen > 0.0){
    printf("mapcorrectApply: trimmed %lu molecule ends by an average of %0.4f kb\n", trimends, trimlen/(double)trimends);
    fflush(stdout);
  }
  if(VERB>=2 && orignumsites != numsites){
    if(colors==2)
      printf("reduced sites from %lu to %lu in %d maps due to minSNR=%0.4f,%0.4f,mres=%0.3f,mresSD=%0.3f:time=%0.6f(wall=%0.6f)\n",
	     orignumsites, numsites, nummaps-orignummaps, minSNR[0], minSNR[1],mres, mresSD, mtime(),wtime());
    else
      printf("reduced sites from %lu to %lu in %d maps due to minSNR=%0.4f,mres=%0.3f,mresSD=%0.3f:time=%0.6f(wall=%0.6f)\n",
	     orignumsites, numsites, nummaps-orignummaps, minSNR[0], mres, mresSD, mtime(),wtime());
    fflush(stdout);
  }
  if(VERB && (VERB>=2 || !BNX) && breakcnt > 0){
    printf("%lu molecules truncated due to -maxInterval %0.3f\n",breakcnt,MaxInterval);
    fflush(stdout);
  }

  if(VERB>=3){/* print map id BIAS_TRACE_ID */
    for(int m = orignummaps; m < nummaps; m++){
      Cmap *pmap = Map[m];
      if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
	int M = pmap->numsite[0];
	double *X = pmap->site[0];
	for(int J = 0; J <= M+1; J++)
	  printf(" J=%d:site[0][J]=%0.6f,SNR[0][J]=%0.3f,\n",J,X[J],(1 <= J && J <= M) ? pmap->SNR[0][J] : -1.0);
	fflush(stdout);
	break;
      }
    }
  }

  colors = origcolors;
}

int oldmaxmaps = 0;

// static int LastLine = -1;

static int idfiltered = -1;

static int RunDataIdCnt= 0, origRunDataIdCnt= 0;/* number of "Run Data" lines read from current file (origRunDataIdCnt are the number in the header rather than in the body of the file) */
static int RunDataIdMax = 0;
static int *RunDataIdmap = NULL;/* mapping from RunData line id in current file to index in RunDataList[] */


/* sort in order of ChipId, FlowCell, date, AM_PM, time (as txt except for numerical FlowCell */
static int ChipCellTimeInc(CRunDataSort *p1, CRunDataSort *p2)
{
  int cmp;
  if(p1->ChipId > p2->ChipId)
    return 1;
  if(p1->ChipId < p2->ChipId)
    return -1;
  if((cmp = p1->FlowCell - p2->FlowCell))
    return cmp;
  if((cmp = strcmp(p1->date,p2->date)))
    return cmp;
  if((cmp = strcmp(p1->AM_PM,p2->AM_PM)))
    return cmp;
  if((cmp = strcmp(p1->time,p2->time)))
    return cmp;
  return (p1 > p2) ? 1 : (p1 < p2) ? -1 : 0;// NEW : stable sort since V3 uses many default fields for ChipId, date, AM_PM, time
}

static FLOAT aNaN = nan("NaN");

static int MapWTpresent = -1;

static int tMaxSites;
static double tMinLen,tMaxLen;

static double fastfilterCnt = 0;// NEW102

static void parse_lines(int numthreads, Cline *lines,int numlines,int lines_per_map,int fileid, char *filename, int &nummaps, int &maxmaps, Cmap **&Map, std::set<long long> *ids, 
			int bnxFXY, int bnxChim)
{
  if(VERB>=2){
    printf("Starting parse_lines():fileid=%d,filename=%s,numthreads=%d,numlines=%d(line %d .. %d),lines_per_map=%d:time=%0.6f(wall=%0.6f)\n",
	   fileid,filename,numthreads,numlines, lines[0].linecnt, numlines ? lines[numlines-1].linecnt : lines[0].linecnt-1, lines_per_map, mtime(),wtime());
    fflush(stdout);
  }
  if(DEBUG) assert((numlines % lines_per_map) == 0);

  if(ids){/* do a quick parse of the first line of each group of lines to find the map id and decide if it should be discarded */
    int remlines = 0;/* number of lines remaining after filtering out lines that don't match the set in "ids" */
    for(int line = 0; line <= numlines - lines_per_map; line += lines_per_map){
      int linecnt = lines[line].linecnt;
      char *buf = lines[line].buf;
      
      if(VERB>=2){
	printf("remlines=%d,line=%d/%d: lines in this group=%d (starting at line %d):\n",remlines,line,numlines,lines_per_map,linecnt);
	for(int t = 0; t < lines_per_map; t++)
	  printf("%s\n",lines[line+t].buf);
	fflush(stdout);
      }

      if(buf[0] != '0'){
	printf("Header line does not start with '0' on line %d of %s (lines_per_map=%d)\n%s\n",linecnt, filename, lines_per_map, buf);
	fflush(stdout);exit(1);
      }
      char *pt = buf, *qt;
      int LabelChannel = strtol(pt,&qt,10);
      if(pt == qt || LabelChannel != 0){
	printf("Expected integer value 0 at start of line %d of %s:\n%s\n",linecnt,filename,buf);
	fflush(stdout);exit(1);
      }
      pt = qt;
      long long id = strtoll(pt,&qt,10);
      if(pt == qt || id <= 0){
	printf("Invalid MapID on line %d of %s:\n%s\n",linecnt,filename,buf);
	fflush(stdout);exit(1);
      }
      if(ids->find(id) == ids->end()){/* this map id is not needed */
	if(VERB>=3){
	  printf("    map id=%lld was filtered out since it is not in the set of contigs being refined\n",id);
	  fflush(stdout);
	}
	idfiltered++;
	continue;
      }
      if(VERB>=2){
	printf("map id=%lld is contigs id set\n", id);
	fflush(stdout);
      }

      if(line > remlines)/* copy lines for this map */
	for(int t = 0; t < lines_per_map; t++)
	  lines[remlines+t] = lines[line+t];
      remlines += lines_per_map;
    }
    numlines = remlines;
  }

  if(DEBUG) assert(origPixelLen > 0.0);
  if(DEBUG) assert(PixelLen > 0.0);
  double C = PixelLen/origPixelLen;// bpp correction to length when applying MinLen,MaxLen

  int oldmaps = nummaps;
  int newmaps = (numlines+lines_per_map-1)/lines_per_map;
  //  oldmaxmaps = maxmaps;
  if(nummaps > (int)MASK(31) - newmaps){
    printf("Too many maps in %s up to line %d: Number of maps limited to 2 billion (%d)\n",filename, lines[numlines-1].linecnt,MASK(31));
    fflush(stdout);exit(1);
  }
  
  int nthreads = min(numthreads, 16 /* WAS 4 */);

  if(VERB>=2){
    long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
    getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
    printf("Starting parse_lines:VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb : CPU time= %0.6f, wall time= %0.6f\n",VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(), wtime());
    printf("nummaps=%d -> %d,maxmaps=%d, minsite=%d, MapSNR=%d, MapIntensit=%d, MapStitched=%d, MapStitchLoc=%d,MapPSFWidth=%d,MapTrueSite=%d,Tracksites=%d\n",
    	   nummaps,nummaps+newmaps,maxmaps, (nthreads >= 4) ? MINSITE : 0 , MapSNR, MapIntensity, MapStitched, MapStitchLoc,MapPSFWidth,MapTrueSite,Tracksites);
    fflush(stdout);
  }

  maxmapalloc(nummaps+newmaps,maxmaps,Map, (nthreads > 4/* WAS101 > 1 */) ? MINSITE : 0, numthreads);
  nummaps += newmaps;

  if(DEBUG>=2)
    for(int i = oldmaps; i < nummaps; i++)
      assert(Map[i] != NULL);

  if(VERB>=2){
    long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
    getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
    printf("Entering parallel section2:VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb : CPU time= %0.6f, wall time= %0.6f\n",VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(), wtime());
    printf("\t numthreads=%d,nthreads=%d,numlines=%d\n",numthreads,nthreads,numlines);
    fflush(stdout);
  }
  
  #pragma omp parallel num_threads(nthreads) if(nthreads > 1)
  {
    double Location[MAXSITES];

    #pragma omp for schedule(dynamic,256/* WAS44 512 */)
    for(int line = 0; line <= numlines - lines_per_map; line += lines_per_map){

      int mapid = oldmaps + line/lines_per_map;

      /* read in next Molecule into Map[mapid] */
      if(VERB>=2){
	#pragma omp critical
	{
	  printf("Parsing %d lines starting at line %d, mapid=%d:\n",lines_per_map,lines[line].linecnt,mapid);
	  for(int t = 0; t < lines_per_map; t++)
	    printf("%s\n",lines[line+t].buf);
	  printf("\n");
	  fflush(stdout);
	}
      }

      int linecnt = lines[line].linecnt;
      char *buf = lines[line].buf;
      if(buf[0] != '0'){
        #pragma omp critical 
	{
	  printf("Header line does not start with '0' on line %d of %s (lines_per_map=%d)\n%s\n",linecnt, filename, lines_per_map, buf);
	  fflush(stdout);exit(1);
	}
      }
      
      assert(mapid < maxmaps);

      Cmap *pmap = Map[mapid];

      pmap->mapid = mapid; 
      pmap->fileid = fileid;
      pmap->linenum = linecnt;

      for(int c = 0; c < colors; c++)
	pmap->Nickase[c] = Nickase[c];

      pmap->startloc = pmap->endloc = 0.0;
      pmap->flipped = 0;

      char *pt = buf, *qt;
      int LabelChannel = strtol(pt,&qt,10);
      if(pt == qt || LabelChannel != 0){
        #pragma omp critical
	{
	  printf("Expected integer value 0 at start of line %d of %s:\n%s\n",linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
      }
    
      pt = qt;
      pmap->id = strtoll(pt,&qt,10);
      if(VERB>=2 /*|| pmap->id==8400128611401000078*/){
        #pragma omp critical
	{
	  printf("mapid=%d:id=%lld (line %d of %s)\n",
		 pmap->mapid,pmap->id,linecnt,filename);
	  fflush(stdout);
	}
      }
      if(pt == qt || pmap->id <= 0){
        #pragma omp critical
	{
	  printf("Invalid MapID on line %d of %s:\n%s\n",linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
      }

      pt = qt;
      double Len = strtod(pt,&qt); 
      if(pt == qt || Len <= 0.0){
        #pragma omp critical
	{
	  printf("Invalid Length on line %d of %s (must be > 0.0):\n%s\n",linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
      }
      pmap->len = Len * 0.001;

      int NumSite = -1;
      int NumSiteLine = -1;

      pmap->RunIndex = pmap->UniqueScanId = 0;// NEW103 : set default to handle BNX without scan information

      if(BNXVersion==0){
        if(MapWTpresent > 0){/* read mapWT on header line */
          pmap->mapWT = 1.0;
	  while(*qt && isspace(*qt))
	    qt++;
	  pt = qt;
	  if(!pt){
	    if(!MapWtWarning){
	      #pragma omp critical
	      {
		printf("WARNING: MapWT value missing on line %d of %s\n",linecnt,filename);
		fflush(stdout);
		MapWtWarning = 1;
	      }
	    }
	    if(bnxmergeerror){
              #pragma omp critical
	      {
		printf("Use -bnxmergeerror 0 if you want this error to be ignored (defaulting to MapWT = 1)\n");
		fflush(stdout);exit(1);
	      }
	    }
	  } else {
	    double MapWT = strtod(pt,&qt);
	    if(pt==qt || !(*qt==0 || isspace(*qt)) || MapWT < 0.0 || MapWT > 1.0){
              #pragma omp critical
	      {
		printf("Invalid MapWT on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->mapWT = MapWT;
	  }
        }
	if(HasLoc){/* check for StartLoc and EndLoc on header line */
	  while(*qt && isspace(*qt))
	    qt++;
	  if(*qt){
	    pt = qt;
	    pmap->startloc = strtod(pt,&qt);
	    if(pt == qt || pmap->startloc < 0.0){
              #pragma omp critical
	      {
		printf("Invalid Startloc on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->startloc *= 0.001;

	    pt = qt;
	    pmap->endloc = strtod(pt,&qt);
	    if(pt == qt || pmap->endloc < 0.0){
              #pragma omp critical
	      {
		printf("Invalid Endloc on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->endloc *= 0.001;

	    pmap->flipped = (pmap->id % 10) & 1;
	    
	    while(*qt && isspace(*qt))
	      qt++;

	    if(*qt){/* read in StartLoc2 and EndLoc2 Flipped ChimFlip */
	      pt = qt;
	      pmap->startloc2 = strtod(pt,&qt);
	      if(pt == qt || pmap->startloc2 < 0.0){
                #pragma omp critical
		{
		  printf("Invalid Startloc2 on line %d of %s:\n%s\n",linecnt,filename,buf);
		  fflush(stdout);exit(1);
	        }  
	      }
	      pmap->startloc2 *= 0.001;

	      pt = qt;
	      pmap->endloc2 = strtod(pt,&qt);
	      if(pt == qt || pmap->endloc2 < 0.0){
                #pragma omp critical
	        {
		  printf("Invalid Endloc2 on line %d of %s:\n%s\n",linecnt,filename,buf);
		  fflush(stdout);exit(1);
	        } 
	      }
	      pmap->endloc2 *= 0.001;
	      
	      pt = qt;
	      pmap->flipped = strtol(pt,&qt,10);
	      if(pt == qt || pmap->flipped < 0){
		#pragma omp critical
		{
		  printf("Invalid Flipped value on line %d of %s:\n%s\n",linecnt,filename,buf);
		  fflush(stdout);exit(1);
	        }
	      }
	      
	      pt = qt;
	      pmap->chimflip = strtol(pt,&qt,10);
	      if(pt == qt || pmap->chimflip < 0){
		#pragma omp critical
		{
		  printf("Invalid ChimFlip value on line %d of %s:\n%s\n",linecnt,filename,buf);
		  fflush(stdout);exit(1);
	        }
	      }
	    }
	  }
	}
      } else {/* read BNX version 1.0 header fields */
	pt = qt;
	pmap->AvgIntensity = strtod(pt,&qt);
	if(pt == qt || !(*qt==0 || isspace(*qt))){
	  #pragma omp critical
	  {
	    printf("Invalid AvgIntensity on line %d of %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	}
	
	pt = qt;
	pmap->bbSNR = strtod(pt,&qt);
	if(pt == qt || !(*qt==0 || isspace(*qt))){
	  #pragma omp critical
	  {
	    printf("Invalid Backbone SNR on line %d of %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	}

	pt = qt;
	NumSite = strtol(pt,&qt,10);
	if(pt==qt || !(*qt==0 || isspace(*qt))){
	  #pragma omp critical
	  {
	    printf("Invalid NumberofLabels on line %d of %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	}
	NumSiteLine = linecnt;

	pt = qt;
	pmap->OriginalMoleculeId = strtol(pt,&qt,10);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || pmap->OriginalMoleculeId < 0){
	  #pragma omp critical
	  {
	    printf("Invalid OriginalMoleculeId on line %d of %s (must be >= 0):\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	}

	pt = qt;
	pmap->ScanNumber = strtol(pt,&qt,10);
	if(pt==qt || !(*qt==0 || isspace(*qt))){
	  #pragma omp critical
	  {
	    printf("Invalid ScanNumber %d on line %d of %s (must be +ve number):\n%s\n",pmap->ScanNumber,linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	}
	if(pmap->ScanNumber <= 0){
	  if(ScanCorrection){
	    #pragma omp critical
	    {
	      printf("Invalid ScanNumber %d on line %d of %s (must be +ve number):\n%s\n",pmap->ScanNumber,linecnt,filename,buf);
	      fflush(stdout);exit(1);
	    }
	  }
	  if(!scannumber_warning){
	    #pragma omp critical
	    {
	      printf("WARNING: Invalid ScanNumber %d on line %d of %s (must be +ve number):\n%s\n",pmap->ScanNumber,linecnt,filename,buf);
	      printf("    Ignoring since -ScanScaling is not being used\n");
	      fflush(stdout);
	      scannumber_warning = 1;
	    }
	  }
	}

	pt = qt;
	pmap->ScanDirection = strtol(pt,&qt,10);
	if(pt==qt || !(*qt==0 || isspace(*qt))){
	  #pragma omp critical
	  {
	    printf("Invalid ScanDirection on line %d of %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	}
	
	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	while(*qt && !isspace(*qt))
	  qt++;
	if(pt == qt || !(*qt==0 || isspace(*qt))){
	  #pragma omp critical
	  {
	    printf("Invalid ChipId=%s on line %d of %s:\n%s\n",pt,linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	}
	char c = *qt;
	*qt = 0;
	// pmap->ChipId = strdup(pt);
	pmap->ChipId = ChipIdAlloc(pt);
	*qt = c;

	pt = qt;
	pmap->FlowCell = strtol(pt,&qt,10);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || pmap->FlowCell < 0){
	  #pragma omp critical
	  {
	    printf("Invalid FlowCell value \"%s\" on line %d of %s (must integer be >= 0):\n%s\n",pt,linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	}

	/* check if RunIndex and mapWT are present at end of line (optional) : If not use most recent Run Data in file (lines[line].RunIndex) (and mapWT=1.0) */
	pmap->RunIndex = lines[line].RunIndex;
	if(DEBUG && RunDataListLen > 0 && !(pmap->RunIndex >= 0 && pmap->RunIndex < RunDataListLen)){
	  #pragma omp critical
	  {
	    printf("line %d in file %s:RunIndex=%d,RunDataListLen=%d\n",lines[line].linecnt,filename,lines[line].RunIndex,RunDataListLen);
	    fflush(stdout);
	    assert(pmap->RunIndex >= 0 && pmap->RunIndex < RunDataListLen);
	  }
	}
	pmap->mapWT = 1.0;

	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	if(!*qt){
	  if(origRunDataIdCnt > 1 || bnxFXY > 0 || bnxChim > 0 || MapWTpresent > 0){
	    if(!RunIndexWarning){
              #pragma omp critical
              {
	        if(!RunIndexWarning){
		  printf("WARNING: Run Index missing on line %d of %s (required since multiple (%d) Run Data line are present in file header)\n", linecnt, filename, origRunDataIdCnt);
		  printf("    BNXheaderFXY=%d,bnxFXY=%d BNXheaderChim=%d, bnxChim=%d, MapWTpresent=%d\n",BNXheaderFXY,bnxFXY,BNXheaderChim,bnxChim,MapWTpresent);
		  printf("    Suggest regenerating %s with newest RefAligner to add those fields : otherwise ScanNumbers will be wrong and map ids may have to be renumbered\n", filename);
		  fflush(stdout);
		  RunIndexWarning = 1;
		}
	      }
	    }
	    if(bnxmergeerror){
              #pragma omp critical
	      {
		printf("Use -bnxmergeerror 0 if you want this error to be ignored\n");
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->RunIndex = 1;
	  }
	} else {
	  int RunIndex = strtol(pt,&qt,10);
	  if(pt==qt || !(*qt==0 || isspace(*qt)) || (RunDataIdCnt > 0 ? (RunIndex < 1 || RunIndex > RunDataIdCnt) : RunIndex > 1)){
	    #pragma omp critical
	    {
	      printf("Invalid Run Index %d on line %d of %s (%d previous Run Data lines in file):\n%s%s\n",RunIndex, linecnt,filename,RunDataIdCnt,pt,buf);
	      fflush(stdout);exit(1);
	    }
	  }
	  if(RunDataIdCnt > 0){
	    pmap->RunIndex = RunDataIdmap[RunIndex-1];
	    if(DEBUG) assert(pmap->RunIndex >= 0 && pmap->RunIndex < RunDataListLen);
	  } else
	    pmap->RunIndex = pmap->UniqueScanId = 0;

	  /* check if Column, FOV & XY fields are present */
	  if(bnxFXY > 0){
	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing Column value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    int Column = strtol(pt,&qt,10);
	    if(pt==qt || !(*qt==0 || isspace(*qt)) || Column < 0){
	      #pragma omp critical
	      {
		printf("Invalid Column value %d on line %d of %s:\n%s\n",Column,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->Column = Column;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing StartFOV value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    int StartFOV = strtol(pt,&qt,10);
	    if(pt==qt || !(*qt==0 || isspace(*qt)) || StartFOV < 0){
	      #pragma omp critical
	      {
		printf("Invalid StartFOV value %d on line %d of %s:\n%s\n",StartFOV,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->StartFOV = StartFOV;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing StartX value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    int StartX = strtol(pt,&qt,10);
	    if(pt==qt || !(*qt==0 || isspace(*qt)) || StartX < 0){
	      #pragma omp critical
	      {
		printf("Invalid StartX value %d on line %d of %s:\n%s\n",StartX,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->StartX = StartX;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing StartY value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    int StartY = strtol(pt,&qt,10);
	    if(pt==qt || !(*qt==0 || isspace(*qt)) || StartY < 0){
	      #pragma omp critical
	      {
		printf("Invalid StartY value %d on line %d of %s:\n%s\n",StartY,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->StartY = StartY;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing EndFOV value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    int EndFOV = strtol(pt,&qt,10);
	    if(pt==qt || !(*qt==0 || isspace(*qt)) || EndFOV < 0){
	      #pragma omp critical
	      {
		printf("Invalid EndFOV value %d on line %d of %s:\n%s\n",EndFOV,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->EndFOV = EndFOV;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing EndX value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    int EndX = strtol(pt,&qt,10);
	    if(pt==qt || !(*qt==0 || isspace(*qt)) || EndX < 0){
	      #pragma omp critical
	      {
		printf("Invalid EndX value %d on line %d of %s:\n%s\n",EndX,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->EndX = EndX;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing EndY value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    int EndY = strtol(pt,&qt,10);
	    if(pt==qt || !(*qt==0 || isspace(*qt)) || EndY < 0){
	      #pragma omp critical
	      {
		printf("Invalid EndY value %d on line %d of %s:\n%s\n",EndY,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->EndY = EndY;
	  } // if bnxFXY > 0

	  /* check if FeatureMinPos .. FeatureMaxScaled are present */
	  if(bnxChim > 0){
	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
              #pragma omp critical
	      {
	        printf("ERROR: Missing FeatureMinPos value on line %d of %s (present in header):\n%s\n",linecnt,filename,buf);
		fflush(stdout);
	      }
	    }
	    int FeatureMinPos = strtol(pt,&qt,10);
	    if(pt==qt || !(*qt==0 || isspace(*qt))){
	      #pragma omp critical
	      {
		printf("Invalid FeatureMinPos value %d on line %d of %s:\n%s\n",FeatureMinPos,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->FeatureMinPos = FeatureMinPos;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing FeatureMinScore value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    double FeatureMinScore = strtod(pt,&qt);
	    if(pt==qt || !(*qt==0 || isspace(*qt))){
	      #pragma omp critical
	      {
		printf("Invalid FeatureMinScore value %0.6f on line %d of %s:\n%s\n",FeatureMinScore,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->FeatureMinScore = FeatureMinScore;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing FeatureMaxPos value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    int FeatureMaxPos = strtol(pt,&qt,10);
	    if(pt==qt || !(*qt==0 || isspace(*qt))){
	      #pragma omp critical
	      {
		printf("Invalid FeatureMaxPos value %d on line %d of %s:\n%s\n",FeatureMaxPos,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->FeatureMaxPos = FeatureMaxPos;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing FeatureMaxScore value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    double FeatureMaxScore = strtod(pt,&qt);
	    if(pt==qt || !(*qt==0 || isspace(*qt))){
	      #pragma omp critical
	      {
		printf("Invalid FeatureMaxScore value %0.6f on line %d of %s:\n%s\n",FeatureMaxScore,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->FeatureMaxScore = FeatureMaxScore;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing FeatureMinScaled value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    double FeatureMinScaled = strtod(pt,&qt);
	    if(pt==qt || !(*qt==0 || isspace(*qt))){
	      #pragma omp critical
	      {
		printf("Invalid FeatureMinScaled value %0.6f on line %d of %s:\n%s\n",FeatureMinScaled,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->FeatureMinScaled = FeatureMinScaled;

	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      #pragma omp critical
	      {
		printf("ERROR: Missing FeatureMaxScaled value on line %d of %s:\n%s\n",linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    double FeatureMaxScaled = strtod(pt,&qt);
	    if(pt==qt || !(*qt==0 || isspace(*qt))){
	      #pragma omp critical
	      {
		printf("Invalid FeatureMaxScaled value %0.6f on line %d of %s:\n%s\n",FeatureMaxScaled,linecnt,filename,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->FeatureMaxScaled = FeatureMaxScaled;
	  }

	  /* check if mapWT is present */
	  if(MapWTpresent > 0){
	    while(*qt && isspace(*qt))
	      qt++;
	    pt = qt;
	    if(!*qt){            
	      if(!MapWtWarning){
		#pragma omp critical
		{
		  printf("WARNING: MapWT value missing on line %d of %s\n",linecnt,filename);
		  fflush(stdout);
		  MapWtWarning = 1;
		}
              }
	      if(bnxmergeerror){
                #pragma omp critical
	        {
                  printf("Use -bnxmergeerror 0 if you want this error to be ignored (defaulting to MapWT = 1)\n");
		  fflush(stdout);exit(1);
	        }
	      }
            }  else {
              double MapWT = strtod(pt,&qt);
	      if(pt==qt || !(*qt==0 || isspace(*qt)) || MapWT < 0.0 || MapWT > 1.0){
                #pragma omp critical
                {
                  printf("Invalid MapWT on line %d of %s:\n%s\n",linecnt,filename,buf);
		  fflush(stdout);exit(1);
                }
              }
	      pmap->mapWT = MapWT;
            }
          }

	}
	/* molecule id will be replaced by a unique number in input_vfixx_duplicates() unless RunIndexWarnin and RunDataWarning is set */
      } // read BNX version 1.0 header fields

      int totsites = 0, color;
      for(color = 1; color <= colors; color++){
	linecnt = lines[line+color].linecnt;
	buf = lines[line+color].buf;

	if(!isdigit(buf[0]) || buf[0] - '0' != color){
          #pragma omp critical
	  {
	    printf("Line for color %d does not start with '%c' on line %d of %s\n%s\n", color,(char)('0'+color),linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	}
    
	/* read in "LabelChannel LabelPosition ... Length" */

	/* scan over the input line to find the site locations (including right end) */
	int numsite = 0;
	for(pt = &buf[1]; *pt; pt = qt){
	  if(numsite >= MAXSITES){
            #pragma omp critical
	    {
	      printf("More than %d sites on line %d of file %s\n",numsite,linecnt,filename);
	      fflush(stdout);exit(1);
	    }
	  }
	  Location[numsite] = strtod(pt,&qt);
	  if(pt==qt){
	    while(*pt && isspace(*pt))
	      pt++;
	    if(!*pt)/* End of Line */
	      break;
            #pragma omp critical
	    {
	      printf("Invalid site value at char %ld of line %d in %s:%s\n%s\n",
		     pt-buf+1,linecnt,filename,pt,buf);
	      fflush(stdout);exit(1);
	    }
	  }
	  if(!(*qt == 0 || isspace(*qt))){
	    while(*pt && isspace(*pt))
	      pt++;
            #pragma omp critical
	    {
	      printf("Invalid site value at char %ld of line %d in %s:%s\n%s\n",
		     pt-buf+1,linecnt,filename,pt,buf);
	      fflush(stdout);exit(1);
	    }
	  }
	  if(Location[numsite] < 0.0 || Location[numsite] > Len){
            #pragma omp critical
  	    {
	      printf("%s:Out of bounds Label Location[%d]=%0.1f on line %d of file %s (map length in header was %0.1f):\n%s\n",
		     bnxerror ? "ERROR" : "WARNING", numsite,Location[numsite],linecnt,filename,Len,buf);
	      fflush(stdout);
	      if(bnxerror) exit(1);
	    }
	  }
	  Location[numsite] *= 0.001;

	  numsite++;
	}

	if(numsite < 1){
          #pragma omp critical
	  {
	    printf("WARNING: no locations found on line %d of %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);
	  }
	} else
	  numsite--;

	if(VERB>=3 && pmap->id==71){
          #pragma omp critical
	  {
	    printf("pmap->id=%lld:pmap->len= %0.4f, numsite= %d, Location[numsite]= %0.4f\n",pmap->id,pmap->len,numsite,Location[numsite]);
	    fflush(stdout);
	  }
	}

	for(int i = 1; i < numsite; i++){
	  if(Location[i] < Location[i-1] - 1e-6){// NEW43
            #pragma omp critical
  	    {
	      printf("ERROR:Out of order Label Locations %0.1f %0.1f on line %d of file %s (map length is %0.1f, numsite=%d):\n%s\n",
		     Location[i-1]*1000.0, Location[i]*1000.0,linecnt,filename,pmap->len*1000.0,numsite,buf);
	      fflush(stdout);
	      exit(1);
	    }
	  }
	}

	if(fabs(Location[numsite]*1000.0 - Len) > 0.1){
          #pragma omp critical
	  {
	    printf("WARNING: map length=%0.1f in header line does not match right end location=%0.1f on line %d of %s\n%s\n",
		   Len, Location[numsite]*1000.0, linecnt,filename,buf);
	    fflush(stdout);
	  }
	}

	for(int i = 0; i <= numsite; i++)// NEW43
	  if(Location[i] > pmap->len)
	    pmap->len = Location[i]; /* extend the right end of the map */
      
	if(VERB>=2 || (BIAS_TRACE && pmap->id==BIAS_TRACE_ID)){
	  #pragma omp critical
	  {
	    printf("pmap->id=%lld:color=%d:pmap->len= %0.6f, numsite=%d from line %d in %s\n%s\n",pmap->id,color,pmap->len, numsite,linecnt,filename,buf);
	    for(int i = 0; i <= numsite; i++)
	      printf("Location[%d]=%0.6f\n",i,Location[i]);
	    fflush(stdout);
	  }
	}

	if(BNX_FASTFILTER && !ids && (pmap->len * C < tMinLen || (tMaxLen > 0.0 && pmap->len * C > tMaxLen) ||
				      (minSNR[color-1] <= 0.0 && (colors==1 || usecolor == color) && 
				       (numsite < MinSites || (tMaxSites && numsite > tMaxSites) ||
					(0 /* WAS102 MaxSiteDensity > 0.0 && !minSNRestimate && numsite * 100.0 > MaxSiteDensity * pmap->len * C */))))){
	  // NOTE : cannot filter based on number of labels if minSNR > 0.0, since minSNR thresholding has not been applied 
	  if(VERB>=2){
	    #pragma omp critical
	    {
	      printf("pmap->id=%lld,mapid=%d: len= %0.3f -> %0.3f, numsite=%d : skipping due to -minlen -maxlen -minsites\n",pmap->id,mapid,pmap->len,pmap->len * C, numsite);
	      fflush(stdout);
	    }
	  }
	  pmap->allfree();
	  pmap->numsite[0] = -1;
	  break;
	}

	int c = color -1;

	if(numsite > (pmap->blockmem ? pmap->blockmem-2 : pmap->numsite[c])){
	  if(VERB>=2){
            #pragma omp critical
	    {
	      printf("pmap->id=%lld,mapid=%d: increasing site[] memory from numsite[%d]=%d to %d (blockmem=%d)\n",
	        pmap->id,mapid,c,(pmap->blockmem ? pmap->blockmem-2 : pmap->numsite[color-1]), numsite, pmap->blockmem);
	      fflush(stdout);
	    }
	  }
	  if(!pmap->blockmem){
	    delete [] pmap->site[c];
	    if(MapTrueSite || Tracksites){
	      delete [] pmap->truesiteL[c];
	      delete [] pmap->truesiteR[c];
	    }
	    if(MapSNR) delete [] pmap->SNR[c];
	    if(MapIntensity) delete [] pmap->Intensity[c];
	    if(MapStitched) delete [] pmap->stitch[c];
	    if(MapPSFWidth) delete [] pmap->PSFWidth[c];
	    if(MapImageLoc){
	      delete [] pmap->ImageFOV[c];
	      delete [] pmap->ImageX[c];
	      delete [] pmap->ImageY[c];
	      if(DEBUG>=2){
		pmap->ImageFOV[c] = NULL;
		pmap->ImageX[c] = NULL;
		pmap->ImageY[c] = NULL;
	      }
            }
	  }
	  pmap->site[c] = new FLOAT[numsite+2];
	  if(VERB>=2){
            #pragma omp critical
	    {
	      printf("re-allocating Map[%d]->site[%d]= %p:id=%lld\n",mapid,c,pmap->site[c],pmap->id);
	      fflush(stdout);
	    }
	  }

	  if(DEBUG>=2)
	    for(int i = 0; i < numsite+2; i++)
	      pmap->site[c][i] = aNaN;
	  if(MapTrueSite || Tracksites){
	    pmap->truesiteL[c] = new long long[numsite+1];
	    pmap->truesiteR[c] = new long long[numsite+1];
	    if(DEBUG >= 2)
	      for(int i = 0; i <= numsite; i++)	      
		pmap->truesiteL[c][i] = pmap->truesiteR[c][i] = -1;
	  }
	  if(MapSNR){
	    pmap->SNR[c] = new FLOAT[numsite+1];
	    if(DEBUG>=2)
	      for(int i = 0; i < numsite+1; i++)
		pmap->SNR[c][i] = aNaN;
	  }
	  if(MapIntensity){
	    pmap->Intensity[c] = new FLOAT[numsite+1];
	    if(DEBUG>=2)
	      for(int i = 0; i < numsite+1; i++)
		pmap->Intensity[c][i] = aNaN;
	  }
	  if(MapStitched){
	    pmap->stitch[c] = new int[numsite+2];
	    if(DEBUG>=2)
	      for(int i = 0; i < numsite+2; i++)
		pmap->stitch[c][i] = -1;
	  }
	  if(MapPSFWidth){
	    pmap->PSFWidth[c] = new FLOAT[numsite+1];
	    if(DEBUG>=2)
	      for(int i = 0; i < numsite+1; i++)
		pmap->PSFWidth[c][i] = aNaN;
	  }
	  if(MapImageLoc){
	    pmap->ImageFOV[c] = new int[numsite+2];
	    pmap->ImageX[c] = new double[numsite+2];
	    pmap->ImageY[c] = new double[numsite+2];
	    if(DEBUG>=2)
	      for(int i = 0; i < numsite+1; i++){
		pmap->ImageFOV[c][i] = -1;
		pmap->ImageX[c][i] = aNaN;
		pmap->ImageY[c][i] = aNaN;
	      }
	  }

	  if(c + 1 < colors) { /* force subsequent allocation for other colors as well */
	    for(int nc = c + 1; nc < colors; nc++){
	      pmap->site[nc] = NULL;
	      if(MapTrueSite || Tracksites) pmap->truesiteL[nc] = pmap->truesiteR[nc] = NULL;
	      if(MapSNR) pmap->SNR[nc] = NULL;
	      if(MapIntensity) pmap->Intensity[nc] = NULL;
	      if(MapStitched) pmap->stitch[nc] = NULL;
	      if(MapPSFWidth) pmap->PSFWidth[nc] = NULL;
	      if(MapImageLoc){
		pmap->ImageFOV[nc] = NULL;
		pmap->ImageX[nc] = NULL;
		pmap->ImageY[nc] = NULL;
	      }
	      pmap->numsite[nc] = -1;
	    }
	  }
	  if(c > 0 && pmap->blockmem) {/* re-allocate non-block memory for previous colors */
	    for(int pc = 0; pc < c; pc++){
	      int numsite = pmap->numsite[pc];
	      if(DEBUG) assert(numsite+2 <= pmap->blockmem);
	      FLOAT *nsite = new FLOAT[numsite + 2];
	      memcpy(nsite, pmap->site[pc], (numsite+2)*sizeof(FLOAT));
	      if(VERB>=2){
		#pragma omp critical
		{
		  printf("re-allocating Map[%d]->site[%d]= %p:id=%lld\n",mapid,pc,nsite,pmap->id);
		  fflush(stdout);
		}
	      }
	      pmap->site[pc] = nsite;
	      if(MapTrueSite || Tracksites){
		long long *nTrueSite = new long long[numsite+1];
		if(MapTrueSite)
		  memcpy(nTrueSite, pmap->truesiteL[pc], (numsite + 1)*sizeof(long long));
		pmap->truesiteL[pc] = nTrueSite;
		
		nTrueSite = new long long[numsite+1];
		if(MapTrueSite)
		  memcpy(nTrueSite, pmap->truesiteR[pc], (numsite + 1) * sizeof(long long));
		pmap->truesiteR[pc] = nTrueSite;
	      }
	      if(MapSNR){
	        double *nSNR = new double[numsite + 1];
		memcpy(nSNR, pmap->SNR[pc], (numsite + 1)*sizeof(double));
		pmap->SNR[pc] = nSNR;
	      }
	      if(MapIntensity){
	        double *nIntensity = new double[numsite + 1];
		memcpy(nIntensity, pmap->Intensity[pc], (numsite + 1)*sizeof(double));
		pmap->Intensity[pc] = nIntensity;
    	      }
	      if(MapStitched){
	        int *nstitch = new int[numsite + 2];
		memcpy(nstitch, pmap->stitch[pc], (numsite + 2)*sizeof(int));
		pmap->stitch[pc] = nstitch;
    	      }
	      if(MapPSFWidth){
		double *nPSFWidth = new double[numsite + 1];
		memcpy(nPSFWidth, pmap->PSFWidth[pc], (numsite + 1)*sizeof(double));
		pmap->PSFWidth[pc] = nPSFWidth;
	      }
	      if(MapImageLoc){
		int *nImageFOV = new int[numsite + 2];
		memcpy(nImageFOV, pmap->ImageFOV[pc], (numsite + 2)*sizeof(int));
		pmap->ImageFOV[pc] = nImageFOV;

		double *nImageX = new double[numsite + 2];
		memcpy(nImageX, pmap->ImageX[pc], (numsite + 2) * sizeof(double));
		pmap->ImageX[pc] = nImageX;

		double *nImageY = new double[numsite + 2];
		memcpy(nImageY, pmap->ImageY[pc], (numsite + 2) * sizeof(double));
		pmap->ImageY[pc] = nImageY;
	      }
  	    }
	  }
	  pmap->blockmem = 0;
	}
	pmap->numsite[c] = numsite;
	pmap->site[c][0] = 0.0;
	for(register int i = 0; i < numsite;i++){
	  if(DEBUG>=2) assert(isnan(pmap->site[c][i+1]));
	  pmap->site[c][i+1] = Location[i];
	}
	if(DEBUG>=2) assert(isnan(pmap->site[c][numsite+1]));
	// WAS43	pmap->site[c][numsite+1] = pmap->len;
	if(VERB>=2 || (BIAS_TRACE && pmap->id==BIAS_TRACE_ID)){
	  printf("pmap->id=%lld:color=%d,c=%d:pmap->len= %0.6f, numsite=%d from line %d in %s\n%s\n",pmap->id,color,c,pmap->len, numsite,linecnt,filename,buf);
	  for(int i = 0; i <= numsite+1; i++)
	    printf("pmap->site[c][%d] = %0.6f\n",i,pmap->site[c][i]);
	  fflush(stdout);
	}
	
	totsites += numsite;

      }// for(color = 1; color <= colors; color++)
      
      if(BNX_FASTFILTER && color <= colors)
	continue;

      for(int c = 0; c < colors; c++){
	pmap->site[c][pmap->numsite[c]+1] = pmap->len;// NEW43

	if(VERB>=2 || (BIAS_TRACE && pmap->id==BIAS_TRACE_ID)){
	  printf("pmap->id=%lld:c=%d:pmap->len= %0.6f, numsite=%d\n",pmap->id,c,pmap->len, pmap->numsite[c]);
	  for(int i = 0; i <= pmap->numsite[c] + 1; i++)
	    printf("pmap->site[c][%d] = %0.6f\n",i,pmap->site[c][i]);
	  fflush(stdout);
	}
      }

      if(Tracksites && !MapTrueSite) {/* generate QX10/QX20 values to track original index values */
	if(colors == 1)
	  for(int i = 1; i <= pmap->numsite[0]; i++)
	    pmap->truesiteL[0][i] = pmap->truesiteR[0][i] = i;
	else {// first convert 2-color map into an interleaved format (code copied from output_xmap.cpp)
	  int *NN = pmap->numsite;
	  FLOAT **YY = pmap->site;
	  Cinterleaved *Ysite = new Cinterleaved[NN[0]+NN[1]+2]; // array of (c,index), corresponding to YY[c][index] and in ascending order of this value
	  int *Yindex[2];
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
	  if(DEBUG) assert(Ysites == NN[0]+NN[1]);
	  // Now Yindex[c][i] is translation from original index i of color c to interleaved index

	  for(int c = 0; c < colors; c++)
	    for(int i = 1; i <= NN[c]; i++)
	      pmap->truesiteL[c][i] = pmap->truesiteR[c][i] = Yindex[c][i];

	  delete [] Ysite;
	  delete [] Yindex[0];
	}
      }

      if(NumSite >= 0 && totsites != NumSite){
        #pragma omp critical
	{
	  if(colors==1)
	    printf("ERROR:NumberofLabels=%d on line %d does not match number of labels=%d on line %d of %s\n%s\n",NumSite,NumSiteLine,totsites,linecnt,filename,buf);
	  else
	    printf("ERROR:NumberofLabels=%d on line %d does not match total number of labels=%d on lines %d to %d of %s\n%s\n",NumSite,NumSiteLine,totsites,NumSiteLine+1,linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
      }

      int subline = colors;

      /* check if we need to read in QX information */

      if(MapTrueSite){
	if(DEBUG) assert(colors <= 2);
	if(VERB>=2){
	   printf("pmap->id=%lld,mapid=%d:numsite[0]=%d:Looking for TrueSite on line %d\n",
		  pmap->id,mapid,pmap->numsite[0],lines[line+subline].linecnt);
	   fflush(stdout);
	}
	for(int color = 1; color <= colors; color++){      
	  linecnt = lines[line+subline+color].linecnt;
	  buf = lines[line+subline+color].buf;
	  int c = color -1;

	  /* check for TrueSite of color=c */
	  const char *Qcode = VQ_TRUESITE[c][BNXVersion];
	  size_t qlen = strlen(Qcode);
	  if(strncmp(buf,Qcode,qlen)){
            #pragma omp critical
	    {
	      printf("TrueSite for color=%d does not start with %s on line %d of %s\n%s\n",
		     color,Qcode,linecnt,filename,buf);
	      fflush(stdout);exit(1);
	    }
	  }

	  if(VERB>=2){
	    printf("pmap->id=%lld,mapid=%d:allocating truesite*[c=%d], numsite[c]=%d\n",pmap->id,mapid,c,pmap->numsite[c]);
	    fflush(stdout);
	  }
	  if(!pmap->blockmem){
	    delete [] pmap->truesiteL[c];
	    delete [] pmap->truesiteR[c];
	    pmap->truesiteL[c] = new long long[pmap->numsite[c]+1];
	    pmap->truesiteR[c] = new long long[pmap->numsite[c]+1];
	    if(DEBUG>=2)
	      for(int i = 0; i <= pmap->numsite[c]; i++)
		pmap->truesiteL[c][i] = pmap->truesiteR[c][i] = -1;
	  }

	  if(DEBUG>=2 && !(pmap->truesiteL[c][0] < 0 && pmap->truesiteR[c][0] < 0)){
	    #pragma omp critical
	    {
	      printf("pmap->id=%lld,mapid=%d:c=%d:numsite[c]=%d,pmap->truesite[c][0]=%lld,%lld,blockmem=%d,MINSITE=%d\n",
		     pmap->id,mapid,c,pmap->numsite[c],pmap->truesiteL[c][0],pmap->truesiteR[c][0],pmap->blockmem,MINSITE);
	      printf("pmap->site[c] = %p, pmap->truesiteL[c]=%p, pmap->truesiteR[c]=%p\n",pmap->site[c],pmap->truesiteL[c],pmap->truesiteR[c]); 
	      for(int i = 0; i <= pmap->numsite[c]+1; i++)
		printf("   pmap->site[c][%d]=%0.5f\n",i, pmap->site[c][i]);
	      for(int i = 0; i < pmap->numsite[c]+1; i++)
		printf("   pmap->truesite[c][%d]=%lld,%lld\n",i, pmap->truesiteL[c][i],pmap->truesiteR[c][i]);
	      fflush(stdout);
	      assert(pmap->truesiteL[c][0] < 0 && pmap->truesiteR[c][0]);
	    }
	  }
	  pmap->truesiteL[c][0] = pmap->truesiteR[c][0] = 0;

	  int site = 1;
	  for(pt = &buf[qlen];*pt; pt = qt){
	    long long truesite = strtoll(pt,&qt,10);
	    if(pt==qt){
	      while(*pt && isspace(*pt))
		pt++;
	      if(!*pt)/* End of Line */
		break;
              #pragma omp critical
	      {
		printf("Invalid %s TrueSite value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(!(*qt == 0 || isspace(*qt))){
	      while(*pt && isspace(*pt))
		pt++;
              #pragma omp critical
	      {
		printf("Invalid %s TrueSite value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(site > pmap->numsite[c]){
              #pragma omp critical
	      {
		printf("Too many %s TrueSite values on line %d in %s: expected %d TrueSite values (numsites=%d):\n%s\n",
		       Qcode, linecnt,filename,pmap->numsite[c]*2,pmap->numsite[c],buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(site > pmap->numsite[c])
	      break;
	    if(DEBUG>=2) assert(pmap->truesiteL[c][site] < 0);
	    pmap->truesiteL[c][site] = truesite;

	    pt = qt;
	    truesite = strtoll(pt,&qt,10);
	    if(pt==qt){
	      while(*pt && isspace(*pt))
		pt++;
	      if(!*pt)/* End of Line */
		break;
              #pragma omp critical
	      {
		printf("Invalid %s TrueSite value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(!(*qt == 0 || isspace(*qt))){
	      while(*pt && isspace(*pt))
		pt++;
              #pragma omp critical
	      {
		printf("Invalid %s TrueSite value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(DEBUG>=2) assert(pmap->truesiteR[c][site] < 0);
	    pmap->truesiteR[c][site++] = truesite;
	  }
	  site--;
	  if(site != pmap->numsite[c]){
	    #pragma omp critical
	    {
	      printf("Incorrect number of %s SNR values on line %d in %s: expected at least %d SNR values, got %d values (numsites=%d):\n%s\n\n",
		     Qcode, linecnt,filename,pmap->numsite[c],pmap->numsite[c],site,buf);
	      fflush(stdout);exit(1);
	    }
	  }
	}
	subline += colors;
      }

      if(MapSNR){
	if(DEBUG) assert(colors <= 2);
	if(VERB>=2){
	   printf("pmap->id=%lld,mapid=%d:numsite[0]=%d:Looking for SNR on line %d\n",
		  pmap->id,mapid,pmap->numsite[0],lines[line+subline].linecnt);
	   fflush(stdout);
	}
	for(int color = 1; color <= colors; color++){      
	  linecnt = lines[line+subline+color].linecnt;
	  buf = lines[line+subline+color].buf;
	  int c = color -1;

	  /* check for SNR of color=c */
	  const char *Qcode = VQ_SNR[c][BNXVersion];
	  size_t qlen = strlen(Qcode);
	  if(strncmp(buf,Qcode,qlen)){
            #pragma omp critical
	    {
	      printf("SNR for color=%d does not start with %s on line %d of %s\n%s\n",
		     color,Qcode,linecnt,filename,buf);
	      fflush(stdout);exit(1);
	    }
	  }

	  if(VERB>=2){
	    printf("pmap->id=%lld,mapid=%d:allocating SNR[c=%d], numsite[c]=%d\n",pmap->id,mapid,c,pmap->numsite[c]);
	    fflush(stdout);
	  }
	  if(!pmap->blockmem){
	    delete [] pmap->SNR[c];
	    pmap->SNR[c] = new double[pmap->numsite[c]+1];
	    if(DEBUG>=2)
	      for(int i = 0; i < pmap->numsite[c]+1; i++)
		pmap->SNR[c][i] = aNaN;
	  }
	  if(DEBUG>=2 && !isnan(pmap->SNR[c][0])){
	    #pragma omp critical
	    {
	      printf("pmap->id=%lld,mapid=%d:c=%d:numsite[c]=%d,pmap->SNR[c][0]=%0.8e,blockmem=%d,MINSITE=%d\n",
		     pmap->id,mapid,c,pmap->numsite[c],pmap->SNR[c][0],pmap->blockmem,MINSITE);
	      printf("pmap->site[c] = %p, pmap->SNR[c]=%p\n",pmap->site[c],pmap->SNR[c]); 
	      for(int i = 0; i <= pmap->numsite[c]+1; i++)
		printf("   pmap->site[c][%d]=%0.5f\n",i, pmap->site[c][i]);
	      for(int i = 0; i < pmap->numsite[c]+1; i++)
		printf("   pmap->SNR[c][%d]=%0.8e\n",i, pmap->SNR[c][i]);
	      fflush(stdout);
	      assert(isnan(pmap->SNR[c][0]));
	    }
	  }
	  pmap->SNR[c][0] = 0.0;

	  int site = 1;
	  for(pt = &buf[qlen];*pt; pt = qt){
	    double SNR = strtod(pt,&qt);
	    if(pt==qt){
	      while(*pt && isspace(*pt))
		pt++;
	      if(!*pt)/* End of Line */
		break;
              #pragma omp critical
	      {
		printf("Invalid %s SNR value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(!(*qt == 0 || isspace(*qt))){
	      while(*pt && isspace(*pt))
		pt++;
              #pragma omp critical
	      {
		printf("Invalid %s SNR value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(site > pmap->numsite[c]+(1-BNXVersion)){
              #pragma omp critical
	      {
		printf("Too many %s SNR values on line %d in %s: expected %d SNR values (numsites=%d):\n%s\n",
		       Qcode, linecnt,filename,pmap->numsite[c]+1-BNXVersion,pmap->numsite[c],buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(site > pmap->numsite[c])
	      break;
	    if(DEBUG>=2) assert(isnan(pmap->SNR[c][site]));
	    pmap->SNR[c][site++] = max(0.001,SNR);
	  }
	  site--;
	  if(site != pmap->numsite[c]){
	    #pragma omp critical
	    {
	      printf("Incorrect number of %s SNR values on line %d in %s: expected at least %d SNR values, got %d values (numsites=%d):\n%s\n\n",
		     Qcode, linecnt,filename,pmap->numsite[c],pmap->numsite[c],site,buf);
	      fflush(stdout);exit(1);
	    }
	  }
	}
	subline += colors;
      }

      if(MapIntensity){
	//	if(DEBUG) assert(BNXVersion==1);
	if(DEBUG>=2) assert(colors <= 2);
	for(int color = 1; color <= colors; color++){      
	  linecnt = lines[line+subline+color].linecnt;
	  buf = lines[line+subline+color].buf;
	  int c = color -1;

	  /* check for Intensity of color=c */
	  const char *Qcode = VQ_INTENSITY[c][BNXVersion];
	  size_t qlen = strlen(Qcode);
	  if(strncmp(buf,Qcode,qlen)){
            #pragma omp critical
	    {
	      printf("Intensity for color=%d does not start with %s on line %d of %s\n%s\n",
		     color,Qcode,linecnt,filename,buf);
	      fflush(stdout);exit(1);
	    }
	  }

	  if(VERB>=2){
	    printf("pmap->id=%lld,mapid=%d:allocating Intensity[c=%d], numsite[c]=%d\n",pmap->id,mapid,c,pmap->numsite[c]);
	    fflush(stdout);
	  }
	  if(!pmap->blockmem){
	    delete [] pmap->Intensity[c];
	    pmap->Intensity[c] = new double[pmap->numsite[c]+1];
	    if(DEBUG>=2)
	      for(int i = 0; i < pmap->numsite[c] + 1; i++)
		pmap->Intensity[c][i] = aNaN;
	  }
	  if(DEBUG>=2) assert(isnan(pmap->Intensity[c][0]));
	  pmap->Intensity[c][0] = 0.0;

	  int site = 1;
	  for(pt = &buf[qlen];*pt; pt = qt){
	    double Intensity = strtod(pt,&qt);
	    if(DEBUG) assert(qt >= pt);
	    if(pt==qt){
	      while(*pt && isspace(*pt))
		pt++;
	      if(!*pt)/* End of Line */
		break;
              #pragma omp critical
	      {
		printf("Invalid %s Intensity value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(!(*qt == 0 || isspace(*qt))){
	      while(*pt && isspace(*pt))
		pt++;
              #pragma omp critical
	      {
		printf("Invalid %s Intensity value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(site > pmap->numsite[c]+(1-BNXVersion)){
              #pragma omp critical
	      {
		printf("Too many %s Intensity values at char %ld..%ld on line %d in %s: expected %d Intensity values (numsites=%d):\n%s\n",
		       Qcode, pt-buf+1, qt-buf, linecnt, filename,pmap->numsite[c]+1-BNXVersion,pmap->numsite[c],buf);
		printf("  Extra Intensity=%0.6f\n",Intensity);
		fflush(stdout);exit(1);
	      }
	    }
	    if(site > pmap->numsite[c])
	      break;
	    if(DEBUG>=2) assert(isnan(pmap->Intensity[c][site]));
	    pmap->Intensity[c][site++] = Intensity;
	  }
	  site--;

	  if(site != pmap->numsite[c]){
	    #pragma omp critical
	    {
	      printf("Incorrect number of %s Intensity values on line %d in %s: expected %d Intensity values, got %d values:\n%s\n",
		     Qcode, linecnt,filename,pmap->numsite[c],site,buf);
	      for(int t = 1; t <= site; t++)
		printf("site[%d]:Intensity=%0.6f\n",t,pmap->Intensity[c][t]);
	      fflush(stdout);exit(1);
	    }
	  }
	}
	subline += colors;
      }

      if(MapStitched){
	if(DEBUG>=2) assert(colors <= 2);
	for(int color = 1; color <= colors; color++){      
	  int c = color -1;
	  linecnt = lines[line+subline+color].linecnt;
	  buf = lines[line+subline+color].buf;
	
	  /* check for stitch of color=c */
	  const char *Qcode = VQ_STITCH[c][BNXVersion];
	  size_t qlen = strlen(Qcode);
	  if(strncmp(buf,Qcode,qlen)){
            #pragma omp critical
	    {
	      printf("Stitch for color=%d does not start with %s on line %d of %s\n%s\n",
		     color,Qcode,linecnt,filename,buf);
	      fflush(stdout);
	      assert(0);
	      fflush(stdout);exit(1);
	    }
	  }
	  if(!pmap->blockmem){
	    delete [] pmap->stitch[c];
	    pmap->stitch[c] = new int[pmap->numsite[c]+2];
	  }
	  pmap->stitch[c][0] = 0;

	  int site = 1;
	  for(pt = &buf[qlen];*pt; pt = qt){
	    int stitch = strtol(pt,&qt,10);
	    if(pt==qt){
	      while(*pt && isspace(*pt))
		pt++;
	      if(!*pt)/* End of Line */
		break;
              #pragma omp critical
	      {
		printf("Invalid %s stitch value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(!(*qt == 0 || isspace(*qt))){
	      while(*pt && isspace(*pt))
		pt++;
              #pragma omp critical
	      {
		printf("Invalid %s stitch value at char %ld of line %d in %s:'%s'\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(site > pmap->numsite[c]+1){
              #pragma omp critical
	      {
		printf("Too many %s stitch values on line %d in %s: expected %d stitch values\n",
		       Qcode,linecnt,filename,pmap->numsite[c]+1);
		fflush(stdout);exit(1);
	      }
	    }
	    pmap->stitch[c][site++] = stitch;
	  }
	  site--;
	  if(site != pmap->numsite[c]+1){
            #pragma omp critical
	    {
	      printf("Incorrect number of %s stitch values on line %d in %s: expected %d stitch values, got %d values:\n%s\n",
		     Qcode,linecnt,filename,pmap->numsite[c]+1,site,buf);
	      fflush(stdout);exit(1);
	    }
	  }
	}
	subline += colors;
      }

      if(MapStitchLoc){
	if(DEBUG>=2) assert(colors <= 2);
	for(int color = 1; color <= colors; color++){
	  int c = color - 1;
	  linecnt = lines[line+subline+color].linecnt;
	  buf = lines[line+subline+color].buf;
	
	  /* check for stitchLoc of color=c */
	  const char *Qcode = VQ_STITCHLOC[c][BNXVersion];
	  size_t qlen = strlen(Qcode);
	  if(strncmp(buf,Qcode,qlen)){
            #pragma omp critical
	    {
	      printf("StitchLoc for color=%d does not start with %s on line %d of %s\n%s\n",
		     color,Qcode,linecnt,filename,buf);
	      //	      printf("line=%d(linecnt=%d),subline=%d\n",line,lines[line].linecnt,subline);
	      fflush(stdout);exit(1);
	    }
	  }
	
	  /* scan over the input line to find the Stitch Locations */
	  int numstitch = 0;
	  for(pt = &buf[qlen]; *pt ; pt = qt){
	    if(numstitch >= MAXSITES){
	      #pragma omp critical
	      {
		printf("More than %d stitch sites on line %d of file %s\n",numstitch,linecnt,filename);
		fflush(stdout);exit(1);
	      }
	    }
	    Location[numstitch] = strtod(pt,&qt);
	    if(pt==qt)
	      break;
	    Location[numstitch] *= 0.001;

	    if(0&&/* HERE */(Location[numstitch] <= 0.0 || Location[numstitch] > pmap->len)){
              #pragma omp critical
	      {
	        printf("ERROR: Out of bounds stitch location=%0.8e on line %d of file %s (map length was %0.1f):\n%s\n",
		  Location[numstitch]*1000.0,linecnt,filename,pmap->len*1000.0,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    /*	 HERE if(Location[numstitch] > pmap->len)
	      Location[numstitch] = pmap->len;
	    if(Location[numstitch] < 0.0)
	    Location[numstitch] = 0.0; */
	    numstitch++;
	  }
	
          if(VERB>=2 && pmap->id == 1){
	    #pragma omp critical
	    {
	      printf("id=%lld:c=%d,numsites=%d,numstitch=%d (line %d of %s):\n",pmap->id,c+1,pmap->numsite[c],numstitch,linecnt,filename);
	      for(int i = 0; i < numstitch; i++)
		printf("stitch[%d]= %0.6f\n",i,Location[i]);
	      printf("%s\n",buf);
	      fflush(stdout);
	    }
	  }

	  pmap->NumStitch[c] = numstitch;
	  delete [] pmap->stitchLocation[c];
	  pmap->stitchLocation[c] = new double[numstitch];

	  for(register int i = 0; i < numstitch; i++)
	    pmap->stitchLocation[c][i] = Location[i];

	}/* color= 1 .. colors */

	subline += colors;
      }// MapStitchLoc

      if(MapPSFWidth){
	if(DEBUG) assert(colors <= 2);
	if(VERB>=2){
	   printf("pmap->id=%lld,mapid=%d:numsite[0]=%d:Looking for MapPSFWidth on line %d\n",
		  pmap->id,mapid,pmap->numsite[0],lines[line+subline].linecnt);
	   fflush(stdout);
	}
	for(int color = 1; color <= colors; color++){      
	  linecnt = lines[line+subline+color].linecnt;
	  buf = lines[line+subline+color].buf;
	  int c = color -1;

	  /* check for PSFWidth for color=c */
	  const char *Qcode = VQ_PSFWIDTH[c][BNXVersion];
	  size_t qlen = strlen(Qcode);
	  if(strncmp(buf,Qcode,qlen)){
            #pragma omp critical
	    {
	      printf("PSFWidth for color=%d does not start with %s on line %d of %s\n%s\n",
		     color,Qcode,linecnt,filename,buf);
	      fflush(stdout);exit(1);
	    }
	  }

	  if(VERB>=2){
	    printf("pmap->id=%lld,mapid=%d:allocating PSFWidth[c=%d], numsite[c]=%d\n",pmap->id,mapid,c,pmap->numsite[c]);
	    fflush(stdout);
	  }
	  if(!pmap->blockmem){
	    delete [] pmap->PSFWidth[c];
	    pmap->PSFWidth[c] = new double[pmap->numsite[c]+1];
	  }
	  if(DEBUG>=2)
	    for(int i = 0; i < pmap->numsite[c]+1; i++)
	      pmap->PSFWidth[c][i] = aNaN;
	  if(DEBUG>=2 && !isnan(pmap->PSFWidth[c][0])){
	    #pragma omp critical
	    {
	      printf("pmap->id=%lld,mapid=%d:c=%d:numsite[c]=%d,pmap->PSFWidth[c][0]=%0.8e,blockmem=%d,MINSITE=%d\n",
		     pmap->id,mapid,c,pmap->numsite[c],pmap->PSFWidth[c][0],pmap->blockmem,MINSITE);
	      printf("pmap->site[c] = %p, pmap->PSFWidth[c]=%p\n",pmap->site[c],pmap->PSFWidth[c]); 
	      for(int i = 0; i <= pmap->numsite[c]+1; i++)
		printf("   pmap->site[c][%d]=%0.5f\n",i, pmap->site[c][i]);
	      for(int i = 0; i < pmap->numsite[c]+1; i++)
		printf("   pmap->PSFWidth[c][%d]=%0.8e\n",i, pmap->PSFWidth[c][i]);
	      fflush(stdout);
	      assert(isnan(pmap->PSFWidth[c][0]));
	    }
	  }
	  pmap->PSFWidth[c][0] = 0.0;

	  int site = 1;
	  for(pt = &buf[qlen];*pt; pt = qt){
	    double PSFWidth = strtod(pt,&qt);
	    if(pt==qt){
	      while(*pt && isspace(*pt))
		pt++;
	      if(!*pt)/* End of Line */
		break;
              #pragma omp critical
	      {
		printf("Invalid %s PSFWidth value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(!(*qt == 0 || isspace(*qt))){
	      while(*pt && isspace(*pt))
		pt++;
              #pragma omp critical
	      {
		printf("Invalid %s PSFWidth value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(site > pmap->numsite[c]+(1-BNXVersion)){
              #pragma omp critical
	      {
		printf("Too many %s PSFWidth values on line %d in %s: expected %d PSFWidth values (numsites=%d):\n%s\n",
		       Qcode, linecnt,filename,pmap->numsite[c]+1-BNXVersion,pmap->numsite[c],buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(site > pmap->numsite[c])
	      break;
	    if(DEBUG>=2) assert(isnan(pmap->PSFWidth[c][site]));
	    pmap->PSFWidth[c][site++] = max(0.001,PSFWidth);
	  }
	  site--;
	  if(site != pmap->numsite[c]){
	    #pragma omp critical
	    {
	      printf("Incorrect number of %s PSFWidth values on line %d in %s: expected at least %d PSFWidth values, got %d values (numsites=%d):\n%s\n\n",
		     Qcode, linecnt,filename,pmap->numsite[c],pmap->numsite[c],site,buf);
	      fflush(stdout);exit(1);
	    }
	  }
	}
	subline += colors;
      } // MapPSFWidth;

      if(MapImageLoc){
	if(DEBUG) assert(colors <= 2);
	if(VERB>=2){
	  printf("pmap->id=%lld,mapid=%d:numsite[0]=%d:Looking for MapImageLoc on line %d\n",
		 pmap->id,mapid,pmap->numsite[0],lines[line+subline].linecnt);
	  fflush(stdout);
	}
	for(int color = 1; color <= colors; color++){
	  linecnt = lines[line+subline+color].linecnt;
	  buf = lines[line+subline+color].buf;
	  int c = color -1;
	  int numsite = pmap->numsite[c];
	  
	  /* check for ImageLoc for color=c */
	  const char *Qcode = VQ_IMAGELOC[c][BNXVersion];
	  size_t qlen = strlen(Qcode);
	  if(strncmp(buf,Qcode,qlen)){
            #pragma omp critical
	    {
	      printf("ImageLoc for color=%d does not start with %s on line %d of %s\n%s\n",
		     color,Qcode,linecnt,filename,buf);
	      fflush(stdout);exit(1);
	    }
	  }
	  
	  if(VERB>=2){
	    printf("pmap->id=%lld,mapid=%d:allocating ImageFOV[c],ImageX[c],ImageY[c] (c=%d], numsite[c]=%d\n",pmap->id,mapid,c,pmap->numsite[c]);
	    fflush(stdout);
	  }
	  if(!pmap->blockmem){
	    delete [] pmap->ImageFOV[c];
	    delete [] pmap->ImageX[c];
	    delete [] pmap->ImageY[c];
	    
	    pmap->ImageFOV[c] = new int[numsite+2];
	    pmap->ImageX[c] = new double[numsite+2];
	    pmap->ImageY[c] = new double[numsite+2];
	  }
	  if(DEBUG>=2){
	    for(int i = 0; i <= numsite+1; i++){
	      pmap->ImageFOV[c][i] = -1;
	      pmap->ImageX[c][i] = aNaN;
	      pmap->ImageY[c][i] = aNaN;
	    }
	  }
	  int site = 1;
	  for(pt = &buf[qlen]; *pt; pt = qt){
	    int FOV = strtol(pt,&qt,10);
	    if(pt==qt){
	      while(*pt && isspace(*pt))
		pt++;
	      if(!*pt)/* End of Line */
		break;
              #pragma omp critical
	      {
		printf("Invalid %s ImageFOV value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(!(*qt == 0 || isspace(*qt))){
	      while(*pt && isspace(*pt))
		pt++;
              #pragma omp critical
	      {
		printf("Invalid %s ImageFOV value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(site > numsite) {
	      #pragma omp critical
	      {
		printf("Too many %s ImageLoc values on line %d in %s: expcected %d sets of values (FOV,X,Y):\n%s\n",
		       Qcode, linecnt, filename, numsite, buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(DEBUG>=2) assert(pmap->ImageFOV[c][site] == -1);
	    if(DEBUG>=2) assert(FOV >= 0);
	    if(FOV > MaxImageFOV){
	      #pragma omp critical
	      {
		if(VERB>=2){
		  printf("ImageFOV=%d at char %ld on line %d in %s\n",FOV,pt-buf+1,linecnt,filename);
		  fflush(stdout);
		}
		MaxImageFOV = max(MaxImageFOV,FOV);
	      }
	    }
	    pmap->ImageFOV[c][site] = FOV;

	    pt = qt;
	    double X = strtod(pt,&qt);
	    if(pt==qt){
	      while(*pt && isspace(*pt))
		pt++;
	      if(!*pt)/* End of Line */
		break;
              #pragma omp critical
	      {
		printf("Invalid %s X coordinate value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(!(*qt == 0 || isspace(*qt))){
	      while(*pt && isspace(*pt))
		pt++;
              #pragma omp critical
	      {
		printf("Invalid %s X-coordinate value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(DEBUG>=2) assert(isnan(pmap->ImageX[c][site]));
	    pmap->ImageX[c][site] = X;

	    pt = qt;
	    double Y = strtod(pt,&qt);
	    if(pt==qt){
	      while(*pt && isspace(*pt))
		pt++;
	      if(!*pt)/* End of Line */
		break;
              #pragma omp critical
	      {
		printf("Invalid %s Y coordinate value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(!(*qt == 0 || isspace(*qt))){
	      while(*pt && isspace(*pt))
		pt++;
              #pragma omp critical
	      {
		printf("Invalid %s Y-coordinate value at char %ld of line %d in %s:%s\n%s\n",
		       Qcode,pt-buf+1,linecnt,filename,pt,buf);
		fflush(stdout);exit(1);
	      }
	    }
	    if(DEBUG>=2) assert(isnan(pmap->ImageY[c][site]));
	    pmap->ImageY[c][site] = Y;

	    site++;
	  }
	  if(--site != numsite){
	    #pragma omp critical
	    {
	      printf("Incorrect number of %s (FOV,X,Y) values on line %d in %s: expected %d sets of values, got %d sets:\n%s\n",
		     Qcode,linecnt,filename,numsite,site,buf);
	      fflush(stdout);exit(1);
	    }
	  }
	}
      }
    } // for
  } // parallel
  
  if(VERB>=2){
    printf("Finished parallel section:mtime=%0.6f(wall=%0.6f)\n",mtime(),wtime());
    fflush(stdout);
  }
  if(BNX_FASTFILTER && !ids){// NEW102
    int i,j;
    for(i = j = oldmaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      if(pmap->numsite[0] < 0)
	continue;

      if(j < i){
	Cmap *tmp = Map[j];
	pmap->mapid = j;
	Map[j] = pmap;
	if(tmp)
	  tmp->mapid = i;
	Map[i] = tmp;
      }
      j++;
    }
    fastfilterCnt += nummaps-j;

    if(VERB>=2 && j != nummaps){
      long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
      getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
      printf("Reduced input maps from %d to %d due to MinLen=%0.3f,MaxLen=%0.3f, MinSites=%d:",  nummaps-oldmaps, j-oldmaps, tMinLen,tMaxLen,MinSites);
      printf("VmRSS= %0.4f, VmHWM= %0.4f : CPU time= %0.6f, wall time= %0.6f\n",VmRSS * 1e-9, VmHWM * 1e-9, mtime(), wtime());      
      fflush(stdout);
    }
    nummaps = j;
  }
}

static int Qcode_warning = 0;

/* If ids != 0 : read in only maps that match the set of id values provided */
int input_bnx(int fileid, int &nummaps, int &maxmaps, Cmap **&Map, std::set<long long> *ids)
{
  if(VERB>=2){
    printf("Calling input_bnx:fileid=%d:time=%0.6f(%0.6f)\n",fileid,mtime(),wtime());
    //    printf("MapStitchLoc=%d, MapStitchLocOrig=%d\n",MapStitchLoc,MapStitchLocOrig);
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  if(!origPixelLen)
    origPixelLen = 0.500;/* default value for .bnx files */
  else if(fabs(origPixelLen - 0.500) > 0.001){
    printf("WARNING: Input BNX file %s assumes bpp = 500, differs from previous input file bpp=%0.6f\n", 
	   vfixx_filename[fileid], origPixelLen*1000.0);
    fflush(stdout);
  }

  RunDataIdCnt = 0;/* number of "Run Data" lines read from current file */
  if(RunDataIdMax < 32){
    delete [] RunDataIdmap;
    RunDataIdMax = 32;
    RunDataIdmap = new int[RunDataIdMax];/* mapping from RunData line id in current file to index in RunDataList[] */
  }

  //  fprintf(stderr, "sizeof(Cmap)=%lu sizeof(Cnanomap)=%lu\n", sizeof(Cmap), sizeof(Cnanomap));
  
  char *filename = vfixx_filename[fileid];
  register char *pt;
  char *qt;
  FILE *fp;

  if((fp = Fopen(filename,"r",NFSdelay))==NULL){
    printf("input_bnx:Failed to read input file %s (%d'th of %d files)\n",filename, fileid+1,num_files);
    fflush(stdout);exit(1);
  }

  if(DEBUG) assert(&Map == &Gmap);
  if(maptype == 1 && MapType == -1){
    printf("input_bnx:fileid=%d, file %s:previous input map type was CMAP which cannot be mixed with .bnx or .vfixx\n", fileid,filename);
    fflush(stdout);exit(1);
  }
  maptype = 0;

  register Cmap *pmap=0;
  int orignummaps = nummaps;

  int NColors= -1;

  int linecnt=1;

  assert(colors >= 1);

  for(int c = 0; c < MAXCOLOR; c++)
    if(!Nickase[c]){
      Nickase[c] = strdup((char*) "unknown");
      if(VERB>=2){
	printf("input_bnx(%s):Setting global Nickase[%d]=%p(%s)\n",filename,c,Nickase[c],Nickase[c]);
	fflush(stdout);
      }
    }

  size_t blen = 10*1024*1024;/* matches scif_io buffer size */
  char *fbuf = (char *)malloc(blen);
  if(!fbuf){
    printf("input_bnx: malloc(%llu) failed\n", (unsigned long long)blen);
    fflush(stdout);exit(1);
  }

  setbuffer(fp, fbuf, blen);
    
  int bnxChim = 0;/* If current file has header with Feature* colums */
  int bnxRunHeaderChim = 0;/* If current file has "rh" header with Feature* fields */
  int bnxFXY = 0;/* If current file has header with "Column .. EndY" */

  int numthreads = 1;
  #ifdef _OPENMP
  numthreads = omp_get_max_threads();
  numthreads = min(numthreads,MaxThreads);
  #endif


  if(!fileid)
    QXmismatch = 0;

  int MapTrueSite_start = MapTrueSite;
  int MapSNR_start = MapSNR;
  int MapStitched_start = MapStitched;
  int MapIntensity_start = MapIntensity;
  int MapStitchLoc_start = MapStitchLoc;
  int MapPSFWidth_start = MapPSFWidth;
  int MapImageLoc_start = MapImageLoc;

  MapTrueSite = 0;
  MapSNR = MapSNROrig;
  MapStitched = MapStitchedOrig;
  MapIntensity = MapIntensityOrig;
  MapStitchLoc = MapStitchLocOrig;
  MapPSFWidth = MapPSFWidthOrig;
  MapImageLoc = MapImageLocOrig;

  int lines_per_map = 1 + colors*(1 + MapTrueSite + MapSNR + MapStitched + MapIntensity + MapStitchLoc + MapPSFWidth + MapImageLoc);

  int maxlines = numthreads * MAPS_PER_THREAD * lines_per_map;
  int numlines = 0;
  int total_lines=0;
  Cline *lines = new Cline[maxlines];/* parse input lines in parallel section */

  for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++) {
    size_t len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      if(len >= LINESIZ-1)
	printf("Line too long: exceeds max length=%d (increase LINESIZE in constants.h) in input file %s on line %d:\n%s\n",LINESIZ-1,filename,linecnt,buf);
      else
	printf("Line not correctly terminated (last char = %c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      fflush(stdout);exit(1);
    }

    if(buf[0] == '#'){/* most comment lines are ignored, except "BNX File Version:","Label Channels:","Run Data" and (if present) "Min Label SNR" and "Min Molecule Length (Kb)" */
      char *key = (char *)"rh";
      if(!strncmp(&buf[1],key,2) && isspace(buf[3])){
	/* check if "FeatureMinSlope FeatureMaxSlope" is present */
	if((pt = strstr(&buf[1],"FeatureMinSlope"))){
	  if(!(qt = strstr(pt,"FeatureMaxSlope"))){
	    printf("ERROR:Missing \"FeatureMaxSlop\" in rh header line after \"FeatureMinSlope\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  bnxRunHeaderChim = 1;
	}
	if(BNXRunHeaderChim < 0){
	  if(VERB/* HERE >=2 */){
	    if(bnxRunHeaderChim)
	      printf("Found extra columns FeatureMinSlope FeatureMaxSlope in Run Data header of %s\n",filename);
	    fflush(stdout);
	  }
	  BNXRunHeaderChim = bnxRunHeaderChim;
	} else if(BNXRunHeaderChim != bnxRunHeaderChim){
	  if(BNXRunHeaderChim == 0)
	    printf("WARNING:Inconsistent BNX Run Data header at line %d of %s: Current file has \"FeatureMinSlope FeatureMaxSlope\", some previous BNX files did not:\n%s\n",linecnt,filename,buf);
	  else {
	    printf("WARNING:Inconsistent BNX Run Data header at line %d of %s: Previous BNX files had \"FeatureMinSlope FeatureMaxSlope\", current file does not:\n%s\n",linecnt,filename,buf);
	    BNXRunHeaderChim = 0;
	  }
	  fflush(stdout);
	}

	continue;
      }

      key = (char *)"0h";
      if(VERB>=2){
	printf("Checking BNX header line %d of %s:\n%s\n",linecnt,filename,buf);
	printf("\tstrstr(&buf[1],\"%s\")= %p, isspace(buf[4])= %d\n",key,strstr(&buf[1],key),isspace(buf[3]) ? 1 : 0);
	fflush(stdout);
      }

      if(!strncmp(&buf[1],key,2) && isspace(buf[3])){
	/* check if Molecule Header includes MapWT */
        int mapWTpresent = (strstr(&buf[1],"MapWT") ? 1 : 0);
	if(MapWTpresent < 0)
	  MapWTpresent = mapWTpresent;
	else if(MapWTpresent != mapWTpresent){
	  if(MapWTpresent)
	    printf("ERROR:Inconsistent BNX Molecule header at line %d of %s: Previous files had \"MapWT\", current file does not:\n%s\n",linecnt,filename,buf);
	  else
	    printf("ERROR:Inconsistent BNX Molecule header at line %d of %s: Current file had \"MapWT\", previous files did not:\n%s\n",linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}

	/* check if Molecule header includes "Column StartFOV .. EndY" */
	if(VERB>=2){
	  printf("Checking BNX Molecule header at line %d of %s for \"Column\"\n%s\n",linecnt,filename,buf);
	  fflush(stdout);
	}
	if((pt = strstr(&buf[1],"Column"))){
	  if(!(qt = strstr(pt,"StartFOV"))){
	    printf("ERROR:Missing \"StartFOV\" in header line after \"Column\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  if(!(pt = strstr(qt,"StartX"))){
	    printf("ERROR:Missing \"StartX\" in header line after \"StartFOV\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  if(!(qt = strstr(pt,"StartY"))){
	    printf("ERROR:Missing \"StartY\" in header line after \"StartX\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  if(!(pt = strstr(qt,"EndFOV"))){
	    printf("ERROR:Missing \"EndFOV\" in header line after \"StartY\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  if(!(qt = strstr(pt,"EndX"))){
	    printf("ERROR:Missing \"EndX\" in header line after \"EndFOV\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  if(!(pt = strstr(qt,"EndY"))){
	    printf("ERROR:Missing \"EndY\" in header line after \"EndX\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  bnxFXY = 1;
	}
	if(BNXheaderFXY < 0){
	  if(VERB/* HERE >=2 */){
	    if(bnxFXY)
	      printf("Found extra columns in BNX Molecule header of %s: Column ... EndY\n",filename);
	    else if(VERB>=2)
	      printf("No extra columns (Column ... EndY) in BNX Molecule header of %s\n",filename);
	    fflush(stdout);
	  }
	  BNXheaderFXY = bnxFXY;
	} else if(BNXheaderFXY != bnxFXY){
	  if(BNXheaderFXY == 0)
	    printf("WARNING:Inconsistent BNX Molecule header at line %d of %s: Current file has \"Column\" etc, some previous files did not:\n%s\n",linecnt,filename,buf);
	  else {
	    printf("WARNING:Inconsistent BNX Molecule header at line %d of %s: Previous files had \"Column\" etc, current file does not:\n%s\n",linecnt,filename,buf);
	    BNXheaderFXY = 0;
	  }
	  fflush(stdout);
	}

	/* check if Molecule header includes "FeatureMinPos ... FeatureMaxScaled" */
	if((pt = strstr(&buf[1],"FeatureMinPos"))){
	  if(!(qt = strstr(pt,"FeatureMinScore"))){
	    printf("ERROR:Missing \"FeatureMinScore\" in header line after \"FeatureMinPos\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  if(!(qt = strstr(pt,"FeatureMaxScore"))){
	    printf("ERROR:Missing \"FeatureMaxScore\" in header line after \"FeatureMinPos\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  if(!(qt = strstr(pt,"FeatureMinScaled"))){
	    printf("ERROR:Missing \"FeatureMinScaled\" in header line after \"FeatureMinPos\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  if(!(qt = strstr(pt,"FeatureMaxScaled")) && !(qt = strstr(pt,"FeatuerMaxScaled"))){
	    printf("ERROR:Missing \"FeatureMaxScaled\" in header line after \"FeatureMinPos\" on line %d of input file %s:\n%s\n",linecnt,filename,buf);
	    fflush(stdout);exit(1);
	  }
	  bnxChim = 1;
	}

        if(BNXheaderChim < 0){
	  if(VERB/* HERE >=2 */){
	    if(bnxChim)
	      printf("Found extra columns in BNX Molecule header of %s: FeatureMinPos ... FeatureMaxScaled\n",filename);
	    fflush(stdout);
	  }
	  BNXheaderChim = bnxChim;
	} else if(BNXheaderChim != bnxChim){
	  if(BNXheaderChim == 0)
	    printf("WARNING:Inconsistent BNX Molecule header at line %d of %s: Current file has \"FeatureMinPos ... FeatureMaxScaled\", some previous files did not:\n%s\n",linecnt,filename,buf);
	  else { // BNXheaderChim==1
	    printf("WARNING:Inconsistent BNX Molecule header at line %d of %s: Previous files had \"FeatureMinPos ... FeatureMaxScaled\", current file does not:\n%s\n",linecnt,filename,buf);

	    BNXheaderChim = 0;
	  }
	  fflush(stdout);
	}

	continue;
      }

      key = (char *)"BNX File Version:";
      if((pt=strstr(&buf[1],key))){
	pt += strlen(key);
	double Version = strtod(pt,&qt);
	if(pt == qt  || !(*qt == 0 || isspace(*qt))){
	  printf("BNX File Version on line %d of %s not recognized:\n%s\n",linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	int majorVersion = (int)floor(Version);
	int minorVersion = ((int)floor(Version*10.0))%10;
	
	if(majorVersion <= 0){/* Version 0.x */
	  if(BNXVersion > 0){
	    printf("Inconsistent BNX File Versions in %s line %d: previous files had version %0.1f vs current file:\n%s\n",filename,linecnt, (BNXVersion==0 ? 0.1 : 1.0), buf);
	    fflush(stdout);exit(1);
	  }
	  BNXVersion = 0;
	} else if(majorVersion==1){/* Version 1.x */
	  if(BNXVersion == 0){
	    printf("Inconsistent BNX File Versions in %s line %d: previous files had version %0.1f vs current file:\n%s\n",filename,linecnt, (BNXVersion==0 ? 0.1 : 1.0), buf);
	    fflush(stdout);exit(1);
	  }
	  BNXVersion = 1;
	  if(BNXMinorVersion < 0)
	    BNXMinorVersion = minorVersion;
	  else {
	    if(BNXMinorVersion != minorVersion){
	      if(min(minorVersion,BNXMinorVersion) >= 2)
		printf("Inconsistent BNX Minor File Versions %d in %s line %d: previous files had minor version %d (upgrading to largest minor version)\n",minorVersion,filename,linecnt,BNXMinorVersion);
	      else
		printf("Inconsistent BNX Minor File Versions %d in %s line %d: previous files had minor version %d (downgrading to lowest minor version)\n",minorVersion,filename,linecnt,BNXMinorVersion);
	      fflush(stdout);
	    }
	    BNXMinorVersion = (min(minorVersion,BNXMinorVersion) >= 2) ? max(minorVersion, BNXMinorVersion) : min(minorVersion,BNXMinorVersion);
	  }
	} else {
	  printf("Unsupported BNX version number %0.1f (only 0.x and 1.x is supported)\n",Version);
	  fflush(stdout);exit(1);
	}

	continue;
      }

      key = (char *)"Label Channels:";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);
	NColors = strtol(pt,&qt,10);
	if(qt == pt || !(*qt == 0 || isspace(*qt)) || NColors <= 0 || NColors > 2){
	  printf("Line %d in %s with %s has invalid integer value (must be 1 or 2):\n%s\n",linecnt,filename, key,buf);
	  fflush(stdout);exit(1);
	}
	if(colors_set){
	  if(NColors != colors){
	    if(fileid > 0){
	      printf("Line %d in %s has %s %d which is different from value %d in previous input file %s:\n%s\n",linecnt,filename,key,NColors,colors,vfixx_filename[fileid-1],buf);
	      printf("All -i input files must have the same number of Label colors\n");
	    } else
	      printf("Line %d in %s has %s %d which is different from value %d found previously:\n%s\n",linecnt,filename,key,NColors,colors,buf);
	    fflush(stdout); exit(1);
	  }
	} else {
	  for(int c = colors; c < NColors; c++)
	    if(!Nickase[c]){
	      Nickase[c] = strdup((char*) "unknown");
	      if(VERB>=2){
		printf("input_bnx(%s):Line %d:Setting global Nickase[%d]=%p(%s)\n",filename,linecnt,c,Nickase[c],Nickase[c]);
		fflush(stdout);
	      }
	    }
	  colors = NColors;
	  colors_set = 1;
	  if(colors == 1 && usecolor==1)
	    usecolor = 0;
	  if(VERB){
	    printf("Line %d in %s: setting colors=%d, colors_set=%d, usecolor= %d due to Label Channels: %d:\n%s\n",linecnt,filename,colors,colors_set,usecolor,NColors,buf);
	    fflush(stdout);
	  }
	}
	if(VERB>=2){
	  printf("Line %d in %s: NColors=%d : colors=%d,colors_set=%d\n",linecnt,filename,NColors,colors,colors_set);
	  fflush(stdout);
	}
	int nlines = 1 + colors*(1 + MapTrueSite + MapSNR + MapIntensity + MapStitched + MapStitchLoc + MapPSFWidth + MapImageLoc);
	if(nlines != lines_per_map){
	  lines_per_map = nlines;
	  maxlines = numthreads * MAPS_PER_THREAD * lines_per_map;
	  Cline *nlines = new Cline[maxlines];
	  for(int k = 0; k < numlines; k++){
	    nlines[k].buf = lines[k].buf;
	    nlines[k].linecnt = lines[k].linecnt;
	  }
	  delete [] lines;
	  lines = nlines;
	}

	continue;
      }
      
      for(int c = 0;  c < colors; c++){
	char keybuf[256];
	sprintf(keybuf,"Nickase Recognition Site %d:",c+1);
	if((pt = strstr(&buf[1],keybuf))){
	  pt += strlen(keybuf);
	  while(*pt && isspace(*pt))
	    pt++;
	  if(!*pt){
	    printf("WARNING: %s on line %d of %s : Nickase string not found at end of line\n",keybuf,linecnt, filename);
	    fflush(stdout);
	    break;// ignore error
	  }
	  /* convert input to lower case */
	  for(int i = strlen(pt)-1; isspace(pt[i]);)
	    pt[i--] = '\0';
	  for(int i = strlen(pt); --i >= 0;)
	    pt[i] = tolower(pt[i]);

	  if(DEBUG) assert(Nickase[c]);
	  if(strcmp(Nickase[c],"unknown")){/* enzyme name was specified in a previous input file */
	    if(strcmp(Nickase[c],pt) && strcmp(pt,"unknown")){/* Enzyme strings are not the same AND neither string is "unknown" */
	      if(BNXMinorVersion >= 3){/* need to check for ; seperated color field */
		char *origColor = strchr(Nickase[c],';');
		char *newColor = strchr(pt,';');
		if(origColor || newColor){
		  if(origColor && newColor){/* If both files specify a color, they must match exactly */
		    printf("ERROR:%s on line %d of %s:\n\t Nickase[%d] string \"%s\" does not match value in previous file of \"%s\"\n",
			   keybuf,linecnt,filename, c, pt, Nickase[c]);
		    fflush(stdout);exit(1);
		  }
		  const char *DefaultColor = c ? "red" : "green";
		  if(origColor && strncmp(&origColor[1], DefaultColor, strlen(DefaultColor))){
		    printf("ERROR:%s on line %d of %s:\n\t Nickase[%d] string \"%s\" is not compactible with value in previous file of \"%s\"\n",
			   keybuf,linecnt,filename, c, pt, Nickase[c]);
		    fflush(stdout);exit(1);
		  }
		  if(newColor && strncmp(&newColor[1], DefaultColor, strlen(DefaultColor))){
		    printf("ERROR:%s on line %d of %s:\n\t Nickase[%d] string \"%s\" is not compactible with value in previous file of \"%s\"\n",
			   keybuf,linecnt,filename, c, pt, Nickase[c]);
		    fflush(stdout);exit(1);
		  }
		  if(newColor) {/* new file specifies color, previous files did NOT : update Nickase[c] to include color */
		    if(VERB/* HERE HERE >=2 */){
		      printf("input_bnx(%s):Line %d:Changing global Nickase[%d] = %s -> %s\n",filename,linecnt,c,Nickase[c],pt);
		      fflush(stdout);
		    }
		    // NOTE : cannot free Nickase[c] since other maps may already point to it. This results in a small memory leakage at end of program
		    Nickase[c] = strdup(pt);
		  }
		} 
	      } else {/* In BNX 1.2 Enzyme strings must match exactly, unless one of them is "unknown" */
		printf("ERROR:%s on line %d of %s:\n\t Nickase[%d] string \"%s\" does not match value in previous file of \"%s\"\n",
		       keybuf,linecnt,filename, c, pt, Nickase[c]);
		fflush(stdout);exit(1);
	      }
	    }
	  } else if(strcmp(pt,"unknown")) {
	    if(VERB/* HERE HERE >=2*/){
	      printf("input_bnx(%s):line %d:Changing global Nickase[%d]=%p(%s) to ",filename,linecnt,c,Nickase[c],Nickase[c]);
	      fflush(stdout);
	    }

	    // NOTE : cannot free Nickase[c] since other maps may already point to it. This results in a small memory leakage at end of program
	    //	    if(Nickase[c]) free(Nickase[c]);
	    Nickase[c] = strdup(pt);
	    if(VERB/* HERE HERE >=2 */){
	      printf("%p(%s)\n",Nickase[c],Nickase[c]);
	      fflush(stdout);
	    }
	  }
	}
      }

      key = (char *)"Run Data";
      if(BNXVersion >= 1 && (pt=strstr(&buf[1],key))){
	if(RunDataIdCnt >= RunDataIdMax){
	  RunDataIdMax = max(2*RunDataIdMax,16);
	  int *oldRunDataIdmap = RunDataIdmap;
	  RunDataIdmap = new int[RunDataIdMax];
	  for(int k = 0; k < RunDataIdCnt; k++)
	    RunDataIdmap[k] = oldRunDataIdmap[k];
	  delete [] oldRunDataIdmap;
	}

	/* check if the line is already in the list */
	while(isspace(buf[len-1]))
	  buf[--len] = '\0';/* strip off white space & CR/LF from end of line */
	int i = 0;
	for(; i < RunDataListLen; i++)
	  if(!strcmp(buf,RunDataList[i]))
	    break;

	if(i < RunDataListLen){/* already at RunDataList[i] */
	  RunDataIdmap[RunDataIdCnt++] = i;
	  if(bnxmergeerror){
	    printf("WARNING:Duplicate Run Data on line %d of %s (present earlier in same file or in a previous -i input file)\n",linecnt,filename);
	    printf("Use -bnxmergeerror 0 if you want this error to be ignored\n");
	    fflush(stdout);
	    fflush(stdout);exit(1);
	  }

	} else {
	  /* copy line to RunDataList[] */
	  if(RunDataListLen >= RunDataListMax){
	    RunDataListMax = max(2*RunDataListLen, 16);
	    char **oldRunDataList = RunDataList;
	    RunDataList = new char*[RunDataListMax];
	    for(int i = 0; i < RunDataListLen; i++)
	      RunDataList[i] = oldRunDataList[i];
	    delete [] oldRunDataList;
	  }
	  RunDataIdmap[RunDataIdCnt++] = RunDataListLen;
	  RunDataList[RunDataListLen++] = strdup(buf);
	}

	continue;
      }

      key = (char *)"Min Label SNR:";
      if((pt = strstr(&buf[1],key))){
	if(!colors_set){
	  printf("Label Channels not specified before line %d in %s:\n%s\n",linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	pt += strlen(key);

	int c = 0;
	for(; c < 1 /* HERE colors */; c++){
	  double MinSNR = strtod(pt,&qt);
	  if(qt == pt || !(*qt == 0 || isspace(*qt))){
	    printf("Line %d in %s with %s has invalid value for Label color %d (must be non-negative number):\n%s\n",linecnt,filename,key,c+1,buf);
	    fflush(stdout);exit(1);
	  }
	  if(MinSNR < 0.0){
	    printf("WARNING:Line %d in %s with %s has invalid value for Label color %d (must be non-negative number):\n%s\n",linecnt,filename,key,c+1,buf);
	    fflush(stdout);
	    MinSNR = 0.0;
	  }
	  pt = qt;
	  BNXminSNR[c] = max(BNXminSNR[c],MinSNR);
	}
	for(; c < MAXCOLOR/* colors */; c++)
	  BNXminSNR[c] = BNXminSNR[0];

	continue;
      }

      key = (char *)"Min Molecule Length (Kb):";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);
	
	double MinLen = strtod(pt,&qt);
	if(qt == pt || !(*qt == 0 || isspace(*qt))){
	  printf("Line %d in %s with %s has invalid value:\n%s\n",linecnt,filename,key,buf);
	  fflush(stdout);exit(1);
	}
	if(MinLen < 0.0){
	  printf("WARNING:Line %d in %s with %s has invalid value (must be non-negative number):\n%s\n",linecnt,filename,key,buf);
	  fflush(stdout);
	  MinLen = 0.0;
	}
	pt = qt;
	BNXminlen = max(MinLen,BNXminlen);
	continue;
      }

      key = (char *)"Original QueryMap:";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);

	while(isspace(buf[len-1]))
	  buf[--len] = '\0';/* strip off white space & CR/LF from end of line */

	if(QueryOrigin && strcmp(QueryOrigin,pt)){
	  printf("Original QueryMap=%s on line %d in %s does not match previous value of %s (in previous -i input files)\n",pt,linecnt,filename,QueryOrigin);
	  fflush(stdout);exit(1);
	}

	QueryOrigin = strdup(pt);
	if(DEBUG) assert(MapTrueSite == 1);
	continue;
      }

      /* Bases per Pixel in bnx file is ignored, using default 500 bpp */

      continue;/* skip all other comment lines */      
    }

    /* done with header comment line */
    if(NColors <= 0){
      printf("lines with \"Label Channels:\" missing or invalid value on line %d of %s:\n%s\n",
	     linecnt,filename,buf);
      fflush(stdout);exit(1);
    }
    if(NColors > MAXCOLOR){
      printf("Lable Channels value %d is too large (MAXCOLOR=%d)\n",NColors,MAXCOLOR);
      fflush(stdout);exit(1);
    }
    if(colors != NColors){
      printf("N Colors=%d in file %s does not match parameter value=%d\n",NColors,filename,colors);
      exit(0);
    }

    if(DEBUG) assert(numlines==0);

    lines[numlines].buf = buf;
    lines[numlines].linecnt = linecnt;
    lines[numlines].RunIndex = 0;
    if(BNXVersion >= 1){
      if(!RunDataIdCnt){
	if(!RunDataWarning){
	  printf("WARNING:Missing Run Data line before line %d in %s (required for BNXVersion=%d)\n",linecnt, vfixx_filename[fileid], BNXVersion);
	  printf("    Suggest regenerating %s with newest RefAligner to add those lines : otherwise ScanNumbers may be wrong and map ids may have to be renumbered\n", filename);
	  fflush(stdout);
	  RunDataWarning = 1;
	}
	lines[numlines].RunIndex = max(1,RunDataListLen)-1;
      } else {
	if(DEBUG) assert(RunDataIdCnt <= RunDataIdMax);
	lines[numlines].RunIndex = RunDataIdmap[RunDataIdCnt - 1];
	if(DEBUG) assert(0 <= lines[numlines].RunIndex && lines[numlines].RunIndex < RunDataListLen);
      }
    }
    numlines++;
    linecnt++;

    break;
  }
  
  idfiltered = 0;
  
  if(BNXVersion == -1){
    printf("No Version specified in %s\n",filename);
    fflush(stdout);exit(1);
  }
  if(BNXVersion > 1){
    printf("Unknown Version id=%d while reading %s\n",BNXVersion,filename);
    fflush(stdout);exit(1);
  }

  tMaxSites = (minSNRestimate ? MaxSites2 : MaxSites);
  tMinLen = (ScanCorrection ? min((ScanCorrectMinLen < 0.0 ? 9999.0 : ScanCorrectMinLen), MinLen * min(0.85, 1.0 - ScaleDelta * ScaleDeltaSize)) : MinLen);
  tMaxLen = (ScanCorrection ? MaxLen * max(1.15, 1.0 + ScaleDelta * ScaleDeltaSize) /* WAS14 1.15 */ : MaxLen);

  origRunDataIdCnt = RunDataIdCnt;

  size_t Q_TRUESITE_LEN[2] = {strlen(VQ_TRUESITE[0][BNXVersion]),strlen(VQ_TRUESITE[1][BNXVersion])};
  size_t Q_SNR_LEN[2] = {strlen(VQ_SNR[0][BNXVersion]),strlen(VQ_SNR[1][BNXVersion])};
  size_t Q_INTENSITY_LEN[2] = {strlen(VQ_INTENSITY[0][BNXVersion]),strlen(VQ_INTENSITY[1][BNXVersion])};
  size_t Q_STITCH_LEN[2] = {strlen(VQ_STITCH[0][BNXVersion]),strlen(VQ_STITCH[1][BNXVersion])};
  size_t Q_STITCHLOC_LEN[2] = {strlen(VQ_STITCHLOC[0][BNXVersion]),strlen(VQ_STITCHLOC[1][BNXVersion])};
  size_t Q_PSFWIDTH_LEN[2] = {strlen(VQ_PSFWIDTH[0][BNXVersion]),strlen(VQ_PSFWIDTH[1][BNXVersion])};
  size_t Q_IMAGELOC_LEN[2] = {strlen(VQ_IMAGELOC[0][BNXVersion]),strlen(VQ_IMAGELOC[1][BNXVersion])};

  char *key = (char *)"Run Data";

  long pos_cur, pos_end, len, plen = 0;
  
  pos_cur = ftell(fp);
  fseek(fp, 0L, SEEK_END);
  pos_end = ftell(fp);
  fseek(fp, pos_cur, SEEK_SET);
  
#if FILE_BUFFER_SIZE==0 // OLD code
  char *buffer = new char[pos_end-pos_cur+2];
  len = fread(buffer, 1, pos_end-pos_cur+1, fp);
  buffer[len] = 0;
#else // NEW code
  long fbsiz = FILE_BUFFER_SIZE;
  long bsiz = min(pos_end - pos_cur + 2, fbsiz);

  int skipped_lines = 0;
  fastfilterCnt = 0;// NEW102

  char *buffer = new char[bsiz];
  plen = 0;// previous content of buffer is at buffer[0..plen-1]

  while(pos_cur < pos_end){
    if(DEBUG>=2){/* verify that current location in file is pos_cur */
      errno = 0;
      long pos = ftell(fp);
      if(pos < 0){
	int eno = errno;
	char *err = strerror(eno);
	printf("input_bnx:ftell(fp) failed on bnx file %s:errno=%d:%s\n",filename, eno,err);
	fflush(stdout);exit(1);
      }
      if(DEBUG && !(pos == pos_cur)){
	printf("ftell(fp) for bnx file %s returned %ld : expected pos_cur= %ld\n",filename,pos,pos_cur);
	fflush(stdout);
	assert(pos == pos_cur);
      }
    }

    long tlen = min(pos_end - pos_cur + 1, fbsiz - 1 - plen);

    len = fread(&buffer[plen], 1, tlen, fp);

    buffer[plen+len] = 0;// null terminate buffer

    if(VERB>=2){
      printf("input_bnx : loaded %ld/%ld bytes from bnx file %s (from pos= %ld,tlen= %ld) into buffer[%ld..%ld],linecnt=%d: time=%0.6f(%0.6f)\n",
	     len, pos_end, vfixx_filename[fileid], pos_cur, tlen, plen, plen+len-1, linecnt, mtime(),wtime());
      fflush(stdout);
    }
    
    long blen = strlen(buffer);
    if(DEBUG) assert(blen <= plen+len);
    if(DEBUG) assert(blen >= plen);
    if(blen < plen + len){
      printf("input_bnx : found NULL character in BNX file %s (at char %ld) : Invalid BNX file\n", filename, pos_cur + blen - plen);
      fflush(stdout);exit(1);
    }

#endif
    pos_cur += len;

    char *p = buffer;

    for(char *q; *p; p = q, linecnt++) {
      if(TSAN){/* avoid calling strstr(), since TSAN calls strlen() with every call! */
	for(int i = 0; ; i++){
	  q = &p[i];
	  if(*q == '\n' || !*q)
	    break;
	}
	if(*q == '\0'){
	  if(pos_cur >= pos_end){
	    printf("Unterminated line %d in %s (length=%lu bytes)\n",linecnt,filename,strlen(p));
	    fflush(stdout);exit(1);
	  }
	  break;	    
	}
      } else {
	if((q = strchr(p,'\n')) == NULL){
	  if(pos_cur >= pos_end){
	    printf("Unterminated line %d in %s (length=%lu bytes)\n",linecnt,filename,strlen(p));
	    fflush(stdout);exit(1);
	  }
	  break;
	}
      }

      size_t len = q - p;
      *q++ = 0;

      if(*p == '#'){
	/* check if we have a "Run Data" comment in the middle of the file, possibly due to multiple BNX files being concatenated.
	   NOTE : This is a harmless extension of the standard to allow multiple raw BNX 1.2 or 1.0 files to be concatenated. */
	if(BNXVersion >=1 && (pt=strstr(&buf[1],key))){
	  if(RunDataIdCnt >= RunDataIdMax){
	    RunDataIdMax = max(2*RunDataIdMax,16);
	    int *oldRunDataIdmap = RunDataIdmap;
	    RunDataIdmap = new int[RunDataIdMax];
	    for(int k = 0; k < RunDataIdCnt; k++)
	      RunDataIdmap[k] = oldRunDataIdmap[k];
	    delete [] oldRunDataIdmap;
	  }

	  char *buf2 = strdup(buf);

	  /* check if the line is already in the list */
	  while(isspace(buf2[len-1]))
	    buf2[--len] = '\0';/* strip off white space & CR/LF from end of line */
	  int i = 0;
	  for(; i < RunDataListLen; i++)
	    if(!strcmp(buf2,RunDataList[i]))
	      break;
	  if(i < RunDataListLen){/* already at RunDataList[i] */
	    RunDataIdmap[RunDataIdCnt++] = i;
	  } else {
	    /* copy line to RunDataList[] */
	    if(RunDataListLen >= RunDataListMax){
	      RunDataListMax = max(2*RunDataListLen, 16);
	      char **oldRunDataList = RunDataList;
	      RunDataList = new char*[RunDataListMax];
	      for(int i = 0; i < RunDataListLen; i++)
		RunDataList[i] = oldRunDataList[i];
	      delete [] oldRunDataList;
	    }
	    RunDataIdmap[RunDataIdCnt++] = RunDataListLen;
	    RunDataList[RunDataListLen++] = strdup(buf2);
	  }
	  free(buf2);

	  skipped_lines++;
	  if(FILE_BUFFER_SIZE > 0)
	    q[-1] = '\n';/* restore EOL from skipped lines in case they end up in the unused buffer end */

	  continue;
	}

	skipped_lines++;
	if(FILE_BUFFER_SIZE > 0)
	  q[-1] = '\n';/* restore EOL from skipped lines in case they end up in the unused buffer end */

	continue;
      }

      if(*p == 'Q') {
	int nlines = 1 + colors*(1 + MapTrueSite + MapSNR + MapIntensity + MapStitched + MapStitchLoc + MapPSFWidth + MapImageLoc);
	if(numlines <= nlines){/* first map : use to detect which Q code are present for each map */
	  for(int c = 0; c < colors; c++){
	    if(!MapTrueSite && !strncmp(p,VQ_TRUESITE[c][BNXVersion],Q_TRUESITE_LEN[c])){
	      if(VERB){
		printf("Detected True Site information (%s) on line %d in %s\n",VQ_TRUESITE[c][BNXVersion],linecnt,vfixx_filename[fileid]);
		fflush(stdout);
	      }
	      if(fileid > 0 && MapTrueSite_start == 0){
		QXmismatch = fileid;
		if(QXerror){
		  printf("ERROR: True Site information (%s) missing in %s (must be present in all BNX input files)\n",VQ_TRUESITE[c][BNXVersion],vfixx_filename[0]);
		  fflush(stdout);exit(1);
		}
	      }
	      MapTrueSite = 1;
	    }
	    if(!MapSNR && !strncmp(p,VQ_SNR[c][BNXVersion],Q_SNR_LEN[c])){
	      if(VERB){
		printf("Detected SNR information (%s) on line %d in %s\n",VQ_SNR[c][BNXVersion],linecnt,vfixx_filename[fileid]);
		fflush(stdout);
	      }
	      if(fileid > 0 && MapSNR_start == 0){
		QXmismatch = fileid;
		if(QXerror){
		  printf("ERROR: SNR information (%s) missing in %s (must be present in all BNX input files)\n",VQ_SNR[c][BNXVersion],vfixx_filename[0]);
		  fflush(stdout);exit(1);
		}
	      }
	      MapSNR = 1;
	    }
	    if(!MapIntensity && !strncmp(p,VQ_INTENSITY[c][BNXVersion],Q_INTENSITY_LEN[c])){
	      if(VERB){
		printf("Detected Intensity information (%s) on line %d in %s\n",VQ_INTENSITY[c][BNXVersion],linecnt,vfixx_filename[fileid]);
		fflush(stdout);
	      }
	      if(fileid > 0 && MapIntensity_start == 0){
		QXmismatch = fileid;
		if(QXerror){
		  printf("ERROR: Intensity information (%s) missing in %s (must be present in all BNX input files)\n",VQ_INTENSITY[c][BNXVersion],vfixx_filename[0]);
		  fflush(stdout);exit(1);
		}
	      }
	      MapIntensity = 1;
	    }
	    if(!MapStitched && !strncmp(p,VQ_STITCH[c][BNXVersion],Q_STITCH_LEN[c])){
	      if(VERB){
		printf("Detected Stitch information (%s) on line %d in %s\n", VQ_STITCH[c][BNXVersion], linecnt, vfixx_filename[fileid]);
		fflush(stdout);
	      }
	      if(fileid > 0 && MapStitched_start == 0){
		QXmismatch = fileid;
		if(QXerror){
		  printf("ERROR: Stitch information (%s) missing in %s (must be present in all BNX input files)\n",VQ_STITCH[c][BNXVersion],vfixx_filename[0]);
		  fflush(stdout);exit(1);
		}
	      }
	      MapStitched = 1;
	    }
	    if(!MapStitchLoc && !strncmp(p,VQ_STITCHLOC[c][BNXVersion],Q_STITCHLOC_LEN[c])){
	      if(VERB){
		printf("Detected Stitch Location information (%s) on line %d in %s\n", VQ_STITCHLOC[c][BNXVersion], linecnt, vfixx_filename[fileid]);
		fflush(stdout);
	      }
	      if(fileid > 0 && MapStitchLoc_start == 0){
		QXmismatch = fileid;
		if(QXerror){
		  printf("ERROR: Stitch Location information (%s) missing in %s (must be present in all BNX input files)\n",VQ_STITCHLOC[c][BNXVersion],vfixx_filename[0]);
		  fflush(stdout);exit(1);
		}
	      }
	      MapStitchLoc = 1;
	    }
	    if(!MapPSFWidth && !strncmp(p,VQ_PSFWIDTH[c][BNXVersion],Q_PSFWIDTH_LEN[c])){
	      if(VERB){
		printf("Detected PSF Width information (%s) on line %d in %s\n",VQ_PSFWIDTH[c][BNXVersion],linecnt,vfixx_filename[fileid]);
		if(VERB>=2)
		  printf("fileid= %d, MapPSFWidth_start= %d\n",fileid, MapPSFWidth_start);
		fflush(stdout);
	      }
	      if(fileid > 0 && MapPSFWidth_start == 0){
		QXmismatch = fileid;
		if(QXerror){
		  printf("ERROR: PSF Width information (%s) missing in %s (must be present in all BNX input files)\n",VQ_PSFWIDTH[c][BNXVersion],vfixx_filename[0]);
		  fflush(stdout);exit(1);
		}
	      }
	      MapPSFWidth = 1;
	    }
	    if(!MapImageLoc && !strncmp(p,VQ_IMAGELOC[c][BNXVersion],Q_IMAGELOC_LEN[c])){
	      if(VERB){
		printf("Detected Image Location information (%s) on line %d in %s\n",VQ_IMAGELOC[c][BNXVersion],linecnt,vfixx_filename[fileid]);
		if(VERB>=2)
		  printf("fileid= %d, QXmismatch=%d, MapImageLoc_start= %d\n",fileid, QXmismatch, MapImageLoc_start);
		fflush(stdout);
	      }
	      if(fileid > 0 && MapImageLoc_start == 0){
		QXmismatch = fileid;
		if(QXerror){
		  printf("ERROR: Image Location information (%s) missing in %s (must be present in all BNX input files)\n",VQ_IMAGELOC[c][BNXVersion],vfixx_filename[0]);
		  fflush(stdout);exit(1);
		}
	      }
	      MapImageLoc = 1;
	    }
	  }
	} // if(numlines <= nlines)

	nlines = 1 + colors*(1 + MapTrueSite + MapSNR + MapIntensity + MapStitched + MapStitchLoc + MapPSFWidth + MapImageLoc);

	if(nlines > lines_per_map){
	  lines_per_map = nlines;
	  maxlines = numthreads * MAPS_PER_THREAD * lines_per_map;
	  Cline *nlines = new Cline[maxlines];
	  for(int k = 0; k < numlines; k++){
	    nlines[k].buf = lines[k].buf;
	    nlines[k].linecnt = lines[k].linecnt;
	    nlines[k].RunIndex = lines[k].RunIndex;
	  }
	  delete [] lines;
	  lines = nlines;
	}
	if(!(MapTrueSite && !(strncmp(p,VQ_TRUESITE[0][BNXVersion],Q_TRUESITE_LEN[0]) || (colors>1 && !strncmp(p,VQ_TRUESITE[1][BNXVersion],Q_TRUESITE_LEN[1]))))
	   && !(MapSNR && (!strncmp(p,VQ_SNR[0][BNXVersion],Q_SNR_LEN[0]) || (colors>1 && !strncmp(p,VQ_SNR[1][BNXVersion],Q_SNR_LEN[1])))) 
	   && !(MapIntensity && (!strncmp(p,VQ_INTENSITY[0][BNXVersion],Q_INTENSITY_LEN[0]) || (colors>1 && !strncmp(p,VQ_INTENSITY[1][BNXVersion],Q_INTENSITY_LEN[1])))) 
	   && !(MapStitched && (!strncmp(p,VQ_STITCH[0][BNXVersion],Q_STITCH_LEN[0]) || (colors>1 && !strncmp(p,VQ_STITCH[1][BNXVersion],Q_STITCH_LEN[1]))))
	   && !(MapStitchLoc && (!strncmp(p,VQ_STITCHLOC[0][BNXVersion],Q_STITCHLOC_LEN[0]) || (colors>1 && !strncmp(p,VQ_STITCHLOC[1][BNXVersion],Q_STITCHLOC_LEN[1]))))
	   && !(MapPSFWidth && (!strncmp(p,VQ_PSFWIDTH[0][BNXVersion],Q_PSFWIDTH_LEN[0]) || (colors>1 && !strncmp(p,VQ_PSFWIDTH[1][BNXVersion],Q_PSFWIDTH_LEN[1]))))
	   && !(MapImageLoc && (!strncmp(p,VQ_IMAGELOC[0][BNXVersion],Q_IMAGELOC_LEN[0]) || (colors>1 && !strncmp(p,VQ_IMAGELOC[1][BNXVersion],Q_IMAGELOC_LEN[1]))))
	   ){
	  if(Qcode_warning < 7){
	    #pragma omp critical
	    {
	      Qcode_warning++;
	      printf("WARNING:Unrecognized Qcode in %s on line %d (OR not present for first molecule in first input file),BNXVersion=%d: %s\n",vfixx_filename[fileid], linecnt, BNXVersion,p);
	      fflush(stdout);
	    }
	  }

	  skipped_lines++;
	  if(FILE_BUFFER_SIZE > 0)
	    q[-1] = '\n';/* restore EOL from skipped lines in case they end up in the unused buffer end */

	  continue;/* skip line */
	}
      } // if (*p == 'Q')

      lines[numlines].buf = p;
      lines[numlines].linecnt = linecnt;
      lines[numlines].RunIndex = 0;
      if(BNXVersion >= 1){
	if(!RunDataIdCnt){
	  if(!RunDataWarning){
	    printf("WARNING:Missing Run Data line before line %d in %s (required for BNXVersion=%d)\n",linecnt, vfixx_filename[fileid], BNXVersion);
	    printf("    Suggest regenerating %s with newest RefAligner to add those lines : otherwise ScanNumbers may be wrong and map ids may have to be renumbered\n", filename);
	    fflush(stdout);
	    RunDataWarning = 1;
	  }
	  lines[numlines].RunIndex = max(1,RunDataListLen)-1;
	} else {
	  if(DEBUG) assert(RunDataIdCnt <= RunDataIdMax);
	  lines[numlines].RunIndex = RunDataIdmap[RunDataIdCnt - 1];
	  if(DEBUG) assert(0 <= lines[numlines].RunIndex && lines[numlines].RunIndex < RunDataListLen);
	}
      }
      numlines++;

      if(QXmismatch && fileid == QXmismatch && MINSITE){/* need to deallocate any unused maps, since their block allocation from previous fileid cannot be used for current fileid */
	maxmaptruncate(nummaps, maxmaps, Map);

	QXmismatch = -1;
      }

      if(numlines >= maxlines){/* parse all lines[0..linecnt-1] in parallel section */
	if(DEBUG) assert((numlines % lines_per_map) == 0);
	parse_lines(numthreads,lines,numlines,lines_per_map,fileid,filename,nummaps,maxmaps,Map,ids, bnxFXY, bnxChim);
	total_lines += numlines;
	numlines = 0;
      }
    } // for(char *q;*p; p = q,linecnt++)
#if FILE_BUFFER_SIZE==0 // OLD CODE
    if(DEBUG) assert(*p == 0);
#else // FILE_BUFFER_SIZE > 0 (NEW CODE)

    if(VERB>=2){
      long rlen = &buffer[plen+len] - p;
      printf("Remaining lines = %d, Remaining bytes = %ld starting at buffer[%ld]\n", numlines,rlen,p-buffer);
      fflush(stdout);
    }

    if(numlines > 0){/* finish processing remaining lines in lines[0..linecnt-1], before buffer is overwritten */
      int nmaps = numlines / lines_per_map;
      parse_lines(numthreads,lines, nmaps * lines_per_map,lines_per_map,fileid,filename,nummaps,maxmaps,Map,ids,bnxFXY, bnxChim);
      total_lines += nmaps * lines_per_map;

      if(numlines > nmaps * lines_per_map){
	char *origp = p;
	long tlen = &buffer[plen+len] - p;// WAS156 size_t tlen = strlen(p);

	if(VERB>=2){
	  printf("buffer[] contains %lu chars from %d lines + %lu chars from a partial line:\n",
		 p - lines[nmaps * lines_per_map].buf, numlines - nmaps*lines_per_map, tlen);
	  if(VERB>=2){
	    for(int t = nmaps*lines_per_map; t < numlines; t++)
	      printf("line %d:%s\n", lines[t].linecnt, lines[t].buf);
	    printf("partial line %d:%s\n", linecnt, p);
	  }
	  fflush(stdout);
	}

	/* add '\n' back to each line */
	for(int t = nmaps * lines_per_map; t < numlines; t++){
	  size_t len = strlen(lines[t].buf);
	  if(DEBUG) assert(lines[t].buf[len] == 0);
	  if(VERB>=2){
	    printf("t=%d:tlen= %lu -> %lu, strlen(lines[t].buf)= %lu, lines[t].buf:\n%s\n",t,tlen,tlen+len+1,len,lines[t].buf);
	    fflush(stdout);
	  }
	  lines[t].buf[len] = '\n';
	  tlen += len + 1;
	}
	p = lines[nmaps * lines_per_map].buf;
	long rlen = &buffer[plen+len] - p;
	if(VERB>=2){
	  printf("Restored surplus lines : buffer[] now contains a total of %ld chars (should be tlen= %ld)\n",rlen,tlen);
	  fflush(stdout);
	}
	if(DEBUG && !skipped_lines && !(rlen == tlen)){
	  printf("numlines=%d, nmaps=%d, lines_per_map=%d, skipped=%d,tlen= %lu, strlen(p)= %lu, strlen(origp)= %lu:origp=\n%s\np=\n%s\n",numlines,nmaps,lines_per_map,skipped_lines,tlen,strlen(p),strlen(origp),origp,p);
	  fflush(stdout);
	  assert(rlen == tlen);
	}
      }
      numlines -= nmaps * lines_per_map;

      if(numlines > 0){  /* remaining lines from partial map will be re-parsed from buffer[] starting at p */
	linecnt -= numlines;
	numlines = 0;
      }
    }
    long rlen = &buffer[plen+len] - p;// WAS157 : strlen(p); // NOTE : cannot assume file does NOT contain NULLs, so need to allow parser to read nulls
    if(VERB>=2){
      printf("Moving remaining %ld bytes from buffer[%ld .. %ld] to start of buffer:\n",rlen, p-buffer, plen+len-1);
      fflush(stdout);
    }
    memmove(buffer, p, rlen);// NOTE : memmove() allows source and destination regions to overlap
    plen = rlen;
  } // while(pos_cur < pos_end)
  if(DEBUG) assert(numlines == 0);
#endif

  if(numlines > 0){
    if(numlines % lines_per_map){
      printf("WARNING:total lines=%d (excluding comments) in %s is not a multiple of %d (%d colors,TrueSite=%d,SNR=%d,stitched=%d,stitchLoc=%d,PSFWidth=%d,ImageLoc=%d): Ignoring last %d lines\n",
	     numlines, filename,lines_per_map,colors,MapTrueSite,MapSNR,MapStitched, MapStitchLoc, MapPSFWidth, MapImageLoc, numlines % lines_per_map);
      fflush(stdout);
      numlines -= numlines % lines_per_map;
    }
    parse_lines(numthreads,lines,numlines,lines_per_map,fileid,filename,nummaps,maxmaps,Map,ids,bnxFXY, bnxChim);
    total_lines += numlines;
  }

  if(VERB>=2){
    printf("Finished reading in %s: %d new maps : wall time= %0.6f secs\n",filename,nummaps-orignummaps,wtime());
    fflush(stdout);
  }

  (void)fclose(fp);
  
  if(fbuf){
    free(fbuf);
    fbuf = 0;
  }
  delete [] lines;
  delete [] buffer;

  if(BNX_FASTFILTER){
    int i,j;
    for(i = j = orignummaps; i < nummaps; i++){
      pmap = Map[i];
      if(pmap->numsite[0] < 0)
	continue;

      if(j < i){
	Cmap *tmp = Map[j];
	pmap->mapid = j;
	Map[j] = pmap;
	if(tmp)
	  tmp->mapid = i;
	Map[i] = tmp;
      }
      j++;
    }

    if(VERB && (j != nummaps || fastfilterCnt > 0)){
      printf("Reduced number of input maps from %d to %d due to MinLen=%0.3f,MaxLen=%0.3f (reflecting bpp= %0.3f), MinSites=%d in file %s (%d of %d)\n",
	     total_lines / lines_per_map, j-orignummaps,tMinLen,tMaxLen,500.0*PixelLen/origPixelLen,MinSites,vfixx_filename[fileid], fileid+1, num_files);
      fflush(stdout);
    }

    nummaps = j;
  }

  if(VERB){
    long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
    getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
    printf("VmRSS= %0.4f, VmHWM= %0.4f : CPU time= %0.6f, wall time= %0.6f\n",VmRSS * 1e-9, VmHWM * 1e-9, mtime(), wtime());      
    fflush(stdout);
  }

  Cmap **maplist = new Cmap*[nummaps-orignummaps];

  if(VERB/* HERE >=2 */ && nummaps > orignummaps){
    double totlen = 0.0, totlen2 = 0.0;
    size_t nsites = 0, nsites2 = 0;

    for(int i= orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      maplist[i-orignummaps] = pmap;
      int M = pmap->numsite[0];
      totlen += pmap->site[0][M+1];
      nsites += M;
      if(colors >= 2){
	int M2 = pmap->numsite[1];
	totlen2 += pmap->site[1][M2 + 1];
	nsites2 += M2;
      }
    }

    qsort(maplist, nummaps-orignummaps, sizeof(Cmap*), (intcmp*)PmapLenDec);/* sort map pointers in descending order of length */

    double N50Sum= 0.0, N50Len = 0.0;
    for(int i = orignummaps; i < nummaps; i++){
      Cmap *pmap = maplist[i-orignummaps];
      int M = pmap->numsite[0];
      if((N50Sum += pmap->site[0][M+1]) >= totlen * 0.5){
	N50Len = pmap->site[0][M+1];
	break;
      }
    }
    if(DEBUG && totlen > 0.0) assert(N50Len > 0.0);

    if(colors >= 2)
      printf("Finished parsing maps in %s: maps=%d, sites=%lu,%lu length= %0.3fkb (avg= %0.3f kb, label density1= %0.3f /100kb, label density2= %0.3f /100kb, N50= %0.4f Mb): wall time= %0.6f secs\n",
	     vfixx_filename[fileid], nummaps-orignummaps,nsites,nsites2,totlen, 0.5 * (totlen+totlen2)/(nummaps-orignummaps), nsites*100.0/max(0.001,totlen), nsites2*100.0/max(0.001,totlen2),
	     N50Len * 0.001, wtime() );
    else
      printf("Finished parsing maps in %s: maps=%d, sites=%lu, length= %0.3fkb (avg= %0.3f kb, label density= %0.3f /100kb, N50= %0.4f Mb): wall time= %0.6f secs\n",
	     vfixx_filename[fileid], nummaps-orignummaps,nsites,totlen,totlen/(nummaps-orignummaps), nsites*100.0/max(0.001,totlen),N50Len * 0.001,wtime());
    fflush(stdout);

    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  if(Tracksites && !MapTrueSite){ /* append current file name to QueryOrigin */
    if(!QueryOrigin)
      QueryOrigin = strdup(filename);
    else {
      size_t len = strlen(QueryOrigin) + strlen(filename) + 2;
      char *nQueryOrigin = (char *)malloc(len * sizeof(char));
      if(!nQueryOrigin){
	printf("malloc(%lu) failed\n",len);
	fflush(stdout);exit(1);
      }
      sprintf(nQueryOrigin,"%s,%s",QueryOrigin,filename);
      free(QueryOrigin);
      QueryOrigin = nQueryOrigin;
    }
  }

  if(BIAS_TRACE){/* print map id BIAS_TRACE_ID */
    for(int m = orignummaps; m < nummaps; m++){
      Cmap *pmap = Map[m];
      for(int c = 0; c < colors; c++){
	if(pmap->id == BIAS_TRACE_ID && pmap->SNR[c]){
	  int N = pmap->numsite[c];
	  double *SNR = pmap->SNR[c];
	  FLOAT *Y =  pmap->site[c];
	  printf("After parsing maps: id=%lld,N=%d:c=%d,pmap->len= %0.4f\n",pmap->id,N,c,pmap->len);
	  for(int J = 1; J <= N+1; J++)
	    printf(" Y[%d]=%0.4f, SNR=%0.4f\n",J,Y[J], (J <= N) ? SNR[J] : -1.0);
	  fflush(stdout);
	}
      }
    }
  }
  
  if(MaxEnd >= 0.0 || (fabs(PixelLen - origPixelLen) > 1e-12 && !NoBpp) || mres >= 0.0 ||
     (MapSNR && (minSNR[0] > 0.0 || (colors >= 2 && minSNR[1] > 0.0)))){/* mresSD = 0 : merge map sites at distance mres or closer
		     mresSD > 0 : merge map sites at distance X with probability 0.5(1-erf((X-mres)/(mresSD*sqrt(2)))) */
    if(VERB>=2){
      printf("Calling mapcorrect with nummaps=%d,orignummaps=%d,fileid=%d\n",nummaps,orignummaps,fileid);
      fflush(stdout);
    }

    mapcorrect(Map, orignummaps, nummaps, maxmaps, fileid, 1);
  } // mres >= 0.0

  if(BIAS_TRACE){/* print map id BIAS_TRACE_ID */
    for(int m = orignummaps; m < nummaps; m++){
      Cmap *pmap = Map[m];
      if(pmap->id == BIAS_TRACE_ID && pmap->SNR[0]){
	int N = pmap->numsite[0];
	double *SNR = pmap->SNR[0];
	FLOAT *Y =  pmap->site[0];
	printf("After applying MaxEnd=%0.4f,mres= %0.6f,mresSD=%0.6f: id=%lld,N=%d:pmap->len=%0.4f\n",MaxEnd, mres,mresSD,pmap->id,N,pmap->len);
	for(int J = 1; J <= N+1; J++)
	  printf(" Y[%d]=%0.4f, SNR=%0.4f\n",J,Y[J], (J <= N) ? SNR[J] : -1.0);
	fflush(stdout);
      }
    }
  }

  if(VERB && ids && idfiltered > 0){
    printf("Filtered %d input maps that did not match id set, %d maps remaining (fileid=%d,filename=%s)\n", idfiltered,nummaps-orignummaps,fileid,vfixx_filename[fileid]);
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  if(DEBUG && colors >= 2){
    for(int i= orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];// WAS12 map[i]
      FLOAT len1 = pmap->site[0][pmap->numsite[0]+1];
      FLOAT len2 = pmap->site[1][pmap->numsite[1]+1];
      if(!(fabs(len1 - len2) < 1e-3)){
	printf("Map[%d] : id=%lld, M1=%d,len1=%0.4f, M2=%d,len2=%0.4f\n",i, pmap->id, pmap->numsite[0], len1, pmap->numsite[1], len2);
	fflush(stdout);
	assert(fabs(len1 - len2) < 1e-3);
      }
    }
  }

  if(VERB/* HERE >=2 && nummaps > orignummaps */){
    double totlen = 0.0, totlen2 = 0.0;
    size_t nsites = 0, nsites2 = 0;
    for(int i= orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      maplist[i-orignummaps] = pmap;      
      int M = pmap->numsite[0];
      totlen += pmap->site[0][M+1];
      nsites += M;
      if(colors >= 2){
	int M2 = pmap->numsite[1];
	totlen2 += pmap->site[1][M2 + 1];
	nsites2 += M2;
      }
    }
    qsort(maplist, nummaps-orignummaps, sizeof(Cmap*), (intcmp*)PmapLenDec);/* sort map pointers in descending order of length */
    double N50Sum= 0.0, N50Len = 0.0;
    for(int i = orignummaps; i < nummaps; i++){
      Cmap *pmap = maplist[i-orignummaps];
      int M = pmap->numsite[0];
      if((N50Sum += pmap->site[0][M+1]) >= totlen * 0.5){
	N50Len = pmap->site[0][M+1];
	break;
      }
    }
    if(DEBUG && totlen > 0.0) assert(N50Len > 0.0);
    if(colors >= 2)
      printf("After applying -maxEnd -bpp -mres -minSNR to %s : maps=%d, sites=%lu,%lu, length= %0.3f kb (avg= %0.3f kb, label density1= %0.3f /100kb, label density2= %0.3f /100kb, N50= %0.4f Mb):wtime=%0.6f\n",
	     vfixx_filename[fileid], nummaps-orignummaps,nsites,nsites2,totlen, 0.5 * (totlen+totlen2)/(nummaps-orignummaps), nsites*100.0/max(0.001,totlen), nsites2*100.0/max(0.001,totlen2),N50Len*0.001,wtime());
    else
      printf("After applying -maxEnd %0.4f -bpp %0.4f -mres %0.8f -minSNR %0.4f to %s : maps=%d, sites=%lu, length= %0.3f kb (avg= %0.3f kb, label density= %0.3f /100kb, N50= %0.4f Mb):wtime=%0.6f\n",
	     MaxEnd, PixelLen * 1e3, mres, minSNR[0]/* NEW4 */, vfixx_filename[fileid], nummaps-orignummaps,nsites,totlen,totlen/(nummaps-orignummaps), nsites*100.0/max(0.001,totlen), N50Len*0.001,wtime());
    fflush(stdout);
  }

  /* filter out maps smaller than MinLen or larger than MaxLen or with fewer sites than MinSites or greater than MaxSites or AvgIntensity > MaxIntensity or any PSFWidth value > maxPSFWidth */
  int i,j = nummaps;

  if(!ids){
    if(VERB>=2){
      printf("Filter maps with MinLen=%0.4f,MaxLen=%0.4f,MinSites=%d,MaxSites=%d,MaxIntensity=%0.6f:orignummaps=%d,nummaps=%d\n",
	     tMinLen,tMaxLen,MinSites,tMaxSites,MaxIntensity,orignummaps,nummaps);
      fflush(stdout);
    }
    int c = 0;
    if(colors == 2 && usecolor==2)
      c = 1;

    for(i=j=orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      if(DEBUG && colors>=2) assert(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3);

      if(VERB>=2 && pmap != NULL){
	printf("Checking Map[%d]:id=%lld,numsites[c=%d]=%d(MinSites=%d,MaxSites=%d),len=%0.3f(MinLen=%0.3f,MaxLen=%0.3f),startloc=%0.3f,endloc=%0.3f,TrueFlip=%d,len=%0.3f\n",
	       i,pmap->id,c,pmap->numsite[c],MinSites,tMaxSites, pmap->site[0][pmap->numsite[0]+1],MinLen,MaxLen,pmap->startloc,pmap->endloc,pmap->flipped,pmap->len);
	fflush(stdout);
      }

      int nsites = pmap->numsite[c];
      if(colors==2 && !usecolor)
	nsites += pmap->numsite[1-c];
      double len = pmap->site[0][pmap->numsite[0]+1];

      if(pmap==NULL || nsites < MinSites || (tMaxSites && nsites > tMaxSites) || 
	 len < tMinLen || (tMaxLen > 0.0 && len > tMaxLen) || 
	 (MaxSiteDensity > 0.0 && !minSNRestimate && nsites * 100.0 > MaxSiteDensity * len * colors) || 
	 (MinSiteDensity > 0.0 && !minSNRestimate && nsites * 100.0 < MinSiteDensity * len * colors)){

	if(VERB>=2 && pmap != NULL){
	  if(colors==2 && !usecolor)
	    printf("Deleting Map[%d]:id=%lld,numsites[0,1]=%d+%d(MinSites=%d,MaxSites=%d),len=%0.3f(MinLen=%0.3f,MaxLen=%0.3f),startloc=%0.3f,endloc=%0.3f,TrueFlip=%d,len=%0.3f\n",
		   i,pmap->id,pmap->numsite[0],pmap->numsite[1],MinSites,tMaxSites, len,tMinLen,tMaxLen,pmap->startloc,pmap->endloc,pmap->flipped,pmap->len);
	  else
	    printf("Deleting Map[%d]:id=%lld,numsites[c=%d]=%d(MinSites=%d,MaxSites=%d),len=%0.3f(MinLen=%0.3f,MaxLen=%0.3f),startloc=%0.3f,endloc=%0.3f,TrueFlip=%d,len=%0.3f\n",
		   i,pmap->id,c,pmap->numsite[c],MinSites,tMaxSites, len,tMinLen,tMaxLen,pmap->startloc,pmap->endloc,pmap->flipped,pmap->len);
	  fflush(stdout);
	}
	if(DEBUG && (pmap!=NULL))
	  pmap->id = -pmap->id;
	continue;
      }

      if(MaxIntensity > 0 && BNXVersion >= 1 && pmap->AvgIntensity > MaxIntensity){
	if(DEBUG && (pmap!=NULL))
	  pmap->id = -pmap->id;
	continue;
      }

      if(MapPSFWidth && maxPSFWidth[c] > 0.0){
	int M = pmap->numsite[c];
	int J = 1;
	for(; J <= M; J++)
	  if(pmap->PSFWidth[c][J] > maxPSFWidth[c])
	    break;
	if(J <= M)
	  continue;
	if(colors==2 && !usecolor && maxPSFWidth[1] > 0.0){
	  M = pmap->numsite[1];
	  for(J=1; J <= M; J++)
	    if(pmap->PSFWidth[1][J] > maxPSFWidth[1])
	      break;
	  if(J <= M)
	    continue;
	}
      }
      
      if(j < i){
	Cmap *tmp = Map[j];
	pmap->mapid = j;
	Map[j] = pmap;
	tmp->mapid = i;
	Map[i] = tmp;
      }
      j++;
    }
  }

  if(VERB && j != nummaps){
    if(MapPSFWidth && maxPSFWidth[0]){
      if(tMaxSites || (!minSNRestimate && (MaxSiteDensity > 0.0 || MinSiteDensity > 0.0)))
	printf("Reduced number of input maps from %d to %d due to MinLen=%0.3f,MaxLen=%0.3f, MinSites=%d,MaxSites=%d,MinDensity=%0.2f,MaxDensity=%0.2f,MaxIntensity=%0.2f,maxPSFWidth[0]=%0.1f in file %s (%d of %d)\n",
	       total_lines / lines_per_map, j-orignummaps,tMinLen,tMaxLen,MinSites,tMaxSites, MinSiteDensity,MaxSiteDensity, MaxIntensity, maxPSFWidth[0], vfixx_filename[fileid], fileid+1, num_files);
      else
	printf("Reduced number of input maps from %d to %d due to MinLen=%0.3f,MaxLen=%0.3f, MinSites=%d,MaxIntensity=%0.2f,maxPSFWidth[0]=%0.1f, in file %s (%d of %d)\n",
	       total_lines / lines_per_map, j-orignummaps,tMinLen,tMaxLen,MinSites,MaxIntensity, maxPSFWidth[0], vfixx_filename[fileid], fileid+1, num_files);
    } else {
      if(tMaxSites || (!minSNRestimate && (MinSiteDensity > 0.0 || MaxSiteDensity > 0.0)))
	printf("Reduced number of input maps from %d to %d due to MinLen=%0.3f,MaxLen=%0.3f, MinSites=%d,MaxSites=%d,MinDensity=%0.2f,MaxDensity=%0.2f,MaxIntensity=%0.2f in file %s (%d of %d)\n",
	       total_lines / lines_per_map, j-orignummaps,tMinLen,tMaxLen,MinSites,tMaxSites, MinSiteDensity,MaxSiteDensity,MaxIntensity, vfixx_filename[fileid], fileid+1, num_files);
      else
	printf("Reduced number of input maps from %d to %d due to MinLen=%0.3f,MaxLen=%0.3f, MinSites=%d,MaxIntensity=%0.2f in file %s (%d of %d)\n",
	       total_lines / lines_per_map, j-orignummaps,tMinLen,tMaxLen,MinSites,MaxIntensity, vfixx_filename[fileid], fileid+1, num_files);
    }
    if(ids){
      printf("ERROR: filtering maps in refineA that were not filtered during Assembly will remove maps that are in the assembled contigs and hence refineA will fail\n");
      fflush(stdout);
      fflush(stdout);exit(1);
    }
    if(VERB/* HERE >=2 */ && j > orignummaps){
      double totlen = 0.0, totlen2 = 0.0;
      size_t nsites = 0, nsites2 = 0;
      for(int i= orignummaps; i < j; i++){
	Cmap *pmap = Map[i];
	maplist[i-orignummaps] = pmap;
	int M = pmap->numsite[0];
	totlen += pmap->site[0][M+1];
	nsites += M;
	if(colors >= 2){
	  int M2 = pmap->numsite[1];
	  totlen2 += pmap->site[1][M2 + 1];
	  nsites2 += M2;
	}
      }
      qsort(maplist, j-orignummaps, sizeof(Cmap*), (intcmp*)PmapLenDec);/* sort map pointers in descending order of length */
      double N50Sum= 0.0, N50Len = 0.0;
      for(int i = orignummaps; i < j; i++){
	Cmap *pmap = maplist[i-orignummaps];
	int M = pmap->numsite[0];
	if((N50Sum += pmap->site[0][M+1]) >= totlen * 0.5){
	  N50Len = pmap->site[0][M+1];
	  break;
	}
      }
      if(DEBUG && totlen > 0.0) assert(N50Len > 0.0);

      if(colors >= 2)
	printf("After applying -minlen etc to %s : maps=%d, sites=%lu,%lu, length= %0.3f kb (avg= %0.3f kb, label density1= %0.3f /100kb, label density2= %0.3f /100kb, N50= %0.4f Mb):wtime=%0.6f\n",
	       vfixx_filename[fileid], j-orignummaps,nsites,nsites2,totlen, 0.5 * (totlen+totlen2)/(j-orignummaps), nsites*100.0/max(0.001,totlen), nsites2*100.0/max(0.001,totlen2), N50Len*0.001,wtime());
      else
	printf("After applying -minlen etc to %s : maps=%d, sites=%lu, length= %0.3f kb (avg= %0.3f kb, label density= %0.3f /100kb, N50= %0.4f Mb):wtime=%0.6f\n",
	       vfixx_filename[fileid], j-orignummaps,nsites,totlen,totlen/(j-orignummaps), nsites*100.0/max(0.001,totlen), N50Len * 0.001,wtime());
    }
    fflush(stdout);
  }
  if(DEBUG>=2){
    for(int k = j; k < nummaps; k++){
      Cmap *pmap = Map[k];
      for(int c = 0; c < colors; c++){
	int numsite = pmap->blockmem ? pmap->blockmem - 2 : pmap->numsite[c];
	for(int i = 0; i <= numsite + 1;  i++)
	  pmap->site[c][i] = aNaN;
	if(MapSNR)
	  for(int i = 0; i <= numsite; i++)
	    pmap->SNR[c][i] = aNaN;
	if(MapIntensity)
	  for(int i = 0; i <= numsite; i++)
	    pmap->Intensity[c][i] = aNaN;
	if(MapPSFWidth)
	  for(int i = 1; i <= numsite; i++)
	    pmap->PSFWidth[c][i] = aNaN;
	if(MapImageLoc)
	  for(int i = 1; i <= numsite; i++){
	    pmap->ImageFOV[c][i] = -1;
	    pmap->ImageX[c][i] = aNaN;
	    pmap->ImageY[c][i] = aNaN;
	  }
      }
    }
  }
  nummaps = j;
  
  if(VERB>=2){
    printf("Finished filtering maps in input_bnx:fileid=%d:time=%0.6f(%0.6f)\n",fileid,mtime(),wtime());
    fflush(stdout);
  }

  if(BIAS_TRACE){/* print map id BIAS_TRACE_ID */
    for(int m = orignummaps; m < nummaps; m++){
      Cmap *pmap = Map[m];
      if(pmap->id == BIAS_TRACE_ID && pmap->SNR[0]){
	int N = pmap->numsite[0];
	double *SNR = pmap->SNR[0];
	FLOAT *Y =  pmap->site[0];
	printf("After filtering maps: id=%lld,N=%d:\n",pmap->id,N);
	for(int J = 1; J <= N+1; J++)
	  printf(" Y[%d]=%0.4f, SNR=%0.4f\n",J,Y[J], (J <= N) ? SNR[J] : -1.0);
	fflush(stdout);
      }
    }
  }

  if(nummaps == orignummaps){/* no maps from this file : reset QX code flags to original value */
    MapTrueSite = MapTrueSite_start;
    MapSNR = MapSNR_start;
    MapStitched = MapStitched_start;
    MapIntensity = MapIntensity_start;
    MapStitchLoc = MapStitchLoc_start;
    MapPSFWidth = MapPSFWidth_start;
    MapImageLoc = MapImageLoc_start;
  } else if(!QXerror){/* keep track of which QX codes are present in any previous file with maps */
    MapTrueSite |= MapTrueSite_start;
    MapSNR |= MapSNR_start;
    MapStitched |= MapStitched_start;
    MapIntensity |= MapIntensity_start;
    MapStitchLoc |= MapStitchLoc_start;
    MapPSFWidth |= MapPSFWidth_start;
    MapImageLoc |= MapImageLoc_start;
  }

  /* check that all maps are monotonic and extend ends so they are at least ENDEXTEND */
  int nthreads = max(1,min(numthreads, nummaps-orignummaps));

  #pragma omp parallel num_threads(nthreads) if(nthreads > 1)
  {  
    #pragma omp for schedule(dynamic,512)
    for(int i=orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      if(DEBUG) assert(pmap->id > 0);

      if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID && colors <= 1){
        printf("id=%lld: pmap->len = %0.6f, site[numsite, numsite+1] = %0.6f, %0.6f\n",pmap->id,pmap->len, pmap->site[0][pmap->numsite[0]], pmap->site[0][pmap->numsite[0]+1]);
	fflush(stdout);
      }

      for(int c=0; c < colors; c++){
	FLOAT *X = pmap->site[c];
	int N = pmap->numsite[c];
	assert(X[0] == 0.0);
	for(int k=1; k < N;k++)
	  if(X[k+1] < X[k]){
	    #pragma omp critical
	    {
	      printf("map %d (id=%lld), color %d:site%d = %0.6f, site%d = %0.6f (out of order) in %s (len=%0.6f,N=%d)\n",
		     i, pmap->id, c+1, k, X[k], k+1, X[k+1], vfixx_filename[pmap->fileid], pmap->len,N);
	      for(int t = 0; t <= N+1; t++)
		printf("X[%d] = %0.6f\n",t,X[t]);
	      fflush(stdout);
	      assert(0);
	      fflush(stdout);exit(1);
	    }
	  }
	/*      if(X[1] < 0.0){
		fprintf(stderr,"map id %d, color %d: site1 = %0.3f (leftmost site is located below 0.0)\n", 
		pmap->id, c+1, X[1]);
		fflush(stdout);exit(1);
		}*/
	if(bnxerror && X[N+1] < X[N]){
	  #pragma omp critical
	  {
	    printf("map %d (id=%lld), color %d: site%d = %0.6f, len = %0.7f (rightmost site is to the right of the end) in %s (len=%0.6f,N=%d)\n",
		   i, pmap->id, c+1, N, X[N], X[N+1], vfixx_filename[pmap->fileid], pmap->len, N);
	    fflush(stdout);
	    assert(0);
	    fflush(stdout);exit(1);
	  }
	}
	if(X[1] < ENDEXTEND){
	  register double delta = ENDEXTEND - X[1];
	  for(register int k=1; k <= N+1; k++)
	    X[k] += delta;
	}
	if(X[N+1] < X[N] + ENDEXTEND)
	  X[N+1] = X[N] + ENDEXTEND;
      } /* c = 0 .. colors-1 */

      if(colors >= 2){/* equalize length for all colors */
	FLOAT len = 0.0;
	for(int c = 0; c < colors; c++){
	  int M = pmap->numsite[c];
	  len = max(len, pmap->site[c][M+1]);
	}
	for(int c = 0; c < colors; c++){
	  int M = pmap->numsite[c];
	  pmap->site[c][M+1] = len;
	}
	pmap->len = len;
      } else
	pmap->len = pmap->site[0][pmap->numsite[0]+1];
      
      if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
        printf("id=%lld: After -minEND %0.3f : pmap->len -> %0.6f, site[0][numsite, numsite+1] = %0.6f, %0.6f\n",pmap->id, ENDEXTEND, pmap->len, pmap->site[0][pmap->numsite[0]], pmap->site[0][pmap->numsite[0]+1]);
	if(colors >= 2)
	  printf("\t site[1][numsite,numsite+1] = %0.6f, %0.6f\n",pmap->site[1][pmap->numsite[1]], pmap->site[1][pmap->numsite[1]+1]);
	fflush(stdout);
      }
    } // for
  }// parallel

  if(VERB>=2){
    printf("input_bnx:Finished checking maps to check for out of order sites:fileid=%d:time=%0.6f(%0.6f)\n",fileid,mtime(),wtime());
    fflush(stdout);
  }

  if(BIAS_TRACE){/* print map id BIAS_TRACE_ID */
    for(int m = orignummaps; m < nummaps; m++){
      Cmap *pmap = Map[m];
      if(pmap->id == BIAS_TRACE_ID && pmap->SNR[0]){
	int N = pmap->numsite[0];
	double *SNR = pmap->SNR[0];
	FLOAT *Y =  pmap->site[0];
	printf("After checking map ends: id=%lld,N=%d:\n",pmap->id,N);
	for(int J = 1; J <= N+1; J++)
	  printf(" Y[%d]=%0.4f, SNR=%0.4f\n",J,Y[J], (J <= N) ? SNR[J] : -1.0);
	fflush(stdout);
      }
    }
  }

  if(DEBUG) assert(&Map != &refmap);

  if(DEBUG>=2){ /* check all maps to make sure all colors are present */
    #pragma omp parallel num_threads(nthreads) if(nthreads > 1)
    {  
      #pragma omp for schedule(static,512)
      for(int i=orignummaps; i < nummaps;i++){
	pmap = Map[i];
	for(register int c=0; c < colors; c++)
	  if(!pmap->site[c]){
            #pragma omp critical
	    {
	      printf("Color %d missing from map[%d] with id=%lld in file %s\n",c+1,i-orignummaps,pmap->id,filename);
	      fflush(stdout);exit(1);
	    }
	  }
      }
    }
    if(VERB>=2){
      printf("input_bnx : Finished checking for missing color data :fileid=%d:time=%0.6f(%0.6f)\n",fileid,mtime(),wtime());
      fflush(stdout);
    }
  }

  if(nummaps - orignummaps > 1)
    SplitMap = 0;

  if(LEAKDEBUG){
    delete [] RunDataIdmap; RunDataIdmap = NULL;
    RunDataIdMax = 0;
  }

  delete [] maplist;

  return nummaps-orignummaps;
}


int merrorexit = 0;// if any input file had an error

int errorexit = 0;

static void cmap_error(long long &contig, long long &lastcontig, FILE * &fp, int &linecnt, char *buf, char *filename, char * &pt, char * &qt, int err = 1)
{
  /* handle error in current map by skipping to next map */
  fflush(stdout);
  if(err)
    errorexit++;

  /* jump to start of next map */
  linecnt++;
  lastcontig = contig;
  for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++) {   int len;
    if(VERB>=2){
      printf("Non-Header line %d in %s:\n%s\n",linecnt,filename,buf);
      fflush(stdout);
    }

    if((len = strlen(buf)) >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      if(len >= LINESIZ-1)
	printf("Line too long (max length=%d, increase LINESIZE in constants.h) in input file %s on line %d:\n%s\n",LINESIZ-1,filename,linecnt,buf);
      else
	printf("Line not correctly terminated (last char= %c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      fflush(stdout);exit(1);
    }
    pt = buf;
    contig = strtoll(pt,&qt,10);
    if(contig >= 0 && contig != lastcontig){
      if(VERB>=2){
	printf("Skipping to next map from contig=%lld to contig=%lld\n", lastcontig,contig);
	fflush(stdout);
      }
      lastcontig = -1;
      return;
    }
  }
  if(errorexit){
    printf("%d maps have errors in %s\n", errorexit, filename);
    fflush(stdout);exit(1);
  }
}

static double cumlen = 0.0;

char *lastfilename = NULL;/* pointer to last input file by input_cmap() */

void check_cmap(int fileid, char *filename, int orignummaps, int &nummaps, int &maxmaps, Cmap **&Map, int HapPresent)
{
  if(VERB>=2){
    if(filename)
      printf("nummaps=%d:After parsing %s\n",nummaps,filename);
    else
      printf("nummaps=%d:After parsing all input files\n",nummaps);
    for(int i = orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      if(colors >= 2)
	printf("Map[%d]:M=%d,%d,Len=%0.3f,%0.3f\n",i,pmap->numsite[0],pmap->numsite[1],pmap->site[0][pmap->numsite[0]+1],pmap->site[1][pmap->numsite[1]+1]);
      else
	printf("Map[%d]:M=%d,Len=%0.3f\n",i,pmap->numsite[0],pmap->site[0][pmap->numsite[0]+1]);
      fflush(stdout);
    }
  }

  int numthreads = 1;

#ifdef _OPENMP
  numthreads = omp_get_max_threads();
  numthreads = min(numthreads,MaxThreads);
  numthreads = max(1,min(numthreads, (nummaps-orignummaps)/1024));
#endif
  
  /* check that all maps sites are monotonic */
  #pragma omp parallel for schedule(dynamic,1024) num_threads(numthreads) if(numthreads > 1)
  for(int i= orignummaps; i < nummaps; i++){
    Cmap *pmap = Map[i];
    if(DEBUG && pmap->contig && pmap->contig->HapSite[0]){
      int m = pmap->contig->numsite[0];
      int *HapSite = pmap->contig->HapSite[0];
      if((DEBUG && !(HapSite[0]==0 /* WAS && HapSite[m+1]==0 */)) || (VERB>=2 && pmap->id==877)){
        printf("WARNING:Map[%d]:id=%lld,contig->numsite[0]=m=%d,HapSite[0]=%d,HapSite[m+1]=%d\n",
          i,pmap->id,m,HapSite[0],HapSite[m+1]);
	fflush(stdout);
	// assert(HapSite[0]==0 && HapSite[m+1]==0);
      }
      HapSite[0] = HapSite[m+1] = 0;
    }
    for(int c=0; c < colors; c++){
      FLOAT *X = pmap->site[c];
      int N = pmap->numsite[c];
      assert(X[0] == 0.0);
      for(int k=1; k < N;k++)
	if(X[k+1] < X[k]){
	  #pragma omp critical
	  {
	    printf("map id %lld, color %d:site%d = %0.3f, site%d = %0.3f (out of order)\n",
		    pmap->id, c+1, k, X[k], k+1, X[k+1]);
	    fflush(stdout);exit(1);
	  }
	}
      if(X[1] < 0.0){
        #pragma omp critical
	{
	  printf("map id %lld, color %d: site1 = %0.3f (leftmost site is located below 0.0)\n", 
		  pmap->id, c+1, X[1]);
	  fflush(stdout);exit(1);
	}
      }
      if(X[N+1] < X[N]){
        #pragma omp critical
	{
	  printf("map id %lld, color %d: site%d = %0.3f, len = %0.3f (rightmost site is to the right of the end)\n",
		  pmap->id, c+1, N, X[N], X[N+1]);
	  fflush(stdout);exit(1);
	}
      }
      if(X[1] < ENDEXTEND){
	if(VERB/* HERE >=2 */){
	  #pragma omp critical
	  {
	    printf("WARNING:map id %lld, color %d: site1 = %0.6f (leftmost site is below %0.6f : extending left end)\n",pmap->id,c+1,X[1],ENDEXTEND);
	    fflush(stdout);
	  }
	}
	register double delta = ENDEXTEND - X[1];
	for(register int k=1; k <= N+1; k++)
	  X[k] += delta;
      }
      if(X[N+1] < X[N] + ENDEXTEND){
	if(VERB/* HERE >=2 */ && X[N+1] < X[N] + ENDEXTEND*0.99){
	  #pragma omp critical
	  {
	    printf("WARNING:map id %lld, color %d: site%d = %0.6f, len = %0.6f (rightmost site is too close to the right end : extending right end at least %0.6f past rightmost site)\n",
		   pmap->id, c+1, N, X[N], X[N+1],ENDEXTEND);
	    fflush(stdout);
	  }
	}

	X[N+1] = X[N] + ENDEXTEND;
      }
    }
    if(VERB>=2 && colors >= 2){
      printf("Map[%d]:M=%d,%d,Len=%0.3f,%0.3f\n",i,pmap->numsite[0],pmap->numsite[1],pmap->site[0][pmap->numsite[0]+1],pmap->site[1][pmap->numsite[1]+1]);
      fflush(stdout);
    }

    if(colors >= 2){  /* equalize lengths of 2 colors */
      FLOAT len = 0.0;
      for(int c = 0; c < colors; c++){
        int M = pmap->numsite[c];
        len = max(len, pmap->site[c][M+1]);
      }
      for(int c = 0; c < colors; c++){
        int M = pmap->numsite[c];
	pmap->site[c][M+1] = len;
      }
    }
  }

  if(VERB>=2 && colors >= 2){
    printf("After checking sites are monotonic\n");
    for(int i = orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      printf("Map[%d]:M=%d,%d,Len=%0.3f,%0.3f\n",i,pmap->numsite[0],pmap->numsite[1],pmap->site[0][pmap->numsite[0]+1],pmap->site[1][pmap->numsite[1]+1]);
      fflush(stdout);
    }
  }

  if(DEBUG && colors >= 2){/* check molecule lengths are the same */
    #pragma omp parallel for schedule(dynamic,1024)
    for(int i= orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      if(VERB>=2){
	printf("Map[%d]:M=%d,%d,Len=%0.3f,%0.3f\n",i,pmap->numsite[0],pmap->numsite[1],pmap->site[0][pmap->numsite[0]+1],pmap->site[1][pmap->numsite[1]+1]);
	fflush(stdout);
      }
      FLOAT len1 = pmap->site[0][pmap->numsite[0]+1];
      FLOAT len2 = pmap->site[1][pmap->numsite[1]+1];
      if(!(fabs(len1 - len2) < 1e-3)){
	#pragma omp critical
	{
	  printf("Map[%d] : id=%lld, M1=%d,len1=%0.4f, M2=%d,len2=%0.4f\n",i, pmap->id, pmap->numsite[0], len1, pmap->numsite[1], len2);
	  fflush(stdout);
	  assert(fabs(len1 - len2) < 1e-3);
	}
      }
    }
  }

  int tMaxSites = (fileid < num_files && minSNRestimate) ? MaxSites2 : MaxSites;    
  if(fileid < num_files){
    if(DEBUG) assert(!minSNRestimate);

    /* filter out maps smaller than MinLen or larger than MaxLen or with fewer sites than MinSites or greater than MaxSites */
    int c = 0;
    if(colors == 2 && usecolor==2)
      c = 1;
    register int i,j;
    for(i=j=orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      int nsites = pmap->numsite[c];
      if(colors==2 && !usecolor)
	nsites += pmap->numsite[1-c];
      double len = pmap->site[0][pmap->numsite[0]+1];
      if(pmap->numsite[0] < MinSites || (tMaxSites && pmap->numsite[0] > tMaxSites) || 
	 len < MinLen || (MaxLen > 0.0 && len > MaxLen) || 
	 (MaxSiteDensity > 0.0 && nsites * 100.0 > MaxSiteDensity * len * colors) ||
	 (MinSiteDensity > 0.0 && nsites * 100.0 < MinSiteDensity * len * colors)){
	if(VERB>=2){
	  printf("Deleting map[%d]:id=%lld,numsites=%d(MinSites=%d,MaxSites=%d),len=%0.3f(MinLen=%0.3f,MaxLen=%0.3f),PixelLen=%0.6f/%0.6f,startloc=%0.3f,endloc=%0.3f,len=%0.3f\n",
		 i,pmap->id,pmap->numsite[0],MinSites,tMaxSites,pmap->site[0][pmap->numsite[0]+1],MinLen,MaxLen,PixelLen,origPixelLen,pmap->startloc,pmap->endloc,pmap->len);
	  fflush(stdout);
	}
	if(DEBUG)
	  pmap->id = -pmap->id;
	continue;
      }
      if(j < i){
	Cmap *tmp = Map[j];
	if(DEBUG) assert(tmp != NULL);
	pmap->mapid = j;
	Map[j] = pmap;
	tmp->mapid = i;
	Map[i] = tmp;
      }
      j++;
    }
    if(VERB && j != nummaps){
      printf("Reduced number of input maps from %d to %d due to MinLen=%0.3f,MaxLen=%0.3f, MinSites=%d,MaxSites=%d in file %s (%d of %d)\n",
	     nummaps-orignummaps,j-orignummaps,MinLen,MaxLen,MinSites,tMaxSites,vfixx_filename[fileid], fileid+1, num_files);
      fflush(stdout);
    }
    nummaps = j;
  }

  if(fileid < num_files && (mres >= 0.0 || MaxInterval > 0.0)){/* mresSD = 0 : merge map sites at distance mres or closer
								  mresSD > 0 : merge map sites at distance X with probability 0.5(1-erf((X-mres)/(mresSD*sqrt(2)))) */
    if(VERB>=2){
      printf("input_cmap.cpp : Calling mapcorrect(BNX=0) : MaxInterval= %0.3f (ref=%d), num_reffiles=%d, mres= %0.4f\n",MaxInterval,refMaxInterval,num_reffiles,mres);
      fflush(stdout);
    }
    if(DEBUG && (refmap || Map)) assert(&Map != &refmap);

    mapcorrect(Map, orignummaps, nummaps, maxmaps, fileid, 0);
  }

  if(fileid >= num_files && MaxInterval > 0.0 && refMaxInterval){
    double origmres = mres;
    mres = 0.0;

    if(VERB>=2){
      printf("input_cmap.cpp : Calling mapcorrect(BNX=0) : MaxInterval= %0.3f (ref=%d), num_reffiles=%d, mres= %0.4f\n",MaxInterval,refMaxInterval,num_reffiles,mres);
      fflush(stdout);
    }
    if(DEBUG && (refmap || Map)) assert(&Map == &refmap);

    mapcorrect(Map, orignummaps, nummaps, maxmaps, fileid, 0);

    mres = origmres;

    for(int i=orignummaps; i < nummaps; i++){/* copy of loop below, but only allocates/computes remap[] and NickCnt[] */
      Cmap *pmap = Map[i];
      /* Initialize refname with filename */
      pmap->refname = strdup(filename);

      /* initialize and copy original site locations (before resolution based condensation) */
      for(int c=0; c < colors; c++){
	if(pmap->remap[c]) delete [] pmap->remap[c];
	pmap->remap[c] = new int[pmap->numsite[c]+1];

	if(pmap->NickCnt[c]) delete [] pmap->NickCnt[c];
	pmap->NickCnt[c] = new int[pmap->numsite[c]+1];

	for(int J=1; J <= pmap->numsite[c]; J++){
	  pmap->remap[c][J] = J;
	  pmap->NickCnt[c][J] = 1;
	}
      }
    }    

  } else if(fileid >= num_files){
    if(DEBUG) assert( &Map == &refmap || svcheck || svcompare1);

    /* HERE : convert to multithreaded version in mapcorrect() so orig* are always block allocated (and reduce memory leak) */
    for(int i=orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      /* Initialize refname with filename */
      pmap->refname = strdup(filename);

      /* initialize and copy original site locations (before resolution based condensation) */
      for(int c=0; c < colors; c++){
	pmap->orignumsite[c] = pmap->numsite[c];

	if(pmap->origsite[c])  delete [] pmap->origsite[c];
	pmap->origsite[c] = new double[pmap->numsite[c]+2];
	for(int J=0; J <= pmap->numsite[c]+1; J++)
	  pmap->origsite[c][J] = pmap->site[c][J];

	if(pmap->remap[c]) delete [] pmap->remap[c];
	pmap->remap[c] = new int[pmap->numsite[c]+1];

	if(pmap->NickCnt[c]) delete [] pmap->NickCnt[c];
	pmap->NickCnt[c] = new int[pmap->numsite[c]+1];

	for(int J=1; J <= pmap->numsite[c]; J++){
	  pmap->remap[c][J] = J;
	  pmap->NickCnt[c][J] = 1;
	}
      }
    }
  }

  if(fileid >= num_files && rres > 0.0){/* merge reference sites closer than rres*0.5 */
    if(DEBUG) assert( &Map == &refmap || svcheck || svcompare1);
    double resKB = rres * 0.500;
    int orignumsites = 0;
    int numsites = 0;
    for(int i=orignummaps; i < nummaps; i++){
      Cmap *pmap = Map[i];
      double reflen = pmap->site[0][pmap->numsite[0]+1];
      for(int c=0; c < colors; c++)
	orignumsites += pmap->numsite[c];

      if(VERB>=2 && &Map == &refmap){
	printf("refmap[%d]: numsite[0]= %d\n",i,pmap->numsite[0]);
	for(int t = 1; t <= pmap->numsite[0]+1; t++)
	  printf("\t t=%d: site[0]= %0.8f (delta= %0.8f)\n", t, pmap->site[0][t], pmap->site[0][t] - pmap->site[0][t-1]);
        fflush(stdout);
      }

      if(!HapPresent)
	deres(pmap,resKB);

      for(int c=0; c < colors; c++){
	numsites += pmap->numsite[c];
	pmap->site[c][pmap->numsite[c]+1] = reflen;
      }
    }
    if(VERB>=2 && orignumsites != numsites){
      printf("reduced sites from %d to %d in %d ref maps due to rres=%0.4f (resKB= %0.4f), SDrange=%0.2f,PixelLen=%0.6f/%0.6f (file = %s), HapPresent=%d\n",
	     orignumsites,numsites, nummaps-orignummaps, rres, resKB, SDrange, PixelLen, origPixelLen, vfixx_filename[fileid],HapPresent);

      if(&Map == &refmap){
	for(int i = orignummaps; i < nummaps; i++){
	  Cmap *pmap = Map[i];
	  printf("refmap[%d]: numsite[0]= %d\n",i,pmap->numsite[0]);
	  for(int t = 1; t <= pmap->numsite[0]+1; t++)
	    printf("\t t=%d: site[0]= %0.8f (delta= %0.8f)\n", t, pmap->site[0][t], pmap->site[0][t] - pmap->site[0][t-1]);
	}
      }
      fflush(stdout);
    }
  }

  /* check all maps to make sure all colors are present */
  for(int i=orignummaps; i < nummaps;i++){
    Cmap *pmap = Map[i];
    for(int c=0; c < colors; c++)
      if(!pmap->site[c]){
	printf("check_cmap:Color %d missing from map[%d] with id=%lld in file %s(fileid=%d,orignummaps=%d,nummaps=%d)\n",c+1,i-orignummaps,pmap->id,filename,fileid,orignummaps,nummaps);
	fflush(stdout);exit(1);
      }
  }

  if(VERB >=2 && colors==1 && (fileid >= num_files || nummaps == orignummaps+1)){
    for(int i = orignummaps; i < nummaps; i++){
      printf("input_cmap:previous maps=%d(total len=%0.3f):fileid=%d/%d,filename=%s:map %d/%d with id=%lld,numsites=%d,len=%0.3f (site density= %0.2f per 100kb)\n",
	     orignummaps,cumlen,fileid,num_files,filename, i+1-orignummaps, nummaps-orignummaps, Map[i]->id, Map[i]->numsite[0],Map[i]->site[0][Map[i]->numsite[0]+1], 
	     Map[i]->numsite[0]*100.0/Map[i]->site[0][Map[i]->numsite[0]+1]);
      cumlen += Map[i]->site[0][Map[i]->numsite[0]+1];
    }
    fflush(stdout);
  }
}

/* NOTE hmap=1 means the filename suffix was .hmap : HapSite information can be present even in .cmap files */
int input_cmap(int fileid,int &nummaps, int &maxmaps, Cmap **&Map, int hmap, bool multithreaded, char *mbuf)
{
  char *filename = vfixx_filename[fileid];
  char *pt, *qt;
  FILE *fp;
  
  int orignummaps = nummaps;// NOTE : not used in multithreaded calling of input_cmap()
  int nmapcnt = 0;

  int MaskPresent = 0;
  int SNRpresent = 0;
  int HapPresent = 0;/* If HapSite & HapDelta are present in header */

  Cmap *pmap=0;

  int NColors= -1,NFragments= -1;

  int linecnt = 1;

  int origcolors = colors;

  assert(colors >= 1);

  for(int c = 0; c < max(2,colors); c++)
    if(!Nickase[c]){
      Nickase[c] = strdup((char *)"unknown");
      if(VERB>=2){
	printf("input_cmap(%s):Setting global Nickase[%d]=%p(%s)\n",filename,c,Nickase[c],Nickase[c]);
	fflush(stdout);
      }
    }

  double lastlen = 0.0;
  long long int lastcontig = -1;
  int lastNumSites = 0;
  int lastSiteID = 0;

  errorexit = 0;/* remember if a fatal error was encountered with any map */

  int ChimQuality = 0;
  int MapTrueSite = 0;

  int len;
  long long int contig = -1;

  int VERSION = 0;

  if( &Map == &Gmap && fileid < num_files){
    if(maptype == 0 && MapType == -1){
      printf("input_cmap:fileid=%d/%d, file %s:previous input map type was .bnx or .vfixx : cannot be mixed with .cmap\n",fileid, num_files, filename);
      fflush(stdout);assert(0);
      fflush(stdout);exit(1);
    }
    maptype = 1;
  }

  char *lastline = NULL;

  if(mbuf == NULL){
    if(DEBUG) assert(!multithreaded);
    mbuf = &buf[0];
  }

  if(fileid < num_files){
    if(DEBUG) assert(&Map != &refmap);

    if(!origPixelLen){
      #pragma omp critical(fileinput) // open file and read header lines in critical section
      {
        if(!origPixelLen){
          origPixelLen = 0.500;/* default value for .cmap files */
	  if(VERB>=2){
            printf("input_cmap:origPixelLen = %0.4f\n",origPixelLen);
	    fflush(stdout);
          }
        }
      } // #pragma omp critical
    } 

    if(origPixelLen && fabs(origPixelLen - 0.500) > 0.001){
      printf("WARNING: Input CMAP file %s assumes PixelLen = 0.500, differs from previous input bpp=%0.6f from %s\n", 
	     filename, origPixelLen*1000.0, lastfilename ? lastfilename : "unknown input file");
      fflush(stdout);
    }
  }

  if(VERB>=2){
    printf("input_cmap:fileid=%d,num_files=%d,num_reffiles=%d,filename=%s(origPixelLen=%0.4f)\n",fileid,num_files,num_reffiles,filename,origPixelLen);
    fflush(stdout);
  }

  errno = 0;
  if((fp = Fopen(filename,"r",NFSdelay))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("input_cmap:Failed to read input file %s (%d'th of %d files):errno=%d:%s\n",filename, fileid+1,num_files,eno,err);
    fflush(stdout);exit(1);
  }

  for(;(lastline = fgets(mbuf,LINESIZ,fp)) != NULL; linecnt++) {
    if((len = strlen(mbuf)) >= LINESIZ-1 || (mbuf[len-1] != '\n' && mbuf[len-1] != '\r')){
      if(len >= LINESIZ-1)
	printf("Line too long (exceeds max length= %d characters, increase LINESIZ in constants.h) in input file %s on line %d:\n%s\n",LINESIZ-1,filename,linecnt,mbuf);
      else
	printf("Line not correctly terminated (last char= %c) in input file %s on line %d:\n%s\n",mbuf[len-1],filename,linecnt,mbuf);
      fflush(stdout);
      (void)fclose(fp);

      merrorexit++;
      break;// will return 0 outside critical section
    }
    
    if(mbuf[0] == '#'){/* most comment lines are ignored, except "Label Channels:", "Number of Consensus Maps:" and "Nickase Recognition Site <color>:" */
      if(VERB>=2){
	printf("Header line %d in %s (colors= %d, origcolors= %d,refcolors= %d, colors_set=%d, usecolor= %d):\n%s\n",linecnt,filename,colors,origcolors,refcolors,colors_set,usecolor,mbuf);
	fflush(stdout);
      }

      char *key = (char *)"CMAP File Version:";
      if((pt=strstr(&mbuf[1],key))){
	pt += strlen(key);
	double Version = strtod(pt,&qt);
	if(qt == pt){
	  printf("Line %d with %s in file %s has invalid Version value:\n%s\n",linecnt,key,filename,mbuf);
	  fflush(stdout);exit(1);
	}
	if(fabs(Version - 0.1) <= 0.00001)
	  VERSION = 1;
	else if(fabs(Version - 0.2) <= 0.00001)
	  VERSION = 2;
	else {
	  printf("CMAP File Version %0.1f in %s not supported (Only 0.1 or 0.2 are supported)\n",Version,filename);
	  fflush(stdout);exit(1);
	}
	continue;
      }
      
      key = (char *)"Label Channels:";
      if((pt=strstr(&mbuf[1],key))){
	pt += strlen(key);
	NColors = strtol(pt,&qt,10);
	if(qt == pt){
	  printf("Line %d with %s in file %s has invalid integer value:\n%s\n",linecnt,key,filename,mbuf);
	  fflush(stdout);exit(1);
	}
	if(fileid >= num_files){
	  if(!refcolors)
	    refcolors = NColors;
	  else if(refcolors != NColors){
	    printf("Line %d in file %s with %s has value %d which is different from value %d in previous -ref input file %s\n%s\n",
		   linecnt,filename,key,NColors,refcolors,vfixx_filename[num_files],mbuf);
	    printf("All -ref input files must have the same number of Label colors\n");
	    fflush(stdout);
	    if(DEBUG) assert(fileid > num_files);
	    exit(1);
	  }
	}

	if(!colors_set){
          #pragma omp critical(fileinput) // open file and read header lines in critical section
	  {
	    if(!colors_set){/* rare condition due to thread race can cause colors_set = 1 : in that case drop to 2nd case below */
	      for(int c = colors; c < NColors; c++)
		if(!Nickase[c]){
		  Nickase[c] = strdup((char *)"unknown");
		  if(VERB>=2){
		    printf("input_cmap(%s):Line %d: Setting global Nickase[%d]=%p(%s)\n",filename,linecnt,c,Nickase[c],Nickase[c]);
		    fflush(stdout);
		  }
		}
	      colors = NColors;
	      colors_set = 1;
	      if(colors == 1 && usecolor==1)
		usecolor = 0;
	      if(VERB){
		printf("Line %d in %s: setting colors=%d, colors_set=%d, usecolor= %d due to Label Channels: %d:\n%s\n",linecnt,filename,colors,colors_set,usecolor,NColors,mbuf);
		fflush(stdout);
	      }
	    }
	  }  // pragma omp critical(fileinput)
	}

	if(colors_set){
	  if(NColors != colors){
	    if(DEBUG) assert(fileid > 0);
	    if(usecolor && NColors == 2 && fileid >= num_files){
	      if(VERB/* HERE HERE >=2 */){
		printf("input_cmap(%s):Line %d:NColors=%d, usecolor=%d: Will read in Reference as 2-color and reduce to 1-color later\n",filename,linecnt,NColors,usecolor);
		fflush(stdout);
	      }
	      colors = 2;/* read in Reference as 2-color then reduce to 1-color later */
	    } else {
	      printf("Line %d in file %s with %s has value %d which is different from value %d in previous input file %s:\n%s\n",linecnt,filename,key,NColors,colors,lastfilename,mbuf);
	      if(colors==2){
		if(fileid == num_files)
		  printf("Please rerun with -usecolor 1 OR -usecolor 2\n");
		else if(fileid > num_files)
		  printf("All -ref input files must have the same number of Label colors\n");
	      }
	      fflush(stdout);exit(1);
	    }
	  }
	}

	continue;
      }

      key = (char *)"h CMapId";
      if((pt = strstr(&mbuf[1],key))){
	pt += strlen(key);
	if((qt = strstr(pt, "ChimQuality"))){
	  if(VERB>=2){
	    printf("Found ChimQuality in header line of %s\n", filename);
	    fflush(stdout);
	  }
	  ChimQuality = 1;
	  if((qt = strstr(pt,"ChimNorm"))){
	    if(VERB>=2 && ChimQuality < 2){
	      printf("Found ChimNorm in header line of %s\n", filename);
	      fflush(stdout);
	    }
	    ChimQuality = 2;
	    if(!(qt = strstr(pt,"SegDupL"))){
	      printf("Found ChimQuality and ChimNorm but NOT SegDupL in header line of %s: This CMAP not supported\n%s\n", filename, pt);
	      fflush(stdout); exit(1);
	    }
	    if(!(qt = strstr(pt,"SegDupR"))){
	      printf("Found ChimQuality and ChimNorm but NOT SegDupR in header line of %s: This CMAP not supported\n%s\n", filename, pt);
	      fflush(stdout); exit(1);
	    }
	    if(!(qt = strstr(pt,"FragileL"))){
	      printf("Found ChimQuality and ChimNorm but NOT FragileL in header line of %s: This CMAP not supported\n%s\n", filename, pt);
	      fflush(stdout); exit(1);
	    }
	    if(!(qt = strstr(pt,"FragileR"))){
	      printf("Found ChimQuality and ChimNorm but NOT FragileR in header line of %s: This CMAP not supported\n%s\n", filename, pt);
	      fflush(stdout); exit(1);
	    }
	    if(!(qt = strstr(pt,"OutlierFrac"))){
	      printf("Found ChimQuality and ChimNorm but NOT OutlierFrac in header line of %s: This CMAP not supported\n%s\n", filename, pt);
	      fflush(stdout); exit(1);
	    }

	    if((qt = strstr(pt,"MolSd"))){
	      if(VERB>=2 && ChimQuality < 3){
		printf("Found MolSd in header line of %s\n",filename);
		fflush(stdout);
	      }
	      ChimQuality = 3;
	      if(!(qt = strstr(pt,"ExpSd"))){
		// HERE HERE 		printf("Found MolSd but NOT ExpSd in header line of %s: This CMAP not supported\n%s\n",filename,pt);
		// HERE HERE 		fflush(stdout);exit(1);
	      } else
		ChimQuality = 4;
	      if(!(qt = strstr(pt,"MolCov"))){
		printf("Found MolSd but NOT MolCov in header line of %s: This CMAP not supported\n%s\n",filename,pt);
		fflush(stdout);exit(1);
	      }
	      if(!(qt = strstr(pt,"ChiSq"))){
		printf("Found MolSd but NOT ChiSq in header line of %s: This CMAP not supported\n%s\n",filename,pt);
		fflush(stdout);exit(1);
	      }
	    }
	  }
	}

	if((qt = strstr(pt,"Mask"))){
	  if(VERB>=2){
	    printf("Found Mask in header line of %s\n", filename);
	    fflush(stdout);
	  }
	  MaskPresent = 1;
	}

	if((qt = strstr(pt,"HapSite"))){
          if(VERB>=2){
	    printf("Found HapSite in header line %d of %s\n", linecnt, filename);
	    fflush(stdout);
	  }
	  HapPresent = 1;
	}

	if((qt = strstr(pt,"SNRgmean"))){
	  if(VERB>=2){
	    printf("Found SNRgmean in header line %d of %s\n", linecnt, filename);
	    fflush(stdout);
	  }
	  SNRpresent = 1;
	}
	continue;
      }

      key = (char *)"Original QueryMap:";
      if((pt=strstr(&mbuf[1],key))){
        pt += strlen(key);
	while(isspace(mbuf[len-1]))
	  mbuf[--len] = '\0';/* strip off white space & CR/LF from end of line */

       #pragma omp critical(fileinput) // open file and read header lines in critical section
       {
	if(VERB>=2){
	    printf("Found Original Querymap on line %d in %s (fileid=%d,num_files=%d)\n%s\n",linecnt,filename,fileid,num_files,pt);
	    fflush(stdout);
	}
	if(fileid < num_files){
	  if(QueryOrigin && strcmp(QueryOrigin,pt)){
	    printf("Original QueryMap=%s on line %d in %s does not match previous value of %s (in previous -i input files)\n",pt,linecnt,filename,QueryOrigin);
	    fflush(stdout);exit(1);
	  }

	  QueryOrigin = strdup(pt);
	  MapTrueSite = 1;
	} else {
	  if(RefOrigin && strcmp(RefOrigin,pt)){
	    printf("Original QueryMap=%s on line %d in %s does not match previous value of %s (in previous -ref input files)\n",pt,linecnt,filename,QueryOrigin);
	    fflush(stdout);exit(1);
	  }

	  RefOrigin = strdup(pt);
	  RefTrueSite = 1;
	  if(VERB>=2){
	    printf("Read Original QueryMap as RefOrigin = %s\n",RefOrigin);
	    fflush(stdout);
	  }
	}
       } // pragma omp critical(fileinput)
	continue;
      }

      key = (char *)"Original RefMap:";
      if((pt=strstr(&mbuf[1],key))){
        pt += strlen(key);
	while(isspace(mbuf[len-1]))
	  mbuf[--len] = '\0';/* strip off white space & CR/LF from end of line */

       #pragma omp critical(fileinput) // open file and read header lines in critical section
       {
	if(fileid < num_files){
	  if(QueryOrigin && strcmp(QueryOrigin,pt)){
	    printf("Original QueryMap=%s on line %d in %s does not match previous value of %s (in previous -i input files)\n",pt,linecnt,filename,QueryOrigin);
	    fflush(stdout);exit(1);
	  }

	  QueryOrigin = strdup(pt);
	  MapTrueSite = 1;
	} else {
	  if(RefOrigin && strcmp(RefOrigin,pt)){
	    printf("Original QueryMap=%s on line %d in %s does not match previous value of %s (in previous -ref input files)\n",pt,linecnt,filename,QueryOrigin);
	    fflush(stdout);exit(1);
	  }

	  RefOrigin = strdup(pt);
	  RefTrueSite = 1;
	}
       } // pragma omp critical(fileinput)
	continue;
      }

      key = (char *)"Number of Consensus Maps:";
      if((pt=strstr(&mbuf[1],key))){
	pt += strlen(key);
	NFragments = strtol(pt,&qt,10);
	if(qt == pt || NFragments < 0){
	  printf("Line with %s has invalid integer value:\n%s\n",key,mbuf);
	  fflush(stdout);exit(1);
	}
	continue;
      }

      for(int c = 0;  c < colors; c++){
	char keybuf[256];
	sprintf(keybuf,"Nickase Recognition Site %d:",c+1);

	if((pt = strstr(&mbuf[1],keybuf))){
	  pt += strlen(keybuf);
	  while(*pt && isspace(*pt))
	    pt++;
	  if(!*pt){
	    printf("WARNING: %s on line %d of %s : Nickase string not found at end of line\n",keybuf,linecnt, filename);
	    fflush(stdout);
	    break;// ignore error
	  }
	  for(int i = strlen(pt)-1; isspace(pt[i]);)
	    pt[i--] = '\0';

	  /* NEW304 : strip off color from CMAP Nickase string */
	  char *pt2;
	  if((pt2 = strchr(pt,';')) || (pt2 = strchr(pt,',')))
	    *pt2 = '\0';

	  if(DEBUG && fileid >= num_files && refcolors <= 0){
	    printf("ERROR: #Label Channels must precede %s on line %d of %s\n",keybuf,linecnt,filename);
	    fflush(stdout);exit(1);
	  }
	  if(DEBUG) assert(colors_set);

	  for(int i = strlen(pt); --i >= 0;)
	    pt[i] = tolower(pt[i]);

	  if(!strcmp(pt,"unknown"))/* ignore color if it is "unknown" */
	    continue;

	  if(DEBUG) assert(Nickase[c]);
	  if(fileid < num_files || refcolors == origcolors){/* reading query OR reference with same number of Labels colors as query */
	    if(!strcmp(Nickase[c],"unknown")){/* enzyme name was NOT previously specified */
              #pragma omp critical(fileinput)
              {
                if(!strcmp(Nickase[c],"unknown")){/* If thread race caused condition to change, drop down to alternate case below */
                  if(VERB/* HERE HERE >=2 */){
		    printf("input_cmap(%s):Changing global Nickase[%d]=%p(%s) to ",filename,c,Nickase[c],Nickase[c]);
		    fflush(stdout);
                  }
		  // NOTE : cannot free Nickase[c] since other maps may already point to it. This results in a small memory leakage at end of program
		  //	    if(Nickase[c]) free(Nickase[c]);
		  Nickase[c] = strdup(pt);
		  if(VERB/* HERE HERE >=2 */){
                    printf("%p(%s)\n",Nickase[c],Nickase[c]);
		    fflush(stdout);
                  }
                }
              } // #pragma omp critical(fileinput)
	    }

	    if(strcmp(Nickase[c],"unknown")){/* enzyme name was specified in a previous input file */
	      char *p1, *p2;
	      if(((p1 = strchr(Nickase[c],';')) || (p1 = strchr(Nickase[c],','))) != ((p2 = strchr(pt,';')) || (p2 = strchr(pt,','))) ?
		 strncmp(Nickase[c],pt,min(p2 ? p2-pt : strlen(pt), p1 ? p1-Nickase[c] : strlen(Nickase[c]))) : /* both strings must match up to ; or , */
		 strcmp(Nickase[c],pt) /* both strings must match completely */){ 
		printf("ERROR1:%s on line %d of %s : Nickase[c=%d] string \"%s\" does not match value in previous file of \"%s\"\n",
		       keybuf,linecnt,filename, c, pt, Nickase[c]);
		if(VERB>=1+RELEASE)
		  printf("strchr(Nickase[c],';')= %d, strchr(Nickase[c],',')= %d, strncmp(Nickase[c],pt,strlen(pt))= %d\n",
		    strchr(Nickase[c],';') ? 1 : 0, strchr(Nickase[c],',') ? 1 : 0, strncmp(Nickase[c],pt,strlen(pt)));
		fflush(stdout); exit(1);
	      } else if(VERB && strcmp(Nickase[c], pt)){
		printf("%s on line %d of %s : Nickase[c=%d] string \"%s\" is consistent with value in previous file of \"%s\"\n",
		       keybuf,linecnt,filename, c, pt, Nickase[c]);
		fflush(stdout);
	      }
	    }

	  } else if(refcolors == 1 && origcolors == 2){// reading reference with 1 color after reading query with 2 colors
	    if(DEBUG) assert(fileid >= num_files);
	    if(!usecolor){
	      printf("Reading reference with 1 color AND query with 2 colors : Please rerun with -usecolor 1 OR -usecolor 2\n");
	      fflush(stdout);exit(1);
	    }
	    char *p1, *p2;
	    int c = usecolor-1;
	    if(((p1 = strchr(Nickase[c],';')) || (p1 = strchr(Nickase[c],','))) != ((p2 = strchr(pt,';')) || (p2 = strchr(pt,','))) ?
	       strncmp(Nickase[c],pt,min(p2 ? p2-pt : strlen(pt), p1 ? p1-Nickase[c] : strlen(Nickase[c]))) : /* both strings must match up to ; or , */
	       strcmp(Nickase[c],pt) /* both strings must match completely */){ 
	      printf("ERROR:%s on line %d of %s : Nickase string \"%s\" does not match Nickase[%d] value in previous file of \"%s\" with -usecolor %d\n",
		     keybuf,linecnt,filename,pt,usecolor-1,Nickase[usecolor-1],usecolor);
	      fflush(stdout);exit(1);
	    } else if(VERB && strcmp(Nickase[c], pt)){
	      printf("%s on line %d of %s : Nickase string \"%s\" is consistent with value in previous file of Nickase[%d] = \"%s\"\n",
		     keybuf,linecnt,filename, pt, c, Nickase[c]);
	      fflush(stdout);
	    }
	  } else if(refcolors == 2 && origcolors == 1){// reading reference with 2 colors after reading query with 1 color (OR query with 2 colors AND usecolor > 0)
	    if(DEBUG) assert(fileid >= num_files);
	    if(DEBUG) assert(colors == 2);
	    if(!usecolor){
	      printf("Reading reference with 2 colors AND query with 1 color : Please rerun with -usecolor 1 OR -usecolor 2\n");
	      fflush(stdout);exit(1);
	    }
	    if(c != usecolor - 1)
	      continue;/* ignore color of reference that will not be used  and check next color c*/
	    else {
	      char *p1, *p2;
	      int c = 0;
	      if(((p1 = strchr(Nickase[c],';')) || (p1 = strchr(Nickase[c],','))) != ((p2 = strchr(pt,';')) || (p2 = strchr(pt,','))) ?
    	            strncmp(Nickase[c],pt,min(p2 ? p2-pt : strlen(pt), p1 ? p1-Nickase[c] : strlen(Nickase[c]))) : /* both strings must match up to ; or , */
	            strcmp(Nickase[c],pt) /* both strings must match completely */){ 
	        printf("ERROR:%s on line %d of %s : Nickase[%d] string \"%s\" does not match Nickase value in previous file of \"%s\" with -usecolor %d\n",
	           keybuf,linecnt,filename,c,pt,Nickase[0],usecolor);
		fflush(stdout);exit(1);
	      } else if(VERB && strcmp(Nickase[c], pt)){
	        printf("%s on line %d of %s : Nickase string \"%s\" is consistent with value in previous file of Nickase[0]= \"%s\"\n",
		     keybuf,linecnt,filename, pt, Nickase[c]);
	        fflush(stdout);
	      }
	    }
	  }
	}
      }// for (int c = 0; c < colors; c++) 

      continue;/* skip all other comment lines */      
    }

    if(hmap && !HapPresent){
      printf("WARNING: HapSite not found in header line of %s (expected with .hmap suffix)\n",filename);
      fflush(stdout);
    }

    /* done with header comment lines */    
    if(CmapChimQuality < 0){
      #pragma omp critical(fileinput)
      {
        if(CmapChimQuality < 0)
          CmapChimQuality = ChimQuality;
      }
    }

    if(CmapChimQuality >= 0 && CmapChimQuality != ChimQuality){
      if(!ChimQuality_warning){
	printf("WARNING: ChimQuality level = %d in %s, was %d in previous CMAP input %s\n",ChimQuality,filename, CmapChimQuality,lastfilename);
	fflush(stdout);
	ChimQuality_warning = 1;
      }
    }
    
    if(CmapMask < 0){
      #pragma omp critical(fileinput)
      {
        if(CmapMask < 0)
          CmapMask = MaskPresent;
      }
    }
    if(CmapMask >= 0 && CmapMask != MaskPresent){
      if(!CmapMask_warning){
	printf("WARNING: MaskPresent = %d in %s, was %d in previous CMAP input %s\n",MaskPresent, filename, CmapMask,lastfilename);
	fflush(stdout);
	CmapMask_warning = 1;
      }
    }

    if(NColors <= 0){
      /*      printf("WARNING : lines with \"N Colors:\" missing or invalid value on line %d of %s:\n%s\n",
	     linecnt,filename,mbuf);
      printf("Assuming N Colors = 1\n");
      fflush(stdout);*/
      NColors = colors;
    }
    if(NColors > MAXCOLOR){
      printf("N Colors value %d is too large (MAXCOLOR=%d)\n",NColors,MAXCOLOR);
      fflush(stdout);exit(1);
    }
    if(colors != NColors){
      printf("N Colors=%d in file %s does not match parameter value=%d\n",NColors,filename,colors);
      fflush(stdout);exit(1);
    }
    /*    if(NFragments < 0){
      printf("lines with \"Number of Consensus Maps:\" missing or invalid value on line %d of %s:\n%s\n",
	     linecnt,filename,mbuf);
      fflush(stdout);exit(1);
      }*/

    // #pragma omp critical
    lastfilename = filename;

    break;// non-header line parsing will be continued outside of critical section
  }

  if(VERB>=2){
    printf("Finished Header lines for fileid=%d in %s: linecnt=%d,lastline:\n%s\n", fileid, filename,linecnt,lastline);
    fflush(stdout);
  }

  if(merrorexit)
    return 0;

  for(;(lastline != NULL) || ((lastline = fgets(mbuf,LINESIZ,fp)) != NULL); linecnt++) {
    if(VERB>=2){
      printf("Non-Header line %d in %s:\n%s\n",linecnt,filename,mbuf);
      fflush(stdout);
    }

    lastline = NULL;/* force reading of next line in next iteration */

    /* read in one site per line */
    pt = mbuf;
    contig = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || contig < 0){
      printf("Invalid long long value %lld for CMapId on line %d of %s (previous line : contig=%lld,NumSites=%d,SiteID=%d):\n%s\n",contig,linecnt,filename,lastcontig,lastNumSites,lastSiteID,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }

  Lnextmap:    
    pt = qt;

    if(contig != lastcontig){/* start of new map */
      if(VERB>=2){
	printf("start of new contig=%lld,lastcontig=%lld\n",contig,lastcontig);
	fflush(stdout);
      }
      if(lastcontig >= 0){ 
	if(DEBUG){/* check if previous contig was completed */
	  register int cnt = 0;
	  for(register int c = 0; c < colors; c++){
	    assert(pmap->numsite[c] <= lastNumSites);
	    cnt += pmap->numsite[c];
	    if(fabs(pmap->site[c][pmap->numsite[c]+1] - pmap->len) >= 2e-4){
	      printf("Right end of Map %lld with %d labels at %0.4f bp does not match map length = %0.4f bp : line %d of %s\n", 
		     pmap->id,pmap->numsite[c], pmap->site[c][pmap->numsite[c]+1]*1000.0, pmap->len*1000.0, linecnt,filename);
	      fflush(stdout);exit(1);
	    }
	    if(DEBUG) assert(fabs(pmap->site[c][pmap->numsite[c]+1] - pmap->len) < 2e-4);
	    if(pmap->len == 0.0){
	      if(MININTERVAL > 0.00011){
		printf("ERROR: Map with zero length and %d labels on line %d of %s (use -minEnd 0.0001 to allow)\n",pmap->numsite[c],linecnt,filename);
		fflush(stdout);exit(1);
	      }
	      printf("WARNING: Map with zero length and %d labels on line %d of %s\n",pmap->numsite[c],linecnt,filename);
	      fflush(stdout);
	    }
	  }
	  if(cnt != lastNumSites){
	    printf("input_cmap:file=%s,id=%lld:cnt=%d,lastNumSites=%d\n",
		   filename,pmap->id,cnt,lastNumSites);
	    for(register int c = 0; c < colors;c++)
	      printf("c=%d:pmap->numsite[c]=%d\n",c,pmap->numsite[c]);
	    fflush(stdout);
	    assert(cnt == lastNumSites);
	  }
	}
      }

      lastcontig = contig;

      /* make sure the array Map[] is large enough */
      #pragma omp critical(fileinput)
      {
	maxmapalloc(nummaps+1,maxmaps,Map,0,1);// NOTE : uses "#pragma omp critical(MAX_BLOCKS)"

        if(VERB/* HERE >=2 */ && nummaps > 0 && (nummaps % 100000)==0){
          printf("%d maps read...\n",nummaps);
	  fflush(stdout);
	}

	pmap = Map[nummaps];
        pmap->mapid = nummaps++;

	if(!(0 < nummaps && nummaps <= MASK(31))){
          printf("Too many maps in %s up to line %d: Total number of maps (including previous input files) limited to 2 billion (%d)\n",filename, linecnt,MASK(31));
	  fflush(stdout);exit(1);
        }

	pmap->id = SEQID ? nummaps : contig;
      }

      nmapcnt++;

      pmap->allfree();
      pmap->linenum = linecnt;
      pmap->fileid = fileid;
      for(int c = 0; c < colors; c++){
	pmap->Nickase[c] = Nickase[c];
	if(VERB>=3 && Map==refmap){
	  printf("refmap[%d]->Nickase[c=%d] -> %s\n", pmap->mapid, c, pmap->Nickase[c]);
	  fflush(stdout);
	}
	pmap->numsite[c] = 0;
      }
      lastlen = 0.0;
      lastSiteID = 0;
      lastNumSites = 0;
    }

    pmap->startloc = pmap->endloc = 0.0;

    pmap->len = strtod(pt,&qt)*0.001;
    if(pt == qt || !(*qt==0 || isspace(*qt))){
      printf("Invalid float value for ContigLength on line %d of %s:\n%s\n",linecnt,filename,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    pt = qt;
    if(lastlen > 0.0 && fabs(pmap->len - lastlen) > 2e-4){
      printf("ContigLength %0.1f on line %d of %s does not match value %0.1f on prevous line:\n%s\n",pmap->len * 1000.0, linecnt,filename, lastlen * 1000.0, mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    lastlen = pmap->len;

    int NumSites = strtol(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || NumSites < 0 ){
      printf("Invalid int value for NumSites on line %d of %s:\n%s\n",linecnt,filename,mbuf);
      printf("Numsites=%d\nqt=\"%s\"\npt=\"%s\"\n",NumSites,qt,pt);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    pt = qt;
    if(lastNumSites > 0 && NumSites != lastNumSites){
      printf("NumSites value on line %d of %s differs from previous line:\n%s\n",linecnt,filename,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    lastNumSites = NumSites;
    for(int c = 0; c < colors; c++){
      if(!pmap->site[c]){
	if(DEBUG>=2) assert(!pmap->blockmem);
	pmap->site[c] = new FLOAT[NumSites+2];
	pmap->siteSD[c] = new double[NumSites+2];
	pmap->sitecov[c] = new float[NumSites+2];
	pmap->sitecnt[c] = new float[NumSites+2];
	if(ChimQuality){
	  pmap->ChimQuality[c] = new float[NumSites+2];
	  if(ChimQuality >= 2){
	    pmap->ChimNorm[c] = new float[NumSites+2];
	    pmap->SegDupL[c] = new float[NumSites+2];
	    pmap->SegDupR[c] = new float[NumSites+2];
	    pmap->FragileEndL[c] = new float[NumSites+2];
	    pmap->FragileEndR[c] = new float[NumSites+2];
	    pmap->OutlierFrac[c] = new float[NumSites+2];
	    if(ChimQuality >= 3){
	      pmap->FragSd[c] = new float[NumSites+2];
	      pmap->ExpSd[c] = new float[NumSites+2];
	      pmap->FragCov[c] = new float[NumSites+2];
	      pmap->FragChiSq[c] = new float[NumSites+2];
	    }
	  }
	}

	if(MaskPresent)
	  pmap->Mask[c] = new size_t[NumSites+2];

	if(HapPresent){
	  if(!pmap->contig){
	    pmap->contig = new Ccontig;
	    pmap->contig->mapid = -1;
	    pmap->contig->id = pmap->id;
	    pmap->contig->nummaps = -1;/* until 2 Alleles have been generated from Hmap information */
	    pmap->contig->numsite[0] = NumSites;
	  }
	  Ccontig *pcontig = pmap->contig;
	  if(!pcontig->site[c]){
	    pcontig->site[c] = new double[NumSites+2];
	    pcontig->HapSite[c] = new int[NumSites+2];
	    pcontig->HapSite[c][0] = pcontig->HapSite[c][NumSites+1] = 0;// NEW
	    pcontig->HapDelta[c] = new double[NumSites+1];
	  }
	  pcontig->site[c][0] = 0.0;/* left end */
	}

	if((fileid < num_files ? MapTrueSite : RefTrueSite) || Tracksites){
          pmap->truesiteL[c] = new long long[NumSites+2];
          pmap->truesiteR[c] = new long long[NumSites+2];
	}

	/* left end */
	pmap->site[c][0] = 0.0;
	pmap->siteSD[c][0] = 0.0;
	pmap->sitecov[c][0] = 1.0;
	pmap->sitecnt[c][0] = 1.0;
	if(ChimQuality){
	  pmap->ChimQuality[c][0] = END_CHIMQUAL;
	  if(ChimQuality >= 2){
	    pmap->ChimNorm[c][0] = END_CHIMQUAL;
	    pmap->SegDupL[c][0] = END_CHIMQUAL;
	    pmap->SegDupR[c][0] = END_CHIMQUAL;
	    pmap->FragileEndL[c][0] = 0.0;
	    pmap->FragileEndR[c][0] = 0.0;
	    pmap->OutlierFrac[c][0] = 0.0;
	    if(ChimQuality >= 3){
	      pmap->FragSd[c][0] = 0.0;
	      pmap->ExpSd[c][0] = 0.0;
	      pmap->FragCov[c][0] = 0.0;
	      pmap->FragCov[c][0] = 0.0;
	    }
	  }
	}
	if(MaskPresent)
	  pmap->Mask[c][0] = 0;

	/* right end : handled later */
	pmap->site[c][1] = 0.0;/* to catch missing right end */
      }
    }

    int SiteID = strtol(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || SiteID <= 0){
      printf("Invalid +ve int value %d for SiteID on line %d of %s:\n%s\n",SiteID, linecnt,filename,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    if(SiteID != lastSiteID + 1){
      printf("Out of order SiteID=%d on line %d of %s (lastSiteID=%d) :\n%s\n",SiteID, linecnt,filename,lastSiteID,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    lastSiteID = SiteID;
    if(SiteID > NumSites + 1){
      printf("SiteID=%d larger than (NumSites=%d)+1 on line %d of %s:\n%s\n",SiteID,NumSites,linecnt,filename,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    pt = qt;
    
    int LabelChannel = strtol(pt,&qt,10);
    if(pt==qt || !(*qt==0 || isspace(*qt)) || LabelChannel < 0){
      printf("Invalid int value for LabelChannel on line %d of %s:\n%s\n",linecnt,filename,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    if(LabelChannel > colors){
      printf("contig=%lld,len=%0.1f,NumSites=%d,SiteID=%d,LabelChannel=%d\n",contig,pmap->len*1000.0, NumSites, SiteID, LabelChannel);
      printf("LabelChannel=%d too large for -colors %d (or \"Label Channels:\" header of %s) on line %d of %s:\n%s\n",LabelChannel,colors,filename,linecnt,filename,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    if(LabelChannel == 0 && !(SiteID == NumSites + 1)){
      printf("While reading line %d of %s: LabelChannel=%d, SiteID=%d, NumSites=%d (Map end with LabelChannel=0 should have SiteID == NumSites + 1)\n%s\n",
	     linecnt,filename,LabelChannel, SiteID, NumSites, mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      //      assert(SiteID == NumSites + 1);
    }
    if(LabelChannel > 0 && SiteID > NumSites){
      printf("While reading line %d of %s: LabelChannel=%d, SiteID=%d, NumSites=%d (Other than Map end with LabelChannel=0, SiteID must be <= NumSites)\n%s\n",
	     linecnt,filename,LabelChannel, SiteID, NumSites, mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      //      assert(SiteID <= NumSites);
    }
    pt = qt;

    double Position = strtod(pt,&qt);
    if(pt==qt || !(*qt==0 || isspace(*qt)) || Position < 0.0){
      printf("Invalid value for Position on line %d of %s:\n%s\n",linecnt,filename,mbuf);
      printf("Position=%0.2f\nqt=\"%s\"\n",Position,qt);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    pt = qt;
    Position *= 0.001;
    if(Position > lastlen + 1e-6){
      printf("Position larger then length on line %d of %s:(Position=%0.1f,map length=%0.1f)\n%s\n",linecnt,filename,Position*1000.0,lastlen*1000.0,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    Position = min(Position,lastlen);

    if(VERB>=3 && pmap->id == 5){
      printf("input_cmap(%s):Map[%d]:id=%lld:site=%d,Position=%0.3f\n", filename, pmap->mapid, pmap->id, pmap->numsite[0]+1, Position);
      fflush(stdout);
    }

    register int color = LabelChannel - 1;
    if(color >= 0){
      int site = pmap->numsite[color];
      if(Position < pmap->site[color][site]){
	printf("Label Position on line %d of %s is to the left of previous Label (or Map Start):(Position=%0.7f, last Position=%0.7f at LabelChannel=%d,SiteID=%d)\n%s\n",
	       linecnt,filename,Position*1000.0,pmap->site[color][site]*1000.0,LabelChannel,site,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }
    } else if(SiteID != NumSites + 1){
      printf("Decreasing SiteID=%d (previous value=%d),LabelChannel=%d on line %d of %s:\n%s\n",SiteID,NumSites,LabelChannel,linecnt,filename,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }

    double StdDev = strtod(pt,&qt);
    if(pt==qt || !(*qt==0 || isspace(*qt))){
      printf("Invalid float value for StdDev on line %d of %s:\n%s\n",linecnt,filename,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    if(StdDev < 0.0){
      if(VERB && !stddev_warning){
	printf("WARNING:Negative float value for StdDev on line %d of %s (changing to 0):\n%s\n",linecnt,filename,mbuf);
	fflush(stdout);
	stddev_warning = 1;
      }
      StdDev = 0.0;
    }
    if(!isfinite(StdDev)){
      if(VERB && !stddev_warning){
	printf("WARNING:Infinite float value for StdDev on line %d of %s (changing to 1e+30):\n%s\n",linecnt,filename,mbuf);
	fflush(stdout);
	stddev_warning = 1;
      }
      StdDev = 1e+30;
    }
    pt = qt;

    double Coverage = strtod(pt,&qt);
    if(pt==qt || !(*qt==0 || isspace(*qt)) || Coverage < 0.0){
      printf("Invalid float or int for Coverage on line %d of %s:\n%s\n",linecnt,filename,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    pt = qt;

    double Occurrence = strtod(pt,&qt);
    if(pt==qt || !(*qt==0 || isspace(*qt)) || Occurrence < 0.0){
      printf("Invalid float or int value for Occurence on line %d of %s:\n%s\n",linecnt,filename,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
    }
    pt = qt;

    double Quality = -1.0, ChimNorm = -1.0, SegDupL = -1.0, SegDupR = -1.0, FragileEndL = 0.0, FragileEndR = 0.0, OutlierFrac = 0.0, FragSd = 0.0, ExpSd = 0.0, FragCov = 0.0, FragChiSq = 1.0;
    size_t Mask = 0;
    long long truesiteL = -1, truesiteR = -1;

    if(!noQuality || ChimQuality){
      /* check this is end of line */
      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	if(ChimQuality)
	  printf("Missing ChimQuality information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	else
	  printf("Missing Quality information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
      }

      Quality = strtod(pt,&qt);
      if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(Quality)){
	if(ChimQuality)
	  printf("Invalid float value for ChimQuality on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	else
	  printf("Invalid float value for Quality on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }
      if(ChimQuality && Quality < -1.0){
	if(!chimqual_warning){
	  printf("WARNING:Invalid negative value for ChimQuality on line %d of %s (replacing value with -1.0):\n%s\n",linecnt,filename,mbuf);
	  fflush(stdout);
	  chimqual_warning = 1;
	}
	Quality = -1.0;
      }
      if(ChimQuality && Quality > 100.0){
	if(!chimqual_warning){
	  printf("WARNING:Invalid large value for ChimQuality on line %d of %s (replacing value with 100.0):\n%s\n",linecnt,filename,mbuf);
	  fflush(stdout);
	  chimqual_warning = 1;
	}
	Quality = 100.0;
      }
      pt = qt;

      if(ChimQuality >= 2){
	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	if(! *pt){
	  printf("Missing SegDupL information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	SegDupL = strtod(pt,&qt);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(SegDupL)){
	  printf("Invalid float value for SegDupL on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	if(SegDupL < -1.0){
	  if(!chimqual_warning){
	    printf("WARNING:Invalid negative value for SegDupL on line %d of %s (replacing value with -1.0):\n%s\n",linecnt,filename,mbuf);
	    fflush(stdout);
	    chimqual_warning = 1;
	  }
	  SegDupL = -1.0;
	}
	if(SegDupL > 100.0){
	  if(!chimqual_warning){
	    printf("WARNING:Invalid large value for SegDupL on line %d of %s (replacing value with 100.0):\n%s\n",linecnt,filename,mbuf);
	    fflush(stdout);
	    chimqual_warning = 1;
	  }
	  SegDupL = 100.0;
	}

	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	if(! *pt){
	  printf("Missing SegDupR information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	SegDupR = strtod(pt,&qt);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(SegDupR)){
	  printf("Invalid float value for SegDupR on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	if(SegDupR < -1.0){
	  if(!chimqual_warning){
	    printf("WARNING:Invalid negative value for SegDupR on line %d of %s (replacing value with -1.0):\n%s\n",linecnt,filename,mbuf);
	    fflush(stdout);
	    chimqual_warning = 1;
	  }
	  SegDupR = -1.0;
	}
	if(SegDupR > 100.0){
	  if(!chimqual_warning){
	    printf("WARNING:Invalid large value for SegDupR on line %d of %s (replacing value with 100.0):\n%s\n",linecnt,filename,mbuf);
	    fflush(stdout);
	    chimqual_warning = 1;
	  }
	  SegDupR = 100.0;
	}

	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	if(! *pt){
	  printf("Missing FragileEndL information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	FragileEndL = strtod(pt,&qt);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(FragileEndL)){
	  printf("Invalid float value for FragileEndL on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	if(FragileEndL < 0.0){
	  if(!chimqual_warning){
	    printf("WARNING:Invalid negative value for FragileEndL on line %d of %s (replacing value with 0.0):\n%s\n",linecnt,filename,mbuf);
	    fflush(stdout);
	    chimqual_warning = 1;
	  }
	  FragileEndL = 0.0;
	}

	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	if(! *pt){
	  printf("Missing FragileEndR information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	FragileEndR = strtod(pt,&qt);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(FragileEndR)){
	  printf("Invalid float value for FragileEndR on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	if(FragileEndR < 0.0){
	  if(!chimqual_warning){
	    printf("WARNING:Invalid negative value for FragileEndR on line %d of %s (replacing value with 0.0):\n%s\n",linecnt,filename,mbuf);
	    fflush(stdout);
	    chimqual_warning = 1;
	  }
	  FragileEndR = 0.0;
	}

	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	if(! *pt){
	  printf("Missing OutlierFrac information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	OutlierFrac = strtod(pt,&qt);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(OutlierFrac)){
	  printf("Invalid float value for OutlierFrac on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	if(OutlierFrac < 0.0){
	  if(!chimqual_warning){
	    printf("WARNING:Invalid negative value for OutlierFrac on line %d of %s (replacing value with 0.0):\n%s\n",linecnt,filename,mbuf);
	    fflush(stdout);
	    chimqual_warning = 1;
	  }
	  OutlierFrac = 0.0;
	}
	if(OutlierFrac > 100.0){
	  if(!chimqual_warning){
	    printf("WARNING:Invalid large value for OutlierFrac on line %d of %s (replacing value with 100.0):\n%s\n",linecnt,filename,mbuf);
	    fflush(stdout);
	    chimqual_warning = 1;
	  }
	  OutlierFrac = 100.0;
	}

	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	if(! *pt){
	  printf("Missing ChimNorm information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	ChimNorm = strtod(pt,&qt);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(ChimNorm)){
	  printf("Invalid float value for ChimNorm on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	if(ChimNorm < -1.0){
	  if(!chimqual_warning){
	    printf("WARNING:Invalid negative value for ChimNorm on line %d of %s (replacing value with -1.0):\n%s\n",linecnt,filename,mbuf);
	    fflush(stdout);
	    chimqual_warning = 1;
	  }
	  ChimNorm = -1.0;
	}
      }
    }

    if(MaskPresent){
      /* check this is end of line */
      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	printf("Missing Mask information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf), linecnt, filename, mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
      }
	
      Mask = strtoul(pt,&qt,16);
      if(pt==qt || !(*qt==0 || isspace(*qt))){
	printf("Invalid Hex value for Mask at char %ld on line %d of %s:\n%s\nstart at:%s\n",(long)(qt-mbuf), linecnt, filename, mbuf, pt);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
      }
    }

    if(ChimQuality >= 3){/* read in MolSd,ExpSd,MolCov,ChiSq */
      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	printf("Missing MolSd information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
      }
      FragSd = strtod(pt,&qt);
      if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(FragSd)){
	printf("Invalid float value for FragSd on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
      }
      if(FragSd < 0.0){
	if(!chimqual_warning){
	  printf("WARNING:Invalid negative value for FragSd on line %d of %s (replacing value with 0.0):\n%s\n",linecnt,filename,mbuf);
	  fflush(stdout);
	  chimqual_warning = 1;
	}
	FragSd = 0.0;
      }

      ExpSd = 0.0;
      if(ChimQuality >= 4){
	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	if(! *pt){
	  printf("Missing ExpSd information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	ExpSd = strtod(pt,&qt);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(ExpSd)){
	  printf("Invalid float value for ExpSd on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	if(ExpSd < 0.0){
	  if(!chimqual_warning){
	    printf("WARNING:Invalid negative value for ExpSd on line %d of %s (replacing value with 0.0):\n%s\n",linecnt,filename,mbuf);
	    fflush(stdout);
	    chimqual_warning = 1;
	  }
	  ExpSd = 0.0;
	}
      }

      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	printf("Missing MolCov information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
      }
      FragCov = strtod(pt,&qt);
      if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(FragCov)){
	printf("Invalid float value for FragCov on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
      }
      if(FragCov < 0.0){
	if(!chimqual_warning){
	  printf("WARNING:Invalid negative value for FragCov on line %d of %s (replacing value with 0.0):\n%s\n",linecnt,filename,mbuf);
	  fflush(stdout);
	  chimqual_warning = 1;
	}
	FragCov = 0.0;
      }

      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	printf("Missing ChiSq information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
      }
      FragChiSq = strtod(pt,&qt);
      if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(FragChiSq)){
	printf("Invalid float value for FragChiSq on line %d of %s:\n%s\n",linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
      }
      if(FragChiSq < 0.0){
	if(!chimqual_warning){
	  printf("WARNING:Invalid negative value for FragChiSq on line %d of %s (replacing value with 0.0):\n%s\n",linecnt,filename,mbuf);
	  fflush(stdout);
	  chimqual_warning = 1;
	}
	FragChiSq = 0.0;
      }
      if(FragChiSq > 1.0){
	if(!chimqual_warning){
	  printf("WARNING:Invalid large value for FragChiSq on line %d of %s (replacing value with 1.0):\n%s\n",linecnt,filename,mbuf);
	  fflush(stdout);
	  chimqual_warning = 1;
	}
	FragChiSq = 1.0;
      }
    }

    int HapSite = 3; 
    double HapDelta = 0.0;
    if(HapPresent){/* for now only saves HapSite and HapDelta information */
      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	printf("Missing HapSite information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }

      HapSite = strtol(pt,&qt,10);
      if(pt==qt || !(*qt==0 || isspace(*qt)) || HapSite < 0 || HapSite > 3){
	printf("Invalid HapSite value %d at char %ld on line %d in %s:\n%s\n", HapSite, (long)(qt-mbuf), linecnt, filename, mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }

      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	printf("Missing HapDelta information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }
      HapDelta = strtod(pt,&qt);
      if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(HapDelta)){
	printf("Invalid HapDelta value=%0.1f at char %ld on line %d in %s:\n%s\n",HapDelta,(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }

      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	printf("Missing SitePhaseConf information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }
      double SitePhaseConf = strtod(pt,&qt);
      if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(SitePhaseConf)){
	printf("Invalid SitePhaseConf value= %0.2f at char %ld on line %d in %s:\n%s\n",SitePhaseConf,(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }

      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	printf("Missing DeltaPhaseConf information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }
      double DeltaPhaseConf = strtod(pt,&qt);
      if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(DeltaPhaseConf)){
	printf("Invalid DeltaPhaseConf= %0.2f value at char %ld on line %d in %s:\n%s\n",DeltaPhaseConf,(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }

      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	printf("Missing HapSiteScore information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }
      double HapSiteScore = strtod(pt,&qt);
      if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(HapSiteScore)){
	printf("Invalid HapSiteScore value= %0.2f at char %ld on line %d in %s:\n%s\n",HapSiteScore,(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }

      while(*qt && isspace(*qt))
	qt++;
      pt = qt;
      if(! *pt){
	printf("Missing HapDeltaScore information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }
      double HapDeltaScore = strtod(pt,&qt);
      if(pt==qt || !(*qt==0 || isspace(*qt)) || !isfinite(HapDeltaScore)){
	printf("Invalid HapDeltaScore value=%0.2f at char %ld on line %d in %s:\n%s\n",
	       HapDeltaScore,(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }
    }

    if(!HapPresent){
      if(fileid < num_files ? MapTrueSite : RefTrueSite){
	/* check this is end of line */
	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	if(!*pt){
	  printf("Missing Left Label Origin information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;		  
	}

	truesiteL = strtoll(pt,&qt,10);
	if(pt==qt || !(*qt==0 || isspace(*qt))){
	  printf("Invalid Left Label Origin (parsed as %lld) at char %ld on line %d in %s:\n%s\n",truesiteL,(long)(pt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
	}

	/* check this is end of line */
	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;
	if(!*pt){
	  printf("Missing Right Label Origin information at char %ld on line %d in %s:\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;		  
        }

	truesiteR = strtoll(pt,&qt,10);
	if(pt==qt || !(*qt==0 || isspace(*qt))){
	  printf("Invalid Right Label Origin (parsed as %lld) at char %ld on line %d in %s:\n%s\n",truesiteR,(long)(pt-mbuf),linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
	}
      }
    }

    /* check if this is end of line */
    while(*qt && isspace(*qt))
      qt++;
    pt = qt;
    if(color >= 0){/* check for optional SNRgmean,lnSNRsd,SNRcnt */
      if(SNRpresent){/* read in SNRgmean,lnSNRsd,SNRcnt */
	if(!pmap->SNRcnt[color]){
	  for(int c = 0; c < colors; c++){
	    pmap->SNRcnt[c] = new int[NumSites+2];
	    pmap->SNRgmean[c] = new double[NumSites+2];
	    pmap->lnSNRsd[c] = new double[NumSites+2];
	    pmap->SNRdist[c] = new double*[NumSites+2];
	    for(int t = 0; t <= NumSites+1; t++){
	      pmap->SNRcnt[c][t] = 0;
	      pmap->SNRdist[c][t] = NULL;
	    }
	  }
	}
	int site = pmap->numsite[color]+1;

	double SNRgmean = strtod(pt,&qt);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || SNRgmean < 0.0 || !isfinite(SNRgmean)){
	  printf("Invalid float value for SNRgmean on line %d in %s:\n%s\n",linecnt,filename,mbuf);
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}
	pt = qt;

	double lnSNRsd = strtod(pt,&qt);
	if(pt==qt || !(*qt==0 || isspace(*qt)) || lnSNRsd < 0.0 || !isfinite(lnSNRsd)){	
	  if(!float_warning){
	    printf("WARNING:Invalid float value \"%s\" for lnSNRsd on line %d in %s:\n%s\n",pt,linecnt,filename,mbuf);
	    fflush(stdout);
	    float_warning = 1;
	  }
	  lnSNRsd = 0.0;
	  cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	}

	pmap->SNRgmean[color][site] = SNRgmean;
	pmap->lnSNRsd[color][site] = lnSNRsd;
	
	/* check if this is end of line */
	while(*qt && isspace(*qt))
	  qt++;
	pt = qt;

	if(*qt){/* read in SNRcnt and SNR values */
	  int SNRcnt = strtol(pt,&qt,10);
	  if(pt==qt || !(*qt==0 || isspace(*qt)) || SNRcnt < 0){
	    printf("Invalid SNRcnt value=%d at char %ld on line %d in %s starting at: %s\nLine %d:%s\n", SNRcnt, (long)(pt-mbuf), linecnt,filename,pt,linecnt,mbuf);
	    cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	  }

	  pt = qt;

	  pmap->SNRcnt[color][site] = SNRcnt;
	  pmap->SNRdist[color][site] = new double[SNRcnt];
	  if(VERB>=3){
	    printf("Allocating Map[%d](%p)->SNRdist[%d][%d] = %p (SNRcnt=%d)\n",pmap->mapid,pmap,color,site,pmap->SNRdist[color][site],SNRcnt);
	    fflush(stdout);
	  }
	
	  for(int t = 0; t < SNRcnt; t++){
	    double SNRval = strtod(pt,&qt);
	    if(pt==qt || !(*qt==0 || isspace(*qt)) || SNRval < 0.0 || !isfinite(SNRval)){
	      if(*pt==0)
		printf("Missing %d'th of %d SNR values at end of line %d in %s:\n%s\n",t+1, SNRcnt, linecnt,filename,mbuf);
	      else
		printf("Invalid %d'th of %d SNR values on line %d in %s:\n%s\n",t+1, SNRcnt, linecnt,filename,mbuf);
	      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;		      
	    }
	    pmap->SNRdist[color][site][t] = SNRval;
	    pt = qt;
	  }

	  /* sort the SNR values in ascending order */
	  qsort(pmap->SNRdist[color][site],SNRcnt,sizeof(double),(intcmp*)doubleInc);

	  /* check if this is end of line */
	  while(*qt && isspace(*qt))
	    qt++;
	  if(*qt && strncmp(qt,"NA",2)){
	    printf("Invalid char(%c) at position %ld on line %d in %s (hmap=%d):\n%s\n",*qt,(long)(qt-mbuf),linecnt,filename,hmap,mbuf);
	    cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
	  }
	}	
      } else if(pmap->SNRcnt[color]){
	printf("Missing SNR information at char %ld on line %d of %s(present for previous sites):\n%s\n",(long)(qt-mbuf),linecnt,filename,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }
    }

    if(0/* HERE HERE */ && color >= 0 && *qt && strncmp(qt,"NA",2)){
      printf("WARNING:Invalid char(%c) at position %ld on line %d in %s (hmap=%d):\n%s\n",*qt,(long)(qt-mbuf),linecnt,filename,hmap,mbuf);
      cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;	
    }

    if(VERB>=2){
      printf("Line %d in %s:mapid %lld,site=%d/%d:color=%d,Position=%0.6f,StdDev=%0.6e,Coverage=%0.1f,Occurence=%0.1f\n",
	     linecnt,filename,pmap->id,pmap->numsite[color]+1,NumSites,color,Position,StdDev,Coverage,Occurrence);
      fflush(stdout);
    }

    /* now add site information to pmap */
    if(color >= 0){
      if(DEBUG) assert(SiteID <= NumSites);
      register int site = ++pmap->numsite[color];
      pmap->site[color][site] = Position;
      pmap->site[color][site+1] = 0.0;/* to catch missing right end */
      if(VERB>=2){
	printf("id=%lld: color=%d,SiteID=%d,NumSites=%d,Position= %0.6f, site=%d,numsite[color]=%d\n",pmap->id,color,SiteID,NumSites,Position,site,pmap->numsite[color]);
	fflush(stdout);
      }
      pmap->siteSD[color][site] = StdDev * 0.001;
      pmap->sitecov[color][site] = Coverage;
      pmap->sitecnt[color][site] = Occurrence;

      if(ChimQuality){
	pmap->ChimQuality[color][site] = Quality;
	if(ChimQuality >= 2){
	  pmap->ChimNorm[color][site] = ChimNorm;
	  pmap->SegDupL[color][site] = SegDupL;
	  pmap->SegDupR[color][site] = SegDupR;
	  pmap->FragileEndL[color][site] = FragileEndL;
	  pmap->FragileEndR[color][site] = FragileEndR;
	  pmap->OutlierFrac[color][site] = OutlierFrac;
	  if(ChimQuality >= 3){
	    pmap->FragSd[color][site] = FragSd;
	    pmap->ExpSd[color][site] = (ChimQuality >= 4) ? ExpSd : 0.0;
	    pmap->FragCov[color][site] = FragCov;
	    pmap->FragChiSq[color][site] = FragChiSq;
	  }
	}
      }

      if(MaskPresent){
	if(VERB>=2 && Mask){
	  printf("Reading from %s, id=%lld: Mask[%d][%d] = %lx\n",filename,pmap->id,color,site,Mask);
	  fflush(stdout);
	}
	pmap->Mask[color][site] = Mask;
      }
      if(HapPresent){
	Ccontig *pcontig = pmap->contig;
	pcontig->site[color][site] = Position;
	if(DEBUG) assert(1 <= site && site <= pcontig->numsite[color]);
	pcontig->HapSite[color][site] = HapSite;
	pcontig->HapDelta[color][site-1] = HapDelta * 0.001;
	// NOTE : currently the other 4 fields are not read in
	if(VERB>=2){
	  printf("id=%lld:i=%d/%d:site[i]=%0.3f,HapSite[i]=%d,HapDelta[i-1]= %0.3f\n",
		 pmap->id, site, NumSites, Position, HapSite, pcontig->HapDelta[0][site-1]);
	  fflush(stdout);
	}
      }
      if(fileid < num_files ? MapTrueSite : RefTrueSite){
	if(DEBUG) assert(truesiteL > 0 && truesiteR > 0);
	pmap->truesiteL[color][site] = truesiteL;
	pmap->truesiteR[color][site] = truesiteR;
      }
    } else {/* right end */
      if(!(fabs(pmap->len - Position) < 2e-4)){
	printf("Last site at %0.4f on line %d of %s does not match map length = %0.4f:\n%s\n",Position*1000.0,linecnt,filename,pmap->len*1000.0,mbuf);
	cmap_error(contig,lastcontig,fp,linecnt,mbuf,filename,pt,qt); goto Lnextmap;
      }
      if(DEBUG) assert(fabs(pmap->len - Position) < 2e-4);
      for(register int color = 0; color < colors; color++){
	register int site = pmap->numsite[color]+1;
	pmap->site[color][site] = Position;
	pmap->siteSD[color][site] = StdDev * 0.001;
	pmap->sitecov[color][site] = Coverage;
	pmap->sitecnt[color][site] = Occurrence;
	if(ChimQuality){
	  pmap->ChimQuality[color][site] = END_CHIMQUAL;
	  if(ChimQuality >= 2){
	    pmap->ChimNorm[color][site] = END_CHIMQUAL;
	    pmap->SegDupL[color][site] = END_CHIMQUAL;
	    pmap->SegDupR[color][site] = END_CHIMQUAL;
	    pmap->FragileEndL[color][site] = 0.0;
	    pmap->FragileEndR[color][site] = 0.0;
	    pmap->OutlierFrac[color][site] = 0.0;
	    if(ChimQuality >= 3){
	      pmap->FragSd[color][site] = 0.0;
	      pmap->ExpSd[color][site] = 0.0;
	      pmap->FragCov[color][site] = 0.0;
	      pmap->FragChiSq[color][site] = 1.0;
	    }
	  }
	}
	if(MaskPresent){
	  if(DEBUG) assert(pmap->Mask[color]);
	  pmap->Mask[color][pmap->numsite[color]+1] = Mask;
	  if(VERB>=2 && Mask){
	    printf("Reading from %s, id=%lld: right end Mask= %lx\n",filename, pmap->id, Mask);
	    fflush(stdout);
	  }
	}
	if(HapPresent){
	  Ccontig *pcontig = pmap->contig;
	  pcontig->site[color][site] = Position;
	  pcontig->HapSite[color][site] = HapSite;
	  pcontig->HapDelta[color][site-1] = HapDelta * 0.001;
	  // NOTE : currently the other 4 fields are not read in
	  if(VERB>=2){
	    printf("id=%lld:i=%d/%d:site[i]=%0.3f,HapSite[i]=%d,HapDelta[i-1]= %0.3f (Right End)\n",
		   pmap->id, site, NumSites, Position, HapSite, pcontig->HapDelta[0][site-1]);
	    fflush(stdout);
	  }
	}
      }

      if(Tracksites && !(fileid < num_files ? MapTrueSite : RefTrueSite)) { /* generate truesite values to track original index values */
	if(colors == 1)
	  for(int i = 1; i <= pmap->numsite[0]; i++)
	    pmap->truesiteL[0][i] = pmap->truesiteR[0][i] = i;
	else {// first convert 2-color map into an interleaved format (code copied from output_xmap.cpp)
	  int *NN = pmap->numsite;
	  FLOAT **YY = pmap->site;
	  Cinterleaved *Ysite = new Cinterleaved[NN[0]+NN[1]+2]; // array of (c,index), corresponding to YY[c][index] and in ascending order of this value
	  int *Yindex[2];
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
	  if(DEBUG) assert(Ysites == NN[0]+NN[1]);
	  // Now Yindex[c][i] is translation from original index i of color c to interleaved index

	  for(int c = 0; c < colors; c++)
	    for(int i = 1; i <= NN[c]; i++)
	      pmap->truesiteL[c][i] = pmap->truesiteR[c][i] = Yindex[c][i];

	  delete [] Ysite;
	  delete [] Yindex[0];
	}
      }
    }
    continue;/* go to next line (no error) */
  }

  if(DEBUG && lastcontig >= 0){/* check if last contig was completed */
    register int cnt = 0;
    for(register int c = 0; c < colors; c++){
      assert(pmap->numsite[c] <= lastNumSites);
      cnt += pmap->numsite[c];
      if(fabs(pmap->site[c][pmap->numsite[c]+1] - pmap->len) >= 2e-4){
	printf("Right end of Map %lld with %d labels at %0.2f bp does not match map length = %0.2f bp : line %d of %s\n", 
	       pmap->id, pmap->numsite[c], pmap->site[c][pmap->numsite[c]+1]*1000.0, pmap->len*1000.0, linecnt,filename);
	fflush(stdout);exit(1);
      }
      if(DEBUG) assert(fabs(pmap->site[c][pmap->numsite[c]+1] - pmap->len) < 2.0e-4);
      if(pmap->len == 0.0){
	if(MININTERVAL > 0.00011){
	  printf("ERROR: Map with zero length and %d labels on line %d of %s (use -minEnd 0.0001 to allow)\n",pmap->numsite[c],linecnt,filename);
	  fflush(stdout);exit(1);
	}
	printf("WARNING: Map with zero length and %d labels on line %d of %s\n",pmap->numsite[c],linecnt,filename);
	fflush(stdout);
      }
    }
    if(cnt != lastNumSites){
      printf("input_cmap:file=%s,id=%lld:cnt=%d,lastNumSites=%d\n",
	     filename,pmap->id,cnt,lastNumSites);
      for(register int c = 0; c < colors;c++)
	printf("c=%d:pmap->numsite[c]=%d\n",c,pmap->numsite[c]);
      fflush(stdout);
      assert(cnt == lastNumSites);
    }
  }

  if(errorexit){
    if(errorexit > 1)
      printf("%d map errors encountered in %s\n", errorexit,filename);
    else
      printf("1 map error encountered in %s\n",filename);
    fflush(stdout);exit(1);
  }
  (void)fclose(fp);

  if(Tracksites && fileid >= num_files && !RefTrueSite){ /* append current file name to RefOrigin */
    if(!RefOrigin)
      RefOrigin = strdup(filename);
    else {
      size_t len = strlen(RefOrigin) + strlen(filename) + 2;
      char *nRefOrigin = (char *)malloc(len * sizeof(char));
      if(!nRefOrigin){
	printf("malloc(%lu) failed\n",len);
	fflush(stdout);exit(1);
      }
      sprintf(nRefOrigin,"%s,%s",RefOrigin,filename);
      free(RefOrigin);
      RefOrigin = nRefOrigin;
    }
  }

  if(Tracksites && fileid < num_files && !MapTrueSite){ /* append current file name to QueryOrigin */
    if(!QueryOrigin)
      QueryOrigin = strdup(filename);
    else {
      size_t len = strlen(QueryOrigin) + strlen(filename) + 2;
      char *nQueryOrigin = (char *)malloc(len * sizeof(char));
      if(!nQueryOrigin){
	printf("malloc(%lu) failed\n",len);
	fflush(stdout);exit(1);
      }
      sprintf(nQueryOrigin,"%s,%s",QueryOrigin,filename);
      free(QueryOrigin);
      QueryOrigin = nQueryOrigin;
    }
  }

  if(fileid >= num_files && nmapcnt > 1)
    SplitRef = 0;

  if(fileid < num_files && nmapcnt > 1){
  //    #pragma omp critical
    SplitMap = 0;
  }

  if(VERB>=2 && nmapcnt <= 0){
    printf("No maps in %s\n",vfixx_filename[fileid]);
    fflush(stdout);
  }

  if(fileid >= num_files && num_files > 0 && origcolors == 1 && colors == 2){
    if(VERB){
      printf("Converting 2-color Reference (%d contigs) in %s to 1-color Reference based on -usecolor %d\n", nmapcnt,filename, usecolor);
      fflush(stdout);
    }
    assert(usecolor >= 1 && usecolor <= 2);
    if(usecolor == 2)
      for(int i = orignummaps; i < nummaps; i++)
	Map[i]->colorswap(usecolor);
    colors = 1;

    if(VERB>=2 && Map[orignummaps]->numsite[0] >= 10){
      printf("Reference Map[%d] (id=%lld) left most 10 labels:\n",orignummaps, Map[orignummaps]->id);
      for(int i = 0; i < 10; i++)
	printf("Map[%d]->site[0][%d]= %0.6f\n",orignummaps,i,Map[orignummaps]->site[0][i]);
      fflush(stdout);
    }
  }

  if(DEBUG && fileid >= num_files) assert(!multithreaded);

  // NOTE : below call to check_cmap() is NOT re-entrant due to dependency on nummaps : call it for all maps after all multithreaded input_cmap() calls are complete

  if(!multithreaded){
    if(VERB>=2){
      printf("nmapcnt=%d, orignummaps=%d, nummaps=%d : calling check_cmap()\n",nmapcnt,orignummaps,nummaps);
      fflush(stdout);
    }
    if(DEBUG) assert(nmapcnt == nummaps-orignummaps);
    
    check_cmap(fileid, filename,orignummaps,nummaps,maxmaps,Map,HapPresent);

    if(VERB>=2){
      printf("After calling check_cmap(): nmapcnt=%d, orignummaps= %d, nummaps=%d\n",nmapcnt,orignummaps,nummaps);
      fflush(stdout);
    }
    //    if(DEBUG) assert(nmapcnt <= nummaps-orignummaps);
    return nummaps - orignummaps;
  }

  return nmapcnt;
}

int input_vfixx(int fileid,int &nummaps, int &maxmaps, Cmap **&Map, std::set<long long> *ids, bool multithreaded, char *mbuf)
{
  char *filename = vfixx_filename[fileid];

  register char *pt;
  char *qt;
  FILE *fp;
  
  size_t len = strlen(filename);

  /* check if input file is actually a .cmap file, if so use seperate input parser */
  if(len >= strlen(".cmap") && (!strcmp(&filename[len-strlen(".cmap")],".cmap") || !strcmp(&filename[len-strlen(".cmap")],".CMAP"))){
    //    printf("Calling input_cmap():filename=%s\n",filename);
    //    fflush(stdout);
    return input_cmap(fileid, nummaps, maxmaps, Map,0,multithreaded,mbuf);
  }

  if(len >= strlen(".hmap") && !strcmp(&filename[len-strlen(".hmap")],".hmap")){
    //    printf("Calling input_cmap():filename=%s\n",filename);
    //    fflush(stdout);
    return input_cmap(fileid, nummaps, maxmaps, Map,1,multithreaded,mbuf);
  }

  if(len >= strlen(".bnx") && !strcmp(&filename[len-strlen(".bnx")],".bnx")){
    int ret;
    #pragma omp critical(fileinput)
    {
      ret = input_bnx(fileid, nummaps, maxmaps, Map, ids);
    }
    return ret;
  }

  if(!(len >= strlen(".vfixx") && !strcmp(&filename[len-strlen(".vfixx")],".vfixx"))){/* Treat any unrecognized file suffix as BNX input */
    printf("WARNING: Will try to read -i input file %s as BNX format file\n", filename);
    fflush(stdout);

    int ret;
    #pragma omp critical(fileinput)
    {
      ret = input_bnx(fileid, nummaps, maxmaps, Map, ids);
    }
    return ret;
  }

  // input VFIXX file
  Cmap *pmap=0;
  int color;
  int orignummaps = nummaps;

  int NColors= -1,NFragments= -1;
  double BasesPerPixel = -1.0;

  int linecnt=1;

  int fpcnt = 0;
  double lensum = 0.0;

  assert(colors >= 1);

  MapTrueSite = 1;/* .vfixx files always contain TrueSite information */

#pragma omp critical(fileinput) // critical section required since VFIXX file input code is not re-entrant
{

  if((fp = Fopen(filename,"r",NFSdelay))==NULL){
    printf("input_vfixx:Failed to read input file %s (%d'th of %d files)\n",filename, fileid+1, num_files);
    fflush(stdout);exit(1);
  }

  if( &Map == &Gmap ){
    if(maptype == 1 && MapType == -1){
      printf("input_vfixx:fileid=%d, file %s:previous input map type was CMAP which cannot be mixed with .bnx or .vfixx\n",fileid,filename);
      fflush(stdout);exit(1);
    }
    maptype = 0;
  }

  for(int c = 0; c < colors; c++)
    if(!Nickase[c]){
      Nickase[c] = strdup((char*) "unknown");
      if(VERB>=2){
	printf("input_vfixx(%s):setting global Nickase[%d]=%p(%s)\n",filename,c,Nickase[c],Nickase[c]);
	fflush(stdout);
      }
    }

  for(;fgets(buf,LINESIZ,fp) != NULL;linecnt++) {
    int len = strlen(buf);
   
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      fflush(stdout);exit(1);
    }
    
    if(buf[0] == '#'){/* most comment lines are ignored, except "N Colors", "N Fragments","Nickase <color>" and "Bases per Pixel" */
      char *key = (char *)"N Colors:";
      if((pt=strstr(&buf[1],key))){
	pt += strlen(key);
	NColors = strtol(pt,&qt,10);
	if(qt == pt){
	  printf("Line with %s has invalid integer value:\n%s\n",key,buf);
	  fflush(stdout);exit(1);
	}
	if(colors_set){
	  if(NColors != colors){
	    if(DEBUG) assert(fileid > 0);
	    printf("Line %d in file %s with %s has value %d which is different from value %d in previous input file %s:\n%s\n",linecnt,filename,key,NColors,colors,vfixx_filename[fileid-1],buf);
	    printf("All -i input files must have the same number of Label colors\n");
	    fflush(stdout);exit(1);
	  }
	} else {
	  for(int c = colors; c < NColors; c++)
	    if(!Nickase[c]){
	      Nickase[c] = strdup((char*) "unknown");
	      if(VERB>=2){
		printf("input_vfixx(%s):setting global Nickase[%d]=%p(%s)\n",filename,c,Nickase[c],Nickase[c]);
		fflush(stdout);
	      }
	    }

	  colors = NColors;
	  colors_set = 1;
	  if(colors == 1 && usecolor==1)
	    usecolor = 0;
	  if(VERB){
	    printf("Line %d in %s: setting colors=%d, colors_set=%d, usecolor= %d due to Label Channels: %d:\n%s\n",linecnt,filename,colors,colors_set,usecolor,NColors,buf);
	    fflush(stdout);
	  }
	}
	continue;
      }

      key = (char *)"N Fragments:";
      if((pt=strstr(&buf[1],key))){
	pt += strlen(key);
	NFragments = strtol(pt,&qt,10);
	if(qt == pt){
	  printf("Line with %s has invalid integer value:\n%s\n",key,buf);
	  fflush(stdout);exit(1);
	}
	continue;
      }

      key = (char *)"Bases per Pixel:";
      if((pt=strstr(&buf[1],key))){
	pt += strlen(key);
	BasesPerPixel = strtod(pt,&qt);
	if(qt == pt){
	  printf("Line with %s has invalid integer value:\n%s\n",key,buf);
	  fflush(stdout);exit(1);
	}
	if(!origPixelLen){
	  origPixelLen = BasesPerPixel * 0.001;
	  if(VERB/* HERE HERE >=2 */){
	    printf("BasesPerPixel = %0.6f, origPixelLen -> %0.9f\n",BasesPerPixel, origPixelLen);
	    fflush(stdout);
	  }
	  if(0 && origPixelLen != 0.500){
	    double invC = PixelLen/0.500;
	    if(VERB){
	      printf("Adjusted MinLen from %0.3f to %0.3f due to Bases per Pixel=%0.3f in %s\n",
		     MinLen, MinLen*invC, BasesPerPixel, filename);
	      fflush(stdout);
	    }
	    MinLen *= invC;
	    if(MaxLen > 0.0){
	      if(VERB){
		printf("Adjusted MaxLen from %0.3f to %0.3f due to Bases per Pixel=%0.3f in %s\n",
		       MaxLen, MaxLen*invC, BasesPerPixel, filename);
		fflush(stdout);
	      }
	      MaxLen *= invC;
	    }
	  }
	} else if(fabs(origPixelLen - BasesPerPixel*0.001) > 0.001){
	  printf("BasesPerPixel=%0.3f on line %d of %s differs from previous PixelLen=%0.3f\n",
		 BasesPerPixel,linecnt,filename,origPixelLen);
	  fflush(stdout);exit(1);
	}
	continue;
      }
      for(int c = 0;  c < colors; c++){
	char keybuf[256];
	sprintf(keybuf,"Nickase %d:",c+1);
	if((pt = strstr(&buf[1],keybuf))){
	  pt += strlen(keybuf);
	  while(*pt && isspace(*pt))
	    pt++;
	  if(!*pt){
	    printf("WARNING: %s on line %d of %s : Nickase string not found at end of line\n",keybuf,linecnt, filename);
	    fflush(stdout);
	    break;// ignore error
	  }
	  for(int i = strlen(pt)-1; isspace(pt[i]);)
	    pt[i--] = '\0';
	  if(DEBUG) assert(Nickase[c]);
	  if(Nickase[c] && strcmp(Nickase[c],"unknown")){
	    if(strcmp(Nickase[c],pt) && strcmp(pt,"unknown")){
	      printf("WARNING:%s on line %d of %s : Nickase[%d] string \"%s\" does not match value in previous file of \"%s\"\n",
		     keybuf,linecnt,filename, c, pt, Nickase[c]);
	      fflush(stdout);
	    }
	  } else if(strcmp(pt,"unknown")) {
	    if(VERB>=2)
	      printf("input_vfixx(%s):Changing global Nickase[%d]=%p(%s) to ",filename,c,Nickase[c],Nickase[c]);
	    // NOTE : cannot free Nickase[c] since other maps may already point to it. This results in a small memory leakage at end of program
	    //	    if(Nickase[c]) free(Nickase[c]);
	    Nickase[c] = strdup(pt);
	    if(VERB>=2){
	      printf("%p(%s)\n",Nickase[c],Nickase[c]);
	      fflush(stdout);
	    }
	  }
	}
      }
      continue;/* skip all other comment lines */
    }

    /* done with header comment lines */
    if(NColors <= 0){
      printf("lines with \"Label Channels:\" missing or invalid value on line %d of %s:\n%s\n",
	     linecnt,filename,buf);
      fflush(stdout);exit(1);
    }
    if(NColors > MAXCOLOR){
      printf("N Colors value %d is too large (MAXCOLOR=%d)\n",NColors,MAXCOLOR);
      fflush(stdout);exit(1);
    }
    if(colors != NColors){
      printf("N Colors=%d in file %s does not match parameter value=%d\n",NColors,filename,colors);
      exit(0);
    }
      
    if(NFragments < 0){
      printf("lines with \"N Fragments:\" missing or invalid value on line %d of %s:\n%s\n",
	     linecnt,filename,buf);
      fflush(stdout);exit(1);
    }
    if(BasesPerPixel <= 0){
      printf("lines with \"Bases per Pixel:\" missing or invalid value on line %d of %s:\n%s\n",
	     linecnt,filename,buf);
      fflush(stdout);exit(1);
    }

    if(buf[0] == '>'){/* Simulated Data information */
      if((color = buf[1] - '0') < 0){
	printf("Invalid line prefix on line %d of file %s:\n%s\n",linecnt,filename,buf);
	fflush(stdout);exit(1);
      }
      long long startloc,endloc,len,startloc2,endloc2,len2;

      switch(color){
      case 0:
	/* read in next Map into Map[nummaps] */
	/* make sure the array map[] is large enough */
	maxmapalloc(nummaps+1,maxmaps,Map,0,1);

	pmap = Map[nummaps];
	pmap->allfree();
	pmap->mapid = nummaps++;
	pmap->fileid = fileid;
	pmap->linenum = linecnt;

	for(register int c = 0; c < colors; c++)
	  pmap->Nickase[c] = Nickase[c];

	/* read in map id and (if present) reference map location information */
	pt = &buf[2];
	pmap->id = strtoll(pt,&qt,10);
	if(pt == qt){
	  printf("Invalid integer value for Mol ID on line %d of %s:\n%s\n",linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	if(VERB>=2){
	  printf("mapid=%d,id=%lld:linecnt=%d,filename=%s\n",pmap->mapid,pmap->id,linecnt,filename);
	  fflush(stdout);
	}

	if(SEQID)
	  pmap->id = nummaps;

	pt = qt;
	startloc = strtoll(pt,&qt,10);
	if(pt == qt){
	  printf("Invalid value %lld for Start (Bp) on line %d of %s:\n%s\n",startloc, linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	pmap->startloc = (double)(startloc) * 0.001;

	pt = qt;
	endloc = strtoll(pt,&qt,10);
	if(pt == qt){
	  printf("Invalid value %lld for End (Bp) on line %d of %s:\n%s\n",endloc,linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	pmap->endloc = (double)(endloc) * 0.001;

	pt = qt;
	len = strtoll(pt,&qt,10);
	if(pt == qt || len < 0){
	  printf("Invalid value %lld for Length (Bp) on line %d of %s:\n%s\n",len,linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	pmap->len = (double)(len) * 0.001;

	pt = qt;
	startloc2 = strtoll(pt,&qt,10);
	if(pt == qt){
	  printf("Invalid value %lld for Start2 (Bp) on line %d of %s:\n%s\n", startloc2, linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	pmap->startloc2 = (double)(startloc2) * 0.001;

	pt = qt;
	endloc2 = strtoll(pt,&qt,10);
	if(pt == qt){
	  printf("Invalid value %lld for End2 (Bp) on line %d of %s:\n%s\n",endloc2,linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	pmap->endloc2 = (double)(endloc2) * 0.001;

	pt = qt;
	len2 = strtoll(pt,&qt,10);
	if(pt == qt || len2 < 0){
	  printf("Invalid value %lld for Length2 (Bp) on line %d of %s:\n%s\n",len2,linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	pmap->len2 = (double)(len2) * 0.001;

	pt = qt;
	pmap->flipped = strtol(pt,&qt,10);
	if(pt == qt || pmap->flipped < 0 || pmap->flipped > 1){
	  printf("Invalid value %d for Flipped on line %d of %s:\n%s\n",pmap->flipped,linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}

	pt = qt;
	pmap->chimflip = strtol(pt,&qt,10);
	if(pt == qt || pmap->chimflip < 0 || pmap->chimflip > 1){
	  printf("Invalid value %d for Chimflip on line %d of %s:\n%s\n",pmap->chimflip,linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}

	pt = qt;
	for(;*pt;pt++)
	  if(!isspace(*pt))
	     break;
	for(qt = pt;*qt;qt++)
	  if(isspace(*qt))
	    break;
	*qt = '\0';
	if(pmap->name) free(pmap->name);
	pmap->name = strdup(pt);

	pmap->startmatch = pmap->endmatch = 0.0;
	pmap->fpcnt = (pmap->startloc || pmap->endloc) ? 0 : -1;

	lensum += pmap->len + pmap->len2;
	continue;

      default:/* read in Reference and actual map site ids for current color channel */
	assert(pmap != 0);
	int origfpcnt = fpcnt;

	/* scan over input line to count how many reference sites (including false positives) are present */
	int numsite = 0;
	long long val;
	for(pt= &buf[2]; *pt;pt = qt){
	  val = strtoll(pt,&qt,10);
	  if(pt == qt)
	    break;
	  if(val < 0)
	    fpcnt++;
	  numsite++;
	}
	while(*qt && isspace(*qt))
	  qt++;
	if(*qt){
	  printf("Invalid last char(%c) at position %ld on line %d in %s:\n%s\n",*qt,(long)(qt-buf),linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}

	if((numsite % 2)){
	  printf("Found %d numbers on line %d in %s (expecting even number of reference sites):\n%s\n",numsite,linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	numsite /= 2;

	pmap->numsite[color-1] = numsite;
	pmap->truesiteL[color-1] = new long long[numsite+1];
	pmap->truesiteR[color-1] = new long long[numsite+1];
	pmap->site[color-1] = new FLOAT[numsite+2];
	pmap->fpcnt += fpcnt - origfpcnt;/* for now, fpcnt is summed over all colors */

	/* rescan line to read in true site locations */
	numsite = 0;
	for(pt= &buf[2]; *pt;pt = qt){
	  val = strtoll(pt,&qt,10);
	  if(pt == qt)
	    break;

	  numsite++;
	  pmap->truesiteL[color-1][numsite] = val;
	  
	  pt = qt;
	  val = strtoll(pt,&qt,10);
	  if(pt == qt)
	    break;
	  pmap->truesiteR[color-1][numsite] = val;
	}
	if(DEBUG) assert(numsite == pmap->numsite[color-1]);

	/* go to next line and read actual map input */
	linecnt++;
	if(fgets(buf,LINESIZ,fp)==NULL || (len = strlen(buf)) >= LINESIZ-1 || buf[len-1] != '\n'){
	  printf("Map line missing, too long or not terminated in input file %s on line %d\n",filename,linecnt);
	  fflush(stdout);exit(1);
	}
	pt = buf;
	double x = strtod(pt,&qt);
	if(pt == qt){
	  printf("Missing size value at start of line %d in %s:\n%s\n",linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	pmap->site[color-1][numsite+1] = x * origPixelLen;
	if(color > 1 && fabs(pmap->site[color-1][numsite+1] - pmap->site[color-2][pmap->numsite[color-2]+1]) >= 1e-6){
	  printf("Length=%0.3f for map on line %d of %s differs from previous length of %0.3f:\n%s\n",x,linecnt,filename,pmap->site[color-2][pmap->numsite[color-2]+1]/origPixelLen,buf);
	  fflush(stdout);exit(1);
	}

	pmap->site[color-1][0] = 0.0;
	for(numsite=0,pt = qt;*pt; pt = qt){
	  x = strtod(pt,&qt);
	  if(pt == qt)
	    break;
	  pmap->site[color-1][++numsite] = x * origPixelLen;
	}
	while(*qt && isspace(*qt))
	  qt++;
	if(*qt){
	  printf("Invalid last character (%c) at position %ld on line %d in %s:\n%s\n",
		 *qt,(long)(qt-buf),linecnt,filename,buf);
	  fflush(stdout);exit(1);
	}
	if(numsite != pmap->numsite[color-1]){
	  printf("%d sites on line %d in %s does not match %d sites on previous line\n%s\n",
		 numsite,linecnt,filename,pmap->numsite[color-1],buf);
	  fflush(stdout);exit(1);
	}

	continue;
      }
      assert(0);
      fflush(stdout);exit(1);
    }
  }

  (void)fclose(fp);

  if(VERB>=2){
    printf("fpcnt=%d,lensum=%0.1f kb : Observed FP=%0.6f\n",
	   fpcnt,lensum,(fpcnt/lensum)*100.0);
    fflush(stdout);
  }

  /* filter out maps smaller than MinLen or larger than MaxLen or with fewer sites than MinSites or greater than MaxSites */
  register int i,j;
  for(i=j=orignummaps; i < nummaps; i++){
    pmap = Map[i];
    if(DEBUG>=2) assert(pmap->mapid == i);
    if(pmap->numsite[0] < MinSites || (MaxSites && pmap->numsite[0] > MaxSites) || pmap->site[0][pmap->numsite[0]+1] < MinLen || (MaxLen > 0.0 && pmap->site[0][pmap->numsite[0]+1] > MaxLen)){
      if(VERB>=2){
	printf("Deleting Map[%d]:id=%lld,numsites=%d(MinSites=%d,MaxSites=%d),len=%0.3f(MinLen=%0.3f,MaxLen=%0.3f),PixelLen=%0.6f,startloc=%0.3f,endloc=%0.3f,len=%0.3f\n",
	       i,pmap->id,pmap->numsite[0],MinSites,MaxSites,pmap->site[0][pmap->numsite[0]+1],MinLen,MaxLen,origPixelLen,pmap->startloc,pmap->endloc,pmap->len);
	fflush(stdout);
      }
      if(DEBUG)
	pmap->id = -pmap->id;
      continue;
    }
    if(j < i){
      Cmap *tmp = Map[j];
      pmap->mapid = j;
      Map[j] = pmap;
      tmp->mapid = i;
      Map[i] = tmp;
    }
    j++;
  }
  if(VERB && j != nummaps){
    printf("Reduced number of input maps in %s from %d to %d due to MinLen=%0.3f,MaxLen=%0.3f MinSites=%d,MaxSites=%d\n",
	   filename, nummaps-orignummaps,j-orignummaps,MinLen,MaxLen,MinSites,MaxSites);
    fflush(stdout);
  }
  nummaps = j;

  /* check that all maps are monotonic and extend ends so they are at least ENDEXTEND */
  for(register int i=orignummaps; i < nummaps; i++){
    pmap = Map[i];
    for(register int c=0; c < colors; c++){
      register FLOAT *X = pmap->site[c];
      register int N = pmap->numsite[c];
      assert(X[0] == 0.0);
      for(register int k=1; k < N;k++)
	if(X[k+1] < X[k]){
	  printf("map id %lld, color %d:site%d = %0.3f, site%d = %0.3f (out of order)\n",
		  pmap->id, c+1, k, X[k]/origPixelLen, k+1, X[k+1]/origPixelLen);
	  fflush(stdout);exit(1);
	}
      if(X[1] < 0.0){
	printf("map id %lld, color %d: site1 = %0.3f (leftmost site is located below 0.0)\n", 
		pmap->id, c+1, X[1]/origPixelLen);
	fflush(stdout);exit(1);
      }
      if(X[N+1] < X[N]){
	printf("map id %lld, color %d: site%d = %0.3f, len = %0.3f (rightmost site is to the right of the end)\n",
		pmap->id, c+1, N, X[N]/origPixelLen, X[N+1]/origPixelLen);
	fflush(stdout);exit(1);
      }
      if(X[1] < ENDEXTEND){
	register double delta = ENDEXTEND-X[1];
	for(register int k=1; k <= N+1; k++)
	  X[k] += delta;
      }
      if(X[N+1] < X[N] + ENDEXTEND)
	X[N+1] = X[N] + ENDEXTEND;
    }
    if(colors >= 2){/* equalize length for all colors */
      FLOAT len = 0.0;
      for(int c = 0; c < colors; c++){
	int M = pmap->numsite[c];
	len = max(len, pmap->site[c][M+1]);
      }
      for(int c = 0; c < colors; c++){
	int M = pmap->numsite[c];
	pmap->site[c][M+1] = len;
      }
      pmap->len = len;
    } else
      pmap->len = pmap->site[0][pmap->numsite[0]+1];
  }
  if(fileid >= num_files){
    if(DEBUG) assert( &Map == &refmap);
    for(register int i=orignummaps; i < nummaps; i++){
      pmap = Map[i];
      /* Initialize refname with filename */
      pmap->refname = strdup(filename);

      /* initialize and copy original site locations (before resolution based condensation) */
      for(register int c=0; c < colors; c++){
	pmap->orignumsite[c] = pmap->numsite[c];
	if(pmap->origsite[c]) delete [] pmap->origsite[c];
	pmap->origsite[c] = new double[pmap->numsite[c]+2];
	if(pmap->remap[c]) delete [] pmap->remap[c];
	pmap->remap[c] = new int[pmap->numsite[c]+1];
	if(pmap->NickCnt[c]) delete [] pmap->NickCnt[c];
	pmap->NickCnt[c] = new int[pmap->numsite[c]+1];
	for(register int J=0; J <= pmap->numsite[c]+1; J++)
	  pmap->origsite[c][J] = pmap->site[c][J];
	for(register int J=1; J <= pmap->numsite[c]; J++){
	  pmap->remap[c][J] = J;
	  pmap->NickCnt[c][J] = 1;
	}
      }
    }
  } else {/* NOT a reference map */
    if(DEBUG && (refmap || Map)) assert(&Map != &refmap);
    if(mres >= 0.0){/* mresSD = 0 : merge map sites at distance mres or closer (if mresSD == 0)
		       mresSD > 0 : merge map sites at distance X with probability 0.5(1-erf((X-mres)/(mresSD*sqrt(2)))) */
      mapcorrect(Map, orignummaps, nummaps, maxmaps, fileid, 0);
    }
  }

  /* check all maps to make sure all colors are present */
  for(register int i=orignummaps; i < nummaps;i++){
    pmap = Map[i];
    for(register int c=0; c < colors; c++)
      if(!pmap->site[c]){
	printf("Color %d missing from Map[%d] with id=%lld in file %s\n",c+1,i-orignummaps,pmap->id,filename);
	fflush(stdout);exit(1);
      }
  }
} // pragma omp critical

  return nummaps-orignummaps;
}

#define MAXREPEATS 3

/* check for duplicate ids in map inputs : this function is called after all input files have been read */
void input_vfixx_duplicates()
{
  int idchanged = 0;/* flag if any molecule id was changed */

  if(VERB>=2){
    printf("input_vfixx_duplicates() entered:time=%0.6f(wall=%0.6f)\n",mtime(),wtime());
    for(int m = 0; m < nummaps; m++){
      Cmap *pmap = Gmap[m];
      if(pmap->id == 1LL && pmap->stitchLocation[0]){
	int N = pmap->NumStitch[0];
	double *StitchLoc = pmap->stitchLocation[0];
	printf("id=%lld,NumStitch[0]=%d:\n",pmap->id,N);
	for(int J = 0; J < N; J++)
	  printf(" StitchLoc[%d]=%0.4f\n",J,StitchLoc[J]);
	fflush(stdout);
      }
    }
    fflush(stdout);
  }

  int numthreads = 1;
#ifdef _OPENMP
  numthreads = omp_get_max_threads();
  numthreads = min(numthreads,MaxThreads);
  numthreads = max(1,min(numthreads, nummaps/1024));
#endif

  if(BNXVersion >= 1 && (RunIndexWarning || RunDataWarning || RunDataListLen <= 0)){/* handle older BNX 1.0 files without RunIndex : set UniqueScanId  & ScanId to ScanNumber-1 and leave id unchanged */
    for(int i = 0; i < nummaps; i++){
      Cmap *pmap = Gmap[i];
      pmap->UniqueScanId = pmap->ScanId = pmap->ScanNumber - 1;
      pmap->RunId = 0;
    }
  }

  if(BNXVersion >= 1 && !(RunIndexWarning || RunDataWarning || RunDataListLen <= 0)){/* compute RunId,ScanId and UniqueScanId for each molecule and compute unique id for each map */

    /* assemble RunDataSort[] with ChipId, FlowCell, date, AM_PM, time fields along with index into RunDataList[] */
    delete [] RunDataSort;
    RunDataSort = new CRunDataSort[RunDataListLen];

    for(int RunIndex = 0; RunIndex < RunDataListLen; RunIndex++){
      CRunDataSort *r = &RunDataSort[RunIndex];
      r->index = RunIndex;
      char *buf = strdup(RunDataList[RunIndex]);
      char *p = buf;

      /* remove spaces from end of string */
      int len = strlen(p);
      while(len > 0 && isspace(p[len-1]))
	p[--len] = 0;

      /* extract tab seperated fields from p */

      /* skip across 3 tabs to get to start of Time field */
      int tabcnt;
      for(tabcnt = 0; *p && tabcnt < 3; p++)
	if(*p == '\t')
	  tabcnt++;
      while(*p && isspace(*p))
	p++;
      if(tabcnt < 3 || *p == 0){
	printf("Failed to find start of 3rd tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      p[-1] = '\0';/* terminate 2nd Field */
      r->F12 = buf;/* Fields 1 and 2 combined */

      /* locate end of 3rd tab seperated field */
      char *q = p;
      while(*q && *q != '\t')
	q++;
      if(!*q || *q != '\t'){
	printf("Failed to find end of 3rd tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      *q++ = '\0';/* end of 3rd Field */

      r->date = p;
      while(*p && !isspace(*p))
	p++;
      if(!*p){
	printf("Failed to find end of 1st space seperated subfield (date) of 3rd tab seperated field (Date/Time) in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      *p++ = '\0';/* end of Field 3.1 (date) */

      while(*p && isspace(*p))
	p++;
      if(!*p){
	printf("Failed to find start of 2nd space seperated subfield (time) of 3rd tab seperated field (Date/Time) in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }

      r->time = p;
      while(*p && !isspace(*p))
	p++;
      if(!*p){
	printf("Failed to find end of 2nd space seperated subfield (time) of 3rd tab seperated field (Date/Time) in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      *p++ = '\0';/* end of Field 3.2 (time) */

      while(*p && isspace(*p))
	p++;
      if(!*p){
	printf("Failed to find start of 3rd space seperated subfield (AM/PM) of 3rd tab seperated field (Date/Time) in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      r->AM_PM = p;

      if(!strcmp(r->AM_PM,"AM") && r->time[0]=='1' && r->time[1]=='2'){/* change 12:xx:xx AM as 00:xx:xx AM, so we can use lexicographical sorting */
	r->time[0] = '0';
	r->time[1] = '0';
      }

      /* end of 3rd subfield is already null terminated at q[-1] */

      r->F456 = p = q;/* start of 4th tab seperated field */

      /* move ahead to start of 7th tab seperated field */
      for(tabcnt = 4; *p && tabcnt < 7; p++)
	if(*p == '\t')
	  tabcnt++;
      while(*p && isspace(*p))
	p++;
      if(tabcnt < 7 || *p == 0){
	printf("Failed to find start of 7th tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      p[-1] = '\0';/* terminate combined Fields 4,5 and 6 */

      r->NumberOfScans = strtol(p,&q,10);
      if(p==q || !(*q==0 || isspace(*q)) || r->NumberOfScans < 0){
	printf("Invalid NumberOfScans %s in %d'th Run Data line in -ii input files:\n%s\n",
	       p, RunIndex+1,RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }

      /* move ahead to start of 8th tab seperated field */
      for(; *p && tabcnt < 8; p++)
	if(*p == '\t')
	  tabcnt++;
      while(*p && isspace(*p))
	p++;
      if(tabcnt < 8 || *p == 0){
	printf("Failed to find start of 8th tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      r->F8 = p;

      /* locate end of 8th tab seperated field */
      q = p;
      while(*q && *q != '\t')
	q++;
      if(!*q || *q != '\t'){
	printf("Failed to find end of 8th tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      *q++ = '\0';

      char *ChipId = p;
      char *pt = strrchr(ChipId,','), *qt;
      if(pt==NULL){
	printf("Invalid ChipId=%s without comma seperated fields in %d'th Run Data line in headers of -i input files:\n%s\n",
	       ChipId,RunIndex+1,RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      r->ChipId = strtoll(&pt[1],&qt,10);
      if(pt==qt || !(*qt==0 || isspace(*qt))){
	printf("Invalid ChipId=%s without valid last numerical field (serial number) in %d'th Run Data line in headers of -i input files:\n%s\n",
	       ChipId,RunIndex+1,RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }

      p = q;/* start of 9th tab seperated field */
      if(!*p){
	printf("Failed to find start of 9th tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      r->FlowCell = strtol(p,&q,10);
      if(p==q || !(*q==0 || isspace(*q)) || r->FlowCell < 0){
	printf("Invalid FlowCell %s in %d'th Run Data line in -i input files:\n%s\n",
	       p,RunIndex+1,RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }

      /* locate end of 9th tab seperated field */
      while(*q && *q != '\t')
	q++;
      if(*q && *q != '\t'){
	printf("Failed to find end of 9th tab seperated field in %d. Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	printf("Start of field at char %ld:%s\n",p-buf,p);
	printf("End of field at char %ld:%s\n",q-buf,q);
	fflush(stdout);exit(1);
      }
      if(*q)
	*q++ = '\0';

      r->F10 = p = q;/* start of 10th tab seperated field */
      if(!*p){/* No 10th field : set default values */
	if(BNXMinorVersion >= 3){
	  printf("Failed to find 10th tab seperated field (SNRFilterType) in %d'th Run Data line in headers of -i input files:\n%s\n",
		 RunIndex+1,RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	r->F10 = NULL;
	r->MinMoleculeLength = 0.0;
	for(int c = 0; c < colors; c++)
	  r->MinLabelSNR[c] = 0.0;
	continue;
      }

      /* skip across 10th tab seperated field */
      while(*q && *q != '\t')
	q++;
      if(*q && *q != '\t'){
	printf("Failed to find end of 10th tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      if(*q)
	*q++ = '\0';
      
      p = q;/* start of 11th tab seperated field */

      if(!*p){/* No 11th field : set default values */
	if(BNXMinorVersion >= 3){
	  printf("Failed to find 11th tab seperated field (MinMoleculeLength) in %d'th Run Data line in headers of -i input files:\n%s\n",
		 RunIndex+1,RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	r->MinMoleculeLength = 0.0;
	for(int c = 0; c < colors; c++)
	  r->MinLabelSNR[c] = 0.0;
	continue;
      }
      r->MinMoleculeLength = strtod(p,&q);
      if(p==q || !(*q==0 || isspace(*q)) || r->MinMoleculeLength < 0.0){
	printf("WARNING:Invalid MinMoleculeLength %s in %d'th Run Data line in -i input files (must be non-negative number):\n%s\n",
	       p,RunIndex+1,RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }

      /* locate end of 11th tab seperated field */
      while(*q && *q != '\t')
	q++;
      if(*q && *q != '\t'){
	printf("Failed to find end of 11th tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
	       RunIndex+1, RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      if(*q)
	*q++ = '\0';
      
      p = q;/* start of 12'th tab seperated field */
      if(!*p){/* No 12'th field */
	if(BNXMinorVersion >= 3){
	  printf("Failed to find 12th tab seperated field (MinLabelSNR) in %d'th Run Data line in headers of -i input files:\n%s\n",
		 RunIndex+1,RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	for(int c = 0; c < colors; c++)
	  r->MinLabelSNR[c] = 0.0;
	continue;
      }

      int c = 0;
      for(; c < (BNXMinorVersion >= 3 ? colors : 1); c++){
	if(c && !*p){
	  printf("Failed to find %d'th tab seperated field (MinLabelSNR for color %d) in %d'th Run Data  line in headers of -i input files:\n%s\n",
		 c+12,c+1,RunIndex+1,RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	r->MinLabelSNR[c] = strtod(p,&q);
	if(p == q || !(*q==0 || isspace(*q))){
	  printf("Invalid MinLabelSNR %s in 12th tab seperated field of %d'th Run Data line in -i input files:\n%s\n",
		 p, RunIndex+1, RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	if(r->MinLabelSNR[c] < 0.0){
	  printf("WARNING:Invalid MinLabelSNR %s in 12th tab seperated field of %d'th Run Data line in -i input files (must be non-negative number):\n%s\n",
		 p, RunIndex+1, RunDataList[RunIndex]);
	  fflush(stdout);
	  r->MinLabelSNR[c] = 0.0;
	}

	/* locate end of c+12th tab seperated field */
	while(*q && *q != '\t')
	  q++;
	if(*q && *q != '\t'){
	  printf("Failed to find end of %d'th tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
		 c+12,RunIndex+1, RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	if(*q)
	  *q++ = '\0';
	
	p = q;
      }
      int col = c+12;
      for(; c < MAXCOLOR; c++)
	r->MinLabelSNR[c] = r->MinLabelSNR[0];
      
      if(BNXRunHeaderChim > 0){/* parse FeatureMinSlope, FeatureMaxSlope */
	if(!*p){
	  printf("Failed to find FeatureMinSlope value in %d'th tab seperated field in %d'th Run Data line in header of -i input files:\n%s\n",
		 col,RunIndex+1,RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	r->FeatureMinSlope = strtod(p,&q);
	if(p == q || !(*q==0 || isspace(*q))){
	  printf("Invalid FeatureMinSlope value %s in %d'th tab seperated field of %d'th Run Data line in -i input files:\n%s\n",
		 p, col, RunIndex+1, RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	while(*q && *q != '\t')
	  q++;
	if(*q && *q != '\t'){
	  printf("Failed to find end of %d'th tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
		 col,RunIndex+1, RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	if(*q)
	  *q++ = '\0';
	p = q;
	col++;

	if(!*p){
	  printf("Failed to find FeatureMaxSlope value in %d'th tab seperated field in %d'th Run Data line in header of -i input files:\n%s\n",
		 col,RunIndex+1,RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	r->FeatureMaxSlope = strtod(p,&q);
	if(p == q || !(*q==0 || isspace(*q))){
	  printf("Invalid FeatureMaxSlope value %s in %d'th tab seperated field of %d'th Run Data line in -i input files:\n%s\n",
		 p, col, RunIndex+1, RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	while(*q && *q != '\t')
	  q++;
	if(*q && *q != '\t'){
	  printf("Failed to find end of %d'th tab seperated field in %d'th Run Data line in headers of -i input files:\n%s\n",
		 col,RunIndex+1, RunDataList[RunIndex]);
	  fflush(stdout);exit(1);
	}
	if(*q)
	  *q++ = '\0';
	p = q;
	col++;
      }

      if(num_files==1){/* check if existing RunId field is correct */
	while(*p && isspace(*p))
	  p++;
	if(*p){
	  int RunId = strtol(p,&q,10);
	  if(p==q || !(*q==0 || isspace(*q)) || RunId != RunIndex+1){
	    printf("WARNING:Invalid RunId %s in %d'th tab seperated field of %d'th Run Data Line in -i input file %s:\n%s\n    RunId should be %d (BNXRunHeaderChim=%d)\n",
		   p, colors+12, RunIndex+1, vfixx_filename[0], RunDataList[RunIndex],RunIndex+1,BNXRunHeaderChim);
	    fflush(stdout);
	  }
	}
      }
    }

    /* regenerate "Run Data" string with correct RunId (last field) */
    for(int RunIndex = 0; RunIndex < RunDataListLen; RunIndex++){
      CRunDataSort *r = &RunDataSort[RunIndex];

      char nbuf[BUFSIZ];
      size_t n;
      /*      if(colors==2)
	n = snprintf(nbuf,BUFSIZ,"%s\t%s %s %s\t%s\t%d\t%s\t%d\t%s\t%0.2f\t%0.6f\t%d",
		r->F12,r->date,r->time,r->AM_PM,r->F456,r->NumberOfScans,r->F8,r->FlowCell,r->F10 ? r->F10 : "static", r->MinMoleculeLength,r->MinLabelSNR[0],r->index+1);
		else */
      if(BNXRunHeaderChim <= 0){
	if(BNXMinorVersion >= 3 && colors >= 2)
	  n = snprintf(nbuf,BUFSIZ,"%s\t%s %s %s\t%s\t%d\t%s\t%d\t%s\t%0.2f\t%0.6f\t%0.6f\t%d",
		       r->F12,r->date,r->time,r->AM_PM,r->F456,r->NumberOfScans,r->F8,r->FlowCell,r->F10 ? r->F10 : "static", r->MinMoleculeLength,r->MinLabelSNR[0],r->MinLabelSNR[1], r->index+1);
	else
	  n = snprintf(nbuf,BUFSIZ,"%s\t%s %s %s\t%s\t%d\t%s\t%d\t%s\t%0.2f\t%0.6f\t%d",
		       r->F12,r->date,r->time,r->AM_PM,r->F456,r->NumberOfScans,r->F8,r->FlowCell,r->F10 ? r->F10 : "static", r->MinMoleculeLength,r->MinLabelSNR[0],r->index+1);
      } else {
	if(BNXMinorVersion >= 3 && colors >= 2)
	  n = snprintf(nbuf,BUFSIZ,"%s\t%s %s %s\t%s\t%d\t%s\t%d\t%s\t%0.2f\t%0.6f\t%0.6f\t%0.5e\t%0.5e\t%d",
		       r->F12,r->date,r->time,r->AM_PM,r->F456,r->NumberOfScans,r->F8,r->FlowCell,r->F10 ? r->F10 : "static", r->MinMoleculeLength,r->MinLabelSNR[0],r->MinLabelSNR[1], 
		       r->FeatureMinSlope, r->FeatureMaxSlope, r->index+1);
	else
	  n = snprintf(nbuf,BUFSIZ,"%s\t%s %s %s\t%s\t%d\t%s\t%d\t%s\t%0.2f\t%0.6f\t%0.5e\t%0.5e\t%d",
		       r->F12,r->date,r->time,r->AM_PM,r->F456,r->NumberOfScans,r->F8,r->FlowCell,r->F10 ? r->F10 : "static", r->MinMoleculeLength,r->MinLabelSNR[0],
		       r->FeatureMinSlope, r->FeatureMaxSlope, r->index+1);
      }
      if(n > BUFSIZ){
	printf("ERROR: %d'th Run Data line in headers of -i input files is larger than %d:\n%s\n",RunIndex+1,BUFSIZ,RunDataList[RunIndex]);
	fflush(stdout);exit(1);
      }
      free(RunDataList[RunIndex]);
      RunDataList[RunIndex] = strdup(nbuf);
    }

    if(RunDataListLen == 1){/* check if global BNXminlen is different from Run Data MinLen : if so reset global value to 0.0 (redundant). Do the same for minSNR[] */
      if(fabs(RunDataSort[0].MinMoleculeLength - BNXminlen) < 0.001)
	BNXminlen = 0.0;
      for(int c = 0; c < colors; c++)
	if(fabs(RunDataSort[0].MinLabelSNR[c] - BNXminSNR[c]) < 0.01)
	  BNXminSNR[c] = 0.0;
    }
    
    /* apply command line values -minlen -minSNR to modify global BNXminlen,BNXminSNR[] */
    BNXminlen = max(BNXminlen,BNXminlen);
    for(int c = 0; c < colors; c++)
      BNXminSNR[c] = max(BNXminSNR[c],minSNR[c]);

    /* sort RunDataSort[0..RunDataListLen-1] in order of ChipId, FlowCell, date, AM_PM, time (as txt except for numerical FlowCell) */
    qsort(RunDataSort, RunDataListLen, sizeof(CRunDataSort), (intcmp *)ChipCellTimeInc);
    
    /* compute RunId, ScanIdStart, UniqueScanIdStart fields */
    int RunId = 0, ScanIdStart = 0, UniqueScanIdStart = 0;
    for(int sindex = 0; sindex < RunDataListLen; sindex++){
      CRunDataSort *r = &RunDataSort[sindex];
      if(!sindex || r[0].ChipId != r[-1].ChipId || r[0].FlowCell != r[-1].FlowCell)
	RunId = ScanIdStart = 0;
      else {
	RunId++;
	ScanIdStart += r[-1].NumberOfScans;
      }
      r->RunId = RunId;
      r->ScanIdStart = ScanIdStart;
      r->UniqueScanIdStart = UniqueScanIdStart;

      /* Restore 00:xx:xx AM to 12:xx:xx AM */
      if(!strcmp(r->AM_PM,"AM") && r->time[0] == '0' && r->time[1] == '0'){
	r->time[0] = '1';
	r->time[1] = '2';
      }

      if(VERB>=2){
	printf("sindex=%d:RunIndex=%d RunId=%d, ScanIdStart=%d, UniqueScanIdStart=%d, NumberOfScans=%d: ChipId=%lld,FlowCell=%d,date=%s,AM_PM=%s,time=%s\n", sindex, r->index, RunId, ScanIdStart, UniqueScanIdStart, r->NumberOfScans, r->ChipId,r->FlowCell,r->date,r->AM_PM,r->time);
	fflush(stdout);
      }

      UniqueScanIdStart += r->NumberOfScans;
    }
    UniqueScans = UniqueScanIdStart;
    if(VERB){
      printf("Input maps contain %d seperate runs and %d unique scans:time=%0.6f(wall=%0.6f)\n",RunDataListLen,UniqueScanIdStart, mtime(), wtime());
      fflush(stdout);
    }

    /* create RunDataIdmap[RunIndex = 0..RunDataListLen-1] as index into RunDataSort[] */
    if(RunDataIdMax < RunDataListLen){
      delete [] RunDataIdmap;
      RunDataIdMax = RunDataListLen;
      RunDataIdmap = new int[RunDataIdMax];
    }
    for(int sindex = 0; sindex < RunDataListLen; sindex++)
      RunDataIdmap[RunDataSort[sindex].index] = sindex;

    /* now compute RunId,ScanId,UniqueScanId for each molecule */
    int maxScanId = 0;
    
    #pragma omp parallel for num_threads(numthreads) schedule(static,1024) if(numthreads > 1)
    for(int i = 0; i < nummaps; i++){
      Cmap *pmap = Gmap[i];
      int RunIndex = pmap->RunIndex;
      int sindex = RunDataIdmap[RunIndex];
      CRunDataSort *r = &RunDataSort[sindex];
      if(DEBUG/* HERE >=2 */) assert(r->index == RunIndex);
      
      if(pmap->ScanNumber > r->NumberOfScans){
	#pragma omp critical
	{
	  printf("Error: Invalid ScanNumber=%d for map id=%lld with RunId=%d (NumberOfScans=%d) in %s\n", pmap->ScanNumber, pmap->id, pmap->RunIndex+1, r->NumberOfScans, vfixx_filename[pmap->fileid]);
	  fflush(stdout);exit(1);
	}
      }
      if(pmap->ScanNumber <= 0){
	if(ScanCorrection){
	  #pragma omp critical
	  {
	    printf("Error: Invalid ScanNumber=%d for map id=%lld with RunId=%d (NumberOfScans=%d) in %s\n", pmap->ScanNumber, pmap->id, pmap->RunIndex+1, r->NumberOfScans, vfixx_filename[pmap->fileid]);
	    fflush(stdout);exit(1);
	  }
	}
	if(!scannumber_warning){
          #pragma omp critical
	  {
	    printf("WARNING: Invalid ScanNumber=%d for map id=%lld with RunId=%d (NumberOfScans=%d) in %s\n", pmap->ScanNumber, pmap->id, pmap->RunIndex+1, r->NumberOfScans, vfixx_filename[pmap->fileid]);
	    printf("    Ignoring since -ScanScaling is not being used\n");
	    fflush(stdout);
	    scannumber_warning = 1;
	  }
	}
	pmap->ScanNumber = 1;// assume this is an old PSFdetect artifact, otherwise we just merged the first 2 scans together
      }

      if(DEBUG>=3 && !(pmap->ScanNumber == 1 && r->RunId == r->UniqueScanIdStart && r->RunId == r->ScanIdStart && r->index == pmap->RunIndex)){
        printf("map id=%lld: pmap->ScanNumber= %d, RunId = %d -> %d, ScanId = %d -> %d, UniqueScanId = %d -> %d, r->index=%d, pmap->RunIndex= %d\n",
           pmap->id,pmap->ScanNumber, pmap->RunId, r->RunId, pmap->ScanId, r->ScanIdStart + pmap->ScanNumber - 1, pmap->UniqueScanId, r->UniqueScanIdStart + pmap->ScanNumber - 1,r->index,pmap->RunIndex);
	fflush(stdout);
      }

      pmap->RunId = r->RunId;
      pmap->ScanId = r->ScanIdStart + pmap->ScanNumber - 1;
      pmap->UniqueScanId = r->UniqueScanIdStart + pmap->ScanNumber - 1;
      maxScanId = max(maxScanId,pmap->ScanId);

      if(bnxBigid){ /* Use 19 digit multi field id number */
	long long newid = ((((r->ChipId % CHIPID_MAX) * FLOWCELL_WIDTH) + pmap->FlowCell) * SCANNUMBER_WIDTH + pmap->ScanId + 1) * MOLID_WIDTH + pmap->OriginalMoleculeId;
	if(newid != pmap->id)
	  idchanged = 1;
	pmap->id = newid;
	if(DEBUG) assert(pmap->id >= 0);
	/* NOTE : occasionally this number may not be unique, eg if some of the fields overflowed : in that case all ids will be sequentialized and RunData information discarded */
      }
      if((VERB>=2 && pmap->id == 21) || (DEBUG && !(pmap->UniqueScanId < UniqueScans))){
	#pragma omp critical
	{
	  printf("pmap->id=%lld: RunIndex=%d, rindex=%d, RunId=%d, ScanId=%d(r->ScanIdStart=%d), UniqueScanId=%d/%d(r->UniqueScanIdStart=%d), ScanNumber=%d (r->NumberOfScans=%d)\n",
		 pmap->id, RunIndex, sindex, pmap->RunId, pmap->ScanId, r->ScanIdStart, pmap->UniqueScanId, UniqueScans, r->UniqueScanIdStart, pmap->ScanNumber, r->NumberOfScans);
	  fflush(stdout);
	  if(DEBUG) assert(pmap->UniqueScanId < UniqueScans);
	}
      }
    }
    if(bnxBigid && maxScanId >= 1000){
      printf("WARNING: maximum value of ScanId = %d : id field for ScanId will overflow\n",maxScanId);
      fflush(stdout);
    }
  }

  if(VERB>=2){
    printf("Before checking for duplicate id's : nummaps=%d,idchanged=%d:time=%0.6f(wall=%0.6f)\n",nummaps,idchanged,mtime(),wtime());
    fflush(stdout);
  }

  if(!SEQID && nummaps > 0){/* check for duplicate or 0 id values */
    int checkcount = 0;
    size_t renumbered = 0;
    int jmax = 0;
    /* set up translation from Gmap[i]->id to Gmap[i]->mapid */
    Cid2mapid *id2mapid = new Cid2mapid[nummaps];

    int idmap_output =  strstr(output_prefix,"/dev/null") ? 0 : 1;
    char filename[PATH_MAX];
    FILE *fp = NULL;
    sprintf(filename,"%s.idmap",output_prefix);

    while(1) {
      int problem = 0;
    
      //      #pragma omp parallel for num_threads(numthreads) schedule(static,1024) if(numthreads > 1)
      for(int i = 0; i < nummaps; i++){
	Cid2mapid *p = &id2mapid[i];
	Cmap *pmap = Gmap[i];
	if(DEBUG>=2) assert(pmap != NULL);
	p->id = pmap->id;
	if(p->id <= 0){
	  problem = 1;
	  if(VERB){
	    printf("WARNING:zero/negative id=%lld numbers in input for %d'th map : replacing with sequential ids starting at 1\n", p->id, i+1);
	    fflush(stdout);
	  }
	  break;
	  //	  continue;
	}
	if(DEBUG/* HERE >=2 */ && !(pmap->mapid == i)){
	  printf("i=%d/%d:Gmap[i]->mapid=%d,id=%lld\n",i,nummaps,pmap->mapid,pmap->id);
	  fflush(stdout);
	  assert(pmap->mapid == i);
	}
	p->mapid = pmap->mapid;
      }

      if(!problem){
	qsort(id2mapid,nummaps,sizeof(Cid2mapid),(intcmp *)idinc);

	//  #pragma omp parallel for num_threads(numthreads) schedule(static,1024) if(numthreads > 1)
	for(int i = 1; i < nummaps; i++){
	  Cid2mapid *p = &id2mapid[i];
	  Cid2mapid *q = &id2mapid[i-1];

	  if(p->id == q->id && ((PairSplit||Cfirst) ? ((p->mapid < Cfirst && q->mapid < Cfirst) || (p->mapid >= Cfirst && q->mapid >= Cfirst)) : 1 /*(CmapMerge || PairMerge || HashGen || Gmap[p->mapid]->fileid == Gmap[q->mapid]->fileid)*/)){
	    //            #pragma omp critical
	    {
	      if(VERB && !problem){
		Cmap *pmap = Gmap[p->mapid];
		Cmap *qmap = Gmap[q->mapid];
		if(PairSplit||Cfirst)
		  printf("WARNING:duplicate id=%lld for %d'th and %d'th maps (-first %d) : ",
			 q->id, p->mapid, id2mapid[i-1].mapid, Cfirst);
		else {
		  if(BNXVersion==1){
		    printf("WARNING:duplicate id=%lld for %d'th and %d'th maps:\n", q->id, q->mapid + 1, p->mapid + 1); 
		    printf("line=%d file=%s (id=%lld,MoleculeId=%d,RunIndex=%d,RunId=%d,ScanId=%d,ChipId=%s,FlowCell=%d,ScanNumber=%d,len=%0.5f)\n",
			   pmap->linenum, vfixx_filename[pmap->fileid], pmap->id, pmap->OriginalMoleculeId, pmap->RunIndex, pmap->RunId, pmap->ScanId, pmap->ChipId, pmap->FlowCell, pmap->ScanNumber, pmap->site[0][pmap->numsite[0]+1]);
		    if(RunDataList)
		      printf("%s\n", RunDataList[min(RunDataListLen-1,pmap->RunIndex)]);
		    printf("line=%d file=%s (id=%lld,MoleculeId=%d,RunIndex=%d,RunId=%d,ScanId=%d,ChipId=%s,FlowCell=%d,ScanNumber=%d,len=%0.5f)\n",
			   qmap->linenum, vfixx_filename[qmap->fileid], qmap->id, qmap->OriginalMoleculeId, qmap->RunIndex, qmap->RunId, qmap->ScanId, qmap->ChipId, qmap->FlowCell, qmap->ScanNumber, qmap->site[0][qmap->numsite[0]+1]);
		    if(RunDataList)
		      printf("%s\n", RunDataList[min(RunDataListLen-1,qmap->RunIndex)]);
		  } else
		    printf("WARNING:duplicate id=%lld for %d'th and %d'th maps:\n\tline=%d file=%s (id=%lld, len=%0.3f)\n\tline=%d file=%s (id=%lld, len=%0.3f)\n",
			   q->id, q->mapid + 1, p->mapid + 1, 
			   pmap->linenum, vfixx_filename[pmap->fileid], pmap->id, pmap->site[0][pmap->numsite[0]+1],
			   qmap->linenum, vfixx_filename[qmap->fileid], qmap->id, qmap->site[0][qmap->numsite[0]+1]);
		}
		//	      if(BNXVersion==1)
		//		printf("Incrementing ScanNumber of duplicate ids by multiples of 100\n"); 
		//	      else
		printf("replacing all ids with sequential ids starting at 1\n"); 
		fflush(stdout);
	      }
	      if(BNXVersion==1 && bnxBigid) 
		RunDataListLen = 0;/* disable output of Run Data information so subsequent ids are not recomputed */
	    }
	    problem = 1;
	    break;
	  }
	}
      }

      if(problem){/* replace id numbers by sequential numbers */
	idchanged = 1;
	checkcount++;// If BNXVersion keep track of how many times we tried resolving (limit to MAXREPEATS using Scan field increments, fall back to sequential numbering)
	if(PairMerge){
	  printf("ERROR: Using -pairmerge with duplicate id numbers in input CMAP(s) will change ids\n");
	  fflush(stdout);
	  exit(1);
	}

	if(idmap_output && fp==NULL && checkFile(filename))
	  idmap_output = 0;
	if(idmap_output && fp==NULL){
	  if((fp = fopen(filename,"w"))==NULL){
	    printf("Cannot write file %s\n",filename);
	    fflush(stdout);exit(1);
	  }
	  printversion(fp);

	  fprintf(fp,"# newID oldID numsites sourceLine sourceFile (NOTE: there may be multiple lines with the same oldID)\n");
	}

	idrenumbered = 1;
	if(0 && BNXVersion == 1 && checkcount <= MAXREPEATS){/* renumber identical ones by adding 100 to the Scan field (assuming repeated Runs from same chipid)  */
	  if(idmap_output)
	    fprintf(fp,"# repeat=%d, renumbered=%lu\n", checkcount, renumbered);
	  for(int i = 0; i < nummaps;i++){
	    long long id = id2mapid[i].id;
	    int j = 1;
	    for(;j < 1000; j++){
	      long long id2 = id2mapid[i+j].id;
	      if(id2 != id)
		break;
	      int mapid = id2mapid[i+j].mapid;
	      long long newid = id2 + (100 * MOLID_WIDTH) * j;
	      Cmap *pmap = Gmap[mapid];
	      if(idmap_output)
		fprintf(fp, "%lld %lld %d %d %s\n", newid, pmap->id, pmap->numsite[0], pmap->linenum, vfixx_filename[pmap->fileid]);
	      /*	      if(newid == 8400128611401000078LL){
		printf("Gmap[%d]->id = %lld -> %lld (j=%d)\n",mapid,pmap->id, newid, j);
		fflush(stdout);
		}*/
	      pmap->id = newid;
	      pmap->ScanNumber += 100 * j;
	      jmax = max(j,jmax);
	      renumbered++;
	    }
	    i += j-1;
	  }
	  continue;/* recheck to make sure none of the new ids collide (can happen if input bnx files were already renumbered) */
	} else {
	  for(int i = 0; i < nummaps;i++){
	    long long newid = i+1;
	    if(idmap_output)
	      fprintf(fp, "%lld %lld %d %d %s\n", newid, Gmap[i]->id, Gmap[i]->numsite[0], Gmap[i]->linenum, vfixx_filename[Gmap[i]->fileid]);
	    Gmap[i]->id = newid;
	  }
	  renumbered = nummaps;
	}
      }
      if(VERB && renumbered){
	if(idmap_output)
	  printf("%lu/%d ID's have been renumbered and a mapping from new to old ids is located in %s\n",renumbered,nummaps,filename);
	else
	  printf("%lu/%d ID's have been renumbered (mapping file %s suppressed)\n",renumbered,nummaps,filename);
	fflush(stdout);
      }
      if(0 && VERB && BNXVersion==1 && checkcount){
	if(checkcount <= MAXREPEATS)
	  printf("NOTE: Only the ScanNumber has been incremented by multiples of 100 : Highest multiple=%d, repeats=%d\n", jmax, checkcount);
	else
	  printf("WARNING: Failed to resolve ids after %d repeats (jmax=%d): replaced all ids with sequential ids starting at 1\n",checkcount-1, jmax);
      }
      break;
    }
    if(idmap_output && fp != NULL)
      (void)fclose(fp);
    delete [] id2mapid;
  } /* !SEQID */

  if(VERB>=2){
    printf("Finished checking for duplicate map ids:time=%0.6f(wall=%0.6f)\n",mtime(),wtime());
    fflush(stdout);
  }

  /* check all maps to locate the largest and smallest id values */
  minid = 1LL << 48;
  maxid = -1;

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {
    long long myminid = 1LL << 48;
    long long mymaxid = -1;
    #pragma omp for schedule(static,512) 
    for(int i= 0; i < nummaps;i++){
      Cmap *pmap = Gmap[i];
      if(pmap->id > mymaxid)
	mymaxid = pmap->id;
      if(pmap->id < myminid)
	myminid = pmap->id;
    }
    
    #pragma omp critical
    {
      minid = min(minid,myminid);
      maxid = max(maxid,mymaxid);
    }
  }

  if(idchanged && !CmapMerge && BNXVersion >= 1){/* need to save <prefix>_renumbered.bnx */
    char buf[PATH_MAX];
    sprintf(buf,"%s_renumbered",output_prefix);
    output_bnx(buf,Gmap,0,nummaps-1, 0, 0, 0, -1);
  }
  if(VERB>=2){
    printf("minid=%lld,maxid=%lld:time=%0.6f(wall=%0.6f)\n",minid,maxid,mtime(),wtime());
    for(int m = 0; m < nummaps; m++){
      Cmap *pmap = Gmap[m];
      if(pmap->id == 1LL && pmap->stitchLocation[0]){
	int N = pmap->NumStitch[0];
	double *StitchLoc = pmap->stitchLocation[0];
	printf("After input_vfixx_duplicates: id=%lld,NumStitch[0]=%d:\n",pmap->id,N);
	for(int J = 0; J < N; J++)
	  printf(" StitchLoc[%d]=%0.4f\n",J,StitchLoc[J]);
	fflush(stdout);
      }
    }
    fflush(stdout);
  }
  
  if(LEAKDEBUG){
    delete [] RunDataIdmap; RunDataIdmap = NULL;
    RunDataIdMax = 0;
  }
}
