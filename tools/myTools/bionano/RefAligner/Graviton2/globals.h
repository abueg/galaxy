#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <set>
#include <map>

#define CALIGN_SMALL 1 /* WAS31 0 */ // replace Calign by CalignA to save memory in Assembler

#ifdef WIN32
#define USE_PFLOAT 1
#define USE_RFLOAT 1
#define PATH_MAX MAX_PATH

#include <io.h>
#include <windows.h>
#include <direct.h>
#define sleep Sleep
#define inline __inline
#define setbuffer(f, buf, s) setvbuf(f, buf, _IOFBF, s)

#ifdef CLANG
#define __extern_always_inline extern always inline
#endif

/* just do regular malloc() under windows, so free() works : this will reduce performance but should work */
static inline int posix_memalign(void** ptr, size_t alignment, size_t siz)
{
  *ptr = malloc(siz);
  /* 	*ptr = (_aligned_malloc(siz, alignment));*/
  return (*ptr==NULL) ? -1 : 0;
}

//Warning these functions are not random enough
#define random rand
#define srandom srand
#define getcwd _getcwd
#if defined(_MSC_VER)
#define strtoll _strtoi64
#endif
#else
#include <unistd.h>
#endif


static inline void linux_tweak(const char *file, const char *setting)
{
#ifndef WIN32
FILE *f;

f=fopen(file, "w");
if(f==NULL) {
	perror(file);
	return;
	}
fprintf(f, "%s", setting);
fclose(f);
#endif
}

//#ifdef _OPENMP
//#include <valgrind/drd.h>
//#define _GLIBCXX_SYNCHRONIZATION_HAPPENS_BEFORE(addr) ANNOTATE_HAPPENS_BEFORE(addr)
//#define _GLIBCXX_SYNCHRONIZATION_HAPPENS_AFTER(addr) ANNOTATE_HAPPENS_AFTER(addr)
//#undef _GLIBCXX_EXTERN_TEMPLATE
//#define _GLIBCXX_EXTERN_TEMPLATE -1
//#endif

#ifdef __ARM_ARCH
#if __ARM_NEON == 1
#include <arm_neon.h>
#endif
#define _mm_prefetch(PT, HINT) // empty, to support ARM

#else
#include <immintrin.h>
#endif


#include <assert.h>
#include <errno.h>
//#include <math.h>
#include <cmath>
#include <iostream>
#include <exception>

#include <string.h>
#include <stdio.h>
#include "ident.h"
#include "constants.h"

extern const char *SVN_ID;
static Ident globals_h_id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/globals.h 11531 2020-08-22 22:23:24Z tanantharaman $");

using namespace std;

#ifdef REPLACE_WITH_AMDLIBM
#include <amdlibm.h>
#endif

#ifdef REPLACE_WITH_INTELM
//#include <mathimf.h>
#endif

#ifdef WIN32
#include <float.h>
#define isfinite(x) (_finite(x))

#define isnan(x) _isnan(x)
#define isinf(x) (!isfinite(x) && !isnan(x))
    #ifndef nan
        static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
        #define nan() (*(const double *) __nan)
    #endif

#else // ! WIN32

#include <float.h>

#if 1 // #ifndef __INTEL_COMPILER // gcc will ignore builtin or "math.h" or <cmath> based isfinite(),isnan() or isinf() if -Ofast or -ffast-math is used, so define our own version

static inline bool isfinite(float x)
{
  union { float f; int I; } u = {x};
  int exponentBits = u.I & 0x7F800000;

  return exponentBits != 0x7F800000;
}

static inline bool isfinite(double x)
{
  union { double f; long long LL; } u = {x};
  long long exponentBits = u.LL & 0x7FF0000000000000;

  return exponentBits != 0x7FF0000000000000;
}

static inline bool isfinite(long double x)
{
  return isfinite((double)x);
}

static inline bool IsNan(float x)
{
  union { float f; int I; } u = {x};
  int exponentBits = u.I & 0x7F800000;
  int mantissaBits = u.I & 0x7FFFFF;

  return (exponentBits == 0x7F800000) && mantissaBits;
}

static inline bool IsNan(double x)
{
  union { double f; long long LL; } u = {x};
  long long exponentBits = u.LL & 0x7FF0000000000000;
  long long mantissaBits = u.LL & 0x000FFFFFFFFFFFFF;

  return (exponentBits == 0x7FF0000000000000) && mantissaBits;
}

static inline bool IsNan(long double x)
{
  return IsNan((double)x);
}

#else // Intel compiler

static inline bool IsNan(float x)
{
  return isnan(x);
}

static inline bool IsNan(double x)
{
  return isnan(x);
}

#endif // Intel Compiler

static inline bool IsInf(float x)
{
  return !isfinite(x) && !IsNan(x);
}

static inline bool IsInf(double x)
{
  return !isfinite(x) && !IsNan(x);
}

static inline bool IsInf(long double x)
{
  return IsInf((double)x);
}

#endif // Linux

/* set default values for Makefile options USE_FLOAT, USE_RFLOAT & USE_SSE */
#ifndef USE_PFLOAT
#define USE_PFLOAT 0
#endif

#ifndef USE_RFLOAT
#define USE_RFLOAT 0
#endif

#ifndef USE_MFLOAT
#define USE_MFLOAT 0
#endif

#ifndef USE_EPOW
#define USE_EPOW 0
#endif

#ifndef USE_VEXP
#define USE_VEXP 0
#endif

#ifndef USE_SEXP
#define USE_SEXP 0
#endif

#ifndef USE_SSE
#ifndef __ARM_ARCH
#define USE_SSE 1     // NOTE : USE_SSE==0 means NEVER use intrinsics. USE_SSE==1 means use intrinsics for SSE2,SSE3 or AVX (if USE_AVX==1) or AVX512 (if USE_AVX512) IF available based on compiler defines of __SSE2__ __SSE3__ and __AVX__ and (__AVX512F etc) respectively.
#else
#define USE_SSE 0
#endif
#endif

#ifndef USE_AVX
#define USE_AVX USE_SSE  // NOTE (currently applies to pairalign only): USE_AVX==0 (only intended for debuggin) means No use of intrinsics for AVX. USE_AVX==1 means use intrinsincs for AVX IF available (AND USE_SSE==1) 
#endif

#ifndef USE_MIC
#define USE_MIC 0 // NOTE: USE_MIC==0 (only intended for debugging) means No use of intrinsics for MIC. USE_MIC==1 means use intrinsincs for AVX IF available (AND USE_SSE==1) 
#endif

#if !defined(USE_AVX512)
#if __AVX512F__ && __AVX512BW__ && __AVX512CD__ && __AVX512DQ__
#define USE_AVX512 USE_SSE
#else
#define USE_AVX512 0
#endif
#endif

#ifndef __SSE__
#undef __SSE2___
#endif

#ifndef __SSE2__
#undef __SSE3__
#endif

#ifndef __SSE3__
#undef USE_SSE
#define USE_SSE 0       // Normally USE_SSE is defined as 1 IFF  __SSE__ and __SSE2__ and __SSE3__ are defined, but can also be forced to 0 in Makefile for debugging.
#endif

#ifndef __AVX__
#undef USE_AVX
#define USE_AVX 0      // Normally USE_AVX is defined as 1 IFF __AVX__ is defined, but can also be forced to 0 in Makefile for debugging.
#endif


#if USE_AVX
#define _mm256_set_m128(hi, lo) \
  _mm256_insertf128_ps(_mm256_castps128_ps256(lo), (hi), 0x1)
#endif

/* Note: FLOAT or PFLOAT are not used if precision is required, PFLOAT is only used in pairalign() */

#if USE_PFLOAT==1
#define PFLOAT float
#define PMAX fmaxf
#define PMIN fminf
#else
#define PFLOAT FLOAT
#define PMAX fmax
#define PMIN fmin
#endif

#undef PMAX
#define PMAX max
#undef PMIN
#define PMIN min
 

#ifdef WIN32
#define MM_HINT int
#else
#ifdef __INTEL_COMPILER
#define MM_HINT int
#define INT32_MIN ((1)<<31)
#else
#ifdef CLANG
#define MM_HINT const int
#else
#define MM_HINT _mm_hint
#endif // not CLANG
#endif // not __INTEL_COMPILER
#endif // not WIN32

/* Moved to memcpy_wrap.c so it can be compiled seperately with gcc 4.6.2 */
#ifndef WIN32
#if (!USE_MIC)
#ifndef __INTEL_COMPILER
#ifndef USE_STATIC
// Magic that makes us compatible to earlier versions of glibc
__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");
#endif
#endif
#endif
#endif
/**/

extern long long OMP_stack_size;

#define INDEL (InDel && (spots_filename[0] || PairSplit)) /* output .indel file of insertions and deletions (same format as .smap version 0.1) */


#if (USE_MIC || USE_SSE || USE_AVX)

#ifndef CLANG

static inline void prefetch_array(void *p, size_t len, MM_HINT hint)
{
  p = (void *)(((size_t)p) & ~0x3f);
  len += 64;
  for(size_t cnt=0;cnt<len;cnt+=64) {
    _mm_prefetch(&(((char *)p)[cnt]), hint);
  }
}

#else // CLANG : does not handle MM_HINT

#define prefetch_array(p,  len, hint)		\
{						\
  p = (typeof(p))(((size_t)p) & ~0x3f);		\
  for(size_t cnt= 0; cnt < len + 64; cnt += 64) {		\
    _mm_prefetch(&(((char *)p)[cnt]), hint);	\
  }						\
}

#endif // CLANG

#else // !(USE_MIC || USE_SSE || USE_AVX)
#define prefetch_array(p,  len, hint) {}
#endif

/* Compute preferred stride for array of M entries of size 4 bytes assuming 8-way associative cache with 64 byte cache lines : (M rounded to nearest multiple of 32) + 16 */
#define COMPUTE_STRIDE4(M)	((((M) + 0x1f) & ~0x1f) + 16)

/* Compute preferred stride for array of M entries of size 4 bytes assuming 8-way associative cache with 64 byte cache lines : (M rounded to nearest multiple of 16) */
// TRY #define COMPUTE_STRIDE4(M)	(((M) + 0xf) & ~0xf)

#if USE_MIC

#include "scif_io.h"
/* remove definition of complex I that RefAligner does not need */
#undef I
#undef complex
#define fopen sic_fopen
#define unlink sic_unlink

#define EXTRACTF0(a)  (_mm512_mask_reduce_gmin_ps(_mm512_int2mask(1), a))

static inline float maxf_mic(float a, float b)
{
float res  __attribute__((aligned(64)));
_mm512_mask_store_ps(&res, _mm512_int2mask(1), _mm512_gmax_ps(_mm512_set1_ps(a), _mm512_set1_ps(b)));
return res;
}

static inline float minf_mic(float a, float b)
{
float res  __attribute__((aligned(64)));
_mm512_mask_store_ps(&res, _mm512_int2mask(1), _mm512_gmin_ps(_mm512_set1_ps(a), _mm512_set1_ps(b)));
return res;
}

static inline float reciprocal_mic(float a)
{
float res  __attribute__((aligned(64)));
_mm512_mask_store_ps(&res, _mm512_int2mask(1), _mm512_rcp23_ps(_mm512_set1_ps(a)));
return res;
//return EXTRACTF0(_mm512_rcp23_ps(_mm512_set1_ps(a)));
}

static inline double maxd_mic(double a, double b)
{
double res  __attribute__((aligned(64)));
_mm512_mask_store_pd(&res, _mm512_int2mask(1), _mm512_gmax_pd(_mm512_set1_pd(a), _mm512_set1_pd(b)));
return res;
}

static inline double mind_mic(double a, double b)
{
double res  __attribute__((aligned(64)));
_mm512_mask_store_pd(&res, _mm512_int2mask(1), _mm512_gmin_pd(_mm512_set1_pd(a), _mm512_set1_pd(b)));
return res;
}

// #ifdef USE_PFLOAT
// #undef PMAX
// #define PMAX maxf_mic
// #undef PMIN
// #define PMIN minf_mic
// #else
// #undef PMAX
// #define PMAX maxd_mic
// #undef PMIN
// #define PMIN mind_mic
// #endif


#endif // USE_MIC

#if USE_RFLOAT==1
#define RFLOAT float
#else
#define RFLOAT FLOAT
#endif

#if USE_MFLOAT==1
#define MFLOAT float
#else
#define MFLOAT FLOAT
#endif

#define TFLOAT RFLOAT // double // floating point precision used in refalign2.cpp

#define FLOAT double

#define PZERO ((PFLOAT)0.0)
#define RZERO ((RFLOAT)0.0)
#define ZERO ((FLOAT)0.0)
#define ZERO128  _mm_setzero_ps()

extern void dumpmemmap();
extern void getmem(long long &VmSize, long long &VmRSS, long long &VmSwap, long long *VmHWM = NULL, long long *VmPeak = NULL);
extern void getswap(long long &MemSize, long long &SwapSize, long long *AvailableMem = NULL);

static inline void POSIX_MEMALIGN_F(void **ptr, size_t alignment, size_t siz)
{
  long long VmSize, VmRSS, VmSwap, VmHWM;

  if(alignment == 0){
    *ptr = (void *)malloc(siz);
    if(! *ptr){
      getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
      printf("\n VmSize= %0.4f GB, VmRSS= %0.4f GB, VmSwap= %0.4f GB, VmHWM= %0.4f Gb\n",
	     VmSize*1e-9,VmRSS*1e-9,VmSwap*1e-9,VmHWM*1e-9);
      printf("ERROR:POSIX_MEMALIGN: malloc of %lu bytes failed\n",siz);
      fflush(stdout);
      assert(0);
    }
  } else {
    int err = posix_memalign(ptr, alignment, siz);  
    if(err){
#pragma omp critical
      {
	getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
	printf("\n VmSize= %0.4f GB, VmRSS= %0.4f GB, VmSwap= %0.4f GB, VmHWM= %0.4f Gb\n",
	  VmSize*1e-9,VmRSS*1e-9,VmSwap*1e-9,VmHWM*1e-9);
	printf("ERROR:posix_memalign() of %lu bytes (alignment= %lu) failed:err=%d:%s\n", siz,alignment,err,strerror(err));
	fflush(stdout);
	//      dumpmemmap();
	assert(0);
      }
    }
  }
}

#define POSIX_MEMALIGN(a,b,c) POSIX_MEMALIGN_F((void **)(&a), b,c)

#define CACHE_LINE_GUARD 16 /* in units of 1 float/int - size of cache line */

#ifndef THREAD_PADDING
#define THREAD_PADDING 4096 /* in bytes - size of page */
#endif

#define MEM_MMAP 0 /* TRY USE_MIC *//* try to MMAP allocations in lightweight heap and for Fmem/Imem in qprobeval/mprobeval */

#if MEM_MMAP
extern void *malloc_huge_pages(size_t size);
extern void free_huge_pages(void *ptr);
#endif

class lightweight_heap {
protected:
#ifdef ASAN
  char **memlist;/* memlist[0..listlen-1] : list of malloc'ed memory regions */
  size_t listlen;/* number of mallocs performed */
  size_t listmax;/* size of allocated memlist[] array */
#else
  char *memory_pool;
  size_t total;/* size of per-thread memory */
  size_t total_hwm;/* total memory_pool size is total_hwm * nthreads */
#endif
  int nthreads;
	
public:

  lightweight_heap(size_t chunk_size, int numthreads) {
    if(VERB>=2){
      printf("Initializing lightweight_heap(%p):chunk_size=%llu,nthreads=%d\n",
	     this,(unsigned long long)chunk_size,numthreads);
      fflush(stdout);
    }
    if(DEBUG) assert(numthreads >= 0);
#ifdef ASAN
    nthreads = numthreads;
    listlen = 0;
    listmax = max(1,nthreads) * 4;
    memlist = (char **)malloc(listmax * sizeof(char *));
    if(!memlist){
      printf("Initializing lightweight_heap(%p):chunk_size=%llu,nthreads=%d\n",
	     this,(unsigned long long)chunk_size,numthreads);
      printf("malloc(%lu) failed\n",listmax * sizeof(char *));
      fflush(stdout);     exit(1);
    }
#else
    nthreads = 0;
    memory_pool = NULL;
    total = 0;
    total_hwm = 0;
    if(chunk_size > 0) reset_lightweight_heap(chunk_size, numthreads);
#endif
  }

#ifndef ASAN
  inline char *thread_block_start(int tid)
  {
    if(DEBUG>=2) assert(tid < nthreads);
    return(&memory_pool[tid*total_hwm]);
  }
	
  inline size_t thread_block_siz()
  {
    return total;
  }
	
  inline size_t *thread_block_free(int tid)
  {
    return((size_t *)(thread_block_start(tid)+(CACHE_LINE_GUARD*4-sizeof(size_t))));
  }
#endif
	
  void reset_lightweight_heap(size_t chunk_size, int &numthreads)
  {
#ifdef ASAN
    nthreads = numthreads;
    for(size_t i = 0; i < listlen; i++)
      free(memlist[i]);
    listlen = 0;
#else

    /* make sure there is one cache line between threads and enough memory for one free pointer.*/
    total = (chunk_size + sizeof(long) + CACHE_LINE_GUARD*4 + THREAD_PADDING-1) & (~ (THREAD_PADDING-1));

    int orig_nthreads = nthreads;
    nthreads = numthreads;

    size_t orig_hwm = total_hwm;
    total_hwm = total;
    
    if((VERB>=2) && (nthreads != orig_nthreads || total_hwm != orig_hwm)){
      printf("reset_lightweight_heap(%p): nthreads=%d -> %d, total_hwm=%llu -> %llu, memory_pool=%p\n",
	     this, orig_nthreads, nthreads, (unsigned long long)orig_hwm, (unsigned long long)total_hwm, memory_pool);
      //      dumpmemmap();
      fflush(stdout);
    }
    if(DEBUG) assert(numthreads >= 0);

    size_t orig_pool_total = orig_hwm * orig_nthreads;
    size_t pool_total = total_hwm * nthreads;

    int T = 0;

    if(pool_total > orig_pool_total || pool_total + (1ul << 28) < orig_pool_total){/* change if total increases OR decreases by at least 0.25 Gbytes */
#if MEM_MMAP
      free_huge_pages(memory_pool);
      memory_pool = (char *) malloc_huge_pages(pool_total);
#else
      free(memory_pool);
      memory_pool = (char *) malloc(pool_total);
#endif
      if(VERB>=2){
	printf("Changing lightweight_heap(%p) from %0.4f GB to %0.4f GB (per thread=%lu->%lu, nthreads=%d->%d):initial used = %d,memory_pool=%p\n",
	       this, orig_pool_total/(1024.0*1024.0*1024.0),pool_total/(1024.0*1024.0*1024.0), orig_hwm, total_hwm, orig_nthreads, nthreads, CACHE_LINE_GUARD*4, memory_pool);
	//	dumpmemmap();
	fflush(stdout);
      }
      const int retries = 10;
      for(; memory_pool==NULL && T < retries; T++){
	printf("lightweight_heap(%p):Could not allocate %d times %lu bytes:\n", this, nthreads, total);
	//	dumpmemmap();
	
	numthreads = max(1, (nthreads * 4) / 5);
	if(numthreads >= nthreads)
	  break;

	sleep(1);

	printf("Reducing number of threads from %d to %d (retry %d of %d)\n",nthreads,numthreads, T+1, retries);
	fflush(stdout);

	nthreads = numthreads;
	pool_total = total_hwm * nthreads;

#if MEM_MMAP
	memory_pool = (char *) malloc_huge_pages(pool_total);
#else
	memory_pool = (char *) malloc(pool_total);
#endif
	if(VERB/* HERE >=2 */){
	  //	  dumpmemmap();
	  fflush(stdout);
	}
      }
      if(VERB>=2){
	printf("lightweight_heap(%p):memory_pool=%p,\n",this, memory_pool);
	fflush(stdout);
      }
      if(memory_pool == NULL){
	printf("lightweight_heap(%p):Cannot allocate %d chunks of %lu bytes: Giving up after %d retries : Not enough system memory (try reduced -maxmem)\n", this,nthreads, total, retries);
	//	dumpmemmap();
	fflush(stdout); exit(1);
      }
    }
    
    for(int i=0; i < nthreads; i++)
      *thread_block_free(i) = CACHE_LINE_GUARD*4;

    if(VERB>=2 && T > 0){
      printf("At end of lighweight_heap():retries= %d\n",T);
      //      dumpmemmap();
      fflush(stdout);
    }

#endif
  }
	
  ~lightweight_heap() {
#ifdef ASAN
    for(size_t i = 0; i < listlen; i++)
      free(memlist[i]);
    free(memlist);
#else
#if MEM_MMAP
    free_huge_pages(memory_pool);
#else
    free(memory_pool);
#endif
    if(VERB>=2){
      printf("freeing lightweight_heap(%p) : pool total= %0.4f Gb, (%llu bytes per thread, nthreads=%d), memory_pool=%p\n",
	     this, total_hwm*nthreads/(1024.0*1024.0*1024.0), (unsigned long long) total_hwm, nthreads, memory_pool);
      // dumpmemmap();
      fflush(stdout);
    }
    total = 0;
    nthreads = 0;
    total_hwm = 0;
    memory_pool = NULL;
#endif
  }
	
  inline char *alloc(size_t piece_size, int tid)
  {
    if(VERB>=2){
      #pragma omp critical
      {
	printf("lightweight_heap(%p)::alloc():piece_size=%lu,tid=%d\n",this,piece_size,tid);
	fflush(stdout);
      }
    }

#ifdef ASAN
    char *ret;

    #pragma omp critical
    {
      if(listlen >= listmax){
	listmax = max(listmax*2, (size_t)(nthreads*4));
	char **new_memlist = (char **)realloc(memlist, listmax * sizeof(char *));
	if(!new_memlist){
	  printf("alloc:realloc(memlist,%llu) failed\n",(unsigned long long)(listmax * sizeof(char *)));
	  fflush(stdout);
	  assert(0);
	  exit(1);
	}
	memlist = new_memlist;
      }
      memlist[listlen] = (char *) malloc(piece_size);
      if(!memlist[listlen]){
	printf("alloc:malloc(%llu) failed\n",(unsigned long long)piece_size);
	exit(1);
      }
      ret = memlist[listlen++];
    }
    return ret;

#else
    if(tid>=nthreads) {
      #pragma omp critical
      {
        printf("lightweight_heap::alloc: INTERNAL ERROR: tid=%d > nthreads=%d\n", tid, nthreads);
	fflush(stdout); assert(0); exit(1);
      }
    }
    size_t *sfree = thread_block_free(tid);
    char *mem = thread_block_start(tid) + *sfree;	
    *sfree += piece_size;
    if(*sfree >= total) {
      #pragma omp critical
      {
	printf("lightweight_heap::alloc: INTERNAL ERROR: heap exhausted: tid=%d, used=%llu available=%llu piece_size=%llu\n", 
	       tid, (unsigned long long)(*sfree), (unsigned long long)total, (unsigned long long)piece_size);
	fflush(stdout); assert(0); exit(1);
      }
    }
    return(mem);
#endif
  }
	
  inline char *alloc(size_t piece_size, int tid, size_t alignment)
  {

    if(VERB>=2){
      #pragma omp critical
      {
	printf("lightweight_heap(%p)::alloc():piece_size=%lu,tid=%d,alignment=%lu\n",this,piece_size,tid,alignment);
	fflush(stdout);
      }
    }

#ifdef ASAN
    char *ret;

    #pragma omp critical
    {
      if(listlen >= listmax){
	listmax = max(listmax*2, (size_t)(nthreads*4));
	char **new_memlist = (char **)realloc(memlist, listmax * sizeof(char *));
	if(!new_memlist){
	  printf("alloc:realloc(memlist,%llu) failed\n",(unsigned long long)(listmax * sizeof(char *)));
	  fflush(stdout);  assert(0);	  exit(1);
	}
	memlist = new_memlist;
      }
      POSIX_MEMALIGN(memlist[listlen], alignment, piece_size);
      ret = memlist[listlen++];
    }
    return ret;
    
#else
    if(tid >= nthreads) {
      #pragma omp critical
      {
        printf("lightweight_heap::alloc(): INTERNAL ERROR: tid=%d > nthreads=%d\n", tid, nthreads);
	fflush(stdout); assert(0); exit(1);
      }
    }

    size_t *sfree = thread_block_free(tid);

    char *mem = thread_block_start(tid) + *sfree;
    int mask = (alignment -1);
    char *p = (char *)((size_t)(mem+mask) & ~mask);

    // Workaround gcc pointer arithmetic bug of some sort.
    //*sfree += (p+piece_size)-mem;
    size_t mem_tail = ((size_t)mem) & mask;
    if(mem_tail) 
      mem_tail = alignment - mem_tail;
    *sfree += piece_size + mem_tail;
    
    if(VERB>=2){
      #pragma omp critical
      {
	printf("alloc(%p):tid=%d/%d,alignment=%llu,tbase=%llu:start_use=%llu,end_use=%llu,total=%llu(hwm=%llu),piece:size=%llu,loc=%llu,tail=%lu\n",
	       this,tid,nthreads,(unsigned long long)alignment,(unsigned long long)(thread_block_start(tid) - memory_pool),(unsigned long long)(mem - thread_block_start(tid)),
	       (unsigned long long)(*sfree), (unsigned long long)total, (unsigned long long)total_hwm,(unsigned long long)piece_size,(unsigned long long)(p-thread_block_start(tid)), mem_tail);
	fflush(stdout);
      }
    }
	  
    if(*sfree >= total) {
      #pragma omp critical
      {
        printf("lightweight_heap::alloc():INTERNAL ERROR: heap exhausted:tid=%d used=%llu total=%llu alignment=%llu mem_tail=%llu piece_size=%llu\n", 
	  tid, (unsigned long long)(*sfree), (unsigned long long)total, 
	  (unsigned long long)alignment, (unsigned long long)mem_tail, (unsigned long long)piece_size);
	fflush(stdout); assert(0); exit(1);
      }
    }
    return(p);
#endif
  }

  inline char **alloc_ptr(size_t piece_size, int tid)	\
  {									\
    return (char **)alloc(piece_size*sizeof(char *), tid, CACHE_LINE_GUARD*4); \
  }
	
  #define ALLOC_FUN(type)						\
    inline type *alloc_ ## type(size_t piece_size, int tid)		\
    {									\
      return (type*)alloc(piece_size*sizeof(type), tid, CACHE_LINE_GUARD*4); \
    }
	

  ALLOC_FUN(int);
  ALLOC_FUN(float);
  ALLOC_FUN(double);
  ALLOC_FUN(long);
  ALLOC_FUN(RFLOAT);
  ALLOC_FUN(PFLOAT);
  ALLOC_FUN(FLOAT);
  ALLOC_FUN(TFLOAT);
  ALLOC_FUN(MFLOAT);
  //ALLOC_FUN(FL128);
#undef ALLOC_FUN
};

#if USE_MIC

static __m512 print_mm512(char *tag, __m512 va)
{
float a[16];
_mm512_store_ps(a, va);
fprintf(stderr, "%s=%.8g %.8g %.8g %.8g : %.8g %.8g %.8g %.8g : %.8g %.8g %.8g %.8g : %.8g %.8g %.8g %.8g\n",
	tag,
	a[0], a[1], a[2], a[3],
	a[4], a[5], a[6], a[7],
	a[8], a[9], a[10], a[11],
	a[12], a[13], a[14], a[15]);
return(va);
}


static void inline _mm512_storeu_ps(void *a, __m512 v)
{
  _mm512_packstorelo_ps(a, v);
  _mm512_packstorehi_ps(a+64, v);
}

static __m512i v512i_asc0123 = _mm512_set_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);

static __m512 inline _mm512_loadu_ps(void const *a)
{
#if 1
	__m512 v1 = _mm512_setzero_ps();
	v1 = _mm512_loadunpacklo_ps(v1, a);
	v1 = _mm512_loadunpackhi_ps(v1, a+64);
	return v1;
#else
	return _mm512_i32gather_ps(v512i_asc0123, a, 4);
#endif
}

static __m512 inline _mm512_loadu_ps(void const *a, char *base, size_t siz, int tid, int callid, int I = -1, int J = -1, int K = -1)
{
#if USE_MIC
	__m512 v1 = _mm512_setzero_ps();

	if(DEBUG>=2){
	  size_t mask = 63;
	  char *a1 = (char *)(((size_t)a) & ~mask);
   	  char *a2 = (char *)(((size_t)a1) + 128);
	  if(DEBUG && !(a1 >= base && a2 <= base+siz)){
            #pragma omp critical
	    {
	      extern lightweight_heap *rs_heap;
	      printf("_mm512_loadu_ps: tid=%d: a=%p,base=%p,siz=%ld:a1=%p,a2=%p,base+siz=%p (rs_heap:start=%p,siz=%ld,start+siz=%p) callid=%d, I=%d, J=%d, K=%d\n",
	        tid, a,base,siz,a1,a2,base+siz, rs_heap->thread_block_start(tid), rs_heap->thread_block_siz(), rs_heap->thread_block_start(tid)+rs_heap->thread_block_siz(),callid, I,J,K);
	      fflush(stdout);
	      assert(a1 >= base);
	      assert(a2 <= base+siz);
	    }
	  }
	}
	v1 = _mm512_loadunpacklo_ps(v1, a);
	v1 = _mm512_loadunpackhi_ps(v1, a+64);
	return v1;
#else
	return _mm512_i32gather_ps(v512i_asc0123, a, 4);
#endif
}

static void inline _mm512_maskstoreu_ps(void *a, __mmask16 mask16,  __m512 v)
{
  _mm512_mask_packstorelo_ps(a, mask16, v);
  _mm512_mask_packstorehi_ps(a+64, mask16, v);
}


/*
 * 
 *  BIG WARNING ! This function only produces correct results for unalinged loads with mask=000000111111b (all bits contigous low order bits)
 *  as it unpacks the elements
 * 
 *  There does not seem to exist a version of generic unaligned mask loads.
 * 
 */

static __m512 inline _mm512_maskunpackloadu_ps(void const*a, __mmask16 mask16)
{
  __m512 v1 = _mm512_setzero_ps();
  v1 = _mm512_mask_loadunpacklo_ps(v1, mask16, a);
  v1 = _mm512_mask_loadunpackhi_ps(v1, mask16, a+64);
  return v1;
}

/* substitute instruction without masking : WARNING this can segfault if reading beyond end of mmap() allocated block of memory */
static __m512 inline _mm512_maskloadu_ps(void const*a, __mmask16 mask16)
{
  return _mm512_loadu_ps(a);
}

/* substitute instruction without masking */
static __m512 inline _mm512_maskloadu_ps(void const*a, __mmask16 mask16, char *base, size_t siz, int tid, int callid, int H = -1, int J = -1, int M = -1)
{
  return _mm512_loadu_ps(a, base, siz, tid, callid, H, J, M);
}

static void inline _mm512_storeu_pd(void *a, __m512d v)
{
  _mm512_packstorelo_pd(a, v);
  _mm512_packstorehi_pd(a+64, v);
}

static __m512d inline _mm512_loadu_pd(void const *a)
{
  __m512d v1=_mm512_setzero_pd();
  v1 = _mm512_loadunpacklo_pd(v1, a);
  v1 = _mm512_loadunpackhi_pd(v1, a+64);
  return v1;
}

static void inline _mm512_maskstoreu_pd(void *a, __mmask8 mask8,  __m512d v)
{
  _mm512_mask_packstorelo_pd(a, mask8, v);
  _mm512_mask_packstorehi_pd(a+64, mask8, v);
}

static __m512d inline _mm512_maskloadu_pd(void const*a, __mmask8 mask8)
{
	
#if 0
  __m512d v1=_mm512_setzero_pd();
  v1 = _mm512_mask_loadunpacklo_pd(v1, mask8, a);
  v1 = _mm512_mask_loadunpackhi_pd(v1, mask8, a+64);
  return v1;
#endif
return _mm512_loadu_pd(a);
}

static void inline _mm512_storeu_epi32(void *a, __m512i v)
{
  _mm512_packstorelo_epi32(a, v);
  _mm512_packstorehi_epi32(a+64, v);
}

static __m512i inline _mm512_loadu_epi32(void const *a)
{
  __m512i v1=_mm512_setzero_epi32();
  v1 = _mm512_loadunpacklo_epi32(v1, a);
  v1 = _mm512_loadunpackhi_epi32(v1, a+64);
  return v1;
}

static void inline _mm512_maskstoreu_epi32(void *a, __mmask16 mask16,  __m512i v)
{
  _mm512_mask_packstorelo_epi32(a, mask16, v);
  _mm512_mask_packstorehi_epi32(a+64, mask16, v);
}

static __m512i inline _mm512_maskloadu_epi32(void const*a, __mmask16 mask16)
{
#if 0
  __m512i v1=_mm512_setzero_epi32();
  v1 = _mm512_mask_loadunpacklo_epi32(v1, mask16, a);
  v1 = _mm512_mask_loadunpackhi_epi32(v1, mask16, a+64);
  return v1;
#endif
	
return _mm512_loadu_epi32(a);
}

static __m512 inline _mm512_alignr_ps(__m512 a, __m512 b, const int s)
{
  return (__m512) _mm512_alignr_epi32((__m512i)a, (__m512i)b, s);
}

static __m512 inline _mm512_mask_alignr_ps(__m512 prev, __mmask16 mask, __m512 a, __m512 b, const int s)
{
return (__m512) _mm512_mask_alignr_epi32((__m512i)prev, mask, (__m512i)a, (__m512i)b, s);
}

static void inline _mm512_copy_few(float *a, float *b, int len)
{
#if 0
  if(len > 16) {
    printf("**** Aieee ! INTERNAL ERROR: _mm512_copy_few called with len=%d > 16\n", len);
    fflush(stdout);assert(0);exit(1);
  }
#endif

  __mmask16 mask = _mm512_int2mask((1<<len)-1);

  //  _mm512_maskstoreu_ps(a, mask, _mm512_maskloadu_ps(b, mask));
  __m512 v = _mm512_maskunpackloadu_ps(b, mask);// NEW
  _mm512_maskstoreu_ps(a, mask, v);// HERE HERE : segfaults in some cases
}

/* unaligned MIC memcpy */
static void inline  _mm512_memcpy(void *a, const void *b, size_t len)
{
  size_t cnt;

  for(cnt= 0; cnt + 64 <= len; cnt += 64) {
    __m512 v = _mm512_loadu_ps(&(((char *)b)[cnt]));
    _mm512_storeu_ps(&(((char *)a)[cnt]), v);
  }

  if(cnt >= len)
    return;

#if 1 // NEW
  for(;cnt + 4 < len; cnt += 4) {
    *(float *)&(((char *)a)[cnt]) = *(float *)&(((char *)b)[cnt]);
  }
#else
  _mm512_copy_few((float *)&(((char *)a)[cnt]), (float *)&(((char *)b)[cnt]), (len-cnt) >> 2);
  cnt = len & ~0x3;
#endif

  for(;cnt < len; cnt++) {
    ((char *)a)[cnt]=((char *)b)[cnt];
  }
}


#endif // USE_MIC


static void inline  memset_int(int *a, int val, size_t len)
{
  size_t cnt = 0;

#if (USE_MIC || USE_AVX512)
  __m512 val512 = _mm512_castsi512_ps(_mm512_set1_epi32(val));
  for(; cnt+16 < len; cnt += 16) {
    _mm512_storeu_ps(&(a[cnt]), val512);
  }
#endif
#if USE_AVX
  __m256 val256 = _mm256_castsi256_ps(_mm256_set1_epi32(val));
  for(; cnt+8 < len; cnt += 8) {
    _mm256_storeu_ps((float *) &(a[cnt]), val256);
  }
#endif
	
  for(; cnt < len; cnt++)
    a[cnt] = val;
}

static void inline  memset_float(float *a, float val, size_t len)
{
  size_t cnt = 0;

#if (USE_MIC || USE_AVX512)
  __m512 val512 = _mm512_set1_ps(val);
  for(; cnt+16 < len; cnt += 16) {
    _mm512_storeu_ps(&(a[cnt]), val512);
  }
#endif
#if USE_AVX
  __m256 val256 = _mm256_set1_ps(val);
  for(; cnt+8 < len; cnt += 8) {
    _mm256_storeu_ps(&(a[cnt]), val256);
  }
#endif
	
  for(; cnt < len; cnt++)
    a[cnt] = val;
}


extern int giter2;/* outer loop variable for -M <A> <B> */

extern int colors;/* number of colors in input files (typically based on -i inputs after processing -usecolor) */

extern pid_t pid;/* value of getpid() */

#include "xmapEntry.h"
#include "bedEntry.h"
#include "simpleRepeat.h"

extern FILE *Fopen(const char *filename, const char *mode, int maxtime);/* calls fopen(filename,mode) and repeats up to maxtime times at 1 second intervals if file cannot be opened) */
extern void FFlush(FILE *fp, int critical = 0);/* calls fflush(fp) and exits with error message, if the call fails */
extern int FILEclose(FILE *fp, int critical = 0);/* calls fclose(fp) after calling FFlush() & fsync() and exits with error message if fclose(fp) fails */

extern size_t fwrite_LOOP(void *buf, size_t size, size_t cnt, FILE *fp);

extern void printversion(FILE *fp);/* output commandline and source code version and compile time flags */

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#define MASK(n) (int)(~((~(0U))<<(n))) // rightmost n bits set to 1 (n <= 31)
#define MAXINT (int)(sizeof(int)==2 ? MASK(15) :  MASK(31))
#define MASK_LL(n) (long long)(~((~(0ULL))<<(n))) // rightmost n bits set to 1 (n <= 63)

extern int gargc;
extern char **gargv;

extern int num_blocks,MAX_BLOCKS;/* number of Cmap & site[] blocks allocated */
extern Cmap **Cmap_blocks;
extern char **site_blocks;/* memory for site[c], SNR[c],Intensity[c],Stitch[c],stitchLocation for a set of maps */
extern Cmap ****id_blocks;/* use of each block identified by a pointer to either &map or &refmap (IF these are initial allocations, to be freed by mapcompact) */
extern size_t *mapcnt_blocks;/* size of Cmap_blocks[] */

extern char **ChipIds; /* set of distinct chipid strings */
extern int numChipIds, maxChipIds;/* number of distinct chipid strings and number allocated */

extern void maxblockalloc(int num);
extern char *ChipIdAlloc(const char *chipid);

extern int startmaps;/* original number of -i input maps */
extern int totalmaps;/* If ScanCorrection > 0: the original total number of maps (may differ from nummaps if -subset was used) */
extern int rawsitemaps;/* number of maps <= totalmaps that have rawsite[] allocated */

extern long long minid,maxid;/* range of mol ids in -i input maps */

extern int nummaps;/* number of nanomaps */
extern int maxmaps;/* number of maps allocated in Gmap[] */
extern Cmap **Gmap;/* array Gmap[i=0..nummaps-1] of pointer to -i input maps */

extern int Hapnummaps;/* the value of nummaps before splitting Haplotype Query maps : the second Allele(s) are placed at refmap[orignummaps .. nummaps-1] */
extern int HapQuerymaps;
extern long long MaxQueryID;/* offset applied to Allele2 query contig id */

extern int idrenumbered;/* 1 IFF map id numbers were replaced by sequential id numbers */
extern int maptype;/* 0 == nanomap (.bnx or .vfixx), 1 == consensus (.cmap) */
extern int UniqueScans;/* For BNX Version 1.0 input : number of unique Scans : Gmap[0..nummaps-1]->UniqueScanId ranges over 0 .. UniqueScans-1 */

extern int RunDataListLen, RunDataListMax;
extern char **RunDataList;/* If != 0 : RunDataList[0..RunDataListLen-1] is the list of RunData header lines from all BNX Version 1.0 input files */

extern int numrefmaps;/* number of ref maps */
extern int maxrefmaps;/* number allocated in refmap[] */
extern Cmap **refmap;/* array refmap[i=0..numrefmaps-1] of pointer to -ref input reference maps */

extern int Hapnumrefmaps;/* the value of numrefmaps before splitting Haplotyed refmaps : the second Allele(s) are placed at refmap[orignumrefmaps .. numrefmaps-1] */
extern int Hapmaps;/* number of original ref maps that were split */
extern long long MaxContigID;/* offset applied to Allele2 contig id */

extern size_t numaligns;/* number of alignments in alignment[] */
extern size_t maxaligns;/* number of alignments allocated in alignment[] */
extern Calign **alignment;/* array alignment[i=0..numaligns-1] of map alignments */

extern int alignment_blockid;/* If >=0, alignment[] array is part of some Calign_block */

extern CalignA **alignmentA;/* see CalignA.h */
extern Calign **alignmentT;/* transposed version of alignmentT[] used by refalign */

extern size_t *numalign_start,*numalign_end;/* see globals.cpp */
extern size_t *numalign_mstart,*numalign_mend;/* see globals.cpp */


extern Cmap ** &gmap;/* Alias for (Cmap **)map in globals.h (due to name alias problem in refine()) */
extern int &gnummaps; /* Alias for (int) nummaps in globals.h (due to name alias problem in Ccontig.cpp) */

extern double draftSD;/* draft map smoothing SD as multiple of PixelLen (typically half of the res value) */
extern double draftTP;/* minimum true positive rate threshold to call draft consensus cut */

extern int numxmapentries;/* number of xmapEntry objects in xmapentries[0 .. numxmapentries-1] */
extern int maxxmapentries;/* allocated size of xmapentries[0.. maxxmapentries-1] */
extern xmapEntry **xmapentries;/* array which corresponds to entries in the xmap--filled in output_xmap */

extern int smapsize;/* number of filtered xmapEntry objects in usexmapentries[0 .. smapsize -1 ] */
extern xmapEntry **usexmapentries;/* array usexmapentries[0 .. smapsize-1] is a re-ordered subset of xmapentries[0..numxmapentries-1] */

extern int numscafentries;/* number of scaffold objects loaded from scaffold file */
class scaffold;/* declared in input_scaffold.cpp */
extern scaffold **scaffolds;/* array of (pointers to) scaffold objects */

extern int numbedentries; /* number of elements in bedentries */
extern bedEntry **bedentries;

extern int numsv; //n elements of sv_array used
extern int maxsv; //size of sv_array
class structuralVariation; //declared in structuralVariation.h
extern structuralVariation *sv_array; //see output_{x,s}map.cpp

extern int nummaskentries; //size is fixed after input, so only one int necessary
class maskSV; //declared in structuralVariation.h
extern maskSV *maskentries; //see input_mask.txt

extern int numrepeats; //number of elements of repeat_array used
extern int maxrepeats; //size of repeat_array
extern repeat *repeat_array;

extern int NumScaleFactor;
extern double *ScaleFactor;/* array of Molecule size Scaling factors tried : ScaleFactor[0..NumScaleFactor-1] */

extern void mapcorrectApply(Cmap **Map, int orignummaps, int &nummaps, int &maxmaps, int BNX, int restore, int numthreads);

extern int input_vfixx(int i, int &nummaps, int &maxmaps, Cmap **&map, std::set<long long> *ids = NULL, bool multithreaded = false, char *mbuf = NULL);
     /* input from .vfixx OR .cmap file with name vfixx_filename[i] and store result in Gmap[nummaps..maxmaps-1] */

extern int input_cmap(int fileid,int &nummaps, int &maxmaps, Cmap **&Map, int hmap, bool multithreaded = false, char *mbuf = NULL);
     /* input from .cmap file with name vfixx_filename[i] and store result in Gmap[nummaps..maxmaps-1] */

extern void check_cmap(int fileid, char *filename, int orignummaps, int &nummaps, int &maxmaps, Cmap **&Map, int HapPresent);

extern void input_vfixx_duplicates();/* check for duplicate map ids */
extern void input_refmap(char *);/* input .map file */
extern void input_spots(char*, int &nummaps, int &maxmaps, Cmap **&map, int fileid);  /* input .spots file and store result in Gmap[nummaps .. maxmaps-1] */
extern void input_align(int,char**);/* input .align file */
extern void input_scaffold(char*);/* input .scaffold file */
extern void input_bed(char*);/* input .bed file (input_bed.cpp) */
extern void input_mask(char*);/* input .mask file (input_mask.cpp) */
extern void input_confidence(char*);/* input confidence file (input_confidence.cpp) */

extern double scaffoldStartPosition(int id, double offset=0);/* see input_scaffold.cpp */
extern double scaffoldGapSize(int id1, int id2);/* see input_scaffold.cpp */
extern double scaffoldEndLength(int id, double startpos=0);/* see input_scaffold.cpp */
extern int input_groups(char *filename, Cmap **refmap, int numrefmaps);/* input group file (see -mapped) */

extern void input_svcheck(char *filename, double svcheckKB, int &numrefmaps, int &maxrefmaps, Cmap **&refmap);/* see -svcheck */
extern void indel_compare(char *filename1, char *filename2, double SVcompareKB, double SVcompareFrac);

extern void refalign_allpairs(int& nummaps, int& numrefmaps, double *pSNRtotLL, double *newPixelLen = NULL);
extern void refalign_allpairs2(int& nummaps,int& numrefmaps, double *pSNRtotLL);/* with 2-colors */

extern int pairalign_allpairs(int mapstart, int mapend, int pairmergeIter, int pairmergeIterLast);
extern int pairalign_allpairs2(int mapstart, int mapend, int pairmergeIter);/* with 2-colors */

extern double pvalue(FLOAT *Y,FLOAT *X,int n, int fn, double fp, int Outliers, double OutlierMN, int EndoutlierNL, double Ylambda, double resKB, int verbose, int Yid, int Xid, int orientation);
extern double KSPvalue(int n1, double *V1, int n2, double *V2,int verb);

extern void output_csites(int refidstart, int refidend, char* prefix, int split_maps, Calign **alignment, size_t *numalign_start, size_t *numalign_end, bool full);
extern void output_refalign(int refidstart, int refidend, char * prefix);
extern void output_refalign2(int refidstart, int refidend, char * prefix);
extern void output_chimmaps(char * prefix);
extern void output_align(char * prefix, int mapstart);
extern void output_xmap(char * prefix, char * &queryfilename, size_t align_start, size_t align_end, bool full);
extern void output_xmap2(char * prefix, size_t align_start, size_t align_end);
extern void output_cmap(char *prefix, Cmap **maps, int start, int end);
extern void output_hmap(char *prefix, Cmap **maps, int start, int end, int hmapsuffix);
extern void output_smap(char * prefix, char * &queryfilename);

#ifdef _OPENMP
#include <omp.h>
#endif

#define SPARSE_MAPWT 1 /* use sparse array (std::map) to represent mapWT[group,mapid] */

#if SPARSE_MAPWT == 1
extern void output_bnx(char *prefix, Cmap **maps, int start, int end, int raw, FILE **pFP, char **filebuffer, int tid, int g = 0, omp_lock_t *lock = NULL, int gmax = 1, int rid = -1, float *mapWT = NULL, std::map<int,float> *smapWT = NULL);
#else
extern void output_bnx(char *prefix, Cmap **maps, int start, int end, int raw, FILE **pFP, char **filebuffer, int tid, int g = 0, omp_lock_t *lock = NULL, int gmax = 1, int rid = -1, float *mapWT = NULL);
#endif


/** Main refinement function. Coded in refine.cpp
    - Perform signal processing on Ccontig.site to identify new sites.
    - Use aligned maps to edit the current Ccontig.
    - Test edit likelihood with qprobeval() before and after changes.
    - Call AddDelete() and SizesEstimate() to update sites*/
class Ccontig;
extern void refine(Ccontig *pcontig, int init, int& Lfrozen, int& Rfrozen);
extern void refineHap(Ccontig *pcontig, int **map1, int **mapK1, int **map2, int **mapK2, int& Lfrozen, int& Rfrozen);

class Cid2mapid {
public:
  long long int id;
  int mapid;
};

class Cinterleaved {/* used to convert between interleaved and non-interleaved 2-color map representation */
public:
  int c;
  int index;
};

#define ID_HASH 1 /* Use hashtable to locate id to mapid translation (faster but uses more memory than the default sort and binary search) */

extern int findid(long long id, register Cid2mapid *p, register int nummaps);
extern int findindex(long long id, register Cid2mapid *p, register int nummaps);
extern int HashFindID(long long id);/* see Assembler_input.cpp */

typedef int intcmp (const void *,const void *);

#ifndef WIN32
static inline double min(register double A,register double B)
{
  return (A <= B) ? A : B;
}

static inline double max(double A, double B)
{
  return (A >= B) ? A : B;
}

static inline float min(float A, float B)
{
  return (A <= B) ? A : B;
}

static inline float max(float A, float B)
{
  return (A >= B) ? A : B;
}

static inline size_t max(size_t A, size_t B)
{
  return (A >= B) ? A : B;
}

static inline size_t min(size_t A, size_t B)
{
  return (A <= B) ? A : B;
}

static inline int max(int A, int B)
{
  return (A >= B) ? A : B;
}

static inline int min(int A, int B)
{
  return (A <= B) ? A : B;
}

static inline long double max(long double A,long double B)
{
  return (A >= B) ? A : B;
}

static inline long long min(long long A, long long B)
{
  return (A <= B) ? A : B;
}

static inline long min(long A, long B)
{
  return (A <= B) ? A : B;
}

#endif

/* sort id in increasing order */
static inline int idinc(const Cid2mapid *p1, const Cid2mapid *p2)
{
  return (p1->id > p2->id) ? 1 : (p1->id < p2->id) ? -1 : 0;
}

#define Y_EXACT 0 /* for K >=2, more exact, by averaging over all possible combinations of internal sites being present/absent */
#define Y_APPROX SCORE_APPROX /* For K == 2, less exact, by using average of range : this approximation is used by the gradient/delta computation in refine() */

extern double mtime(),wtime();

/* compute location of condensed site formed from the K+1 consecutive sites Y[I-K ... I] */
static inline double Yc(register double *Y,
			 register int I,
			 register int K)
{
  if(Y_APPROX)
    return 0.5*(Y[I] + Y[I-K]);
  else {
    extern double FN[MAXCOLOR];
    switch(K){
    case 0:
      return Y[I];
    case 1:
      return 0.5*(Y[I]+Y[I-1]);
    case 2:
      if(Y_EXACT){
	register double fn = FN[0];
	register double p = 1.0 - fn;
	return p * Y[I-1] + fn*0.5*(Y[I]+Y[I-1]);
      } 
      return Y[I-1];
    case 3:
      if(Y_EXACT){
	register double fn = FN[0];
	register double p = 1.0 - fn;
	return p*(1.0 - p * 0.5) * (Y[I-1]+Y[I-K+1]) + fn*fn*0.5*(Y[I]+Y[I-K]);
      }
    case 4:
      if(Y_EXACT){
	register double fn = FN[0];
	register double p = 1.0 - fn;
	return p*(p*(p-2.0)+2.0)*0.5*(Y[I-1]+Y[I-K+1]) + p*fn*Y[I-2] + fn*fn*fn*0.5*(Y[I]+Y[I-K]);
      }
    default:
      if(DEBUG >= 2) assert(K>2);
      return 0.5*(Y[I-1]+Y[I-K+1]);
    }
  }
}

static inline float Yc(register float *Y,
		       register int I,
		       register int K)
{
  if(Y_APPROX)
    return 0.5f *(Y[I] + Y[I-K]);
  else {
    extern double FN[MAXCOLOR];
    switch(K){
    case 0:
      return Y[I];
    case 1:
      return 0.5f *(Y[I]+Y[I-1]);
    case 2:
      if(Y_EXACT){
	register float fn = FN[0];
	register float p = 1.0f - fn;
	return p * Y[I-1] + fn*0.5*(Y[I]+Y[I-1]);
      } 
      return Y[I-1];
    case 3:
      if(Y_EXACT){
	register float fn = FN[0];
	register float p = 1.0f - fn;
	return p*(1.0f - p * 0.5f) * (Y[I-1]+Y[I-K+1]) + fn*fn*0.5f*(Y[I]+Y[I-K]);
      }
    case 4:
      if(Y_EXACT){
	register float fn = FN[0];
	register float p = 1.0f - fn;
	return p*(p*(p-2.0f)+2.0f)*0.5f*(Y[I-1]+Y[I-K+1]) + p*fn*Y[I-2] + fn*fn*fn*0.5f*(Y[I]+Y[I-K]);
      }
    default:
      if(DEBUG >= 2) assert(K>2);
      return 0.5f*(Y[I-1]+Y[I-K+1]);
    }
  }
}

inline double Pr(double y, double resKB, double IresSD)
{
  return 0.5*(1.0+erf((y - resKB)*IresSD));
}

extern const char *VQ_TRUESITE[2][2];
extern const char *VQ_SNR[2][2];
extern const char *VQ_INTENSITY[2][2];
extern const char *VQ_STITCH[2][2];
extern const char *VQ_STITCHLOC[2][2];
extern const char *VQ_PSFWIDTH[2][2];
extern const char *VQ_IMAGELOC[2][2];

class Cmerge {
public:
  int color;
  int newid;/* index into site[color][] */
  int origid;/* index into orig[color][] */
  double site;/* (orig) location */
};

extern int siteinc(Cmerge *p1, Cmerge *p2);

extern int checkFile(char *filename);

extern long long nextContigId(char *filename);/* see output_draft.cpp */

extern char *currentDateTime();


#ifdef WIN32
// error function as found at http://stackoverflow.com/questions/6281020/error-function-erfx-not-found-in-math-h-for-visual-studio-2005
#include <float.h>
static double erf(register double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}

// the complement of the error function
static inline double erfc(double x)
{
	return 1.0 - erf(x);
}

#define isfinite(x) _finite(x)
#endif

#endif // GLOBALS_H
