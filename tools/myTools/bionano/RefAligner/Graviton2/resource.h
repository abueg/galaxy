#ifndef __RESOURCE_H__
#define __RESOURCE_H__

#include <time.h>

static Ident Resource_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/resource.h 10616 2020-02-16 10:02:32Z tanantharaman $");

#define WAIT_TIME 60 /* print message after this many seconds of waiting for resource */

#define WAIT_MAX 1 /* wait time as multiple of WAIT_INC before swapping out MatchIndex[] & MatchOverflow[] of current thread to disk */

#define WAIT_INC (USE_MIC ? 1000 : 100) /* wait time increment in milliseconds (lower values increase critical section overhead, high values add wake up delay when there is contention) */

#define TIMEOUT 180000 /* time out and over-allocate after this many seconds of waiting for resource (risk of swapping) : only for critical section */
#define TIMEOUT2 360000 /* time out and over-allocate after this many seconds of waiting for resource (risk of swapping) : all allocations  */

extern double wtime();

#ifdef WIN32

struct timespec {
  int tv_sec;
  int tv_nsec;
};

void nanosleep(timespec *tm, void *)
{
  Sleep(tm->tv_sec * 1000 + (tm->tv_nsec / 1000000));
}

#endif

// #include <mcheck.h>

/* read back MatchIndex[] & MatchOverflow[] memory of current thread from disk : called inside a critical section so VmRSS is raised back before next thread checks it */
template<class OFFSET_T, class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
static void MemRestore(FILE *FPswap, size_t memsiz, CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *MT, int tid, size_t IndexSiz, size_t Imemsiz,
		       int i, struct timespec& tm, char *filename /* these args are only used for debugging or verbose output */   
)
{
  double wt = wtime();
  if(VERB>=2){
    printf("tid=%d: After waiting %0.2f secs swapping back in %lu+%lu bytes from %s: wt= %0.6f\n",
	   tid, 1e-9 * i * tm.tv_nsec, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * IndexSiz, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * (size_t)MT->MatchOverflowNext, filename,wt);
    fflush(stdout);
  }
  //    mcheck_check_all();

  fseek(FPswap, 0L, SEEK_SET);

  size_t r;
  if((r = fread(MT->MatchIndex, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>), IndexSiz, FPswap)) != IndexSiz){
    printf("tid=%d:Failed to read %lu bytes from %s (r= %lu, IndexSiz= %lu)\n",tid, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * IndexSiz, filename, r, IndexSiz);
    fflush(stdout);exit(1);
  }

  if(MT->MatchOverflowNext > 0){
    if((r = fread(MT->MatchOverflow, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>), (size_t)MT->MatchOverflowNext, FPswap)) != (size_t)MT->MatchOverflowNext){
      printf("tid=%d:Failed to read %lu bytes from %s (r= %lu, MatchOverflowNext= %u)\n",tid, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * (size_t) MT->MatchOverflowNext, filename, r, MT->MatchOverflowNext);
      fflush(stdout);exit(1);
    }
  }

  fclose(FPswap);
  unlink(filename);

  if(VERB>=2){
    double nwt = wtime();
    printf("tid=%d: After waiting %0.2f secs swapped back in %0.4f + %0.4f G (%0.2f Mb/sec) from %s: wall time= %0.6f (Elapsed= %0.6f)\n",
	   tid, 1e-9 * i * tm.tv_nsec, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * IndexSiz * 1e-9, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * (size_t)MT->MatchOverflowNext * 1e-9, 
	   (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * IndexSiz + sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * (size_t)MT->MatchOverflowNext) * 1e-6 / (nwt-wt), filename,nwt-wt,nwt);
    fflush(stdout);
  }

  //    mcheck_check_all();
}

template<class OFFSET_T,class OFFSET1_T, int HASHWIN, int RANGE, int HASH_SCALE>
class Resource {
  long long max_amount;
  long long amount;
  int allocations, active_threads;
  int max_threads, *thread_allocations;
  long long *thread_total;
  size_t IndexSiz, Imemsiz;
  long long VmSize, VmRSS, VmSwap, VmHWM;
  size_t blen;

  struct timespec tm;

 public:
  Resource(long long m, int numthreads)
    {
      blen = USE_MIC ? 10*1024*1024 : 64*1024; // If USE_MIC : must match scif_io.h buffer size to avoid triggering bug (failed fwrite())

      max_amount = m;
      amount = m;
      allocations = active_threads = 0;
      tm.tv_sec = 0;
      tm.tv_nsec = WAIT_INC * 1000 * 1000;/* WAIT_INC ms */
      max_threads = numthreads;
      thread_allocations = new int[numthreads];
      thread_total = new long long[numthreads];
      for(int i = 0; i < numthreads; i++){
	thread_allocations[i] = 0;
	thread_total[i] = 0;
      }
      IndexSiz = (1UL << MHASH_BITS);
      Imemsiz = (sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>) * IndexSiz + PAGE-1) & ~(PAGE-1);
      VmSize = VmRSS = VmSwap = VmHWM = 0;
      getmem(VmSize,VmRSS,VmSwap,&VmHWM);
    }
	
  ~Resource()
   {
     if(VERB>=2 || max_amount > amount){
       if(max_amount > amount)
	 printf("WARNING:Deleting Resource: amount=%0.4f,max_amount=%0.4f G(used=%0.4f G),allocations=%d\n",amount*1e-9,max_amount*1e-9,(max_amount-amount)*1e-9,allocations);
       else
	 printf("Deleting Resource: amount=%0.4f,max_amount=%0.4f G(used=%0.4f G),allocations=%d\n",amount*1e-9,max_amount*1e-9,(max_amount-amount)*1e-9,allocations);
       if(max_amount > amount){
	 long long total = 0;
	 for(int i = 0; i < max_threads; i++){
	   total += thread_total[i];
	   printf("tid=%d/%d: thread_total[tid]= %0.4f G, cum= %0.4f G\n",i, max_threads, thread_total[i]*1e-9, total*1e-9);
	 }
       }
       fflush(stdout);
     }
     delete [] thread_allocations;
     delete [] thread_total;
   }
	
  /* The last allocation (in MatchFind()) of each thread is marked as critical==1 : At least one thread must be past this allocation (active_threads > 0), before any allocation is blocked */
  /* If MatchOverflow != NULL and wait time has exceeded WAIT_MAX seconds, save MatchOverflow[0..MatchOverflowNext-1] to disk and apply madvise(DONT_NEED) to MatchOverflow. When memory becomes available,
     read the data back in from disk */
     
  void allocate(long long amt, int critical, size_t memsiz, CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *MT) // unsigned int MatchOverflowNext, CMatchLine* MatchOverflow)
  {
    FILE *FPswap = NULL;
    char *write_buffer = NULL;
    long long memswapped = 0;
    char filename[PATH_MAX];

    int tid = MULTIMATCH_TWOTHREADS ? MT->tid * 2 + (MT->matchtable ? 1 : 0) : MT->tid;

    if((DEBUG && !(0 <= tid && tid < max_threads)) || VERB>=2){
      #pragma omp critical
      {
	printf("allocate:amt=%0.4f G,critical=%d,memsiz=%0.4f G, allocated= %0.4f G, MT->tid=%d, MT->matchtable= %p : tid=%d, max_threads=%d, remaining mem=%0.4f G:wtime= %0.6f\n",
	       amt*1e-9,critical,memsiz*1e-9,thread_total[tid]*1e-9, MT->tid,MT->matchtable, tid, max_threads,amount *1e-9,wtime());
	fflush(stdout);
	assert(0 <= tid && tid < max_threads);
      }
    }

    if(amt > max_amount) {
      #pragma omp critical
      {
	printf("WARNING: amount requested (%0.4f G) is larger than maximum amount (%0.4f G). Thread %d may block until it is the only thread\n", amt*1e-9, max_amount*1e-9, tid);
	fflush(stderr);
      }
    }

    int ret = 0;

    long long Hmaxmem = HMaxMem * (1000LL * 1000LL * 1000LL);

    for(long long i = 0;;i++){

      #pragma omp critical(resource)
      {
	if(amount > amt || active_threads <= 0 || VmRSS + VmSwap + amt + memswapped <= Hmaxmem){
	  amount -= amt;
	  thread_total[tid] += amt;
	  allocations++;
	  if(critical){
	    if(DEBUG) assert(thread_allocations[tid] >= 0);
	    if(!thread_allocations[tid])
	      active_threads++;
	    thread_allocations[tid]++;
	  }
	  if(VERB && i > 0){
	    printf("tid=%d:%0.4f/%0.4f G received after %0.2f secs (remaining=%0.4f,max=%0.4f G,allocations=%d,%d,active threads=%d,cr=%d,VmRSS=%0.3f,VmSwap=%0.3f,VmHWM=%0.3f,Max=%0.3f,Disk=%0.3f): wall time=%0.6f \n",
	      tid,amt*1e-9,thread_total[tid]*1e-9,1e-9 * i * tm.tv_nsec,amount*1e-9,max_amount*1e-9,thread_allocations[tid],allocations,active_threads,critical,VmRSS*1e-9,VmSwap*1e-9,VmHWM*1e-9,Hmaxmem*1e-9,memswapped*1e-9,wtime());
	    fflush(stdout);
	  }
	  VmRSS += amt + memswapped;/* conservative : forces other threads to recheck VmRSS, if needed */
	  ret = 1;
	  if(FPswap != NULL) {/* read back MatchOverflow memory from disk */
            MemRestore(FPswap, memsiz, MT, tid, IndexSiz, Imemsiz, i, tm, filename);
	    free(write_buffer);
          }
        }
      }
      if(ret)
	return;

      // NEW3 : before checking VmRSS, see if MatchOverflow[] & MatchIndex[] array can be swapped out to disk
      if(MT->MatchOverflow && !FPswap && i >= WAIT_MAX){

	sprintf(filename,"%s_p%lld_T%d.Hswap", output_prefix, (long long)pid, tid);

	double wt = wtime();
	if(VERB>=2){
	  printf("tid=%d: After waiting %0.2f secs for %0.4f G, will swap out %0.4f G to %s and free up %0.4f G + %0.4f G of memory:wt= %0.6f\n",
		 tid, 1e-9 * i * tm.tv_nsec, amt*1e-9, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>)*(size_t)MT->MatchOverflowNext*1e-9,filename,Imemsiz*1e-9,memsiz*1e-9,wt);
	  fflush(stdout);
	}
	// mcheck_check_all();

	if(DEBUG && strlen(filename) >= PATH_MAX){
	  printf("filename = %s is too long (%lu chars), must not exceed %d chars\n",filename, strlen(filename),PATH_MAX-1);
	  fflush(stdout);exit(1);
	}
	unlink(filename);/* in case it is present from previous run */

	if((FPswap = fopen(filename,"w+")) == NULL){
	  int eno = errno;
	  char *err = strerror(eno);
	  #pragma omp critical
	  {
	    printf("tid=%d:failed to open file %s for writing:errno=%d:%s\n",tid,filename,eno,err);
	    fflush(stdout);exit(1);
	  }
	}
	if(!(write_buffer = (char *)malloc(blen))){
          printf("allocate:malloc(%llu) failed\n",(unsigned long long)blen);
          fflush(stdout);exit(1);
        }
        setbuffer(FPswap, write_buffer, blen);

	size_t n = fwrite_LOOP(MT->MatchIndex, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>), IndexSiz, FPswap);
	if(n != IndexSiz){
	  #pragma omp critical
	  {
	    printf("tid=%d:fwrite= %lu: write of %lu bytes to %s failed\n", tid, n, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>)*IndexSiz, filename);
	    fflush(stdout);exit(1);
	  }
	}
	memswapped += IndexSiz * sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>);
	
	if(DEBUG) assert(MT->MatchOverflow != NULL);

	if(MT->MatchOverflowNext > 0){
	  n = fwrite_LOOP(MT->MatchOverflow, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>), (size_t)MT->MatchOverflowNext, FPswap);
	  if(n != (size_t)MT->MatchOverflowNext){
	    #pragma omp critical
	    {
	      printf("tid=%d:fwrite= %lu: write of %lu bytes to %s failed\n", tid, n, sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>)*(size_t)MT->MatchOverflowNext, filename);
	      fflush(stdout);exit(1);
	    }
	  }
	  memswapped += (size_t) MT->MatchOverflowNext * sizeof(CMatchLine<OFFSET_T,RANGE,HASH_SCALE>);
	}
	fflush(FPswap);
	
	if(madvise(MT->MatchIndex, Imemsiz, MADV_DONTNEED)){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("madvise(%p,%lu,MADV_DONTNEED) failed:errno=%d:%s\n",MT->MatchIndex,memsiz,eno,err);
	  fflush(stdout);exit(1);
	}
	if(madvise(MT->MatchOverflow, memsiz, MADV_DONTNEED)){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("madvise(%p,%lu,MADV_DONTNEED) failed:errno=%d:%s\n",MT->MatchOverflow,memsiz,eno,err);
	  fflush(stdout);exit(1);
	}
	
        if(VERB>=2){
          #pragma omp critical
          {
	    double nwt = wtime();

            long long rVmSize, rVmRSS, rVmSwap, rVmHWM;
	    getmem(rVmSize,rVmRSS,rVmSwap,&rVmHWM);

            printf("tid=%d: After waiting %0.2f secs for %0.4f G, swapped out %0.4f G (%0.2f Mb/sec) to %s, freed %0.4f G + %0.4f G: rVmRSS=%0.3f,rVmSwap=%0.3f,rmHWM=%0.3f,wall time= %0.6f (elapsed= %0.6f)\n",
	       tid, 1e-9 * i * tm.tv_nsec, amt*1e-9, memswapped*1e-9, memswapped * 1e-6 /(nwt-wt), filename,Imemsiz*1e-9,memsiz*1e-9,rVmRSS*1e-9,rVmSwap*1e-9,rVmHWM*1e-9,nwt-wt,nwt);
	    fflush(stdout);
          }
        }
      }

      // mcheck_check_all();

      // NEW3 : check VmRSS etc to see if there is enough memory available
      #pragma omp critical(resource)
      {
	getmem(VmSize,VmRSS,VmSwap,&VmHWM);
	if(DEBUG) assert(VmHWM >= VmRSS);
	if(VmRSS + VmSwap + amt + memswapped <= Hmaxmem){
	  amount -= amt;
	  thread_total[tid] += amt;
	  allocations++;
	  if(critical){
	    if(DEBUG) assert(thread_allocations[tid] >= 0);
	    if(!thread_allocations[tid])
	      active_threads++;
	    thread_allocations[tid]++;
	  }
	  if(VERB && i > 0 && i >= (WAIT_TIME * 1000 / WAIT_INC)){
	    printf("tid=%d:%0.4f/%0.4f G received after %0.2f secs (remaining=%0.4f,max=%0.4f G,allocations=%d,%d,active threads=%d,cr=%d,VmRSS=%0.3f,VmSwap=%0.3f,VmHWM=%0.3f,Max=%0.3f,Disk=%0.3f Gb): wall time=%0.6f \n",
  	        tid,amt*1e-9,thread_total[tid]*1e-9,1e-9 * i * tm.tv_nsec,amount*1e-9,max_amount*1e-9,thread_allocations[tid],allocations,active_threads,critical,VmRSS*1e-9,VmSwap*1e-9,VmHWM*1e-9,Hmaxmem*1e-9,memswapped*1e-9,wtime());
	    fflush(stdout);
	  }
	  VmRSS += amt + memswapped;/* conservative : forces other threads to recheck VmRSS, if needed */
	  ret = 1;
	  if(FPswap != NULL){/* read back MatchIndex[] & MatchOverflow[] memory of current thread from disk */
	    MemRestore(FPswap, memsiz, MT, tid, IndexSiz, Imemsiz, i, tm, filename);
	    free(write_buffer);
          }
	}
      }
      if(ret)
	return; 

      if(VERB && (VERB>=2 || !(i % (WAIT_TIME * 1000 / WAIT_INC)))){
        #pragma omp critical(stdout_crit)
	{
	  if(i==0)
	    printf("tid=%d: Starting wait for %0.4f G (available=%0.4f,max=%0.4f,swap=%0.4f G,allocations=%d,%d,active threads=%d,cr=%d, VmRSS=%0.3f,VmSwap=%0.3f,VmHWM=%0.3f,Max=%0.3f Gb): wall time=%0.6f\n",
	  tid,amt*1e-9,amount*1e-9,max_amount*1e-9,memswapped*1e-9,thread_allocations[tid],allocations,active_threads,critical,VmRSS*1e-9,VmSwap*1e-9,VmHWM*1e-9,Hmaxmem*1e-9,wtime());
	  else
	    printf("tid=%d: Waited %0.2f secs for %0.4f G (available=%0.4f,max=%0.4f,swap=%0.4f G,allocations=%d,%d,active threads=%d,cr=%d, VmRSS=%0.3f,VmSwap=%0.3f,VmHWM=%0.3f,Max=%0.3f Gb): wall time=%0.6f\n",
	  tid, 1e-9 * i * tm.tv_nsec, amt*1e-9,amount*1e-9,max_amount*1e-9,memswapped*1e-9,thread_allocations[tid],allocations,active_threads,critical,VmRSS*1e-9,VmSwap*1e-9,VmHWM*1e-9,Hmaxmem*1e-9,wtime());
	  fflush(stdout);
	}
      }

      if(i >= (critical ? TIMEOUT : TIMEOUT2) * 1000 / WAIT_INC){
        #pragma omp critical(resource)
	{
	  if(VERB){
	    printf("WARNING: tid=%d:Waited %0.2f secs for %0.4f G (available=%0.4f/%0.4f G,allocations=%d,%d,active threads=%d,cr=%d) : over allocating : wall time=%0.6f\n",
		   tid, 1e-9 * i * tm.tv_nsec, amt*1e-9, amount*1e-9, max_amount*1e-9,thread_allocations[tid],allocations,active_threads,critical,wtime());
	    fflush(stdout);
	  }

	  amount -= amt;
	  thread_total[tid] += amt;
	  allocations++;
	  if(critical){
	    if(DEBUG) assert(thread_allocations[tid] >= 0);
	    if(DEBUG) assert(active_threads >= 0);
	    if(thread_allocations[tid] <= 0)
	      active_threads++;
	    thread_allocations[tid]++;
	  }
	  VmRSS += amt + memswapped;/* conservative : forces other threads to recheck VmRSS, if needed */
	  if(FPswap != NULL){/* read back MatchIndex[] & MatchOverflow[] memory of current thread from disk */
	    MemRestore(FPswap, memsiz, MT, tid, IndexSiz, Imemsiz, i, tm, filename);
	    free(write_buffer);
          }
        }

	return;
      }
      struct timespec rem;
      nanosleep(&tm,&rem);/* sleep for WAIT_INC ms */
    }/* wait loop : for(;;i++) */
  }
	
  void ResFree(long long amt, int critical, CMatchTable<OFFSET_T,OFFSET1_T,HASHWIN,RANGE,HASH_SCALE> *MT)
  {
    int tid = MULTIMATCH_TWOTHREADS ? MT->tid * 2 + (MT->matchtable ? 1 : 0) : MT->tid;
    if(DEBUG/* HERE >= 2*/) assert(0 <= tid && tid < max_threads);

    #pragma omp critical(resource)
    {
      amount += amt;
      allocations--;
      if(critical){
	thread_allocations[tid]--;
	if(thread_allocations[tid] <= 0)
	  active_threads--;
	if(DEBUG) assert(thread_allocations[tid] >= 0);
	if(DEBUG) assert(active_threads >= 0);
      }
      if(VERB>=2){
        long long rVmSize, rVmRSS, rVmSwap, rVmHWM;
	getmem(rVmSize,rVmRSS,rVmSwap,&rVmHWM);

      	printf("tid=%d:released %0.4f/%0.4f G: (available=%0.4f max=%0.4f G,allocations=%d,%d,active threads=%d,cr=%d):VmRSS=%0.3f,VmSwap=%0.3f G:wall time=%0.6f\n",
          tid,amt*1e-9,thread_total[tid]*1e-9,amount*1e-9,max_amount*1e-9,thread_allocations[tid],allocations,active_threads,critical,rVmRSS*1e-9,rVmSwap*1e-9,wtime());
      	fflush(stdout);
      }
      thread_total[tid] -= amt;
    }
  }
	
};


#endif
