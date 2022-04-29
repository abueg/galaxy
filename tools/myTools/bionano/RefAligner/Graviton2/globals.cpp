#include "globals.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/globals.cpp 11531 2020-08-22 22:23:24Z tanantharaman $");

#include "version.h"

#if USE_MIC
#include "scif_io.c"
#endif

#include <time.h>

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
static char   timebuf[80];

char *currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(timebuf, sizeof(timebuf), "%Y-%m-%d.%X", &tstruct);

    return timebuf;
}

/* same as fopen() but loop up to maxcnt (typically 65) seconds in case NFS client directory attribute cache is out of date (should not be needed if acdirmax=0) */

FILE *Fopen(const char *filename,  const char *mode, int maxcnt)
{
  double starttime = wtime(),maxtime = maxcnt+0.999;
  errno = 0;
  FILE *fp = NULL;
  for(int cnt = 0; cnt < maxcnt && (fp = fopen(filename,mode)) == NULL; cnt++){
    int eno = errno;
    char *err = strerror(eno);
    double w = wtime();
    printf("Cannot open file %s, retry=%d/%d:errno=%d:%s (Current Time= %s), elapsed time= %0.3f secs\n",filename,cnt,maxcnt,eno,err,currentDateTime(),w-starttime);
    fflush(stdout);

    if(w - starttime >= maxtime){
      printf("Elapsed time exceeded limit of %0.1f secs\n",maxtime);
      fflush(stdout);
      break;
    }

    errno = 0;
    sleep(1);
  }
  
  return fp;
}

/* same as fflush() and checks return code and prints error message and exits, if needed */
void FFlush(FILE *fp, int critical)
{
  int eno, r;
  if(!USE_MIC && critical){
    #pragma omp critical(FFlushlock1)
    {
      errno = 0;
      r = fflush(fp);
      eno = errno;
    }
  } else {
    errno = 0;
    r = fflush(fp);
    eno = errno;
  }

  if(r){
    char *err = strerror(eno);
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num ();
#endif
    #pragma omp critical(FFflushlock2)
    {
      printf("ERROR:tid=%d,fflush(fp=%p)=%d(failed),critical=%d,errno=%d: %s\n", tid, fp, r, critical, eno, err);
      fflush(stdout); exit(1);
    }
  }
  return;
}

#define SCIF_MAX 10*1024*1024 // need to limit fwrite() calls to no more than this many bytes to avoid errors in the scif implementation of fwrite()

size_t fwrite_LOOP(void *buf, size_t size, size_t num, FILE *fp)
{
  char *q = (char *)buf;
  size_t maxreq = SCIF_MAX / size;
  size_t cnt = 0;
  while(cnt < num){
    errno = 0;
    char *p = &q[cnt * size];
    size_t req = min(num - cnt, maxreq);
    size_t r = fwrite((void *)p, size, req, fp);
    if(r <= 0){
      int eno = errno;
      char *err = strerror(eno);
      printf("fwrite(&buf[%lu],%lu,%lu) failed (ret=%lu):errno=%d:%s\n",p-q,size,req,r,eno,err);
      fflush(stdout);
      return cnt;
    }
    if(r < req){
      int eno = errno;
      char *err = strerror(eno);
      printf("WARNING: fwrite(&buf[%lu],%lu,%lu) came up short (ret=%lu):errno=%d:%s\n",p-q,size,req,r,eno,err);
      printf("    Will continue: completed %lu, remaining=%lu elements of %lu bytes\n",cnt+r,num-(cnt+r),size);
      fflush(stdout);
    } else if(VERB>=2){
      printf("fwrite(&buf[%lu],%lu,%lu) succeeeded (ret=%lu): completed %lu elements\n",p-q,size,req,r,cnt+r);
      fflush(stdout);
    }
    cnt += r;
  }
  if(DEBUG) assert(cnt == num);
  return cnt;
}

int FILEclose(FILE *fp, int critical)
{
  FFlush(fp, critical);

#ifndef WIN32
#if USE_MIC == 0
  errno = 0;
  int fd = fileno(fp);
  if(fd < 0){
    char *err = strerror(errno);    
    fprintf(stdout,"ERROR:fileno(fp) returned -1: %s\n",err);
    fflush(stdout);   exit(1);
  } else {
    errno = 0;
    if(fsync(fd) < 0){
      char *err = strerror(errno);
      printf("ERROR:fsync(fd=%d) returned -1:%s\n",fd,err);
      fflush(stdout);   exit(1);
    }
  }
#endif
#endif

  errno = 0;
  int ret = fclose(fp);
  if(ret){
    char *err = strerror(errno);    
    #pragma omp critical
    {
      printf("ERROR:fclose(fp=%p) returned %d: %s\n",fp,ret,err);
      fflush(stdout);    exit(1);
    }
  }
  return ret;
}

void printversion(FILE *fp)
{
  char cwd[BUFSIZ];

#ifndef WIN32
  gethostname(cwd,BUFSIZ);

  fprintf(fp,"# hostname=%s\n", cwd);

#endif

  fprintf(fp,"# $ cd %s; %s", getcwd(cwd,BUFSIZ), gargv[0]);
  for(register int i = 1; i < gargc; i++)
    fprintf(fp," %s", gargv[i]);
  fprintf(fp,"\n");


#if 1 // OLD version
  (void)fprintf(fp,"# %s %s\n", VERSION, SVN_ID);
  (void)fprintf(fp,"# FLAGS: USE_SSE=%d USE_AVX=%d USE_MIC=%d USE_PFLOAT=%d USE_RFLOAT=%d USE_MFLOAT=%d USE_EPOW=%d DEBUG=%d VERB=%d\n", USE_SSE, USE_AVX, USE_MIC, USE_PFLOAT, USE_RFLOAT, USE_MFLOAT, USE_EPOW, DEBUG, VERB);
#else // NEW version
  (void)fprintf(fp,"# %s\n", COMPILE_COMMAND);
  (void)fprintf(fp,"# FLAGS: USE_SSE=%d USE_AVX=%d USE_MIC=%d USE_PFLOAT=%d USE_RFLOAT=%d USE_MFLOAT=%d USE_EPOW=%d DEBUG=%d VERB=%d\n", USE_SSE, USE_AVX, USE_MIC, USE_PFLOAT, USE_RFLOAT, USE_MFLOAT, USE_EPOW, DEBUG, VERB);
  (void)fprintf(fp,"# %s %s\n", VERSION, SVN_ID);
#endif
  fflush(fp);
}

void dumpmemmap()
{
#ifndef WIN32
  if(!(RELEASE && USE_MIC)){
    char filename[BUFSIZ];
    pid_t pid = getpid();
    sprintf(filename,"/proc/%lld/status",(long long)pid);

    printf("---------Start of /proc/%lld/status ----------- wall time=%0.6f\n", (long long)pid, wtime());
    fflush(stdout);

    errno = 0;

#if USE_MIC
#undef fopen
#endif
    FILE *fp = fopen(filename,"r");// Need to use real fopen() NOT sic_fopen in MIC
#if USE_MIC
#define fopen sic_fopen
#endif

    if(fp == NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("dumpmem(): Failed to read %s: errno=%d:%s\n",filename,eno,err);
      fflush(stdout);
    } else {
      char buf[LINESIZ];
      int linecnt = 1;
      for(;fgets(buf,LINESIZ,fp) != NULL; linecnt++){
	size_t len = strlen(buf);
	if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	  if(len >= LINESIZ-1)
	    printf("Line too long: exceeds max length=%d (increase LINESIZE in constants.h) in input file %s on line %d:\n%s\n",LINESIZ-1,filename,linecnt,buf);
	  else
	    printf("Line not correctly terminated (last char = %c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
	  fflush(stdout); 
	  break;
	}
	buf[len-1] = '\0';
	printf("%s\n",buf);
      }
      fclose(fp);
    }

    printf("---------End of /proc/%lld/status -------------\n", (long long)pid);
    fflush(stdout);

    sprintf(filename,"/proc/meminfo");
  
    printf("---------Start of /proc/meminfo -----------\n");
    fflush(stdout);

    errno = 0;

#if USE_MIC
#undef fopen
#endif
    fp = fopen(filename,"r");// Need to use real fopen() NOT sic_fopen in MIC
#if USE_MIC
#define fopen sic_fopen
#endif

    if(fp == NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("dumpmem(): Failed to read %s: errno=%d:%s\n",filename,eno,err);
      fflush(stdout);
    } else {
      char buf[LINESIZ];
      int linecnt = 1;
      for(;fgets(buf,LINESIZ,fp) != NULL; linecnt++){
	size_t len = strlen(buf);
	if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	  if(len >= LINESIZ-1)
	    printf("Line too long: exceeds max length=%d (increase LINESIZE in constants.h) in input file %s on line %d:\n%s\n",LINESIZ-1,filename,linecnt,buf);
	  else
	    printf("Line not correctly terminated (last char = %c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
	  fflush(stdout); 
	  break;
	}
	buf[len-1] = '\0';
	printf("%s\n",buf);
      }
      fclose(fp);
    }

    printf("---------End of /proc/meminfo -----------\n");
    fflush(stdout);
  }// if(!(RELEASE && USE_MIC))

#if 0
  if(!RELEASE){
    sprintf(buf,"ps auxwww");

    err = system(buf);
    if(err){
      printf("system(%s) failed: return code = %d\n", buf, err);
      fflush(stdout);
    }
  }
#endif
    
#endif // ifndef WIN32
}

#if 0 // historical Linux code to find out current malloc free space
long alloc_ttl, alloc_max;

void available_alloc (void)
{
  struct __mp__ MTYP *q;

  alloc_ttl = 0;
  alloc_max = 0;

  q = (struct __mp__ MTYP *) &__mp__;
  q = q->next;

  while (q) {
    alloc_ttl += q->len;
    if (alloc_max < q->len) allox_max = q->len;
    q = q->next;
  }
}
#endif

long long OMP_stack_size = 2*1024*1024;/* this value should be updated at run time by setting OMP_STACKSIZE */
static int stack_size_init = 0;

void getmem(long long &VmSize, long long &VmRSS, long long &VmSwap, long long *pVmHWM /* If != 0 */, long long *pVmPeak /* If != 0 */)
{

#ifdef WIN32
  // make up some reasonable values
  VmSize = 2*1024*1024*1024;
  VmRSS =     128*1024*1024;
  VmSwap = 0;
#else
  if(!stack_size_init){/* try to read OMP_STACKSIZE */
    char *OMP_STACKSIZE = getenv("OMP_STACKSIZE");
    if(OMP_STACKSIZE){
      char *pt = OMP_STACKSIZE,*qt;
      long long size = strtol(pt,&qt,10);
      if(pt == qt){
	printf("Invalid OMP_STACKSIZE=%s : Assuming 2M stacksize\n",pt);
	fflush(stdout);
      } else {
	size *= 1024;
	if(*qt == 'M' || *qt == 'm')
	  size *= 1024;
	else if(*qt == 'G' || *qt == 'g')
	  size *= 1024*1024;
	else if(*qt == 'B' || *qt == 'b')
	  size /= 1024;
	if(VERB/* HERE >= 2*/){
	  printf("OMP_STACKSIZE = %lld bytes\n",size);
	  fflush(stdout);
	}
	OMP_stack_size = size;
      }
      stack_size_init = 1;
    }
  }

  char filename[PATH_MAX];
  pid_t pid = getpid();
  sprintf(filename,"/proc/%lld/status",(long long)pid);

#if USE_MIC
#undef fopen
#endif
  FILE *fp = fopen(filename,"r");// Need to use real fopen() NOT sic_fopen in MIC
#if USE_MIC
#define fopen sic_fopen
#endif

  if(fp==NULL){
    char *err = strerror(errno);
    printf("getmem(): Failed to read %s: %s\n",filename,err);
    fflush(stdout);
    char buf[BUFSIZ];
    sprintf(buf,"ls -l /proc/%lld",(long long)pid);
    printf("$ %s\n",buf);
    fflush(stdout);
    int ret = system(buf);
    if(ret){
      printf("system(%s) failed: return code = %d\n", buf, ret);
      fflush(stdout);
    }
    sprintf(buf,"cat /proc/%lld/status",(long long)pid);
    printf("$ %s\n",buf);
    fflush(stdout);
    ret = system(buf);
    if(ret){
      printf("system(%s) failed: return code = %d\n", buf, ret);
      fflush(stdout);
    }
    dumpmemmap();
    fflush(stdout);  exit(1);
  }

  char buf[LINESIZ];
  int linecnt = 1;
  
  long long nVmRSS = -1, nVmSize = -1, nVmSwap = -1, nVmPeak = -1, nVmHWM = -1;// initialize as not read

  for(;fgets(buf,LINESIZ,fp) != NULL;linecnt++){
    size_t len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      if(len >= LINESIZ-1)
	printf("Line too long: exceeds max length=%d (increase LINESIZE in constants.h) in input file %s on line %d:\n%s\n",LINESIZ-1,filename,linecnt,buf);
      else
	printf("Line not correctly terminated (last char = %c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      fflush(stdout); exit(1);
    }

    char *pt, *qt, *key;

    if(pVmPeak != NULL){
      key = (char *)"VmPeak:";
      if((pt = strstr(buf,key))){
	pt += strlen(key);
	nVmPeak = strtoll(pt,&qt,10);
	if(pt == qt || !(*qt == 0 || isspace(*qt))){
	  printf("value of %s on line %d of %s not recognized:\n%s\n",key,linecnt,filename,buf);
	  fflush(stdout); exit(1);
	}
	pt = qt;
	while(*pt && isspace(*pt))
	  pt++;
	if(strncmp(pt,"kB",2)){
	  printf("WARNING:Unknown memory size unit (%s) on line %d of %s (expected kB):\n%s\n",pt,linecnt,filename,buf);
	  fflush(stdout);
	}
	nVmPeak *= 1024LL;
	if(VERB>=2){
	  printf("VmSize=%lld:%s",nVmSize,buf);
	  fflush(stdout);
	}
	continue;// next line
      }
    }

    key = (char *)"VmSize:";
    if((pt = strstr(buf,key))){
      pt += strlen(key);
      nVmSize = strtoll(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt))){
	printf("value of %s on line %d of %s not recognized:\n%s\n",key,linecnt,filename,buf);
	fflush(stdout); exit(1);
      }
      pt = qt;
      while(*pt && isspace(*pt))
	pt++;
      if(strncmp(pt,"kB",2)){
	printf("WARNING:Unknown memory size unit (%s) on line %d of %s (expected kB):\n%s\n",pt,linecnt,filename,buf);
	fflush(stdout);
      }
      nVmSize *= 1024LL;
      if(VERB>=2){
	printf("VmSize=%lld:%s",nVmSize,buf);
	fflush(stdout);
      }
      continue;// next line
    }

    if(pVmHWM != NULL){
      key = (char *)"VmHWM:";
      if((pt = strstr(buf,key))){
	pt += strlen(key);
	nVmHWM = strtoll(pt,&qt,10);
	if(pt == qt || !(*qt == 0 || isspace(*qt))){
	  printf("value of %s on line %d of %s not recognized:\n%s\n",key,linecnt,filename,buf);
	  fflush(stdout); exit(1);
	}
	pt = qt;
	while(*pt && isspace(*pt))
	  pt++;
	if(strncmp(pt,"kB",2)){
	  printf("WARNING:Unknown memory size unit (%s) on line %d of %s (expected kB):\n%s\n",pt,linecnt,filename,buf);
	  fflush(stdout);
	}
	nVmHWM *= 1024LL;
	if(VERB>=2){
	  printf("VmSize=%lld:%s",nVmSize,buf);
	  fflush(stdout);
	}
	continue;// next line
      }
    }

    key = (char *)"VmRSS:";
    if((pt = strstr(buf,key))){
      pt += strlen(key);
      nVmRSS = strtoll(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt))){
	printf("value of %s on line %d of %s no recognized:\n%s\n",key,linecnt,filename,buf);
	fflush(stdout);	exit(1);
      }
      pt = qt;
      while(*pt && isspace(*pt))
	pt++;
      if(strncmp(pt,"kB",2)){
	printf("WARNING:Unknown memory size unit (%s) on line %d of %s (expected kB):\n%s\n",pt,linecnt,filename,buf);
	fflush(stdout);
      }
      nVmRSS *= 1024LL;
      if(VERB>=2){
	printf("VmRSS=%lld:%s",nVmRSS,buf);
	fflush(stdout);
      }
      continue;// next line
    }

    key = (char *)"VmSwap:";
    if((pt = strstr(buf,key))){
      pt += strlen(key);
      nVmSwap = strtoll(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt))){
	printf("value of %s on line %d of %s not recognized:\n%s\n",key,linecnt,filename,buf);
	fflush(stdout);
	exit(1);
      }
      pt = qt;
      while(*pt && isspace(*pt))
	pt++;
      if(strncmp(pt,"kB",2)){
	printf("WARNING:Unknown memory size unit (%s) on line %d of %s (expected kB):\n%s\n",pt,linecnt,filename,buf);
	fflush(stdout);
      }
      nVmSwap *= 1024LL;
      if(VERB>=2){
	printf("VmSwap=%lld:%s",nVmSwap,buf);
	fflush(stdout);
      }
      //      continue;// next line
    }
    if(nVmRSS >= 0 && nVmSize >= 0 && nVmSwap >= 0)
      break;
  }
  if(fclose(fp)){
    char *err = strerror(errno);
    printf("fclose() for /proc/self/status failed: %s\n",err);
    fflush(stdout);
  }
  if(pVmHWM){
    if(DEBUG>=2 && !(nVmHWM >= nVmRSS)){
      printf("VmRSS = %lld, VmHWM = %lld\n",nVmRSS, nVmHWM);
      fflush(stdout);
      assert(nVmHWM >= nVmRSS);
    }
  }

  if(pVmPeak) *pVmPeak = nVmPeak;
  VmSize = nVmSize;
  if(pVmHWM)
    *pVmHWM = nVmHWM;
  VmRSS = nVmRSS;
  VmSwap = nVmSwap;
#endif  
}

void getswap(long long &MemSize, long long &SwapSize, long long *AvailableMem)
{
  SwapSize = 0;

#ifndef WIN32

  FILE *fp = popen("free -m","r");

  if(fp==NULL){
    char *err = strerror(errno);
    printf("getswap(): Failed to read \"free -m\": %s\n",err);
    dumpmemmap();
    fflush(stdout); exit(1);
  }
  
  char buf[LINESIZ];
  int linecnt = 1;
  
  MemSize = SwapSize = -1;

  for(;fgets(buf,LINESIZ,fp) != NULL;linecnt++){
    size_t len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      if(len >= LINESIZ-1)
	printf("Line too long: exceeds max length=%d (increase LINESIZE in constants.h) in output of \"free -m\" on line %d:\n%s\n",LINESIZ-1,linecnt,buf);
      else
	printf("Line not correctly terminated (last char = %c) in output of \"free -m\" on line %d:\n%s\n",buf[len-1],linecnt,buf);
      fflush(stdout); exit(1);
    }

    char *pt, *qt, *key;

    key = (char *)"Mem:";
    if((pt = strstr(buf,key))){
      pt += strlen(key);

      MemSize = strtoll(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt))){
	printf("1st value after %s on line %d from output of \"free -m\" not recognized:\n%s\n",key,linecnt,buf);
	fflush(stdout);	exit(1);
      }
      MemSize *= 1024LL * 1024LL;

      if(AvailableMem != 0){
	/* skip next 5 fields */
	while(*pt && isspace(*pt))
	  pt++;
	while(*pt && !isspace(*pt))
	  pt++;
	
	while(*pt && isspace(*pt))
	  pt++;
	while(*pt && !isspace(*pt))
	  pt++;
	
	while(*pt && isspace(*pt))
	  pt++;
	while(*pt && !isspace(*pt))
	  pt++;
	
	while(*pt && isspace(*pt))
	  pt++;
	while(*pt && !isspace(*pt))
	  pt++;
	
	while(*pt && isspace(*pt))
	  pt++;
	while(*pt && !isspace(*pt))
	  pt++;

	*AvailableMem = strtoll(pt,&qt,10);
	if(pt == qt || !(*qt == 0 || isspace(*qt))){
	  printf("6th value after %s on line %d from output of \"free -m\" not recognized:\n%s\n",key,linecnt,buf);
	  fflush(stdout);	exit(1);
	}
	*AvailableMem *= 1024LL * 1024LL;
      }

      continue;// next line
    }

    key = (char *)"Swap:";
    if((pt = strstr(buf,key))){
      pt += strlen(key);

      /* skip next 2 fields */
      while(*pt && isspace(*pt))
	pt++;
      while(*pt && !isspace(*pt))
	pt++;
      
      while(*pt && isspace(*pt))
	pt++;
      while(*pt && !isspace(*pt))
	pt++;
      
      SwapSize = strtoll(pt,&qt,10);
      if(pt == qt || !(*qt == 0 || isspace(*qt))){
	printf("3rd value after %s on line %d from output of \"free -m\" not recognized:\n%s\n",key,linecnt,buf);
	fflush(stdout);	exit(1);
      }
      SwapSize *= 1024LL * 1024LL;

      continue;// next line
    }

    if(MemSize >= 0 && SwapSize >= 0)
      break;
  }

  if(SwapSize < 0){
    printf("Could not find Swap size using \"free -m\": Assuming no swap\n");
    fflush(stdout);

    SwapSize = 0;
  }

  if(pclose(fp)){
    char *err = strerror(errno);
    printf("pclose() for \"free -m\" failed: %s\n",err);
    fflush(stdout);
  }

#endif
}

const char *global_id = 0;

int num_blocks = 0;/* number of Cmap & site[] blocks allocated */
int MAX_BLOCKS = 0;
Cmap **Cmap_blocks = 0;
char **site_blocks = 0;/* memory for site[c], SNR[c],Intensity[c],Stitch[c],stitchLocation for a set of maps */
Cmap ****id_blocks = 0;/* use of each block identified by a pointer to either &map or &refmap (IF these are initial allocations, to be freed by mapcompact) */
size_t *mapcnt_blocks = 0;;/* size of Cmap_blocks[] */

char **ChipIds = 0; /* set of distinct chipid strings */
int numChipIds = 0, maxChipIds = 0;/* number of distinct chipid strings and number allocated */

void maxblockalloc(int num)
{
  if(num < MAX_BLOCKS)
    return;
  if(num < 2*MAX_BLOCKS + 1024)
    num = 2*MAX_BLOCKS + 1024;
  if(VERB>=2 || LEAKDEBUG){
    printf("MAX_BLOCKS= %d -> %d\n", MAX_BLOCKS, num);
    fflush(stdout);
  }

  #pragma omp critical(MAX_BLOCKS)
  {
    if(MAX_BLOCKS <= 0){
      Cmap_blocks = (Cmap **)malloc(num * sizeof(Cmap *));
      site_blocks = (char **)malloc(num * sizeof(char *));
      id_blocks = (Cmap ****)malloc(num * sizeof(Cmap ***));
      mapcnt_blocks = (size_t *)malloc(num * sizeof(size_t));
    } else {
      Cmap_blocks = (Cmap **)realloc(Cmap_blocks, num * sizeof(Cmap *));
      site_blocks = (char **)realloc(site_blocks, num * sizeof(char *));
      id_blocks = (Cmap ****)realloc(id_blocks, num * sizeof(Cmap ***));
      mapcnt_blocks = (size_t *)realloc(mapcnt_blocks, num * sizeof(size_t));
    }
    if(!Cmap_blocks || !site_blocks || !id_blocks || !mapcnt_blocks){
      printf("Failed to reallocate %d->%d blocks of memory\n",MAX_BLOCKS, num);
      fflush(stdout); exit(1);
    }
    MAX_BLOCKS = num;
  }
}

/* return a pointer to a copy of chipid in ChipIds[0..numChipIds-1] */
char *ChipIdAlloc(const char *chipid)
{
  char *ret = 0;

  #pragma omp critical(ChipIdslock)
  {
    /* first recheck if chipid is already in ChipIds[0..numChipIds-1] */
    for(int i = 0; i < numChipIds; i++)
      if(!strcmp(chipid,ChipIds[i])){
        ret = ChipIds[i];
        break;
      }

    if(!ret){    /* allocate a copy of chipid and store it in ChipIds[0..numChipIds-1] */
      if(numChipIds >= maxChipIds){
        maxChipIds = max(maxChipIds * 2, 64);
        ChipIds = (char **)realloc(ChipIds, sizeof(char *) * maxChipIds);
      }
      if(VERB>=2){
        printf("Adding chipid=%s as ChipIds[%d]\n",chipid,numChipIds);
	fflush(stdout);
      }
      ChipIds[numChipIds++] = ret = strdup(chipid);
    }
  }

  return ret;
}

int giter2 = 0;/* outer loop variable for -M <A> <B> */

pid_t pid = 0;/* value of getpid() */

int startmaps = 0;/* number of maps read in. If ScanCorrection > 0 : Gmap[0..startmaps-1] will be output in _rescaled.bnx (NOT used in pairalign) */
int nummaps = 0;/* number of maps passed to refalign or pairalign (may be less than startmaps if -subset etc was used) */
int totalmaps = 0;/* In refalign : totalmaps >= startmaps, including map splits (NOT used in pairalign) */
int rawsitemaps = 0;/* number of maps <= totalmaps that have rawsite[] allocated */


long long minid= 0,maxid = 0;

int maxmaps = 0;
Cmap **Gmap = 0;
int Hapnummaps = 0;/* the value of nummaps before splitting Haplotype Query maps : the second Allele(s) are placed at refmap[orignummaps .. nummaps-1] */
int HapQuerymaps = 0;/* number of split Haplotype Query maps */
long long MaxQueryID = 0;/* offset applied to Allele2 query contig id */

int idrenumbered = SEQID;/* 1 IFF map id numbers were replaced by sequential id numbers */
int maptype= -1;/* 0 == nanomap (.bnx or .vfixx), 1 == consensus (.cmap) */
int UniqueScans = 0;/* For BNX Version 1.0 input : number of unique Scans : Gmap[0..nummaps-1]->UniqueScanId ranges over 0 .. UniqueScans-1 */

int RunDataListLen = 0, RunDataListMax = 0;
char **RunDataList = 0;/* If != 0 : RunDataList[0..RunDataListLen-1] is the list of RunData header lines from all BNX Version 1.0 input files */

int numrefmaps = 0;
int maxrefmaps = 0;
Cmap **refmap = 0;
int Hapnumrefmaps = 0;/* the value of numrefmaps before splitting Haplotyed refmaps : the second Allele(s) are placed at refmap[orignumrefmaps .. numrefmaps-1] */
int Hapmaps = 0;/* number of original ref maps that were split */
long long MaxContigID=0;/* offset applied to Allele2 contig id */

size_t numaligns = 0;/* number of alignments in alignment[] */
size_t maxaligns = 0;
Calign **alignment = 0;/* array alignment[i=0..numaligns-1] of map alignments */
int alignment_blockid = -1;/* If >=0, alignment[] array is part of some Calign_block */

CalignA **alignmentA = 0;
Calign **alignmentT = 0;/* transposed version of alignmentT[] used by refalign */

/* The alignments for each refmap i = 0..numrefmaps-1 are located in alignment[numalign_start[i]..numalign_end[i]-1] : useful when there are multiple reference maps */
size_t *numalign_start = 0;
size_t *numalign_end = 0;

size_t *numalign_mstart = 0;
size_t *numalign_mend = 0;

Cmap ** &gmap = Gmap;/* Alias for (Cmap **)map in globals.h (due to name alias problem in refine()) */
int &gnummaps = nummaps; /* Alias for (int) nummaps in globals.h (due to name alias problem in Ccontig.cpp) */

double draftSD = 0.40;/* draft map smoothing SD as multiple of res */
double draftTP = 0.36;/* minimum true positive rate threshold to call draft consensus cut */

int numxmapentries = 0;/* number of active xmapEntry objects in xmapentries[0..numxmapentries-1] */
int maxxmapentries = 0;/* allocated size of xmapentries[0.. maxxmapentries-1] */
xmapEntry **xmapentries = NULL;/* array which corresponds to entries in the xmap--filled in output_xmap */
int smapsize = 0;/* number of filtered xmapEntry objects in usexmapentries[0 .. smapsize -1 ] */
xmapEntry **usexmapentries = NULL;/* array usexmapentries[0 .. smapsize-1] is a re-ordered subset of xmapentries[0..numxmapentries-1] */

int numsv = 0; //n elements of sv_array used
int maxsv = 0; //size of sv_array
structuralVariation *sv_array = 0;

int nummaskentries = 0; //size is fixed after input, so only one int necessary
maskSV *maskentries = 0; //see input_mask.txt

int numscafentries = 0;/* number of scaffold objects loaded from scaffold file */
scaffold **scaffolds = 0;/* array of (pointers to) scaffold objects */

int numbedentries = 0; /* number of elements in bedentries */
bedEntry **bedentries = 0;

int numrepeats = 0; //number of elements of repeat_array used
int maxrepeats = 0; //size of repeat_array
repeat *repeat_array = 0;

int NumScaleFactor = 0;
double *ScaleFactor = 0;/* array of Molecule size Scaling factors tried : ScaleFactor[0..NumScaleFactor-1] */

const char *VQ_TRUESITE[2][2] = {{"QX10","QX10"},   /* Q-score code for TrueSite for labels for Channel/color 1 (for simulated data) */
				 {"QX20","QX20"}};   /* Q-score code for TrueSite for labels for Channel/color 2 (for simulated data) */
const char *VQ_SNR[2][2] = {{"QX02","QX11"}, /* Q-score code of SNR signal for Channel/color 1 (For BNX version 0.1, 1.0 respectively) */
                            {"QY02","QX21"}}; /* Q-score code of SNR signal for Channel/color 2 */
const char *VQ_INTENSITY[2][2] = {{"QX12","QX12"}, /* Q-score code of Intensity signal for Channel/color 1 */
				  {"QX22","QX22"}};/* Q-score code of Intensity signal for Channel/color 2 */
const char *VQ_STITCH[2][2] = {{"QX04","QX13"},  /* Q-score code of Stitch and/or Knot signal for Channel/color 1 */
			       {"QY04","QX23"}}; /* Q-score code of Stitch and/or Knot signal for Channel/color 2 */
const char *VQ_STITCHLOC[2][2] = {{"QX05","QX14"}, /* Q-score code for Stitch location for Channel/color 1 */
				  {"QY05","QX24"}};/* Q-score code for Stitch location for Channel/color 2 */

const char *VQ_PSFWIDTH[2][2] = {{"QX15","QX15"},  /* Q-score code for PSF width for labels for Channel/color 1 */
				 {"QX25","QX25"}};  /* Q-score code for PSF width for labels for Channel/color 2 */
const char *VQ_IMAGELOC[2][2] = {{"QX16","QX16"},  /* Q-score code for Image Locations (FOV1 X1 Y1 FOV2 X2 Y2 ... ) for labels for Channel/color 1 */
                                {"QX26","QX26"}}; /* Q-score code for Image Locations (FOV1 X1 Y1 FOV2 X2 Y2 ... ) for labels for Channel/color 2 */


/* find id in sorted array p[0..nummaps-1] and return corresponding mapid OR -1 if not found */
int findid(long long id, Cid2mapid *p, int nummaps)
{
  if(VERB>=2){
    printf("findid(id=%lld,nummaps=%d):id range=%lld..%lld:",id,nummaps,p[0].id,p[nummaps-1].id);
    fflush(stdout);
  }

  int low = 0;
  int high = nummaps-1;
  while(high > low){
    int mid = (low+high)/2;
    long long int midID = p[mid].id;
    _mm_prefetch((const char *) &p[(low+mid)>>1].id, _MM_HINT_T0);
    _mm_prefetch((const char *) &p[(mid+high)>>1].id, _MM_HINT_T0);
    if(midID < id)
      low = mid+1;
    else if(midID > id)
      high = mid-1;
    else
      low = high = mid;
  }
  if(DEBUG && !(low>=high)){
    printf("\nfindid(id=%lld,nummaps=%d):id range=%lld..%lld:low=%d,high=%d\n",id,nummaps,p[0].id,p[nummaps-1].id,low,high);
    fflush(stdout);
    assert(low>=high);
  }
  if(high < low || p[low].id != id)
    return -1;
  if(VERB>=2){
    printf("mapid=%d\n",p[low].mapid);
    fflush(stdout);
  }
  return p[low].mapid;
}

/* unlike findid() it returns the index : this allows nearby id values (including duplicate id values) to be located */
int findindex(long long int id, Cid2mapid *p, int nummaps)
{
  if(VERB>=2){
    printf("findid(id=%lld,nummaps=%d):id range=%lld..%lld:",id,nummaps,p[0].id,p[nummaps-1].id);
    fflush(stdout);
  }

  int low = 0;
  int high = nummaps-1;
  while(high > low){
    int mid = (low+high)/2;
    long long int midID = p[mid].id;
    _mm_prefetch((const char *) &p[(low+mid)>>1].id, _MM_HINT_T0);
    _mm_prefetch((const char *) &p[(mid+high)>>1].id, _MM_HINT_T0);
    if(midID < id)
      low = mid+1;
    else if(midID > id)
      high = mid-1;
    else
      low = high = mid;
  }
  if(DEBUG && !(low>=high)){
    printf("\nfindid(id=%lld,nummaps=%d):id range=%lld..%lld:low=%d,high=%d\n",id,nummaps,p[0].id,p[nummaps-1].id,low,high);
    fflush(stdout);
    assert(low>=high);
  }
  if(high < low || p[low].id != id)
    return -1;
  if(VERB>=2){
    printf("mapid=%d\n",p[low].mapid);
    fflush(stdout);
  }
  return low;
}

int siteinc(register Cmerge *p1, register Cmerge *p2)
{
  return (p1->site > p2->site) ? 1 : (p1->site < p2->site) ? -1 : 0;
}

#include "parameters.h"

int checkFile(char *filename)
{
  if(!output_filename_matched(filename))
    return 1;

  FILE *fp;
  if((fp = fopen(filename,"r")) == NULL)
    return 0;

  if(ForceOverwrite){
    #pragma omp critical
    {
      int tid = -1;
#ifdef _OPENMP
      tid = omp_get_thread_num ();
#endif

      printf("WARNING:Output file %s already exists(tid=%d)\n",filename,tid);
      fflush(stdout);

#if 0 // DEBUG code
      char buf[BUFSIZ];
#if MIC
      sprintf(buf,"ssh -y -y -oStrictHostKeyChecking=no -oBatchMode=yes host ls -l %s",filename);
#else
      sprintf(buf,"ls -l %s",filename);
#endif
      int err = system(buf);
      if(err){
	printf("system(%s) failed: return code = %d\n", buf,err);
	fflush(stdout);
      }
#endif

      fclose(fp);
    
      /* unlink the filename, in case it is a hard link to another file : otherwise the other linked files get changed as well */
      int ret, retries = 3;
      errno = 0;
      for(int i = 0; i < retries && (ret = unlink(filename)); i++){
	int eno = errno;
	char *err = strerror(eno);
	if(i + 1 >= retries){
	  printf("unlink(%s) failed:errno=%d:%s (Giving up : Will try to overwrite file):tid=%d\n",filename,eno,err,tid);
	  fflush(stdout);
	  break;
	}
	printf("unlink(%s) failed:errno=%d:%s (Retry attempt %d/%d in 10 seconds):tid=%d\n",filename,eno,err,i+1, retries,tid);
	fflush(stdout);

	sleep(10);
	errno = 0;
      }

#if 0 // DEBUG code
      //      char buf[BUFSIZ];
#if MIC
      sprintf(buf,"ssh -y -y -oStrictHostKeyChecking=no -oBatchMode=yes host ls -l %s",filename);
#else
      sprintf(buf,"ls -l %s",filename);
#endif
      err = system(buf);
      if(err){
	printf("system(%s) failed: return code = %d\n", buf,err);
	fflush(stdout);
      }
#endif
    }
    return 0;
  }

  #pragma omp critical
  {

    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num ();
#endif

    printf("Output file %s already exists (use -f to overwrite):tid=%d,fp=%p\n",filename,tid,fp);
    fflush(stdout);

    char buf[BUFSIZ];
#if MIC
    sprintf(buf,"ssh -y -y -oStrictHostKeyChecking=no -oBatchMode=yes host ls -l %s",filename);
#else
    sprintf(buf,"ls -l %s",filename);
#endif
    int err = system(buf);
    if(err){
      printf("system(%s) failed: return code = %d\n", buf,err);
      fflush(stdout);
    }

    if(0){/* print out start of file */
      const size_t tlen = BUFSIZ;
      char rbuf[tlen];
      errno = 0;
      size_t len = fread(rbuf, 1, tlen - 1, fp);
      int eno = errno;
      char *err = strerror(eno);
      printf("fread(BUFSIZ=%lu)=%lu,feof()=%d,ferror()=%d,errno=%d:%s\n",tlen,len,feof(fp),ferror(fp),eno,err);
      if(len > 0){
	rbuf[len] = '\0';
	for(size_t i = 0; i < len; i++)
	  printf("%c",rbuf[i]);
	printf("\n");
      }
      fflush(stdout);
    }

    exit(1);
  }

  return 0;
}

#if USE_MMAP

#include <assert.h>
#include <stdlib.h>
#include <sys/mman.h>
#define HUGE_PAGE_SIZE (2 * 1024 * 1024)
#define ALIGN_TO_PAGE_SIZE(x) ((((x) + HUGE_PAGE_SIZE -1) / HUGE_PAGE_SIZE) * HUGE_PAGE_SIZE)

void *malloc_huge_pages(size_t size)
{
  // Use 1 extra page to store allocation metadata
  // (libhugetlbfs is more efficient in this regard)
  size_t real_size = ALIGN_TO_PAGE_SIZE(size + HUGE_PAGE_SIZE);
  //  char *ptr = (char *)mmap(NULL, real_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB, -1, 0);
  char *ptr = (char *)mmap(NULL, real_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
  if (ptr == MAP_FAILED) {  // The mmap() call failed. Try to regular mmap() then malloc instead
    ptr = (char *)malloc(real_size);
    if (ptr == NULL) return NULL;
    real_size = 0;
  }
  // Save real_size since mmunmap() requires a size parameter
  *((size_t *)ptr) = real_size;

  return ptr + HUGE_PAGE_SIZE;  // Skip the page with metadata
}

void free_huge_pages(void *ptr)
{
  if (ptr == NULL) return;
  // Jump back to the page with metadata
  void *real_ptr = (char *)ptr - HUGE_PAGE_SIZE;

  // Read the original allocation size
  size_t real_size = *((size_t *)real_ptr);
  assert(real_size % HUGE_PAGE_SIZE == 0);
  if (real_size != 0) // The memory was allocated via mmap() and must be deallocated via munmap()
    munmap(real_ptr, real_size);
  else // The memory was allocated via malloc() and must be deallocated via free()
    free(real_ptr);
}

#endif // USE_MMAP
