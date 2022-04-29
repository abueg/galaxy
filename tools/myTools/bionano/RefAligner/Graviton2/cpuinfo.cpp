#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <math.h>
#include <ctype.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <fcntl.h>
#include <malloc.h>

#ifndef WIN32
#include <unistd.h>
#endif

#include <sys/stat.h>

using namespace std;

#include "globals.h"
#include "parameters.h"
#include "constants.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/trunk/cpuinfo.cpp 11489 2020-08-14 23:07:00Z tanantharaman $");

/* return number of avx512 units and/or AMD Zen2 info from /proc/cpuinfo */
static void cpuinfo(int &avx512, int &znver2)
{
  avx512 = 0;
  znver2 = 0;

#ifndef WIN32
#if USE_MIC==0
  char filename[BUFSIZ];
  sprintf(filename,"/proc/cpuinfo");
  
  errno = 0;
  FILE *fp = fopen(filename,"r");
  if(fp == NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("cpuinfo(): Failed to read %s: errno=%d:%s\n",filename,eno,err);
    fflush(stdout);exit(1);
  }
  
  char buf[LINESIZ];
  int linecnt = 1;

  for(;fgets(buf,LINESIZ,fp) != NULL; linecnt++){
    size_t len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      if(len >= LINESIZ-1)
	printf("Line too long: exceeds max length=%d (increase LINESIZE in constants.h) in output of \"free -m\" on line %d:\n%s\n",LINESIZ-1,linecnt,buf);
      else
	printf("Line not correctly terminated (last char = %c) in output of \"free -m\" on line %d:\n%s\n",buf[len-1],linecnt,buf);
      fflush(stdout); exit(1);
    }

    if(strstr(buf,"model name")){
      if(strstr(buf,"Xeon")){
	if(strstr(buf,"Silver") || strstr(buf,"Gold 5"))
	  avx512 = 1;
	else if(strstr(buf,"Gold 6") || strstr(buf,"Platinum"))
	  avx512 = 2;
      } else if(strstr(buf,"AMD EPYC 7"))
	znver2 = 1;

      printf("/proc/cpuinfo: avx512= %d, znver2= %d:\n%s", avx512,znver2, buf);
      fflush(stdout);

      break;
    }
  }

  fclose(fp);
#endif // !USE_MIC
#endif // !WIN32

  return;
}

/* check if current cpu supports current binary : If not try to switch to alternate binary, by using execvp().
   SSE3 means current binary requires sse3 instructions (same as __SSE3__ in other source files)
   AVX means current binary requires avx instructions (same as __AVX__ in other source files)
   AVX2 means current binary requires avx2 instructions (same as __AVX2__ in other source files)
   AVX512 means current binary requires avx512f + avx512bw + avx512cd + avx512dq instructions (same as __AVX512F__ && AVX512BW__ && AVX512CD__ && AVX512DE__ in other source files)
 */
void check_cpu(int argc, char**argv, int SSE3, int AVX, int AVX2, int AVX512)
{
#if !USE_MIC
#ifndef __ARM_ARCH
  if(VERB){
    printf("Current binary %s compiled with SSE3=%d, AVX=%d, AVX2=%d, AVX512=%d\n",argv[0],SSE3,AVX,AVX2,AVX512);
    fflush(stdout);
  }

  int avx512 = 0, znver2 = 0, Refine = 0;
  cpuinfo(avx512, znver2);

  if(DEBUG>=1+RELEASE && znver2) assert(avx512 == 0);
  if(znver2)
    avx512 = 0;

  /* check if -refine 2 or -refine 3 is being used */
  for(int i = 1; i < argc - 1; i++){
    if(!strcmp(argv[i],"-refine")){
      char *qt;
      Refine = strtol(argv[i+1],&qt,10);
      if(qt == argv[i+1] || !(*qt==0 || isspace(*qt))){
	printf("%s value  of %s is invalid\n",argv[i], argv[i+1]);
	fflush(stdout);exit(1);
      }
    }
  }

  __builtin_cpu_init();// needed to enable __builtin_cpu_supports() function to work

  char cwd[BUFSIZ];
  gethostname(cwd,BUFSIZ);

  int mismatch = 0;/* 1 IFF the current binary cannot run on current hardware OR a better binary may exist for current hardware */
  int revert = 0;/* 1 IFF it is OK to continue with current binary, if alternate binary does not exist */

#if __GNUC__ >= 8
  if(AVX512){
    printf("__builtin_cpu_supports(avx512f)= %d\n", __builtin_cpu_supports("avx512f"));
    printf("__builtin_cpu_supports(avx512bw)= %d\n", __builtin_cpu_supports("avx512bw"));
    printf("__builtin_cpu_supports(avx512cd)= %d\n", __builtin_cpu_supports("avx512cd"));
    printf("__builtin_cpu_supports(avx512dq)= %d\n", __builtin_cpu_supports("avx512dq"));
    printf("__builtin_cpu_supports(avx512vl)= %d\n", __builtin_cpu_supports("avx512vl"));
#if __GNUC__ >= 9
    printf("__builtin_cpu_is(znver2)= %d\n", __builtin_cpu_is("znver2"));
#endif
    fflush(stdout);

    if( avx512 <= 1 || Refine <= 1 || ! ( __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512cd") && __builtin_cpu_supports("avx512dq") /* && __builtin_cpu_supports("avx512vl")*/)){
      mismatch = 1;
      if ( __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512cd") && __builtin_cpu_supports("avx512dq") )
	revert = 1;
      if ( avx512 >= 2)
	printf("WARNING: host %s has dual avx512 FP units: switching to Skylake avx2 binary anyway due to -refine %d\n",cwd,Refine);
      else if ( avx512 == 1)
	printf("WARNING: host %s does not have dual avx512 FP units: switching to Skylake avx2 binary\n",cwd);
      else
	printf("WARNING: host %s does not support Intel avx512f, avx512bw, avx512cd, avx512dq instruction sets\n",cwd);
      fflush(stdout);
    }
  } else
#endif // __GNUC__ >= 8
  if(AVX2){
    printf("__builtin_cpu_supports(avx2)= %d, GNUC= %d\n", __builtin_cpu_supports("avx2"), __GNUC__);
#if __GNUC__ >= 9
    printf("__builtin_cpu_is(znver2)= %d\n", __builtin_cpu_is("znver2"));
#endif
    fflush(stdout);

    if( ! __builtin_cpu_supports("avx2") || (znver2 != (strstr(argv[0],"/amd2/")?1:0)) || ((avx512 >= 1) != (strstr(argv[0],"/skylake_avx2/") ? 1 : 0)) || (!znver2 && !avx512 && !strstr(argv[0],"/avx2/"))){
      mismatch = 1;
      if( __builtin_cpu_supports("avx2"))
	revert = 1;
      if ( znver2 && !strstr(argv[0],"/amd2/"))
	printf("WARNING: host %s supports AMD Zen2 : switching to better binary\n",cwd);
      else if ( !znver2 && strstr(argv[0],"/amd2/"))
	printf("WARNING: host %s does not support AMD Zen2 : switching to better bineary\n",cwd);
      else if( avx512 && !strstr(argv[0],"/skylake_avx2/"))
	printf("WARNING: host %s supports Intel Skylake avx2 instruction set: switching to better binary\n",cwd);      
      else if( !avx512 && strstr(argv[0],"/skylake_avx2/"))
	printf("WARNING: host %s does not support Intel Skylake avx2 instruction set: switching to better binary\n",cwd);      
      else if( !znver2 && !avx512 && !strstr(argv[0],"/avx2/"))
	printf("WARNING: host %s does not support AMD Zen2 or Intel Skylake avx2: switching to Haswell avx2 binary\n",cwd);
      else
	printf("WARNING: host %s does not support Intel avx2 instruction set\n",cwd);
      fflush(stdout);
    }
  } else if(AVX){
    printf("__builtin_cpu_supports(avx)= %d\n", __builtin_cpu_supports("avx"));
    fflush(stdout);
    if( ! __builtin_cpu_supports("avx")){
      mismatch = 1;
      printf("WARNING: host %s does not support Intel avx instruction set\n",cwd);
      fflush(stdout);
    }
  } else if(SSE3){
    printf("__builtin_cpu_supports(sse3)= %d\n", __builtin_cpu_supports("sse3"));
    fflush(stdout);

    if( ! __builtin_cpu_supports("sse3")){
      mismatch = 1;
      printf("WARNING: host %s does not support Intel sse3 instruction set\n",cwd);      
      fflush(stdout);
    }
  }

  if(mismatch){
    /* see if we can locate the correct binary */
    int found = 0;
    char execname[BUFSIZ], *p;
    if(argv[0][0] == '/')
      (void)strcpy(execname,argv[0]);
    else {
      getcwd(cwd,BUFSIZ);
      (void)sprintf(execname,"%s/%s",cwd,argv[0]);
    }
    p = strrchr(execname,'/');
    char *q = strdup(p);

    /* check if execname already incorporates /avx2/ or /skylake_avx2/ or /amd2/ or /avx/ : if so remove that from that path name */
    *p = '\0';
    char *r;
    if((r = strstr(execname,"/avx2")) && strlen(r) == strlen("/avx2"))
      p = r;
    else if((r = strstr(execname,"/skylake_avx2")) && strlen(r) == strlen("/skylake_avx2"))
      p = r;
    else if((r = strstr(execname,"/amd2")) && strlen(r) == strlen("/amd2"))
      p = r;
    else if((r = strstr(execname,"/avx")) && strlen(r) == strlen("/avx"))
      p = r;

    if(VERB/* HERE >=2 */){
      printf("execname=%s,q=%s,p=%s\n",execname,q,p);
      fflush(stdout);
    }

    if( __builtin_cpu_supports("avx2")){
      if( avx512 >= 1){
	(void)sprintf(p,"/skylake_avx2/%s",q+1);
	found = 1;
      } else if( znver2 ) {
	(void)sprintf(p,"/amd2/%s",q+1);
	found = 1;
      } else {
	(void)sprintf(p,"/avx2/%s",q+1);
	found = 1;
      }
    } else if( __builtin_cpu_supports("avx")){
      (void)sprintf(p,"/avx/%s",q+1);
      found = 1;
    } else if( __builtin_cpu_supports("sse3")){
      (void)sprintf(p,"/sse/%s",q+1);
      found = 1;
    }

    if(found && strcmp(execname,argv[0])){
      printf("Trying alternate binary: %s\n",execname);
      fflush(stdout);

      char *argv0 = argv[0];

      argv[0] = strdup(execname);
      errno = 0;
      int e = execvp(argv[0], argv);
      int eno = errno;
      char *err = strerror(eno);	
      printf("ERROR: exevp(%s,...) returned %d: errno=%d : %s\n",argv[0],e,eno,err);
      fflush(stdout);
      if(!revert)
	exit(1);

      printf("Using original binary %s\n", argv0);
      fflush(stdout);

      free(argv[0]);
      argv[0] = argv0;
    } else {
      if(!revert){
	printf("Unsupported CPU architecture (must support avx512, avx2, avx or sse3)\n");
	fflush(stdout);	exit(1);
      }
      if(found)
	printf("WARNING: build misconfiguration, proposed alternate binary= %s is same as original!\n", execname);
      printf("Using original binary %s\n", argv[0]);
      fflush(stdout);
    }
  }

#endif // !USE_MIC
#endif // !__ARM_ARCH
}
