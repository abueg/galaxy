#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#ifndef WIN32
#include <unistd.h>
#else
#include <direct.h>
#define getcwd _getcwd
#endif

#include "constants.h"
#include "globals.h"
#include "parameters.h"
#include "Assembler_parameters.h"
#include "Ccontig.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/output_contigs.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

extern void printversion(FILE *fp);
extern double mtime(),wtime();

void output_contigs(Ccontig *contigs, int numcontigs, char *basename)
{
  if(strstr(basename,"/dev/null"))
    return;

  /* output in .contigs file format */
  char filename[PATH_MAX];
  strcpy(filename,basename);

  /* no need to strip suffix : for default basename this was already done in RefAligner.cpp */
  int i = strlen(filename);
  sprintf(&filename[i],".contigs");

  if(checkFile(filename))
    return;

  double start = mtime(), wstart = wtime();
  
  if(VERB){
    printf("Generating %s (Assembly contigs data structure)\n",filename);
    fflush(stdout);
  }

  FILE *fp;
  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("failed to open file %s for writing Assembly contigs data structure:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  char *write_buffer = NULL;
  size_t blen = contig_format ? 1024*1024 : 10*1024*1024;
  write_buffer = (char *)malloc(blen);
  setbuffer(fp, write_buffer, blen);

  FILE *fpbin = NULL;
  char filenamebin[PATH_MAX];
  char *write_bufferbin = NULL;
  size_t blen2 = 10*1024*1024;
  if(contig_format){
    sprintf(filenamebin,"%s.bin",filename);
    if((fpbin = fopen(filenamebin,"wb"))==NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("failed to open file %s for writing Assembly contigs binary data:errno=%d:%s\n", filenamebin,eno,err);
      exit(1);
    }
    write_bufferbin = (char *)malloc(blen2);
    setbuffer(fpbin, write_bufferbin, blen2);
  }

  /* write out commandline */
  printversion(fp);

  fprintf(fp,"colors=%d,numcontigs=%d\n",colors,numcontigs);
  if(colors != 1){
    printf("colors=%d: multiple colors not implemented\n",colors);
    exit(1);
  }

  size_t offset = 0;/* binary file offset at start of each contig */
  for(register int C = 0; C < numcontigs; C++){
    register Ccontig *Ycontig = &contigs[C];

    register int n = Ycontig->numsite[0];
    register double *Y = Ycontig->site[0];
    register float *sitecnt = Ycontig->sitecnt[0];
    register float *sitecntFN = Ycontig->sitecntFN[0];
    register float *fragcnt = Ycontig->fragcnt[0];

    register int MD = Ycontig->nummaps;
    register Ccontig *Xcontig = Ycontig->contig;
    register int *Xflip = Ycontig->flip;
    register int **Xsitemap = Ycontig->sitemap[0];
    register double **X = Ycontig->X[0];

    if(DEBUG) assert(Ycontig -> mapid < 0 && Ycontig->nummaps > 0);
    if(DEBUG) assert(Y[0] == 0.0);

    fprintf(fp,"Contig=%d numsite=%d,nummaps=%d\n", C, n, MD);
    if(contig_format==0){
      for(register int i = 1; i <= n; i++)
	fprintf(fp,"i=%d:site[i]=%0.6f,fragcnt[i-1]=%0.2f,sitecnt[i]=%0.2f,sitecntFN[i]=%0.2f\n",
		i,Y[i],fragcnt[i-1],sitecnt[i],sitecntFN[i]);
      fprintf(fp,"i=%d:site[i]=%0.6f,fragcnt[i-1]=%0.2f\n",
	      n+1,Y[n+1],fragcnt[n]);
    
      for(register int m = 0; m < MD; m++){
	if(DEBUG) assert(Xcontig[m].mapid >= 0);
	register int M = Xcontig[m].numsite[0];
	if(DEBUG) assert(Gmap[Xcontig[m].mapid]->mapid == Xcontig[m].mapid);/* If this fails, the input Gmap[] array be have been re-sorted since it was read in */
	fprintf(fp,"map=%d:mapid=%d,id=%lld,trimL=%d,trimR=%d,flip=%d,numsite=%d\n",
		m,Xcontig[m].mapid,Gmap[Xcontig[m].mapid]->id,Xcontig[m].trimL[0],Xcontig[m].trimR[0],Xflip[m],M);
	for(register int j = 0; j <= M+1; j++)
	  fprintf(fp,"  j=%d:sitemap=%d,X=%0.6f\n",j,Xsitemap[m][j],X[m][j]);
      }
    } else {/* use binary format : output entire contig as binary chunk followed by CONTIG_MAGIC and "\n"*/
      fprintf(fp,"binary offset=%llu\n", (unsigned long long)offset);
      if(VERB>=2){
	printf("Contig=%d numsite=%d,nummaps=%d\n",C,n,MD);
	for(int i = 1; i <= min(5,n+1); i++)
	  printf("  i=%d:site[i]=%lf,fragcnt[i-1]=%f,sitecnt[i]=%f,sitecntFN[i]=%f\n",
		 i,Y[i],fragcnt[i-1],sitecnt[i],sitecntFN[i]);
	fflush(stdout);
      }

      size_t r;
      if((r = fwrite_LOOP(&Y[1],sizeof(double), n+1, fpbin)) != (size_t)(n+1)){
	printf("fwrite_LOOP of %d doubles from &Y[1] to %s failed (ret=%llu)\n",n+1,filenamebin,(unsigned long long)r);
	exit(1);
      }
      offset += r*sizeof(double);
      if((r = fwrite_LOOP(&fragcnt[0],sizeof(float), n+1, fpbin)) != (size_t)(n+1)){
	printf("fwrite_LOOP of %d floats from &fragcnt[0] to %s failed (ret=%llu)\n",n+1,filenamebin,(unsigned long long)r);
	exit(1);
      }
      offset += r*sizeof(float);
      if((r = fwrite_LOOP(&sitecnt[1],sizeof(float), n, fpbin)) != (size_t)n){
	printf("fwrite_LOOP of %d floats from &sitecnt[1] to %s failed (ret=%llu)\n",n,filenamebin,(unsigned long long)r);
	exit(1);
      }
      offset += r*sizeof(float);
      if((r = fwrite_LOOP(&sitecntFN[1],sizeof(float), n, fpbin)) != (size_t)n){
	printf("fwrite_LOOP of %d floats from &sitecntFN[1] to %s failed (ret=%llu)\n",n,filenamebin,(unsigned long long)r);
	exit(1);
      }
      offset += r*sizeof(float);
      long long *id = new long long[MD];
      int *mapid = new int[MD*5];
      int *trimL = &mapid[MD];
      int *trimR = &mapid[2*MD];
      int *flip = &mapid[3*MD];
      int *numsite = &mapid[4*MD];
      for(int m = 0; m < MD; m++){
	id[m] = Gmap[Xcontig[m].mapid]->id;
	mapid[m] = Xcontig[m].mapid;
	trimL[m] = Xcontig[m].trimL[0];
	trimR[m] = Xcontig[m].trimR[0];
	flip[m] = Xflip[m];
	numsite[m] = Xcontig[m].numsite[0];
      }
      if((r = fwrite_LOOP(id, sizeof(long long), MD, fpbin)) != (size_t)MD){
	printf("fwrite_LOOP of %d long longs of id[] to %s failed (ret=%llu)\n",MD,filenamebin,(unsigned long long)r);
	exit(1);
      }
      offset += r*sizeof(long long);
      if((r = fwrite_LOOP(mapid, sizeof(int), MD*5, fpbin)) != (size_t)(MD*5)){
	printf("fwrite_LOOP of 5*%d ints to %s failed (ret=%llu)\n",MD,filenamebin,(unsigned long long)r);
	exit(1);
      }
      offset += r*sizeof(int);
      if(VERB>=2){
	for(int m = 0; m < min(5,MD); m++)
	  printf("  m=%d/%d:id=%lld,mapid=%d,trimL=%d,trimR=%d,flip=%d,numsite=%d\n",
		 m,MD,id[m],mapid[m],trimL[m],trimR[m],flip[m],numsite[m]);
	fflush(stdout);
      }

      for(int m = 0; m < MD; m++){
	if(VERB>=2){
	  for(int j = 0; j <= numsite[m]+1; j++)
	    printf("  m=%d/%d:j=%d/%d:sitemap[m][j]=%d(n=%d),X[m][j]=%0.3f\n",m,MD,j,numsite[m],Xsitemap[m][j],n, X[m][j]);
	  fflush(stdout);
	}
	if((r = fwrite_LOOP(&Xsitemap[m][0], sizeof(int), numsite[m]+2, fpbin)) != (size_t)(numsite[m]+2)){
	  printf("fwrite of %d ints to %s failed (ret=%llu,m=%d)\n",numsite[m],filenamebin,(unsigned long long)r,m);
	  exit(1);
	}
	offset += r*sizeof(int);
	if((r = fwrite_LOOP(&X[m][0], sizeof(double), numsite[m]+2, fpbin)) != (size_t)(numsite[m]+2)){
	  printf("fwrite of %d doubles to %s failed (ret=%llu,m=%d)\n",numsite[m]+2,filenamebin,(unsigned long long)r,m);
	  exit(1);
	}
	offset += r*sizeof(double);
      }

      delete [] mapid;
      delete [] id;
    }
  }

  if(contig_format){
    fprintf(fp,"total binary output = %llu bytes\n",(unsigned long long)offset);

    FILEclose(fpbin);
    if(write_bufferbin)
      free(write_bufferbin);
  }
  FILEclose(fp);
  if(write_buffer)
    free(write_buffer);
  if(VERB){
    printf("Finished writing %s: time=%0.6f, wall time=%0.6f secs\n",filename, mtime()-start,wtime()-wstart);
    fflush(stdout);
  }
}

