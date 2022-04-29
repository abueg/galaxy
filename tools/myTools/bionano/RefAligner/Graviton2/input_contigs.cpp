#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include <math.h>
#include <omp.h>

#include "constants.h"
#include "globals.h"
#include "parameters.h"
#include "Assembler_parameters.h"
#include "Ccontig.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/input_contigs.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

static  char buf[LINESIZ];

extern double mtime(), wtime();

#define PARTIAL_READ 1 // read only part of the binary contig file at a time to minimize memory required

typedef struct {
  char *buf;
  long start;
  long length;
  size_t binstart;
  int linecnt;
} PARSE_CHUNK;

void input_contigs(Ccontig * &contigs, int &numcontigs, char *filename)
{
  FILE *fp;
  if((fp = fopen(filename,"r"))==NULL){
    printf("failed to open file %s for reading Assembly contigs data structure\n",filename);
    exit(1);
  }
  if(VERB){
    printf("Reading %s (Assembly contigs data structure): wall time=%0.6f\n",filename, wtime());
    fflush(stdout);
  }
  FILE *fpbin = NULL;
  char filenamebin[PATH_MAX];
  size_t blen = 10*1024*1024;/* matches scif_io buffer size */
  char *fbufbin = 0;
  if(contig_format){
    sprintf(filenamebin,"%s.bin",filename);
    if((fpbin = fopen(filenamebin,"rb"))==NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("failed to open file %s for reading Assembly contigs binary data:errno=%d:%s\n",filenamebin,eno,err);
      exit(1);
    }
    fbufbin = (char *)malloc(blen);
    if(!fbufbin){
      printf("input_align:malloc(%llu) failed\n",(unsigned long long)blen);
      exit(1);
    }
    setbuffer(fpbin, fbufbin, blen);  
  }

  int linecnt = 1;
  int NColors, input;

  while(fgets(buf,LINESIZ,fp) != NULL){
    int len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      printf("ERROR:Line too long (len=%d) or not terminated (last char='%c') in input file %s on line %d (increase LINESIZE):\n%s\n",len, buf[len-1],filename,linecnt,buf);
      exit(1);
    }
    // WAS70    linecnt++;
    if(buf[0] != '#')
      break;
    linecnt++;// NEW70
  }

  if((input = sscanf(buf,"colors=%d,numcontigs=%d\n",&NColors,&numcontigs)) != 2){
    printf("ERRROR:reading line %d of %s : sscanf returned %d (expected 2):\n%s\n",linecnt,filename,input,buf);
    exit(1);
  }
  if(NColors != colors){
    if(colors_set){
      printf("ERROR:colors=%d on line %d of %s: must match -colors or -i input value or %d\n",NColors,linecnt,filename,colors);
      exit(1);
    }
    colors = NColors;
    colors_set = 1;
  }
  if(colors != 1){
    printf("ERROR:colors=%d: multiple colors not implemented\n",colors);
    exit(1);
  }
  linecnt++;

  size_t pos_cur = ftell(fp);
  fseek(fp, 0L, SEEK_END);
  size_t pos_end = ftell(fp);
  fseek(fp, pos_cur, SEEK_SET);
  
  char *buffer = (char *)malloc(pos_end-pos_cur+1);
  if(!buffer){
    printf("input_contigs: malloc(%llu) failed\n",(unsigned long long)(pos_end-pos_cur+1));
    exit(1);
  }
  size_t len = fread(buffer, 1, pos_end-pos_cur, fp);

  if(DEBUG && !(len == pos_end-pos_cur)){
    printf("pos_cur=%llu, pos_end=%llu, len=%llu\n",(unsigned long long)pos_cur,(unsigned long long)pos_end,(unsigned long long)len);
    fflush(stdout);
    assert(len == pos_end-pos_cur);
  }
  buffer[len] = 0;
  fclose(fp);
  fp = NULL;
  
  if(VERB){
    printf("%s read into memory buffer: len=%llu, contig_start=%d, contig_end=%d: time=%0.6f, wall time=%0.6f\n", filename, (unsigned long long)len, contig_start,contig_end, mtime(), wtime());
    fflush(stdout);
  }

  /* allocate contigs */
  contigs = new Ccontig[numcontigs];

  PARSE_CHUNK *parse_chunks = new PARSE_CHUNK[numcontigs];
  int contig_count= 0, chunk_linecnt= linecnt;
  parse_chunks[contig_count].buf = buffer;
  parse_chunks[contig_count].start = 0;
  parse_chunks[contig_count].linecnt = chunk_linecnt;

  /* initialize newmapid = -1 for all maps : Maps to be saved will be set to >= 0 */
  for(int i = 0; i < nummaps; i++)
    Gmap[i]->newmapid = -1;
  
  size_t i;
  for(i = 0; i < len; i++) {
    if(buffer[i]=='\n') {
      chunk_linecnt++;
      if((i+8 < len) && !strncmp(&buffer[i+1], "Contig=", 7)) {
	buffer[i] = 0;
	parse_chunks[contig_count].length = i - parse_chunks[contig_count].start;
	if(VERB>=2){
	  printf("Contig[%d]:start=%ld,length=%ld bytes in %s\n",contig_count,parse_chunks[contig_count].start,parse_chunks[contig_count].length,filename);
	  fflush(stdout);
	}
	contig_count++;

	parse_chunks[contig_count].buf = &buffer[i+1];
	parse_chunks[contig_count].start = i+1;
	parse_chunks[contig_count].linecnt = chunk_linecnt;
      }
    }
  }
  buffer[i] = 0;
  parse_chunks[contig_count].length = i - parse_chunks[contig_count].start;
  if(VERB>=2){
    printf("Contig[%d]:start=%ld,length=%ld bytes in %s\n",contig_count,parse_chunks[contig_count].start,parse_chunks[contig_count].length,filename);
    fflush(stdout);
  }
  contig_count++;
  if(contig_count != numcontigs){
    printf("ERROR:Found %d contigs in %s: expected %d contigs\n",contig_count,filename,numcontigs);
    exit(1);
  }
  
  int contig_max = min(contig_end+1,numcontigs);

  int numthreads = 1;
  #ifdef _OPENMP
  numthreads = omp_get_max_threads();
  numthreads = min(numthreads,MaxThreads);
  numthreads = max(1, min(contig_max-contig_start,numthreads));
  if(DEBUG) assert(numthreads >= 1);
  #endif

  if(VERB){
    printf("Parse chunks complete. Parsing contigs %d .. %d (total=%d) using %d threads(contig_format=%d): time=%0.6f, wall time=%0.6f secs\n", 
	   contig_start+1,contig_max, numcontigs, numthreads, contig_format,mtime(), wtime());
    if(VERB>=2){
      printf("Before copying content of contig[]: contig_start=%d,contig_end=%d:\n",contig_start,contig_end);
      for(int i = 0; i < numcontigs; i++){
	Ccontig *p = &contigs[i];
	if((max(0,contig_start-1) <= i && i <= contig_max+1) || p->nummaps > 0)
	  printf("contig[%d]=%p: mapid=%d, nummaps=%d, contig= %p\n",i, p, p->mapid, p->nummaps, p->contig);
      }
      fflush(stdout);
    }

    fflush(stdout);
  }

  int Cmax = min(contig_max, numcontigs - 1);// NOTE : need to read in binary byte range for 1 extra contig

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1 && contig_format)
  {

    #pragma omp for schedule(dynamic,1)
    for(int C = contig_start; C <= Cmax; C++){
      int linecnt = parse_chunks[C].linecnt;
      char *buf = parse_chunks[C].buf;
      long length = parse_chunks[C].length;
      long pos = 0;
      int input;
      int tmpC;
      Ccontig *Ycontig = &contigs[C];
      Ycontig->mapid = -1;
    
      if((input = sscanf(&buf[pos],"Contig=%d numsite=%d,nummaps=%d\n", &tmpC, &Ycontig->numsite[0], &Ycontig->nummaps))==EOF || input != 3){
	printf("ERROR:reading line %d of %s : sscanf returned %d (expected 3)\n",linecnt,filename,input);
	exit(1);
      }
      if(tmpC != C){
	printf("ERROR:Expected Contig=%d (got Contig=%d) on line %d of %s\n",C,tmpC,linecnt,filename);
	exit(1);
      }
      if(Ycontig->numsite[0] < 0 || Ycontig->nummaps <= 0){
	printf("ERROR:Invalid values for numsite=%d,nummaps=%d on line %d of %s\n",Ycontig->numsite[0],Ycontig->nummaps,linecnt,filename);
	exit(1);
      }
      while(buf[pos] && buf[pos] != '\n')
	pos++;
      pos++;
      if(VERB>=2){
	printf("Contig=%d/%d:numsite=%d,nummaps=%d on line %d of %s\n",C,contig_max,Ycontig->numsite[0],Ycontig->nummaps,linecnt,filename);
	fflush(stdout);
      }
      linecnt++;

      /* just read in binary offset */
      size_t offset;
      if((input = sscanf(&buf[pos], "binary offset=%llu\n", (unsigned long long *)&offset)) == EOF || input != 1){
	printf("ERROR: reading line %d of %s : scanf returned %d (expected 1)\n", linecnt, filename, input);
	exit(1);
      }
	
      while(buf[pos] && buf[pos] != '\n')
	pos++;
      //      pos++; /* NOTE : The final \n was not including in the parse_chunks[C] */
      if(VERB>=2){
	printf("Contig=%d/%d:binary offset = %llu on line %d of %s\n",C,contig_max,(unsigned long long)offset,linecnt,filename);
	fflush(stdout);
      }
      linecnt++;
      if(DEBUG) assert(pos <= length);
      
      parse_chunks[C].binstart = offset;
    }
  }

  if(Cmax >= contig_max){ /* NEW70 : reset nummaps & numsite[0] for last contig (that will NOT be read in) */
    Ccontig *lastYcontig = &contigs[Cmax];
    lastYcontig->nummaps = 0;
    lastYcontig->numsite[0] = 0;
  }

  if(contig_format){
    pos_cur = ftell(fpbin);
    fseek(fpbin, 0L, SEEK_END);
    pos_end = ftell(fpbin);
    fseek(fpbin, pos_cur, SEEK_SET);
  }

#if !PARTIAL_READ // original complete file read code
  char *bufferbin = 0;
  size_t lenbin = 0;

  if(contig_format){
    if(!(bufferbin = (char *)malloc(pos_end-pos_cur+1))){
      printf("input_contigs: malloc(%llu) failed\n", (unsigned long long)(pos_end-pos_cur+1));
      exit(1);
    }
    lenbin = fread(bufferbin, 1, pos_end - pos_cur, fpbin);
    if(lenbin != pos_end - pos_cur){
      printf("input_contigs:fread of %llu bytes from %s failed (return value= %llu)\n", (unsigned long long)(pos_end - pos_cur), filenamebin, (unsigned long long)lenbin);
      exit(1);
    }
    bufferbin[lenbin] = 0;
    fclose(fpbin);
    fpbin = NULL;

    if(fbufbin){
      free(fbufbin);
      fbufbin = 0;
    }
  }
  if(VERB){
    printf("%s read into binary memory buffer: lenbin=%llu, contig_start=%d, contig_end=%d: time=%0.6f, wall time=%0.6f\n", filenamebin, (unsigned long long)lenbin, contig_start,contig_end, mtime(), wtime());
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }
#endif

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {

#if PARTIAL_READ
   char *bufferbin = 0;
   size_t lenbin = 0;

   int tid = 0;
   #ifdef _OPENMP
   tid = omp_get_thread_num ();
   #endif
#endif

    #pragma omp for schedule(dynamic,1)
    for(int C = contig_start; C < contig_max; C++){
      int linecnt = parse_chunks[C].linecnt;
      char *buf = parse_chunks[C].buf;
      //      long length = parse_chunks[C].length;
      long pos = 0;
      int input;
      int tmpC;
      Ccontig *Ycontig = &contigs[C];
      Ycontig->mapid = -1;
    
      if(!contig_format){
	if((input = sscanf(&buf[pos],"Contig=%d numsite=%d,nummaps=%d\n", &tmpC, &Ycontig->numsite[0], &Ycontig->nummaps))==EOF || input != 3){
	  printf("ERROR:reading line %d of %s : sscanf returned %d (expected 3)\n",linecnt,filename,input);
	  exit(1);
	}
	if(tmpC != C){
	  printf("ERROR:Expected Contig=%d (got Contig=%d) on line %d of %s\n",C,tmpC,linecnt,filename);
	  exit(1);
	}
	if(Ycontig->numsite[0] < 0 || Ycontig->nummaps <= 0){
	  printf("ERROR:Invalid values for numsite=%d,nummaps=%d on line %d of %s\n",Ycontig->numsite[0],Ycontig->nummaps,linecnt,filename);
	  exit(1);
	}
	while(buf[pos] && buf[pos] != '\n')
	  pos++;
	pos++;
	if(VERB>=2){
	  printf("Contig=%d/%d:numsite=%d,nummaps=%d on line %d of %s\n",C,contig_max,Ycontig->numsite[0],Ycontig->nummaps,linecnt,filename);
	  fflush(stdout);
	}
	linecnt++;
      }

      int n = Ycontig->numsite[0];
      int MD = Ycontig->nummaps;

      double *Y = Ycontig->site[0] = new double[n+2];

      float *sitecnt = Ycontig->sitecnt[0] = new float[n+1];
      float *sitecntFN = Ycontig->sitecntFN[0] = new float[n+2];
      float *fragcnt = Ycontig->fragcnt[0] = new float[n+2];
      Ycontig->fragcntT[0] = new float[n+2];
      if(CovNorm)
	Ycontig->fragcntTnorm[0] = new float[n+2];
      if(DEBUG) assert(TrimNorm < 0);

      Ccontig *Xcontig = Ycontig->contig = new Ccontig[MD];
      int *Xflip = Ycontig->flip = new int[MD];
      int **Xsitemap = Ycontig->sitemap[0] = new int*[MD];
      double **X = Ycontig->X[0] = new double*[MD];

      int tmp;
    
      Y[0] = 0.0;

      if(contig_format){/* new binary format */
	size_t offset = parse_chunks[C].binstart;

#if PARTIAL_READ // NEW partial file read code
	size_t offset2 = (C == numcontigs-1) ? pos_end : parse_chunks[C+1].binstart;
	if(VERB>=2){
	  printf("C=%d,numcontigs=%d:offset= %lu, offset2= %lu, pos_end= %lu\n", C,numcontigs,offset,offset2,pos_end);
	  fflush(stdout);
	}
	if(bufferbin) free(bufferbin);
	lenbin = offset2 - offset;
	if(!(bufferbin = (char *)malloc(lenbin+1))){
	  printf("input_contigs: malloc(%lu) failed\n", lenbin+1);
	  exit(1);
	}

        size_t rlen = 0;

        #pragma omp critical(fpbin)
	{
	  fseek(fpbin, offset, SEEK_SET);
	  rlen = fread(bufferbin, 1, lenbin, fpbin);
	  if(rlen != lenbin){
	    printf("input_contigs:fread of %lu bytes from %s failed (return value= %lu, offset=%lu,tid=%d)\n", lenbin, filenamebin, rlen, offset,tid);
	    exit(1);
	  }
	}

	if(VERB){
          printf("C=%d:read %lu bytes from %s into binary memory buffer: offset=%lu,lenbin=%lu,end=%lu,tid=%d: time=%0.6f, wall time=%0.6f\n", C, rlen, filenamebin, offset, lenbin, pos_end, tid,mtime(), wtime());
	  fflush(stdout);
	}

	offset = 0;
	bufferbin[lenbin] = 0;
#endif

	if(DEBUG>=2) assert(offset + sizeof(double)*(n+1) <= lenbin);
	memcpy(&Y[1], &bufferbin[offset], sizeof(double)*(n+1)); 
	offset += sizeof(double)*(n+1);

	if(DEBUG>=2) assert(offset + sizeof(float)*(n+1) <= lenbin);
	memcpy(&fragcnt[0], &bufferbin[offset], sizeof(float)*(n+1)); 
	offset += sizeof(float)*(n+1);
	
	if(DEBUG>=2) assert(offset + sizeof(float)*n <= lenbin);
	memcpy(&sitecnt[1], &bufferbin[offset], sizeof(float)*n);
	offset += sizeof(float)*n;
	
	if(DEBUG>=2) assert(offset + sizeof(float)*n <= lenbin);
	memcpy(&sitecntFN[1], &bufferbin[offset], sizeof(float)*n);
	offset += sizeof(float)*n;

	long long *id = (long long *)&bufferbin[offset];
	offset += sizeof(long long) * MD;
	if(DEBUG) assert(offset <= lenbin);

	if(VERB>=2){
	  for(int i = 1; i <= n; i++)
	    printf("i=%d:site[i]=%lf,fragcnt[i-1]=%f,sitecnt[i]=%f,sitecntFN[i]=%f\n",i,Y[i],fragcnt[i-1],sitecnt[i],sitecntFN[i]);
	  printf("i=%d:site[i]=%lf,fragcnt[i-1]=%f\n",n+1,Y[n+1],fragcnt[n]);
	  fflush(stdout);
	}

	int *mapid = (int *)&bufferbin[offset];
	offset += sizeof(int)*MD*5;

	int *trimL = &mapid[MD];
	int *trimR = &mapid[2*MD];
	int *flip = &mapid[3*MD];
	int *numsite = &mapid[4*MD];
	int *CumSiteCnt = new int[MD+1];
	int TotSiteCnt = CumSiteCnt[0] = 0;

	if(VERB>=2){
	  for(int m = 0; m < MD; m++)
	    printf("m=%d/%d:id=%lld,mapid=%d,trimL=%d,trimR=%d,flip=%d,numsite=%d\n",
		   m,MD,id[m],mapid[m],trimL[m],trimR[m],flip[m],numsite[m]);
	  fflush(stdout); 
	}

	for(int m = 0; m < MD; m++){
	  memcpy(&Xcontig[m].id, &id[m], sizeof(long long));// HERE HERE HERE : id may NOT be aligned to 8 byte boundary in output_contig.cpp (fix by moving id array to earlier in binary file)
	  Xcontig[m].mapid = mapid[m];
	  if(!AVOID_BNX){
	    if(mapid[m] < 0 || mapid[m] >= nummaps){
	      printf("ERROR:Contig=%d/%d:m=%d/%d:mapid[m]=%d in %s is out of range (nummaps=%d)\n",C,contig_max,m,MD,mapid[m],filename,nummaps);
	      
	      exit(1);
	    }
	    if(Xcontig[m].id != Gmap[mapid[m]]->id){
	      printf("ERROR:Contig=%d/%d:m=%d/%d:Expected id=%lld (got id=%lld) on line %d of %s for mapid=%d (Molecule file or -minsites -maxsites -minlen -maxlen -mres options were different during original Assembly ?)\n",
		     C,contig_max,m, MD, Gmap[Xcontig[m].mapid]->id, Xcontig[m].id, linecnt,filename,Xcontig[m].mapid);
	      exit(1);
	    }
	    if(Gmap[mapid[m]]->newmapid < 0)
	      Gmap[mapid[m]]->newmapid = 0;/* mark map as needed */
	  } else {
	    if(mapid[m] < 0 || (DEBUG && mapid[m] > 100000000)){
	      printf("ERROR:Contig=%d/%d:m=%d/%d:mapid[m]=%d on line %d of %s is out of range (nummaps=%d)\n",C,contig_max,m,MD,mapid[m],linecnt,filename,nummaps);
	      exit(1);
	    }
	    if(Xcontig[m].id < 0 || (0 && DEBUG && id[m] > 2500LL * 1000LL*1000LL*1000LL*1000LL*1000LL)){
	      printf("ERROR:Contig=%d/%d:m=%d/%d:id[m]=%lld on line %d of %s is out of range for simulated data\n",C,contig_max,m,MD,Xcontig[m].id,linecnt,filename);
	      exit(1);
	    }
	  }

	  Xcontig[m].trimL[0] = trimL[m];
	  Xcontig[m].trimR[0] = trimR[m];;
	  Xflip[m] = flip[m];
	  if(numsite[m] <= 0 || (DEBUG && numsite[m] > 10000)){
	    printf("ERROR:m=%d/%d:numsite[m]=%d on line %d of %s is out of range\n",m,MD,numsite[m],linecnt,filename);
	    exit(1);
	  }
	  CumSiteCnt[m+1] = TotSiteCnt += Xcontig[m].numsite[0] = numsite[m];
	} 
	Ycontig->blockmem = 1;
	int *SiteMem = new int[TotSiteCnt + MD*3];
	double *Xmem = new double[TotSiteCnt + MD*2];
	for(int m = 0; m < MD; m++){// this loop can be multithreaded (2nd level)
	  Xsitemap[m] = &SiteMem[1 + 3*m + CumSiteCnt[m]];
	  Xsitemap[m][-1] = -1;
	  if(DEBUG) assert(&Xsitemap[m][numsite[m]+1] <= &SiteMem[TotSiteCnt + MD*3-1]);
	  size_t SitePos = offset + (CumSiteCnt[m]+2*m)*(sizeof(int) + sizeof(double));
	  memcpy(&Xsitemap[m][0], &bufferbin[SitePos], sizeof(int) * (numsite[m]+2));

	  X[m] = &Xmem[CumSiteCnt[m]+2*m];
	  if(DEBUG) assert(&X[m][numsite[m]+1] <= &Xmem[TotSiteCnt + MD*2-1]);
	  size_t XPos = SitePos + (numsite[m]+2)*sizeof(int);
	  memcpy(&X[m][0], &bufferbin[XPos], sizeof(double) * (numsite[m]+2));
	  if(DEBUG) assert(XPos + numsite[m]*sizeof(double) <= lenbin);

	  if(VERB>=2){
	    for(int j = 0; j <= numsite[m]+1; j++)
	      printf("  m=%d/%d:j=%d/%d:sitemap[m][j]=%d(n=%d),X[m][j]=%0.3f\n",m,MD,j, numsite[m], Xsitemap[m][j],n,X[m][j]);
	    fflush(stdout);
	  }
	}

	delete [] CumSiteCnt;

      } else {/* original text format */

	int next = -1;

	for(register int i = 1; i <= n; i++){
	  if((input = sscanf(&buf[pos],"i=%d:site[i]=%lf,fragcnt[i-1]=%f,sitecnt[i]=%f,sitecntFN[i]=%f\n%n",
			     &tmp,&Y[i],&fragcnt[i-1],&sitecnt[i],&sitecntFN[i],&next)) == EOF || input != 5){
	    printf("ERROR: reading line %d of %s : fscanf returned %d (expected 5)\n",linecnt,filename,input);
	    exit(1);
	  }
	  if(tmp != i){
	    printf("ERROR: expected i=%d (got i=%d) on line %d of %s\n",i,tmp,linecnt,filename);
	    exit(1);
	  }
	  pos += next;
	  linecnt++;
	}
	if((input = sscanf(&buf[pos],"i=%d:site[i]=%lf,fragcnt[i-1]=%f\n%n", &tmp,&Y[n+1],&fragcnt[n],&next)) == EOF || input != 3){
	  printf("ERROR: reading line %d of %s : fscanf returned %d (expected 3)\n",linecnt,filename,input);
	  exit(1);
	}
	if(tmp != n+1){
	  printf("ERROR: expected i=%d (got i=%d) on line %d of %s\n",n+1,tmp,linecnt,filename);
	  exit(1);
	}
	pos += next;
	linecnt++;
    
	for(register int m = 0; m < MD; m++){
	  int tmp;
	  if((input = sscanf(&buf[pos],"map=%d:mapid=%d,id=%lld,trimL=%d,trimR=%d,flip=%d,numsite=%d\n%n",
			     &tmp,&Xcontig[m].mapid,&Xcontig[m].id,&Xcontig[m].trimL[0],&Xcontig[m].trimR[0],&Xflip[m],&Xcontig[m].numsite[0],&next))==EOF || input != 7){
	    printf("ERROR: reading line %d of %s : fscanf returned %d (expected 7)\n",linecnt,filename,input);
	    exit(1);
	  }
	  if(tmp != m){
	    printf("ERROR: Expected m=%d (got m=%d) on line %d of %s\n", m, tmp, linecnt, filename);
	    exit(1);
	  }
      
	  int mapid = Xcontig[m].mapid;
	  if(!AVOID_BNX){
	    if(mapid < 0 || mapid >= nummaps){
	      printf("ERROR: mapid=%d on line %d of %s is out of range (nummaps=%d)\n",mapid,linecnt,filename,nummaps);
	      exit(1);
	    }
	  
	    if(Xcontig[m].id != Gmap[mapid]->id){
	      printf("ERROR: Expected id=%lld (got id=%lld) on line %d of %s for mapid=%d (Molecule file or -minsites -maxsites -minlen -maxlen -mres options were different during original Assembly ?)\n",
		     Gmap[Xcontig[m].mapid]->id, Xcontig[m].id, linecnt,filename,Xcontig[m].mapid);
	      exit(1);
	    }
	    if(Gmap[mapid]->newmapid < 0)
	      Gmap[mapid]->newmapid = 0;/* mark map as needed */
	  }

	  pos += next;
	  linecnt++;

	  register int M = Xcontig[m].numsite[0];
	  Xsitemap[m] = new int[M+3];
	  *Xsitemap[m]++ = -1;
	  X[m] = new double[M+2];

	  for(register int j = 0; j <= M+1; j++){
	    if((input = sscanf(&buf[pos],"  j=%d:sitemap=%d,X=%lf\n%n",&tmp,&Xsitemap[m][j],&X[m][j],&next)) == EOF || input != 3){
	      printf("ERROR: reading line %d of %s : fscanf returned %d (expected 3)\n",linecnt,filename,input);
	      exit(1);
	    }
	    if(tmp != j){
	      printf("ERROR: Expected j=%d (got j=%d) on line %d of %s\n", j, tmp, linecnt, filename);
	      exit(1);
	    }
	    if(Xsitemap[m][j] < -1 || Xsitemap[m][j] > n+1){
	      printf("ERROR: out of range sitemap=%d on line %d of %s (m=%d,j=%d,n=%d)\n",Xsitemap[m][j],linecnt,filename,m,j,n);
	      exit(1);
	    }
	    pos+=next;
	    linecnt++;
	  }

	  Xsitemap[m][0] = -1;// Can delete this fix after debugging 
	}
      } // contig_format==0
    } // for(int C = contig_start; C < contig_max; C++)
#if PARTIAL_READ
    if(bufferbin){
      free(bufferbin);
      bufferbin = 0;
    }
#endif
  } // omp parallel

  if(VERB){
    printf("Parsing contigs complete: time=%0.6f, wall time=%0.6f secs\n", mtime(), wtime());
    fflush(stdout);
  }
  delete [] parse_chunks;
  if(buffer){
    free(buffer);
    buffer = 0;
  }
  if(contig_format){
#if PARTIAL_READ
    fclose(fpbin);
    fpbin = NULL;
    if(fbufbin){
      free(fbufbin);
      fbufbin = 0;
    }
#else
    if(bufferbin){
      free(bufferbin);
      bufferbin = 0;
    }
#endif
  }

  if(!AVOID_BNX){
    /* renumber mapids in Gmap[] and contig[contig_start..contig_max-1] */
    int newmapid = 0;
    for(int i = 0; i < nummaps; i++)
      if(Gmap[i]->newmapid >= 0)
	Gmap[i]->newmapid = newmapid++;
    for(int C = contig_start; C < contig_max; C++){
      register Ccontig *Ycontig = &contigs[C];
      int MD = Ycontig->nummaps;
      register Ccontig *Xcontig = Ycontig->contig;
      for(register int m = 0; m < MD; m++){
	int mapid = Xcontig[m].mapid;
	Xcontig[m].mapid = Gmap[mapid]->newmapid;
      }
    }

    int j = 0;
    for(int i = 0; i < nummaps; i++){
      register Cmap *pmap = Gmap[i];
      if(pmap->newmapid < 0)
	continue;
      if(j < i){
	Cmap *tmp = Gmap[j];
	if(DEBUG) assert(pmap->newmapid == j);
	pmap->mapid = j;
	Gmap[j] = pmap;
	tmp->mapid = i;
	Gmap[i] = tmp;
      }
      j++;
    }
    if(VERB && j != nummaps){
      printf("Reduced number of input maps from %d to %d due to contig range %d .. %d: wall time=%0.6f secs\n",
	     nummaps,j, contig_start, contig_max-1, wtime());
      fflush(stdout);
    }
    nummaps = j;
  }
}

