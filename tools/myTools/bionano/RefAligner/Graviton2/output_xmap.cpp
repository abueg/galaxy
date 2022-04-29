#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include <math.h>
#ifndef WIN32
#include <unistd.h>
#else
#include <direct.h>
#define getcwd _getcwd
#endif

#include "constants.h"
#include "globals.h"
#include "parameters.h"
#include "Calign.h"
#include "structuralVariation.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/output_xmap.cpp 11194 2020-06-18 00:25:28Z tanantharaman $");

//#undef DEBUG
//#define DEBUG 2

#define INDEL_TRACE 0
#define REF_TRACE 21LL
#define MAP_TRACE 291LL

extern double zscore(double pvalue);

extern Cmap **YYmap,**XXmap;
extern int numY,numX;

class Cindel {
public:
  long long Xmap_id;
  long long Ymap_id;
  double qrystartpos;
  double qryendpos;
  double refstartpos;
  double refendpos;
  double confidence;
  int orientation;
  int c;/* color */
  size_t aligncnt;
};

static int IndelRefInc(Cindel *p1, Cindel *p2)
{
  return (p1->refstartpos > p2->refstartpos) ? 1 : (p1->refstartpos < p2->refstartpos) ? -1 :
    (p1->refendpos > p2->refendpos) ? 1 : (p1->refendpos < p2->refendpos) ? -1 : 0;
}

static int Id1SiteInc(register Calign **p1, register Calign **p2)
{
  Calign *align1 = p1[0];
  if(!align1) return 1;
  Calign *align2 = p2[0];
  if(!align2) return -1;
  if(DEBUG/* HERE >=2 */) assert(0 <= align1->mapid2 && align1->mapid2 < numX);
  if(DEBUG/* HERE >=2 */) assert(0 <= align2->mapid2 && align2->mapid2 < numX);

  Cmap *Ymap1 = YYmap[align1->mapid1];
  Cmap *Ymap2 = YYmap[align2->mapid1];
  if (Ymap1->id > Ymap2->id) 
    return 1;
  if (Ymap1->id < Ymap2->id)
    return -1;

  double LogPV1 = align1->logPV;
  double LogPV2 = align2->logPV;
  if(LogPV1 < LogPV2) 
    return 1;
  if(LogPV1 > LogPV2)
    return -1;

  return 0;
}

static int Id2SiteInc(register Calign **p1, register Calign **p2)
{
  Calign *align1 = p1[0];
  if(!align1) return 1;
  Calign *align2 = p2[0];
  if(!align2) return -1;
  if(DEBUG/* HERE >=2 */) assert(0 <= align1->mapid2 && align1->mapid2 < numX);
  if(DEBUG/* HERE >=2 */) assert(0 <= align2->mapid2 && align2->mapid2 < numX);
  Cmap *Xmap1 = XXmap[align1->mapid2];
  Cmap *Xmap2 = XXmap[align2->mapid2];
  if (Xmap1->id > Xmap2->id) 
    return 1;
  if (Xmap1->id < Xmap2->id)
    return -1;

  double LogPV1 = align1->logPV; // WAS77 RefSplit ? align1->logPV2 : align1->logPV;
  double LogPV2 = align2->logPV; // WAS77 RefSplit ? align2->logPV2 : align2->logPV;
  if(LogPV1 < LogPV2) 
    return 1;
  if(LogPV1 > LogPV2)
    return -1;

  Cmap *Ymap1 = YYmap[align1->mapid1];
  Cmap *Ymap2 = YYmap[align2->mapid1];
  if (Ymap1->id > Ymap2->id) 
    return 1;
  if (Ymap1->id < Ymap2->id)
    return -1;

  return 0;
}

/* sort alignments in order of logPV, id1, sites1[0], id2, sites2[0] (null alignments last) */
int LogPVSiteInc(register Calign **p1, register Calign **p2)
{
  register Calign *align1 = p1[0];
  register Calign *align2 = p2[0];
  if(!align1) return 1;
  if(!align2) return -1;

  register double LogPV1 = align1->logPV;
  register double LogPV2 = align2->logPV;
  if(LogPV1 < LogPV2) return 1;
  if(LogPV1 > LogPV2) return -1;

  Cmap *Ymap1 = YYmap[align1->mapid1];
  Cmap *Ymap2 = YYmap[align2->mapid1];
  int left1 = align1->sites1 ? (align1->sites1[0] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int left2 = align2->sites1 ? (align2->sites1[0] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  int ret = (Ymap1->id > Ymap2->id) ? 1 : (Ymap1->id < Ymap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  int numpairs1 = align1->numpairs;
  int numpairs2 = align2->numpairs;
  int right1 = align1->sites1 ? (align1->sites1[numpairs1-1] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int right2 = align2->sites1 ? (align2->sites1[numpairs2-1] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  ret = (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
  if(ret)
    return ret;

  Cmap *Xmap1 = XXmap[align1->mapid2];
  Cmap *Xmap2 = XXmap[align2->mapid2];
  int Xshift1 = Xmap1->origmap ? max(1,Xmap1->left[0]) - 1 : 0;
  int Xshift2 = Xmap2->origmap ? max(1,Xmap2->left[0]) - 1 : 0;
  int LJ1 = align1->sites2[0];
  int LJ2 = align2->sites2[0];
  int RJ1 = align1->sites2[numpairs1 - 1];
  int RJ2 = align2->sites2[numpairs2 - 1];
  left1 = align1->orientation ? Xshift1 + Xmap1->numsite[0] + 1 - RJ1 : Xshift1 + LJ1;
  left2 = align2->orientation ? Xshift2 + Xmap2->numsite[0] + 1 - RJ2 : Xshift2 + LJ2;
  ret = (Xmap1->id > Xmap2->id) ? 1 : (Xmap1->id < Xmap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  right1 = align1->orientation ? Xshift1 + Xmap1->numsite[0] + 1 - LJ1 : Xshift1 + RJ1;
  right2 = align2->orientation ? Xshift2 + Xmap2->numsite[0] + 1 - LJ2 : Xshift1 + RJ2;

  return (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
}

/* sort alignments in order of id2, id1, sites1[0] (null alignments last) */
int Site2Inc(register Calign **p1, register Calign **p2)
{
  register Calign *align1 = p1[0];
  register Calign *align2 = p2[0];
  if(!align1) return 1;
  if(!align2) return -1;

  Cmap *Xmap1 = XXmap[align1->mapid2];
  Cmap *Xmap2 = XXmap[align2->mapid2];
  int ret = (Xmap1->id > Xmap2->id) ? 1 : (Xmap1->id < Xmap2->id) ? -1 : 0;
  if(ret)
    return ret;

  Cmap *Ymap1 = YYmap[align1->mapid1];
  Cmap *Ymap2 = YYmap[align2->mapid1];
  int left1 = align1->sites1 ? (align1->sites1[0] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int left2 = align2->sites1 ? (align2->sites1[0] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  ret = (Ymap1->id > Ymap2->id) ? 1 : (Ymap1->id < Ymap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  int numpairs1 = align1->numpairs;
  int numpairs2 = align2->numpairs;
  int right1 = align1->sites1 ? (align1->sites1[numpairs1-1] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int right2 = align2->sites1 ? (align2->sites1[numpairs2-1] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  ret = (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
  return ret;
}

/* sort alignments in order of id1, sites1[0], sites1[N+1], id2, sites2[0],sites2[M+1] */
static int SiteInc(Calign **p1, Calign **p2)
{
  Calign *align1 = p1[0];
  Calign *align2 = p2[0];
  if(!align1) return 1;
  if(!align2) return -1;

  Cmap *Ymap1 = YYmap[align1->mapid1];
  Cmap *Ymap2 = YYmap[align2->mapid1];
  int left1 = align1->sites1 ? (align1->sites1[0] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int left2 = align2->sites1 ? (align2->sites1[0] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  int ret = (Ymap1->id > Ymap2->id) ? 1 : (Ymap1->id < Ymap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  int numpairs1 = align1->numpairs;
  int numpairs2 = align2->numpairs;
  int right1 = align1->sites1 ? (align1->sites1[numpairs1-1] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int right2 = align2->sites1 ? (align2->sites1[numpairs2-1] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  ret = (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
  if(ret)
    return ret;

  Cmap *Xmap1 = XXmap[align1->mapid2];
  Cmap *Xmap2 = XXmap[align2->mapid2];
  int Xshift1 = Xmap1->origmap ? max(1,Xmap1->left[0]) - 1 : 0;
  int Xshift2 = Xmap2->origmap ? max(1,Xmap2->left[0]) - 1 : 0;
  int LJ1 = align1->sites2[0];
  int LJ2 = align2->sites2[0];
  int RJ1 = align1->sites2[numpairs1 - 1];
  int RJ2 = align2->sites2[numpairs2 - 1];
  left1 = align1->orientation ? Xshift1 + Xmap1->numsite[0] + 1 - RJ1 : Xshift1 + LJ1;
  left2 = align2->orientation ? Xshift2 + Xmap2->numsite[0] + 1 - RJ2 : Xshift2 + LJ2;
  ret = (Xmap1->id > Xmap2->id) ? 1 : (Xmap1->id < Xmap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  right1 = align1->orientation ? Xshift1 + Xmap1->numsite[0] + 1 - LJ1 : Xshift1 + RJ1;
  right2 = align2->orientation ? Xshift2 + Xmap2->numsite[0] + 1 - LJ2 : Xshift1 + RJ2;

  return (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
}

/* 2-color version of LogPVSiteInc, only uses sites of first color to order */
int LogPVSiteInc2(register Calign **p1, register Calign **p2)
{
  register Calign *align1 = p1[0];
  register Calign *align2 = p2[0];
  if(!align1) return 1;
  if(!align2) return -1;

  register double LogPV1 = align1->logPV;
  register double LogPV2 = align2->logPV;
  if(LogPV1 < LogPV2) return 1;
  if(LogPV1 > LogPV2) return -1;

  Cmap *Ymap1 = YYmap[align1->mapid1];
  Cmap *Ymap2 = YYmap[align2->mapid1];
  int left1 = align1[1].sites1 ? (align1[1].sites1[0] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int left2 = align2[1].sites1 ? (align2[1].sites1[0] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  int ret = (Ymap1->id > Ymap2->id) ? 1 : (Ymap1->id < Ymap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  int numpairs1 = align1[1].numpairs;
  int numpairs2 = align2[1].numpairs;
  int right1 = align1[1].sites1 ? (align1[1].sites1[numpairs1-1] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int right2 = align2[1].sites1 ? (align2[1].sites1[numpairs2-1] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  ret = (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
  if(ret)
    return ret;

  Cmap *Xmap1 = XXmap[align1->mapid2];
  Cmap *Xmap2 = XXmap[align2->mapid2];
  int Xshift1 = Xmap1->origmap ? max(1,Xmap1->left[0]) - 1 : 0;
  int Xshift2 = Xmap2->origmap ? max(1,Xmap2->left[0]) - 1 : 0;
  int LJ1 = align1[1].sites2[0];
  int LJ2 = align2[1].sites2[0];
  int RJ1 = align1[1].sites2[numpairs1 - 1];
  int RJ2 = align2[1].sites2[numpairs2 - 1];
  left1 = align1->orientation ? Xshift1 + Xmap1->numsite[0] + 1 - RJ1 : Xshift1 + LJ1;
  left2 = align2->orientation ? Xshift2 + Xmap2->numsite[0] + 1 - RJ2 : Xshift2 + LJ2;
  ret = (Xmap1->id > Xmap2->id) ? 1 : (Xmap1->id < Xmap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  right1 = align1->orientation ? Xshift1 + Xmap1->numsite[0] + 1 - LJ1 : Xshift1 + RJ1;
  right2 = align2->orientation ? Xshift2 + Xmap2->numsite[0] + 1 - LJ2 : Xshift1 + RJ2;

  return (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
}

/* 2-color version of SiteInc, only uses sites of first color to order */
static int SiteInc2(Calign **p1, Calign **p2)
{
  Calign *align1 = p1[0];
  Calign *align2 = p2[0];
  if(!align1) return 1;
  if(!align2) return -1;

  Cmap *Ymap1 = YYmap[align1->mapid1];
  Cmap *Ymap2 = YYmap[align2->mapid1];
  int left1 = align1[1].sites1[0] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0);
  int left2 = align2[1].sites1[0] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0);
  int ret = (Ymap1->id > Ymap2->id) ? 1 : (Ymap1->id < Ymap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  int numpairs1 = align1[1].numpairs;
  int numpairs2 = align2[1].numpairs;
  int right1 = align1[1].sites1[numpairs1-1] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0);
  int right2 = align2[1].sites1[numpairs2-1] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0);
  ret = (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
  if(ret)
    return ret;

  Cmap *Xmap1 = XXmap[align1->mapid2];
  Cmap *Xmap2 = XXmap[align2->mapid2];
  int Xshift1 = Xmap1->origmap ? max(1,Xmap1->left[0]) - 1 : 0;
  int Xshift2 = Xmap2->origmap ? max(1,Xmap2->left[0]) - 1 : 0;
  int LJ1 = align1[1].sites2[0];
  int LJ2 = align2[1].sites2[0];
  int RJ1 = align1[1].sites2[numpairs1 - 1];
  int RJ2 = align2[1].sites2[numpairs2 - 1];
  left1 = align1->orientation ? Xshift1 + Xmap1->numsite[0] + 1 - RJ1 : Xshift1 + LJ1;
  left2 = align2->orientation ? Xshift2 + Xmap2->numsite[0] + 1 - RJ2 : Xshift2 + LJ2;
  ret = (Xmap1->id > Xmap2->id) ? 1 : (Xmap1->id < Xmap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 0;
  if(ret)
    return ret;

  right1 = align1->orientation ? Xshift1 + Xmap1->numsite[0] + 1 - LJ1 : Xshift1 + RJ1;
  right2 = align2->orientation ? Xshift2 + Xmap2->numsite[0] + 1 - LJ2 : Xshift1 + RJ2;

  return (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
}

class Cinterval {
public:
  long long id;/* reference/query id */
  FLOAT left, right;/* alignment interval */
  int mapid;/* index into YYmap[] or XXmap[] */
};

/* sort in ascending order of id and left ends */
static int IntervalIdSiteInc(Cinterval *p1, Cinterval *p2)
{
  if(p1->id > p2->id) return 1;
  if(p1->id < p2->id) return -1;
  return (p1->left > p2->left) ? 1 : (p1->left < p2->left) ? -1 : 0;
}

/* perform XMAP output to fp of Lalignment[start..end-1] 
   Use output prefix = basename (absolute path = absbasename, current working directory = cwd)
   
   Xscale and Yscale are use to scale alignment locations in the output : typically 1 except for pairwise alignment (or pairsplit), when this undoes -bpp in output

   If Indel != 0, also perform .indel output. If (dosvdetect) also appends indels to global sv_array[0..maxsv-1]
   if xmapentries != NULL, add xmap information to xmapentries[0 .. (end-start) -1]
 */

void output_xmapIO(char *filename, Calign **Lalignment, size_t start, size_t end, char *basename, char * &queryfilename, int Indel, xmapEntry **xmapentries, double Xscale, double Yscale, bool full)
{
  if(VERB>=2){
    printf("output_xmapIO:filename=%s,start=%lu,end=%lu,basename=%s,Indel=%d,full=%d\n",filename,start,end,basename,Indel,full);
    fflush(stdout);
  }

  /* get current directory pathname */
  char cwd[PATH_MAX];
  if(getcwd(cwd,PATH_MAX) == NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("getcwd() failed (errno=%d):%s\n",eno, err);
    exit(1);
  }

  /* convert basename into absolute pathname */
#ifndef WIN32
  char absbasename[PATH_MAX];
  if(basename[0]=='/')
    strcpy(absbasename,basename);
  else
    sprintf(absbasename,"%s/%s",cwd,basename);
  if(queryfilename && queryfilename[0] != '/'){/* Also convert queryfilename into an absolute pathname */
    char absname[PATH_MAX];
    sprintf(absname,"%s/%s",cwd,queryfilename);
    free(queryfilename);
    queryfilename = strdup(absname);
  }
#else
  char *absbasename = basename;
#endif

  if(checkFile(filename))
    return;

  if(VERB){
    printf("Generating %s\n",filename);
    fflush(stdout);
  }

  char *origfilename = strdup(filename);
  if(1){/* write to .tmp file first, then rename it later so -RefineOverwrite works correctly */
    int len = strlen(filename);
    sprintf(&filename[len], ".tmp");
  }
  char *tmpfilename = strdup(filename);

  FILE *fp;
  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("failed to open file %s for writing xmap file:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  size_t blen = 10*1024*1024;/* matches scif_io buffer size */
  char *write_buffer = (char *)malloc(blen);
  if(!write_buffer){
    printf("output_xmap:malloc(%llu) failed\n",(unsigned long long)blen);
    exit(1);
  }
  setbuffer(fp, write_buffer, blen);

  /* Also output indel SVs (if -indel) */
  FILE *fpSV = NULL;
  char filenameSV[PATH_MAX];
  char *write_bufferSV = NULL;
  if(Indel){
    sprintf(filenameSV,"%s.indel",basename);
    if(!checkFile(filenameSV)){
      if(VERB){
	printf("Generating %s\n",filenameSV);
	fflush(stdout);
      }
      if((fpSV = fopen(filenameSV,"w"))==NULL){
	int eno = errno;
	char *err = strerror(eno);
	printf("failed to open file %s for writing indel file:errno=%d:%s\n",filenameSV,eno,err);
	exit(1);
      }
    }
  }

  if(fpSV != NULL){
    write_bufferSV = (char *)malloc(blen);
    if(!write_bufferSV){
      printf("output_xmap:malloc(%llu) for write_bufferSV failed\n",(unsigned long long)blen);
      exit(1);
    }
    setbuffer(fpSV, write_bufferSV, blen);

    printversion(fpSV);

    fprintf(fpSV,"# SMAP File Version:\t0.0\n"); //  (0.0 was indels only; 0.1 has newer types)
    if(spots_filename[0] || (svcheck && numrefmaps > 0)){/* Reference Aligner : refer to condensed version of reference CMAP */
      if(full)
	fprintf(fpSV,"# Reference Maps From:\t%s_full_r.cmap\n", absbasename);
      else
	fprintf(fpSV,"# Reference Maps From:\t%s_r.cmap\n", absbasename);
      fprintf(fpSV,"# Query Maps From:\t%s\n",queryfilename);/* This may actually be a modified cmap file */
    } else {
      if(mres > 0.001 || mresSD > 0.0 || fabs(PixelLen - origPixelLen) > 1e-12){
	if(full)
	  fprintf(fpSV,"# Reference Maps From:\t%s_full_r.cmap\n", absbasename);
	else
	  fprintf(fpSV,"# Reference Maps From:\t%s_r.cmap\n", absbasename);
	fprintf(fpSV,"# Query Maps From:\t%s_q.cmap\n", absbasename);
      } else {
	if(full)
	  fprintf(fpSV,"# Reference Maps From:\t%s_full.cmap", vfixx_filename[0]);
	else
	  fprintf(fpSV,"# Reference Maps From:\t%s.cmap", vfixx_filename[0]);
	for(register int i = 1; i < RefIndex; i++)
#ifndef WIN32
	  if(vfixx_filename[i][0] == '/')
	    fprintf(fpSV,",%s",vfixx_filename[i]);
	  else
	    fprintf(fpSV,",%s/%s",cwd,vfixx_filename[i]);
#else
	    fprintf(fpSV,",%s",vfixx_filename[i]);
#endif
	fprintf(fpSV,"\n");
#ifndef WIN32
	if(vfixx_filename[QueryIndex][0] == '/')
	  fprintf(fpSV,"# Query Maps From:\t%s",vfixx_filename[QueryIndex]);
	else
	  fprintf(fpSV,"# Query Maps From:\t%s/%s",cwd,vfixx_filename[QueryIndex]);
#else
	fprintf(fpSV,"# Query Maps From:\t%s",vfixx_filename[QueryIndex]);	
#endif
	for(register int i = QueryIndex+1; i < num_files;i++)
#ifndef WIN32
	  if(vfixx_filename[i][0] == '/')
	    fprintf(fpSV,",%s",vfixx_filename[i]);
	  else
	    fprintf(fpSV,",%s/%s",cwd,vfixx_filename[i]);
#else
	    fprintf(fpSV,",%s",vfixx_filename[i]);
#endif
	fprintf(fpSV,"\n");
      }
    }
#ifndef WIN32
    if(filename[0] == '/')
      fprintf(fpSV,"# Xmap Entries From:\t%s\n",filename);    
    else
      fprintf(fpSV,"# Xmap Entries From:\t%s/%s\n",cwd,filename);    
#else
      fprintf(fpSV,"# Xmap Entries From:\t%s\n",filename);    
#endif
    fprintf(fpSV,"#h SmapEntryID\tQryContigID\tRefcontigID1\tRefcontigID2\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tType\tXmapID1\tXmapID2\n");
    fprintf(fpSV,"#f int        \tint        \tint         \tint         \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat  \tstring \tint    \tint    \n");
  }

  /* write out commandline */
  printversion(fp);

  fprintf(fp,"# XMAP File Version:\t0.2\n");
  fprintf(fp,"# Label Channels:\t%d\n", colors);

  if(spots_filename[0] || (svcheck && numrefmaps > 0)){/* Reference Aligner : refer to condensed version of reference CMAP */
    if(full)
      fprintf(fp,"# Reference Maps From:\t%s_full_r.cmap\n", absbasename);
    else
      fprintf(fp,"# Reference Maps From:\t%s_r.cmap\n", absbasename);
    fprintf(fp,"# Query Maps From:\t%s",queryfilename);/* This may actually be a modified cmap file */
  } else {
    if(mres > 0.001 || mresSD > 0.0 || fabs(PixelLen - origPixelLen) > 1e-12){
      if(full)
	fprintf(fp,"# Reference Maps From:\t%s_full_r.cmap\n", absbasename);
      else
	fprintf(fp,"# Reference Maps From:\t%s_r.cmap\n", absbasename);
      fprintf(fp,"# Query Maps From:\t%s_q.cmap",absbasename);
    } else {
      fprintf(fp,"# Reference Maps From:\t%s", vfixx_filename[0]);
      for(register int i = 1; i < RefIndex; i++)
#ifndef WIN32
	if(vfixx_filename[i][0] == '/')
	  fprintf(fp,",%s",vfixx_filename[i]);
	else
	  fprintf(fp,",%s/%s",cwd,vfixx_filename[i]);
#else
        fprintf(fp,",%s",vfixx_filename[i]);
#endif
      fprintf(fp,"\n");

#ifndef WIN32
      if(vfixx_filename[QueryIndex][0] == '/')
	fprintf(fp,"# Query Maps From:\t%s",vfixx_filename[QueryIndex]);
      else
	fprintf(fp,"# Query Maps From:\t%s/%s",cwd,vfixx_filename[QueryIndex]);
#else
      fprintf(fp,"# Query Maps From:\t%s",vfixx_filename[QueryIndex]);
#endif
      for(register int i = QueryIndex+1; i < num_files;i++)
#ifndef WIN32
	if(vfixx_filename[i][0] == '/')
	  fprintf(fp,",%s",vfixx_filename[i]);
	else
	  fprintf(fp,",%s/%s",cwd,vfixx_filename[i]);
#else
	  fprintf(fp,",%s",vfixx_filename[i]);
#endif
    }
  }
  fprintf(fp,"\n");

  if((PairSplit || (MultiMatches && RefSplit)) && dosvdetect){ /* see parameters.{h,cpp} */
    fprintf(fp,"# Smap Entries From:\t%s.smap\n",basename);
  }

  if(SecondBest && SinglelineXmap){
    fprintf(fp,"#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment\tQryEndPos2\tRefEndPos2\tOrientation2\tConfidence2\tNumMatch2\n");
    fprintf(fp,"#f int        \tint        \tint        \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat     \tstring \tflat  \tfloat \tint         \tstring   \tfloat     \tfloat     \tstring      \tfloat      \tint      \n");
  } else if(SecondBest && !SinglelineXmap){
    fprintf(fp,"#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment\tSecondID\n");
    fprintf(fp,"#f int        \tint        \tint        \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat     \tstring \tfloat \tfloat \tint         \tstring   \tint\n");
    fprintf(fp,"# SecondID refers to XmapEntryID of Second Best alignment for same query and reference (typically next line): 0 if none, -1 if current line is Second Best alignment\n");
  } else if(0 && MultiMatches){
    fprintf(fp,"#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment\tBestID\n");
    fprintf(fp,"#f int        \tint        \tint        \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat     \tstring \tfloat \tfloat \tint         \tstring   \tint\n");
  } else if(XMAPmapWT){
    fprintf(fp,"#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment\tMapWt\n");
    fprintf(fp,"#f int        \tint        \tint        \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat     \tstring \tfloat \tfloat \tint         \tstring   \tfloat\n");
  } else {
    fprintf(fp,"#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment\n");
    fprintf(fp,"#f int        \tint        \tint        \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat     \tstring \tfloat \tfloat \tint         \tstring   \n");
  }

  /* allocate interval array to keep track of total and unique alignment lengths */
  Cinterval *Rinterval = new Cinterval[(end - start + 1)*2];
  Cinterval *Qinterval = &Rinterval[end - start + 1];
  int RintervalCnt = 0, QintervalCnt = 0;
  
  size_t aligncnt = 0;
  size_t SmapCnt = 0;

  double Ilog10 = 1.0/log(10.0);
  double EndOutlier = log(max(1.0e-300,PoutlierEnd)) + log(max(0.0001,1.0-PoutlierEnd));;
  double Outlier = log(max(1.0e-300,Poutlier));

  if(INDEL_TRACE){
    printf("Ilog10=%0.6f,EndOutlier=%0.6f,Outlier=%0.6f,Lalignment[%lu..%lu]\n",Ilog10,EndOutlier,Outlier,start,end-1);
    fflush(stdout);
  }

  for(size_t i = start; i < end; i++){
    Calign *p = Lalignment[i];
    if(!p)
      continue;

    if(INDEL_TRACE >= 2){
      printf("i=%lu(%lu..%lu),aligncnt=%lu:Lalignment[i]=%p,Lalignment[i]->align2=%p,refid=%d,mapid=%d,or=%d,chimpair=%d:numpairs=%d,score=%0.3f,logPV=%0.2f(at loop start)\n",
	     i,start,end-1,aligncnt,p,p->align2,p->mapid1,p->mapid2,p->orientation,p->chimpair,p->numpairs,p->score,p->logPV);
      fflush(stdout);
    }

    int U = p->numpairs;

    int SecondLine = 0;
    if(SecondBest && !SinglelineXmap && p->chimpair < 0){
      if(DEBUG/* HERE >=2 */) assert(!(U<=0 || p->score <= ScoreThreshold2 || p->logPV <= LogPvThreshold2));
      SecondLine = 1;
    }

    /* Below p->chimpair is changed to the XmapEntryID in the output file (0 means no output) */

    if(!SecondLine){
      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold))){
	if(VERB>=2){
	  printf("Lalignment[%lu]= %p: numpairs=%d,score=%0.3f,logPV=%0.2f, align2= %p, chimpair=%d: skipping\n",i, p, U, p->score, p->logPV, p->align2, p->chimpair);
	  fflush(stdout);
	}
	p->chimpair = 0;
	continue;
      }
    } else {
      if(DEBUG) assert(U > 0);
    }

    int rid = p->mapid1;
    int mid = p->mapid2;
    if(DEBUG && !(0 <= rid && rid < numY)){
      printf("i=%lu(%lu..%lu),aligncnt=%lu:Lalignment[i]=%p,Lalignment[i]->align2=%p,refid=%d,mapid=%d,or=%d,chimpair=%d:numpairs=%d,score=%0.3f,logPV=%0.2f(at loop start)\n",
	     i,start,end-1,aligncnt,p,p->align2,p->mapid1,p->mapid2,p->orientation,p->chimpair,p->numpairs,p->score,p->logPV);
      printf("\n rid=%d,numY=%d,mid=%d,numX=%d\n",rid,numY,mid,numX);
      fflush(stdout);
      assert(0 <= rid && rid < numY);
    }
    if(DEBUG && !(0 <= mid && mid <= numX)){
      printf("i=%lu(%lu..%lu),aligncnt=%lu:Lalignment[i]=%p,Lalignment[i]->align2=%p,refid=%d,mapid=%d,or=%d,chimpair=%d:numpairs=%d,score=%0.3f,logPV=%0.2f(at loop start)\n",
	     i,start,end-1,aligncnt,p,p->align2,p->mapid1,p->mapid2,p->orientation,p->chimpair,p->numpairs,p->score,p->logPV);
      printf("\n rid=%d,numY=%d,mid=%d,numX=%d\n",rid,numY,mid,numX);
      fflush(stdout);
      
      assert(0 <= mid && mid < numX);
    }
    //    Cmap *Ymap = YYmap[rid];
    Cmap *Xmap = XXmap[mid];

    Cmap *origmap = Xmap;
    while(origmap->origmap){
      Xmap->origmap = origmap = origmap->origmap;
      if(origmap->origmap){/* otherwise origmap->left[] is not defined */
	Xmap->left[0] += max(0,origmap->left[0]-1);
	Xmap->right[0] += max(0,origmap->left[0]-1);
      }
    }
    if(BestRef && origmap->align->mapid1 != rid){
      p->chimpair = 0;
      continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
    }

    aligncnt++;
    p->chimpair = (p->chimpair < 0) ? -aligncnt : aligncnt;
  }

  if(VERB>=2){
    printf("output_xmapIO: aligncnt=%lu, start=%lu, end=%lu\n",aligncnt,start,end);
    fflush(stdout);
  }

  size_t totaligncnt = aligncnt;

  /* now all |p->chimpair| are the planned XmapEntryID in output file (or 0 if alignment is not output) */

  int IndelRangeExpandCnt = 0, IndelTotCnt = 0;

  extern double raYlambda;// see RGentigRefScore.h
  double AvgInterval = raYlambda;

  for(size_t i = start; i < end; i++){
    Calign *p = Lalignment[i];
    if(!p || p->chimpair==0)
      continue;

    int U = p->numpairs;

    int SecondLine = 0;
    if(SecondBest && !SinglelineXmap && p->chimpair < 0)
      SecondLine = 1;

    aligncnt = abs(p->chimpair);

    int rid = p->mapid1;
    if(DEBUG) assert(0 <= rid && rid < numY);
    int mid = p->mapid2;
    if(DEBUG && !(0 <= mid && mid <= numX)){
      printf("i=%lu(%lu..%lu),aligncnt=%lu:Lalignment[i]=%p,Lalignment[i]->align2=%p,refid=%d,mapid=%d,or=%d,chimpair=%d:numpairs=%d,score=%0.3f,logPV=%0.2f\n",
	     i,start,end-1,aligncnt,p,p->align2,p->mapid1,p->mapid2,p->orientation,p->chimpair,p->numpairs,p->score,p->logPV);
      printf("\n rid=%d,numY=%d,mid=%d,numX=%d\n",rid,numY,mid,numX);
      fflush(stdout);
      
      assert(0 <= mid && mid < numX);
    }
    Cmap *Ymap = YYmap[rid];
    Cmap *Xmap = XXmap[mid];

    Cmap *origmap = Xmap;
    while(origmap->origmap){
      Xmap->origmap = origmap = origmap->origmap;
      if(origmap->origmap){/* otherwise origmap->left[] is not defined */
	Xmap->left[0] += max(0,origmap->left[0]-1);
	Xmap->right[0] += max(0,origmap->left[0]-1);
      }
    }

    if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
      printf("aligncnt=%lu: Xmap Lalignment[%lu]:mid=%lld,rid=%lld,or=%d:numpairs=%d,score=%0.6f,logPV=%0.2f,chimpair=%d\n",
	     aligncnt,i,Xmap->id,Ymap->id,p->orientation,U,p->score,p->logPV,p->chimpair);
      fflush(stdout);
    }

    int M = Xmap->numsite[0];
    int N = Ymap->numsite[0];
    FLOAT *Y = Ymap->site[0];
    FLOAT *X = Xmap->site[0];
    Cmap *origYmap = Ymap->origmap ? Ymap->origmap : Ymap;
    Cmap *origXmap = Xmap->origmap ? Xmap->origmap : Xmap;
    FLOAT *origY = Ymap->origmap ? Ymap->origmap->site[0] : Y;
    FLOAT *origX = Xmap->origmap ? Xmap->origmap->site[0] : X;
    int origN = Ymap->origmap ? Ymap->origmap->numsite[0] : N;
    int origM = Xmap->origmap ? Xmap->origmap->numsite[0] : M;
    if(DEBUG) assert(Xmap->origmap==0 || Xmap->origmap->origmap==0);/* If this fails fix origmap link in pairalign() or refalign() to point to root map */
    if(DEBUG) assert(Ymap->origmap==0 || Ymap->origmap->origmap==0);/* If this fails fix origmap link (same as for Xmap above) to point to root map */
    int Yshift = Ymap->origmap ? max(1,Ymap->left[0]) - 1 : 0;
    int Xshift = Xmap->origmap ? max(1,Xmap->left[0]) - 1 : 0;
    
    int LI = p->sites1[0];
    int LK = (resSD[0] > 0.0) ? p->sitesK1[0] : 0;
    if(DEBUG) assert(1 <= LI && LI <= N);
    int LJ = p->sites2[0];
    int RI = p->sites1[U-1];
    //    register int RK = (resSD > 0.0 ) ? p->sitesK1[0] : 0;
    if(DEBUG) assert(LI <= RI && RI <= N);
    int RJ = p->sites2[U-1];
    if(DEBUG) assert(1 <= LJ && LJ <= RJ);// NEW9

    // Get the indices, subsequently use to compute number of un-aligned labels 
    /* NOTE : here qrystartidx is the absolute label index that is aligned with LI, while in structuralVariation.h, qrystartidx is the smaller of absolute label index aligned with LI or RI */
    int qrystartidx = p->orientation ? Xshift + M + 1 - LJ : Xshift + LJ;
    int qryendidx   = p->orientation ? Xshift + M + 1 - RJ : Xshift + RJ;
    double qrystartpos = origX[qrystartidx]*1000.0*Xscale; /* QryStartPos */
    double qryendpos   = origX[qryendidx  ]*1000.0*Xscale; /* QryEndPos */
    double qryLen = origX[origM+1]*1000.0*Xscale; /* QryLen */
    int refstartidx = Yshift + LI - LK;
    int refendidx   = Yshift + RI;
    double refstartpos = origY[refstartidx]*1000.0*Yscale; /* RefStartPos */
    double refendpos   = origY[refendidx  ]*1000.0*Yscale; /* RefEndPos */
    double refLen = origY[origN+1]*1000.0*Yscale;/* RefLen */
    
    if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
      printf("    LI=%d,LK=%d,LJ=%d,RI=%d,RJ=%d,Xshift=%d,Yshift=%d,Xscale=%0.8f,Yscale=%0.8f:qrystartidx=%d,qryendidx=%d,qrystartpos=%0.6f,qryendpos=%0.6f\n",
	     LI,LK,LJ,RI,RJ,Xshift,Yshift,Xscale,Yscale,qrystartidx,qryendidx,qrystartpos,qryendpos);
      fflush(stdout);
    }

    /* add interval information to Rinterval[],Qinterval[] */
    Cinterval *ri = &Rinterval[RintervalCnt++];
    Cinterval *qi = &Qinterval[QintervalCnt++];
    ri->id = YYmap[rid]->id;
    ri->mapid = rid;
    ri->left = refstartpos;
    ri->right = refendpos;
    qi->id = XXmap[mid]->id;
    qi->mapid = mid;
    qi->left = min(qrystartpos,qryendpos);
    qi->right = max(qrystartpos,qryendpos);

    fprintf(fp,"%llu\t%lld\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%c\t%0.2f\t",
	    (unsigned long long)aligncnt, /* XmapEntryID */
	    Xmap->id, /* QryContigID */
	    Ymap->id, /* RefContigID */
	    qrystartpos,
	    qryendpos,
	    refstartpos,
	    refendpos,
	    p->orientation ? '-' : '+', /* Orientation */
	    p->logPV/*zscore(pow(10.0,-max(0.0,p->logPV)))*/ /* Confidence */
	    );

    int nindel = 0; //count n indels added for use in making xmapentry objects
    if(Indel){/* output internal outliers of current alignment as indels */
      int lastI = p->sites1[0];
      int lastK = (resSD[0] > 0.0) ? p->sitesK1[0] : 0;
      int lastJ = p->sites2[0];
      int IL = lastI, KL = lastK, JL = lastJ;
      int I,K,J;

      int psI = 0, psK = 0, psJ = 0, pslastI = 0, pslastK = 0, pslastJ = 0;// last output Indel boundaries : used to avoid duplicates

      for(int k = 1; k < U; lastI = I, lastK = K, lastJ = J, k++){
	I = p->sites1[k];
	K = (resSD[0] > 0.0) ? p->sitesK1[k] : 0;
	J = p->sites2[k];

	int outlier = (p->iscore[k] > p->outscore[k] + 0.01 || (outlierExtend && (lastJ >= J || lastI >= I-K)) || p->outscore[k] < 0.0) ? 1 : 0;
	if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE && outlier){
	  int qrystartidx = p->orientation ? Xshift + M + 1 - lastJ : Xshift + lastJ;
	  int qryendidx   = p->orientation ? Xshift + M + 1 - J : Xshift + J;
	  double qrystartpos = origX[qrystartidx]*1000.0*Xscale; /* QryStartPos */
	  double qryendpos   = origX[qryendidx  ]*1000.0*Xscale; /* QryEndPos */
	  int refstartidx    = Yshift+lastI-lastK;
	  int refendidx      = Yshift+I;
	  double refstartpos = origY[refstartidx]*1000.0*Yscale; /* RefStartPos */
	  double refendpos   = origY[refendidx  ]*1000.0*Yscale; /* RefEndPos */

	  printf("xid=%lu:rid=%lld,mid=%lld,or=%d:k=%d/%d:I=%d,%d,K=%d,%d,J=%d,%d: iscore[k]= %0.6f, outscore[k]= %0.6f, outlier=%d, qry=%d..%d(%0.3f..%0.3f),ref=%d..%d(%0.3f..%0.3f)\n",
		 aligncnt,Ymap->id,Xmap->id,p->orientation, k,U,lastI,I,lastK,K,lastJ,J,p->iscore[k],p->outscore[k],outlier,qrystartidx,qryendidx,qrystartpos,qryendpos,refstartidx,refendidx,refstartpos,refendpos);
	  printf("\t CmapChimQuality=%d,IndelChiSq= %0.3e, BreakChiSq= %0.3f, NoBreak=%d\n",CmapChimQuality,IndelChiSq,BreakChiSq,NoBreak);
	  fflush(stdout);
	}

	if(!outlier)/* not an outlier */
	  continue;

	/*  merge nearby outliers into a single SV ? This is currently done by -outlierExtend or -refSplit (and subsequently by -svcheck) */

	/* outlier lies between sites lastI .. I-K on reference and lastJ .. J on query */
	int qrystartidx = p->orientation ? Xshift + M + 1 - lastJ : Xshift + lastJ;
	int qryendidx   = p->orientation ? Xshift + M + 1 - J : Xshift + J;
	double qrystartpos = origX[qrystartidx]*1000.0*Xscale; /* QryStartPos */
	double qryendpos   = origX[qryendidx  ]*1000.0*Xscale; /* QryEndPos */
	int refstartidx    = Yshift+lastI-lastK;
	int refendidx      = Yshift+I;
	double refstartpos = origY[refstartidx]*1000.0*Yscale; /* RefStartPos */
	double refendpos   = origY[refendidx  ]*1000.0*Yscale; /* RefEndPos */

	double origqsiz = p->orientation ? (qrystartpos - qryendpos) : (qryendpos - qrystartpos);
	double origqrystartpos = qrystartpos;
	double origqryendpos = qryendpos;
	double origrefstartpos = (origY[Yshift+lastI]*1000.0*Yscale + refstartpos)*0.5;
	double origrefendpos = (refendpos + origY[Yshift+I-K]*1000.0*Yscale)*0.5;

	if(DEBUG/* HERE >=2 */) assert(origqsiz >= 0.0);

	/* include half of de-res interval in query interval, so query interval can be replaced with ref interval without introducing a length bias */
	if(!p->orientation){
	  if(lastK > 0)
	    qrystartpos -= 0.5*(origY[Yshift+lastI] - origY[Yshift+lastI-lastK]) * 1000.0 * Yscale;
	  if(K > 0)
	    qryendpos += 0.5*(origY[Yshift+I] - origY[Yshift+I-K]) * 1000.0 * Yscale;
	} else {/* query runs from high to low positions along alignment */ 
	  if(lastK > 0)
	    qrystartpos += 0.5*(origY[Yshift+lastI] - origY[Yshift+lastI-lastK]) * 1000.0 * Yscale;
	  if(K > 0)
	    qryendpos -= 0.5*(origY[Yshift+I] - origY[Yshift+I-K]) * 1000.0 * Yscale;
	}
	if(DEBUG/* HERE >=2 */) assert(refstartpos <= refendpos);
	double newqsiz =  p->orientation ? (qrystartpos - qryendpos) : (qryendpos - qrystartpos);
	if(DEBUG/* HERE >=2 */ && !(newqsiz >= origqsiz)){
	  printf("Yid=%lld,Xid=%lld,or=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f:query=%0.3f..%0.3f,ref=%0.3f..%0.3f,qryindex=%d..%d(M=%d)\n",
		 Ymap->id,Xmap->id,p->orientation,I,K,J,lastI,lastK,lastJ,p->outscore[k],p->iscore[k],qrystartpos,qryendpos,refstartpos,refendpos,qrystartidx,qryendidx,M);
	  printf("Xmap->origmap=%p,Xmap=%p,Xshift=%d,origX=%p,X=%p,Xscale=%0.8f\n",Xmap->origmap,Xmap,Xshift,origX,X,Xscale);
	  for(int I = 0; I <= M+1; I++)
	    printf("X[%d]=%0.6f %c\n",I,X[I], (I > 0 && X[I] < X[I-1]) ? '!' : ' ');
	  fflush(stdout);
	  assert(newqsiz >= origqsiz);
	}
	if(DEBUG/* HERE >=2 */) assert(J > lastJ);
	if(DEBUG/* HERE >=2 */) assert(I-K > lastI);

	if(DEBUG && IndelChiSq > 0.0) assert(!(BreakChiSq > 0.0 && NoBreak));

	if((IndelChiSq > 0.0 || (BreakChiSq > 0.0 && NoBreak)) && CmapChimQuality >= 3 && spots_filename[0] && qryendpos - qrystartpos < refendpos - refstartpos){// ignore outlier based in -indelChiSq
	  double OutlierKB = 0.001 * fabs(refendpos - refstartpos - qryendpos + qrystartpos);// indel size in KB

	  float *OutlierFrac = Xmap->OutlierFrac[0];
	  float *FragSd = Xmap->FragSd[0];
	  float *ExpSd = Xmap->ExpSd[0];
	  float *FragCov = Xmap->FragCov[0];
	  float *FragChiSq = Xmap->FragChiSq[0];
	  if(DEBUG/* HERE HERE >=2 */) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);
	  if(DEBUG/* HERE HERE >=2 */) assert(AvgInterval > 0.0);

	  float *ChimQuality = Xmap->ChimQuality[0];
	  float *ChimNorm = Xmap->ChimNorm[0];
	  float *SegDupL = Xmap->SegDupL[0];
	  float *SegDupR = Xmap->SegDupR[0];
	  float *FragileEndL = Xmap->FragileEndL[0];
	  float *FragileEndR = Xmap->FragileEndR[0];

	  double LRange = (origY[Yshift + I] - origY[Yshift + lastI]) * Yscale;
	  double LRangeRatio = LRange / AvgInterval;

	  int jmin = qrystartidx;
	  int jmax = qryendidx;
	  if(VERB/* HERE HERE >=2 && Ymap->id == 1 && Xmap->id == 5298 */){
	    printf("xid=%lu:Yid=%lld,Xid=%lld,or=%d:I=%d..%d,K=%d..%d,J=%d..%d(j=%d..%d):outscore=%0.6f,iscore=%0.6f:query=%0.3f..%0.3f,ref=%0.3f..%0.3f,qryindex=%d..%d(M=%d),OutlierKB= %0.3f\n",
		   aligncnt,Ymap->id,Xmap->id,p->orientation,lastI,I,lastK,K,lastJ,J,jmin,jmax,p->outscore[k],p->iscore[k],qrystartpos,qryendpos,refstartpos,refendpos,qrystartidx,qryendidx,M,OutlierKB);
	    printf("\t    MaxSize= %0.1f, MaxOutlierFrac= %0.1f, MaxCov= %0.2f, ChiSq= %0.4e, MinSf= %0.4f, MinSr= %0.4f, SdRatio= %0.4f, MinInterval= %0.4f, MinCov= %0.2f, MinCovR= %0.2f(%0.2f), NoBreak=%d\n",
		   IndelMaxSize, IndelFrac,IndelMaxCov,IndelChiSq,IndelMinSf,IndelMinSr,IndelSdRatio,IndelMinInterval, IndelMinCov, IndelMinCovScale, IndelMinCovRatio, (BreakChiSq > 0.0) ? NoBreak : 0);
	    fflush(stdout);
	  }	    

	  if(IndelMinCovRatio > 0 && I-K-lastI != J-lastJ)
	    LRangeRatio = IndelMinCovRatio;

	  int j = jmin;
	  for(; j < jmax; j++){
	    double SRange = (origX[j+1] - origX[j]) * Xscale;
	    if(DEBUG && !NoBreak) assert(FragChiSq[j] >= 0.0);

	    double FragCovj = FragCov[j];
	    if(IndelChimNorm > 0.0 && ChimNorm[j] > 0.0 && ChimNorm[j+1] > 0.0)
	      FragCovj = min(FragCovj, min(ChimNorm[j]*ChimQuality[j],ChimNorm[j+1]*ChimQuality[j+1]) * IndelChimNorm * 0.01);

	    if((BreakChiSq > 0.0 && NoBreak) ? (FragChiSq[j] < 0.0 && refendpos - refstartpos > qryendpos - qrystartpos && OutlierKB < IndelMaxSize) :
	       (OutlierKB < IndelMaxSize && refendpos - refstartpos > qryendpos - qrystartpos && LRange > IndelMinInterval && 
		(OutlierFrac[j] > IndelFrac || 
		 max(FragileEndL[j]+FragileEndR[j],FragileEndL[j+1]+FragileEndR[j+1]) > IndelFragile ||
		 max(SegDupL[j]+SegDupR[j],SegDupL[j+1]+SegDupR[j+1]) > IndelSegDup ||
		 (IndelWin > 0 && FragChiSq[j] < 0.0) || 
		 (FragCovj < IndelMaxCov && 
		  ((0.0 <= FragChiSq[j] && FragChiSq[j] < IndelChiSq &&
		    FragSd[j] > ExpSd[j] + IndelMinSf + IndelMinSr * SRange &&
		    FragSd[j] > ExpSd[j] * IndelSdRatio) ||		     
		   FragCovj < IndelMinCov + ((IndelMinCovRatio > 0 && I-K-lastI != J-lastJ) ? min(IndelMinCovRatio,LRangeRatio) : LRangeRatio) *  IndelMinCovScale))))){
	      if(VERB/* HERE HERE >=2 && Ymap->id == 1 && Xmap->id == 5298 */ && !(BreakChiSq > 0.0 && NoBreak)){
		/*		printf("Yid=%lld,Xid=%lld,or=%d:I=%d..%d,K=%d..%d,J=%d..%d:outscore=%0.6f,iscore=%0.6f:query=%0.3f..%0.3f,ref=%0.3f..%0.3f,qryindex=%d..%d(M=%d)\n",
				Ymap->id,Xmap->id,p->orientation,lastI,I,lastK,K,lastJ,J,p->outscore[k],p->iscore[k],qrystartpos,qryendpos,refstartpos,refendpos,qrystartidx,qryendidx,M);*/
		printf("\t j=%d(%d..%d):J=%d..%d:LRange= %0.4f(Ratio= %0.3f), SRange= %0.4f, OutlierFrac= %0.1f, FragChiSq= %0.4e, OutlierKB= %0.4f, FragCov= %0.2f(N1=%0.2f,%0.2f), FragSd= %0.1f, ExpSd= %0.1f : skipping outlier\n",
		       j,jmin,jmax,lastJ,J,LRange,LRangeRatio,SRange, OutlierFrac[j], FragChiSq[j], OutlierKB, FragCov[j], ChimNorm[j]*ChimQuality[j]*0.01, ChimNorm[j+1]*ChimQuality[j+1]*0.01, FragSd[j],ExpSd[j]);
		/*		printf("\t    MaxSize= %0.1f, MaxOutlierFrac= %0.1f, MaxCov= %0.2f, ChiSq= %0.4e, MinSf= %0.4f, MinSr= %0.4f, SdRatio= %0.4f, MinInterval= %0.4f, MinCov= %0.2f, MinCovR= %0.2f\n",
				IndelMaxSize, IndelFrac,IndelMaxCov,IndelChiSq,IndelMinSf,IndelMinSr,IndelSdRatio,IndelMinInterval, IndelMinCov, IndelMinCovScale);*/
		fflush(stdout);
	      }

	      break;
	    }
	    if(VERB/* HERE HERE >=3 && Ymap->id == 1 && (Xmap->id == 409 || Xmap->id==638) */ && !(BreakChiSq > 0.0 && NoBreak)){
	      printf("\t j=%d(%d..%d):J=%d..%d:LRange= %0.4f(Ratio= %0.3f), SRange= %0.4f, OutlierFrac= %0.1f, FragChiSq= %0.4e, OutlierKB= %0.4f, FragCov= %0.2f(N1=%0.2f,%0.2f), FragSd= %0.1f, ExpSd= %0.1f\n",
		     j,jmin,jmax,lastJ,J,LRange,LRangeRatio,SRange, OutlierFrac[j], FragChiSq[j], OutlierKB, FragCov[j], ChimNorm[j]*ChimQuality[j]*0.01,ChimNorm[j+1]*ChimQuality[j]*0.01,FragSd[j],ExpSd[j]);
	      //	      printf("\t\t (OutlierKB < IndelMaxSize) = %d\n",OutlierKB < IndelMaxSize ? 1 : 0);
	      fflush(stdout);
	    }	    
	  }

	  if(j < jmax)
	    continue;// outlier will be ignored
	}

	/* compute confidence in outlier based on lesser of:
	   1. Confidence in entire alignment
	   2. Confidence that this region is an outlier (vs a regular alignment)
	   3. Confidence that this region is an outlier (vs an EndOutlier in either direction).
	   NOTE : should really be adding the Pvalues, but this is approximately correct and quicker
	*/
	double centerConfidence, outlierConfidence,PenSm = 0.0;
	int sI = I, sK = K, sJ = J, slastI = lastI, slastK = lastK, slastJ = lastJ;// values to output to .indel or .smap : can be updated without modifying loop variables I,lastI etc 
	double SVsizeFix = 0.0;/* change in SVsize caused by -svSmallDelConfirm and/or -svIndelRangeExpand : If +ve deltaY will be expanded by this value, else deltaX will be expanded by -SVsizeFix */

	if(spots_filename[0] || (svcheck && numrefmaps > 0)){/* Reference Aligner : need to recompute interval score to exclude misresolved label penalty */
	  // declared in RGentigRefScore.h
	  extern void SintDetailRA(double X, double Y, int m, int n, int J, int I, int K /* right */, int T /* left */, FLOAT *Yref, double &Bias, double &Pen, double &Gauss, double &PenSm, int verb);
	  
	  extern RFLOAT Frate;
	  double  Bias,Pen,Gauss;
	  double deltaX = fabs(origX[qryendidx] - origX[qrystartidx]);
	  double deltaY = Yc(origY, I, K) - Yc(origY, lastI, lastK);
	  SintDetailRA(deltaX, deltaY, J - lastJ, I-K-lastI, J, I, K, lastK, Y, Bias,Pen,Gauss,PenSm, 
			   (VERB>=2 && Ymap->id == 10 && Xmap->id == 111) ? 1 : 0);
	  if(FRATE_FIX1 && smap_IndelConfBySize){
	    Gauss += deltaX * Frate;
	    Pen -= deltaX * Frate;
	  }
	  double outscore = smap_IndelConfBySize ? Gauss : Bias + Pen + Gauss; // NOTE : excludes PenSm which is the misresolved penalty for the right end of the interval
	  double SVsize = deltaY - deltaX;

	  if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
	    double OutlierKB = 0.001 * fabs(refendpos - refstartpos - qryendpos + qrystartpos);// indel size
	    printf("Yid=%lld,Xid=%lld,or=%d:I=%d..%d,K=%d..%d,J=%d..%d,OutlierKB=%0.4f: checking outscore[%d]/L10= %0.4f -> %0.4f (Bias= %0.4f, Gauss= %0.4f, Pen= %0.4f, PenSm= %0.4f), ref= %0.3f..%0.3f, qry= %0.3f..%0.3f\n",
		   Ymap->id,Xmap->id,p->orientation,lastI,I,lastK,K,lastJ,J,OutlierKB, k,p->outscore[k]*Ilog10, outscore*Ilog10, Bias, Gauss, Pen, PenSm, refstartpos*1e-3, refendpos*1e-3, qrystartpos*1e-3, qryendpos*1e-3);
	    fflush(stdout);
	  }

	  if(smap_SmallDelConfirm > 0.0 && (0 < k-1 || k+1 < U)){
	    int nI = sI, nK = sK, nJ = sJ, nlastI = slastI, nlastK = slastK, nlastJ = slastJ;
	    int pk = k-1, nk = k;

	    /* For deletions in ref intervals < smap_SmallDelConfirm : try to include neighboring aligned intervals and see if confidence goes down */
	    if(pk > 0 && deltaY < smap_SmallDelConfirm && deltaY > deltaX){
	      pk--;
	      nlastI = p->sites1[pk];
	      nlastK = (resSD[0] > 0.0) ? p->sitesK1[pk] : 0;
	      nlastJ = p->sites2[pk];
	    }
	    if(nk+1 < U && deltaY < smap_SmallDelConfirm && deltaY > deltaX){
	      nk++;
	      nI = p->sites1[nk];
	      nK = (resSD[0] > 0.0) ? p->sitesK1[nk] : 0;
	      nJ = p->sites2[nk];
	    }

	    /* compute score for alternate boundaries */
	    int nqrystartidx = p->orientation ? Xshift + M + 1 - nlastJ : Xshift + nlastJ;
	    int nqryendidx   = p->orientation ? Xshift + M + 1 - nJ : Xshift + nJ;
	    double ndeltaX = fabs(origX[nqryendidx] - origX[nqrystartidx]);
	    double ndeltaY = Yc(origY, nI, nK) - Yc(origY, nlastI, nlastK);
	    double nSVsize = ndeltaY - ndeltaX;
	    double nSVsizeFix = nSVsize - SVsize;

	    double origBias = Bias, origPen = Pen, origGauss = Gauss, origPenSm = PenSm;
	    if(DEBUG) assert(ndeltaY >= deltaY);
	    double ndeltaY3 = nSVsizeFix > 0.0 ? deltaY + nSVsizeFix : deltaY; // WAS41 min(ndeltaY, deltaY);
	    double ndeltaX3 = nSVsizeFix < 0.0 ? deltaX - nSVsizeFix : deltaX; // WAS41 max(rres*0.500, ndeltaY3 - nSVsize);
	    SintDetailRA(ndeltaX3, ndeltaY3, nJ - nlastJ, nI-nK-nlastI, nJ, nI, nK, nlastK, Y, Bias,Pen,Gauss,PenSm, 
			   (VERB>=2 && Ymap->id == 10 && Xmap->id == 111) ? 1 : 0);
	    if(FRATE_FIX1 && smap_IndelConfBySize){
	      Gauss += ndeltaX * Frate;
	      Pen -= ndeltaX * Frate;
	    }
	    double noutscore = smap_IndelConfBySize ? Gauss : Bias + Pen + Gauss; // NOTE : excludes PenSm which is the misresolved penalty for the right end of the interval

	    if(noutscore > outscore){// NOTE : outscore values are typically -ve, so the largest value results in the smallest confidence value
	      if((SVsize < 0.0) == (nSVsize < 0.0) && fabs(SVsize - nSVsize) <= fabs((SVsize + nSVsize)*0.5) * smap_IndelRangeMinSVchange){/* keep original ranges, just update SVsizeFix */
		SVsizeFix = nSVsizeFix;

		if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
		  printf("\t Yid=%lld,Xid=%lld: outscore[%d]/L10= %0.4f -> %0.4f -> %0.4f (Bias= %0.4f, Gauss= %0.4f, Pen= %0.4f, PenSm= %0.4f): keeping original range with SVsize= %0.3f -> %0.3f(SVsizeFix=%0.3f)\n",
			 Ymap->id,Xmap->id, k,p->outscore[k]*Ilog10, outscore*Ilog10, noutscore*Ilog10, Bias, Gauss, Pen, PenSm, SVsize,nSVsize,SVsizeFix);
		  fflush(stdout);
		}
	      } else {/* switch to new ranges and update coordinates */
		SVsize = nSVsize;
		SVsizeFix = 0.0;

		/* update coordinates */
		sI = nI, sK = nK, sJ = nJ, slastI = nlastI, slastK = nlastK, slastJ = nlastJ;

		qrystartidx = p->orientation ? Xshift + M + 1 - nlastJ : Xshift + nlastJ;
		qryendidx   = p->orientation ? Xshift + M + 1 - nJ : Xshift + nJ;
		refstartidx    = Yshift + nlastI - nlastK;
		refendidx      = Yshift + nI;

		if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
		  printf("\t Yid=%lld,Xid=%lld: outscore[%d]/L10= %0.4f -> %0.4f -> %0.4f (Bias= %0.4f, Gauss= %0.4f, Pen= %0.4f, PenSm= %0.4f), ref= %0.3f..%0.3f -> %0.3f..%0.3f, qry= %0.3f..%0.3f -> %0.3f..%0.3f\n",
			 Ymap->id,Xmap->id, k,p->outscore[k]*Ilog10, outscore*Ilog10, noutscore*Ilog10, Bias, Gauss, Pen, PenSm, 
			 refstartpos, refendpos, origY[refstartidx]*1000.0*Yscale, origY[refendidx]*1000.0*Yscale,
			 qrystartpos, qryendpos, origX[qrystartidx]*100.0*Xscale, origX[qryendidx]*1000.0*Xscale);
		  printf("\t I=%d..%d -> %d..%d, K=%d..%d -> %d..%d, J=%d..%d -> %d..%d\n",
			 lastI,I,nlastI,nI,lastK,K,nlastK,nK,lastJ,J,nlastJ,nJ);
		  fflush(stdout);
		}

		qrystartpos = origX[qrystartidx]*1000.0*Xscale; /* QryStartPos */
		qryendpos   = origX[qryendidx  ]*1000.0*Xscale; /* QryEndPos */
		refstartpos = origY[refstartidx]*1000.0*Yscale; /* RefStartPos */
		refendpos   = origY[refendidx  ]*1000.0*Yscale; /* RefEndPos */

		/* include half of de-res interval in query interval, so query interval can be replaced with ref interval without introducing a length bias */
		origqsiz = p->orientation ? (qrystartpos - qryendpos) : (qryendpos - qrystartpos);
		origqrystartpos = qrystartpos;
		origqryendpos = qryendpos;
		origrefstartpos = (origY[Yshift + nlastI]*1000.0*Yscale + refstartpos)*0.5;
		origrefendpos = (refendpos + origY[Yshift + nI -nK]*1000.0*Yscale)*0.5;

		if(!p->orientation){
		  if(nlastK > 0)
		    qrystartpos -= 0.5*(origY[Yshift+nlastI] - origY[Yshift+nlastI-nlastK]) * 1000.0 * Yscale;
		  if(nK > 0)
		    qryendpos += 0.5*(origY[Yshift+nI] - origY[Yshift+nI-nK]) * 1000.0 * Yscale;
		} else {/* query runs from high to low positions along alignment */ 
		  if(nlastK > 0)
		    qrystartpos += 0.5*(origY[Yshift+nlastI] - origY[Yshift+nlastI-nlastK]) * 1000.0 * Yscale;
		  if(nK > 0)
		    qryendpos -= 0.5*(origY[Yshift+nI] - origY[Yshift+nI-nK]) * 1000.0 * Yscale;
		}
		if(DEBUG/* HERE >=2 */) assert(refstartpos <= refendpos);
	      }

	      outscore = noutscore;
	    }// if(noutscore > outscore)
	    else {
	      if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
		printf("\t Yid=%lld,Xid=%lld: outscore[%d]/L10= %0.4f -> %0.4f -> %0.4f (Bias= %0.4f, Gauss= %0.4f, Pen= %0.4f, PenSm= %0.4f), deltaX= %0.3f -> %0.3f, deltaY= %0.3f -> %0.3f :reverting\n",
		       Ymap->id,Xmap->id, k,p->outscore[k]*Ilog10, outscore*Ilog10, noutscore*Ilog10, Bias, Gauss, Pen, PenSm, deltaX, ndeltaX, deltaY, ndeltaY);
		printf("\t\t I=%d..%d -> %d..%d (reverting), K=%d..%d -> %d..%d (reverting), J=%d..%d -> %d..%d (reverting)\n",
		       lastI,I,nlastI,nI,lastK,K,nlastK,nK,lastJ,J,nlastJ,nJ);
		fflush(stdout);
	      }
	      // NEW37 : undo changes to nI,nK,nJ etc since score did NOT get worse
	      nI = I, nK = K, nJ = J, nlastI = lastI, nlastK = lastK, nlastJ = lastJ;
	      pk = k-1, nk = k;

	      ndeltaX = deltaX;
	      ndeltaY = deltaY;

	      Bias = origBias;
	      Gauss = origGauss;
	      Pen = origPen;
	      PenSm = origPenSm;
	    }

	    if(SMALL_CONFIRM_FIX){/* check if reference interval is bounded by deres'd labels (OR nearby labels) : if so expand to nearest non-deres'd (AND isolated) aligned label */
	      int nI2 = nI, nK2 = nK, nJ2 = nJ, nlastI2 = nlastI, nlastK2 = nlastK, nlastJ2 = nlastJ;
	      int pk2 = pk, nk2 = nk;

	      while(pk2 > 0){
		if(DEBUG) assert(0 <= nlastK2 && 1 < nlastI2 - nlastK2 && nlastI2 < N);
		if((nlastK2 <= 0 || origY[nlastI2] - origY[nlastI2 - nlastK2] <= smap_IndelRangeExpandRes) && 
		   Yc(origY,nlastI2,nlastK2) - origY[nlastI2 - nlastK2 - 1] > smap_IndelRangeExpandSep && origY[nlastI2 + 1] - origY[nlastI2] > smap_IndelRangeExpandSep)
		  break;

		pk2--;
		nlastI2 = p->sites1[pk2];
		nlastK2 = (resSD[0] > 0.0) ? p->sitesK1[pk2] : 0;
		nlastJ2 = p->sites2[pk2];
	      }
	      while(nk2 + 1 < U){
		if(DEBUG) assert(0 <= nK2 && 1 < nI2 - nK2 && nI2 < N);
		if((nK2 <= 0 || origY[nI2] - origY[nI2 - nK2] <= smap_IndelRangeExpandRes) 
		   && Yc(origY,nI2,nK2) - origY[nI2 - nK2 - 1] > smap_IndelRangeExpandSep/*NEW37*/ && origY[nI2 + 1] - origY[nI2] > smap_IndelRangeExpandSep/*NEW37*/)
		  break;

		nk2++;
		nI2 = p->sites1[nk2];
		nK2 = (resSD[0] > 0.0) ? p->sitesK1[nk2] : 0;
		nJ2 = p->sites2[nk2];
	      }
	      
	      int nqrystartidx2 = p->orientation ? Xshift + M + 1 - nlastJ2 : Xshift + nlastJ2;
	      int nqryendidx2   = p->orientation ? Xshift + M + 1 - nJ2 : Xshift + nJ2;
	      double ndeltaX2 = fabs(origX[nqryendidx2] - origX[nqrystartidx2]);
	      double ndeltaY2 = Yc(origY, nI2, nK2) - Yc(origY, nlastI2, nlastK2);
	      double nSVsize2 = ndeltaY2 - ndeltaX2;
	      double nSVsizeFix2 = nSVsize2 - SVsize;

	      if(nk2 > nk || pk2 < pk){
		double Bias2,Pen2,Gauss2,PenSm2;
		if(DEBUG) assert(ndeltaY2 >= ndeltaY);
		double ndeltaY3 = nSVsizeFix2 > 0.0 ? ndeltaY + nSVsizeFix2 : ndeltaY; // WAS41 min(ndeltaY2, ndeltaY);
		double ndeltaX3 = nSVsizeFix2 < 0.0 ? ndeltaX - nSVsizeFix2 : ndeltaX; // WAS41 max(rres*0.5, ndeltaY3 - nSVsize2);
		SintDetailRA(ndeltaX3, ndeltaY3, nJ2 - nlastJ2, nI2-nK2-nlastI2, nJ2, nI2, nK2, nlastK2, Y, Bias2,Pen2,Gauss2,PenSm2, 
			     (VERB>=2 && Ymap->id==10 && Xmap->id==111) ? 1 : 0);
		if(FRATE_FIX1 && smap_IndelConfBySize){
		  Gauss2 += ndeltaX3 * Frate;
		  Pen2 -= ndeltaX3 * Frate;
		}
		double noutscore2 = smap_IndelConfBySize ? Gauss2 : Bias2 + Pen2 + Gauss2; // NOTE : excludes PenSm which is the misresolved penalty for the right end of the interval

		if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
		  int nrefstartidx2    = Yshift + nlastI2 - nlastK2;
		  int nrefendidx2      = Yshift + nI2;

		  printf("\t Yid=%lld,Xid=%lld: outscore[%d]/L10= %0.4f -> %0.4f -> %0.4f (Bias= %0.4f -> %0.4f, Gauss= %0.4f -> %0.4f, Pen= %0.4f -> %0.4f, PenSm= %0.4f -> %0.4f)\n\t    ref= %0.3f..%0.3f -> %0.3f..%0.3f, qry= %0.3f..%0.3f -> %0.3f..%0.3f, SVsize= %0.3f -> %0.3f, deltaY= %0.3f -> %0.3f -> %0.3f, deltaX= %0.3f -> %0.3f -> %0.3f\n\t    I=%d..%d -> %d..%d, K=%d..%d -> %d..%d, J=%d..%d -> %d..%d\n",
			 Ymap->id,Xmap->id, k,p->outscore[k]*Ilog10, outscore*Ilog10, noutscore2*Ilog10, Bias, Bias2, Gauss, Gauss2, Pen, Pen2, PenSm, PenSm2,
			 refstartpos*1e-3, refendpos*1e-3, origY[nrefstartidx2]*Yscale, origY[nrefendidx2]*Yscale,
			 qrystartpos*1e-3, qryendpos*1e-3, origX[nqrystartidx2]*Xscale, origX[nqryendidx2]*Xscale, SVsize, nSVsize2, ndeltaY, ndeltaY2, ndeltaY3, ndeltaX, ndeltaX2, ndeltaX3,
			 nlastI,nI,nlastI2,nI2,nlastK,nK,nlastK2,nK2,nlastJ,nJ,nlastJ2,nJ2);
		  fflush(stdout);
		}

		if(noutscore2 < outscore /* HERE HERE && fabs(nSVsize2) < 0.700*/){// Why is this still needed ? To avoid a doubling of FP calls for small deletions 300-500bp
		  double origBias2 = Bias2, origPen2 = Pen2, origGauss2 = Gauss2, origPenSm2 = PenSm2;

		  SintDetailRA(ndeltaX2, ndeltaY2, nJ2 - nlastJ2, nI2-nK2-nlastI2, nJ2, nI2, nK2, nlastK2, Y, Bias2,Pen2,Gauss2,PenSm2, 
			       (VERB>=2 && Ymap->id==10 && Xmap->id==111) ? 1 : 0);
		  if(FRATE_FIX1 && smap_IndelConfBySize){
		    Gauss2 += ndeltaX2 * Frate;
		    Pen2 -= ndeltaX2 * Frate;
		  }
		  double noutscore3 = smap_IndelConfBySize ? Gauss2 : Bias2 + Pen2 + Gauss2; // NOTE : excludes PenSm which is the misresolved penalty for the right end of the interval	 

		  if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
		    int nrefstartidx2    = Yshift + nlastI2 - nlastK2;
		    int nrefendidx2      = Yshift + nI2;

		    printf("\t Yid=%lld,Xid=%lld: outscore[%d]/L10= %0.4f -> %0.4f -> %0.4f (Bias= %0.4f -> %0.4f -> %0.4f, Gauss= %0.4f -> %0.4f -> %0.4f, Pen= %0.4f -> %0.4f -> %0.4f, PenSm= %0.4f -> %0.4f -> %0.4f)\n\t    ref= %0.3f..%0.3f -> %0.3f..%0.3f, qry= %0.3f..%0.3f -> %0.3f..%0.3f, SVsize= %0.3f -> %0.3f, deltaY= %0.3f -> %0.3f, deltaX= %0.3f -> %0.3f\n",
			   Ymap->id,Xmap->id, k, outscore*Ilog10, noutscore2*Ilog10, noutscore3*Ilog10, Bias, origBias2, Bias2, Gauss, origGauss2, Gauss2, Pen, origPen2, Pen2, PenSm, origPenSm2, PenSm2,
			   refstartpos*1e-3, refendpos*1e-3, origY[nrefstartidx2]*Yscale, origY[nrefendidx2]*Yscale,
			   qrystartpos*1e-3, qryendpos*1e-3, origX[nqrystartidx2]*Xscale, origX[nqryendidx2]*Xscale, SVsize, nSVsize2, ndeltaY, ndeltaY2, ndeltaX, ndeltaX2);
		    fflush(stdout);
		  }
		  noutscore2 = max(noutscore3, noutscore2);
		}

		/* NEW380 : If SVsize sign is reversed OR SVsize more than doubles, skip this indel since it will be rediscovered from the other larger indel */
                /* (alternately revert to original range and SVsize) */
		if((SVsize < 0.0) != (nSVsize2 < 0.0) || fabs(nSVsize2) > 2.0 * fabs(SVsize)){
#if 1
		  if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
		    printf("\t Yid=%lld,Xid=%lld: skipping Indel sinze larger Indel is nearby (SVsize= %0.3f -> %0.3f kb)\n",
			   Ymap->id,Xmap->id, SVsize, nSVsize2);
		    fflush(stdout);
		  }		  
		  continue;// skip this indel
#else
		  if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
		    printf("\t Yid=%lld,Xid=%lld: reverting to original range and SVsize= %0.3f\n",
			   Ymap->id,Xmap->id, SVsize);
		    fflush(stdout);
		  }		  
		  SVsizeFix = 0.0;/* revert to original range and size */
#endif
		} else if((SVsize < 0.0) == (nSVsize2 < 0.0) && fabs(SVsize - nSVsize2) <= fabs((SVsize + nSVsize2)*0.5) * smap_IndelRangeMinSVchange){/* keep original ranges, just update SVsizeFix */
		  SVsizeFix = nSVsizeFix2;

		  if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
		    printf("\t Yid=%lld,Xid=%lld: keeping original range with SVsize= %0.3f -> %0.3f(SVsizeFix= %0.3f) kb\n",
			   Ymap->id,Xmap->id, SVsize, nSVsize2, SVsizeFix);
		    fflush(stdout);
		  }
		} else {/* switch to new ranges and update coordinates */
		  SVsize = nSVsize2;
		  SVsizeFix = 0.0;

		  /* update coordinates */
		  sI = nI2, sK = nK2, sJ = nJ2, slastI = nlastI2, slastK = nlastK2, slastJ = nlastJ2;

		  qrystartidx = p->orientation ? Xshift + M + 1 - nlastJ2 : Xshift + nlastJ2;
		  qryendidx   = p->orientation ? Xshift + M + 1 - nJ2 : Xshift + nJ2;
		  refstartidx    = Yshift + nlastI2 - nlastK2;
		  refendidx      = Yshift + nI2;

		  qrystartpos = origX[qrystartidx]*1000.0*Xscale; /* QryStartPos */
		  qryendpos   = origX[qryendidx  ]*1000.0*Xscale; /* QryEndPos */
		  refstartpos = origY[refstartidx]*1000.0*Yscale; /* RefStartPos */
		  refendpos   = origY[refendidx  ]*1000.0*Yscale; /* RefEndPos */

		  /* include half of de-res interval in query interval, so query interval can be replaced with ref interval without introducing a length bias */
		  origqsiz = p->orientation ? (qrystartpos - qryendpos) : (qryendpos - qrystartpos);
		  origqrystartpos = qrystartpos;
		  origqryendpos = qryendpos;
		  origrefstartpos = (origY[Yshift + nlastI2]*1000.0*Yscale + refstartpos)*0.5;
		  origrefendpos = (refendpos + origY[Yshift + nI2 -nK2]*1000.0*Yscale)*0.5;

		  if(!p->orientation){
		    if(nlastK2 > 0)
		      qrystartpos -= 0.5*(origY[Yshift + nlastI2] - origY[Yshift + nlastI2 - nlastK2]) * 1000.0 * Yscale;
		    if(nK2 > 0)
		      qryendpos += 0.5*(origY[Yshift + nI2] - origY[Yshift + nI2 - nK2]) * 1000.0 * Yscale;
		  } else {/* query runs from high to low positions along alignment */ 
		    if(nlastK2 > 0)
		      qrystartpos += 0.5*(origY[Yshift + nlastI2] - origY[Yshift + nlastI2 - nlastK2]) * 1000.0 * Yscale;
		    if(nK2 > 0)
		      qryendpos -= 0.5*(origY[Yshift + nI2] - origY[Yshift + nI2 - nK2]) * 1000.0 * Yscale;
		  }
		  if(DEBUG/* HERE >=2 */) assert(refstartpos <= refendpos);
		}

		outscore = noutscore2;// NEW37 : accept new score regardless of whether it is better or worse
	      } // if(nk2 > nk || pk2 < pk)
	      else if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
		double ndeltaY3 = nSVsizeFix2 > 0.0 ? ndeltaY + nSVsizeFix2 : ndeltaY; // WAS41 min(ndeltaY2, ndeltaY);
		double ndeltaX3 = nSVsizeFix2 < 0.0 ? ndeltaX - nSVsizeFix2 : ndeltaX; // WAS41 max(rres*0.5, ndeltaY3 - nSVsize2);

		int nrefstartidx2    = Yshift + nlastI2 - nlastK2;
		int nrefendidx2      = Yshift + nI2;

		printf("\t Yid=%lld,Xid=%lld:k=%d: ref= %0.3f..%0.3f -> %0.3f..%0.3f, qry= %0.3f..%0.3f -> %0.3f..%0.3f, SVsize= %0.3f -> %0.3f, deltaY= %0.3f -> %0.3f -> %0.3f, deltaX= %0.3f -> %0.3f -> %0.3f:reverting\n",
		       Ymap->id,Xmap->id, k,
		       refstartpos*1e-3, refendpos*1e-3, origY[nrefstartidx2]*Yscale, origY[nrefendidx2]*Yscale,
		       qrystartpos*1e-3, qryendpos*1e-3, origX[nqrystartidx2]*Xscale, origX[nqryendidx2]*Xscale, SVsize, nSVsize2, ndeltaY, ndeltaY2, ndeltaY3, ndeltaX, ndeltaX2, ndeltaX3);
		fflush(stdout);
	      }
	    } // if(SMALL_CONFIRM_FIX)
          } // if(smap_SmallDelConfirm > 0.0)

	  centerConfidence = -outscore*Ilog10;
	  outlierConfidence = min(centerConfidence,p->logPV);
	  if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
	    printf("\t outscore[%d]= %0.4f -> %0.4f (Bias= %0.4f,Pen= %0.4f, Gauss= %0.4f, PenSm= %0.4f): PenoutlierConfidence= %0.2f, smap_min_conf= %0.2f\n",
		   k, p->outscore[k], outscore, Bias, Pen, Gauss, PenSm, outlierConfidence, smap_min_conf);
	    fflush(stdout);
	  }
	} else {
	  centerConfidence = (-p->outscore[k]) * Ilog10;
	  outlierConfidence = min(centerConfidence,p->logPV);
	  if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
	    printf("\t outscore[%d]= %0.4f: outlierConfidence= %0.2f, smap_min_conf= %0.2f\n",k, p->outscore[k], outlierConfidence, smap_min_conf);
	    fflush(stdout);
	  }
	}
	if(DEBUG/* HERE >=2 */) assert(isfinite(outlierConfidence));

	// NOTE : .indel includes indels below smap_min_conf as long as center confidence is above smap_min_conf
	if(outlierConfidence <= smap_min_conf)
	  continue;

	/* compare with end outlier to the left */
	extern double SmRA(int J, int I, int K, FLOAT *Yref);
	double Leftscore = SmRA(JL,IL,KL,origY);
	for(int t = 1; t < k; t++)
	  Leftscore += p->iscore[t];
	if(OUTLIER_CONF_FIX && p->Lend >= -1 && k >= OUTLIER_CONF_FIX)
	  Leftscore += p->iscore[0];
	//if((Leftscore - EndOutlier)*Ilog10 < IndelConfidence) //old way
	//IndelConfidence = (Leftscore - EndOutlier)*Ilog10;
	double LeftConfidence = Leftscore*Ilog10;
	if(DEBUG/* HERE >=2 */) assert(isfinite(Leftscore));

	/* compare with end outlier to the right */
	double Rightscore = PenSm;// Sm(JR,IR,KR,origY)
	for(int t = k+1; t < U; t++)
	  Rightscore += p->iscore[t];
	if(OUTLIER_CONF_FIX && p->Rend >= -1 && U-k >= OUTLIER_CONF_FIX)
	  Rightscore += p->iscore[U];
	//if((Rightscore - EndOutlier)*Ilog10 < IndelConfidence) //old way
	//IndelConfidence = (Rightscore - EndOutlier)*Ilog10;
	double RightConfidence = Rightscore*Ilog10;
	if(DEBUG/* HERE >=2 */) assert(isfinite(Rightscore));

	//new way: fn in structuralVariation.{h,cpp}
	double IndelConfidence = makeConfidence(outlierConfidence, LeftConfidence, RightConfidence);
	if(DEBUG/* HERE >=2 */) assert(isfinite(IndelConfidence));

	if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
	  printf("Yid=%lld,Xid=%lld,or=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f:CenterConf=%0.6f,LeftConf=%0.6f,RightConf=%0.6f,logPV=%0.6f:Confidence=%0.6f\n",
		 Ymap->id,Xmap->id,p->orientation,I,K,J,lastI,lastK,lastJ,p->outscore[k],p->iscore[k],centerConfidence,LeftConfidence,RightConfidence,p->logPV, IndelConfidence);
	  printf("\t iscore[0]= %0.6f, iscore[U]= %0.6f, Lend=%d, Rend=%d\n",p->iscore[0],p->iscore[U],p->Lend,p->Rend);
	  if(VERB>=2){
	    printf("Leftscore components:\n");
	    for(int t = 1; t < k; t++)
	      printf("t=%d,k=%d:iscore[t]=%0.6f\n",t,k,p->iscore[t]);
	    if(OUTLIER_CONF_FIX && p->Lend >= -1 && k >= OUTLIER_CONF_FIX)
	      printf("Lend=%d:iscore[0]=%0.6f\n",p->Lend,p->iscore[0]);
	  }
	  if(VERB>=2){
	    printf("Rightscore components:\n");
	    for(int t = k+1; t < U; t++)
	      printf("t=%d,U=%d:iscore[t]=%0.6f\n",t,U,p->iscore[t]);
	    if(OUTLIER_CONF_FIX && p->Rend >= -1 && U-k >= OUTLIER_CONF_FIX)
	      printf("Rend=%d,U=%d:iscore[U]=%0.6f\n",p->Rend,U,p->iscore[U]);
	  }
	  fflush(stdout);
	}

	if(InDelMinConf > 0.0 && max(RightConfidence,LeftConfidence) < InDelMinConf)
	  continue;

	/* make sure there is no other alignment of this Query (Xmap) overlapping this indel location for which the
	   basic alignment has a better logPV : If they exist, the overlap region was previously computed in the alignment as sites p->suppress[1..M] and are non-zero in range p->suppress[p->start..p->end] */
	int start,end;
	if(p->orientation){
	  if(DEBUG) assert(qryendidx <= qrystartidx);
	  start = max(qryendidx,p->start);
	  end = min(qrystartidx,p->end);
	} else {
	  if(DEBUG) assert(qryendidx >= qrystartidx);
	  start = max(qrystartidx,p->start);
	  end = min(qryendidx,p->end);
	}
	if(start < end){
	  int i;
	  for(i = start; i <= end; i++)
	    if(p->suppress[i])
	      break;
	  if(i <= end){
	    if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
	      int j = i+1;
	      for(; j <= end; j++)
		if(!p->suppress[i])
		  break;
	      printf("Yid=%lld,Xid=%lld,or=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f:query=%d..%d(%0.3f..%0.3f),ref=%0.3f..%0.3f:outlier skipped due to higher scoring alignment overlap at i=%d..%d\n",
		     Ymap->id,Xmap->id,p->orientation,I,K,J,lastI,lastK,lastJ,p->outscore[k],p->iscore[k],qrystartidx,qryendidx,qrystartpos,qryendpos,refstartpos,refendpos,i,j-1);
	      fflush(stdout);
	    }
	    continue;
	  }
	}

	if(sI == psI && sK == psK && sJ == psJ && slastI == pslastI && slastK == pslastK && slastJ == pslastJ){
	  if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
	    printf("Yid=%lld,Xid=%lld,or=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f:query=%0.3f..%0.3f,ref=%0.3f..%0.3f:duplicate indel found(Xshift=%d,Yshift=%d)\n",
		   Ymap->id,Xmap->id,p->orientation,I,K,J,lastI,lastK,lastJ,p->outscore[k],p->iscore[k],qrystartpos,qryendpos,refstartpos,refendpos,Xshift,Yshift);
	    fflush(stdout);
	  }
	  continue;
	}
	double delY = (SVsizeFix > 0.0) ? SVsizeFix * 500.0 /* WAS380 0.5 */ : 0.0;/* In bp ! */
	double delX = (SVsizeFix < 0.0) ? (p->orientation ? SVsizeFix : -SVsizeFix) * 500.0 /* WAS380 0.5 */ : 0.0;/* In bp ! */

	if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE){
	  if(p->start > 0){
	    printf("Yid=%lld,Xid=%lld,or=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f:query=%0.3f..%0.3f,ref=%0.3f..%0.3f,Fix=%0.3f:indel found(Xshift=%d,Yshift=%d),qovr=%d..%d\n",
		   Ymap->id,Xmap->id,p->orientation,I,K,J,lastI,lastK,lastJ,p->outscore[k],p->iscore[k],qrystartpos*1e-3,qryendpos*1e-3,refstartpos*1e-3,refendpos*1e-3,SVsizeFix,Xshift,Yshift,p->start,p->end);
	  } else
	    printf("Yid=%lld,Xid=%lld,or=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f:query=%0.3f..%0.3f,ref=%0.3f..%0.3f,Fix=%0.3f:indel found(Xshift=%d,Yshift=%d)\n",
		   Ymap->id,Xmap->id,p->orientation,I,K,J,lastI,lastK,lastJ,p->outscore[k],p->iscore[k],qrystartpos*1e-3,qryendpos*1e-3,refstartpos*1e-3,refendpos*1e-3,SVsizeFix,Xshift,Yshift);
	  printf("    sI= %d..%d -> %d..%d, sK= %d..%d -> %d..%d, sJ= %d..%d -> %d..%d (previous -> current indel) : delY= %0.1f, delX= %0.1f\n",
		 pslastI,psI,slastI,sI,pslastK,psK,slastK,sK,pslastJ,psJ,slastJ,sJ, delY, delX);
	  fflush(stdout);
	}

	psI = sI, psK = sK, psJ = sJ, pslastI = slastI, pslastK = slastK, pslastJ = slastJ;


	SmapCnt++;
	bool isdeletion = (refendpos-refstartpos > fabs(qryendpos-qrystartpos));
	if(fpSV)
	  fprintf(fpSV,"%llu\t%lld\t%lld\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%c\t%0.2f\t%s%dI%dD%d\t%lu\t%lu\n",
		(unsigned long long)SmapCnt,    // SmapEntryID
		Xmap->id,   // QryContigID
		Ymap->id,   // RefContigID1
		Ymap->id,   // RefContigID2
		  qrystartpos - delX,
		  qryendpos + delX,
		refstartpos - delY,
		refendpos + delY,
		p->orientation ? '-' : '+', // Orientation
		IndelConfidence, 
		isdeletion ? "deletion_" : "insertion_", // Type
		  (int)fabs(fabs(qryendpos-qrystartpos+2.0*delX) - (refendpos-refstartpos+2.0*delY)),// indel size in bp
		sJ - slastJ - 1 /* + (lastK > 0 ? 1 : 0) + (K > 0 ? 1 : 0)*/, // inserted sites
		(sI-sK) - slastI /* I - (lastI-lastK)*/ - 1, // deleted sites
		aligncnt, aligncnt); // XmapID1, XmapID2

	// NOTE : .indel includes indels below smap_min_conf as long as center confidence is above smap_min_conf
	if(IndelConfidence <= smap_min_conf)
	  continue;

	if( dosvdetect ) {
	  //create structuralVariation object to store this data
	  //structuralVariation(int ind1, int ind2, double qsize, double rsize, double qstart, double qstop, int refid1, int refid2, double rstart, double rstop, const char* svtd, sv_type_t svt, int qsi=-1, int qei=-1, int rsi=-1, int rei=-1) {
	  if( numsv+1 > maxsv ) //must increase size of sv_array
	    maxsv = growArray(sv_array, maxsv, 256);

	  double indelLeft = min(origqrystartpos,origqryendpos) - min(qrystartpos,qryendpos) ;
	  double indelRight = max(qrystartpos,qryendpos) - max(origqrystartpos,origqryendpos);

	  if(INDEL_TRACE && Ymap->id == REF_TRACE && Xmap->id == MAP_TRACE)
	    printf("Indel %2i: %9s: chr %lld,qryid=%lld, qry=%0.3f .. %0.3f, ref=%0.3f .. %0.3f: Conf= %5.2f: outlierConf= %0.2f, Leftconf= %.2f, Rightconf= %0.2f,LeftRes=%0.3f,RightRes=%0.3f\n", 
		   numsv+1, isdeletion ? "deletion" : "insertion", Ymap->id, Xmap->id, (qrystartpos-delX)*1e-3, (qryendpos+delX)*1e-3,(refstartpos-delY)*1e-3, (refendpos+delY)*1e-3, 
		   IndelConfidence, outlierConfidence, LeftConfidence, RightConfidence, indelLeft*1e-3,indelRight*1e-3); //debug

	  sv_array[numsv] = structuralVariation(numxmapentries, numxmapentries, fabs(qryendpos-qrystartpos) - indelLeft - indelRight + delX*2.0, (refendpos-refstartpos) -indelLeft - indelRight + delY * 2.0,
						min(qrystartpos-delX, qryendpos+delX) + indelLeft,  max(qrystartpos-delX, qryendpos+delX) - indelRight, Ymap->id, Ymap->id,
						origrefstartpos - delY, origrefendpos + delY, 0, isdeletion ? deletion : insertion,
						min(qrystartidx, qryendidx), max(qrystartidx, qryendidx), refstartidx, refendidx, outlierConfidence, LeftConfidence, RightConfidence, sv_indel, aligncnt, aligncnt);
	  numsv++;
	  nindel++;
	}
      } /* lastJ, J loop */
      
      if(InDelEnds){/* check for end outliers to add to .indel file */
	if(p->Lend <= -2){
	  SmapCnt++;
	  if(fpSV)
	    fprintf(fpSV,"%llu\t%lld\t%lld\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%c\t%0.2f\t%s\t%lu\t%lu\n",
		  (unsigned long long)SmapCnt,    // SmapEntryID
		  Xmap->id,   // QryContigID
		  Ymap->id,   // RefContigID1
		  Ymap->id,   // RefContigID2
		  p->orientation ? qryLen : 0.0, // Left End of Query region 
		  qrystartpos, // Right End of Query region
		  refstartpos, // Right End of Reference region
		  -1.0,         // Not Defined
		  p->orientation ? '-' : '+', // Orientation
		  -1.0, // Not computed
		  "end", // Type
		  aligncnt, aligncnt); // XmapID1, XmapID2
	}
	if(p->Rend <= -2){
	  SmapCnt++;
	  if(fpSV)
	    fprintf(fpSV,"%llu\t%lld\t%lld\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%c\t%0.2f\t%s\t%lu\t%lu\n",
		  (unsigned long long)SmapCnt,    // SmapEntryID
		  Xmap->id,   // QryContigID
		  Ymap->id,   // RefContigID1
		  Ymap->id,   // RefContigID2
		  qryendpos,   // Left End of Query region 
		  p->orientation ? 0.0 : qryLen, // Right End of Query region
		  refendpos, // Left End of Reference region
		  -1.0,         // Not Defined
		  p->orientation ? '-' : '+', // Orientation
		  -1.0, // Not computed
		  "end", // Type
		  aligncnt, aligncnt); // XmapID1, XmapID2
	}
      }
      if(VERB && end - start <= 10000)
	fflush(stdout);
    } /* if(Indel) */

    if( xmapentries != NULL /* (dosvdetect || doscaffold) */ &&  !SecondLine) { /* see parameters.{h,cpp} */
      int orignsite = (Xmap->origmap ? Xmap->origmap->numsite[0] : Xmap->numsite[0]); //need origmap bc of pairsplit (?)
      double qlen = origX[orignsite+1];
      //remember, orientation is backwards... true is reverse
      int startua = p->orientation ? qryendidx - 1 : qrystartidx - 1;
      int endua   = p->orientation ? orignsite - qrystartidx : orignsite - qryendidx;
      int numrefsites = refendidx - refstartidx + 1; // WAS44 RI - LI + 1; //+1 bc both ends are included
      if(DEBUG>=1+RELEASE) assert(numrefsites == RI - LI + LK + 1);
      if(VERB>=2 || (DEBUG && !Xmap->origmap && !(orignsite == XXmap[p->mapid2]->numsite[0]))){
	printf("xmapentries[%d]:aligncnt=%lu,refid=%lld,qryid=%lld,or=%d:qrystartpos=%0.1f,qryendpos=%0.1f,refstartpos=%0.1f,refendpos=%0.1f,logPV=%0.2f,orignsite=%d,N=%d,Yscale= %0.8f, Yshift=%d, RefSplit=%d\n",
	       numxmapentries,aligncnt,Ymap->id,Xmap->id,p->orientation,qrystartpos,qryendpos,refstartpos,refendpos,p->logPV,orignsite,N, Yscale,Yshift,RefSplit);
	fflush(stdout);
	if(DEBUG && !Xmap->origmap) assert(orignsite == XXmap[p->mapid2]->numsite[0]);
      }
      if(DEBUG>=1+RELEASE && RefSplit) assert(Yscale == 1.0 && Yshift == 0);// If this is not true during SVstage need to add Yscale and/or Yshift in output_smap.cpp

      xmapentries[numxmapentries] = new xmapEntry(p, aligncnt, Xmap->id, Ymap->id, qrystartpos, qryendpos, qrystartidx, qryendidx, refstartpos, refendpos, refstartidx, refendidx, !bool(p->orientation), p->logPV, startua, endua, qlen, numrefsites, orignsite, N); //arry index is one less than xmap id
      xmapentries[numxmapentries] -> RefMap = Ymap;// Used by Small Inversions
      xmapentries[numxmapentries] -> QryMap = Xmap;// Used by Small Inversions
      //array refsites is allocated in xmapEntry constructor and filled here
      //printf("align %i: %i\n", aligncnt, numrefsites);
      for(int i= refstartidx /* WAS130 Yshift+LI*/; i <= refendidx/*WAS130 Yshift+RI */; i++) { //include site at Yshift+RI: this is RefEndPos
	xmapentries[numxmapentries]->refsites[i-refstartidx/* Yshift-LI*/] = origY[i]*1000.*Yscale; //refsites is allocated in xmapEntry constructor
	//printf("%2i  %f\n", i, xmapentries[aligncnt-1]->refsites[i-Yshift-LI]); //origY[i]*1000.*Yscale);
      }
      if(p->orientation){
	if(DEBUG) assert(qrystartidx >= qryendidx);
	if(DEBUG) assert(qrystartpos >= qryendpos);
	for(int t = qryendidx; t <= qrystartidx; t++){
	  xmapentries[numxmapentries]->qrysites[t - qryendidx] = origX[t] * 1000.0;
	  if(DEBUG>=2) assert(qryendpos <= origX[t] * 1000.0  && origX[t] * 1000.0 <= qrystartpos);
	}
      } else {
	if(DEBUG) assert(qrystartidx <= qryendidx);
	if(DEBUG) assert(qrystartpos <= qryendpos);
	for(int t = qrystartidx; t <= qryendidx; t++){
	  xmapentries[numxmapentries]->qrysites[t - qrystartidx] = origX[t] * 1000.0;
	  if(DEBUG>=2) assert(qrystartpos <= origX[t] * 1000.0  && origX[t] * 1000.0 <= qryendpos);
	}
      }

      if(PairSplit && p->align2 && p->align2->logPV > p->logPV - log(10.0)*Palindromic_Delta)
	xmapentries[numxmapentries]->Palindromic = true;

      numxmapentries++;
      if(DEBUG) assert(0 <= numxmapentries && (size_t)numxmapentries <= end - start); 
      if(DEBUG) assert(aligncnt == (size_t)numxmapentries);
    }

    /* Output Alignment Cigar */
    int lastI = p->sites1[0];
    int lastJ = p->sites2[0];
    int matchcnt = 1;

    const bool DEBUG_CIGAR = false; //if true, output all aligned sites (1 row per match)
    if (DEBUG_CIGAR) {
      printf("==============================================================================\n");
      printf("XmapEntryID=%lu\tQryContigID=%lld\tRefContigID=%lld\n", aligncnt, Xmap->id, Ymap->id);
      printf("RefLab\tRefRes\tQryLab\n"); //"I\tK\tJ\n");
      int K = (resSD[0] > 0.0) ? p->sitesK1[0] : 0;
      int origJ = p->orientation ? Xshift + M + 1 - lastJ : Xshift + lastJ;
      printf("%i\t%i\t%i\n", Yshift + lastI, K, origJ); //first match is not in loop below
    } // if DEBUG_CIGAR

    for(int k = 1; k < U; k++){
      int I = p->sites1[k];
      int K = (resSD[0] > 0.0) ? p->sitesK1[k] : 0;
      if(DEBUG) assert(1 <= I-K && I <= N);
      int J = p->sites2[k];

      if(DEBUG && !(I-K > lastI && J > lastJ)){
	printf("QryID=%lld,RefID=%lld,or=%d:qrystart=%0.3f,qryend=%0.3f,refstart=%0.3f,refend=%0.3f,logPV=%0.2f:k=%d/%d:lastI=%d,lastJ=%d,I=%d,K=%d,J=%d\n",
	       Xmap->id,Ymap->id,p->orientation,qrystartpos,qryendpos,refstartpos,refendpos,p->logPV,k,U,lastI,lastJ,I,K,J);
	for(int t = 0; t < U; t++)
	  printf("t=%d:sites1[t]=%d,sites2[t]=%d\n",t,p->sites1[t],p->sites2[t]);
	fflush(stdout);
	assert(I-K > lastI && J > lastJ);
      }

      if (DEBUG_CIGAR) {      
	/* locations on complete maps origX[],origY[] */
	int origI = Yshift + I;
	int origJ = p->orientation ? Xshift + M + 1 - J : Xshift + J;
	printf("%i\t%i\t%i\n", origI, K, origJ); //origK = K
      } // if DEBUG_CIGAR

      //      printf("i=%d:mapid=%d(id=%d):k=%d/%d:I=%d,%d,K=%d,%d,J=%d,%d,matchcnt=%d\n",i,p->mapid2, XXmap[p->mapid2]->id, k,K,lastI,I,lastK,K,lastJ,J,matchcnt);
      //      fflush(stdout);

      if(K==0 && I == lastI + 1 && J == lastJ + 1){
	matchcnt++;
      } else {
	/* output match so far */
	if(matchcnt > 0)
	  fprintf(fp,"%dM",matchcnt);
	if(I == lastI + 1)
	  fprintf(fp,"%dI",J - lastJ -1);
	else if(J == lastJ + 1)
	  fprintf(fp,"%dD", I - lastI -1);
	else /* both an insertion and a deletion */
	  fprintf(fp,"%dI%dD", J - lastJ - 1, I - lastI - 1);
	matchcnt = 1;
      }
      lastI = I;
      lastJ = J;
    }
    if(matchcnt > 0)
      fprintf(fp,"%dM", matchcnt);
    /* done with cigar */

    fprintf(fp,"\t%0.1f\t%0.1f\t%d\t",qryLen,refLen, (usecolor==2) ? 2 : 1);
    
    /* Output Alignment doublet string */
    for(int k = 0; k < U; k++){
      int I = p->sites1[k];
      int K = (resSD[0] > 0.0) ? p->sitesK1[k] : 0;
      if(DEBUG) assert(1 <= I-K && I <= N);
      int J = p->sites2[k];
      
      int origI = Yshift + I;
      int origJ = p->orientation ? Xshift + M + 1 - J : Xshift + J;
      if(DEBUG && !(origI-K >= 1 && origI <= origN)){
	printf("i=%lu/(%lu..%lu):rid=%d(id=%lld,orig=%lld,Yshift=%d),mid=%d(id=%lld,orig=%lld,Xshift=%d),or=%d:k=%d/%d:I=%d,K=%d,J=%d,origI=%d,origJ=%d,origN=%d,origM=%d\n",
	       i,start,end,rid,Ymap->id,origYmap->id,Yshift,mid,Xmap->id,origXmap->id,Xshift,p->orientation,k,U,I,K,J,origI,origJ,origN,origM);
	fflush(stdout);
	fflush(fp);
	assert(origI-K >= 1 && origI <= origN);
      }
      if(DEBUG) assert(origJ >= 1 && origJ <= origM);
      if(K>0)
	fprintf(fp,"(%d,%d)",origI-K,origJ);
      fprintf(fp,"(%d,%d)",origI,origJ);
    }

    if(SecondBest && SinglelineXmap){
      Calign *p2 = p->align2;
      if(p2){
	int U2 = p2->numpairs;
	int RI2 = p2->sites1[U2-1];
	//      register int LK2 = p2->sitesK1[U2-1];
	int RJ2 = p2->sites2[U2-1];
	fprintf(fp,"\t%0.1f\t%0.1f\t%c\t%0.2f\t%d",
		(p2->orientation ? origX[Xshift + M + 1 - RJ2] : origX[Xshift + RJ2])*1000.0*Xscale,/* QryEndPos2 */
		origY[Yshift + RI2]*1000.0*Yscale,/* RefEndPos2 */
		p2->orientation ? '-' : '+', /* Orientation2 */
		p2->logPV,/*zscore(pow(10.0,-max(0.0,p->logPV2)))*/ /* Confidence2 */
		p2->numpairs /* NumMatch2 */
		);
      } else
	fprintf(fp,"\t%0.1f\t%0.1f\t%c\t%0.2f\t%d",
		0.0,/* QryEndPos2 */
		0.0,/* RefEndPos2 */
		'+', /* Orientation2 */
		0.0,/*zscore(pow(10.0,-max(0.0,p->logPV2)))*/ /* Confidence2 */
		0 /* NumMatch2 */
		);
    }

    if(SecondBest && !SinglelineXmap){
      if(VERB>=2){
	printf("i=%lu,aligncnt=%lu:Lalignment[i]=%p,Lalignment[i]->align2=%p,SecondLine=%d\n",i,aligncnt,Lalignment[i],Lalignment[i]->align2,SecondLine);
	fflush(stdout);
      }
      if(SecondLine){
	fprintf(fp,"\t%d", -1);
      } else {
	if(DEBUG>=2) assert(p == Lalignment[i]);
	size_t t = -1;
	if(p->align2){/* locate 2nd best alignment in Lalignment[] and outputs its XmapEntry value (==chimpair)*/
	  for(t = start; t < end; t++)
	    if(Lalignment[t] == p->align2)
	      break;
	  if(DEBUG) assert(start <= t && t < end);
	  assert(Lalignment[t]->align2 == NULL);
	}
	fprintf(fp,"\t%lu", p->align2 ? (size_t)abs(Lalignment[t]->chimpair) : 0);
      }
    }

    if(!SecondBest && XMAPmapWT)
      fprintf(fp,"\t%0.6e", p->mapWT);

    fprintf(fp,"\n");
  }
  
  if(VERB){
    printf("Generated %s with %lu match groups\n",filename, aligncnt);
    fflush(stdout);
  }
  if(VERB && fpSV){
    printf("Generated %s with %lu indels\n",filenameSV, SmapCnt);
    fflush(stdout);
  }

  if(VERB>=2){
    printf("RintervalCnt=%d,QintervalCnt=%d\n",RintervalCnt,QintervalCnt);
    fflush(stdout);
  }

  /* merge Rinterval[0..RintervalCnt-1] and Qinterval[0..QintervalCnt-1] to report total and unique interval size totals */
  if(RintervalCnt > 1)
    qsort(Rinterval,RintervalCnt, sizeof(Cinterval), (intcmp*)IntervalIdSiteInc);
  if(QintervalCnt > 1)
    qsort(Qinterval,QintervalCnt, sizeof(Cinterval), (intcmp*)IntervalIdSiteInc);

  int i = 0;
  FLOAT Rtotlen = (RintervalCnt <= 0) ? 0.0 : Rinterval[0].right - Rinterval[0].left;
  FLOAT Runiqlen = 0.0;
  if(RintervalCnt > 0){
    for(int j = 1; j < RintervalCnt; j++){
      Cinterval *pi = &Rinterval[i];
      Cinterval *pj = &Rinterval[j];
      Rtotlen += pj->right - pj->left;
      
      if(pi->id == pj->id && pj->left <= pi->right){/* merge intervals */
	if(DEBUG) assert(pi->left <= pj->left);
	pi->right = max(pi->right, pj->right);
	continue;
      }
      
      Runiqlen += pi->right - pi->left;
      
      i++;
      if(i < j)
	Rinterval[i] = Rinterval[j];
    }
    Runiqlen += Rinterval[i].right - Rinterval[i].left;
    RintervalCnt = i+1;
  }
  
  FLOAT RuniqRange = 0.0, RlinkLen = 0.0;
  for(int i = 0; i < RintervalCnt; i++){
    Cinterval *pi = &Rinterval[i];
    long long id = pi->id;

    if(i > 0 && Rinterval[i-1].id == id)
      continue;
    
    Cmap *pmap = YYmap[pi->mapid];
    int N = pmap->numsite[0];
    RlinkLen += pmap->site[0][N+1];

    FLOAT minleft = pi->left, maxright = pi->right;
    for(int j = i+1; j < RintervalCnt; j++){
      Cinterval *pj = &Rinterval[j];
      if(pj->id != id)
	break;
      minleft = min(minleft,pj->left);
      maxright = max(maxright,pj->right);
    }
    
    RuniqRange += maxright - minleft;
  }

  i = 0;
  FLOAT Qtotlen = (QintervalCnt <= 0) ? 0.0 : Qinterval[0].right - Qinterval[0].left;
  FLOAT Quniqlen = 0.0;
  if(QintervalCnt > 0){
    for(int j = 1; j < QintervalCnt; j++){
      Cinterval *pi = &Qinterval[i];
      Cinterval *pj = &Qinterval[j];
      Qtotlen +=  pj->right - pj->left;

      if(pi->id == pj->id && pj->left <= pi->right){/* merge intervals */
	if(DEBUG) assert(pi->left <= pj->left);
	pi->right = max(pi->right, pj->right);
	continue;
      }
    
      Quniqlen += pi->right - pi->left;

      if(++i < j)
	Qinterval[i] = Qinterval[j];
    }
    QintervalCnt = i+1;
  }

  if(VERB>=2){
    printf("After merging overlapped intervals : RintervalCnt=%d,QintervalCnt=%d\n",RintervalCnt,QintervalCnt);
    fflush(stdout);
  }

  FLOAT QuniqRange = 0.0, QlinkLen = 0.0;
  for(int i = 0; i < QintervalCnt; i++){
    Cinterval *pi = &Qinterval[i];
    long long id = pi->id;

    if(i > 0 && Qinterval[i-1].id == id)
      continue;
    
    if(DEBUG && !(0 <= pi->mapid && pi->mapid < numX && XXmap && XXmap[pi->mapid])){
      printf("i=%d,QintervalCnt=%d:pi->id=%lld,pi->mapid=%d,numX=%d\n",i,QintervalCnt,pi->id,pi->mapid,numX);
      fflush(stdout);
      assert(0 <= pi->mapid && pi->mapid < numX);
      assert(XXmap);
      assert(XXmap[pi->mapid]);
    }
    Cmap *pmap = XXmap[pi->mapid];
    int N = pmap->numsite[0];
    QlinkLen += pmap->site[0][N+1];

    FLOAT minleft = pi->left, maxright = pi->right;
    for(int j = i+1; j < QintervalCnt; j++){
      Cinterval *pj = &Qinterval[j];
      if(pj->id != id)
	break;
      minleft = min(minleft,pj->left);
      maxright = max(maxright,pj->right);
    }
    
    QuniqRange += maxright - minleft;
  }

  if(VERB){
    printf("Reference : Total Aligned length = %0.3f Mb, Unique Aligned Length = %0.3f Mb, Aligned Ranges = %0.3f Mb, Aligned Maps Len = %0.3f Mb\n", 
	   Rtotlen * 1e-6, Runiqlen * 1e-6, RuniqRange * 1e-6, RlinkLen * 1e-3);
    printf("Query : Total Aligned length = %0.3f Mb, Unique Aligned Length = %0.3f Mb, Aligned Ranges = %0.3f Mb, Aligned Maps Len = %0.3f Mb\n", 
	   Qtotlen * 1e-6, Quniqlen * 1e-6, QuniqRange * 1e-6, QlinkLen * 1e-3);
    fflush(stdout);
  }

  delete [] Rinterval;

  FILEclose(fp);
  free(write_buffer);

  if(fpSV){
    FILEclose(fpSV);
    free(write_bufferSV);
  }

  /* strip .tmp from output file name */
  errno = 0;
  if(rename(tmpfilename,origfilename)){
    int eno = errno;
    char *err = strerror(eno);
    printf("failed to rename %s to %s:errno=%d:%s\n",tmpfilename,origfilename,eno,err);
    fflush(stdout);  exit(1);
  }
  free(tmpfilename);
  free(origfilename);
}

/* output alignment[align_start .. align_end-1] that are above threshold.
   output .xmap in <basename>.xmap
   Reference maps already exist in <basename>_r.cmap
   Query maps already exist in queryfilename
*/
void output_xmap(char *basename, char * &queryfilename, size_t align_start, size_t align_end, bool full)
{
  if(DEBUG) assert(colors==1);

  if(VERB>=2){
    printf("output_xmap:basename=%s,align_start=%lu,align_end=%lu,full=%d\n",basename,align_start,align_end,full);
    fflush(stdout);
  }

  if(strstr(basename,"/dev/null"))
    return;

  char filename[PATH_MAX];
  strcpy(filename,basename);
  /* no need to strip suffix : for default basename this was already done in RefAligner.cpp */
  if(1){
    int i = strlen(filename);
    strcpy(&filename[i],".xmap");
  }
  if(checkFile(filename))
    return;

  /* get current directory pathname */
  char cwd[PATH_MAX];
  if(getcwd(cwd,PATH_MAX) == NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("getcwd() failed (errno=%d):%s\n",eno, err);
    exit(1);
  }

#ifndef WIN32
  /* convert basename into absolute pathname */
  char absbasename[PATH_MAX];
  if(basename[0]=='/')
    strcpy(absbasename,basename);
  else
    sprintf(absbasename,"%s/%s",cwd,basename);
  if(queryfilename && queryfilename[0] != '/'){/* convert queryfilename into an absolute pathname */
    char absname[PATH_MAX];
    sprintf(absname,"%s/%s",cwd,queryfilename);
    free(queryfilename);
    queryfilename = strdup(absname);
  }
#else
  char *absbasename = basename;
#endif

  if(VERB>=2){
    printf("output_xmap(%s): align_start=%lu,align_end=%lu\n",filename,align_start,align_end);
    fflush(stdout);
  }

  /* remove below threshold or deallocated entries */
  if(RepeatShift > 0.0 || FirstAlignments > 0){
    if(DEBUG) assert(align_start == 0 && align_end == numaligns);
    size_t aligncnt = 0;
    for(register size_t i=0; i < numaligns; i++){
      register Calign *p = alignment[i];
      if(!p)
	continue;
      alignment[i] = 0;
      if(RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || p->numpairs < AlignedSiteThreshold2) 
	 : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold)){
	delete [] p;
	continue;
      }
	
      if(DEBUG) assert(alignment[i] == 0);
      alignment[aligncnt++] = p;
    }
    align_end = numaligns = aligncnt;
  }

  Calign **Lalignment = 0;
  size_t start = align_start;
  size_t end = align_end;

  if(SecondBest && !SinglelineXmap){/* interleave 2nd best alignments after best alignments */
    /* NOTE : the 2nd best alignments are marked with chimpair = -1, others with chimpair = 0 */
    delete [] Lalignment;
    Lalignment = new Calign*[(align_end - align_start)*2];
    start = 0;
    end = 0;
    size_t SecondCnt = 0;
    size_t SplitCnt = 0;
    if(VERB>=2){
      printf("ScoreThreshold=%0.6f,LogPvThreshold=%0.2f,AlignedSiteThreshold=%d\n",ScoreThreshold,LogPvThreshold,AlignedSiteThreshold);
      printf("SecondBest alignments (ScoreThreshold2=%0.6f,LogPvThreshold2=%0.2f,AlignedSiteThreshold2=%d):\n", ScoreThreshold2,LogPvThreshold2,AlignedSiteThreshold2);
      fflush(stdout);
    }

    if(DEBUG) assert(!RefSplit);

    for(size_t i = align_start; i < align_end; i++){
      Calign *p = alignment[i];
      if(!p)
	continue;
      int U = p->numpairs;
      p->chimpair = 0;
      if(U<=0 || (!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold){
	if(VERB>=2){
	  printf("alignment[%lu/%lu]= %p:refid=%d,mapid=%d,or=%d:numpairs=%d,score=%0.3f,logPV=%0.2f, align2= %p: skipping\n",i, align_end, p, 
		 p->mapid1,p->mapid2,p->orientation,U, p->score, p->logPV, p->align2);
	  fflush(stdout);
	}
	continue;
      }

      if(Gmap[p->mapid2]->origmap){
	SplitCnt++;
	if(VERB>=2){
	  printf("Lalignment[%lu]= %p: numpairs=%d,score=%0.3f,logPV=%0.2f: refid=%d, split mapid=%d(orig=%d,id=%lld),or=%d\n",
		 end, p, U, p->score, p->logPV, p->mapid1,p->mapid2,Gmap[p->mapid2]->origmap->mapid,Gmap[p->mapid2]->id,p->orientation);
	  fflush(stdout);
	}
      }
      Lalignment[end++] = p;

      if(p->align2){
	Calign *q = p->align2;
	if(DEBUG) assert(q->numpairs > 0);
	if((!REFSPLIT_FIX && q->score <= ScoreThreshold2) || q->logPV <= LogPvThreshold2){
	  if(VERB>=2){
	    printf("Lalignment[%lu]= %p:refid=%d,mapid=%d,or=%d:numpairs=%d,score=%0.3f,logPV=%0.2f, p->align2= %p: numpairs=%d,score=%0.3f,logPV=%0.2f (skipping align2)\n",
		   end-1, p, p->mapid1, p->mapid2, p->orientation, U, p->score, p->logPV, p->align2, q->numpairs, q->score, q->logPV);
	    fflush(stdout);
	  }
	  delete [] p->align2; p->align2 = NULL;
	  continue;
	}
	if(DEBUG) assert(q->align2 == NULL);
	if(VERB>=2){
	  printf("Lalignment[%lu]= %p: numpairs=%d,score=%0.3f,logPV=%0.2f: Lalignment[%lu]->p->align2= %p: refid=%d,mapid=%d(%lld),or=%d,numpairs=%d,score=%0.3f,logPV=%0.2f\n",
		 end-1, p, U, p->score, p->logPV, end, q, q->mapid1,q->mapid2,Gmap[q->mapid2]->id, q->orientation, q->numpairs, q->score, q->logPV);
	  fflush(stdout);
	}
	q->chimpair = -1;
	Lalignment[end++] = q;
	SecondCnt++;
      }
    }
    if(VERB){
      printf("%lu Xmap alignments (including %lu splits) + %lu SecondBest alignments\n", end - SecondCnt, SplitCnt, SecondCnt);
      fflush(stdout);
    }
  }

  if(MultiMatches){/* output all matchgroups and tag the best one for each map pair */
    /* NOTE : the 2nd best alignments are marked with chimpair = -1, best one with chimpair = 0 */

    // first count how many total matchgroups there are
    size_t cnt = 0;
    for(size_t i = align_start; i < align_end; i++){
      Calign *p = alignment[i];
      if(!p)
	continue;
      int U = p->numpairs;
      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || U < AlignedSiteThreshold2)
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || U < AlignedSiteThreshold))){
	if(VERB>=2 && YYmap[p->mapid1]->id == 25 && XXmap[p->mapid2]->id == 6170){
	  printf("alignment[%lu]:pairs=%d,score=%0.6f,logPV=%0.6f,logPV2=%0.6f,multicnt=%d:refid=%d(id=%lld),mapid=%d(id=%lld),or=%d: skipping\n",
		 i,p->numpairs,p->score,p->logPV,RefSplit ? p->logPV2 : p->logPV, p->multicnt, p->mapid1,YYmap[p->mapid1]->id,p->mapid2,XXmap[p->mapid2]->id,p->orientation);
	  fflush(stdout);
	}
	continue;
      }
      if((DEBUG>=1+RELEASE && RefSplit && !(p->logPV2 >= p->logPV)) ||VERB>=2){
	int rid = p->mapid1;
	int mid = p->mapid2;

	#pragma omp critical
	{
	  printf("alignment[i=%lu]:p->score=%0.6f,numpairs=%d,logPV=%0.6f,logPV2=%0.6f,multicnt=%d(cnt=%lu):mid=%lld,rid=%lld,or=%d\n",
		 i,p->score,p->numpairs,p->logPV,RefSplit ? p->logPV2 : p->logPV, p->multicnt,cnt,XXmap[mid]->id,YYmap[rid]->id,p->orientation);
	  fflush(stdout);

	  if(DEBUG>=1+RELEASE && RefSplit) assert(p->logPV2 >= p->logPV);
	}
      }
      if(RefSplit ? p->logPV2 <= LogPvThreshold2 : p->logPV <= LogPvThreshold){
	if(VERB>=2 && YYmap[p->mapid1]->id == 25 && XXmap[p->mapid2]->id == 6170){
	  printf("alignment[%lu]:pairs=%d,score=%0.6f,logPV=%0.6f,logPV2=%0.6f,multicnt=%d:refid=%d(id=%lld),mapid=%d(id=%lld),or=%d: skipping due to LogPvThreshold= %0.2f\n",
		 i,p->numpairs,p->score,p->logPV,RefSplit ? p->logPV2 : p->logPV, p->multicnt, p->mapid1,YYmap[p->mapid1]->id,p->mapid2,XXmap[p->mapid2]->id,p->orientation,LogPvThreshold);
	  fflush(stdout);
	}
	continue;
      }
      if(DEBUG>=2/* HERE HERE >= 1+RELEASE */ && !extendSplit && !RefSplit && !(p->Malign && p->multicnt > 0)){
	int rid = p->mapid1;
	int mid = p->mapid2;
	#pragma omp critical
	{
	  printf("alignment[i=%lu]:p->score=%0.6f,numpairs=%d,logPV=%0.6f,logPV2=%0.6f,multicnt=%d(cnt=%lu):mid=%lld,rid=%lld,or=%d:p->Malign= %p\n",
		 i,p->score,p->numpairs,p->logPV,RefSplit ? p->logPV2 : p->logPV, p->multicnt,cnt,XXmap[mid]->id,YYmap[rid]->id,p->orientation,p->Malign);
	  fflush(stdout);
	
	  if(DEBUG>=1+RELEASE) assert(p->Malign && p->multicnt > 0);
	}
      }
      cnt += p->Malign ? p->multicnt : 1;
    }
    Lalignment = new Calign*[cnt];
    start = 0;
    end = 0;
    size_t SecondCnt = 0;

    for(size_t i = align_start; i < align_end; i++){
      Calign *p = alignment[i];
      if(!p)
	continue;
      int U = p->numpairs;
      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold)))
	continue;

      if(VERB>=2 && YYmap[p->mapid1]->id == 25285LL){
	int rid = p->mapid1;
	int mid = p->mapid2;

	#pragma omp critical
	{
	  printf("alignment[i=%lu]:p->score=%0.6f,numpairs=%d,logPV=%0.6f,logPV2=%0.6f,multicnt=%d(SecondCnt=%lu):mid=%lld,rid=%lld,or=%d,Malign=%p,end=%lu\n",
		 i,p->score,p->numpairs,p->logPV,RefSplit ? p->logPV2 : p->logPV,p->multicnt,SecondCnt,XXmap[mid]->id,YYmap[rid]->id,p->orientation,p->Malign,end);
	  if(p->Malign){
	    for(int t = 0; t < p->multicnt; t++){
	      Calign *q = p->Malign[t];
	      printf("\t Malign[%d]:rid=%lld,mid=%lld,or=%d:pairs=%d,score=%0.6f,logPV=%0.2f,logPV2=%0.2f\n", t, YYmap[q->mapid1]->id,XXmap[q->mapid2]->id,q->orientation,q->numpairs,q->score,q->logPV,q->logPV2);
	    }
	  }
	  fflush(stdout);
	}
      }

      if(!p->Malign){
	if((DEBUG>=2 /* HERE HERE >=1+RELEASE*/ && !(extendSplit || RefSplit)) ||
	   (VERB>=2 && YYmap[p->mapid1]->id == 25 && XXmap[p->mapid2]->id == 6170)){
	  printf("alignment[i=%lu]:refid=%d(id=%lld),mapid=%d(id=%lld),p->score=%0.6f,pairs=%d,logPV=%0.6f,logPV2=%0.6f,multicnt=%d(SecondCnt=%lu),or=%d,Malign=%p,end=%lu\n",
		 i,p->mapid1,YYmap[p->mapid1]->id,p->mapid2,XXmap[p->mapid2]->id,p->score,p->numpairs,p->logPV,RefSplit ? p->logPV2 : p->logPV,p->multicnt,SecondCnt,p->orientation,p->Malign,end);
	  fflush(stdout);
	  assert(extendSplit || RefSplit);
	}
	p->chimpair = 0;
	Lalignment[end++] = p;
      } else {
	for(int t = 0; t < p->multicnt; t++){
	  Calign *q = p->Malign[t];
	  if(DEBUG) assert(0 <= q->mapid1 && q->mapid1 < numY);
	  if(DEBUG) assert(0 <= q->mapid2 && q->mapid2 < numX);
	  if((DEBUG>=1+RELEASE && RefSplit && !(q->logPV2 >= q->logPV)) ||
	   (VERB>=2 && YYmap[p->mapid1]->id == 25 && XXmap[p->mapid2]->id == 6170)){
	    printf("alignment[i=%lu]:refid=%d(id=%lld),mapid%d(id=%lld),or=%d:pairs=%d,score=%0.6f,logPV=%0.6f,logPV2=%0.6f,multicnt=%d:t=%d:pairs=%d,score=%0.6f,logPV=%0.6f,logPV2=%0.6f,or=%d\n",
		   i,p->mapid1,YYmap[p->mapid1]->id,p->mapid2,XXmap[p->mapid2]->id,p->orientation,p->numpairs,p->score,p->logPV,p->logPV2,p->multicnt,t,q->numpairs,q->score,q->logPV,q->logPV2,q->orientation);
	    fflush(stdout);
	    assert(q->logPV2 >= q->logPV);
	  }

	  q->chimpair = (t==0) ? 0 : -1;
	  Lalignment[end++] = q;
	}
	SecondCnt += p->multicnt - 1;
      }
    }
    //    fprintf(fp,"# BestID refers to XmapEntryID of Best alignment for same query and reference\n");
    if(VERB){
      printf("%lu Xmap best alignments + %lu alternate alignments\n", end - SecondCnt, SecondCnt);
      fflush(stdout);
    }
  }

  if(!Lalignment){/* make a local copy so we can re-order the alignments without affecting the calling function */
    Lalignment = new Calign*[end - start];
    for(size_t i = start; i < end; i++)
      Lalignment[i-start] = alignment[i];
    end -= start;
    start = 0;
  }

  if(!SecondBest){  /* filter out alignments that are NULL, below threshold or with mapWT < UnrefineMinMapwt */
    size_t aligncnt = start;
    for(size_t i = start; i < end; i++){
      Calign *p = Lalignment[i];
      if(!p)
	continue;
      int U = p->numpairs;
      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold)))
	continue;
      if(p->mapWT < UnrefineMinMapwt)
	continue;
      Lalignment[aligncnt++] = p;
    }
    if(VERB>=2 && aligncnt < end && aligncnt > start){
      int rid = Lalignment[start]->mapid1;
      printf("RefId=%lld:Reduced number of alignments from %lu to %lu due to -MinMapwt %0.2f %0.2f etc\n",YYmap[rid]->id,end-start,aligncnt-start, RefineMinMapwt,UnrefineMinMapwt);
      fflush(stdout);
    }
    end = aligncnt;
  }
  

  double Xscale = 1.0;
  double Yscale = 1.0;
  /* NOTE if -ref or -mres -bpp were used, rescaled _r.cmap and _q.cmap files will be used, so no need to rescale the alignment locations */
  if(!(spots_filename[0] || (svcheck && numrefmaps > 0)) && !(mres >= 0.001 || mresSD > 0.0 || fabs(PixelLen - origPixelLen) > 1e-12)){
    Yscale = Xscale = origPixelLen/PixelLen;
    if(PairSplit && PairsplitRef)/* -i inputs that correspond to reference is not scaled by bpp/500 */
      Yscale = 1.0;
  }

  double *RefFilteredConf = new double[numrefmaps];

  if(XmapNoQueryOverlap){
    if(RefSplit){
      printf("-XmapNoQueryOverlap not implemented with -RefSplit\n");
      fflush(stdout);exit(1);
    }

    /* sort alignments in order of id2, logPV  (null alignments last) */
    qsort(&Lalignment[start],end - start,sizeof(Calign *), (intcmp*)Id2SiteInc);

    for(size_t i = start; i < end; i++){
      Calign *p = Lalignment[i];
      if(!p)
	continue;
      RefFilteredConf[p->mapid1] = 0.0;
    }

    for(size_t i = start; i+1  < end; i++){
      Calign *p = Lalignment[i];
      if(!p)
	continue;
      int U = p->numpairs;
      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold)))
	continue;
      if(SecondBest && p->chimpair)
	continue;/* skip 2nd best alignments */

      int rid = p->mapid1;
      int mid = p->mapid2;
      Cmap *Ymap = YYmap[rid];
      Cmap *Xmap = XXmap[mid];
      if(DEBUG) assert(Xmap->origmap == NULL);
      int N = Ymap->numsite[0];
      int M = Xmap->numsite[0];
      FLOAT *Y = Ymap->site[0];
      FLOAT *X = Xmap->site[0];
      int LI = p->sites1[0];
      int LJ = p->sites2[0];
      int RI = p->sites1[U-1];
      int RJ = p->sites2[U-1];
      if(DEBUG) assert(LJ <= RJ);
      int qrystart = p->orientation ? M + 1 - RJ : LJ;
      int qryend   = p->orientation ? M + 1 - LJ : RJ;
      if(DEBUG) assert(qrystart <= qryend);
      double qrystartpos = X[qrystart]*Xscale;/* QryStartPos */
      double qryendpos = X[qryend]*Xscale;/* QryStartPos */
      double refstartpos = Y[LI] * Yscale; /* RefStartPos */
      double refendpos   = Y[RI] * Yscale; /* RefEndPos */

      for(size_t k = i+1; k < end; k++){
	Calign *q = Lalignment[k];
	if(!q)
	  continue;
	if(mid != q->mapid2)
	  break;

	int U2 = q->numpairs;
	if(U2 <= 0 || (RefSplit ? ((!REFSPLIT_FIX && q->score <= ScoreThreshold2) || q->logPV2 <= LogPvThreshold2 || U2 < AlignedSiteThreshold2) 
		       : ((!REFSPLIT_FIX && q->score <= ScoreThreshold) || q->logPV <= LogPvThreshold || U2 < AlignedSiteThreshold)))
	  continue;
	if(SecondBest && q->chimpair < 0)
	  continue;/* skip 2nd best alignments */

	int rid2 = q->mapid1;
	Cmap *Ymap2 = YYmap[rid2];
	int N2 = Ymap2->numsite[0];
	FLOAT *Y2 = Ymap2->site[0];
	if(DEBUG) assert(Ymap2->origmap == NULL);
	int LI2 = q->sites1[0];
	int LJ2 = q->sites2[0];
	int RJ2 = q->sites2[q->numpairs-1];
	int RI2 = q->sites1[q->numpairs-1];
	if(DEBUG) assert(LJ2 <= RJ2);
	int qrystart2 = q->orientation ? M + 1 - RJ2 : LJ2;
	int qryend2   = q->orientation ? M + 1 - LJ2 : RJ2;
	if(DEBUG) assert(qrystart2 <= qryend2);	
	double qrystartpos2 = X[qrystart2]*Xscale;/* QryStartPos2 */
	double qryendpos2 = X[qryend2]*Xscale;/* QryStartPos2 */
	double refstartpos2 = Y2[LI2]*Yscale; /* RefStartPos2 */
	double refendpos2   = Y2[RI2]*Yscale; /* RefEndPos2 */
	
	/* now check if alignment q should be filtered out due to XmapNoQueryOverlap : If so update RefFilteredConf[rid2] */

	/* check if p & q overlap by at least XmapNoQueryOverlap kb in query map */
	if(max(qrystartpos,qrystartpos2) + XmapNoQueryOverlap < min(qryendpos,qryendpos2)){
	  if(VERB>=2){
	    printf("Overlapped QueryId=%lld(M=%d,len=%0.3f) with RefId=%lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f and %0.2f):Igoring 2nd alignment\n",
		   Xmap->id,M,X[M+1]*Xscale,YYmap[rid]->id,N,Y[N+1]*Yscale,LI,RI,refstartpos,refendpos,qrystart,qryend,p->orientation,qrystartpos,qryendpos,
		   YYmap[q->mapid1]->id,N2,Y2[N2+1]*Yscale,LI2,RI2,refstartpos2,refendpos2,qrystart2,qryend2,q->orientation,qrystartpos2,qryendpos2,p->logPV,q->logPV);
	    fflush(stdout);
	  }
	  RefFilteredConf[rid2] = max(RefFilteredConf[rid2], q->logPV);

	  // NOTE : cannot actually delete the alignment since -resEstimate may need it to estimate res & resSD during last iteration in refalign() : just set the logPV below the threshold
	  Lalignment[k]->logPV = LogPvThreshold2 - 1.0;
	  continue;
	}
      }/* for k = i+1 .. end-1 */
    } /* for i = 0 .. end-2 */
  }

  if(XmapFilterThresh1 > 0.0){
    if(RefSplit){
      printf("-XmapFilterThresh not implemented with -RefSplit\n");
      fflush(stdout);exit(1);
    }
    if(DEBUG) assert(XmapFilterThresh1 >= XmapFilterThresh2 && XmapFilterThresh2 > XmapFilterDelta);

    /* sort alignments in order of id1, logPV (null alignments last) */
    qsort(&Lalignment[start],end - start,sizeof(Calign *), (intcmp*)Id1SiteInc);

    int inext = -1;
    for(size_t i = start; i+1  < end; i = inext){
      inext = i+1;

      Calign *p = Lalignment[i];
      if(!p)
	continue;
      int U = p->numpairs;
      if(U<=0 || p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold)
	continue;
      if(SecondBest && p->chimpair)
	continue;/* skip 2nd best alignments */

      int rid = p->mapid1;
      int mid = p->mapid2;

      for(size_t k = i+1; k < end; k++){
	Calign *q = Lalignment[k];
	if(!q)
	  continue;
	if(rid != q->mapid1){
	  inext = k;
	  break;
	}

	if(q->numpairs <= 0 || q->score <= ScoreThreshold || q->logPV <= LogPvThreshold || q->numpairs < AlignedSiteThreshold)
	  continue;
	if(SecondBest && q->chimpair < 0)
	  continue;/* skip 2nd best alignments */

	int mid2 = q->mapid2;
	if(!((p->logPV > q->logPV + XmapFilterDelta && p->logPV > XmapFilterThresh2 && RefFilteredConf[rid] <= p->logPV) || p->logPV > XmapFilterThresh1)){/* remove all matchgroups with mapid1 == rid */
	  // NOTE : cannot actually delete the alignment since -resEstimate may need it to estimate res & resSD during last iteration in refalign() : just set the logPV below the threshold
	  if(VERB/* HERE >=2 */){
	    if(DEBUG) assert(p->logPV <= XmapFilterThresh1);

	    Cmap *Ymap = YYmap[rid];
	    int N = Ymap->numsite[0];
	    FLOAT *Y = Ymap->site[0];

	    Cmap *Xmap = XXmap[mid];
	    int M = Xmap->numsite[0];
	    FLOAT *X = Xmap->site[0];

	    Cmap *Xmap2 = XXmap[mid2];
	    int M2 = Xmap2->numsite[0];
	    FLOAT *X2 = Xmap2->site[0];

	    if(p->logPV <= q->logPV + XmapFilterDelta){

	      printf("Deleting ambiguous alignment of RefId=%lld(N=%d,len=%0.3f) with QueryId=%lld(M=%d,len=%0.3f,or=%d,logPV=%0.2f) and %lld(M=%d,len=%0.3f:or=%d,log10PV=%0.2f),Delta=%0.2f,ThreshPV1=%0.2f\n",
		     Ymap->id,N,Y[N+1],Xmap->id,M,X[M+1],p->orientation,p->logPV,Xmap2->id,M2,X2[M2+1],q->orientation,q->logPV,XmapFilterDelta,XmapFilterThresh1);
	    } else if(p->logPV <= RefFilteredConf[rid]){
	      printf("Deleting non-best alignment of RefId=%lld(N=%d,len=%0.3f) with QueryId=%lld(M=%d,len=%0.3f,or=%d,logPV=%0.2f) and %lld(M=%d,len=%0.3f:or=%d,log10PV=%0.2f),BestPV=%0.2f,ThreshPV1=%0.2f\n",
		     Ymap->id,N,Y[N+1],Xmap->id,M,X[M+1],p->orientation,p->logPV,Xmap2->id,M2,X2[M2+1],q->orientation,q->logPV,RefFilteredConf[rid],XmapFilterThresh1);
	    } else // HERE : this should never happen
	      printf("Deleting unambiguously best alignment of RefId=%lld(N=%d,len=%0.3f) with QueryId=%lld(M=%d,len=%0.3f,or=%d,logPV=%0.2f) and %lld(M=%d,len=%0.3f:or=%d,log10PV=%0.2f) due to ThreshPV2=%0.2f\n",
		     Ymap->id,N,Y[N+1],Xmap->id,M,X[M+1],p->orientation,p->logPV,Xmap2->id,M2,X2[M2+1],q->orientation,q->logPV,XmapFilterThresh2);
	    fflush(stdout);
	  }
	  Lalignment[i]->logPV = LogPvThreshold2 - 1.0;
	}

	/* NOTE : We only care about the best alignment with same refid, all others will be deleted */
	Lalignment[k]->logPV = LogPvThreshold2 - 1.0;

	/* scan ahead to next mapid1 (deleting all alignments with same refid) */
	for(k++; k < end; k++){
	  Calign *q = Lalignment[k];
	  if(!q)
	    continue;
	  if(rid != q->mapid1){
	    inext = k;
	    break;
	  }
	  Lalignment[k]->logPV = LogPvThreshold2 - 1.0;
	}
	break;
      }
    }    
  }

  delete [] RefFilteredConf;

  FILE *fpCH = NULL;
  char filenameCH[PATH_MAX];
  //  char *write_bufferCH = NULL;

  if(XmapChim || XmapUnique)  /* sort alignments in ascending order of id2, descending order of logPV (WAS77 logPV2 if RefSplit)  (null alignments last) : break ties in ascending order of id1 (prefers chr X over Y) */
    qsort(&Lalignment[start],end - start,sizeof(Calign *), (intcmp*)Id2SiteInc);

  /* keep track of query map regions with overlapped alignments : suppress the overlapped part of the alignment with the worse logPV when reporting indels (HERE HERE : should compare logPV in overlap region only) */
  if(INDEL || XmapUnique)
    for(size_t i = start; i < end; i++){
      Calign *p = Lalignment[i];
      if(!p) {
	if(DEBUG>=1+RELEASE) assert(!(RefSplit && MultiMatches));
	continue;
      }
      int U = p->numpairs;

      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold))){
	if(DEBUG>=1+RELEASE && (RefSplit && MultiMatches)){
	  printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:pairs=%d,score=%0.6f,logPV=%0.2f(logPV2=%0.2f),S2=%0.6f,T2=%0.2f,A2=%d\n",
		 p->mapid1,YYmap[p->mapid1]->id,p->mapid2,XXmap[p->mapid2]->id,p->orientation,U,p->score,p->logPV,p->logPV2,ScoreThreshold2,LogPvThreshold2,AlignedSiteThreshold2);
	  fflush(stdout);
	  assert(!(RefSplit && MultiMatches));
	}
	continue;
      }
      int mid = p->mapid2;
      Cmap *Xmap = XXmap[mid];
      int M = Xmap->numsite[0];

      if(p->suppress)
	delete [] p->suppress;
      p->suppress = new int[M+2];
      for(int I = 0; I <= M+1; I++)
	p->suppress[I] = 0;
      /* NOTE: p->suppress[I] is a union of various ranges and will be computed as the change between two consecutive values of the number of ranges overlapping any label
	 and then the actual values accumulated at the end : All non-zero values will represent regions for which indel output will be suppressed for alignment p */
    }

  if(XmapUnique){  /* delete alignments if there is another alignment of the same query map with better logPV and the first alignment does not have at least XmapUnique unique sites on the query map */
    int delcnt = 0;

    for(size_t i = start; i+1 < end; i++){
      Calign *p = Lalignment[i];
      if(!p){
	if(DEBUG>=1+RELEASE) assert(!(RefSplit && MultiMatches));
	continue;
      }
      int U = p->numpairs;
      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold))){
	if(DEBUG) assert(!(RefSplit && MultiMatches));
	continue;
      }

      if(SecondBest && p->chimpair){
	if(DEBUG>=1+RELEASE) assert(!(RefSplit && MultiMatches));
	continue;/* skip 2nd best alignments */
      }
      
      int rid = p->mapid1;
      int mid = p->mapid2;

      Cmap *Ymap = YYmap[rid];
      Cmap *Xmap = XXmap[mid];
      if(VERB>=2 && Ymap->id == 9 && Xmap->id == 2183){
	printf("Checking Lalignment[%lu]:mid=%lld,rid=%lld:numpairs=%d,score=%0.6f,logPV=%0.2f(%0.2f):XmapUnique=%d\n",i,Xmap->id,Ymap->id,U,p->score,p->logPV,RefSplit ? p->logPV2 : p->logPV,XmapUnique);
	fflush(stdout);
      }
      Cmap *origmap = Xmap;
      while(origmap->origmap){
	Xmap->origmap = origmap = origmap->origmap;
	if(origmap->origmap){/* otherwise origmap->left[] is not defined */
	  Xmap->left[0] += max(0,origmap->left[0]-1);
	  Xmap->right[0] += max(0,origmap->left[0]-1);
	}
      }
      int N = Ymap->numsite[0];
      int M = Xmap->numsite[0];
      FLOAT *Y = Ymap->site[0];
      FLOAT *X = Xmap->site[0];
      FLOAT *origX = Xmap->origmap ? Xmap->origmap->site[0] : X;
      FLOAT *origY = Ymap->origmap ? Ymap->origmap->site[0] : Y;
      if(DEBUG>=2) assert(Xmap->origmap==0 || Xmap->origmap->origmap==0);/* If this fails fix origmap link in refalign() to point to root map */
      if(DEBUG>=2) assert(Ymap->origmap==0 || Ymap->origmap->origmap==0);/* If this fails fix origmap link (same as for Xmap above) to point to root map */
      int Yshift = Ymap->origmap ? max(1,Ymap->left[0]) - 1 : 0;
      int Xshift = Xmap->origmap ? max(1,Xmap->left[0]) - 1 : 0;
      int LI = p->sites1[0];
      int LJ = p->sites2[0];
      int RI = p->sites1[U-1];
      int RJ = p->sites2[U-1];
      if(DEBUG) assert(LJ <= RJ);
      int qrystart = p->orientation ? Xshift + M + 1 - RJ : Xshift + LJ;
      int qryend   = p->orientation ? Xshift + M + 1 - LJ : Xshift + RJ;
      double refstartpos = origY[Yshift + LI]*Yscale; /* RefStartPos */
      double refendpos   = origY[Yshift + RI]*Yscale; /* RefEndPos */
      if(DEBUG) assert(qrystart <= qryend);

      for(size_t k = i+1; k < end; k++){
	Calign *q = Lalignment[k];
	if(!q){
	  if(DEBUG>=1+RELEASE) assert(!(RefSplit && MultiMatches));
	  continue;
	}
	if(mid != q->mapid2)
	  break;
	if(VERB>=2 && Ymap->id == 9 && Xmap->id == 2183){
	  printf("Checking alignments i=%lu,k=%lu:mid1=%lld,mid2=%lld,rid1=%lld,rid2=%lld:q->numpairs=%d,score=%0.6f,logPV=%0.2f\n",
		 i,k,Xmap->id,XXmap[q->mapid2]->id,Ymap->id,YYmap[q->mapid1]->id,q->numpairs,q->score,q->logPV);
	  fflush(stdout);
	}
	int U2 = q->numpairs;
	if(U2 <= 0 || (RefSplit ? ((!REFSPLIT_FIX && q->score <= ScoreThreshold2) || q->logPV2 <= LogPvThreshold2 || U2 < AlignedSiteThreshold2) 
		       : ((!REFSPLIT_FIX && q->score <= ScoreThreshold) || q->logPV <= LogPvThreshold || U2 < AlignedSiteThreshold))){
	  if(DEBUG>=1+RELEASE) assert(!(RefSplit && MultiMatches));
	  continue;
	}
	if(SecondBest && q->chimpair < 0){
	  if(DEBUG>=1+RELEASE) assert(!(RefSplit && MultiMatches));
	  continue;/* skip 2nd best alignments */
	}
	Cmap *Ymap2 = YYmap[q->mapid1];
	int N2 = Ymap2->numsite[0];
	FLOAT *Y2 = Ymap2->site[0];
	FLOAT *origY2 = Ymap2->origmap ? Ymap2->origmap->site[0] : Y2;
	if(DEBUG>=2) assert(Ymap2->origmap==0 || Ymap2->origmap->origmap==0);/* If this fails fix origmap link (same as for Xmap above) to point to root map */
	int Yshift2 = Ymap2->origmap ? max(1,Ymap2->left[0]) - 1 : 0;

	/* check if alignments p & q have at least XmapUnique distinct sites on both p->mapid2 and q->mapid2 */
	int LI2 = q->sites1[0];
	int LJ2 = q->sites2[0];
	int RJ2 = q->sites2[U2 - 1];
	int RI2 = q->sites1[U2-1];
	if(DEBUG && !(LJ2 <= RJ2)){
	  printf("i=%lu,k=%lu:rid=%d,mid=%d,q->mapid1=%d,q->mapid2=%d:LJ=%d,RJ=%d,LJ2=%d,RJ2=%d\n",
		 i,k,rid,mid,q->mapid1,q->mapid2,LJ,RJ,LJ2,RJ2);
	  fflush(stdout);
	  assert(LJ2 <= RJ2);
	}
	int qrystart2 = q->orientation ? Xshift + M + 1 - RJ2 : Xshift + LJ2;
	int qryend2   = q->orientation ? Xshift + M + 1 - LJ2 : Xshift + RJ2;
	double refstartpos2 = origY2[Yshift2 + LI2]*Yscale; /* RefStartPos */
	double refendpos2   = origY2[Yshift2 + RI2]*Yscale; /* RefEndPos */
	if(DEBUG) assert(qrystart2 <= qryend2);
	
	if(VERB>=3 && Ymap->id == 9 && Xmap->id == 2183){
	  printf("   LI=%d,RI=%d,LJ=%d,RJ=%d,M=%d,or=%d;LI2=%d,RI2=%d,LJ2=%d,RJ2=%d,M2=%d,or2=%d;qrystart=%d(%0.3f),qryend=%d(%0.3f),qrystart2=%d(%0.3f),qryend=%d(%0.3f),XmapChim=%d\n",
		 LI,RI,LJ,RJ,M,p->orientation,LI2,RI2,LJ2,RJ2,XXmap[q->mapid2]->numsite[0],q->orientation,
		 qrystart,origX[qrystart]*Xscale,
		 qryend,origX[qryend]*Xscale,
		 qrystart2,origX[qrystart2]*Xscale,
		 qryend2,origX[qryend2]*Xscale,
		 XmapChim);
	  fflush(stdout);
	}

	if(!MultiMatches &&
	   !(qryend2 >= qryend + XmapUnique /* && qrystart2 >= qrystart + XmapChim */) &&
	   !(/* qryend >= qryend2 + XmapChim && */ qrystart >= qrystart2 + XmapUnique)){/* Not enough unique query part in 2nd alignment : delete the 2nd alignment (with worse logPV) */
	  delcnt++;
	  if(VERB>=2 && Ymap->id == 9 && Xmap->id == 2183){
	    printf("Overlapped QueryId=%lld(M=%d,len=%0.3f) with RefId=%lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f and %0.2f):Igoring 2nd alignment\n",
		   Xmap->id,M,origX[M+1]*Xscale,YYmap[rid]->id,N,origY[N+1]*Yscale,LI,RI,refstartpos,refendpos,qrystart,qryend,p->orientation,origX[qrystart],origX[qryend],
		   YYmap[q->mapid1]->id,N2,origY2[N2+1]*Yscale,LI2,RI2,refstartpos2,refendpos2,qrystart2,qryend2,q->orientation,origX[qrystart2],origX[qryend2],p->logPV,q->logPV);
	    fflush(stdout);
	  }
	  // NOTE : cannot actually delete the alignment since -resEstimate may need it to estimate res & resSD during last iteration in refalign() : just set the logPV below the threshold
	  Lalignment[k]->logPV = LogPvThreshold2 - 1.0;
	  continue;
	}
	int overlapstart = max(qrystart,qrystart2);
	int overlapend = min(qryend,qryend2);

	if(VERB>=2 && Ymap->id == 9 && Xmap->id == 2183  && overlapstart < overlapend && INDEL
	   && ((qrystart2 >= qrystart && qryend2 <= qryend) || (qrystart >= qrystart2 && qryend <= qryend2))){
	  printf("Fully overlapped Alignment of QueryId=%lld(len=%0.3f) with RefId=%lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f,%0.2f and %0.2f,%0.2f):overlap=%d..%d\n",
		 Xmap->id,origX[M+1]*Xscale,YYmap[rid]->id,N,origY[N+1]*Yscale,LI,RI,refstartpos,refendpos,qrystart,qryend,p->orientation,origX[qrystart],origX[qryend],
		 YYmap[q->mapid1]->id,N2,origY2[N2+1]*Yscale,LI2,RI2,refstartpos2,refendpos2,qrystart2,qryend2,q->orientation,origX[qrystart2],origX[qryend2],
		 p->logPV,p->logPV2,q->logPV,q->logPV2, overlapstart,overlapend);
	  fflush(stdout);
	}
	if(overlapstart < overlapend && INDEL){
	  /* flag the lower scoring alignment (q) query region that overlaps the higher scoring alignment's query region so indels are NOT output from that region */
	  /* HERE HERE : should check which alignment (p or q) is better in overlap region : either rescore logPV for this region OR just approximate with count of number of aligned labels. Not needed with InDelNoOverlap */
	  
	  q->suppress[overlapstart]++;
	  q->suppress[overlapend+1]--;
	  if(InDelNoOverlap && (RefSplit ? q->logPV2 : q->logPV) >= max(InDelNoOverlapT,LogPvThreshold)){/* also suppress output from higher scoring alignment, provided lower scoring alignment is above -T threshold and has at least XmapUnique
										       unique labels. NOTE : q->logPV could be lower than LogPvThreshold due to MultiMatch's -T2 threshold or -CutFlip */
	    if(InDelNoOverlap <= 1 || ((qryend2 >= qryend + XmapUnique || qrystart >= qrystart2 + XmapUnique) && 
				       (InDelNoOverlap <= 2 || !(q->orientation != p->orientation && max(refstartpos,refstartpos2) <= min(refendpos,refendpos2))))){
	      p->suppress[overlapstart]++;
	      p->suppress[overlapend+1]--;
	    }
	  }
	}
	if(VERB>=2 && Ymap->id == 9 && Xmap->id == 2183  && overlapstart < overlapend && INDEL){
	  if(InDelNoOverlap && (RefSplit ? q->logPV2 : q->logPV) >= max(InDelNoOverlapT,LogPvThreshold) && 
	     (InDelNoOverlap <= 1 || ((qryend2 >= qryend + XmapUnique || qrystart >= qrystart2 + XmapUnique) &&
				      (InDelNoOverlap <= 2 || !(q->orientation != p->orientation && max(refstartpos,refstartpos2) <= min(refendpos,refendpos2))))))
	    printf("Overlap Region in QueryId=%lld(len=%0.3f) with RefId=%lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f,%0.2f and %0.2f,%0.2f):suppress[%lu,%lu]=%d..%d(%0.3f..%0.3f):Lalignment[%lu,%lu]\n",
		 Xmap->id,origX[M+1]*Xscale,YYmap[rid]->id,N,origY[N+1]*Yscale,LI,RI,refstartpos,refendpos,qrystart,qryend,p->orientation,origX[qrystart],origX[qryend],
		   YYmap[q->mapid1]->id,N2,origY2[N2+1]*Yscale,LI2,RI2,refstartpos2,refendpos2,qrystart2,qryend2,q->orientation,origX[qrystart2],origX[qryend2],p->logPV, p->logPV2, q->logPV, q->logPV2,
		 i,k,overlapstart,overlapend,origX[overlapstart],origX[overlapend],i,k);
	  else
	    printf("Overlap Region in QueryId=%lld(len=%0.3f) with RefId=%lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f,%0.2f and %0.2f,%0.2f):suppress[%lu]=%d..%d(%0.3f..%0.3f):Lalignment[%lu,%lu]\n",
		   Xmap->id,origX[M+1]*Xscale,YYmap[rid]->id,N,origY[N+1]*Yscale,LI,RI,refstartpos,refendpos,qrystart,qryend,p->orientation,origX[qrystart],origX[qryend],
		   YYmap[q->mapid1]->id,N2,origY2[N2+1]*Yscale,LI2,RI2,refstartpos2,refendpos2,qrystart2,qryend2,q->orientation,origX[qrystart2],origX[qryend2],p->logPV, p->logPV2, q->logPV, q->logPV2,
		   k,overlapstart,overlapend,origX[overlapstart],origX[overlapend],i,k);
	  fflush(stdout);
	}
      } /* k = i+1 .. */
    } /* i = start .. end -1 */

    if(VERB){
      printf("Deleted %d fully overlapped alignments\n",delcnt);
      if(VERB>=2){
	printf("Alignments with partial overlaps:\n");
	fflush(stdout);
      }
    }

    int lasti = -1;
    for(size_t i = start; i < end; i++){
      Calign *p = Lalignment[i];
      if(!p){
	if(DEBUG>=1+RELEASE) assert(!(RefSplit && MultiMatches));
	continue;
      }
      int U = p->numpairs;
      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold))){
	if(DEBUG>=1+RELEASE) assert(!(RefSplit && MultiMatches));
	continue;
      }
      int mid = p->mapid2;
      Cmap *Xmap = XXmap[mid];
      int M = Xmap->numsite[0];
      double *X = Xmap->site[0];

      int cnt = 0;
      p->start = p->end = 0;
      for(int I = 1; I <= M+1; I++){
	p->suppress[I] += p->suppress[I-1];
	if(p->suppress[I]){
	  if(!p->start)
	    p->start = I;
	  p->end = I;
	  cnt++;
	}
      }
      if(cnt > 0){
	if(VERB>=2){
	  printf("Lalignment[%lu]:rid=%lld,qid=%lld,or=%d,M=%d,logPV=%0.2f(%0.2f): suppressed intervals include %d labels in X[%d..%d]=%0.3f..%0.3f:\n",
		 i,YYmap[p->mapid1]->id, Xmap->id,p->orientation,M, RefSplit ? p->logPV2 : p->logPV, p->logPV, cnt, p->start,p->end,X[p->start],X[p->end]);
	  if(VERB>=2){
	    for(int I = 1; I <= M; I++)
	      if(p->suppress[I])
		printf("\t X[%d]= %0.3f\n",I,X[I]);
	  }
	  fflush(stdout);
	}
	if(DEBUG>=1+RELEASE && RefSplit && !(p->logPV2 >= p->logPV)){
	  printf("Lalignment[%lu]:rid=%lld,qid=%lld,or=%d,M=%d,logPV=%0.2f(%0.2f): suppressed intervals include %d labels in X[%d..%d]=%0.3f..%0.3f:\n",
		 i,YYmap[p->mapid1]->id, Xmap->id,p->orientation,M, RefSplit ? p->logPV2 : p->logPV, p->logPV, cnt, p->start,p->end,X[p->start],X[p->end]);
	  printf("p->logPV=%0.6f,p->logPV2=%0.6f\n",p->logPV,p->logPV2);
	  fflush(stdout);
	  assert(p->logPV2 >= p->logPV);
	}
      }
      if(DEBUG && !InDelNoOverlap && (i==0 || (lasti >= 0 && p->mapid2 != Lalignment[lasti]->mapid2))) assert(cnt <= 0);// highest scoring matchgroup for each contig should not have any suppression regions 
      lasti = i;
    }
  }

  if(XmapChim){  /* check for chimeric query alignments */
    sprintf(filenameCH,"%s.chimeric", basename);
    if(!checkFile(filenameCH)){
      if(VERB){
	printf("Generating %s\n",filenameCH);
	fflush(stdout);
      }
      if((fpCH = fopen(filenameCH,"w"))==NULL){
	int eno = errno;
	char *err = strerror(eno);
	printf("failed to open file %s for writing chimeric contigs file:errno=%d:%s\n",filenameCH,eno,err);
	exit(1);
      }
      printversion(fpCH);
      //    if(DEBUG) assert(spots_filename[0]);
      if(full)
	fprintf(fpCH,"# Reference Maps From:\t%s_full_r.cmap\n", absbasename);
      else
	fprintf(fpCH,"# Reference Maps From:\t%s_r.cmap\n", absbasename);
      fprintf(fpCH,"# Query Maps From:\t%s\n",queryfilename);/* This may actually be a modified cmap file */
#ifndef WIN32
      if(filename[0] == '/')
	fprintf(fpCH,"# Xmap Entries From:\t%s\n",filename);    
      else
	fprintf(fpCH,"# Xmap Entries From:\t%s/%s\n",cwd,filename);    
#else
      fprintf(fpCH,"# Xmap Entries From:\t%s\n",filename);    
#endif
      
      fprintf(fpCH,"#h EntryID\tQryContigID\tRefContigID1\tRefContigID2\tQryStartPos1\tQryEndPos1\tRefStartPos1\tRefEndPos1\tOrientation1\tConfidence1\tQryStartPos2\tQryEndPos2\tRefStartPos2\tRefEndPos2\tOrientation2\tConfidence2\n");
      fprintf(fpCH,"#f int    \tint        \tint         \tint         \tfloat       \tfloat     \tfloat       \tfloat     \tstring      \tfloat      \tfloat       \tfloat     \tfloat       \tfloat     \tstring      \tfloat      \n");
    }

    int chimcnt = 0;/* total chimeric contigs */
    int chimcntS = 0;/* chimeric contigs within same reference contig : Requires -MultiMatches and may be spurious due to insufficiently large -hashdelta splitting the same alignment into 2 pieces */

    for(size_t i = start; i+1 < end; i++){
      Calign *p = Lalignment[i];
      if(!p)
	continue;
      int U = p->numpairs;
      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold)))
	continue;
      
      int rid = p->mapid1;
      int mid = p->mapid2;

      if(!((i <= start || mid != Lalignment[i-1]->mapid2)))
	continue;/* NOTE : only check chimeric alignments where the first alignment has the highest log10PV, so we only count one chimerism per query contig */

      Cmap *Ymap = YYmap[rid];
      Cmap *Xmap = XXmap[mid];
      Cmap *origmap = Xmap;
      while(origmap->origmap){
	Xmap->origmap = origmap = origmap->origmap;
	if(origmap->origmap){
	  Xmap->left[0] += max(0,origmap->left[0]-1);
	  Xmap->right[0] += max(0,origmap->left[0]-1);
	}
      }
      int N = Ymap->numsite[0];
      int M = Xmap->numsite[0];
      FLOAT *Y = Ymap->site[0];
      FLOAT *X = Xmap->site[0];
      FLOAT *origX = Xmap->origmap ? Xmap->origmap->site[0] : X;
      FLOAT *origY = Ymap->origmap ? Ymap->origmap->site[0] : Y;
      if(DEBUG>=2) assert(Xmap->origmap==0 || Xmap->origmap->origmap==0);/* If this fails fix origmap link in refalign() to point to root map */
      if(DEBUG>=2) assert(Ymap->origmap==0 || Ymap->origmap->origmap==0);/* If this fails fix origmap link (same as for Xmap above) to point to root map */
      int Yshift = Ymap->origmap ? max(1,Ymap->left[0]) - 1 : 0;
      int Xshift = Xmap->origmap ? max(1,Xmap->left[0]) - 1 : 0;
      int LI = p->sites1[0];
      int LJ = p->sites2[0];
      int RI = p->sites1[U-1];
      int RJ = p->sites2[U-1];
      if(DEBUG) assert(LJ <= RJ);
      int qrystart = p->orientation ? Xshift + M + 1 - RJ : Xshift + LJ;
      int qryend   = p->orientation ? Xshift + M + 1 - LJ : Xshift + RJ;
      double refstartpos = origY[Yshift + LI]*Yscale; /* RefStartPos */
      double refendpos   = origY[Yshift + RI]*Yscale; /* RefEndPos */
      if(DEBUG) assert(qrystart <= qryend);

      int typeMax = (MultiMatches && XmapChimSpacing) ? 1 : 0;
      for(int type = 0; type <= typeMax; type++){/* first look for cross chromosome chimeric contigs, then same chromosome */
	int found = 0;/* flag if we found a chimeric match to alignment p : in that case break out of type = 0 .. typeMax loop */
	for(size_t k = i+1; k < end; k++){
	  Calign *q = Lalignment[k];
	  if(DEBUG) assert(q);
	  if(!q)
	    break;
	  if(mid != q->mapid2)
	    break;
	  if(q->mapid1 == rid && !type){	  // If MultiMatches need to check for chimeric matches with same rid
	    if(DEBUG) assert(NoSplit <= 1 || SecondBest || MultiMatches);
	    continue;
	  }
	  if(q->mapid1 != rid && type)
	    continue;
	  int U2 = q->numpairs;
	  if(U2 <= 0 || (RefSplit ? ((!REFSPLIT_FIX && q->score <= ScoreThreshold2) || q->logPV2 <= LogPvThreshold2 || U2 < AlignedSiteThreshold2) 
			 : ((!REFSPLIT_FIX && q->score <= ScoreThreshold) || q->logPV <= LogPvThreshold || U2 < AlignedSiteThreshold)))
	    continue;
	  Cmap *Ymap2 = YYmap[q->mapid1];
	  int N2 = Ymap2->numsite[0];
	  FLOAT *Y2 = Ymap2->site[0];
	  FLOAT *origY2 = Ymap2->origmap ? Ymap2->origmap->site[0] : Y2;
	  if(DEBUG>=2) assert(Ymap2->origmap==0 || Ymap2->origmap->origmap==0);/* If this fails fix origmap link (same as for Xmap above) to point to root map */
	  int Yshift2 = Ymap2->origmap ? max(1,Ymap2->left[0]) - 1 : 0;

	  /* check if alignments p & q have at least XmapChim distinct sites on q->mapid2 */
	  int LI2 = q->sites1[0];
	  int LJ2 = q->sites2[0];
	  int RJ2 = q->sites2[q->numpairs-1];
	  int RI2 = q->sites1[q->numpairs-1];
	  if(DEBUG && !(LJ2 <= RJ2)){
	    printf("i=%lu,k=%lu:rid=%d,mid=%d,q->mapid1=%d,q->mapid2=%d:LJ=%d,RJ=%d,LJ2=%d,RJ2=%d\n",
		   i,k,rid,mid,q->mapid1,q->mapid2,LJ,RJ,LJ2,RJ2);
	    fflush(stdout);
	    assert(LJ2 <= RJ2);
	  }
	  int qrystart2 = q->orientation ? Xshift + M + 1 - RJ2 : Xshift + LJ2;
	  int qryend2   = q->orientation ? Xshift + M + 1 - LJ2 : Xshift + RJ2;
	  double refstartpos2 = origY2[Yshift2 + LI2]*Yscale; /* RefStartPos */
	  double refendpos2   = origY2[Yshift2 + RI2]*Yscale; /* RefEndPos */
	  if(DEBUG) assert(qrystart2 <= qryend2);
	
	  if(q->mapid1 == rid){/* make sure alignments are separated in reference or Query by at least XmapChimSpacing */
	    if(DEBUG) assert(MultiMatches && type && XmapChimSpacing);
	    /* minimum gap in either must be at least XmapChim * (XmapLength[0] * PixelLen / origPixelLen) / XmapSites[0] */
	    double refgap = max(refstartpos2 - refendpos, refstartpos - refendpos2);
	    double querygap = max(origX[qrystart2] - origX[qryend], origX[qrystart] - origX[qryend2]);
	    if(max(refgap,querygap) < XmapChimSpacing){
	      if(VERB>=2 && (qryend2 >= qryend + XmapChim || qrystart >= qrystart2 + XmapChim)){
		printf("Skipping Chimeric QueryId=%lld(len=%0.3f) with RefId=%lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f and %0.2f):\n",
		       Xmap->id,origX[M+1]*Xscale,YYmap[rid]->id,N,origY[N+1]*Yscale,LI,RI,refstartpos,refendpos,qrystart,qryend,p->orientation,origX[qrystart],origX[qryend],
		       YYmap[q->mapid1]->id,N2,origY2[N2+1]*Yscale,LI2,RI2,refstartpos2,refendpos2,qrystart2,qryend2,q->orientation,origX[qrystart2],origX[qryend2],p->logPV,q->logPV);
		fflush(stdout);
	      }
	      continue;
	    }
	    if(VERB>=2 && (qryend2 >= qryend + XmapChim || qrystart >= qrystart2 + XmapChim)){
	      printf("QueryId=%lld,RefId=%lld,or=%d,%d:refgap=%0.3f,querygap=%0.3f,XmapChimSpacing=%0.3f\n",Xmap->id,YYmap[rid]->id,p->orientation,q->orientation,refgap,querygap,XmapChimSpacing);
	      fflush(stdout);
	    }
	  }

	  if(qryend2 >= qryend + XmapChim || qrystart >= qrystart2 + XmapChim){
	    chimcnt++;
	    if(q->mapid1 == rid)
	      chimcntS++;
	    if(VERB && end - start < 100000){
	      printf("Chimeric QueryId=%lld(len=%0.3f) with RefId=%lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f and %0.2f):\n",
		     Xmap->id,origX[M+1]*Xscale,YYmap[rid]->id,N,origY[N+1]*Yscale,LI,RI,refstartpos,refendpos,qrystart,qryend,p->orientation,origX[qrystart],origX[qryend],
		     YYmap[q->mapid1]->id,N2,origY2[N2+1]*Yscale,LI2,RI2,refstartpos2,refendpos2,qrystart2,qryend2,q->orientation,origX[qrystart2],origX[qryend2],p->logPV,q->logPV);
	      //	      printf("    i=%lu,k=%lu:rid1=%d,rid2=%d\n",i,k,p->mapid1,q->mapid1);
	      fflush(stdout);
	    }
	    if(fpCH)
	      fprintf(fpCH,"%6d    \t%10lld\t%12lld\t%12lld\t%12.1f\t%10.1f\t%12.1f\t%10.1f\t%12c\t%11.2f\t%12.1f\t%10.1f\t%12.1f\t%10.1f\t%12c\t%10.2f\n",
		      chimcnt,Xmap->id,YYmap[p->mapid1]->id,YYmap[q->mapid1]->id,origX[qrystart],origX[qryend],refstartpos,refendpos,p->orientation ? '-' : '+',p->logPV,
		      origX[qrystart2],origX[qryend2],refstartpos2,refendpos2,q->orientation ? '-' : '+', q->logPV);

	    /* skip ahead i until last alignment with mapid2 equals p->mapid2 == mid */
	    for(k++;k < end;k++)
	      if(Lalignment[k]->mapid2 != mid)
		break;
	    if(DEBUG) assert(Lalignment[k-1]->mapid2 == mid);
	    i = k - 1;
	    found = 1;
	    break;
	  }
	  if(VERB>=2 && ((qrystart2 >= qrystart && qryend2 <= qryend) ||
			 (qrystart >= qrystart2 && qryend <= qryend2))){
	    printf("Duplicate Region in QueryId=%lld(len=%0.3f) with RefId=%lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f and %0.2f):\n",
		   Xmap->id,origX[M+1]*Xscale,YYmap[rid]->id,N,origY[N+1]*Yscale,LI,RI,refstartpos,refendpos,qrystart,qryend,p->orientation,origX[qrystart],origX[qryend],
		   YYmap[q->mapid1]->id,N2,origY2[N2+1]*Yscale,LI2,RI2,refstartpos2,refendpos2,qrystart2,qryend2,q->orientation,origX[qrystart2],origX[qryend2],p->logPV,q->logPV);
	    fflush(stdout);
	  }
	}// for k = i+1 .. end-1
	if(found)/* If type==0 : already found a chimeric contig across chromosomes, no need to look for same chromosome chimeric contig */
	  break;
      }// for type = 0 .. 1
    }
    if(VERB){
      printf("Found %d Chimeric Query contigs (including %d across multiple chromosomes)\n",chimcnt, chimcnt - chimcntS);
      fflush(stdout);
    }
  }

  if(fpCH)
    FILEclose(fpCH);
  
  if(VERB>=2){
    printf("Alignments before sorting:\n");
    for(size_t i = start; i < end; i++){
      Calign *p = Lalignment[i];
      if(!p)
	continue;
      int U = p->numpairs;
      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold)))
	continue;
      printf("Lalignment[%lu]:refid=%d,mapid=%d(id=%lld),or=%d:numpairs=%d,score=%0.3f,logPV=%0.2f,chimpair=%d\n",
	     i, p->mapid1, p->mapid2, Gmap[p->mapid2]->id, p->orientation, U, p->score, p->logPV,p->chimpair);
    }
    fflush(stdout);
  }

  //  qsort(&Lalignment[start],end - start,sizeof(Calign *), (intcmp*)Id2SiteInc);// HERE HERE : sort by id2,logPV,id1

  /* sort alignments in order of id1, sites1[0],sites1[N+1], id2, sites2[0],sites2[M+1] (null alignments last) */
  qsort(&Lalignment[start],end - start,sizeof(Calign *), (intcmp*)SiteInc);


  if(VERB>=2){
    printf("After sorting alignments for xmap output:start=%lu,end=%lu\n",start,end);
    for(size_t i = start; i < end; i++){
      Calign *p = Lalignment[i];
      if(!p)
	continue;
      int U = p->numpairs;

      if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		  : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold)))
	continue;
      
      int rid = p->mapid1;
      int mid = p->mapid2;

      printf("Lalignment[%lu]: rid=%d(id=%lld),mid=%d(id=%lld):logPV=%0.2f\n",i,rid,YYmap[rid]->id,mid,XXmap[mid]->id,p->logPV);
    }
    fflush(stdout);
  }

  /* sort alignments in order of LogPV, id2, id1, sites1[0],sites1[numpairs-1] (null alignments last) */
  if(0)// for debugging
    qsort(&Lalignment[start],end - start,sizeof(Calign *), (intcmp*)LogPVSiteInc);
  
  /* sort alignments in order of id2, id1, sites1[0],sites1[numpairs-1] (null alignments last) */
  if(0)// for debugging
    qsort(&Lalignment[start],end - start,sizeof(Calign *), (intcmp*)Site2Inc);

  if((PairSplit || (MultiMatches && RefSplit)) && (dosvdetect || doscaffold)){ /* see parameters.{h,cpp} */
    if(DEBUG) assert(end - start <= MASK(31));
    delete [] xmapentries;
    maxxmapentries = end - start;
    xmapentries = new xmapEntry*[maxxmapentries]; /* allocate for numaligns entries, though actual number will be <= this */
    numxmapentries = 0;
  }

  if(VERB>=2){
    printf("output_xmap(%s): align_start=%lu,align_end=%lu (before output_xmapIO)\n",filename,align_start,align_end);
    fflush(stdout);
  }

  /* Perform IO to .xmap (and .indel, if needed) */
  output_xmapIO(filename, Lalignment, start, end, basename, queryfilename, INDEL, xmapentries, Xscale, Yscale, full);

  if(DEBUG) assert(Lalignment != 0 && Lalignment != alignment);
  delete [] Lalignment;
}

void output_xmap2(char *basename, size_t align_start, size_t align_end)
{
  if(DEBUG) assert(colors==2);

  if(strstr(basename,"/dev/null"))
    return;

  char filename[PATH_MAX];
  strcpy(filename,basename);
  /* no need to strip suffix : for default basename this was already done in RefAligner.cpp */
  if(1){
    int i = strlen(filename);
    strcpy(&filename[i],".xmap");
  }
  if(checkFile(filename))
    return;

  if(VERB){
    printf("Generating %s\n",filename);
    fflush(stdout);
  }
  FILE *fp;
  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("failed to open file %s for writing xmap file:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  size_t blen = 10*1024*1024;/* matches scif_io buffer size */
  char *write_buffer = (char *)malloc(blen);
  if(!write_buffer){
    printf("output_xmap2:malloc(%llu) failed\n",(unsigned long long)blen);
    exit(1);
  }
  setbuffer(fp, write_buffer, blen);

  /* write out commandline */
  printversion(fp);

  if(RefSplit){
    printf("RefSplit not yet implemented with 2 colors\n");
    fflush(stdout);exit(1);
  }

  /* HERE HERE : update changes from output_xmap */

  /* remove below-threshold alignments */
  for(size_t i= align_start; i < align_end; i++){
    register Calign *p = alignment[i];
    if(!p)
      continue;
    if(p->numpairs <= 0 || (!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold){
      delete [] alignment[i];
      alignment[i] = 0;
    }
  }

  /* get current directory pathname */
  char cwd[PATH_MAX];
  if(getcwd(cwd,PATH_MAX) == NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("getcwd() failed (errno=%d):%s\n",eno,err);
    exit(1);
  }

  /* convert basename into absolute pathname */
#ifndef WIN32
  char absbasename[PATH_MAX];
  if(basename[0]=='/')
    strcpy(absbasename,basename);
  else
    sprintf(absbasename,"%s/%s",cwd,basename);
  if(vfixx_filename[0] && vfixx_filename[0][0] != '/'){/* convert vfixx_filename[0] into an absolute pathname */
    char absname[PATH_MAX];
    sprintf(absname,"%s/%s",cwd,vfixx_filename[0]);
    free(vfixx_filename[0]);
    vfixx_filename[0] = strdup(absname);
  }
#else
  char *absbasename = basename;
#endif

  /* Also output indel SVs (if -indel) */
  FILE *fpSV = NULL;
  char filenameSV[PATH_MAX];
  char *write_bufferSV = NULL;
  if(INDEL){
    sprintf(filenameSV,"%s.indel",basename);
    if(!checkFile(filenameSV)){
      if(VERB){
	printf("Generating %s\n",filenameSV);
	fflush(stdout);
      }
      if((fpSV = fopen(filenameSV,"w"))==NULL){
	int eno = errno;
	char *err = strerror(eno);
	printf("failed to open file %s for writing indel file:errno=%d:%s\n",filenameSV,errno,err);
	exit(1);
      }
    }
  }

  if(fpSV != NULL){
    write_bufferSV = (char *)malloc(blen);
    if(!write_bufferSV){
      printf("output_xmap2:malloc(%llu) for write_bufferSV failed\n",(unsigned long long)blen);
      exit(1);
    }
    setbuffer(fpSV, write_bufferSV, blen);

    printversion(fpSV);

    fprintf(fpSV,"# SMAP File Version:\t0.0\n"); //  (0.0 was indels only; 0.1 has newer types)
    if(spots_filename[0] || (svcheck && numrefmaps > 0)){/* Reference Aligner : refer to condensed version of reference CMAP */
      fprintf(fpSV,"# Reference Maps From:\t%s_r.cmap\n", absbasename);
      fprintf(fpSV,"# Query Maps From:\t%s\n",vfixx_filename[0]);/* This may actually be a modified cmap file */
    } else {
      if(mres > 0.001 || mresSD > 0.0 || fabs(PixelLen - origPixelLen) > 1e-12){
	fprintf(fpSV,"# Reference Maps From:\t%s_r.cmap\n", absbasename);
	fprintf(fpSV,"# Query Maps From:\t%s_q.cmap\n", absbasename);
      } else {
	fprintf(fpSV,"# Reference Maps From:\t%s_r.cmap", vfixx_filename[0]);
	for(register int i = 1; i < RefIndex; i++)
#ifndef WIN32
	  if(vfixx_filename[i][0] == '/')
	    fprintf(fpSV,",%s",vfixx_filename[i]);
	  else
	    fprintf(fpSV,",%s/%s",cwd, vfixx_filename[i]);
#else
	    fprintf(fpSV,",%s",vfixx_filename[i]);
#endif
	fprintf(fpSV,"\n");
#ifndef WIN32
	if(vfixx_filename[QueryIndex][0] == '/')
	  fprintf(fpSV,"# Query Maps From:\t%s",vfixx_filename[QueryIndex]);
	else
	  fprintf(fpSV,"# Query Maps From:\t%s/%s",cwd,vfixx_filename[QueryIndex]);
#else
	  fprintf(fpSV,"# Query Maps From:\t%s",vfixx_filename[QueryIndex]);
#endif
	for(register int i = QueryIndex+1; i < num_files;i++)
#ifndef WIN32
	  if(vfixx_filename[i][0] == '/')
	    fprintf(fpSV,",%s",vfixx_filename[i]);
	  else
	    fprintf(fpSV,",%s/%s",cwd, vfixx_filename[i]);
#else
	    fprintf(fpSV,",%s",vfixx_filename[i]);
#endif
	fprintf(fpSV,"\n");
      }
    }
    fprintf(fpSV,"# Xmap Entries From:\t%s\n",filename);    
    fprintf(fpSV,"#h SmapEntryID\tQryContigID\tRefcontigID1\tRefcontigID2\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tType\tXmapID1\tXmapID2\n");
    fprintf(fpSV,"#f int        \tint        \tint         \tint         \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat \tstring \tint\tint\n");
  }

  /* remove below-threshold or deallocated alignments */
  if(RepeatShift > 0.0 || FirstAlignments > 0){
    if(DEBUG) assert(align_start == 0 && align_end == numaligns);
    size_t aligncnt = 0;
    for(register size_t i=0; i < numaligns; i++){
      register Calign *p = alignment[i];
      if(!p)
	continue;
      alignment[i] = 0;
      if(p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold || p->numpairs <= 0){
	delete [] p;
	continue;
      }
      if(DEBUG) assert(alignment[i] == 0);
      alignment[aligncnt++] = p;
    }
    align_end = numaligns = aligncnt;
  }

  Calign **Lalignment = alignment;
  size_t start = align_start;
  size_t end = align_end;

  if(SecondBest && !SinglelineXmap){/* interleave 2nd best alignments after best alignments */
    /* NOTE : the 2nd best alignments are marked with chimpair = -1, others with chimpair = 0 */
    Lalignment = new Calign*[(align_end - align_start)*2];
    start = 0;
    end = 0;
    size_t SecondCnt = 0;
    size_t SplitCnt = 0;
    if(VERB>=2){
      printf("ScoreThreshold=%0.6f,LogPvThreshold=%0.2f,AlignedSiteThreshold=%d\n",ScoreThreshold,LogPvThreshold,AlignedSiteThreshold);
      printf("SecondBest alignments:\n");
      fflush(stdout);
    }

    for(size_t i = align_start; i < align_end; i++){
      Calign *p = alignment[i];
      if(!p)
	continue;
      int U = p->numpairs;
      p->chimpair = 0;
      if(U<=0 || p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold){
	if(VERB>=2){
	  printf("alignment[%lu/%lu]= %p:refid=%d,mapid=%d,or=%d:numpairs=%d,score=%0.3f,logPV=%0.2f, align2= %p: skipping\n",i, align_end, p, 
		 p->mapid1,p->mapid2,p->orientation,U, p->score, p->logPV, p->align2);
	  fflush(stdout);
	}
	continue;
      }
      if(Gmap[p->mapid2]->origmap){
	SplitCnt++;
	if(VERB>=2){
	  printf("Lalignment[%lu]= %p: numpairs=%d,score=%0.3f,logPV=%0.2f: refid=%d, split mapid=%d(orig=%d,id=%lld),or=%d\n",
		 end, p, U, p->score, p->logPV, p->mapid1,p->mapid2,Gmap[p->mapid2]->origmap->mapid,Gmap[p->mapid2]->id,p->orientation);
	  fflush(stdout);
	}
      }
      Lalignment[end++] = p;

      if(p->align2){
	Calign *q = p->align2;
	if(DEBUG) assert(q->numpairs > 0);
	if(q->score <= ScoreThreshold2 || q->logPV <= LogPvThreshold2){
	  if(VERB>=2){
	    printf("Lalignment[%lu]= %p:refid=%d,mapid=%d,or=%d:numpairs=%d,score=%0.3f,logPV=%0.2f, p->align2= %p: numpairs=%d,score=%0.3f,logPV=%0.2f (skipping align2)\n",
		   end-1, p, p->mapid1, p->mapid2, p->orientation, U, p->score, p->logPV, p->align2, q->numpairs, q->score, q->logPV);
	    fflush(stdout);
	  }
	  delete [] p->align2; p->align2 = NULL;
	  continue;
	}
	if(DEBUG) assert(q->align2 == NULL);
	if(VERB>=2){
	  printf("Lalignment[%lu]= %p: numpairs=%d,score=%0.3f,logPV=%0.2f: Lalignment[%lu]->p->align2= %p: refid=%d,mapid=%d(%lld),or=%d,numpairs=%d,score=%0.3f,logPV=%0.2f\n",
		 end-1, p, U, p->score, p->logPV, end, q, q->mapid1,q->mapid2,Gmap[q->mapid2]->id, q->orientation, q->numpairs, q->score, q->logPV);
	  fflush(stdout);
	}
	q->chimpair = -1;
	Lalignment[end++] = q;
	SecondCnt++;
      }
    }
    //    fprintf(fp,"# SecondID refers to XmapEntryID of Second Best alignment for same query and reference (typically next line): 0 if none, -1 if current line is Second Best alignment\n");

    if(VERB){
      printf("%lu Xmap alignments (including %lu splits) + %lu SecondBest alignments\n", end - SecondCnt, SplitCnt, SecondCnt);
      fflush(stdout);
    }
  }

  double Xscale = 1.0;
  double Yscale = 1.0;
  if(!(spots_filename[0] || (svcheck && numrefmaps > 0)) && fabs(PixelLen - origPixelLen) > 1e-12 /* NEW */){
    Yscale = Xscale = origPixelLen/PixelLen;
    if(PairSplit)/* -i inputs that correspond to reference is not scaled by bpp/500 */
      Yscale = 1.0;
  }

  // FILE *fpCH = NULL;
  //   char filenameCH[PATH_MAX];
  //  char *write_bufferCH = NULL;

  if(XmapChim || XmapUnique)  /* sort alignments in order of id2, logPV  (null alignments last) */
    qsort(&Lalignment[start],end - start,sizeof(Calign *), (intcmp*)Id2SiteInc);

  if(XmapChim){  /* check for chimeric query alignments */
    printf("-xmapChim not yet implemented for 2 colors\n");
    exit(1);
  }

  /* keep track of query map regions with overlapped alignments : suppress part of the alignment with the worse logPV in the overlap region when reporting indels */
  if(INDEL || XmapUnique)
    for(size_t i = start; i < end; i++){
      Calign *p = Lalignment[i];
      if(!p)
	continue;
      // HERE HERE : update to handle multiple regions as in 1-color code
      for(int c = 0; c < colors; c++)
	p[1+c].start = p[1+c].end = 0;
    }

  if(XmapUnique){  /* delete alignments if there is another alignment of the same query map with better logPV and the first alignment does not have at least XmapUnique unique sites on the query map */
    int delcnt = 0;

    for(size_t i = start; i+1 < end; i++){
      Calign *p = Lalignment[i];
      if(!p)
	continue;
      if(p->numpairs <=0 || p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold)
	continue;
      if(SecondBest && p->chimpair)
	continue;/* skip 2nd best alignments */
      
      int rid = p->mapid1;
      int mid = p->mapid2;
      Cmap *Ymap = YYmap[rid];
      Cmap *Xmap = XXmap[mid];
      if(VERB>=2){
	printf("Checking Lalignment[%lu]:mid=%lld,rid=%lld:numpairs=%d,score=%0.6f,logPV=%0.2f\n",i,Xmap->id,Ymap->id,p->numpairs,p->score,p->logPV);
	fflush(stdout);
      }

      Cmap *origmap = Xmap;
      while(origmap->origmap){
	Xmap->origmap = origmap = origmap->origmap;
	if(origmap->origmap){/* otherwise origmap->left[] is not defined */
	  for(int c = 0; c < colors; c++){
	    Xmap->left[c] += max(0,origmap->left[c]-1);
	    Xmap->right[c] += max(0,origmap->left[c]-1);
	  }
	}
      }
      int *NN = Ymap->numsite;
      int *MM = Xmap->numsite;
      FLOAT **YY = Ymap->site;
      FLOAT **XX = Xmap->site;
      if(DEBUG>=2) assert(Xmap->origmap==0 || Xmap->origmap->origmap==0);/* If this fails fix origmap link in refalign() to point to root map */
      if(DEBUG>=2) assert(Ymap->origmap==0 || Ymap->origmap->origmap==0);/* If this fails fix origmap link (same as for Xmap above) to point to root map */
      FLOAT **origXX = Xmap->origmap ? Xmap->origmap->site : XX;
      FLOAT **origYY = Ymap->origmap ? Ymap->origmap->site : YY;
      int Yshift[2],Xshift[2],LI[2],LJ[2],RI[2],RJ[2],qrystart[2],qryend[2],U[2];
      for(int c = 0; c < 2; c++){
	if(DEBUG) assert(p[1+c].numpairs > 0);/* HERE : if this fails, need to update the code */
	U[c] = p[1+c].numpairs;
	Yshift[c] = Ymap->origmap ? max(1,Ymap->left[c]) - 1 : 0;
	Xshift[c] = Xmap->origmap ? max(1,Xmap->left[c]) - 1 : 0;
	LI[c] = p[1+c].sites1[0];
	LJ[c] = p[1+c].sites2[0];
	RI[c] = p[1+c].sites1[U[c]-1];
	RJ[c] = p[1+c].sites2[U[c]-1];
	if(DEBUG) assert(LJ[c] <= RJ[c]);
	qrystart[c] = p->orientation ? Xshift[c] + MM[c] + 1 - RJ[c] : Xshift[c] + LJ[c];
	qryend[c]   = p->orientation ? Xshift[c] + MM[c] + 1 - LJ[c] : Xshift[c] + RJ[c];
	if(DEBUG) assert(qrystart[c] <= qryend[c]);
      }
      double qrystartpos = min(origXX[0][qrystart[0]],origXX[1][qrystart[1]] * Xscale);
      double qryendpos = max(origXX[0][qryend[0]],origXX[1][qryend[1]]) * Xscale;
      double refstartpos = min(origYY[0][Yshift[0] + LI[0]], origYY[1][Yshift[1] + LI[1]]) * Yscale; /* RefStartPos */
      double refendpos   = max(origYY[0][Yshift[0] + RI[0]], origYY[1][Yshift[1] + RI[1]]) * Yscale; /* RefEndPos */

      for(size_t k = i+1; k < end; k++){
	Calign *q = Lalignment[k];
	if(!q)
	  continue;
	if(mid != q->mapid2)
	  break;
	if(VERB>=3){
	  printf("Checking alignments i=%lu,k=%lu:mid1=%lld,mid2=%lld,rid1=%lld,rid2=%lld:q->numpairs=%d,score=%0.6f,logPV=%0.2f\n",
		 i,k,Xmap->id,XXmap[q->mapid2]->id,Ymap->id,YYmap[q->mapid1]->id,q->numpairs,q->score,q->logPV);
	  fflush(stdout);
	}
	if(q->numpairs <= 0 || q->score <= ScoreThreshold || q->logPV <= LogPvThreshold || q->numpairs < AlignedSiteThreshold)
	  continue;
	if(SecondBest && q->chimpair < 0)
	  continue;/* skip 2nd best alignments */
	Cmap *Ymap2 = YYmap[q->mapid1];
	int *NN2 = Ymap2->numsite;
	FLOAT **YY2 = Ymap2->site;
	if(DEBUG>=2) assert(Ymap2->origmap==0 || Ymap2->origmap->origmap==0);/* If this fails fix origmap link (same as for Xmap above) to point to root map */
	FLOAT **origYY2 = Ymap2->origmap ? Ymap2->origmap->site : YY2;

	/* check if alignments p & q have at least XmapUnique distinct sites on both p->mapid1 and q->mapid1 */

	int Yshift2[2],LI2[2],LJ2[2],RJ2[2],RI2[2],qrystart2[2],qryend2[2];
	for(int c = 0; c < colors; c++){
	  if(DEBUG) assert(q[1+c].numpairs > 0);// HERE : if this fails, the code must be changed 
	  Yshift[c] = Ymap2->origmap ? max(1,Ymap2->left[c]) - 1 : 0;
	  LI2[c] = q[1+c].sites1[0];
	  LJ2[c] = q[1+c].sites2[0];
	  RJ2[c] = q[1+c].sites2[q[1+c].numpairs-1];
	  RI2[c] = q[1+c].sites1[q[1+c].numpairs-1];
	  if(DEBUG && !(LJ2[c] <= RJ2[c])){
	    printf("i=%lu,k=%lu:rid=%d,mid=%d,q->mapid1=%d,q->mapid2=%d,c=%d:LJ[c]=%d,RJ[c]=%d,LJ2[c]=%d,RJ2[c]=%d\n",
		   i,k,rid,mid,q->mapid1,q->mapid2,c,LJ[c],RJ[c],LJ2[c],RJ2[c]);
	    fflush(stdout);
	    assert(LJ2[c] <= RJ2[c]);
	  }
	  qrystart2[c] = q->orientation ? Xshift[c] + MM[c] + 1 - RJ2[c] : Xshift[c] + LJ2[c];
	  qryend2[c]   = q->orientation ? Xshift[c] + MM[c] + 1 - LJ2[c] : Xshift[c] + RJ2[c];
	  if(DEBUG) assert(qrystart2[c] <= qryend2[c]);
	}
	double qrystartpos2 = min(origXX[0][qrystart2[0]],origXX[1][qrystart2[1]]) * Xscale;
	double qryendpos2 = max(origXX[0][qryend2[0]],origXX[1][qryend2[1]]) * Xscale;
	double refstartpos2 = min(origYY2[0][Yshift2[0] + LI2[0]], origYY2[1][Yshift2[1] + LI2[1]]) * Yscale; /* RefStartPos */
	double refendpos2   = max(origYY2[0][Yshift2[0] + RI2[0]], origYY2[1][Yshift2[1] + RI2[1]]) * Yscale; /* RefEndPos */
	
	if(VERB>=3){
	  for(int c = 0; c < colors; c++)
	    printf("   c=%d:LI=%d,RI=%d,LJ=%d,RJ=%d,M=%d,or=%d;LI2=%d,RI2=%d,LJ2=%d,RJ2=%d,M2=%d,or2=%d;qrystart=%d(%0.3f),qryend=%d(%0.3f),qrystart2=%d(%0.3f),qryend=%d(%0.3f),XmapChim=%d\n",
		   c,LI[c],RI[c],LJ[c],RJ[c],MM[c],p->orientation,LI2[c],RI2[c],LJ2[c],RJ2[c],XXmap[q->mapid2]->numsite[c],q->orientation,
		   qrystart[c],origXX[c][qrystart[c]]*Xscale,
		   qryend[c],origXX[c][qryend[c]]*Xscale,
		   qrystart2[c],origXX[c][qrystart2[c]]*Xscale,
		   qryend2[c],origXX[c][qryend2[c]]*Xscale,
		   XmapChim);
	  fflush(stdout);
	}

	if(!(qryend2[0]+qryend2[1] >= qryend[0]+qryend[1] + XmapChim /* && qrystart2[0]+qrystart2[1] >= qrystart[0]+qrystart[1] + XmapChim */) &&
	   !(/* qryend[0]+qryend[1] >= qryend2[0]+qryend2[1] + XmapChim && */ qrystart[0]+qrystart[1] >= qrystart2[0]+qrystart[1] + XmapChim)){/* Not enough unique query part in 2nd alignment : delete the 2nd alignment (with worse logPV) */
	  delcnt++;
	  if(VERB>=2){
	    printf("Overlapped QueryId=%lld(M=%d,%d,len=%0.3f) with RefId=%lld(N=%d,%d,len=%0.3f:R=%d,%d..%d,%d:%0.3f..%0.3f;Q=%d,%d..%d,%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,%d,len=%0.3f:R=%d,%d..%d,%d:%0.3f..%0.3f;Q=%d,%d..%d,%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f and %0.2f):Igoring 2nd alignment\n",
		   Xmap->id,MM[0],MM[1],origXX[0][MM[0]+1]*Xscale,YYmap[rid]->id,NN[0],NN[0],origYY[0][NN[0]+1]*Yscale,LI[0],LI[1],RI[0],RI[1],refstartpos,refendpos,
		   qrystart[0],qrystart[1],qryend[0],qryend[1],p->orientation,qrystartpos,qryendpos,
		   YYmap[q->mapid1]->id,NN2[0],NN2[1],origYY2[0][NN2[0]+1]*Yscale,LI2[0],LI2[1],RI2[0],RI2[1],refstartpos2,refendpos2,
		   qrystart2[0],qrystart2[1],qryend2[0],qryend2[1],q->orientation,qrystartpos2,qryendpos2,p->logPV,q->logPV);
	    fflush(stdout);
	  }
	  // NOTE : cannot actually delete the alignment since -resEstimate may need it to estimate res & resSD during last iteration : just set the logPV below the threshold
	  //	  delete [] Lalignment[k];
	  //	  Lalignment[k] = q = 0;
	  Lalignment[k]->logPV = LogPvThreshold - 1.0;
	  continue;
	}

	if(VERB>=2 && ((qrystartpos2 >= qrystartpos - 1e-6 && qryendpos2 <= qryendpos + 1e-6) ||
		       (qrystartpos >= qrystartpos2 - 1e-6 && qryendpos <= qryendpos2 + 1e-6))){
	  printf("Duplicate Region in QueryId=%lld(len=%0.3f) with RefId=%lld(N=%d,%d,len=%0.3f:R=%d,%d..%d,%d:%0.3f..%0.3f;Q=%d,%d..%d,%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,%d,len=%0.3f:R=%d,%d..%d,%d:%0.3f..%0.3f;Q=%d,%d..%d,%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f and %0.2f):\n",
		 Xmap->id,origXX[0][MM[0]+1]*Xscale,YYmap[rid]->id,NN[0],NN[1],origYY[0][NN[0]+1]*Yscale,LI[0],LI[1],RI[0],RI[1],refstartpos,refendpos,
		 qrystart[0],qrystart[1],qryend[0],qryend[1],p->orientation,qrystartpos,qryendpos,
		 YYmap[q->mapid1]->id,NN2[0],NN2[1],origYY2[0][NN2[0]+1]*Yscale,LI2[0],LI2[1],RI2[0],RI2[1],refstartpos2,refendpos2,
		 qrystart2[0],qrystart2[1],qryend2[0],qryend2[1],q->orientation,qrystartpos2,qryendpos2,p->logPV,q->logPV);
	  fflush(stdout);
	}
	for(int c = 0; c < colors; c++){
	  int overlapstart = max(qrystart[c],qrystart2[c]);
	  int overlapend = min(qryend[c],qryend2[c]);
	  if(overlapstart < overlapend - 1e-6 && INDEL){
	    /* flag the lower scoring alignment (q) query region that overlaps the higher scoring alignment's query region so indels are NOT output from that region */
	    // HERE HERE : update to handle multiple regions as in 1-color code
	    if(q[1+c].start < q[1+c].end){/* combine new region with old region */
	      q[1+c].start = min(q[1+c].start,overlapstart);
	      q[1+c].end = max(q[1+c].end,overlapend);
	    } else {
	      q[1+c].start = overlapstart;
	      q[1+c].end = overlapend;
	    }
	  }
	  if(VERB>=2 && overlapstart < overlapend){
	    printf("c=%d:Overlap Region in QueryId=%lld(len=%0.3f) with RefId=%lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) and %lld(N=%d,len=%0.3f:R=%d..%d:%0.3f..%0.3f;Q=%d..%d,or=%d:%0.3f..%0.3f) (Log10PV=%0.2f and %0.2f):overlap=%d..%d(%0.3f..%0.3f):alignment[%lu]:start=%d,end=%d\n",
		   c,Xmap->id,origXX[c][MM[c]+1]*Xscale,YYmap[rid]->id,NN[c],origYY[c][NN[c]+1]*Yscale,LI[c],RI[c],refstartpos,refendpos,qrystart[c],qryend[c],p->orientation,origXX[c][qrystart[c]],origXX[c][qryend[c]],
		   YYmap[q->mapid1]->id,NN2[c],origYY2[c][NN2[c]+1]*Yscale,LI2[c],RI2[c],refstartpos2,refendpos2,qrystart2[c],qryend2[c],q->orientation,origXX[c][qrystart2[c]],origXX[c][qryend2[c]],p->logPV,q->logPV,
		   overlapstart,overlapend,origXX[c][overlapstart],origXX[c][overlapend],k,q[1+c].start,q[1+c].end);
	    fflush(stdout);
	  }
	} /* c = 0..colors-1 */
      } /* k = i+1 .. */
    } /* i = start .. end -1 */
    if(VERB){
      printf("Deleted %d Overlapped alignments\n",delcnt);
      if(VERB/* HERE >=2 */){
	printf("Alignments with partial overlaps:\n");
	for(size_t i = start; i < end; i++){
	  Calign *p = Lalignment[i];
	  if(p){
	    for(int c = 0; c < colors; c++){
	      if(p[1+c].start < p[1+c].end){
		Cmap *qmap = XXmap[p->mapid2];
		printf("Lalignment[%lu]:rid=%lld,qid=%lld,or=%d,M=%d,c=%d: suppressed interval %d .. %d (%0.3f .. %0.3f)\n",
		       i,YYmap[p->mapid1]->id, qmap->id,p->orientation,qmap->numsite[c],c,p[1+c].start,p[1+c].end,qmap->site[c][p[1+c].start],qmap->site[c][p[1+c].end]);
	      }
	    }
	  }		   
	}
      }
      fflush(stdout);
    }
  }

  /* sort alignments in order of id1, sites1[0],sites1[N+1], id2, sites2[0],sites2[M+1] (null alignments last) */
  qsort(&Lalignment[start], end - start, sizeof(Calign *), (intcmp*)SiteInc2);

  /* sort alignments in order of LogPV, id2, id1, sites1[0],sites1[numpairs-1] (null alignments last) */
  if(0)// for debugging
    qsort(&Lalignment[start],end - start,sizeof(Calign *), (intcmp*)LogPVSiteInc2);

  fprintf(fp,"# XMAP File Version:\t3.0\n");
  fprintf(fp,"# Label Channels:\t%d\n", colors);
  fprintf(fp,"# XmapEntryID's for Nth matchgroup are N on two lines with LabelChannels 1 and 2 respectively\n");

  if(spots_filename[0] || (svcheck && numrefmaps > 0)){/* Reference Aligner : refer to condensed version of reference CMAP */
    fprintf(fp,"# Reference Maps From:\t%s_r.cmap\n", absbasename);
    fprintf(fp,"# Query Maps From:\t%s", vfixx_filename[0]);/* This may actually be a modified cmap file */
    if(DEBUG) assert(num_files==1);
  }  else {
    if(mres > 0.001 || mresSD > 0.0 || fabs(PixelLen - origPixelLen) > 1e-12){
      fprintf(fp,"# Reference Maps From:\t%s_r.cmap\n", absbasename);
      fprintf(fp,"# Query Maps From:\t%s_q.cmap", absbasename);
    } else {
      fprintf(fp,"# Reference Maps From:\t%s", vfixx_filename[0]);
      for(register int i = 1; i < RefIndex; i++)
#ifndef WIN32
	if(vfixx_filename[i][0] == '/')
	  fprintf(fp,",%s",vfixx_filename[i]);
	else
	  fprintf(fp,",%s/%s",cwd,vfixx_filename[i]);
#else
	  fprintf(fp,",%s",vfixx_filename[i]);
#endif
      fprintf(fp,"\n");
#ifndef WIN32
      if(vfixx_filename[QueryIndex][0] == '/')
	fprintf(fp,"# Query Maps From:\t%s",vfixx_filename[QueryIndex]);
      else
	fprintf(fp,"# Query Maps From:\t%s/%s",cwd,vfixx_filename[QueryIndex]);
#else
	fprintf(fp,"# Query Maps From:\t%s",vfixx_filename[QueryIndex]);
#endif
      for(register int i = QueryIndex+1; i < num_files;i++)
#ifndef WIN32
	if(vfixx_filename[i][0] == '/')
	  fprintf(fp,",%s",vfixx_filename[i]);
	else
	  fprintf(fp,",%s/%s",cwd,vfixx_filename[i]);
#else
	  fprintf(fp,",%s",vfixx_filename[i]);
#endif
    }
  }
  fprintf(fp,"\n");

  if((PairSplit || (MultiMatches && RefSplit)) && dosvdetect){ /* see parameters.{h,cpp} */
    fprintf(fp,"# Smap Entries In:\t%s.smap\n",basename);
  }

  if(SecondBest && SinglelineXmap){
    fprintf(fp,"#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment\tQryEndPos2\tRefEndPos2\tOrientation2\tConfidence2\tNumMatch2\n");
    fprintf(fp,"#f int        \tint        \tint        \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat     \tstring \tflat  \tfloat \tint         \tstring   \tfloat     \tfloat     \tstring      \tfloat      \tint      \n");
  } else {
    fprintf(fp,"#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment\n");
    fprintf(fp,"#f int        \tint        \tint        \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat     \tstring \tfloat \tfloat \tint         \tstring   \n");
  }


  if((PairSplit || (MultiMatches && RefSplit)) && (dosvdetect || doscaffold)){ /* see parameters.{h,cpp} */
    if(DEBUG) assert(end - start <= MASK(31));
    delete [] xmapentries;
    maxxmapentries = end - start;
    xmapentries = new xmapEntry*[maxxmapentries]; /* allocate for numaligns entries, though actual number will be <= this */
    numxmapentries = 0;
  }

  /* allocate interval array to keep track of total and unique alignment lengths */
  Cinterval *Rinterval = new Cinterval[(end - start + 1)*2];
  Cinterval *Qinterval = &Rinterval[end - start + 1];
  int RintervalCnt = 0, QintervalCnt = 0;

  size_t aligncnt = 0;
  size_t SmapCnt = 0;

  double Ilog10 = 1.0/log(10.0);
  double EndOutlier = log(max(1.0e-300,PoutlierEnd));
  double Outlier = log(max(1.0e-300,Poutlier));

  for(size_t i= start; i < end; i++){
    Calign *p = Lalignment[i];
    if(!p)
      continue;
    if(VERB>=2){
      printf("i=%lu,aligncnt=%lu:Lalignment[i]=%p,Lalignment[i]->align2=%p,chimpair=%d:numpairs=%d,score=%0.3f,logPV=%0.2f(at loop start)\n",
	     i,aligncnt,p,p->align2,p->chimpair,p->numpairs,p->score,p->logPV);
      fflush(stdout);
    }

    int U = p->numpairs;

    int SecondLine = 0;
    if(SecondBest && !SinglelineXmap && p->chimpair < 0){
      if(DEBUG/* HERE >=2 */) assert(!(U<=0 || p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold));
      SecondLine = 1;
    }

    /* Below p->chimpair is changed to the XmapEntryID in the output file (0 means no output) */

    if(!SecondLine){
      if(U<=0 || p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold){
	if(VERB>=2){
	  printf("Lalignment[%lu]= %p: numpairs=%d,score=%0.3f,logPV=%0.2f, align2= %p, chimpair=%d: skipping\n",i, p, U, p->score, p->logPV, p->align2, p->chimpair);
	  fflush(stdout);
	}
	p->chimpair = 0;
	continue;
      }
    } else {
      if(DEBUG) assert(U > 0);
    }

    int rid = p->mapid1;
    if(DEBUG) assert(0 <= rid && rid < numY);
    int mid = p->mapid2;
    if(DEBUG) assert(0 <= mid && mid < numX);
    //    Cmap *Ymap = YYmap[rid];
    Cmap *Xmap = XXmap[mid];

    Cmap *origmap = Xmap;
    while(origmap->origmap){
      Xmap->origmap = origmap = origmap->origmap;
      if(origmap->origmap){/* otherwise origmap->left[] is not defined */
	for(int c = 0; c < colors; c++){
	  Xmap->left[c] += max(0,origmap->left[c]-1);
	  Xmap->right[c] += max(0,origmap->left[c]-1);
	}
      }
    }
    if(BestRef && origmap->align->mapid1 != rid){
      p->chimpair = 0;
      continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */
    }

    aligncnt++;
    p->chimpair = (p->chimpair < 0) ? -aligncnt : aligncnt;
  }

  extern double raYlambda;// see RGentigRefScore.h
  double AvgInterval = raYlambda;

  int numindel = 0, maxindel = 1024;
  Cindel *indel = new Cindel[maxindel];

  /* now all |p->chimpair| are the planned XmapEntryID in output file (or 0 if alignment is not output) */

  for(size_t i = start; i < end; i++){
    Calign *p = Lalignment[i];
    if(!p || p->chimpair==0)
      continue;

    int U = p->numpairs;

    int SecondLine = 0;
    if(SecondBest && !SinglelineXmap && p->chimpair < 0)
      SecondLine = 1;

    aligncnt = abs(p->chimpair);

    int rid = p->mapid1;
    if(DEBUG) assert(0 <= rid && rid < numY);
    int mid = p->mapid2;
    if(DEBUG) assert(0 <= mid && mid < numX);
    Cmap *Ymap = YYmap[rid];
    Cmap *Xmap = XXmap[mid];

    Cmap *origmap = Xmap;
    while(origmap->origmap){
      Xmap->origmap = origmap = origmap->origmap;
      if(origmap->origmap){/* otherwise origmap->left[] is not defined */
	Xmap->left[0] += max(0,origmap->left[0]-1);
	Xmap->right[0] += max(0,origmap->left[0]-1);
      }
    }

    if(VERB>=2){
      printf("aligncnt=%lu: Xmap Lalignment[%lu]:mid=%lld,rid=%lld:numpairs=%d,score=%0.6f,logPV=%0.2f,chimpair=%d,start=%d,end=%d\n",aligncnt,i,Xmap->id,Ymap->id,U,p->score,p->logPV,p->chimpair,p->start,p->end);
      fflush(stdout);
    }

    int *MM = Xmap->numsite;
    int *NN = Ymap->numsite;
    FLOAT **YY = Ymap->site;
    FLOAT **XX = Xmap->site;
    FLOAT **origYY = Ymap->origmap ? Ymap->origmap->site : YY;
    FLOAT **origXX = Xmap->origmap ? Xmap->origmap->site : XX;
    int *origNN = Ymap->origmap ? Ymap->origmap->numsite : NN;
    int *origMM = Xmap->origmap ? Xmap->origmap->numsite : MM;

    if(DEBUG) assert(Xmap->origmap==0 || Xmap->origmap->origmap==0);/* If this fails fix origmap link in pairalign() or refalign() to point to root map */
    if(DEBUG) assert(Ymap->origmap==0 || Ymap->origmap->origmap==0);/* If this fails fix origmap link (same as for Xmap above) to point to root map */

    Cinterleaved *origYsite = NULL, *origXsite = NULL;
    int *Yindex[2] = {NULL,NULL}, Ysites=0;/* Yindex[c][I=1..origNN[c]] is the interleaved index into origYsite[] corresponding to origYY[c][I] */
    int *Xindex[2] = {NULL,NULL}, Xsites=0;/* Xindex[c][I=1..origMM[c]] is the interleaved index into origXsite[] corresponding to origXX[c][I] */

    if(XmapInterleaved){  // convert each 2-color map into an interleaved format : array of (c,index), corresponding to origYY[c][index] and in ascending order of this value
      origYsite = new Cinterleaved[origNN[0]+origNN[1]+2];
      Yindex[0] = new int[origNN[0]+origNN[1]+2];
      Yindex[1] = &Yindex[0][origNN[0]+1];

      origYsite[0].c = 0;
      origYsite[0].index = 0;
      Ysites = 1;
      int I[2] = {1,1};
      for(;I[0] <= origNN[0] || I[1] <= origNN[1];){
	if(I[0] > origNN[0] || (I[1] <= origNN[1] && origYY[0][I[0]] > origYY[1][I[1]])){/* next site is origYY[1][I[1]] */
	  Yindex[1][I[1]] = Ysites;
	  origYsite[Ysites].c = 1;
	  origYsite[Ysites++].index = I[1]++;
	} else {/* next site is origYY[0][I[0]] */
	  if(DEBUG) assert(I[1] > origNN[1] || (I[0] <= origNN[0] && origYY[1][I[1]] >= origYY[0][I[0]]));
	  Yindex[0][I[0]] = Ysites;
	  origYsite[Ysites].c = 0;
	  origYsite[Ysites++].index = I[0]++;
	}
      }
      origYsite[Ysites].c = 0;
      origYsite[Ysites--].index = origNN[0]+1;
      if(DEBUG) assert(Ysites == origNN[0]+origNN[1]);
      
      origXsite = new Cinterleaved[origMM[0]+origMM[1]+2];
      Xindex[0] = new int[origMM[0]+origMM[1]+2];
      Xindex[1] = &Xindex[0][origMM[0]+1];

      origXsite[0].c = 0;
      origXsite[0].index = 0;
      Xsites = 1;
      I[0] = I[1] = 1;
      for(;I[0] <= origMM[0] || I[1] <= origMM[1];){
	if(I[0] > origMM[0] || (I[1] <= origMM[1] && origXX[0][I[0]] > origXX[1][I[1]])){/* next site is origXX[1][I[1]] */
	  Xindex[1][I[1]] = Xsites;
	  origXsite[Xsites].c = 1;
	  origXsite[Xsites++].index = I[1]++;
	} else {/* next site is origXX[0][I[0]] */
	  if(DEBUG) assert(I[1] > origMM[1] || (I[0] <= origMM[0] && origXX[1][I[1]] >= origXX[0][I[0]]));
	  Xindex[0][I[0]] = Xsites;
	  origXsite[Xsites].c = 0;
	  origXsite[Xsites++].index = I[0]++;
	}
      }
      origXsite[Xsites].c = 0;
      origXsite[Xsites--].index = origMM[0]+1;
      if(DEBUG) assert(Xsites == origMM[0]+origMM[1]);
    }

    /*    if(Ymap->origmap || Xmap->origmap){
      printf("output_xmap2: split maps not handled with 2 colors (try -nosplit 2)\n");
      exit(1);
      }*/

    /* add interval information to Rinterval[],Qinterval[] */
    Cinterval *ri = &Rinterval[RintervalCnt++];
    Cinterval *qi = &Qinterval[QintervalCnt++];
    ri->id = YYmap[rid]->id;
    ri->left = origYY[0][origNN[0]+1];
    ri->right = 0.0;
    qi->id = XXmap[mid]->id;
    qi->left = origXX[0][origMM[0]+1];
    qi->right = 0.0;

    /* NOTE : indels will be accumulated for both colors at indel[0 .. numindel-1] then combined if their intervals intersect */
    numindel = 0;

    for(int c = 0; c < colors; c++){/* output one line per Label Channel */
      Calign *q = &p[c+1];

      int M = MM[c];
      int N = NN[c];
      //      FLOAT *Y = YY[c];
      FLOAT *X = XX[c];
      FLOAT *origX = origXX[c];
      FLOAT *origY = origYY[c];
      int Yshift = Ymap->origmap ? max(1,Ymap->left[c]) - 1 : 0;
      int Xshift = Xmap->origmap ? max(1,Xmap->left[c]) - 1 : 0;
      //      int origM = origMM[c];
      //      int origN = origNN[c];

      int U = q->numpairs;
      int LI = q->sites1[0];
      int LJ = q->sites2[0];
      int RI = q->sites1[U-1];
      int RJ = q->sites2[U-1];

      // Get the indices, subsequently use to compute number of un-aligned labels 
      int qrystartidx = p->orientation ? Xshift + M + 1 - LJ : Xshift + LJ;
      int qryendidx = p->orientation ? Xshift + M + 1 - RJ : Xshift + RJ;
      double qrystartpos = origX[qrystartidx]*1000.0*Xscale; /* QryStartPos */
      double qryendpos   = origX[qryendidx  ]*1000.0*Xscale; /* QryEndPos */
      double qryLen = origX[M+1]*1000.0*Xscale; /* QryLen */
      double refstartpos = origY[Yshift + LI]*1000.0*Yscale; /* RefStartPos */
      double refendpos   = origY[Yshift + RI]*1000.0*Yscale; /* RefEndPos */
      double refLen = origY[N+1]*1000.0*Yscale;/* RefLen */

      /* update interval information to Rinterval[],Qinterval[] */
      ri->left = min(ri->left,refstartpos);
      ri->right = max(ri->right,refendpos);
      qi->left = min(qi->left,min(qrystartpos,qryendpos));
      qi->right = max(qi->right,max(qrystartpos,qryendpos));

      fprintf(fp,"%lu\t%lld\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%c\t%0.2f\t",
	    aligncnt, /* XmapEntryID */
	    Xmap->id, /* QryContigID */
	    Ymap->id, /* RefContigID */
	    qrystartpos,
	    qryendpos,
	    refstartpos,
	    refendpos,
	    p->orientation ? '-' : '+', /* Orientation */
	    p->logPV/*zscore(pow(10.0,-max(0.0,p->logPV)))*/ /* Confidence */
	    );

      if(INDEL){/* append indels for current color to indel[0..numindel-1] */
	int lastI = q->sites1[0];
	int lastK = (resSD[0] > 0.0) ? q->sitesK1[0] : 0;
	int lastJ = q->sites2[0];
	int I,K,J;

	for(int k = 1; k < U; lastI = I, lastK = K, lastJ = J, k++){
	  I = q->sites1[k];
	  K = (resSD[0] > 0.0) ? q->sitesK1[k] : 0;
	  J = q->sites2[k];

	  int outlier = (q->iscore[k] > q->outscore[k] + (RFLOAT)0.01 || (outlierExtend && (lastJ >= J || lastI >= I-K))) ? 1 : 0;
	  if(!outlier /* q->outscore[k] >= q->iscore[k] - Poutlier - (RFLOAT)0.01*/)/* not an outlier */
	    continue;

	  /*  merge nearby outliers into a single SV ? This is currently done by -outlierExtend (and subsequently by -svcheck) */

	  /* outlier lies between sites lastI .. I-K on reference and lastJ .. J on query */
	  int qrystartidx = q->orientation ? Xshift + M + 1 - lastJ : Xshift + lastJ;
	  int qryendidx   = q->orientation ? Xshift + M + 1 - J : Xshift + J;
	  double qrystartpos = origX[qrystartidx]*1000.0*Xscale; /* QryStartPos */
	  double qryendpos   = origX[qryendidx  ]*1000.0*Xscale; /* QryEndPos */
	  double refstartpos = origY[Yshift+lastI-lastK]*1000.0*Yscale; /* RefStartPos */
	  double refendpos   = origY[Yshift+I]*1000.0*Yscale; /* RefEndPos */

	  /* include half of de-res interval in query interval, so query interval can be replaced with ref interval without introducing a length bias */
	  double origqsiz = p->orientation ? (qrystartpos - qryendpos) : (qryendpos - qrystartpos);
	  if(DEBUG/* HERE >=2 */) assert(origqsiz >= 0.0);
	  if(!p->orientation){
	    if(lastK > 0)
	      qrystartpos -= 0.5*(origY[Yshift+lastI] - origY[Yshift+lastI-lastK]) * 1000.0 * Yscale;
	    if(K > 0)
	      qryendpos += 0.5*(origY[Yshift+I] - origY[Yshift+I-K]) * 1000.0 * Yscale;
	  } else {/* query runs from high to low positions along alignment */ 
	    if(lastK > 0)
	      qrystartpos += 0.5*(origY[Yshift+lastI] - origY[Yshift+lastI-lastK]) * 1000.0 * Yscale;
	    if(K > 0)
	      qryendpos -= 0.5*(origY[Yshift+I] - origY[Yshift+I-K]) * 1000.0 * Yscale;
	  }
	  if(DEBUG/* HERE >=2 */) assert(refstartpos <= refendpos);
	  double newqsiz =  p->orientation ? (qrystartpos - qryendpos) : (qryendpos - qrystartpos);
	  if(DEBUG/* HERE >=2 */ && !(newqsiz >= origqsiz)){
	    printf("Yid=%lld,Xid=%lld,or=%d,c=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f:query=%0.3f..%0.3f,ref=%0.3f..%0.3f,qryindex=%d..%d(M=%d)\n",
		   Ymap->id,Xmap->id,p->orientation,c,I,K,J,lastI,lastK,lastJ,q->outscore[k],q->iscore[k],qrystartpos,qryendpos,refstartpos,refendpos,qrystartidx,qryendidx,M);
	    printf("Xmap->origmap=%p,Xmap=%p,Xshift=%d,origX=%p,X=%p,Xscale=%0.8f\n",Xmap->origmap,Xmap,Xshift,origX,X,Xscale);
	    for(int I = 0; I <= M+1; I++)
	      printf("X[%d]=%0.6f %c\n",I,X[I], (I > 0 && X[I] < X[I-1]) ? '!' : ' ');
	    fflush(stdout);
	    assert(newqsiz >= origqsiz);
	  }
	  if(DEBUG/* HERE >=2 */) assert(J > lastJ);
	  if(DEBUG/* HERE >=2 */) assert(I-K > lastI);

	  if(IndelChiSq > 0.0 && CmapChimQuality >= 3 && spots_filename[0] && qryendpos - qrystartpos < refendpos - refstartpos){// ignore outlier based in -indelChiSq
	    double Outlier = fabs(refendpos - refstartpos - qryendpos + qrystartpos);// indel size

	    float *OutlierFrac = Xmap->OutlierFrac[c];
	    float *FragSd = Xmap->FragSd[c];
	    float *ExpSd = Xmap->ExpSd[c];
	    float *FragCov = Xmap->FragCov[c];
	    float *FragChiSq = Xmap->FragChiSq[c];
	    if(DEBUG/* HERE HERE >=2 */) assert(OutlierFrac != NULL && FragSd != NULL && ExpSd != NULL && FragCov != NULL && FragChiSq != NULL);
	    double LRange = (origY[Yshift + I] - origY[Yshift + lastI]) * Yscale;
	    double LRangeRatio = LRange / AvgInterval;

	    int jmin = qrystartidx;
	    int jmax = qryendidx;
	    int j = jmin;
	    for(; j < jmax; j++){
	      double SRange = (origX[j+1] - origX[j]) * Xscale;
	      if(Outlier < IndelMaxSize &&
		 (OutlierFrac[j] > IndelFrac || 
		  (FragCov[j] < IndelMaxCov && 
		   LRange > IndelMinInterval &&
		   ((FragChiSq[j] < IndelChiSq &&
		     FragSd[j] > ExpSd[j] + IndelMinSf + IndelMinSr * SRange &&
		     FragSd[j] > ExpSd[j] * IndelSdRatio) ||		     
		    FragCov[j] < IndelMinCov + LRangeRatio *  IndelMinCovScale)))){
		if(VERB>=2){
		  printf("\t j=%d(%d..%d):J=%d..%d:LRange= %0.4f(Ratio= %0.3f), SRange= %0.4f, OutlierFrac= %0.1f, FragChiSq= %0.4e, Outlier= %0.4f, FragCov= %0.2f, FragSd= %0.1f, ExpSd= %0.1f : skipping outlier\n",
			 j,jmin,jmax,lastJ,J,LRange,LRangeRatio,SRange, OutlierFrac[j], FragChiSq[j], Outlier, FragCov[j], FragSd[j],ExpSd[j]);
		  printf("\t    MaxSize= %0.1f, MaxOutlierFrac= %0.1f, MaxCov= %0.2f, ChiSq= %0.4e, MinSf= %0.4f, MinSr= %0.4f, SdRatio= %0.4f, MinInterval= %0.4f, MinCov= %0.2f, MinCovR= %0.2f\n",
			 IndelMaxSize, IndelFrac,IndelMaxCov,IndelChiSq,IndelMinSf,IndelMinSr,IndelSdRatio,IndelMinInterval, IndelMinCov, IndelMinCovScale);
		  fflush(stdout);
		}

		break;
	      }
	    }

	    if(j < jmax)
	      continue;// outlier will be ignored
	  }


	  /* compute confidence in outlier based on lesser of:
	     1. Confidence in entire alignment
	     2. Confidence that this region is an outlier (vs a regular alignment)
	     3. Confidence that this region is an outlier (vs an EndOutlier in either direction).
	     NOTE : should really be adding the Pvalues, but this is approximately correct and quicker
	  */
	  double IndelConfidence = min((q->iscore[k] - q->outscore[k]-Outlier)*Ilog10, p->logPV);
	  /* compare with end outlier to the left */
	  double Leftscore = 0.0;
	  for(int t = 0; t <= k; t++)
	    Leftscore += q->iscore[t];

	  if((Leftscore - EndOutlier)*Ilog10 < IndelConfidence)
	    IndelConfidence = (Leftscore - EndOutlier)*Ilog10;
	  /* compare with end outlier to the right */
	  double Rightscore = 0.0;
	  for(int t = k; t <= U; t++)
	    Rightscore += q->iscore[t];
	  if((Rightscore - EndOutlier)*Ilog10 < IndelConfidence)
	    IndelConfidence = (Rightscore - EndOutlier)*Ilog10;

	  if(VERB>=2 && Ymap->id==1 && Xmap->id==318){
	    printf("Yid=%lld,Xid=%lld,or=%d,c=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f,Outlier=%0.6f:LeftScore=%0.6f,RightScore=%0.6f,EndOutlier=%0.6f,Confidence=%0.6f\n",
		   Ymap->id,Xmap->id,p->orientation,c,I,K,J,lastI,lastK,lastJ,q->outscore[k],q->iscore[k],Outlier,Leftscore,Rightscore,EndOutlier,IndelConfidence);
	    fflush(stdout);
	  }

	  if(IndelConfidence < -Outlier * Ilog10)
	    continue;

	  /* make sure there is no other alignment of this Query (Xmap) overlapping this indel location for which the
	     basic alignment has a better logPV : If they exist, the overlap region is specified in the alignment as sites q->start .. q->end */
	  if(DEBUG) assert(q->start <= q->end);
	  if(p->orientation){
	    if(DEBUG) assert(qryendidx <= qrystartidx);
	    if(max(qryendidx,q->start) < min(qrystartidx,q->end)){
	      if(VERB>=2 && Ymap->id==23 && Xmap->id==60){
		printf("Yid=%lld,Xid=%lld,or=%d,c=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f:query=%0.3f..%0.3f,ref=%0.3f..%0.3f:outlier skipped due to higher scoring alignment overlap\n",
		       Ymap->id,Xmap->id,p->orientation,c,I,K,J,lastI,lastK,lastJ,q->outscore[k],q->iscore[k],qrystartpos,qryendpos,refstartpos,refendpos);
		fflush(stdout);
	      }
	      continue;
	    }
	  } else {
	    if(DEBUG) assert(qryendidx >= qrystartidx);
	    if(max(qrystartidx,q->start) < min(qryendidx,q->end)){
	      if(VERB>=2 && Ymap->id==23 && Xmap->id==60){
		printf("Yid=%lld,Xid=%lld,or=%d,c=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f:query=%0.3f..%0.3f,ref=%0.3f..%0.3f:outlier skipped due to higher scoring alignment overlap\n",
		       Ymap->id,Xmap->id,p->orientation,c,I,K,J,lastI,lastK,lastJ,q->outscore[k],q->iscore[k],qrystartpos,qryendpos,refstartpos,refendpos);
		fflush(stdout);
	      }
	      continue;
	    }
	  }

	  if(VERB>=2 && Ymap->id==23 && Xmap->id==60){
	    printf("Yid=%lld,Xid=%lld,or=%d,c=%d:I=%d,K=%d,J=%d,lastI=%d,lastK=%d,lastJ=%d:outscore=%0.6f,iscore=%0.6f:query=%0.3f..%0.3f,ref=%0.3f..%0.3f:indel found(Xshift=%d,Yshift=%d)\n",
		   Ymap->id,Xmap->id,p->orientation,c,I,K,J,lastI,lastK,lastJ,q->outscore[k],q->iscore[k],qrystartpos,qryendpos,refstartpos,refendpos,Xshift,Yshift);
	    fflush(stdout);
	  }

	  if(numindel >= maxindel){
	    maxindel = max(1024,2*maxindel);
	    Cindel *nindel = new Cindel[maxindel];
	    memcpy(nindel,indel,numindel*sizeof(Cindel));
	    delete [] indel;
	    indel = nindel;
	  }
	  Cindel *pindel = &indel[numindel++];
	  pindel->Xmap_id = Xmap->id;
	  pindel->Ymap_id = Ymap->id;
	  pindel->qrystartpos = qrystartpos;
	  pindel->qryendpos = qryendpos;
	  pindel->refstartpos = refstartpos;
	  pindel->refendpos = refendpos;
	  if(DEBUG) assert(pindel->refstartpos < pindel->refendpos);
	  pindel->orientation = p->orientation;
	  pindel->confidence = IndelConfidence;
	  pindel->c = c;
	  pindel->aligncnt = aligncnt;

#if 0
	  SmapCnt++;
	  if(fpSV)
	    fprintf(fpSV,"%llu\t%lld\t%lld\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%c\t%0.2f\t%s%dI%dD%d\t%lu\t%lu\n",
		  (unsigned long long)SmapCnt,    // SmapEntryID
		  Xmap->id,   // QryContigID
		  Ymap->id,   // RefContigID1
		  Ymap->id,   // RefContigID2
		  qrystartpos,
		  qryendpos,
		  refstartpos,
		  refendpos,
		  p->orientation ? '-' : '+', // Orientation
		  IndelConfidence, 
		  (refendpos-refstartpos > fabs(qryendpos-qrystartpos)) ? "deletion_" : "insertion_", // Type
		  (int)fabs(fabs(qryendpos-qrystartpos) - (refendpos-refstartpos)),// indel size in bp
		  J - lastJ - 1 + (lastK > 0 ? 1 : 0) + (K > 0 ? 1 : 0), // inserted sites
		  I - (lastI-lastK) - 1, // deleted sites
		  aligncnt, aligncnt); // XmapID1, XmapID2
#endif
	} /* lastJ, J loop */
      }

      if((dosvdetect || doscaffold) && !SecondLine){
	printf("svdetect or doscaffold not implemented for 2 colors\n");
	exit(1);
      }

      /* Output Alignment Cigar */
      int lastI = q->sites1[0];
      int lastJ = q->sites2[0];
      int matchcnt = 1;

      for(int k = 1; k < U; k++){
	int I = q->sites1[k];
	int K = (resSD[0] > 0.0) ? q->sitesK1[k] : 0;
	int J = q->sites2[k];

	if(DEBUG) assert(I-K > lastI && J > lastJ);

	if(K==0 && I == lastI + 1 && J == lastJ + 1){
	  matchcnt++;
	} else {
	  /* output match so far */
	  if(matchcnt > 0)
	    fprintf(fp,"%dM",matchcnt);
	  if(I == lastI + 1)
	    fprintf(fp,"%dI",J - lastJ -1);
	  else if(J == lastJ + 1)
	    fprintf(fp,"%dD", I - lastI -1);
	  else /* both an insertion and a deletion */
	    fprintf(fp,"%dI%dD", J - lastJ - 1, I - lastI - 1);
	  matchcnt = 1;
	}
	lastI = I;
	lastJ = J;
      }
      if(matchcnt > 0)
	fprintf(fp,"%dM", matchcnt);
      /* done with cigar */

      fprintf(fp,"\t%0.1f\t%0.1f\t%d\t",qryLen,refLen,c+1);
  
      /* Output Alignment doublet string */
      if(!XmapInterleaved){
	for(int k = 0; k < U; k++){
	  int I = q->sites1[k];
	  int K = (resSD[0] > 0.0) ? q->sitesK1[k] : 0;
	  int J = q->sites2[k];
	  
	  int origI = Yshift + I;
	  int origJ = p->orientation ? Xshift + M + 1 - J : Xshift + J;
	  if(K > 0 )
	    fprintf(fp,"(%d,%d)",origI-K,origJ);
	  fprintf(fp,"(%d,%d)",origI,origJ);
	}
      } else {
	for(int k = 0; k < U; k++){
	  int I = q->sites1[k];
	  int K = (resSD[0] > 0.0) ? q->sitesK1[k] : 0;
	  int J = q->sites2[k];
	  
	  int origI = Yshift + I;
	  int origK = origI-K;
	  int origJ = p->orientation ? Xshift + M + 1 - J : Xshift + J;
	  
	  /* convert per-color index values (origI,origK,origJ) into interleaved index values */
	  origI = Yindex[c][origI];
	  origK = Yindex[c][origK];
	  origJ = Xindex[c][origJ];

	  if(K > 0 )
	    fprintf(fp,"(%d,%d)",origK,origJ);
	  fprintf(fp,"(%d,%d)",origI,origJ);
	}
      }

      if(SecondBest && SinglelineXmap){
	Calign *p2 = p->align2;
	if(p2){
	  if(DEBUG) assert(!p2->align2);
	  Calign *q2 = &p2[c+1];
	  int U2 = q2->numpairs;
	  int RI2 = q2->sites1[U2-1];
	  //      register int LK2 = q2->sitesK1[U2-1];
	  int RJ2 = q2->sites2[U2-1];
	  fprintf(fp,"\t%0.1f\t%0.1f\t%c\t%0.2f\t%d",
		  (p2->orientation ? origX[Xshift + M + 1 - RJ2] : origX[Xshift + RJ2])*1000.0*Xscale,/* QryEndPos2 */
		  origY[Yshift + RI2]*1000.0*Yscale,/* RefEndPos2 */
		  p2->orientation ? '-' : '+', /* Orientation2 */
		  p2->logPV, /* Confidence2 */
		  q2->numpairs /* NumMatch2 */
		  );
	} else
	  fprintf(fp,"\t%0.1f\t%0.1f\t%c\t%0.2f\t%d",
		  0.0,/* QryEndPos2 */
		  0.0,/* RefEndPos2 */
		  '+', /* Orientation2 */
		  0.0, /* Confidence2 */
		  0 /* NumMatch2 */
		  );
      }
      
      fprintf(fp,"\n");
    }/* c = 0 .. colors-1 */

    if(XmapInterleaved){
      delete [] origYsite;
      delete [] origXsite;
      delete [] Yindex[0];
      delete [] Xindex[0];
      origYsite = origXsite = NULL;
      Yindex[0] = Yindex[1] = Xindex[0] = Xindex[1] = NULL;
    }

    if(INDEL){/* combine indels of opposite color in indel[0..numindel-1] that overlap into a single indel, then output all remaining indels */
      for(int i = 0; i < numindel-1; i++){
	Cindel *pindel = &indel[i];
	for(int j = i+1; j < numindel; j++){
	  Cindel *qindel = &indel[j];
	  if(pindel->c != qindel->c){
	    double qstart = pindel->orientation ? max(pindel->qryendpos, qindel->qryendpos) : max(pindel->qrystartpos, qindel->qrystartpos);
	    double qend = pindel->orientation ? min(pindel->qrystartpos, qindel->qrystartpos) : min(pindel->qryendpos, qindel->qryendpos);
	    double rstart = max(pindel->refstartpos,qindel->refstartpos);
	    double rend = min(pindel->refendpos,qindel->refendpos);
	    double qoverlap = qend - qstart;
	    double roverlap = rend - rstart;
	    if(qoverlap > 0.0 && roverlap > 0.0){/* merge the two indels into a single smaller one */
	      pindel->qrystartpos = pindel->orientation ? qend : qstart;
	      pindel->qryendpos = pindel->orientation ? qstart : qend;
	      pindel->refstartpos = rstart;
	      pindel->refendpos = rend;
	      if(DEBUG) assert(pindel->refstartpos < pindel->refendpos);
	      pindel->confidence = max(pindel->confidence, qindel->confidence);

	      for(int k = j+1; k < numindel; k++)
		indel[k-1] = indel[k];
	      numindel--;

	      break;
	    }
	  }
	}
      }
      
      qsort(indel, numindel, sizeof(Cindel), (intcmp*) IndelRefInc);

      for(int i = 0; i < numindel; i++){
	Cindel *pindel = &indel[i];
	SmapCnt++;
	
	/* compute number of inserted and deleted sites */
	int insertcnt = 0, deletecnt = 0;
	if(pindel->orientation){
	  if(DEBUG) assert(pindel->qrystartpos > pindel->qryendpos);
	  for(int c = 0; c < colors; c++){
	    int J = 1;
	    for(; J <= origMM[c]; J++)
	      if(origXX[c][J] >= pindel->qryendpos)
		break;
	    for(;J <= origMM[c]+1; J++){
	      if(origXX[c][J] > pindel->qrystartpos)
		break;
	      insertcnt++;
	    }
	  }
	} else {
	  if(DEBUG) assert(pindel->qrystartpos < pindel->qryendpos);
	  for(int c = 0; c < colors; c++){
	    int J = 1;
	    for(; J <= origMM[c]; J++)
	      if(origXX[c][J] >= pindel->qrystartpos)
		break;
	    for(;J <= origMM[c]+1; J++){
	      if(origXX[c][J] > pindel->qryendpos)
		break;
	      insertcnt++;
	    }
	  }
	}
	if(DEBUG && !(pindel->refstartpos < pindel->refendpos)){
	  printf("Xid=%lld,Yid=%lld,indel[%d]:refstartpos= %0.4f, refendpos= %0.4f\n",Xmap->id,Ymap->id,i,pindel->refstartpos,pindel->refendpos);
	  fflush(stdout);
	  assert(pindel->refstartpos < pindel->refendpos);
	}
	for(int c = 0; c < colors; c++){
	  int I = 1;
	  for(; I <= origNN[c]; I++)
	    if(origYY[c][I] >= pindel->refstartpos)
	      break;
	  for(;I <= origNN[c]+1; I++){
	    if(origYY[c][I] > pindel->refendpos)
	      break;
	    deletecnt++;
	  }
	}

	if(fpSV)
	  fprintf(fpSV,"%llu\t%lld\t%lld\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%c\t%0.2f\t%s%dI%dD%d\t%lu\t%lu\n",
		(unsigned long long)SmapCnt,    // SmapEntryID
		Xmap->id,   // QryContigID
		Ymap->id,   // RefContigID1
		Ymap->id,   // RefContigID2
		pindel->qrystartpos,
		pindel->qryendpos,
		pindel->refstartpos,
		pindel->refendpos,
		p->orientation ? '-' : '+', // Orientation
		pindel->confidence,
		(pindel->refendpos - pindel->refstartpos > fabs(pindel->qryendpos - pindel->qrystartpos)) ? "deletion_" : "insertion_", // Type
		(int)fabs(fabs(pindel->qryendpos - pindel->qrystartpos) - (pindel->refendpos - pindel->refstartpos)),// indel size in bp
		insertcnt, // inserted sites
		deletecnt, // deleted sites
		pindel->aligncnt, pindel->aligncnt); // XmapID1, XmapID2
      }
    }

  }/* i = start .. end-1 */
  
  if(SecondBest && !SinglelineXmap)
    delete [] Lalignment;
  delete [] Rinterval;
  delete [] indel;

  FILEclose(fp);
  free(write_buffer);

  if(fpSV){
    FILEclose(fpSV);
    free(write_bufferSV);
  }
}
