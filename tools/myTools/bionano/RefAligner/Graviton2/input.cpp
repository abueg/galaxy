#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include <unistd.h>
//#include <math.h>

#include "constants.h"
#include "globals.h"
#include "parameters.h"
#include "Ccontig.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/input.cpp 10137 2019-10-15 01:13:53Z tanantharaman $");

static Cmap **grefmap = refmap;

static int doubleInc(double *p1, double *p2)
{
  return (p1[0] > p2[0]) ? 1 : (p1[0] < p2[0]) ? -1 : 0;
}

/* merge reference sites closer than resKB */
void deres(Cmap *pmap, double resKB)
{
  Ccontig *pcontig = pmap->contig;

  for(int c=0; c < colors; c++){
    FLOAT *Y =  pmap->site[c];
    int N = pmap->numsite[c];
    if(VERB>=2){
      int cnt = 0;
      for(int I=1;I<N;I++)
	if(Y[I+1]-Y[I] <= resKB)
	  cnt++;
      printf("refid=%lld,color=%d,resKB=%0.3f(res=%0.3f,resSD=%0.3f,bpp=%0.3f):%d/%d intervals are below resKB\n",
	     pmap->id,c,resKB,res[c], resSD[c], PixelLen, cnt,N-1);
      fflush(stdout);
    }
    int shift = 0;
    for(int I=1;I <= N;I++){
      int J,K;
      for(J=I+1; J <= N; J++)
	if(Y[J]-Y[J-1] > resKB)
	  break;
      --J;

      if(VERB>= 2)
	printf("c=%d,I=%d,shift=%d,J=%d:NickLocation=%d,LabelId=%d->%d,NickCnt=%d\n",
	       c,I,shift,J,(int)floor(pmap->origsite[c][I+shift]*1000.0+0.5),
	       pmap->remap[c][I+shift],pmap->remap[c][I+shift]-shift,J-I+1);

      int rm = pmap->remap[c][I+shift] -= shift;
      if(J > I){/* merge sites Y[I] .. Y[J] */
	if(BIAS_TRACE>=2 && pmap->id == BIAS_TRACE_ID){
          #pragma omp critical
	  {
	    printf("map id=%lld:merging sites %d .. %d (at %0.10f ... %0.10f) into site %0.10f\n",
		   pmap->id,I,J,Y[I],Y[J],Yc(Y,J,J-I));
	    fflush(stdout);
	  }
	}
	//	FLOAT sum = 0.0;
	for(K=I; K <= J; K++){
	  //	  sum += Y[K];
	  pmap->NickCnt[c][K+shift] = J-I+1;
	}
	Y[I] = Yc(Y,J,J-I); // WAS sum/(J-I+1);
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
	if(pmap->siteSD[c]){/* average siteSD,sitecov,sitecnt */
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
	    pmap->ChimQuality[c][I] = -1.0;
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
	  if(cnt > 0){
	    pmap->ChimNorm[c][I] = sum / cnt;
	    if(DEBUG/* HERE >=2 */) assert(pmap->ChimNorm[c][I] >= 0.0);
	  } else
	    pmap->ChimNorm[c][I] = -1.0;
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
	  if(cnt > 0){
	    pmap->SegDupL[c][I] = sum / cnt;
	    if(DEBUG/* HERE >=2 */) assert(pmap->SegDupL[c][I] >= 0.0);
	  } else
	    pmap->SegDupL[c][I] = -1.0;
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
	  if(cnt > 0){
	    pmap->SegDupR[c][I] = sum / cnt;
	    if(DEBUG/* HERE >=2 */) assert(pmap->SegDupR[c][I] >= 0.0);
	  } else
	    pmap->SegDupR[c][I] = -1.0;
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
	  if(cnt > 0){
	    pmap->FragileEndL[c][I] = sum / cnt;
	    if(DEBUG/* HERE >=2 */) assert(pmap->FragileEndL[c][I] >= 0.0);
	  } else
	    pmap->FragileEndL[c][I] = -1.0;
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
	  if(cnt > 0){
	    pmap->FragileEndR[c][I] = sum / cnt;
	    if(DEBUG/* HERE >=2 */) assert(pmap->FragileEndR[c][I] >= 0.0);
	  } else
	    pmap->FragileEndR[c][I] = -1.0;
	}
	if(pmap->OutlierFrac[c]){/* average OutlierFrac[c] : skip -ve values */
	  float sum = 0.0;
	  int  cnt = 0;
	  for(int k = I; k <= J; k++){
	    if(pmap->OutlierFrac[c][k] >= 0.0){
	      sum += pmap->OutlierFrac[c][k];
	      cnt++;
	    }
	  }
	  if(cnt > 0){
	    pmap->OutlierFrac[c][I] = sum / cnt;
	    if(DEBUG/* HERE >=2 */) assert(pmap->OutlierFrac[c][I] >= 0.0);
	  } else
	    pmap->OutlierFrac[c][I] = -1.0;
	}

	if(pmap->Mask[c]){/* OR Mask[c] values */
	  size_t mask = 0;
	  for(int k = I; k <= J; k++)
	    mask |= pmap->Mask[c][k];
	  pmap->Mask[c][I] = mask;
	}

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
	    cnt += pmap->SNRcnt[c][k];
	    if(DEBUG && pmap->SNRcnt[c][k] && !pmap->SNRdist[c][k]){
	      printf("pmap->id=%lld:c=%d,k=%d,SNRcnt[c][k]=%d,SNRdist[c][k]=%p\n",
		     pmap->id,c,k,pmap->SNRcnt[c][k],pmap->SNRdist[c][k]);
	      fflush(stdout);
	      assert(pmap->SNRdist[c][k]);
	    }
	  }
	  double *newSNRdist = new double[cnt];/* NOTE: this memory will be leaked, but this doesn't waste much memory, since SNRcnt[] is only used for consensus maps */
	  int ncnt = 0;
	  for(int k = I; k <= J; k++)
	    for(int t = 0; t < pmap->SNRcnt[c][k]; t++)
	      newSNRdist[ncnt++] = pmap->SNRdist[c][k][t];
	  qsort(newSNRdist,cnt,sizeof(double),(intcmp*)doubleInc);	      /* sort the SNR values in ascending order */

	  if(VERB>=3 && pmap->id == 934 && c==0 && (I <= 229 && 229 <= J)){
	    printf("Replacing SNRdist[c=%d][I=%d]=%p -> %p,SNRcnt[c][I]=%d->%d,J=%d,id=%lld\n",c,I,pmap->SNRdist[c][I],newSNRdist,pmap->SNRcnt[c][I],cnt,J,pmap->id);
	    fflush(stdout);
	  }
	  if(0 /* HERE HERE blockmem >= 0 */)
	    for(int k = I; k <= J; k++)
	      delete pmap->SNRdist[c][k];
	  pmap->SNRcnt[c][I] = cnt;
	  pmap->SNRdist[c][I] = newSNRdist;
	  for(int k = I+1; k <= J; k++){
	    pmap->SNRcnt[c][k] = 0;
	    pmap->SNRdist[c][k] = NULL;
	  }
	}

	for(K=I+1;K <= J; K++){
	  if(VERB >= 2)
	    printf("c=%d,K=%d,shift=%d,J=%d:NickLocation=%d,LabelId=%d->%d,NickCnt=%d\n",
		   c,K,shift,J,(int)floor(pmap->origsite[c][K+shift]*1000.0+0.5),
		   pmap->remap[c][K+shift],pmap->remap[c][I+shift],pmap->NickCnt[c][K+shift]);
	  pmap->remap[c][K+shift] = rm;
	}

	/* shift remaining sites */
	shift += J-I;
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

	  if(pmap->siteSD[c]){
	    pmap->siteSD[c][K] = pmap->siteSD[c][J];
	    pmap->sitecov[c][K] = pmap->sitecov[c][J];
	    pmap->sitecnt[c][K] = pmap->sitecnt[c][J];
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

	  if(pmap->Mask[c])
	    pmap->Mask[c][K] = pmap->Mask[c][J];

	  if(pmap->SNRcnt[c]){
	    pmap->SNRgmean[c][K] = pmap->SNRgmean[c][J];
	    pmap->lnSNRsd[c][K] = pmap->lnSNRsd[c][J];

	    pmap->SNRcnt[c][K] = pmap->SNRcnt[c][J];
	    if(DEBUG) assert(pmap->SNRdist[c][K] == NULL);
	    pmap->SNRdist[c][K] = pmap->SNRdist[c][J];
	    pmap->SNRcnt[c][J] = 0;
	    pmap->SNRdist[c][J] = NULL;
	  }
	}

	/* copy last site at N+1 to K+1 */
	Y[K+1] = Y[N+1];
	if(pcontig){
	  pcontig->site[c][K+1] = pcontig->site[c][N+1];
	  if(pcontig->HapSite[c]){
	    if(DEBUG>=1+RELEASE) assert(pcontig->HapSite[c][N+1]==0);
	    pcontig->HapSite[c][K+1] = 0;
	  }
	}
	if(pmap->Mask[c])
	  pmap->Mask[c][K+1] = pmap->Mask[c][N+1];
	if(pmap->SNRcnt[c]){
	  pmap->SNRgmean[c][K+1] = pmap->SNRgmean[c][N+1];
	  pmap->lnSNRsd[c][K+1] = pmap->lnSNRsd[c][N+1];
	  pmap->SNRcnt[c][K+1] = pmap->SNRcnt[c][N+1];
	  if(VERB>=3){
	    printf("Updating pmap->id=%lld:SNRdist[c=%d][K+1=%d]=%p -> %p,N+1=%d\n",pmap->id,c,K+1,pmap->SNRdist[c][K+1],pmap->SNRdist[c][N+1],N+1);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(pmap->SNRdist[c][K+1] == NULL);
	  pmap->SNRdist[c][K+1] = pmap->SNRdist[c][N+1];
	  pmap->SNRcnt[c][N+1] = 0;
	  pmap->SNRdist[c][N+1] = NULL;
	}
	N = K;
	if(DEBUG)
	  assert(shift == pmap->numsite[c] - N);
      }
    }

    if(VERB>=2 && pmap->numsite[c] != N){
      printf("refid=%lld,color=%d:reduced sites from %d to %d due to res=%0.4f,resSD=%0.4f,rres=%0.2f,PixelLen=%0.3f(resKB=%0.4f)\n",
	     pmap->id,c,pmap->numsite[c],N,res[c],resSD[c],rres,PixelLen,resKB);
      fflush(stdout);
    }
    if(DEBUG) assert(Y[N+1] == Y[pmap->numsite[c] + 1]);
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

    if(DEBUG>=2){/* check that site Y[0..N+1] are monotonic */
      for(int I = 0; I <= N; I++)
	if(DEBUG && !(Y[I+1] >= Y[I])){
	  printf("refid=%lld:I=%d,N=%d:Y[I]=%0.8f,Y[I+1]=%0.8f\n",pmap->id,I,N,Y[I],Y[I+1]);
	  fflush(stdout);
	  assert(Y[I+1] >= Y[I]);
	}
    }
  }/* c = 0 .. pmap->numcolor -1 */
}

static char buf[LINESIZ];

/** input spots file and append reference map to map[0..nummaps-1] and increment refmaps count */
void input_spots(char *filename, int &numrefmaps, int &maxrefmaps, Cmap **&refmap, int spots_fileid)
{
  register char *pt;
  char *qt;
  FILE *fp;

  if(num_files + spots_fileid >= MAXFILES){
    printf("Increase MAXFILES to accomodate both -i and -ref input filenames\n");
    exit(1);
  }
  vfixx_filename[num_files+spots_fileid] = filename;

  size_t len = strlen(filename);

  /* check if input file is a .hmap file : if so call input_hmap() instead */
  if(len >= strlen(".hmap") && !strcmp(&filename[len-strlen(".hmap")],".hmap")){
    (void)input_cmap(num_files+spots_fileid,numrefmaps,maxrefmaps,refmap,1);
    return;
  }

  if(!(len >= strlen(".spots") && !strcmp(&filename[len-strlen(".spots")],".spots"))){// anything other than .spots suffix is treated as CMAP format input
    if(!(len >= strlen(".cmap") && !strcmp(&filename[len-strlen(".cmap")],".cmap"))){
      printf("WARNING: Will try to read -ref input file %s as CMAP format file\n",filename);
      fflush(stdout);
    }

    (void)input_cmap(num_files+spots_fileid,numrefmaps,maxrefmaps,refmap,0);
    return;
  }

  if((fp = Fopen(filename,"r",NFSdelay))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("Failed to read input file %s:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  int origcolors = colors;

  register int nickid,color;
  double location, reflen = -1.0;

  int maxsite = 1024;// initialize allocation size : will be adjusted dynamically, if needed

  maxmapalloc(numrefmaps+1,maxrefmaps,refmap,0,1);
  grefmap = refmap;
  register Cmap *pmap = refmap[numrefmaps++];
  pmap->mapid = pmap->id = numrefmaps;
  pmap->fileid = num_files + spots_fileid;
  for(register int c=0; c < colors; c++)
    pmap->Nickase[c] = 0;

  int linecnt=1;
  for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++){
    size_t len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      exit(1);
    }
    if(buf[0] == '#'){/* most comment lines are ignored, except "Reference Size (Bp)", 
			 "Reference Name", "N Colors" */
      char *key = (char *)"Reference Size (Bp):";
      if((pt=strstr(&buf[1],key))){
	pt += strlen(key);
	reflen = strtod(pt,&qt) * 0.001;
	if(qt == pt){
	  printf("Line with %s has invalid numerical value:\n%s\n",key,buf);
	  exit(1);
	}
	continue;
      }
      key = (char *)"Reference Name:";
      if((pt=strstr(&buf[1],key))){
	pt += strlen(key);
	while(*pt && isspace(*pt))
	  pt++;
	qt = &buf[len-1];
	while(qt > pt && isspace(*qt))
	  qt--;
	if(qt <= pt){
	  printf("Cannot find Reference Name on line %d of %s:\n%s\n", linecnt, filename, buf);
	  exit(1);
	}
	qt[1] = '\0';
	pmap->refname = strdup(pt);
	continue;
      }
      key = (char *)"N Colors:";
      if((pt=strstr(&buf[1],key))){
	pt += strlen(key);
	int Ncolors = strtol(pt,&qt,10);
	if(qt == pt){
	  printf("Line with %s has invalid integer value:\n%s\n",key,buf);
	  exit(1);
	}
	if(colors != Ncolors){
	  if(usecolor && Ncolors == 2){
	    colors = 2;/* read in Reference as 2-color than reduce to 1-color later */
	  } else {
	    printf("N Colors=%d in file %s does not match parameter value=%d\n",Ncolors,filename,colors);
	    exit(1);	  
	  }
	}
	if(Ncolors > MAXCOLOR){
	  printf("N Colors value %d is too large (MAXCOLOR=%d)\n",Ncolors,MAXCOLOR);
	  exit(1);
	}

	for(register int c=0; c < colors;c++){
	  pmap->numsite[c] = 0;
	  pmap->site[c] = new FLOAT[maxsite+2];
	  if(DEBUG>=2){
	    double init = nan("NaN");
	    for(register int i = 0; i < maxsite+2; i++)
	      pmap->site[c][i] = init;
	  }
	}
	continue;
      }

      for(register int c=0; c < colors; c++){
	char name[128];
	sprintf(name,"Nickase%d:",c+1);
	key = name;
	if((pt=strstr(&buf[1],key))){
	  pt += strlen(key);
	  while(*pt && isspace(*pt))
	    pt++;
	  for(qt=pt; *qt;qt++)
	    if(isspace(*qt)){
	      *qt = '\0';
	      break;
	    }
	  if(qt == pt){
	    printf("Line with %s is missing recognition sequence:\n%s\n",key,buf);
	    exit(1);
	  }
	  pmap->Nickase[c] = strdup(pt);
	  break;
	}

	sprintf(name,"Nickase or Nick Period%d:",c+1);
	key = name;
	if((pt=strstr(&buf[1],key))){
	  pt += strlen(key);
	  while(*pt && isspace(*pt))
	    pt++;
	  for(qt=pt; *qt;qt++)
	    if(isspace(*qt)){
	      *qt = '\0';
	      break;
	    }
	  if(qt == pt){
	    printf("Line with %s is missing recognition sequence:\n%s\n",key,buf);
	    exit(1);
	  }
	  pmap->Nickase[c] = strdup(pt);
	  break;
	}
      }
      continue;
    }
    for(register int c=0; c < colors; c++)
      if(!pmap->Nickase[c]){
	//	printf("Failed to find \"#\tNickase%d:\" in first %d lines of file %s: setting to \"unknown\"\n",c+1,linecnt,filename);
	pmap->Nickase[c] = strdup("unknown");
      }
    if(strstr(buf,"NickID") == buf){/* header line:"NickID Color Location" */
      if(!pmap->refname){
	//	printf("Failed to find \"Reference Name\" in %s: setting to \"unknown\"\n",filename);
	pmap->refname = strdup("unknown");
      }
      if(reflen < 0.0){
	printf("Failed to find \"Reference Size\" in %s\n", filename);
	exit(1);
      }
      continue;
    }

    /* read each subsequent line with NickID,Color,Location */
    pt = buf;
    nickid = strtol(pt,&qt,10);
    if(pt==qt){
      printf("unable to read NickID on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    color = strtol(pt,&qt,10);
    if(pt==qt){
      printf("unable to read Color on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    location = strtod(pt,&qt);
    if(pt==qt){
      printf("unable to read Location on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    
    if(color <= 0 || color > colors){
      printf("Invalid Color=%d on line %d of %s (for -colors %d)\n:%s\n",
	     color,linecnt,filename,colors,buf);
      exit(1);
    }

    location *= 0.001;

    if(location >= reflen - 0.0005){
      if(!(location > reflen - 0.021)){
	printf("location = %0.1f on line %d of %s is too large given that length = %0.1f\n",
	       location*1000.0,linecnt,filename,reflen*1000.0);
	exit(1);
      }
      continue;/* skip site at end : this is not really a site and end will be added later */
    }

    if(pmap->numsite[color-1] >= maxsite){/* increase allocation */
      int newmax = maxsite * 2;
      for(register int c=0; c < colors; c++){
	FLOAT *origsite = pmap->site[c];
	pmap->site[c] = new FLOAT[newmax+2];
	for(register int i = 1;i <= pmap->numsite[c];i++)
	  pmap->site[c][i] = origsite[i];
	if(DEBUG>=2){
	  double init = nan("NaN");
	  for(register int i = pmap->numsite[c]+1; i < newmax+2; i++)
	    pmap->site[c][i] = init;
	}
	delete [] origsite;
      }
      maxsite = newmax;
    }
    
    register int siteid = pmap->numsite[color-1]++;
    if(siteid > 0 && location < pmap->site[color-1][siteid] - 1e-6){
      printf("location = %0.1f on line %d of %s is larger than previous location = %0.1f\n",
	     location*1000.0, linecnt, filename, pmap->site[color-1][siteid] * 1000.0);
      exit(1);
    }
    if(VERB>=2)
      printf("nickid=%d,color=%d,location=%0.4f:siteid=%d,maxsite=%d\n",
	     nickid,color,location,siteid,maxsite);
    pmap->site[color-1][siteid+1] = location;
  }

  (void)fclose(fp);

  if(DEBUG) assert(refmap == grefmap);

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

  if(VERB){
    printf("refmaps=%d:refname=%s,colors=%d,reflen=%0.3f from %s\n",numrefmaps,pmap->refname,colors,reflen,filename);
    fflush(stdout);
  }

  for(register int c=0; c < colors; c++)
    pmap->site[c][0] = 0.0;

  if(DEBUG>=2){/* check that site Y[0..N+1] are monotonic */
    for(register int c = 0; c < colors; c++){
      register FLOAT *Y =  pmap->site[c];
      register int N = pmap->numsite[c];
      for(register int I = 0; I < N; I++){
	if(DEBUG && !(Y[I+1] >= Y[I])){
	  printf("refid=%lld:I=%d,N=%d:Y[I]=%0.8f,Y[I+1]=%0.8f\n",pmap->id,I,N,Y[I],Y[I+1]);
	  fflush(stdout);
	  assert(Y[I+1] >= Y[I]);
	}
      }
    }
  }

  if(rres > 0.0){
    double resKB = rres * 0.500;
    deres(pmap,resKB);/* merge reference sites closer than resKB */
  }

  /* append start and end locations of reference map for each color */
  for(register int c=0; c < colors; c++){
    pmap->site[c][0] = 0.0;
    pmap->site[c][pmap->numsite[c]+1] = reflen;
    if(DEBUG) assert(reflen > pmap->site[c][pmap->numsite[c]]);
  }
  pmap->startloc = 0.0;
  pmap->endloc = reflen;

  if(VERB>=2){
    for(register int c=0;c < colors; c++){
      printf("Color=%d:numsites=%d\n",c+1,pmap->numsite[c]);
      for(register int i = 0; i < pmap->numsite[c];i++)
	printf("  site[%d]=%0.3f\n",i+1,pmap->site[c][i+1]);
    }
  }

  if(origcolors == 1 && colors == 2){
    assert(usecolor >= 1 && usecolor <= 2);
    if(VERB){
      printf("Converting 2-color Reference in %s to 1-color Reference based on -usecolor %d\n", filename, usecolor);
      fflush(stdout);
    }
    if(usecolor == 2)
      pmap->colorswap(usecolor);
    colors = 1;
  }
}

/* refmap[0..numrefmaps-1]->paired are set to group ids in groupfile and the number of groups returned */
int input_groups(char *filename, Cmap **refmap, int numrefmaps)
{
  Cid2mapid *id2mapid = new Cid2mapid[numrefmaps];
  for(int i = numrefmaps; --i >= 0;){
    Cid2mapid *p = &id2mapid[i];
    p->id = refmap[i]->id;
    p->mapid = i;
  }
  qsort(id2mapid,numrefmaps,sizeof(Cid2mapid),(intcmp *)idinc);

  for(int i = numrefmaps; --i >= 0;)
    refmap[i]->paired = 0;

  FILE *fp;
  if((fp = Fopen(filename,"r",NFSdelay))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("Failed to read input file %s:errno=%d:%s\n",filename,eno,err);
    fflush(stdout);exit(1);
  }

  int linecnt = 1;
  int group = 0;
  for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++){
    int len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      exit(1);
    }
    if(buf[0] == '#')/* comment lines are ignored */
      continue;

    group++;
    if(!(0 < group && group <= MASK(31))){// NEW27
      printf("Too many groups in %s by line %d : cannot exceed max 31-bit value of %d (group %d)\n",filename,linecnt,MASK(31),group);
      fflush(stdout);exit(1);
    }

    char *pt = buf, *qt;
    for(; *pt; pt = qt){
      long long id = strtoll(pt,&qt,10);
      if((*qt && !isspace(*qt)) || id < 0){
	printf("Invalid id on line %d of %s:\n%s\n",linecnt,filename,pt);
	fflush(stdout);exit(1);
      }
      if(pt == qt)
	break;
      
      int index = findindex(id,id2mapid,numrefmaps);/* translate id to index value in id2mapid[0..numrefmaps-1] */
      if(index < 0){
	printf("ERROR:id=%lld on line %d of %s not found in %d reference maps from %s ... (%d files)\n",id,linecnt,groupfile,numrefmaps,spots_filename[0], num_reffiles);
	fflush(stdout);exit(1);// HERE : Was commented out
	continue;
      }
      int rid = id2mapid[index].mapid;
      refmap[rid]->paired = group;
      if(VERB>=2){
	printf("refmap[%d] (id=%lld) placed in group%d\n",rid,refmap[rid]->id,group);
	fflush(stdout);
      }

      if(Hapmaps){/* check for duplicate entries in id2mapid : can happen with Hapmap input and SplitHmap == 1*/
	for(int i = index+1; i < numrefmaps; i++){
	  if(id2mapid[i].id != id)
	    break;
	  int rid = id2mapid[i].mapid;
	  refmap[rid]->paired = group;
	  if(VERB>=2){
	    printf("refmap[%d] (id=%lld) placed in group%d\n",rid,refmap[rid]->id,group);
	    fflush(stdout);
	  }
	}
	for(int i = index-1; i >= 0; i--){
	  if(id2mapid[i].id != id)
	    break;
	  int rid = id2mapid[i].mapid;
	  refmap[rid]->paired = group;
	  if(VERB>=2){
	    printf("refmap[%d] (id=%lld) placed in group%d\n",rid,refmap[rid]->id,group);
	    fflush(stdout);
	  }
	}
      }
    }
  }
  fclose(fp);

  delete [] id2mapid;

  for(int i = numrefmaps; --i >= 0;)
    if(refmap[i]->paired <= 0){
      printf("refmap[%d/%d]: id=%lld not found in group file %s\n",i,numrefmaps,refmap[i]->id, filename);
      exit(1);
    }

  return group;
}
