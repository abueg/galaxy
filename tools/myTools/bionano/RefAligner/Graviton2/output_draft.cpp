#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#ifndef WIN32
#include <unistd.h>
#else

typedef int ssize_t;
#include <sys/locking.h>
#endif

#include "constants.h"
#include "globals.h"
#include "parameters.h"
#include "Ccontig.h"

#define COVFIX 1 /* Fix CMAP output coverage to represet coverage of current site */

#define MIN_MAPWT RefineMinMapwt /* minimum map weight to output in refined .xmap */

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/output_draft.cpp 11297 2020-07-09 19:48:45Z tanantharaman $");
extern void printversion(FILE *fp);

typedef int intcmp (const void *,const void *);

#define BVERB 0 /* Verbose display of coverage profile and breakpoint decisions */

// SplitRev  :  TEST 0 , 1 , 2 or 3  : New experimental version of Contig Split
// SplitSite : Only check coverage at consensus sites (and ignores less reliable values in between consensus sites)

extern int MaxImageFOV;

struct Cmaporder {
  int m;/* index on contig->contig[] etc */
  double left;/* leftmost aligned site location on reference (any color) */
};

static int MaporderSiteInc(Cmaporder *p1, Cmaporder *p2)
{
  return (p1->left > p2->left) ? 1 : (p1->left < p2->left) ? -1 : 0;
}

/* sort doubles in increasing order */
static int doubleInc(double *p1, double *p2)
{
  return (p1[0] > p2[0]) ? 1 : (p1[0] < p2[0]) ? -1 : 0;
}

/* sort floats in decreasing order */
static int floatdec(register float *p1, register float *p2)
{
  register float id1 = *p1;
  register float id2 = *p2;
  
  return (id1 < id2) ? 1 : (id1 > id2) ? -1 : 0;
}

/* given an entry in Ccontig->sitemap, its final site index is only known by iterating through the peak array
   So loop over the non-zero entries of peak, returning the index which corresponds to peakidx */
int getSiteIdx(int *peakind, int peaksize, int peakidx) 
{
  for(int i=0; i<=peaksize; i++)
    if( peakind[i] == peakidx )
      return i+1; //plus one becuase site numbers start at 1
  return 0; //not found--I think this is ok if the right end is trimmed, but it's still not a site in the refined map, so return 0
}

extern int noQuality; //parameters.h

/* output specified region (sites = left .. right) of a refined contig */
static void output_draftIO(Ccontig *contig, long long contigid, char *basename, int left, int right, int n, int *peak, double *Y, int Allele, double newleft, double newright)
{
  //  double leftend = (left > 0 ? MININTERVAL : 0.0);
  //  double rightend = (right <= n ? MININTERVAL : 0.0);
  if(Refine == 1 /* WAS && extend <= 1*/){/* trim Y[left .. Yright] down to real sites */
    while(left < right && !peak[left])
      left++;
    while(left < right && !peak[right])
      right--;
  }

  double leftend = (left > 0) ? (Refine == 1 /* WAS && extend <= 1*/ && EndTrimCov <= 0.0) ? max(MININTERVAL,Y[left]-contig->left) : MININTERVAL : 0.0;/* extension on left end beyond Y[left] */
  double rightend = (right <= n) ? (Refine == 1 /* WAS && extend <= 1*/ && EndTrimCov <= 0.0) ? max(MININTERVAL,contig->right - Y[right]) : MININTERVAL : 0.0;/* extension on right end beyond Y[right] */
  if(DEBUG) assert(leftend >= 0.0);
  if(DEBUG) assert(rightend >= 0.0);

  if(DEBUG>=2){/* check that Y[0..n+1] is monotonic */
    for(int i = 1; i <= n+1; i++){
      if(!(Y[i] >= Y[i-1])){
	printf("output_draftIO: i=%d,n=%d:Y[i-1]= %0.4f, Y[i]= %0.4f, Y[n+1]= %0.4f\n",i,n,Y[i-1],Y[i],Y[n+1]);
	fflush(stdout);
	assert(Y[i] >= Y[i-1]);
      }
    }
  }

  if(VERB>=2 && Refine==1){
    printf("Y[left]=%0.6f,contig->left=%0.6f,leftend=%0.6f, Y[right]=%0.6f,contig->right=%0.6f,rightend=%0.6f\n",
	   Y[left],contig->left,leftend,Y[right],contig->right,rightend);
    fflush(stdout);
  }

  if(strstr(basename,"/dev/null"))
    return;

  char filename[PATH_MAX];
  FILE *fp;

  if(VERB && Refine){
    printf("output_draftIO:CMapID=%lld,contigid=%lld,left=%d(%0.3f),right=%d(%0.3f),n=%d,leftend=%0.3f,rightend=%0.3f,Allele=%d,HapSite[0]=%p,MaskL=0x%lx,MaskR=0x%lx,EndTrim=%0.2f\n",
	   CMapID,contigid,left,Y[left],right,Y[right],n,leftend,rightend,Allele,contig->HapSite[0],contig->MaskL,contig->MaskR,EndTrimCov);
    fflush(stdout);
  }

  long long origCMapID = CMapID;
  if(CMapID > 0 && HapSitePvalue > 0.0)
    CMapID = CMapID * 10 + Allele;
  
  /* output consensus map as .cmap file */
  strcpy(filename,basename);
  /* no need to strip suffix : for default basename this was already done in RefAligner.cpp */
  if(1){
    int i = strlen(filename);
    if(contigid > 0)
      sprintf(&filename[i],"_contig%lld.cmap", contigid);
    else if(CMapID > 0){
      if(contig->HapSite[0])
	sprintf(&filename[i],"_contig%lld.hmap",CMapID);
      else
	sprintf(&filename[i],"_contig%lld.cmap",CMapID);
    } else
      sprintf(&filename[i],".cmap");
  }

  if(checkFile(filename)){
    CMapID = origCMapID;
    return;
  }

  if(VERB){
    if(contigid > 0) {
      if(Refine)
	printf("Generating %s (Assembler Refined contig%lld in .cmap format)\n",filename, contigid);
      else
	printf("Generating %s (draft contig%lld in .cmap format)\n",filename, contigid);
    } else if(CMapID > 0) {
      if(Allele)
	printf("Generating %s (refined reference of contig%lld Allele %d in .cmap format):CmapSNR=%d,MapSNR=%d\n",filename, origCMapID,Allele, CmapSNR,MapSNR);
      else if(contig->HapSite[0])
	printf("Generating %s (haplotyped refined reference of contig%lld in .hmap format):CmapSNR=%d,MapSNR=%d\n",filename,CMapID,CmapSNR,MapSNR);
      else
	printf("Generating %s (refined reference of contig%lld in .cmap format):CmapSNR=%d,MapSNR=%d\n",filename,CMapID,CmapSNR,MapSNR);
    } else
      printf("Generating %s (refined reference %s in .cmap format)\n",filename, basename);
    fflush(stdout);
  }

  char *origfilename = strdup(filename);
  if(1 /* NEW contig->HapSite[0] */){
    int len = strlen(filename);
    sprintf(&filename[len],".tmp");
  }
  char *tmpfilename = strdup(filename);

  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("failed to open file %s for writing contig%lld draft consensus(as .cmap):errno=%d:%s\n",filename,contigid,eno,err);
    fflush(stdout);exit(1);
  }

  if(Refine /* WAS175 >= 2 */ && EndTrimRatio > 0.0 && contigid <= 0){
    // Trim back the last K labels at either end if highest label count for all of these K labels is <= EndTrimRatio times the label cnt of the next label.

    if(colors >= 2){
      printf("ERROR : -EndTrimRatio not implemented for %d colors\n",colors);
      fflush(stdout);exit(1);
    }

    int nleft = left;
    float maxocc = contig->sitecnt[0][left], maxleft = 0.0, maxright = 0.0;

    for(int i = left;  ++i <= right;){
      if(!peak[i])
	continue;
      float occ = contig->sitecnt[0][i];
      if(maxocc <= occ * EndTrimRatio){
	nleft = i;
	maxleft = maxocc;
      }
      maxocc = max(occ, maxocc);
    }

    if(nleft > left){
      if(VERB){
	printf("Trimmed off labels %d..%d on left end leaving first label at site[%d]= %0.3f (max occ[%d..%d]= %0.2f, occ[%d]= %0.2f)\n",
	       left,nleft-1,nleft,contig->site[0][nleft],left,nleft-1,maxleft,nleft,contig->sitecnt[0][nleft]);
	fflush(stdout);
      }

      leftend += contig->site[0][nleft] - contig->site[0][left];// NEW281
      left = nleft;
    }

    int nright = right;
    maxocc = contig->sitecnt[0][right];

    for(int i = right; --i >= left;){
      if(!peak[i])
	continue;
      float occ = contig->sitecnt[0][i];
      if(VERB>=2){
	printf("i=%d,right=%d:occ[i]= %0.2f, maxocc[%d..%d]= %0.2f, EndTrimRatio= %0.3f: nright= %d -> %d\n",
	       i,right,occ,i+1,right,maxocc,EndTrimRatio, nright, (maxocc <= occ * EndTrimRatio) ? i : nright);
	fflush(stdout);
      }
      if(maxocc <= occ * EndTrimRatio){
	nright = i;
	maxright = maxocc;
      }
      maxocc = max(occ, maxocc);
    }
    
    if(nright < right){
      if(VERB){
	printf("Trimmed off labels %d..%d on right end leaving last label at Y[%d]= %0.3f : max occ[%d..%d]= %0.2f, occ[%d]= %0.2f\n",
	       nright+1,right,nright,contig->site[0][nright],nright+1,right,maxright,nright,contig->sitecnt[0][nright]);
	fflush(stdout);
      }
      rightend += contig->site[0][right] - contig->site[0][nright];// NEW280
      right = nright;
    }
  }

  register int NumSites = 0;
  for(register int k = max(1,left); k <= min(n,right); k++)
    if(peak[k] > 0)
      NumSites++;

  if(!noQuality && contig->HapSite[0]){
    printf("Warning: -contigQual not implemented with -Haplotype\n");
    fflush(stdout);
  }
  if(ChimQuality && contig->HapSite[0]){
    printf("Warning: -ChimQuality not implemented with -Haplotype for .hmap output\n");
    fflush(stdout);
  }

  //Compute quality score
  double* avgmolpv    = 0; //this will store the sum of mol pvalues for each molecule which aligns to this site, weighted
  double* avgmolpvden = 0; //denominator to get average--divide avgmolpv by this
  double totavg = 0;
  if(!noQuality && contig->mapWT && contig->logPV && contig->LP && contig->sitemapL[0] && !contig->HapSite[0]) {
    int* peakindices = new int[NumSites+1]; //store peak indices by looping through peak again
    int totalign = 0; //count indicies in peak, reuse below
    //printf("NUMSITES = %i  left = %i  right = %i\n", NumSites, left, right);
    for(register int k = max(1,left); k <= min(n,right); k++)
      if(peak[k] > 0) {
	//printf("peakindices[%i] = %i\n", totalign, k);
	peakindices[totalign] = k; //note--start at zero, so indices in this list are off by 1 from indices in cmap
	totalign++;
      }
    avgmolpv    = new double[NumSites+1]; //because site numbers go from 1 to N, get one extra, and don't use 0 
    avgmolpvden = new double[NumSites+1];
    for(int i=0; i<=NumSites; i++) { 
      avgmolpv[i] = 0;
      avgmolpvden[i] = 0;
    }
    for(int i= 0; i < contig->nummaps; i++) { //loop over all the molecules which align to this contig
      //long long id = Gmap[contig->contig[i].mapid]->id; //printf below
      int moltot = 0; //count how many sites are in peak for this molecule -- need for filtering below
      for(int j= 0; j <= contig->contig[i].numsite[0] + 1; j++) {
	int peakidx = contig->sitemap[0][i][j]; //also idx in fragcnt, sitecnt, etc
	if(peakidx > -1) { //-1 is not aligned--ignore those
	  int hit = (peak[peakidx] > 0 ? getSiteIdx(peakindices, NumSites, peakidx) : 0); //getSiteIdx is index in final contig--need to check if trimmed
	  moltot += (hit > 0 ? 1 : 0);
	  //printf("SITEMAP id %lld idx %i (%2i) cov %5.2f wcov %5.2f occ %5.2f\n", id, peakidx, hit, contig->fragcnt[0][peakidx], contig->sitecntFN[0][peakidx], contig->sitecnt[0][peakidx]); //printing here prints unused mols
	}
      }
      //check pvalue, score, numsites after refinement
      double confwt = contig->logPV[i] * contig->mapWT[i] / moltot;
      if(contig->LP[i] <= ScoreThreshold || contig->logPV[i] <= LogPvThreshold || moltot < AlignedSiteThreshold || moltot <= 0) 
	continue;
      //count sites in above loop, populate avgmolpv{,den} in this loop
      for(int j= 0; j <= contig->contig[i].numsite[0] + 1; j++) {
	int peakidx = contig->sitemap[0][i][j]; //also idx in fragcnt, sitecnt, etc
	if(peakidx < 0 || peak[peakidx] <= 0)  //see above loop
	  continue;
	int hit = getSiteIdx(peakindices, NumSites, peakidx); //getSiteIdx is index in final contig--need to check if trimmed
	int peakidxl = contig->sitemapL[0][i][j]; //check for another site (don't count it in denominator, only numerator)
	int hitl = (peakidxl != peakidx ? getSiteIdx(peakindices, NumSites, peakidxl) : 0);
	//printf("SITEMAP id %lld idx %3i (%2i) (%i) cov %5.2f wcov %5.2f occ %5.2f\n", id, peakidx, hit, hitl, contig->fragcnt[0][peakidx], contig->sitecntFN[0][peakidx], contig->sitecnt[0][peakidx]); //no unused mols
	if(hit > 0) {
	  avgmolpv[hit] += confwt;
	  avgmolpvden[hit] += contig->mapWT[i];
	}
	if(hitl > 0) {
	  avgmolpv[hitl] += confwt;
	  avgmolpvden[hitl] += contig->mapWT[i];
	}
      }
      //printf(" id %10lld : %2i  weight %.3f  score %5.2f  logPV %5.2f  /n  %4.2f  *wt  %4.2f\n", id, moltot, contig->mapWT[i], contig->LP[i], contig->logPV[i], contig->logPV[i]/moltot, confwt);
    }
  } //end if noQuality
  
  size_t blen = 10*1024*1024;/* matches scif_io buffer size */
  char *write_buffer = (char *)malloc(blen);
  if(!write_buffer){
    printf("output_draftIO:malloc(%llu) failed\n",(unsigned long long)blen);
    fflush(stdout);exit(1);
  }
  setbuffer(fp, write_buffer, blen);

  /* write out commandline */
  printversion(fp);

  size_t *Mask = new size_t[n+2];
  memset(Mask,0,(n+2)*sizeof(size_t));

  if(contig->SegDupCnt > 0 && SegDupMask){/* recompute Mask[0..n+1] */
    double len = newright - newleft;

    for(int S = 0; S < contig->SegDupCnt; S++){
      double start = newleft + contig->SegDupStart[S] * len;
      double end = newleft + contig->SegDupEnd[S] * len;
      
      if(VERB/* HERE >=2 */){
	printf("output_draftIO:CMapID=%lld:SegDupCnt=%d,S=%d:SegDupStart=%0.6f,SegDupEnd=%0.6f:len=%0.4f,start=%0.4f(%0.4f),end=%0.4f(%0.4f),newleft=%0.4f,newright=%0.4f\n",
	       CMapID,contig->SegDupCnt,S,contig->SegDupStart[S],contig->SegDupEnd[S],len,start,start-newleft,end,end-newleft,newleft,newright);
	fflush(stdout);
      }

      for(int i = 1; i <= n; i++){
	if(Y[i] < start)
	  continue;
	if(Y[i] > end)
	  break;
	Mask[i] = SegDupMask;
      }
    }
  }

  int mapSNR = MapSNR ? 1 : 0;

  /* compute SNR mean,sd etc from alignments to this reference */
  double *SNRwtsum[MAXCOLOR];/* NOTE : this is only used for SNR computation : it is not the Occurrence value output in .cmap */
  double *SNRsum[MAXCOLOR];
  double *SNRsumsq[MAXCOLOR];
  int *SNRcnt[MAXCOLOR];
  int *SNRmax[MAXCOLOR];
  double **SNR[MAXCOLOR];
  int *nn = contig->numsite;
  if(DEBUG) assert(nn[0] == n);
  if(DEBUG && CovNorm) assert(!REPLACE_COV);

  double CovLambdaInv = 1.0;
  if(CovNorm){/* compute adjustment to fragcntT[] and (sitecnt[]) based on fragcntTNorm[] */
    CovLambdaInv = 1.0/CovLambda;
    for(int t = 0; t <= n; t++){
      float cov = contig->fragcntT[0][t];
      if(DEBUG && !isfinite(cov)){
	printf("CMapID=%lld:t=%d,n=%d,contig->fragcntT[0][t]=%0.6f\n",CMapID,t,n,contig->fragcntT[0][t]);
	fflush(stdout);
	assert(isfinite(contig->fragcntT[0][t]));
      }
      if(cov <= 0.001){
	contig->fragcntTnorm[0][t] = 1.0;
	continue;
      }
      double mean = contig->fragcntTnorm[0][t]/cov;/* mean size of maps contributing to contig->fragcntT[0][t] */
      if(DEBUG && !isfinite(mean)){
	printf("CMapID=%lld:t=%d,n=%d,contig->fragcntT[0][t]=%0.6f,contig->fragcntTnorm[0][t]=%0.6f,mean=%0.6f\n",CMapID,t,n,contig->fragcntT[0][t],contig->fragcntTnorm[0][t],mean);
	fflush(stdout);
      }
      if(mean > 1000.0){
	printf("WARNING:CMapID=%lld:t=%d,n=%d,contig->fragcntT[0][t]=%0.6f,mean=%0.6f,MinLen=%0.6f,CovLambda=%0.6f,fragcntTnorm[0][t]=%0.6f\n",
	       CMapID,t,n,contig->fragcntT[0][t],mean,MinLen,CovLambda,contig->fragcntTnorm[0][t]);
	fflush(stdout);
      }
      mean = min(mean,1000.0);// HERE : stopgap fix until we have time to track down cause of large values above 1000.0 in run5.log32
      contig->fragcntTnorm[0][t] = exp((mean-MinLen-CovLambda)*CovLambdaInv);
      if(DEBUG && !isfinite(contig->fragcntTnorm[0][t])){
	printf("CMapID=%lld:t=%d,n=%d,contig->fragcntT[0][t]=%0.6f,mean=%0.6f,MinLen=%0.6f,CovLambda=%0.6f,fragcntTnorm[0][t]->%0.6f\n",
	       CMapID,t,n,contig->fragcntT[0][t],mean,MinLen,CovLambda,contig->fragcntTnorm[0][t]);
	fflush(stdout);
	assert(isfinite(contig->fragcntTnorm[0][t]));
      }
    }
  }

  if(mapSNR){
    for(int c = 0; c < colors; c++){
      SNRwtsum[c] = new double[nn[c]+1];/* SNRwtsum[1..nn[c]] weighted sum (site based) */
      if(CmapSNR){
	SNRcnt[c] = new int[nn[c]+1];/* SNR count per site */
	SNRmax[c] = new int[nn[c]+1];/* SNR count max allocation per site */
	SNR[c] = new double*[nn[c]+1];/* SNR values array per site */
      }
      SNRsum[c] = new double[nn[c]+1];
      SNRsumsq[c] = new double[nn[c]+1];
      for(register int t = 1 ; t <= nn[c]; t++){
	SNRwtsum[c][t] = 0.0;
	SNRsum[c][t] = SNRsumsq[c][t] = 0.0;
	if(CmapSNR){
	  SNR[c][t] = new double[SNRmax[c][t] = 64];
	  SNRcnt[c][t] = 0;
	}
      }
    }

    /* loop over all maps in contig using alignment information in sitemap[] to compute SNR*[] values */

    for(int m = 0; m < contig->nummaps; m++){
      double mapwt = contig->mapWT ? contig->mapWT[m] : 1.0;
      if(DEBUG/* >=2 */ && !isfinite(mapwt)){
	printf("m=%d,contig->nummaps=%d,contig->mapWT[m]= %0.6e\n",m,contig->nummaps,contig->mapWT[m]);
	fflush(stdout);
	assert(isfinite(mapwt));
      }
      if(mapwt <= 0.001)
	continue;
      Ccontig *mcontig = &contig->contig[m];
      int flip = contig->flip[m];
      int mapid = mcontig->mapid;
      Cmap *pmap = Gmap[mapid];

      for(int c = 0; c < colors; c++){
	int M = mcontig->numsite[c];
	int left  = mcontig->trimL[c];/* left end of map sub-region in use for map m (0 if complete map is in use) */
	int right = mcontig->trimR[c];/* right end of map sub-region in use for map m (pmap->numsites[0]+1 if complete map is in use) */
	int *mapR = contig->sitemap[c][m];/* Gmap[j=0..M+1] maps site j of color c of map m to one of contig->numsite[c] sites : 0 refers to contig's left end, contig->numsite[c]+1 to the right end */
	int *mapL = contig->sitemapL[c] ? contig->sitemapL[c][m] : mapR;/* when site j of color c of map m aligns to multiple consensus sites, mapL[j] is the leftmost consensus site index, while mapR[j] is the rightmost */

	for(int J = 1; J <= M; J++){
	  int R = mapR[J];
	  if(R < 0)
	    continue;
	  if(DEBUG) assert(0 < R && R <= contig->numsite[c]);
	  int L = mapL[J];

	  if(DEBUG/* HERE >=2 */) assert(0 < L && L <= R);
	  if(DEBUG/* HERE >=2 */ && flip==0) assert(J+left <= pmap->numsite[c]);
	  if(DEBUG/* HERE >=2 */ && flip==1) assert(1 <= right-J && right-J <= pmap->numsite[c]);
	  double wt = mapwt;// mapwt/(1.0 + (resSD > 0 ? R-L : 0));
	  for(int t = L; t <= R; t++){
	    SNRwtsum[c][t] += wt;
	    double snr = pmap->SNR[c][flip ? right-J : J+left];
	    if(DEBUG>=2 && ! (isfinite(snr) && snr > 0.0)){
	      printf("m=%d/%d,flip=%d,mapid=%d(id=%lld):c=%d,left=%d,right=%d,M=%d,numsite[c]=%d,J=%d:L=%d,R=%d,wt=%0.6f:t=%d,snr=%g\n",
		     m, contig->nummaps,flip,mapid, pmap->id,c,left,right,M,pmap->numsite[c],J,L,R,wt,t,snr);
	      for(int j = 1; j <= pmap->numsite[c]; j++)
		printf("pmap->SNR[c][%d] = %0.8e\n",j,pmap->SNR[c][j]);
	      fflush(stdout);
	      assert(isfinite(snr) && snr > 0.0);
	    }
	    double lnsnr = log(snr);
	    if(DEBUG>=2) assert(isfinite(lnsnr));
	    if(CmapSNR && wt > 0.5){
	      if(SNRcnt[c][t] >= SNRmax[c][t]){
		int newmax = SNRmax[c][t]*2;
		double *oldSNR = SNR[c][t];
		SNR[c][t] = new double[newmax];
		for(int u = SNRmax[c][t]; --u >= 0; )
		  SNR[c][t][u] = oldSNR[u];
		delete [] oldSNR;
		SNRmax[c][t] = newmax;
	      }
	      SNR[c][t][SNRcnt[c][t]++] = snr;
	    }
	      
	    SNRsum[c][t] += wt * lnsnr;
	    if(DEBUG/* HERE >=2 */) assert(isfinite(SNRsum[c][t]));
	    SNRsumsq[c][t] += wt * lnsnr*lnsnr;
	    if(DEBUG/* HERE >=2 */) assert(isfinite(SNRsumsq[c][t]));
	  }
	} // J = 1 .. M
      }
    }
  }

  if(VERB>=2){
    printf("Refine= %d, EndTrimRatio= %0.3f, contigid=%lld: left=%d,right=%d\n", Refine, EndTrimRatio, contigid, left, right);
    fflush(stdout);
  }

  if(Refine/* WAS172 >= 2*/ && EndTrimCov > 0.0 && contigid <= 0){ /* Extend ends based on unaligned molecule ends  : EndTrimFrac'th percentile or EndTrimCov'th length in list sorted in descending order */
    double *leftList = new double[contig->nummaps * 4];
    double *rightList = &leftList[contig->nummaps];
    double *leftListFull = &leftList[contig->nummaps*2];
    double *rightListFull = &leftList[contig->nummaps*3];

    int leftCnt = 0, rightCnt = 0;

    double minwt = (EndTrimWt > 0.0) ? EndTrimWt : MIN_MAPWT;

    for(int m = 0; m < contig->nummaps; m++){
      double mapwt = contig->mapWT ? contig->mapWT[m] : 1.0;
      if(DEBUG/* >=2 */ && !isfinite(mapwt)){
	printf("m=%d,contig->nummaps=%d,contig->mapWT[m]= %0.6e\n",m,contig->nummaps,contig->mapWT[m]);
	fflush(stdout);
	assert(isfinite(mapwt));
      }
      if(VERB>=2){
	printf("m=%d/%d: mapwt= %0.6f, minwt= %0.6f\n",m, contig->nummaps, mapwt, minwt);
	fflush(stdout);
      }
      if(mapwt < minwt)
	continue;

      Ccontig *mcontig = &contig->contig[m];
      int flip = contig->flip[m];
      int mapid = mcontig->mapid;
      Cmap *pmap = Gmap[mapid];
      double leftlen = -1000.0, rightlen = -1000.0;
      double leftlenFull = -1000.0, rightlenFull = -1000.0;
    
      for(int c = 0; c < colors; c++){
	int M = pmap->numsite[c];
	FLOAT *X = pmap->site[c];
	int *mapR = contig->sitemap[c][m];/* map[j=0..M+1] maps site j of color c of map m to one of contig->numsite[c] sites : 0 refers to contig's left end, contig->numsite[c]+1 to the right end */
	int *mapL = contig->sitemapL[c] ? contig->sitemapL[c][m] : mapR;/* when site j of color c of map m aligns to multiple consensus sites, mapL[j] is the leftmost consensus site index, while mapR[j] is the rightmost */
      
	int JL = M, JR = 0;
      
	for(int J = 1; J <= M; J++){
	  int R = mapR[J];
	  if(R < 0)
	    continue;
	  int L = mapL[J];
	  if(DEBUG) assert(0 < L && L <= R && R <= contig->numsite[c]);
	  
	  if(left <= R && R <= right){
	    JL = min(JL, J);
	    JR = J;
	  }
	}

	if(JL <= JR){
	  if(DEBUG/* HERE >=2 */) assert(0 < JL && JL <= JR && JR <= M);

	  int IL = mapR[JL], KL = mapL[JL];
	  int IR = mapR[JR], KR = mapL[JR];
	
	  /* lengths to end or label after last aligned label (whichever comes first) */
	  leftlen = max(leftlen, (flip ? X[M+1-JL+1] - X[M+1-JL] : X[JL] - X[JL-1]) - (Yc(Y,IL,IL-KL) - Y[left]));
	  rightlen = max(rightlen, (flip ? X[M+1-JR] - X[M-JR] : X[JR+1] - X[JR]) - (Y[right] - Yc(Y,IR,IR-KR)));

	  /* lengths to end */
	  leftlenFull = max(leftlenFull, (flip ? X[M+1] - X[M+1-JL] : X[JL]) - (Yc(Y,IL,IL-KL) - Y[left]));
	  rightlenFull = max(rightlenFull, (flip ? X[M+1-JR] : X[M+1] - X[JR]) - (Y[right] - Yc(Y,IR,IR-KR)));

	  if(VERB>=2){
	    printf("m=%d/%d(id=%lld),wt=%0.6f,c=%d,M=%d,n=%d,flip=%d: IL=%d,KL=%d,JL=%d,IR=%d,KR=%d,JR=%d, left=%d, right=%d : leftlen=%0.3f, rightlen=%0.3f(full= %0.3f,%0.3f),len=%0.3f\n",
		   m,contig->nummaps,pmap->id,mapwt,c,M,n,flip,IL,IL-KL,JL,IR,IR-KR,JR,left,right,leftlen,rightlen,leftlenFull,rightlenFull,X[M+1]);
	    fflush(stdout);
	  }
	}
      }

      if(leftlen > 0.0){
	leftListFull[leftCnt] = leftlenFull;
	leftList[leftCnt++] = leftlen;
      }

      if(rightlen > 0.0){
	rightListFull[rightCnt] = rightlenFull;
	rightList[rightCnt++] = rightlen;
      }
    }

    if(leftCnt > EndTrimCov){/* adjust leftend size */
      qsort(leftList, leftCnt, sizeof(double), (intcmp*) doubleInc);
    
      int rank = max(floor(EndTrimCov + 0.999999), floor(leftCnt * EndTrimFrac + 0.5)); // WAS11 max(floor(EndTrimCov + 0.9), floor(leftCnt * 0.5 + 0.5));
      if(DEBUG) assert(1 <= rank && rank <= leftCnt);
      double end = leftList[leftCnt - rank];
      if(VERB/* HERE >=2 */){
	printf("EndTrimCov = %0.3f(R=%0.4f), leftCnt = %d : rank = %d, new leftend= %0.3f\n",EndTrimCov,EndTrimFrac,leftCnt,rank, end);
	fflush(stdout);
      }
      double origleftend = leftend;
      leftend = max(end, leftend);

      if(EndTrimCnt > 0.0){
	int rank50 = max(floor(EndTrimCov + 0.999999), floor(leftCnt * 0.5 + 0.5));
	if(DEBUG) assert(1 <= rank50 && rank50 <= leftCnt);
	double end50 = leftList[leftCnt - rank50];
	double leftend50 = max(end50, origleftend);

	double maxerr = 0.0;/* determine sizing error = sqrt(sf^2 + (sr * leftend)^2), for whichever color has larger error */
	for(int c = 0; c < colors; c++){
	  double err2 = SR[c] * leftend50;
	  maxerr = max(maxerr, sqrt(SF[c]*SF[c] + SD[c]*fabs(SD[c])*leftend50 + err2*err2));
	}
	maxerr *= 3.0;/* use 3 x SD as error bound */

	qsort(leftListFull, leftCnt, sizeof(double), (intcmp*)doubleInc);
	
	int cnt = 0;
	for(int i = 0; i < leftCnt; i++){
	  double err = fabs(leftend50 - leftListFull[i]);
	  if(err < maxerr)
	    cnt++;
	  if(VERB>=3){
	    printf("i=%d/%d: leftList[i]= %0.4f, leftListFull[i]= %0.4f, leftend50= %0.4f, err= %0.4f, maxerr= %0.4f: cnt= %d\n",
		   i,leftCnt,leftList[i], leftListFull[i], leftend50, err, maxerr, cnt);
	    fflush(stdout);
	  }
	}
	if(cnt >= EndTrimCnt){
	  contig->MaskL |= (END_NOEXT | END_CHR);
	  leftend = min(leftend, leftend50);
	}
	if(VERB/* HERE >=2 */){
	  printf("EndTrimCnt= %0.1f : left maxerr= %0.3f, cnt= %d, MaskL -> %lu, leftend -> %0.4f\n", EndTrimCnt,maxerr, cnt, contig->MaskL, leftend);
	  fflush(stdout);
	}
      }
    }

    if(rightCnt > EndTrimCov){/* adjust rightend size */
      qsort(rightList, rightCnt, sizeof(double), (intcmp*) doubleInc);
    
      int rank = max(floor(EndTrimCov + 0.999999), floor(rightCnt * EndTrimFrac + 0.5));// WAS11 max(floor(EndTrimCov + 0.9), floor(rightCnt * 0.5 + 0.5));
      if(DEBUG) assert(1 <= rank && rank <= rightCnt);
      double end = rightList[rightCnt - rank];
      if(VERB/* HERE >=2 */){
	printf("EndTrimCov = %0.3f(R=%0.4f), rightCnt = %d : rank = %d, new rightend= %0.3f\n",EndTrimCov,EndTrimFrac,rightCnt,rank, end);
	fflush(stdout);
      }
      double origrightend = rightend;
      rightend = max(end, rightend);

      if(EndTrimCnt > 0.0){
	int rank50 = max(floor(EndTrimCov + 0.999999), floor(rightCnt * 0.5 + 0.5));
	if(DEBUG) assert(1 <= rank50 && rank50 <= rightCnt);
	double end50 = rightList[rightCnt - rank50];
	double rightend50 = max(end50, origrightend);

	double maxerr = 0.0;/* determine sizing error = sqrt(sf^2 + (sr * rightend50)^2), for whichever color has larger error */
	for(int c = 0; c < colors; c++){
	  double err2 = SR[c] * rightend50;
	  maxerr = max(maxerr, sqrt(SF[c]*SF[c] + SD[c]*fabs(SD[c])*rightend50 + err2*err2));
	}
	maxerr *= 3.0;/* use 3 x SD as error bound */

	qsort(rightListFull, rightCnt, sizeof(double), (intcmp*)doubleInc);
	
	int cnt = 0;
	for(int i = 0; i < rightCnt; i++){
	  double err = fabs(rightend50 - rightListFull[i]);
	  if(err < maxerr)
	    cnt++;
	  if(VERB>=3){
	    printf("i=%d/%d: rightList[i]= %0.4f, rightListFull[i]= %0.4f, rightend50= %0.4f, err= %0.4f, maxerr= %0.4f: cnt= %d\n",
		   i,rightCnt,rightList[i], rightListFull[i], rightend50,err,maxerr,cnt);
	    fflush(stdout);
	  }
	}
	if(cnt >= EndTrimCnt){
	  contig->MaskR |= (END_NOEXT | END_CHR);
	  rightend = min(rightend, rightend50);
	}
	if(VERB/* HERE >=2 */){
	  printf("EndTrimCnt= %0.1f : right maxerr= %0.3f, cnt= %d, MaskR -> %lu, rightend -> %0.4f\n",
		 EndTrimCnt,maxerr, cnt, contig->MaskR, rightend);
	  fflush(stdout);
	}
      }
    }

    delete [] leftList;
    leftList = rightList = NULL;
    leftCnt = rightCnt = 0;
  }

  /* remove MaskL or MaskR if more than 50kb has been trimmed off original ends */
  if(contig->MaskL && Y[left]-leftend > newleft + 50.0){
    if(VERB/* HERE >=2 */){
      printf(" MaskL= 0x%lx -> 0 (Y[left=%d]=%0.4f,leftend=%0.4f,newleft=%0.4f,contig->left=%0.4f): More than 50kb trimmed\n",contig->MaskL,left,Y[left],leftend,newleft,contig->left);
      fflush(stdout);
    }
    contig->MaskL = 0;
  }
  if(contig->MaskR && Y[right]+rightend < newright - 50.0){
    if(VERB/* HERE >=2 */){
      printf(" MaskR= 0x%lx -> 0 (Y[right=%d]=%0.4f,rightend=%0.4f,newright=%0.4f,contig->right=%0.4f): More than 50kb trimmed\n",contig->MaskR,right,Y[right],rightend,newright,contig->right);
      fflush(stdout);
    }
    contig->MaskR = 0;
  }

  /* output key header lines */
  if(ChimQuality >= 2 && contig->sitecntN2[0] && contig->sitecntN3[0] && contig->sitecntN4[0] && contig->sitecntN5[0])
    fprintf(fp,"# CMAP File Version:\t0.2\n");
  else
    fprintf(fp,"# CMAP File Version:\t0.1\n");
  fprintf(fp,"# Label Channels:\t%d\n", colors);
  if(nummaps > 0 && Gmap[0]){
    if(usecolor >= 2 && colors == 1){// NEW71
      if(Gmap[0]->Nickase[usecolor-1])
	fprintf(fp,"# Nickase Recognition Site %d:\t%s\n", usecolor, Gmap[0]->Nickase[usecolor-1]);
      else
	fprintf(fp,"# Nickase Recognition Site %d:\tunknown\n", usecolor);
    } else {
      for(register int c = 0; c < colors; c++)
	if(Gmap[0]->Nickase[c])
	  fprintf(fp,"# Nickase Recognition Site %d:\t%s\n", c+1, Gmap[0]->Nickase[c]);
	else
	  fprintf(fp,"# Nickase Recognition Site %d:\tunknown\n", c+1);
    }
  }
  fprintf(fp,"# Number of Consensus Maps:\t1\n");
  if(VERB/* HERE >=2 */){
    printf("mapSNR=%d,noQuality=%d,CmapSNR=%d\n",mapSNR,noQuality,CmapSNR);
    fflush(stdout);
  }
  if(ChimQuality && contig->sitecntFN[0] && !contig->HapSite[0]){
    const char *OutlierString = (OutlierQuality < 2) ? "OutlierFrac" : "OutlierFrac & MolSd & MolCov & ChiSq";
    if(FRAGCOV_FIX && Refine && !(!CovNorm && REPLACE_COV))
      fprintf(fp,"# StdDev & Coverage & %s refer to the interval between the current site and the next site\n",OutlierString);
    else
      fprintf(fp,"# StdDev & %s refer to the interval between the current site and the next site\n",OutlierString);
    if(!CovNorm && !REPLACE_COV){
      if(COVERAGE_FIX)
	fprintf(fp,"# Coverage includes any ends of molecules beyond the first or last aligned site (excluding endoutliers)\n");
      else
	fprintf(fp,"# Coverage does NOT include the ends of molecules beyond the first orlast aligned site\n");
    }
  } else
    fprintf(fp,"# StdDev refers to the interval between the current site and the next site\n");
  if(contig->HapSite[0] && contig->HapDelta[0])
    fprintf(fp,"# HapDelta, HapDeltaPhase and HapDeltaScore refer to the interval between the current site and the previous site\n");
  fprintf(fp,"#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence");
  if(ChimQuality && contig->sitecntFN[0] && !contig->HapSite[0]){
    fprintf(fp,"\tChimQuality");
    if(ChimQuality >= 2 && contig->sitecntN2[0] && contig->sitecntN3[0] && contig->sitecntN4[0] && contig->sitecntN5[0]){
      fprintf(fp,"\tSegDupL\tSegDupR\tFragileL\tFragileR\tOutlierFrac\tChimNorm");
    }
  } else if(!noQuality && !contig->HapSite[0])
    fprintf(fp,"\tQuality");
  if(contig->MaskL || contig->MaskR || (contig->SegDupCnt && SegDupMask))
    fprintf(fp,"\tMask");
  if(OutlierQuality >= 2 && contig->fragSd[0] && contig->expSd[0] && contig->fragCov[0] && contig->fragChiSq[0])
    fprintf(fp,"\tMolSd\tExpSd\tMolCov\tChiSq");
  if(contig->HapSite[0])
    fprintf(fp,"\tHapSite\tHapDelta\tSitePhaseConf\tDeltaPhaseConf\tHapSiteScore\tHapDeltaScore");
  if(mapSNR){
    if(CmapSNR)
      fprintf(fp,"\tGmeanSNR\tlnSNRsd\tSNRcount\tSNR ...");
    else
      fprintf(fp,"\tGmeanSNR\tlnSNRsd");
  }
  fprintf(fp,"\n");

  fprintf(fp,"#f int\tfloat\tint\tint\tint\tfloat\tfloat\tfloat\tfloat");
  if(ChimQuality && contig->sitecntFN[0] && !contig->HapSite[0]){
    fprintf(fp,"\tfloat");
    if(ChimQuality >= 2 && contig->sitecntN2[0] && contig->sitecntN3[0] && contig->sitecntN4[0] && contig->sitecntN5[0]){
      fprintf(fp,"\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat");
    }
  } else if(!noQuality && !contig->HapSite[0])
    fprintf(fp,"\tfloat");
  if(contig->MaskL || contig->MaskR || (contig->SegDupCnt && SegDupMask))
    fprintf(fp,"\tHex");
  if(OutlierQuality >= 2)
    fprintf(fp,"\tfloat\tfloat\tfloat\tfloat");
  if(contig->HapSite[0])
    fprintf(fp,"\tint\tfloat\tfloat\tfloat\tfloat\tfloat");
  if(mapSNR){
    if(CmapSNR)
      fprintf(fp,"\tfloat\tfloat\tint\tfloat ...");
    else
      fprintf(fp,"\tfloat\tfloat");
  }
  fprintf(fp,"\n");

#if 0 // multiple color support : not yet implemented 
  /* interleave colors/channels so sites are in increasing order */
  int newtotsites = 0;
  for(register int c = 0; c < colors; c++)
    newtotsites += contig->numsite[0];
  Cmerge *pn = new Cmerge[newtotsites+1];
  int cnt = 0;
  for(int c = 0; c < colors; c++){
    for(register int j = 1; j <= contig->numsite[c]; j++){
      pn[cnt].color = c;
      pn[cnt].newid = j;
      pn[cnt].site = contig->site[c][j];
      cnt++;
    }
  }
  assert(cnt == newtotsites);
  qsort(pn, newtotsites, sizeof(Cmerge), (intcmp*) siteinc);
  /* Add right end */
  pn[newtotsites].color = pn[newtotsites].newid = -1;
  pn[newtotsites].site = contig->site[0][contig->numsite[0]+1];
#endif

  if(BVERB){
    printf("output_draftIO:CovNorm=%d,REPLACE_COV=%d,left=%d(%0.3f),right=%d(%0.3f),n=%d(%0.3f),leftend=%0.5f,rightend=%0.5f\n",CovNorm,REPLACE_COV,left,Y[left],right,Y[right],n,Y[n],leftend,rightend);
    fflush(stdout);
  }

  double invLog10 = 1.0/log(10.0);
  double HapSiteOffset = -log(max(1e-300, HapSitePvalue)) * invLog10;
  double HapDeltaOffset = -log(max(1e-300, HapIndelPvalue)) * invLog10;

  int SiteID = 0;
  float *chimConf = contig->sitecntFN[0];
  float *SegDupL = contig->sitecntN2[0];
  float *SegDupR = contig->sitecntN3[0];
  float *FragileL = contig->sitecntN4[0];
  float *FragileR = contig->sitecntN5[0];

  float *OutlierQ = contig->sitecntN6[0];
  float *fragSd = contig->fragSd[0];
  float *expSd = contig->expSd[0];
  float *fragCov = contig->fragCov[0];
  float *fragChiSq = contig->fragChiSq[0];

  float *ChimNorm = contig->sitecntFNnorm[0];
  float *fragcnt = CovNorm ? contig->fragcntT[0] : (REPLACE_COV ? contig->sitecntFN[0] : COVERAGE_FIX ? contig->fragcntT[0] : contig->fragcnt[0]);
  float *CovScale = contig->fragcntTnorm[0];
  double lastY = 0.0;
  int prevk = 0;
  if(CMapID <= 0)
    CMapID = 1;

  double Ylen =  (Y[right]-Y[left]+leftend+rightend)*1000.0;

  if(VERB>=2){
    printf("contigid=%lld,llabs(contigid)=%lld,CMapID=%lld,Ylen=%0.4f\n",contigid,llabs(contigid),CMapID,Ylen*0.001);
    fflush(stdout);
  }

  for(int k = max(1,left); k <= min(n,right); k++){
    if(peak[k] > 0){
      int nextk = k+1;
      for(; nextk <= min(n,right); nextk++)
	if(peak[nextk])
	  break;

      if(DEBUG && !(Y[k] >= lastY && (nextk > min(n,right) || Y[nextk] >= Y[k])) && Y[k] <= Y[right]){
	printf("output_draftIO:k=%d,peak[k]=%d,Y[k]=%0.6f,peak[prevk=%d]=%d,Y[prevk]=%0.6f, lastY=%0.6f, peak[nextk=%d]=%d, Y[nextk]=%0.6f,n=%d,left=%d,right=%d,Y[left]=%0.6f,Y[right]=%0.6f\n", 
	       k,peak[k],Y[k],prevk,peak[prevk],Y[prevk], lastY,nextk,peak[nextk],Y[nextk],n,left,right,Y[left],Y[right]);
	if(contig->HapSite[0])
	  printf(" \t HapSite[prevk]=%d,HapSite[k]=%d,HapSite[nextk]=%d\n",contig->HapSite[0][prevk],contig->HapSite[0][k],contig->HapSite[0][nextk]);
	fflush(stdout);
	assert(Y[k] >= lastY);
	if(nextk <= min(n,right)) assert(Y[nextk] >= Y[k]);
	assert(Y[k] <= Y[right]);
      }
      lastY = Y[k];

      ++SiteID;

      double fragcntk = fragcnt[k];
      if(DEBUG && !REPLACE_COV && !(fragcntk >= 0.0)){
	printf("SiteID=%d:k=%d,Y[k]=%0.3f,peak[k]=%d,sitecnt[0][k]=%0.2f,fragcnt[0][k]=%0.2f,fragcntT[0][k]=%0.2f,fragcntk=%0.2f,CovNorm=%d,REPLACE_COV=%d,,fragcntTnorm[0][k]=%0.2f\n",
	       SiteID,k,Y[k],peak[k],contig->sitecnt[0][k],contig->fragcnt[0][k],contig->fragcntT[0][k],fragcntk,CovNorm,REPLACE_COV, CovNorm ? contig->fragcntTnorm[0][k] : 1.0);
	fflush(stdout);

	assert(fragcntk >= 0.0);
      }

      float cov = max(0.0,fragcntk);// NOTE : if -ReplaceCov is used, fragcntk is ChimQuality, which can be -1.0, which is NOT a valid coverage value
      float occ = COVFIX ? contig->sitecnt[0][k] : min(contig->sitecnt[0][k],cov);
      if(VERB>=2 || !isfinite(cov) || !isfinite(occ) || !(CovNorm==0 || isfinite(CovScale[k]))){
	printf("SiteID=%d:k=%d,Y[k]=%0.3f,peak[k]=%d,sitecnt[0][k]=%0.2f,fragcnt[0][k]=%0.2f,fragcntT[0][k]=%0.2f,fragcnt[k]=%0.2f,CovNorm=%d,REPLACE_COV=%d,cov=%0.2f,occ=%0.2f,fragcntTnorm[0][k]=%0.2f\n",
	       SiteID,k,Y[k],peak[k],contig->sitecnt[0][k],contig->fragcnt[0][k],contig->fragcntT[0][k],fragcntk,CovNorm,REPLACE_COV,cov,occ,CovNorm ? contig->fragcntTnorm[0][k] : 1.0);
	fflush(stdout);
      }
      if(CovNorm){
	float scale = CovScale[k];
	cov *= scale;
	occ *= scale;
      }
      if(DEBUG) assert(isfinite(cov));
      if(DEBUG) assert(isfinite(occ));
      // HERE HERE : update how StdDev is computed, so it corresponds to actual weighted molecule interval SD 
      double SDcov = min(16.0f,max(1.0f,cov));
      double StdDev = 0.0;
      if(nextk <= min(n,right) && peak[nextk])
	StdDev = 1000.0*sqrt((SF[0]*SF[0]+SD[0]*SD[0]*(Y[nextk]-Y[k]))/SDcov)/(1.0-FN[0]);// HERE HERE : apply QUADRATIC_VARIANCE
      if(DEBUG && !(isfinite(StdDev))){
	printf("SiteID=%d:k=%d,n=%d,peak[k]=%d,sitecnt[0][k]=%0.1f,fragcnt[0][k]=%0.1f,CovNorm=%d,REPLACE_COV=%d,cov=%0.1f,occ=%0.1f,fragcntTnorm[0][k]=%0.1f\n",
	       SiteID,k,n,peak[k],contig->sitecnt[0][k],contig->fragcnt[0][k],CovNorm,REPLACE_COV,cov,occ,CovNorm ? contig->fragcntTnorm[0][k] : 1.0);
	printf("SDcov=%0.6e,StdDev=%0.6e,FN[0]=%0.6f,SF[0]=%0.3f,SD[0]=%0.6f,SR[0]=%0.6f:Y[nextk=%d]=%0.3f,Y[k]=%0.3f,Y[n]=%0.3f,Y[n+1]=%0.3f,peak[nextk]=%d,left=%d,right=%d\n",
	       SDcov,StdDev,FN[0],SF[0],SD[0],SR[0],nextk,Y[nextk],Y[k],Y[n],Y[n+1],peak[nextk],left,right);
	for(int t = k; t <= nextk; t++){
	  float cov = fragcnt[t];
	  float occ = COVFIX ? contig->sitecnt[0][t] : min(contig->sitecnt[0][t],cov);
	  printf("t=%d,peak[t]=%d,Y[t]=%0.3f:sitecnt[0][t]=%0.1f,fragcnt[0][t]=%0.1f,CovNorm=%d,REPLACE_COV=%d,cov=%0.1f,occ=%0.1f,fragcntTnorm[0][t]=%0.1f\n",
		 t,peak[t],Y[t],contig->sitecnt[0][t],contig->fragcnt[0][t],CovNorm,REPLACE_COV,cov,occ,CovNorm ? contig->fragcntTnorm[0][t] : 1.0);
	}
	fflush(stdout);
	assert(isfinite(StdDev));
      }
      fprintf(fp,"%lld\t%0.1f\t%d\t%d\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.1f",
	      (contigid > 0) ? llabs(contigid) : CMapID, // CMapID
	      Ylen, // ContigLength
	      NumSites, // NumSites
	      SiteID, // SiteID
	      1, // LabelChannel
	      (Y[k]-Y[left]+leftend)*1000.0, // Position
	      StdDev,// StdDev
	      cov, // Coverage
	      occ  // Occurrence
	      );
      if(ChimQuality && chimConf && !contig->HapSite[0]){
	fprintf(fp,"\t%0.2f", chimConf[k]);// ChimQuality

	if(DEBUG && !(isfinite(chimConf[k]) && -1.0 <= chimConf[k] && chimConf[k] <= 100.001f)){
	  printf("k=%d/%d:chimConf[k]= %0.3f\n",k,n,chimConf[k]);
	  fflush(stdout);
	  assert(isfinite(chimConf[k]) && -1.0 <= chimConf[k] && chimConf[k] <= 100.001f);
	}

	if(ChimQuality >= 2 && SegDupL && SegDupR && FragileL && FragileR && OutlierQ && ChimNorm){
	  fprintf(fp,"\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f",SegDupL[k], SegDupR[k], FragileL[k], FragileR[k], OutlierQ[k], ChimNorm[k]);
	  
	  if(DEBUG && !(isfinite(SegDupL[k]) && -1.0 <= SegDupL[k] && SegDupL[k] <= 100.001f)){
	    printf("WARNING:k=%d(left=%d,right=%d,n=%d):SegDupL[k]= %0.6e (%0.6e)\n",k,left,right,n,SegDupL[k], contig->sitecntN2[0][k]);
	    fflush(stdout);fflush(fp);
	    assert(isfinite(SegDupL[k]) && -1.0 <= SegDupL[k] && SegDupL[k] <= 100.001f);
	  }
	  if(DEBUG && !(isfinite(SegDupR[k]) && -1.0 <= SegDupR[k] && SegDupR[k] <= 100.001f)){
	    printf("WARNING:k=%d(left=%d,right=%d,n=%d):SegDupR[k]= %0.6e (%0.6e)\n",k,left,right,n,SegDupR[k], contig->sitecntN3[0][k]);
	    fflush(stdout);fflush(fp);
	    assert(isfinite(SegDupR[k]) && -1.0 <= SegDupR[k] && SegDupR[k] <= 100.001f);
	  }
	  if(DEBUG && !(isfinite(FragileL[k]) && -1.0 <= FragileL[k])){
	    printf("WARNING:k=%d(left=%d,right=%d,n=%d):FragileL[k]= %0.6e (%0.6e)\n",k,left,right,n,FragileL[k], contig->sitecntN4[0][k]);
	    fflush(stdout);fflush(fp);
	    assert(isfinite(FragileL[k]) && -1.0 <= FragileL[k]);
	  }
	  if(DEBUG && !(isfinite(FragileR[k]) && -1.0 <= FragileR[k])){
	    printf("WARNING:k=%d(left=%d,right=%d,n=%d):FragileR[k]= %0.6e (%0.6e)\n",k,left,right,n,FragileR[k], contig->sitecntN5[0][k]);
	    fflush(stdout);fflush(fp);
	    assert(isfinite(FragileR[k]) && -1.0 <= FragileR[k]);
	  }
	  if(DEBUG && !(isfinite(OutlierQ[k]) && 0.0 <= OutlierQ[k] && OutlierQ[k] <= 100.001f)){
	    printf("WARNING:k=%d(left=%d,right=%d,n=%d):OutlierQ[k]= %0.6e (%0.6e)\n",k,left,right,n,OutlierQ[k], contig->sitecntN6[0][k]);
	    fflush(stdout);fflush(fp);
	    assert(isfinite(OutlierQ[k]) && 0.0 <= OutlierQ[k] && OutlierQ[k] <= 100.001f);
	  }
	  if(DEBUG && !(isfinite(ChimNorm[k]) && -1.0 <= ChimNorm[k])){
	    printf("WARNING:k=%d(left=%d,right=%d,n=%d):ChimNorm[k]= %0.6e (%0.6e)\n",k,left,right,n,ChimNorm[k], contig->sitecntFNnorm[0][k]);
	    fflush(stdout);fflush(fp);
	    assert(isfinite(ChimNorm[k]) && -1.0 <= ChimNorm[k]);
	  }
	}
      } else if(!noQuality && !contig->HapSite[0]){
	if(DEBUG) assert(fabs(avgmolpvden[SiteID])>0.0); //protect against nans in output files
	fprintf(fp, "\t%.2f", avgmolpv[SiteID]/avgmolpvden[SiteID]);
	totavg += avgmolpv[SiteID]/avgmolpvden[SiteID]; //total for right end
      }
      if(contig->MaskL || contig->MaskR || (contig->SegDupCnt && SegDupMask))
	fprintf(fp,"\t%lx", (SiteID==1) ? contig->MaskL | Mask[k] : Mask[k]);
      if(OutlierQuality >= 2){
	fprintf(fp,"\t%0.1f\t%0.1f\t%0.2f\t%0.2e", fragSd[k]*1000.0, expSd[k]*1000.0, fragCov[k], fragChiSq[k]);

	if(DEBUG && !(isfinite(fragSd[k]) && 0.0 <= fragSd[k])){
	  printf("WARNING:k=%d(left=%d,right=%d,n=%d):fragSd[k]= %0.6e (%0.6e)\n",k,left,right,n,fragSd[k], contig->fragSd[0][k]);
	  fflush(stdout);fflush(fp);
	  assert(isfinite(fragSd[k]) && 0.0 <= fragSd[k]);
	}
	if(DEBUG && !(isfinite(expSd[k]) && 0.0 <= expSd[k])){
	  printf("WARNING:k=%d(left=%d,right=%d,n=%d):expSd[k]= %0.6e (%0.6e)\n",k,left,right,n,expSd[k], contig->expSd[0][k]);
	  fflush(stdout);fflush(fp);
	  assert(isfinite(expSd[k]) && 0.0 <= expSd[k]);
	}
	if(DEBUG && !(isfinite(fragCov[k]) && 0.0 <= fragCov[k])){
	  printf("WARNING:k=%d(left=%d,right=%d,n=%d):fragCov[k]= %0.6e (%0.6e)\n",k,left,right,n,fragCov[k], contig->fragCov[0][k]);
	  fflush(stdout);fflush(fp);
	  assert(isfinite(fragCov[k]) && 0.0 <= fragCov[k]);
	}
	if(DEBUG && !(isfinite(fragChiSq[k]) && 0.0 <= fragChiSq[k] && fragChiSq[k] <= 1.0)){
	  printf("WARNING:k=%d(left=%d,right=%d,n=%d):fragChiSq[k]= %0.6e (%0.6e)\n",k,left,right,n,fragChiSq[k], contig->fragChiSq[0][k]);
	  fflush(stdout);fflush(fp);
	  assert(isfinite(fragChiSq[k]) && 0.0 <= fragChiSq[k] && fragChiSq[k] <= 1.0);
	}
      }
      if(contig->HapSite[0]){
	int *HapSite = contig->HapSite[0];
	double HapDeltaScore = contig->HapDeltaScore[0][prevk];
	double delta = contig->HapDelta[0][prevk];
	for(int t = prevk+1; t < k; t++){
	  delta += contig->HapDelta[0][t];
	  HapDeltaScore = max(HapDeltaScore,contig->HapDeltaScore[0][t]);
	}
	HapDeltaScore = (fabs(delta) > 0.001) ? HapDeltaOffset + HapDeltaScore * invLog10 : 0.0;

	double HapSiteScore = HapSiteOffset + contig->HapSiteScore[0][k] * invLog10;
	double HapSitePhase = (1 <= HapSite[k] && HapSite[k] <= 2) ? contig->HapSitePhase[0][k]*invLog10 : -1.0;
	double HapDeltaPhase = (fabs(delta) > 0.001) ? contig->HapDeltaPhase[0][prevk] * invLog10 : -1.0;
	fprintf(fp,"\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.2f\t%0.2f", HapSite[k], delta * 1000.0, HapSitePhase, HapDeltaPhase, HapSiteScore, HapDeltaScore);
      }
      if(mapSNR){
	double wt = max(0.1,SNRwtsum[0][k]);
	if(DEBUG) assert(isfinite(wt));

	double sum = SNRsum[0][k];
	if(DEBUG && !isfinite(sum)){
	  printf("ERROR: non-finite number: sum=%g from SNRsum[0][k=%d] (n=%d,left=%d,right=%d)\n", sum, k,n,left,right);
	  fflush(stdout);
	  assert(isfinite(sum));
	}
	if(DEBUG && !isfinite(exp(sum/wt))){
	  printf("SNRsum[0][k=%d]= sum = %0.8e, wt=%0.8e, exp(sum/wt)= %g\n",k, SNRsum[0][k], wt, exp(sum/wt));
	  fflush(stdout);
	  assert(isfinite(exp(sum)));
	}

	double sumsq = SNRsumsq[0][k];
	sum /= wt;
	sumsq /= wt;
	double sd = sqrt(max(0.0, (sumsq - sum*sum)));
	if(DEBUG && !isfinite(sd)){
	  printf("ERROR: non-finite number: sd=%g from sumsq=%.15g sum=%.15g wt=%.15g\n", sd, sumsq, sum, wt);
	  fflush(stdout);
	  assert(isfinite(sd));
	}

	fprintf(fp,"\t%0.4f\t%0.4f",
		exp(sum),                                   // GmeanSNR
		sd                                          // lnSNRsd
		);
	if(CmapSNR){

	  int cnt = SNRcnt[0][k];
	  double *pSNR = SNR[0][k];

	  qsort(pSNR,cnt,sizeof(double),(intcmp*)doubleInc);
	    
	  fprintf(fp,"\t%d",cnt);
	  for(int u = 0; u < cnt; u++)
	    fprintf(fp,"\t%0.4f",pSNR[u]);
	}
      }
      fprintf(fp,"\n");

      prevk = k;
    }
  }

  /* right end */
  int kR = min(n,right);
  double fragcntkR = fragcnt[kR];
  float cov = max(0.0,fragcntkR);
  double SDcov = max(1.0f,min(16.0f,contig->fragcnt[0][kR]));
  totavg /= NumSites; //now average
  if(VERB>=2){
    printf("SDcov=%0.6e:n=%d,right=%d,kR=%d,contig->fragcnt[kR]=%0.1f,contig->fragcnt[kR-1]=%0.1f,cov=%0.1f\n",SDcov,n,right,kR,contig->fragcnt[0][kR],contig->fragcnt[0][kR-1],cov);
    fflush(stdout);
  }
  fprintf(fp,"%lld\t%0.1f\t%d\t%d\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.1f",
	  (contigid > 0) ? llabs(contigid) : CMapID, // CMapID
	  Ylen, // NEW : ContigLength
	  NumSites, // NumSites
	  NumSites + 1, // SiteID = right end
	  0, // LabelChannel
	  Ylen, // NEW : Position (same as ContigLength)
	  0.0 /* 1000.0*sqrt((SF[0]*SF[0]+SD[0]*SD[0]*(Y[kR+1]-Y[kR]))/SDcov)/(1.0-FN[0])*/, // StdDev
	  1.0 /* cov */, // Coverage
	  0.0/* WAS 1.0 */   // Occurrence
	  );
  if(ChimQuality && chimConf && !contig->HapSite[0]){
    fprintf(fp,"\t%0.2f", END_CHIMQUAL);// ChimQuality
    if(ChimQuality >= 2 && SegDupL && SegDupR && FragileL && FragileR && OutlierQ && ChimNorm){
      fprintf(fp,"\t%0.2f\t%0.2f\t0.00\t0.00\t0.00\t%0.2f",END_CHIMQUAL,END_CHIMQUAL,END_CHIMQUAL);
    }
  } else if(!noQuality && !contig->HapSite[0])
    fprintf(fp,"\t%0.2f", totavg);
  if(contig->MaskL || contig->MaskR || (contig->SegDupCnt && SegDupMask))
    fprintf(fp,"\t%lx", contig->MaskR | Mask[n+1]);
  if(OutlierQuality >= 2)
    fprintf(fp,"\t%0.1f\t%0.1f\t%0.2f\t%0.2e",0.0,0.0,0.0,1.0);
  if(contig->HapSite[0]){
    double HapDeltaScore = contig->HapDeltaScore[0][prevk];
    double delta = contig->HapDelta[0][prevk];
    for(int t = prevk+1; t <= n; t++){
      delta += contig->HapDelta[0][t];
      HapDeltaScore = max(HapDeltaScore,contig->HapDeltaScore[0][t]);
    }
    HapDeltaScore = (fabs(delta) > 0.001) ? HapDeltaOffset + HapDeltaScore * invLog10 : 0.0;
    double HapDeltaPhase = (fabs(delta) > 0.001) ? contig->HapDeltaPhase[0][prevk] * invLog10 : -1.0;

    fprintf(fp,"\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.2f\t%0.2f",0, delta * 1000.0, -1.0, HapDeltaPhase, -1.0, HapDeltaScore);
  }
  if(mapSNR){
    fprintf(fp,"\t%0.3f\t%0.3f",
	    0.0, // SNRgmean
	    0.0  // lnSNRsd
	    );
    if(CmapSNR)
      fprintf(fp,"\t0");
  }
  fprintf(fp,"\n");

  FILEclose(fp);

  if(mapSNR){/* free up memory */
    for(int c = 0; c < colors; c++){
      delete [] SNRwtsum[c];
      delete [] SNRsum[c];
      delete [] SNRsumsq[c];
      if(CmapSNR){
	delete [] SNRcnt[c];
	delete [] SNRmax[c];
	for(int t = 1; t <= nn[c]; t++)
	  delete [] SNR[c][t];
	delete [] SNR[c];
      }
    }
  }

  if(!Refine || !contigid /* WAS !(CMapID > 0 && contigid < 0) */){
    printf("Skipping output of refined .xmap:CMapID=%lld,contigid=%lld\n",CMapID,contigid);
    fflush(stdout);

    if(1 /* NEW contig->HapSite[0] */){
      errno = 0;
      if(rename(tmpfilename,origfilename)){
	int eno = errno;
	char *err = strerror(eno);
	printf("failed to rename %s to %s:errno=%d:%s\n",tmpfilename,origfilename,eno,err);
	fflush(stdout);  exit(1);
      }
    }

    CMapID = origCMapID;
    free(write_buffer);

    if(VERB){
      printf("Completed output of %s: cum wall time= %0.6f\n",origfilename,wtime());
      fflush(stdout);
    }

    free(tmpfilename);
    free(origfilename);
    
    delete [] Mask;

    return;
  }

  if(contigid > 0)
    CMapID = contigid;

  int split_maps = (MappedUnsplit >= 0) ? !MappedUnsplit : (SplitRef && CMapID > 0) ? 1 : 0;
  if(!split_maps){
    printf("WARNING:Cannot output refined .xmap into a single file for all contigs: input each -ref map from a seperate file OR use -mapped-unsplit 0 to generate separate .xmap's for each contig\n");
    fflush(stdout);
    CMapID = origCMapID;
    free(write_buffer);

    if(1 /* NEW contig->HapSite[0] */){
      errno = 0;
      if(rename(tmpfilename,origfilename)){
	int eno = errno;
	char *err = strerror(eno);
	printf("failed to rename %s to %s:errno=%d:%s\n",tmpfilename,origfilename,eno,err);
	fflush(stdout);  exit(1);
      }
    }

    if(VERB){
      printf("Completed output of %s: cum wall time= %0.6f\n",origfilename,wtime());
      fflush(stdout);
    }

    free(tmpfilename);
    free(origfilename);

    delete [] Mask;

    return;
  }

  Cmap **goodmaps = new Cmap *[contig->nummaps];
  int *goodflip = new int[contig->nummaps];
  int numgoodmaps = 0;
  char cmapprefix[PATH_MAX];

  /* use array to order alignments as in output_xmap() */
  if(VERB>=2 && CMapID == 19842){
    printf("%s_refined_contig%lld_q.cmap: maps to be output:\n",basename,CMapID);
    fflush(stdout);
  }

  Cmaporder *maporder = new Cmaporder[contig->nummaps];
  for(int m = 0; m < contig->nummaps; m++){
    maporder[m].m = m;
    maporder[m].left = Y[n+1];

    Ccontig *mcontig = &contig->contig[m];
    int mapid = mcontig->mapid;
    Cmap *pmap = Gmap[mapid];
    double mapwt = 1.0;

    if(contig->mapWT && (mapwt = contig->mapWT[m]) >= RefineMinMapwt){
      if(VERB>=2 && CMapID == 19842){
	printf("m=%d: id=%lld, mapWT=%0.4f\n",m,pmap->id,mapwt);
	fflush(stdout);
      }
      if(1){/* duplicate maps are possible : remove them from _q.cmap */
	int t = 0;
	for(; t < numgoodmaps; t++)
	  if(goodmaps[t]->id == pmap->id)
	    break;
	if(t >= numgoodmaps){
	  goodflip[numgoodmaps] = contig->flip[m];
	  goodmaps[numgoodmaps++] = pmap;
	} else if(!(MultiMatches || (CMapID > 0 && contig->HapSite[0]) || goodflip[t] != contig->flip[m])){
	  printf("WARNING:Duplicate map ids in %s_refined_contig%lld_q.cmap : goodmaps[%d]:mapid=%d,id=%lld,flip=%d; goodmaps[%d]:mapid=%d,id=%lld,flip=%d\n",
		 basename,CMapID,t,goodmaps[t]->mapid,goodmaps[t]->id,goodflip[t],m,pmap->mapid,pmap->id,contig->flip[m]);
	  printf("\t contig->nummaps= %d\n",contig->nummaps);
	  for(int m = 0; m < contig->nummaps; m++)
	    printf("m=%d/%d: mapid=%d, id=%lld, flip=%d, mapWT= %0.6f, logPV= %0.2f\n",
		   m,contig->nummaps,contig->contig[m].mapid,gmap[contig->contig[m].mapid]->id,contig->flip[m],contig->mapWT[m],contig->logPV[m]);
	  fflush(stdout);
	  assert(goodmaps[t]->id != pmap->id);
	}
      }
    }

    for(int c = 0; c < colors; c++){
      int M = mcontig->numsite[c];
      int *map = contig->sitemap[c][m];/* Gmap[j=0..M+1] maps site j of color c of map m to one of contig->numsite[c] sites : 0 refers to contig's left end, contig->numsite[c]+1 to the right end */
      int *mapL = contig->sitemapL[c] ? contig->sitemapL[c][m] : map;
      for(int J = 1; J <= M; J++){
	int R = map[J];
	if(R < 0)
	  continue;
	if(DEBUG) assert(0 < R && R <= contig->numsite[c]);
	int L = mapL[J];
	if(DEBUG/* HERE >=2 */) assert(0 < L && L <= R);
	if(Y[L] < maporder[m].left)
	  maporder[m].left = Y[L];
	break;
      }
    }
  }

  if(DEBUG>=2){
    for(int m = 1; m < numgoodmaps; m++){
      Cmap *pmap = goodmaps[m];
      for(int t = 0; t < m; t++){
	if(goodmaps[t]->id == pmap->id){
	  printf("WARNING:Duplicate map ids in %s_refined_contig%lld_q.cmap : goodmaps[%d]:mapid=%d,id=%lld, goodmaps[%d]:mapid=%d,id=%lld\n",
		 basename,CMapID,t,goodmaps[t]->mapid,goodmaps[t]->id,m,pmap->mapid,pmap->id);
	  fflush(stdout);
	  assert(goodmaps[t]->id != pmap->id);
	}
      }
    }
  }

  // HERE HERE : wait until xmap has been output so molecules in _q.cmap can be limited to those in the .xmap 

  /* output <basename>_refined_contig<id>_q.cmap */
  sprintf(cmapprefix,"%s_refined_contig%lld_q",basename, CMapID);

  if(usecolor)
    colors = 2;
  if(usecolor==2){
    for(int i = 0; i < numgoodmaps; i++)
      goodmaps[i]->colorswap(usecolor);
    char *tmp = Nickase[0]; Nickase[0] = Nickase[1]; Nickase[1] = tmp;
  }

  if(CmapMergeBnx)
    output_bnx(cmapprefix,goodmaps,0,numgoodmaps-1, 0, NULL, NULL, -1);
  else
    output_cmap(cmapprefix,goodmaps,0,numgoodmaps-1);

  if(usecolor==2){
    for(int i = 0; i < numgoodmaps; i++)
      goodmaps[i]->colorswap(usecolor);
    char *tmp = Nickase[0]; Nickase[0] = Nickase[1]; Nickase[1] = tmp;
  }
  if(usecolor)
    colors = 1;

  /* output <basename>_refined_contig%lld.xmap */
  sprintf(filename,"%s_refined_contig%lld.xmap",basename, CMapID);

  if(!checkFile(filename)) {
    if(VERB){
      printf("Generating %s\n",filename);
      fflush(stdout);
    }
    if((fp = fopen(filename,"w"))==NULL){
      printf("failed to open file %s for writing xmap file\n",filename);
      fflush(stdout);exit(1);
    }

    setbuffer(fp, write_buffer, blen);

    /* write out commandline */
    printversion(fp);

    /* Also output indel SVs */
    FILE *fpSV = NULL;
    char filenameSV[PATH_MAX];
    char *write_bufferSV = NULL;
    if(INDEL){
      sprintf(filenameSV,"%s_refined_contig%lld.indel",basename, CMapID);
      if(!checkFile(filenameSV)){
	if(VERB){
	  printf("Generating %s\n",filenameSV);
	  fflush(stdout);
	}
	if((fpSV = fopen(filenameSV,"w"))==NULL){
	  int eno = errno;
	  char *err = strerror(errno);
	  printf("failed to open file %s for writing indel file:errno=%d:%s\n",filenameSV,eno,err);
	  fflush(stdout);exit(1);
	}
	size_t blen = 10*1024*1024;/* matches scif_io buffer size */
	write_bufferSV = (char *)malloc(blen);
	if(!write_bufferSV){
	  printf("output_draftIO:malloc(%llu) failed\n",(unsigned long long)blen);
	  fflush(stdout);exit(1);
	}
	setbuffer(fpSV, write_bufferSV, blen);

	printversion(fpSV);

	fprintf(fpSV,"# SMAP File Version:\t0.0\n"); //  (0.0 was indels only; 0.1 has newer types)
	if(contig->HapSite[0])
	  fprintf(fpSV,"# Reference Maps From:\t%s_contig%lld.hmap\n", basename, CMapID);/* output in this function above */
	else
	  fprintf(fpSV,"# Reference Maps From:\t%s_contig%lld.cmap\n", basename, CMapID);/* output in this function above */
	fprintf(fpSV,"# Query Maps From:\t%s.%s\n",cmapprefix, CmapMergeBnx ? "bnx" : "cmap");/* output in this function above */
	fprintf(fpSV,"# Xmap Entries From:\t%s\n",filename);    
	fprintf(fpSV,"#h SmapEntryID\tQryContigID\tRefcontigID1\tRefcontigID2\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tType\tXmapID1\tXmapID2\tMapWt\n");
	fprintf(fpSV,"#f int        \tint        \tint         \tint         \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat \tstring \tint\tint\tfloat\n");
      }
    }

    if(XmapChim || XmapUnique){
      printf("WARNING: -xmapchim and -xmapUnique ignored for output %s\n",filename);
      fflush(stdout);
    }
  
    /* sort alignments in order of id1, sites1[0] (null alignments last) */
    qsort(maporder, contig->nummaps,sizeof(Cmaporder), (intcmp*)MaporderSiteInc);

    if(DEBUG && contigid <= 0) assert(spots_filename[0]);

    fprintf(fp,"# XMAP File Version:\t0.2\n");
    fprintf(fp,"# Label Channels:\t%d\n", colors);

    if(contig->HapSite[0])
      fprintf(fp,"# Reference Maps From:\t%s_contig%lld.hmap\n", basename, CMapID);/* output in this function above */
    else
      fprintf(fp,"# Reference Maps From:\t%s_contig%lld.cmap\n", basename, CMapID);/* output in this function above */
    fprintf(fp,"# Query Maps From:\t%s.%s\n",cmapprefix, CmapMergeBnx ? "bnx" : "cmap");/* output in this function above */
    if(contig->HapSite[0]){
      fprintf(fp,"#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tAllele\tQryLen\tRefLen\tLabelChannel\tAlignment\tMapWt\tMaxOutlier\n");
      fprintf(fp,"#f int        \tint        \tint        \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat     \tstring \tint   \tfloat \tfloat \tint         \tstring   \tfloat\tfloat\n");
    } else {
      fprintf(fp,"#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment\tMapWt\tMaxOutlierKb\n");
      fprintf(fp,"#f int        \tint        \tint        \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat     \tstring \tfloat \tfloat \tint         \tstring   \tfloat\tfloat\n");
    }

    int aligncnt = 0;
    int SmapCnt = 0;

    double refLen = (Y[right] - Y[left] +leftend + rightend) * 1000.0;            /* RefLen */
    /* create mapping for Y[I] index to consensus site count (the sites with peak[I] > 0 in the range peak[left .. right) */
    int *refsite = new int[n+2];
    int refcnt = 0;
    for(int I = 0; I < left; I++)
      refsite[I] = -1;
    for(int I = left; I <= right; I++){
      if(peak[I] > 0)
	refsite[I] = ++refcnt;
      else
	refsite[I] = -1;
    }
    for(int I = right+1; I <= n+1; I++)
      refsite[I] = -1;

    if(VERB>=2){
      printf("created mapping from Hcuts[] index to consensus index : n=%d,N=refcnt=%d\n",n,refcnt);
      for(int I = left; I <= right; I++)
	if(peak[I]){
	  if(contig->HapSite[0])
	    printf("i=%d:Y[i]=%0.4f:refsite[i]=%d,contig->Hdel[0][i]=%d,peak[i]=%d,contig->HapSite[0][i]=%d\n",I,Y[I],refsite[I],contig->Hdel[0][I],peak[I],contig->HapSite[0][I]);
	  else
	    printf("i=%d:Y[i]=%0.4f:refsite[i]=%d,contig->Hdel[0][i]=%d,peak[i]=%d\n",I,Y[I],refsite[I],contig->Hdel[0][I],peak[I]);
	}
      fflush(stdout);
    }

    /* locate leftmost and rightmost real site with peak[] != 0*/
    int leftpeak = left, rightpeak = right;
    for(; leftpeak < rightpeak; leftpeak++)
      if(peak[leftpeak])
	break;
    for(; leftpeak < rightpeak; rightpeak--)
      if(peak[rightpeak])
	break;

    double Ilog10 = 1.0/log(10.0);

    for(int t = 0; t < contig->nummaps; t++){
      int m = maporder[t].m;
      Ccontig *mcontig = &contig->contig[m];
      int flip = contig->flip[m];
      int mapid = mcontig->mapid;
      double mapwt = 1.0;
      if(VERB>=2 && Gmap[mapid]->id == 158077LL){
	printf("t=%d/%d:mapid=%d,id=%lld,mapWT=%0.6f,RefineMinMapwt= %0.6f\n",
	       t,contig->nummaps,mapid,Gmap[mapid]->id,contig->mapWT ? contig->mapWT[m] : 1.0, RefineMinMapwt);
	fflush(stdout);
      }

      if(contig->mapWT && (mapwt = contig->mapWT[m]) < RefineMinMapwt)
	continue;
    
      aligncnt++;

      Cmap *Xmap = Gmap[mapid];
      int M = mcontig->numsite[0];
      int Xleft = mcontig->trimL[0];
      int Xright = mcontig->trimR[0];
      int *map = contig->sitemap[0][m];
      int *mapL = contig->sitemapL[0] ? contig->sitemapL[0][m] : map;
      int *outlier = contig->outlier[0][m];
      FLOAT *X = Xmap->site[0];
      //      double **XX = contig->X[0];
    
      /* locate leftmost and rightmost aligned sites */
      int LI = -1, LJ = -1, RI = -1, RJ = -1, LK= -1;// RK = -1;
      for(int J = 1; J <= M; J++){
	int R = map[J];
	if(R < 0 || R < left || R > right)
	  continue;
	LI = R;
	LK = R - mapL[J];
	LJ = J;
	break;
      }
      if(LI < 0)/* molecule m does not overlap reference in interval Y[left..right] */
	continue;

      for(int J = M; J >= 1; J--){
	int R = map[J];
	if(R < 0 || R < left || R > right)
	  continue;
	RI = R;
	//	RK = R - mapL[J];
	RJ = J;
	break;
      }
      if(DEBUG) assert(RI >= left);

      int qrystartidx = flip ? Xright - LJ : Xleft + LJ;
      int qryendidx = flip ? Xright - RJ : Xleft + RJ;
      double qrystartpos = X[qrystartidx]*1000.0; /* QryStartPos */
      double qryendpos = X[qryendidx]*1000.0;   /* QryEndPos */
      double qryLen = X[M+1]*1000.0;              /* QryLen */
      double refstartpos = (Y[LI] - Y[left] + leftend) * 1000.0;          /* RefStartPos */
      double refendpos    = (Y[RI] - Y[left] + leftend) * 1000.0;         /* RefEndPos */

      if(VERB>=2 && Xmap->id == 158077LL){
	printf("m=%d,Xmap->id=%lld:LJ=%d,RJ=%d,M=%d\n",m, Xmap->id,LJ,RJ,M);
	fflush(stdout);
      }
      if(DEBUG && !(refendpos >= refstartpos - 1e-6)){
	printf("LI=%d,LJ=%d,RI=%d,RJ=%d,n=%d,M=%d,left=%d,right=%d\n",LI,LJ,RI,RJ,n,M,left,right);
	printf("qrystartpos= %0.4f, qryendpos= %0.4f, refstartpos= %0.4f, refendpos= %0.1f,Y[LI]= %0.7f, Y[RI]= %0.7f\n",
	       qrystartpos,qryendpos, refstartpos, refendpos, Y[LI], Y[RI]);
	fflush(stdout);
	assert(refendpos >= refstartpos - 1e-6);
      }

      fprintf(fp,"%llu\t%lld\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%c\t%0.2f\t",
	      (unsigned long long)aligncnt, /* XmapEntryID */
	      Xmap->id, /* QryContigID */
	      CMapID /* -contigid */, /* RefContigID */
	      qrystartpos,
	      qryendpos,
	      refstartpos,
	      refendpos,
	      flip ? '-' : '+', /* Orientation */
	      contig->logPV[m]  /* Confidence */
	      );

      if(INDEL){/* output internal outliers of current alignment as indels */
	int lastI = LI;
	int lastJ = LJ;
	int lastK = LK;
	int I = -1,K,J;

	for(J = LJ+1; J <= RJ; lastI = I, lastK = K, lastJ = J, J++){
	  while(J <= RJ && (I = map[J]) < 0)
	    J++;
	  if(DEBUG) assert(J <= RJ);
	  K = mapL[J];

	  if(VERB>=2 && Xmap->id == 158077LL){
	    printf("Between J=%d..%d and I=%d..%d (M=%d,left=%d,right=%d,Y[left,right]=%0.3f,%0.3f):outlier=%d\n",lastJ,J,lastI,K,M,left,right,Y[left],Y[right],outlier[lastJ]);
	    fflush(stdout);
	  }
	  if(!outlier[lastJ])/* not an outlier at left end of current interval lastJ .. J */
	    continue;

	  /* HERE : merge nearby outliers into a single SV */

	  /* outlier lies between sites lastI .. K on reference and lastJ .. J on query */
	  int qrystartidx = flip ? Xright - lastJ : Xleft + lastJ;
	  int qryendidx = flip ? Xright - J : Xleft + J;
	  double qrystartpos = X[qrystartidx] * 1000.0; /* QryStartPos */
	  double qryendpos = X[qryendidx] * 1000.0; /* QryEndPos */
	  double refstartpos = (Y[lastI] - Y[left] + leftend) * 1000.0;/* RefStartPos */
	  double refendpos = (Y[K] - Y[left] + leftend) * 1000.0;/* RefEndPos */

	  if(DEBUG/* HERE >=2 */) assert(refstartpos <= refendpos);
	  if(DEBUG/* HERE >=2 */ && !(flip ? (qryendpos <= qrystartpos) : (qrystartpos <= qryendpos))){
	    printf("m=%d:id=%lld:J=%d..%d,I=%d..%d,outlier=%d:flip=%d,Xleft=%d,Xright=%d,qrystartidx=%d,qryendidx=%d,M=%d,qrystartpos=%0.3f,qryendpos=%0.3f\n",
		   m,Xmap->id,lastJ,J,lastI,K,outlier[lastJ],flip,Xleft,Xright,qrystartidx,qryendidx,M,qrystartpos,qryendpos);
	    fflush(stdout);
	    assert(flip ? (qryendpos <= qrystartpos) : (qrystartpos <= qryendpos));
	  }
	  if(DEBUG/* HERE >=2 */) assert(J > lastJ);
	  if(DEBUG/* HERE >=2 */) assert(K > lastI);

	  /* HERE : xmapUnique not yet implemented */

	  // compute outlier confidence value as in output_xmap.cpp, but ignore left and right end confidence
	  // declared in RGentigRefScore.h
	  extern void SintDetailRA(double X, double Y, int m, int n, int J, int I, int K /* right */, int T /* left */, FLOAT *Yref, double &Bias, double &Pen, double &Gauss, double &PenSm, int verb);
	  
	  double PenSm = 0.0, Bias,Pen,Gauss;
	  SintDetailRA(fabs(X[qryendidx] - X[qrystartidx]), Yc(Y, I, I-K) - Yc(Y, lastI, lastI-lastK), J - lastJ, K-lastI, J, I, I-K, lastI-lastK, Y, Bias,Pen,Gauss,PenSm, 0) ;
	  double outscore = Bias + Pen + Gauss; // NOTE : excludes PenSm which is the misresolved penalty for the right end of the interval
	  double centerConfidence = -outscore * Ilog10;
	  double outlierConfidence = min(centerConfidence, contig->logPV[m]);

	  if(VERB>=2 && CMapID==9992 && Xmap->id == 355200){
	    printf("\t I=%d,K=%d,J=%d,lasti=%d,lastk=%d,lastJ=%d:Yc(Y,I,K)=%0.4f,Yc(Y,lastI,lastK)=%0.4f\n",I,I-K,J,lastI,lastI-lastK,lastJ,Yc(Y,I,I-K),Yc(Y,lastI,lastI-lastK));
	    printf("\t centerConf= %0.4f, logPV= %0.4f, OutlierConf= %0.4f, min_conf= %0.4e: x=%0.4f,y=%0.4f,m=%d,n=%d:Bias=%0.6e,Pen=%0.6e,Gauss=%0.6e\n",
		   centerConfidence, contig->logPV[m], outlierConfidence, smap_min_conf, fabs(X[qryendidx] - X[qrystartidx]), Yc(Y,I,I-K) - Yc(Y, lastI, lastI-lastK), J - lastJ, K-lastI, Bias,Pen,Gauss);
	    fflush(stdout);
	  }

	  if(outlierConfidence < smap_min_conf)
	    continue;

	  // HERE HERE : include endoutlier confidence correction from output_xmap.cpp

	  SmapCnt++;
	  if(fpSV)
	    fprintf(fpSV,"%d\t%lld\t%lld\t%lld\t%10.1f\t%10.1f\t%10.1f\t%10.1f\t%c\t%10.2f\t%s%dI%dD%d\t%d\t%d\t%0.6f\n",
		  SmapCnt,    // SmapEntryID
		  Xmap->id,   // QryContigID
		  CMapID /* -contigid */,   // RefContigID1
		  CMapID /* -contigid */,   // RefContigID2
		  qrystartpos,
		  qryendpos,
		  refstartpos,
		  refendpos,
		  flip ? '-' : '+', // Orientation
		  outlierConfidence /* WAS contig->logPV[m] */, // Confidence
		  (refendpos-refstartpos > fabs(qryendpos-qrystartpos)) ? "deletion_" : "insertion_", // Type
		  (int)fabs(fabs(qryendpos-qrystartpos) - (refendpos-refstartpos)),// indel size in bp
		  J - lastJ - 1, // inserted sites
		  K - lastI - 1, // deleted sites
		  aligncnt, aligncnt, // XmapID1, XmapID2
		  mapwt);
	}/* lastJ, J loop */
	if(InDelEnds){/* check for end outliers to add to .indel file */
	  if(VERB>=2 && Xmap->id == 1574383LL){
	    printf("Xmap->id=%lld:M=%d,outlier[0]=%d,outlier[M]=%d\n",Xmap->id,M,outlier[0],outlier[M]);
	    fflush(stdout);
	  }
	  if(outlier[0] & 2){/* left endoutlier */
	    SmapCnt++;
	    if(fpSV)
	      fprintf(fpSV,"%d\t%lld\t%lld\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%c\t%0.2f\t%s\t%d\t%d\t%0.6e\n",
		    SmapCnt,    // SmapEntryID
		    Xmap->id,   // QryContigID
		    CMapID /* -contigid */,   // RefContigID1
		    CMapID /* -contigid */,   // RefContigID2
		    flip ? qryLen : 0.0, // Left End of Query region 
		    qrystartpos, // Right End of Query region
		    refstartpos, // Right End of Reference region
		    -1.0,         // Not Defined
		    flip ? '-' : '+', // Orientation
		    -1.0, // Not computed
		    "end", // Type
		    aligncnt, aligncnt, // XmapID1, XmapID2
		    mapwt); 
	  }
	  if(outlier[M] & 2){/* right end outlier */
	    if(fpSV)
	      fprintf(fpSV,"%d\t%lld\t%lld\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%c\t%0.2f\t%s\t%d\t%d\t%0.6e\n",
		    SmapCnt,    // SmapEntryID
		    Xmap->id,   // QryContigID
		    CMapID /* -contigid */,   // RefContigID1
		    CMapID /* -contigid */,   // RefContigID2
		    qryendpos,   // Left End of Query region 
		    flip ? 0.0 : qryLen, // Right End of Query region
		    refendpos, // Left End of Reference region
		    -1.0,         // Not Defined
		    flip ? '-' : '+', // Orientation
		    -1.0, // Not computed
		    "end", // Type
		    aligncnt, aligncnt, // XmapID1, XmapID2
		    mapwt);
	  }
	}
      } /* if(INDEL) */

      /* Output Alignment Cigar */
      fprintf(fp,"%8s"," ");
      int lastI = LI;
      int lastJ = LJ;
      int I= -1,J;
      int matchcnt = 1;
      for(J = LJ+1; J <= RJ; lastI = I, lastJ = J, J++){
	while(J <= RJ && (I = map[J]) < 0)
	  J++;
	if(DEBUG) assert(J <= RJ);
	int K = I - mapL[J];

	/* interval lies between sites lastI .. I-K on reference and lastJ .. J on query */
	if(DEBUG && !(refsite[I] > 0 && refsite[I-K] > 0 && refsite[lastI] > 0)){
	  printf("m=%d(id=%lld),Between J=%d..%d and I=%d..%d,K=%d:refsite[%d]=%d,refsite[%d]=%d,refsite[%d]=%d,outlier=%d\n",
		 m,Xmap->id,lastJ,J,lastI,I-K,K,lastI,refsite[lastI],I-K,refsite[I-K],I,refsite[I],outlier[lastJ]);
	  printf("Hdel[%d]=%d,Hdel[%d]=%d,Hdel[%d]=%d\n",
		 lastI,contig->Hdel[0][lastI],I-K,contig->Hdel[0][I-K],I,contig->Hdel[0][I]);
	  printf("map[%d]=%d(mapL=%d),map[%d]=%d(mapL=%d)\n",
		 lastJ,map[lastJ],mapL[lastJ],J,map[J],mapL[J]);
	  fflush(stdout);
	  assert(refsite[I] > 0 && refsite[lastI] > 0);
	}
	if(VERB>=2){
	  printf("m=%d(id=%lld),Between J=%d..%d and I=%d..%d,K=%d:refsite[%d]=%d,refsite[%d]=%d,refsite[%d]=%d,outlier=%d\n",
		 m,Xmap->id,lastJ,J,lastI,I-K,K,lastI,refsite[lastI],I-K,refsite[I-K],I,refsite[I],outlier[lastJ]);
	  fflush(stdout);
	}

	if(K==0 && refsite[I] == refsite[lastI] + 1 && J == lastJ + 1){
	  matchcnt++;
	} else {
	  /* output match so far */
	  if(matchcnt > 0)
	    fprintf(fp,"%dM",matchcnt);
	  if(refsite[I] == refsite[lastI] + 1)
	    fprintf(fp,"%dI",J - lastJ -1);
	  else if(J == lastJ + 1)
	    fprintf(fp,"%dD", refsite[I] - refsite[lastI] -1);
	  else /* both an insertion and a deletion */
	    fprintf(fp,"%dI%dD", J - lastJ - 1, refsite[I] - refsite[lastI] - 1);
	  matchcnt = 1;
	}
      }/* lastJ, J loop */
      if(matchcnt > 0)
	fprintf(fp,"%dM", matchcnt);
      /* done with cigar */

      if(contig->HapSite[0]){
	double pA = contig->MapPhase[m];
	fprintf(fp,"\t%d", (pA >= HapAlleleProb) ? 1 : (pA <= 1.0-HapAlleleProb) ? 2 : 0);//  Xmap->RunIndex + 1, Xmap->ScanNumber);
      }
      fprintf(fp,"\t%0.1f\t%0.1f\t%d\t",qryLen,refLen, (usecolor==2) ? 2 : 1);

      /* Output Alignment doublet string */
      lastI = LI;
      lastJ = LJ;
      I = -1;
      for(J = LJ; J <= RJ; J++){
	while(J <= RJ && (I = map[J]) < 0)
	  J++;
	if(I < left)
	  continue;
	if(DEBUG) assert(J <= RJ);
	int K = I - max(leftpeak,mapL[J]);
	if(I-K > rightpeak)
	  break;
	if(I > rightpeak){
	  I = min(rightpeak,I);
	  K = I - max(leftpeak,mapL[J]);
	  if(DEBUG) assert(I <= rightpeak && K >= 0);
	}

	if(K > 0){
	  if(DEBUG && !(0 < I-K && I-K <= n && 0 < refsite[I-K] && refsite[I-K] <= refcnt)){
	    if(contig->HapSite[0])
	      printf("J=%d,I=map[J]=%d,mapL[J]=%d,K=%d,n=%d:refsite[I-K]=%d,peak[I-K]=%d,HapSite[0][I-K]=%d,refcnt=%d(left=%d,right=%d,leftpeak=%d,rightpeak=%d)\n",
		     J,I,mapL[J],K,n,refsite[I-K],peak[I-K],contig->HapSite[0][I-K],refcnt,left,right,leftpeak,rightpeak);
	    else
	      printf("J=%d,I=map[J]=%d,mapL[J]=%d,K=%d,n=%d:refsite[I-K]=%d,peak[I-K]=%d,Hdel[0][I-K]=%d,refcnt=%d(left=%d,right=%d,leftpeak=%d,rightpeak=%d)\n",
		     J,I,mapL[J],K,n,refsite[I-K],peak[I-K],contig->Hdel[0][I-K],refcnt,left,right,leftpeak,rightpeak);
	    fflush(stdout);
	    assert(0 < I-K && I-K <= n && 0 < refsite[I-K] && refsite[I-K] <= refcnt);
	  }

	  fprintf(fp,"(%d,%d)",refsite[I-K], flip ? (M+1-J) : J);
	}
        if(DEBUG) assert(0 < I && I <= n && 0 < refsite[I] && refsite[I] <= refcnt);
	fprintf(fp,"(%d,%d)",refsite[I], flip ? (M+1-J) : J);
      }

      //      fprintf(fp,"\t%0.6f\t%0.3f\n", mapwt, (contig->align && contig->align[m]) ? contig->align[m]->maxoutlier : -1.0);
      fprintf(fp,"\t%0.6e\t%0.3f\n", mapwt, (contig->maxoutlier) ? contig->maxoutlier[m] : -1.0);
    } /* loop over maps */

    FILEclose(fp);

    if(INDEL && fpSV){
      FILEclose(fpSV);
      free(write_bufferSV);
    }
    delete [] refsite;
  }

  free(write_buffer);

  delete [] maporder;
  delete [] goodmaps;
  delete [] goodflip;
  CMapID = origCMapID;

  if(1 /* NEW contig->HapSite[0] */){
    errno = 0;
    if(rename(tmpfilename,origfilename)){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to rename %s to %s:errno=%d:%s\n",tmpfilename,origfilename,eno,err);
      fflush(stdout);  exit(1);
    }
  }

  if(VERB){
    printf("Completed output of %s: cum wall time= %0.6f\n",origfilename,wtime());
    fflush(stdout);
  }

  free(origfilename);
  free(tmpfilename);

  delete [] Mask;

  return;
}


#include <sys/types.h>
#include <sys/stat.h>
#ifdef WIN32
#include <io.h>
//#include <stdio.h>
#include <stddef.h>
//#include <iostream>
#define flock(fd, mode) _locking(fd, mode, 1000)
#define LOCK_EX _LK_LOCK
#define LOCK_UN _LK_UNLCK
#else
#include <sys/file.h>
//#include <unistd.h>
#endif
#include <fcntl.h>
//#include <errno.h>

static char buf[LINESIZ];

/* If SplitCnt == 0 : lock the file, read the single integer in the file, increment it and unlock the file

   If SplitCnt != 0 : Just read the value in the file and return it (no file locking is used)
                      If file does not exist try to create it with the largest contig ID in the group_manifest file
 */
long long nextContigId(char *filename)
{
  char absfilename[BUFSIZ], *r;
  if(filename[0] == '/')
    strcpy(absfilename,filename);
  else if(draft_prefix[0] == '/'){
    strcpy(absfilename,draft_prefix);
    r = strrchr(absfilename,'/');
    if(DEBUG) assert(r != NULL);
    strcpy(&r[1], filename);
  } else 
    strcpy(absfilename,filename);

  int repeatcnt = 0;
 Lrepeat:

  errno = 0;
  int fd = open(absfilename, SplitCnt ? O_RDONLY : O_RDWR);
  int eno = errno;
  char *err = strerror(errno);
  if(VERB>=1+RELEASE){
    printf("nextContigID:filename=%s,absfilename=%s,SplitCnt=%lld,fd=%d\n",filename,absfilename,SplitCnt,fd);
    fflush(stdout);
  }
  if(fd < 0 && absfilename[0] == '/' && SplitCnt){  /* locate largest contig id in /<B>/<A><stage>0<N>/<A><stage>0<N>_group_manifest and write the value to absfilename */
    const char *stage;
    if(extendSplit)
      stage = "extension";
    else if(strstr(absfilename,"refineFinal"))
      stage = "refineFinal";
    else if(strstr(absfilename,"refineB"))
      stage = "refineB";
    else {
      printf("Unable to identify refinement stage (refineB, extension or refineFinal) from %s\n",absfilename);
      printf("open(\"%s\") for read/write failed: errno=%d(%s)\n", absfilename, eno, err);
      fflush(stdout); exit(1);
    }

    /* Assume absfilename matches /<B>/<A><stage>1<N>/ID */
   
    int len = strlen(stage);
    char manifest[BUFSIZ];
    strcpy(manifest,absfilename);
    char *E = strstr(manifest,stage);
    while(!E){/* make sure we have the last (rightmost) match in absfilename */
      char *F = strstr(&E[len],stage);
      if(!F)
	break;
      E = F;
    }
    if(!E || (E[len] != '1' && E[len] != '0')){
      if(E)
	printf("Cannot find %s1 or %s0 embedded in %s (E=%s): needed to generate missing file %s\n",stage, stage,absfilename,E, filename);
      else
	printf("Cannot find %s1 or %s0 embedded in %s: needed to generate missing file %s\n",stage, stage, absfilename, filename);
      fflush(stdout);exit(1);
    }
    
    E[len] = '0';// replace "<stage>1" by "<stage>0"

    char *N = strchr(E,'/');
    if(!N){
      printf("Cannot find folder of %s: needed to generate missing file %s\n", absfilename, filename);
      fflush(stdout); exit(1);
    }
    if(strcmp(N,"/ID")){
      printf("absfilename= %s: must end in '/ID', instead found %s\n",absfilename, N);
      fflush(stdout); exit(1);
    }
    *N = '\0';/* manifest should now look like /<B>/<A><stage>0<N> */

    char *A = strrchr(manifest,'/');
    if(!A){
      printf("Cannot find parent folder of %s: needed to generate missing file %s (manifest=%s)\n", absfilename, filename, manifest);
      fflush(stdout); exit(1);
    }
    A++;

    char *Astage0_N = strdup(A);/* should match <A><stage>0<N> */
    if(DEBUG) assert(Astage0_N != NULL);
    sprintf(N,"/%s_group_manifest",Astage0_N);/* manifest should now look like /<B>/<A><stage>0<N>/<A><stage>0<N>_group_manifest */
	
    FILE *fp;
    errno = 0;
    if((fp = fopen(manifest,"r"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("Failed to read group manifest file %s:errno=%d:%s\n",manifest,eno,err);
      printf("absfilename=%s,Astage0_N=%s\n",absfilename,Astage0_N);
      fflush(stdout);exit(1);
    }

    free(Astage0_N);

    long long maxID = 1;
    int linecnt = 1, maxline = -1;
    for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++){
      int len = strlen(buf);
      if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
	printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],manifest,linecnt,buf);
	fflush(stdout);exit(1);
      }
      if(buf[0] == '#')/* comment lines are ignored */
	continue;
      char *pt = buf, *qt;
      for(; *pt; pt = qt){
	long long id = strtoll(pt,&qt,10);
	if((*qt && !isspace(*qt)) || id < 0){
	  printf("Invalid id on line %d of %s:\n%s\n",linecnt,manifest,pt);
	  fflush(stdout);exit(1);
	}
	if(pt == qt)
	  break;
	if(id > maxID){
	  maxID = id;
	  maxline = linecnt;
	}
      }
    }
    (void)fclose(fp);

    if(VERB){
      printf("Largest value in group manifest was %lld on line %d of %s\n",maxID, maxline, manifest);
      printf("Writing value %lld to %s\n",maxID,absfilename);
      fflush(stdout);
    }
    sprintf(manifest,"%s.tmp%lld",absfilename,startCMapID);
    errno = 0;
    if((fp = fopen(manifest,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("Failed to write ID file %s:errno=%d:%s\n",manifest,eno,err);
      fflush(stdout);exit(1);
    }
    fprintf(fp,"%lld\n",maxID);
    fflush(fp);
    (void)fclose(fp);
	
    errno = 0;
    if(rename(manifest,absfilename)){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to rename %s to %s:errno=%d:%s\n",manifest,absfilename,eno,err);
      fflush(stdout);  exit(1);
    }

    return maxID;
  }
  if(fd < 0){
    if(absfilename[0] == '/')
      printf("open(\"%s\") for %s failed: errno=%d:%s\n", absfilename, SplitCnt ? "read" : "read/write", eno, err);
    else {
      char cwd[BUFSIZ];
      printf("open(\"%s\") for %s failed (cwd=%s): errno=%d:%s\n", absfilename, SplitCnt ? "read" : "read/write", getcwd(cwd,BUFSIZ),eno,err);
    }
    fflush(stdout); exit(1);
  }

  if(!SplitCnt){ /* lock the file for exclusive read */

#ifdef WIN32 // OLD : only works if all jobs run on same server
    if(flock(fd, LOCK_EX) < 0){
      char *err = strerror(errno);
      printf("flock of %s failed: errno=%d(%s)\n", absfilename, errno, err);
      fflush(stdout);exit(1);
    }
#else // use fcntl() to lock fd, which should work on NFS mounted files shared across multiple servers
      // Won't work if NFS lockd is not running on v3 NFS or any version of NFS before Linux 2.26.12 (such as on Xeon Phi).
    struct flock lock;
    lock.l_type = F_WRLCK;
    lock.l_whence = SEEK_SET;
    lock.l_start = 0;
    lock.l_len = 16;

    errno = 0;
    if(fcntl(fd,F_SETLKW,&lock) < 0){
      int eno = errno;
      char *err = strerror(errno);
      printf("fcntl(SETLKW,F_WRLCK) of %s failed: errno=%d(%s)\n",absfilename,eno,err);
      fflush(stdout);exit(1);
    }
#endif
  }

  /* read/update integer in file */
  FILE *fp = fdopen(fd, SplitCnt ? "r" : "r+");
  errno = 0;
  if(fp==NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("fdopen() for file %s failed:errno=%d:%s\n",absfilename,eno,err);
    fflush(stdout);exit(1);
  }
  char buf[LINESIZ];
  errno = 0;
  if(fgets(buf,LINESIZ,fp) == NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("Failed to read from %s:errno=%d:%s\n",absfilename,eno,err);
    fflush(stdout);
    if(errno == 116 && repeatcnt < 5){/* stale file handle : another job may have overwritten the file after we opened it : retry up to 5 times */
      printf("Assuming other job overwrote file : retrying %d more times\n",5-repeatcnt);
      fflush(stdout);
      repeatcnt++;
      sleep(1);
      goto Lrepeat;
    }
    exit(1);
  }
  ssize_t cnt = strlen(buf);
  if(cnt >= LINESIZ-1){
    printf("Line too long (LINESIZ=%d) in input file %s\n%s\n",LINESIZ,absfilename,buf);
    fflush(stdout);exit(1);
  }

  char *qt;
  long long ID = strtoll(buf,&qt,10);
  if(qt == buf || !(*qt==0 || isspace(*qt))){
    printf("ID value in %s has invalid syntax :cnt=%lld,buf=%s\n", absfilename,(long long)cnt,buf);
    fflush(stdout);exit(1);
  }
  if(ID <= 0){
    printf("ID value in %s of %lld is invalid (must be positive): cnt=%lld,buf=%s\n",absfilename,ID,(long long)cnt,buf);
    fflush(stdout);exit(1);
  }
  if(VERB>=2 && SplitCnt){
    printf("Read SplitCnt=%lld from %s:%s\n",ID,absfilename,buf);
    fflush(stdout);
  }

  if(!SplitCnt){
    ID++;

    off_t offset = lseek(fd,0,SEEK_SET);
    if(offset < 0){
      char *err = strerror(errno);
      printf("lseek of %s to start of file failed: errno=%d(%s)\n", absfilename, errno, err);
      fflush(stdout);exit(1);
    }
    sprintf(buf,"%lld\n",ID);
    cnt = write(fd,buf,strlen(buf));
    if(cnt < 0){
      char *err = strerror(errno);
      printf("write of %lld to %s failed: errno=%d(%s)\n", ID, absfilename, errno, err);
      fflush(stdout);exit(1);
    }
    if((size_t)cnt != (size_t)strlen(buf)){
      printf("write to %s failed (bytes written=%ld):%s\n", absfilename, cnt, buf);
      fflush(stdout);exit(1);
    }

    /* unlock file */
#if WIN32 // OLD method, only works with all jobs on same server
    if(flock(fd,LOCK_UN) < 0){
      char *err = strerror(errno);
      printf("unlock of %s failed: errno=%d(%s)\n", absfilename, errno, err);
      fflush(stdout);exit(1);
    }
#else
    // No need to unlock, since close(fd) releases all fcntl locks
#endif
  }

  if(FILEclose(fp) != 0){
    char *err = strerror(errno);
    printf("%s:FILEclose() failed for file %s\n",err,absfilename);
    fflush(stdout);exit(1);
  }
#if 0
  if(close(fd) < 0){
    char *err = strerror(errno);
    printf("close of file descriptor to %s failed: errno=%d(%s)\n", absfilename, errno, err);
    fflush(stdout);exit(1);
  }
#endif  

  return ID;
}


void output_draft(Ccontig *contig, long long contigid, char *basename, int Lfrozen, int Rfrozen, int Allele)
{
  if(DEBUG) assert(Rfrozen >= Lfrozen);
  if(colors > 1){
    printf("output_draft: not implemented for colors=%d (usecolor=%d)\n",colors,usecolor);
    fflush(stdout);
    if(DEBUG) assert(!usecolor);
    fflush(stdout);exit(1);
  }
  
  if(VERB && Refine){
    printf("output_draft:CMapID=%lld,contigid=%lld,Lfrozen=%d,Rfrozen=%d,Allele=%d,HapSite[0]=%p\n",CMapID,contigid,Lfrozen,Rfrozen,Allele,contig->HapSite[0]);
    int n = contig->numsite[0];
    printf("\t MaskL = 0x%lx, MaskR = 0x%lx, n=%d\n", contig->MaskL, contig->MaskR, n);

    fflush(stdout);
  }

  /* first compute consensus map to output */
  int n = contig->numsite[0];
  double *Y = contig->site[0];
  double *sitecnt = new double[n+2];
  double *smoothcnt = new double[n+2];
  int *peak = new int[n+2];

  /* copy raw sitecnt */
  sitecnt[0] = sitecnt[n+1] = 0.0;
  for(register int i = 1; i <= n; i++)
    sitecnt[i] = contig->sitecnt[0][i];

  register double reskb = res[0] * PixelLen;
  register double sd = draftSD * reskb;

  if(!Refine){/* must be called from Assembler for draft contig */
    if(COVERAGE_FIX)
      for(int i = 0; i <= n+1; i++)
	contig->fragcntT[0][i] = contig->fragcnt[0][i];// Assembler does not compute fragcntT

    if(VERB>=2){
      for(register int m = 0; m < contig->nummaps; m++){
	printf("m=%d:mapid=%d,id=%lld:M=%d,N=%d,flip=%d:\n",m, contig->contig[m].mapid, Gmap[contig->contig[m].mapid]->id, contig->contig[m].numsite[0],contig->numsite[0], contig->flip[m]);
	int trimL = contig->contig[m].trimL[0];
	int trimR = contig->contig[m].trimR[0];
	for(register int j = 1; j <= contig->contig[m].numsite[0] + 1; j++)
	  if(contig->flip[m])
	    printf("  X'[%d]=%0.3f, sitemap[m][%d] = %d\n",j, Gmap[contig->contig[m].mapid]->site[0][trimR] - Gmap[contig->contig[m].mapid]->site[0][trimR - j], j, contig->sitemap[0][m][j]);
	  else
	    printf("  X[%d]=%0.3f, sitemap[m][%d] = %d\n",j, Gmap[contig->contig[m].mapid]->site[0][j+trimL], j, contig->sitemap[0][m][j]);
      }
      fflush(stdout);
    }

    /* first compute a smoothed the sitecnt value and determine its significant peaks */

    double W = 6.0*sd;/* window half width */
    double Ivar = 1.0/(2.0*sd*sd);

    smoothcnt[0] = smoothcnt[n+1] = 0.0;
    int left = 1, right = 1;/* left and right end of window */
    for(int i = 1; i <= n; i++){
      /* move ahead left end of window until it is W to the left of Y[i] */
      while(Y[left] < Y[i] - W)
	left++;

      /* move ahead right end of window until it is W to the right of Y[i] */
      while(right <= n && Y[right] < Y[i] + W)
	right++;
    
      /* now smooth over window left ... right-1 */
      register double sum = 0.0;
      for(register int k = left; k < right; k++){
	register double err = Y[i] - Y[k];
	sum += sitecnt[k] * exp(-err*err*Ivar);
      }
      smoothcnt[i] = sum;
    }
  
    if(VERB>=2)
      printf("sd=%0.4f,reskb=%0.3f\n",sd,reskb);

    /* locate local peaks in smoothcnt[1..n] */
    peak[0] = peak[n+1] = 0;
    for(register int i = 1; i <= n; i++){
      /* scan forward to first site that differs */
      for(right = i+1; right <= n+1;right++){
	if(Y[right] - Y[i] > sd)
	  break;
	if(fabs(smoothcnt[right] - smoothcnt[i]) >= 0.0001)
	  break;
      }
      /*scan backward to first site that differs */
      for(left = i; --left >= 0; ){
	if(Y[i] - Y[left] > sd)
	  break;
	if(fabs(smoothcnt[left] - smoothcnt[i]) >= 0.0001)
	  break;
      }


      if(left < 0 || right > n+1 || smoothcnt[i] <= 0.0)
	peak[i] = 0;
      else
	peak[i] = ((Y[right] - Y[i] > sd || smoothcnt[i] > smoothcnt[right]) && 
		   (Y[i] - Y[left] > sd || smoothcnt[i] > smoothcnt[left])) ? 1 : 0;

      if(VERB>=2)
	printf("i=%d:site[i]=%0.3f,smoothcnt[i]=%0.4f,smoothcnt[left=%d]=%0.4f,smoothcnt[right=%d]=%0.4f,peak=%d\n",
	       i,Y[i],smoothcnt[i],left,smoothcnt[left],right,smoothcnt[right],peak[i]);
    }


    /* remove peaks that are within res of a larger peak */

    for(register int i = 1; i <= n; i++){
      if(peak[i]){
	/* scan left */
	for(register int k = i; --k >= 1;){
	  if(Y[i] - Y[k] > reskb)
	    break;
	  if(peak[k] && smoothcnt[k] > smoothcnt[i]){
	    peak[i] = -1;
	    break;
	  }
	}
	if(peak[i] < 0)
	  continue;
	/* scan right */
	for(register int k = i; ++k <= n; ){
	  if(Y[k] - Y[i] > reskb)
	    break;
	  if(peak[k] && smoothcnt[k] >= smoothcnt[i]){/* priority to the rightmost peak with the same value */
	    peak[i] = -1;
	    break;
	  }
	}
      }
    }
  } else {/* contig has been refined */
    if(contig->Hdel[0]){/* new method */
      if(contig->HapSite[0]){/* Haplotype map : consider any site with HapSite[0][i] != 0 as a confirmed site */
	for(int k=1;k <= n; k++){
	  smoothcnt[k] = contig->HapSite[0][k] ? contig->fragcnt[0][k-1] + 1.0 : 0.0;
	  peak[k] = contig->HapSite[0][k] ? 1 : 0;
	}
      } else {/* consider any site with Hdel[0][i] == 0 as a confirmed site */
	for(register int k=1;k <= n; k++){
	  smoothcnt[k] = (contig->Hdel[0][k]==0) ? contig->fragcnt[0][k-1] + 1.0 : 0.0;
	  peak[k] = (contig->Hdel[0][k]==0) ? 1 : 0;
	}
      }
    } else { /* old method just consider any site with sitecnt[i] > 0 as a confirmed site */
      for(register int k=1;k <= n; k++){
	smoothcnt[k] = sitecnt[k] ? contig->fragcnt[0][k-1] + 1.0 : 0.0;
	peak[k] = (sitecnt[k] ? 1 : 0);
      }
    }
  }

  /* estimate current location of contig->left and contig->right based on distance from original Y[Lfrozen,Rfrozen] stored in contig->Lfrozen, contig->Rfrozen */
  /*  NOTE : This should match contig->left .. contig->right, but using Y[Lfrozen .. Rfrozen] is more rebust since refinement tries hard to keep Y[Lfrozen,Rfrozen] valid after stretching of consensus map due to refinement */
  //  if(DEBUG && SegDupMask && contig->SegDupCnt > 0) assert(extend <= 1);// since only refineFinal should propagate SegDupMask
  if(DEBUG && SegDupMask && contig->SegDupCnt && extend >= 2) assert(Rfrozen == 0);


  double newleft = (Rfrozen == 0)/*NEW139*/ ? contig->left : Y[Lfrozen] + (contig->left - contig->Lfrozen);
  double newright = (Rfrozen == 0)/*NEW139*/ ? contig->right : Y[Rfrozen] - (contig->Rfrozen - contig->right);
  if(VERB/* HERE >=2 */){
    printf("contig->left= %0.4f, contig->right= %0.4f remapped to newleft= %0.4f, newright= %0.4f, Lfrozen=%d, Rfrozen=%d,n=%d,extendonly=%d,MaskL=%lu,MaskR=%lu(contig->Lfrozen=%0.4f,contig->Rfrozen=%0.4f,Y[n+1]=%0.4f)\n", 
	   contig->left,contig->right,newleft,newright,Lfrozen,Rfrozen,n,extendonly,contig->MaskL,contig->MaskR,contig->Lfrozen,contig->Rfrozen,Y[n+1]);
    fflush(stdout);
  }
  
  /* determined how much of Consensus ends should be trimmed */
  int left = 1, right = n;
  if(/* WAS 0 && */Refine <= 1 && EndTrimCov <= 0.0){/* avoid any trimming */
    left = 1;
    right = n;
  } else {
    if(EndTrimCov > 0.0){
      int origleft = left, origright = right;
      if(extendonly && Rfrozen > Lfrozen){ /* limit trimming to exclude sites Y[Lfrozen .. Rfrozen] UNLESS MaskL (or MaskR) AND CloneTrim is true */
	while(left < right && (((Y[left] < Y[Lfrozen] || (contig->MaskL && CloneTrim)) && contig->fragcnt[0][left] < EndTrimCov) || !peak[left])){
	  if(VERB>=2){
	    printf("left=%d,right=%d: Y[left,Lfrozen]= %0.4f,%0.4f : peak[left]=%d, contig->fragcnt[0][left]= %0.1f, EndTrimCov= %0.1f: trimming left end\n",
		   left,right,Y[left], Y[Lfrozen], peak[left],contig->fragcnt[0][left],EndTrimCov);
	    fflush(stdout);
	  }
	  left++;
	}

	while(right > left && (((Y[right] > Y[Rfrozen] || (contig->MaskR && CloneTrim)) && contig->fragcnt[0][right - ((FRAGCOV_FIX&&Refine)?1:0)] < EndTrimCov) || !peak[right])){
	  if(VERB>=2){
	    printf("left=%d,right=%d: Y[right,Rfrozen]= %0.4f,%0.4f : peak[right]=%d, contig->fragcnt[0][right]= %0.1f, contig->fragcnt[0][right-1]= %0.1f, EndTrimCov= %0.1f: trimming right end\n",
		   left,right,Y[right], Y[Rfrozen],peak[right],contig->fragcnt[0][right],contig->fragcnt[0][right-1],EndTrimCov);
	    fflush(stdout);
	  }

	  right--;
	}
	
	if(EndTrimCov2 > 0.0){/* further trim back to any internal region with cov < EndTrimCov2, if beyond original ends (or within EndTrimLen2 of original ends) */
	  //	  int origleft2 = left, origright2 = right;
	  double leftlim = newleft + EndTrimLen2;
	  double rightlim = newright - EndTrimLen2;
	  for(int t = left+1; t < right; t++){
	    if(Y[t] >= rightlim && peak[t] && contig->fragcnt[0][t] < EndTrimCov2){
	      right = t;
	      break;
	    }
	  }

	  for(int t = right-1; t > left; t--){
	    if(Y[t-1] <= leftlim && peak[t-1] && contig->fragcnt[0][t-1] < EndTrimCov2){
	      left = t;
	      break;
	    }
	  }
	}

	if(contig->MaskL && CloneTrim && left > 1 && Y[left] > Y[Lfrozen])
	  contig->MaskL = 0;

	if(contig->MaskR && CloneTrim && right < n && Y[right] < Y[Rfrozen])
	  contig->MaskR = 0;

      } else {
	while(left < right && (contig->fragcnt[0][left] < EndTrimCov || !peak[left])){
	  if(VERB>=2){
	    printf("left=%d,right=%d: Y[left,Lfrozen]= %0.4f,%0.4f : peak[left]=%d, contig->fragcnt[0][left]= %0.1f, EndTrimCov= %0.1f: trimming left end\n",
		   left,right,Y[left], Y[Lfrozen], peak[left],contig->fragcnt[0][left],EndTrimCov);
	    fflush(stdout);
	  }
	  left++;
	}
	while(right > left && (contig->fragcnt[0][right - ((FRAGCOV_FIX&&Refine)?1:0)] < EndTrimCov || !peak[right])){
	  if(VERB>=2){
	    printf("left=%d,right=%d: Y[right]= %0.4f : peak[right]=%d, contig->fragcnt[0][right]= %0.1f, contig->fragcnt[0][right-1]= %0.1f, EndTrimCov= %0.1f: trimming right end\n",
		   left,right,Y[right], peak[right],contig->fragcnt[0][right],contig->fragcnt[0][right-1],EndTrimCov);
	    fflush(stdout);
	  }
	  right--;
	}

	if(EndTrimCov2 > 0.0){/* further trim back to any internal region with cov < EndTrimCov2, if beyond original ends (or within EndTrimLen2 of original ends) */
	  //	  int origleft2 = left, origright2 = right;
	  double leftlim = newleft + EndTrimLen2;
	  double rightlim = newright - EndTrimLen2;
	  for(int t = left+1; t < right; t++){
	    if(Y[t] >= rightlim && peak[t] && contig->fragcnt[0][t] < EndTrimCov2){
	      right = t;
	      break;
	    }
	  }

	  for(int t = right-1; t > left; t--){
	    if(Y[t-1] <= leftlim && peak[t-1] && contig->fragcnt[0][t-1] < EndTrimCov2){
	      left = t;
	      break;
	    }
	  }
	}
      }
      if(VERB/* HERE >=2 */ && (left != origleft || right != origright)){
	printf("Trimmed ends due to EndTrimCov= %0.2f, CloneTrim=%d, EndTrimInternal=%0.2f,%0.2f: left= %d -> %d (%0.3f -> %0.3f), right= %d -> %d (%0.3f -> %0.3f), MaskL= %lx, MaskR= %lx\n",
	       EndTrimCov,CloneTrim, EndTrimCov2,EndTrimLen2,origleft,left,Y[origleft],Y[left],origright,right,Y[origright],Y[right],contig->MaskL,contig->MaskR);
	fflush(stdout);
      }
    }
    if(!extendonly && extend<=1 && Rfrozen > Lfrozen){/* trim off extension region outside of Y[Lfrozen .. Rfrozen] */
      while(left < right && (left < Lfrozen || !peak[left]))
	left++;
      while(right > left && (right > Rfrozen || !peak[right]))
	right--;
    }
  }
  //  double leftend = (left > 0 ? MININTERVAL : 0.0);
  //  double rightend = (right <= n ? MININTERVAL : 0.0);
  double leftend = (left > 0) ? (Refine == 1 /* WAS && extend <= 1 */ && EndTrimCov <= 0.0) ? max(MININTERVAL, Y[left] - contig->left) : MININTERVAL : 0.0;/* extension on left end beyond Y[left] */
  double rightend = (right <= n) ? (Refine == 1 /* WAS && extend <= 1 */ && EndTrimCov <= 0.0) ? max(MININTERVAL, contig->right - Y[right]) : MININTERVAL : 0.0;/* extension on right end beyond Y[right] */
  
  if(VERB>=2 && Refine==1){
    printf("Y[left]=%0.6f,contig->left=%0.6f,leftend=%0.6f, Y[right]=%0.6f,contig->right=%0.6f,rightend=%0.6f\n",
	   Y[left],contig->left,leftend,Y[right],contig->right,rightend);
    fflush(stdout);
  }

  if(VERB/* >=2 */){
     if(Refine){
       if(Rfrozen > Lfrozen){
	 if(extend>=2)
	   printf("    numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f,cov=%0.1f),right=%d(%0.3f,cov=%0.1f),Raw Length=%0.3f,Trimmed Length=%0.3f):original ends: left=%0.3f,right=%0.3f,Lfrozen=%d(%0.3f),Rfrozen=%d(%0.3f)\n",
		  n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],contig->fragcnt[0][left],right,Y[right], contig->fragcnt[0][right - ((FRAGCOV_FIX&&Refine)?1:0)],
		  Y[n+1], Y[right] - Y[left]+leftend+rightend, contig->left,contig->right,Lfrozen,Y[Lfrozen],Rfrozen,Y[Rfrozen]);
	 else
	   printf("    numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f,cov=%0.1f),right=%d(%0.3f,cov=%0.1f),Raw Length=%0.3f,Trimmed Length=%0.3f),Lfrozen=%d(%0.3f),Rfrozen=%d(%0.3f)\n",
		  n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],contig->fragcnt[0][left],right,Y[right], contig->fragcnt[0][right- ((FRAGCOV_FIX&&Refine)?1:0)],
		  Y[n+1], Y[right] - Y[left] + leftend+rightend, Lfrozen,Y[Lfrozen],Rfrozen,Y[Rfrozen]);
       } else {
	 if(extend>=2)
	   printf("    numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f,cov=%0.1f),right=%d(%0.3f,cov=%0.1f),Raw Length=%0.3f,Trimmed Length=%0.3f):original ends: left=%0.3f,right=%0.3f\n",
		  n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],contig->fragcnt[0][left],right,Y[right], contig->fragcnt[0][right- ((FRAGCOV_FIX&&Refine)?1:0)],
		  Y[n+1], Y[right] - Y[left] +leftend+rightend, contig->left,contig->right);
	 else
	   printf("    numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f,cov=%0.1f),right=%d(%0.3f,cov=%0.1f),Raw Length=%0.3f,Trimmed Length=%0.3f)\n",
		  n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],contig->fragcnt[0][left],right,Y[right], contig->fragcnt[0][right- ((FRAGCOV_FIX&&Refine)?1:0)],
		  Y[n+1], Y[right]-Y[left]+leftend+rightend);
       }
     } else
       printf("contig%lld:numsite=%d(tot=%d),maps=%d,sd=%0.3f(EndTrimCov=%0.2f,left=%d(%0.3f,cov=%0.1f),right=%d(%0.3f,cov=%0.1f), Raw Length=%0.3f,Trimmed Length=%0.3f)\n",
	      contigid,n,contig->totsites(),contig->nummaps,sd,EndTrimCov,left, Y[left], contig->fragcnt[0][left], right, Y[right], contig->fragcnt[0][right- ((FRAGCOV_FIX&&Refine)?1:0)],
	      Y[n+1], Y[right]-Y[left]+leftend+rightend);
  }
  if(BVERB && Refine){
    register int k;
    printf("n=%d,left=%d,right=%d,site[left]=%0.3f,site[right]=%0.3f\n",n,left,right,Y[left],Y[right]);
    for(k = 1; k <= n; k++){
      if(SplitSite && !peak[k] /* contig->sitecnt[0][k] <= 0.0 */)
	continue;
      if(contig->fragcntT[0]){
	if(contig->sitecntFNnorm[0]){
	  printf("k=%d:site[k]=%0.3f,sitecnt[k]=%0.1f,sitecntFN[k]=%0.1f(norm=%0.1f),fragcnt[k]=%0.1f,fragcntT[k]=%0.1f,smoothcnt[k]=%0.4f,peak[k]=%d%c\n",
		 k,Y[k],contig->sitecnt[0][k],contig->sitecntFN[0][k],contig->sitecntFNnorm[0][k],contig->fragcnt[0][k],contig->fragcntT[0][k],smoothcnt[k],peak[k],(peak[k]>0 && smoothcnt[k] > contig->fragcnt[0][k-1]*draftTP) ? '!':' ');
	  if(ChimQuality >= 2){
	    printf("\t sitecntN2[k]=%0.1f,sitecntN3[k]=%0.1f,sitecntN4[k]=%0.1f,sitecntN5[k]=%0.1f,sitecntN6[k]=%0.1f",
		   contig->sitecntN2[0][k], contig->sitecntN3[0][k], contig->sitecntN4[0][k], contig->sitecntN5[0][k], contig->sitecntN6[0][k]);
	    if(OutlierQuality >= 2)
	      printf(",FragSd=%0.1f,ExpSd=%0.1f,Cov=%0.2f,ChiSq=%0.2e",contig->fragSd[0][k]*1000.0,contig->expSd[0][k]*1000.0,contig->fragCov[0][k],contig->fragChiSq[0][k]);
	    printf("\n");
	  }
	} else
	  printf("k=%d:site[k]=%0.3f,sitecnt[k]=%0.1f,sitecntFN[k]=%0.1f,fragcnt[k]=%0.1f,fragcntT[k]=%0.1f,smoothcnt[k]=%0.4f,peak[k]=%d%c\n",
		 k,Y[k],contig->sitecnt[0][k],contig->sitecntFN[0][k],contig->fragcnt[0][k],contig->fragcntT[0][k],smoothcnt[k],peak[k],(peak[k]>0 && smoothcnt[k] > contig->fragcnt[0][k-1]*draftTP) ? '!':' ');
      } else
	printf("k=%d:site[k]=%0.3f,sitecnt[k]=%0.1f,sitecntFN[k]=%0.1f,fragcnt[k]=%0.1f,smoothcnt[k]=%0.4f,peak[k]=%d%c\n",
	       k,Y[k],contig->sitecnt[0][k],contig->sitecntFN[0][k],contig->fragcnt[0][k],smoothcnt[k],peak[k],(peak[k] > 0 && smoothcnt[k] > contig->fragcnt[0][k-1]*draftTP) ? '!':' ');
    }
    printf("k=%d:site[k]=%0.3f,fragcnt[k-1]=%0.1f\n", k,Y[k],contig->fragcnt[0][k-1]);
  }

  if(!Refine){ /* reset peak[k] = 1 to 0 unless smoothcnt[k] > contig->fragcnt[k-1]*draftTP */
    for(register int k = 1; k <= n; k++)
      if(peak[k] > 0 && !(smoothcnt[k] > contig->fragcnt[0][k-1]*draftTP))
	peak[k] = 0;
  }

  if(VERB>=2){
    printf("Refine=%d,ContigSplitRatio= %0.2f, ContigSplitRatio2= %0.2f,extendonly=%d,Rfrozen=%d,CMapID=%lld,contig->HapSite[0]=%p, EndTrimCQ=%d\n",
	   Refine,ContigSplitRatio,ContigSplitRatio2,extendonly,Rfrozen,CMapID,contig->HapSite[0],EndTrimCQ);
    fflush(stdout);
  }

  if(Refine && ContigSplitRatio > 0.0 && ContigSplitRatio2 > 0.0 && !(extendonly && Rfrozen > 0) && /* NEW */ !(CMapID > 0 && contig->HapSite[0])){
    if(EndTrimCQ){ /* NEW89 : trim ends (Y[left..right]) up to COVERAGE_TRIM_LEN from start of contig->sitecntFN[0][k] != END_CHIMQUAL */
      int left1 = left, right1 = right;
      int left2 = left, right2 = right;
      while(left2 < right2 && contig->sitecntFN[0][left2] == END_CHIMQUAL)
	left2++;
      while(left2 < right2 && contig->sitecntFN[0][right2] == END_CHIMQUAL)
	right2--;
      

      int origleft2 = left2, origright2 = right2;

      /* NEW99 : back off by 1 real label */
      for(left2--; left2 > left1; left2--)
	if(peak[left2])
	  break;

      for(right2++; right2 < right1; right2++)
	if(peak[right2])
	  break;

      while(left < left2 && Y[left2] - Y[left+1] >= COVERAGE_TRIM_LEN)
	left++;
      while(right2 < right && Y[right-1] - Y[right2] >= COVERAGE_TRIM_LEN)
	right--;

      /* Backoff trimming to nearest real label */
      while(left > left1 && !peak[left])
	left--;
      while(right < right1 && !peak[right])
	right++;
      
      if(VERB){
	printf("Trimmed ends from Y[%d..%d](%0.4f..%0.4f) to Y[%d..%d](%0.4f..%0.4f) due to ChimQual range Y[%d..%d](%0.4f..%0.4f) and range+1 Y[%d..%d](%0.4f..%0.4f)\n",
	       left1,right1,Y[left1],Y[right1],left,right,Y[left],Y[right],origleft2,origright2,Y[origleft2],Y[origright2],left2,right2,Y[left2],Y[right2]);
	fflush(stdout);
      }
    }

    /* check for internal region with local minimum in sitecntFN[] that is below ContigSplitRatio times average value of sitecntFN */
    if(DEBUG) assert(contigid <= 0); /* since ContigSplitRatio can only be set by RefAligner */
    
    if(SplitCnt < 0){
      SplitCnt = nextContigId(ContigCntFile);
      if(DEBUG)assert(SplitCnt > 0);
      if(VERB/* HERE >=2 */){
	printf("nextContigId(%s) returned SplitCnt=%lld\n",ContigCntFile,SplitCnt);
	fflush(stdout);
      }
    }
    if(DEBUG && CMapID > SplitCnt){
      printf("ERROR: CMapID = %lld is larger than value of %lld in ID file %s (Invalid Pipeline restart ?)\n",CMapID, SplitCnt, ContigCntFile);
      fflush(stdout);
      exit(1);
    }

    int NumSplits = 0;

    /* compute median value of sitecntFN and fragcnt (excluding any 0 values) */
    float *sitecntFN = new float[right-left+1];
    float *fragcnt = new float[right-left+1];
    int lenFN = 0, lenFC = 0;
    
    if(SplitSite){
      float PsitecntFN = -1.0, NsitecntFN = -1.0;
      for(int i = left; i <= right; i++){
	if(!peak[i] /* contig->sitecnt[0][i] <= 0.0*/){
	  if(PsitecntFN < 0.0){
	    contig->sitecntFN[0][i] = END_CHIMQUAL;/* left end region */
	    continue;
	  }
	  if(NsitecntFN < 0.0){
	    for(int j = i+1; j <= n; j++)
	      if(contig->sitecnt[0][j] > 0.0){
		NsitecntFN = contig->sitecntFN[0][j];
		break;
	      }
	  }
	  if(NsitecntFN >= 0.0)
	    contig->sitecntFN[0][i] = 0.5*(PsitecntFN + NsitecntFN);
	  else
	    contig->sitecntFN[0][i] = END_CHIMQUAL;/* right end region */
	  continue;
	}
	sitecntFN[lenFN++] = PsitecntFN = contig->sitecntFN[0][i];
	fragcnt[lenFC++] = FRAGCOV_FIX ? max(contig->fragcnt[0][i],contig->fragcnt[0][max(0,i-1)]) : contig->fragcnt[0][i];
	NsitecntFN = -1.0;
      }
    } else {
      for(register int i = left; i <= right; i++){
	sitecntFN[i-left] = contig->sitecntFN[0][i];
	fragcnt[i-left] = FRAGCOV_FIX ? max(contig->fragcnt[0][i],contig->fragcnt[0][max(0,i-1)]) : contig->fragcnt[0][i];
      }
      lenFN = right-left+1;
      lenFC = right-left+1;
    }

    qsort(sitecntFN,lenFN,sizeof(float),(intcmp *)floatdec);
    qsort(fragcnt,lenFC,sizeof(float),(intcmp *)floatdec);
    while(lenFN > 0 && sitecntFN[lenFN-1] <= 0.5)
      lenFN--;
    while(lenFC > 0 && fragcnt[lenFC-1] <= 0.5)
      lenFC--;
    float medianFN = (lenFN%2) ? sitecntFN[(lenFN-1)/2] : (lenFN==0) ? 50.0 : (sitecntFN[(lenFN-1)/2] + sitecntFN[lenFN/2])/2;
    float medianFC = (lenFC%2) ? fragcnt[(lenFC-1)/2] : (lenFC==0) ? 2.0 : (fragcnt[(lenFC-1)/2] + fragcnt[lenFC/2])/2;
    delete [] sitecntFN;
    delete [] fragcnt;

    int L = left, R = right;
    float threshold = SplitRev ? (TrimNorm >= 0 && TrimNormMed ? TrimNormMed : medianFN) * ContigSplitRatio2 : min(medianFC * ContigSplitRatio, medianFN * ContigSplitRatio2);
    if(1 /* HERE BVERB */){
      printf("median:sitecntFN=%0.1f,fragcnt=%0.1f,left=%d,right=%d,lenFN=%d,lenFC=%d,TrimNormMed=%d:threshold=%0.1f\n",medianFN,medianFC,left,right,lenFN,lenFC,TrimNormMed,threshold);
      fflush(stdout);
    }
    if(DEBUG && !isfinite(threshold)){
      printf("n=%d,left=%d,right=%d,site[left]=%0.3f,site[right]=%0.3f\n",n,left,right,Y[left],Y[right]);
      for(int k = 1; k <= n; k++){
	if(contig->fragcntT[0]){
	  if(contig->sitecntFNnorm[0])
	    printf("k=%d:site[k]=%0.3f,sitecnt[k]=%0.1f,sitecntFN[k]=%0.1f(norm=%0.1f),fragcnt[k-1]=%0.1f,fragcntT[k-1]=%0.1f,smoothcnt[k]=%0.4f,peak[k]=%d%c\n",
		   k,Y[k],contig->sitecnt[0][k],contig->sitecntFN[0][k],contig->sitecntFNnorm[0][k],contig->fragcnt[0][k-1],contig->fragcntT[0][k-1],smoothcnt[k],peak[k],(peak[k]>0 && smoothcnt[k] > contig->fragcnt[0][k-1]*draftTP) ? '!':' ');
	  else
	    printf("k=%d:site[k]=%0.3f,sitecnt[k]=%0.1f,sitecntFN[k]=%0.1f,fragcnt[k-1]=%0.1f,fragcntT[k-1]=%0.1f,smoothcnt[k]=%0.4f,peak[k]=%d%c\n",
		   k,Y[k],contig->sitecnt[0][k],contig->sitecntFN[0][k],contig->fragcnt[0][k-1],contig->fragcntT[0][k-1],smoothcnt[k],peak[k],(peak[k]>0 && smoothcnt[k] > contig->fragcnt[0][k-1]*draftTP) ? '!':' ');
	} else
	  printf("k=%d:site[k]=%0.3f,sitecnt[k]=%0.1f,sitecntFN[k]=%0.1f,fragcnt[k-1]=%0.1f,smoothcnt[k]=%0.4f,peak[k]=%d%c\n",
		 k,Y[k],contig->sitecnt[0][k],contig->sitecntFN[0][k],contig->fragcnt[0][k-1],smoothcnt[k],peak[k],(peak[k]>0 && smoothcnt[k] > contig->fragcnt[0][k-1]*draftTP) ? '!':' ');
      }
      printf("median:sitecntFN=%0.1f,fragcnt=%0.1f,left=%d,right=%d,lenFN=%d,lenFC=%d,threshold=%0.1f\n",medianFN,medianFC,left,right,lenFN,lenFC,threshold);
      fflush(stdout);
      assert(isfinite(threshold));
    }
    if(SplitRev <= 1){    /* skip ends, until sitecntFN > threshold */
      while(L < R && contig->sitecntFN[0][L] <= threshold)
	L++;
      while(L < R && contig->sitecntFN[0][R] <= threshold)
	R--;
    }
    if(BVERB>=2 && Refine){
      printf("After SplitSite and SplitRev adjustments: threadshold=%0.1f,L=%d,R=%d:\n",threshold,L,R);
      for(int k = 1; k <= n; k++)
	printf("k=%d:site[k]=%0.3f,sitecnt[k]=%0.1f,sitecntFN[k]=%0.1f,peak[k]=%d\n",
	       k,Y[k],contig->sitecnt[0][k],contig->sitecntFN[0][k],peak[k]);
      fflush(stdout);
    }
    if(L < R-1){
      /* Check for internal sitecntFN minima between L and R where sitecntFN < threshold.
	 If SplitFragileD > 0, also check for one sides Fragile sites (chromosome ends) in the middle of a contig (see -splitFragile) 
       */
      int start, end = L;
      for(register int k = L+1; k < R; k = end+1){
	start = k-1;

	/* locate start and end of next local minima sitecntFN[start .. end] */
	for(end = k; end < R;end++){
	  if(contig->sitecntFN[0][end+1] < contig->sitecntFN[0][end]){/* possible start of local minima region */
	    start = end+1;
	    continue;
	  }
	  if(start >= k && contig->sitecntFN[0][end+1] > contig->sitecntFN[0][end])/* end of local minima region */
	    break;
	}

	if(start < k || end >= R)
	  continue;

	/* NEW74 : locate k s.t. Y[k] is closest to (Y[start] + Y[end])/2.0 */
	double Ymid = (Y[start] + Y[end]) * 0.5, err = Ymid - Y[start];

	//	k = (start+end)/2;/* midpoint of minima region */

	k = start;
	for(int t = start+1; t <= end; t++){
	  if(fabs(Y[t] - Ymid) < err){
	    err = fabs(Y[t] - Ymid);
	    k = t;
	  }
	}
	if(DEBUG && !(start <= k && k <= end)){
	  printf("start=%d,end=%d:Y[start]=%0.3f,Ymid=%0.3f,Y[end]=%0.3f,k=%d,Y[k]=%0.3f,err=%0.3f\n",start,end,Y[start],Ymid,Y[end],k,Y[k],err);
	  fflush(stdout);
	  assert(start <= k && k <= end);
	}


	/* scan beyond local minima to make sure FN increases by at least ContigSplitRatio3 on either side before dropping below this level */
	float minFN = (SplitRev ? contig->sitecntFN[0][k] : max(contig->sitecntFN[0][k],threshold)) * ContigSplitRatio3;
	int kleft = start, kright = end;
	float rightFN = contig->sitecntFN[0][end];
	for(register int t = end+1; t <= R; t++){
	  if(contig->sitecntFN[0][t] > rightFN){
	    rightFN = contig->sitecntFN[0][t];
	    kright = t;
	    if(rightFN > minFN)
	      break;
	  }
	  if(contig->sitecntFN[0][t] < contig->sitecntFN[0][end])
	    break;
	}

	float leftFN = contig->sitecntFN[0][start];
	for(register int t = start-1; t >= L; t--){
	  if(contig->sitecntFN[0][t] > leftFN){
	    leftFN = contig->sitecntFN[0][t];
	    kleft = t;
	    if(leftFN > minFN)
	      break;
	  }
	  if(contig->sitecntFN[0][t] < contig->sitecntFN[0][end])
	    break;
	}

	if(leftFN < minFN)/* invalid minima */
	  continue;

	if(rightFN < minFN)/* invalid minima */
	  continue;

	if(contig->sitecntFN[0][k] >= threshold) /* minima is too high */
	  continue;

	if(BVERB){
	  printf("Found minima at k=%d(range=%d..%d),site[k]=%0.3f:sitecntFN[k]=%0.1f,fragcnt[k-1]=%0.1f,fragcnt[k]=%0.1f,leftFN=%0.1f(at site[%d]=%0.3f),rightFN=%0.1f(at site[%d]=%0.3f),minFN=%0.1f,ContigSplitRatio=%0.4f,%0.4f,%0.4f:medianFN=%0.1f,medianFC=%0.1f,threshold=%0.1f\n",
		 k,start,end,Y[k],contig->sitecntFN[0][k],contig->fragcnt[0][k-1],contig->fragcnt[0][k],leftFN,kleft,Y[kleft],rightFN,kright,Y[kright],minFN,ContigSplitRatio,ContigSplitRatio2,ContigSplitRatio3,medianFN,medianFC,threshold);
	  fflush(stdout);
	}
	if(DEBUG) assert(start <= k && k <= end);

	/* check midpoint of local minima to see if a chimeric split is warranted */

	/* scan ahead to end of chimeric junction */
	int t, cend = (SplitRev >= 3 ? k+1 : end), cstart = (SplitRev >= 3 ? k : start);
	if(SplitRev <= 1){
	  for(t = k+1; t < R; t++)
	    if(contig->sitecntFN[0][t] >= threshold)
	      break;
	  cend = t-1;
	}

	/* scan back to start of chimeric junction */
	if(SplitRev <= 1){
	  for(t = k-1; t > L; t--)
	    if(contig->sitecntFN[0][t] >= threshold)
	      break;
	  cstart = t+1;
	}

	/* trim back and forward from ends of non-chimeric region based on EndTrimCov */
	int Lbreak = cstart-1, Rbreak = cend+1;
	while(Rbreak < right && contig->fragcnt[0][Rbreak - (FRAGCOV_FIX?1:0)] < EndTrimCov)
	  Rbreak++;
	while(Lbreak > left && contig->fragcnt[0][Lbreak-1] < EndTrimCov)
	  Lbreak--;
	  
	if(!SplitRev){
	  /* check if left or right region outside of the chimeric region is too small (less than 2 sites) ie to near the shallow ends */
	  int Lcnt = 0, Rcnt = 0;
	  for(t = left; t <= Lbreak; t++)
	    Lcnt += (peak[t] && contig->sitecntFN[0][t] >= threshold) ? 1 : 0;
	  for(t = Rbreak; t <= right; t++)
	    Rcnt += (peak[t] && contig->sitecntFN[0][t] >= threshold) ? 1 : 0;
	      
	  if(Lcnt < 2 || Rcnt < 2){
	    if(BVERB){
	      printf("      Lbreak=%d,Rbreak=%d,Lcnt=%d,Rcnt=%d: skipping this minima\n",Lbreak,Rbreak,Lcnt,Rcnt);
	      fflush(stdout);
	    }
	    continue;
	  }
	}

	int Ftype = -1, Fleft, Fright;
	if(ChimQuality >= 2 && SplitFragileD > 0.0){/* check if one sided Fragile site start is present between left and Lbreak (can end up R-1)  */
	  int t = left+1;
	  for(; t < Lbreak; t++){
	    //	    if(SplitSite && !peak[t]) continue;

	    double fragileLeft = contig->sitecntN4[0][t];
	    double fragileRight = contig->sitecntN5[0][t];
	    if(fragileLeft - fragileRight > SplitFragileD && fragileLeft > fragileRight * SplitFragileR){
	      Ftype = 1;
	      break;
	    }
	    if(fragileRight - fragileLeft > SplitFragileD && fragileRight > fragileLeft * SplitFragileR){
	      Ftype = 0;
	      break;
	    }
	  }

	  if(DEBUG) assert((Ftype >= 0) == (t < Lbreak));

	  if(t < Lbreak){
	    Fright = Fleft = t;
	    for(t++; t < R; t++){
	      //	      if(SplitSite && !peak[t]) continue;

	      double fragileLeft = contig->sitecntN4[0][t];
	      double fragileRight = contig->sitecntN5[0][t];
	      if(Ftype ? (fragileLeft - fragileRight > SplitFragileD && fragileLeft > fragileRight * SplitFragileR) :
		 (fragileRight - fragileLeft > SplitFragileD && fragileRight > fragileLeft * SplitFragileR) ){
		Fright = t;
		continue;
	      }

	      break;
	    }

	    int origk = k;

	    k = Lbreak = (Fleft + Fright)/2;
	    Rbreak = max(Fright, Lbreak+1);

	    if(VERB){
	      printf("Found %s Contig End at k=%d(range=%d..%d), site[k]=%0.3f:peak[k]=%d(%d..%d)sitecntFN[k]=%0.1f,segdupLeft=%0.1f,segdupRight=%0.1f,fragileLeft=%0.1f,%0.1f,fragileRight=%0.1f,%0.1f,Outlier=%0.1f,Norm=%0.1f,Sd=%0.1f,ChiSq=%0.2e\n",
		     Ftype ? "Left" : "Right", k,Fleft,Fright,Y[k], peak[k],peak[Fleft],peak[Fright],contig->sitecntFN[0][k],contig->sitecntN2[0][k],contig->sitecntN3[0][k],contig->sitecntN4[0][Fleft],contig->sitecntN4[0][Fright],
		     contig->sitecntN5[0][Fleft],contig->sitecntN5[0][Fright],contig->sitecntN6[0][k],contig->sitecntFNnorm[0][k],contig->fragSd[0][k]*1000.0,contig->fragChiSq[0][k]);
	      fflush(stdout);
	    }
	    if(BVERB>=2) k = origk;
	  }
	}

	if(VERB && Ftype < 0){
	  printf("Found Chimeric junction at k=%d(range=%d..%d),site[k]=%0.3f:sitecntFN[k]=%0.1f,fragcnt[k-1]=%0.1f,fragcnt[k]=%0.1f,leftFN=%0.1f(at site[%d]=%0.3f),rightFN=%0.1f(at site[%d]=%0.3f),minFN=%0.1f,ContigSplitRatio=%0.4f,%0.4f,%0.4f:Lbreak=%d,Rbreak=%d,medianFN=%0.1f,medianFC=%0.1f,threshold=%0.1f\n",
		 k,cstart,cend,Y[k],contig->sitecntFN[0][k],contig->fragcnt[0][k-1],contig->fragcnt[0][k],leftFN,kleft,Y[kleft],rightFN,kright,Y[kright],minFN,ContigSplitRatio,ContigSplitRatio2,ContigSplitRatio3,
		 Lbreak,Rbreak,medianFN,medianFC,threshold);
	  fflush(stdout);
	}
	if(BVERB>=2) {
	  if(VERB && Ftype < 0){
	    printf("\t Ignoring Chimeric junction due to BVERB=%d\n",BVERB);
	    fflush(stdout);
	  }
	  if(VERB && Ftype >= 0){
	    printf("\t Ignoring Fragile split junction due to BVERB=%d\n",BVERB);
	    fflush(stdout);
	  }
	  continue;
	}

	if(Lbreak > left){  /* output region between left and Lbreak as a seperate contig */
	  /* check that region has at least AlignedSiteThreshold sites */
	  int cnt = 0;
	  for(t = left; t <= Lbreak; t++)
	    cnt += (peak[t] ? 1 : 0);

	  if(cnt >= AlignedSiteThreshold && Y[Lbreak] - Y[left] >= MinSplitLen){
	    if(DEBUG) assert(ContigCntFile != 0);
	    long long oldCMapID = CMapID;
	    NumSplits++;
	    if(SplitCnt)
	      CMapID += NumSplits * SplitCnt;
	    else
	      CMapID = nextContigId(ContigCntFile);
	    if(CMapID < 0){
	      printf("CMapID = %lld wrapped around 63 bit integer limit (orig CMapID= %lld, NumSplits=%d, SplitCnt = %lld)\n",
		     CMapID, oldCMapID, NumSplits, SplitCnt);
	      fflush(stdout);
	      assert(CMapID > 0);
	    }
	    if(1 /* HERE BVERB */){
	      if(extend>=2){
		if(DEBUG) assert(Refine && contigid<=0);
		printf("%s_contig%lld:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f,sites=%d),MinTP=%d:original ends: left=%0.3f,right=%0.3f,NumSplits=%d,SplitCnt=%lld\n",
		       basename, CMapID,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],Lbreak,Y[Lbreak],Y[n+1], 
		       Y[Lbreak] - Y[left]+leftend+rightend, cnt, AlignedSiteThreshold,contig->left,contig->right,NumSplits,SplitCnt);
	      } else if(contigid<=0){
		if(DEBUG) assert(Refine);
		printf("%s_contig%lld:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f,sites=%d),MinTP=%d,MinSplitLen=%0.3f,NumSplits=%d,SplitCnt=%lld\n",
		       basename,CMapID,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],Lbreak,Y[Lbreak],Y[n+1], Y[Lbreak]-Y[left]+leftend+rightend,
		       cnt,AlignedSiteThreshold,MinSplitLen,NumSplits,SplitCnt);
	      }
	      fflush(stdout);
	    }

	    // HERE HERE : If Ftype >= 0 && SplitFragileE : mark Fragile ends as non-extendable

	    output_draftIO(contig,contigid,basename,left,Lbreak,n,peak, Y, Allele, newleft, newright);

	    CMapID = oldCMapID;
	  } else if(1/* BVERB */){
	    if(extend>=2){
	      if(DEBUG) assert(Refine && contigid <= 0);
	      printf("%s_contigN_refined:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f):original ends: left=%0.3f,right=%0.3f\n",
		     basename,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],Lbreak,Y[Lbreak],Y[n+1], Y[Lbreak]-Y[left]+leftend+rightend, contig->left,contig->right);
	    } else if(contigid <= 0){
	      if(DEBUG) assert(Refine);
	      printf("%s_contigN_refined:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f),TP=%d,MinTP=%d,MinSplitLen=%0.3f\n",
		     basename,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],Lbreak,Y[Lbreak],Y[n+1], Y[Lbreak]-Y[left]+leftend+rightend,
		     cnt,AlignedSiteThreshold,MinSplitLen);
	    }
	    fflush(stdout);
	  }
	}

	/* truncate the remaining contig to between Rbreak and right */
	L = left = Rbreak;
	if(SplitRev <= 1)
	  while(L < R && contig->sitecntFN[0][L] < threshold)
	    L++;
	k = end = L;/* continue scanning remaining contig from L+1 */
      } /* advance k to next minima */

      if(ChimQuality >= 2 && SplitFragileD > 0.0){/* check if one sided Fragile site is present between left and right */
	int origleft = left;

	for(int k = left + 1; k < right; k = left + 1){
	  int Ftype = -1, Fleft, Fright, Lbreak = left,Rbreak;

	  int t = left+1;
	  for(; t < right; t++){
	    //	    if(SplitSite && !peak[t]) continue;

	    double fragileLeft = contig->sitecntN4[0][t];
	    double fragileRight = contig->sitecntN5[0][t];
	    if(fragileLeft - fragileRight > SplitFragileD && fragileLeft > fragileRight * SplitFragileR){
	      Ftype = 1;
	      break;
	    }
	    if(fragileRight - fragileLeft > SplitFragileD && fragileRight > fragileLeft * SplitFragileR){
	      Ftype = 0;
	      break;
	    }
	  }

	  if(DEBUG) assert((Ftype >= 0) == (t < right));

	  if(t >= right)
	    break;

	  Fright = Fleft = t;
	  for(t++; t < right; t++){
	    //	    if(SplitSite && !peak[t])	      continue;
	    double fragileLeft = contig->sitecntN4[0][t];
	    double fragileRight = contig->sitecntN5[0][t];
	    if(Ftype ? (fragileLeft - fragileRight > SplitFragileD && fragileLeft > fragileRight * SplitFragileR) :
	       (fragileRight - fragileLeft > SplitFragileD && fragileRight > fragileLeft * SplitFragileR) ){
	      Fright = t;
	      continue;
	    }
	    break;
	  }
	  
	  Lbreak = k = (Fleft + Fright)/2;
	  if(DEBUG) assert(Lbreak > left);
	  Rbreak = max(Fright, Lbreak+1);
	  if(VERB){
	    printf("Found %s Contig End at k=%d(range=%d..%d), site[k]=%0.3f:sitecntFN[k]=%0.1f,segdupLeft=%0.1f,segdupRight=%0.1f,fragileLeft=%0.1f,%0.1f,fragileRight=%0.1f,%0.1f,Outlier=%0.1f,Norm=%0.1f,Sd=%0.1f,expSd=%0.1f,Cov=%0.1f,ChiSq=%0.2e\n",
		   Ftype ? "Left" : "Right", k,Fleft,Fright,Y[k], contig->sitecntFN[0][k],contig->sitecntN2[0][k],contig->sitecntN3[0][k],
		   contig->sitecntN4[0][Fleft],contig->sitecntN4[0][Fright],contig->sitecntN5[0][Fleft],contig->sitecntN5[0][Fright],contig->sitecntN6[0][k],contig->sitecntFNnorm[0][k],contig->fragSd[0][k]*1000.0,contig->expSd[0][k]*1000.0, contig->fragCov[0][k],contig->fragChiSq[0][k]);
	    fflush(stdout);
	  }

	  /* check that region has at least AlignedSiteThreshold sites */
	  int cnt = 0;
	  for(t = left; t <= Lbreak; t++)
	    cnt += (peak[t] ? 1 : 0);

	  if(cnt >= AlignedSiteThreshold && Y[Lbreak] - Y[left] >= MinSplitLen && BVERB <= 1){
	    if(DEBUG) assert(ContigCntFile != 0);
	    long long oldCMapID = CMapID;
	    NumSplits++;
	    if(SplitCnt)
	      CMapID += NumSplits * SplitCnt;
	    else
	      CMapID = nextContigId(ContigCntFile);
	    if(CMapID < 0){
	      printf("CMapID = %lld wrapped around 63 bit integer limit (orig CMapID= %lld, NumSplits=%d, SplitCnt = %lld)\n",
		     CMapID, oldCMapID, NumSplits, SplitCnt);
	      fflush(stdout);
	      assert(CMapID > 0);
	    }
	    if(1 /* HERE BVERB */){
	      if(extend>=2){
		if(DEBUG) assert(Refine && contigid<=0);
		printf("%s_contig%lld:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f,sites=%d),MinTP=%d:original ends: left=%0.3f,right=%0.3f,NumSplits=%d,SplitCnt=%lld\n",
		       basename, CMapID,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],Lbreak,Y[Lbreak],Y[n+1], 
		       Y[Lbreak] - Y[left]+leftend+rightend, cnt, AlignedSiteThreshold,contig->left,contig->right,NumSplits,SplitCnt);
	      } else if(contigid<=0){
		if(DEBUG) assert(Refine);
		printf("%s_contig%lld:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f,sites=%d),MinTP=%d,MinSplitLen=%0.3f,NumSplits=%d,SplitCnt=%lld\n",
		       basename,CMapID,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],Lbreak,Y[Lbreak],Y[n+1], Y[Lbreak]-Y[left]+leftend+rightend,
		       cnt,AlignedSiteThreshold,MinSplitLen,NumSplits,SplitCnt);
	      }
	      fflush(stdout);
	    }

	    // HERE HERE : If Ftype >= 0 && SplitFragileE : mark Fragile ends as non-extendable

	    output_draftIO(contig,contigid,basename,left,Lbreak,n,peak, Y, Allele, newleft, newright);

	    CMapID = oldCMapID;
	  } else if(1/* BVERB */){
	    if(extend>=2){
	      if(DEBUG) assert(Refine && contigid <= 0);
	      printf("%s_contigN_refined:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f):original ends: left=%0.3f,right=%0.3f\n",
		     basename,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],Lbreak,Y[Lbreak],Y[n+1], Y[Lbreak]-Y[left]+leftend+rightend, contig->left,contig->right);
	    } else if(contigid <= 0){
	      if(DEBUG) assert(Refine);
	      printf("%s_contigN_refined:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f),TP=%d,MinTP=%d,MinSplitLen=%0.3f\n",
		     basename,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],Lbreak,Y[Lbreak],Y[n+1], Y[Lbreak]-Y[left]+leftend+rightend,
		     cnt,AlignedSiteThreshold,MinSplitLen);
	    }
	    fflush(stdout);
	  }

	  /* truncate the remaining contig to between Rbreak and right */
	  left = Rbreak;/* continue scanning remaining contig from left + 1 */
	}
	if(BVERB>=2) left = origleft;
      }
    }
  } // if(Refine && ContigSplitRatio > 0.0 ... )

  /* check that (remaining) region Hcuts[left .. right] has at least AlignedSitedThreshold sites */
  int cnt = 0;
  for(register int t = left; t <= right; t++)
    cnt += (peak[t] ? 1 : 0);
  if(1 /* BVERB */){
    if(extend>=2){
      if(DEBUG) assert(Refine && contigid <= 0 && CMapID > 0);
      printf("%s_contig%lld:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f):original ends: left=%0.3f,right=%0.3f,MinSplitLen=%0.3f,cnt=%d,%d\n",
	     basename,CMapID, n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],right,Y[right],Y[n+1], Y[right]-Y[left]+leftend+rightend, 
	     contig->left,contig->right,MinSplitLen,cnt,AlignedSiteThreshold);
    } else if(contigid <= 0){
      if(DEBUG) assert(Refine);
      if(CMapID > 0)
	printf("%s_contig%lld:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f),TP=%d,MinTP=%d,MinSplitLen=%0.3f,cnt=%d,%d\n",
	       basename,CMapID,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],right,Y[right],Y[n+1], Y[right]-Y[left]+leftend+rightend,
	       cnt,AlignedSiteThreshold,MinSplitLen,cnt,AlignedSiteThreshold);
      else
	printf("%s.cmap:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f),Raw Length=%0.3f,Trimmed Length=%0.3f),TP=%d,MinTP=%d,MinSplitLen=%0.3f,cnt=%d,%d\n",
	       basename,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],right,Y[right],Y[n+1], Y[right]-Y[left]+leftend+rightend,
	       cnt,AlignedSiteThreshold,MinSplitLen,cnt,AlignedSiteThreshold);
    } else 
      printf("%s_contig%lld.cmap:numsite=%d(tot=%d),maps=%d,(EndTrimCov=%0.2f,left=%d(%0.3f),right=%d(%0.3f), Raw Length=%0.3f,Trimmed Length=%0.3f),TP=%d,MinTP=%d,MinSplitLen=%0.3f,cnt=%d,%d\n",
    	     basename,contigid,n,contig->totsites(),contig->nummaps,EndTrimCov,left,Y[left],right,Y[right],Y[n+1], Y[right]-Y[left]+leftend+rightend,
	     cnt,AlignedSiteThreshold,MinSplitLen,cnt,AlignedSiteThreshold);
    fflush(stdout);
  }
  if(/*ContigSplitRatio <= 0.0 || */(cnt >= AlignedSiteThreshold && Y[right] - Y[left] >= MinSplitLen)){
    output_draftIO(contig,contigid,basename,left,right,n,peak,Y, Allele, newleft, newright);
  } else if(!strstr(basename,"/dev/null")){/* create an empty .cmap file */

    long long origCMapID = CMapID;
    if(CMapID > 0 && HapSitePvalue > 0.0)
      CMapID = CMapID * 10 + Allele;

    char filename[PATH_MAX];
    strcpy(filename,basename);
    int i = strlen(filename);
    if(contigid > 0)
      sprintf(&filename[i],"_contig%lld.cmap", contigid);
    else if(CMapID > 0){
      if(contig->HapSite[0])
	sprintf(&filename[i],"_contig%lld.hmap",CMapID);
      else
	sprintf(&filename[i],"_contig%lld.cmap",CMapID);
    } else
      sprintf(&filename[i],".cmap");

    if(checkFile(filename)){
      CMapID = origCMapID;
      return;
    }

    if(VERB){
      if(contigid > 0){
	if(Refine)
	  printf("Generating empty %s (Assembler Refined contig%lld in .cmap format)\n",filename, contigid);
	else
	  printf("Generating empty %s (draft contig%lld in .cmap format)\n",filename, contigid);
      } else if(CMapID > 0){
	if(Allele)
	  printf("Generating empty %s (refined contig%lld Allele %d in .cmap format)\n",filename, origCMapID,Allele);
	else if(contig->HapSite[0])
	  printf("Generating empty %s (haplotyped refined contig%lld in .hmap format)\n",filename,CMapID);
	else
	  printf("Generating empty %s (refined contig%lld in .cmap format)\n",filename,CMapID);
      } else
	printf("Generating empty %s (refined contig%lld in .cmap format)\n",filename, CMapID);
      fflush(stdout);
    }

    FILE *fp;
    if((fp = fopen(filename,"w"))==NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("failed to open file %s for writing contig%lld draft consensus(as .cmap):errno=%d:%s\n",filename, contigid ? contigid : CMapID,eno,err);
      fflush(stdout);exit(1);
    }
    FILEclose(fp);

    CMapID = origCMapID;
  }

  delete [] sitecnt;
  delete [] smoothcnt;
  delete [] peak;
}

/* output maps in .hmap format based on maps[i]->contig (ignores maps[i]->site etc ) */
/* use .hmap suffix if hmapsuffix = 1, otherwise use .cmap suffix for output filename */
/* this function is only called from pairalign() when -pairmergeHmap is used */
void output_hmap(char *basename, Cmap *maps[], int start, int end, int hmapsuffix)
{
  if(end > start){
    printf("output_hmap(): Only supports a single contig output: start=%d, end=%d\n", start, end);
    fflush(stdout); exit(1);
  }

  if(DEBUG) assert(maps[start]->contig && maps[start]->contig->HapSite[0]);
  if(DEBUG) assert(maps[start]->id == maps[start]->contig->id);

  Ccontig *contig = maps[start]->contig;
  int *HapSite = contig->HapSite[0];

  if(VERB>=2 && contig->id == 2875){
    printf("output_hmap:id=%lld,N=%d,Len=%0.3f\n",contig->id,contig->numsite[0], contig->site[0][contig->numsite[0]+1]);
    fflush(stdout);
  }

  if(strstr(basename,"/dev/null"))
    return;

  char filename[PATH_MAX];
  FILE *fp;
  if(hmapsuffix)
    sprintf(filename,"%s.hmap",basename);
  else
    sprintf(filename,"%s.cmap",basename);

  if(checkFile(filename))
    return;
  
  if(VERB){
    printf("Generating %s (Haplotyped map produced by merging 2 or more Cmaps): id=%lld, N=%d Labels, Len= %0.3f kb\n",filename,contig->id,contig->numsite[0],contig->site[0][contig->numsite[0]+1]);
    fflush(stdout);
  }

  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("failed to open file %s for writing Haplotype Map:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  int n = contig->numsite[0];
  if(DEBUG>=1+RELEASE) assert(HapSite[0]==0 && HapSite[n+1]==0);

  int NumSites = 0;
  for(int k = 1; k <= n; k++)
    if(HapSite[k])
      NumSites++;
  if(DEBUG && PairMergeHmap > 0) assert(NumSites == n);
  
  size_t blen = 1024*1024;/* must be less than scif_io buffer size */
  char *write_buffer = (char *)malloc(blen);
  if(!write_buffer){
    printf("output_draftIO:malloc(%llu) failed\n",(unsigned long long)blen);
    exit(1);
  }
  setbuffer(fp, write_buffer, blen);

  /* write out commandline */
  printversion(fp);
  
  /* output key header lines */
  fprintf(fp,"# CMAP File Version:\t0.1\n");
  fprintf(fp,"# Label Channels:\t%d\n", colors);
  if(nummaps > 0 && Gmap[0]){
    for(register int c = 0; c < colors; c++)
      if(Gmap[0]->Nickase[c])
	fprintf(fp,"# Nickase Recognition Site %d:\t%s\n", c+1, Gmap[0]->Nickase[c]);
      else
	fprintf(fp,"# Nickase Recognition Site %d:\tunknown\n", c+1);
  }
  fprintf(fp,"# Number of Consensus Maps:\t1\n");

  fprintf(fp,"# StdDev & Coverage refer to the interval between the current site and the next site\n");
  if(!contig->fragcnt[0])
    fprintf(fp,"# Coverage currently not computed and set to 0\n");
  if(contig->HapSite[0] && contig->HapDelta[0])
    fprintf(fp,"# HapDelta, HapDeltaPhase and HapDeltaScore refer to the interval between the current site and the previous site\n");

  fprintf(fp,"#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence");
  if(VERB>=3 && contig->id == 125716166596LL){
    printf("output_hmap():contig->id=%lld, MaskL = %lx, MaskR = %lx\n",contig->id, contig->MaskL, contig->MaskR);
    fflush(stdout);
  }
  if(contig->MaskL || contig->MaskR)
    fprintf(fp,"\tMask");
  if(contig->HapSite[0])
    fprintf(fp,"\tHapSite\tHapDelta\tSitePhaseConf\tDeltaPhaseConf\tHapSiteScore\tHapDeltaScore");
  fprintf(fp,"\n");

  fprintf(fp,"#f int\tfloat\tint\tint\tint\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat");
  if(contig->MaskL || contig->MaskR)
    fprintf(fp,"\tHex");
  if(contig->HapSite[0])
    fprintf(fp,"\tint\tfloat\tfloat\tfloat\tfloat\tfloat");
  fprintf(fp,"\n");

  double invLog10 = 1.0/log(10.0);

  double *Y = contig->site[0];
  int SiteID = 0;
  int prevk = 0;

  for(int k = 1; k <= n; k++){
    if(HapSite[k]){
      int nextk = k+1;
      for(; nextk <= n; nextk++)
	if(HapSite[k])
	  break;

      if(DEBUG) assert(Y[k] >= Y[prevk]);
      if(DEBUG && nextk <= n) assert(Y[nextk] >= Y[k]);

      ++SiteID;
      float cov = contig->fragcnt[0] ? contig->fragcnt[0][k] : 0.0;/* fragcnt not yet implemented in merge_contigs() */
      float occ = contig->sitecnt[0] ? (COVFIX ? contig->sitecnt[0][k] : min(contig->sitecnt[0][k], cov)) : 0.0;/* sitecnt not yet implemented in merge_contigs() */
      if(DEBUG) assert(isfinite(cov));
      if(DEBUG) assert(isfinite(occ));
      double SDcov = min(16.0f,max(1.0f,cov));
      double StdDev = 0.0;
      if(nextk <= n && HapSite[nextk])
	StdDev = 1000.0*sqrt((SF[0]*SF[0]+SD[0]*SD[0]*(Y[nextk]-Y[k]))/SDcov)/(1.0-FN[0]);// HERE HERE : apply QUADRATIC_VARIANCE
      if(DEBUG) assert(isfinite(StdDev));
      
      fprintf(fp,"%lld\t%0.1f\t%d\t%d\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.1f",
	      contig->id, // CMapID
	      Y[n+1]*1000.0, // ContigLength
	      NumSites, // NumSites
	      SiteID, // SiteID
	      1, // LabelChannel
	      Y[k]*1000.0, // Position
	      StdDev,// StdDev
	      cov, // Coverage
	      occ  // Occurrence
	      );

      if(contig->MaskL || contig->MaskR)
	fprintf(fp,"\t%lx", (SiteID==1) ? contig->MaskL : (size_t)0);

      double HapDeltaScore =  contig->HapDeltaScore[0] ? contig->HapDeltaScore[0][prevk] : 0.0;
      double HapDelta = contig->HapDelta[0][prevk];
      for(int t = prevk+1; t < k; t++){
	HapDelta += contig->HapDelta[0][t];
	if(contig->HapDeltaScore[0])
	  HapDeltaScore = max(HapDeltaScore,contig->HapDeltaScore[0][t]);
      }
      HapDeltaScore = (fabs(HapDelta) > 0.001) ? HapDeltaScore * invLog10 : 0.0;
      
      double HapSiteScore = contig->HapSiteScore[0] ? contig->HapSiteScore[0][k] * invLog10 : 0.0;
      double HapSitePhase = 0.0, HapDeltaPhase = 0.0;
      if(contig->HapSitePhase[0])
	HapSitePhase = (1 <= HapSite[k] && HapSite[k] <= 2) ? contig->HapSitePhase[0][k]*invLog10 : -1.0;
      if(contig->HapDeltaPhase[0])
	HapDeltaPhase = (fabs(HapDelta) > 0.001) ? contig->HapDeltaPhase[0][prevk]*invLog10 : -1.0;
      
      if(VERB>=2 && contig->id == 2875){
	printf("ID=%lld,Len=%0.4f,NumSites=%d,SiteID=%d,Position=%0.4f,HapSite=%d,HapDelta=%0.4f (k=%d,prevk=%d,n=%d)\n",contig->id,Y[n+1],NumSites,SiteID,Y[k],HapSite[k],HapDelta,k,prevk,n);
	fflush(stdout);
      }

      fprintf(fp,"\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.2f\t%0.2f", HapSite[k], HapDelta * 1000.0, HapSitePhase, HapDeltaPhase, HapSiteScore, HapDeltaScore);

      fprintf(fp,"\n");

      prevk = k;
    }
  }

  /* right end Y[n+1] */
  fprintf(fp,"%lld\t%0.1f\t%d\t%d\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.1f",
	  contig->id, // CMapID
	  Y[n+1]*1000.0, // ContigLength
	  NumSites, // NumSites
	  NumSites + 1, // SiteID = right end
	  0, // LabelChannel
	  Y[n+1]*1000.0, // Position
	  0.0, // StdDev
	  1.0 /* cov */, // Coverage
	  0.0/* WAS 1.0 */   // Occurrence
	  );

  if(contig->MaskL || contig->MaskR)
    fprintf(fp,"\t%lx", contig->MaskR);

  if(contig->HapSite[0]){
    double HapDeltaScore = contig->HapDeltaScore[0] ? contig->HapDeltaScore[0][prevk] : 0.0;
    double HapDelta = contig->HapDelta[0][prevk];
    for(int t = prevk+1; t <= n; t++){
      HapDelta += contig->HapDelta[0][t];
      if(contig->HapDeltaScore[0])
	HapDeltaScore = max(HapDeltaScore,contig->HapDeltaScore[0][t]);
    }
    HapDeltaScore = (fabs(HapDelta) > 0.001) ? HapDeltaScore * invLog10 : 0.0;
      
    double HapDeltaPhase = 0.0;
    if(contig->HapDeltaPhase[0])
      HapDeltaPhase = (fabs(HapDelta) > 0.001) ? contig->HapDeltaPhase[0][prevk]*invLog10 : -1.0;

    if(VERB>=2 && contig->id == 2875){
      printf("ID=%lld,Len=%0.4f,NumSites=%d,SiteID=%d,Position=%0.4f,HapDelta=%0.4f (prevk=%d,n=%d):Right End\n",contig->id,Y[n+1],NumSites,NumSites+1,Y[n+1],HapDelta,prevk,n);
      fflush(stdout);
    }

    fprintf(fp,"\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.2f\t%0.2f", 0 , HapDelta * 1000.0,-1.0,HapDeltaPhase,-1.0,HapDeltaScore);
  }

  fprintf(fp,"\n");

  FILEclose(fp);

  free(write_buffer);
}

/* output maps in .cmap format (see also older output_csites() in output.cpp, only used for output of _r.cmap with coverage profile computed from alignments) */
void output_cmap(char *basename, Cmap *maps[], int start, int end)
{
  if(strstr(basename,"/dev/null"))
    return;

  char filename[PATH_MAX];
  sprintf(filename,"%s.cmap",basename);
  
  if(checkFile(filename))
    return;

  if(VERB && !PairMerge){
    #pragma omp critical
    {
      printf("Generating %s (with %d contigs or maps)\n",filename, end - start + 1);
      fflush(stdout);
    }
  }

  char *origfilename = strdup(filename);
  if(1){
    int len = strlen(filename);
    sprintf(&filename[len],".tmp");
  }
  char *tmpfilename = strdup(filename);

  FILE *fp;
  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("failed to open file %s for writing .cmap file of %d maps: errno=%d: %s\n",filename, end - start + 1, eno, err);
    fflush(stdout);exit(1);
  }
  
  size_t blen = 10*1024*1024;/* matches scif_io buffer size */
  char *write_buffer = (char *)malloc(blen);
  if(!write_buffer){
    printf("output_cmap:malloc(%llu) failed\n",(unsigned long long)blen);
    fflush(stdout);exit(1);
  }
  setbuffer(fp, write_buffer, blen);

  /* write out commandline */
  printversion(fp);

  /* output key header lines */
  if(CmapChimQuality >= 2)
    fprintf(fp,"# CMAP File Version:\t0.2\n");
  else
    fprintf(fp,"# CMAP File Version:\t0.1\n");
  fprintf(fp,"# Label Channels:\t%d\n", colors);
  if(end >= start && maps[start]){
    for(int c = 0; c < colors; c++)
      if(maps[start]->Nickase[c]){
	int L = strlen(maps[start]->Nickase[c]);
	while(L > 0 && maps[start]->Nickase[c][L-1] == '\n'){
	  //	  printf("maps[start]->Nickase[%d] has newline at end: removed\n",c);
	  maps[start]->Nickase[c][--L] = '\0';
	}
	fprintf(fp,"# Nickase Recognition Site %d:\t%s\n", c+1, maps[start]->Nickase[c]);
      }
  } else
    for(int c = 0; c < colors; c++)
      fprintf(fp,"# Nickase Recognition Site %d:\t%s\n", c+1, "unknown");

  /* count number of maps excluding map splits */
  int mapcnt = 0;
  for(int i = start; i <= end; i++)
    if(!maps[i]->origmap)
      mapcnt++;

  if(VERB>=2){
    printf("start=%d,end=%d,mapcnt=%d\n",start,end,mapcnt);
    fflush(stdout);
  }

  int SNRpresent = 0;
  if(CmapSNR){
    for(int i = start; i <= end; i++)
      if(maps[i]->SNRcnt[0]){
	SNRpresent = 1;
	break;
      }
  }
  // if SNRgmean[] & lnSNRsd[] (but NOT SNRcnt) are present, output just SNRgmean and lnSNRsd */

  int SNRval = 0;
  for(int i = start; i <= end; i++)
    if(maps[i]->SNR[0]){
      SNRval = 1;
      break;
    }

  int MaskPresent = (CmapMask > 0) ? 1 : 0;/* output Mask if it is present in some input CMAP or is present in any output CMAP in current output file */
  for(int i = start; i <= end; i++)
    if(maps[i]->Mask[0]){
      MaskPresent = 1;
      break;
    }

  int HapSitePresent = 0;
  for(int i = start; i <= end; i++)
    if(maps[i]->contig && maps[i]->contig->HapSite[0]){
      HapSitePresent = 1;
      break;
    }

  int FragSdPresent = 0;
  for(int i = start; i <= end; i++)
    if(maps[i]->FragSd[0] && maps[i]->ExpSd[0] && maps[i]->FragCov[0] && maps[i]->FragChiSq[0]){
      if(VERB>=2){
	printf("output_cmap:Found FragSd,ExpSd,FragCov,FragChiSq in Gmap[%d]:filename= %s\n",i,filename);
	fflush(stdout);
      }
      FragSdPresent = 1;
      break;
    }

  if(QueryOrigin)
    fprintf(fp,"# Original Query Maps from: %s\n", QueryOrigin);
  if(RefOrigin && !QueryOrigin && num_files <= 0)
    fprintf(fp,"# Original Reference Maps from: %s\n", RefOrigin);

  const char *OutlierString = (CmapChimQuality >= 2) ? (FragSdPresent ? " & OutlierFrac & MolSd & MolCov & ChiSq" : " & OutlierFrac") : 
    FragSdPresent ? " & MolSd & MolCov & ChiSq" : "";// NEW372
  fprintf(fp,"# Number of Consensus Maps:\t%d\n", mapcnt);
  if(FRAGCOV_FIX)
    fprintf(fp,"# StdDev & Coverage%s refer to the interval between current site and next site\n",OutlierString);
  else
    fprintf(fp,"# StdDev%s refer%s to the interval between current site and next site\n",OutlierString,CmapChimQuality >= 2 ? "" : "s");
  if(HapSitePresent)
    fprintf(fp,"# HapDelta, HapDeltaPhase and HapDeltaScore refer to the interval between the current site and the previous site\n");

  fprintf(fp,"#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence");

  if(CmapChimQuality >= 1 && !HapSitePresent){
    fprintf(fp,"\tChimQuality");
    if(CmapChimQuality >= 2){
      fprintf(fp,"\tSegDupL\tSegDupR\tFragileL\tFragileR\tOutlierFrac\tChimNorm");
    }
  }
  if(MaskPresent)
    fprintf(fp,"\tMask");
  if(FragSdPresent)
    fprintf(fp,"\tMolSd\tExpSd\tMolCov\tChiSq");
  if(HapSitePresent)
    fprintf(fp,"\tHapSite\tHapDelta\tSitePhaseConf\tDeltaPhaseConf\tHapSiteScore\tHapDeltaScore");    
  if(Tracksites)
    fprintf(fp,"\tOriginLeft\tOriginRight");
  if(SNRpresent)
    fprintf(fp,"\tGmeanSNR\tlnSNRsd\tSNRcount\tSNR ...");
  else if(SNRval)
    fprintf(fp,"\tGmeanSNR\tlnSNRsd");
  fprintf(fp,"\n");

  if(floatcov)
    fprintf(fp,"#f int\tfloat\tint\tint\tint\tfloat\tfloat\tfloat\tfloat");
  else
    fprintf(fp,"#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint");
  if(CmapChimQuality >= 1 && !HapSitePresent){
    fprintf(fp,"\tfloat");
    if(CmapChimQuality >= 2){
      fprintf(fp,"\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat");
    }
  }
  if(MaskPresent)
    fprintf(fp,"\tHex");
  if(FragSdPresent)
    fprintf(fp,"\tfloat\tfloat\tfloat\tfloat");
  if(HapSitePresent)
    fprintf(fp,"\tint\tfloat\tfloat\tfloat\tfloat\tfloat");
  if(Tracksites)
    fprintf(fp,"\tint\tint");
  if(SNRpresent)
    fprintf(fp,"\tfloat\tfloat\tint\tfloat ...");
  else if(SNRval)
    fprintf(fp,"\tfloat\tfloat");
  fprintf(fp,"\n");

  for(int i = start; i <= end; i++){
    Cmap *pmap = maps[i];
    if(pmap->origmap)
      continue;/* don't output map fragments, since they have the same ID as the complete map */

    if(VERB>=3 && pmap->id == 125716166596LL){
      int N = pmap->numsite[0];
      printf("output_cmap():pmap->id=%lld, N=%d, MaskL = %lx, MaskR = %lx\n",
	     pmap->id, N, pmap->Mask[0] ? pmap->Mask[0][1] : (size_t)0, pmap->Mask[0] ? pmap->Mask[0][N+1] : (size_t)0);
      fflush(stdout);
    }


    if(DEBUG)
      for(int c = 0; c < colors; c++)
	if(pmap->siteSD[c]){
	  assert(pmap->sitecov[c]);
	  assert(pmap->sitecnt[c]);
	}

    int N = 0;
    for(int c = 0; c < colors; c++)
      N += pmap->numsite[c];
    if(VERB>=3 && pmap->id == 438){
      printf("output_cmap:i=%d,pmap->id=%lld, N=%d: contig=%p\n", i, pmap->id, N, pmap->contig);
      if(pmap->contig && pmap->contig->HapSite[0] && pmap->contig->HapDelta[0])
	for(int i = 1; i <= N+1; i++)
	  printf("\t i=%d: HapSite[i]=%d,HapDelta[i-1]=%0.4f\n",i,pmap->contig->HapSite[0][i],pmap->contig->HapDelta[0][i-1]);
      fflush(stdout);
    }

    Cmerge *pn = new Cmerge[N+1];
    int SiteID = 0;
    for(int c = 0; c < colors; c++){
      for(int j = 1;  j <= pmap->numsite[c]; j++){
	if(DEBUG/* HERE >=2 */ && !(pmap->site[c][j+1] >= pmap->site[c][j])){
	  printf("pmap->id=%lld,OriginalMoleculeId=%d,ChipId=%s,FlowCell=%d,ScanNumber=%d,RunIndex=%d:c=%d:numsite=%d,j=%d,site[c][j]= %0.10f,site[c][j+1]= %0.10f\n",
		 pmap->id,pmap->OriginalMoleculeId,pmap->ChipId,pmap->FlowCell,pmap->ScanNumber,pmap->RunIndex,c,pmap->numsite[c],j,pmap->site[c][j],pmap->site[c][j+1]);
	  for(int t = 1; t <= pmap->numsite[c]+1; t++)
	    printf("  site[c][%d]= %0.6f\n", t, pmap->site[c][t]);
	  fflush(stdout);
	  assert(pmap->site[c][j+1] >= pmap->site[c][j]);
	}
	pn[SiteID].color = c;
	pn[SiteID].newid = j;
	pn[SiteID].site = pmap->site[c][j];
	SiteID++;
      }
    }
    if(DEBUG) assert(SiteID == N);
    qsort(pn, N, sizeof(Cmerge), (intcmp *)siteinc);    
    if(DEBUG/* HERE >=2 */)
      for(int i = 0; i < N-1;i++)
	if(!(pn[i].site <= pn[i+1].site)){
	  printf("i=%d,N=%d:pn[i].site=%0.6f (color=%d,siteid=%d), pn[i+1].site=%0.6f(color=%d,siteid=%d)\n",
		 i,N,pn[i].site,pn[i].color,pn[i].newid, pn[i+1].site,pn[i+1].color, pn[i].newid);
	  fflush(stdout);
	  assert(pn[i].site <= pn[i+1].site);
	}

    /* Add right end */
    pn[N].color = pn[N].newid = -1;
    pn[N].site = pmap->site[0][pmap->numsite[0]+1];
    if((DEBUG && N > 0 && !(pn[N].site >= pn[N-1].site)) || (VERB>=2 && pmap->id==8)){
      printf("i=%d(%d..%d):colors=%d:pmap->id=%lld:\n",i,start,end,colors,pmap->id);
      for(int c = 0; c < colors; c++){
	printf("  c=%d:pmap->numsite[c]=%d,pmap->site[c][pmap->numsite[c]+1] = %0.6f\n",
	       c,pmap->numsite[c],pmap->site[c][pmap->numsite[c]+1]);
	for(int t = 1; t <= pmap->numsite[c]+1; t++)
	  printf("  pmap->site[c][%d] = %0.6f\n", t,pmap->site[c][t]);
      }
      printf("N=%d:pn[N-1].site=%0.6f (color=%d, siteid=%d), numsite[0]=%d, site[0][numsite[0]+1]= pn[N].site = %0.6f\n",
	     N, pn[N-1].site, pn[N-1].color, pn[N-1].newid, pmap->numsite[0], pmap->site[0][pmap->numsite[0]+1]);
      for(int t = 0; t <= N; t++)
	printf("   pn[%d]:site=%0.6f, color=%d\n",t, pn[t].site,pn[t].color);
      printf("   pmap->site[%d][%d]=%0.6f\n",pn[N-1].color,pn[N-1].newid,pmap->site[pn[N-1].color][pn[N-1].newid]);
      fflush(stdout);
      assert(pn[N].site >= pn[N-1].site);
    }
    if(DEBUG)
      for(int c = 1; c < colors; c++){
	if(DEBUG && !(fabs(pmap->site[c-1][pmap->numsite[c-1]+1] - pmap->site[c][pmap->numsite[c]+1]) < 1e-6)){
	  printf("pmap->id=%lld:c=%d,numsite[c-1]=%d,numsite[c]=%d,site[0][1]=%0.6f,len[0]=%0.6f,site[1][1]=%0.6f,len[1]=%0.6f\n",
		 pmap->id,c,pmap->numsite[c-1],pmap->numsite[c],pmap->site[0][1],pmap->site[0][pmap->numsite[0]+1],pmap->site[1][1],pmap->site[1][pmap->numsite[1]+1]);
	  fflush(stdout);
	  assert(fabs(pmap->site[c-1][pmap->numsite[c-1]+1] - pmap->site[c][pmap->numsite[c]+1]) < 1e-6);
	}
      }

    double invLog10 = 1.0/log(10.0);

    SiteID = 0;
    for(int k = 1; k <= N+1; k++){
      int c = pn[k-1].color;
      int newid = pn[k-1].newid;
      ++SiteID;
      double Coverage = (k <= N && pmap->siteSD[c]) ? pmap->sitecov[c][newid] : /* NEW372 */(FRAGCOV_FIX ? (k <= N /* && newid < pmap->numsite[c]*/) : (k <= N+1)) ? 1.0 : 0.0;
      double Occurrence = (k <= N && pmap->siteSD[c]) ? pmap->sitecnt[c][newid] : /* NEW372 */(k <= N) ? 1.0 : 0.0;
      if(floatcov)
	fprintf(fp, "%lld\t%0.1f\t%d\t%d\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.1f",
		pmap->id, // CMapID
		pmap->site[0][pmap->numsite[0]+1] * 1000.0, //ContigLength
		N, // NumSites
		SiteID, //SiteID
		c + 1, // LabelChannel
		pn[k-1].site * 1000.0, // Position
		(k <= N && pmap->siteSD[c]) ? pmap->siteSD[c][newid] * 1000.0 : 0.0, // StdDev
		Coverage,
		Occurrence);
      else
	fprintf(fp, "%lld\t%0.1f\t%d\t%d\t%d\t%0.1f\t%0.1f\t%d\t%d", 
		pmap->id, // CMapID
		pmap->site[0][pmap->numsite[0]+1] * 1000.0, //ContigLength
		N, // NumSites
		SiteID, //SiteID
		c + 1, // LabelChannel
		pn[k-1].site * 1000.0, // Position
		(k <= N && pmap->siteSD[c]) ? pmap->siteSD[c][newid] * 1000.0 : 0.0, // StdDev
		(int)floor(Coverage + 0.5),
		(int)floor(Occurrence + 0.5));

      if(CmapChimQuality >= 1 && !HapSitePresent){
	if(pmap->ChimQuality[c]){
	  if(DEBUG && k <= N) assert(isfinite(pmap->ChimQuality[c][newid]) && -1.0 <= pmap->ChimQuality[c][newid] && pmap->ChimQuality[c][newid] <= 100.001f);
	  fprintf(fp, "\t%0.2f", (k <= N) ? pmap->ChimQuality[c][newid] : 0.0);
	} else
	  fprintf(fp, "\t%0.2f", 0.0);
	if(CmapChimQuality >= 2){
	  if(pmap->SegDupL[c]){
	    if(DEBUG && k <= N) assert(isfinite(pmap->SegDupL[c][newid]) && -1.0 <= pmap->SegDupL[c][newid] && pmap->SegDupL[c][newid] <= 100.001f);
	    fprintf(fp, "\t%0.2f", (k <= N) ? pmap->SegDupL[c][newid] : 0.0);
	  } else
	    fprintf(fp, "\t%0.2f", 0.0);
	  if(pmap->SegDupR[c]){
	    if(DEBUG && k <= N) assert(isfinite(pmap->SegDupR[c][newid]) && -1.0 <= pmap->SegDupR[c][newid] && pmap->SegDupR[c][newid] <= 100.001f);
	    fprintf(fp, "\t%0.2f", (k <= N) ? pmap->SegDupR[c][newid] : 0.0);
	  } else
	    fprintf(fp, "\t%0.2f", 0.0);

	  if(pmap->FragileEndL[c]){
	    if(DEBUG && k <= N) assert(isfinite(pmap->FragileEndL[c][newid]) && 0.0 <= pmap->FragileEndL[c][newid] /* && pmap->FragileEndL[c][newid] <= 1000.001f*/);
	    fprintf(fp, "\t%0.2f", (k <= N) ? pmap->FragileEndL[c][newid] : 0.0);
	  } else
	    fprintf(fp, "\t%0.2f", 0.0);
	  if(pmap->FragileEndR[c]){
	    fprintf(fp, "\t%0.2f", (k <= N) ? pmap->FragileEndR[c][newid] : 0.0);
	    if(DEBUG && k <= N && !(isfinite(pmap->FragileEndR[c][newid]) && 0.0 <= pmap->FragileEndR[c][newid] /* && pmap->FragileEndR[c][newid] <= 1000.001f */)){
	      fflush(fp);
	      printf("pmap->id= %lld: k=%d,N=%d,c=%d,newid=%d,pmap->numsite[c]=%d:pmap->FragileEndR[c][newid]= %0.8e\n",pmap->id, k,N,c,newid,pmap->numsite[c],pmap->FragileEndR[c][newid]);
	      fflush(stdout);
	      assert(isfinite(pmap->FragileEndR[c][newid]) && 0.0 <= pmap->FragileEndR[c][newid] /* && pmap->FragileEndR[c][newid] <= 1000.001f*/);
	    }
	  } else
	    fprintf(fp, "\t%0.2f", 0.0);

	  if(pmap->OutlierFrac[c]){
	    if(DEBUG && k <= N) assert(isfinite(pmap->OutlierFrac[c][newid]) && 0.0 <= pmap->OutlierFrac[c][newid] && pmap->OutlierFrac[c][newid] <= 100.001f);
	    fprintf(fp, "\t%0.2f", (k <= N) ? pmap->OutlierFrac[c][newid] : 0.0);
	  } else
	    fprintf(fp, "\t%0.2f", 0.0);

	  if(pmap->ChimNorm[c]){
	    if(DEBUG && k <= N) assert(isfinite(pmap->ChimNorm[c][newid]) && -1.0 <= pmap->ChimNorm[c][newid]);
	    fprintf(fp, "\t%0.2f", (k <= N) ? pmap->ChimNorm[c][newid] : 0.0);
	  } else
	    fprintf(fp, "\t%0.2f", 0.0);
	}
      }
      if(MaskPresent){
	size_t Mask = 0;
	if(k <= N){
	  if(pmap->Mask[c])
	    Mask = pmap->Mask[c][newid];
	} else {
	  for(int lc = 0; lc < colors; lc++)
	    if(pmap->Mask[lc])
	      Mask |= pmap->Mask[lc][pmap->numsite[lc]+1];
	}
	if(VERB>=2 && Mask){
	  printf("Writing to %s, id=%lld, k=%d/%d: Mask= %lx\n",filename,pmap->id,k,N,Mask);
	  fflush(stdout);
	}
	fprintf(fp, "\t%lx", Mask);
      }

      if(FragSdPresent){
	if(pmap->FragSd[c]){
	  if(DEBUG && k <= N && !(isfinite(pmap->FragSd[c][newid]) && 0.0 <= pmap->FragSd[c][newid])){
	    printf("output_smap:filename=%s:maps[%d]:id=%lld,k=%d,N=%d,c=%d,newid=%d,FragSd[c]= %p, FragSd[c][newid]= %0.6e\n",
		   filename,i,pmap->id,k,N,c,newid,pmap->FragSd[c],pmap->FragSd[c][newid]);
	    printf("\t CmapChimQuality= %d, OutlierQuality= %d\n",CmapChimQuality, OutlierQuality);
	    fflush(stdout);fflush(fp);
	    assert(isfinite(pmap->FragSd[c][newid]) && 0.0 <= pmap->FragSd[c][newid]);
	  }
	  fprintf(fp,"\t%0.1f", (k <= N) ? pmap->FragSd[c][newid] : 0.0);
	} else
	  fprintf(fp,"\t%0.1f", 0.0);

	if(pmap->ExpSd[c]){
	  if(DEBUG && k <= N && !(isfinite(pmap->ExpSd[c][newid]) && 0.0 <= pmap->ExpSd[c][newid])){
	    printf("output_smap:filename=%s:maps[%d]:id=%lld,k=%d,N=%d,c=%d,newid=%d,ExpSd[c]= %p, ExpSd[c][newid]= %0.6e\n",
		   filename,i,pmap->id,k,N,c,newid,pmap->ExpSd[c],pmap->ExpSd[c][newid]);
	    printf("\t CmapChimQuality= %d, OutlierQuality= %d\n",CmapChimQuality, OutlierQuality);
	    fflush(stdout);fflush(fp);
	    assert(isfinite(pmap->ExpSd[c][newid]) && 0.0 <= pmap->ExpSd[c][newid]);
	  }
	  fprintf(fp,"\t%0.1f", (k <= N) ? pmap->ExpSd[c][newid] : 0.0);
	} else
	  fprintf(fp,"\t%0.1f", 0.0);

	if(pmap->FragCov[c]){
	  if(DEBUG && k <= N && !(isfinite(pmap->FragCov[c][newid]) && 0.0 <= pmap->FragCov[c][newid])){
	    printf("output_smap:filename=%s:maps[%d]:id=%lld,k=%d,N=%d:c=%d,newid=%d,FragCov[c]= %p, FragCov[c][newid] = %0.6e\n",
		   filename,i,pmap->id,k,N,c,newid,pmap->FragCov[c],pmap->FragCov[c][newid]);
	    printf("\t CmapChimQuality= %d, OutlierQuality= %d\n",CmapChimQuality, OutlierQuality);
	    fflush(stdout);fflush(fp);
	    assert(isfinite(pmap->FragCov[c][newid]) && 0.0 <= pmap->FragCov[c][newid]);
	  }
	  fprintf(fp,"\t%0.2f", (k <= N) ? pmap->FragCov[c][newid] : 0.0);
	} else
	  fprintf(fp,"\t%0.2f", 0.0);

	if(pmap->FragChiSq[c]){
	  if(DEBUG && k <= N && !(isfinite(pmap->FragChiSq[c][newid]) && (NoBreak ? -2.0 : 0.0) <= pmap->FragChiSq[c][newid])){
	    printf("output_smap:filename=%s:maps[%d]:id=%lld,k=%d,N=%d:c=%d,newid=%d,FragChiSq[c]= %p, FragChiSq[c][newid] = %0.6e\n",
		   filename,i,pmap->id,k,N,c,newid,pmap->FragChiSq[c],pmap->FragChiSq[c][newid]);
	    printf("\t CmapChimQuality= %d, OutlierQuality= %d\n",CmapChimQuality, OutlierQuality);
	    fflush(stdout);fflush(fp);
	    assert(isfinite(pmap->FragChiSq[c][newid]) && (NoBreak ? -2.0 : 0.0) <= pmap->FragChiSq[c][newid]);
	  }
	  fprintf(fp,"\t%0.2e", (k <= N) ? pmap->FragChiSq[c][newid] : 1.0);
	} else
	  fprintf(fp,"\t%0.2e", 1.0);
      }

      if(HapSitePresent){
	Ccontig *contig = pmap->contig;
	if(!(contig && contig->HapSite[0]))
	  fprintf(fp,"\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.2f\t%0.2f", 3, 0.0, -1.0, 0.0, -1.0, 0.0);
	else if(k <= N){
	  int HapSite = contig->HapSite[c] ? contig->HapSite[c][newid] : 3;
	  double HapDelta = contig->HapDelta[c] ? contig->HapDelta[c][newid-1] : 0.0;
	  double HapSiteScore = contig->HapSiteScore[c] ? contig->HapSiteScore[c][newid] * invLog10 : 0.0;
	  double HapDeltaScore = contig->HapDeltaScore[c] ? contig->HapDeltaScore[c][newid-1] * invLog10 : 0.0;
	  double HapSitePhase = (contig->HapSitePhase[c] && 1 <= HapSite && HapSite <= 2) ? contig->HapSitePhase[c][newid] * invLog10 : -1.0;
	  double HapDeltaPhase = contig->HapDeltaPhase[c] ? contig->HapDeltaPhase[c][newid-1] * invLog10 : -1.0;
	  fprintf(fp,"\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.2f\t%0.2f",HapSite,HapDelta*1000.0,HapSitePhase,HapDeltaPhase,HapSiteScore,HapDeltaScore);
	} else if(k >= N+1){/* right end */
	  double HapDelta = contig->HapDelta[0] ? contig->HapDelta[0][N] : 0.0;
	  double HapDeltaScore = contig->HapDeltaScore[0] ? contig->HapDeltaScore[0][N] * invLog10 : 0.0;
	  double HapDeltaPhase = contig->HapDeltaPhase[0] ? contig->HapDeltaPhase[0][N] * invLog10 : -1.0;
	  if(DEBUG && !(fabs(HapDelta) <= pmap->site[0][N+1] - pmap->site[0][N])){
	    printf("id=%lld:k=%d,N=%d:site[N]= %0.6f, site[N+1]= %0.6f, HapDelta[N]= %0.6f\n", pmap->id,k,N,pmap->site[0][N],pmap->site[0][N+1],contig->HapDelta[0][N]);
	    fflush(stdout);
	    assert(fabs(HapDelta) <= pmap->site[0][N+1] - pmap->site[0][N]);
	  }
	  fprintf(fp,"\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.2f\t%0.2f",0,HapDelta*1000.0, -1.0, HapDeltaPhase, -1.0, HapDeltaScore);
	} 
      }
      if(Tracksites){
	if(DEBUG && k <= N && !(pmap->truesiteL[c] && pmap->truesiteR[c])){
	  printf("output_cmap:id=%lld:k=%d,N=%d,c=%d,newid=%d,SiteID=%d:pmap->truesiteL[c]=%p,pmap->truesiteR[c]=%p\n",
		 pmap->id,k,N,c,newid,SiteID,pmap->truesiteL[c],pmap->truesiteR[c]);
	  fflush(stdout);
	  assert(pmap->truesiteL[c] && pmap->truesiteR[c]);
	}
	fprintf(fp, "\t%lld\t%lld", (k <= N) ? pmap->truesiteL[c][newid] : -1, (k <= N) ? pmap->truesiteR[c][newid] : -1);
      }

      if(DEBUG && k > 1 && !((pn[k-2].color != c) ? pn[k-1].site >= pn[k-2].site : pn[k-1].site > pn[k-2].site)){
	fflush(fp);
	printf("Out of order sites in .cmap output:pmap->id=%lld: pn[%d]:Site=%d(color=%d,newid=%d),pos=%0.1f vs pn[%d]:Site=%d(color=%d,newid=%d),pos=%0.1f,N=%d\n",
	       pmap->id,k-2,SiteID-1,pn[k-2].color,pn[k-2].newid,pn[k-2].site*1000.0,k-1,SiteID,pn[k-1].color,pn[k-1].newid,pn[k-1].site*1000.0,N);
	fflush(stdout);
	assert(pn[k-1].site > pn[k-2].site);
      }
      if(DEBUG && k <= N) assert(c >= 0 && c < colors);
      if(VERB>=2 && !strcmp(filename,"test_q.cmap")){
	printf("pmap->id=%lld:k=%d,N=%d,c=%d:CmapSNR=%d,SNRpresent=%d,pmap->SNRcnt[c]=%p\n", pmap->id, k,N,c,CmapSNR,SNRpresent, pmap->SNRcnt[c]);
	fflush(stdout);
      }
      if(SNRpresent){
	if(k <= N && pmap->SNRcnt[c]){
	  if(DEBUG) assert(isfinite(pmap->lnSNRsd[c][newid]));
	  fprintf(fp,"\t%0.4f\t%0.4f\t%d", pmap->SNRgmean[c][newid], pmap->lnSNRsd[c][newid], pmap->SNRcnt[c][newid]);
	  for(int u = 0; u < pmap->SNRcnt[c][newid]; u++)
	    fprintf(fp,"\t%0.4f",pmap->SNRdist[c][newid][u]);
	} else 
	  fprintf(fp,"\t%0.4f\t%0.4f\t%d", 0.0, 0.0, 0);
      } else if(SNRval){
	if(k <= N){
	  if(DEBUG) assert(pmap->SNR[c] != NULL);
	  fprintf(fp,"\t%0.4f\t%0.4f", pmap->SNR[c][newid], 0.0);
	} else
	  fprintf(fp,"\t%0.4f\t%0.4f", 0.0, 0.0);
      }
      
      fprintf(fp,"\n");
    }
    delete [] pn;
  }
  (void)FILEclose(fp);
  free(write_buffer);

  if(1){
    errno = 0;
    if(rename(tmpfilename,origfilename)){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to rename %s to %s:errno=%d:%s\n",tmpfilename,origfilename,eno,err);
      fflush(stdout);  exit(1);
    }
  }
  if(VERB){
    printf("Completed output of %s: cum wall time= %0.6f\n",origfilename,wtime());
    fflush(stdout);
  }

  free(origfilename);
  free(tmpfilename);
}

extern char *Nickase[MAXCOLOR];

extern int RunIndexWarning, RunDataWarning;/* see input_vfixx.cpp */
static int BNXversionWARNING = 0;

/* output maps in .bnx format : If pFP != 0 : continue output to open File (if pFP[g] != NULL && != -1) OR open file (append if pFP[g]== -1) and allocate pwrite_buffer[g] and save FILE pointer to pFP[g] */
/* The caller is responsible for mutex locking the file pointer (if needed) */
/* If raw != 0 : output the original map (without -bpp or -resbias correction) : currently only used with -ScanCorrection to output *_rescaled.bnx */
/* If tid >= 0 : the thread that called output_bnx() from multithreaded section */
/* If mapWT != 0 : mapWT[maps[i]->mapid] is a WT <= 1.0 added to the header of each molecule (otherwise a value of 1.0 is implied) */
/* If smapWT != 0 : smapWT[maps[i]->mapid][g] is a WT <= 1.0 added to the header of each molecule (otherwise a value of 1.0 is implied) */
#if SPARSE_MAPWT == 1
void output_bnx(char *basename, Cmap *maps[], int start, int end, int raw, FILE **pFP, char **pwrite_buffer, int tid, int g, omp_lock_t *lock, int gmax, int refid, float *mapWT, std::map<int,float> *smapWT)
#else
void output_bnx(char *basename, Cmap *maps[], int start, int end, int raw, FILE **pFP, char **pwrite_buffer, int tid, int g, omp_lock_t *lock, int gmax, int refid, float *mapWT)
#endif
{
  if(strstr(basename,"/dev/null"))
    return;

  FILE *fp = NULL;
  char *write_buffer = NULL;

  if(pFP){
    if(DEBUG) assert(pwrite_buffer != NULL && 1 <= g && g <= gmax && lock != NULL && tid >= 0);
    #pragma omp critical(pFPlock) /* critical section required to insure read is not just in this thread's register */
    {
      fp = pFP[g];
      write_buffer = pwrite_buffer[g];
    }
  }

  if(VERB>=2){
    printf("output_xmap(%s):BNXVersion=%d\n",basename,BNXVersion);
    fflush(stdout);
  }

  if(TSAN || BNXVersion < 0){
    #pragma omp critical(BNXVersion)
    BNXVersion = max(0,BNXVersion);
  }

  /* count number of maps exluding map splits */
  int mapcnt = 0;
  for(int i = start; i <= end; i++)
    if(!maps[i]->origmap)
      mapcnt++;

  if(VERB>=2){
    printf("output_bnx:    start=%d,end=%d,mapcnt=%d\n",start,end,mapcnt);
    fflush(stdout);
  }

  int HasStartLoc = 0, HasTrueSite = MapTrueSite;

  FILE *origfp = fp;

  if(fp==NULL || fp== (FILE *)-1){
    if(DEBUG && fp == (FILE *) -1) assert(pFP && pFP[g] == (FILE *)-1);

    char filename[PATH_MAX];
    strcpy(filename,basename);
    /* no need to strip suffix : for default basename this was already done in RefAligner.cpp */
    if(1){
      int i = strlen(filename);
      sprintf(&filename[i],".bnx");
    }

    if(VERB){
      #pragma omp critical
      {
	if(fp==NULL){
	  if(tid >= 0)
	    printf("Generating %s (starting with %d maps), group=%d, refid=%d: tid=%d",filename, mapcnt, g, refid, tid);
	  else
	    printf("Generating %s (with %d maps), refid=%d",filename, mapcnt, refid);
	  if(VERB>=2){
	    printf(":BNXVersion=%d",BNXVersion);
	    if(BNXVersion >= 1)
	      printf(":colors=%d,MapSNR=%d,MapIntensity=%d,MapStitched=%d,MapStitchLoc=%d,MapPSFWidth=%d",colors,MapSNR,MapIntensity,MapStitched,MapStitchLoc,MapPSFWidth);
	  }
	  printf("\n");
	} else if(VERB>=2)
	  printf("Re-Opening %s (%d maps), group=%d, refid=%d: tid=%d\n",filename, mapcnt, g, refid, tid);
	fflush(stdout);
      }
    }
    
    if(fp==NULL && checkFile(filename))
      return;

    const char *mode = (fp==NULL) ? "w" : "a";

    if(fp==NULL){/* make sure parent folder exists : this is only needed for stage0 -mapped BNX output where parent folder is a sub-directory of the main output folder */
      char *r = strrchr(filename,'/');
      if(r != NULL){
	char dirname[PATH_MAX];
	strcpy(dirname,filename);
	dirname[r - filename] = 0;
	
	errno = 0;
	if(mkdir(dirname,S_IRWXU|S_IRGRP|S_IROTH) < 0 && errno != EEXIST){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("Failed to create directory %s for output file %s:%s\n",dirname,filename,err);
	  fflush(stdout);exit(1);
	}
      }
    }

    for(int cnt = 0; (fp = fopen(filename,mode)) == NULL; cnt++){
      int eno = errno;
      char *err = strerror(errno);
      int gcnt = 0;

      #pragma omp critical(filelock)
      {
	if(tid >= 0)
	  printf("tid=%d:",tid);
	if(mode[0] == 'a')
	  printf("failed to re-open");
	else
	  printf("failed to open");
	printf(" file %s for writing merged .bnx file of %d maps:errno=%d:%s",filename, mapcnt,eno,err);
	if(eno != 13){// errno == 13 : permission denied
          if(cnt <= 100)
	    printf(" (retry %d of 100 after 1 second)",cnt+1);
	  else if(eno != 28)// errno == 28 means no space left on device
	    printf(" (retry %d after 60 seconds)",cnt+1);
        } else if(cnt <= 0){
          printf(" (retry once after 1 second)");
        }
	printf("\n");

#ifndef WIN32
        char buf[BUFSIZ];

	if(!USE_MIC && cnt > 0 && eno == 24){ /* errno==24 means too many open file handles : output list of open file handles for current process */
	  pid_t pid = getpid();
	  printf("List of open files of current pid=%lld\n",(long long)pid);
	  fflush(stdout);

	  sprintf(buf,"ls -l /proc/%lld/fd", (long long)pid);
	  int err = system(buf);
	  if(err)
	    printf("system(%s) failed: return code = %d\n", buf, err);
	}
	if(!USE_MIC && cnt >= 100  && eno==24){  /* output list of all open file handles */
	  printf("List of all open files (lsof)\n");
	  fflush(stdout);
	  
	  pid_t pid = getpid();
#if USE_MIC
	    sprintf(buf,"ssh -t -t -oStrictHostKeyChecking=no -oBatchMode=yes host /usr/sbin/lsof -p %lld",(long long)pid);
#else
	    sprintf(buf,"/usr/sbin/lsof -p %lld", (long long)pid);
#endif
	  int err = system(buf);
	  if(err)
	    printf("system(%s) failed: return code = %d\n", buf, err);
	}
#endif
	fflush(stdout);

#ifdef _OPENMP
	if(eno==24) {
          /* close some of the open files that are currently not locked by other threads */
          // NOTE : this does NOT seem to work on some MICs : subsequent newly opened files all fail to append data
	  //             The problem is caused by using outdated scif_io.c that did not correctly handle the "a" mode flag
          for(int group = 1; group <= gmax; group++){
            if(group == g)
	      continue;
	    if(omp_test_lock(&lock[group])){
              #pragma omp critical(pFPlock) /* critical section required to insure read/update is not just in this thread's register */
	      {
                if(pFP[group] && pFP[group] != (FILE *)-1){
	          if(VERB>=2){
                    printf("tid=%d:calling FILEclose(FP[%d]= %p -> -1)\n",tid,group,pFP[group]);
		    fflush(stdout);
                  }

		  (void)FILEclose(pFP[group], 1);
		  pFP[group] = (FILE *)-1;
		    
		  if(DEBUG) assert(pwrite_buffer[group] != NULL);
		  free(pwrite_buffer[group]);// WAS128 delete [] pwrite_buffer[group]
		  pwrite_buffer[group] = NULL;
		  
		  gcnt++;
                }
	      }
	      omp_unset_lock(&lock[group]);
	    }
	  }
        }// if(eno==24)
#endif
      }
      if(!gcnt){
        if(eno != 13){// errno == 13 : permission denied
          if(cnt <= 100){
	    sleep(1);
	    continue;
          }
	  if(eno != 28) {// errno == 28 means no space left on device
	    sleep(60);
	    continue;
          }
        }
	if(cnt <= 0){// typically errno == 13 : permission denied
	  sleep(1);
	  continue;
	}
	printf("ERROR:Too many failures with errno == %d\n",eno);
	fflush(stdout); exit(1);
      } else {
	printf("Closed %d open files: retrying after 1 second\n",gcnt);
	fflush(stdout);

	sleep(1);
      }
    }

    //    size_t blen = USE_MIC ? 128*1024 : 10*1024*1024;/* Cannot exceed scif_io buffer size = 10M */
    size_t blen = 10*1024*1024;/* matches scif_io buffer size */
    if(DEBUG>=2) assert(write_buffer == NULL);
    if(pFP && blen * gmax > (unsigned long long)(0.5 * MaxMem * 1024.0 * 1024.0 * 1024.0)){/* adjust buffer size so no more than half of MaxMem is used for output buffers */
      blen = (unsigned long long)(0.5 * MaxMem * 1024.0 * 1024.0 * 1024.0) / gmax;
      blen = max(1ull, blen/4096ull);
      blen *= 4096ull;
    }

    write_buffer = (char *)malloc(blen);
    if(!write_buffer){
      printf("output_bnx:malloc(%llu) failed\n",(unsigned long long)blen);
      fflush(stdout);exit(1);
    }
    if(VERB>=3){
      #pragma omp critical
      {
        printf("tid=%d:group=%d, fp=%p,buffer=%p,blen=%lu:calling setbuffer\n",tid,g,fp,write_buffer,blen);
        fflush(stdout);
      }
    }
    setbuffer(fp, write_buffer, blen);

    if(VERB>=2){
      #pragma omp critical
      {
        printf("tid=%d:group=%d, fp -> %p, buffer -> %p, filename= %s\n",tid,g,fp,write_buffer,filename);
        fflush(stdout);
      }
    }

    if(DEBUG/* HERE >=2 */ && origfp == (FILE *)(-1)){
      errno = 0;
      long origpos_cur = ftell(fp);
      if(origpos_cur < 0 || errno){
	int eno = errno;
	char *err = strerror(eno);
	#pragma omp critical
	{
	  printf("ERROR:tid=%d:group=%d,fp=%p,buffer=%p: orig ftell(fp)=%ld,eno=%d:%s\n",tid,g,fp,write_buffer,origpos_cur,eno,err);
	  fflush(stdout);exit(1);
	}
      }

      errno = 0;
      int ret = fprintf(fp,"# file reopened\n");
      if(ret < 0){
	int eno = errno;
	char *err = strerror(eno);
	#pragma omp critical
	{
	  printf("ERROR:tid=%d:group=%d,fp=%p,buffer=%p: fprintf()=%d,eno=%d:%s\n",tid,g,fp,write_buffer,ret,eno,err);
	  fflush(stdout);exit(1);
	}
      }

      errno = 0;
      long pos_cur = ftell(fp);
      if(pos_cur < 0 || errno){
	int eno = errno;
	char *err = strerror(eno);
	#pragma omp critical
	{
	  printf("ERROR:tid=%d:group=%d,fp=%p,buffer=%p: ftell(fp)=%ld,eno=%d:%s\n",tid,g,fp,write_buffer,origpos_cur,eno,err);
	  fflush(stdout);exit(1);
	}
      }

      if(VERB>=2){
	printf("tid=%d:group=%d,fp=%p,buffer=%p,filename=%s:wrote %d bytes:pos_cur=%ld->%ld\n",tid,g,fp,write_buffer,filename,ret,origpos_cur,pos_cur);
	fflush(stdout);
      }
      FFlush(fp, 1);
    }
  }

  if(origfp == NULL){/* write to new file */

    /* write out commandline */
    printversion(fp);

    if(BNXVersion == 1){/* check that all maps being output have Version 1.0 information (may be missing for .cmap, or .spot file inputs) */
      int maxScanNumber = 1;
      for(int i = start; i <= end; i++) {
	if(maps[i]->origmap)
	  continue;
	if(maps[i]->OriginalMoleculeId < 0){
	  if(!BNXversionWARNING){
	    if(VERB){
	      printf("maps[%d]:id=%lld,OriginalMoleculeId=%d (should be >= 0)\n",i,maps[i]->id,maps[i]->OriginalMoleculeId);
	      printf("   BNX input will be treated as 0.1 version and 1.x output will be generated using default values for all post BNX 0.1 header fields\n");
	      fflush(stdout);
	    }
	    BNXversionWARNING = 1;
	  }
	  maps[i]->OriginalMoleculeId = maps[i]->mapid + 1;
	  maps[i]->ScanNumber = 1;
	  maps[i]->ScanDirection = -1;
	  maps[i]->ChipId = ChipIdAlloc("unknown");
	  maps[i]->FlowCell = 1;
	  maps[i]->AvgIntensity = 0.0;
	  maps[i]->bbSNR = 0.0;
	  if(MapIntensity || MapStitchLoc || MapPSFWidth){
	    if(VERB){
	      printf("WARNING: Intensity, stitchLocation or PSFWidth information cannot be output from BNX 0.1 or CMAP input\n");
	      fflush(stdout);
	    }
	    MapIntensity = MapStitchLoc = MapPSFWidth = 0;
	    for(int c = 0; c < colors; c++)
	      maps[i]->NumStitch[c] = 0;
	  }
	  maps[i]->RunIndex = 0;
	  maps[i]->UniqueScanId = 0;
	  maps[i]->RunId = 0;
	  maps[i]->ScanId = 0;
	} else if(maps[i]->ScanNumber > maxScanNumber){
	  maxScanNumber = maps[i]->ScanNumber;
	  if(BNXversionWARNING){
	    printf("WARNING:maps[%d]:id=%lld,ScanNumber=%d: This is inconsistent with previous molecule with OriginalMoleculeId= -1\n",
		   i,maps[i]->id,maps[i]->ScanNumber);
	    fflush(stdout);
	  }
	}
      }
      if(BNXversionWARNING && maxScanNumber > 1){
	printf("ERROR: max value of maps[i]->ScanNumber = %d, which means some input is BNX 1.x, but some molecules have no OriginalMoleculeID suggesting BNX 0.1 or CMAP input\n",
	       maxScanNumber);
	fflush(stdout);exit(1);
      }
      if(BNXversionWARNING){/* create a RunData string */
	if(RunDataListLen > 0){
	  printf("ERROR: While converting CMAP to BNX 1.x found previous RunData information of unknown origin:\n");
	  for(int i = 0; i < RunDataListLen; i++)
	    printf("%s\n",RunDataList[i]);
	  fflush(stdout);exit(1);
	}
	if(num_files != 1){
	  printf("ERROR: Cannot convert multiple CMAP input files to BNX 1.x : convert one file at a time\n");
	  fflush(stdout);exit(1);
	}
	if(!RunDataListMax){
	  RunDataListMax = 16;
	  delete [] RunDataList;
	  RunDataList = new char*[RunDataListMax];
	}
	char RunData[1024];
	sprintf(RunData,"# Run Data\t%s\tUnknown\t01/01/2000 01:00:00 AM\t0\t1.0\t500\t1\tUnknown,1\t1",vfixx_filename[0]);
	RunDataList[RunDataListLen++] = strdup(RunData);
      }
    }

    /* output key header lines */
    if(BNXVersion <= 0)
      fprintf(fp,"# BNX File Version:\t0.1\n");
    else if(RunIndexWarning || RunDataWarning || RunDataListLen <= 0 || BNXnoRunData)
      fprintf(fp,"# BNX File Version:\t1.0\n");
    else if(BNXMinorVersion <= 2)
      fprintf(fp,"# BNX File Version:\t1.2\n");
    else
      fprintf(fp,"# BNX File Version:\t1.3\n");
    fprintf(fp,"# Label Channels:\t%d\n", colors);
    if(end >= start && maps[start]->Nickase[0]){
      if(VERB>=2){
        printf("output_bnx:basename=%s,start=%d,end=%d:global Nickase[0]=%p(%s),maps[start]:id=%lld,Nickase=%p,Nickase[0]=%p(%s)\n",
	       basename,start,end,Nickase[0],Nickase[0],maps[start]->id,maps[start]->Nickase,maps[start]->Nickase[0],maps[start]->Nickase[0]);
        fflush(stdout);
      }
      for(register int c = 0; c < colors; c++)
        if(maps[start]->Nickase[c]){
	  int L = strlen(maps[start]->Nickase[c]);
	  while(L > 0 && maps[start]->Nickase[c][L-1] == '\n'){
	    //	  printf("maps[start]->Nickase[%d] has newline at end: removed\n",c);
	    maps[start]->Nickase[c][--L] = '\0';
	  }
	  fprintf(fp,"# Nickase Recognition Site %d:\t%s\n", c+1, maps[start]->Nickase[c]);
        }
    } else
      for(register int c = 0; c < colors; c++)
        fprintf(fp,"# Nickase Recognition Site %d:\t%s\n", c+1, "unknown");

    fprintf(fp,"# Bases per Pixel:\t500\n");

    if(BNXVersion <= 0){
      fprintf(fp,"# Number of Molecules:\t%d\n", mapcnt);
      for(int i = start; i <= end; i++){
        Cmap *pmap = maps[i];
	if(pmap->origmap)
	  continue;
	if(pmap->startloc > 0.0 && pmap->endloc > 0.0){
	  HasStartLoc = 1;
	  break;
	}
      }
      for(int i = start; i <= end; i++){
        Cmap *pmap = maps[i];
	if(pmap->origmap)
	  continue;
	if(pmap->truesiteL[0] != NULL)
	  HasTrueSite = 1;
	break;
      }
      
      fprintf(fp,"#0h LabelChannel\tMapID\tLength");
#if SPARSE_MAPWT == 1
      if(mapWT || smapWT)
#else
      if(mapWT)
#endif
	fprintf(fp,"\tMapWT");
      if(HasStartLoc)
	fprintf(fp,"\tStartLoc\tEndLoc\tStartLoc2\tEndLoc2\tFlipped\tChimFlip");
      fprintf(fp,"\n");

      fprintf(fp,"#0f int\tint\tfloat");
#if SPARSE_MAPWT == 1
      if(mapWT || smapWT)
#else
      if(mapWT)
#endif
	fprintf(fp,"\tfloat");
      if(HasStartLoc)
	fprintf(fp,"\tfloat\tfloat\tfloat\tfloat\tint\tint");
      fprintf(fp,"\n");

      fprintf(fp,"#1h LabelChannel\tLabelPositions[N]\n");
      fprintf(fp,"#1f int\tfloat\n");
      if(colors > 1){
	fprintf(fp,"#2h LabelChannel\tLabelPositions[N]\n");
	fprintf(fp,"#2f int\tfloat\n");
      }
    } else {/* BNXVersion >= 1 */
      if(!BNXnoRunData){
	if(BNXRunHeaderChim > 0){
	  if(colors==2)
	    fprintf(fp,"#rh SourceFolder\tInstrumentSerial\tTime\tNanoChannelPixelsPerScan\tStretchFactor\tBasesPerPixel\tNumberofScans\tChipId\tFlowCell\tSNRFilterType\tMinMoleculeLength\tMinLabelSNR1\tMinLabelSNR2\tFeatureMinSlope\tFeatureMaxSlope\tRunId\n");
	  else
	    fprintf(fp,"#rh SourceFolder\tInstrumentSerial\tTime\tNanoChannelPixelsPerScan\tStretchFactor\tBasesPerPixel\tNumberofScans\tChipId\tFlowCell\tSNRFilterType\tMinMoleculeLength\tMinLabelSNR\tFeatureMinSlope\tFeatureMaxSlope\tRunId\n");
	} else {
	  if(colors==2)
	    fprintf(fp,"#rh SourceFolder\tInstrumentSerial\tTime\tNanoChannelPixelsPerScan\tStretchFactor\tBasesPerPixel\tNumberofScans\tChipId\tFlowCell\tSNRFilterType\tMinMoleculeLength\tMinLabelSNR1\tMinLabelSNR2\tRunId\n");
	  else
	    fprintf(fp,"#rh SourceFolder\tInstrumentSerial\tTime\tNanoChannelPixelsPerScan\tStretchFactor\tBasesPerPixel\tNumberofScans\tChipId\tFlowCell\tSNRFilterType\tMinMoleculeLength\tMinLabelSNR\tRunId\n");
	}
	for(int i = 0; i < RunDataListLen; i++)
	  fprintf(fp,"%s\n",RunDataList[i]);
      }
      if(tid >= 0)
	fprintf(fp,"# Number of Molecules for refid= %d:\t%d\n", refid, mapcnt);
      else
	fprintf(fp,"# Number of Molecules:\t%d\n", mapcnt);
      if(BNXminlen > 0.0)
	fprintf(fp,"# Min Molecule Length (Kb):\t%0.3f\n",BNXminlen);
      if(BNXminSNR[0] > 0.0 || (colors >= 2 && BNXminSNR[1] > 0.0)){
	fprintf(fp,"# Min Label SNR:");
	for(int c = 0; c < colors; c++)
	  fprintf(fp,"\t%0.2f",BNXminSNR[c]);
	fprintf(fp,"\n");
      }
      if(VERB>=2){
        printf("mapWT= %p, BNXheaderFXY=%d\n",mapWT,BNXheaderFXY);
	fflush(stdout);
      }

      fprintf(fp,"#0h LabelChannel\tMoleculeID\tLength\tAvgIntensity\tSNR\tNumberofLabels\tOriginalMoleculeId\tScanNumber\tScanDirection\tChipId\tFlowcell\tRunId");
      if(BNXheaderFXY > 0)
	fprintf(fp,"\tColumn\tStartFOV\tStartX\tStartY\tEndFOV\tEndX\tEndY");
      if(BNXheaderChim > 0)
	fprintf(fp,"\tFeatureMinPos\tFeatureMinScore\tFeatureMaxPos\tFeatureMaxScore\tFeatureMinScaled\tFeatureMaxScaled");
#if SPARSE_MAPWT == 1
      if(mapWT || smapWT)
#else
      if(mapWT)
#endif
	fprintf(fp,"\tMapWT");
      if(!(RunIndexWarning || RunDataWarning || RunDataListLen<=0 || BNXnoRunData))
	fprintf(fp,"\tGlobalScanNumber");
      fprintf(fp,"\n");

      fprintf(fp,"#0f int\t int\t float\tfloat\tfloat\tint\tint\tint\tint\tstring\tint\tint");
      if(BNXheaderFXY > 0)
	fprintf(fp,"\tint\tint\tint\tint\tint\tint\tint");
      if(BNXheaderChim > 0)
	fprintf(fp,"\tint\tfloat\tint\tfloat\tfloat\tfloat");
#if SPARSE_MAPWT == 1
      if(mapWT || smapWT)
#else
      if(mapWT)
#endif
	fprintf(fp,"\tfloat");
      if(!(RunIndexWarning || RunDataWarning || RunDataListLen<=0 || BNXnoRunData))
	fprintf(fp,"\tint");
      fprintf(fp,"\n");

      fprintf(fp,"#1h LabelChannel\tLabelPositions[N]\n");
      fprintf(fp,"#1f int\tfloat\n");
      if(colors > 1){
	fprintf(fp,"#2h LabelChannel\tLabelPositions[N]\n");
	fprintf(fp,"#2f int\tfloat\n");
      }
      fprintf(fp,"#Qh QualityScoreID\tQualityScores[N]\n");
      fprintf(fp,"#Qf string\tfloat[N]\n");

      if(HasTrueSite)
	for(int c = 0; c < colors; c++)
	  fprintf(fp,"# Quality Score %s: True Sites (intervals, two values) for channel %d\n",VQ_TRUESITE[c][BNXVersion], c+1);
      if(MapSNR)
	for(int c = 0; c < colors; c++)
	  fprintf(fp,"# Quality Score %s: Label SNR for channel %d\n",VQ_SNR[c][BNXVersion], c+1);
      if(MapIntensity)
	for(int c = 0; c < colors; c++)
	  fprintf(fp,"# Quality Score %s: Label Intensity for channel %d\n",VQ_INTENSITY[c][BNXVersion], c+1);
      if(MapStitched)
	for(int c = 0; c < colors; c++)
	  fprintf(fp,"# Quality Score %s: Stitched Intervals for channel %d\n",VQ_STITCH[c][BNXVersion], c+1);
      if(MapStitchLoc)
	for(int c = 0; c < colors; c++)
	  fprintf(fp,"# Quality Score %s: Stitch Locations for channel %d\n",VQ_STITCHLOC[c][BNXVersion], c+1);
      if(MapPSFWidth)
	for(int c = 0; c < colors; c++)
	  fprintf(fp,"# Quality Score %s: PSF Width for channel %d\n", VQ_PSFWIDTH[c][BNXVersion], c+1);
      if(MapImageLoc)
	for(int c = 0; c < colors; c++)
	  fprintf(fp,"# Quality Score %s: Image Location (FOV,X,Y) for channel %d\n", VQ_IMAGELOC[c][BNXVersion], c+1);
    }

    if(QueryOrigin)
      fprintf(fp,"# Original QueryMap: %s\n", QueryOrigin);

  } else {/* appending to an open FILE pointer */
    if(DEBUG) assert(pwrite_buffer != NULL && 1 <= g && g <= gmax && lock != NULL && tid >= 0);

    if(VERB>=2){
      char filename[PATH_MAX];
      strcpy(filename,basename);
      /* no need to strip suffix : for default basename this was already done in RefAligner.cpp */
      if(1){
	int i = strlen(filename);
	sprintf(&filename[i],".bnx");
      }
      
      #pragma omp critical
      {
        printf("Appending to %s (%d maps), group=%d, refid=%d: tid=%d\n",filename, end - start + 1, g, refid, tid);
	fflush(stdout);
      }
    }
    if(BNXVersion == 1){/* check that all maps being output have Version 1.0 information (may be missing for .cmap, or .spot file inputs) */
      for(int i = start; i <= end; i++) {
	if(maps[i]->origmap)
	  continue;
	if(maps[i]->OriginalMoleculeId < 0){
	  printf("BNX output version %d not possible due to input map=%lld missing Version 1.x information\n",
		 BNXVersion, maps[i]->id);
	  fflush(stdout);  exit(1);
        }
      }
    }

    fprintf(fp,"# Number of Additional Molecules for refid= %d:\t%d\n", refid, mapcnt);
  }

  if(pFP && fp != pFP[g]){
    #pragma omp critical(pFPlock) /* critical section required to insure read/update is not just in this thread's register */
    {
      pFP[g] = fp;
      pwrite_buffer[g] = write_buffer;
      if(VERB>=2){
        printf("tid=%d:Writing %d maps to %s.bnx: group=%d, fp=%p, FP[group]=%p, buffer[group]=%p\n",tid,end-start+1, basename, g, fp, pFP[g], pwrite_buffer[g]);
	fflush(stdout);
      }
    }
  }

  if(tid >= 0)
    fprintf(fp,"# Thread ID=%d\n",tid);

  int linecnt = 0;
  int lines_per_map = 1 + colors*(1 + HasTrueSite + MapSNR + MapStitched + MapIntensity + MapStitchLoc + MapPSFWidth + MapImageLoc);

  /* HERE : if !pFP, try using multithreaded output to per-thread buffers, then sequentially write out the buffers as binary output */

  for(int i = start; i <= end; i++){
    Cmap *pmap = maps[i];
    if(pmap->origmap)
      continue;/* don't output map fragments, since they have the same ID as the complete map */

    int *numsite = pmap->numsite;
    if(DEBUG && !NoStat && maxresbias > mres * 0.5) assert(pmap->rawsite != NULL);
    double **site = (raw && pmap->rawsite[0] != NULL) ? pmap->rawsite : pmap->site;
    double scale = raw ? origPixelLen * 1000.0 / PixelLen : 1000.0;
    if(raw && UniqueScans > 1 && fabs(origPixelLen - 0.500) > 1e-6){
      printf("WARNING: BNX 1.x input has origPixelLen = %0.4f (expected 0.5)\n", origPixelLen);
      fflush(stdout);
    }

    if(BNXVersion <= 0){
      fprintf(fp,"0\t%lld\t%0.2f", pmap->id, site[0][numsite[0]+1]*scale);

      if(DEBUG>=2) assert(0 <= pmap->mapid && pmap->mapid < nummaps);
      if(mapWT){
	if(DEBUG>=2) assert(0.0 <= mapWT[pmap->mapid] && mapWT[pmap->mapid] <= 1.0 + 1e-7);
	fprintf(fp,"\t%0.6e", min(1.0,(double)mapWT[pmap->mapid])/* NEW157*/);
      }
#if SPARSE_MAPWT == 1
      if(smapWT){
	// NOTE : smapWT[pmap->mapid][g] is not thread-safe : if index g is not present in map, an element of value 0 is inserted (and returned) and the insertion is not thread-safe
	std::map<int,float> *p = &smapWT[pmap->mapid];
	double val = 0.0;
	if(p->find(g) != p->end()) 
	  val = p[0][g];
	if(DEBUG>=2) assert(0.0 <= val && val <= 1.0 + 1e-7);
	fprintf(fp,"\t%0.6e", min(1.0, val));
      }
#endif

      if(pmap->startloc > 0.0 && pmap->endloc > 0.0)
	fprintf(fp,"\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%d\t%d", pmap->startloc*1000.0, pmap->endloc*1000.0, pmap->startloc2*1000.0,pmap->endloc2*1000.0, pmap->flipped,pmap->chimflip);
      fprintf(fp,"\n");
    } else {// BNXVersion > 0
      int NumberofLabels = 0;
      for(int c = 0; c < colors; c++)
	NumberofLabels += numsite[c];
      
      fprintf(fp,"0\t%lld\t%0.2f\t%0.2f\t%0.2f\t%d\t%d\t%d\t%d\t%s\t%d",
          pmap->id, site[0][numsite[0]+1]*scale, pmap->AvgIntensity, pmap->bbSNR,  NumberofLabels, pmap->OriginalMoleculeId, pmap->ScanNumber, pmap->ScanDirection, pmap->ChipId, pmap->FlowCell); 
      if(RunIndexWarning || RunDataWarning || RunDataListLen <= 0 || BNXnoRunData)/* input was missing RunIndex or Run Data information */
	fprintf(fp,"\t%d", 1);
      else
	fprintf(fp,"\t%d", pmap->RunIndex+1);

      if(BNXheaderFXY > 0)// add Column etc 
	fprintf(fp,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d",  pmap->Column, pmap->StartFOV, pmap->StartX, pmap->StartY, pmap->EndFOV, pmap->EndX, pmap->EndY); 
      if(BNXheaderChim > 0)// add FeatureMinPos .. FeatureMaxScaled
	fprintf(fp,"\t%d\t%0.6f\t%d\t%0.6f\t%0.6f\t%0.6f", pmap->FeatureMinPos,pmap->FeatureMinScore,pmap->FeatureMaxPos,pmap->FeatureMaxScore,pmap->FeatureMinScaled,pmap->FeatureMaxScaled);

      if(DEBUG>=2) assert(0 <= pmap->mapid && pmap->mapid < nummaps);  
      if(mapWT) { // add mapWT value 
	if(DEBUG>=2) assert(0.0 <= mapWT[pmap->mapid] && mapWT[pmap->mapid] <= 1.0 + 1e-7);
	fprintf(fp,"\t%0.6e", min(1.0,(double)mapWT[pmap->mapid])); 
      }
#if SPARSE_MAPWT == 1
      if(smapWT){
	// NOTE : smapWT[pmap->mapid][g] is not thread-safe : if index g is not present in map, an element of value 0 is inserted (and returned) and the insertion is not thread-safe
	std::map<int,float> *p = &smapWT[pmap->mapid];
	double val = 0.0;
	if(p->find(g) != p->end()) 
	  val = p[0][g];
	if(DEBUG>=2) assert(0.0 <= val && val <= 1.0 + 1e-7);
	fprintf(fp,"\t%0.6e", min(1.0,val));
      }
#endif
      if(RunIndexWarning || RunDataWarning || RunDataListLen <= 0 || BNXnoRunData)
	fprintf(fp,"\t%d\n", 1);
      else
	fprintf(fp,"\t%d\n", pmap->UniqueScanId + 1);
    }
    linecnt++;
    if(DEBUG/* HERE >=2 */) assert((linecnt % lines_per_map) == 1);

    for(int c = 0; c < colors; c++){
      fprintf(fp,"%d",c+1);
      for(int k = 1; k <= numsite[c]+1; k++)
	fprintf(fp,"\t%0.2f", site[c][k]*scale);
      fprintf(fp,"\n");
      linecnt++;
    }
    if(DEBUG/* HERE >=2 */ && !((linecnt % lines_per_map) == ((1+colors) % lines_per_map))){
      #pragma omp critical
      {
	printf("While writing map lines to %s.bnx:\n",basename);
	printf("BNXVersion=%d",BNXVersion);
	if(BNXVersion >= 1)
	  printf(":colors=%d,HasTrueSite=%d,MapSNR=%d,MapIntensity=%d,MapStitched=%d,MapStitchLoc=%d,MapPSFWidth=%d,MapImageLoc=%d\n",
		 colors,HasTrueSite,MapSNR,MapIntensity,MapStitched,MapStitchLoc,MapPSFWidth,MapImageLoc);
	printf("i=%d,start=%d,end=%d:linecnt=%d, lines_per_map=%d\n",i,start,end,linecnt,lines_per_map);
	fflush(stdout);
	fflush(fp);
	assert((linecnt % lines_per_map) == ((1+colors) % lines_per_map));
      }
    }
    
    if(HasTrueSite){
      for(int c = 0; c < colors; c++){
	if(DEBUG) assert(pmap->truesiteL[c] && pmap->truesiteR[c]);
	fprintf(fp,"%s", VQ_TRUESITE[c][BNXVersion]);
	for(int k = 1; k <= numsite[c]; k++)
	  fprintf(fp,"\t%lld\t%lld", pmap->truesiteL[c][k],pmap->truesiteR[c][k]);
	fprintf(fp,"\n");
	linecnt++;
      }
    }

    if(MapSNR){
      for(int c = 0; c < colors; c++){
	if(DEBUG) assert(pmap->SNR[c]);
	fprintf(fp,"%s", VQ_SNR[c][BNXVersion]);
	for(register int k = 1; k <= numsite[c]; k++)
	  fprintf(fp,"\t%0.3f"/* NEW11 */, pmap->SNR[c][k]);
	if(BNXVersion <= 0)
	  fprintf(fp,"\t0.0");
	fprintf(fp,"\n");
	linecnt++;
      }
      if(DEBUG/* HERE >=2 */) assert((linecnt % lines_per_map) == ((1+2*colors) % lines_per_map));
    }

    if(MapIntensity){
      //      if(DEBUG) assert(BNXVersion >= 1);
      for(int c = 0; c < colors; c++){
	if(DEBUG) assert(pmap->Intensity[c]);
	fprintf(fp,"%s", VQ_INTENSITY[c][BNXVersion]);
	for(register int k = 1; k <= numsite[c]; k++)
	  fprintf(fp,"\t%0.3f" /* NEW11 */, pmap->Intensity[c][k]);
	if(BNXVersion <= 0)
	  fprintf(fp,"\t0.0");
	fprintf(fp,"\n");
	linecnt++;
      }
      if(DEBUG/* HERE >=2 */ && MapSNR) assert((linecnt % lines_per_map) == ((1+3*colors) % lines_per_map));
    }
    if(MapStitched){
      for(int c = 0; c < colors; c++){
	if(DEBUG) assert(pmap->stitch[c]);
	fprintf(fp,"%s", VQ_STITCH[c][BNXVersion]);
	for(register int k = 1; k <= pmap->numsite[c]+1; k++)
	  fprintf(fp,"\t%d", pmap->stitch[c][k]);
	fprintf(fp,"\n");
	linecnt++;
      }
    }
    if(MapStitchLoc){
      if(DEBUG) assert(BNXVersion >= 1);
      for(int c = 0; c < colors; c++){
	if(DEBUG) assert(pmap->stitchLocation[c]);
	fprintf(fp,"%s", VQ_STITCHLOC[c][BNXVersion]);
	for(int k = 0; k < pmap->NumStitch[c]; k++)
	  fprintf(fp,"\t%0.3f", pmap->stitchLocation[c][k]);
	fprintf(fp,"\n");
	linecnt++;
      }
    }
    if(MapPSFWidth){
      if(DEBUG) assert(BNXVersion >= 1);
      for(int c = 0; c < colors; c++){
	if(DEBUG) assert(pmap->PSFWidth[c]);
	fprintf(fp,"%s", VQ_PSFWIDTH[c][BNXVersion]);
	for(int k = 1; k <= numsite[c]; k++){
	  //	  assert(pmap->PSFWidth[c][k] > 0.0 && pmap->PSFWidth[c][k] < 50000.0);
	  fprintf(fp,"\t%0.1f", pmap->PSFWidth[c][k]);
	}
	fprintf(fp,"\n");
	linecnt++;
      }
    }
    if(MapImageLoc){
      if(DEBUG) assert(BNXVersion >= 1);
      for(int c = 0; c < colors; c++){
	if(DEBUG) assert(pmap->ImageFOV[c] && pmap->ImageX[c] && pmap->ImageY[c]);
	fprintf(fp,"%s", VQ_IMAGELOC[c][BNXVersion]);
	for(int k = 1; k <= numsite[c]; k++){
	  if(DEBUG>=2){
            #pragma omp critical 
            {
              if(!(pmap->ImageFOV[c][k] <= MaxImageFOV)){
                printf("pmap->id=%lld,c=%d,k=%d/%d:pmap->ImageFOV[c][k]=%d,MaxImageFOV=%d\n",pmap->id,c,k,numsite[c],pmap->ImageFOV[c][k],MaxImageFOV);
                fflush(stdout);
                fflush(fp);
		assert(pmap->ImageFOV[c][k] <= MaxImageFOV);
              }
            }
          }
	  fprintf(fp,"\t%d\t%0.2f\t%0.2f", pmap->ImageFOV[c][k],pmap->ImageX[c][k],pmap->ImageY[c][k]);
	}
	fprintf(fp,"\n");
	linecnt++;
      }
    }
    
    if(DEBUG && !(linecnt == lines_per_map * (i+1-start))){
      #pragma omp critical
      {
	printf("While writing map lines to %s.bnx:\n",basename);
	printf("BNXVersion=%d:",BNXVersion);
	if(BNXVersion >= 1)
	  printf("colors=%d,HasTrueSite=%d,MapSNR=%d,MapIntensity=%d,MapStitched=%d,MapStitchLoc=%d,MapPSFWidth=%d,MapImageLoc=%d\n",
		 colors,HasTrueSite,MapSNR,MapIntensity,MapStitched,MapStitchLoc,MapPSFWidth,MapImageLoc);
	printf("i=%d,start=%d,end=%d:linecnt=%d(expected %d),lines_per_map=%d\n",
	       i,start,end,linecnt,lines_per_map * (i+1-start),lines_per_map);
	fflush(stdout);
	assert(linecnt == lines_per_map * (i+1-start));
      }
    }
  }

  if(!pFP){
    (void)FILEclose(fp);
    if(DEBUG) assert(write_buffer != NULL);
    free(write_buffer);
  } else {
    if(VERB>=2 && linecnt > 0){
      #pragma omp critical
      {
        printf("    Wrote %d map lines to %s.bnx:FP[%d]= %p, buffer=%p, tid=%d\n",linecnt, basename, g, pFP[g], pwrite_buffer[g], tid);
        fflush(stdout);
      }
    }

#if 0 // NOTE : only used for testing or until MIC bug with handling of exess open file handles is fixed
    if(VERB>=2){
      #pragma omp critical
      {
        printf("tid=%d:calling FILEclose(FP[%d]=%p,fp=%p)\n",tid,g,pFP[g],fp);
	fflush(stdout);
      }
    }

    #pragma omp critical(pFPlock) /* critical section required to insure read/update is not just in this thread's register */
    {
      (void)FILEclose(fp, 1);
      pFP[g] = (FILE *)-1;
      if(DEBUG) assert(pwrite_buffer[g] != NULL);
      free(pwrite_buffer[g]);// WAS128 delete [] pwrite_buffer[g];
      pwrite_buffer[g] = NULL;
    }
#endif
  }
}  
