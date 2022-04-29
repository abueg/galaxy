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

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/output.cpp 11060 2020-05-21 23:16:28Z tanantharaman $");


/* sort doubles in increasing order */
static int doubleInc(double *p1, double *p2)
{
  return (p1[0] > p2[0]) ? 1 : (p1[0] < p2[0]) ? -1 : 0;
}

/* Output <prefix>_r.cmap
   Uses global alignment[] array to compute local coverage profile.
   If !ForceOverwrite, exit if file exists */
void output_csites(int refstart, int refend, char *prefix, int split_maps, Calign **alignment, size_t *numalign_start, size_t *numalign_end, bool full)
{
  FILE *fp;

  if(strstr(prefix,"/dev/null"))
    return;

  // HERE HERE : multithread output of _r.cmap files (see output of _q.cmap in refalign.cpp)

  /* output a .cmap file with the condensed reference (for each reference contig) */
  char cfilename[PATH_MAX];

  if(DEBUG && !(refend - refstart + 1 == numrefmaps)){
    printf("output_csites: refstart=%d, refend=%d, numrefmaps=%d\n", refstart,refend, numrefmaps);
    fflush(stdout);
    assert(refend - refstart + 1 == numrefmaps);
  }

  if(VERB>=2){
    printf("output_csites:colors=%d,usecolor=%d,refstart=%d,refend=%d,prefix=%s,SplitRef=%d,CMapID=%lld:split_maps=%d\n",colors,usecolor,refstart,refend,prefix,SplitRef,CMapID,split_maps);
    fflush(stdout);
  }

  double NoCoverage[MAXCOLOR], TotalLen[MAXCOLOR];
  for(int c = 0; c < colors; c++)
    NoCoverage[c] = TotalLen[c] = 0.0;

  for(int rid = refstart; rid <= refend; rid++){
    int start = rid;
    int end = rid;
    if(!split_maps){/* output only a single file */
      rid = refend;
      start = refstart;
      end = refend;
      if(full)
	sprintf(cfilename,"%s_full_r.cmap", prefix);
      else
	sprintf(cfilename,"%s_r.cmap", prefix);
    } else {
      if(full)
	sprintf(cfilename,"%s_contig%lld_full_r.cmap", prefix, refmap[rid]->id);
      else
	sprintf(cfilename,"%s_contig%lld_r.cmap", prefix, refmap[rid]->id);      
    }

    if(checkFile(cfilename))
      continue;

    if(VERB){
      printf("Generating %s (condensed reference)\n",cfilename);
      /*      printf("Generating %s (condensed reference):SplitRef=%d,numrefmaps=%d,num_reffiles=%d,CMapID=%lld,split=%d,rid=%d,start=%d,end=%d,refstart=%d,refend=%d\n",
	      cfilename,SplitRef,numrefmaps,num_reffiles,CMapID,split_maps,rid,start,end,refstart,refend);*/
      fflush(stdout);
    }

    char *origfilename = strdup(cfilename);
    if(1){/* write to .tmp file first, then rename it later */
      int len = strlen(cfilename);
      sprintf(&cfilename[len], ".tmp");
    }
    char *tmpfilename = strdup(cfilename);

    if((fp = fopen(cfilename,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open tmp file %s for writing:errno=%d:%s\n", cfilename,eno,err);
      exit(1);
    }

    /* write out commandline */
    printversion(fp);

    /* output key header lines */
    fprintf(fp,"# CMAP File Version:\t0.1\n");
    if(RefOrigin)
      fprintf(fp,"# Original RefMap: %s\n", RefOrigin);

    fprintf(fp,"# Label Channels:\t%d\n", colors);
    
    if(/* WAS156 */0 && usecolor >= 2 && colors == 1){/* NEW71 : This is incorrect since ref map Nickase[] values will always match query Nickase colors AFTER applying usecolor (if needed). 
							 Also if Query had 2 colors that were reduced to 1 by -usecolor, this function is called before query output and hence before switching back to 2 colors */
      if(refmap[refstart]->Nickase[0]){
	if(VERB/* HERE HERE >=2 */){
	  printf("rid=%d,refstart=%d,refend=%d:c=%d,colors=%d:refmap[refstart]->Nickase[c]= %s\n",rid,refstart,refend,0,colors,refmap[refstart]->Nickase[0]);
	  printf("rid=%d,refstart=%d,refend=%d:c=%d,colors=%d:refmap[refstart]->Nickase[c]= %s\n",rid,refstart,refend,usecolor-1,colors,refmap[refstart]->Nickase[usecolor-1]);
	  fflush(stdout);
	}
	  
	fprintf(fp,"# Nickase Recognition Site %d:\t%s\n", usecolor, refmap[refstart]->Nickase[usecolor - 1]);
      }
    } else {
      for(int c = 0; c < colors; c++){
	if(refmap[refstart]->Nickase[c]){
	  if(VERB/* HERE HERE >=2*/){
	    printf("rid=%d,refstart=%d,refend=%d:c=%d,colors=%d:refmap[refstart]->Nickase[c]= %s\n",rid,refstart,refend,c,colors,refmap[refstart]->Nickase[c]);
	    fflush(stdout);
	  }
	  
	  fprintf(fp,"# Nickase Recognition Site %d:\t%s\n", c+1, refmap[refstart]->Nickase[c]);
	}
      }
    }

    int MaskPresent = 0;
    for(int i = start; i <= end; i++)
      if(refmap[i]->Mask[0]){
	MaskPresent = 1;
	break;
      }

    const char *OutlierString = (CmapChimQuality >= 2) ? " & OutlierFrac" : "";

    fprintf(fp,"# Number of Consensus Maps:\t%d\n", end-start+1 /* WAS156 refend-refstart+1 */);
    if(FRAGCOV_FIX && FRAGCOV_FIX_CSITES && !CovNorm){
      fprintf(fp,"# StdDev & Coverage%s refer to the interval between the current site and the next site\n",OutlierString);
      if(!COVERAGE_FIX)
	fprintf(fp,"# Coverage does NOT include the ends of Query Map beyond the first or last aligned site\n");
      else
	fprintf(fp,"# Coverage includes any ends of Query Map beyond the first or last aligned site (excluding endoutliers)\n");
    } else 
      fprintf(fp,"# StdDev%s refers to the interval between the current site and the next site\n",OutlierString);

    fprintf(fp,"#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence");
    if(CmapChimQuality >= 1){
      fprintf(fp,"\tChimQuality");
      if(CmapChimQuality >= 2)
	fprintf(fp,"\tSegDupL\tSegDupR\tFragileL\tFragileR\tOutlierFrac\tChimNorm");
    }
    if(MaskPresent)
      fprintf(fp,"\tMask");
    if(Tracksites)
      fprintf(fp,"\tOriginLeft\tOriginRight");
    if(MapSNR){
      fprintf(fp,"\tGmeanSNR\tlnSNRsd");
      if(CmapSNR)
	fprintf(fp,"\tSNRcount\tSNR ...");
    }
    fprintf(fp,"\n");

    if(floatcov)
      fprintf(fp,"#f int\tfloat\tint\tint\tint\tfloat\tfloat\tfloat\tfloat");
    else
      fprintf(fp,"#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint");
    if(CmapChimQuality >= 1){
      fprintf(fp,"\tfloat");
      if(CmapChimQuality >= 2)
	fprintf(fp,"\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat");
    }
    if(MaskPresent)
      fprintf(fp,"\tHex");
    if(Tracksites)
      fprintf(fp,"\tint\tint");
    if(MapSNR){
      fprintf(fp,"\tfloat\tfloat");
      if(CmapSNR)
	fprintf(fp,"\tint\tfloat ...");
    }
    fprintf(fp,"\n");

    if(DEBUG) assert(refmap == YYmap);
    if(DEBUG) assert(Gmap == XXmap);

    for(int refid = start; refid <= end; refid++){
      Cmap *pmap = refmap[refid];
      int *NN = pmap->numsite;

      /* compute coverage profile and SNR mean,sd etc from alignments to this reference */
      double *Coverage[MAXCOLOR];/* site based coverge */
      double *ICoverage[MAXCOLOR];/* interval based coverage (used to determined unmapped region */
      double *CoverageNorm[MAXCOLOR];/* if CovNorm : used to compute map length normalized Coverage */
      double *Occurence[MAXCOLOR];
      // HERE HERE HERE : add arrays to compute SD of label intervals
      double *SNRsum[MAXCOLOR];
      double *SNRsumsq[MAXCOLOR];
      int *SNRcnt[MAXCOLOR];
      int *SNRmax[MAXCOLOR];
      double **SNR[MAXCOLOR];
      for(int c = 0; c < colors; c++){
	Coverage[c] = new double[NN[c]+1];/* Coverage[i=1..N] (based on site i) */
	ICoverage[c] = new double[NN[c]+1];/* Coverage[i=0..N] (based on interval i..i+1) */
	if(CovNorm)
	  CoverageNorm[c] = new double[NN[c]+1];
	Occurence[c] = new double[NN[c]+1];/* Occurence[1..N] weight sum (site based) */
	for(register int t = 1; t <= NN[c]; t++)
	  Coverage[c][t] = Occurence[c][t] = 0.0;
	for(register int t = 0; t <= NN[c]; t++)
	  ICoverage[c][t] = 0.0;
	if(CovNorm)
	  for(register int t = 1; t <= NN[c]; t++)
	    CoverageNorm[c][t] = 0.0;
	if(MapSNR){
	  if(CmapSNR){
	    SNRcnt[c] = new int[NN[c]+1];/* SNR count per site */
	    SNRmax[c] = new int[NN[c]+1];/* SNR count max allocation per site */
	    SNR[c] = new double *[NN[c]+1];/* SNR values array per site */
	  }
	  SNRsum[c] = new double[NN[c]+1];
	  SNRsumsq[c] = new double[NN[c]+1];
	  for(register int t = 1; t <= NN[c]; t++){
	    SNRsum[c][t] = SNRsumsq[c][t] = 0.0;
	    if(CmapSNR){
	      SNR[c][t] = new double[SNRmax[c][t] = 64];
	      SNRcnt[c][t] = 0;
	    }
	  }
	}
      }

      int coverage_trim = (REPLACE_COV ? COVERAGE_TRIM : 0);

      if(VERB>=2){
	printf("Appending refmap[%d], id=%lld, N=%d to %s (coverage_trim= %d)\n",refid,pmap->id,NN[0],cfilename,coverage_trim);
	printf("\t Computing Occurence & Coverage based on alignment[%lu .. %lu]\n",numalign_start[refid],numalign_end[refid]-1);
	fflush(stdout);
      }

      Calign **Lalignment = 0;
      size_t align_start = numalign_start[refid];
      size_t align_end = numalign_end[refid];
      size_t start = align_start;
      size_t end = align_end;
      size_t MGcnt = 0;
      if(MultiMatches){
	// first count how many total matchgroups there are
	for(size_t i = start; i < end; i++){
	  Calign *p = alignment[i];
	  if(DEBUG) assert(p->mapid1 == refid);
	  if(!p)
	    continue;
	  int U = p->numpairs;
	  if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || U < AlignedSiteThreshold2)
		      : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || U < AlignedSiteThreshold)))
	    continue;

	  if(DEBUG>=1+RELEASE && RefSplit && !(p->logPV2 >= p->logPV)){
	    printf("alignment[%lu]: xid= %lld, yid=%lld, or=%d:logPV=%0.6f,logPV2=%0.6f,score=%0.6f,np=%d\n",
		   i,XXmap[p->mapid2]->id,YYmap[p->mapid1]->id,p->orientation,p->logPV,p->logPV2,p->score,p->numpairs);
	    fflush(stdout);
	    assert(p->logPV2 >= p->logPV);
	  }

	  if(RefSplit ? p->logPV2 <= LogPvThreshold2 : p->logPV <= LogPvThreshold)
	    continue;

	  if(DEBUG>=1+RELEASE && !extendSplit && !RefSplit && colors == 1 ) assert(p->Malign && p->multicnt > 0);

	  MGcnt += p->Malign ? p->multicnt : 1;
	}

	Lalignment = new Calign*[MGcnt];
	start = 0;
	end = 0;

	for(size_t i = align_start; i < align_end; i++){
	  Calign *p = alignment[i];
	  if(!p)
	    continue;
	  int U = p->numpairs;
	  if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		      : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold)))
	    continue;
	  
	  if(!p->Malign){
	    //	    if(DEBUG>=1+RELEASE) assert(extendSplit || RefSplit);
	    p->chimpair = 0;
	    Lalignment[end++] = p;
	  } else {
	    for(int t = 0; t < p->multicnt; t++){
	      Calign *q = p->Malign[t];
	      if(DEBUG) assert(0 <= q->mapid1 && q->mapid1 < numY);
	      if(DEBUG) assert(0 <= q->mapid2 && q->mapid2 < numX);
	      if(DEBUG) assert(q->mapid1 == p->mapid1);
	      if(DEBUG) assert(q->mapid2 == p->mapid2);
	      if(DEBUG>=1+RELEASE && RefSplit) assert(q->logPV2 >= q->logPV);
	      
	      q->chimpair = (t==0) ? 0 : -1;
	      Lalignment[end++] = q;
	    }	    
	  }
	}
	if(DEBUG) assert(end-start <= MGcnt);
      }
      
      if(!Lalignment){/* make a local copy so we can re-order the alignments without affecting the calling function */
	Lalignment = new Calign*[end - start];
	for(size_t i = start; i < end; i++)
	  Lalignment[i-start] = alignment[i];
	end -= start;
	start = 0;
      }

      for(size_t i = start; i < end; i++){
	Calign *p = Lalignment[i];
	if(DEBUG && p && !(p->mapid1==refid)){
	  printf("refid=%d(%d..%d):i=%llu(%llu..%llu),numaligns=%llu:p->mapid1=%d,p->mapid2=%d\n",
		 refid,refstart,refend,(unsigned long long)i,(unsigned long long)numalign_start[refid],
		 (unsigned long long)numalign_end[refid]-1,(unsigned long long)numaligns,p->mapid1,p->mapid2);
	  for(int r = 0; r < numrefmaps; r++)
	    printf("refid=%d:numalign_start=%llu,numalign_end=%llu,numsites=%d\n",
		   r,(unsigned long long)numalign_start[r],(unsigned long long)numalign_end[r],refmap[r]->numsite[0]);
	  fflush(stdout);
	  assert(p->mapid1 == refid);
	}
	/*	if(!p || p->numpairs <= 0)
	  continue;
	if(p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold)
	continue; */
	if(!p)
	  continue;
	int U = p->numpairs;
	if(U<=0 || (RefSplit ? ((!REFSPLIT_FIX && p->score <= ScoreThreshold2) || p->logPV2 <= LogPvThreshold2 || U < AlignedSiteThreshold2) 
		    : ((!REFSPLIT_FIX && p->score <= ScoreThreshold) || p->logPV <= LogPvThreshold || U < AlignedSiteThreshold)))
	  continue;

	Cmap *Xmap = XXmap[p->mapid2];
	int *MM = Xmap->numsite;
	double Xlen = Xmap->site[0][MM[0]+1];
	Cmap *origmap = Xmap;
	int left = 0;
	//	int right = M+1;
	if(Xmap->origmap){
	  while(origmap->origmap){
	    Xmap->origmap = origmap = origmap->origmap;
	    Xmap->left[0] += max(0,origmap->left[0]-1);
	    Xmap->right[0] += max(0,origmap->left[0]-1);
	  }
	  left = max(0,Xmap->left[0] - 1);
	  //	  right = Xmap->origmap->numsite[0]+1;
	}

	if(colors > 1)
	  p++;

	for(int c = 0; c < colors; c++){
	  int K = p[c].numpairs;
	  int LI = p[c].sites1[0];
	  if(DEBUG) assert(1 <= LI && LI <= NN[c]);
	  int LK = LI - (resSD[c] > 0 ? p[c].sitesK1[0] : 0);
	  if(DEBUG) assert(1 <= LK && LK <= NN[c]);
	  int RI = p[c].sites1[K-1];
	  if(DEBUG) assert(1 <= p[c].sites1[K-1] && p[c].sites1[K-1] <= NN[c]);
	  if(coverage_trim > 0){
	    LK += coverage_trim;
	    RI -= coverage_trim;
	  } else if(COVERAGE_FIX){/* include ref labels overlapped by query map ends */
	    if(DEBUG &&  !(0 <= p[c].Lij1 && p[c].Lij1 <= LK && RI <= p[c].Rij1 && p[c].Rij1 <= NN[c]+1)){ 
	      printf("alignment[i=%lu]:mapid1=%d(id=%lld),mapid2=%d(id=%lld),or=%d,c=%d:MM[c]=%d,NN[c]=%d,Lij1=%d,Rij1=%d,LI=%d,LK=%d,RI=%d,RK=%d\n",
		     i,p->mapid1,refmap[p->mapid1]->id,p->mapid2,Gmap[p->mapid2]->id,p->orientation,c,MM[c],NN[c],p[c].Lij1,p[c].Rij1,LI,LK,RI,RI-p[c].sitesK1[K-1]);
	      fflush(stdout);
	      assert(0 <= p[c].Lij1 && p[c].Lij1 <= LK);
	      assert(RI <= p[c].Rij1 && p[c].Rij1 <= NN[c]+1);
	    }
	    LK = max(1,p[c].Lij1);
	    if(DEBUG && !(RI <= p[c].Rij1 && p[c].Rij1 <= NN[c]+1))
	    RI = min(NN[c],p[c].Rij1);
	  }

	  if(VERB>=3){
	    printf("\t\t alignment[%lu]: xid= %lld, M= %d, Xlen= %0.4f: logPV= %0.2f, np= %d, WT= %0.6f, Coverage[%d][%d..%d] += WT, Occurrence[%d][%d..%d] += WT\n",
		   i, Xmap->id, MM[0],Xlen,p->logPV, p->numpairs, p->mapWT,c,LK,RI,c,p[c].sites1[0] - (resSD[c] > 0 ? p[c].sitesK1[0] : 0), p[c].sites1[K-1]);
	    fflush(stdout);
	  }

	  if(DEBUG) assert(0 <= p->mapWT && p->mapWT <= 1.0);
	  for(int k = LK; k <= RI; k++)
	    Coverage[c][k] += p->mapWT /* WAS47 1.0 */;
	  for(int k = LK; k < RI; k++)
	    ICoverage[c][k] += p->mapWT /* WAS47 1.0 */;
	  if(CovNorm)
	    for(int k = LK; k <= RI; k++)
	      CoverageNorm[c][k] += Xlen;
	  for(int k = 0; k < K; k++){
	    int R = p[c].sites1[k];
	    int L = R - (resSD[c] > 0 ? p[c].sitesK1[k] : 0);
	    double wt = p->mapWT;
	    int J = p[c].sites2[k];
	    for(register int t = L; t <= R; t++){
	      Occurence[c][t] += wt;
	      if(MapSNR){
		double snr = origmap->SNR[c][p->orientation ? (MM[c] + 1-J)+left : J+left];
		double lnsnr = log(snr);

		if(CmapSNR){
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

		if(DEBUG/* HERE >=2 */ && !(isfinite(wt) && isfinite(lnsnr))){
		  int index = p->orientation ? (MM[c]+1-J)+left : J+left;
		  printf("output_csites:refid=%d(id=%lld),N=%d:alignment[%llu],mapid=%d(id=%lld):c=%d,numpairs=%d:k=%d,L=%d,R=%d,M=%d:t=%d,SNR[c][%d]=snr=%0.6e,lnsnr=%0.6e,wt=%0.6e\n",
			 refid,pmap->id,pmap->numsite[c],(unsigned long long)i,p->mapid2,Gmap[p->mapid2]->id,c,K,k,L,R,MM[c],t,index,snr,lnsnr,wt);
		  for(int t = 1; t <= MM[c]; t++)
		    printf("   SNR[c][%d]=%0.6e\n",t,origmap->SNR[c][t]);
		  
		  fflush(stdout);
		  assert(isfinite(wt));
		  assert(isfinite(lnsnr));
		}

		SNRsum[c][t] += wt * lnsnr;
		SNRsumsq[c][t] += wt * lnsnr*lnsnr;

		if(DEBUG/* HERE >=2 */) assert(isfinite(SNRsum[c][t]));
		if(DEBUG/* HERE >=2 */) assert(isfinite(SNRsumsq[c][t]));
	      }
	    }
	  }
	}
      }

      if(CovNorm){/* compute normalized coverage scaling factor CoverageNorm[c][t] : both Coverage[c][t] and Occurence[c][t] will be scaled by this factor */
	double CovLambdaInv = 1.0/CovLambda;
	for(int c = 0; c < colors; c++){
	  for(int t = 1; t <= NN[c]; t++){
	    int cov = Coverage[c][t];
	    if(!cov){
	      CoverageNorm[c][t] = 1.0;
	      continue;
	    }
	    double mean = CoverageNorm[c][t]/cov;/* mean size of maps contributing to Coverage[c][t] */
	    CoverageNorm[c][t] = exp((mean-MinLen-CovLambda)*CovLambdaInv);
	  }
	}
      }

      /* update NoCoverage and TotalLen */
      for(int c = 0; c < colors; c++){
	TotalLen[c] += pmap->site[c][NN[c]+1];
	for(int t = 1; t < NN[c]; t++)
	  if(!ICoverage[c][t]){
	    /* scan forward to find end of zero coverage region */
	    int u = t+1;
	    for(;u < NN[c]; u++)
	      if(Coverage[c][u])
		break;
	    NoCoverage[c] += pmap->site[c][u] - pmap->site[c][t];
	    if(VERB>=2){
	      printf("Unmapped Region refid=%d:%0.3f to %0.3f kb (Sum=%0.3f kb)\n",refid,pmap->site[c][t],pmap->site[c][u],NoCoverage[c]);
	      fflush(stdout);
	    }
	    t = u;
	  }
      }

      /* interleave colors/channels so sites are in increasing order */
      int newtotsites = 0;
      for(register int c = 0; c < colors; c++)
	newtotsites += pmap->numsite[c];
      Cmerge *pn = new Cmerge[newtotsites+1];
      int cnt = 0;
      for(int c = 0; c < colors; c++){
	for(register int j = 1;  j <= pmap->numsite[c]; j++){
	  pn[cnt].color = c;
	  pn[cnt].newid = j;
	  pn[cnt].site = pmap->site[c][j];
	  cnt++;
	}
      }
      assert(cnt == newtotsites);
      qsort(pn, newtotsites, sizeof(Cmerge), (intcmp *)siteinc);
      /* Add right end */
      pn[newtotsites].color = pn[newtotsites].newid = -1;
      pn[newtotsites].site = pmap->site[0][pmap->numsite[0]+1];
      if(DEBUG)
	for(int c = 1; c < colors; c++)
	  assert(fabs(pmap->site[c-1][pmap->numsite[c-1]+1] - pmap->site[c][pmap->numsite[c]+1]) < 1e-6);

      int N = newtotsites;
      for(register int k = 1; k <= N + 1;k++){
	double pos = pn[k-1].site * 1000.0; // pmap->site[0][k] * 1000.0;
	int channel = pn[k-1].color + 1;// (k <= N ? 1 : 0);
	int c = channel - 1;
	int newid = pn[k-1].newid;
	
	double cov = 1.0, occ = 0.0;
	if(k <= N){
	  cov = (FRAGCOV_FIX && FRAGCOV_FIX_CSITES && !CovNorm) ? ICoverage[c][newid] : Coverage[c][newid];
	  occ = Occurence[c][newid];
	  if(CovNorm){
	    double scale = CoverageNorm[c][newid];
	    cov *= scale;
	    occ *= scale;
	  }
	  // WAS  cov = max(1.0,cov);
	}
	
	if(floatcov)
	  fprintf(fp, "%lld\t%0.1f\t%d\t%d\t%d\t%0.1f\t%0.1f\t%0.1f\t%0.1f",
		  pmap->id, // CMapID
		  pmap->site[0][NN[0]+1]*1000.0, // ContigLength
		  N, // NumSites
		  k, // SiteID
		  channel, // LabelChannel
		  pos, // Position
		  0.0, // StdDev
		  (k <= N) ? cov : 1.0, // Coverage
		  (k <= N) ? occ : 1.0 // Occurence
		  );
	else
	  fprintf(fp, "%lld\t%0.1f\t%d\t%d\t%d\t%0.1f\t%0.1f\t%d\t%d",
		  pmap->id, // CMapID
		  pmap->site[0][NN[0]+1]*1000.0, // ContigLength
		  N, // NumSites
		  k, // SiteID
		  channel, // LabelChannel
		  pos, // Position
		  0.0, // StdDev
		  (k <= N) ? (int)floor(cov+0.5) : 1, // Coverage
		  (k <= N) ? (int)floor(occ+0.5) : 1 // Occurence
		  );

	if(CmapChimQuality >= 1){
	  if(pmap->ChimQuality[0]){
	    if(DEBUG && k <= N) assert(isfinite(pmap->ChimQuality[c][newid]) && -1.0 <= pmap->ChimQuality[c][newid] && pmap->ChimQuality[c][newid] <= 100.001f);
	    fprintf(fp, "\t%0.2f", (k <= N) ? pmap->ChimQuality[c][newid] : 0.0);
	  } else
	    fprintf(fp, "\t%0.2f", 0.0);
	  if(CmapChimQuality >= 2){
	    if(pmap->SegDupL[0]){
	      if(DEBUG && k <= N) assert(isfinite(pmap->SegDupL[c][newid]) && -1.0 <= pmap->SegDupL[c][newid] && pmap->SegDupL[c][newid] <= 100.001f);
	      fprintf(fp, "\t%0.2f", (k <= N) ? pmap->SegDupL[c][newid] : 0.0);
	    } else
	      fprintf(fp, "\t%0.2f", 0.0);
	    if(pmap->SegDupR[0]){
	      if(DEBUG && k <= N) assert(isfinite(pmap->SegDupR[c][newid]) && -1.0 <= pmap->SegDupR[c][newid] && pmap->SegDupR[c][newid] <= 100.001f);
	      fprintf(fp, "\t%0.2f", (k <= N) ? pmap->SegDupR[c][newid] : 0.0);
	    } else
	      fprintf(fp, "\t%0.2f", 0.0);

	    if(pmap->FragileEndL[0]){
	      if(DEBUG && k <= N) assert(isfinite(pmap->FragileEndL[c][newid]) && 0.0 <= pmap->FragileEndL[c][newid] /* && pmap->FragileEndL[c][newid] <= 1000.001f */);
	      fprintf(fp, "\t%0.2f", (k <= N) ? pmap->FragileEndL[c][newid] : 0.0);
	    } else
	      fprintf(fp, "\t%0.2f", 0.0);
	    if(pmap->FragileEndR[0]){
	      if(DEBUG && k <= N) assert(isfinite(pmap->FragileEndR[c][newid]) && 0.0 <= pmap->FragileEndR[c][newid] /* && pmap->FragileEndR[c][newid] <= 1000.001f */);
	      fprintf(fp, "\t%0.2f", (k <= N) ? pmap->FragileEndR[c][newid] : 0.0);
	    } else
	      fprintf(fp, "\t%0.2f", 0.0);

	    if(pmap->OutlierFrac[0]){
	      if(DEBUG && k <= N) assert(isfinite(pmap->OutlierFrac[c][newid]) && 0.0 <= pmap->OutlierFrac[c][newid] && pmap->OutlierFrac[c][newid] <= 100.001f);
	      fprintf(fp, "\t%0.2f", (k <= N) ? pmap->OutlierFrac[c][newid] : 0.0);
	    } else
	      fprintf(fp, "\t%0.2f", 0.0);

	    if(pmap->ChimNorm[0]){
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
	    printf("Writing to %s, id=%lld, k=%d/%d: Mask= %lx\n",cfilename,pmap->id,k,N,Mask);
	    fflush(stdout);
	  }
	  fprintf(fp, "\t%lx", Mask);
	}

	if(Tracksites){
	  if(DEBUG && k <= N) assert(pmap->truesiteL[c] && pmap->truesiteR[c]);	  
	  fprintf(fp, "\t%lld\t%lld", (k <= N) ? pmap->truesiteL[c][newid] : -1, (k <= N) ? pmap->truesiteR[c][newid] : -1);
	}

	if(MapSNR){
	  if(k <= N){
	    double wt = max(0.001,Occurence[c][newid]);
	    double sum = SNRsum[c][newid];
	    double sumsq = SNRsumsq[c][newid];
	    sum /= wt;
	    sumsq /= wt;
	    if(DEBUG && !isfinite(sum)){
	      printf("non-finite number: sumsq=%.15g sum=%.15g wt=%.15g\n", sumsq, sum, wt);
	      fflush(stdout);
	      assert(isfinite(sum));
	    }
	    double sd = sqrt(max(0.0, (sumsq - sum*sum)));
	    if(DEBUG && !isfinite(sd)) {
	      printf("non-finite number: sd=%g from sumsq=%.15g sum=%.15g wt=%.15g\n", sd, sumsq, sum, wt);
	      fflush(stdout);
	      assert(isfinite(sd));
	    }
	    if(DEBUG) assert(isfinite(exp(sum)));

	    fprintf(fp, "\t%0.4f\t%0.4f",
		    exp(sum), // GmeanSNR
		    sd   // lnSNRsd
		    );

	    if(CmapSNR){
	      int cnt = SNRcnt[c][newid];
	      double *pSNR = SNR[c][newid];

	      qsort(pSNR,cnt,sizeof(double),(intcmp*)doubleInc);

	      fprintf(fp,"\t%d",cnt);
	      for(int u = 0; u < cnt; u++)
		fprintf(fp,"\t%0.4f",pSNR[u]);
	    }
	  } else { /* right end */
	    fprintf(fp, "\t%0.4f\t%0.4f",
		    0.0, // SNRmean
		    0.0  // SNRsd
		    );
	    if(CmapSNR)
	      fprintf(fp,"\t0");
	  }
	}
	fprintf(fp, "\n");
      }
      for(int c = 0; c < colors; c++){
	delete [] Coverage[c];
	delete [] ICoverage[c];
	if(CovNorm)
	  delete [] CoverageNorm[c];
	delete [] Occurence[c];
	if(MapSNR){
	  delete [] SNRsum[c];
	  delete [] SNRsumsq[c];
	  if(CmapSNR){
	    delete [] SNRcnt[c];
	    delete [] SNRmax[c];
	    for(register int t = 1; t <= NN[c]; t++)
	      delete [] SNR[c][t];
	    delete [] SNR[c];
	  }
	}
      }
      delete [] pn;
      delete [] Lalignment;
    }
    FILEclose(fp);

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
  
  if(VERB){/* display Total length of reference and reference region without any coverage from query maps */
    for(int c = 0; c < colors; c++)
      if(TotalLen[c] > 0.0)
	printf("color=%d:Total Reference Length= %0.3f Kb, Unmapped Length= %0.3f Kb (%0.1f%%)\n",c+1,TotalLen[c],NoCoverage[c], NoCoverage[c]*100.0/TotalLen[c]);
      else
	printf("color=%d:Total Reference Length= %0.3f Kb, Unmapped Length= %0.3f Kb\n",c+1,TotalLen[c],NoCoverage[c]);
    fflush(stdout);
  }
  if(VERB){
    printf("Completed output of %d _r.cmap files: cum wall time= %0.6f\n", split_maps ? (refend-refstart+1) : 1, wtime());
    fflush(stdout);
  }
}

extern double zscore(double pvalue);

/* output error parameters from refalign iterations */
void output_err(char *basename, Cparameters *parameters, int iterations, double *MappingRatePV, double CoverageMult, int MIN_PV, int MAX_PV)
{
  if(strstr(output_prefix,"/dev/null"))
    return;

  char filename[PATH_MAX];
  if(CMapID > 0)
    sprintf(filename,"%s_id%lld.err",basename,CMapID);
  else
    sprintf(filename,"%s.err",basename);

  int origForceOverwrite = ForceOverwrite;
  if(giter2 >= 1)
    ForceOverwrite = 1;

  if(giter2 <= 0 && checkFile(filename)){
    ForceOverwrite = origForceOverwrite;
    return;
  }

  if(VERB){
    printf("Generating %s (error file)\n",filename);
    fflush(stdout);
  }

  if(colors > 2){
    fprintf(stderr,"output_err: colors=%d not supported\n",colors);
    exit(1);
  }

  int origcolors = colors;

  if(usecolor && refcolors >= 2){/* output 2-color error parameters : only one color's parameters will be changed (based on usecolor) */
    if(DEBUG) assert(colors == 1);
    if(usecolor >= 2){/* swap parameters for color "usecolor" with parameters for 1st color */
      if(VERB){
	printf("Swapping colors 1 and %d due to -usecolor %d including parameters[0..%d]\n",usecolor,usecolor,iterations);
	fflush(stdout);
      }
      swapparam(usecolor,parameters,iterations+1);
    }
    colors = 2;
  }

  if(DEBUG) assert(iterations <= RefRepeats);

  FILE *fp;
  if(giter2 > 0){/* append to file */
    if((fp = fopen(filename,"a"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for appending .err file:errno=%d:%s\n",filename,eno,err);
      exit(1);
    }
    fprintf(fp,"# set %d/%d with %d iterations\n", giter2 + 1, RefRepeats2, RefRepeats);

  } else {/* start new file */
    if((fp = fopen(filename,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for writing .err file:errno=%d:%s\n",filename,eno,err);
      exit(1);
    }
  
    /* write out commandline */
    printversion(fp);

    if(ResBins[0] > 0)
      fprintf(fp,"# See %s.errbias for bias parameters\n", output_prefix);
    if(MappingRatePV && MAX_PV)
      fprintf(fp,"# See %s.maprate for final mapping rate distribution\n", output_prefix);
    fprintf(fp,"# Label Channels:\t%d\n", colors);
  }

  if(colors==2){
    if(giter2 <= 0){
      fprintf(fp, "Iteration\tFP(/100kb)\tFN(rate)\tsf\tsd\tbpp\tres(pixels)\tMaps\tLog10LR(/Maps)\tGoodMaps\tlog10LR(/GoodMaps)\tbppSD\tFPrate\tsr\tse\tLabelDensity(/100kb)\tresSD(pixels)\tmres(x500bp)\tmresSD(x500bp)\tminSNR");
      if(outlierRate)
	fprintf(fp,"\tOutlierRate\tEndoutlierRate");
      if(perLabelScore)
	fprintf(fp,"\tLog10LR(/Label)");
      fprintf(fp, "\tLabelChannel\tMAmean(kb)\tMAsf\tMAsd");
      fprintf(fp,"\n");
    }

    for(int i = 0; i <= iterations; i++){
      for(int c = 0; c < colors; c++){
	fprintf(fp,"%d\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.2f\t%0.3f\t%d\t%0.6f\t%d\t%0.6f\t%0.2f\t%0.6f\t%0.6f\t%0.3f\t%0.3f\t%0.6f\t%0.3f\t%0.3f\t%0.3f",
		i,
		parameters[i].FP[c],      /* FalsePositives(per100kb) */
		parameters[i].FN[c],      /* FalseNegativeRate */
		parameters[i].SF[c],      /* SF[c] */
		parameters[i].SD[c],      /* SD[c] */
		parameters[i].PixelLen*1000.0, /* BasesPerPixel */
		parameters[i].res[c],        /* reference resolution in pixels */
		nummaps /* parameters[i].mapcnt*/,     /* number of Aligned Maps (usually same as total maps) */
		parameters[i].logLR/log(10.0), /* log10(LikelihoodRatio) per Aligned Map */
		parameters[i].ATmapcnt,        /* number of Maps that aligned above all 3 thresholds */
		parameters[i].ATlogLR/log(10.0),/* log10(LikelihoodRatio) per Map that aligned above all 3 threshols */
		parameters[i].PixelLenSD*1000.0, /* SD(BasesPerPixel) */
		parameters[i].FPfrac[c],    /* FalsePositive Rate */
		parameters[i].SR[c],        /* Relative Error */
		parameters[i].SE[c],        /* Resolution Error */
		parameters[i].LabelDensity[c], /* Total Label Density (per 100kb) */
		parameters[i].resSD[c],
		parameters[i].mres,
		parameters[i].mresSD,
		parameters[i].minSNR[c]);
	if(outlierRate)
	  fprintf(fp,"\t%0.5f\t%0.5f", parameters[i].outlierRate, parameters[i].EndoutlierRate);
	if(perLabelScore)
	  fprintf(fp,"\t%0.6f", parameters[i].ATlogLR/(log(10.0)*max(1,parameters[i].sitecnt - parameters[i].ATmapcnt)));
	fprintf(fp,"\t%d\t%0.6f\t%0.6f\t%0.6f",
		c+1,                                /* LabelChannel */
		c ? parameters[i].MA_mean : 0.0,    /* MAbias(Kb) */
		c ? parameters[i].SF[2] : 0.0,      /* SiteMA(Kb) */
		c ? parameters[i].SD[2] : 0.0);      /* ScalingMA(Kb^1/2) */
	fprintf(fp,"\n");
      }
    }
  } else {/* single color */
    if(giter2 <= 0){
      //      fprintf(fp, "Iteration\tFP(/100kb)\tFNrate\tSiteSD(Kb)\tScalingSD(Kb^1/2)\tbpp\tres(pixels)\tMaps\tLog10LR(/Maps)\tGoodMaps\tlog10LR(/GoodMaps)\tbppSD\tFPrate\tRelativeSD\tResolutionSD\tLabelDensity(/100kb)\tresSD\tmres\tmresSD");
      fprintf(fp,"# MapRate is by length\n");
      fprintf(fp, "Iter\tFP/100k\tFNrate \tsf     \tsd     \tbpp   \tres(px)\tMaps  \tMapRate\tAligned\tlog10LR\tbppSD\tFPrate \tsr     \tse   \tYLab/100k\tresSD\tmres \tmresSD\tminSNR\tXlab/100k");
      if(outlierRate)
	fprintf(fp,"\tOutlierRate\tEndoutlierRate");
      if(perLabelScore)
	fprintf(fp,"\tLog10LR(/Label)");
      fprintf(fp,"\n");
    }
    for(register int i = 0; i <= iterations; i++){
      fprintf(fp, "%d\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%6.2f\t%0.3f\t%-6d\t%7.5f\t%-6d\t%7.4f\t%5.2f\t%0.5f\t%0.5f\t%0.3f\t%8.5f \t%0.3f\t%0.3f\t%0.3f \t%0.3f\t%8.5f ",
	      i,
	      parameters[i].FP[0],      /* FalsePositives(per100kb) */
	      parameters[i].FN[0],      /* FalseNegative Rate */
	      parameters[i].SF[0],      /* SiteErr(Kb) */
	      parameters[i].SD[0],      /* ScalingError(Kb^1/2) */
	      parameters[i].PixelLen*1000.0, /* BasesPerPixel */
	      parameters[i].res[0],        /* reference resolution in pixels */
	      nummaps/* parameters[i].mapcnt */,     /* number of Aligned Maps (usually same as total maps) */
	      parameters[i].sumX/max(0.001,parameters[i].totlen), /* Effective Coverage Fraction by aligned length (including internal outliers) */
	      parameters[i].ATmapcnt,        /* number of Maps that aligned above all 3 thresholds */
	      parameters[i].ATlogLR/log(10.0),/* log10(LikelihoodRatio) per Map that aligned above all 3 threshols */
	      parameters[i].PixelLenSD*1000.0, /* SD(BasesPerPixel) */
	      parameters[i].FPfrac[0],  /* FalsePositive Rate */
	      parameters[i].SR[0],     /* RelativeError */
	      parameters[i].SE[0],     /* ResolutionError */
	      parameters[i].LabelDensity[0], /* Aligned Ref Label Density (per 100kb) */
	      parameters[i].resSD[0],     
	      parameters[i].mres,
	      parameters[i].mresSD,
	      parameters[i].minSNR[0] ,
	      parameters[i].LabelDensityX[0] /* Aligned Qry label Density (per 100kb) */
	      );
      if(outlierRate)
	fprintf(fp,"\t%0.5f\t%0.5f", parameters[i].outlierRate, parameters[i].EndoutlierRate);
      if(perLabelScore)
	fprintf(fp,"\t%0.6f", parameters[i].ATlogLR/(log(10.0)*max(1,parameters[i].sitecnt - parameters[i].ATmapcnt)));
      fprintf(fp,"\n");
    }
  }
  FILEclose(fp);
  fp = NULL;

  if(MappingRatePV && MAX_PV){  /* output the mapping rates */
    sprintf(filename,"%s.maprate", basename);

    if(checkFile(filename)){
      ForceOverwrite = origForceOverwrite;
      return;
    }
    
    if(VERB){
      printf("Generating %s (mapping rate file)\n", filename);
      fflush(stdout);
    }

    if((fp = fopen(filename,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for writing mapping rate:errno=%d:%s\n",filename,eno,err);
      exit(1);
    }

    printversion(fp);

    fprintf(fp,"# LogPV\tMappingRate\tCoverage\n");
    for(int i = MIN_PV; i <= MAX_PV; i++)
      fprintf(fp,"%0.1f\t%0.5f\t%0.3f\n",(double)i, MappingRatePV[i], MappingRatePV[i]*CoverageMult);

    FILEclose(fp);
    fp = NULL;
  }

  /* output the bias parameters as a <prefix>.errbias file */
  if(ResBins[0] + (colors >= 2 ? ResBins[1] : 0) > 0){
    sprintf(filename,"%s.errbias%d",basename,giter2);

    if(checkFile(filename)){
      ForceOverwrite = origForceOverwrite;
      return;
    }

    if(VERB){
      printf("Generating %s (error bias file for giter2=%d)\n",filename,giter2);
      fflush(stdout);
    }
    
    if((fp = fopen(filename,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for writing .errbias file:errno=%d:%s\n",filename,eno,err);
      exit(1);
    }

    printversion(fp);

    if(colors > 1)
      fprintf(fp,"# Number of PieceWise Linear segments = %d %d\n", ResBins[0], ResBins[1]);
    else
      fprintf(fp,"# Number of PieceWise Linear segments = %d\n", ResBins[0]);
    fprintf(fp,"# Number of Iterations = %d\n", iterations);
    for(int c = 0; c < colors; c++){
      fprintf(fp,"#color%d:\t",c+1);
      for(int i = 0; i <= iterations; i++)
	fprintf(fp,"Size%d\tBias%d\t", i,i);
      fprintf(fp,"\n");
    
      int maxResBins = 0;
      for(int i = 0; i <= iterations; i++)
	maxResBins = max(maxResBins, parameters[i].ResBins[c]);
      if(DEBUG) assert(maxResBins <= RESBINS);

      for(int Bin = 0; Bin <= maxResBins; Bin++){
	for(int i = 0; i <= iterations; i++)
	  if(Bin >= parameters[i].ResBins[c])
	    fprintf(fp,"%0.3f\t%0.3f\t", 0.0,0.0);
	  else
	    fprintf(fp,"%0.3f\t%0.3f\t", parameters[i].resbiasX[c][Bin], parameters[i].resbias[c][Bin]);
	fprintf(fp,"\n");
      }
    }

    FILEclose(fp);
    fp = NULL;
  }

  if(VERB>=2){
    printf("giter2=%d,RefRepeats2=%d,RefRepeats=%d\n",giter2,RefRepeats2,RefRepeats);
    fflush(stdout);
  }

  if(giter2 == RefRepeats2 - 1){  /* output the last iteration as a binary file */
    if(CMapID > 0)
      sprintf(filename,"%s_id%lld.errbin",basename,CMapID);
    else
      sprintf(filename,"%s.errbin",basename);

    if(checkFile(filename)){
      ForceOverwrite = origForceOverwrite;
      return;
    }

    if(VERB){
      printf("Generating %s (binary error file)\n",filename);
      fflush(stdout);
    }

    if((fp = fopen(filename,"wb"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for writing .errbin file:errno=%d:%s\n",filename,eno,err);
      exit(1);
    }

    int iter = (SaveRepeat >= 0 && giter2 == RefRepeats2-1) ? min(iterations,SaveRepeat) : iterations ; // - 1;
    parameters[iter].version = Cparameters_Version;
    parameters[iter].colors = colors;
    parameters[iter].MaxResBins = RESBINS;
    if(usecolor){
      parameters[iter].SF[colors] = MA_SF;
      parameters[iter].SD[colors] = MA_SD;
      parameters[iter].MA_mean = MA_mean;
    }
    if(DEBUG){
      for(int c = 0; c < colors; c++){
	if(DEBUG && parameters[iter].ResBins[c] > 0 &&  ! (parameters[iter].ResBins[c] <= RESBINS && parameters[iter].maxresbias >= parameters[iter].resbiasX[c][parameters[iter].ResBins[c]] - 1e-6)){
	  printf("c=%d:parameters[%d]:ResBins[c]=%d/%d,maxresbias=%0.3f,resbiasX[ResBins]=%0.3f\n",iter,
		 c,parameters[iter].ResBins[c],RESBINS,parameters[iter].maxresbias,parameters[iter].resbiasX[c][parameters[iter].ResBins[c]]);
	  for(int t = 0; t <= parameters[iter].ResBins[c]; t++)
	    printf("    t=%d:resbiasX[c][t]=%0.3f,resbias[c][t]=%0.6f\n",t,parameters[iter].resbiasX[c][t],parameters[iter].resbias[c][t]);
	  fflush(stdout);
	  assert(parameters[iter].maxresbias >= parameters[iter].resbiasX[c][parameters[iter].ResBins[c]] - 1e-6);
	}
      }
    }

    if(VERB>=2){
      printf("writing %s : parameters[%d].MaxResBins=%d,RESBINS=%d,sizeof(Cparameters)=%lu,",filename, iter,parameters[iter].MaxResBins, RESBINS, sizeof(Cparameters));
      printf("mapcnt=%d,aligned=%d\n",parameters[iter].mapcnt,parameters[iter].ATmapcnt);
      for(int c = 0; c < colors; c++)
	printf("c=%d:  FP=%0.8f, FN=%0.8f, sf=%0.8f, sd=%0.8f,sr=%0.8f,se=%0.8f,res=%0.8f,resSD=%0.8f\n", 
	       c,parameters[iter].FP[c],parameters[iter].FN[c],parameters[iter].SF[c],parameters[iter].SD[c],parameters[iter].SR[c],parameters[iter].SE[c],parameters[iter].res[c],parameters[iter].resSD[c]);
      fflush(stdout);
    }

    size_t i = fwrite(&parameters[iter],sizeof(Cparameters),1,fp);
    if(i != 1){
      char *err = strerror(errno);
      printf("failed to write to %s (i=%llu):%s\n",filename,(unsigned long long)i,err);
      exit(1);
    }
    FILEclose(fp);
  }

  if(usecolor >= 2 && colors >= 2){/* swap parameters for color "usecolor" with parameters for 1st color */
    if(VERB){
      printf("Swapping colors 1 and %d due to -usecolor %d including parameters[0..%d]\n",usecolor,usecolor,iterations);
      fflush(stdout);
    }
    if(DEBUG) assert(refcolors >= 2);
    swapparam(usecolor,parameters, iterations+1);
  }
  colors = origcolors;
  ForceOverwrite = origForceOverwrite;
}

void output_refalign(int refstart, int refend, char *basename)
{
  if(strstr(basename,"/dev/null"))
    return;

  if(DEBUG) assert(colors == 1);

  char filename[PATH_MAX];
  size_t blen = 10*1024*1024;/* matches scif_io buffer size */
  char *write_buffer =  (char *)malloc(blen);
  if(!write_buffer){
    printf("output_refalign:malloc(%llu) failed\n",(unsigned long long)blen);
    exit(1);
  }

  if(DEBUG) assert(refend - refstart + 1 == numrefmaps);
  
  for(int refid = refstart; refid <= refend; refid++){
    int start = refid;
    int end = refid;
    if(SplitRef || CMapID <= 0){/* output only a single file */
      refid = refend;
      start = refstart;
      end = refend;
      sprintf(filename,"%s.map",basename);
    } else
      sprintf(filename,"%s_contig%lld.map", basename, refmap[refid]->id);

    if(checkFile(filename))
      continue;

    if(VERB){
      printf("Generating %s\n",filename);
      fflush(stdout);
    }

    FILE *fp;
    if((fp = fopen(filename,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for writing .map file:errno=%d:%s\n",filename,eno,err);
      exit(1);
    }

    setbuffer(fp, write_buffer, blen);
  
    /* write out commandline */
    printversion(fp);

    //    fprintf(fp,"# Source: %s\n",basename);
    fprintf(fp,"# N Channels: %d\n", colors);
    if(!NoSplit)
      fprintf(fp,"#WARNING: Molecule sites may not be correct for maps with more than 1 chimeric split due to SPLIT_LIMIT=0\n");

    //    fprintf(fp,"Software version: %s\n",SVN_ID);
    if(scale_output)
      fprintf(fp,"MappedMoleculeId\tMoleculeId\tMoleculeIndex\tContigId\tScore\tZscore   \tDirection\tStartLocation\tEndLocation\tStartMatchLocation\tEndMatchLocation\tDetectedLabelCount\tTruePositiveLabelCount\tFalsePositiveLabelCount\tFalseNegativeLabelCount\tLog10Pvalue\tLeftEndChim\tRightEndChim\tBPP");
    else
      fprintf(fp,"MappedMoleculeId\tMoleculeId\tMoleculeIndex\tContigId\tScore\tZscore   \tDirection\tStartLocation\tEndLocation\tStartMatchLocation\tEndMatchLocation\tDetectedLabelCount\tTruePositiveLabelCount\tFalsePositiveLabelCount\tFalseNegativeLabelCount\tLog10Pvalue\tLeftEndChim\tRightEndChim");

    if(VERB){
      double pvalue = pow(10.0,-max(0.0,LogPvThreshold));
      register double Zscore = zscore(pvalue);
      printf("Pvalue threshold=%0.6e corresponds to Z-score threshold of %0.6f\n",
	     pvalue, Zscore);
      fflush(stdout);
    }

    /* locate longest aligned map */
    register int maxsites = 0;
    for(register size_t i = numalign_start[start]; i < numalign_end[end]; i++){
      register Calign *p = alignment[i];
      if(!p)
	continue;
      register Cmap *pmap = Gmap[p->mapid2];
      if(pmap->numsite[0] > maxsites)
	maxsites = pmap->numsite[0];
    }

    for(register int i = 1; i <= maxsites; i++)
      if(resSD[0] > 0.0)
	fprintf(fp,"\tNanoLabelID%d\tRefLabelIDleft%d\tRefLabelIDright%d",i,i,i);
      else
	fprintf(fp,"\tNanoLabelID%d\tRefLabelID%d",i,i);
    fprintf(fp,"\n");
  
    int aligncnt = 0;
    for(register size_t i = numalign_start[start]; i < numalign_end[end]; i++){
      register Calign *p = alignment[i];
      if(!p || p->numpairs <= 0)
	continue;
      if(p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold)
	continue;

      int rid = p->mapid1;
      int mid = p->mapid2;
      if(DEBUG) assert(mid < nummaps);
      if(DEBUG) assert(start <= rid && rid <= end && rid < numrefmaps);

      register Cmap *Xmap = Gmap[mid];

      Cmap *origmap = Xmap;
      while(origmap->origmap)
	origmap = origmap->origmap;

      if(BestRef && origmap->align->mapid1 != rid)
	continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */

      aligncnt++;
      register Cmap *Ymap = refmap[rid];
      register double pvalue = pow(10.0,-max(0.0,p->logPV));
      register double Zscore = zscore(pvalue);

      fprintf(fp,"%d\t%lld\t%d\t%lld",aligncnt,Xmap->id,mid,Ymap->id);
      fprintf(fp,"\t%0.3f", p->score);
      fprintf(fp,"\t%0.6f", Zscore);
      fprintf(fp,"\t%c", p->orientation ? '-' : '+');/* Direction */
      
      register FLOAT *Y = Ymap->site[0];
      register FLOAT *X = Xmap->site[0];
      register int N = Ymap->numsite[0];
      register int M = Xmap->numsite[0];

      register int left = 0;
      register int right = M+1;
      if(Xmap->origmap){
	Cmap *origmap = Xmap;
	while(origmap->origmap){
	  Xmap->origmap = origmap = origmap->origmap;
	  Xmap->left[0] += max(0,origmap->left[0]-1);
	  Xmap->right[0] += max(0,origmap->left[0]-1);
	}
	left = max(0,Xmap->left[0] - 1);
	right = Xmap->origmap->numsite[0]+1;
      }

      register int numpairs = p->numpairs;
      if(DEBUG)	assert(numpairs >= 1);
      int LI = p->sites1[0];
      int LK = resSD[0] > 0.0 ? p->sitesK1[0] : 0;
      int LJ = p->sites2[0];
      int RI = p->sites1[numpairs-1];
      int RK = resSD[0] > 0.0 ? p->sitesK1[numpairs-1] : 0;
      int RJ = p->sites2[numpairs-1];
      if(DEBUG) assert(LI >= 1 && LI <= N);
      if(DEBUG) assert(RI >= 1 && RI <= N);
      if(DEBUG) assert(LJ >= 1 && LJ <= M);
      if(DEBUG) assert(RJ >= 1 && RJ <= M);
      int Lend = -p->Lend;
      int Rend = -p->Rend;
      if(Xmap->origmap){/* chimeric map fragment : check which end is the local alignment end */
	if(Xmap->left[0] <= 0 && Xmap->right[0] <= Xmap->origmap->numsite[0]){
	  if(p->orientation){
	    Lend = 2;
	  } else {
	    Rend = 2;
	  }
	}
	if(Xmap->left[0] > 0 && Xmap->right[0] > Xmap->origmap->numsite[0]){
	  if(p->orientation){
	    Rend = 2;
	  } else {
	    Lend = 2;
	  }
	}
      }
      
      double StartLocation = Y[LI] - (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]);
      double EndLocation = Y[RI] + (p->orientation ? X[M+1-RJ] : X[M+1]-X[RJ]);
      double StartMatchLocation = Y[LI];
      double EndMatchLocation = Y[RI];
      
      if(VERB >=2 && Xmap->id==85 && Ymap->id == 1){
	printf("Xmap->id=%lld:mapid=%d,LI=%d,RI=%d,N=%d,LJ=%d,RJ=%d,M=%d,left=%d,right=%d,Lend=%d,Rend=%d,orientation=%d\n",
	       Xmap->id,p->mapid2,LI,RI,N,LJ,RJ,M,left,right,Lend,Rend,p->orientation);
	printf("  StartLocation=%0.3f,EndLocation=%0.3f,StartMatchLocation=%0.3f,EndMatchLocation=%0.3f\n",
	       StartLocation,EndLocation,StartMatchLocation,EndMatchLocation);
	fflush(stdout);
      }

      int Ksum = 0;
      if(resSD[0] > 0.0)
	for(register int k = 0; k < numpairs; k++)
	  Ksum += p->sitesK1[k];

      fprintf(fp,"\t%lld", (long long int)floor(StartLocation*1000.0 + 0.5));/* StartLocation */
      fprintf(fp,"\t%lld", (long long int)floor(EndLocation*1000.0 + 0.5));/* EndLocation */
      fprintf(fp,"\t%lld", (long long int)floor(StartMatchLocation*1000.0 + 0.5));/* StartMatchLocation */
      fprintf(fp,"\t%lld", (long long int)floor(EndMatchLocation*1000.0 + 0.5));/* EndMatchLocation */
      fprintf(fp,"\t%d", M);/* DetectedLabelCount */
      fprintf(fp,"\t%d", numpairs); /* TruePositiveLabelCount */
      fprintf(fp,"\t%d", M-numpairs - (Lend >= 2 ? LJ-1 : 0) - (Rend >= 2 ? M-RJ : 0));/* FalsePositiveLabelCount */
      fprintf(fp,"\t%d", (Rend >= 2 ? RI : p->Rij1) - (Lend >= 2 ? LI-LK : p->Lij1) + 1 - numpairs - Ksum); /* FalseNegativeLabelCount */
      fprintf(fp,"\t%0.3f",p->logPV);/* Log10Pvalue */
      fprintf(fp,"\t%d\t%d", Lend-1, Rend-1);    /* left and right alignment end types (-1 = Chr End, 0 = Global, 1 = Chimeric/Local) */
      if(scale_output)
	fprintf(fp, "\t%0.2f", Xmap->incscale * PixelLen * 1000.0);

      register int I,K,J,H,T,G;

      /* output left end of alignment left of aligned pair (LI,LJ) */
      if(Lend < 2){
	for(J = 1; J < LJ; J++){
	  if(resSD[0] > 0.0)
	    fprintf(fp,"\t%4d\t%4d\t%4d",(p->orientation ? (M+1-J)+left : J+left),-1,-1);/* NanoLabelID, RefLabelIDleft(-1),RefLabelIDright(-1) */
	  else
	    fprintf(fp,"\t%4d\t%4d",(p->orientation ? (M+1-J)+left : J+left),-1);/* NanoLabelID, RefLabelID(-1) */
	}
      }

      /* output leftmost aligned pair (LI,LK,LJ) */
      if(resSD[0] > 0.0)
	fprintf(fp,"\t%4d\t%4d\t%4d", (p->orientation ? (M+1-LJ)+left : LJ+left), LI-LK,LI);  /* NanoLabelID, RefLabelIDleft,RefLabelIDright */
      else
	fprintf(fp,"\t%4d\t%4d", (p->orientation ? (M+1-LJ)+left : LJ+left), LI);  /* NanoLabelID, RefLabelID */
	      
      I = G = LI;
      K = T = LK;
      J = H = LJ;
      for(register int k = 1; k < numpairs;k++, G = I, T = K, H = J){
	I = p->sites1[k];
	K = resSD[0] > 0.0 ? p->sitesK1[k] : 0;
	J = p->sites2[k];
	/* output false positive sites between H and J (based on linear interpolation) */
	for(register int j = H+1; j < J; j++)
	  if(resSD[0] > 0.0)
	    fprintf(fp,"\t%4d\t%4d\t%4d",  (p->orientation ? (M+1-j)+left : j+left),-1,-1);/* NanoLabelID, RefLabelIDleft, RefLabelIDright */      
	  else
	    fprintf(fp,"\t%4d\t%4d",  (p->orientation ? (M+1-j)+left : j+left),-1);/* NanoLabelID, RefLabelID */      
      
	/* output aligned pair (I,J) (OR (I-K,I,J) if resSD > 0.0)*/
	if(resSD[0] > 0.0)
	  fprintf(fp,"\t%4d\t%4d\t%4d", (p->orientation ? (M+1-J)+left : J+left), I-K, I);  /* NanoLabelID, RefLabelIDleft, RefLabelIDright */      
	else
	  fprintf(fp,"\t%4d\t%4d", (p->orientation ? (M+1-J)+left : J+left), I);  /* NanoLabelID, RefLabelID */      
      }
      if(DEBUG)
	assert(I==RI && J == RJ && K == RK);

      /* output right end of alignment right of alignment (RI,LJ) */
      if(Rend < 2){
	for(J=RJ+1;J <= M; J++){
	  if(resSD[0] > 0.0)
	    fprintf(fp,"\t%4d\t%4d\t%4d",(p->orientation ? (M+1-J)+left : J+left), -1, -1); /* NanoLabelID, RefLabelIDleft, RefLabelIDright */      
	  else
	    fprintf(fp,"\t%4d\t%4d",(p->orientation ? (M+1-J)+left : J+left), -1); /* NanoLabelID, RefLabelID */      
	}
      }

      fprintf(fp,"\n");
    }

    FILEclose(fp);
  } /* refid = refstart ... refend */

  free(write_buffer);
}

/** Details for single aligned sites pair for output. Two-color. */
class Calignpair {
public:
  int color;
  int site1,siteIK1,site2;/**< site1 is -1 for unaligned sites in map2 */
  double loc;/**< sum of distance from left end for both maps */
};

#if 0
static int CalignpairLocInc(Calignpair *p1, Calignpair *p2)
{
  if(DEBUG>=2) assert(p1->loc > -1e9);
  if(DEBUG>=2) assert(p2->loc > -1e9);
  return (p1->loc > p2->loc) ? 1 : (p1->loc < p2->loc) ? -1 : 0;
}
#endif


#if 0 // obsolete .map file no longer supported
/** Two-color version of .map output */
void output_refalign2(int refstart, int refend, char *basename)
{
  if(strstr(basename,"/dev/null"))
    return;

  if(DEBUG) assert(colors==2);

  char filename[PATH_MAX];
  strcpy(filename,basename);
  /* no need to strip suffix : for default basename this was already done in RefAligner.cpp */
  if(1){
    int i = strlen(filename);
    strcpy(&filename[i],".map");
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

    printf("failed to open file %s for writing .map file:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  /* write out commandline */
  printversion(fp);

  //  fprintf(fp,"# Source: %s\n",basename);
  fprintf(fp,"# N Channels: %d\n", colors);
  if(!NoSplit)
    fprintf(fp,"#WARNING: Nanomap sites may not be correct for maps with more than 1 chimeric split due to SPLIT_LIMIT=0\n");

  //  fprintf(fp,"Software version: %s\n",SVN_ID);
  if(scale_output)
    fprintf(fp,"MappedMoleculeId\tMoleculeId\tMoleculeIndex\tContigId\tScore\tConfidence\tDirection\tStartLocation\tEndLocation\tStartMatchLocation1\tStartMatchLocation2\tEndMatchLocation1\tEndMatchLocation2\tDetectedLabelCount1\tDetectedLabelCount2\tTruePositiveLabelCount1\tTruePositiveLabelCnt2\tFalsePositiveLabelCount1\tFalsePositiveLabelCount2\tFalseNegativeLabelCount1\tFalseNegativeLabelCount1\tLog10Pvalue\tLeftEndChim1\tLeftEndChim2\tRightEndChim1\tRightEndChim2\tBPP");
  else
    fprintf(fp,"MappedMoleculeId\tMoleculeId\tMoleculeIndex\tContigId\tScore\tConfidence\tDirection\tStartLocation\tEndLocation\tStartMatchLocation1\tStartMatchLocation2\tEndMatchLocation1\tEndMatchLocation2\tDetectedLabelCount1\tDetectedLabelCount2\tTruePositiveLabelCount1\tTruePositiveLabelCnt2\tFalsePositiveLabelCount1\tFalsePositiveLabelCount2\tFalseNegativeLabelCount1\tFalseNegativeLabelCount1\tLog10Pvalue\tLeftEndChim1\tLeftEndChim2\tRightEndChim1\tRightEndChim2");

  if(VERB){
    double pvalue = pow(10.0,-max(0.0,LogPvThreshold));
    register double Zscore = zscore(pvalue);
    printf("Pvalue threshold=%0.6e corresponds to Z-score threshold of %0.6f\n",
	   pvalue, Zscore);
    fflush(stdout);
  }

  /* locate longest aligned map */
  register int maxsites = 0;
  for(register size_t i = refstart; i < refend; i++){
    register Calign *p = alignment[i];
    if(!p)
      continue;
    register Cmap *pmap = Gmap[p->mapid2];
    int numsite = 0;
    for(int c = 0; c < colors; c++)
      numsite += pmap->numsite[0];
    if(numsite > maxsites)
      maxsites = numsite;
  }

  for(register int i = 1; i <= maxsites; i++)
    if(colors > 1){
      if(DEBUG) assert(resSD[0] > 0.0);
      fprintf(fp,"\tChannelID%d\tNanoLabelID%d\tRefLabelIDleft%d\tRefLabelIDright%d",i,i,i,i);
    } else if(resSD[0] > 0.0)
      fprintf(fp,"\tNanoLabelID%d\tRefLabelIDleft%d\tRefLabelIDright%d",i,i,i);
    else
      fprintf(fp,"\tNanoLabelID%d\tRefLabelID%d",i,i);
  fprintf(fp,"\n");
  
  int aligncnt = 0;
  for(register size_t i = 0; i < numaligns; i++){
    register Calign *p = alignment[i];
    if(!p || !p->numpairs)
      continue;
    if(p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold)
      continue;

    int rid = p->mapid1;
    int mid = p->mapid2;
    if(DEBUG) assert(refstart <= rid && rid < refend);
    if(DEBUG) assert(mid < nummaps);

    register Cmap *Xmap = Gmap[mid];
    Cmap *origmap = Xmap;
    while(origmap->origmap)
      origmap = origmap->origmap;
    if(BestRef && origmap->align->mapid1 != rid)
      continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */

    aligncnt++;
    register Cmap *Ymap = refmap[rid];
    register double pvalue = pow(10.0,-max(0.0,p->logPV));
    register double Zscore = zscore(pvalue);

    fprintf(fp,"%d\t%lld\t%d\t%lld",aligncnt,Xmap->id,mid,Ymap->id);
    fprintf(fp,"\t%0.3f", p->score);
    fprintf(fp,"\t%0.6f", Zscore);
    fprintf(fp,"\t%c", p->orientation ? '-' : '+');/* Direction */

    register FLOAT **Y = Ymap->site;
    register FLOAT **X = Xmap->site;
    register int *N = Ymap->numsite;
    register int *M = Xmap->numsite;
    int left[MAXCOLOR],right[MAXCOLOR];
    for(int c = 0; c < colors; c++){
      left[c] = 0;
      right[c] = M[c]+1;

      if(Xmap->origmap){
	left[c] = max(0,Xmap->left[c] - 1);
	right[c] = Xmap->origmap->numsite[c]+1;
      }
    }

    if(DEBUG) assert(p->numpairs >= 1);

    double StartLocation[2],EndLocation[2],StartMatchLocation[2],EndMatchLocation[2];
    int LI[2],LK[2],LJ[2],RI[2],RJ[2],Lend[2],Rend[2],Ksum[2];//RK[2]

    p++;

    for(int c = 0; c < colors; c++){
      register int U = p[c].numpairs;
      LI[c] = p[c].sites1[0];
      LK[c] = p[c].sitesK1[0];
      LJ[c] = p[c].sites2[0];
      RI[c] = p[c].sites1[U-1];
      //      RK[c] = p[c].sitesK1[U-1];
      RJ[c] = p[c].sites2[U-1];
      assert(LI[c] >= 1 && LI[c] <= N[c]);
      assert(RI[c] >= 1 && RI[c] <= N[c]);
      assert(LJ[c] >= 1 && LJ[c] <= M[c]);
      assert(RJ[c] >= 1 && RJ[c] <= M[c]);
      Lend[c] = -p[c].Lend;
      Rend[c] = -p[c].Rend;
      if(Xmap->origmap){/* chimeric map fragment : check which end is the local alignment end */
	if(Xmap->left[c] <= 0 && Xmap->right[c] <= Xmap->origmap->numsite[c]){
	  if(p[-1].orientation){
	    Lend[c] = 2;
	  } else {
	    Rend[c] = 2;
	  }
	}
	if(Xmap->left[c] > 0 && Xmap->right[c] > Xmap->origmap->numsite[c]){
	  if(p[-1].orientation){
	    Rend[c] = 2;
	  } else {
	    Lend[c] = 2;
	  }
	}
      }

      StartLocation[c] = Y[c][LI[c]] - (p[-1].orientation ? X[c][M[c]+1] - X[c][M[c]+1-LJ[c]] : X[c][LJ[c]]);
      StartMatchLocation[c] = Y[c][LI[c]];
      EndLocation[c] = Y[c][RI[c]] + (p[-1].orientation ? X[c][M[c]+1-RJ[c]] : X[c][M[c]+1]-X[c][RJ[c]]);
      EndMatchLocation[c] = Y[c][RI[c]];

      if(VERB >= 2){
	printf("Xmap->id=%lld,c=%d:mapid=%d:LI=%d,RI=%d,N=%d,LJ=%d,RJ=%d,M=%d,left=%d,right=%d,Lend=%d,Rend=%d,orientation=%d\n",
	       Xmap->id,c,p[-1].mapid2,LI[c],RI[c],N[c],LJ[c],RJ[c],M[c],left[c],right[c],Lend[c],Rend[c],p[-1].orientation);
	printf("  StartLocation=%0.3f,EndLocation=%0.3f,StartMatchLocation=%0.3f,EndMatchLocation=%0.3f\n",
	       StartLocation[c],EndLocation[c],StartMatchLocation[c],EndMatchLocation[c]);
	fflush(stdout);
      }

      Ksum[c] = 0;
      for(register int k = 0; k < U; k++)
	Ksum[c] += p[c].sitesK1[k];
    }

    fprintf(fp,"\t%lld", (long long int)floor(min(StartLocation[0],StartLocation[1])*1000.0 + 0.5));/* StartLocation */
    fprintf(fp,"\t%lld", (long long int)floor(max(EndLocation[0],EndLocation[1])*1000.0 + 0.5));/* EndLocation */
    fprintf(fp,"\t%lld", (long long int)floor(StartMatchLocation[0]*1000.0 + 0.5));/* StartMatchLocation1 */
    fprintf(fp,"\t%lld", (long long int)floor(StartMatchLocation[1]*1000.0 + 0.5));/* StartMatchLocation2 */
    fprintf(fp,"\t%lld", (long long int)floor(EndMatchLocation[0]*1000.0 + 0.5));/* EndMatchLocation1 */
    fprintf(fp,"\t%lld", (long long int)floor(EndMatchLocation[1]*1000.0 + 0.5));/* EndMatchLocation2 */
    fprintf(fp,"\t%d", M[0]);/* DetectedLabelCount1 */
    fprintf(fp,"\t%d", M[1]);/* DetectedLabelCount2 */
    fprintf(fp,"\t%d", p[0].numpairs); /* TruePositiveLabelCount1 */
    fprintf(fp,"\t%d", p[1].numpairs); /* TruePositiveLabelCount2 */
    fprintf(fp,"\t%d", M[0]-p[0].numpairs - (Lend[0] >= 2 ? LJ[0]-1 : 0) - (Rend[0] >= 2 ? M[0]-RJ[0] : 0));/* FalsePositiveLabelCount1 */
    fprintf(fp,"\t%d", M[1]-p[1].numpairs - (Lend[1] >= 2 ? LJ[1]-1 : 0) - (Rend[1] >= 2 ? M[1]-RJ[1] : 0));/* FalsePositiveLabelCount2 */
    fprintf(fp,"\t%d", (Rend[0] >= 2 ? RI[0] : p[0].Rij1) - (Lend[0] >= 2 ? LI[0]-LK[0] : p[0].Lij1) + 1 - p[0].numpairs - Ksum[0]); /* FalseNegativeLabelCount1 */
    fprintf(fp,"\t%d", (Rend[1] >= 2 ? RI[1] : p[1].Rij1) - (Lend[1] >= 2 ? LI[1]-LK[1] : p[1].Lij1) + 1 - p[1].numpairs - Ksum[1]); /* FalseNegativeLabelCount2 */
    fprintf(fp,"\t%0.3f",p->logPV);/* Log10Pvalue */
    fprintf(fp,"\t%d\t%d\t%d\t%d", Lend[0]-1, Lend[1]-1,Rend[0]-1,Rend[1]-1);    /* left and right alignment end types for Channel 1 & 2 (-1 = Chr End, 0 = Global, 1 = Chimeric/Local) */
    if(scale_output)
      fprintf(fp, "\t%0.2f", Xmap->incscale * PixelLen * 1000.0);

    /* sort alignment pairs based on total distance from left end of both maps and both channels */
    int maxpairs = M[0]+M[1];
    Calignpair *alignpairs = new Calignpair[maxpairs];
    for(int i = 0; i < maxpairs; i++)
      alignpairs[i].loc = -1e+10;

    int pairs = 0;
    for(int c = 0; c < colors; c++){
      int I,K,J, JL = (Lend[c] < 2) ? 0 : p[c].sites2[0];
      for(int k = 0; k < p[c].numpairs; k++, JL = J){
	alignpairs[pairs].color = c;
	alignpairs[pairs].site1 = I = p[c].sites1[k];
	K = p[c].sitesK1[k];
	alignpairs[pairs].siteIK1 = I-K;
	J = p[c].sites2[k];
	alignpairs[pairs].site2 = p[-1].orientation ? M[c]+1-J : J;

	double Yik = Yc(Y[c],I,K);
	double Xj = p[-1].orientation ? X[c][M[c]+1] - X[c][M[c]+1-J] : X[c][J];
	alignpairs[pairs++].loc = Yik + Xj;

	/* save unaligned sites between JL and J */
	for(int j = JL+1; j < J; j++){
	  alignpairs[pairs].color = c;
	  alignpairs[pairs].site1 = alignpairs[pairs].siteIK1 = -1;
	  alignpairs[pairs].site2 = p[-1].orientation ? M[c]+1-j : j;
	  alignpairs[pairs++].loc = Yik + (p[-1].orientation ? X[c][M[c]+1] - X[c][M[c]+1-j] : X[c][j]);
	}
      }
      if(Rend[c] < 2){/* save unaligned sites between JL and M+1 */
	double Yik = Y[c][N[c]+1];
	for(int j = JL+1; j <= M[c]; j++){
	  alignpairs[pairs].color = c;
	  alignpairs[pairs].site1 = alignpairs[pairs].siteIK1 = -1;
	  alignpairs[pairs].site2 = p[-1].orientation ? M[c]+1-j : j;
	  alignpairs[pairs++].loc = Yik + (p[-1].orientation ? X[c][M[c]+1] - X[c][M[c]+1-j] : X[c][j]);
	}
      }
    }
    if(DEBUG) assert(pairs <= maxpairs);
    qsort(alignpairs,pairs,sizeof(Calignpair),(intcmp*) CalignpairLocInc);

    /* output alignment quads (color,Nanolabel,RefLabelLeft,RefLabelRight) */
    for(int i = 0; i < pairs; i++){
      Calignpair *pa = &alignpairs[i];
      fprintf(fp,"\t%1d\t%4d\t%4d\t%4d", pa->color+1, pa->site2, pa->siteIK1, pa->site1);/* Color, NanoLabelID, RefLabelIDleft(or -1),RefLabelIDright(or -1) */
    }
    fprintf(fp,"\n");
    
    delete [] alignpairs;
  }

  FILEclose(fp);
}
#endif // obsolete 2-color .map file output

void output_chimmaps(char *basename)
{
  if(strstr(basename,"/dev/null"))
    return;

  if(DEBUG) assert(!PairMerge);/* since minid,maxid would not be defined */

  char filename[PATH_MAX];
  strcpy(filename,basename);
  /* no need to strip suffix : for default basename this was already done in RefAligner.cpp */
  if(1){
    int i = strlen(filename);
    strcpy(&filename[i],".chim");
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

    printf("failed to open file %s for writing .chim file:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  /* write out commandline */
  printversion(fp);

  //  fprintf(fp,"# Source: %s\n",basename);
  //  fprintf(fp,"# Software Version: %s\n", SVN_ID);
  fprintf(fp,"MoleculeID\tHitCnt\tDistantHits\n");
  
  int *HitCnt = new int[startmaps+1];
  int *DistantHits = new int[startmaps+1];
  for(register int i=0; i < startmaps; i++)
    HitCnt[i] = DistantHits[i] = 0;

  for(register int i=numaligns;--i>=0;){
    register Calign *p = alignment[i];
    if(!p || p->numpairs <= 0)
      continue;
    if(p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold)
      continue;

    register Cmap *pmap = Gmap[p->mapid2];
    Cmap *origmap = pmap;
    while(origmap->origmap)
      origmap = origmap->origmap;
    if(BestRef && origmap->align->mapid1 != p->mapid1)
      continue;/* skip alignment unless it is with the reference with the best alignment score with Xmap */

    while(pmap->origmap)
      pmap = pmap->origmap;
    register int mapid = pmap->mapid;
    if(DEBUG) assert(mapid >= 0 && mapid < startmaps);

    HitCnt[mapid]++;
    /*    if(VERB && pmap->id==96 && Gmap[p->mapid2]->origmap){
      printf("pmap->id=%d:mapid=%d:chimpair=%d\n",
	     pmap->id,p->mapid2,p->chimpair);
	     fflush(stdout);
	     }*/

    if(Gmap[p->mapid2]->origmap && !p->chimpair)
      DistantHits[mapid]++;
  }
  for(register int i = 0; i < startmaps; i++){
    long long int id = Gmap[i]->id;
    int mapid = Gmap[i]->mapid;
    if(DEBUG)  assert(id >= minid && id <= maxid);
    fprintf(fp,"%lld\t%d\t%d\n", id, HitCnt[mapid],DistantHits[mapid]);
  }
  FILEclose(fp);

  delete [] HitCnt;
  delete [] DistantHits;
}

extern Cmap **YYmap,**XXmap;
extern int numY,numX;

#if 0
/* sort alignments that are above threshold before those that are below threshold and then in reducing order of logPV */
int PVdec(register Calign **pp1, register Calign **pp2)
{
  register Calign *p1 = *pp1, *p2 = *pp2;
  register int pass1 = (p1<=0) ? -1 : (p1->score > ScoreThreshold && p1->logPV > LogPvThreshold && p1->numpairs >= AlignedSiteThreshold) ? 1 : 0;
  register int pass2 = (p2<=0) ? -1 : (p2->score > ScoreThreshold && p2->logPV > LogPvThreshold && p2->numpairs >= AlignedSiteThreshold) ? 1 : 0;

  if(DEBUG>=2 && pass1 == pass2 && p1 <= 0) assert(p2 <= 0);

  return (pass1 < pass2) ? 1 : (pass1 > pass2) ? -1 : (p1 <= 0) ? 0 : 
    (p1->logPV < p2->logPV) ? 1 : (p1->logPV > p2->logPV) ? -1 : 0;
}
#endif

void output_align(char *basename, int mapstart)
{
  if(colors > 2){
    printf("Multiple channel data for %d colors not yet implemented for .align output\n",colors);
    exit(1);
  }

  if(strstr(basename,"/dev/null"))
    return;

  char filename[PATH_MAX];
  strcpy(filename,basename);
  /* no need to strip suffix : for default basename this was already done in RefAligner.cpp */
  if(1){
    int i = strlen(filename);
    strcpy(&filename[i],".align");
  }

  if(checkFile(filename))
    return;

  char *origfilename = strdup(filename);
  if(1){/* write to .tmp file first, then rename it later so -RefineOverwrite works correctly */
    int len = strlen(filename);
    sprintf(&filename[len], ".tmp");
  }
  char *tmpfilename = strdup(filename);

  /* compact alignment array */
  size_t aligncnt = 0;
  for(size_t i = 0; i < numaligns; i++){
    Calign *p = alignment[i];
    if(!p || p->numpairs <= 0)
      continue;
    if(p->score <= ScoreThreshold || p->logPV <= LogPvThreshold || p->numpairs < AlignedSiteThreshold)
      continue;

    if(aligncnt < i){
      Calign *tmp = alignment[aligncnt];
      alignment[aligncnt] = p;
      alignment[i] = tmp;
    }
    aligncnt++;
  }
  numaligns = aligncnt;

  if(VERB){
    printf("Generating %s (%lu alignments)\n",origfilename, numaligns);
    fflush(stdout);
  }

  FILE *fp;
  if((fp = fopen(filename,align_format ? "wb" : "w"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("failed to open file %s for writing .align file:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  char *write_buffer = NULL;
  size_t blen = 10*1024*1024; // must match scif_io.h buffer size to avoid triggering bug (failed fwrite())
  if(!(write_buffer = (char *)malloc(blen))){
    printf("output_align:malloc(%llu) failed\n",(unsigned long long)blen);
    exit(1);
  }
  setbuffer(fp, write_buffer, blen);

  FILE *fpbin = NULL;
  char filenamebin[PATH_MAX];
  char *write_bufferbin = NULL;
  if(align_format){
    sprintf(filenamebin,"%s.bin",origfilename);
    if((fpbin = fopen(filenamebin,"wb"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for writing alignment binary data:errno=%d:%s\n",filenamebin,eno,err);
      exit(1);
    }
    size_t b2len = 10*1024*1024;
    if(!(write_bufferbin = (char *)malloc(b2len))){
      int eno = errno;
      char *err = strerror(eno);
      printf("output_align: malloc(%llu) failed:errno=%d:%s\n", (unsigned long long)b2len,eno,err);
      exit(1);
    }
    setbuffer(fpbin, write_bufferbin, b2len);
  }
  
  /* write out commandline */
  printversion(fp);

  fprintf(fp,"# Map File Input:\t%s", vfixx_filename[0]);
  for(register int i = 1; i < num_files;i++)
    fprintf(fp,",%s",vfixx_filename[i]);
  fprintf(fp,"\n");

  /*  fprintf(fp,"#\tScoring Method:\t%s\n",
	  PAIRSCORETYPE==1 ? "Thomas's Method" :
	  PAIRSCORETYPE==2 ? "Nguyen's Method" :
	  PAIRSCORETYPE==3 ? "John's refactored Method" :
	  PAIRSCORETYPE==4 ? "Thomas's Modified Method" : 
	  PAIRSCORETYPE==5 ? "Thomas's Refactored Method" : "Unknown");*/
  fprintf(fp,"# Number of Maps:\t%d\n", startmaps);
  fprintf(fp,"# Label Channels:\t%d\n", colors);
  if(align_format){
    fprintf(fp,"# binary format with AlignScore=%d\n", AlignScore);
    fprintf(fp,"# alignments= %lu\n",numaligns);
  } else {
    fprintf(fp,"#>0\tAlignmentID\tMol0ID\tMol1ID\tScore\tCenterOffset\tOverlap\tOrientation\tPvalueLog10\tTrueOffset\tTrueOverlapFraction\tFileID0\tFileID1\tOffset\n");
    if(colors==1){
      fprintf(fp,"#Mol0 Participating Label IDs\n");
      fprintf(fp,"#Mol1 Participating Label IDs\n");
      if(AlignScore){
	fprintf(fp,"#Aligned Segment Scores left of each label and right end\n");
	fprintf(fp,"#Aligned Segment Raw Scores without outlier adjustments\n");
      }
    } else {
      for(int c = 0; c < colors; c++){
	fprintf(fp,"#color=%d: Mol0 Participating Label IDs\n",c+1);
	fprintf(fp,"#color=%d: Mol1 Participating Label IDs\n",c+1);
	if(AlignScore){
	  fprintf(fp,"#color=%d: Aligned Segment Scores left of each label and right end\n",c+1);
	  fprintf(fp,"#color=%d: Aligned Segment Raw Scores without outlier adjustments\n", c+1);
	}
      }
    }
  }

  size_t offset = 0;/* binary file offset */
  if(align_format){
    long long *Yid_vec = new long long[numaligns];
    long long *Xid_vec = new long long[numaligns];
    int *numpairs_vec = new int[numaligns];

    /* first output Yid,Xid and numpairs as vectors to binary file */
    for(size_t i=0; i < numaligns; i++){
      Calign *p = alignment[i];
      numpairs_vec[i] = p->numpairs;
      Yid_vec[i] = YYmap[p->mapid1]->id;
      Xid_vec[i] = XXmap[p->mapid2]->id;
    }

    size_t r;

    /* output Yid vector */
    if((r = fwrite_LOOP(Yid_vec,sizeof(long long), numaligns, fpbin)) != numaligns){
      printf("fwrite_LOOP of %lu long long of Ymap->id to %s failed (ret=%lu)\n",numaligns,filenamebin,r);
      fflush(stdout);  exit(1);
    }
    offset += r*sizeof(long long);

    /* output Xid vector */
    if((r = fwrite_LOOP(Xid_vec,sizeof(long long), numaligns, fpbin)) != numaligns){
      printf("fwrite_LOOP of %llu long long of Xmap->id to %s failed (ret=%llu)\n",(unsigned long long)numaligns,filenamebin,(unsigned long long)r);
      exit(1);
    }
    offset += r*sizeof(long long);

    /* output numpairs vector */
    if((r = fwrite_LOOP(numpairs_vec,sizeof(int), numaligns, fpbin)) != numaligns){
      printf("fwrite_LOOP of %llu int(s) of numpairs to %s failed (ret=%llu)\n",(unsigned long long)numaligns,filenamebin,(unsigned long long)r);
      exit(1);
    }
    offset += r*sizeof(int);
    
    delete [] Yid_vec;
    delete [] Xid_vec;
    delete [] numpairs_vec;
  }
  size_t *offset_vec = new size_t[numaligns];/* binary file offset at start of each contig */

  if(colors==1){
    for(register size_t i=0; i < numaligns; i++){
      register Calign *p = alignment[i];

      register int K = p->numpairs;
      register Cmap *Ymap = YYmap[p->mapid1];
      register Cmap *Xmap = XXmap[p->mapid2];
      register int N = Ymap->numsite[0];
      register int M = Xmap->numsite[0];
      register FLOAT *Y = Ymap->site[0];
      register FLOAT *X = Xmap->site[0];
      register int Yshift = Ymap->origmap ? max(1,Ymap->left[0]) - 1 : 0;
      register int Xshift = Xmap->origmap ? max(1,Xmap->left[0]) - 1 : 0;
      register int LI = p->sites1[0];
      register int LJ = p->sites2[0];
      register int RI = p->sites1[K-1];
      register int RJ = p->sites2[K-1];
      double trueoffset = Xmap->startloc - Ymap->startloc;
      if(VERB>=2 && Ymap->id==69 && Xmap->id==2742){
	printf("Yid=%lld,Xid=%lld,or=%d:Xmap->startloc=%0.1f,Ymap->startloc=%0.1f,trueoffset=%0.3f kb\n",Ymap->id,Xmap->id,p->orientation,Xmap->startloc,Ymap->startloc,trueoffset);
	fflush(stdout);
      }
      double trueoverlap = (trueoffset>=0.0) ? min(Y[N+1]-trueoffset,X[M+1]) : min(X[M+1]+trueoffset,Y[N+1]);
      /* average center to center distance over all aligned sites */
      register double distanceYtoX = 0;
      /* adjust location to be relative to origmap */
      register Cmap *origYmap = Ymap->origmap ? Ymap->origmap : Ymap;
      register Cmap *origXmap = Xmap->origmap ? Xmap->origmap : Xmap;
      register FLOAT *origY = origYmap->site[0];
      register FLOAT *origX = origXmap->site[0];
      register int origN = origYmap->numsite[0];
      register int origM = origXmap->numsite[0];
      register double centerY = origY[origN+1]*0.5;
      register double centerX = origX[origM+1]*0.5;
      for(register int index = 0; index < K; index++){
	register int I = p->sites1[index];
	register int J = p->sites2[index];
	if(VERB>=3){
	  printf("Yid=%lld,Xid=%lld,or=%d,N=%d(%d),M=%d(%d),centerX=%0.4f,centerY=%0.4f,index=%d/%d,I=%d,J=%d: YtoX= %0.4f\n",
		 Ymap->id,Xmap->id,p->orientation,N,origN,M,origM,centerX,centerY,index,K,I,J,
		 (origY[I+Yshift]-centerY) + (centerX - (p->orientation ? (origX[origM+1]-origX[M+1-J+Xshift]) : origX[J+Xshift])));
	  fflush(stdout);
	}
	distanceYtoX += (origY[I+Yshift]-centerY) + (centerX - (p->orientation ? (origX[origM+1]-origX[M+1-J+Xshift]) : origX[J+Xshift]));
      }
      distanceYtoX /= K;

      double overlap = (Y[RI]-Y[LI] + (p->orientation ? X[M+1-LJ] - X[M+1-RJ] : X[RJ]-X[LJ]))*0.5;
      //      double leftEnd = (p->Lend >= -1) ? 0.0 : min(Y[LI], (p->orientation ? X[M+1] - X[M+1-LJ] : X[LJ]));
      //      double rightEnd = (p->Rend >= -1) ? 0.0 : min(Y[N+1]-Y[RI], (p->orientation ? X[M+1-RJ] : X[M+1] - X[RJ]));

      if(align_format){
	offset_vec[i] = offset;/* accumulate vector of bytes offsets into binary file, to allow multithreaded read */

	if(VERB>=3 && i < 10){
	  printf("i=%lu/%lu: offset=%lu\n", i, numaligns, offset);
	  fflush(stdout);
	}

	size_t r;

	if((r = fwrite(&p->score,sizeof(double), 1, fpbin)) != 1){
	  printf("fwrite of 1 double from alignment[%llu]->score to %s failed (ret=%llu)\n",
		 (unsigned long long)i,filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(double);

	if((r = fwrite(&p->logPV,sizeof(double), 1, fpbin)) != 1){
	  printf("fwrite of 1 double from alignment[%llu]->logPV to %s failed (ret=%llu)\n",
		 (unsigned long long)i,filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(double);

	double TrueOffset = trueoffset/PixelLen;
	if((r = fwrite(&TrueOffset,sizeof(double), 1, fpbin)) != 1){
	  printf("fwrite of 1 double from alignment[%llu] TrueOffset to %s failed (ret=%llu)\n",
		 (unsigned long long)i,filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(double);

	double TrueOverlapFraction = (Xmap->startloc||Ymap->startloc) ? max(0.0,trueoverlap)*2.0/(X[M+1]+Y[N+1]) : 0.0;
	if((r = fwrite(&TrueOverlapFraction,sizeof(double), 1, fpbin)) != 1){
	  printf("fwrite of 1 double from alignment[%llu] TrueOverlapFraction to %s failed (ret=%llu)\n",
		 (unsigned long long)i, filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(double);

	if((r = fwrite(&p->orientation,sizeof(int), 1, fpbin)) != 1){
	  printf("fwrite of 1 int from alignment[%llu]->orientation to %s failed (ret=%llu)\n",
		 (unsigned long long)i, filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(int);

	if((r = fwrite(&p->scaleID,sizeof(int), 1, fpbin)) != 1){
	  printf("fwrite of 1 int from alignment[%llu]->orientation to %s failed (ret=%llu)\n",
		 (unsigned long long)i, filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(int);

	if((r = fwrite_LOOP(&p->sites1[0],sizeof(int), p->numpairs, fpbin)) != (size_t)p->numpairs){
	  printf("fwrite_LOOP of %d int(s) from alignment[%llu]->sites1[] to %s failed (ret=%llu)\n",
		 p->numpairs,(unsigned long long)i, filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(int);

	if((r = fwrite_LOOP(&p->sites2[0],sizeof(int), p->numpairs, fpbin)) != (size_t)p->numpairs){
	  printf("fwrite_LOOP of %d int(s) from alignment[%llu]->sites1[] to %s failed (ret=%llu)\n",
		 p->numpairs,(unsigned long long)i, filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(int);

	if(AlignScore){
	  if((r = fwrite_LOOP(&p->iscore[0],sizeof(FLOAT), p->numpairs, fpbin)) != (size_t)p->numpairs){
	    printf("fwrite_LOOP of %d doubles from alignment[%llu]->iscore[] to %s failed (ret=%llu)\n",
		   p->numpairs,(unsigned long long)i, filenamebin,(unsigned long long)r);
	    exit(1);
	  }
	  offset += r*sizeof(FLOAT);

	  if((r = fwrite_LOOP(&p->outscore[0],sizeof(FLOAT), p->numpairs, fpbin)) != (size_t)p->numpairs){
	    printf("fwrite_LOOP of %d doubles from alignment[%llu]->iscore[] to %s failed (ret=%llu)\n",
		   p->numpairs,(unsigned long long)i, filenamebin,(unsigned long long)r);
	    exit(1);
	  }
	  offset += r*sizeof(FLOAT);
	}

	//      if(DEBUG>=2) fflush(fpbin);

      } else {/* text format */
	fprintf(fp,">0\t%llu\t%lld\t%lld\t%0.3f\t%0.3f\t%0.3f\t%d\t%0.3f\t%0.3f\t%0.4f\t%d\t%d\t%0.3f\n", 
		(unsigned long long)(i+1), /* AlignmentID */
		Ymap->id, /* Mol0ID */
		Xmap->id, /* Mol1ID */
		p->score, /* Score */
		distanceYtoX/PixelLen, /* CenterOffset (pixels) */
		overlap/PixelLen, /* Overlap (pixels) */
		p->orientation ? -1 : 1, /* Orientation */
		p->logPV,                 /* PvalueLog10 */
		trueoffset/PixelLen, /* TrueOffset */
		(Xmap->startloc||Ymap->startloc) ? max(0.0,trueoverlap)*2.0/(X[M+1]+Y[N+1]) : 0.0,/* TrueOverlapFraction */
		Ymap->fileid+1, /* File0ID */
		Xmap->fileid+1, /* File1ID */
		(distanceYtoX+centerY-centerX)/PixelLen /* OffSet (pixels) */
		);

	/* output Mol0 label IDs */
	for(register int k = 0; k < K; k++)
	  fprintf(fp,"%d\t", Yshift + p->sites1[k] - 1); /* Label ID (zero-based) */
	fprintf(fp,"\n");

	/* output Mol1 label IDs */
	for(register int k = 0; k < K; k++)
	  if(p->orientation)
	    fprintf(fp,"%d\t", Xshift + M - p->sites2[k]);/* Label ID (zero-based) */
	  else
	    fprintf(fp,"%d\t", Xshift + p->sites2[k] - 1);
	fprintf(fp,"\n");

	if(AlignScore){
	  /* output aligned segment scores */
	  for(register int k = 0; k < K; k++)
	    fprintf(fp,"%0.3f\t", p->iscore[k]);
	  fprintf(fp,"%0.3f\n", p->iscore[K]);
	  
	  /* output aligned segment scores without outlier correction */
	  for(register int k = 0; k < K; k++)
	    fprintf(fp,"%0.3f\t", p->outscore[k]);
	  fprintf(fp,"%0.3f\n", p->outscore[K]);
	}
      }
    } // for i = 0 .. numaligns - 1

  } else { // colors==2
    if(CmapChimQuality >= 1){
      printf("CMAP ChimQuality not yet supported for 2 colors\n");
      exit(1);
    }

    for(size_t i=0; i < numaligns; i++){
      Calign *p = &alignment[i][1];

      Cmap *Ymap = YYmap[p[-1].mapid1];
      Cmap *Xmap = XXmap[p[-1].mapid2];
      int *NN = Ymap->numsite;
      int *MM = Xmap->numsite;
      FLOAT **YY = Ymap->site;
      FLOAT **XX = Xmap->site;

      int KK[2], Yshift[2], Xshift[2], LI[2],LJ[2],RI[2],RJ[2];
      for(int c = 0; c < colors; c++){
	KK[c] = p[c].numpairs;
	Yshift[c] = Ymap->origmap ? max(1,Ymap->left[c]) - 1 : 0;
	Xshift[c] = Xmap->origmap ? max(1,Xmap->left[c]) - 1 : 0;
	LI[c] = p[c].sites1[0];
	LJ[c] = p[c].sites2[0];
	RI[c] = p[c].sites1[KK[c]-1];
	RJ[c] = p[c].sites2[KK[c]-1];
      }

      double trueoffset = Xmap->startloc - Ymap->startloc;
      double trueoverlap = (trueoffset>=0.0) ? YY[0][NN[0]+1] - trueoffset : XX[0][MM[0]+1] + trueoffset;

      /* average center to center distance over all aligned sites */
      double distanceYtoX = 0;
      /* adjust location to be relative to origmap */
      Cmap *origYmap = Ymap->origmap ? Ymap->origmap : Ymap;
      Cmap *origXmap = Xmap->origmap ? Xmap->origmap : Xmap;
      FLOAT **origYY = origYmap->site;
      FLOAT **origXX = origXmap->site;
      int *origNN = origYmap->numsite;
      int *origMM = origXmap->numsite;
      double centerY = origYY[0][origNN[0]+1]*0.5;
      double centerX = origXX[0][origMM[0]+1]*0.5;
      for(int c = 0; c < colors; c++){
	for(int index = KK[c]; --index >= 0;){
	  int I = p[c].sites1[index];
	  int J = p[c].sites2[index];
	  distanceYtoX += (origYY[c][I+Yshift[c]] - centerY) + (centerX - (p[c].orientation ? (origXX[c][origMM[c]+1] - origXX[c][MM[c]+1-J+Xshift[c]]) : origXX[c][J+Xshift[c]]));
	}
      }
      distanceYtoX /= max(1,KK[0]+KK[1]);

      //double leftEnd[2],rightEnd[2], maxOutlier = 0.0;
      double overlap[2];
      for(int c = 0; c < colors; c++){
	overlap[c] = (YY[c][RI[c]] - YY[c][LI[c]] + (p[c].orientation ? XX[c][MM[c]+1-LJ[c]] - XX[c][MM[c]+1-RJ[c]] : XX[c][RJ[c]] - XX[c][LJ[c]]))*0.5;
	//	leftEnd[c] = min(YY[c][LI[c]], (p[c].orientation ? XX[c][MM[c]+1] - XX[c][MM[c]+1-LJ[c]] : XX[c][LJ[c]]));
	//	rightEnd[c] = min(YY[c][NN[c]+1]-YY[c][RI[c]], (p->orientation ? XX[c][MM[c]+1-RJ[c]] : XX[c][MM[c]+1] - XX[c][RJ[c]]));
      }

      if(align_format){
	printf("binary .align not yet implemented for 2 colors\n");// HERE
	exit(1);

#if 0
	offset_vec[i] = offset;/* accumulate vector of bytes offsets into binary file, to allow multithreaded read */

	if(VERB>=3 && i < 10){
	  printf("i=%lu/%lu: offset=%lu\n", i, numaligns, offset);
	  fflush(stdout);
	}

	size_t r;

	if((r = fwrite(&p->score,sizeof(double), 1, fpbin)) != 1){
	  printf("fwrite of 1 double from alignment[%llu]->score to %s failed (ret=%llu)\n",
		 (unsigned long long)i,filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(double);

	if((r = fwrite(&p->logPV,sizeof(double), 1, fpbin)) != 1){
	  printf("fwrite of 1 double from alignment[%llu]->logPV to %s failed (ret=%llu)\n",
		 (unsigned long long)i,filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(double);

	double TrueOffset = trueoffset/PixelLen;
	if((r = fwrite(&TrueOffset,sizeof(double), 1, fpbin)) != 1){
	  printf("fwrite of 1 double from alignment[%llu] TrueOffset to %s failed (ret=%llu)\n",
		 (unsigned long long)i,filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(double);

	double TrueOverlapFraction = (Xmap->startloc||Ymap->startloc) ? max(0.0,trueoverlap)*2.0/(X[M+1]+Y[N+1]) : 0.0;
	if((r = fwrite(&TrueOverlapFraction,sizeof(double), 1, fpbin)) != 1){
	  printf("fwrite of 1 double from alignment[%llu] TrueOverlapFraction to %s failed (ret=%llu)\n",
		 (unsigned long long)i, filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(double);

	if((r = fwrite(&p->orientation,sizeof(int), 1, fpbin)) != 1){
	  printf("fwrite of 1 int from alignment[%llu]->orientation to %s failed (ret=%llu)\n",
		 (unsigned long long)i, filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(int);

	if((r = fwrite(&p->scaleID,sizeof(int), 1, fpbin)) != 1){
	  printf("fwrite of 1 int from alignment[%llu]->orientation to %s failed (ret=%llu)\n",
		 (unsigned long long)i, filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(int);

	if((r = fwrite_LOOP(&p->sites1[0],sizeof(int), p->numpairs, fpbin)) != p->numpairs){
	  printf("fwrite of %d int(s) from alignment[%llu]->sites1[] to %s failed (ret=%llu)\n",
		 p->numpairs,(unsigned long long)i, filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(int);

	if((r = fwrite_LOOP(&p->sites2[0],sizeof(int), p->numpairs, fpbin)) != p->numpairs){
	  printf("fwrite of %d int(s) from alignment[%llu]->sites1[] to %s failed (ret=%llu)\n",
		 p->numpairs,(unsigned long long)i, filenamebin,(unsigned long long)r);
	  exit(1);
	}
	offset += r*sizeof(int);

	if(AlignScore){
	  if((r = fwrite_LOOP(&p->iscore[0],sizeof(FLOAT), p->numpairs, fpbin)) != p->numpairs){
	    printf("fwrite of %d doubles from alignment[%llu]->iscore[] to %s failed (ret=%llu)\n",
		   p->numpairs,(unsigned long long)i, filenamebin,(unsigned long long)r);
	    exit(1);
	  }
	  offset += r*sizeof(FLOAT);

	  if((r = fwrite_LOOP(&p->outscore[0],sizeof(FLOAT), p->numpairs, fpbin)) != p->numpairs){
	    printf("fwrite of %d doubles from alignment[%llu]->iscore[] to %s failed (ret=%llu)\n",
		   p->numpairs,(unsigned long long)i, filenamebin,(unsigned long long)r);
	    exit(1);
	  }
	  offset += r*sizeof(FLOAT);
	}

	//      if(DEBUG>=2) fflush(fpbin);
#endif

      } else {/* text format */
	fprintf(fp,">0\t%llu\t%lld\t%lld\t%0.3f\t%0.3f\t%0.3f\t%d\t%0.3f\t%0.3f\t%0.4f\t%d\t%d\t%0.3f\n", 
		(unsigned long long)(i+1), /* AlignmentID */
		Ymap->id, /* Mol0ID */
		Xmap->id, /* Mol1ID */
		p[-1].score, /* Score */
		distanceYtoX/PixelLen, /* CenterOffset (pixels) */
		(overlap[0]+overlap[1])*0.5/PixelLen, /* Overlap (pixels) */
		p[-1].orientation ? -1 : 1, /* Orientation */
		p[-1].logPV,                 /* PvalueLog10 */
		trueoffset/PixelLen, /* TrueOffset */
		(Xmap->startloc||Ymap->startloc) ? max(0.0,trueoverlap)*2.0/(XX[0][MM[0]+1]+YY[0][NN[0]+1]) : 0.0,/* TrueOverlapFraction */
		Ymap->fileid+1, /* File0ID */
		Xmap->fileid+1, /* File1ID */
		(distanceYtoX+centerY-centerX)/PixelLen /* OffSet (pixels) */
		);

	for(int c = 0; c < colors; c++){
	  /* output Mol0 label IDs */
	  for(int k = 0; k < KK[c]; k++)
	    fprintf(fp,"%d\t", Yshift[c] + p[c].sites1[k] - 1); /* Label ID (zero-based) */
	  fprintf(fp,"\n");

	  /* output Mol1 label IDs */
	  for(int k = 0; k < KK[c]; k++)
	    if(p->orientation)
	      fprintf(fp,"%d\t", Xshift[c] + MM[c] - p[c].sites2[k]);/* Label ID (zero-based) */
	    else
	      fprintf(fp,"%d\t", Xshift[c] + p[c].sites2[k] - 1);
	  fprintf(fp,"\n");

	  if(AlignScore){
	    /* output aligned segment scores */
	    for(int k = 0; k < KK[c]; k++)
	      fprintf(fp,"%0.3f\t", p[c].iscore[k]);
	    fprintf(fp,"%0.3f\n", p[c].iscore[KK[c]]);
	  
	    /* output aligned segment scores without outlier correction */
	    for(register int k = 0; k < KK[c]; k++)
	      fprintf(fp,"%0.3f\t", p[c].outscore[k]);
	    fprintf(fp,"%0.3f\n", p[c].outscore[KK[c]]);
	  }
	}
      }
    } // for i = 0 .. numaligns-1
  } // colors==2

  if(align_format){
    if(VERB>=2){
      printf("Calling FILEclose(fpbin): filenamebin=%s\n",filenamebin);
      fflush(stdout);
    }
    FILEclose(fpbin);
  }

  if(align_format){/* save offset_vec[0..numaligns-1] to .align file as binary array */
    fflush(fp);
    errno = 0;
    clearerr(fp);
    size_t r;
    if((r = fwrite_LOOP(offset_vec, sizeof(size_t), numaligns, fp)) != numaligns){
      char *err = strerror(errno);
      printf("fwrite_LOOP of %llu size_t's (byte offsets) to %s failed (ret= %llu):ferror()=%d,feof()=%d\n%s\n",
	     (unsigned long long)numaligns, filename, (unsigned long long)r, ferror(fp),feof(fp),err);
      fflush(stdout);
      printf("\t offset_vec[0]=%lu, offset_vec[%lu]=%lu\n",offset_vec[0],numaligns-1,offset_vec[numaligns-1]);
      fflush(stdout);
      printf("\t Calling exit(1)\n");
      fflush(stdout); exit(1);
    }
  }

  if(VERB>=2){
    printf("Calling FILEclose(fpbin): filenamebin=%s\n",filenamebin);
    fflush(stdout);
  }
  FILEclose(fp);

  if(write_buffer){
    free(write_buffer);
    write_buffer = NULL;
  }
  if(write_bufferbin){
    free(write_bufferbin);
    write_bufferbin = NULL;
  }
  if(offset_vec){
    delete [] offset_vec;
    offset_vec = NULL;
  }

  if(1 /* NEW contig->HapSite[0] */){
    errno = 0;
    if(rename(tmpfilename,origfilename)){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to rename %s to %s:errno=%d:%s\n",tmpfilename,origfilename,eno,err);
      fflush(stdout);  exit(1);
    }
  }
  free(tmpfilename);
  free(origfilename);
}
