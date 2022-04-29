#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include <math.h>
#ifndef WIN32
#include <unistd.h>
#include <sys/stat.h>
#else
#include <direct.h>
#define getcwd _getcwd
#endif

#include "globals.h"
#include "parameters.h"
#include "Calign.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/input_svcheck.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

char *xmap_filename = 0, *query_filename = 0, *ref_filename = 0;
int StartRef = -1,NummapsRef = -1,StartQuery = -1,NummapsQuery = -1;
xmapEntry *SVtocheck = 0;/* The SVs to be checked : SVtocheck[i = 0 .. numrefmaps/2 - 1] corresponds to refmap[2*i .. 2i+1] */
int svMax = 0;
int svCnt = 0, svCntOrig = 0;


static int Id1Site1Id2Inc(xmapEntry *p1, xmapEntry *p2)
{
  Cmap *Ymap1 = p1->RefMap;
  Cmap *Ymap2 = p2->RefMap;
  Cmap *Xmap1 = p1->QryMap;
  Cmap *Xmap2 = p2->QryMap;
  double left1 = p1->refstartpos;
  double left2 = p2->refstartpos;
  double right1 = p1->refendpos;
  double right2 = p2->refendpos;

  return (Ymap1->id > Ymap2->id) ? 1 : (Ymap1->id < Ymap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 
    (right1 > right2) ? 1 : (right1 < right2) ? -1 : 
    (Xmap1->id > Xmap2->id) ? 1 : (Xmap1->id < Xmap2->id) ? -1 : 0;
}

static int Id1Id2Site1Inc(xmapEntry *p1, xmapEntry *p2)
{
  Cmap *Ymap1 = p1->RefMap;
  Cmap *Ymap2 = p2->RefMap;
  Cmap *Xmap1 = p1->QryMap;
  Cmap *Xmap2 = p2->QryMap;
  double left1 = p1->refstartpos;
  double left2 = p2->refstartpos;
  double right1 = p1->refendpos;
  double right2 = p2->refendpos;

  return (Ymap1->id > Ymap2->id) ? 1 : (Ymap1->id < Ymap2->id) ? -1 :
    (Xmap1->id > Xmap2->id) ? 1 : (Xmap1->id < Xmap2->id) ? -1 :
    (left1 > left2) ? 1 : (left1 < left2) ? -1 : 
    (right1 > right2) ? 1 : (right1 < right2) ? -1 : 0;
}

static int CmapIdInc(Cmap **p1, Cmap **p2)
{
  /* Cannot use p1[0]->id - p2[0]->id since that may be out of range for integer */
  return (p1[0]->id > p2[0]->id) ? 1 : (p1[0]->id < p2[0]->id) ? -1 : 0;
}

static int xmapEntryRefInc(xmapEntry *p1, xmapEntry *p2)
{
  if(p1->refcontigid > p2->refcontigid)
    return 1;
  if(p1->refcontigid < p2->refcontigid)
    return -1;

  if(p1->refstartpos > p2->refstartpos)
    return 1;
  if(p1->refstartpos < p2->refstartpos)
    return 0;

  return (p1->refendpos > p2->refendpos) ? 1 : (p1->refendpos < p2->refendpos) ? -1 : 0;
}

/* use binary search to locate fmap[i=0..num-1]->id == id */
Cmap *findidBS(Cmap **fmap, int num, long long id)
{
  int low = 0;
  int high = num-1;
  while(high > low){
    int mid = (low+high)/2;
    long long midID = fmap[mid]->id;
    if(midID < id)
      low = mid + 1;
    else if(midID > id)
      high = mid - 1;
    else
      low = high = mid;
  }
  if(high < low || fmap[low]->id != id)
    return NULL;
  return fmap[low];
}

/* Allocate and create a map in pmap[0] (with specified id) :
   If corrected==0 : Just copy a region of QryMap including the SV and +- svcheckKB around the SV
   If corrected==1 : Same, but replace the SV region of the QryMap with the corresponding region of RefMap (possibly flipped).

   NOTE : If pSV == 0 : just copy the entire QryMap
   NOTE : If svcheckKB is very large, the entire QryMap is copied (modified at the SV if corrected==1)

   pmap may be NULL (or will be first deallocated) and is fully allocated here.
 */
void mapsplice(Cmap * &pmap, long long id, xmapEntry *pSV, double svcheckKB, int corrected, Cmap *QryMap, Cmap *RefMap)
{
  if(DEBUG) assert(QryMap != NULL);
  if(DEBUG && corrected) assert(RefMap != NULL);

  if(pmap == NULL)
    pmap = new Cmap;
  if(id > MASK(31)){
    printf("mapsplice:id=%lld exceeds 31 bit limit\n",id);
    exit(1);
  }
  pmap->mapid = pmap->id = id;
  pmap->fileid = -1;
  for(int c = 0; c < colors; c++){
    pmap->Nickase[c] = strdup(QryMap->Nickase[c]);
    pmap->numsite[c] = 0;
    if(pmap->site[c]) {
      if(!pmap->blockmem) {
	delete [] pmap->site[c];
	delete [] pmap->siteSD[c];
	delete [] pmap->sitecov[c];
	delete [] pmap->sitecnt[c];
      }
      pmap->site[c] = NULL;
      pmap->siteSD[c] = NULL;
      pmap->sitecov[c] = NULL;
      pmap->sitecnt[c] = NULL;
    }
    if(pmap->SNRcnt[c]){
      delete [] pmap->SNRcnt[c];
      pmap->SNRcnt[c] = 0;
    }
  }
  pmap->startloc = pmap->endloc = 0.0;
  pmap->blockmem = 0;

  int *Nqry = QryMap->numsite;
  double qrylen = QryMap->site[0][Nqry[0]+1];

  double qryleft,qryright,qrystart,qryend;
  if(pSV){
    qryleft = min(pSV->qrystartpos,pSV->qryendpos);
    qryright = max(pSV->qrystartpos,pSV->qryendpos);
    qrystart = max(0.0, qryleft - svcheckKB);
    qryend = min(qrylen, qryright + svcheckKB);
    if(DEBUG) assert(qryleft < qryright - 1e-6);
  } else {
    qrystart = qryleft = qryright = 0.0;
    qryend = qrylen;
  }

  if(!corrected || !pSV){/* copy QryMap in Kb range qrystart to qryend */
    for(int c = 0; c < colors; c++){
      /* compute number of sites and allocate memory */
      int qstart = 1;
      while(qstart < Nqry[c] && QryMap->site[c][qstart] <= qrystart)
	qstart++;
      int qend = Nqry[c];
      while(qend >= qstart && qend > 1 && QryMap->site[c][qend] >= qryend)
	qend--;

      int NumSites = max(0,qend - qstart + 1);
      pmap->numsite[c] = NumSites;
      pmap->site[c] = new FLOAT[NumSites+2];
      pmap->siteSD[c] = new double[NumSites+2];
      pmap->sitecov[c] = new float[NumSites+2];
      pmap->sitecnt[c] = new float[NumSites+2];
      pmap->SNRcnt[c] = 0;
      if(QryMap->SNRcnt[c] || (corrected && RefMap->SNRcnt[c])){
	pmap->SNRcnt[c] = new int[NumSites+2];
	pmap->SNRgmean[c] = new double[NumSites+2];
	pmap->lnSNRsd[c] = new double[NumSites+2];
	pmap->SNRdist[c] = new double*[NumSites+2];
      }

      /* left end */
      pmap->site[c][0] = 0.0;
      pmap->siteSD[c][0] = 0.0;
      pmap->sitecov[c][0] = 1;
      pmap->sitecnt[c][0] = 1;
      if(pmap->SNRcnt[c]){
	pmap->SNRcnt[c][0] = 0;
	pmap->SNRgmean[c][0] = 0.0;
	pmap->lnSNRsd[c][0] = 0.0;
	pmap->SNRdist[c][0] = NULL;
      }

      /* right end */
      pmap->site[c][NumSites+1] = qryend - qrystart;
      pmap->siteSD[c][NumSites+1] = 0.0;
      pmap->sitecov[c][NumSites+1] = 1;
      pmap->sitecnt[c][NumSites+1] = 1;
      if(pmap->SNRcnt[c]){
	pmap->SNRcnt[c][NumSites+1] = 0;
	pmap->SNRgmean[c][NumSites+1] = 0.0;
	pmap->lnSNRsd[c][NumSites+1] = 0.0;
	pmap->SNRdist[c][NumSites+1] = NULL;
      }

      /* copy QryMap from sites qstart ... qend with left end at qrystart and right end at qryend */
      FLOAT *Y = QryMap->site[c];
      FLOAT *Z = pmap->site[c];
      int i = 1;
      for(int j = qstart;j <= qend; j++, i++){
	if(DEBUG) assert(i <= NumSites);
	Z[i] = Y[j] - qrystart;
	if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
	if(DEBUG && !(Z[i] < Z[NumSites+1])){
	  printf("mapsplice:id=%lld,pSV=%p,svcheckKB=%0.4f,corrected=%d,QryMap->id=%lld,RefMap->id=%lld,pSV=%p:qrystartpos=%0.6f,qryendpos=%0.6f\n",
		 id,pSV,svcheckKB,corrected,QryMap->id,RefMap->id,pSV,pSV->qrystartpos,pSV->qryendpos);
	  printf("  qrystart=%0.6f,qryend=%0.6f,qryleft=%0.6f,qryright=%0.6f,qrylen=%0.6f:c=%d,qstart=%d,qend=%d,NumSites=%d,Y[qstart]=%0.6f,Y[qend]=%0.6f\n",
		 qrystart,qryend,qryleft,qryright,qrylen, c, qstart, qend, NumSites, Y[qstart],Y[qend]);
	  printf("Z[NumSites+1]=%0.4f (should match qryend - qrystart), i=%d,j=%d:Y[j]=%0.6f, Z[i]=%0.6f\n",
		 Z[NumSites+1],i,j,Y[j],Z[i]);
	  fflush(stdout);
	  assert(Z[i] < Z[NumSites+1]);
	}
	pmap->siteSD[c][i] = QryMap->siteSD[c][j];
	pmap->sitecov[c][i] = QryMap->sitecov[c][j];
	pmap->sitecnt[c][i] = QryMap->sitecnt[c][j];
	if(pmap->SNRcnt[c]){
	  pmap->SNRdist[c][i] = NULL;
	  if(QryMap->SNRcnt[c]){
	    pmap->SNRcnt[c][i] = QryMap->SNRcnt[c][j];
	    pmap->SNRgmean[c][i] = QryMap->SNRgmean[c][j];
	    pmap->lnSNRsd[c][i] = QryMap->lnSNRsd[c][j];
	    if(pmap->SNRcnt[c][i]){
	      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
	      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
		pmap->SNRdist[c][i][t] = QryMap->SNRdist[c][j][t];
	    }
	  } else {/* use default filler values */
	    pmap->SNRcnt[c][i] = 0;
	    pmap->SNRgmean[c][i] = Z[i];
	    pmap->lnSNRsd[c][i] = 0.0;
	  }
	}
      }
    }
  } else {/* copy QryMap from qrystart to qryleft (inclusive), then the RefMap SV region (possibly flipped), then QryMap from qryright(inclusive) to qryend */
    double refleft = pSV->refstartpos;
    double refright = pSV->refendpos;
    int *Nref = RefMap->numsite;
    for(int c = 0; c < colors; c++){
      /* compute total number of sites and allocate memory */

      /* locate first site beyond qrystart */
      int qstart = 1;
      while(qstart < Nqry[c] && QryMap->site[c][qstart] <= qrystart + 1e-6)
	qstart++;

      /* locate first site at OR before qryleft */
      int qleft = qstart - 1;
      while(qleft < Nqry[c] && QryMap->site[c][qleft+1] <= qryleft + 1e-6)
	qleft++;
      /*      if(DEBUG && colors==1 && !(fabs(QryMap->site[c][qleft] - qryleft) < 1e-6)){
	printf("qryid=%d,refid=%d:qrystart=%0.3f,qryleft=%0.3f,qryright=%0.3f,qryend=%0.3f:qleft=%d(%0.3f)\n",
	       pSV->qrycontigid,pSV->refcontigid,qrystart,qryleft,qryright,qryend,qleft,QryMap->site[c][qleft]);
	printf("Unable To find site on QryMap corresponding to qryleft:QryMap->site[%d][%d]=%0.3f,QryMap->site[%d][%d]=%0.3f\n",
	       c,qleft-1,QryMap->site[0][qleft-1],c,qleft+1,QryMap->site[0][qleft+1]);
	fflush(stdout);
	assert(fabs(QryMap->site[c][qleft] - qryleft) < 1e-6);
	}*/

      /* locate first site at OR beyond qryright */
      int qright = qleft;
      while(qright <= Nqry[c] && QryMap->site[c][qright] < qryright - 1e-6)
	qright++;
      if(DEBUG && !(qright > qleft)){
	printf("qryid=%lld,refid=%lld:qrystart=%0.3f,qryleft=%0.3f,qryright=%0.3f,qryend=%0.3f:c=%d:qleft=%d(%0.3f),qright=%d(%0.3f),M[c]=%d\n",
	       pSV->qrycontigid,pSV->refcontigid,qrystart,qryleft,qryright,qryend,c,qleft,QryMap->site[c][qleft],qright,QryMap->site[c][qright],Nqry[c]);
	fflush(stdout);
	assert(qright > qleft);
      }
      /*      if(DEBUG && colors==1) assert(fabs(QryMap->site[c][qright] - qryright) < 1e-6);*/

      /* locate last site before qryend */
      int qend = Nqry[c];
      while(qend >= qstart && qend > 1 && QryMap->site[c][qend] >= qryend - 1e-6)
	qend--;

      /* Locate first site at OR beyond refleft */
      int rleft = 1;
      while(rleft < Nref[c] && RefMap->site[c][rleft] < refleft - 1e-6)
	rleft++;

      /* Locate last site at OR before refright */
      int rright = Nref[c];
      while(rright >= 1  && RefMap->site[c][rright] > refright + 1e-6)
	rright--;

      if(pSV->orientforward){
	/* Do NOT include both sites at qryleft AND at refleft : discard the site at refleft, if needed */
	if(fabs(QryMap->site[c][qleft] - qryleft) < 1e-6 && fabs(RefMap->site[c][rleft] - refleft) < 1e-6)
	  rleft++;
	/* Do NOT include both sites at qryright AND at refright : discard the site at refright, if needed */
	if(fabs(QryMap->site[c][qright] - qryright) < 1e-6 && fabs(RefMap->site[c][rright] - refright) < 1e-6)
	  rright--;
      } else {
	/* Do NOT include both sites at qryright AND at refleft : discard the site at refleft, if needed */
	if(fabs(QryMap->site[c][qright] - qryright) < 1e-6 && fabs(RefMap->site[c][rleft] - refleft) < 1e-6)
	  rleft++;
	/* Do NOT include both sites at qryleft AND at refright : discard the site at refright, if needed */
	if(fabs(QryMap->site[c][qleft] - qryleft) < 1e-6 && fabs(RefMap->site[c][rright] - refright) < 1e-6)
	  rright--;
      }

      /* NOTE qryleft may be < qrystart (due to misresolved labels on left). Similarly qryright may be > qryend */
      /* NOTE rright may be < rleft if the SV interval has no internal sites on Reference AND query interval coincides with sites in query */

      int NumSites = max(0,qleft-qstart+1) + max(0,rright - rleft + 1) + max(0,qend-qright+1);
      pmap->numsite[c] = NumSites;
      pmap->site[c] = new FLOAT[NumSites+2];
      pmap->siteSD[c] = new double[NumSites+2];
      pmap->sitecov[c] = new float[NumSites+2];
      pmap->sitecnt[c] = new float[NumSites+2];
      pmap->SNRcnt[c] = 0;
      if(QryMap->SNRcnt[c] || RefMap->SNRcnt[c]){
	pmap->SNRcnt[c] = new int[NumSites+2];
	pmap->SNRgmean[c] = new double[NumSites+2];
	pmap->lnSNRsd[c] = new double[NumSites+2];
	pmap->SNRdist[c] = new double*[NumSites+2];
      }

      /* left end */
      pmap->site[c][0] = 0.0;
      pmap->siteSD[c][0] = 0.0;
      pmap->sitecov[c][0] = 1;
      pmap->sitecnt[c][0] = 1;
      if(pmap->SNRcnt[c]){
	pmap->SNRcnt[c][0] = 0;
	pmap->SNRgmean[c][0] = 0.0;
	pmap->lnSNRsd[c][0] = 0.0;
	pmap->SNRdist[c][0] = NULL;
      }

      /* right end */
      pmap->site[c][NumSites+1] = max(MININTERVAL, qryleft - qrystart) + (refright - refleft) + max(MININTERVAL, qryend - qryright);
      pmap->siteSD[c][NumSites+1] = 0.0;
      pmap->sitecov[c][NumSites+1] = 1;
      pmap->sitecnt[c][NumSites+1] = 1;
      if(pmap->SNRcnt[c]){
	pmap->SNRcnt[c][NumSites+1] = 0;
	pmap->SNRgmean[c][NumSites+1] = 0.0;
	pmap->lnSNRsd[c][NumSites+1] = 0.0;
	pmap->SNRdist[c][NumSites+1] = NULL;
      }

      /* copy QryMap including sites qstart to qleft and locations qrystart to qryleft (If qrystart < qryleft) */
      if(DEBUG && qrystart >= qryleft) assert(qstart > qleft);

      FLOAT *Y = QryMap->site[c];
      FLOAT *Z = pmap->site[c];
      int i = 1;
      for(int j = qstart;j <= qleft; j++, i++){
	if(DEBUG) assert(i <= NumSites);
	Z[i] = Y[j] - qrystart;
	if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
	if(DEBUG) assert(Z[i] < Z[NumSites+1]);
	pmap->siteSD[c][i] = QryMap->siteSD[c][j];
	pmap->sitecov[c][i] = QryMap->sitecov[c][j];
	pmap->sitecnt[c][i] = QryMap->sitecnt[c][j];
	if(pmap->SNRcnt[c]){
	  pmap->SNRdist[c][i] = NULL;
	  if(QryMap->SNRcnt[c]){
	    pmap->SNRcnt[c][i] = QryMap->SNRcnt[c][j];
	    pmap->SNRgmean[c][i] = QryMap->SNRgmean[c][j];
	    pmap->lnSNRsd[c][i] = QryMap->lnSNRsd[c][j];
	    if(pmap->SNRcnt[c][i]){
	      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
	      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
		pmap->SNRdist[c][i][t] = QryMap->SNRdist[c][j][t];
	    }
	  } else {/* use default filler values */
	    pmap->SNRcnt[c][i] = 0;
	    pmap->SNRgmean[c][i] = Z[i];
	    pmap->lnSNRsd[c][i] = 0.0;
	  }
	}
      }
      
      FLOAT *X = RefMap->site[c];
      if(pSV->orientforward){ /* copy RefMap from sites rleft to rright and locations refleft to refright */
	for(int k = rleft; k <= rright; k++, i++){
	  if(DEBUG) assert(i <= NumSites);
	  if(DEBUG) assert(refright - X[k] < refright - refleft + 1e-6);
	  Z[i] = max(MININTERVAL, qryleft-qrystart) + X[k] - refleft;
	  if(DEBUG && i > 0 && !(Z[i] > Z[i-1])){
	    printf("qryid=%lld,refid=%lld,or=0:qrystart=%0.3f,qryleft=%0.3f,qryright=%0.3f,qryend=%0.3f:c=%d:qstart=%d(%0.3f),qleft=%d(%0.3f),qright=%d(%0.3f),qend=%d(%0.3f),M[c]=%d\n",
		   pSV->qrycontigid,pSV->refcontigid,qrystart,qryleft,qryright,qryend,c,qstart,QryMap->site[c][qstart],qleft,QryMap->site[c][qleft],qright,QryMap->site[c][qright],qend,QryMap->site[c][qend],Nqry[c]);
	    printf("     refleft=%0.3f,refright=%0.3f,refend=%0.3f:rleft=%d(%0.3f),rright=%d(%0.3f),N[c]=%d\n",
		   refleft,refright,X[Nref[c]+1],rleft,X[rleft],rright,X[rright],Nref[c]);
	    printf("i=%d,k=%d:Z[i-1]=%0.6f, X[k]=%0.6f,Z[i]=%0.6f\n",i,k,Z[i-1],X[k],Z[i]);
	    fflush(stdout);
	    assert(Z[i] > Z[i-1]);
	  }
	  if(DEBUG) assert(Z[i] < Z[NumSites+1]);
	  pmap->siteSD[c][i] = RefMap->siteSD[c][k];
	  pmap->sitecov[c][i] = RefMap->sitecov[c][k];
	  pmap->sitecnt[c][i] = RefMap->sitecnt[c][k];
	  if(pmap->SNRcnt[c]){
	    pmap->SNRdist[c][i] = NULL;
	    if(RefMap->SNRcnt[c]){
	      pmap->SNRcnt[c][i] = RefMap->SNRcnt[c][k];
	      pmap->SNRgmean[c][i] = RefMap->SNRgmean[c][k];
	      pmap->lnSNRsd[c][i] = RefMap->lnSNRsd[c][k];
	      if(pmap->SNRcnt[c][i]){
		pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
		  pmap->SNRdist[c][i][t] = RefMap->SNRdist[c][k][t];
	      }
	    } else {/* use default filler values */
	      pmap->SNRcnt[c][i] = 0;
	      pmap->SNRgmean[c][i] = Z[i];
	      pmap->lnSNRsd[c][i] = 0.0;
	    }
	  }
	}
      } else {/* copy RefMap from sites rright to rleft and locations refright to refleft (reversed orientatation) */
	for(int k = rright; k >= rleft; k--, i++){
	  if(DEBUG) assert(i <= NumSites);
	  if(DEBUG) assert(refright - X[k] < refright - refleft + 1e-6);
	  Z[i] = max(MININTERVAL, qryleft-qrystart) + refright - X[k];
	  if(DEBUG && i > 0) assert(Z[i] > Z[i-1]);
	  if(DEBUG) assert(Z[i] < Z[NumSites+1]);
	  pmap->siteSD[c][i] = RefMap->siteSD[c][k];
	  pmap->sitecov[c][i] = RefMap->sitecov[c][k];
	  pmap->sitecnt[c][i] = RefMap->sitecnt[c][k];
	  if(pmap->SNRcnt[c]){
	    pmap->SNRdist[c][i] = NULL;
	    if(RefMap->SNRcnt[c]){
	      pmap->SNRcnt[c][i] = RefMap->SNRcnt[c][k];
	      pmap->SNRgmean[c][i] = RefMap->SNRgmean[c][k];
	      pmap->lnSNRsd[c][i] = RefMap->lnSNRsd[c][k];
	      if(pmap->SNRcnt[c][i]){
		pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
		for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
		  pmap->SNRdist[c][i][t] = RefMap->SNRdist[c][k][t];
	      }
	    } else {/* use default filler values */
	      pmap->SNRcnt[c][i] = 0;
	      pmap->SNRgmean[c][i] = Z[i];
	      pmap->lnSNRsd[c][i] = 0.0;
	    }
	  }
	}
      }

      if(DEBUG) assert(refright > refleft);

      /* copy QryMap from sites qright to qend and locations qryright to qryend */
      if(DEBUG && qryright >= qryend) assert(qright > qend);
      for(int j = qright;j <= qend; j++, i++){
	if(DEBUG) assert(i <= NumSites);
	if(DEBUG && !(Y[j] >= qryright - 1e-6)){
	  printf("QryMap->id=%lld,M=%d,RefMap->id=%lld,N=%d,pSV=%p:qrystartpos=%0.6f,qryendpos=%0.6f\n",
		 QryMap->id,QryMap->numsite[0],RefMap->id,RefMap->numsite[0], pSV, pSV->qrystartpos, pSV->qryendpos);
	  printf("j=%d:qrystart=%0.6f,qryleft=%0.6f,qryright=%0.6f,qryend=%0.6f:qstart=%d(%0.6f),qleft=%d(%0.6f),qright=%d(%0.6f),qend=%d(%0.6f):Y[j]=%0.6f\n",
		 j,qrystart,qryleft,qryright,qryend,qstart,Y[qstart],qleft,Y[qleft],qright,Y[qright],qend,Y[qend],Y[j]);
	  printf("  refleft=%0.3f,refright=%0.3f,rleft=%d(%0.3f),rright=%d(%0.3f)\n",refleft,refright,rleft,X[rleft],rright,X[rright]);
	  printf("  i=%d:Z[i-1]=%0.3f,orientation=%d\n",i,Z[i-1],pSV->orientforward);
	  for(int t = qright; t <= qend; t++)
	    printf("t=%d,Y[t]=%0.6f\n",t,Y[t]);
	  fflush(stdout);
	  assert(Y[j] >= qryright - 1e-6);
	}
	Z[i] = max(MININTERVAL, qryleft - qrystart) + (refright - refleft) + Y[j] - qryright;
	if(DEBUG && i > 0 && !(Z[i] > Z[i-1])){
	  printf("j=%d:qrystart=%0.3f,qryleft=%0.3f,qryright=%0.3f,qryend=%0.3f:qstart=%d(%0.3f),qleft=%d(%0.3f),qright=%d(%0.3f),qend=%d(%0.3f):Y[j]=%0.3f\n",
		 j,qrystart,qryleft,qryright,qryend,qstart,Y[qstart],qleft,Y[qleft],qright,Y[qright],qend,Y[qend],Y[j]);
	  printf("  refleft=%0.3f,refright=%0.3f,rleft=%d(%0.3f),rright=%d(%0.3f)\n",refleft,refright,rleft,X[rleft],rright,X[rright]);
	  printf("  i=%d:Z[i]=%0.3f,Z[i-1]=%0.3f,orientation=%d\n",i,Z[i],Z[i-1],pSV->orientforward);
	  fflush(stdout);
	  assert(Z[i] > Z[i-1]);
	}
	if(DEBUG) assert(Z[i] < Z[NumSites+1]);
	pmap->siteSD[c][i] = QryMap->siteSD[c][j];
	pmap->sitecov[c][i] = QryMap->sitecov[c][j];
	pmap->sitecnt[c][i] = QryMap->sitecnt[c][j];
	if(pmap->SNRcnt[c]){
	  pmap->SNRdist[c][i] = NULL;
	  if(QryMap->SNRcnt[c]){
	    pmap->SNRcnt[c][i] = QryMap->SNRcnt[c][j];
	    pmap->SNRgmean[c][i] = QryMap->SNRgmean[c][j];
	    pmap->lnSNRsd[c][i] = QryMap->lnSNRsd[c][j];
	    if(pmap->SNRcnt[c][i]){
	      pmap->SNRdist[c][i] = new double[pmap->SNRcnt[c][i]];
	      for(int t = 0; t < pmap->SNRcnt[c][i]; t++)
		pmap->SNRdist[c][i][t] = QryMap->SNRdist[c][j][t];
	    }
	  } else {/* use default filler values */
	    pmap->SNRcnt[c][i] = 0;
	    pmap->SNRgmean[c][i] = Z[i];
	    pmap->lnSNRsd[c][i] = 0.0;
	  }
	}
      }
    }
  }
}

static char buf_input_indel[LINESIZ];

/* read in the SVs and append query & ref maps to map[0..nummaps-1] */
void input_indel(char *filename, xmapEntry * &SVtocheck, int &svCnt, int &svMax, int skipTwoMatchgroup = 0)
{
  char *buf = buf_input_indel;
  char *pt;
  char *qt;
  FILE *fp;

  if(rres > 0.0){
    printf("WARNING: input_indel: rres = %0.3f ignored, since reference map in _r.cmap should already be de-resed\n",rres);
    fflush(stdout);

    rres = 0;
  }

  int SMAP_SUFFIX = 0;
  if((pt = strstr(filename,".smap")) && !strcmp(pt,".smap"))
    SMAP_SUFFIX = 1;

  int SMAP = 0;/* SMAP revision */
  int smapCnt = 0;

  if((fp = fopen(filename,"r"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("Failed to read input file %s:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  int SV2cnt = 0;/* number 2-matchgroup SVs skipped from .smap */

  xmapEntry *xmap = 0;
  int xmapCnt = 0, xmapMax = 0;

  int linecnt = 1;
  for(;fgets(buf,LINESIZ,fp) != NULL; linecnt++){
    int len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      exit(1);
    }
    
    if(buf[0] == '#'){/* most comment lines are ignored, except "SMAP File Version", "Reference Maps From:" , "Query Maps From:" and "Label Channels:" */
      char *key = (char *)"SMAP File Version:";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);
	double version = strtod(pt,&qt);
	if(qt == pt || !(*qt == 0 || isspace(*qt)) || version < 0.0){
	  printf("Line %d in %s with %s has invalid version number (must be non-negative):\n%s\n",linecnt,filename,key,buf);
	  exit(1);
	}
	if(SMAP_SUFFIX){
	  if(fabs(version - 0.5) <= 1e-6)
	    SMAP = 5;
	  else if(fabs(version - 0.6) <= 1e-6)
	    SMAP = 6;
	  else if(fabs(version - 0.7) <= 1e-6)
	    SMAP = 7;
	  else {
	    printf("WARNING: Line %d in %s with %s has unsupported version (must be 0.5, 0.6 or 0.7):\n%s\n",linecnt,filename,key,buf);
	    fflush(stdout);
	  }
	} else {
	  if(version != 0.0){
	    printf("Line %d in %s with %s has unsupported version (must be 0.0):\n%s\n",linecnt,filename,key,buf);
	    exit(1);
	  }
	  SMAP = 0;
	}
	continue;
      }

      key = (char *)"Label Channels:";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);
	int NColors = strtol(pt,&qt,10);
	if(qt == pt || !(*qt == 0 || isspace(*qt)) || NColors <= 0 || NColors > 2){
	  printf("Line %d in %s with %s has invalid integer value (must be 1 or 2):\n%s\n",linecnt, filename, key, buf);
	  exit(1);
	}
	if(NColors != colors){
	  printf("Label Channels on Line %d in %s does not match -colors %d\n",linecnt, filename, colors);
	  exit(1);
	}
	continue;
      }

      key = (char *)"Reference Maps From:";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);
	char name[PATH_MAX];
	int i;
	errno = 0;
	if(1 != (i=sscanf(pt,"%s",name))){
	  char *err  = strerror(errno);
	  printf("sscanf returned %d (expected 1) on line %d of %s:%s\n",i,linecnt,filename,err);
	  exit(1);
	}
	ref_filename = strdup(name);
	continue;
      }

      key = (char *)"Query Maps From:";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);
	char name[PATH_MAX];
	int i;
	errno = 0;
	if(1 != (i=sscanf(pt,"%s",name))){
	  char *err  = strerror(errno);
	  printf("sscanf returned %d (expected 1) on line %d of %s:%s\n",i,linecnt,filename,err);
	  exit(1);
	}
	query_filename = strdup(name);
	continue;
      }
      
      key = (char *)"Xmap Entries From:";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);
	char name[PATH_MAX];
	int i;
	errno = 0;
	if(1 != (i=sscanf(pt,"%s",name))){
	  char *err  = strerror(errno);
	  printf("sscanf returned %d (expected 1) on line %d of %s:%s\n",i,linecnt,filename,err);
	  exit(1);
	}
	xmap_filename = strdup(name);
	continue;
      }

      continue;/* skip all other comment lines */
    }

    /* done with header comment lines */
    if(!ref_filename){
      printf("Reference Maps name not found in header lines of %s\n",filename);
      exit(1);
    }
    if(!query_filename){
      printf("Query Maps name not found in header lines of %s\n",filename);
      exit(1);
    }

    if(StartRef < 0){
      /* disable all filtering or rescaling of input maps */
      double origPixelLen = PixelLen;
      double origmres = mres;
      double origMinLen = MinLen;
      double origMaxLen = MaxLen;
      int origMinSites = MinSites;
      int origMaxSites = MaxSites;
      int origResBins[MAXCOLOR] = {0,0};
      double origminSNR[MAXCOLOR] = {0,0};

      /* suppress checking for QX codes, since input maps are typically CMAPS */
      int origMapTrueSite = MapTrueSite;
      int origMapSNR = MapSNR;
      int origMapIntensity = MapIntensity;
      int origMapStitched = MapStitched;
      int origMapPSFWidth = MapPSFWidth;
      int origMapImageLoc = MapImageLoc;

      PixelLen = 0.500;/* Read in Reference Maps and Query maps without rescaling */
      mres = 0.0;
      MinLen = MaxLen = 0.0;
      MaxSites = MinSites = 0;
      for(int c = 0; c < colors; c++){
	origResBins[c] = ResBins[c];
	origminSNR[c] = minSNR[c];

	ResBins[c] = 0;/* Don't apply bias correction to Reference Maps and Query maps */
	minSNR[c] = 0.0;
      }
      MapTrueSite = MapSNR = MapIntensity = MapStitched = MapPSFWidth = MapImageLoc = 0;

#ifndef WIN32
      /* extract .indel directory as absolute pathname from .indel filename */
      char basename[PATH_MAX], *qt;
      if(filename[0] == '/')/* absolute path name */
	strcpy(basename,filename);
      else {
	if(getcwd(basename,PATH_MAX) == NULL){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("getcwd() failed (errno=%d):%s\n",eno, err);
	  exit(1);
	}
	pt = basename + strlen(basename);
	sprintf(pt,"/%s",filename);
      }

      if((pt = strrchr(basename,'/')))
	pt++;
      else {
	printf("Could not find '/' in absolute .indel filename:%s\n",basename);
	exit(1);
      }

      struct stat statbuf;

      if(xmap_filename[0] != '/' || stat(xmap_filename, &statbuf)) { /* change xmap_filename into absolute pathname */
	if((qt = strrchr(xmap_filename,'/')))
	  strcpy(pt,qt+1);
	else
	  strcpy(pt,xmap_filename);
	if(VERB/* HERE >=2 */ && strcmp(xmap_filename,basename)){
	  printf("Changing xmap_filename from \"%s\" to \"%s\"\n",xmap_filename,basename);
	  fflush(stdout);
	}
	free(xmap_filename);
	xmap_filename = strdup(basename);
      }
      if(ref_filename[0] != '/' || stat(ref_filename, &statbuf)){/* change ref_filename into absolute pathname */
	if((qt = strrchr(ref_filename,'/')))
	  strcpy(pt,qt+1);
	else
	  strcpy(pt, ref_filename);
	if(VERB/* HERE >=2 */ && strcmp(ref_filename,basename)){
	  printf("Changing ref_filename from \"%s\" to \"%s\"\n",ref_filename,basename);
	  fflush(stdout);
	}
	free(ref_filename);
	ref_filename = strdup(basename);
      }
      if(query_filename[0] != '/' || stat(query_filename, &statbuf)){/* change query_filename into abolute pathname */
	if((qt = strrchr(query_filename,'/')))
	  strcpy(pt,qt+1);
	else
	  strcpy(pt, query_filename);
	if(VERB/* HERE >=2 */ && strcmp(query_filename,basename)){
	  printf("Changing query_filename from \"%s\" to \"%s\"\n", query_filename,basename);
	  fflush(stdout);
	}
	free(query_filename);
	query_filename = strdup(basename);
      }
#endif

      int origmaptype = maptype;

      /* read in Reference Maps to map[nummaps ...] (Process like with -ref) */
      if(num_files+1 >= MAXFILES){
	printf("Too many input files : increase MAXFILES=%d to at least %d\n", MAXFILES, num_files+2);
	exit(1);
      }
      vfixx_filename[num_files] = strdup(ref_filename);
      if(VERB>=2){
	printf("Reading Reference maps from %s(ref_filename=%s)\n",vfixx_filename[num_files],ref_filename);
	fflush(stdout);
      }

      StartRef = nummaps;
      NummapsRef = input_vfixx(num_files, nummaps, maxmaps, Gmap);
      if(DEBUG) assert(nummaps == StartRef + NummapsRef);
      num_files++;

      qsort(&Gmap[StartRef],NummapsRef,sizeof(Cmap *),(intcmp*)CmapIdInc);      

      if(VERB/* HERE >=2 */){
	printf("Read %d Reference maps from %s\n",NummapsRef,vfixx_filename[num_files-1]);
	fflush(stdout);
      }

      /* append Query Maps to map[nummaps ... ] */
      if(++num_files >= MAXFILES){
	printf("Too many input files : increase MAXFILES=%d to at least %d\n", MAXFILES, num_files+1);
	exit(1);
      }
      vfixx_filename[num_files-1] = strdup(query_filename);
      if(VERB>=2){
	printf("Reading Query maps from %s(query_filename=%s)\n",vfixx_filename[num_files-1],query_filename);
	fflush(stdout);
      }
      StartQuery = nummaps;
      NummapsQuery = input_vfixx(num_files-1,nummaps,maxmaps, Gmap);
      if(DEBUG) assert(nummaps == StartQuery + NummapsQuery);
      qsort(&Gmap[StartQuery],NummapsQuery,sizeof(Cmap *),(intcmp*)CmapIdInc);      

      if(VERB/* HERE >=2 */){
	printf("Read %d Query maps from %s\n",NummapsQuery,vfixx_filename[num_files-1]);
	if(VERB>=3){
	  long long qryid = 5LL;
	  Cmap *QryMap = findidBS(&Gmap[StartQuery], NummapsQuery,qryid);
	  printf("qryid=%lld:numsite=%d\n",qryid, QryMap->numsite[0]);
	  for(int i = 1; i <= QryMap->numsite[0]+1; i++)
	    printf("QryMap->site[0][%d] = %0.3f\n",i,QryMap->site[0][i]);
	}
	fflush(stdout);
      }
      
      if(SMAP){/* input the xmaps, since we need the orientation of matchgroups to validate querystart,queryend values */
	input_xmap(xmap_filename,xmap,xmapCnt,xmapMax, StartQuery, NummapsQuery, StartRef, NummapsRef, query_filename, ref_filename);
	if(DEBUG/* HERE >=2 */){
	  for(int i = 0; i < xmapCnt; i++){
	    xmapEntry *pXmap = &xmap[i];
	    assert(pXmap->smapentryid == i+1);
	  }
	}
      }

      maptype = origmaptype;
      
      for(int c = 0; c < colors; c++){
	ResBins[c] = origResBins[c];
	minSNR[c]= origminSNR[c];
      }
      PixelLen = origPixelLen;
      mres = origmres;
      MinLen = origMinLen;
      MaxLen = origMaxLen;
      MinSites = origMinSites;
      MaxSites = origMaxSites;

      MapTrueSite = origMapTrueSite;
      MapSNR = origMapSNR;
      MapIntensity = origMapIntensity;
      MapStitched = origMapStitched;
      MapPSFWidth = origMapPSFWidth;
      MapImageLoc = origMapImageLoc;
    }

    if(svCnt >= svMax){
      svMax = max(2*svMax,1024);
      if(DEBUG) assert(svMax > 0);// check for 31 bit overflow
      xmapEntry *newSVs = new xmapEntry[svMax];
      if(svCnt)
	memcpy(newSVs,SVtocheck,sizeof(xmapEntry)*svCnt);
      delete [] SVtocheck;
      SVtocheck = newSVs;
    }
    xmapEntry *pSV = &SVtocheck[svCnt];
    pSV->mergecnt = 1;

    /* parse .smap lines with insertion or deletion */
    pt = buf;
    int SmapEntryId = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || SmapEntryId <= 0){
      printf("Invalid SmapEntryId on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pSV->smapentryid = SmapEntryId;
    if(DEBUG && SmapEntryId != (SMAP ? smapCnt : svCnt)+1){
      printf("Invalid SmapEntryId=%d on line %d of %s (should match %d'th SV line)\n%s\n",SmapEntryId,linecnt,filename,(SMAP ? smapCnt : svCnt)+1,buf);
      printf("\t pt=%s\n",pt);
      printf("\t qt=%s\n",qt);
      exit(1);
    }
    smapCnt++;

    pt = qt;
    pSV->qrycontigid = strtol(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pSV->qrycontigid <= 0){
      printf("Invalid QryContigID on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    Cmap *QryMap = findidBS(&Gmap[StartQuery],NummapsQuery,(long long)pSV->qrycontigid);
    if(!QryMap){
      printf("QryContigID=%lld on line %d of %s not found amongst %d Query maps in %s\n",
	     pSV->qrycontigid, linecnt, filename, NummapsQuery, vfixx_filename[num_files-1]);
      exit(1);
    }
    
    if(DEBUG) assert(QryMap != NULL && QryMap->id == pSV->qrycontigid);
    int *Nqry = QryMap->numsite;
    
    pt = qt;
    pSV->refcontigid = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pSV->refcontigid <= 0){
      printf("Invalid RefContigID1 on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    Cmap *RefMap = findidBS(&Gmap[StartRef],NummapsRef,(long long)pSV->refcontigid);
    if(!RefMap){
      printf("RefContigID1=%lld on line %d of %s not found amongst %d Reference maps in %s\n",
	     pSV->refcontigid, linecnt, filename, NummapsRef, vfixx_filename[num_files-2]);
      exit(1);
    }
    if(DEBUG) assert(RefMap != NULL && RefMap->id == pSV->refcontigid);
    int *NRef = RefMap->numsite;

    pt = qt;
    long long RefContigID2 = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt))){
      printf("Invalid RefContigID2 on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    /* check later that pSV->refcontigid == RefContigID2 for insertions or deletions */
    
    pt = qt;
    pSV->qrystartpos = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt))/* || pSV->qrystartpos < (SMAP ? -1.0001 : -10000.0)*/){
      printf("Invalid QryStartPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    
    if(pSV->qrystartpos < -1.0001)
      pSV->qrystartpos = 0.001;

    if(!SMAP && pSV->qrystartpos > 0.0){
      pSV->qrystartpos *= 0.001;
      if(pSV->qrystartpos > QryMap->site[0][Nqry[0]+1] + 10.0){
	printf("QryStartPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s (+10kb): pSV->qrycontigid=%lld\n",
	       pSV->qrystartpos,linecnt,filename,QryMap->site[0][Nqry[0]+1], QryMap->id, vfixx_filename[num_files-1], pSV->qrycontigid);
	exit(1);
      }
    }

    pt = qt;
    pSV->qryendpos = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pSV->qryendpos < (SMAP ? -1.0001 : -10000.0)){
      printf("Invalid QryEndPos= %0.2f on line %d of %s:\n%s\n",pSV->qryendpos, linecnt, filename, buf);
      exit(1);
    }
    
    if(!SMAP && pSV->qryendpos > 0.0){
      pSV->qryendpos *= 0.001;
      if(pSV->qryendpos > QryMap->site[0][Nqry[0]+1] + 10.0){
	printf("QryEndPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s (+10kb)\n",
	       pSV->qryendpos,linecnt,filename,QryMap->site[0][Nqry[0]+1], QryMap->id, vfixx_filename[num_files-1]);
	exit(1);
      }
    }

    pt = qt;
    pSV->refstartpos = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pSV->refstartpos < (SMAP ? -1.000001 : 0.0)){
      printf("Invalid RefStartPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
#if 0 // SMAP case is handled later after verifying that this is a deletion or insertion and NOT a dual matchgroup case */
    if(pSV->refstartpos > 0.0){
      pSV->refstartpos *= 0.001;
      if(pSV->refstartpos > RefMap->site[0][NRef[0]+1]){
	printf("RefStartPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s\n",
	       pSV->refstartpos,linecnt,filename,RefMap->site[0][NRef[0]+1], RefMap->id, vfixx_filename[num_files-2]);
	exit(1);
      }
    }
#endif

    pt = qt;
    pSV->refendpos = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pSV->refendpos < (SMAP ? -1.00001 : 0.0)){
      printf("Invalid RefEndPos= %0.2f on line %d of %s:\n%s\n",pSV->refendpos, linecnt,filename,buf);
      exit(1);
    }
#if 0 // SMAP case is handled later after verifying that this is a deletion or insertion */
    if(pSV->refendpos > 0.0){
      pSV->refendpos *= 0.001;
      if(pSV->refendpos > RefMap->site[0][NRef[0]+1]){
	printf("RefEndPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s\n",
	       pSV->refendpos,linecnt,filename,RefMap->site[0][NRef[0]+1], RefMap->id, vfixx_filename[num_files-2]);
	exit(1);
      }
    }
#endif
    while(*qt && isspace(*qt))
      qt++;

    if(SMAP==0){  /* read orientation character */
      if(*qt != '+' && *qt != '-'){
	printf("Invalid Orientation=%c on line %d of %s:\n%s\n",*qt,linecnt,filename,buf);
	exit(1);
      }
      pSV->orientforward = (*qt == '+');
      while(*qt && !isspace(*qt))
	qt++;

      while(*qt && isspace(*qt))
	qt++;
    }
    
    pt = qt;
    pSV->confidence = strtod(pt,&qt);// NOTE : if SMAP==6, then this is PPV
    if(pt == qt || !(*qt==0 || isspace(*qt)) || !isfinite(pSV->confidence)){
      printf("Invalid Confidence %0.3f on line %d of %s:%s\n%s\n",pSV->confidence,linecnt,filename,pt,buf);
      exit(1);
    }

    /* check Type field : if it isn't insertion or deletion skip this SV */
    while(*qt && isspace(*qt))
      qt++;
    pt = qt;
    if(strncmp(qt,"insertion",strlen("insertion")) && strncmp(pt,"deletion",strlen("deletion")))
      continue;/* skip this SV */
    while(*qt && !isspace(*qt))
      qt++;
    char c = *qt;
    *qt = '\0';
    pSV->type = strdup(pt);
    *qt = c;

    if(pSV->refcontigid != RefContigID2){
      printf("RefcontigID1=%lld does not match RefcontigID2=%lld on line %d of %s\n",pSV->refcontigid,RefContigID2,linecnt,filename);
      fflush(stdout);exit(1);
    }

    pt = qt;
    pSV->xmapid = strtol(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pSV->xmapid <= 0){
      printf("Invalid XmapID1 on line %d of %s:\n%s\n",linecnt,filename,buf);
      fflush(stdout);exit(1);
    }

    pt = qt;
    int XmapID2 = strtol(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || XmapID2 <= 0){
      printf("Invalid XmapID2 on line %d of %s:\n%s\n",linecnt,filename,buf);
      fflush(stdout);exit(1);
    }
    if(SMAP==0 && pSV->xmapid != XmapID2){/* with .smap, the two values need not match for large indels */
      printf("XmapID1=%d and XmapID2=%d on line %d of %s do not match:\n%s\n",pSV->xmapid,XmapID2,linecnt,filename,buf);
      fflush(stdout);exit(1);
    }
    if(SMAP && skipTwoMatchgroup && pSV->xmapid != XmapID2){/* for now, just skip 2 matchgroup calls */
      if(VERB>=2){
	printf("Skipping Two matchgroup SV: XmapID1=%d and XmapID2=%d on line %d of %s\n",pSV->xmapid,XmapID2,linecnt,filename);
	fflush(stdout);
      }
      SV2cnt++;
      continue;
    }
    
    if(SMAP && pSV->qrystartpos > 0.0){
      pSV->qrystartpos *= 0.001;
      if(pSV->qrystartpos > QryMap->site[0][Nqry[0]+1] + 10.0){
	printf("QryStartPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s (+10kb): pSV->qrycontigid=%lld\n",
	       pSV->qrystartpos,linecnt,filename,QryMap->site[0][Nqry[0]+1], QryMap->id, vfixx_filename[num_files-1], pSV->qrycontigid);
	exit(1);
      }
    }
    if(SMAP && pSV->qryendpos > 0.0){
      pSV->qryendpos *= 0.001;
      if(pSV->qryendpos > QryMap->site[0][Nqry[0]+1] + 10.0){
	printf("QryEndPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s (+10kb)\n",
	       pSV->qryendpos,linecnt,filename,QryMap->site[0][Nqry[0]+1], QryMap->id, vfixx_filename[num_files-1]);
	exit(1);
      }
    }
    if(pSV->refstartpos > 0.0){
      pSV->refstartpos *= 0.001;
      if(pSV->refstartpos > RefMap->site[0][NRef[0]+1]){
	printf("RefStartPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s\n",
	       pSV->refstartpos,linecnt,filename,RefMap->site[0][NRef[0]+1], RefMap->id, vfixx_filename[num_files-2]);
	exit(1);
      }
    }
    if(pSV->refendpos > 0.0){
      pSV->refendpos *= 0.001;
      if(pSV->refendpos > RefMap->site[0][NRef[0]+1]){
	printf("RefEndPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s\n",
	       pSV->refendpos,linecnt,filename,RefMap->site[0][NRef[0]+1], RefMap->id, vfixx_filename[num_files-2]);
	exit(1);
      }
    }

    pSV->QryMap = QryMap;
    pSV->RefMap = RefMap;

    if(SMAP >= 6){/* scan ahead to RawConfidence column */
      pt = qt;
      while(*pt && isspace(*pt))
	pt++;

      /* skip over LinkID QryStartIdx QryEndIdx RefStartIdx RefEndIdx Zygosity Genotype GenotypeGroup */
      for(int skip = 0; skip < 8; skip++){
	while(*qt && !isspace(*pt))
	  pt++;
	while(*pt && isspace(*pt))
	  pt++;
      }
      
      double RawConfidence = strtod(pt,&qt);
      if(pt == qt || !(*qt == 0 || isspace(*qt))){
	printf("Invalid RawConfidence on line %d of %s: %s\n%s\n",linecnt,filename,pt,buf);
	fflush(stdout);exit(1);
      }
      
      //      pSV->confidence_scaled = pSV->confidence;
      pSV->confidence = RawConfidence;
    }

    if(SMAP){/* set pSV->orientforward based on xmapid */
      if(DEBUG) assert(1 <= pSV->xmapid && pSV->xmapid <= xmapCnt);
      pSV->orientforward = xmap[pSV->xmapid - 1].orientforward;
    }

    svCnt++;
  }
  (void)fclose(fp);
  if(VERB && SV2cnt > 0){
    printf("WARNING: Skipped %d Indels with two matchgroups\n",SV2cnt);
    fflush(stdout);
  }
}

static char buf_input_xmap[LINESIZ];

/* NOTE : If StartRef < 0 : read in Reference Maps and Query Maps (get their names from the .xmap), otherwise assume they are already read in from ref_filename and query_filename.
   Input maps are appended to map[0.. nummaps-1] : Query maps will be in map[StartQuery .. StartQuery + NummapsQuery-1] and Reference maps in map[StartRef .. StartRef + NummapsRef-1]. 

   Append xmap information to xmap[0..xmapMax-1] starting at xmap[xmapCnt] :  xmap[] is reallocated if needed and xmapCnt incremented by the number of xmaps read (Normally xmapCnt should be 0 initially)
] */

void input_xmap(char *filename, xmapEntry * &xmap, int &xmapCnt, int &xmapMax, int &StartQuery, int &NummapsQuery, int &StartRef, int &NummapsRef, char* &query_filename, char* &ref_filename)
{
  char *buf = buf_input_xmap;

  if(colors != 1){
    printf("input_xmap(): not yet implemented for -colors %d\n",colors);
    exit(1);
  }

  char *pt, *qt;
  FILE *fp;
  
  if((fp = fopen(filename,"r"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("Failed to read input file %s:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  int linecnt = 1;
  for(;fgets(buf,LINESIZ,fp) != NULL; linecnt++){
    int len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      exit(1);
    }
    
    if(buf[0] == '#'){/* most comment lines are ignored */
      char *key = (char *)"Label Channels:";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);
	int NColors = strtol(pt,&qt,10);
	if(qt == pt || !(*qt == 0 || isspace(*qt)) || NColors <= 0 || NColors > 2){
	  printf("Line %d in %s with %s has invalid integer value (must be 1 or 2):\n%s\n",linecnt, filename, key, buf);
	  exit(1);
	}
	if(NColors != colors){
	  printf("Label Channels on Line %d in %s does not match -colors %d\n",linecnt, filename, colors);
	  exit(1);
	}
	continue;
      }

      key = (char *)"Reference Maps From:";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);
	char name[PATH_MAX];
	int i;
	errno = 0;
	if(1 != (i=sscanf(pt,"%s",name))){
	  char *err  = strerror(errno);
	  printf("sscanf returned %d (expected 1) on line %d of %s:%s\n",i,linecnt,filename,err);
	  exit(1);
	}
	if(ref_filename != NULL){
	  char *pt = strrchr(ref_filename,'/');
	  char *qt = strrchr(name,'/');
	  if(!pt)
	    pt = ref_filename;
	  if(!qt)
	    qt = name;
	  if(strcmp(pt,qt)){
	    printf("Reference Maps filename %s on line %d in %s does not agree with previous filename %s (from .indel or .smap)\n",qt+1,linecnt,filename,ref_filename);
	    exit(1);
	  }
	} else
	  ref_filename = strdup(name);
	continue;
      }

      key = (char *)"Query Maps From:";
      if((pt = strstr(&buf[1],key))){
	pt += strlen(key);
	char name[PATH_MAX];
	int i;
	errno = 0;
	if(1 != (i=sscanf(pt,"%s",name))){
	  char *err  = strerror(errno);
	  printf("sscanf returned %d (expected 1) on line %d of %s:%s\n",i,linecnt,filename,err);
	  exit(1);
	}
	if(query_filename != NULL){
	  char *pt = strrchr(query_filename,'/');
	  char *qt = strrchr(name,'/');
	  if(!pt)
	    pt = query_filename;
	  if(!qt)
	    qt = name;
	  if(strcmp(pt,qt)){
	    printf("Query Maps filename %s on line %d in %s does not agree with previous filename %s (from .indel or .smap)\n",qt,linecnt,filename,query_filename);
	    exit(1);
	  }
	} else
	  query_filename = strdup(name);
	continue;
      }

      continue;
    }

    /* done with header comment lines */
    if(!ref_filename){
      printf("Reference Maps name not found in header lines of %s\n",filename);
      exit(1);
    }
    if(!query_filename){
      printf("Query Maps name not found in header lines of %s\n",filename);
      exit(1);
    }

    if(StartRef < 0){/* need to read in Reference and Query Maps */
      /* temporarily change global parameters so Maps are read in without rescaling or bias correction */
      double origPixelLen = PixelLen;
      double origMinLen = MinLen;
      double origMaxLen = MaxLen;
      int origMinSites = MinSites;
      int origMaxSites = MaxSites;
      int origResBins[MAXCOLOR] = {0,0};
      MinLen = MaxLen = 0.0;
      MaxSites = MinSites = 0;
      PixelLen = 0.500;/* Read in Reference Maps and Query maps without rescaling */
      for(int c = 0; c < colors; c++){
	origResBins[c] = ResBins[c];
	ResBins[c] = 0;/* Don't apply bias correction to Reference Maps and Query maps */
      }

#ifndef WIN32
      /* extract directory from absolute .xmap filename */
      char basename[PATH_MAX], *qt;
      if(filename[0] == '/')/* absolute path name */
	strcpy(basename,filename);
      else {
	if(getcwd(basename,PATH_MAX) == NULL){
	  int eno = errno;
	  char *err = strerror(eno);
	  printf("getcwd() failed (errno=%d):%s\n",eno, err);
	  exit(1);
	}
	pt = basename + strlen(basename);
	sprintf(pt,"/%s",filename);
      }
      if((pt = strrchr(basename,'/')))
	pt++;
      else {
	printf("Could not find '/' in absolute .xmap filename:%s\n",basename);
	exit(1);
      }

      struct stat statbuf;

      if(ref_filename[0] != '/' || stat(ref_filename, &statbuf)){/* change ref_filename into absolute pathname */
	if((qt = strrchr(ref_filename,'/')))
	  strcpy(pt,qt+1);
	else
	  strcpy(pt, ref_filename);
	if(VERB/* HERE >=2 */ && strcmp(ref_filename,basename)){
	  printf("Changing ref_filename from \"%s\" to \"%s\"\n",ref_filename,basename);
	  fflush(stdout);
	}
	free(ref_filename);
	ref_filename = strdup(basename);
      }
      if(query_filename[0] != '/' || stat(query_filename, &statbuf)){/* change query_filename into abolute pathname */
	if((qt = strrchr(query_filename,'/')))
	  strcpy(pt,qt+1);
	else
	  strcpy(pt, query_filename);
	if(VERB/* HERE >=2 */ && strcmp(query_filename,basename)){
	  printf("Changing query_filename from \"%s\" to \"%s\"\n", query_filename,basename);
	  fflush(stdout);
	}
	free(query_filename);
	query_filename = strdup(basename);
      }
#endif

      /* read in Reference Maps to map[nummaps ...] */
      if(num_files+1 >= MAXFILES){
	printf("Too many input files : increase MAXFILES=%d to at least %d\n", MAXFILES, num_files+2);
	exit(1);
      }
      vfixx_filename[num_files] = strdup(ref_filename);
      if(VERB>=2){
	printf("Reading Reference maps from %s(ref_filename=%s)\n",vfixx_filename[num_files],ref_filename);
	fflush(stdout);
      }
      StartRef = nummaps;
      NummapsRef = input_vfixx(num_files, nummaps, maxmaps, Gmap);
      if(DEBUG) assert(nummaps == StartRef + NummapsRef);
      num_files++;

      qsort(&Gmap[StartRef],NummapsRef,sizeof(Cmap *),(intcmp*)CmapIdInc);      

      if(VERB/* HERE >=2 */){
	printf("Read %d Reference maps from %s\n",NummapsRef,vfixx_filename[num_files-1]);
	fflush(stdout);
      }

      /* append Query Maps to map[nummaps ... ] */
      if(++num_files >= MAXFILES){
	printf("Too many input files : increase MAXFILES=%d to at least %d\n", MAXFILES, num_files+1);
	exit(1);
      }
      vfixx_filename[num_files-1] = strdup(query_filename);
      if(VERB>=2){
	printf("Reading Query maps from %s(query_filename=%s)\n",vfixx_filename[num_files-1],query_filename);
	fflush(stdout);
      }
      StartQuery = nummaps;
      NummapsQuery = input_vfixx(num_files-1,nummaps,maxmaps, Gmap);
      if(DEBUG) assert(nummaps == StartQuery + NummapsQuery);
      qsort(&Gmap[StartQuery],NummapsQuery,sizeof(Cmap *),(intcmp*)CmapIdInc);      

      if(VERB/* HERE >=2 */){
	printf("Read %d Query maps from %s\n",NummapsQuery,vfixx_filename[num_files-1]);
	fflush(stdout);
      }
      
      /* restore global parameters */
      for(int c = 0; c < colors; c++)
	ResBins[c] = origResBins[c];
      PixelLen = origPixelLen;
      MinLen = origMinLen;
      MaxLen = origMaxLen;
      MinSites = origMinSites;
      MaxSites = origMaxSites;
    }

    if(xmapCnt >= xmapMax){
      xmapMax = max(2*xmapMax,1024);
      xmapEntry *newxmap = new xmapEntry[xmapMax];
      if(xmapCnt)
	memcpy(newxmap,xmap,sizeof(xmapEntry)*xmapCnt);
      delete [] xmap;
      xmap = newxmap;
    }
    xmapEntry *pXmap = &xmap[xmapCnt];
    
    /* parse .xmap lines */
    pt = buf;
    int XmapEntryId = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || XmapEntryId <= 0){
      printf("Invalid XmapEntryId on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pXmap->smapentryid = XmapEntryId;
    if(DEBUG && XmapEntryId != xmapCnt+1){
      printf("Invalid XmapEntryId=%d on line %d of %s (should match %d'th xmap line)\n",XmapEntryId,linecnt,filename,xmapCnt+1);
      exit(1);
    }

    pt = qt;
    pXmap->qrycontigid = strtol(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pXmap->qrycontigid <= 0){
      printf("Invalid QryContigID on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    Cmap *QryMap = findidBS(&Gmap[StartQuery],NummapsQuery,(long long)pXmap->qrycontigid);
    if(!QryMap){
      printf("QryContigID=%lld on line %d of %s not found amoungs %d Query maps in %s\n",
	     pXmap->qrycontigid, linecnt, filename, NummapsQuery, query_filename);
      exit(1);
    }
    if(DEBUG) assert(QryMap != NULL && QryMap->id == pXmap->qrycontigid);
    int *Nqry = QryMap->numsite;

    pt = qt;
    pXmap->refcontigid = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pXmap->refcontigid <= 0){
      printf("Invalid RefContigID1 on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    Cmap *RefMap = findidBS(&Gmap[StartRef],NummapsRef,(long long)pXmap->refcontigid);
    if(!RefMap){
      printf("RefContigID1=%lld on line %d of %s not found amongst %d Reference maps in %s\n",
	     pXmap->refcontigid, linecnt, filename, NummapsRef, ref_filename);
      exit(1);
    }
    if(DEBUG) assert(RefMap != NULL && RefMap->id == pXmap->refcontigid);
    int *Nref = RefMap->numsite;

    pt = qt;
    pXmap->qrystartpos = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pXmap->qrystartpos < 0.0){
      printf("Invalid QryStartPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pXmap->qrystartpos *= 0.001;
    if(pXmap->qrystartpos > QryMap->site[0][Nqry[0]+1] + 1e-6){
      printf("QryStartPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s: qrycontigid=%lld\n",
	     pXmap->qrystartpos,linecnt,filename,QryMap->site[0][Nqry[0]+1], QryMap->id, query_filename, pXmap->qrycontigid);
      exit(1);
    }

    pt = qt;
    pXmap->qryendpos = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pXmap->qryendpos < 0.0){
      printf("Invalid QryEndPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pXmap->qryendpos *= 0.001;
    if(pXmap->qryendpos > QryMap->site[0][Nqry[0]+1]){
      printf("QryEndPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s: qrycontigid=%lld\n",
	     pXmap->qryendpos,linecnt,filename,QryMap->site[0][Nqry[0]+1], QryMap->id, query_filename, pXmap->qrycontigid);
      exit(1);
    }

    pt = qt;
    pXmap->refstartpos = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pXmap->refstartpos < 0.0){
      printf("Invalid RefStartPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pXmap->refstartpos *= 0.001;
    if(pXmap->refstartpos > RefMap->site[0][Nref[0]+1]){
      printf("RefStartPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s\n",
	     pXmap->refstartpos, linecnt, filename, RefMap->site[0][Nref[0]+1], RefMap->id, ref_filename);
      exit(1);
    }

    pt = qt;
    pXmap->refendpos = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt)) || pXmap->refendpos < 0.0){
      printf("Invalid RefEndPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pXmap->refendpos *= 0.001;
    if(pXmap->refendpos > RefMap->site[0][Nref[0]+1]){
      printf("RefEndPos=%0.3f kB on line %d of %s is larger than size=%0.3f kB of map id=%lld in %s\n",
	     pXmap->refendpos,linecnt,filename,RefMap->site[0][Nref[0]+1], RefMap->id, ref_filename);
      exit(1);
    }

    /* read orientation character */
    while(*qt && isspace(*qt))
      qt++;
    if(*qt != '+' && *qt != '-'){
      printf("Invalid Orientation=%c on line %d of %s:\n%s\n",*qt,linecnt,filename,buf);
      exit(1);
    }
    pXmap->orientforward = (*qt == '+');
    while(*qt && !isspace(*qt))
      qt++;
    
    pt = qt;
    pXmap->confidence = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt))){
      printf("Invalid Confidence on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }

    /* skip CIGAR string */
    while(*qt && isspace(*qt))
      qt++;
    while(*qt && !isspace(*qt))
      qt++;
    
    pt = qt;
    pXmap->querycontiglen = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt))){
      printf("Invalid QryLen on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    
    pt = qt;
    pXmap->refcontiglen = strtod(pt,&qt);
    if(pt == qt || !(*qt==0 || isspace(*qt))){
      printf("Invalid RefLen on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    int LabelChannel = strtol(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || LabelChannel != 1){
      printf("Invalid LabelChannel on line %d of %s (must be 1):\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    
    /* Allocate alignment */
    Calign *palign = pXmap->align = new Calign[1];
    palign->mapid1 = RefMap->mapid;
    palign->mapid2 = QryMap->mapid;
    palign->orientation = (pXmap->orientforward ? 0 : 1);
    palign->allocated_size = 32;
    palign->numpairs = 0;
    palign->sites1 = new int[palign->allocated_size];
    palign->sitesK1 = new int[palign->allocated_size];
    palign->sites2 = new int[palign->allocated_size];
    
    /* parse alignment */
    pt = qt;
    int lastI = -1, lastJ = -1;
    while((qt = strchr(pt,'('))){
      pt = qt+1;
      int I = strtol(pt,&qt,10);
      if(pt==qt || !(*qt==0 || isspace(*qt) || *qt==',')){
	printf("Invalid Ref index on line %d of %s starting at %s:\n%s\n",linecnt,filename,pt,buf);
	exit(1);
      }
      if(I <= 0 || I > Nref[0]){
	printf("Out of range Ref index %d on line %d of %s (N=%d):\n%s\n",I,linecnt,filename,Nref[0],buf);
	exit(1);
      }

      while(*qt && isspace(*qt))
	qt++;
      if(*qt != ','){
	printf("Expected ',' on line %d of %s starting at %s:\n%s\n",linecnt,filename,qt,buf);
	exit(1);
      }

      pt = qt+1;
      int J = strtol(pt,&qt,10);
      if(pt==qt || !(*qt==0 || isspace(*qt) || *qt==')')){
	printf("Invalid Query index on line %d of %s starting at %s:\n%s\n",linecnt,filename,&pt[-1],buf);
	exit(1);
      }
      if(J <= 0 || J > Nqry[0]){
	printf("Out of range Query index %d on line %d of %s (M=%d):\n%s\n",J,linecnt,filename,Nqry[0],buf);
	exit(1);
      }

      while(*qt && isspace(*qt))
	qt++;
      if(*qt != ')'){
	printf("Expected ')' on line %d of %s starting at %s:\n%s\n",linecnt,filename,qt,buf);
	exit(1);
      }
      qt++;
      while(*qt && isspace(*qt))
	qt++;
      pt = qt;

      if(J == lastJ){
	if(DEBUG) assert(palign->numpairs >= 1);
	int U = palign->numpairs - 1;
	if(DEBUG) assert(palign->sites1[U] == lastI);
	if(DEBUG) assert(palign->sitesK1[U] == 0);
	if(DEBUG && !(palign->sites2[U] == (palign->orientation ? (Nqry[0]+1 - lastJ) : lastJ))){
	  printf("XmapID=%d:QryID=%lld,RefID=%lld,or=%d:QryStart=%0.3f,QryEnd=%0.3f,RefStart=%0.3f,RefEnd=%0.3f:numpairs=%d:I=%d,J=%d,lastI=%d,lastJ=%d,U=%d,palign->sites2[U]=%d,Nqry[0]=%d\n%s\n", 
		 XmapEntryId,pXmap->qrycontigid,pXmap->refcontigid,palign->orientation, pXmap->qrystartpos, pXmap->qryendpos, pXmap->refstartpos, pXmap->refendpos, palign->numpairs, I,J,lastI,lastJ,
		 U,palign->sites2[U],Nqry[0],buf);
	  fflush(stdout);

	  assert(palign->sites2[U] == (palign->orientation ? (Nqry[0]+1 - lastJ) : lastJ));
	}
	if(DEBUG) assert(I > lastI);
	
	palign->sitesK1[U] = I - lastI;
	palign->sites1[U] = I;
      } else {
	int U = palign->numpairs;
	if(++palign->numpairs >= palign->allocated_size){
	  palign->allocated_size *= 2;
	  int *sites1 = new int[palign->allocated_size];
	  int *sitesK1 = new int[palign->allocated_size];
	  int *sites2 = new int[palign->allocated_size];
	  memcpy(sites1,palign->sites1,U*sizeof(int));
	  memcpy(sitesK1,palign->sitesK1,U*sizeof(int));
	  memcpy(sites2,palign->sites2,U*sizeof(int));
	  delete [] palign->sites1;
	  delete [] palign->sitesK1;
	  delete [] palign->sites2;
	  palign->sites1 = sites1;
	  palign->sitesK1 = sitesK1;
	  palign->sites2 = sites2;
	}
	if(DEBUG) assert(I > lastI);
	if(DEBUG && lastJ > 0 && !(pXmap->orientforward ? (J > lastJ) : (J < lastJ))){
	  printf("XmapID=%d:QryID=%lld,RefID=%lld,or=%d:QryStart=%0.3f,QryEnd=%0.3f,RefStart=%0.3f,RefEnd=%0.3f:numpairs=%d:I=%d,J=%d,lastI=%d,lastJ=%d\n%s\n", 
		 XmapEntryId,pXmap->qrycontigid,pXmap->refcontigid,palign->orientation, pXmap->qrystartpos, pXmap->qryendpos, pXmap->refstartpos, pXmap->refendpos, palign->numpairs, I,J,lastI,lastJ,buf);
	  fflush(stdout);
	  assert(pXmap->orientforward ? (J > lastJ) : (J < lastJ));
	}
	palign->sites1[U] = I;
	palign->sitesK1[U] = 0;
	palign->sites2[U] = palign->orientation ? (Nqry[0]+1 - J) : J;
      }
      lastI = I;
      lastJ = J;
    }
    pXmap->QryMap = QryMap;
    pXmap->RefMap = RefMap;
    xmapCnt++;
  }
  if(VERB){
    printf("Read in %d alignments from %s\n",xmapCnt, filename);
    fflush(stdout);
  }
}

/* return 1 IFF there is an xmap entry in xmap[0..xmapCnt-1] that completely overlaps the range refstart .. refend in reference contig refcontigid */
/* Requires that xmap[] is sorted in ascending order of refcontigid, refstartpos and refendpos */

int findXmapOverlap(xmapEntry *xmap, int xmapCnt, int refcontigid, double refstartpos, double refendpos)
{
  /* use binary search to locate xmap entry (if it exists) */
  int low = 0;
  int high = xmapCnt - 1;
  while(low <= high){
    int mid = (low+high)/2;
    xmapEntry *p = &xmap[mid];

    /* first check refcontigid */
    if(p->refcontigid < refcontigid){
      low = mid + 1;
      continue;
    }
    if(p->refcontigid > refcontigid){
      high = mid - 1;
      continue;
    }

    /* next check refstartpos,refendpos */
    if(p->refendpos < refendpos){
      low = mid + 1;
      continue;
    }
    if(p->refstartpos > refstartpos){
      high = mid - 1;
      continue;
    }

    if(DEBUG/* HERE HERE >=2 */) assert(p->refcontigid == refcontigid && p->refendpos >= refendpos && p->refstartpos <= refstartpos);
    return 1;
  }

  return 0;// not found
}

#define SIZEBINS 29
static double SizeValues[SIZEBINS] = {0.0,0.1,0.2,0.3,0.5,0.7,1.0,1.4,2.0,3.0,5.0,7.0,10.0,14.0,20.0,30.0,50.0,70.0,100.0,140.0,200.0,300.0,500.0,700.0,1000.0,1400.0,2000.0,3000.0,5000.0};

/* encode indel size (+ve for insertions -ve for deletions) as a small integer (+ve for insertions, -ve for deletions) */
/* An absolute bin value of N=1..SIZEBINS-1 corresponds to the interval SizeValues[N-1] .. SizeValues[N] and value N=SIZEBINS correponds to all values > SizeValues[SIZEBINS-1] */
static int sizebin(double size)
{
  int sn = (size > 0.0) ? 1 : -1;
  double asiz = fabs(size);
  
  if(asiz >= SizeValues[SIZEBINS-1])
    return SIZEBINS * sn;

  /* locate bin using binary search */
  int high = SIZEBINS-1;/* value is below SizeValues[high] */
  int low = 0;/* value is at or above SizeValues[low] */
  while(high - low > 1){
    int mid = (low + high) / 2;
    if(asiz >= SizeValues[mid])
      low = mid;
    else
      high = mid;
  }
  if(DEBUG) assert(1 <= high && high <= SIZEBINS && abs(sn) == 1);
  return high * sn;
}

/* return 1 IFF the two SVs are considered a match. Assume SVs are in ascending order of refcontigid,refstartpos,refendpos */
static int checkmatch(xmapEntry *SV1, int j1, int svCnt1,
		      xmapEntry *SV2, int j2, int svCnt2,
		      double& overlapstart,
		      double& overlapend)
{
  long long refcontigid1 = SV1[j1].refcontigid;
  long long refcontigid2 = SV2[j2].refcontigid;

  if(DEBUG/* HERE >=2 */) assert(refcontigid1 == refcontigid2);

  double refstartpos1 = SV1[j1].refstartpos;
  double refendpos1 = SV1[j1].refendpos;
  double qrystartpos1 = SV1[j1].qrystartpos;
  double qryendpos1 = SV1[j1].qryendpos;

  double refstartpos2 = SV2[j2].refstartpos;
  double refendpos2 = SV2[j2].refendpos;
  double qrystartpos2 = SV2[j2].qrystartpos;
  double qryendpos2 = SV2[j2].qryendpos;
    
  double delta1 = (refendpos1 - refstartpos1) - fabs(qryendpos1 - qrystartpos1);
  double delta2 = (refendpos2 - refstartpos2) - fabs(qryendpos2 - qrystartpos2);
  double len1 = max(fabs(qryendpos1 - qrystartpos1),refendpos1 - refstartpos1);
  double len2 = max(fabs(qryendpos2 - qryendpos1), refendpos2 - refstartpos2);
  double expand1 = max(0.0, fabs(qryendpos1-qrystartpos1) - (refendpos1-refstartpos1));
  double expand2 = max(0.0, fabs(qryendpos2-qrystartpos2) - (refendpos2-refstartpos2));
  overlapstart = max(refstartpos1 - expand1*0.5, refstartpos2 - expand2*0.5);
  overlapend = min(refendpos1 + expand1*0.5, refendpos2 + expand2*0.5);

  int match = 0;

  if(overlapend - overlapstart >= SVcompareKB && overlapend - overlapstart >= min(len1,len2) * SVcompareFrac && 
     (SVDeltaMatch <= 0.0 || fabs(delta1-delta2) <= max(fabs(delta1),fabs(delta2)) * SVDeltaMatch)){
    match = 1;
    if(overlapend < overlapstart){
      if(SV1 == SV2)
	match = 0;
      else {/* check if next SV is closer */
	if(DEBUG) assert(min(refendpos1,refendpos2) < max(refstartpos1,refstartpos2));
	if(refendpos1 < refstartpos2){
	  /* check if SV1[j1+1] is closer to SV2[j2] (than SV1[j1]) */
	  if(j1+1 < svCnt1 && SV1[j1+1].refcontigid == refcontigid1){
	    double nrefstartpos1 = SV1[j1+1].refstartpos;
	    if(nrefstartpos1 - refendpos2 < refstartpos2 - refendpos1)
	      match = 0;
	  }
	  /* check if SV2[j2-1] is closer to SV1[j1] ( than SV2[j2]) */
	  if(j2 > 0 && SV1[j2-1].refcontigid == refcontigid2){
	    double prefendpos2 = SV2[j2-1].refendpos;
	    if(refstartpos1 - prefendpos2 < refstartpos2 - refendpos1)
	      match = 0;
	  }
	} else {
	  if(DEBUG) assert(refstartpos1 > refendpos2);

	  /* check if SV2[j2+1] is closer to SV1[j1] (than SV2[j2]) */
	  if(j2+1 < svCnt2 && SV2[j2+1].refcontigid == refcontigid2){
	    double nrefstartpos2 = SV2[j2+1].refstartpos;
	    if(nrefstartpos2 - refendpos1 < refstartpos1 - refendpos2)
	      match = 0;
	  }

	  /* check if SV1[j1-1] is closer to SV2[j2] ( than SV1[j1]) */
	  if(j1 > 0 && SV1[j1-1].refcontigid == refcontigid1){
	    double prefendpos1 = SV1[j1-1].refendpos;
	    if(refstartpos2 - prefendpos1 < refstartpos1 - refendpos2)
	      match = 0;
	  }
	}
      }
    }
  }

  if(VERB>=3 && refcontigid1 == 2){
    printf("j1=%d,j2=%d,match=%d:\n",j1,j2,match);
    printf("  SV1[%d]:queryid=%lld refid=%lld refstart=%0.3f,refend=%0.3f,qrystart=%0.3f,qryend=%0.3f,type=%s,confidence=%0.2f\n",
	   j1, SV1[j1].qrycontigid,SV1[j1].refcontigid,SV1[j1].refstartpos,SV1[j1].refendpos,SV1[j1].qrystartpos,SV1[j1].qryendpos,SV1[j1].type,SV1[j1].confidence);
    printf("  SV2[%d]:queryid=%lld refid=%lld refstart=%0.3f,refend=%0.3f,qrystart=%0.3f,qryend=%0.3f,type=%s,confidence=%0.2f\n",
	   j2, SV2[j2].qrycontigid,SV2[j2].refcontigid,SV2[j2].refstartpos,SV2[j2].refendpos,SV2[j2].qrystartpos,SV2[j2].qryendpos,SV2[j2].type,SV2[j2].confidence);
	     
  }

  return match;
}

/* display comparison of 2 indel sets (must be for same reference mapset) */
void indel_compare(char *filename1, char *filename2, double SVcompareKB, double SVcompareFrac)
{
  if(SVSiteMatch != 0){
    printf("-svcompare : SiteMatch=%d not yet implemented\n", SVSiteMatch);
    exit(1);
  }
  /*  if(SVXmapOverlap > 0.0){
    printf("-svcompare : SVXmapOverlap=%0.3f not yet implemented\n", SVXmapOverlap);
    exit(1);
    }*/

  svMax = svCnt = 0;
  SVtocheck = 0;
  input_indel(filename1, SVtocheck, svCnt, svMax);/* read in the SVs and append query & ref maps to map[0..nummaps-1] */

  if(VERB){
    printf("Read in %d SVs from %s\n", svCnt, filename1);
    fflush(stdout);
  }

  /* save 1st set of SVs */
  char *xmap_filename1 = xmap_filename;
  char *query_filename1 = query_filename;
  char *ref_filename1 = ref_filename;
  int StartRef1 = StartRef;
  int NummapsRef1 = NummapsRef;
  int StartQuery1 = StartQuery;
  int NummapsQuery1 = NummapsQuery;
  //  int svMax1 = svMax;
  int svCnt1 = svCnt;
  xmapEntry *SV1 = SVtocheck;

  // input the xmaps as well: helps categorize unique SVs AND detect mismatched sites
  xmapEntry *xmap1 = 0;
  int xmapCnt1 = 0, xmapMax1 = 0;
  input_xmap(xmap_filename1, xmap1, xmapCnt1, xmapMax1, StartQuery1, NummapsQuery1, StartRef1, NummapsRef1, query_filename1, ref_filename1);

  xmap_filename = query_filename = ref_filename = 0;
  StartRef = NummapsRef = StartQuery = NummapsQuery = -1;
  svMax = svCnt = 0;
  SVtocheck = 0;
  input_indel(filename2, SVtocheck, svCnt, svMax);/* read in the SVs and append query & ref maps  to map[0..nummaps-1] */
  if(VERB){
    printf("Read in %d SVs from %s\n", svCnt, filename2);
    fflush(stdout);
  }

  char *xmap_filename2 = xmap_filename;
  char *query_filename2 = query_filename;
  char *ref_filename2 = ref_filename;
  int StartRef2 = StartRef;
  int NummapsRef2 = NummapsRef;
  int StartQuery2 = StartQuery;
  int NummapsQuery2 = NummapsQuery;
  //  int svMax2 = svMax;
  int svCnt2 = svCnt;
  xmapEntry *SV2 = SVtocheck;

  // input the xmaps as well: helps categorize unique SVs AND detect mismatched sites
  xmapEntry *xmap2 = 0;
  int xmapCnt2 = 0, xmapMax2 = 0;
  input_xmap(xmap_filename2, xmap2, xmapCnt2, xmapMax2, StartQuery2, NummapsQuery2, StartRef2, NummapsRef2, query_filename2, ref_filename2);

  if(strcmp(ref_filename2,ref_filename1)){
    printf("WARNING: -svcompare %s %s : _r.cmaps are not the same (%s and %s respectively)\n",
	   filename1, filename2, ref_filename1,ref_filename2);
    printf("   (Assuming the files have the same reference maps)\n");
    fflush(stdout);
  }

  /* sort both sets of SVs by reference id and location then query id (id1, sites1[0], sites1[N], id2) */ 
  qsort(SV1,svCnt1,sizeof(xmapEntry), (intcmp*)Id1Site1Id2Inc);  
  qsort(SV2,svCnt2,sizeof(xmapEntry), (intcmp*)Id1Site1Id2Inc);  

  double overlapstart,overlapend;

  /* merge overlapped SVs in each set and remove SVs below LogPvThreshold */
  int i = 0;
  for(int j1 = 0; j1 < svCnt1; j1++){
    if(DEBUG) assert(isfinite(SV1[j1].confidence));
    if(SV1[j1].confidence < LogPvThreshold)
      continue;

    int k = i;
    for(; --k >= 0;){
      if(SV1[k].refcontigid != SV1[j1].refcontigid)
	break;
      /* check if SV[k] is overlapped with SV[j1] */
      if(checkmatch(SV1,k,svCnt1,SV1,j1,svCnt1,overlapstart, overlapend)){
	if(VERB>=2 && SV1[k].refcontigid == 1 && (SV1[k].qrycontigid == 14720 || SV1[j1].qrycontigid == 14720)){
	  printf("SV1[%d]:queryid=%lld,refid=%lld,fwd=%d,conf=%0.2f,refstart=%0.3f,refend=%0.3f; SV1[%d]:queryid=%lld,refid=%lld,fwd=%d,conf=%0.2f,refstart=%0.3f,refend=%0.3f: Overlaped(start=%0.3f,end=%0.3f)\n",
		 k,SV1[k].qrycontigid,SV1[k].refcontigid,SV1[k].orientforward?1:0,SV1[k].confidence,SV1[k].refstartpos,SV1[k].refendpos,
		 j1,SV1[j1].qrycontigid,SV1[j1].refcontigid,SV1[j1].orientforward?1:0,SV1[j1].confidence,SV1[j1].refstartpos,SV1[j1].refendpos,
		 overlapstart, overlapend);
	  fflush(stdout);
	}
	if(SV1[j1].confidence > SV1[k].confidence) /* save SV with larger confidence */
	  SV1[k] = SV1[j1];
	break;
      }
    }
    if(k >= 0 && SV1[k].refcontigid == SV1[j1].refcontigid)
      continue;/* match found */

    if(i < j1)
      SV1[i] = SV1[j1];
    i++;
  }
  if(VERB/* HERE >=2 */ && i < svCnt1){
    printf("%d SVs from %s merged and reduced to %d\n",svCnt1,filename1, i);
    fflush(stdout);
  }
  svCnt1 = i;

  for(int i = 0; i < svCnt1; i++){
    /* compute indel size bin for SV1[i] */
    double refstartpos = SV1[i].refstartpos;
    double refendpos = SV1[i].refendpos;
    double qrystartpos = SV1[i].qrystartpos;
    double qryendpos = SV1[i].qryendpos;
    double delta = (refendpos - refstartpos) - fabs(qryendpos - qrystartpos);
    SV1[i].sbin = sizebin(-delta);
  }
  
  qsort(SV1,svCnt1,sizeof(xmapEntry), (intcmp*)Id1Site1Id2Inc);  /* sort again since order may no longer be correct */

  if(VERB>=2){
    for(int j1 = 0; j1 < svCnt1; j1++)
      printf("  SV1[%d]:queryid=%lld refid=%lld refstart=%0.3f,refend=%0.3f,qrystart=%0.3f,qryend=%0.3f,type=%s,confidence=%0.2f,sbin=%d\n",
	     j1, SV1[j1].qrycontigid,SV1[j1].refcontigid,SV1[j1].refstartpos,SV1[j1].refendpos,SV1[j1].qrystartpos,SV1[j1].qryendpos,SV1[j1].type,SV1[j1].confidence,SV1[j1].sbin);
  }

  i = 0;
  for(int j2 = 0; j2 < svCnt2; j2++){
    if(DEBUG) assert(isfinite(SV2[j2].confidence));

    if(SV2[j2].confidence < LogPvThreshold)
      continue;

    int k = i;
    for(; --k >= 0;){
      if(SV2[k].refcontigid != SV2[j2].refcontigid)
	break;
      /* check if SV[k] is overlapped with SV[j2] */
      if(checkmatch(SV2, k, svCnt2, SV2, j2, svCnt2, overlapstart, overlapend)){
	if(VERB>=2 && SV1[k].refcontigid == 1 && (SV1[k].qrycontigid == 14720 || SV1[j2].qrycontigid == 14720)){
	  printf("SV2[%d]:queryid=%lld,refid=%lld,fwd=%d,conf=%0.2f,refstart=%0.3f,refend=%0.3f; SV2[%d]:queryid=%lld,refid=%lld,fwd=%d,conf=%0.2f,refstart=%0.3f,refend=%0.3f: Overlaped(start=%0.3f,end=%0.3f)\n",
		 k,SV2[k].qrycontigid,SV2[k].refcontigid,SV2[k].orientforward?1:0,SV2[k].confidence,SV2[k].refstartpos,SV2[k].refendpos,
		 j2,SV2[j2].qrycontigid,SV2[j2].refcontigid,SV2[j2].orientforward?1:0,SV2[j2].confidence,SV2[j2].refstartpos,SV2[j2].refendpos,
		 overlapstart, overlapend);
	  fflush(stdout);
	}
	if(SV2[j2].confidence > SV2[k].confidence) /* save SV with higher confidence */
	  SV2[k] = SV2[j2];
	break;
      }
    }
    if(k >= 0 && SV2[k].refcontigid == SV2[j2].refcontigid)
      continue;/* match found */

    if(i < j2)
      SV2[i] = SV2[j2];
    i++;
  }

  if(VERB/* HERE >=2 */ && i < svCnt2){
    printf("%d SVs from %s merged and reduced to %d\n",svCnt2,filename2, i);
    fflush(stdout);
  }
  svCnt2 = i;

  for(int i = 0; i < svCnt2; i++){
    /* compute indel size bin for SV2[i] */
    double refstartpos = SV2[i].refstartpos;
    double refendpos = SV2[i].refendpos;
    double qrystartpos = SV2[i].qrystartpos;
    double qryendpos = SV2[i].qryendpos;
    double delta = (refendpos - refstartpos) - fabs(qryendpos - qrystartpos);
    SV2[i].sbin = sizebin(-delta);
  }

  qsort(SV2,svCnt2,sizeof(xmapEntry), (intcmp*)Id1Site1Id2Inc);  /* sort again since order may no longer be correct */

  if(VERB>=2)
    for(int j2 = 0; j2 < svCnt2; j2++)
      printf("  SV2[%d]:queryid=%lld refid=%lld refstart=%0.3f,refend=%0.3f,qrystart=%0.3f,qryend=%0.3f,type=%s,confidence=%0.2f,sbin=%d\n",
	     j2, SV2[j2].qrycontigid, SV2[j2].refcontigid, SV2[j2].refstartpos, SV2[j2].refendpos, SV2[j2].qrystartpos, SV2[j2].qryendpos, SV2[j2].type, SV2[j2].confidence, SV2[j2].sbin);

  double matchlen = 0.0;  /* also compute kb intersection in reference intervals between two sets of SVs : first requires merging overlapped SVs in each set */
  for(int j1 = 0; j1 < svCnt1; j1++)
    SV1[j1].sv_overlap = false;
  for(int j2 = 0; j2 < svCnt2; j2++)
    SV2[j2].sv_overlap = false;

  int matchcnt = 0;/* number of matching SVs with both confidence values >= LogPvThreshold */
  int tcnt1 = 0;/* number of SV1[] with confidence < LogPvThreshold */
  int tcntPV1 = 0;/* number of SV1[] with confidence < max(LogPvThreshold, SVConf) */
  int tcnt2 = 0;/* number of SV2[] with confidence < LogPvThreshold */
  int tcntPV2 = 0;/* number of SV2[] with confidence < max(LogPvThreshold, SVConf) */

  for(int j1 = 0; j1 < svCnt1; j1++)
    if(SV1[j1].confidence < max(LogPvThreshold, SVConf))
      tcntPV1++;
  for(int j2 = 0; j2 < svCnt2; j2++)
    if(SV2[j2].confidence < max(LogPvThreshold, SVConf))
      tcntPV2++;

  if(VERB>=2){
    printf("SVcompareKB=%0.3f, SVcompareFrac=%0.6f, SVConf=%0.2f, SVDeltaMatch=%0.6f, SVSiteMatch=%d, SVXmapOverlap= %0.3f\n",SVcompareKB,SVcompareFrac,SVConf,SVDeltaMatch,SVSiteMatch, SVXmapOverlap);
    fflush(stdout);
  }

  int k2 = 0;
  for(int j1 = 0; j1 < svCnt1; j1++){
    if(DEBUG) assert(SV1[j1].confidence >= LogPvThreshold);

    while(k2 < svCnt2 && SV2[k2].refcontigid < SV1[j1].refcontigid)
      k2++;
    
    for(int j2 = k2; j2 < svCnt2; j2++){
      if(DEBUG) assert(SV2[j2].confidence >= LogPvThreshold);

      if(SV2[j2].refcontigid > SV1[j1].refcontigid)
	break;

      if(checkmatch(SV1,j1, svCnt1, SV2,j2,svCnt2, overlapstart, overlapend)){
	matchcnt++;
	SV1[j1].sv_overlap = true;
	SV2[j2].sv_overlap = true;

	if(overlapstart < overlapend)
	  matchlen += overlapend - overlapstart;
      }
    }
  }
  
  int unique1 = 0;/* number of unmatched SV1[] with confidence >= LogPvThreshold */
  int unique2 = 0;/* number of unmatched SV2[] with confidence >= LogPvThreshold */
  int uniqueFP1 = 0;/* subset of unique1 that has a xmap2[] entry that overlaps SV by at least SVXmapOverlap on both sides */
  int uniqueFP2 = 0;/* subset of unique2 that has a xmap1[] entry that overlaps SV by at least SVXmapOverlap on both sides */

  double ConfSum1[2] = {0.0,0.0};/* sum of SV1[] confidence values for matched(index 1) and unmatched(index 0) SVs with confidence >= LogPvThreshold */
  double ConfSum2[2] = {0.0,0.0};/* sum of SV2[] confidence values for matched(index 1) and unmatched(index 0) SVs with confidence >= LogPvThreshold */

  double totlen1 = 0.0;/* sum of SV1[] ref lengths for SVs with confidence >= LogPvThreshold */
  double totlen2 = 0.0;/* sum of SV2[] ref lengths for SVs with confidence >= LogPvThreshold */

  int uniquePV1 = 0;/* number of unmatched SV1[] with confidence >= max(LogPvThreshold,SVConf) */
  int uniquePV2 = 0;/* number of unmatched SV2[] with confidence >= max(LogPvThreshold,SVConf) */
  int totPV1 = 0;/* number of SV1[] with confidence >= max(LogPvThreshold,SVConf) */
  int totPV2 = 0;/* number of SV2[] with confidence >= max(LogPvThreshold,SVConf) */
  int uniquePVFP1 = 0;/* subset of uniquePV1 that has a xmap2[] entry that overlaps SV by at least SVXmapOverlap on both sides */
  int uniquePVFP2 = 0;/* subset of uniquePV2 that has a xmap1[] entry that overlaps SV by at least SVXmapOverlap on both sides */
  int uniquePV1siz[2*SIZEBINS+1];/* uniquePV1siz[i] is number of unmatched SV1[] with confidence >= max(LogPvThreshold,SVConf) for size bin i-(SIZEBINS+1) */
  int uniquePV2siz[2*SIZEBINS+1];/* uniquePV2siz[i] is number of unmatched SV2[] with confidence >= max(LogPvThreshold,SVConf) for size bin i-(SIZEBINS+1) */
  int totPV1siz[2*SIZEBINS+1];/* totPV1siz[i] is number of total SV1[] with confidence >= max(LogPvThreshold,SVConf) for size bin i-(SIZEBINS+1) */
  int totPV2siz[2*SIZEBINS+1];/* totPV2siz[i] is number of total SV2[] with confidence >= max(LogPvThreshold,SVConf) for size bin i-(SIZEBINS+1) */
  double uniqueRlenPV1siz[2*SIZEBINS+1];/* uniqueRlenPV1siz[i] is the sum of Reference intervals lengths for unmatched SV1[] with confidence >= max(LogPvThreshold,SVConf) for size bin i-(SIZEBINS+1) */
  double uniqueRlenPV2siz[2*SIZEBINS+1];/* uniqueRlenPV2siz[i] is the sum of Reference intervals lengths for unmatched SV2[] with confidence >= max(LogPvThreshold,SVConf) for size bin i-(SIZEBINS+1) */
  double totRlenPV1siz[2*SIZEBINS+1];/* totRlenPV1siz[i] is the sum of Reference intervals lengths for total SV1[] with confidence >= max(LogPvThreshold,SVConf) for size bin i-(SIZEBINS+1) */
  double totRlenPV2siz[2*SIZEBINS+1];/* totRlenPV2siz[i] is the sum of Reference intervals lengths for total SV2[] with confidence >= max(LogPvThreshold,SVConf) for size bin i-(SIZEBINS+1) */

  for(int i = 0; i <= 2*SIZEBINS; i++){
    uniquePV1siz[i] = uniquePV2siz[i] = totPV1siz[i] = totPV2siz[i] = 0;
    uniqueRlenPV1siz[i] = uniqueRlenPV2siz[i] = totRlenPV1siz[i] = totRlenPV2siz[i] = 0.0;
  }

  double ConfSumPV1[2] = {0.0,0.0};/* sum of SV1[] confidence values for matched(index 1) and unmatched(index 0) SVs with confidence >= max(LogPvThreshold,SVConf) */
  double ConfSumPV2[2] = {0.0,0.0};/* sum of SV2[] confidence values for matched(index 1) and unmatched(index 0) SVs with confidence >= max(LogPvThreshold,SVConf) */

  double totlenPV1 = 0.0;/* sum of SV1[] ref lengths for SVs with confidence >= max(LogPvThreshold,SVConf) */
  double totlenPV2 = 0.0;/* sum of SV2[] ref lengths for SVs with confidence >= max(LogPvThreshold,SVConf) */

  qsort(xmap1,xmapCnt1,sizeof(xmapEntry),(intcmp*)xmapEntryRefInc);
  qsort(xmap2,xmapCnt2,sizeof(xmapEntry),(intcmp*)xmapEntryRefInc);

  if(SVXmapOverlap > 0.0){/* set field have_ref_repeat == true to indicate if SV is confirmed by presence of xmap in other assembly */
    for(int j1 = 0; j1 < svCnt1; j1++){
      SV1[j1].have_ref_repeat = false;
      if(findXmapOverlap(xmap2, xmapCnt2, SV1[j1].refcontigid, SV1[j1].refstartpos - SVXmapOverlap, SV1[j1].refendpos + SVXmapOverlap))
	SV1[j1].have_ref_repeat = true;
    }
    for(int j2 = 0; j2 < svCnt2; j2++){
      SV2[j2].have_ref_repeat = false;
      if(findXmapOverlap(xmap1, xmapCnt1, SV2[j2].refcontigid, SV2[j2].refstartpos - SVXmapOverlap, SV2[j2].refendpos + SVXmapOverlap))
	SV2[j2].have_ref_repeat = true;
    }
  }


  for(int j1 = 0; j1 < svCnt1; j1++){
    if(SV1[j1].confidence >= LogPvThreshold){
      unique1 += SV1[j1].sv_overlap ? 0 : 1;
      if(SVXmapOverlap > 0.0 && !SV1[j1].sv_overlap && SV1[j1].have_ref_repeat)
	uniqueFP1++;
      ConfSum1[SV1[j1].sv_overlap ? 1 : 0] += SV1[j1].confidence;
      if(DEBUG && !(isfinite(ConfSum1[1]))){
	printf("j1=%d/%d:SV1[j1].confidence= %0.6e, SV1[j1].sv_overlap=%d,ConfSum1[1]=%0.6e\n",
	       j1,svCnt1,SV1[j1].confidence, SV2[j1].sv_overlap ? 1 : 0, ConfSum1[1]);
	fflush(stdout);
	assert(isfinite(ConfSum1[1]));
      }
      totlen1 += SV1[j1].refendpos - SV1[j1].refstartpos;

      if(SV1[j1].confidence >= SVConf){
	uniquePV1 += SV1[j1].sv_overlap ? 0 : 1;
	totPV1++;
	if(SIZEBINS > 0){
	  if(DEBUG) assert(SV1[j1].sbin != 0 && abs(SV1[j1].sbin) <= SIZEBINS);
	  int bin = SV1[j1].sbin + SIZEBINS;
	  double Rlen = SV1[j1].refendpos - SV1[j1].refstartpos;
	  totPV1siz[bin]++;
	  totRlenPV1siz[bin] += Rlen;
	  if(!SV1[j1].sv_overlap){
	    uniquePV1siz[bin]++;
	    uniqueRlenPV1siz[bin] += Rlen;
	  }
	}

	if(SVXmapOverlap > 0.0 && !SV1[j1].sv_overlap && SV1[j1].have_ref_repeat)
	  uniquePVFP1++;
	ConfSumPV1[SV1[j1].sv_overlap ? 1 : 0] += SV1[j1].confidence;
	if(DEBUG && !(isfinite(ConfSumPV1[1]))){
	  printf("j1=%d/%d:SV1[j1].confidence= %0.6e, SV1[j1].sv_overlap=%d,ConfSumPV1[1]=%0.6e\n",
		 j1,svCnt1,SV1[j1].confidence, SV1[j1].sv_overlap ? 1 : 0, ConfSumPV1[1]);
	  fflush(stdout);
	  assert(isfinite(ConfSumPV1[1]));
	}
	totlenPV1 += SV1[j1].refendpos - SV1[j1].refstartpos;
      }
    }
  }
  for(int j2 = 0; j2 < svCnt2; j2++){
    if(SV2[j2].confidence >= LogPvThreshold){
      unique2 += SV2[j2].sv_overlap ? 0 : 1;
      if(SVXmapOverlap > 0.0 && !SV2[j2].sv_overlap && SV2[j2].have_ref_repeat)
	uniqueFP2++;
      ConfSum2[SV2[j2].sv_overlap ? 1 : 0] += SV2[j2].confidence;
      if(DEBUG && !(isfinite(ConfSum2[1]))){
	printf("j2=%d/%d:SV2[j2].confidence= %0.6e, SV2[j2].sv_overlap=%d,ConfSum2[1]=%0.6e\n",
	       j2,svCnt2,SV2[j2].confidence, SV2[j2].sv_overlap ? 1 : 0, ConfSum2[1]);
	fflush(stdout);
	assert(isfinite(ConfSum2[1]));
      }
      totlen2 += SV2[j2].refendpos - SV2[j2].refstartpos;

      if(SV2[j2].confidence >= SVConf){
	uniquePV2 += SV2[j2].sv_overlap ? 0 : 1;
	totPV2++;
	if(SIZEBINS > 0){
	  if(DEBUG && !(SV2[j2].sbin != 0 && abs(SV2[j2].sbin) <= SIZEBINS)){
	    printf("j2=%d/%d:SV2[j2].sbin=%d,SIZEBINS=%d\n",j2,svCnt2,SV2[j2].sbin,SIZEBINS);
	    fflush(stdout);
	    assert(SV2[j2].sbin != 0 && abs(SV2[j2].sbin) <= SIZEBINS);
	  }
	  int bin = SV2[j2].sbin + SIZEBINS;
	  double Rlen = SV2[j2].refendpos - SV2[j2].refstartpos;
	  totPV2siz[bin]++;
	  totRlenPV2siz[bin] += Rlen;
	  if(!SV2[j2].sv_overlap){
	    uniquePV2siz[bin]++;
	    uniqueRlenPV2siz[bin] += Rlen;
	  }
	}

	if(SVXmapOverlap > 0.0 && !SV2[j2].sv_overlap && SV2[j2].have_ref_repeat)
	  uniquePVFP2++;
	ConfSumPV2[SV2[j2].sv_overlap ? 1 : 0] += SV2[j2].confidence;
	if(DEBUG && !(isfinite(ConfSumPV2[1]))){
	  printf("j2=%d/%d:SV2[j2].confidence= %0.6e, SV2[j2].sv_overlap=%d,ConfSumPV2[1]=%0.6e\n",
		 j2,svCnt2,SV2[j2].confidence, SV2[j2].sv_overlap ? 1 : 0, ConfSumPV2[1]);
	  fflush(stdout);
	  assert(isfinite(ConfSumPV2[1]));
	}
	totlenPV2 += SV2[j2].refendpos - SV2[j2].refstartpos;
      }
    }
  }

  printf("Indels with confidence above %0.2f:\n",LogPvThreshold);
  if(SVXmapOverlap > 0.0){
    printf("  %s:indels=%d (%0.3f kb), overlapped=%d (%0.3f kb, avg conf=%0.2f), unique=%d(avg conf=%0.2f),confirmed=%d\n",
	   filename1, svCnt1-tcnt1, totlen1, svCnt1-tcnt1-unique1, matchlen, ConfSum1[1]/max(1,svCnt1-tcnt1-unique1), unique1, ConfSum1[0]/max(1,unique1),uniqueFP1);
    printf("  %s:indels=%d (%0.3f kb), overlapped=%d (%0.3f kb, avg conf=%0.2f), unique=%d(avg conf=%0.2f),confirmed=%d\n",
	   filename2, svCnt2-tcnt2, totlen2, svCnt2-tcnt2-unique2, matchlen, ConfSum2[1]/max(1,svCnt2-tcnt2-unique2), unique2, ConfSum2[0]/max(1,unique2),uniqueFP2);
  } else {
    printf("  %s:indels=%d (%0.3f kb), overlapped=%d (%0.3f kb, avg conf=%0.2f), unique=%d(avg conf=%0.2f)\n",
	   filename1, svCnt1-tcnt1, totlen1, svCnt1-tcnt1-unique1, matchlen, ConfSum1[1]/max(1,svCnt1-tcnt1-unique1), unique1, ConfSum1[0]/max(1,unique1));
    printf("  %s:indels=%d (%0.3f kb), overlapped=%d (%0.3f kb, avg conf=%0.2f), unique=%d(avg conf=%0.2f)\n",
	   filename2, svCnt2-tcnt2, totlen2, svCnt2-tcnt2-unique2, matchlen, ConfSum2[1]/max(1,svCnt2-tcnt2-unique2), unique2, ConfSum2[0]/max(1,unique2));
  }
  fflush(stdout);
  if(DEBUG) assert(isfinite(ConfSum1[1]));
  if(DEBUG) assert(isfinite(ConfSum2[1]));

  FILE *SVfp = NULL;
  if(SVfile){
    if((SVfp = fopen(SVfile,"a"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("failed to open file %s for append of indel statistics:errno=%d:%s\n",SVfile,eno,err);
      exit(1);
    }
    if(VERB/* HERE HERE >=2 */){
      printf("Appending statistics to %s\n",SVfile);
      fflush(stdout);
    }
  }

  if(SVConf > LogPvThreshold){
    printf("Indels with confidence above %0.2f:\n",SVConf);
    if(SVXmapOverlap > 0.0){
      printf("  %s:indels=%d (%0.3f kb), overlapped=%d (avg conf=%0.2f), unique=%d(avg conf=%0.2f),confirmed=%d\n",
	     filename1, svCnt1-tcntPV1, totlenPV1, svCnt1-tcntPV1-uniquePV1, ConfSumPV1[1]/max(1,svCnt1-tcntPV1-uniquePV1), uniquePV1, ConfSumPV1[0]/max(1,uniquePV1),uniquePVFP1);
      printf("  %s:indels=%d (%0.3f kb), overlapped=%d (avg conf=%0.2f), unique=%d(avg conf=%0.2f),confirmed=%d\n",
	     filename2, svCnt2-tcntPV2, totlenPV2, svCnt2-tcntPV2-uniquePV2, ConfSumPV2[1]/max(1,svCnt2-tcntPV2-uniquePV2), uniquePV2, ConfSumPV2[0]/max(1,uniquePV2),uniquePVFP2);

      if(DEBUG && !(totPV1 == svCnt1-tcntPV1)){
	printf("WARNING:totPV1=%d,svCnt1=%d,tcntPV1=%d: expected totPV1 == svCnt1-tcntPV1 (err=%d)\n",
	       totPV1,svCnt1,tcntPV1,totPV1 - (svCnt1-tcntPV1));
	fflush(stdout);
      }
      if(DEBUG && !(totPV2 == svCnt2-tcntPV2)){
	printf("WARNING:totPV2=%d,svCnt2=%d,tcntPV2=%d: expected totPV2 == svCnt2-tcntPV2 (err=%d)\n",
	       totPV2,svCnt2,tcntPV2,totPV2 - (svCnt2-tcntPV2));
	fflush(stdout);
      }

      printf("# SVConf=%0.2f TotCnt1=%d,Uniq1=%0.2f%%,ExclEnds1=%0.2f%% TotCnt2=%d,Uniq2=%0.2f%%,ExclEnds2=%0.2f%%\n",
	     SVConf, svCnt1-tcntPV1,uniquePV1*100.0/((double)(svCnt1-tcntPV1)),uniquePVFP1*100.0/((double)(svCnt1-tcntPV1)),
	     svCnt2-tcntPV2, uniquePV2*100.0/((double)(svCnt2-tcntPV2)), uniquePVFP2*100.0/((double)(svCnt2-tcntPV2)));
      if(SIZEBINS > 0){
	printf("#h Size\tDelCnt1(RefL)\tDelUniq1(RefL)\tDelCnt2(RefL)\tDelUniq2(RefL)\tInsCnt1(RefL)\tInsUniq1(RefL)\tInsCnt2(RefL)\tInsUniq2(RefL)\n");
	int totcnt1 = 0, totcnt2 = 0;
	int uniqcnt1 = 0, uniqcnt2 = 0;
	for(int i = 1; i <= SIZEBINS; i++){
	  printf("%0.1f",SizeValues[i-1]);
	  int sbin = -i;
	  int bin = sbin + SIZEBINS;
	  if(totPV1siz[bin] <= 0)
	    printf("\t%3d(%4.1f)\t%6.2f%%(%4.1f)", 0, -1.0, 100.0, -1.0);
	  else
	    printf("\t%3d(%4.1f)\t%6.2f%%(%4.1f)", totPV1siz[bin], totRlenPV1siz[bin]/totPV1siz[bin], 100.0*uniquePV1siz[bin]/(double)totPV1siz[bin], uniquePV1siz[bin] ? uniqueRlenPV1siz[bin]/uniquePV1siz[bin] : -1.0);

	  if(totPV2siz[bin] <= 0)
	    printf("\t%3d(%4.1f)\t%6.2f%%(%4.1f)", 0, -1.0, 100.0, -1.0);
	  else
	    printf("\t%3d(%4.1f)\t%6.2f%%(%4.1f)", totPV2siz[bin], totRlenPV2siz[bin]/totPV2siz[bin], 100.0*uniquePV2siz[bin]/(double)totPV2siz[bin], uniquePV2siz[bin] ? uniqueRlenPV2siz[bin]/uniquePV2siz[bin] : -1.0);

	  totcnt1 += totPV1siz[bin];
	  totcnt2 += totPV2siz[bin];
	  uniqcnt1 += uniquePV1siz[bin];
	  uniqcnt2 += uniquePV2siz[bin];

	  sbin = i;
	  bin = sbin + SIZEBINS;

	  if(totPV1siz[bin] <= 0)
	    printf("\t%3d(%4.1f)\t%6.2f%%(%4.1f)", 0, -1.0, 100.0, -1.0);
	  else
	    printf("\t%3d(%4.1f)\t%6.2f%%(%4.1f)", totPV1siz[bin], totRlenPV1siz[bin]/totPV1siz[bin], 100.0*uniquePV1siz[bin]/(double)totPV1siz[bin], uniquePV1siz[bin] ? uniqueRlenPV1siz[bin]/uniquePV1siz[bin] : -1.0);

	  if(totPV2siz[bin] <= 0)
	    printf("\t%3d(%4.1f)\t%6.2f%%(%4.1f)", 0, -1.0, 100.0, -1.0);
	  else
	    printf("\t%3d(%4.1f)\t%6.2f%%(%4.1f)", totPV2siz[bin], totRlenPV2siz[bin]/totPV2siz[bin], 100.0*uniquePV2siz[bin]/(double)totPV2siz[bin], uniquePV2siz[bin] ? uniqueRlenPV2siz[bin]/uniquePV2siz[bin] : -1.0);

	  printf("\n");

	  totcnt1 += totPV1siz[bin];
	  totcnt2 += totPV2siz[bin];
	  uniqcnt1 += uniquePV1siz[bin];
	  uniqcnt2 += uniquePV2siz[bin];
	}

	if(DEBUG && !(totcnt1 == totPV1)){
	  printf("WARNING:totcnt1=%d,totPV1=%d (should match)\n", totcnt1,totPV1);
	  fflush(stdout);
	}
	if(DEBUG && !(totcnt2 == totPV2)){
	  printf("WARNING:totcnt2=%d,totPV2=%d (should match)\n", totcnt2,totPV2);
	  fflush(stdout);
	}
	if(DEBUG && !(uniqcnt1 == uniquePV1)){
	  printf("WARNING:uniqcnt1=%d,uniquePV1=%d (should match)\n", uniqcnt1,uniquePV1);
	  fflush(stdout);
	}
	if(DEBUG && !(uniqcnt2 == uniquePV2)){
	  printf("WARNING:uniqcnt2=%d,uniquePV2=%d (should match) : totcnt2=%d,totPV2=%d\n", uniqcnt2,uniquePV2, totcnt2, totPV2);
	  fflush(stdout);
	}
      }
      printf("\n");
      fflush(stdout);

      if(SVfp){
	printversion(SVfp);
	fprintf(SVfp, "# SVConf=%0.2f TotCnt1=%d,Uniq1=%0.2f%%,ExclEnds1=%0.2f%% TotCnt2=%d,Uniq2=%0.2f%%,ExclEnds2=%0.2f%%\n",
		SVConf, svCnt1-tcntPV1,uniquePV1*100.0/((double)(svCnt1-tcntPV1)),uniquePVFP1*100.0/((double)(svCnt1-tcntPV1)),
		svCnt2-tcntPV2, uniquePV2*100.0/((double)(svCnt2-tcntPV2)), uniquePVFP2*100.0/((double)(svCnt2-tcntPV2)));
	if(SIZEBINS > 0){
	  fprintf(SVfp,"#h Size\tDelCnt1(RefL)\tDelUniq1(RefL)\tDelCnt2(RefL)\tDelUniq2(RefL)\tInsCnt1(RefL)\tInsUniq1(RefL)\tInsCnt2(RefL)\tInsUniq2(RefL)\n");
	  for(int i = 1; i <= SIZEBINS; i++){
	    fprintf(SVfp,"%0.1f",SizeValues[i-1]);
	    int sbin = -i;
	    int bin = sbin + SIZEBINS;
	    if(totPV1siz[bin] <= 0)
	      fprintf(SVfp,"\t%3d(%4.1f)\t%6.2f%%(%4.1f)", 0, -1.0, 100.0, -1.0);
	    else
	      fprintf(SVfp,"\t%3d(%4.1f)\t%6.2f%%(%4.1f)", totPV1siz[bin], totRlenPV1siz[bin]/totPV1siz[bin], 100.0*uniquePV1siz[bin]/(double)totPV1siz[bin], uniquePV1siz[bin] ? uniqueRlenPV1siz[bin]/uniquePV1siz[bin] : -1.0);

	    if(totPV2siz[bin] <= 0)
	      fprintf(SVfp,"\t%3d(%4.1f)\t%6.2f%%(%4.1f)", 0, -1.0, 100.0, -1.0);
	    else
	      fprintf(SVfp,"\t%3d(%4.1f)\t%6.2f%%(%4.1f)", totPV2siz[bin], totRlenPV2siz[bin]/totPV2siz[bin], 100.0*uniquePV2siz[bin]/(double)totPV2siz[bin], uniquePV2siz[bin] ? uniqueRlenPV2siz[bin]/uniquePV2siz[bin] : -1.0);


	    sbin = i;
	    bin = sbin + SIZEBINS;

	    if(totPV1siz[bin] <= 0)
	      fprintf(SVfp,"\t%3d(%4.1f)\t%6.2f%%(%4.1f)", 0, -1.0, 100.0, -1.0);
	    else
	      fprintf(SVfp,"\t%3d(%4.1f)\t%6.2f%%(%4.1f)", totPV1siz[bin], totRlenPV1siz[bin]/totPV1siz[bin], 100.0*uniquePV1siz[bin]/(double)totPV1siz[bin], uniquePV1siz[bin] ? uniqueRlenPV1siz[bin]/uniquePV1siz[bin] : -1.0);

	    if(totPV2siz[bin] <= 0)
	      fprintf(SVfp,"\t%3d(%4.1f)\t%6.2f%%(%4.1f)", 0, -1.0, 100.0, -1.0);
	    else
	      fprintf(SVfp,"\t%3d(%4.1f)\t%6.2f%%(%4.1f)", totPV2siz[bin], totRlenPV2siz[bin]/totPV2siz[bin], 100.0*uniquePV2siz[bin]/(double)totPV2siz[bin], uniquePV2siz[bin] ? uniqueRlenPV2siz[bin]/uniquePV2siz[bin] : -1.0);

	    fprintf(SVfp,"\n");
	  }
	  fflush(SVfp);
	}
      }
    } else {
      printf("  %s:indels=%d (%0.3f kb), overlapped=%d (avg conf=%0.2f), unique=%d(avg conf=%0.2f)\n",
	     filename1, svCnt1-tcntPV1, totlenPV1, svCnt1-tcntPV1-uniquePV1, ConfSumPV1[1]/max(1,svCnt1-tcntPV1-uniquePV1), uniquePV1, ConfSumPV1[0]/max(1,uniquePV1));
      printf("  %s:indels=%d (%0.3f kb), overlapped=%d (avg conf=%0.2f), unique=%d(avg conf=%0.2f)\n",
	     filename2, svCnt2-tcntPV2, totlenPV2, svCnt2-tcntPV2-uniquePV2, ConfSumPV2[1]/max(1,svCnt2-tcntPV2-uniquePV2), uniquePV2, ConfSumPV2[0]/max(1,uniquePV2));
      printf("# SVConf=%0.2f TotCnt=%0.1f Uniq%%=%0.2f\n",SVConf, ((double)(svCnt1-tcntPV1+svCnt2-tcntPV2))/2, (uniquePV1+uniquePV2)*100.0/((double)(svCnt1-tcntPV1+svCnt2-tcntPV2))); 
      if(SVfp){
	fprintf(SVfp,"# SVConf=%0.2f TotCnt=%0.1f Uniq%%=%0.2f\n",SVConf, ((double)(svCnt1-tcntPV1+svCnt2-tcntPV2))/2, (uniquePV1+uniquePV2)*100.0/((double)(svCnt1-tcntPV1+svCnt2-tcntPV2))); 
	fflush(SVfp);
      }
    }
    fflush(stdout);

    if(DEBUG) assert(isfinite(ConfSumPV1[1]));
    if(DEBUG) assert(isfinite(ConfSumPV2[1]));
  }

  if(SVfp){
    (void)fclose(SVfp);
    SVfp = NULL;
  }

  if(VERB/* HERE >=2 */){
    printf("\nUnique SVs from %s with confidence above %0.2f:\n",filename1,max(LogPvThreshold,SVConf));
    for(int j1 = 0; j1 < svCnt1; j1++){
      if(SV1[j1].sv_overlap == false && SV1[j1].confidence >= max(LogPvThreshold,SVConf)){
	if(SVXmapOverlap > 0.0 && SV1[j1].have_ref_repeat)
	  printf("queryid=%lld refid=%lld refstart=%0.3f,refend=%0.3f,type=%s,confidence=%0.2f: confirmed\n", SV1[j1].qrycontigid,SV1[j1].refcontigid,SV1[j1].refstartpos,SV1[j1].refendpos,SV1[j1].type,SV1[j1].confidence);
	else
	  printf("queryid=%lld refid=%lld refstart=%0.3f,refend=%0.3f,type=%s,confidence=%0.2f\n", SV1[j1].qrycontigid,SV1[j1].refcontigid,SV1[j1].refstartpos,SV1[j1].refendpos,SV1[j1].type,SV1[j1].confidence);
      }
    }

    printf("\nUnique SVs from %s with confidence above %0.2f:\n",filename2,max(LogPvThreshold,SVConf));
    for(int j2 = 0; j2 < svCnt2; j2++){
      if(SV2[j2].sv_overlap == false && SV2[j2].confidence >= max(LogPvThreshold,SVConf)){
	if(SVXmapOverlap > 0.0 && SV2[j2].have_ref_repeat)
	  printf("queryid=%lld refid=%lld refstart=%0.3f,refend=%0.3f,type=%s,confidence=%0.2f: confirmed\n", SV2[j2].qrycontigid,SV2[j2].refcontigid,SV2[j2].refstartpos,SV2[j2].refendpos,SV2[j2].type,SV2[j2].confidence);
	else
	  printf("queryid=%lld refid=%lld refstart=%0.3f,refend=%0.3f,type=%s,confidence=%0.2f\n", SV2[j2].qrycontigid,SV2[j2].refcontigid,SV2[j2].refstartpos,SV2[j2].refendpos,SV2[j2].type,SV2[j2].confidence);
      }
    }
  }

  fflush(stdout);
}

void input_svcheck(char *filename, double svcheckKB, int &numrefmaps, int &maxrefmaps, Cmap **&refmap)
{
  svMax = svCnt = 0;
  SVtocheck = 0;
  input_indel(filename, SVtocheck, svCnt, svMax, 1);/* read in the SVs and append query & ref maps  to map[0..nummaps-1] */
  if(VERB){
    printf("Read %d Indels from %s\n",svCnt,filename);
    fflush(stdout);
  }

  if(smap_min_conf > 0.0){/* filter output SVs below -svMinConf level */
    int cnt = 0;
    for(int S = 0; S < svCnt; S++){
      if(SVtocheck[S].confidence < smap_min_conf)
	continue;
      if(cnt < S)
	SVtocheck[cnt] = SVtocheck[S];
      cnt++;
    }
    if(VERB){
      printf("Filtered out %d of %d SVs due to -svMinConf %0.3f : %d SVs remaining\n",svCnt - cnt, svCnt, smap_min_conf, cnt);
      fflush(stdout);
    }
    svCnt = cnt;
  }

  svCntOrig = svCnt;

  /* sort SVs in order of id1, id2, sites1[0] so mergeable SVs are next to each other */
  qsort(SVtocheck,svCnt,sizeof(xmapEntry), (intcmp*)Id1Id2Site1Inc);  
  for(int S = 0; S < svCntOrig; S++){
    xmapEntry *pSV = &SVtocheck[S];
    pSV->smapentryid = S+1;
    pSV->mergecnt = 1;
    pSV->mergestart = S;
  }
  for(int S = 0; S < svCntOrig; S++){
    xmapEntry *pSV = &SVtocheck[S];

    /* check if multiple SVs starting at S can be combined */
    for(int R = S+1; R < svCntOrig; R++){
      xmapEntry *qSV = &SVtocheck[R-1];
      xmapEntry *rSV = &SVtocheck[R];
      if(qSV->xmapid != rSV->xmapid)
	break;
      if(DEBUG/* HERE >=2 */) assert(qSV->qrycontigid == rSV->qrycontigid);
      if(DEBUG/* HERE >=2 */) assert(qSV->refcontigid == rSV->refcontigid);
      //      if(DEBUG/* HERE >=2 */) assert(qSV->refendpos <= rSV->refstartpos + 1e-6);// May not be true due to outlier region including the mis-resolved regions at both ends
      if(DEBUG/* HERE >=2 */) assert(qSV->orientforward == rSV->orientforward);
      if(qSV->refendpos < rSV->refstartpos - svcheckMerge)
	break;

      /* create a merged SV combining SVtocheck[S..R] */
      if(svCnt >= svMax){
	svMax *= 2;
	xmapEntry *newSVs = new xmapEntry[svMax];
	memcpy(newSVs,SVtocheck,sizeof(xmapEntry)*svCnt);
	delete [] SVtocheck;
	SVtocheck = newSVs;
	pSV = &SVtocheck[S];
	qSV = &SVtocheck[R-1];
	rSV = &SVtocheck[R];
      }
      xmapEntry *tSV = &SVtocheck[svCnt++];
      *tSV = *pSV;/* copy most fields from SVtocheck[S] */
      tSV->smapentryid = svCnt;
      tSV->mergestart = S;
      tSV->mergecnt = R-S+1;

      tSV->refendpos = rSV->refendpos;
      if(pSV->orientforward){
	tSV->qrystartpos = pSV->qrystartpos;
	tSV->qryendpos = rSV->qryendpos;
	if(DEBUG && !(tSV->qrystartpos < tSV->qryendpos - 1e-6)){
	  printf("While merging SVtocheck[S..R],S=%d,R=%d:SVtocheck[S].qrystartpos=%0.7f,SVtocheck[R].qryendpos=%0.7f\n",
		 S,R,tSV->qrystartpos,tSV->qryendpos);
	  for(int t = S; t <= R; t++){
	    xmapEntry *nSV = &SVtocheck[t];
	    printf("\t SVtocheck[%d]: xmapid=%d,qrycontigid=%lld,refcontigid=%lld,or=%d:qrystartpos=%0.6f,qryendpos=%0.6f,refstartpos=%0.6f,refendpos=%0.6f\n",
		   t,nSV->xmapid,nSV->qrycontigid,nSV->refcontigid,nSV->orientforward ? 0 : 1,nSV->qrystartpos,nSV->qryendpos,nSV->refstartpos,nSV->refendpos);
	  }
	  fflush(stdout);
	  assert(tSV->qrystartpos < tSV->qryendpos - 1e-6);
	}
      } else {
	tSV->qrystartpos = pSV->qrystartpos;
	tSV->qryendpos = rSV->qryendpos;
	if(DEBUG && !(tSV->qryendpos < tSV->qrystartpos - 1e-6)){
	  printf("SV[%d]:qrystartpos=%0.3f,qryendpos=%0.3f\n", svCnt-1,tSV->qrystartpos,tSV->qryendpos);
	  printf("SV[%d]:qrystartpos=%0.3f,qryendpos=%0.3f\n", S,pSV->qrystartpos,pSV->qryendpos);
	  printf("SV[%d]:qrystartpos=%0.3f,qryendpos=%0.3f\n", R,rSV->qrystartpos,rSV->qryendpos);
	  fflush(stdout);

	  assert(tSV->qryendpos < tSV->qrystartpos - 1e-6);
	}
      }
      for(int t = S+1; t <= R; t++)
	tSV->confidence += SVtocheck[t].confidence;
      double reflen = tSV->refendpos - tSV->refstartpos;
      double qrylen = fabs(tSV->qryendpos - tSV->qrystartpos);
      char buf[BUFSIZ];
      sprintf(buf,"%s_%d",reflen > qrylen ? "deletion" : "insertion", (int)floor(1000.0*fabs(qrylen-reflen)+0.5));
      tSV->type = strdup(buf);

      if(VERB/* HERE >=2 */){
	printf("Merging SV[%d..%d] to create SV[%d]:qid=%lld,rid=%lld:QueryRange=%0.3f..%0.3f(%0.3f),RefRange=%0.3f..%0.3f(%0.3f),or=%d,type=%s,confidence=%0.2f\n",
	       S,R,svCnt-1,tSV->qrycontigid,tSV->refcontigid, tSV->qrystartpos, tSV->qryendpos, tSV->QryMap->site[0][tSV->QryMap->numsite[0]+1],
	       tSV->refstartpos,tSV->refendpos,tSV->RefMap->site[0][tSV->RefMap->numsite[0]+1],tSV->orientforward ? 0 : 1, tSV->type, tSV->confidence);
	if(VERB>=3){
	  for(int t = S; t <= R; t++){
	    xmapEntry *pSV = &SVtocheck[t];
	    printf("\t SV[%d/%d]:qid=%lld,rid=%lld:QueryRange=%0.3f..%0.3f(%0.3f),RefRange=%0.3f..%0.3f(%0.3f),or=%d,type=%s,confidence=%0.2f\n",
		   t, svCntOrig, pSV->qrycontigid,pSV->refcontigid,pSV->qrystartpos,pSV->qryendpos, pSV->QryMap->site[0][pSV->QryMap->numsite[0]+1],
		   pSV->refstartpos,pSV->refendpos,pSV->RefMap->site[0][pSV->RefMap->numsite[0]+1], pSV->orientforward ? 0 : 1, pSV->type, pSV->confidence);
	  }
	}
	fflush(stdout);
      }
    } /* R = S+1 .. svCntOrig-1 */
  } /* S = 0 .. svCntOrig-1 */

  for(int S = 0; S < svCnt; S++){
    xmapEntry *pSV = &SVtocheck[S];
    if(DEBUG && !(pSV->smapentryid == S+1)){
      printf("SV[%d]:smapentryid=%d (expected %d)\n",S,pSV->smapentryid,S+1);
      fflush(stdout);
      assert(pSV->smapentryid == S+1);
    }
    Cmap *QryMap = pSV->QryMap;
    Cmap *RefMap = pSV->RefMap;

    if(VERB/* HERE >=2 */){
      printf("SV[%d]:qid=%lld,rid=%lld:QueryRange=%0.3f..%0.3f,RefRange=%0.3f..%0.3f,or=%d,type=%s,conf=%0.2f:Testing with contig id=%d and %d\n",
	     S,pSV->qrycontigid,pSV->refcontigid, pSV->qrystartpos, pSV->qryendpos, pSV->refstartpos, pSV->refendpos, pSV->orientforward ? 0 : 1, pSV->type, pSV->confidence, numrefmaps+1, numrefmaps+2);
      fflush(stdout);
    }

    /* generate reference map pair from SV */
    maxmapalloc(numrefmaps+2,maxrefmaps,refmap,0,1);

    mapsplice(refmap[numrefmaps],(long long)(numrefmaps+1),pSV,svcheckKB,0,QryMap,RefMap);
    refmap[numrefmaps++]->startloc = pSV->refstartpos + RefMap->id * CHR_OFFSET * 0.001;

    mapsplice(refmap[numrefmaps],(long long)(numrefmaps+1),pSV,svcheckKB,1,QryMap,RefMap);
    refmap[numrefmaps++]->startloc = pSV->refstartpos + RefMap->id * CHR_OFFSET * 0.001;

    if(DEBUG) assert(numrefmaps > 0);// check for 31 bit overflow
  }

  if(VERB){
    printf("Read %d indels from %s and added %d Merged indels\n",svCntOrig,filename,svCnt-svCntOrig);
    printf("Created %d reference maps (one pair per indel)\n",numrefmaps);
    fflush(stdout);
  }

  nummaps = startmaps;/* extra Query and Ref maps are still at maps[startmaps ... ] */
  num_files -= 2;
}
