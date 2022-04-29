#ifndef STRUCTURALVARIATION_H
#define STRUCTURALVARIATION_H

#include "globals.h"
#include "parameters.h"
//#include "xmapEntry.h"
#include "conf_params.h"

#define LATE_ZYGOSITY 1 // compute zygosity only after final xmap list is computed

#define SV_SIZE smapSize /* add SVsize column to .smap output and change smap version to 0.8 */
#define SV_FREQ smapFreq /* add SVsize AND SVfreq columns to .smap output and change smap version to 0.9 */
#define SV_ORIENT smapTransOrient /* add SVsize, SVfreq AND orientation colums to .smap output and change smap version to 0.9 */

#define SMALL_CONFIRM_FIX 1 /* modify -svSmallDelConfirm L to expand ref interval until a label is found that is NOT within L of next label on either side */

//used in simpleRepeat.cpp and output_{x,s}map.cpp
template<class T> 
int growArray(T*& ina, const int oldsize, int newsize=0) {
  //printf("growArray start: oldsize %i newsize %i\n", oldsize, newsize); fflush(stdout); //debug
  T* tmp = ina;
  //if newsize > oldsize, take newsize, otherwise, try to get oldsize*2, but be sure you get at least one ele
  newsize = (newsize > oldsize ? newsize : max(oldsize+1, oldsize*2) );
  assert(newsize > 0);
  ina = new T[newsize];
  for(int i=0; i<oldsize; i++) //copy eles of old into new
    ina[i] = tmp[i];
  if( oldsize ) //if no elements, can't delete
    delete [] tmp; //tmp is old array--container no longer necessary
  //printf("growArray end  : oldsize %i newsize %i\n", oldsize, newsize); //debug
  return newsize;
}


extern const char* transintra_name          ;
extern const char* transinter_name          ;
extern const char* transintra_overlap_name  ;
extern const char* transinter_overlap_name  ;
extern const char* transintra_repeat_name   ;
extern const char* transinter_repeat_name   ;
extern const char* transintra_duplicate_name;
extern const char* transinter_duplicate_name;
extern const char* transintra_common_name;
extern const char* transinter_common_name;
extern const char* transintra_segdup_name;
extern const char* transinter_segdup_name;
extern const bool smap_xmapid_col; //include the xmap entry ids in smap; could be moved to parameters.cpp if promote to command line arg
extern const bool smap_linkid_col; //add column after XmapIDs for LinkID
extern const bool smap_labelind_col; //add four columns for label indices after XmapIDs and LinkID
extern const bool smap_xmap_allpairs; //don't just take neighboring xmap entries, take all (unique) pairs
extern const float sv_overlap_frac; //fraction two svs must overlap to be called 'the same'
//extern int smap_genotype_ngroup; //count number of genotypeGroups -- no need to be global
//extern const double translocation_overlap_buf; //buffer around translocation breapoint for overlap (alignmentOverlapTranslocation and svOverlap only, NOT bed overlap, so doesn't need to be global)

//write all data to first argument--a file handle
extern void writeToSmap(FILE *fp, int nsv, long long qrycontigid, long long refcontigid1, long long refcontigid2, double qrystart, double qrystop, double refstart, double refstop, const char *svtype, int xmapid1, int xmapid2, int linkid, int qsi, int qei, int rsi, int rei, double conf, double confl, double confr, double confc, const char* zyg, int gen, int ggrp, float confs, double SVsize, double qrycov, double refcov, double SVfreq, char *orientation);

//enum for algo: currently only pairsplit -sv and Thomas' -indel
enum sv_algo_t {
  sv_unknown = -1,
  sv_ps = 0, //pairsplit : really means that indices[2] are index into usexmapentries[] (rather than xmapentries[]) and includes Indels with 2 different matchgroups
  sv_indel = 1 // really means indices[2] are index in xmapentries[] and includes Indels with 1 single matchgroup and small Duplications & Inversions with 2 different matchgroups
};

//for type_enum data member of structuralVariation
enum sv_type_t {
  unknown_type = 0,
  insertion = 1,
  deletion = 2,
  intrachr_translocation = 3,
  interchr_translocation = 4,
  inversion_type = 5,
  inversion_small_type = 6,// will be paired with another sv_array[] entry (also happens to have algo == sv_indel) : (NOTE : in rare cases the paired inversion will be filtered out and this type changed to inversion_type)
  duplication_type = 7, // some of these could be Small Duplications : in that case algo == sv_indel
  duplicationinv_type = 8,// some of these could be Small Inverted Duplications : in that case algo == sv_indel
  duplicationsplit_type = 9,
  compression = 10, // low confidence large deletion with Qry overlap larger than smap_del_bigqryoverlap
  complex_type = 11 
};


extern double makeConfidence(double a, double b, double c);
extern int compare_svs(const void *a, const void* b); //for qsort'ing structuralVariation objects
extern bool overlapPoint(double min1, double max1, double p);
extern double overlapLength(double min1, double max1, double min2, double max2);
extern void alignmentOverlapTranslocation(xmapEntry **usexmapentries, int numxmapent);
extern void alignmentOverlapInversion(xmapEntry **usexmapentries, int numxmapent);
extern void alignmentOverlap(xmapEntry **xmapent, int numxmapent); //check alignment overlapping SV
extern void svOverlap(xmapEntry **xmapent, int numxmapent); //check overlapping SVs (indels)
//extern void doTransConf(xmapEntry **usexmapentries);
extern void doTransInvConf(xmapEntry **usexmapentries, bool dotrans);
extern void doPairInversionBP(xmapEntry **usexmapentries, structuralVariation *sv_array, int numsvIndel, int numsv, bool verbose = false);

//utility fn for overlap btwn two segments min1-max1 and min2-max2 (must be sorted)
//return true if there is overlap, false if none
//if minA == maxB, this is NOT considered overlap
//originally just had doubles, but want ints too: use template
template<class T> 
bool overlap(T min1, T max1, T min2, T max2) {
  return( ( min1 < max2 ) && ( max1 > min2 ) );
}


//store the results of input_mask
class maskSV {
public:
  int chr1, chr2;
  double chr1_start, chr2_start;
  //for now, don't bother storing the SD/range values
  maskSV() {
    chr1 = chr2 = 0;
    chr1_start = chr2_start = 0;
  }
  maskSV(int c1, int c2, double c1s, double c2s) {
    chr1 = c1;
    chr2 = c2;
    chr1_start = c1s;
    chr2_start = c2s;
  }
}; //end class maskSV


extern Cmap **YYmap,**XXmap;

#include "Calign.h"

// Given one ref site idxs and one alignment, find corresponding (absolute) query location (by interpolating, if needed, based on nearest aligned labels) and also return nearest aligned query labels
//   If prefdir == 0 (default) : If interpolating query location, pick nearest label on qry in either direction AND don't adjust ref label or location
//   If prefdir == 1, prefer nearest aligned ref label >= refidx : adjust both ref and qry labels and locations (only applied if aligned labels are present at or on both sides of refidx)
//   If prefdir == -1, prefer nearest aligned ref label <= refidx : adjust both ref and qry labels and locations (only applied if aligned labels are present at or on both sides of refidx)
// NOTE: refloc and qryloc are in bp (NOT kb) 
static void getQry(Calign* p1, int& refidx, double& qryloc, int& qryidx, int prefdir = 0, double *refloc = NULL) {

  Cmap *Ymap = YYmap[p1->mapid1];
  Cmap *Xmap = XXmap[p1->mapid2];
  FLOAT *Y = Ymap->site[0];
  FLOAT *X = Xmap->site[0];
  int N = Ymap->numsite[0];
  int M = Xmap->numsite[0];
  int U1 = p1->numpairs; //n aligned sites in first matchgroup

  if(DEBUG) assert(0 < refidx && refidx <= N);

  int Lqryidx=-1, Rqryidx=-1; //find bounding index(s) in p1->sites2 which matches refidx
  int Lrefidx=-1,Rrefidx=-1;
  for(int t = 0; t < U1; t++) {// NOTE : avoid matches with misresolved labels
    if( p1->sites1[t] <= refidx && (p1->sitesK1[t] <= 0 || t==U1-1 || p1->sites1[t+1] > refidx) ){
      Lrefidx = p1->sites1[t];
      Lqryidx = p1->orientation ? M+1 - p1->sites2[t] : p1->sites2[t];
      if(DEBUG>=2 && !(1 <= Lqryidx && Lqryidx <= M)){
	printf("getQry:Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d:M=%d,Lqryidx=%d,Rqryidx=%d\n",p1->mapid1,Ymap->id,p1->mapid2,Xmap->id,p1->orientation,M,Lqryidx,Rqryidx);
	printf("\t t=%d,U1=%d,p1->sites1[t]=%d,p1->sitesK1[t]=%d,p1->sites2[t]=%d,refidx=%d\n",t,U1,p1->sites1[t],p1->sitesK1[t],p1->sites2[t],refidx);
	fflush(stdout);
	assert(1 <= Lqryidx && Lqryidx <= M);
      }
    }
    if( p1->sites1[t] - p1->sitesK1[t] >= refidx && (p1->sitesK1[t] <= 0 || t==0 || p1->sites1[t-1] - p1->sitesK1[t-1] < refidx) ){
      Rrefidx = p1->sites1[t];
      Rqryidx = p1->orientation ? M+1 - p1->sites2[t] : p1->sites2[t];
      if(DEBUG>=2 && !(1 <= Rqryidx && Rqryidx <= M)){
	printf("getQry:Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d:M=%d,Lqryidx=%d,Rqryidx=%d\n",p1->mapid1,Ymap->id,p1->mapid2,Xmap->id,p1->orientation,M,Lqryidx,Rqryidx);
	printf("\t t=%d,U1=%d,p1->sites1[t]=%d,p1->sitesK1[t]=%d,p1->sites2[t]=%d,refidx=%d\n",t,U1,p1->sites1[t],p1->sitesK1[t],p1->sites2[t],refidx);
	fflush(stdout);
	assert(1 <= Rqryidx && Rqryidx <= M);
      }
      if(p1->sitesK1[t] <= 0)
	break;
    }
  }
  if(DEBUG && !((1 <= Lqryidx && Lqryidx <= M) || (1 <= Rqryidx && Rqryidx <= M))){
    printf("getQry:Yid=%d(id=%lld),Xid=%d(id=%lld),or=%d:M=%d,Lqryidx=%d,Rqryidx=%d,N=%d,refidx=%d\n",p1->mapid1,Ymap->id,p1->mapid2,Xmap->id,p1->orientation,M,Lqryidx,Rqryidx,N,refidx);
    printf("\t U1=%d, I=%d..%d,K=%d..%d,J=%d..%d\n",U1,p1->sites1[0],p1->sites1[U1-1],p1->sitesK1[0],p1->sitesK1[U1-1],p1->sites2[0],p1->sites2[U1-1]);
    for(int t = 0; t < U1; t++)
      printf("t=%d: I=%d,K=%d,J=%d\n",t, p1->sites1[t],p1->sitesK1[t],p1->sites2[t]);
    fflush(stdout);
    assert((1 <= Lqryidx && Lqryidx <= M) || (1 <= Rqryidx && Rqryidx <= M));
  }

  if(1 <= Lqryidx && Lqryidx <= M && 1 <= Rqryidx && Rqryidx <= M){/* 1 <= Lrefidx <= refidx <= Rrefidx <= N */
    if(Lrefidx < Rrefidx){/* interpolate */
      if(prefdir == 0){
	double ratio = (Y[refidx] - Y[Lrefidx]) / (Y[Rrefidx] - Y[Lrefidx]);
	qryloc = X[Lqryidx] + ratio * (X[Rqryidx] - X[Lqryidx]);
	qryidx = fabs(qryloc - X[Lqryidx]) < fabs(qryloc - X[Rqryidx]) ? Lqryidx : Rqryidx;// nearest
	qryloc *= 1e3;
      } else {
	if(prefdir == 1){
	  refidx = Rrefidx;
	  qryidx = Rqryidx;
	} else {
	  refidx = Lrefidx;
	  qryidx = Lqryidx;
	}
	if(refloc)
	  *refloc = Y[refidx] * 1e3;
	qryloc = X[qryidx] * 1e3;
	if(VERB>=2){
	  printf("getQry:prefdir=%d:refidx -> %d, qryidx -> %d, Y[refidx]= %0.4f kb, X[qryidx] = %0.4f kb\n",prefdir,refidx,qryidx,Y[refidx],X[qryidx]);
	  fflush(stdout);
	}
      }
    } else {/* perfect match */
      if(DEBUG) assert(Lrefidx == refidx && refidx == Rrefidx && Lqryidx == Rqryidx);
      if(DEBUG && refloc) assert(fabs(*refloc - Y[refidx] * 1e3) <= 1e-3);
      qryloc = X[Lqryidx] * 1e3;
      qryidx = Lqryidx;
      if(VERB>=2){
	printf("getQry:prefdir=%d:refidx= %d, qryidx -> %d, Y[refidx]= %0.4f kb, X[qryidx] = %0.4f kb\n",prefdir,refidx,qryidx,Y[refidx],X[qryidx]);
	fflush(stdout);
      }
    }
  } else if(1 <= Lqryidx && Lqryidx <= M){/* All aligned labels (except possibly, misresolved labels) were left of refidx : extrapolate from left */
    if(p1->sites1[U1-1] >= refidx && prefdir == -1){
      refidx = Lrefidx;
      qryidx = Lqryidx;
      if(refloc)
	*refloc = Y[refidx] * 1e3;
      qryloc = X[qryidx] * 1e3;
      if(VERB>=2){
	printf("getQry:prefdir=%d:refidx -> %d, qryidx -> %d, Y[refidx]= %0.4f kb, X[qryidx] = %0.4f kb\n",prefdir,refidx,qryidx,Y[refidx],X[qryidx]);
	fflush(stdout);
      }
    } else {/* extrapolate from left */
      double shift = Y[refidx] - Y[Lrefidx];
      qryloc = (X[Lqryidx] + shift) * 1e3;
      qryidx   = Lqryidx;
    }
  } else {/* extrapolate from right */
    if(DEBUG) assert(1 <= Rqryidx && Rqryidx <= M);
    if(p1->sites1[0] - p1->sitesK1[0] <= refidx && prefdir == 1){
      refidx = Rrefidx;
      qryidx = Rqryidx;
      if(refloc)
	*refloc = Y[refidx] * 1e3;
      qryloc = X[qryidx] * 1e3;
      if(VERB>=2){
	printf("getQry:prefdir=%d:refidx -> %d, qryidx -> %d, Y[refidx]= %0.4f kb, X[qryidx] = %0.4f kb\n",prefdir,refidx,qryidx,Y[refidx],X[qryidx]);
	fflush(stdout);
      }
    } else { /* extrapolate from right */
      double shift = Y[Rrefidx] - Y[refidx];
      qryloc = (X[Rqryidx] - shift) * 1e3;
      qryidx = Rqryidx;
      if(VERB>=2){
        printf("getQry:prefdir=%d:refidx= %d,Rrefidx= %d:shift= %0.4f kb, qryidx -> %d, Y[refidx]= %0.4f kb, X[qryidx] = %0.4f kb\n",prefdir,refidx,Rrefidx,shift,qryidx,Y[refidx],X[qryidx]);
        fflush(stdout);
      }
    }
  }
} //end getQry

// Given one qryidx (absolute index on query) and one alignment, find corresponding ref location (by interpolating, if needed, based on nearest aligned labels) and also return nearest aligned ref labels
//   If prefdir == 0 (default) : If interpolating ref location, pick nearest label on ref in either direction AND don't adjust qry label or location
//   If prefdir == 1, prefer nearest aligned qry label >= qryidx : adjust both ref and qry labels and locations (only applied if aligned labels are present on both sides of qryidxx)
//   If prefdir == -1, prefer nearest aligned qry label <= qryidx : adjust both ref and qry labels and locations (only applied if aligned labels are present on both sides of qryidx)
// Returns index into p1->sites1[] & p1->sites2[] that corresponds to ref label refidx returned
// NOTE: refloc and qryloc are in bp (NOT kb) 
static int getRef(Calign* p1, int& qryidx, double& refloc, int& refidx, int prefdir = 0, double *qryloc = NULL, bool verbose = false) {

  Cmap *Ymap = YYmap[p1->mapid1];
  Cmap *Xmap = XXmap[p1->mapid2];
  FLOAT *Y = Ymap->site[0];
  FLOAT *X = Xmap->site[0];
  int N = Ymap->numsite[0];
  int M = Xmap->numsite[0];
  int U1 = p1->numpairs; //n aligned sites in first matchgroup

  if(DEBUG) assert(0 < qryidx && qryidx <= M);

  if(verbose){
    printf("getRef: initial qryidx= %d, M=%d, N=%d, prefdir= %d, orientation= %d\n",qryidx,M,N, prefdir,p1->orientation);
    fflush(stdout);
  }

  int Lqryidx=-1, Rqryidx=-1; //find bounding labels in p1->sites2 which matches qryidx
  int Lrefidx=-1, Rrefidx=-1;
  int Lindex= -1, Rindex= -1;/* index into p1->sites1[] that matches Lrefidx & Rrefidx respectively */

  if(p1->orientation){
    for(int t = U1; --t >= 0;) {
      if( M+1 - p1->sites2[t] <= qryidx){
	Rrefidx = p1->sites1[Rindex = t];
	Lqryidx = M+1 - p1->sites2[t];
      }
      if( M+1 - p1->sites2[t] >= qryidx){
	Lrefidx = p1->sites1[Lindex = t];
	Rqryidx = M+1 - p1->sites2[t];
	break;
      }
    }
  } else {
    for(int t = 0; t < U1; t++) {
      if( p1->sites2[t] <= qryidx){
	Lrefidx = p1->sites1[Lindex = t];
	Lqryidx = p1->sites2[t];
      }
      if( p1->sites2[t] >= qryidx ){
	Rrefidx = p1->sites1[Rindex = t];
	Rqryidx = p1->sites2[t];
	break;
      }
    }
  }
  if(DEBUG && Lqryidx >= 1 && Rqryidx >= 1) assert(Lqryidx <= qryidx && qryidx <= Rqryidx && Rqryidx <= M);
  if(DEBUG) assert((1 <= Lrefidx && Lrefidx <= N) || (1 <= Rrefidx && Rrefidx <= N));

  if(verbose){
    printf("getRef:qryidx=%d (qrypos=%0.3f) : Lqryidx=%d(%0.3f),Rqryidx=%d(%0.3f),Lrefidx=%d,Rrefidx=%d\n",qryidx,X[qryidx],Lqryidx,X[Lqryidx],Rqryidx,X[Rqryidx],Lrefidx,Rrefidx);
    fflush(stdout);
  }

  int ret = -1;

  if(1 <= Lrefidx && Lrefidx <= N && 1 <= Rrefidx && Rrefidx <= N){/* 1 <= Lqryidx <= qryidx <= Rqryidx <= M */
    if(Lqryidx < Rqryidx){/* interpolate */
      if(prefdir == 0){
	double ratio = (X[qryidx] - X[Lqryidx]) / (X[Rqryidx] - X[Lqryidx]);
	refloc = (Y[Lrefidx] + (p1->orientation ? 1.0 - ratio : ratio) * (Y[Rrefidx] - Y[Lrefidx]));
	refidx = fabs(refloc - Y[Lrefidx]) < fabs(refloc - Y[Rrefidx]) ? Lrefidx : Rrefidx;// nearest
	ret = fabs(refloc - Y[Lrefidx]) < fabs(refloc - Y[Rrefidx]) ? Lindex : Rindex;
	refloc *= 1e3;
      } else {
	if(prefdir == 1){
	  qryidx = Rqryidx;
	  refidx = (p1->orientation ? Lrefidx : Rrefidx);
	  ret = (p1->orientation ? Lindex : Rindex);
	} else {
	  qryidx = Lqryidx;
	  refidx = (p1->orientation ? Rrefidx : Lrefidx);
	  ret = (p1->orientation ? Rindex : Lindex);
	}
	refloc = Y[refidx] * 1e3;
	if(qryloc){
	  *qryloc = X[qryidx] * 1e3;
	  if(verbose){
	    printf("qryidx -> %d, qryloc -> %0.3f\n", qryidx, *qryloc * 1e-3);
	    fflush(stdout);
	  }
	}
      }
    } else {/* exact match */
      if(DEBUG && !(Lqryidx == qryidx && qryidx == Rqryidx && Lrefidx == Rrefidx)){
	printf("getRef: refid=%d(id=%lld),mapid=%d(id=%lld),or=%d,M=%d,N=%d:qryidx=%d,Lrefidx=%d,Rrefidx=%d,Lqryidx=%d,Rqryidx=%d\n",
	       p1->mapid1,Ymap->id,p1->mapid2,Xmap->id,p1->orientation,M,N,qryidx,Lrefidx,Rrefidx,Lqryidx,Rqryidx);
	fflush(stdout);
	assert(Lqryidx == qryidx && qryidx == Rqryidx && Lrefidx == Rrefidx);
      }
      refloc = Y[Lrefidx] * 1e3;
      refidx = Lrefidx;
      ret = Lindex;
      if(qryloc){
	*qryloc = X[qryidx] * 1e3;
	if(verbose){
	  printf("qryidx -> %d, qryloc -> %0.3f\n", qryidx, *qryloc * 1e-3);
	  fflush(stdout);
	}
      }
    }
  } else if(1 <= Lrefidx && Lrefidx <= N){/* extrapolate from left */
    if(p1->orientation){
      double shift = X[Rqryidx] - X[qryidx];
      refloc = (Y[Lrefidx] + shift) * 1e3;
    } else {
      double shift = X[qryidx] - X[Lqryidx];
      refloc = (Y[Lrefidx] + shift) * 1e3;
    }
    refidx = Lrefidx;
    ret = Lindex;
    if(qryloc){
      *qryloc = X[qryidx] * 1e3;
      if(verbose){
	printf("qryidx -> %d, qryloc -> %0.3f\n", qryidx, *qryloc * 1e-3);
	fflush(stdout);
      }
    }
  } else {/* extrapolate from right */
    if(DEBUG) assert(1 <= Rrefidx && Rrefidx <= N);
    if(p1->orientation){
      double shift = X[qryidx] - X[Lqryidx];
      refloc = (Y[Rrefidx] - shift) * 1e3;
    } else {
      double shift = X[Rqryidx] - X[qryidx];
      refloc = (Y[Rrefidx] - shift) * 1e3;
    }
    refidx   = Rrefidx;
    ret = Rindex;
    if(qryloc){
      *qryloc = X[qryidx] * 1e3;
      if(verbose){
	printf("qryidx -> %d, qryloc -> %0.3f\n", qryidx, *qryloc * 1e-3);
	fflush(stdout);
      }
    }
  }
  return ret;
} //end getRef


/* After the introduction of filtering criteria for translocations, two loops
   are required. One to make the SVs, one to filter them. Store the information
   from the first loop in this object for use in the second. */
class structuralVariation {
public:
  int indices[2]; //indices in array of xmapEntry objects -- note that which array this is is only defined when using these
  //whether or not to apply filtering of overlap events of the same type
  // (ie, if two translocations conflict, one will be called 'overlap') 
  int link_index; // xmapid of inverted MG (currently only for inversion_partial entries, see makeInversionPartial)
  int pair_index; //for paired svs: index into sv_array; currently inversions only, see pairInversionBP
  bool exclude; //exclude means overlap
  bool duplicate; //cominatorial duplicates (1-2, 1-3, 2-3: 1-3 is duplicate)
  bool repeat; //repeat label pattern, see simpleRepeat.{cpp,h}
  bool common; //matches entry in mask file
  double querysize, refsize;/* gap size : in kb */
  double querystart, querystop, refstart, refstop;/* SV coordinates (indel or breakpoint) in bp*/
  int querystartidx, querystopidx, refstartidx, refstopidx;/* SV coordinates (indel or breakpoint gap) : absolute coordinates for query, with querystartidx < querystopidx */
  double frequency, qrycov, refcov;/* number between 0.0 and 1.0 : fraction of query coverage vs total coverage of ref based on average of previous 4 coordinates */
  char *orientation;/* NULL or Translocation orientation string */
  long long refcontigid1, refcontigid2;
  int xmapid1, xmapid2; //store for sorting
  char* svtype; //type in smap
  const char* svtype_debug; //debugging string only--encodes explanation of why complex event failed to be any other type
  sv_type_t type_enum;
  int gap_overlap; //overlap with gap from bed
  double confidence, confidence_left, confidence_right, confidence_center;
  float confidence_scaled;
  sv_algo_t algo_enum;
  bool alignment_overlap;
  int indel_overlap;
  long long indel_overlap_id;
  int genotype, genotypeGroup;
  bool link_coord; // boolean condition on coordinates of link_index's MG: see makeInversionPartial (NOT USED ANYMORE)
  int zygosity; // -1=unassigned, 0=hom, 1=het, 2=unk
  int smap_index;/* index in .smap output : precomputed during no-write loop, so pair_index can be translated to location in .smap */
  bool inv_MGi;// true IFF MG i (indices[0]) is the inverted matchgroup (only valid for inversion OR inverted duplication, otherwise set to false)
  bool MGi_qryleft;// true IFF MG i (indices[0]) is left of MG j (indices[1]) on query (only computed for inversion, otherwise set to false)
  bool trans_overlapped;// true IFF the translocation was called based on small matchgroup overlapping insertion in larger matchgroup
  
  double refstart2,refstop2;/* for small inversion paired SV coordinates */

  ~structuralVariation() {
    if( svtype )
      delete [] svtype;
    if(orientation != NULL){
      free(orientation);
      orientation = NULL;
    }
    //if( svtype_debug ) //may get double-free error
    //delete [] svtype_debug;
  }
  structuralVariation() {
    indices[0] = indices[1] = link_index = pair_index = -1;
    querysize = refsize = querystart = querystop = refstart = refstop = 0;
    querystartidx = querystopidx = refstartidx = refstopidx = 0;
    svtype = 0;
    svtype_debug = 0;
    algo_enum = sv_unknown;
    refcontigid1 = refcontigid2 = 0;
    exclude = duplicate = repeat = common = alignment_overlap = false;
    gap_overlap = 0;
    type_enum = unknown_type;
    xmapid1 = xmapid2 = indel_overlap = indel_overlap_id = 0;
    confidence_left = confidence_right = confidence_center = confidence = -1000.0;
    confidence_scaled = -1.;
    genotype = genotypeGroup = -1;
    link_coord = false;
    smap_index = zygosity = -1;
    frequency = refcov = qrycov = -1.0;
    orientation = NULL;
    trans_overlapped = false;
    refstop2 = refstart2 = -1.0;
  }

  structuralVariation(int ind1, int ind2, double qsize, double rsize, double qstart, double qstop, long long refid1, long long refid2, double rstart, double rstop, const char* svtd, sv_type_t svt, int qsi=-1, int qei=-1, int rsi=-1, int rei=-1, double conf=-1., double lconf=-1., double rconf=-1., sv_algo_t al=sv_indel, int xid1=-1, int xid2=-1, bool trans_ovr = false) {
    indices[0]	  = ind1;
    indices[1]	  = ind2;
    link_index    = -1;
    pair_index    = -1;
    querysize	  = qsize;
    refsize	  = rsize;
    if(DEBUG) assert(qstart >= 0.0);
    querystart	  = qstart;
    querystop	  = qstop;
    querystartidx = qsi;
    querystopidx  = qei;
    refcontigid1  = refid1;
    refcontigid2  = refid2;
    if(DEBUG) assert(rstart >= 0.0);
    refstart	  = rstart;
    refstop	  = rstop;
    refstartidx   = rsi;
    refstopidx    = rei;
    exclude       = duplicate = repeat = common = alignment_overlap = false;
    svtype        = 0; //svtyp;
    svtype_debug  = svtd;
    type_enum     = svt;
    //printf("constructor: type=%d\n", type_enum); //debug
    gap_overlap   = 0;
    confidence_center = conf;
    confidence_left   = lconf;
    confidence_right  = rconf;
    confidence        = (conf > 0 && lconf > 0 && rconf > 0) ? makeConfidence(conf, lconf, rconf) : -1.;
    confidence_scaled = -1.;
    algo_enum     = al;
    xmapid1       = xid1;
    xmapid2       = xid2;
    indel_overlap = indel_overlap_id = 0;
    genotype = genotypeGroup = -1;
    link_coord = false;
    smap_index = zygosity = -1;
    inv_MGi = false;
    MGi_qryleft = false;
    frequency = refcov = qrycov = -1.0;
    orientation = NULL;
    trans_overlapped = trans_ovr;

    refstart2 = refstart;
    refstop2 = refstop;

#if 0 // NOTE : not valid for indels ???
    if(VERB>=2 && indices[0] != indices[1]){
      xmapEntry *entry1 = (algo_enum == sv_ps) ? usexmapentries[indices[0]] : xmapentries[indices[0]];//copy from end of output_smap
      xmapEntry *entry2 = (algo_enum == sv_ps) ? usexmapentries[indices[1]] : xmapentries[indices[1]];
      printf("structuralVariation constructor: algo=%d,type=%d,indices[0]=%d,indices[1]=%d,xmapid1=%d(%d),xmapid2=%d(%d)\n",
	     algo_enum, type_enum, indices[0],indices[1],xmapid1,entry1->xmapid,xmapid2,entry2->xmapid);
      fflush(stdout);
    }
#endif
  }

  //if -bed argument is used, numbedentries is the number of lines read from
  //bed file, stored as bedEntry objects
  //for indels and inversions,
  //check if this SV overlaps with any 'gap' type entries, set gap_overlap
  //add flagging translocations also if bed entry type is 'common' or 'segdupe'
  void checkBedGapOverlap() {
    if( duplicate || exclude ) //do not check duplicates, overlaps
      return;
    //if false, compare translocations to 'common' and 'segdupe', if true, compare indels and inversions to gaps
    bool checkgap = false; 
    if(type_enum == insertion || type_enum == deletion || type_enum == compression || type_enum == inversion_type) //types compared to gaps
      checkgap = true;
    //no other types considered
    else if(type_enum != interchr_translocation && type_enum != intrachr_translocation) 
      return;

    //if no bed file supplied, numbedentries == 0
    for(int bi=0; bi < numbedentries; bi++) {
      //be sure bed entry type is 'gap'--same reference contig (only check 1 bc they're equal for indels and inversions)
      if( checkgap && (strcmp(bedentries[bi]->type,"gap") != 0 || bedentries[bi]->contigid != refcontigid1) )
	continue;
      if( !checkgap && strcmp(bedentries[bi]->type,"common") != 0 && strcmp(bedentries[bi]->type,"segdupe") != 0 )
	continue;
      //because of intra_chr translocations, cannot assume single ref contig id
      /* int checkid = 0; //1 = check refcontigid1, 2 = refcontigid2
      if( !checkgap && bedentries[bi]->contigid == refcontigid1 )
	checkid = 1;
      else if( !checkgap && bedentries[bi]->contigid == refcontigid2 )
	checkid = 2;
      else if( !checkgap ) // && checkid == 0 )
	continue; //neither chromosome matches
      */

      if( !checkgap ) {
	//translocations: no buffer is necessary: just compare the right chromosome
	//for translocations, set common to true if bed entry is 'common' and set gap_overlap = 3 for 'segdupe'
	bool match = false;
	for( int i=0; i<2; i++ ) {
	  //0/1 doesn't matter, but refstart is refcontigid1, so must match ternary below
	  int refid = ( i == 1 ? refcontigid1 : refcontigid2 ); 
	  if( bedentries[bi]->contigid != refid )
	    continue;
	  double rpos = ( i == 1 ? refstart : refstop );
	  if( bedentries[bi]->start < rpos && bedentries[bi]->stop > rpos ) {
	    match = true;
	    if( strcmp(bedentries[bi]->type,"common") == 0 )
	      common = 1;
	    else if( strcmp(bedentries[bi]->type,"segdupe") == 0 )
	      gap_overlap = 3;
	    if(VERB>=2){
	      xmapEntry *entry1 = (algo_enum == sv_ps) ? usexmapentries[indices[0]] : xmapentries[indices[0]];//copy from end of output_smap
	      printf("type=%d:xmapid=%d,%d,refid=%lld,%lld,qryid=%lld,%s=%0.1f:bedentries[%d]:start=%0.1f,stop=%0.1f:common=%d,gap_overlap=%d\n",
		     type_enum, xmapid1,xmapid2,refcontigid1,refcontigid2,entry1->qrycontigid,(i==1) ? "refstart":"refstop", rpos,bi,bedentries[bi]->start,bedentries[bi]->stop,common,gap_overlap);
	      fflush(stdout);
	    }
	    break; //first match is sufficient
	  }
	}
	if( match )
	  break; //take only first match
	continue; //do not execute code below (that is for checkgap == true)
      }

      //overlap criteria: since labels inside gaps aren't possible,
      // the SV can neither start nor end inside the gap. Lets record such events as an error (NOTE : often the Bed file is the default human file, so this should NOT be an error)
      // encode the error as gap_overlap == 2
      //therefore, a correct 'gap' sv must have one (or more) gaps inside of it
      if( bedentries[bi]->start >= refstart && bedentries[bi]->stop <= refstop ) {
	gap_overlap = 1; //correct gap
	break; //no need to check for others
      }
      else if( (bedentries[bi]->start >= refstart && bedentries[bi]->start < refstop) || //gap starts inside sv
	       (bedentries[bi]->stop  > refstart && bedentries[bi]->stop  <= refstop) || //gap stops inside sv
	       (bedentries[bi]->start <  refstart && bedentries[bi]->stop  >  refstop) ) { //gap spans sv
	gap_overlap = 1 /* WAS81 2*/; //error--the gap must be in the wrong place (or there's a bug somewhere)
	if(DEBUG>=1+RELEASE) {
	  printf("checkBedGapOverlap:bedentries[bi=%d]->start= %0.1f, stop= %0.1f, refstart= %0.1f, refstop= %0.1f bp : invalid gap overlap condition (incorrect bed file ?)\n",
		 bi,bedentries[bi]->start,bedentries[bi]->stop, refstart,refstop);
	  fflush(stdout);
	  //	  assert(0);/* investigate if this is a valid case */
	}
	break;
      }
    } //end for numbedentries
  } //end checkBedGapOverlap

  //svtype is not set in either constructor above
  //set it here based on type_enum, exclude, and duplicate
  //also set svtype_debug if not set

  void setTypeStr() {
    size_t retLen = strlen(transintra_name) + strlen("_nbaseerror");// 22 + 11 == 33 byte maximum string length
    char* ret = new char[retLen + 1]; // + 1 for string termination
    ret[0] = '\0';// initialize to empty string for debugging to catch unhandled cases
    assert(type_enum != unknown_type); //these should not exist
    if(type_enum == insertion)
      strcpy(ret, "insertion"); //9
    else if(type_enum == deletion || type_enum == compression)
      strcpy(ret, "deletion");
    else if(type_enum == complex_type)
      strcpy(ret, "complex");
    else if(type_enum == inversion_type && repeat && !duplicate) //if duplicate, don't use repeat suffix
      strcpy(ret, "inversion_repeat"); //16
    else if(type_enum == inversion_type && (!repeat || duplicate)) //if duplicate, get suffix below
      strcpy(ret, "inversion");
    else if(type_enum == inversion_small_type && (!repeat || duplicate))
      strcpy(ret, "inversion_paired");
    //for translocations, the names are not "overlap" appended to above
    else if(type_enum == intrachr_translocation) {
      if(duplicate) //duplicate is worse than repeat
	strcpy(ret, transintra_duplicate_name);
      else if(repeat) //repeat is first, ie, 'worst' or 'most important'
	strcpy(ret, transintra_repeat_name);
      else if(exclude) //exclude means overlap
	strcpy(ret, transintra_overlap_name);
      else if(common) //common means masked
	strcpy(ret, transintra_common_name);
      else if(gap_overlap == 3) //1 and 2 are not used for translocations, 3 is segdupe
	strcpy(ret, transintra_segdup_name);
      else //none of above is normal event
	strcpy(ret, transintra_name);
    }
    else if(type_enum == interchr_translocation) {
      if(duplicate)
	strcpy(ret, transinter_duplicate_name);	
      else if(repeat)
	strcpy(ret, transinter_repeat_name);	
      else if(exclude)
	strcpy(ret, transinter_overlap_name);
      else if(common) //common means masked
	strcpy(ret, transinter_common_name);
      else if(gap_overlap == 3)
	strcpy(ret, transinter_segdup_name);
      else
	strcpy(ret, transinter_name);
    }
    else if(type_enum == duplication_type)
      strcpy(ret, "duplication");
    else if(type_enum == duplicationinv_type)
      strcpy(ret, "duplication_inverted");
    else if(type_enum == duplicationsplit_type)
      strcpy(ret, "duplication_split");

    //below are annotations on top of above types
    //duplicate takes precedence over overlap
    // rational is that if it's duplicate, you know it's garbage, whereas overlap still has some possibility of being correct
    if(duplicate && (type_enum == insertion || type_enum == deletion || type_enum == compression || type_enum == inversion_type))
      strcat(ret, "_duplicate"); //22+10 < 32
    //now do overlap for indels + inversions, but not if repeat
    else if(exclude && !(type_enum == inversion_type && repeat) &&
	    type_enum != complex_type && 
	    type_enum != intrachr_translocation && 
	    type_enum != interchr_translocation )
      strcat(ret, "_overlap"); //8+9 < 22
    //n-base gap overlap, only indel and inversion (not repeat): above two are 'worse', call those first
    else if( gap_overlap == 1 && !(type_enum == inversion_type && repeat) && 
	     type_enum != intrachr_translocation && type_enum != interchr_translocation )
      strcat(ret, "_nbase"); //6+9 < 22
    else if( gap_overlap == 2 && !(type_enum == inversion_type && repeat) &&
	     type_enum != intrachr_translocation && type_enum != interchr_translocation ) //see checkBedGapOverlap
      strcat(ret, "_nbaseerror"); //11+9 < 22
    //introduce 'tiny' type for indels: size is less than smap_indel_tiny_maxsize, default 1kb
    else if( (type_enum == insertion || type_enum == deletion || type_enum == compression) &&
	     fabs(querysize - refsize) < smap_indel_tiny_maxsize )
      strcat(ret, "_tiny"); //5+9 < 22

    //printf("ret %s (%zd, %zd), type %d\n", ret, strlen(ret), retLen, int(type_enum)); //debug
    if(DEBUG && (type_enum == unknown_type || ret[0]=='\0' || strlen(ret) > retLen)){
      printf("ERROR in setTypeStr: type_enum=%d,repeat=%d,duplicate=%d, ret=%s\n", type_enum, repeat,duplicate?1:0,ret);
      fflush(stdout);
      assert(ret[0] && strlen(ret) <= retLen);
    }
    //svtype = ret; //rather than assign, copy and delete orig
    svtype = new char[retLen+1];
    strcpy(svtype, ret);
    delete [] ret;
    //also set debug str if not already set
    if( svtype_debug == 0 ) {
      svtype_debug = svtype;
    }
  } //end setTypeStr

  //for inversion SVs, need to add inversion_partial to smap
  // these are the best guess at the other breakpoint of the inversion
  // do this only if svtype == "inversion", and write the new smap
  // entry to FILE* fp. Return 0 for no SV added, 1 for added.
  void makeInversionPartial(xmapEntry **usexmapentries, bool verbose=true) {
    if( !doInversionPartial() )
      return;

    xmapEntry *entry1 = (algo_enum == sv_ps) ? usexmapentries[indices[0]] : xmapentries[indices[0]];//copy from end of output_smap
    xmapEntry *entry2 = (algo_enum == sv_ps) ? usexmapentries[indices[1]] : xmapentries[indices[1]];

    // qry orientation is based on which MG is the inverted matchgroup (previously computed, see inv_MGi) : The MG without the inversion determines the orientation of qry. 
    // The MG with the inversion has the inversion_partial coordinates.
    bool qry_orientf = inv_MGi ? entry2->orientforward : entry1->orientforward;

    //xmapEntry* entry_inv = (entry1->orientforward != qry_orientf ? entry1 : entry2);
    if(DEBUG/* HERE >= 2*/) assert(inv_MGi == (entry1->orientforward != qry_orientf));
    link_index = inv_MGi ? indices[0] : indices[1];

    xmapEntry* entry_inv = (algo_enum == sv_ps) ? usexmapentries[link_index] : xmapentries[link_index];
    link_coord = (entry_inv->qrystartpos != querystart && entry_inv->qrystartpos != querystop);
    if( VERB>=2 && verbose ){
      printf("Inv part: xmapids=%d,%d, inv_MG= %d, qry_orientf %i, link_index %i(xmapid=%d), link_coord %i\n", 
	     entry1->xmapid,entry2->xmapid, inv_MGi ? entry1->xmapid : entry2->xmapid, qry_orientf, link_index, entry_inv->xmapid,link_coord);
      fflush(stdout);
    }
  } //end makeInversionPartial

  int writeInversionPartial(xmapEntry **usexmapentries, int nsv, FILE *fp, const char* zyg, double SVsize) {
    assert(svtype != 0); //setTypeStr must have assigned this char* a non-zero value

    if(DEBUG) assert(link_index >= 0);// link_index should already have been computed
    xmapEntry* entry_inv = (algo_enum == sv_ps) ? usexmapentries[link_index] : xmapentries[link_index];
    bool qry_orientf = !entry_inv->orientforward;

    double thisqry = 0, thisref = 0; // coords for inversion_partial smap entry
    int thisqryidx = 0, thisrefidx = 0; // corresponding label indices

    if(link_index == indices[0]){/* inverted matchgroup is MG i (left matchgroup on reference) */
      /* partial locations are obtained by mapping refstart and (qry_orientf ? querystart : querystop) via inverted matchgroup alignment */
      getQry(entry_inv->align, refstartidx, thisqry, thisqryidx);
      
      if(MGi_qryleft)/* MG i is also the left matchgroup on query */
	getRef(entry_inv->align, qry_orientf ? querystartidx : querystopidx , thisref, thisrefidx);
      else /* MG i is the right matchgroup on query : this is rare and typically means that MG i completely overlaps MG j on reference */
	getRef(entry_inv->align, qry_orientf ? querystopidx : querystartidx , thisref, thisrefidx);

    } else { /* inverted matchgroup is MG j (right matchgroup on reference) */
      if(DEBUG) assert(link_index == indices[1]);
      /* partial locations are obtained by mapping refstop and (qry_orientf ? querystop : querystart) via inverted matchgroup alignment */
      getQry(entry_inv->align, refstopidx, thisqry, thisqryidx);

      /* MG i is also the left matchgroup on query  */
      if(MGi_qryleft)
	getRef(entry_inv->align, qry_orientf ? querystopidx : querystartidx , thisref, thisrefidx);
      else /* NEW101 : this is rare : the overlapped part of the inverted matchgroup is assumed to be trimmed away */
	getRef(entry_inv->align, qry_orientf ? querystopidx : querystartidx , thisref, thisrefidx, qry_orientf ? 1 : -1, qry_orientf ? &querystop : &querystart );	    
    }

    if(VERB>=2){
      xmapEntry *entry1 = (algo_enum == sv_ps) ? usexmapentries[indices[0]] : xmapentries[indices[0]];
      xmapEntry *entry2 = (algo_enum == sv_ps) ? usexmapentries[indices[1]] : xmapentries[indices[1]];
      printf("Inv part: xmapid=%d,%d, refid=%lld,qryid=%lld: qry_orientf=%d, MGi_qryleft=%d, refgap=%d,%d (%0.3f,%0.3f), querygap=%d,%d (%0.3f,%0.3f): part: thisrefidx=%d,thisref=%0.3f,thisqryidx=%d,thisqry=%0.3f,SVsize= %0.3f\n",
	     entry1->xmapid,entry2->xmapid,entry_inv->refcontigid,entry_inv->qrycontigid, qry_orientf ? 1 : 0, MGi_qryleft ? 1 : 0,
	     refstartidx,refstopidx,refstart*1e-3,refstop*1e-3,querystartidx,querystopidx,querystart*1e-3,querystop*1e-3, thisrefidx,thisref*1e-3,thisqryidx,thisqry*1e-3,SVsize);
      fflush(stdout);
    }

    //write each entry to file--there is only one MG, so put -1 in all the second MG place
    writeToSmap(fp, nsv, entry_inv->qrycontigid, entry_inv->refcontigid, -1, thisqry, -1, thisref, -1, "inversion_partial", entry_inv->xmapidfilt, -1, nsv-1, thisqryidx, -1, thisrefidx, -1, -1, -1, -1., -1., zyg, -1, -1, confidence_scaled,SVsize,frequency,qrycov,refcov,orientation);
    return 1;
  } //end writeInversionPartial

  //utility method for makeInversionPartial and writeSVToSmap: return whether or not partial entry should be made
  bool doInversionPartial() {
    return type_enum == inversion_type && !duplicate && !exclude;
    //above condition is equivalent to below after setTypeStr is called (use above so that this doesn't need to be called)
    // return !( strcmp(svtype, "inversion") != 0 &&
    // 	      strcmp(svtype, "inversion_repeat") != 0 &&
    // 	      strcmp(svtype, "inversion_nbase") != 0 &&
    // 	      strcmp(svtype, "inversion_nbaseerror") != 0 );
  }

  //wrapper to makeInversionPartial which also calls fn writeToSmap
  // If dowrite== false just mark corresponding xmapentry->output = true

  int writeSVToSmap(xmapEntry **usexmapentries, int nsv, FILE *fp, bool dowrite, bool verbose=true) {
    //duplicate SV calls are artifact of pairsplit (or, less often, refsplit) and therefore uninteresting
    if ( duplicate ) {
      if (VERB>=3 && (type_enum == inversion_small_type || this - sv_array == 5640)){
	printf("writeSVToSmap: sv_array[%ld] : nsv=%d, pair_index=%d, smap_index=%d, type_enum=%d : duplicate= true\n",this - sv_array, nsv, pair_index, smap_index, type_enum);
	fflush(stdout);
      }
      return 0;
    }
    //confidence filter: confidence -1 is ok for non-indel;
    // for indel, confidence == 0 is not an indel, so these should always be omitted
    if( (confidence > 0.0 && confidence < smap_min_conf) || confidence == 0.0 ){
      if (VERB>=3 && (type_enum == inversion_small_type || this - sv_array == 5640)){
	printf("writeSVToSmap: sv_array[%ld] : nsv=%d, pair_index=%d, smap_index=%d, type_enum=%d : confidence = %0.4f (min= %0.4f)\n",
	       this - sv_array, nsv, pair_index, smap_index, type_enum, confidence, smap_min_conf);
	fflush(stdout);
      }
      return 0;
    }
    //this is criteria for which an inversion_partial event will be created
    //inversion partial is only created for these types, not inversion_duplicate, inversion_overlap (I think those are the only two left)
    bool invpart = doInversionPartial();

    //for -indel calls, the indices refer to the original xmap entry array, not usexmapentries (Same for inversion_small_type of Small Duplications)
    xmapEntry *entry1 = (algo_enum == sv_ps) ? usexmapentries[indices[0]] : xmapentries[indices[0]];
    xmapEntry *entry2 = (algo_enum == sv_ps) ? usexmapentries[indices[1]] : xmapentries[indices[1]];
    if(DEBUG/* HERE >=2 */ && type_enum == inversion_small_type )  assert( algo_enum == sv_indel );
    entry1->output = true;
    entry2->output = true;
    if(VERB>=2){
      printf("writeSVToSmap:sv_array[%ld]:dowrite=%d,nsv=%d:algo=%d,type=%d,indices[0]=%d,indices[1]=%d,xmapid1=%d(%d),xmapid2=%d(%d),output1=%d,output2=%d,invpart=%d,link_index=%d\n",
	     this - sv_array, dowrite, nsv,algo_enum, type_enum, indices[0],indices[1],xmapid1,entry1->xmapid,xmapid2,entry2->xmapid,entry1->output?1:0,entry2->output?1:0,invpart,link_index);
      fflush(stdout);
    }
    if(DEBUG && !(xmapid1 == entry1->xmapid && xmapid2 == entry2->xmapid)){
      printf("writeSVToSmap:sv_array[%ld]:dowrite=%d,nsv=%d:algo=%d,type=%d,indices[0]=%d,indices[1]=%d,xmapid1=%d(%d),xmapid2=%d(%d),output1=%d,output2=%d,invpart=%d,link_index=%d\n",
	     this - sv_array, dowrite, nsv,algo_enum, type_enum, indices[0],indices[1],xmapid1,entry1->xmapid,xmapid2,entry2->xmapid,entry1->output?1:0,entry2->output?1:0,invpart,link_index);
      printf("\t entry1->xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid, entry1->qrycontigid,
	     entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	     entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	     entry1->orientforward, entry1->confidence, entry1->align->numpairs);
      printf("\t entry2->xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry2->xmapid, entry2->refcontigid,entry2->qrycontigid,
	     entry2->qrystartpos * 1e-3, entry2->qryendpos * 1e-3, entry2->qrystartidx, entry2->qryendidx, 
	     entry2->refstartpos * 1e-3, entry2->refendpos * 1e-3, entry2->refstartidx, entry2->refendidx,
	     entry2->orientforward, entry2->confidence, entry2->align->numpairs);
      fflush(stdout);

      assert(xmapid1 == entry1->xmapid);
      assert(xmapid2 == entry2->xmapid);
    }

    //dummy check for ref contig id
    assert(entry1->refcontigid == refcontigid1 && entry2->refcontigid == refcontigid2);

    if( !dowrite ){
      smap_index = nsv /* WAS (type_enum == inversion_small_type) ? nsv : -1 */;// only used by type_enum == inversion_small_type
      
      if(VERB>=2){
	printf("writeSVToSmap: sv_array[%ld] : nsv=%d, pair_index=%d, smap_index=%d, type_enum=%d, dowrite=false\n", this - sv_array, nsv, pair_index, smap_index, type_enum);
	fflush(stdout);
      }

#if 0
      if(invpart && algo_enum != sv_ps){/* also mark the xmapids that makeInversionPartial would use */
	xmapEntry *entry1 = (algo_enum == sv_ps) ? usexmapentries[indices[0]] : xmapentries[indices[0]];
	xmapEntry *entry2 = (algo_enum == sv_ps) ? usexmapentries[indices[1]] : xmapentries[indices[1]];
	printf("invpart && !ps: ids: %i %i\n", entry1->xmapid, entry2->xmapid);
	fflush(stdout);
	assert(false);
	entry1->output = true;
	entry2->output = true;
      }
#endif

      return (invpart ? 2 : 1);// number of smap entries that would have been added
    }

    const char* zyg = getZygosity(); //this can change value of genotype, so call before next line

    int SVlink_index = invpart ? nsv+1 : -1;// index into .smap of partial inversion (NOTE : not to be confused with SV link_index field, which is the index into usexmapentries[] for inverted matchgroup)

    if(DEBUG /* WAS && type_enum == inversion_small_type */) assert(smap_index == nsv);/* smap_index was set during the no-write sv_array[] traversal (see above) */

    if(type_enum == inversion_small_type){
      if(VERB>=3){
	printf("writeSVToSmap: sv_array[%ld] : nsv=%d, pair_index=%d, smap_index=%d, type_enum=%d, type=%s, confidence= %0.4f, dowrite=true\n", 
	       this - sv_array, nsv, pair_index, smap_index, type_enum, svtype, confidence);
	fflush(stdout);
      }

      if(DEBUG && !(pair_index >= 0)){
	printf("writeSVToSmap: sv_array[%ld] : nsv=%d, pair_index=%d, smap_index=%d, type_enum=%d, dowrite=false\n", this - sv_array, nsv, pair_index, smap_index, type_enum);
	fflush(stdout);
	assert(pair_index >= 0);
      }

      double pconfidence = sv_array[pair_index].confidence;
      if(DEBUG && !(!sv_array[pair_index].duplicate && !((pconfidence > 0.0 && pconfidence < smap_min_conf) || pconfidence == 0.0))){
	printf("WARNING: SVlink_index=%d : sv_array[pair_index=%d]: duplicate=%d,confidence=%0.4f (min_conf= %0.4f)\n",
	       SVlink_index,pair_index,sv_array[pair_index].duplicate,sv_array[pair_index].confidence,smap_min_conf);
	fflush(stdout);
      }

      SVlink_index = sv_array[pair_index].smap_index;
      
      if(DEBUG && !(SVlink_index >= 1)){
	printf("writeSVToSmap: sv_array[%ld] : nsv=%d, pair_index=%d, smap_index=%d, type_enum=%d, type=%s, dowrite=true\n", this - sv_array, nsv, pair_index, smap_index, type_enum, svtype);
	printf("\t xmapid1=%d,xmapid2=%d,SVlink_index = %d\n", xmapid1,xmapid2,SVlink_index);
	fflush(stdout);
	assert(SVlink_index >= 1);
      }
    }

    double SVsize = -1.0;
    if(SV_SIZE){
      if(type_enum == insertion) {
	SVsize = fabs(querystop - querystart) - (refstop - refstart);
	if(DEBUG>=1+RELEASE && !(SVsize >= 0.0)){
	  printf("writeSVToSmap:sv_array[%ld]:dowrite=%d,nsv=%d:algo_enum=%d,indices[0]=%d,indices[1]=%d,xmapid1=%d(%d),xmapid2=%d(%d),output1=%d,output2=%d,invpart=%d,link_index=%d\n",
		 this - sv_array, dowrite, nsv,algo_enum, indices[0],indices[1],xmapid1,entry1->xmapid,xmapid2,entry2->xmapid,entry1->output,entry2->output,invpart,link_index);
	  printf("\t type_enum=%d:querystop= %0.1f, querystart= %0.1f, refstop= %0.1f, refstart= %0.1f bp\n",type_enum,querystop,querystart,refstop,refstart);
	  fflush(stdout);
	  assert(SVsize >= 0.0);
	}
	SVsize = fabs(SVsize);
      } else if(type_enum == deletion || type_enum == compression) {
	SVsize = (refstop - refstart) - fabs(querystop - querystart);
	if(DEBUG>=1+RELEASE) assert(SVsize >= 0.0);
	SVsize = fabs(SVsize);
      } else if(type_enum == duplication_type || type_enum == duplicationinv_type || type_enum == duplicationsplit_type){
	SVsize = refstop - refstart;
	if(DEBUG>=1+RELEASE && !(SVsize >= 0.0)){
	  printf("writeSVToSmap:sv_array[%ld]:refid=%lld,qryid=%lld:dowrite=%d,nsv=%d:algo_enum=%d,indices[0]=%d,indices[1]=%d,xmapid1=%d(%d),xmapid2=%d(%d),output1=%d,output2=%d,invpart=%d,link_index=%d\n",
		 this - sv_array, entry1->refcontigid,entry1->qrycontigid,dowrite, nsv,algo_enum, indices[0],indices[1],
		 xmapid1,entry1->xmapid,xmapid2,entry2->xmapid,entry1->output,entry2->output,invpart,link_index);
	  printf("\t type_enum=%d:querystop= %0.1f, querystart= %0.1f, refstop= %0.1f, refstart= %0.1f bp\n",type_enum,querystop,querystart,refstop,refstart);
	  fflush(stdout);
	  assert(SVsize >= 0.0);
	}
	SVsize = fabs(SVsize);
      } else if(type_enum == inversion_type){
	/* NEW38 : duplicate code from writeInversionPartial to compute inversion_partial coordinates */
	if(DEBUG) assert(link_index >= 0);// link_index should already have been computed
	xmapEntry* entry_inv = (algo_enum == sv_ps) ? usexmapentries[link_index] : xmapentries[link_index];/* inverted matchgroup */
	bool qry_orientf = !entry_inv->orientforward;

	double thisqry = 0, thisref = 0; // coords for inversion_partial smap entry
	int thisqryidx = 0, thisrefidx = 0; // corresponding label indices

	if(link_index == indices[0]){/* inverted matchgroup is MG i (left matchgroup on reference) */
	  /* partial locations are obtained by mapping refstart and (qry_orientf ? querystart : querystop) via inverted matchgroup alignment */
	  getQry(entry_inv->align, refstartidx, thisqry, thisqryidx);
      
	  if(MGi_qryleft) { /* MG i is also the left matchgroup on query */
	    getRef(entry_inv->align, qry_orientf ? querystartidx : querystopidx , thisref, thisrefidx);
	    SVsize = refstart - thisref;
	  } else { /* MG i is the right matchgroup on query : this is rare and typically means that MG i completely overlaps MG j on reference */
	    getRef(entry_inv->align, qry_orientf ? querystopidx : querystartidx , thisref, thisrefidx);
	    SVsize = thisref - refstart;
	  }

	  if(DEBUG) assert(refstart <= refstop);
	  if(DEBUG>=1+RELEASE && !(SVsize >= 0.0)){
	    printf("writeSVToSmap:sv_array[%ld]:dowrite=%d,nsv=%d:algo=%d,indices[0]=%d,indices[1]=%d,xmapid1=%d(%d),xmapid2=%d(%d),output1=%d,output2=%d,invpart=%d,link_index=%d\n",
		   this - sv_array, dowrite, nsv,algo_enum, indices[0],indices[1],xmapid1,entry1->xmapid,xmapid2,entry2->xmapid,entry1->output?1:0,entry2->output?1:0,invpart,link_index);
	    printf("\t xmapid1=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
		   entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
		   entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
		   entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	    printf("\t xmapid2=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry2->xmapid, entry2->refcontigid,entry2->qrycontigid,
		   entry2->qrystartpos * 1e-3, entry2->qryendpos * 1e-3, entry2->qrystartidx, entry2->qryendidx, 
		   entry2->refstartpos * 1e-3, entry2->refendpos * 1e-3, entry2->refstartidx, entry2->refendidx,
		   entry2->orientforward, entry2->confidence, entry2->align->numpairs);

	    printf("\t type_enum=%d:entry_inv->xmapid=%d,fwd=%d:inv_MGi=%d,MGi_qryleft=%d:refstart=%0.1f,refstop=%0.1f,thisref=%0.1f,SVsize=%0.1f bp\n",
		   type_enum,entry_inv->xmapid,entry_inv->orientforward,inv_MGi,MGi_qryleft,refstart,refstop, thisref,SVsize);
	    fflush(stdout);
	    assert(SVsize >= 0.0);
	  }
	  SVsize = fabs(SVsize);

	} else { /* inverted matchgroup is MG j (right matchgroup on reference) */
	  if(DEBUG) assert(link_index == indices[1]);
	  /* partial locations are obtained by mapping refstop and (qry_orientf ? querystop : querystart) via inverted matchgroup alignment */
	  getQry(entry_inv->align, refstopidx, thisqry, thisqryidx);

	  if(MGi_qryleft) 	  /* MG i is also the left matchgroup on query  */
	    getRef(entry_inv->align, qry_orientf ? querystopidx : querystartidx , thisref, thisrefidx, 0);
	  else {/* NEW101 : This is rare : the  overlapped part of the inverted matchgroup is assumed to be trimmed away */
	    getRef(entry_inv->align, qry_orientf ? querystopidx : querystartidx , thisref, thisrefidx, qry_orientf ? 1 : -1, qry_orientf ? &querystop : &querystart );	    
	    
	    if(VERB/* HERE HERE >=2 */){
	      xmapEntry *entry1 = (algo_enum == sv_ps) ? usexmapentries[indices[0]] : xmapentries[indices[0]];
	      xmapEntry *entry2 = (algo_enum == sv_ps) ? usexmapentries[indices[1]] : xmapentries[indices[1]];
	      printf("Inversion: xmapid=%d,%d, refid=%lld,qryid=%lld: qry_orientf=%d, MGi_qryleft=%d, refgap=%d,%d (%0.3f,%0.3f), querygap=%d,%d (%0.3f,%0.3f): part: thisrefidx=%d,thisref=%0.3f,thisqryidx=%d,thisqry=%0.3f,SVsize= %0.3f\n",
		     entry1->xmapid,entry2->xmapid,entry_inv->refcontigid,entry_inv->qrycontigid, qry_orientf ? 1 : 0, MGi_qryleft ? 1 : 0,
		     refstartidx,refstopidx,refstart*1e-3,refstop*1e-3,querystartidx,querystopidx,querystart*1e-3,querystop*1e-3, thisrefidx,thisref*1e-3,thisqryidx,thisqry*1e-3,thisref - refstop);
	      fflush(stdout);
	    }
	  }

	  if(DEBUG) assert(refstart <= refstop);
	  SVsize = thisref - refstop;
	  if(DEBUG>=1+RELEASE) assert(SVsize >= 0.0);
	  SVsize = fabs(SVsize);
	}

	if(VERB>=2){
	  xmapEntry *entry1 = (algo_enum == sv_ps) ? usexmapentries[indices[0]] : xmapentries[indices[0]];
	  xmapEntry *entry2 = (algo_enum == sv_ps) ? usexmapentries[indices[1]] : xmapentries[indices[1]];

	  printf("Inv SVsize: xmapid=%d,%d, refid=%lld,qryid=%lld: qry_orientf=%d, MGi_qryleft=%d, refgap=%d,%d (%0.3f,%0.3f), querygap=%d,%d (%0.3f,%0.3f): part: thisrefidx=%d,thisref=%0.3f,thisqryidx=%d,thisqry=%0.3f,SVsize= %0.3f\n",
	     entry1->xmapid,entry2->xmapid,entry_inv->refcontigid,entry_inv->qrycontigid, qry_orientf ? 1 : 0, MGi_qryleft ? 1 : 0,
	     refstartidx,refstopidx,refstart*1e-3,refstop*1e-3,querystartidx,querystopidx,querystart*1e-3,querystop*1e-3, thisrefidx,thisref*1e-3,thisqryidx,thisqry*1e-3,SVsize);
	  fflush(stdout);
	}
      } else if(type_enum == inversion_small_type){/* NOTE : the inverted (small) matchgroup is always xmapid2 */
	structuralVariation *that = &sv_array[pair_index];
	if(that->refstop > refstart){
	  if(DEBUG && !(refstop <= that->refstart)){
	    printf("this=&sv_array[%ld], pair_index= %d: refstart= %0.3f, refstop= %0.3f, that->refstart= %0.3f, that->refstop= %0.3f: refid= %lld,%lld\n",
		   this - sv_array, pair_index, refstart, refstop, that->refstart, that->refstop,refcontigid1,refcontigid2);
	    fflush(stdout);
	    assert(refstop <= that->refstart);
	  }
	  SVsize = ((that->refstop - refstart) + (that->refstart - refstop)) * 0.5;// NOTE : includes 1/2 the gap regions
	  if(DEBUG) assert(SVsize > 0.0);
	} else {
	  if(DEBUG) assert(that->refstop <= refstart);
	  SVsize = ((refstop - that->refstart) + (refstart - that->refstop)) * 0.5;// NOTE : includes 1/2 the gap regions
	  if(DEBUG) assert(SVsize > 0.0);
	}
      }
    }

    writeToSmap(fp, nsv, entry1->qrycontigid, entry1->refcontigid, entry2->refcontigid, querystart, querystop, refstart, refstop, svtype, entry1->xmapidfilt, entry2->xmapidfilt, 
		SVlink_index, querystartidx, querystopidx, refstartidx, refstopidx, confidence, confidence_left, confidence_right, confidence_center, zyg, genotype, genotypeGroup, confidence_scaled, SVsize, frequency, qrycov,refcov, orientation);
    if( invpart ) {
      writeInversionPartial(usexmapentries, nsv+1, fp, zyg, SVsize); //will return 0 if not inversion
    }
    return (invpart ? 2 : 1); //number of smap entries added
  } //end writeSVToSmap

  void maskTranslocation() {
    //if this is a translocation, scan global maskentries, if both chrs match and position is
    //within +- smap_translocation_mask_tol (parameters.h/cpp, in kb), set flag common to true
    //because of chr match requirement, intra vs inter is automatic
    if(type_enum != intrachr_translocation && type_enum != interchr_translocation)
      return;
    //printf("checking mask for trans: %2i %12.2f, %2i %12.2f\n", refcontigid1, refstart, refcontigid2, refstop); //debug
    for(int mi=0; mi<nummaskentries; mi++) {
      //need 0. otherwise declaration in globals is ambiguous; smap_translocation_mask_tol is in kb, convert to bp
      double chr1min = max(0., maskentries[mi].chr1_start - smap_translocation_mask_tol*1e3); 
      double chr1max = maskentries[mi].chr1_start + smap_translocation_mask_tol*1e3;
      double chr2min = max(0., maskentries[mi].chr2_start - smap_translocation_mask_tol*1e3); 
      double chr2max = maskentries[mi].chr2_start + smap_translocation_mask_tol*1e3;
      //it should be the case that refstart is always on refcontigid1 and refstop is always on refcontigid2
      //do not require that chrs are sorted here or in mask
      //note that there is an emacs bug which requries the < condition before the > condition or indentation is broken
      if( ( maskentries[mi].chr1 == refcontigid1 && refstart < chr1max && refstart > chr1min &&
	    maskentries[mi].chr2 == refcontigid2 && refstop  < chr2max && refstop  > chr2min ) ||
	  ( maskentries[mi].chr1 == refcontigid2 && refstart < chr2max && refstart > chr2min &&
	    maskentries[mi].chr2 == refcontigid1 && refstop  < chr1max && refstop  > chr1min ) ) {
	if(VERB){
	  printf("masked translocation: %2i  %2i  %2i  %12.2f  %12.2f\n", mi+1, maskentries[mi].chr1, maskentries[mi].chr2, maskentries[mi].chr1_start, maskentries[mi].chr2_start);
	  fflush(stdout);
	}
	common = true;
      }
    }
  } //end maskTranslocation

  void wrapperIndelConfidence(xmapEntry **usexmapentries, int verbose) {
    if( type_enum == compression ) {
      confidence = smap_min_conf;
      confidence_scaled = 0.01;
      return;
    }

    if( type_enum != insertion && type_enum != deletion ) //indel only
      return;
    if( duplicate ) //skip duplicates
      return;
    //confidence filter (for confidence set in output_xmap.cpp)
    if( confidence >= 0 && confidence < smap_min_conf )
      return;
    setIndelConfidence(usexmapentries, verbose);
    //now, confidence must be set and > smap_min_conf (sometimes it may be 0)
    if( confidence >= smap_min_conf )
      setIndelConfidenceScaled();
  }

  void setIndelConfidence(xmapEntry **usexmapentries, int verbose);

  // int sv_ppvbin_conf[sv_conflen];		 
  // int sv_ppvbin_sizekb[sv_sizelen];		 
  // float sv_ppvbin_del[sv_conflen][sv_sizelen];
  // float sv_ppvbin_ins[sv_conflen][sv_sizelen];

  int getConfBin(double conf) {
    for( int i=sv_conflen-1; i >= 0; i-- ) {
      if( conf >= (double)sv_ppvbin_conf[i] ) //stored as int
	return i;
    }
    printf("ERROR in getConfBin: size=%.3f, conf=%.2f\n", refsize - querysize, conf);
    fflush(stdout);
    assert(0); //dummy check
  }

  int getSizeBin(double size) {
    for( int i=sv_sizelen-1; i >= 0; i-- ) {
      if( size >= (double)sv_ppvbin_sizekb[i] ) //stored as int
	return i;
    }
    assert(0); //dummy check
  }

  void setIndelConfidenceScaled() {
    int cb = getConfBin(confidence);
    //NOTE: refsize and querysize are in bp for *internal outliers* (sv_indel) and in kb for *split indels* (sv_ps): needs to be kb for getSizeBin
    double size = fabs(refsize - querysize) * (algo_enum == sv_ps ? 1.0 : 1e-3);
    int sb = getSizeBin(size);
    if(size < svConfMinsize) /* below smallest value simulated : set to value that mean NA */
      confidence_scaled = svConfNAvalue;
    else if( querysize > refsize ) //insertion
      confidence_scaled = sv_ppvbin_ins[cb][sb];
    else
      confidence_scaled = sv_ppvbin_del[cb][sb];
    
    float raw_confidence_scaled = 1.0 - pow(0.1,confidence);// NEW66
    confidence_scaled = min(confidence_scaled, raw_confidence_scaled);

    //    printf("%9s: size=%7.3f (%10.1f %10.1f; %10.1f %10.1f : %10.1f %10.1f), conf=%6.2f, scaled=%.2f\n", type_enum == insertion ? "insertion" : "deletion", size, refstart, refstop, querystart, querystop, querysize, refsize, confidence, confidence_scaled); //debug--size is in kb, but start/stop are in bp
    fflush(stdout);
  }

  void setZygosity() {
    //use alignment_overlap and indel_overlap data members to determine zygosity interpretation
    if( !smap_zygosity || //put zygosity flag in smap
	(type_enum != insertion && type_enum != deletion && type_enum != compression &&
	 type_enum != intrachr_translocation && type_enum != interchr_translocation &&
	 type_enum != inversion_type && type_enum != inversion_small_type) )
      return;
    //zygosity: -1=unassigned, 0=hom, 1=het, 2=unk
    //see svOverlap (structuralVariation.cpp)
    if( !alignment_overlap ) {
      if( indel_overlap == 0 || indel_overlap == 1 ) {
	zygosity = 0; //ret = "homozygous";
      } else if( indel_overlap == 2 || indel_overlap == 3 ) {
	zygosity = 1; //ret = "heterozygous";
      }
    }
    else { //alignment_overlap
      if( indel_overlap > 1 && indel_overlap < 4 ) {
	genotype = -1; //for unknown zygosity, do not call genotype
	zygosity = 2; //ret = "unknown";
      }
      else if( indel_overlap == 0 || indel_overlap == 1 ) {
	zygosity = 1; //ret = "heterozygous";
      }
    }
    if(VERB>=2 && refcontigid1==refcontigid2 && refcontigid1 == 21){
      printf("setZygosity: svtype=%10s alignment_overlap= %i, indel_overlap= %i, zygosity= %d\n", svtype, alignment_overlap, indel_overlap, zygosity);      
      fflush(stdout);
    }

    if( (zygosity < 0 || zygosity > 2) && VERB && DEBUG ){
      printf("Error in getZygosity: alignment_overlap=%i, indel_overlap=%i\n", alignment_overlap, indel_overlap);
      fflush(stdout); assert(false);
    }
  } //setZygosity

  const char* getZygosity() {
    if( zygosity < 0 || zygosity > 2 ){
      return 0; //this is not an error: setZygosity does not change zygosity if it's not a type for which we compute zygosity
      //printf("Error in getZygosity: alignment_overlap=%i, indel_overlap=%i\n", alignment_overlap, indel_overlap);
      //fflush(stdout); assert(false);
    }

    //promote zygosity to zygosity of pair if 'greater': -1=unassigned, 0=hom, 1=het, 2=unk
    if( (type_enum == inversion_type || type_enum == inversion_small_type) && zygosity != 2 && pair_index > -1 ) {
      structuralVariation* sv = &sv_array[pair_index];
      if( DEBUG && (sv->zygosity < 0 || sv->zygosity > 2) ){
	printf("Error in getZygosity: invalid zygosity= %i in sv_array[pair_index= %i]: type=%d,contigid=%lld,%lld, ref= %0.3f .. %0.3f kb\n", 
	       sv->zygosity, pair_index, sv->type_enum, sv->refcontigid1,sv->refcontigid2,sv->refstart*1e-3,sv->refstop*1e-3);
	
	printf("\t Needed to compute zygosity of this sv=%p: type=%d, contigid=%lld,%lld, ref= %0.3f .. %0.3f kb\n",
	       this, type_enum, refcontigid1, refcontigid2, refstart*1e-3, refstop*1e-3);
	fflush(stdout); 
	assert(!(sv->zygosity < 0 || sv->zygosity > 2));
      }
      if( sv->zygosity > zygosity )
	zygosity = sv->zygosity;
    }

    if( zygosity == 0 )
      return "homozygous";
    if( zygosity == 1 )
      return "heterozygous";
    //if( zygosity == 2 )
    return "unknown";
  }
}; //end class structuralVariation

static Ident structuralVariation_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/structuralVariation.h 11194 2020-06-18 00:25:28Z tanantharaman $");

#endif
