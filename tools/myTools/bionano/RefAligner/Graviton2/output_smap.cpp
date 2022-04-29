#include <stdlib.h>
#include <stdio.h>
#ifndef WIN32
#include <unistd.h>
#else
#include <direct.h>
#define getcwd _getcwd
#endif

#include "globals.h"
#include "parameters.h"
#include "bedEntry.h"

#include "structuralVariation.h"
#include "Calign.h"

//#undef DEBUG
//#define DEBUG 2

#define MATCHGROUP_TRIM smap_MATCHGROUP_TRIM /* If matchgroups are fully overlapped on query, trim overlapped region of smaller matchgroup down to outlier regions : this does NOT affect _full.xmap
						>= 1 : Only if fully overlapped on both query and reference AND both matchgroups have same reference contig
						>= 2 : Even if NOT overlapped on reference
						>= 3 : Even if matchgroups have different reference contigs
					     */

#define MATCHGROUP_PARTIAL_TRIM smap_MATCHGROUP_PARTIAL_TRIM
                                  /* If matchgroups are partly overlapped on query, trim smaller matchgroup (for 2-MG duplications or overlapped Indels) OR 
				     inverted matchgroup (for 2-MG inversions or inverted duplications) until overlap in query is eliminated : does NOT affect _full.xmap.

				     Also trim matchgroups that have already been eliminated but will be output (due to 1-MG Indels and Small Inversions/Duplications) to exclude query regions overlapped by
				     any larger matchgroup that is NOT eliminated, except for outlier regions corresponding to Small Inversions/Duplications

				     Exception : For deletion with qry overlap exceeding smap_del_maxqryoverlap, don't trim the query overlap so deletion is suppressed ( and we can see why !)
				                 For translocations with qry overlap exceeding smap_translocation_ovlmx, don't trim the qry overlap so translocations is suppressed (and we can see why !)
						 For inversions with qry overlap exceeding smap_InvMaxOverlapSize, don't trim the qry overlap so inversion is suppressed (and we can see why !)

				                 Also if MATCHGROUP_PARTIAL_TRIM < MATCHGROUP_PARTIAL_TRIM_COMPLEX, do not trim query overlap for matchgroup pairs that do NOT result in any SV call
				  */

#define MATCHGROUP_PARTIAL_TRIM_COMPLEX 2 /* WAS151 1 */

#define MATCHGROUP_PAIR_FIX 1 /* NEW30 : Don't generate MG pairs (same qry and ref) if there is an intermediate MG that overlapps the gap of both ref and qry by at least 
				 smap_GapOverlapRatio x min(gap size, intermediate MG size) AND the intermediate MG overlaps the gap + end MGs by at least smap_GapOverlapRatio 
				 times the intermediate MG size AND the intermediate MG orientation matches the orientation of either of the outer matchgroups */

#define MATCHGROUP_PAIR_QFIX 1 /* NEW56 : temporary fix until dominant orientation based calling of neighboring inverted matchgroups is done : 
				  apply MATCHGROUP_PAIR_FIX only if all 3 MGs have the same orientation */

#define INVERTED_DUP_FIX 1 /* fix inversion duplicate filter to avoid eliminating inversions that are not definitely duplicates of each other */
#define TRANS_DUP_FIX 1 /* disables translocation duplicate filter for overlapped translocations */

// HERE HERE : need to suppress deletion call if most of the ref deleted region is actually aligned to the same qry via a different matchgroup
// HERE HERE : need to suppress insertion call if most of the qry inserted region is actually aligned to the same ref via a different matchgroup

#define INDEL_FILTER_FIX 1 /* NEW30 : include duplication SV in analysis of redundant indel with MGs 1-3 based on SV 1-2 and SV 2-3 being indel OR duplication */


#define QRY_PRIMARY_ORIENTATION 1 /* If query has MGs in both orientations for a given ref and the size of MGs for one orientation is at least smap_PrimaryOrientation (2.0) times larger
				       than the size for the other orientation then don't consider 2-MG indels that are not in the primary orientation */

#define SMALL_INVERSION_FIX 1 /* Try to trim inverted matchgroup on either end to remove palindromic extensions */
#define SMALL_DUPLICATION_FIX 1 /* Try to trim smaller matchgroup on either end to reduced query overlap with larger matchgroup */

#define INVERSION_OVERLAP_FIX 1 /* NEW8 : Don't call inversion if inversion matchgroup barely overlaps in query and mainly overlaps in reference */

#define INV_OVERLAP_FILTER_FIX 1 /* see OverlapFilterMG : also waives Qry Label overlap requirement if smaller MG is fully overlapped in both ref and qry by larger MG */

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/output_smap.cpp 11394 2020-07-31 17:18:50Z tanantharaman $");

static int xmapidMax= 0;

//usexmapentries are all the xmapEntries which pass the confidence threshold, and smapsize is the length
//find pairs which have the same query contig--store in indices{1,2}

//the old way is this:
//  for each entry, find the next closest entry only
//    this avoids unnecessary duplications like A->B, A->C, B->C. Instead, just do A->B, B->C.

//the new way is to generate all pairs and filter redundant calls (eg indels from A->C) later

// However this may not always work if the matchgroups are not all the same orientation or with the same ref contig, and hence the breakpoints correspond to different SV types
//      Simple fix with MATCHGROUP_PAIR_FIX : see above (also MATCHGROUP_PAIR_QFIX until dominant orientation is implemented correctly)

// NOTE : assumes usexmapentries[] are sorted in ascending order of qrycontigid, refcontigid, refstartpos, refendpos

//return len of indices arrays
int getSameQueryIndices(xmapEntry **usexmapentries, int smapsize, int*& indices1, int*& indices2, int indsize) {
  int indlen = 0;

  if(QRY_PRIMARY_ORIENTATION && smap_PrimaryOrientation <= 1.0)
    for(int i=0; i < smapsize; i++)
      usexmapentries[i]->dominant_orientation = 0;

  if(QRY_PRIMARY_ORIENTATION && smap_PrimaryOrientation > 1.0){/* compute dominant_orientation of query for all MGs, based on maximizing number of neighboring pairs of MGs on ref being in same
								  order as on query (not crossed) */
    for(int i=0; i < smapsize; i++){
      xmapEntry *pMGi = usexmapentries[i];
      if(i==0 || pMGi->qrycontigid != usexmapentries[i-1]->qrycontigid || pMGi->refcontigid != usexmapentries[i-1]->refcontigid){
	double forward = 0.0, backward = 0.0;

	// HERE HERE : should consider all pairs i,j and sum up crossed regions (in kb on qry or ref, whichever is less)

	int j = i;
	for(; j < smapsize; j++){
	  xmapEntry *pMGj = usexmapentries[j];
	  if(pMGi->qrycontigid != pMGj->qrycontigid || pMGi->refcontigid != pMGj->refcontigid)
	    break;
	  if(pMGj->refendpos <= pMGi->refendpos)
	    continue;
	  double minq1 = min(pMGi->qrystartpos, pMGi->qryendpos);
	  double maxq1 = max(pMGi->qrystartpos, pMGi->qryendpos);
	  double minq2 = min(pMGj->qrystartpos, pMGj->qryendpos);
	  double maxq2 = max(pMGj->qrystartpos, pMGj->qryendpos);

	  if(minq1 < minq2 && maxq1 < maxq2)
	    forward += min(maxq1-minq1, maxq2-minq2);
	  else if(minq1 > minq2 && maxq1 > maxq2)
	    backward += min(maxq1-minq1, maxq2-minq2);
	}
	int dominant = 0;
	if(forward > backward * smap_PrimaryOrientation)
	  dominant = 1;
	else if(backward > forward * smap_PrimaryOrientation)
	  dominant = 2;

	if(VERB>=2){
	  printf("usexmapentries[%d..%d] : dominant_orientation -> %d, forward= %0.3f, backward= %0.3f kb\n",i,j-1,dominant, forward*1e-3,backward*1e-3);
	  fflush(stdout);
	}

	for(int t = i; t < j; t++)
	  usexmapentries[t]->dominant_orientation = dominant;
      }
    }
  }

  /* HERE HERE HERE : should also compute local orientation for each MG pair based on flanking (surrounding) 2 MGs : should replace indices1[],indices2[] by array of class with indices1,indices2,orientation */

  for(int i=0; i < smapsize - 1; i++) { //first idx doesn't need to include last ele
    xmapEntry *pMGi = usexmapentries[i];

    for(int j= i+1; j < smapsize; j++) { //start at next ele
      xmapEntry *pMGj = usexmapentries[j];
      if( pMGi->qrycontigid != pMGj->qrycontigid )
	break;// since qrycontigid are in ascending order

      if(VERB>=3 && MATCHGROUP_PAIR_FIX && pMGi->refcontigid == 5 && pMGj->refcontigid == 19 && pMGi->qrycontigid== 3461 && pMGj->qrycontigid== 3461){
	double minq1 = min(pMGi->qrystartpos, pMGi->qryendpos);
	double maxq1 = max(pMGi->qrystartpos, pMGi->qryendpos);
	double minq2 = min(pMGj->qrystartpos, pMGj->qryendpos);
	double maxq2 = max(pMGj->qrystartpos, pMGj->qryendpos);
	printf("i=%d,j=%d: xmapids=%d,%d,fwd=%d,%d: qryi= %0.3f .. %0.3f, qryj= %0.3f .. %0.3f, refi= %0.3f .. %0.3f, refj= %0.3f .. %0.3f:smap_GapOverlapRatio= %0.4f\n",
	       i,j,pMGi->xmapid,pMGj->xmapid,pMGi->orientforward,pMGj->orientforward,
	       minq1*1e-3,maxq1*1e-3,minq2*1e-3,maxq2*1e-3,pMGi->refstartpos*1e-3,pMGi->refendpos*1e-3,pMGj->refstartpos*1e-3,pMGj->refendpos*1e-3,smap_GapOverlapRatio);
	fflush(stdout);
      }

      if(MATCHGROUP_PAIR_FIX && smap_GapOverlapRatio > 0.0 && pMGi->refcontigid == pMGj->refcontigid && pMGi->refendpos < pMGj->refstartpos && 
	 !(MATCHGROUP_PAIR_QFIX && pMGi->orientforward != pMGj->orientforward)){
	double minq1 = min(pMGi->qrystartpos, pMGi->qryendpos);
	double maxq1 = max(pMGi->qrystartpos, pMGi->qryendpos);
	double minq2 = min(pMGj->qrystartpos, pMGj->qryendpos);
	double maxq2 = max(pMGj->qrystartpos, pMGj->qryendpos);

	if(DEBUG>=2 && !((!overlap(minq1,maxq1,minq2,maxq2)) == (min(maxq1,maxq2) < max(minq1,minq2)))){
	  printf("i=%d,j=%d: xmapids=%d,%d: qryi= %0.3f .. %0.3f, qryj= %0.3f .. %0.3f, refi= %0.3f .. %0.3f, refj= %0.3f .. %0.3f:smap_GapOverlapRatio= %0.4f\n",
		 i,j,pMGi->xmapid,pMGj->xmapid,minq1*1e-3,maxq1*1e-3,minq2*1e-3,maxq2*1e-3,pMGi->refstartpos*1e-3,pMGi->refendpos*1e-3,pMGj->refstartpos*1e-3,pMGj->refendpos*1e-3,smap_GapOverlapRatio);
	  printf("\toverlap(minq1,maxq1,minq2,maxq2)= %d\n",overlap(minq1,maxq1,minq2,maxq2));
	  fflush(stdout);
	  assert((!overlap(minq1,maxq1,minq2,maxq2)) == (min(maxq1,maxq2) <= max(minq1,minq2)));
	}

	if(min(maxq1,maxq2) < max(minq1,minq2)){ // WAS40 !overlap(minq1,maxq1,minq2,maxq2)

	  double refstart = pMGi->refstartpos, qrystart;
	  double refgapstart = pMGi->refendpos, qrygapstart;
	  double refgapstop = pMGj->refstartpos, qrygapstop;
	  double refstop = pMGj->refendpos, qrystop;

	  //	  bool query_forward;

	  if( maxq1 < minq2){
	    //	    query_forward = true;
	    qrystart = minq1;
	    qrygapstart = maxq1;
	    qrygapstop = minq2;
	    qrystop = maxq2;
	  } else {
	    if(DEBUG) assert(maxq2 < minq1);
	    //	    query_forward = false;
	    qrystart = minq2;
	    qrygapstart = maxq2;
	    qrygapstop = minq1;
	    qrystop = maxq1;// WAS67 maxq2
	  }

	  if(VERB>=3 && MATCHGROUP_PAIR_FIX && pMGi->refcontigid == 17 && pMGj->refcontigid == 17 && pMGi->qrycontigid==531 && pMGj->qrycontigid==531){
	    printf("\t qrystart= %0.3f, qrystop= %0.3f, qrygapstart= %0.3f, qrygapstop= %0.3f\n", qrystart,qrystop, qrygapstart,qrygapstop);
	    fflush(stdout);
	  }

	  int k = i+1;
	  for(; k < j; k++){
	    xmapEntry *pMGk = usexmapentries[k];
	    if(VERB>=3 && MATCHGROUP_PAIR_FIX && pMGi->refcontigid == 17 && pMGj->refcontigid == 17 && pMGi->qrycontigid==531 && pMGj->qrycontigid==531){
	      double minq3 = min(pMGk->qrystartpos, pMGk->qryendpos);
	      double maxq3 = max(pMGk->qrystartpos, pMGk->qryendpos);
	      printf("\t k=%d/%d: xmapids=%d,%d,%d,qryk= %0.3f .. %0.3f, refk= %0.3f .. %0.3f kb, fwdi=%d,fwdj=%d,fwdk=%d(refidk=%lld,qryidk=%lld)\n",
		     k, j, pMGi->xmapid,pMGj->xmapid,pMGk->xmapid,minq3*1e-3,maxq3*1e-3,pMGk->refstartpos*1e-3,pMGk->refendpos*1e-3, pMGi->orientforward,pMGj->orientforward,pMGk->orientforward,pMGk->refcontigid,pMGk->qrycontigid);
	      fflush(stdout);
	    }
	    if(pMGk->refcontigid != pMGi->refcontigid || pMGk->qrycontigid != pMGi->qrycontigid){
	      k = smapsize;
	      break;
	    }
	    if(MATCHGROUP_PAIR_QFIX  ? 
	       (pMGk->orientforward != pMGi->orientforward || pMGk->orientforward != pMGj->orientforward) :
	       (pMGk->orientforward != pMGi->orientforward && pMGk->orientforward != pMGj->orientforward))
	      continue;
	    double minq3 = min(pMGk->qrystartpos, pMGk->qryendpos);
	    double maxq3 = max(pMGk->qrystartpos, pMGk->qryendpos);
	    if(VERB>=3 && MATCHGROUP_PAIR_FIX && pMGi->refcontigid == 17 && pMGj->refcontigid == 17 && pMGi->qrycontigid==531 && pMGj->qrycontigid==531){
	      printf("i=%d,j=%d,k=%d:xmapids=%d,%d,%d, qryi= %0.3f .. %0.3f, qryj= %0.3f .. %0.3f, qryk= %0.3f .. %0.3f, refi= %0.3f .. %0.3f, refj= %0.3f .. %0.3f, refk= %0.3f .. %0.3f kb\n",
		     i,j,k,pMGi->xmapid,pMGj->xmapid,pMGk->xmapid,minq1*1e-3,maxq1*1e-3,minq2*1e-3,maxq2*1e-3,minq3*1e-3,maxq3*1e-3,
		     pMGi->refstartpos*1e-3,pMGi->refendpos*1e-3,pMGj->refstartpos*1e-3,pMGj->refendpos*1e-3,pMGk->refstartpos*1e-3,pMGk->refendpos*1e-3);
	      printf("\toverlapLength(qrygapij,qryk)=%0.3f,overlapLength(qryij,qryk)=%0.3f,overlapLength(refgapij,refk)=%0.3f,overlapLength(refij,refk)=%0.3f\n",
		     overlapLength(qrygapstart,qrygapstop, minq3, maxq3)*1e-3,overlapLength(qrystart, qrystop, minq3, maxq3)*1e-3,
		     overlapLength(refgapstart,refgapstop, pMGk->refstartpos, pMGk->refendpos)*1e-3,overlapLength(refstart,refstop, pMGk->refstartpos, pMGk->refendpos));
	      fflush(stdout);
	    }
	    if(overlapLength(qrygapstart,qrygapstop, minq3, maxq3) >= min(maxq3 - minq3, qrygapstop - qrygapstart) * smap_GapOverlapRatio &&
	       overlapLength(qrystart, qrystop, minq3, maxq3) >= (maxq3 - minq3) * smap_GapOverlapRatio && 
	       overlapLength(refgapstart,refgapstop, pMGk->refstartpos, pMGk->refendpos) >= min(pMGk->refendpos - pMGk->refstartpos, refgapstop - refgapstart) * smap_GapOverlapRatio &&
	       overlapLength(refstart,refstop, pMGk->refstartpos, pMGk->refendpos) >= (pMGk->refendpos - pMGk->refstartpos) * smap_GapOverlapRatio)
	      break;
	  }
	  if(k < j){
	    if(VERB>=2 && pMGi->refcontigid==17 && pMGj->refcontigid==17 && pMGi->qrycontigid==531 && pMGj->qrycontigid==531){
	      xmapEntry *pMGk = usexmapentries[k];
	      double minq3 = min(pMGk->qrystartpos, pMGk->qryendpos);
	      double maxq3 = max(pMGk->qrystartpos, pMGk->qryendpos);
	      printf("Skipping MG pair i=%d, j=%d (xmapid=%d,%d) due to MG k=%d: qryi= %0.3f .. %0.3f, qryj= %0.3f .. %0.3f, qryk= %0.3f .. %0.3f, refi= %0.3f .. %0.3f, refj= %0.3f .. %0.3f, refk= %0.3f .. %0.3f kb,fwd=%d,%d,%d\n",
		     i,j,usexmapentries[i]->xmapid,usexmapentries[j]->xmapid,k,minq1*1e-3,maxq1*1e-3,minq2*1e-3,maxq2*1e-3,minq3*1e-3,maxq3*1e-3,pMGi->refstartpos*1e-3,pMGi->refendpos*1e-3,pMGj->refstartpos*1e-3,pMGj->refendpos*1e-3,
		     pMGk->refstartpos*1e-3,pMGk->refendpos*1e-3,pMGi->orientforward,pMGk->orientforward,pMGj->orientforward);
	      fflush(stdout);
	    }
	    continue;
	  }

	  k = i-1;
	  for(;k >= 0; k--){
	    xmapEntry *pMGk = usexmapentries[k];
	    if(VERB>=3 && MATCHGROUP_PAIR_FIX && pMGi->refcontigid == 17 && pMGj->refcontigid == 17 && pMGi->qrycontigid==531 && pMGj->qrycontigid==531){
	      double minq3 = min(pMGk->qrystartpos, pMGk->qryendpos);
	      double maxq3 = max(pMGk->qrystartpos, pMGk->qryendpos);
	      printf("\t k=%d: xmapids=%d,%d,%d,qryk= %0.3f .. %0.3f, refk= %0.3f .. %0.3f kb, fwdi=%d,fwdj=%d,fwdk=%d(refidk=%lld,qryidk=%lld)\n",
		     k, pMGi->xmapid,pMGj->xmapid,pMGk->xmapid,minq3*1e-3,maxq3*1e-3,pMGk->refstartpos*1e-3,pMGk->refendpos*1e-3, pMGi->orientforward,pMGj->orientforward,pMGk->orientforward,pMGk->refcontigid,pMGk->qrycontigid);
	      fflush(stdout);
	    }
	    if(pMGk->refcontigid != pMGi->refcontigid || pMGk->qrycontigid != pMGi->qrycontigid){
	      k = -1;
	      break;
	    }
	    //	    if(pMGk->refendpos <= refgapstart)
	    //	      continue;
	    if(MATCHGROUP_PAIR_QFIX  ? 
	       (pMGk->orientforward != pMGi->orientforward || pMGk->orientforward != pMGj->orientforward) :
	       (pMGk->orientforward != pMGi->orientforward && pMGk->orientforward != pMGj->orientforward))
	      continue;
	    double maxq3 = max(pMGk->qrystartpos, pMGk->qryendpos);
	    double minq3 = min(pMGk->qrystartpos, pMGk->qryendpos);
	    if(VERB>=2 && MATCHGROUP_PAIR_FIX && pMGi->refcontigid == 17 && pMGj->refcontigid == 17 && pMGi->qrycontigid==531 && pMGj->qrycontigid==531){
	      printf("i=%d,j=%d,k=%d:xmapids=%d,%d,%d: qryi= %0.3f .. %0.3f, qryj= %0.3f .. %0.3f, qryk= %0.3f .. %0.3f, refi= %0.3f .. %0.3f, refj= %0.3f .. %0.3f, refk= %0.3f .. %0.3f kb\n",
		     i,j,k,pMGi->xmapid,pMGj->xmapid,pMGk->xmapid,minq1*1e-3,maxq1*1e-3,minq2*1e-3,maxq2*1e-3,minq3*1e-3,maxq3*1e-3,pMGi->refstartpos*1e-3,pMGi->refendpos*1e-3,pMGj->refstartpos*1e-3,pMGj->refendpos*1e-3,
		     pMGk->refstartpos*1e-3,pMGk->refendpos*1e-3);
	      printf("\toverlapLength(qrygapij,qryk)=%0.3f,overlapLength(qryij,qryk)=%0.3f,overlapLength(refgapij,refk)=%0.3f,overlapLength(refij,refk)=%0.3f\n",
		     overlapLength(qrygapstart,qrygapstop, minq3, maxq3),overlapLength(qrystart, qrystop, minq3, maxq3),overlapLength(refgapstart,refgapstop, pMGk->refstartpos, pMGk->refendpos),
		     overlapLength(refstart,refstop, pMGk->refstartpos, pMGk->refendpos));
	      fflush(stdout);
	    }
	    if(overlapLength(qrygapstart,qrygapstop, minq3, maxq3) >= min(maxq3 - minq3, qrygapstop - qrygapstart) * smap_GapOverlapRatio &&
	       overlapLength(qrystart, qrystop, minq3, maxq3) >= (maxq3 - minq3) * smap_GapOverlapRatio && 
	       overlapLength(refgapstart,refgapstop, pMGk->refstartpos, pMGk->refendpos) >= min(pMGk->refendpos - pMGk->refstartpos, refgapstop - refgapstart) * smap_GapOverlapRatio &&
	       overlapLength(refstart,refstop, pMGk->refstartpos, pMGk->refendpos) >= (pMGk->refendpos - pMGk->refstartpos) * smap_GapOverlapRatio)
	      break;
	  }
          if(k >= 0){
	    if(VERB>=3 && pMGi->refcontigid == 17 && pMGj->refcontigid == 17 && pMGi->qrycontigid==531 && pMGj->qrycontigid==531){
	      xmapEntry *pMGk = usexmapentries[k];
	      double minq3 = min(pMGk->qrystartpos, pMGk->qryendpos);
	      double maxq3 = max(pMGk->qrystartpos, pMGk->qryendpos);
	      printf("Skipping MG pair i=%d,j=%d(xmapid=%d,%d) due to MG k=%d: qryi= %0.3f .. %0.3f, qryj= %0.3f .. %0.3f, qryk= %0.3f .. %0.3f, refi= %0.3f .. %0.3f, refj= %0.3f .. %0.3f, refk= %0.3f .. %0.3f kb, fwd=%d,%d,%d\n",
		     i,j,usexmapentries[i]->xmapid,usexmapentries[j]->xmapid,k,minq1*1e-3,maxq1*1e-3,minq2*1e-3,maxq2*1e-3,minq3*1e-3,maxq3*1e-3,pMGi->refstartpos*1e-3,pMGi->refendpos*1e-3,pMGj->refstartpos*1e-3,pMGj->refendpos*1e-3,
		     pMGk->refstartpos*1e-3,pMGk->refendpos*1e-3,pMGi->orientforward,pMGk->orientforward,pMGj->orientforward);
	      fflush(stdout);
	    }
	    continue;
	  }
	}
      }

      if( indlen >= indsize ) { // need to increase size
	growArray(indices1, indsize); //double size (but indsize is used on next line, so don't change it yet)
	indsize = growArray(indices2, indsize); //double size
      }
      indices1[indlen] = i;
      indices2[indlen] = j;	
      indlen++;
      if( ! smap_xmap_allpairs ) //don't break to keep all pairs (goes from N-1 to N(N-1)/2)
	  break; //for each entry, find the next closest entry only
    }
  }

  if(VERB>=2){
    printf("At end of getSameQueryIndices:\n");
    for(int t = 0; t < indlen; t++)
      printf("  t=%d/%d:indices1[t]=%d,indices2[t]=%d (xmapid=%d,%d)\n",t,indlen,indices1[t],indices2[t],usexmapentries[indices1[t]]->xmapid,usexmapentries[indices2[t]]->xmapid);
    fflush(stdout);
  }

  return indlen;
}


//struct for dealing with end-type SVs--keep track of number of un-aligned labels and length at start and end of a query contig
//then make an array of these for all query contigs
class endalign {
public:
  endalign(long long, int, int, double, double, long long, double, bool, int, int, int);
  endalign( const endalign& );
  endalign operator=( const endalign& orig );
  //data members
  long long contigid;
  int mapid;// index into Gmap[] 
  int startunalignlab; // number of unaligned labels
  double startunalignlen; // unaligned length (in kb)
  double contiglen; //contig length in kb
  int  contignsite; //contig n sites
  bool startsvused;
  bool isstart;
  int refstartlen; //length of all arrays below
  long long* refstartids;
  double* refstartpos;
  int* refidx; //index for position above
  int* xids; //indices in xmap of matches store in ref arrays above

  void addRefStart(double, long long, int, int);
  void replaceRefStart(double, long long, int, int);
  void printData(int);
  endalign() { 
    refstartpos = NULL;
    refstartids = NULL;
    refidx = xids = NULL;
  };
  ~endalign() {
    delete [] refstartids;
    delete [] refstartpos;
    delete [] refidx;
    delete [] xids;
  }
};

endalign::endalign(long long cid, int mid, int startu=int(1e4), double startl=1e6, double clen=0, long long startid=0, double startp=0, bool isrt=false, int cns=0, int ridx=0, int xi=0) {
  contigid = cid;
  mapid = mid;
  startunalignlab = startu;
  startunalignlen = startl;
  contiglen       = clen;
  contignsite     = cns;
  startsvused     = false;
  isstart         = isrt;
  if( startid > 0 ) { //if startid, have single element
    refstartids = new long long[1];
    refstartids[0] = startid;
    refstartpos = new double[1];
    refstartpos[0] = startp;
    refidx      = new int[1];
    refidx[0]   = ridx;
    xids        = new int[1];
    xids[0]     = xi;
  }
  else {
    refstartids = 0;
    refstartpos = 0;
    refidx      = 0;
    xids        = 0;
  }
  refstartlen     = (startid > 0 ? 1 : 0); //if startid, have single element
}

//Need this bc of array data members.
endalign::endalign( const endalign& orig ) {
  contigid = orig.contigid;
  mapid = orig.mapid;
  startunalignlab = orig.startunalignlab;
  startunalignlen = orig.startunalignlen;
  contiglen       = orig.contiglen      ;
  contignsite     = orig.contignsite    ;
  startsvused     = orig.startsvused    ;
  isstart         = orig.isstart;
  refstartlen     = orig.refstartlen;
  if( refstartlen > 0 ) {
    refstartids     = new long long[refstartlen];
    refstartpos     = new double[refstartlen];
    refidx          = new int   [refstartlen];
    xids            = new int   [refstartlen];
    for( int i=0; i<refstartlen; i++ ) {
      refstartids[i] = orig.refstartids[i];
      refstartpos[i] = orig.refstartpos[i];
      refidx[i]      = orig.refidx[i];
      xids[i]        = orig.xids[i];
    }
  }
  else {
    refstartids = 0;
    refstartpos = 0;
    refidx = 0;
    xids = 0;
  }
  //printf("Copy constructor\n");
  //this->printData(-1);
}

//Note: need this when you assign to memory already allocated
//"Any time you need to write your own custom copy constructor, you also need to write a custom assignment operator."
endalign endalign::operator=( const endalign& orig ) {
  contigid = orig.contigid;
  mapid = orig.mapid;
  startunalignlab = orig.startunalignlab;
  startunalignlen = orig.startunalignlen ;
  contiglen       = orig.contiglen       ;
  contignsite     = orig.contignsite    ;
  startsvused     = orig.startsvused     ;
  isstart         = orig.isstart;
  refstartlen     = orig.refstartlen;
  if( refstartlen > 0 ) {
    refstartids     = new long long[refstartlen];
    refstartpos     = new double[refstartlen];
    refidx          = new int   [refstartlen];
    xids            = new int   [refstartlen];
    for( int i=0; i<refstartlen; i++ ) {
      refstartids[i] = orig.refstartids[i];
      refstartpos[i] = orig.refstartpos[i];
      refidx[i]      = orig.refidx[i];
      xids[i]        = orig.xids[i];
    }
  }
  else {
    refstartids = 0;
    refstartpos = 0;
    refidx = 0;
    xids = 0;
  }
  return *this;
}

//if new refstart is found, replace existing data with new
void endalign::replaceRefStart(double refstart, long long refcontigid, int ridx, int xi) {
  if( refstartlen > 0 ) {
    delete [] refstartids;
    delete [] refstartpos;
    delete [] refidx;
    delete [] xids;
  }
  refstartlen = 1;
  refstartids = new long long[1];
  refstartids[0] = refcontigid;
  refstartpos = new double[1];
  refstartpos[0] = refstart;
  refidx = new int[1];
  refidx[0] = ridx;
  xids = new int[1];
  xids[0] = xi;
}

//grow the refstartids and refstartpos lists--these will only have a few entries, so just grow by 1 each time
void endalign::addRefStart(double refstart, long long refcontigid, int ridx, int xi) {
  long long* newstartids = new long long[refstartlen+1];
  double* newstartpos = new double[refstartlen+1];
  int*    newidx      = new int[refstartlen+1];
  int*    newxi       = new int[refstartlen+1];
  for(int i=0; i<refstartlen; i++) {
    newstartids[i] = refstartids[i];
    newstartpos[i] = refstartpos[i];
    newidx[i]      = refidx[i];
    newxi[i]       = xids[i];
  }
  newstartids[refstartlen] = refcontigid;
  newstartpos[refstartlen] = refstart;
  newidx[refstartlen] = ridx;
  newxi[refstartlen]  = xi;
  if( refstartlen > 0 ) { //if previously allocated, delete
    delete [] refstartids;
    delete [] refstartpos;
    delete [] refidx;
    delete [] xids;
  }
  refstartids = newstartids;
  refstartpos = newstartpos;
  refidx = newidx;
  xids   = newxi;
  refstartlen += 1;
}


void endalign::printData(int idx=-1) {
  printf("%2i id %4lld  labs: %3i  lens: %7.3f;  %i\n", idx, contigid, startunalignlab, startunalignlen, refstartlen );
  for(int k=0; k<refstartlen; k++) 
    printf("  %s: %2lld  %.3f\n", (isstart ? "start" : "end"), refstartids[k], refstartpos[k]);
  fflush(stdout);
}


/* Hashtable to quickly map query id (64 bits) to index into endalign[0..nummaps*2-1] */

class Cid2indexHash {/* NOTE : struct is 64 bytes long to fit in 1 or 2 cache lines */
public:
  long long id[5];
  int endindex[5];
  int cnt;

  inline void init(int C,long long Id, int EndIndex){
    if(DEBUG>=2) assert(0 <= C && C < 5);
    id[C] = Id;
    endindex[C] = EndIndex;
    cnt = C+1;
  };
};

static Cid2indexHash *HashTable = 0;
static int *HashIndex = 0;
static int Next = 5, HashMax = 0;/* Next is last used entry in HashTable[6..HashMax-1] */

#define H_BITS 20 // static size of hashtable index is (1 << H_BITS)

static inline int HashKey(long long id)
{
  return (id ^ (id >> H_BITS)/* ^ (id >> (2*H_BITS))*/) & MASK(H_BITS);
}

static inline int HashNext()
{
  if(++Next >= HashMax){
    int origHashMax = HashMax;
    HashMax = (HashMax * 3) / 2;
    Cid2indexHash *newHashTable = (Cid2indexHash *)malloc(HashMax * sizeof(Cid2indexHash));
    if(!newHashTable){
      printf("HashNext:malloc(%llu) failed\n", (unsigned long long)(HashMax * sizeof(Cid2indexHash)));
      exit(1);
    }
    memcpy(&newHashTable[6],&HashTable[6], (origHashMax-6) * sizeof(Cid2indexHash));
    free(HashTable);
    HashTable = newHashTable;
  }
  return Next;
}

static void HashInit()
{
  if(HashIndex)
    free(HashIndex);
  if(!(HashIndex = (int *)calloc(1 << H_BITS, sizeof(int)))){
    printf("HashInit:calloc(%llu) failed\n",(unsigned long long)((1 << H_BITS)*sizeof(int)));
    exit(1);
  }

  if(HashTable)
    free(HashTable);
  HashMax = (1 << H_BITS);
  if(!(HashTable = (Cid2indexHash *)malloc(HashMax*sizeof(Cid2indexHash)))){
    printf("HashInit:malloc(%llu) failed\n", (unsigned long long)(HashMax * sizeof(Cid2indexHash)));
    exit(1);
  }
  Next = 5;
}

/* NOTE : it is illegal to insert the same id twice */
static void HashInsert(long long id, int EndIndex)
{
  int key = HashKey(id);
  int index = HashIndex[key];
  if(!index){
    HashIndex[key] = index = HashNext();
    if(DEBUG>=2) assert(6 <= index && index <= Next);
    Cid2indexHash *p = &HashTable[index];
    p->init(0,id,EndIndex);
  } else {
    if(DEBUG>=2) assert(6 <= index && index <= Next);
    Cid2indexHash *p = &HashTable[index];
    int cnt = p->cnt;
    if(DEBUG>=2)/* check that id is not already in HashTable */
      for(int i = min(5,cnt); --i >= 0;)
	assert(p->id[i] != id);
    while(cnt > 5){/* overflowed */
      if(DEBUG>=2) assert(cnt <= Next);
      p = &HashTable[index = cnt];
      cnt = p->cnt;
      if(DEBUG>=2)/* check that id is not already in HashTable */
	for(int i = min(5,cnt); --i >= 0;)
	  assert(p->id[i] != id);
    }
    if(DEBUG>=2) assert(p == &HashTable[index] && cnt == p->cnt);

    if(cnt >= 5){/* full */
      int nindex = HashNext();
      HashTable[index].cnt = nindex;// HashTable may have been reallocated by HashNext(), invalidating previous value of p
      p = &HashTable[nindex];
      p->init(0,id,EndIndex);
    } else
      p->init(cnt,id,EndIndex);
  }
}

/* returns EndIndex or -1 if not found */
int HashFindID(long long id)
{
  int key = HashKey(id);
  int index = HashIndex[key];
  if(!index)
    return -1;
  if(DEBUG>=2) assert(6 <= index && index <= Next);
  Cid2indexHash *p = &HashTable[index];
  int cnt = p->cnt;
  for(int i = min(5,cnt); --i >= 0;)
    if(p->id[i] == id)
      return p->endindex[i];
  while(cnt > 5){/* overflowed */
    if(DEBUG>=2) assert(cnt <= Next);
    p = &HashTable[cnt];
    cnt = p->cnt;
    for(int i = min(5,cnt); --i >= 0;)
      if(p->id[i] == id)
	return p->endindex[i];
  }
  return -1;
}

//when looping over contigs, either update the endalign object for that contig, 
//or add a new one if it doesn't exist yet
//return next index to be used in next call of this fn (size of ea_array is fixed)
int mergeEndaligns(endalign*& ea_array, int easize, xmapEntry* xmapentry) {

  long long cid = xmapentry->qrycontigid;
  int mapid = xmapentry->align->mapid2;
  int startu = xmapentry->startunalign; //number of unaligned labels at start of query (id == cid)
  int endu   = xmapentry->endunalign; //number of unaligned labels at end of query (id == cid)
  if(DEBUG>=1+RELEASE/* HERE HERE >=2 */) assert(startu == min(xmapentry->qrystartidx, xmapentry->qryendidx) - 1);
  if(DEBUG>=1+RELEASE/* HERE HERE >=2 */) assert(endu == xmapentry->querycontignsite  - max(xmapentry->qrystartidx, xmapentry->qryendidx));
  double startl = min(xmapentry->qrystartpos, xmapentry->qryendpos) * 1e-3;
  double endl   = xmapentry->querycontiglen - max(xmapentry->qrystartpos, xmapentry->qryendpos) * 1e-3;
  double refstart = (xmapentry->qrystartpos < xmapentry->qryendpos ? xmapentry->refstartpos : xmapentry->refendpos);
  double refend   = (xmapentry->qrystartpos > xmapentry->qryendpos ? xmapentry->refstartpos : xmapentry->refendpos);
  int    refsidx  = (xmapentry->qrystartpos < xmapentry->qryendpos ? xmapentry->refstartidx : xmapentry->refendidx);
  int    refeidx  = (xmapentry->qrystartpos > xmapentry->qryendpos ? xmapentry->refstartidx : xmapentry->refendidx);

  //first, check if the contig id is already there
  // NOTE : since this function is called with qrycontigid in ascending order we only need to check the last 2 entries
  bool havecontig = false; //set to true if endalign already exists for this contig
  for(int i= /* WAS 0 */ max(0, easize - 2); i < easize; i++) {
    if( ea_array[i].contigid == cid ) {
      if(DEBUG>=1+RELEASE) assert(ea_array[i].mapid == mapid);/* If cid matches, so should mapid */
      //      if(DEBUG) assert(i >= easize - 2);
      havecontig = true;

      //check for update for each start/end
      if( startu < ea_array[i].startunalignlab && ea_array[i].isstart ) {
	if(DEBUG) assert( startl < ea_array[i].startunalignlen ); //if n lab is <, length must be also
	ea_array[i].startunalignlab = startu; 
	ea_array[i].startunalignlen = startl;
	ea_array[i].replaceRefStart(refstart, xmapentry->refcontigid, refsidx, xmapentry->xmapid); //remove old, replace with current
      }
      else if( startu == ea_array[i].startunalignlab && ea_array[i].isstart ) { //same end occurs twice--add ref info
	if(DEBUG) assert( startl == ea_array[i].startunalignlen );
	ea_array[i].addRefStart(refstart, xmapentry->refcontigid, refsidx, xmapentry->xmapid);
      }
      if( endu < ea_array[i].startunalignlab && !ea_array[i].isstart ) {
	if(DEBUG) assert( endl < ea_array[i].startunalignlen ); //if n lab is <, length must be also
	ea_array[i].startunalignlab   = endu; 
	ea_array[i].startunalignlen   = endl;
	ea_array[i].replaceRefStart(refend, xmapentry->refcontigid, refeidx, xmapentry->xmapid);
      }
      else if( endu == ea_array[i].startunalignlab && !ea_array[i].isstart ) {
	if(DEBUG) assert( endl == ea_array[i].startunalignlen );
	ea_array[i].addRefStart(refend, xmapentry->refcontigid, refeidx, xmapentry->xmapid);
      }
      //can't return bc need to check start and end
    }
  }
  if( havecontig )
    return easize;

  // all allocation is done ahead of time: still need to insert at end
  //if not yet returned, it's a new contig--allocate space for it

  HashInsert(cid, easize);// remember mapping from cid to index into ea_array[]

  ea_array[easize  ] = endalign(cid, mapid, startu, startl, xmapentry->querycontiglen, xmapentry->refcontigid, refstart, true, xmapentry->querycontignsite, refsidx, xmapentry->xmapid);
  ea_array[easize+1] = endalign(cid, mapid, endu  , endl  , xmapentry->querycontiglen, xmapentry->refcontigid, refend , false, xmapentry->querycontignsite, refeidx, xmapentry->xmapid);

  return easize+2;
}


//given a query contig id, check if it has a valid 'end' sv at the beginning, 
// and if so, write it to file and set its startsvused data member false
//bool checkstart is needed so this fn can work for 'end' svs at the beginning and the end of contigs
//Thresholds are: smap_min_end_nlabels and smap_min_end_lenkb -- see parameters.cpp
//  If dowrite==false, don't write anything to smap (fp), just mark corresponding xmapentry->output to true
//  xmapidMap[xmapid = 1 .. xmapidMax] is map from xmapid to xmapidfilt (NOT valid if dowrite==false)
//  xmapindexMap[xmapid = 1 .. xmapidMax] is map from xmapid to xmapentries[] index (used for debugging)
//  global array xmapentries[i=0..numxmapentries-1] is the complete list of xmap entries (used for debugging)
//return number of sv added so that nsv can be incremented

int checkEndSV(FILE *fp, int nsv, long long qrycontigid, endalign* endaligns, int easize, bool checkstart, bool dowrite, int* xmapidMap, int *xmapindexMap) {

  //first, find the ele of endaligns array which corresponds to qrycontigid
  endalign* enda = 0; //need to use pointer because need the startsvused assignment to carry over to next use

  /* Use Hashtable to find first of two entries in endaligns[] with contigid == qrycontigid */
  int endindex = HashFindID(qrycontigid);// returns -1 if qrycontigid was not found
  if(endindex >= 0){
    if(DEBUG) assert(endaligns[endindex].contigid == qrycontigid && endaligns[endindex+1].contigid == qrycontigid);

    //  for(int i= 0; i < easize; i++) {
    for(int i= endindex; i <= endindex + 1; i++) {
      if(endaligns[i].contigid == qrycontigid && checkstart == endaligns[i].isstart){
	enda = &endaligns[i];
	break;
      }
    }
  }

  //if no element of endaligns matches this qrycontigid, below makes no sense
  //assert or return? May as well just return, though I think this should never happen--print a warning too
  if(enda /* WAS enda->contigid */ == 0) {
    printf("Warning in output_smap, checkEndSV: qry contigid %lld not found in endaligns array\n", qrycontigid);
    fflush(stdout);
    return 0;
  }

  if( enda->startunalignlab >= smap_min_end_nlabels && //doc says 'at least this', which is >=
      enda->startunalignlen > smap_min_end_lenkb &&
      !enda->startsvused) { //make sure you don't report the same sv twice
    //valid end sv at start of contig
    enda->startsvused = dowrite; //should be redundant
    //because these are all at the beginning, use '20' as the qrystart
    //and the end is the start unalign len, but need to convert from kb to b
    assert(enda->refstartlen > 0);
    for( int i=0; i < enda->refstartlen; i++ ) {
      int xmapid = enda->xids[i];
      if(DEBUG) assert(1 <= xmapid && xmapid <= xmapidMax);
      int index = xmapindexMap[xmapid];
      if(DEBUG) assert(0 <= index && index < numxmapentries);      
      
      // NEW90 : filter out if confidence is below -T or failing smap_confbylen or smap_ConfByInterval
      xmapEntry *entry = xmapentries[index];
      double Len = min( fabs(entry->qrystartpos - entry->qryendpos),
			fabs(entry->refstartpos - entry->refendpos) )  / 1000.;
      if(entry->confidence < LogPvThreshold || entry->confidence < smap_confbylen * Len || entry->confidence < smap_ConfByInterval * (entry->align->numpairs - 1.0))
	continue;

      if(!dowrite){
	xmapentries[index]->output = true;
	continue;
      }
      if(DEBUG) assert(xmapentries[index]->output == true);

      int xmapidfilt = xmapidMap[xmapid];
      if(DEBUG && !(xmapidfilt > 0)){
	printf("checkEndSV:i=%d,enda->xids[i]=xmapid=%d,xmapidMax=%d,xmapidMap[xmapid]=xmapidfilt=%d,xmapindexMap[xmapid]=index=%d,numxmapentries=%d\n", i,xmapid,xmapidMax,xmapidfilt,index,numxmapentries);
	printf("\t xmapentries[index]=(xmapid=%d,xmapidfilt=%d,output=%d)\n",xmapentries[index]->xmapid,xmapentries[index]->xmapidfilt,xmapentries[index]->output ? 1 : 0);
	fflush(stdout);
	assert(xmapidfilt > 0);
      }

      if( enda->isstart )
	writeToSmap(fp, nsv+i, qrycontigid, enda->refstartids[i], -1, 20., enda->startunalignlen*1e3, -1., enda->refstartpos[i], "end", 0, xmapidfilt, -1, 1, enda->startunalignlab+1, -1, enda->refidx[i], -1., -1., -1., -1., 0, -1, -1, -1,-1.0, -1.0, -1.0, -1.0, NULL);
      else
	writeToSmap(fp, nsv+i, qrycontigid, enda->refstartids[i], -1, (enda->contiglen-enda->startunalignlen)*1e3, enda->contiglen*1e3, enda->refstartpos[i], -1., "end", xmapidfilt, 0, -1, enda->contignsite - enda->startunalignlab, enda->contignsite, enda->refidx[i], -1, -1., -1., -1., -1.,  0, -1, -1, -1, -1.0, -1.0, -1.0, -1.0, NULL);
    }
    return enda->refstartlen;
  }
  return 0;
} //end checkEndSV



//check for duplicate inversions and translocations
// these are similar to duplicate indels in the sense that two matchgroups are each called to make an event with a third while only one of them makes sense.
// Simplest example: 1-2 is indel, 1-3 is duplicate inversion, 2-3 is true inversion 
//
// For inversions, the common MG (3 in above example) must be the inverted (or non-inverted) MG for both inversions, otherwise neither is marked as duplicate
//
// For translocations the orientation of matchgroups does NOT matter
//
//see comments below for details
//also note that because the svs are only created on the same query contig, their xmapentries
// must also be on the same query contig, so no change is necessary for multiple contigs

void checkDuplicates( structuralVariation*& sv_array, int numsvIndel, int indiceslen, xmapEntry **usexmapentries, sv_type_t exclude_type, int verbose=1) {
  if(VERB>=2){
    printf("checkDuplicates:type= %d,verbose=%d\n",exclude_type,verbose);
    fflush(stdout);
  }

  //these types are the only ones supported
  if(DEBUG) assert(exclude_type == intrachr_translocation ||
		   exclude_type == interchr_translocation ||
		   exclude_type == inversion_type);

  if(DEBUG>=2)
    for(int i = 0; i < numsvIndel; i++)
      assert(sv_array[i].type_enum != exclude_type);

  for(int i= numsvIndel; i < indiceslen; i++) {
    //printf("checkDuplicates: type=%d\n", sv_array[i].type_enum); //debug
    if (sv_array[i].type_enum != exclude_type )
      continue;
    if (TRANS_DUP_FIX && sv_array[i].type_enum != inversion_type && sv_array[i].trans_overlapped)
      continue;

    int qrycontigid1 = usexmapentries[sv_array[i].indices[0]]->qrycontigid;

    for(int j= i+1; j < indiceslen; j++) { //only need to compare to subsequent SVs
      if( sv_array[j].type_enum != exclude_type || sv_array[j].duplicate )
	continue;
      if (TRANS_DUP_FIX && sv_array[i].type_enum != inversion_type && sv_array[i].trans_overlapped)
	continue;
      int qrycontigid2 = usexmapentries[sv_array[j].indices[0]]->qrycontigid;
      if(qrycontigid2 != qrycontigid1)
	break;// NEW : rely on qrycontig being in ascending order, hence no need to check any larger j values

      //      const char* common = "No common xmapEntry";
      bool common = false;
      int dupe = -1;

      // below rule should only apply when common matchgroup is inverted OR if non-common matchgroup is further on query but closer on ref (see INVERTED_DUP_FIX)

      //check if the same xmapEntry is used in both SVs
      for(int k= 0; k < 2; k++) { // k goes over two xmapEntries which make up SV i
	for(int l= 0; l < 2; l++) { // l goes over two xmapEntries which make up SV j
	  if( sv_array[i].indices[k] == sv_array[j].indices[l] ) {
	    common = true;// "Common xmapEntry";
	    
	    /* first decide if common xmapEntry is the inverted MG or non-inverted MG for both SVs : If neither cannot apply this filter (HERE HERE : handle this case using dominant orientation) */
	    
	    bool inverted_i = (k==0 ? sv_array[i].inv_MGi : !sv_array[i].inv_MGi);
	    bool inverted_j = (l==0 ? sv_array[j].inv_MGi : !sv_array[j].inv_MGi);
	    if(INVERTED_DUP_FIX && inverted_i != inverted_j && exclude_type == inversion_type/*NEW377*/)
	      break;// this case not handled

	    //decide which is duplicate based on reference position: if both on same side, further one is duplicate
	    xmapEntry *entry1 = usexmapentries[sv_array[i].indices[k==0 ? 1 : 0]]; // k is common MG, need other one
	    xmapEntry *entry2 = usexmapentries[sv_array[j].indices[l==0 ? 1 : 0]]; // l is common MG, need other one
	    xmapEntry *entryc = usexmapentries[sv_array[j].indices[l]]; // common MG

	    // for ref
	    double entry1_rstart = entry1->refstartpos;
	    double entry1_rend   = entry1->refendpos;
	    double entry2_rstart = entry2->refstartpos;
	    double entry2_rend   = entry2->refendpos;
	    //	    double entryc_rstart = entryc->refstartpos;
	    //	    double entryc_rend   = entryc->refendpos; 

	    // for qry
	    double entry1_qstart = min(entry1->qrystartpos, entry1->qryendpos);
	    double entry1_qend   = max(entry1->qrystartpos, entry1->qryendpos);
	    double entry2_qstart = min(entry2->qrystartpos, entry2->qryendpos);
	    double entry2_qend   = max(entry2->qrystartpos, entry2->qryendpos);
	    double entryc_qstart = min(entryc->qrystartpos, entryc->qryendpos);
	    double entryc_qend   = max(entryc->qrystartpos, entryc->qryendpos);
	
	    if( entry1_qend < entryc_qend && entry2_qend < entryc_qend ) { // both on left -- allow overlap
	      if( entry1_qend > entry2_qend) { // entry1 ends to right of entry2, exclude SV which uses entry2
		if ( exclude_type != inversion_type/*NEW377*/ || inverted_i || entry1_rend < entry2_rend){
		  dupe = j;
		  sv_array[j].duplicate = true;
		}
	      } else if( entry1_qend < entry2_qend){
		if( exclude_type != inversion_type/*NEW377*/ || inverted_i || entry1_rend > entry2_rend){
		  dupe = i;
		  sv_array[i].duplicate = true;
		}
	      }
	    } //end if both on left
	    else if( entry1_qstart > entryc_qstart && entry2_qstart > entryc_qstart ) { //both on right
	      if( entry1_qstart > entry2_qstart ) { //entry1 starts to right of entry2, exclude SV which uses entry1
		if ( exclude_type != inversion_type/*NEW377*/|| !INVERTED_DUP_FIX || inverted_i || entry1_rstart < entry2_rstart){
		  dupe = i;
		  sv_array[i].duplicate = true;
		}
	      } else if( entry1_qstart < entry2_qstart ) {
		if ( exclude_type != inversion_type/*NEW377*/ || !INVERTED_DUP_FIX || inverted_i || entry1_rstart > entry2_rstart){
		  dupe = j;
		  sv_array[j].duplicate = true;
		}
	      }	      
	    } //end if both on right
	  } //end if indices match
	} //end for loop on xmapentries of inner SV
      } //end for loop on xmapentries of outter SV
      if( verbose > 1){ // debug 
	printf("%-19s : %2i %2i (%2i %2i, %2i %2i) ; ", common ? "Common xmapEntry" : "No common xmapEntry", i+1, j+1, 
	       usexmapentries[sv_array[i].indices[0]]->xmapid, usexmapentries[sv_array[i].indices[1]]->xmapid,
	       usexmapentries[sv_array[j].indices[0]]->xmapid, usexmapentries[sv_array[j].indices[1]]->xmapid);
	fflush(stdout);
      }
      if( verbose > 0 && dupe != -1 ) { //default is to print just dupes
	const char* typ = (exclude_type == inversion_type ? "inversion" : "translocation");
	int excl = (dupe == i ? j : i); //opposite mg from duplicated one
	printf("Duplicate %s: %2i (%2i %2i), excluded by %2i (%2i %2i)", typ, dupe+1, usexmapentries[sv_array[dupe].indices[0]]->xmapid, usexmapentries[sv_array[dupe].indices[1]]->xmapid, excl+1, usexmapentries[sv_array[excl].indices[0]]->xmapid, usexmapentries[sv_array[excl].indices[1]]->xmapid);
	fflush(stdout);
      }
      if( (verbose > 0 && dupe != -1) || verbose > 1 ){
	printf("\n");
	fflush(stdout);
      }
    } //end for inner loop on SVs (indiceslen)
  } //end for outter loop on SVs (indiceslen) for duplicate inversions
} //end checkDuplicates


//given two ref site idxs and one alignment, find correspond query location (by interpolating, if needed, based on nearest aligned query labels) and also return surrounding aligned query labels
void getQryCoord(Calign* p1, int refstartidx, int refstopidx, double& startlocQ, int& startidx, double& stoplocQ, int& stopidx) {

  Cmap *Ymap = YYmap[p1->mapid1];
  Cmap *Xmap = XXmap[p1->mapid2];
  FLOAT *Y = Ymap->site[0];
  FLOAT *X = Xmap->site[0];
  int N = Ymap->numsite[0];
  int M = Xmap->numsite[0];
  int U1 = p1->numpairs; //n aligned sites in first matchgroup

  if(DEBUG && !(0 < refstartidx && refstartidx < refstopidx && refstopidx <= N)){
    printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:refstartidx=%d,refstopidx=%d,N=%d\n",p1->mapid1,YYmap[p1->mapid1]->id,p1->mapid2,XXmap[p1->mapid2]->id,p1->orientation,refstartidx,refstopidx,N);
    fflush(stdout);
    assert(0 < refstartidx && refstartidx < refstopidx && refstopidx <= N);
  }

  int LstartidxQ=-1, RstartidxQ=-1, LstopidxQ=-1, RstopidxQ=-1; //find bounding index(s) in p1->sites2 which matches refstart{idx}
  int Lstartidx=-1,Rstartidx=-1,Lstopidx=-1,Rstopidx=-1;
  for(int t = 0; t < U1; t++) {
    if( p1->sites1[t] <= refstartidx ){
      Lstartidx = p1->sites1[t];
      LstartidxQ = p1->orientation ? M+1 - p1->sites2[t] : p1->sites2[t];
    }
    if( p1->sites1[t] >= refstartidx ){
      Rstartidx = p1->sites1[t];
      RstartidxQ = p1->orientation ? M+1 - p1->sites2[t] : p1->sites2[t];
      break;
    }
  }
  if(DEBUG && !(1 <= LstartidxQ && LstartidxQ <= M) && !(1 <= RstartidxQ && RstartidxQ <= M)){
    printf("getQryCoord: refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:refstartidx=%d,refstopidx=%d,N=%d: LstartidxQ=%d(I=%d),RstartidxQ=%d(I=%d),M=%d\n",
	   p1->mapid1,Ymap->id,p1->mapid2,Xmap->id,p1->orientation,refstartidx,refstopidx,N,LstartidxQ,Lstartidx,RstartidxQ,Rstartidx,M);
    fflush(stdout);
    assert((1 <= LstartidxQ && LstartidxQ <= M) || (1 <= RstartidxQ && RstartidxQ <= M));
  }

  for(int t = 0; t < U1; t++) {
    if( p1->sites1[t] <= refstopidx ){
      Lstopidx = p1->sites1[t];
      LstopidxQ = p1->orientation ? M+1 - p1->sites2[t] : p1->sites2[t];
    }
    if( p1->sites1[t] >= refstopidx ){    
      Rstopidx = p1->sites1[t];
      RstopidxQ = p1->orientation ? M+1 - p1->sites2[t] : p1->sites2[t];
      break;
    }
  }
  if(DEBUG && !(1 <= LstopidxQ && LstopidxQ <= M) && !(1 <= RstopidxQ && RstopidxQ <= M)){
    printf("getQryCoord: refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:refstartidx=%d,refstopidx=%d,N=%d: LstopidxQ=%d(I=%d),RstopidxQ=%d(I=%d),M=%d\n",
	   p1->mapid1,Ymap->id,p1->mapid2,Xmap->id,p1->orientation,refstartidx,refstopidx,N,LstopidxQ,Lstopidx,RstopidxQ,Rstopidx,M);
    fflush(stdout);
    assert((1 <= LstopidxQ && LstopidxQ <= M) || (1 <= RstopidxQ && RstopidxQ <= M));
  }

  if(1 <= LstartidxQ && LstartidxQ <= M && 1 <= RstartidxQ && RstartidxQ <= M){/* interpolate */
    if(Lstartidx < Rstartidx){/* interpolate */
      double ratio = (Y[refstartidx] - Y[Lstartidx]) / (Y[Rstartidx] - Y[Lstartidx]);
      startlocQ = X[LstartidxQ] + ratio * (X[RstartidxQ] - X[LstartidxQ]);
   } else {
      if(DEBUG) assert(Lstartidx == refstartidx && refstartidx == Rstartidx && LstartidxQ == RstartidxQ);
      startlocQ = X[LstartidxQ];
    }
    startidx   = LstartidxQ;// outermost (could be smaller or larger than stopidx)
  } else if(1 <= LstartidxQ && LstartidxQ <= M){/* extrapolate from left */
    double shift = Y[refstartidx] - Y[Lstartidx];
    startlocQ = X[LstartidxQ] + shift;
    startidx   = LstartidxQ;// outermost (could smaller or larger than stopidx)
  } else {/* extrapolate from right */
    if(DEBUG) assert(1 <= RstartidxQ && RstartidxQ <= M);
    double shift = Y[Rstartidx] - Y[refstartidx];
    startlocQ = X[RstartidxQ] - shift;
    startidx = RstartidxQ;
  }

  if(1 <= LstopidxQ && LstopidxQ <= M && 1 <= RstopidxQ && RstopidxQ <= M){/* interpolate */
    if(Lstopidx < Rstopidx){/* interpolate */
      double ratio = (Y[refstopidx] - Y[Lstopidx]) / (Y[Rstopidx] - Y[Lstopidx]);
      stoplocQ = X[LstopidxQ] + ratio * (X[RstopidxQ] - X[LstopidxQ]);
   } else {
      if(DEBUG) assert(Lstopidx == refstopidx && refstopidx == Rstopidx && LstopidxQ == RstopidxQ);
      stoplocQ = X[LstopidxQ];
    }
    stopidx    = RstopidxQ;// outermost
  } else if(1 <= LstopidxQ && LstopidxQ <= M){/* extrapolate from left */
    double shift = Y[refstopidx] - Y[Lstopidx];
    stoplocQ = X[LstopidxQ] + shift;
    stopidx   = LstopidxQ;
  } else {/* extrapolate from right */
    if(DEBUG) assert(1 <= RstopidxQ && RstopidxQ <= M);
    double shift = Y[Rstopidx] - Y[refstopidx];
    stoplocQ = X[RstopidxQ] - shift;
    stopidx = RstopidxQ;// outermost
  }
  
  /* convert from kb to bp */
  startlocQ *= 1e3;
  stoplocQ *= 1e3;

} //end getQryCoord

/* return label I in 0..N+1 so that Y[I] is closest to position (NOTE : position can be outside the range Y[0] .. Y[N+1]) */
static int findlabel(double position, FLOAT *Y, int N)
{
  // use binary search
  int low = 0, high = N+1;

  while(high > low + 1){
    int mid = (high + low) / 2;
    if(position < Y[mid])
      high = mid;
    else
      low = mid;
  }

  return fabs(position - Y[low]) < fabs(position - Y[high]) ? low : high;
}

/* return number of overlaped query labels between two matchgroups */
static int QryLabelOverlap(Calign *p, Calign *q, int &U1, int &U2)
{
  int cnt = 0;
  U1 = p->numpairs;
  U2 = q->numpairs;

  int t1 = 0, t2 = 0;// number of matching aligned query labels

  if(p->orientation == q->orientation){
    for(; t1 < U1 && t2 < U2; ){
      if(p->sites2[t1] < q->sites2[t2]){
	t1++;
	continue;
      }
      if(p->sites2[t1] > q->sites2[t2]){
	t2++;
	continue;
      }
      if(DEBUG) assert(p->sites2[t1] == q->sites2[t2]);
      cnt++;
      t1++;
      t2++;
    }
  } else {// opposite orientation matchgroups
    int M = XXmap[p->mapid2]->numsite[0];
    if(p->orientation){
      for(; t1 < U1 && t2 < U2; ){
	if(M+1 - p->sites2[U1-1-t1] < q->sites2[t2]){
	  t1++;
	  continue;
	}
	if(M+1 - p->sites2[U1-1-t1] > q->sites2[t2]){
	  t2++;
	  continue;
	}
	if(DEBUG) assert(M+1 - p->sites2[U1-1-t1] == q->sites2[t2]);
	cnt++;
	t1++;
	t2++;
      }
    } else {
      for(; t1 < U1 && t2 < U2; ){
	if(VERB>=3){
	  printf("t1=%d,t2=%d:U1=%d,U2=%d,M=%d:p->sites2[t1]=%d,q->sites2[U2-1-t2]=%d: cnt=%d\n",
		 t2,t2,U1,U2,M,p->sites2[t1],q->sites2[U2-1-t2],cnt);
	  fflush(stdout);
	}
	if(p->sites2[t1] < M+1 - q->sites2[U2-1-t2]){
	  t1++;
	  continue;
	}
	if(p->sites2[t1] > M+1 - q->sites2[U2-1-t2]){
	  t2++;
	  continue;
	}
	if(DEBUG) assert(p->sites2[t1] == M+1- q->sites2[U2-1-t2]);
	cnt++;
	t1++;
	t2++;
      }
    }
  }

  return cnt;
}

extern RFLOAT OutlierEndBias, OutlierEndPenalty;// see RGentigRefScore.h
extern double Sm(int J, int I, int K, double *Y);// see RGentigRefScore.h
extern double alignFP(Calign *p, FLOAT *Y, FLOAT *X, int N, int M, int orientation, int Yid, int Xid, double lscore, int verb);// see refalign.cpp

static int MaxXmapid;
static int *firstIndel = NULL;

#define CVERBOSE 0

/* check for outliers of smap_Outlier_minQrySites sites in qry of larger matchgroup (usexmapentries[j]) that overlaps smaller matchgroup in query.

A : fulloverlap == true :
   1. If no such outliers, mark usexmapentries[i] for deletion (unless MATCHGROUP_TRIM <=2 && usexmapentries[i]->refcontigid != usexmapentries[j]->refcontigid) and return.
   2. For largest such outlier, trim usexmapentries[i] to only include the outlier region on query. 
   3. IF saveoutliers == true : For each additional outlier, add a new xmapEntry to end of xmapentries[0..numxmapentries-1] that matches that part of the original usexmapentries[i]. 
   4. Set xmapid_overlapped == usexmapentries[j]->xmapid for all MGs in steps 2 & 3.
B : fulloverlap == false : (HERE HERE HERE : should determine which map has higher confidence in overlap region and trim the other one instead of always trimming usexmapentries[i])
   1. IF saveoutliers == true : For each outlier add a new xmapEntry to end of xmapentries[0..numxmapentries-1] that matches that part of the original usexmapentries[i]
      Also set xmapid_overlapped == usexmapentries[j]->xmapid for each such new MG.
   2. Trim away the entire overlapped region of usexmapentries[i]
      If usexmapentries[i] is larger and fully overlaps usexmapentries[j] prefer the non-overlapped region towards right end of query (with query oriented as in usexmapentries[j])
      unless keepleft = true is specified

return 1 IFF any changes were made to usexmapentries[i] or new MGs appended to xmapentries[]
*/

// NOTE : usexmapentries[j] and usexmapentries[i] have the same querycontigid but NOT necessarily the same refcontigid

static __attribute__ ((noinline)) int checkOverlap(xmapEntry **&usexmapentries,
			 int &smapsize,
                         xmapEntry **&xmapentries,
			 int &numxmapentries,
			 int &maxxmapentries,
			 int j, // usexmapentries[j] is the typically larger (higher confidence) matchgroup 
			 int i, // usexmapentries[i] is the typically smaller (lower confidence) matchgroup of (possibly) opposite orientation : this matchgroup may be truncated/filtered/broken
			 bool fulloverlap, // true IFF usexmapentries[i] is fully overlapped by usexmapentries[j] on query
			 bool saveoutliers, // If false (AND fulloverlap == false) : Don't save outliers as extra xmapEntry (step B1 above)
   		         int numsvIndel = -1, // If >= 0 : i is the larger MG and if any portion of it is filtered, any Indel sv_array[0..numsvIndel-1] from that region needs to be filtered out
			 bool keepleft = false // If true AND !fulloverlap: retain the nonoverlapped region of usexmapentries[i] towards left end of query (oriented as in usexmapentries[j])		       
			 )
{
  int ret = 0;

  Calign *q = usexmapentries[i]->align;
  int matchcnt = q->numpairs;

  Calign *p = usexmapentries[j]->align;
  int U = p->numpairs;

  Cmap *Xmap = usexmapentries[j]->QryMap;
  Cmap *Yjmap = usexmapentries[j]->RefMap;
  Cmap *Yimap = usexmapentries[i]->RefMap;
  FLOAT *Yj = Yjmap->site[0];
  FLOAT *Yi = Yimap->site[0];
  int Nj = usexmapentries[j]->refcontignsite;
  int Ni = usexmapentries[i]->refcontignsite;
  FLOAT *X = Xmap->site[0];
  int M = usexmapentries[i]->querycontignsite;

  int orientation = usexmapentries[j]->orientforward ? 0 : 1;

  int refstartidx = usexmapentries[i]->refstartidx;
  int refendidx = usexmapentries[i]->refendidx;

  bool inverted = (usexmapentries[i]->orientforward != usexmapentries[j]->orientforward);
  int qrystartidx = inverted ? usexmapentries[i]->qryendidx : usexmapentries[i]->qrystartidx; // absolute index of "leftmost" MGi qry label with query oriented same as in usexmapentries[j]
  int qryendidx = inverted ? usexmapentries[i]->qrystartidx : usexmapentries[i]->qryendidx; // absolute index of "rightmost" MGi qry label with query oriented same as in usexmapentries[j]

  int qrystartidx1 = orientation ? (M+1 - qrystartidx) : qrystartidx;// relative index of leftmost MGi qry label with qry oriented as in usexmapentries[j] (larger matchgroup)
  int qryendidx1 = orientation ? (M+1 - qryendidx) : qryendidx;// relative index rightmost MGi qry label with qry oriented as in usexmapentries[j] (larger matchgroup)

  if(CVERBOSE){
    printf("Checking %sqry overlap : refid=%lld,%lld,qryid=%lld,j=%d,xmapid1=%d(conf=%0.2f,np=%d),i=%d,xmapid2=%d(conf=%0.2f,np=%d),fwd=%d,%d:RefSmall= %d..%d(%0.3f..%0.3f),QrySmall= %d..%d(%0.3f..%0.3f), M=%d,Ni=%d,Nj=%d\n",
	   fulloverlap ? "full " : "", usexmapentries[j]->refcontigid,usexmapentries[i]->refcontigid,usexmapentries[j]->qrycontigid,j,usexmapentries[j]->xmapid,usexmapentries[j]->confidence,U,
	   i,usexmapentries[i]->xmapid,usexmapentries[i]->confidence,matchcnt, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[i]->orientforward ? 1 : 0, 
	   refstartidx,refendidx, Yi[refstartidx],Yi[refendidx], qrystartidx,qryendidx,X[qrystartidx],X[qryendidx],M,Ni,Nj);
    printf("\t qry1=%d..%d,M=%d,inverted=%d,orientation=%d:xmapid1:Ref= %d..%d(%0.3f..%0.3f),Qry= %d..%d(%0.3f..%0.3f),xmapid_overlapped=%d,%d,saveoutliers=%d,fulloverlap=%d,numsvIndel=%d\n",
	   qrystartidx1,qryendidx1,M,inverted,orientation,
	   usexmapentries[j]->refstartidx,usexmapentries[j]->refendidx, Yj[usexmapentries[j]->refstartidx],Yj[usexmapentries[j]->refendidx], 
	   usexmapentries[j]->qrystartidx,usexmapentries[j]->qryendidx, X[usexmapentries[j]->qrystartidx],X[usexmapentries[j]->qryendidx],
	   usexmapentries[j]->xmapid_overlapped,usexmapentries[i]->xmapid_overlapped,saveoutliers,fulloverlap,numsvIndel);
    fflush(stdout);
  }

  //  if(DEBUG>=1+RELEASE && !numsvIndel) assert(X[qryendidx]-X[qryendidx] <= fabs(X[usexmapentries[j]->qryendidx]-X[usexmapentries[j]->qrystartidx]));
  if(DEBUG>=1+RELEASE && numsvIndel < 0 && !(usexmapentries[j]->confidence >= usexmapentries[i]->confidence)){
    printf("numsvIndel= %d\n",numsvIndel);
    fflush(stdout);
    assert(usexmapentries[j]->confidence >= usexmapentries[i]->confidence);
  }

  if(DEBUG>=1+RELEASE) assert(usexmapentries[i]->numrefsites >= usexmapentries[i]->refendidx - usexmapentries[i]->refstartidx + 1);
  if(DEBUG>=1+RELEASE) assert(usexmapentries[j]->numrefsites >= usexmapentries[j]->refendidx - usexmapentries[j]->refstartidx + 1);

  if(!fulloverlap && !saveoutliers){/* check if this is a full overlap of MGi by MGj */
    int JL = p->sites2[0], JR = p->sites2[p->numpairs-1];
    if(min(JR,qryendidx1) - max(JL,qrystartidx1) >= qryendidx1 - qrystartidx1 /* WAS42 min(JR-JL, qryendidx1-qrystartidx1)*/){
#if 1 // NEW 
      // Switch to fulloverlap == true : need to be able to handle deletion of smaller matchgroup j with sv_overlap == true in main SV loop
      fulloverlap = true;
      if(CVERBOSE){
	printf("\t Full qry overlap : continuing with fulloverlap=true: JL=%d, JR=%d, qrystartdx1=%d, qryendidx1=%d\n",JL,JR,qrystartidx1,qryendidx1);
	fflush(stdout);
      }
#else // OLD
      /* there is nothing to do */
      if(CVERBOSE){
	printf("\t Full qry overlap : Nothing trimmed: JL=%d, JR=%d, qrystartdx1=%d, qryendidx1=%d\n",JL,JR,qrystartidx1,qryendidx1);
	fflush(stdout);
      }
      if(DEBUG>=1+RELEASE) assert(usexmapentries[i]->xmapid_overlapped == usexmapentries[j]->xmapid);
      return 0;
#endif
    }
  }

  if(DEBUG) assert(usexmapentries[j]->qrycontigid == usexmapentries[i]->qrycontigid);
  if(DEBUG) assert(usexmapentries[j]->querycontignsite == M);
  if(DEBUG) assert(1 <= refstartidx && refstartidx <= refendidx && refendidx <= Ni);

  if(DEBUG) assert(refstartidx == q->sites1[0] - q->sitesK1[0]);
  if(DEBUG && !(refendidx == q->sites1[matchcnt-1])){
    printf("\t refendidx= %d, matchcnt=%d, q->sites1[matchcnt-1]= %d, q->sitesK1[matchcnt-1]= %d\n",
	   refendidx, matchcnt, q->sites1[matchcnt-1], q->sitesK1[matchcnt-1]);
    for(int t = 0; t < q->numpairs; t++)
      printf("\t\t t=%d: I=%d,K=%d,J=%d (q->sites2[t]=%d)\n",
	     t, q->sites1[t],q->sitesK1[t], q->orientation ? M+1 - q->sites2[t] : q->sites2[t], q->sites2[t]);
    fflush(stdout);
    assert(refendidx == q->sites1[matchcnt-1]);
  }
  if(DEBUG) assert(usexmapentries[i]->qrystartidx == (q->orientation ? M+1 - q->sites2[0] : q->sites2[0]));
  if(DEBUG) assert(usexmapentries[i]->qryendidx == (q->orientation ? M+1 - q->sites2[matchcnt-1] : q->sites2[matchcnt-1]));

  if(DEBUG) assert(refstartidx <= refendidx);
  if(DEBUG && !(qrystartidx1 <= qryendidx1)){
    #pragma omp critical
    {
      printf("\t refid=%lld,qryid=%lld:qry1=%d..%d\n",usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, qrystartidx1,qryendidx1);
      fflush(stdout);
      assert(qrystartidx1 <= qryendidx1);
    }
  }

  int localtype = (PoutlierEnd > 0.0) ? -3 : -2;

  /* locate all outliers in qry overlap region of larger matchgroup j that are smap_Outlier_minQrySites labels or larger */
  int outliers = 0;
  int qL[16],qR[16];
  int tL[16],tR[16];
  int maxOverlap = 0, maxoutlier = -1;

  if(fulloverlap || saveoutliers){
    int lastJ = p->sites2[0];
    for(int t = 1; t < U; t++){
      int J = p->sites2[t];
      int Overlap = min(J, qryendidx1) - max(lastJ, qrystartidx1);
      if(Overlap >= smap_Outlier_minQrySites - 1){
	if(outliers >= 16){
	  printf("WARNING: more than 16 outliers in qry overlap region of refid=%lld,qryid=%lld (qry1=%d..%d)\n",usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,qrystartidx1,qryendidx1);
	  fflush(stdout);
	  break;
	}

	/* WAS83	if(Overlap > maxOverlap){
	  maxOverlap = Overlap;
	  maxoutlier = outliers;
	  } */

	qL[outliers] = max(lastJ,qrystartidx1);
	qR[outliers] = min(J,qryendidx1);

#if 0 // OLD CODE
	/* determinate largest range q->sites1[tL2 .. tR2] that falls within qL[outlier] .. qR[outlier] (NOTE : also include any interval that partly overlaps qL[outlier] .. qR[outlier]) */
	int tL2 = 0, tR2 = q->numpairs - 1;
	int firstJ2 = inverted ? M+1 - q->sites2[tL2] : q->sites2[tL2];
	for(tL2 = 1; tL2 <= tR2; tL2++){
	  int J2 = inverted ? M+1 - q->sites2[tL2] : q->sites2[tL2];
	  int JL = min(firstJ2,J2), JR = max(firstJ2,J2);
	  if(max(qL[outliers], JL) <  min(qR[outliers],JR))
	    break;
	  firstJ2 = J2;
	}
	//	if(DEBUG) assert(tL2 <= tR2);
	tL2--;

	int lastJ2 = inverted ? M+1 - q->sites2[tR2] : q->sites2[tR2];
	for(tR2--; tL2 <= tR2; tR2--){
	  int J2 = inverted ? M+1 - q->sites2[tR2] : q->sites2[tR2];
	  int JL = min(lastJ2,J2), JR = max(lastJ2,J2);
	  if(max(qL[outliers], JL) <  min(qR[outliers],JR))
	    break;
	  lastJ2 = J2;
	}
	//	if(DEBUG) assert(tL2 <= tR2);
	tR2++;

	if(CVERBOSE){
	  if(tL2 < tR2)
	    printf("\t outlier %d: qry1=%d..%d, outlierJ= %d..%d, t=%d..%d/%d(J=%d..%d)\n",outliers,qrystartidx1,qryendidx1,lastJ,J, tL2,tR2,q->numpairs,
		   inverted ? M+1 - q->sites2[tL2] : q->sites2[tL2], inverted ? M+1 - q->sites2[tR2] : q->sites2[tR2]);
	  else {
	    printf("\t outlier %d: qry1=%d..%d, outlierJ= %d..%d, t=%d..%d/%d\n",outliers,qrystartidx1,qryendidx1,lastJ,J, tL2,tR2,q->numpairs);
	    for(int t = 0; t < q->numpairs; t++)
	      printf("\t\t t=%d:q->sites2[t]=%d, J=%d\n",t, q->sites2[t], inverted ? M+1 - q->sites2[t] : q->sites2[t]);
	  }
	  fflush(stdout);
	}
	if(DEBUG) assert(tL2 < tR2);

#else // NEW CODE (NEW82)
	/* determine largest range q->sites2[tL2 .. tR2] that falls within qL[outlier] .. qR[outlier] (NOTE : also include any interval that is at least 50% overlaped with qL[outlier] .. qR[outlier]) */
	if(DEBUG) assert(qL[outliers] < qR[outliers]);
	int absJL, absJR;
	FLOAT L1 = X[absJL = (usexmapentries[j]->orientforward ? qL[outliers] : M+1 - qR[outliers])];
	FLOAT R1 = X[absJR = (usexmapentries[j]->orientforward ? qR[outliers] : M+1 - qL[outliers])];
	if(DEBUG) assert(L1 < R1);
	if(CVERBOSE){
	  printf("\t found outlier %d: qL=%d,qR=%d(abs J1=%d..%d), L1= %0.4f, R1= %0.4f, Len1= %0.4f kb\n",outliers,qL[outliers],qR[outliers],absJL,absJR,L1,R1,X[M+1]);
	  fflush(stdout);
	}

	int tL2 = 0, tR2 = q->numpairs - 1;
	int firstJ2 = usexmapentries[i]->orientforward ? q->sites2[tL2] : M+1 - q->sites2[tL2];
	for(tL2 = 1; tL2 <= tR2; tL2++){
	  int J2 = usexmapentries[i]->orientforward ? q->sites2[tL2] : M+1 - q->sites2[tL2];
	  FLOAT L2 = min(X[J2],X[firstJ2]);
	  FLOAT R2 = max(X[J2],X[firstJ2]);
	  if(CVERBOSE && VERB>=2){
	    printf("\t tL2=%d/%d:firstJ2=%d,J2=%d:L1=%0.4f,R1=%0.4f,L2=%0.4f,R2=%0.4f:overlap= %0.4f, R2-L2= %0.4f kb\n",
		   tL2,tR2,firstJ2,J2,L1,R1,L2,R2,min(R1,R2) - max(L1,L2), R2-L2);
	    fflush(stdout);
	  }
	  if(DEBUG>=1+RELEASE) assert(L2 < R2);
	  if(min(R1,R2) > max(L1,L2)){
	    if(min(R1,R2) - max(L1,L2) <= 0.5 * (R2-L2))
	      tL2++;
	    break;
	  }
	  firstJ2 = J2;
	}
	tL2--;

	int lastJ2 = usexmapentries[i]->orientforward ? q->sites2[tR2] : M+1 - q->sites2[tR2];
	for(tR2--; tL2 <= tR2; tR2--){
	  int J2 = usexmapentries[i]->orientforward ? q->sites2[tR2] : M+1 - q->sites2[tR2];
	  double L2 = min(X[J2],X[lastJ2]);
	  double R2 = max(X[J2],X[lastJ2]);
	  if(CVERBOSE && VERB>=2){
	    printf("\t tR2=%d(tL2=%d):lastJ2=%d,J2=%d:L1=%0.4f,R1=%0.4f,L2=%0.4f,R2=%0.4f:overlap= %0.4f, R2-L2= %0.4f kb\n",
		   tR2,tL2,lastJ2,J2,L1,R1,L2,R2,min(R1,R2) - max(L1,L2), R2-L2);
	    fflush(stdout);
	  }
	  if(DEBUG>=1+RELEASE) assert(L2 < R2);
	  if(min(R1,R2) > max(L1,L2)){
	    if(min(R1,R2) - max(L1,L2) <= 0.5 * (R2-L2))
	      tR2--;
	    break;
	  }
	  lastJ2 = J2;
	}
	tR2++;

	if(CVERBOSE){
	  if(tL2 < tR2)
	    printf("\t outlier %d: qry1=%d..%d, qL=%d,qR=%d(abs J1=%d..%d): t=%d..%d/%d(abs J2=%d..%d),maxoutlier=%d\n",outliers,qrystartidx1,qryendidx1,qL[outliers],qR[outliers],absJL,absJR,tL2,tR2,q->numpairs,
		   usexmapentries[i]->orientforward ? q->sites2[tL2] : M+1 - q->sites2[tL2], usexmapentries[i]->orientforward ? q->sites2[tR2] : M+1 - q->sites2[tR2],maxoutlier);
	  else {
	    printf("\t Failed outlier %d: qry1=%d..%d, qL=%d,qR=%d(absJ1=%d..%d): t=%d..%d/%d\n",outliers,qrystartidx1,qryendidx1,qL[outliers],qR[outliers],absJL,absJR, tL2,tR2,q->numpairs);
	    if(CVERBOSE>=2){
	      for(int t = 0; t < q->numpairs; t++){
		int J2 = usexmapentries[i]->orientforward ? q->sites2[t] : M+1 - q->sites2[t];
		printf("\t\t t=%d:q->sites2[t]=%d, abs J2=%d (X[J2]= %0.4f kb)\n",t, q->sites2[t], J2, X[J2]);
	      }
	    }
	  }
	  fflush(stdout);
	}
	if(DEBUG) assert(0 <= tL2 && tL2 <= tR2 && tR2 < q->numpairs);
#endif

	if(tL2 < tR2 && abs(q->sites2[tR2] - q->sites2[tL2] + 1) >= smap_Outlier_minQrySites - 1){
	  if(DEBUG) assert(0 <= tL2 && tL2 < tR2 && tR2 < q->numpairs);
	  tL[outliers] = tL2;
	  tR[outliers] = tR2;

	  if(tR2 - tL2 + 1 > maxOverlap){// NEW83
	    maxOverlap = tR2 - tL2;
	    maxoutlier = outliers;
	  }

	  outliers++;
	}
      }
      lastJ = J;
    }

    if(fulloverlap && !outliers){
      if(MATCHGROUP_TRIM >= 3 || usexmapentries[i]->refcontigid == usexmapentries[j]->refcontigid){
	if(CVERBOSE){
	  printf("\t No outliers of %d labels or more found in qry overlap region of matchgroup j: filtering out matchgroup i\n",smap_Outlier_minQrySites);
	  fflush(stdout);
	}

	usexmapentries[i]->sv_overlap = true;// delete smaller matchgroup
      } else {
	if(CVERBOSE){
	  printf("\t No outliers of %d labels or more found in qry overlap region of matchgroup j\n",smap_Outlier_minQrySites);
	  fflush(stdout);
	}
      }
      return 1;
    }
  }

  /* recompute local X[] as Xrev[], which may depend on ScaleID and p->orientation */
  double *Xrev = new double[M+2];
  FLOAT scale = ScaleFactor[q->scaleID];
  if(!q->orientation){
    for(int J=0;J <= M+1;J++)
      Xrev[J] = X[J]*scale;
  } else {
    for(int j=0;j <= M+1;j++)
      Xrev[j] = (X[M+1]-X[M+1-j])*scale;
  }

  double ChimScore = OutlierEndBias + OutlierEndPenalty;

  if(outliers && (outliers >= 2 || !fulloverlap) && saveoutliers){/* for each outlier other than the largest one, create a new xmapentries[numxmapentries++] */
    /*    xmapEntry **newxmapentries = new xmapEntry*[numxmapentries + (outliers - 1)];
    for(int i = 0; i < numxmapentries; i++)
      newxmapentries[i] = xmapentries[i];
    
      delete [] xmapentries;  xmapentries = newxmapentries;*/

    if(DEBUG) assert(maxoutlier >= 0 && maxoutlier < outliers);
    for(int i = 0; i < outliers; i++){
      if(fulloverlap && i == maxoutlier)
	continue;

      /* create alignment subset q->sites*[tL[i] .. tR[i]] */
      ret = 1;
      maxalignallocNULL(numaligns+1, alignment, numaligns, maxaligns);
      Calign *r = alignment[numaligns++] = new Calign[1];
      copy(r,q,1,1);
      int S = tL[i], E = tR[i];
      int nU = E - S + 1;

      if(S > 0){/* new left endoutlier */
	int I = q->sites1[S];
	int K = q->sitesK1[S];
	int J = q->sites2[S];
	r->iscore[0] = r->outscore[0] = ChimScore + Sm(0,I,K,Yi);
	r->Lij1 = I-K;
	r->Lij2 = J;
	r->Lend = localtype;
      }

      int t = S;
      int lastI = r->sites1[0] = q->sites1[t];
      int lastK = r->sitesK1[0] = q->sitesK1[t];
      int lastJ = r->sites2[0] = q->sites2[t];
      r->noutliers = 0;
      r->maxoutlier = 0.0;
      r->maxoutlierLabels = 0;
      for(int u = 1; ++t <= E; u++){
	r->iscore[u] = q->iscore[t];
	r->outscore[u] = q->outscore[t];
	int I = r->sites1[u] = q->sites1[t];
	int K = r->sitesK1[u] = q->sitesK1[t];
	int J = r->sites2[u] = q->sites2[t];
	if(r->outscore[u] + (FLOAT)0.01 < r->iscore[u]){
	  r->noutliers++;
	  double deltaY = Yc(Yi,I,K) - Yc(Yi,lastI,lastK);
	  double deltaX = Xrev[J] - Xrev[lastJ];
	  if(DEBUG) assert(deltaX > 0.0);
	  double delta = fabs(deltaY-deltaX);
	  r->maxoutlier = max(delta, r->maxoutlier);
	  r->maxoutlierLabels = max(I-K-lastI + J-lastJ - 2, r->maxoutlierLabels);
	}
	lastI = I;
	lastK = K;
	lastJ = J;
      }
      
      if(E < q->numpairs - 1){/* new right end */
	int I = q->sites1[E];
	//	int K = q->sitesK1[E];
	int J = q->sites2[E];
	r->iscore[nU] = r->outscore[nU] = ChimScore;
	r->Rij1 = I;
	r->Rij2 = J;
	r->Rend = localtype;
      } else if(S > 0)
	r->iscore[nU] = r->outscore[nU] = r->iscore[r->numpairs];

      r->numpairs = nU;

      r->score = 0;
      for(int t = 0; t <= nU; t++)
	r->score += r->iscore[t];
      
      r->logPV = alignFP(r, Yi, Xrev, Ni, M, r->orientation, r->mapid1, r->mapid2, r->score, 0);
      r->chimpair = 1;// flag as split matchgroup so it isn't filtered out by LogPvThreshold
      r->logPV2 = max(r->logPV, q->logPV2);
      if(DEBUG>=2) assert(r->logPV2 >= r->logPV);

      int qrystartidx0 = r->orientation ? M + 1 - r->sites2[0] : r->sites2[0];
      int qryendidx0 = r->orientation ? M + 1 - r->sites2[r->numpairs - 1] : r->sites2[r->numpairs - 1];
      double qrystartpos0 = X[qrystartidx0] * 1000.0;
      double qryendpos0 = X[qryendidx0] * 1000.0;
      int refstartidx0 = r->sites1[0] - r->sitesK1[0];
      int refendidx0 = r->sites1[r->numpairs - 1];
      double refstartpos0 = Yi[refstartidx0] * 1000.0;
      double refendpos0 = Yi[refendidx0] * 1000.0;
      int startua = r->orientation ? qryendidx0 - 1 : qrystartidx0 - 1;
      int endua   = r->orientation ? M - qrystartidx0 : M - qryendidx0;
      if(DEBUG) assert(startua == min(qrystartidx0,qryendidx0) - 1);
      if(DEBUG) assert(endua == M - max(qrystartidx0,qryendidx0));

      int k = numxmapentries;
      if(k >= maxxmapentries){
	if(VERB>=3){
	  for(int i = 0; i < numxmapentries; i++){
	    xmapEntry *p = usexmapentries[i];
	    if(p->xmapid == 868 || p->xmapid == 866){
	      printf("\t usexmapentries[%d]: refid=%lld,qryid=%lld, xmapid=%d, xmapid_overlapped=%d\n",
		     i, p->refcontigid,p->qrycontigid,p->xmapid,p->xmapid_overlapped);
	      fflush(stdout);
	    }
	  }
	}
	long long nmaxxmapentries = max((k * 5LL) >> 2LL, 64LL);
	if(nmaxxmapentries > MASK_LL(31) || nmaxxmapentries < 0){
	  printf("31-bit overflow of numxmapentries= %d -> %lld\n",numxmapentries,nmaxxmapentries);
	  fflush(stdout);exit(1);
	}
	xmapEntry **newxmapentries = new xmapEntry*[nmaxxmapentries];
	xmapEntry **newusexmapentries = new xmapEntry*[nmaxxmapentries];
	if(CVERBOSE){
	  printf("re-allocating xmapentries[],usexmapentries[] : maxxmapentries = %d -> %lld, usexmapentries -> %p, xmapentries -> %p\n", numxmapentries,nmaxxmapentries, newusexmapentries, newxmapentries);
	  fflush(stdout); 
	}
	maxxmapentries = nmaxxmapentries;
	memcpy(newxmapentries, xmapentries, k * sizeof(xmapEntry *));
	memcpy(newusexmapentries, usexmapentries, smapsize * sizeof(xmapEntry *));
	delete [] xmapentries; xmapentries = newxmapentries;
	delete [] usexmapentries; usexmapentries = newusexmapentries;

	if(VERB>=3){
	  for(int i = 0; i < numxmapentries; i++){
	    xmapEntry *p = usexmapentries[i];
	    if(p->xmapid == 868 || p->xmapid == 866){
	      printf("\t usexmapentries[%d]: refid=%lld,qryid=%lld, xmapid=%d, xmapid_overlapped=%d\n",
		     i, p->refcontigid,p->qrycontigid,p->xmapid,p->xmapid_overlapped);
	      fflush(stdout);
	    }
	  }
	}
      }

      usexmapentries[smapsize] = xmapentries[k] = new xmapEntry(r, k+1, Xmap->id, Yimap->id, qrystartpos0, qryendpos0, qrystartidx0, qryendidx0, refstartpos0, refendpos0,refstartidx0,refendidx0,!bool(r->orientation), r->logPV /* WAS176 r->logPV2 */, startua, endua, X[M+1], refendidx0 - refstartidx0 + 1, M, Ni) ;
      xmapentries[k]->RefMap = Yimap;
      xmapentries[k]->QryMap = Xmap;
      for(int t = refstartidx0; t <= refendidx0; t++){
	xmapentries[k]->refsites[t-refstartidx0] = Yi[t] * 1000.0;
	if(VERB>=2){
	  printf("usexmapentries[%d]->refsites[%d] = %0.1f (xmapid=%d)\n", smapsize, t-refstartidx, xmapentries[k]->refsites[t-refstartidx], xmapentries[k]->xmapid);
	  fflush(stdout);
	}
	if(DEBUG>=2) assert(xmapentries[k]->refstartpos <= Yi[t] * 1000.0  && Yi[t] * 1000.0 <= xmapentries[k]->qryendpos);
      }
      if(r->orientation){/* reverse orientation */
	if(DEBUG) assert(xmapentries[k]->qrystartidx >= xmapentries[k]->qryendidx);
	if(DEBUG) assert(xmapentries[k]->qrystartpos >= xmapentries[k]->qryendpos);
	for(int t = xmapentries[k]->qryendidx; t <= xmapentries[k]->qrystartidx; t++){
	  xmapentries[k]->qrysites[t - xmapentries[k]->qryendidx] = X[t] * 1000.0;
	  if(DEBUG>=2) assert(xmapentries[k]->qryendpos <= X[t] * 1000.0  && X[t] * 1000.0 <= xmapentries[k]->qrystartpos);
	}
      } else {
	if(DEBUG) assert(xmapentries[k]->qrystartidx <= xmapentries[k]->qryendidx);
	if(DEBUG) assert(xmapentries[k]->qrystartpos <= xmapentries[k]->qryendpos);
	for(int t = xmapentries[k]->qrystartidx; t <= xmapentries[k]->qryendidx; t++){
	  xmapentries[k]->qrysites[t - xmapentries[k]->qrystartidx] = X[t] * 1000.0;
	  if(DEBUG>=2) assert(xmapentries[k]->qrystartpos <= X[t] * 1000.0  && X[t] * 1000.0 <= xmapentries[k]->qryendpos);
	}
      }

      if(!fulloverlap){
	xmapentries[k]->sv_overlap = true;
	xmapentries[k]->output = true;
      }
      xmapentries[k]->xmapid_overlapped = usexmapentries[j]->xmapid;

      if(CVERBOSE)
	printf("Splitting qry overlap region:smapsize=%d: refid=%lld,qryid=%lld,xmapid=%d(conf=%0.2f,np=%d),fwd=%d:Ref= %d..%d(%0.3f..%0.3f),Qry= %d..%d(%0.3f..%0.3f),M=%d,Ni=%d(xmapid_overlapped -> %d)\n",
	       smapsize, xmapentries[k]->refcontigid,xmapentries[k]->qrycontigid,k+1,xmapentries[k]->confidence,r->numpairs,!r->orientation,refstartidx0,refendidx0,Yi[refstartidx0],Yi[refendidx0],
	       qrystartidx0,qryendidx0,X[qrystartidx0],X[qryendidx0], M, Ni, xmapentries[k]->xmapid_overlapped);

      if(VERB>=3){
	for(int i = 0; i < numxmapentries; i++){
	  xmapEntry *p = usexmapentries[i];
	  if(p->xmapid == 868 || p->xmapid == 866){
	    printf("\t usexmapentries[%d]: refid=%lld,qryid=%lld, xmapid=%d, xmapid_overlapped=%d\n",
		   i, p->refcontigid,p->qrycontigid,p->xmapid,p->xmapid_overlapped);
	    fflush(stdout);
	  }
	}
      }

      numxmapentries++;
      smapsize++;
    }// for(int i = 0; i < outliers; i++)
  }// if(outliers >= 2)

  int t2, r2;
  if(!fulloverlap){
    /* locate region of usexmapentries[i] that is NOT overlapped on query by larger matchgroup usexmapentries[j] */
    int JL = p->sites2[0], JR = p->sites2[p->numpairs-1];

    t2 = 0;
    r2 = q->numpairs - 1;

    if((DEBUG && !(qryendidx1 > JR || qrystartidx1 < JL)) || (VERB>=3)){
      printf("JL=%d,JR=%d,qrystartidx1=%d,qryendidx1=%d,q->sites2[0,%d]= %d,%d\n",JL,JR,qrystartidx1,qryendidx1,r2,q->sites2[0],q->sites2[r2]);
      fflush(stdout);
      assert(qryendidx1 > JR || qrystartidx1 < JL);
    }

    if(qryendidx1 > JR && !(keepleft && qrystartidx1 < JL)){/* include label JR (or first aligned label <= JR) in NOT overlapped range */
      if(!inverted){
	for(; t2 <= r2; t2++){
	  int J = q->sites2[t2];
	  if(!(/* WAS80 JL <= J && */ J <= JR)){
	    if(t2 > 0)
	      t2--;
	    break;
	  }
	}
	if(DEBUG) assert(qryendidx1 == q->sites2[r2]);
	if(DEBUG) assert(t2 < r2);
      } else {/* inverted */
	if(DEBUG) assert(qryendidx1 == M+1 - q->sites2[t2]);
	for(; t2 <= r2; r2--){
	  int J = M+1 - q->sites2[r2];
	  if(!(/* WAS80 JL <= J && */ J <= JR)){
	    if(r2 < q->numpairs - 1)
	      r2++;
	    break;
	  }
	}
	if(DEBUG) assert(t2 < r2);
      }
    } else {/* include label JL (or first aligned label >= JL) in NOT overlapped range */
      if(DEBUG) assert(qrystartidx1 < JL);
      if(!inverted){
	if(DEBUG) assert(qrystartidx1 == q->sites2[t2]);
	for(; t2 <= r2; r2--){
	  int J = q->sites2[r2];
	  if(!(JL <= J /* WAS80 && J <= JR */)){
	    if(r2 < q->numpairs - 1)
	      r2++;
	    break;
	  }
	}
	if(DEBUG) assert(t2 < r2);
      } else {// inverted 
	for(; t2 <= r2; t2++){
	  int J = M+1 - q->sites2[t2];
	  if(!(JL <= J /* WAS80 && J <= JR */)){
	    if(t2 > 0)
	      t2--;
	    break;
	  }
	}
	if(DEBUG) assert(qrystartidx1 == M+1 - q->sites2[r2]);
	if(DEBUG) assert(t2 < r2);
      }
    }
    if(CVERBOSE){
      printf("\t Non-overlaped region of xmapid2: t2=%d(I=%d,K=%d,J=%d) .. r2=%d(I=%d,K=%d,J=%d), np=%d\n",
	     t2,q->sites1[t2],q->sitesK1[t2], q->orientation ? M+1 - q->sites2[t2] : q->sites2[t2],
	     r2,q->sites1[r2],q->sitesK1[r2], q->orientation ? M+1 - q->sites2[r2] : q->sites2[r2], q->numpairs);
      fflush(stdout);
    }
  } else {
    t2 = tL[maxoutlier];
    r2 = tR[maxoutlier];

    if(CVERBOSE){
      printf("\t maxoutlier=%d/%d:usexmapentries[i]->xmapid_overlapped -> %d: t2=%d,r2=%d, np= %d\n",maxoutlier,outliers,usexmapentries[j]->xmapid,t2,r2,q->numpairs);
      fflush(stdout);
    }
    usexmapentries[i]->xmapid_overlapped = usexmapentries[j]->xmapid;

    if(VERB>=3){
      for(int i = 0; i < numxmapentries; i++){
	xmapEntry *p = usexmapentries[i];
	if(p->xmapid == 868 || p->xmapid == 866){
	  printf("\t usexmapentries[%d]: refid=%lld,qryid=%lld, xmapid=%d, xmapid_overlapped=%d\n",
		 i, p->refcontigid,p->qrycontigid,p->xmapid,p->xmapid_overlapped);
	  fflush(stdout);
	}
      }
    }
  }

  /* trim back usexmapentries[i] & q to aligned pairs t2 .. r2 */
  //  int origrefstartidx = refstartidx;
  //  int origrefendidx = refendidx;

  if(t2 > 0){/* trim left end */
    refstartidx = q->sites1[t2] - q->sitesK1[t2];
    if(inverted){
      qryendidx = q->orientation ? M+1 - q->sites2[t2] : q->sites2[t2];// absolute index of new "rightmost" MGi qry label with query oriented same as in usexmapentries[j]
      qryendidx1 = orientation ? (M+1 - qryendidx) : qryendidx;// relative index of rightmost MGi qry label with qry as oriented same as in usexmapentries[j] 
    } else {
      qrystartidx = q->orientation ? M+1 - q->sites2[t2] : q->sites2[t2];// absolute index of new "leftmost" MGi qry label with query oriented same as usexmapentries[j]
      qrystartidx1 = orientation ? (M+1 - qrystartidx) : qrystartidx;// leftmost MGi qry label with qry oriented as in usexmapentries[j]
    }

    matchcnt = q->numpairs - t2;

    if(CVERBOSE){
      printf("\t Trimming left end of xmapid2 : new Ref= %d..%d(%0.3f..%0.3f),Qry= %d..%d(%0.3f..%0.3f),np=%d\n", 
	     refstartidx,refendidx,Yi[refstartidx],Yi[refendidx],qrystartidx,qryendidx,X[qrystartidx],X[qryendidx],matchcnt);
      fflush(stdout);
    }
  }
  if(DEBUG) assert(0 <= t2 && t2 < r2 && r2 < q->numpairs);

  if(r2 < q->numpairs - 1){
    refendidx = q->sites1[r2];
    if(inverted){
      qrystartidx = q->orientation ? M+1 - q->sites2[r2] : q->sites2[r2];// absolute index of new "leftmost" MGi qry label with query oriented same as usexmapentries[j]
      qrystartidx1 = orientation ? (M+1 - qrystartidx) : qrystartidx;// leftmost MGi qry label with qry oriented as in usexmapentries[j] 
    } else {
      qryendidx = q->orientation ? M+1 - q->sites2[r2] : q->sites2[r2];// absolute index of new "rightmost" MGi qry label with query oriented same as usexmapentries[j]
      qryendidx1 = orientation ? (M+1 - qryendidx) : qryendidx;// rightmost MGi qry label with qry as oriented as in usexmapentries[j] 
    }
    matchcnt = r2 - t2 + 1;

    if(CVERBOSE){
      printf("\t Trimming right end of xmapid2 : new Ref= %d..%d(%0.3f..%0.3f),Qry= %d..%d(%0.3f..%0.3f),np=%d\n", 
	     refstartidx,refendidx,Yi[refstartidx],Yi[refendidx],qrystartidx,qryendidx,X[qrystartidx],X[qryendidx],matchcnt);
      fflush(stdout);
    }
  }

  if(t2 > 0 || r2 < q->numpairs - 1){/* confirm trimming of q and usexmapentries[i] */

    if(numsvIndel > 0){/* NEW82 : remove Indels with usexmapentries[i]->xmapid that no longer map to trimmed ref site interval refstartidx .. refstartidx */
      if(DEBUG)assert(usexmapentries[j]->refcontigid == usexmapentries[i]->refcontigid);
      long long rid = usexmapentries[i]->refcontigid;
      long long qid = usexmapentries[i]->refcontigid;
      int xmapid1 = usexmapentries[i]->xmapid;
		 
      if(DEBUG) assert(firstIndel != 0);
      if(DEBUG) assert(0 < xmapid1 && xmapid1 <= MaxXmapid);

      for(int sv = firstIndel[xmapid1]; sv < numsvIndel; sv++){
	structuralVariation *pSV = &sv_array[sv];
	if(DEBUG>=2) assert(pSV->algo_enum != sv_ps);
	xmapEntry *entry1 = xmapentries[pSV->indices[0]];
	if(entry1->xmapid != xmapid1)
	  break;

	if(DEBUG>=2) assert(pSV->refcontigid1 == pSV->refcontigid2);
	if(pSV->refcontigid1 != rid)
	  continue;
	if(pSV->confidence == 0.0)
	  continue;
	if(pSV->refstartidx <= refstartidx || pSV->refstopidx >= refendidx){/* single label overlap would have no confidence */
	  if(CVERBOSE){
	    printf("\t Filtered out Indel sv_array[%d]:rid=%lld,qid=%lld,fwd=%d:ref=%d..%d(%0.3f..%0.3f),qry=%d..%d(%0.3f..%0.3f):conf=%0.2f\n",
		   sv,rid,qid,entry1->orientforward,pSV->refstartidx,pSV->refstopidx,pSV->refstart*1e-3,pSV->refstop*1e-3,
		   pSV->querystartidx,pSV->querystopidx,pSV->querystart*1e-3,pSV->querystop*1e-3, pSV->confidence);
	    fflush(stdout);
	  }
	  pSV->confidence = 0.0;
	  continue;
	}
      }
    }

    ret = 1;

    if(t2 > 0){/* new left endoutlier */
      int I = q->sites1[t2];
      int K = q->sitesK1[t2];
      int J = q->sites2[t2];
      q->iscore[0] = q->outscore[0] = ChimScore + Sm(0,I,K,Yi);
      q->Lij1 = I-K;
      q->Lij2 = J;
      q->Lend = localtype;
    }

    int t = t2;
    int lastI = q->sites1[0] = q->sites1[t2];
    int lastK = q->sitesK1[0] = q->sitesK1[t2];
    int lastJ = q->sites2[0] = q->sites2[t2];
    q->noutliers = 0;
    q->maxoutlier = 0.0;
    q->maxoutlierLabels = 0;
    for(int u = 1; ++t <= r2; u++){
      q->iscore[u] = q->iscore[t];
      q->outscore[u] = q->outscore[t];
      int I = q->sites1[u] = q->sites1[t];
      int K = q->sitesK1[u] = q->sitesK1[t];
      int J = q->sites2[u] = q->sites2[t];
      if(q->outscore[u] + (FLOAT)0.01 < q->iscore[u]){
	q->noutliers++;
	double deltaY = Yc(Yi,I,K) - Yc(Yi,lastI,lastK);
	double deltaX = Xrev[J] - Xrev[lastJ];
	if(DEBUG) assert(deltaX > 0.0);
	double delta = fabs(deltaY-deltaX);
	q->maxoutlier = max(delta, q->maxoutlier);
	q->maxoutlierLabels = max(I-K-lastI + J-lastJ - 2, q->maxoutlierLabels);
      }
      lastI = I;
      lastK = K;
      lastJ = J;
    }

    int qU = r2 - t2 + 1;
    if(r2 < q->numpairs - 1){/* new right end */
      int I = q->sites1[r2];
      //      int K = q->sitesK1[r2];
      int J = q->sites2[r2];
      q->iscore[qU] = q->outscore[qU] = ChimScore;
      q->Rij1 = I;
      q->Rij2 = J;
      q->Rend = localtype;
    } else if(t2 > 0)
      q->iscore[qU] = q->outscore[qU] = q->iscore[q->numpairs];

    q->numpairs = qU;
    
    q->score = 0;
    for(int t = 0; t <= qU; t++)
      q->score += q->iscore[t];
      
    q->logPV = alignFP(q, Yi, Xrev, Ni, M, q->orientation, q->mapid1, q->mapid2, q->score, 0);
    q->logPV2 = max(q->logPV,q->logPV2);// NEW176
    if(DEBUG>=2) assert(q->logPV2 >= q->logPV);
    q->chimpair = 1;// flag as split matchgroup so it isn't filtered out by LogPvThreshold
    if(CVERBOSE){
      printf("\t\t conf -> %0.2f\n",q->logPV);
      fflush(stdout);
    }

    usexmapentries[i]->confidence = q->logPV;

    if(DEBUG) assert(matchcnt == q->numpairs);
    if(DEBUG) assert(refstartidx == q->sites1[0] - q->sitesK1[0]);
    if(DEBUG) assert(refendidx == q->sites1[qU-1]);
    if(DEBUG) assert(refstartidx <= refendidx);
    if(DEBUG && !(refendidx - refstartidx + 1 <= usexmapentries[i]->numrefsites)){
      printf("refstartidx= %d -> %d, refendidx= %d -> %d, numrefsites = %d -> %d\n",
	     usexmapentries[i]->refstartidx,refstartidx,usexmapentries[i]->refendidx,refendidx,usexmapentries[i]->numrefsites, refendidx - refstartidx + 1);
      fflush(stdout);
      assert(refendidx - refstartidx + 1 <= usexmapentries[i]->numrefsites);
    }
    usexmapentries[i]->refstartidx = refstartidx;
    usexmapentries[i]->refendidx = refendidx;

    if(inverted){
      usexmapentries[i]->qrystartidx = qryendidx;
      usexmapentries[i]->qryendidx = qrystartidx;
    } else {
      usexmapentries[i]->qrystartidx = qrystartidx;
      usexmapentries[i]->qryendidx = qryendidx;
    }
    if(DEBUG) assert(usexmapentries[i]->qrystartidx == (q->orientation ? M+1 - q->sites2[0] : q->sites2[0]));
    if(DEBUG) assert(usexmapentries[i]->qryendidx == (q->orientation ? M+1 - q->sites2[qU-1] : q->sites2[qU-1]));
    if(DEBUG) assert(usexmapentries[i]->orientforward != (q->orientation ? true : false));
    if(DEBUG/* HERE >=2 */) assert(usexmapentries[i]->orientforward ? (usexmapentries[i]->qrystartidx <= usexmapentries[i]->qryendidx) : 
				   (usexmapentries[i]->qrystartidx >= usexmapentries[i]->qryendidx) );

    usexmapentries[i]->refstartpos = Yi[refstartidx] * 1000.0;
    usexmapentries[i]->qrystartpos = X[usexmapentries[i]->qrystartidx] * 1000.0;
    usexmapentries[i]->refendpos = Yi[refendidx] * 1000.0;
    usexmapentries[i]->qryendpos = X[usexmapentries[i]->qryendidx] * 1000.0;

    usexmapentries[i]->startunalign = min(usexmapentries[i]->qrystartidx, usexmapentries[i]->qryendidx) - 1;
    usexmapentries[i]->endunalign = usexmapentries[i]->querycontignsite - max(usexmapentries[i]->qrystartidx, usexmapentries[i]->qryendidx);
    if(DEBUG) assert(M == usexmapentries[i]->querycontignsite);

    if(DEBUG) assert(usexmapentries[i]->numrefsites >= refendidx - refstartidx + 1);
    for(int t = refstartidx; t <= refendidx; t++){
      usexmapentries[i]->refsites[t-refstartidx] = Yi[t] * 1000.0;
      if(VERB>=2){
	printf("usexmapentries[%d]->refsites[%d] = %0.1f (xmapid=%d)\n", i, t-refstartidx, usexmapentries[i]->refsites[t-refstartidx], usexmapentries[i]->xmapid);
	fflush(stdout);
      }
      if(DEBUG>=2) assert(usexmapentries[i]->refstartpos <= Yi[t] * 1000.0  && Yi[t] * 1000.0 <= usexmapentries[i]->refendpos);
    }
    if(DEBUG) assert(usexmapentries[i]->numqrysites >= abs(usexmapentries[i]->qryendidx - usexmapentries[i]->qrystartidx));
    if(q->orientation){
      if(DEBUG) assert(usexmapentries[i]->qrystartidx >= usexmapentries[i]->qryendidx);
      if(DEBUG) assert(usexmapentries[i]->qrystartpos >= usexmapentries[i]->qryendpos);
      for(int t = usexmapentries[i]->qryendidx; t <= usexmapentries[i]->qrystartidx; t++){
	usexmapentries[i]->qrysites[t - usexmapentries[i]->qryendidx] = X[t] * 1000.0;
	if(DEBUG>=2) assert(usexmapentries[i]->qryendpos <= X[t] * 1000.0  && X[t] * 1000.0 <= usexmapentries[i]->qrystartpos);
      }
    } else {
      if(DEBUG) assert(usexmapentries[i]->qrystartidx <= usexmapentries[i]->qryendidx);
      if(DEBUG) assert(usexmapentries[i]->qrystartpos <= usexmapentries[i]->qryendpos);
      for(int t = usexmapentries[i]->qrystartidx; t <= usexmapentries[i]->qryendidx; t++){
	usexmapentries[i]->qrysites[t - usexmapentries[i]->qrystartidx] = X[t] * 1000.0;
	if(DEBUG>=2) assert(usexmapentries[i]->qrystartpos <= X[t] * 1000.0  && X[t] * 1000.0 <= usexmapentries[i]->qryendpos);
      }
    }
  }

  delete [] Xrev;

  return ret;
}

class CSmallDuplication
{
  /* NOTE : refers to a large matchgroup with a smaller non-inverted (or inverted) matchgroup that overlaps mostly Outliers in the Query and non-outlier region of the Reference */

public:
  long long refcontigid;/* reference contig id */
  long long qrycontigid;/* query contig id */
  xmapEntry *Large;/* large matchgroup with the Outlier */
  xmapEntry *Small;/* small matchgroup with the duplicate region */
  bool orientforward;/* orientation of Query in larger matchgroup (with Outlier) : same for both matchgroups */
  bool inverted;/* If smaller matchgroup has opposite orientation of larger matchgroup (Small Inverted duplication) */

  FLOAT *Y;
  int N;
  FLOAT *X;
  int M;

  int RefDupStart;/* start label of Duplication in Reference (from Small matchgroup) */
  int RefDupEnd;/* end label of Duplication in Reference (from Small matchgroup) */
  int QryDupStart;/* start label of Duplication in Query (from Small or Large matchgroup, whichever is leftmost, and if unaligned in Large matchgroup, the nearest surrounding label) */
  int QryDupEnd;/* end label of Duplication in Query (from Small or Large matchgroup, whichever is rightmost, and if unaligned in Large matchgroup, the nearest surrounding label) */

  double confBC;/* Bonferroni corrected confidence of smaller matchgroup */
  double LeftConf;/* confidence (LogPV) of left side of large matchgroup (excluding regions in Query or Reference of the smaller matchgroup) */
  double RightConf; /* confidence (LogPV) of right side of large matchgroup (excluding regions in Query or Reference of the smaller matchgroup) */
};

// Small Duplications will be appended to SmallDuplications[0..NumSmallDuplications-1] : return 1 if anything was appended

#define SVERBOSE 0 // 1 or 2

static int checkSmallDuplication(CSmallDuplication *&SmallDuplications, 
			       int &NumSmallDuplications, 
			       int &MaxSmallDuplications,
			       xmapEntry **usexmapentries, 
			       int j, // usexmapentries[j] is the larger (higher confidence) matchgroup 
			       int i // usexmapentries[i] is the smaller (lower confidence) matchgroup of opposite orientation
			       )
{
#if 0
  int refstartidx, refendidx, qrystartidx,qryendidx,qrystartidx1,qryendidx1,matchcnt,matchcnt1,matchcnt2;
  int refstartidx2 = -1,refendidx2 = -1,qrystartidx2 = -1,qryendidx2 = -1;/* reference and query label ranges for matchgroup j that most tightly bound matchgroup i on Query */
  int M,N,U;
  FLOAT *X,*Y;
  Calign *p,*q;
#endif  

  //  verify Small Duplication (with outlier in usexmapentries[j]) and save the Duplication information (4 boundary locations)
  int refstartidx = usexmapentries[i]->refstartidx;
  int refendidx = usexmapentries[i]->refendidx;
  if(DEBUG) assert(refstartidx <= refendidx);

  bool inverted = (usexmapentries[i]->orientforward != usexmapentries[j]->orientforward);
  int qrystartidx = inverted ? usexmapentries[i]->qryendidx : usexmapentries[i]->qrystartidx; // absolute index, with qry oriented as in usexmapentries[j] (larger matchgroup)
  int qryendidx = inverted ? usexmapentries[i]->qrystartidx : usexmapentries[i]->qryendidx; // absolute index, with qry oriented as in usexmapentries[j] (larger matchgroup)

  Calign *q = usexmapentries[i]->align;
  int matchcnt = q->numpairs;

  Calign *p = usexmapentries[j]->align;
  int U = p->numpairs;

  FLOAT *Y = usexmapentries[j]->RefMap->site[0];
  int N = usexmapentries[j]->refcontignsite;
  FLOAT *X = usexmapentries[j]->QryMap->site[0];
  int M = usexmapentries[i]->querycontignsite;

  int orientation = usexmapentries[j]->orientforward ? 0 : 1;

  if(SVERBOSE){
    printf("Checking for Small Duplication : refid=%lld,qryid=%lld,j=%d,xmapid1=%d(conf=%0.2f,np=%d),i=%d,xmapid2=%d(conf=%0.2f,U=%d),fwd=%d,%d:RefDup= %d..%d(%0.3f..%0.3f),QryDup= %d..%d(%0.3f..%0.3f), M=%d,N=%d\n",
	   usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,j,usexmapentries[j]->xmapid,usexmapentries[j]->confidence,U,i,usexmapentries[i]->xmapid,usexmapentries[i]->confidence,matchcnt,
	   usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[i]->orientforward ? 1 : 0, refstartidx,refendidx, Y[refstartidx],Y[refendidx],
	   qrystartidx,qryendidx,X[qrystartidx],X[qryendidx],M,N);
    fflush(stdout);
  }

  if(DEBUG/* HERE HERE >=2 */) assert(usexmapentries[i]->orientforward ? (usexmapentries[i]->qrystartidx <= usexmapentries[i]->qryendidx) : 
					   (usexmapentries[i]->qrystartidx >= usexmapentries[i]->qryendidx) );
  if(DEBUG/* HERE HERE >=2 */) assert(usexmapentries[j]->orientforward ? (usexmapentries[j]->qrystartidx <= usexmapentries[j]->qryendidx) : 
					   (usexmapentries[j]->qrystartidx >= usexmapentries[j]->qryendidx) );
  if(DEBUG/* HERE HERE >=2 */) assert(usexmapentries[i]->align->orientation == (usexmapentries[i]->orientforward ? 0 : 1));
  if(DEBUG/* HERE HERE >=2 */) assert(usexmapentries[j]->align->orientation == (usexmapentries[j]->orientforward ? 0 : 1));

  int qrystartidx1 = orientation ? (M+1 - qrystartidx) : qrystartidx;// leftmost Duplication label with qry oriented as in usexmapentries[j] (larger matchgroup)
  int qryendidx1 = orientation ? (M+1 - qryendidx) : qryendidx;// rightmost Duplication label with qry as oriented as in usexmapentries[j] (larger matchgroup)

  if(SVERBOSE){
    printf("\t refid=%lld,qryid=%lld:qry1=%d..%d,M=%d\n",usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, qrystartidx1,qryendidx1,M);
    fflush(stdout);
  }

  if(DEBUG && !(qrystartidx1 <= qryendidx1)){
    #pragma omp critical
    {
      printf("\t refid=%lld,qryid=%lld:qry1=%d..%d\n",usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, qrystartidx1,qryendidx1);
      fflush(stdout);
      assert(qrystartidx1 <= qryendidx1);
    }
  }

  int trimL = 0, trimR = matchcnt;/* range will be reduced if temporary trimming occurs : Trimming will be made permanant if Duplication is confirmed AND MATCHGROUP_TRIM > 0 */

  if(!inverted){/* trim smaller matchgroup to remove identical alignment pairs (in case this was not fixed in refalign) */
    /* locate leftmost matching alignments in p & q */
    int t1= 0, t2 = 0;
    for(; t1 < U && t2 < matchcnt; ){
      if(p->sites1[t1] < q->sites1[t2]){
	t1++;
	continue;
      }
      if(p->sites1[t1] > q->sites1[t2]){
	t2++;
	continue;
      }
      if(p->sites2[t1] < q->sites2[t2]){
	t1++;
	continue;
      }
      if(p->sites2[t1] > q->sites2[t2]){
	t2++;
	continue;
      }
      if(p->sitesK1[t1] == q->sitesK1[t2]) /* overlap match found */
	break;
	    
      if(t1 + 1 < U)
	t1++;
      else
	t2++;
    }

    if(SVERBOSE){
      printf("t1=%d/%d,t2=%d/%d\n",t1,U,t2,matchcnt);
      fflush(stdout);
    }

    if(t1 < U && t2 < matchcnt){/* overlap match found */
      /* locate rightmost contiguous match */
      int r1 = t1, r2 = t2;
      while(r1 < U - 1 && r2 < matchcnt - 1){
	if(p->sites1[r1+1] != q->sites1[r2+1])
	  break;
	if(p->sites2[r1+1] != q->sites2[r2+1])
	  break;
	if(p->sitesK1[r1+1] != q->sitesK1[r2+1])
	  break;
	r1++;
	r2++;
      }

      if(SVERBOSE){
	printf("t1=%d/%d,t2=%d/%d, r1=%d,r2=%d\n",t1,U,t2,matchcnt,r1,r2);
	fflush(stdout);
      }

      if(t2 == 0 && r2 < matchcnt - 1){/* trim left end locally */
	refstartidx = q->sites1[r2+1] - q->sitesK1[r2+1];
	qrystartidx1 = q->sites2[r2+1];
	qrystartidx = !usexmapentries[i]->orientforward /* WAS78 orientation */ ? (M+1 - qrystartidx1) : qrystartidx1;// absolute index
	matchcnt -= r2+1;

	trimL = max(trimL, r2+1);

	if(SVERBOSE){
	  printf("\t Trimming left end of xmapid2 : new RefDup= %d..%d(%0.3f..%0.3f),QryDup= %d..%d(%0.3f..%0.3f),np=%d\n", refstartidx,refendidx,Y[refstartidx],Y[refendidx],qrystartidx,qryendidx,X[qrystartidx],X[qryendidx],matchcnt);
	  fflush(stdout);
	}

      } else if(t2 > 0 && r2 == matchcnt - 1){/* trim right end locally */

	refendidx = q->sites1[t2 - 1];
	qryendidx1 = q->sites2[t2 - 1];
	qryendidx = !usexmapentries[i]->orientforward /* WAS78 orientation */ ? (M+1 - qryendidx1) : qryendidx1;// absolute index
	matchcnt = t2;

	trimR = min(trimR, t2);

	if(SVERBOSE){
	  printf("\t Trimming right end of xmapid2 : new RefDup= %d..%d(%0.3f..%0.3f),QryDup= %d..%d(%0.3f..%0.3f),np=%d\n", refstartidx,refendidx,Y[refstartidx],Y[refendidx],qrystartidx,qryendidx,X[qrystartidx],X[qryendidx],matchcnt);
	  fflush(stdout);
	}
      }
    }
  }

  if(SMALL_DUPLICATION_FIX){
    /* need to trim smaller matchgroup to remove query overlap, if possible, for both inverted or non-inverted case */
    /* First locate largest outlier in larger matchgroup j that overlaps matchgroup i on qry */

    int qL = M+1, qR = 0;/* range of labels in Query of smaller matchgroup i that are overlapped by outlier of larger matchgroup j (with query orientation matching larger matchgroup) */
    int lastJ = p->sites2[0];
    int maxOverlap = 0;/*number of sites in largest Outlier overlap with qry range of matchgroup i */
    for(int t = 1; t < U; t++){
      int J = p->sites2[t];
      int Overlap = min(J, qryendidx1) - max(lastJ, qrystartidx1);
      if(Overlap > maxOverlap){
	maxOverlap = Overlap;
	qL = max(lastJ,qrystartidx1);
	qR = min(J,qryendidx1);
      }
      lastJ = J;
    }

    if(qL < qR){/* the largest Outlier Qry range is qL .. qR (with qry oriented like larger matchgroup j) : try to trim matchgroup i outside of this range */
      int t2 = trimL, r2 = trimR - 1;
      if(qrystartidx1 < qL){/* Trim off any part of matchgroup i that is left of qL on Qry */
	while(t2 < r2){
	  int J = inverted ? M+1 - q->sites2[t2] : q->sites2[t2];
	  if(J < qL){
	    t2++;
	    continue;
	  }
	  J = inverted ? M+1 - q->sites2[r2] : q->sites2[r2];
	  if(J < qL){
	    r2--;
	    continue;
	  }
	  break;
	}
      }
      if(qryendidx1 > qR){/* Trim off any part of matchgroup i that is right of qR on Qry */
	while(t2 < r2){
	  int J = inverted ? M+1 - q->sites2[t2] : q->sites2[t2];
	  if(J > qR){
	    t2++;
	    continue;
	  }
	  J = inverted ? M+1 - q->sites2[r2] : q->sites2[r2];
	  if(J > qR){
	    r2--;
	    continue;
	  }
	  break;
	}
      }

      if(t2 > trimL){/* trim left end */
	refstartidx = q->sites1[t2] - q->sitesK1[t2];
	if(inverted){// NEW141
	  qryendidx1 = M+1 - q->sites2[t2];
	  qryendidx =  orientation ? (M+1 - qryendidx1) : qryendidx1;// absolute index
	} else {
	  qrystartidx1 = q->sites2[t2];
	  qrystartidx = orientation ? (M+1 - qrystartidx1) : qrystartidx1;// absolute index
	}
	trimL = t2;
	matchcnt = trimR - trimL;

	if(SVERBOSE){
	  printf("\t Trimming left end of xmapid2 : new RefDup= %d..%d(%0.3f..%0.3f),QryDup= %d..%d(%0.3f..%0.3f),trimL -> %d, trimR=%d,np=%d->%d\n", 
		 refstartidx,refendidx,Y[refstartidx],Y[refendidx],qrystartidx,qryendidx,X[qrystartidx],X[qryendidx],trimL,trimR,q->numpairs,matchcnt);
	  fflush(stdout);
	}
      }
      if(r2 < trimR - 1){/* trim right end */
	refendidx = q->sites1[r2];
	if(inverted){// NEW141
	  qrystartidx1 = M+1 - q->sites2[r2];
	  qrystartidx = orientation ? (M+1 - qrystartidx1) : qrystartidx1;// absolute index	  
	} else {
	  qryendidx1 = q->sites2[r2];
	  qryendidx = orientation ? (M+1 - qryendidx1) : qryendidx1;// absolute index
	}
	trimR = r2+1;
	matchcnt = trimR - trimL;

	if(SVERBOSE){
	  printf("\t Trimming right end of xmapid2 : new RefDup= %d..%d(%0.3f..%0.3f),QryDup= %d..%d(%0.3f..%0.3f),np=%d\n", refstartidx,refendidx,Y[refstartidx],Y[refendidx],qrystartidx,qryendidx,X[qrystartidx],X[qryendidx],matchcnt);
	  fflush(stdout);
	}
      }
    }
  }

  /* find the two aligned pairs in the Large matchgroup that bound the smaller matchgroup in the Query */
  int tL = -1, tR = -1;/* alignment index values that most closely bound the smaller matchgroup on the Query (not including overlapped identical alignment pairs) */
  int refstartidx2 = -1,refendidx2 = -1,qrystartidx2 = -1,qryendidx2 = -1;/* reference and query label ids corresponding to tL and tR */
  int matchcnt1 = 0;// Count of aligned pairs in the Large matchgroup (outside bounding range, including bounding pairs) that overlap the smaller matchgroup Reference Range : this should be high fraction of matchcnt
  int matchcnt2 = 0;// Count of aligned pairs in Large matchgroup (within bounding range, including bounding pairs) that overlap the smaller matchgroup Query Range : this should be small fraction of matchcnt

  for(int t = 0; t < U; t++){
    int I = p->sites1[t];
    int K = p->sitesK1[t];
    int J = p->sites2[t];

    // update tL,tR,refstartidx2,refendidx2,qrystartidx2,qryendidx2
    if(J <= qrystartidx1){/* update left bound */
      tL = t;
      if(DEBUG) assert(I >= refstartidx2);
      refstartidx2 = I;
      if(DEBUG) assert(J >= qrystartidx2);
      qrystartidx2 = J;
      if(J >= qrystartidx1)
	matchcnt2++;

      if(refstartidx <= I && I <= refendidx)
	matchcnt1++;
    } else if(J >= qryendidx1){/* at or beyond right bound */
      if(refendidx2 < 0){/* found right bound */
	tR = t;
	refendidx2 = I-K;
	qryendidx2 = J;
      }
      if(J <= qryendidx1)
	matchcnt2++;
      if(refstartidx <= I && I <= refendidx)
	matchcnt1++;
    } else /* between bounds */
      //    if(J >= qrystartidx1 && J <= qryendidx1)
      matchcnt2++;

    if(SVERBOSE>=2)
      printf("\t t=%d/%d:I=%d,K=%d,J=%d : ref2=%d..%d,qry2=%d..%d,matchcnt=%d, matchcnt1=%d,matchcnt2=%d,tL=%d,tR=%d\n",t, U, I,K,J, refstartidx2,refendidx2,qrystartidx2,qryendidx2,matchcnt, matchcnt1,matchcnt2,tL,tR);// HERE
  }

  if(SVERBOSE){
    printf("\t tL=%d,tR=%d,or=%d,inv=%d:refstartidx=%d,refendidx=%d,refstartidx2=%d,refendidx2=%d,qrystartidx2=%d,qryendidx2=%d,matchcnt1=%d,matchcnt2=%d,matchcnt=%d\n",
	   tL,tR,orientation,inverted?1:0,refstartidx,refendidx,refstartidx2,refendidx2,qrystartidx2,qryendidx2,matchcnt1,matchcnt2,matchcnt);
    fflush(stdout);
  }
  if(DEBUG && refstartidx2 >= 1 && refendidx2 >= 1){
    assert(refstartidx2 <= refendidx2);
    assert(qrystartidx2 <= qryendidx2);
  }

  if(refstartidx >= refendidx || matchcnt1 < smap_Dup_minRefOverlapSites)
    return 0;

  if(!(refstartidx2 >= 1 && refendidx2 >= 1 && matchcnt2 <= floor(matchcnt * SmallDupMaxOverlap + 0.5) /* overlap in query labels <= SmallDupMaxOverlap * matchcnt */ 
       && matchcnt1 >= floor(matchcnt * SmallDupMinOverlap + 0.5) /* overlap in ref labels  >= SmallDupMinOverlap * matchcnt */))
    return 0;

  int misscnt = 0;/* count number of aligned labels of small matchgroup in reference that do NOT match an aligned label in the larger matchgroup */
  int t1 = trimL, t2 = 0;
  for(; t1 < trimR && t2 < U; ){
    int I1 = q->sites1[t1];
    int K1 = I1 - q->sitesK1[t1];
    int I2 = p->sites1[t2];
    int K2 = I2 - p->sitesK1[t2];
    
    if(SVERBOSE)
      printf("\t t1=%d/%d, t2=%d/%d: I1=%d,K1=%d,J1=%d; I2=%d,K2=%d,J2=%d : matchcnt=%d, misscnt=%d\n", t1,matchcnt,t2,U, I1,I1-K1,q->sites2[t1],I2,I2-K2,p->sites2[t2], matchcnt, misscnt);
    
    if(I1 < K2){
      misscnt++;
      t1++;
      continue;
    }
    if(I2 < K1){
      t2++;
      continue;
    }
    /* at least a touching match of q with p */
    t1++;
    t2++;
  }
 
  if(SVERBOSE){
    printf("\t Final: matchcnt=%d, misscnt=%d (SmallDupMaxMiss= %0.4f, %d)\n",   matchcnt, misscnt, SmallDupMaxMiss, (int)floor(matchcnt * SmallDupMaxMiss + 0.5));  
    fflush(stdout);
  }

  if(misscnt > floor(matchcnt * SmallDupMaxMiss + 0.5))
    return 0;

  /* remap qrystartidx2 and qryendidx2 to absolute coordinates */
  if(orientation){
    qrystartidx2 = M+1 - qrystartidx2;
    qryendidx2 = M+1 - qryendidx2;
  }

  /* map refstartidx,refendidx to query using larger matchgroup */
  double qrystartloc3, qryendloc3;
  int qrystartidx3, qryendidx3;
  if(DEBUG && !(refstartidx < refendidx)){
    printf("WARNING: refid=%lld,qryid=%lld,or=%d,inv=%d:xmapid1=%d,xmapid2=%d,fwd=%d,%d:refstartidx=%d,refendidx=%d,qrystartidx=%d,qryendidx=%d: refstartidx2=%d,refendidx2=%d; qrystartidx2=%d,qryendidx2=%d\n",
	   usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,orientation,inverted?1:0,usexmapentries[i]->xmapid,usexmapentries[j]->xmapid,usexmapentries[i]->orientforward,usexmapentries[j]->orientforward,
	   refstartidx,refendidx,qrystartidx,qryendidx,refstartidx2,refendidx2,qrystartidx2,qryendidx2);
    fflush(stdout);
    assert(refstartidx < refendidx);
  }
  getQryCoord(p,refstartidx,refendidx, qrystartloc3, qrystartidx3, qryendloc3, qryendidx3);

  /* compute gap between both (absolute label) intervals in query : (qrystartidx,qryendidx), (qrystartidx3,qryendidx3) */
  int minidx1 = min(qrystartidx,qryendidx);
  int maxidx1 = max(qrystartidx,qryendidx);
  int minidx3 = min(qrystartidx3,qryendidx3);
  int maxidx3 = max(qrystartidx3,qryendidx3);

  int gapL = min(maxidx1,maxidx3);
  int gapR = max(minidx1,minidx3);

  int gap = max(minidx1 - maxidx3, minidx3 - maxidx1);
  if(SVERBOSE){
    printf("\t qryidx=%d..%d,qryidx2=%d..%d,qryidx3=%d..%d, gapL=%d,gapR=%d,gap=%d\n",qrystartidx,qryendidx,qrystartidx2,qryendidx2,qrystartidx3,qryendidx3,gapL,gapR,gap);
    fflush(stdout);
  }

  if(DEBUG && !inverted && !(gap == gapR - gapL)){// following assertion can fail if either matchgroup has large internal outliers
    printf("WARNING: refid=%lld,qryid=%lld,or=%d,inv=%d:xmapid1=%d,xmapid2=%d,fwd=%d,%d:refstartidx=%d,refendidx=%d,qrystartidx=%d,qryendidx=%d; qrystartidx3=%d,qryendidx3=%d: gap=%d, gapL=%d,gapR=%d : refstartidx2=%d,refendidx2=%d; qrystartidx2=%d,qryendidx2=%d\n",
	   usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,orientation,inverted?1:0,usexmapentries[i]->xmapid,usexmapentries[j]->xmapid,usexmapentries[i]->orientforward,usexmapentries[j]->orientforward,
	   refstartidx,refendidx,qrystartidx,qryendidx,qrystartidx3,qryendidx3,gap,gapL,gapR,refstartidx2,refendidx2,qrystartidx2,qryendidx2);
    fflush(stdout);
    //    assert(gap == gapR - gapL);
  }

  if(/* WAS36 */0&& gap > 0 && !inverted){ /* expand the reference interval and one of the query intervals by either (refstartidx2,qrystartidx2) OR (refendidx2,qryendidx2), whichever falls inside gap on query */
    if(gapL <= qrystartidx2 && qrystartidx2 <= gapR){
      if(DEBUG) assert(refstartidx2 >= refendidx);
      gapL = qrystartidx2;
      refendidx = refstartidx2;
    } else {
      if(DEBUG && !(gapL <= qryendidx2 && qryendidx2 <= gapR)){
	printf("refid=%lld,qryid=%lld,or=%d,inv=%d:refstartidx=%d,refendidx=%d,qrystartidx=%d,qryendidx=%d; qrystartidx3=%d,qryendidx3=%d: gapL=%d,gapR=%d : refstartidx2=%d,refendidx2=%d; qrystartidx2=%d,qryendidx2=%d\n",
	       usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,orientation,inverted?1:0,refstartidx,refendidx,qrystartidx,qryendidx,qrystartidx3,qryendidx3,gapL,gapR,refstartidx2,refendidx2,qrystartidx2,qryendidx2);
	fflush(stdout);
	assert(gapL <= qryendidx2 && qryendidx2 <= gapR);
      }
      if(DEBUG) assert(refendidx2 <= refstartidx);
      gapR = qryendidx2;
      refstartidx = refendidx2;
    }
  }

  double confBC = q->logPV - log((double)max(1,gap))/log(10.0);
  if(SVERBOSE){
    printf("\t qrystartidx=%d,qryendidx=%d,qrystartidx3=%d,qryendidx3=%d:gap=%d,confBC= %0.2f\n", qrystartidx,qryendidx,qrystartidx3,qryendidx3,gap, confBC);
    fflush(stdout);
  }

  if(confBC < SmallDupConfBC)
    return 0;

  /* Small Duplication is confirmed */
  if(MATCHGROUP_TRIM && (trimL > 0 || trimR < q->numpairs - 1)){/* trim the smaller matchgroup permanantly to q->sites1[trimL .. trimR-1] */
    for(int t = trimL; t < trimR; t++){
      q->sites1[t - trimL] = q->sites1[t];
      q->sitesK1[t - trimL] = q->sitesK1[t];
      q->sites2[t - trimL] = q->sites2[t];
    }
    if(DEBUG) assert(matchcnt == trimR - trimL);
    q->numpairs = trimR - trimL;
    
    usexmapentries[i]->refstartidx = refstartidx;
    usexmapentries[i]->refendidx = refendidx;
    if(inverted){// NEW141
      usexmapentries[i]->qryendidx = qrystartidx;
      usexmapentries[i]->qrystartidx = qryendidx;
    } else {
      usexmapentries[i]->qrystartidx = qrystartidx;
      usexmapentries[i]->qryendidx = qryendidx;
    }

    if(DEBUG/* HERE HERE >=2 */ && !(usexmapentries[i]->orientforward ? (usexmapentries[i]->qrystartidx <= usexmapentries[i]->qryendidx) : 
					  (usexmapentries[i]->qrystartidx >= usexmapentries[i]->qryendidx) )){
      printf("While checking for Small Duplication : refid=%lld,qryid=%lld,j=%d,xmapid1=%d(conf=%0.2f,np=%d),i=%d,xmapid2=%d(conf=%0.2f,U=%d),fwd=%d,%d\n",
	     usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,j,usexmapentries[j]->xmapid,usexmapentries[j]->confidence,U,i,usexmapentries[i]->xmapid,usexmapentries[i]->confidence,matchcnt,
	     usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[i]->orientforward ? 1 : 0);
      printf(" trimL=%d,trimR=%d: q->sites1[trimL]= %d, q->sites1[trimR-1]= %d, N= %d, q->sites2[trimL]= %d, q->sites2[trimR]= %d, M= %d\n",
	     trimL,trimR,q->sites1[trimL],q->sites1[trimR-1],N,q->sites2[trimL],q->sites2[trimR-1],M);
      printf(" refstartidx -> %d, refendidx -> %d, qrystartidx -> %d, qryendidx -> %d, orientforward = %d\n",refstartidx,refendidx,qrystartidx,qryendidx,usexmapentries[i]->orientforward);
      fflush(stdout);

      assert(usexmapentries[i]->orientforward ? (usexmapentries[i]->qrystartidx <= usexmapentries[i]->qryendidx) : 
	     (usexmapentries[i]->qrystartidx >= usexmapentries[i]->qryendidx) );
    }

    usexmapentries[i]->refstartpos = usexmapentries[i]->RefMap->site[0][refstartidx];
    usexmapentries[i]->qrystartpos = usexmapentries[i]->QryMap->site[0][qrystartidx];
    usexmapentries[i]->refendpos = usexmapentries[i]->RefMap->site[0][refendidx];
    usexmapentries[i]->qryendpos = usexmapentries[i]->QryMap->site[0][qryendidx];

    usexmapentries[i]->startunalign = min(usexmapentries[i]->qrystartidx, usexmapentries[i]->qryendidx) - 1;
    usexmapentries[i]->endunalign = usexmapentries[i]->querycontignsite - max(usexmapentries[i]->qrystartidx, usexmapentries[i]->qryendidx);

    if(SVERBOSE){
      printf("\t Confirming trimming of xmapid2 : new RefDup= %d..%d(%0.3f..%0.3f),QryDup= %d..%d(%0.3f..%0.3f),np=%d\n", refstartidx,refendidx,Y[refstartidx],Y[refendidx],qrystartidx,qryendidx,X[qrystartidx],X[qryendidx],matchcnt);
      fflush(stdout);
    }
  }

  if(NumSmallDuplications >= MaxSmallDuplications){
    MaxSmallDuplications *= 2;
    CSmallDuplication *newSmallDuplications = new CSmallDuplication[MaxSmallDuplications];
    memcpy(newSmallDuplications,SmallDuplications,NumSmallDuplications * sizeof(CSmallDuplication));
    delete [] SmallDuplications;
    SmallDuplications = newSmallDuplications;
  }

  CSmallDuplication *pDup = &SmallDuplications[NumSmallDuplications++];
  pDup->refcontigid = usexmapentries[j]->refcontigid;
  pDup->qrycontigid = usexmapentries[j]->qrycontigid;
  pDup->orientforward = usexmapentries[j]->orientforward;
  pDup->inverted = inverted;
	      
  pDup->Large = usexmapentries[j];
  pDup->Small = usexmapentries[i];

  if(DEBUG) assert(usexmapentries[j]->QryMap && usexmapentries[j]->QryMap == usexmapentries[i]->QryMap);
  if(DEBUG) assert(usexmapentries[j]->RefMap && usexmapentries[j]->RefMap == usexmapentries[i]->RefMap);
  pDup->Y = Y;
  pDup->N = N;
  pDup->X = X;
  pDup->M = M;

  pDup->confBC = confBC;

  pDup->RefDupStart = refstartidx;
  pDup->RefDupEnd = refendidx;
  pDup->QryDupStart = pDup->orientforward ? min(minidx1,minidx3) : max(maxidx1,maxidx3);
  pDup->QryDupEnd = pDup->orientforward ? max(maxidx1,maxidx3) : min(minidx1,minidx3);

  /* left confidence is based on alignment p->sites1[t = 0..tL], but skip Duplication region */
  extern double SmRA(int J, int I, int K, FLOAT *Y);
  double Leftscore = SmRA(p->sites2[0], p->sites1[0], p->sitesK1[0], Y);
  for(int t = 1; t <= tL; t++){
    if(p->sites1[t] > refstartidx)
      break;
    Leftscore += p->iscore[t];
  }
  pDup->LeftConf = Leftscore / log(10.0);
  
  /* right confidence is based on alignment p->sites1[t = tR .. U-1] but skip Duplication region */
  double Rightscore = SmRA(p->sites2[U-1], p->sites1[U-1], p->sitesK1[U-1], Y);
  for(int t = tR+1; t < U; t++){
    if(p->sites1[t] <= refendidx)
      continue;
    Rightscore += p->iscore[t];
  }
  pDup->RightConf = Rightscore / log(10.0);

  int QryOutlierStart = orientation ? M + 1 - qrystartidx2 : qrystartidx2;
  int QryOutlierEnd = orientation ? M + 1 - qryendidx2 : qryendidx2;

  if(VERB){
    if(inverted)
      printf("Small Inverted Duplication: refid=%lld,qryid=%lld, xmapid1= %d, xmapid2= %d, forwd= %d,inv=%d:RefOutlier= %d .. %d(%0.3f..%0.3f), QryOutlier= %d .. %d(%0.3f..%0.3f), RefDup= %d .. %d(%0.3f..%0.3f), QryInvDup= %d .. %d(%0.3f..%0.3f),cnt=%d,%d,%d,%d ConfBC= %0.2f (conf= %0.2f, %0.2f)\n",
	     pDup->refcontigid,pDup->qrycontigid,usexmapentries[j]->xmapid,usexmapentries[i]->xmapid,pDup->orientforward ? 1 : 0, pDup->inverted ? 1 : 0, refstartidx2,refendidx2,Y[refstartidx2],Y[refendidx2],
	     QryOutlierStart,QryOutlierEnd,X[QryOutlierStart],X[QryOutlierEnd],
	     pDup->RefDupStart, pDup->RefDupEnd, Y[pDup->RefDupStart], Y[pDup->RefDupEnd], pDup->QryDupStart,pDup->QryDupEnd, X[pDup->QryDupStart], X[pDup->QryDupEnd], 
	     matchcnt, matchcnt1,matchcnt2,misscnt,confBC,pDup->LeftConf,pDup->RightConf);
    else
      printf("Small Duplication: refid=%lld,qryid=%lld, xmapid1= %d, xmapid2= %d, forwd= %d,inv=%d:RefOutlier= %d .. %d(%0.3f..%0.3f), QryOutlier= %d .. %d(%0.3f..%0.3f), RefDup= %d .. %d(%0.3f..%0.3f), QryDup= %d .. %d(%0.3f..%0.3f),cnt=%d,%d,%d,%d ConfBC= %0.2f (conf= %0.2f, %0.2f)\n",
	     pDup->refcontigid,pDup->qrycontigid,usexmapentries[j]->xmapid,usexmapentries[i]->xmapid,pDup->orientforward ? 1 : 0, pDup->inverted ? 1 : 0, refstartidx2,refendidx2,Y[refstartidx2],Y[refendidx2],
	     QryOutlierStart,QryOutlierEnd,X[QryOutlierStart],X[QryOutlierEnd],
	     pDup->RefDupStart, pDup->RefDupEnd, Y[pDup->RefDupStart], Y[pDup->RefDupEnd], pDup->QryDupStart,pDup->QryDupEnd, X[pDup->QryDupStart], X[pDup->QryDupEnd], 
	     matchcnt, matchcnt1,matchcnt2,misscnt,confBC,pDup->LeftConf,pDup->RightConf);
    fflush(stdout);
  }

  return 1;
}

class CSmallInversion
{
  /* NOTE : refers to a large matchgroup with a smaller inverted matchgroup that overlaps mostly Outliers in both Reference and Query */

public:
  long long refcontigid;/* reference contig id */
  long long qrycontigid;/* query contig id */
  xmapEntry *Large;/* large matchgroup with the Outlier */
  xmapEntry *Small;/* small matchgroup with the Inversion */
  bool orientforward;/* orientation of Query in larger matchgroup (with Outlier). For Inversion Matchgroup Query orientation is the opposite */

  FLOAT *Y;
  int N;
  FLOAT *X;
  int M;

  int RefOutlierStart;/* start label of Outlier in Reference */
  int RefOutlierEnd;/* end label of outlier in Reference */
  int QryOutlierStart;/* start label of Outlier in Query */
  int QryOutlierEnd; /* end label of Outlier in Query */

  int RefInversionStart;/* start label of Inversion in Reference */
  int RefInversionEnd;/* end label of Inversion in Reference */
  int QryInversionStart;/* start label of Inversion in Query */
  int QryInversionEnd;/* end label of Inversion in Query */

  double LeftConf;/* confidence (LogPV) of left side of large matchgroup */
  double RightConf; /* confidence (LogPV) of right side of large matchgroup */
};

#define INV_VERB 0

/* Computes qry aligned label overlap matchcnt2/matchcnt after trimming any inverted palindromic ends of smaller matchgroup i that overlap larger matchgroup in qry
   Return 0, if inversion has been ruled out.
   Return 1, if inversion is still possible.
*/
static int SmallInversionOverlap(xmapEntry **usexmapentries,
				 int j, // usexmapentries[j] is the larger (higher confidence) matchgroup 
				 int i, // usexmapentries[i] is the smaller (lower confidence) matchgroup of opposite orientation
				 Calign* &p,
				 Calign* &q,
				 int &M, int &N, 
				 FLOAT* &X, FLOAT* &Y,
				 int &U,
				 int &refstartidx, int &refendidx,/* ref label range of matchgroup i */
				 int &qrystartidx, int &qryendidx,/* qry label range of matchgroup i : qryendidx is aligned with refstartidx
								     AND qrystartidx < qryendidx IFF matchgroup i has query oriented backwards */
				 int &qrystartidx1, int &qryendidx1, /* qry label range of matchgroup i with qry oriented same as matchgroup j : 
									qryendidx1 is aligned with refstartidx AND qrystartidx1 < qryendidx1 */
				 int &matchcnt, /* aligned qry labels on matchgroup i (after trimming Palindromic ends) */
				 int &matchcnt1, /* If InvDup : aligned labels on matchgroup j that overlap matchgroup i Ref Range  */
				 int &matchcnt2, /* aligned labels on matchgroup j that overlap matchgroup i Query and Ref ranges (or just Query range if InvDup) */
				 int &refstartidx2, int &refendidx2, /* ref label sub-range of matchgroup j that bound matchgroup i */
				 int &qrystartidx2, int &qryendidx2, /* qry label sub-range of matchgroup j that correspond to alignment range tL..tR and bounds qry range of match group i
									(with qry oriented as in matchgroup j) */
				 int &tL, int &tR,
				 bool InvDup = false
				 )
{
  //  verify inversion (with outlier in usexmapentries[j]) and save the inversion information (8 boundary locations)
  refstartidx = usexmapentries[i]->refstartidx;
  refendidx = usexmapentries[i]->refendidx;
  if(DEBUG) assert(refstartidx <= refendidx);

  qrystartidx = usexmapentries[i]->qryendidx; // absolute index // WAS min(usexmapentries[i]->qrystartidx, usexmapentries[i]->qryendidx);
  qryendidx = usexmapentries[i]->qrystartidx; // absolute index // WAS max(usexmapentries[i]->qrystartidx, usexmapentries[i]->qryendidx);
  // WAS  if(DEBUG) assert(qrystartidx <= qryendidx);

  q = usexmapentries[i]->align;
  matchcnt = q->numpairs;

  p = usexmapentries[j]->align;
  U = p->numpairs;

  Y = usexmapentries[j]->RefMap->site[0];
  N = usexmapentries[j]->refcontignsite;
  X = usexmapentries[j]->QryMap->site[0];
  M = usexmapentries[i]->querycontignsite;

  int orientation = usexmapentries[j]->orientforward ? 0 : 1;

  if(INV_VERB){
    printf("Checking for Small Inversion: refid=%lld(%lld),qryid=%lld(%lld),j=%d,xmapid1=%d(conf=%0.2f),i=%d,xmapid2=%d(conf=%0.2f),fwd1=%d:RefInversion= %d..%d(%0.3f..%0.3f),QryInversion= %d..%d(%0.3f..%0.3f), M=%d,N=%d,numpairs=%d\n",
	   usexmapentries[j]->refcontigid,usexmapentries[i]->refcontigid,usexmapentries[j]->qrycontigid,usexmapentries[i]->qrycontigid, j,usexmapentries[j]->xmapid,usexmapentries[j]->confidence,
	   i,usexmapentries[i]->xmapid,usexmapentries[i]->confidence,usexmapentries[j]->orientforward ? 1 : 0, refstartidx,refendidx, Y[refstartidx],Y[refendidx],
	   qrystartidx,qryendidx,X[qrystartidx],X[qryendidx],M,N,U);
    if(INV_VERB>=2){
      printf("\t p->mapid1=%d(id=%lld),p->mapid2=%d(id=%lld),q->mapid1=%d(id=%lld),q->mapid2=%d(id=%lld)\n",
	     p->mapid1,YYmap[p->mapid1]->id,p->mapid2,XXmap[p->mapid2]->id,q->mapid1,YYmap[q->mapid1]->id,q->mapid2,XXmap[q->mapid2]->id);
      printf("\t usexmapentries[j,i]->RefMap->mapid=%d,%d, usexmapentries[j,i]->QryMap->mapid=%d,%d\n",
	     usexmapentries[j]->RefMap->mapid, usexmapentries[i]->RefMap->mapid, usexmapentries[j]->QryMap->mapid, usexmapentries[i]->QryMap->mapid);
    }
    fflush(stdout);
  }

  if(DEBUG/* HERE HERE >=2 */) assert(usexmapentries[i]->orientforward ? (usexmapentries[i]->qrystartidx <= usexmapentries[i]->qryendidx) : 
					   (usexmapentries[i]->qrystartidx >= usexmapentries[i]->qryendidx) );
  if(DEBUG/* HERE HERE >=2 */) assert(usexmapentries[j]->orientforward ? (usexmapentries[j]->qrystartidx <= usexmapentries[j]->qryendidx) : 
					   (usexmapentries[j]->qrystartidx >= usexmapentries[j]->qryendidx) );
  if(DEBUG/* HERE HERE >=2 */) assert(usexmapentries[i]->align->orientation == (usexmapentries[i]->orientforward ? 0 : 1));
  if(DEBUG/* HERE HERE >=2 */) assert(usexmapentries[j]->align->orientation == (usexmapentries[j]->orientforward ? 0 : 1));

  qrystartidx1 = orientation ? (M+1 - qrystartidx /* WAS qryendidx */) : qrystartidx;// leftmost inversion label with qry oriented as in usexmapentries[j] (larger matchgroup)
  qryendidx1 = orientation ? (M+1 - qryendidx /* WAS qrystartidx */) : qryendidx;// rightmost inversion label with qry oriented as in usexmapentries[j] (larger matchgroup)

  if(DEBUG/* HERE >=2 */ && !(usexmapentries[i]->refcontigid == usexmapentries[j]->refcontigid && usexmapentries[i]->qrycontigid == usexmapentries[j]->qrycontigid &&
					refstartidx == q->sites1[0]  - q->sitesK1[0] && refendidx == q->sites1[matchcnt-1] && 
					qrystartidx1 == M+1 - q->sites2[matchcnt-1] && qryendidx1 == M+1 - q->sites2[0])){
    printf("refid=%lld,qryid=%lld,xmapid1=%d(conf=%0.2f),xmapid2=%d(conf=%0.2f),fwd=%d,%d:U2=%d:RefInv=%d..%d,QryInv=%d..%d(qry1=%d..%d),q->sites1[0]=%d,q->sitesK1[0]=%d,M=%d\n",
	   usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,j,usexmapentries[j]->confidence,i,usexmapentries[i]->confidence,usexmapentries[j]->orientforward,usexmapentries[i]->orientforward,
	   matchcnt,refstartidx,refendidx,qrystartidx,qryendidx,qrystartidx1,qryendidx1,q->sites1[0],q->sitesK1[0],M);
    printf("\t refstartidx=%d, refendidx= %d, matchcnt=%d, q->sites1[%d]= %d\n",refstartidx, refendidx, matchcnt, matchcnt-1, q->sites1[matchcnt-1]);
    printf("\t qrystartidx1= %d, M=%d, numpairs=%d, matchcnt= %d, q->sites2[%d]= %d\n",qrystartidx, M, q->numpairs, matchcnt, matchcnt-1, q->sites2[matchcnt-1]);
    printf("\t qryendidx1= %d, M=%d, q->sites2[0]= %d\n", qryendidx1, M, q->sites2[0]);
    fflush(stdout);
    
    assert(usexmapentries[i]->refcontigid == usexmapentries[j]->refcontigid);
    assert(usexmapentries[i]->qrycontigid == usexmapentries[j]->qrycontigid);
    assert(refstartidx == q->sites1[0]  - q->sitesK1[0]);
    assert(refendidx == q->sites1[matchcnt-1]);
    assert(qrystartidx1 == M+1 - q->sites2[matchcnt-1]);
    assert(qryendidx1 == M+1 - q->sites2[0]);
  }

  if(INV_VERB>=2){
    printf("\t refid=%lld,qryid=%lld:qry1=%d..%d\n",usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, qrystartidx1,qryendidx1);
    fflush(stdout);
  }

  if(DEBUG && !(qrystartidx1 <= qryendidx1)){// NEW
    #pragma omp critical
    {
      printf("\t refid=%lld,qryid=%lld:qry1=%d..%d\n",usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, qrystartidx1,qryendidx1);
      fflush(stdout);
      assert(qrystartidx1 <= qryendidx1);
    }
  }

  if(SMALL_DUPLICATION_FIX && InvDup){/* Try to trim inverted matchgroup locally on either end to reduce query overlap with larger matchgroup */
    // HERE HERE 
  }

  if(SMALL_INVERSION_FIX && !InvDup){  /* Try to trim inverted matchgroup locally on either end to remove palindromic extensions.
					  NOTE : Just adjust index ranges (refstartidx,refendidx,qrystartidx,qryendidx,qrystartidx1,qryendidx1) and matchcnt (remaining aligned labels in matchgroup i)
					  since alignment q == usematchgroup[i]->align is not used below this block */
    /* NOTE : the following algorithm is not optimal since it is greedy and may fail to remove the entire palindromic extension if the larger matchgroup has a missing alignment outside the
       outlier region : Should first identify the largest outlier region in the larger matchgroup that overlaps smaller matchgroup in both ref and qry. See MATCHGROUP_TRIM for a solution */

    int U2 = matchcnt; // q->numpairs

    int tL2 = 0, tR2 = U2 - 1;
    int backtrack = 0;/* remember if last change to tL2,tR2 involved changing both simultaneously : last such change will be undone to maximize size of inversion region */
    for(int t2 = tL2; t2 <= tR2; t2++){
      int I2 = q->sites1[t2];
      int K2 = q->sitesK1[t2];
      int J2 = M+1 - q->sites2[t2];/* orient same as larger matchgroup */

      int I1= -1, K1= -1, J1 = -1;
      /* check if there is an alignment from I2-K2 .. I2 to query via larger matchgroup, that is outside the range qrystartidx1 .. qryendidx1 */
      for(int t = 0; t < U; t++){
	int I = p->sites1[t];
	int K = p->sitesK1[t];
	int J = p->sites2[t];

	if(max(I-K,I2-K2) <= min(I,I2)){/* found aligned pair */
	  J1 = J;
	  break;
	}
	if(I-K >= I2)
	  break;
      }

      if(J1 >= 0 && J1 < qrystartidx1 && (matchcnt >= 2 || !backtrack)){/* delete aligned pair t2 since it it will reduce matchcnt by 1 and matchcnt2 by at least 1 */
	if(INV_VERB>=2){
	  printf("t2=%d/%d:I2=%d,K2=%d,J2=%d,tL2=%d,tR2=%d,:J1=%d,qrystartidx1= %d: deleting aligned pair t2 from left end of Inverted matchgroup\n",
		 t2,U2,I2,K2,J2,tL2,tR2,J1,qrystartidx1);
	  fflush(stdout);
	}

	if(DEBUG) assert(tL2 == t2);
	tL2 = t2 + 1;
	backtrack = 0;
	if(tL2 < U2){// NEW11
	  refstartidx = q->sites1[tL2] - q->sitesK1[tL2];
	  qryendidx1 = M+1 - q->sites2[tL2];
	} else {
	  refstartidx = q->sites1[U2-1] - q->sitesK1[U2-1] + 1;
	  qryendidx1 = M+1 - (q->sites2[U2-1] + 1);
	}
	if(DEBUG) assert(refstartidx <= N+1 && qryendidx1 >= 0);
	qryendidx = orientation ? (M+1 - qryendidx1) : qryendidx1;

	matchcnt--;
	continue;/* try to delete additional aligned pairs t2+1 etc */
      }

      /* check if there is an alignment from ref to J2 via larger matchgroup, that is outside the range refstartidx .. refendidx */
      for(int t = U; --t >= 0;){
	int I = p->sites1[t];
	int K = p->sitesK1[t];
	int J = p->sites2[t];
      
	if(J == J2){/* found aligned pair */
	  I1 = I;
	  K1 = K;
	  break;
	}
	if(J < J2)
	  break;
      }

      if(I1 - K1 > refendidx && (matchcnt >= 2 || !backtrack)){ /* delete aligned pair t2 since it it will reduce matchcnt by 1 and matchcnt2 by at least 1 */
	if(INV_VERB>=2){
	  printf("t2=%d/%d:I2=%d,K2=%d,J2=%d,tL2=%d,tR2=%d:I1=%d,refendidx= %d: deleting aligned pair t2 from left end Inverted matchgroup\n",
		 t2,U2,I2,K2,J2,tL2,tR2,I1,refendidx);
	  fflush(stdout);
	}
	if(DEBUG) assert(tL2 == t2);
	tL2 = t2 + 1;
	backtrack = 0;
	if(tL2 < U2){// NEW11
	  refstartidx = q->sites1[tL2] - q->sitesK1[tL2];
	  qryendidx1 = M+1 - q->sites2[tL2];
	} else {
	  refstartidx = q->sites1[U2-1] - q->sitesK1[U2-1] + 1;
	  qryendidx1 = M+1 - (q->sites2[U2-1] + 1);
	}
	if(DEBUG) assert(refstartidx <= N+1 && qryendidx1 >= 0);
	qryendidx = orientation ? (M+1 - qryendidx1) : qryendidx1;

	matchcnt--;
	continue;/* try to delete additional aligned pairs t2+1 etc */
      }

      /* check if we can delete both t2 and tR2 reducing both matchcnt and matchcnt2 by 2 : will revert last such change to keep inversion region as large as plausible */
      if(J1 >= 0 && J1 == qrystartidx1 && I1 >= refendidx && I1-K1 <= refendidx && tR2 >= tL2 + 1){
	if(INV_VERB>=2){
	  printf("t2=%d/%d,tL2=%d,tR2=%d:I2=%d,K2=%d,J2=%d:J1=%d,I1=%d,K1=%d,qryidx1=%d..%d, refidx= %d..%d: deleting aligned pairs tL2,tR2 from both ends of Inverted matchgroup\n",
		 t2,U2,tL2,tR2,I2,K2,J2,J1,I1,K1,qrystartidx1,qryendidx1,refstartidx,refendidx);
	  fflush(stdout);
	}
	if(DEBUG) assert(tL2 == t2);
	backtrack = 1;
	tL2 = tL2 + 1;// NEW11
	tR2 = tR2 - 1;
        if(DEBUG)assert(0 <= tR2 && tL2 < U2);
	refstartidx = q->sites1[tL2] - q->sitesK1[tL2];
	qryendidx1 = M+1 - q->sites2[tL2];
	qryendidx = orientation ? (M+1 - qryendidx1) : qryendidx1;
	refendidx = q->sites1[tR2];
	qrystartidx1 = M+1 - q->sites2[tR2];
	qrystartidx = orientation ? (M+1 - qrystartidx1) : qrystartidx1;

	matchcnt -= 2;
	continue;/* try to delete additional aligned pairs */
      }
      

      if(DEBUG) assert(tL2 == t2);
      // WAS11      tL2 = t2;
      break;/* leftmost aligned pair in inverted matchgroup is confirmed at tL2 as (I2,K2,J2)  */
    }

    for(int t2 = tR2; t2 >= tL2; t2--){
      int I2 = q->sites1[t2];
      int K2 = q->sitesK1[t2];
      int J2 = M+1 - q->sites2[t2];/* orient same as larger matchgroup */

      /* check if there is an alignment from I2-K2 .. I2 to query via larger matchgroup that is outside the range qrystartidx1 .. qryendidx1 */
      int J1 = -1;
      for(int t = U; --t >= 0; t--){
	int I = p->sites1[t];
	int K = p->sitesK1[t];
	int J = p->sites2[t];
	if(max(I-K,I2-K2) <= min(I,I2)){/* found aligned pair */
	  J1 = J;
	  break;
	}
	if(I <= I2-K2)
	  break;
      }
      if(J1 > qryendidx1 && (matchcnt >= 2 || !backtrack)){/* delete aligned pair t2 since it it will reduce matchcnt by 1 and matchcnt2 by at least 1 */
	if(INV_VERB>=2){
	  printf("t2=%d:I2=%d,K2=%d,J2=%d,tL2=%d,tR2=%d:J1=%d,qryendidx1= %d: deleting aligned pair t2 from right end of Inverted matchgroup\n",
		 t2,I2,K2,J2,tL2,tR2,J1,qryendidx1);
	  fflush(stdout);
	}
	if(DEBUG) assert(tR2 == t2);
	tR2 = t2-1;
	backtrack = 0;
        if(tR2 >= 0){// NEW11
	  refendidx = q->sites1[tR2];
	  qrystartidx1 = M+1 - q->sites2[tR2];
	} else {
          refendidx = q->sites1[0] - 1;
          qrystartidx1 = M+1 - (q->sites2[0] - 1);
        }
        if(DEBUG) assert(refendidx >= 0 && qrystartidx1 <= M+1);
	qrystartidx = orientation ? (M+1 - qrystartidx1) : qrystartidx1;

	matchcnt--;
	continue;/* try to delete additional aligned pairs t2-1 etc */
      }

      /* check if there is an alignment from ref to J2 via larger matchgroup, that is outside the range refstartidx .. refendidx */
      int I1 = -1, K1= -1;
      for(int t = 0; t < U; t++){
	int I = p->sites1[t];
	int K = p->sitesK1[t];
	int J = p->sites2[t];

	if(J == J2){/* found aligned pair */
	  I1 = I;
	  K1 = K;
	  break;
	}
	if(J > J2)
	  break;
      }
      if(I1 >= 0 && I1 < refstartidx && (matchcnt >= 2 || !backtrack)) { /* delete aligned pair t2 since it it will reduce matchcnt by 1 and matchcnt2 by at least 1 */
	if(INV_VERB>=2){
	  printf("t2=%d:I2=%d,K2=%d,J2=%d,tL2=%d,tR2=%d:I1=%d,refstartidx= %d: deleting aligned pair t2 from right end of Inverted matchgroup\n",
		 t2,I2,K2,J2,tL2,tR2,I1,refstartidx);
	  fflush(stdout);
	}

	if(DEBUG) assert(tR2 == t2);
	tR2 = t2 - 1;
	backtrack = 0;
        if(tR2 >= 0){// NEW11
	  refendidx = q->sites1[tR2];
	  qrystartidx1 = M+1 - q->sites2[tR2];
	} else {
          refendidx = q->sites1[0] - 1;
          qrystartidx1 = M+1 - (q->sites2[0] - 1);
        }
        if(DEBUG) assert(refendidx >= 0 && qrystartidx1 <= M+1);
	qrystartidx = orientation ? (M+1 - qrystartidx1) : qrystartidx1;

	matchcnt--;
	continue;/* try to delete additional aligned pairs t2 - 1 etc */
      }

      /* check if we can delete both t2 and tR2 reducing both matchcnt and matchcnt2 by 2 : will revert last such change to keep inversion region as large as plausible */
      if(I1 >= 0 && I1 - K1 <= refstartidx && I1 >= refstartidx && J1 == qryendidx1 && tR2 >= tL2+1){
	if(INV_VERB>=2){
	  printf("t2=%d,tL2=%d,tR2=%d,U2=%d:I2=%d,K2=%d,J2=%d:J1=%d,I1=%d,K1=%d,qryendidx1=%d,refstartidx=%d: deleting aligned pair tL2,tR2 from both ends of Inverted matchgroup\n",
		 t2,tL2,tR2,U2,I2,K2,J2,J1,I1,K1,qryendidx1,refstartidx);
	  fflush(stdout);
	}
	if(DEBUG) assert(tR2 == t2);
	backtrack = 1;
	tL2 = tL2 + 1;// NEW11
	tR2 = tR2 - 1;// NEW11
        if(DEBUG)assert(0 <= tR2 && tL2 < U2);
	refstartidx = q->sites1[tL2] - q->sitesK1[tL2];
	qryendidx1 = M+1 - q->sites2[tL2];
	qryendidx = orientation ? (M+1 - qryendidx1) : qryendidx1;
	refendidx = q->sites1[tR2];
	qrystartidx1 = M+1 - q->sites2[tR2];
	qrystartidx = orientation ? (M+1 - qrystartidx1) : qrystartidx1;

	matchcnt -= 2;
	continue;/* try to delete additional aligned pairs */
      }

      if(DEBUG) assert(tR2 == t2);
      // WAS11      tR2 = t2;
      break;/* rightmost aligned pair in inverted matchgroup is confirmed at tR2 as (I2,K2,J2) */
    }

    if(backtrack){
      if(DEBUG) assert(tL2 > 0 && tR2 < U2-1);
      tL2--;
      tR2++;

      refstartidx = q->sites1[tL2] - q->sitesK1[tL2];
      qryendidx1 = M+1 - q->sites2[tL2];
      qryendidx = orientation ? (M+1 - qryendidx1) : qryendidx1;
    
      refendidx = q->sites1[tR2];
      qrystartidx1 = M+1 - q->sites2[tR2];
      qrystartidx = orientation ? (M+1 - qrystartidx1) : qrystartidx1;

      matchcnt = tR2 - tL2 + 1;

      if(INV_VERB>=2){
	printf("Backtracked truncating inverted matchgroup: tL2=%d, tR2=%d : refidx = %d..%d, qryidx = %d .. %d, qry1 = %d .. %d, matchcnt = %d -> %d\n",
	       tL2,tR2,refstartidx,refendidx,qrystartidx,qryendidx,qrystartidx1,qryendidx1, U2, matchcnt);
	fflush(stdout);
      }

    } else if(DEBUG>=1+RELEASE/* HERE HERE >=2 */ && tL2 <= tR2){

      if(!(refstartidx == q->sites1[tL2] - q->sitesK1[tL2])){
	printf("refid=%lld,qryid=%lld,xmapid1=%d(conf=%0.2f),xmapid2=%d(conf=%0.2f),fwd1=%d:tL2=%d,tR2=%d,U2=%d:RefInv=%d..%d,QryInv=%d..%d(qry1=%d..%d),q->sites1[tL2]=%d,q->sitesK1[tL2]=%d\n",
	       usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,j,usexmapentries[j]->confidence,i,usexmapentries[i]->confidence,usexmapentries[j]->orientforward ? 1 : 0,
	       tL2,tR2,U2,refstartidx,refendidx,qrystartidx,qryendidx,qrystartidx1,qryendidx1,q->sites1[tL2],q->sitesK1[tL2]);
	fflush(stdout);
	assert(refstartidx == q->sites1[tL2] - q->sitesK1[tL2]);
      }
      assert(qryendidx1 == M+1 - q->sites2[tL2]);
      assert(refendidx == q->sites1[tR2]);
      if(!(qrystartidx1 == M+1 - q->sites2[tR2])){
	printf("refid=%lld,qryid=%lld,xmapid1=%d(conf=%0.2f),xmapid2=%d(conf=%0.2f),fwd1=%d:tL2=%d,tR2=%d,U2=%d:RefInv=%d..%d,QryInv=%d..%d(qrystartidx1=%d,qryendidx1=%d),q->sites2[0]=%d,q->sites2[tR2]=%d\n",
	       usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,j,usexmapentries[j]->confidence,i,usexmapentries[i]->confidence,usexmapentries[j]->orientforward ? 1 : 0,
	       tL2,tR2,U2,refstartidx,refendidx,qrystartidx,qryendidx,qrystartidx1,qryendidx1,q->sites2[0],q->sites2[tR2]);
	fflush(stdout);
	assert(qrystartidx1 == M+1 - q->sites2[tR2]);
      }
	
      assert(matchcnt == tR2 - tL2 + 1);
    }

    if(DEBUG) assert(tL2 >= 0 && tR2 < U2);

    if(tL2 > tR2)
      return 0;/* inverted matchgroup must not overlap larger matchgroup on either query or reference, OR all alignments overlap alignments of larger matchgroup */
  }

  if(qrystartidx == qryendidx) 
    return 0;// NEW83

  /* find the two aligned pairs in the Large matchgroup that bound the smaller matchgroup in both the Query and Reference (or just Query if InvDup is true) */
  tL = tR = -1;/* alignment index values that most closely bound the smaller (inverted matchgroup) on both Query and Reference (or just Query if InvDup) */
  refstartidx2 = refendidx2 = qrystartidx2 = qryendidx2 = -1;/* reference and query label ranges for matchgroup j corresponding to tL and tR */
  matchcnt1 = 0;/* If InvDup : count of aligned pairs in Large matchgroup (outside bounding range, including bounding pairs) that overlap the smaller matchgroup Reference Range : this should be a high fraction of matchcnt */
  matchcnt2 = 0;/* Count of aligned pairs in Large matchgroup (within bounding range, including the bounding pairs) that overlap the smaller matchgroup in either Query or Reference ranges 
		   (or just Query ranges if InvDup) : this should be small fraction of matchcnt */

  for(int t = 0; t < U; t++){
    int I = p->sites1[t];
    int K = p->sitesK1[t];
    int J = p->sites2[t];

    // update tL, tR, refstartidx2,refendidx2,qrystartidx2,qryendidx2
    if((InvDup || I <= refstartidx) && J <= qrystartidx1){/* update left bound */
      tL = t;
      if(DEBUG) assert(I >= refstartidx2);
      refstartidx2 = I;
      if(DEBUG) assert(J >= qrystartidx2);
      qrystartidx2 = J;
      if((!InvDup && I >= refstartidx) || J >= qrystartidx1)
	matchcnt2++;

      if(InvDup && refstartidx <= I && I <= refendidx)
	matchcnt1++;
    } else if((InvDup || I-K >= refendidx) && J >= qryendidx1){/* at or beyond right bound */
      if(refendidx2 < 0){/* found right bound */
	tR = t;
	refendidx2 = I-K;
	qryendidx2 = J;
      }
      if((!InvDup && I-K <= refendidx) || J <= qryendidx1)
	matchcnt2++;
      if(InvDup && refstartidx <= I && I <= refendidx)
	matchcnt1++;
    } else /* between bounds */
      //    if((InvDup || (I >= refstartidx && I-K <= refendidx)) && J >= qrystartidx1 && J <= qryendidx1)
      matchcnt2++;

    if(INV_VERB>=3)
       printf("\t t=%d/%d:I=%d,K=%d,J=%d : ref2=%d..%d,qry2=%d..%d,matchcnt=%d, matchcnt2=%d\n",t, U, I,K,J, refstartidx2,refendidx2,qrystartidx2,qryendidx2,matchcnt, matchcnt2);
  }
  
  if(DEBUG && refstartidx2 >= 1 && refendidx2 >= 1){
    if(!(refstartidx2 <= refendidx2 && qrystartidx2 <= qryendidx2)){
      printf("refstartidx2=%d,refendidx2=%d,qrystartidx2=%d,qryendidx2=%d\n",
	     refstartidx2,refendidx2,qrystartidx2,qryendidx2);
      fflush(stdout);
      assert(refstartidx2 <= refendidx2);
      assert(qrystartidx2 <= qryendidx2);
    }
  }

  return 1;
}

// New inversions will be appended to SmallInversions[0..NumSmallInversions-1] : return 1 if anything was appended
static int checkSmallInversion(CSmallInversion *&SmallInversions, 
			       int &NumSmallInversions, 
			       int &MaxSmallInversions,
			       xmapEntry **usexmapentries, 
			       int j, // usexmapentries[j] is the larger (higher confidence) matchgroup 
			       int i // usexmapentries[i] is the smaller (lower confidence) matchgroup of opposite orientation
			       )
{
  int refstartidx, refendidx, qrystartidx, qryendidx, qrystartidx1, qryendidx1, matchcnt, matchcnt1, matchcnt2;
  int refstartidx2 = -1,refendidx2 = -1,qrystartidx2 = -1,qryendidx2 = -1;/* reference and query label ranges for matchgroup j that most tightly bound matchgroup i (after Palindromic end trimming) */
  int tL = -1, tR = -1;/* alignment index values that most closely bound the smaller (inverted matchgroup) on both Query and Reference */
  int M,N,U;
  FLOAT *X,*Y;
  Calign *p, *q;

  if(!SmallInversionOverlap(usexmapentries, j, i, p, q, M, N, X, Y, U,
			    refstartidx, refendidx, qrystartidx, qryendidx, qrystartidx1,qryendidx1, matchcnt, matchcnt1, matchcnt2,
			    refstartidx2, refendidx2, qrystartidx2, qryendidx2,tL,tR, false))
    return 0;// corner case : inversion ruled out since all qry aligned labels were trimmed away as part of an inverted palindrome


  if(INV_VERB>=2){
    printf("refstartidx2=%d,refendidx2=%d, qrystartidx2=%d,qryendidx2=%d,M=%d: refstartidx=%d, refendidx=%d, qrystartidx1=%d,qryendidx1=%d:matchcnt2=%d, matchcnt=%d, CutFlipMaxOverlap= %0.4f, CutFlipMaxMiss= %0.4f\n",
	   refstartidx2,refendidx2,qrystartidx2,qryendidx2,M,refstartidx,refendidx,qrystartidx1,qryendidx1,matchcnt2,matchcnt,CutFlipMaxOverlap, CutFlipMaxMiss);
    printf("\t matchcnt2= %d, matchcnt=%d\n",matchcnt2,matchcnt);
    fflush(stdout);
  }

  if(refstartidx2 >= 1 && refendidx2 >= 1 && matchcnt2 <= floor(matchcnt * CutFlipMaxOverlap + 0.5) /* overlap in labels / matchcnt <= CutFlipMaxOverlap (default = 50%) */
     && (CutFlipMaxMiss <= 0.0 || max(refstartidx - refstartidx2 + qrystartidx1 - qrystartidx2, refendidx2 - refendidx + qryendidx2 - qryendidx1) - 2 <= floor(matchcnt * CutFlipMaxMiss + 0.5))){
    if(NumSmallInversions >= MaxSmallInversions){
      MaxSmallInversions *= 2;
      CSmallInversion *newSmallInversions = new CSmallInversion[MaxSmallInversions];
      memcpy(newSmallInversions,SmallInversions,NumSmallInversions * sizeof(CSmallInversion));
      delete [] SmallInversions;
      SmallInversions = newSmallInversions;
    }

    CSmallInversion *pInv = &SmallInversions[NumSmallInversions++];
    pInv->refcontigid = usexmapentries[j]->refcontigid;
    pInv->qrycontigid = usexmapentries[j]->qrycontigid;
    pInv->orientforward = usexmapentries[j]->orientforward;
	      
    pInv->Large = usexmapentries[j];
    pInv->Small = usexmapentries[i];

    if(DEBUG) assert(usexmapentries[j]->QryMap && usexmapentries[j]->QryMap == usexmapentries[i]->QryMap);
    if(DEBUG) assert(usexmapentries[j]->RefMap && usexmapentries[j]->RefMap == usexmapentries[i]->RefMap);
    pInv->Y = Y;
    pInv->N = N;
    pInv->X = X;
    pInv->M = M;

    pInv->RefOutlierStart = refstartidx2;
    pInv->RefOutlierEnd = refendidx2;
    pInv->QryOutlierStart = pInv->orientforward ? qrystartidx2 : M+1 - qrystartidx2 /* WAS qryendidx2 */;
    pInv->QryOutlierEnd = pInv->orientforward ? qryendidx2 : M+1 - qryendidx2 /* WAS qrystartidx2 */;

    pInv->RefInversionStart = refstartidx;
    pInv->RefInversionEnd = refendidx;
    pInv->QryInversionStart = qrystartidx; // aligned to refendidx // WAS usexmapentries[i]->qrystartidx;
    pInv->QryInversionEnd = qryendidx; // aligned to refstartidx // WAS usexmapentries[i]->qryendidx;

    /* left confidence is based on alignment p->sites1[t = 0..tL] etc */
    extern double SmRA(int J, int I, int K, FLOAT *Y);
    double Leftscore = SmRA(p->sites2[0], p->sites1[0], p->sitesK1[0], Y);
    for(int t = 1; t <= tL; t++)
      Leftscore += p->iscore[t];
    pInv->LeftConf = Leftscore / log(10.0);

    double Rightscore = SmRA(p->sites2[U-1], p->sites1[U-1], p->sitesK1[U-1], Y);
    for(int t = tR+1; t < U; t++)
      Rightscore += p->iscore[t];
    pInv->RightConf = Rightscore / log(10.0);

    if(VERB){
      printf("Found Small Inversion: refid=%lld,qryid=%lld, xmapid1= %d, xmapid2= %d, forwd1= %d:RefOutlier= %d .. %d(%0.3f..%0.3f), QryOutlier= %d .. %d(%0.3f..%0.3f), RefInversion= %d .. %d(%0.3f..%0.3f), QryInversion= %d .. %d(%0.3f..%0.3f),cnt=%d,%d, Conf:Left= %0.2f, Right= %0.2f, Inv= %0.2f\n",
	     pInv->refcontigid,pInv->qrycontigid,usexmapentries[j]->xmapid,usexmapentries[i]->xmapid,
	     pInv->orientforward ? 1 : 0, pInv->RefOutlierStart, pInv->RefOutlierEnd, Y[pInv->RefOutlierStart],Y[pInv->RefOutlierEnd],
	     pInv->QryOutlierStart, pInv->QryOutlierEnd,X[pInv->QryOutlierStart],X[pInv->QryOutlierEnd],
	     pInv->RefInversionStart, pInv->RefInversionEnd, Y[pInv->RefInversionStart], Y[pInv->RefInversionEnd], 
	     pInv->QryInversionStart, pInv->QryInversionEnd, X[pInv->QryInversionStart], X[pInv->QryInversionEnd],matchcnt,matchcnt2, pInv->LeftConf, pInv->RightConf, pInv->Small->confidence);
      fflush(stdout);
    }

    return 1;
  }

  return 0;
}

#if 0 // only used for debugging
/* use binary search to locate index i such that xmapentries[i]->xmapid == xmapid. 
   Returns index i if found.
   Exit with error if not found
*/
static int findXmapid(int xmapid, xmapEntry **xmapentries, int numxmapentries)
{
  int low = 0, high = numxmapentries-1;
  while(high > low + 1){
    int mid = (high + low) / 2;
    int val = xmapentries[mid]->xmapid;
    if(val > xmapid)
      high = mid;
    else if(val < xmapid)
      low = mid;
    else
      return mid;
  }
  if(xmapid == xmapentries[low]->xmapid)
    return low;
  if(xmapid == xmapentries[high]->xmapid)
    return high;
  
  printf("findXmapid: Cannot find xmapid=%d in xmapentries[0..%d]\n",xmapid,numxmapentries-1);
  fflush(stdout);exit(1);
}
#endif

static int findXmapEntry(xmapEntry *p, xmapEntry **xmapentries, int numxmapentries)
{
  for(int i = 0; i < numxmapentries; i++)
    if(p == xmapentries[i]){
      if(DEBUG>=2){/* make sure there is not a 2nd instance */
	for(int j = i+1; j < numxmapentries; j++){
	  if(p == xmapentries[j]){
	    printf("WARNING: 2 copies in xmapentries[0..%d] of xmapEntry=%p at i=%d,j=%d:xmapid=%d,refid=%lld,qryid=%lld,fwd=%d:qrystartpos=%0.1f,qryendpos=%0.1f,refstartpos=%0.1f,refendpos=%0.1f,conf=%0.2f\n",
		   numxmapentries-1,p,i,j,p->xmapid,p->refcontigid,p->qrycontigid,p->orientforward ? 1 : 0, p->qrystartpos,p->qryendpos,p->refstartpos, p->refendpos, p->confidence);
	    fflush(stdout);
	    assert(xmapentries[i] != xmapentries[j]);
	  }
	}
      }
      return i;
    }
  return -1;
}

/* sort alignments in order of id1, sites1[0], sites1[N+1], id2, sites2[0],sites2[M+1] (null alignments last) */
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

/* sort alignments in order of mapid1, sites1[0], sites1[N+1], id2, sites2[0],sites2[M+1] (null alignments last) */
static int RefSiteInc(Calign **p1, Calign **p2)
{
  Calign *align1 = p1[0];
  Calign *align2 = p2[0];
  if(!align1) return 1;
  if(!align2) return -1;

  Cmap *Ymap1 = YYmap[align1->mapid1];
  Cmap *Ymap2 = YYmap[align2->mapid1];
  int left1 = align1->sites1 ? (align1->sites1[0] + (Ymap1->origmap ? max(1,Ymap1->left[0]) -1 : 0)) : 0 ;
  int left2 = align2->sites1 ? (align2->sites1[0] + (Ymap2->origmap ? max(1,Ymap2->left[0]) -1 : 0)) : 0 ;
  int ret = (align1->mapid1 > align2->mapid1) ? 1 : (align1->mapid1 < align2->mapid1) ? -1 :
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

/* qsort() intcmp function to sort sv_array[] in ascending order of type_enum and SVsize (computed locally) */
static inline int SVenumSVsizeInc(structuralVariation *x, structuralVariation *y)
{
  if(x->type_enum < y->type_enum)
    return -1;
  if(x->type_enum > y->type_enum)
    return 1;

  double xSVsize = -1.0;
  if(x->type_enum == insertion)
    xSVsize = fabs(x->querystop - x->querystart) - (x->refstop - x->refstart);
  else if(x->type_enum == deletion || x->type_enum == compression)
    xSVsize = (x->refstop - x->refstart) - fabs(x->querystop - x->querystart);
  else if(x->type_enum == duplication_type || x->type_enum == duplicationinv_type)
    xSVsize = x->refstop - x->refstart;
  
  double ySVsize = -1.0;
  if(y->type_enum == insertion)
    ySVsize = fabs(y->querystop - y->querystart) - (y->refstop - y->refstart);
  else if(y->type_enum == deletion || y->type_enum == compression)
    ySVsize = (y->refstop - y->refstart) - fabs(y->querystop - y->querystart);
  else if(y->type_enum == duplication_type || y->type_enum == duplicationinv_type)
    ySVsize = y->refstop - y->refstart;
  
  if(xSVsize < ySVsize)
    return -1;
  if(xSVsize > ySVsize)
    return 1;

  // for ties maintain original order
  if(x < y)
    return -1;
  if(x > y)
    return 1;
  return 0;
}

/* qsort() intcmp function to sort usexmapentries[] array in ascending order of qrycontigid, refcontigid, refstartpos,refendpos (since indel duplicate check relies on order by refcontigid,refstartpos,refendpos): 
   NOTE : For ties preserve original order
 */
static inline int xmapEntryQryidInc(xmapEntry **p1, xmapEntry **p2)
{
  long long qrycontigid1 = p1[0]->qrycontigid;
  long long qrycontigid2 = p2[0]->qrycontigid;
  return (qrycontigid1 > qrycontigid2) ? 1 : (qrycontigid1 < qrycontigid2) ? -1 : 
    (p1[0]->refcontigid > p2[0]->refcontigid) ? 1 : (p1[0]->refcontigid < p2[0]->refcontigid) ? -1 :
    (p1[0]->refstartpos > p2[0]->refstartpos) ? 1 : (p1[0]->refstartpos < p2[0]->refstartpos) ? -1 :
    (p1[0]->refendpos > p2[0]->refendpos) ? 1 : (p1[0]->refendpos < p2[0]->refendpos) ? -1 :
    (!p1[0]->orientforward && p2[0]->orientforward) ? 1 : (p1[0]->orientforward && !p2[0]->orientforward) ? -1 :
    (p1[0]->orientforward ? (p1[0]->qrystartpos > p2[0]->qrystartpos) : (p1[0]->qrystartpos < p2[0]->qrystartpos)) ? 1 :
    (p1[0]->orientforward ? (p1[0]->qrystartpos < p2[0]->qrystartpos) : (p1[0]->qrystartpos > p2[0]->qrystartpos)) ? -1 :
    (p1[0]->orientforward ? (p1[0]->qryendpos > p2[0]->qryendpos) : (p1[0]->qryendpos < p2[0]->qryendpos)) ? 1 :
    (p1[0]->orientforward ? (p1[0]->qryendpos < p2[0]->qryendpos) : (p1[0]->qryendpos > p2[0]->qryendpos)) ? -1 :
    (p1 > p2) ? 1 : (p1 < p2) ? -1 : 0;
}

#if 0
/* qsort() intcmp function to sort xmapentries[] array in ascending order of xmapid (there should be no ties) */
static inline int xmapEntryXmapidInc(xmapEntry **p1, xmapEntry **p2)
{
  return (p1[0]->xmapid > p2[0]->xmapid) ? 1 : (p1[0]->xmapid < p2[0]->xmapid) ? -1 : 0;
}
#endif

/* compute duplication query coords.
   Also compute a safe lower bound on the spacing between the duplicate regions on the query (in kb) */
void getDuplicationCoords(xmapEntry* xe1, xmapEntry* xe2, int refstartidx, int refstopidx,
			  double& qrystart, int& qrystartidx, double& qrystop, int& qrystopidx, 
			  double &DupSpacing, bool verbose)
{
  Calign *p1 = xe1->align;
  Calign *p2 = xe2->align;
  if(DEBUG) assert(p1->mapid1 == p2->mapid1 && p1->mapid2 == p2->mapid2);
  Cmap *Ymap = YYmap[p1->mapid1];
  Cmap *Xmap = XXmap[p1->mapid2];
  FLOAT *Y = Ymap->site[0];
  FLOAT *X = Xmap->site[0];
  //  int N = Ymap->numsite[0];
  //  int M = Xmap->numsite[0];


  double startloc1Q=0, stoploc1Q=0, startloc2Q=0, stoploc2Q=0;
  int startidx1=0, stopidx1=0, startidx2=0, stopidx2=0;
  if(DEBUG && !(refstartidx < refstopidx)){
    printf("  dupCoord: refid=%lld,qryid=%lld: refstartidx= %d(%0.3f), refstopidx= %d (%0.3f)\n", YYmap[p1->mapid1]->id,XXmap[p1->mapid2]->id, refstartidx, Y[refstartidx], refstopidx, Y[refstopidx]);
    fflush(stdout);
    assert(refstartidx < refstopidx);
  }
  getQryCoord(p1, refstartidx, refstopidx, startloc1Q, startidx1, stoploc1Q, stopidx1);
  getQryCoord(p2, refstartidx, refstopidx, startloc2Q, startidx2, stoploc2Q, stopidx2);

  if(verbose){
    printf("  dupCoord: refid=%lld,qryid=%lld: refidx= %i,%i (%0.3f,%0.3f) qryloc1= %.1f, %.1f qidx1= %i %i (%.3f,%.3f) qryloc2= %.1f, %.1f qidx2= %i %i (%0.3f,%0.3f)\n", YYmap[p1->mapid1]->id,XXmap[p1->mapid2]->id,
	   refstartidx, refstopidx, Y[refstartidx], Y[refstopidx], startloc1Q, stoploc1Q, startidx1, stopidx1, X[startidx1],X[stopidx1], startloc2Q, stoploc2Q, startidx2, stopidx2, X[startidx2], X[stopidx2]);
    fflush(stdout);
  }

  int min1 = min(startidx1,stopidx1);
  int max1 = max(startidx1,stopidx1);
  int min2 = min(startidx2,stopidx2);
  int max2 = max(startidx2,stopidx2);

  qrystart = min(min(startloc1Q, stoploc1Q), min(startloc2Q, stoploc2Q));
  qrystop  = max(max(startloc1Q, stoploc1Q), max(startloc2Q, stoploc2Q));
  qrystartidx = min(min1,min2);
  qrystopidx  = max(max1,max2);

  if(max1 < min2){/* qry1 is left of qry2 with a gap between them : subtract first label interval from each end of gap, so DupSpacing is not overestimated */
    DupSpacing = max(0.0, X[min2-1] - X[max1+1]);
  } else if(max2 < min1){/* qry2 is left of qyr1 */
    DupSpacing = max(0.0, X[min1-1] - X[max2+1]);
  } else /* qry regions have no gap or may even be overlapped */
    DupSpacing = 0.0;
}


//figure out query coordinates for duplications:
// though there are four coordinates, we currently output only two, so for now just store those two in
// querystart, querystop. They are min(first copy coords), max(second copy coords) where first/second
// is order on query.
void getDuplicationCoords_orig(xmapEntry* xe1, xmapEntry* xe2, double refstart, double refstop, double minrefalignlen,
			       double refoverlapsize, bool qryoverlap, double minq1, double minq2, double maxq1, double maxq2,
			       double qrysize, double& qrystart, double& qrystop) {
  bool verbose = true;
  //hack some locals so that I don't have to change below
  xmapEntry **usexmapentries = new xmapEntry*[2]; //copy elements of xmapentries
  int i=0, j=1;
  usexmapentries[i] = xe1;
  usexmapentries[j] = xe2;

  if(verbose)
    printf("  duplication: "); //  dupli-fixed: ");
  //query coords are a bit tricky...
  if( minrefalignlen == refoverlapsize && //complete ref overlap
      usexmapentries[i]->refstartpos != usexmapentries[j]->refstartpos &&
      usexmapentries[i]->refendpos   != usexmapentries[j]->refendpos ) {
    if( fabs(usexmapentries[i]->refstartpos - usexmapentries[i]->refendpos) <=
	fabs(usexmapentries[j]->refstartpos - usexmapentries[j]->refendpos) ) {
      //turns out that this is impossible, especially with the requirement above that the start coordinates are not the same
      printf("\nUnhandled duplication: i in j: ref i: %.1f %.1f, ref j: %.1f  %.1f\n", usexmapentries[i]->refstartpos, usexmapentries[i]->refendpos, usexmapentries[j]->refstartpos, usexmapentries[j]->refendpos); fflush(stdout);
      fflush(stdout);
      assert( !(fabs(usexmapentries[i]->refstartpos - usexmapentries[i]->refendpos) <=
		fabs(usexmapentries[j]->refstartpos - usexmapentries[j]->refendpos)) );
    } else if( qryoverlap ) {
      //find qrystart/end on i which is inside j (ql is opposite)
      bool lr = ((usexmapentries[i]->qrystartpos > usexmapentries[j]->qrystartpos &&
		  usexmapentries[i]->qrystartpos < usexmapentries[j]->qryendpos) ||
		 (usexmapentries[i]->qrystartpos > usexmapentries[j]->qryendpos &&
		  usexmapentries[i]->qrystartpos < usexmapentries[j]->qrystartpos));
      double ql = ( lr ? usexmapentries[i]->qryendpos : usexmapentries[i]->qrystartpos );
      double r1 = ( lr ? usexmapentries[i]->refendpos : usexmapentries[i]->refstartpos );
      double r2 = ( fabs(r1-usexmapentries[j]->refstartpos) < fabs(r1-usexmapentries[j]->refendpos) ?
		    usexmapentries[j]->refstartpos : usexmapentries[j]->refendpos );
      if(verbose)
	printf("ji ovl, ql %.1f r1 %.1f r2 %.1f ", ql, r1, r2);
      ql = ql + ( lr == usexmapentries[i]->orientforward ? -1. : 1. ) * fabs(r1-r2);
      qrystart = ( minq1 <= minq2 ? ql : minq2 );
      qrystop  = ( minq1 <= minq2 ? maxq2 : ql );
    } else { //j is inside i (should be same as above after changing i <-> j and >/< in lr)
      bool lr = (usexmapentries[i]->qrystartpos - minq2 > usexmapentries[i]->qryendpos - minq2);
      double ql = ( lr ? usexmapentries[i]->qrystartpos : usexmapentries[i]->qryendpos );
      double r1 = ( lr ? usexmapentries[i]->refstartpos : usexmapentries[i]->refendpos );
      double r2 = ( lr ? usexmapentries[j]->refstartpos : usexmapentries[j]->refendpos );
      if(verbose)
	printf("j in i, ql %.1f r1 %.1f r2 %.1f ", ql, r1, r2);
      float s = ( lr == usexmapentries[i]->orientforward ? 1. : -1. ); //sign in next expression
      ql = ql + s*fabs(r1-r2) + s*fabs(usexmapentries[j]->qrystartpos - usexmapentries[j]->qryendpos);
      qrystart = ( minq1 <= minq2 ? ql : minq2 );
      qrystop  = ( minq1 <= minq2 ? maxq2 : ql );
    }
  } else { //partial ref overlap of both MGs
    if( qrysize < 0 ) { // qrysize > 0 means the MGs are in order
      if(verbose)
	printf("p out ");
      //for this case, it's always the min/max of qry coords
      qrystart = min(minq1, minq2);
      qrystop  = max(maxq1, maxq2);
    } else {
      if(verbose)
	printf("p in  ");
      //simple formula for remaining cases: add/subtract refoverlapsize from qrystart/stop based on direction of overlap:
      // because refstart i is always <= j, it is always refend of i and refstart of j which overlap because complete 
      // overlap is handled above (right?); but add/subtract from that depends on orientation
      double q1 = usexmapentries[i]->qryendpos   + (usexmapentries[i]->orientforward ? -1 : 1)*refoverlapsize*1e3;
      double q2 = usexmapentries[j]->qrystartpos + (usexmapentries[j]->orientforward ? 1 : -1)*refoverlapsize*1e3;
      qrystart = min(q1,q2);
      qrystop  = max(q1,q2);
    }
  }
  //hack: if qrystart is < 20, set it to 20: this could happen, eg, due to stretch or indels. It happened in two duplications in the sample I'm testing on, positions were -50 and - 160 (bp).
  qrystart = max(20., qrystart);
  if(verbose)
    printf(": refstart/stop: %.1f %.1f, qrystart/stop: %.1f %.1f\n", refstart, refstop, qrystart, qrystop);

  delete[] usexmapentries;
}

class CmatchgroupPair {/* pair of matchgroups that are overlapped */
public:
  int i;/* 1st matchgroup is usexmapentries[i] */
  int j;/* 2nd matchgroup is usexmapentries[j] */
  int iLabels;/* usexmapentries[i]->align->numpairs */
  int jLabels;/* usexmapentries[j]->align->numpairs */
};

static inline int MGlabelsDec(CmatchgroupPair *p1, CmatchgroupPair *p2)
{
  int min1 = min(p1->iLabels,p1->jLabels);
  int max1 = max(p1->iLabels,p1->jLabels);
  int min2 = min(p2->iLabels,p2->jLabels);
  int max2 = max(p2->iLabels,p2->jLabels);
  
  return (min1 < min2) ? 1 : (min1 > min2) ? -1 :
    (max1 < max2) ? 1 : (max1 > max2) ? -1 : 
    (p1 > p2) ? 1 : (p1 < p2) ? -1 : 0;
}

/* All matchgroup pairs in MGlist[0..MGlistlen-1] have same qryid and refid AND query overlap (in kb) that exceeds both CutFlipOverlapFilter AND OverlapFilterQryRatio.

   1. Sort matchgroup pairs in descending order of min(iLabels,jLabels), max(iLabels,jLabels)
   2. Delete smaller of each matchgroup pair unless it needs to be retained due to possible (non-inverted) Small Duplication or small (CutFlip) Inversion (or inverted Duplication).
      (See below for detailed conditions based  on -CutFlipFilter, -SmallDupFilter, -OverlapFilter and -InversionOverlapFilter)

*/
static void FilterMG(CmatchgroupPair *MGlist, 
		     int MGlistlen,
		     xmapEntry **usexmapentries,
		     int smapsize)
{
  if(DEBUG) assert(MGlistlen > 0);
  long long qid = usexmapentries[MGlist[0].i]->qrycontigid;
  long long rid = usexmapentries[MGlist[0].i]->refcontigid;

  qsort(MGlist,MGlistlen, sizeof(CmatchgroupPair), (intcmp *)MGlabelsDec);

  if(VERB>=2){
    printf("FilterMG : checking overlap of %d matchgroup pairs with rid=%lld,qid=%lld:\n",MGlistlen, rid,qid);
    fflush(stdout);
  }

  for(int t = 0; t < MGlistlen; t++){
    CmatchgroupPair *pMG = &MGlist[t];
    int i = pMG->i;
    int j = pMG->j;
    if(DEBUG && ! (0 <= i && i < smapsize)){
      printf("MGlist[%d]:i=%d,j=%d,smapsize=%d\n",t,i,j,smapsize);
      fflush(stdout);
      assert(0 <= i && i < smapsize);
    }
    if(DEBUG) assert(0 <= j && j < smapsize);
    if(usexmapentries[i]->confidence > usexmapentries[j]->confidence){/* swap i,j, so i is always the smaller matchgroup */
      int tmp = i;
      i = j;
      j = tmp;
    }
    if(DEBUG>=2) assert(usexmapentries[i]->qrycontigid == usexmapentries[j]->qrycontigid);
    if(DEBUG>=2) assert(usexmapentries[i]->refcontigid == usexmapentries[j]->refcontigid);
    if(DEBUG) assert(usexmapentries[i]->qrycontigid == qid);
    if(DEBUG) assert(usexmapentries[j]->qrycontigid == qid);
    if(DEBUG) assert(usexmapentries[i]->refcontigid == rid);
    if(DEBUG) assert(usexmapentries[j]->refcontigid == rid);
    if(DEBUG>=2){
      xmapEntry *p = usexmapentries[i];
      if(!XXmap[p->align->mapid2]->origmap) assert(p->querycontignsite == XXmap[p->align->mapid2]->numsite[0]);

      p = usexmapentries[j];
      if(!XXmap[p->align->mapid2]->origmap) assert(p->querycontignsite == XXmap[p->align->mapid2]->numsite[0]);
    }
    
#if 0
    if(usexmapentries[i]->sv_overlap)
      continue;
    if(usexmapentries[j]->sv_overlap)
      continue;
#endif

    double mini = min(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
    double minj = min(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
    double maxi = max(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
    double maxj = max(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
    double overlap = min(maxi,maxj) - max(mini,minj);
    double minlen = min(maxi-mini, maxj-minj);
    
    Calign *p = usexmapentries[i]->align;
    Calign *q = usexmapentries[j]->align;

    if(VERB>=2 && overlap > 0.0){
      printf("FilterMG:t=%d:xmapEntry i=%d (xid=%d, rid=%lld,qid=%lld,fwd=%d,ov=%d,conf=%0.2f,np=%d,qry=%0.3f..%0.3f), j=%d (xid=%d,rid=%lld,qid=%lld,fwd=%d,ov=%d,conf=%0.2f,np=%d,qry=%0.3f..%0.3f) : overlap= %0.3f, minlen= %0.3f\n",
             t,i, usexmapentries[i]->xmapid,usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid,usexmapentries[i]->orientforward ? 1 : 0, 
	     usexmapentries[i]->sv_overlap ? 1 : 0, usexmapentries[i]->confidence, usexmapentries[i]->align->numpairs, mini*1e-3, maxi * 1e-3,
	     j, usexmapentries[j]->xmapid,usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,usexmapentries[j]->orientforward ? 1 : 0, 
	     usexmapentries[j]->sv_overlap ? 1 : 0, usexmapentries[j]->confidence, usexmapentries[j]->align->numpairs, minj*1e-3, maxj * 1e-3, overlap * 1e-3, minlen * 1e-3);
      fflush(stdout);
    }

#if 1
    if(usexmapentries[i]->sv_overlap)
      continue;
    if(usexmapentries[j]->sv_overlap)
      continue;
#endif

    /* For same orientation matchgroups, protect small matchgroup if query aligned label overlap fraction does NOT exceed max(max(0.5,SmallDupMaxOverlap),OverlapFilterQryLabRatio) */
    if(usexmapentries[i]->orientforward == usexmapentries[j]->orientforward &&
       SmallDupMaxOverlap > 0.0 && SmallDupMinOverlap > 0.0){
      if(DEBUG && !(p->orientation == q->orientation)){
	printf("xmapEntry i=%d (xid=%d, rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d,qry=%0.3f..%0.3f), j=%d (xid=%d,rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d,qry=%0.3f..%0.3f) : overlap= %0.3f, minlen= %0.3f\n",
	       i, usexmapentries[i]->xmapid,usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid,usexmapentries[i]->orientforward ? 1 : 0, 
	       usexmapentries[i]->confidence, usexmapentries[i]->align->numpairs, mini*1e-3, maxi * 1e-3,
	       j, usexmapentries[j]->xmapid,usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,usexmapentries[j]->orientforward ? 1 : 0, 
	       usexmapentries[j]->confidence, usexmapentries[j]->align->numpairs, minj*1e-3, maxj * 1e-3, overlap * 1e-3, minlen * 1e-3);
	printf("\n usexmapentries[i]->align->orientation=%d,usexmapentries[j]->align->orientation=%d\n",p->orientation,q->orientation);
	fflush(stdout);
	assert(p->orientation == q->orientation);
      }

      int U1,U2;
      int cnt = QryLabelOverlap(p,q, U1, U2);
      
      if(VERB>=2){
	printf(" i=%i,xid=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d) overlaps in len by %.4f with xid=%i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d), and in aligned query labels by %0.4f (%d/%d)\n", 
	       i,usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,usexmapentries[i]->align->numpairs,
	       overlap/minlen, usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
	       usexmapentries[j]->align->numpairs, (double)cnt / min(U1,U2), cnt, min(U1,U2));
	fflush(stdout);
      }

      if(cnt <= max(max(0.5, SmallDupMaxOverlap),OverlapFilterQryLabRatio) * min(U1,U2))
	continue;
    }

    /* For opposite orientation matchgroups, protect smaller matchgroup unless one of the following conditions (A OR B) apply :
       A. Smaller matchgroups is NOT fully overlapped on reference (May still be CutFlip candidate with a 3rd matchgroup) AND all of the following are true :
	   1. Larger matchgroup confidence exceeds smaller one by InvOverlapFilterConfDelta.
           2. Query overlap (fraction of kb length) exceeds InvOverlapFilterQryRatio.
           3. Query Overlap (fraction of aligned labels) exceeds max(0.5, CutFlipMaxOverlap,OverlapFilterQryLabRatio).
       B. Smaller matchgroup is fully overlapped on reference (CutFlip candidate) AND all of the following are true:
	   1. Larger matchgroup confidence exceeds smaller one by CutFlipFilterConf.
           2. (Temporarily) Trim palindromic ends of smaller matchgroup overlapping larger matchgroup, then check if Query Overlap (fraction of aligned labels) exceeds 
	      max(0.5, CutFlipMaxOverlap,OverlapFilterQryLabRatio).

       Note that typically the conditions for A are less stringent (less likely to protect smaller matchgroup) than the conditions for B. After checking for Small SegDups and CutFlip inversions, condition A
             is rechecked for all opposite orientation matchgroups (and the protected matchgroups are marked for use in inversion calls only)
    */
    if(usexmapentries[i]->orientforward != usexmapentries[j]->orientforward && CutFlipMaxOverlap > 0.0){
      if(DEBUG && !(p->orientation != q->orientation)){
	printf("xmapEntry i=%d (xid=%d, rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d,qry=%0.3f..%0.3f), j=%d (xid=%d,rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d,qry=%0.3f..%0.3f) : overlap= %0.3f, minlen= %0.3f\n",
	       i, usexmapentries[i]->xmapid,usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid,usexmapentries[i]->orientforward ? 1 : 0, 
	       usexmapentries[i]->confidence, usexmapentries[i]->align->numpairs, mini*1e-3, maxi * 1e-3,
	       j, usexmapentries[j]->xmapid,usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,usexmapentries[j]->orientforward ? 1 : 0, 
	       usexmapentries[j]->confidence, usexmapentries[j]->align->numpairs, minj*1e-3, maxj * 1e-3, overlap * 1e-3, minlen * 1e-3);
	printf("\n usexmapentries[i]->align->orientation=%d,usexmapentries[j]->align->orientation=%d\n",p->orientation,q->orientation);
	fflush(stdout);
	assert(p->orientation != q->orientation);
      }

      double refoverlap = min(usexmapentries[i]->refendpos, usexmapentries[j]->refendpos) - max(usexmapentries[i]->refstartpos, usexmapentries[j]->refstartpos);
      double refminlen = min(usexmapentries[i]->refendpos - usexmapentries[i]->refstartpos, usexmapentries[j]->refendpos - usexmapentries[j]->refstartpos);
      if(refoverlap >= refminlen){// Case B
	if(fabs(usexmapentries[i]->confidence - usexmapentries[j]->confidence) <= CutFlipFilterConf)
	  continue;

	/* to protect inverted duplications check Query Label Overlap first : If Label Overlap is <= max(0.5,CutFlipMaxOverlap,SmallDupMaxOverlap,OverlapFilterQryLabRatio)
	   no need to check Palindromic inversion ends */
	int U1,U2;
	int cnt = QryLabelOverlap(p,q, U1, U2);
	if(VERB>=2){
	  printf(" i=%i,xid=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d) overlaps in len by %.4f with xid=%i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d), and in aligned query labels by %0.4f (%d/%d)\n", 
		 i,usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,usexmapentries[i]->align->numpairs,
		 overlap/minlen, usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
		 usexmapentries[j]->align->numpairs, (double)cnt / min(U1,U2), cnt, min(U1,U2));
	  fflush(stdout);
	}

	if(cnt <= max(max(0.5, SmallDupMaxOverlap),max(SmallDupMaxOverlap,OverlapFilterQryLabRatio)) * min(U1,U2))
	  continue;

	/* remove palindromic ends of smaller matchgroup  and compute overlap in aligned qry labels */
	int refstartidx, refendidx, qrystartidx, qryendidx, qrystartidx1, qryendidx1, matchcnt, matchcnt1, matchcnt2;
	int refstartidx2 = -1,refendidx2 = -1,qrystartidx2 = -1,qryendidx2 = -1;/* reference and query label ranges for matchgroup j that most tightly bound matchgroup i (after Palindromic end trimming) */
	int tL = -1, tR = -1;/* alignment index values that most closely bound the smaller (inverted matchgroup) on both Query and Reference */
	int M,N,U;
	FLOAT *X,*Y;
	Calign *p, *q;
	if(SmallInversionOverlap(usexmapentries, j, i, p, q, M, N, X, Y, U,
				 refstartidx, refendidx, qrystartidx, qryendidx, qrystartidx1,qryendidx1, matchcnt, matchcnt1, matchcnt2,
				 refstartidx2, refendidx2, qrystartidx2, qryendidx2,tL,tR,false)){
	  if(INV_VERB>=2){
	    printf("refstartidx2=%d,refendidx2=%d, qrystartidx2=%d,qryendidx2=%d,M=%d,N=%d,U=%d: refstartidx=%d, refendidx=%d, qrystartidx1=%d,qryendidx1=%d:matchcnt2=%d, matchcnt=%d, CutFlipMaxOverlap= %0.4f, CutFlipMaxMiss= %0.4f\n",
		   refstartidx2,refendidx2,qrystartidx2,qryendidx2,M,N,U,refstartidx,refendidx,qrystartidx1,qryendidx1,matchcnt2,matchcnt,CutFlipMaxOverlap, CutFlipMaxMiss);
	    printf("\t tL=%d,tR=%d\n",tL,tR);
	    fflush(stdout);
	  }

	  /* check if Query overlap (fraction of aligned labels) exceeds max(0.5, CutFlipMaxOverlap) , ignoring parlindromic ends */
	  if(refstartidx2 >= 1 && refendidx2 >= 1 && matchcnt2 <= floor( matchcnt * max(max(0.5, CutFlipMaxOverlap),OverlapFilterQryLabRatio) + 0.5))/* possible CutFlip Inversion */
	    continue;
	}

      } /* case B */ else {// Case A

	if(fabs(usexmapentries[i]->confidence - usexmapentries[j]->confidence) <= InvOverlapFilterConfDelta)
	  continue;

	if(overlap/minlen < InvOverlapFilterQryRatio)
	  continue;

	int U1,U2;
	int cnt = QryLabelOverlap(p,q,U1,U2);

	if(DEBUG>=2){
	  /* check if Query overlap (fraction of aligned labels) exceeds max(max(0.5, CutFlipMaxOverlap),OverlapFilterQryLabRatio) */
	  int M1 = usexmapentries[i]->querycontignsite;
	  int M2 = usexmapentries[j]->querycontignsite;
	  if(DEBUG && !(M1==M2)){
	    printf(" i=%d(xid=%d,rid=%lld,qid=%lld,M1=%d,mapid1=%d,mapid2=%d,id1=%lld,id2=%lld,N=%d,M=%d)\n j=%d(xid=%d,rid=%lld,qid=%lld,M2=%d,mapid1=%d,mapid2=%d,id1=%lld,id2=%lld,N=%d,M=%d)\n",
		   i,usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, M1, p->mapid1,p->mapid2,
		   YYmap[p->mapid1]->id,XXmap[p->mapid2]->id, YYmap[p->mapid1]->numsite[0],XXmap[p->mapid2]->numsite[0],
		   j,usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, M2, q->mapid1,q->mapid2,
		   YYmap[q->mapid1]->id,XXmap[q->mapid2]->id, YYmap[q->mapid1]->numsite[0],XXmap[q->mapid2]->numsite[0]);
	    fflush(stdout);
	    assert(M1==M2);
	  }
	  
	  assert(U1 == p->numpairs);
	  assert(U2 == q->numpairs);
	  int cnt2 = 0;// number of matching aligned query labels
	  int t1,t2;
	  if(p->orientation){/* matchgroup p is inverted */
	    for(t1 = U1-1, t2 = 0; t1 >= 0 && t2 < U2; ){
	      if(M1+1 - p->sites2[t1] < q->sites2[t2]){
		t1--;
		continue;
	      }
	      if(M1+1 - p->sites2[t1] > q->sites2[t2]){
		t2++;
		continue;
	      }
	      if(DEBUG) assert(M1+1 - p->sites2[t1] == q->sites2[t2]);
	      cnt2++;
	      t1--;
	      t2++;
	    }
	  } else {/* matchgroup q is inverted */
	    for(t1 = 0, t2 = U2-1; t1 < U1 && t2 >= 0; ){
	      if(VERB>=3 ){
		printf("t1=%d,t2=%d:U1=%d,U2=%d,M2=%d:p->sites2[t1]=%d,q->sites2[t2]=%d: cnt=%d\n",
		       t1,t2,U1,U2,M2,p->sites2[t1],q->sites2[t2],cnt2);
		fflush(stdout);
	      }
	      if(p->sites2[t1] < M2+1 - q->sites2[t2]){
		t1++;
		continue;
	      }
	      if(p->sites2[t1] > M2+1 - q->sites2[t2]){
		t2--;
		continue;
	      }
	      if(DEBUG) assert(p->sites2[t1] == M2+1 - q->sites2[t2]);
	      cnt2++;
	      t1++;
	      t2--;
	    }
	  }
      
	  if(!(cnt == cnt2)){
	    printf("cnt2=%d: U1=%d,U2=%d,cnt=%d\n",cnt2,U1,U2,cnt);
	    fflush(stdout);
	    assert(cnt == cnt2);
	  }
	}

	if(VERB>=2){
	  printf("xmapid1=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d) overlaps in qry len by %.4f with xmapid2=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d), and in aligned query labels by %0.4f (%d/%d)\n", 
		 usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,usexmapentries[i]->align->numpairs,
		 overlap/minlen, usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
		 usexmapentries[j]->align->numpairs, (double)cnt / min(U1,U2), cnt, min(U1,U2));
	  if(SMALL_INVERSION_FIX && refminlen > 0.0)
	    printf("\t Ref len overlap fraction is %0.6f (refoverlap= %0.4f kb, refminlen= %0.4f kb)\n",max(0.0,refoverlap)/refminlen, refoverlap, refminlen);
	  fflush(stdout);
	}
	if(cnt <= max(max(0.5, CutFlipMaxOverlap),OverlapFilterQryLabRatio) * min(U1,U2))
	  continue;// protect possible CutFlip inversion
      } // Case A
    }
	     
    if( usexmapentries[i]->confidence < usexmapentries[j]->confidence ) {
      usexmapentries[i]->sv_overlap = true;
      if(VERB>=2){
	printf("xmapid1=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d) is excluded due to overlap (%.2f) with xmapid2=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d)\n", 
	       usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,usexmapentries[i]->align->numpairs,
	       overlap/minlen, usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence, usexmapentries[j]->align->numpairs);
	fflush(stdout);
      }
    } else if( usexmapentries[j]->confidence < usexmapentries[i]->confidence ) {
      usexmapentries[j]->sv_overlap = true;
      if(VERB>=2){
	printf("xmapid1=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f) is excluded due to overlap (%.2f) with xmapid2=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f)\n", 
	       usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
	       overlap/minlen, usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence);
	fflush(stdout);
      }
    }
  } // for t = 0 .. MGlistlen - 1 
}

/* return true of any matching SV with same xmapid1,xmapid2 is present in sv_array[0..numsv-1] or false if not found 
   NOTE : xmapid1 is the MG with larger confidence value than xmapid2 */
bool findSV(int xmapid1, int xmapid2, structuralVariation *sv, int numsv)
{
  if(VERB>=2){
    printf("findSV: xmapid1=%d,xmapid2=%d:sv== &sv_array[%ld], numsv=%d\n", xmapid1,xmapid2,sv-sv_array,numsv);
    fflush(stdout);
  }

  int low = 0, high = numsv -1;
  while(high > low){
    int mid = (low + high)/2;
    structuralVariation *p = &sv[mid];
    if(DEBUG>=2){
      assert(p->xmapid1 != p->xmapid2);
      assert(!(p->xmapid1 == xmapid2 && p->xmapid2 == xmapid1));
    }
    if(p->xmapid1 < xmapid1)
      low = mid+1;
    else if(p->xmapid1 > xmapid1)
      high = mid-1;
    else {
      if(p->xmapid2 < xmapid2)
	low = mid+1;
      else if(p->xmapid2 > xmapid2)
	high = mid-1;
      else {
	return true;
      }
    }
  }

  return low <= high && sv[low].xmapid1 == xmapid1 && sv[low].xmapid2 == xmapid2;
}

  // Filter out overlapping matchgroups based on -OverlapFilter, which requires the following three conditions to be satisfied for the smaller confidence matchgroup :
  // A1. Query overlap ratio in kb must exceed OverlapFilterQryRatio (0.7 by default)
  // A2. Query overlap ratio in Labels must exceed OverlapFilterQryLabRatio (0 by default) 
  // A3. Confidence of larger matchgroup must exceed that of smaller matchgroup by at least OverlapFilterConfDelta (0 by default)
  // A4. IF OverlapFilterSameRef == 1 : refcontigid's must be the same.
  //
  // However if matchgroups have opposite orientation (or have different refcontigid's) the smaller matchgroup is NOT deleted and flagged with inversion_only=true, so it is only available for inversion calls, 
  //               unless the following 3 conditions are also satisfied :
  // B1. Query overlap ratio in kb must exceed InvOverlapFilterQryRatio (99.9% with different refcontid's)
  // B2. Query overlap ratio in Labels must exceed InvOverlapFilterQryLabRatio (0 by default)
  // B3. Confidence of larger matchgroup must exceed that of smaller matchgroup by at least InvOverlapFilterConfDelta.
  //
  // HERE HERE : InvOverlapFilter should not protect matchgroups with same orientation as dominant (largest) matchgroup of query (with same ref), if this call to OverlapFilterMG() is after applying 
  //	       MG merges (AfterMergeMG == true) since they cannot call a valid inversion, unless it is a nested inversions, which should already have been called by Small Inversion or inverted Duplication.
  //
  // HERE HERE : Order overlapped matchgroup pairs (satisfying conditions A1 & A3) from large to small (cf FilterMG()).
  //
  // If two matchgroups of opposite orientation already call a Small Inversion or Inverted Duplication, then flag the smaller MG with largeMG_only=true, so it can only be used subsequently for
  //    non-indel SV calls (or only inversions & inverted duplication calls if inversion_only is also set) as the larger MG (which will happen rarely). Note that this effectively waives the 
  //    Qry Label Overlap conditions (A2 & B2). Only applied if !AfterMergeMG.
  //
  // If INV_OVERLAP_FILTER_FIX : Alternately (even if AfterMergeMG) and smaller matchgroup is fully overlapped in both Ref and Qry, then at least set smaller MG to inversion_only
  //
  // NOTE : Small Inversions & Duplications can be found at sv_array[numsvIndel1..numsvIndel-1] and are sorted by xmapid1 and xmapid2
  //
  // Only query overlaps are considered here -- reference positions are considered later (Heng's criteria)
  //
  // uses the flag xmapEntry.sv_overlap to flag matchgroups for deletion
  //
  // returns number of matchgroups marked for deletion

int OverlapFilterMG(xmapEntry **usexmapentries, int smapsize, int numsvIndel1, int numsvIndel, bool AfterMergeMG, bool verbose)
{
  int filtercnt = 0;

  for(int i = 0; i < smapsize; i++) {
    if(VERB>=2){
      printf("Checking xmapid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv=%d) for overlap : sv_overlap=%d\n",
	     usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
	     usexmapentries[i]->sv_overlap ? 1 : 0, usexmapentries[i]->inversion_only ? 1 : 0);
      fflush(stdout);
    }

    if( usexmapentries[i]->sv_overlap ) //already excluded (by previous i), no need to recheck
      continue;
    for(int j = i+1; j < smapsize; j++) {
      if( usexmapentries[i]->qrycontigid != usexmapentries[j]->qrycontigid ) //only exclude if same query contig
	break;// since qrycontigid are in ascending order
      if( usexmapentries[j]->sv_overlap ) //already excluded (by previous i), no need to recheck
	continue;

      if(OverlapFilterSameRef && usexmapentries[i]->refcontigid != usexmapentries[j]->refcontigid)
	continue;

      double mini = min(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
      double minj = min(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
      double maxi = max(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
      double maxj = max(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
      double overlap = min(maxi,maxj) - max(mini,minj);

      double minlen = min(maxi-mini, maxj-minj);

      double refoverlap = min(usexmapentries[i]->refendpos, usexmapentries[j]->refendpos) - max(usexmapentries[i]->refstartpos, usexmapentries[j]->refstartpos);
      double refminlen = min(usexmapentries[i]->refendpos - usexmapentries[i]->refstartpos, usexmapentries[j]->refendpos - usexmapentries[j]->refstartpos);

      if(DEBUG>=2){
	double ol = overlapLength(mini, maxi, minj, maxj); //returns 0 for no overlap
	if( ol > 0 && overlap != ol ) {
	  printf("overlap: %f ovlplen: %f; i= %f %f, j= %f %f\n", overlap, ol, mini, maxi, minj, maxj); 
	  fflush(stdout);
	  assert( overlap == overlapLength(mini, maxi, minj, maxj) );
	}
      }
      if( overlap <= 0.0 )
	continue;

      /* check if MG i & j already calls a Small Inversion or Duplications */
      int xmapid1 = usexmapentries[i]->xmapid;
      int xmapid2 = usexmapentries[j]->xmapid;
      if(usexmapentries[i]->confidence < usexmapentries[j]->confidence){
	int tmp = xmapid1;
	xmapid1 = xmapid2;
	xmapid2 = tmp;
      }
      bool SmallSV = !AfterMergeMG && findSV(xmapid1,xmapid2,&sv_array[numsvIndel1],numsvIndel-numsvIndel1);
      if(SmallSV){/* mark smaller MG with largeMG_only = true : only matters if MG is not filtered */
	//	if(DEBUG>=1+RELEASE) assert(fabs(usexmapentries[i]->confidence - usexmapentries[j]->confidence) >= 2.0);

	if(usexmapentries[i]->confidence < usexmapentries[j]->confidence)
	  usexmapentries[i]->largeMG_only = true;
	else
	  usexmapentries[j]->largeMG_only = true;
      }

      if(VERB>=2){
	Calign *p = usexmapentries[i]->align;
	Calign *q = usexmapentries[j]->align;
	int U1,U2;
	int cnt = QryLabelOverlap(p,q, U1, U2);

	printf("Checking xmapid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv=%d,lMG=%d) for overlap with xmapid %i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv=%d,lMG=%d): overlap= %0.3f(%d), minlen= %0.3f(%d) (ratio= %0.4f,%0.4f),SmallSV=%d,AfterMG=%d\n", 
	       usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
	       usexmapentries[i]->inversion_only ? 1 : 0, usexmapentries[i]->largeMG_only ? 1 : 0,
	       usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
	       usexmapentries[j]->inversion_only ? 1 : 0, usexmapentries[j]->largeMG_only ? 1 : 0,
	       overlap*1e-3,cnt,minlen*1e-3,min(U1,U2),overlap/minlen,((double)cnt)/min(U1,U2),SmallSV ? 1 : 0, AfterMergeMG?1:0);
	fflush(stdout);
      }

#if 0
      if(DEBUG && !( overlap/minlen > OverlapFilterQryRatio && fabs(usexmapentries[i]->confidence - usexmapentries[j]->confidence) >= OverlapFilterConfDelta && 
		     usexmapentries[i]->orientfoward == usexmapentries[i]->orientforward) && SmallSV){/* should not be possible to get a Small Duplication in this case */
	Calign *p = usexmapentries[i]->align;
	Calign *q = usexmapentries[j]->align;
	int U1,U2;
	int cnt = QryLabelOverlap(p,q, U1, U2);

	printf("Checking xmapid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv=%d,lMG=%d) for overlap with xmapid %i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv=%d,lMG=%d): overlap= %0.3f(%d), minlen= %0.3f(%d) (ratio= %0.4f,%0.4f),SmallSV=%d,AfterMG=%d\n", 
	       usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
	       usexmapentries[i]->inversion_only ? 1 : 0, usexmapentries[i]->largeMG_only ? 1 : 0,
	       usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
	       usexmapentries[j]->inversion_only ? 1 : 0, usexmapentries[j]->largeMG_only ? 1 : 0,
	       overlap*1e-3,cnt,minlen*1e-3,min(U1,U2),overlap/minlen,((double)cnt)/min(U1,U2),SmallSV ? 1 : 0, AfterMergeMG?1:0);
	printf("SmallSV = %d, AfterMergeMG=%d\n",SmallSV,AfterMergeMG);
	fflush(stdout);
      }
#endif

      if( overlap/minlen > OverlapFilterQryRatio && fabs(usexmapentries[i]->confidence - usexmapentries[j]->confidence) >= OverlapFilterConfDelta ) { //parameters.cpp
	Calign *p = usexmapentries[i]->align;
	Calign *q = usexmapentries[j]->align;
	int U1,U2;
	int cnt = QryLabelOverlap(p,q, U1, U2);

	if(cnt > min(U1,U2) * OverlapFilterQryLabRatio || 
	   (INV_OVERLAP_FILTER_FIX && usexmapentries[i]->orientforward != usexmapentries[j]->orientforward && overlap/minlen >= 0.99 && refoverlap/refminlen >= 0.99)){

	  if(usexmapentries[i]->orientforward != usexmapentries[j]->orientforward || usexmapentries[i]->refcontigid != usexmapentries[j]->refcontigid){
	    /* protect smaller matchgroup overlapped in query by larger opposite orientation matchgroup with same ref contig (or with a diff ref contig), unless each of the following is true :
	       B1. Query Overlap ratio (in kb) exceeds max(OverlapFilterQryRatio,InvOverlapFilterQryRatio) (99.9% for diff ref contig).
	       B2. Query Overlap ratio (in labels) exceeds max(OverlapFilterQryLabRatio,InvOverlapFilterQryLabRatio).
	       B3. Larger matchgroup confidence exceeds smaller matchgroup confidence by max(OverlapFilterConfDelta,InvOverlapFilterConfDelta).

	       This is done because the smaller matchgroup might be involved in inversion calls with some other opposite orientation matchgroup with no significant overlap in query or reference
	       However, flag the smaller matchgroup with inversion = true, so that it is only used for inversion calls with other opposite orientation matchgroups that do NOT overlap on the query by 
	       OverlapFilterQryRatio
	       
	       HERE HERE : InvOverlapFilter should not protect matchgroups with same orientation as dominant (largest) matchgroup of query (with same ref), if AfterMergeMG==true
	       (see MergeMG) since they cannot call a valid inversion, unless it is a nested inversions, which should already have been called by Small Inversion or inverted Duplication : to be sure
	       just mark such MGs with largeMG_only = true.
	    */
	  
	    if(!(overlap / minlen >= ((usexmapentries[i]->refcontigid == usexmapentries[j]->refcontigid) ? max(OverlapFilterQryRatio,InvOverlapFilterQryRatio) : 0.999) && 
		 cnt > min(U1,U2) * max(OverlapFilterQryLabRatio,InvOverlapFilterQryLabRatio) && 
		 fabs(usexmapentries[i]->confidence - usexmapentries[j]->confidence) >= max(OverlapFilterConfDelta,InvOverlapFilterConfDelta))){

	      //printf("inversion_only: i=%i j=%i qryi: %10.1f %10.1f, qryj: %10.1f %10.1f\n", i, j, usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos, usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
	      // flag lower confidence matchgroup as inversion_only, unless one of the matchgroups is already flagged 
	      if( verbose ){
		if(usexmapentries[i]->inversion_only || usexmapentries[j]->inversion_only){
		  if( usexmapentries[i]->confidence < usexmapentries[j]->confidence)
		    printf("xid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv_only=%d) is %s with overlap (%.4f) with xid %i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv_only=%d)\n", 
			   usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
			   usexmapentries[i]->inversion_only ? 1 : 0,  usexmapentries[i]->inversion_only ? "ALREADY excluded" : "NOT excluded even", overlap/minlen, 
			   usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
			   usexmapentries[j]->inversion_only ? 1 : 0);
		  else
		    printf("xid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv_only=%d) is %s with overlap (%.4f) with xid %i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv_only=%d)\n", 
			   usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
			   usexmapentries[j]->inversion_only ? 1 : 0,  usexmapentries[j]->inversion_only ? "ALREADY excluded" : "NOT excluded even",  overlap/minlen, 
			   usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
			   usexmapentries[i]->inversion_only ? 1 : 0);		  
		} else if( usexmapentries[i]->confidence < usexmapentries[j]->confidence)
		  printf("xid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv_only=%d) is excluded (except for inversions) due to overlap (%.4f) with xid %i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv_only=%d)\n", 
			 usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
			 usexmapentries[i]->inversion_only ? 1 : 0,
			 overlap/minlen, usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
			 usexmapentries[j]->inversion_only ? 1 : 0);
		else
		  printf("xid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv_only=%d) is excluded (except for inversions) due to overlap (%.4f) with xid %i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inv_only=%d)\n", 
			 usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
			 usexmapentries[j]->inversion_only ? 1 : 0,
			 overlap/minlen, usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
			 usexmapentries[i]->inversion_only ? 1 : 0);
		fflush(stdout);
	      }

	      if(!usexmapentries[i]->inversion_only && !usexmapentries[j]->inversion_only){// NEW11
		if( usexmapentries[i]->confidence < usexmapentries[j]->confidence){
		  if(DEBUG>=2 && U1 >= U2 + 10){/* larger matchgroup is being deleted ! */
		    printf("WARNING:Larger overlapped matchgroup is being marked inversion_only : U1= %d, U2= %d, qrylen1= %0.4f, qrylen2= %0.4f, conf1= %0.2f, conf2= %0.2f\n", 
			   U1,U2,(maxi-mini)*1e-3,(maxj-minj)*1e-3,usexmapentries[i]->confidence, usexmapentries[j]->confidence);
		    fflush(stdout);
		    //		    assert(U1 < U2 + 10);
		  }
		  usexmapentries[i]->inversion_only = true;
		} else {
		  if(DEBUG>=2 && U2 >= U1 + 5){/* larger matchgroup is being deleted ! */
		    printf("WARNING:Larger overlapped matchgroup is being marked inversion_only : U1= %d, U2= %d, qrylen1= %0.4f, qrylen2= %0.4f, conf1= %0.2f, conf2= %0.2f\n", 
			   U1,U2,(maxi-mini)*1e-3,(maxj-minj)*1e-3,usexmapentries[i]->confidence, usexmapentries[j]->confidence);
		    fflush(stdout);
		    //		    assert(U2 < U1 + 5);
		  }
		  usexmapentries[j]->inversion_only = true;
		}
	      }
	      continue; 
	    }
	  }

	  //equal confidence means neither excluded
	  if( usexmapentries[i]->confidence < usexmapentries[j]->confidence ) {
	    usexmapentries[i]->sv_overlap = true;
	    filtercnt++;
	    if( verbose ){
	      printf("xmapid1=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,ovrl=%d) is excluded due to overlap (%.4f) with xmapid2=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,ovrl=%d)\n", 
		     usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
		     usexmapentries[i]->xmapid_overlapped, overlap/minlen, usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, 
		     usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence, usexmapentries[j]->xmapid_overlapped);
	      fflush(stdout);
	    }
	    if(DEBUG>=2 && U1 > U2){/* larger matchgroup is being deleted ! */
	      printf("WARNING:Larger overlapped matchgroup is being deleted : U1= %d, U2= %d, qrylen1= %0.4f, qrylen2= %0.4f, conf1= %0.2f, conf2= %0.2f\n", 
		     U1,U2,(maxi-mini)*1e-3,(maxj-minj)*1e-3,usexmapentries[i]->confidence, usexmapentries[j]->confidence);
	      fflush(stdout);
	      //	      assert(U1 < U2 + 10);
	    }
	  }
	  else if( usexmapentries[j]->confidence < usexmapentries[i]->confidence ) {
	    usexmapentries[j]->sv_overlap = true;
	    filtercnt++;
	    if( verbose ){
	      printf("xmapid1=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,ovrl=%d) is excluded due to overlap (%.2f) with xmapid2=%d (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,ovrl=%d)\n", 
		     usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence, 
		     usexmapentries[j]->xmapid_overlapped, overlap/minlen, 
		     usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
		     usexmapentries[i]->xmapid_overlapped);
	      fflush(stdout);
	    }
	    if(DEBUG>=2 && U2 > U1){/* larger matchgroup is being deleted ! */
	      printf("WARNING:Larger overlapped matchgroup is being deleted : U1= %d, U2= %d, qrylen1= %0.4f, qrylen2= %0.4f, conf1= %0.2f, conf2= %0.2f\n", 
		     U1,U2,(maxi-mini)*1e-3,(maxj-minj)*1e-3,usexmapentries[i]->confidence, usexmapentries[j]->confidence);
	      fflush(stdout);
	      //	      assert(U2 < U1 + 10);
	    }
	  }
	}
      }
      //printf("%2i %2i, %i %i, %7.3f %7.3f, %.2f (%9.1f %9.1f)\n", usexmapentries[i]->xmapid, usexmapentries[j]->xmapid, usexmapentries[i]->sv_overlap, usexmapentries[j]->sv_overlap, usexmapentries[i]->confidence, usexmapentries[j]->confidence, overlap/minlen, overlap, minlen); //debug -- comment out above break
      if(usexmapentries[i]->sv_overlap)
	break;
    } //end inner loop for simple overlap filtering 
  } //end outer loop for simple overlap filtering 

  return filtercnt;
}


// Compute confidence of a subrange of an existing alignment
double SubrangeAlignFP(int S, int E, Calign *q, Cmap *Ymap, Cmap *Xmap)
{
  if(E <= S)
    return 0.0;

  int N = Ymap->numsite[0];
  FLOAT *Y = Ymap->site[0];
  int M = Xmap->numsite[0];
  FLOAT *X = Xmap->site[0];

  /* recompute local X[] as Xrev[], which may depend on ScaleID and p->orientation */
  double *Xrev = new double[M+2];
  FLOAT scale = ScaleFactor[q->scaleID];
  if(!q->orientation){
    for(int J=0;J <= M+1;J++)
      Xrev[J] = X[J]*scale;
  } else {
    for(int j=0;j <= M+1;j++)
      Xrev[j] = (X[M+1]-X[M+1-j])*scale;
  }

  double ChimScore = OutlierEndBias + OutlierEndPenalty;
  int localtype = (PoutlierEnd > 0.0) ? -3 : -2;

  /* create alignment subrange q->sites*[S..E] */
  Calign *r = new Calign[1];
  copy(r,q,1,1);
  int nU = E - S + 1;
  
  if(S > 0){/* new left endoutlier */
    int I = q->sites1[S];
    int K = q->sitesK1[S];
    int J = q->sites2[S];
    r->iscore[0] = r->outscore[0] = ChimScore + Sm(0,I,K,Y);
    r->Lij1 = I-K;
    r->Lij2 = J;
    r->Lend = localtype;
  }
  
  int t = S;
  
  int lastI = r->sites1[0] = q->sites1[t];
  int lastK = r->sitesK1[0] = q->sitesK1[t];
  int lastJ = r->sites2[0] = q->sites2[t];
  r->noutliers = 0;
  r->maxoutlier = 0.0;
  r->maxoutlierLabels = 0;
  for(int u = 1; ++t <= E; u++){
    r->iscore[u] = q->iscore[t];
    r->outscore[u] = q->outscore[t];
    int I = r->sites1[u] = q->sites1[t];
    int K = r->sitesK1[u] = q->sitesK1[t];
    int J = r->sites2[u] = q->sites2[t];
    if(r->outscore[u] + (FLOAT)0.01 < r->iscore[u]){
      r->noutliers++;
      double deltaY = Yc(Y,I,K) - Yc(Y,lastI,lastK);
      double deltaX = Xrev[J] - Xrev[lastJ];
      if(DEBUG) assert(deltaX > 0.0);
      double delta = fabs(deltaY-deltaX);
      r->maxoutlier = max(delta, r->maxoutlier);
      r->maxoutlierLabels = max(I-K-lastI + J-lastJ - 2, r->maxoutlierLabels);
    }
    lastI = I;
    lastK = K;
    lastJ = J;
  }
      
  if(E < q->numpairs - 1){/* new right end */
    int I = q->sites1[E];
    //	int K = q->sitesK1[E];
    int J = q->sites2[E];
    r->iscore[nU] = r->outscore[nU] = ChimScore;
    r->Rij1 = I;
    r->Rij2 = J;
    r->Rend = localtype;
  } else if(S > 0)
    r->iscore[nU] = r->outscore[nU] = r->iscore[r->numpairs];
  
  r->numpairs = nU;
  
  r->score = 0;
  for(int t = 0; t <= nU; t++)
    r->score += r->iscore[t];
      
  double conf = alignFP(r, Y, Xrev, N, M, r->orientation, r->mapid1, r->mapid2, r->score, 0);
  
  delete [] Xrev;
  delete [] r;

  return conf;
}


void output_smap(char *basename, char * &queryfilename)
{
   if(colors>=2){
     printf("output_smap() not yet implemented for colors=%d\n",colors);
     exit(1);
   }

   if(strstr(basename,"/dev/null"))
     return;

   bool verbose = (numxmapentries <= 50000) ? true : false; //local for debugging of this file only

   char filename[PATH_MAX];
   strcpy(filename,basename);
   int i = strlen(filename);
   strcpy(&filename[i],".smap");

   //also make the xmap file name (this will be the filtered xmap)
   char xfilename[PATH_MAX];
   strcpy(xfilename,basename);
   strcpy(&xfilename[i],".xmap"); //i is same as above

   if(checkFile(filename))
     return;

   if(VERB){
     printf("Starting SV analysis, numxmapentries= %d (wall time=%.6f)\n", numxmapentries, wtime());
     fflush(stdout);
   }

   if(DEBUG>=2){/* verify that xmapentries[] are in ascending order of xmapid */
     for(int i = 0; i < numxmapentries - 1; i++){
       if(!(xmapentries[i]->xmapid < xmapentries[i+1]->xmapid)){
	 printf("i=%d:xmapentries[i]->xmapid=%d,xmapentries[i+1]->xmapid=%d\n",i,xmapentries[i]->xmapid,xmapentries[i+1]->xmapid);
	 fflush(stdout);
	 assert(xmapentries[i]->xmapid < xmapentries[i+1]->xmapid);
       }
     }
   }

   // process the alignments: alignment and numaligns are declared in globals.h
   if(DEBUG) assert(numxmapentries <= maxxmapentries);
   if(DEBUG) assert(usexmapentries == NULL);
   usexmapentries = new xmapEntry*[maxxmapentries];
   if(VERB/* HERE HERE >=2 */){
     printf("Allocated usexmapentries[] : maxxmapentries = %d, usexmapentries = %p, xmapentries = %p: wall time= %0.6f secs\n", maxxmapentries, usexmapentries, xmapentries,wtime());
     fflush(stdout); 
   }

   /* initialize usexmapentries[0..smapsiz-1] to include all xmap entries */ 
   smapsize = numxmapentries; //this is not actually the smap size, it's the len of usexmapentries
   for(int i=0; i < numxmapentries; i++){
     if(DEBUG) assert(xmapentries[i]->output == false);
     usexmapentries[i] = xmapentries[i];
     if(DEBUG) assert(usexmapentries[i]->sv_overlap == false);
     if(DEBUG) assert(usexmapentries[i]->inversion_only == false);

     usexmapentries[i]->callRepeat(verbose,REPEAT_TOLERANCE,REPEAT_MINELE,REPEAT_CONF_PENALTY); //check this entry for repeats--DEBUG ONLY
     if(VERB/* HERE HERE HERE >=2 */){
       xmapEntry *p = xmapentries[i];
       printf("i=%d/%d: xmapentries[i]:xmapid=%d:refid=%lld,qryid=%lld,fwd=%d:qrystartpos=%0.1f,qryendpos=%0.1f,refstartpos=%0.1f,refendpos=%0.1f,conf=%0.2f,repeat=%d,simple_repeat=%d\n",
	      i,numxmapentries,p->xmapid,p->refcontigid,p->qrycontigid,p->orientforward ? 1 : 0, p->qrystartpos,p->qryendpos,p->refstartpos, p->refendpos, p->confidence,p->align->repeat,p->have_ref_repeat);
       fflush(stdout);
     }
   }

   if(DEBUG>=2){
     for(int i = 0; i < smapsize;i++){
       xmapEntry *p = usexmapentries[i];

       /* check p->refsites[] (bp) vs refendpos (bp) */
       assert(p->numrefsites >= p->refendidx - p->refstartidx + 1);
       int Istart = p->refstartidx;
       for(int I = Istart; I <= p->refendidx; I++){
	 double sitesI = p->refsites[I - Istart];
	 if(DEBUG && !(p->refstartpos <= sitesI && sitesI <= p->refendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refstartpos= %0.1f,refendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart, p->refsites[I-Istart],p->refstartpos,p->refendpos);
	   fflush(stdout);
	   assert(p->refstartpos <= sitesI && sitesI <= p->refendpos);
	 }
	 if(I > p->refstartidx){
	   if(DEBUG && !(p->refsites[I-1 - Istart] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refsites[%d]= %0.1f, refendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart,p->refsites[I - Istart],I-1-Istart,p->refsites[I-1-Istart],p->refendpos);
	     fflush(stdout);
	     assert(p->refsites[I-1-Istart] <= sitesI);
	   }
	 }
       }
       /* check p->qrysites[0..numqrysites-1] (bp) vs refendpos (bp) */
       double qrystartpos = min(p->qrystartpos,p->qryendpos);
       double qryendpos = max(p->qrystartpos,p->qryendpos);
       int Imin = min(p->qrystartidx,p->qryendidx);
       int Imax = max(p->qrystartidx,p->qryendidx);
       if(DEBUG) assert(p->numqrysites >= Imax-Imin+1);
       for(int I = Imin; I <= Imax; I++){
	 double sitesI = p->qrysites[I - Imin];
	 if(DEBUG && !(qrystartpos <= sitesI && sitesI <= qryendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites, I, I-Imin,p->qrysites[I-Imin],qrystartpos,qryendpos);
	   fflush(stdout);
	   assert(qrystartpos <= sitesI && sitesI <= qryendpos);
	 }
	 if(I > Imin){
	   if(DEBUG && !(p->qrysites[I-1-Imin] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites,I,I-Imin,p->qrysites[I-Imin],I-1-Imin,p->qrysites[I-1-Imin],qrystartpos,qryendpos);
	     fflush(stdout);
	     assert(p->qrysites[I-1-Imin] <= sitesI);
	   }
	 }
       }
     }
   }

   if(DEBUG>=2){
     for(int i = 0; i < numxmapentries;i++){
       xmapEntry *p = usexmapentries[i];
       if(!XXmap[p->align->mapid2]->origmap) assert(p->querycontignsite == XXmap[p->align->mapid2]->numsite[0]);
     }
   }
   if(DEBUG>=2){
     for(int i = 0; i < numsv; i++){
       structuralVariation *p = &sv_array[i];
       if(p->duplicate || p->confidence == 0.0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	 continue;

       xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
       xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
       if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	 printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
		i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	 fflush(stdout);
	 assert(entry1->xmapid == p->xmapid1);
	 assert(entry2->xmapid == p->xmapid2);
       }
     }
   }

   int numsvIndel1 = numsv;/* original number of SVs with Indels with single Matchgroup */

   // sort Indel SVs in order of xmapid1,xmapid2 etc in smap : used to speed up Indel removal in checkOverlap() etc
   qsort(sv_array, numsv, sizeof(*sv_array), compare_svs);

   MaxXmapid = (numxmapentries <= 0) ? 0 : xmapentries[numxmapentries-1]->xmapid;
   delete [] firstIndel;
   firstIndel = new int[MaxXmapid+1];
   for(int i = 0; i <= MaxXmapid; i++)
     firstIndel[i] = MaxXmapid + 1;

   int Xmapid = 0;
   for(int i = 0; i < numsv; i++){
     structuralVariation *pSV = &sv_array[i];
     if(DEBUG>=2) assert(pSV->algo_enum == sv_indel);
     if(DEBUG>=2) assert(pSV->indices[0] == pSV->indices[1]);
     xmapEntry *entry1 = xmapentries[pSV->indices[0]];
     if(entry1->xmapid > Xmapid){
       Xmapid = entry1->xmapid;
       firstIndel[Xmapid] = i;
     }
   }

   /* first mark matchgroups with 1-MG Indels with output == true */
   for(int i = 0; i < numsvIndel1; i++){
     structuralVariation *p = &sv_array[i];
     if(DEBUG) assert(p->algo_enum != sv_ps);
     if(DEBUG) assert(p->indices[0] == p->indices[1]);
     if(DEBUG) assert(p->duplicate == false);

     if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
       continue;

     xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
     entry1->output = true;
   }

   /* sort usexmapentries[0.. smapsize-1] by qrycontigid, refcontigid, refstartpos, refendpos in ascending order
      NOTE : this does not change the order of xmapentries[] which remains in ascending order of xmapid, so output is not affected */
   qsort(usexmapentries, smapsize, sizeof(usexmapentries[0]), (intcmp *)xmapEntryQryidInc);

   if(VERB){
     printf("Finished sorting usexmapentries[]:smapsiz=%d: wall time= %0.6f secs\n",smapsize,wtime());
     fflush(stdout);
   }

   if(DEBUG>=2){
     for(int i = 0; i < numxmapentries;i++){
       xmapEntry *p = usexmapentries[i];
       if(!XXmap[p->align->mapid2]->origmap) assert(p->querycontignsite == XXmap[p->align->mapid2]->numsite[0]);
     }
   }

   if(VERB>=2){
     printf("\nAfter sorting usexmapentries[]: numxmapentries=%d\n",numxmapentries);
     for(int i = 0; i < numxmapentries; i++){
       xmapEntry *p = usexmapentries[i];
       printf("i=%d/%d: xmapentries[i]:xmapid=%d:refid=%lld,qryid=%lld,fwd=%d:qrystartpos=%0.1f,qryendpos=%0.1f,refstartpos=%0.1f,refendpos=%0.1f,conf=%0.2f:sv_overlap=%d\n",
	      i,numxmapentries,p->xmapid,p->refcontigid,p->qrycontigid,p->orientforward ? 1 : 0, p->qrystartpos,p->qryendpos,p->refstartpos, p->refendpos, p->confidence,p->sv_overlap);
       fflush(stdout);
     }
   }
   if(DEBUG>=2){
     for(int i = 0; i < numsv; i++){
       structuralVariation *p = &sv_array[i];
       if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	 continue;

       xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
       xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
       if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	 printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
		i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	 fflush(stdout);
	 assert(entry1->xmapid == p->xmapid1);
	 assert(entry2->xmapid == p->xmapid2);
       }
     }
   }

   if(DEBUG>=2)
     for(int i = 0; i < smapsize; i++)
       assert(usexmapentries[i]->numrefsites >= usexmapentries[i]->refendidx - usexmapentries[i]->refstartidx + 1);

   /* If smaller matchgroup is fully overlapped on qry (and ref if MATCHGROUP_TRIM <= 1 and has same ref contig if MATCHGROUP_TRIM <= 2),
      first trim the smaller matchgroup so it only overlaps qry at outliers of smap_Outlier_minQrySites labels or more in
      larger matchgroup. May need to split smaller MG if there are multiple such outliers in larger MG. Delete smaller MG if there are no such outliers.
      Finally resort matchgroups in usexmapentries[] if there are new or deleted MGs */
   if(MATCHGROUP_TRIM && XmapUnique){
     size_t overlapcnt = 0;
     int orignumxmapentries = numxmapentries;

     for(int i = 0; i < smapsize; i++) {    
       xmapEntry *pMGi = usexmapentries[i];
       if( pMGi->sv_overlap ) //already excluded (by previous i), no need to recheck
	 continue;

       for(int j = i+1; j < smapsize; j++) {
	 xmapEntry *pMGj = usexmapentries[j];
	 if( pMGj->sv_overlap ) //already excluded (by previous i), no need to recheck
	   continue;

	 if( pMGi->qrycontigid != pMGj->qrycontigid ) //only exclude if same query contig
	   break;// since qrycontigid are in ascending order

	 /* NOTE : If MATCHGROUP_TRIM <= 2 && pMGi->refcontigid != pMGj->refcontigid : checkOverlap will not filter out the smaller MG if it does not overlap an outlier
		   However, if it does overlap an outlier it will be trimmed to just the outlier region */
	 /*	if(MATCHGROUP_TRIM <= 2 && pMGi->refcontigid != pMGj->refcontigid)
		 continue;*/

	 // HERE HERE : if translocation, avoid filtering MGs that are repeats

	 double mini = min(pMGi->qrystartpos, pMGi->qryendpos);
	 double minj = min(pMGj->qrystartpos, pMGj->qryendpos);
	 double maxi = max(pMGi->qrystartpos, pMGi->qryendpos);
	 double maxj = max(pMGj->qrystartpos, pMGj->qryendpos);
	 double overlap = min(maxi,maxj) - max(mini,minj);
	 double minlen = min(maxi-mini, maxj-minj);
	 double maxlen = max(maxi-mini, maxj-minj);

	 if(VERB>=2){
	   printf("Checking qry overlap for i=%d,j=%d:xmapid=%d,%d:refid=%lld,%lld,qryid=%lld:conf=%0.2f,%0.2f,qryi= %0.4f..%0.4f, qryj= %0.4f..%0.4f, overlap= %0.4f, minlen= %0.4f, maxlen= %0.4f kb\n",
		  i,j,pMGi->xmapid,pMGj->xmapid,pMGi->refcontigid,pMGj->refcontigid,pMGi->qrycontigid,pMGi->confidence,pMGj->confidence,
		  mini*1e-3,maxi*1e-3,minj*1e-3,maxj*1e-3,overlap*1e-3,minlen*1e-3,maxlen*1e-3);
	   fflush(stdout);
	 }

	 /* HERE HERE HERE : if both matchgroups are same size and opposite orientation (palindromic inversion or repeat) keep orientation that matches dominant orientation of nearby matchgroups 
	    For same orientation, for matchgroups of same size, keep matchgroup that is closest on reference to other matchgroups
	  */
	 if(overlap < 0.999 * minlen /* WAS88  || (overlap >= 0.99 * maxlen) */)
	   continue;

	 /* check reference overlap */
	 double roverlap = min(pMGi->refendpos, pMGj->refendpos) - max(pMGi->refstartpos, pMGj->refstartpos);

	 /* NEW77 : confidence and size must agree, otherwise not safe to trim/remove overlapped MG */
	 if(pMGi->confidence < pMGj->confidence && maxi - mini <= /* WAS88 < */ maxj - minj){
	   if(MATCHGROUP_TRIM <= 1 && roverlap < 0.999 * (pMGi->refendpos - pMGi->refstartpos))
	     continue;

	   overlapcnt++;
	   (void)checkOverlap(usexmapentries,smapsize,xmapentries,numxmapentries,maxxmapentries,j,i,true,true);
	 } else if(pMGi->confidence > pMGj->confidence && maxi - mini >= /* WAS88 > */ maxj - minj) {
	   if(MATCHGROUP_TRIM <= 1 && roverlap < 0.999 * (pMGj->refendpos - pMGj->refstartpos))
	     continue; // Small Duplication should only apply when smaller matchgroup is fully overlapped on reference by larger matchgroup	    

	   overlapcnt++;
	   (void)checkOverlap(usexmapentries,smapsize,xmapentries,numxmapentries,maxxmapentries,i,j,true,true);
	 }

	 if( pMGi->sv_overlap ) // NEW104 : no need to check remaining values of j
	   break;
       }
     }

     int newsmapsize = 0;
     for( int i=0; i < smapsize; i++) {
       if(usexmapentries[i]->sv_overlap)
	 continue;
       usexmapentries[newsmapsize++] = usexmapentries[i];
     }
     if( VERB ){
       printf("After query overlap filter (overlapcnt= %lu): numxmapentries= %d -> %d -> %d (wall time= %0.6f)\n", overlapcnt, orignumxmapentries, smapsize, newsmapsize, wtime());
       fflush(stdout);
     }
     smapsize = newsmapsize; //this is len of usexmapentries

     qsort(usexmapentries, smapsize, sizeof(usexmapentries[0]), (intcmp *)xmapEntryQryidInc);
     if(DEBUG>=2){
       for(int i = 0; i < numxmapentries;i++){
	 xmapEntry *p = usexmapentries[i];
	 if(!XXmap[p->align->mapid2]->origmap) assert(p->querycontignsite == XXmap[p->align->mapid2]->numsite[0]);
       }
     }

     if(VERB){
       printf("After resorting usexmapentries[]: smapsize=%d: wall time= %0.6f secs\n",smapsize,wtime());
       fflush(stdout);
     }

     if(VERB>=3){
       for(int i = 0; i < smapsize; i++){
	 xmapEntry *p = usexmapentries[i];
	 if(1 /* p->xmapid == 868 || p->xmapid == 866 */){
	   printf("usexmapentries[%d]: refid=%lld,qryid=%lld, xmapid=%d, xmapid_overlapped=%d : fwd=%d,conf=%0.2f,np=%d,qry=%0.3f..%0.3f\n",
		  i, p->refcontigid,p->qrycontigid,p->xmapid,p->xmapid_overlapped,p->orientforward ? 1 : 0, p->confidence, p->align->numpairs,
		  min(p->qrystartpos,p->qryendpos),max(p->qrystartpos,p->qryendpos));
	   fflush(stdout);
	 }
       }
     }
   }

   if(DEBUG>=2){/* check if any 1-MG Indels no longer overlap the MGs they refer to */
     for(int i = 0; i < numsvIndel1; i++){
       structuralVariation *p = &sv_array[i];
       if(DEBUG) assert(p->algo_enum != sv_ps);
       if(DEBUG && i < numsvIndel1) assert(p->indices[0] == p->indices[1]);

       if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	 continue;

       xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
       if(DEBUG) assert(entry1->output == true);
       if(entry1->output == false)
	 continue;

       double qmin = min(entry1->qrystartpos,entry1->qryendpos);
       double qmax = max(entry1->qrystartpos,entry1->qryendpos);
       int iqmin = min(entry1->qrystartidx,entry1->qryendidx);
       int iqmax = max(entry1->qrystartidx,entry1->qryendidx);

       if(!((qmin <= p->querystart && p->querystart <= qmax) &&
	    (qmin <= p->querystop && p->querystop <= qmax) &&
	    (iqmin <= p->querystartidx && p->querystartidx <= iqmax) &&
	    (iqmin <= p->querystopidx && p->querystopidx <= iqmax) &&
	    (entry1->refstartpos <= p->refstart && p->refstart <= entry1->refendpos) &&
	    (entry1->refstartpos <= p->refstop && p->refstop <= entry1->refendpos) &&
	    (entry1->refstartidx <= p->refstartidx && p->refstartidx <= entry1->refendidx) &&
	    (entry1->refstartidx <= p->refstopidx && p->refstopidx <= entry1->refendidx))){
	 printf("sv_array[%d]:Single MG indel's MG was trimmed and no longer includes the indel:\n",i);
	 printf("sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
		p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
		p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	 printf("xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
		entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
		entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
		entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	 fflush(stdout);
	 assert(0);
       }
     }
   }

   if(DEBUG>=2)
     for(int i = 0; i < smapsize; i++)
       assert(usexmapentries[i]->numrefsites >= usexmapentries[i]->refendidx - usexmapentries[i]->refstartidx + 1);

   if(DEBUG>=2){
     for(int i = 0; i < smapsize;i++){
       xmapEntry *p = usexmapentries[i];

       /* check p->refsites[] (bp) vs refendpos (bp) */
       assert(p->numrefsites >= p->refendidx - p->refstartidx + 1);
       int Istart = p->refstartidx;
       for(int I = Istart; I <= p->refendidx; I++){
	 double sitesI = p->refsites[I - Istart];
	 if(DEBUG && !(p->refstartpos <= sitesI && sitesI <= p->refendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refstartpos= %0.1f,refendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart, p->refsites[I-Istart],p->refstartpos,p->refendpos);
	   fflush(stdout);
	   assert(p->refstartpos <= sitesI && sitesI <= p->refendpos);
	 }
	 if(I > p->refstartidx){
	   if(DEBUG && !(p->refsites[I-1 - Istart] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refsites[%d]= %0.1f, refendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart,p->refsites[I - Istart],I-1-Istart,p->refsites[I-1-Istart],p->refendpos);
	     fflush(stdout);
	     assert(p->refsites[I-1-Istart] <= sitesI);
	   }
	 }
       }
       /* check p->qrysites[0..numqrysites-1] (bp) vs refendpos (bp) */
       double qrystartpos = min(p->qrystartpos,p->qryendpos);
       double qryendpos = max(p->qrystartpos,p->qryendpos);
       int Imin = min(p->qrystartidx,p->qryendidx);
       int Imax = max(p->qrystartidx,p->qryendidx);
       if(DEBUG) assert(p->numqrysites >= Imax-Imin+1);
       for(int I = Imin; I <= Imax; I++){
	 double sitesI = p->qrysites[I - Imin];
	 if(DEBUG && !(qrystartpos <= sitesI && sitesI <= qryendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites, I, I-Imin,p->qrysites[I-Imin],qrystartpos,qryendpos);
	   fflush(stdout);
	   assert(qrystartpos <= sitesI && sitesI <= qryendpos);
	 }
	 if(I > Imin){
	   if(DEBUG && !(p->qrysites[I-1-Imin] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites,I,I-Imin,p->qrysites[I-Imin],I-1-Imin,p->qrysites[I-1-Imin],qrystartpos,qryendpos);
	     fflush(stdout);
	     assert(p->qrysites[I-1-Imin] <= sitesI);
	   }
	 }
       }
     }
   }

   /* initialize contigendalign[0..endalignsize-1] */
   endalign *contigendalign = new endalign[nummaps*2]; //fixed size
   int endalignsize = 0; //number of entries which are used

   /* Use hashtable to map from qrycontigid (64 bit) to index into contigendalign[0..nummaps*2-1] : will point to first of two consecutive entries in contigendalign[] */
   HashInit();

   for(int i=0; i < smapsize; i++){
     //printf("Align %2i %4i %.1f: %3i %3i %7.3f %7.3f\n", xmapentries[i]->xmapid, xmapentries[i]->qrycontigid, xmapentries[i]->querycontiglen, xmapentries[i]->startunalign, xmapentries[i]->endunalign, startunalignlen, endunalignlen); 

     endalignsize = mergeEndaligns(contigendalign, endalignsize, usexmapentries[i]);
     if(DEBUG && !(endalignsize <= nummaps*2) ) {
       printf("ERROR in mergeEndaligns: endalignsize=%i nummaps=%i\n", endalignsize, nummaps);
       fflush(stdout);
       assert( endalignsize <= nummaps*2 ); 
     }
   }

   if(VERB/* HERE >=2 */){
     printf("finished merging end alignments from %d maps (smapsiz=%d) : wall time= %0.6f secs\n", nummaps, smapsize, wtime());
     fflush(stdout);
   }

   if(DEBUG>=2){
     for(int i = 0; i < smapsize;i++){
       xmapEntry *p = usexmapentries[i];

       /* check p->refsites[] (bp) vs refendpos (bp) */
       assert(p->numrefsites >= p->refendidx - p->refstartidx + 1);
       int Istart = p->refstartidx;
       for(int I = Istart; I <= p->refendidx; I++){
	 double sitesI = p->refsites[I - Istart];
	 if(DEBUG && !(p->refstartpos <= sitesI && sitesI <= p->refendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refstartpos= %0.1f,refendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart, p->refsites[I-Istart],p->refstartpos,p->refendpos);
	   fflush(stdout);
	   assert(p->refstartpos <= sitesI && sitesI <= p->refendpos);
	 }
	 if(I > p->refstartidx){
	   if(DEBUG && !(p->refsites[I-1 - Istart] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refsites[%d]= %0.1f, refendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart,p->refsites[I - Istart],I-1-Istart,p->refsites[I-1-Istart],p->refendpos);
	     fflush(stdout);
	     assert(p->refsites[I-1-Istart] <= sitesI);
	   }
	 }
       }
       /* check p->qrysites[0..numqrysites-1] (bp) vs refendpos (bp) */
       double qrystartpos = min(p->qrystartpos,p->qryendpos);
       double qryendpos = max(p->qrystartpos,p->qryendpos);
       int Imin = min(p->qrystartidx,p->qryendidx);
       int Imax = max(p->qrystartidx,p->qryendidx);
       if(DEBUG) assert(p->numqrysites >= Imax-Imin+1);
       for(int I = Imin; I <= Imax; I++){
	 double sitesI = p->qrysites[I - Imin];
	 if(DEBUG && !(qrystartpos <= sitesI && sitesI <= qryendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites, I, I-Imin,p->qrysites[I-Imin],qrystartpos,qryendpos);
	   fflush(stdout);
	   assert(qrystartpos <= sitesI && sitesI <= qryendpos);
	 }
	 if(I > Imin){
	   if(DEBUG && !(p->qrysites[I-1-Imin] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites,I,I-Imin,p->qrysites[I-Imin],I-1-Imin,p->qrysites[I-1-Imin],qrystartpos,qryendpos);
	     fflush(stdout);
	     assert(p->qrysites[I-1-Imin] <= sitesI);
	   }
	 }
       }
     }
   }

   if(CutFlipOverlapFilter > 0.0){/* Apply partial filtering of overlapped matchgroups : Similar to less stringent loop after Small Inversions & Duplications, 
				     but some matchgroups are protected (not filtered) and overlap fraction required may be higher */

     // HERE HERE : If Query Overlap ratio (in kb) is 99.9%, keep lower confidence matchgroup only if it overlaps an outlier of 4 (Small_minQrySites) qry label intervals or more in larger matchgroup 

     if(VERB){
       printf("Apply partial matchgroup overlap filter before checking for small SegDup & CutFlip inversions : Query Overlap ratio must be %0.4f (in kb) AND %0.4f (in labels)\n",
	      max(CutFlipOverlapFilter,OverlapFilterQryRatio), max(max(0.5,SmallDupMaxOverlap),OverlapFilterQryLabRatio));
       printf("\t (but %0.4f for inversion qry labels, and ignoring palindromic part of inversion with full ref overlap) AND confidence delta must be %0.2f (%0.2f for inversion with full ref overlap)\n",
	      max(max(0.5,CutFlipMaxOverlap),OverlapFilterQryLabRatio), max(OverlapFilterConfDelta, InvOverlapFilterConfDelta), 
	      max(CutFlipFilterConf,max(OverlapFilterConfDelta, InvOverlapFilterConfDelta)));
       fflush(stdout);
     }

     int MGlistmax = 64;
     CmatchgroupPair *MGlist = new CmatchgroupPair[MGlistmax];

     /* NOTE : matchgroups are handled in groups that correspond to the same refcontigid and qrycontigid : MGlist[0..MGlistlen-1]
	Overlapped pairs in the same groups are processed in descending order of smaller,larger matchgroups, so largest matchgroups are filtered out first (avoiding filtering out smaller matchgroups) */

     long long MGrefid = -1, MGqryid = -1;/* refcontigid and qrycontigid of current group */
     int MGlistlen = 0;

     for(int i = 0; i < smapsize; i++) {
       if(DEBUG>=2){
	 xmapEntry *p = usexmapentries[i];
	 if(!XXmap[p->align->mapid2]->origmap) assert(p->querycontignsite == XXmap[p->align->mapid2]->numsite[0]);
       }
       if( usexmapentries[i]->sv_overlap ) //already excluded (by previous i), no need to recheck
	 continue;

       if(!(usexmapentries[i]->refcontigid == MGrefid && usexmapentries[i]->qrycontigid == MGqryid)){/* start of new refcontigid,qrycontigid group */
	 if(MGlistlen > 0)/* process current group of overlapped matchgroup pairs in descending order of sizes */
	   FilterMG(MGlist,MGlistlen,usexmapentries,smapsize);

	 MGlistlen = 0;
	 MGrefid = usexmapentries[i]->refcontigid;
	 MGqryid = usexmapentries[i]->qrycontigid;
       }

       for(int j = i+1; j < smapsize; j++) {
	 if(DEBUG>=2){
	   xmapEntry *p = usexmapentries[j];
	   if(!XXmap[p->align->mapid2]->origmap) assert(p->querycontignsite == XXmap[p->align->mapid2]->numsite[0]);
	 }

	 if( usexmapentries[j]->sv_overlap ) //already excluded (by previous i), no need to recheck
	   continue;

	 if( usexmapentries[i]->qrycontigid != usexmapentries[j]->qrycontigid ) //only exclude if same query contig
	   break;// since qrycontigid are in ascending order

	 if(usexmapentries[i]->refcontigid != usexmapentries[j]->refcontigid)
	   continue;

	 /* protect smaller matchgroup overlapped by larger opposite orientation matchgroup for same query and ref contig : might contain CutFlip inversion : This match group will be filtered out again later */

	 double mini = min(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
	 double minj = min(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
	 double maxi = max(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
	 double maxj = max(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
	 double overlap = min(maxi,maxj) - max(mini,minj);
	 double minlen = min(maxi-mini, maxj-minj);

	 if(VERB>=2 && overlap > 0.0){
	   printf(" i=%d (xid=%d, rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d,qry=%0.3f..%0.3f), j=%d (xid=%d,rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,np=%d,qry=%0.3f..%0.3f) : overlap= %0.3f, minlen= %0.3f\n",
		  i, usexmapentries[i]->xmapid,usexmapentries[i]->refcontigid,usexmapentries[i]->qrycontigid,usexmapentries[i]->orientforward ? 1 : 0, 
		  usexmapentries[i]->confidence, usexmapentries[i]->align->numpairs, mini*1e-3, maxi * 1e-3,
		  j, usexmapentries[j]->xmapid,usexmapentries[j]->refcontigid,usexmapentries[j]->qrycontigid,usexmapentries[j]->orientforward ? 1 : 0, 
		  usexmapentries[j]->confidence, usexmapentries[j]->align->numpairs, minj*1e-3, maxj * 1e-3, overlap * 1e-3, minlen * 1e-3);
	   fflush(stdout);
	 }

	 if( overlap <= 0.0 )
	   continue;

	 if( overlap/minlen >= max(CutFlipOverlapFilter, OverlapFilterQryRatio) &&
	     fabs(usexmapentries[i]->confidence - usexmapentries[j]->confidence) > max(InvOverlapFilterConfDelta, OverlapFilterConfDelta) ) { // see parameters.cpp
	   if(MGlistlen >= MGlistmax){
	     CmatchgroupPair *newMGlist = new CmatchgroupPair[MGlistmax *= 2];
	     memcpy(newMGlist,MGlist,MGlistlen * sizeof(CmatchgroupPair));
	     delete [] MGlist;
	     MGlist = newMGlist;
	   }
	   CmatchgroupPair *pMG = &MGlist[MGlistlen++];
	   pMG->i = i;
	   pMG->j = j;
	   pMG->iLabels = usexmapentries[i]->align->numpairs;
	   pMG->jLabels = usexmapentries[j]->align->numpairs;

	   continue;
	 }
       } // end inner loop for simple overlap filtering 
     } // end outer loop for simple overlap filtering 

     if(VERB/* HERE HERE >=2 */){
       printf("Final MGlistlen=%d\n",MGlistlen);
       fflush(stdout);
     }

     if(MGlistlen > 0)
       FilterMG(MGlist,MGlistlen,usexmapentries,smapsize);      
     delete [] MGlist; MGlist = NULL;

     // filter usexmapentries based on sv_overlap
     int newsmapsize=0;
     for( int i=0; i < smapsize; i++) {
       if(usexmapentries[i]->sv_overlap)
	 continue;
       usexmapentries[newsmapsize++] = usexmapentries[i];
     }
     if( VERB ){
       printf("Xmap entries removed due to pre-inversion overlap: %i, Xmaps remaining: %i (wall time= %0.6f)\n", smapsize-newsmapsize, newsmapsize, wtime());
       fflush(stdout);
     }
     smapsize = newsmapsize; //this is len of usexmapentries
   }

   if(DEBUG>=2){
     for(int i = 0; i < smapsize;i++){
       xmapEntry *p = usexmapentries[i];

       /* check p->refsites[] (bp) vs refendpos (bp) */
       assert(p->numrefsites >= p->refendidx - p->refstartidx + 1);
       int Istart = p->refstartidx;
       for(int I = Istart; I <= p->refendidx; I++){
	 double sitesI = p->refsites[I - Istart];
	 if(DEBUG && !(p->refstartpos <= sitesI && sitesI <= p->refendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refstartpos= %0.1f,refendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart, p->refsites[I-Istart],p->refstartpos,p->refendpos);
	   fflush(stdout);
	   assert(p->refstartpos <= sitesI && sitesI <= p->refendpos);
	 }
	 if(I > p->refstartidx){
	   if(DEBUG && !(p->refsites[I-1 - Istart] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refsites[%d]= %0.1f, refendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart,p->refsites[I - Istart],I-1-Istart,p->refsites[I-1-Istart],p->refendpos);
	     fflush(stdout);
	     assert(p->refsites[I-1-Istart] <= sitesI);
	   }
	 }
       }
       /* check p->qrysites[0..numqrysites-1] (bp) vs refendpos (bp) */
       double qrystartpos = min(p->qrystartpos,p->qryendpos);
       double qryendpos = max(p->qrystartpos,p->qryendpos);
       int Imin = min(p->qrystartidx,p->qryendidx);
       int Imax = max(p->qrystartidx,p->qryendidx);
       if(DEBUG) assert(p->numqrysites >= Imax-Imin+1);
       for(int I = Imin; I <= Imax; I++){
	 double sitesI = p->qrysites[I - Imin];
	 if(DEBUG && !(qrystartpos <= sitesI && sitesI <= qryendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites, I, I-Imin,p->qrysites[I-Imin],qrystartpos,qryendpos);
	   fflush(stdout);
	   assert(qrystartpos <= sitesI && sitesI <= qryendpos);
	 }
	 if(I > Imin){
	   if(DEBUG && !(p->qrysites[I-1-Imin] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites,I,I-Imin,p->qrysites[I-Imin],I-1-Imin,p->qrysites[I-1-Imin],qrystartpos,qryendpos);
	     fflush(stdout);
	     assert(p->qrysites[I-1-Imin] <= sitesI);
	   }
	 }
       }
     }
   }

   if(DEBUG>=2){
     for(int i = 0; i < numsv; i++){
       structuralVariation *p = &sv_array[i];
       if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	 continue;

       xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
       xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
       if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	 printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
		i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	 printf("\t numsv=%d\n", numsv);
	 fflush(stdout);
	 assert(entry1->xmapid == p->xmapid1);
	 assert(entry2->xmapid == p->xmapid2);
       }
       if(p->algo_enum == sv_ps && !(p->xmapid1 != p->xmapid2)){
	 printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
		i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	 printf("\t numsv=%d\n", numsv);
	 fflush(stdout);
	 assert(p->xmapid1 != p->xmapid2);
       }
     }
   }

   // Locate Small Duplications based on 2 matchgroups
   int NumSmallDuplications = 0, MaxSmallDuplications = 16;
   CSmallDuplication *SmallDuplications = new CSmallDuplication[MaxSmallDuplications];// will grow dynamically as needed 

   if(SmallDupMaxOverlap > 0.0 && SmallDupMinOverlap > 0.0){/* check for small duplications based on overlap of small matchgroup with insertion outlier in large matchgroup */
     if(VERB/* HERE >=2 */){
       printf("checking for Small Duplications (inverted OR non-inverted): smapsize=%d, numxmapentries=%d:wall time=%0.3f secs\n",smapsize, numxmapentries,wtime());
       fflush(stdout);
     }

     for(int i=0; i < smapsize; i++){
       xmapEntry *pMGi = usexmapentries[i];
       if(VERB>=3){
	 printf("usexmap[i=%d]: qrycontigid=%lld,refcontigid=%lld,or=%d,conf= %0.2f:qrystartpos= %0.3f, qryendpos= %0.3f, refstartpos= %0.3f, refendpos= %0.3f\n",
		i,pMGi->qrycontigid,pMGi->refcontigid, (pMGi->orientforward ? 0 : 1), pMGi->confidence, pMGi->qrystartpos*1e-3, pMGi->qryendpos*1e-3, pMGi->refstartpos*1e-3, pMGi->refendpos*1e-3);
	 fflush(stdout);
       }
       for(int j = i+1; j < smapsize; j++){
	 xmapEntry *pMGj = usexmapentries[j];
	 if(VERB>=3){
	   printf("    usexmap[j=%d]: qrycontigid=%lld,refcontigid=%lld,or=%d,conf= %0.2f:qrystartpos= %0.3f, qryendpos= %0.3f, refstartpos= %0.3f, refendpos= %0.3f\n",
		  j,pMGj->qrycontigid,pMGj->refcontigid, (pMGj->orientforward ? 0 : 1), pMGj->confidence, pMGj->qrystartpos*1e-3, pMGj->qryendpos*1e-3, pMGj->refstartpos*1e-3, pMGj->refendpos*1e-3);
	   fflush(stdout);
	 }

	 if(pMGi->qrycontigid != pMGj->qrycontigid)
	   break;// since qrycontigid are in ascending order
	 if(pMGi->refcontigid != pMGj->refcontigid)
	   continue;

	 double mini = min(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
	 double minj = min(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
	 double maxi = max(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
	 double maxj = max(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
	 double overlap = min(maxi,maxj) - max(mini,minj);
	 double minlen = min(maxi-mini, maxj-minj);
	 double maxlen = max(maxi-mini, maxj-minj);

	 if(overlap <= 0.0)
	   continue;

	 if( overlap < 0.999 * minlen || /* NEW74 */overlap >= 0.99 * maxlen)// Small Duplication should only apply when smaller matchgroup is fully overlapped on qry by larger matchgroup
	   continue;

	 if(pMGi->orientforward == pMGj->orientforward && min(pMGi->confidence,pMGj->confidence) < LogPvThreshold)
	   continue;// NEW88 : Non-Inverted Duplication should satisfy -T threshold for both matchgroups

	 if(VERB>=2){
	   printf("i=%d,j=%d:minlen= %0.3f, overlap= %0.3f, OverlapFilterQryRatio=%0.6f\n", i,j,minlen*1e-3,overlap*1e-3,OverlapFilterQryRatio);
	   fflush(stdout);
	 }

	 /* check reference overlap */
	 double roverlap = min(pMGi->refendpos, pMGj->refendpos) - max(pMGi->refstartpos, pMGj->refstartpos);

	 // HERE HERE : Save the paired inversion breakpoint  ?

	 if(pMGi->confidence < pMGj->confidence){
	   if( overlap < 0.999 * (maxi - mini))// Small Duplication should only apply when smaller matchgroup is fully overlapped on qry by larger matchgroup
	     continue;

	   if(roverlap < 0.999 * (pMGi->refendpos - pMGi->refstartpos))
	     continue; // Small Duplication should only apply when smaller matchgroup is fully overlapped on reference by larger matchgroup

	   if(checkSmallDuplication(SmallDuplications,NumSmallDuplications,MaxSmallDuplications,usexmapentries,j,i)){
	     usexmapentries[i]->output = true;
	     usexmapentries[j]->output = true;
	   }
	 } else {
	   if( overlap < 0.999 * (maxj - minj))// Small Duplication should only apply when smaller matchgroup is fully overlapped on qry by larger matchgroup
	     continue;

	   if(roverlap < 0.999 * (usexmapentries[j]->refendpos - usexmapentries[j]->refstartpos))
	     continue; // Small Duplication should only apply when smaller matchgroup is fully overlapped on reference by larger matchgroup	    

	   if(checkSmallDuplication(SmallDuplications,NumSmallDuplications,MaxSmallDuplications,usexmapentries,i,j)){
	     usexmapentries[i]->output = true;
	     usexmapentries[j]->output = true;
	   }
	 }
       }
     }

     if(VERB/* HERE >=2 */){
       printf("finished checking for small Duplications (inverted OR non-inverted) : %d Duplications found : wall time= %0.6f secs\n", NumSmallDuplications, wtime());
       fflush(stdout);
     }
   }

   if(DEBUG>=2){
     for(int i = 0; i < smapsize;i++){
       xmapEntry *p = usexmapentries[i];

       /* check p->refsites[] (bp) vs refendpos (bp) */
       assert(p->numrefsites >= p->refendidx - p->refstartidx + 1);
       int Istart = p->refstartidx;
       for(int I = Istart; I <= p->refendidx; I++){
	 double sitesI = p->refsites[I - Istart];
	 if(DEBUG && !(p->refstartpos <= sitesI && sitesI <= p->refendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refstartpos= %0.1f,refendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart, p->refsites[I-Istart],p->refstartpos,p->refendpos);
	   fflush(stdout);
	   assert(p->refstartpos <= sitesI && sitesI <= p->refendpos);
	 }
	 if(I > p->refstartidx){
	   if(DEBUG && !(p->refsites[I-1 - Istart] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refsites[%d]= %0.1f, refendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart,p->refsites[I - Istart],I-1-Istart,p->refsites[I-1-Istart],p->refendpos);
	     fflush(stdout);
	     assert(p->refsites[I-1-Istart] <= sitesI);
	   }
	 }
       }
       /* check p->qrysites[0..numqrysites-1] (bp) vs refendpos (bp) */
       double qrystartpos = min(p->qrystartpos,p->qryendpos);
       double qryendpos = max(p->qrystartpos,p->qryendpos);
       int Imin = min(p->qrystartidx,p->qryendidx);
       int Imax = max(p->qrystartidx,p->qryendidx);
       if(DEBUG) assert(p->numqrysites >= Imax-Imin+1);
       for(int I = Imin; I <= Imax; I++){
	 double sitesI = p->qrysites[I - Imin];
	 if(DEBUG && !(qrystartpos <= sitesI && sitesI <= qryendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites, I, I-Imin,p->qrysites[I-Imin],qrystartpos,qryendpos);
	   fflush(stdout);
	   assert(qrystartpos <= sitesI && sitesI <= qryendpos);
	 }
	 if(I > Imin){
	   if(DEBUG && !(p->qrysites[I-1-Imin] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites,I,I-Imin,p->qrysites[I-Imin],I-1-Imin,p->qrysites[I-1-Imin],qrystartpos,qryendpos);
	     fflush(stdout);
	     assert(p->qrysites[I-1-Imin] <= sitesI);
	   }
	 }
       }
     }
   }

   // Locate Small Inversions (see -CutFlip) based on 2 matchgroups
   int NumSmallInversions = 0, MaxSmallInversions = 16;
   CSmallInversion *SmallInversions = new CSmallInversion[MaxSmallInversions];// will grow dynamically as needed 

   if(VERB/* HERE >=2 */){
     printf("checking for Small Inversions: smapsize=%d, numxmapentries=%d: wall time= %0.6f secs\n",smapsize, numxmapentries,wtime());
     fflush(stdout);
   }

   for(int i=0; i < smapsize; i++){
     if(VERB>=3){
       xmapEntry *p = usexmapentries[i];
       printf("usexmap[i=%d]: xmapid=%d: qrycontigid=%lld,refcontigid=%lld,or=%d,conf= %0.2f:qrystartpos= %0.3f, qryendpos= %0.3f, refstartpos= %0.3f, refendpos= %0.3f\n",
	      i,p->xmapid,p->qrycontigid,p->refcontigid, (p->orientforward ? 0 : 1), p->confidence, p->qrystartpos*1e-3, p->qryendpos*1e-3, p->refstartpos*1e-3, p->refendpos*1e-3);
       fflush(stdout);
     }
     for(int j = i+1; j < smapsize; j++){
       if(VERB>=3){
	 xmapEntry *p = usexmapentries[j];
	 printf("    usexmap[j=%d]: qrycontigid=%lld,refcontigid=%lld,or=%d,conf= %0.2f:qrystartpos= %0.3f, qryendpos= %0.3f, refstartpos= %0.3f, refendpos= %0.3f\n",
		j,p->qrycontigid,p->refcontigid, (p->orientforward ? 0 : 1), p->confidence, p->qrystartpos*1e-3, p->qryendpos*1e-3, p->refstartpos*1e-3, p->refendpos*1e-3);
	 fflush(stdout);
       }

       if(usexmapentries[i]->qrycontigid != usexmapentries[j]->qrycontigid)
	 break;// since qrycontigid are in ascending order
       if(usexmapentries[i]->refcontigid != usexmapentries[j]->refcontigid)
	 continue;
       if(usexmapentries[i]->orientforward == usexmapentries[j]->orientforward)
	 continue;

       double mini = min(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
       double minj = min(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
       double maxi = max(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
       double maxj = max(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
       double overlap = min(maxi,maxj) - max(mini,minj);
       double minlen = min(maxi-mini, maxj-minj);
       double maxlen = max(maxi-mini, maxj-minj);

       if(VERB>=2){
	 printf("i=%d,j=%d:minlen= %0.3f, overlap= %0.3f, InvOverlapFilterQryRatio=%0.6f\n", i,j,minlen*1e-3,overlap*1e-3,InvOverlapFilterQryRatio);
	 fflush(stdout);
       }

       if(overlap <= 0.0)
	 continue;

       if( overlap < 0.999 * minlen || /* NEW74 */overlap >= 0.99 * maxlen)// Small Inversion should only apply when smaller matchgroup is fully overlapped on qry by larger matchgroup
	 continue;

       /* check reference overlap */
       double roverlap = min(usexmapentries[i]->refendpos, usexmapentries[j]->refendpos) - max(usexmapentries[i]->refstartpos, usexmapentries[j]->refstartpos);

       if(usexmapentries[i]->confidence < usexmapentries[j]->confidence){
	 if( overlap < 0.999 * (maxi - mini))// Small Inversion should only apply when smaller matchgroup is fully overlapped on qry by larger matchgroup
	   continue;

	 if(roverlap < 0.999 * (usexmapentries[i]->refendpos - usexmapentries[i]->refstartpos))
	   continue; // Small Inversion should only apply when smaller matchgroup is fully overlapped on reference by larger matchgroup

	 if(checkSmallInversion(SmallInversions,NumSmallInversions,MaxSmallInversions,usexmapentries,j,i)){
	   usexmapentries[i]->output = true;
	   usexmapentries[j]->output = true;
	 }
       } else {
	 if( overlap < 0.999 * (maxj - minj))// Small Inversion should only apply when smaller matchgroup is fully overlapped on qry by larger matchgroup
	   continue;

	 if(roverlap < 0.999 * (usexmapentries[j]->refendpos - usexmapentries[j]->refstartpos))
	   continue; // Small Inversion should only apply when smaller matchgroup is fully overlapped on reference by larger matchgroup	    

	 if(checkSmallInversion(SmallInversions,NumSmallInversions,MaxSmallInversions,usexmapentries,i,j)){
	   usexmapentries[i]->output = true;
	   usexmapentries[j]->output = true;
	 }
       }
     }
   }

   if(VERB/* HERE >=2 */){
     printf("finished checking for small Inversions : %d inversions found : wall time= %0.6f secs\n", NumSmallInversions, wtime());
     fflush(stdout);
   }

   if(DEBUG>=2){
     for(int i = 0; i < smapsize;i++){
       xmapEntry *p = usexmapentries[i];

       /* check p->refsites[] (bp) vs refendpos (bp) */
       assert(p->numrefsites >= p->refendidx - p->refstartidx + 1);
       int Istart = p->refstartidx;
       for(int I = Istart; I <= p->refendidx; I++){
	 double sitesI = p->refsites[I - Istart];
	 if(DEBUG && !(p->refstartpos <= sitesI && sitesI <= p->refendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refstartpos= %0.1f,refendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart, p->refsites[I-Istart],p->refstartpos,p->refendpos);
	   fflush(stdout);
	   assert(p->refstartpos <= sitesI && sitesI <= p->refendpos);
	 }
	 if(I > p->refstartidx){
	   if(DEBUG && !(p->refsites[I-1 - Istart] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refsites[%d]= %0.1f, refendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart,p->refsites[I - Istart],I-1-Istart,p->refsites[I-1-Istart],p->refendpos);
	     fflush(stdout);
	     assert(p->refsites[I-1-Istart] <= sitesI);
	   }
	 }
       }
       /* check p->qrysites[0..numqrysites-1] (bp) vs refendpos (bp) */
       double qrystartpos = min(p->qrystartpos,p->qryendpos);
       double qryendpos = max(p->qrystartpos,p->qryendpos);
       int Imin = min(p->qrystartidx,p->qryendidx);
       int Imax = max(p->qrystartidx,p->qryendidx);
       if(DEBUG) assert(p->numqrysites >= Imax-Imin+1);
       for(int I = Imin; I <= Imax; I++){
	 double sitesI = p->qrysites[I - Imin];
	 if(DEBUG && !(qrystartpos <= sitesI && sitesI <= qryendpos)){
	   printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		  i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites, I, I-Imin,p->qrysites[I-Imin],qrystartpos,qryendpos);
	   fflush(stdout);
	   assert(qrystartpos <= sitesI && sitesI <= qryendpos);
	 }
	 if(I > Imin){
	   if(DEBUG && !(p->qrysites[I-1-Imin] <= sitesI)){
	     printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		    i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites,I,I-Imin,p->qrysites[I-Imin],I-1-Imin,p->qrysites[I-1-Imin],qrystartpos,qryendpos);
	     fflush(stdout);
	     assert(p->qrysites[I-1-Imin] <= sitesI);
	   }
	 }
       }
     }
   }

   printf("sv array: numsv= %d, smallInversions= %d maxsv= %d: wall time= %0.6f secs\n", numsv, NumSmallInversions, maxsv,wtime()); fflush(stdout); //debug
   if( numsv + NumSmallInversions*2 > maxsv ) //if not, no need to increase size
     maxsv = growArray(sv_array, maxsv, max(numsv * 3 / 2 , numsv + NumSmallInversions*2)); //newsize is third arg  

   if(NumSmallInversions > 0){    // append small inversions in SmallInversions[0 .. NumSmallInversions - 1] to sv_array[] : 2 entries per inversion
     for(int t = 0; t < NumSmallInversions; t++){
       CSmallInversion *p = &SmallInversions[t];
       int Large = findXmapEntry(p->Large,xmapentries,numxmapentries);
      if(DEBUG && !(0 <= Large && Large < numxmapentries && xmapentries[Large] == p->Large)){
	printf("SmallInversions[%d]: findXmapEntry(p->Large = %p) returned Large = %d (numxmapentries=%d)\n",
	       t, p->Large, Large, numxmapentries);
	fflush(stdout);
	assert(0 <= Large && Large < numxmapentries && xmapentries[Large] == p->Large);
      }
      int Small = findXmapEntry(p->Small,xmapentries,numxmapentries);
      if(DEBUG && !(0 <= Small && Small < numxmapentries && xmapentries[Small] == p->Small)){
	printf("SmallInversions[%d]: findXmapEntry(p->Small = %p) returned Small = %d (numxmapentries=%d)\n",
	       t, p->Small, Small, numxmapentries);
	fflush(stdout);
	assert(0 <= Small && Small < numxmapentries && xmapentries[Small] == p->Small);
      }

      FLOAT *Y = p->Y;
      FLOAT *X = p->X;

      p->Small->smallDupInv = true;// mark smaller matchgroup so it is not use to call an intra-chromosomal translocation
      if(MATCHGROUP_TRIM && XmapUnique && !(p->Small->xmapid_overlapped > 0)){
	printf("SmallInversion[%d]: refid=%lld,qryid=%lld,xmapid1=%d,xmapid2=%d (Small=%d,p->Small->xmapid_overlapped=%d)\n",
	       t,p->refcontigid,p->qrycontigid,p->Large->xmapid,p->Small->xmapid,Small,p->Small->xmapid_overlapped);

	if(VERB>=2){
	  for(int i = 0; i < numxmapentries; i++){
	    xmapEntry *p = usexmapentries[i];
	    if(p->xmapid == 9337 || p->xmapid == 1385){
	      printf("\t usexmapentries[%d]: refid=%lld,qryid=%lld, xmapid=%d, xmapid_overlapped=%d\n",
		     i, p->refcontigid,p->qrycontigid,p->xmapid,p->xmapid_overlapped);
	      fflush(stdout);
	    }
	  }
	}

	fflush(stdout);
	assert(p->Small->xmapid_overlapped > 0);
      }

      if(VERB>=2){
	printf("sv_array[%d,%d] : added small inversion pair with xmapid1=%d,xmapid2=%d\n",numsv,numsv+1,p->Large->xmapid,p->Small->xmapid);
	fflush(stdout);
      }

      //note that the algo_enum is sv_indel because the xmapids used (Small & Large) may no longer be in usexmapentries (even though the actual algo is RefSplit)
      //and *1e3 is for units: X/Y are in kb

      /* first breakpoint is the "Start"(left) side of inversion */
      sv_array[numsv] = structuralVariation(Large,Small,0.0,0.0,X[p->QryOutlierStart]*1e3,X[p->QryInversionStart]*1e3,
					    p->refcontigid,p->refcontigid,Y[p->RefOutlierStart]*1e3,Y[p->RefInversionStart]*1e3, "small_inversion",
					    inversion_small_type, p->QryOutlierStart, p->QryInversionStart,
					    p->RefOutlierStart,p->RefInversionStart,-1.0,p->LeftConf,p->Small->confidence, sv_indel, p->Large->xmapid, p->Small->xmapid);
      //if(DEBUG) assert(xmapentries[Large]->qrycontigid == xmapentries[Small]->qrycontigid);

      /* second breakpoint is the "End"(right) side of inversion */
      sv_array[numsv+1] = structuralVariation(Large,Small,0.0,0.0,X[p->QryInversionEnd]*1e3,X[p->QryOutlierEnd]*1e3,
					      p->refcontigid,p->refcontigid,Y[p->RefInversionEnd]*1e3,Y[p->RefOutlierEnd]*1e3, "small_inversion",
					      inversion_small_type, p->QryInversionEnd, p->QryOutlierEnd,
					      p->RefInversionEnd,p->RefOutlierEnd,-1.0, p->Small->confidence, p->RightConf, sv_indel, p->Large->xmapid, p->Small->xmapid);

      sv_array[numsv+1].refstart2 = sv_array[numsv].refstart;
      sv_array[numsv+1].refstop2 = sv_array[numsv].refstop;

      sv_array[numsv].refstart2 = sv_array[numsv+1].refstart;
      sv_array[numsv].refstop2 = sv_array[numsv+1].refstop;

      numsv += 2;
    }
  }

  delete [] SmallInversions; SmallInversions = NULL;

  int numsvIndel2 = numsv;/* original number of single Matchgroup Indels + small Inversions */

  if(DEBUG>=2){/* check if any 1-MG Indels no longer overlap the MGs they refer to */
    for(int i = 0; i < numsvIndel2; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum != sv_ps);
      if(DEBUG && i < numsvIndel1) assert(p->indices[0] == p->indices[1]);

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;
      
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      if(DEBUG) assert(entry1->output == true);
      if(entry1->output == false)
	continue;

      double qmin = min(entry1->qrystartpos,entry1->qryendpos);
      double qmax = max(entry1->qrystartpos,entry1->qryendpos);
      int iqmin = min(entry1->qrystartidx,entry1->qryendidx);
      int iqmax = max(entry1->qrystartidx,entry1->qryendidx);

      if(!((qmin <= p->querystart && p->querystart <= qmax) &&
	   (qmin <= p->querystop && p->querystop <= qmax) &&
	   (iqmin <= p->querystartidx && p->querystartidx <= iqmax) &&
	   (iqmin <= p->querystopidx && p->querystopidx <= iqmax) &&
	   (entry1->refstartpos <= p->refstart && p->refstart <= entry1->refendpos) &&
	   (entry1->refstartpos <= p->refstop && p->refstop <= entry1->refendpos) &&
	   (entry1->refstartidx <= p->refstartidx && p->refstartidx <= entry1->refendidx) &&
	   (entry1->refstartidx <= p->refstopidx && p->refstopidx <= entry1->refendidx))){
	if(i < numsvIndel1)
	  printf("sv_array[%d]:Single MG indel's MG was trimmed and no longer includes the SV coordinates:\n",i);
	else 
	  printf("sv_array[%d]: 2-MG small inversion's large MG was trimmed and no longer includes the SV coordinates:\n",i);
	printf("sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	printf("numsvIndel1=%d,numsvIndel2=%d\n",numsvIndel1,numsvIndel2);
	fflush(stdout);
	assert(0);
      }
      if(i < numsvIndel1) assert(p->type_enum == insertion || p->type_enum == deletion);
      else assert(p->type_enum == inversion_small_type);
    }
  }

  if(DEBUG>=2){
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsv=%d\n", numsv);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
      if(p->algo_enum == sv_ps && !(p->xmapid1 != p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsv=%d\n", numsv);
	fflush(stdout);
	assert(p->xmapid1 != p->xmapid2);
      }
    }
  }

  printf("sv array: numsv= %d, smallDuplications= %d maxsv= %d: wall time= %0.6f secs\n", numsv, NumSmallDuplications, maxsv, wtime()); fflush(stdout); //debug
  if( numsv + NumSmallDuplications > maxsv ) //if not, no need to increase size
    maxsv = growArray(sv_array, maxsv, max(numsv * 3 / 2, numsv + NumSmallDuplications)); //newsize is third arg  

  if(NumSmallDuplications > 0){ //  append small duplications in SmallDuplications[0 .. NumSmallDuplications - 1] to sv_array[] : 1 entry per duplication
    for(int t = 0; t < NumSmallDuplications; t++){
      CSmallDuplication *p = &SmallDuplications[t];
      int Large = findXmapEntry(p->Large,xmapentries,numxmapentries);
      if(DEBUG && !(0 <= Large && Large < numxmapentries && xmapentries[Large] == p->Large)){
	printf("SmallDuplications[%d]: findXmapEntry(p->Large = %p) returned Large = %d (numxmapentries=%d)\n",
	       t, p->Large, Large, numxmapentries);
	fflush(stdout);
	assert(0 <= Large && Large < numxmapentries && xmapentries[Large] == p->Large);
      }
      int Small = findXmapEntry(p->Small,xmapentries,numxmapentries);
      if(DEBUG && !(0 <= Small && Small < numxmapentries && xmapentries[Small] == p->Small)){
	printf("SmallDuplications[%d]: findXmapEntry(p->Small = %p) returned Small = %d (numxmapentries=%d)\n",
	       t, p->Small, Small, numxmapentries);
	fflush(stdout);
	assert(0 <= Small && Small < numxmapentries && xmapentries[Small] == p->Small);
      }

      FLOAT *Y = p->Y;
      FLOAT *X = p->X;

      p->Small->smallDupInv = true;// mark smaller matchgroup so it is not use to call an intra-chromosomal translocation
      if(MATCHGROUP_TRIM && XmapUnique) assert(p->Small->xmapid_overlapped > 0);

      //note that the algo_enum is sv_indel because the xmapids used (Small & Large) may no longer be in usexmapentries (even though the actual algo is RefSplit)
      //and *1e3 is for units: X/Y are in kb

      sv_array[numsv] = structuralVariation(Large, Small, 0.0, 0.0, X[p->QryDupStart]*1e3, X[p->QryDupEnd]*1e3, 
					    p->refcontigid, p->refcontigid, Y[p->RefDupStart]*1e3, Y[p->RefDupEnd]*1e3, 
					    p->inverted ? "inverted_duplication" : "duplication", p->inverted ? duplicationinv_type : duplication_type, 
					    p->QryDupStart, p->QryDupEnd, p->RefDupStart, p->RefDupEnd, p->confBC, p->LeftConf, p->RightConf, sv_indel, p->Large->xmapid, p->Small->xmapid);
      //if(DEBUG) assert(xmapentries[Large]->qrycontigid == xmapentries[Small]->qrycontigid);
      numsv++;
      if(DEBUG) assert(numsv <= maxsv);

      // HERE HERE : also output paired inversion breakpoint
    }
  }

  delete [] SmallDuplications; SmallDuplications = NULL;

  if(DEBUG>=2){
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsv=%d\n", numsv);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
      if(p->algo_enum == sv_ps && !(p->xmapid1 != p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsv=%d\n", numsv);
	fflush(stdout);
	assert(p->xmapid1 != p->xmapid2);
      }

    }
  }

  if(DEBUG>=2){
    for(int i = 0; i < smapsize;i++){
      xmapEntry *p = usexmapentries[i];

      /* check p->refsites[] (bp) vs refendpos (bp) */
      assert(p->numrefsites >= p->refendidx - p->refstartidx + 1);
      int Istart = p->refstartidx;
      for(int I = Istart; I <= p->refendidx; I++){
	double sitesI = p->refsites[I - Istart];
	if(DEBUG && !(p->refstartpos <= sitesI && sitesI <= p->refendpos)){
	  printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refstartpos= %0.1f,refendpos= %0.1f\n",
		 i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart, p->refsites[I-Istart],p->refstartpos,p->refendpos);
	  fflush(stdout);
	  assert(p->refstartpos <= sitesI && sitesI <= p->refendpos);
	}
	if(I > p->refstartidx){
	  if(DEBUG && !(p->refsites[I-1 - Istart] <= sitesI)){
	    printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refsites[%d]= %0.1f, refendpos= %0.1f\n",
		   i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart,p->refsites[I - Istart],I-1-Istart,p->refsites[I-1-Istart],p->refendpos);
	    fflush(stdout);
	    assert(p->refsites[I-1-Istart] <= sitesI);
	  }
	}
      }
      /* check p->qrysites[0..numqrysites-1] (bp) vs refendpos (bp) */
      double qrystartpos = min(p->qrystartpos,p->qryendpos);
      double qryendpos = max(p->qrystartpos,p->qryendpos);
      int Imin = min(p->qrystartidx,p->qryendidx);
      int Imax = max(p->qrystartidx,p->qryendidx);
      if(DEBUG) assert(p->numqrysites >= Imax-Imin+1);
      for(int I = Imin; I <= Imax; I++){
	double sitesI = p->qrysites[I - Imin];
	if(DEBUG && !(qrystartpos <= sitesI && sitesI <= qryendpos)){
	  printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		 i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites, I, I-Imin,p->qrysites[I-Imin],qrystartpos,qryendpos);
	  fflush(stdout);
	  assert(qrystartpos <= sitesI && sitesI <= qryendpos);
	}
	if(I > Imin){
	  if(DEBUG && !(p->qrysites[I-1-Imin] <= sitesI)){
	    printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		   i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites,I,I-Imin,p->qrysites[I-Imin],I-1-Imin,p->qrysites[I-1-Imin],qrystartpos,qryendpos);
	    fflush(stdout);
	    assert(p->qrysites[I-1-Imin] <= sitesI);
	  }
	}
      }
    }
  }

  int numsvIndel = numsv;/* original number of SVs with Indels with single matchgroup created in output_xmap() or small Inversion or Duplication (all with algo_enum == sv_indel) */

  if(DEBUG>=2){/* check if any 1-MG Indels no longer overlap the MGs they refer to */
    for(int i = 0; i < numsvIndel; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum != sv_ps);
      if(DEBUG && i < numsvIndel1) assert(p->indices[0] == p->indices[1]);

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;
      
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      if(DEBUG) assert(entry1->output == true);
      if(entry1->output == false)
	continue;

      double qmin = min(entry1->qrystartpos,entry1->qryendpos);
      double qmax = max(entry1->qrystartpos,entry1->qryendpos);
      int iqmin = min(entry1->qrystartidx,entry1->qryendidx);
      int iqmax = max(entry1->qrystartidx,entry1->qryendidx);

      if(!((qmin <= p->querystart && p->querystart <= qmax) &&
	   (qmin <= p->querystop && p->querystop <= qmax) &&
	   (iqmin <= p->querystartidx && p->querystartidx <= iqmax) &&
	   (iqmin <= p->querystopidx && p->querystopidx <= iqmax) &&
	   (entry1->refstartpos <= p->refstart && p->refstart <= entry1->refendpos) &&
	   (entry1->refstartpos <= p->refstop && p->refstop <= entry1->refendpos) &&
	   (entry1->refstartidx <= p->refstartidx && p->refstartidx <= entry1->refendidx) &&
	   (entry1->refstartidx <= p->refstopidx && p->refstopidx <= entry1->refendidx))){
	if(i < numsvIndel1)
	  printf("sv_array[%d]:Single MG indel's MG was trimmed and no longer includes the SV coordinates:\n",i);
	else if(i < numsvIndel2)
	  printf("sv_array[%d]: 2-MG small inversion's large MG was trimmed and no longer includes the SV coordinates:\n",i);
	else
	  printf("sv_array[%d]: 2-MG small duplication's large MG was trimmed and no longer includes the SV coordinates:\n",i);
	printf("sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	printf("numsvIndel1=%d,numsvIndel2=%d,numsvIndel=%d\n",numsvIndel1,numsvIndel2,numsvIndel);
	fflush(stdout);
	assert(0);
      }
      if(i < numsvIndel1 ? !(p->type_enum == insertion || p->type_enum == deletion) :
	 i < numsvIndel2 ? !(p->type_enum == inversion_small_type) :
	 i < numsvIndel ? !(p->type_enum == duplication_type || p->type_enum == duplicationinv_type) : 0){
	printf("sv_array[%d]: type_enum=%d, numsvIndel1=%d,numsvIndel2=%d, numsvIndel=%d\n",i,p->type_enum,numsvIndel1,numsvIndel2,numsvIndel);
	printf("\t sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("\t xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	fflush(stdout);
	assert(0);
      }
    }
  }

  // apply xmap filter based on  -svSideMinLen && -svAlignConfByLen (Also filter out matchgroups that don't satisfy the original -S -T thresholds)
  // NOTE : for RefSplit matchgroups the full unsplit matchgroup confidence (logPV2) is used
  // NOTE : for CutFlip matchgroups the confidence logPV2 is raised to at least the -T2 threshold.
  //        Mark the MGs to be filtered due to -S or -T by inversion_only, so Endoutlier CutFlip matchgroups can still be used
  int nsmapsize = 0;
  for(register int i=0; i < smapsize; i++){

    Calign *p = usexmapentries[i]->align;
    if(RefSplit && (p->score <= ScoreThreshold || p->logPV2 <= LogPvThreshold)){
      if(verbose){
    	printf("Xmapid %2i (rid=%lld,qid=%lld,fwd=%d) fails -S %0.2f -T %0.2f thresholds : score= %0.6f logPV= %0.2f, logPV2= %0.2f : marking as inversion_only\n", 
	       usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid,  usexmapentries[i]->orientforward ? 1 : 0, ScoreThreshold, LogPvThreshold, p->score, p->logPV, p->logPV2);
	fflush(stdout);
      }
      usexmapentries[i]->inversion_only = true;
    }

    // require aligned regions are at least smap_side_alignkb (commandline option -svSideMinLen)
    // this requirement should be applied to each alignment, not to each potential SV.
    // the latter can potentially prevent a valid sv from being found.
    // example: alignments 1 and 3 form an SV, but alignment 2 is too small
    //  when considering 1-2 and 2-3, both are discarded the old way, and 1-3 is never considered.
    if( fabs(usexmapentries[i]->qrystartpos - usexmapentries[i]->qryendpos)/1000. < smap_side_alignkb || //current query
    	fabs(usexmapentries[i]->refstartpos - usexmapentries[i]->refendpos)/1000. < smap_side_alignkb ) { //current ref
      if( verbose ){
    	printf("xmapid %2i fails align len: qryid=%lld,refid=%lld, qry= %f ... %f, ref= %f ... %f\n", usexmapentries[i]->xmapid, usexmapentries[i]->qrycontigid, usexmapentries[i]->refcontigid,
	       usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos, usexmapentries[i]->refstartpos, usexmapentries[i]->refendpos );
	fflush(stdout);
      }
      continue;
    }
    
    //minalign is min of qry and ref--these are in bases, convert to kb
    double minalign = min( fabs(usexmapentries[i]->qrystartpos - usexmapentries[i]->qryendpos),
			   fabs(usexmapentries[i]->refstartpos - usexmapentries[i]->refendpos) ) / 1000.; 

    /* smap: minimum confidence (in xmap) over length (kb) of aligned region -- from parameters.h/cpp */
    //float smap_confbylen = 0.07; 
    if( usexmapentries[i]->confidence / minalign < smap_confbylen ) {
      if( verbose ){
	printf("xmapid %2i fails conf by len: %f, %f, %f\n", usexmapentries[i]->xmapid, usexmapentries[i]->confidence, minalign, usexmapentries[i]->confidence/minalign );
	fflush(stdout);
      }
      continue;
    }

    usexmapentries[nsmapsize] = usexmapentries[i]; //no need to copy, just assign pointers
    usexmapentries[nsmapsize]->callRepeat(false, REPEAT_TOLERANCE,REPEAT_MINELE,REPEAT_CONF_PENALTY); //check this entry for repeats: first arg is verbose: this is too verbose: better to report excluded SVs rather than every single repeat (default false)
    nsmapsize++; //different from xmapEntry.xmapid bc of above continue
  } //end for numxmapentries

  if( VERB ){
    printf("Xmap entries removed= %d, Surviving xmap entries= %d (wall time= %0.6f)\n", smapsize - nsmapsize, nsmapsize, wtime());
    fflush(stdout);
  }
  smapsize = nsmapsize;

  if( smapsize == 0 )
    verbose = false; //no point in printing anything else

  if(VERB>=2){
    for(int i = 0; i < smapsize; i++){
      xmapEntry *p = usexmapentries[i];
      printf("usexmapentries[i=%d]: xmapid=%d, rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inversion_only=%d\n",
	     i,p->xmapid,p->refcontigid,p->qrycontigid,p->orientforward ? 1 : 0, p->confidence,
	     p->inversion_only ? 1 : 0);
    }
    fflush(stdout);
  }

  if(VERB){
    printf("Applying matchgroup overlap filter : Query Overlap ratio must be %0.4f (in kb) AND %0.4f (in labels) AND confidence delta must be %0.2f\n",
	   OverlapFilterQryRatio, OverlapFilterQryLabRatio, OverlapFilterConfDelta);
    printf("\t Any opposite orientation matchgroup is saved for inversions only, unless Query Overlap ratios exceed %0.4f (in kb) AND %0.4f (in labels) AND conf delta exceeds %0.2f\n",
	   max(OverlapFilterQryRatio,InvOverlapFilterQryRatio), max(OverlapFilterQryLabRatio, InvOverlapFilterQryLabRatio), max(OverlapFilterConfDelta, InvOverlapFilterConfDelta));
    printf("\t Qry Label overlap requirement is waived for MG pairs that have already generated a Small Inversion or Inverted Duplication call : numsvIndel1=%d,numsvIndel2=%d,numsvIndel=%d\n",
	   numsvIndel1,numsvIndel2,numsvIndel);
    if(INV_OVERLAP_FILTER_FIX)
      printf("\t Qry Label overlap requirement is also waived for inverted MG pairs where smaller MG is fully overlapped by larger MG in both ref and qry\n");
    if(OverlapFilterSameRef)
      printf("\t Filter is disabled if overlapped MGs have different ref contigs (-OverlapFilterSameRef 1)\n");
    fflush(stdout);
  }

  qsort(&sv_array[numsvIndel1], numsvIndel-numsvIndel1, sizeof(*sv_array), compare_svs);

  if(VERB>=2){
    printf("Small Inversions and Duplications:\n");
    for(int s = numsvIndel1;  s < numsvIndel; s++){
      structuralVariation *p = &sv_array[s];
      printf("sv_array[%d]: xmapid1=%d,xmapid2=%d, type_enum=%d\n", s,p->xmapid1,p->xmapid2,p->type_enum);
    }
    fflush(stdout);
  }

  int filtercnt = OverlapFilterMG(usexmapentries, smapsize, numsvIndel1, numsvIndel, false, verbose);

  // filter usexmapentries based on sv_overlap
  if( OverlapFilterQryRatio > 0.0 ) {
    int newsmapsize=0;
    for( int i=0; i < smapsize; i++) {
      if(usexmapentries[i]->sv_overlap)
	continue;
      if(usexmapentries[i]->largeMG_only) continue;// HERE HERE : remove after implementing largeMG_only restriction
      usexmapentries[newsmapsize++] = usexmapentries[i];
    }
    if( VERB ){
      printf("%d xmap entries removed due to overlap: %d xmap entries remaining (time = %0.6f secs)\n", smapsize-newsmapsize, newsmapsize, wtime());
      fflush(stdout);
    }
    // HERE HERE : uncomment after implementing largeMG_only restriction    if(DEBUG) assert(filtercnt == smapsize-newsmapsize);
    smapsize = newsmapsize; //this is len of usexmapentries
  } //end if dosv_simpleoverlap for filtering usexmapentries

  if(VERB>=2){
    printf("After OverlapFilter:smapsize=%d\n",smapsize);
    for(int i = 0; i < smapsize; i++){
      xmapEntry *p = usexmapentries[i];
      printf("usexmapentries[%d]: xmapid=%d, rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inversion_only=%d\n",
	     i,p->xmapid,p->refcontigid,p->qrycontigid,p->orientforward ? 1 : 0, p->confidence,
	     p->inversion_only ? 1 : 0);
    }
    fflush(stdout);
  }

  if(MergeMG) {
    for(int i = 0; i < smapsize; i++){
      usexmapentries[i]->orig = NULL;
      if(DEBUG) usexmapentries[i]->jMerge = -1;
    }

    int filtercnt = 0;/* number of MGs marked for deletion with sv_overlap = true */

    while(1){/* repeat MG merging and filtering until there is no more progress */
      int origfiltercnt = filtercnt;

      if(VERB/* HERE HERE >=2 */){
	printf("Trying to temporarily merge matchgroups with gaps <= %0.2f kb and no other matchgroup overlapping the gap in qry or ref: wall time= %0.6f secs\n", MergeMG,wtime());
	fflush(stdout);
      }

      /* merge non-overlapping matchgroup pairs with same orientation and same query and refid that are seperated by the smallest qry+ref gap (in labels) for either matchgroup 
	 (NOTE : No need to save indel implied by gap, since it should later be called as 2 matchgroup indel)
	 This allows additional superfluous matchgroups to be eliminated in segdup regions and acts as a very conservative parsimony filter */
      for(int i = 0; i < smapsize - 1; i++){
	if( usexmapentries[i]->sv_overlap )
	  continue;
	for(int j = i+1; j < smapsize; j++) {
	  if(usexmapentries[i]->qrycontigid != usexmapentries[j]->qrycontigid)
	    break;
	  if(usexmapentries[i]->refcontigid != usexmapentries[j]->refcontigid)
	    break;
	  if( usexmapentries[j]->sv_overlap )
	    continue;
	  if( usexmapentries[i]->orientforward != usexmapentries[j]->orientforward) // need same orientation
	    continue;

	  /* try to merge MG j with MG i */

	  if(DEBUG>=2) assert(usexmapentries[i]->refstartpos <= usexmapentries[j]->refstartpos);
	  double refoverlap = usexmapentries[i]->refendpos - usexmapentries[j]->refstartpos;	
	  //	  double mini = min(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
	  //	  double minj = min(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
	  //	  double maxi = max(usexmapentries[i]->qrystartpos, usexmapentries[i]->qryendpos);
	  //	  double maxj = max(usexmapentries[j]->qrystartpos, usexmapentries[j]->qryendpos);
	  double overlap = usexmapentries[i]->orientforward ? usexmapentries[i]->qryendpos - usexmapentries[j]->qrystartpos :  
	    usexmapentries[j]->qrystartpos - usexmapentries[i]->qryendpos; // NOT min(maxi,maxj) - max(mini,minj);
	
	  if(VERB>=2){
	    printf("Checking xid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f) for merge with xid %i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f): overlap= %0.3f, refoverlap= %0.3f kb\n", 
		   usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
		   usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
		   overlap*1e-3,refoverlap*1e-3);
	    fflush(stdout);
	  }
	
	  if(overlap >= 0.0 || refoverlap >= 0.0)
	    continue;

	  /* verify that no other MG k  = j + 1 .. end - 1 is closer : this should always be the case on reference but may not be true for query (if MG j and MG k are "crossed) */
	  double gap = -(overlap + refoverlap);/* ref + qry gap size as measure of "closeness" : could use number of labels instead */

	  for(int k = j+1; k < smapsize; k++){
	    if(usexmapentries[i]->qrycontigid != usexmapentries[k]->qrycontigid)
	      break;
	    if(usexmapentries[i]->refcontigid != usexmapentries[k]->refcontigid)
	      break;
	    if( usexmapentries[k]->sv_overlap ) //already excluded (by previous i), no need to recheck
	      continue;
	    if( usexmapentries[i]->orientforward != usexmapentries[k]->orientforward) // need same orientation
	      continue;

	    /* verify the ref gap is larger for MG k than for MG j (vs MG i) */
	    if(DEBUG>=2) assert(usexmapentries[i]->refstartpos <= usexmapentries[k]->refstartpos);
	    double refoverlapk = usexmapentries[i]->refendpos - usexmapentries[k]->refstartpos;
	    if(DEBUG)  assert(refoverlapk <= refoverlap);

	    //	    double mink = min(usexmapentries[k]->qrystartpos, usexmapentries[k]->qryendpos);
	    //	    double maxk = max(usexmapentries[k]->qrystartpos, usexmapentries[k]->qryendpos);
	    double overlapk = usexmapentries[i]->orientforward ? usexmapentries[i]->qryendpos - usexmapentries[k]->qrystartpos :  
	      usexmapentries[k]->qrystartpos - usexmapentries[i]->qryendpos; // NOT min(maxi,maxk) - max(mini,mink);

	    if(overlapk >= 0.0)
	      continue;

	    if(VERB>=2 && -(overlapk + refoverlapk) < gap){
	      printf("Checking xid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f) for merge with xid %i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f): overlap= %0.3f, refoverlap= %0.3f kb (best xid=%i)\n", 
		     usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
		     usexmapentries[k]->xmapid, usexmapentries[k]->refcontigid, usexmapentries[k]->qrycontigid, usexmapentries[k]->orientforward ? 1 : 0, usexmapentries[k]->confidence,
		     overlapk*1e-3,refoverlapk*1e-3, usexmapentries[j]->xmapid);
	      fflush(stdout);
	    }

	    if(-refoverlapk >= gap || -refoverlapk > MergeMG)
	      break;/* no possibility of find a closer MG with gap sizes <= MergeMG in usexmapentries[k .. end-1], since MGs are sorted in ascending order of reference location */

	    double gapk = -(overlapk + refoverlapk);

	    if(gapk < gap){
	      j = k;
	      gap = gapk;

	      //	      minj = mink;
	      //	      maxj = maxk;
	      overlap = overlapk;
	      refoverlap = refoverlapk;
	    }
	  }

	  if(max(-overlap,-refoverlap) * 1e-3 > MergeMG)
	    continue;

	  if(VERB>=2){
	    printf("Merging xid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f) with xid %i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f): qrygap= %0.3f, refgap= %0.3f kb (MaxGap= %0.3f)\n", 
		   usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence,
		   usexmapentries[j]->xmapid, usexmapentries[j]->refcontigid, usexmapentries[j]->qrycontigid, usexmapentries[j]->orientforward ? 1 : 0, usexmapentries[j]->confidence,
		   -overlap*1e-3,-refoverlap*1e-3,MergeMG);
	    fflush(stdout);
	  }

	  /* create temp copy of MG i so it can be restored later */
	  xmapEntry *pxmap = new xmapEntry(*usexmapentries[i]);// NOTE: only pointers to refsites[] and qrysites[] are copied by default copy constructor
	  pxmap->jMerge = j;/* remember which matchgroup was merged into MG i */
	  pxmap->orig = usexmapentries[i];
	  usexmapentries[i] = pxmap;

	  /* temporarily merge MG j into MG i */ 
	  usexmapentries[i]->refendpos = usexmapentries[j]->refendpos;
	  usexmapentries[i]->refendidx = usexmapentries[j]->refendidx;
	  usexmapentries[i]->qryendpos = usexmapentries[j]->qryendpos;
	  usexmapentries[i]->qryendidx = usexmapentries[j]->qryendidx;
	  usexmapentries[i]->confidence += usexmapentries[j]->confidence;// approximation
	  usexmapentries[i]->numrefsites = usexmapentries[i]->refendidx - usexmapentries[i]->refstartidx + 1;
	  usexmapentries[i]->refsites = NULL;/* refsites[] & qrysites[] will be restored later from pxmap->orig */
	  usexmapentries[i]->qrysites = NULL;
	
	  /* compute merged alignment */
	  Calign *p = usexmapentries[i]->align;
	  Calign *q = usexmapentries[j]->align;
	  int U1 = p->numpairs;
	  int U2 = q->numpairs;
	  Cmap *Ymap = YYmap[p->mapid1];
	  Cmap *Xmap = XXmap[p->mapid2];
	  FLOAT *Y = Ymap->site[0];
	  FLOAT *X = Xmap->site[0];
	  int N = Ymap->numsite[0];
	  int M = Xmap->numsite[0];

	  Calign *r = new Calign;
	  r->numpairs = U1+U2;
	  r->expand_arrays(r->numpairs);
	  r->mapid1 = p->mapid1;
	  r->mapid2 = p->mapid2;
	  r->orientation = p->orientation;
	  r->rev = p->rev;
	  r->scaleID = p->scaleID;
	  r->repeat = 0;
	  r->maxoutlier = max(p->maxoutlier, q->maxoutlier);
	  r->maxoutlier = max(r->maxoutlier, fabs(overlap-refoverlap));
	  r->noutliers = p->noutliers + q->noutliers + 1;
	  r->maxoutlierLabels = max(p->maxoutlierLabels, q->maxoutlierLabels);
	  // HERE HERE  r->maxoutlierLabels = max(r->maxoutlierLabels, I-K-lastK + J-lastJ - 2);
	  
	  if(CutFlip){
	    if(DEBUG) assert(1 <= p->Nrange && p->Nrange <= N && 1 <= p->Mrange && p->Mrange <= M);
	    if(DEBUG) assert(1 <= q->Nrange && q->Nrange <= N && 1 <= q->Mrange && q->Mrange <= M);
	    r->Nrange = min(p->Nrange + q->Nrange + q->sites1[0] - p->sites1[U1-1], N);/* upper bound */
	    r->Mrange = min(p->Mrange + q->Mrange + abs(q->sites2[0] - p->sites2[U1-1]), M);/* upper bound */
	  } else {
	    r->Nrange = N;
	    r->Mrange = M;
	  }

	  /* left end */
	  r->Lend = p->Lend;
	  r->Lij1 = p->Lij1;
	  r->Lij2 = p->Lij2;
	  
	  /* copy p */
	  for(int t = 0; t < U1; t++){
	    r->iscore[t] = p->iscore[t];
	    r->outscore[t] = p->outscore[t];
	    r->sites1[t] = p->sites1[t];
	    r->sitesK1[t] = p->sitesK1[t];
	    r->sites2[t] = p->sites2[t];
	  }

	  /* compute gap scores */
	  double Ygap = Yc(Y,q->sites1[0],q->sitesK1[0]) - Yc(Y,p->sites1[U1-1],p->sitesK1[U1-1]);
	  double Xgap = r->orientation ? X[M+1-p->sites2[U1-1]] - X[M+1-q->sites2[0]] : X[q->sites2[0]] - X[p->sites2[U1-1]];
	  // NOTE : Ygap and Xgap may differ from -refoverlap and -overlap if K values are > 0, since xmapEntry ignores K values
	  if(DEBUG && !q->sitesK1[0] && !p->sitesK1[U1-1] && !(fabs(Ygap + refoverlap*1e-3) <= 1e-6)){
	    printf("q->sites1[0]=IR=%d, q->sitesK1[0]=KR=%d,p->sites1[U1-1]=IL=%d, p->sitesK1[U1-1]=KL=%d: Yc(Y,IR,KR)= %0.4f, Yc(Y,IL,KL)= %0.4f, refoverlap= %0.6f, Ygap= %0.6f\n",
		   q->sites1[0],q->sitesK1[0],p->sites1[U1-1], p->sitesK1[U1-1], Yc(Y,q->sites1[0],q->sitesK1[0]),Yc(Y,p->sites1[U1-1],p->sitesK1[U1-1]),refoverlap*1e-3, Ygap);
	    fflush(stdout);
	    assert(fabs(Ygap + refoverlap*1e-3) <= 1e-6);
	  }
	  if(DEBUG && !q->sitesK1[0] && !p->sitesK1[U1-1] && !(fabs(Xgap + overlap*1e-3) <= 1e-6)){
	    printf("q->sites1[0]=IR=%d, q->sitesK1[0]=KR=%d,p->sites1[U1-1]=IL=%d, p->sitesK1[U1-1]=KL=%d: Yc(Y,IR,KR)= %0.4f, Yc(Y,IL,KL)= %0.4f, refoverlap= %0.6f, Ygap= %0.6f\n",
		   q->sites1[0],q->sitesK1[0],p->sites1[U1-1], p->sitesK1[U1-1], Yc(Y,q->sites1[0],q->sitesK1[0]),Yc(Y,p->sites1[U1-1],p->sitesK1[U1-1]),refoverlap*1e-3,Ygap);
	    printf("q->sites2[0]=JR=%d,p->sites2[U1-1]=JL=%d,or=%d, M=%d: X[JL]= %0.6f, X[JR]= %0.6f : Xgap= %0.6f, -overlap= %0.6f\n",
		   q->sites2[0],p->sites2[U1-1],r->orientation,M,
		   r->orientation ? X[M+1]-X[M+1-p->sites2[U1-1]] : X[p->sites2[U1-1]],
		   r->orientation ? X[M+1]-X[M+1-q->sites2[0]] : X[q->sites2[0]],
		   Xgap,overlap*1e-3);
	    fflush(stdout);
	    assert(fabs(Xgap + overlap*1e-3) <= 1e-6);
	  }
	  extern void SintDetail(double X, double Y, int m, int n, int J, int I, int K, int T, FLOAT *Ya, 
				 double &Bias, double &Pen, double &Gauss, double &PenSm, double &OutPen, double &iscore, int verb);/* see RGentigRefScore.h */
	  double Bias,Pen,Gauss,PenSm, OutPen, iscore;
	  SintDetail(Xgap,Ygap,q->sites2[0] - p->sites2[U1-1], q->sites1[0]-q->sitesK1[0]-p->sites1[U1-1], q->sites2[0], q->sites1[0], q->sitesK1[0], p->sitesK1[U1-1], Y, 
		     Bias,Pen,Gauss,PenSm,OutPen,iscore,0);
	  extern RFLOAT OutlierBias;/* see RGentigRefScore.h */

	  r->iscore[U1] = iscore;
	  r->outscore[U1] = OutlierBias + Bias  + Pen + Gauss + PenSm;
	  r->sites1[U1] = q->sites1[0];
	  r->sitesK1[U1] = q->sitesK1[0];
	  r->sites2[U1] = q->sites2[0];
	  
	  if(VERB>=2){
	    printf("Gap: x=%0.4f,y=%0.6f,I=%d,K=%d,J=%d,T=%d:m=%d,n=%d:iscore= %0.6f, outscore= %0.6f (OutlierBias= %0.6f, Bias=%0.6f, Pen=%0.6f, Gauss= %0.6f, PenSm= %0.6f)\n",
		   Xgap,Ygap,q->sites1[0],q->sitesK1[0],q->sites2[0],p->sitesK1[U1-1],q->sites2[0] - p->sites2[U1-1],q->sites1[0]-q->sitesK1[0]-p->sites1[U1-1],
		   iscore,r->outscore[U1],OutlierBias, Bias, Pen, Gauss, PenSm);
	    fflush(stdout);
	  }

	  /* copy q right of first label */
	  for(int t = 1; t < U2; t++){
	    r->iscore[U1 + t] = q->iscore[t];
	    r->outscore[U1 + t] = q->outscore[t];
	    r->sites1[U1 + t] = q->sites1[t];
	    r->sitesK1[U1 + t] = q->sitesK1[t];
	    r->sites2[U1 + t] = q->sites2[t];
	  }

	  /* score of right end */
	  r->iscore[U1 + U2] = q->iscore[U2];
	  r->outscore[U1 + U2] = q->outscore[U2];

	  /* right end */
	  r->Rend = q->Rend;
	  r->Rij1 = q->Rij1;
	  r->Rij2 = q->Rij2;

	  /* compute total score */
	  r->score = 0.0;
	  for(int t = 0; t <= U1+U2; t++)
	    r->score += r->iscore[t];

	  /* compute confidence */
	  double *Xrev = new double[M+2];/* copy of correctly oriented and scaled X[] */
	  FLOAT scale = ScaleFactor[p->scaleID];
	  if(!p->orientation){
	    for(int J=0;J <= M+1;J++)
	      Xrev[J] = X[J]*scale;
	  } else {
	    for(int j=0;j <= M+1;j++)
	      Xrev[j] = (X[M+1]-X[M+1-j])*scale;
	  }

	  extern double alignFP(Calign *r, FLOAT *Y, FLOAT *X, int N, int M, int orientation, int Yid, int Xid, double lscore, int verb);// see refalign.cpp
	  r->logPV = alignFP(r, Y, Xrev, N, M, r->orientation, r->mapid1,r->mapid2, r->score, 0);
	  r->logPV2 = max(r->logPV,r->logPV2);// NEW176
	  if(DEBUG>=2) assert(r->logPV2 >= r->logPV);
	  delete [] Xrev;

	  usexmapentries[i]->align = r;/* this overwrites the pointer : original value is in usexmapentries[i]->orig->align */
	  usexmapentries[i]->confidence = r->logPV;

	  if(VERB>=2){
	    printf("Created temporary merged xid %2i (rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,score=%0.6f)\n",
		   usexmapentries[i]->xmapid, usexmapentries[i]->refcontigid, usexmapentries[i]->qrycontigid, usexmapentries[i]->orientforward ? 1 : 0, usexmapentries[i]->confidence, 
		   usexmapentries[i]->align->score);
	    fflush(stdout);
	  }

	  /* mark MG j for deletion */
	  usexmapentries[j]->sv_overlap = true;
	  filtercnt++;
	}// j = i+1 .. smapsize - 1
      }// i = 0 .. smapsize - 2

      if(DEBUG) assert(filtercnt >= origfiltercnt);
      if(filtercnt == origfiltercnt)/* no more MGs were merged */
	break;// while(1) 

      if(VERB/* HERE >=2 */){
	printf("%d matchgroup merges were performed with gaps <= %0.2f kb : re-running overlap filter: wall time= %0.6f secs\n", filtercnt - origfiltercnt, MergeMG,wtime());
	fflush(stdout);
      }

      /* rerun Overlap Filter */
      int startfiltercnt = filtercnt;
      filtercnt += OverlapFilterMG(usexmapentries, smapsize, numsvIndel1, numsvIndel2, true, verbose);

      /* repeat merging maps and re-running filter loop until no progress is made */
      if(DEBUG) assert(filtercnt >= startfiltercnt);
      if(startfiltercnt == filtercnt)
	break;
    } // while(1)

    /* restore original split matchgroups */
    for(int i = 0; i < smapsize; i++){
      if(usexmapentries[i]->sv_overlap)
	continue;
      
      xmapEntry *p;
      while((p = usexmapentries[i]->orig) != NULL){
	if(DEBUG) assert(p->sv_overlap == false);
	int j = usexmapentries[i]->jMerge;
	if(DEBUG) assert(j > i);
	if(DEBUG) assert(usexmapentries[j]->sv_overlap);
	usexmapentries[j]->sv_overlap = false;
	filtercnt--;
	if(DEBUG) assert(usexmapentries[i]->align != p->align);
	delete usexmapentries[i]->align; 
	delete usexmapentries[i];// only refsites[],qrysites[] arrays are deallocated by xmapEntry destructor, but these have been set to NULL above
	usexmapentries[i] = p;
      }
    }

    /* remove xmapEntries with sv_overlap == true */
    int newsmapsize=0;
    for( int i=0; i < smapsize; i++) {
      if(usexmapentries[i]->sv_overlap)
	continue;
      // HERE      if(usexmapentries[i]->largeMG_only) continue;// HERE HERE : remove after implementing largeMG_only restriction
      usexmapentries[newsmapsize++] = usexmapentries[i];
    }
    if( VERB ){
      printf("%i xmap entries removed due to -MergeMG %0.3f(filtercnt=%d): %i xmap entries remaining (time = %0.6f secs)\n", smapsize-newsmapsize, MergeMG, filtercnt,newsmapsize, wtime());
      fflush(stdout);
    }
    // HERE HERE : uncomment after implementing largeMG_only restriction    if(DEBUG) assert(filtercnt == smapsize-newsmapsize);
    smapsize = newsmapsize;

    if(VERB>=2){
      printf("After MergeMG OverlapFilter: smapsize=%d\n",smapsize);
      for(int i = 0; i < smapsize; i++){
	xmapEntry *p = usexmapentries[i];
	printf("usexmapentries[%d]: xmapid=%d, rid=%lld,qid=%lld,fwd=%d,conf=%0.2f,inversion_only=%d\n",
	       i,p->xmapid,p->refcontigid,p->qrycontigid,p->orientforward ? 1 : 0, p->confidence,
	       p->inversion_only ? 1 : 0);
      }
      fflush(stdout);
    }
  }// if (MergeMG)

  if(DEBUG>=2)
    for(int i = 0; i < smapsize; i++)
      assert(usexmapentries[i]->numrefsites >= usexmapentries[i]->refendidx - usexmapentries[i]->refstartidx + 1);

  if(DEBUG>=2){/* check if any 1-MG Indels no longer overlap the MGs they refer to */
    for(int i = 0; i < numsvIndel; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum != sv_ps);
      if(DEBUG && i < numsvIndel1) assert(p->indices[0] == p->indices[1]);

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;
      
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      if(DEBUG) assert(entry1->output == true);
      if(entry1->output == false)
	continue;

      double qmin = min(entry1->qrystartpos,entry1->qryendpos);
      double qmax = max(entry1->qrystartpos,entry1->qryendpos);
      int iqmin = min(entry1->qrystartidx,entry1->qryendidx);
      int iqmax = max(entry1->qrystartidx,entry1->qryendidx);

      if(!((qmin <= p->querystart && p->querystart <= qmax) &&
	   (qmin <= p->querystop && p->querystop <= qmax) &&
	   (iqmin <= p->querystartidx && p->querystartidx <= iqmax) &&
	   (iqmin <= p->querystopidx && p->querystopidx <= iqmax) &&
	   (entry1->refstartpos <= p->refstart && p->refstart <= entry1->refendpos) &&
	   (entry1->refstartpos <= p->refstop && p->refstop <= entry1->refendpos) &&
	   (entry1->refstartidx <= p->refstartidx && p->refstartidx <= entry1->refendidx) &&
	   (entry1->refstartidx <= p->refstopidx && p->refstopidx <= entry1->refendidx))){
	if(i < numsvIndel1)
	  printf("sv_array[%d]:Single MG indel's MG was trimmed and no longer includes the SV coordinates:\n",i);
	else 
	  printf("sv_array[%d]: 2-MG small inversion or duplication's large MG was trimmed and no longer includes the SV coordinates:\n",i);
	printf("sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	fflush(stdout);
	assert(0);
      }
      if(i < numsvIndel1 ? !(p->type_enum == insertion || p->type_enum == deletion) :
	 i < numsvIndel ? !(p->type_enum == inversion_small_type || p->type_enum == duplication_type || p->type_enum == duplicationinv_type) : 0){
	printf("sv_array[%d]: type_enum=%d, numsvIndel1=%d,numsvIndel2=%d, numsvIndel=%d\n",i,p->type_enum,numsvIndel1,numsvIndel2,numsvIndel);
	printf("\t sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("\t xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	fflush(stdout);
	assert(0);
      }
    }
  }

  if(MATCHGROUP_PARTIAL_TRIM && XmapUnique){

    /* For matchgroups that have been filtered out (sv_overlap == true) but will be output due to Small Inversions/Duplications or 1-MG Indels (output == true), check if they partly overlap another
       larger matchgroup that is NOT filtered out (sv_overlap == false) : If so split off the parts of the overlapped query region corresponding to outliers in the larger MG (if MATCHGROUP_TRIM), then
       trim away the overlapped query region (which should not contain any 1-MG Indels, if XmapUnique, but if there are any filter them out) */

    int orignumxmapentries = numxmapentries;
    //    int origsmapsize = smapsize;
    
    for(int i = 0; i < smapsize; i++){
      xmapEntry *pMGi = usexmapentries[i];
      for(int j = i+1; j < smapsize; j++){
	xmapEntry *pMGj = usexmapentries[j];
	if( pMGi->qrycontigid != pMGj->qrycontigid )
	  break;// since qrycontigid are in ascending order

	double mini = min(pMGi->qrystartpos, pMGi->qryendpos);
	double minj = min(pMGj->qrystartpos, pMGj->qryendpos);
	double maxi = max(pMGi->qrystartpos, pMGi->qryendpos);
	double maxj = max(pMGj->qrystartpos, pMGj->qryendpos);
	double overlap = min(maxi,maxj) - max(mini,minj);
	double minlen = min(maxi-mini, maxj-minj);

	if(overlap >= 0.9999 * minlen)
	  continue;

	if(pMGi->confidence < pMGj->confidence /* WAS73 maxi - mini < maxj - minj*/){
	  if(pMGi->sv_overlap == true && pMGi->output && pMGj->sv_overlap == false){
	    (void)checkOverlap(usexmapentries,smapsize,xmapentries,numxmapentries,maxxmapentries,j,i,false,true);
	  }
	} else {/* pMGj is smaller matchgroup */
	  if(pMGj->sv_overlap == true && pMGj->output && pMGi->sv_overlap == false){	  
	    (void)checkOverlap(usexmapentries,smapsize,xmapentries,numxmapentries,maxxmapentries,i,j,false,true);
	  }
	}
      }
    }

    if(DEBUG && numxmapentries > orignumxmapentries){
      for(int i = orignumxmapentries; i < numxmapentries; i++){
	assert(xmapentries[i]->sv_overlap == true);
	assert(xmapentries[i]->output == true);
      }
    }

    if(numxmapentries > orignumxmapentries){
      if(DEBUG) assert(smapsize == orignumxmapentries);
      smapsize = numxmapentries;
    }

    int newsmapsize = 0;
    for( int i=0; i < smapsize; i++) {
      if(usexmapentries[i]->sv_overlap)
	continue;
      usexmapentries[newsmapsize++] = usexmapentries[i];
    }
    if( VERB ){
      printf("After partial query overlap filter: smapsize= %d -> %d (wall time= %0.6f)\n",smapsize,newsmapsize,wtime());
      fflush(stdout);
    }
    smapsize = newsmapsize;

    if(numxmapentries > orignumxmapentries)
      qsort(usexmapentries, smapsize, sizeof(usexmapentries[0]), (intcmp *)xmapEntryQryidInc);

    if(VERB>=3){
      for(int i = 0; i < numxmapentries; i++){
	xmapEntry *p = usexmapentries[i];
	if(1 /* p->xmapid == 868 || p->xmapid == 866*/){
	  printf("usexmapentries[%d]: refid=%lld,qryid=%lld, xmapid=%d, xmapid_overlapped=%d: logPV= %0.2f\n",
		 i, p->refcontigid,p->qrycontigid,p->xmapid,p->xmapid_overlapped, p->confidence);
	  fflush(stdout);
	}
      }
    }
  }

  if(DEBUG>=2){/* check if any 1-MG Indels no longer overlap the MGs they refer to */
    for(int i = 0; i < numsvIndel; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum != sv_ps);
      if(DEBUG && i < numsvIndel1) assert(p->indices[0] == p->indices[1]);

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;
      
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      if(DEBUG) assert(entry1->output == true);
      if(entry1->output == false)
	continue;

      double qmin = min(entry1->qrystartpos,entry1->qryendpos);
      double qmax = max(entry1->qrystartpos,entry1->qryendpos);
      int iqmin = min(entry1->qrystartidx,entry1->qryendidx);
      int iqmax = max(entry1->qrystartidx,entry1->qryendidx);

      if(!((qmin <= p->querystart && p->querystart <= qmax) &&
	   (qmin <= p->querystop && p->querystop <= qmax) &&
	   (iqmin <= p->querystartidx && p->querystartidx <= iqmax) &&
	   (iqmin <= p->querystopidx && p->querystopidx <= iqmax) &&
	   (entry1->refstartpos <= p->refstart && p->refstart <= entry1->refendpos) &&
	   (entry1->refstartpos <= p->refstop && p->refstop <= entry1->refendpos) &&
	   (entry1->refstartidx <= p->refstartidx && p->refstartidx <= entry1->refendidx) &&
	   (entry1->refstartidx <= p->refstopidx && p->refstopidx <= entry1->refendidx))){
	if(i < numsvIndel1)
	  printf("sv_array[%d]:Single MG indel's MG was trimmed and no longer includes the SV coordinates:\n",i);
	else 
	  printf("sv_array[%d]: 2-MG small inversion or duplication's large MG was trimmed and no longer includes the SV coordinates:\n",i);
	printf("sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	printf("numsvIndel1=%d,numsvIndel2=%d,numsvIndel=%d\n",numsvIndel1,numsvIndel2,numsvIndel);
	fflush(stdout);
	assert(0);
      }
      if(i < numsvIndel1 ? !(p->type_enum == insertion || p->type_enum == deletion) :
	 i < numsvIndel ? !(p->type_enum == inversion_small_type || p->type_enum == duplication_type || p->type_enum == duplicationinv_type) : 0){
	printf("sv_array[%d]: type_enum=%d, numsvIndel1=%d,numsvIndel2=%d, numsvIndel=%d\n",i,p->type_enum,numsvIndel1,numsvIndel2,numsvIndel);
	printf("\t sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("\t xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	fflush(stdout);
	assert(0);
      }
    }
  }

  if(DEBUG>=2)
    for(int i = 0; i < smapsize; i++)
      assert(usexmapentries[i]->numrefsites >= usexmapentries[i]->refendidx - usexmapentries[i]->refstartidx + 1);

  if(DEBUG>=2){
    for(int i = 0; i < smapsize;i++){
      xmapEntry *p = usexmapentries[i];

      /* check p->refsites[] (bp) vs refendpos (bp) */
      assert(p->numrefsites >= p->refendidx - p->refstartidx + 1);
      int Istart = p->refstartidx;
      for(int I = Istart; I <= p->refendidx; I++){
	double sitesI = p->refsites[I - Istart];
	if(DEBUG && !(p->refstartpos <= sitesI && sitesI <= p->refendpos)){
	  printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refstartpos= %0.1f,refendpos= %0.1f\n",
		 i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart, p->refsites[I-Istart],p->refstartpos,p->refendpos);
	  fflush(stdout);
	  assert(p->refstartpos <= sitesI && sitesI <= p->refendpos);
	}
	if(I > p->refstartidx){
	  if(DEBUG && !(p->refsites[I-1 - Istart] <= sitesI)){
	    printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refsites[%d]= %0.1f, refendpos= %0.1f\n",
		   i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart,p->refsites[I - Istart],I-1-Istart,p->refsites[I-1-Istart],p->refendpos);
	    fflush(stdout);
	    assert(p->refsites[I-1-Istart] <= sitesI);
	  }
	}
      }
      /* check p->qrysites[0..numqrysites-1] (bp) vs refendpos (bp) */
      double qrystartpos = min(p->qrystartpos,p->qryendpos);
      double qryendpos = max(p->qrystartpos,p->qryendpos);
      int Imin = min(p->qrystartidx,p->qryendidx);
      int Imax = max(p->qrystartidx,p->qryendidx);
      if(DEBUG) assert(p->numqrysites >= Imax-Imin+1);
      for(int I = Imin; I <= Imax; I++){
	double sitesI = p->qrysites[I - Imin];
	if(DEBUG && !(qrystartpos <= sitesI && sitesI <= qryendpos)){
	  printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		 i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites, I, I-Imin,p->qrysites[I-Imin],qrystartpos,qryendpos);
	  fflush(stdout);
	  assert(qrystartpos <= sitesI && sitesI <= qryendpos);
	}
	if(I > Imin){
	  if(DEBUG && !(p->qrysites[I-1-Imin] <= sitesI)){
	    printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		   i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites,I,I-Imin,p->qrysites[I-Imin],I-1-Imin,p->qrysites[I-1-Imin],qrystartpos,qryendpos);
	    fflush(stdout);
	    assert(p->qrysites[I-1-Imin] <= sitesI);
	  }
	}
      }
    }
  }

  // Store pairs of indices from usexmapentries[0..smapsize-1] which have the same query contig;
  // NOTE : usexmapentries[] was previously sorted by qrycontigid, refcontigid, refstartpos, refendpos in ascending order
  int indices_size = smapsize*2; //dynamically grow in getSameQueryIndices
  int *indices1 = new int[indices_size], *indices2 = new int[indices_size]; 
  int indiceslen = getSameQueryIndices(usexmapentries, smapsize, indices1, indices2, indices_size);

  //must increase size of global array by fixed value indiceslen
  //note in loop below, index in sv_array is no longer loop index due to -indels added
  printf("sv array: numsv=%i, indiceslen=%i, maxsv=%i: wall time= %0.6f secs\n", numsv, indiceslen, maxsv,wtime()); fflush(stdout); //debug
  if(DEBUG) assert(indiceslen >= 0);
  if( numsv + indiceslen > maxsv ) //if not, no need to increase size
    maxsv = growArray(sv_array, maxsv, max(numsv * 3 / 2, numsv + indiceslen + 16)); //newsize is third arg

  if( VERB ){
    printf("Starting main sv loop (wall time=%.6f)\n", wtime());
    printf("Query pairs: %i (smapsize=%d, indices_size=%d, NumSmallInversions=%d, maxsv=%d)\n", indiceslen, smapsize, indices_size, NumSmallInversions, maxsv);
    fflush(stdout);
  }
  if(DEBUG>=2){
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
    }
  }

  //printf("Translocation params: maxsep:%.0f ovlfr:%.3f ovlmx:%.0f\n", smap_translocation_maxsep, smap_translocation_ovlfr, smap_translocation_ovlmx); //DEBUG

  int repeatcnt = 0, repeatcnt2 = 0, swapcnt = 0, trimcnt = 0;
  size_t cumrepeatcnt2 = 0;

 LmainRepeat:

  //for each element surviving above requirement (usexmapentries), find SVs and smap_sv_sizekb requirements
  for(int indp = 0; indp < indiceslen; indp++) {
    if(VERB && (indp > 0 && !(indp % 100000))){
      printf("sv main loop indp=%d/%d: numsv= %d/%d: wall time= %0.6f\n",indp,indiceslen,numsv,maxsv,wtime());
      fflush(stdout);
    }

    int i = indices1[indp]; //'previous' index
    int j = indices2[indp]; //'current' index

    if(DEBUG) assert(0 <= i && i < j && j < smapsize);

    xmapEntry *pMGi = usexmapentries[i];
    xmapEntry *pMGj = usexmapentries[j];

    if(repeatcnt2 > 0)
      cumrepeatcnt2 += repeatcnt2;
    repeatcnt2 = 0;

  Lrepeat:
    if(MATCHGROUP_PARTIAL_TRIM && (pMGi->sv_overlap || pMGj->sv_overlap)) // check if checkOverlap() filtered out either MG
      continue;

    if(MATCHGROUP_PARTIAL_TRIM && pMGi->refcontigid == pMGj->refcontigid && !(pMGi->refstartpos <= pMGj->refstartpos && (pMGi->refstartpos < pMGj->refstartpos || pMGi->refendpos <= pMGj->refendpos))){
    Lrepeat2:
      if(VERB>=2){
	printf("SV main loop indp=%d/%d: swapping i=%d,j=%d (xmapids=%d,%d): pMGi= %0.3f .. %0.3f kb, pMGj= %0.3f .. %0.3f kb (on reference): repeatcnt=%d, numsv=%d\n",
	       indp,indiceslen,i,j,pMGi->xmapid,pMGj->xmapid,pMGi->refstartpos*1e-3,pMGi->refendpos*1e-3,pMGj->refstartpos*1e-3,pMGj->refendpos*1e-3,repeatcnt, numsv);
	fflush(stdout);
      }
      if(DEBUG) assert(repeatcnt==0);
      swapcnt++;

      /* swap i,j */
      indices1[indp] = j;
      indices2[indp] = i;

      i = indices1[indp]; //'previous' index
      j = indices2[indp]; //'current' index

      pMGi = usexmapentries[i];
      pMGj = usexmapentries[j];
    }

    if(VERB>=2){
      printf("SV main loop indp=%d/%d:i=%d,j=%d (xmapids=%d,%d), numsv=%d\n",indp,indiceslen,i,j,pMGi->xmapid,pMGj->xmapid, numsv);
      fflush(stdout);
    }

    if(MATCHGROUP_TRIM && XmapUnique && pMGi->xmapid_overlapped > 0 && pMGj->xmapid != pMGi->xmapid_overlapped){
      if(VERB>=2){
	printf("SV main loop indp=%d/%d:i=%d,j=%d:xmapids=%d,%d: skipping since pMGi->xmapid_overlapped=%d\n",indp,indiceslen,i,j,pMGi->xmapid,pMGj->xmapid,pMGi->xmapid_overlapped);
	fflush(stdout);
      }
      continue;
    }
    if(MATCHGROUP_TRIM && XmapUnique && pMGj->xmapid_overlapped > 0 && pMGi->xmapid != pMGj->xmapid_overlapped){
      if(VERB>=2){
	printf("SV main loop indp=%d/%d:i=%d,j=%d:xmapids=%d,%d: skipping since pMGj->xmapid_overlapped=%d\n",indp,indiceslen,i,j,pMGi->xmapid,pMGj->xmapid,pMGj->xmapid_overlapped);
	fflush(stdout);
      }
      continue;
    }

    //svs are only allowed to be within a single query contig--check current entry has same ids as prev
    //this should be satisfied by indices arrays--change if to assert
    if(DEBUG) assert( pMGi->qrycontigid == pMGj->qrycontigid );
    long long contigid = pMGi->qrycontigid;
    Cmap *Xmap = XXmap[pMGi->align->mapid2];
    FLOAT *X = Xmap->site[0];
    int M = Xmap->numsite[0];

    long long refid1 = pMGi->refcontigid;
    Cmap *Ymap1 = YYmap[pMGi->align->mapid1];
    FLOAT *Y1 = Ymap1->site[0];
    int N1 = Ymap1->numsite[0];

    //be sure neither xmap entry was flagged as overlap above--these should have been removed from usexmapentries already
    if(DEBUG) assert( !pMGi->sv_overlap && !pMGj->sv_overlap );

    /* first compute some useful values */
    bool refoverlap = overlap(pMGi->refstartpos, pMGi->refendpos,
			      pMGj->refstartpos, pMGj->refendpos);
    double refoverlapsize = overlapLength(pMGi->refstartpos, pMGi->refendpos,
					  pMGj->refstartpos, pMGj->refendpos) * 1e-3; //kb!
    int refoverlapSites = min(pMGi->refendidx, pMGj->refendidx) - max(pMGi->refstartidx, pMGj->refstartidx);

    double minq1 = min(pMGi->qrystartpos, pMGi->qryendpos);
    double maxq1 = max(pMGi->qrystartpos, pMGi->qryendpos);
    double minq2 = min(pMGj->qrystartpos, pMGj->qryendpos);
    double maxq2 = max(pMGj->qrystartpos, pMGj->qryendpos);

    bool qryoverlap = overlap( minq1, maxq1, minq2, maxq2 );
    double qryoverlapsize = overlapLength( minq1, maxq1, minq2, maxq2 ) * 1e-3; //kb!

    int iminq1 = min(pMGi->qrystartidx, pMGi->qryendidx);
    int imaxq1 = max(pMGi->qrystartidx, pMGi->qryendidx);
    int iminq2 = min(pMGj->qrystartidx, pMGj->qryendidx);
    int imaxq2 = max(pMGj->qrystartidx, pMGj->qryendidx);
    int qryoverlapSites = min(imaxq1, imaxq2) - max(iminq1, iminq2);

    // locate range of usexmapentries[] with same qrycontigid as i, j 
    int kmin = min(i,j), kmax = max(i,j);
    for(;--kmin >= 0;){
      xmapEntry *pMGk = usexmapentries[kmin];	
      if(pMGk->qrycontigid != contigid){
	kmin++;
	break;
      }
    }
    kmin = max(0,kmin);
    if(kmin == min(i,j))
      kmin++;

    for(;++kmax < smapsize;){
      xmapEntry *pMGk = usexmapentries[kmax];
      if(pMGk->qrycontigid != contigid){
	kmax--;
	break;
      }
    }
    kmax = min(smapsize-1, kmax);
    if(kmax == max(i,j))
      kmax--;

    if(VERB>=2){
      printf("SV main loop indp=%d/%d: i=%d,j=%d: xmapids= %d,%d qryid= %lld,%lld, refid= %lld,%lld, fwd=%d,%d, inv_only=%d,%d:\n",
	     indp,indiceslen,i,j,pMGi->xmapid,pMGj->xmapid,pMGi->qrycontigid,pMGj->qrycontigid,pMGi->refcontigid,pMGj->refcontigid,pMGi->orientforward,pMGj->orientforward,
	     pMGi->inversion_only,pMGj->inversion_only);
      printf("\t ri=%d..%d(%0.3f..%0.3f),rj=%d..%d(%0.3f..%0.3f),qi=%d..%d(%0.3f..%0.3f),qj=%d..%d(%0.3f..%0.3f),kmin=%d,kmax=%d\n",
	     pMGi->refstartidx,pMGi->refendidx, pMGi->refstartpos * 1e-3, pMGi->refendpos * 1e-3,
	     pMGj->refstartidx,pMGj->refendidx, pMGj->refstartpos * 1e-3, pMGj->refendpos * 1e-3,
	     iminq1,imaxq1, minq1 * 1e-3, maxq1 * 1e-3,  iminq2,imaxq2, minq2 * 1e-3, maxq2 * 1e-3,kmin,kmax);
      fflush(stdout);
    }

    //a rearrangement is when one query contig aligns to two reference contigs
    //in this case, do not make any restrictions on refsize 
    const bool rearrangement = ( pMGi->refcontigid != pMGj->refcontigid );

    if(rearrangement){
      // NEW79 : check if both matchgroups have confidence >= LogPvThreshold (and not just LogPvThreshold2, which could be the small inversion/duplication threshold)
      if(pMGi->confidence < LogPvThreshold || pMGj->confidence < LogPvThreshold){
	if(verbose){
	  printf("xmapid=%d,%d: qryid=%lld,refid=%lld,%lld, conf= %0.2f, %0.2f : skipping translocation (confidence of both MGs must be >= %0.2f)\n",pMGi->xmapid,pMGj->xmapid,
		 contigid,pMGi->refcontigid,pMGj->refcontigid,pMGi->confidence,pMGj->confidence,LogPvThreshold);
	  fflush(stdout);
	}
	continue;
      }

      // NEW87 : check if either matchgroup is flagged as a repeat : If so skip translocation calls 
      if(pMGi->align->repeat > 0 || pMGj->align->repeat > 0){
	if(verbose){
	  printf("xmapid=%d,%d: qryid=%lld,refid=%lld,%lld, conf= %0.2f, %0.2f, Repeat Mask= %d, %d : skipping translocation (both MGs must not have repeat pattern > 0)\n",
		 pMGi->xmapid,pMGj->xmapid,
		 contigid,pMGi->refcontigid,pMGj->refcontigid,pMGi->confidence,pMGj->confidence, pMGi->align->repeat, pMGj->align->repeat);
	  fflush(stdout);
	}
	continue;
      }

      // NEW79 : check if another same qry MG overlaps most of one of original 2 MGs and maps back to the other MGs's ref contig : if so disallow translocation call
      // NEW150 : Also check if another same qry MG overlaps most of the gap between the original 2 MGs : if so disallow translocation call

      int k = kmin;
      for(; k <= kmax; k++){
	if(k==i || k==j)
	   continue;
	xmapEntry *pMGk = usexmapentries[k];
	if(DEBUG>=2) assert(pMGk->qrycontigid == contigid);
	
	double minq3 = min(pMGk->qrystartpos, pMGk->qryendpos);
	double maxq3 = max(pMGk->qrystartpos, pMGk->qryendpos);
	
	if(VERB>=3){
	  printf("i=%d,j=%d,k=%d:xmapid=%d,%d: qryid=%lld,refid=%lld,%lld, qry1= %0.3f..%0.3f kb, qry2= %0.3f .. %0.3f kb : xmapid=%d(qryid=%lld,refid=%lld), qry3= %0.3f..%0.3f kb\n",
		 i,j,k,pMGi->xmapid,pMGj->xmapid,contigid,pMGi->refcontigid,pMGj->refcontigid,minq1*1e-3,maxq1*1e-3,minq2*1e-3,maxq2*1e-3,
		 pMGk->xmapid,pMGk->qrycontigid,pMGk->refcontigid,minq3*1e-3,maxq3*1e-3);
	  fflush(stdout);
	}

	if(smap_TransMaxOverlap > 0.0 && pMGk->refcontigid == pMGj->refcontigid && overlapLength(minq1,maxq1, minq3, maxq3) >= (maxq1-minq1) * smap_TransMaxOverlap){
	  if(verbose){
	    printf("xmapid=%d,%d: qryid=%lld,refid=%lld,%lld : xmapid=%d(qryid=%lld,refid=%lld) overlaps first MG on qry by %0.1f%% (-svTransMaxOverlap): skipping translocation\n",pMGi->xmapid,pMGj->xmapid,
		   contigid,pMGi->refcontigid,pMGj->refcontigid,pMGk->xmapid,pMGk->qrycontigid,pMGk->refcontigid,smap_TransMaxOverlap*100.0);
	    fflush(stdout);
	  }
	  break;
	}
	if(smap_TransMaxOverlap > 0.0 && pMGk->refcontigid == pMGj->refcontigid && pMGk->confidence >= LogPvThreshold &&
	   overlapLength(min(maxq2,minq1),max(maxq1,minq2), minq3, maxq3) >= (maxq1-minq1) * smap_TransMaxOverlap){
	  if(verbose){
	    printf("xmapid=%d,%d: qryid=%lld,refid=%lld,%lld : xmapid=%d(qryid=%lld,refid=%lld,conf=%0.2f) overlaps first MG + gap on qry by %0.1f%% (-svTransMaxOverlap): skipping translocation\n",pMGi->xmapid,pMGj->xmapid,
		   contigid,pMGi->refcontigid,pMGj->refcontigid,pMGk->xmapid,pMGk->qrycontigid,pMGk->refcontigid,pMGk->confidence,smap_TransMaxOverlap*100.0);
	    fflush(stdout);
	  }
	  break;
	}

	if(smap_TransMaxOverlap > 0.0 && pMGk->refcontigid == pMGi->refcontigid && overlapLength(minq2,maxq2, minq3, maxq3) >= (maxq2-minq2) * smap_TransMaxOverlap){
	  if(verbose){
	    printf("xmapid=%d,%d: qryid=%lld,refid=%lld,%lld : xmapid=%d(qryid=%lld,refid=%lld) overlaps second MG on qry by %0.1f%% (-svTransMaxOverlap): skipping translocation\n",pMGi->xmapid,pMGj->xmapid,
		   contigid,pMGi->refcontigid,pMGj->refcontigid,pMGk->xmapid,pMGk->qrycontigid,pMGk->refcontigid,smap_TransMaxOverlap*100.0);
	    fflush(stdout);
	  }
	  break;
	}
	if(smap_TransMaxOverlap > 0.0 && pMGk->refcontigid == pMGi->refcontigid && pMGk->confidence >= LogPvThreshold && 
	   overlapLength(min(maxq1,minq2),max(maxq2,minq2), minq3, maxq3) >= (maxq2-minq2) * smap_TransMaxOverlap){
	  if(verbose){
	    printf("xmapid=%d,%d: qryid=%lld,refid=%lld,%lld : xmapid=%d(qryid=%lld,refid=%lld) overlaps second MG + gap on qry by %0.1f%% (-svTransMaxOverlap): skipping translocation\n",pMGi->xmapid,pMGj->xmapid,
		   contigid,pMGi->refcontigid,pMGj->refcontigid,pMGk->xmapid,pMGk->qrycontigid,pMGk->refcontigid,smap_TransMaxOverlap*100.0);
	    fflush(stdout);
	  }
	  break;
	}

	if(smap_TransMaxGapOverlap > 0.0 && pMGk->confidence >= LogPvThreshold){
	  if(minq2 > maxq1 ? overlapLength(maxq1, minq2, minq3, maxq3) >= (maxq3-minq3) * smap_TransMaxGapOverlap :
	     minq1 > maxq2/*WAS405 maxq3 */ ? overlapLength(maxq2, minq1, minq3, maxq3) >= (maxq3-minq3) * smap_TransMaxGapOverlap : false){
	    if(verbose){
	      printf("xmapid=%d,%d:qryid=%lld,refid=%lld,%lld : xmapid=%d(qryid=%lld,refid=%lld) overlaps qry gap by %0.1f%% (-svTransMaxGapOverlap): skipping translocation\n",pMGi->xmapid,pMGj->xmapid,
		     contigid,pMGi->refcontigid,pMGj->refcontigid,pMGk->xmapid,pMGk->qrycontigid,pMGk->refcontigid,smap_TransMaxGapOverlap*100.0);
	      fflush(stdout);
	    }
	    break;
	  }
	}

      }
      if(k <= kmax)
	continue;
    }

    double refsize = (pMGj->refstartpos - pMGi->refendpos) * 1e-3; //convert to kb

    bool inversion_only = (pMGi->inversion_only || pMGj->inversion_only);

    // New code to handle overlapped matchgroups that were NOT deleted due to -InversionOverlapFilter AND MG pair cannot be inversion
    if(inversion_only && (rearrangement || pMGi->orientforward == pMGj->orientforward)){
      if(VERB>=2){
	printf("\txmapid=%d,%d:inversion_only=%d,%d,rearrangement=%d,fwd=%d,%d, qryoverlap=%d,repeatcnt=%d: skipping pair due to inversion_only\n",
	       pMGi->xmapid,pMGj->xmapid, pMGi->inversion_only, pMGj->inversion_only,
	       rearrangement?1:0,pMGi->orientforward?1:0,pMGj->orientforward?1:0, qryoverlap, repeatcnt);
	fflush(stdout);
      }

      if(MATCHGROUP_PARTIAL_TRIM >= MATCHGROUP_PARTIAL_TRIM_COMPLEX/* NEW151 */ && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	int origsmapsize = smapsize;
	if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
	  trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries, j, i, false, false);
	} else {/* trim usexmapentries[j] where it overlaps on query */
	  trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries, i, j, false, false);
	}
	if(DEBUG) assert(origsmapsize == smapsize);
      }

      continue;/* inversion_only matchgroups can only be used to call inversions */
    }

    if(inversion_only && max(pMGi->confidence,pMGj->confidence) <= LogPvThreshold){/* At least one matchgroup must satisfy -T for inversion_only calls */
      if(VERB>=2){
	printf("\txmapid=%d,%d:inversion_only=%d,%d,conf=%0.2f,%0.2f, qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld : skipping pair since they are both below -T threshold\n",
	       pMGi->xmapid,pMGj->xmapid, pMGi->inversion_only, pMGj->inversion_only,
	       pMGi->confidence,pMGj->confidence,qryoverlap,repeatcnt,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
	fflush(stdout);
      }

      if(MATCHGROUP_PARTIAL_TRIM >= MATCHGROUP_PARTIAL_TRIM_COMPLEX /* NEW151 */ && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	int origsmapsize = smapsize;
	if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
	  trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	} else {/* trim usexmapentries[j] where it overlaps on query */
	  trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	}
	if(DEBUG) assert(origsmapsize == smapsize);
      }

      continue;
    }

    if(inversion_only){/* Make sure current 2 matchgroups are not themselves overlapped in Query by smap_InvMaxOverlapSize (or based on OverlapFilterQryRatio) */
      double minlen = min(maxq1 - minq1, maxq2 - minq2) * 1e-3;      
      if(VERB>=2){
        if(qryoverlapsize > (smap_InvMaxOverlapSize ? smap_InvMaxOverlapSize : minlen * OverlapFilterQryRatio))
	  printf("\t inversion_only=%d,%d: qryoverlapsize= %0.3f,minlen= %0.3f, svInvMaxOverlapSize= %0.3f,OverlapFilterQryRatio= %0.6f, qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld: skipping due to overlap\n",
		 pMGi->inversion_only,pMGj->inversion_only,qryoverlapsize,minlen,smap_InvMaxOverlapSize,OverlapFilterQryRatio,qryoverlap,repeatcnt,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
	else
	  printf("\t inversion_only=%d,%d: qryoverlapsize= %0.3f,minlen= %0.3f, svInvMaxOverlapSize= %0.3f,OverlapFilterQryRatio= %0.6f\n",
		 pMGi->inversion_only,pMGj->inversion_only,qryoverlapsize,minlen,smap_InvMaxOverlapSize,OverlapFilterQryRatio);
	fflush(stdout);
      }

      if(qryoverlapsize > (smap_InvMaxOverlapSize ? smap_InvMaxOverlapSize : minlen * OverlapFilterQryRatio)){
	if(smap_InvMaxOverlapSize <= 0.0/*NEW151*/ && MATCHGROUP_PARTIAL_TRIM >= MATCHGROUP_PARTIAL_TRIM_COMPLEX && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	  int origsmapsize = smapsize;
	  if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
	    trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	  } else {/* trim usexmapentries[j] where it overlaps on query */
	    trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	  }
	  if(DEBUG) assert(origsmapsize == smapsize);
	}

	continue;
      }
    }

    // ref entries are always sorted, so refstart < refstop, for each entry, but there is NO requirement on the end
    // second (j) can be shorter than first (i), and therefore end before first.
    // This guarantees that refsize < 0 implies overlap IF chromosomes are same.
    if( !( pMGi->refcontigid != pMGj->refcontigid ||
	   pMGj->refstartpos >= pMGi->refstartpos ) ) {
      printf("ERROR: xmap entries out of order: %i %i\n", pMGi->xmapid, pMGj->xmapid);
      fflush(stdout);
      assert( pMGi->refcontigid != pMGj->refcontigid ||
	      pMGj->refstartpos >= pMGi->refstartpos );
    }

    if(VERB>=2){
      printf("SV main loop indp=%d/%d:rearrangement=%d, refsize=%0.3f\n",indp,indiceslen, rearrangement ? 1 : 0, refsize);

      printf("\tpMGi:fwd=%d,refstartpos=%0.3f,refendpos=%0.3f,qryendpos=%0.3f,qrystartpos=%0.3f kb, conf= %0.2f\n",
	     pMGi->orientforward, pMGi->refstartpos*1e-3,pMGi->refendpos*1e-3,pMGi->qryendpos*1e-3, pMGi->qrystartpos*1e-3,pMGi->confidence);
      printf("\tpMGj:fwd=%d,refstartpos=%0.3f,refendpos=%0.3f,qryendpos=%0.3f,qrystartpos=%0.3f kb, conf= %0.2f\n",
	     pMGj->orientforward ? 1 : 0, pMGj->refstartpos*1e-3, pMGj->refendpos*1e-3, pMGj->qryendpos*1e-3, pMGj->qrystartpos*1e-3,pMGj->confidence);
      fflush(stdout);
    }

    /* HERE HERE : determine locally optimal orientation of query based on MGi, MGj AND usexmapentries[i-1,i+1,j-1,j+1]. Current method based on only MGi, MGj handles some cases incorrectly
                        and misses (eg) case with both MGi and MGj being inverted matchgroups with an Indel between them : this indel is currently never called
     */

    // since entries are sorted by ref pos, check if MGj is before or after MGi on query
    double qrysize = 0, qrystart = 0, qrystop = 0;
    int qrystartidx = 0, qrystopidx = 0;

    bool inversion = false; //opposite orientation, no overlap
    bool opporient = false; //opposite orientation, regardless of overlap--must distinguish from indel
    //new convention: j = current, i = previous (confusing? Yes.)
    //now, classify based on orientation
    //orientation here is flipped from the align struct--here, true is forward, false is backwards

    if( pMGi->orientforward && pMGj->orientforward ) { //both forward
      qrystart    = pMGi->qryendpos;
      qrystartidx = pMGi->qryendidx;
      qrystop     = pMGj->qrystartpos;
      qrystopidx  = pMGj->qrystartidx;
    }
    else if( !pMGi->orientforward && !pMGj->orientforward ) { //both backward
      qrystart    = pMGj->qrystartpos;
      qrystartidx = pMGj->qrystartidx;
      qrystop     = pMGi->qryendpos;
      qrystopidx  = pMGi->qryendidx;
    }
    else { //opposite orientation
      opporient = true;
      if(minq1 < minq2) { //previous query (MG i) starts before current query (MG j) starts
	qrystart = maxq1;
	qrystartidx = imaxq1;
	qrystop  = minq2;
	qrystopidx  = iminq2;
      } else if( minq2 < minq1) { //current query (MG j) starts before previous query (MG i) starts
	qrystart = maxq2;
	qrystartidx = imaxq2;
	qrystop  = minq1;
	qrystopidx  = iminq1;
      } else { // minq1 == minq2 (NOTE : this case is not possible if OverlapFilterQryRatio < 1.0, since then the smaller matchgroups in the query would have been filtered out)
	if(DEBUG && /* NEW109 */ !MATCHGROUP_PARTIAL_TRIM && !(OverlapFilterQryRatio >= 0.9999 || OverlapFilterConfDelta > 0.0)){
          printf("SV indp=%d/%d:i=%d,j=%d:refid=%lld,%lld,qryid=%lld,%lld,fwd=%d,%d\n",indp,indiceslen,i,j,pMGi->refcontigid,pMGj->refcontigid,
		 pMGi->qrycontigid,pMGj->qrycontigid,pMGi->orientforward ? 1 : 0, pMGj->orientforward ? 1 : 0);
	  printf("\t qrypos[i] = %0.4f .. %0.4f, qrypos[j] = %0.4f .. %0.4f, refpos[i] = %0.4f .. %0.4f, refpos[j] = %0.4f .. %0.4f\n",
		 pMGi->qrystartpos, pMGi->qryendpos, pMGj->qrystartpos, pMGj->qryendpos,
		 pMGi->refstartpos, pMGi->refendpos, pMGj->refstartpos, pMGj->refendpos);
	  printf("\t OverlapFilterQryRatio= %0.6f, OverlapFilterConfDelta= %0.3f\n",OverlapFilterQryRatio, OverlapFilterConfDelta);
	  fflush(stdout);
          assert(OverlapFilterQryRatio >= 0.9999 || OverlapFilterConfDelta > 0.0);
        }
#if 1 // NEW CODE
	if(maxq1 > maxq2) {/* treat qry MG j as left of qry MG i */
	  qrystart = maxq2;
	  qrystartidx = imaxq2;
	  qrystop  = minq1;
	  qrystopidx  = iminq1;
	} else { /* maxq1 <= maxq2 : treat qry MG i as left of qry MG j */
	  qrystart = maxq1;
	  qrystartidx = imaxq1;
	  qrystop  = minq2;
	  qrystopidx  = iminq2;
	}
#else // OLD CODE
	qrystart = minq1;
	qrystartidx = iminq1; // WAS30 (pMGi->qrystartpos < pMGi->qryendpos) ? pMGi->qrystartidx : pMGi->qryendidx;
	if(maxq1 > maxq2) {/* treat qry MG j as left of qry MG i */
	  qrystop = maxq2;
	  qrystopidx = imaxq2; // WAS30 (pMGj->qrystartpos > pMGj->qryendpos) ? pMGj->qrystartidx : pMGj->qryendidx;
	} else { /* maxq1 <= maxq2 : treat qry MG i as left of qry MG j */
	  qrystop = maxq1;
	  qrystopidx = imaxq1; // WAS30 (pMGi->qrystartpos > pMGi->qryendpos) ? pMGi->qrystartidx : pMGi->qryendidx;
	}
#endif // OLD CODE
      }
      if(VERB>=2){
	printf("i=%d,j=%d:xmapid1=%d,xmapid2=%d:qryid1=%lld,qryid2=%lld,refid1=%lld,refid2=%lld, qryoverlap= %d, rearrangement= %d, inversion_only= %d, qrystart=%0.3f,qrystop=%0.3f\n",
	       i,j,pMGi->xmapid,pMGj->xmapid,pMGi->qrycontigid,pMGj->qrycontigid,refid1,pMGj->refcontigid,
	       qryoverlap ? 1 : 0, rearrangement ? 1 : 0, inversion_only ? 1 : 0, qrystart*1e-3,qrystop*1e-3);
	fflush(stdout);
      }
    } /* opposite orientation */


    /* determine which matchgroup should be treated as inverted matchgroup (for inversion and inverted duplication) :
       Assume OverlapFilterQryRatio < 1.0 : This means matchgroups cannot overlap completely on query and hence ordering on query is always unambigous.
       1. If MG j is not completely overlapped by MG i, orient query so MG i is left of MG j.
          If MG i is completely overlapped on reference by larger MG j (left end on reference are the same label), special handling may be required :
	      1A) If MG i is inverted, then inverted duplication is the only possibility : no special handling required, but need to rule out inversion.
	      1B) If MG j is inverted, can call both inversion (of non-overlapped part of MG j) and inverted duplication (of MG i) :
	          Treat MG i as inverted MG, unless Duplication Spacing (measured on ref(???)) is larger than smap_InvDup_maxSpacing, since we prefer to call inverted duplications over regular inversions
       2. If MG j is completely overlapped on reference by larger (or equal size) MG i, then MG j is the inverted MG for calling inverted duplications, but MG i is the inverted MG for calling regular inversions!
          Treat MG j as inverted MG, unless the Duplication Spacing (measured via reference) is larger than smap_InvDup_maxSpacing, since we prefer to call inverted duplications over regular inversions
          (Could also call both by Treating MG i as inverted MG, but inversion call may be very small or zero size)
    */
    double inv_refalignlen = 0;/* ref size of MG that is treated as inverted (only computed if opporient) */
    int inv_refalignSites = 0;/* ref sites in MG that is treated as inverted ( only computed if opporient) */
    double inv_qryalignlen = 0;/* query size of MG to be treated as inverted (only computed if opporient) */ 
    int inv_qryalignSites = 0;/* query sites of MG to be treated as inverted (only computed if opporient) */
    int noninv_refalignSites = 0;/* ref sites in MG that is treated as non-inverted (only computed if opporient) */
    bool inv_MGi = true;/* If MG i is treated as the inverted matchgroup for calling inversions (only valid if opporient==1) */
    bool NoInvDup = false;/* set to true if (MG j is completely overlapped on reference by MG i and inv_MGi == true) OR (MG i is completely overlapped on reference by MG j and inv_MGi == false) */
    bool NoInv = false;/* set to true in case 1A above, to rule out Inversion if Inverted Duplication fails because of Duplication Spacing */
    bool MGi_qryleft = true;/* If Matchgroup i is left of Matchgroup j on query (only computed if opporient==1 AND an inversion call is made) */

    if(opporient && !rearrangement) {
      if(OverlapFilterQryRatio >= 0.99){
	printf("Inversion handling requires OverlapFilterQryRatio <= 0.99 (current value is %0.6f)\n",OverlapFilterQryRatio);
	fflush(stdout);exit(1);
      }

      if(DEBUG && !(pMGj->refstartpos >= pMGi->refstartpos)){
	printf("i=%d,j=%d:refstartpos= %0.1f,%0.1f (refid=%lld,%lld,qryid=%lld,%lld)\n",
	       i,j,pMGi->refstartpos,pMGj->refstartpos,pMGi->refcontigid,pMGi->refcontigid,pMGi->qrycontigid,pMGj->qrycontigid);
	fflush(stdout);
	assert(pMGj->refstartpos >= pMGi->refstartpos);
      }

      if (pMGi->refendpos < pMGj->refendpos){/* MG j is NOT completely overlapped by MG i on reference : orient query so MG i is left of MG j */
	if (min(pMGi->qrystartpos, pMGi->qryendpos) < min(pMGj->qrystartpos, pMGj->qryendpos) ) { // MG i starts before MG j on query : typically use normal orientation of query
	  if (max(pMGi->qrystartpos, pMGi->qryendpos) - max(pMGj->qrystartpos, pMGj->qryendpos) > min(pMGj->qrystartpos, pMGj->qryendpos) - min(pMGi->qrystartpos, pMGi->qryendpos)){
	    // NEW80 : exception : MG i ends after MG j on query AND non-overlapped region of MG i on right side of query is larger than on left side : use reversed orientation of query
	    inv_MGi = pMGi->orientforward;
	  } else
	    inv_MGi = !pMGi->orientforward;
	} else {// MG j starts before MG i on query : typically use reversed orientation of query
	  if (max(pMGj->qrystartpos, pMGj->qryendpos) - max(pMGi->qrystartpos, pMGi->qryendpos) > min(pMGi->qrystartpos, pMGi->qryendpos) - min(pMGj->qrystartpos, pMGj->qryendpos)) {
	    // NEW80 : exception : MG j ends after MG i on query AND non-overlapped region on MG j on right side of query is larger than on left side : use normal orientation of query
	    inv_MGi = !pMGi->orientforward;
	  } else
	    inv_MGi = pMGi->orientforward;
	}

	if (pMGi->refstartpos >= pMGj->refstartpos /* && pMGi->refendpos <= pMGj->refendpos */){/* MG i is completely overlapped by MG j*/
	  if(DEBUG/* HERE >=2 */) assert(pMGi->refendpos <= pMGj->refendpos);

	  if(inv_MGi) { /* MG i is inverted : just rule out inversion */
	    NoInv = true;
	  } else { /* MG j is inverted */
	    /* special handling is requred : Treat MG i as inverted matchgroup and check if spacing between implied inverted duplication exceeds smap_InvDup_maxSpacing (excluding first and last label interval in query gap) */
	    double Spacing = 0.0;
	    if(pMGi->orientforward){/* use inverted query to compute Inverted Duplication Spacing */
	      if(DEBUG) assert(pMGi->qrystartpos < pMGi->qryendpos);
	      if(DEBUG) assert(pMGj->qrystartpos > pMGj->qryendpos);
	      if( pMGj->qrystartpos < pMGi->qryendpos ) { // MG i starts before MG j on inverted query
		Spacing = (pMGi->refstartpos - pMGj->refstartpos) * 1e-3;
		if(pMGi->qrystartidx - pMGj->qrystartidx > 2)
		  Spacing += X[pMGi->qrystartidx - 1] - X[pMGj->qrystartidx + 1];
	      } else { // MG j starts before MG i on inverted query
		Spacing = (pMGj->refendpos - pMGi->refendpos) * 1e-3;
		if(pMGj->qryendidx - pMGi->qryendidx > 2)
		  Spacing += X[pMGj->qryendidx - 1] - X[pMGi->qryendidx + 1];
	      }
	    } else {/* use un-inverted query to compute Inverted Duplication Spacing */
	      if(DEBUG) assert(pMGi->qrystartpos > pMGi->qryendpos);
	      if(DEBUG) assert(pMGj->qrystartpos < pMGj->qryendpos);
	      if( pMGj->qrystartpos < pMGi->qryendpos ) { // MG j starts before MG i on non-inverted query
		Spacing = (pMGj->refendpos - pMGi->refendpos) * 1e-3;
		if(pMGi->qryendidx - pMGj->qryendidx > 2)
		  Spacing += X[pMGi->qryendidx - 1] - X[pMGj->qryendidx + 1];
	      } else {// MG i starts before MG j on non-inverted query
		Spacing = (pMGi->refstartpos - pMGj->refstartpos) * 1e-3;
		if(pMGj->qrystartidx - pMGi->qrystartidx > 2)
		  Spacing += X[pMGj->qrystartidx - 1] - X[pMGi->qrystartidx + 1];
	      }
	    }
	    if(Spacing > smap_InvDup_maxSpacing /* NEW */ || qrystopidx - qrystartidx > refoverlapSites * smap_InvDup_maxQryGapMult){
	      if(DEBUG) assert(inv_MGi == false);
	      NoInvDup = true;/* rules out Inverted Duplication */
	    } else
	      inv_MGi = true;
	  }/* MG j is inverted */
	  
	}/* MG i is completely overlapped by MG j */

      } else {/* MG j is completely overlapped by MG i on reference (or both are the same interval), hence cannot be an inversion with MG j inverted */
	if(DEBUG && !(pMGj->refendpos <= pMGi->refendpos)){
	  printf("i=%d,j=%d:refstartpos= %0.1f,%0.1f (refid=%lld,qryid=%lld)\n",i,j,pMGi->refstartpos,pMGj->refstartpos,pMGi->refcontigid,pMGi->qrycontigid);
	  fflush(stdout);
          assert(pMGj->refendpos <= pMGi->refendpos);
	}

	/* Treat MG j as inverted matchgroup and check if spacing between implied inverted duplications exceeds smap_InvDup_maxSpacing (but don't include first and last label interval of query gap in spacing, since it might be part of the duplication region that is unlabelled). 
	   If so rule out inverted duplication by setting inv_MGi = NoInvDup = true UNLESS inversion is also not possible (more than half of MGi would have to be trimmed away)
	 */
	
	double Spacing = 0.0;/* inverted duplication spacing with MGj as inverted matchgroup */
	double TrimmedMGi = 0.0;/* size of trimmed MGi with MGi as inverted matchgroup */
	if(pMGj->orientforward){/* Use inverted query to compute Inverted Duplication Spacing */
	  if(DEBUG) assert(pMGj->qrystartpos < pMGj->qryendpos);
	  if(DEBUG) assert(pMGi->qrystartpos > pMGi->qryendpos);
	  if( pMGi->qrystartpos < pMGj->qryendpos ) { // MG j starts before MG i on inverted query
	    TrimmedMGi = Spacing = (pMGj->refstartpos - pMGi->refstartpos) * 1e-3;
	    if(pMGj->qrystartidx - pMGi->qrystartidx > 2)
	      Spacing += X[pMGj->qrystartidx - 1] - X[pMGi->qrystartidx + 1];
	  } else { // MG i starts before MG j on inverted query
	    TrimmedMGi = Spacing = (pMGi->refendpos - pMGj->refendpos) * 1e-3;
	    if(pMGi->qryendidx - pMGj->qryendidx > 2)
	      Spacing += X[pMGi->qryendidx - 1] - X[pMGj->qryendidx + 1];
	  }
	} else {/* Use un-inverted query to compute Inverted Duplication Spacing */
	  if(DEBUG && !(pMGj->qrystartpos > pMGj->qryendpos)){
	    printf("i=%d,j=%d:xmapid=%d,%d, qryid=%lld,refid=%lld: qryi= %0.3f .. %0.3f, qryj= %0.3f .. %0.3f, refi= %0.3f .. %0.3f kb, refj= %0.3f .. %0.3f kb\n",
		   i,j,pMGi->xmapid,pMGj->xmapid, pMGi->qrycontigid, pMGi->refcontigid, 
		   pMGi->qrystartpos*1e-3, pMGi->qryendpos*1e-3, pMGj->qrystartpos*1e-3, pMGj->qryendpos*1e-3,
		   pMGi->refstartpos*1e-3, pMGi->refendpos*1e-3, pMGj->refstartpos*1e-3, pMGj->refendpos*1e-3);
	    fflush(stdout);
	    assert(pMGj->qrystartpos > pMGj->qryendpos);
	  }
	  if(DEBUG) assert(pMGi->qrystartpos < pMGi->qryendpos);
	  if( pMGi->qrystartpos < pMGj->qryendpos ) { // MG i starts before inverted MG j on non-inverted query
	    TrimmedMGi = Spacing = (pMGi->refendpos - pMGj->refendpos) * 1e-3;
	    if(pMGj->qryendidx - pMGi->qryendidx > 2)
	      Spacing += X[pMGj->qryendidx - 1] - X[pMGi->qryendidx + 1];
	  } else {// MG j starts before MG i on non-inverted query
	    TrimmedMGi = Spacing = (pMGj->refstartpos - pMGi->refstartpos) * 1e-3;
	    if(pMGi->qrystartidx - pMGj->qrystartidx > 2)
	      Spacing += X[pMGi->qrystartidx - 1] - X[pMGj->qrystartidx + 1];
	  }
	}
	
	if(((Spacing > smap_InvDup_maxSpacing /* NEW */ || qrystopidx - qrystartidx > refoverlapSites * smap_InvDup_maxQryGapMult) &&
	    /*NEW50 */ TrimmedMGi >= (pMGi->refendpos - pMGi->refstartpos) * 0.001 * (1.0 - smap_Inv_MaxTrimFrac)) ||
	   TrimmedMGi >= (pMGi->refendpos - pMGi->refstartpos) * 0.001 * (1.0 - smap_Inv_MinTrimFrac)){
	  inv_MGi = true;
	  NoInvDup = true;/* rules out Inverted Duplication */
	} else
	  inv_MGi = false;

	if(VERB>=2){
	  printf("i=%d,j=%d:xmapid=%d,%d, qryid=%lld,refid=%lld: qryi= %0.3f .. %0.3f, qryj= %0.3f .. %0.3f, refi= %0.3f .. %0.3f kb, refj= %0.3f .. %0.3f kb\n",
		 i,j,pMGi->xmapid,pMGj->xmapid, pMGi->qrycontigid, pMGi->refcontigid, 
		 pMGi->qrystartpos*1e-3, pMGi->qryendpos*1e-3, pMGj->qrystartpos*1e-3, pMGj->qryendpos*1e-3,
		 pMGi->refstartpos*1e-3, pMGi->refendpos*1e-3, pMGj->refstartpos*1e-3, pMGj->refendpos*1e-3);
	  printf("\tSpacing= %0.3f, TrimmedMGi= %0.3f (smap_InvDup_maxSpacing= %0.3f, smap_Inv_TrimFrac= %0.3f .. %0.3f): inv_MGi=%d,NoInvDup=%d\n", 
		 Spacing, TrimmedMGi, smap_InvDup_maxSpacing, smap_Inv_MinTrimFrac, smap_Inv_MaxTrimFrac, inv_MGi, NoInvDup);
	  fflush(stdout);
	}
      }

      if(inv_MGi){
	inv_refalignlen = (pMGi->refendpos - pMGi->refstartpos) * 1e-3 ;
	inv_refalignSites = pMGi->refendidx - pMGi->refstartidx;
	noninv_refalignSites = pMGj->refendidx - pMGj->refstartidx;

	inv_qryalignlen = fabs(pMGi->qryendpos - pMGi->qrystartpos) * 1e-3;
	inv_qryalignSites = abs(pMGi->qryendidx - pMGi->qrystartidx);
      } else {
	inv_refalignlen = (pMGj->refendpos - pMGj->refstartpos) * 1e-3;/* kb */
	inv_refalignSites = pMGj->refendidx - pMGj->refstartidx;
	noninv_refalignSites = pMGi->refendidx - pMGi->refstartidx;

	inv_qryalignlen = fabs(pMGj->qryendpos - pMGj->qrystartpos) * 1e-3; /* kb */
	inv_qryalignSites = abs(pMGj->qryendidx - pMGj->qrystartidx);
      }

      /* NEW81 : Heuristic : disable inversion if inverted matchgroup is 5x larger than non-inverted matchgroup */
      if(inv_MGi ? (pMGi->refendpos - pMGi->refstartpos > (pMGj->refendpos - pMGj->refstartpos) * 5.0) : 
	 (pMGj->refendpos - pMGj->refstartpos > (pMGi->refendpos - pMGi->refstartpos)*5.0)){
	if((inv_MGi ? pMGj->confidence : pMGi->confidence) < max(InversionNormalLogPv,InversionNormalLogPv5)){
	  if(VERB/* HERE HERE >=2 */){
	    printf("i=%d,j=%d:xmapid=%d,%d:qryid=%lld,refid=%lld,fwd=%d,%d:Suppressing inversion since inverted matchgroup is 5x larger than non-inverted matchgroup and non-inv conf= %0.2f/%0.2f\n",
		   i,j,pMGi->xmapid,pMGj->xmapid,pMGi->qrycontigid,pMGi->refcontigid,pMGi->orientforward,pMGj->orientforward, inv_MGi ? pMGj->confidence : pMGi->confidence, 
		   max(InversionNormalLogPv,InversionNormalLogPv5));
	    fflush(stdout);
	  }
	  NoInv = true;
	}
      }
    }

    qrysize = (qrystop - qrystart) * 1e-3; //convert to kb
    sv_type_t svtype_enum = (qrysize > refsize ? insertion : deletion);// default to indel for each pair of MGs, unless something more specific is detected

    if(VERB>=2){
      printf("SV main loop indp=%d/%d : qrysize= %0.4f, refsize= %0.4f kb\n",indp,indiceslen,qrysize, refsize);
      fflush(stdout);
    }


    //will apply Heng's translocation criteria to both inter-chromosomal and intra-chromosomal events
    // the former are 'rearrangement', the latter are 'intra_rearr'
    bool intra_rearr = !rearrangement && ( qrysize > smap_sv_sizekb || refsize > smap_sv_sizekb ) && /* NEW30 */ fabs(qrysize - refsize) > smap_sv_sizekb;
    if(opporient && intra_rearr && smap_inv_MaxGap > 0)
      intra_rearr = (qrysize > smap_inv_MaxGap || refsize > smap_inv_MaxGap);

    //need minqryalignlen for translocation overlap, and need that and minrefalignlen for overlapped indels
    double minqryalignlen = min(fabs(pMGi->qrystartpos - pMGi->qryendpos),
				fabs(pMGj->qrystartpos - pMGj->qryendpos)) * 1e-3; //kb!
    double minrefalignlen = min(fabs(pMGi->refstartpos - pMGi->refendpos),
				fabs(pMGj->refstartpos - pMGj->refendpos)) * 1e-3; //kb!
    int minrefalignSites = min(abs(pMGi->refstartidx - pMGi->refendidx),
			       abs(pMGj->refstartidx - pMGj->refendidx));
 
    int LeftSites = min(iminq1,iminq2) - 1;/* Number of unaligned labels on left end of query, taking into account only MG's i & j */
    int RightSites = M - max(imaxq1,imaxq2);/* Number of unaligned labels on right end of query, taking into account only MG's i & j */

    // HERE HERE : should check if qry contig ends are fragile sites : If so skip smap_DupSplit_MinQryLen and smap_DupSplit_MinQrySites tests
    bool LeftExtend = false, RightExtend = false;/* If left/right end of query is likely to have been truncated by -extsplit and hence can be extrapolated (only computed if opporient==0) */
    if(!opporient && !rearrangement){
      if(LeftSites <= smap_DupSplit_maxend && (minq1 < minq2 ? ((maxq1-minq1)*0.001 >= smap_DupSplit_MinQryLen && imaxq1-iminq1+1 >= smap_DupSplit_MinQrySites)
					      : ((maxq2-minq2)*0.001 >= smap_DupSplit_MinQryLen && imaxq2-iminq2+1 >= smap_DupSplit_MinQrySites)))
	LeftExtend = true;
      if(RightSites <= smap_DupSplit_maxend && (maxq1 > maxq2 ? ((maxq1-minq1)*0.001 >= smap_DupSplit_MinQryLen && imaxq1-iminq1+1 >= smap_DupSplit_MinQrySites)
					       : ((maxq2-minq2)*0.001 >= smap_DupSplit_MinQryLen && imaxq2-iminq2+1 >= smap_DupSplit_MinQrySites)))
	RightExtend = true;
    }
    if(VERB>=2){
      printf("iminq1=%d,imaxq1=%d,iminq2=%d,imaxq2=%d,minq1=%0.3f,maxq1=%0.3f,minq2=%0.3f,maxq2=%0.3f:LeftSites=%d,RightSites=%d,LeftExtend=%dd,RightExtend=%d\n",
	     iminq1,imaxq1,iminq2,imaxq2,minq1*0.001,maxq1*0.001,minq2*0.001,maxq2*0.001,LeftSites,RightSites,LeftExtend,RightExtend);
      printf("\tsmap_DupSplit_MinQryLen= %0.3f, smap_DupSplit_MinQrySites=%d, smap_DupSplit_maxend= %d\n",smap_DupSplit_MinQryLen, smap_DupSplit_MinQrySites,smap_DupSplit_maxend);
      fflush(stdout); 
    }

    if(VERB/* HERE HERE >=2 */ && verbose){
      printf("i=%d,j=%d: xmapids=%d,%d:qrycontig=%lld,%lld refid=%lld,%lld fwd=%d,%d, conf=%0.2f,%0.2f: qryoverlap=%d(Sites=%d,size=%0.3f),refoverlap=%d(Sites=%d,size=%0.3f)\n",
	     i,j, pMGi->xmapid, pMGj->xmapid,pMGi->qrycontigid,pMGj->qrycontigid,pMGi->refcontigid,
	     pMGj->refcontigid,pMGi->orientforward,pMGj->orientforward,
	     pMGi->confidence, pMGj->confidence, qryoverlap, qryoverlapSites, qryoverlapsize, refoverlap,refoverlapSites,refoverlapsize);
      printf("\t gap qryidx=%d..%d,qry=%0.3f .. %0.3f : ri=%d..%d(%0.3f..%0.3f),rj=%d..%d(%0.3f..%0.3f)\n\t LeftSites=%d,RightSites=%d,qi=%d..%d(%0.3f..%0.3f),qj=%d..%d(%0.3f..%0.3f),qrysize= %0.3f,Extend=%d,%d,rearr=%d,intra-rearr=%d,opporient=%d\n",
	     qrystartidx, qrystopidx, qrystart * 0.001, qrystop *0.001, 
	     pMGi->refstartidx, pMGi->refendidx, pMGi->refstartpos * 0.001, pMGi->refendpos * 0.001, 
	     pMGj->refstartidx, pMGj->refendidx, pMGj->refstartpos * 0.001, pMGj->refendpos * 0.001,
	     LeftSites, RightSites, iminq1, imaxq1, minq1 * 0.001, maxq1 * 0.001, iminq2, imaxq2, minq2 * 0.001,maxq2 * 0.001, qrysize,
	     LeftExtend,RightExtend, rearrangement,intra_rearr, opporient);
      printf("\t xmapid_overlapped= %d,%d, smallDupInv= %d,%d, inversion_only=%d(%d,%d)\n",
	     pMGi->xmapid_overlapped, pMGj->xmapid_overlapped, pMGi->smallDupInv, pMGj->smallDupInv, inversion_only,pMGi->inversion_only, pMGj->inversion_only);
      fflush(stdout);
    }

    if( DEBUG && !rearrangement && (refoverlap != (refsize < 0)) ) {
      printf("SV main loop indp=%d/%d: qid=%4lld, rearrangement=%d, qryoverlap=%d, qrysize=%0.4f, qryoverlapsize=%.3f, refsize=%0.4f, refoverlap=%d/%.3f\n", indp, indiceslen, contigid, rearrangement, qryoverlap, qrysize, qryoverlapsize, refsize, refoverlap, refoverlapsize);
      printf(" usexmapentries[i=%d]: orientforward=%d, qryendpos=%0.3f, qrystartpos=%0.3f\n",i,pMGi->orientforward ? 1 : 0, pMGi->qryendpos, pMGi->qrystartpos);
      printf(" usexmapentries[j=%d]: orientforward=%d, qryendpos=%0.3f, qrystartpos=%0.3f\n",j,pMGj->orientforward ? 1 : 0, pMGj->qryendpos, pMGj->qrystartpos);
      fflush(stdout);
      assert(0);
    }

    bool translocation = false; //translocation is a rearrangement that passes Heng's criteria
    bool overlap_indel = false, duplication = false, duplication_split = false;
    bool complex_sv = false; //bug in intel compiler creates problems using 'complex' as variable nam

    double refstart = 0.0, refstop = 0.0;
    int refstartidx = 0, refstopidx = 0;

    if(MATCHGROUP_TRIM && XmapUnique && pMGi->xmapid_overlapped == pMGj->xmapid && qryoverlapSites < imaxq1 - iminq1){
      if(VERB>=2){
	printf(" i=%d,j=%d (xmapids=%d,%d) : pMGi->xmapid_overlapped = %d -> -1 : qryoverlapSites=%d, qryi= %d..%d, qryj= %d..%d\n",
	       i,j,pMGi->xmapid,pMGj->xmapid, pMGi->xmapid_overlapped,qryoverlapSites,iminq1,imaxq1,iminq2,imaxq2);
	fflush(stdout);
      }
      pMGi->xmapid_overlapped = -1;/* MGj must have been trimmed since analyzing overlap with MGi */
    }

    if(MATCHGROUP_TRIM && XmapUnique && pMGj->xmapid_overlapped == pMGi->xmapid && qryoverlapSites < imaxq2 - iminq2){
      if(VERB>=2){
	printf(" i=%d,j=%d (xmapids=%d,%d) : pMGj->xmapid_overlapped = %d -> -1 : qryoverlapSites=%d, qryi= %d..%d, qryj= %d..%d\n",
	       i,j,pMGi->xmapid,pMGj->xmapid, pMGj->xmapid_overlapped,qryoverlapSites,iminq1,imaxq1,iminq2,imaxq2);
	fflush(stdout);
      }
      pMGj->xmapid_overlapped = -1;/* MGi must have been trimmed since analyzing overlap with MGj */
    }

    bool qryiFullOverlapped = false;
    bool qryjFullOverlapped = false;
    if(MATCHGROUP_TRIM <= 2){// NEW162
      double mini = min(pMGi->qrystartpos, pMGi->qryendpos);
      double minj = min(pMGj->qrystartpos, pMGj->qryendpos);
      double maxi = max(pMGi->qrystartpos, pMGi->qryendpos);
      double maxj = max(pMGj->qrystartpos, pMGj->qryendpos);
      double overlap = min(maxi,maxj) - max(mini,minj);
      double minlen = min(maxi-mini, maxj-minj);
      
      if(overlap >= 0.999 * minlen){
	if(maxi-mini > maxj-minj)
	  qryjFullOverlapped = true;

	if(maxi-mini < maxj-minj)
	  qryiFullOverlapped = true;
      }
    }
      
    if(MATCHGROUP_TRIM && XmapUnique && (pMGi->xmapid_overlapped > 0 || (MATCHGROUP_TRIM <= 2 && qryiFullOverlapped))
       && !pMGi->smallDupInv && !inversion_only){/* check for translocations */
      if(DEBUG>=1+RELEASE && MATCHGROUP_TRIM >= 3) assert(pMGi->xmapid_overlapped == pMGj->xmapid);

      // NOTE : MGi must be left of MGj on reference, but can be fully overlapped on reference if left ends coincide
      if(pMGi->refcontigid == pMGj->refcontigid && pMGi->refstartpos >= pMGj->refstartpos){// call duplication, provided gap is smaller than smap_Dup_maxSpacing (or smap_InvDup_maxSpacing)
	if((DEBUG && !(pMGi->refstartpos == pMGj->refstartpos)) || (VERB>=2)){
	  printf("i=%d,j=%d:xmapid=%d,%d,qryid=%lld,%lld,refid=%lld,%lld: refi= %0.3f .. %0.3f, refj= %0.3f .. %0.3f, qryi= %0.3f .. %0.3f, qryj = %0.3f .. %0.3f kb\n",
		 i,j,pMGi->xmapid,pMGj->xmapid, pMGi->qrycontigid, pMGj->qrycontigid, pMGi->refcontigid, pMGj->refcontigid,
		 pMGi->refstartpos * 1e-3, pMGi->refendpos * 1e-3, pMGj->refstartpos * 1e-3, pMGj->refendpos * 1e-3,
		 pMGi->qrystartpos * 1e-3, pMGi->qryendpos * 1e-3, pMGj->qrystartpos * 1e-3, pMGj->qryendpos * 1e-3);
	  fflush(stdout);
		 
	  assert(pMGi->refstartpos == pMGj->refstartpos);
	}
	if(DEBUG) assert(pMGi->refendpos <= pMGj->refendpos);

	// query orientation will be based on larger MGj 
	if(pMGi->orientforward != pMGj->orientforward){/* check if Inverted-Duplication spacing is greater than smap_InvDup_maxSpacing */
	  //	  double Spacing = fabs(pMGj->qrystartpos - pMGi->qryendpos);
	  // HERE HERE  	  
	} else {/* check if Duplication spacing is greater than smap_Dup_maxSpacing */
	  //	  double Spacing = fabs(pMGj->qrystartpos - pMGi->qrystartpos);
	  // HERE HERE  	  
	}

	// HERE HERE 
	if(VERB>=1+RELEASE){
	  printf("Possible unhandled Duplication\n");
	  fflush(stdout); // exit(1);
	}
	continue;
      }

      // call translocation with MGi overlapping Insertion in MGj 
      if(VERB>=1+RELEASE/* HERE HERE >=2 */){
	printf("i=%d,j=%d:xmapid=%d,%d,qryid=%lld,%lld,refid=%lld,%lld,fwd=%d,%d: refi= %0.3f .. %0.3f, refj= %0.3f .. %0.3f, qryi= %0.3f .. %0.3f, qryj = %0.3f .. %0.3f kb\n",
	       i,j,pMGi->xmapid,pMGj->xmapid, pMGi->qrycontigid, pMGj->qrycontigid, pMGi->refcontigid, pMGj->refcontigid,
	       pMGi->orientforward,pMGj->orientforward,
	       pMGi->refstartpos * 1e-3, pMGi->refendpos * 1e-3, pMGj->refstartpos * 1e-3, pMGj->refendpos * 1e-3,
	       pMGi->qrystartpos * 1e-3, pMGi->qryendpos * 1e-3, pMGj->qrystartpos * 1e-3, pMGj->qryendpos * 1e-3);
	fflush(stdout);
      }

      if(pMGi->refcontigid == pMGj->refcontigid){/* Case A: intra-chromosomal translocation : may have partial reference overlap */
	if(VERB>=2){
	  printf("\t Possible Overlap intra-Translocation\n");
	  fflush(stdout);
	}

	/* map MGi query coordinates to reference via MGj */
	int QryLeftMGi, QryRightMGi;/* Left and Right ends of MGi on query with query orientated as in MGj */
	double QryLeftLocMGi, QryRightLocMGi;
	if(pMGj->orientforward){
	  QryLeftMGi = min(pMGi->qrystartidx, pMGi->qryendidx);
	  QryRightMGi = max(pMGi->qrystartidx, pMGi->qryendidx);
	  QryLeftLocMGi = min(pMGi->qrystartpos, pMGi->qryendpos);
	  QryRightLocMGi = max(pMGi->qrystartpos, pMGi->qryendpos);
	} else {
	  QryLeftMGi = max(pMGi->qrystartidx, pMGi->qryendidx);
	  QryRightMGi = min(pMGi->qrystartidx, pMGi->qryendidx);
	  QryLeftLocMGi = max(pMGi->qrystartpos, pMGi->qryendpos);
	  QryRightLocMGi = min(pMGi->qrystartpos, pMGi->qryendpos);
	}

	int QryLeftMGj = QryLeftMGi, QryRightMGj = QryRightMGi, RefLeftMGj, RefRightMGj;/* coordinates of boundary of insertion (with query oriented as in MGj) */
	double QryLeftLocMGj, QryRightLocMGj, RefLeftLocMGj, RefRightLocMGj;/* corresponding locations in bp */
	int Lindex = getRef(pMGj->align, QryLeftMGj, RefLeftLocMGj, RefLeftMGj, pMGj->orientforward ? -1 : 1, &QryLeftLocMGj);
	int Rindex = getRef(pMGj->align, QryRightMGj, RefRightLocMGj, RefRightMGj, pMGj->orientforward ? 1 : -1, &QryRightLocMGj);

	double ConfLeftLogPV = SubrangeAlignFP(0,Lindex,pMGj->align,YYmap[pMGj->align->mapid1],Xmap);
	double ConfRightLogPV = SubrangeAlignFP(Rindex,pMGj->align->numpairs - 1,pMGj->align,YYmap[pMGj->align->mapid1],Xmap);

	if(DEBUG && !(RefLeftLocMGj > pMGi->refstartpos)){
	  printf("QryLeftMGj= %d, QryRightMGj= %d: QryLeftLocMGj= %0.1f, QryRightLocMGj= %0.1f,RefLeftLocMGj= %0.1f, RefRightLocMGj= %0.1f,pMGi->refstartpos= %0.1f\n",
		 QryLeftMGj,QryRightMGj,QryLeftLocMGj,QryRightLocMGj,RefLeftLocMGj,RefRightLocMGj,pMGi->refstartpos);
	  fflush(stdout);
	  assert(RefLeftLocMGj > pMGi->refstartpos);
	}

#if 0 // WAS405
	/* NOTE : the following code requires intra-translocations to span 5Mb : is this needed ? */
	double trans_sizekb = (RefLeftLocMGj - pMGi->refstartpos) * 1e-3;
	if(trans_sizekb < smap_sv_sizekb){/* mark as complex */
	  if(VERB>=1+RELEASE){
	    printf("RefLeftLocMGj= %0.4f, pMGi->refendpos = %0.4f: trans_sizekb= %0.4f, smap_sv_sizekb= %0.4f: marking as complex\n",
		   RefLeftLocMGj*1e-3, pMGi->refendpos*1e-3,trans_sizekb,smap_sv_sizekb);
	    fflush(stdout);
	  }
	  complex_sv = true;
	  continue;
	}
#endif

	translocation = true;
	
	/* compute coordinates for left translocation breakpoint (on ref) */
	qrystart = QryLeftLocMGj;
	qrystop = QryLeftLocMGi;
	qrystartidx = QryLeftMGj;
	qrystopidx = QryLeftMGi;

	refstart = RefLeftLocMGj;
	refstartidx = RefLeftMGj;
	if(pMGj->orientforward == pMGi->orientforward){
	  refstop = pMGi->refstartpos;
	  refstopidx = pMGi->refstartidx;
	} else {
	  refstop = pMGi->refendpos;
	  refstopidx = pMGi->refendidx;
	}
	
	if(verbose){
	  printf("  overlap intra-translocation left breakpoint: refid=%lld,qryid=%lld: qry= %0.3f .. %0.3f, ref= %0.3f .. %0.3f kb, ConfLeft= %0.2f -> %0.2f, numsv=%d\n",
		 pMGj->refcontigid, pMGi->qrycontigid, qrystart,qrystop,refstart,refstop, pMGj->confidence, ConfLeftLogPV,numsv);
	  fflush(stdout);
	}
	
	svtype_enum = intrachr_translocation;
	const char *tyst = "transintra"; //for debugging--more info than 'complex'

	if(verbose && VERB /* HERE HERE >= 2 */){
	  printf("SV main loop indp=%d/%d: xmapids=%d,%d, svtype_enum=%d, tyst=%s, numsv=%d\n",
		 indp,indiceslen,pMGj->xmapid, pMGi->xmapid, svtype_enum, tyst ? tyst : "null", numsv);
	  fflush(stdout);
	}

	if(min(ConfLeftLogPV, pMGi->confidence) >= LogPvThreshold){//create structuralVariation object, store indices
	  if(numsv >= maxsv)
	    maxsv = growArray(sv_array,maxsv, max(numsv * 3 / 2, numsv + 17));
	  sv_array[numsv] = structuralVariation(j, i, qrysize, refsize, qrystart, qrystop, pMGj->refcontigid, pMGi->refcontigid, refstart, refstop, tyst, svtype_enum, qrystartidx, qrystopidx, refstartidx, refstopidx, -1., ConfLeftLogPV, pMGi->confidence, sv_ps, pMGj->xmapid, pMGi->xmapid, true/* trans_overlapped*/);
	  sv_array[numsv].inv_MGi = inv_MGi;
	  sv_array[numsv].MGi_qryleft = MGi_qryleft;
	  sv_array[numsv].repeat = (pMGi->have_ref_repeat || pMGj->have_ref_repeat); //move to method of sv class
	
	  if(SV_ORIENT){
	    if(pMGi->orientforward){
	      if(pMGj->orientforward){
		sv_array[numsv].orientation = strdup("+/+");
	      } else {
		sv_array[numsv].orientation = strdup("+/-");	    
	      }
	    } else {
	      if(pMGj->orientforward){
		sv_array[numsv].orientation = strdup("+/-");	    
	      } else {
		sv_array[numsv].orientation = strdup("-/-");
	      }
	    }
	  }
	  numsv++;
	  if(DEBUG) assert(numsv <= maxsv);
	}

	/* compute coordinate for right translocation breakpoint */
	qrystart = QryRightLocMGi;
	qrystartidx = QryRightMGi;
	qrystop = QryRightLocMGj;
	qrystopidx = QryRightMGj;
	if(pMGj->orientforward == pMGi->orientforward){
	  refstart = pMGi->refendpos;
	  refstartidx = pMGi->refendidx;
	} else {
	  refstart = pMGi->refstartpos;
	  refstartidx = pMGi->refstartidx;
	}
	refstop = RefRightLocMGj;
	refstopidx = RefRightMGj;

	if(verbose){
	  printf("  overlap intra-translocation right breakpoint: refid=%lld,qryid=%lld: qry= %0.3f .. %0.3f, ref= %0.3f .. %0.3f kb, ConfRight= %0.2f -> %0.2f,numsv=%d\n",
		 pMGi->refcontigid, pMGi->qrycontigid, qrystart,qrystop,refstart,refstop,pMGj->confidence,ConfRightLogPV,numsv);
	  fflush(stdout);
	}

	svtype_enum = intrachr_translocation;
	tyst = "transintra"; //for debugging--more info than 'complex'

	if(verbose && VERB /* HERE HERE >= 2 */){
	  printf("SV main loop indp=%d/%d: xmapids=%d,%d, svtype_enum=%d, tyst=%s, numsv=%d\n",
		 indp,indiceslen,pMGi->xmapid, pMGj->xmapid, svtype_enum, tyst ? tyst : "null", numsv);
	  fflush(stdout);
	}

	if(min(ConfRightLogPV, pMGi->confidence) >= LogPvThreshold){//create structuralVariation object, store indices
	  if(numsv >= maxsv)
	    maxsv = growArray(sv_array,maxsv, max(numsv * 3 / 2, numsv + 17));
	  sv_array[numsv] = structuralVariation(i, j, qrysize, refsize, qrystart, qrystop, pMGi->refcontigid, pMGj->refcontigid, refstart, refstop, tyst, svtype_enum, qrystartidx, qrystopidx, refstartidx, refstopidx, -1., pMGi->confidence, ConfRightLogPV, sv_ps, pMGi->xmapid, pMGj->xmapid, true/* trans_overlapped*/);
	  sv_array[numsv].inv_MGi = inv_MGi;
	  sv_array[numsv].MGi_qryleft = MGi_qryleft;
	  sv_array[numsv].repeat = (pMGi->have_ref_repeat || pMGj->have_ref_repeat); //move to method of sv class
	
	  if(SV_ORIENT){
	    if(pMGi->orientforward){
	      if(pMGj->orientforward){
		sv_array[numsv].orientation = strdup("+/+");
	      } else {
		sv_array[numsv].orientation = strdup("-/+");	    
	      }
	    } else {
	      if(pMGj->orientforward){
		sv_array[numsv].orientation = strdup("-/+");	    
	      } else {
		sv_array[numsv].orientation = strdup("-/-");
	      }
	    }
	  }
	  numsv++;
	  if(DEBUG) assert(numsv <= maxsv);
	}

	continue;
      } // if ( pMGi->refcontigid == pMGj->refcontigid)

      /* caseB : inter-chromosomal translocation */
      translocation = true;

      /* map MGi query coordinates to reference via MGj */
      int QryLeftMGi, QryRightMGi;/* Left and Right ends of MGi on query with query orientated as in MGj */
      double QryLeftLocMGi, QryRightLocMGi;
      if(pMGj->orientforward){
	QryLeftMGi = min(pMGi->qrystartidx, pMGi->qryendidx);
	QryLeftLocMGi = min(pMGi->qrystartpos, pMGi->qryendpos);
	QryRightMGi = max(pMGi->qrystartidx, pMGi->qryendidx);
	QryRightLocMGi = max(pMGi->qrystartpos, pMGi->qryendpos);
      } else {
	QryLeftMGi = max(pMGi->qrystartidx, pMGi->qryendidx);
	QryLeftLocMGi = max(pMGi->qrystartpos, pMGi->qryendpos);
	QryRightMGi = min(pMGi->qrystartidx, pMGi->qryendidx);
	QryRightLocMGi = min(pMGi->qrystartpos, pMGi->qryendpos);
      }

      int QryLeftMGj = QryLeftMGi, QryRightMGj = QryRightMGi, RefLeftMGj, RefRightMGj;/* coordinates of boundary of insertion (with query oriented as in MGj) */
      double QryLeftLocMGj, QryRightLocMGj, RefLeftLocMGj, RefRightLocMGj;/* corresponding locations in bp */
      int Lindex = getRef(pMGj->align, QryLeftMGj, RefLeftLocMGj, RefLeftMGj, pMGj->orientforward ? -1 : 1, &QryLeftLocMGj); // (pMGi->refcontigid==2 && pMGj->refcontigid==9 && pMGi->qrycontigid==82)
      int Rindex = getRef(pMGj->align, QryRightMGj, RefRightLocMGj, RefRightMGj, pMGj->orientforward ? 1 : -1, &QryRightLocMGj);

      double ConfLeftLogPV = SubrangeAlignFP(0,Lindex,pMGj->align,YYmap[pMGj->align->mapid1],Xmap);
      double ConfRightLogPV = SubrangeAlignFP(Rindex,pMGj->align->numpairs - 1,pMGj->align,YYmap[pMGj->align->mapid1],Xmap);

      if(VERB>=1+RELEASE/* HERE HERE >=2 */){
	printf("  refid=%lld,%lld,qryid=%lld:qryMGi=%d..%d(%0.3f..%0.3f),refMGi=%d..%d(%0.3f..%0.3f),qryMGj=%d..%d(%0.3f..%0.3f),refMGj=%d..%d(%0.3f..%0.3f)\n",
	       pMGi->refcontigid,pMGj->refcontigid,pMGi->qrycontigid,QryLeftMGi,QryRightMGi,QryLeftLocMGi*1e-3,QryRightLocMGi*1e-3,
	       pMGi->refstartidx,pMGi->refendidx,pMGi->refstartpos*1e-3,pMGi->refendpos*1e-3,
	       QryLeftMGj,QryRightMGj,QryLeftLocMGj*1e-3,QryRightLocMGj*1e-3,
	       RefLeftMGj,RefRightMGj,RefLeftLocMGj*1e-3,RefRightLocMGj*1e-3);
	fflush(stdout);
      }

      /* compute coordinate for left translocation breakpoint (on ref for MGj) */
      qrystart = QryLeftLocMGj;
      qrystartidx = QryLeftMGj;
      qrystop = QryLeftLocMGi;
      qrystopidx = QryLeftMGi;

      refstart = RefLeftLocMGj;
      refstartidx = RefLeftMGj;
      if(pMGi->orientforward == pMGj->orientforward){
	refstop = pMGi->refstartpos;
	refstopidx = pMGi->refstartidx;
      } else {
	refstop = pMGi->refendpos;
	refstopidx = pMGi->refendidx;
      }

      if(verbose){
	printf("  overlap translocation left breakpoint: refid=%lld,%lld,qryid=%lld: qry= %0.3f .. %0.3f, ref= %0.3f (MGj), %0.3f kb (MGi), ConfLeft= %0.2f -> %0.2f,numsv=%d\n",
	       pMGj->refcontigid, pMGi->refcontigid, pMGi->qrycontigid, qrystart*1e-3,qrystop*1e-3,refstart*1e-3,refstop*1e-3,pMGj->confidence, ConfLeftLogPV,numsv);
	fflush(stdout);
      }

      svtype_enum = interchr_translocation;
      const char *tyst = "transinter"; //for debugging--more info than 'complex'

      if(verbose && VERB /* HERE HERE >= 2 */){
	printf("SV main loop indp=%d/%d: xmapids=%d,%d, svtype_enum=%d, tyst=%s, numsv=%d\n",
	       indp,indiceslen,pMGi->xmapid, pMGj->xmapid, svtype_enum, tyst ? tyst : "null", numsv);
	fflush(stdout);
      }

      if(min(ConfLeftLogPV, pMGi->confidence) >= LogPvThreshold){//create structuralVariation object, store indices
	if(numsv >= maxsv)
	  maxsv = growArray(sv_array,maxsv, max(numsv * 3 / 2, numsv + 17));

	sv_array[numsv] = structuralVariation(j, i, qrysize, refsize, qrystart, qrystop, pMGj->refcontigid, pMGi->refcontigid, refstart, refstop, tyst, svtype_enum, qrystartidx, qrystopidx, refstartidx, refstopidx, -1., ConfLeftLogPV, pMGi->confidence, sv_ps, pMGj->xmapid, pMGi->xmapid, true/* trans_overlapped */);
	sv_array[numsv].inv_MGi = inv_MGi;
	sv_array[numsv].MGi_qryleft = MGi_qryleft;
	sv_array[numsv].repeat = (pMGi->have_ref_repeat || pMGj->have_ref_repeat); //move to method of sv class

	if(SV_ORIENT){
	  if(pMGi->orientforward){
	    if(pMGj->orientforward){
	      sv_array[numsv].orientation = strdup("+/+");
	    } else {
	      sv_array[numsv].orientation = strdup("+/-");	    
	    }
	  } else {
	    if(pMGj->orientforward){
	      sv_array[numsv].orientation = strdup("+/-");	    
	    } else {
	      sv_array[numsv].orientation = strdup("-/-");
	    }
	  }
	}
	numsv++;
	if(DEBUG) assert(numsv <= maxsv);
      }

      /* compute coordinate for right translocation breakpoint */
      qrystart = QryRightLocMGi;
      qrystartidx = QryRightMGi;
      qrystop = QryRightLocMGj;
      qrystopidx = QryRightMGj;
      if(pMGi->orientforward == pMGj->orientforward){
	refstart = pMGi->refendpos;
	refstartidx = pMGi->refendidx;
      } else {
	refstart = pMGi->refstartpos;
	refstartidx = pMGi->refstartidx;
      }
      refstop = RefRightLocMGj;
      refstopidx = RefRightMGj;
      
      if(verbose){
	printf("  overlap translocation right breakpoint: refid=%lld,%lld,qryid=%lld: qry= %0.3f .. %0.3f, ref= %0.3f (MGi) .. %0.3f kb (MGj), ConfRight= %0.2f -> %0.2f,numsv=%d\n",
	       pMGi->refcontigid, pMGj->refcontigid, pMGi->qrycontigid, qrystart*1e-3, qrystop*1e-3, refstart*1e-3, refstop*1e-3, pMGj->confidence, ConfRightLogPV,numsv);
	fflush(stdout);
      }

      svtype_enum = interchr_translocation;
      tyst = "transinter"; //for debugging--more info than 'complex'

      if(verbose && VERB /* HERE HERE >= 2 */){
	printf("SV main loop indp=%d/%d: xmapids=%d,%d, svtype_enum=%d, tyst=%s, numsv=%d\n",
	       indp,indiceslen,pMGi->xmapid, pMGj->xmapid, svtype_enum, tyst ? tyst : "null", numsv);
	fflush(stdout);
      }

      if(min(ConfRightLogPV, pMGi->confidence) >= LogPvThreshold){//create structuralVariation object, store indices
	if(numsv >= maxsv)
	  maxsv = growArray(sv_array,maxsv, max(numsv * 3 / 2, numsv + 17));

	sv_array[numsv] = structuralVariation(i, j, qrysize, refsize, qrystart, qrystop, pMGi->refcontigid, pMGj->refcontigid, refstart, refstop, tyst, svtype_enum, qrystartidx, qrystopidx, refstartidx, refstopidx, -1., pMGi->confidence, ConfRightLogPV, sv_ps, pMGi->xmapid, pMGj->xmapid, true /* trans_overlapped */);
	sv_array[numsv].inv_MGi = inv_MGi;
	sv_array[numsv].MGi_qryleft = MGi_qryleft;
	sv_array[numsv].repeat = (pMGi->have_ref_repeat || pMGj->have_ref_repeat); //move to method of sv class

	if(SV_ORIENT){
	  if(pMGi->orientforward){
	    if(pMGj->orientforward){
	      sv_array[numsv].orientation = strdup("+/+");
	    } else {
	      sv_array[numsv].orientation = strdup("-/+");	    
	    }
	  } else {
	    if(pMGj->orientforward){
	      sv_array[numsv].orientation = strdup("-/+");	    
	    } else {
	      sv_array[numsv].orientation = strdup("-/-");
	    }
	  }
	}
	numsv++;
	if(DEBUG) assert(numsv <= maxsv);
      }

      continue;
    }

    if(MATCHGROUP_TRIM && XmapUnique && (pMGj->xmapid_overlapped > 0 || (MATCHGROUP_TRIM <= 2 && qryjFullOverlapped)) 
       && !pMGj->smallDupInv && !inversion_only){/* check for translocations */
      if(DEBUG>=1+RELEASE && MATCHGROUP_TRIM >= 3) assert(pMGj->xmapid_overlapped == pMGi->xmapid);

      // NOTE : MGi must be left of MGj on ref, but can be fully overlapped on ref 
      if(pMGi->refcontigid == pMGj->refcontigid && pMGi->refendpos >= pMGj->refendpos){// call duplication, provided gap is smaller than smap_Dup_maxSpacing (or smap_InvDup_maxSpacing)
	if((DEBUG && !(pMGi->refstartpos <= pMGj->refstartpos)) || (VERB>=2)){
	  printf("i=%d,j=%d:xmapid=%d,%d,qryid=%lld,%lld,refid=%lld,%lld: refi= %0.3f .. %0.3f, refj= %0.3f .. %0.3f, qryi= %0.3f .. %0.3f, qryj = %0.3f .. %0.3f kb\n",
		 i,j,pMGi->xmapid,pMGj->xmapid, pMGi->qrycontigid, pMGj->qrycontigid, pMGi->refcontigid, pMGj->refcontigid,
		 pMGi->refstartpos * 1e-3, pMGi->refendpos * 1e-3, pMGj->refstartpos * 1e-3, pMGj->refendpos * 1e-3,
		 pMGi->qrystartpos * 1e-3, pMGi->qryendpos * 1e-3, pMGj->qrystartpos * 1e-3, pMGj->qryendpos * 1e-3);
	  fflush(stdout);
		 
	  assert(pMGi->refstartpos <= pMGj->refstartpos);
	}
	
	// query orientation will be based on larger MGi
	if(pMGi->orientforward != pMGj->orientforward){/* check if Inverted-Duplication spacing is greater than smap_InvDup_maxSpacing */
	  //	  double Spacing = fabs(pMGi->qrystartpos - pMGj->qryendpos);
	  // HERE HERE  	  
	} else {/* check if Duplication spacing is greater than smap_Dup_maxSpacing */
	  //	  double Spacing = fabs(pMGi->qrystartpos - pMGj->qrystartpos);
	  // HERE HERE  	  
	}

	// HERE HERE  
	if(VERB >=1+RELEASE){
	  printf("Possible unhandled Duplication\n");
	  fflush(stdout);// exit(1);
	}
	continue;
      }
      
      // call translocation with MGj overlapping Insertion in MGi 
      if(VERB>=1+RELEASE/* HERE HERE >=2 */){
	printf("i=%d,j=%d:xmapid=%d,%d,qryid=%lld,%lld,refid=%lld,%lld: refi= %0.3f .. %0.3f, refj= %0.3f .. %0.3f, qryi= %0.3f .. %0.3f, qryj = %0.3f .. %0.3f kb\n",
	       i,j,pMGi->xmapid,pMGj->xmapid, pMGi->qrycontigid, pMGj->qrycontigid, pMGi->refcontigid, pMGj->refcontigid,
	       pMGi->refstartpos * 1e-3, pMGi->refendpos * 1e-3, pMGj->refstartpos * 1e-3, pMGj->refendpos * 1e-3,
	       pMGi->qrystartpos * 1e-3, pMGi->qryendpos * 1e-3, pMGj->qrystartpos * 1e-3, pMGj->qryendpos * 1e-3);
	printf("\t xmapid_overlapped= %d, %d, SmalldupInv=%d, %d, inversion_only= %d\n",pMGi->xmapid_overlapped, pMGj->xmapid_overlapped, pMGi->smallDupInv, pMGj->smallDupInv, inversion_only);
	fflush(stdout);
      }

      if(pMGi->refcontigid == pMGj->refcontigid){/* Case C: intra-chromosomal translocation : may have partial reference overlap */
	if(VERB>=2){
	  printf("\t Possible Overlap intra-Translocation\n");
	  fflush(stdout);
	}

	/* map MGj query coordinates to reference via MGi */
	int QryLeftMGj, QryRightMGj;/* Left and Right ends of MGj on query with query orientated as in MGi */
	double QryLeftLocMGj, QryRightLocMGj;
	if(pMGi->orientforward){
	  QryLeftMGj = min(pMGj->qrystartidx, pMGj->qryendidx);
	  QryRightMGj = max(pMGj->qrystartidx, pMGj->qryendidx);
	  QryLeftLocMGj = min(pMGj->qrystartpos, pMGj->qryendpos);
	  QryRightLocMGj = max(pMGj->qrystartpos, pMGj->qryendpos);
	} else {
	  QryLeftMGj = max(pMGj->qrystartidx, pMGj->qryendidx);
	  QryRightMGj = min(pMGj->qrystartidx, pMGj->qryendidx);
	  QryLeftLocMGj = max(pMGj->qrystartpos, pMGj->qryendpos);
	  QryRightLocMGj = min(pMGj->qrystartpos, pMGj->qryendpos);
	}

	int QryLeftMGi = QryLeftMGj, QryRightMGi = QryRightMGj, RefLeftMGi, RefRightMGi;/* coordinates of boundary of insertion (with query oriented as in MGi) */
	double QryLeftLocMGi, QryRightLocMGi, RefLeftLocMGi, RefRightLocMGi;/* corresponding locations in bp */
	int Lindex = getRef(pMGi->align, QryLeftMGi, RefLeftLocMGi, RefLeftMGi, pMGi->orientforward ? -1 : 1, &QryLeftLocMGi);
	int Rindex = getRef(pMGi->align, QryRightMGi, RefRightLocMGi, RefRightMGi, pMGi->orientforward ? 1 : -1, &QryRightLocMGi);

	double ConfLeftLogPV = SubrangeAlignFP(0,Lindex,pMGi->align,YYmap[pMGi->align->mapid1],Xmap);
	double ConfRightLogPV = SubrangeAlignFP(Rindex,pMGi->align->numpairs - 1,pMGi->align,YYmap[pMGi->align->mapid1],Xmap);

	if(DEBUG && !(RefLeftLocMGi < pMGj->refendpos)){
	  printf("\tQryMGj=%d..%d(%0.3f..%0.3f), RefMGj=%d..%d(%0.3f..%0.3f) : QryMGi=%d..%d(%0.3f..%0.3f), RefMGi=%d..%d(%0.3f..%0.3f)\n",
		 QryLeftMGj,QryRightMGj,QryLeftLocMGj*1e-3,QryRightLocMGj*1e-3, pMGj->refstartidx, pMGj->refendidx, pMGj->refstartpos*1e-3,pMGj->refendpos*1e-3,
		 QryLeftMGi,QryRightMGi,QryLeftLocMGi*1e-3,QryRightLocMGi*1e-3, RefLeftMGi, RefRightMGi, RefLeftLocMGi*1e-3, RefRightLocMGi*1e-3);
	  printf("\t RefLeftLocMGi= %0.3f, pMGj->refendpos= %0.3f kb\n", RefLeftLocMGi * 1e-3, pMGj->refendpos * 1e-3);
	  fflush(stdout);
	  assert(RefLeftLocMGi < pMGj->refendpos);
	}

#if 0 // WAS405
	/* NOTE : the following code requires intra-translocations to span 5Mb : is this needed ? */
	double trans_sizekb = (pMGj->refstartpos - RefLeftLocMGi) * 1e-3;
	if(trans_sizekb < smap_sv_sizekb){/* mark as complex */
	  if(VERB>=1+RELEASE){
	    printf("RefLeftLocMGi= %0.4f, pMGj->refstartpos = %0.4f: trans_sizekb= %0.4f, smap_sv_sizekb= %0.4f: marking as complex\n",
		   RefLeftLocMGi*1e-3, pMGj->refstartpos*1e-3,trans_sizekb,smap_sv_sizekb);
	    fflush(stdout);
	  }
	  complex_sv = true;
	  continue;
	}
#endif

	translocation = true;

	/* compute coordinates for left translocation breakpoint (on ref) */
	qrystart = QryLeftLocMGi;
	qrystop = QryLeftLocMGj;
	qrystartidx = QryLeftMGi;
	qrystopidx = QryLeftMGj;

	refstart = RefLeftLocMGi;
	refstartidx = RefLeftMGi;
	if(pMGj->orientforward == pMGi->orientforward){
	  refstop = pMGj->refstartpos;
	  refstopidx = pMGj->refstartidx;
	} else {
	  refstop = pMGj->refendpos;
	  refstopidx = pMGj->refendidx;
	}
	
	if(verbose){
	  printf("  overlap intra-translocation left breakpoint: refid=%lld,qryid=%lld: qry= %0.3f .. %0.3f, ref= %0.3f .. %0.3f kb, ConfLeft= %0.2f -> %0.2f, numsv=%d\n",
		 pMGi->refcontigid, pMGj->qrycontigid, qrystart,qrystop,refstart,refstop, pMGi->confidence, ConfLeftLogPV, numsv);
	  fflush(stdout);
	}

	svtype_enum = intrachr_translocation;
	const char *tyst = "transintra"; //for debugging--more info than 'complex'

	if(verbose && VERB /* HERE HERE >= 2 */){
	  printf("SV main loop indp=%d/%d: xmapids=%d,%d, svtype_enum=%d, tyst=%s, numsv=%d\n",
		 indp,indiceslen,pMGi->xmapid, pMGj->xmapid, svtype_enum, tyst ? tyst : "null", numsv);
	  fflush(stdout);
	}

	if(min(ConfLeftLogPV, pMGj->confidence) >= LogPvThreshold){//create structuralVariation object, store indices
	  if(numsv >= maxsv)
	    maxsv = growArray(sv_array,maxsv, max(numsv * 3 / 2, numsv + 17));

	  sv_array[numsv] = structuralVariation(i, j, qrysize, refsize, qrystart, qrystop, pMGi->refcontigid, pMGj->refcontigid, refstart, refstop, tyst, svtype_enum, qrystartidx, qrystopidx, refstartidx, refstopidx, -1., ConfLeftLogPV, pMGj->confidence, sv_ps, pMGi->xmapid, pMGj->xmapid, true /* trans_overlapped */);
	  sv_array[numsv].inv_MGi = inv_MGi;
	  sv_array[numsv].MGi_qryleft = MGi_qryleft;
	  sv_array[numsv].repeat = (pMGi->have_ref_repeat || pMGj->have_ref_repeat); //move to method of sv class
	
	  if(SV_ORIENT){
	    if(pMGj->orientforward){
	      if(pMGi->orientforward){
		sv_array[numsv].orientation = strdup("+/+");
	      } else {
		sv_array[numsv].orientation = strdup("+/-");	    
	      }
	    } else {
	      if(pMGi->orientforward){
		sv_array[numsv].orientation = strdup("+/-");	    
	      } else {
		sv_array[numsv].orientation = strdup("-/-");
	      }
	    }
	  }
	  numsv++;
	  if(DEBUG) assert(numsv <= maxsv);
	}

	/* compute coordinate for right translocation breakpoint */
	qrystart = QryRightLocMGj;
	qrystartidx = QryRightMGj;
	qrystop = QryRightLocMGi;
	qrystopidx = QryRightMGi;
	if(pMGj->orientforward == pMGi->orientforward){
	  refstart = pMGj->refendpos;
	  refstartidx = pMGj->refendidx;
	} else {
	  refstart = pMGj->refstartpos;
	  refstartidx = pMGj->refstartidx;
	}
	refstop = RefRightLocMGi;
	refstopidx = RefRightMGi;

	if(verbose){
	  printf("  overlap intra-translocation right breakpoint: refid=%lld,qryid=%lld: qry= %0.3f .. %0.3f, ref= %0.3f .. %0.3f kb, ConfRight= %0.2f -> %0.2f,numsv=%d\n",
		 pMGi->refcontigid, pMGi->qrycontigid, qrystart,qrystop,refstart,refstop, pMGi->confidence, ConfRightLogPV, numsv);
	  fflush(stdout);
	}

	svtype_enum = intrachr_translocation;
	tyst = "transintra"; //for debugging--more info than 'complex'

	if(verbose && VERB /* HERE HERE >= 2 */){
	  printf("SV main loop indp=%d/%d: xmapids=%d,%d, svtype_enum=%d, tyst=%s, numsv=%d\n",
		 indp,indiceslen,pMGj->xmapid, pMGi->xmapid, svtype_enum, tyst ? tyst : "null", numsv);
	  fflush(stdout);
	}

	if(min(ConfRightLogPV, pMGj->confidence) >= LogPvThreshold){//create structuralVariation object, store indices
	  if(numsv >= maxsv)
	    maxsv = growArray(sv_array,maxsv, max(numsv * 3 / 2, numsv + 17));

	  sv_array[numsv] = structuralVariation(j, i, qrysize, refsize, qrystart, qrystop, pMGj->refcontigid, pMGi->refcontigid, refstart, refstop, tyst, svtype_enum, qrystartidx, qrystopidx, refstartidx, refstopidx, -1., pMGj->confidence, ConfRightLogPV, sv_ps, pMGj->xmapid, pMGi->xmapid, true /* trans_overlapped */);
	  sv_array[numsv].inv_MGi = inv_MGi;
	  sv_array[numsv].MGi_qryleft = MGi_qryleft;
	  sv_array[numsv].repeat = (pMGi->have_ref_repeat || pMGj->have_ref_repeat); //move to method of sv class
	  
	  if(SV_ORIENT){
	    if(pMGj->orientforward){
	      if(pMGi->orientforward){
		sv_array[numsv].orientation = strdup("+/+");
	      } else {
		sv_array[numsv].orientation = strdup("-/+");	    
	      }
	    } else {
	      if(pMGi->orientforward){
		sv_array[numsv].orientation = strdup("-/+");	    
	      } else {
		sv_array[numsv].orientation = strdup("-/-");
	      }
	    }
	  }
	  numsv++;
	  if(DEBUG) assert(numsv <= maxsv);
	}

	continue;
      } // if(pMGi->refcontigid == pMGj->refcontigid)

      /* case D: inter chromosomal translocation */
      translocation = true;

      /* map MGj query coordinates to reference via MGi */
      int QryLeftMGj, QryRightMGj;/* Left and Right ends of MGj on query with query orientated as in MGi */
      double QryLeftLocMGj, QryRightLocMGj;
      if(pMGi->orientforward){
	QryLeftMGj = min(pMGj->qrystartidx, pMGj->qryendidx);
	QryLeftLocMGj = min(pMGj->qrystartpos, pMGj->qryendpos);
	QryRightMGj = max(pMGj->qrystartidx, pMGj->qryendidx);
	QryRightLocMGj = max(pMGj->qrystartpos, pMGj->qryendpos);
      } else {
	QryLeftMGj = max(pMGj->qrystartidx, pMGj->qryendidx);
	QryLeftLocMGj = max(pMGj->qrystartpos, pMGj->qryendpos);
	QryRightMGj = min(pMGj->qrystartidx, pMGj->qryendidx);
	QryRightLocMGj = min(pMGj->qrystartpos, pMGj->qryendpos);
      }
      int QryLeftMGi = QryLeftMGj, QryRightMGi = QryRightMGj, RefLeftMGi, RefRightMGi;/* coordinates of boundary of insertion (with query oriented as in MGj) */
      double QryLeftLocMGi, QryRightLocMGi, RefLeftLocMGi, RefRightLocMGi;/* corresponding locations in bp */
      int Lindex = getRef(pMGi->align, QryLeftMGi, RefLeftLocMGi, RefLeftMGi, pMGi->orientforward ? -1 : 1, &QryLeftLocMGi);
      int Rindex = getRef(pMGi->align, QryRightMGi, RefRightLocMGi, RefRightMGi, pMGi->orientforward ? 1 : -1, &QryRightLocMGi);

      double ConfLeftLogPV = SubrangeAlignFP(0,Lindex,pMGi->align,YYmap[pMGi->align->mapid1],Xmap);
      double ConfRightLogPV = SubrangeAlignFP(Rindex,pMGi->align->numpairs - 1,pMGi->align,YYmap[pMGi->align->mapid1],Xmap);

      if(VERB>=1+RELEASE/* HERE HERE >=2 */){
	printf("  refid=%lld,%lld,qryid=%lld:qryMGj=%d..%d(%0.3f..%0.3f),refMGj=%d..%d(%0.3f..%0.3f),qryMGi=%d..%d(%0.3f..%0.3f),refMGi=%d..%d(%0.3f..%0.3f)\n",
	       pMGj->refcontigid,pMGi->refcontigid,pMGj->qrycontigid,QryLeftMGj,QryRightMGj,QryLeftLocMGj*1e-3,QryRightLocMGj*1e-3,
	       pMGj->refstartidx,pMGj->refendidx,pMGj->refstartpos*1e-3,pMGj->refendpos*1e-3,
	       QryLeftMGi,QryRightMGi,QryLeftLocMGi*1e-3,QryRightLocMGi*1e-3,
	       RefLeftMGi,RefRightMGi,RefLeftLocMGi*1e-3,RefRightLocMGi*1e-3);
	fflush(stdout);
      }

      /* compute coordinate for left translocation breakpoint (on ref for MGi) */
      qrystart = QryLeftLocMGi;
      qrystartidx = QryLeftMGi;
      qrystop = QryLeftLocMGj;
      qrystopidx = QryLeftMGj;
      refstart = RefLeftLocMGi;
      refstartidx = RefLeftMGi;
      if(pMGi->orientforward == pMGj->orientforward){
	refstop = pMGj->refstartpos;
	refstopidx = pMGj->refstartidx;
      } else {
	refstop = pMGj->refendpos;
	refstopidx = pMGj->refendidx;
      }

      if(verbose){
	printf("  overlap translocation left breakpoint: refid=%lld,%lld,qryid=%lld: qry= %0.3f .. %0.3f, ref= %0.3f (MGi), %0.3f kb (MGj), ConfLeft= %0.2f -> %0.2f, numsv=%d\n",
	       pMGi->refcontigid, pMGj->refcontigid, pMGj->qrycontigid, qrystart*1e-3,qrystop*1e-3,refstart*1e-3,refstop*1e-3, pMGi->confidence, ConfLeftLogPV, numsv);
	fflush(stdout);
      }

      svtype_enum = interchr_translocation;
      const char *tyst = "transinter"; //for debugging--more info than 'complex'

      if(verbose && VERB /* HERE HERE >= 2 */){
	printf("SV main loop indp=%d/%d: xmapids=%d,%d, svtype_enum=%d, tyst=%s, numsv=%d\n",
	       indp,indiceslen,pMGj->xmapid, pMGi->xmapid, svtype_enum, tyst ? tyst : "null", numsv);
	fflush(stdout);
      }

      if(min(ConfLeftLogPV, pMGj->confidence) >= LogPvThreshold){//create structuralVariation object, store indices
	if(numsv >= maxsv)
	  maxsv = growArray(sv_array,maxsv, max(numsv * 3 / 2, numsv + 17));

	sv_array[numsv] = structuralVariation(i, j, qrysize, refsize, qrystart, qrystop, pMGi->refcontigid, pMGj->refcontigid, refstart, refstop, tyst, svtype_enum, qrystartidx, qrystopidx, refstartidx, refstopidx, -1., ConfLeftLogPV, pMGj->confidence, sv_ps, pMGi->xmapid, pMGj->xmapid, true /* trans_overlapped */);
	sv_array[numsv].inv_MGi = inv_MGi;
	sv_array[numsv].MGi_qryleft = MGi_qryleft;
	sv_array[numsv].repeat = (pMGi->have_ref_repeat || pMGj->have_ref_repeat); //move to method of sv class

	if(SV_ORIENT){
	  if(pMGj->orientforward){
	    if(pMGi->orientforward){
	      sv_array[numsv].orientation = strdup("+/+");
	    } else {
	      sv_array[numsv].orientation = strdup("+/-");	    
	    }
	  } else {
	    if(pMGi->orientforward){
	      sv_array[numsv].orientation = strdup("+/-");	    
	    } else {
	      sv_array[numsv].orientation = strdup("-/-");
	    }
	  }
	}
	numsv++;
	if(DEBUG) assert(numsv <= maxsv);
      }

      /* compute coordinate for right translocation breakpoint */
      qrystart = QryRightLocMGj;
      qrystartidx = QryRightMGj;
      qrystop = QryRightLocMGi;
      qrystopidx = QryRightMGi;
      if(pMGi->orientforward == pMGj->orientforward){
	refstart = pMGj->refendpos;
	refstartidx = pMGj->refendidx;
      } else {
	refstart = pMGj->refstartpos;
	refstartidx = pMGj->refstartidx;
      }
      refstop = RefRightLocMGi;
      refstopidx = RefRightMGi;
      
      if(verbose){
	printf("  overlap translocation right breakpoint: refid=%lld,%lld,qryid=%lld: qry= %0.3f .. %0.3f, ref= %0.3f (MGj) .. %0.3f kb (MGi), ConfRight= %0.2f -> %0.2f, numsv=%d\n",
	       pMGj->refcontigid, pMGi->refcontigid, pMGj->qrycontigid, qrystart*1e-3, qrystop*1e-3, refstart*1e-3, refstop*1e-3, pMGi->confidence, ConfRightLogPV, numsv);
	fflush(stdout);
      }

      svtype_enum = interchr_translocation;
      tyst = "transinter"; //for debugging--more info than 'complex'

      if(verbose && VERB /* HERE HERE >= 2 */){
	printf("SV main loop indp=%d/%d: xmapids=%d,%d, svtype_enum=%d, tyst=%s, numsv=%d\n",
	       indp,indiceslen,pMGj->xmapid, pMGi->xmapid, svtype_enum, tyst ? tyst : "null", numsv);
	fflush(stdout);
      }

      //create structuralVariation object, store indices
      if(min(ConfRightLogPV, pMGj->confidence) >= LogPvThreshold){//create structuralVariation object, store indices
	if(numsv >= maxsv)
	  maxsv = growArray(sv_array,maxsv, max(numsv * 3 / 2, numsv + 17));

	sv_array[numsv] = structuralVariation(j, i, qrysize, refsize, qrystart, qrystop, pMGj->refcontigid, pMGi->refcontigid, refstart, refstop, tyst, svtype_enum, qrystartidx, qrystopidx, refstartidx, refstopidx, -1., pMGj->confidence, ConfRightLogPV, sv_ps, pMGj->xmapid, pMGi->xmapid, true /* trans_overlapped */);
	sv_array[numsv].inv_MGi = inv_MGi;
	sv_array[numsv].MGi_qryleft = MGi_qryleft;
	sv_array[numsv].repeat = (pMGi->have_ref_repeat || pMGj->have_ref_repeat); //move to method of sv class
	
	if(SV_ORIENT){
	  if(pMGj->orientforward){
	    if(pMGi->orientforward){
	      sv_array[numsv].orientation = strdup("+/+");
	    } else {
	      sv_array[numsv].orientation = strdup("-/+");	    
	    }
	  } else {
	    if(pMGi->orientforward){
	      sv_array[numsv].orientation = strdup("-/+");	    
	    } else {
	      sv_array[numsv].orientation = strdup("-/-");
	    }
	  }
	}
	numsv++;
	if(DEBUG) assert(numsv <= maxsv);
      }

      continue;
    }

    if( rearrangement || intra_rearr ){ // possible translocations, otherwise defaults to complex_sv
      //Heng's translocation criteria: (see parameters.{cpp,h})
      //double smap_translocation_maxsep = 500.; /* smap: for transloactions with no overlap, maximum separation on query (in kb) */
      //double smap_translocation_ovlfr  = 0.3; /* smap: for translocations with overlap, maximum overlap length/matchgroup length */
      //double smap_translocation_ovlmx  = 200.; /* smap: for translocations with overlap, maximum overlap length (in kb) */
      
      if( qryoverlap ) { 
	if( qryoverlapsize < smap_translocation_ovlmx &&
	    qryoverlapsize / minqryalignlen < smap_translocation_ovlfr ) { //these are Heng's two criteria for the overlap case
	  translocation = true;
	} else if( qryoverlapsize >= minqryalignlen){// handle case where smaller MG overlaps an internal outlier in the larger MG : most of these cases were handled previously
	  if(VERB>=1+RELEASE){
	    printf("WARNING: Unhandled overlap of smaller MG by larger MG : possible translocation\n");
	    fflush(stdout);// exit(1);
	  }
	  /* defaults to complex_sv below */
	}
      } // else (no overlap), -qryoverlapsize is really the gap size
      else if( -qryoverlapsize < smap_translocation_maxsep ) {
	translocation = true;
      }

      if(VERB && verbose){
	printf("\tqryoverlap=%d,qryoverlapsize= %0.3f,smap_translocation_ovlmx= %0.3f, smap_translocation_ovlfr= %0.6f, smap_translocation_maxsep= %0.3f: translocation=%d\n",
	       qryoverlap ? 1:0, qryoverlapsize, smap_translocation_ovlmx,smap_translocation_ovlfr,smap_translocation_maxsep,translocation ? 1 : 0);
	fflush(stdout);
      }

      if(!translocation) complex_sv = true;// defaults to complex_sv
    }
    else if (opporient){/* handle all opposite orientation cases : inversion, inverted duplication (otherwise default to complex) */

      double nonovsz = inv_refalignlen - max(0.0,refoverlapsize); //non-overlap size in reference inverted MG
      int nonovSites = inv_refalignSites - max(0,refoverlapSites);// non-overlap sites in reference inverted MG

      double qrynonovsz = inv_qryalignlen - max(0.0,qryoverlapsize);// non-overlap size of query inverted MG
      int qrynonovSites = inv_qryalignSites - max(0,qryoverlapSites);// non-overlap sites in query inverted MG
      
      if(VERB/* HERE HERE >=2 */ && verbose){
	printf("opporient=%d: xmapids %d %d(conf=%0.2f,%0.2f,fwd=%d,%d,inv_only=%d,%d),refid=%lld,qryid=%lld:qrygap=%d,%d(%0.3f,%0.3f)\n\t inv_MG=%d,NoInvDup=%d,NoInv=%d,nonovsz= %0.3f(%d),qrynonovsz=%0.3f(%d),qrysize= %0.3f,qryoverlap=%d(%0.3f kb, %d Sites),refoverlap=%d(%0.3f kb, %d Sites):inv_refalignlen=%0.3f,inv_refalignSites=%d(noninv=%d),qrystopidx=%d,qrystartidx=%d: checking for Inverted Duplication\n",
	       opporient ? 1 : 0, pMGi->xmapid, pMGj->xmapid, pMGi->confidence, pMGj->confidence, 
	       pMGi->orientforward, pMGj->orientforward,pMGi->inversion_only,pMGj->inversion_only,
	       pMGi->refcontigid, pMGi->qrycontigid, qrystartidx, qrystopidx, qrystart*1e-3,qrystop*1e-3,
	       usexmapentries[inv_MGi ? i : j]->xmapid, NoInvDup ? 1 : 0, NoInv ? 1 : 0, nonovsz,nonovSites,qrynonovsz,qrynonovSites,qrysize, 
	       qryoverlap ? 1 : 0, qryoverlapsize, qryoverlapSites, refoverlap ? 1 : 0, refoverlapsize, refoverlapSites, inv_refalignlen, inv_refalignSites, noninv_refalignSites, qrystopidx, qrystartidx);
	fflush(stdout);
      }
      if(DEBUG && qryoverlap) assert(qryoverlapSites > 0);
      if(DEBUG && !qryoverlap) assert(qrystopidx >= qrystartidx);

      /* NEW375 : check that inversion_only is consistent with inv_MGi */
      if(1 && (inv_MGi ? pMGj->inversion_only : pMGi->inversion_only)){/* Probably a missed Cut-Flip-Paste inversion */
	if(VERB/* HERE HERE >=2 */ && verbose){
	  printf("\t non-inverted MG has inversion_only=%d : ruling out Inversion or Inverted duplication (possibly missed small Inversion)\n", usexmapentries[inv_MGi ? j : i]->inversion_only);
	  fflush(stdout);
	}
	if(MATCHGROUP_PARTIAL_TRIM >= MATCHGROUP_PARTIAL_TRIM_COMPLEX /* NEW151 */ && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	  int origsmapsize = smapsize;
	  if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
	    trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	  } else {/* trim usexmapentries[j] where it overlaps on query */
	    trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	  }
	  if(DEBUG) assert(origsmapsize == smapsize);
	}
	
	continue;
      }

      if(InversionNormalLogPv > 0.0){       /* The non-inverted MG  needs to be above this LogPV threshold */ 
	if(usexmapentries[inv_MGi ? j : i]->confidence <= InversionNormalLogPv) {
	  if(VERB/* HERE HERE >=2 */ && verbose){
	    printf("\t non-inverted MG has conf=%0.2f,qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld : ruling out Inversion or Inverted duplication\n", 
		   usexmapentries[inv_MGi ? j : i]->confidence,qryoverlap,repeatcnt,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
	    fflush(stdout);
	  }

	  if(MATCHGROUP_PARTIAL_TRIM >= MATCHGROUP_PARTIAL_TRIM_COMPLEX /* NEW151 */ && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	    int origsmapsize = smapsize;
	    if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
	      trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	    } else {/* trim usexmapentries[j] where it overlaps on query */
	      trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	    }
	    if(DEBUG) assert(origsmapsize == smapsize);
	  }

	  continue;
	}
      }

      /* The inverted matchgroup needs to either be above -T threshold or satisfy the -CutFlip & -CutFlipBC rules */

      // Check if low confidence inverted matchgroup fails the -CutFlip or -CutFlipBC rules
      if(usexmapentries[inv_MGi ? i : j]->confidence <= (InversionNormalLogPv > 0.0 ? InversionNormalLogPv : LogPvThreshold)){
	if(inv_MGi){
	  int qrydist = abs(pMGi->qrystartidx - pMGi->qryendidx) + qrystopidx - qrystartidx;
	  int refdist = pMGj->refstartidx - pMGi->refstartidx;
	  if(0 && refdist > qrydist + CutFlipEndExpand){
	    if(VERB/* HERE HERE >=2 */ && verbose){
	      printf("\t inverted MG has conf=%0.2f: qrydist=%d, refdist=%d, CutFlipEndExpand=%d,qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld: Violates -CutFlip, skipping\n", 
		     pMGi->confidence, qrydist,refdist,CutFlipEndExpand,qryoverlap,repeatcnt,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
	      fflush(stdout);
	    }

	    if(MATCHGROUP_PARTIAL_TRIM>=2 && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	      int origsmapsize = smapsize;
	      if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
		trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	      } else {/* trim usexmapentries[j] where it overlaps on query */
		trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	      }
	      if(DEBUG) assert(origsmapsize == smapsize);
	    }

	    continue;
	  }

	  if(CutFlipBC){
	    double myLogPvThreshold4 = LogPvThreshold4 + log(max(1.0, qrydist+1.0-CutFlip) * max(1.0, refdist + 1.0 - CutFlip))/log(10.0);
	    if(pMGi->confidence <= myLogPvThreshold4){
	      if(VERB/* HERE HERE >=2 */ && verbose){
		printf("\t inverted MG has conf=%0.2f: BC threshold = %0.2f (inv_MGi=%d),qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld: skipping\n", 
		       pMGi->confidence, LogPvThreshold4, inv_MGi, qryoverlap, repeatcnt,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
		fflush(stdout);
	      }

	      if(MATCHGROUP_PARTIAL_TRIM >= MATCHGROUP_PARTIAL_TRIM_COMPLEX /* NEW151 */ && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
		if(DEBUG) assert(pMGi == usexmapentries[i]);
		if(DEBUG) assert(pMGj == usexmapentries[j]);

		int origsmapsize = smapsize;
		// HERE HERE : should defer (like complex pairs) to after SV loop
		if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
		  trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
		} else {/* trim usexmapentries[j] where it overlaps on query */
		  trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
		}
		if(DEBUG) assert(origsmapsize == smapsize);
	      }

	      continue;
	    }
	  }
	} else {// invMGi == false
	  int qrydist = abs(pMGj->qrystartidx - pMGj->qryendidx) + qrystopidx - qrystartidx;
	  int refdist = pMGj->refendidx - pMGi->refendidx;
	  if(0 && refdist > qrydist + CutFlipEndExpand){
	    if(VERB/* HERE HERE >=2 */ && verbose){
	      printf("\t inverted MG has conf=%0.2f: qrydist=%d, refdist=%d, CutFlipEndExpand=%d, qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld: Violates -CutFlip, skipping\n",
		     pMGj->confidence, qrydist,refdist,CutFlipEndExpand,qryoverlap,repeatcnt,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
	      fflush(stdout);
	    }

	    if(MATCHGROUP_PARTIAL_TRIM>=2 && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	      int origsmapsize = smapsize;
	      // HERE HERE : should defer (like complex pairs) to after SV loop
	      if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
		trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	      } else {/* trim usexmapentries[j] where it overlaps on query */
		trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	      }
	      if(DEBUG) assert(origsmapsize == smapsize);
	    }

	    continue;
	  }
	  if(CutFlipBC){
	    double myLogPvThreshold4 = LogPvThreshold4 + log(max(1.0, qrydist+1.0-CutFlip) * max(1.0, refdist + 1.0 - CutFlip))/log(10.0);
	    if(pMGj->confidence <= myLogPvThreshold4){
	      if(VERB/* HERE HERE >=2 */ && verbose){
		printf("\t inverted MG has conf=%0.2f: BC threshold = %0.2f (inv_MGi=%d),qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld: skipping\n", 
		       pMGj->confidence, LogPvThreshold4, inv_MGi, qryoverlap, repeatcnt, pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
		fflush(stdout);
	      }

	      if(MATCHGROUP_PARTIAL_TRIM >= MATCHGROUP_PARTIAL_TRIM_COMPLEX /* NEW151 */ && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
		if(DEBUG) assert(pMGi == usexmapentries[i]);
		if(DEBUG) assert(pMGj == usexmapentries[j]);

		int origsmapsize = smapsize;
		// HERE HERE : should defer (like complex pairs) to after SV loop
		if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
		  trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
		} else {/* trim usexmapentries[j] where it overlaps on query */
		  trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
		}
		if(DEBUG) assert(origsmapsize == smapsize);
	      }

	      continue;
	    }
	  }
	} // invMGi == false
      } // inverted MG confidence fails threshold for inversions
      
      // NOTE : Will later disable inversion or inverted duplication calls if inverted matchgroup has same orientation as dominant query orientation
      if(DEBUG && !qryoverlap) assert(qryoverlapSites <= 0);

      // NEW : Inverted duplication is now based on the number of labels in the ref overlap region vs inverted MG reference region AND query gap (should be small relative to ref overlap)

      // NEW11 : Inverted duplication is disabled if query overlap >= OverlapFilterQryRatio : If valid, this cases should have been called as Small Duplication
      double queryoverlap, queryminlen;
      if(1){
	double mini = min(pMGi->qrystartpos, pMGi->qryendpos);
	double minj = min(pMGj->qrystartpos, pMGj->qryendpos);
	double maxi = max(pMGi->qrystartpos, pMGi->qryendpos);
	double maxj = max(pMGj->qrystartpos, pMGj->qryendpos);
	queryoverlap = min(maxi,maxj) - max(mini,minj);
	queryminlen = min(maxi-mini, maxj-minj);      
      }

      // NOTE : Duplication Spacing <= smap_InvDup_maxSpacing is checked after coordinates are know 
      if( !NoInvDup && refoverlap && refoverlapSites >= inv_refalignSites * smap_InvDup_minRefOverlapFrac && qrystopidx - qrystartidx  <= refoverlapSites * smap_InvDup_maxQryGapMult 
	  && queryoverlap < max(smap_InvMaxOverlapSize, queryminlen * OverlapFilterQryRatio)){
	duplication = true;
	inversion = true; //for inverted duplication

	double DupSpacing = -1.0,qrystart,qrystop;// compute DupSpacing
	int refstartidx = max(pMGi->refstartidx, pMGj->refstartidx);// WAS11 pMGj->refstartidx;
	int refstopidx  = min(pMGi->refendidx, pMGj->refendidx);
	if(DEBUG && !(refstartidx <= refstopidx)){
	  printf("refstartidx=%d, refstopidx=%d\n",refstartidx,refstopidx);
	  fflush(stdout);
	  assert(refstartidx <= refstopidx);
	}

	if(qryoverlap){// NEW4 : trim away query overlap region from refstartidx..refstopidx
	  int refidx1,refidx2;
	  double refloc1,refloc2;
	  int mini = min(pMGi->qrystartidx, pMGi->qryendidx);
	  int maxi = max(pMGi->qrystartidx, pMGi->qryendidx);
	  int minj = min(pMGj->qrystartidx, pMGj->qryendidx);
	  int maxj = max(pMGj->qrystartidx, pMGj->qryendidx);
	  if(DEBUG) assert(min(maxi,maxj) >= max(mini,minj));
	  int qryidx1 = max(mini,minj), qryidx2 = min(maxi,maxj);// WAS11 qryidx1 = qrystartidx, qrydix2 = qrystopidx;
	  if(inv_MGi){	  /* map query overlap region to reference via inverted matchgroup i */
	    getRef(pMGi->align, qryidx1, refloc1, refidx1, -1);
	    getRef(pMGi->align, qryidx2, refloc2, refidx2, +1);
	  } else { /* map query overlap region to reference via inverted matchgroup j */
	    getRef(pMGj->align, qryidx1, refloc1,refidx1, -1);
	    getRef(pMGj->align, qryidx2, refloc2,refidx2, +1);
	  }
	  if(max(refidx1,refidx2) > refstartidx && min(refidx1,refidx2) <= refstartidx)
	    refstartidx = max(refidx1,refidx2);
	  if(min(refidx1,refidx2) < refstopidx && max(refidx1,refidx2) >= refstopidx)
	    refstopidx = min(refidx1,refidx2);
	}

	int qrystartidx = -1, qrystopidx = -1;
	if(refstartidx <= refstopidx)
	  getDuplicationCoords(pMGi, pMGj, refstartidx, refstopidx, qrystart, qrystartidx, qrystop, qrystopidx, DupSpacing, 0 /* verbose */);

	if(verbose){/* NOTE : these are not the final coordinates, which are computed later : these are just to determinate if this is a valid Inverted duplication */
	  printf("  inverted duplication ? ref interval: %d %d, qry interval: %d %d (%.3f %.3f), Dup Spacing= %0.3f kb\n", refstartidx, refstopidx, qrystartidx, qrystopidx, qrystart*1e-3, qrystop*1e-3, DupSpacing);
	  fflush(stdout);
	}

	if(!(refstartidx <= refstopidx) || DupSpacing > smap_InvDup_maxSpacing){
	  NoInvDup = true;
	  duplication = false;
	  inversion = false;
	  if(!NoInv && !NoInv && refoverlapSites < inv_refalignSites * (smap_Inv_maxRefOverlapFrac > 0.0 ? smap_Inv_maxRefOverlapFrac : smap_InvDup_minRefOverlapFrac)// NEW151
	     && qryoverlapSites < inv_qryalignSites * (smap_Inv_maxQryOverlapFrac > 0.0 ? smap_Inv_maxQryOverlapFrac : smap_InvDup_minRefOverlapFrac) 
	     && nonovSites >= smap_Inv_minRefSites && qrynonovSites >= smap_Inv_minQrySites 
	     && qryoverlapsize < /* WAS155 qryoverlap > */ (smap_InvMaxOverlapSize ? smap_InvMaxOverlapSize : queryminlen * OverlapFilterQryRatio)/* NEW151 */)/* same test as below */	  
	    inversion = true;
	}
      } else if(!NoInv && refoverlapSites < inv_refalignSites * (smap_Inv_maxRefOverlapFrac > 0.0 ? smap_Inv_maxRefOverlapFrac : smap_InvDup_minRefOverlapFrac)// NEW151
		&& qryoverlapSites < inv_qryalignSites * (smap_Inv_maxQryOverlapFrac > 0.0 ? smap_Inv_maxQryOverlapFrac : smap_InvDup_minRefOverlapFrac)
		&& nonovSites >= smap_Inv_minRefSites && qrynonovSites >= smap_Inv_minQrySites
		&& qryoverlapsize < /* WAS155 qryoverlap > */ (smap_InvMaxOverlapSize ? smap_InvMaxOverlapSize : queryminlen * OverlapFilterQryRatio)/* NEW151 */){
	inversion = true;
      } else { //this is complex
	if(verbose)
	  printf("opporient else\n"); //can this happen? Yes, rarely, but it can.
      }
      if(verbose)
	printf("opporient %5lld: xmapids %3i %3i; gapsize=%.3f %.3f; align: %.3f %.3f; overlapsize: qry= %.3f ref= %.3f; nonovsz: %.3f; inv %i, dup %i\n", 
	       contigid, pMGi->xmapid, pMGj->xmapid, qrysize, refsize, qrysize+minqryalignlen, refsize+minrefalignlen, qryoverlapsize, refoverlapsize, nonovsz, inversion, duplication);

      if(!inversion){
	if(VERB>=2){
	  printf("complex: qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld: skipping\n",qryoverlap,repeatcnt,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
	  fflush(stdout);
	}
	if(MATCHGROUP_PARTIAL_TRIM>=2/* NEW151 */ && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	  int origsmapsize = smapsize;
	  // HERE HERE : should defer (like complex pairs) to after SV loop
	  if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
	    trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	  } else {/* trim usexmapentries[j] where it overlaps on query */
	    trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	  }
	  if(DEBUG) assert(origsmapsize == smapsize);
	}

	continue;//	complex_sv = true;// NOTE : complex SV should not be output
      }

    } else { // all remaining cases have !rearrangement && !intra_rearr && !opporient (Thomas)

      int qrygapSites = abs(pMGi->qrystartidx - pMGj->qryendidx); 

      if(VERB/* HERE HERE >=2 */ && verbose){
	printf("\t opporient=%d: xmapids %d %d;refid=%lld,qryid=%lld:qrysize= %0.3f:qryoverlapSites=%d,refoverlapSites=%d,minrefalignSites=%d,qrystopidx=%d,qrystartidx=%d,Split_qrygapSites=%d,maxq1=%0.3f,minq1=%0.3f,Extend=%d,%d: checking for non-Inverted Duplication\n",
	       opporient ? 1 : 0, pMGi->xmapid, pMGj->xmapid, 
	       pMGi->refcontigid, pMGi->qrycontigid, qrysize, qryoverlapSites, refoverlapSites, minrefalignSites, qrystopidx, qrystartidx,qrygapSites,maxq1*0.001,minq1*0.001,LeftExtend?1:0,RightExtend?1:0);
	printf("\t MG[i]->refendidx=%d, MG[j]->refstartidx=%d\n",pMGi->refendidx,pMGj->refstartidx);
	fflush(stdout);
      }

      if(DEBUG) assert(!rearrangement && !intra_rearr && !opporient);
      if(DEBUG) assert(pMGi->refcontigid == pMGj->refcontigid);

      if(qrysize < -(maxq1-minq1) * 0.001 ) { // Left ends of matchgroups are crossed : check for duplication_split OR duplication (otherwise defaults to complex)
	if(DEBUG) assert(!opporient);
	if( qryoverlapSites <= smap_DupSplit_maxoverlap && qrygapSites <= smap_DupSplit_maxgap && (LeftExtend || RightExtend) && 
	    ((LeftExtend && RightExtend) || smap_DupSplit_singleEnd >= 3 || 
	     (smap_DupSplit_singleEnd >= 1 && refoverlapSites >= 1) ||
	     (smap_DupSplit_singleEnd >= 2 && ((!LeftExtend ? (minq1 < maxq2) : (maxq1 > minq2)) ? (maxq1-minq1) : (maxq2-minq2)) > 0.5 * (pMGj->refendpos - pMGi->refstartpos)))){
	  duplication = duplication_split = true;
	  if(LeftSites <= smap_DupSplit_maxend)
	    LeftExtend = true;// NEW25
	  if(RightSites <= smap_DupSplit_maxend)
	    RightExtend = true;// NEW25
	  if(verbose){
	    printf("\t dupli-split %lld,qryid=%lld: xmapids= %3i %3i: qrysize=%.3f, qryoverlapSites=%d, qrygapSites=%d, LeftSites=%d,RightSites=%d,LeftExtend=%d,RightExtend=%d\n", 
		   contigid, pMGi->qrycontigid, pMGi->xmapid, pMGj->xmapid, qrysize, qryoverlapSites,qrygapSites,LeftSites,RightSites,LeftExtend?1:0,RightExtend?1:0);
	    fflush(stdout);
	  }
	} else if(refoverlap && refoverlapSites >= smap_Dup_minRefOverlapSites + max(0,qryoverlapSites) && refoverlapsize > qryoverlapsize &&
		  (pMGi->orientforward ? (pMGi->qrystartidx - pMGj->qryendidx) : (pMGj->qryendidx - pMGi->qrystartidx)) <=
		  ((pMGi->orientforward ? LeftExtend : RightExtend) ? pMGj->refendidx - pMGi->refstartidx : refoverlapSites) * smap_Dup_maxQryGapMult){// NEW278 : Left ends of matchgroups are Crossed and MGs overlap on reference AND gap is not too large : duplication 
	  duplication = true;
	  if(verbose){
	    printf("\t ref overlp3: qryid=%lld: xmapids= %d,%d:qryoverlapSites=%d,refoverlapSites=%d,qrystartidx=%d,qrystopidx=%d,fwd=%d,LeftExtend=%d,RightExtend=%d\n",
		   pMGi->qrycontigid, pMGi->xmapid, pMGj->xmapid, qryoverlapSites,refoverlapSites, pMGi->qrystartidx, pMGj->qryendidx, 
		   pMGi->orientforward ? 1 : 0, LeftExtend ? 1 : 0, RightExtend ? 1 : 0);
	    fflush(stdout);
	  }
	} else {
	  // HERE HERE HERE : should call transposition (or intra-tranlocation pair), IF the 2 flanking MGs (with query orientation like current 2 MGs) don't cross either of these 2 MGs
	  complex_sv = true;// crossed with large query overlap or MG not near end of contigs : defaults to complex-sv

	  if(VERB>=2){
	    printf("complex: qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld: skipping\n",qryoverlap,repeatcnt,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
	    fflush(stdout);
	  }

	  if(MATCHGROUP_PARTIAL_TRIM >= MATCHGROUP_PARTIAL_TRIM/* NEW151 */ && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	    int origsmapsize = smapsize;
	    // HERE HERE : should defer (like complex pairs) to after SV loop
	    if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
	      trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	    } else {/* trim usexmapentries[j] where it overlaps on query */
	      trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	    }
	    if(DEBUG) assert(origsmapsize == smapsize);
	  }
	}
      }

      else if( refoverlap && refoverlapSites >= smap_Dup_minRefOverlapSites + max(0, qryoverlapSites) && /* NEW2 */ refoverlapsize > qryoverlapsize && 
	       qrystopidx - qrystartidx <= ((pMGi->orientforward ? RightExtend : LeftExtend) ? pMGi->refendidx - pMGj->refstartidx : refoverlapSites) * smap_Dup_maxQryGapMult){// Left ends of matchgroups are NOT crossed and MGs overlap on reference AND gap is not too large : duplication
	if(DEBUG) assert(inversion == false);

	duplication = true;
	if(verbose){
	  printf("\t ref overlp2: qryid=%lld: xmapids= %d,%d:refoverlapSites=%d,qrystartidx=%d,qrystopidx=%d,fwd=%d,LeftExtend=%d,RightExtend=%d\n",
		 pMGi->qrycontigid, pMGi->xmapid, pMGj->xmapid, refoverlapSites, qrystartidx, qrystopidx, 
		 pMGi->orientforward ? 1 : 0, LeftExtend ? 1 : 0, RightExtend ? 1 : 0);
	  fflush(stdout);
	}
      }

      // remaining cases below handle 2-matchgroup indels using Warren's original code

      else if((qryoverlap && minqryalignlen - qryoverlapsize < smap_side_alignkb) ||
	      (refoverlap && minrefalignlen - refoverlapsize < smap_side_alignkb)) {
	if(verbose){
	  printf("\t flank too small: qryid= %lld, xmapids %3i %3i; %.3f %.3f; align: %.3f %.3f; overlapsize: %.3f; new ovlpsz q: %.3f r: %.3f\n", contigid, pMGi->xmapid, pMGj->xmapid, qrysize, refsize, qrysize+minqryalignlen, refsize+minrefalignlen, qryoverlapsize, qryoverlapsize, refoverlapsize);
	  fflush(stdout);
	}

	// NOTE : any case with qryoverlap or refoverlap that is NOT inversion or duplication or overlap_indel will be output as a complex (or not at all)
	if(VERB>=2){
	  printf("complex: qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld,: skipping\n",qryoverlap,repeatcnt,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
	  fflush(stdout);
	}
	if(MATCHGROUP_PARTIAL_TRIM >= MATCHGROUP_PARTIAL_TRIM_COMPLEX/* NEW151 */ && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	  int origsmapsize = smapsize;
	  // HERE HERE : should defer (like complex pairs) to after SV loop
	  if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
	    trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	  } else {/* trim usexmapentries[j] where it overlaps on query */
	    trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	  }
	  if(DEBUG) assert(origsmapsize == smapsize);
	}

	continue;// WAS33 complex_sv = true;// default to complex SV
      }

      else if( qryoverlap /* WAS29 && qryoverlapsize < sv_indel_maxoverlap && fabs(qrysize + qryoverlapsize) <= 0.001 */ && (!refoverlap /* NEW29 */|| qryoverlapsize > refoverlapsize)) {/* overlapped deletion */
	if(qryoverlapsize > smap_del_maxqryoverlap){// -svDelMaxQryOverlap 140 (WAS405 110) : assume repeat compression (assembly error)
	  if(verbose){
	    printf("\t qry overlap too large:qryid=%lld,refid=%lld,%lld: qry gap= %0.3f .. %0.3f (%i .. %i), qrysize= %.3f, ref= %0.3f .. %0.3f (%d .. %d), refsize=%0.3f, qryoverlapsize= %0.3f, refoverlapsize=%0.3f\n", 
		   contigid, pMGi->refcontigid,pMGj->refcontigid,qrystart * 1e-3, qrystop * 1e-3, qrystartidx, qrystopidx, qrysize, pMGi->refendpos * 1e-3, pMGj->refstartpos * 1e-3, 
		   pMGi->refendidx, pMGj->refstartidx, refsize, qryoverlapsize, refoverlapsize);
	    fflush(stdout);
	  }

	  // NOTE : any case with qryoverlap or refoverlap that is NOT inversion or duplication or overlap_indel will be output as a complex SV (or not at all)
	  // NOTE : Avoid trimming query overlap so user can see why deletion was NOT called (due to large query overlap)
	  if(0&&/*NEW48*/ MATCHGROUP_PARTIAL_TRIM && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	    int origsmapsize = smapsize;
	    // HERE HERE : should defer (like complex pairs) to after SV loop
	    if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
	      trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	    } else {/* trim usexmapentries[j] where it overlaps on query */
	      trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	    }
	    if(DEBUG) assert(origsmapsize == smapsize);
	  }

	  continue;// WAS33 complex_sv = true;// default to complex SV
	} else { // qryoverlapsize <= smap_del_maxqryoverlap

	  overlap_indel = true;
	  //start/stops are in bp, sizes are in kb: qrysize is the overlap
	  if(VERB/* HERE HERE >=2 */ && verbose){
	    printf("\t qry overlap:qryid=%lld, qry gap= %0.3f .. %0.3f (%i .. %i), qrysize= %.3f, ref= %0.3f .. %0.3f (%d .. %d), refsize=%0.3f, qryoverlapsize= %0.3f, refoverlapsize=%0.3f\n", 
		   contigid, qrystart * 1e-3, qrystop * 1e-3, qrystartidx, qrystopidx, qrysize, pMGi->refendpos * 1e-3, pMGj->refstartpos * 1e-3, 
		   pMGi->refendidx, pMGj->refstartidx, refsize, qryoverlapsize, refoverlapsize);
	    fflush(stdout);
	  }

	  if(DEBUG) assert( fabs(qrysize + qryoverlapsize) <= 0.001);

	  // NEW29 : for query remove the overlap region from both alignments to reflect the uncertainty about which part of reference was deleted in query 
	  if(VERB>=2){
	    printf("\t  Before getRef: ref = %d..%d(%0.3f..%0.3f), qry = %d..%d(%0.3f..%0.3f)\n",
		   refstartidx, refstopidx, refstart * 1e-3, refstop * 1e-3,qrystartidx,qrystopidx,qrystart*1e-3,qrystop*1e-3);
	    fflush(stdout);
	  }

	  double tmp = qrystart; qrystart = qrystop; qrystop = tmp;
	  int itmp = qrystartidx; qrystartidx = qrystopidx; qrystopidx = itmp;
	  if(pMGi->orientforward){
	    getRef(pMGi->align, qrystartidx, refstart, refstartidx, -1, &qrystart);
	    getRef(pMGj->align, qrystopidx, refstop, refstopidx, 1,  &qrystop);
	  } else {
	    getRef(pMGi->align, qrystopidx, refstart, refstartidx, 1, &qrystop );
	    getRef(pMGj->align, qrystartidx, refstop, refstopidx, -1, &qrystart );
	  }

	  refsize = (refstop - refstart) * 1e-3; //size is in kb, positions are bp
	  qrysize = (qrystop - qrystart) * 1e-3;

	  if(VERB>=2){
	    printf("\t  After getRef: ref = %d..%d(%0.3f..%0.3f), qry = %d..%d(%0.3f..%0.3f)\n",
		   refstartidx, refstopidx, refstart * 1e-3, refstop * 1e-3,qrystartidx,qrystopidx,qrystart*1e-3,qrystop*1e-3);
	    fflush(stdout);
	  }


	  sv_type_t orig_svtype_enum = svtype_enum;
	  if(qrysize < refsize && qryoverlapsize > smap_del_bigqryoverlap) // -svDelBigQryOverlap 110 (WAS404 60) kb
	    svtype_enum = compression;/* low confidence deletion */
	  else
	    svtype_enum = (qrysize > refsize ? insertion : deletion);// may no longer be deletion

	  if(verbose){
	    printf("\t qry overlap indel: type_enum= %d -> %d : qryid= %lld, qry= %0.3f ..%0.3f(%i ..%i), qrysize= %.3f, ref= %0.3f .. %.3f(%i..%i),refsize=%.3f, qryoverlapsize= %0.3f\n", 
		   orig_svtype_enum, svtype_enum,
		   contigid, qrystart*1e-3, qrystop*1e-3, qrystartidx, qrystopidx, qrysize, refstart*1e-3, refstop*1e-3, refstartidx, refstopidx, refsize, qryoverlapsize);
	    printf("\t minqryalignlen= %0.4f kb, minrefalignlen= %0.4f kb\n",minqryalignlen, minrefalignlen);
	    fflush(stdout);
	  }
	}
      }

      else if( refoverlap /* WAS29 && refoverlapsize < sv_indel_maxoverlap && qrysize > 0 */ /* NEW29*/&& (!qryoverlap || refoverlapsize > qryoverlapsize)) {/* overlapped insertion */
	overlap_indel = true;
	if(verbose){
	  printf("\t ref overlap:qryid=%lld, qry gap= %0.3f .. %0.3f (%i .. %i), qrysize= %.3f, ref= %0.3f .. %0.3f (%d .. %d), refsize=%0.3f, qryoverlapsize= %0.3f, refoverlapsize=%0.3f\n", 
		 contigid, qrystart * 1e-3, qrystop * 1e-3, qrystartidx, qrystopidx, qrysize, pMGi->refendpos * 1e-3, pMGj->refstartpos * 1e-3, 
		 pMGi->refendidx, pMGj->refstartidx, refsize, qryoverlapsize, refoverlapsize);
	  fflush(stdout);
	}
	if(DEBUG) assert( fabs(refsize + refoverlapsize) <= 0.001);

	// NEW29 : for ref remove the overlap region from both alignments to reflect the uncertainty about which part of qry was inserted into reference
	refstart = pMGj->refstartpos;
	refstop = pMGi->refendpos;
	refstartidx = pMGj->refstartidx;
	refstopidx = pMGi->refendidx;
	if(DEBUG) assert(refstart < refstop);
	if(DEBUG) assert(refstartidx < refstopidx);

	if(VERB>=2){
	  printf("\t  Before getQry: ref -> %d..%d(%0.3f..%0.3f)\n",refstartidx, refstopidx, refstart * 1e-3, refstop * 1e-3);
	  fflush(stdout);
	}

	// NOTE : qrystart <= qrystop is always true
	if(pMGi->orientforward){
	  getQry(pMGi->align, refstartidx, qrystart, qrystartidx, -1, &refstart);
	  getQry(pMGj->align, refstopidx, qrystop, qrystopidx, 1, &refstop);
	} else {
	  getQry(pMGi->align, refstartidx, qrystop, qrystopidx, -1, &refstart);
	  getQry(pMGj->align, refstopidx, qrystart, qrystartidx, 1, &refstop);
	}

	if(VERB>=2){
	  printf("\t  After getQry: ref -> %d..%d(%0.3f..%0.3f), qry -> %d..%d(%0.3f..%0.3f)\n",
		 refstartidx, refstopidx, refstart * 1e-3, refstop * 1e-3,qrystartidx,qrystopidx,qrystart*1e-3,qrystop*1e-3);
	  fflush(stdout);
	}

	if(qrystart > qrystop){/* NEW39 : the MGs must have embedded internal outliers : remove qry overlap region as well from both alignments */
	  int itmp = qrystartidx;
	  double tmp = qrystart;
	  qrystartidx = qrystopidx;
	  qrystart = qrystop;
	  qrystopidx = itmp;
	  qrystop = tmp;

	  if(pMGi->orientforward){
	    getRef(pMGi->align, qrystartidx, refstart, refstartidx, -1, &qrystart);
	    getRef(pMGj->align, qrystopidx, refstop, refstopidx, 1,  &qrystop);
	  } else {
	    getRef(pMGi->align, qrystopidx, refstart, refstartidx, 1, &qrystop );
	    getRef(pMGj->align, qrystartidx, refstop, refstopidx, -1, &qrystart );
	  }
	  	
	  if(VERB>=2){
	    printf("\t  After getRef: ref -> %d..%d(%0.3f..%0.3f), qry -> %d..%d(%0.3f..%0.3f)\n",
		   refstartidx, refstopidx, refstart * 1e-3, refstop * 1e-3,qrystartidx,qrystopidx,qrystart*1e-3,qrystop*1e-3);
	    fflush(stdout);
	  }

	  if(DEBUG && !(refstart <= refstop)){
	    printf("\t ref overlap insertion: contigid= %lld,qry= %.3f .. %.3f(%i .. %i) qrysize=%.3f, ref= %.3f .. %.3f(%i .. %i) refsize=%.3f\n", 
		   contigid, qrystart * 1e-3, qrystop * 1e-3, qrystartidx, qrystopidx, qrysize, refstart * 1e-3, refstop * 1e-3, refstartidx, refstopidx, refsize);
	    fflush(stdout);
	    assert(refstart <= refstop);
	  }
	}

	if(DEBUG && !(qrystart <= qrystop)){
	  printf("\t ref overlap insertion: contigid= %lld,qry= %.3f .. %.3f(%i .. %i) qrysize=%.3f, ref= %.3f .. %.3f(%i .. %i) refsize=%.3f\n", 
		 contigid, qrystart * 1e-3, qrystop * 1e-3, qrystartidx, qrystopidx, qrysize, refstart * 1e-3, refstop * 1e-3, refstartidx, refstopidx, refsize);
	  printf("qrystart= %0.1f, qrystop= %0.1f bp\n",qrystart,qrystop);
	  fflush(stdout);
	  assert(qrystart <= qrystop);
	}

	refsize = (refstop - refstart) * 1e-3; //size is in kb, positions are bp
	qrysize = (qrystop - qrystart) * 1e-3;

	if(verbose)
	  printf("\t ref overlap indel: type_enum= %d -> %d, contigid= %lld,qry= %.3f .. %.3f(%i .. %i) qrysize=%.3f, ref= %.3f .. %.3f(%i .. %i) refsize=%.3f\n", 
		 svtype_enum, (qrysize > refsize ? insertion : deletion), contigid, 
		 qrystart * 1e-3, qrystop * 1e-3, qrystartidx, qrystopidx, qrysize, refstart * 1e-3, refstop * 1e-3, refstartidx, refstopidx, refsize);

	svtype_enum = (qrysize > refsize ? insertion : deletion);// may no longer be insertion
      }
      else {/* regular 2-MG Indel without any overlaps */
	if(verbose){
	  printf("\t no overlap Indel: qryid= %lld, qry= %0.3f .. %0.3f (%i .. %i), qrysize= %.3f, ref= %0.3f .. %0.3f (%d .. %d), refsize=%0.3f, xmapids= %i %i : inv=%d,dupl=%d\n", 
		 contigid, qrystart * 1e-3, qrystop * 1e-3, qrystartidx, qrystopidx, qrysize, pMGi->refendpos * 1e-3, pMGj->refstartpos * 1e-3, 
		 pMGi->refendidx, pMGj->refstartidx, refsize, 
		 pMGi->xmapid, pMGj->xmapid, inversion,duplication);
	  fflush(stdout);
	}
	if(DEBUG) assert(qrysize >= 0.0);
	if(DEBUG) assert(refsize >= 0.0);
      }
    }// end of !rearrangement && !instr_rearr && !opporient


    if(0&& QRY_PRIMARY_ORIENTATION && inversion){/* filter out inversions or inverted duplications IFF dominant qry orientation does not match the non-inverted matchgroup's qry orientation
						    AND the inversion matchgroup is larger than the non-inverted matchgroup
						    AND the inversion matchgroup is NOT ( flanked by a non-inverted matchgroup on the other side  AND both non-inverted matchgroups are some minimum fraction of the inverted matchgroup AND the gaps are not too large relative to the non-inverted matchgroups )
 */
      // HERE HERE

      if(DEBUG && !(pMGi->dominant_orientation == pMGj->dominant_orientation)){
        printf("\t inversion=%d,duplication=%d,complex_sv=%d: dominant_orientation=%d,%d, fwd=%d,%d\n",
	       inversion,duplication,complex_sv, pMGi->dominant_orientation, pMGj->dominant_orientation, pMGi->orientforward, pMGj->orientforward);
        fflush(stdout);
        assert(pMGi->dominant_orientation == pMGj->dominant_orientation);/* since qrycontigid and refcontigid are same */
      }
      if(DEBUG) assert(!complex_sv);

      bool wrong_orientation = false;

      if(pMGi->dominant_orientation == 1){/* dominant qry orientation is forward */
	if(inv_MGi){/* MG i is the inverted matchgroup for calling inversion or inverted duplication */
	  if(!pMGj->orientforward)
	    wrong_orientation = true;
	} else {/* MG j is the inverted matchgroup for calling inversion or inverted duplication */
	  if(!pMGi->orientforward)
	    wrong_orientation = true;
	}
      } else if(pMGi->dominant_orientation == 2){/* dominant qry orientation is backwards */
	if(inv_MGi){/* MG i is the inverted matchgroup for calling inversion or inverted duplication */
	  if(pMGj->orientforward)
	    wrong_orientation = true;
	} else {/* MG j is the inverted matchgroup for calling inversion or inverted duplication */
	  if(pMGi->orientforward)
	    wrong_orientation = true;
	}
      }

      if(wrong_orientation){
	if(verbose){
	  printf("\t suppressed inversion or inverted-duplication due to invalid dominant qry orientation = %d : inv_MGi=%d, fwd= %d,%d, qryoverlap=%d,repeatcnt=%d,qryid=%lld,refid=%lld,%lld\n",
		 pMGi->dominant_orientation, inv_MGi, pMGi->orientforward, pMGj->orientforward,qryoverlap,repeatcnt,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid);
	  fflush(stdout);
	}

	if(MATCHGROUP_PARTIAL_TRIM >= MATCHGROUP_PARTIAL_TRIM_COMPLEX/* NEW151 */ && !repeatcnt && qryoverlap && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
	  int origsmapsize = smapsize;
	  // HERE HERE : should defer (like complex pairs) to after SV loop
	  if(pMGi->confidence < pMGj->confidence){/* trim pMGi where it overlaps on query */
	    trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false);
	  } else {/* trim usexmapentries[j] where it overlaps on query */
	    trimcnt += checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false);
	  }
	  if(DEBUG) assert(origsmapsize == smapsize);
	}

	continue;
      }
    }

    if(DEBUG && ((rearrangement || intra_rearr) ? !translocation :
		 opporient ? !inversion :
		 (duplication || duplication_split) ? 0 :
		 (qryoverlap || qrysize < 0 || refsize < 0) ? !overlap_indel :
		 0)){
      if(complex_sv != true){
	printf("\t rearrangement=%d, intra_rearr=%d, translocation=%d, opporient=%d, inversion=%d, duplication=%d, duplication_split=%d, qryoverlap=%d,qrysize=%0.3f,refsize=%0.3f,overlap_indel=%d,complex_sv=%d\n",
	       rearrangement,intra_rearr,translocation,opporient,inversion,duplication,duplication_split,qryoverlap,qrysize,refsize,overlap_indel,complex_sv);
	fflush(stdout);
	assert(complex_sv == true);
      }
    }
    // NOTE : complex_sv can also be true for crossed matchgroups that are not duplication_split

    if(MATCHGROUP_PARTIAL_TRIM && !repeatcnt && qryoverlap && !complex_sv && (MATCHGROUP_PARTIAL_TRIM >= 2 || !rearrangement)){
      // NOTE : This may change ordering of matchgroups : need to flip MG pairs when needed (including for current MG pair)
      // NOTE : should prioritize trimming from highest confidence calls to lowest confidence calls : for now just skip complex calls and do them later (if needed)
      bool trim_MGi = false, trim_big = false;
      if(!inversion /* duplication_split || (duplication && !inversion) || overlap_indel */){/* trim smaller matchgroup to remove query overlap region */
	trim_MGi = (pMGi->confidence < pMGj->confidence); // WAS41 */ (maxq1 - minq1 < maxq2 - minq2);
      } else /* if(inversion) */{/* trim inverted matchgroup (see inv_MGi) to remove query overlap region */
	if(DEBUG) assert(opporient);
	trim_MGi = inv_MGi;

	// NOTE : in this case the larger MG may be trimmed : If so need to remove any indels from trimmed region. Also try to avoid this (see where inv_MGi is determined) if MGi is much larger than MGj
	trim_big = (trim_MGi != (pMGi->confidence < pMGj->confidence /* WAS78 maxq1 - minq1 < maxq2 - minq2*/));
      }

      if(VERB>=2){
	printf("qryoverlap=%d,repeatcnt=%d,complex_sv=%d,qryid=%lld,refid=%lld,%lld: trimming partial query overlap (trim_MGi=%d,trim_big=%d)\n",
	       qryoverlap,repeatcnt,complex_sv,pMGi->qrycontigid,pMGi->refcontigid,pMGj->refcontigid,trim_MGi,trim_big);
	fflush(stdout);
      }

      int origsmapsize = smapsize;
      int trimmed = 0;
      if(trim_MGi){/* trim pMGi where it overlaps on query */
	trimcnt += trimmed = checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,j, i, false, false, trim_big ? numsvIndel1 : -1, true);
      } else {/* trim usexmapentries[j] where it overlaps on query */
	trimcnt += trimmed = checkOverlap(usexmapentries, smapsize,xmapentries, numxmapentries, maxxmapentries,i, j, false, false, trim_big ? numsvIndel1 : -1, false);
      }
      if(DEBUG) assert(origsmapsize == smapsize);

      if(trimmed /* WAS42 && overlap_indel*/){/* Coordinates for overlap_indel have already been computed and may no longer be valid */
	if(++repeatcnt2 < 100)
	  goto Lrepeat;
	  
	if(VERB/* HERE HERE >=2 */){
	  printf("WARNING: indp=%d/%d: repeatcnt= %d,%d: skipping further repeats for this MG pair\n",indp,indiceslen,repeatcnt,repeatcnt2);
	  fflush(stdout);
	}
      }

      if(0 && pMGi->refcontigid == pMGj->refcontigid && !(pMGi->refstartpos <= pMGj->refstartpos && (pMGi->refstartpos < pMGj->refstartpos || pMGi->refendpos <= pMGj->refendpos))){
	if(DEBUG) assert(trimmed);
	goto Lrepeat2;
      }

      if(0 && trimmed){ // update local variables depending coordinates of usexmapentries[i,j] : minq1, maxq1, minq2, maxq2, qrystartidx, qrystopidx etc

	refoverlap = overlap(pMGi->refstartpos, pMGi->refendpos, pMGj->refstartpos, pMGj->refendpos);
	refoverlapsize = overlapLength(pMGi->refstartpos, pMGi->refendpos, pMGj->refstartpos, pMGj->refendpos) * 1e-3; //kb!
	refoverlapSites = min(pMGi->refendidx, pMGj->refendidx) - max(pMGi->refstartidx, pMGj->refstartidx);

	minq1 = min(pMGi->qrystartpos, pMGi->qryendpos);
	maxq1 = max(pMGi->qrystartpos, pMGi->qryendpos);
	minq2 = min(pMGj->qrystartpos, pMGj->qryendpos);
	maxq2 = max(pMGj->qrystartpos, pMGj->qryendpos);

	qryoverlap = overlap( minq1, maxq1, minq2, maxq2 );
	qryoverlapsize = overlapLength( minq1, maxq1, minq2, maxq2 ) * 1e-3; //kb!

	iminq1 = min(pMGi->qrystartidx, pMGi->qryendidx);
	imaxq1 = max(pMGi->qrystartidx, pMGi->qryendidx);
	iminq2 = min(pMGj->qrystartidx, pMGj->qryendidx);
	imaxq2 = max(pMGj->qrystartidx, pMGj->qryendidx);
	qryoverlapSites = min(imaxq1, imaxq2) - max(iminq1, iminq2);

	if( pMGi->orientforward && pMGj->orientforward ) { //both forward
	  qrystart    = pMGi->qryendpos;
	  qrystartidx = pMGi->qryendidx;
	  qrystop     = pMGj->qrystartpos;
	  qrystopidx  = pMGj->qrystartidx;
	} else if( !pMGi->orientforward && !pMGj->orientforward ) { //both backward
	  qrystart    = pMGj->qrystartpos;
	  qrystartidx = pMGj->qrystartidx;
	  qrystop     = pMGi->qryendpos;
	  qrystopidx  = pMGi->qryendidx;
	} else { //opposite orientation
	  opporient = true;
	  if(minq1 < minq2) { //previous query (MG i) starts before current query (MG j) starts
	    qrystart = maxq1;
	    qrystartidx = imaxq1;
	    qrystop  = minq2;
	    qrystopidx  = iminq2;
	  } else if( minq2 < minq1) { //current query (MG j) starts before previous query (MG i) starts
	    qrystart = maxq2;
	    qrystartidx = imaxq2;
	    qrystop  = minq1;
	    qrystopidx  = iminq1;
	  } else { // minq1 == minq2 (NOTE : this case is not possible if OverlapFilterQryRatio < 1.0, since then the smaller matchgroups in the query would have been filtered out)
	    if(DEBUG) assert(OverlapFilterQryRatio >= 0.9999 || OverlapFilterConfDelta > 0.0);
#if 1 // NEW CODE	    
	    if(maxq1 > maxq2) {/* treat qry MG j as left of qry MG i */
	      qrystart = maxq2;
	      qrystartidx = imaxq2;
	      qrystop  = minq1;
	      qrystopidx  = iminq1;
	    } else { /* maxq1 <= maxq2 : treat qry MG i as left of qry MG j */
	      qrystart = maxq1;
	      qrystartidx = imaxq1;
	      qrystop  = minq2;
	      qrystopidx  = iminq2;
	    }
#else // OLD CODE
	    qrystart = minq1;
	    qrystartidx = iminq1; // WAS30 (pMGi->qrystartpos < pMGi->qryendpos) ? pMGi->qrystartidx : pMGi->qryendidx;
	    if(maxq1 > maxq2) {/* treat qry MG j as left of qry MG i */
	      qrystop = maxq2;
	      qrystopidx = imaxq2; // WAS30 (pMGj->qrystartpos > pMGj->qryendpos) ? pMGj->qrystartidx : pMGj->qryendidx;
	    } else { /* maxq1 <= maxq2 : treat qry MG i as left of qry MG j */
	      qrystop = maxq1;
	      qrystopidx = imaxq1; // WAS30 (pMGi->qrystartpos > pMGi->qryendpos) ? pMGi->qrystartidx : pMGi->qryendidx;
	    }
#endif // OLD CODE
	  }
	} /* opposite orientation */

	qrysize = (qrystop - qrystart) * 1e-3;
	refsize = (refstop - refstart) * 1e-3; //size is in kb, positions are bp
      }// if(trimmed)
    }

    // Now figure out the coordinates to output in SMAP for each SV type
    if( duplication_split  || (duplication && !inversion && 
			       !(pMGj->refendpos > pMGi->refendpos && pMGj->refstartpos > pMGi->refstartpos))) { //figure out the coords using duplication_split logic

      /* determine which matchgroup is left on the query (oriented with MG not twisted relative to ref). Since overlap on query must be minimal, this should always be unambigous (Thomas) */
      /* Unless mentioned, all positions are in bp */
      bool QueryIleftofJ;
      if(pMGi->orientforward ? (pMGj->qryendpos > pMGi->qryendpos) : 
	 (pMGj->qryendpos < pMGi->qryendpos)){/* query i is left of query j : typically this means there is more overlap on reference */
	QueryIleftofJ = true;
	refstart = pMGj->refstartpos;
	refstartidx = pMGj->refstartidx;
	if(pMGj->orientforward ? RightExtend : LeftExtend){
	  refstop = pMGi->refendpos;
	  refstopidx = pMGi->refendidx;
	} else {
	  refstop = min(pMGi->refendpos, pMGj->refendpos);
	  refstopidx = min(pMGi->refendidx, pMGj->refendidx);
	}
      } else {/* query j is left of query i : typically this means there is little or no overlap on reference (includes original duplication_split with no overlap on reference) */
	QueryIleftofJ = false;
	if(pMGj->orientforward ? LeftExtend : RightExtend){
	  refstart    = pMGi->refstartpos;
	  refstartidx = pMGi->refstartidx;
	} else {
	  refstart    = max(pMGi->refstartpos, pMGj->refstartpos);
	  refstartidx = max(pMGi->refstartidx, pMGj->refstartidx);
	}
 	if(pMGj->orientforward ? RightExtend : LeftExtend){
	  refstop     = pMGj->refendpos;
	  refstopidx  = pMGj->refendidx;
	} else {
	  refstop = min(pMGi->refendpos, pMGj->refendpos);
	  refstopidx = min(pMGi->refendidx, pMGj->refendidx);
	}
      }

      // NOTE : qrygap is NOT assumed to be part of tandem duplication region 
      qrystart = min(maxq1,maxq2) - (refstop - refstart);
      qrystop = max(minq1,minq2) + (refstop - refstart);
      
      /* Find nearest label on reference and query for each position */
      refstartidx = findlabel(refstart * 1e-3, Y1, N1);
      refstopidx = findlabel(refstop * 1e-3, Y1, N1);
      qrystartidx = findlabel(qrystart * 1e-3, X, M);
      qrystopidx = findlabel(qrystop * 1e-3, X, M);

      double qrygapstart = min(maxq1,maxq2);
      double qrygapstop = max(minq1,minq2);
      int qrygapstartidx = findlabel(qrygapstart * 1e-3, X , M);
      int qrygapstopidx = findlabel(qrygapstop * 1e-3, X, M);
      if(DEBUG && !qryoverlap) assert(qrygapstart <= qrygapstop);
      if(DEBUG && !qryoverlap) assert(qrygapstartidx <= qrygapstopidx);
      double DupSpacing =  max(0.0, X[qrygapstopidx - 1] - X[qrygapstartidx + 1]);// in kb

      if(verbose){
	printf("duplication: IleftofJ=%d,i=%d(ref=%0.3f..%0.3f,qry=%0.3f..%0.3f),j=%d(ref=%0.3f..%0.3f,qry=%0.3f..%0.3f): refstart=%0.3f,refend=%0.3f,qrystart=%0.3f,qryend=%0.3f(idx=%d,%d,%d,%d), SpaceQryIdx=%d,%d, Spacing= %0.3f kb\n", QueryIleftofJ,
	       i,pMGi->refstartpos,pMGi->refendpos,pMGi->qrystartpos,pMGi->qryendpos,
	       j,pMGj->refstartpos,pMGj->refendpos,pMGj->qrystartpos,pMGj->qryendpos,
	       refstart,refstop,qrystart,qrystop,refstartidx,refstopidx,qrystartidx,qrystopidx,qrygapstartidx,qrygapstopidx,DupSpacing);
	fflush(stdout);
      }
      
      qrystart = max(0.0, qrystart);
      qrystop = min(X[M+1]*1e3, qrystop);

      if(DupSpacing > smap_Dup_maxSpacing){
	duplication = false;
	complex_sv = true;// defaults to complex_sv : really a duplication with large spacing exceeding smap_Dup_maxSpacing
      }
    }
    else if( duplication ) { //figure out the coords for regular duplication (or inverted duplication)
      refstart    = pMGj->refstartpos;
      refstartidx = pMGj->refstartidx;
      refstop     = min(pMGi->refendpos, pMGj->refendpos);
      refstopidx  = min(pMGi->refendidx, pMGj->refendidx);

      if(qryoverlap){/* need to trim duplication region */
	int refidx1,refidx2;
	double refloc1,refloc2;
	int qryidx1 = qrystartidx, qryidx2 = qrystopidx;

	if(!inversion){
	  /* map overlap query region to reference via right matchgroup j */
	  getRef(pMGj->align, qryidx1, refloc1, refidx1);
	  getRef(pMGj->align, qryidx2, refloc2,refidx2);

	  if(max(refidx1,refidx2) >= refstopidx && min(refidx1,refidx2) <= refstartidx){/* NOT a duplication after all */
	    if(VERB/* HERE HERE >=2 */){
	      printf("  duplication invalid: qryidx=%d,%d, refidx=%d,%d, refstartidx=%d, refstopidx=%d\n",
		     qryidx1,qryidx2, refidx1,refidx2, refstartidx, refstopidx);
	      fflush(stdout);
	    }
	    duplication = false;
	    continue;// WAS33 complex_sv = true;// defaults to complex-sv : not really an SV, just indication of same duplication in both reference and query
	  } else {/* trim duplication region to exclude the query overlap region mapped to ref by matchgroup j*/
	    if(max(refidx1,refidx2) > refstartidx && min(refidx1,refidx2) <= refstartidx){
	      if(VERB/* HERE HERE >=2 */){
		printf("  modified dup coords: qryidx=%d,%d, refidx=%d,%d: refstartidx=%d -> %d, refstopidx=%d, refstart= %0.3f -> %0.3f, refstop= %0.3f\n",
		       qryidx1,qryidx2, refidx1,refidx2, refstartidx, max(refidx1,refidx2), refstopidx, refstart * 1e-3, max(refloc1,refloc2) * 1e-3, refstop * 1e-3);
		fflush(stdout);
	      }
	      refstartidx = max(refidx1,refidx2);
	      refstart = max(refloc1,refloc2);
	    }
	    if(min(refidx1,refidx2) < refstopidx && max(refidx1,refidx2) >= refstopidx){
	      if(VERB/* HERE HERE >=2 */){
		printf("  modified dup coords: qryidx=%d,%d, refidx=%d,%d: refstartidx=%d, refstopidx=%d -> %d, refstart= %0.3f, refstop= %0.3f -> %0.3f\n",
		       qryidx1,qryidx2, refidx1,refidx2, refstartidx, refstopidx, min(refidx1,refidx2), refstart * 1e-3, refstop * 1e-3, min(refloc1,refloc2) * 1e-3);
		fflush(stdout);
	      }
	      refstopidx = min(refidx1,refidx2);
	      refstop = min(refloc1,refloc2);
	    }
	  }
	} else {/* inverted duplication with query overlap */

	  if(inv_MGi){	  /* map query overlap region to reference via inverted matchgroup i */
	    getRef(pMGi->align, qryidx1, refloc1,refidx1);
	    getRef(pMGi->align, qryidx2, refloc2,refidx2);
	  } else { /* map query overlap region to reference via inverted matchgroup j */
	    getRef(pMGj->align, qryidx1, refloc1,refidx1);
	    getRef(pMGj->align, qryidx2, refloc2,refidx2);
	  }

	  if(max(refidx1,refidx2) > refstartidx && min(refidx1,refidx2) <= refstartidx){
	    if(VERB/* HERE HERE >=2 */){
	      printf("  modified inverted dup coords: qryidx=%d,%d, refidx=%d,%d: refstartidx=%d -> %d, refstopidx=%d, refstart= %0.3f -> %0.3f, refstop= %0.3f\n",
		     qryidx1,qryidx2, refidx1,refidx2, refstartidx, max(refidx1,refidx2), refstopidx, refstart * 1e-3, max(refloc1,refloc2) * 1e-3, refstop * 1e-3);
	      fflush(stdout);
	    }
	    refstartidx = max(refidx1,refidx2);
	    refstart = max(refloc1,refloc2);
	  }
	  if(min(refidx1,refidx2) < refstopidx && max(refidx1,refidx2) >= refstopidx){
	    if(VERB/* HERE HERE >=2 */){
	      printf("  modified inverted dup coords: qryidx=%d,%d, refidx=%d,%d: refstartidx=%d, refstopidx=%d -> %d, refstart= %0.3f, refstop= %0.3f -> %0.3f\n",
		     qryidx1,qryidx2, refidx1,refidx2, refstartidx, refstopidx, min(refidx1,refidx2), refstart * 1e-3, refstop * 1e-3, min(refloc1,refloc2) * 1e-3);
	      fflush(stdout);
	    }
	    refstopidx = min(refidx1,refidx2);
	    refstop = min(refloc1,refloc2);
	  }
	}
      }
      if(DEBUG) assert( (pMGj->qrystartpos < pMGj->qryendpos) == pMGj->orientforward &&
			(pMGi->qrystartpos < pMGi->qryendpos) == pMGi->orientforward );

      //getDuplicationCoords_orig(pMGi, pMGj, refstart, refstop, minrefalignlen,
      //			refoverlapsize, qryoverlap, minq1, minq2, maxq1, maxq2, qrysize, qrystart, qrystop);
      //void getDuplicationCoords(xmapEntry* xe1, xmapEntry* xe2, int refstartidx, int refstopidx,
      //		  double& qrystart, int& qrystartidx, double& qrystop, int& qrystopidx) {
      double DupSpacing;
      if(DEBUG) assert(refstartidx <= refstopidx);
      getDuplicationCoords(pMGi, pMGj, refstartidx, refstopidx, qrystart, qrystartidx, qrystop, qrystopidx, DupSpacing, verbose);

      if(verbose){
        printf("  duplication (inv=%d): ref interval: %.3f %.3f, qry interval: %.3f %.3f, Dup Spacing= %0.3f kb(max=%0.3f)\n", 
	       inversion ? 1 : 0, refstart * 1e-3, refstop * 1e-3, qrystart*1e-3, qrystop*1e-3, DupSpacing, (inversion ? smap_InvDup_maxSpacing : smap_Dup_maxSpacing));
	fflush(stdout);
      }
      if(DEBUG) assert(qrystop > qrystart && refstop > refstart);
      //      if(DEBUG) assert(DupSpacing <= max(0.0, (qrystop - qrystart - 2.0 * (refstop - refstart))*1e-3 + 1e-3));

      qrystart = max(0.0, qrystart);
      qrystop = min(X[M+1]*1e3, qrystop);

      if(!qryoverlap && DupSpacing > (inversion ? smap_InvDup_maxSpacing : smap_Dup_maxSpacing)){
	if(DEBUG>=1+RELEASE) assert(!(inversion && duplication));/* this case should have been ruled out previously by checking DupSpacing */
	inversion = false;/* NEW : the inversion case was ruled out previously */
	duplication = false;
	complex_sv = true;// defaults to complex-sv : really a duplication with large spacing exceeding smap_Dup_maxSpacing
      }
    }

    else if(!duplication && inversion) {// set coords for inversion 
      if(inv_MGi){/* MG i is the inverted matchgroup */
	if(pMGi->orientforward){/* check order of matchgroups i & j on inverted query */
	  MGi_qryleft = (pMGi->qryendpos > pMGj->qrystartpos);/* inverted MG i starts before MG j on inverted query */
	  // MGi_qryleft = (pMGi->qryendpos - pMGj->qrystartpos) > (pMGj->qryendpos - pMGi->qrystartpos);/* inverted MG i has more query range before (than after) MG j on inverted query */
	} else {/* check order of matchgroups i & j on un-inverted query */
	  MGi_qryleft = (pMGi->qryendpos < pMGj->qrystartpos);/* inverted MG i starts before MG j on non-inverted query */
	  // MGi_qryleft = (pMGj->qrystartpos - pMGi->qryendpos) > (pMGi->qrystartpos - pMGj->qryendpos);/* inverted MG i has more query range before (than after) MG j on non-inverted query */
	}

	if(MGi_qryleft){
	  refstop     = pMGj->refstartpos;
	  refstopidx  = pMGj->refstartidx;
	  if(pMGj->refstartidx < pMGi->refendidx){/* matchgroups are overlaped on reference : don't include overlapped region of inverted matchgroup i in inversion breakpoint interval */
	    refstart    = refstop;
	    refstartidx = refstopidx;
	  } else {
	    refstart    = pMGi->refendpos;
	    refstartidx = pMGi->refendidx;
	  }
	  if(qryoverlap){/* adjust query gap interval to match the boundary of non-inverted matchgroup j */
	    if(pMGj->orientforward){
	      qrystart = qrystop = minq2;
	      qrystartidx = qrystopidx = pMGj->qrystartidx;
	    } else {
	      qrystop = qrystart = maxq2;
	      qrystopidx = qrystartidx = pMGj->qrystartidx;
	    }
	    qrysize = 0.0;
	  }
	} else { // MGi_qryleft == false 
	  if(DEBUG) assert(pMGj->refendidx <= pMGi->refendidx);/* matchgroup j is completely overlapped by matchgroup i on ref : exclude overlapped ref region of inverted matchgroup i */
	  // NOTE : if above assertion fails, there are additional cases to handle

	  // Effectively trimmed inverted matchgroup i is treated as being right of non-inverted matchgroup j

	  refstart = pMGj->refendpos;
	  refstartidx = pMGj->refendidx;

	  refstop = refstart;
	  refstopidx = refstartidx;

	  if(qryoverlap){/* adjust query interval to match the boundary of the non-inverted matchgroup j */
	    if(pMGj->orientforward){
	      qrystop = qrystart = maxq2;
	      qrystopidx = qrystartidx = pMGj->qryendidx;
	    } else {
	      qrystart = qrystop = minq2;
	      qrystartidx = qrystopidx = pMGj->qryendidx;
	    }
	    qrysize = 0.0;
	  }
	}

      } else {/* MG j is the inverted matchgroup */
	if(pMGi->orientforward){/* check order of matchgroups i & j on un-inverted query */
	  MGi_qryleft = (pMGi->qrystartpos < pMGj->qryendpos);/* MG i starts before inverted MG j on non-inverted query */
	} else {/* check order of matchgroups i & j on inverted query */
	  MGi_qryleft = (pMGi->qrystartpos > pMGj->qryendpos);/* MG i starts before inverted MG j on inverted query */
	}

	if(MGi_qryleft){
	  refstart    = pMGi->refendpos;
	  refstartidx = pMGi->refendidx;
	  if(pMGj->refstartidx < pMGi->refendidx){/* matchgroups are overlaped on reference : don't include overlapped region of inverted matchgroup j in inversion breakpoint interval */
	    refstop     = refstart;
	    refstopidx  = refstartidx;
	  } else {
	    refstop     = pMGj->refstartpos;
	    refstopidx  = pMGj->refstartidx;
	  }
	  if(qryoverlap){/* adjust query interval to match the boundary of non-inverted matchgroup i */
	    if(pMGi->orientforward){
	      qrystop = qrystart = maxq1;
	      qrystopidx = qrystartidx = pMGi->qryendidx;
	    } else {
	      qrystart = qrystop = minq1;
	      qrystartidx = qrystopidx = pMGi->qryendidx;
	    }
	    qrysize = 0.0;
	  }
	} else { // NEW18 : MGi_qryleft == false : matchgroup i is fully overlapped by matchgroup j on qry
	  if(DEBUG) assert(pMGi->refendidx < pMGj->refendidx); /* matchgroup j cannot be completely overlapped by matchgroup i on both ref AND qry ! */
	  
	  refstart = pMGi->refendpos;
	  refstartidx = pMGi->refendidx;
	  if(pMGi->refendidx > pMGj->refstartidx){/* matchgroups are overlapped on reference (should be rare) : trim off part of inverted matchgroup j that is overlapped on ref */
	    refstop = refstart;
	    refstopidx = refstartidx;
	  } else {
	    refstop = pMGj->refstartpos;
	    refstopidx = pMGj->refstartidx;
	  }

	  if(DEBUG) assert(qryoverlap);
	  if(DEBUG && !(qrystart > qrystop)){
	    printf("qrystart= %0.4f, qrystop= %0.4f kb\n",qrystart*1e-3,qrystop*1e-3);
	    fflush(stdout);
	    assert(qrystart > qrystop);
	  }

	  /* trim off part of inverted matchgroup j corresponding to query overlap interval */
	  if(pMGi->orientforward){
	    qrystop = qrystart = maxq1;// NEW101
	    qrystopidx = qrystartidx = imaxq1;// NEW101
	  } else {
	    qrystart = qrystop = minq1;// NEW101
	    qrystartidx = qrystopidx = iminq1;// NEW101
	  }
	  qrysize = 0.0;
	}
      }
      if(verbose){
	printf("  inversion: inv_MGi=%d, MGi_qryleft=%d: ref interval=%d,%d(%.3f %.3f), qry interval=%d,%d(%.3f %.3f)\n", 
	       inv_MGi, MGi_qryleft, refstartidx,refstopidx,refstart * 1e-3, refstop * 1e-3, qrystartidx,qrystopidx,qrystart*1e-3, qrystop*1e-3);
	fflush(stdout);
      }
    } // if(!duplication && inversion)

    else  if( /* WAS rearrangement || */ translocation ) { //same for both intra- and inter-chromosomal translocations
      // case where smaller MG overlaps indel in larger MG was handled previously

      //note <=, >= required to prevent assert in else, but these can't be translocations because of overlap criteria
      if( ( minq1 <= minq2 && pMGi->orientforward && pMGj->orientforward) || 
	  ( minq1 >= minq2 && !pMGi->orientforward && !pMGj->orientforward) ) {
	qrystart    = pMGi->qryendpos;
	qrystop     = pMGj->qrystartpos;
	qrystartidx = pMGi->qryendidx;
	qrystopidx  = pMGj->qrystartidx;
	refstart    = pMGi->refendpos;
	refstop     = pMGj->refstartpos;
	refstartidx = pMGi->refendidx;
	refstopidx  = pMGj->refstartidx;
      }
      else if( ( minq1 <= minq2 && !pMGi->orientforward && !pMGj->orientforward ) || 
	       ( minq1 >= minq2 && pMGi->orientforward && pMGj->orientforward ) ) {
	qrystart    = pMGi->qrystartpos;
	qrystop     = pMGj->qryendpos;
	qrystartidx = pMGi->qrystartidx;
	qrystopidx  = pMGj->qryendidx;
	refstart    = pMGi->refstartpos;
	refstop     = pMGj->refendpos;	
	refstartidx = pMGi->refstartidx;
	refstopidx  = pMGj->refendidx;	
      }
      else if( ( minq1 <= minq2 && pMGi->orientforward && !pMGj->orientforward ) ||
	       ( minq1 >= minq2 && !pMGi->orientforward && pMGj->orientforward ) ) {
	qrystart    = pMGi->qryendpos;
	qrystop     = pMGj->qryendpos;
	qrystartidx = pMGi->qryendidx;
	qrystopidx  = pMGj->qryendidx;
	refstart    = pMGi->refendpos;
	refstop     = pMGj->refendpos; 
	refstartidx = pMGi->refendidx;
	refstopidx  = pMGj->refendidx; 
      }
      else if( ( minq1 <= minq2 && !pMGi->orientforward && pMGj->orientforward ) ||
	       ( minq1 >= minq2 && pMGi->orientforward && !pMGj->orientforward ) ) {
	qrystart    = pMGi->qrystartpos;
	qrystop     = pMGj->qrystartpos;
	qrystartidx = pMGi->qrystartidx;
	qrystopidx  = pMGj->qrystartidx;
	refstart    = pMGi->refstartpos;
	refstop     = pMGj->refstartpos; 
	refstartidx = pMGi->refstartidx;
	refstopidx  = pMGj->refstartidx; 
      }
      else {
	printf("ERROR: breakpoint error for translocation\n\n");
	assert(false);
      }
      if(verbose){
	printf("  translocation: refid=%lld,%lld,qryid=%lld,%lld: qry = %0.3f .. %0.3f, ref= %0.3f .. %0.3f kb, translocation=%d\n",
	       pMGi->refcontigid,pMGj->refcontigid,pMGi->qrycontigid,pMGj->qrycontigid,
	       qrystart,qrystop,refstart,refstop, translocation);
	fflush(stdout);
      }
    } // if(rearrangement || translocation)
    else if(complex_sv){
      // coordinates are computed later, since some of the other cases also revert to complex_sv
    }
    else if (overlap_indel){// coordinates for this case were computed previously 
      if(verbose){
	printf("  Overlapped Indel: refid=%lld,%lld,qryid=%lld,%lld: qry = %0.3f .. %0.3f, ref= %0.3f .. %0.3f kb, qrysiz= %0.3f, refsiz= %0.3f\n",
	       pMGi->refcontigid,pMGj->refcontigid,pMGi->qrycontigid,pMGj->qrycontigid,
	       qrystart*1e-3,qrystop*1e-3,refstart*1e-3,refstop*1e-3,qrysize,refsize);
 	fflush(stdout);
      }
      if(DEBUG) assert(refstart <= refstop);
      if(DEBUG) assert(qrystart <= qrystop);
      if(DEBUG) assert(fabs((qrystop - qrystart)*1e-3 - qrysize) < 1e-6);
      if(DEBUG) assert(fabs((refstop - refstart)*1e-3 - refsize) < 1e-6);
    }
    else { // default is 2 MG indel : all other cases must be handled previously
      refstart    = pMGi->refendpos;
      refstartidx = pMGi->refendidx;
      refstop     = pMGj->refstartpos;
      refstopidx  = pMGj->refstartidx;
      //	refsize = (refstop - refstart) * 1e-3; //size is in kb, positions are bp
      if(verbose){
	printf("  2 MG Indel: refid=%lld,%lld,qryid=%lld,%lld: qry = %0.3f .. %0.3f, ref= %0.3f .. %0.3f kb, qrysize= %0.3f, refsize= %0.3f: logPV= %0.2f,%0.2f, logPV2= %0.2f,%0.2f\n",
	       pMGi->refcontigid,pMGj->refcontigid,pMGi->qrycontigid,pMGj->qrycontigid,
	       qrystart*1e-3,qrystop*1e-3,refstart*1e-3,refstop*1e-3,qrysize,refsize, pMGi->align->logPV,pMGj->align->logPV,pMGi->align->logPV2,pMGj->align->logPV2);
	fflush(stdout);
      }
      if(DEBUG) assert(!qryoverlap);
      if(DEBUG) assert(!refoverlap);
      if(DEBUG) assert(refstart <= refstop);
      if(DEBUG) assert(refstartidx <= refstopidx);
      if(DEBUG) assert(qrystart <= qrystop);
      if(DEBUG) assert(fabs((qrystop - qrystart)*1e-3 - qrysize) < 1e-6);
      if(DEBUG) assert(fabs((refstop - refstart)*1e-3 - refsize) < 1e-6);
    }

    // NOTE : complex_sv may have been set to true in above if(){} else {} block

    if(complex_sv){

      //for complex, call from start of first MG to max of ends of both matchgroups
      refstart    = pMGi->refstartpos; //starts are guaranteed sorted
      refstartidx = pMGi->refstartidx;
      refstop  = max(pMGi->refendpos,  pMGj->refendpos);
      refstopidx  = (pMGi->refendpos > pMGj->refendpos ?
		     pMGi->refendidx : pMGj->refendidx);

      if(verbose){
	printf("  complex: refid=%lld,%lld,qryid=%lld,%lld: qry = %0.3f .. %0.3f, ref= %0.3f .. %0.3f kb\n",
	       pMGi->refcontigid,pMGj->refcontigid,pMGi->qrycontigid,pMGj->qrycontigid,
	       qrystart,qrystop,refstart,refstop);
	fflush(stdout);
      }
    } // complex types
    
    // Next determine svtype_enum : default of insertion or deletion/compression was set previously
    if( translocation ) {
      if( intra_rearr )
	svtype_enum = intrachr_translocation;
      else
	svtype_enum = interchr_translocation;
    } else if(rearrangement || intra_rearr){
      svtype_enum = complex_type;
    } else if( duplication /* NEW33 */|| duplication_split) { 
      svtype_enum = (inversion ? duplicationinv_type : duplication_type);
      svtype_enum = (duplication_split ? duplicationsplit_type : svtype_enum);
    } //if none of above and inversion found, it's true inversion (except the !intra_rearr is not checked above)
    else if( inversion /* WAS30 && !intra_rearr */ ) { 
      if(DEBUG) assert(!intra_rearr);
      svtype_enum = inversion_type;
    } else if((qryoverlap || refsize < 0 || qrysize < 0) && !overlap_indel){
      if(DEBUG) assert(complex_sv);
      svtype_enum = complex_type;
    } else if(opporient){
      if(DEBUG) assert(!inversion);
      if(DEBUG) assert(complex_sv);
      svtype_enum = complex_type;
    } else if( complex_sv){
      svtype_enum = complex_type;
    }

    if(DEBUG) assert(svtype_enum != unknown_type);

    const char *tyst = 0; //for debugging--more info than 'complex'
    if(verbose) { //only printed for verbose
      //note that this if structure sets precendence for these naming conventions
      // ie, overlaps will be called before bad inversions
      if( translocation && intra_rearr )
	tyst = "transintra"; //translocation, shortened for printout
      else if( translocation && !intra_rearr )
	tyst = "transinter"; //translocation, shortened for printout
      else if( rearrangement && !translocation )
	tyst = "inter-complex"; //inter-chromosomal complex
      else if( intra_rearr && !translocation )
	tyst = "intra-complex"; //intra-chromosomal complex
      else if(duplication_split)
	tyst = "duplication-split";
      else if(duplication && !inversion)
	tyst = "duplication";
      else if(duplication && inversion)
	tyst = "inverted-duplication";
      else if(!duplication && inversion)
	tyst = "inversion";
      else if( qryoverlap && !overlap_indel )
	tyst = "qryoverlap-complex"; //overlap on qry is complex
      else if( refsize < 0 && !overlap_indel )
	tyst = "refoverlap-complex"; //overlap on ref is complex
      else if( qrysize < 0 && !overlap_indel )
	tyst = "qrynegsiz-complex"; //this happens for out-of-order matches with no overlap
      else if( opporient && !inversion )
	tyst = "badinv-complex"; //inversion which fails smap_inversion_ratio (is complex)
#if 0
      else if( qrysize > smap_sv_sizekb )
	tyst = "qrylarge-complex";
      else if( refsize > smap_sv_sizekb )
	tyst = "reflarge-complex";
#endif
    }

    if(verbose && VERB /* HERE HERE >= 2 */){
      printf("SV main loop indp=%d/%d: inv=%d(inv_only=%d),dup=%d,dup_split=%d,overlap_indel=%d,complex_sv=%d : xmapids=%d,%d, svtype_enum=%d, tyst=%s: sv_array[%d]\n",
	     indp,indiceslen,inversion,inversion_only,duplication,duplication_split,overlap_indel, complex_sv, pMGi->xmapid, pMGj->xmapid, svtype_enum, tyst ? tyst : "null", numsv);
      fflush(stdout);
    }

    if(inversion_only && !inversion)
      continue;// One of the two matchgroups would have been filtered due to overlap except for -InversionOverlapFilter, which is only intended to save inversion calls

    /* NEW88 : filter out calls other than small inversion or small inverted duplication where both matchgroups don't satisfy -T (for 2-MG indels only based on logPV2 or -indelMinConf) */
    if(svtype_enum == insertion || svtype_enum == deletion){
      if(min(pMGi->align->logPV2, pMGj->align->logPV2) < LogPvThreshold || max(pMGi->align->logPV,pMGj->align->logPV) < InDelMinConf){
	if(verbose){
	  printf("i=%d,j=%d:xmapids=%d,%d:svtype=%d:logPV=%0.2f,%0.2f(indelMinConf=%0.2f),logPV2=%0.2f,%0.2f(LogPvThreshold=%0.2f): filtered out 2-MG indel due to final MG confidence\n",
		 i,j,pMGi->xmapid,pMGj->xmapid,svtype_enum,pMGi->align->logPV,pMGj->align->logPV,InDelMinConf,pMGi->align->logPV2,pMGj->align->logPV2,LogPvThreshold);
	  fflush(stdout);
	}
	continue;
      }
    } else if((svtype_enum == inversion_small_type || svtype_enum == inversion_type || svtype_enum == duplicationinv_type) && sv_ps == sv_indel){
    } else {
      if(min(pMGi->confidence, pMGj->confidence) < LogPvThreshold){
	if(verbose){
	  printf("i=%d,j=%d:xmapids=%d,%d:svtype=%d:conf=%0.2f,%0.2f,LogPvThreshold=%0.2f : filtered out SV due to final MG confidence\n",
		 i,j,pMGi->xmapid,pMGj->xmapid,svtype_enum,pMGi->confidence,pMGj->confidence,LogPvThreshold);
	  fflush(stdout);
	}
	continue;
      }
    }

    //printf("bot loop %2i %2i (%2i %2i): qrystart/stop: %.4f %.4f; refstart/top: %.4f %.4f\n", pMGi->xmapid, pMGj->xmapid, i, j, qrystart, qrystop, refstart, refstop); //debug -- always gave same result as top loop (but may not for inter-chr) -- now enabling above logic for translocations

    //create structuralVariation object, store indices
    //old index was indp, which is also index of this loop
    //new index is numsv, global for size of sv_arry
    //always assign left, right conf since those are valid for all types; center is only valid for indels, so exclude that and set later
    if(numsv >= maxsv)
      maxsv = growArray(sv_array,maxsv, max(numsv * 3 / 2, numsv + 17));

    sv_array[numsv] = structuralVariation(i, j, qrysize, refsize, qrystart, qrystop, pMGi->refcontigid, pMGj->refcontigid, refstart, refstop, tyst, svtype_enum, qrystartidx, qrystopidx, refstartidx, refstopidx, -1., pMGi->confidence, pMGj->confidence, sv_ps, pMGi->xmapid, pMGj->xmapid);
    sv_array[numsv].inv_MGi = inv_MGi;
    sv_array[numsv].MGi_qryleft = MGi_qryleft;
    sv_array[numsv].repeat = (pMGi->have_ref_repeat || pMGj->have_ref_repeat); //move to method of sv class

    if(SV_ORIENT && translocation){
      if(pMGi->orientforward){
	if(pMGj->orientforward){
	  sv_array[numsv].orientation = strdup("+/+");
	} else {
	  if(minq1 < minq2)
	    sv_array[numsv].orientation = strdup("+/-");	    
	  else
	    sv_array[numsv].orientation = strdup("-/+");	    
	}
      } else {
	if(pMGj->orientforward){
	  if(minq1 < minq2)
	    sv_array[numsv].orientation = strdup("-/+");	    
	  else
	    sv_array[numsv].orientation = strdup("+/-");	    
	} else {
	  sv_array[numsv].orientation = strdup("-/-");
	}
      }
    }

    //if(DEBUG) assert(pMGi->qrycontigid == pMGj->qrycontigid);
    numsv++;
    if(DEBUG) assert(numsv <= maxsv);

    if(VERB>=2){
      printf("SV main loop indp=%d/%d: end of loop, numsv=%d,maxsv=%d\n",indp,indiceslen,numsv,maxsv);
      fflush(stdout);
    }
  } // END for indp = 0 .. indiceslen-1 (main sv loop)

  if(VERB){
    printf("Finished main sv loop: indiceslen=%d, numsvIndel=%d, numsv=%d,maxsv=%d, repeatcnt=%d, repeatcnt2=%lu, swapcnt=%d,trimcnt=%d (wall time=%.6f)\n", 
	   indiceslen, numsvIndel,numsv,maxsv,repeatcnt,cumrepeatcnt2, swapcnt,trimcnt,wtime());
    fflush(stdout);
  }

  if(DEBUG>=2){
    for(int i = 0; i < smapsize;i++){
      xmapEntry *p = usexmapentries[i];

      /* check p->refsites[] (bp) vs refendpos (bp) */
      assert(p->numrefsites >= p->refendidx - p->refstartidx + 1);
      int Istart = p->refstartidx;
      for(int I = Istart; I <= p->refendidx; I++){
	double sitesI = p->refsites[I - Istart];
	if(DEBUG && !(p->refstartpos <= sitesI && sitesI <= p->refendpos)){
	  printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refstartpos= %0.1f,refendpos= %0.1f\n",
		 i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart, p->refsites[I-Istart],p->refstartpos,p->refendpos);
	  fflush(stdout);
	  assert(p->refstartpos <= sitesI && sitesI <= p->refendpos);
	}
	if(I > p->refstartidx){
	  if(DEBUG && !(p->refsites[I-1 - Istart] <= sitesI)){
	    printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numrefsites=%d:I=%d,refsites[%d]= %0.1f, refsites[%d]= %0.1f, refendpos= %0.1f\n",
		   i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numrefsites,I,I-Istart,p->refsites[I - Istart],I-1-Istart,p->refsites[I-1-Istart],p->refendpos);
	    fflush(stdout);
	    assert(p->refsites[I-1-Istart] <= sitesI);
	  }
	}
      }
      /* check p->qrysites[0..numqrysites-1] (bp) vs refendpos (bp) */
      double qrystartpos = min(p->qrystartpos,p->qryendpos);
      double qryendpos = max(p->qrystartpos,p->qryendpos);
      int Imin = min(p->qrystartidx,p->qryendidx);
      int Imax = max(p->qrystartidx,p->qryendidx);
      if(DEBUG) assert(p->numqrysites >= Imax-Imin+1);
      for(int I = Imin; I <= Imax; I++){
	double sitesI = p->qrysites[I - Imin];
	if(DEBUG && !(qrystartpos <= sitesI && sitesI <= qryendpos)){
	  printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		 i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites, I, I-Imin,p->qrysites[I-Imin],qrystartpos,qryendpos);
	  fflush(stdout);
	  assert(qrystartpos <= sitesI && sitesI <= qryendpos);
	}
	if(I > Imin){
	  if(DEBUG && !(p->qrysites[I-1-Imin] <= sitesI)){
	    printf("usexmapentries[i=%d/%d]:refid=%lld,qryid=%lld,fwd=%d,numqrysites=%d:I=%d,qrysites[%d]= %0.1f, qrysites[%d]= %0.1f, qrystartpos= %0.1f, qryendpos= %0.1f\n",
		   i,smapsize, p->refcontigid,p->qrycontigid,p->orientforward,p->numqrysites,I,I-Imin,p->qrysites[I-Imin],I-1-Imin,p->qrysites[I-1-Imin],qrystartpos,qryendpos);
	    fflush(stdout);
	    assert(p->qrysites[I-1-Imin] <= sitesI);
	  }
	}
      }
    }
  }

  if(DEBUG>=2){/* check if any 1-MG Indels no longer overlap the MGs they refer to */
    //    int cnt = 0;
    for(int i = 0; i < numsvIndel; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum != sv_ps);
      if(DEBUG && i < numsvIndel1) assert(p->indices[0] == p->indices[1]);

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;
      
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      if(DEBUG) assert(entry1->output == true);
      if(entry1->output == false)
	continue;

      double qmin = min(entry1->qrystartpos,entry1->qryendpos);
      double qmax = max(entry1->qrystartpos,entry1->qryendpos);
      int iqmin = min(entry1->qrystartidx,entry1->qryendidx);
      int iqmax = max(entry1->qrystartidx,entry1->qryendidx);

      if(DEBUG>=2 && !((qmin <= p->querystart && p->querystart <= qmax) &&
	   (qmin <= p->querystop && p->querystop <= qmax) &&
	   (iqmin <= p->querystartidx && p->querystartidx <= iqmax) &&
	   (iqmin <= p->querystopidx && p->querystopidx <= iqmax) &&
	   (entry1->refstartpos <= p->refstart && p->refstart <= entry1->refendpos) &&
	   (entry1->refstartpos <= p->refstop && p->refstop <= entry1->refendpos) &&
	   (entry1->refstartidx <= p->refstartidx && p->refstartidx <= entry1->refendidx) &&
	   (entry1->refstartidx <= p->refstopidx && p->refstopidx <= entry1->refendidx))){
	if(i < numsvIndel1)
	  printf("WARNING:sv_array[%d]:Single MG indel's MG was trimmed and no longer includes the SV coordinates:\n",i);
	else 
	  printf("WARNING:sv_array[%d]: 2-MG small inversion or duplication's large MG was trimmed and no longer includes the SV coordinates:\n",i);
	printf("\t sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d) : deleting SV\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("\t xmapid1=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	if(i >= numsvIndel1){
	  xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
	  printf("\t xmapid2=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry2->xmapid, entry2->refcontigid,entry2->qrycontigid,
	       entry2->qrystartpos * 1e-3, entry2->qryendpos * 1e-3, entry2->qrystartidx, entry2->qryendidx, 
	       entry2->refstartpos * 1e-3, entry2->refendpos * 1e-3, entry2->refstartidx, entry2->refendidx,
	       entry2->orientforward, entry2->confidence, entry2->align->numpairs);
	}
	printf("\t numsvIndel1=%d,numsvIndel2=%d,numsvIndel=%d: deleting SV\n",numsvIndel1,numsvIndel2,numsvIndel);
	fflush(stdout);

	//	if(DEBUG>=1+RELEASE && p->indices[0] != p->indices[1]) assert(0);
      }
      if(i < numsvIndel1 ? !(p->type_enum == insertion || p->type_enum == deletion) :
	 i < numsvIndel ? !(p->type_enum == inversion_small_type || p->type_enum == duplication_type || p->type_enum == duplicationinv_type) : 0){
	printf("sv_array[%d]: type_enum=%d, numsvIndel1=%d,numsvIndel2=%d, numsvIndel=%d\n",i,p->type_enum,numsvIndel1,numsvIndel2,numsvIndel);
	printf("\t sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("\t xmapid=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	fflush(stdout);
	assert(0);
      }
    }
  }

  if(DEBUG>=2)
    for(int i = 0; i < smapsize; i++)
      assert(usexmapentries[i]->numrefsites >= usexmapentries[i]->refendidx - usexmapentries[i]->refstartidx + 1);

  if(MATCHGROUP_PARTIAL_TRIM && !repeatcnt && (swapcnt > 0 || trimcnt > 0)){/* repeat main SV loop one more time, since matchgroups have been trimmed */
    repeatcnt++;

    qsort(usexmapentries, smapsize, sizeof(usexmapentries[0]), (intcmp *)xmapEntryQryidInc);
    indiceslen = getSameQueryIndices(usexmapentries, smapsize, indices1, indices2, indices_size);    

    if(VERB){
      printf("sv array: numsv=%d -> %d (numsvIndel), indiceslen=%i, maxsv=%i : repeating main SV loop: wall time= %0.6f secs\n", numsv, numsvIndel, indiceslen, maxsv,wtime()); fflush(stdout); //debug
      fflush(stdout);
    }

    numsv = numsvIndel;

    swapcnt = 0;

    goto LmainRepeat;
  }

  if(DEBUG>=2){
    for(int i = numsvIndel; i < numsv - 1; i++){
      xmapEntry *entry1 = usexmapentries[sv_array[i].indices[0]];
      xmapEntry *entry2 = usexmapentries[sv_array[i+1].indices[0]];
      if(!(entry1->qrycontigid <= entry2->qrycontigid && entry1->qrycontigid == usexmapentries[sv_array[i].indices[1]]->qrycontigid)){
	for(int j = i; j <= i+1; j++)
	  printf("i=%d,sv_array[i].indices[0,1]= %d,%d, usexmapentries[sv_array[i].indices[0,1]]->qrycontigid= %lld,%lld\n",
		 j,sv_array[j].indices[0],sv_array[j].indices[1],usexmapentries[sv_array[j].indices[0]]->qrycontigid,usexmapentries[sv_array[j].indices[1]]->qrycontigid);
	fflush(stdout);
	assert(entry1->qrycontigid <= entry2->qrycontigid);
	assert(entry1->qrycontigid == usexmapentries[sv_array[i].indices[1]]->qrycontigid);
      }
      assert(entry2->qrycontigid == usexmapentries[sv_array[i+1].indices[1]]->qrycontigid);
    }
  }
  if(DEBUG>=2){
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum == sv_ps || p->algo_enum == sv_indel);
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsvIndel=%d, numsv=%d\n", numsvIndel,numsv);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
      if(p->algo_enum == sv_ps && !(p->xmapid1 != p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsv=%d\n", numsv);
	fflush(stdout);
	assert(p->xmapid1 != p->xmapid2);
      }
    }
  }

  delete [] firstIndel; firstIndel = NULL;

  /* NEW88 : filter out small duplications where both matchgroups don't satisfy -T */
  /*         Also filter out small inversion or small inverted duplication if outer matchgroup does not satisfy -T or inner matchgroups does not satisfy LogPvThreshold3 */
  //  for(int t = numsvIndel1; t < numsvIndel; t++){
  for(int t = 0; t < numsv; t++){
    structuralVariation *p = &sv_array[t];
    if(p->algo_enum != sv_indel)
      continue;
    if(p->indices[0] == p->indices[1])
      continue;
    if(DEBUG) assert(numsvIndel1 <= t && t < numsvIndel);

    if(DEBUG) assert(p->algo_enum == sv_indel);
    if(DEBUG) assert(p->indices[0] != p->indices[1]);
    xmapEntry *Large = xmapentries[p->indices[0]];
    xmapEntry *Small = xmapentries[p->indices[1]];
    if(p->type_enum == duplication_type){
      if(min(Large->confidence, Small->confidence) < LogPvThreshold){
	if(verbose){
	  printf("sv_array[%d]:xmapids=%d,%d:svtype=%d:conf=%0.2f,%0.2f,LogPvThreshold=%0.2f : filtered out small duplication due to final MG confidence\n",
		 t,Large->xmapid,Small->xmapid,p->type_enum,Large->confidence,Small->confidence,LogPvThreshold);
	  fflush(stdout);
	}
	p->duplicate = true;
      }
      continue;
    }
    
    if(Large->confidence < LogPvThreshold || Small->confidence < LogPvThreshold3){
      if(verbose){
	printf("sv_array[%d]:xmapids=%d,%d:svtype=%d:conf=%0.2f,%0.2f,LogPvThreshold=%0.2f,%0.2f : filtered out SV due to final MG confidence\n",
	       t,Large->xmapid,Small->xmapid,p->type_enum,Large->confidence,Small->confidence,LogPvThreshold,LogPvThreshold3);
	fflush(stdout);
      }
      p->duplicate = true;
    }
  }

  /* NEW89 : filter output any SV based on matchgroups that fail smap_confbylen or smap_ConfByInterval */
  if(verbose && (smap_confbylen > 0.0 || smap_ConfByInterval > 0.0)){
    printf("Re-checking all %d SVs for matchgroup confidence values based on ConfByLen=% 0.4f, ConfByInterval= %0.2f\n",numsv, smap_confbylen, smap_ConfByInterval);
    fflush(stdout);
  }
  for(int t = 0; t < numsv; t++){
    structuralVariation *p = &sv_array[t];

    xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
    xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];

    double Len1 = min( fabs(entry1->qrystartpos - entry1->qryendpos),
		       fabs(entry1->refstartpos - entry1->refendpos) ) / 1000.; 

    double Len2 = min( fabs(entry2->qrystartpos - entry2->qryendpos),
		       fabs(entry2->refstartpos - entry2->refendpos) ) / 1000.; 
    if(entry1->confidence < smap_confbylen * Len1 || entry1->confidence < smap_ConfByInterval * (entry1->align->numpairs - 1.0) ||
       entry2->confidence < smap_confbylen * Len2 || entry2->confidence < smap_ConfByInterval * (entry2->align->numpairs - 1.0)){
      if(verbose){
	if(entry1 == entry2)
	  printf("sv_array[%d]:xmapids=%d:svtype=%d:conf=%0.2f,Len=%0.4f kb, pairs=%d : ConfByLen=%0.4f,ConfByInterval=%0.2f : filtered out SV due to final MG confidence\n",
		 t,entry1->xmapid,p->type_enum,entry1->confidence,Len1, entry1->align->numpairs, smap_confbylen, smap_ConfByInterval);
	else
	  printf("sv_array[%d]:xmapids=%d,%d:svtype=%d:conf=%0.2f,%0.2f,Len=%0.4f,%0.4f kb, pairs=%d,%d : ConfByLen=%0.4f,ConfByInterval=%0.2f : filtered out SV due to final MG confidence\n",
		 t,entry1->xmapid,entry2->xmapid,p->type_enum,entry1->confidence,entry2->confidence,Len1,Len2,entry1->align->numpairs,entry2->align->numpairs,smap_confbylen,smap_ConfByInterval);
	fflush(stdout);
      }
      p->duplicate = true;
    }
  }

  /* NEW406 : Filter out Insertion if a 3rd MG overlaps query gap by at least svInsMaxGapOverlap (default 50%)
     Filter out Deletion if a 3rd MG overlaps reference gap by at least svDelMaxGapOverlap (default 50%) */
  int dupcnt = 0, totcnt = 0;
  long long lastqrycontigid = -1;
  int kmin = 0,kmax = 0;

  for(int i= numsvIndel; i < numsv; i++) {
    if(sv_array[i].type_enum != insertion && sv_array[i].type_enum != deletion && sv_array[i].type_enum != compression)
      continue; //only indels

    totcnt++;
    if(sv_array[i].type_enum == insertion && svInsMaxGapOverlap <= 0.0)
      continue;

    if(sv_array[i].type_enum != insertion && svDelMaxGapOverlap <= 0.0)
      continue;

    if(DEBUG>=2) assert(sv_array[i].algo_enum == sv_ps);
    if(DEBUG>=2) assert(sv_array[i].indices[0] < sv_array[i].indices[1]);

    xmapEntry *p = usexmapentries[sv_array[i].indices[0]];
    xmapEntry *q = usexmapentries[sv_array[i].indices[1]];
    long long qrycontigid = p->qrycontigid;
    long long refcontigid = p->refcontigid;
    if(DEBUG>=2) assert(q->qrycontigid == qrycontigid);
    if(DEBUG>=2) assert(q->refcontigid == refcontigid);

    double minr1 = min(p->refstartpos, p->refendpos);
    double maxr1 = max(p->refstartpos, p->refendpos);
    double minr2 = min(q->refstartpos, q->refendpos);
    double maxr2 = max(q->refstartpos, q->refendpos);

    double minq1 = min(p->qrystartpos, p->qryendpos);
    double maxq1 = max(p->qrystartpos, p->qryendpos);
    double minq2 = min(q->qrystartpos, q->qryendpos);
    double maxq2 = max(q->qrystartpos, q->qryendpos);

    if(qrycontigid != lastqrycontigid){/* locate range of usexmapentries[kmin..kmax] with same qrycontigid, to speed up search for 3rd MG*/
      kmin = min(sv_array[i].indices[0],sv_array[i].indices[1]);
      kmax = max(sv_array[i].indices[0],sv_array[i].indices[1]);
      for(;--kmin >= 0;){
	xmapEntry *pMGk = usexmapentries[kmin];	
	if(pMGk->qrycontigid != qrycontigid){
	  kmin++;
	  break;
	}
      }
      kmin = max(0,kmin);

      for(;++kmax < smapsize;){
	xmapEntry *pMGk = usexmapentries[kmax];
	if(pMGk->qrycontigid != qrycontigid){
	  kmax--;
	  break;
	}
      }
      kmax = min(smapsize-1, kmax);
    }
    
    int k = kmin;
    for(; k <= kmax; k++){
      xmapEntry *r = usexmapentries[k];
      if(r==p || r==q)
	continue;
      if(r->refcontigid != refcontigid)
	continue;
      if(r->confidence < LogPvThreshold)
	continue;


      if(sv_array[i].type_enum == deletion || sv_array[i].type_enum == compression){/* check for overlap with reference gap */
	double minr3 = min(r->refstartpos, r->refendpos);
	double maxr3 = max(r->refstartpos, r->refendpos);
	if(VERB>=2){
	  printf("xmapid=%d,%d:qryid=%lld,refid=%lld (ref=%0.3f..%0.3f) : xmapid=%d(ref=%0.3f..%0.3f) overlaps ref gap by %0.1f%% (-svDelMaxGapOverlap): ignoring Deletion\n",
		 p->xmapid,q->xmapid,  qrycontigid,refcontigid,maxr1*1e-3,minr2*1e-3,r->xmapid,minr3*1e-3,maxr3*1e-3,svDelMaxGapOverlap*100.0);
	  fflush(stdout);
	}
	if(minr2 > maxr1 ? overlapLength(maxr1,minr2,minr3,maxr3) >= (minr2-maxr1) * svDelMaxGapOverlap :
	   minr1 > maxr2 ? overlapLength(maxr2,minr1,minr3,maxr3) >= (minr1-maxr2) * svDelMaxGapOverlap : false){
	  if(verbose){
	    printf("xmapid=%d,%d:qryid=%lld,refid=%lld : xmapid=%d(ref=%0.3f..%0.3f) overlaps ref gap by %0.1f%% (-svDelMaxGapOverlap): ignoring Deletion\n",
		   p->xmapid,q->xmapid,  qrycontigid,refcontigid,r->xmapid,minr3*1e-3,maxr3*1e-3,svDelMaxGapOverlap*100.0);
	    fflush(stdout);
	  }
	  break;
	}
      } else {/* check for overlap with qry gap */
	double minq3 = min(r->qrystartpos, r->qryendpos);
	double maxq3 = max(r->qrystartpos, r->qryendpos);
	if(minq2 > maxq1 ? overlapLength(maxq1,minq2,minq3,maxq3) >= (minq2-maxq1) * svInsMaxGapOverlap :
	   minq1 > maxq2 ? overlapLength(maxq2,minq1,minq3,maxq3) >= (minq1-maxq2) * svInsMaxGapOverlap : false){
	  if(verbose){
	    printf("xmapid=%d,%d:qryid=%lld,refid=%lld : xmapid=%d(qry=%0.3f..%0.3f) overlaps qry gap by %0.1f%% (-svInsMaxGapOverlap): ignoring Insertion\n",
		   p->xmapid,q->xmapid,  qrycontigid,refcontigid,r->xmapid,minq3*1e-3,maxq3*1e-3,svInsMaxGapOverlap*100.0);
	    fflush(stdout);
	  }
	  break;
	}
      }
    }// for(k = kmin; k <= kmax; k++)
    if(k <= kmax){
      dupcnt++;
      sv_array[i].duplicate = true;
    }
  }// int i = numsvIndel; i < numsv; i++

  if(VERB/* HERE >=2 */){
    printf("finished Indel gap overlap check : found %d Indels with 3rd MG overlapping Indel gap : wall time= %0.6f secs\n",dupcnt,wtime());
    fflush(stdout);
  }

  // check for duplicate indels ie, cominatorial duplicates: 1-2, 1-3, 2-3: 1-3 is duplicate
  // because sv_array only compares same query & ref contigs, the xmapentries pointed to by indices
  // must also, so this is guaranteed to be same query & ref contig (same for other duplicates)

  /* HERE HERE : 1-2, 1-3 without 2-3 should also filter out 1-3 IFF 2 & 3 are NOT crossed (assume 2-3 is a different SV type, eg duplication)
     HERE HERE : 1-3, 2-3 without 1-2 should also filter out 1-3 IFF 1 & 2 are NOT crossed (assume 1-3 is a different SV type, eg duplication)
     (should not be needed if -svGapOverlapRatio 0.7 is used) */

  // NOTE : Only checks sv_array[i = numsvIndel .. numsv -1], since 0..numsvIndel -1 are Indels with algo_enum == sv_indel
  if(DEBUG>=2)
    for(int i = 0; i < numsvIndel; i++)
      assert(sv_array[i].algo_enum == sv_indel);

  // NOTE : Assumes that usexmapentries[] for same qrycontigid & refcontigid is sorted in ascending order of refstartpos, refendpos, otherwise it is not clear which of the 3 Indels is redundant
  // The usexmapentries[] qsort() above sorts qrycontigid, then refcontigid then refstartpos then refendpos
  // NEW faster code relies on fact that sv_array[numsvIndel .. numsv-1] is already sorted in ascending order of indices[0], then indices[1] due to getSameQueryIndices() using usexmapentries[] ordering 

  dupcnt = 0;

  for(int i= /* WAS 1 */ numsvIndel + 1; i < numsv - 1; i++) { //first and last can't be duplicate
    bool firstused = false;
    bool secondused = false;
    if( sv_array[i].type_enum != insertion && sv_array[i].type_enum != deletion && sv_array[i].type_enum != compression )
      continue; //only indels

    if(DEBUG>=2) assert(sv_array[i].algo_enum == sv_ps);
    //    if( sv_array[j].algo_enum != sv_ps ) //duplicate checking only applies to pairsplit results
    //	continue;

    if(DEBUG>=2) assert(sv_array[i].indices[0] < sv_array[i].indices[1]);

    // compare first indices to previous first indices
    for(int j = i; --j >= numsvIndel;){
      if( sv_array[j].type_enum != insertion && sv_array[j].type_enum != deletion && sv_array[j].type_enum != compression && (!INDEL_FILTER_FIX || sv_array[j].type_enum != duplication_type))
	continue; //only indels or non-inverted duplication

      if(DEBUG>=2) assert(sv_array[j].algo_enum == sv_ps);
      if(DEBUG>=2) assert(sv_array[j].indices[0] < sv_array[j].indices[1]);
      if(DEBUG>=2) assert(sv_array[j].indices[0] <= sv_array[i].indices[0]);

      if(sv_array[i].indices[0] == sv_array[j].indices[0]){
	if(DEBUG>=2) assert(sv_array[j].indices[1] < sv_array[i].indices[1]);
	firstused = true;
      }
      break;
    }
    
    if(!firstused)
      continue;

    // compare second index to following second
    for(int j = i+1; j < numsv; j++){
      if( sv_array[j].type_enum != insertion && sv_array[j].type_enum != deletion && sv_array[j].type_enum != compression && (!INDEL_FILTER_FIX || sv_array[j].type_enum != duplication_type))
	continue; //only indels or non-inverted duplication

      if(DEBUG>=2) assert(sv_array[j].algo_enum == sv_ps);
      if(DEBUG>=2) assert(sv_array[j].indices[0] >= sv_array[i].indices[0]);
      if(DEBUG>=2) assert(sv_array[j].indices[0] < sv_array[j].indices[1]);

      if(sv_array[j].indices[0] >= sv_array[i].indices[1])
	break;/* since (with condition in assertion in previous line) it is no longer possible for sv_array[j].indices[1] to equal sv_array[i].indices[1] */

      if(sv_array[i].indices[1] == sv_array[j].indices[1] ){
	if(DEBUG>=2) assert(sv_array[i].indices[0] < sv_array[j].indices[0]);
	secondused = true;
	break;
      }
    }

    if( firstused && secondused ) {
      if( VERB && verbose ){
	printf("Duplicate indel sv_array[%2i] (xmapid= %2i %2i)\n", i,
	       usexmapentries[sv_array[i].indices[0]]->xmapid,
	       usexmapentries[sv_array[i].indices[1]]->xmapid);
	fflush(stdout);
      }
      if(!sv_array[i].duplicate) dupcnt++;

      sv_array[i].duplicate = true;
    }

    if(0 && firstused){// found 1-2 (j) and 1-3 (i) : If the matchgroups 2 & 3 are not crossed, delete SV 1-3
      if(DEBUG) assert(!secondused);
      // HERE HERE : should not be needed if -svGapOverlapRatio 0.7 is used
    }
    if(0 && secondused){// found 1-3 (i) and 2-3 (j) : If the matchgroups 1 & 2 are not crossed, delete SV 1-3
      if(DEBUG) assert(!firstused);
      // HERE HERE : should not be needed if -svGapOverlapRatio 0.7 is used
    }
  } //end duplicate indel outer loop

  if(VERB/* HERE >=2 */){
    printf("finished duplicate indel loop : found %d duplicate indels : wall time= %0.6f secs\n", dupcnt, wtime());
    fflush(stdout);
  }
  if(DEBUG>=2){
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum == sv_ps || p->algo_enum == sv_indel);
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsvIndel=%d, numsv=%d\n", numsvIndel,numsv);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
      if(p->algo_enum == sv_ps && !(p->xmapid1 != p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsv=%d\n", numsv);
	fflush(stdout);
	assert(p->xmapid1 != p->xmapid2);
      }
    }
  }

  //call duplicate function for inversions and both kinds of translocations
  checkDuplicates( sv_array, numsvIndel, numsv, usexmapentries, inversion_type, 0 /* verbose */ );
  if(VERB/* HERE >=2 */){
    printf("finished checkDuplicates for inversions : wall time= %0.6f secs\n", wtime());
    fflush(stdout);
  }
  if(DEBUG>=2){
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum == sv_ps || p->algo_enum == sv_indel);

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsvIndel=%d, numsv=%d\n", numsvIndel,numsv);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
      if(p->algo_enum == sv_ps && !(p->xmapid1 != p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsv=%d\n", numsv);
	fflush(stdout);
	assert(p->xmapid1 != p->xmapid2);
      }
    }
  }

  checkDuplicates( sv_array, numsvIndel, numsv, usexmapentries, intrachr_translocation, 0 );
  if(VERB/* HERE >=2 */){
    printf("finished checkDuplicates for intrachr translocations : wall time= %0.6f secs\n", wtime());
    fflush(stdout);
  }
  if(DEBUG>=2){
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum == sv_ps || p->algo_enum == sv_indel);
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsvIndel=%d, numsv=%d\n", numsvIndel,numsv);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
      if(p->algo_enum == sv_ps && !(p->xmapid1 != p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsv=%d\n", numsv);
	fflush(stdout);
	assert(p->xmapid1 != p->xmapid2);
      }
    }
  }

  checkDuplicates( sv_array, numsvIndel, numsv, usexmapentries, interchr_translocation, 0 );
  if(VERB/* HERE >=2 */){
    printf("finished checkDuplicates for interchr translocations : wall time= %0.6f secs\n", wtime());
    fflush(stdout);
  }

  bool verbose_trans = false; //print all pairs of translocations and filtering results--very long
  int ntrans = 0; //get number of translocations
  for(int i=0; i < numsv; i++) {
    if( sv_array[i].type_enum == interchr_translocation || sv_array[i].type_enum == intrachr_translocation )
      ntrans++;
  }
  if( verbose && ntrans > 0 ){ //don't bother if zero bc it's redundant with more verbose list below
    printf("N translocations: %i\n", ntrans);
    fflush(stdout);
  }

  // Thomas : Speed up by only checking i = numsvIndel .. numsv -1, since 0..numsvIndel -1 are single match group Indels with algo_enum == sv_indel (either Indel or small Inversion)
  // NOTE : SV's (other than sv_indel) are currently sorted by the indices[0] then indices[1]. Since usexmapentries[] are sorted by qrycontid
  // before the sv_array[numsvindel .. numsv -1] were generated, this means that these SV's are sorted in ascending order of qrycontigid (assuming each SV refers to only one qrycontigid)
  if(DEBUG>=2){
    for(int i = numsvIndel; i < numsv - 1; i++){
      xmapEntry *entry1 = usexmapentries[sv_array[i].indices[0]];
      xmapEntry *entry2 = usexmapentries[sv_array[i+1].indices[0]];
      if(!(entry1->qrycontigid <= entry2->qrycontigid && entry1->qrycontigid == usexmapentries[sv_array[i].indices[1]]->qrycontigid)){
	for(int j = i; j <= i+1; j++)
	  printf("i=%d,sv_array[i].indices[0,1]= %d,%d, usexmapentries[sv_array[i].indices[0,1]]->qrycontigid= %lld,%lld\n",
		 j,sv_array[j].indices[0],sv_array[j].indices[1],usexmapentries[sv_array[j].indices[0]]->qrycontigid,usexmapentries[sv_array[j].indices[1]]->qrycontigid);
	fflush(stdout);
	assert(entry1->qrycontigid <= entry2->qrycontigid);
	assert(entry1->qrycontigid == usexmapentries[sv_array[i].indices[1]]->qrycontigid);
      }
      assert(entry2->qrycontigid == usexmapentries[sv_array[i+1].indices[1]]->qrycontigid);
    }
  }
  if(DEBUG>=2){
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum == sv_ps || p->algo_enum == sv_indel);

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsvIndel=%d, numsv=%d\n", numsvIndel,numsv);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
      if(p->algo_enum == sv_ps && !(p->xmapid1 != p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsv=%d\n", numsv);
	fflush(stdout);
	assert(p->xmapid1 != p->xmapid2);
      }
    }
  }

  //double loop over all SVs for overlap -- exclude those with different xmapEntries which overlap
  //note: this block is to check for multiple SVs on the same query contig of the same type which overlap on the query
  // Thomas : this seems to lose perfectly valid indel calls : suppressed for indels for now
  // Thomas : this seems to lose perfectly valid trans calls : suppressed for translocations for now
  for(int i= numsvIndel /* WAS 0 */; i < numsv; i++) {
    if(DEBUG>=2) assert(sv_array[i].algo_enum == sv_ps);
    //    if( sv_array[i].algo_enum != sv_ps ) // overlap checking only applies to "pairsplit" results (meaning SVs that use indices[2] into usexmapentries[])
    //      continue;
    bool indeli = bool(sv_array[i].type_enum == insertion || sv_array[i].type_enum == deletion || sv_array[i].type_enum == compression);
    if(indeli)
      continue;// Thomas
    bool transi = bool(sv_array[i].type_enum == intrachr_translocation ||
		       sv_array[i].type_enum == interchr_translocation);
    if(transi)
      continue;// Thomas

    //bool inveri = bool(sv_array[i].type_enum == inversion_type);
    bool inveri = false;
    if(sv_array[i].type_enum == complex_type || sv_array[i].type_enum == inversion_type || 
       sv_array[i].type_enum == duplication_type || sv_array[i].type_enum == duplicationinv_type ||
       sv_array[i].type_enum == duplicationsplit_type )
      continue;
    assert( indeli || transi || inveri || sv_array[i].type_enum == complex_type ); //be sure all types are checked

    int qrycontigid1 = usexmapentries[sv_array[i].indices[0]]->qrycontigid;

    for(int j=i+1; j < numsv; j++) {
      if(DEBUG>=2) assert(sv_array[j].algo_enum == sv_ps);
      //      if( sv_array[j].algo_enum != sv_ps ) // overlap checking only applies to "pairsplit" results (meaning SVs that use indices[2] into usexmapentries[])
      //	continue;

      int qrycontigid2 = usexmapentries[sv_array[j].indices[0]]->qrycontigid;
      if(qrycontigid2 != qrycontigid1)
	break;// NEW : rely on qrycontig being in ascending order, hence no need to check any larger j values

      //instead of checking overlap only for translocations, check for any two events of the same type
      // skip complex since there's no need for them
      // treat insertions and deletions as the same type, also for two types of translocations
      // and only remaining type which gets overlap is inversions
      bool indelj = bool(sv_array[j].type_enum == insertion || sv_array[j].type_enum == deletion || sv_array[j].type_enum == compression);
      if(indelj)
	continue;// Thomas
      bool transj = bool(sv_array[j].type_enum == intrachr_translocation ||
			 sv_array[j].type_enum == interchr_translocation);
      if(transj)
	continue;// Thomas

      bool inverj = bool(sv_array[j].type_enum == inversion_type);	  
      if( (indeli && (transj || inverj)) ||
	  (transi && (indelj || inverj)) ||
	  (inveri && (indelj || transj)) ||
	  sv_array[j].type_enum == complex_type ||
	  sv_array[j].type_enum == duplication_type || sv_array[j].type_enum == duplicationinv_type ||
	  sv_array[j].type_enum == duplicationsplit_type)
	continue;
      assert( indelj || transj || inverj || sv_array[j].type_enum == complex_type ); 
      //due to -indel svs, these are no longer equivalent, and the indices arrays should not be used
      //assert(indices1[i] == sv_array[i].indices[0]); //sv_array was filled this way--dummy check
      //assert(indices2[i] == sv_array[i].indices[1]);
      //assert(indices1[j] == sv_array[j].indices[0]);
      //assert(indices2[j] == sv_array[j].indices[1]);
      double overlap = 0; //4 potential overlaps, keep first which causes exclusion
      double conf1 = 0, conf2 = 0; //need to get sum of confidence for each SV
      bool exclude = false;

      for(int k=0; k < 2; k++) { //k goes over two xmapEntries which make up SV i
	xmapEntry *entry1 = usexmapentries[sv_array[i].indices[k]];
	double min1 = min(entry1->qrystartpos, entry1->qryendpos);
	double max1 = max(entry1->qrystartpos, entry1->qryendpos);
	conf1 += entry1->confidence;
	for(int l=0; l < 2; l++) { //l goes over two xmapEntries which make up SV j
	  xmapEntry *entry2 = usexmapentries[sv_array[j].indices[l]];
	  if( k == 0 ) //don't double count
	    conf2 += entry2->confidence;
	  if( exclude ) //just take first exclusion 
	    continue; //continue rather than break in order to properly sum confidences
	  //an xmap entry cannot exclude itself; it's allowed to participate in multiple translocations
	  if( entry1->xmapid == entry2->xmapid ) 
	    continue;

	  //for multi-contig case, must check qrycontigid of entries because their indices are not checked
	  if(DEBUG>=2) assert(entry1->qrycontigid == entry2->qrycontigid);
	  //	  if( entry1->qrycontigid != entry2->qrycontigid )
	  //	    continue;
	  double min2 = min(entry2->qrystartpos, entry2->qryendpos);
	  double max2 = max(entry2->qrystartpos, entry2->qryendpos);
	  //check for overlap
	  if( min1 <= min2 && max1 >= min2 ) {
	    if( max1 >= max2 ) //2 is within 1
	      overlap = max2 - min2;
	    else
	      overlap = max1 - min2;
	  }
	  else if( min2 <= min1 && max2 >= min1 ) {
	    if( max2 >= max1 ) //1 is within 2
	      overlap = max1 - min1;
	    else
	      overlap = max2 - min1;
	  }
	  exclude = bool(overlap/min(max1-min1,max2-min2) > smap_translocation_ovlfr); //default 0.3 -- see parameters.cpp
	} //end inner indices loop
      } //end outer indices loop

      //use confs to get index in sv_array to exclude -- note, if conf1 == conf2, do not exclude either
      int excludeidx = 0;
      if( exclude && conf1 > conf2 ) {
	excludeidx = j;
	sv_array[j].exclude = true;
      }
      else if( exclude && conf2 > conf1 ) {
	excludeidx = i;
	sv_array[i].exclude = true;
      }
      else 
	exclude = false;
      
      //report for this SV-pair
      if( verbose_trans ) {
	printf("SV_array[%2i, %2i] : exclude=%d,%d : (xmapids= %2i,%2i, %2i %2i) exclude=%i overlap=%8.1f, conf=%6.2f %6.2f, excludeidx=%2i\n", i, j,
	       sv_array[i].exclude ? 1 : 0, sv_array[j].exclude ? 1 : 0,
	       usexmapentries[sv_array[i].indices[0]]->xmapid,
	       usexmapentries[sv_array[i].indices[1]]->xmapid,
	       usexmapentries[sv_array[j].indices[0]]->xmapid,
	       usexmapentries[sv_array[j].indices[1]]->xmapid,
	       exclude ? 1 : 0, overlap, conf1, conf2, excludeidx);
	fflush(stdout);
      }
    } //end inner numsv loop
  } //end outer numsv loop

  if( verbose ) { //report on excluded
    int nexclude = 0;
    for(int indp=0; indp < numsv; indp++) {
      if(sv_array[indp].exclude) {
	nexclude++;
	if(VERB){
	  if( nexclude == 1 ) //first one only
	    printf("Indices of overlapping SVs:  ");
	  printf("%i  ", indp+1);
	  fflush(stdout);
	}
      }
    }
    if(VERB){
      if( nexclude > 0 )
	printf("\n");
      printf("N SVs excluded due to overlap: %i\n", nexclude);
      fflush(stdout);
    }
  }

  if(VERB/* HERE >=2 */){
    printf("Completed SV duplicate checking, numsv=%d: wall time= %0.6f secs\n", numsv, wtime());
    fflush(stdout);
  }

  if(VERB>=2){
    printf("\nBefore sorting sv_array: numsv=%d,numxmapentries=%d\n",numsv,numxmapentries);
    for(int i = 0; i < numxmapentries; i++){
      xmapEntry *p = usexmapentries[i];
      printf("i=%d/%d: xmapentries[i]:xmapid=%d:refid=%lld,qryid=%lld,fwd=%d:qrystartpos=%0.1f,qryendpos=%0.1f,refstartpos=%0.1f,refendpos=%0.1f,conf=%0.2f\n",
	     i,numxmapentries,p->xmapid,p->refcontigid,p->qrycontigid,p->orientforward ? 1 : 0, p->qrystartpos,p->qryendpos,p->refstartpos, p->refendpos, p->confidence);
      fflush(stdout);
    }
  }
  if(DEBUG>=2){
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum == sv_ps || p->algo_enum == sv_indel);
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];

      if(p->type_enum == unknown_type){
	printf("sv_array[%d]:refid=%lld,%lld,qryid=%lld,%lld:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d:numsvIndel=%d\n",
	       i,entry1->refcontigid,entry2->refcontigid,entry1->qrycontigid,entry2->qrycontigid,
	       p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid,numsvIndel);
	fflush(stdout);
	assert(p->type_enum != unknown_type);
      }

      if(!(p->algo_enum == (i < numsvIndel ? sv_indel : sv_ps))){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d:numsvIndel=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid,numsvIndel);
	fflush(stdout);
	assert(p->algo_enum == (i < numsvIndel ? sv_indel : sv_ps));
      }

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
      if(p->algo_enum == sv_ps && !(p->xmapid1 != p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,type_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,p->type_enum, entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	printf("\t numsv=%d\n", numsv);
	fflush(stdout);
	assert(p->xmapid1 != p->xmapid2);
      }
    }
  }

  if(DEBUG>=2){/* check if any 1-MG Indels no longer overlap the MGs they refer to */
    for(int i = 0; i < numsvIndel; i++){
      structuralVariation *p = &sv_array[i];

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];

      if((DEBUG && !(p->algo_enum != sv_ps && (i >= numsvIndel1 || p->indices[0] == p->indices[1]))) ||
	 (VERB>=3 && (entry1->xmapid == 53590 || entry1->xmapid == 59946 || entry2->xmapid==53590 || entry2->xmapid==59946))){
	printf("sv_array[%d]: type_enum=%d, query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n", i, p->type_enum,
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("\t indices= %d, %d : xmapid= %d,%d\n", p->indices[0], p->indices[1], entry1->xmapid, entry2->xmapid);
	printf("\t xmapid1=%d: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, 
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	printf("\t xmapid2=%d: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry2->xmapid, 
	       entry2->qrystartpos * 1e-3, entry2->qryendpos * 1e-3, entry2->qrystartidx, entry2->qryendidx, 
	       entry2->refstartpos * 1e-3, entry2->refendpos * 1e-3, entry2->refstartidx, entry2->refendidx,
	       entry2->orientforward, entry2->confidence, entry2->align->numpairs);

	fflush(stdout);

	assert(p->algo_enum != sv_ps);
	if(i < numsvIndel1) assert(p->indices[0] == p->indices[1]);
      }

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;
      
      //      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      if(DEBUG) assert(entry1->output == true);
      if(entry1->output == false)
	continue;

      double qmin = min(entry1->qrystartpos,entry1->qryendpos);
      double qmax = max(entry1->qrystartpos,entry1->qryendpos);
      int iqmin = min(entry1->qrystartidx,entry1->qryendidx);
      int iqmax = max(entry1->qrystartidx,entry1->qryendidx);

      if(DEBUG>=2 && !((qmin <= p->querystart && p->querystart <= qmax) &&
	   (qmin <= p->querystop && p->querystop <= qmax) &&
	   (iqmin <= p->querystartidx && p->querystartidx <= iqmax) &&
	   (iqmin <= p->querystopidx && p->querystopidx <= iqmax) &&
	   (entry1->refstartpos <= p->refstart && p->refstart <= entry1->refendpos) &&
	   (entry1->refstartpos <= p->refstop && p->refstop <= entry1->refendpos) &&
	   (entry1->refstartidx <= p->refstartidx && p->refstartidx <= entry1->refendidx) &&
	   (entry1->refstartidx <= p->refstopidx && p->refstopidx <= entry1->refendidx))){
	if(i < numsvIndel1)
	  printf("WARNING:sv_array[%d]:Single MG indel's MG was trimmed and no longer includes the SV coordinates:\n",i);
	else 
	  printf("WARNING:sv_array[%d]: 2-MG small inversion or duplication's large MG was trimmed and no longer includes the SV coordinates:\n",i);
	printf("\t sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("\t xmapid1=%d: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, 
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	if(!(p->type_enum==insertion || p->type_enum==deletion)){
	  xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
	  printf("\t xmapid2=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry2->xmapid, entry2->refcontigid,entry2->qrycontigid,
	       entry2->qrystartpos * 1e-3, entry2->qryendpos * 1e-3, entry2->qrystartidx, entry2->qryendidx, 
	       entry2->refstartpos * 1e-3, entry2->refendpos * 1e-3, entry2->refstartidx, entry2->refendidx,
	       entry2->orientforward, entry2->confidence, entry2->align->numpairs);
	}
	fflush(stdout);
	//	if(DEBUG>=1+RELEASE && p->indices[0] != p->indices[1]) assert(0);
      }
      if(i < numsvIndel1 ? !(p->type_enum == insertion || p->type_enum == deletion) :
	 i < numsvIndel ? !(p->type_enum == inversion_small_type || p->type_enum == duplication_type || p->type_enum == duplicationinv_type) : 0){
	printf("sv_array[%d]: type_enum=%d, numsvIndel1=%d,numsvIndel2=%d, numsvIndel=%d\n",i,p->type_enum,numsvIndel1,numsvIndel2,numsvIndel);
	printf("\t sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n",p->type_enum, 
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("\t xmapid1=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, entry1->refcontigid,entry1->qrycontigid,
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	if(!(p->type_enum==insertion || p->type_enum==deletion)){
	  xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
	  printf("\t xmapid2=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry2->xmapid, entry2->refcontigid,entry2->qrycontigid,
	       entry2->qrystartpos * 1e-3, entry2->qryendpos * 1e-3, entry2->qrystartidx, entry2->qryendidx, 
	       entry2->refstartpos * 1e-3, entry2->refendpos * 1e-3, entry2->refstartidx, entry2->refendidx,
	       entry2->orientforward, entry2->confidence, entry2->align->numpairs);
	}
	fflush(stdout);
	assert(0);
      }
    }
  }

  if(VERB>=3){
    printf("Before sorting SVs in order of xmapid1,xmapid2 etc: numsvIndel=%d, numsv=%d, maxsv=%d\n",numsvIndel,numsv,maxsv);
    for(int i = numsvIndel; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      //      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      
      printf("sv_array[%d]: type=%d, xmapid=%d,%d, rid=%lld,%lld,qryid=%lld\n", i, p->type_enum, p->xmapid1,p->xmapid2,p->refcontigid1,p->refcontigid2,entry1->qrycontigid);
    }
    fflush(stdout);
  }

  // sort SVs in order of xmapid1,xmapid2 etc in smap: also used in svOverlap (???)
  qsort(sv_array, numsv, sizeof(*sv_array), compare_svs);

  if(VERB>=3){
    printf("After sorting SVs in order of xmapid1,xmapid2 etc: numsv=%d, maxsv=%d\n",numsv,maxsv);
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      //      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      
      printf("sv_array[%d]: type=%d, xmapid=%d,%d, rid=%lld,%lld,qryid=%lld\n", i, p->type_enum, p->xmapid1,p->xmapid2,p->refcontigid1,p->refcontigid2,entry1->qrycontigid);
    }
    fflush(stdout);
  }

  if(1){/* check if any 1-MG Indels no longer overlap the MGs they refer to (and filter them out, if needed) */
    int cnt = 0;
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum == sv_ps || p->algo_enum == sv_indel);
      if(p->algo_enum == sv_ps)
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];

      if(DEBUG) assert(p->algo_enum != sv_ps);
      if((DEBUG && (p->type_enum==insertion || p->type_enum==deletion) && !(p->indices[0] == p->indices[1])) ||
	 (VERB>=3 && (entry1->xmapid == 53590 || entry1->xmapid == 59946 || entry2->xmapid==53590 || entry2->xmapid==59946))){
	printf("sv_array[%d]: type_enum=%d, query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n", i, p->type_enum,
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("\t indices= %d, %d : xmapid= %d,%d\n", p->indices[0], p->indices[1], entry1->xmapid, entry2->xmapid);
	printf("\t xmapid1=%d: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, 
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	printf("\t xmapid2=%d: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry2->xmapid, 
	       entry2->qrystartpos * 1e-3, entry2->qryendpos * 1e-3, entry2->qrystartidx, entry2->qryendidx, 
	       entry2->refstartpos * 1e-3, entry2->refendpos * 1e-3, entry2->refstartidx, entry2->refendidx,
	       entry2->orientforward, entry2->confidence, entry2->align->numpairs);

	fflush(stdout);
	if(DEBUG && (p->type_enum==insertion || p->type_enum==deletion))
	  assert(p->indices[0] == p->indices[1]);
      }

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;
      
      //      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      if(DEBUG) assert(entry1->output == true);
      if(entry1->output == false)
	continue;

      double qmin = min(entry1->qrystartpos,entry1->qryendpos);
      double qmax = max(entry1->qrystartpos,entry1->qryendpos);
      int iqmin = min(entry1->qrystartidx,entry1->qryendidx);
      int iqmax = max(entry1->qrystartidx,entry1->qryendidx);

      if(!((qmin <= p->querystart && p->querystart <= qmax) &&
	   (qmin <= p->querystop && p->querystop <= qmax) &&
	   (iqmin <= p->querystartidx && p->querystartidx <= iqmax) &&
	   (iqmin <= p->querystopidx && p->querystopidx <= iqmax) &&
	   (entry1->refstartpos <= p->refstart && p->refstart <= entry1->refendpos) &&
	   (entry1->refstartpos <= p->refstop && p->refstop <= entry1->refendpos) &&
	   (entry1->refstartidx <= p->refstartidx && p->refstartidx <= entry1->refendidx) &&
	   (entry1->refstartidx <= p->refstopidx && p->refstopidx <= entry1->refendidx))){

	if(VERB>=2){
	  if(p->type_enum==insertion || p->type_enum==deletion)
	    printf("WARNING:sv_array[%d]:Single MG indel's MG was trimmed and no longer includes the SV coordinates:\n",i);
	  else 
	    printf("WARNING:sv_array[%d]: 2-MG small inversion or duplication's large MG was trimmed and no longer includes the SV coordinates:\n",i);
	  printf("\t sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d):deleting SV\n",p->type_enum, 
		 p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
		 p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	  printf("\t xmapid1=%d: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, 
		 entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
		 entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
		 entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	  if(!(p->type_enum==insertion || p->type_enum==deletion)){
	    //	  xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
	    printf("\t xmapid2=%d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry2->xmapid, entry2->refcontigid,entry2->qrycontigid,
		   entry2->qrystartpos * 1e-3, entry2->qryendpos * 1e-3, entry2->qrystartidx, entry2->qryendidx, 
		   entry2->refstartpos * 1e-3, entry2->refendpos * 1e-3, entry2->refstartidx, entry2->refendidx,
		   entry2->orientforward, entry2->confidence, entry2->align->numpairs);
	  }
	  fflush(stdout);
	}

	p->duplicate = true;
	cnt++;
	if(p->type_enum == inversion_small_type){/* NEW79 : must change type of paired inversion to inversion_type : this will allow the partial_inversion to be added later */
	  /* located the paired inversion */
	  int jmatch = -1;
	  int matchcnt = 0;
	  for(int j = 0; j < numsv; j++){
	    if(j == i)
	      continue;
	    structuralVariation *q = &sv_array[j];
	    if(q->type_enum != inversion_small_type)
	      continue;
	    if(p->indices[0] == q->indices[0] && p->indices[1] == q->indices[1]){
	      matchcnt++;
	      jmatch = j;
	      if(DEBUG) assert(p->refcontigid1 == q->refcontigid1);
	      if(DEBUG) assert(p->refcontigid2 == q->refcontigid2);
	    }
	  }

	  if(DEBUG && matchcnt > 1){
	    printf("WARNING: Found %d matches for small inversion sv_array[%d] : xmapids=%d,%d: refid= %lld: refstart= %0.1f, refend= %0.1f, qrystart= %0.1f, qryend= %0.1f\n",
		   matchcnt, i, p->indices[0],p->indices[1],p->refcontigid1,p->refstart,p->refstop,p->querystart,p->querystop);
	    for(int j = 0; j < numsv; j++){
	      if(i==j)
		continue;
	      structuralVariation *q = &sv_array[j];
	      if(q->type_enum != inversion_small_type)
		continue;
	      if(p->indices[0] == q->indices[0] && p->indices[1] == q->indices[1])
		printf("\t Match at sv_array[%d] : xmapids=%d,%d, refid=%lld: refstart= %0.1f, refend= %01.f, qrystart= %0.1f, qryend= %0.1f\n",
		       j, q->indices[0],q->indices[1],q->refcontigid1,q->refstart,q->refstop,q->querystart,q->querystop);
	    }
	    fflush(stdout);
	    assert(matchcnt <= 1);
	  }

	  if(matchcnt == 1){
	    structuralVariation *q = &sv_array[jmatch];
	    xmapEntry *pMGi = ((q->algo_enum == sv_ps) ? usexmapentries : xmapentries)[q->indices[0]];
	    xmapEntry *pMGj = ((q->algo_enum == sv_ps) ? usexmapentries : xmapentries)[q->indices[1]];

	    /* make sure pMGi is left of pMGj on reference : this may not be true */
	    if(xmapEntryQryidInc(&pMGi,&pMGj) > 0){
	      int tmp = q->indices[0];
	      q->indices[0] = q->indices[1];
	      q->indices[1] = tmp;

	      pMGi = ((q->algo_enum == sv_ps) ? usexmapentries : xmapentries)[q->indices[0]];
	      pMGj = ((q->algo_enum == sv_ps) ? usexmapentries : xmapentries)[q->indices[1]];
	      q->xmapid1 = pMGi->xmapid;// NEW104
	      q->xmapid2 = pMGj->xmapid;// NEW104

	      if(VERB){
		printf("\t Swapping xmapids : xmapid1-> %d, xmapid2-> %d\n", pMGi->xmapid,pMGj->xmapid);
		fflush(stdout);
	      }
	    }

	    /* need to compute inv_MGi and MGi_qryleft for this inversion */
	    bool inv_MGi = true;
	    if(pMGi->refendpos < pMGj->refendpos){/* MG j is NOT completely overlapped by MG i on reference : orient query so MG i is left of MG j */
	      if (min(pMGi->qrystartpos, pMGi->qryendpos) < min(pMGj->qrystartpos, pMGj->qryendpos) ) { // MG i starts before MG j on query : typically use normal orientation of query
		if (max(pMGi->qrystartpos, pMGi->qryendpos) - max(pMGj->qrystartpos, pMGj->qryendpos) > min(pMGj->qrystartpos, pMGj->qryendpos) - min(pMGi->qrystartpos, pMGi->qryendpos)){
		  // NEW80 : exception : MG i ends after MG j on query AND non-overlapped region of MG i on right side of query is larger than on left side : use reversed orientation of query
		  inv_MGi = pMGi->orientforward;
		} else
		  inv_MGi = !pMGi->orientforward;
	      } else {// MG j starts before MG i on query : typically use reversed orientation of query
		if (max(pMGj->qrystartpos, pMGj->qryendpos) - max(pMGi->qrystartpos, pMGi->qryendpos) > min(pMGi->qrystartpos, pMGi->qryendpos) - min(pMGj->qrystartpos, pMGj->qryendpos)) {
		  // NEW80 : exception : MG j ends after MG i on query AND non-overlapped region on MG j on right side of query is larger than on left side : use normal orientation of query
		  inv_MGi = !pMGi->orientforward;
		} else
		  inv_MGi = pMGi->orientforward;
	      }
	    } /* MG j is completely overlapped by MG i on reference (or both are the same interval), hence cannot be an inversion with MG j inverted */

	    bool MGi_qryleft = true;/* If Matchgroup i is left of Matchgroup j on query */
	    if(inv_MGi){/* MG i is the inverted matchgroup */
	      if(pMGi->orientforward){/* check order of matchgroups i & j on inverted query */
		MGi_qryleft = (pMGi->qryendpos > pMGj->qrystartpos);/* inverted MG i starts before MG j on inverted query */
	      } else {/* check order of matchgroups i & j on un-inverted query */
		MGi_qryleft = (pMGi->qryendpos < pMGj->qrystartpos);/* inverted MG i starts before MG j on non-inverted query */
	      }
	    } else {/* MG j is the inverted matchgroup */
	      if(pMGi->orientforward){/* check order of matchgroups i & j on un-inverted query */
		MGi_qryleft = (pMGi->qrystartpos < pMGj->qryendpos);/* MG i starts before inverted MG j on non-inverted query */
	      } else {/* check order of matchgroups i & j on inverted query */
		MGi_qryleft = (pMGi->qrystartpos > pMGj->qryendpos);/* MG i starts before inverted MG j on inverted query */
	      }
	    }
	    q->inv_MGi = inv_MGi;
	    q->MGi_qryleft = MGi_qryleft;

	    if(VERB){
	      printf("\t sv_array[%d].type_enum = %d -> %d, inv_MGi -> %d, MGi_qryleft -> %d, since paired small inversion sv_array[%d] had to be deleted\n",
		     jmatch,q->type_enum,inversion_type, inv_MGi, MGi_qryleft, i);
	      fflush(stdout);
	    }

	    q->type_enum = inversion_type;

	  }
	}

	//	if(DEBUG>=1+RELEASE && p->indices[0] != p->indices[1]) assert(0);
      }
    }
    if(cnt > 0){
      printf("Deleted %d 1-MG indels or small inversion/duplications\n",cnt);
      fflush(stdout);
    }
  }

  if(DEBUG>=2){
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(DEBUG) assert(p->algo_enum == sv_ps || p->algo_enum == sv_indel);
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];

      if(p->type_enum == unknown_type){
	printf("sv_array[%d]:refid=%lld,%lld,qryid=%lld,%lld:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d:numsvIndel=%d\n",
	       i,entry1->refcontigid,entry2->refcontigid,entry1->qrycontigid,entry2->qrycontigid,
	       p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid,numsvIndel);
	fflush(stdout);
	assert(p->type_enum != unknown_type);
      }
    }
  }

  if(DEBUG>=2){
    for(int i = 0; i < numsv - 1; i++){
      int k;
      if((k = compare_svs(&sv_array[i],&sv_array[i+1])) != -1){
	printf("compare_svs(&sv_array[%d],&sv_array[%d]) = %d : expected -1 (numsv=%d)\n",i,i+1,k,numsv);
	fflush(stdout);
	assert(k== -1);
      }
      if(!(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2))){
	printf("sv_array[%d,%d]:xmapid1=%d,%d,xmapid2=%d,%d (not ordered in ascending order)\n",i,i+1,sv_array[i].xmapid1,sv_array[i+1].xmapid1,sv_array[i].xmapid2,sv_array[i+1].xmapid2);
	fflush(stdout);
	assert(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2));
      }
    }
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];

#if 0
      if(!(p->algo_enum == (i < numsvIndel ? sv_indel : sv_ps))){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,algo_enum=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d:numsvIndel=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid,numsvIndel);
	fflush(stdout);
	assert(p->algo_enum == (i < numsvIndel ? sv_indel : sv_ps));
      }
#endif

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
    }
  }

  /* HERE HERE : eliminate duplicate calls based on the same region of the query using a parsimony based priority ordering
     1. For balanced Indels prefer calls in this order :
        A. (Small) Inversions
	B. Balanced (paired) translocation.
	C. Unbalanced translocation.
	B. (tiny) Indels.
     2. For Insertion prefer calls in this order :
        A. Duplication or Inverted Duplication (Gap size under XXX)
	B. Unbalanced translocation (inter-chromosomal OR Gap size over XXX)
	C. Insertion
     
   */

  for(int indp=0; indp < numsv; indp++) {    //call this in preparation for svOverlap below
    if(VERB>=3){
      structuralVariation *p = &sv_array[indp];
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];

      printf("sv_array[%d]=%p:refid=%lld,%lld,qryid=%lld,%lld:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d:sv_ps=%d, type_enum=%d: computing scaled confidence\n",
	     indp,p,entry1->refcontigid,entry2->refcontigid,entry1->qrycontigid,entry2->qrycontigid, p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->algo_enum == sv_ps ? 1 : 0, p->type_enum);
      fflush(stdout);
    }
    sv_array[indp].wrapperIndelConfidence(usexmapentries, verbose ? 1 : 0); //calls setIndelConfidence
    if(VERB>=3){
      structuralVariation *p = &sv_array[indp];
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];

      printf("sv_array[%d]=%p:refid=%lld,%lld,qryid=%lld,%lld:duplicate=%d,confidence=%0.2f,%0.2f,xmapid1=%d,xmapid2=%d:sv_ps=%d, type_enum=%d: checking for partial Inversions\n",
	     indp,p,entry1->refcontigid,entry2->refcontigid,entry1->qrycontigid,entry2->qrycontigid, p->duplicate?1:0,p->confidence,p->confidence_scaled,p->xmapid1,p->xmapid2,p->algo_enum == sv_ps ? 1 : 0, p->type_enum);
      fflush(stdout);
    }
    sv_array[indp].makeInversionPartial(usexmapentries);
  }
  if(VERB/* HERE >=2 */){
    printf("Completed SV sorting and overlap prep, numsv=%d: wall time= %0.6f secs\n", numsv,wtime());
    fflush(stdout);
  }
  if(VERB>=3){
    printf("SVs after SV sorting:\n");
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      //      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      
      printf("sv_array[%d]: type=%d, xmapid=%d,%d, rid=%lld,%lld,qryid=%lld, ref= %0.3f .. %0.3f kb\n", 
	     i, p->type_enum, p->xmapid1,p->xmapid2,p->refcontigid1,p->refcontigid2,entry1->qrycontigid,p->refstart*1e-3,p->refstop*1e-3);
    }
    printf("\nusexmapentries[0..%d]:\n", smapsize-1);
    for(int i = 0; i < smapsize; i++){
      xmapEntry *pxmap = usexmapentries[i];
      printf("usexmapentries[%d]: qryid=%lld,refid=%lld: qry= %0.3f .. %0.3f, ref= %0.3f .. %0.3f kb\n",
	     i,pxmap->qrycontigid,pxmap->refcontigid,pxmap->qrystartpos*1e-3,pxmap->qryendpos*1e-3,pxmap->refstartpos*1e-3,pxmap->refendpos*1e-3);
    }
    fflush(stdout);
  }

#if LATE_ZYGOSITY == 0  // NOTE :  zygosity calling should be based on all MGs that will be output in final XMAP, not just on currently active MGs in usexmapentries[0..smapsize-1] 

  //call fns in structuralVariation.h/cpp to check for alignment, sv overlaps for inversion zygosity calling
  doPairInversionBP(usexmapentries, sv_array, numsvIndel, numsv, false);
  if(VERB/* HERE >=2 */){
    printf("Completed doPairInversionBP: wall time= %0.6f secs\n", wtime());
    fflush(stdout);
  }
  alignmentOverlap(usexmapentries, smapsize);
  alignmentOverlapInversion(usexmapentries, smapsize);
  alignmentOverlapTranslocation(usexmapentries, smapsize);
  svOverlap(usexmapentries, smapsize);

  if(VERB){
    printf("Finished sv zygosity, numsv=%d (wall time=%.6f)\n", numsv,wtime());
    fflush(stdout);
  }
#endif

  if(DEBUG>=2){
    for(int i = 0; i < numsv - 1; i++){
      int k;
      if((k = compare_svs(&sv_array[i],&sv_array[i+1])) != -1){
	printf("compare_svs(&sv_array[%d],&sv_array[%d]) = %d : expected -1 (numsv=%d)\n",i,i+1,k,numsv);
	fflush(stdout);
	assert(k== -1);
      }
      if(!(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2))){
	printf("sv_array[%d,%d]:xmapid1=%d,%d,xmapid2=%d,%d (not ordered in ascending order)\n",i,i+1,sv_array[i].xmapid1,sv_array[i+1].xmapid1,sv_array[i].xmapid2,sv_array[i+1].xmapid2);
	fflush(stdout);
	assert(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2));
      }
    }
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
    }
  }

  doTransInvConf(usexmapentries, true); //translocations
  doTransInvConf(usexmapentries, false); //inversions

  if(VERB){
    printf("Finished Translocation & Inversion confidence (wall time= %0.6f)\n",wtime());
    fflush(stdout);
  }
  if(DEBUG>=2){
    for(int i = 0; i < numsv - 1; i++){
      int k;
      if((k = compare_svs(&sv_array[i],&sv_array[i+1])) != -1){
	printf("compare_svs(&sv_array[%d],&sv_array[%d]) = %d : expected -1 (numsv=%d)\n",i,i+1,k,numsv);
	fflush(stdout);
	assert(k== -1);
      }
      if(!(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2))){
	printf("sv_array[%d,%d]:xmapid1=%d,%d,xmapid2=%d,%d (not ordered in ascending order)\n",i,i+1,sv_array[i].xmapid1,sv_array[i+1].xmapid1,sv_array[i].xmapid2,sv_array[i+1].xmapid2);
	fflush(stdout);
	assert(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2));
      }
    }
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
    }
  }

  if(DEBUG/* HERE >=2 */){/* check if any 1-MG Indels no longer overlap the MGs they refer to */
    int cnt = 0;
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(p->algo_enum == sv_ps)
	continue;

      if(DEBUG) assert(p->algo_enum != sv_ps);
      if(DEBUG && (p->type_enum==insertion || p->type_enum==deletion) && !(p->indices[0] == p->indices[1])){
	xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
	xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
	printf("sv_array[%d]: type_enum=%d, query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d)\n", i, p->type_enum,
	       p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
	       p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	printf("\t indices= %d, %d : xmapid= %d,%d\n", p->indices[0], p->indices[1], entry1->xmapid, entry2->xmapid);
	printf("xmapid=%d: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, 
	       entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
	       entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
	       entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	printf("xmapid=%d: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry2->xmapid, 
	       entry2->qrystartpos * 1e-3, entry2->qryendpos * 1e-3, entry2->qrystartidx, entry2->qryendidx, 
	       entry2->refstartpos * 1e-3, entry2->refendpos * 1e-3, entry2->refstartidx, entry2->refendidx,
	       entry2->orientforward, entry2->confidence, entry2->align->numpairs);

	fflush(stdout);
	assert(p->indices[0] == p->indices[1]);
      }

      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;
      
      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      if(DEBUG) assert(entry1->output == true);
      if(entry1->output == false)
	continue;

      double qmin = min(entry1->qrystartpos,entry1->qryendpos);
      double qmax = max(entry1->qrystartpos,entry1->qryendpos);
      int iqmin = min(entry1->qrystartidx,entry1->qryendidx);
      int iqmax = max(entry1->qrystartidx,entry1->qryendidx);

      if(!((qmin <= p->querystart && p->querystart <= qmax) &&
	   (qmin <= p->querystop && p->querystop <= qmax) &&
	   (iqmin <= p->querystartidx && p->querystartidx <= iqmax) &&
	   (iqmin <= p->querystopidx && p->querystopidx <= iqmax) &&
	   (entry1->refstartpos <= p->refstart && p->refstart <= entry1->refendpos) &&
	   (entry1->refstartpos <= p->refstop && p->refstop <= entry1->refendpos) &&
	   (entry1->refstartidx <= p->refstartidx && p->refstartidx <= entry1->refendidx) &&
	   (entry1->refstartidx <= p->refstopidx && p->refstopidx <= entry1->refendidx))){
	if(VERB>=2){
	  if(p->type_enum==insertion || p->type_enum==deletion)
	    printf("WARNING:sv_array[%d]:Single MG indel's MG was trimmed and no longer includes the SV coordinates:\n",i);
	  else 
	    printf("WARNING:sv_array[%d]: 2-MG small inversion or duplication's large MG was trimmed and no longer includes the SV coordinates:\n",i);
	  printf("\t sv: type_enum=%d,query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d) : deleting SV\n",p->type_enum, 
		 p->querystart * 1e-3, p->querystop * 1e-3, p->querystartidx, p->querystopidx,
		 p->refstart * 1e-3, p->refstop * 1e-3, p->refstartidx, p->refstopidx);
	  printf("\t xmapid1= %d: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry1->xmapid, 
		 entry1->qrystartpos * 1e-3, entry1->qryendpos * 1e-3, entry1->qrystartidx, entry1->qryendidx, 
		 entry1->refstartpos * 1e-3, entry1->refendpos * 1e-3, entry1->refstartidx, entry1->refendidx,
		 entry1->orientforward, entry1->confidence, entry1->align->numpairs);
	  if(!(p->type_enum==insertion || p->type_enum==deletion)){
	    xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
	    printf("\t xmapid2= %d:rid=%lld,qid=%lld: query= %0.4f .. %0.4f kb (%d..%d), ref= %0.4f .. %0.4f kb (%d..%d): fwd=%d, conf= %0.2f, np=%d\n", entry2->xmapid, entry2->refcontigid,entry2->qrycontigid,
		   entry2->qrystartpos * 1e-3, entry2->qryendpos * 1e-3, entry2->qrystartidx, entry2->qryendidx, 
		   entry2->refstartpos * 1e-3, entry2->refendpos * 1e-3, entry2->refstartidx, entry2->refendidx,
		   entry2->orientforward, entry2->confidence, entry2->align->numpairs);
	  }
	  fflush(stdout);
	}

	p->duplicate = true;
	cnt++;

	if(DEBUG && p->indices[0] != p->indices[1] && p->type_enum == inversion_small_type) assert(0);
      }
    }
    if(cnt > 0){
      printf("Deleted %d 1-MG indels or small inversion/duplications\n",cnt);
      fflush(stdout);
    }
  }

  // sort SVs in ascending order of type_enum and SVsize
  //  qsort(sv_array, numsv, sizeof(*sv_array), (intcmp *)SVenumSVsizeInc);

  xmapidMax = 0;
  for(int i = 0; i < numxmapentries; i++){
    if(DEBUG) assert(xmapentries[i]->xmapid > 0);
    xmapidMax = max(xmapidMax, xmapentries[i]->xmapid);
  }
  int *xmapindexMap = new int[xmapidMax+1];/* map from xmapid to xmapentries[] index (used for debugging) */
  memset(xmapindexMap, -1, xmapidMax * sizeof(int));
  for(int i = 0; i < numxmapentries; i++)
    xmapindexMap[xmapentries[i]->xmapid] = i;

  /* set pair_index for inversion_small_type */
  for(int i = 0; i < numsv; i++){
    structuralVariation *pSV = &sv_array[i];
    if(pSV->type_enum != inversion_small_type)
      continue;
    if(pSV->duplicate)
      continue;
    if(pSV->pair_index >= 0)
      continue;

    int jmatch = -1;
    int matchcnt = 0;// make sure there is only 1 matching small inversion 

    for(int j = i+1; j < numsv; j++){
      structuralVariation *qSV = &sv_array[j];
      if(qSV->type_enum != inversion_small_type)
	continue;
      if(qSV->pair_index >= 0)
	continue;
      if(pSV->indices[0] == qSV->indices[0] && pSV->indices[1] == qSV->indices[1]){
	matchcnt++;
	jmatch = j;
	if(DEBUG) assert(pSV->refcontigid1 == qSV->refcontigid1);
	if(DEBUG) assert(pSV->refcontigid2 == qSV->refcontigid2);
      }
    }
    
    if(DEBUG && matchcnt != 1){
      printf("WARNING: Found %d matches for small inversion sv_array[%d] : xmapids=%d,%d: refid= %lld: refstart= %0.1f, refend= %0.1f, qrystart= %0.1f, qryend= %0.1f, duplicate=%d\n",
	     matchcnt, i, pSV->indices[0],pSV->indices[1],pSV->refcontigid1,pSV->refstart,pSV->refstop,pSV->querystart,pSV->querystop,pSV->duplicate);
      for(int j = i+1; j < numsv; j++){
	structuralVariation *qSV = &sv_array[j];
	if(qSV->type_enum != inversion_small_type)
	  continue;
	if(qSV->pair_index >= 0)
	  continue;
	if(pSV->indices[0] == qSV->indices[0] && pSV->indices[1] == qSV->indices[1])
	  printf("\t Match at sv_array[%d] : xmapids=%d,%d, refid=%lld: refstart= %0.1f, refend= %01.f, qrystart= %0.1f, qryend= %0.1f,duplicate=%d\n",
		 j, qSV->indices[0],qSV->indices[1],qSV->refcontigid1,qSV->refstart,qSV->refstop,qSV->querystart,qSV->querystop,qSV->duplicate);
      }
      fflush(stdout);
      if(MATCHGROUP_PARTIAL_TRIM && matchcnt==0){
	printf("\t Deleted small inversion\n");
	fflush(stdout); 
	sv_array[i].duplicate = true;
	continue;
      }
      assert(matchcnt == 1);
    }
    
    sv_array[jmatch].pair_index = i;
    sv_array[i].pair_index = jmatch;

    if(DEBUG && abs(i-jmatch) > 1){
      printf("WARNING: paired inversions are not located next to each other in sv_array:\n");
      printf("\t sv_array[%d].pair_index = %d, duplicate=%d, confidence= %0.4f\n",i, jmatch, sv_array[i].duplicate,sv_array[i].confidence);
      printf("\t sv_array[%d].pair_index = %d, duplicate=%d, confidence= %0.4f\n",jmatch, i, sv_array[jmatch].duplicate,sv_array[jmatch].confidence);
    }

    if( 0 || verbose) { //don't think we need this one
      printf("\t sv_array[%d].pair_index = %d, duplicate=%d, confidence= %0.4f\n",i, jmatch, sv_array[i].duplicate,sv_array[i].confidence);
      printf("\t sv_array[%d].pair_index = %d, duplicate=%d, confidence= %0.4f\n",jmatch, i, sv_array[jmatch].duplicate,sv_array[jmatch].confidence);
    }
  }

  if(DEBUG>=2){
    for(int i = 0; i < numsv - 1; i++){
      int k;
      if((k = compare_svs(&sv_array[i],&sv_array[i+1])) != -1){
	printf("compare_svs(&sv_array[%d],&sv_array[%d]) = %d : expected -1 (numsv=%d)\n",i,i+1,k,numsv);
	fflush(stdout);
	assert(k== -1);
      }
      if(!(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2))){
	printf("sv_array[%d,%d]:xmapid1=%d,%d,xmapid2=%d,%d (not ordered in ascending order)\n",i,i+1,sv_array[i].xmapid1,sv_array[i+1].xmapid1,sv_array[i].xmapid2,sv_array[i+1].xmapid2);
	fflush(stdout);
	assert(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2));
      }
    }
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	fflush(stdout);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
    }
  }

  //open output file
  FILE *fp;
  printf("Generating %s\n",filename);
  if((fp = fopen(filename,"w"))==NULL){
    printf("failed to open file %s for writing smap file\n",filename);
    exit(1);
  }

  /* write out commandline */
  printversion(fp);

  char absbasename[BUFSIZ];

#ifndef WIN32
  /* get current directory pathname */
  char cwd[PATH_MAX];
  if(getcwd(cwd,PATH_MAX) == NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("getcwd() failed (errno=%d):%s\n",eno,err);
    exit(1);
  }

  /* convert basename into absolute pathname */
  if(basename[0]=='/')
    strcpy(absbasename,basename);
  else
    sprintf(absbasename,"%s/%s",cwd,basename);
  if(queryfilename && queryfilename[0] != '/'){/* convert queryfilename into an absolute pathname */
    char absname[BUFSIZ];
    sprintf(absname,"%s/%s",cwd,queryfilename);
    free(queryfilename);
    queryfilename = strdup(absname);
  }
#else
  strcpy(absbasename,basename);
#endif

  //sorting of alignment should remain after output_xmap call

  //versions:
  //0.0 was indels only; 0.1 has newer types
  //0.2 removes orientation and has many other changes
  //0.3 adds LinkID (bool smap_linkid_col must be true)
  //0.4 adds indices in cmaps of positions reported in {Qry,Ref}{Start,End}Pos columns
  //0.5 adds zygosity
  //0.6 adds RawConfidence column at end
  //0.7 extends RawConfidence to RawConfidenceLeft RawConfidenceRight and RawConfidenceCenter
  //0.8 adds SVsize (for Indels, Inversions and Duplications, value is -1 for all other SVs)
  //0.9 adds SVfreq (-1 for complex SVs)
  //fprintf(fp,"# SMAP File Version:\t%s\n", (smap_labelind_col ? "0.4" : "0.3")); 
  //fprintf(fp,"# SMAP File Version:\t%s\n", (smap_zygosity ? "0.5" : "0.4")); 

  if(SV_FREQ || SV_ORIENT)
    fprintf(fp,"# SMAP File Version:\t%s\n", "0.9"); 
  else if(SV_SIZE)
    fprintf(fp,"# SMAP File Version:\t%s\n", "0.8"); 
  else
    fprintf(fp,"# SMAP File Version:\t%s\n", "0.7"); 

  //take all this from output_xmap even though currently it only runs for pairsplit, which is only for pairwise, ie, no spots
  if(spots_filename[0]){/* Reference Aligner : refer to condensed version of reference CMAP */
#ifndef WIN32
    if(output_prefix[0] == '/')
      fprintf(fp,"# Reference Maps From:\t%s_r.cmap\n", output_prefix);
    else
      fprintf(fp,"# Reference Maps From:\t%s/%s_r.cmap\n", cwd, output_prefix);
#else
      fprintf(fp,"# Reference Maps From:\t%s_r.cmap\n", output_prefix);
#endif
    fprintf(fp,"# Query Maps From:\t%s",queryfilename);/* This is actually a modified cmap file */
    //if(DEBUG) assert(num_files==1); //this fails when using -if, but that should be handled correctly
    /*    for(register int i = 1; i < num_files;i++)
	  fprintf(fp,",%s",vfixx_filename[i]);*/
  }  else {
    if(mres >= 0.001 || mresSD > 0.0 || fabs(PixelLen - origPixelLen) > 1e-12){
      fprintf(fp,"# Reference Maps From:\t%s_r.cmap\n", absbasename);
      fprintf(fp,"# Query Maps From:\t%s_q.cmap", absbasename);
    } else {
#ifndef WIN32
      if(vfixx_filename[0][0]=='/')
	fprintf(fp,"# Reference Maps From:\t%s", vfixx_filename[0]);
      else
	fprintf(fp,",%s/%s",cwd,vfixx_filename[0]);
#else
	fprintf(fp,"# Reference Maps From:\t%s", vfixx_filename[0]);
#endif
      for(register int i = 1; i < RefIndex; i++){
#ifndef WIN32
	if(vfixx_filename[i][0]=='/')
	  fprintf(fp,",%s",vfixx_filename[i]);
	else
	  fprintf(fp,",%s/%s",cwd,vfixx_filename[i]);
#else
        fprintf(fp,",%s",vfixx_filename[i]);
#endif
      }
      fprintf(fp,"\n");
#ifndef WIN32
      if(vfixx_filename[QueryIndex][0] == '/')
	fprintf(fp,"# Query Maps From:\t%s",vfixx_filename[QueryIndex]);
      else
	fprintf(fp,"# Query Maps From:\t%s/%s", cwd, vfixx_filename[QueryIndex]);
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
  fprintf(fp,"# Xmap Entries From:\t%s\n",xfilename);

  //header (similar to xmap header)
  fprintf(fp,"#h SmapEntryID\tQryContigID\tRefcontigID1\tRefcontigID2\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tConfidence\tType");
  if( smap_xmapid_col )
    fprintf(fp,"\tXmapID1\tXmapID2");
  if( smap_linkid_col )
    fprintf(fp,"\tLinkID");
  if( smap_labelind_col )
    fprintf(fp,"\tQryStartIdx\tQryEndIdx\tRefStartIdx\tRefEndIdx");
  if( smap_zygosity )
    fprintf(fp,"\tZygosity\tGenotype\tGenotypeGroup");
  fprintf(fp,"\tRawConfidence\tRawConfidenceLeft\tRawConfidenceRight\tRawConfidenceCenter");

  if(SV_SIZE || SV_FREQ || SV_ORIENT)
    fprintf(fp,"\tSVsize");
  if(SV_FREQ || SV_ORIENT){
    fprintf(fp,"\tSVfreq");
    if(SV_FREQ>=2)
      fprintf(fp,"\tSVcov\tTotcov");
  }
  if(SV_ORIENT)
    fprintf(fp,"\torientation");
  fprintf(fp,"\n");
  fprintf(fp,"#f int        \tint        \tint         \tint         \tfloat      \tfloat    \tfloat      \tfloat     \tfloat \tstring ");
  if( smap_xmapid_col )
    fprintf(fp,"\tint\tint");
  if( smap_linkid_col )
    fprintf(fp,"\tint");
  if( smap_labelind_col )
    fprintf(fp,"\tint\tint\tint\tint");
  if( smap_zygosity )
    fprintf(fp,"\tstring\tint\tint");
  fprintf(fp,"\tfloat\tfloat\tfloat\tfloat"); //RawConfidence x4
  if(SV_SIZE || SV_FREQ || SV_ORIENT)
    fprintf(fp,"\tfloat");// -svSize
  if(SV_FREQ || SV_ORIENT){
    fprintf(fp,"\tfloat");// -svFreq
    if(SV_FREQ>=2)
      fprintf(fp,"\tfloat\tfloat");// -svFreq 2
  }
  if(SV_ORIENT)
    fprintf(fp,"\tstring");// -svTransOrient
  fprintf(fp,"\n");

  /* NEW88 : reset all xmapentries[i]->output to false then mark any matchgroups with 1-MG Indels with output = true */
  for(int i = 0; i < numxmapentries; i++)
    xmapentries[i]->output = false;

  for(int i = 0; i < numsv; i++){
    structuralVariation *p = &sv_array[i];
    if(p->algo_enum == sv_ps || p->indices[0] != p->indices[1])
      continue;

    if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
      continue;
    
    xmapentries[p->indices[0]] -> output = true;
  }

  // NOTE All SVs and 'end' svs are first traversed without writing to smap, just to mark all used xmapentries[]->output = true
  //  Now loop over all SVs and "write" them all to the smap (printing detailed info also) : The write is not actually done, just smap index (line number) is computed
  if(VERB){
    printf("Calling writeSVToSmap() with dowrite=false : numsv = %d\n",numsv);
    fflush(stdout);
  }

  int nsv2 = 1;
  for(int indp=0; indp < numsv; indp++) {
    structuralVariation* sv = &sv_array[indp];
    sv->checkBedGapOverlap(); //check for overlap with bed gap
    sv->maskTranslocation();
    sv->setTypeStr(); //set svtype and svtype_debug data members
    if(!LATE_ZYGOSITY)
      sv->setZygosity(); //must call for all entries before first final writeSVToSmap (ie, the one which actually writes)
    if(DEBUG>=2 && sv->algo_enum == sv_ps) assert(sv->type_enum != inversion_small_type);
    xmapEntry *entry1 = (sv->algo_enum == sv_ps) ? usexmapentries[sv->indices[0]] : xmapentries[sv->indices[0]];
    xmapEntry *entry2 = (sv->algo_enum == sv_ps) ? usexmapentries[sv->indices[1]] : xmapentries[sv->indices[1]];
    if( verbose && (sv->confidence < 0 || sv->confidence > smap_min_conf) ) { //conf filter
      printf("SV %d : sv_array[%d](xmapid= %d,%d) is %10s (%19s): qryid=%lld, %lld, refid=%lld, %lld, qrysiz=%10.4f, refsiz=%10.4f kb, gap=%i,conf=%.2f,dup=%d\n", 
	     nsv2, indp, entry1->xmapid, entry2->xmapid, sv->svtype, sv->svtype_debug, entry1->qrycontigid, entry2->qrycontigid, entry1->refcontigid, entry2->refcontigid, sv->querysize, sv->refsize, sv->gap_overlap, sv->confidence,sv->duplicate);
      fflush(stdout);

      if(DEBUG) assert(entry1->refcontigid == sv->refcontigid1);
      if(DEBUG) assert(entry2->refcontigid == sv->refcontigid2);
    }

    nsv2 += sv->writeSVToSmap(usexmapentries, nsv2, fp, false); //false for no write to fp
  }

  // all the indels are done--check to make sure no more 'end' svs at the begining of contigs that were skipped (because that contig has no indels)
  if(VERB>=2){
    printf("Calling checkEndSV() for contigs without indels : smapsize=%d,dowrite=false\n",smapsize);
    fflush(stdout);
  }
  for(int i = 0; i < smapsize; i++)  
    nsv2 += checkEndSV(fp, nsv2, usexmapentries[i]->qrycontigid, contigendalign, endalignsize, true, false, NULL, xmapindexMap);

  // do the 'end' svs at the starts of contigs with indels
  if(VERB>=2){
    printf("Calling checkEndSV() for the starts of contigs with indels: endalignsize=%d,dowrite=false\n",endalignsize);
    fflush(stdout);
  }
  for(int i=0; i < endalignsize; i++) 
    nsv2 += checkEndSV(fp, nsv2, contigendalign[i].contigid, contigendalign, endalignsize, true, false, NULL, xmapindexMap);

  // lastly, do the 'end' svs at the ends of contigs with indels
  if(VERB>=2){
    printf("Calling checkEndSV() for the ends of contigs with indels: endalignsize=%d,dowrite=false\n",endalignsize);
    fflush(stdout);
  }
  for(int i=0; i < endalignsize; i++)
    nsv2 += checkEndSV(fp, nsv2, contigendalign[i].contigid, contigendalign, endalignsize, false, false, NULL, xmapindexMap);

  if(DEBUG>=2){
    for(int i = 0; i < numsv - 1; i++){
      int k;
      if((k = compare_svs(&sv_array[i],&sv_array[i+1])) != -1){
	printf("compare_svs(&sv_array[%d],&sv_array[%d]) = %d : expected -1 (numsv=%d)\n",i,i+1,k,numsv);
	fflush(stdout);
	assert(k== -1);
      }
      if(!(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2))){
	printf("sv_array[%d,%d]:xmapid1=%d,%d,xmapid2=%d,%d (not ordered in ascending order)\n",i,i+1,sv_array[i].xmapid1,sv_array[i+1].xmapid1,sv_array[i].xmapid2,sv_array[i+1].xmapid2);
	fflush(stdout);
	assert(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2));
      }
    }

    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->output && entry2->output && entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum == sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	fflush(stdout);
	assert(entry1->output);
	assert(entry2->output);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
    }
  }

  double Xscale = 1.0;
  double Yscale = 1.0;
  /* NOTE if -ref or -mres -bpp were used, rescaled _r.cmap and _q.cmap files will be used, so no need to rescale the alignment locations */
  if(!(spots_filename[0] || (svcheck && numrefmaps > 0)) && fabs(PixelLen - origPixelLen) > 1e-12 /* NEW */){
    Yscale = Xscale = origPixelLen/PixelLen;
    if(PairSplit && PairsplitRef)/* -i inputs that correspond to reference is not scaled by bpp/500 */
      Yscale = 1.0;
  }

  /* mark all remaining entries in filtered usexmapentries[]->output as true (unless marked for removal by inversion_only=true), so they will be output even if not involved in any SV call */
  for(int i=0; i < smapsize; i++){
    if(VERB>=2){
      xmapEntry *p = usexmapentries[i];
      printf("usexmapentries[%d/%d]: xmapid=%d,qid=%lld,rid=%lld, inversion_only= %d, confidence= %0.2f, LogPvThreshold= %0.2f\n",
	     i,smapsize,p->xmapid,p->qrycontigid,p->refcontigid,p->inversion_only,p->confidence,LogPvThreshold);
      fflush(stdout);
    }
    if(!usexmapentries[i]->inversion_only /* NEW88 */ && usexmapentries[i]->confidence >= LogPvThreshold)
      usexmapentries[i]->output = true;
  }

  if(DEBUG>=2){/* verify that xmapentries[] are in ascending order of xmapid */
    for(int i = 0; i < numxmapentries - 1; i++){
      if(!(xmapentries[i]->xmapid < xmapentries[i+1]->xmapid)){
	printf("i=%d:xmapentries[i]->xmapid=%d,xmapentries[i+1]->xmapid=%d\n",i,xmapentries[i]->xmapid,xmapentries[i+1]->xmapid);
	fflush(stdout);
	assert(xmapentries[i]->xmapid < xmapentries[i+1]->xmapid);
      }
    }
  }

#if LATE_ZYGOSITY // NEW148 : Compute Zygosity of SVs based on MGs that will actually be output 

  //call fns in structuralVariation.h/cpp to check for alignment, sv overlaps for inversion zygosity calling
  doPairInversionBP(usexmapentries, sv_array, numsvIndel, numsv, false);
  if(VERB/* HERE >=2 */){
    printf("Completed doPairInversionBP: wall time= %0.6f secs\n", wtime());
    fflush(stdout);
  }
  alignmentOverlap(xmapentries, numxmapentries);
  alignmentOverlapInversion(xmapentries, numxmapentries);
  alignmentOverlapTranslocation(xmapentries, numxmapentries);
  svOverlap(xmapentries, numxmapentries);

  if(VERB){
    printf("Setting zygosity flag\n");
    fflush(stdout);
  }
  for(int indp=0; indp < numsv; indp++){
    structuralVariation* sv = &sv_array[indp];
    sv->setZygosity(); //must call for all entries before first final writeSVToSmap (ie, the one which actually writes)
  }

  if(VERB){
    printf("Finished sv zygosity, numsv=%d (wall time=%.6f)\n", numsv,wtime());
    fflush(stdout);
  }

#endif

  /* sort numxmapentries by refid, so corrected _r.cmap can be output */
  
  /* map each xmapentry to its output line number in the filtered xmap file (xmapidfilt) */
  int xmapcnt = 0;
  for(int i = 0; i < numxmapentries; i++){
    xmapEntry *xmap = xmapentries[i];
    xmapcnt += xmap->output ? 1 : 0;
    xmap->xmapidfilt = xmap->output ? xmapcnt : 0;
    if(VERB>=2){
      Calign *p = xmap->align;
      int rid = p->mapid1;
      int qid = p->mapid2;
      Cmap *Ymap = refmap[rid];
      Cmap *Xmap = Gmap[qid];
      printf("i=%d/%d: output=%d, xmapid=%d, xmapcnt= %d, xmapidfilt=%d: rid=%d(id=%lld),qid=%d(id=%lld)\n", 
	     i, numxmapentries, xmap->output ? 1 : 0, xmap->xmapid, xmapcnt, xmap->xmapidfilt, rid,Ymap->id,qid,Xmap->id);
      fflush(stdout);
    }
    if(DEBUG && xmapentries[i]->output) assert(xmapentries[i]->xmapidfilt > 0);
  }
  int *xmapidMap = new int[xmapidMax+1];/* map from original xmapid (location in _full.xmap) to xmapidfilt, which is 0 for entries that will not be output OR xmapid location in the new XMAP output */
  memset(xmapidMap, 0, xmapidMax * sizeof(int));
  for(int i = 0; i < numxmapentries; i++)
    xmapidMap[xmapentries[i]->xmapid] = xmapentries[i]->xmapidfilt;

  if(SV_FREQ){/* compute SV call frequency value */

    /* Update total reference BNX (weighted) coverage to reflect sum of BNX (weighted) coverage from consensus maps of filtered xmap entries */
    for(int n = 0; n < numrefmaps; n++){/* initialize refmap[n].sitecov[] = 0 */
      Cmap *p = refmap[n];
      int N = p->numsite[0];
      if(DEBUG) assert(p->sitecov[0]);
      for(int i = 0; i <= N+1; i++)
	p->sitecov[0][i] = 0.0;
    }

    for(int i = 0; i < numxmapentries; i++){
      xmapEntry *xmap = xmapentries[i];
      if(!xmap->output)
	continue;
    
      Calign *p = xmap->align;
      int rid = p->mapid1;
      int qid = p->mapid2;
      int U = p->numpairs; 
      Cmap *Ymap = refmap[rid];
      Cmap *Xmap = Gmap[qid];
      if(DEBUG) assert(Ymap->id == xmap->refcontigid);
      if(DEBUG) assert(Xmap->id == xmap->qrycontigid);
      //    int N = Ymap->numsite[0];
      int M = Xmap->numsite[0];
      FLOAT *Y = Ymap->site[0];
      FLOAT *X = Xmap->site[0];

      int LI = p->sites1[0];
      int LK = LI - p->sitesK1[0];
      int LJ = p->sites2[0];
      int RI = p->sites1[U-1];
      int RK = RI - p->sitesK1[U-1];
      int RJ = p->sites2[U-1];

      int Lij = p->Lend <= -2 ? LK : p->Lij1;
      int Rij = p->Rend <= -2 ? RI : p->Rij1;

      if(VERB>=2 && Ymap->id == 8 && max(Lij,23099) <= min(Rij,23169)){
	printf("xmapentries[i=%d]: xmapid=%d -> %d: rid=%d(id=%lld),qid=%d(id=%lld),Lij=%d,Rij=%d,U=%d,M=%d,N=%d:LI=%d,LK=%d,LJ=%d(cov=%0.3f),RI=%d,RK=%d,RJ=%d(cov=%0.3f)\n",
	       i,xmap->xmapid, xmap->xmapidfilt,rid,Ymap->id,qid,Xmap->id,Lij,Rij,U,M,Ymap->numsite[0],LI,LK,LJ,Xmap->sitecov[0][p->orientation ? M+1-LJ : LJ],   RI,RK,RJ,Xmap->sitecov[0][p->orientation ? M+1-RJ : RJ]);
	fflush(stdout);
      }

      if(LJ >= 1){ /* process left end */
	int J = LJ - 1;
	double xmax = p->orientation ? X[M+1-J] - X[M+1-LJ] : X[LJ] - X[J];
	for(int I = LK; --I >= Lij; ){
	  double y = Y[LK] - Y[I];
	  while(y > xmax){
	    if(--J < 0)
	      break;
	    xmax = p->orientation ? X[M+1-J] - X[M+1-LJ] : X[LJ] - X[J];
	  }
	  if(y > xmax || J < 0)
	    break;
	  double covX;
	  Ymap->sitecov[0][I] += covX = Xmap->sitecov[0]==NULL ? 1.0f : Xmap->sitecov[0][p->orientation ? M+(FRAGCOV_FIX?0:1)-J : J];
	  if(VERB>=2 && Ymap->id == 8 && max(Lij,23099) <= min(Rij,23169)){
	    printf("\t i=%d(Lij=%d,LK=%d,LI=%d): Ymap->sitecov[0][i] += %0.3f -> %0.3f (LJ= %d,J=%d)\n",i,Lij,LK,LI,covX,Ymap->sitecov[0][i],LJ,J);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(isfinite(Ymap->sitecov[0][I]));
	}
      }

      /* handle label (LI,LK,LJ) */
      float covL = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? Xmap->sitecov[0][p->orientation ? M-LJ : LJ] : 
	min(Xmap->sitecov[0][p->orientation ? M+1-LJ : LJ], Xmap->sitecov[0][p->orientation ? M-LJ : LJ+1]);
      for(int i = LK; i < LI; i++){
	Ymap->sitecov[0][i] += covL;
	if(VERB>=2 && Ymap->id == 8 && max(Lij,23099) <= min(Rij,23169)){
	  printf("\t i=%d(LK=%d,LI=%d,LJ=%d): Ymap->sitecov[0][i] += %0.3f -> %0.3f (LJ= %d)\n",i,LK,LI,LJ,covL,Ymap->sitecov[0][i],LJ);
	  fflush(stdout);
	}
	if(DEBUG) assert(isfinite(Ymap->sitecov[0][i]));
      }

      for(int u = 1; u < U; u++){
	int RI = p->sites1[u];
	int RK = RI - p->sitesK1[u];
	int RJ = p->sites2[u];

	/* interpolate between (LI,LK,LJ) and (RI,RK,RJ) */
	float cov = 0.0, wt = 0.0;
	for(int J = LJ; J < RJ; J++){
	  cov += Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? Xmap->sitecov[0][p->orientation ? M-J : J] : min(Xmap->sitecov[0][p->orientation ? M+1-J : J], Xmap->sitecov[0][p->orientation ? M-J : J+1]);
	  wt += 1.0;
	}
	if(DEBUG) assert(wt > 0.0);
	cov /= wt;
	
	float covR = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? Xmap->sitecov[0][p->orientation ? M-RJ : RJ] :
	  min(Xmap->sitecov[0][p->orientation ? M+1-RJ : RJ], Xmap->sitecov[0][p->orientation ? M-RJ : RJ+1]);

	for(int i = LI; i < RK; i++){
	  Ymap->sitecov[0][i] += cov;
	  if(VERB>=2 && Ymap->id == 8 && max(Lij,23099) <= min(Rij,23169)){
	    printf("\t i=%d(LI=%d,RK=%d,LJ=%d,RJ=%d): Ymap->sitecov[0][i] += %0.3f -> %0.3f (covR= %0.3f,covL= %0.3f)\n",i,LI,RK,LJ,RJ,cov,Ymap->sitecov[0][i], covR,covL);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(isfinite(Ymap->sitecov[0][i]));
	}

	/* handle label (RI,RK,RJ) */
	for(int i = RK; i < RI; i++){
	  Ymap->sitecov[0][i] += covR;
	  if(VERB>=2 && Ymap->id == 8 && max(Lij,23099) <= min(Rij,23169)){
	    printf("\t i=%d(RK=%d,RI=%d,RJ=%d): Ymap->sitecov[0][i] += %0.3f -> %0.3f\n",i,RK,RI,RJ,covR,Ymap->sitecov[0][i]);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(isfinite(Ymap->sitecov[0][i]));
	}

	LI = RI;
	LK = RK;
	LJ = RJ;
	covL = covR;
      }

      if(LJ <= M){/* process right end */
	int J = LJ;
	double xmax = p->orientation ? X[M+1-LJ] - X[M-J] : X[J+1] - X[LJ];
	for(int I = LI; ++I <= Rij; ){
	  double y = Y[I] - Y[LI];
	  while(y > xmax){
	    if(++J > M)
	      break;
	    xmax = p->orientation ? X[M+1-LJ] - X[M-J] : X[J+1] - X[LJ];
	  }
	  if(y > xmax || J > M)
	    break;

	  Ymap->sitecov[0][I] += Xmap->sitecov[0]==NULL ? 1.0f : Xmap->sitecov[0][p->orientation ? M+(FRAGCOV_FIX?0:1)-J : J];
	  if(VERB>=2 && Ymap->id == 8 && max(Lij,23099) <= min(Rij,23169)){
	    printf("\t I=%d(LI=%d,Rij=%d): Ymap->sitecov[0][I] += %0.3f -> %0.3f (LJ= %d,M=%d)\n",I,LI,Rij,covL,Ymap->sitecov[0][I],LJ,M);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(isfinite(Ymap->sitecov[0][I]));
	}
      }
    }

    /* NEW175: Update reference coverage to include gaps between MGs for 2-MG SVs (same reference contig), provided no 3rd MG spans part of the qry gap */
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf) || p->type_enum == complex_type){
	p->frequency = -1.0;
	continue;
      }

      if(p->algo_enum != sv_ps)
	continue;/* 1-MG indels and small Duplications & Inversions don't have regions not covered by any MG */

      if(p->refcontigid1 != p->refcontigid2 || p->type_enum == intrachr_translocation)
	continue;/* This should only be translocations, which use only the 2 ends of the breakpoint interval */

      int k1 = p->indices[0];
      int k2 = p->indices[1];
      xmapEntry *entry1 = usexmapentries[k1];
      xmapEntry *entry2 = usexmapentries[k2];
      if(DEBUG) assert(entry1->qrycontigid == entry2->qrycontigid);
      if(DEBUG) assert(entry1->refcontigid == p->refcontigid1 && entry2->refcontigid == p->refcontigid2);
      Calign *p1 = entry1->align;
      Calign *p2 = entry2->align;
      if(DEBUG) assert(p2->mapid2 == p1->mapid2);
      int qid = p1->mapid2;
      
      double minq1 = min(entry1->qrystartpos, entry1->qryendpos);
      double maxq1 = max(entry1->qrystartpos, entry1->qryendpos);
      double minq2 = min(entry2->qrystartpos, entry2->qryendpos);
      double maxq2 = max(entry2->qrystartpos, entry2->qryendpos);

      if(overlapLength(minq1,maxq1,minq2,maxq2) < 0.0){/* check if qrygap is overlapped by 3rd MG */
	double qrystart = min(maxq1,maxq2);
	double qryend = max(minq1,minq2);
	if(DEBUG>=2) assert(qryend > qrystart);

	// locate range of usexmapentries with same qrycontigid 
	long long qrycontigid = entry1->qrycontigid;
	int kmin = min(k1,k2), kmax = max(k1,k2);
	for(;--kmin >= 0;){
	  xmapEntry *pMGk = usexmapentries[kmin];
	  if(pMGk->qrycontigid != qrycontigid){
	    kmin++;
	    break;
	  }
	}
	kmin = max(0,kmin);
      
	for(;++kmax < smapsize;){
	  xmapEntry *pMGk = usexmapentries[kmax];
	  if(pMGk->qrycontigid != qrycontigid){
	    kmax--;
	    break;
	  }
	}
	kmax = min(smapsize-1, kmax);
	if(DEBUG>=2){
	  for(int k = 0; k < kmin; k++)
	    assert(usexmapentries[k]->qrycontigid < qrycontigid);
	  for(int k = kmax; ++k < smapsize;)
	    assert(usexmapentries[k]->qrycontigid > qrycontigid);
	}

	if(kmin == min(k1,k2))
	  kmin++;
	if(kmax == max(k1,k2))
	  kmax--;
      
	int k = kmin;
	for(;k <= kmax; k++){
	  if(k == k1 || k == k2)
	    continue;
	  xmapEntry *pMGk = usexmapentries[k];
	  if(DEBUG>=2) assert(pMGk->qrycontigid == qrycontigid);
	  
	  double minq3 = min(pMGk->qrystartpos, pMGk->qryendpos);
	  double maxq3 = max(pMGk->qrystartpos, pMGk->qryendpos);

	  if(overlapLength(qrystart,qryend, minq3,maxq3) > 0.0)
	    break;// suppress this SV
	}
	if(k <= kmax)
	  continue;// suppress this SV : do not add qrygap coverge to refgap
      }

      if(DEBUG) assert(p2->mapid1 == p1->mapid1);
      int rid = p1->mapid1;
      Cmap *Ymap = refmap[rid];
      Cmap *Xmap = Gmap[qid];
      int N = Ymap->numsite[0];
      int M = Xmap->numsite[0];
      if(DEBUG) assert(Ymap->id == p->refcontigid1);
      if(DEBUG) assert(Ymap->sitecov[0]);

      int LI = p->refstartidx;
      int RI = p->refstopidx;
      int LJ = p->querystartidx;
      int RJ = p->querystopidx;
      if(DEBUG>=1+RELEASE && !(1 <= LI && LI <= RI && RI <= N)){
	printf("sv_array[i=%d/%d]:querystartidx=%d,querystopidx=%d,refstartidx=%d,refstopidx=%d, refcontigid=%lld,%lld(rid=%d),qid=%d(id=%lld),M=%d,N=%d\n",
	       i,numsv,p->querystartidx,p->querystopidx,p->refstartidx,p->refstopidx,p->refcontigid1,p->refcontigid2,rid,qid,Xmap->id,Xmap->numsite[0],Ymap->numsite[0]);
	printf("\t type_num=%d,svtype=%s, querystart= %0.4f, querystop= %0.4f, refstart= %0.4f, refstop= %0.4f, xmapid1= %d -> %d, xmapid2= %d -> %d\n",
	       p->type_enum,p->svtype,p->querystart * 1e-3, p->querystop * 1e-3, p->refstart * 1e-3, p->refstop * 1e-3,
	       entry1->xmapid,entry1->xmapidfilt,entry2->xmapid,entry2->xmapidfilt); 
	printf("\t LI=%d,RI=%d,N=%d\n",LI,RI,N);
	fflush(stdout);
	assert(1 <= LI && LI <= RI && RI <= N);
      }

      if(DEBUG && !(0 <= LJ && 0 <= RJ)){
	printf("sv_array[i=%d/%d]:querystartidx=%d,querystopidx=%d,refstartidx=%d,refstopidx=%d, refcontigid=%lld,%lld(rid=%d),qid=%d(id=%lld),M=%d,N=%d\n",
	       i,numsv,p->querystartidx,p->querystopidx,p->refstartidx,p->refstopidx,p->refcontigid1,p->refcontigid2,rid,qid,Xmap->id,Xmap->numsite[0],Ymap->numsite[0]);
	printf("\t type_num=%d,svtype=%s, querystart= %0.4f, querystop= %0.4f, refstart= %0.4f, refstop= %0.4f, xmapid1= %d -> %d, xmapid2= %d -> %d\n",
	       p->type_enum,p->svtype,p->querystart * 1e-3, p->querystop * 1e-3, p->refstart * 1e-3, p->refstop * 1e-3,
	       entry1->xmapid,entry1->xmapidfilt,entry2->xmapid,entry2->xmapidfilt); 
	
	fflush(stdout);
	assert(0 <= p->querystartidx && 0 <= p->querystopidx);
      }
      if(DEBUG) assert(LJ <= M+1);
      if(DEBUG) assert(RJ <= M+1);
      LI = min(N,max(1,LI));
      RI = min(N,max(1,RI));
      LJ = min(M,max(1,LJ));
      RJ = min(M,max(1,RJ));
      
      if(RJ < LJ){/* swap LJ,RJ */
	int tmp = RJ;
	RJ = LJ;
	LJ = tmp;
      }

      if(VERB>=2 && Ymap->id == 8 && max(LI,23099) <= min(RI,23169)){
	printf("sv_array[i=%d/%d]:querystartidx=%d,querystopidx=%d,refstartidx=%d,refstopidx=%d, rid=%d(id=%lld),qid=%d(id=%lld),M=%d,N=%d\n",
	       i,numsv,p->querystartidx,p->querystopidx,p->refstartidx,p->refstopidx,rid,p->refcontigid1,qid,Xmap->id,Xmap->numsite[0],Ymap->numsite[0]);
	printf("\t type_num=%d,svtype=%s, querystart= %0.4f, querystop= %0.4f, refstart= %0.4f, refstop= %0.4f\n",
	       p->type_enum,p->svtype,p->querystart * 1e-3, p->querystop * 1e-3, p->refstart * 1e-3, p->refstop * 1e-3); 
	printf("\t LI=%d,RI=%d,LJ=%d,RJ=%d:\n",
	       LI,RI,LJ,RJ);
	fflush(stdout);
      }
      
      /* interpolate between (LI,LJ) and (RI,RJ) */
      float cov = 0.0, wt = 0.0;
      if(RJ <= LJ){
	cov = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? min(Xmap->sitecov[0][LJ], Xmap->sitecov[0][max(0,LJ-1)]) : 
	  min(Xmap->sitecov[0][LJ], min(Xmap->sitecov[0][max(0,LJ-1)], Xmap->sitecov[0][min(M,LJ+1)]));
	if(DEBUG) assert(isfinite(cov));
	wt = 1.0;
      } else {
	for(int J = LJ; J <= RJ - (FRAGCOV_FIX ? 1 : 0); J++){
	  float xcovJ = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? Xmap->sitecov[0][J] : min(Xmap->sitecov[0][J], Xmap->sitecov[0][J+1]);
	  if(VERB>=2 && Ymap->id == 8 && max(LI,23099) <= min(RI,23169)){
	    printf("\t J=%d:xcovJ= %0.3f (sitecov[J]=%0.3f,sitecov[J+1]=%0.3f)\n",J,xcovJ,Xmap->sitecov[0][J],Xmap->sitecov[0][J+1]);
	    fflush(stdout);
	  }
	  if(xcovJ > 0.0){
	    cov += xcovJ;
	    wt += 1.0;
	  }
	  if(DEBUG) assert(isfinite(cov));
	}
      }
      if(wt > 0.0)
	cov /= wt;

      if(VERB>=3 && Ymap->id == 8){
	printf("\t LI=%d,RI=%d:cov = %0.3f\n",LI,RI,cov);
	fflush(stdout);
      }

      for(int I = LI; I <= max(LI,RI - (FRAGCOV_FIX?1:0)); I++){
	Ymap->sitecov[0][I] += cov;
	if((VERB>=2 && Ymap->id == 8 && max(LI,23099) <= min(RI,23169)) || (DEBUG && !isfinite(Ymap->sitecov[0][I]))){
	  float covL = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? Xmap->sitecov[0][p->orientation ? M-LJ : LJ] : 
	    min(Xmap->sitecov[0][p->orientation ? M+1-LJ : LJ], Xmap->sitecov[0][p->orientation ? M-LJ : LJ+1]);
	  float covR = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? Xmap->sitecov[0][p->orientation ? M-RJ : RJ] :
	    min(Xmap->sitecov[0][p->orientation ? M+1-RJ : RJ], Xmap->sitecov[0][p->orientation ? M-RJ : RJ+1]);
	  printf("\t I=%d(LI=%d,RI=%d): Ymap->sitecov[0][I] += %0.3f -> %0.3f (RJ= %d, covR= %0.3f,covL= %0.3f)\n",I,LI,RI,cov,Ymap->sitecov[0][I], RJ,covR,covL);
	  fflush(stdout);
	}
	if(DEBUG && !isfinite(Ymap->sitecov[0][I])){
	  printf("Ymap->sitecov[0][I]= %0.6e, I=%d\n",Ymap->sitecov[0][I],I);
	  fflush(stdout);
	  assert(isfinite(Ymap->sitecov[0][I]));
	}
      }
    }

    /* compute SV call frequency based on ratio of weighted BNX coverage (SV call average query contig BNX coverage divided by reference contig BNX coverage) */
    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf) || p->type_enum == complex_type){
	p->frequency = -1.0;
	continue;
      }

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(DEBUG) assert(entry1->qrycontigid == entry2->qrycontigid);
      if(DEBUG) assert(entry1->refcontigid == p->refcontigid1 && entry2->refcontigid == p->refcontigid2);
      Calign *p1 = entry1->align;
      Calign *p2 = entry2->align;
      int rid1 = p1->mapid1;
      int rid2 = p2->mapid1;
      if(DEBUG) assert(p2->mapid2 == p1->mapid2);
      int qid = p1->mapid2;
      Cmap *Ymap1 = refmap[rid1];
      Cmap *Ymap2 = refmap[rid2];
      Cmap *Xmap = Gmap[qid];
      int N1 = Ymap1->numsite[0];
      int N2 = Ymap2->numsite[0];
      int M = Xmap->numsite[0];
      if(DEBUG) assert(Ymap1->id == p->refcontigid1);
      if(DEBUG) assert(Ymap2->id == p->refcontigid2);
      if(DEBUG) assert(Ymap1->sitecov[0] && Ymap2->sitecov[0]);
      int LI = p->refstartidx;
      int RI = p->refstopidx;
      int LJ = p->querystartidx;
      int RJ = p->querystopidx;
      if(DEBUG>=1+RELEASE) assert(1 <= LI && 1 <= RI);
      if(DEBUG>=1+RELEASE) assert(LI <= N1 && RI <= N2);
      if(DEBUG>=1+RELEASE && p->refcontigid1 == p->refcontigid2 && p->type_enum != intrachr_translocation) assert(1 <= LI && LI <= RI && RI <= N1 && N1==N2);

      if(DEBUG && !(0 <= LJ && 0 <= RJ)){
	printf("sv_array[i=%d/%d]:querystartidx=%d,querystopidx=%d,refstartidx=%d,refstopidx=%d, refcontigid=%lld,%lld(rid=%d,%d),qid=%d(id=%lld),M=%d,N=%d,%d\n",
	       i,numsv,p->querystartidx,p->querystopidx,p->refstartidx,p->refstopidx,p->refcontigid1,p->refcontigid2,rid1,rid2,qid,Xmap->id,Xmap->numsite[0],Ymap1->numsite[0],Ymap2->numsite[0]);
	printf("\t type_num=%d,svtype=%s, querystart= %0.4f, querystop= %0.4f, refstart= %0.4f, refstop= %0.4f, xmapid1= %d -> %d, xmapid2= %d -> %d\n",
	       p->type_enum,p->svtype,p->querystart * 1e-3, p->querystop * 1e-3, p->refstart * 1e-3, p->refstop * 1e-3,
	       entry1->xmapid,entry1->xmapidfilt,entry2->xmapid,entry2->xmapidfilt); 
	
	fflush(stdout);
	assert(0 <= p->querystartidx && 0 <= p->querystopidx);
      }
      if(DEBUG) assert(LJ <= M+1);
      if(DEBUG) assert(RJ <= M+1);
      LI = min(N1,max(1,LI));
      RI = min(N2,max(1,RI));
      LJ = min(M,max(1,LJ));
      RJ = min(M,max(1,RJ));

      if(RJ < LJ){/* swap LJ,RJ */
	int tmp = RJ;
	RJ = LJ;
	LJ = tmp;
      }

      if(VERB>=2 && (Xmap->id == 11824 || Xmap->id==351)){
	printf("sv_array[i=%d/%d]:querystartidx=%d,querystopidx=%d,refstartidx=%d,refstopidx=%d, rid=%d,%d(id=%lld,%lld),qid=%d(id=%lld),M=%d,N=%d,%d\n",
	       i,numsv,p->querystartidx,p->querystopidx,p->refstartidx,p->refstopidx,rid1,rid2,p->refcontigid1,p->refcontigid2,qid,Xmap->id,Xmap->numsite[0],Ymap1->numsite[0],Ymap2->numsite[0]);
	printf("\t type_num=%d,svtype=%s, querystart= %0.4f, querystop= %0.4f, refstart= %0.4f, refstop= %0.4f\n",
	       p->type_enum,p->svtype,p->querystart * 1e-3, p->querystop * 1e-3, p->refstart * 1e-3, p->refstop * 1e-3); 
	printf("\t LI=%d,RI=%d,LJ=%d,RJ=%d:\n",
	       LI,RI,LJ,RJ);
	fflush(stdout);
      }

      float refcov = 0.0;
      if(p->refcontigid1 == p->refcontigid2 && p->type_enum != intrachr_translocation){
	if(DEBUG) assert(Ymap1 == Ymap2);
	refcov = 0.0;
	float wt = 0.0;
	for(int I = LI; I <= max(LI, RI - (FRAGCOV_FIX?1:0)); I++){
	  if(Ymap1->sitecov[0][I] > 0.0){
	    if(VERB>=2 && (Xmap->id == 11824 || Xmap->id==351)){
	      printf("\t I=%d:ycovI= %0.3f\n",I,Ymap1->sitecov[0][I]);
	      fflush(stdout);
	    }
	    refcov += Ymap1->sitecov[0][I];
	    wt += 1.0;
	  }
	}
	if(wt > 0.0)
	  refcov /= wt;
      } else { /* just average the ends of the breakpoint interval */
	if(FRAGCOV_FIX)
	  refcov = 0.5 * (max(Ymap1->sitecov[0][LI],Ymap1->sitecov[0][LI-1]) + max(Ymap2->sitecov[0][RI],Ymap2->sitecov[0][RI-1]));
	else
	  refcov = 0.5 * (Ymap1->sitecov[0][LI] + Ymap2->sitecov[0][RI]);
      }

      float qrycov = 0.0, wt = 0.0;
      if(p->refcontigid1 == p->refcontigid2 && p->type_enum != intrachr_translocation){
	if(RJ <= LJ){
	  qrycov = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? min(Xmap->sitecov[0][LJ], Xmap->sitecov[0][max(0,LJ-1)]) : 
	    min(Xmap->sitecov[0][LJ], min(Xmap->sitecov[0][max(0,LJ-1)], Xmap->sitecov[0][min(M,LJ+1)]));
	  wt = 1.0;
	} else {
	  for(int J = LJ; J <= RJ - (FRAGCOV_FIX ? 1 : 0); J++){
	    float xcovJ = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? Xmap->sitecov[0][J] : min(Xmap->sitecov[0][J], Xmap->sitecov[0][J+1]);
	    if(VERB>=2 && (Xmap->id == 11824 || Xmap->id==351)){
	      printf("\t J=%d:xcovJ= %0.3f (sitecov[J]=%0.3f,sitecov[J+1]=%0.3f)\n",J,xcovJ,Xmap->sitecov[0][J],Xmap->sitecov[0][J+1]);
	      fflush(stdout);
	    }
	    if(xcovJ > 0.0){
	      qrycov += xcovJ;
	      wt += 1.0;
	    }
	  }
	}
      } else { /* just average the ends of the breakpoint interval */
	if(RJ <= LJ){
	  int J = LJ;
	  float xcovJ = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? min(Xmap->sitecov[0][J-1],Xmap->sitecov[0][J]) : Xmap->sitecov[0][J];
	  if(xcovJ > 0.0){
	    qrycov += xcovJ;
	    wt += 1.0;
	  }
	} else {
	  int J = LJ;
	  float xcovJ = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? min(Xmap->sitecov[0][J-1],Xmap->sitecov[0][J]) : Xmap->sitecov[0][J];
	  if(xcovJ > 0.0){
	    qrycov += xcovJ;
	    wt += 1.0;
	  }
	  
	  J = RJ - (FRAGCOV_FIX ? 1 : 0);
	  xcovJ = Xmap->sitecov[0]==NULL ? 1.0f : FRAGCOV_FIX ? min(Xmap->sitecov[0][J-1],Xmap->sitecov[0][J]) : Xmap->sitecov[0][J];
	  if(xcovJ > 0.0){
	    qrycov += xcovJ;
	    wt += 1.0;
	  }
	}
      }
      if(wt > 0.0)
	qrycov /= wt;
    
      if((DEBUG>=1+RELEASE && !(qrycov >= 0.0 && refcov > 0.0 && refcov >= qrycov * 0.9)) || (VERB>=2 && (Xmap->id == 11824 || Xmap->id == 351))){
	printf("WARNING:sv_array[i=%d/%d]:querystartidx=%d,querystopidx=%d,refstartidx=%d,refstopidx=%d, rid=%d,%d(id=%lld,%lld),qid=%d(id=%lld),M=%d,N=%d,%d\n",
	       i,numsv,p->querystartidx,p->querystopidx,p->refstartidx,p->refstopidx,rid1,rid2,p->refcontigid1,p->refcontigid2,qid,Xmap->id,Xmap->numsite[0],Ymap1->numsite[0],Ymap2->numsite[0]);
	printf("\t type_num=%d,svtype=%s, querystart= %0.4f, querystop= %0.4f, refstart= %0.4f, refstop= %0.4f\n",
	       p->type_enum,p->svtype,p->querystart * 1e-3, p->querystop * 1e-3, p->refstart * 1e-3, p->refstop * 1e-3); 
	printf("\t LI=%d,RI=%d,LJ=%d,RJ=%d: refcov= %0.3f, qrycov= %0.3f,xmapid1= %d -> %d, xmapid2= %d -> %d\n",
	       LI,RI,LJ,RJ,refcov,qrycov,	       entry1->xmapid,entry1->xmapidfilt,entry2->xmapid,entry2->xmapidfilt); 
	printf("\t Ymap1->sitecov[0][LI]= %0.3f, Ymap2->sitecov[0][RI]= %0.3f\n", Ymap1->sitecov[0][LI], Ymap2->sitecov[0][RI]);
	fflush(stdout);

	//  HERE HERE	if(DEBUG>=1+RELEASE) assert(qrycov >= 0.0 && refcov > 0.0 && refcov >= qrycov * 0.9);
      }

      p->frequency = (refcov > 0.0f) ? min(1.0f, qrycov / refcov) : -1.0;
      p->qrycov = qrycov;
      p->refcov = max(qrycov,refcov);
    }
  }// if(SV_FREQ)

  if(DEBUG>=2){
    for(int i = 0; i < numsv - 1; i++){
      int k;
      if((k = compare_svs(&sv_array[i],&sv_array[i+1])) != -1){
	printf("compare_svs(&sv_array[%d],&sv_array[%d]) = %d : expected -1 (numsv=%d)\n",i,i+1,k,numsv);
	fflush(stdout);
	assert(k== -1);
      }
      if(!(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2))){
	printf("sv_array[%d,%d]:xmapid1=%d,%d,xmapid2=%d,%d (not ordered in ascending order)\n",i,i+1,sv_array[i].xmapid1,sv_array[i+1].xmapid1,sv_array[i].xmapid2,sv_array[i+1].xmapid2);
	fflush(stdout);
	assert(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2));
      }

      /* make sure xmapid is still ordered when replaced by the new filtered (compacted) xmapid values, but only for those SV's that will be output */
      if(sv_array[i].duplicate || sv_array[i].confidence == 0 || (sv_array[i].confidence > 0 && sv_array[i].confidence < smap_min_conf))
	continue;
      if(sv_array[i+1].duplicate || sv_array[i+1].confidence == 0 || (sv_array[i+1].confidence > 0 && sv_array[i+1].confidence < smap_min_conf))
	continue;

      if(!(xmapidMap[sv_array[i].xmapid1] <= xmapidMap[sv_array[i+1].xmapid1] && (xmapidMap[sv_array[i].xmapid1] < xmapidMap[sv_array[i+1].xmapid1] || xmapidMap[sv_array[i].xmapid2] <= xmapidMap[sv_array[i+1].xmapid2]))){
	printf("sv_array[%d,%d]:xmapid1=%d,%d,xmapid2=%d,%d, xmapid1filt=%d,%d,xmapid2filt=%d,%d, output1=%d,%d output2=%d,%d(not ordered in ascending order)\n",i,i+1,
	       sv_array[i].xmapid1,sv_array[i+1].xmapid1,sv_array[i].xmapid2,sv_array[i+1].xmapid2,
	       xmapidMap[sv_array[i].xmapid1],xmapidMap[sv_array[i+1].xmapid1],xmapidMap[sv_array[i].xmapid2],xmapidMap[sv_array[i+1].xmapid2],
	       usexmapentries[sv_array[i].indices[0]]->output,usexmapentries[sv_array[i+1].indices[0]]->output,usexmapentries[sv_array[i].indices[1]]->output,usexmapentries[sv_array[i+1].indices[1]]->output);
	fflush(stdout);

	if(xmapidMap[sv_array[i].xmapid1] > 0 && xmapidMap[sv_array[i+1].xmapid1] > 0 && xmapidMap[sv_array[i].xmapid2] > 0 && xmapidMap[sv_array[i+1].xmapid2] > 0)
	  assert(xmapidMap[sv_array[i].xmapid1] <= xmapidMap[sv_array[i+1].xmapid1] && (xmapidMap[sv_array[i].xmapid1] < xmapidMap[sv_array[i+1].xmapid1] || xmapidMap[sv_array[i].xmapid2] <= xmapidMap[sv_array[i+1].xmapid2]));
      }
    }

    for(int i = 0; i < numsv; i++){
      structuralVariation *p = &sv_array[i];
      if(p->duplicate || p->confidence == 0 || (p->confidence > 0 && p->confidence < smap_min_conf))
	continue;

      xmapEntry *entry1 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[0]];
      xmapEntry *entry2 = ((p->algo_enum == sv_ps) ? usexmapentries : xmapentries)[p->indices[1]];
      if(!(entry1->output && entry2->output && entry1->xmapid == p->xmapid1 && entry2->xmapid == p->xmapid2)){
	printf("sv_array[%d]:duplicate=%d,confidence=%0.2f,xmapid1=%d,xmapid2=%d, indices[0]=%d,indices[1]=%d,sv_ps=%d,entry1->output=%d,xmapid=%d; entry2->output=%d,xmapid=%d\n",
	       i,p->duplicate?1:0,p->confidence,p->xmapid1,p->xmapid2,p->indices[0],p->indices[1],p->algo_enum==sv_ps ? 1 : 0,entry1->output,entry1->xmapid,entry2->output,entry2->xmapid);
	fflush(stdout);
	assert(entry1->output);
	assert(entry2->output);
	assert(entry1->xmapid == p->xmapid1);
	assert(entry2->xmapid == p->xmapid2);
      }
    }
  }

  Calign **Lalignment = new Calign*[xmapcnt];
  size_t aligncnt = 0;
  for(int i = 0; i < numxmapentries; i++){
    if(xmapentries[i]->output){
      Lalignment[aligncnt++] = xmapentries[i]->align;
      if(DEBUG>=1+RELEASE && !(aligncnt == (size_t)xmapentries[i]->xmapidfilt)){
	printf("WARNING: i=%d/%d:xmapentries[i]:xmapid=%d,xmapidfilt=%d,aligncnt=%lu\n",i,numxmapentries,xmapentries[i]->xmapid,xmapentries[i]->xmapidfilt,aligncnt);
	fflush(stdout);
	assert(aligncnt == (size_t)xmapentries[i]->xmapidfilt);
      }
    }
  }

  /* rename original XMAP and replace it with filtered version */
  char xfilenamefull[PATH_MAX]; // unfiltered xmap name
  sprintf(xfilenamefull,"%s_full.xmap",basename);

  printf("Renaming %s to %s\n", xfilename, xfilenamefull);
  rename(xfilename,xfilenamefull);

  output_xmapIO(xfilename, Lalignment, (size_t)0, aligncnt, basename, queryfilename, 0, NULL, Xscale, Yscale, false);

  /* write condensed reference (for filtered xmap) */

  /* sort Lalignment[] array : this does NOT affect .smap or .xmap output since Lalignment[] array is not used any more and .xmap output already is completed */
  /* sort alignments in order of mapid1, sites1[0],sites1[N+1], id2, sites2[0],sites2[M+1] (null alignments last) */
  qsort(Lalignment, aligncnt, sizeof(Calign *), (intcmp*)RefSiteInc);

  // NOTE : assumes reference maps YYmap[0..numrefmaps-1] are in ascending order of id value

  /* compute numalign_start[0..numrefmaps-1], numalign_end[0..numrefmaps-1] */
  size_t *numalign_start = new size_t[numrefmaps*2];
  size_t *numalign_end = &numalign_start[numrefmaps];

  int refstart = (aligncnt <= 0) ? 0 : Lalignment[0]->mapid1;
  long long pid = YYmap[refstart]->id;
  int refid = refstart;
  for(int i = 0; i < refid; i++)
    numalign_start[i] = numalign_end[i] = 0;
  numalign_start[refid] = 0;

  for(size_t i = 1; i < aligncnt; i++){
    Calign *p = Lalignment[i];
    if(DEBUG && !(refid <= p->mapid1 && p->mapid1 < numrefmaps)){
      printf("i=%lu,aligncnt=%lu:refid=%d,pid=%lld,mapid1=%d,numrefmaps=%d\n",
	     i,aligncnt,refid,pid,p->mapid1,numrefmaps);
      for(int I = 0; I < numrefmaps; I++)
	printf("YYmap[%d]: id=%lld\n",I,YYmap[I]->id);
      fflush(stdout);
      assert(refid <= p->mapid1 && p->mapid1 < numrefmaps);
    }

    long long nid = YYmap[p->mapid1]->id;

    if(nid != pid){
      numalign_end[refid] = i;
      if(VERB>=2){
	printf("refid=%d:numalign_start[refid]=%lu,numalign_end[refid]=%lu (i=%lu,p->mapid1=%d,nid=%lld,pid=%lld)\n",
	       refid,numalign_start[refid],numalign_end[refid],i,p->mapid1,nid,pid);
	fflush(stdout);
      }

      while(++refid < p->mapid1)
	numalign_start[refid] = numalign_end[refid] = i;

      numalign_start[refid = p->mapid1] = i;
      pid = nid;
    }
    if(DEBUG && !(p->mapid1 == refid)){
      printf("i=%lu/%lu: pid=%lld,nid=%lld,p->mapid1=%d,refid=%d\n",
	     i,aligncnt,pid,nid,p->mapid1,refid);
      fflush(stdout);
      assert(p->mapid1 == refid);
    }
  }
  numalign_end[refid] = aligncnt;
  if(VERB>=2){
    printf("refid=%d:numalign_start[refid]=%lu,numalign_end[refid]=%lu\n",  refid,numalign_start[refid],numalign_end[refid]);
    fflush(stdout);
  }
  if(DEBUG) assert(refid < numrefmaps);
  for(int i = refid+1; i < numrefmaps; i++)
    numalign_start[i] = numalign_end[i] = aligncnt;

  /* output condensed reference (for filtered xmap) */
  int split_maps = (MappedUnsplit >= 0) ? !MappedUnsplit : (SplitRef && CMapID > 0) ? 1 : 0;
  output_csites(0, numrefmaps-1, output_prefix, split_maps, Lalignment, numalign_start, numalign_end, false);

  delete [] Lalignment;
  delete [] numalign_start;


  if(DEBUG>=2){
    for(int i = 0; i < numsv - 1; i++){
      int k;
      if((k = compare_svs(&sv_array[i],&sv_array[i+1])) != -1){
	printf("compare_svs(&sv_array[%d],&sv_array[%d]) = %d : expected -1 (numsv=%d)\n",i,i+1,k,numsv);
	fflush(stdout);
	assert(k== -1);
      }
      if(!(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2))){
	printf("sv_array[%d,%d]:xmapid1=%d,%d,xmapid2=%d,%d (not ordered in ascending order)\n",i,i+1,sv_array[i].xmapid1,sv_array[i+1].xmapid1,sv_array[i].xmapid2,sv_array[i+1].xmapid2);
	fflush(stdout);
	assert(sv_array[i].xmapid1 <= sv_array[i+1].xmapid1 && (sv_array[i].xmapid1 < sv_array[i+1].xmapid1 || sv_array[i].xmapid2 <= sv_array[i+1].xmapid2));
      }

      /* make sure xmapid is still filtered when replaced by the new filtered (compacted) xmapid values, but only for those SV's that will be output */
      if(sv_array[i].duplicate || sv_array[i].confidence == 0 || (sv_array[i].confidence > 0 && sv_array[i].confidence < smap_min_conf))
	continue;
      if(sv_array[i+1].duplicate || sv_array[i+1].confidence == 0 || (sv_array[i+1].confidence > 0 && sv_array[i+1].confidence < smap_min_conf))
	continue;

      if(!(xmapidMap[sv_array[i].xmapid1] <= xmapidMap[sv_array[i+1].xmapid1] && (xmapidMap[sv_array[i].xmapid1] < xmapidMap[sv_array[i+1].xmapid1] || xmapidMap[sv_array[i].xmapid2] <= xmapidMap[sv_array[i+1].xmapid2]))){
	printf("sv_array[%d,%d]:xmapid1=%d,%d,xmapid2=%d,%d, xmapid1filt=%d,%d,xmapid2filt=%d,%d, output1=%d,%d output2=%d,%d(not ordered in ascending order)\n",i,i+1,
	       sv_array[i].xmapid1,sv_array[i+1].xmapid1,sv_array[i].xmapid2,sv_array[i+1].xmapid2,
	       xmapidMap[sv_array[i].xmapid1],xmapidMap[sv_array[i+1].xmapid1],xmapidMap[sv_array[i].xmapid2],xmapidMap[sv_array[i+1].xmapid2],
	       usexmapentries[sv_array[i].indices[0]]->output,usexmapentries[sv_array[i+1].indices[0]]->output,usexmapentries[sv_array[i].indices[1]]->output,usexmapentries[sv_array[i+1].indices[1]]->output);
	fflush(stdout);

	if(xmapidMap[sv_array[i].xmapid1] > 0 && xmapidMap[sv_array[i+1].xmapid1] > 0 && xmapidMap[sv_array[i].xmapid2] > 0 && xmapidMap[sv_array[i+1].xmapid2] > 0)
 	  assert(xmapidMap[sv_array[i].xmapid1] <= xmapidMap[sv_array[i+1].xmapid1] && (xmapidMap[sv_array[i].xmapid1] < xmapidMap[sv_array[i+1].xmapid1] || xmapidMap[sv_array[i].xmapid2] <= xmapidMap[sv_array[i+1].xmapid2]));
      }
    }
  }

  // NOTE All SVs and 'end' svs are again traversed to actually do the writing to smap, with the correct filtered xmapid values
  if(VERB){
    printf("writing sv to smap: numsv=%i\n", numsv);
    fflush(stdout);
  }
  if(VERB){
    printf("Calling writeSVToSmap() with dowrite=true : numsv = %d\n",numsv);
    fflush(stdout);
  }

  int nsv = 1; //increment only if pass all requirements
  for(int indp=0; indp < numsv; indp++)
    nsv += sv_array[indp].writeSVToSmap(usexmapentries, nsv, fp, true);


  //all the indels are done--check to make sure no more 'end' svs at the begining of contigs were skipped (because that contig has no indels)
  if(VERB>=2){
    printf("Calling checkEndSV() for contigs without indels : smapsize=%d,dowrite=true\n",smapsize);
    fflush(stdout);
  }
  for(int i = 0; i < smapsize; i++)  
    nsv += checkEndSV(fp, nsv, usexmapentries[i]->qrycontigid, contigendalign, endalignsize, true, true, xmapidMap, xmapindexMap);

  // do the 'end' svs at the starts of contigs with indels
  if(VERB>=2){
    printf("Calling checkEndSV() for the starts of contigs with indels: endalignsize=%d,dowrite=true\n",endalignsize);
    fflush(stdout);
  }
  for(int i=0; i < endalignsize; i++) 
    nsv += checkEndSV(fp, nsv, contigendalign[i].contigid, contigendalign, endalignsize, true, true, xmapidMap, xmapindexMap);

  // lastly, do the 'end' svs at the ends of contigs with indels
  if(VERB>=2){
    printf("Calling checkEndSV() for the ends of contigs with indels: endalignsize=%d,dowrite=true\n",endalignsize);
    fflush(stdout);
  }
  for(int i=0; i < endalignsize; i++)
    nsv += checkEndSV(fp, nsv, contigendalign[i].contigid, contigendalign, endalignsize, false, true, xmapidMap, xmapindexMap);

  FILEclose(fp);
  if(VERB){
    printf("Finished all sv (wall time=%.6f)\n", wtime());
    printf("Generated %s with %d SVs\n",filename, nsv-1);
    fflush(stdout);
  }

  // usexmapentries points to global xmapentries which we may want to use subsequently, so don't delete the actual xmapEntries just the pointer array usexmapentries[]
  delete [] usexmapentries; usexmapentries = NULL;

  delete [] indices1;
  delete [] indices2;
  if( endalignsize ) //if allocated
    delete [] contigendalign;
  delete [] xmapidMap;
  delete [] xmapindexMap;

} //end output_smap

extern char *xmap_filename;

/* output a list of SVs using the .smap format to <prefix>.indel */
void output_indel(char *basename, xmapEntry *SVs, int SVcnt)
{
  if(strstr(basename,"/dev/null"))
    return;

  char filename[BUFSIZ];
  sprintf(filename,"%s.indel",basename);
   
  if(checkFile(filename))
    return;

  if(VERB){
    printf("Generating %s (with %d indels)\n",filename, SVcnt);
    fflush(stdout);
  }

  /* open output file */
  FILE *fp;
  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("failed to open file %s for writing %d indels:errno=%d:%s\n",filename, SVcnt,eno,err);
    exit(1);
  }

  /* write out commandline */
  printversion(fp);

#ifndef WIN32
  /* get current directory pathname */
  char cwd[BUFSIZ];
  if(getcwd(cwd,BUFSIZ) == NULL){
    char error[BUFSIZ];
    char* err = strerror_r(errno, error, BUFSIZ);
    printf("getcwd() failed (errno=%d):%s\n",errno,err);
    exit(1);
  }
#endif

  fprintf(fp,"# SMAP file Version:\t0.0\n");
#ifndef WIN32
  if(vfixx_filename[num_files][0]=='/')
    fprintf(fp,"# Reference Maps From:\t%s\n", vfixx_filename[num_files]);
  else
    fprintf(fp,"# Reference Maps From:\t%s/%s\n", cwd, vfixx_filename[num_files]);
#else
    fprintf(fp,"# Reference Maps From:\t%s\n", vfixx_filename[num_files]);
#endif
#ifndef WIN32
  if(vfixx_filename[num_files+1][0] == '/')
    fprintf(fp,"# Query Maps From:\t%s\n", vfixx_filename[num_files+1]);
  else
    fprintf(fp,"# Query Maps From:\t%s/%s\n", cwd, vfixx_filename[num_files+1]);
#else
    fprintf(fp,"# Query Maps From:\t%s\n", vfixx_filename[num_files+1]);
#endif
#ifndef WIN32
  if(xmap_filename[0] == '/')
    fprintf(fp,"# Xmap Entries From:\t%s\n", xmap_filename);
  else
    fprintf(fp,"# Xmap Entries From:\t%s/%s\n", cwd,xmap_filename);
#else
    fprintf(fp,"# Xmap Entries From:\t%s\n", xmap_filename);
#endif
  fprintf(fp,"#h SmapEntryID\tQryContigID\tRefcontigID1\tRefcontigID2\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tType\tXmapID1\tXmapID2\n");
  fprintf(fp,"#f int        \tint        \tint         \tint         \tfloat      \tfloat    \tfloat      \tfloat    \tstring     \tfloat \tstring \tint\tint\n");

  for(int i = 0; i < SVcnt; i++)
    fprintf(fp,"%d\t%lld\t%lld\t%lld\t%10.1f\t%10.1f\t%10.1f\t%10.1f\t%c\t%10.2f\t%s\t%d\t%d\n",
	    i+1, // SmapEntryID
	    SVs[i].qrycontigid, // QryContigID
	    SVs[i].refcontigid, // RefContigID1
	    SVs[i].refcontigid, // RefContigID2
	    SVs[i].qrystartpos*1000.0,
	    SVs[i].qryendpos*1000.0,
	    SVs[i].refstartpos*1000.0,
	    SVs[i].refendpos*1000.0,
	    SVs[i].orientforward ? '+' : '-', // Orientation
	    SVs[i].confidence,                // Corrected Confidence
	    SVs[i].type,                      // Type
	    SVs[i].xmapid, SVs[i].xmapid);    // XmapID1, XmapID2
  
  FILEclose(fp);
}
