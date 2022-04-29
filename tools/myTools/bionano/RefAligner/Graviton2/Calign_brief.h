#ifndef CALIGN_BRIEF_H
#define CALIGN_BRIEF_H

class Calign;// forward declaration

/* abreviated version of Calign : only refid,orientation,scaleID,score,logPV,sites1[0],sites1[numpairs-1] are required for align and nalign and could be saved without calling new/delete */
class Calign_brief {
public:
  double score; ///< score of best alignment 
  double logPV; ///< -log10(Pvalue), Pvalue without Bonferroni correction
  int sites1_0;
  int sites1_M;
  int numpairs;
  int orientation; ///< orientation of 2nd map in best alignment (0 is normal, 1 = reversed) 
  int mapid1; ///< index into map[] of 1st map (or index into refmap[] of reference map) 
  int mapid2; ///< index into map[] of 1st map (or index into refmap[] of reference map) 
  int scaleID; ///< 2nd map was scaled by ScaleFactor[scaleID] (ScaleFactor[0] == 1.0) 
	
  Calign_brief() {
    reset();
  }
	
  void inline reset() {
    score = MINSCORE;
    numpairs = -1;
    mapid1 = -1;
    mapid2 = -1;
    scaleID = 0;
  }
	
  void inline update(const Calign_brief &align) {
    score = align.score;
    numpairs = align.numpairs;
    orientation = align.orientation;
    mapid1 = align.mapid1;
    mapid2 = align.mapid2;
    scaleID = align.scaleID;
    logPV = align.logPV;
    sites1_0 = align.sites1_0;
    sites1_M = align.sites1_M;		
  }
  void update(const Calign &align, int M);
};

#include "ident.h"
static Ident Calign_brief_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Calign_brief.h 5409 2016-10-03 19:32:08Z tanantharaman $");

#endif
