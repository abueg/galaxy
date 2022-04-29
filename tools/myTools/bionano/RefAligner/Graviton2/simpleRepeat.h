#ifndef SIMPLEREPEAT_H
#define SIMPLEREPEAT_H

#include "Cmap.h"

extern float stdDev(double* inlist, int nin, double mean);

//class to store repeat data
//instances will be created by findRepeatSimple function
//so repeat data is defined relative to its arguments
class repeat {
public:
  double startpos; //first element which is repeat
  double endpos;
  int    nele;    //number of repeat elements
  double avgsize; //average size of all repeat elements
  float stddev; //standard deviation of repeat units
  long long id; //molecule ID: 0 if not using molecules
  FLOAT length; //molecule length: 0 if not using molecules (see Cmap.h)
  //data members for FP/FN
  int  fpIndices_len;
  int* fpIndices    ;
  int  fnIndices_len;
  int* fnIndices    ;


  repeat() {
    startpos = 0;
    endpos = 0;  
    nele   = 0;    
    avgsize = 0;
    stddev = 0;
    id = 0;
    length = 0;
    fpIndices_len = 0;
    fpIndices     = 0;
    fnIndices_len = 0;
    fnIndices     = 0;
  }
  repeat(double start, double end, int ne, double avg) {
    startpos = start;
    endpos   = end;  
    nele     = ne;    
    avgsize  = avg; 
    stddev = 0;
    id = 0;
    length = 0;
    fpIndices_len = 0;
    fpIndices     = 0;
    fnIndices_len = 0;
    fnIndices     = 0;
  }
  //constructor for FP/FN: copy input FP/FN indices arrays to data members
  repeat(double start, double end, double* intervals, int nintervals, int* fpind, int nfpind, int* fnind, int nfnind) {
    bool debug = false; //true
    startpos = start;
    endpos   = end;  
    nele     = nintervals+nfnind; //fn count as two intervals, but only one has been accounted for in nintervals
    avgsize  = 0; 
    for( int i=0; i<nintervals; i++ ) {
      if(debug) printf("  %.4f", intervals[i]);
      avgsize += intervals[i];
    }
    avgsize /= nintervals;
    stddev = stdDev(intervals, nintervals, avgsize); 
    if(debug) printf("\n  mean=%f, stddev=%f\n", avgsize, stddev);

    id = 0;
    length = 0;
    fpIndices_len = nfpind;
    fpIndices = (int*)malloc(fpIndices_len*sizeof(double));
    for( int i=0; i<fpIndices_len; i++ ) {
      if(debug && i==0) printf("  FP: ");
      if(debug) printf("%i ", fpind[i]);
      fpIndices[i] = fpind[i];
    }
    if(debug && fpIndices_len ) printf("\n");
    fnIndices_len = nfnind;
    fnIndices = (int*)malloc(fnIndices_len*sizeof(double));
    for( int i=0; i<fnIndices_len; i++ ) {
      if(debug && i==0 ) printf("  FN: ");
      if(debug) printf("%i ", fnind[i]);
      fnIndices[i] = fnind[i];
    }
    if(debug && fnIndices_len ) printf("\n");
  } //end constructor for FP/FN
}; //end class repeat


extern void findRepeatSimpleMaps(int filtertype, float tolerance, int minele, bool verbose=false);
extern int  findRepeatSimple(double* inlist, int inlist_len, float tolerance=0.1, int minele=3 );
extern int  findSimpleRepeatFPFN(double* labelPosList, int labelPosList_len, float tolerance=0.1, int minele=3 );
extern void output_repeat(char *basename);
extern void output_repeatStats(char *basename);
extern void output_repeatAll(int fitertype=0, Cmap **partmaps=0, int numpartmaps=0);
extern void maskRepeatLabels(Cmap*& inmap, int repeat_start, int repeat_end, bool verbose=false);

static Ident simpleRepeat_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/simpleRepeat.h 3381 2014-11-06 02:44:41Z tanantharaman $");

#endif
