#ifndef XMAPENTRY_H
#define XMAPENTRY_H

#define MASK(n) (int)(~((~(0U))<<(n))) // rightmost n bits set to 1 

//#include "Cmap.h"
//#include "Calign.h"

extern double smap_Repeat_Tolerance;
extern int smap_Repeat_Minele;
extern double smap_Repeat_ConfPenalty;

#define QRY_REPEAT 1 /* HERE HERE : make commandline parameter */
#define REPEAT_TOLERANCE smap_Repeat_Tolerance /* see -svRepeat */
#define REPEAT_MINELE smap_Repeat_Minele /* see -svRepeat */
#define REPEAT_CONF_PENALTY smap_Repeat_ConfPenalty /* see -svRepeat */

class Cmap;// forward declaration
class Calign;// forward declaration
class CalignA;// forward declaration

extern Cmap **YYmap,**XXmap;
extern int numY,numX;

#ifdef _OPENMP
//#include <omp.h>
#endif


/** intermediate class to store xmap entries since there's a lot of processing between the 
 * alignment object (globals.h) and the xmap data, but the latter is used for the sv analysis. */
class xmapEntry {
public:
  int smapentryid; //used in input_svcheck.cpp. Also used in first part of output_smap() as index into usexmapentries[] 
  int xmapid; // xmapid - 1 is index into xmapentries[0..numxmapentries-1] : xmapid is the line number in the _full.xmap output file
  int xmapidfilt; //new id for filtered xmap : xmapidfilt is the line number in the final .xmap output file
  long long qrycontigid;
  long long refcontigid;
  double qrystartpos;/* absolute position, aligned with refstartpos */
  double qryendpos;/* absolute position, aligned with refendpos */
  int qrystartidx;// absolute index, aligned with refstartidx : If !orientforward then qrystartidx > qryendidex 
  int qryendidx;// absolute index, aligned with refendidx : If !orientforward then qrystartidx > qryendidex 
  double refstartpos;
  double refendpos;
  int refstartidx;
  int refendidx;
  bool orientforward; /**< true means orientforward == (!align->orientation) */
  int dominant_orientation;/* set in getSameQueryIndices() : Based on all MGs between same querycontigid,refcontigid, the dominant query orientation (based on smap_PrimaryOrientation ratio)
			      0 : No dominant orientation
			      1 : Dominant query orientation is forward
			      -1 : Dominant query orientation is backwards */
  double confidence; /**< log10 based confidence */
  char *type;/* Type field (read in by input_svcheck() from .indel or .smap) */
  //cigar string is not stored because it is not used in sv analysis (output_smap), but see field align below that contains full alignment
  //extra information for finding unaligned regions at start and end of contigs
  int startunalign; /**< number of un-aligned labels from beginning of query contig to first aligned site of this alignment */
  int endunalign; /**< number of un-aligned labels from end of query contig to last aligned site of this alignment */
  double querycontiglen; /**< length of query contig (kb) */
  double refcontiglen;/**< length of ref contig (kb) */
  int querycontignsite; // n sites in query contig
  int refcontignsite;// n sites in ref contig
  int numrefsites;// number of ref sites in alignment region : >= (refendidx - refstartidx + 1)
  int numqrysites;// number of qry sites in alignment region : >= abs(qryendidx - qrystartidx) + 1
  double* refsites; // refsites[0..numrefsites-1] : for checking repeats
  double* qrysites; // qrysites[0..numqrysites-1] : for checking repeats
  bool have_ref_repeat; //whether self.callRepeat found any repeats in ref (or qry if QRY_REPEAT)
  bool sv_overlap; //output_smap: overlap with another matchgroup (scheduled to be filtered out)
  bool inversion_only;/* this matchgroup is flagged to be filtered (similar to sv_overlap) and can be used for inversion calls only */
  bool largeMG_only;/* this matchgroup is flagged to be filtered (similar to sv_overlap) and can only be used in SV calls as the larger matchgroup */
  int xmapid_overlapped;/* If > 0 : this matchgroup is completely overlapped by xmapid == xmapid_overlapped : only consider SV call with this other MG xmapid_overlapped */
  bool smallDupInv;/* If xmapid_overlapped > 0 : true IFF a small duplication,inverted-duplication or inversion was called with this MG */

  int mergecnt, mergestart;
  
  Cmap *QryMap;
  Cmap *RefMap;

  int sbin;/* indel size bin : > 0 for insertions, < 0 for deletions */
  bool Palindromic;// If reverse alignment is almost as good
  bool output;// If this matchgroup has an SV that was output

  Calign *align;/* full alignment : used when reading in XMAP 0.2 or for small inversions or duplication */

  int jMerge;/* temporary index into usexmapentries[jMerge] to MG that was merged with the current MG */
  xmapEntry *orig;/* temporary pointer to original MG before merging with usexmapentries[jMerge] */

  /* member functions */
  //constructors
  xmapEntry() {
    xmapid      = 0;
    xmapidfilt  = 0;
    qrycontigid = 0;
    refcontigid = 0;
    qrystartpos = 0;
    qryendpos   = 0;
    qrystartidx = 0;
    qryendidx   = 0;
    refstartpos = 0;
    refendpos   = 0;
    refstartidx = 0;
    refendidx   = 0;
    orientforward = true;
    confidence  = 0;
    startunalign = 0;
    endunalign = 0;
    querycontiglen = 0;
    numrefsites = 0;
    numqrysites = 0;
    refsites = NULL;
    qrysites = NULL;
    have_ref_repeat = false;
    sv_overlap = false;
    Palindromic = false;

    output = false;
    inversion_only = false;
    largeMG_only = false;// NEW11
    xmapid_overlapped = -1;
    smallDupInv = false;

    align = NULL;
    QryMap = RefMap = NULL;
  }
  xmapEntry(Calign* clgn, size_t xid, long long qid, long long rid, double qspos, double qepos, int qsi, int qei, double rspos, double repos, int rsi, int rei, bool ofor, double conf, int startua=0, int endua=0, double qlen=0, int nrs=0, int nqs=0, int N = 0) {
    if(xid > MASK(31)){
      printf("ERROR in xmapEntry():xid=%lu exceeds range of int = %d\n",xid,MASK(31));
      exit(1);
    }
    align = clgn;
    xmapid      = xid;
    xmapidfilt  = 0;
    qrycontigid = qid;
    refcontigid = rid;
    qrystartpos = qspos;
    qryendpos   = qepos;
    qrystartidx = qsi;
    qryendidx   = qei;
    refstartpos = rspos;
    refendpos   = repos;
    if(DEBUG/* HERE HERE >=2 */) assert(refstartpos <= refendpos);
    refstartidx = rsi;
    refendidx   = rei;
    if(DEBUG/* HERE HERE >=2 */) assert(refstartidx <= refendidx);
    orientforward  = ofor;
    if(DEBUG/* HERE HERE >=2 */) assert(orientforward ? (qrystartidx <= qryendidx) : (qrystartidx >= qryendidx));
    //    if(DEBUG/* HERE HERE >=2 */) assert(align!=NULL && align->orientation == (orientforward ? 0 : 1));
    confidence     = conf;
    startunalign   = startua;
    endunalign     = endua;
    querycontiglen = qlen;
    querycontignsite = nqs;
    refcontignsite = N;
    //    if(DEBUG>=3 && !XXmap[align->mapid2]->origmap) assert(querycontignsite == XXmap[align->mapid2]->numsite[0]);
    numrefsites    = nrs;
    if(DEBUG/* HERE HERE >=2 */) assert(numrefsites >= refendidx - refstartidx + 1);
    refsites       = new double[numrefsites]; //best to allocate in constructor and fill outside rather than both outside
    numqrysites = abs(qryendidx - qrystartidx) + 1;
    qrysites = new double[numqrysites];
    if(DEBUG>=2){
      double aNaN = nan("NaN");
      for(int t = 0; t < numrefsites; t++)
	refsites[t] = aNaN;
      for(int t = 0; t < numqrysites; t++)
	qrysites[t] = aNaN;
    }
    have_ref_repeat = false;
    sv_overlap = false;
    output = false;
    inversion_only = false;// NEW
    largeMG_only = false;// NEW11

    xmapid_overlapped = -1;
    if(VERB>=2 && (xmapid == 866 || xmapid == 868)){
      printf("xmapEntry():xmapid=%d, refid=%lld, qryid=%lld, xmapid_overlapped = %d\n",xmapid, refcontigid, qrycontigid, xmapid_overlapped);
      fflush(stdout);
    }
    
    smallDupInv = false;

    Palindromic = false;
    QryMap = RefMap = NULL;
  }
  ~xmapEntry() {
    delete [] refsites; refsites = NULL;
    delete [] qrysites; qrysites = NULL;
  }

  void callRepeat(bool verbose=false, double tolerance=REPEAT_TOLERANCE, int minele=REPEAT_MINELE, double confPenalty= REPEAT_CONF_PENALTY ); //defined in simpleRepeat.cpp
};

extern void input_xmap(char *filename, xmapEntry *&xmap, int &xmapCnt, int &xmapMax, int &StartQuery, int &NummapsQuery, int &StartRef, int &NummapsRef, char* &query_filename, char * &ref_filename);/* see input_svcheck.cpp */
extern void output_xmapIO(char *filename, Calign **Lalignment, size_t start, size_t end, char *basename, char * &queryfilename, int Indel, xmapEntry **xmapentries, double Xscale, double Yscale, bool full);

extern xmapEntry **xmapentries;/* array which corresponds to entries in the xmap--filled in output_xmap */

#include "ident.h"

static Ident xmapEntry_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/xmapEntry.h 11317 2020-07-15 20:27:50Z tanantharaman $");

#endif
