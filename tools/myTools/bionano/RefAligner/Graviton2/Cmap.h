#ifndef CMAP_H
#define CMAP_H

#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "Calign_brief.h"

class Ccontig;/* forward declaration : cannot include "Ccontig.h" until after Cmap.h is declared */

static Ident Cmap_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Cmap.h 10584 2020-02-12 02:25:46Z tanantharaman $");

#define CHR_OFFSET 1000000000LL /* bp location offset as multiple of chromosome number (used when computing absolute reference location : starts each chr at chr * CHR_OFFSET) */

/** Primary consensus nanomap object */
class Cmap {
 public:
  int fileid;/**< file from which Map was read in is vfixx_filename[fileid] */
  int linenum;/**< line number at which Map was read in vfixx_filename[fileid] */
  int mapid;/**< original index into array map[] (or refmap[]) after maps were read in and after Assembler_input (or as modified by pairalign) */

  long long id;/**< id of Map (must be unique accross all input files) */

  char *name;/**< id name : only used for .vfixx files */

  // NOTE number of colors is given by global "colors"
  char *Nickase[MAXCOLOR];/**< Nickase[c] recognition sequence for color c */
  int numsite[MAXCOLOR];/**< numsite[c=0..colors-1] : number of sites for color c in this map (or -1 if site[c] is not allocated) */
  FLOAT *site[MAXCOLOR];/**< site[c=0..colors-1][i=0..numsite[c]+1] : the location of site i in color c (in kb from the left end site[0])
		  * site[0] is always 0, and site[numsite+1] is the length of this map */
  int blockmem;/**< If >= 1 : site[] is allocated in block memory (see maxmapalloc() and mapcompact() in Cmap.cpp) and has memory allocated for blockmem-2 sites
		       <= -1 : site[],siteSD[],sitecov[],ChimQuality[],sitecnt[],SNRcnt[],SNRdist[],SNRgmean[],lnSNRsd[] are (copies of) (possibly) block-allocated (orig*) arrays */
  int newmapid;/**< Filler, Used in input_contigs() to keep only maps that will actually be used */
  FLOAT *rawsite[MAXCOLOR];/* site locations before -resbias bias removal : only defined if NoStat==0 && maxresbias > mres * 0.5 */
  int rblockmem;/**< 1 IFF rawsite[] is allocated in block memory (see rawsitealloc() in Cmap.cpp) */

  /* The following fields are only defined for .bnx input (see also SNR distribution below, which is read in from .cmap)*/
  long long *truesiteL[MAXCOLOR];/**< leftmost true location truesiteL[c][1..numsite[c]] of nicking sites in Reference Genome (-1 for false positives) */
  long long *truesiteR[MAXCOLOR];/**< rightmost true locations truesiteR[c][1..numsite[c]] of nicking sites in Reference Genome (-1 for false positives) */
  // NOTE : If Tracksite, truesiteL[],truesiteR[] refer to interleaved label index in the original Query or Ref map (QueryOrigin or RefOrigin) */                         

  double *SNR[MAXCOLOR];/**< SNR[c][i=1..numsite[c]] Signal to Noise ratio of site i : Only valid if MapSNR==1 */
  double *Intensity[MAXCOLOR];/* Nearby Average Intensity for each site : Only valid if MapIntensity==1 */
  int *stitch[MAXCOLOR];/**< stitch[c][i=1..numsite[c]+1] : If image was stitched between site i and i-1 : Only valid if MapStitched==1 */
  double *stitchLocation[MAXCOLOR];/* backbone stitch point locations stitchLocation[c][0..NumStitch[c]-1] : Only valid if MapStitchLoc==1 */
  double *PSFWidth[MAXCOLOR];/* PSF Width for each site : Only valid if MapPSFWidth==1 */
  int *ImageFOV[MAXCOLOR];/* Image Field Of View for each site : Only valid if MapImageLoc==1 */
  double *ImageX[MAXCOLOR];/* Image X coordinate for each site : Only valid if MapImageLoc==1 */
  double *ImageY[MAXCOLOR];/* Image Y coordinate for each site : Only valid if MapImageLoc==1 */

  /* The following fileds are only defined for Version 1.0 .bnx input (Only valid if BNXVersion >= 1) */
  int OriginalMoleculeId;/* Molecule ID from original MOL file (-1 if unknown) */
  int ScanNumber;/* Scan number from run (1 or larger) */
  int ScanDirection;/* Scan Direction (-1 if unknown) */
  char *ChipId;/* Serial number & date of the chip of the run */
  int FlowCell;/* Flow Cell number */
  double AvgIntensity;/* average Backbone intensity */
  double bbSNR;/* average backbone SNR */
  int NumStitch[MAXCOLOR];

  int RunIndex;/* index into RunDataList[0..RunDataListLen-1] (RunIndex+1 is present in column "RunId" for extended version of BNX 1.0, optional if only one "Run Data" field is present implying RunIndex==0) */

  /* following 7 optional scalars are only present if BNXVersion = 1, BNXMinorVersion = 1 (1.3) and BNXheaderFXY = 1. They will be present after the "RunId" field (which has RunIndex + 1) and before MapWT */
  int Column;
  int StartFOV, StartX, StartY;
  int EndFOV, EndX, EndY;

  /* following 6 optional scalars are only present if BNXVersion = 1, BNXMinorVersion = 1 and BNXheaderChim = 1. They will be present after the "RunId" field (and after optional BNXheaderFXY fields) and before optional MapWT */
  int FeatureMinPos;/* location in bp */
  int FeatureMaxPos;
  double FeatureMinScore;
  double FeatureMaxScore;
  double FeatureMinScaled;
  double FeatureMaxScaled;

  int UniqueScanId;/* unique mapping from (RunIndex,ScanNumber), 0 based. UniqueScanId+1 present in extended version of BNX1.0 as GlobalScanNumber, but never read in since it can be computed from other fields and changes after -merge. */
  int RunId;/* Time Order of this map's RunIndex amongst all RunDataList[] entries with same ChipId,FlowCell as this map (0 based)
	       NOTE : there is a unique mapping from RunIndex to RunId, but different RunIndex can have the same RunId (if they have different ChipId or FlowCell)
	       This is only used to create a unique map id and is computed from other fields */
  int ScanId;/* unique mapping from (RunId,ScanNumber) : only used to create a unique id and computed from other fields (0 based) */

  /* The following fields are only defined for .cmap input */
  double *siteSD[MAXCOLOR];/**< siteSD[c][i = 1 .. numsite[c]] : StdDev of site interval site[c][i..i+1] */
  float *sitecov[MAXCOLOR];/**< sitecov[c][i = 1 .. numsite[c]] : underlying total coverage of site site[c][i] (may be weighted sum) 
			    NOTE : If FRAGCOV_FIX, then sitecov[c][i] is the interval coverage between label i and i+1 */
  float *sitecnt[MAXCOLOR];/**< sitecnt[c][i = 1 .. numsite[c]] : underlying number of expressed labels of site location site[c][i] (may be weighted sum) */

  /* the following field is only defined for .cmap consensus map input with ChimQuality column (not defined if quality[c] == 0) */
  float *ChimQuality[MAXCOLOR];/**< If quality[c] != 0 : quality[c][i = 1 .. numsite[c]] : Refine quality scores at site location site[c][i] */
  float *ChimNorm[MAXCOLOR];/* N1+N2+N3 */
  float *SegDupL[MAXCOLOR];
  float *SegDupR[MAXCOLOR];
  float *FragileEndL[MAXCOLOR];
  float *FragileEndR[MAXCOLOR];
  float *OutlierFrac[MAXCOLOR];

  float *FragSd[MAXCOLOR];
  float *ExpSd[MAXCOLOR];
  float *FragCov[MAXCOLOR];
  float *FragChiSq[MAXCOLOR];
  
  /* the following field is only defined for .cmap consensus map input with Mask column. The following bit-wise values are defined:
     END_NOEXT : For 1st and End label this denotes that the corresponding contig end will not be extended during refinement (generated by -CloneLimitLen)
 */
  size_t *Mask[MAXCOLOR];// NOTE : Not always block allocated

  /* The following fields are only defined for .cmap consensus map input with underlying SNR values : 
     If SNRcnt[c] == 0 the corresponding SNRgmean[c][i] lnSNRsd[c][i] and lnSNRdist[c][i] are invalid (not present in .cmap input).
     if SNRcnt[c][i=0..numsite[c]+1] == 0, only SNRgmean[c][i] and lnSNRsd[c][i] are valid  */
  int *SNRcnt[MAXCOLOR];/* number of SNR values underlying each consensus map site */
  double *SNRgmean[MAXCOLOR];/* geometric mean of SNR values underlying each consensus map site */
  double *lnSNRsd[MAXCOLOR];/* sd of log(SNR) values underlying each consensus map site */
  double **SNRdist[MAXCOLOR];/* actual SNR values underlying each consensus map site */

  /* The following fields are used for map fragments */
  Cmap *origmap;/**< pointer to map from which this map fragment is derived (may be another map fragment with origmap != 0) :
		   * (0 if this map is not a map fragment) */
  int left[MAXCOLOR],right[MAXCOLOR];/**< leftmost and rightmost site of origmap included in this map :
					This maps has intervals of size MININTERVAL on either end (from the first/last site of either color, if colors > 1),
					unless left/right corresponds to the ends of the original maps (0 or origmap->numsite[]+1 respectively)
					NOTE : If colors > 1, one of the colors may have no sites in this fragment and left[c] == 1 and right[c] == 0
				     */

  int centercloned;/* > 0 : If aligned center of this map was cloned : Used with PairSplit to ignore this original map.
		            Used with -pairmerge to mark maps whose orig* fields are invalid/undefined */
  int paired;/* If > 0 : If this map has been aligned : Used with PairSplit to prevent further splitting of this map with lower scoring alignments
		Used with -pairmerge to represet the 1 + index into map[] of the new map into which this map has been merged. Also prevents output of this map */
  
  long long origid;/* used with -pairmerge to keep track of id or original map : normally this is the same as id unless -SplitSegDup was used */
  int Xid;/* mapid of root query map (root of origmap->origmap-> chain) of this map (or aligned with this map for reference maps) : Used with PairSplit to handle multiple query maps together */

  double y2sum, xysum;/* used to compute inscale = y2sum / xysum */
  double incscale;/**< incremental scaling factor (to be multiplied into scale and site[][] values later) : same for all colors */
  double incwt;/**< relative weight of this map scaling factor (relative to other nanomaps) (sum of y for non-outlier intervals) */
#if 0
  double xsum;/* sum of x (non-outlier intervals) */
  double Xlen;/* sum of x for all intervals */
  double Ylen;/* sum of y for all intervals */
#endif

  double cumscale;/* cumulative scaling factor : used by -MapScale */

  /* The following fields are used for simulated data */
  double startloc;/**< leftmost location in specified contig of Reference Genome (in kb) (0.0 if unknown) */
  double endloc;/**< rightmost location in specified contig of Reference Genome (in kb) */
  double len;/**< length : this should equal endloc-startloc */
  double startloc2;/* If != 0.0 : leftmost location of 2nd chimeric fragment */
  double endloc2;/* If != 0.0 : rightmost location of 2nd chimeric fragment */
  double len2;/**< If != 0.0 : length of 2nd chimeric fragment : should equal endloc2-startloc2 */
  int flipped;/**< if true orientation is flipped */
  int chimflip;/**< if 2nd chimeric fragment is flipped relative to 1st */

  int refid;/**< contig id in Reference Genome */
  int fpcnt;/**< number of false positives in this molecule (-1 if unknown) */
  double Qscore;/**< reliability of startloc/endloc information (used during .map input 
		   * to track best scoring local alignment based on alignment Score) */
  double startmatch,endmatch;/**< leftmost/rightmost location that aligned in Reference Genome (in kb) (0.0 if unknown) */
  
  /* The following fields are only defined for reference maps (or if mres > 0, only fields starting at orignumsite[] are defined, IFF origsite[c] != NULL)
     These values are only used for .cmap output, hence do not include SNR[],Intensity[] or stitch[] information (unless minSNRestimate is true) */
  char *refname;/**< NULL if not a reference map */
  int *remap[MAXCOLOR];/**< remap[c][j=1..orignumsite[c]] : the new location of original site j in site[c][] */
  int *NickCnt[MAXCOLOR];/**< NickCnt[c][j=1..orignumsite[c]]] : the number of original locations mapped to same location */
  int orignumsite[MAXCOLOR];/**< orignumsite[c] : original number for sites for color c in this Reference amp */
  double *origsite[MAXCOLOR];/**< origsite[c][j=0..orignumsite[c]+1] : original location of site j in color c */
  double *origsiteSD[MAXCOLOR];
  float *origsitecov[MAXCOLOR];
  float *origsitecnt[MAXCOLOR];

  float *origChimQuality[MAXCOLOR];
  float *origChimNorm[MAXCOLOR];
  float *origSegDupL[MAXCOLOR];
  float *origSegDupR[MAXCOLOR];
  float *origFragileEndL[MAXCOLOR];
  float *origFragileEndR[MAXCOLOR];
  float *origOutlierFrac[MAXCOLOR];

  float *origFragSd[MAXCOLOR];
  float *origExpSd[MAXCOLOR];
  float *origFragCov[MAXCOLOR];
  float *origFragChiSq[MAXCOLOR];

  size_t *origMask[MAXCOLOR];

  int *origSNRcnt[MAXCOLOR];
  /*     If origSNRcnt[c] == 0 the corresponding origSNRgmean[c][i] origlnSNRsd[c][i] and origlnSNRdist[c][i] are invalid (not present in .cmap input).
        if origSNRcnt[c][i=0..orignumsite[c]+1] == 0, only origSNRgmean[c][i] and origlnSNRsd[c][i] are valid  */
  double *origSNRgmean[MAXCOLOR];
  double *origlnSNRsd[MAXCOLOR];
  double **origSNRdist[MAXCOLOR];

  long long *origtruesiteL[MAXCOLOR];
  long long *origtruesiteR[MAXCOLOR];
  double *origSNR[MAXCOLOR];
  double *origIntensity[MAXCOLOR];
  int *origstitch[MAXCOLOR];
  // NOTE : origstitchLocation[] is not needed since stitchLocation is NOT modified by mapcorrect()
  double *origPSFWidth[MAXCOLOR];
  int *origImageFOV[MAXCOLOR];
  double *origImageX[MAXCOLOR];
  double *origImageY[MAXCOLOR];

  /* The following field is only defined for query maps */
  Calign_brief align[MAXCOLOR];/**< NULL if not defined, otherwise best known alignment(s) with refid = align->mapid1 (best over all reference maps if there are multiple reference maps) */
  Calign_brief nalign[MAXCOLOR];/**< value of align for next iteration */
  /* only refid,orientation,scaleID,score,LogPV,sites1[0],sites1[numpairs-1] are required for align and nalign and can be saved without calling new/delete */

  double Likelihood;/* used by BestRefWT to sum up Likelihoods of all alignments of this map */

  /* The following field is only used by pairalign to precompute site intervals up to DELTA sites apart */
  PFLOAT **deltaX[MAXCOLOR]; /**< deltaX[c][J=2..M][n = max(0,DELTA+1-J) .. DELTA-1] = XX[c][J] - XX[c][J-DELTA+n] */
  PFLOAT **deltaR[MAXCOLOR]; /**< deltaR[c][J=2..M][n = max(0,DELTA+1-J) .. DELTA-1] = XX[c][M+1 - (J-DELTA+n)] - XX[c][M+1-J]  */
  /* NOTE : if (USE_SSE && USE_AVX || USE_MIC) && USE_PFLOAT then deltaX and deltaR are defined differently
      deltaX[c][n=1..DELTA][J=1+n...M] = XX[c][J] - XX[c][J-n]
      deltaR[c][n=1..DELTA][J=1+n...M] = XX[c][M+1 - (J-n)] - XX[c][M+1-J]
   */
  PFLOAT **deltaY[MAXCOLOR];/* same as deltaX[] assuming !(USE_SSE==1 && USE_AVX==1) */
  
  double mapWT;/* used with -BestRefWT -mapped to output map weight to BNX file */

  // leftoverlap,rightoverlap : Used in pairmerge only to mark if left or right end has an end overlap with another map and points to the corresponding alignment or NULL otherwise
  Calign *leftoverlap;
  Calign *rightoverlap;

  Ccontig *contig;/* Used with -pairmergeHmap to create an Hmap : Typically this Cmap will be the primary Allele that is used to find bubbles to add to the alternate Allele and the
		     contig will list all the component contigs as maps aligned to the Hmap (including this Cmap as the 1st map contig[0]).

		     Also used to read in .hmap : In that case this Cmap will match the contig map and contig.nummaps will be -1. 
		     If -HaploType is used, contig will be used to compute 2 Allele maps  with the first Allele saved in this Cmap and the 2nd Allele located in refmap[Allele].
		     Also nummaps will be changed to 2 and contig->contig[A=0..1] will refer to the 2 Allele Cmaps (The current Cmap and refmap[Allele]) and all of the following will be true : 
		           this == refmap[contig->contig[0].mapid]
			   id == contig->contig[0].id
			   numsite[0] == contig->contig[0].numsite[0]
			   Allele == contig->contig[1].mapid
			   refmap[Allele]->id == contig->contig[1].id (If RefineHap==1, this is to id, so -mapped -grouped stores mappings for both alleles together and is otherwise not used
			                                               If RefineHap==0, this is set to a unique id, using a mechanism like -contigsplit)
			   refmap[Allele]->numsite[0] == contig->contig[1].numsite[0]
			   contig->flip[0] == contig->flip[1] == 0

		    The mapping from Allele 1 (this Cmap's site[0][I=0..numsite[0]+1]]) to contig->site[0][i=0..n+1] is given by i == contig->sitemap[0][0][I], where n == contig->numsite[0]
		    The mapping from Allele 2 (refmap[Allele]->site[0][I=0..N2+1]) to contig->site[0][i=0..n+1] is given by i == contig->sitemap[0][1][I], where N2 == refmap[Allele]->numsite[0]
		  */
  int Allele;/* If contig != NULL : points to the 2nd Allele refmap[Allele] for which this Cmap is the 1st Allele : Both alleles are derived from this contig.
	        If contig == NULL : points to the 1st Allele refmap[Allele] for which this Cmap is the 2nd Allele (or -1 if NA) : Both alleles are derived from refmap[Allele]->contig */

  double *pMapWT;/* Used with extendWTdup == 0 : If != 0, pointer to last alignment weight of this map that was raised to extendWT (from original value = origMapWT) */
  double origMapWT;

  /**< member functions */
  Cmap()/**< : align(), nalign()*/ {
    init();
  }
  void init() {
    /*    printf("Cmap():this=x%p\n", this);
	  fflush(stdout);*/
    name = 0;
    refname = 0;

    Column = -1;

    for(int c = 0; c < MAXCOLOR; c++){
      origsite[c] = 0;

      remap[c] = 0;
      NickCnt[c] = 0;
      site[c] = 0;
      numsite[c] = -1;
//       if(sizeof(FLOAT) != sizeof(float))
      rawsite[c] = 0;
      origsiteSD[c] = siteSD[c] = 0;
      origsitecov[c] = sitecov[c] = 0;
      origsitecnt[c] = sitecnt[c] = 0;

      origChimQuality[c] = ChimQuality[c] = NULL;
      origChimNorm[c] = ChimNorm[c] = NULL;
      origSegDupL[c] = SegDupL[c] = NULL;
      origSegDupR[c] = SegDupR[c] = NULL;
      origFragileEndL[c] = FragileEndL[c] = NULL;
      origFragileEndR[c] = FragileEndR[c] = NULL;
      origOutlierFrac[c] = OutlierFrac[c] = NULL;
      
      origFragSd[c] = FragSd[c] = NULL;
      origExpSd[c] = ExpSd[c] = NULL;
      origFragCov[c] = FragCov[c] = NULL;
      origFragChiSq[c] = FragChiSq[c] = NULL;

      origMask[c] = Mask[c] = NULL;

      Nickase[c] = NULL;

      SNR[c] = NULL;
      stitch[c] = NULL;
      Intensity[c] = NULL;
      stitchLocation[c] = NULL;
      PSFWidth[c] = NULL;
      ImageFOV[c] = NULL;
      ImageX[c] = ImageY[c] = NULL;

      deltaY[c] = deltaX[c] = deltaR[c] = 0;

      origSNRcnt[c] = SNRcnt[c] = 0;

      origtruesiteL[c] = origtruesiteR[c] = truesiteL[c] = truesiteR[c] = NULL;
      origSNR[c] = origIntensity[c] = origPSFWidth[c] = NULL;
      origstitch[c] = NULL;
      origImageFOV[c] = NULL;
      origImageX[c] = origImageY[c] = NULL;
      align[c].reset();
      nalign[c].reset();
    }

    if(DEBUG) RunIndex = -1;

    blockmem = rblockmem = 0;
    id = mapid = -1;
    fileid = -1;
    flipped = 0;
    //align = nalign = 0;
    origmap = NULL;
    OriginalMoleculeId = -1;
    ScanDirection = -1;
    ScanNumber = 1;
    UniqueScanId = ScanId = 0;
    refid = 0;
    startloc = endloc = len = 0.0;
    startloc2 = endloc2 = len2 = 0.0;
    mapWT = 1.0;
    cumscale = 1.0;
    contig = NULL;
    Allele = -1;
    leftoverlap = rightoverlap = NULL;
  }
  ~Cmap(){
    allfree();
  }

  void allfree();

  void allfree2();/* free memory of all except 1st color */
  
  void colorswap(int U) {/* swap 1st and Uth color information */
    U--;
    if(DEBUG) assert(U > 0);

    char *p = Nickase[0];
    Nickase[0] = Nickase[U];
    Nickase[U] = p;

    int tmp = numsite[0];
    numsite[0] = numsite[U];
    numsite[U] = tmp;

    if(1){
      FLOAT *ptmpF = site[0];
      site[0] = site[U];
      site[U] = ptmpF;
    }

    if(rawsite[0] || rawsite[U]){
      FLOAT *ptmpf = rawsite[0];
      rawsite[0] = rawsite[U];
      rawsite[U] = ptmpf;
    }
    
    if(truesiteL[0] || truesiteL[U]){
      long long *ptmp = truesiteL[0];
      truesiteL[0] = truesiteL[U];
      truesiteL[U] = ptmp;
    }
    if(truesiteR[0] || truesiteR[U]){
      long long *ptmp = truesiteR[0];
      truesiteR[0] = truesiteR[U];
      truesiteR[U] = ptmp;
    }

    if(SNR[0] || SNR[U]){
      double *ptmpd = SNR[0];
      SNR[0] = SNR[U];
      SNR[U] = ptmpd;
    }
    if(Intensity[0] || Intensity[U]){
      double *ptmpd = Intensity[0];
      Intensity[0] = Intensity[U];
      Intensity[U] = ptmpd;
    }
    if(stitch[0] || stitch[U]){
      int *ptmp = stitch[0];
      stitch[0] = stitch[U];
      stitch[U] = ptmp;
    }
    if(stitchLocation[0] || stitchLocation[U]){
      double *ptmpd = stitchLocation[0];
      stitchLocation[0] = stitchLocation[U];
      stitchLocation[U] = ptmpd;
    }
    if(PSFWidth[0] || PSFWidth[U]){
      double *ptmpd = PSFWidth[0];
      PSFWidth[0] = PSFWidth[U];
      PSFWidth[U] = ptmpd;
    }

    if(ImageFOV[0] || ImageFOV[U]){
      int *ptmp = ImageFOV[0];
      ImageFOV[0] = ImageFOV[U];
      ImageFOV[U] = ptmp;
    }
    if(ImageX[0] || ImageX[U]){
      double *ptmp = ImageX[0];
      ImageX[0] = ImageX[U];
      ImageX[U] = ptmp;
    }
    if(ImageY[0] || ImageY[U]){
      double *ptmp = ImageY[0];
      ImageY[0] = ImageY[U];
      ImageY[U] = ptmp;
    }

    if(siteSD[0] || siteSD[U]){
      double *ptmpd = siteSD[0];
      siteSD[0] = siteSD[U];
      siteSD[U] = ptmpd;
    }
    if(sitecov[0] || sitecov[U]){
      float *ptmp = sitecov[0];
      sitecov[0] = sitecov[U];
      sitecov[U] = ptmp;
    }
    if(sitecnt[0] || sitecnt[U]){
      float *ptmp = sitecnt[0];
      sitecnt[0] = sitecnt[U];
      sitecnt[U] = ptmp;
    }
    
    if(ChimQuality[0] || ChimQuality[U]){
      float *ptmp = ChimQuality[0];
      ChimQuality[0] = ChimQuality[U];
      ChimQuality[U] = ptmp;
    }
    if(ChimNorm[0] || ChimNorm[U]){
      float *ptmp = ChimNorm[0];
      ChimNorm[0] = ChimNorm[U];
      ChimNorm[U] = ptmp;
    }
    if(SegDupL[0] || SegDupL[U]){
      float *ptmp = SegDupL[0];
      SegDupL[0] = SegDupL[U];
      SegDupL[U] = ptmp;
    }
    if(SegDupR[0] || SegDupR[U]){
      float *ptmp = SegDupR[0];
      SegDupR[0] = SegDupR[U];
      SegDupR[U] = ptmp;
    }
    if(FragileEndL[0] || FragileEndL[U]){
      float *ptmp = FragileEndL[0];
      FragileEndL[0] = FragileEndL[U];
      FragileEndL[U] = ptmp;
    }
    if(FragileEndR[0] || FragileEndR[U]){
      float *ptmp = FragileEndR[0];
      FragileEndR[0] = FragileEndR[U];
      FragileEndR[U] = ptmp;
    }
    if(OutlierFrac[0] || OutlierFrac[U]){
      float *ptmp = OutlierFrac[0];
      OutlierFrac[0] = OutlierFrac[U];
      OutlierFrac[U] = ptmp;
    }
    if(FragSd[0] || FragSd[U]){
      float *ptmp = FragSd[0];
      FragSd[0] = FragSd[U];
      FragSd[U] = ptmp;
    }
    if(ExpSd[0] || ExpSd[U]){
      float *ptmp = ExpSd[0];
      ExpSd[0] = ExpSd[U];
      ExpSd[U] = ptmp;
    }
    if(FragCov[0] || FragCov[U]){
      float *ptmp = FragCov[0];
      FragCov[0] = FragCov[U];
      FragCov[U] = ptmp;
    }
    if(FragChiSq[0] || FragChiSq[U]){
      float *ptmp = FragChiSq[0];
      FragChiSq[0] = FragChiSq[U];
      FragChiSq[U] = ptmp;
    }

    if(Mask[0] || Mask[U]){
      size_t *ptmp = Mask[0];
      Mask[0] = Mask[U];
      Mask[U] = ptmp;
    }

    if(SNRcnt[0] || SNRcnt[U]){
      int *ptmp = SNRcnt[0];
      SNRcnt[0] = SNRcnt[U];
      SNRcnt[U] = ptmp;
      double *ptmpd = SNRgmean[0];
      SNRgmean[0] = SNRgmean[U];
      SNRgmean[U] = ptmpd;
      ptmpd = lnSNRsd[0];
      lnSNRsd[0] = lnSNRsd[U];
      lnSNRsd[U] = ptmpd;
      
      double **pptmpd = SNRdist[0];
      SNRdist[0] = SNRdist[U];
      SNRdist[U] = pptmpd;
    }
    int tmpd = left[0];
    left[0] = left[U];
    left[U] = tmpd;
    tmpd = right[0];
    right[0] = right[U];
    right[U] = tmpd;

    if(remap[0] || remap[U]){
      int *ptmp = remap[0];
      remap[0] = remap[U];
      remap[U] = ptmp;
    }
    if(NickCnt[0] || NickCnt[U]){
      int *ptmp = NickCnt[0];
      NickCnt[0] = NickCnt[U];
      NickCnt[U] = ptmp;
    }

    tmp = orignumsite[0];
    orignumsite[0] = orignumsite[U];
    orignumsite[U] = tmp;

    if(origsite[0] || origsite[U]){
      double *ptmpd = origsite[0];
      origsite[0] = origsite[U];
      origsite[U] = ptmpd;
    }
    if(origsiteSD[0] || origsiteSD[U]){
      double *ptmpd = origsiteSD[0];
      origsiteSD[0] = origsiteSD[U];
      origsiteSD[U] = ptmpd;
    }
    if(origsitecov[0] || origsitecov[U]){
      float *ptmpd = origsitecov[0];
      origsitecov[0] = origsitecov[U];
      origsitecov[U] = ptmpd;
    }
    if(origsitecnt[0] || origsitecnt[U]){
      float *ptmpd = origsitecnt[0];
      origsitecnt[0] = origsitecnt[U];
      origsitecnt[U] = ptmpd;
    }
    if(origChimQuality[0] || origChimQuality[U]){
      float *ptmp = origChimQuality[0];
      origChimQuality[0] = origChimQuality[1];
      origChimQuality[U] = ptmp;
    }
    if(origChimNorm[0] || origChimNorm[U]){
      float *ptmp = origChimNorm[0];
      origChimNorm[0] = origChimNorm[1];
      origChimNorm[U] = ptmp;
    }
    if(origSegDupL[0] || origSegDupL[U]){
      float *ptmp = origSegDupL[0];
      origSegDupL[0] = origSegDupL[1];
      origSegDupL[U] = ptmp;
    }
    if(origSegDupR[0] || origSegDupR[U]){
      float *ptmp = origSegDupR[0];
      origSegDupR[0] = origSegDupR[1];
      origSegDupR[U] = ptmp;
    }
    if(origFragileEndL[0] || origFragileEndL[U]){
      float *ptmp = origFragileEndL[0];
      origFragileEndL[0] = origFragileEndL[1];
      origFragileEndL[U] = ptmp;
    }
    if(origFragileEndR[0] || origFragileEndR[U]){
      float *ptmp = origFragileEndR[0];
      origFragileEndR[0] = origFragileEndR[1];
      origFragileEndR[U] = ptmp;
    }
    if(origOutlierFrac[0] || origOutlierFrac[U]){
      float *ptmp = origOutlierFrac[0];
      origOutlierFrac[0] = origOutlierFrac[1];
      origOutlierFrac[U] = ptmp;
    }
    if(origFragSd[0] || origFragSd[U]){
      float *ptmp = origFragSd[0];
      origFragSd[0] = origFragSd[U];
      origFragSd[U] = ptmp;
    }
    if(origExpSd[0] || origExpSd[U]){
      float *ptmp = origExpSd[0];
      origExpSd[0] = origExpSd[U];
      origExpSd[U] = ptmp;
    }
    if(origFragCov[0] || origFragCov[U]){
      float *ptmp = origFragCov[0];
      origFragCov[0] = origFragCov[U];
      origFragCov[U] = ptmp;
    }
    if(origFragChiSq[0] || origFragChiSq[U]){
      float *ptmp = origFragChiSq[0];
      origFragChiSq[0] = origFragChiSq[U];
      origFragChiSq[U] = ptmp;
    }

    if(origMask[0] || origMask[U]){
      size_t *ptmp = origMask[0];
      origMask[0] = origMask[1];
      origMask[U] = ptmp;
    }

    if(origSNRcnt[0] || origSNRcnt[U]){
      if(1){
	int *ptmpd = origSNRcnt[0];
	origSNRcnt[0] = origSNRcnt[U];
	origSNRcnt[U] = ptmpd;
      }

      if(origSNRgmean[0] || origSNRgmean[U]){
	double *ptmpd = origSNRgmean[0];
	origSNRgmean[0] = origSNRgmean[U];
	origSNRgmean[U] = ptmpd;
      }
      if(origlnSNRsd[0] || origlnSNRsd[U]){
	double *ptmpd = origlnSNRsd[0];
	origlnSNRsd[0] = origlnSNRsd[U];
	origlnSNRsd[U] = ptmpd;
      }
      if(origSNRdist[0] || origSNRdist[U]){
	double **pptmpd = origSNRdist[0];
	origSNRdist[0] = origSNRdist[U];
	origSNRdist[U] = pptmpd;
      }
    }
    if(origtruesiteL[0] || origtruesiteL[U]){
      long long *ptmp = origtruesiteL[0];
      origtruesiteL[0] = origtruesiteL[U];
      origtruesiteL[U] = ptmp;
    }
    if(origtruesiteR[0] || origtruesiteR[U]){
      long long *ptmp = origtruesiteR[0];
      origtruesiteR[0] = origtruesiteR[U];
      origtruesiteR[U] = ptmp;
    }

    if(origSNR[0] || origSNR[U]){
      double *ptmpd = origSNR[0];
      origSNR[0] = origSNR[U];
      origSNR[U] = ptmpd;
    }
    if(origIntensity[0] || origIntensity[U]){
      double *ptmpd = origIntensity[0];
      origIntensity[0] = origIntensity[U];
      origIntensity[U] = ptmpd;
    }
    if(origstitch[0] || origstitch[U]){
      int *ptmp = origstitch[0];
      origstitch[0] = origstitch[U];
      origstitch[U] = ptmp;
    }

    if(origPSFWidth[0] || origPSFWidth[U]){
      double *ptmpd = origPSFWidth[0];
      origPSFWidth[0] = origPSFWidth[U];
      origPSFWidth[U] = ptmpd;
    }

    if(origImageFOV[0] || origImageFOV[U]){
      int *ptmp = origImageFOV[0];
      origImageFOV[0] = origImageFOV[U];
      origImageFOV[U] = ptmp;
    }
    if(origImageX[0] || origImageX[U]){
      double *ptmp = origImageX[0];
      origImageX[0] = origImageX[U];
      origImageX[U] = ptmp;
    }
    if(origImageY[0] || origImageY[U]){
      double *ptmp = origImageY[0];
      origImageY[0] = origImageY[U];
      origImageY[U] = ptmp;
    }

    //    if(DEBUG) assert(deltaX[0]==0 && deltaX[U]==0);
    //    if(DEBUG) assert(deltaR[0]==0 && deltaR[U]==0);
    //    if(DEBUG) assert(deltaY[0]==0 && deltaY[U]==0);
  }

  /* NOTE: deltaXfreeBlock() and deltaXinitBlock() operate over a range of maps and are multithreaded/faster */
  void deltaXfree()
  {
    for(int c = 0; c < colors; c++){
      if(deltaX[c]){
	if((USE_SSE && USE_AVX) || USE_MIC){
	  if(VERB>=2 || LEAKDEBUG>=2){
	    printf("Cmap::deltaXfree(): this=%p, c=%d, &deltaX[c][1][2]=%p, origmap=%p, deltaX[c]=%p, deltaY[c]=%p\n", this, c, &deltaX[c][1][2], origmap, deltaX[c],deltaY[c]);
	    fflush(stdout);
	  }
	  delete [] &deltaX[c][1][2];
	  delete [] deltaY[c];
	} else { 
	  delete [] deltaX[c][2];
	}
	delete [] deltaX[c];
	deltaX[c] = deltaR[c] = 0;
      }
    }
  }

  void deltaXinit()
  {
    if((VERB>=2 && id== 101000011LL) || LEAKDEBUG>=2){
      printf("deltaXinit(): mapid=%d,id=%lld\n",mapid,id);
      fflush(stdout);
    }

    for(int c = 0; c < colors; c++){
      if(deltaX[c] != 0)
	continue;
      int M = numsite[c];
      FLOAT *X = site[c];
      long long cnt = 0;
      if((USE_SSE && USE_AVX && USE_PFLOAT) || (USE_MIC && USE_PFLOAT)){
	int M_stride = COMPUTE_STRIDE4(M-1);
	PFLOAT *mem = new PFLOAT[M_stride*DELTA*3LL];
	deltaX[c] = new PFLOAT*[2*(DELTA+1)];
	deltaY[c] = new PFLOAT*[M+1];
	if(VERB>=2 || LEAKDEBUG>=2){
	  extern int numrefmaps, nummaps;
	  extern Cmap **map, **refmap;
	  if(numrefmaps > 0 && nummaps > 0)
	    printf("Cmap::deltaXinit(): this=%p,c=%d,mem=%p, origmap=%p, deltaX[c]=%p,deltaY[c]=%p (refmap[%d]=%p, map[%d]=%p)\n",
		   this, c, mem, origmap, deltaX[c], deltaY[c], numrefmaps-1,refmap[numrefmaps-1],nummaps-1,map[nummaps-1]);
	  else
	    printf("Cmap::deltaXinit(): this=%p,c=%d,mem=%p, origmap=%p, deltaX[c]=%p, deltaY[c]=%p\n", 
		   this, c, mem, origmap, deltaX[c], deltaY[c]);
	  fflush(stdout);
	}
	deltaR[c] = &deltaX[c][DELTA+1];
	for(register int n = 1; n <= DELTA; n++){
	  deltaX[c][n] = &mem[cnt-2]; cnt += M_stride;
	  for(register int J = n+1; J <= M; J++)
	    deltaX[c][n][J] = X[J] - X[J-n];
	}
	for(register int n = 1; n <= DELTA; n++){
	  deltaR[c][n] = &mem[cnt-2]; cnt += M_stride;
	  for(register int J = n+1; J <= M; J++)
	    deltaR[c][n][J] = X[M+1-J+n] - X[M+1-J];
	}

	for(register int J = 2; J <= M; J++){
	  deltaY[c][J] = &mem[cnt]; cnt += DELTA;
	  FLOAT Xj = X[J];
	  FLOAT *pX = &X[J-DELTA];
	  int s = 0;
	  if(J < DELTA+1)
	    s = DELTA+1-J;
	  for(register int n = s; n < DELTA; n++)
	    deltaY[c][J][n] = Xj - pX[n];
	}
	if(DEBUG)assert(cnt == (M-1+2*M_stride)*DELTA);
      } else {
	PFLOAT *mem = new PFLOAT[(M-1)*DELTA*2LL];
	deltaY[c] = deltaX[c] = new PFLOAT*[2LL*(M+1)];
	if(VERB>=2 || LEAKDEBUG>=2){
	  extern int numrefmaps, nummaps;
	  extern Cmap **map, **refmap;
	  if(numrefmaps > 0 && nummaps > 0)
	    printf("Cmap::deltaXinit(): this=%p,c=%d,mem=%p, origmap=%p, deltaX[c]=%p,deltaY[c]=%p (refmap[%d]=%p, map[%d]=%p)\n",
		   this, c, mem, origmap, deltaX[c], deltaY[c], numrefmaps-1,refmap[numrefmaps-1],nummaps-1,map[nummaps-1]);
	  else
	    printf("Cmap::deltaXinit(): this=%p,c=%d,mem=%p, origmap=%p, deltaX[c]=%p, deltaY[c]=%p\n", 
		   this, c, mem, origmap, deltaX[c], deltaY[c]);
	  fflush(stdout);
	}
	deltaR[c] = &deltaX[c][M+1];
	for(register int J = 2; J <= M; J++){
	  deltaX[c][J] = &mem[cnt]; cnt += DELTA;
	  FLOAT Xj = X[J];
	  FLOAT *pX = &X[J-DELTA];
	  int s = 0;
	  if(J < DELTA+1)
	    s = DELTA+1-J;
	  for(register int n = s; n < DELTA; n++)
	    deltaX[c][J][n] = Xj - pX[n];
	  deltaR[c][J] = &mem[cnt];
	  cnt += DELTA;
	  Xj = X[M+1-J];
	  pX = &X[M+1+DELTA-J];
	  for(register int n = s; n < DELTA; n++)
	    deltaR[c][J][n] = pX[-n] - Xj;
	}
	if(DEBUG)assert(cnt == (M-1)*DELTA*2LL);
      }
    }
  }

  Cmap(Ccontig *contig);/**< convert a contig into a single cmap representing the consensus map of the contig */

  Cmap(Cmap *pmap);/**< creat a copy of a single cmap : assumes .cmap input (only site[],siteSD[],sitecov[],sitecnt[] are copied)*/
  void trim(double Left, double Right);/**< trim map to include only the range Left ... Right (in kb) : assumes .cmap input (only site[],siteSD[],sitecov[],sitecnt[] are copied)*/
  void inversion(double Left, double Right); /**< Invert the map interval Left .. Right (in kb) : assumes .cmap input (only site[],siteSD[],sitecov[],sitecnt[] (ChimQuality etc if present) are copied)*/
  Cmap(Cmap *Ymap, Cmap *Xmap, int Yor, int Xor, int I, int J);/**< create pairmerge map from Ymap,Xmap in orientation Yor,Xor with crossover point at labels I,J respectively
								assumes .cmap input (only site[],siteSD[],sitecov[],sitecnt[],(ChimQuality[], if present) are updated) */
  
 private:

};

extern void maxmapalloc(int minmaps, int &maxmaps, Cmap  ** &map, int minsites, int numthreads);
extern void maxmaptruncate(int nummaps, int &maxmaps, Cmap **&map);
extern void mapcompact(int nummaps, int &maxmaps, Cmap ** &map);
extern void maxmap_free();

extern void deltaXinitBlock(Cmap **rmap, int mapstart, int mapend);
extern void deltaXfreeBlock(Cmap **rmap, int mapstart, int mapend);

extern void rawsitealloc(Cmap **Map, int start, int end);
extern void BiasCorrectSlope(double **Offset, double **Slope);
extern void BiasCorrect(Cmap **Map, int orignummaps, int nummaps, int rawalloc);
extern void BiasCorrect2(Cmap **Map, int orignummaps, int nummaps, int nthreads, int **SizeToBin, int Imaxresbias, double **Offset, double **Slope);
#endif
