#ifndef CCONTIG_H
#define CCONTIG_H

#include <stdlib.h>
#include <stdio.h>

#include "constants.h"
#include "Cnode.h"

class Cmap;/* forward declaration : cannot include "Cmap.h" until after Ccontig.h is declared */

static Ident Ccontig_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Ccontig.h 10447 2019-12-19 22:06:41Z tanantharaman $");

/** data structure to hold a Contig/Nanomap and the alignment 
   of each component Nanomap(or another Contig consensus) with the consensus map */

class Ccontig {
 public:
  int mapid; /**< If >=0, the mapid of the map of which this contig is a copy (normal orientation implied) */
  long long id;/* If mapid >=0, the id of the map of which this contig is a copy (normal orientation implied) : same as map[mapid]->id
		  If mapid < 0, the id of the consensus map of which this contig is a copy : should match refmap[i]->id for some i = 0 .. numrefmaps -1 (Not valid for Assembler) */
  int numsite[MAXCOLOR]; /**< number of sites for each color */
  int trimL[MAXCOLOR],trimR[MAXCOLOR]; /**< If mapid >= 0 the region of map[mapid] that was copied (trimL to trimR) */
  int contigid; /**< can be used to ID this contid : see Assembler_graph and refalign */

  /*  NOTE : Either (mapid>=0) OR (nummaps != 0) :
   *    If nummaps > 0 && mapid = -1 : This is a contig.
   *    If nummaps = -1 && mapid= -1 : This is a contig, but component maps are not present (this is used when reading in an Hmap until the 2 Alleles have been generated)
   *    If nummaps = 0 && mapid >= 0 : This is a Nanomap and ALL the following fields are not defined */

  int nummaps; ///< total number of maps/contigs that make up this contig (0 if this contig is a map). NOTE : Each map can typically only occur once, unless -BestRefWt or -MultiMatch is used
  int totalsites; /**< tmp : sum of numsites for all component contigs (OR -1) */
  double left,right; /**< original consensus ends (if extend > 0) */
  double Y1,YN;/* original first and last label */
  double Lfrozen, Rfrozen;/* original location of site[0][Lfrozen], site[0][Rfrozen] : used with SegDupMask to keep track of change in location of left,right during refinement, since change final location of 
			     site[0][Lfrozen],site[0][Rfrozen] is known */
  double *site[MAXCOLOR]; /**< site[c][i=0..numsite[c]+1] is the location of site i for color c */

  float *sitecnt[MAXCOLOR]; ///< sitecnt[c][i=1..numsite] is a count of the number of aligned maps aligned to site i (If sitecnt[i]==0, then site[i] is the end of a map and NOT an actual site) 
  float *sitecntFN[MAXCOLOR]; /**< sitecntFN[i=0..numsite+1] is a count of the number of aligned maps spanning site i, but NOT necessarily aligned to site i 
			         NOTE: In refine(), sitecntFN[] is instead used to store ChimQuality */
  float *fragcnt[MAXCOLOR]; /**< fragcnt[i=0..numsite+1] is the count of the number of aligned map spanning the contig site interval i .. i+1 (excludes the unaligned ends of each map) 
			       NOTE: In refine(), unless FRAGCOV_FIX, fragcnt[] instead is the number of aligned maps (excluding unaligned ends) crossing each site i */

  float *fragcntT[MAXCOLOR]; ///< Only computed by refine() : same as fragcnt[] but includes the unaligned ends of each map up to the last ref site overlapped
  float *fragcntTnorm[MAXCOLOR];///< Only computeed by refine() if -CovNorm : average map size to generate map size normalized version of fragcntT[] 
  float *sitecntFNnorm[MAXCOLOR];///< Only computed by refine() if -TrimNorm : map population count to generate sitecntFN[] 

  float *sitecntN2[MAXCOLOR];/* Only computed by refine() : sitecntN2[0][k] is count of molecules that aligned on right side but overlap on left side of Hcuts[k] */
  float *sitecntN3[MAXCOLOR];/* Only computed by refine() : sitecntN3[0][k] is count of molecules that aligned on left side but overlap on right side of Hcuts[k] */
  float *sitecntN4[MAXCOLOR];/* Only computed by refine() : sitecntN4[0][k] is count of molecules that aligned on right side of Hcuts[k] but end within FragileQualityEnd (kb) of Hcuts[k] */
  float *sitecntN5[MAXCOLOR];/* Only computed by refine() : sitecntN5[0][k] is count of molecules that aligned on left side of Hcuts[k] but end within FragileQualityEnd (kb) of Hcuts[k] */

  float *sitecntN6[MAXCOLOR];/* Only computed by refine() : sitecntN6[0][k] is fraction of molecules that aligned with an outlier across interval Hcuts[k]..Hcuts[k+1] (only valid for Hdel[k] == 0) */

  float *fragSd[MAXCOLOR];/* Only computed by refine() : sample sd for non-outlier aligned intervals across Hcuts[k]..Hcuts[k+1] (only valid for Hdel[k] == 0) */
  float *expSd[MAXCOLOR];/* Only computed by refine() : Expected sample sd for non-outlier aligned intervals across Hcuts[k]..Hcuts[k+1] (only valid for Hdel[k] == 0) */
  float *fragBias[MAXCOLOR];/* Only computed by refine() : sample bias for non-outlier aligned intervals across Hcuts[k] .. Hcuts[k+1] (only valid for Hdel[k] == 0) */
  float *fragCov[MAXCOLOR];/* Only computed by refine() : sample coverage for non-outlier aligned intervals across Hcuts[k]..Hcuts[k+1] (only valid for Hdel[k] == 0) */
  float *fragChiSq[MAXCOLOR];/* Only computed by refine() : ChiSq statistic for fragSd relative to expected value for Hcuts[k]..Hcuts[k+1] (only valid for Hdel[k] == 0) */

  float *sitecntN1[MAXCOLOR][MAX_CHIM];/* alternate values of ChimQuality : sitecntN1[c][0][i] == sitecntFN[c][i] */

  Ccontig *contig; ///< contig[i=0..nummaps-1] is the i'th component of this contig (either a map or a contig) 
  int *flip; ///< flip[i=0..nummaps-1] is the orientation of contig[i] relative to the consensus map 
  int **sitemap[MAXCOLOR]; /* sitemap[c][i=0..nummaps-1][j=0..contig[i].numsite[c]+1] maps site j of color c of map i (in orientation flip[i]) to a site of this contig (same color) :
			      all values lie in the range 0 ... numsite[c]+1, where 0 refers to the contig's left end and numsite[c]+1 to the contig's right end (value can be -1 if
			      a site is not present/aligned in the consensus). For convenience sitemap[c][i][-1] always returns -1.
			    */
  int **sitemapL[MAXCOLOR]; /* sitemapL[c][i=0..nummaps-1][j=0..contig[i].numsite[c]+1] is same as sitemap[c][i][j] except when a map site j maps to multiple consensus sites that failed to resolve :
			       in that case sitemapL[c][i][j] points to the leftmost (smallest) consensus site index, while sitemap[c][i][j] points to the rightmost (largest) consensus site index */
  double *mapWT; /**< pre-refinement : mapWT[i=0..nummaps-1] is input Map weight (typically 1.0 unless -Erefine or -BestRefWT was used) 
                  * post-refinement : mapWT[i=0..nummaps-1] is Map weight (based on -LRbias -TB and -BestRefWT) used for QualityScore and coverage computation
		  */
  double *origmapWT; /* The pre-refinement mapWT[] value based on -BestRefWT only (excluding -Erefine) to be used to compute the post-refinement mapWT[] values */
  double *maxoutlier;/* pre-refinement : maxoutlier[i=0..nummaps-1] is maxoutlier value of original alignment for matchgroup/map i as computed by refalign */

  double *scale; ///< scale[i=0..nummaps] is the Scaling factor applied to contig[i] to get X[i] 
  double **X[MAXCOLOR]; ///< If X[c] != 0 : X[c][m=0..nummaps-1][j=0..contig[i].numsite[c]+1] are the Nanomap sites of map i for color c in orientation flip[i] 

  MFLOAT ***delX[MAXCOLOR];/* Used in qprobeval.cpp & mprobeval.cpp :
	     delX[c][m=0..nummaps-1][J=2..M][n = max(0,deltaExtXRef+1-J) .. deltaExtXRef-1] = X[c][m][J] - X[c][m][J - deltaExtXRef + n], where M == contig[m].numsite[c] */

  int blockmem;/* If X[][] and sitemap[][] are each allocated as a single contiguous block (so only sitemap[c][0] and X[c][0] need be deleted) */
  Calign **align;/* If != 0 : align[i=0..nummaps-1] is the pointer to the original alignment for matchgroup/map i as computed by refalign 
		   (currently only used in refalign to compute maxoutlier[] and mapWT[]) */
  
  /* The following are only valid (and allocated) post-refinement : actual sites are now identified by sitecnt[i] > 0 */
  double *LP; ///< LP[i=0..nummaps-1] : post-refinement log-likelihood score (always assumes biaswt=0) 
  double *logPV; ///< logPV[i=0..nummaps-1] : post-refinement log10(Pvalue) of each alignment 
  int **outlier[MAXCOLOR]; /**< outlier[c][i=0..nummaps-1][j=0..contig[i].numsite[c]] is 1 IFF the interval (j..j+1) of map i, color c, is part of an alignment outlier region
			      outlier[c][0] & 2 is used to flag left endoutlier and outlier[c][numsite[c]] & 2 is used to flag right endoutlier */
  int *Hdel[MAXCOLOR];///< Hdel[c][i=1..numsite] == 0 IFF site[c][i] is present (computed during refinement) */

  /* The following are only valid (and allocated) if HapSite[0] != NULL. 
     This only happens post-refinement if -Haplotype was specified or for Hmaps created by pairmerge or read in from an .hmap. If so:

     1. Each component map m = 0.. nummaps-1 (if nummaps > 0) must correspond to a specific Allele, identified by MapPhase[m]
   */
  int *HapSite[MAXCOLOR];   /**< HapSite[c][i=1..numsite[c]] is the haplotype status of site[c][i] for color c :
			       0 : site[c][i] is not present (always true IFF sitecnt[i]==0)
			       1 : site[c][i] is present in Allele A only (implies sitecnt[c][i] > 0)
			       2 : site[c][i] is present in Allele B only (implies sitecnt[c][i] > 0)
			       3 : site[c][i] is present for both alleles (implies sitecnt[c][i] > 0)
			    */
  double *HapDelta[MAXCOLOR];/**< HapDelta[c][i=0..numsite[c]] is the haplotype delta value for interval between site[i] and site[i+1] for color c : 
				==0 : interval has no indel : This is required to be the case if Hapsite[i] == 0 to avoid duplicate representations.
				> 0 : A heterozygous indel overlaps the interval site[i..i+1] and Allele A is larger by HapDelta[i] while Allele B is smaller by HapDelta[i]
				< 0 : A heterizygous indel overlaps the interval site[i..i+1] and Allele A is smaller by -HapDelta[i] while Allele B is larger by -HapDelta[i].
				NOTE : large indels may be split into several pieces spanning multiple intervals : they can be recognized since the 
				internal intervals will be (approximately) completely deleted in one Allele (ie |HapDelta[i]| ~ site[i+1] - site[i]) 
				along the internal sites (see HapSite[i] & HapSite[i+1]) for that same Allele. */
  double *HapSiteScore[MAXCOLOR];/**< HapSiteScore[c][i=1..numsite[c]] is only defined if HapSite[c][i] == 1 or 2 : It is the log-likelihood score improvement due adding a Haplotype Site */
  double *HapDeltaScore[MAXCOLOR];/**< HapDeltaScore[c][i=0..numsite[c]] is only defined if HapDelta[c][i] != 0.0 : It is the log-likelihood score improvement due to the Haplotype Indel */
                                  /* NOTE : HapSiteScore[] and HapDeltaScore[] are natural log here but log10 in the .hmap or .cmap files and must be converted during file IO */
  double *HapSitePhase[MAXCOLOR];/**< HapSitePhase[c][i=1..numsite[c]] is only defined if HapSite[c][i] == 1 or 2 :
				    It is the LP drop if all haplotype differences starting at HapSite[c][i] (including HapDelta[c][i]) and to the right are reversed in orientation : 
				    small or 0 value signals start of phase contig */
  double *HapDeltaPhase[MAXCOLOR];/**< HapDeltahase[c][i=0..numsite[c]] is only defined if HapDelta[c][i] != 0.0 :
				     It is the LP drop if all haplotype differences starting at HapDelta[c][i] (excluding HapSite[c][i]) and to the right are reversed in orientation : 
				     small or 0 value signals start of phase contig */				     

  double *MapPhase;/**< MapPhase[m=0..MD-1] is the "probability that map m belongs to Allele A" ( same as 1 - "probability map m belongs to Allele B") */


  /* Mask information to be passed through refinement */
  size_t *Mask[MAXCOLOR];/* If Mask[c] != 0 : The full Mask information Mask[c][i=0..numsite[c]] : typically only present for Haplotype contigs */

  /* NOTE : The index values of Mask[c][i] are misleading after refinement since labels may have been added/deleted/shifted hence the following summary is generated before refinement and used
    after refinement to regenerate reasonable Mask values */

  size_t MaskL,MaskR;/* end Mask values */
  int SegDupCnt;/* number of SegDup regions (with Mask & SegDupMask) */
  int SegDupMax;/* allocation size of SegDupStart[], SegDupEnd[] */
  double *SegDupStart;/* SegDupStart[0..SegDupCnt-1] is start of SegDup regions as a fraction of total contig length (including extension regions, before trimming) */
  double *SegDupEnd;/* SegDupEnd[0..SegDupCnt-1] is end of SegDup regions as a fraction of total contig length (including extension regions, before trimming) */

  //  int HapContigs;///< number Haplotype Phase contigs (Not Yet Implemented) */
  //  int *HapFragId[MAXCOLOR];/*< HapFragId[c][i=0..numsite] identifies the Haplotype Phase contig (0..HapContigs-1) for interval site[c][i..i+1] (IFF HapDelta[c][i] != 0.0) */
  //  int *HapSiteId[MAXCOLOR];/*< HapSiteId[c][i=1..numsite] identifies the Haplotype Phase contig (0..HapContigs-1) for site[c][c] (IFF HapSite[c][i] == 1 or 2) */
  
  // Ccontig *Haplotype;/* Haplotype[0..1] represent the two Haplotypes of the current contig as independent contigs (Phasing is arbitrary between Phase Contigs) */
      

  Ccontig(){
    init();
  };
  Ccontig(Cnode *Nnode);/**< convert a graph node to a contig with a single map (nested 1 deep, with nummaps=1) */
  Ccontig(Cmap *pmap);/* convert a Cmap to a Hap contig with a single map (nexted 1 deep, with nummaps = 1), where every label is Homozygous (HapSite[1..numsite] = 3) and the single map is Allele A
			 (MapPhase[0] = 1.0)
			 If the Cmap is already an Hmap (Cmap->contig != 0), this is an error : this function is used to initialize Cmap->contig for a regular Cmap with Cmap->contig=0
			 NOTE : sitecnt[],sitecntFN[] ... sitecntN6[] or sitemapL[] ... Hdel[] or HapSiteScore[] .. HapDeltaPhase[] are not allocated or initialized
		      */

  ~Ccontig(){
    allfree();
  }
  void init(){
    if(VERB>=2){
      printf("Initializing empty Ccontig: this= %p\n",this);
      fflush(stdout);
    }
    blockmem = 0;
    mapid = -1;
    nummaps = 0;
    flip = NULL;
    contig = NULL;
    scale = NULL;
    LP = logPV = mapWT = origmapWT = NULL;
    align = NULL;
    maxoutlier = NULL;
    //    Haplotype = 0;
    MapPhase = NULL;
    MaskL = MaskR = 0;
    SegDupCnt = SegDupMax = 0;
    SegDupStart = SegDupEnd = NULL;
    for(int c = colors; --c >= 0;){
      Mask[c] = NULL;

      numsite[0] = 0;
      X[c] = 0;
      delX[c] = 0;
      site[c] = 0;
      sitecnt[c] = sitecntFN[c] = fragcnt[c] = 0;
      sitecntN1[c][0] = 0;

      fragcntT[c] = fragcntTnorm[c] = sitecntFNnorm[c] = 0;
      sitecntN6[c] = sitecntN5[c] = sitecntN2[c] = sitecntN3[c] = sitecntN4[c] = 0;
      fragSd[c] = expSd[c] = fragBias[c] = fragCov[c] = fragChiSq[c] = 0;
      sitemap[c] = sitemapL[c] = 0;
      outlier[c] = 0;
      Hdel[c] = 0;

      HapDelta[c] = 0;
      HapSite[c] = 0;
      HapSiteScore[c] = 0;
      HapDeltaScore[c] = 0;
      HapSitePhase[c] = 0;
      HapDeltaPhase[c] = 0;
      //HapFragId[c] = HapSiteId[c] = 0;
    }
  }

  void allfree(){
    if(VERB>=2){
      printf("Ccontig::allfree():mapid=%d,id=%lld,nummaps=%d,colors=%d,blockmem=%d\n",mapid,id,nummaps,colors,blockmem);
      fflush(stdout);
    }
    if(mapid < 0) {
      if(nummaps){
	delete [] flip;
	if(DEBUG>=1+RELEASE && nummaps > 0) assert(contig != NULL);
	for(int i= 0; i < nummaps; i++){
	  if(contig != NULL)
	    contig[i].allfree();
	  for(int c = colors; --c >= 0;){
	    if(VERB>=2){
	      printf("\t i=%d/%d,c=%d:sitemap[c]=%p,sitemap[c][i]=%p\n",
		     i,nummaps,c,sitemap[c],sitemap[c] ? sitemap[c][i] : NULL);
	      fflush(stdout);
	    }
	    if(Mask[c]){
	      delete [] Mask[c];
	      Mask[c] = NULL;
	    }
	    if(sitemap[c] && sitemap[c][i]){
	      if(!blockmem || i==0)
		delete [] &sitemap[c][i][-1];
	      sitemap[c][i] = 0;
	    }
	    if(sitemapL[c] && sitemapL[c][i]){
	      delete [] &sitemapL[c][i][-1];
	      sitemapL[c][i] = 0;
	    }
	    if(outlier[c] && outlier[c][i]){
	      delete [] outlier[c][i];
	      outlier[c][i] = 0;
	    }
	  }
	}
	for(int c = colors; --c >= 0;){	
	  delete [] sitemap[c];
	  if(sitemapL[c]){
	    delete [] sitemapL[c];
	    sitemapL[c] = 0;
	  }
	  if(outlier[c]){
	    delete [] outlier[c];
	    outlier[c] = 0;
	  }
	  if(Hdel[c]){
	    delete [] Hdel[c];
	    Hdel[c] = 0;
	  }
	  if(HapDelta[c]){
	    delete [] HapDelta[c];
	    HapDelta[c] = 0;
	  }
	  if(HapSite[c]){
	    delete [] HapSite[c];   HapSite[c] = 0;
	    delete [] HapSiteScore[c]; HapSiteScore[c] = 0;
	    delete [] HapDeltaScore[c]; HapDeltaScore[c]= 0;
	    delete [] HapSitePhase[c]; HapSitePhase[c] = 0;
	    delete [] HapDeltaPhase[c]; HapDeltaPhase[c] = 0;
	  }
	  
#if 0
	  if(HapFragId[c]){
	    delete [] HapFragId[c];
	    HapFragId[c] = 0;
	  }
	  if(HapSiteId[c]){
	    delete [] HapSiteId[c];
	    HapSiteId[c] = 0;
	  }
#endif
	}

	delete [] contig;
	if(scale){
	  delete [] scale;
	  scale = 0;
	}
	for(int c = 0; c < colors; c++){
	  if(X[c]){
	    if(!blockmem)
	      for(register int i = 0; i < nummaps; i++)
		delete [] X[c][i];
	    else if(nummaps > 0)
	      delete [] X[c][0];
	    delete [] X[c];
	    X[c] = 0;
	  }
	}
	if(LP){
	  delete [] LP;
	  LP = 0;
	}
	if(logPV){
	  delete [] logPV;
	  logPV = 0;
	}
	if(mapWT){
	  delete [] mapWT;
	  mapWT = NULL;
	}
	if(origmapWT){
	  delete [] origmapWT;
	  origmapWT = NULL;
	}
	if(MapPhase){
	  delete [] MapPhase;
	  MapPhase = NULL;
	}
	if(SegDupStart){
	  delete [] SegDupStart;
	  SegDupStart = NULL;
	}
	if(SegDupEnd){
	  delete [] SegDupEnd;
	  SegDupEnd = NULL;
	}
	SegDupCnt = 0;
	if(align){
	  delete [] align;
	  align = NULL;
	}
	if(maxoutlier){
	  delete [] maxoutlier;
	  align = NULL;
	}
	nummaps = 0;
      }
      for(int c = colors; --c >= 0;){	
	if(numsite[c] > 0){
	  delete [] site[c];
	  //	delete [] color;
	  delete [] sitecnt[c];
	  delete [] sitecntFN[c];
	  delete [] sitecntN1[c][0];
	  delete [] fragcnt[c];
	  if(sitecntFNnorm[c]){
	    delete [] sitecntFNnorm[c];
	    sitecntFNnorm[c] = 0;
	  }
	  if(sitecntN2[c]){
	    delete [] sitecntN2[c];
	    sitecntN2[c] = 0;
	  }
	  if(sitecntN3[c]){
	    delete [] sitecntN3[c];
	    sitecntN3[c] = 0;
	  }
	  if(sitecntN4[c]){
	    delete [] sitecntN4[c];
	    sitecntN4[c] = 0;
	  }
	  if(sitecntN5[c]){
	    delete [] sitecntN5[c];
	    sitecntN5[c] = 0;
	  }
	  if(sitecntN6[c]){
	    delete [] sitecntN6[c];
	    sitecntN6[c] = 0;
	  }
	  if(fragSd[c]){
	    delete [] fragSd[c];
	    fragSd[c] = 0;
	  }
	  if(expSd[c]){
	    delete [] expSd[c];
	    expSd[c] = 0;
	  }
	  if(fragBias[c]){
	    delete [] fragBias[c];
	    fragBias[c] = 0;
	  }
	  if(fragCov[c]){
	    delete [] fragCov[c];
	    fragCov[c] = 0;
	  }
	  if(fragChiSq[c]){
	    delete [] fragChiSq[c];
	    fragChiSq[c] = 0;
	  }

	  if(fragcntT[c]){
	    delete [] fragcntT[c];
	    fragcntT[c] = 0;
	  }
	  if(fragcntTnorm[c]){
	    delete [] fragcntTnorm[c];
	    fragcntTnorm[c] = 0;
	  }

	  numsite[c] = 0;
	}
      }
    }
  };

  int totsites() {
    if(mapid >= 0){
      totalsites = 0;
      for(register int c = colors; --c >= 0;)
	totalsites += numsite[c];
      return totalsites;
    }
    if(totalsites >= 0)
      return totalsites;
    totalsites = 0;
    for(register int i = nummaps; --i >= 0;)
      totalsites += contig[i].totsites();
    return totalsites;
  }

  void contigflip();/**< flip the entire contig (along with all subcontigs) */

  /** create a new contig by merging contig1 and contig2 using alignment between one map of contig1 and another
     map of contig2 provided by align[]. Creates a contig with 2 components: contig1 and contig2 
     NOTE that the alignment between contig1 and contig will be incomplete, since only sites from the two aligned maps are aligned.
     The unaligned sites (of Contig1 or Contig2) can either be interpolated into the combined contig (with sites closer than
     resolution "res" are combined) OR left unaligned. Either way this may produce FP and FN errors in the consensus,
     that must be fixed later by Refinement.

     contig1 and contig2 are incorporated into the new contig and the original contigs are left empty with numsite and nummaps set to 0.

     For efficiency the merging should be done in balanced (binary tree) fashion for a G*log(G) run time. Otherwise
     if the merging is done left to right the total run time will be G^2, where G is the contig size in sites.
  */

#if CALIGN_SMALL
  Ccontig(Ccontig *contig1, Ccontig *contig2, CalignA *align, int rverb);
#else
  Ccontig(Ccontig *contig1, Ccontig *contig2, Calign *align, int rverb);
#endif

  void flatten();/**< flatten the current contig so all components in contig[] are maps (nexted 1-deep) */
  int findmap(int id) {/** return index to component map/contig with given mapid, -1 if no match (Does NOT check if current mapid == id) */
    if(DEBUG) assert(mapid < 0);
    for(register int i=nummaps;--i>=0;){
      if(contig[i].mapid == id)
	return i;
    }
    return -1;
  }
  
  void Xmap();/**< allocate and compute X[c][i=0..nummaps-1][0..contig[i].numsite[c]+1] provided X[c] == 0  AND nummaps > 0 */
  void display();

  void delXinit();
  void delXfree();
};

extern void output_draft(Ccontig *contig, long long contigid, char *prefix, int Lfrozen, int Rfrozen, int Allele);
extern void output_contigs(Ccontig *contig, int numcontigs, char *prefix);
extern void input_contigs(Ccontig * &contig, int &numcontigs, char *prefix);

extern int merge_contigs(Ccontig *contig1, Ccontig *contig2, Calign *palign, int outlierExtend);
// NOTE : sitecnt[],sitecntFN[] ... sitecntN6[] or fragSd[],expSd[],fragBias[],fragCov[],fragChiSq[] or sitemapL[] ... Hdel[] or HapSiteScore[] .. HapDeltaPhase[] are not allocated or initialized

extern void split_hapmap(int k, long long newid, Cmap **& refmap, int &refmaps, int &maxrefmaps, int Hapnumrefmaps);

#endif
