#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

#include "Ccontig.h"
#include "Cmap.h"

// NOTE : unlike Ccontig.cpp this version does not include functions that depend on Cnode.h or Assembler_parameters.h and can be linked with RefAligner

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/trunk/CcontigR.cpp 5211 2016-08-02 02:46:28Z tanantharaman $");

/* allocate and compute X[c][i=0..nummaps-1][0..contig[i].numsite[c]+1] provided X[c] == 0  AND nummaps > 0 */
void Ccontig::Xmap()
{
  if(nummaps <= 0)
    return;

  for(int i= 0; i < nummaps; i++)
    if(contig[i].nummaps > 0)
      return;

  for(int c = 0; c < colors; c++){
    if(X[c])
      continue;

    X[c] = new double*[nummaps];
    for(int i = 0;i < nummaps;i++){
      int F = flip[i];
      double Scale = (scale ? scale[i] : 1.0);
      Ccontig *p = &contig[i];
      Cmap *pmap = Gmap[p->mapid];
      X[c][i] = new double[p->numsite[c]+2];
      int trimL = p->trimL[c];
      int trimR = p->trimR[c];
      if(DEBUG) assert(p->numsite[c]+1 == trimR-trimL);
      double left = pmap->site[c][trimL];
      double right = pmap->site[c][trimR];
      for(int j= 0; j <= p->numsite[c]+1; j++)
	X[c][i][j] = Scale * (F ? right - pmap->site[c][trimR - j] : pmap->site[c][trimL+j] - left);
    }
  }
}


// NOTE : sitecnt[],sitecntFN[] ... sitecntN6[] or sitemapL[] ... Hdel[] or HapSiteScore[] .. HapDeltaPhase[] are not allocated or initialized
Ccontig::Ccontig(Cmap *pmap)
{
  if(colors > 1){
    printf("Ccontig::Ccontig(Cmap): more than 1 color not supported\n");
    exit(1);
  }

  if(pmap->contig != NULL){
    printf("Attempting to converting a Haplotype Cmap (id=%lld) into a Haplotype Contig : not implemented (just use Cmap->contig instead)\n",
	   pmap->id);
    fflush(stdout);exit(1);
  }

  int M = pmap->numsite[0];
  FLOAT *X = pmap->site[0];
  if(VERB>=2){
    printf("Ccontig:converting pmap(mapid=%d,id=%lld,M=%d,Len=%0.3f) into a single map contig\n",pmap->mapid,pmap->id,M,X[M+1]);
    fflush(stdout);
  }

  /* set default initialization */
  init();

  numsite[0] = M;
  site[0] = new double[M+2];
  //  color = new char[M+1];

  for(int i = 0; i <= M; i++){
    site[0][i] = pmap->site[0][i];
    //    color[i] = 0;
  }
  site[0][M+1] = pmap->site[0][M+1];

  nummaps = 1;
  contig = new Ccontig[1];
  flip = new int[1];
  sitemap[0] = new int*[1];
  contig[0].mapid = pmap->mapid;
  if(DEBUG)assert(pmap->mapid >= 0);
  contig[0].id = pmap->id;
  contig[0].trimL[0] = 0;
  contig[0].trimR[0] = M+1;
  contig[0].numsite[0] = M;
  flip[0] = 0;
  sitemap[0][0] = new int[M+3];
  *sitemap[0][0]++ = -1;
  for(int i= 0; i <= M+1; i++)
    sitemap[0][0][i] = i;

  HapSite[0] = new int[M+2];
  HapDelta[0] = new double[M+2];
  MapPhase = new double[nummaps];

  for(int i = 0; i <= M+1; i++){
    HapSite[0][i] = 3;
    HapDelta[0][i] = 0.0;
  }
  HapSite[0][0] = HapSite[0][M+1] = 0;
  MapPhase[0] = 1.0;
  
  if(pmap->Mask[0]){
    MaskL = pmap->Mask[0][1];
    MaskR = pmap->Mask[0][pmap->numsite[0]+1];
  }
}

class CLabel {
public:
  int i; /* index in original contig1 or contig2 (normal orientation) */
  int Allele;/* 1 for unaligned label from contig1, 2 for unaligned label from contig2, 3 for aligned label (same as new HapSite value) */
  double loc;/* distance from left end of aligned interval of current Allele */
  double frac;/* fraction of aligned interval covered along current Allele */
};

/* qsort() intcmp function to sort CLabel array in increasing order of frac */
static inline int CLabelFracInc(CLabel *p1, CLabel *p2)
{
  return (p1->frac > p2->frac) ? 1 : (p1->frac < p2->frac) ? -1 : 0;
}

/* Create a new Hap contig by merging Hap contig1 and Hap contig2 based on the pairwise alignment align between the 1st map of each contig (contig1->contig[0] and contig2->contig[0])

   The resulting merged contig replaces (overwrites) contig1, while contig2 remains unchanged.

   Typically contig2->contig[0] will be completely overlapped by contig1->contig[0], so the new contig will match contig1 except where the alignment has internal outliers.
   In general, beyond the ends of the alignment between contig1 and contig2, the shorter contig is ignored and the new consensus will match the larger contig (typically contig1).

   In the region of the alignment, aligned labels correspond to homozygous labels (HapSite=3) while unaligned labels correspond to Het labels with labels from contig1 resulting in HapSite=1
   and labels from contig2 resulting in HapSite = 2. 

   It is an error if any of the labels in the aligned region of either contig1 or contig2 are already Hap labels OR any aligned interval of either contig1 or contig2 already has a Hap Indel
       return 1 and don't change contig1.

   Maps of contig2 are incorporated with their Allele reversed so that contig2->contig[0] is incorporated with the opposite Allele from contig1->contig[0] (To reverse  the Allele
   map m of contig2 means that MapPhase[m] becomes 1.0 - MapPhase[m], HapSite[i] becomes (HapSite[i]==3 ? 3 : 3-HapSite[i]), and HapDelta[i] becomes -HapDelta[i] in the merged contig1)

   contig1->contig[0] remains the 1st map of the merged contig.

   Beyond the ends of alignments between contig2->contig[0] and contig1->contig[0], the shorter contig is ignored and the consensus assumed to match the larger contig.

   0 is returned on success, 1 if the contigs could not be merged (eg due to Het labels in the aligned regions), in which case contig1 will remain unchanged.

   NOTE : sitecnt[],sitecntFN[] ... sitecntN6[] or sitemapL[] .. Hdel[] are not allocated or initialized

   HERE HERE : Generalize by allowing alignment to be between any component map m1 of contig1 and any component map m2 of contig2. Will fail (return 1) unless :
       1. Map m1 fully overlaps Map m2 and Map m2's ends (almost) match up with contig2's ends (or vice versa), OR
       2. One end of Map m1 overlaps an end of Map m2 and the overlapped ends of m1 and m2 (almost) match the corresponding ends of contig1 and contig2.

 */
int merge_contigs(Ccontig *contig1, Ccontig *contig2, Calign *align, int outlierExtend)
{
  if(colors != 1){
    printf("merge_contigs() not implemented for colors=%d\n",colors);
    fflush(stdout);exit(1);
  }
  if(!contig1->HapSite[0] || !contig2->HapSite[0] || !contig1->nummaps || !contig2->nummaps){
    printf("merge_contigs() only implemented for Haplotype contigs\n");
    fflush(stdout);exit(1);
  }
  if(contig1->MapPhase[0] < 0.999 || contig2->MapPhase[0] < 0.999){
    printf("merge_contigs(id1=%lld,id2=%lld) requires that the first map correspond to Allele A: contig1->MapPhase[0]= %0.6f,contig2->MapPhase[0]= %0.6f\n",
	   contig1->contig[0].id, contig2->contig[0].id, contig1->MapPhase[0], contig2->MapPhase[0]);
    fflush(stdout); exit(1);
  }
  if(DEBUG){/* verify that all sub-contigs are simple maps (with mapid >= 0 and nummaps = 0) */
    for(int m = 0; m < contig1->nummaps; m++){
      assert(contig1->contig[m].mapid >= 0);
      assert(contig1->contig[m].nummaps == 0);
    }
    for(int m = 0; m < contig2->nummaps; m++){
      assert(contig2->contig[m].mapid >= 0);
      assert(contig2->contig[m].nummaps == 0);
    }
  }

  int N1 = contig1->numsite[0];
  int N2 = contig2->numsite[0];
  int M1 = contig1->contig[0].numsite[0];
  int M2 = contig2->contig[0].numsite[0];
  int orientation = align->orientation;
  int U = align->numpairs;

  if(DEBUG){/* verify that all input HapSite[] values are non-zero and <= 3 */
    for(int I = 1; I <= N1; I++)
      assert(1 <= contig1->HapSite[0][I] && contig1->HapSite[0][I] <= 3);
    for(int I = 1; I <= N2; I++)
      assert(1 <= contig2->HapSite[0][I] && contig2->HapSite[0][I] <= 3);
    if(DEBUG>=1+RELEASE){
      if(!(contig1->HapSite[0][0] == 0 && contig1->HapSite[0][N1+1] == 0)){
	printf("merge_contigs(id1=%lld,id2=%lld): N1=%d, contig1->HapSite[0][0]=%d,contig1->HapSite[0][N1+1]=%d\n",
	       contig1->contig[0].id,contig2->contig[0].id, N1, contig1->HapSite[0][0],contig1->HapSite[0][N1+1]);
	fflush(stdout);
	assert(contig1->HapSite[0][0] == 0 && contig1->HapSite[0][N1+1] == 0);
      }
      assert(contig2->HapSite[0][0] == 0 && contig2->HapSite[0][N2+1] == 0);
    }
  }

  /* locate start/end of aligned regions in contig1->contig[0] */
  int iL = align->sites1[0];
  int iR = align->sites1[U-1];
  if(DEBUG) assert(0 < iL && iL <= iR && iR <= M1);

  /* locate start/end of aligned regions in contig2->contig[0] */
  int jL = align->sites2[0];// in flipped orientation of contig2
  int jR = align->sites2[U-1];// in flipped orientation of contig2
  if(DEBUG) assert(0 < jL && jL <= jR && jR <= M2);

  /* locate start/end of aligned regions in contig1 */
  int IL = contig1->sitemap[0][0][iL];  
  int IR = contig1->sitemap[0][0][iR];

  /* locate start/end of aligned regions in and contig2 */
  int JR = contig2->sitemap[0][0][orientation ? M2+1-jL : jR];// in normal orientation of contig2
  int JL = contig2->sitemap[0][0][orientation ? M2+1-jR : jL];// in normal orientation of contig2

  /* NOTE : the alignment may have to be trimmed since some of the labels from the original "larger" contig may no longer be aligned to the consensus */
  while(U > 1 && (IL < 0 || ((orientation ? JR : JL) < 0))){
    for(int t = 1; t < U; t++){
      align->sites1[t-1] = align->sites1[t];
      align->sites2[t-1] = align->sites2[t];
    }
    U--;

    iL = align->sites1[0];
    iR = align->sites1[U-1];
    if(DEBUG) assert(0 < iL && iL <= iR && iR <= M1);

    jL = align->sites2[0];// in flipped orientation of contig2
    jR = align->sites2[U-1];// in flipped orientation of contig2
    if(DEBUG) assert(0 < jL && jL <= jR && jR <= M2);

    IL = contig1->sitemap[0][0][iL];  
    IR = contig1->sitemap[0][0][iR];
    JR = contig2->sitemap[0][0][orientation ? M2+1-jL : jR];// in normal orientation of contig2
    JL = contig2->sitemap[0][0][orientation ? M2+1-jR : jL];// in normal orientation of contig2
  }

  while(U > 1 && (IR < 0 || (orientation ? JL : JR) < 0)){
    U--;

    iR = align->sites1[U-1];
    if(DEBUG) assert(0 < iL && iL <= iR && iR <= M1);

    jR = align->sites2[U-1];// in flipped orientation of contig2
    if(DEBUG) assert(0 < jL && jL <= jR && jR <= M2);

    IR = contig1->sitemap[0][0][iR];
    JR = contig2->sitemap[0][0][orientation ? M2+1-jL : jR];// in normal orientation of contig2
    JL = contig2->sitemap[0][0][orientation ? M2+1-jR : jL];// in normal orientation of contig2
  }

  if(VERB>=2 && contig1->contig[0].id == 1790 && contig2->contig[0].id==977 && orientation){
    printf("contig1->contig[0].id= %lld,contig2->contig[0].id=%lld,or=%d:\n",contig1->contig[0].id,contig2->contig[0].id,orientation);
    for(int t = 0; t < U; t++)
      printf("t=%d:sites1[t]=%d,sites2[t]=%d\n",t,align->sites1[t],align->sites2[t]);
    printf("\n");

    printf("U=%d,orientation=%d:iL=%d,iR=%d,M1=%d,jL=%d,jR=%d,M2=%d,IL=%d,IR=%d,N1=%d,JL=%d,JR=%d,N2=%d\n",
	   U,orientation,iL,iR,M1,jL,jR,M2,IL,IR,N1,JL,JR,N2);
    fflush(stdout);
  }

  if(DEBUG && !(0 < IL && IL <= IR && IR <= N1)){
    printf("U=%d,orientation=%d:iL=%d,iR=%d,M1=%d,jL=%d,jR=%d,M2=%d,IL=%d,IR=%d,N1=%d,JL=%d,JR=%d,N2=%d\n",
	   U,orientation,iL,iR,M1,jL,jR,M2,IL,IR,N1,JL,JR,N2);
    fflush(stdout);
    if(DEBUG >= 1+RELEASE) assert(0 < IL && IL <= IR && IR <= N1);
    return 1;
  }
  if(DEBUG && !(0 < JL && JL <= JR && JR <= N2)){  
    printf("U=%d,orientation=%d:iL=%d,iR=%d,M1=%d,jL=%d,jR=%d,M2=%d,IL=%d,IR=%d,N1=%d,JL=%d,JR=%d,N2=%d\n",
	   U,orientation,iL,iR,M1,jL,jR,M2,IL,IR,N1,JL,JR,N2);
    fflush(stdout);
    if(DEBUG >= 1+RELEASE) assert(0 < JL && JL <= JR && JR <= N2);  
    return 1;
  }

  /* verify that the aligned labels in both contig1 and contig2 are all homozygous labels (else return 1) */
  for(int I = IL; I <= IR; I++)
    if(contig1->HapSite[0][I] != 3 || (I < IR && contig1->HapDelta[0][I] != 0.0)){
      if(VERB/* HERE >=2 */){
	printf("merge_contigs(id=%lld,id2=%lld) failed due to label label I=%d/%d in aligned region %d..%d of contig1 being already Het (HapSite[I]=%d,HapDelta[I]=%0.4f)\n",
	       contig1->contig[0].id,contig2->contig[0].id,I,N1,IL,IR,contig1->HapSite[0][I], (I < IR ? contig1->HapDelta[0][I] : 0.0));
	fflush(stdout);
      }
      return 1;
    }

  for(int J = JL; J <= JR; J++){
    if(contig2->HapSite[0][J] != 3 || (J < JR && contig2->HapDelta[0][J] != 0.0)){
      if(VERB/* HERE >=2 */){
	printf("merge_contigs(id=%lld,id2=%lld) failed due to label label J=%d/%d in aligned region %d..%d of contig2 being already Het (HapSite[J]=%d,HapDelta[J]=%0.4f)\n",
	       contig1->contig[0].id,contig2->contig[0].id,J,N2,JL,JR,contig2->HapSite[0][J],(J < JR ? contig2->HapDelta[0][J] : 0.0));
	fflush(stdout);
      }
      return 1;
    }
  }

  /* decide which left and right ends will be copied to the result */
  int LeftContig = (contig1->site[0][IL] > (orientation ? (contig2->site[0][N2+1] - contig2->site[0][JR]) : contig2->site[0][JL])) ? 1 : 2;
  int RightContig = (contig1->site[0][N1+1] - contig1->site[0][IR] > (orientation ? contig2->site[0][JL] : contig2->site[0][N2+1] - contig2->site[0][JR])) ? 1 : 2;
  //  int LeftContig = 1;
  //  int RightContig = 1; 

  /* count number of consensus labels in new contig1 */
  int n = (LeftContig==1 ? IL : (orientation ? N2+1-JR : JL)) + (RightContig==1 ? N1+1-IR : (orientation ? JL : N2+1-JR));/* Unaligned ends including first and last aligned labels */
  n += IR-IL-1; /* adds Middle labels in contig1, includes remaining aligned labels (remain Hom) and unaligned labels (will change to Het) */

  /* add unaligned labels within alignment of contig2, which will be added as new Het sites */
  int lastJ = (orientation ? JR : JL) /* WAS JL */, J;
  for(int t = 1; t < U; t++, lastJ = J){
    int j = align->sites2[t];// in flipped orientation of contig2
    J = contig2->sitemap[0][0][orientation ? M2+1-j : j];// in normal orientation of contig2
    if(DEBUG && !(orientation ? J < lastJ : J > lastJ)){
      printf("t=%d,U=%d,orientation=%d:lastJ=%d,j=%d,M2=%d,J=%d (JL=%d,JR=%d)\n",
	     t,U,orientation,lastJ,j,M2, J, JR, JR);
      fflush(stdout);
      assert(orientation ? J < lastJ : J > lastJ);
    }
    n += (orientation ? lastJ - J : J - lastJ) - 1;
  }
  int nummaps = contig1->nummaps + contig2->nummaps;
  
  /* count total sites in all input maps */
  int totalsites = 0;
  for(int m = 0; m < contig1->nummaps; m++)
    totalsites += contig1->contig[m].numsite[0];
  for(int m = 0; m < contig2->nummaps; m++)
    totalsites += contig2->contig[m].numsite[0];

  // allocate arrays for new contig (will replace arrays in contig1[] when we are done)
  double *site = new double[n+2];
  Ccontig *contig = new Ccontig[nummaps];
  int *flip = new int[nummaps];
  int **sitemap = new int*[nummaps];
  int sitecnt = 0, mcnt = 0;
  for(int m = 0; m < contig1->nummaps; m++,mcnt++){
    sitemap[mcnt] = new int[contig1->contig[m].numsite[0]+3];
    *sitemap[mcnt]++ = -1;
    sitecnt += contig1->contig[m].numsite[0]+3;
  }
  for(int m = 0; m < contig2->nummaps; m++,mcnt++){
    sitemap[mcnt] = new int[contig2->contig[m].numsite[0]+3];
    *sitemap[mcnt]++ =-1;
    sitecnt += contig2->contig[m].numsite[0]+3;
  }
  if(DEBUG) assert(sitecnt == totalsites + 3*nummaps);
  if(DEBUG) assert(mcnt == nummaps);

  int *HapSite = new int[n+2];
  double *HapDelta = new double[n+2];
  double *MapPhase = new double[nummaps];
  HapSite[0] = HapSite[n+1] = 0;

  /* copy the component maps contig[] and orientation flip[] and Allele MapPhase[] */
  mcnt = 0;
  for(int m = 0; m < contig1->nummaps; m++, mcnt++){
    contig[mcnt] = contig1->contig[m];
    flip[mcnt] = contig1->flip[m];
    MapPhase[mcnt] = contig1->MapPhase[m];
  }
  for(int m = 0; m < contig2->nummaps; m++, mcnt++){
    contig[mcnt] = contig2->contig[m];
    flip[mcnt] = orientation ^ contig2->flip[m];
    MapPhase[mcnt] = 1.0 - contig2->MapPhase[m];
  }
  

  int *sitemap1 = new int[N1+3];/* temporary sitemap1[I=0..N1+1] maps site I of contig1 to a site on the new consensus */
  int *sitemap2 = new int[N2+3];/* temporary sitemap2[J=0..N1+1] maps site J of contig2 (in normal orientation) to a site on the new consensus */
  *sitemap1++ = -1;
  *sitemap2++ = -1;
  if(DEBUG){
    for(int I = 0; I <= N1+1; I++)
      sitemap1[I] = -1;
    for(int J = 0; J <= N2+1; J++)
      sitemap2[J] = -1;
  }

  /* next build site[0..n+1],HapSite[1..n] and HapDelta[0..n] together with temporary sitemap1[0..N1+1], sitemap2[0..N2+1] */
  site[0] = 0.0;/* left end */

  /* first copy consensus sites of Left End from (LeftContig==1 ? contig1 : contig2), updating site[], HapSite[] and HapDelta[] */
  sitecnt = 1;
  if(LeftContig==1){/* copy contig1 sites from 1 .. IL */
    sitemap1[0] = 0;
    for(int I = 1; I <= IL; I++, sitecnt++){
      sitemap1[I] = sitecnt;
      site[sitecnt] = contig1->site[0][I];
      HapSite[sitecnt] = contig1->HapSite[0][I];
      HapDelta[sitecnt - 1] = contig1->HapDelta[0][I-1];
    }
    sitemap2[orientation ? JR : JL] = sitecnt - 1;/* first aligned site */
  } else {/* copy contig2 sites from either 1 .. JL or N2 .. JR */
    if(orientation){
      sitemap2[N2+1] = 0;
      for(int J = N2; J >= JR; J--, sitecnt++){
	sitemap2[J] = sitecnt;
	site[sitecnt] = contig2->site[0][N2+1] - contig2->site[0][J];
	int H = contig2->HapSite[0][J];
	HapSite[sitecnt] = (H==3 ? 3 : 3-H);
	HapDelta[sitecnt-1] = -contig2->HapDelta[0][J];
      }
    } else {
      sitemap2[0] = 0;
      for(int J = 1; J <= JL; J++, sitecnt++){
	sitemap2[J] = sitecnt;
	site[sitecnt] = contig2->site[0][J];
	int H = contig2->HapSite[0][J];
	HapSite[sitecnt] = (H==3 ? 3 : 3-H);
	HapDelta[sitecnt - 1] = -contig2->HapDelta[0][J-1];
      }
    }
    sitemap1[IL] = sitecnt - 1;/* first aligned site */
  }
  
  /* next handle each aligned interval from left to right */
  int lasti = align->sites1[0],i;
  int lastj = align->sites2[0],j;
  int lastI = IL, I;
  lastJ = (orientation ? JR : JL)/* WAS JL */;
  CLabel *Labels = new CLabel[N1+N2+1];

  if(DEBUG) assert(align->iscore != NULL && align->outscore != NULL);

  for(int t = 1; t < U; t++, lastI = I, lastJ = J, lasti = i, lastj = j){
    i = align->sites1[t];
    j = align->sites2[t];// in flipped orientation of contig2
    int outlier = (align->iscore[t] > align->outscore[t] + (RFLOAT)0.01 || (outlierExtend && (lastj >= j || lasti >= i))) ? 1 : 0;

    I = contig1->sitemap[0][0][i];
    if(DEBUG) assert(I > lastI);
    double len1 = contig1->site[0][I] - contig1->site[0][lastI], len2;
    /* build Clabel array for all labels in the current aligned interval from (lastI..I) and (lastJ..J) (excluding the aligned labels) */
    int cnt = 0;
    for(int L = lastI + 1; L < I; L++, cnt++){
      CLabel *p = &Labels[cnt];
      p->Allele = 1;
      p->loc = contig1->site[0][L] - contig1->site[0][lastI];
      p->i = L;
      p->frac = p->loc / len1;
    }
    if(orientation){
      J = contig2->sitemap[0][0][M2+1-j];// in normal orientation of contig2    
      if(DEBUG) assert(J < lastJ);
      len2 = contig2->site[0][lastJ] - contig2->site[0][J];
      for(int L = lastJ - 1; L > J; L--, cnt++){
	CLabel *p = &Labels[cnt];
	p->Allele = 2;
	p->loc = contig2->site[0][lastJ] - contig2->site[0][L];
	p->i = L;
	p->frac = p->loc / len2;
      }
    } else {
      J = contig2->sitemap[0][0][j];// in normal orientation of contig2    
      if(DEBUG) assert(J > lastJ);
      len2 = contig2->site[0][J] - contig2->site[0][lastJ];
      for(int L = lastJ + 1; L < J; L++, cnt++){
	CLabel *p = &Labels[cnt];
	p->Allele = 2;
	p->loc = contig2->site[0][L] - contig2->site[0][lastJ];
	p->i = L;
	p->frac = p->loc / len2;
      }
    }
    qsort(Labels, cnt, sizeof(CLabel), (intcmp *)CLabelFracInc);

    double Hlen = (len1+len2)*0.5;
    double Hdelta = (len1-len2)*0.5;
    if(!outlier) Hdelta = 0.0;
    double lastF = 0.0, F;
    if(DEBUG) assert(sitecnt > 1);
    for(int L = 0; L < cnt; L++, lastF = F){
      CLabel *p = &Labels[L];
      F = p->frac;
      site[sitecnt] = site[sitecnt-1] + Hlen * (F - lastF);
      HapSite[sitecnt] = p->Allele;
      if(p->Allele == 1)
	sitemap1[p->i] = sitecnt;
      else
	sitemap2[p->i] = sitecnt;
      HapDelta[sitecnt - 1] = Hdelta * (F - lastF);
      sitecnt++;
    }
    /* add aligned label at (I,J) */
    F = 1.0;
    site[sitecnt] = site[sitecnt-1] + Hlen * (F - lastF);
    HapSite[sitecnt] = 3;
    HapDelta[sitecnt - 1] = Hdelta * (F - lastF);
    sitemap1[I] = sitecnt;
    sitemap2[J] = sitecnt;
    sitecnt++;
  }
  
  /* Finally copy consensus sites of Right End from (RightContig==1 ? contig1 : contig2), excluding last aligned site, updating site[], HapSite[] and HapDelta[] */
  double SiteStart = site[sitecnt-1];
  if(RightContig==1){/* copy contig1 sites from IR + 1 .. N1+1 */
    for(int I = IR + 1; I <= N1; I++, sitecnt++){
      if(DEBUG) assert(sitecnt <= n);
      site[sitecnt] = SiteStart + contig1->site[0][I] - contig1->site[0][IR];
      sitemap1[I] = sitecnt;
      HapSite[sitecnt] = contig1->HapSite[0][I];
      HapDelta[sitecnt - 1] = contig1->HapDelta[0][I-1];
    }
    if(DEBUG) assert(sitecnt == n+1);
    site[sitecnt] = SiteStart + contig1->site[0][N1+1] - contig1->site[0][IR];
    HapDelta[sitecnt - 1] = contig1->HapDelta[0][N1];
    sitemap1[N1+1] = sitecnt;
  } else {/* copy contig2 sites from either JR + 1 .. N2 + 1 or JL-1 ... 0 */
    if(orientation){
      for(int J = JL - 1; J > 0; J--, sitecnt++){
	if(DEBUG) assert(sitecnt <= n);
	site[sitecnt] = SiteStart + contig2->site[0][JL] - contig2->site[0][J];
	int H = contig2->HapSite[0][J];
	HapSite[sitecnt] = (H==3 ? 3 : 3-H);
	HapDelta[sitecnt-1] = -contig2->HapDelta[0][J];
      }
      if(DEBUG) assert(sitecnt == n+1);
      site[sitecnt] = SiteStart + contig2->site[0][JL];
      HapDelta[sitecnt - 1] = -contig2->HapDelta[0][0];
      sitemap2[0] = sitecnt;
    } else {
      for(int J = JR + 1; J <= N2; J++, sitecnt++){
	if(DEBUG) assert(sitecnt <= n);
	site[sitecnt] = SiteStart + contig2->site[0][J] - contig2->site[0][JR];
	int H = contig2->HapSite[0][J];
	HapSite[sitecnt] = (H==3 ? 3 : 3-H);
	HapDelta[sitecnt - 1] = -contig2->HapDelta[0][J-1];
      }
      if(DEBUG) assert(sitecnt == n+1);
      site[sitecnt] = SiteStart + contig2->site[0][N2+1] - contig2->site[0][JR];
      HapDelta[sitecnt-1] = -contig2->HapDelta[0][N2];
      sitemap2[N2+1] = sitecnt;
    }
  }

  if(VERB>=3){
    printf("After merging:id=%lld: n=%d\n",min(contig1->id, contig2->id), n);
    for(int i = 0; i <= n; i++)
      printf("i=%d:site[i]=%0.4f,HapSite[i]=%d,HapDelta[i]=%0.4f\n",i,site[i], HapSite[i], HapDelta[i]);
    fflush(stdout);
  }

  /* now remap the labels of each component map m using sitemap1[] and sitemap2[] from the original sitemap[m][] values  */
  mcnt = 0;
  for(int m = 0; m  < contig1->nummaps; m++, mcnt++){
    int cnt = contig1->contig[m].numsite[0];
    for(int j = 0; j <= cnt + 1; j++){
      int J = contig1->sitemap[0][m][j];
      if(DEBUG) assert(-1 <= J && J <= N1+1);
      if(DEBUG) assert(-1 <= sitemap1[J] && sitemap1[J] <= n+1);
      sitemap[mcnt][j] = sitemap1[J];
    }
  }
  for(int m = 0; m < contig2->nummaps; m++, mcnt++){
    int cnt = contig2->contig[m].numsite[0];
    for(int j = 0; j <= cnt + 1; j++){
      int J = contig2->sitemap[0][m][j];
      if(DEBUG) assert(-1 <= J && J <= N2+1);
      if(DEBUG) assert(-1 <= sitemap2[J] && sitemap2[J] <= n+1);
      sitemap[mcnt][j] = sitemap2[J];
    }
  }
  size_t MaskL, MaskR;
  MaskL = (LeftContig==1) ? contig1->MaskL : orientation ? contig2->MaskR : contig2->MaskL;
  MaskR = (RightContig==1) ? contig1->MaskR : orientation ? contig2->MaskL : contig2->MaskR;

  if(DEBUG>=1+RELEASE) assert(HapSite[0] == 0 && HapSite[n+1] == 0);

  /* replace contig1[0] arrays with merged contig */
  delete [] contig1->site[0];  contig1->site[0] = site;
  delete [] contig1->contig; contig1->contig = contig;
  delete [] contig1->flip; contig1->flip = flip;
  for(int i = 0; i < contig1->nummaps; i++)
    if(!contig1->blockmem || i==0)
      delete [] &contig1->sitemap[0][i][-1];
  delete [] contig1->sitemap[0]; contig1->sitemap[0] = sitemap;
  delete [] contig1->HapSite[0]; contig1->HapSite[0] = HapSite;
  delete [] contig1->HapDelta[0]; contig1->HapDelta[0] = HapDelta;
  delete [] contig1->MapPhase; contig1->MapPhase = MapPhase;
  contig1->nummaps = nummaps;
  contig1->numsite[0] = n;
  contig1->totalsites = totalsites;
  contig1->MaskL = MaskL;
  contig1->MaskR = MaskR;
  contig1->blockmem = 0;

  int m = contig1->numsite[0];
  if(DEBUG>=1+RELEASE && !(contig1->HapSite[0][0] == 0 && contig1->HapSite[0][m+1] == 0)){
    printf("m=%d,n=%d:contig1->HapSite[0][0]=%d,contig1->HapSite[0][m+1]=%d\n",m,n,contig1->HapSite[0][0],contig1->HapSite[0][m+1]);
    fflush(stdout);
    assert(contig1->HapSite[0][0] == 0 && contig1->HapSite[0][m+1] == 0);
  }

  delete [] Labels;
  delete [] &sitemap1[-1];
  delete [] &sitemap2[-1];

  return 0;
}

/* split haplotype map refmap[k]->contig into two alleles stored in refmap[k] and refmap[numrefmaps++]

  For now, only the site[0][] array of the two Alleles are computed or updated. All other arrays are deallocated (if present).
  Eventually this will have to be extended to:  siteSD[], sitecov[], sitecnt[], ChimQuality[] (etc) and orig*[]

  NOTE : Mask[0][0..m+1] information is left unchanged and index refers to hapmap index, same as contig->site[0][0..m+1] index and NOT the final refmap[k]->site[0][] index
*/
   

void split_hapmap(int k, long long newid, Cmap **&refmap, int &numrefmaps, int &maxrefmaps, int Hapnumrefmaps)
{
  if(colors==2){
    printf("split_hmap() not implemented for colors=%d\n", colors);
    fflush(stdout); exit(1);
  }

  if(DEBUG) assert(0 <= k && k < Hapnumrefmaps);
  Cmap *pmap = refmap[k];
  Ccontig *pcontig = pmap->contig;
  if(DEBUG) assert(pcontig != NULL && pmap->id == pcontig->id && pcontig->mapid == -1);

  /* locate 2nd allele at end of refmap[0..numrefmaps-1] */
  maxmapalloc(numrefmaps + 1, maxrefmaps, refmap, 0, 1);
  pmap->Allele = numrefmaps++;
  Cmap *qmap = refmap[pmap->Allele];
  qmap->Allele = k;

  /* set up 2 sub-contigs in pcontig */
  if(DEBUG) assert(pcontig->nummaps == -1);
  pcontig->nummaps = 2;
  pcontig->contig = new Ccontig[2];
  pcontig->flip = new int[2];
  pcontig->sitemap[0] = new int*[2];
  pcontig->contig[0].mapid = k;
  pcontig->contig[0].id = pmap->id;
  pcontig->contig[1].mapid = qmap->mapid = pmap->Allele;
  pcontig->contig[1].id = qmap->id = newid;
  pcontig->flip[0] = pcontig->flip[1] = 0;

  /* delete unsupported arrays in pmap and qmap  : SNRdist[], siteSD[],sitecov[],sitecnt[],ChimQuality[] (etc) and orig*[]  */
  if(pmap->SNRcnt[0]){
    if(pmap->SNRdist[0]){
      for(int i = 0; i <= pmap->numsite[0]+1; i++){
	if(pmap->blockmem >= 0) delete [] pmap->SNRdist[0][i]; 
	pmap->SNRdist[0][i] = NULL;
      }
      if(pmap->blockmem >= 0) delete [] pmap->SNRdist[0]; 
      pmap->SNRdist[0] = NULL;
    }
    if(pmap->blockmem >= 0) delete [] pmap->SNRcnt[0]; 
    pmap->SNRcnt[0] = NULL;
    if(pmap->blockmem >= 0) delete [] pmap->SNRgmean[0]; 
    pmap->SNRgmean[0] = NULL;
    if(pmap->blockmem >= 0) delete [] pmap->lnSNRsd[0]; 
    pmap->lnSNRsd[0] = NULL;
  }
  
  if(pmap->blockmem >= 0) delete [] pmap->siteSD[0];
  pmap->siteSD[0] = NULL;
  if(pmap->blockmem >= 0) delete [] pmap->sitecov[0];
  pmap->sitecov[0] = NULL;
  if(pmap->blockmem >= 0) delete [] pmap->sitecnt[0]; 
  pmap->sitecnt[0] = NULL;
  if(pmap->blockmem >= 0) delete [] pmap->ChimQuality[0]; 
  pmap->ChimQuality[0] = NULL;
  if(pmap->blockmem >= 0) delete [] pmap->ChimNorm[0]; 
  pmap->ChimNorm[0] = NULL;
  if(pmap->blockmem >= 0) delete [] pmap->SegDupL[0]; 
  pmap->SegDupL[0] = NULL;
  if(pmap->blockmem >= 0) delete [] pmap->SegDupR[0]; 
  pmap->SegDupR[0] = NULL;
  if(pmap->blockmem >= 0) delete [] pmap->FragileEndL[0]; 
  pmap->FragileEndL[0] = NULL;
  if(pmap->blockmem >= 0) delete [] pmap->FragileEndR[0]; 
  pmap->FragileEndR[0] = NULL;
  if(pmap->blockmem >= 0) delete [] pmap->OutlierFrac[0]; 
  pmap->OutlierFrac[0] = NULL;

  if(pmap->origSNRcnt[0]){
    if(pmap->origSNRdist[0]){
      for(int i = 0; i <= pmap->orignumsite[0]+1; i++)
	pmap->origSNRdist[0][i] = NULL;
      pmap->origSNRdist[0] = NULL;
    }
    pmap->origSNRcnt[0] = NULL;
    pmap->origSNRgmean[0] = NULL;
    pmap->origlnSNRsd[0] = NULL;
  }
  pmap->origsite[0] = NULL;
  pmap->origsiteSD[0] = NULL;
  pmap->origsitecov[0] = NULL;
  pmap->origsitecnt[0] = NULL;
  pmap->origChimQuality[0] = NULL;
  pmap->origChimNorm[0] = NULL;
  pmap->origSegDupL[0] = NULL;
  pmap->origSegDupR[0] = NULL;
  pmap->origFragileEndL[0] = NULL;
  pmap->origFragileEndR[0] = NULL;
  pmap->origOutlierFrac[0] = NULL;

  if(qmap->SNRcnt[0]){
    if(qmap->SNRdist[0]){
      for(int i = 0; i <= qmap->numsite[0]+1; i++){
	if(qmap->blockmem >= 0) delete [] qmap->SNRdist[0][i]; 
	qmap->SNRdist[0][i] = NULL;
      }
      if(qmap->blockmem >= 0) delete [] qmap->SNRdist[0]; 
      qmap->SNRdist[0] = NULL;
    }
    if(qmap->blockmem >= 0) delete [] qmap->SNRcnt[0]; 
    qmap->SNRcnt[0] = NULL;
    if(qmap->blockmem >= 0) delete [] qmap->SNRgmean[0]; 
    qmap->SNRgmean[0] = NULL;
    if(qmap->blockmem >= 0) delete [] qmap->lnSNRsd[0]; 
    qmap->lnSNRsd[0] = NULL;
  }
  
  if(qmap->blockmem >= 0) delete [] qmap->siteSD[0]; 
  qmap->siteSD[0] = NULL;
  if(qmap->blockmem >= 0) delete [] qmap->sitecov[0];
  qmap->sitecov[0] = NULL;
  if(qmap->blockmem >= 0) delete [] qmap->sitecnt[0]; 
  qmap->sitecnt[0] = NULL;
  if(qmap->blockmem >= 0) delete [] qmap->ChimQuality[0]; 
  qmap->ChimQuality[0] = NULL;
  if(qmap->blockmem >= 0) delete [] qmap->ChimNorm[0]; 
  qmap->ChimNorm[0] = NULL;
  if(qmap->blockmem >= 0) delete [] qmap->SegDupL[0]; 
  qmap->SegDupL[0] = NULL;
  if(qmap->blockmem >= 0) delete [] qmap->SegDupR[0]; 
  qmap->SegDupR[0] = NULL;
  if(qmap->blockmem >= 0) delete [] qmap->FragileEndL[0]; 
  qmap->FragileEndL[0] = NULL;
  if(qmap->blockmem >= 0) delete [] qmap->FragileEndR[0]; 
  qmap->FragileEndR[0] = NULL;
  if(qmap->blockmem >= 0) delete [] qmap->OutlierFrac[0]; 
  qmap->OutlierFrac[0] = NULL;

  if(qmap->origSNRcnt[0]){
    if(qmap->origSNRdist[0]){
      for(int i = 0; i <= qmap->orignumsite[0]+1; i++)
	qmap->origSNRdist[0][i] = NULL;
      qmap->origSNRdist[0] = NULL;
    }
    qmap->origSNRcnt[0] = NULL;
    qmap->origSNRgmean[0] = NULL;
    qmap->origlnSNRsd[0] = NULL;
  }
  qmap->origsite[0] = NULL;
  qmap->origsiteSD[0] = NULL;
  qmap->origsitecov[0] = NULL;
  qmap->origsitecnt[0] = NULL;
  qmap->origChimQuality[0] = NULL;
  qmap->origChimNorm[0] = NULL;
  qmap->origSegDupL[0] = NULL;
  qmap->origSegDupR[0] = NULL;
  qmap->origFragileEndL[0] = NULL;
  qmap->origFragileEndR[0] = NULL;
  qmap->origOutlierFrac[0] = NULL;

  /* retrieve Haplotype Map information from pcontig */
  int m = pcontig->numsite[0];/* total labels present in either Allele1 OR Allele2 */
  int *HapSite = pcontig->HapSite[0];
  if(DEBUG>=1+RELEASE && !(HapSite[0]==0 && HapSite[m+1]==0)){
    printf("WARNING in split_hmap:refmap[%d]->id=%lld, contig:numsite[0]=m=%d,HapSite[0]=%d,HapSite[m+1]=%d\n",
	   k,refmap[k]->id, m,HapSite[0],HapSite[m+1]);
    fflush(stdout);
    assert(HapSite[0]==0 && HapSite[m+1]==0);
    HapSite[0] = HapSite[m+1] = 0;
  }

  double *HapDelta = pcontig->HapDelta[0];
  double *site = pcontig->site[0];

  /* compute number of labels in each Allele */
  int N1 = 0, N2 = 0;

  for(int i = 1; i <= m; i++){
    if(DEBUG) assert(0 < HapSite[i]  && HapSite[i] <= 3);
    N1 += HapSite[i] & 1;
    N2 += HapSite[i] >> 1;
  }
  pcontig->contig[0].numsite[0] = N1;
  pcontig->contig[1].numsite[0] = N2;

  /* allocate sitemap and site arrays for each Allele */
  int *sitemap1 = new int[N1+3];
  int *sitemap2 = new int[N2+3];
  *sitemap1++ = -1;
  *sitemap2++ = -1;
  pcontig->sitemap[0][0] = sitemap1;
  pcontig->sitemap[0][1] = sitemap2;
  
  delete [] pmap->site[0];
  if(DEBUG) assert(qmap->site[0] == NULL);
  FLOAT *site1 = new FLOAT[N1+2];
  FLOAT *site2 = new FLOAT[N2+2];
  pmap->site[0] = site1;
  qmap->site[0] = site2;
  pmap->numsite[0] = N1;
  qmap->numsite[0] = N2;

  /* Use site[0..m+1], HapSite[1..m] and HapDelta[0..m] to compute site1[0..N1+1], site2[0..N2+1] and sitemap1[0..N1+1],sitemap2[0..N2+1] */
  if(DEBUG) assert(site[0] == 0.0);
  HapSite[0] = HapSite[m+1] = 3;// temporary
  double cum1 = 0.0, cum2 = 0.0;
  int cnt1 = 0, cnt2 = 0;
  site1[0] = site2[0] = 0.0;
  sitemap1[0] = sitemap2[0] = 0;
  for(int i = 0; i <= m; i++){
    double Hdelta = HapDelta[i];
    double frag = site[i+1] - site[i];
    if(DEBUG && !(frag >= fabs(Hdelta) - 1e-7)){
      printf("Haplotype refmap[%d](id=%lld):i=%d,m=%d,HapDelta[i]=%0.10f,site[i]=%0.10f,site[i+1]=%0.10f,frag=%0.10f(frag - fabs(HapDelta[i])= %0.10e\n",
	     k,pmap->id,i,m,HapDelta[i],site[i],site[i+1],frag, frag - fabs(Hdelta));
      fflush(stdout);
      assert(frag >= fabs(Hdelta) - 1e-7);
    }
    Hdelta = copysign(min(frag,fabs(Hdelta)), Hdelta);
    if(DEBUG && !(frag >= fabs(Hdelta))){
      printf("Haplotype refmap[%d](id=%lld):i=%d,m=%d,HapDelta[i]=%0.10f->%0.10f,site[i]=%0.10f,site[i+1]=%0.10f,frag=%0.10f(frag - fabs(HapDelta[i])= %0.10e\n",
	     k,pmap->id,i,m,HapDelta[i],Hdelta, site[i],site[i+1],frag, frag - fabs(Hdelta));
      fflush(stdout);
      assert(frag >= fabs(Hdelta));
    }
    HapDelta[i] = Hdelta;

    cum1 += frag + Hdelta;
    cum2 += frag - Hdelta;
    if(HapSite[i+1] & 1){
      site1[++cnt1] = cum1;
      sitemap1[cnt1] = i+1;
    }
    if(HapSite[i+1] & 2){
      site2[++cnt2] = cum2;
      sitemap2[cnt2] = i+1;
    }
  }
  if(DEBUG) assert(cnt1 == N1+1);
  if(DEBUG) assert(cnt2 == N2+1);
  if(DEBUG) assert(sitemap1[N1+1] == m+1);
  if(DEBUG) assert(sitemap2[N2+1] == m+1);
  HapSite[0] = HapSite[m+1] = 0;// NEW ???

  if(VERB/* HERE >=2 */){
    printf("Split Haplotype refmap[%d],id=%lld with m=%d into 2nd Allele at refmap[%d] with id=%lld: N1=%d,N2=%d,L1=%0.3f,L2=%0.3f (numrefmaps=%d->%d)\n",
	   k, pmap->id, m, pmap->Allele, qmap->id, N1,N2, pmap->site[0][N1+1], qmap->site[0][N2+1], Hapnumrefmaps, numrefmaps);
    fflush(stdout);
  }
  if(DEBUG) assert(qmap->Allele == k && qmap->contig==NULL);
  if(DEBUG>=2){
    assert(N2 == qmap->numsite[0]);
    int N = qmap->numsite[0];
    FLOAT *Y = qmap->site[0];
    /* check that site Y[0..N+1] are monotonic */
    for(int I = 0; I <= N; I++){
      if(DEBUG && !(Y[I+1] >= Y[I])){
	printf("refid=%lld:I=%d,N=%d:Y[I]=%0.8f,Y[I+1]=%0.8f\n",pmap->id,I,N,Y[I],Y[I+1]);
	fflush(stdout);
	assert(Y[I+1] >= Y[I]);
      }
    }
  }
}

// delX is similar to deltaX but used by qprobeval/mprobeval instead of pairalign

#define DELTA_BLOCK 16
static MFLOAT *delXmem[DELTA_BLOCK];
static MFLOAT **delXpmem[DELTA_BLOCK];
static int numblock = 0;

#include "parameters.h"

void Ccontig::delXinit()
{
  if(VERB>=2){
    printf("delXinit:nummaps=%d,deltaExtXRef=%d,numblock=%d (starting):wall time=%0.6f\n",nummaps, deltaExtXRef, numblock,wtime());
    fflush(stdout);
  }

  long long *fmem_start = new long long[(nummaps + 1)*2];
  long long *pmem_start = &fmem_start[nummaps+1];

  long long fmem = 0, pmem = 0 ;
  for(int m = 0; m < nummaps; m++){
    fmem_start[m] = fmem;
    pmem_start[m] = pmem;
    for(int c = 0; c < colors; c++){
      int M = contig[m].numsite[c];
      fmem += max(0,M-1) * deltaExtXRef;
      pmem += (M+1);
    }    
#if USE_AVX && !USE_AVX512 && USE_MFLOAT
    fmem += 7;
#elif (USE_MIC || USE_AVX512) && USE_MFLOAT
    fmem += 15;
#endif
  }
  fmem_start[nummaps] = fmem;
  pmem_start[nummaps] = pmem;

  for(int c = 0; c < colors; c++)
    delX[c] = new MFLOAT**[nummaps];

  if(numblock >= DELTA_BLOCK){
    printf("delXinit:numblock=%d : increase DELTA_BLOCK in CcontigR.cpp\n",numblock);
    fflush(stdout);
  }
  MFLOAT *fmemblock = delXmem[numblock] = new MFLOAT[fmem];
  MFLOAT **pmemblock = delXpmem[numblock] = new MFLOAT*[pmem];

  if(VERB>=2){
    printf("delXinit:numblock=%d,fmem=%lld,pmem=%lld: delXmem[numblock]= %p, delXpmem[numblock]= %p\n",
	   numblock,fmem,pmem,delXmem[numblock],delXpmem[numblock]);
    fflush(stdout);
  }

  int nthreads = 1;
  #ifdef _OPENMP
  nthreads = omp_get_max_threads();
  nthreads = min(nthreads,MaxThreads);
  #endif

  //  #pragma omp parallel for num_threads(nthreads) schedule(dynamic,16) if(nthreads > 1)
  for(int m = 0; m < nummaps; m++){
    long long fcnt = fmem_start[m];
    long long pcnt = pmem_start[m];

#if (USE_AVX || USE_MIC) && USE_MFLOAT
    if(VERB>=3 && m==183){
      printf("zeroing delX memory for m=%d : %lld floats from %p to just before %p\n",
	     m, fmem_start[m+1]-fcnt, &fmemblock[fcnt], &fmemblock[fmem_start[m+1]]);
      fflush(stdout);
    }
    memset(&fmemblock[fcnt],0,(fmem_start[m+1]-fcnt) * sizeof(float));
#endif

    for(int c = 0; c < colors; c++){
      int M = contig[m].numsite[c];
      FLOAT *Xm = X[c][m];

      MFLOAT **pdelXc = delX[c][m] = &pmemblock[pcnt]; pcnt += M+1;
      for(int J = 2; J <= M; J++){
        MFLOAT *pdelXcJ = pdelXc[J] = &fmemblock[fcnt]; fcnt += deltaExtXRef;
	FLOAT Xj = Xm[J];
	FLOAT *pX = &Xm[J - deltaExtXRef];
	int s = 0;
	if(J < deltaExtXRef + 1)
	  s = deltaExtXRef + 1 - J;
	for(int n = s; n < deltaExtXRef; n++)
	  pdelXcJ[n] = Xj - pX[n];
      }
    }/* c = 0 .. colors-1 */

    if(DEBUG) assert(fcnt <= fmem_start[m+1]);
    if(DEBUG) assert(pcnt <= pmem_start[m+1]);
  }/* i = mapstart .. mapend */

  numblock++;

  delete [] fmem_start;

  if(VERB>=2){
    printf("delXinit:nummaps=%d,deltaExtXRef=%d,numblock=%d (completed):wall time=%0.6f\n",nummaps,deltaExtXRef,numblock,wtime());
    fflush(stdout);
  }
}

void Ccontig :: delXfree()
{
  if(VERB>=2 || LEAKDEBUG>=2){
    printf("delXfree:nummaps=%d,deltaExtXRef=%d,numblock=%d (starting):wall time=%0.6f\n",nummaps,deltaExtXRef,numblock,wtime());
    fflush(stdout);
  }

  while(--numblock >= 0){
    delete [] delXmem[numblock];
    delete [] delXpmem[numblock];
  }
  numblock = 0;

  for(int c = 0;c < colors;c++){
    delete [] delX[c];
    delX[c] = 0;
  }

  if(VERB>=2 || LEAKDEBUG>=2){
    printf("delXfree:nummaps=%d,deltaExtXRef=%d,numblock=%d (completed):wall time=%0.6f\n",nummaps,deltaExtXRef,numblock,wtime());
    fflush(stdout);
  }
}
