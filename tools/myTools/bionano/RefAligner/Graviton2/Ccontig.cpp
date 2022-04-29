#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

#include "Ccontig.h"
#include "Cmap.h"
#include "Assembler_parameters.h"

extern Cnode *node;
extern int *nodeorder;

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Ccontig.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

#define INTERLEAVE 1 /* 0 : Keep the longer unaligned end, leave the other one unaligned \
			1 : Interleave both unaligned ends */

Ccontig::Ccontig(Cnode *Node)
{
  if(colors > 1){
    printf("Ccontig::Ccontig(Cnode): more than 1 color not supported\n");
    exit(1);
  }

  if(VERB>=2){
    printf("Ccontig:converting nodeid=%d,mapid=%d to contig\n",Node->nodeid,Node->mapid);
    fflush(stdout);
  }

  /* set default initialization */
  init();
  /*  numsite = 0;
  mapid = -1;
  nummaps = 0;
  X = 0;
  scale = 0;
  LP = logPV = mapWT = 0;*/

  register Cmap *pmap = Gmap[Node->mapid];
  register int M = pmap->numsite[0];
  register int left = Node->trimL;
  register int right = Node->trimR;
  if(DEBUG)assert(left >= 0 && right <= M+1);

  numsite[0] = right - left - 1;
  site[0] = new double[numsite[0]+2];
  //  color = new char[numsite[0]+1];
  sitecnt[0] = new float[numsite[0]+1];
  sitecntFN[0] = new float[numsite[0]+2];
  fragcnt[0] = new float[numsite[0]+2];
  fragcntT[0] = new float[numsite[0]+2];
  if(CovNorm)
    fragcntTnorm[0] = new float[numsite[0]+2];
  if(TrimNorm >= 0){
    sitecntFNnorm[0] = new float[numsite[0]+2];
    sitecntN2[0] = new float[numsite[0]+2];
    sitecntN3[0] = new float[numsite[0]+2];
    sitecntN4[0] = new float[numsite[0]+2];
    sitecntN5[0] = new float[numsite[0]+2];
    sitecntN6[0] = new float[numsite[0]+2];

    fragSd[0] = new float[numsite[0]+2];
    expSd[0] = new float[numsite[0]+2];
    fragBias[0] = new float[numsite[0]+2];
    fragCov[0] = new float[numsite[0]+2];
    fragChiSq[0] = new float[numsite[0]+2];
  }
  site[0][0] = 0.0;
  fragcnt[0][0] = 1;
  for(register int i = 1; i <= numsite[0]; i++){
    site[0][i] = pmap->site[0][i+left] - pmap->site[0][left];
    //    color[i] = 0;
    sitecnt[0][i] = 1;
    sitecntFN[0][i] = 0;
    fragcnt[0][i] = 1;
  }
  site[0][numsite[0]+1] = pmap->site[0][numsite[0]+1+left] - pmap->site[0][left];

  nummaps = 1;
  contig = new Ccontig[1];
  flip = new int[1];
  sitemap[0] = new int*[1];
  contig[0].mapid = Node->mapid;
  contig[0].trimL[0] = Node->trimL;
  contig[0].trimR[0] = Node->trimR;
  if(DEBUG)assert(pmap->mapid >= 0);
  contig[0].numsite[0] = numsite[0];
  flip[0] = 0;
  sitemap[0][0] = new int[numsite[0]+3];
  *sitemap[0][0]++ = -1;
  for(register int i= 0; i <= numsite[0]+1; i++)
    sitemap[0][0][i] = i;
}

class CAlignOrder {
public:
  int leftsite;
  int m;
  Ccontig *contig;
};

static int AlignedSiteInc(CAlignOrder *p1, CAlignOrder *p2)
{
  return (p1->leftsite > p2->leftsite) ? 1 : (p1->leftsite < p2->leftsite) ? -1 : 0;
}

void Ccontig::display()
{
  if(mapid >= 0)
    printf("mapid=%d,trimL=%d,trimR=%d,id=%lld,numsite=%d\n",mapid,trimL[0],trimR[0],Gmap[mapid]->id,numsite[0]);
  else {
    printf("numsite=%d:\n",numsite[0]);
    for(register int t = 0; t <= numsite[0]+1;t++){
      if(t<=0)
	printf(" t=0:site[t]=%0.3f\n",site[0][0]);
      else if(t > numsite[0])
	printf(" t=%d:site[t]=%0.3f,fragcnt[t-1]=%0.2f\n",t,site[0][t],fragcnt[0][t-1]);
      else
	printf(" t=%d:site[t]=%0.3f,fragcnt[t-1]=%0.2f,sitecnt[t]=%0.2f,sitecntFN[t]=%0.2f\n",t,site[0][t],fragcnt[0][t-1],sitecnt[0][t],sitecntFN[0][t]);
    }
    fflush(stdout);

    /* Sort maps by leftmost aligned site */
    CAlignOrder *maporder = new CAlignOrder[nummaps];
    for(int m = 0; m < nummaps; m++){
      maporder[m].m = m;
      maporder[m].leftsite = sitemap[0][m][0];
      maporder[m].contig = &contig[m];
    }
    if(DEBUG/* HERE >=2 */){
      for(int t = 0; t < nummaps; t++){
	int m = maporder[t].m;
	if(DEBUG && !(maporder[t].contig == &contig[m])){
	  printf("Before qsort: t=%d,m=%d, maporder[t].contig - contig = %ld, maporder[t].contig = %p, &contig[m]= %p\n",
		 t,m,maporder[t].contig - &contig[0],maporder[t].contig,&contig[m]);
	  fflush(stdout);
	  assert(maporder[t].contig == &contig[m]);
	}
      }
    }
    qsort(maporder,nummaps,sizeof(CAlignOrder),(intcmp*)AlignedSiteInc);

    if(DEBUG/* HERE >=2 */){
      for(int t = 0; t < nummaps; t++){
	int m = maporder[t].m;
	if(DEBUG && !(maporder[t].contig == &contig[m])){
	  printf("After qsort: t=%d,m=%d, maporder[t].contig - contig = %ld, maporder[t].contig = %p, &contig[m]= %p\n",
		 t,m,maporder[t].contig - &contig[0],maporder[t].contig,&contig[m]);
	  fflush(stdout);
	  for(int k = 0; k < nummaps; k++)
	    printf("    t=%d,m=%d, maporder[t].contig - contig = %ld, maporder[t].contig = %p, &contig[m]= %p\n",
		   k,maporder[k].m,maporder[k].contig - &contig[0],maporder[k].contig, &contig[maporder[k].m]);
	    
	  assert(maporder[t].contig == &contig[m]);
	}
      }
    }

    for(int t = 0; t < nummaps; t++){
      int m = maporder[t].m;
      if(DEBUG) assert(maporder[t].contig == &contig[m]);
      if(contig[m].mapid>=0){
	register Cmap *pmap = Gmap[contig[m].mapid]; 
	printf("  map=%d(mapid=%d,id=%lld,trimL=%d,trimR=%d),flip=%d,numsite=%d:\n",
	       m,contig[m].mapid,Gmap[contig[m].mapid]->id,contig[m].trimL[0],contig[m].trimR[0],flip[m],contig[m].numsite[0]);
	if(!flip[m])
	  for(register int k = 0; k <= contig[m].numsite[0]+1; k++)
	    printf("    k=%d:site[%d]=%0.3f,sitemap[k]=%d\n",k,k+contig[m].trimL[0],pmap->site[0][k+contig[m].trimL[0]],sitemap[0][m][k]);
	else
	  for(register int k = 0; k <= contig[m].numsite[0]+1; k++)
	    printf("    k=%d:site[%d]=%0.3f,sitemap[k]=%d\n",k,contig[m].trimR[0]-k,pmap->site[0][contig[m].trimR[0]-k],sitemap[0][m][k]);
      } else {
	printf("  map=%d,flip=%d,numsite=%d:\n",m,flip[m],contig[m].numsite[0]);
	if(!flip[m])
	  for(register int k = 0; k <= contig[m].numsite[0]+1; k++)
	    printf("    k=%d:site[k]=%0.3f,sitemap[k]=%d\n",k,contig[m].site[0][k],sitemap[0][m][k]);
	else
	  for(register int k = 0; k <= contig[m].numsite[0]+1; k++)
	    printf("    k=%d:site[%d]=%0.3f,sitemap[k]=%d\n",k,contig[m].numsite[0]+1-k,contig[m].site[0][contig[m].numsite[0]+1-k],sitemap[0][m][k]);
      }
    }
    delete [] maporder;
  }
  fflush(stdout);
}

void Ccontig::flatten()
{
  if(mapid >= 0)
    return;

  /* first flatten all sub-contigs and count total number of maps */
  int newnummaps = 0;
  int numcontigs = 0;
  for(register int i = 0; i < nummaps;i++){
    contig[i].flatten();
    if(contig[i].mapid >= 0)
      newnummaps++;
    else {
      numcontigs++;
      newnummaps += contig[i].nummaps;
    }
  }
  if(!numcontigs)/* all sub-contigs are simple maps : nothing left to flatten */
    return;

  if(VERB>=2){
    printf("Flattening contig:nummaps=%d -> %d, numcontigs=%d\n", nummaps,newnummaps,numcontigs);
    fflush(stdout);
  }

  /* allocate new arrays for sub-contigs */
  register Ccontig *newcontig = new Ccontig[newnummaps];  
  register int *newflip = new int[newnummaps];
  register int **newsitemap = new int*[newnummaps];

  /* move each map from each sub-contig into newcontig[] and translate its sitemap into newsitemap */
  register int mapcnt = 0;
  for(register int i = 0; i < nummaps; i++){
    register Ccontig *pi = &contig[i];
    if(pi->mapid >= 0){/* just move the single map contig */
      newcontig[mapcnt] = *pi;
      newflip[mapcnt] = flip[i];
      newsitemap[mapcnt] =  sitemap[0][i];
      sitemap[0][i] = 0;
      mapcnt++;
      continue;
    }

    register int N = pi->numsite[0];
    register int inummaps = pi->nummaps;
    for(register int j = 0; j < inummaps; j++){
      if(DEBUG) assert(pi->contig[j].mapid >= 0);
      register int M = pi->contig[j].numsite[0];
      newcontig[mapcnt] = pi->contig[j];
      register int flipped = flip[i];
      newflip[mapcnt] = pi->flip[j] ^ flipped;
      newsitemap[mapcnt] = new int[M+3];
      *newsitemap[mapcnt]++ = -1;
      /* now translate pi->sitemap[j][0..M+1] using sitemap[i][0..N+1] and flipped 
	 and generate newsitemap[mapcnt][0..M+1] */
      register int *np = newsitemap[mapcnt];
      register int *ip = sitemap[0][i];
      register int *jp = pi->sitemap[0][j];
      if(flipped)
	for(register int k=0; k <= M+1; k++)
	  np[k] = ip[N + 1 - jp[k]];
      else
	for(register int k=0; k <= M+1; k++)
	  np[k] = ip[jp[k]];
      mapcnt++;
    }
  }
  
  /* delete previous contig[],flip[] and sitemap[] */
  for(register int i=0; i < nummaps; i++)
    if(sitemap[0][i])
      delete [] &sitemap[0][i][-1];
  delete [] sitemap[0];
  delete [] flip;
  delete [] contig;
  contig = newcontig;
  flip = newflip;
  sitemap[0] = newsitemap;
  nummaps = newnummaps;
  if(DEBUG)
    for(register int i = 0; i < nummaps; i++)
      assert(contig[i].mapid >= 0);
}

/** flip the entire contig along with all sub-contigs */
void Ccontig::contigflip()
{
  if(VERB>=2){
    printf("flipping contig:mapid=%d,numsite=%d,nummaps=%d\n",mapid,numsite[0],nummaps);
    fflush(stdout);
  }

  if(DEBUG) assert(mapid < 0);// cannot flip trivial contig

  Ccontig *p = new Ccontig();
  p->mapid = mapid;

  register int N = numsite[0];

  p->numsite[0] = N;
  p->site[0] = new double[N+2];
  //  p->color = new char[N+1];
  p->sitecnt[0] = new float[N+1];
  p->sitecntFN[0] = new float[N+2];
  p->fragcnt[0] = new float[N+2];
  p->fragcntT[0] = new float[N+2];
  if(CovNorm)
    p->fragcntTnorm[0] = new float[N+2];
  if(TrimNorm >= 0){
    p->sitecntFNnorm[0] = new float[N+2];
    p->sitecntN2[0] = new float[N+2];
    p->sitecntN3[0] = new float[N+2];
    p->sitecntN4[0] = new float[N+2];
    p->sitecntN5[0] = new float[N+2];
    p->sitecntN6[0] = new float[N+2];

    p->fragSd[0] = new float[N+2];
    p->expSd[0] = new float[N+2];
    p->fragBias[0] = new float[N+2];
    p->fragCov[0] = new float[N+2];
    p->fragChiSq[0] = new float[N+2];
  }
  p->site[0][0] = 0.0;
  p->site[0][N+1] = site[0][N+1];
  for(register int i = 1; i <= N; i++){
    p->site[0][i] = site[0][N+1] - site[0][N+1-i];
    //    p->color[i] = color[N+1-i];
    p->sitecnt[0][i] = sitecnt[0][N+1-i];
    p->sitecntFN[0][i] = sitecntFN[0][N+1-i];
  }
  for(register int i = 0; i <= N; i++)
    p->fragcnt[0][i] = fragcnt[0][N-i];

  p->nummaps = nummaps;
  p->contig = new Ccontig[nummaps];
  p->flip = new int[nummaps];
  p->sitemap[0] = new int*[nummaps];
  for(register int k=0; k < nummaps; k++){
    register int M = contig[k].numsite[0];
    p->flip[k] = flip[k] ^ 1;
    p->sitemap[0][k] = new int[M+3];
    *(p->sitemap[0][k])++ = -1;
    p->contig[k] = contig[k];
    contig[k].nummaps = contig[k].numsite[0] = 0;

    for(register int i = 0; i <= M+1; i++)
      p->sitemap[0][k][M+1-i] = N+1 - sitemap[0][k][i];
  }

  allfree();
  *this = *p;
  p->numsite[0] = p->nummaps = 0;
  delete p;
}

/** Merge two contigs into a single contig based on the alignment between one map in each contig */
Ccontig::Ccontig(register Ccontig *contig1, 
		 register Ccontig *contig2,
#if CALIGN_SMALL
		 register CalignA *align,
#else
		 register Calign *align,
#endif
		 int rverb)
{
  if(DEBUG) assert(align!=0 && contig1 != 0 && contig2 != 0);
  if(VERB && (VERB>=2 || rverb)){
    printf("Merging contig1=x%p(nummaps=%d,numsite=%d) and contig2=x%p(nummaps=%d,numsite=%d) using align(mapid1=N%d,id1=%lld,mapid2=N%d,id2=%lld)\n",
	   contig1,contig1->nummaps,contig1->numsite[0],contig2,contig2->nummaps,contig2->numsite[0],nodeorder[align->mapid1],Gmap[align->mapid1]->id,nodeorder[align->mapid2],Gmap[align->mapid2]->id);
    fflush(stdout);
  }
  if(DEBUG) assert(contig1->nummaps > 0 && contig1->numsite[0] > 0);
  if(DEBUG) assert(contig2->nummaps > 0 && contig2->numsite[0] > 0);

  /* set default initialization */
  init();

  contig1->flatten();
  contig2->flatten();

  if(contig1->mapid>=0 || contig2->mapid >= 0){
    printf("Contig Merge not implemented for trivial contigs with mapid>=0\n");
    fflush(stdout);
    assert(0);
  }

  register int mapid1 = align->mapid1;
  register int mapid2 = align->mapid2;
  register int swapped = 0;

  register int index1 = contig1->findmap(mapid1);
  register int index2 = contig2->findmap(mapid2);
  if(index1 < 0){/* swap contig1 and contig2 */
    register Ccontig *contig = contig1;
    contig1 = contig2;
    contig2 = contig;
    swapped = 1;
    index1 = contig1->findmap(mapid1);
    index2 = contig2->findmap(mapid2);
  }
  if(index1 < 0){
    printf("Contig Merge: Cannot find index for mapid1=%d(N%d) in contig1(nummaps=%d) or contig2(nummaps=%d):swapped=%d\n",
	   mapid1,nodeorder[mapid1],contig1->nummaps,contig2->nummaps,swapped);
    printf("  align:(mapid1=%d(N%d),id1=%lld,mapid2=%d(N%d),id2=%lld,or=%d,numpairs=%d,logPV=%0.2f,score=%0.3f)\n",
	   mapid1,nodeorder[mapid1],Gmap[mapid1]->id,mapid2,nodeorder[mapid2],Gmap[mapid2]->id,align->orientation,align->numpairs,align->logPV,align->score);
    for(register int i = 0; i < contig1->nummaps;i++)
      printf("contig1->contig[%d]->mapid=%d(N%d)\n",i,contig1->contig[i].mapid,nodeorder[contig1->contig[i].mapid]);
    for(register int i = 0; i < contig2->nummaps;i++)
      printf("contig2->contig[%d]->mapid=%d(N%d)\n",i,contig2->contig[i].mapid,nodeorder[contig2->contig[i].mapid]);
    fflush(stdout);
    exit(1);
  }
  if(index2 < 0){
    printf("Contig Merge: Cannot find index for mapid2=N%d in contig1(nummaps=%d) or contig2(nummaps=%d):swapped=%d\n",
	   nodeorder[mapid2],contig1->nummaps,contig2->nummaps,swapped);
    for(register int i = 0; i < contig1->nummaps;i++)
      printf("contig1->contig[%d]->mapid=N%d\n",i,nodeorder[contig1->contig[i].mapid]);
    for(register int i = 0; i < contig2->nummaps;i++)
      printf("contig2->contig[%d]->mapid=N%d\n",i,nodeorder[contig2->contig[i].mapid]);
    fflush(stdout);
    exit(1);
  }

  if(VERB && (VERB>=2 || rverb)){
    printf("After adjusting order : Merging contig1(nummaps=%d,numsite=%d) and contig2(nummaps=%d,numsite=%d) using align(mapid1=N%d,mapid2=N%d),index1=%d,index2=%d,flip1=%d,flip2=%d\n",
	   contig1->nummaps,contig1->numsite[0],contig2->nummaps,contig2->numsite[0],nodeorder[align->mapid1],nodeorder[align->mapid2],index1,index2,contig1->flip[index1],contig2->flip[index2]);
    fflush(stdout);
  }

  mapid = -1;
  nummaps = 2;
  contig = new Ccontig[2];
  flip = new int[2];
  sitemap[0] = new int*[2];
  contig[0] = *contig1;
  contig[1] = *contig2;
  contig1->nummaps = contig1->numsite[0] = 0;/* mark as free */
  contig2->nummaps = contig2->numsite[0] = 0;/* mark as free */
  for(register int i = 0; i <= 1; i++){
    sitemap[0][i] = new int[contig[i].numsite[0]+3];
    *sitemap[0][i]++ = -1;
    if(DEBUG)
      for(register int k = 0; k <= contig[i].numsite[0]+1; k++)
	sitemap[0][i][k] = -2;
    flip[i] = 0;
  }

  if(DEBUG) assert(align->numpairs > 0);

  /* force mapid1 into orientation = 0 */
  if(contig[0].flip[index1]){
    if(VERB && (VERB>=2 || rverb)){
      printf("Flipping contig1 (contig[0].flip[index1]=%d)\n",contig[0].flip[index1]);
      fflush(stdout);
    }
    contig[0].contigflip();
    if(DEBUG) assert(!contig[0].flip[index1]);
    if(VERB && (VERB>=2 || rverb)){
      printf("contig1 = contig[0]:numsite=%d:\n",contig[0].numsite[0]);
      for(register int t = 0; t <= contig[0].numsite[0]+1; t++){
	if(t<=0)
	  printf(" t=%d:site[t]=%0.3f\n",t,contig[0].site[0][0]);
	else if(t > contig[0].numsite[0])
	  printf(" t=%d:site[t]=%0.3f,fragcnt[t-1]=%0.2f\n",t,contig[0].site[0][t],contig[0].fragcnt[0][t-1]);
	else
	  printf(" t=%d:site[t]=%0.3f,fragcnt[t-1]=%0.2f,sitecnt[t]=%0.2f,sitecntFN[t]=%0.2f\n",t,contig[0].site[0][t],contig[0].fragcnt[0][t-1],contig[0].sitecnt[0][t],contig[0].sitecntFN[0][t]);
      }
      fflush(stdout);
    }
  }

  /* forces mapid2 into orientation == align->orientation */
  if((contig[1].flip[index2] ^ align->orientation)){
    if(VERB && (VERB>=2 || rverb)){
      printf("Flipping contig2 (contig[1].flip[index2]=%d,align->orientation=%d)\n",contig[1].flip[index2],align->orientation);
      fflush(stdout);
    }
    contig[1].contigflip();
    if(DEBUG) assert(!(contig[1].flip[index2] ^ align->orientation));
    if(VERB && (VERB>=2 || rverb)){
      printf("contig2 = contig[1]:numsite=%d:\n",contig[1].numsite[0]);
      for(register int t = 0; t <= contig[1].numsite[0]+1; t++){
	if(t<=0)
	  printf(" t=%d:site[t]=%0.3f\n",t,contig[1].site[0][0]);
	else if(t > contig[1].numsite[0])
	  printf(" t=%d:site[t]=%0.3f,fragcnt[t-1]=%0.2f\n",t,contig[1].site[0][t],contig[1].fragcnt[0][t-1]);
	else
	  printf(" t=%d:site[t]=%0.3f,fragcnt[t-1]=%0.2f,sitecnt[t]=%0.2f,sitecntFN[t]=%0.2f\n",t,contig[1].site[0][t],contig[1].fragcnt[0][t-1],contig[1].sitecnt[0][t],contig[1].sitecntFN[0][t]);
      }
      fflush(stdout);
    }
  }

  /* first compute number of sites in new contig */
  int N[2];/* number of sites in each contig */
  N[0] = contig[0].numsite[0];
  N[1] = contig[1].numsite[0];
  int M[2];/* number of sites in the aligned maps of each contig */
  M[0] = contig[0].contig[index1].numsite[0];
  M[1] = contig[1].contig[index2].numsite[0];
  register int K = align->numpairs;

  extern int *map2nodeid;

  Cnode *pnode1, *pnode2;
  if(map2nodeid){
    pnode1 = &node[map2nodeid[mapid1]];
    pnode2 = &node[map2nodeid[mapid2]];
  } else {
    pnode1 = &node[mapid1];
    pnode2 = &node[mapid2];
  }
  if(DEBUG) assert(pnode1->mapid == mapid1);
  if(DEBUG) assert(pnode2->mapid == mapid2);

  int trimL1 = pnode1->trimL;
  int trimR1 = pnode1->trimR;
  int trimL2 = align->orientation ? Gmap[mapid2]->numsite[0]+1 - pnode2->trimR : pnode2->trimL;
  int trimR2 = align->orientation ? Gmap[mapid2]->numsite[0]+1 - pnode2->trimL : pnode2->trimR;

  int JL[2];/* left most aligned site of each contig */
  JL[0] = JL[1] = -1;
  int Gmap = -1, Hmap = -1, Imap = -1,Jmap = -1;// Save index into maps for debugging
  for(register int k=0; k < K; k++){
    register int i = align->sites1[k];
    if(i <= trimL1 || i >= trimR1)
      continue;
    register int j = align->sites2[k];
    if(j <= trimL2 || j >= trimR2)
      continue;
    if(DEBUG) assert(0 <= i-trimL1 && i-trimL1 <= contig[0].contig[index1].numsite[0]+1);
    if(DEBUG) assert(0 <= j-trimL2 && j-trimL2 <= contig[1].contig[index2].numsite[0]+1);
    JL[0] = contig[0].sitemap[0][index1][i-trimL1];
    JL[1] = contig[1].sitemap[0][index2][j-trimL2];
    if(DEBUG) assert(JL[0] >= -1 && JL[1] >= -1);
    if(JL[0] >= 0 && JL[1] >= 0){
      Imap = i - trimL1;
      Jmap = j - trimL2;
      break;
    }
  }
  if(DEBUG && !(JL[0] >= 0 && JL[1] >= 0)){
    printf("mapid1=%d,mapid2=%d:trimL1=%d,trimR1=%d,numsite1=%d,trimL2=%d,trimR2=%d,numsite2=%d:numpairs=%d,JL[0]=%d,JL[1]=%d\n",
	   mapid1,mapid2,trimL1,trimR1,contig[0].numsite[0],trimL2,trimR2,contig[1].numsite[0],align->numpairs,JL[0],JL[1]);
    for(register int k=0;k < K; k++)
      printf("k=%d:align->sites1[k]=%d,align->sites2[k]=%d\n",k,align->sites1[k],align->sites2[k]);
    printf("Unable to find leftmost aligned site within range of trimL..trimR !\n");
    fflush(stdout);
    assert(JL[0] >= 0 && JL[1] >= 0);
  }
  int L[2];/* left end of aligned maps on each contig (can be 0) */
  if(align->Lend <= -1){
    L[0] = contig[0].sitemap[0][index1][0];
    L[1] = contig[1].sitemap[0][index2][0];
    if(DEBUG) assert(L[0] >= 0 && L[1] >= 0);
  } else {
    L[0] = JL[0];
    L[1] = JL[1];
  }

#if INTERLEAVE==0
  /* decide which left end to copy */
  register int left = (contig[0].site[0][JL[0]] >= contig[1].site[0][JL[1]]) ? 0 : 1;
  numsite[0] = JL[left]; /* sites of the left end that will be copied up to (including) the first aligned site */
  numsite[0] += JL[1-left] - L[1-left];/* add the sites of the other contig to the left of the first aligned site up to the left end of the aligned map */
#else
  /* interleave both left ends */
  numsite[0] = JL[0] + JL[1];/* includes left end of shorter end and merged first aligned site */
#endif

  register int G = JL[0];
  register int H = JL[1];
  register int I= -1,J = -1;
  for(register int k = 1; k < K; k++, G = I, H = J){
    for(;k < K; k++){
      register int i = align->sites1[k];
      if(i <= trimL1 || i >= trimR1)
	continue;
      register int j = align->sites2[k];
      if(j <= trimL2 || j >= trimR2)
	continue;
      if(DEBUG) assert(0 <= i-trimL1 && i-trimL1 <= contig[0].contig[index1].numsite[0]+1);
      if(DEBUG) assert(0 <= j-trimL2 && j-trimL2 <= contig[1].contig[index2].numsite[0]+1);
      I = contig[0].sitemap[0][index1][i-trimL1];
      J = contig[1].sitemap[0][index2][j-trimL2];
      if(DEBUG) assert(I >= -1 && J >= -1);
      if(I > G && J > H)
	break;
    }
    if(k >= K)
      break;
    if(DEBUG) assert(I > G && J > H);
    numsite[0] += (I-G)+(J-H)-1;
  }
  int JR[2];/* right most aligned site of each contig */
  JR[0] = G;
  JR[1] = H;
  if(DEBUG) assert(JR[0] <= N[0] && JR[1] <= N[1]);
  int R[2];/* right end of aligned maps on each contig (can be the right end of the contig at N[]+1) */
  if(align->Rend <= -1){
    R[0] = contig[0].sitemap[0][index1][M[0]+1];
    R[1] = contig[1].sitemap[0][index2][M[1]+1];
    if(DEBUG) assert(R[0] >= 0 && R[1] >= 0 && R[0] <= N[0]+1 && R[1] <= N[1]+1);
  } else {
    R[0] = JR[0];
    R[1] = JR[1];
  }

#if INTERLEAVE==0
  /* decide which contig right end to copy */
  register int right = (contig[0].site[0][N[0]+1] - contig[0].site[0][JR[0]] >= 
			contig[1].site[0][N[1]+1] - contig[1].site[0][JR[1]]) ? 0 : 1;
  numsite[0] += N[right] - JR[right];
  numsite[0] += R[1-right] - JR[1-right];/* add sites of other contig to the right of the last aligned site up to the right end of the aligned map */
#else
  /* interleave both right ends */
  numsite[0] += N[0]-JR[0] + N[1]-JR[1] + 1;/* includes right end of shorter end */
#endif

  /* now allocate arrays */
  site[0] = new double[numsite[0]+2];
  //  color = new char[numsite+1];
  sitecnt[0] = new float[numsite[0]+1];
  sitecntFN[0] = new float[numsite[0]+2];
  fragcnt[0] = new float[numsite[0]+2];
  fragcntT[0] = new float[numsite[0]+2];
  if(CovNorm)
    fragcntTnorm[0] = new float[numsite[0]+2];
  if(TrimNorm >= 0){
    sitecntFNnorm[0] = new float[numsite[0]+2];
    sitecntN2[0] = new float[numsite[0]+2];
    sitecntN3[0] = new float[numsite[0]+2];
    sitecntN4[0] = new float[numsite[0]+2];
    sitecntN5[0] = new float[numsite[0]+2];
    sitecntN6[0] = new float[numsite[0]+2];

    fragSd[0] = new float[numsite[0]+2];
    expSd[0] = new float[numsite[0]+2];
    fragBias[0] = new float[numsite[0]+2];
    fragCov[0] = new float[numsite[0]+2];
    fragChiSq[0] = new float[numsite[0]+2];
  }

#if INTERLEAVE==0
  J = L[1-left];/* next site on contig[1-left] to copy */
  /* copy contig[left] up to site J on contig[1-left] mapped to contig[left] */
  /* NOTE : unaligned parts of contig[1-left] have no effect on sitecntFN[] or fragcnt[] and sitecnt[1-left] are set to -1 */
  site[0][0] = 0.0;
  sitemap[0][left][0] = 0;
  for(I = 1; contig[left].site[0][JL[left]] - contig[left].site[0][I] > contig[1-left].site[0][JL[1-left]] - contig[1-left].site[0][J]; I++){
    if(DEBUG) assert(I <= JL[left]);
    site[0][I] = contig[left].site[0][I];
    if(DEBUG)assert(site[0][I] >= site[0][I-1]);
    //    color[I] = contig[left].color[I];
    sitecnt[0][I] = contig[left].sitecnt[0][I];
    sitecntFN[0][I] = contig[left].sitecntFN[0][I];
    fragcnt[0][I-1] = contig[left].fragcnt[0][I-1];
    sitemap[0][left][I] = I;
  }
  if(DEBUG) assert(J <= JL[1-left]);
  if(DEBUG) assert(I > 0);
  for(register int i = 0; i < J;i++)
    sitemap[0][1-left][i] = -1;
  register int next = I;/* next site on merged contig */
  /* copy contig[left] from site I to JL[left] and contig[1-left] from site J to JL[1-left] (the last sites are merged, all others interleaved) */
  while(I < JL[left] || J < JL[1-left]){
    if(J >= JL[1-left] || (I < JL[left] && contig[left].site[0][JL[left]] - contig[left].site[0][I] >
			   contig[1-left].site[0][JL[1-left]] - contig[1-left].site[0][J])){/* copy site I on contig[left] */
      if(DEBUG) assert(I < JL[left]);
      site[0][next] = contig[left].site[0][I];
      if(DEBUG) assert(site[0][next] >= site[0][next-1]);
      //      color[next] = contig[left].color[I];
      sitecnt[0][next] = contig[left].sitecnt[0][I];
      sitecntFN[0][next] = contig[left].sitecntFN[0][I] + (J > 0 ? contig[1-left].fragcnt[0][J-1] : 0.0);
      fragcnt[0][next-1] = contig[left].fragcnt[0][I-1] + (J > 0 ? contig[1-left].fragcnt[0][J-1] : 0.0);
      sitemap[0][left][I++] = next++ ;
    } else { /* copy site J on contig[1-left] */
      if(DEBUG) assert(J < JL[1-left]);
      if(DEBUG) assert(I>0 && next > 0);
      site[next] = contig[left].site[JL[left]] - (contig[1-left].site[JL[1-left]] - contig[1-left].site[J]);
      if(DEBUG) assert(site[next] >= site[next-1]);
      if(J == 0){
	//	color[next] = contig[1-left].color[J+1];
	sitecnt[0][next] = 0;
	sitecntFN[0][next] = contig[left].fragcnt[0][I-1];
	fragcnt[0][next-1] = contig[left].fragcnt[0][I-1];
      } else {
	//	color[next] = contig[1-left].color[J];
	sitecnt[0][next] = contig[1-left].sitecnt[0][J];
	sitecntFN[0][next] = contig[1-left].sitecntFN[0][J] + contig[left].fragcnt[0][I-1];
	fragcnt[0][next-1] = contig[1-left].fragcnt[0][J-1] + contig[left].fragcnt[0][I-1];
      }
      sitemap[0][1-left][J++] = next++;
    }
  }
  if(DEBUG) assert(I==JL[left] && J == JL[1-left]);
  /* next is first merged site */
  site[next] = contig[left].site[I];
  if(DEBUG) assert(site[next] >= site[next-1]);
  //  color[next] = contig[left].color[I];
  //  if(DEBUG) assert(contig[left].color[I] == contig[1-left].color[J]);
  sitecnt[0][next] = contig[left].sitecnt[0][I] + contig[1-left].sitecnt[0][J];
  sitecntFN[0][next] = contig[left].sitecntFN[0][I] + contig[1-left].sitecntFN[0][J];
  fragcnt[0][next-1] = contig[left].fragcnt[0][I-1] + contig[1-left].fragcnt[0][J-1];
  sitemap[0][left][I] = next;
  sitemap[0][1-left][J] = next++;
#else
  /* interleave left ends up to JL[0] and JL[1] */
  register int next = 0;
  register double y = contig[0].site[0][JL[0]];
  register double x = contig[1].site[0][JL[1]];
  register double z = max(y,x);
  for(I=J=0; I < JL[0] || J < JL[1];){
    if(J >= JL[1] || (I < JL[0] && y - contig[0].site[0][I] > x - contig[1].site[0][J])){/* copy site I on contig[left] */
      if(DEBUG) assert(I < JL[0]);
      site[0][next] = z - (y - contig[0].site[0][I]);
      if(next > 0){
	if(DEBUG) assert(site[0][next] >= site[0][next-1] - 1.0e-8);
	if(I == 0){
	  //	  color[next] = contig[0].color[I+1];
	  sitecnt[0][next] = 0;
	  sitecntFN[0][next] = (J > 0 ? contig[1].fragcnt[0][J-1] : 0.0);
	  fragcnt[0][next-1] = (J > 0 ? contig[1].fragcnt[0][J-1] : 0.0);
	} else {
	  //	  color[next] = contig[0].color[I];
	  sitecnt[0][next] = contig[0].sitecnt[0][I];
	  sitecntFN[0][next] = contig[0].sitecntFN[0][I] + (J > 0 ? contig[1].fragcnt[0][J-1] : 0.0);
	  fragcnt[0][next-1] = contig[0].fragcnt[0][I-1] + (J > 0 ? contig[1].fragcnt[0][J-1] : 0.0);
	}
      }
      sitemap[0][0][I++] = next++ ;
    } else { /* copy site J on contig[1-left] */
      if(DEBUG) assert(J < JL[1]);
      site[0][next] = z - (x - contig[1].site[0][J]);
      if(next > 0){
	if(DEBUG) assert(site[0][next] >= site[0][next-1] - 1.0e-8);
	if(J == 0){
	  //	  color[next] = contig[1].color[J+1];
	  sitecnt[0][next] = 0;
	  sitecntFN[0][next] = (I > 0 ? contig[0].fragcnt[0][I-1] : 0.0);
	  fragcnt[0][next-1] = (I > 0 ? contig[0].fragcnt[0][I-1] : 0.0);
	} else {
	  //	  color[next] = contig[1].color[J];
	  sitecnt[0][next] = contig[1].sitecnt[0][J];
	  sitecntFN[0][next] = contig[1].sitecntFN[0][J] + (I > 0 ? contig[0].fragcnt[0][I-1] : 0.0);
	  fragcnt[0][next-1] = contig[1].fragcnt[0][J-1] + (I > 0 ? contig[0].fragcnt[0][I-1] : 0.0);
	}
      }
      sitemap[0][1][J++] = next++;
    }
  }
  if(DEBUG) assert(I==JL[0] && J == JL[1]);
  if(DEBUG) assert(I>0 && J > 0);
  /* next is first merged site */
  site[0][next] = z;
  if(DEBUG) assert(site[0][next] >= site[0][next-1] - 1.0e-8);
  //  color[next] = contig[0].color[I];
  //  if(DEBUG) assert(contig[0].color[I] == contig[1].color[J]);
  sitecnt[0][next] = contig[0].sitecnt[0][I] + contig[1].sitecnt[0][J];
  sitecntFN[0][next] = contig[0].sitecntFN[0][I] + contig[1].sitecntFN[0][J];
  fragcnt[0][next-1] = contig[0].fragcnt[0][I-1] + contig[1].fragcnt[0][J-1];
  sitemap[0][0][I] = next;
  sitemap[0][1][J] = next++;
#endif

  /* Now copy the overlapped regions of contig[0] and contig[1] */
  G = JL[0];
  H = JL[1];
  Gmap = Imap;
  Hmap = Jmap;
  for(register int k = 1; k < K; k++, G = I, H = J, Gmap = Imap, Hmap = Jmap){
    for(;k < K; k++){
      register int i = align->sites1[k];
      if(i <= trimL1 || i >= trimR1)
	continue;
      register int j = align->sites2[k];
      if(j <= trimL2 || j >= trimR2)
	continue;
      if(DEBUG) assert(0 <= i-trimL1 && i-trimL1 <= contig[0].contig[index1].numsite[0]+1);
      if(DEBUG) assert(0 <= j-trimL2 && j-trimL2 <= contig[1].contig[index2].numsite[0]+1);
      I = contig[0].sitemap[0][index1][i-trimL1];
      J = contig[1].sitemap[0][index2][j-trimL2];
      if(DEBUG) assert(I >= -1 && J >= -1);
      if(I > G && J > H){
	Imap = i - trimL1;
	Jmap = j - trimL2;
	break;
      }
    }
    if(k >= K)
      break;
    if(DEBUG) assert(I > G && J > H);

    /* first compute averaged interval size z in new contig : assumes var(y) = SD[0]*SD[0]*y
       SizeType == 0 : Old code (not a good approximation & fails to detect outliers)
       SizeType >= 1 : Better approximation
       SizeType >= 2 : Also tries to detect outliers (not reliable and order dependent)
    */
    double varI = 0.0;
    for(register int i = G; i < I; i++){
      double mincnt = contig[0].fragcnt[0][i];
      if(DEBUG) assert(mincnt > 0.0);
      varI += (contig[0].site[0][i+1] - contig[0].site[0][i]) / mincnt; // *(SizeType ? mincnt : 1.0/mincnt);
    }
    if(DEBUG) assert(varI > 0.0);

    double varJ = 0.0;
    for(register int j = H; j < J; j++){
      double mincnt = contig[1].fragcnt[0][j];
      if(DEBUG) assert(mincnt > 0.0);
      varJ += (contig[1].site[0][j+1] - contig[1].site[0][j]) / mincnt; // *(SizeType ? mincnt : 1.0/mincnt);
    }
    if(DEBUG) assert(varJ > 0.0);

    register double y = contig[0].site[0][I] - contig[0].site[0][G];/* contig[0] interval size */
    register double x = contig[1].site[0][J] - contig[1].site[0][H];/* contig[1] interval size */
    double wtI = SizeType ? y / varI : 1.0/varI;
    double wtJ = SizeType ? x / varJ : 1.0/varJ;
    if(DEBUG) assert(wtI > 0.0);
    if(DEBUG) assert(wtJ > 0.0);
    register double z = (y*wtI + x *wtJ)/(wtI+wtJ);/* new contig interval size */
    if(DEBUG) assert(isfinite(z));

    /* compute weighting of new size estimate based on Poutlier */
    double zwt = 1.0;
    if(SizeType>=2  && Poutlier > 0.0){
      double var = 2.0 * SF[0]*SF[0] + SD[0]*SD[0]*(x+y);
      double Ksd = fabs(x-y)/sqrt(2.0*var);
      double pv = 0.5 * erfc(Ksd);
      zwt = pv/(pv + Poutlier);
      if(DEBUG) assert(isfinite(zwt));
      if(SizeType>=3 && zwt * (wtI + wtJ) < max(wtI,wtJ) && min(wtI,wtJ) * 4.0 <= max(wtI,wtJ)){/* assume higher weight interval is correct */
	zwt = max(wtI,wtJ)/(wtI+wtJ);
	z = (wtI > wtJ) ? y : x;
      }
      if(DEBUG) assert(zwt >= 0.0);
      zwt = max(0.01,zwt);/* avoid problems caused by zero fragcnt[] values */
    }

    register double start = site[0][next-1];
    if(VERB && (VERB>=2 || rverb)){
      int mapid1 = contig[0].contig[index1].mapid;
      if(DEBUG) assert(mapid1 >= 0 && mapid1 < gnummaps);
      int numsite1 = contig[0].contig[index1].numsite[0];
      if(DEBUG) assert(Imap >= 1 && Imap <= numsite1);
      if(DEBUG && !(1 <= Gmap && Gmap <= numsite1)){
	printf("  k=%d:G=%d,H=%d,I=%d,J=%d:y=%0.3f,x=%0.3f,wtI=%0.3f,wtJ=%0.3f,z=%0.3f,consensus left site=%d,loc=%0.3f\n",k,G,H,I,J,y,x,wtI,wtJ,z,next-1,start);
	printf("  Imap=%d,Gmap=%d,contig[0].contig[index1].numsite=%d,contig[0].flip[index1]=%d\n",
	       Imap,Gmap,contig[0].contig[index1].numsite[0], contig[0].flip[index1]);
	fflush(stdout);
	assert(1 <= Gmap && Gmap <= numsite1);
      }
      double ymap = contig[0].flip[index1] ? 
	gmap[mapid1]->site[0][contig[0].contig[index1].trimR[0] - Gmap] - gmap[mapid1]->site[0][contig[0].contig[index1].trimR[0] - Imap] :
	gmap[mapid1]->site[0][contig[0].contig[index1].trimL[0] + Imap] - gmap[mapid1]->site[0][contig[0].contig[index1].trimL[0] + Gmap];

      int mapid2 = contig[1].contig[index2].mapid;
      if(DEBUG) assert(mapid2 >= 0 && mapid2 < gnummaps);
      int numsite2 = contig[1].contig[index2].numsite[0];
      if(DEBUG) assert(Jmap >= 1 && Jmap <= numsite2);
      if(DEBUG) assert(1 <= Hmap && Hmap <= numsite2);
      double xmap = contig[1].flip[index2] ?
	gmap[mapid2]->site[0][contig[1].contig[index2].trimR[0] - Hmap] - gmap[mapid2]->site[0][contig[1].contig[index2].trimR[0] - Jmap] :
	gmap[mapid2]->site[0][contig[1].contig[index2].trimL[0] + Jmap] - gmap[mapid2]->site[0][contig[1].contig[index2].trimL[0] + Hmap];
      if(fabs(x-y) > 4.0 * sqrt(2.0*SF[0]*SF[0]+SD[0]*SD[0]*(x+y)))
	printf("  k=%d:G=%d,H=%d,I=%d,J=%d:y=%0.3f(ymap=%0.3f),x=%0.3f(xmap=%0.3f),wtI=%0.3f,wtJ=%0.3f,z=%0.3f(wt=%0.3f),consensus left site=%d,loc=%0.3f ??\n",k,G,H,I,J,y,ymap,x,xmap,wtI,wtJ,z,zwt,next-1,start);
      else
	printf("  k=%d:G=%d,H=%d,I=%d,J=%d:y=%0.3f(ymap=%0.3f),x=%0.3f(xmap=%0.3f),wtI=%0.3f,wtJ=%0.3f,z=%0.3f(wt=%0.3f),consensus left site=%d,loc=%0.3f\n",k,G,H,I,J,y,ymap,x,xmap,wtI,wtJ,z,zwt,next-1,start);
      fflush(stdout);
    }
    register int i = G+1;
    register int j = H+1;
    while(i < I || j < J){/* copy either contig[0].site[i] or contig[1].site[j], whichever is leftmost */
      if(DEBUG) assert(next <= numsite[0]);
      if(DEBUG) assert(i <= I && j <= J);
      register double zi = (contig[0].site[0][i] - contig[0].site[0][G])*z/y;
      register double zj = (contig[1].site[0][j] - contig[1].site[0][H])*z/x;
      if(j >= J || (i < I && zi < zj)){/* copy contig[0].site[i] */
	if(DEBUG) assert(i < I);
	site[0][next] = start + zi;
	if(DEBUG) assert(site[0][next] >= site[0][next-1] - 1.0e-8);
	//	color[next] = contig[0].color[i];
	sitecnt[0][next] = contig[0].sitecnt[0][i];
	sitecntFN[0][next] = contig[0].sitecntFN[0][i] + contig[1].fragcnt[0][j-1];
	fragcnt[0][next-1] = zwt*(contig[0].fragcnt[0][i-1] + contig[1].fragcnt[0][j-1]);
	sitemap[0][0][i++] = next++;
      } else {/* copy contig[1].site[j] */
	if(DEBUG) assert(j < J);
	site[0][next] = start + zj;
	if(DEBUG) assert(site[0][next] >= site[0][next-1] - 1.0e-8);
	//	color[next] = contig[1].color[j];
	sitecnt[0][next] = contig[1].sitecnt[0][j];
	sitecntFN[0][next] = contig[1].sitecntFN[0][j] + contig[0].fragcnt[0][i-1];
	fragcnt[0][next-1] = zwt*(contig[1].fragcnt[0][j-1] + contig[0].fragcnt[0][i-1]);
	sitemap[0][1][j++] = next++;
      }
    }
    if(DEBUG) assert(i==I && j==J);
    /* copy aligned sites contig[0].site[0][I] and contig[0].site[0][J] */
    site[0][next] = start + z;
    if(DEBUG) assert(site[0][next] >= site[0][next-1] - 1.0e-8);
    //    color[next] = contig[0].color[I];
    //    if(DEBUG) assert(contig[0].color[I] == contig[1].color[J]);
    sitecnt[0][next] = contig[0].sitecnt[0][I] + contig[1].sitecnt[0][J];
    sitecntFN[0][next] = contig[0].sitecntFN[0][I] + contig[1].sitecntFN[0][J];
    fragcnt[0][next-1] = zwt*(contig[0].fragcnt[0][I-1] + contig[1].fragcnt[0][J-1]);
    sitemap[0][0][I] = next;
    sitemap[0][1][J] = next++;
  }
  if(DEBUG) assert(G == JR[0] && H == JR[1]);
  if(DEBUG) assert(next <= numsite[0]+1);

#if INTERLEAVE==0
  J = R[1-right];/* last site on contig[1-right] to copy (could be N[1-right] + 1) */
  /* copy overlapped right end of aligned maps till contig[1-right].site[J] */
  register int i = JR[right] + 1;
  register int j = JR[1-right] + 1;
  register double y = contig[right].site[0][JR[right]];
  register double x = contig[1-right].site[0][JR[1-right]];
  register double start = site[0][next-1];
  while(j <= J || contig[right].site[0][i] - y < contig[1-right].site[0][J] - x){
    if(DEBUG) assert(j <= J+1);
    if(j > J || contig[right].site[0][i] - y < contig[1-right].site[0][j] - x){/* copy contig[right].site[i] */
      site[0][next] = start + contig[right].site[0][i] - y;
      /*      printf("site[next=%d]=%0.4f,start=%0.4f,contig[right=%d].site[i=%d]=%0.4f,y=%0.4f(JR[0]=%d,JR[1]=%d,J=%d,i=%d,j=%d)\n",
	     next,site[next],start,right,i,contig[right].site[i],y,JR[0],JR[1],J,i,j);
	     fflush(stdout);*/
      if(DEBUG) assert(site[0][next] >= site[0][next-1]);
      //      color[next] = contig[right].color[i];
      sitecnt[0][next] = contig[right].sitecnt[0][i];
      sitecntFN[0][next] = contig[right].sitecntFN[0][i] + contig[1-right].fragcnt[0][j-1];
      fragcnt[0][next-1] = contig[right].fragcnt[0][i-1] + contig[1-right].fragcnt[0][j-1];
      sitemap[0][right][i++] = next++;
    } else {/* copy contig[1-right].site[j] */
      site[0][next] = start + contig[1-right].site[0][j] - x;
      /*      printf("site[next=%d]=%0.4f,start=%0.4f,contig[1-right=%d].site[j=%d]=%0.4f,x=%0.4f(JR[0]=%d,JR[1]=%d,J=%d,i=%d,j=%d)\n",
	     next,site[0][next],start,1-right,j,contig[1-right].site[0][j],x,JR[0],JR[1],J,i,j);
	     fflush(stdout);*/
      if(DEBUG) assert(site[0][next] >= site[0][next-1]); 
   }
  }
  if(DEBUG) assert(next <= numsite[0]+1);
  if(DEBUG) assert(i <= N[right]+1);
  if(DEBUG) assert(j == J+1);
  for(;j <= N[1-right] + 1; j++)
    sitemap[0][1-right][j] = -1;
  /* copy remaining sites on contig[right].site[i..N[right]+1] */
  for(;i <= N[right];){
    site[0][next] = start + contig[right].site[0][i] - y;
    /*    printf("site[next=%d]=%0.4f,start=%0.4f,contig[right=%d].site[i=%d]=%0.4f,y=%0.4f(JR[right]=%d,N[right]=%d)\n",
	   next,site[0][next],start,right,i,contig[right].site[0][i],y,JR[right],N[right]);
	   fflush(stdout);*/
    if(DEBUG) assert(site[0][next] >= site[0][next-1]);
    //    color[next] = contig[right].color[i];
    sitecnt[0][next] = contig[right].sitecnt[0][i];
    sitecntFN[0][next] = contig[right].sitecntFN[0][i];
    fragcnt[0][next-1] = contig[right].fragcnt[0][i-1];
    sitemap[0][right][i++] = next++;
  }
  if(DEBUG) assert(i == N[right]+1);
  if(DEBUG) assert(next == numsite[0]+1);
  site[0][next] = start + contig[right].site[0][i] - y;
  if(DEBUG) assert(site[0][next] >= site[0][next-1]);
  fragcnt[0][next-1] = contig[right].fragcnt[0][i-1];
  sitemap[0][right][i] = next;
#else
  /* interleave all remaining site right of JR[0] and JR[1] */
  register int i = JR[0]+1;
  register int j = JR[1]+1;
  y = contig[0].site[0][JR[0]];
  x = contig[1].site[0][JR[1]];
  register double start = site[0][next-1];
  while(i <= N[0]+1 || j <= N[1]+1){
    if(j > N[1]+1 || (i <= N[0] + 1 && contig[0].site[0][i] - y < contig[1].site[0][j] - x)){/* copy contig[0].site[i] */
      if(DEBUG) assert(i <= N[0]+1);
      site[0][next] = start + contig[0].site[0][i] - y;
      if(DEBUG) assert(site[0][next] >= site[0][next-1] - 1.0e-8);
      if(next <= numsite[0]){
	if(i > N[0]){
	  //	  color[next] = contig[0].color[i-1];
	  sitecnt[0][next] = 0;
	  sitecntFN[0][next] = (j <= N[1]+1) ? contig[1].fragcnt[0][j-1] : 0.0;
	} else {
	  //	  color[next] = contig[0].color[i];
	  sitecnt[0][next] = contig[0].sitecnt[0][i];
	  sitecntFN[0][next] = contig[0].sitecntFN[0][i] + (j <= N[1]+1 ? contig[1].fragcnt[0][j-1] : 0.0);
	}
      }
      fragcnt[0][next-1] = contig[0].fragcnt[0][i-1] + (j <= N[1]+1 ? contig[1].fragcnt[0][j-1] : 0.0);
      sitemap[0][0][i++] = next++;
    } else {/* copy contig[1].site[j] */
      if(DEBUG) assert(j <= N[1]+1);
      site[0][next] = start + contig[1].site[0][j] - x;      
      if(DEBUG) assert(site[0][next] >= site[0][next-1] - 1.0e-8);
      if(next <= numsite[0]){
	if(j > N[1]){
	  //	  color[next] = contig[1].color[j-1];
	  sitecnt[0][next] = 0;
	  sitecntFN[0][next] = (i <= N[0]+1 ? contig[0].fragcnt[0][i-1] : 0.0);
	} else {
	  //	  color[next] = contig[1].color[j];
	  sitecnt[0][next] = contig[1].sitecnt[0][j];
	  sitecntFN[0][next] = contig[1].sitecntFN[0][j] + (i <= N[0]+1 ? contig[0].fragcnt[0][i-1] : 0.0);
	}
      }
      fragcnt[0][next-1] = contig[1].fragcnt[0][j-1] + (i <= N[0]+1 ? contig[0].fragcnt[0][i-1] : 0.0);
      sitemap[0][1][j++] = next++;
    }
  }
  if(DEBUG) assert(next == numsite[0]+2);
#endif

  if(DEBUG){
    /* check that sitemap was fully initialized */
    for(register int i=0; i <= 1; i++)
      for(register int k = -1; k <= N[i]+1;k++)
	assert(sitemap[0][i][k] >= -1);
    /* check that sites are monotonic */
    for(register int i = 1; i <= numsite[0]+1; i++){
      assert(site[0][i] >= site[0][i-1] - 1.0e-8);
      if(site[0][i-1] > site[0][i])
	site[0][i] = site[0][i-1];
    }

    float cnt = contig[0].nummaps + contig[1].nummaps + 0.0001;
    /* check that fragcount is > 0 and does not exceed the number of maps */
    for(register int i = 0; i <= numsite[0]; i++){
      if(DEBUG && !(fragcnt[0][i] > 0.0 && fragcnt[0][i] <= cnt)){
	printf("fragcnt[i=%d]= %0.8e,cnt=%0.3f\n",i,fragcnt[0][i],cnt);
	fflush(stdout);
	assert(fragcnt[0][i] > 0.0 && fragcnt[0][i] <= cnt);
      }
    }
  }
}
