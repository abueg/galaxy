#include <stdlib.h>
#include <stdio.h>
//#include <math.h>
#include <omp.h>

#include <vector>
#include <stack>
#include <queue>

#ifndef __ARM_ARCH
#include <emmintrin.h>

#ifdef __SSE3__
#include <pmmintrin.h>

#ifdef __SSSE3__
#include <tmmintrin.h>

#endif
#endif
#endif 

#include <sys/mman.h>

#include <malloc.h>

//#undef DEBUG
//#define DEBUG 2


#include "Assembler.h"
#include "Assembler_parameters.h"

#include "Cmap.h"

#if CALIGN_SMALL
#include "CalignA.h"
#define Calign CalignA
#define alignment alignmentA
#else
#include "Calign.h"
#endif


#include "Cnode.h"
#include "Ccontig.h"

#define GRAPH_SMALL 1 /* WAS30 0 */ // use CnodeA and CedgeA in graph_build() and graph_nodetrim() to save memory, then convert to Cnode/Cedge after swapping out alignments

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Assembler_graph.cpp 11381 2020-07-28 23:26:39Z tanantharaman $");

#define HUGE_COV 1e+300 /* largest allowed coverage value for edges */

#define FAST 0  /* 0 : Full implementation of redundancy counting in graph_redundancy(): can be slow due to cycles in graph that have to be checked. \
                   2 : fastest, but least accurate redundancy counting (Each deleted edge is added to only the first other pair of edges it aligns with) */
#define CFAST 2 /* 0 : Full implementation \
		   1 : faster, but less accurate redundancy counting (In first iteration, Each deleted edge is only added to pairs of edges it aligns with if they are both shorter edges OR \
		      this is the first pair of edges for this deleted edge) \
		   >= 2 : In first CFAST-1 iterations, Each deleted edge is only added to pairs of edges it aligns with if they are both shorter edges. */

#define AVCOV 1 /* 0 : use the minimum edge coverage for collapsed chains \
		   1 : Use the average edge coverage for collapsed chains (average coverage = sum(distance)/sum(distance/cov)) */

#define FORKBREAK 0 /* TRY 1 */ /* before CROSSBREAK, add cross links to Forking nodes : may cause loss of maps in final consensus and produce too many cross edges */

#define CROSSBREAK 1 /* At end of graph_redudundancy, for each node check if there are more than two dominant edges and check if any pair of them does not have a strong direct edge (with valid < 0) : \
			  If so replace the weaker of the two edges with the direct edge to any of the other edges with a strong direct edge (change valid from -1 to 1) */
#define CROSSBREAK_FIX 1 /* fix updated linked list */

#define CYCLEBREAK  1 /* 0 : break cycles based on lowest coverage < 2nd lowest coverage * CycleCoverageRatio (No longer supported) \
		         1 : break cycles based on lowest coverage < average coverage *  CycleCoverageRatio (provided lowest coverage <= BulgeBreakMaxCov) */
#define MAXCYCLE 0 /* WAS 1 */ /* limit cycle length to MaxBulge (even if cycle is already found) : should not be needed */
#define CYCLESTEPS 1 /* WAS 0 */ /* limit cycle length to BFSdepth*2 steps (in addition to MaxBulge kb) during search */

#define CHIMTRIM 1  /* 0 : For nodes connected to disjoint sets of nodes, break the smaller set of edges  \
 		       1 : Also for almost disjoint sets of nodes (number of cross edges <= "expected number" * ChimTrimRatio) */

#define EDGETRIM 2  /* 0 : trim edge if max(1,Ccnt) is smaller than max(Lcntmax,Rcntmax)*EdgeTrimRatio for both nodes \
				   1 : trim edge if max(1,Ccnt) is smaller than min(Lcntmax,Rcntmax)*EdgeTrimRatio for at least one node  \
				   2 : trim edge if max(1,Ccnt) is smaller than min(Lcntmax,Rcntmax)*EdgeTrimRatio for both nodes \
				*/
#define EDGECHECK 1 /* in graph_edgecnt() don't count cases where the orientation and distances are not consistent (as in graph_redundancy) */
#define EDGEBRIDGE 0  /* >= 1 : Compute Bcnt (the number of briding edges) for each edge \
			 >= 2 : Add Bcnt to Ccnt (causes chimtrim() to delete too many edges due to fixed threshold used : perhaps it would work if threshold is adjusted) */

#define BREAK_BULGE 0 /* 0 : Only break weaker mismatched path if coverage is weaker by MinCoverageRatio \
			 1 : Break weaker mismatched path regardless of coverage */
#define MAXBULGE 0 /* WAS 1 */ /* limit bulge length (sum of 2 paths) to MaxBulge (even if bulge has already been found) : should not be needed*/

#define BULGE_FAST BulgeFast /* Speed up graph_bulge() by skipping Dysktra shortest path check if BFSdepth <= BULGE_FAST OR paths are same length and will be merged */ 

#define CHIMTRIM_FAST 1 /* try to speed up chimtrim() by keeping track of changes in edgecnt of nodes : nodes with unchanged edgecnt for themselves and all their neighbors \
			   during the last call to chimtrim() don't need to be checked.	*/

#define CHIMTRIM_AVOID FragilePreserve /* 0 : assertion failure line 5010 */ /* avoid deleting chimeric edges near fragile sites by re-estimating edge counts after every edge deletion (cannot be repeatably multithreaded) */

#define FAST_PATH FastPaths /* Try -FastPaths 1,2 */ /* faster (slightly approximate) way to enumerate the shortest paths between pairs of nodes in a graph component in order of longest to shortest distance */

#define CC_FAST 1  /* speed up connected component traversal by creating a sublist of topnodes[] for each connected component */
#define CC_COMPRESS 1 /* WAS 0 */ /* compress SP[] array to include only the sublist of topnodes in the current connected component */
#define CYCLE_PAR 1 /* Multithreaded cycle check in graph_bulge() */
#define CYCLE_DEFER 1 /* Defer graph updates to break cycles to after multithreaded section (which is used only to locate graph bulges) : 
			 This ensures graph updates are done in the same order as with a single thread */

#define BULGE_PAR 1 /* TRY 2 */ /* Multithreaded graph bulge check in graph_bulge() (If 2 do this even if BFSdepth > BULGE_FAST : not tested) 
			       NOTE : currently causes nonrepeatable behaviour due to order in which bulges with common edges are handled (Can be fixed with BULGE_DEFER=1) */


#define BULGE_AVOID 0 /* TRY 1 */ /* If > 0 : In General bulge avoid deleting edge that is the only forward or backward edge from/to a node : 
				     instead delete another edge (even if it is not the lowest coverage edge) and preserve the lowest coverage edge with a minimum coverage value of 1.
				  */
#define BULGE_DEFER 1  /* Defer graph updates to remove General bulges to after multithreaded section (which is used only to locate graph bulges) : 
			  This ensures graph updates are done in the same order as with a single thread */

#define BFSINIT 1 /* WAS 0 */ /* Use Bread First Search (instead of DFS) to initialize ShortestPaths when there is a distance limit : theoretically faster (fewer steps), but can be slower (each step takes more time than with DFS) */

#define PVERB 1 /* ((refmap_filename || Gmap[0]->startloc > 0.0) ? 1 : 2) */ /* to display best paths set PVERB = 1, otherwise PVERB=2 */
#define PDEBUG DEBUG

#define MAXERR 4.0 /* maximum distance discrepancy as multiple of SD */
#define NMAXERR 4.0 /* used for complex bulges in place of MAXERR */

typedef short SINT;

extern double mtime();
extern double wtime();

int maxnodes = 0;
int numnodes = 0;
Cnode *node = 0;
CnodeA *nodeA = 0;/* nodeA is used in Assembler_graph() and Assembler_nodetrim() to save memory */
int *map2nodeid = 0;/* If != 0 : maps mapid to node id (if nodes are cloned, map is to first node id). Used in Ccontig.cpp */

unsigned int Maxedges = 0;
unsigned int Numedges = 0;
Cedge *edge = 0;/* edge[0..Maxedges-1] are allocated : edge[0] is not used, edge[i=1..Numedges] are used */
CedgeA *edgeA = 0;/* edgeA is used in Assembler_graph() and Assembler_nodetrim() to save memory */
int maxedges = -1;

int *nodemap = 0;/* map ids in order of location (if known, otherwise nodemap[i] == i) */
int *nodeorder = 0;/* mapping from mapid to location (if known) */

Cnodedistance *edgelist = 0;

int *EnodeidMem = 0;/* block allocation memory for node[i].Enodeid[] */
unsigned int *edgeidMem = 0;/* block allocation memory for node[i].edgeid[] */
size_t edgeMem = 0;/* size of block allocation for EnodeidMem[] and edgeidMem[] : also indicates that graph_compact() has been called (after which graph becomes more complicated) */

int numtopnodes = 0;/* number of top level nodes (all other nodes should be reachable from these nodes via Super-Edges) */
Cnode **topnodes = 0;/* array of pointers to top level nodes in the graph with Super-Edges 
			Note that if topnodes==0, all nodes in node[] are topnodes */

/* NOTE : need to update nodeorder[i] to make sure i is mapid and NOT a nodeid :

   Replace nodeorder[nodeid] by nodeorder[node[nodeid].mapid]
   Replace node[mapid] by node[nodeid]
*/

int PosEdgeInit = 0;/* If PosEdge fields are valid */

static int numthreads = -1;

typedef int intcmp (const void *,const void *);



static int qsort_nodeid = -1;

/** sort edges in increasing order of nodedistance */
static int Cnodedistinc(const Cnodedistance* p1, const Cnodedistance* p2)
{
  double distance1 = p1->nodedistance;
  double distance2 = p2->nodedistance;

  return (distance1 > distance2) ? 1 : (distance1 < distance2) ? -1 : 0;
}
/** sort edges in increasing order of nodedistance */
static int CnodedistAinc(const CnodedistanceA* p1, const CnodedistanceA* p2) // for use with CnodeA : use template
{
  double distance1 = p1->nodedistance;
  double distance2 = p2->nodedistance;

  return (distance1 > distance2) ? 1 : (distance1 < distance2) ? -1 : 0;
}

#if 0
/** sort edges in decreasing order of coverage */
static int edgecovdec(const Cedge** p1, const Cedge** p2)
{
  double cov1 = p1[0]->cov;
  double cov2 = p2[0]->cov;

  return (cov1 < cov2) ? 1 : (cov1 > cov2) ? -1 : 0;
}
#endif

/** sort map ids in increasing order of startloc + len/2 (provided startloc > -1000) */
static int startlocinc(const int *pid1, const int *pid2)
{
  Cmap *pmap1 = Gmap[*pid1];
  Cmap *pmap2 = Gmap[*pid2];
  double startloc1 = pmap1->startloc;
  double startloc2 = pmap2->startloc;
  if(startloc1 > 0.0 /* -999.9 */)
    startloc1 += pmap1->site[0][pmap1->numsite[0]+1]*0.5;
  if(startloc2 > 0.0 /* -999.9 */)
    startloc2 += pmap2->site[0][pmap2->numsite[0]+1]*0.5;
  
  return (pmap1->refid > pmap2->refid) ? 1 : (pmap1->refid < pmap2->refid) ? -1 : (startloc1 > startloc2) ? 1 : (startloc1 < startloc2) ? -1 : 0;
}

/** sort Ccontig array in decreasing order of totsites() */
static int ContigTotsitesDec(Ccontig *p1, Ccontig *p2)
{
  int totsites1 = p1->totsites();
  int totsites2 = p2->totsites();

  return (totsites1 < totsites2) ? 1 : (totsites1 > totsites2) ? -1 : 0;
}

static inline void edgereorder(Cnode *Vnode, Cnodedistance *edgelist)
{
  unsigned int *pedgeid = Vnode->edgeid;
  for(int j = Vnode->edgecnt; --j >= 0;){
    Cedge *pedge = &edge[pedgeid[j]];
    edgelist[j].pedge = pedge;
    if(pedge->Lmap == Vnode->nodeid)
      edgelist[j].nodedistance = (pedge->Lflipped ? -pedge->distance : pedge->distance);
    else {
      if(DEBUG>=2 && !(pedge->Rmap == Vnode->nodeid)){
	printf("edgereorder(Vnode->nodeid=%d):j=%d/%d:edgeid[j]=%u,Lmap=%d,Rmap=%d\n",Vnode->nodeid,j,Vnode->edgecnt,pedgeid[j],pedge->Lmap,pedge->Rmap);
	fflush(stdout);
	assert(pedge->Rmap == Vnode->nodeid);
      }
      edgelist[j].nodedistance = (pedge->Rflipped ? pedge->distance : -pedge->distance);
    }
  }

  /* sort edges by nodedistance */
  qsort(edgelist,Vnode->edgecnt,sizeof(Cnodedistance),(intcmp *)Cnodedistinc);

  /* update Vnode->edgeid[], Vnode->Enodeid[] */
  int *Enodeid = Vnode->Enodeid;
  for(int vw = Vnode->edgecnt; --vw >= 0;){
    Cedge *VWedge = edgelist[vw].pedge;
    pedgeid[vw] = VWedge->edgeid;
    if(DEBUG>=2) assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
    Enodeid[vw] = (VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap;
  }

  /* locate first node with positive nodedistance */
  Vnode->PosEdge = Vnode->edgecnt;/* in case there are no postive edges */
  for(int j=0;j < Vnode->edgecnt; j++)
    if(edgelist[j].nodedistance > 0.0){
      Vnode->PosEdge = j;
      break;
    }
}

static inline void edgereorder(CnodeA *Vnode, CnodedistanceA *edgelist) // for use with CnodeA,CedgeA : use template
{
  unsigned int *pedgeid = Vnode->edgeid;
  for(int j = Vnode->edgecnt; --j >= 0;){
    CedgeA *pedge = &edgeA[pedgeid[j]];
    edgelist[j].pedge = pedge;
    if(pedge->Lmap == Vnode->nodeid)
      edgelist[j].nodedistance = (pedge->Lflipped ? -pedge->distance : pedge->distance);
    else {
      if(DEBUG>=2) assert(pedge->Rmap == Vnode->nodeid);
      edgelist[j].nodedistance = (pedge->Rflipped ? pedge->distance : -pedge->distance);
    }
  }

  /* sort edges by nodedistance */
  qsort(edgelist,Vnode->edgecnt,sizeof(CnodedistanceA),(intcmp *)CnodedistAinc);

  /* update Vnode->edgeid[], Vnode->Enodeid[] */
  int *Enodeid = Vnode->Enodeid;
  for(int vw = Vnode->edgecnt; --vw >= 0;){
    CedgeA *VWedge = edgelist[vw].pedge;
    pedgeid[vw] = VWedge->edgeid;
    Enodeid[vw] = (VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap;
  }

  /* locate first node with positive nodedistance */
  Vnode->PosEdge = Vnode->edgecnt;/* in case there are no postive edges */
  for(int j=0;j < Vnode->edgecnt; j++)
    if(edgelist[j].nodedistance > 0.0){
      Vnode->PosEdge = j;
      break;
    }
}

/** compute nodedistance and PosEdge for all active nodes */
static void PosEdge()
{
  /* Sort Edges from most negative distance to most positive distance and remember index of first edge with positive distance */
  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {
    Cnodedistance *edgelist = new Cnodedistance[maxedges];

#pragma omp for schedule(dynamic,16)
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      if(DEBUG>=2) assert(0 <= Vnode->edgecnt && Vnode->edgecnt < maxedges);

      /* re-order edges for current node (Vnode->mapid) */
      edgereorder(Vnode,edgelist);
 
      /* compute forward and backward edge cov sum */
      Vnode->forward = 0;
      for(int vw = Vnode->PosEdge; vw < Vnode->edgecnt; vw++){    
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	Vnode->forward += VWedge->cov;
      }
      Vnode->backward = 0;
      for(int vw = 0; vw < Vnode->PosEdge; vw++){    
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	Vnode->backward += VWedge->cov;
      }
    }
    delete [] edgelist;
  }// parallel
  if(VERB>=2){
    printf("computed forward and backward edge counts\n");
    fflush(stdout);
  }
  PosEdgeInit = 1;
}

#if GRAPH_SMALL
static void PosEdgeA() /* for use with CedgeA,CnodeA */
{
  /* Sort Edges from most negative distance to most positive distance and remember index of first edge with positive distance */
#pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {
    CnodedistanceA *edgelist = new CnodedistanceA[maxedges];

#pragma omp for schedule(dynamic,16)
    for(int i=0; i < numtopnodes; i++){
      CnodeA *Vnode = (CnodeA *)topnodes[i];

      /* re-order edges for current node (Vnode->mapid) */
      edgereorder(Vnode,edgelist);
 
      /* compute forward and backward edge cov sum */
      Vnode->forward = 0;
      for(int vw = Vnode->PosEdge; vw < Vnode->edgecnt; vw++){    
	CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	Vnode->forward += VWedge->cov;
      }
      Vnode->backward = 0;
      for(int vw = 0; vw < Vnode->PosEdge; vw++){    
	CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	Vnode->backward += VWedge->cov;
      }
    }
    delete [] edgelist;
  }// parallel
  if(VERB>=2){
    printf("computed forward and backward edge counts\n");
    fflush(stdout);
  }
  PosEdgeInit = 1;
}
#endif

static int edge_delete_verb = 0;

/** mark VXedge for deletion (adding any super-Edge (chain of nodes) to side chain(s) of Vnode) */
static inline void edge_delete(Cedge *VXedge, Cnode *Vnode)
{
  if(VERB>=2 && edge_delete_verb){
    printf("edge_delete:deleting VXedge=(N%d,N%d,%0.3f%c,%d,%d),cov=%0.3f,Vnode=N%d\n",
	   nodeorder[node[VXedge->Lmap].mapid],nodeorder[node[VXedge->Rmap].mapid],VXedge->distance,VXedge->alignid ? ' ':'*',
	   VXedge->Lflipped,VXedge->Rflipped,VXedge->cov,nodeorder[Vnode->mapid]);
    fflush(stdout);
  }
  if(DEBUG>=2 && !(VXedge->cov >= -0.01)){
    printf("edge_delete:deleting VXedge=(N%d,N%d,%0.3f%c,%d,%d),cov=%0.3f,Vnode=N%d\n",
	   nodeorder[node[VXedge->Lmap].mapid],nodeorder[node[VXedge->Rmap].mapid],VXedge->distance,VXedge->alignid ? ' ':'*',
	   VXedge->Lflipped,VXedge->Rflipped,VXedge->cov,nodeorder[Vnode->mapid]);
    fflush(stdout);
    assert(VXedge->cov >= -0.01);
  }

  VXedge->valid = -1;

  if(!VXedge->alignid && VXedge->cov > 0.0){/* Super-Edge */
    if(DEBUG>=2) assert(VXedge->Lchainid && VXedge->Rchainid);
    if(VXedge->Lmap == Vnode->nodeid){/* Add Edge VXedge->Lchainid to Vnode->chain (list of side chains) */
      if(DEBUG>=2) assert(1 <= VXedge->Lchainid && VXedge->Lchainid <= Numedges);
      Vnode->chain.push_back(VXedge->Lchainid);
    } else {/* Add Edge VXedge->Rchain to Vnode->chain (list of side chains) */
      if(DEBUG>=2) assert(VXedge->Rmap == Vnode->nodeid);
      if(DEBUG>=2) assert(1 <= VXedge->Rchainid && VXedge->Rchainid <= Numedges);
      Vnode->chain.push_back(VXedge->Rchainid);
    }
    /* NOTE : no need to add restored chain nodes to topnodes[] since they are only needed during construction of consensus map */
  }
}

/** mark VXedge for deletion (adding any super-Edge (chain of nodes) to side chain(s) of Vnode) */
static inline void edge_delete(CedgeA *VXedge, CnodeA *Vnode) // for use with CedgeA, CnodeA : use template
{
  if(VERB>=2 && edge_delete_verb){
    printf("edge_delete:deleting VXedge=(N%d,N%d,%0.3f%c,%d,%d),cov=%0.3f,Vnode=N%d\n",
	   nodeorder[node[VXedge->Lmap].mapid],nodeorder[node[VXedge->Rmap].mapid],VXedge->distance,VXedge->alignid ? ' ':'*',
	   VXedge->Lflipped,VXedge->Rflipped,VXedge->cov,nodeorder[Vnode->mapid]);
    fflush(stdout);
  }
  if(DEBUG>=2 && !(VXedge->cov >= -0.01)){
    printf("edge_delete:deleting VXedge=(N%d,N%d,%0.3f%c,%d,%d),cov=%0.3f,Vnode=N%d\n",
	   nodeorder[node[VXedge->Lmap].mapid],nodeorder[node[VXedge->Rmap].mapid],VXedge->distance,VXedge->alignid ? ' ':'*',
	   VXedge->Lflipped,VXedge->Rflipped,VXedge->cov,nodeorder[Vnode->mapid]);
    fflush(stdout);
    assert(VXedge->cov >= -0.01);
  }

  VXedge->valid = -1;

  if(DEBUG>=2) assert(VXedge->alignid);
}

/** just count number of internal nodes (maps) in super-edge (return 0 for regular edge) */

static int edgenodecnt(Cedge *Nedge, Cnode *Lnode, int Lflipped)
{
  int mapcnt = 0;
  if(DEBUG>=2) assert(Nedge->alignid >= 0);

  if(Nedge->alignid)
    return mapcnt;

  if(DEBUG>=2) assert(Nedge->valid >= 0);
  if(DEBUG>=2) assert(Lnode->nodeid == Lnode - node && 0 <= Lnode->nodeid && Lnode->nodeid < numnodes);

  if(Nedge->Lmap == Lnode->nodeid){/* go from Lmap to Rmap */
    Cedge *LRedge = &edge[Nedge->Lchainid];

    if(DEBUG>=2) assert(LRedge->edgeid == LRedge - edge && 1 <= LRedge->edgeid && LRedge->edgeid <= Numedges);
    if(DEBUG>=2) assert(LRedge->Lmap == Lnode->nodeid || LRedge->Rmap == Lnode->nodeid);

    Cnode *Rnode = &node[(LRedge->Lmap == Lnode->nodeid) ? LRedge->Rmap : LRedge->Lmap];

    if(DEBUG>=2) assert(Rnode->nodeid == Rnode - node && 0 <= Rnode->nodeid && Rnode->nodeid < numnodes);

    int Rflipped = Lflipped ^ (LRedge->Lflipped | LRedge->Rflipped);

    while(Rnode->nodeid != Nedge->Rmap){
      mapcnt += edgenodecnt(LRedge,Lnode,Lflipped);
      mapcnt++;
      Lnode = Rnode;
      Lflipped = Rflipped;

      if(DEBUG>=2) assert(Lnode->edgecnt == 2);
      if(DEBUG>=2) assert(&edge[Lnode->edgeid[0]] == LRedge || &edge[Lnode->edgeid[1]] == LRedge);

      LRedge = (&edge[Lnode->edgeid[0]] == LRedge ? &edge[Lnode->edgeid[1]] : &edge[Lnode->edgeid[0]]);

      if(DEBUG>=2) assert(LRedge->edgeid == LRedge - edge && 1 <= LRedge->edgeid && LRedge->edgeid <= Numedges);
      if(DEBUG>=2) assert(LRedge->Lmap == Lnode->nodeid || LRedge->Rmap == Lnode->nodeid);

      int Rnodeid = (LRedge->Lmap == Lnode->nodeid) ? LRedge->Rmap : LRedge->Lmap;

      if(DEBUG>=2) assert(0 <= Rnodeid && Rnodeid < numnodes);

      Rnode = &node[Rnodeid];

      if(DEBUG>=2) assert(Rnode->nodeid == Rnodeid && Rnode->nodeid == Rnode - node && 0 <= Rnode->nodeid && Rnode->nodeid < numnodes);

      Rflipped = Lflipped ^ (LRedge->Lflipped | LRedge->Rflipped);
    }
    mapcnt += edgenodecnt(LRedge,Lnode,Lflipped);
  } else { /* go from Rmap to Lmap */

    if(DEBUG>=2) assert(Nedge->Rmap == Lnode->nodeid);

    Cedge *LRedge = &edge[Nedge->Rchainid];

    if(DEBUG>=2) assert(LRedge->edgeid == LRedge - edge && 1 <= LRedge->edgeid && LRedge->edgeid <= Numedges);
    if(DEBUG>=2) assert(LRedge->Lmap == Lnode->nodeid || LRedge->Rmap == Lnode->nodeid);

    Cnode *Rnode = &node[(LRedge->Lmap == Lnode->nodeid) ? LRedge->Rmap : LRedge->Lmap];

    if(DEBUG>=2) assert(Rnode->nodeid == Rnode - node && 0 <= Rnode->nodeid && Rnode->nodeid < numnodes);

    int Rflipped = Lflipped ^ (LRedge->Lflipped | LRedge->Rflipped);

    while(Rnode->nodeid != Nedge->Lmap){
      mapcnt += edgenodecnt(LRedge,Lnode,Lflipped);
      mapcnt++;
      Lnode = Rnode;
      Lflipped = Rflipped;

      if(DEBUG>=2) assert(Lnode->edgecnt == 2);
      if(DEBUG>=2) assert(&edge[Lnode->edgeid[0]] == LRedge || &edge[Lnode->edgeid[1]] == LRedge);

      LRedge = (&edge[Lnode->edgeid[0]] == LRedge ? &edge[Lnode->edgeid[1]] : &edge[Lnode->edgeid[0]]);

      if(DEBUG>=2) assert(LRedge->edgeid == LRedge - edge && 1 <= LRedge->edgeid && LRedge->edgeid <= Numedges);
      if(DEBUG>=2) assert(LRedge->Lmap == Lnode->nodeid || LRedge->Rmap == Lnode->nodeid);

      int Rnodeid = (LRedge->Lmap == Lnode->nodeid) ? LRedge->Rmap : LRedge->Lmap;

      if(DEBUG>=2) assert(0 <= Rnodeid && Rnodeid < numnodes);

      Rnode = &node[Rnodeid];

      if(DEBUG>=2) assert(Rnode->nodeid == Rnodeid && Rnode->nodeid == Rnode - node && 0 <= Rnode->nodeid && Rnode->nodeid < numnodes);

      Rflipped = Lflipped ^ (LRedge->Lflipped | LRedge->Rflipped);
    }
    mapcnt += edgenodecnt(LRedge,Lnode,Lflipped);
  }
  return mapcnt;
}

/** assemble (and display) path underlying super-edge on a single line starting at the end with node Lnode (but don't display nodes corresponding to Lmap or Rmap)
 returns number of nodes embedded in super-edge */

static int edgedisplay(Cedge *pedge, Cnode *Lnode, int Lflipped,
		       Ccontigpath *path, int cnt, int verb)
{
  int mapcnt = 0;

  if(DEBUG>=2) assert(pedge->alignid >= 0);

  if(pedge->alignid){
    if(verb >= PVERB)
      printf("  (D=%0.3f,C=%0.1f) ", pedge->distance, pedge->cov);
    if(path)
      path->align[cnt+mapcnt] = alignment[pedge->alignid - 1];
    return mapcnt;
  }

  if(DEBUG>=2) assert(pedge->valid >= 0);

  if(pedge->Lmap == Lnode->nodeid){/* go from Lmap to Rmap */
    Cedge *LRedge = &edge[pedge->Lchainid];
    Cnode *Rnode = &node[(LRedge->Lmap == Lnode->nodeid) ? LRedge->Rmap : LRedge->Lmap];
    int Rflipped = Lflipped ^ (LRedge->Lflipped | LRedge->Rflipped);
    while(Rnode->nodeid != pedge->Rmap){
      mapcnt += edgedisplay(LRedge,Lnode,Lflipped,path,cnt+mapcnt,verb);
      mapcnt++;
      if(verb >= PVERB)
	printf("N%d%c %lld\n", nodeorder[Rnode->mapid],Rflipped ? '\'':' ', Gmap[Rnode->mapid]->id);
      if(path){
	path->nodeid[cnt+mapcnt] = Rnode->nodeid;
	path->flip[cnt+mapcnt] = Rflipped;
      }
      Lnode = Rnode;
      Lflipped = Rflipped;
      if(DEBUG>=2) assert(Lnode->edgecnt == 2);
      LRedge = (&edge[Lnode->edgeid[0]] == LRedge ? &edge[Lnode->edgeid[1]] : &edge[Lnode->edgeid[0]]);
      Rnode = &node[(LRedge->Lmap == Lnode->nodeid) ? LRedge->Rmap : LRedge->Lmap];
      Rflipped = Lflipped ^ (LRedge->Lflipped | LRedge->Rflipped);
    }
    mapcnt += edgedisplay(LRedge,Lnode,Lflipped,path,cnt+mapcnt,verb);
  } else { /* go from Rmap to Lmap */
    if(DEBUG>=2) assert(pedge->Rmap == Lnode->nodeid);
    Cedge *LRedge = &edge[pedge->Rchainid];
    Cnode *Rnode = &node[(LRedge->Lmap == Lnode->nodeid) ? LRedge->Rmap : LRedge->Lmap];
    int Rflipped = Lflipped ^ (LRedge->Lflipped | LRedge->Rflipped);
    while(Rnode->nodeid != pedge->Lmap){
      mapcnt += edgedisplay(LRedge,Lnode,Lflipped,path,cnt+mapcnt,verb);
      mapcnt++;
      if(verb >= PVERB)
	printf("N%d%c %lld\n", nodeorder[Rnode->mapid],Rflipped ? '\'':' ', Gmap[Rnode->mapid]->id);
      if(path){
	path->nodeid[cnt+mapcnt] = Rnode->nodeid;
	path->flip[cnt+mapcnt] = Rflipped;
      }
      Lnode = Rnode;
      Lflipped = Rflipped;
      if(DEBUG>=2) assert(Lnode->edgecnt == 2);
      LRedge = (&edge[Lnode->edgeid[0]] == LRedge ? &edge[Lnode->edgeid[1]] : &edge[Lnode->edgeid[0]]);
      Rnode = &node[(LRedge->Lmap == Lnode->nodeid) ? LRedge->Rmap : LRedge->Lmap];
      Rflipped = Lflipped ^ (LRedge->Lflipped | LRedge->Rflipped);
    }
    mapcnt += edgedisplay(LRedge,Lnode,Lflipped,path,cnt+mapcnt,verb);
  }
  return mapcnt;
}

/* remove edges with valid <= -1 (or NULL edgeids) from edge lists of all current nodes
   returns number of edges removed */

static long long GCedges(unsigned int &edgecnt, Cedge* &edge)
{
  long long origcnt = 0;
  long long finalcnt = 0;

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {
    long long myorigcnt = 0;
    long long myfinalcnt = 0;

    #pragma omp for schedule(dynamic,16) nowait
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      myorigcnt += Vnode->edgecnt;
      int k = 0;
      for(int vw = 0; vw < Vnode->edgecnt;vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	if(!VWedgeid)
	  continue;
	Cedge *VWedge = &edge[VWedgeid];
	if(VWedge->valid <= -1)
	  continue;
	if(DEBUG>=2 && !(VWedge->valid == 1)){
	  printf("GCedges():i=%d,Vnode=N%d,vw=%d,VWedge=(N%d,N%d,%0.3f,%d,%d):valid=%d\n",
		 i,nodeorder[Vnode->mapid],vw,nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],
		 VWedge->distance,VWedge->Lflipped,VWedge->Rflipped,VWedge->valid);
	  assert(VWedge->valid == 1);
	}
	Vnode->Enodeid[k] = Vnode->Enodeid[vw];
	Vnode->edgeid[k++] = VWedgeid;
      }
      if(VERB>=2 && Vnode->nodeid == 0){
	printf("topnodes[%d]: nodeid=%d, edgecnt=%d -> %d\n",i,Vnode->nodeid,Vnode->edgecnt, k);
	fflush(stdout);
      }
      Vnode->edgecnt = k;
      myfinalcnt += k;
    }// omp for

    #pragma omp atomic
    origcnt += myorigcnt;

    #pragma omp atomic
    finalcnt += myfinalcnt;
  } // omp parallel

  PosEdgeInit = 0;

  if(DEBUG>=2 && !((origcnt % 2) == 0 && (finalcnt % 2) == 0)){
    unsigned int edgecnt = 0;
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edge[i].valid < 0)
	continue;
      assert(edge[i].mark == 0);
      edgecnt++;
    }

    printf("GCedges: origcnt= %lld, finalcnt= %lld, edgecnt= %u\n",origcnt,finalcnt,edgecnt);
    fflush(stdout);


    assert((origcnt % 2) == 0 && (finalcnt % 2) == 0);
  }

  if(VERB>=2 && origcnt > finalcnt){
    printf("GCedges: numedges = %lld -> %lld: wtime= %0.6f\n",origcnt/2,finalcnt/2, wtime());
    fflush(stdout);
  }

  edgecnt = finalcnt/2;

  return origcnt-finalcnt;
}

#if GRAPH_SMALL
static long long GCedges(unsigned int &edgecnt, CedgeA* &edgeA) // for use with CedgeA,CnodeA : use template
{

  /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
  if(DEBUG>=2){
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edgeA[i].valid < 0)
	continue;
      assert(edgeA[i].mark == 0);
    }
    size_t edgecnt = 0;
    for(int i=0; i < numtopnodes; i++){
      CnodeA *Vnode = (CnodeA *) topnodes[i];
      if(Vnode->mark >= 2){
	assert(Vnode->edgecnt <= 0);
	continue;
      }
      edgecnt += Vnode->edgecnt;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	if(!(1 <= VWedgeid && VWedgeid <= Numedges)){
	  printf("i=%d/%d:node[i]:mapid=%d,nodeid=%d,mark=%d,edgecnt=%d/%d,vw=%d:node[i].edgeid[vw]=%u,VWedgeid=%u,Numedges=%u\n",
		 i,numnodes,Vnode->mapid,Vnode->nodeid,Vnode->mark,Vnode->edgecnt,Vnode->maxedge,vw,Vnode->edgeid[vw],VWedgeid,Numedges);
	  fflush(stdout);
	  assert(1 <= VWedgeid && VWedgeid <= Numedges);
	}
	assert(edge[VWedgeid].mark == 0);
      }
    }    
    if(DEBUG){
      for(int i=0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];

	  VWedge->mark++;

	  if(VWedge->valid <= 0)
	    continue;
	  assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
	  assert(VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap));
	}
      }    
      for(int i=0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->mark != 2){
	    printf("i=%d/%d: mapid=%d: vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d\n",
		   i,numnodes,Vnode->mapid, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark, VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped);
	    fflush(stdout);
	    assert(VWedge->mark == 2);
	  }
	}
      }    
      /* reset edge mark values to 0 */
      for(int i=0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edge[Vnode->edgeid[vw]].mark = 0;
      }    
    }
    if((edgecnt % 2) != 0){
      printf("GCedges:sum of node[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",edgecnt,numnodes,Numedges);
      fflush(stdout);
      assert((edgecnt%2) == 0);
    }
  }

  long long origcnt = 0;
  long long finalcnt = 0;

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {
    long long myorigcnt = 0;
    long long myfinalcnt = 0;

    #pragma omp for schedule(dynamic,16) nowait
    for(int i=0; i < numtopnodes; i++){
      CnodeA *Vnode = (CnodeA *) topnodes[i];
      myorigcnt += Vnode->edgecnt;
      int k = 0;
      for(int vw = 0; vw < Vnode->edgecnt;vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	if(!VWedgeid)
	  continue;
	CedgeA *VWedge = &edgeA[VWedgeid];
	if(VWedge->valid <= -1)
	  continue;
	if(DEBUG>=2 && !(VWedge->valid == 1)){
	  printf("GCedges():i=%d,Vnode=N%d,vw=%d,VWedge=(N%d,N%d,%0.3f,%d,%d):valid=%d\n",
		 i,nodeorder[Vnode->mapid],vw,nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],
		 VWedge->distance,VWedge->Lflipped,VWedge->Rflipped,VWedge->valid);
	  assert(VWedge->valid == 1);
	}
	Vnode->Enodeid[k] = Vnode->Enodeid[vw];
	Vnode->edgeid[k++] = VWedgeid;
      }
      if(VERB>=2 && Vnode->nodeid == 0){
	printf("topnodes[%d]: nodeid=%d, edgecnt=%d -> %d\n",i,Vnode->nodeid,Vnode->edgecnt, k);
	fflush(stdout);
      }
      Vnode->edgecnt = k;
      myfinalcnt += k;
    }// omp for

    #pragma omp atomic
    origcnt += myorigcnt;

    #pragma omp atomic
    finalcnt += myfinalcnt;
  } // omp parallel

  PosEdgeInit = 0;

  /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
  if(DEBUG>=2){
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edgeA[i].valid < 0)
	continue;
      assert(edgeA[i].mark == 0);
    }
    size_t edgecnt = 0;
    for(int i=0; i < numtopnodes; i++){
      CnodeA *Vnode = (CnodeA *) topnodes[i];
      if(Vnode->mark >= 2){
	assert(Vnode->edgecnt <= 0);
	continue;
      }
      edgecnt += Vnode->edgecnt;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	if(!(1 <= VWedgeid && VWedgeid <= Numedges)){
	  printf("i=%d/%d:node[i]:mapid=%d,nodeid=%d,mark=%d,edgecnt=%d/%d,vw=%d:node[i].edgeid[vw]=%u,VWedgeid=%u,Numedges=%u\n",
		 i,numnodes,Vnode->mapid,Vnode->nodeid,Vnode->mark,Vnode->edgecnt,Vnode->maxedge,vw,Vnode->edgeid[vw],VWedgeid,Numedges);
	  fflush(stdout);
	  assert(1 <= VWedgeid && VWedgeid <= Numedges);
	}
	assert(edge[VWedgeid].mark == 0);
      }
    }    
    if(DEBUG){
      for(int i=0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];

	  VWedge->mark++;

	  if(VWedge->valid <= 0)
	    continue;
	  assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
	  assert(VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap));
	}
      }    
      for(int i=0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->mark != 2){
	    printf("i=%d/%d: mapid=%d: vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d\n",
		   i,numnodes,Vnode->mapid, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark, VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped);
	    fflush(stdout);
	    assert(VWedge->mark == 2);
	  }
	}
      }    
      /* reset edge mark values to 0 */
      for(int i=0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edge[Vnode->edgeid[vw]].mark = 0;
      }    
    }
    if((edgecnt % 2) != 0){
      printf("GCedges:sum of node[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",edgecnt,numnodes,Numedges);
      fflush(stdout);
      assert((edgecnt%2) == 0);
    }
  }

  if(DEBUG && !((origcnt % 2) == 0 && (finalcnt % 2) == 0)){
    unsigned int edgecnt = 0;
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edgeA[i].valid < 0)
	continue;
      assert(edgeA[i].mark == 0);
      edgecnt++;
    }

    printf("GCedges: origcnt= %lld, finalcnt= %lld, edgecnt= %u\n",origcnt,finalcnt,edgecnt);
    fflush(stdout);

    assert((origcnt % 2) == 0 && (finalcnt % 2) == 0);
  }

  if(VERB>=2 && origcnt > finalcnt){
    printf("GCedges: numedges = %lld -> %lld: wtime= %0.6f\n",origcnt/2,finalcnt/2, wtime());
    fflush(stdout);
  }

  edgecnt = finalcnt/2;

  return origcnt-finalcnt;
}
#endif

/** remove nodes with mark==2 from topnodes[]. 
   Also remove nodes with no edges
   return count of nodes removed */
static int GCnodes()
{
  int k=0;
  for(int i=0; i < numtopnodes; i++){
    Cnode *Vnode = topnodes[i];
    if(VERB>=2 && (Vnode->mark >= 2 || Vnode->edgecnt <= 0) && i==0){
      printf("GCnodes():topnodes[%d]: mark=%d, edgecnt=%d : deleting node\n",i,Vnode->mark,Vnode->edgecnt);
      fflush(stdout);
    }
    if(VERB >= 2 && Vnode->edgecnt <= 0)
      printf("Deleted isolated node N%d\n",nodeorder[Vnode->mapid]);
    if(Vnode->mark >= 2){
      //      if(DEBUG>=2) assert(Vnode->edgecnt <= 0);
      continue;
    }
    if(Vnode->edgecnt <= 0){
      Vnode->mark = 2;
      continue;
    }
    topnodes[k++] = topnodes[i];
  }
  int delta = numtopnodes - k;
  if(VERB>=2 && delta > 0){
    printf("GCnodes: numtopnodes = %d -> %d (delta=%d): wtime= %0.6f\n",numtopnodes, k, delta, wtime());
    fflush(stdout);
  }
  numtopnodes = k;

  if(DEBUG>=3 && edgeMem){
    for(int i = 0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(VWedge->valid <= 0)
	  continue;
	if(!((VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid) && (VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap)))){
          printf("i=%d/%d:nodeid=%d(N%d), vw=%d/%d: valid=%d,edgeid=%u,Lmap=%d(N%d),Rmap=%d(N%d),node[i].Enodeid[vw]= %d\n",
                 i,numnodes,Vnode->nodeid,Vnode->mapid,vw,Vnode->edgeid[vw],Vnode->edgecnt,VWedge->valid,VWedge->Lmap,node[VWedge->Lmap].mapid,VWedge->Rmap,node[VWedge->Rmap].mapid,Vnode->Enodeid[vw]);
	  fflush(stdout);
          assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
	  assert(VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap));
        }
      }
    }    
  }

  /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
  if(DEBUG>=3 && !edgeMem){
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edge[i].valid < 0)
	continue;
      assert(edge[i].mark == 0);
    }
    size_t edgecnt = 0;
    for(int i = 0; i < numnodes; i++){
      Cnode *Vnode = &node[i];
      if(Vnode->mark >= 2){
	assert(Vnode->edgecnt <= 0);
	continue;
      }
      edgecnt += Vnode->edgecnt;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	if(!(1 <= VWedgeid && VWedgeid <= Numedges)){
	  printf("i=%d/%d:node[i]:mapid=%d,nodeid=%d,mark=%d,edgecnt=%d/%d,vw=%d:node[i].edgeid[vw]=%u,VWedgeid=%u,Numedges=%u\n",
		 i,numnodes,Vnode->mapid,Vnode->nodeid,Vnode->mark,Vnode->edgecnt,Vnode->maxedge,vw,Vnode->edgeid[vw],VWedgeid,Numedges);
	  fflush(stdout);
	  assert(1 <= VWedgeid && VWedgeid <= Numedges);
	}
	assert(edge[VWedgeid].mark == 0);
      }
    }    
    if(DEBUG){
      for(int i = 0; i < numnodes; i++){
	Cnode *Vnode = &node[i];
	assert(Vnode->nodeid == i);
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  VWedge->mark++;

	  assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
	  assert(VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap));
	}
      }    
      for(int i = 0; i < numnodes; i++){
	Cnode *Vnode = &node[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->mark != 2){
	    printf("i=%d/%d: mapid=%d : vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d\n",
		   i,numnodes,Vnode->mapid, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark, VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped);
	    fflush(stdout);
	    assert(VWedge->mark == 2);
	  }
	}
      }    
      for(int i = 0; i < numnodes; i++){
	Cnode *Vnode = &node[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edge[Vnode->edgeid[vw]].mark = 0;
      }    
      /* check that topnodes[0..numtopnodes-1] reach all nodeA[] that do not have mark >= 2 */
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	assert(Vnode->mark == 0);
	Vnode->mark = 2;
      }
      for(int i = 0; i < numnodes; i++)
	assert(node[i].mark >= 2);
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	Vnode->mark = 0;
      }
    }
    if((edgecnt % 2) != 0){
      printf("GCnodes():sum of node[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",edgecnt,numnodes,Numedges);
      fflush(stdout);
      assert((edgecnt%2) == 0);
    }
  }

  return delta;
}

#if GRAPH_SMALL
static int GCnodesA()/* for use with CnodeA */
{
  int k=0;
  for(int i=0; i < numtopnodes; i++){
    CnodeA *Vnode = (CnodeA *) topnodes[i];
    if(VERB>=2 && (Vnode->mark >= 2 || Vnode->edgecnt <= 0) && i==0){
      printf("GCnodesA():topnodes[%d]: mark=%d, edgecnt=%d : deleting node\n",i,Vnode->mark,Vnode->edgecnt);
      fflush(stdout);
    }
    if(VERB >= PVERB+1 && Vnode->edgecnt <= 0)
      printf("Deleted isolated node N%d\n",nodeorder[Vnode->mapid]);
    if(Vnode->mark >= 2){
      if(DEBUG>=2) assert(Vnode->edgecnt <= 0);
      continue;
    }
    if(Vnode->edgecnt <= 0){
      Vnode->mark = 2;
      continue;
    }
    topnodes[k++] = topnodes[i];
  }
  int delta = numtopnodes - k;
  if(VERB>=2 && delta > 0){
    printf("GCnodesA: numtopnodes = %d -> %d (delta=%d): wtime= %0.6f\n",numtopnodes, k, delta, wtime());
    fflush(stdout);
  }
  numtopnodes = k;

  /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
  if(DEBUG>=3){
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edgeA[i].valid < 0)
	continue;
      assert(edgeA[i].mark == 0);
    }
    size_t edgecnt = 0;
    for(int i = 0; i < numnodes; i++){
      CnodeA *Vnode = &nodeA[i];
      if(Vnode->mark >= 2){
	assert(Vnode->edgecnt <= 0);
	continue;
      }
      edgecnt += Vnode->edgecnt;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	assert(1 <= VWedgeid && VWedgeid <= Numedges);
	assert(edgeA[VWedgeid].mark == 0);
      }
    }    
    if(DEBUG){
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edgeA[Vnode->edgeid[vw]].mark++;
      }    
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	  if(VWedge->mark != 2){
	    printf("i=%d/%d: mapid=%d: vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d:align->logPV= %0.2f,or=%d\n",
		   i,numnodes,Vnode->mapid, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark,
		   VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped,alignment[VWedge->alignid - 1]->logPV, alignment[VWedge->alignid - 1]->orientation);
	    fflush(stdout);
	    assert(VWedge->mark == 2);
	  }
	}
      }    
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edgeA[Vnode->edgeid[vw]].mark = 0;
      }    
      /* check that topnodes[0..numtopnodes-1] reach all nodeA[] that do not have mark >= 2 */
      for(int i = 0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	assert(Vnode->mark == 0);
	Vnode->mark = 2;
      }
      for(int i = 0; i < numnodes; i++)
	assert(nodeA[i].mark >= 2);
      for(int i = 0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	Vnode->mark = 0;
      }
    }
    if((edgecnt % 2) != 0){
      printf("GCnodesA():sum of nodeA[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",edgecnt,numnodes,Numedges);
      fflush(stdout);
      assert((edgecnt%2) == 0);
    }
  }

  return delta;
}
#endif

#if GRAPH_SMALL==0
#define CnodeA Cnode
#define CedgeA Cedge
#define nodeA node
#define edgeA edge

#define GCnodesA GCnodes
#define PosEdgeA PosEdge
#endif

/** sort edges in decreasing order of alignment logPV */
static int edgePVdec(const unsigned int* p1, const unsigned int* p2)
{
  double logPV1 = alignment[edgeA[p1[0]].alignid - 1]->logPV;
  double logPV2 = alignment[edgeA[p2[0]].alignid - 1]->logPV;

  return (logPV1 < logPV2) ? 1 : (logPV1 > logPV2) ? -1 : 0;
}

/** sort edges in increasing order of nodeid */
static int edgeidnodeinc(const unsigned int* p1, const unsigned int* p2)
{
  CedgeA *edge1 = &edgeA[*p1];
  CedgeA *edge2 = &edgeA[*p2];
  int nodeid1 = (edge1->Lmap == qsort_nodeid) ? edge1->Rmap : edge1->Lmap;
  int nodeid2 = (edge2->Lmap == qsort_nodeid) ? edge2->Rmap : edge2->Lmap;

  return nodeid1 - nodeid2;
}

/** Write all nodes and edges (passing filter) to the graph object. Calculate center-to-center
offsets for edges. Link nodes by edges.
 */
void graph_build(Cmap **map, int nummaps, Calign **alignment, size_t &numaligns)
{
  if(VERB){
    printf("Completed Data input,PVERB=%d(Gmap[0]->id=%lld,mapid=%d,startloc=%0.1f,endloc=%0.1f):time=%0.6f(total=%0.6f)\n",PVERB,Gmap[0]->id,Gmap[0]->mapid,Gmap[0]->startloc,Gmap[0]->endloc,wtime(),mtime());
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  if(colors != 1){
    printf("graph_build not implemented for colors=%d\n",colors);
    exit(1);
  }

  numthreads = 1;
  #ifdef _OPENMP
  numthreads = omp_get_max_threads();
  if(DEBUG) assert(numthreads >= 1);
  if(MaxThreads > 0)
    numthreads = min(numthreads,MaxThreads);
  #endif

  if(nummaps >= (int)MASK(31) || nummaps < 0 || numaligns >= MASK(31) /* || numaligns < 0 */){
    printf("graph_build cannot support more than 2 billion maps or alignments : nummaps=%d,numaligns=%lu\n",
	   nummaps,numaligns);
    exit(1);
  }

  nodeorder = new int[nummaps];
  nodemap = new int[nummaps];
  for(int i = nummaps; --i >= 0;)
    nodemap[i] = i;
  if(Gmap[0]->startloc)  /* sort mapid's in nodemap[0..nummaps-1] in ascending order of startloc */
    qsort(nodemap,nummaps,sizeof(int),(intcmp *)startlocinc);
    
  if(VERB/* HERE HERE */){
    printf("Initialized and sorted %d nodes:time=%0.6f\n",nummaps,wtime());
    fflush(stdout);
  }

  if(DEBUG>=2)
    for(int i = nummaps; --i >= 0;)
      nodeorder[i] = -1;
  for(int i = nummaps; --i >= 0;)
    nodeorder[nodemap[i]] = i;
  if(DEBUG>=2)
    for(int i = nummaps; --i >= 0;)
      assert(nodeorder[i] >= 0 && nodeorder[i] < nummaps);
  if(VERB>=PVERB && Gmap[0]->startloc){
    printf("%d maps in order of startloc:\n",nummaps);
    for(int i = 0; i < nummaps;i++)
      printf("N%d: mapid=%d : id=%lld, startloc=%0.3f, endloc=%0.3f, flipped=%d\n",
	     i, nodemap[i], Gmap[nodemap[i]]->id, Gmap[nodemap[i]]->startloc, Gmap[nodemap[i]]->endloc, Gmap[nodemap[i]]->flipped);
    fflush(stdout);
  }

  /* create one graph node for each map, so that initially Gmap[] index is same as node[] index (this can change if CROSSBREAK>=1) */
  numnodes = nummaps;
  maxnodes = numnodes;/* NOTE : If node cloning is implemented, maxnodes may need to be set larger than numnodes */

  nodeA = new CnodeA[maxnodes];
  CnodeA *pnode = nodeA;

  for(int i = 0; i < nummaps; i++, pnode++){
    pnode->mapid = pnode->nodeid = i;
    pnode->cov = 1;
    pnode->forward = pnode->backward = 0;
    pnode->len = Gmap[i]->site[0][Gmap[i]->numsite[0]+1];
  }

  /* Now add one edge to the graph for each alignment */
  if(numaligns > MASK(31)){
    printf("Too many alignments = %lld : Cannot exceed 31 bits (%d)\n",(long long)numaligns,MASK(31));
    fflush(stdout);
  }
  for(size_t i = 0; i < numaligns; i++){
    if(DEBUG && alignment[i]==0){
      printf("i=%lu,numaligns=%lu,alignment[i]=0 !\n",i,numaligns);
      fflush(stdout);
      assert(alignment[i] != 0);
    }
    if(DEBUG>=2) assert(alignment[i]->mapid1 != alignment[i]->mapid2);
    //    alignment[i]->chimpair = i;/* used as index into alignment[] : helps with garbage collection */
  }

  if(numaligns + 1024 >= MASK_LL(32)){
    printf("numaligns = %lu : Number of alignments cannot exeed %lld\n", numaligns, MASK_LL(32) - 1024);
    fflush(stdout); exit(1);
  }
  Numedges = numaligns;
  maxedgealloc((unsigned int)(numaligns + 1024), Maxedges, edgeA);

  if(VERB/* HERE HERE >=2 */){
    printf("edgeA= %p, numaligns= %lu, Numedges=%u, Maxedges=%u\n",edgeA, numaligns, Numedges, Maxedges);
    fflush(stdout);
  }

  if(DEBUG>=2)
    for(unsigned int i = 0; i < Maxedges; i++)
      if(!(edgeA[i].mark == 0)){
	printf("edge=%p, &edge[i=%d]=%p, edge[i].mark=%d, sizeof(Cedge)= %lu\n",edgeA,i,&edgeA[i],edgeA[i].mark,sizeof(Cedge));
	fflush(stdout);
	assert(edgeA[i].mark == 0);
      }

  double CovSum = 0.0, MaxCov = 0.0;
  int CovCnt = 0;

  CedgeA *pedge = &edgeA[1];// NOTE : edgeA[0] is not used, so edgeid == 0 can be used as NULL index
  for(size_t i = 0; i < numaligns; i++, pedge++){
    pedge->edgeid = i+1;

    Calign *palign = alignment[i];
    if(DEBUG>=2) assert(palign->mapid1 != palign->mapid2);

    pedge->alignid = i+1;
    pedge->valid = 1;
    pedge->cov = 1.0;
    if(GraphWeight > 0.0){/* adjust initial pedge->cov based on align->logPV */
      double X = palign->logPV - LogPvThreshold;
      if(GraphWeight == 1.0)
	pedge->cov = X;
      else
	pedge->cov = (pow(GraphWeight,X)-1)/(GraphWeight-1.0);
      pedge->cov = min(GraphMaxWeight, pedge->cov);

      if(VERB){
	MaxCov = max(MaxCov,pedge->cov);
	CovSum += pedge->cov;
	CovCnt++;
      }
    }

    //    pedge->ucnt = 0; // now initialized in Cedge constructor
    
    /* compute distance of center of mapid2 from center of mapid1 */
    Cmap *pmap1 = Gmap[palign->mapid1];
    Cmap *pmap2 = Gmap[palign->mapid2];

    int N = pmap1->numsite[0];
    int M = pmap2->numsite[0];
    int K = palign->numpairs;

    FLOAT *Y = pmap1->site[0];
    FLOAT *X = pmap2->site[0];

    double centerY = Y[N+1]*0.5;
    double centerX = X[M+1]*0.5;

    double distanceYtoX;

    if(1){/* average over all aligned sites, for more robust estimate of center to center distance */
      distanceYtoX = 0.0;
      for(int index = 0; index < K; index++){
	int I = palign->sites1[index];
	int J = palign->sites2[index];
	double dist = (Y[I]-centerY) + (centerX-(palign->orientation ? (X[M+1]-X[M+1-J]) : X[J]));
	distanceYtoX += dist;
      }
      distanceYtoX /= K;
    } else {/* old less robust method */
      /* locate aligned site just to the right of each center point */
      int index1;
      for(index1 = 0;index1 < K;index1++){
	int i = palign->sites1[index1];
	if(Y[i] > centerY)
	  break;
      } 
      if(DEBUG>=2) assert(index1 >= K || Y[palign->sites1[index1]] > centerY);
      if(DEBUG>=2) assert(index1 <= 0 || Y[palign->sites1[index1-1]] <= centerY);

      int index2 = 0;
      for(;index2 < K; index2++){
	int i = palign->sites2[index2];
	if((palign->orientation ? (X[M+1] - X[M+1-i]) : X[i]) > centerX)
	  break;
      }
      if(DEBUG>=2) assert(index2 >= K || (palign->orientation ? X[M+1]-X[M+1-palign->sites2[index2]] : X[palign->sites2[index2]]) > centerX);
      if(DEBUG>=2 && !(index2 <= 0 || (palign->orientation ? X[M+1] - X[M+1-palign->sites2[index2-1]] : X[palign->sites2[index2-1]]) <= centerX)){
	printf("alignment[%lu]:Yid=%d,Xid=%d,N=%d,M=%d,centerY=%0.3f,centerX=%0.3f:K=%d,index2=%d\n",
	       i,palign->mapid1,palign->mapid2,N,M,centerY,centerX,K,index2);
	if(index2 > 0)
	  printf("   sites1[index2-1]=%d,sites2[index2-1]=%d:X[%d]=%0.3f\n",
		 palign->sites1[index2-1],palign->sites2[index2-1],
		 palign->sites2[index2-1],X[palign->sites2[index2-1]]);
	assert(palign->sites1[index2-1] >= 0);
	assert(palign->sites1[index2-1] < N);
	assert(index2 <= 0 || (palign->orientation ? X[M+1] - X[M+1-palign->sites2[index2-1]] : X[palign->sites2[index2-1]]) <= centerX);
      }

      /* map centers on to other map */
      double centerYonX, centerXonY;
    
      if(!palign->orientation){
	if(index1 >= K){/* translate to right end of X */
	  centerYonX = X[palign->sites2[K-1]] + (centerY - Y[palign->sites1[K-1]]);
	} else if(index1 <= 0){/* translate to left end of X */
	  centerYonX = X[palign->sites2[0]] + (centerY - Y[palign->sites1[0]]);
	} else {/* interpolate between aligned sites */
	  double scale = (centerY - Y[palign->sites1[index1-1]])/(Y[palign->sites1[index1]] - Y[palign->sites1[index1-1]]);
	  centerYonX = X[palign->sites2[index1-1]] + scale * (X[palign->sites2[index1]] - X[palign->sites2[index1-1]]);
	}

	if(index2 >= K){/* translate to right end of Y */
	  centerXonY = Y[palign->sites1[K-1]] + (centerX - X[palign->sites2[K-1]]);
	} else if(index2 <= 0){/* translate to left end of Y */
	  centerXonY = Y[palign->sites1[0]] + (centerX - X[palign->sites2[0]]);
	} else {/* interpolate between aligned sites */
	  double scale = (centerX - X[palign->sites2[index2-1]])/(X[palign->sites2[index2]] - X[palign->sites2[index2-1]]);
	  centerXonY = Y[palign->sites1[index2-1]] + scale * (Y[palign->sites1[index2]] - Y[palign->sites1[index2-1]]);
	}
      } else {
	if(index1 >= K){/* translate to right end of X */
	  centerYonX = (X[M+1] - X[M+1-palign->sites2[K-1]]) + (centerY - Y[palign->sites1[K-1]]);
	} else if(index1 <= 0){/* translate to left end of X */
	  centerYonX = (X[M+1] - X[M+1-palign->sites2[0]]) + (centerY - Y[palign->sites1[0]]);
	} else {/* interpolate between aligned sites */
	  double scale = (centerY - Y[palign->sites1[index1-1]])/(Y[palign->sites1[index1]] - Y[palign->sites1[index1-1]]);
	  centerYonX = (X[M+1] - X[M+1-palign->sites2[index1-1]]) + 
	    scale * (X[M+1-palign->sites2[index1-1]] - X[M+1-palign->sites2[index1]]);
	}

	if(index2 >= K){/* translate to right end of Y */
	  centerXonY = Y[palign->sites1[K-1]] + (centerX - (X[M+1]-X[M+1-palign->sites2[K-1]]));
	} else if(index2 <= 0){/* translate to left end of Y */
	  centerXonY = Y[palign->sites1[0]] + (centerX - (X[M+1]-X[M+1-palign->sites2[0]]));
	} else {/* interpolate between aligned sites */
	  double scale = (centerX - (X[M+1]-X[M+1-palign->sites2[index2-1]]))/
	    (X[M+1-palign->sites2[index2-1]] - X[M+1-palign->sites2[index2]]);
	  centerXonY = Y[palign->sites1[index2-1]] + scale * (Y[palign->sites1[index2]] - Y[palign->sites1[index2-1]]);
	}
      }
       
      distanceYtoX = ((centerX - centerYonX)+(centerXonY - centerY))*0.5;
      /*    if(fabs(distanceYtoX) > centerX+centerY ||
	    (palign->mapid1==303 && palign->mapid2==2277) || 
	    (palign->mapid1==2277 && palign->mapid2==303)){
	    printf("palign:mapid1(Y)=%d,mapid2(X)=%d,distanceYtoX=%0.3f,centerX=%0.3f,centerY=%0.3f,centerYonX=%0.3f,centerXonY=%0.3f\n",
	    palign->mapid1,palign->mapid2,distanceYtoX,centerX,centerY,centerYonX,centerXonY);
	    printf("orientation=%d,Y[N+1]=%0.3f,X[M+1]=%0.3f\n",palign->orientation,Y[N+1],X[M+1]);
	    fflush(stdout);
	    }*/
      if(DEBUG) assert(fabs(distanceYtoX) <= centerX+centerY+0.00001);
    }

    if(distanceYtoX == 0.0)
      distanceYtoX = 0.0001;
      
    if(distanceYtoX >= 0.0){/* X is right of Y */
      pedge->Lmap = palign->mapid1;
      pedge->Rmap = palign->mapid2;
      pedge->Lflipped = 0;
      pedge->Rflipped = palign->orientation;
      pedge->distance = distanceYtoX;
    } else {/* Y is right of X */
      pedge->Lmap = palign->mapid2;
      pedge->Rmap = palign->mapid1;
      pedge->Lflipped = palign->orientation;
      pedge->Rflipped = 0;
      pedge->distance = -distanceYtoX;
    }
    if(DEBUG) assert(pedge->distance >= 0.0);
    if(DEBUG>=2) assert(pedge->Lmap != pedge->Rmap);
  }
  
  if(VERB && GraphWeight > 0.0 && CovCnt > 0){
    printf("Average Graph Edge Weight = %0.3f, Max = %0.3f\n", CovSum/CovCnt, MaxCov);
    fflush(stdout);
  }

  if(VERB/* HERE HERE */){
    printf("Initialized %u/%u edges:time=%0.6f\n",Numedges,Maxedges,wtime());
    fflush(stdout);
  }

  /* link the edges to the nodes */
  /* first count how many edges there are at each node */
  for(int i = 0; i < numnodes; i++)
    nodeA[i].edgecnt = 0;

  pedge = &edgeA[1];
  for(unsigned int i = 1; i <= Numedges; i++, pedge++){
    /* initially Gmap[] index is the same as the node[] index */
    if(DEBUG>=2){
      assert(nodeA[pedge->Lmap].nodeid == pedge->Lmap);
      assert(nodeA[pedge->Rmap].nodeid == pedge->Rmap);
    }

    nodeA[pedge->Lmap].edgecnt++;
    nodeA[pedge->Rmap].edgecnt++;

    if(VERB>=3 && (pedge->Lmap == 0 || pedge->Rmap == 0)){
      printf("edgeA[i=%d]:Lmap=%d,Rmap=%d:nodeA[Lmap].edgecnt -> %d, nodeA[Rmap].edgecnt -> %d\n",i,pedge->Lmap,pedge->Rmap,nodeA[pedge->Lmap].edgecnt,nodeA[pedge->Rmap].edgecnt);
      fflush(stdout);
    }
    if(DEBUG>=2) assert(pedge->Lmap != pedge->Rmap);
  }
    
  /* now allocate the right number of edgeid's at each node */
  pnode = nodeA;
  for(int i = 0; i < numnodes; i++, pnode++){
    pnode->maxedgealloc(pnode->edgecnt);
    pnode->edgecnt = 0;
  }
  
  /* Now add the edge index values to the nodes */
  for(unsigned int i = 1; i <= Numedges; i++){
    int Lmap = edgeA[i].Lmap;
    pnode = &nodeA[Lmap];
    if(DEBUG>=2) assert(pnode->edgecnt < pnode->maxedge);
    pnode->edgeid[pnode->edgecnt++] = i;

    int Rmap = edgeA[i].Rmap;
    pnode = &nodeA[Rmap];
    if(DEBUG>=2) assert(pnode->edgecnt < pnode->maxedge);
    pnode->edgeid[pnode->edgecnt++] = i;
  }  

  if(!topnodes){ /* initialize topnodes[] to all nodes */
    numtopnodes = numnodes;
    topnodes = new Cnode *[maxnodes];
    for(int i=0;i < numnodes; i++)
      topnodes[i] = (Cnode *) &nodeA[i];
  }

  /* find largest edgecnt */
  int maxedgecnt = 0;
  int maxedgecnt_node = -1;
  for(int i = 0; i < numnodes; i++)
    if(nodeA[i].edgecnt > maxedgecnt){
      maxedgecnt = nodeA[i].edgecnt;
      maxedgecnt_node = i;
    }

  if(VERB){
    printf("Initialized graph with %d nodes and %d edges (max edgecnt=%d):time=%0.6f(total=%0.6f)\n",numnodes, Numedges, maxedgecnt,wtime(),mtime());
    printf("  Node size= %lu bytes, Edge size= %lu bytes\n", sizeof(CnodeA), sizeof(CedgeA) + sizeof(int));
    if(VERB>=2)
      dumpmemmap();
    if(VERB >=2){
      int node1 = maxedgecnt_node;
      CnodeA *Vnode = &nodeA[node1];
      int map1 = Vnode->mapid;
      double startloc1 = Gmap[map1]->startloc;
      double reflen = startloc1;
      if(startloc1 > -1000.0){
	reflen = Gmap[map1]->endloc - startloc1;
	startloc1 += 0.5*Vnode->len;
      }
      printf("N%d:",nodeorder[map1]);
      printf("nodeid=%d,",node1);
      printf("mapid=%d", map1);
      printf("(id=%lld),",Gmap[map1]->id);
      printf("cov=%0.2f,",Vnode->cov);
      printf("PosEdge=%d,",Vnode->PosEdge);
      printf("edgecnt=%d,",Vnode->edgecnt);
      printf("loc=%0.3f",startloc1);
      printf("(%0.3f..",Gmap[map1]->startloc);
      printf("%0.3f,",Gmap[map1]->endloc);
      printf("%3f..",Gmap[map1]->startmatch);
      printf("%3f),",Gmap[map1]->endmatch);
      printf("flip=%d,",Gmap[map1]->flipped);
      printf("len=%0.3f",Vnode->len);
      printf("(%0.3f),",reflen);
      printf("left=%0.1f,right=%0.1f,",Vnode->backward,Vnode->forward);
      printf("trim=%d,%d,",Vnode->trimL,Vnode->trimR);
      printf("N=%d,",Gmap[map1]->numsite[0]);
      printf("Qscore=%0.3f\n",Gmap[map1]->Qscore);

      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	CnodeA *Wnode = &nodeA[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
	int map2 = nodeA[(VWedge->Lmap == Vnode->mapid) ? VWedge->Rmap : VWedge->Lmap].mapid;
	double startloc2 = Gmap[map2]->startloc;
	if(startloc2 > -1000.0)
	  startloc2 += 0.5*Wnode->len;
	printf(" %cEdge%d = (%d,%d,%0.3f,%d,%d) to N%d(id=%lld):valid=%d,cov=%0.2f,wt=%0.3f,loc=%0.3f(delta=%0.3f)",
	       VWedge->alignid ? ' ' : 'S', vw, VWedge->Lmap, VWedge->Rmap, VWedge->distance, VWedge->Lflipped,VWedge->Rflipped,
	       nodeorder[map2], Gmap[map2]->id, VWedge->valid, VWedge->cov, 0.0, startloc2, startloc2-startloc1);
	if(VWedge->alignid){
	  Calign *p = alignment[VWedge->alignid - 1];
	  int LI = p->sites1[0];
	  int LJ = p->sites2[0];
	  int RI = p->sites1[p->numpairs - 1];
	  int RJ = p->sites2[p->numpairs - 1];
	  int N = Gmap[p->mapid1]->numsite[0];
	  int M = Gmap[p->mapid2]->numsite[0];
	  int mis = RI-LI+1+RJ-LJ+1-2*p->numpairs;
	  double len1 = Gmap[p->mapid1]->site[0][RI] - Gmap[p->mapid1]->site[0][LI];
	  double len2 = p->orientation ? 
	    Gmap[p->mapid2]->site[0][M+1-LJ] - Gmap[p->mapid2]->site[0][M+1-RJ] :
	    Gmap[p->mapid2]->site[0][RJ] - Gmap[p->mapid2]->site[0][LJ];
	  printf(",Yid=%d,Xid=%d:LI=%d,RI=%d,N=%d:LJ=%d,RJ=%d,M=%d,mis=%d,R=%0.3f(/100kb),score=%0.3f,logPV=%0.2f\n",
		 p->mapid1,p->mapid2,LI,RI,N,LJ,RJ,M,mis,100.0*mis/(len1+len2),p->score,p->logPV);
	} else
	  printf("\n");
      }
    }
    fflush(stdout);
  }

  /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
  if(DEBUG>=2){
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edgeA[i].valid < 0)
	continue;
      assert(edgeA[i].mark == 0);
    }
    size_t edgecnt = 0;
    for(int i = 0; i < numnodes; i++){
      CnodeA *Vnode = &nodeA[i];
      if(Vnode->mark >= 2){
	assert(Vnode->edgecnt <= 0);
	continue;
      }
      edgecnt += Vnode->edgecnt;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	assert(1 <= VWedgeid && VWedgeid <= Numedges);
	assert(edgeA[VWedgeid].mark == 0);
      }
    }    
    if(DEBUG>=3){
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edgeA[Vnode->edgeid[vw]].mark++;
      }    
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	  if(VWedge->mark != 2){
	    printf("i=%d/%d: mapid=%d, id=%lld : vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d:align->logPV= %0.2f,or=%d\n",
		   i,numnodes,Vnode->mapid,gmap[Vnode->mapid]->id, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark,
		   VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped,alignment[VWedge->alignid - 1]->logPV, alignment[VWedge->alignid - 1]->orientation);
	    fflush(stdout);
	    assert(VWedge->mark == 2);
	  }
	}
      }    
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edgeA[Vnode->edgeid[vw]].mark = 0;
      }    
      /* check that topnodes[0..numtopnodes-1] reach all nodeA[] that do not have mark >= 2 */
      for(int i = 0; i < numtopnodes; i++){
        CnodeA *Vnode = (CnodeA *) topnodes[i];
	assert(Vnode->mark == 0);
	Vnode->mark = 2;
      }
      for(int i = 0; i < numnodes; i++)
	assert(nodeA[i].mark >= 2);
      for(int i = 0; i < numtopnodes; i++){
        CnodeA *Vnode = (CnodeA *) topnodes[i];
	Vnode->mark = 0;
      }
    }
    if((edgecnt % 2) != 0){
      printf("sum of nodeA[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",edgecnt,numnodes,Numedges);
      fflush(stdout);
    }
  }

  /* Sort edges by nodeid */
  for(int i = 0; i < numnodes; i++){
    CnodeA *Vnode = &nodeA[i];
    /* sort edges in ascending order of nodeid */
    qsort_nodeid = Vnode->nodeid;/* global arg to edgenodeinc() */

    if(Vnode->edgeid)
      qsort(Vnode->edgeid,Vnode->edgecnt,sizeof(unsigned int),(intcmp *)edgeidnodeinc);
    else     
      if(DEBUG) assert(Vnode->edgecnt == 0);
  }

  if(VERB){
    printf("Sorted edges at nodes by nodeid:time=%0.6f(total=%0.6f)\n",wtime(),mtime());
    fflush(stdout);
  }

  /* initialize node[i].Enodeid[0..node[i].maxedge-1] */
  for(int i = 0; i < numnodes; i++){
    CnodeA *Vnode = &nodeA[i];
    for(int vw = 0; vw < Vnode->edgecnt; vw++){
      CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
      Vnode->Enodeid[vw] = (VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap;
      if(DEBUG>=2 && vw > 0) assert(Vnode->Enodeid[vw] >= Vnode->Enodeid[vw-1]);
    }
  }

  /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
  if(DEBUG>=2){
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edgeA[i].valid < 0)
	continue;
      assert(edgeA[i].mark == 0);
    }
    size_t edgecnt = 0;
    for(int i = 0; i < numnodes; i++){
      CnodeA *Vnode = &nodeA[i];
      if(Vnode->mark >= 2){
	assert(Vnode->edgecnt <= 0);
	continue;
      }
      edgecnt += Vnode->edgecnt;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	assert(1 <= VWedgeid && VWedgeid <= Numedges);
	assert(edgeA[VWedgeid].mark == 0);
      }
    }    
    if(DEBUG>=3){
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edgeA[Vnode->edgeid[vw]].mark++;
      }    
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	  if(VWedge->mark != 2){
	    printf("i=%d/%d: mapid=%d, id=%lld : vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d:align->logPV= %0.2f,or=%d\n",
		   i,numnodes,Vnode->mapid,gmap[Vnode->mapid]->id, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark,
		   VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped,alignment[VWedge->alignid - 1]->logPV, alignment[VWedge->alignid - 1]->orientation);
	    fflush(stdout);
	    assert(VWedge->mark == 2);
	  }
	}
      }    
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edgeA[Vnode->edgeid[vw]].mark = 0;
      }    
      /* check that topnodes[0..numtopnodes-1] reach all nodeA[] that do not have mark >= 2 */
      for(int i = 0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	assert(Vnode->mark == 0);
	Vnode->mark = 2;
      }
      for(int i = 0; i < numnodes; i++)
	assert(nodeA[i].mark >= 2);
      for(int i = 0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	Vnode->mark = 0;
      }
    }
    if((edgecnt % 2) != 0){
      printf("sum of nodeA[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",edgecnt,numnodes,Numedges);
      fflush(stdout);
    }
  }

  unsigned int edgecnt = Numedges;

  if(!SkipFilter){  /* check for duplicate alignments */
    int dupcnt = 0;
    for(int i = 0; i < numtopnodes;i++){
      CnodeA *Vnode = (CnodeA *) topnodes[i];
      for(int vw = 1; vw < Vnode->edgecnt; vw++){
	CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	if(DEBUG>=2) assert(!(VWedge->valid < 0));
	//      if(VWedge->valid < 0)
	//	continue;
	if(DEBUG>=2) assert(((VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap) == Vnode->Enodeid[vw]);
	if(Vnode->Enodeid[vw] == Vnode->Enodeid[vw-1]){
	  if(VERB && !dupcnt){
	    printf("WARNING:Duplicate edge found from node=%d(mapid=%d,id=%lld):\n", Vnode->nodeid,Vnode->mapid,Gmap[Vnode->mapid]->id);
	    fflush(stdout);
	    CedgeA *edge1 = &edgeA[Vnode->edgeid[vw-1]];
	    Calign *align1 = alignment[edge1->alignid - 1];
	    printf("\t 1. To node %d(mapid=%d,id=%lld):edgeid=%d,distance=%0.3f,align:score=%0.6f,Pvalue=%0.3f,numpairs=%d,filename=%s\n",
		   Vnode->Enodeid[vw-1],node[Vnode->Enodeid[vw-1]].mapid,Gmap[node[Vnode->Enodeid[vw-1]].mapid]->id,edge1->edgeid,edge1->distance,
		   align1->score,align1->logPV,align1->numpairs,align_filename[max(0,align1->fileid)]);
	    fflush(stdout);
	    CedgeA *edge2 = &edgeA[Vnode->edgeid[vw]];
	    Calign *align2 = alignment[edge2->alignid - 1];
	    printf("\t 2. To node %d(mapid=%d,id=%lld):edgeid=%d,distance=%0.3f,align:score=%0.6f,Pvalue=%0.3f,numpairs=%d,filename=%s\n",
		   Vnode->Enodeid[vw],node[Vnode->Enodeid[vw]].mapid,Gmap[node[Vnode->Enodeid[vw]].mapid]->id,edge2->edgeid,edge2->distance,
		   align2->score,align2->logPV,align2->numpairs,align_filename[max(0,align2->fileid)]);
	    fflush(stdout);
	  }
	  dupcnt++;
	  edge_delete(VWedge,Vnode);
	}
      }
    }
    if(dupcnt > 0){
      GCedges(edgecnt,edgeA);
      GCnodesA();
    }
    if(VERB){
      if(dupcnt > 0)
	printf("Warning:Deleted %d duplicate alignments(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",dupcnt,edgecnt,numtopnodes,wtime(),mtime());
      else
	printf("Checked for  duplicate alignments(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",edgecnt,numtopnodes,wtime(),mtime());
      fflush(stdout);
    }

    /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
    if(DEBUG>=2){
      for(unsigned int i = 1; i <= Numedges; i++){
	if(edgeA[i].valid < 0)
	  continue;
	assert(edgeA[i].mark == 0);
      }
      size_t edgecnt = 0;
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	if(Vnode->mark >= 2){
	  assert(Vnode->edgecnt <= 0);
	  continue;
	}
	edgecnt += Vnode->edgecnt;
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  unsigned int VWedgeid = Vnode->edgeid[vw];
	  assert(1 <= VWedgeid && VWedgeid <= Numedges);
	  assert(edgeA[VWedgeid].mark == 0);
	}
      }    
      if(DEBUG>=3){
	for(int i = 0; i < numnodes; i++){
	  CnodeA *Vnode = &nodeA[i];
	  for(int vw = 0; vw < Vnode->edgecnt; vw++)
	    edgeA[Vnode->edgeid[vw]].mark++;
	}    
	for(int i = 0; i < numnodes; i++){
	  CnodeA *Vnode = &nodeA[i];
	  for(int vw = 0; vw < Vnode->edgecnt; vw++){
	    CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	    if(VWedge->mark != 2){
	      printf("i=%d/%d: mapid=%d, id=%lld : vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d:align->logPV= %0.2f,or=%d\n",
		     i,numnodes,Vnode->mapid,gmap[Vnode->mapid]->id, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark,
		     VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped,alignment[VWedge->alignid - 1]->logPV, alignment[VWedge->alignid - 1]->orientation);
	      fflush(stdout);
	      assert(VWedge->mark == 2);
	    }
	  }
	}    
	for(int i = 0; i < numnodes; i++){
	  CnodeA *Vnode = &nodeA[i];
	  for(int vw = 0; vw < Vnode->edgecnt; vw++)
	    edgeA[Vnode->edgeid[vw]].mark = 0;
	}    
	/* check that topnodes[0..numtopnodes-1] reach all nodeA[] that do not have mark >= 2 */
	for(int i = 0; i < numtopnodes; i++){
	  CnodeA *Vnode = (CnodeA *) topnodes[i];
	  assert(Vnode->mark == 0);
	  Vnode->mark = 2;
	}
	for(int i = 0; i < numnodes; i++)
	  assert(nodeA[i].mark >= 2);
	for(int i = 0; i < numtopnodes; i++){
	  CnodeA *Vnode = (CnodeA *) topnodes[i];
	  Vnode->mark = 0;
	}
      }
      if((edgecnt % 2) != 0){
	printf("sum of nodeA[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",edgecnt,numnodes,Numedges);
	fflush(stdout);
      }
    }

    int smallcnt = 0;
    /* filter out maps smaller than MinLen */
    for(int i=0; i < numtopnodes; i++){
      CnodeA *Vnode = (CnodeA *) topnodes[i];
      if(Vnode->len < MinLen){
	Vnode->mark = 2;
	smallcnt++;
	/* mark all edges out of this node for deletion */
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edge_delete(&edgeA[Vnode->edgeid[vw]],Vnode);
      }
    }
    
    if(smallcnt > 0){
      GCedges(edgecnt,edgeA);
      GCnodesA();
    }
    if(VERB && (VERB/*>=2*/ || smallcnt > 0)){
      printf("Deleted %d Maps smaller than %0.3f kb(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",smallcnt,MinLen,edgecnt,numtopnodes,wtime(),mtime());
      fflush(stdout);
    }

    if(VERB>=2){
      GCedges(edgecnt,edgeA);
      GCnodesA();
      printf("After calling GCedges,GCnodes: edges=%d,nodes=%d\n",edgecnt,numtopnodes);
      fflush(stdout);
    }

    /* check if any nodes have excessively high edgecnt (see MaxRelCoverage) : if so keep only the best edges */
    int highedgecnt = floor(0.5 + MaxRelCoverage * (edgecnt*2.0)/numnodes);

    if(MaxAbsCoverage > 0 && MaxAbsCoverage < highedgecnt)
      highedgecnt = MaxAbsCoverage;
  
    int highedgecnt2 = (MaxAbsCoverage2 > 0) ? min(highedgecnt,MaxAbsCoverage2) : highedgecnt;

    int highedgedelcnt = 0;
    if(maxedgecnt > highedgecnt){
      /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
      if(DEBUG>=2){
	for(int i = 0; i < numnodes; i++){      /* mark all edges with mark = 0 */
	  CnodeA *Vnode = &nodeA[i];
	  for(int vw = 0; vw < Vnode->edgecnt; vw++)
	    assert(edgeA[Vnode->edgeid[vw]].mark == 0);
	}    
	if(DEBUG>=3){
	  for(int i = 0; i < numnodes; i++){
	    CnodeA *Vnode = &nodeA[i];
	    for(int vw = 0; vw < Vnode->edgecnt; vw++)
	      edgeA[Vnode->edgeid[vw]].mark++;
	  }    
	  for(int i = 0; i < numnodes; i++){
	    CnodeA *Vnode = &nodeA[i];
	    for(int vw = 0; vw < Vnode->edgecnt; vw++){
	      CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	      if(VWedge->mark != 2){
		Calign *p = alignment[VWedge->alignid - 1];
		printf("i=%d/%d: mapid=%d, id=%lld : vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d:align->logPV= %0.2f,or=%d\n",
		       i,numnodes,Vnode->mapid,gmap[Vnode->mapid]->id, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark,
		       VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped,p->logPV, p->orientation);
		fflush(stdout);
		assert(VWedge->mark == 2);
	      }
	    }
	  }    
	  for(int i = 0; i < numnodes; i++){
	    CnodeA *Vnode = &nodeA[i];
	    for(int vw = 0; vw < Vnode->edgecnt; vw++)
	      edgeA[Vnode->edgeid[vw]].mark = 0;
	  }    
	}
      }

      /* Sort edges by logPV */
      //    #pragma omp parallel for schedule(dynamic,64)
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	if(Vnode->edgecnt > highedgecnt2){
	  /* sort edges in descending order of logPV */
	  qsort(Vnode->edgeid,Vnode->edgecnt,sizeof(unsigned int),(intcmp *)edgePVdec);
	  if(DEBUG>=2){
	    for(int vw = 1; vw < Vnode->edgecnt; vw++){
	      CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	      CedgeA *VWm1edge = &edgeA[Vnode->edgeid[vw-1]];
	      if(DEBUG && !(alignment[VWedge->alignid - 1]->logPV <= alignment[VWm1edge->alignid - 1]->logPV)){
		printf("node[%d](id=%lld):vw=%d:edgeA[vw]->logPV=%0.2f,edgeA[vw-1]->logPV=%0.2f\n",
		       i,gmap[Vnode->mapid]->id,vw,alignment[VWedge->alignid - 1]->logPV,alignment[VWm1edge->alignid - 1]->logPV);
		fflush(stdout);
		assert(alignment[VWedge->alignid - 1]->logPV <= alignment[VWm1edge->alignid - 1]->logPV);
	      }
	    }
	  }
	}
      }

      if(VERB){
	printf("Sorted edges at nodes with %d edges by logPV:time=%0.6f(total=%0.6f)\n",highedgecnt2,wtime(),mtime());
	fflush(stdout);
      }

      /* increment the edge->mark for each end that is superfluous */
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	if(Vnode->edgecnt > highedgecnt2){ /* mark all edges above the best highedgecnt2 once and above highedgecnt twice for deletion */
	  if(DEBUG>=2){
	    for(int vw = 1; vw < Vnode->edgecnt; vw++){
	      CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	      CedgeA *VWm1edge = &edgeA[Vnode->edgeid[vw-1]];
	      if(DEBUG && !(alignment[VWedge->alignid - 1]->logPV <= alignment[VWm1edge->alignid - 1]->logPV)){
		printf("nodeA[%d](id=%lld):vw=%d:edgeA[vw]->logPV=%0.2f,edgeA[vw-1]->logPV=%0.2f\n",
		       i,gmap[Vnode->mapid]->id,vw,alignment[VWedge->alignid - 1]->logPV,alignment[VWm1edge->alignid - 1]->logPV);
		fflush(stdout);
		assert(alignment[VWedge->alignid - 1]->logPV <= alignment[VWm1edge->alignid - 1]->logPV);
	      }
	    }
	  }

	  int vw = highedgecnt2;
	  int vwmax = min(highedgecnt,Vnode->edgecnt);
	  for(; vw < vwmax; vw++)
	    edgeA[Vnode->edgeid[vw]].mark++;
	  for(vwmax = Vnode->edgecnt; vw < vwmax; vw++)
	    edgeA[Vnode->edgeid[vw]].mark += 2;
	}
      }

      if(VERB){
	if(highedgecnt2 < highedgecnt)
	  printf("Marked edges with rank above %d on both nodes and above %d on one:time=%0.6f(total=%0.6f)\n",highedgecnt2,highedgecnt,wtime(),mtime());
	else
	  printf("Marked edges with rank above %d on both nodes:time=%0.6f(total=%0.6f)\n",highedgecnt,wtime(),mtime());
	fflush(stdout);
      }

      if(VERB>=2){/* count how many topnodes have no edges, no edges with mark < 3, and how many have node mark >= 2 */
	int nodecnt0 = 0, nodecnt1 = 0, nodecnt2 = 0;
	for(int i = 0; i < numtopnodes; i++){
	  CnodeA *Vnode = (CnodeA *)topnodes[i];
	  if(Vnode->edgecnt <= 0){
	    nodecnt0++;
	    continue;
	  }

	  int edgecnt = 0;// how many edges have mark < 3
	  for(int vw = 0; vw < Vnode->edgecnt; vw++){
	    if(DEBUG) assert(edgeA[Vnode->edgeid[vw]].mark <= 4);

	    if(edgeA[Vnode->edgeid[vw]].mark < 3)
	      edgecnt++;
	  }

	  if(!edgecnt){
	    if(!nodecnt1){
	      printf("Node[%d](mapid=%d, id=%lld): edgecnt=%d : none with mark < 3:\n", i, Vnode->mapid, gmap[Vnode->mapid]->id,Vnode->edgecnt);
	      for(int vw = 0; vw < Vnode->edgecnt; vw++){
		CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
		Calign *p = alignment[VWedge->alignid - 1];
		printf("\t edgeA[%d] : edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,distance=%0.3f,Lflipped=%d,Rflipped=%d:align->logPV= %0.2f\n",
		       vw, VWedge->edgeid, VWedge->mark, VWedge->Lmap, VWedge->Rmap, VWedge->distance, VWedge->Lflipped, VWedge->Rflipped, p->logPV);
	      }
	      fflush(stdout);
	    }
	    nodecnt1++;
	    continue;
	  }
	  if(Vnode->mark >= 2)
	    nodecnt2++;
	}
	
	printf("topnodes= %d: Includes %d nodes with no edges, %d nodes with all edges marked for deletion + %d other nodes with mark >= 2\n", numtopnodes, nodecnt0, nodecnt1, nodecnt2);
	fflush(stdout);
      }

      /* now delete edges marked at least 3 times for deletion (ie below highedgecnt2 for both edges and below highedgecnt for at least one edge ) */
      int NodeIcnt = 0;/* count of nodes with at least highedgecnt alignments */
      long long EdgeIcnt = 0;
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  if(edgeA[Vnode->edgeid[vw]].mark >= 3){
	    highedgedelcnt++;
	    edge_delete(&edgeA[Vnode->edgeid[vw]],Vnode);
	  }
	  edgeA[Vnode->edgeid[vw]].mark = 0;
	}
	if(Vnode->edgecnt >= highedgecnt){
	  NodeIcnt++;
	  EdgeIcnt += Vnode->edgecnt;
	}
      }
      GCedges(edgecnt,edgeA);
      GCnodesA();
      if(VERB){
	if(highedgecnt2 < highedgecnt)
	  printf("Deleted %d edges with both nodes that had excessive edgecnt above %d and one above %d (involving %d nodes with an average of %0.1f edges each) (edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
		 highedgedelcnt,highedgecnt2,highedgecnt,NodeIcnt,((double)EdgeIcnt)/max(1,NodeIcnt),edgecnt,numtopnodes,wtime(),mtime());
	else
	  printf("Deleted %d edges with both nodes that had excessive edgecnt above %d (involving %d nodes with an average of %0.1f edges each) (edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
		 highedgedelcnt,highedgecnt,NodeIcnt,((double)EdgeIcnt)/max(1,NodeIcnt), edgecnt,numtopnodes,wtime(),mtime());
	fflush(stdout);
      }

      if(MaxRelCoverage <= 2){
	int new_highedgecnt;
	while((new_highedgecnt = floor(0.5 + MaxRelCoverage * (edgecnt*2.0)/numnodes)) < highedgecnt){
	  highedgecnt = new_highedgecnt;
	  highedgedelcnt = 0;
	  for(int i = 0; i < numnodes; i++){
	    CnodeA *Vnode = &nodeA[i];
	    if(Vnode->edgecnt > highedgecnt){
	      for(int vw = highedgecnt; vw < Vnode->edgecnt; vw++)
		edge_delete(&edgeA[Vnode->edgeid[vw]],Vnode);
	      highedgedelcnt += Vnode->edgecnt - highedgecnt;
	    }
	  }
	  GCedges(edgecnt,edgeA);
	  GCnodesA();
	  if(VERB){
	    printf("Deleted %d edge ends at nodes with excessive edgecnt above %d(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",highedgedelcnt,highedgecnt,edgecnt,numtopnodes,wtime(),mtime());
	    fflush(stdout);
	  }
	}
      }
    } else if(VERB){
      printf("MaxRelCoverage=%d:highedgecnt=%d,maxedgecnt=%d\n",MaxRelCoverage,highedgecnt,maxedgecnt);
      fflush(stdout);
    }
  }

  if(maxedges < 0){
    maxedges = 0;
    for(int i = 0; i < numtopnodes; i++) {
      CnodeA *Vnode = (CnodeA *) topnodes[i];
      int edgecnt = Vnode->edgecnt;
      if(edgecnt > maxedges)
	maxedges = edgecnt;
    }
    if(VERB>=2){
      printf("max alignments for any map =%d\n",maxedges);
      fflush(stdout);
    }
    if(maxedges >= 32*1024){
      printf("maxedges=%d exceeds largest 16 bit value (use -MaxRelCoverage to limit value)\n",maxedges);
      exit(1);
    }
  }

  if(!PosEdgeInit)
    PosEdgeA();

  if(VERB>=2 && (refmap_filename || Gmap[0]->startloc>0.0)){/* display node locations */
    printf("Known map left end locations (mapid = index in input vfixx file(s)):\n");

    for(int i=0; i < numnodes;i++){
      int node1;
      if(i < nummaps)
	node1 = nodemap[i];
      else {
	if(DEBUG) assert(CROSSBREAK);
	node1 = i;
      }
      CnodeA *Vnode = &nodeA[node1];
      int map1 = Vnode->mapid;
      if(DEBUG>=2) assert(nodemap[nodeorder[map1]] == map1);
      if(DEBUG>=2) assert(Vnode->mapid == map1);
      if(DEBUG>=2 && i < nummaps){
	assert(node1 == map1);
	assert(nodeorder[map1] == i);
      } 
      double startloc1 = Gmap[map1]->startloc;      
      if(VERB && startloc1 > 0.0)
      	printf("N%d: mapid=%d(id=%lld) : loc=%0.3f\n", nodeorder[map1], map1, Gmap[map1]->id, startloc1);
    }    
    fflush(stdout);
  }

  if(VERB){
    printf("graph build completed:numedges=%d/%d,numnodes=%d/%d,time=%0.6f(total=%0.6f)\n",edgecnt,Numedges,numtopnodes,numnodes,wtime(),mtime());
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
  if(DEBUG>=2){
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edgeA[i].valid < 0)
	continue;
      assert(edgeA[i].mark == 0);
    }
    size_t edgecnt = 0;
    for(int i = 0; i < numnodes; i++){
      CnodeA *Vnode = &nodeA[i];
      if(Vnode->mark >= 2){
	assert(Vnode->edgecnt <= 0);
	continue;
      }
      edgecnt += Vnode->edgecnt;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	assert(1 <= VWedgeid && VWedgeid <= Numedges);
	assert(edgeA[VWedgeid].mark == 0);
      }
    }    
    if(DEBUG>=3){
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edgeA[Vnode->edgeid[vw]].mark++;
      }    
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	  if(VWedge->mark != 2){
	    printf("i=%d/%d: mapid=%d, id=%lld : vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d:align->logPV= %0.2f,or=%d\n",
		   i,numnodes,Vnode->mapid,gmap[Vnode->mapid]->id, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark,
		   VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped,alignment[VWedge->alignid - 1]->logPV, alignment[VWedge->alignid - 1]->orientation);
	    fflush(stdout);
	    assert(VWedge->mark == 2);
	  }
	}
      }    
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edgeA[Vnode->edgeid[vw]].mark = 0;
      }    
      /* check that topnodes[0..numtopnodes-1] reach all nodeA[] that do not have mark >= 2 */
      for(int i = 0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	assert(Vnode->mark == 0);
	Vnode->mark = 2;
      }
      for(int i = 0; i < numnodes; i++)
	assert(nodeA[i].mark >= 2);
      for(int i = 0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	Vnode->mark = 0;
      }
    }
    if((edgecnt % 2) != 0){
      printf("sum of nodeA[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",edgecnt,numnodes,Numedges);
      fflush(stdout);
    }
  }
}

/** display all remaining graph nodes and edges with cov >= mincov. If suffix != 0, output is to file, otherwise to stdout */

extern void graph_displayA(int mincov, int verb, char *suffix);

void graph_display(int mincov, int verb, char *suffix)
{
  if(node == NULL && nodeA != NULL)
    return graph_displayA(mincov, verb, suffix);

  FILE *fp = stdout;
  if(suffix){
    char filename[1024];
    sprintf(filename,"%s_%s.graph",output_prefix,suffix);
    if((fp = fopen(filename,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("Unable to write to %s:errno=%d:%s\n",filename,eno,err);
      fflush(stdout);
    }
    if(VERB){
      printf("Generating %s (graph dump)\n",filename);
      fflush(stdout);
    }
  }

  int nodecnt=0;/* count of active nodes */
  int edgecnt=0;/* count of active edges */
  int highcnt=0;/* count of how many times edges with cov >= MinCoverage occured */
  int errcnt=0;/* count how many times edges with cov >= MinCoverage had true offset >= 20kb */
  double errsiz=0.0;
  int maxEcov = 0;
  int maxNcov = 0;

  if(!PosEdgeInit)
    PosEdge();

  for(int i=0; i < numnodes;i++){
    Cnode *Vnode;
    int map1,node1;
    if(i < nummaps)
      node1 = nodemap[i];
    else {
      if(DEBUG) assert(CROSSBREAK);
      node1 = i;
    }
    Vnode = &node[node1];
    map1 = Vnode->mapid;
    if(DEBUG>=2) assert(nodemap[nodeorder[map1]] == map1);
    if(DEBUG>=2) assert(Vnode->mapid == map1);
    if(DEBUG>=2 && i < nummaps){
      assert(node1 == map1);
      assert(nodeorder[map1] == i);
    } 

    if(Vnode->mark >= 2)/* no longer in topnodes[] */
      continue;
    nodecnt++;
    if(Vnode->cov > maxNcov)
      maxNcov = Vnode->cov;

    double startloc1 = Gmap[map1]->startloc;
    double reflen = startloc1;
    if(startloc1 > -1000.0){
      reflen = Gmap[map1]->endloc - startloc1;
      startloc1 += 0.5*Vnode->len;
    }
    if(VERB >= verb){
	//	fprintf(fp,"node%d:nodeid=%d,mapid=%d(id=%d),cov=%0.2f,PosEdge=%d,edgecnt=%d,loc=%0.3f(%0.3f..%0.3f,%0.3f..%0.3f),flip=%d,len=%0.3f(%0.3f),left=%0.1f,right=%0.1f,Ccntmax=%d,%d,trim=%d,%d,N=%d,Qscore=%0.3f",
	//	       nodeorder[map1],node1,map1,Gmap[map1]->id,Vnode->cov,Vnode->PosEdge,Vnode->edgecnt,startloc1,Gmap[map1]->startloc,Gmap[map1]->endloc,Gmap[map1]->startmatch,Gmap[map1]->endmatch,Gmap[map1]->flipped,
	//	       Vnode->len,reflen,Vnode->backward,Vnode->forward,Vnode->Lcntmax,Vnode->Rcntmax,Vnode->trimL,Vnode->trimR,Gmap[map1]->numsite[0],Gmap[map1]->Qscore);
      fprintf(fp,"N%d:",nodeorder[map1]);
      fprintf(fp,"nodeid=%d,",node1);
      fprintf(fp,"(mapid=%d,id=%lld),",map1,Gmap[map1]->id);
      fprintf(fp,"cov=%0.2f,",Vnode->cov);
      fprintf(fp,"PosEdge=%d,",Vnode->PosEdge);
      fprintf(fp,"edgecnt=%d,",Vnode->edgecnt);
      fprintf(fp,"loc=%0.3f",startloc1);
      fprintf(fp,"(%0.3f..",Gmap[map1]->startloc);
      fprintf(fp,"%0.3f,",Gmap[map1]->endloc);
      fprintf(fp,"%3f..",Gmap[map1]->startmatch);
      fprintf(fp,"%3f),",Gmap[map1]->endmatch);
      fprintf(fp,"flip=%d,",Gmap[map1]->flipped);
      fprintf(fp,"len=%0.3f",Vnode->len);
      fprintf(fp,"(%0.3f),",reflen);
      fprintf(fp,"left=%0.1f,right=%0.1f,",Vnode->backward,Vnode->forward);
      fprintf(fp,"Ccntmax=%d,%d,",Vnode->Lcntmax,Vnode->Rcntmax);
      fprintf(fp,"trim=%d,%d,",Vnode->trimL,Vnode->trimR);
      fprintf(fp,"N=%d,",Gmap[map1]->numsite[0]);
      fprintf(fp,"Qscore=%0.3f",Gmap[map1]->Qscore);

      if(Gmap[map1]->name)
	fprintf(fp,":%s\n",Gmap[map1]->name);
      else
	fprintf(fp,"\n");
    }
    for(int vw = 0; vw < Vnode->edgecnt; vw++){
      Cedge *VWedge = &edge[Vnode->edgeid[vw]];
      Cnode *Wnode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
      if(VWedge->valid < 0 || VWedge->cov <= 0.0)
	continue;
      if(VWedge->cov > maxEcov)
	maxEcov = VWedge->cov;

      int map2 = node[(VWedge->Lmap == Vnode->mapid) ? VWedge->Rmap : VWedge->Lmap].mapid;
      if(DEBUG>=2) assert(nodemap[nodeorder[map2]] == map2);
      if(DEBUG>=2 && map1 == map2){
	printf("map1=%d,Vnode=%d(N%d),map2=%d,VWedge=(%d,%d,%0.3f,%d,%d)\n",
	       map1,Vnode->mapid,nodeorder[Vnode->mapid],map2,VWedge->Lmap,VWedge->Rmap,VWedge->distance,VWedge->Lflipped,VWedge->Rflipped);
	fflush(stdout);
	assert(map1 != map2);
      }

      if(map1 < map2)
	edgecnt++;
      if(floor(VWedge->cov+0.5) < mincov)
	continue;
      if(map1 < map2)
	highcnt++;
      double startloc2 = Gmap[map2]->startloc;
      if(startloc2 > -1000.0)
	startloc2 += 0.5*Wnode->len;
      if(VERB >= verb){
	fprintf(fp," %cEdge%d = (N%d,N%d,%0.3f,%d,%d) to N%d(mapid=%d,id=%lld):valid=%d,cov=%0.2f,wt=%0.3f,loc=%0.3f(delta=%0.3f),Lcnt=%d,Rcnt=%d,Bcnt=%d",
	       VWedge->alignid ? ' ' : 'S', vw, nodeorder[node[VWedge->Lmap].mapid], nodeorder[node[VWedge->Rmap].mapid], VWedge->distance, VWedge->Lflipped,VWedge->Rflipped,
	       nodeorder[map2], map2, Gmap[map2]->id,VWedge->valid, VWedge->cov, 0.0, startloc2, startloc2-startloc1,VWedge->Lcnt,VWedge->Rcnt,VWedge->Bcnt);
	if(VWedge->alignid){
	  Calign *p = alignment[VWedge->alignid - 1];
	  int LI = p->sites1[0];
	  int LJ = p->sites2[0];
	  int RI = p->sites1[p->numpairs - 1];
	  int RJ = p->sites2[p->numpairs - 1];
	  int N = Gmap[p->mapid1]->numsite[0];
	  int M = Gmap[p->mapid2]->numsite[0];
	  int mis = RI-LI+1+RJ-LJ+1 - 2 * p->numpairs;
	  double len1 = Gmap[p->mapid1]->site[0][RI] - Gmap[p->mapid1]->site[0][LI];
	  double len2 = p->orientation ? 
	    Gmap[p->mapid2]->site[0][M+1-LJ] - Gmap[p->mapid2]->site[0][M+1-RJ] :
	    Gmap[p->mapid2]->site[0][RJ] - Gmap[p->mapid2]->site[0][LJ];
	  fprintf(fp,",Yid=%d,Xid=%d:LI=%d,RI=%d,N=%d:LJ=%d,RJ=%d,M=%d,mis=%d,R=%0.3f(/100kb),score=%0.3f,logPV=%0.2f\n",
		 p->mapid1,p->mapid2,LI,RI,N,LJ,RJ,M,mis,100.0*mis/(len1+len2),p->score,p->logPV);
	} else
	  fprintf(fp,"\n");
      }
      if(map1 < map2){
	if(fabs(startloc2-startloc1) >= 20.0)
	  errcnt++;
	if(fabs(startloc2-startloc1) >= errsiz)
	  errsiz = fabs(startloc2-startloc1);
      }
    }
  }
  if(VERB){
    fprintf(fp,"%d nodes, %d edges: %d with cov>=%d\n",
	   nodecnt, edgecnt, highcnt, mincov);
    fflush(fp);
  }
  if(fp != stdout)
    FILEclose(fp);
}

void graph_displayA(int mincov, int verb, char *suffix) // for use with CnodeA,CedgeA : called from graph_display() if node == NULL && nodeA != NULL
{
  FILE *fp = stdout;
  if(suffix){
    char filename[1024];
    sprintf(filename,"%s_%s.graph",output_prefix,suffix);
    if((fp = fopen(filename,"w"))==NULL){
      int eno = errno;
      char *err = strerror(eno);
      printf("Unable to write to %s:errno=%d:%s\n",filename,eno,err);
      fflush(stdout);
    }
    if(VERB){
      printf("Generating %s (graph dump)\n",filename);
      fflush(stdout);
    }
  }

  int nodecnt=0;/* count of active nodes */
  int edgecnt=0;/* count of active edges */
  int highcnt=0;/* count of how many times edges with cov >= MinCoverage occured */
  int errcnt=0;/* count how many times edges with cov >= MinCoverage had true offset >= 20kb */
  double errsiz=0.0;
  int maxEcov = 0;
  int maxNcov = 0;

  if(!PosEdgeInit)
    PosEdgeA();

  for(int i=0; i < numnodes;i++){
    CnodeA *Vnode;
    int map1,node1;
    if(i < nummaps)
      node1 = nodemap[i];
    else {
      if(DEBUG) assert(CROSSBREAK);
      node1 = i;
    }
    Vnode = &nodeA[node1];
    map1 = Vnode->mapid;
    if(DEBUG>=2) assert(nodemap[nodeorder[map1]] == map1);
    if(DEBUG>=2) assert(Vnode->mapid == map1);
    if(DEBUG>=2 && i < nummaps){
      assert(node1 == map1);
      assert(nodeorder[map1] == i);
    } 

    if(Vnode->mark >= 2)/* no longer in topnodes[] */
      continue;
    nodecnt++;
    if(Vnode->cov > maxNcov)
      maxNcov = Vnode->cov;

    double startloc1 = Gmap[map1]->startloc;
    double reflen = startloc1;
    if(startloc1 > -1000.0){
      reflen = Gmap[map1]->endloc - startloc1;
      startloc1 += 0.5*Vnode->len;
    }
    if(VERB >= verb){
	//	fprintf(fp,"node%d:nodeid=%d,mapid=%d(id=%d),cov=%0.2f,PosEdge=%d,edgecnt=%d,loc=%0.3f(%0.3f..%0.3f,%0.3f..%0.3f),flip=%d,len=%0.3f(%0.3f),left=%0.1f,right=%0.1f,Ccntmax=%d,%d,trim=%d,%d,N=%d,Qscore=%0.3f",
	//	       nodeorder[map1],node1,map1,Gmap[map1]->id,Vnode->cov,Vnode->PosEdge,Vnode->edgecnt,startloc1,Gmap[map1]->startloc,Gmap[map1]->endloc,Gmap[map1]->startmatch,Gmap[map1]->endmatch,Gmap[map1]->flipped,
	//	       Vnode->len,reflen,Vnode->backward,Vnode->forward,Vnode->Lcntmax,Vnode->Rcntmax,Vnode->trimL,Vnode->trimR,Gmap[map1]->numsite[0],Gmap[map1]->Qscore);
      fprintf(fp,"N%d:",nodeorder[map1]);
      fprintf(fp,"nodeid=%d,",node1);
      fprintf(fp,"(mapid=%d,id=%lld),",map1,Gmap[map1]->id);
      fprintf(fp,"cov=%0.2f,",Vnode->cov);
      fprintf(fp,"PosEdge=%d,",Vnode->PosEdge);
      fprintf(fp,"edgecnt=%d,",Vnode->edgecnt);
      fprintf(fp,"loc=%0.3f",startloc1);
      fprintf(fp,"(%0.3f..",Gmap[map1]->startloc);
      fprintf(fp,"%0.3f,",Gmap[map1]->endloc);
      fprintf(fp,"%3f..",Gmap[map1]->startmatch);
      fprintf(fp,"%3f),",Gmap[map1]->endmatch);
      fprintf(fp,"flip=%d,",Gmap[map1]->flipped);
      fprintf(fp,"len=%0.3f",Vnode->len);
      fprintf(fp,"(%0.3f),",reflen);
      fprintf(fp,"left=%0.1f,right=%0.1f,",Vnode->backward,Vnode->forward);
      fprintf(fp,"trim=%d,%d,",Vnode->trimL,Vnode->trimR);
      fprintf(fp,"N=%d,",Gmap[map1]->numsite[0]);
      fprintf(fp,"Qscore=%0.3f",Gmap[map1]->Qscore);

      if(Gmap[map1]->name)
	fprintf(fp,":%s\n",Gmap[map1]->name);
      else
	fprintf(fp,"\n");
    }
    for(int vw = 0; vw < Vnode->edgecnt; vw++){
      CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
      CnodeA *Wnode = &nodeA[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
      if(VWedge->valid < 0 || VWedge->cov <= 0.0)
	continue;
      if(VWedge->cov > maxEcov)
	maxEcov = VWedge->cov;

      int map2 = node[(VWedge->Lmap == Vnode->mapid) ? VWedge->Rmap : VWedge->Lmap].mapid;
      if(DEBUG>=2) assert(nodemap[nodeorder[map2]] == map2);
      if(DEBUG>=2 && map1 == map2){
	printf("map1=%d,Vnode=%d(N%d),map2=%d,VWedge=(%d,%d,%0.3f,%d,%d)\n",
	       map1,Vnode->mapid,nodeorder[Vnode->mapid],map2,VWedge->Lmap,VWedge->Rmap,VWedge->distance,VWedge->Lflipped,VWedge->Rflipped);
	fflush(stdout);
	assert(map1 != map2);
      }

      if(map1 < map2)
	edgecnt++;
      if(floor(VWedge->cov+0.5) < mincov)
	continue;
      if(map1 < map2)
	highcnt++;
      double startloc2 = Gmap[map2]->startloc;
      if(startloc2 > -1000.0)
	startloc2 += 0.5*Wnode->len;
      if(VERB >= verb){
	fprintf(fp," %cEdge%d = (N%d,N%d,%0.3f,%d,%d) to N%d(mapid=%d,id=%lld):valid=%d,cov=%0.2f,wt=%0.3f,loc=%0.3f(delta=%0.3f)",
	       VWedge->alignid ? ' ' : 'S', vw, nodeorder[node[VWedge->Lmap].mapid], nodeorder[node[VWedge->Rmap].mapid], VWedge->distance, VWedge->Lflipped,VWedge->Rflipped,
		nodeorder[map2], map2, Gmap[map2]->id,VWedge->valid, VWedge->cov, 0.0, startloc2, startloc2-startloc1);
	if(VWedge->alignid){
	  Calign *p = alignment[VWedge->alignid - 1];
	  int LI = p->sites1[0];
	  int LJ = p->sites2[0];
	  int RI = p->sites1[p->numpairs - 1];
	  int RJ = p->sites2[p->numpairs - 1];
	  int N = Gmap[p->mapid1]->numsite[0];
	  int M = Gmap[p->mapid2]->numsite[0];
	  int mis = RI-LI+1+RJ-LJ+1 - 2 * p->numpairs;
	  double len1 = Gmap[p->mapid1]->site[0][RI] - Gmap[p->mapid1]->site[0][LI];
	  double len2 = p->orientation ? 
	    Gmap[p->mapid2]->site[0][M+1-LJ] - Gmap[p->mapid2]->site[0][M+1-RJ] :
	    Gmap[p->mapid2]->site[0][RJ] - Gmap[p->mapid2]->site[0][LJ];
	  fprintf(fp,",Yid=%d,Xid=%d:LI=%d,RI=%d,N=%d:LJ=%d,RJ=%d,M=%d,mis=%d,R=%0.3f(/100kb),score=%0.3f,logPV=%0.2f\n",
		 p->mapid1,p->mapid2,LI,RI,N,LJ,RJ,M,mis,100.0*mis/(len1+len2),p->score,p->logPV);
	} else
	  fprintf(fp,"\n");
      }
      if(map1 < map2){
	if(fabs(startloc2-startloc1) >= 20.0)
	  errcnt++;
	if(fabs(startloc2-startloc1) >= errsiz)
	  errsiz = fabs(startloc2-startloc1);
      }
    }
  }
  if(VERB){
    fprintf(fp,"%d nodes, %d edges: %d with cov>=%d\n",
	   nodecnt, edgecnt, highcnt, mincov);
    fflush(fp);
  }
  if(fp != stdout)
    FILEclose(fp);
}

#define SKIP 4 /* ignore last SKIP aligned fragments on either end of an alignment when computing fragment coverage to detect chimeric molecules */

/** compute local coverage profiles of each node (map) using all other maps that align with it (works best with local alignments).
   Then look for:
   1. Shallow coverage at either end : trim the poor quality ends of the maps (just mark with trimL and trimR).
   2. Shallow internal region or site : Chimeric maps are marked for deletion (mark = 2).
   3. Shallow average coverage below TrimCoverage : mark for deletion (mark = 2). 

   4. Delete any edges with fewer than 1+ MinSites/2 aligned sites that are within trimL .. trimR (for both maps)

 */

long long graph_nodetrim(int verb)
{
  int shallowcnt = 0;// count of maps with low average coverage
  int shallowchimcnt = 0;// count of maps with low average coverage that were chimeric 
  int chimcnt = 0;//count of maps deleted as chimeric maps

  int maxsites = 0;
  for(int i=0; i < numtopnodes; i++){
    CnodeA *Vnode = (CnodeA *) topnodes[i];
    int M = Gmap[Vnode->mapid]->numsite[0];
    if(M > maxsites)
      maxsites = M;
  }

  int *fragcov = new int[maxsites+1];
  int *sitecov = new int[maxsites+1];

  for(int i = 0; i < numtopnodes; i++){
    CnodeA *Vnode = (CnodeA *) topnodes[i];
    int M = Gmap[Vnode->mapid]->numsite[0];
    //    int *sitecov = Vnode->sitecov = new int[M+1];
    int k;
    for(k=0; k <= M; k++)
      fragcov[k] = sitecov[k] = 0;

    for(int vw = 0; vw < Vnode->edgecnt; vw++){
      CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
      Calign *VWalign = alignment[VWedge->alignid - 1];
      int L,R;
      if(VWalign->mapid1 == Vnode->mapid){
	L = VWalign->sites1[0];
	R = VWalign->sites1[VWalign->numpairs - 1];
	for(int t = 0; t < VWalign->numpairs; t++)
	  sitecov[VWalign->sites1[t]]++;
      } else {
	if(DEBUG>=2) assert(VWalign->mapid2 == Vnode->mapid);
	if(VWalign->orientation){
	  L = M+1 - VWalign->sites2[VWalign->numpairs - 1];
	  R = M+1 - VWalign->sites2[0];
	  for(int t = 0; t < VWalign->numpairs; t++)
	    sitecov[M+1 - VWalign->sites2[t]]++;
	} else {
	  L = VWalign->sites2[0];
	  R = VWalign->sites2[VWalign->numpairs - 1];
	  for(int t = 0; t < VWalign->numpairs; t++)
	    sitecov[VWalign->sites2[t]]++; 
	}
      }
      if(DEBUG>=2) assert(R >= L);
      for(k = L+SKIP; k < R-SKIP; k++)
	fragcov[k]++;
    }

    /* compute average site coverage */
    int sum=0;
    for(k = 1; k <= M; k++)
      sum += sitecov[k];
    double meancov = ((double)sum)/M;

    /* adjust average site coverage by ignoring sites below TrimCoverageRatio * meancov */
    int origcnt = M;
    while(1){
      int cnt;
      for(sum=0.0,cnt =0, k=1; k <= M; k++)
	if(sitecov[k] >= TrimCoverageRatio * meancov){
	  sum += sitecov[k];
	  cnt++;
	}
      if(cnt >= origcnt)
	break;
      if(DEBUG>=2) assert(cnt > 0);
      meancov = ((double)sum)/cnt;
      origcnt = cnt;
    }

    if(meancov < TrimCoverage){
      if(VERB>=2 && verb){
	if(Gmap[Vnode->mapid]->name)
	  printf("N%d (%s) deleted due to mean coverage = %0.1f < TrimCoverage(=%0.1f)\n",
		 nodeorder[Vnode->mapid], Gmap[Vnode->mapid]->name, meancov, TrimCoverage);
	else
	  printf("N%d deleted due to mean coverage = %0.1f < TrimCoverage(=%0.1f)\n",
		 nodeorder[Vnode->mapid], meancov, TrimCoverage);
	fflush(stdout);
      }
      shallowcnt++;
      if(Gmap[Vnode->mapid]->name && (strstr(Gmap[Vnode->mapid]->name,"R_") || strstr(Gmap[Vnode->mapid]->name,"C_")))
	shallowchimcnt++;
      Vnode->mark = 2;
      /* mark all edges out of this node for deletion */
      for(int vw = 0; vw < Vnode->edgecnt; vw++)
	edge_delete(&edgeA[Vnode->edgeid[vw]],Vnode);
      continue;
    }

    /* trim off left end with coverage below TrimCoverageRatio * meancov */
    Vnode->trimL = 0;
    for(k = 1; k <= M; k++)
      if(sitecov[k] >= TrimCoverageRatio * meancov)
	break;
    if(DEBUG>=2) assert(k <= M);
    Vnode->trimL = k-1;
    /* trim off right end with coverage below TrimCoverage */
    for(k = M; k >= 1; k--)
      if(sitecov[k] >= TrimCoverageRatio * meancov)
	break;
    if(DEBUG>=2) assert(k >= 1);
    Vnode->trimR = k+1;

    if(DEBUG>=2) assert(Vnode->trimL < Vnode->trimR);
    
    /* analyze fragcov[] to see if any value is below ChimCoverageRatio times largest values on either side */
    /* first find the largest peak */
    int Peak = 1;
    for(k = 2; k < M; k++)
      if(fragcov[k] > fragcov[Peak])
	Peak = k;
    if(fragcov[Peak] >= TrimCoverage){
      if(Peak < M-2){ /* look for valley to right of peak */
	int LWM= Peak+1, HWM = Peak+1;
	for(k = LWM+1; k < M; k++){
	  if(fragcov[k] < fragcov[LWM])
	    LWM = HWM = k;
	  else if(fragcov[k] > fragcov[HWM])
	    HWM = k;
	  if(fragcov[HWM] >= TrimCoverage &&
	     fragcov[LWM] < ChimCoverageRatio * fragcov[Peak] &&
	     fragcov[LWM] < ChimCoverageRatio * fragcov[HWM])
	    break;
	}
	if(k < M){/* found chimeric valley */
	  if(VERB>=2 && verb){
	    for(;++k < M;)
	      if(fragcov[k] > fragcov[HWM])
		HWM = k;
	    if(Gmap[Vnode->mapid]->name)
	      printf("N%d (%s) deleted as chimeric due to fragcov[%d]=%d,fragcov[%d]=%d,fragcov[%d]=%d,M=%d\n",
		     nodeorder[Vnode->mapid], Gmap[Vnode->mapid]->name, Peak, fragcov[Peak], LWM, fragcov[LWM], HWM, fragcov[HWM], M);
	    else
	      printf("N%d deleted as chimeric due to fragcov[%d]=%d,fragcov[%d]=%d,fragcov[%d]=%d,M=%d\n",
		     nodeorder[Vnode->mapid], Peak, fragcov[Peak], LWM, fragcov[LWM], HWM, fragcov[HWM], M);
	    fflush(stdout);
	  }
	
	  Vnode->mark = 2;
	}
      }
      if(!Vnode->mark && Peak > 2){/* look for valley to left of peak */
	int LWM= Peak-1, HWM = Peak-1;
	for(k = LWM-1; k >= 1; k--){
	  if(fragcov[k] < fragcov[LWM])
	    LWM = HWM = k;
	  else if(fragcov[k] > fragcov[HWM])
	    HWM = k;
	  if(fragcov[HWM] >= TrimCoverage &&
	     fragcov[LWM] < ChimCoverageRatio * fragcov[Peak] &&
	     fragcov[LWM] < ChimCoverageRatio * fragcov[HWM])
	    break;
	}
	if(k >= 1){/* found chimeric valley */
	  if(VERB>=2 && verb){
	    for(;--k >= 1;)
	      if(fragcov[k] > fragcov[HWM])
		HWM = k;
	    if(Gmap[Vnode->mapid]->name)
	      printf("N%d (%s) deleted as chimeric due to fragcov[%d]=%d,fragcov[%d]=%d,fragcov[%d]=%d,M=%d\n",
		     nodeorder[Vnode->mapid], Gmap[Vnode->mapid]->name, HWM, fragcov[HWM], LWM, fragcov[LWM], Peak, fragcov[Peak], M);
	    else
	      printf("N%d deleted as chimeric due to fragcov[%d]=%d,fragcov[%d]=%d,fragcov[%d]=%d,M=%d\n",
		     nodeorder[Vnode->mapid], HWM, fragcov[HWM], LWM, fragcov[LWM], Peak, fragcov[Peak], M);
	    fflush(stdout);
	  }
	  Vnode->mark = 2;
	}
      }
    }
    if(Vnode->mark == 2){
      chimcnt++;

      /* mark all edges out of this node for deletion */
      for(int vw = 0; vw < Vnode->edgecnt; vw++)
	edge_delete(&edgeA[Vnode->edgeid[vw]],Vnode);
    } else  if(VERB>=2 && verb && ((Gmap[Vnode->mapid]->name && (strstr(Gmap[Vnode->mapid]->name,"R_") || strstr(Gmap[Vnode->mapid]->name,"C_"))) ||
				nodeorder[Vnode->mapid] == 788 /* 1255 */ /* 462*/)){
      printf("nodetrim:failed to detect chimeric map N%d (%s): average site coverage=%0.1f, peak fragcov[%d]=%d:M=%d,trimL=%d,trimR=%d\n",
	     nodeorder[Vnode->mapid], Gmap[Vnode->mapid]->name, meancov, Peak, fragcov[Peak], M, Vnode->trimL, Vnode->trimR);
      for(int t = 1; t <= M; t++)
	printf("t=%d:fragcov[t]=%d,sitecov[t]=%d\n",t,fragcov[t],sitecov[t]);
      //	  exit(1);
      fflush(stdout);
    }
    
    //    delete [] Vnode->sitecov; 
    //    Vnode->sitecov = 0;
  }

  delete [] fragcov;
  delete [] sitecov;

  int EdgeTrimCnt = 0;/* number of edges trimmed due to trimL..trimR */
  int minsites = 1 + MinSites/2;
  for(int i = 0; i < numtopnodes; i++){
    CnodeA *Vnode = (CnodeA *) topnodes[i];
    if(Vnode->mark >= 2)
      continue;
    int M = Gmap[Vnode->mapid]->numsite[0];
    for(int vw = 0; vw < Vnode->edgecnt; vw++){
      CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
      if(VWedge->valid < 0)
	continue;
      CnodeA *Wnode = &nodeA[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
      int N = Gmap[Wnode->mapid]->numsite[0];
      if(Vnode->trimL <= 0 && Vnode->trimR >= M+1 && Wnode->trimL <= 0 && Wnode->trimR >= N+1)
	continue;
      if(!VWedge->alignid)
	continue;
      Calign *align = alignment[VWedge->alignid - 1];

      int K = align->numpairs;
      int cnt = 0;
      int trimL1,trimR1,trimL2,trimR2;
      if(align->mapid1 == Vnode->mapid) {
	if(DEBUG>=2) assert(align->mapid2 == Wnode->mapid);
	trimL1 = Vnode->trimL;
	trimR1 = Vnode->trimR;
	trimL2 = align->orientation ? N+1-Wnode->trimR : Wnode->trimL;
	trimR2 = align->orientation ? N+1-Wnode->trimL : Wnode->trimR;
      } else {
	if(DEBUG>=2) assert(align->mapid2 == Vnode->mapid);
	if(DEBUG>=2) assert(align->mapid1 == Wnode->mapid);
	trimL2 = align->orientation ? M+1-Vnode->trimR : Vnode->trimL;
	trimR2 = align->orientation ? M+1-Vnode->trimL : Vnode->trimR;
	trimL1 = Wnode->trimL;
	trimR1 = Wnode->trimR;
      }

      for(int k = 0; k < K; k++){
	int i = align->sites1[k];
	if(i <= trimL1 || i >= trimR1)
	  continue;
	int j = align->sites2[k];
	if(j <= trimL2 || j >= trimR2)
	  continue;
	cnt++;
      }
      if(cnt < minsites){
	VWedge->valid = -1;
	EdgeTrimCnt++;
      }
    }
  }

  int edgecnt = 0;// count of edges trimmed due to Ccnt 

  unsigned int finalcnt;
  long long Ecnt = GCedges(finalcnt,edgeA);
  int Ncnt = GCnodesA();

  if(VERB /*&& (chimcnt > 0 || shallowcnt>0 || edgecnt>0)*/){
    printf("graph_nodetrim(%d): deleted %d shallow maps(chimeric=%d) and %d chimeric maps and %d edges with Ccnt=0, %d edges due to trimL,trimR(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
	   verb,shallowcnt,shallowchimcnt,chimcnt,edgecnt,EdgeTrimCnt,finalcnt,numtopnodes,wtime(),mtime());
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
  if(DEBUG>=2){
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edgeA[i].valid < 0)
	continue;
      assert(edgeA[i].mark == 0);
    }
    size_t edgecnt = 0;
    for(int i = 0; i < numnodes; i++){
      CnodeA *Vnode = &nodeA[i];
      if(Vnode->mark >= 2){
	assert(Vnode->edgecnt <= 0);
	continue;
      }
      edgecnt += Vnode->edgecnt;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	assert(1 <= VWedgeid && VWedgeid <= Numedges);
	assert(edgeA[VWedgeid].mark == 0);
      }
    }    
    if(DEBUG>=3){
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edgeA[Vnode->edgeid[vw]].mark++;
      }    
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  CedgeA *VWedge = &edgeA[Vnode->edgeid[vw]];
	  if(VWedge->mark != 2){
	    printf("i=%d/%d: mapid=%d, id=%lld : vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d:align->logPV= %0.2f,or=%d\n",
		   i,numnodes,Vnode->mapid,gmap[Vnode->mapid]->id, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark,
		   VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped,alignment[VWedge->alignid - 1]->logPV, alignment[VWedge->alignid - 1]->orientation);
	    fflush(stdout);
	    assert(VWedge->mark == 2);
	  }
	}
      }    
      for(int i = 0; i < numnodes; i++){
	CnodeA *Vnode = &nodeA[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edgeA[Vnode->edgeid[vw]].mark = 0;
      }    
      /* check that topnodes[0..numtopnodes-1] reach all nodeA[] that do not have mark >= 2 */
      for(int i = 0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	assert(Vnode->mark == 0);
	Vnode->mark = 2;
      }
      for(int i = 0; i < numnodes; i++)
	assert(nodeA[i].mark >= 2);
      for(int i = 0; i < numtopnodes; i++){
	CnodeA *Vnode = (CnodeA *) topnodes[i];
	Vnode->mark = 0;
      }
    }
    if((edgecnt % 2) != 0){
      printf("sum of nodeA[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",
	     edgecnt,numnodes,Numedges);
      fflush(stdout);
    }
  }

  return Ncnt+Ecnt;
}

/** see graph_edgetrim() for details */
class Cneighbor {
 public:
  int nodeid;
  int flip;
  double Edistance;/* short for edge->distance */
};

/* sort Cneighbor in increasing order of nodeid */
int neighborinc(Cneighbor *p1, Cneighbor *p2)
{
  int id1 = p1->nodeid;
  int id2 = p2->nodeid;
  
  return (id1 > id2) ? 1 : (id1 < id2) ? -1 : 0;
}

Cneighbor *VLmem = 0;
Cneighbor **VL = 0;
Cneighbor **VR = 0;
int *VLcnt = 0;
int *VRcnt = 0;

static void compute_neighbors(Cnode *Vnode, Cneighbor *VL, Cneighbor *VR, int &VLcnt, int &VRcnt)
{
  VLcnt = VRcnt = 0;
  
  for(register int vl = 0; vl < Vnode->PosEdge; vl++){
    register Cedge *VLedge = &edge[Vnode->edgeid[vl]];
    register Cneighbor *pn = &VR[VRcnt++];
    pn->nodeid = (VLedge->Lmap == Vnode->nodeid) ? VLedge->Rmap : VLedge->Lmap;
    pn->flip = (VLedge->Lflipped | VLedge->Rflipped);
    pn->Edistance = VLedge->distance;
  }
  for(register int vr = Vnode->PosEdge; vr < Vnode->edgecnt; vr++){
    register Cedge *VRedge = &edge[Vnode->edgeid[vr]];
    register Cneighbor *pn = &VL[VLcnt++];
    pn->nodeid = (VRedge->Lmap == Vnode->nodeid) ? VRedge->Rmap : VRedge->Lmap;
    pn->flip = (VRedge->Lflipped | VRedge->Rflipped);
    pn->Edistance = VRedge->distance;
  }

  if(DEBUG>=2) assert(VLcnt <= Vnode->edgecnt);
  if(DEBUG>=2) assert(VRcnt <= Vnode->edgecnt);

  /* sort all the node lists so we can compute the number of common nodes */
  qsort(VL,VLcnt,sizeof(Cneighbor),(intcmp *)neighborinc);
  qsort(VR,VRcnt,sizeof(Cneighbor),(intcmp *)neighborinc);
}

static double varF3,varF4,varSD;

static void compute_edgecnt(Cnode *Vnode, int Vflipped, 
			    Cnode *Wnode, int Wflipped, 
			    double VWdistance,
			    int &VWLcnt, int &VWRcnt, int &VWBcnt, int &VWCcnt,
			    int edgebridge, int verb,
			    Cneighbor *VL, Cneighbor *VR, Cneighbor *WL, Cneighbor *WR,
			    int VLcnt, int VRcnt, int WLcnt, int WRcnt)
{

  /* to compute Lcnt,Rcnt,Bcnt collect all nodes to the left(right) of V and W (in their current orientation) */
  if(!Vflipped){/* swap VL & VR */
    Cneighbor *tmp = VL; VL = VR; VR = tmp;
    int t = VLcnt; VLcnt = VRcnt; VRcnt = t;
  }
  if(!Wflipped){/* swap WL & WR */
    Cneighbor *tmp = WL; WL = WR; WR = tmp;
    int t = WLcnt; WLcnt = WRcnt; WRcnt = t;
  }

  int v,w;
  int Lcnt = 0;
  int VWflip = Vflipped ^ Wflipped;
  for(v=w=0;v < VLcnt && w < WLcnt;){
    if(VL[v].nodeid == WL[w].nodeid){
      if(!EDGECHECK)
	Lcnt++;
      else if((VL[v].flip ^ WL[w].flip) == VWflip){
	double x = VWdistance + VL[v].Edistance;
	double y = WL[w].Edistance;
	double var = varF3 + varSD*(x + y);
	double err = x - y;
	if(err * err  < var)
	  Lcnt++;
      }
      v++;
      w++;
    } else if(VL[v].nodeid < WL[w].nodeid)
      v++;
    else
      w++;
  }
      
  int Rcnt = 0;
  for(v=w=0;v < VRcnt && w < WRcnt;){
    if(VR[v].nodeid == WR[w].nodeid){
      if(!EDGECHECK)
	Rcnt++;
      else if((VR[v].flip ^ WR[w].flip) == VWflip){
	double x = VWdistance + WR[w].Edistance;
	double y = VR[v].Edistance;
	double var = varF3 + varSD*(x + y);
	double err = x - y;
	if(err * err < var)
	  Rcnt++;
      }
      v++;
      w++;
    } else if(VR[v].nodeid < WR[w].nodeid)
      v++;
    else
      w++;
  }
      
  int Bcnt = 0;
  if(edgebridge){
    /* now count the number of edges between VL[] and WR[] that have consistent orientation and length */
    /* first mark nodes on WR[w] with mark = w+1 */
    if(DEBUG>=2)
      for(w = 0; w < WRcnt; w++)
	assert(node[WR[w].nodeid].mark == 0);

    for(w = 0; w < WRcnt; w++){
      if(DEBUG>=2) assert(node[WR[w].nodeid].mark == 0);
      node[WR[w].nodeid].mark = w+1;
    }
    for(v = 0; v < VLcnt; v++){
      Cneighbor *pn = &VL[v];
      Cnode *Xnode = &node[pn->nodeid];
      int Xflipped = pn->flip ^ Vflipped;
      /* If Xflipped==1: look for neighbors of X with -ve node distance in increasing order of distance
	 If Xflipped==0: look for neighers of X with +ve node distance in increasing order of distance */
      for(int xt = (Xflipped ? 0 : Xnode->PosEdge); xt < (Xflipped ? Xnode->PosEdge : Xnode->edgecnt); xt++) {
	Cedge *XTedge = &edge[Xnode->edgeid[xt]];
	Cnode *Tnode = &node[(XTedge->Lmap == Xnode->nodeid) ? XTedge->Rmap : XTedge->Lmap];
	int Tflipped = Xflipped ^ (XTedge->Lflipped | XTedge->Rflipped);
	if(Tnode->mark){
	  int w = Tnode->mark-1;
	  if(DEBUG>=2 && Tflipped == (WR[w].flip ^ Wflipped) && !(Xflipped ? (xt < Xnode->PosEdge) : (xt >= Xnode->PosEdge))){
	    printf("Vnode=N%d%c,Wnode=N%d%c,Xnode=N%d%c,Tnode=N%d%c,XTedge=(N%d,N%d,%0.3f,%d,%d),xt=%d,Xnode->PosEdge=%d,Tnode->mark=%d,WR[%d].flip=%d,WR[%d].nodeid=N%d\n",
		   nodeorder[Vnode->mapid],Vflipped ? '\'':' ',
		   nodeorder[Wnode->mapid],Wflipped ? '\'':' ',
		   nodeorder[Xnode->mapid],Xflipped ? '\'':' ',
		   nodeorder[Tnode->mapid],Tflipped ? '\'':' ',
		   nodeorder[node[XTedge->Lmap].mapid],nodeorder[node[XTedge->Rmap].mapid],XTedge->distance,XTedge->Lflipped,XTedge->Rflipped,xt,Xnode->PosEdge,
		   Tnode->mark,w,WR[w].flip, w,nodeorder[node[WR[w].nodeid].mapid]);
	    fflush(stdout);
	    assert(Xflipped ? (xt < Xnode->PosEdge) : (xt >= Xnode->PosEdge));
	  }
	  if(!EDGECHECK)
	    Bcnt++;
	  else if(Tflipped == (WR[w].flip ^ Wflipped)){
	    double x = VWdistance + VL[v].Edistance + WR[w].Edistance;
	    double y = XTedge->distance;
	    double var = varF4 + varSD*(x + y);
	    double err = x - y;
	    if(err * err < var){
	      Bcnt++;
	    }
	  }
	}
      }
    }
    for(w = 0; w < WRcnt; w++)
      node[WR[w].nodeid].mark = 0;

    if(DEBUG>=2){
      for(int k=0; k < numtopnodes; k++){
	Cnode *Tnode = topnodes[k];
	assert(Tnode->mark == 0);
      }
    }
  }

  if(VERB>=2 && verb>=2 && ((nodeorder[Vnode->mapid]==958 && nodeorder[Wnode->mapid] == 481)|| 
			 (nodeorder[Vnode->mapid]==481 && nodeorder[Wnode->mapid] == 958))){
    printf("Vnode=N%d(flip=%d),Wnode=N%d(flip=%d),VLcnt=%d,VRcnt=%d,WLcnt=%d,WRcnt=%d,Lcnt=%d,Rcnt=%d,Bcnt=%d\n",
	   nodeorder[Vnode->mapid],Vflipped,nodeorder[Wnode->mapid],Wflipped,VLcnt,VRcnt,WLcnt,WRcnt,Lcnt,Rcnt,Bcnt);
    for(int k=0; k < VLcnt; k++)
      printf("VL[%d]=%d(N%d)\n",k,VL[k].nodeid,nodeorder[node[VL[k].nodeid].mapid]);
    for(int k=0; k < VRcnt; k++)
      printf("VR[%d]=%d(N%d)\n",k,VR[k].nodeid,nodeorder[node[VR[k].nodeid].mapid]);
    for(int k=0; k < WLcnt; k++)
      printf("WL[%d]=%d(N%d)\n",k,WL[k].nodeid,nodeorder[node[WL[k].nodeid].mapid]);
    for(int k=0; k < WRcnt; k++)
      printf("WR[%d]=%d(N%d)\n",k,WR[k].nodeid,nodeorder[node[WR[k].nodeid].mapid]);
    fflush(stdout);
  }

  VWLcnt = Lcnt;
  VWRcnt = Rcnt;
  VWBcnt = Bcnt;
  VWCcnt = min(Lcnt,Rcnt);
  if(edgebridge >= 2)
    VWCcnt += Bcnt;
}

static void compute_edgecnt_alloc()
{
  // allocate global VL,VR,WL,WR for each node (if not already done) */
  if(VL==0){
    VL = new Cneighbor*[2*numnodes];
    //    VR = new Cneighbor*[numnodes];
    VR = &VL[numnodes];
    VLcnt = new int[2*numnodes];
    //    VRcnt = new int[numnodes];
    VRcnt = &VLcnt[numnodes];

    for(int i = 0; i < numnodes; i++)
      VL[i] = 0;
    
    size_t VLmemsiz = 0;
    for(int i = 0; i < numtopnodes; i++)   
      VLmemsiz += 2 * topnodes[i]->edgecnt;
    VLmem = new Cneighbor[VLmemsiz];

    size_t VLmemcnt = 0;

    for(int i = 0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      int V = Vnode->nodeid;
      VL[V] = &VLmem[VLmemcnt];
      VLmemcnt += Vnode->edgecnt;
      VR[V] = &VLmem[VLmemcnt];
      VLmemcnt += Vnode->edgecnt;
    }
    
    if(DEBUG) assert(VLmemcnt == VLmemsiz);
  }
}

static int graph_edgecnt(int verb)
{
  
  varF3 = 3.0*SM*SM * (MAXERR*MAXERR);
  varSD = SD[0]*SD[0] * (MAXERR*MAXERR);
  varF4 = 4.0*SM*SM * (MAXERR*MAXERR);

  if(!PosEdgeInit)
    PosEdge();

  if(DEBUG>=2){
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      assert(Vnode->mark == 0);
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	VWedge->Ccnt = VWedge->Lcnt = VWedge->Rcnt = VWedge->Bcnt = -1;
	assert(VWedge->valid == 1);
      }
    }
  }

  compute_edgecnt_alloc();
  int block=16;

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {
    #pragma omp for schedule(dynamic,block)
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      int V = Vnode->nodeid;
      compute_neighbors(Vnode,VL[V],VR[V],VLcnt[V],VRcnt[V]);
    }

    #pragma omp for schedule(dynamic,block)
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      int V = Vnode->nodeid;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(VWedge->Lmap != V)/* will be handled when Vnode matches Lmap */
	  continue;
	int W = VWedge->Rmap;
	Cnode *Wnode = &node[W];
	int Vflipped = VWedge->Lflipped;
	int Wflipped = VWedge->Rflipped;

	compute_edgecnt(Vnode,Vflipped,Wnode,Wflipped, VWedge->distance, VWedge->Lcnt, VWedge->Rcnt, VWedge->Bcnt, VWedge->Ccnt, EDGEBRIDGE>= 2 ? 2 : 0, verb,VL[V],VR[V],VL[W],VR[W],VLcnt[V],VRcnt[V],VLcnt[W],VRcnt[W]);
      }
    }

    /* now compute largest value of Ccnt for each node's left & right edges */
    #pragma omp for schedule(dynamic,block)
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      Vnode->Lcntmax = Vnode->Rcntmax = 0;
      for(int vw = 0; vw < Vnode->PosEdge; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(DEBUG>=2) assert(VWedge->Ccnt >= 0);
	if(VWedge->Ccnt > Vnode->Lcntmax)
	  Vnode->Lcntmax = VWedge->Ccnt;
      }
      for(int vw = Vnode->PosEdge; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(DEBUG>=2) assert(VWedge->Ccnt >= 0);
	if(VWedge->Ccnt > Vnode->Rcntmax)
	  Vnode->Rcntmax = VWedge->Ccnt;
      }
    }
  }// omp parallel

  /* NOTE global VL,VR,WL,WR are deleted in compute_redundancy(), the last function to call compute_edgecnt() */

  return 0;
}

/** For each edge connecting nodes A -> B, check if it is spurious by counting (in graph_edgecnt()):
   1. The number of shared left nodes (C->A and C->B) : Lcnt
   2. The number of shared right nodes (A->D and B->D) : Rcnt
   3. The number of bridged edges (C->A and B->D and C->D) : Bcnt (If EDGEBRIDGE >= 1)

   The value of min(Lcnt,Rcnt) is saved as Ccnt.

   The largest value of Ccnt for any left edge connected to node Vnode is saved as Vnode->Lcntmax
   The largest value of Ccnt for any right edge connected to node Vnode is saved as Vnode->Rcntmax

   If Ccnt is small compared to min(Lcntmax,Rcntmax) for EITHER/BOTH (see EDGETRIM=1,2) nodes, then this edge it will be deleted. 
   This should handle chimeric edges (if EDGETRIM=1) since the one non-chimeric end should have a large min(Lcntmax,Rcntmax) in most cases), while avoiding
   avoiding edges spanning fragile nicking sites (which should have min(Lcntmax,Rcntmax) small on both sides). EDGETRIM=2 is better at avoiding edges spanning
   fragile nicking sites but will not handle chimeric edges (but these are handled more accurately in chimtrim())

  NOTE : Will trim away true edges that have the longest distance since there are no longer edges to confirm them. This is a reasonable loss
      since the dominant path should not be affected. However this means this functions cannot be called more than once, otherwise
      each time the next most distant true edges (at each node) will be trimmed away.

   1. Delete all nodes with max(Vnode->Lcntmax,Vnode->Rcntmax) < min(MinCcntmax,min(Vnode->forward,Vnode->backward)*EdgeTrimRatio)
      This deletes a surprisingly large number of nodes in the complete drosophila data set (even after filtering out maps under 175kb, 14 sites). 
      Probably caused by MaxRelCoverage

      Also may have the problem that it delete nodes near the ends of linear contigs.

		
*/
int graph_edgetrim(int verb)
{
  
  if(VERB){
    printf("Using %d threads\n",numthreads);
    fflush(stdout);
  }

  graph_edgecnt(verb);

  unsigned int edgecnt = 0;// count of edges trimmed 
  int chimcnt = 0;

  if(DEBUG>=2)
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	assert(VWedge->mark == 0);
      }
    }

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {
    unsigned int myedgecnt = 0;// count of edges trimmed 
    int mychimcnt = 0;

    /* now threshold edges by comparing Ccnt with Lcntmax,Rcntmax of the connected nodes based on EDDGETRIM method */
    #pragma omp for schedule(dynamic,16)
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(VWedge->Lmap != Vnode->mapid)/* will be handled when Vnode matches Lmap */
	  continue;
	if(DEBUG>=2) assert(VWedge->valid == 1);
	Cnode *Wnode = &node[VWedge->Rmap];
	/* Note : Edges with Ccnt==0 are probably chimeric edges, but removing them based on a ratio will remove them
	   whenever Lcntmax,Rcntmax >= 1. Hence we use max(1,Ccnt) instead */
	if(EDGETRIM==0 ? (max(1,VWedge->Ccnt) < min(max(Vnode->Lcntmax,Vnode->Rcntmax),max(Wnode->Lcntmax,Wnode->Rcntmax)) * EdgeTrimRatio) :
	   EDGETRIM==1 ? (max(1,VWedge->Ccnt) < min(Vnode->Lcntmax,Vnode->Rcntmax) * EdgeTrimRatio || max(1,VWedge->Ccnt) < min(Wnode->Lcntmax,Wnode->Rcntmax) * EdgeTrimRatio) :
	   (max(1,VWedge->Ccnt) < min(Vnode->Lcntmax,Vnode->Rcntmax) * EdgeTrimRatio && max(1,VWedge->Ccnt) < min(Wnode->Lcntmax,Wnode->Rcntmax) * EdgeTrimRatio) ){

          #pragma omp critical(valid)
          {
	    if(VERB>=2 &&  (verb /* || nodeorder[Vnode->mapid]==462 || nodeorder[Wnode->mapid]==462 ||
				 nodeorder[Vnode->mapid] == 856 || nodeorder[Wnode->mapid]==856 ||
				 (nodeorder[Vnode->mapid] == 1528 && nodeorder[Wnode->mapid] == 134) ||
				 (nodeorder[Vnode->mapid] == 134 && nodeorder[Wnode->mapid] == 1528)*/)){
	      printf("Trimming edge=(N%d,N%d,%0.3f,%d,%d):Lcnt=%d,Rcnt=%d,Ccnt=%d,Vcntmax=(%d,%d),Wcntmax=(%d,%d),mark=%d\n",
		     nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,
		     VWedge->Lflipped,VWedge->Rflipped,VWedge->Lcnt,VWedge->Rcnt,
		     VWedge->Ccnt,Vnode->Lcntmax,Vnode->Rcntmax,Wnode->Lcntmax,Wnode->Rcntmax,VWedge->mark);
	      fflush(stdout);
	    }
	    VWedge->valid = -1;
	  }

	  if(DEBUG>=2) assert(VWedge->mark == 0);
	  myedgecnt++;
	}
	if(DEBUG>=2) VWedge->mark = 1;/* check that all edges were checked */
      }
    }
  
    if(DEBUG>=2){
      #pragma omp master
      {
	for(int i=0; i < numtopnodes; i++){
	  Cnode *Vnode = topnodes[i];
	  for(int vw = 0; vw < Vnode->edgecnt; vw++){
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	    assert(VWedge->mark == 1);
	  }
	}
	for(int i=0; i < numtopnodes; i++){
	  Cnode *Vnode = topnodes[i];
	  for(int vw = 0; vw < Vnode->edgecnt; vw++){
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	    VWedge->mark = 0;
	  }
	}
      }
    }

    if(1){

      #pragma omp for schedule(dynamic,16)
      for(int i=0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	if(max(Vnode->Lcntmax,Vnode->Rcntmax) < MinCcntmax &&
	   max(Vnode->Lcntmax,Vnode->Rcntmax) < floor(min(Vnode->forward,Vnode->backward)*EdgeTrimRatio)){
          #pragma omp critical(mark)
	  {
	    if(VERB>=2 && verb){
	      if(Gmap[Vnode->mapid]->name)
		printf("N%d (%s, id=%lld) deleted as chimeric due to Ccntmax=%d,%d less than MinCcntmax=%d(edgecnt=%d)\n",
		       nodeorder[Vnode->mapid], Gmap[Vnode->mapid]->name, Gmap[Vnode->mapid]->id, Vnode->Lcntmax, Vnode->Rcntmax, MinCcntmax, Vnode->edgecnt);
	      else
		printf("N%d (id=%lld) deleted as chimeric due to Ccntmax=%d,%d less than MinCcntmax=%d(edgecnt=%d)\n",
		       nodeorder[Vnode->mapid], Gmap[Vnode->mapid]->id, Vnode->Lcntmax,Vnode->Rcntmax, MinCcntmax, Vnode->edgecnt);
	    }
	    Vnode->mark = 2;
	    /* mark all edges out of this node for deletion */
	    for(int vw = 0; vw < Vnode->edgecnt; vw++)
	      edge_delete(&edge[Vnode->edgeid[vw]],Vnode);
	  }

	  mychimcnt++;
	}

	if(MinCcnt > 0){  /* delete edges with Ccnt < MinCcnt, provided Ccnt < min(forward,backward)*EdgeTrimRatio for both nodes */
	  for(int vw = 0; vw < Vnode->edgecnt; vw++){
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	    if(VWedge->valid < 0)
	      continue;
	    Cnode *Lnode = &node[VWedge->Lmap];
	    Cnode *Rnode = &node[VWedge->Rmap];
	    if(VWedge->Ccnt < MinCcnt && 
	       VWedge->Ccnt < floor(min(Lnode->forward,Lnode->backward)*EdgeTrimRatio) &&
	       VWedge->Ccnt < floor(min(Rnode->forward,Rnode->backward)*EdgeTrimRatio)){
              #pragma omp critical(valid)
	      {
		if(VERB>=2 &&  (verb /* || nodeorder[Vnode->mapid]==462 || nodeorder[Wnode->mapid]==462 ||
				     nodeorder[Vnode->mapid] == 856 || nodeorder[Wnode->mapid]==856 ||
				     (nodeorder[Vnode->mapid] == 1528 && nodeorder[Wnode->mapid] == 134) ||
				     (nodeorder[Vnode->mapid] == 134 && nodeorder[Wnode->mapid] == 1528)*/)){
		  printf("Trimming edge=(N%d,N%d,%0.3f,%d,%d):Lcnt=%d,Rcnt=%d,Ccnt=%d,Lmap(forward=%0.1f,backward=%0.1f),Rmap(forward=%0.1f,backward=%0.1f)\n",
			 nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,
			 VWedge->Lflipped,VWedge->Rflipped,VWedge->Lcnt,VWedge->Rcnt,VWedge->Ccnt,
			 Lnode->forward,Lnode->backward,Rnode->forward,Rnode->backward);
		  fflush(stdout);
		}
		VWedge->valid = -1;
	      }

	      myedgecnt++;
	    }
	  }    
	}

	if(MaxCcnt > 0){/* delete edges with max(Lcnt,Rcnt) < MaxCcnt */
	  Cnode *Vnode = topnodes[i];
	  for(int vw = 0; vw < Vnode->edgecnt; vw++){
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	    if(VWedge->valid < 0)
	      continue;
	    if(max(VWedge->Lcnt,VWedge->Rcnt) < MaxCcnt){
              #pragma omp critical(valid)
	      {
		if(VERB>=2 &&  (verb /* || nodeorder[Vnode->mapid]==462 || nodeorder[Wnode->mapid]==462 ||
				     nodeorder[Vnode->mapid] == 856 || nodeorder[Wnode->mapid]==856 ||
				     (nodeorder[Vnode->mapid] == 1528 && nodeorder[Wnode->mapid] == 134) ||
				     (nodeorder[Vnode->mapid] == 134 && nodeorder[Wnode->mapid] == 1528)*/)){
		  printf("Trimming edge=(N%d,N%d,%0.3f,%d,%d):Lcnt=%d,Rcnt=%d,Ccnt=%d,MaxCcnt=%d\n",
			 nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,
			 VWedge->Lflipped,VWedge->Rflipped,VWedge->Lcnt,VWedge->Rcnt,VWedge->Ccnt,MaxCcnt);
		  fflush(stdout);
		}
		VWedge->valid = -1;
	      }

	      myedgecnt++;
	    }
	  }
	}
      }// omp for
    } // if(1)
    
    #pragma omp critical(edgecnt)
    {
      edgecnt += myedgecnt;
      chimcnt += mychimcnt;
    }
  } // pragma omp parallel

  unsigned int finalcnt;
  long long Ecnt = GCedges(finalcnt,edge);
  int Ncnt = GCnodes();

  if(VERB){
    printf("edgetrim(%d): removed %d unconfirmed edges and %d chimeric nodes (edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
	   verb,edgecnt,chimcnt, finalcnt, numtopnodes,wtime(),mtime());
    if(VERB>=2)
      printf("Ecnt=%lld,Ncnt=%d\n",Ecnt,Ncnt);
    if(VERB>=2)
      dumpmemmap();
    fflush(stdout);
  }

  return 0;
}

void DumpMapAlignments()
{
  /* Dump filtered alignments and maps and exit, if -dumpalign was specified */

  /* mark all alignments with trueTP = 0 */
  for(size_t i = 0; i < numaligns; i++)
    alignment[i]->trueTP = 0;

  /* mark all alignments still in Graph with trueTP = 1 */
  for(int i = 0; i < numtopnodes; i++){
    Cnode *Vnode = topnodes[i];
    for(int vw = 0; vw < Vnode->edgecnt; vw++){
      Cedge *VWedge = &edge[Vnode->edgeid[vw]];
      if(DEBUG) assert(VWedge->alignid > 0);
      alignment[VWedge->alignid - 1]->trueTP = 1;
    }
  }
    
  /* now remove all alignments with trueTP = 0 from alignment[0..numaligns-1] */
  size_t cnt = 0;
  for(size_t i = 0; i < numaligns; i++){
    Calign *p = alignment[i];
    if(!p->trueTP)
      continue;
    if(cnt < i){/* swap alignment[cnt] with alignment[i] */
      alignment[i] = alignment[cnt];
      alignment[cnt] = p;
    }
    cnt++;
  }
  if(VERB){
    printf("Reduced number of alignments from %lu to %lu\n",numaligns, cnt);
    fflush(stdout);
  }
  numaligns = cnt;

  /* output maps remaining */
  int *newmapid = new int[nummaps];/* mapping from old to new locations in Gmap[] (-1 if mapid was deleted) */

  for(int i = 0; i < nummaps; i++){
    if(DEBUG) assert(Gmap[i]);
    Gmap[i]->paired = 0;
  }
  for(int i = 0; i < numtopnodes; i++){
    int mapid = topnodes[i]->mapid;
    if(DEBUG) assert(Gmap[mapid]);
    Gmap[mapid]->paired = 1;
  }
  int j = 0;
  for(int i = 0; i < nummaps; i++){
    Cmap *pmap = Gmap[i];
    if(DEBUG) assert(pmap);
    if(!pmap->paired){
      newmapid[i] = -1;
      continue;
    }
    newmapid[i] = j;
    if(j < i){
      Cmap *tmp = Gmap[j];
      Gmap[j] = pmap;
      Gmap[i] = tmp;
    }
    j++;
  }
  if(VERB){
    printf("Reduced number of maps from %d to %d\n",nummaps,j);
    fflush(stdout);
  }
  nummaps = j;

  if(usecolor)
    colors = 2;

  output_bnx(output_prefix,Gmap,0,nummaps-1, 1, NULL, NULL, -1);

  /* renumber alignment mapids */
  for(size_t i = 0; i < numaligns; i++){
    Calign *p = alignment[i];
    if(DEBUG) assert(p!=0);
    p->mapid1 = newmapid[p->mapid1];
    p->mapid2 = newmapid[p->mapid2];
  }
    
  delete [] newmapid;

  /* output alignments remaining */
  if(DEBUG) assert(!Cfirst);
  output_align(output_prefix, 0);
    
  if(VERB>=2){
    dumpmemmap();      
    fflush(stdout);
  }
}

static Cnode *Find(Cnode *X, Cnode **parent)
{
  Cnode **p = &parent[X->nodeid];
  if(*p != X)
    *p = Find(*p,parent);
  return *p;
}

static void Union(Cnode *X, Cnode *Y, Cnode **parent, int *mark)
{
  Cnode *Xroot = Find(X,parent);
  Cnode *Yroot = Find(Y,parent);
  if(Xroot == Yroot)
    return;
  
  if(mark[Xroot->nodeid] < mark[Yroot->nodeid])
    parent[Xroot->nodeid] = Yroot;
  else if(mark[Xroot->nodeid] > mark[Yroot->nodeid])
    parent[Yroot->nodeid] = Xroot;
  else {
    parent[Yroot->nodeid] = Xroot;
    mark[Xroot->nodeid]++;
  }
}


class CNodeSet {
public:
  Cnode *node;
  Cedge *edge;
  int nodeid;/* shortcut for node->nodeid */
  int breakable;
};

/** inline function to pull out critical loop (to be vectorized) */
static inline void matrixAB(SINT *EmatrixW, SINT *nodesetMarked, int E16, int E, int *wAB)
{
  // NOTE : extend to case where SINT = char : requires E <= 127*16 (should be checked in chimtrim) : potentially 2x faster
#if USE_SSE==1
  if(sizeof(SINT)==2){/* vectorize with 16 bit signed int vectors */
    __m128i marked = _mm_load_si128((__m128i *)nodesetMarked);
    __m128i EmatrixWX = _mm_load_si128((__m128i *)EmatrixW);
    __m128i vA = _mm_andnot_si128(marked,EmatrixWX);
    __m128i vB = _mm_and_si128(marked,EmatrixWX);
    for(int x = 8; x < E16; x += 8){/* critical loop takes 19% of run time (51% exluding openMP overhead) */
      marked = _mm_load_si128((__m128i *)&nodesetMarked[x]);
      EmatrixWX = _mm_load_si128((__m128i *)&EmatrixW[x]);
      vA = _mm_add_epi16(vA, _mm_andnot_si128(marked,EmatrixWX));
      vB = _mm_add_epi16(vB, _mm_and_si128(marked,EmatrixWX));
    }

#ifdef __SSSE3__
    /* use horizontal add of vA and vB */
    vA = _mm_hadds_epi16(vA,vA);
    vA = _mm_hadds_epi16(vA,vA);
    vA = _mm_hadds_epi16(vA,vA);
    short wA = _mm_extract_epi16(vA,0);

    vB = _mm_hadds_epi16(vB,vB);
    vB = _mm_hadds_epi16(vB,vB);
    vB = _mm_hadds_epi16(vB,vB);
    short wB = _mm_extract_epi16(vB,0);

#else // No __SSSE3__

    short tmpA[8],tmpB[8];
    _mm_store_si128((__m128i*)tmpA,vA);
    _mm_store_si128((__m128i*)tmpB,vB);
    short wA = 0, wB = 0;
    for(int i = 0; i < 8; i++){
      wA += tmpA[i];
      wB += tmpB[i];
    }
#endif // No __SSSE3__

    if(DEBUG>=2){
      short newA = wA, newB = wB;
      wA = wB = 0;
      for(int x = 0; x < E; x++){
	SINT marked = nodesetMarked[x];
	SINT EmatrixWX = EmatrixW[x];
	wA += (~marked) & EmatrixWX;
	wB += marked & EmatrixWX;
      }
      assert(wA == newA);
      assert(wB == newB);
    }
    wAB[0] = wA;
    wAB[1] = wB;
  } else 
#endif // USE_SSE==1
  {/* non-intrinsic code */
    int wA = 0, wB = 0;
    for(int x = 0; x < E; x++){
      SINT marked = nodesetMarked[x];
      SINT EmatrixWX = EmatrixW[x];
      wA += (~marked) & EmatrixWX;
      wB += marked & EmatrixWX;
    }
    wAB[0] = wA;
    wAB[1] = wB;
  }
}

static inline int check_chimtrim(Cnode  *Vnode)
{
  return (0 /* nodeorder[Vnode->mapid] == 1211*//* || nodeorder[Vnode->mapid] == 247 || nodeorder[Vnode->mapid] == 245 || nodeorder[Vnode->mapid] == 243 */ /*nodeorder[Vnode->mapid] ==371 || nodeorder[Vnode->mapid]==857*/ /*1207*//*1250*/ /*246*//*243*/ /* 788 */ /*|| nodeorder[Vnode->mapid] == 1238*//* 1255*/ /*nodeorder[Vnode->mapid] == 462 ||*/ /*nodeorder[Vnode->mapid] == 958*//* nodeorder[Vnode->mapid]==243 || nodeorder[Vnode->mapid] == 856 */) ? 1 : 0;
}



int *edgecnt_old = 0;/* edgecnt[nodeid] = number of edges at start of last call to graph_chimtrim() */
int *edgecnt_new = 0;/* used to update edgecnt[] at end of call to graph_chimtrim() */
int maxpid = 0;
int **Pmark= 0;/* Pmark[0..maxpid-1] : preallocated (and zero'd) memory for each thread */
Cnode ***Pparent = 0;/* Pparent[0..maxpid-1] : preallocated (and zero'd) memory for each thread */

/** check if node is chimeric based on its connected nodes falling into disjoint sets based on their pairwise connections */
long long graph_chimtrim(int verb)
{
  /* NOTE : first call to graph_chimtrim() should compact Gmap[],node[] and edge[] arrays to include only surviving maps/nodes & edges */

  if(!PosEdgeInit){
    if(VERB>=2){
      printf("graph_chimtrim: Calling PosEdge():time=%0.6f(total=%0.6f)\n",wtime(),mtime());
      fflush(stdout);
    }
    PosEdge();
  }

  if(VERB>=2){
    printf("graph_chimtrim: Calling graph_edgecnt():time=%0.6f(total=%0.6f)\n",wtime(),mtime());
    fflush(stdout);
  }
  graph_edgecnt(0);/* update edge counts */


  int edgecnt = 0;// count of edges trimmed 
  int unchanged_cnt = 0;// number of cases where little work was needed
  int partition_cnt = 0;// number of cases where graph partitioning was required
  double partitionE_cnt = 0.00;// sum of E (node edges) when graph partitioning was required
  double partitionE4_cnt = 0.0;// sum of E^4 when graph partitioning was required
  double totalE4_cnt = 0.0;// sum of E^4 of all nodes

  if(DEBUG >=2){/* reset mark and parent */
    int mcnt = 0,pcnt=0;
    for(int i=0; i < numtopnodes; i++){
      Cnode *Wnode = topnodes[i];
      mcnt += (Wnode->mark ? 1 : 0);
      pcnt += (Wnode->parent ? 1 : 0);
      Wnode->mark = 0;
      Wnode->parent = 0;
    }
    if(mcnt > 0 || pcnt > 0){
      printf("Warning:graph_chimtrim found %d nodes with mark and %d nodes with parent not reset\n",mcnt,pcnt);
      fflush(stdout);
    }
  }

  if(CHIMTRIM_FAST){
    if(!edgecnt_old){
      edgecnt_old = new int[numnodes];
      edgecnt_new = new int[numnodes];
      for(int i = numnodes;--i>=0;)
	edgecnt_new[i] = -1;
    }
  }

  if(VERB>=2){
    printf("Starting parallel chimtrim:time=%0.6f(total=%0.6f)\n",wtime(),mtime());
    fflush(stdout);
  }

  #if USE_SSE==1
  #ifdef __SSSE3__
  if(VERB>=2){
    printf("Using SSSE3:");
    fflush(stdout);
  }
  #endif
  #endif 

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {

    int tid = 0;
    int numpid = 1;
    #ifdef _OPENMP
    tid = omp_get_thread_num ();
    numpid = omp_get_num_threads();
    if(DEBUG>=2) assert(tid < numpid);
    #endif

    if(CHIMTRIM_FAST)
    #pragma omp for schedule(static,512)
    for(int i = 0; i < numtopnodes; i++)
    {
      Cnode *Vnode = topnodes[i];
      int edgecnt = Vnode->edgecnt;
      int k = Vnode->nodeid;
      edgecnt_old[k] = edgecnt_new[k];
      edgecnt_new[k] = edgecnt;
      if(DEBUG>=2 && edgecnt_old[k] >= 0 && !(edgecnt_new[k] <= edgecnt_old[k])){
	#pragma omp critical
	{
	  printf("i=%d,k=nodeid=%d:edgecnt_old[k]=%d,edgecnt=edgecnt_new[k]=%d\n",
		 i,k,edgecnt_old[k],edgecnt);
	  fflush(stdout);
	  assert(edgecnt_new[k] <= edgecnt_old[k]);
	}
      }
    } // omp for

    // local per-thread even counters
    int myedgecnt = 0;// count of edges trimmed 
    int myunchanged_cnt = 0;// number of cases where little work was needed
    int mypartition_cnt = 0;// number of cases where graph partitioning was required
    double mypartitionE_cnt = 0;// sum of E value when graph partitioning was required
    double mypartitionE4_cnt = 0;// sum of E^4 when graph partitioning was required
    double mytotalE4_cnt = 0;// sum of E^4 of all nodes

    /* allocate heap memory local to each thread : 
       nodesetMarked[] and Ematrix[k][] are allocated as multiples of 16 bytes for vectorization */
    int boolblock = 16/sizeof(SINT);// number of SINTS to fit in 16 bytes
    int maxedges16 = (USE_SSE && (maxedges % boolblock)) ? (maxedges/boolblock + 1)*boolblock : maxedges;/* round up SINT array size to 16 bytes multiple */

    /* NOTE : allocate the following memory once (see allocation of mark[],parent[] below) */
    CNodeSet *nodeset = new CNodeSet[maxedges];
    SINT *nodesetMarked = new SINT[maxedges16]; /* nodesetMarked is shortcut for (mark[node->nodeid] ? 1 : 0) : not always updated */
    int *subset = new int[maxedges];
    SINT **Ematrix = new SINT *[maxedges];
    SINT *Ematrix_mem = new SINT[maxedges*maxedges16];
    for(int k = maxedges; --k >= 0;)
      Ematrix[k] = &Ematrix_mem[k*maxedges16];

    /* need to allocate thread-local mark[0..numnodes-1] and parent[0..numnodes-1] */
    /* NOTE: index parent[] by nodeset index= 0..E-1 instead of nodeid=0..numnodes. Unfortunately mark[] needs to be indexed by nodeid during initialization of matrix[][], but can use a SINT type to save memory  */ 
    if(TSAN || !Pmark){
      #pragma omp critical
      {
	if(!Pmark){
	  maxpid = numpid;
	  int ** LPmark = new int*[maxpid];
	  Cnode *** LPparent = new Cnode**[maxpid];
	  for(int i = maxpid; --i >= 0;)
	    LPmark[i] = 0;
	  Pparent = LPparent;
	  Pmark = LPmark;
	}
      }
    }

    if(DEBUG) assert(numpid <= maxpid);
    if(!Pmark[tid]){/* this memory will be allocated once and reused across multiple calls to graph_chimtrim() */
      Pmark[tid] = new int[numnodes];
      Pparent[tid] = new Cnode*[numnodes];
      int *mark = Pmark[tid];
      Cnode **parent = Pparent[tid];
      for(int i = numnodes; --i >= 0;){
	mark[i] = 0;
	parent[i] = 0;
      }
    }
    int *mark = Pmark[tid];
    Cnode **parent = Pparent[tid];
    if(DEBUG>=2)
      for(int i = numnodes; --i >= 0;){
	assert(mark[i] == 0);
	assert(parent[i] == 0);
      }

    #pragma omp for schedule(dynamic,16)
    for(int i=0; i < numtopnodes; i++) {
      Cnode *Vnode = topnodes[i];

      /* NOTE : if Vnode->edgecnt is large , randomly sub-sample a set of edges */

      /* assemble nodeset[0..E-1], the set of nodes connected to Vnode */
      int E = 0;
      int changed = (Vnode->edgecnt != edgecnt_old[Vnode->nodeid]) ? 1 : 0;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(DEBUG>=2) assert(((VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap) == Vnode->Enodeid[vw]);	
	Cnode *Wnode = &node[Vnode->Enodeid[vw]];
	if(CHIMTRIM_FAST && Wnode->edgecnt != edgecnt_old[Wnode->nodeid])
	  changed = 1;
	if(DEBUG>=2 && mark[Wnode->nodeid] != 0){
	  printf("Vnode=N%d(id=%lld):vw=%d,Wnode=N%d(id=%lld),Wnode->mark=%d\n",nodeorder[Vnode->mapid],Gmap[Vnode->mapid]->id,
		 vw,nodeorder[Wnode->mapid],Gmap[Wnode->mapid]->id,Wnode->mark);
	  for(int vw = 0; vw < Vnode->edgecnt; vw++){
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	    Cnode *Znode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
	    printf("  vw=%d/%d:Wnode=N%d[id=%lld],mark[Wnode->nodeid]=%d\n",vw,Vnode->edgecnt,nodeorder[Znode->mapid],Gmap[Znode->mapid]->id,mark[Znode->nodeid]);
	  }
	  fflush(stdout);
	  assert(mark[Wnode->nodeid] == 0);
	}
	nodeset[E].node = Wnode;
	nodeset[E].nodeid = Wnode->nodeid;
	nodeset[E].edge = VWedge;

	/* If connection to Wnode cannot be broken because min(Wnode->Lcntmax,Wnode->Rcntmax) < MinChimCnt) skip this node.
	   This additional condition is to avoid deleting edges in shallow regions near fragile nicking sites. */
	nodeset[E].breakable = (min(Wnode->Lcntmax,Wnode->Rcntmax) >= MinChimCcnt) ? 1 : 0;

	parent[Wnode->nodeid] = Wnode;
	mark[Wnode->nodeid] = ++E;
      }
      if(DEBUG>=2)assert(E == Vnode->edgecnt);
      double E1 = E;
      mytotalE4_cnt += E1*E1*E1*E1;

      if(E <= 0){
	if(VERB>=2 && verb){
	  printf("Vnode=N%d:edgecnt=%d,E=%d (No flexible edges)\n",
		 nodeorder[Vnode->mapid],Vnode->edgecnt,E);
	  fflush(stdout);
	}
	continue;
      }
      if(CHIMTRIM_FAST && !changed){// Why is the effect of CHIMTRIM_FAST on run time so small even when unchanged_cnt is 90% of total nodes ?
	myunchanged_cnt++;

	/* reset mark, parent */
	for(int vw = E; --vw >= 0;){
	  if(DEBUG>=2) assert(nodeset[vw].node->nodeid == nodeset[vw].nodeid);
	  int nodeid = nodeset[vw].nodeid;
	  mark[nodeid] = 0;
	  parent[nodeid] = 0;
	}
	continue;
      }

      /* assemble the edge matrix Ematrix[i][j]==1 IFF there is an edge between nodeset[i] and nodeset[j] */
      for(int w = E; --w >= 0;)
	for(int x = E; --x >= 0;)
	  Ematrix[w][x] = 0;

      for(int w = 0; w < E; w++){
	Cnode *Wnode = nodeset[w].node;
	for(int wx = Wnode->edgecnt; --wx >= 0;){
          if(DEBUG>=2){
	    Cedge *WXedge = &edge[Wnode->edgeid[wx]];
	    assert(((WXedge->Lmap == Wnode->nodeid) ? WXedge->Rmap : WXedge->Lmap) == Wnode->Enodeid[wx]);
	  }
	  int Xmark = mark[Wnode->Enodeid[wx]];
	  if(!Xmark)
	    continue;// NOTE : avoid branch stall by allowing -1 index in Ematrix[][] (but need to maintain 16-byte alignment of 2nd index 0, does not require extra space unless E==E16, if padded elements are reset to 0 later)
	  int x = Xmark - 1;
	  Ematrix[w][x] = Ematrix[x][w] = 1;
	}
      }

      /* NOTE : from here on mark[] index 0..numnodes could be replaced by the nodeset index 0..E-1 (or use seperate Nmark[]) */

      /* first do a fast disjoint set check using Union-Find operations */

      /* now use edges between nodes in nodeset[0..E-1] to do a union of the nodes into sub-sets. 
	 Use (parent and mark) as (parent and rank+1) for the Union-Find operations */
      for(int w = 0; w < E; w++){
	Cnode *Wnode = nodeset[w].node;
	for(int x = 0; x < E; x++){
	  if(!Ematrix[w][x])
	    continue;
	  Cnode *Xnode = nodeset[x].node;
	  /* perform a union over Wnode and Xnode */
	  Union(Wnode,Xnode,parent,mark);
	}
      }

      /* reset mark[] */
      for(int vw = E; --vw >= 0;){
	if(DEBUG>=2) assert(nodeset[vw].node->nodeid == nodeset[vw].nodeid);
	mark[nodeset[vw].nodeid] = 0;
      }

      /* now check if we ended up with more than one sub-set and compute the average Ccnt value of each subset */

      int rootcnt = 0;
      for(int vw = E; --vw >= 0;){
	Cnode *Wnode = nodeset[vw].node;
	Cnode *Wroot = Find(Wnode,parent);
	if(Wnode == Wroot)
	  rootcnt++;
	mark[Wroot->nodeid]++;
      }
      int maxcnt = 0;
      Cnode *Troot = 0;
      for(int vw = E; --vw >= 0;){
	Cnode *Wroot = Find(nodeset[vw].node,parent);
	if(mark[Wroot->nodeid] > maxcnt){
	  maxcnt = mark[Wroot->nodeid];
	  Troot = Wroot;
	}
      }
      int origedgecnt = myedgecnt;

      if(rootcnt > 1 || (VERB>=2 && verb && check_chimtrim(Vnode))){/* found disjoint sets */
	if(rootcnt > 1){/* delete edges that don't belong to largest disjoint set (unless min(Wnode->Lcntmax,Wnode->Rcntmax) < MinChimCcnt, ie breakable==0). */
	  int Tmark = mark[Troot->nodeid];
	  for(int vw = E; --vw >= 0;){
	    Cnode *Wnode = nodeset[vw].node;
	    Cnode *Wroot = Find(Wnode,parent);
	    Cedge *VWedge = nodeset[vw].edge;
	    if(mark[Wroot->nodeid] <  Tmark && nodeset[vw].breakable){
	      int valid;
	      if(TSAN || (valid = VWedge->valid) != -1){
		if(!CHIMTRIM_AVOID){
                  #pragma omp critical (valid) // HERE HERE HERE : try atomic ?
		  {
		    valid = VWedge->valid;
		    VWedge->valid = -1;
                  }
                } else { // CHIMTRIM_AVOID
                  #pragma omp critical (valid)
		  {
		    valid = VWedge->valid;
		    if(valid != -1){
		      VWedge->valid = -1;

		      /* reduce Lcntmax,Rcntmax at Vnode,Wnode to avoid deleting multiple edges near fragile nicking sites before re-computing edge counts */
		      if(vw < Vnode->PosEdge){/* negative edge : from Wnode to Vnode */
		        if(Wnode->Rcntmax > 0)
			  Wnode->Rcntmax--;
			if(Vnode->Lcntmax > 0)
			  Vnode->Lcntmax--;
	              } else {/* positive edge : from Vnode to Wnode */
		        if(Vnode->Rcntmax > 0)
			  Vnode->Rcntmax--;
			if(Wnode->Lcntmax > 0)
			  Wnode->Lcntmax--;
	              }
		    }
	          }
		} // CHIMTRIM_AVOID
              } // if(TSAN || (valid = VWedge->valid) != -1)

	      if(/* WAS10 !CHIMTRIM_AVOID || */ valid != -1)
	        myedgecnt++;
	    }
	  }
	}

	if(VERB>=2 && verb) {
	  printf("N%d (Lcntmax=%d,Rcntmax=%d): %d connected nodes form %d disjoint sets, largest set has %d nodes",
		 nodeorder[Vnode->mapid], Vnode->Lcntmax,Vnode->Rcntmax, E, rootcnt, maxcnt);
	  printf(" (edgecnt=%d->%d)\n",origedgecnt,myedgecnt);
	  for(int i=0; i < E; i++){
	    Cnode *Tnode = nodeset[i].node;
	    Cnode *Troot = Find(Tnode,parent);
	    Cedge *VWedge = nodeset[i].edge;
	    printf(" edge%d to N%d: root=N%d, cnt=%d (N%d:Lcntmax=%d,Rcntmax=%d; N%d:Lcntmax=%d,Rcntmax=%d) valid=%d\n",
		   i, nodeorder[Tnode->mapid], nodeorder[Troot->mapid], mark[Troot->nodeid],
		   nodeorder[Vnode->mapid], Vnode->Lcntmax, Vnode->Rcntmax, nodeorder[Tnode->mapid], Tnode->Lcntmax, Tnode->Rcntmax,
		   VWedge->valid);
	  }
	  fflush(stdout);
	}
      }

      /* reset mark, parent  and initialize nodeset[].marked */
      for(int vw = E; --vw >= 0;){
	if(DEBUG>=2) assert(nodeset[vw].node->nodeid == nodeset[vw].nodeid);
	int nodeid = nodeset[vw].nodeid;
	mark[nodeid] = 0;
	parent[nodeid] = 0;
      }

      if(CHIMTRIM && myedgecnt <= origedgecnt){/* try a greedy two way graph-partitioning approach for "almost" disjoint sets */
	mypartition_cnt++;
	double dE = E;
	mypartitionE_cnt += dE;
	mypartitionE4_cnt += dE*dE*dE*dE;

	int subsetcnt = 0;/* smaller node subset subset[0..subsetcnt-1] */
	int smallcnt = 0;/* number of internal edges within smaller node subset */
	int largecnt = 0;/* number of internal edges within larger node subset (all nodes in nodeset[], except those with mark=1 that are in subset[]) */
	int crosscnt = 0;/* number of edges between two subsets */

	int E16 = (USE_SSE==1 && (E % boolblock)) ? (E/boolblock + 1)*boolblock : E;/* round up E to 16 bytes multiple */

	#if USE_SSE==1	/* pad nodesetMarked, Ematrix[w] with 0 up to E16 */
	/* initialize nodeset[].marked */
	for(int vw = E16; --vw >= 0;)
	  nodesetMarked[vw] = 0;
	for(int w = E; --w >= 0;)
	  for(int x = E16; --x >= E;)
	    Ematrix[w][x] = 0;
	#else
	/* initialize nodeset[].marked */
	for(int vw = E; --vw >= 0;)
	  nodesetMarked[vw] = 0;
	#endif

	/* initialize largecnt */
	for(int w = 0; w < E; w++)
	  for(int x = 0; x < w; x++)
	    largecnt += Ematrix[w][x];

	int bestsubsetcnt = 0;/* number of nodes in subset[] for best partitioning */
	double bestratio = 1.0;/* ratio of crosscnt / (subsetcnt * (E-subsetcnt))  for best partitioning */
	int bestsmallcnt = 0,bestlargecnt = largecnt,bestcrosscnt=0;/* number of internal/cross edges in the best two subsets */

	int breakable = 0 /* Vmax */;/* smaller subset number of breakable edges with min(Lcntmax,Rcntmax) >= MinChimCcnt */
	int bestbreakable = breakable;

	/* keep adding one node to smaller subset: pick the node that minimizes crosscnt */
	while(subsetcnt < E){
	  /* find node w in larger set with smallest value if wA - wB, where
	     wA = number of connections between node w and larger set (A)
	     wB = number of connections between node w and smaller set (B) */
	  int bestnode = -1, bestA = E, bestB = 0;
	  double bratio = 1.0;
	  for(int w = 0; w < E; w++){
	    if(DEBUG>=2)assert((mark[nodeset[w].node->nodeid]?1:0) == nodesetMarked[w]);
	    if(nodesetMarked[w])
	      continue;

	    int wAB[2];
	    matrixAB(Ematrix[w],nodesetMarked,E16,E,wAB);

	    if(smallcnt > 0 /* && (nodeorder[Vnode->mapid]==371 || nodeorder[Vnode->mapid]==857)*/){
	      int nsubsetcnt = subsetcnt+1;
	      int nlargecnt = largecnt - wAB[0];
	      int nsmallcnt = smallcnt + wAB[1];
	      int ncrosscnt = crosscnt + wAB[0] - wAB[1];
	      int nbreakable = breakable + nodeset[w].breakable;
	      double ratio = (nsubsetcnt >= E-1 || nsmallcnt <= 0 || nlargecnt <= 0) ? 1.0 : (0.5 * ncrosscnt)/sqrt((double)(nsmallcnt*nlargecnt));	    
	      if(bestnode < 0 || (nbreakable && ratio <= bratio && (ratio < bratio - 1e-8 || wAB[0]-wAB[1] < bestA-bestB))){
		bestnode = w;
		bratio = (nbreakable ? ratio : 1.0);
		bestA = wAB[0];
		bestB = wAB[1];
	      }
	    } else if(wAB[0] - wAB[1] < bestA - bestB){
	      bestnode = w;
	      bestA = wAB[0];
	      bestB = wAB[1];
	    }
	  }
	  /* move node bestnode to subset */
	  if(DEBUG>=2) assert(bestnode >= 0);
	  breakable += nodeset[bestnode].breakable;
	  //	  Cnode *Bnode = nodeset[bestnode].node;
	  subset[subsetcnt++] = bestnode;
	  if(DEBUG>=2) assert(nodeset[bestnode].node->nodeid == nodeset[bestnode].nodeid);
	  if(DEBUG>=2) assert(subsetcnt != 0);
	  mark[nodeset[bestnode].nodeid] = subsetcnt;
	  nodesetMarked[bestnode] = 1;
	  largecnt -= bestA;
	  smallcnt += bestB;
	  crosscnt += bestA - bestB;
	  double ratio = (subsetcnt<=1 || subsetcnt >= E-1 || smallcnt <= 0 || largecnt <= 0) ? 1.0 : (0.5 * crosscnt)/sqrt((double)(smallcnt*largecnt));
	
	  if(ratio < bestratio - 1e-8 && breakable){
	    bestratio = ratio;
	    bestsubsetcnt = subsetcnt;
	    bestsmallcnt = smallcnt;
	    bestlargecnt = largecnt;
	    bestcrosscnt = crosscnt;
	    bestbreakable = breakable;
	  }
	  if(VERB>=2 && verb && (verb >=2 || check_chimtrim(Vnode))){
	    printf("Vnode=N%d:E=%d,subsetcnt=%d:bestnode=N%d,smallcnt=%d,largecnt=%d,crosscnt=%d,ratio=%0.6f(bestratio=%0.6f,bestcnt=%d, smallcnt=%d,largecnt=%d,crosscnt=%d),b=%d\n",
		   nodeorder[Vnode->mapid],E,subsetcnt,nodeorder[nodeset[bestnode].node->mapid],smallcnt,largecnt,crosscnt,ratio,bestratio,bestsubsetcnt,bestsmallcnt,bestlargecnt,bestcrosscnt,breakable);
	    fflush(stdout);
	  }
	}

	/* remove sub optimal nodes from subset[] */
	for(int w = bestsubsetcnt;w < subsetcnt;w++){
	  int s = subset[w];
	  if(DEBUG>=2) assert(nodeset[s].node->nodeid == nodeset[s].nodeid);
	  mark[nodeset[s].nodeid] = 0 ;
	  nodesetMarked[s] = 0;
	}
	subsetcnt = bestsubsetcnt;
	smallcnt = bestsmallcnt;
	largecnt = bestlargecnt;
	crosscnt = bestcrosscnt;
	breakable = bestbreakable;

	/* try moving nodes between the two sets until no further progress is possible */
	int progress = /* (nodeorder[Vnode->mapid]==371 || nodeorder[Vnode->mapid]==857) ? 1 : 0*/ 1;
	int cnt = 0;
	for(; cnt < 1000 && progress; cnt++){
	  progress = 0;

	  /* try deleting a node from marked subset */
	  double bratio = bestratio;
	  double origratio = bratio;
	  int bestnode = -1;
	  int bbreakable = breakable;
	  int bestA = 0, bestB = 0;
	  for(int w = 0; w < E; w++){
	    if(DEBUG>=2) assert((mark[nodeset[w].node->nodeid]?1:0) == nodesetMarked[w]);
	    if(!nodesetMarked[w])
	      continue;

	    int wAB[2];
	    matrixAB(Ematrix[w],nodesetMarked,E16,E,wAB);

	    int nsubsetcnt = subsetcnt-1;
	    int nlargecnt = largecnt + wAB[0];
	    int nsmallcnt = smallcnt - wAB[1];
	    int ncrosscnt = crosscnt - wAB[0] + wAB[1];
	    int nbreakable = breakable - nodeset[w].breakable;
	    double ratio = (nsubsetcnt<=1 || nsmallcnt <= 0 || nlargecnt <= 0) ? 1.0 : (0.5 * ncrosscnt)/sqrt((double)(nsmallcnt*nlargecnt));
	    if(ratio < bratio - 1e-8 && nbreakable){
	      bestnode = w;
	      bratio = ratio;
	      bestA = wAB[0];
	      bestB = wAB[1];
	      bbreakable = nbreakable;
	    }
	  }
	  if(bestnode != -1){
	    if(DEBUG>=2) assert(bratio < bestratio - 1e-9);
	    if(DEBUG>=2) assert(bbreakable);
	    progress++;
	    subsetcnt--;
	    if(DEBUG>=2) assert(nodeset[bestnode].node->nodeid == nodeset[bestnode].nodeid);
	    nodesetMarked[bestnode] = mark[nodeset[bestnode].nodeid] = 0;
	    largecnt += bestA;
	    smallcnt -= bestB;
	    crosscnt -= bestA-bestB;
	    bestratio = bratio;
	    breakable = bbreakable;
	    if(DEBUG>=2)assert(bratio < origratio - 1e-9);
	    if(VERB>=2 && verb && (verb >=2 || check_chimtrim(Vnode))){
	      printf("Vnode=N%d:E=%d,subsetcnt=%d:deleted N%d,smallcnt=%d,largecnt=%d,crosscnt=%d,ratio=%0.8f->%0.8f,b=%d\n",
		     nodeorder[Vnode->mapid],E,subsetcnt,nodeorder[nodeset[bestnode].node->mapid],smallcnt,largecnt,crosscnt,origratio,bestratio,breakable);
	      fflush(stdout);
	    }
	  }

	  /* try adding a node to subset[] */
	  bestnode = -1;
	  origratio = bratio;
	  for(int w = 0; w < E; w++){
	    if(DEBUG>=2) assert((mark[nodeset[w].node->nodeid]?1:0) == nodesetMarked[w]);
	    if(nodesetMarked[w])
	      continue;

	    int wAB[2];
	    matrixAB(Ematrix[w],nodesetMarked,E16,E,wAB);

	    int nsubsetcnt = subsetcnt+1;
	    int nlargecnt = largecnt - wAB[0];
	    int nsmallcnt = smallcnt + wAB[1];
	    int ncrosscnt = crosscnt + wAB[0] - wAB[1];
	    int nbreakable = breakable + nodeset[w].breakable;
	    double ratio = (nsubsetcnt >= E-1 || nsmallcnt <= 0 || nlargecnt <= 0) ? 1.0 : (0.5 * ncrosscnt)/sqrt((double)(nsmallcnt*nlargecnt));
	    if(ratio < bratio - 1e-8 && nbreakable){
	      bestnode = w;
	      bratio = ratio;
	      bestA = wAB[0];
	      bestB = wAB[1];
	      bbreakable = nbreakable;
	    }
	  }

	  if(bestnode != -1){
	    if(DEBUG>=2) assert(bratio < bestratio - 1e-9);
	    if(DEBUG>=2) assert(bbreakable);
	    progress++;
	    subsetcnt++;
	    if(DEBUG>=2) assert(nodeset[bestnode].node->nodeid == nodeset[bestnode].nodeid);
	    if(DEBUG>=2) assert(subsetcnt != 1);
	    mark[nodeset[bestnode].nodeid] = subsetcnt;
	    nodesetMarked[bestnode] = 1;
	    largecnt -= bestA;
	    smallcnt += bestB;
	    crosscnt += bestA-bestB;
	    bestratio = bratio;
	    breakable = bbreakable;
	    if(DEBUG>=2)assert(bratio < origratio - 1e-9);
	    if(VERB>=2 && verb && (verb >=2 || check_chimtrim(Vnode))){
	      printf("Vnode=N%d:E=%d,subsetcnt=%d:added N%d,smallcnt=%d,largecnt=%d,crosscnt=%d,ratio=%0.8f->%0.8f,b=%d\n",
		     nodeorder[Vnode->mapid],E,subsetcnt,nodeorder[nodeset[bestnode].node->mapid],smallcnt,largecnt,crosscnt,origratio,bestratio,breakable);
	      fflush(stdout);
	    }
	  }
	}
	if(cnt >= 1000){
	  printf("Loop time out: cnt=%d, Vnode=N%d\n",cnt,nodeorder[Vnode->mapid]);
	  exit(1);
	}

	if(bestratio < ChimTrimRatio){
	  /* delete edges to smaller subset, provided breakable==1 */
	  if(VERB>=2 && verb){
	    printf("Vnode=N%d(C=%d,%d):E=%d,bestratio=%0.6f(min=%0.3f),bestcnt=%d,smallcnt=%d,largecnt=%d,crosscnt=%d):\n",
		   nodeorder[Vnode->mapid],Vnode->Lcntmax,Vnode->Rcntmax,E,bestratio,ChimTrimRatio,subsetcnt,smallcnt,largecnt,crosscnt);
	    printf("%s", smallcnt <= largecnt ? "small subset:" : "large subset:");
	    for(int w = 0; w < E; w++){
	      Cnode *Wnode  = nodeset[w].node;
	      if(mark[Wnode->nodeid])
		printf("  N%d(C=%d,%d)", nodeorder[Wnode->mapid], Wnode->Lcntmax,Wnode->Rcntmax);
	    }
	    printf("%s", smallcnt <= largecnt ? "\nlarge subset:" : "\nsmall subset:");
	    for(int w = 0; w < E; w++){
	      Cnode *Wnode  = nodeset[w].node;
	      if(!mark[Wnode->nodeid])
		printf(" N%d(C=%d,%d)", nodeorder[Wnode->mapid], Wnode->Lcntmax, Wnode->Rcntmax );
	    }
	    printf("\n");
	    fflush(stdout);
	  }

	  for(int vw = 0; vw < E; vw++){
	    Cnode *Wnode = nodeset[vw].node;
	    if(DEBUG>=2) assert((mark[Wnode->nodeid] ? 1 : 0) == nodesetMarked[vw]);
	    if((smallcnt <= largecnt) ? !nodesetMarked[vw] : nodesetMarked[vw])
	      continue;
	    if(!nodeset[vw].breakable)
	      continue;
	    Cedge *VWedge = nodeset[vw].edge;

#if 1
	      int valid;
	      if(TSAN || (valid = VWedge->valid) != -1){
		if(!CHIMTRIM_AVOID){
                  #pragma omp critical (valid) // HERE HERE : try atomic ?
		  {
		    valid = VWedge->valid;
		    VWedge->valid = -1;
                  }
                } else { // CHIMTRIM_AVOID
                  #pragma omp critical (valid)
		  {
		    valid = VWedge->valid;
		    if(valid != -1){
		      VWedge->valid = -1;

		      /* reduce Lcntmax,Rcntmax at Vnode,Wnode to avoid deleting multiple edges near fragile nicking sites before re-computing edge counts */
		      if(vw < Vnode->PosEdge){/* negative edge : from Wnode to Vnode */
		        if(Wnode->Rcntmax > 0)
			  Wnode->Rcntmax--;
			if(Vnode->Lcntmax > 0)
			  Vnode->Lcntmax--;
	              } else {/* positive edge : from Vnode to Wnode */
		        if(Vnode->Rcntmax > 0)
			  Vnode->Rcntmax--;
			if(Wnode->Lcntmax > 0)
			  Wnode->Lcntmax--;
	              }
		    }
	          }
		} // CHIMTRIM_AVOID
              } // if(TSAN || (valid = VWedge->valid) != -1)

	      if(/* WAS10 !CHIMTRIM_AVOID || */ valid != -1)
	        myedgecnt++;
#else // OLD CODE triggered TSAN
	    int valid = VWedge->valid;
	    if(valid != -1){
	      if(!CHIMTRIM_AVOID){
		VWedge->valid = -1;
	      } else {
                #pragma omp critical (valid)
		{
		  valid = VWedge->valid;
		  if(valid != -1){
		    VWedge->valid = -1;

		    if(CHIMTRIM_AVOID){
		      /* reduce Lcntmax,Rcntmax at Vnode,Wnode to avoid deleting multiple edges near fragile nicking sites before re-computing edge counts */
		      if(vw < Vnode->PosEdge){/* negative edge : from Wnode to Vnode */
			if(Wnode->Rcntmax > 0)
			  Wnode->Rcntmax--;
			if(Vnode->Lcntmax > 0)
			  Vnode->Lcntmax--;
		      } else {/* positive edge : from Vnode to Wnode */
			if(Vnode->Rcntmax > 0)
			  Vnode->Rcntmax--;
			if(Wnode->Lcntmax > 0)
			  Wnode->Lcntmax--;
		      }
		    }
		  }
		}
	      }

	      if(!CHIMTRIM_AVOID || valid != -1)
		myedgecnt++;
	    }
#endif // OLD CODE
	  }
	}

	/* reset mark */
	for(int w = E; --w >= 0;){
	  if(DEBUG>=2) assert(nodeset[w].node->nodeid == nodeset[w].nodeid);
	  mark[nodeset[w].nodeid] = 0;
	}
      }
    }

    /* free up heap memory local to each thread */
    delete [] nodeset;
    delete [] nodesetMarked;
    delete [] subset;
    delete [] Ematrix_mem;
    delete [] Ematrix;

    /* update event counters */
    #pragma omp critical (counters)
    {
       edgecnt += myedgecnt;
       unchanged_cnt += myunchanged_cnt;
       partition_cnt += mypartition_cnt;
       partitionE_cnt += mypartitionE_cnt;
       partitionE4_cnt += mypartitionE4_cnt;
       totalE4_cnt += mytotalE4_cnt;
    }
  } // pragma omp parallel

  if(DEBUG>=2){
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      assert(Vnode->mark == 0);
      assert(Vnode->parent == 0);
    }
  }

  if(DEBUG>=2 && CHIMTRIM_FAST){
    #pragma omp parallel for num_threads(numthreads) schedule(static,32) if(numthreads > 1)
    for(int i = 0; i < numtopnodes; i++)
    {
      Cnode *Vnode = topnodes[i];
      int edgecnt = Vnode->edgecnt;
      int k = Vnode->nodeid;
      if(DEBUG && !(edgecnt == edgecnt_new[k])){
	#pragma omp critical
	{
	  printf("i=%d,k=nodeid=%d:edgecnt_new[k]=%d,edgecnt=%d\n",
		 i,k,edgecnt_new[k],edgecnt);
	  fflush(stdout);
	  assert(edgecnt == edgecnt_new[k]);
	}
      }
    }
  }

  unsigned int finalcnt;
  long long Ecnt = GCedges(finalcnt,edge);
  int Ncnt = GCnodes();

  if(VERB /*&& edgecnt > 0 */){
    printf("graph_chimtrim(%d): deleted %d chimeric edges,unchanged=%d,partition=%d(Eav=%0.1f,E4sum=%8.6e),totalE4sum=%8.6e,(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
	   verb,edgecnt,unchanged_cnt,partition_cnt,partitionE_cnt/partition_cnt,partitionE4_cnt,totalE4_cnt,finalcnt,numtopnodes,wtime(),mtime());
    fflush(stdout);
  }
  if(DEBUG>=2 && CHIMTRIM_FAST){
    #pragma omp parallel for num_threads(numthreads) schedule(static,32) if(numthreads > 1)
    for(int i = 0; i < numtopnodes; i++)
    {
      Cnode *Vnode = topnodes[i];
      int edgecnt = Vnode->edgecnt;
      int k = Vnode->nodeid;
      if(DEBUG && !(edgecnt <= edgecnt_new[k])){
	#pragma omp critical
	{
	  printf("i=%d,k=nodeid=%d:edgecnt_new[k]=%d,edgecnt=%d\n",
		 i,k,edgecnt_new[k],edgecnt);
	  fflush(stdout);
	  assert(edgecnt <= edgecnt_new[k]);
	}
      }
    }
  }

  return Ecnt+Ncnt;
}

static int cdepth = 0;/* depth of recursion of cov_update */
#if DEBUG>=2
static int depthwarning = 0;
static int *edgeid_stack = 0;
#endif
static Cedge *cov_root = 0;/* used for debugging : root value of pedge */
static int redudancy_iter = -1;/* used for debugging : graph_redudancy value of iter */

/* recursively update "cov" value of an Edge */
/* CROSSBREAK==0 : free memory of update lists
   CROSSBREAK==1 : copy edge->update list to edge->updated */
static void cov_update(Cedge *pedge)
{
#if DEBUG>=2
  if(!edgeid_stack)
    edgeid_stack = new int[100000];
  edgeid_stack[cdepth] = pedge->edgeid;
  if(++cdepth >= 100000 && !depthwarning){
    depthwarning = 1;
    printf("Warning: cov_update recursion depth=%d : probably a cycle in DAG, which will cause eventual stack overlow\n",cdepth);
    printf("pedge=(%d,%d,%0.3f,%d,%d,mark=%d,id=%d), root=(%d,%d,%0.3f,%d,%d,mark=%d,id=%d)\n",
	   pedge->Lmap,pedge->Rmap,pedge->distance,pedge->Lflipped,pedge->Rflipped,pedge->updatemark,pedge->edgeid,
	   cov_root->Lmap,cov_root->Rmap,cov_root->distance,cov_root->Lflipped,cov_root->Rflipped,cov_root->updatemark,cov_root->edgeid);
    fflush(stdout);
    exit(1);
  }
#endif // DEBUG>=2

#if DEBUG >= 2
  size_t origsize,origsized;
  if(CROSSBREAK){
    origsize = pedge->update.size();
    origsized = pedge->updated.size();
  }
#endif
#if DAG_CHECK
  if(pedge->updatemark){
    printf("cov_update:updatemark=%d:depth=%d pedge already on call stack.\nCycle in update graph!\n", pedge->updatemark,cdepth);
    printf("pedge=(%d,%d,%0.3f,%d,%d,mark=%d,id=%d), root=(%d,%d,%0.3f,%d,%d,mark=%d,id=%d)\n",
	   pedge->Lmap,pedge->Rmap,pedge->distance,pedge->Lflipped,pedge->Rflipped,pedge->updatemark,pedge->edgeid,
	   cov_root->Lmap,cov_root->Rmap,cov_root->distance,cov_root->Lflipped,cov_root->Rflipped,cov_root->updatemark,cov_root->edgeid);
    dumpmemmap();
    for(int d = 0; d < cdepth; d++){
      Cedge *tedge = &edge[edgeid_stack[d]];
      printf("d=%2d,edgeid=%d:edge=(%d,%d,%0.3f,%d,%d,mark=%d,id=%d)\n",
	     d,edgeid_stack[d],tedge->Lmap,tedge->Rmap,tedge->distance,tedge->Lflipped,tedge->Rflipped,tedge->updatemark,tedge->edgeid);
    }
    fflush(stdout);
    exit(1);
  }
  pedge->updatemark = max(1,cdepth);/* to catch cycles in the tree */
#endif

  std::vector<int> origupdate;
  origupdate.swap(pedge->update);

  for(std::vector<int>::reverse_iterator it = origupdate.rbegin() ; it != origupdate.rend(); ++it){ int edgeid = *it;
    Cedge *uedge = &edge[edgeid];
#if DEBUG >=2 
    if(cov_root->Lmap == 0 && cov_root->Rmap==695732){
      printf("%2d:",cdepth);
      for(int t = 0; t < cdepth; t++)
	printf(" ");
      printf("uedge=(%7d,%7d,%7.3f,%7d,%7d,mark=%2d,id=%9d), pedge=(%7d,%7d,%7.3f,%7d,%7d,mark=%2d,id=%9d), root=(%7d,%7d,%7.3f,%7d,%7d,mark=%2d,id=%9d)\n",
	     uedge->Lmap,uedge->Rmap,uedge->distance,uedge->Lflipped,uedge->Rflipped,uedge->updatemark,uedge->edgeid,
	     pedge->Lmap,pedge->Rmap,pedge->distance,pedge->Lflipped,pedge->Rflipped,pedge->updatemark,pedge->edgeid,
	     cov_root->Lmap,cov_root->Rmap,cov_root->distance,cov_root->Lflipped,cov_root->Rflipped,cov_root->updatemark,cov_root->edgeid);
      fflush(stdout);
    }
#endif

    if(DEBUG>=2) assert(pedge->Lmap == uedge->Lmap ||
					pedge->Lmap == uedge->Rmap ||
					pedge->Rmap == uedge->Lmap ||
					pedge->Rmap == uedge->Rmap);
#if DEBUG>=2
    if(CFAST >= 2 && redudancy_iter < CFAST-1){
      if(!(uedge->distance > pedge->distance)){
	printf("pedge=(%d,%d,%0.3f,%d,%d,mark=%d,id=%d), root=(%d,%d,%0.3f,%d,%d,mark=%d,id=%d)\n",
	       pedge->Lmap,pedge->Rmap,pedge->distance,pedge->Lflipped,pedge->Rflipped,pedge->updatemark,pedge->edgeid,
	       cov_root->Lmap,cov_root->Rmap,cov_root->distance,cov_root->Lflipped,cov_root->Rflipped,cov_root->updatemark,cov_root->edgeid);
	printf("uedge=(%d,%d,%0.3f,%d,%d,mark=%d,id=%d): it= %ld/%ld\n", 
	       uedge->Lmap,pedge->Rmap,pedge->distance,pedge->Lflipped,pedge->Rflipped,pedge->updatemark,pedge->edgeid,
	       it - origupdate.rbegin(), origupdate.rend()-origupdate.rbegin());
	dumpmemmap();
	fflush(stdout);
	assert(uedge->distance > pedge->distance);
      }
    }
#endif
    double Iucnt = 1.0/uedge->ucnt;
    cov_update(uedge);
    pedge->cov += uedge->cov * Iucnt;

    if(CROSSBREAK) /* append edgeid to pedge->updated */
      pedge->updated.push_back(edgeid);
  }

#if DEBUG>=2
  if(CROSSBREAK){
    if(pedge->update.size() != 0 || pedge->updated.size() != origsize + origsized){
      printf("cov_update:pedge:id=%d,origsize=%lld,origsized=%lld,size=%lld,sized=%lld\n",
	     pedge->edgeid,(long long)origsize,(long long)origsized,
	     (long long)pedge->update.size(),(long long)pedge->updated.size());
      fflush(stdout);

      assert(pedge->update.size() == 0);
      assert(pedge->updated.size() == origsize + origsized);
    }
  }
#endif

  if(DEBUG>=2) cdepth--;

#if DAG_CHECK
  pedge->updatemark = 0;
#endif
}

#if 0
static int DFSpath(int parent, int dest)
{
  Cnode *Pnode = &node[parent];
  if(Pnode->mark)
    return 0;
  Pnode->mark = 1;
  for(int k = 0; k < Pnode->edgecnt; k++){
    Cedge *pedge = &edge[Pnode->edgeid[k]];
    if(pedge->valid <=0)
      continue;
    int child = (pedge->Lmap == parent) ? pedge->Rmap : pedge->Lmap;
    if(child == dest)
      return 1;
    if(DFSpath(child,dest))
      return 1;
  }
  return 0;
}
#endif

#if 0
/* return 1 if there is path from nodeid=source to nodeid=dest */
static int findpath(int source, int dest)
{
  if(source == dest)
    return 1;

  int result = DFSpath(source,dest);

  Cnode *Wnode = node;
  for(int k=0; k < numnodes; k++, Wnode++)
    Wnode->mark = 0;
  
  return result;
}
#endif

static int DFSverb = 0;

/** If inc = 1, mark will be incremented from 0 to 1 and 2.
   If inc = -1, mark will be decremented from 2 to 1 and 0. 
   return 0 if no cycle OR path to edge with mark == -1 was found, 
   return 1 if path to edge with mark == -1 was found
   return -1 if cycle was found
 */
static int updateDFS(Cedge *source, int inc, int depth, char *Emark)
{
  int edgeid = source->edgeid;
  int mark = Emark[edgeid];

  if(DEBUG>=2 && mark == 1){
    printf("updateDFS:source=(N%d,N%d,%0.3f,%d,%d):mark=%d(already a cycle in update graph)\n",
	   nodeorder[node[source->Lmap].mapid],nodeorder[node[source->Rmap].mapid],source->distance,
	   source->Lflipped,source->Rflipped,mark);
    fflush(stdout);
    assert(mark != 1);
  }

  if(DEBUG>=2) assert(mark >= -1);
  if(VERB>=2){
    printf("updateDFS:source=x%p(N%d,N%d,%0.3f,%d,%d),inc=%d:mark=%d,depth=%d\n",
	   (void *)source, nodeorder[node[source->Lmap].mapid],nodeorder[node[source->Rmap].mapid],source->distance,
	   source->Lflipped,source->Rflipped,inc,mark,depth);
    fflush(stdout);
  }
  if(mark == -1)
    return 1;

  if(mark == inc+1/* (inc > 0 ? 2 : 0)*/)/* this node already explored and finished */
    return 0;

  int result = 0;

  Emark[edgeid] += inc;/* mark source as active */
  if(DEBUG>=2) assert(Emark[edgeid] == 1);

  /* recurse over all children */
  for(std::vector<int>::reverse_iterator it = source->update.rbegin() ; it != source->update.rend(); ++it){ int uedgeid = *it;
//for(Cupdate *up = source->update;up; up = up->next){
    if((result = updateDFS(&edge[uedgeid], inc, depth+1, Emark))){
      if(DEBUG>=2 && result < 0){
	printf("updateDFS:source=(N%d,N%d,%0.3f,%d,%d):mark=%d (returning from cycle in update graph)\n",
	       nodeorder[node[source->Lmap].mapid],nodeorder[node[source->Rmap].mapid],source->distance,
	       source->Lflipped,source->Rflipped,mark);
	fflush(stdout);
      }
      break;
    }
  }

  Emark[edgeid] += inc;/* mark source as finished */
  if(DEBUG>=2) assert(Emark[edgeid] == (inc > 0 ? 2 : 0));
  
  return result;
}

/** find directed path from source to any target with mark == -1 using update links: 
   return 1 if path is found, -1 if cycle was found, 0 otherwise.

   Uses the mark field (assumed to be initialized at 0 except for targets)
 */
static int FindpathUpdate(Cedge *source, char *Emark)
{
  if(VERB>=2){
    printf("FindpathUpdate:source=x%p(N%d,N%d,%0.3f,%d,%d)\n",
	   (void *)source, nodeorder[node[source->Lmap].mapid],nodeorder[node[source->Rmap].mapid],source->distance, source->Lflipped,source->Rflipped);
    fflush(stdout);
  }

  int result = updateDFS(source,1,1,Emark);
  (void)updateDFS(source,-1,1,Emark);/* to reset mark field to 0 */
  return result;
}

/** recursively update updistance value of edge and its update successors to at least newdistance */
static void distanceUpdate(double newdistance, Cedge *pedge)
{
  if(newdistance <= pedge->updistance)
    return;
  pedge->updistance = newdistance;

  for(std::vector<int>::reverse_iterator it = pedge->update.rbegin() ; it != pedge->update.rend(); ++it)
    distanceUpdate(newdistance, &edge[*it]);
}

#if 0
/** If inc = 1, mark will be incremented from 0 to 1 and 2.
   If inc = -1, mark will be decremented from 2 to 1 and 0.
   Rdist is the distance from the root node (Rnode) to the current node (Snode)
   Only forward (+ve) edges are explored.

   Nodes with mark== -1 are skipped (disallowed nodes)
   Nodes with mark== -2 result in a failure value of -1.0 being returned.
       (If nodes with mark= -2 were set during the reverse DFS, then this means the forward and reverse
       DFS collide and the graph has a cycle with an odd number of map reversals)
*/
static double DFSforward(Cnode *Vnode, int Vflipped, int inc, double Rdist)
{
  if(DFSverb && inc > 0){
    printf("DFSforward(Vnode=N%d,flip=%d,inc=%d,Rdist=%0.3f):Vnode->mark=%d\n",
	   nodeorder[Vnode->mapid],Vflipped,inc,Rdist,Vnode->mark);
    fflush(stdout);
  }

  if(Vnode->mark == 1)/* cycle */
    return 0.0;
  if(Vnode->mark == inc+1/*(inc>0 ? 2 : 0)*/)/* this node already explored and finished */
    return 0.0;

  if(Vnode->mark == -1)/* disallowed node */
    return 0.0;
  if(Vnode->mark == -2)/* cycle with odd number of reversals detected */
    return 0.0;

  Vnode->flipDFS = Vflipped;
  Vnode->mark += inc;/* mark Vnode as active */
  if(DEBUG>=2) assert(Vnode->mark == 1);

  /* recurse over all +ve children (relative to Vnode being Vflipped) */
  /* If Vflipped==1: look for neighbors of Vnode with -ve nodedistance 
     If Vflipped==0: look for neighbors of Vnode with +ve node distance */

  double maxdist = Rdist;
  Cedge *Bedge = 0;
  for(int vw = (Vflipped ? 0 : Vnode->PosEdge); vw < (Vflipped ? Vnode->PosEdge : Vnode->edgecnt); vw++){
    Cedge *VWedge = &edge[Vnode->edgeid[vw]]; 
    Cnode *Wnode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
    int Wflipped = (Vflipped ^ (VWedge->Lflipped | VWedge->Rflipped));
    double Ndist = DFSforward(Wnode,Wflipped, inc, Rdist + VWedge->distance);
    if(Ndist < 0.0){
      maxdist = Ndist;
      Bedge = VWedge;
      break;
    }
    if(Ndist > maxdist){
      maxdist = Ndist;
      Bedge = VWedge;
    }
  }
  
  Vnode->mark += inc;/* mark Vnode as finished */
  if(DEBUG>=2) assert(Vnode->mark == (inc > 0 ? 2 : 0));
  Vnode->bestDFS = Bedge;

  if(DFSverb && inc > 0){
    if(Bedge)
      printf("DFSforward(Vnode=N%d,flip=%d,inc=%d,Rdist=%0.3f):final maxdist=%0.3f(mark=%d),best edge=(N%d,N%d,%0.3f,%d,%d)\n",
	     nodeorder[Vnode->mapid],Vflipped,inc,Rdist,maxdist,Vnode->mark,nodeorder[node[Bedge->Lmap].mapid],nodeorder[node[Bedge->Rmap].mapid],
	     Bedge->distance, Bedge->Lflipped,Bedge->Rflipped);
    else
      printf("DFSforward(Vnode=N%d,flip=%d,inc=%d,Rdist=%0.3f):final maxdist=%0.3f(mark=%d),best edge = NULL\n",
	     nodeorder[Vnode->mapid],Vflipped,inc,Rdist,maxdist,Vnode->mark);
    fflush(stdout);
  }
  return maxdist;
}
#endif

#if 0
/** find the longest path distance from node Rnode using +ve edges only (skip Nodes with mark= -1) */
static double FindForward(Cnode *Rnode,int reset)
{
  if(DFSverb){
    printf("FindForward(Rnode=N%d,reset=%d)\n",nodeorder[Rnode->mapid],reset);
    fflush(stdout);
  }
  if(reset)
    return DFSforward(Rnode,0,-1,0.0);/* to reset mark fields to 0 */
  return DFSforward(Rnode, 0, 1,0.0);
}


/** find the longest path distance from node Rnode using -ve edges only (skip Nodes with mark= -1) 
* return -1.0 if Nodes with mark = -2 are encountered : these are nodes marked during the forward sweep */
static double FindBackward(Cnode *Rnode, int reset)
{
  if(DFSverb){
    printf("FindBackward(Rnode=N%d,reset=%d)\n",nodeorder[Rnode->mapid],reset);
    fflush(stdout);
  }
  if(reset)
    return DFSforward(Rnode,1,-1,0.0);/* to reset mark fields to 0 */
  return DFSforward(Rnode, 1, 1,0.0);
}
#endif

#define DATA(SP, node, field) SP[node->nodeid].field

/** Here Heap means a Min-heap wrt heap[j].node->SPcost[dir][heap[j].flip]
   Pre : heap[1..heapcnt] is a Heap except for heap[i].node->SPcost[dir][flip] being (possibly) too large
   Post : heap[1..heapcnt] is a Heap.
   Backpointer heap[j].node->heapindex[heap[j].flip] must be maintained to equal j
   Children of heap[j] are heap[2j] and heap[2j+1].
 */
static inline void sift_down(Cnodeflip *heap,
			     int heapcnt,
			     int i,
			     int dir,
			     CSPdata *SP) /* Get per node data from SP[node->nodeid].heapindex[] instead of node->heapindex[] etc : used for multithreading */
{
  Cnodeflip heapI = heap[i];
  double Idist = SP[heapI.nodeid].SPcost[dir][heapI.flip];
  int j;
  while((j = 2*i) <= heapcnt){
    if(j < heapcnt && 
       SP[heap[j].nodeid].SPcost[dir][heap[j].flip] > SP[heap[j+1].nodeid].SPcost[dir][heap[j+1].flip]){
      if(Idist <= SP[heap[j+1].nodeid].SPcost[dir][heap[j+1].flip])
	break;
      heap[i] = heap[j+1];
      SP[heap[i].nodeid].heapindex[heap[i].flip] = i;
      i = j+1;
    } else {
      if(Idist <= SP[heap[j].nodeid].SPcost[dir][heap[j].flip])
	break;
      heap[i] = heap[j];
      SP[heap[i].nodeid].heapindex[heap[i].flip] = i;
      i = j;
    }
  }
  SP[heapI.nodeid].heapindex[heapI.flip] = i;
  heap[i] = heapI;
}

/** Here Heap means a Min-heap wrt heap[j]->hdist
   Pre : heap[1..heapcnt] is a Heap except for heap[i]->hdist being (possibly) too large
   Post : heap[1..heapcnt] is a Heap.
   Backpointer heap[j]->hindex must be maintained to equal j
   Children of heap[j] are heap[2j] and heap[2j+1].
 */
static inline void sift_down(Cnode **heap,
			     int heapcnt,
			     int i,
			     CSPdata *SP) /* Get per node data from SP[node->nodeid].hindex[] instead of node->hindex[] etc : used for multithreading */
{
  Cnode *heapI = heap[i];
  double Idist = DATA(SP,heapI,hdist);
  int j;
  while((j = 2*i) <= heapcnt){
    if(j < heapcnt && 
       DATA(SP,heap[j],hdist) > DATA(SP,heap[j+1],hdist)){
      if(Idist <= DATA(SP,heap[j+1],hdist))
	break;
      heap[i] = heap[j+1];
      DATA(SP,heap[i],hindex) = i;
      i = j+1;
    } else {
      if(Idist <= DATA(SP,heap[j],hdist))
	break;
      heap[i] = heap[j];
      DATA(SP,heap[i],hindex) = i;
      i = j;
    }
  }
  DATA(SP,heapI,hindex) = i;
  heap[i] = heapI;
}

/** Here Heap means a Min-heap wrt heap[i].node->SPcost[dir][heap[i].flip]

   Pre : heap[1..heapcnt] is a Heap except for heap[i].node->SPcost[dir][flip] being (possibly) too small
   Post : heap[1..heapcnt] is a Heap.
   Backpointer heap[j].node->heapindex[flip] must be maintained to equal j
   Parent of heap[j] is heap[j/2] 
 */

static void sift_up(Cnodeflip *heap,
		    int heapcnt,
		    int i,
		    int dir,
		    CSPdata *SP) /* IF != 0 : get per node data from SP[node->nodeid].heapindex[] instead of node->heapindex[] etc : used for multithreading */
{
  if(DEBUG>=2) assert(1 <= i && i <= heapcnt);
  Cnodeflip heapI = heap[i];
  double distI = SP[heapI.nodeid].SPcost[dir][heapI.flip];
  for(int j;i > 1; i = j){
    j = i/2;
    if(distI >= SP[heap[j].nodeid].SPcost[dir][heap[j].flip])
      break;
    heap[i] = heap[j];
    SP[heap[i].nodeid].heapindex[heap[i].flip] = i;
  }
  SP[heapI.nodeid].heapindex[heapI.flip] = i;
  heap[i] = heapI;
}

/** Here Heap means a Min-heap wrt heap[i]->SPdist[0][0]

   Pre : heap[1..heapcnt] is a Heap except for heap[i]->hdist being (possibly) too small
   Post : heap[1..heapcnt] is a Heap.
   Backpointer heap[j]->hindex[0] must be maintained to equal j
   Parent of heap[j] is heap[j/2] 
 */

static inline void sift_up(Cnode **heap,
		    int heapcnt,
		    int i,
		    CSPdata *SP) /* Get per node data from SP[node->nodeid].hindex[] instead of node->hindex[] etc : used for multithreading */
{
  if(PDEBUG>=2) assert(1 <= i && i <= heapcnt);
  Cnode *heapI = heap[i];
  double distI = DATA(SP,heapI,hdist);
  for(int j;i > 1; i = j){
    j = i/2;
    if(distI >= DATA(SP,heap[j],hdist))
      break;
    heap[i] = heap[j];
    DATA(SP,heap[i],hindex) = i;
  }
  DATA(SP,heapI,hindex) = i;
  heap[i] = heapI;
}



/** used to modify all nodes reachable from Vnode based on dir and rev args. 
   Also modifies the SPmark field as a side effect based on inc.
   Returns number of nodes encountered.
   Used by ShortestPaths()

   Gets per node data from TSP[node->nodeid].<field> instead of node-><field> to support multithreading
*/

/** this function takes 49% of runtime for human, and is over 98% of the cost of ShortestPaths :

   Has been modified to allow function cloning with only one int arg (V) remaining.
 */

/** This function takes over 25% of runtime for human, and is over 98% of the cost of ShortestPaths: see if DFS without recursion OR BFS is faster */

template<int inc, int dir, int rev>
//     const int inc; /* +1 or -1 depending on whether SPmark is to be incremented from 0 to 2 or decremented from 2 to 0 */
//     const int dir; /* 1 : reset SPcost[1] to -1.0 for all nodes reached
//                       0 : reset SPcost[0] to -1.0 for all nodes reached
//			 -1 (or rev=1) : no action. */
//     const int rev; /* If rev != 0, reverse the orientation of (node,flip) pairs by swapping SPcost[dir][0,1],SPdist[dir][0,1] and SPedge[dir][0,1] for all nodes reached */
static inline int SPinitDFS_orig(Cnode *Vnode, CSPdata *SP)
{
  CSPdata *pV = &SP[Vnode->nodeid];

  int SPmark = pV->SPmark & 3;
  if(SPmark == 1)/* cycle */
    return 0;
  if(SPmark == inc+1)/* this node already explored and finished */
    return 0;

  int cnt = 1;/* node count of all newly explored nodes */

  pV->SPmark += inc;/* SPmark Vnode as active */
  if(PDEBUG>=2) assert((pV->SPmark & 3) == 1);

  if(rev==0 && dir >= 0){
    double *p = pV->SPcost[dir];
    p[0] = p[1] = -1.0;
    if(VERB>=2 && DFSverb>=2){
      printf("SPinitDFS:Vnode=N%d,dir=%d,inc=%d,rev=%d:SPcost[dir][0]=%0.3f,SPcost[dir][1]=%0.3f\n",
	     nodeorder[Vnode->mapid],dir,inc,rev,pV->SPcost[dir][0],pV->SPcost[dir][1]);
      fflush(stdout);
    }
  }
  if(rev){
    double *pcost = pV->SPcost[dir];
    double cost = pcost[0];
    pcost[0] = pcost[1];
    pcost[1] = cost;

    double *pdist = pV->SPdist[dir];
    double dist = pdist[0];
    pdist[0] = pdist[1];
    pdist[1] = dist;

    Cedge **pedge = pV->SPedge[dir];
    Cedge *edge0 = pedge[0];
    pedge[0] = pedge[1];
    pedge[1] = edge0;
  }

  /* recurse over all children */
  int *Enodeid = Vnode->Enodeid;
  for(int vw = 0; vw < Vnode->edgecnt; vw++){
    if(PDEBUG>=2){
      Cedge *VWedge = &edge[Vnode->edgeid[vw]];
      assert(((VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap) == Enodeid[vw]);	
    }
    Cnode *Wnode = &node[Enodeid[vw]];
    cnt += SPinitDFS_orig<inc,dir,rev>(Wnode,SP);
  }

  pV->SPmark += inc;/* mark source as finished */
  if(PDEBUG>=2) assert((pV->SPmark & 3) == ((inc>0) ? 2 : 0));
  
  return cnt;
}

template<int inc, int dir, int rev>
//     const int inc; /* +1 or -1 depending on whether SPmark is to be incremented from 0 to 2 or decremented from 2 to 0 */
//     const int dir; /* 1 : reset SPcost[1] to -1.0 for all nodes reached
//                       0 : reset SPcost[0] to -1.0 for all nodes reached
//			 -1 (or rev=1) : no action. */
//     const int rev; /* If rev != 0, reverse the orientation of (node,flip) pairs by swapping SPcost[dir][0,1],SPdist[dir][0,1] and SPedge[dir][0,1] for all nodes reached */
static inline int SPinitDFS_node(int nodeid, CSPdata *SP)
{
  CSPdata *pV = &SP[nodeid];

  int SPmark = pV->SPmark & 3;
  if(SPmark == 1)/* cycle */
    return 0;
  if(SPmark == inc+1)/* this node already explored and finished */
    return 0;

  int cnt = 1;/* node count of all newly explored nodes */

  pV->SPmark += inc;/* mark Vnode as active */
  if(DEBUG>=2) assert((pV->SPmark & 3) == 1);

  if(rev){
    double *pcost = pV->SPcost[dir];
    double cost = pcost[0];
    pcost[0] = pcost[1];
    pcost[1] = cost;

    double *pdist = pV->SPdist[dir];
    double dist = pdist[0];
    pdist[0] = pdist[1];
    pdist[1] = dist;

    Cedge **pedge = pV->SPedge[dir];
    Cedge *edge0 = pedge[0];
    pedge[0] = pedge[1];
    pedge[1] = edge0;
  } else if(dir >= 0){
    double *p = pV->SPcost[dir];
    p[0] = p[1] = -1.0;
  }

  return cnt;
}

template<int inc, int dir, int rev>
//     const int inc; /* +1 or -1 depending on whether SPmark is to be incremented from 0 to 2 or decremented from 2 to 0 */
//     const int dir; /* 1 : reset SPcost[1] to -1.0 for all nodes reached
//                       0 : reset SPcost[0] to -1.0 for all nodes reached
//			 -1 (or rev=1) : no action. */
//     const int rev; /* If rev != 0, reverse the orientation of (node,flip) pairs by swapping SPcost[dir][0,1],SPdist[dir][0,1] and SPedge[dir][0,1] for all nodes reached */
static int SPinitDFS(int nodeid, CSPdata *SP, int nodes)// If nodes > 0 : only this many nodes are reachable (otherwise up to numnodes nodes may be reachable)
{
  /* non-recursive DFS function that uses Vnode[] and Vw[] to stack recursive variables nodeid, w */
  if(nodes <= 0) 
    nodes = numnodes;

  size_t memsiz = nodes * (sizeof(int) + sizeof(int));
  char *mem = (char *)malloc(memsiz);
  if(!mem){
    printf("Error allocating %d structs of %lu bytes totalling %lu bytes\n",nodes, sizeof(int) + sizeof(int), memsiz);
    fflush(stdout);exit(1);
  }
  int *Vnodeid = (int *)mem;
  int *Vw = (int *)&mem[nodes * sizeof(int)];

  int last = 0;
  Vnodeid[last] = nodeid;
  int w = Vw[last] = 0;

  int cnt = SPinitDFS_node<inc,dir,rev>(nodeid,SP);

  while(1) {
  Lcall:
    CSPdata *pW = &SP[nodeid];
    int *Enodeid = pW->Enodeid;

    for(; w < pW->hindex;) {
      int nodeid2 = Enodeid[w++];
      int node_cnt = SPinitDFS_node<inc,dir,rev>(nodeid2,SP);
      if(node_cnt == 0)
	continue;

      cnt += node_cnt;

      /* simulate recursive call */
      Vw[last++] = w;

      Vnodeid[last] = nodeid = nodeid2;
      Vw[last] = w = 0;
      goto Lcall;
    }
		
    /* Roll back stack when we are done */
    do {
      pW = &SP[nodeid];
      pW->SPmark += inc;
		
      if(--last < 0) { /* All done */
	free(mem);
	return cnt;
      }

      nodeid = Vnodeid[last];
      w = Vw[last];
      pW = &SP[nodeid];
    } while(w >= pW->hindex);
  }
  
  /* The following code should never run */
  free(mem);
  return cnt;
}

template<int inc, int dir, int rev>
//     const int inc; /* +1 or -1 depending on whether SPmark is to be incremented from 0 to 2 or decremented from 2 to 0 for all nodes reached */
//     const int dir; /* 1 : reset SPcost[1][0..1] to -1.0 for all nodes reached
//                       0 : reset SPcost[0][0..1] to -1.0 for all nodes reached
//			 -1 (or rev=1) : no action. */
//     const int rev; /* If rev != 0, reverse the orientation of (node,flip) pairs by swapping SPcost[dir][0,1],SPdist[dir][0,1] and SPedge[dir][0,1] for all nodes reached */
static int SPinitBFS(Cnode *Rnode, CSPdata *SP, double maxdist)// only initialize nodes until cumulative distance exceeds maxdist
{
  if(PDEBUG) assert(maxdist > 0.0);

  Cnode **heap = new Cnode*[numnodes+1];
  int heapcnt = DATA(SP,Rnode,hindex) = 1;
  heap[1] = Rnode;
  DATA(SP,Rnode,hdist) = 0.0;
  SP[Rnode->nodeid].SPmark += inc;/* mark Rnode as active (in heap) */

  int cnt = 0;/* node count of all nodes reached (and updated) */

  while(heapcnt > 0){
    cnt++;

    /* remove top node from heap */
    Cnode **pUnode = &heap[1];
    Cnode *Unode = *pUnode;
    
    if(PDEBUG>=2) assert(DATA(SP,Unode,hindex) == 1);
    DATA(SP,Unode,hindex) = -1;

    if(--heapcnt > 0){
      heap[1] = heap[heapcnt+1];
      DATA(SP,heap[1],hindex) = 1;
      if(heapcnt > 1)
	sift_down(heap,heapcnt,1,SP);
    }

    /* mark Unode as finished (removed from heap) */
    CSPdata *pU = &SP[Unode->nodeid];
    pU->SPmark += inc;
    if(PDEBUG>=2) assert((pU->SPmark & 3) == ((inc>0) ? 2 : 0));

    /* update Unode depending on rev && dir */
    if(rev==0 && dir >= 0){
      double *p = pU->SPcost[dir];
      p[0] = p[1] = -1.0;
      if(VERB>=2 && DFSverb>=2){
	printf("SPinitBFS:Unode=N%d,dir=%d,inc=%d,rev=%d:SPcost[dir][0]=%0.3f,SPcost[dir][1]=%0.3f\n",
	       nodeorder[Unode->mapid],dir,inc,rev,pU->SPcost[dir][0],pU->SPcost[dir][1]);
	fflush(stdout);
      }
      if(DEBUG>=3){/* mark SPdist[dir][0,1] as uninitialized */
	p = pU->SPdist[dir];
	p[0] = p[1] = -1.0;
      }
    }
    if(rev){
      double *pcost = pU->SPcost[dir];
      double cost = pcost[0];
      pcost[0] = pcost[1];
      pcost[1] = cost;
      
      double *pdist = pU->SPdist[dir];
      double dist = pdist[0];
      pdist[0] = pdist[1];
      pdist[1] = dist;

      Cedge **pedge = pU->SPedge[dir];
      Cedge *edge0 = pedge[0];
      pedge[0] = pedge[1];
      pedge[1] = edge0;
    }
    
    if(DATA(SP,Unode,hdist) > maxdist)
      continue;
    
    double Ucost = pU->hdist;
    if(PDEBUG>=2) assert(Ucost >= 0.0);

    /* find all successors of Unode */
    int *Enodeid = Unode->Enodeid;
    for(int v = 0; v < Unode->edgecnt; v++){
      Cedge *UVedge = &edge[Unode->edgeid[v]];
      //      if(UVedge->valid <= 0)
      //	continue;
      if(PDEBUG>=2) assert(((UVedge->Lmap == Unode->nodeid) ? UVedge->Rmap : UVedge->Lmap) == Enodeid[v]);
      Cnode *Vnode = &node[Enodeid[v]];
      CSPdata *pV = &SP[Vnode->nodeid];
      int SPmark = pV->SPmark & 3;
      if(SPmark == inc+1)/* this node already explored and finished */
	continue;

      double cost = Ucost + UVedge->distance;
      if(PDEBUG>=2) assert(cost >= 0.0);
      if(SPmark != 1){/* Vnode was previously unreached : add it to the heap */
	pV->hdist = cost;
	pV->SPmark += inc;/* mark Vnode as active (in heap) */
	if(PDEBUG>=2) assert((pV->SPmark & 3) ==1);
	Cnode **pheap = &heap[++heapcnt];
	if(PDEBUG>=2) assert(heapcnt <= numnodes);
	pheap[0] = Vnode;
	pV->hindex = heapcnt;
	sift_up(heap,heapcnt,heapcnt,SP);
      } else if(cost < pV->hdist){/* update Vnode in heap */
	pV->hdist = cost;
	int hindex = pV->hindex;
	if(PDEBUG>=2) assert(hindex >= 1);/* Vnode should still be in heap */
	Cnode **pheap = &heap[hindex];
	if(PDEBUG>=2) assert(pheap[0] == Vnode);
	sift_up(heap,heapcnt,hindex,SP);
      }
    } /* for all neighbors Vnode of Unode */
  } // while(heapcnt > 0)
  
  delete [] heap;

  return cnt;
}

/** If Vnodeid >= 0 : check if path back from (Unode,Uflipped) to source node traverses Vnode (any orientation).
                      If so return (Vnode,Vflipped), else return (0,0)
   If Vnodeid < 0 : check if path back from (Unode,Uflipped) to source node traverses any node twice.
                     If so return (node,flip) for occurance closest to source node, else return (0,0).
   sets SPmark |= 4 when a node is traversed with flip=0 and SPmark |= 8 when node is traversed with flip=1.
   SPmark is reset when this recursive function unwinds */

static Cnodeflip loopback(int Unodeid,
			  int Uflipped,
			  int Vnodeid,
			  int dir, /**< dir index of path edges : path is always traversed in the forward direction from source to Unode */
			  CSPdata *SP, /**< Get per node data from SP[nodeid].heapindex[] instead of node[nodeid] : used for multithreading */
			  int nodes = 0, /**< If != 0 : number of nodes reachable from Rnode (will be <= numnodes) */
			  Cnode **CCnodes = 0, /**< If != 0 : array of node pointers indexed by nodeid : Used when index into SP[] is NOT same as index into node[], in which case SP[] index range is 0 .. nodes-1 */
			  int *CCnodeid = 0 /** < If != 0 : CCnodeid[Vnode->nodeid] is the index of Vnode in SP[] : Used when index into SP[] is NOT same index into node[], in which case SP[] index range is 0 .. nodes-1 */
			  )
{
  if(DEBUG>=2 && CCnodes){
    assert(nodes > 0);
    assert(CCnodeid != NULL);
    assert(0 <= Unodeid && Unodeid < nodes);
    if(Vnodeid >= 0)
      assert(Vnodeid < nodes);
  }

  Cnode *Unode = CCnodes ? CCnodes[Unodeid] : &node[Unodeid];
  Cnode *Vnode = (Vnodeid < 0) ? 0 : CCnodes ? CCnodes[Vnodeid] : &node[Vnodeid];

  Cnodeflip ret;
  if(DEBUG>=2 && Vnode) assert(!(SP[Unodeid].SPmark & (4|8)));

  if(VERB>=2 && DFSverb>=2 && !Vnode){
    Cedge *pedge = SP[Unodeid].SPedge[dir][Uflipped];
    if(pedge)
      printf("loopback(Unode=N%d%c,dir=%d):SPmark=%d,cost=%0.3f,pedge=(N%d,N%d,%0.3f,%d,%d)\n",
	     nodeorder[Unode->mapid],Uflipped ? '\'':' ',dir,SP[Unodeid].SPmark & (4|8), SP[Unodeid].SPcost[dir][Uflipped],
	     nodeorder[node[pedge->Lmap].mapid],nodeorder[node[pedge->Rmap].mapid],pedge->distance,pedge->Lflipped,pedge->Rflipped);
    else
      printf("loopback(Unode=N%d%c,dir=%d):SPmark=%d,cost=%0.3f,pedge=0\n",
	     nodeorder[Unode->mapid],Uflipped ? '\'':' ',dir, SP[Unodeid].SPmark & (4|8),  SP[Unodeid].SPcost[dir][Uflipped]);
    fflush(stdout);
  }

  if(Unodeid == Vnodeid){
    ret.nodeid = Unodeid;
    ret.flip = Uflipped;
    return ret;
  }

  if(!Vnode && (SP[Unodeid].SPmark & (4|8))){/* this node already traversed */
    ret.nodeid = Unodeid;
    ret.flip = Uflipped;
    if(VERB>=2 && DFSverb>=2 && !Vnode){
      printf("loopback(Unode=N%d%c,dir=%d):returning Unode due to SPmark\n",
	     nodeorder[Unode->mapid],Uflipped ? '\'':' ',dir);
      fflush(stdout);
    }
    return ret;
  }

  Cedge *Uedge = SP[Unodeid].SPedge[dir][Uflipped];
  if(!Uedge){/* Unode is the source node */
    if(DEBUG>=2 && !(SP[Unodeid].SPcost[dir][Uflipped] == 0.0)){
      printf("loopback(Unode=N%d%c,dir=%d):SPmark=%d,cost=%0.3f,pedge=0\n",
	     nodeorder[Unode->mapid],Uflipped ? '\'':' ',dir,DATA(SP,Unode,SPmark) & (4|8),  DATA(SP,Unode,SPcost)[dir][Uflipped]);
      printf("cost should be 0 when Unode == source node\n");
      fflush(stdout);
      assert(SP[Unodeid].SPcost[dir][Uflipped] == 0.0);
    }
    ret.nodeid = -1;
    ret.flip = 0;
    if(VERB>=2 && DFSverb>=2 && !Vnode){
      printf("loopback(Unode=N%d%c,dir=%d):returning 0 at source node\n",
	     nodeorder[Unode->mapid],Uflipped ? '\'':' ',dir);
      fflush(stdout);
    }
    return ret;
  }
  int Pnodeid = (Uedge->Lmap == Unode->nodeid) ? Uedge->Rmap : Uedge->Lmap;
  Cnode *Pnode = &node[Pnodeid];
  if(CCnodeid){
    int nPnodeid = CCnodeid[Pnodeid];
    if(DEBUG>=2 && !(0 <= nPnodeid && nPnodeid < nodes)){
      printf("loopback(nodes=%d): Unodeid=%d, Unode->nodeid=%d,Uedge->Lmap=%d,Uedge->Rmap=%d,Pnodeid=%d,nPnodeid=%d\n",
	     nodes,Unodeid,Unode->nodeid,Uedge->Lmap,Uedge->Rmap,Pnodeid,nPnodeid);
      fflush(stdout);
      assert(0 <= nPnodeid  && nPnodeid < nodes);
    }
    Pnodeid = nPnodeid;
  }
  /* If Uflipped==0 : look for neighbors with -ve nodedistance (backwards)
     If Uflipped==1 : look for neighbors with +ve nodedistance (backwards) */
  int Pflipped = (Uflipped ^ (Uedge->Lflipped | Uedge->Rflipped));

  if(VERB>=2 && DFSverb>=2 && !Vnode){
    printf("  loopback(Unode=N%d%c,dir=%d):Pnode=N%d%c,Uedge=(N%d,N%d,%0.3f,%d,%d),dist=%0.3f\n",
	   nodeorder[Unode->mapid],Uflipped ? '\'':' ',dir,
	   nodeorder[Pnode->mapid],Pflipped ? '\'':' ',
	   nodeorder[node[Uedge->Lmap].mapid],nodeorder[node[Uedge->Rmap].mapid],Uedge->distance,Uedge->Lflipped,Uedge->Rflipped,
	   SP[Unodeid].SPcost[dir][Uflipped]);
    fflush(stdout);
  }

  if(Vnode)
    return loopback(Pnodeid,Pflipped,Vnodeid,dir,SP,nodes,CCnodes,CCnodeid);

  SP[Unodeid].SPmark |= Uflipped ? 4 : 8;

  ret = loopback(Pnodeid,Pflipped,-1,dir,SP,nodes,CCnodes,CCnodeid);

  SP[Unodeid].SPmark &= ~(Uflipped ? 4 : 8);

  if(VERB>=2 && DFSverb>=2 && !Vnode){
    if(ret.nodeid >= 0){
      Cnode *Rnode = CCnodes ? CCnodes[ret.nodeid] : &node[ret.nodeid];
      printf("loopback(Unode=N%d%c,dir=%d):ret = N%d%c,SPmark=%d\n",
	     nodeorder[Unode->mapid],Uflipped ? '\'':' ',dir,
	     nodeorder[Rnode->mapid],ret.flip ? '\'':' ', SP[Unodeid].SPmark & (4|8));
    } else
      printf("loopback(Unode=N%d%c,dir=%d):ret = 0,SPmark=%d\n",
	     nodeorder[Unode->mapid],Uflipped ? '\'':' ',dir, SP[Unodeid].SPmark & (4|8));
    fflush(stdout);
  }
  
  return ret;
}

static int cyclecnt = 0;
static int pathcrosscnt = 0;

static int cyclewarning = 1;
static int pathcrosswarning = 1;

/** mark all nodes reachable by forward (or backward see dir) links from Rnode 
    (with orientation Rflipped) by setting the SPcost,SPedge and SPdist fields
   of all reachable nodes using Disjkstra's shortest path (lowest cost) algorithm. 
   Each Cnode corresponds to 2 path nodes (one for each orientation (flip status)). 
   A map may be reachable in one orientation, but not the other.
   If reachable in both orientations, the lowest cost path may be different for each orientation.
   (In general if a map is reachable in both orientations there is a bad cycle in the graph and a Warning will be printed)
   maps/nodes marked with (SPmark & 16) and edges with valid <=0 may not be traversed.
   
   Also keep track of the actual path distance (in case the cost used is not the same as distance) and returns the most distant node from Rnode along the lowest cost paths.

   returns most distant reachable (node,flip) pair 
*/

Cnodeflip ShortestPaths(int Rnodeid, /**< source node index */
			int Rflipped,/**< orientation of source node */
			int dir,/**< 0 for forward links, 1 for backward links */
			int metric,/**< cost/distance metric to use : 
					     metric = 0 : cost= Sum(UVedge->distance)    distance = Sum(UVedge->distance)
					     metric = 1 : cost= Sum(1000.0/UVedge->cov)  (Don't check for cycles or track max distance)
					     metric = 2 : cost= Sum((AVCOV ? UVedge->distance : 1000.0)/UVedge->cov)  distance = Sum(UVedge->distance * (AVCOV ? 1.0 : UVedge->cov)) */
			double mincov,/**< minimum coverage of edges to be traversed */
			double maxcost, /* If != 0 && metric==0 : stop after reaching node with path distance from Rnode exceeding this value */
			CSPdata *SP, /**< Modify/Write per node data to SP[node->nodeid].heapindex[] instead of node-> : used for multithreading */
			Cnodeflip *heapmem = 0, /* If != 0 : preallocated memory heap[0..2*numnodes] */
			int nodes = 0, /**< If != 0 : number of nodes reachable from Rnode (will be <= numnodes) */
			Cnode **CCnodes = 0, /**< If != 0 : array of node pointers indexed by nodeid : Used when index into SP[] is NOT same as index into node[], in which case SP[] index range is 0 .. nodes-1 */
			int *CCnodeid = 0/** < If != 0 : CCnodeid[Vnode->nodeid] is the index of Vnode in SP[] : Used when index into SP[] is NOT same as index into node[], in which case SP[] index range is 0 .. nodes-1 */
			)
{
  if(DEBUG && CCnodes) assert(nodes > 0 && CCnodeid != 0 && 0 <= Rnodeid && Rnodeid < nodes);

  //  int tid = 0;
#ifdef _OPENMP
  //  tid = omp_get_thread_num ();
#endif

  Cnode *Rnode = CCnodes ? CCnodes[Rnodeid] : &node[Rnodeid];

  if(VERB>=2 && (Rnodeid==148965 || Rnodeid==1240440) && Rflipped==0 && dir==1 && metric==0 /*  && tid==0 && DFSverb*/){
    printf("Starting ShortestPaths from Rnodeid=%d,node=%d(N%d), Rflipped=%d, dir=%d, metric=%d, mincov= %0.2f, maxcost= %0.2f: wtime=%0.6f\n",
	   Rnodeid,Rnode->nodeid,nodeorder[Rnode->mapid],Rflipped,dir,metric,mincov,maxcost,wtime());
    fflush(stdout);
  }
  if(DEBUG>=2) assert(Rnode->nodeid == Rnode - node);

  /* A (node,flip) pair is in one of three states:
       1. Unreached : node->SPcost[dir][flip] = -1.0 (logically same as nodes in min-Heap with node->SPcost[dir][flip] = Inf)
       2. Inserted into min-Heap : node->SPcost[dir][flip] >= 0.0, node->heapindex[flip] >= 1
       3. Removed from min-Heap : node->SPcost[dir][flip] >= 0.0, node->heapindex[flip] = -1  */

  if(DEBUG>=2 && (SP[Rnodeid].SPmark & 16)){
    printf("ShortestPaths:Root node %d(N%d) has SPmark=%d\n",Rnode->nodeid, nodeorder[Rnode->mapid],Rnode->SPmark);
    fflush(stdout);
    assert(0);
  }
  if(DEBUG && maxcost > 0.0 && metric != 0){
    printf("ShortestPaths not supported with metric=%d maxcost=%0.1f\n",metric,maxcost);
    fflush(stdout); exit(1);
  }

  /* first reset all reachable nodes (in either direction) to SPcost[dir][0..1] = -1.0 */
  //  if(DEBUG>=2 && DFSverb)
  //    DFSverb = 2;
  int nodecnt = (BFSINIT && maxcost > 0.0) ? (dir ? SPinitBFS<1,1,0>(Rnode,SP,maxcost) : SPinitBFS<1,0,0>(Rnode,SP,maxcost)) :
    (dir ? SPinitDFS<1,1,0>(Rnodeid,SP,nodes) : SPinitDFS<1,0,0>(Rnodeid,SP,nodes));
  if(DEBUG>=2) assert(nodecnt >= 1);
  if(DEBUG>=2 && nodes > 0) assert(nodecnt <= nodes);

  /* decide on source node orientation : If dir = 1, we invert the source node orientation,
     traverse in the forward direction, then reverse all SPcost[1][0,1] values */
  int flip = (dir ? (!Rflipped) : Rflipped);

  /* initialize the heap[1..heapcnt] with the Root node  (Rnode,flip) */
  Cnodeflip *heap = heapmem ? heapmem : new Cnodeflip[2*nodecnt+1];
  int heapcnt = SP[Rnodeid].heapindex[flip] = 1;
  heap[1].nodeid = Rnodeid;
  heap[1].flip = flip;
  SP[Rnodeid].SPcost[dir][flip] = 0.0;
  SP[Rnodeid].SPdist[dir][flip] = 0.0;
  SP[Rnodeid].SPedge[dir][flip] = 0;
  if(VERB>=2 && (Rnodeid==148965 || Rnodeid==1240440) && Rflipped==0 && dir==1 && metric==0 /*  && tid==0 && DFSverb*/){
    printf("source node %d(N%d%c),flip=%d(dir=%d):nodecnt=%d,wtime= %0.6f\n",
	   Rnode->nodeid, nodeorder[Rnode->mapid],Rflipped ? '\'': ' ',flip,dir,nodecnt,wtime());
    fflush(stdout);
  }

  /* keep track of most distant completed nodes (exlude paths that loopback on themselves) */
  double maxdist = 0.0;
  Cnodeflip ret;
  ret.nodeid = Rnodeid;
  ret.flip = flip;

  while(heapcnt > 0){
    /* remove top node from heap */
    Cnodeflip *Unodeflip = &heap[1];
    int Unodeid = Unodeflip->nodeid;
    int Uflipped = Unodeflip->flip;
    Cnode *Unode = CCnodes ? CCnodes[Unodeid] : &node[Unodeid];
    if(DEBUG>=2) assert(Unode->nodeid == Unode - node);
    if(DEBUG>=2) assert(SP[Unodeid].heapindex[Uflipped] == 1);
    SP[Unodeid].heapindex[Uflipped] = -1;
    if(--heapcnt > 0){
      heap[1] = heap[heapcnt+1];
      SP[heap[1].nodeid].heapindex[heap[1].flip] = 1;
      if(heapcnt > 1)
	sift_down(heap,heapcnt,1,dir,SP);
    }

    if(VERB>=2 && (Rnodeid==148965 || Rnodeid==1240440) && Rflipped==0 && dir==1 && metric==0 /*  && tid==0 && DFSverb*/){
      Cedge *Uedge = SP[Unodeid].SPedge[dir][Uflipped];
      if(Uedge)
	printf("removed from heap:node%d(N%d%c),cost=%0.3f,dist=%0.3f,Uedge=(N%d,N%d,%0.3f,%d,%d,cov=%0.3f)(Unode->heapindex[%d]=%d):wtime= %0.6f\n",
	       Unode->nodeid,nodeorder[Unode->mapid],Uflipped ? '\'':' ',SP[Unodeid].SPcost[dir][Uflipped],SP[Unodeid].SPdist[dir][Uflipped],
	       nodeorder[node[Uedge->Lmap].mapid], nodeorder[node[Uedge->Rmap].mapid],Uedge->distance,Uedge->Lflipped,Uedge->Rflipped,Uedge->cov,
	       Uflipped,SP[Unodeid].heapindex[Uflipped],wtime());
      else
	printf("removed from heap:node%d(N%d%c),cost=%0.3f,dist=%0.3f,edge=0:wtime= %0.6f\n",
	       Unode->nodeid,nodeorder[Unode->mapid],Uflipped ? '\'':' ',SP[Unodeid].SPcost[dir][Uflipped],SP[Unodeid].SPdist[dir][Uflipped],wtime());
      fflush(stdout);
    }

    /* check if Unode is most distant completed node so far */
    if(metric != 1 && SP[Unodeid].SPdist[dir][Uflipped] > maxdist){
      Cnodeflip loop = loopback(Unodeid,Uflipped,-1,dir,SP,nodes,CCnodes,CCnodeid);
      if(loop.nodeid >= 0){
	if(DEBUG>=2 && CCnodes) assert(loop.nodeid < nodes);
	Cnode *Lnode = CCnodes ? CCnodes[loop.nodeid] : &node[loop.nodeid];
	if(DEBUG>=2) assert(Lnode->nodeid == Lnode - node);
	pathcrosscnt = 1;
	if(!pathcrosswarning){
	  printf("Warning: path from N%d%c to N%d%c (dist=%0.3f,cost=%0.3f) traverses N%d twice (first flip=%d)\n",
		 nodeorder[Rnode->mapid],flip ? '\'':' ', nodeorder[Unode->mapid],Uflipped ? '\'':' ',
		 SP[Unodeid].SPdist[dir][Uflipped], SP[Unodeid].SPcost[dir][Uflipped], nodeorder[Lnode->mapid],loop.flip);
	  fflush(stdout);
	  pathcrosswarning = 1;
	}
      } else {
	maxdist = SP[Unodeid].SPdist[dir][Uflipped];
	ret.nodeid = Unodeid;
	ret.flip = Uflipped;
	if(VERB>=2 && (Rnodeid==148965 || Rnodeid==1240440) && Rflipped==0 && dir==1 && metric==0 /*  && tid==0 && DFSverb*/){
	  printf("updating most distant node to %d(N%d%c), distance=%0.3f:wtime= %0.6f\n",
		 Unode->nodeid,nodeorder[Unode->mapid],Uflipped ? '\'':' ',maxdist,wtime());
	  fflush(stdout);
	}
	if(metric==0 && maxcost > 0.0 && maxdist > maxcost)
	  break;
      }
    }

    if(metric==0 && maxcost > 0.0 && SP[Unodeid].SPdist[dir][Uflipped] > maxcost)
      continue;

    /* find all successors of Unode (except those marked with (SPmark & 16) */
    /* If Uflipped==1 : look for neighbors with -ve nodedistance
       If Uflipped==0 : look for neighbors with +ve nodedistance */
    int v;
    if(VERB>=2 && (Rnodeid==148965 || Rnodeid==1240440) && Rflipped==0 && dir==1 && metric==0 /*  && tid==0 && DFSverb*/){
      printf("Unode->PosEdge=%d,Unode->edgecnt=%d,Uflipped=%d:wtime=%0.6f\n",Unode->PosEdge,Unode->edgecnt,Uflipped,wtime());
      for(int v=0; v < Unode->edgecnt; v++){
	Cedge *UVedge = &edge[Unode->edgeid[v]];
	printf("\t v=%d,UVedge=%d,%d(N%d,N%d,%d,%d):valid=%d\n",
	       v,UVedge->Lmap,UVedge->Rmap,nodeorder[node[UVedge->Lmap].mapid],nodeorder[node[UVedge->Rmap].mapid],UVedge->Lflipped,UVedge->Rflipped,UVedge->valid);
      }
      fflush(stdout);
    }
    if(DEBUG>=2){
      for(int v = 0; v < Unode->edgecnt; v++){
	Cedge *UVedge = &edge[Unode->edgeid[v]];
	if(UVedge->valid <= 0)
	  continue;
	if(!((UVedge->Lmap == Unode->nodeid || UVedge->Rmap == Unode->nodeid) && (UVedge->Lmap == Unode->nodeid ? (Unode->Enodeid[v] == UVedge->Rmap) : (Unode->Enodeid[v] == UVedge->Lmap)))){
          printf("Unode->nodeid=%d, v=%d/%d: valid=%d,Lmap=%d,Rmap=%d,node[i].Enodeid[v]= %d\n",
		 Unode->nodeid,v,Unode->edgecnt,UVedge->valid,UVedge->Lmap,UVedge->Rmap,Unode->Enodeid[v]);
	  fflush(stdout);

          assert(UVedge->Lmap == Unode->nodeid || UVedge->Rmap == Unode->nodeid);
	  assert(UVedge->Lmap == Unode->nodeid ? (Unode->Enodeid[v] == UVedge->Rmap) : (Unode->Enodeid[v] == UVedge->Lmap));
        }
      }
    }

    for(v = (Uflipped ? 0 : Unode->PosEdge); v < (Uflipped ? Unode->PosEdge : Unode->edgecnt); v++){
      Cedge *UVedge = &edge[Unode->edgeid[v]];
      if(VERB>=2 && (Rnodeid==148965 || Rnodeid==1240440) && Rflipped==0 && dir==1 && metric==0 /*  && tid==0 && DFSverb*/){
	printf("v=%d,UVedge=%d,%d(N%d,N%d,%d,%d):valid=%d,cov=%0.1f(mincov=%0.1f):wtime=%0.6f\n",
	       v,UVedge->Lmap,UVedge->Rmap,nodeorder[node[UVedge->Lmap].mapid],nodeorder[node[UVedge->Rmap].mapid],UVedge->Lflipped,UVedge->Rflipped,UVedge->valid,UVedge->cov,mincov,wtime());
	fflush(stdout);
      }

      if(UVedge->valid <= 0)
	continue;
      if(UVedge->cov < mincov)
	continue;
      Cnode *Vnode = &node[(UVedge->Lmap == Unode->nodeid) ? UVedge->Rmap : UVedge->Lmap];
      if(DEBUG>=2) assert(Vnode->nodeid == Vnode - node);
      int Vnodeid = CCnodeid ? CCnodeid[Vnode->nodeid] : Vnode->nodeid;
      if(DEBUG>=2 && CCnodeid) assert(0 <= Vnodeid && Vnodeid < nodes);
      if((SP[Vnodeid].SPmark & 16))
	continue;
      int Vflipped = (Uflipped ^ (UVedge->Lflipped | UVedge->Rflipped));
      double distance = SP[Unodeid].SPdist[dir][Uflipped];
      double cost = SP[Unodeid].SPcost[dir][Uflipped];
      if(PDEBUG>=2) assert(cost >= 0.0);

      switch(metric){
      case 0:
	cost +=  UVedge->distance;
	distance +=  UVedge->distance;
	break;
      case 1:
	cost += 1000.0/UVedge->cov;
	distance +=  UVedge->distance;
	break;
      case 2:
	cost += (AVCOV ? UVedge->distance : 1000.0)/UVedge->cov;
	distance +=  AVCOV ? UVedge->distance : UVedge->distance * UVedge->cov;
	break;
      default:
	printf("ShortestPaths: Invalid cost metric=%d(only 0,1 or 2 supported)\n",metric);
	exit(1);
      }
      if(DEBUG>=2 && !(cost >= 0.0)){/* Dijkstra's algorithm only supports non-negative edge costs */
	printf("Unode=%d(N%d%c),Vnode=%d(N%d%c),dir=%d:Old cost=%0.3f,UVedge=(cov=%0.3f,dist=%0.3f,valid=%d,mark=%d)total cost=%0.3f\n",
	       Unode->nodeid,nodeorder[Unode->mapid],Uflipped ? '\'':' ', Vnode->nodeid,nodeorder[Vnode->mapid],Vflipped ? '\'':' ',
	       dir, SP[Unodeid].SPcost[dir][Uflipped], UVedge->cov,UVedge->distance, UVedge->valid, UVedge->mark, cost);
	fflush(stdout);
	assert(cost >= 0.0);
      }
      if(metric != 1 && (TSAN || !cyclecnt)){      /* check if current path loops back through node Vnode : if so print warning, if not done already */
	int mycyclecnt;

        #pragma omp atomic capture
	mycyclecnt = cyclecnt++;// NOTE : global cyclecnt

	if(!mycyclecnt){
	  Cnodeflip loop = loopback(Unodeid,Uflipped,Vnodeid,dir,SP,nodes,CCnodes,CCnodeid);
	  if(loop.nodeid == Vnodeid){
            #pragma omp critical(cyclewarning)
	    {
	      if(!cyclewarning){
		printf("Warning:graph contains cycle from N%d%c to (%0.3f kb) N%d%c to (%0.3f kb) N%d%c (Rnode=N%d%c)\n",
		     nodeorder[Vnode->mapid], loop.flip ? '\'':' ', 
		     SP[Unodeid].SPdist[dir][Uflipped] - SP[Vnodeid].SPdist[dir][loop.flip],
		     nodeorder[Unode->mapid], Uflipped ? '\'':' ',
		     UVedge->distance,
		     nodeorder[Vnode->mapid], Vflipped ? '\'':' ',
		     nodeorder[Rnode->mapid], flip ? '\'':' ');
		fflush(stdout);
		cyclewarning = 1;
	      }
	    }
	  } else {
	    #pragma omp atomic
	    cyclecnt--;//NOTE : global cyclecnt
	  }
	}
      }
      if(VERB>=2 && (Rnodeid==148965 || Rnodeid==1240440) && Rflipped==0 && dir==1 && metric==0 /*  && tid==0 && DFSverb*/){
	printf("%s heap:node%d(N%d%c),cost=%0.3f,dist=%0.3f,edge=%d,%d(N%d,N%d,%0.3f,%d,%d,cov=%0.3f):v=%d,PosEdge=%d,edgecnt=%d(Vnode->heapindex[%d]=%d,old cost=%0.3f,dist=%0.3f):wtime=%0.6f\n",
	       SP[Vnodeid].SPcost[dir][Vflipped] < 0.0 ? "inserted into" : 
	       SP[Vnodeid].SPcost[dir][Vflipped] > cost ? "updated in" : "not updated in",
	       Vnode->nodeid,nodeorder[Vnode->mapid],Vflipped ? '\'':' ',cost,distance,
	       UVedge->Lmap,UVedge->Rmap,nodeorder[node[UVedge->Lmap].mapid],nodeorder[node[UVedge->Rmap].mapid],UVedge->distance,UVedge->Lflipped,UVedge->Rflipped,UVedge->cov,
	       v, Unode->PosEdge, Unode->edgecnt, Vflipped, SP[Vnodeid].heapindex[Vflipped], SP[Vnodeid].SPcost[dir][Vflipped], SP[Vnodeid].SPdist[dir][Vflipped],wtime());
	fflush(stdout);
      }
      if(DEBUG>=2) assert((SP[Vnodeid].SPmark & 3) == 2);
      if(DEBUG>=3) assert(SP[Vnodeid].SPcost[dir][Vflipped] >= -1.000001);
      if(SP[Vnodeid].SPcost[dir][Vflipped] < 0.0){/* (Vnode,Vflipped) was previously unreached : add it to the heap */
	SP[Vnodeid].SPcost[dir][Vflipped] = cost;
	SP[Vnodeid].SPdist[dir][Vflipped] = distance;
	SP[Vnodeid].SPedge[dir][Vflipped] = UVedge;
	Cnodeflip *pheap = &heap[++heapcnt];
	if(DEBUG>=2 && !(heapcnt <= nodecnt*2)){
	  printf("heapcnt=%d,nodecnt=%d\n",heapcnt,nodecnt);
	  fflush(stdout);
	  assert(heapcnt <= 2*nodecnt);
	}
	pheap->nodeid = Vnodeid;
	pheap->flip = Vflipped;
	SP[Vnodeid].heapindex[Vflipped] = heapcnt;
	sift_up(heap,heapcnt,heapcnt,dir,SP);
      } else if(cost < SP[Vnodeid].SPcost[dir][Vflipped]){/* update (Vnode,Vflipped) in heap */
	int heapindex = SP[Vnodeid].heapindex[Vflipped];
	if(DEBUG>=2 && !(heapindex >= 1 && heapindex <= heapcnt)){
	  printf("Rnode=%d(N%d%c),dir=%d:Unode=%d(N%d%c), Vnode=%d(N%d%c), Vnode->SPcost[dir][Vflipped]=%0.3f, Vnode->heapindex[Vflipped]=%d\n",
	      Rnode->nodeid,nodeorder[Rnode->mapid],flip ? '\'':' ',dir,Unode->nodeid,nodeorder[Unode->mapid],Uflipped ? '\'':' ',
	      Vnode->nodeid,nodeorder[Vnode->mapid],Vflipped ? '\'':' ', SP[Vnodeid].SPcost[dir][Vflipped],SP[Vnodeid].heapindex[Vflipped]);
	  fflush(stdout);
	  assert(heapindex >= 1);/* Vnode should still be in heap */
	}
	Cnodeflip *pheap = &heap[heapindex];
	if(DEBUG>=2) assert(pheap->nodeid == Vnodeid);
	if(DEBUG>=2) assert(pheap->flip == Vflipped);
	SP[Vnodeid].SPcost[dir][Vflipped] = cost;
	SP[Vnodeid].SPdist[dir][Vflipped] = distance;
	SP[Vnodeid].SPedge[dir][Vflipped] = UVedge;
	sift_up(heap,heapcnt,heapindex,dir,SP);
      }
    }
  } // while(heapcnt > 0)

  if(VERB>=2 && (Rnodeid==148965 || Rnodeid==1240440) && Rflipped==0 && dir==1 && metric==0 /*  && tid==0 && DFSverb*/){
    Cnode *Lnode = CCnodes ? CCnodes[ret.nodeid] : &node[ret.nodeid];
    printf("longest distance forward from node %d(N%d%c) is to node %d(N%d%c),dist=%0.3f(cost=%0.3f) %s: wtime=%0.6f\n",
	      Rnode->nodeid,nodeorder[Rnode->mapid],flip ? '\'':' ', Lnode->nodeid,nodeorder[Lnode->mapid],ret.flip ? '\'':' ',
	      SP[ret.nodeid].SPdist[dir][ret.flip], SP[ret.nodeid].SPcost[dir][ret.flip], dir ? "(path will be reversed)" : "",wtime());
    fflush(stdout);
  }

  /* reset SPmark and (if dir==1) reverse SPflip[1][0..1] for all nodes */
  if(BFSINIT &&  maxcost > 0.0){
    if(dir>0)
      SPinitBFS<-1,1,1>(Rnode,SP,maxcost);
    else
      SPinitBFS<-1,-1,0>(Rnode,SP,maxcost);
  } else {
    int nodecnt;
    if(dir>0)
      nodecnt = SPinitDFS<-1,1,1>(Rnodeid,SP,nodes);
    else
      nodecnt = SPinitDFS<-1,-1,0>(Rnodeid,SP,nodes);
    if(DEBUG>=2 && nodes > 0) assert(1 <= nodecnt && nodecnt <= nodes);
  }
  if(dir)
    ret.flip = 1 - ret.flip;

  if(!heapmem)
    delete [] heap;
  if(VERB>=2 && (Rnodeid==148965 || Rnodeid==1240440) && Rflipped==0 && dir==1 && metric==0 /*  && tid==0 && DFSverb*/){
	    printf("Completed ShortestPaths from Rnodeid=%d(nodeid=%d): wtime=%0.6f\n",Rnodeid,Rnode->nodeid,wtime());
    fflush(stdout);
  }
  return ret; 
}

/** If inc = 1, mark will be incremented from 0 to 1 and 2.
   If inc = -1, mark will be decremented from 2 to 1 and 0. */
static void displayDFS(Cedge *pedge, int inc, int ucnt)
{
  if(pedge->mark == 1){/* cycle in existing DAG detected */
    printf("Cycle in existing update graph detected\n");
    exit(1);
  }
  if(pedge->mark == inc+1/*(inc>0 ? 2 : 0)*/){/* this Edge already explored and finished */
    if(inc>0 && cdepth+1 <= 10000)
      printf("%3d:N%d to N%d:cov=%0.1f,ucnt=%d(cum=%d)(updated)\n",cdepth+1,nodeorder[node[pedge->Lmap].mapid],nodeorder[node[pedge->Rmap].mapid],pedge->cov,pedge->ucnt,ucnt*pedge->ucnt);
    return;
  }

  pedge->mark += inc;/* mark this Edge as active */
  if(DEBUG>=2) assert(pedge->mark == 1);

  int origcdepth = cdepth;
  cdepth++;
  if(inc>0 && cdepth <= 10000)
    printf("%3d:N%d to N%d:cov=%0.1f,ucnt=%d(cum=%d)\n",cdepth,nodeorder[node[pedge->Lmap].mapid],nodeorder[node[pedge->Rmap].mapid],pedge->cov,pedge->ucnt,ucnt*pedge->ucnt);

  /* recurse over all children */
  for(std::vector<int>::reverse_iterator it = pedge->update.rbegin() ; it != pedge->update.rend(); ++it)
    displayDFS(&edge[*it], inc, ucnt * pedge->ucnt);
  
  pedge->mark += inc;/* mark this Edge as finished */
  if(DEBUG>=2) assert(pedge->mark == (inc > 0 ? 2 : 0));
  cdepth--;
  
  if(DEBUG>=2) assert(cdepth == origcdepth);
}

/** display the coverage value of all other edges that will be added to the specified edges coverage */
static void covdisplay(Cedge *pedge)
{
  printf("coverage updates scheduled for edge=(%d,%d,%0.3f,%d,%d):Initial cov=%0.1f(N%d to N%d),ucnt=%d\n",
	 pedge->Lmap,pedge->Rmap,pedge->distance,pedge->Lflipped,pedge->Rflipped,pedge->cov,
	 nodeorder[node[pedge->Lmap].mapid],nodeorder[node[pedge->Rmap].mapid],pedge->ucnt);
  fflush(stdout);
  
  cdepth = 0;
  displayDFS(pedge,1,pedge->ucnt);
  if(DEBUG) assert(cdepth==0);
  displayDFS(pedge,-1,pedge->ucnt);
  if(DEBUG) assert(cdepth==0);
}

extern int verb;

/** See Algorithm 4.1.1 in Nguyen's thesis : redundancy is initialized at cov=1 rather than cov=0 in Nguyen's thesis
* Basic step is to find edges that span (are in parallel to) 2 other edges with consistent orientation and distance 
* and eliminate the single edge and add its
   coverage to the other 2 edges */

void graph_redundancy(int verb)
{
  compute_edgecnt_alloc();/* compute maxedge and allocate global VL,VR,WL,WR (if not already done) */

  /* initialization */
  /* Edge's already set to valid=1, cov >= 1, mark=0, update=0 */
  /* Node's already set to valid=1, cov >= 1, mark=0, update=0 */
  if(DEBUG>=2){
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      assert(Vnode->update.size() == 0);
      assert(Vnode->mark == 0);
      assert(Vnode->cov >= 1);
      for(int k = 0; k < Vnode->edgecnt; k++){
	Cedge *pedge = &edge[Vnode->edgeid[k]];
	assert(pedge->update.size() == 0);
	assert(pedge->mark == 0);
	assert(pedge->valid == 1);
	assert(pedge->cov >= 1);
      }
    }
  }

  if(!PosEdgeInit)
    PosEdge();

  double varF3 = 3.0*SM*SM * (MAXERR*MAXERR);
  double varSD = SD[0]*SD[0] * (MAXERR*MAXERR);

  if(DEBUG>=2){
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      assert(Vnode == &node[Vnode->nodeid]);
      for(int k = 0; k < Vnode->edgecnt; k++){
	Cedge *pedge = &edge[Vnode->edgeid[k]];
	assert(pedge->mark == 0);
	assert(pedge == &edge[pedge->edgeid]);
      }
    }
  }

  if(VERB>=2){
    printf("At Start of graph_redundancy:\n");
    dumpmemmap();
    fflush(stdout);
  }

  int progress = 1;
  for(int iter=0;progress; iter++){/* repeat until no progress */
    progress = 0; 
    int conflictcheck = 0;/* number of conflicts checked */
    int scheduled = 0;/* deletions scheduled */
    long long updates = 0;/* updates scheduled */
    int edgecnt = 0;
    
    if(VERB>=2){
      printf("graph_redundancy:starting iter=%d\n",iter);
      fflush(stdout);
    }

    /* initialize updistance = distance for all edges */
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      for(int vw= 0; vw < Vnode->edgecnt;vw++){
	Cedge *Vedge = &edge[Vnode->edgeid[vw]];
	Vedge->updistance = Vedge->distance;
      }
    }      

    if(1) {
      int myprogress = 0;
      int myconflictcheck = 0;
      int myscheduled = 0;
      long long myupdates = 0;
      
      /* allocate per-thread heap memory */
      char *Vmark = new char[numnodes];
      for(int i = numnodes; --i >= 0;)
	Vmark[i] = 0;
      char *Emark = new char[Numedges];
      for(int i = Numedges; --i >= 0;)
	Emark[i] = 0;

      /* NOTE : valid = -1 for nodes that were marked for deletion in previous iterations over "iter".
	 valid will be changed from 1 to 0 in the following loop for nodes deleted in this iteration of "iter" */

      /* Main loop over all nodes */
      for(int i=0; i < numtopnodes; i++){
	if(VERB>=2 && !(i%1000)){
	  printf("i=%d/%d\n",i,numtopnodes);
	  fflush(stdout);
	}


	Cnode *Vnode = topnodes[i];
	if(!Vnode->edgecnt)
	  continue;

	/* mark all nodes reachable from Vnode */
	for(int vw= 0; vw < Vnode->edgecnt;vw++){
	  Cedge *Vedge = &edge[Vnode->edgeid[vw]];
	  if(Vedge->valid < 0)
	    continue;
	  Vmark[node[(Vedge->Lmap == Vnode->nodeid) ? Vedge->Rmap : Vedge->Lmap].nodeid] = 1;/* Marked */
	}
      
	for(int Vflipped = 0; Vflipped <= 1; Vflipped++){
	  /* Vflipped == 0 : look for +ve neighbors of +ve neighbors that are marked nodes 
	     Vflipped == 1 : look for -ve neighbors of -ve neighbors that are marked nodes  */
	  for(int vw = Vnode->PosEdge; Vflipped ? (--vw >= 0) : (vw < Vnode->edgecnt); vw += (Vflipped ? 0 : 1)){    
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	    if(VWedge->valid < 0)
	      continue;
	    Cnode *Wnode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
	    int Wflipped = (Vflipped ^ (VWedge->Lflipped | VWedge->Rflipped));

	    /* If Wflipped==1: look for neighbors of W with -ve node distance in increasing order of distance
	       If Wflipped==0: look for neighers of W with +ve node distance in increasing order of distance */
	    for(int wx = Wnode->PosEdge; Wflipped ? (--wx >= 0) : (wx < Wnode->edgecnt); wx += (Wflipped ? 0 : 1)){
	      Cedge *WXedge = &edge[Wnode->edgeid[wx]];
	      if(WXedge->valid < 0)
		continue;

	      Cnode *Xnode = &node[(WXedge->Lmap == Wnode->nodeid) ? WXedge->Rmap : WXedge->Lmap];

	      //	      if(DEBUG) assert(Xnode->mark == Vmark[Xnode->nodeid]);
	      if(!Vmark[Xnode->nodeid])
		continue;
	      if(Xnode->nodeid < Vnode->nodeid)
		continue;/* avoid duplicate work */
	      int Xflipped = (Wflipped ^ (WXedge->Lflipped | WXedge->Rflipped));

	      /* locate edge from V to X  with correct orientation of X */
	      int vx;
	      Cedge *VXedge = 0;
	      for(vx = 0; vx < Vnode->edgecnt; vx++){
		VXedge = &edge[Vnode->edgeid[vx]];
		if(VXedge->valid < 0)
		  continue;
		if((VXedge->Lmap == Vnode->nodeid ? VXedge->Rmap : VXedge->Lmap) == Xnode->nodeid)
		  break;
	      }
	      if(DEBUG>=2 && !(vx < Vnode->edgecnt)){
		printf("V=N%d%c,W=N%d%c,X=N%d%c:Xnode=>mark=%d,Vnode->PosEdge=%d,Vnode->edgecnt=%d\n",
		       nodeorder[Vnode->mapid],Vflipped ? '\'':' ',
		       nodeorder[Wnode->mapid],Wflipped ? '\'':' ',
		       nodeorder[Xnode->mapid],Xflipped ? '\'':' ',
		       Vmark[Xnode->nodeid],Wnode->PosEdge,Wnode->edgecnt);
		for(vx = 0; vx < Vnode->edgecnt; vx++){
		  VXedge = &edge[Vnode->edgeid[vx]];
		  printf("vx=%d:(N%d,N%d,%0.3f,%d,%d,cov=%0.1f)\n",vx,nodeorder[node[VXedge->Lmap].mapid],nodeorder[node[VXedge->Rmap].mapid],VXedge->distance,VXedge->Lflipped,VXedge->Rflipped,VXedge->cov);
		}	      
		fflush(stdout);
		assert(vx < Vnode->edgecnt);
	      }
	      if(DEBUG>=2) assert(VXedge); 

	      /* check that orientation of X is correct */
	      if(Xflipped != (Vflipped ^ (VXedge->Lflipped | VXedge->Rflipped)))
		continue;

	      if(DEBUG>=2) assert(VXedge->valid >= 0 && VWedge->valid >= 0 && WXedge->valid >= 0);

	      double VXdistance = VXedge->distance;
	      if(Vflipped ? (vx >= Vnode->PosEdge) : (vx < Vnode->PosEdge))
		VXdistance = -VXdistance;

	      if(VXdistance < 0.0)/* triangle inequality : can happen for values close to 0.0 */
		continue;/* This case is best handled by graph_bulge(), to determine which edge is erroneous */

	      /* check that VXedge->distance ~= VWedge->distance + WXedge->distance */
	      double sum = VWedge->distance + WXedge->distance;
	      double var = varF3 + varSD*(sum + VXdistance);
	      double err = sum - VXdistance;
	      if(err*err >= var)/* path distances do NOT match */
		continue;	      /* This is best handled by graph_bulge() */

	      if(DEBUG/* HERE >=2 */) assert(VXedge->alignid || VXedge->cov <= 0.0);// Super-Edges not supported in graph_redundancy()

	      if(FAST>=2 && VXedge->valid <= 0)
		continue;
	    
	      double newdistance = max(VWedge->updistance,WXedge->updistance);
	      int checkconflict = (newdistance >= VXedge->updistance) ? 1 : 0;/* If 1, need to check for conflicts */
	      if(checkconflict){
		if(CFAST && iter < CFAST-1 && (CFAST>=2 || VXedge->valid <= 0)){
		  if(!myprogress)
		    myprogress = 1;
		  continue;
		}
		if(DEBUG>=2 && CFAST>=2) assert(iter >= CFAST-1);
		myconflictcheck++;

		/* (Perform DFS search of existing "update" DAG to see if there is already a path from VXedge to VWedge or WXedge) */
		if(DEBUG>=2) assert(Emark[VWedge->edgeid] == 0);
		if(DEBUG>=2) assert(Emark[WXedge->edgeid] == 0);
		if(DEBUG>=2) assert(Emark[VXedge->edgeid] == 0);
		Emark[VWedge->edgeid] = -1;
		Emark[WXedge->edgeid] = -1;
		int conflict = FindpathUpdate(VXedge,Emark);

		Emark[VWedge->edgeid] = 0;
		Emark[WXedge->edgeid] = 0;
		if(DEBUG>=2) assert(Emark[VXedge->edgeid] == 0);

		if(conflict){
		  myprogress++;/* make sure it is checked again later  */
		  continue;
		}

		/* propagate max(VWedge->updistance,WXedge->updistance) to VXedge and its update successors */
		distanceUpdate(newdistance,VXedge);
	      }

	      if(VXedge->valid > 0)
		myscheduled++;
	      myupdates += 2;

	      if(VERB>=2 /* && ((nodeorder[Vnode->mapid] == 1536 && nodeorder[Xnode->mapid] == 63) ||
			    (nodeorder[Vnode->mapid] == 63 && nodeorder[Xnode->mapid] == 1536) ||
			    (nodeorder[Vnode->mapid] == 70 && nodeorder[Xnode->mapid] == 63) ||
			    (nodeorder[Vnode->mapid] == 63 && nodeorder[Xnode->mapid] == 70)
			    )*/){
		printf("V=N%d(Vflipped=%d),W=N%d(Wflipped=%d),X=N%d(Xflipped=%d),VW=%0.3f,WX=%0.3f,VX=%0.3f:sum=%0.3f,sd=%0.3f,err=%0.3f\n",
		       nodeorder[Vnode->mapid],Vflipped,nodeorder[Wnode->mapid],Wflipped,nodeorder[Xnode->mapid],Xflipped,VWedge->distance,WXedge->distance,
		       VXdistance,sum,sqrt(var),err);
		printf(" VWedge=(N%d,N%d,%0.3f,%d,%d),cov=%0.1f:valid=%d,ucnt=%d\n",
		       nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,VWedge->Lflipped,VWedge->Rflipped,VWedge->cov,VWedge->valid,VWedge->ucnt);
		printf(" WXedge=(N%d,N%d,%0.3f,%d,%d),cov=%0.1f:valid=%d,ucnt=%d\n",
		       nodeorder[node[WXedge->Lmap].mapid],nodeorder[node[WXedge->Rmap].mapid],WXedge->distance,WXedge->Lflipped,WXedge->Rflipped,WXedge->cov,WXedge->valid,WXedge->ucnt);
		printf(" VXedge=(N%d,N%d,%0.3f,%d,%d),cov=%0.1f:valid=%d,ucnt=%d\n",
		       nodeorder[node[VXedge->Lmap].mapid],nodeorder[node[VXedge->Rmap].mapid],VXedge->distance,VXedge->Lflipped,VXedge->Rflipped,VXedge->cov,VXedge->valid,VXedge->ucnt);
		fflush(stdout);
	      }

	      if(DEBUG>=2) assert(VXedge->ucnt >= ((VXedge->valid <= 0) ? 1 : 0));

	      /* mark VX for deletion */
	      VXedge->valid = 0;
	      VXedge->ucnt++;

	      /* schedule updates of "cov" values of VWedge,WXedge and Wnode from final "cov" value of VXedge */
	      if(DEBUG>=2 && CFAST >=2 && iter < CFAST-1){
		assert(VXedge->distance > VWedge->distance);
		assert(VXedge->distance > WXedge->distance);
	      }
	      VWedge->update.push_back(VXedge->edgeid);// = new Cupdate(VXedge,VWedge->update);
	      WXedge->update.push_back(VXedge->edgeid);// = new Cupdate(VXedge,WXedge->update);
	      Wnode->update.push_back(VXedge->edgeid); // = new Cupdate(VXedge,Wnode->update);
	    } /* wx loop */
	  } /* vw loop */
	}/* Vflipped  loop */

	/* reset mark from all neighbors of V */      
	for(int vw = 0; vw < Vnode->edgecnt; vw++){    
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  Vmark[node[VWedge->Lmap == Vnode->nodeid ? VWedge->Rmap : VWedge->Lmap].nodeid] = 0;
	}
      }  /* end Main loop over all nodes */

      progress += myprogress;
      conflictcheck += myconflictcheck;
      scheduled += myscheduled;
      updates += myupdates;

      if(VERB >= 3){
	for(int i=0; i < numtopnodes; i++){
	  Cnode *Wnode = topnodes[i];
	  if(nodeorder[Wnode->mapid] != 1536)
	    continue;
	  printf("iter=%d,Wnode=N%d:Node coverage updates\n",iter,nodeorder[Wnode->mapid]);
	
	  for(std::vector<int>::reverse_iterator it = Wnode->update.rbegin() ; it != Wnode->update.rend(); ++it)
	    covdisplay(&edge[*it]);

	  fflush(stdout);
	}
	for(int i=0; i < numtopnodes; i++){
	  Cnode *Vnode = topnodes[i];
	  if(!Vnode->edgecnt)
	    continue;
	
	  for(int vw = 0; vw < Vnode->edgecnt; vw++){    
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	    Cnode *Wnode = &node[VWedge->Lmap == Vnode->nodeid ? VWedge->Rmap : VWedge->Lmap];
	    if((nodeorder[Vnode->mapid] == 1536 && nodeorder[Wnode->mapid] == 155) ||
	       (nodeorder[Vnode->mapid] == 1536 && nodeorder[Wnode->mapid] == 63)
	       ){
	      printf("iter=%d,Vnode=N%d,Wnode=N%d:Edge coverage updates:\n",iter,nodeorder[Vnode->mapid],nodeorder[Wnode->mapid]);

	      for(std::vector<int>::reverse_iterator it = VWedge->update.rbegin() ; it != VWedge->update.rend(); ++it){
		Cedge *uedge = &edge[*it];
		printf("    update from edge=(N%d,N%d)\n", nodeorder[node[uedge->Lmap].mapid],nodeorder[node[uedge->Rmap].mapid]);
	      }	  	

	      //	    covdisplay(VWedge);
	    }
	  }
	}
      }// VERB>=3

      /* free per-thread heap memory */
      if(DEBUG>=2){
	for(int i = numnodes; --i >= 0;)
	  assert(Vmark[i] == 0);
	for(int i = Numedges; --i >= 0;)
	  assert(Emark[i] == 0);
      }
      delete [] Vmark;
      delete [] Emark;
    } /* if(1) */

    if(VERB>=2){
      printf("graph_redundancy:iter=%d: parallel code completed:time=%0.6f(total=%0.6f)\n", iter, wtime(),mtime());
      fflush(stdout);
    }
    if(VERB>=2){/* count number of update,updated links */
      long long cnt = 0,cntd = 0,ncnt = 0, ncntd = 0,totedgecnt = 0;
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	ncnt += Vnode->update.size();
	ncntd += Vnode->updated.size();

	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  totedgecnt++;

	  cnt += VWedge->update.size();
	  cntd += VWedge->updated.size();
#if DEBUG>=2
	  VWedge->size = VWedge->update.size();
	  VWedge->sized = VWedge->updated.size();
#endif

#if DAG_CHECK
	  VWedge->updatemark = 0;
#endif
	}
      }
      printf("graph_redundancy:iter=%d:Edge update links=%lld(updated=%lld),Node update links=%lld(updated=%lld),total eges=%lld\n",iter,cnt/2,cntd/2,ncnt,ncntd,totedgecnt);
      fflush(stdout);
    }

    if(DEBUG) redudancy_iter = iter;

    /* perform all coverage updates over entire graph */
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      if(!Vnode->edgecnt)
	continue;

      /* schedule edge deletion */
      for(int vw = 0; vw < Vnode->edgecnt; vw++){    
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(VWedge->valid == 0){
	  VWedge->valid = -1;
	  edgecnt++;
	}
      }

      for(int Vflipped = 0; Vflipped <= 1; Vflipped++){
	/* propagate the scheduled coverage updates of +ve (-ve) neighbors of V (recursively) */
	for(int vw = Vnode->PosEdge; Vflipped ? (--vw >= 0) : (vw < Vnode->edgecnt); vw += (Vflipped ? 0 : 1)){    
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(DEBUG>=2){
	    assert(cdepth == 0);
	    cov_root = VWedge;
	  }
	  cov_update(VWedge);

	  Cnode *Wnode = &node[VWedge->Lmap == Vnode->nodeid ? VWedge->Rmap : VWedge->Lmap];
	  int Wflipped = (Vflipped ^ (VWedge->Lflipped | VWedge->Rflipped));
	  for(int wx = (Wflipped ? 0 : Wnode->PosEdge); wx < (Wflipped ? Wnode->PosEdge : Wnode->edgecnt); wx++){
	    Cedge *WXedge = &edge[Wnode->edgeid[wx]];
	    if(DEBUG>=2){
	      assert(cdepth == 0);
	      cov_root = WXedge;
	    }
	    cov_update(WXedge);
	  }

	  for(std::vector<int>::reverse_iterator it = Wnode->update.rbegin() ; it != Wnode->update.rend(); ++it){
	    int uedgeid = *it;
	    Cedge *uedge = &edge[uedgeid];
	    if(DEBUG>=2){
	      assert(cdepth == 0);
	      cov_root = uedge;
	    }
	    cov_update(uedge);
	    Wnode->cov += uedge->cov / uedge->ucnt;
	    if(CROSSBREAK)
	      Wnode->updated.push_back(uedgeid);/* append Wnode->update elements to Wnode->updated[] */
	  }
	  std::vector<int>().swap(Wnode->update);/* empty Wnode->update[] */
	}
      }
    }

    if(DEBUG>=2){
      for(int i=0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	if(DEBUG>=2) assert(Vnode->update.size() == 0);
	if(DEBUG>=2) assert(Vnode->mark == 0);
	for(int vw = 0; vw < Vnode->edgecnt;vw++){
	  Cedge *Vedge = &edge[Vnode->edgeid[vw]];
	  if(DEBUG>=2 && !(Vedge->update.size() == 0)){
	    printf("After graph_redundancy(before GCedges):Vnode=N%d,Vedge=(N%d,N%d,%0.3f,%d,%d):update.size()= %lld\n",
		   nodeorder[Vnode->mapid],nodeorder[node[Vedge->Lmap].mapid],nodeorder[node[Vedge->Rmap].mapid],Vedge->distance,Vedge->Lflipped,Vedge->Rflipped,(long long)Vedge->update.size());
	    fflush(stdout);
	    assert(Vedge->update.size() == 0);
	  }
	  if(DEBUG>=2) assert(Vedge->mark == 0);
	}
      }
    }

    if(VERB /* edgecnt > 0*/){
      printf("graph_redundancy:iter=%d: reduced number of edges by %d,conflicts=%d/%d,deletions=%d, updates=%lld:time=%0.6f(total=%0.6f)\n",
	     iter,edgecnt, progress,conflictcheck, scheduled, updates, wtime(),mtime());
      fflush(stdout);
    }

    if(VERB>=2){/* count number of update,updated links */
      long long cnt = 0,cntd = 0,ncnt = 0, ncntd = 0, totedgecnt = 0;
#if DEBUG>=2
      long long origcnt = 0, origcntd = 0;
#endif
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	ncnt += Vnode->update.size();
	ncntd += Vnode->updated.size();

	for (int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  totedgecnt++;

	  cnt += VWedge->update.size();
	  cntd += VWedge->updated.size();
#if DEBUG>=2
	  if (CROSSBREAK){
	    assert(VWedge->update.size() == 0);
	    assert(VWedge->updated.size() == VWedge->size + VWedge->sized);
	  }
	  origcnt += VWedge->size;
	  origcntd += VWedge->sized;
#endif
	}
      }
#if DEBUG>=2
      printf("graph_redundancy:iter=%d:Edge update links=%lld->%lld(updated=%lld->%lld),Node update links=%lld(updated=%lld),total edges=%lld\n",iter,origcnt/2,cnt/2,origcntd/2,cntd/2,ncnt,ncntd,totedgecnt);
#else
      printf("graph_redundancy:iter=%d:Edge update links=%lld(updated=%lld),Node update links=%lld(updated=%lld),total edges=%lld\n",iter,cnt/2,cntd/2,ncnt,ncntd,totedgecnt);
#endif
      fflush(stdout);
    }
  }/* while loop */

  if(VERB>=2){
    dumpmemmap();
    fflush(stdout);
  }

  if(GraphWeight > 0.0){/* make sure all cov values are <= HUGE_COV */
    double maxCov = 0.0;
    Cedge *maxEdge = 0;
    for(int t=0; t < numtopnodes; t++){
      Cnode *Tnode = topnodes[t];
      for(int u = 0; u < Tnode->edgecnt; u++){
	Cedge *TUedge = &edge[Tnode->edgeid[u]];
	if(DEBUG) assert(TUedge->cov >= 0.0);
	if(TUedge->cov > maxCov){
	  maxCov = TUedge->cov;
	  maxEdge = TUedge;
	}
	if(TUedge->cov > HUGE_COV)
	  TUedge->cov = HUGE_COV;
      }	
    }
    if(VERB){
      printf("graph_redundancy: GraphWeight=%0.1f:Largest coverage = %0.2f : (N%d,N%d,%0.3f,%d,%d)\n",
	     GraphWeight,maxCov, nodeorder[node[maxEdge->Lmap].mapid],nodeorder[node[maxEdge->Rmap].mapid],maxEdge->distance,maxEdge->Lflipped,maxEdge->Rflipped);
      fflush(stdout);
    }
  }

  if(CROSSBREAK){
    int bcntmax = 0, fcntmax = 0;
    for(int i = 0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      if(Vnode->PosEdge > bcntmax)
	bcntmax = Vnode->PosEdge;
      if(Vnode->edgecnt - Vnode->PosEdge > fcntmax)
	fcntmax = Vnode->edgecnt - Vnode->PosEdge;
    }
    int *bconnect = new int[bcntmax];
    int *fconnect = new int[fcntmax];

    if(DEBUG>=2){
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  assert(VWedge->valid != 0);
	}
      }
    }

    int strongcnt = 0;
    int undocnt = 0;
    if(FORKBREAK){
      /* Before breaking crosses directly, try to undo edge merges at forked nodes :
	 If X is a node with more than 1 incoming edge but only one outgoing edge (or vise-versa), check for each pair of incomming node A and outgoing node B
	 if a direct edge A to B existed (valid = -1 or 1) with valid length and restore it with cov = min(cov(AX),cov(XB)) and reduce the coverage of AX and XB by the same amount.

	 This is a risk free operation since graph_bulge is able to undo it at a later time, but will help with break crosses when two chimeric nodes are
	 connected to each other.
      */

      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	if(Vnode->PosEdge <= 1 && Vnode->edgecnt <= Vnode->PosEdge+1)
	  continue;/* only one positive or negative edge */
	/* locate dominant forward and backward edge */
	double fcov = 0.0, bcov = 0.0;
	for(int vw = 0; vw < Vnode->PosEdge; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->valid <= 0)
	    continue;
	  if(VWedge->cov > bcov)
	    bcov = VWedge->cov;
	}
	if(bcov < CrossbreakMincov)
	  continue;

	for(int vw = Vnode->PosEdge; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->valid <= 0)
	    continue;
	  if(VWedge->cov > fcov)
	    fcov = VWedge->cov;
	}
	if(fcov < CrossbreakMincov)
	  continue;

	/* count number of strong backward and forward edges */
	int fcnt=0,bcnt=0;
	for(int vx = 0; vx < Vnode->PosEdge; vx++){
	  Cedge *VXedge = &edge[Vnode->edgeid[vx]];
	  if(VXedge->valid <= 0)
	    continue;
	  if(VXedge->cov < CrossbreakMincov)
	    continue;
	  bconnect[bcnt++] = vx;
	}
	for(int vw = Vnode->PosEdge; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];	
	  if(VWedge->valid <= 0)
	    continue;
	  if(VWedge->cov < CrossbreakMincov)
	    continue;
	  fconnect[fcnt++] = vw;
	}
	if(bcnt <=1 && fcnt <= 1)
	  continue;
	if(bcnt > 1 && fcnt > 1)
	  continue;

	double F=0.0,B=0.0;
	for(int x = 0; x < bcnt; x++){
	  Cedge *VXedge = &edge[Vnode->edgeid[bconnect[x]]];
	  if(VXedge->valid > 0)
	    B += VXedge->cov;
	}
	for(int w= 0; w < fcnt; w++){
	  Cedge *VWedge = &edge[Vnode->edgeid[fconnect[w]]];
	  if(VWedge->valid > 0)
	    F += VWedge->cov;
	}
	if(DEBUG>=2) assert(F > 0.0 && B > 0.0);
	double Bratio = min(1.0,F/B);
	double Fratio = min(1.0,B/F);

	/* locate pairs of strong backward/forward edges that have strong direct connections that can be restored */
	for(int x = 0; x < bcnt; x++){  /* scan over strong backward edges */
	  int vx = bconnect[x];
	  Cedge *VXedge = &edge[Vnode->edgeid[vx]];
	  Cnode *Xnode = &node[(VXedge->Lmap == Vnode->nodeid) ? VXedge->Rmap : VXedge->Lmap];
	  //	  int Xflipped = VXedge->Lflipped | VXedge->Rflipped;
	  if(DEBUG>=2) assert(VXedge->valid == 1);

	  for(int w = 0; w < fcnt; w++){	/* scan over strong forward edges */
	    int vw = fconnect[w];
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];	
	    Cnode *Wnode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
	    //	    int Wflipped = VWedge->Lflipped | VWedge->Rflipped;
	    if(DEBUG>=2) assert(VWedge->valid == 1);

	    /* check if there is a direct edge from Xnode to Wnode */
	    int x = 0;
	    Cedge *XTedge = 0;
	    for(; x < Xnode->edgecnt; x++){
	      XTedge = &edge[Xnode->edgeid[x]];
	      Cnode *Tnode = &node[(XTedge->Lmap == Xnode->nodeid) ? XTedge->Rmap : XTedge->Lmap]; 
	      if(Tnode == Wnode)
		break;
	    }
	    if(x < Xnode->edgecnt){/* found direct edge : now check if the edge was signficant and has the correct distance */
	      if(DEBUG>=2) assert(XTedge != 0);
	      double x = VXedge->distance + VWedge->distance;
	      double y = XTedge->distance;
	      double var = 2.0*SM*SM + SD[0]*SD[0]*(x+y);
	      double err = x - y;
	      if(err * err < var*(MAXERR*MAXERR) /* && XTedge->Ccnt >= CrossbreakMincnt */){/* restore/strengthen this edge */
		strongcnt++;

		/* make sure side with single edge does not have coverage reduced by more than pro-rata amount (Bratio or Fratio) */
		double covDelta = min(VXedge->cov*Bratio,VWedge->cov*Fratio) - 1.0;
		if(XTedge->valid == -1){
		  undocnt++;
		  XTedge->valid = 0;/* avoid propagation of this edge's coverage, which would be order dependent */
		  XTedge->cov = max(1.0,covDelta);
		} else {
		  XTedge->cov += covDelta;
		}
		VXedge->cov -= covDelta;
		VWedge->cov -= covDelta;
	      }
	    }
	  }
	}
      }

      /* change any edge valid==0 to 1 */
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->valid == 0)
	    VWedge->valid = 1;
	  if(DEBUG>=2 && VWedge->valid > 0 && !(VWedge->cov >= 0.99)){
	    printf("VWedge=(N%d,N%d,%0.3f,%d,%d):cov=%0.3f,valid=%d\n",
		   nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,VWedge->Lflipped,VWedge->Rflipped,VWedge->cov,VWedge->valid);
	    fflush(stdout);
	    assert(VWedge->cov >= 0.99);
	  }
	}
      }

      if(VERB && undocnt > 0){
	printf("graph_redundancy: restored %d cross edges and strengthended %d additional edges\n",undocnt, strongcnt-undocnt);
	fflush(stdout);
      }
    }

    int bypasscnt= 0;

    #pragma omp parallel for num_threads(numthreads) schedule(dynamic,3) if(numthreads > 1)
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      int V = Vnode->nodeid;
      compute_neighbors(Vnode,VL[V],VR[V],VLcnt[V],VRcnt[V]);
    }

    if(VERB/* HERE >=2 */){
      printf("Finished compute_neighbors() for all %d nodes:time=%0.6f(total=%0.6f)\n", numtopnodes,wtime(),mtime());
      fflush(stdout);
    }

    for(int i = 0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      if(Vnode->PosEdge <= 1 && Vnode->edgecnt <= Vnode->PosEdge+1)
	continue;/* only one positive or negative edge */
      
      /* locate dominant forward and backward edge */
      double fcov = 0.0, bcov = 0.0;
      for(int vw = 0; vw < Vnode->PosEdge; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(VWedge->valid < 0)
	  continue;
	if(VWedge->cov > bcov)
	  bcov = VWedge->cov;
      }
      if(bcov < CrossbreakMincov)
	continue;
      for(int vw = Vnode->PosEdge; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];	
	if(VWedge->valid < 0)
	  continue;
	if(VWedge->cov > fcov)
	  fcov = VWedge->cov;
      }
      if(fcov < CrossbreakMincov)
	continue;

      // NOTE : try to avoid thresholding with EdgeCoverageRatio : only needed in original implementation that did not clone/split chimeric nodes 
      fcov *= EdgeCoverageRatio;
      bcov *= EdgeCoverageRatio;

      /* count number of strong backward and forward edges */
      int fcnt=0,bcnt=0;
      for(int vx = 0; vx < Vnode->PosEdge; vx++){
	Cedge *VXedge = &edge[Vnode->edgeid[vx]];
	if(VXedge->valid < 0)
	  continue;
	if(VXedge->cov < CrossbreakMincov)
	  continue;
	bcnt++;
      }
      for(int vw = Vnode->PosEdge; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];	
	if(VWedge->valid < 0)
	  continue;
	if(VWedge->cov < CrossbreakMincov)
	  continue;
	fcnt++;
      }
      if(bcnt <=1 || fcnt <= 1)
	continue;

      if(DEBUG>=2) assert(bcnt <= bcntmax);
      if(DEBUG>=2) assert(fcnt <= fcntmax);

      Vnode->bcnt = bcnt;
      Vnode->fcnt = fcnt;
      Vnode->pairs = new Ccrossbreak*[bcnt];
      Vnode->pairsmem = new Ccrossbreak[bcnt * fcnt];
      for(int i = 0; i < bcnt; i++)
	Vnode->pairs[i] = &Vnode->pairsmem[i*fcnt];

      /* locate pairs of strong backward/forward edges that have strong direct connections or shared updates */
      /* Try using min(Lcnt,Rcnt) >= CrossbreakMincnt to pick strong connections 
	 Test case : node N247 with strong connections on left to (N244,N245) N1302 and the right to N249,N1310. There is no direct edge between N1302 and either N249 or N1310.
	 Try adding ratio based test as in edgetrim : min(1,Ccnt) >= min(Lcntmax,Rcntmax)*EdgeTrimRatio for both nodes (see N858 for example)
      */

      int backward = -1;
      for(int vx = 0; vx < Vnode->PosEdge; vx++){  /* scan over strong backward edges */
	Cedge *VXedge = &edge[Vnode->edgeid[vx]];
	if(VXedge->valid < 0)
	  continue;
	if(VXedge->cov < CrossbreakMincov)
	  continue;
	backward++;
	if(DEBUG>=2) assert(backward < bcnt);

	Cnode *Xnode = &node[(VXedge->Lmap == Vnode->nodeid) ? VXedge->Rmap : VXedge->Lmap];
	int Xid = Xnode->nodeid;
	int Xflipped = VXedge->Lflipped | VXedge->Rflipped;

	int forward = -1;
	for(int vw = Vnode->PosEdge; vw < Vnode->edgecnt; vw++){	/* scan over strong forward edges */
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];	
	  if(VWedge->valid < 0)
	    continue;
	  if(VWedge->cov < CrossbreakMincov)
	    continue;
	  forward++;
	  if(DEBUG>=2) assert(forward < fcnt);

	  Cnode *Wnode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
	  int Wid = Wnode->nodeid;
	  int Wflipped = VWedge->Lflipped | VWedge->Rflipped;

	  Ccrossbreak *p = &Vnode->pairs[backward][forward];
	  /* compute the edge counts Lcnt,Rcnt,Bcnt for the (virtual) edge from Xnode to Wnode */
	  int Lcnt,Rcnt,Bcnt,Ccnt;
	  double XWdistance = VXedge->distance + VWedge->distance;

	  compute_edgecnt(Xnode,Xflipped,Wnode,Wflipped,XWdistance,Lcnt,Rcnt,Bcnt,Ccnt,EDGEBRIDGE?1:0,verb, VL[Xid],VR[Xid],VL[Wid],VR[Wid],VLcnt[Xid],VRcnt[Xid],VLcnt[Wid],VRcnt[Wid]);

	  if(VXedge->cov < bcov || VWedge->cov < fcov){
	    p->connected = 0;
	    if(VERB>=2){
	      printf("V=N%d,back edge to X=N%d(Cmax=%d,cov=%0.1f,dist=%0.3f), forward edge to W=N%d(Cmax=%d,cov=%0.1f,dist=%0.3f) (thresholded by EdgeCoverageRatio=%0.4f,bcov=%0.1f,fcov=%0.1f),Lcnt=%d,Rcnt=%d,Bcnt=%d\n",
		     nodeorder[Vnode->mapid],nodeorder[Xnode->mapid],min(Xnode->Lcntmax,Xnode->Rcntmax),
		     VXedge->cov,VXedge->distance,nodeorder[Wnode->mapid],min(Wnode->Lcntmax,Wnode->Rcntmax),
		     VWedge->cov,VWedge->distance,EdgeCoverageRatio,bcov,fcov,Lcnt,Rcnt,Bcnt);
	      fflush(stdout);
	    }
	    continue;
	  }

	  /* check if there is a strong direct edge from Xnode to Wnode */
	  int x = 0;
	  Cedge *XTedge = 0;
	  for(; x < Xnode->edgecnt; x++){
	    XTedge = &edge[Xnode->edgeid[x]];
	    Cnode *Tnode = &node[(XTedge->Lmap == Xnode->nodeid) ? XTedge->Rmap : XTedge->Lmap]; 
	    if(Tnode == Wnode)
	      break;
	  }
	  /* compute shared cov fraction between VXedge and VWedge based on updated links */
	  double bsum=0.0,fsum=0.0,shared=0.0;
	  for(std::vector<int>::reverse_iterator it = VXedge->updated.rbegin() ; it != VXedge->updated.rend(); ++it){
	    Cedge *bedge = &edge[*it];
	    double delta = bedge->cov / bedge->ucnt;
	    bsum += delta;
	    /* check if same edge (bedge) was used to update VWedge->cov */
	    for(std::vector<int>::reverse_iterator it = VWedge->updated.rbegin() ; it != VWedge->updated.rend(); ++it){
	      if(&edge[*it] == bedge){
		shared += delta;
		break;
	      }
	    }
	  }
	  for(std::vector<int>::reverse_iterator it = VWedge->updated.rbegin() ; it != VWedge->updated.rend(); ++it){
	    Cedge *fedge = &edge[*it];
	    fsum += fedge->cov / fedge->ucnt;
	  }
	  
	  if(VERB>=2){
	    if(x >= Xnode->edgecnt)
	      printf("V=N%d,back edge to X=N%d(Cmax=%d,cov=%0.1f,dist=%0.3f), forward edge to W=N%d(Cmax=%d,cov=%0.1f,dist=%0.3f): no direct edge! (bcov=%0.1f,fcov=%0.1f,shared=%0.1f),Lcnt=%d,Rcnt=%d,Bcnt=%d,Ccnt=%d\n",
		     nodeorder[Vnode->mapid],nodeorder[Xnode->mapid],min(Xnode->Lcntmax,Xnode->Rcntmax),
		     VXedge->cov,VXedge->distance,nodeorder[Wnode->mapid],min(Wnode->Lcntmax,Wnode->Rcntmax),
		     VWedge->cov,VWedge->distance,bsum,fsum,shared,Lcnt,Rcnt,Bcnt,Ccnt);
	    else {
	      printf("V=N%d,back edge to X=N%d(Cmax=%d,cov=%0.1f,dist=%0.3f), forward edge to W=N%d(Cmax=%d,cov=%0.1f,dist=%0.3f): direct edge=(N%d,N%d,%0.3f,%d,%d) cov=%0.1f(ucnt=%d),valid=%d (bcov=%0.1f,fcov=%0.1f,shared=%0.1f),Lcnt=%d,Rcnt=%d,Bcnt=%d,Ccnt=%d\n",
		     nodeorder[Vnode->mapid],nodeorder[Xnode->mapid],min(Xnode->Lcntmax,Xnode->Rcntmax),
		     VXedge->cov,VXedge->distance,nodeorder[Wnode->mapid],min(Wnode->Lcntmax,Wnode->Rcntmax),
		     VWedge->cov, VWedge->distance, nodeorder[node[XTedge->Lmap].mapid],nodeorder[node[XTedge->Rmap].mapid],XTedge->distance,XTedge->Lflipped,XTedge->Rflipped,
		     XTedge->cov, XTedge->ucnt, XTedge->valid,bsum,fsum,shared,Lcnt,Rcnt,Bcnt,Ccnt);
	    }
	    fflush(stdout);
	  }

	  if(x < Xnode->edgecnt){
	    if(DEBUG>=2) assert(XTedge != 0);
	    if(XTedge->cov > shared)
	      shared = XTedge->cov;
	  } else
	    XTedge = 0;
	  if(0){/* old method */
	    double ratio = shared/min(bsum,fsum);
	    p->connected = (ratio >= CrossbreakMinratio) ? 1 : 0;
	  } else {/* new method */
	    p->connected = 0;
	    if(Ccnt >= CrossbreakMincnt)
	      p->connected = 1;
	    else if(max(1,Ccnt) >= min(Xnode->Lcntmax,Xnode->Rcntmax)*EdgeTrimRatio &&
		    max(1,Ccnt) >= min(Wnode->Lcntmax,Wnode->Rcntmax)*EdgeTrimRatio)
	      p->connected = 1;
	  }
	  p->directEdge = XTedge;
	  p->backEdge = VXedge;
	  p->backNode = Xnode;
	  p->foreEdge = VWedge;
	  p->foreNode = Wnode;
	}
	if(DEBUG>=2) assert(forward == fcnt-1);
      }
      if(DEBUG>=2) assert(backward == bcnt-1);
    }

    if(VERB/* HERE >=2 */){
      printf("Finished locating strong backward/foward edges with strong direct connections: time=%0.6f(total=%0.6f)\n",
	     wtime(),mtime());
      fflush(stdout);
    }

    for(int i = 0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      if(!Vnode->pairs)
	continue;
      int bcnt = Vnode->bcnt;
      int fcnt = Vnode->fcnt;
      /* check if a subset of the backward edges is strongly connected to a subset of the forward edges : If so seperate them from Vnode and connect them directly */

      for(int b = 0; b < bcnt; b++){
	Ccrossbreak *p = Vnode->pairs[b];
	int fconnect_cnt = 0;
	for(int f = 0; f < fcnt; f++, p++){
	  if(p->connected){
	    if(p->backEdge->valid < 0 || p->foreEdge->valid < 0)/* these edges have been deleted due to a neighboring chimeric node */
	      p->connected = 0;
	    else
	      fconnect[fconnect_cnt++] = f;
	  }
	}
	if(!fconnect_cnt)
	  continue;
	if(fconnect_cnt >= fcnt)
	  break;/* all forward edges are connected to backward edge b */
	
	for(int bb = 0; bb < bcnt; bb++)
	  bconnect[bb] = 0;
	bconnect[b] = 1;
	/* check if fconnect[0..fconnect_cnt-1] are connected to any edges other than b */
	for(int fi = 0; fi < fconnect_cnt; fi++){
	  int f = fconnect[fi];
	  for(int bb = 0; bb < bcnt; bb++){
	    if(bb == b)
	      continue;
	    p = &Vnode->pairs[bb][f];
	    if(p->connected){
	      if(p->backEdge->valid < 0 || p->foreEdge->valid < 0)/* these edges have been deleted due to a neighboring chimeric node */
		p->connected = 0;
	      else
		bconnect[bb] = 1;
	    }
	  }
	}
	
	/* count how many back edges were connected */
	int bconnect_cnt = 0;
	for(int bb = 0;bb < bcnt; bb++)
	  if(bconnect[bb])
	    bconnect[bconnect_cnt++] = bb;
	if(DEBUG>=2) assert(bconnect_cnt >= 1);/* must at least include edge b */
	if(bconnect_cnt >= bcnt)
	  break;/* all backward edges are connected to forward edges that are connected to backward edge b */

	/* NOTE : use node cloning instead of using existing direct edges only : this is more flexible and allows each chimeric clones trimL,trimR to be re-optimized */

	if(bconnect_cnt > 1 || fconnect_cnt > 1)
	  continue;/* This case is more difficult to handle  : most chimeric nodes can be handled without it */

	/* break out forward edges fconnect[0..fconnect_cnt-1] and backward edges bconnect[0..bconnect_cnt-1] and replace them with direct edges */
	for(int bi = 0; bi < bconnect_cnt; bi++){
	  int bb = bconnect[bi];
	  for(int fi = 0; fi < fconnect_cnt; fi++){
	    int ff = fconnect[fi];
	    Ccrossbreak *p = &Vnode->pairs[bb][ff];
	    if(p->connected && p->directEdge){
	      p->connected = 0;/* so we don't find this set of connections again */
	      Cedge *nedge = p->directEdge;
	      if(nedge->valid==1)
		continue;/* presumably this edge did not match the distance of the foreEdge + backEdge */

	      bypasscnt++;

	      Cedge *backEdge = p->backEdge;
	      Cedge *foreEdge = p->foreEdge;
	      Cnode *Bnode = p->backNode;
	      int Bflipped = (backEdge->Lflipped | backEdge->Rflipped);
	      Cnode *Fnode = p->foreNode;
	      int Fflipped = (foreEdge->Lflipped | foreEdge->Rflipped);

	      nedge->valid = 1;
	      nedge->ucnt = 0;
	      nedge->cov = min(backEdge->cov,foreEdge->cov);

	      if(VERB>=2) {
		printf("Replacing N%d%c to N%d to N%d%c with direct connection bypassing middle node:(N%d,N%d,%0.3f,%d,%d),cov=%0.1f (valid=%d,%d)\n",
		       nodeorder[Bnode->mapid],Bflipped ? '\'':' ',
		       nodeorder[Vnode->mapid],
		       nodeorder[Fnode->mapid],Fflipped ? '\'':' ',
		       nodeorder[node[nedge->Lmap].mapid],nodeorder[node[nedge->Rmap].mapid],nedge->distance,nedge->Lflipped,nedge->Rflipped,nedge->cov,
		       backEdge->valid,foreEdge->valid);
		fflush(stdout);
	      }
	      if(DEBUG>=2) assert(backEdge->valid > 0);
	      if(DEBUG>=2) assert(foreEdge->valid > 0);
	      backEdge->valid = -1;
	      foreEdge->valid = -1;
	    }
	  } /* fi loop */
	} /* bi loop */
      } /* b loop */
    }/* Vnode loop */

    delete [] bconnect;
    delete [] fconnect;

    if(VERB && bypasscnt > 0){
      printf("graph_redundancy: bypassed %d crossed nodes:time=%0.6f(total=%0.6f)\n", bypasscnt,wtime(),mtime());
      fflush(stdout);
    }

    if(VERB>=2){/* count number of update,updated links */
      long long cnt = 0,cntd = 0,ncnt = 0, ncntd = 0;
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	ncnt += Vnode->update.size();
	ncntd += Vnode->updated.size();

	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  cnt += VWedge->update.size();
	  cntd += VWedge->updated.size();
	}
      }
      printf("graph_redundancy:Before freeing updated:Edge update links=%lld(updated=%lld),Node update links=%lld(updated=%lld):time=%0.6f\n",cnt/2,cntd/2,ncnt,ncntd, wtime());
      fflush(stdout);
    }

    /* free up updated, pairs, pairsmem */
    for(int i = 0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	VWedge->update_free();
      }
      Vnode->update_free();
    }

    if(VERB/* HERE >=2 */){
      printf("Finshed freeing update links:time=%0.6f(total=%0.6f)\n",wtime(),mtime());
      fflush(stdout);
    }

    if(VERB>=2){/* count number of update,updated links */
      long long cnt = 0,cntd = 0,ncnt = 0, ncntd = 0;
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	ncnt += Vnode->update.size();
	ncntd += Vnode->updated.size();
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  cnt += VWedge->update.size();
	  cntd += VWedge->updated.size();
	}
      }
      printf("graph_redundancy:After freeing updated:Edge update links=%lld(updated=%lld),Node update links=%lld(updated=%lld)\n",cnt/2,cntd/2,ncnt,ncntd);
      fflush(stdout);
    }
  }

  /* compact Gmap[],node[] and edge[] arrays to include only surviving maps/nodes & edges (at this point there are no compound edges but 2 nodes may have the same mapid */


  int totcnt= 0;/* total surviving edge pointers */
  int highcnt = 0;/* edge pointer with cov >= MinCoverage */
  int endcnt= 0;/* nodes with all +ve or all -ve edges */
  int chaincnt = 0;/* nodes with exactly 1 -ve and 1 -ve edge */

  unsigned int finalcnt;
  long long delta = GCedges(finalcnt,edge);    /* remove edges with valid <= -1 from edge lists of all nodes */
  PosEdge();    /* reconstruct PosEdges and nodedistances */
    
  if(DEBUG>=2 || delta > 0){
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      if(DEBUG>=2) assert(Vnode->update.size() == 0);
      if(DEBUG>=2) assert(Vnode->mark == 0);
      for(int vw = 0; vw < Vnode->edgecnt;vw++){
	Cedge *Vedge = &edge[Vnode->edgeid[vw]];
	if(DEBUG>=2 && !(Vedge->update.size() == 0)){
	  printf("After graph_redundancy:Vnode=%d,Vedge=(%d,%d,%0.3f,%d,%d):update.size() = %lld\n",
		 Vnode->mapid,Vedge->Lmap,Vedge->Rmap,Vedge->distance,Vedge->Lflipped,Vedge->Rflipped,(long long)Vedge->update.size());
	  fflush(stdout);
	  assert(Vedge->update.size() == 0);
	}
	if(DEBUG>=2) assert(Vedge->valid == 1);
	if(DEBUG>=2) assert(Vedge->mark == 0);
	if(Vedge->cov >= MinCoverage)
	  highcnt++;
      }
      totcnt += Vnode->edgecnt;

      if(Vnode->edgecnt && (Vnode->PosEdge <= 0 || Vnode->PosEdge >= Vnode->edgecnt))
	endcnt++;
      if(Vnode->edgecnt==2 && Vnode->PosEdge == 1)
	chaincnt++;
    }
  }
  if(VERB && delta > 0){
    printf("graph_redundancy: reduced number of edges from %lld to %d(%d with cov >= %d), %d end nodes, %d chain nodes(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
	   (delta+totcnt)/2,totcnt/2,highcnt/2,MinCoverage,endcnt,chaincnt,finalcnt,numtopnodes,wtime(),mtime());
    fflush(stdout);
  }

  if(VL){
    delete [] VLmem;
    delete [] VLcnt;
    delete [] VL;
    VL = 0;
  }
  /*
  if(VERB){
    printf("Calling malloc_trim(0): wtime= %0.6f\n",wtime());
    fflush(stdout);
  }
  malloc_trim(0);
  */
}

static int chaincycle_warning = 0;

/** Represent linear chains as single Super-Edges (with pointers to the start/end of the original chain)
   (Chain Collapsal in Nygen's thesis) Coverage of Super-Edges is based on adding distance/coverage ratio and distance values.

   If one edge in the chain has coverage >= MaxCoverage, all edges must have coverage >= MaxCoverage (node coverage must be >= MaxCoverage/2)

   Otherwise if one edge lies between 2^n and below 2^(n+1) all edges must lie in the same range (node coverage must be between 2^(n-1) and 2^(n+1))

   Also the end nodes cannot be chimeric (low node coverage), even though they are not part of the chain.

   Returns number of Edges collapsed (0 if nothing was changed)
*/

int graph_collapse(int MinCoverage, int verb)
{
  if(DEBUG>=2) assert(topnodes);

  if((Numedges + 1024LL) * 5 >= Maxedges * 4LL)/* allocate more edges */
    maxedgealloc(min(MASK_LL(32), (Numedges + 1024LL) * 5LL/4LL), Maxedges, edge);

  if(VERB>=2){
    printf("graph_collapse(): Numedges=%u, Maxedges=%u\n",Numedges,Maxedges);
    fflush(stdout);
  }

  if(!PosEdgeInit)    /* update PosEdge of remaining topnodes */
    PosEdge();

  int numchains = 0;/* number of chains found */
  int chaincnt = 0;/* Count of nodes collapsed */

  if(DEBUG>=2){
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      assert(Vnode->mark == 0);
      for(int vw= 0; vw < Vnode->edgecnt;vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	assert(VWedge->mark == 0);
	if(VWedge->valid <= 0)
	  continue;
	if(!((VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid) && (VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap)))
	   /* || (Vnode->nodeid == 2414381) || (Vnode->nodeid == 766410) || (Vnode->nodeid == 928430) || (Vnode->nodeid == 301017)*/){
          printf("i=%d/%d:nodeid=%d, vw=%d/%d: valid=%d,Lmap=%d,Rmap=%d,node[i].Enodeid[vw]= %d\n",
		 i,numnodes,Vnode->nodeid,vw,Vnode->edgecnt,VWedge->valid,VWedge->Lmap,VWedge->Rmap,Vnode->Enodeid[vw]);
	  fflush(stdout);
          assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
	  assert(VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap));
        }
      }
    }
  }

  stack<Cnode*> ChimericNodes;/* stack of pointers to nodes marked as chimeric (so they can be recognized later even if their characteristics change) */

  /* Main loop over all current nodes */
  for(int i=0; i < numtopnodes; i++){
    Cnode *Vnode = topnodes[i];
    double mincov = MaxCoverage, maxcov = 1000000.0;
    if(Vnode->mark == 2)/* already part of a collapsed chain */
      continue;
    if(Vnode->edgecnt != 2)
      continue;
    if(VERB>=2 && verb){
      printf("Vnode=N%d:PosEdge=%d,edgecnt=%d,edge[0]->distance=%0.3f,edge[1]->distance=%0.3f,Vnode->cov=%0.1f,edge[0]->cov=%0.1f,edge[1]->cov=%0.1f\n",
	     nodeorder[Vnode->mapid],Vnode->PosEdge,Vnode->edgecnt,edge[Vnode->edgeid[0]].distance,edge[Vnode->edgeid[1]].distance,Vnode->cov,edge[Vnode->edgeid[0]].cov,edge[Vnode->edgeid[1]].cov);
      fflush(stdout);
    }
    if(DEBUG>=2){
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(VWedge->valid <= 0)
	  continue;
	if(!((VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid) && (VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap)))){
          printf("Vnode->nodeid=%d(N%d), vw=%d/%d: valid=%d,edgeid=%u,Lmap=%d(N%d),Rmap=%d(N%d),Vnode->Enodeid[vw]= %d\n",
                 Vnode->nodeid,Vnode->mapid,vw,Vnode->edgeid[vw],Vnode->edgecnt,VWedge->valid,VWedge->Lmap,node[VWedge->Lmap].mapid,VWedge->Rmap,node[VWedge->Rmap].mapid,Vnode->Enodeid[vw]);
	  fflush(stdout);
          assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
	  assert(VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap));
        }
      }
    }

    if(min(edge[Vnode->edgeid[0]].cov,edge[Vnode->edgeid[1]].cov) < MaxCoverage){
      double Xmincov = min(edge[Vnode->edgeid[0]].cov,edge[Vnode->edgeid[1]].cov);
      double Xmaxcov = max(edge[Vnode->edgeid[0]].cov,edge[Vnode->edgeid[1]].cov);
      mincov = 0;
      maxcov = 1;
      while(Xmincov >= maxcov){
	mincov = maxcov;
	maxcov = 2.0 * mincov;
      }
      if(Xmaxcov > maxcov)
	continue;
    }
    if(edge[Vnode->edgeid[0]].cov < mincov || edge[Vnode->edgeid[1]].cov < mincov ||
       edge[Vnode->edgeid[0]].cov >= maxcov || edge[Vnode->edgeid[1]].cov >= maxcov)
      continue;
    if(Vnode->cov >= maxcov || Vnode->cov < mincov*0.5)
      continue;

    /* Even if Vnode->PosEdge != 1, one edge might be close to 0 */
    if(Vnode->PosEdge == 0 && edge[Vnode->edgeid[0]].distance > MAXERR*SM)
      continue;
    if(Vnode->PosEdge == 2 && edge[Vnode->edgeid[1]].distance > MAXERR*SM)
      continue;

    if(edge[Vnode->edgeid[0]].mark || edge[Vnode->edgeid[1]].mark)
      continue;/* Vnode is chimeric or next to chimeric node */

    int origchaincnt = chaincnt;

    if(DEBUG>=2) assert(Vnode->mark==0);
    Vnode->mark = 2;

    double distance = 0.0;/* sum of distance along the chain */
    double Adist = 0.0;/* sum of abs(distance) along the chain */
    double AdistBYcov = 0.0;/* sum of abs(distance)/cov along the chain */
    double cov = min(edge[Vnode->edgeid[0]].cov,edge[Vnode->edgeid[1]].cov);/* cov will be the minimum cov value along the chain, if AVCOV==0 */
    int chainlen = 1;

    /* traverse left over chain nodes */
    Cnode *Lnode = Vnode;
    int Lflipped = 0;
    Cedge *Pedge = &edge[Vnode->edgeid[1]];/* R-side edge */
    Cnode *Wnode;
    int Wflipped;
    Cedge *Ledge;

    while(1){
      int Lindex = (&edge[Lnode->edgeid[0]] == Pedge) ? 1 : 0;
      Ledge = &edge[Lnode->edgeid[Lindex]];
      double delta = (Lindex >= Lnode->PosEdge) ? -Ledge->distance : Ledge->distance;
      if(Lflipped)
	delta = -delta;
      if(DEBUG>=2) assert(fabs(delta) > 0.0);
      if(DEBUG>=2) assert(delta  >= -MAXERR*SM);

      distance += delta;
      if(AVCOV){
	Adist += fabs(delta);
	AdistBYcov += fabs(delta) / Ledge->cov;
      } else
	cov = min(cov,Ledge->cov);

      Wnode = &node[(Ledge->Lmap == Lnode->nodeid) ? Ledge->Rmap : Ledge->Lmap];
      Wflipped = (Lflipped ^ (Ledge->Lflipped | Ledge->Rflipped));
      if(VERB>=2 && verb){
	printf("Wnode=%d(N%d%c),delta=%0.3f,Lnode=%d(N%d%c),...,Vnode=%d(N%d):Wnode:edgecnt=%d,PosEdge=%d,cov=%0.1f,mark=%d:dist=%0.3f,ecov=%0.2f,Adist=%0.3f,AdistBYcov=%0.8f(scov=%0.2f)\n",
	       Wnode->nodeid,nodeorder[Wnode->mapid],Wflipped ? '\'':' ', delta,
	       Lnode->nodeid,nodeorder[Lnode->mapid],Lflipped ? '\'':' ', Vnode->nodeid,nodeorder[Vnode->mapid],
	       Wnode->edgecnt,Wnode->PosEdge,Wnode->cov,Wnode->mark,delta,Ledge->cov,Adist,AdistBYcov, Adist/AdistBYcov);
	fflush(stdout);
      }
      if(DEBUG>=2){
	for(int wl = 0; wl < Wnode->edgecnt; wl++){
	  Cedge *WLedge = &edge[Wnode->edgeid[wl]];
	  if(WLedge->valid <= 0)
	    continue;
	  if(!((WLedge->Lmap == Wnode->nodeid || WLedge->Rmap == Wnode->nodeid) && (WLedge->Lmap == Wnode->nodeid ? (Wnode->Enodeid[wl] == WLedge->Rmap) : (Wnode->Enodeid[wl] == WLedge->Lmap)))){
	    printf("Wnode->nodeid=%d(N%d), wl=%d/%d: valid=%d,edgeid=%u/%u,Lmap=%d(N%d),Rmap=%d(N%d),Wnode->Enodeid[wl]= %d\n",
		   Wnode->nodeid,Wnode->mapid,wl,Wnode->edgecnt,WLedge->valid,Wnode->edgeid[wl],Numedges,WLedge->Lmap,node[WLedge->Lmap].mapid,WLedge->Rmap,node[WLedge->Rmap].mapid,Wnode->Enodeid[wl]);
	    fflush(stdout);
	    assert(WLedge->Lmap == Wnode->nodeid || WLedge->Rmap == Wnode->nodeid);
	    assert(WLedge->Lmap == Wnode->nodeid ? (Wnode->Enodeid[wl] == WLedge->Rmap) : (Wnode->Enodeid[wl] == WLedge->Lmap));
	  }
	}
      }

      if(Wnode->edgecnt != 2)
	break;
      int Windex = (&edge[Wnode->edgeid[0]] == Ledge) ? 1 : 0;
      if(DEBUG>=2) assert(Wnode->edgeid[1-Windex] == Ledge->edgeid);
      if(DEBUG>=2) assert(Wnode->Enodeid[1-Windex] == Lnode->nodeid);
      double Wdelta = (Windex >= Wnode->PosEdge) ? -edge[Wnode->edgeid[Windex]].distance : edge[Wnode->edgeid[Windex]].distance;
      if(Wflipped)
	Wdelta = -Wdelta;
      if(Wdelta < -MAXERR*SM)
	break;
      if(Wnode->cov < mincov*0.5 || Wnode->cov >= maxcov)
      	break;
      if(edge[Wnode->edgeid[0]].cov < mincov || edge[Wnode->edgeid[1]].cov < mincov ||
	 edge[Wnode->edgeid[0]].cov >= maxcov || edge[Wnode->edgeid[1]].cov >= maxcov)
	break;
      if(Wnode->mark == 2)
	break; /* cycle */
      Cedge *Nedge = &edge[Wnode->edgeid[Windex]];
      if(Nedge->mark)/* node next to Wnode is chimeric : terminate chain at Wnode */
	break;
      if(DEBUG>=2 && Wnode->mark != 0){
	printf("Wnode->mark=%d:W=N%d(flip=%d),L=N%d(flip=%d)(possible chain cycle to left of V=N%d)\n",
	       Wnode->mark,nodeorder[Wnode->mapid],Wflipped,nodeorder[Lnode->mapid],Lflipped,nodeorder[Vnode->mapid]);
	printf("N%d%c,edge=(N%d,N%d,%0.2f,%d,%d,cov=%0.1f),N%d%c\n",nodeorder[Wnode->mapid],Wflipped ? '\'':' ',
	       nodeorder[node[Ledge->Lmap].mapid],nodeorder[node[Ledge->Rmap].mapid],Ledge->distance,Ledge->Lflipped,Ledge->Rflipped,Ledge->cov,
	       nodeorder[Lnode->mapid],Lflipped?'\'':' ');
	while(Lnode != Vnode){
	  int Lindex = (&edge[Lnode->edgeid[0]] == Ledge) ? 1 : 0;
	  Ledge = &edge[Lnode->edgeid[Lindex]];
	  double delta = (Lindex >= Lnode->PosEdge) ? Ledge->distance : -Ledge->distance;
	  if(Lflipped)
	    delta = -delta;
	  if(DEBUG) assert(fabs(delta)>0.0);
	  if(DEBUG) assert(delta  >= -MAXERR*SM);
	  Cnode *Nnode = &node[(Ledge->Lmap == Lnode->nodeid) ? Ledge->Rmap : Ledge->Lmap];
	  int Nflipped = (Lflipped ^ ( Ledge->Lflipped | Ledge->Rflipped));
	  printf("N%d%c,edge=(N%d,N%d,%0.2f,%d,%d,cov=%0.1f)N%d%c\n",nodeorder[Lnode->mapid],Lflipped?'\'':' ',
		 nodeorder[node[Ledge->Lmap].mapid],nodeorder[node[Ledge->Rmap].mapid],Ledge->distance,Ledge->Lflipped,Ledge->Rflipped,Ledge->cov,
		 nodeorder[Nnode->mapid],Nflipped ? '\'':' ');
	  Lnode = Nnode;
	  Lflipped = Nflipped;
	}
	fflush(stdout);
	assert(Wnode->mark == 0);
      }
      Wnode->mark = 2;

      Lnode = Wnode;
      Lflipped = Wflipped;
      Pedge = Ledge;
      chainlen++;
    }

    Cnode *Rnode = Vnode;
    int Rflipped = 0;
    Pedge = &edge[Vnode->edgeid[0]];/* L-side edge */
    Cnode *Xnode;
    int Xflipped;
    Cedge *Redge;
    while(1){
      int Rindex = (&edge[Rnode->edgeid[0]] == Pedge) ? 1 : 0;
      Redge = &edge[Rnode->edgeid[Rindex]];
      double delta = (Rindex >= Rnode->PosEdge) ? Redge->distance : -Redge->distance;
      if(Rflipped)
	delta = -delta;
      if(DEBUG>=2) assert(fabs(delta)>0.0);
      if(DEBUG>=2) assert(delta  >= -MAXERR*SM);

      distance += delta;
      if(AVCOV){
	Adist += fabs(delta);
	AdistBYcov += fabs(delta)/Redge->cov;
      } else
	cov = min(cov,Redge->cov);

      Xnode = &node[(Redge->Lmap == Rnode->nodeid) ? Redge->Rmap : Redge->Lmap];
      Xflipped = (Rflipped ^ (Redge->Lflipped | Redge->Rflipped));
      if(VERB>=2 && verb){
	printf("Wnode=%d(N%d%c),Lnode=%d(N%d%c),...,Vnode=%d(N%d),...,Rnode=%d(N%d%c),delta=%0.3f,Xnode=%d(N%d%c):Xnode:edgecnt=%d,PosEdge=%d,cov=%0.1f,mark=%d:dist=%0.3f,ecov=%0.2f,Adist=%0.3f,AdistBYcov=%0.8f(scov=%0.2f)\n",
	       Wnode->nodeid,nodeorder[Wnode->mapid],Wflipped ? '\'':' ',
	       Lnode->nodeid,nodeorder[Lnode->mapid],Lflipped ? '\'':' ', Vnode->nodeid,nodeorder[Vnode->mapid],
	       Rnode->nodeid,nodeorder[Rnode->mapid],Rflipped ? '\'':' ', delta,
	       Xnode->nodeid,nodeorder[Xnode->mapid],Xflipped ? '\'':' ', Xnode->edgecnt,Xnode->PosEdge,Xnode->cov,Xnode->mark, delta, Redge->cov, Adist, AdistBYcov, Adist/AdistBYcov);
	fflush(stdout);
      }
      if(DEBUG>=2){
	for(int xr = 0; xr < Xnode->edgecnt; xr++){
	  Cedge *XRedge = &edge[Xnode->edgeid[xr]];
	  if(XRedge->valid <= 0)
	    continue;
	  if(!((XRedge->Lmap == Xnode->nodeid || XRedge->Rmap == Xnode->nodeid) && (XRedge->Lmap == Xnode->nodeid ? (Xnode->Enodeid[xr] == XRedge->Rmap) : (Xnode->Enodeid[xr] == XRedge->Lmap)))){
	    printf("Xnode->nodeid=%d(N%d), xr=%d/%d: valid=%d,edgeid=%u/%u,Lmap=%d(N%d),Rmap=%d(N%d),Xnode->Enodeid[xr]= %d\n",
		   Xnode->nodeid,Xnode->mapid,xr,Xnode->edgecnt,XRedge->valid,Xnode->edgeid[xr],Numedges,XRedge->Lmap,node[XRedge->Lmap].mapid,XRedge->Rmap,node[XRedge->Rmap].mapid,Xnode->Enodeid[xr]);
	    printf("\tnumnodes=%d,numtopnodes=%d\n",numnodes,numtopnodes);

	    for(int v = 0; v < Xnode->edgecnt; v++){
	      if(v==xr)
		continue;
	      Cedge *Vedge = &edge[Xnode->edgeid[v]];
	      printf("Xnode->nodeid=%d(N%d), v=%d/%d: valid=%d,edgeid=%u/%u,Lmap=%d(N%d),Rmap=%d(N%d),Xnode->Enodeid[v]= %d\n",
		     Xnode->nodeid,Xnode->mapid,v,Xnode->edgecnt,Vedge->valid,Xnode->edgeid[v],Numedges,Vedge->Lmap,node[Vedge->Lmap].mapid,Vedge->Rmap,node[Vedge->Rmap].mapid,Xnode->Enodeid[v]);
	    }
	    fflush(stdout);

	    assert(XRedge->Lmap == Xnode->nodeid || XRedge->Rmap == Xnode->nodeid);
	    assert(XRedge->Lmap == Xnode->nodeid ? (Xnode->Enodeid[xr] == XRedge->Rmap) : (Xnode->Enodeid[xr] == XRedge->Lmap));
	  }
	}
      }

      if(Xnode->edgecnt != 2)
	break;
      int Xindex = (&edge[Xnode->edgeid[0]] == Redge) ? 1 : 0;
      if(DEBUG>=2) assert(Xnode->edgeid[1-Xindex] == Redge->edgeid);
      if(DEBUG>=2) assert(Xnode->Enodeid[1-Xindex] == Rnode->nodeid);
      double Xdelta = (Xindex >= Xnode->PosEdge) ? 
	edge[Xnode->edgeid[Xindex]].distance : -edge[Xnode->edgeid[Xindex]].distance;
      if(Xflipped)
	Xdelta = -Xdelta;
      if(Xdelta < -MAXERR*SM)
	break;
      if(Xnode->cov < mincov*0.5 || Xnode->cov >= maxcov)
	break;
      if(edge[Xnode->edgeid[0]].cov < mincov || edge[Xnode->edgeid[1]].cov < mincov ||
	 edge[Xnode->edgeid[0]].cov >= maxcov || edge[Xnode->edgeid[1]].cov >= maxcov)
	break;
      if(Xnode->mark == 2)
      	break;/* cycle */
      Cedge *Nedge = &edge[Xnode->edgeid[Xindex]];
      if(Nedge->mark)/* node next to Xnode is chimeric : terminate chain at Xnode */
	break;
      if(DEBUG>=2 && !(Xnode->mark == 0)){
	printf("Wnode=N%d(flip=%d),Lnode=N%d(flip=%d)...Vnode=N%d...Rnode=N%d(flip=%d),Xnode=N%d(flip=%d),Xnode->mark=%d\n",
	       nodeorder[Wnode->mapid],Wflipped,nodeorder[Lnode->mapid],Lflipped,nodeorder[Vnode->mapid],
	       nodeorder[Rnode->mapid],Rflipped,nodeorder[Xnode->mapid],Xflipped,Xnode->mark);
	fflush(stdout);
	assert(Xnode->mark == 0);
      }
      Xnode->mark = 2;

      Rnode = Xnode;
      Rflipped = Xflipped;
      Pedge = Redge;
      chainlen++;
    }
    
    if(VERB>=2 && verb){
      printf("Final chain: Wnode=%d(N%d%c),Lnode=%d(N%d%c),...,Vnode=%d(N%d),...,Rnode=%d(N%d%c),Xnode=%d(N%d%c):chainlen=%d\n",
	     Wnode->nodeid, nodeorder[Wnode->mapid], Wflipped ? '\'':' ',
	     Lnode->nodeid, nodeorder[Lnode->mapid], Lflipped ? '\'':' ', 
	     Vnode->nodeid, nodeorder[Vnode->mapid],
	     Rnode->nodeid, nodeorder[Rnode->mapid], Rflipped ? '\'':' ',
	     Xnode->nodeid, nodeorder[Xnode->mapid], Xflipped ? '\'':' ', chainlen);
      fflush(stdout);
    }

    /* check for cycle in chain */
    if(Xnode->mark || Wnode->mark || Xnode == Wnode){
      /* unmark chain from Lnode ... Vnode ... Rnode */
      if(VERB && !chaincycle_warning){
	printf("Warning: circular chain > N%d%c- N%d%c- ... - N%d - ... - N%d%c- N%d%c< (%d nodes,len=%0.3f,cov=%0.1f)\n",
	       nodeorder[Wnode->mapid],Wflipped ? '\'' : ' ',
	       nodeorder[Lnode->mapid],Lflipped ? '\'' : ' ',
	       nodeorder[Vnode->mapid],
	       nodeorder[Rnode->mapid],Rflipped ? '\'' : ' ',
	       nodeorder[Xnode->mapid],Xflipped ? '\'' : ' ', chainlen,distance,cov);
	fflush(stdout);
	chaincycle_warning = 1;
      }
      Pedge = Ledge;
      while(Lnode != Rnode){
	Lnode->mark = 0;
	int Lindex = (&edge[Lnode->edgeid[0]] == Pedge) ? 1 : 0;
	Cedge *Ledge = &edge[Lnode->edgeid[Lindex]];
	double delta = (Lindex >= Lnode->PosEdge) ? Ledge->distance : -Ledge->distance;
	if(Lflipped)
	  delta = -delta;
	Cnode *Nnode = &node[(Ledge->Lmap == Lnode->nodeid) ? Ledge->Rmap : Ledge->Lmap];
	int Nflipped = (Lflipped ^ ( Ledge->Lflipped | Ledge->Rflipped));
	if(DEBUG>=2) assert(fabs(delta)>0.0);

	if(DEBUG>=2 && !(delta >= -MAXERR*SM)){
	  printf("Pedge=(N%d,N%d,%0.3f,%d,%d):Lnode=N%d%c, Ledge=(N%d,N%d,%0.3f,%d,%d) Nnode=N%d%c ... Rnode=N%d%c\n",
		 nodeorder[node[Pedge->Lmap].mapid],nodeorder[node[Pedge->Rmap].mapid],Pedge->distance,Pedge->Lflipped,Pedge->Rflipped,
		 nodeorder[Lnode->mapid],Lflipped ? '\'':' ',
		 nodeorder[node[Ledge->Lmap].mapid],nodeorder[node[Ledge->Rmap].mapid],Ledge->distance,Ledge->Lflipped,Ledge->Rflipped,
		 nodeorder[Nnode->mapid],Nflipped ? '\'':' ',
		 nodeorder[Rnode->mapid],Rflipped ? '\'':' ');
	  printf("Lindex=%d,Lflipped=%d,delta=%0.3f,SM=%0.3f\n",Lindex,Lflipped,delta,SM);
	  fflush(stdout);

	  assert(delta  >= -MAXERR*SM);
	} 
	if(VERB>=2 && verb)
	  printf("Pedge=%d,%d(N%d,N%d,%0.3f,%d,%d):Lnode=N%d%c, Ledge=%d,%d,(N%d,N%d,%0.3f,%d,%d) Nnode=N%d%c ... Rnode=N%d%c\n",
                 Pedge->Lmap,Pedge->Rmap,nodeorder[node[Pedge->Lmap].mapid],nodeorder[node[Pedge->Rmap].mapid],Pedge->distance,Pedge->Lflipped,Pedge->Rflipped,
		 nodeorder[Lnode->mapid],Lflipped ? '\'':' ',
                 Ledge->Lmap,Ledge->Rmap,nodeorder[node[Ledge->Lmap].mapid],nodeorder[node[Ledge->Rmap].mapid],Ledge->distance,Ledge->Lflipped,Ledge->Rflipped,
		 nodeorder[Nnode->mapid],Nflipped ? '\'':' ',
		 nodeorder[Rnode->mapid],Rflipped ? '\'':' ');

	//	assert(Nnode != Wnode);
	Lnode = Nnode;
	Lflipped = Nflipped;
	Pedge = Ledge;
      }
      Lnode->mark = 0;
      chaincnt = origchaincnt;
      continue;
    }

    numchains++;

    /* found a linear chain >Wnode-Lnode- ... -Vnode- ... -Rnode-Xnode< */
    /* replace -Lnode-...-Vnode-...-Rnode- by a Super-Edge with computed cov, distance */
    if(Numedges+1 >= Maxedges){/* allocate more edges */
      printf("Numedges = %u, Maxedges = %u : edge[] array overflow : increase allocation at start of graph_collapse()\n",Numedges,Maxedges);
      fflush(stdout);exit(1);
    }

    Cedge *Sedge = &edge[++Numedges];// WAS32 &edge[Numedges];
    Sedge->edgeid = Numedges;// WAS32 Numedges++;
    
    Sedge->cov = AVCOV ? (Adist/AdistBYcov) : cov;
    Sedge->valid = 1;
    Sedge->alignid = 0;
    if(distance >= 0.0){
      Sedge->distance = max(0.0001,distance);
      if(Wflipped && Xflipped){/* Make Super-Edge go from Xnode to Wnode */
	Sedge->Lmap = Xnode->nodeid;
	Sedge->Rmap = Wnode->nodeid;
	Sedge->Lflipped = Sedge->Rflipped = 0;
	
	Sedge->Lchainid = Redge->edgeid;
	Sedge->Rchainid = Ledge->edgeid;
      } else {/* make Super-Edge go from Wnode to Xnode */
	Sedge->Lmap = Wnode->nodeid;
	Sedge->Lflipped = Wflipped;
	Sedge->Rmap = Xnode->nodeid;
	Sedge->Rflipped = Xflipped;
	
	Sedge->Lchainid = Ledge->edgeid;
	Sedge->Rchainid = Redge->edgeid;
      }
    } else {
      Sedge->distance = -distance;
      if(Wflipped && Xflipped){/* Make Super-Edge go from Wnode to Xnode */
	Sedge->Lmap = Wnode->nodeid;
	Sedge->Rmap = Xnode->nodeid;
	Sedge->Lflipped = Sedge->Rflipped = 0;
	
	Sedge->Lchainid = Ledge->edgeid;
	Sedge->Rchainid = Redge->edgeid;
      } else {/* Make Super-Edge go from Xnode to Wnode */
	Sedge->Lmap = Xnode->nodeid;
	Sedge->Lflipped = Xflipped;
	Sedge->Rmap = Wnode->nodeid;
	Sedge->Rflipped = Wflipped;
	
	Sedge->Lchainid = Redge->edgeid;
	Sedge->Rchainid = Ledge->edgeid;
      }
    }
    if(DEBUG>=2)assert(Sedge->distance > 0.0);

    /* In Wnode replace Ledge by Sedge */
    int k;
    for(k = 0; k < Wnode->edgecnt; k++){
      Cedge *pedge = &edge[Wnode->edgeid[k]];
      if(pedge == Ledge){
	if(DEBUG>=2) assert(pedge->valid > 0);
	Wnode->edgeid[k] = Sedge->edgeid;
	if(DEBUG>=2) assert(Sedge->Lmap == Wnode->nodeid || Sedge->Rmap == Wnode->nodeid);
	Wnode->Enodeid[k] = (Sedge->Lmap == Wnode->nodeid) ? Sedge->Rmap : Sedge->Lmap;
	break;
      }
    }
    if(DEBUG>=2) assert(k < Wnode->edgecnt);
    if(DEBUG>=2){
      for(int wl = 0; wl < Wnode->edgecnt; wl++){
	Cedge *WLedge = &edge[Wnode->edgeid[wl]];
	if(WLedge->valid <= 0)
	  continue;
	if(!((WLedge->Lmap == Wnode->nodeid || WLedge->Rmap == Wnode->nodeid) && (WLedge->Lmap == Wnode->nodeid ? (Wnode->Enodeid[wl] == WLedge->Rmap) : (Wnode->Enodeid[wl] == WLedge->Lmap)))){
	  printf("Wnode->nodeid=%d(N%d), wl=%d/%d: valid=%d,edgeid=%u/%u,Lmap=%d(N%d),Rmap=%d(N%d),Wnode->Enodeid[wl]= %d\n",
		 Wnode->nodeid,Wnode->mapid,wl,Wnode->edgecnt,WLedge->valid,Wnode->edgeid[wl],Numedges,WLedge->Lmap,node[WLedge->Lmap].mapid,WLedge->Rmap,node[WLedge->Rmap].mapid,Wnode->Enodeid[wl]);
	  fflush(stdout);
	  assert(WLedge->Lmap == Wnode->nodeid || WLedge->Rmap == Wnode->nodeid);
	  assert(WLedge->Lmap == Wnode->nodeid ? (Wnode->Enodeid[wl] == WLedge->Rmap) : (Wnode->Enodeid[wl] == WLedge->Lmap));
	}
      }
    }

    /* In Xnode replace Redge by Sedge */
    for(k = 0; k < Xnode->edgecnt; k++){
      Cedge *pedge = &edge[Xnode->edgeid[k]];
      if(pedge == Redge){
	Xnode->edgeid[k] = Sedge->edgeid;
	if(DEBUG>=2) assert(Sedge->Lmap == Xnode->nodeid || Sedge->Rmap == Xnode->nodeid);
	Xnode->Enodeid[k] = (Sedge->Lmap == Xnode->nodeid) ? Sedge->Rmap : Sedge->Lmap;
	break;
      }
    }
    if(DEBUG>=2) assert(k < Xnode->edgecnt);
    if(DEBUG>=2){
      for(int xr = 0; xr < Xnode->edgecnt; xr++){
	Cedge *XRedge = &edge[Xnode->edgeid[xr]];
	if(XRedge->valid <= 0)
	  continue;
	if(!((XRedge->Lmap == Xnode->nodeid || XRedge->Rmap == Xnode->nodeid) && (XRedge->Lmap == Xnode->nodeid ? (Xnode->Enodeid[xr] == XRedge->Rmap) : (Xnode->Enodeid[xr] == XRedge->Lmap)))){
	  printf("Xnode->nodeid=%d(N%d), xr=%d/%d: valid=%d,edgeid=%u/%u,Lmap=%d(N%d),Rmap=%d(N%d),node[i].Enodeid[xr]= %d\n",
		 Xnode->nodeid,Xnode->mapid,xr,Xnode->edgecnt,XRedge->valid,Xnode->edgeid[xr],Numedges,XRedge->Lmap,node[XRedge->Lmap].mapid,XRedge->Rmap,node[XRedge->Rmap].mapid,Xnode->Enodeid[xr]);
	  fflush(stdout);
	  assert(XRedge->Lmap == Xnode->nodeid || XRedge->Rmap == Xnode->nodeid);
	  assert(XRedge->Lmap == Xnode->nodeid ? (Xnode->Enodeid[xr] == XRedge->Rmap) : (Xnode->Enodeid[xr] == XRedge->Lmap));
	}
      }
    }

    /* nodes Lnode ... Rnode are already scheduled to be removed from topnodes (with mark=2) */

    if(VERB>=2 && verb){
      printf("Collapsed chain Wnode->nodeid=%d : N%d - N%d - ... - N%d - ... - N%d - N%d : Xnode->nodeid=%d (%d nodes/edges removed, len=%0.3f,cov=%0.1f)\n",
	     Wnode->nodeid,nodeorder[Wnode->mapid],nodeorder[Lnode->mapid],nodeorder[Vnode->mapid],nodeorder[Rnode->mapid],nodeorder[Xnode->mapid], Xnode->nodeid, chainlen, Sedge->distance, Sedge->cov);
      printf("Sedge=%d,%d(N%d,N%d,%0.3f,%d,%d)\n",Sedge->Lmap,Sedge->Rmap,nodeorder[node[Sedge->Lmap].mapid],nodeorder[node[Sedge->Rmap].mapid],Sedge->distance,Sedge->Lflipped,Sedge->Rflipped);
      fflush(stdout);
    }
    
    chaincnt += chainlen;
    
    if(!edgelist)
      edgelist = new Cnodedistance[maxedges];

    /* resort edges of Xnode and Wnode */
    for(Vnode = Wnode;;Vnode = Xnode){
      /* re-order edges for current node (Vnode->mapid) */
      edgereorder(Vnode,edgelist);
      if(Vnode == Xnode)
	break;
    }
  }

  /* unmark all edges to possibly chimeric nodes (low node coverage)*/
  while(!ChimericNodes.empty()){
    Cnode *Vnode = ChimericNodes.top(); ChimericNodes.pop();
    for(int vw= 0; vw < Vnode->edgecnt;vw++){
      Cedge *VWedge = &edge[Vnode->edgeid[vw]];
      VWedge->mark = 0;
    }
  }

  if(DEBUG>=3){
    for(int i = 0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(VWedge->valid <= 0)
	  continue;
	if(!((VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid) && (VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap)))){
          printf("i=%d/%d:nodeid=%d(N%d), vw=%d/%d: valid=%d,edgeid=%u,Lmap=%d(N%d),Rmap=%d(N%d),node[i].Enodeid[vw]= %d\n",
                 i,numnodes,Vnode->nodeid,Vnode->mapid,vw,Vnode->edgeid[vw],Vnode->edgecnt,VWedge->valid,VWedge->Lmap,node[VWedge->Lmap].mapid,VWedge->Rmap,node[VWedge->Rmap].mapid,Vnode->Enodeid[vw]);
	  fflush(stdout);
          assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
	  assert(VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap));
        }
      }
    }    
  }

  (void)GCnodes();  /* remove nodes with mark=2 from topnodes[] */
  unsigned int finalcnt;
  (void)GCedges(finalcnt,edge);

  if(VERB /* && chaincnt > 0*/){
    printf("graph_collapse(%d,%d): %d nodes removed in %d collapsed chains (edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
	   MinCoverage,verb,chaincnt,numchains, finalcnt, numtopnodes, wtime(),mtime());
    fflush(stdout);
  }

  if(DEBUG>=3){
    for(int i = 0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(VWedge->valid <= 0)
	  continue;
	if(!((VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid) && (VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap)))){
          printf("i=%d/%d:nodeid=%d(N%d), vw=%d/%d: valid=%d,edgeid=%u,Lmap=%d(N%d),Rmap=%d(N%d),node[i].Enodeid[vw]= %d\n",
                 i,numnodes,Vnode->nodeid,Vnode->mapid,vw,Vnode->edgeid[vw],Vnode->edgecnt,VWedge->valid,VWedge->Lmap,node[VWedge->Lmap].mapid,VWedge->Rmap,node[VWedge->Rmap].mapid,Vnode->Enodeid[vw]);
	  fflush(stdout);
          assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
	  assert(VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap));
        }
      }
    }    
  }
  if(DEBUG>=2)
    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      assert(Vnode->mark == 0);
      for(int vw= 0; vw < Vnode->edgecnt;vw++){
	Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	if(!(VWedge->mark == 0)){
	  printf("VWedge=(N%d,N%d,%0.3f,%d,%d):mark=%d still set at end of graph_collapse()\n",
		 nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,
		 VWedge->Lflipped,VWedge->Rflipped,VWedge->mark);
	  assert(VWedge->mark == 0);
	}
      }
    }

  PosEdgeInit = 0;/* Is this needed ? */

  return chaincnt;
}


class CCycle {/* Cycle information for a specific Vnode,Vflipped combination */
public:
  Cnode *Vnode,*Pnode;
  int Vflipped,Pflipped;

  Cedge *PNedge;/* back edge from Pnode to Nnode==Vnode to complete cycle */

  int PathLen;/* number of nodes in path from Vnode to Pnode (excluding Vnode) : <= 2*BFSdepth */
  Cedge **PathEdge;/* edges of path from Vnode to Pnode : PathEdge[k] is edge between PathNode[k] and its parent node */

  Cedge *DelEdge;/* edge scheduled for deletion (if previous cycle did not already break this cycle) */
  Cnode *DelNode;/* node of edge scheduled for deletion */
};

/* sort graph Cycles by Vnode->nodeid and preserve order of Vflipped etc */
static int CCycleNodeInc(const CCycle *p1, const CCycle *p2)
{
  return (p1->Vnode->nodeid > p2->Vnode->nodeid) ? 1 : (p1->Vnode->nodeid < p2->Vnode->nodeid) ? -1 :
    (p1 > p2) ? 1 : (p1 < p2) ? -1 : 0;
    //    (p1->Vflipped > p2->Vflipped) ? 1 : (p1->Vflipped < p2->Vflipped) ? -1 :
    //    (p1->Pnode->nodeid > p2->Pnode->nodeid) ? 1 : (p1->Pnode->nodeid < p2->Pnode->nodeid) ? -1 :
    //    p1->Pflipped - p2->Pflipped;
}

static int cyclebreak(Cnode *Vnode,
		      int Vflipped,
		      Cnode *Pnode,
		      int Pflipped,
		      Cedge *PNedge,
		      /* Nnode == Vnode */
		      int Nflipped,
		      CCycle *Cycles,
		      int &cntCycles,
		      int maxCycles,
		      int BFSdepth,
		      int verb,
		      int &cyclecnt,
		      int &cyclefix,
		      BFSdata *BFS) /**< if BFS != 0, use thread specific node memory BFS[node->nodeid].field instead of node->field. See DATA(BFS,node,field) */
{
  if(VERB>=2 && verb){
    printf("Found cycle from Vnode=N%d(flip=%d) to Pnode=N%d(flip=%d) to Vnode=N%d(flip=%d):V to P:cov=%0.3f,dist=%0.3f,wt=%0.3f:P to V:cov=%0.3f,dist=%0.3f,wt=%0.3f\n",
	   nodeorder[Vnode->mapid],Vflipped,nodeorder[Pnode->mapid],Pflipped,nodeorder[Vnode->mapid],Nflipped,DATA(BFS,Pnode,mincov),DATA(BFS,Pnode,Rdist),0.0,
	   PNedge->cov,PNedge->distance,0.0);
    fflush(stdout);
  }
  if(DEBUG>=2) assert(PNedge->cov > 0.0 && DATA(BFS,Pnode,mincov) > 0.0);

  #pragma omp atomic // NOTE : replace with local counter in parent
  cyclecnt++;

  int mycntCycles = 0;// NEW10

  if(PNedge->cov < DATA(BFS,Pnode,mincov) * CycleCoverageRatio){
    if(CYCLE_DEFER){/* save Cycle information and handle the actual edge deletion outside the multithreaded section (after check if this cycle still exists)  */
      #pragma omp atomic capture
      mycntCycles = cntCycles++;// NEW10

      if(mycntCycles >= maxCycles){
	#pragma omp atomic
	cntCycles--;
      } else {
	CCycle *pc = &Cycles[mycntCycles];
	pc->Vnode = Vnode;
	pc->Pnode = Pnode;
	pc->Vflipped = Vflipped;
	pc->Pflipped = Pflipped;
	pc->PNedge = PNedge;
	pc->DelEdge = PNedge;
	pc->DelNode = Vnode;
	int pathlen = 0;
	for(Cnode *Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
	  pc->PathEdge[pathlen] = DATA(BFS,Tnode,pedge);
	  pathlen++;
	}
	if(DEBUG) assert(pathlen <= BFSdepth*2);
	pc->PathLen = pathlen;
      }
    } else {
      #pragma omp critical
      {
        cyclefix++;
        edge_delete(PNedge,Vnode);
        if(VERB>=2 && verb){
  	  for(Cnode *Nnode = Pnode; Nnode != Vnode; Nnode = DATA(BFS,Nnode,parent))
	    printf("Nnode=N%d%c,mincov=%0.3f,pedge=(N%d,N%d,%0.3f,%d,%d)\n",
		   nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'':' ', DATA(BFS,Nnode,mincov), nodeorder[node[DATA(BFS,Nnode,pedge)->Lmap].mapid], 
		   nodeorder[node[DATA(BFS,Nnode,pedge)->Rmap].mapid], DATA(BFS,Nnode,pedge)->cov,
		   DATA(BFS,Nnode,pedge)->Lflipped,DATA(BFS,Nnode,pedge)->Rflipped);
	  printf("  (Deleted cycle link (N%d,N%d) with cov=%0.3f,valid=%d): cyclecnt=%d,cyclefix=%d)\n",
		 nodeorder[node[PNedge->Lmap].mapid],nodeorder[node[PNedge->Rmap].mapid], PNedge->cov, PNedge->valid, cyclecnt,cyclefix);
	}
      }
    }
    return 1;
  }

  /* locate smallest coverage link in cycle (excluding PNedge) and longest distance. Also determine average coverage weighted by distance */
  Cnode *Tnode = 0;
  double covsum = 0.0, distsum = 0.0;
  double mincov = 1.1*HUGE_COV;
  double maxdistance = PNedge->distance;
  for(Cnode *Nnode = Pnode; Nnode != Vnode; Nnode = DATA(BFS,Nnode,parent)){
    double cov = DATA(BFS,Nnode,pedge)->cov;
    double dist = DATA(BFS,Nnode,pedge)->distance;
    if(cov <= mincov){
      mincov = cov;
      Tnode = Nnode;
    }
    if(dist > maxdistance)
      maxdistance = dist;
    covsum += cov*dist;
    distsum += dist;
  }
  if(DEBUG>=2 && Tnode == 0){
    if(!(VERB>=2 && verb))
      printf("Found cycle from Vnode=N%d(flip=%d) to Pnode=N%d(flip=%d) to Vnode=N%d(flip=%d):V to P:cov=%0.1f,dist=%0.3f,wt=%0.3f:P to V:cov=%0.1f,dist=%0.3f,wt=%0.3f\n",
	     nodeorder[Vnode->mapid],Vflipped,nodeorder[Pnode->mapid],Pflipped,nodeorder[Vnode->mapid],Nflipped,DATA(BFS,Pnode,mincov),DATA(BFS,Pnode,Rdist),0.0,
	     PNedge->cov,PNedge->distance,0.0);

    printf("failed to find weakest link in cycle\n");
    fflush(stdout);
    assert(Tnode != 0);
  }
  if(DEBUG>=2 && !BULGE_PAR && !(fabs(mincov - DATA(BFS,Pnode,mincov)) <= 0.001)){
    #pragma omp critical
    {
      if(!(VERB>=2 && verb))
	printf("Found cycle from Vnode=N%d(flip=%d) to Pnode=N%d(flip=%d) to Vnode=N%d(flip=%d):V to P:mincov=%0.1f,Rist=%0.3f:P to V:cov=%0.1f,dist=%0.3f\n",
	       nodeorder[Vnode->mapid],Vflipped,nodeorder[Pnode->mapid],Pflipped,nodeorder[Vnode->mapid],Nflipped,DATA(BFS,Pnode,mincov),DATA(BFS,Pnode,Rdist),PNedge->cov,PNedge->distance);
      printf("Pnode->mincov=%0.3f, mincov=%0.3f, Tnode=N%d, Tnode->pedge=(N%d,N%d,cov=%0.3f,%d,%d,dist=%0.3f),verb=%d\n",
	     DATA(BFS,Pnode,mincov), mincov, nodeorder[Tnode->mapid], nodeorder[node[DATA(BFS,Tnode,pedge)->Lmap].mapid],nodeorder[node[DATA(BFS,Tnode,pedge)->Rmap].mapid],DATA(BFS,Tnode,pedge)->cov,
	     DATA(BFS,Tnode,pedge)->Lflipped,DATA(BFS,Tnode,pedge)->Rflipped,DATA(BFS,Tnode,pedge)->distance,verb);
      for(Cnode *Nnode = Pnode; Nnode != Vnode; Nnode = DATA(BFS,Nnode,parent))
	printf("Nnode=N%d%c,mincov=%0.3f,pedge=(id=%d,N%d,N%d,cov=%0.3f,%d,%d,dist=%0.3f)\n",
	       nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'':' ',DATA(BFS,Nnode,mincov), DATA(BFS,Nnode,pedge)->edgeid,
	       nodeorder[node[DATA(BFS,Nnode,pedge)->Lmap].mapid], nodeorder[node[DATA(BFS,Nnode,pedge)->Rmap].mapid], DATA(BFS,Nnode,pedge)->cov,
	       DATA(BFS,Nnode,pedge)->Lflipped,DATA(BFS,Nnode,pedge)->Rflipped,DATA(BFS,Nnode,pedge)->distance);
      fflush(stdout);
      assert(fabs(mincov - DATA(BFS,Pnode,mincov)) <= 0.001);
    }
  }
  if(0 && DEBUG>=2 && BULGE_PAR && !(mincov >= DATA(BFS,Pnode,mincov) - 0.001)){
    #pragma omp critical
    {
      if(!(VERB>=2 && verb))
	printf("Found cycle from Vnode=N%d(flip=%d) to Pnode=N%d(flip=%d) to Vnode=N%d(flip=%d):V to P:mincov=%0.1f,Rist=%0.3f:P to V:cov=%0.1f,dist=%0.3f\n",
	       nodeorder[Vnode->mapid],Vflipped,nodeorder[Pnode->mapid],Pflipped,nodeorder[Vnode->mapid],Nflipped,DATA(BFS,Pnode,mincov),DATA(BFS,Pnode,Rdist),PNedge->cov,PNedge->distance);
      printf("Pnode->mincov=%0.3f, mincov=%0.3f, Tnode=N%d, Tnode->pedge=(N%d,N%d,cov=%0.3f,%d,%d,dist=%0.3f),verb=%d\n",
	     DATA(BFS,Pnode,mincov), mincov, nodeorder[Tnode->mapid], nodeorder[node[DATA(BFS,Tnode,pedge)->Lmap].mapid],nodeorder[node[DATA(BFS,Tnode,pedge)->Rmap].mapid],DATA(BFS,Tnode,pedge)->cov,
	     DATA(BFS,Tnode,pedge)->Lflipped,DATA(BFS,Tnode,pedge)->Rflipped,DATA(BFS,Tnode,pedge)->distance,verb);
      for(Cnode *Nnode = Pnode; Nnode != Vnode; Nnode = DATA(BFS,Nnode,parent))
	printf("Nnode=N%d%c,mincov=%0.3f,pedge=(id=%d,N%d,N%d,cov=%0.3f,%d,%d,dist=%0.3f)\n",
	       nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'':' ',DATA(BFS,Nnode,mincov), DATA(BFS,Nnode,pedge)->edgeid,
	       nodeorder[node[DATA(BFS,Nnode,pedge)->Lmap].mapid], nodeorder[node[DATA(BFS,Nnode,pedge)->Rmap].mapid], DATA(BFS,Nnode,pedge)->cov,
	       DATA(BFS,Nnode,pedge)->Lflipped,DATA(BFS,Nnode,pedge)->Rflipped,DATA(BFS,Nnode,pedge)->distance);
      fflush(stdout);
      assert(mincov >= DATA(BFS,Pnode,mincov) - 0.001);/* cov  may have increased due to another thread since Pnode->mincov was computed */
    }
  }

  /* check if all edge lengths are < MAXERR*SM : if so delete the weakest edge */
  if(maxdistance < MAXERR*SM){
    if(CYCLE_DEFER){
      #pragma omp atomic capture
      mycntCycles = cntCycles++;// NEW10

      if(mycntCycles >= maxCycles){
	#pragma omp atomic
	cntCycles--;
      } else {
	CCycle *pc = &Cycles[mycntCycles];
	pc->Vnode = Vnode;
	pc->Pnode = Pnode;
	pc->Vflipped = Vflipped;
	pc->Pflipped = Pflipped;
	pc->PNedge = PNedge;
	if(PNedge->cov < mincov){
	  pc->DelEdge = PNedge;
	  pc->DelNode = Vnode;
	} else {
	  pc->DelEdge = DATA(BFS,Tnode,pedge);
	  pc->DelNode = Tnode;
	}
	int pathlen = 0;
	for(Cnode *Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
	  pc->PathEdge[pathlen] = DATA(BFS,Tnode,pedge);
	  pathlen++;
	}
	if(DEBUG) assert(pathlen <= BFSdepth*2);
	pc->PathLen = pathlen;
      }
    } else {
      #pragma omp critical
      {
	cyclefix++;
	if(PNedge->cov < mincov)
	  edge_delete(PNedge,Vnode);
	else
	  edge_delete(DATA(BFS,Tnode,pedge),Tnode);
      }
    }
    return 1;
  }
  
  if(MAXCYCLE && DATA(BFS,Pnode,Rdist) + PNedge->distance > MaxBulge){/* recheck when MaxBulge is increased */
    if(VERB>=2 && verb){
      printf("   cycle length=%0.3f is too large (MaxBulge=%0.3f)\n",DATA(BFS,Pnode,Rdist)+PNedge->distance,MaxBulge);
      for(Cnode *Nnode = Pnode; Nnode != Vnode; Nnode = DATA(BFS,Nnode,parent))
	printf("Nnode=N%d%c,mincov=%0.3f,pedge=(N%d,N%d,%0.3f,%d,%d)\n",
	       nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'':' ', DATA(BFS,Nnode,mincov), nodeorder[node[DATA(BFS,Nnode,pedge)->Lmap].mapid], 
	       nodeorder[node[DATA(BFS,Nnode,pedge)->Rmap].mapid], DATA(BFS,Nnode,pedge)->cov, DATA(BFS,Nnode,pedge)->Lflipped,DATA(BFS,Nnode,pedge)->Rflipped);
      fflush(stdout);
    }
    return 0;
  }

  /* Possible improvement : look for alternate path with highest coverage from Tnode->parent to Tnode AND Pnode to Tnode (excluding edge Tnode->pedge) */

  if(CYCLEBREAK){
    if(mincov < PNedge->cov){
      covsum -= mincov * DATA(BFS,Tnode,pedge)->distance;
      distsum -= DATA(BFS,Tnode,pedge)->distance;
      covsum += PNedge->cov * PNedge->distance;
      distsum += PNedge->distance;
    }
    covsum /= distsum;

    if(mincov < PNedge->cov){
      if(mincov < covsum * CycleCoverageRatio && mincov <= CycleBreakMaxCov) {
	if(CYCLE_DEFER){/* save Cycle information and handle the actual edge deletion outside the multithreaded section (after check if this cycle still exists)  */
          #pragma omp atomic capture
	  mycntCycles = cntCycles++;// NEW10

	  if(mycntCycles >= maxCycles){
	    #pragma omp atomic
	    cntCycles--;
	  } else {
	    CCycle *pc = &Cycles[mycntCycles];
	    pc->Vnode = Vnode;
	    pc->Pnode = Pnode;
	    pc->Vflipped = Vflipped;
	    pc->Pflipped = Pflipped;
	    pc->PNedge = PNedge;
	    pc->DelEdge = DATA(BFS,Tnode,pedge);
	    pc->DelNode = Tnode;
	    int pathlen = 0;
	    for(Cnode *Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
	      pc->PathEdge[pathlen] = DATA(BFS,Tnode,pedge);
	      pathlen++;
	    }
	    if(DEBUG) assert(pathlen <= BFSdepth*2);
	    pc->PathLen = pathlen;
	  }
	} else {
          #pragma omp critical
	  {
	    cyclefix++;
	    edge_delete(DATA(BFS,Tnode,pedge),Tnode);
	    if(VERB>=2 && verb){
	      for(Cnode *Nnode = Pnode; Nnode != Vnode; Nnode = DATA(BFS,Nnode,parent))
		printf("Nnode=N%d%c,mincov=%0.3f,pedge=(N%d,N%d,%0.3f,%d,%d)\n",
		       nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'':' ', DATA(BFS,Nnode,mincov), nodeorder[node[DATA(BFS,Nnode,pedge)->Lmap].mapid],
		       nodeorder[node[DATA(BFS,Nnode,pedge)->Rmap].mapid], DATA(BFS,Nnode,pedge)->cov, DATA(BFS,Nnode,pedge)->Lflipped,DATA(BFS,Nnode,pedge)->Rflipped);
	      printf("  (Deleted cycle link (N%d,N%d) with cov min=%0.1f,av=%0.1f ): cyclecnt=%d,cyclefix=%d)\n",
		     nodeorder[node[DATA(BFS,Tnode,pedge)->Lmap].mapid],nodeorder[node[DATA(BFS,Tnode,pedge)->Rmap].mapid], mincov, covsum, cyclecnt,cyclefix);
	    }
	  }
	}
	return 1;
      } else {
	if(VERB>=2 && verb){
	  if(Nflipped != Vflipped)
	    printf("  link (N%d to N%d with cov=%0.1f) is not sufficiently weaker than average coverage of other links = %0.1f and MaxCov=%0.1f (cycle with reversal!)\n",
		   nodeorder[node[DATA(BFS,Tnode,pedge)->Lmap].mapid], nodeorder[node[DATA(BFS,Tnode,pedge)->Rmap].mapid], mincov, covsum, CycleBreakMaxCov);
	  else
	    printf("  link (N%d to N%d with cov=%0.1f) is not sufficiently weaker than average coverage of other links = %0.1f and MaxCov=%0.1f\n",
		   nodeorder[node[DATA(BFS,Tnode,pedge)->Lmap].mapid], nodeorder[node[DATA(BFS,Tnode,pedge)->Rmap].mapid], mincov, covsum, CycleBreakMaxCov);
	  for(Cnode *Nnode = Pnode; Nnode != Vnode; Nnode = DATA(BFS,Nnode,parent))
	    printf("Nnode=N%d%c,mincov=%0.3f,pedge=(N%d,N%d,%0.3f,%d,%d,cov=%0.1f)\n",
		   nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'':' ',DATA(BFS,Nnode,mincov), nodeorder[node[DATA(BFS,Nnode,pedge)->Lmap].mapid], 
		   nodeorder[node[DATA(BFS,Nnode,pedge)->Rmap].mapid], DATA(BFS,Nnode,pedge)->distance, DATA(BFS,Nnode,pedge)->Lflipped,DATA(BFS,Nnode,pedge)->Rflipped, DATA(BFS,Nnode,pedge)->cov);
	  fflush(stdout);
	}
      
	return 0;
      }
    } else {
      if(PNedge->cov < covsum * CycleCoverageRatio && PNedge->cov <= CycleBreakMaxCov){
	if(CYCLE_DEFER){/* save Cycle information and handle the actual edge deletion outside the multithreaded section (after check if this cycle still exists)  */
          #pragma omp atomic capture
	  mycntCycles = cntCycles++;// NEW10

	  if(mycntCycles >= maxCycles){
            #pragma omp atomic
	    cntCycles--;
	  } else {
	    CCycle *pc = &Cycles[mycntCycles];
	    pc->Vnode = Vnode;
	    pc->Pnode = Pnode;
	    pc->Vflipped = Vflipped;
	    pc->Pflipped = Pflipped;
	    pc->PNedge = PNedge;
	    pc->DelEdge = PNedge;
	    pc->DelNode = Vnode;
	    int pathlen = 0;
	    for(Cnode *Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
	      pc->PathEdge[pathlen] = DATA(BFS,Tnode,pedge);
	      pathlen++;
	    }
	    if(DEBUG) assert(pathlen <= BFSdepth*2);
	    pc->PathLen = pathlen;
	  }
	} else {
	  #pragma omp critical
	  {
	    cyclefix++;
	    edge_delete(PNedge, Vnode);
	    if(VERB>=2 && verb)
	      printf("  (Deleted cycle link (N%d,N%d) with cov min=%0.1f,av=%0.1f ): cyclecnt=%d,cyclefix=%d)\n",
		     nodeorder[node[PNedge->Lmap].mapid],nodeorder[node[PNedge->Rmap].mapid], PNedge->cov, covsum, cyclecnt,cyclefix);
	  }
	}
	return 1;
      } else {
	if(VERB>=2 && verb){
	  if(Nflipped != Vflipped)
	    printf("  link (N%d to N%d with cov=%0.1f) is not sufficiently weaker than average coverage of other links = %0.1f (cycle with reversal!)\n",
		   nodeorder[Pnode->mapid], nodeorder[Vnode->mapid], PNedge->cov, covsum);
	  else
	    printf("  link (N%d to N%d with cov=%0.1f) is not sufficiently weaker than average coverage of other links = %0.1f\n",
		   nodeorder[Pnode->mapid], nodeorder[Vnode->mapid], PNedge->cov, covsum);
	  for(Cnode *Nnode = Pnode; Nnode != Vnode; Nnode = DATA(BFS,Nnode,parent))
	    printf("Nnode=N%d%c,mincov=%0.3f,pedge=(N%d,N%d,%0.3f,%d,%d,cov=%0.1f)\n",
		   nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'':' ',DATA(BFS,Nnode,mincov), nodeorder[node[DATA(BFS,Nnode,pedge)->Lmap].mapid],
		   nodeorder[node[DATA(BFS,Nnode,pedge)->Rmap].mapid], DATA(BFS,Nnode,pedge)->distance,
		   DATA(BFS,Nnode,pedge)->Lflipped,DATA(BFS,Nnode,pedge)->Rflipped,DATA(BFS,Nnode,pedge)->cov);
	  fflush(stdout);
	}
      
	return 0;
      }
    }
  } 

  if(DEBUG) assert(CYCLEBREAK==0);

  printf("CYCLEBREAK=0 no longer supported\n");
  fflush(stdout);
  exit(1);

  return 0;
}

class CBulge {/* Bulge information for a specific Vnode,Vflipped combination : PathLen1 == 0 if there is no bulge information */
public:
  Cnode *Vnode, *Xnode, *Wnode;
  int Vflipped, Xflipped, Wflipped;

  int PathLen1;/* number of nodes in path from Vnode to Xnode to Nnode (excluding Vnode) : <= BFSdepth */
  int *PathNode1;/* node ids of path from Vnode to Nnode (excluding Vnode) : PathNode1[0] == Nnode->nodeid, PathNode1[k+1] is the parent of PathNode1[k] */
  int *PathFlipped1;/* orientation of nodes in path from Vnode to Nnode (excluding Vnode) */
  Cedge **PathEdge1;/* edges of path from Vnode to Nnode (excluding Vnode) : PathEdge1[k] is edge between PathNode1[k] and its parent node */

  int PathLen2;/* number of nodes in path from Vnode to Pnode (excluding Vnode) : <= BFSdepth */
  int *PathNode2;/* node ids of path from Vnode to Pnode (excluding Vnode) : PathNode2[0] == Pnode->nodeid, PathNode1[k+1] is the parent of PathNode1[k] */
  int *PathFlipped2;/* orientation of nodes in path from Vnode to Pnode (excluding Vnode) */
  Cedge **PathEdge2;/* edges of path from Vnode to Pnode (excluding Vnode) : PathEdge2[k] is edge between PathNode2[k] and its parent node */

  Cedge *PNedge;
  
  double cov1,cov2;/* planned changes to path1 and path2 respectively */
};

/* sort graph bulges by Vnode->nodeid and preserve order of Vflipped, Xnode->nodeid, Xflipped, Wnode->nodeid, Wflipped */
static int CBulgeNodeInc(const CBulge *p1, const CBulge *p2)
{
  return (p1->Vnode->nodeid > p2->Vnode->nodeid) ? 1 : (p1->Vnode->nodeid < p2->Vnode->nodeid) ? -1 :
    (p1 > p2) ? 1 : (p1 < p2) ? -1 : 0;
    //    (p1->Vflipped > p2->Vflipped) ? 1 : (p1->Vflipped < p2->Vflipped) ? -1 :
    //    (p1->Xnode->nodeid > p2->Xnode->nodeid) ? 1 : (p1->Xnode->nodeid < p2->Xnode->nodeid) ? -1 :
    //    (p1->Xflipped > p2->Xflipped) ? 1 : (p1->Xflipped < p2->Xflipped) ? -1 :
    //    (p1->Wnode->nodeid > p2->Wnode->nodeid) ? 1 : ( p1->Wnode->nodeid < p2->Wnode->nodeid) ? -1 :
    //    p1->Wflipped - p2->Wflipped;
}

/** Looks for graph bulges:

   1. Simple bulges (two (Super)Edges connecting the same two nodes): 
         If distances match, merge lower coverage path. 
         If distances don't match and lower coverage path is well below coverage of other path, delete it.

   2. General Bulge (two paths which connecting the same two nodes, found by BFS with limited number of hops):
         If distances match, delete lowest coverage edge and add this coverage to other path and subtract it from this path.
	 If distances don't match and lower coverage path is well below coverage of other path (by MinCoverageRatio) 
    	    delete lowest coverage edge (provided the coverage is not more than BulgeBreakMaxCov)
         This can create duplicate contigs when the edge deleted is the only (or highest coverage) forward or backward edge from a node : 
	 If BULGE_AVOID try to avoid deleting these kinds of edges by finding another possible edge to delete (even it is not the lowest coverage edge).
	 

   3. Cycles : break weakest link based, if possible (see CYCLEBREAK)

   NOTE : choice of weakest link to break when path lengths do not match should be based on more than just coverage. For example if an edge
       with cov=44 is going to be subsequently thresholded by EdgeCoverageRatio, while another edge with cov=20 is not going to be threshold by
       EdgeCoverageRatio, we should delete the edge with cov=44 rather than the edge with cov=20 (which may be near a fragile site or low coverage ratio).
       Other measures to consider:
       1. Ratio to highest coverage edge in same direction (from either node), as used by EdgeCoverageRatio
       2. Ratio of max(1,Ccnt) to min(Lcntmax,Rcntmax) of EITHER/BOTH (see EDGETRIM=1,2) nodes, as used by graph_edgetrim().
*/
int graph_bulge(int BFSdepth, int verb)/* Maximum depth of BFS used to find complex bulges */
{
  if(VERB>=2 || (DEBUG>=2 && GraphWeight > 0.0)){/* locate highest coverage edge */
    double maxCov = 0.0;
    Cedge *maxEdge = 0;
    for(int t=0; t < numtopnodes; t++){
      Cnode *Tnode = topnodes[t];
      for(int u = 0; u < Tnode->edgecnt; u++){
	Cedge *TUedge = &edge[Tnode->edgeid[u]];
	if(DEBUG) assert(TUedge->cov >= 0.0);
	if(TUedge->cov > maxCov){
	  maxCov = TUedge->cov;
	  maxEdge = TUedge;
	}
      }	
    }
    printf("graph_bulge: Largest coverage = %0.2f : (N%d,N%d,%0.3f,%d,%d)\n",
	   maxCov, nodeorder[node[maxEdge->Lmap].mapid],nodeorder[node[maxEdge->Rmap].mapid],maxEdge->distance,
	   maxEdge->Lflipped,maxEdge->Rflipped);
    fflush(stdout);
    if(DEBUG) assert(maxCov < HUGE_COV);
  }

  int bulgecnt = 0;  /* count how many simple bulges exist */
  int bulgefix = 0;  /* count how many bulges could be fixed (removed) */
  int cyclecnt = 0; /* count number of cycles */
  int cyclefix = 0; /* count number of cycles that could be fixed */

  /* first check for simple graph bulges */
  for(int i=0; i < numtopnodes; i++){
    Cnode *Vnode = topnodes[i];
    if(DEBUG>=2) assert(Vnode->mark != 2);
    if(Vnode->edgecnt < 2)
      continue;
    int V = Vnode->nodeid;
    for(int vw = 0; vw < Vnode->edgecnt; vw++){
      Cedge *VWedge = &edge[Vnode->edgeid[vw]];
      if(VWedge->valid < 0)
	continue;
      int W = (VWedge->Lmap == V) ? VWedge->Rmap : VWedge->Lmap;
      if(DEBUG>=2) assert(V != W);
      if(V > W)
	continue;/* Will be checked with V,W reversed */
      for(int vx = vw+1; vx < Vnode->edgecnt; vx++){
	Cedge *VXedge = &edge[Vnode->edgeid[vx]];
	if(VXedge->valid < 0)
	  continue;
	int X = (VXedge->Lmap == V) ? VXedge->Rmap : VXedge->Lmap;
	if(X==W){
	  bulgecnt++;
	  /* check if distance of VW and VX are consistant */
	  Cnode *Xnode = &node[X];
	  double var = 2.0*SM*SM + SD[0]*SD[0]*(Vnode->len + Xnode->len);
	  double err = VXedge->distance - VWedge->distance;
	  if(VERB>=2 && verb)
	    printf("bulge%d:V=N%d,W=X=N%d:VWedge=(N%d,N%d,%0.3f,%d,%d,C=%0.1f,S=%d),VXedge=(N%d,N%d,%0.3f,%d,%d,C=%0.1f,S=%d),sd=%0.3f,err=%0.3f\n",
		   bulgecnt,nodeorder[Vnode->mapid],nodeorder[node[W].mapid],
		   nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,VWedge->Lflipped,VWedge->Rflipped,VWedge->cov,VWedge->alignid ? 0 : 1,
		   nodeorder[node[VXedge->Lmap].mapid],nodeorder[node[VXedge->Rmap].mapid],VXedge->distance,VXedge->Lflipped,VXedge->Rflipped,VXedge->cov,VXedge->alignid ? 0 : 1,sqrt(var),err);
	  if(err*err >= var*(NMAXERR*NMAXERR)) {
	    if(VXedge->cov > VWedge->cov * MinCoverageRatio && 
	       VWedge->cov > VXedge->cov * MinCoverageRatio && 
	       (VWedge->Lflipped||VWedge->Rflipped) == (VXedge->Lflipped||VXedge->Rflipped))
	      continue;
	    bulgefix++;
	    if(VXedge->cov <= VWedge->cov){
	      if(DEBUG>=2) assert(VXedge->valid > 0);
	      edge_delete(VXedge,Vnode);/* mark VXedge for deletion (adding any super-Edge nodes to side chain of Vnode) */
	      continue;
	    }
	    if(VWedge->cov < VXedge->cov){
	      if(DEBUG>=2) assert(VWedge->valid > 0);
	      edge_delete(VWedge,Vnode);
	      continue;
	    }
	    printf("bulge%d:V=N%d,W=X=N%d:VWedge=(N%d,N%d,%0.3f,%d,%d,C=%0.1f,S=%d),VXedge=(N%d,N%d,%0.3f,%d,%d,C=%0.1f,S=%d):SD=%0.3f,err=%0.3f(not implemented)\n",
		   bulgecnt,nodeorder[Vnode->mapid],nodeorder[node[W].mapid],
		   nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,VWedge->Lflipped,VWedge->Rflipped,VWedge->cov,VWedge->alignid ? 0 : 1,
		   nodeorder[node[VXedge->Lmap].mapid],nodeorder[node[VXedge->Rmap].mapid],VXedge->distance,VXedge->Lflipped,VXedge->Rflipped,VXedge->cov,VXedge->alignid ? 0 : 1,
		   sqrt(var),fabs(err));
	    assert(0);
	  }

	  /* If distances are larger than 100% of map size, we also need to check if the maps are consistent (Not implemented: just print a warning) */
	  if(VXedge->distance + VWedge->distance > 2.0 * min(Vnode->len,Xnode->len)){
	    printf("graph_bulge:Warning:V=N%d,W=X=N%d:length = %0.3f,%0.3f exceed 100%% of map lengths = %0.3f,%0.3f: Not checking entire path for matching maps\n",
		   nodeorder[Vnode->mapid],nodeorder[node[W].mapid],VXedge->distance,VWedge->distance,Vnode->len,Xnode->len);
	    fflush(stdout);
	  }

	  bulgefix++;
	  if(VWedge->cov >= VXedge->cov){/* remove VXedge and add its coverage to other edge */
	    VWedge->cov += VXedge->cov;
	    if(DEBUG>=2) assert(VXedge->valid > 0);
	    edge_delete(VXedge,Vnode);/* mark VXedge for deletion (adding any super-Edge nodes to side chain of Vnode) */
	  } else {/* remove VWedge and add its coverage to other edge */
	    VXedge->cov += VWedge->cov;
	    if(DEBUG>=2) assert(VWedge->valid > 0);
	    edge_delete(VWedge,Vnode);

	    break; /* Since VWedge is no longer valid, break out of inner loop  and advance to next VWedge */
	  }
	}
      }
    }
  }
  
  /* remove edges with valid <= -1 from edge lists of all nodes */
  unsigned int finalcnt;
  (void)GCedges(finalcnt,edge);

  if(VERB && (bulgecnt > 0 || verb)){
    printf("graph_bulge(%d,%0.1f): %d simple graph bulges found(%d bulge edges deleted)(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
	   BFSdepth, MaxBulge, bulgecnt, bulgefix, finalcnt, numtopnodes,wtime(),mtime());
    fflush(stdout);
  }
  
  int Vcnt = 0;
  int BFScnt1= 0, BFScnt2 = 0;

  bulgecnt = 0;/* count number of complex bulges */
  bulgefix = 0;/* count number of edges deleted */

  if(BFSdepth >= 1){/* perform BFS from each node up to depth BFSdepth + 1 */
    int origPosEdgeInit= PosEdgeInit;
    if(!PosEdgeInit)    /* update PosEdge of remaining topnodes */
      PosEdge();

    if(DEBUG>=2){/* check that node distances are correctly ordered */
      for(int i=0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	for(int Vflipped = 0;Vflipped <= 1; Vflipped++){/* try both orientations of V */      
	  for(int vw= (Vflipped ? 0 : Vnode->PosEdge); vw < (Vflipped ? Vnode->PosEdge : Vnode->edgecnt); vw++){
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];	
	    int NDsign = (VWedge->Lmap == Vnode->nodeid) ? (VWedge->Lflipped ? -1 : 1) : (VWedge->Rflipped ? 1 : -1);
	    if(DEBUG && !(NDsign == (Vflipped ? -1 : 1))){
	      printf("Error in order of edges for Vnode=N%d,Vflipped=%d,vw=%d,PosEdge=%d:NDsign=%d,distance=%0.3f\n",
		     nodeorder[Vnode->mapid],Vflipped,vw,Vnode->PosEdge,NDsign,VWedge->distance);
	      printf("PosEdgeInit=%d->%d\n",origPosEdgeInit,PosEdgeInit);
	      fflush(stdout);
	      assert(NDsign == (Vflipped ? -1 : 1));
	    }
	  }
	}
      }
    }

    if(DEBUG>=2){/* all nodes should have parent == 0, mark == 0, SPmark == 0 */
      for(int t=0; t < numtopnodes; t++){
	Cnode *Tnode = topnodes[t];
	assert(Tnode->parent == 0);
	assert(Tnode->mark == 0);
	assert(Tnode->SPmark == 0);
	/* all edges should have valid == 1 */
	for(int u = 0; u < Tnode->edgecnt; u++){
	  Cedge *TUedge = &edge[Tnode->edgeid[u]];
	  assert(TUedge->valid == 1);
	  assert(TUedge->cov > -0.01);
	}
      }
    }

    /* Perform BFS based on distance up to distance MaxBulge (OR BFSdepth*2 steps) to locate cycles */

    /* allocate memory for cycle information for each Vnode & Vflipped */
    CCycle *Cycles = 0;
    Cedge **CEdgeMem = 0;
    int maxCycles = 0;
    int cntCycles = 0;
    if(CYCLE_DEFER){
      maxCycles = numtopnodes*4;/* Note : If more cycles are found the excess will be ignored till the next call */
      Cycles = new CCycle[maxCycles];
      CEdgeMem = new Cedge*[maxCycles*(BFSdepth*2)];
      int EdgeCnt = 0;
      for(int i = 0; i < maxCycles; i++){
	CCycle *p = &Cycles[i];
	p->PathLen = 0;
	p->PathEdge = &CEdgeMem[EdgeCnt]; EdgeCnt += 2*BFSdepth;
      }
      if(DEBUG) assert(EdgeCnt == maxCycles*BFSdepth*2);
    }

    int nthreads = min(numthreads, max(1,numtopnodes/8));
    int block = 1;
    while(block < 16 && numtopnodes > numthreads * block * 8)
      block *= 2;

    #pragma omp parallel num_threads(nthreads) if(CYCLE_PAR && nthreads > 1)
    { 
      /* allocate per-thread memory */
      BFSdata *BFS = new BFSdata[numnodes];
      for(int i=0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	DATA(BFS,Vnode,parent) = 0;
	DATA(BFS,Vnode,mark) = 0;
	DATA(BFS,Vnode,SPmark) = 0;
      }

      //      int mycntCycles = 0;

      #pragma omp for schedule(dynamic,block)
      for(int i=0; i < numtopnodes; i++){
        if(CYCLE_DEFER && !TSAN){
      //          #pragma omp critical
      //          mycntCycles = cntCycles;

	  if(cntCycles >= maxCycles)
	    continue;//NOTE : cannot break out of multithreaded loop, so just skip all subsequent iterations
        }

	Cnode *Vnode = topnodes[i];

	if(DEBUG>=2) assert(DATA(BFS,Vnode,parent) == 0);
	DATA(BFS,Vnode,Rdist) = 0.0;
	DATA(BFS,Vnode,mincov) = 1000000000;
	if(VERB>=2 && !(i%1000)){
	  #pragma omp critical
	  {
	    printf("i=%d/%d: BFSdepth=%d,MaxBulge=%0.1f: %d cycles found (%d edges deleted):time=%0.6f(total=%0.6f)\n",i,numtopnodes,BFSdepth,MaxBulge,cyclecnt,cyclefix,wtime(),mtime());
	    fflush(stdout);
	  }
	}

	/* look for cycles : due to reversal edges need any two edges */
	if(Vnode->edgecnt >= 2){
	  int found = 0;
	  DATA(BFS,Vnode,parent) = Vnode;
	  for(int Vflipped = 0;Vflipped <= 1; Vflipped++){/* try both orientations of V */
	    DATA(BFS,Vnode,pflip) = Vflipped;
	    /* If Vflipped==1: look for neighbors of Vnode with -ve nodedistance 
	       If Vflipped==0: look for neighbors of Vnode with +ve node distance */
	    for(int vw= (Vflipped ? 0 : Vnode->PosEdge); vw < (Vflipped ? Vnode->PosEdge : Vnode->edgecnt); vw++){
	      Cedge *VWedge = &edge[Vnode->edgeid[vw]];	
	      if(VWedge->valid < 0)
		continue;
	      Cnode *Wnode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
	      int Wflipped = (Vflipped ^ (VWedge->Lflipped ^ VWedge->Rflipped));
	  
	      if(VERB>=2 && verb>=2){
		printf("Vnode=N%d(flip=%d),Wnode=N%d(flip=%d)\n",
		       nodeorder[Vnode->mapid],Vflipped,nodeorder[Wnode->mapid],Wflipped);
		fflush(stdout);
	      }

	      if(DEBUG>=2)assert(DATA(BFS,Wnode,parent) == 0);
	      DATA(BFS,Wnode,parent) = Vnode;
	      DATA(BFS,Wnode,pedge) = VWedge;
	      DATA(BFS,Wnode,Rdist) = VWedge->distance;
	      DATA(BFS,Wnode,mincov) = VWedge->cov;
	      DATA(BFS,Wnode,mark) = 1;
	      DATA(BFS,Wnode,pflip) = Wflipped;

	      priority_queue< Cpath, vector<Cpath> , less<Cpath> > BFSqueue;/* NATE : less (<) is implemented as greater (>) in Cpath! */
	      stack<Cnode*> FinishedNodes;/* stack of nodes removed from BFSqueue : used to reset parent field later */
	      BFSqueue.push(Cpath(Wnode,Wflipped,CYCLESTEPS ? 1.0 : DATA(BFS,Wnode,Rdist)));
	      if(VERB>=2 && verb>=2) printf("queue.push:priority=%0.1f,Wnode=N%d(flip=%d),parent=N%d(flip=%d),mincov=%0.1f,Rdist=%0.3f,Rwt=%0.3f,mark=%d\n",
					 1.0,nodeorder[Wnode->mapid],DATA(BFS,Wnode,pflip),nodeorder[Vnode->mapid],Vflipped,
					 DATA(BFS,Wnode,mincov),DATA(BFS,Wnode,Rdist),0.0,DATA(BFS,Wnode,mark));

	      while(BFSqueue.size() > 0){
		Cpath top = BFSqueue.top(); BFSqueue.pop();
		FinishedNodes.push(top.node);
		if(VERB>=2)
		  BFScnt1++;
		Cnode *Pnode = top.node;
		int Pflipped = top.flipped;
		if(DEBUG>=2) assert(DATA(BFS,Pnode,parent));
		if(VERB>=2 && verb>=2) printf("queue.pop():priority=%0.1f,Pnode=N%d(flip=%d),parent=N%d(flip=%d),mincov=%0.1f,Rdist=%0.3f,Rwt=%0.3f,mark=%d\n",
					   top.priority,nodeorder[Pnode->mapid],Pflipped,nodeorder[DATA(BFS,Pnode,parent)->mapid],DATA(BFS,DATA(BFS,Pnode,parent),pflip),
					   DATA(BFS,Pnode,mincov),DATA(BFS,Pnode,Rdist),0.0,DATA(BFS,Pnode,mark));
		if(top.priority >= (CYCLESTEPS ? BFSdepth*2.0 : MaxBulge))
		  break;
		if(DATA(BFS,Pnode,Rdist) >= MaxBulge)
		  continue;

		/* If Pflipped==1: look for neighbors of Pnode with -ve nodedistance 
		   If Pflipped==0: look for neighbors of Pnode with +ve node distance */
		int pn;
		for(pn = (Pflipped ? 0 : Pnode->PosEdge); pn < (Pflipped ? Pnode->PosEdge : Pnode->edgecnt); pn++){
		  Cedge *PNedge = &edge[Pnode->edgeid[pn]];
		  Cnode *Nnode = &node[(PNedge->Lmap == Pnode->nodeid) ? PNedge->Rmap : PNedge->Lmap];
		  int Nflipped = (Pflipped ^ (PNedge->Lflipped | PNedge->Rflipped));
		  if(VERB>=2 && verb>=2)
		    printf("pn=%d/%d:Nnode=N%d(flip=%d):PNedge->valid=%d,cov=%0.3f,Nnode->parent=N%d,Nnode->mark=%d\n",
			   pn,(Pflipped ? Pnode->PosEdge : Pnode->edgecnt),
			   nodeorder[Nnode->mapid],Nflipped,PNedge->valid,PNedge->cov,
			   DATA(BFS,Nnode,parent) ? nodeorder[DATA(BFS,Nnode,parent)->mapid] : -1, DATA(BFS,Nnode,mark));
		  if(PNedge->valid < 0)
		    continue;

		  double distance2 = DATA(BFS,Pnode,Rdist) + PNedge->distance;
		  double mincov2 = min(DATA(BFS,Pnode,mincov),PNedge->cov);
	      
		  if(Nnode == Vnode){/* a cycle : try to break the weakest link(s) */
		    if(cyclebreak(Vnode,Vflipped,Pnode,Pflipped,PNedge,Nflipped,Cycles,cntCycles,maxCycles,BFSdepth,verb,cyclecnt,cyclefix,BFS)){
		      found = 1;
		      break;/* terminate BFS */
		    }
		    continue;
		  }
		  if(DATA(BFS,Nnode,mark) == DATA(BFS,Pnode,mark))
		    continue;/* self-bulge (will be found later when Vnode is the common ancestor of Nnode and Pnode) */
		  if(DEBUG>=2) assert(0 <= Pnode - node && Pnode - node < numnodes && Pnode-node == Pnode->nodeid);
		  DATA(BFS,Nnode,parent) = Pnode; 
		  DATA(BFS,Nnode,pedge) = PNedge;
		  DATA(BFS,Nnode,Rdist) = distance2;//Pnode->Rdist + PNedge->distance;
		  if(DEBUG>=2) assert(mincov2 == min(DATA(BFS,Pnode,mincov),PNedge->cov));
		  DATA(BFS,Nnode,mincov) = mincov2;//min(Pnode->mincov,PNedge->cov);
		  DATA(BFS,Nnode,mark) = DATA(BFS,Pnode,mark);
		  DATA(BFS,Nnode,pflip) = Nflipped;
		  if(VERB>=2 && verb>=2){
		    printf("queue.push:priority=%0.1f,Node=N%d(flip=%d),parent=N%d,mincov=%0.1f,Rdist=%0.3f,wt=%0.3f,mark=%d\n",
			   top.priority+1.0,nodeorder[Nnode->mapid],Nflipped,nodeorder[Pnode->mapid],
			   DATA(BFS,Nnode,mincov),DATA(BFS,Nnode,Rdist),0.0,DATA(BFS,Nnode,mark));
		    fflush(stdout);
		  }
		  BFSqueue.push(Cpath(Nnode,Nflipped,CYCLESTEPS ? top.priority + 1.0 : DATA(BFS,Nnode,Rdist)));	      
		}
		if(pn < (Pflipped ? Pnode->PosEdge : Pnode->edgecnt))/* early termination */
		  break;/* terminate BFS */
	      }
	      /* clean up : reset parent fields in FinishedNodes or BFSqueue */
	      while(!FinishedNodes.empty()){
		Cnode *Pnode = FinishedNodes.top(); FinishedNodes.pop();
		if(DEBUG>=2 && !DATA(BFS,Pnode,parent)){
		  printf("Finished node N%d :parent=0\n",nodeorder[Pnode->mapid]);
		  fflush(stdout);
		  assert(DATA(BFS,Pnode,parent));
		}
		if(VERB>=2 && verb>=2){
		  printf("reseting Finished node N%d: parent=N%d,mark=%d\n",
			 nodeorder[Pnode->mapid],nodeorder[DATA(BFS,Pnode,parent)->mapid],DATA(BFS,Pnode,mark));
		  fflush(stdout);
		}
		DATA(BFS,Pnode,parent) = 0;
		DATA(BFS,Pnode,mark) = 0;
	      }
	      while(BFSqueue.size() > 0){
		Cpath top = BFSqueue.top(); BFSqueue.pop();
		Cnode *Pnode = top.node;
		if(DEBUG>=2)assert(DATA(BFS,Pnode,parent));
		if(VERB>=2 && verb>=2)
		  printf("queue.pop():priority=%0.1f,Node=N%d,parent=N%d,mincov=%0.1f,Rdist=%0.3f,mark=%d(cleanup)\n",
			 top.priority,nodeorder[Pnode->mapid],nodeorder[DATA(BFS,Pnode,parent)->mapid],
			 DATA(BFS,Pnode,mincov),DATA(BFS,Pnode,Rdist),DATA(BFS,Pnode,mark));
		DATA(BFS,Pnode,parent) = 0;
		DATA(BFS,Pnode,mark) = 0;
	      }
	      if(DEBUG>=3){/* all nodes except Vnode should have parent field == 0 */
		for(int t=0; t < numtopnodes; t++){
		  Cnode *Tnode = topnodes[t];
		  if(Tnode != Vnode)
		    assert(DATA(BFS,Tnode,parent) == 0);
		  assert(DATA(BFS,Tnode,mark) == 0);
		  assert(DATA(BFS,Tnode,SPmark) == 0);
		}
	      }
	      if(found)
		break;
	    }
	    if(found)
	      break;
	  }
	  DATA(BFS,Vnode,parent) = 0;
	}
      } // omp for i = 0 .. numtopnodes-1
      /* free up per thread memory */
      if(BFS) delete [] BFS;
    } // omp parallel

    if(CYCLE_DEFER){/* perform graph updates to break cyles in a deterministic order */
      if(cntCycles >= maxCycles - numthreads){
	printf("WARNING:increase maxCycles=%d (numtopnodes=%d,cntCycles=%d)\n",maxCycles,numtopnodes,cntCycles);
	fflush(stdout);
      }

      qsort(Cycles,cntCycles,sizeof(CCycle),(intcmp*)CCycleNodeInc);/* make order deterministic */

      for(int i = 0; i < cntCycles; i++){
	CCycle *pc = &Cycles[i];

	/* check that all cycle edges still have valid >= 1 */
	if(pc->DelEdge->valid <= 0)
	  continue;
	if(DEBUG>=2) assert(pc->PathLen > 0);
	int valid = 1;
	for(int pathlen = 0; pathlen < pc->PathLen; pathlen++){
	  Cedge *Tpedge = pc->PathEdge[pathlen];
	  if(!(Tpedge->valid > 0)){
	    valid = 0;
	    break;
	  }
	}
	if(!valid)
	  continue;
	cyclefix++;
	edge_delete(pc->DelEdge,pc->DelNode);
      }

      /* reset Cycles[] memory */
      cntCycles = 0;
    }

    /* remove edges with valid <= -1 from edge lists of all nodes */
    (void)GCedges(finalcnt,edge);
    if(VERB){
      if(bulgecnt > 0 || verb)
	printf("graph_bulge(%d,%0.1f): %d complex bulges found(%d edges deleted)(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n", BFSdepth, MaxBulge, bulgecnt,bulgefix, finalcnt,numtopnodes,wtime(),mtime());
      if(cyclecnt >= 0 || verb)
	printf("graph_bulge(%d,%0.1f): %d cycles found(%d edges deleted)(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n", BFSdepth, MaxBulge, cyclecnt, cyclefix, finalcnt,numtopnodes,wtime(),mtime());
      if(VERB>=2)
	printf("BFScnt=%d\n",BFScnt1);
      fflush(stdout);
    }

    if(VERB>=2){/* locate highest coverage edge */
      double maxCov = 0.0;
      Cedge *maxEdge = 0;
      for(int t=0; t < numtopnodes; t++){
	Cnode *Tnode = topnodes[t];
	for(int u = 0; u < Tnode->edgecnt; u++){
	  Cedge *TUedge = &edge[Tnode->edgeid[u]];
	  if(DEBUG) assert(TUedge->cov >= 0.0);
	  if(TUedge->cov > maxCov){
	    maxCov = TUedge->cov;
	    maxEdge = TUedge;
	  }
	}	
      }
      printf("graph_bulge(before complex bulges): Largest coverage = %0.2f : (N%d,N%d,%0.3f,%d,%d)\n",
	     maxCov, nodeorder[node[maxEdge->Lmap].mapid],nodeorder[node[maxEdge->Rmap].mapid],maxEdge->distance,
	     maxEdge->Lflipped,maxEdge->Rflipped);
      fflush(stdout);
    }

    if(DEBUG>=2){/* all nodes should have parent == 0, mark == SPmark = 0 */
      for(int t=0; t < numtopnodes; t++){
	Cnode *Tnode = topnodes[t];
	assert(Tnode->parent == 0);
	assert(Tnode->mark == 0);
	assert(Tnode->SPmark == 0);
	/* all edges should have valid == 1  && cov >= 0.0 */
	for(int u = 0; u < Tnode->edgecnt; u++){
	  Cedge *TUedge = &edge[Tnode->edgeid[u]];
	  assert(TUedge->valid == 1);
	  assert(TUedge->cov > -0.01);
	}
      }
    }

    /* Perform BFS based on number of hops, up to BFSdepth hops (or distance MaxBulge) to find 2 paths in same direction that converge (or a cycle within the same path) */
    bulgecnt = bulgefix = cyclecnt = cyclefix = 0;

    if(!PosEdgeInit)    /* update PosEdge of remaining topnodes */
      PosEdge();

    if(DEBUG>=2){/* check that node distances are correctly ordered */
      for(int i=0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	for(int Vflipped = 0;Vflipped <= 1; Vflipped++){/* try both orientations of V */      
	  for(int vw= (Vflipped ? 0 : Vnode->PosEdge); vw < (Vflipped ? Vnode->PosEdge : Vnode->edgecnt); vw++){
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];	
	    int NDsign = (VWedge->Lmap == Vnode->nodeid) ? (VWedge->Lflipped ? -1 : 1) : (VWedge->Rflipped ? 1 : -1);
	    if(DEBUG)assert(NDsign == (Vflipped ? -1 : 1));
	  }
	}
      }
    }

    if(BULGE_AVOID){/* mark all edges with mark=1 IFF the edge is the last forward/backward edge from a node */
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	for(int vw = Vnode->edgecnt; --vw >= 0;){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  VWedge->mark = 0;
	}
      }

      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	if(Vnode->PosEdge == 1){/* mark single backward edge as critical */
	  Cedge *VWedge = &edge[Vnode->edgeid[0]];
	  VWedge->mark = 1;
	}
	if(Vnode->edgecnt == Vnode->PosEdge+1){/* mark single forward edge as critical */
	  Cedge *VWedge = &edge[Vnode->edgeid[Vnode->edgecnt-1]];
	  VWedge->mark = 1;
	}
      }
    }

    /* allocate memory for bulge information for each Vnode = topnodes[i=0..numtopnodes-1] and Vfliped = 0,1 */
    CBulge *bulges = 0;
    int *NodeMem = 0;
    Cedge **EdgeMem = 0;
    int maxbulges = 0;
    int cntbulges = 0;
    if(BULGE_DEFER){
      maxbulges = numtopnodes*4;/* NOTE : If more bulges are found the excess will be ignored till the next call */
      bulges = new CBulge[maxbulges];/* bulges[i] has PathLen1 == 0 OR contains a graph bulge at Vnode == bulges[i].Vnode */
      NodeMem = new int[maxbulges*(BFSdepth*4+2)];
      EdgeMem = new Cedge*[maxbulges*(BFSdepth*2+1)];
      int NodeCnt = 0;
      int EdgeCnt = 0;
      for(int i = 0; i < maxbulges; i++){
	CBulge *p = &bulges[i];
	p->PathLen1 = 0;
	p->PathNode1 = &NodeMem[NodeCnt]; NodeCnt += BFSdepth+1;
	p->PathFlipped1 = &NodeMem[NodeCnt]; NodeCnt += BFSdepth+1;
	p->PathNode2 = &NodeMem[NodeCnt]; NodeCnt += BFSdepth;
	p->PathFlipped2 = &NodeMem[NodeCnt]; NodeCnt += BFSdepth;
	p->PathEdge1 = &EdgeMem[EdgeCnt]; EdgeCnt += BFSdepth+1;
	p->PathEdge2 = &EdgeMem[EdgeCnt]; EdgeCnt += BFSdepth;
      }
      if(DEBUG) assert(NodeCnt == maxbulges*(BFSdepth*4+2));
      if(DEBUG) assert(EdgeCnt == maxbulges*(BFSdepth*2+1));
    }

    #pragma omp parallel num_threads(nthreads) if(BULGE_PAR && nthreads > 1 && (BULGE_PAR >=2 || BFSdepth <= BULGE_FAST))
    {
      /* allocate per-thread memory */
      int tid = 0;
      #ifdef _OPENMP
      tid = omp_get_thread_num ();
      #endif

      BFSdata *BFS = new BFSdata[numnodes];
      CSPdata *SP = new CSPdata[numnodes];
      for(int i=0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	if(DEBUG>=2) assert(0 <= Vnode->nodeid  && Vnode->nodeid < numnodes && Vnode-node == Vnode->nodeid);
	DATA(BFS,Vnode,parent) = 0;
	DATA(BFS,Vnode,mark) = 0;
	DATA(BFS,Vnode,SPmark) = 0;
	DATA(SP,Vnode,SPmark) = 0;
	DATA(SP,Vnode,hindex) = Vnode->edgecnt;
	DATA(SP,Vnode,Enodeid) = Vnode->Enodeid;
      }

      int mycntbulges = 0;

      #pragma omp for schedule(dynamic,16)
      for(int i=0; i < numtopnodes; i++){

	if(BULGE_DEFER && mycntbulges >= maxbulges)
	  continue;//NOTE : cannot break out of multithreaded loop, so just skip all subsequent iterations

	Cnode *Vnode = topnodes[i];
	if(DEBUG>=2) assert(DATA(BFS,Vnode,parent) == 0);
	DATA(BFS,Vnode,Rdist) = 0.0;
	DATA(BFS,Vnode,mincov) = 1.1*HUGE_COV;

	for(int Vflipped = 0;Vflipped <= 1; Vflipped++){/* try both orientations of V */      
	  /* need 2 forward(backward if Vflipped) edges for start of complex bulge */

	  if(Vflipped ? (Vnode->PosEdge < 2) : (Vnode->edgecnt < Vnode->PosEdge + 2))
	    continue;

	  if(VERB>=2 && verb>=2){
	    #pragma omp critical
	    {
	      Vcnt++;
	      printf("tid=%d:i=%d/%d:Vnode=N%d%c:PosEdge=%d,edgecnt=%d(Vcnt=%d,BFScnt=%d)\n",
		     tid,i,numtopnodes,nodeorder[Vnode->mapid],Vflipped ? '\'':' ', Vnode->PosEdge,Vnode->edgecnt,Vcnt,BFScnt2);
	      fflush(stdout);
	    }
	  }

	  if(DEBUG>=2) assert(0 <= Vnode->nodeid && Vnode->nodeid < numnodes && Vnode-node == Vnode->nodeid);
	  DATA(BFS,Vnode,parent) = Vnode;
	  DATA(BFS,Vnode,pflip) = Vflipped;
	  for(int vw= (Vflipped ? 0 : Vnode->PosEdge); vw < (Vflipped ? Vnode->PosEdge : Vnode->edgecnt); vw++){
	    Cedge *VWedge = &edge[Vnode->edgeid[vw]];	
	    if(VWedge->valid < 0)
	      continue;
	    if(DEBUG>=2) {/* check that nodedistance has correct sign */
	      int NDsign = (VWedge->Lmap == Vnode->nodeid) ? (VWedge->Lflipped ? -1 : 1) : (VWedge->Rflipped ? 1 : -1);
	      assert(NDsign == (Vflipped ? -1 : 1));
	    }
	    Cnode *Wnode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
	    int Wflipped = (Vflipped ^ (VWedge->Lflipped | VWedge->Rflipped));

	    for(int vx= vw+1; vx < (Vflipped ? Vnode->PosEdge : Vnode->edgecnt); vx++){
	      Cedge *VXedge = &edge[Vnode->edgeid[vx]];
	      if(VXedge->valid < 0)
		continue;
	      if(DEBUG>=2) {/* check that nodedistance has correct sign */
		int NDsign = (VXedge->Lmap == Vnode->nodeid) ? (VXedge->Lflipped ? -1 : 1) : (VXedge->Rflipped ? 1 : -1);
		assert(NDsign == (Vflipped ? -1 : 1));
	      }
	      Cnode *Xnode = &node[(VXedge->Lmap == Vnode->nodeid) ? VXedge->Rmap : VXedge->Lmap];
	      int Xflipped = (Vflipped ^ (VXedge->Lflipped ^ VXedge->Rflipped));
	      if(Xnode == Wnode)/* simple bulge : handled seperately */
		continue;
	  
	      if(VERB>=2 && verb>=2){
		#pragma omp critical
		{
		  printf("tid=%d:Vnode=N%d%c,edgecnt=%d,PosEdge=%d,Wnode=N%d%c,Xnode=N%d%c,vw=%d,VWedge=(id=%d,N%d,N%d,cov=%0.3f,f%d,f%d,dist=%0.3f),vx=%d,VXedge=(id=%d,N%d,N%d,cov=%0.3f,f%d,f%d,dist=%0.3f)\n",
			 tid,nodeorder[Vnode->mapid],Vflipped ? '\'':' ',Vnode->edgecnt,Vnode->PosEdge,
			 nodeorder[Wnode->mapid],Wflipped ? '\'':' ',
			 nodeorder[Xnode->mapid],Xflipped ? '\'':' ',
			 vw,VWedge->edgeid,nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->cov,VWedge->Lflipped,VWedge->Rflipped,VWedge->distance,
			 vx,VXedge->edgeid,nodeorder[node[VXedge->Lmap].mapid],nodeorder[node[VXedge->Rmap].mapid],VXedge->cov,VXedge->Lflipped,VXedge->Rflipped,VXedge->distance);
		  fflush(stdout);
		}
	      }

	      /* try to find two paths starting from Wnode and Xnode that lead to a common node */
	      /* need to color nodes with 2 colors (mark=1 and mark=2) and don't count cases 
		 where the same color paths converge */
	      int found = 0;
	      if(DEBUG>=2)assert(DATA(BFS,Wnode,parent) == 0);
	      if(DEBUG>=2) assert(0 <= Wnode->nodeid && Wnode->nodeid < numnodes && Wnode - node == Wnode->nodeid);
	      if(DEBUG>=2) assert(0 <= Vnode->nodeid && Vnode->nodeid < numnodes && Vnode - node == Vnode->nodeid);
	      DATA(BFS,Wnode,parent) = Vnode;
	      DATA(BFS,Wnode,pedge) = VWedge;
	      DATA(BFS,Wnode,Rdist) = VWedge->distance;
	      DATA(BFS,Wnode,mincov) = VWedge->cov;
	      DATA(BFS,Wnode,mark) = 1;
	      DATA(BFS,Wnode,pflip) = Wflipped;

	      if(DEBUG>=2) assert(DATA(BFS,Xnode,parent) == 0);
	      if(DEBUG>=2) assert(0 <= Xnode->nodeid && Xnode->nodeid < numnodes && Xnode - node == Xnode->nodeid);
	      DATA(BFS,Xnode,parent) = Vnode;
	      DATA(BFS,Xnode,pedge) = VXedge;
	      DATA(BFS,Xnode,Rdist) = VXedge->distance;
	      DATA(BFS,Xnode,mincov) = VXedge->cov;
	      DATA(BFS,Xnode,mark) = 2;
	      DATA(BFS,Xnode,pflip) = Xflipped;

	      priority_queue< Cpath, vector<Cpath> , less<Cpath> > BFSqueue;/* Note less (<) is implemented as greater (>) in Cpath! */
	      stack<Cnode*> FinishedNodes;/* stack of nodes removed from BFSqueue : used to reset parent field later */
	      BFSqueue.push(Cpath(Wnode,Wflipped,1.0));
	      if(VERB>=2 && verb>=2){
                #pragma omp critical
		{
		  printf("tid=%d:queue.push:priority=%0.1f,Wnode=N%d%c,parent=N%d%c,mincov=%0.1f,Rdist=%0.3f,Rwt=%0.3f,mark=%d\n",
			 tid,1.0,nodeorder[Wnode->mapid],Wflipped ? '\'':' ',nodeorder[Vnode->mapid],Vflipped ? '\'':' ',
			 DATA(BFS,Wnode,mincov),DATA(BFS,Wnode,Rdist),0.0, DATA(BFS,Wnode,mark));
		  fflush(stdout);
		}
	      }
	      BFSqueue.push(Cpath(Xnode,Xflipped,1.0));
	      if(VERB>=2 && verb>=2){
                #pragma omp critical
		{
		  printf("tid=%d:queue.push:priority=%0.1f,Xnode=N%d%c,parent=N%d%c,mincov=%0.1f,Rdist=%0.3f,Rwt=%0.3f,mark=%d\n",
			 tid,1.0,nodeorder[Xnode->mapid],Xflipped ? '\'':' ',nodeorder[Vnode->mapid],Vflipped ? '\'':' ',
			 DATA(BFS,Xnode,mincov),DATA(BFS,Xnode,Rdist), 0.0,DATA(BFS,Xnode,mark));
		  fflush(stdout);
		}
	      }
	      while(BFSqueue.size() > 0){
		Cpath top = BFSqueue.top(); BFSqueue.pop();
		FinishedNodes.push(top.node);
		if(VERB>=2) BFScnt2++;
		if(top.priority > BFSdepth)
		  break;
		Cnode *Pnode = top.node;
		if(DATA(BFS,Pnode,Rdist) > MaxBulge)
		  continue;
		int Pflipped = top.flipped;
		if(DEBUG>=2) assert(DATA(BFS,Pnode,parent));
		if(VERB>=2 && verb>=2){
                  #pragma omp critical
		  {
		    printf("tid=%d:queue.pop():priority=%0.1f,Pnode=N%d%c,parent=N%d%c,mincov=%0.1f,Rdist=%0.3f,Rwt=%0.3f,mark=%d\n",
			   tid,top.priority,nodeorder[Pnode->mapid],Pflipped ? '\'':' ',nodeorder[DATA(BFS,Pnode,parent)->mapid],DATA(BFS,DATA(BFS,Pnode,parent),pflip) ? '\'':' ',
			   DATA(BFS,Pnode,mincov),DATA(BFS,Pnode,Rdist),0.0,DATA(BFS,Pnode,mark));
		    fflush(stdout);
		  }
		}

		/* If Pflipped==1: look for neighbors of Pnode with -ve nodedistance 
		   If Pflipped==0: look for neighbors of Pnode with +ve node distance */
		int pn;
		for(pn = (Pflipped ? 0 : Pnode->PosEdge); pn < (Pflipped ? Pnode->PosEdge : Pnode->edgecnt); pn++){
		  Cedge *PNedge = &edge[Pnode->edgeid[pn]];
		  if(PNedge->valid < 0)
		    continue;
		  Cnode *Nnode = &node[(PNedge->Lmap == Pnode->nodeid) ? PNedge->Rmap : PNedge->Lmap];
		  int Nflipped = (Pflipped ^ (PNedge->Lflipped | PNedge->Rflipped));
		  if(VERB>=2 && verb>=2)
		    printf("pn=%d:Pnode=N%d%c,PosEdge=%d,edgecnt=%d),Nnode=N%d%c:valid=%d,parent=N%d%c,mark=%d\n",
			   pn,nodeorder[Pnode->mapid],Pflipped?'\'':' ',Pnode->PosEdge,Pnode->edgecnt,nodeorder[Nnode->mapid],Nflipped ? '\'':' ',
			   PNedge->valid, DATA(BFS,Nnode,parent) ? nodeorder[DATA(BFS,Nnode,parent)->mapid] : -1, (DATA(BFS,Nnode,parent) && DATA(BFS,DATA(BFS,Nnode,parent),pflip)) ? '\'':' ',
			   DATA(BFS,Nnode,mark));
		  double distance2 = DATA(BFS,Pnode,Rdist) + PNedge->distance;
		  double mincov2 = min(DATA(BFS,Pnode,mincov),PNedge->cov);

		  if(DEBUG>=2) {/* check that nodedistance has correct sign */
		    int NDsign = (PNedge->Lmap == Pnode->nodeid) ? (PNedge->Lflipped ? -1 : 1) : (PNedge->Rflipped ? 1 : -1);
		    if(DEBUG && !(NDsign == (Pflipped ? -1 : 1))){
		      #pragma omp critical
		      {		      
			printf("tid=%d:pn=%d:Pnode=N%d%c,PosEdge=%d,edgecnt=%d,PNedge=(id=%d,N%d,N%d,cov=%0.3f,f%d,f%d,distance=%0.3f),Nnode=N%d%c:valid=%d,parent=N%d%c,mark=%d:mincov->%0.3f,Rdist->%0.3f\n",
			       tid,pn,nodeorder[Pnode->mapid],Pflipped?'\'':' ',Pnode->PosEdge,Pnode->edgecnt,PNedge->edgeid,nodeorder[node[PNedge->Lmap].mapid],nodeorder[node[PNedge->Rmap].mapid],
			       PNedge->cov,PNedge->Lflipped,PNedge->Rflipped,PNedge->distance,nodeorder[Nnode->mapid],Nflipped ? '\'':' ',
			       PNedge->valid, DATA(BFS,Nnode,parent) ? nodeorder[DATA(BFS,Nnode,parent)->mapid] : -1, (DATA(BFS,Nnode,parent) && DATA(BFS,DATA(BFS,Nnode,parent),pflip)) ? '\'':' ',
			       DATA(BFS,Nnode,mark),mincov2,distance2);
			printf("NDsign = %d, Pflipped=%d\n",NDsign,Pflipped);
			fflush(stdout);

			assert(NDsign == (Pflipped ? -1 : 1));
		      }
		    }
		  }

		  if(VERB>=2 && verb>=2){
                    #pragma omp critical
		    {
		      printf("tid=%d:pn=%d:Pnode=N%d%c,PosEdge=%d,edgecnt=%d,PNedge=(id=%d,N%d,N%d,cov=%0.3f,f%d,f%d,distance=%0.3f),Nnode=N%d%c:valid=%d,parent=N%d%c,mark=%d:mincov->%0.3f,Rdist->%0.3f\n",
			     tid,pn,nodeorder[Pnode->mapid],Pflipped?'\'':' ',Pnode->PosEdge,Pnode->edgecnt,PNedge->edgeid,nodeorder[node[PNedge->Lmap].mapid],nodeorder[node[PNedge->Rmap].mapid],
			     PNedge->cov,PNedge->Lflipped,PNedge->Rflipped,PNedge->distance,nodeorder[Nnode->mapid],Nflipped ? '\'':' ',
			     PNedge->valid, DATA(BFS,Nnode,parent) ? nodeorder[DATA(BFS,Nnode,parent)->mapid] : -1, (DATA(BFS,Nnode,parent) && DATA(BFS,DATA(BFS,Nnode,parent),pflip)) ? '\'':' ',
			     DATA(BFS,Nnode,mark),mincov2,distance2);
		      fflush(stdout);
		    }
		  }


		  if(DATA(BFS,Nnode,parent)){/* found two paths from Vnode to Nnode */
		    if(DEBUG>=2) assert(node <= DATA(BFS,Nnode,parent) && DATA(BFS,Nnode,parent) < &node[numnodes]);
		    if(Nnode == Vnode){/* a cycle : try to break the weakest link(s) */
		      if(VERB>=2 && verb >=2){
			#pragma omp critical
			{
			  printf("tid=%d:Found cycle from Vnode=N%d(flip=%d) to Pnode=N%d(flip=%d) to Vnode=N%d(flip=%d):V to P:mincov=%0.3f,Rdist=%0.3f:P to V:cov=%0.3f,dist=%0.3f\n",
				 tid,nodeorder[Vnode->mapid],Vflipped,nodeorder[Pnode->mapid],Pflipped,nodeorder[Vnode->mapid],Nflipped,DATA(BFS,Pnode,mincov),DATA(BFS,Pnode,Rdist),
				 PNedge->cov,PNedge->distance);
			  fflush(stdout);
			}
		      }
		      if(DEBUG>=2 && ! (DATA(BFS,Pnode,pedge) != PNedge)){
			#pragma omp critical
			{
			  printf("tid=%d:Vnode=N%d(flip=%d)\n",tid,nodeorder[Vnode->mapid],Vflipped);
			  fflush(stdout);
			  assert(DATA(BFS,Pnode,pedge) != PNedge);
			}
		      }
		      if(cyclebreak(Vnode,Vflipped,Pnode,Pflipped,PNedge,Nflipped,Cycles,cntCycles,maxCycles,BFSdepth,verb,cyclecnt,cyclefix,BFS)){
			found = 1;/* force termination of BFS for (Vnode,Wnode,Xnode) later */
			break;/* terminate BFS */
		      }
		      continue;
		    }
		    if(DATA(BFS,Nnode,mark) == DATA(BFS,Pnode,mark))
		      continue;/* self-bulge (will be found later when Vnode is the common ancestor of Nnode and Pnode) */
		  
		    if(DATA(BFS,Nnode,Rdist) + DATA(BFS,Pnode,Rdist) + PNedge->distance > MaxBulge){
		      if(VERB>=2 && verb){
			printf("tid=%d:Nnode->Rdist=%0.3f,Pnode->Rdist=%0.3f,PNedge->distance=%0.3f: Total distance exceeds MaxBulge=%0.3f\n",
			       tid,DATA(BFS,Nnode,Rdist),DATA(BFS,Pnode,Rdist),PNedge->distance,MaxBulge);
			fflush(stdout);
		      }
		      continue;
		    }

		    found = 1;/* force termination of BFS for (Vnode,Wnode,Xnode) later */

		    double var = 2.0*SM*SM + SD[0]*SD[0]*(distance2 + DATA(BFS,Nnode,Rdist));
		    double err = distance2 - DATA(BFS,Nnode,Rdist);
		    if(VERB>=2 && verb){
                      #pragma omp critical
		      {
			printf("tid=%d:Found two paths from Vnode=N%d%c to Nnode=N%d (Pnode=N%d%c)\nPath1:distance=%8.3f,cov=%5.3f,mark=%d: N%d%c(cov=%0.3f)",
			       tid,nodeorder[Vnode->mapid], Vflipped ? '\'':' ',nodeorder[Nnode->mapid],nodeorder[Pnode->mapid],Pflipped?'\'':' ',
			       distance2, mincov2, Pnode->mark, nodeorder[Nnode->mapid], Nflipped ? '\'' : ' ', PNedge->cov);
			for(Cnode *Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			  assert(Tnode != 0);
			  printf(" <- N%d%c(cov=%5.1f)", nodeorder[Tnode->mapid], DATA(BFS,Tnode,pflip) ? '\'' : ' ', DATA(BFS,Tnode,pedge)->cov);
			}
			printf(" <- N%d%c\nPath2:distance=%8.3f,cov=%5.3f,mark=%d: N%d%c(cov=%0.3f)", nodeorder[Vnode->mapid], Vflipped ? '\'':' ',
			       DATA(BFS,Nnode,Rdist), DATA(BFS,Nnode,mincov), DATA(BFS,Nnode,mark),nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'' : ' ', DATA(BFS,Nnode,pedge)->cov);
			for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			  if(DEBUG>=2 && Tnode == 0){
			    printf("Tnode = 0!\n");
			    fflush(stdout);
			    printf("Vnode=N%d, Nnode=N%d%c, Pnode=N%d%c\n", nodeorder[Vnode->mapid],
				   nodeorder[Nnode->mapid],Nflipped ? '\'':' ',
				   nodeorder[Pnode->mapid],Pflipped ? '\'':' ');
			    fflush(stdout);
			    assert(Tnode != 0);
			  }
			  printf(" <- N%d%c(cov=%5.3f)", nodeorder[Tnode->mapid], DATA(BFS,Tnode,pflip) ? '\'' : ' ', DATA(BFS,Tnode,pedge)->cov);
			}
			printf(" <- N%d%c\n", nodeorder[Vnode->mapid],Vflipped ? '\'':' ');
			printf("path distance err=%0.3f,sd=%0.3f\n",fabs(err),sqrt(var));
			fflush(stdout);
		      }
		    }

		    Cedge *newPNedge;
		    Cnode *newPnode;
		    int newPflipped;

		    /* check if paths are compatible (orientation and distance) */
		    int compatible = 1;
		    if(DATA(BFS,Nnode,pflip) != Nflipped)
		      compatible = 0;
		    else if(err*err >= var*(NMAXERR*NMAXERR))
		      compatible = 0;

		    if(BULGE_FAST && (BFSdepth <= BULGE_FAST || compatible)){/* just use original paths */
		      newPNedge = PNedge;
		      newPnode = Pnode;
		      newPflipped = Pflipped;
		    } else {
		      /* First try to find the strongest (highest coverage) version of each path, so we keep the strongest possible path */

		      /* Try to strengthen the path via P and PNedge (1st path displayed, but using mincov2,cov2 etc) */
		      for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent))
			DATA(SP,Tnode,SPmark) |= 16; /* mark the internal nodes of the other path with SPmark = 16 */

		      if(DEBUG>=2) assert(DATA(BFS,Nnode,pedge)->valid == 1);
		      if(!BULGE_PAR || nthreads <= 1)/* the pedge->valid field is shared by all threads so cannot be changed without a critical section */
			DATA(BFS,Nnode,pedge)->valid = 0;

		      (void) ShortestPaths(Vnode->nodeid,Vflipped,0,1,0.0,0.0,SP);

		      if(!BULGE_PAR || nthreads <= 1)
			DATA(BFS,Nnode,pedge)->valid = 1;
		      else if(DATA(SP,Nnode,SPedge)[0][Nflipped] == DATA(BFS,Nnode,pedge)){/* check if new best path includes Nnode->pedge */
			#pragma omp critical
			{			/* repeat ShortestPaths() in critical section */
			  DATA(BFS,Nnode,pedge)->valid = 0;
			  (void) ShortestPaths(Vnode->nodeid,Vflipped,0,1,0.0,0.0,SP);
			  DATA(BFS,Nnode,pedge)->valid = 1;
			}
		      }

		      for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent))
			DATA(SP,Tnode,SPmark) &= ~16; /* reset SPmark */

		      if(DEBUG>=2 && DATA(SP,Nnode,SPcost)[0][Nflipped] < 0.0){
			printf("Failed to find shortest path from Vnode=N%d%c to Nnode=N%d%c (updated 1st path)!\n",
			       nodeorder[Vnode->mapid],Vflipped ? '\'':' ',nodeorder[Nnode->mapid],Nflipped ? '\'':' ');
			printf("Nnode->SPcost[0][Nflipped] = %0.3f\n", DATA(SP,Nnode,SPcost)[0][Nflipped]);
			printf("Nnode->pedge=(N%d,N%d,%0.3f,%d,%d)\n",
			       nodeorder[node[DATA(BFS,Nnode,pedge)->Lmap].mapid],nodeorder[node[DATA(BFS,Nnode,pedge)->Rmap].mapid],
			       DATA(BFS,Nnode,pedge)->distance, DATA(BFS,Nnode,pedge)->Lflipped, DATA(BFS,Nnode,pedge)->Rflipped);

			#pragma omp critical
			{
			  printf("Calling ShortestPaths(Vnode,Vflipped,0,1) again with DFSverb=1\n");
			  fflush(stdout);
			  DFSverb = 1;

			  for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent))
			    DATA(SP,Tnode,SPmark) |= 16; /* mark the internal nodes of the other path with SPmark = 16 */
			  DATA(BFS,Nnode,pedge)->valid = 0;
			  (void) ShortestPaths(Vnode->nodeid,Vflipped,0,1,0.0,0.0,SP);
			  printf("Nnode->SPcost[0][Nflipped] = %0.3f\n", DATA(SP,Nnode,SPcost)[0][Nflipped]);
			  fflush(stdout);
			  assert(DATA(SP,Nnode,SPcost)[0][Nflipped] >= 0.0);
			}
			exit(1);
		      }

		      Cnodeflip backflip = loopback(Nnode->nodeid,Nflipped,-1,0,SP);
		      if(backflip.nodeid >= 0){
			if(VERB>=2 && verb){
			  Cnode *Bnode = &node[backflip.nodeid];
			  printf("Best alternate to 1st path traverses N%d%c twice : using original path\n",
				 nodeorder[Bnode->mapid],backflip.flip ? '\'':' ');
			  fflush(stdout);
			}
			newPNedge = PNedge;
			newPnode = Pnode;
			newPflipped = Pflipped;
		      } else {
			/* now replace the path PNedge, Pnode ... Vnode : only the fields parent,pedge,pflip,Rdist are updated and the overal mincov2 computed 
			   (Rdwt and mincov fields are left unchanged and will be invalid) */
			newPNedge = DATA(SP,Nnode,SPedge)[0][Nflipped];
			newPnode = &node[(newPNedge->Lmap == Nnode->nodeid) ? newPNedge->Rmap : newPNedge->Lmap];
			if(DEBUG>=2) assert(newPnode != Vnode);
			newPflipped = (Nflipped ^ (newPNedge->Lflipped | newPNedge->Rflipped));
			if(VERB>=2 && verb)
			  printf("Found best alternate for 1st path from N%d%c(cov=%0.3f) <- N%d%c",
				 nodeorder[Nnode->mapid],Nflipped ? '\'':' ', newPNedge->cov,
				 nodeorder[newPnode->mapid],newPflipped ? '\'':' ');

			mincov2 = newPNedge->cov;
			Cnode *Tnode = newPnode;
			int Tflipped = newPflipped;
			for(; Tnode != Vnode;){
			  Cedge *TUedge = DATA(SP,Tnode,SPedge)[0][Tflipped];
			  Cnode *Unode = &node[(TUedge->Lmap == Tnode->nodeid) ? TUedge->Rmap : TUedge->Lmap];
			  int Uflipped = (Tflipped ^ (TUedge->Lflipped | TUedge->Rflipped));
			  if(!DATA(BFS,Tnode,parent))
			    FinishedNodes.push(Tnode);
			  DATA(BFS,Tnode,pedge) = TUedge;
			  if(DEBUG>=2) assert(0 <= Tnode->nodeid && Tnode->nodeid < numnodes && Tnode-node == Tnode->nodeid);
			  if(DEBUG>=2) assert(0 <= Unode->nodeid && Unode->nodeid < numnodes && Unode-node == Unode->nodeid);
			  DATA(BFS,Tnode,parent) = Unode;
			  DATA(BFS,Tnode,pflip) = Tflipped;
			  DATA(BFS,Tnode,Rdist) = DATA(SP,Tnode,SPdist)[0][Tflipped];
			  if(TUedge->cov < mincov2)
			    mincov2 = TUedge->cov;

			  Tnode = Unode;
			  Tflipped = Uflipped;
			  if(VERB>=2 && verb)
			    printf("(cov=%0.3f) <- N%d%c", TUedge->cov, nodeorder[Tnode->mapid], Tflipped ? '\'':' ');
			}
			if(DEBUG>=2)assert(Tflipped == Vflipped);
			distance2 = DATA(BFS,newPnode,Rdist) + newPNedge->distance;
			if(VERB>=2 && verb)
			  printf(" (mincov=%0.3f,distance=%0.3f)\n", mincov2, distance2);

			/* Try to strengthen the path via N->parent (2nd path displayed) */
			for(Cnode *Tnode = newPnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent))
			  DATA(SP,Tnode,SPmark) = 16; /* mark the internal nodes of the other path with SPmark = 16 */

			if(DEBUG>=2) assert(newPNedge->valid == 1);
			if(!BULGE_PAR || nthreads <= 1)/* the valid field is shared by all threads so cannot be changed without a critical section */
			  newPNedge->valid = 0;
			(void)ShortestPaths(Vnode->nodeid,Vflipped,0,1,0.0,0.0,SP);
			if(!BULGE_PAR || nthreads <= 1)
			  newPNedge->valid = 1;
			else if(DATA(SP,Nnode,SPedge)[0][DATA(BFS,Nnode,pflip)] == newPNedge){
			  #pragma omp critical
			  {
			    newPNedge->valid = 0;
			    (void)ShortestPaths(Vnode->nodeid,Vflipped,0,1,0.0,0.0,SP);
			    newPNedge->valid = 1;
			  }
			}

			for(Cnode *Tnode = newPnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent))
			  DATA(SP,Tnode,SPmark) = 0;

			if(DEBUG>=2 && DATA(SP,Nnode,SPcost)[0][DATA(BFS,Nnode,pflip)] < 0.0){
			  printf("Failed to find shortest path from Vnode=N%d%c to Nnode=N%d%c (updated 2nd path)!\n",
				 nodeorder[Vnode->mapid],Vflipped ? '\'':' ',nodeorder[Nnode->mapid],DATA(BFS,Nnode,pflip) ? '\'':' ');
			  printf(" Nnode->SPcost[0][Nnode->pflip] = %0.3f\n",DATA(SP,Nnode,SPcost)[0][DATA(BFS,Nnode,pflip)]);
			  printf("Calling ShortestPaths(Vnode,Vflipped,0,1) again with DFSverb=1\n");
			  fflush(stdout);
			  DFSverb = 1;
			  for(Cnode *Tnode = newPnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent))
			    DATA(SP,Tnode,SPmark) = 16; /* mark the internal nodes of the other path with SPmark = 16 */
			  #pragma omp critical
			  {
			    newPNedge->valid = 0;
			    (void)ShortestPaths(Vnode->nodeid,Vflipped,0,1,0.0,0.0,SP);
			    printf(" Nnode->SPcost[0][Nnode->pflip] = %0.3f\n",DATA(SP,Nnode,SPcost)[0][DATA(BFS,Nnode,pflip)]);
			    fflush(stdout);
			    assert(DATA(SP,Nnode,SPcost)[0][DATA(BFS,Nnode,pflip)] >= 0.0);
			  }
			  exit(1);
			}
			backflip = loopback(Nnode->nodeid,DATA(BFS,Nnode,pflip),-1,0,SP);
			if(backflip.nodeid >= 0){
			  if(VERB>=2 && verb){
			    Cnode *Bnode = &node[backflip.nodeid];
			    printf("Best alternate to 2nd path traverses N%d%c twice : using original path\n",
				   nodeorder[Bnode->mapid],backflip.flip ? '\'':' ');
			    fflush(stdout);
			  }
			} else {
			  /* now replace the path Nnode->parent ... Vnode : only the fields parent,pedge,pflip,Rdist are updated and the overal mincov1 computed 
			     (Rdwt and mincov fields are left unchanged and will be invalid) */
			  if(VERB>=2 && verb)
			    printf("Found best alternate for 2nd path from N%d%c",
				   nodeorder[Nnode->mapid],DATA(BFS,Nnode,pflip) ? '\'':' ');
			  double mincov1 = 1.1*HUGE_COV;
			  Cnode *Tnode = Nnode;
			  int Tflipped = DATA(BFS,Nnode,pflip);
			  for(;Tnode != Vnode;){
			    Cedge *TUedge = DATA(SP,Tnode,SPedge)[0][Tflipped];
			    Cnode *Unode = &node[(TUedge->Lmap == Tnode->nodeid) ? TUedge->Rmap : TUedge->Lmap];
			    int Uflipped = (Tflipped ^ (TUedge->Lflipped | TUedge->Rflipped));
			    if(!DATA(BFS,Tnode,parent))
			      FinishedNodes.push(Tnode);
			    DATA(BFS,Tnode,pedge) = TUedge;
			    if(DEBUG>=2) assert(0 <= Tnode->nodeid && Tnode->nodeid < numnodes && Tnode-node == Tnode->nodeid);
			    if(DEBUG>=2) assert(0 <= Unode->nodeid && Unode->nodeid < numnodes && Unode-node == Unode->nodeid);
			    DATA(BFS,Tnode,parent) = Unode;
			    DATA(BFS,Tnode,pflip) = Tflipped;
			    DATA(BFS,Tnode,Rdist) = DATA(SP,Tnode,SPdist)[0][Tflipped];
			    if(TUedge->cov < mincov1)
			      mincov1 = TUedge->cov;
			    Tnode = Unode;
			    Tflipped = Uflipped;

			    if(VERB>=2 && verb)
			      printf("(cov=%0.3f) <- N%d%c", TUedge->cov, nodeorder[Tnode->mapid], Tflipped ? '\'':' ');
			  }
			  if(DEBUG>=2)assert(Tflipped == Vflipped);
			  if(VERB>=2 && verb)
			    printf(" (mincov=%0.3f,distance=%0.3f)\n", mincov1,DATA(BFS,Nnode,Rdist));
			  DATA(BFS,Nnode,mincov) = mincov1;
			}		

			/* recompute path errors */
			var = 2.0*SM*SM + SD[0]*SD[0]*(distance2 + DATA(BFS,Nnode,Rdist));
			err = distance2 - DATA(BFS,Nnode,Rdist);
			if(VERB>=2 && verb){
			  printf("updated path distance err=%0.3f,sd=%0.3f\n",fabs(err),sqrt(var));
			  fflush(stdout);
			}
		      }
		    }

		    if(BULGE_PAR && !(mincov2 > 0.0 && DATA(BFS,Nnode,mincov) > 0.0)){/* another thread already took care of this bulge */
      		      if(DEBUG>=1+RELEASE) assert(!BULGE_DEFER);
		      if(DEBUG>=2){
			int valid = 1;
			if(!(DATA(BFS,Nnode,pedge)->valid > 0))
			  valid = 0;
			if(valid)
			  for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent))
			    if(!(DATA(BFS,Tnode,pedge)->valid > 0)){
			      valid = 0;
			      break;
			    }
			if(!(newPNedge->valid > 0))
			  valid = 0;
			if(valid)
			  for(Cnode *Tnode = newPnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent))
			    if(!(DATA(BFS,Tnode,pedge)->valid > 0)){
			      valid = 0;
			      break;
			    }
			assert(!valid);
			assert(found);
		      }
		      continue;
		    }

		    /* check if orientation of two paths are compatible */
		    if(DEBUG>=2) assert(mincov2 > 0.0);
		    if(DEBUG>=2) assert(DATA(BFS,Nnode,mincov) > 0.0);
		    double cov1=0.0,cov2=0.0;

		    if(DATA(BFS,Nnode,pflip) != Nflipped) {  /* paths have incompatible orientation */
		      /* check if distance of two paths are compatible */
		      if(MAXBULGE && distance2 + DATA(BFS,Nnode,Rdist) > MaxBulge){
			if(VERB>=2 && verb)
			  printf("  V=%d(flip=%d):(paths have incompatible orientation AND paths are too long (ignoring):distance=%0.3f+%0.3f,MaxBulge=%0.3f)\n",
				 Vnode->nodeid,Vflipped,distance2,DATA(BFS,Nnode,Rdist),MaxBulge);
		      } else {/* just delete the weakest link (even if not significantly weaker) provided cov < BulgeBreakMaxCov */
			if(min(mincov2,DATA(BFS,Nnode,mincov)) > BulgeBreakMaxCov){
			  if(VERB>=2 && verb)
			    printf("  (paths have incompatible orientation but weakest path is too strong to break:minvo2=%0.1f,Nnode->mincov=%0.1f)\n",mincov2,DATA(BFS,Nnode,mincov));
			} else if(mincov2 <= DATA(BFS,Nnode,mincov)){
			  cov1 = 0;
			  cov2 = -mincov2;
			  if(VERB>=2 && verb)
			    printf("  (paths have incompatible orientation: breaking 1st path:cov1=%0.1f,cov2=%0.1f)\n",cov1,cov2);
			} else if(DATA(BFS,Nnode,mincov) < mincov2){
			  cov1 = -DATA(BFS,Nnode,mincov);
			  cov2 = 0;
			  if(VERB>=2 && verb)
			    printf("  (paths have incompatible orientation: breaking 2nd path:cov1=%0.1f,cov2=%0.1f)\n",cov1,cov2);
			}
		      }
		    } else if(err*err >= var*(NMAXERR*NMAXERR)){/* just delete the weakest link of the weaker path (if it is significantly weaker than other path AND not stronger than BulgeBreakMaxCov) */
		      if(MAXBULGE && distance2 + DATA(BFS,Nnode,Rdist) > MaxBulge){
			if(VERB>=2 && verb)
			  printf("   V=%d(flip=%d):(paths do not match AND paths are too long (ignoring): distance=%0.3f+%0.3f,MaxBulge=%0.3f)\n",
				 Vnode->nodeid,Vflipped,distance2,DATA(BFS,Nnode,Rdist),MaxBulge);			       
		      } else {
			/* NOTE : we should check that the weakest link is unique : otherwise it is not clear which link needs to be deleted */
			if(min(mincov2,DATA(BFS,Nnode,mincov)) > BulgeBreakMaxCov){
			  if(VERB>=2 && verb)
			    printf("  V=%d(flip=%d):(paths do not match AND weakest path is too strong to break: mincov2=%0.1f,Nnode->mincov=%0.1f)\n",Vnode->nodeid,Vflipped,mincov2,DATA(BFS,Nnode,mincov));
			} else if(mincov2 < DATA(BFS,Nnode,mincov)*MinCoverageRatio){
			  cov1 = 0;
			  cov2 = -mincov2;
			  if(VERB>=2 && verb)
			    printf("  (paths do not match: breaking 1st path:cov1=%0.1f,cov2=%0.1f)\n",cov1,cov2);
			} else if(DATA(BFS,Nnode,mincov) < mincov2*MinCoverageRatio){
			  cov1 = -DATA(BFS,Nnode,mincov);
			  cov2 = 0;
			  if(VERB>=2 && verb)
			    printf("  (paths do not match: breaking 2nd path:cov1=%0.1f,cov2=%0.1f)\n",cov1,cov2);
			} else {
			  if(BREAK_BULGE){
			    if(mincov2 > DATA(BFS,Nnode,mincov)){
			      if(mincov2*MinCoverageRatio > DATA(BFS,Nnode,mincov)){
				cov1 = -DATA(BFS,Nnode,mincov);
				cov2 = 0;
				if(VERB>=2 && verb)
				  printf("  (paths do not match: breaking 2nd path:cov1=%0.1f,cov2=%0.1f)\n",cov1,cov2);
			      }
			    } else {
			      if(mincov2 < DATA(BFS,Nnode,mincov) * MinCoverageRatio){
				cov1 = 0;
				cov2 = -mincov2;
				if(VERB>=2 && verb)
				  printf("  (paths do not match: breaking 1st path:cov1=%0.1f,cov2=%0.1f)\n",cov1,cov2);
			      }
			    }
			  } else {
			    if(VERB>=2 && verb)
			      printf("  V=%d(flip=%d):(paths do not match but neither is clearly weaker:mincov2=%0.1f,Nnode->mincov=%0.1f)\n",Vnode->nodeid,Vflipped,mincov2,DATA(BFS,Nnode,mincov));
			  }
			}
		      }
		    } else {/* add the coverage of the weakest link of the weaker path to the other path and delete the weakest link */
		      double mincov1 = DATA(BFS,Nnode,mincov);
		      if(BULGE_AVOID){ /* determine weakest link excluding critical edges that will be protected (to ensure at least one non-critical link gets deleted) */
			double minC1 = 1.1*HUGE_COV, minC2 = 1.1*HUGE_COV;/* estimates of mincov1,mincov1 igoring critical edges */
			for(Cnode *Tnode = Nnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			  Cedge *Tpedge = DATA(BFS,Tnode,pedge);
			  if(!Tpedge->mark && Tpedge->cov < minC1)
			    minC1 = Tpedge->cov;
			}
			if(DEBUG && !(minC1 >= mincov1 - 1e-6)){
			  #pragma omp critical
			  {
			    printf("tid=%d:Found two paths from Vnode=N%d%c to Nnode=N%d (Pnode=N%d%c)\nPath1:distance=%8.3f,cov=%5.3f,mark=%d: N%d%c(cov=%0.3f)",
				   tid,nodeorder[Vnode->mapid], Vflipped ? '\'':' ',nodeorder[Nnode->mapid],nodeorder[Pnode->mapid],Pflipped?'\'':' ',
				   distance2, mincov2, Pnode->mark, nodeorder[Nnode->mapid], Nflipped ? '\'' : ' ', PNedge->cov);
			    for(Cnode *Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      assert(Tnode != 0);
			      printf(" <- N%d%c(cov=%5.1f)", nodeorder[Tnode->mapid], DATA(BFS,Tnode,pflip) ? '\'' : ' ', DATA(BFS,Tnode,pedge)->cov);
			    }
			    printf(" <- N%d%c\nPath2:distance=%8.3f,cov=%5.3f,mark=%d: N%d%c(cov=%0.3f)", nodeorder[Vnode->mapid], Vflipped ? '\'':' ',
				   DATA(BFS,Nnode,Rdist), DATA(BFS,Nnode,mincov), DATA(BFS,Nnode,mark),nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'' : ' ', DATA(BFS,Nnode,pedge)->cov);
			    for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      if(DEBUG>=2 && Tnode == 0){
				printf("Tnode = 0!\n");
				fflush(stdout);
				printf("Vnode=N%d, Nnode=N%d%c, Pnode=N%d%c\n", nodeorder[Vnode->mapid],
				       nodeorder[Nnode->mapid],Nflipped ? '\'':' ',
				       nodeorder[Pnode->mapid],Pflipped ? '\'':' ');
				fflush(stdout);
				assert(Tnode != 0);
			      }
			      printf(" <- N%d%c(cov=%5.3f)", nodeorder[Tnode->mapid], DATA(BFS,Tnode,pflip) ? '\'' : ' ', DATA(BFS,Tnode,pedge)->cov);
			    }
			    printf(" <- N%d%c\n", nodeorder[Vnode->mapid],Vflipped ? '\'':' ');
			    printf("path distance err=%0.3f,sd=%0.3f\n",fabs(err),sqrt(var));

			    printf("mincov1 without critical edges = %0.6f, mincov1=%0.6f\n",minC1,mincov1);
			    fflush(stdout);

			    assert(minC1 >= mincov1 - 1e-6);
			  }
			}
			for(Cnode *Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			  Cedge *Tpedge = DATA(BFS,Tnode,pedge);
			  if(!Tpedge->mark && Tpedge->cov < minC2)
			    minC2 = Tpedge->cov;
			}
			if(!PNedge->mark && PNedge->cov < minC2)
			  minC2 = PNedge->cov;
			if(DEBUG && !(minC2 >= mincov2 - 1e-6)){
			  #pragma omp critical
			  {
			    printf("tid=%d:Found two paths from Vnode=N%d%c to Nnode=N%d (Pnode=N%d%c)\nPath1:distance=%8.3f,cov=%5.3f,mark=%d: N%d%c(cov=%0.3f)",
				   tid,nodeorder[Vnode->mapid], Vflipped ? '\'':' ',nodeorder[Nnode->mapid],nodeorder[Pnode->mapid],Pflipped?'\'':' ',
				   distance2, mincov2, Pnode->mark, nodeorder[Nnode->mapid], Nflipped ? '\'' : ' ', PNedge->cov);
			    for(Cnode *Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      assert(Tnode != 0);
			      printf(" <- N%d%c(cov=%5.1f)", nodeorder[Tnode->mapid], DATA(BFS,Tnode,pflip) ? '\'' : ' ', DATA(BFS,Tnode,pedge)->cov);
			    }
			    printf(" <- N%d%c\nPath2:distance=%8.3f,cov=%5.3f,mark=%d: N%d%c(cov=%0.3f)", nodeorder[Vnode->mapid], Vflipped ? '\'':' ',
				   DATA(BFS,Nnode,Rdist), DATA(BFS,Nnode,mincov), DATA(BFS,Nnode,mark),nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'' : ' ', DATA(BFS,Nnode,pedge)->cov);
			    for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      if(DEBUG>=2 && Tnode == 0){
				printf("Tnode = 0!\n");
				fflush(stdout);
				printf("Vnode=N%d, Nnode=N%d%c, Pnode=N%d%c\n", nodeorder[Vnode->mapid],
				       nodeorder[Nnode->mapid],Nflipped ? '\'':' ',
				       nodeorder[Pnode->mapid],Pflipped ? '\'':' ');
				fflush(stdout);
				assert(Tnode != 0);
			      }
			      printf(" <- N%d%c(cov=%5.3f)", nodeorder[Tnode->mapid], DATA(BFS,Tnode,pflip) ? '\'' : ' ', DATA(BFS,Tnode,pedge)->cov);
			    }
			    printf(" <- N%d%c\n", nodeorder[Vnode->mapid],Vflipped ? '\'':' ');
			    printf("path distance err=%0.3f,sd=%0.3f\n",fabs(err),sqrt(var));

			    printf("mincov2 without critical edges = %0.6f, mincov2=%0.6f\n",minC2,mincov2);
			    fflush(stdout);

			    assert(minC2 >= mincov2 - 1e-6);
			  }// assertion failure
			}
			mincov1 = minC1;
			mincov2 = minC2;
			if(mincov1 > HUGE_COV && mincov2 > HUGE_COV)/* all edges are critical : don't delete any edge */
			  mincov1 = mincov2 = 0.0;
		      }
		      if(mincov2 < mincov1){/* delete 2nd path (Nnode to PNedge to Pnode ...) */
			cov1 = mincov2;
			cov2 = -mincov2;
		      } else { /* delete 1st path (Nnode to Nnode->pedge to Nnode->parent ... ) */
			cov1 = -mincov1;
			cov2 = mincov1;
		      }
		    }
		    if(VERB>=2 && verb){
		      printf("cov1=%0.3f(path2 delta),cov2=%0.3f(path1 delta)\n",cov1,cov2);
		      fflush(stdout);
		    }

		    if(BULGE_DEFER){ /* save bulge information for Vnode and handle the actual graph update outside the multithreaded section, since this graph update is serialized by the critical section anyway */
                      if(mycntbulges < maxbulges){
		        #pragma omp atomic capture
                        mycntbulges = cntbulges++;
		      
			if(mycntbulges >= maxbulges){
			  #pragma omp atomic
                          cntbulges--;
		        } else {
			  CBulge *pb = &bulges[mycntbulges];
			  pb->Vnode = Vnode;
			  pb->Xnode = Xnode;
			  pb->Wnode = Wnode;
			  pb->Vflipped = Vflipped;
			  pb->Xflipped = Xflipped;
			  pb->Wflipped = Wflipped;

			  int pathlen = 0;
			  int Tflipped = Nflipped;
			  Cnode *Tnode = Nnode;
			  for(; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent), Tflipped = DATA(BFS,Tnode,pflip) ){
                            pb->PathNode1[pathlen] = Tnode->nodeid;
			    pb->PathFlipped1[pathlen] = Tflipped;
			    pb->PathEdge1[pathlen] = DATA(BFS,Tnode,pedge);
			    pathlen++;
			  }
			  if(DEBUG) assert(pathlen <= BFSdepth+1);
			  pb->PathLen1 = pathlen;
		      
			  pathlen = 0;
			  for(Tflipped = Pflipped, Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent), Tflipped = DATA(BFS,Tnode,pflip) ){
			    pb->PathNode2[pathlen] = Tnode->nodeid;
			    pb->PathFlipped2[pathlen] = Tflipped;
			    pb->PathEdge2[pathlen] = DATA(BFS,Tnode,pedge);
			    pathlen++;
			  }
			  if(DEBUG) assert(pathlen <= BFSdepth);
			  pb->PathLen2 = pathlen;
		      
			  pb->PNedge = PNedge;
			  pb->cov1 = cov1;
			  pb->cov2 = cov2;
		        } /* else : just skip this bulge till next call to graph_bulge() */
                      }
	            } else /* if(!BULGE_DEFER) */ { /* perform graph update */
		    
		      // NOTE : replace with thread-local counter
                      #pragma omp atomic
		      bulgecnt++;

		      #pragma omp critical 
		      {

			/* check that both paths still have all edges with valid >= 1 and that cov1 and cov2 do not exceed any current ->cov values (except for critical edges if BULGE_AVOID) */
			int valid = 1;
			if(BULGE_PAR && nthreads > 1){
			  if(!(DATA(BFS,Nnode,pedge)->valid > 0))
			    valid = 0;
			  if(valid){
			    if(DATA(BFS,Nnode,pedge)->cov + cov1 < 0.0 && !(BULGE_AVOID && DATA(BFS,Nnode,pedge)->mark))
			      cov1 = -DATA(BFS,Nnode,pedge)->cov;
			    for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      if(!(DATA(BFS,Tnode,pedge)->valid > 0)){
				valid = 0;
				if(BULGE_PAR)
				  break;
			      }
			      if(DATA(BFS,Tnode,pedge)->cov + cov1 < 0.0 && !(BULGE_AVOID && DATA(BFS,Tnode,pedge)->mark))
				cov1 = -DATA(BFS,Tnode,pedge)->cov;
			    }
			  }
			  if(!(newPNedge->valid > 0))
			    valid = 0;
			  if(valid) {
			    if(newPNedge->cov + cov2 < 0.0 && !(BULGE_AVOID && newPNedge->mark))
			      cov2 = -newPNedge->cov;
			    for(Cnode *Tnode = newPnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      if(!(DATA(BFS,Tnode,pedge)->valid > 0)){
				valid = 0;
				if(BULGE_PAR)
				  break;
			      }
			      if(DATA(BFS,Tnode,pedge)->cov + cov2 < 0.0 && !(BULGE_AVOID && DATA(BFS,Tnode,pedge)->mark))
				cov2 = -DATA(BFS,Tnode,pedge)->cov;
			    }
			  }
			}

			if(VERB>=2 && verb){
			  printf(" V=N%d%c,W=N%d%c,X=N%d%c:\n",
				 nodeorder[Vnode->nodeid],Vflipped ? '\'' : ' ', nodeorder[Wnode->nodeid], Wflipped ? '\'' : ' ', nodeorder[Xnode->nodeid], Xflipped ? '\'' : ' ');
			  printf("  path1 from N to V:\n");
			  int Tflipped = Nflipped;
			  Cnode *Tnode = Nnode;
			  int pathlen = 0;
			  for(; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent), Tflipped = DATA(BFS,Tnode,pflip), pathlen++ ){
			    Cedge *Tpedge = DATA(BFS,Tnode,pedge);
			    printf("    %d:Tnode=N%d%c, Tedge=(N%d,N%d,%0.3f%c,%d,%d),cov=%0.3f,valid=%d,mark=%d\n",
				   pathlen, nodeorder[Tnode->nodeid], Tflipped ? '\'' : ' ', nodeorder[node[Tpedge->Lmap].mapid], nodeorder[node[Tpedge->Rmap].mapid], 
				   Tpedge->distance, Tpedge->alignid ? ' ':'*', Tpedge->Lflipped,Tpedge->Rflipped,Tpedge->cov,Tpedge->valid,Tpedge->mark);
			  }
			  printf("  path2 from P to V:\n");
			  for(pathlen= 0, Tflipped = Pflipped, Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent), Tflipped = DATA(BFS,Tnode,pflip), pathlen++ ){
			    Cedge *Tpedge = DATA(BFS,Tnode,pedge);
			    printf("    %d:Tnode=N%d%c, Tedge=(N%d,N%d,%0.3f%c,%d,%d),cov=%0.3f,valid=%d,mark=%d\n",
				   pathlen, nodeorder[Tnode->nodeid], Tflipped ? '\'' : ' ', nodeorder[node[Tpedge->Lmap].mapid], nodeorder[node[Tpedge->Rmap].mapid], 
				   Tpedge->distance, Tpedge->alignid ? ' ':'*', Tpedge->Lflipped,Tpedge->Rflipped,Tpedge->cov,Tpedge->valid,Tpedge->mark);
			  }
			  printf("  PNedge=(N%d,N%d,%0.3f%c,%d,%d),cov=%0.3f,valid=%d,mark=%d\n", nodeorder[node[PNedge->Lmap].mapid], nodeorder[node[PNedge->Rmap].mapid],
				 PNedge->distance, PNedge->alignid ? ' ':'*', PNedge->Lflipped,PNedge->Rflipped,PNedge->cov,PNedge->valid,PNedge->mark);
			  fflush(stdout);
			}

			if(DEBUG && !BULGE_PAR) assert(valid);

			if(!BULGE_PAR || valid){/* Update graph */
			  if(cov1 == 0.0 && cov2 != 0.0){	  /* break weakest link(s) (with coverage matching -cov2) on path Nnode PNedge Pnode ... Vnode */
			    if(DEBUG>=2) assert(cov2 < 0.0);
			    if((newPNedge->cov + cov2) <= 0.01 && (!BULGE_PAR || newPNedge->valid > 0)){
			      if(DEBUG>=2)  assert(newPNedge->cov >= -0.01);
			      edge_delete(newPNedge,Nnode);
			      if(DEBUG>=2)  assert(newPNedge->valid < 0);
			      bulgefix++;
			    }
			    for(Cnode *Tnode = newPnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      if(DEBUG>=2) assert(!BULGE_PAR || DATA(BFS,Tnode,pedge)->valid > 0);
			      if((DATA(BFS,Tnode,pedge)->cov + cov2) <= 0.01 /* && (!BULGE_PAR || DATA(BFS,Tnode,pedge)->valid > 0) */){
				if(DEBUG>=2)  assert(DATA(BFS,Tnode,pedge)->cov >= -0.01);
				edge_delete(DATA(BFS,Tnode,pedge),Tnode);
				if(DEBUG>=2)  assert(DATA(BFS,Tnode,pedge)->valid < 0);
				bulgefix++;
			      }
			    }
			  } else if(cov1 != 0.0 && cov2 == 0.0){	/* break weakest link (with coverage matching -cov1) on path Nnode Nnode->pedge Nnode->parent ... Vnode */
			    if(DEBUG>=2) assert(cov1 < 0.0);
			    if(DEBUG>=2) assert(!BULGE_PAR || DATA(BFS,Nnode,pedge)->valid > 0);
			    if((DATA(BFS,Nnode,pedge)->cov + cov1) <= 0.01 /* && (!BULGE_PAR || DATA(BFS,Nnode,pedge)->valid > 0)*/){
			      if(DEBUG>=2) assert(DATA(BFS,Nnode,pedge)->cov >= -0.01);
			      edge_delete(DATA(BFS,Nnode,pedge),Nnode);
			      if(DEBUG>=2) assert(DATA(BFS,Nnode,pedge)->valid < 0);
			      bulgefix++;
			    }
			    for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      if(DEBUG>=2) assert(!BULGE_PAR || DATA(BFS,Tnode,pedge)->valid > 0);
			      if((DATA(BFS,Tnode,pedge)->cov + cov1) <= 0.01 /* && (!BULGE_PAR || DATA(BFS,Tnode,pedge)->valid > 0)*/){
				if(DEBUG>=2) assert(DATA(BFS,Tnode,pedge)->cov >= -0.01);
				edge_delete(DATA(BFS,Tnode,pedge),Tnode);
				if(DEBUG>=2) assert(DATA(BFS,Tnode,pedge)->valid < 0);
				bulgefix++;
			      }
			    }
			  } else if(cov1 != 0.0 && cov2 != 0.0) {
			    /* add cov1 to Nnode->mincov from path Nnode Nnode->pedge Nnode->parent ... Vnode and
			       add cov2 to Nnode->mincov from path Nnode PNedge Pnode ... Vnode and
			       Delete any edges whose cov values drops to 0, unless (BULGE_AVOID and mark==1) : in that case keep cov to at least 1.0. */

			    
			    Cedge *Tpedge = DATA(BFS,Nnode,pedge);
			    Tpedge->cov += cov1;
			    if(BULGE_AVOID && Tpedge->mark && Tpedge->cov < 1.0)
			      Tpedge->cov = 1.0;
			    else if(Tpedge->cov <= 0.01){
			      if(DEBUG>=2) assert(Tpedge->cov >= -0.01);
			      edge_delete(Tpedge,Nnode);
			      if(DEBUG>=2) assert(Tpedge->valid < 0);
			      bulgefix++;
			    }
			    for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      if((Tnode->cov += cov1) <= 0.99)
				Tnode->cov = 1;
			      Tpedge = DATA(BFS,Tnode,pedge);
			      Tpedge->cov += cov1;
			      if(BULGE_AVOID && Tpedge->mark && Tpedge->cov < 1.0)
				Tpedge->cov = 1.0;
			      else if(Tpedge->cov <= 0.01){
				if(DEBUG>=2) assert(Tpedge->cov >= -0.01);
				edge_delete(Tpedge,Tnode);
				if(DEBUG>=2) assert(Tpedge->valid < 0);
				bulgefix++;
			      }
			    }

			    newPNedge->cov += cov2;
			    if(BULGE_AVOID && newPNedge->mark && newPNedge->cov < 1.0)
			      newPNedge->cov = 1.0;
			    else if(newPNedge->cov <= 0.01){
			      if(DEBUG>=2) assert(newPNedge->cov >= -0.01);
			      edge_delete(newPNedge,Nnode);
			      if(DEBUG>=2) assert(newPNedge->valid < 0);
			      bulgefix++;
			    }
			    for(Cnode *Tnode = newPnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      if((Tnode->cov += cov2) <= 0.99)
				Tnode->cov = 1;
			      Tpedge = DATA(BFS,Tnode,pedge);
			      Tpedge->cov += cov2;
			      if(BULGE_AVOID && Tpedge->mark && Tpedge->cov < 1.0)
				Tpedge->cov = 1.0;
			      else if(Tpedge->cov <= 0.01){
				if(DEBUG>=2 && !(Tpedge->cov >= -0.01)){
				  Cedge *VXedge = Tpedge;
				  printf("edge_delete:deleting VXedge=(N%d,N%d,%0.3f%c,%d,%d),valid=%d,cov=%0.3f,Vnode=N%d\n",
					 nodeorder[node[VXedge->Lmap].mapid],nodeorder[node[VXedge->Rmap].mapid],VXedge->distance,VXedge->alignid ? ' ':'*',
					 VXedge->Lflipped,VXedge->Rflipped,VXedge->valid,VXedge->cov,nodeorder[Tnode->mapid]);
				  printf("cov2=%0.3f,cov1=%0.3f\n",cov2,cov1);
				  fflush(stdout);
				  
				  assert(Tpedge->cov >= -0.01);
				}
				edge_delete(Tpedge,Tnode);
				if(DEBUG>=2) assert(Tpedge->valid < 0);
				bulgefix++;
			      }
			    }
			  }
			  if(VERB>=2 && verb && (cov1 || cov2)){
			    if(cov1==0.0 || cov2 == 0.0)
			      printf("  V=%d(flip=%d):(Deleted weakest link(cov1=%0.3f,cov2=%0.3f): bulgecnt=%d,bulgefix=%d,i=%d/%d)\n",
				     Vnode->nodeid,Vflipped,cov1,cov2,bulgecnt,bulgefix,i,numtopnodes);
			    else
			      printf("  V=%d(flip=%d):(Merged weakest path(cov1=%0.3f,cov2=%0.3f): bulgecnt=%d,bulgefix=%d,i=%d/%d)\n",
				     Vnode->nodeid,Vflipped,cov1,cov2,bulgecnt,bulgefix,i,numtopnodes);
			  } else if(VERB>=2 && verb){
			    printf("  V=%d(flip=%d):(Nothing changed (cov1=%0.3f,cov2=%0.3f): bulgecnt=%d,bulgefix=%d,i=%d/%d)\n",
				   Vnode->nodeid,Vflipped,cov1,cov2,bulgecnt,bulgefix,i,numtopnodes);
			    printf("mincov2=%0.1f,mincov1=%0.1f\n",mincov2,DATA(BFS,Nnode,mincov));
			    printf("tid=%d:Found two paths from Vnode=N%d%c to Nnode=N%d (Pnode=N%d%c)\nPath1:distance=%8.3f,cov=%5.3f,mark=%d: N%d%c(cov=%0.3f)",
				   tid,nodeorder[Vnode->mapid], Vflipped ? '\'':' ',nodeorder[Nnode->mapid],nodeorder[Pnode->mapid],Pflipped?'\'':' ',
				   distance2, mincov2, Pnode->mark, nodeorder[Nnode->mapid], Nflipped ? '\'' : ' ', PNedge->cov);
			    for(Cnode *Tnode = Pnode; Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      assert(Tnode != 0);
			      printf(" <- N%d%c(cov=%5.1f)", nodeorder[Tnode->mapid], DATA(BFS,Tnode,pflip) ? '\'' : ' ', DATA(BFS,Tnode,pedge)->cov);
			    }
			    printf(" <- N%d%c\nPath2:distance=%8.3f,cov=%5.3f,mark=%d: N%d%c(cov=%0.3f)", nodeorder[Vnode->mapid], Vflipped ? '\'':' ',
				   DATA(BFS,Nnode,Rdist), DATA(BFS,Nnode,mincov), DATA(BFS,Nnode,mark),nodeorder[Nnode->mapid], DATA(BFS,Nnode,pflip) ? '\'' : ' ', DATA(BFS,Nnode,pedge)->cov);
			    for(Cnode *Tnode = DATA(BFS,Nnode,parent); Tnode != Vnode; Tnode = DATA(BFS,Tnode,parent)){
			      if(DEBUG>=2 && Tnode == 0){
				printf("Tnode = 0!\n");
				fflush(stdout);
				printf("Vnode=N%d, Nnode=N%d%c, Pnode=N%d%c\n", nodeorder[Vnode->mapid],
				       nodeorder[Nnode->mapid],Nflipped ? '\'':' ',
				       nodeorder[Pnode->mapid],Pflipped ? '\'':' ');
				fflush(stdout);
				assert(Tnode != 0);
			      }
			      printf(" <- N%d%c(cov=%5.3f)", nodeorder[Tnode->mapid], DATA(BFS,Tnode,pflip) ? '\'' : ' ', DATA(BFS,Tnode,pedge)->cov);
			    }
			    printf(" <- N%d%c\n", nodeorder[Vnode->mapid],Vflipped ? '\'':' ');
			    printf("path distance err=%0.3f,sd=%0.3f\n",fabs(err),sqrt(var));
			  }
			  if(DEBUG>=2) assert(found);
			} // if (!BULGE_PAR || valid)
		      } // critical section 

		    } // (!BULGE_DEFER)

		    break;/* terminate BFS */
		  } // if(Nnode->parent)

		  if(DEBUG>=2) assert(0 <= Nnode->nodeid && Nnode->nodeid < numnodes && Nnode-node == Nnode->nodeid);
		  if(DEBUG>=2) assert(0 <= Pnode->nodeid && Pnode->nodeid < numnodes && Pnode-node == Pnode->nodeid);

		  DATA(BFS,Nnode,parent) = Pnode; 
		  DATA(BFS,Nnode,pedge) = PNedge;
		  DATA(BFS,Nnode,Rdist) = distance2;//Pnode->Rdist + PNedge->distance;
		  if(BULGE_PAR)
		    DATA(BFS,Nnode,mincov) = min(DATA(BFS,Pnode,mincov),PNedge->cov);
		  else {
		    if(DEBUG>=2) assert(mincov2 == min(DATA(BFS,Pnode,mincov),PNedge->cov));
		    DATA(BFS,Nnode,mincov) = mincov2;//min(Pnode->mincov,PNedge->cov);
		  }
		  DATA(BFS,Nnode,mark) = DATA(BFS,Pnode,mark);
		  DATA(BFS,Nnode,pflip) = Nflipped;
		  if(VERB>=2 && verb>=2){
                    #pragma omp critical
		    {
		      printf("tid=%d:queue.push:priority=%0.1f,Node=N%d%c,parent=N%d,mincov=%0.1f,Rdist=%0.3f,wt=%0.3f,mark=%d\n",
			     tid,top.priority+1.0,nodeorder[Nnode->mapid],DATA(BFS,Nnode,pflip) ? '\'':' ',nodeorder[Pnode->mapid],
			     DATA(BFS,Nnode,mincov),DATA(BFS,Nnode,Rdist),0.0,DATA(BFS,Nnode,mark));
		      fflush(stdout);
		    }
		  }
		  BFSqueue.push(Cpath(Nnode,Nflipped,top.priority+1.0));
		}
		if(pn < (Pflipped ? Pnode->PosEdge : Pnode->edgecnt))/* early termination */
		  break;/* terminate BFS */
		if(DEBUG>=2 && !(found == 0)){
		  printf("Vnode=N%d%c,Wnode=N%d%c,Xnode=N%d%c:Pnode=N%d%c,pn=%d,Pnode->PosEdge=%d,Pnode->edgecnt=%d,found=%d\n",
			 nodeorder[Vnode->mapid],Vflipped ? '\'':' ',nodeorder[Wnode->mapid],Wflipped ? '\'':' ',nodeorder[Xnode->mapid],Xflipped ? '\'':' ',
			 nodeorder[Pnode->mapid],Pflipped ?'\'':' ',pn,Pnode->PosEdge,Pnode->edgecnt,found);
		  printf("failed to break out of pn loop\n");
		  fflush(stdout);
		  assert(found == 0);
		}
	      }
	      /* clean up : reset parent/mark fields in FinishedNodes or BFSqueue */
	      while(!FinishedNodes.empty()){
		Cnode *Pnode = FinishedNodes.top(); FinishedNodes.pop();
		if(DEBUG>=2)assert(DATA(BFS,Pnode,parent));
		DATA(BFS,Pnode,parent) = 0;
		DATA(BFS,Pnode,mark) = 0;
	      }
	      while(BFSqueue.size() > 0){
		Cpath top = BFSqueue.top(); BFSqueue.pop();
		Cnode *Pnode = top.node;
		if(DEBUG>=2)assert(DATA(BFS,Pnode,parent));
		if(VERB>=2 && verb>=2)
		  printf("queue.pop():priority=%0.1f,Node=N%d,parent=N%d,mincov=%0.1f,Rdist=%0.3f,mark=%d(cleanup)\n",
			 top.priority,nodeorder[Pnode->mapid],nodeorder[DATA(BFS,Pnode,parent)->mapid],
			 DATA(BFS,Pnode,mincov),DATA(BFS,Pnode,Rdist),DATA(BFS,Pnode,mark));
		DATA(BFS,Pnode,parent) = 0;
		DATA(BFS,Pnode,mark) = 0;
	      }
	      if(DEBUG>=4){/* all nodes except Vnode should have parent field == 0 */
		for(int t=0; t < numtopnodes; t++){
		  Cnode *Tnode = topnodes[t];
		  if(Tnode != Vnode && DATA(BFS,Tnode,parent)){
		    printf("topnode[%d]=N%d,Vnode=N%d,Tnode->parent=N%d(should be 0), mark=%d\n",
			   t,nodeorder[Tnode->mapid],nodeorder[Vnode->mapid],nodeorder[DATA(BFS,Tnode,parent)->mapid],DATA(BFS,Tnode,mark));
		    fflush(stdout);
		    assert(DATA(BFS,Tnode,parent) == 0);
		  }
		  assert(DATA(BFS,Tnode,mark) == 0);
		  assert(DATA(BFS,Tnode,SPmark) == 0);
		}
	      }
	      if(found && VWedge->valid < 0)/* outer loop edge deleted : break out of loop */
		break;
	    }
	  }
	  DATA(BFS,Vnode,parent) = 0;
	  if(DEBUG>=3){/* all nodes should have parent field == 0, mark == SPmark = 0 */
	    for(int t=0; t < numtopnodes; t++){
	      Cnode *Tnode = topnodes[t];
	      assert(DATA(BFS,Tnode,parent) == 0);
	      assert(DATA(BFS,Tnode,mark) == 0);
	      assert(DATA(BFS,Tnode,SPmark) == 0);
	      assert(DATA(SP,Tnode,SPmark) == 0);
	    }
	  }
	} // Vflipped = 0 .. 1
      } // omp for : Vnode = topnodes[i = 0 .. numtopnodes-1]
      if(DEBUG>=2){/* all nodes should have parent field == 0, mark == SPmark = 0 */
	for(int t=0; t < numtopnodes; t++){
	  Cnode *Tnode = topnodes[t];
	  assert(0 <= Tnode->nodeid && Tnode->nodeid < numnodes);
	  assert(DATA(BFS,Tnode,parent) == 0);
	  assert(DATA(BFS,Tnode,mark) == 0);
	  assert(DATA(BFS,Tnode,SPmark) == 0);
	  assert(DATA(SP,Tnode,SPmark) == 0);
	}
      }
      /* free per-thread memory */
      if(BFS) delete [] BFS;      
      if(SP) delete [] SP;
    } // omp parallel

    if(CYCLE_DEFER){/* perform graph updates to break cyles in a deterministic order */
      if(cntCycles >= maxCycles - numthreads){
	printf("WARNING:increase maxCycles=%d (numtopnodes=%d,cntCycles=%d)\n",maxCycles,numtopnodes,cntCycles);
	fflush(stdout);
      }

      qsort(Cycles,cntCycles,sizeof(CCycle),(intcmp*)CCycleNodeInc);/* make order deterministic */

      for(int i = 0; i < cntCycles; i++){
	CCycle *pc = &Cycles[i];

	/* check that all cycle edges still have valid >= 1 */
	if(pc->DelEdge->valid <= 0)
	  continue;
	if(DEBUG>=2) assert(pc->PathLen > 0);
	int valid = 1;
	for(int pathlen = 0; pathlen < pc->PathLen; pathlen++){
	  Cedge *Tpedge = pc->PathEdge[pathlen];
	  if(!(Tpedge->valid > 0)){
	    valid = 0;
	    break;
	  }
	}
	if(!valid)
	  continue;
	cyclefix++;
	edge_delete(pc->DelEdge,pc->DelNode);
      }

      /* free up memory */
      delete [] Cycles;
      delete [] CEdgeMem;
    }

    if(BULGE_DEFER){/* perform graph updates to break/merge graph bulges */
      if(VERB>=2){/* locate highest coverage edge */
	double maxCov = 0.0;
	Cedge *maxEdge = 0;
	for(int t=0; t < numtopnodes; t++){
	  Cnode *Tnode = topnodes[t];
	  for(int u = 0; u < Tnode->edgecnt; u++){
	    Cedge *TUedge = &edge[Tnode->edgeid[u]];
	    if(DEBUG) assert(TUedge->cov >= 0.0);
	    if(TUedge->cov > maxCov){
	      maxCov = TUedge->cov;
	      maxEdge = TUedge;
	    }
	  }	
	}
	printf("graph_bulge:BULGE_DEFER: Largest coverage = %0.2f : (N%d,N%d,%0.3f,%d,%d)\n",
	       maxCov, nodeorder[node[maxEdge->Lmap].mapid],nodeorder[node[maxEdge->Rmap].mapid],maxEdge->distance,
	       maxEdge->Lflipped,maxEdge->Rflipped);
	fflush(stdout);
      }

      if(VERB>= 2 && verb){
	printf("bulges=%d/%d:\n",cntbulges,maxbulges);
	fflush(stdout);
      }
      qsort(bulges,cntbulges,sizeof(CBulge),(intcmp*)CBulgeNodeInc);/* make order deterministic */

      for(int i = 0; i < cntbulges; i++){
	CBulge *pb = &bulges[i];

	Cnode *Vnode = pb->Vnode;
	int Vflipped = pb->Vflipped;
	//	Cnode *Xnode = pb->Xnode;
	//	Cnode *Wnode = pb->Wnode;
	if(DEBUG>=2) assert(pb->PathLen1 > 0);

	bulgecnt++;

	double cov1 = pb->cov1;
	double cov2 = pb->cov2;

	/* check that both paths still have all edges with valid >= 1 and that cov1 and cov2 do not exceed any current ->cov values (except for critical edges if BULGE_AVOID) */
	int valid = 1;
	Cedge *PNedge = pb->PNedge;

	for(int pathlen = 0; pathlen < pb->PathLen1; pathlen++){
	  Cedge *Tpedge = pb->PathEdge1[pathlen];
	  if(!(Tpedge->valid > 0)){
	    valid = 0;
	    break;
	  }
	  if(Tpedge->cov + cov1 < 0.0 && !(BULGE_AVOID && Tpedge->mark))
	    cov1 = -Tpedge->cov;
	}

	if(!(PNedge->valid > 0))
	  valid = 0;
	if(valid){
	  if(PNedge->cov + cov2 < 0.0 && !(BULGE_AVOID && PNedge->mark))
	    cov2 = -PNedge->cov;
	  for(int pathlen = 0; pathlen < pb->PathLen2; pathlen++){
	    Cedge *Tpedge = pb->PathEdge2[pathlen];
	    if(!(Tpedge->valid > 0)){
	      valid = 0;
	      break;
	    }
	    if(Tpedge->cov + cov2 < 0.0 && !(BULGE_AVOID && Tpedge->mark))
	      cov2 = -Tpedge->cov;
	  }
	}

	if((VERB>=2 && verb) || cov1 > HUGE_COV || cov2 > HUGE_COV){
	  printf(" V=N%d%c,W=N%d%c,X=N%d%c:\n",
		 nodeorder[Vnode->mapid],Vflipped ? '\'' : ' ', nodeorder[pb->Wnode->mapid], pb->Wflipped ? '\'' : ' ', nodeorder[pb->Xnode->mapid], pb->Xflipped ? '\'' : ' ');
	  printf("  path1 from N to V (pathlen=%d):\n",pb->PathLen1);
	  for(int pathlen = 0; pathlen < pb->PathLen1; pathlen++){
	    Cnode *Tnode = &node[pb->PathNode1[pathlen]];
	    Cedge *Tpedge = pb->PathEdge1[pathlen];
	    printf("    %d:Tnode=N%d%c, Tedge=(N%d,N%d,%0.3f%c,%d,%d),cov=%0.3f,valid=%d,mark=%d\n",
		   pathlen, nodeorder[Tnode->mapid], pb->PathFlipped1[pathlen] ? '\'' : ' ', nodeorder[node[Tpedge->Lmap].mapid], nodeorder[node[Tpedge->Rmap].mapid], Tpedge->distance, Tpedge->alignid ? ' ':'*',
		   Tpedge->Lflipped,Tpedge->Rflipped,Tpedge->cov,Tpedge->valid,Tpedge->mark);
	  }
	  printf("  path2 from P to V (pathlen=%d):\n",pb->PathLen2);
	  for(int pathlen = 0; pathlen < pb->PathLen2; pathlen++){
	    Cnode *Tnode = &node[pb->PathNode2[pathlen]];
	    Cedge *Tpedge = pb->PathEdge2[pathlen];
	    printf("    %d:Tnode=N%d%c, Tedge=(N%d,N%d,%0.3f%c,%d,%d),cov=%0.3f,valid=%d,mark=%d\n",
		   pathlen, nodeorder[Tnode->mapid], pb->PathFlipped2[pathlen] ? '\'' : ' ', nodeorder[node[Tpedge->Lmap].mapid], nodeorder[node[Tpedge->Rmap].mapid], Tpedge->distance, Tpedge->alignid ? ' ':'*',
		   Tpedge->Lflipped,Tpedge->Rflipped,Tpedge->cov,Tpedge->valid,Tpedge->mark);
	  }
	  printf("  PNedge=(N%d,N%d,%0.3f%c,%d,%d),cov=%0.3f,valid=%d,mark=%d\n", nodeorder[node[PNedge->Lmap].mapid], nodeorder[node[PNedge->Rmap].mapid], PNedge->distance, PNedge->alignid ? ' ':'*',
		 PNedge->Lflipped,PNedge->Rflipped,PNedge->cov,PNedge->valid,PNedge->mark);
	  if(cov1 > HUGE_COV || cov2 > HUGE_COV)
	    printf("cov1=%0.3f,cov2=%0.3f\n",cov1,cov2);
	  fflush(stdout);
	}

	if(DEBUG) assert(cov1 <= HUGE_COV);
	if(DEBUG) assert(cov2 <= HUGE_COV);

	if(!valid){
	  if(VERB>=2 && verb){
	    printf("  V=%d(flip=%d):(Nothing changed due to missing edges (cov1=%0.3f,cov2=%0.3f): bulgecnt=%d,bulgefix=%d,i=%d/%d:\n",
		   Vnode->nodeid,Vflipped,cov1,cov2,bulgecnt,bulgefix,i,numtopnodes);
	    fflush(stdout);
	  }
	  continue;
	}

	/* Update graph */

	if(cov1 == 0.0 && cov2 != 0.0){	  /* break weakest link(s) (with coverage matching -cov2) on path Nnode PNedge Pnode ... Vnode */
	  if(DEBUG>=2) assert(cov2 < 0.0);
	  if(PNedge->cov + cov2 <= 0.01){
	    if(DEBUG>=2)  assert(PNedge->cov >= -0.01);
	    edge_delete(PNedge,&node[pb->PathNode1[0]]);
	    if(DEBUG>=2)  assert(PNedge->valid < 0);
	    bulgefix++;
	  }
	  for(int pathlen = 0; pathlen < pb->PathLen2; pathlen++){
	    Cnode *Tnode = &node[pb->PathNode2[pathlen]];
	    Cedge *Tpedge = pb->PathEdge2[pathlen];
	    if(DEBUG>=2) assert(Tpedge->valid > 0);
	    if((Tpedge->cov + cov2) <= 0.01){
	      if(DEBUG>=2)  assert(Tpedge->cov >= -0.01);
	      edge_delete(Tpedge,Tnode);
	      if(DEBUG>=2)  assert(Tpedge->valid < 0);
	      bulgefix++;
	    }
	  }
	} else if(cov1 != 0.0 && cov2 == 0.0){	/* break weakest link (with coverage matching -cov1) on path Nnode Nnode->pedge Nnode->parent ... Vnode */
	  if(DEBUG>=2) assert(cov1 < 0.0);
	  for(int pathlen = 0; pathlen < pb->PathLen1; pathlen++){
	    Cnode *Tnode = &node[pb->PathNode1[pathlen]];
	    Cedge *Tpedge = pb->PathEdge1[pathlen];
	    if(DEBUG>=2) assert(Tpedge->valid > 0);
	    if((Tpedge->cov + cov1) <= 0.01){
	      if(DEBUG>=2) assert(Tpedge->cov >= -0.01);
	      edge_delete(Tpedge,Tnode);
	      if(DEBUG>=2) assert(Tpedge->valid < 0);
	      bulgefix++;
	    }
	  }
	} else if(cov1 != 0.0 && cov2 != 0.0) {
	  /* add cov1 to Nnode->mincov from path Nnode Nnode->pedge Nnode->parent ... Vnode and
	     add cov2 to Nnode->mincov from path Nnode PNedge Pnode ... Vnode and
	     Delete any edges whose cov values drops to 0, unless (BULGE_AVOID and mark==1) : in that case keep cov = 1.0 */
	      
	  for(int pathlen = 0; pathlen < pb->PathLen1; pathlen++){
	    Cnode *Tnode = &node[pb->PathNode1[pathlen]];
	    Cedge *Tpedge = pb->PathEdge1[pathlen];
	    double origcov = Tpedge->cov;
	    if(pathlen > 0 && (Tnode->cov += cov1) <= 0.99)
	      Tnode->cov = 1;
	    Tpedge->cov += cov1;
	    if(BULGE_AVOID && Tpedge->mark && Tpedge->cov < 1.0)
	      Tpedge->cov = 1.0;
	    else if(Tpedge->cov <= 0.01){
	      if(DEBUG>=2) assert(Tpedge->cov >= -0.01);
	      edge_delete(Tpedge,Tnode);
	      if(DEBUG>=2) assert(Tpedge->valid < 0);
	      bulgefix++;
	    }
	    if(DEBUG>=2 && !(Tpedge->cov <= HUGE_COV)){
	      #pragma omp critical
	      {
		printf("cov1=%0.1f,cov2=%0.1f,Tpedge->cov=%0.1f -> %0.1f\n",cov1,cov2,origcov,Tpedge->cov);
		fflush(stdout);
		assert(Tpedge->cov <= HUGE_COV);
	      }
	    }
	  }
	      
	  PNedge->cov += cov2;
	  if(BULGE_AVOID && PNedge->mark && PNedge->cov < 1.0)
	    PNedge->cov = 1.0;
	  else if(PNedge->cov <= 0.01){
	    if(DEBUG>=2) assert(PNedge->cov >= -0.01);
	    edge_delete(PNedge,&node[pb->PathNode1[0]]);
	    if(DEBUG>=2) assert(PNedge->valid < 0);
	    bulgefix++;
	  }
	  if(DEBUG>=2) assert(PNedge->cov <= HUGE_COV);
	  for(int pathlen = 0; pathlen < pb->PathLen2; pathlen++){
	    Cnode *Tnode = &node[pb->PathNode2[pathlen]];
	    Cedge *Tpedge = pb->PathEdge2[pathlen];
	    if((Tnode->cov += cov2) <= 0.99)
	      Tnode->cov = 1;
	    Tpedge->cov += cov2;
	    if(BULGE_AVOID && Tpedge->mark && Tpedge->cov < 1.0)
	      Tpedge->cov = 1.0;
	    else if(Tpedge->cov <= 0.01){
	      if(DEBUG>=2) assert(Tpedge->cov >= -0.01);
	      edge_delete(Tpedge,Tnode);
	      if(DEBUG>=2) assert(Tpedge->valid < 0);
	      bulgefix++;
	    }
	    if(DEBUG>=2) assert(Tpedge->cov <= HUGE_COV);
	  }
	} // if (cov1 != 0.0 && cov2 != 0.0)

	if(VERB>=2 && verb && (cov1 || cov2)){
	  if(cov1==0.0 || cov2 == 0.0)
	    printf("  V=%d(flip=%d):(Deleted weakest link(cov1=%0.3f,cov2=%0.3f): bulgecnt=%d,bulgefix=%d,i=%d/%d)\n",
		   Vnode->nodeid,Vflipped,cov1,cov2,bulgecnt,bulgefix,i,numtopnodes);
	  else
	    printf("  V=%d(flip=%d):(Merged weakest path(cov1=%0.3f,cov2=%0.3f): bulgecnt=%d,bulgefix=%d,i=%d/%d)\n",
		   Vnode->nodeid,Vflipped,cov1,cov2,bulgecnt,bulgefix,i,numtopnodes);
	} else if(VERB>=2 && verb)
	  printf("  V=%d(flip=%d):(Nothing changed (cov1=%0.3f,cov2=%0.3f): bulgecnt=%d,bulgefix=%d,i=%d/%d)\n",
		 Vnode->nodeid,Vflipped,cov1,cov2,bulgecnt,bulgefix,i,numtopnodes);
      } // for i = 0 .. cntbulges -1 
      
      /* delete graph bulge memory */
      delete [] bulges;
      delete [] NodeMem;
      delete [] EdgeMem;
    }

    if(DEBUG>=2){/* all nodes should have parent == 0, mark == SPmark = 0 */
      for(int t=0; t < numtopnodes; t++){
	Cnode *Tnode = topnodes[t];
	assert(Tnode->parent == 0);
	assert(Tnode->mark == 0);
	assert(Tnode->SPmark == 0);
      }
    }

    /* remove edges with valid <= -1 from edge lists of all nodes */
    (void)GCedges(finalcnt,edge);
    if(VERB){
      if(bulgecnt > 0 || verb)
	printf("graph_bulge(%d,%0.3f): %d complex bulges found(%d edges deleted)(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n", BFSdepth, MaxBulge, bulgecnt,bulgefix, finalcnt,numtopnodes,wtime(),mtime());
      if(cyclecnt > 0 || verb)
	printf("graph_bulge(%d,%0.3f): %d cycles found(%d edges deleted)(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n", BFSdepth, MaxBulge, cyclecnt, cyclefix, finalcnt,numtopnodes,wtime(),mtime());
      if(VERB>=2)
	printf("BFScnt=%d\n",BFScnt2);
      fflush(stdout);
    }
  }

  if(VERB>=2){/* locate highest coverage edge */
    double maxCov = 0.0;
    Cedge *maxEdge = 0;
    for(int t=0; t < numtopnodes; t++){
      Cnode *Tnode = topnodes[t];
      for(int u = 0; u < Tnode->edgecnt; u++){
	Cedge *TUedge = &edge[Tnode->edgeid[u]];
	if(DEBUG) assert(TUedge->cov >= 0.0);
	if(TUedge->cov > maxCov){
	  maxCov = TUedge->cov;
	  maxEdge = TUedge;
	}
      }	
    }
    printf("graph_bulge(End): Largest coverage = %0.2f : (N%d,N%d,%0.3f,%d,%d)\n",
	   maxCov, nodeorder[node[maxEdge->Lmap].mapid],nodeorder[node[maxEdge->Rmap].mapid],maxEdge->distance,
	   maxEdge->Lflipped,maxEdge->Rflipped);
    fflush(stdout);
  }

  return bulgefix + cyclefix;
}

/** remove edges with coverage < MinCoverage.
 Delete (chimeric) nodes with low node coverage AND total coverage of both "left" OR "right" edges above MinCoverage (see ChimRatio) (Done once) 
 Delete edges with coverage well below that of the top two edges for each node by a ratio of EdgeCoverageRatio 
 If sidechain != 0 : Delete side nodes connected by a single edge to the middle of a long chain (This is harmful if side chains are not recovered later) 
 returns number of edges + nodes removed */
int graph_prune(int sidechain,int verb)
{
  if(VERB>=2){/* locate highest coverage edge */
    double maxCov = 0.0;
    Cedge *maxEdge = 0;
    for(int t=0; t < numtopnodes; t++){
      Cnode *Tnode = topnodes[t];
      for(int u = 0; u < Tnode->edgecnt; u++){
	Cedge *TUedge = &edge[Tnode->edgeid[u]];
	if(DEBUG) assert(TUedge->cov >= 0.0);
	if(TUedge->cov > maxCov){
	  maxCov = TUedge->cov;
	  maxEdge = TUedge;
	}
      }	
    }
    printf("graph_prune(Start): Largest coverage = %0.2f : (N%d,N%d,%0.3f,%d,%d)\n",
	   maxCov, nodeorder[node[maxEdge->Lmap].mapid],nodeorder[node[maxEdge->Rmap].mapid],maxEdge->distance,
	   maxEdge->Lflipped,maxEdge->Rflipped);
    fflush(stdout);
  }

  numthreads = 1;
  #ifdef _OPENMP
  numthreads = omp_get_max_threads();
  if(DEBUG) assert(numthreads >= 1);
  if(MaxThreads > 0)
    numthreads = min(numthreads,MaxThreads);
  #endif

  if(!PosEdgeInit)
    PosEdge();

  int mincnt = 0;/* edges below MinCoverage */
  int weakcnt = 0;/* edges below EdgeCoverageRatio * 2nd highest edge coverage */
  int chimcnt=0;
  int chimFPcnt=0;

  for(int i=0; i < numtopnodes; i++){
    Cnode *Vnode = topnodes[i];
    for(int vw = 0; vw < Vnode->edgecnt; vw++){
      Cedge *VWedge = &edge[Vnode->edgeid[vw]];
      if(VWedge->valid < 0)
	continue;
      if(floor(VWedge->cov + 0.5) < MinCoverage){
	if(VERB>=2 && verb){
	  printf("graph_prune(%d,%d):VWedge=(N%d,N%d,%0.3f,%d,%d):cov=%0.2f,MinCoverage=%d(Deleted)\n",
		 sidechain,verb,nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,
		 VWedge->Lflipped,VWedge->Rflipped,VWedge->cov,MinCoverage);
	  fflush(stdout); 
	}
	mincnt++;
	edge_delete(VWedge,Vnode);
      }
    }

    if(EdgeCoverageRatio > 0.0 && Vnode->edgecnt >= 2){
      if(Vnode->PosEdge >= 2){
	/* locate dominant backward edge */
	double cov1 = 0.0;
	for(int vw = 0; vw < Vnode->PosEdge; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->cov > cov1)
	    cov1 = VWedge->cov;
	}
	cov1 *= EdgeCoverageRatio;
	for(int vw = 0; vw < Vnode->PosEdge; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->valid < 0)
	    continue;
	  //	    Cnode *Wnode = &node[(VWedge->Lmap == Vnode->mapid) ? VWedge->Rmap : VWedge->Lmap];
	  if(VWedge->cov < cov1/* && (Wnode->edgecnt >= 2 || max(Wnode->Lcntmax,Wnode->Rcntmax) * EdgeCoverageRatio >= 1.0)*/){
	    if(VERB>=2 && verb){
	      printf("Vnode=N%d,vw=%d,VWedge=(N%d,N%d,%0.3f,%d,%d,cov=%0.1f),cov1=%0.1f (deleting backward edge VWedge)\n",
		     nodeorder[Vnode->mapid],vw,nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,
		     VWedge->Lflipped,VWedge->Rflipped,VWedge->cov,cov1);
	      fflush(stdout);
	    }
	    weakcnt++;
	    edge_delete(VWedge,Vnode);
	  } else if(VERB>=2 && verb){
	    printf("Vnode=N%d,vw=%d,VWedge=(N%d,N%d,%0.3f,%d,%d,cov=%0.1f),cov1=%0.1f (backward edge VWedge not deleted)\n",
		   nodeorder[Vnode->mapid],vw,nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,
		   VWedge->Lflipped,VWedge->Rflipped,VWedge->cov,cov1);
	    fflush(stdout);
	  }
	}
      }
      if(Vnode->edgecnt >= Vnode->PosEdge + 2){/* locate dominant forward edge */
	double cov2 = 0.0;
	for(int vw = Vnode->PosEdge; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->cov > cov2)
	    cov2 = VWedge->cov;
	}
	cov2 *= EdgeCoverageRatio;
	for(int vw = Vnode->PosEdge; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->valid < 0)
	    continue;
	  //	    Cnode *Wnode = &node[(VWedge->Lmap == Vnode->mapid) ? VWedge->Rmap : VWedge->Lmap];
	  if(VWedge->cov < cov2/* && (Wnode->edgecnt >= 2 || max(Wnode->Lcntmax,Wnode->Rcntmax) * EdgeCoverageRatio >= 1.0)*/){
	    if(VERB>=2 && verb){
	      printf("Vnode=N%d,vw=%d,VWedge=(N%d,N%d,%0.3f,%d,%d,cov=%0.1f),cov2=%0.1f (deleting forward edge VWedge)\n",
		     nodeorder[Vnode->mapid],vw,nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,
		     VWedge->Lflipped,VWedge->Rflipped,VWedge->cov,cov2);
	      fflush(stdout);
	    }
	    weakcnt++;
	    edge_delete(VWedge,Vnode);
	  } else 	      if(VERB>=2 && verb){
	    printf("Vnode=N%d,vw=%d,VWedge=(N%d,N%d,%0.3f,%d,%d,cov=%0.1f),cov2=%0.1f (forward edge VWedge not deleted)\n",
		   nodeorder[Vnode->mapid],vw,nodeorder[node[VWedge->Lmap].mapid],nodeorder[node[VWedge->Rmap].mapid],VWedge->distance,
		   VWedge->Lflipped,VWedge->Rflipped,VWedge->cov,cov2);
	    fflush(stdout);
	  }
	}
      }
    }
  }


  unsigned int finalcnt;
  long long Ecnt = GCedges(finalcnt,edge);
  int Ncnt = GCnodes();

  if(VERB && (Ecnt > 0 || Ncnt > 0) ){
    printf("graph_prune(%d,%d): Pruned %d mincov edges(MinCoverage=%d), %d weak edges and %d chimeric nodes (FP=%d,ChimRatio=%0.3f) and %d isolated nodes(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
	   sidechain, verb, mincnt, MinCoverage, weakcnt,chimcnt,chimFPcnt,ChimRatio,Ncnt-chimcnt,finalcnt,numtopnodes,wtime(),mtime());
    fflush(stdout);
  }

  /* look for two nodes with only a single edge between them : if the edge is not a Super-edge, delete the edge (and both nodes) */
  for(int i=0; i < numtopnodes; i++){
    Cnode *Vnode = topnodes[i];
    if(Vnode->edgecnt != 1)
      continue;
    Cedge *VWedge = &edge[Vnode->edgeid[0]];
    Cnode *Wnode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
    if(Wnode->edgecnt == 1 && VWedge->alignid){
      if(VERB && (VERB >=PVERB+1 || verb)){
	printf("deleting N%d and N%d connected by a single edge\n",
	       nodeorder[Vnode->mapid],nodeorder[Wnode->mapid]);
	fflush(stdout);
      }
      edge_delete(VWedge,Vnode);
    }
  }

  Ecnt = GCedges(finalcnt,edge);
  Ncnt = GCnodes();

  if(VERB && (Ecnt > 0 || Ncnt > 0) ){
    printf("graph_prune(%d,%d): Pruned %lld simple edges connecting doubled nodes along with %d nodes(edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
	   sidechain, verb, Ecnt/2, Ncnt, finalcnt,numtopnodes,wtime(),mtime());
    fflush(stdout);
  }

  /* look for edges connecting a single edge node to a longer chain of connected nodes : move the single node
     to a side chain to allow further analysis on the main graph */
  if(sidechain){
    if(!PosEdgeInit)
      PosEdge();

    int nthreads = min(numthreads, max(1,numtopnodes/8));

    int test = 0;
    //    for(test = 8; --test >= 0; ){

    if(VERB>=2){
      printf("Scanning topnodes for single edge nodes: nthreads=%d,numtopnodes=%d:mtime= %0.6f, wtime= %0.6f(test=%d)\n",nthreads,numtopnodes,mtime(),wtime(),test);
      fflush(stdout);
    }

    int *singlenodes = new int[numtopnodes];
    int numsinglenodes = 0;

    for(int i=0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      if(Vnode->edgecnt != 1)
	continue;
      
      Cedge *VWedge = &edge[Vnode->edgeid[0]];
      if(VWedge->distance > MaxSidebranch)
	continue;
      
      singlenodes[numsinglenodes++] = i;
    }
    
    int block = 16;
    /*    while(block < 16 && numsinglenodes > numthreads * block * 8)
	  block *= 2;*/
    //    nthreads = 1;

    if(VERB>=2){
      printf("Starting multithreaded section for sidechain analysis: nthreads=%d,block=%d,numtopnodes=%d,numsinglenodes=%d:mtime= %0.6f, wtime= %0.6f(test=%d)\n",nthreads,block,numtopnodes,numsinglenodes,mtime(),wtime(),test);
      fflush(stdout);
    }

    CSPdata *GSPT[nthreads];
    Cnodeflip *heapT[nthreads];

    //    #pragma omp parallel num_threads(nthreads)
    {
      int tid = 0;
      #ifdef _OPENMP
      tid = omp_get_thread_num ();
      #endif


      CSPdata *GSP = new CSPdata[numnodes];
      for(int i=0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	DATA(GSP,Vnode,SPmark) = 0;
        DATA(GSP,Vnode,hindex) = Vnode->edgecnt;
        DATA(GSP,Vnode,Enodeid)= Vnode->Enodeid;
      }
      GSPT[tid] = GSP;
      heapT[tid] = new Cnodeflip[2*numnodes+1];    
    }

    if(VERB>=2){
      printf("Initialized multithreaded section for sidechain analysis: nthreads=%d,block=%d,numsinglenodes=%d,numtopnodes=%d,numnodes=%d:mtime= %0.6f, wtime= %0.6f\n",
	     nthreads,block,numsinglenodes,numtopnodes,numnodes,mtime(),wtime());
      fflush(stdout);
    }

    // NOTE : for some strange reason multithreading slows down and spends most of the time in some unknown system call (does not show up in perf trace)

    if(VERB>=2 && verb)
      cyclewarning = 0;

    cyclecnt = 0;// NOTE : global cyclecnt

    // HERE HERE    #pragma omp parallel num_threads(nthreads)
    {
      int tid = 0;
      #ifdef _OPENMP
      tid = omp_get_thread_num ();
      #endif


      CSPdata *GSP = GSPT[tid];

      int pathcnt = 0;/* number of times ShortestPaths() had to be computed (in this thread) */

      #pragma omp for schedule(dynamic,block)
      for(int t=0; t < numsinglenodes; t++){
	int i = singlenodes[t];
	Cnode *Vnode = topnodes[i];
	if(DEBUG>=2) assert(Vnode->edgecnt==1);
	/*	if(Vnode->edgecnt != 1)
		continue;*/

	Cedge *VWedge = &edge[Vnode->edgeid[0]];
	if(DEBUG>=2) assert(VWedge->distance <= MaxSidebranch);
	/*	if(VWedge->distance > MaxSidebranch)
		continue;*/

	Cnode *Wnode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];

	/* Now compute the longest forward and reverse path from Wnode (not going over VWedge or any edge with coverage less than min(MaxCoverage,VWedge->cov)) */
	DATA(GSP,Vnode,SPmark) = 16;/* block this node in ShortestPaths() */
	Cnodeflip Fflip = ShortestPaths(Wnode->nodeid,0,0,0,min((double)MaxCoverage,VWedge->cov),VWedge->distance,GSP,heapT[tid]);
	Cnodeflip Bflip = ShortestPaths(Wnode->nodeid,0,1,0,min((double)MaxCoverage,VWedge->cov),VWedge->distance,GSP,heapT[tid]);
	DATA(GSP,Vnode,SPmark) = 0;
	pathcnt++;

	double forward = GSP[Fflip.nodeid].SPcost[0][Fflip.flip];
	double backward = GSP[Bflip.nodeid].SPcost[1][Bflip.flip];

	//	if(test > 0)	  continue;

	if(VWedge->distance < min(forward,backward) && VWedge->distance <= MaxSidebranch){
          #pragma omp critical
	  {
	    if(VERB>=2 && verb){
	      printf("t=%d/%d,i=%d: Moving single connected node N%d to side chain of N%d(fdist=%0.3f(N%d%c),bdist=%0.3f(N%d%c),edge=%0.3f(%d maps),cycle=%d):tid=%d,paths=%d,wtime=%0.6f\n",
		     t,numsinglenodes,i,nodeorder[Vnode->mapid],nodeorder[Wnode->mapid],forward,nodeorder[node[Fflip.nodeid].mapid],Fflip.flip ? '\'':' ',
		     backward,nodeorder[node[Bflip.nodeid].mapid],Bflip.flip ? '\'':' ',VWedge->distance,edgenodecnt(VWedge,Vnode,0)+2,cyclecnt,tid,pathcnt,wtime());
	      fflush(stdout);
	    }
	    edge_delete(VWedge,Wnode);
	    Vnode->mark = 2;
	    if(VWedge->alignid){/* need to explictly add VWedge to side chain of Wnode */
	      Wnode->chain.push_back(VWedge->edgeid);
            }
	  }
	} else if(VERB>=2 && verb){
          #pragma omp critical
	  {
	    printf("t=%d/%d,i=%d:Single connected node N%d (connected to N%d):fdist=%0.3f(N%d%c),bdist=%0.3f(N%d%c),edge=%0.3f,cycle=%d:tid=%d,paths=%d,wtime=%0.6f\n",
		   t,numsinglenodes,i,nodeorder[Vnode->mapid],nodeorder[Wnode->mapid],forward,nodeorder[node[Fflip.nodeid].mapid],Fflip.flip ? '\'':' ',
		   backward,nodeorder[node[Bflip.nodeid].mapid],Bflip.flip ? '\'':' ',VWedge->distance,cyclecnt,tid,pathcnt,wtime());
	    fflush(stdout);
	  }
	}
      } // for i = 0 .. numtopnodes-1

      delete [] GSPT[tid];
      delete [] heapT[tid];
    } // parallel

    if(VERB>=2){
      printf("Completed multithreaded section for sidechain analysis: nthreads=%d,block=%d,numtopnodes=%d,numsinglenodes=%d:mtime= %0.6f, wtime= %0.6f\n",nthreads,block,numtopnodes,numsinglenodes,mtime(),wtime());
      fflush(stdout);
    }
    delete [] singlenodes;

    //    } // for test = 
    //    exit(1);

    long long Ecnt2 = GCedges(finalcnt,edge);
    int Ncnt2 = GCnodes();
    if(VERB && (Ecnt2 > 0 || Ncnt2 > 0)){
      printf("graph_prune(%d,%d): moved %d single connected nodes to side chain (removed %lld edges) (edges=%d,nodes=%d):time=%0.6f(total=%0.6f)\n",
	     sidechain, verb, Ncnt2, Ecnt2/2, finalcnt, numtopnodes,wtime(),mtime());
      fflush(stdout);
    }
    Ecnt += Ecnt2;
    Ncnt += Ncnt2;
  }

  if(VERB>=2){/* locate highest coverage edge */
    double maxCov = 0.0;
    Cedge *maxEdge = 0;
    for(int t=0; t < numtopnodes; t++){
      Cnode *Tnode = topnodes[t];
      for(int u = 0; u < Tnode->edgecnt; u++){
	Cedge *TUedge = &edge[Tnode->edgeid[u]];
	if(DEBUG) assert(TUedge->cov >= 0.0);
	if(TUedge->cov > maxCov){
	  maxCov = TUedge->cov;
	  maxEdge = TUedge;
	}
      }	
    }
    printf("graph_prune(End): Largest coverage = %0.2f : (N%d,N%d,%0.3f,%d,%d)\n",
	   maxCov, nodeorder[node[maxEdge->Lmap].mapid],nodeorder[node[maxEdge->Rmap].mapid],maxEdge->distance,
	   maxEdge->Lflipped,maxEdge->Rflipped);
    fflush(stdout);
  }

  return Ecnt+Ncnt;
}

/** mark all nodes reachable from current node (any direction) 
 mark will be incremented from 0 to 1 and 2 
 returns number of nodes found */
static int connectedDFS(Cnode *Vnode)
{
  if(Vnode->mark)/* this node already explored (either ancestor or finished) */
    return 0;

  int result = 1;

  Vnode->mark += 1;/* mark Vnode as active */
  if(DEBUG>=2) assert(Vnode->mark == 1);

  /* recurse over all children */
  for(int vw = 0; vw < Vnode->edgecnt; vw++){
    Cedge *VWedge = &edge[Vnode->edgeid[vw]];
    Cnode *Wnode = &node[(VWedge->Lmap == Vnode->nodeid) ? VWedge->Rmap : VWedge->Lmap];
    result += connectedDFS(Wnode);
  }

  Vnode->mark = 2;/* mark source as finished */
  if(DEBUG>=2) assert(Vnode->mark == 2);
  
  return result;
}

static int graph_components_called = 0;

/** find connected components : delete all but one representative node from each connected component 
  Calls connectedDFS() from nodes then searches for ShortestPaths(), finally output_contigs().
*/
void graph_components()
{
  graph_components_called++;

  /* reset all node's mark (and SPmark) to 0 */
  for(int i = 0; i < numtopnodes; i++){
    Cnode *Vnode = topnodes[i];
    Vnode->mark = Vnode->SPmark = 0;
    if(DEBUG>=2){
      unsigned int *pedgeid = Vnode->edgeid;

      for(int j = Vnode->edgecnt; --j >= 0; pedgeid++)
	assert(edge[*pedgeid].valid == 1);
    }
  }

  cyclewarning = 0;
  pathcrosswarning = 0;

  if(!PosEdgeInit)
    PosEdge();

  int numcontigpath = 0;
  Ccontigpath *contigpath = new Ccontigpath[numtopnodes];

  int components = 0;/* number of graph components */
  int singletons = 0;/* number of singleton node graph components */
  CSPresult *LongestPath = new CSPresult[numtopnodes*2];
  if(DEBUG>=2)
    for(int k = numtopnodes*2; --k >= 0;)
      LongestPath[k].SourceNode = -1;

  Cnode **CCnodes = topnodes;
  int CCcnt = numtopnodes;
  int *CCnodeid = NULL;
  int *CCEnodeid = NULL, CCedgecnt = 0, CCedgesiz = 0;
  if(CC_COMPRESS){
    CCnodeid = new int[numnodes];
    if(DEBUG>=2)
      for(int i = 0; i < numnodes; i++)
	CCnodeid[i] = -1;
    int edgecnt = 0;
    for(int k = 0; k < numtopnodes; k++)
      edgecnt += topnodes[k]->edgecnt;
    CCEnodeid = new int[edgecnt];
    CCedgesiz = edgecnt;
  }
  if(CC_FAST){
    CCnodes = new Cnode*[numtopnodes];
    CCcnt = 0;
  }

  CSPdata *GSP = new CSPdata[CC_COMPRESS ? numtopnodes : numnodes];

  CSPdata *allSP[numthreads];
  for(int i = 0; i < numthreads; i++)
    allSP[i] = 0;

  for(int i = 0; i < numtopnodes; i++){
    Cnode *Vnode = topnodes[i];
    if(Vnode->mark >= 2)
      continue;
    if(DEBUG>=2) assert(Vnode->mark == 0);

    int cnt = connectedDFS(Vnode);/* locate connected component with Vnode (marks all nodes in component with mark==2) */

    if(DEBUG) assert(cnt >= 1);
    components++;
    if(cnt == 1)
      singletons++;
    if(DEBUG) assert(Vnode->mark == 2);

    /* collect a list of nodes with mark==2, so this component can be traversed quickly below */
    if(CC_FAST){
      CCcnt = 0;
      if(DEBUG>=2) assert(cnt <= numtopnodes);
      for(int k = 0; k < numtopnodes;k++){
	Cnode *Rnode = topnodes[k];
	if(Rnode->mark != 2)
	  continue;/* Rnode is not in this graph component */
	if(DEBUG>=2) assert(!(Rnode->SPmark & 16));
	CCnodes[CCcnt++] = Rnode;
      }    
      if(DEBUG>=2) assert(CCcnt == cnt);
#if 0 // This block displays the connected graph component sizes and exits
	printf("component=%d,Vnode=%d:cumulative nodes=%d/%d\n",components,Vnode->nodeid,CCcnt,numtopnodes);
	fflush(stdout);
	continue;
#endif

    }

    if(CC_COMPRESS)     /* create mapping from node[] index to CCnodes[] index */
      for(int t = 0; t < CCcnt; t++)
	CCnodeid[CCnodes[t]->nodeid] = t;

    Cnode **pn = CCnodes;
    for(int k = 0; k < CCcnt; k++){
      Cnode *Rnode = pn[k];
      if(CC_COMPRESS){
        GSP[k].SPmark = 0;
	GSP[k].hindex = Rnode->edgecnt;
	GSP[k].Enodeid = &CCEnodeid[CCedgecnt];
	CCedgecnt += Rnode->edgecnt;
	if(DEBUG) assert(CCedgecnt <= CCedgesiz);
	for(int t = 0; t < Rnode->edgecnt; t++){
          int nodeid = GSP[k].Enodeid[t] = CCnodeid[Rnode->Enodeid[t]];
	  if(DEBUG>=2) assert(0 <= nodeid && nodeid < CCcnt && CCnodes[nodeid]->nodeid == Rnode->Enodeid[t]);
	}
      } else {
        GSP[Rnode->nodeid].SPmark = 0;
	GSP[Rnode->nodeid].hindex = Rnode->edgecnt;
	GSP[Rnode->nodeid].Enodeid = Rnode->Enodeid;
      }
    }

    int order = 1;/* count of paths extracted from current graph component */

    int nthreads = numthreads;
    if(CCcnt < nthreads)
      nthreads = max(1,CCcnt);

    #pragma omp parallel num_threads(nthreads)
    {
      int tid = 0;
      #ifdef _OPENMP
      tid = omp_get_thread_num ();
      #endif

      /* allocate/initialize thread-local memory */
      if(!allSP[tid])
        allSP[tid] = new CSPdata[CC_COMPRESS ? numtopnodes : numnodes];

      CSPdata *SP = allSP[tid];
      for(int k = 0; k < CCcnt; k++){
	Cnode *Rnode = pn[k];
	if(CC_COMPRESS){
	  SP[k].hindex = Rnode->edgecnt;
	  SP[k].Enodeid = GSP[k].Enodeid;
	} else {
	  SP[Rnode->nodeid].hindex = Rnode->edgecnt;
	  SP[Rnode->nodeid].Enodeid = Rnode->Enodeid;
	}
      }
    }

    while(1){
      /* try to find the pair of nodes for which the shortest path is the largest of all possible node pairs in this component, 
	 excluding edges in previously generated paths (edges with valid <= 0) and (if FAST_PATH) exclude nodes with (GSP[nodeid].SPmark & 16) */
      double pairdist = 0.0;
      double paircost = -1.0;
      Cnode *pairNode = 0;
      Cnodeflip Endnode;
      Endnode.nodeid = -1;
      Endnode.flip = 0;
      int pairdir = 0, pairk = -1;
      
      int block = 1;
      //      int block = min(16,max(1,CCcnt/(8*nthreads)));
      while(block < 16 && CCcnt > numthreads * block * 8)
	block *= 2;

      if(VERB>=2){
	printf("component=%d(paths=%d):Starting parallel section using %d threads,block=%d(CCcnt=%d,i=%d,numtopnodes=%d):wtime= %0.6f, mtime= %0.6f\n",components,order,nthreads,block,CCcnt,i,numtopnodes, wtime(),mtime());
	fflush(stdout);
      }
      //      if(components >= 55) exit(1);

      cyclecnt = 0;// NOTE : global cyclecnt

      #pragma omp parallel num_threads(nthreads)
      {
	int tid = 0;
	#ifdef _OPENMP
	tid = omp_get_thread_num ();
	#endif

	/* allocate/initialize thread-local memory */
	if(DEBUG) assert(allSP[tid] != NULL);

	CSPdata *SP = allSP[tid];
	for(int k = 0; k < CCcnt; k++){
	  Cnode *Rnode = pn[k];
	  if(CC_COMPRESS)
	    SP[k].SPmark = 0;
	  else
	    SP[Rnode->nodeid].SPmark = 0;
	}

        #pragma omp for schedule(dynamic,block)
	for(int k = 0; k < CCcnt; k++){
	  Cnode *Rnode = pn[k];

	  if(!CC_FAST && Rnode->mark != 2)
	    continue;/* Rnode is not in this graph component */

	  int Rnodeid = CC_COMPRESS ? k : Rnode->nodeid;

	  if(DEBUG>=2 && CC_FAST && !(Rnode->mark==2)){
	    #pragma omp critical
	    {
	      printf("k=%d/%d:Rnode->mark=%d,CCnodes[%d]->mark=%d,CCcnt=%d,tid=%d/%d\n",
		     k,CCcnt,Rnode->mark,k,CCnodes[k]->mark,CCcnt,tid,nthreads);
	      printf("components=%d:Rnode=x%p, CCnodes[k]=x%p,pn=x%p,CCnodes=x%p\n",components,Rnode,CCnodes[k],pn,CCnodes);
	      fflush(stdout);
	      assert(Rnode->mark==2);
	    }
	  }
	  if(DEBUG>=2) assert(!(Rnode->SPmark & 16));
	  if(DEBUG>=2 && (SP[Rnodeid].SPmark & 16)){
	    printf("k=%d/%d:Rnode->mark=%d,Rnode->nodeid=%d,SP[Rnodeid=%d].SPmark=%d\n",
	       k,CCcnt,Rnode->mark,Rnode->nodeid,Rnodeid,SP[Rnodeid].SPmark);
	    fflush(stdout);
	    assert(!(SP[Rnodeid].SPmark & 16));
	  }

	  /* locate most distant node forward/backwards from Rnode */
	  for(int dir = 0; dir <= 1; dir++){
            Cnodeflip Lnodeflip = ShortestPaths(Rnodeid, 0, dir, 2, 0.0, 0.0, SP, 0, CCcnt, CC_COMPRESS ? CCnodes : 0 , CC_COMPRESS ? CCnodeid : 0);

	    CSPresult *pSP = &LongestPath[k*2 + dir];
	    pSP->SourceNode = k;
	    pSP->dir = dir;
	    pSP->EndNode = Lnodeflip;
	    pSP->pairdist = SP[Lnodeflip.nodeid].SPdist[dir][Lnodeflip.flip];
	    pSP->paircost = SP[Lnodeflip.nodeid].SPcost[dir][Lnodeflip.flip];

	    if(DEBUG>=3){
	      for(int t = 0; t < CCcnt; t++){
		Cnode *Tnode = pn[t];
		if(DEBUG && (SP[CC_COMPRESS ? t : Tnode->nodeid].SPmark & 16)){
		  printf("k=%d,dir=%d:t=%d:Tnode->mark=%d,Tnode->nodeid=%d,SP[Tnode->nodeid].SPmark=%d\n",
			 k,dir,t,Tnode->mark,Tnode->nodeid,SP[CC_COMPRESS ? t : Tnode->nodeid].SPmark);
		  fflush(stdout);
		  assert(!(SP[CC_COMPRESS ? t : Tnode->nodeid].SPmark & 16));
		}
	      }
	    }
	  } // dir = 0 .. 1
	  if(VERB>=3){
	    printf("component=%d(paths=%d):k=%d/%d,tid=%d/%d: wtime= %0.6f\n",components,order,k,CCcnt,tid,nthreads,wtime());
	    fflush(stdout);
	  }
	} // omp for k = 0 ... CCcnt - 1

	/* free-up thread-local memory */
	if(VERB>=2){
#pragma omp critical 
	  {
	    printf("tid=%d:thread termination: SP=x%p will be deleted\n",tid,SP);
	    fflush(stdout);
	  }
	}
      } // omp parallel

      if(VERB>=2){
	printf("component=%d(paths=%d):Completed parallel section using %d threads:wtime= %0.6f, mtime= %0.6f\n",components,order,nthreads,wtime(),mtime());
	fflush(stdout);
      }

      /* locate best overall path */
      for(int k = 0; k < CCcnt; k++){
	Cnode *Rnode = pn[k];
	if(!CC_FAST && Rnode->mark != 2)
	  continue;/* Rnode is not in this graph component */
	for(int dir = 0; dir <= 1; dir++){
	  CSPresult *pSP = &LongestPath[k*2 + dir];
	  if(DEBUG>=2) assert(pSP->SourceNode >= 0 && pSP->SourceNode < numtopnodes);
	  if(pSP->pairdist > pairdist){
	    pairdist = pSP->pairdist;
	    paircost = pSP->paircost;
	    pairk = pSP->SourceNode;
	    pairNode = Rnode;
	    Endnode = pSP->EndNode;
	    pairdir = pSP->dir;
	  }
	}
      }      

      if(DEBUG>=2 && order <= 1 && !(pairdist > 0.0)){
	printf("Component %d(%d nodes), path #%d: no path found!\n",  components, cnt, order);
	for(int k = 0; k < CCcnt; k++){
	  Cnode *Rnode = pn[k];
	  if(!CC_FAST && Rnode->mark != 2)
	    continue;/* Rnode is not in this graph component */
	  printf("Rnode=N%d:mark=%d,SPmark=%d\n", nodeorder[Rnode->mapid],Rnode->mark,Rnode->SPmark);
	  for(int rp = 0; rp < Rnode->edgecnt; rp++){
	    Cedge *RPedge = &edge[Rnode->edgeid[rp]];
	    printf("  rp=%d:edge=(N%d,N%d,%0.3f,%d,%d),cov=%0.1f\n", rp,
		   nodeorder[node[RPedge->Lmap].mapid],nodeorder[node[RPedge->Rmap].mapid],
		   RPedge->distance, RPedge->Lflipped,RPedge->Rflipped, RPedge->cov);
	  }
	}
	fflush(stdout);
	assert(pairdist > 0.0);
      }

      if(VERB>=2){
	printf("component=%d(paths=%d):Best remaining path: pairdist= %0.3f, ContigMinLen= %0.3f:wtime= %0.6f, mtime= %0.6f\n",components,order,pairdist,ContigMinLen,wtime(),mtime());
	fflush(stdout);
      }

      if(pairdist <= 0.0 || (order > 1 && pairdist <= ContigMinLen))
	break;
      
      double distance = 0.0;
      int mapcnt = 0;
      if(DEBUG>=2) assert(!(pairNode->SPmark & 16));

      if(DEBUG>=2)
	for(int k = 0; k < CCcnt; k++){
	  Cnode *Rnode = pn[k];
	  assert(GSP[CC_COMPRESS ? k : Rnode->nodeid].SPmark == 0);
	}

      Cnodeflip nodeflip = ShortestPaths(CC_COMPRESS ? pairk : pairNode->nodeid, 0, pairdir, 2, 0.0, 0.0, GSP, 0, CCcnt, CC_COMPRESS ? CCnodes : 0, CC_COMPRESS ? CCnodeid : 0);
      if(DEBUG>=2) assert(Endnode.nodeid == nodeflip.nodeid && Endnode.flip == nodeflip.flip);
      double origpairdist = pairdist;
      pairdist = GSP[nodeflip.nodeid].SPdist[pairdir][nodeflip.flip];
      if(DEBUG>=2) assert(fabs(pairdist - origpairdist) <= 0.01);

      int origorder = order;
      int SPcnt = 0;/* count of calls to ShortestPaths() */

      if(VERB>=2){
	printf("component=%d(paths=%d):Confirmed best remaining path: pairdist= %0.3f, ContigMinLen= %0.3f:wtime= %0.6f, mtime= %0.6f\n",components,order,pairdist,ContigMinLen,wtime(),mtime());
	fflush(stdout);
      }

      while(1){/* This loop has only one iteration, unless FAST_PATH > 0 */

	/* now look for next longest path : first mark the current best path edges with valid = 0 so it is ignored in the future */
	/* also save the current path to the output contiglist if it involves at least ContigMinMaps different maps */

	/* initialize next contig list */
	Ccontigpath *path = &contigpath[numcontigpath];
	if(!path->nodeid){/* NOTE: cannot use CCcnt as array size since compound edges may include additional maps */
	  path->nodeid = new int[numnodes];
	  path->align = new Calign *[numnodes];
	  path->flip = new int[numnodes];
	}
	
	Cnode *Lnode = CC_COMPRESS ? CCnodes[nodeflip.nodeid] : &node[nodeflip.nodeid];

	path->nodeid[mapcnt] = Lnode->nodeid;
	path->flip[mapcnt] = nodeflip.flip;

	int rcnt = 0;/* count of top level edges (not including sub-edges of top level compound edges) */

	for(;Lnode->nodeid != pairNode->nodeid;){
	  if(FAST_PATH)
	    GSP[nodeflip.nodeid].SPmark |= 16;
	  Cedge *Sedge = GSP[nodeflip.nodeid].SPedge[pairdir][nodeflip.flip];
	  if(DEBUG>=2) assert(Sedge != 0);
	  if(DEBUG>=2) assert(Sedge->valid == 1);
	  Sedge->valid = 0;
	  rcnt++;

	  distance += Sedge->distance;

	  if(VERB >= PVERB && !(cnt <= 2 && pairNode->edgecnt==1 && edge[pairNode->edgeid[0]].alignid))
	    printf("N%d%c %lld(D=%0.3f%c,C=%0.1f):Sedge=(N%d,N%d,%d,%d)%c\n", nodeorder[Lnode->mapid], nodeflip.flip ? '\'':' ',Gmap[Lnode->mapid]->id,
		   Sedge->distance,Sedge->alignid ? ' ' : '*', Sedge->cov,
		   nodeorder[node[Sedge->Lmap].mapid],nodeorder[node[Sedge->Rmap].mapid],Sedge->Lflipped,Sedge->Rflipped,Sedge->alignid ? ' ':':');
	  mapcnt += edgedisplay(Sedge,Lnode,nodeflip.flip,path,mapcnt,VERB);
	  mapcnt++;

	  Cnode *Nnode = &node[(Sedge->Lmap == Lnode->nodeid) ? Sedge->Rmap : Sedge->Lmap];
	  int Nflipped;
	  if(pairdir)
	    Nflipped = nodeflip.flip ? (!Sedge->Lflipped) : Sedge->Rflipped;
	  else
	    Nflipped = nodeflip.flip ? (!Sedge->Rflipped) : Sedge->Lflipped;
	  if(DEBUG>=2 && !(Nflipped == (nodeflip.flip ^ (Sedge->Lflipped | Sedge->Rflipped)))){
	    printf("\nN%d%c(D=%0.3f%c,C=%0.1f):edge=(N%d,N%d,%d,%d)%c\n", nodeorder[Lnode->mapid], nodeflip.flip ? '\'':' ',
		   Sedge->distance,Sedge->alignid ? ' ' : '*', Sedge->cov,
		   nodeorder[node[Sedge->Lmap].mapid],nodeorder[node[Sedge->Rmap].mapid],Sedge->Lflipped,Sedge->Rflipped,Sedge->alignid ? ' ':':');
	    printf("  Nflipped=%d,nodeflip.flip=%d,Sedge->Lflipped=%d,Sedge->Rflipped=%d,pairdir=%d\n",
		   Nflipped,nodeflip.flip,Sedge->Lflipped,Sedge->Rflipped,pairdir);
	    fflush(stdout);
	    assert(Nflipped == (nodeflip.flip ^ (Sedge->Lflipped | Sedge->Rflipped)));
	  }
	
	  /* update contig list */
	  path->nodeid[mapcnt] = Nnode->nodeid;
	  path->flip[mapcnt] = Nflipped;

	  nodeflip.nodeid = CC_COMPRESS ? CCnodeid[Nnode->nodeid] : Nnode->nodeid;
	  nodeflip.flip = Nflipped;
	  if(DEBUG>=2){
	    Lnode = CC_COMPRESS ? CCnodes[nodeflip.nodeid] : &node[nodeflip.nodeid];
	    assert(Lnode == Nnode);
	  }
	  Lnode = Nnode;
        }
	path->align[mapcnt++] = 0;

	if(VERB >= PVERB && !(cnt <= 2 && pairNode->edgecnt==1 && edge[pairNode->edgeid[0]].alignid)){
	  printf("N%d ( total distance = %0.3f kb, maps=%d,block=%d,threads=%d)\n\n",nodeorder[pairNode->mapid],distance,mapcnt,block,nthreads);
	  fflush(stdout);
	}

	if(mapcnt >= ContigMinMaps && (!AVCOV || pairdist/paircost >= ContigMinCov) && distance > ContigMinLen){  /* Finalize/confirm contig list */
	  if(DEBUG>=2) assert(Lnode->nodeid == pairNode->nodeid);
	  path->nummaps = mapcnt;
	  numcontigpath++;
	}

	if(VERB && !(cnt <= 2 && pairNode->edgecnt==1 && edge[pairNode->edgeid[0]].alignid)){
	  int ncnt = 0; /* count remaining nodes */
	  int ecnt = 0; /* count remaining edges */

	  Cnode *Enode = CC_COMPRESS ? CCnodes[Endnode.nodeid] : &node[Endnode.nodeid];

	  if(VERB >= PVERB){
	    for(int k = 0; k < CCcnt; k++){
	      Cnode *Rnode = pn[k];
	      if(!CC_FAST && Rnode->mark != 2)
		continue;/* Rnode is not in this graph component */
	      if((GSP[CC_COMPRESS ? k : Rnode->nodeid].SPmark & 16))
		continue;/* Rnode is in a longer path */
	      ncnt++;

	      for(int rw = 0; rw < Rnode->edgecnt; rw++)
		if(edge[Rnode->edgeid[rw]].valid)
		  ecnt++;
	    }      

	    if(mapcnt >= ContigMinMaps && (!AVCOV || pairdist/paircost >= ContigMinCov) && distance > ContigMinLen)
	      printf("Component %d(%d nodes), path #%d (Consensus Path %d): length=%0.3f(%0.3f) from N%d to N%d%c, maps=%d,rcnt=%d (dir=%d),cost=%0.3f,dist=%0.3f (ncnt=%d,ecnt=%d),i=%d/%d:time=%0.6f(total=%0.6f)\n",
		     components, cnt, order, numcontigpath-1,distance,distance + (pairNode->len + Enode->len)*0.5,
		     nodeorder[pairNode->mapid], nodeorder[Enode->mapid], Endnode.flip ? '\'': ' ', 
		     mapcnt, rcnt, pairdir, paircost,pairdist, ncnt, ecnt, i, numtopnodes, wtime(),mtime());
	    else if(VERB/* >=2 */)
	      printf("Component %d(%d nodes), path #%d: length=%0.3f(%0.3f) from N%d to N%d%c, maps=%d,rcnt=%d (dir=%d),cost=%0.3f,dist=%0.3f: Discarded (ncnt=%d,ecnt=%d), i=%d/%d:time=%0.6f(total=%0.6f)\n",
		     components, cnt, order, distance, distance + (pairNode->len + Enode->len)*0.5,
		     nodeorder[pairNode->mapid], nodeorder[Enode->mapid], Endnode.flip ? '\'': ' ', 
		     mapcnt, rcnt, pairdir, paircost,pairdist, ncnt, ecnt, i, numtopnodes,wtime(),mtime());
	  } else {
	    if(mapcnt >= ContigMinMaps && (!AVCOV || pairdist/paircost >= ContigMinCov) && distance > ContigMinLen)
	      printf("Component %d(%d nodes), path #%d (Consensus Path %d): length=%0.3f(%0.3f) from N%d to N%d%c, maps=%d (dir=%d),cost=%0.3f,dist=%0.3f,i=%d/%d:time=%0.6f(total=%0.6f)\n",
		     components, cnt, order, numcontigpath-1,distance,distance + (pairNode->len + Enode->len)*0.5,
		     nodeorder[pairNode->mapid], nodeorder[Enode->mapid], Endnode.flip ? '\'': ' ', 
		     mapcnt, pairdir, paircost,pairdist, i, numtopnodes, wtime(),mtime());
	    else if(VERB/* >=2 */)
	      printf("Component %d(%d nodes), path #%d: length=%0.3f(%0.3f) from N%d to N%d%c, maps=%d (dir=%d),cost=%0.3f,dist=%0.3f,i=%d/%d: Discarded:time=%0.6f(total=%0.6f)\n",
		     components, cnt, order, distance, distance + (pairNode->len + Enode->len)*0.5,
		     nodeorder[pairNode->mapid], nodeorder[Enode->mapid], Endnode.flip ? '\'': ' ', 
		     mapcnt, pairdir, paircost,pairdist, i, numtopnodes,wtime(),mtime());
	    
	  }
	  fflush(stdout);
	}

	if(order > 1 && distance <= ContigMinLen)
	  break;

	order++;
	if(!FAST_PATH)
	  break;/* recompute longest path from scratch */

	LongestPath[pairk*2 + pairdir].SourceNode = -1;/* mark this source node as checked */

	/* locate the next shortest path from CSPresult[2*(k=0..numnodes-1) + (dir=0..1)] : skip paths that start or end on previous paths with GSP[nodeid].SPmark & 16 */
	/* NOTE : the complete path will be recomputed to skip all nodes with SPmark & 16, and if this path is too short, this loop will terminate */
	pairdist = 0.0;
	while(1){/* repeatedly look for next best path and try to confirm it */
	  for(int k = 0; k < CCcnt; k++){
	    Cnode *Rnode = pn[k];
	    if((!CC_FAST && Rnode->mark != 2) || (GSP[CC_COMPRESS ? k : Rnode->nodeid].SPmark & 16))
	      continue;
	    for(int dir = 0; dir <= 1; dir++){
	      CSPresult *pSP = &LongestPath[k*2 + dir];
	      if(pSP->SourceNode < 0)
		continue;
	      if((GSP[pSP->EndNode.nodeid].SPmark & 16))
		continue;
	      if(pSP->pairdist > pairdist){
		pairdist = pSP->pairdist;
		paircost = pSP->paircost;
		pairk = pSP->SourceNode;
		pairNode = Rnode;
		Endnode = pSP->EndNode;
		pairdir = pSP->dir;
	      }
	    }
	  }
	  origpairdist = pairdist;
	  if(pairdist <= ContigMinLen)
	    break;

	  /* try to confirm if path still exists while avoiding SPmark & 16 && valid==0 */
	  distance = 0.0;
	  mapcnt = 0;

	  SPcnt++;
	  if(DEBUG>=2) assert(!(pairNode->SPmark & 16));
 	  nodeflip = ShortestPaths(CC_COMPRESS ? pairk : pairNode->nodeid, 0, pairdir, 2, 0.0, 0.0, GSP, 0, CCcnt, CC_COMPRESS ? CCnodes : 0, CC_COMPRESS ? CCnodeid : 0);/* recheck path */

	  pairdist = GSP[nodeflip.nodeid].SPdist[pairdir][nodeflip.flip];
	  paircost = GSP[nodeflip.nodeid].SPcost[pairdir][nodeflip.flip];
	  if (FAST_PATH <= 1 || pairdist > ContigMinLen)
	    break;
	  LongestPath[pairk*2 + pairdir].SourceNode = -1;/* mark this source node as checked */
	}/* while(1) : find next best path */

	if(pairdist <= ContigMinLen)
	  break;
      }/* while(1) : FAST_PATH loop */

      if(FAST_PATH && VERB>=2 && order > origorder+1){
	printf("Generated %d paths without recomputing all paths (last pairdist=%0.3f->%0.3f,ContigMinLen=%0.3f,SPcnt=%d):wtime= %0.6f\n",order-origorder,origpairdist,pairdist,ContigMinLen,SPcnt,wtime());
	fflush(stdout);
      }
      if(FAST_PATH){/* remove SPmark & 16 from all nodes */
	for(int k = 0; k < CCcnt; k++){
	  if(CC_COMPRESS)
	    GSP[k].SPmark &= ~16;
	  else {
	    Cnode *Rnode = pn[k];
	    GSP[Rnode->nodeid].SPmark &= ~16;
	  }
	}
      }
    }/* while(1) : enumerate all paths of current connect component */

    Vnode->mark = 0;/* Vnode is the representative node of the current connected component */
    for(int k = numtopnodes; --k >= 0;){
      Cnode *Rnode = topnodes[k];
      if(Rnode->mark == 2)
	Rnode->mark = 3;
    }
    /* find next component */
  }/* for(int i = 0; i < numtopnodes; i++) */
	
  for(int i = 0; i < numthreads; i++)
    delete [] allSP[i];

  delete [] GSP;

  if(CC_FAST){
    delete [] CCnodes;
    CCnodes = topnodes;
    CCcnt = numtopnodes;
  }
  delete [] CCnodeid;
  delete [] CCEnodeid;
  delete [] LongestPath;

  for(int k = numtopnodes; --k >= 0;){
    Cnode *Rnode = topnodes[k];
    if(Rnode->mark == 3)
      Rnode->mark = 2;
  }

  if(VERB){
    printf("graph_components: Found %d graph components (%d singleton nodes)\n",
	   components,singletons);
    fflush(stdout);
  }

  if(MinSidebranch > 0.0){/* output side chains larger than MinSidebranch kb and with number of maps >= ContigMinMaps */
    if(VERB){
      printf("Locating side chaings larger than %0.3f kb and at least %d maps\n",MinSidebranch,ContigMinMaps);
      fflush(stdout);
    }

    int sidechain_cnt = 0;
    int sidechain_output = 0;
    for(int k = 0; k < numnodes; k++){
      Cnode *Rnode = &node[k];
      int Rflipped = 0;
      for(std::vector<int>::reverse_iterator it = Rnode->chain.rbegin() ; it != Rnode->chain.rend(); ++it){ Cedge *Redge = &edge[*it];
	/* sidechains of Rnode */

	Ccontigpath *path = &contigpath[numcontigpath];
	if(!path->nodeid){
	  path->nodeid = new int[numnodes];
	  path->align = new Calign *[numnodes];
	  path->flip = new int[numnodes];
	}
	int cnt = edgenodecnt(Redge,Rnode,Rflipped);
	int mapcnt = 0;
	path->nodeid[mapcnt] = Rnode->nodeid;
	path->flip[mapcnt] = Rflipped;
	
	double distance = Redge->distance;
	if(VERB >= PVERB && cnt >= ContigMinMaps && distance >= MinSidebranch)
	  printf("N%d%c %lld(D=%0.3f%c,C=%0.1f):edge=(N%d,N%d,%d,%d)%c (sidechain)\n", nodeorder[Rnode->mapid], Rflipped ? '\'':' ', Gmap[Rnode->mapid]->id,
		 Redge->distance,Redge->alignid ? ' ' : '*', Redge->cov,
		 nodeorder[node[Redge->Lmap].mapid],nodeorder[node[Redge->Rmap].mapid],Redge->Lflipped,edge->Rflipped,edge->alignid ? ' ':':');
	
	mapcnt = edgedisplay(Redge,Rnode,Rflipped,path,mapcnt, (cnt >= ContigMinMaps && distance >= MinSidebranch) ? VERB : 0);
	
	/* process end node at other end of edge */
	mapcnt++;
	
	Cnode *Nnode = &node[(Redge->Lmap == Rnode->nodeid) ? Redge->Rmap : Redge->Lmap];
	int Nflipped = Redge->Lflipped ^ Redge->Rflipped;
	
	/* update contig list with Nnode */
	path->nodeid[mapcnt] = Nnode->nodeid;
	path->flip[mapcnt] = Nflipped;
	path->align[mapcnt++] = 0;

	if(VERB >= PVERB && cnt >= ContigMinMaps && distance >= MinSidebranch){
	  printf("N%d (total distance = %0.3f kb, maps=%d) (sidechain)\n\n", nodeorder[Nnode->mapid],distance,mapcnt);
	  fflush(stdout);
	}

	sidechain_cnt++;

	if(mapcnt >= ContigMinMaps && distance >= MinSidebranch){/* Finalize/confirm contig list */
	  sidechain_output++;

	  path->nummaps = mapcnt;
	  numcontigpath++;
	  
	  if(VERB){
	    printf("Component %d(%d nodes),Rnode=N%d,sidechain=(N%d,N%d,%d,%d) (Consensus Path %d): length=%0.3f from N%d to N%d%c, maps=%d:time=%0.6f(total=%0.6f)\n",
		   components, mapcnt, nodeorder[Rnode->mapid], nodeorder[node[Redge->Lmap].mapid],nodeorder[node[Redge->Rmap].mapid],Redge->Lflipped,Redge->Rflipped,
		   numcontigpath-1,distance,nodeorder[node[path->nodeid[0]].mapid],nodeorder[node[path->nodeid[mapcnt-1]].mapid],path->flip[mapcnt-1] ? '\'':' ',mapcnt,wtime(),mtime());
	    if(VERB>=PVERB){
	      for(int t = 0; t < mapcnt; t++){
		if(t < mapcnt-1)
		  printf("N%d%c(id=%lld):align=(N%d,N%d,%d)\n",
			 nodeorder[node[path->nodeid[t]].mapid],path->flip[t] ? '\'':' ', Gmap[node[path->nodeid[t]].mapid]->id,
			 nodeorder[path->align[t]->mapid1],nodeorder[path->align[t]->mapid2],path->align[t]->orientation);
		
		else
		  printf("N%d%c(id=%lld)\n", nodeorder[node[path->nodeid[t]].mapid],path->flip[t] ? '\'':' ', Gmap[node[path->nodeid[t]].mapid]->id);
	      }
	    }
	    fflush(stdout);
	  }
	}
      } /* sidechains of Rnode */
    } /* Rnode = &node[k = 0 .. numnodes-1] */
    if(VERB){
      printf("graph_components: Found %d sidechains (%d with at least %d maps and distance >= %0.1f kb)\n",sidechain_cnt,sidechain_output,ContigMinMaps,MinSidebranch);
      fflush(stdout);
    }
  } /* MinSidebranch */

  if(VERB && numcontigpath <= 0){
    printf("WARNING: no contigs generated : possibly due to insufficient data\n");
    fflush(stdout);
  }
  if(numcontigpath >= 1){
    /* build a contig contig[i] for each contigpath[i=0..numcontigpath-1] */
    Ccontig *contig = new Ccontig[numcontigpath];
    int maxmaps = 0;
    for(int i = 0; i < numcontigpath; i++){
      int nummaps = contigpath[i].nummaps;
      if(nummaps > maxmaps)
	maxmaps = nummaps;
    }
    Ccontig *contiglist = new Ccontig[maxmaps];
    for(int i = 0; i < numcontigpath; i++){
      //      FILE *fp = NULL;/* debugging fp */
      //      char filename[PATH_MAX];/* debugging file name */
      int rverb = 0;
      /*      if(i==2617 || i==2615)
	      rverb = 1;*/
      if(rverb){/* open debugging file */

      }

      /* copy the maps into an array of contigs */
      Ccontigpath *path = &contigpath[i];
      int mapcnt = path->nummaps;

#if 0
#define MAXMAP 4
      if(DEBUG>=2 && VERB>=2 && i==1 && mapcnt > MAXMAP){// copy only maps 1..MAXMAP */
	printf("Assembling Consensus Path %d:mapcnt was %d\n",i,mapcnt);
	for(int t = 0; t < MAXMAP; t++){
	  path->nodeid[t] = path->nodeid[t+1];
	  path->align[t] = path->align[t+1];
	  path->flip[t] = path->flip[t+1];
	}
	path->nummaps = mapcnt = MAXMAP;
      }
#endif

      if(VERB){
	printf("Assembling Consensus Path %d:mapcnt=%d:\n",i,mapcnt);
	fflush(stdout);
      }
      for(int cnt = 0; cnt < mapcnt;cnt++){
	if(VERB>=PVERB){
	  if(cnt < mapcnt-1)
	    printf("N%d%c(id=%lld):align=(N%d,N%d,%d,PV=%0.1f,score=%0.2f,A=%d)\n",
		   nodeorder[node[path->nodeid[cnt]].mapid],path->flip[cnt] ? '\'':' ', Gmap[node[path->nodeid[cnt]].mapid]->id,
		   nodeorder[path->align[cnt]->mapid1],nodeorder[path->align[cnt]->mapid2],path->align[cnt]->orientation,path->align[cnt]->logPV,path->align[cnt]->score,path->align[cnt]->numpairs);
	  else
	    printf("N%d%c(id=%lld)\n", nodeorder[node[path->nodeid[cnt]].mapid],path->flip[cnt] ? '\'':' ', Gmap[node[path->nodeid[cnt]].mapid]->id);
	  fflush(stdout);

	  /*	  if(DEBUG && Gmap[node[path->nodeid[cnt]].mapid]->id == 3221000210LL)
		  rverb = 1;*/
	}
	Ccontig *p = new Ccontig(&node[contigpath[i].nodeid[cnt]]);
	if(DEBUG>=2) assert(p->X[0] == 0);
	if(DEBUG) assert(p->scale == 0);
	if(DEBUG) assert(p->mapid < 0);
	contiglist[cnt] = *p;
	p->numsite[0] = p->nummaps = 0;
	delete p;
      }

      /* repeatedly merge pairs of contigs in contigslist[0..mapcnt-1], until only 1 contig is left */
      while(mapcnt > 1){
	int j,k;
	for(k=j=0; j < mapcnt; k++, j += 2){
	  if(j+1 < mapcnt){ /* combine contig j and j+1 and save the result as contig i */
	    if(VERB && (VERB>=2 || rverb)){
	      printf("i=%d/%d:j=%d/%d,k=%d:Merging contiglist[j] and contiglist[j+1] to contiglist[k] using align[j](N%d(%lld),N%d(%lld),%d):",
		     i,numcontigpath,j,mapcnt,k,nodeorder[path->align[j]->mapid1],Gmap[path->align[j]->mapid1]->id,nodeorder[path->align[j]->mapid2],Gmap[path->align[j]->mapid2]->id,path->align[j]->orientation);
	      if(j+2 < mapcnt)
		printf("align[j+1]=(N%d(%lld),N%d(%lld),%d)",
		       nodeorder[path->align[j+1]->mapid1],Gmap[path->align[j+1]->mapid1]->id,nodeorder[path->align[j+1]->mapid2],Gmap[path->align[j+1]->mapid2]->id,path->align[j+1]->orientation);
	      printf("\n");
	      fflush(stdout);
	    }
	    Ccontig *merged = new Ccontig(&contiglist[j],&contiglist[j+1],path->align[j],rverb);
	    if(DEBUG>=2) assert(merged->X[0] == 0);
	    contiglist[k] = *merged;
	    merged->numsite[0] = merged->nummaps = 0;
	    delete merged;

	    if(VERB && (VERB>=2 || rverb)){
	      printf("contiglist[k=%d]:\n",k);
	      fflush(stdout);
	      contiglist[k].flatten();

	      contiglist[k].display();
	    }

	    path->align[k] = path->align[j+1];/* link between j+1 and j+2 */
	  } else {/* just copy contig j */
	    if(VERB && (VERB>=2 || rverb)){
	      printf("i=%d/%d:j=%d/%d,k=%d:moving contiglist[j] to contiglist[k]\n",i,numcontigpath,j,mapcnt,k);
	      fflush(stdout);
	    }
	    contiglist[k] = contiglist[j];
	    contiglist[j].numsite[0] = contiglist[j].nummaps = 0;// mark as empty
	  }
	}
	mapcnt = k;
      }
      contig[i] = contiglist[0];
      contiglist[0].numsite[0] = contiglist[0].nummaps = 0;// mark as empty 

      contig[i].flatten();

      if(VERB && (VERB>=2 || rverb)){
	//	contig[i].contigflip();
	printf("Assembled Consensus Path %d:\n",i);
	contig[i].display();
      }

      /* pre-compute X[][] for contig[i] */
      if(DEBUG>=2) assert(contig[i].X[0] ==0);
      contig[i].Xmap();
    }

    /* sort contig[0..numcontigpath-1] from largest to smallest contig (total number of sites in contig) */
    for(int i = 0; i < numcontigpath; i++){
      contig[i].totalsites = -1;/* actual value is computed (if not already done) by ContigTotsitesDec() */
      contig[i].contigid = i;
    }
    qsort(contig,numcontigpath,sizeof(Ccontig),(intcmp*)ContigTotsitesDec);

    /* output contig[0..numcontigpath-1] */
    if(VERB){
      for(int i = 0; i < numcontigpath; i++){
	printf("Unrefined contig[%d] corresponds to filename suffix unrefined_contig%d.cmap and Consensus Path %d:\n",
	       i, i+1, contig[i].contigid);

	if(VERB>=2 && (contig[i].contigid==2617 || contig[i].contigid==2615)){
	  printf("Final Consensus Path %d:\n",contig[i].contigid);
	  contig[i].display();
	}
      }
      fflush(stdout);
    }
    output_contigs(contig,numcontigpath,output_prefix);
    /*    if(contig_format){// output text format as well
      char prefix[PATH_MAX];
      sprintf(prefix,"%s.txt",output_prefix);

      contig_format = 0;
      output_contigs(contig,numcontigpath,prefix);
      contig_format = 1;
      }*/

    for(int i = 0; i < numcontigpath; i++){
      if(Refine){ /* Refine consensus map of contig[i] */
	if(VERB){
	  printf("Refining contig[%d] (Consensus Path %d):\n",i, contig[i].contigid);
	  fflush(stdout);
	}
	if(SAVE_MAPSET){/* save mapset as <PREFIX>_contig<N>.bnx */
	  Cmap **goodmaps = new Cmap *[contig[i].nummaps];
	  for(int k = 0; k < contig[i].nummaps; k++)
	    goodmaps[k] = Gmap[contig[i].contig[k].mapid];
	  char prefix[PATH_MAX];
	  sprintf(prefix,"%s_contig%d",output_prefix,i+1);
	  output_bnx(prefix,goodmaps,0,contig[i].nummaps-1,0,NULL,NULL,-1);
	  delete [] goodmaps;
	}
	int Lfrozen = 0,Rfrozen = 0;
	contig[i].id = i+1;
	refine(&contig[i],1,Lfrozen,Rfrozen);
      }

      if(VERB>=2 && contig[i].nummaps > 0){/* display contig[i] */
	printf("contig %d:numsite=%d:\n",i,contig[i].numsite[0]);
	int k;
	for(k = 1; k <= contig[i].numsite[0]; k++)
	  printf("   k=%d:site[k]=%0.3f,fragcnt[k-1]=%0.1f,sitecnt[k]=%0.1f,sitecntFN[k]=%0.1f\n",
		 k,contig[i].site[0][k],contig[i].fragcnt[0][k-1],contig[i].sitecnt[0][k],contig[i].sitecntFN[0][k]);
	printf("   k=%d:site[k]=%0.3f,fragcnt[k-1]=%0.1f\n", k, contig[i].site[0][k], contig[i].fragcnt[0][k-1]);
	fflush(stdout);
      }

      /* output draft consensus for contig[i] (as prefix_contig<i+1>.cmap file) */
      output_draft(&contig[i],(long long)(i+1),output_prefix,0,0,0);
    }

    delete [] contiglist;

    /* output the trimmed Maps */
    //  output_TrimmedMaps(nummaps,Gmap,node,output_prefix);

    /*    for(int i = 0; i < numcontigpath; i++)
	  output_draft(&contig[i],(long long)(i+1),output_prefix);*/

    delete [] contig;
  }
  delete [] contigpath;

  if(DEBUG){
    /* reset all edge's valid to 1 : so this function can be called again */
    for(int i = 0; i < numtopnodes; i++){
      Cnode *Vnode = topnodes[i];
      unsigned int *pedgeid = Vnode->edgeid;
      for(int j = Vnode->edgecnt; --j >= 0; pedgeid++)
	edge[*pedgeid].valid = 1;
    }
  }
}

/** free up all memory in graph */
void graph_free()
{
  if(topnodes){
    delete [] topnodes;
    topnodes = NULL;
    numtopnodes = 0;
  }

  delete [] node;

#if 0
  Cedge *pedge = &edge[1];
  for(unsigned int i= 1; i <= Numedges; i++, pedge++)
    pedge->allfree();
#endif

  delete [] edge;
  delete [] nodemap;
  delete [] nodeorder;
  if(edgelist)
    delete [] edgelist;
  if(EnodeidMem)
    delete [] EnodeidMem;
  if(edgeidMem)
    delete [] edgeidMem;
  if(map2nodeid)
    delete [] map2nodeid;
}

void graph_compact()
{
  double stime = wtime();

  size_t maxedgeMem = 0;
  for(int i = 0; i < numnodes; i++)
    maxedgeMem += node[i].maxedge;

  int newnumnodes = numtopnodes;
  
  size_t newedgeMem = 0;
  for(int i = 0; i < numtopnodes; i++)
    newedgeMem += topnodes[i]->edgecnt;

  if(DEBUG && !((newedgeMem % 2) == 0)){
    size_t edgeCnt = 0;
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edge[i].valid < 0)
	continue;
      assert(edge[i].mark == 0);
      edgeCnt++;
    }
    printf("newedgeMem = %lu, edgeCnt = %lu\n", newedgeMem, edgeCnt);
    fflush(stdout);

    for(int i = 0; i < numtopnodes; i++){
      Cnode *Vnode  = topnodes[i];
      for(int k = 0; k < Vnode->edgecnt; k++){
        unsigned int edgeid = Vnode->edgeid[k];
	assert(1 <= edgeid && edgeid <= Numedges);
	edge[edgeid].mark++;
      }
    }
    for(unsigned int i = 1; i <= Numedges;i++){
      if(edge[i].valid < 0)
	continue;
      if(edge[i].mark != 2){
	printf("edge[%u]:edgeid=%d,mark=%d,valid=%d, Lmap=%d, Rmap=%d\n",i,edge[i].edgeid,edge[i].mark,edge[i].valid,edge[i].Lmap,edge[i].Rmap);
	fflush(stdout);
	
	assert(0 <= edge[i].Lmap && edge[i].Lmap < numnodes);
	assert(0 <= edge[i].Rmap && edge[i].Rmap < numnodes);
	
	Cnode *Lnode = &node[edge[i].Lmap];
	printf("node[Lmap].edgecnt=%d/%d, mark=%d:\n",Lnode->edgecnt,Lnode->maxedge,Lnode->mark);
	for(int k = 0; k < Lnode->edgecnt; k++)
	  printf("k=%d: Lnode->edgeid[k] = %u\n", k, Lnode->edgeid[k]);
	fflush(stdout);

	Cnode *Rnode = &node[edge[i].Rmap];
	printf("node[Rmap].edgecnt=%d/%d, mark=%d:\n",Rnode->edgecnt,Rnode->maxedge,Rnode->mark);
	for(int k = 0; k < Rnode->edgecnt; k++)
	  printf("k=%d: Rnode->edgeid[k] = %u\n", k, Rnode->edgeid[k]);
	fflush(stdout);

	assert(edge[i].mark == 2);
      }
    }

    assert((newedgeMem % 2) == 0);
  }
  if(DEBUG)assert((newedgeMem/2) < MASK_LL(32));

  unsigned int newNumedges = newedgeMem / 2;
  unsigned int newMaxedges = min(MASK_LL(32), newNumedges * 5LL / 4LL);

  int *newEnodeidMem = new int[newedgeMem];
  unsigned int *newedgeidMem = new unsigned int[newedgeMem];
  Cnode *newnode = new Cnode[newnumnodes];
  Cedge *newedge = new Cedge[newMaxedges];
  
  size_t EnodeidCnt = 0;
  size_t edgeidCnt = 0;

  map2nodeid = new int[numnodes];/* maps old node id (== mapid) to new node id */
  if(DEBUG>=2)
    for(int i = 0; i < numnodes; i++)
      map2nodeid[i] = -1;
  for(int i = 0; i < numtopnodes; i++)
    map2nodeid[topnodes[i]->nodeid] = i;

  unsigned int *mapedgeid = new unsigned int[Numedges+1];/* maps old edge id to new edge id */
  unsigned int edgeCnt = 0;
  for(unsigned int i = 1; i <= Numedges; i++){
    if(DEBUG>=2) mapedgeid[i] = 0;
    if(DEBUG>=2) assert(edge[i].edgeid == i);
    if(edge[i].valid < 0)
      continue;
    mapedgeid[i] = ++edgeCnt;
  }
  assert(edgeCnt == newNumedges);

  // Free up update and sitecov memory of node[0..numnodes-1] : should not be needed, but play it safe
  for(int i = 0; i < numnodes; i++){
    Cnode *pNode = &node[i];
    pNode->update_free();
  }

  /* copy and update node data */
  for(int i = 0; i < numtopnodes; i++){
    Cnode *pNode = &newnode[i];
    *pNode = *topnodes[i];// copy node data

    // reallocate Enodeid[] and edgeid[] as part of block memory EnodeidMem[], edgeidMem[]
    int *newEnodeid = &newEnodeidMem[EnodeidCnt]; EnodeidCnt += pNode->edgecnt;
    unsigned int *newedgeid = &newedgeidMem[edgeidCnt]; edgeidCnt += pNode->edgecnt;

    // update node data
    pNode->nodeid = i;
    pNode->maxedge = pNode->edgecnt;
    for(int k = 0; k < pNode->edgecnt; k++){
      int nodeid = pNode->Enodeid[k];
      if(DEBUG>=2) assert(0 <= map2nodeid[nodeid] && map2nodeid[nodeid] < numtopnodes);
      newEnodeid[k] = map2nodeid[nodeid];

      unsigned int edgeid = pNode->edgeid[k];
      if(DEBUG>=2) assert(1 <= edgeid && edgeid <= Numedges);
      if(DEBUG>=2) assert(1 <= mapedgeid[edgeid] && mapedgeid[edgeid] <= newNumedges);
      newedgeid[k] = mapedgeid[edgeid];
    }

    pNode->Enodeid = newEnodeid;
    pNode->edgeid = newedgeid;
    topnodes[i] = pNode;
  }

  /* free original node array and replace it with new array */
  delete [] node;
  node = newnode;

  // Free up update vectors for all edges : should not be needed, but play it safe 
  for(unsigned int i = 1; i <= Numedges; i++)
    edge[i].allfree();

  /* copy and update edge data */  
  edgeCnt = 0;
  for(unsigned int i = 1; i <= Numedges; i++){
    if(edge[i].valid < 0)
      continue;
    ++edgeCnt;

    /* copy edge data */
    newedge[edgeCnt] = edge[i];

    /* update edge data */
    Cedge *pedge = &newedge[edgeCnt];
    pedge->edgeid = edgeCnt;
    pedge->Lmap = map2nodeid[pedge->Lmap];
    if(DEBUG>=2) assert(0 <= pedge->Lmap && pedge->Lmap < newnumnodes);
    pedge->Rmap = map2nodeid[pedge->Rmap];
    if(DEBUG>=2) assert(0 <= pedge->Rmap && pedge->Rmap < newnumnodes);
  }

  /* free original edge array and replace it with new array */
  delete [] edge;
  edge = newedge;

  if(VERB){
    double etime = wtime();
    printf("Compacted graph : nodes= %d/%d -> %d, edges= %u/%u -> %u/%u, links= %lu -> %lu, Mem = %0.4f -> %0.4f Gb: time= %0.4f (wall= %0.4f)\n",
	   numnodes,maxnodes,newnumnodes,Numedges,Maxedges,newNumedges, newMaxedges,maxedgeMem, newedgeMem,
	   (maxnodes * sizeof(Cnode) + Maxedges * sizeof(Cedge) + maxedgeMem * 2ul * sizeof(int)) * 1e-9, 
	   (newnumnodes * sizeof(Cnode) + newMaxedges * sizeof(Cedge) + newedgeMem * 2ul * sizeof(int)) * 1e-9,
	   etime - stime, etime);
    fflush(stdout);
  }

  /* update node/edge counts */
  numnodes = maxnodes = newnumnodes;
  Numedges = newNumedges;
  Maxedges = newMaxedges;

  /* update EnodeidMem and edgeidMem */
  EnodeidMem = newEnodeidMem;
  edgeidMem = newedgeidMem;
  edgeMem = newedgeMem;

  delete [] mapedgeid;

  /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
  if(DEBUG>=2){
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edge[i].valid < 0)
	continue;
      assert(edge[i].mark == 0);
    }
    size_t edgecnt = 0;
    for(int i = 0; i < numnodes; i++){
      Cnode *Vnode = &node[i];
      if(Vnode->mark >= 2){
	assert(Vnode->edgecnt <= 0);
	continue;
      }
      edgecnt += Vnode->edgecnt;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	if(!(1 <= VWedgeid && VWedgeid <= Numedges)){
	  printf("i=%d/%d:node[i]:mapid=%d,nodeid=%d,mark=%d,edgecnt=%d/%d,vw=%d:node[i].edgeid[vw]=%u,VWedgeid=%u,Numedges=%u\n",
		 i,numnodes,Vnode->mapid,Vnode->nodeid,Vnode->mark,Vnode->edgecnt,Vnode->maxedge,vw,Vnode->edgeid[vw],VWedgeid,Numedges);
	  fflush(stdout);
	  assert(1 <= VWedgeid && VWedgeid <= Numedges);
	}
	assert(edge[VWedgeid].mark == 0);
      }
    }    
    if(DEBUG){
      for(int i = 0; i < numnodes; i++){
	Cnode *Vnode = &node[i];
	assert(Vnode->nodeid == i);
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];

	  VWedge->mark++;

	  if(VWedge->valid <= 0)
	    continue;
	  assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
	  assert(VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap));
	}
      }    
      for(int i = 0; i < numnodes; i++){
	Cnode *Vnode = &node[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->mark != 2){
	    printf("i=%d/%d: mapid=%d: vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d\n",
		   i,numnodes,Vnode->mapid, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark, VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped);
	    fflush(stdout);
	    assert(VWedge->mark == 2);
	  }
	}
      }    
      for(int i = 0; i < numnodes; i++){
	Cnode *Vnode = &node[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edge[Vnode->edgeid[vw]].mark = 0;
      }    
      /* check that topnodes[0..numtopnodes-1] reach all node[] that do not have mark >= 2 */
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	assert(Vnode->mark == 0);
	Vnode->mark = 2;
      }
      for(int i = 0; i < numnodes; i++)
	assert(node[i].mark >= 2);
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	Vnode->mark = 0;
      }
    }
    if((edgecnt % 2) != 0){
      printf("graph_compact():sum of node[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",edgecnt,numnodes,Numedges);
      fflush(stdout);
      assert((edgecnt%2) == 0);
    }
  }

  if(VERB){
    printf("Calling malloc_trim(0): wtime= %0.6f\n",wtime());
    fflush(stdout);
  }
  malloc_trim(0);
}

#if GRAPH_SMALL==0
#undef CnodeA
#undef CedgeA
#undef nodeA
#undef edgeA
#endif

// convert graph by expanding CnodeA -> Cnode, CedgeA -> Cedge, nodeA[] -> node[], edgeA[] -> edge[]
void graph_convert()
{
  if(!nodeA){
    assert(!edgeA);
    return;
  }

  if(VERB){
    printf("graph_convert(): expanding CnodeA -> Cnode, CedgeA -> Cedge\n");
    fflush(stdout);
  }

  edge = new Cedge[Maxedges];
  for(unsigned int i = 1; i <= Numedges; i++){
    Cedge *p = &edge[i];
    ((CedgeA *)p)[0] = edgeA[i];
  }
  delete [] edgeA;
  edgeA = NULL;

  node = new Cnode[maxnodes];
  for(int i = 0; i < numnodes; i++){
    if(DEBUG>=2) assert(0 <= nodeA[i].edgecnt && nodeA[i].edgecnt <= maxedges);

    Cnode *p = &node[i];
    ((CnodeA *)p)[0] = nodeA[i];

    if(DEBUG>=2) assert(0 <= node[i].edgecnt && node[i].edgecnt <= maxedges);
  }

  /* fix topnodes[] to point to node[] instead of nodeA[] */
  for(int i = 0; i < numtopnodes; i++){
    CnodeA *p = (CnodeA *)topnodes[i];
    int index = p - nodeA;
    if(DEBUG>=2) assert(0 <= index && index < numnodes);
    if(DEBUG>=2 && !(0 <= node[index].edgecnt && node[index].edgecnt <= maxedges)){
      Cnode *q = &node[index];
      printf("i=%d:topnodes[i]->mapid=%d,nodeid=%d,edgecnt=%d/%d, topnodes[i]-nodeA= index= %d:\n",
	     i,p->mapid,p->nodeid,p->edgecnt,p->maxedge,index);
      fflush(stdout);
      printf("\t node[index].mapid=%d,nodeid=%d,edgecnt=%d/%d, maxedges=%d\n",
	     q->mapid,q->nodeid,q->edgecnt,q->maxedge,maxedges);
      fflush(stdout);
      assert(0 <= node[index].edgecnt && node[index].edgecnt <= maxedges);
    }

    topnodes[i] = &node[index];
  }

  for(int i = 0; i < numnodes; i++){
    nodeA[i].edgeid = NULL;// so ~CnodeA() will not delete[] this array below
    nodeA[i].Enodeid = NULL;// so ~CnodeA() will not delete[] this array below
  }
  delete [] nodeA;
  nodeA = NULL;

  /* make sure each edge is connected to only 2 nodes : increment edge mark from each node that can reach it (should add up to 2 nodes/edge) */
  if(DEBUG>=2){
    for(unsigned int i = 1; i <= Numedges; i++){
      if(edge[i].valid < 0)
	continue;
      assert(edge[i].mark == 0);
    }
    size_t edgecnt = 0;
    for(int i = 0; i < numnodes; i++){
      Cnode *Vnode = &node[i];
      if(Vnode->mark >= 2){
	assert(Vnode->edgecnt <= 0);
	continue;
      }
      edgecnt += Vnode->edgecnt;
      for(int vw = 0; vw < Vnode->edgecnt; vw++){
	unsigned int VWedgeid = Vnode->edgeid[vw];
	if(!(1 <= VWedgeid && VWedgeid <= Numedges)){
	  printf("i=%d/%d:node[i]:mapid=%d,nodeid=%d,mark=%d,edgecnt=%d/%d,vw=%d:node[i].edgeid[vw]=%u,VWedgeid=%u,Numedges=%u\n",
		 i,numnodes,Vnode->mapid,Vnode->nodeid,Vnode->mark,Vnode->edgecnt,Vnode->maxedge,vw,Vnode->edgeid[vw],VWedgeid,Numedges);
	  fflush(stdout);
	  assert(1 <= VWedgeid && VWedgeid <= Numedges);
	}
	assert(edge[VWedgeid].mark == 0);
      }
    }    
    if(DEBUG){
      for(int i = 0; i < numnodes; i++){
	Cnode *Vnode = &node[i];
	assert(Vnode->nodeid == i);
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  VWedge->mark++;
	  if(VWedge->valid <= 0)
	    continue;
	  assert(VWedge->Lmap == Vnode->nodeid || VWedge->Rmap == Vnode->nodeid);
	  assert(VWedge->Lmap == Vnode->nodeid ? (Vnode->Enodeid[vw] == VWedge->Rmap) : (Vnode->Enodeid[vw] == VWedge->Lmap));
	}
      }    
      for(int i = 0; i < numnodes; i++){
	Cnode *Vnode = &node[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++){
	  Cedge *VWedge = &edge[Vnode->edgeid[vw]];
	  if(VWedge->mark != 2){
	    printf("i=%d/%d: mapid=%d: vw= %d/%d: edgeid=%d, mark=%d, Lmap=%d,Rmap=%d,Lflipped=%d,Rflipped=%d\n",
		   i,numnodes,Vnode->mapid, vw, Vnode->edgecnt, VWedge->edgeid, VWedge->mark, VWedge->Lmap, VWedge->Rmap, VWedge->Lflipped,VWedge->Rflipped);
	    fflush(stdout);
	    assert(VWedge->mark == 2);
	  }
	}
      }    
      for(int i = 0; i < numnodes; i++){
	Cnode *Vnode = &node[i];
	for(int vw = 0; vw < Vnode->edgecnt; vw++)
	  edge[Vnode->edgeid[vw]].mark = 0;
      }    
      /* check that topnodes[0..numtopnodes-1] reach all node[] that do not have mark >= 2 */
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	assert(Vnode->mark == 0);
	Vnode->mark = 2;
      }
      for(int i = 0; i < numnodes; i++)
	assert(node[i].mark >= 2);
      for(int i = 0; i < numtopnodes; i++){
	Cnode *Vnode = topnodes[i];
	Vnode->mark = 0;
      }
    }
    if((edgecnt % 2) != 0){
      printf("graph_convert():sum of node[0..numnodes-1].edgecnt = %lu : should be multiple of 2 (numnodes=%d,Numedges=%u)\n",edgecnt,numnodes,Numedges);
      fflush(stdout);
      assert((edgecnt%2) == 0);
    }
  }

  if(VERB){
    printf("Calling malloc_trim(0): wtime= %0.6f\n",wtime());
    fflush(stdout);
  }
  malloc_trim(0);
}
