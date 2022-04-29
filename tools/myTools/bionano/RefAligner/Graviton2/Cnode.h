#ifndef CNODE_H
#define CNODE_H

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "globals.h"

#if CALIGN_SMALL
#include "CalignA.h"
#endif

#include "Calign.h"


extern int *EnodeidMem;
extern unsigned int *edgeidMem;

class Cedge;
class Cnode;

// static long long Cupdate_cnt = 0, Cupdate_hwm = 0;

/** Class to store significantly confirmed edges w detail.  */
class Ccrossbreak {
public:
  int connected;
  Cedge *directEdge;

  Cedge *backEdge;
  Cnode *backNode;

  Cedge *foreEdge;
  Cnode *foreNode;
};

/** Graph Edge represents the overlap between a map pair */
class CedgeA {
 public:

#if DEBUG>=2
  size_t size,sized;
#endif

#if DAG_CHECK || DEBUG>=2
  int updatemark;
#endif

  unsigned int edgeid;/**< index into edge[1..NumEdges] */

  int Lmap;/**< index into node[] of leftmost map (relative to map centers and with Lflipped/Rflipped) */
  int Rmap;/**< index into node[] of rightmost map (relative to map centers and with Lflipped/Rflipped) */

  int Lflipped;/**< If left map is flipped */
  int Rflipped;/**< If right map is flipped */
  /* Note : Lflipped and Rflipped are never both 1 : whenever this happens, both maps are flipped and Lmap,Rmap exchanged. */

  double distance;/**< Distance between Lmap and Rmap (>= 0.0) : initially the distance between map centers (= distsum/wt) */
  double cov;/**< coverage(redudancy) count for this edge (initially 1 (OR log10(PV)/log10(PVthreshold) if REDUNDANCY_WEIGHTED)) */

  long long alignid;/* alignment[alignid-1] is pointer to overlap alignment (0 for Super-Edges, -1 for Uninitialized). For deallocated edges valid == -1 */

  int mark; /**< initialized to 0 */

  int valid;/**< 1 : Normal Edge
	       0 : Edge scheduled for removal after propagation of coverage updates are complete.
	       -1 : Edge marked for removal from Graph (It is treated as removed, until it is actually removed by GCedge) */

  CedgeA(){
#if DAG_CHECK
    updatemark = -1;
#endif
    alignid = -1;
    mark = 0;
  };

 private:
};

class Cedge : public CedgeA {
 public:
  /* NOTE all fields after this are only used after graph_nodetrim() */
  
  double wt;/**< weight(quality) of distance for all merged edges in this edge : initially the number of aligned site pairs */
  double distsum;/**< sum of distance * wt for all merged edges in this edge */

  int Lcnt;/**< number shared left neighbors (computed in graph_edgecnt()) */
  int Rcnt;/**< number of shared right neighbors (computed in graph_edgecnt()) */
  int Bcnt;/**< number of edges bridging this edge (computed in graph_edgecnt()) */
  int Ccnt;/**< min(Lcnt,Rcnt) (computed in graph_edgecnt()) */
  int ucnt;/**< number of other paths that will add cov/ucnt to their edge/node's cov value 
	      (If ucnt > 0 then this edge will be eventually deleted and valid==0) */

  std::vector<int> update;/**< vector of edge ids whose final "cov/ucnt" value will be added to this edge */
  std::vector<int> updated;/**< IF CROSSBREAK : vector of edge ids whose final "cov/ucnt" value was already added to this edge */
  double updistance;/**< maximum distance value of any ancestor edge across *update links (or current distance, whichever
		       is larger). Whenever a new update link is added, the current edge's distance value is propagated
		       forward (across update links) as long as the next edge'e updistance value is not already as large.
		       Conversely, there is no need to check for cycles if the new update edge is going to an edge
		       with larger updistance value. Most of the time update links point from smaller to larger updistance
		       values, so this should be less time consuming than checking for cycles each time a new update edge
		       is added.
		    */

  unsigned int Lchainid, Rchainid;
  //  Cedge *Lchain,*Rchain;
  /**< If valid >= 0 && alignid == 0, this edge is a Super-Edge that represets a chain of nodes
		       starting at edge Lchain and ending at edge Rchain. 

		       If Lmap and Rmap are oriented the same way in the chain, this edge is oriented
		       to make Lflipped == Rflipped == 0.

		       Otherwise this Super-Edge is oriented arbitrarily.
		       
		       Note that(Lchain->Lmap == Lmap || Lchain->Rmap == Lmap) && 
		                (Rchain->Lmap == Rmap || Rchain->Rmap == Rmap)

		       This Super-Edge can be expanded by :

		       1. Adding Edge Lchain to Lmap (Lchain is still connected to the other end)
		            (in place of this edge).
		       2. Adding Edge Rchain to Rmap (Rchain is still connected to the other end)
		            (in place of this edge).
		       3. Add the nodes in the chain back to the topnodes list.
		       
		       Note that edges in the chain (including Lchain,Rchain) can themselves be Super-Edges,
		       since multiple cycles of chain-collapsel and bulge removal take place.
		    */

  Cedge(){
    ucnt = 0;
    Ccnt = Lcnt = Rcnt = Bcnt = -1;
  };
  ~Cedge(){
    allfree();
  };
  void update_free(){
    std::vector<int>().swap(update);
    std::vector<int>().swap(updated);
  }
  void allfree(){
    update_free();
  };
 private:
};


/** class type to save result of a ShortestPath() call */
/** (node,flip) pair (oriented node) type used for Shortest Path algorithm */
class Cnodeflip {
 public:
  int nodeid;
  int flip;
};

class CSPresult {
 public:
  int SourceNode;/**< source node as index into node[] */
  int dir;/**< orientation of source node */
  Cnodeflip EndNode;/**< most distant node from (SourceNode,dir) */
  double pairdist;/**< distance of EndNode from (SourceNode,dir) */
  double paircost;/**< heuristic cost of least cost path to EndNode from (SourceNode,dir) */
};

/** class type for ShortestPath() per node data (each thread has its own copy with multi-threading) */
class CSPdata {
 public:
  int SPmark;/**< used for DFS/BFS */
  int hindex;/* used by BFS (used as edgecnt by DFS) */
  int *Enodeid;/* copy of node[nodeid].Enodeid[0..edgecnt-1] */
  double hdist;/* used by BFS */
  int heapindex[2];/**< heapindex[flip] is the index into heap used by Dijkstra's algorithm for the node in orientation = flip */
  double SPcost[2][2];/**< SPcost[dir][flip] is the lowest cost distance from root node to this node in orientation=flip using forward (backward IFF dir=1) links only (-1 if not reachable) */
  double SPdist[2][2];/**< SPdist[dir][flip] is the actual distance of the lowest cost path (see SPcost[][]) (If cost == distance, then SPdist[][] == SPcost[][]) */
  Cedge *SPedge[2][2];/**< SPedge[dir][flip] is the edge to next node on shortest path from root to this node in orientation=flip using forward (backward IFF dir=1) links only */
};

/** class type for Breadth First Search per node data (each thread has its own copy with multi-threading) */
class BFSdata {
 public:
  int mark;
  int SPmark;
  Cnode *parent;/**< previous node on path from root */
  Cedge *pedge;/**< Edge between parent and current node */
  int pflip;/**< orientation of current node on path from root (root may be flipped as well) */
  double Rdist;/**< net distance from root node (via parent node) */
  double Rdwt;/**< sum of distance/wt from root node (via parent node) */
  double mincov;/**< minimum coverage of all edges on path from root node */
};

extern void maxedgealloc(unsigned int num, unsigned int &maxedges, Cedge* &edge);
extern void maxedgealloc(unsigned int num, unsigned int &maxedges, CedgeA* &edgeA);

/** Graph Node represents a map and a set of Edges (overlaps with other maps) */
/** Each Edge has two pointers : one from each Node at either end of the Edge */
class CnodeA {
 public:
  int nodeid;/**< index into node[] */
  int mapid; /**< index into map[] : two nodes may have the same mapid if the map was chimeric (see CROSSBREAK, once node cloning is implemented) */

  int edgecnt;/**< number of Edges */
  int *Enodeid;/**< If != 0 : Enodeid[i= 0..edgecnt-1] : nodeid of neighboring node : short for ((edge[edgeid[i]]->Lmap == nodeid  ? edge[edgeid[i]]->Rmap : edge[edgeid[i]]->Lmap) */
  unsigned int *edgeid;/**< If != 0 : edgeid[i= 0..edgecnt-1] : edgeid (index into edge[]) of link to neighboring nodes */

  int PosEdge;/**< for use by graph algorithms : first edge with positive distance from this node (unflipped) */

  int maxedge;/**< Number of Edges allocated for edgeid[] and Enodeid[] */
  
  /* NOTE : all fields below only needed after graph_build() */

  int mark;/** for use by graph algorithms : 0 = Unmarked, 1 = Marked, 2 = Eliminated (not in topnodes[]) */

  //  int *sitecov;/** coverage profile relative to all other maps that aligned with it (can be used to clean up map data) */
  int trimL,trimR;/** trimmed sites : left end sites 0..trimL and right end sites trimR .. M+1 (M == map[mapid]->numsite[0])
		     are treated as trimmed (NOTE : trimL==0 and trimR==M+1 nothing was trimmed) */


  double len;/**< map length in kb */
  double cov;/**< coverage (redudancy) count for this node (initially 1) : Only computed if valid == 1*/
  double forward,backward;/** total coverage of forward and backward edges (computed when graph was first constructed and in PosEdge()) */

  CnodeA() {
    maxedge = 0;
    edgeid = NULL;
    Enodeid = NULL;
    mark = 0;
    trimL = trimR = 0;
  };
  
  ~CnodeA() {
    allfreeA();
  }

  void allfreeA(){
    if(edgeid){
      if(!edgeidMem) delete [] edgeid;
      edgeid = 0;
    }
    if(Enodeid){
      if(!EnodeidMem) delete [] Enodeid;
      Enodeid = 0;
    }
    maxedge = 0;
  };

  void maxedgealloc(int numedges);/**< allocate/expand Edge pointer array : new pointers are initialized to 0 
				     (no actual Cedge allocated) */
 private:
};

class Cnode : public CnodeA {
 public:

  /* NOTE : all field below only needed after graph_nodetrim() */
  int Lcntmax,Rcntmax;/** largest value for Ccnt for left/right edges respectively : see graph_edgecnt() */

  std::vector<int> update;/**< In graph_redundancy : list of edges whose final "cov" value will be added to this Node  */
  std::vector<int> updated;/**< list of edges already added to this Node (used only if CROSSBREAK) */
  std::vector<int> updateL;
  std::vector<int> updateR;/**< to track coverage from edges to the right : If seperate coverage from "left" and "right" is tracked */
             
  std::vector<int> chain; /**<  list of side chain edges, to be used in final consensus map construction */

  /* used in BFS (along with mark) */
  Cnode *parent;/**< previous node on path from root */
  Cedge *pedge;/**< Edge between parent and current node */
  int pflip;/**< orientation of current node on path from root (root may be flipped as well) */
  double Rdist;/**< net distance from root node (via parent node) */
  double Rdwt;/**< sum of distance/wt from root node (via parent node) */
  double mincov;/**< minimum coverage of all edges on path from root node */

  Cedge *bestDFS;/**< used by DFSforward() */
  int flipDFS;/**< used by DFSforward() */

  /** used by graph_redudundancy if CROSSBREAK */
  int bcnt,fcnt;
  Ccrossbreak **pairs;/**< pairs[0..bcnt-1][0..fcnt-1] for pairs of edges with significant join confirmation */
  Ccrossbreak *pairsmem;

  /** fields used by ShortestPaths() */
  int SPmark;/**< used for DFS */
  //  int heapindex[2];/* heapindex[flip] is the index into heap used by Dijkstra's algorithm for the node in orientation = flip */
  //  double SPcost[2][2];/* SPcost[dir][flip] is the lowest cost distance from root node to this node in orientation=flip using forward (backward IFF dir=1) links only (-1 if not reachable) */
  //  double SPdist[2][2];/* SPdist[dir][flip] is the actual distance of the lowest cost path (see SPcost[][]) (If cost == distance, then SPdist[][] == SPcost[][]) */
  //  Cedge *SPedge[2][2];/* SPedge[dir][flip] is the edge to next node on shortest path from root to this node in orientation=flip using forward (backward IFF dir=1) links only */

  Cnode() {
    SPmark = 0;
    parent = 0;
    pairs = 0;
    pairsmem = 0;
  };
  
  ~Cnode() {
    update_free();
    std::vector<int>().swap(chain);
  }

  void update_free(){
    std::vector<int>().swap(update);
    std::vector<int>().swap(updated);
    std::vector<int>().swap(updateL);
    std::vector<int>().swap(updateR);
    if(pairs){
      delete [] pairsmem;
      delete [] pairs;
      pairs = 0;
      pairsmem = 0;
    }
  };

  void allfree(){
    allfreeA();
    update_free();
    std::vector<int>().swap(chain);
  };

 private:
};

using namespace std;

/** data type to hold the end of a path from a shared source(root) node */
class Cpath {
 public:
  Cnode *node;/**< end of path */
  int flipped;/**< if map node->mapid is flipped (relative to shared source node not being flipped) */
  double priority;/**< smaller values have higher priority : priority corresponds to distance
		   from source node (either number of hops or actual edge distance) */
  Cpath(Cnode *N, int flip, double prio){
    node = N; flipped = flip; priority = prio;
  };
  bool operator<(const Cpath &a) const {
    return priority > a.priority;/* this is deliberately backwards so we can use less<Cpath> as the comparison function !!! */
  };
};

/** Path through graph defined by array of nodeids (N), alignments (N-1), and orientations (N). */
class Ccontigpath {
 public:
  int nummaps;
  Cnode *start;/**< linked list of nodes: each node Vnode is linked to Vnode->parent (orientation Vnode->pflip)
		  via edge Vnode->pedge and alignment Vnode->pedge->align */
  int flipstart;/* orientation of first map */

  int *nodeid;/**< array of node ids nodeid[i=0..nummaps-1] : corresponding Cmap is map[node[nodeid[i]]->mapid] */
#if CALIGN_SMALL
  CalignA **align;/**< align[i = 0..nummaps-2] are alignments between maps of nodeid[i] and nodeid[i+1]
 */
#else
  Calign **align;/**< align[i = 0..nummaps-2] are alignments between maps of nodeid[i] and nodeid[i+1] */
#endif
  int *flip;/**< array of orientation of maps */

  Ccontigpath(){
    nodeid = 0;
    align = 0;
    flip = 0;
  };
  ~Ccontigpath(){
    if(nodeid)
      delete [] nodeid;
    if(align)
      delete [] align;
    if(flip)
      delete [] flip;
  }
};

class CnodedistanceA {
public:
  CedgeA *pedge;
  double nodedistance;
};

class Cnodedistance {
public:
  Cedge *pedge;
  double nodedistance;
};

#include "ident.h" // "globals.h"
static Ident Cnode_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Cnode.h 10473 2020-01-08 08:18:14Z tanantharaman $");

#endif

