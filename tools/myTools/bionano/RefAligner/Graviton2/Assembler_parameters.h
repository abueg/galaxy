#ifndef ASSEMBLER_PARAMETERS_H
#define ASSEMBLER_PARAMETERS_H

#include "globals.h"

static Ident Assembler_parameters_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Assembler_parameters.h 5405 2016-10-02 03:19:46Z tanantharaman $");

#define AVOID_BNX 1 /* With -contigs  avoid reading in maps until after .contigs file and read in only maps that are needed by .contigs file */

#include "parameters.h"

/* parameters used only in Assembler */
extern char *align_filename[MAXFILES];/* input .align file name(s) */
extern int num_align_files;/* number of input .align files */
extern char *refmap_filename;/* input .map file name (if available) */
extern char *contig_filename;/* input .contig file name (if available) */
extern int SAVE_MAPSET; /* save mapset used for each refined contig <N> in <PREFIX>_contig<N>.bnx */

extern int contig_start,contig_end;/* contig index range to refine (if -contigs is specified) */
extern int contig_format;/* 0 : old text format
			    1 : new compact (mostly binary, except for headers) format */
#define CONTIG_MAGIC "End of Binary Chunk" /* Magic string to seperate binary contig chunks : it is assumed this byte sequence followed by "\n" never occurs in a binary chunk */

extern double Pchim;/* chimeric alignments (low score in the middle with at least Pchim_minlen segments) pvalue */
extern int Pchim_minlen;/* minimum number of segments for Pchim */

extern int ENDSITES; /* filter out chimeric alignments with misaligned ends with at least this many misaligned sites on X or Y */

extern double SM;/* fixed error sigma for map distances */

#define NODETRIM_FIRST 1 /* call graph_nodetrim() before graph_edgetrim() */

extern double EdgeTrimRatio;/* minimum ratio of edge min(Lcnt,Rcnt)==Cnt to node Ccntmax (otherwise edge is trimmed away as spurious provided Ccnt > 0) */
extern int MinCcnt;/* miminum edge value of min(Lcnt,Rcnt)==Ccnt (otherwise edge is trimmed away as spurious, 
		      provided that also Ccnt < min(forward,backward)*EdgeTrimRatio) */
extern int MaxCcnt;/* minimum edge value of max(Lcnt,Rcnt) (otherwise edge is trimmed away as spurious) : should not but somes can trim away edges near fragile nicking sites */

extern double TrimCoverage;/* should be a fraction of average over all nodes */ /* minimum average coverage : nodes below this are trimmed away */
extern double TrimCoverageRatio;/* Minimum pairwise site coverage relative to the mean value : parts of maps below this are trimmed away */
extern double ChimCoverageRatio;/* minimum site interval coverage relative to highest coverage each side of an internal region : If internal coverage
				   drops below this threshold, node is deleted as chimeric */
extern int MinCcntmax;/* If largest Ccnt for any edge at a node is below this value, delete the node as chimeric */

extern int MinCoverage;/* Minimum coverage for assembly path (edges and nodes) */
extern int MaxCoverage;/* Treat all coverage values above this as equivalent (Used by graph_collapse())*/

extern double ChimRatio;	/* Ratio of node coverage to forward or backward edge coverage for Chimeric nodes
				   (nodes with ratios below this are assumed Chimeric maps) */
extern double EdgeCoverageRatio; /* Minimum ratio of coverage for 3rd edge vs 2nd edge (sorted by coverage) 
				    for 3rd edge to be retained */
extern double MinCoverageRatio; /* Minimum ratio of edge coverage for alternate path vs top path 
				   (otherwise alternate path is deleted) */
extern double CycleCoverageRatio;/* Minimum ratio of edge coverage for weakest edge vs average edge in cycle (otherwise weakest edge is deleted) */

extern double BulgeBreakMaxCov; /* try to increase to avoid absolute threshold */  /* maximum coverage link to break when 2 paths do not match 
										     and the lower coverage link < average coverage * MinCoverageRatio */
extern double CycleBreakMaxCov;/* maximum coverage link to break in cycles with at least one edge >  MAXERR*SM */

extern double ChimTrimRatio; /* Minimum Ratio of crossedges*0.5/sqrt(set1edges*set2edges) : otherwise smaller set's edges are removed */
extern int MinChimCcnt; /* Mimimum Ccnt average value for edges that may be deleted in graph_chimtrim() from nodes that appear to join 
			   two different regions of the genome .  This avoids disconnecting connections crossing a weak nick site. */

extern double CrossbreakMincov;/* ignore edges with coverage less than this amount during CROSSBREAK analysis */
extern int CrossbreakMincnt;/* confirm connections with min(Lcnt,Rcnt) >= CrossbreakMincnt during CROSSBREAK analysis (New method) */
extern double CrossbreakMinratio;/* minimum independent coverage confirmation for pairs of edges connecting through a node (Old method) */

extern double MaxBulge;/* maximum length (in kb) of a bulge or cycle (larger values could correspond to actual genome circular chromosome length) */
extern double MinSidebranch;/* If > 0.0 : Output side branches larger than this size (and having more than ContigMinMaps maps) */
extern double MaxSidebranch; /* hide side branches up to this size in side chains */
extern int SideChain;/* If != 0 : Enable hiding side branches with single super-edge to help simplify graph */
extern int SidebranchExt;/* If != 0 : Extend side branches and alternate main paths in both directions up to the next branch point */

extern int ContigMinMaps;/* minimum number of maps in any Contig (path) to output */
extern double ContigMinCov;/* minimum average "cov" value for any Contig (path) to output */
extern double ContigMinLen;/* minimum contig length (in excess of average map length) */

extern int MaxRelCoverage;/* maximum input edgecnt as multiple of average edgecnt : if higher, keep only edges with highest logPV */
extern int MaxAbsCoverage;/* If > 0 : maximum input edgecnt per map : if higher, keep only only edges with highest logPV */
extern int MaxAbsCoverage2;/* If > 0 : keep only edges with higher logPV rank (on both nodes) */
extern int MaxCoverageDump;/* If > 1 : output alignments remaining after applying -MaxRelCoverage & exit (useful for computing more accurate Hashtable ROC curve) */
extern int SkipFilter;/* If > 1 : assume input is a result of -dumpalign and skip -MaxRelCoverage, -TrimCoverage, -MinCcntmax filters */

extern int BulgeFast;/* Speed up graph_bulge() by skipping Sykstra shortest path check if BFSdepth <= BulgeFast */
extern int FastPaths;/* Speed up final graph component analysis : 1 = conservative, 2 = agressive */

extern int SizeType;/* 0 : Assumes var(y) = SD[0]*SD[0]*y
		       1 : Assumes var(y) = SF[0]*SF[0] */

extern int FragilePreserve; /* Try to preserve fragile sites during graph_chimtrim() : Will make it non-deterministic with multithreading */

extern double GraphWeight; /* If > 0.0 : weight input graph edges (alignments) by Weight = (pow(GraphWeight,log10(T)+logPV)-1)/(GraphWeight-1.0)
			     NOTE : If GraphWeight = 1.0 + Eps : Weight ~ log10(T)+logPV as Eps -> 0
			            If GraphWeight = Eps : Weight ~ 1 as Eps -> 0 (The default behaviour with GraphWeight == 0.0) */
extern double GraphMaxWeight;/* If > 0.0 : limit Weight of graph edges (when adjusted due to GraphWeight) to no more than this value */

extern int AlignmentFilter;/* If > 0 : If number of alignments per map exceeds AlignmentFilter, increment MinLen by minlenDelta and LogPvThreshold by LogPvDelta */
extern double minlenDelta;
extern double LogPvDelta;
extern double AlignmentsPerGb;

#endif
