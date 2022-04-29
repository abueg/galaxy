#include "Assembler_parameters.h"

static Ident Assembler_parameters_cpp_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Assembler_parameters.cpp 7414 2018-02-26 05:14:39Z tanantharaman $");

#include "parameters.cpp" /* moved shared parameters to parameters.cpp */

/* parameters used only in Assembler */

char *align_filename[MAXFILES];/* input .align file name(s) */
int num_align_files = 0;/* number of input .align files */
char *refmap_filename = 0;/* input .map file name (if available for debugging) */
char *contig_filename = 0;/* input .contig file name (if available) */
int SAVE_MAPSET= 0; /* save mapset used for each refined contig <N> in <PREFIX>_contig<N>.bnx */

int contig_start = 0,contig_end = 0;/* contig numbers to refine (if -contigs is specified) */
int contig_format = 0;/* 0 : old text format
			 1 : new compact format (mostly binary, except for headers) */

double Pchim = 0.001;/* chimeric alignments pvalue : if low score in middle of Pchim_minlen or more segments exceeds -log(pvalue) the
			alignment is considered chimeric and ignored */
int Pchim_minlen = 3;

int ENDSITES = 6; /* If more than this many sites at either end are unaligned (on both molecules) the alignment 
		     is considered chimeric and ignored/deleted */

double SM = 2.0;/* kb sizing error in pairwise overlap offset (midpoint distance) */

double EdgeTrimRatio = 0.20;
int MinCcnt = 0;
int MaxCcnt = 0/* TRY 1,2 */;

double TrimCoverage = 2.99;
double TrimCoverageRatio = 0.20;
double ChimCoverageRatio = 0.20;
int MinCcntmax = 2;

int MinCoverage = 2;
int MaxCoverage = 64 ;

double ChimRatio = 0.20;
double EdgeCoverageRatio = 0.20;
double MinCoverageRatio = 0.20;
double CycleCoverageRatio = 0.20;/* Minimum ratio of edge coverage for weakest edge vs average edge in cycle (otherwise weakest edge is deleted) */

double BulgeBreakMaxCov = 3.9; /* maximum coverage link to break when 2 paths do not match */
double CycleBreakMaxCov = 19.9;/* maximum coverage link to break in cycles with at least one edge >  MAXERR*SM */

double ChimTrimRatio = 0.15 /* TRY 0.20 : may trim away fragile regions unless MinChimCcnt is increased */;
int MinChimCcnt = 6; /* Using an absolute threshold is ugly but seems to work to avoid trimming regions near fragile sites */

double CrossbreakMincov = 7.5;/* ignore edges with coverage less than this amount during CROSSBREAK analysis */
int CrossbreakMincnt = 5;/* confirm connections with min(Lcnt,Rcnt) >= CrossbreakMincnt during CROSSBREAK analysis (New method) */
double CrossbreakMinratio = 0.20;/* minimum independent coverage confirmation for pairs of edges connecting through a node (Old method) */

double MaxBulge = 2000.0;/* The effective value is increased gradually from 100.0 to this value (in kb) */
double MinSidebranch = 0.0;/* If > 0.0 : Output side branches larger than this size (if branch has more than ContigMinMaps maps) */
double MaxSidebranch = 60.0; /* hide side branches up to this size in side chains (typically set to a large value if MinSidebranch > 0.0) */
int SideChain = 0;/* If != 0 : Enable hiding side branches with single super-edge to help simplify graph */
int SidebranchExt = 0;/* If != 0 : Extend side branches and alternate main paths in both directions up to the next branch point */

int ContigMinMaps = 50;/* minimum maps in output contigs */
double ContigMinCov = 4.0;/* minimum average depth of output contigs */
double ContigMinLen = 50.0;/* minimum contig length (in excess of average map length) */

int MaxRelCoverage = 100;/* maximum input edgecnt as multiple of average edgecnt : if higher, keep only edges with highest logPV */
int MaxAbsCoverage = 200;/* If > 0 : maximum input edgecnt per map : if higher, keep only only edges with higher logPV rank (on either node) */
int MaxAbsCoverage2 = 0;/* If > 0 : keep only edges with higher logPV rank (on both nodes) */
int MaxCoverageDump = 0;/* If > 1 : output alignments remaining after applying -MaxRelCoverage, -TrimCoverage, -MinCcntmax exit (useful for computing more accurate Hashtable ROC curve OR saving memory) */
int SkipFilter = 0;/* If > 1 : assume input is a result of -dumpalign and skip -MaxRelCoverage, -TrimCoverage, -MinCcntmax filters */

int BulgeFast = 8;/* Speed up graph_bulge() by skipping Dykstra shortest path check if BFSdepth <= BulgeFast */
int FastPaths = 0;/* Speed up final graph component analysis : 1 = conservative, 2 = agressive */

int SizeType = 0;/* How to average sizes during draft contig computation :
		    0 : Assumes var(y) = SD[0]*SD[0]*y
		    1 : Assumes var(y) = SD[0]*SF[0]  (And downweights size mismatches based on -Poutlier) */
int FragilePreserve = 1; /* Try to preserve fragile sites during graph_chimtrim() : Will make it non-deterministic with multithreading */

double GraphWeight = 0.0; /* If > 0.0 : weight input graph edges (alignments) by Weight = (pow(GraphWeight,log10(T)+logPV)-1)/(GraphWeight-1.0)
			     NOTE : If GraphWeight = 1.0 + Eps : Weight ~ log10(T)+logPV as Eps -> 0
			            If GraphWeight = Eps : Weight ~ 1 as Eps -> 0 (The default behaviour with GraphWeight == 0.0) */
double GraphMaxWeight = 0.0;/* If > 0.0 : limit Weight of graph edges (when adjusted due to GraphWeight) to no more than this value */

int AlignmentFilter = 0;/* If > 0 : If number of alignments per map exceeds AlignmentFilter, increment MinLen by minlenDelta and LogPvThreshold by LogPvDelta */
double minlenDelta = 0.0;
double LogPvDelta = 0.0;
double AlignmentsPerGb = 1000000.0;
