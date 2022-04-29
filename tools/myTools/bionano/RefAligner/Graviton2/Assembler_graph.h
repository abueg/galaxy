# pragma once /*http://connect.microsoft.com/VisualStudio/feedback/details/651405/pch-warning-after-installing-sp1*/
#ifndef ASSEMBLER_GRAPH_H
#define ASSEMBLER_GRAPH_H

#include "Cnode.h"

/** \file Assembler_graph.h \brief Prototypes for Assembler_graph.cpp. 
*
* Mostly functions, minimal class structure.
*/

static Ident Assembler_graph_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Assembler_graph.h 7926 2018-09-27 01:35:42Z tanantharaman $");

extern int numnodes;
extern Cnode *node;

extern int numedges;
extern Cedge *edge;

#if CALIGN_SMALL
extern void graph_build(Cmap **map, int nummaps, CalignA **alignmentA, size_t &numaligns);/* build graph from maps and alignments */
#else
extern void graph_build(Cmap **map, int nummaps, Calign **alignment, size_t &numaligns);/* build graph from maps and alignments */
#endif

extern void graph_display(int mincov, int verb, char *suffix);/* display graph showing only edges with coverage >= mincov */

extern void graph_edgecnt(int verb);

extern int graph_edgetrim(int verb);

extern long long graph_nodetrim(int verb);

extern void graph_convert();

extern long long graph_chimtrim(int verb);

extern void graph_redundancy(int verb);/* Algorithm 4.1.1 in Nguyen's thesis */

extern void graph_compact();

extern int graph_collapse(int mincov, int verb);/* Chain Collapsal in Nguyen's thesis */

extern int graph_bulge(int BFSdepth, int verb);/* Graph Bulges in Nguyen's thesis */

extern int graph_prune(int side, int verb);/* prune edges with coverage < MinCoverage
	      If side>0: If single node connects to middle of graph, move it to side chain */
extern void graph_components();/* find connected components using DFS */

extern void graph_free();

#endif
