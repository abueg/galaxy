#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "constants.h"
#include "Cmap.h"
#include "Cnode.h"
#include "Ccontig.h"

/** \file Assembler.h \brief Simple Prototyper for Assembler.cpp. 
*
*/


static Ident Assembler_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Assembler.h 2946 2014-07-15 00:13:26Z tanantharaman $");

extern void output_TrimmedMaps(int nummaps, Cmap **map, Cnode *node, char *filename);


#endif
