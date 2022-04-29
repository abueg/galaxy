
#include <stdlib.h>
#include <stdio.h>

//#undef DEBUG
//#define DEBUG 2

#include "globals.h"
#include "Cnode.h"


/** \file Cnode.cpp increase size of array edge[0..maxedge-1] to at least edge[0..num-1] */

static Ident Id("$Header: $");


void CnodeA::maxedgealloc(int num)
{
  int newedges = num;
  if(newedges <= maxedge)
    return;
  if(newedges < (maxedge+1024)*5/4)
    newedges = (maxedge+1024)*5/4;
  if(maxedge <= 0){
    edgeid = new unsigned int[newedges];
    if(DEBUG>=2)
      for(int i=0;i < newedges;i++)
	edgeid[i] = 0;
    Enodeid = new int[newedges];
  } else {
    unsigned int *origedgeid = edgeid;
    int *origEnodeid = Enodeid;
    edgeid = new unsigned int[newedges];
    Enodeid = new int[newedges];
    for(register int i=0; i < maxedge; i++){
      edgeid[i] = origedgeid[i];
      Enodeid[i] = origEnodeid[i];
    }
    if(DEBUG>=2)
      for(int i=maxedge; i < newedges;i++)
	edgeid[i] = 0;
    delete [] origedgeid;
    delete [] origEnodeid;
  }
  maxedge = newedges;
}

void maxedgealloc(unsigned int num, unsigned int &maxedges, Cedge* &edge)
{
  unsigned int newedges = num + 1;
  if(newedges <= maxedges)
    return;

  if(newedges < 1024)
    newedges = 1024;
  if(maxedges > 0 && newedges < min(MASK_LL(32), (maxedges+1024)*5LL/4))
    newedges = min(MASK_LL(32), (maxedges+1024)*5LL/4);

  if(VERB/* HERE HERE >=2 */){
    printf("maxedgealloc: num=%u,maxedges=%u: newedges= %u\n",num,maxedges,newedges);
    fflush(stdout);
  }

  if(maxedges <= 0){
    edge = new Cedge[newedges];

    if(DEBUG>=2)
      for(unsigned int i = 0; i < newedges; i++)
	assert(edge[i].mark == 0);

  } else {
    Cedge *origedge = edge;
    edge = new Cedge[newedges];
    if(DEBUG>=2)
      for(unsigned int i = 0; i < newedges; i++)
	assert(edge[i].mark == 0);

    for(register unsigned int i=0; i < maxedges;i++)
      edge[i] = origedge[i];
    delete [] origedge;
  }
  maxedges = newedges;
}

void maxedgealloc(unsigned int num, unsigned int &maxedges, CedgeA* &edge)
{
  unsigned int newedges = num + 1;
  if(newedges <= maxedges)
    return;

  if(newedges < 1024)
    newedges = 1024;
  if(maxedges > 0 && newedges < min(MASK_LL(32), (maxedges+1024)*5LL/4))
    newedges = min(MASK_LL(32), (maxedges+1024)*5LL/4);

  if(VERB/* HERE HERE >=2 */){
    printf("maxedgealloc: num=%u,maxedges=%u: newedges= %u\n",num,maxedges,newedges);
    fflush(stdout);
  }

  if(maxedges <= 0){
    edge = new CedgeA[newedges];

    if(DEBUG>=2)
      for(unsigned int i = 0; i < newedges; i++)
	assert(edge[i].mark == 0);

  } else {
    CedgeA *origedge = edge;
    edge = new CedgeA[newedges];
    if(DEBUG>=2)
      for(unsigned int i = 0; i < newedges; i++)
	assert(edge[i].mark == 0);

    for(register unsigned int i=0; i < maxedges;i++)
      edge[i] = origedge[i];
    delete [] origedge;
  }
  maxedges = newedges;
}

