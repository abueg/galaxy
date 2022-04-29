#ifndef BEDENTRY_H
#define BEDENTRY_H

#include "ident.h"

static Ident bedEntry_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/bedEntry.h 5409 2016-10-03 19:32:08Z tanantharaman $");

/* Bed files have up to 12 columns,
   but for now, I only need the first four */
class bedEntry {
public:
  bedEntry() {
    contigid = 0;
    start = 0;
    stop = 0;
    type = 0;
  }

  bedEntry(int cid, double st, double sp, char* tp) {
    contigid = cid;
    start = st;
    stop = sp;
    type = tp; /**< assign pointer--assume fhead arg is already copied */
  }

  int contigid; //contig id (or chromosome): first column of bed file
  double start; //although stored in bed as int, read as double to compare to positions of labels/SVs
  double stop; //third column
  char* type; //fourth column
};


#endif
