#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#ifndef WIN32
#include <unistd.h>
#else
#include <direct.h>
#define getcwd _getcwd
#endif

#include "globals.h"
#include "parameters.h"


static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/output_amap.cpp 3896 2015-06-19 22:31:31Z tanantharaman $");

//just two ints--basically the stl pair, but ints only (not templated)
class pairints {
public:
  pairints() {
    first = 0;
    second = 0;
  }

  pairints(int f, int s) {
    first = f;
    second = s;
  }

  int first;
  int second;
};


//double the size of the argument array, return the new size
int growPairints(pairints*& ina, int oldsize) {
  pairints* tmp = ina; //copy pointer to start of old array
  int newsize = oldsize*2; //new size is twice old size
  ina = new pairints[newsize]; //allocate new size in old object ina
  for(int i=0; i<oldsize; i++) //copy pointers to existing scaffold objects
    ina[i] = tmp[i];
  delete [] tmp; //tmp is old array--container no longer necessary
  return newsize;
}


//find all sizable gaps
//use global xmapentries which is populated in output_xmap
//integer argument is number of sized gaps--it ends up being length of returned array
pairints* findGapPairs(int& npairs) {
  npairs = 0; //start at zero and add as found
  int newpairsize = 10; //arbitrary number
  pairints* newpairs = new pairints[newpairsize];

  for(int i=0; i<numxmapentries; i++) {
    for(int j=i+1; j<numxmapentries; j++) {
      if( xmapentries[i]->qrycontigid == xmapentries[j]->qrycontigid ||
	  xmapentries[i]->refcontigid == xmapentries[j]->refcontigid ) {
	newpairs[npairs] = pairints(i,j);
	npairs++;
	if( npairs > newpairsize )
	  newpairsize = growPairints(newpairs, newpairsize);
      }  
    }
  }
  return newpairs;
} //end findGapPairs


//process xmapentries object, derive the amap from it, and write it to disk
void output_amap(char *basename) {

  if(strstr(basename,"/dev/null"))
    return;

  //output amap name
  char filename[PATH_MAX];
  strcpy(filename,basename);
  int i = strlen(filename);
  strcpy(&filename[i],".amap");
  //also make the xmap file name
  char xfilename[PATH_MAX];
  strcpy(xfilename,basename);
  strcpy(&xfilename[i],".xmap"); //i is same as above

  if(checkFile(filename))
    return;

  //open output file
  if(VERB){
    printf("Generating %s\n",filename);
    fflush(stdout);
  }

  FILE *fp;
  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("failed to open file %s for writing amap file:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  /* write out commandline */
  printversion(fp);

  //sorting of alignment should remain after output_xmap call

  fprintf(fp,"# AMAP File Version:\t0.0\n"); //use 0.1 once this is fully described

  //take all this from output_xmap even though currently it only runs for pairsplit, which is only for pairwise, ie, no spots
  if(spots_filename[0]){/* Reference Aligner : refer to condensed version of reference CMAP */
    fprintf(fp,"# Reference Maps From:\t%s_r.cmap\n", basename);
    fprintf(fp,"# Query Maps From:\t%s",vfixx_filename[0]);/* This is actually a modified cmap file */
    if(DEBUG) assert(num_files==1);
  } else {
    fprintf(fp,"# Reference Maps From:\t%s", vfixx_filename[0]);
    for(register int i = 1; i < RefIndex; i++)
      fprintf(fp,",%s",vfixx_filename[i]);
    fprintf(fp,"\n");
    fprintf(fp,"# Query Maps From:\t%s",vfixx_filename[QueryIndex]);
    for(register int i = QueryIndex+1; i < num_files;i++)
      fprintf(fp,",%s",vfixx_filename[i]);
  }
  fprintf(fp,"\n");
  fprintf(fp,"# Xmap Entries From:\t%s\n",xfilename);

  //header
  fprintf(fp,"# GapIndex  QryContig1  QryContig2  RefContig1  RefContig2  QryStart  QryStop  QryGapSize  QryGapSizeErr  RefStart  RefStop  RefGapSize  RefGapSizeN  Confidence\n");

  //process the xmap data to derive the amap data

  const int colw = 15; //actually width is probably less--use this for char arrays
  int npairs = 0;
  pairints* newpairs = findGapPairs(npairs);

  for(int i=0; i<npairs; i++) {
    int id1 = newpairs[i].first;
    int id2 = newpairs[i].second;

    int qrycontig1 = xmapentries[id1]->qrycontigid;
    int qrycontig2 = xmapentries[id2]->qrycontigid;
    int refcontig1 = xmapentries[id1]->refcontigid;
    int refcontig2 = xmapentries[id2]->refcontigid;

    int gaptype = 0; //gaptype determines which (query or ref) gaps get sized
    if( qrycontig1 == qrycontig2 )
      gaptype = 1; //1 here means 1 or 2--have qry gap size (ref gap size in case 1, but not 2)
    else if( refcontig1 == refcontig2 )
      gaptype = 3; //3 means that the query gap has no size bc different contigs;
    else {
      printf("ERROR in output_amap: unknown gap type\n");
      break;
    }

    double prevqrystart = xmapentries[id1]->qrystartpos;
    double prevqrystop  = xmapentries[id1]->qryendpos;
    double currqrystart = xmapentries[id2]->qrystartpos;
    double currqrystop  = xmapentries[id2]->qryendpos;
    double qrystart = 0, qrystop = 0;
    char qrygap[colw];

    if( gaptype == 1 ) {
      if( max(prevqrystart, prevqrystop) > max(currqrystart, currqrystop) ) {
	qrystart = min(prevqrystart, prevqrystop);
	qrystop  = max(currqrystart, currqrystop);
      }
      else {
	qrystart = max(prevqrystart, prevqrystop);
	qrystop  = min(currqrystart, currqrystop);
      }
      sprintf(qrygap, "%9.1f", fabs(qrystop - qrystart));
    }
    else { //gaptype == 3
      qrystart = max(prevqrystart, prevqrystop);
      qrystop  = min(currqrystart, currqrystop);
      sprintf(qrygap, "%9s", "N");
    }

    //now process ref positions

    double refstart = scaffoldStartPosition(refcontig1, xmapentries[id1]->refendpos);
    double refstop  = scaffoldStartPosition(refcontig2, xmapentries[id2]->refstartpos);
    double refsizen = scaffoldGapSize(refcontig1, refcontig2);
    if( refsizen == -1 ) //unsized reference gap means different scaffolds means gap type 2
      gaptype = 2;

    char refgap[colw];
    char refsizenstr[colw];
    if( gaptype == 1 ) {
      //RefGapSize is refGapSizeN _plus_ the un-aligned portion between the last alignment and the end of the contig
      //since the xmap is sorted by reference position, you always go to the end of the first and beginning of second
      double gap1 = scaffoldEndLength(refcontig1, xmapentries[id1]->refendpos);
      double gap2 = xmapentries[id2]->refstartpos;
      sprintf(refgap     , "%9.1f", gap1 + gap2 + refsizen);
      sprintf(refsizenstr, "%9.1f", refsizen);
    }
    else if( gaptype == 2 ) { //here, we know nothing bc different fasta ids
      sprintf(refgap     , "%9s", "N");
      sprintf(refsizenstr, "%9s", "N");
    }
    else if( gaptype == 3 ) { //here, same contig, so simply difference
      sprintf(refgap     , "%9.1f", refstop-refstart);
      sprintf(refsizenstr, "%9s", "N");
    }
    else { //this can never happen
      printf("ERROR in output_amap: invalid gaptype %i", gaptype);
      exit(1);
    }

    double conf = xmapentries[id1]->confidence + xmapentries[id2]->confidence;

    //printf("%3i  %3i  %3i  %3i  %3i  %9.1f  %9.1f  %s   -1  %10.1f  %10.1f  %s  %s\n", i+1, qrycontig1, qrycontig2, refcontig1, refcontig2, qrystart, qrystop, qrygap, refstart, refstop, refgap, refsizenstr); //debug
    fprintf(fp,"%3i  %3i  %3i  %3i  %3i  %9.1f  %9.1f  %s   -1  %10.1f  %10.1f  %s  %s  %6.2f\n", i+1, qrycontig1, qrycontig2, refcontig1, refcontig2, qrystart, qrystop, qrygap, refstart, refstop, refgap, refsizenstr, conf);
  }

  FILEclose(fp);
}
