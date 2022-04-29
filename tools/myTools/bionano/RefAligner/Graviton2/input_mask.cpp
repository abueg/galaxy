#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

//#include "globals.h"
//#include "parameters.h"
//#include "bedEntry.h"
#include "structuralVariation.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/input_mask.cpp 3927 2015-07-08 00:54:12Z tanantharaman $");

static char buf[LINESIZ];


void input_mask(char* filename) {
  int currsize = 60; //arbitrary number
  //bedentries = new bedEntry*[currsize]; //initialize global pointer to pointer of scaffolds with size currscafsize
  maskentries = new maskSV[currsize];

  register char *pt;
  char *qt;
  FILE *fp;

  if( !(qt = strstr(filename,".mask")) || strlen(qt) != strlen(".mask") ){ //either .mask
    if( !(qt = strstr(filename,".txt")) || strlen(qt) != strlen(".txt") ){ //or .txt
      printf("Input mask file must end in \".mask\" or \".txt\" %s\n",filename);
      exit(1);
    }
  }

  if((fp = fopen(filename,"r"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("Failed to read input file %s:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  //start reading file
  //note: LINESIZ is in constants.h, and is 256*1024; buf is in RefAligner
  int linecnt=1;
  for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++){
    if(buf[0] == '#') /* all comment lines are ignored */
      continue;

    pt = &buf[0];
    while(*pt && isspace(*pt)) //pt should be beginning of first column, stripped of opening whitespace
      pt++;
    qt = pt;

    //first column is index (not used)
    (void)strtol(pt,&qt,10);
    if(pt==qt){
      printf("unable to read contigID on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    //columns 2 and 3 are chromosomes
    pt = qt;
    int chr1 = strtol(pt,&qt,10);
    if(pt==qt){
      printf("unable to read Chr1 on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;
    int chr2 = strtol(pt,&qt,10);
    if(pt==qt){
      printf("unable to read Chr2 on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    //columns 4 and 5 are positions
    pt = qt;
    double chr1_start = strtod(pt,&qt);
    if(pt==qt){
      printf("unable to read Chr1Pos on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;
    double chr2_start = strtod(pt,&qt);
    if(pt==qt){
      printf("unable to read Chr2Pos on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    //columns 6 and 7 are SDs for positions (not required)
    /*
    pt = qt;
    double chr1_startsd = strtod(pt,&qt);
    if(pt==qt){
      printf("unable to read Chr1PosSD on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;
    double chr2_startsd = strtod(pt,&qt);
    if(pt==qt){
      printf("unable to read Chr2PosSD on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    //column 8 is N samples
    pt = qt;
    int nsamples = strtol(pt,&qt);
    if(pt==qt){
      printf("unable to read N on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    while(*pt && isspace(*pt)) //pt should be beginning of first column, stripped of opening whitespace
      pt++;
    qt = pt;
    while(!isspace(*qt))
      qt++;
    if(qt == pt){
      printf("Invalid bed line %d of %s:\n%s\n", linecnt, filename, buf);
      exit(1);
    }
    qt--; //you get an extra one with this
    char* name = strdup(pt); //copies pt--use when pt is a string to keep
    name[strlen(pt)-strlen(qt)+1] = '\0'; //need to terminate after first 
    pt += strlen(name);
    */

    //chrs must be > 0, starts must be >= 0
    if( chr1 < 1 ) {
      printf("Invalid Chr1, must be > 0, on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    if( chr2 < 1 ) {
      printf("Invalid Chr2, must be > 0, on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    if( chr1_start < 0 ) {
      printf("Invalid Chr1Pos, must be >= 0, on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    if( chr2_start < 0 ) {
      printf("Invalid Chr2Pos, must be >= 0, on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }

    //make new object, insert in array
    if( nummaskentries == currsize ) {
      currsize = growArray(maskentries, nummaskentries);
    }
    maskentries[nummaskentries] = maskSV(chr1, chr2, chr1_start, chr2_start);
    nummaskentries++; //n ele of array used
  } //end for line in fp  

  //print result
  printf("Read %i records from %s\n", nummaskentries, filename);

  if(0) { //debug: print all maskSV objects
    for(int bi=0; bi<nummaskentries; bi++)
      printf("maskSV %2i:  %2i  %2i  %12.2f  %12.2f\n", bi+1,
	     maskentries[bi].chr1,
	     maskentries[bi].chr2,
	     maskentries[bi].chr1_start,
	     maskentries[bi].chr2_start);
  }
} //end input_bed
