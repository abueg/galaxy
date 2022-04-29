#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "globals.h"
#include "parameters.h"
#include "bedEntry.h"


static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/input_bed.cpp 3896 2015-06-19 22:31:31Z tanantharaman $");

static char buf[LINESIZ];


//if the input bed file is larger than the initial bedentries size, you need to allocate more memory
//double the size, return the new size
int growBedEntries() {
  bedEntry** tmp = bedentries; //copy pointer to start of old bedentries
  int newsize = numbedentries*2; //new size is twice old size
  bedentries = new bedEntry*[newsize]; //allocate new size in scaffolds
  for(int i=0; i<numbedentries; i++) //copy pointers to existing bed objects, now at tmp
    bedentries[i] = tmp[i];
  delete [] tmp; //tmp is old--container no longer necessary
  return newsize;
}


void input_bed(char* filename) {
  int currsize = 100; //arbitrary number
  bedentries = new bedEntry*[currsize]; //initialize global pointer to pointer of scaffolds with size currscafsize

  register char *pt;
  char *qt;
  FILE *fp;

  if( !(qt = strstr(filename,".bed")) || strlen(qt) != strlen(".bed") ){
    printf("Input scaffold file must end in \".bed\" %s\n",filename);
    exit(1);    
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

    //first column is contig id
    int contigid = strtol(pt,&qt,10);
    if(pt==qt){
      printf("unable to read contigID on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    //columns 2 and 3 are start and stop
    pt = qt;
    double start = strtod(pt,&qt);
    if(pt==qt){
      printf("unable to read start on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;
    double stop = strtod(pt,&qt);
    if(pt==qt){
      printf("unable to read stop on line %d of %s\n:%s\n",linecnt,filename,buf);
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
    char* name = strdup(pt); //copies pt
    name[strlen(pt)-strlen(qt)+1] = '\0'; //need to terminate after first 
    pt += strlen(name);

    //contigid must be > 0, start and stop must be >= 0
    if( contigid < 1 ) {
      printf("Invalid contigid, must be > 0, on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    if( start < 0 ) {
      printf("Invalid start, must be >= 0, on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    if( stop < 0 ) {
      printf("Invalid stop, must be >= 0, on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }
    //fstop must always be > fstart
    if( stop < start ) {
      printf("Invalid stop, must be >= start, on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }

    //make new bedEntry object, insert in bedentries
    if( numbedentries == currsize ) {
      currsize = growBedEntries();
    }
    bedentries[numbedentries] = new bedEntry(contigid, start, stop, name);
    numbedentries++; //global size of array
  } //end for line in fp  

  //print result
  printf("Read %i records from %s\n", numbedentries, filename);

  if(0) { //debug: print all bedEntries
    for(int bi=0; bi<numbedentries; bi++)
      printf("bedEntry %2i:  %2i  %9.0f  %9.0f  %s\n", bi+1,
	     bedentries[bi]->contigid,
	     bedentries[bi]->start,
	     bedentries[bi]->stop,
	     bedentries[bi]->type);
  }
} //end input_bed
