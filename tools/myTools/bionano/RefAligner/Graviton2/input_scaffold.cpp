#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "globals.h"
#include "parameters.h"


static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/input_scaffold.cpp 3896 2015-06-19 22:31:31Z tanantharaman $");
static char buf[LINESIZ];


/** class which stores the contents of the scaffold file */
class scaffold {
public:

  scaffold(){
    fastahead  = 0;
    fastaid    = 0;
    contigid   = 0;
    scaffoldid = 0;
    fastastart = 0;
    fastastop  = 0;
  }

  scaffold(char* fhead, int fid, int cid, int sid, double fstart, double fstop){
    fastahead  = fhead; /**< assign pointer--assume fhead arg is already copied */
    fastaid    = fid;
    contigid   = cid;
    scaffoldid = sid;
    fastastart = fstart;
    fastastop  = fstop;
  }

  char*  fastahead;
  int    fastaid;
  int    contigid;
  int    scaffoldid;
  double fastastart;
  double fastastop;
};


//remember, the contigid differs from the index in the scaffold array by 1

//get the gap size between two scaffold contigs
//this is the fastastart of the second less the fastastop of the first
double scaffoldGapSize(int id1, int id2) {
  if( id1 == id2 )
    return 0; //no gap if same ids
  assert(id1 <= numscafentries && id2 <= numscafentries); //protect against seg fault
  //this is the condition to size the gap: must be on same scaffold, ie, same fastaid
  if( scaffolds[id1-1]->fastaid == scaffolds[id2-1]->fastaid ) {
    if( id1 > id2 )
      return scaffolds[id1-1]->fastastart - scaffolds[id2-1]->fastastop ;
    else
      return scaffolds[id2-1]->fastastart - scaffolds[id1-1]->fastastop ;
  }
  else //this means no gap sized bc different scaffolds--above can't be negative
    return -1;
}


//get the start position of the entry in the global scaffolds list
// (which is index one less than the id), and add offset to it
double scaffoldStartPosition(int id, double offset) {
  assert(id <= numscafentries); //protect against seg fault in next line
  return scaffolds[id-1]->fastastart + offset;
}


//get the distance from a position on a contig to its end
//the end is the total length, which is just fastastop - fastastart; then do length - startpos
double scaffoldEndLength(int id, double startpos) {
  assert(id <= numscafentries); //protect against seg fault in next line
  scaffold* sc = scaffolds[id-1];
  return sc->fastastop - sc->fastastart - startpos;
}


//if the input scaffold file is larger than the initial scaffolds size, you need to allocate more memory
//double the size, return the new size
int growScaffolds() {
  scaffold** tmpscaf = scaffolds; //copy pointer to start of old scaffolds
  int newsize = numscafentries*2; //new size is twice old size
  scaffolds = new scaffold*[newsize]; //allocate new size in scaffolds
  for(int i=0; i<numscafentries; i++) //copy pointers to existing scaffold objects
    scaffolds[i] = tmpscaf[i];
  delete [] tmpscaf; //tmpscaf is old scaffolds--container no longer necessary
  return newsize;
}


//fn which reads scaffold file and populates <global array of scaffold objects>
void input_scaffold(char *filename)
{
  int currscafsize = 100; //arbitrary number
  scaffolds = new scaffold*[currscafsize]; //initialize global pointer to pointer of scaffolds with size currscafsize

  register char *pt;
  char *qt;
  FILE *fp;

  if( !(qt = strstr(filename,".scaffold")) || strlen(qt) != strlen(".scaffold") ){
    printf("Input scaffold file must end in \".scaffold\" %s\n",filename);
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

    //first column is fastaHead, which is a string
    pt = &buf[0];
    while(*pt && isspace(*pt)) //pt should be beginning of first column, stripped of opening whitespace
      pt++;
    qt = pt;
    while(!isspace(*qt))
      qt++;
    if(qt == pt){
      printf("Invalid scaffold line %d of %s:\n%s\n", linecnt, filename, buf);
      exit(1);
    }
    qt--; //you get an extra one with this
    char* fastahead = strdup(pt); //copies pt
    fastahead[strlen(pt)-strlen(qt)+1] = '\0'; //need to terminate after first 
    pt += strlen(fastahead);

    //second column is fastaID
    int fastaid = strtol(pt,&qt,10);
    if(pt==qt){
      printf("unable to read fastaID on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    int contigid = strtol(pt,&qt,10);
    if(pt==qt){
      printf("unable to read contigID on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    int scafid = strtol(pt,&qt,10);
    if(pt==qt){
      printf("unable to read scafID on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }

    //for fasta start and stop, use double even though they're normally integers
    pt = qt;
    double fstart = strtod(pt,&qt);
    if(pt==qt){
      printf("unable to read fstart on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    double fstop = strtod(pt,&qt);
    if(pt==qt){
      printf("unable to read fstop on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }

    //make new scaffold object from above
    // scaffold(char* fhead, int fid, int cid, int sid, double fstart, double fstop){
    int si = numscafentries; //just for readability
    if( si == currscafsize ) { //index must be < size
      currscafsize = growScaffolds();
    }
    scaffolds[si] = new scaffold(fastahead, fastaid, contigid, scafid, fstart, fstop);

    //consistency checks

    //fstop must always be > fstart: current entry
    if( fstop < fstart ) {
      printf("Invalid fastaStop, must be < fastaStart, on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1);
    }

    //first contigID must be 1--strictly speaking, this isn't necessary, however,
    // it is necessary if the index == contigID-1, and I've designed the file format this way
    if( si == 0 && scaffolds[si]->contigid != 1 ) {
      printf("First contigID must be 1 on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1); 
    }

    //contigID is always one more than previous
    if( si > 0 && scaffolds[si]->contigid != scaffolds[si-1]->contigid+1 ) {
      printf("ContigIDs are not sequential on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1); 
    }

    //if same fastaHead, fastaID must also be same (strcmp == 0 for equal)
    if( si > 0 && !strcmp(scaffolds[si]->fastahead, scaffolds[si-1]->fastahead) &&
	scaffolds[si]->fastaid != scaffolds[si-1]->fastaid ) {
      printf("fastaHead same, but fastaID different on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1); 
    }

    //and vice-versa
    if( si > 0 && strcmp(scaffolds[si]->fastahead, scaffolds[si-1]->fastahead) &&
	scaffolds[si]->fastaid == scaffolds[si-1]->fastaid ) {
      printf("fastaID same, but fastaHead different on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1); 
    }

    //if same fastaID, scaffoldID must be one more than previous
    if( si > 0 && scaffolds[si]->fastaid == scaffolds[si-1]->fastaid &&
	scaffolds[si]->scaffoldid != scaffolds[si-1]->scaffoldid+1 ) {
      printf("fastaID same, but non-sequential scaffoldIDs on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1); 
    }

    //if same fastaID, fastastop must be > previous fasta start (allow equal)
    if( si > 0 && scaffolds[si]->fastaid == scaffolds[si-1]->fastaid &&
	scaffolds[si]->fastastart < scaffolds[si-1]->fastastop ) {
      printf("fastaID same, but out of order fastaStart/fastaStop on line %d of %s\n:%s\n",linecnt,filename,buf);
      exit(1); 
    }

    numscafentries++; //global; size of scaffolds array

    //if( linecnt > 100 ) break; //debug--only read first N lines
  } //end loop on lines in scaffold file

  if(0) { //debug--print every entry of scaffolds object
    for( int si=0; si<numscafentries; si++ ) {
      printf("entry %3i = %s %i %i %i %9.1f %9.1f;\n", si+1,
	     scaffolds[si]->fastahead,
	     scaffolds[si]->fastaid,
	     scaffolds[si]->contigid,
	     scaffolds[si]->scaffoldid,
	     scaffolds[si]->fastastart,
	     scaffolds[si]->fastastop);
    }
  }
}
