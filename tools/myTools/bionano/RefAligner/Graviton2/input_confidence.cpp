#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>

#include "constants.h"
#include "conf_params.h" //for globals

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/trunk/input_confidence.cpp 3927 2015-07-08 00:54:12Z tanantharaman $");

static char buf[LINESIZ];


void input_confidence(char* filename) {

  register char *pt;
  char *qt;
  FILE *fp;

  if( !(qt = strstr(filename,".txt")) || strlen(qt) != strlen(".txt") ){ //or .txt
    printf("Error: Input file must end in \".txt\" %s\n",filename);
    fflush(stdout);
    exit(1);
  }

  if((fp = fopen(filename,"r"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("Error: Failed to read input file %s:errno=%d:%s\n",filename,eno,err);
    fflush(stdout);
    exit(1);
  }

  //first null all the data in globals in conf_params.h because it needs to be replaced
  for( int i=0; i<sv_conflen; i++ ) {
    sv_ppvbin_conf[i] = 0; 
    for( int j=0; j<sv_sizelen; j++ ) {
      if( i == 0 )
	sv_ppvbin_sizekb[j] = 0;		 
      sv_ppvbin_del[i][j] = 0;
      sv_ppvbin_ins[i][j] = 0;
    }
  }

  //the way you tell where to put the data is by the order:
  // first is conf bins, then size bins, then ppv table for deletions then insertions
  // for the table, the number of rows is the nubmer of conf bins and
  //  the number of columns is the number of size bins
  //also note blank lines are not allowed

  int currsize = SVCONFLEN; //this is max size in conf_params.h
  int tmplen = 0;
  int   *tmpi = new int[currsize];
  float *tmpf = new float[currsize];
  int rowi = -1, datai = 0; //which data we are currently reading

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

    rowi++; //move here from below
    tmplen = 0;
    int v = 0;
    float w = 0;
    while(qt) { //number of columns is variable, so use while
      if( datai < 2 ) { //first two lines are integers for bin--NOTE: may want these to be floats
	v = strtol(pt,&qt,10);
	//printf("v=%i  ", v); //debug
      } else {
	w = strtod(pt,&qt);
	//printf("w=%.2f  ", w); //debug
      }
      if(pt==qt){ //this happens at end of line, so it's not an error
      	//printf("Error: Invalid value on line %d of %s\n:%s\n",linecnt,filename,buf);
      	//exit(1);
	break;
      }
      if( datai < 2 ) {
	assert(tmplen < currsize);
	tmpi[tmplen] = v;
	tmplen++;
      } else {
	assert(tmplen < currsize);
	tmpf[tmplen] = w;
	tmplen++;
      }
      pt = qt;
    } //end while qt
    //printf("\n"); //debug

    if( datai == 0 ) { //copy tmpi into conf bins
      assert(tmplen < SVCONFLEN); //array max size
      sv_conflen = tmplen;
      for( int i=0; i<tmplen; i++ )
	sv_ppvbin_conf[i] = tmpi[i];
      datai++;
      rowi = -1; //need when move rowi++ above
    }
    else if( datai == 1 ) { //size bins
      assert(tmplen < SVSIZELEN); //array max size
      sv_sizelen = tmplen;
      for( int i=0; i<tmplen; i++ )
	sv_ppvbin_sizekb[i] = tmpi[i];
      datai++;
      rowi = -1; //need when move rowi++ above
    }
    else if( datai == 2 || datai == 3 ) { //deletion, insertion ppv
      if(tmplen != sv_sizelen) { //each row is the size bins
	printf("ERROR: number of columns is %i, but must be number of size bins (%i): on line %i of %s\n", tmplen, sv_sizelen, linecnt, filename);
	fflush(stdout);
	exit(1);
      }
      assert(rowi < SVCONFLEN); //check row number is in bounds--shouldn't be possible due to below if
      for( int i=0; i<tmplen; i++ ) {
	if( datai == 2 )
	  sv_ppvbin_del[rowi][i] = tmpf[i];
	else if( datai == 3 )
	  sv_ppvbin_ins[rowi][i] = tmpf[i];
      }
      if( rowi == sv_conflen-1 ) { //go to deletions after all rows of insertions are read
	if( datai == 3 ) //this means all data is done being read
	  break;
	datai++;
	rowi = -1; //was 0
      }
      //else 
      //rowi++;
    }
    else //this should not happen
      assert(0); //dummy check
  } //end for line in fp  

  //must check that number of rows read is correct, but must do this after done, outside loop on file lines
  if( rowi != sv_conflen-1 ) {
    printf("ERROR: number of rows expected (for last table) is %i, but end of file reached with only %i rows found\n", sv_conflen, rowi+1);
    fflush(stdout);
    //assert(rowi == sv_conflen-1);
    exit(1); //exit is nicer than assert
  }

  //print result
  printf("Parameters read for confidence in %i size bins and %i conf bins from %s\n", sv_sizelen, sv_conflen, filename);

  if(0) { //debug: print all data
    //sv_conflen //printed above
    printf("conf vals: ");
    for( int i=0; i<sv_conflen; i++ )
      printf("%i ", sv_ppvbin_conf[i]);
    printf("\nsize vals: ");
    //sv_sizelen //printed above
    for( int i=0; i<sv_sizelen; i++ )
      printf("%i ", sv_ppvbin_sizekb[i]);
    printf("\ndeletion ppv:\n");
    for( int i=0; i<sv_conflen; i++ ) {
      for( int j=0; j<sv_sizelen; j++ )
	printf("%.2f ", sv_ppvbin_del[i][j]);
      printf("\n");
    }
    printf("insertion ppv:\n");
    for( int i=0; i<sv_conflen; i++ ) {
      for( int j=0; j<sv_sizelen; j++ )
	printf("%.2f ", sv_ppvbin_ins[i][j]);
      printf("\n");
    }
  }

  fflush(stdout);
} //end input_confidence
