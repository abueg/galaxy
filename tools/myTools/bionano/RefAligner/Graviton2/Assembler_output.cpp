#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include <math.h>

#include "constants.h"
#include "globals.h"
#include "Assembler_parameters.h"
#include "Assembler.h"
#include "Cmap.h"
#include "Cnode.h"
#include "Ccontig.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Assembler_output.cpp 3896 2015-06-19 22:31:31Z tanantharaman $");

void output_TrimmedMaps(int nummaps, Cmap **map, Cnode *node, char *basename)
{
  char filename[PATH_MAX];
  strcpy(filename,basename);
  int i = strlen(filename);
  sprintf(&filename[i],"_trim.vfixx");

  FILE *fp;
  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("failed to open file %s for writing trimmed vfixx file:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }
  if(VERB){
    printf("Generating %s (trimmed Maps)\n",filename);
    fflush(stdout);
  }
  
  if(colors > 1){
    printf("output_TrimmedMaps: not implemented for multiple colors (%d)\n", colors);
    exit(1);
  }

  /* output key header lines */
  fprintf(fp,"#\tN Colors:\t%d\n", colors);
  fprintf(fp,"#\tN Fragments:\t%d\n", nummaps);
  fprintf(fp,"#\tBases per Pixel:\t%d\n", (int)floor(PixelLen*1000.0+0.5));
  fprintf(fp,"#\t>0	Mol ID	Start (Bp)	End (Bp)	Length (Bp)    trimL    trimR\n");
  for(register int c = 0; c < colors; c++)
    fprintf(fp,"#\t>%d	Color %d	Nick Id (In untrimmed Nanomap)\n",c+1,c+1);
  fprintf(fp,"#\tFrag Len (px)	Nick Locs (px)...\n");
  for(register int i = 0; i < nummaps;i++){
    register Cmap *pmap = map[i];
    register Cnode *pnode = &node[i];
    //    register int M = pmap->numsite[0];
    register FLOAT *X = pmap->site[0];
    register int trimL = pnode->trimL;
    register int trimR = pnode->trimR;
    /* Output reference location line */
    fprintf(fp,">0\t%lld\t%d\t%d\t%d\t%d\t%d\n", 
	    pmap->id,  /* Mol ID */
	    (int)floor(pmap->startloc*1000.0+0.5), /* Start (Bp) */
	    (int)floor(pmap->endloc*1000.0+0.5),   /* End (Bp) */
	    (int)floor(pmap->len*1000.0+0.5),      /* Length (Bp) */
	    trimL, trimR);
    /* output original Nick ID (before trimming) */
    fprintf(fp,">1");
    for(register int k = trimL+1; k < trimR; k++)
      fprintf(fp,"\t%d", k);
    fprintf(fp,"\n");
    /* output Nick locations (in Pixels) */
    fprintf(fp,"%0.3f", (X[trimR] - X[trimL])/PixelLen);
    for(register int k = trimL+1; k < trimR; k++)
      fprintf(fp,"\t%0.3f", (X[k] - X[trimL])/PixelLen);
    fprintf(fp,"\n");
  }
  FILEclose(fp);
}

