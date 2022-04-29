#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include <math.h>

#include "globals.h"
#include "parameters.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/input_refmap.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

static char buf[LINESIZ];

static int warning = 0;

#define CONTIG_SIZE 1000000.0 /* add this size offset (in kb) for each contig to seperate them */

void input_refmap(char *filename)
{
  register char *pt;
  char *qt;
  FILE *fp;
  
  if((fp = fopen(filename,"r"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("Failed to read input file %s:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  /* set up translation table from Gmap[i]->id to i */
  int maxmapid = 0;
  for(register int i=nummaps; --i >= 0;)
    if(Gmap[i]->id > maxmapid)
      maxmapid = Gmap[i]->id;
  int *id2index = new int[maxmapid+1];
  for(register int i=maxmapid; i >= 0; i--)
    id2index[i] = -1;
  for(register int i=nummaps; --i >= 0;)
    id2index[Gmap[i]->id] = i;

  /* reset known ground truth to "unknown" */
  for(register int i=nummaps; --i >= 0;){
    register Cmap *p = Gmap[i];
    p->startloc = p->endloc = p->startmatch = p->endmatch = -1000.0;
    p->fpcnt = p->flipped = -1;
    p->Qscore = -9999.999;
  }

  int linecnt=1;
  int loccnt = 0;
  for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++){
    int len = strlen(buf);
   
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      exit(1);
    }

    if(buf[0] == '#')// comment lines are skipped
      continue;

    if(!strncmp(buf,"Software version:",strlen("Software version")))
      continue;

    if(!strncmp(buf,"MappedMoleculeId",strlen("MappedMoleculeId")))
      continue;

    pt = buf;
    int MappedMoleculeId = strtol(pt,&qt,10);
    if(pt == qt){
      printf("Invalid integer value for MappedMoleculeId=%d on line %d of %s\n%s\n",
	     MappedMoleculeId,linecnt,filename,buf);
      exit(1);
    }
    
    pt = qt;
    int MoleculeId = strtol(pt,&qt,10);
    if(pt == qt || MoleculeId < 0){
      printf("Invalid integer value for MoleculeId on line %d of %s\n%s\n",linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    int MoleculeIndex = strtol(pt,&qt,10);
    if(pt == qt){
      printf("Invalid integer value for MoleculeIndex=%d on line %d of %s\n%s\n",
	     MoleculeIndex,linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    int ContigId = strtol(pt,&qt,10);
    if(pt == qt){
      printf("Invalid integer value for ContigId=%d on line %d of %s\n%s\n",
	     ContigId,linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    double Score = strtod(pt,&qt);
    if(pt == qt){
      printf("Invalid float value for Score on line %d of %s\n%s\n",linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    double Confidence = strtod(pt,&qt);
    if(pt == qt){
      printf("Invalid float value for Confidence=%0.3f on line %d of %s\n%s\n",
	     Confidence,linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    while(*pt && isspace(*pt))
      pt++;
    int orientation = (*pt == '-') ? 1 : (*pt == '+') ? 0 : -1;
    if(orientation < 0){
      printf("Invalid Direction \'%c\' on line %d of %s\n%s\n", *pt, linecnt, filename, buf);
      exit(1);
    }

    pt++;
    int StartLocation = strtol(pt,&qt,10);
    if(pt == qt){
      printf("Invalid integer value %d for StartLocation on line %d of %s\n%s\n",
	     StartLocation,linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    int EndLocation = strtol(pt,&qt,10);
    if(pt == qt){
      printf("Invalid integer value for EndLocation on line %d of %s\n%s\n",linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    int StartMatchLocation = strtol(pt,&qt,10);
    if(pt == qt){
      printf("Invalid integer value for StartMatchLocation on line %d of %s\n%s\n",linecnt,filename,buf);
      exit(1);
    }

    pt = qt;
    int EndMatchLocation = strtol(pt,&qt,10);
    if(pt == qt){
      printf("Invalid integer value for EndMatchLocation on line %d of %s\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    
    if(MoleculeId > maxmapid){
      printf("MoleculeId=%d is too large on line %d of %s : (largest value is %d)\n",
	     MoleculeId,linecnt,filename,maxmapid);
      exit(1);
    }
    int mapid = id2index[MoleculeId];
    if(mapid < 0){
      if(!warning){
	printf("WARNING:MoleculeId=%d(MappedMoleculeId=%d,MoleculeIndex=%d,ContigId=%d) on line %d of %s is not present in input file\n",MoleculeId,MappedMoleculeId,MoleculeIndex,ContigId,linecnt,filename);
	fflush(stdout);
	warning = 1;
      }
      continue;
    }
    register Cmap *Xmap = Gmap[mapid];
	   
    /*    printf("MappedId=%d,MoleculeId=%d,MoleculeIndex=%d,StartLocation=%d,EndLocation=%d:mapid=%d,Xmap->Qscore=%8.3f,Score=%0.3f,loccnt=%d\n",
	   MappedMoleculeId,MoleculeId,MoleculeIndex,StartLocation,EndLocation,mapid,Xmap->Qscore < -99999.0 ? -99999.0 : Xmap->Qscore,Score,loccnt);
	   fflush(stdout);*/

    if(Xmap->Qscore < Score){/* update startloc,endloc,flipped,Qscore */
      loccnt++;
      //      Xmap->startloc = (SCORES ? EndLocation*0.001 - Xmap->site[0][Xmap->numsite[0]] : StartLocation*0.001);
      Xmap->refid = ContigId;
      Xmap->startloc = StartLocation*0.001 + ContigId*CONTIG_SIZE;
      Xmap->endloc = EndLocation*0.001 + ContigId*CONTIG_SIZE;
      Xmap->flipped = orientation;
      Xmap->Qscore = Score;
      Xmap->startmatch = StartMatchLocation*0.001;
      Xmap->endmatch = EndMatchLocation*0.001;
    }

  }

  register int cnt=0;
  for(register int i=nummaps; --i>=0;)
    if(Gmap[i]->startloc > -1000.0)
      cnt++;
  if(VERB){
    printf("Read in true locations for %d/%d maps\n",cnt,nummaps);
    fflush(stdout);
  }
      
  delete [] id2index;

  fclose(fp);
}

