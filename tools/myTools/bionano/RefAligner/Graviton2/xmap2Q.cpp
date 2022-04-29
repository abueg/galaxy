#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#ifndef WIN32
#include <unistd.h>
#else
#include <direct.h>
#define getcwd _getcwd
#endif

#include "constants.h"

const char *SVN_ID = (char *)"$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/xmap2Q.cpp 1264 2013-03-13 19:37:30Z tanantharaman $";
#include "version.h"


int main(int argc, char **argv)
{
  if(argc < 2 || argc > 3){
    printf("Usage:\n%s assembly.xmap [refmap.xmap]\n",argv[0]);
    printf("%s %s\n",VERSION,SVN_ID);
    exit(1);
  }
  char *filename = argv[1];
  char *refname = 0;
  if(argc >= 3)
    refname = argv[2];

  FILE *fp;
  char buf[LINESIZ];  
  if((fp = fopen(filename,"r"))==NULL){
    printf("Cannot read file %s\n",filename);
    exit(1);
  }
  

  int mapcnt = 0;
  double lenSum = 0.0;/* sum of all aligned parts of maps */
  double refStart = 1e+30,refEnd = 0.0;/* range of aligned region in reference (over which average coverage is computed) */
  double confSum = 0.0;/* sum of confidence values for all maps */
  
  long long RefContigID = -1;

  int linecnt=1;
  for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++) {
    int len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      exit(1);
    }

    if(buf[0] == '#'){/* ignore comment lines */
      continue;
    }
    
    int XmapEntryID;
    long long QryContigID, RefcontigID;
    double QryStartPos,QryEndPos,RefStartPos,RefEndPos,Confidence;
    char orientation;
    
    char *pt = buf, *qt;

    XmapEntryID = strtol(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || XmapEntryID < 0){
      printf("Invalid integer value for XmapEntryID on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    QryContigID = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || QryContigID < 0){    
      printf("Invalid long long value for QryContigID on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    RefcontigID = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || RefcontigID < 0){    
      printf("Invalid long long value for RefcontigID on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;
    if(RefContigID < 0)
      RefContigID = RefcontigID;
    else if(RefContigID != RefcontigID){
      printf("Inconsistent RefcontigID=%lld on line %d of %s (previously %lld):\n%s\n",RefcontigID,linecnt,filename,RefContigID,buf);
      exit(1);
    }

    QryStartPos = strtod(pt,&qt);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || QryStartPos < 0){    
      printf("Invalid float value for QryStartPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    QryEndPos = strtod(pt,&qt);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || QryEndPos < 0){    
      printf("Invalid float value for QryEndPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;
    
    RefStartPos = strtod(pt,&qt);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || RefStartPos < 0){    
      printf("Invalid float value for RefStartPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    RefEndPos = strtod(pt,&qt);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || RefEndPos < 0){    
      printf("Invalid float value for RefEndPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    while(*pt && isspace(*pt))
      pt++;
    orientation = *pt;
    if(!orientation || !(orientation == '-' || orientation == '+') || !(pt[1] == 0 || isspace(pt[1]))){
      printf("Invalid Orientation on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt++;
    
    Confidence = strtod(pt,&qt);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || Confidence < 0){    
      printf("Invalid float value for Confidence on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;
    
    mapcnt++;
    lenSum += fabs(QryStartPos - QryEndPos);
    if(min(RefStartPos,RefEndPos) < refStart)
      refStart = min(RefStartPos,RefEndPos);
    if(max(RefStartPos,RefEndPos) > refEnd)
      refEnd = max(RefStartPos,RefEndPos);
    confSum += Confidence;
  }
  fclose(fp);
  
  if(argc < 3){
    if(mapcnt <= 0 || refStart >= refEnd)
      printf("%lld %0.1f 0.00 0.00\n", RefContigID, refEnd-refStart);
    else
      printf("%lld %0.1f %0.2f %0.2f\n", RefContigID, refEnd-refStart, lenSum/(refEnd-refStart), confSum/mapcnt);

    return 0;
  }

  /* locate RefcontigID in file refname */
  if((fp = fopen(refname,"r"))==NULL){
    printf("Cannot read file %s\n", filename);
    exit(1);
  }
  linecnt=1;
  for(;fgets(buf,LINESIZ,fp)!=NULL;linecnt++) {
    int len = strlen(buf);
    if(len >= LINESIZ-1 || (buf[len-1] != '\n' && buf[len-1] != '\r')){
      printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",buf[len-1],filename,linecnt,buf);
      exit(1);
    }

    if(buf[0] == '#')/* ignore comment lines */
      continue;
    
    int XmapEntryID;
    long long QryContigID, RefcontigID = -1;
    double QryStartPos,QryEndPos,RefStartPos,RefEndPos;
    
    char *pt = buf, *qt;

    XmapEntryID = strtol(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || XmapEntryID < 0){
      printf("Invalid integer value for XmapEntryID on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    QryContigID = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || QryContigID < 0){    
      printf("Invalid long long value for QryContigID on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    RefcontigID = strtoll(pt,&qt,10);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || RefcontigID < 0){    
      printf("Invalid long long value for RefcontigID on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    QryStartPos = strtod(pt,&qt);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || QryStartPos < 0){    
      printf("Invalid float value for QryStartPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    QryEndPos = strtod(pt,&qt);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || QryEndPos < 0){    
      printf("Invalid float value for QryEndPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;
    
    RefStartPos = strtod(pt,&qt);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || RefStartPos < 0){    
      printf("Invalid float value for RefStartPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    RefEndPos = strtod(pt,&qt);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || RefEndPos < 0){    
      printf("Invalid float value for RefEndPos on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;

    while(*pt && isspace(*pt))
      pt++;
    char orientation = *pt;
    if(!orientation || !(orientation == '-' || orientation == '+') || !(pt[1] == 0 || isspace(pt[1]))){
      printf("Invalid Orientation on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt++;
    
    double Confidence = strtod(pt,&qt);
    if(pt == qt || !(*qt == 0 || isspace(*qt)) || Confidence < 0){    
      printf("Invalid float value for Confidence on line %d of %s:\n%s\n",linecnt,filename,buf);
      exit(1);
    }
    pt = qt;
    //    printf("line=%d:QryContigID=%lld, target id = %lld\n",linecnt,QryContigID,RefContigID);

    if(QryContigID == RefContigID){
      if(mapcnt <= 0 || refStart >= refEnd)
	printf("%lld %0.1f 0.00 0.00 Mapped %0.1f %0.1f %0.1f %0.1f %c %0.2f\n", 
	       RefContigID, refEnd-refStart, RefStartPos, RefEndPos, min(QryStartPos, QryEndPos),max(QryStartPos,QryEndPos),orientation,Confidence);
      else
	printf("%lld %0.1f %0.2f %0.2f Mapped %0.1f %0.1f %0.1f %0.f %c %0.2f\n", 
	       RefContigID, refEnd-refStart, lenSum/(refEnd-refStart), confSum/mapcnt, RefStartPos, RefEndPos, min(QryStartPos, QryEndPos),max(QryStartPos,QryEndPos),orientation,Confidence);
      
      fclose(fp);
      return 0;
    }
    
  }
  fclose(fp);

  if(mapcnt <= 0 || refStart >= refEnd)
    printf("%lld %0.1f 0.00 0.00 Unmapped\n", RefContigID, refEnd-refStart);
  else
    printf("%lld %0.1f %0.2f %0.2f Unmapped\n", RefContigID, refEnd-refStart, lenSum/(refEnd-refStart), confSum/mapcnt);

}
