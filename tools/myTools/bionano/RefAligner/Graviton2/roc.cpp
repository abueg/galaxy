#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <ctype.h>

const char *SVN_ID =  (char *)"$Header: $";

char *usage = (char *)"roc input.align\n";

class Calign {
  int id1;
  int id2;
  double score;
  double offset;
  double overlap;
  int orientation;
  double logPV;
};

int main(int argc, char **argv)
{

  if(argc < 2){
    printf("usage:\n%s\n",usage);
    exit(1);
  }

  char *filename = argv[1];

  FILE *fp = fopen(filename,"r");
  if(fp==NULL){
    printf("Failed to read input file %s\n",filename);
    exit(1);
  }

  //  int maxalign = 1024;
  //  Calign *align = new Calign[maxalign];
  //  int numalign = 0;

  char buf[BUFSIZ];

  for(int linecnt=1;fgets(buf,BUFSIZ,fp)!=NULL;linecnt++){
    int len = strlen(buf);
    if(len >= BUFSIZ-1 || (buf[len-1] != '\n'  && buf[len-1] != '\r')){
      printf("Line too long or not terminated (%c) in input file %s on line %d:\n%s\n",
	     buf[len-1],filename,linecnt,buf);
      exit(1);
    }
    if(buf[0] ==  '#')
      continue;
    /* HERE */
  }

  fclose(fp);
  exit(0);
  return 0;
}
