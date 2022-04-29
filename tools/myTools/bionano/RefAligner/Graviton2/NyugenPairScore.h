#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "globals.h"
#include "parameters.h"

#ifndef SCORE_H
#define SCORE_H

static Ident NyugenPairScore_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/NyugenPairScore.h 5409 2016-10-03 19:32:08Z tanantharaman $");

/* Nguyen's Pairwise Alignment scoring functions */
/* The Alpha value is taken from the global Kbias value
   The Gamma value is taken from the glabal Gbias value
   The Beta value is a command line parameter */

#define ALPHA Kbias
#define BETA Beta
#define GAMMA Gbias

/* precomputed values used by pairalign.cpp */
static double MinScore,Ylambda,Ytheta;

/* precomputed values to speed up likelihood score computation */
static double Isd;/* used locally : 1.0/(SD[0]*sqrt(2.0)) */

static void score_init(int mapstart, int mapend)
{
  if(mapend < mapstart)
    return;

  if(map[mapstart]->numcolor != 1){
    printf("alignment score for multi-color data not yet implemented\n");
    exit(1);
  }

  register double pFN = FN[0];
  register double pTP = 1.0 - pFN;
  register double F = FP[0]*0.01;
  register double var = SD[0]*SD[0];

  /* compute theta, the average observed interval size */
  register double theta = 0.0;
  register int sitecnt = 0;
  for(register int i=mapstart; i <= mapend;i++){
    register Cmap *Xmap = map[i];
    register int M = Xmap->numsite[0];
    sitecnt += M;
    theta += Xmap->site[0][M+1];
  }
  if(sitecnt <= 0 || theta <= 0.0){
    fprintf(stderr,"score_init():Insufficient data:mapstart=%d,mapend=%d,sitecnt=%d\n",
	    mapstart,mapend,sitecnt);
    exit(1);
  }
  theta /= sitecnt;

  register double lambda = pTP/(1.0/(theta - var*0.5) - F);

  Ylambda = lambda;
  Ytheta = theta;

  if(VERB){
    printf("Alpha=%0.4f,Beta=%0.4f,Gamma=%0.4f,theta=%0.3f,lambda=%0.3f\n",
	   ALPHA,BETA,GAMMA,theta,lambda);
    fflush(stdout);
  }


  MinScore = -2*(DELTA-1)*GAMMA;
  Isd = sqrt(1.0/(2.0*var));

  if(VERB){
    printf("var=%0.4f,minscore=%0.3f\n", var,MinScore);
    fflush(stdout);
  }
}

static void score_free()
{
}

/* Sint() includes Sm(J,I) */
static inline double Sint(register double X,
			  register double Y,
			  register int m,
			  register int n,
			  register int J,
			  register int I)
  
{
  register double XYsum = X+Y;
  register double err = fabs(X-Y);
  register double d = err*Isd/sqrt(XYsum);

  return (d < ALPHA) ? (ALPHA-d)-GAMMA*(m+n-2.0) 
    : BETA*(ALPHA-d)-GAMMA*(m+n-2.0);
}

static inline double Sbnd(register double X,
			  register double Y,
			  register int m,
			  register int n)
{
  register double err = fabs(X-Y);
  register double XYsum = X+Y;
  register double d = err*Isd/sqrt(XYsum);

  return (d < ALPHA) ? -d - GAMMA*(m+n-2.0)
    : BETA*(ALPHA-d)-ALPHA - GAMMA*(m+n-2.0);
}

static inline double Send(register double X,
			  register int m,
			  register int n)
{
  return - GAMMA*(m+n-2.0);
}

static inline double Sm(register int J,register int I)
{
  return 0.0;
}

#endif
