#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "globals.h"
#include "parameters.h"

#ifndef SCORE_H
#define SCORE_H

static Ident PGentigPairScore_h_Id("$Header: $");

/* Gentig Pairwise Alignment scoring functions */

/* precompute values to speed up likelihood score computation */
static double LBias, FBias, MissPenalty, MatchScore, Ylambda, Ytheta, MinScore, Ivar;

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

  if(F*theta >= 1.0){
    printf("FP density=%0.4f, theta=%0.6f: FP is too large for meaningful alignment scoring\n",F,theta);
    exit(1);
  }

  register double Pmiss = pTP*pFN/(pTP*pFN + F*lambda);

  MissPenalty = -Pmiss*log(pFN) + (1.0-Pmiss)*log(2*pTP/(F*sqrt(2*M_PI*var*theta)));

  if(GAMMA != 1.0)
    MissPenalty *= GAMMA;

  FBias = Kbias * ( pTP*0.5 + pFN*MissPenalty );
  LBias = F*Kbias*MissPenalty;

  if(VERB){
    printf("Kbias=%0.2f,Gbias=%0.2f,pTP=%0.6f,F=%0.6f,theta=%0.6f,lambda=%0.6f\n",
	   Kbias,Gbias,pTP,F,theta,lambda);
    fflush(stdout);
  }

  MatchScore = FBias;
  FBias *= 0.5;

  MinScore = -2*(DELTA-1)*MissPenalty;
  Ivar = 0.5/var;

  if(VERB){
    printf("LBias=%0.6f,FBias=%0.6f,MissPenalty=%0.6f,MatchScore=%0.6f,var=%0.4f,minscore=%0.3f\n",
	   LBias,FBias,MissPenalty,MatchScore,var,MinScore);
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
  register double err = X-Y;
  register double XYsum = X+Y;

  /*  printf("S(x=%0.3f,y=%0.3f,m=%d,n=%d)+S(%d,%d)=%0.6f\n",
      X,Y,m,n,J,I,MatchScore + XYsum*LBias - err*err*Ivar/XYsum +(FBias - MissPenalty)*(m+n-2.0));*/
  return MatchScore + XYsum*LBias - err*err*Ivar/XYsum +(FBias - MissPenalty)*(m+n-2.0);
}

/* short for Send() + S(X,Y,1,1) - S((X+Y)/2,(X+Y)/2,1,1) */
static inline double Sbnd(register double X,
			  register double Y,
			  register int m,
			  register int n)
{
  register double err = X-Y;
  register double XYsum = X+Y;

  return 2.0*X*LBias - err*err*Ivar/XYsum + (FBias - MissPenalty)*(m+n-2.0) - 0.25;
}

static inline double Send(register double X,
			  register int m,
			  register int n)
{
  /*  printf("Send(%0.3f,%d,%d)=%0.6f\n",X,m,n,
      2.0*X*LBias + (FBias - MissPenalty)*(m+n-2.0) - 0.25);*/
  return 2.0*X*LBias + (FBias - MissPenalty)*(m+n-2.0) - 0.25;
}

static inline double Sm(register int J,register int I)
{
  
  /*  printf("Sm(%d,%d)=%0.6f\n",J,I,MatchScore);*/
  return MatchScore;
}

#endif
