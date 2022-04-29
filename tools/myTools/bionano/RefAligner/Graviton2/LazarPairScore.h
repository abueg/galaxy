#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "globals.h"
#include "parameters.h"

#ifndef SCORE_H
#define SCORE_H

static Ident LazarPairScore_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/LazarPairScore.h 5409 2016-10-03 19:32:08Z tanantharaman $");

/* John Lazar's refactored Pairwise Alignment scoring functions */
/* The Missaligned cut penalty is scaled by the global Gbias factor (normally Gbias=1) */

/* precomputed values used by pairalign.cpp */
static double MinScore,Ylambda,Ytheta;

/* precomputed values to speed up likelihood score computation */
static double MissPenalty, MatchScore, InvLambda,Ivar;

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

  MissPenalty = log((lambda/theta)/(pTP*pFN + F*lambda));
  if(GAMMA != 1.0)
    MissPenalty *= GAMMA;


  if(VERB){
    printf("Gbias=%0.2f,pTP=%0.6f,F=%0.6f,theta=%0.3f,lambda=%0.3f\n",
	   Gbias,pTP,F,theta,lambda);
    fflush(stdout);
  }

  InvLambda = 1.0/lambda;
  MatchScore = log(theta/sqrt(2.0*M_PI*var));
  Ivar = 0.5/var;
  MinScore = -2*(DELTA-1)*MissPenalty;

  if(VERB){
    printf("MissPenalty=%0.4f,MatchScore=%0.4f,var=%0.4f,minscore=%0.3f\n",
	   MissPenalty,MatchScore,var,MinScore);
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
  register double err = fabs(X-Y);
  register double XYsum = X+Y;

  return MatchScore + 0.5*log(XYsum) - err*err*Ivar/XYsum + err*InvLambda - MissPenalty*(m+n-2.0);
}

static inline double Sbnd(register double X,
			  register double Y,
			  register int m,
			  register int n)
{
  register double err = fabs(X-Y);
  register double XYsum = X+Y;

  return -err*err*Ivar/XYsum + err*InvLambda - MissPenalty*(m+n-2.0);
}

static inline double Send(register double X,
			  register int m,
			  register int n)
{
  return - MissPenalty*(m+n-2.0);
}

static inline double Sm(register int J,register int I)
{
  return 0.0;
}

#endif
