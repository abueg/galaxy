#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

#include "globals.h"
#include "parameters.h"

#ifndef SCORE_H
#define SCORE_H

static Ident GentigPairScore_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/GentigPairScore.h 5409 2016-10-03 19:32:08Z tanantharaman $");

/* Gentig Pairwise Alignment scoring functions */

#define RES_BND 0 /* Apply LBias only to interval size above resKB x 2 */

#define KLBIAS KLbias 
#define KFBIAS KFbias

/* precompute values to speed up likelihood score computation */
static double LBias, FBias, MissPenalty, MatchScore, Ylambda, Ytheta, MinScore, Ivar, resKB2;

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
    printf("FP density=%0.4f, theta=%0.3f: FP is too large for meaningful alignment scoring\n",F,theta);
    exit(1);
  }

  register double Pmiss = pTP*pFN/(pTP*pFN + F*lambda);

  MissPenalty = -Pmiss*log(pFN) + (1.0-Pmiss)*log(2*pTP*pTP/(F*sqrt(2*M_PI*var*lambda)));
  if(GAMMA != 1.0)
    MissPenalty *= GAMMA;

  //  register double FBias = Kbias * 2.0*pFN*MissPenalty/pTP;
  //  register double LBias = F*Kbias*MissPenalty;
  FBias = KFBIAS * 2.0*pFN*MissPenalty/pTP;
  LBias = KLBIAS * F * MissPenalty;

  if(VERB){
    printf("Kbias=%0.2f,KLbias=%0.2f,KFbias=%0.2f,Gbias=%0.2f,pTP=%0.6f,F=%0.6f,theta=%0.6f,lambda=%0.6f\n",
	   Kbias,KLbias,KFbias,Gbias,pTP,F,theta,lambda);
    fflush(stdout);
  }

  MatchScore = FBias + Kbias*0.5;

  MinScore = -2*(DELTA-1)*MissPenalty;
  MinScore *= 2.0;

  Ivar = 0.5/var;

  if(RES_BND)
    resKB2 = 2.0 * mres * PixelLen;

  if(VERB){
    printf("LBias=%0.6f,FBias=%0.4f,MissPenalty=%0.5f,MatchScore=%0.4f,var=%0.4f,minscore=%0.3f\n",
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

  if(RES_BND)
    return (max(0.0,X-resKB2)+max(0.0,Y-resKB2))*LBias + MatchScore - err*err*Ivar/XYsum - MissPenalty*(m+n-2.0);
  else
    return XYsum*LBias + MatchScore - err*err*Ivar/XYsum - MissPenalty*(m+n-2.0);
}

static inline double Sbnd(register double X,
			  register double Y,
			  register int m,
			  register int n)
{
  register double err = X-Y;
  register double XYsum = X+Y;

  if(RES_BND)
    return (max(0.0,X-resKB2)+max(0.0,Y-resKB2))*LBias + FBias*0.5 - err*err*Ivar/XYsum - MissPenalty*(m+n-2.0);
  else
    return XYsum*LBias + FBias*0.5 - err*err*Ivar/XYsum - MissPenalty*(m+n-2.0);
}

static inline double Send(register double X,
			  register int m,
			  register int n)
{
  if(RES_BND)
    return 2.0*max(0.0,X-resKB2)*LBias + FBias*0.5 - MissPenalty*(m+n-2.0);
  else
    return 2.0*X*LBias + FBias*0.5 - MissPenalty*(m+n-2.0);
}

static inline double Sm(register int J,register int I)
{
  return 0.0;
}

#endif
