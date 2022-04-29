#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

#include "globals.h"
#include "parameters.h"

#ifndef SCORE_H
#define SCORE_H

static Ident GentigRefScore_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/GentigRefScore.h 5409 2016-10-03 19:32:08Z tanantharaman $");

class Ypen {
 public:
  double VarPenalty;
  double RefBias;
} **Ypenalty = 0;

/* Gentig Guided Alignment scoring functions */

/* precompute values to speed up likelihood score computation */
static double LBias, FnPenalty, FpPenalty, MatchScore, Ylambda, MinScore;

static void score_init(register Cmap *refmap)
{
  if(refmap->numcolor != 1){
    printf("alignment score for multi-color data not yet implemented\n");
    exit(1);
  }
  register int N = refmap->numsite[0];
  register double *Y = refmap->site[0];/* Y[0..N+1] */
  
  register double pFN = FN[0];
  register double pTP = 1.0 - pFN;
  register double F = FP[0]*0.01;
  register double var = SD[0]*SD[0];

  register double lambda = Y[N+1]/((N>0) ? N : 1);
  register double theta = 1.0/(F + pTP/lambda) + var*0.5;
  register double Ygm = log(Y[1]+Y[N+1]-Y[N]);/* to compute Geometric Mean of Y intervals */
  for(register int i=1;i<N;i++)
    Ygm += log(Y[i+1]-Y[i]);
  Ygm /= ((N>0) ? N : 1);
  Ygm = exp(Ygm);

  register double Gcorr = 0.5*log(Ygm/lambda);

  if(SCORES)
    Ylambda = lambda;

  if(F*theta >= 1.0){
    printf("FP density=%0.4f, theta=%0.3f: FP is too large for meaningful alignment scoring\n",F,theta);
    exit(1);
  }

  register double FBias = Kbias * (-(pFN/pTP)*log(pFN) + 0.5*(1.0-log(1-F*theta)));
  register double Bias = F*(Kbias*log(pTP/(F*sqrt(2*M_PI*var*theta))) - FBias);

  if(Bias <= 0.0 || FBias <= 0.0){
    printf("FP density=%0.4f is too large for meaningful alignment scoring (Fbias=%e,Bias=%e)\n",
	   F,FBias,Bias);
    printf("Kbias=%0.2f,pFN=%0.4f,pTP=%0.4f,sd=%0.4f,F=%0.6f,theta=%0.3f,lambda=%0.3f\n",
	   Kbias,pFN,pTP,SD[0],F,theta,lambda);
    exit(1);
  }

  if(VERB){
    printf("Kbias=%0.2f,pTP=%0.4f,F=%0.6f,theta=%0.3f,lambda=%0.3f\n",
	   Kbias,pTP,F,theta,lambda);
    fflush(stdout);
  }

  /* allocate Ypenalty[1..DELTA_Y][0..N] */
  assert(Ypenalty==0);
  Ypenalty = new Ypen*[DELTA_Y+1];
  for(register int delta = 1; delta <= DELTA_Y; delta++)
    Ypenalty[delta] = new Ypen[N+1];

  LBias = Bias;/* X length bias */
  FnPenalty = -log(pFN);/* missing cut penalty*/
  FpPenalty = -log(F*sqrt(2*M_PI*var*theta)) - FBias;/* false cut penalty */
  MatchScore = FBias + log(pTP);

  if(GAMMA != 1.0){
    FnPenalty *= GAMMA;
    FpPenalty *= GAMMA;
  }

  register double Ivar = 0.5/var;
  if(VERB){
    printf("LBias=%0.6f,FnPenalty=%0.4f,FpPenalty=%0.4f,MatchScore=%0.3f,Ivar=%0.3f,Gcorr=%0.3f,minscore=%0.3f\n",
	   LBias,FnPenalty,FpPenalty,MatchScore,Ivar,Gcorr,-(DELTA_Y-1)*FnPenalty-(DELTA_X-1)*FpPenalty);
    fflush(stdout);
  }

  register double offset = MatchScore + Gcorr;

  for(register int delta = 1; delta <= DELTA_Y; delta++)
    for(register int I = 0; I <= N+1-delta; I++){/* precompute terms for reference interval Y[I..I+delta] */
      register Ypen *p = &Ypenalty[delta][I];
      register double y = Y[I+delta] - Y[I];
      register double yInv = 1.0/y;
      p->VarPenalty = Ivar*yInv;/* penalty factor for error squared for y = Y[I..I+delta] */
      p->RefBias = 0.5*log(theta*yInv) + offset;/* Bias factor for y = Y[I..I+delta] + Sm(J,I) constants */
    }

  MinScore = -(DELTA_Y-1)*FnPenalty-(DELTA_X-1)*FpPenalty;
}

static void score_free(register Cmap *refmap)
{
  if(Ypenalty){
    for(register int delta = 1; delta <= DELTA_Y;delta++)
      if(Ypenalty[delta])
	delete [] Ypenalty[delta];
    delete [] Ypenalty;
    Ypenalty = 0;
  }
}

/* includes Sm(J,I) */
static inline double Sint(register double X,
			  register double Y,
			  register int m,
			  register int n,
			  register int J,
			  register int I)
  
{
  register Ypen *p = &Ypenalty[n][I-n];/* penalty terms for Y[I-n .. I] */
  register double err = X-Y;

  /* MatchScore is included in RefBias */

  return X*LBias + p->RefBias - p->VarPenalty * err*err - (n-1)*FnPenalty - (m-1)*FpPenalty;
}

static inline double Send(register double X,
			  register int m,
			  register int n)
{
  return X*LBias - (n-1)*FnPenalty -(m-1)*FpPenalty;
}

static inline double Sm(register int J,register int I)
{
  return MatchScore;
}

#endif
