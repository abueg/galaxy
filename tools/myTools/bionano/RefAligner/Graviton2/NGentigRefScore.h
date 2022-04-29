#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "parameters.h"

#ifndef SCORE_H
#define SCORE_H

static Ident NGentigRefScore_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/NGentigRefScore.h 337 2012-04-10 23:41:03Z tanantharaman $");

#define RES_BND 1 /* Apply LBias only to interval size above resKB x 2 */

#define KLBIAS KLbias /* to scale LBias */
#define KFBIAS KFbias /* to scale EBias */
#define BETA Beta /* scale  log(G/py) term */

class Ypen {
 public:
  double VarPenalty;/* 1.0/(2*(SF[0]*SF[0]+SD[0]*SD[0]*y)) */
  double RefBias;/* 0.5*log(Ygm/(y*pTP)) */
} **Ypenalty = 0;

/* Gentig Guided Alignment scoring functions */

/* precompute values to speed up likelihood score computation */
static double LBias, EBias, FnPenalty, FpPenalty, Ylambda, MinScore, resKB2;
//static int Verb=0;

static void score_init(register Cmap *refmap)
{
  if(refmap->numcolor != 1){
    printf("alignment score for multi-color data not yet implemented\n");
    exit(1);
  }
  register int N = refmap->numsite[0];
  register double *Y = refmap->site[0];/* Y[0..N+1] */
  
  register double resKB = res * PixelLen;
  resKB2 = (RES_BND ? 2.0 * resKB : 0.0);

  register double pFN = FN[0];
  register double pTP = 1.0 - pFN;
  register double F = FP[0]*0.01;
  register double var = SD[0]*SD[0];
  register double varF = SF[0]*SF[0];

  register double lambda = Y[N+1]/((N>0) ? N : 1);
  register double theta = 1.0/(F + pTP/lambda) + var*0.5 + resKB2*0.5;
  register double Ygm = log(Y[1]+Y[N+1]-Y[N]);/* to compute Geometric Mean of Y intervals */
  for(register int i=1;i<N;i++)
    Ygm += log(Y[i+1]-Y[i]);
  Ygm /= ((N>0) ? N : 1);
  Ygm = exp(Ygm);

  Ylambda = lambda;

  assert(FN[0] > 0.0);
  assert(FN[0] < 1.0);
  assert(theta > 0.0);
  assert(Ylambda > 0.0);

  if(F*theta >= 1.0){
    printf("FP density=%0.4f, theta=%0.3f: FP is too large for meaningful alignment scoring\n",F,theta);
    exit(1);
  }

  LBias = -KLBIAS * F*log(F*(theta-resKB2*0.5));/* Bias per aligned kb */
  EBias = -0.5*KFBIAS*pFN*log(pFN)/pTP;/* Bias per aligned end */
  register double FBias = 0.5*Kbias + 2.0*EBias;/* Bias per aligned internal interval */

  if(VERB){
    printf("Kbias=%0.2f,KLbias=%0.2f,KFbias=%0.2f,G=%0.2f,B=%0.2f,pTP=%0.4f,F=%0.6f,theta=%0.3f,lambda=%0.3f,Ygm=%0.3f,resKB2=%0.3f\n",
	   Kbias,KLBIAS,KFBIAS,GAMMA,BETA,pTP,F,theta,lambda,Ygm,resKB2);
    fflush(stdout);
  }

  /* allocate Ypenalty[1..DELTA_YEND][0..N] */
  assert(Ypenalty==0);
  Ypenalty = new Ypen*[DELTA_YEND+1];
  for(register int delta = 1; delta <= DELTA_YEND; delta++)
    Ypenalty[delta] = new Ypen[N+1];

  FnPenalty = -log(pFN);/* missing cut penalty*/
  FpPenalty = -log(F*(theta - resKB2*0.5));/* false cut penalty */

  if(GAMMA != 1.0){
    FnPenalty *= GAMMA;
    FpPenalty *= GAMMA;
  }

  MinScore = -(DELTA_Y-1)*FnPenalty-(DELTA_X-1)*FpPenalty;
  MinScore *= 2.0;


  if(VERB){
    printf("LBias=%0.6f,EBias=%0.3f,FBias=%0.3f,FnPenalty=%0.4f,FpPenalty=%0.4f,sf=%0.4f,sd=%0.4f,minscore=%0.3f\n",
	   LBias,EBias,FBias,FnPenalty,FpPenalty,SF[0],SD[0],MinScore);
    fflush(stdout);
  }

  register double YgmBpTP = Ygm/pTP;
  for(register int delta = 1; delta <= DELTA_YEND; delta++){
    register int I;
    for(I = 0; I <= N+1-delta; I++){/* precompute terms for reference interval Y[I..I+delta] */
      register Ypen *p = &Ypenalty[delta][I];
      register double y = Y[I+delta] - Y[I];
      p->VarPenalty = 0.5/(varF+var*y);/* penalty factor for error squared for y = Y[I..I+delta] */
      p->RefBias = BETA * 0.5*log((varF + var*YgmBpTP)/(varF+var*y)) + FBias ;/* per aligned interval Bias factor for y = Y[I..I+delta] + Sm(J,I) constants */
      if(RES_BND)
	p->RefBias += max(0.0,y-resKB2)*LBias;
      else
	p->RefBias += y*LBias;
    }
    if(DEBUG)
      for(;I <= N;I++){
	register Ypen *p = &Ypenalty[delta][I];
#ifdef WIN32
	double denom = 0;
	p->VarPenalty = 0/denom;
	p->RefBias = 0/denom;
#else
	p->VarPenalty = nan("NaN");
	p->RefBias = nan("NaN");
#endif
      }
  }
}

static void score_free(register Cmap *refmap)
{
  if(Ypenalty){
    for(register int delta = 1; delta <= DELTA_YEND;delta++)
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
			  register int I)/* right end of interval! */
  
{
  assert(n>=1 && n <= DELTA_Y);
  assert(I>=n);
  register Ypen *p = &Ypenalty[n][I-n];/* penalty terms for Y[I-n .. I] */
  register double err = X-Y;

  /* LBias & FBias are included in p->RefBias */

  /*  if(Verb){
    printf("Sint(X=%0.3f,Y=%0.3f,m=%d,n=%d,J=%d,I=%d)\n",X,Y,m,n,J,I);
    printf("     VarPenalty=%0.4f,err=%0.3f(delta=%0.4f),RefBias=%0.4f\n",
	   p->VarPenalty,err,p->VarPenalty*err*err,p->RefBias);
	   }*/

  return p->RefBias - p->VarPenalty * err*err - (n-1)*FnPenalty - (m-1)*FpPenalty;
}

/* short for Send(Y,m,n) + S(X,Y,1,1) - S((X+Y)/2,(X+Y)/2,1,1) */
static inline double Sbnd(register double X,
			  register double Y,
			  register int m,
			  register int n,
			  register int I)/* right end of interval! */
{
  assert(n>=1 && n <= DELTA_YEND);
  assert(I>=n);
  register Ypen *p = &Ypenalty[n][I-n];/* penalty terms for Y[I-n .. I] */
  register double err = X-Y;

  if(RES_BND)
    return max(0.0,Y-resKB2)*LBias + EBias - (n-1)*FnPenalty - (m-1)*FpPenalty -p->VarPenalty*err*err;
  else
    return Y*LBias + EBias - (n-1)*FnPenalty - (m-1)*FpPenalty - p->VarPenalty*err*err;
}

static inline double Send(register double X,
			  register int m,
			  register int n)
{
  if(RES_BND)
    return max(0.0,X-resKB2)*LBias + EBias - (n-1)*FnPenalty -(m-1)*FpPenalty;
  else
    return X*LBias + EBias - (n-1)*FnPenalty -(m-1)*FpPenalty;
}

static inline double Sm(register int J,register int I)
{
  return 0.0;
}

#endif
