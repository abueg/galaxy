#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

#include "globals.h"
#include "parameters.h"

#ifndef RGENTIGREFSCORE2_H
#define RGENTIGREFSCORE2_H

static Ident RGentigRefScore2_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/trunk/RGentigRefScore2.h 1126 2013-01-23 22:15:43Z tanantharaman $");

#define COLOR_SHIFT 0 /* 0 : assume color shift (at right end) depends on absolute distance between rightmost aligned sites of each color
			 1 : assume color shift (at right end) depends on Molecule location of rightmost site (approximate Molecule length) 
			 2 : assume color shift (at right end) depends on Reference location of rightmost site (approximate Reference length) (Required for REFDEBUG) */

#define OUTLIER_TYPE 0 /* OUTLIER type :
			  0 : Include both misaligned sites and sizing error in outlier (requires DELTA_X=DELTA_Y >= n to evaluate exactly)
			  1 : Include only sizing error in outlier 
			  NOTE : Only implemented for maptype=0 : If maptype==1, always behaves like OUTLIER_TYPE==0
		       */

#define OUTLIER_LTYPE 1 /* OUTLIER exponential decay type :
			   0 : outlier likelihood is proportional to exp(-(x+y)/outlierLambda)
			   1 : outlier likelihood is proportional to exp(-|x-y|/outlierLambda)
			   NOTE : Only implemented with maptype=1 : If maptype==0, always behaves like OUTLIER_LTYPE==1.
			*/

#define OUTLIER_DELTA(err) (!OUTLIER_MAX || fabs(err) < ((TFLOAT)outlierMax))

#define OUTLIER(maptype,x,y) (maptype || OUTLIER_DELTA(x-y)) /* only test alignments for outliers if they satisfy this condition : x = Nanomap, y = reference (interval size), maptype = 0(vfixx or bnx) or 1(Cmap) */

#define SCORE_Y 0 /* make score independent of X (less accurate Theta estimate and depends on Y, but result does not change when mapset is filtered) */


#define FAST_OUTLIER 1 /* WAS 0 */ /* approximate OutlierBias = log(1-Poutlier)  with 0.0, since Poutlier is typically 0.0003 is less */

#define KBIAS Kbias /* to scale Gaussian penalty bias */
#define KLBIAS KLbias /* to scale False cut penalty bias */
#define KFBIAS KFbias /* to scale Missing cut penalty bias */
#define BETA 1.0 /* Beta */ /* scale  log(G/py) term */

static inline double mylog(double x)
{
  return log(x);
}

static inline float mylog(float x)
{
  return logf(x);
}

static TFLOAT *OutlierPenaltyBC = NULL;

static class Cprtab {
 public:
  FLOAT LogPr;/* log(Pr(Y[c][I]-Y[c][I-K])) */
  TFLOAT Sm;/* Sm(0,I,K,Y[c],c) */
} ***PRtab[MAXCOLOR];/* PRtab[c][refid][K][I] has precomputed penalties that depend on refmap[refid]->site[c][I-K ... I] */

static class CYPen {
 public:
  TFLOAT deltaY;/**< Yc(YY[refid],I,K)-Yc(YY[refid],I-K-n,T) */
  TFLOAT Pen;/**< RefPen - 0.5*log(var) - (n-1)*FnPenalty + FpPenalty + PRtab[refid][n-1][I-K].LogPr + (FIX_OUTLIER_MIS ? 0.0 : Sm(0,I,K,Y)) */
  TFLOAT Sm;/**< Sm(0,I,K,Y) : only initialized if FIX_OUTLIER_MIS */
  TFLOAT Ivar;/**< 1.0/var where var = F2var + S2var*deltaY + R2var*deltaY*deltaY + E2var*(resL*resL+resR*resR) */
} *****YPen[MAXCOLOR];/**< YPen[c][refid][K=0..min(N-1,KMAX)][T=0..KMAX][I=K+1..N][n=1..DELTA_Y] */ 

class CMemCnt {
 public:
  long long YPmemPPP, YPmemPP, YPmemP, YPmem;
  long long PRmemP,PRmem;
};

static Cprtab **cprtab_p_pool[MAXCOLOR];
static Cprtab *cprtab_pool[MAXCOLOR];
static CYPen **** cypen_ppp_pool[MAXCOLOR];
static CYPen *** cypen_pp_pool[MAXCOLOR];
static CYPen **  cypen_p_pool[MAXCOLOR];
static CYPen *   cypen_pool[MAXCOLOR];
static size_t score_init_mem = 0;/* total memory that will be freed by score_free() */

class CXPen {
public:
  TFLOAT **deltaX[MAXCOLOR];/**< X[J] - X[J-m] */
  TFLOAT **Pen[MAXCOLOR];/**< m*FpPenalty + deltaX*Frate */
  TFLOAT **OutlierBias[MAXCOLOR]; /**< OutlierBias + deltaX*FrateBias + m*FBiasWT */
  TFLOAT **PenBias[MAXCOLOR]; /**< deltaX*FrateBias + m*FBiasWT - Pen */
  char *mem_pool;
  size_t mem_pool_size;
  CXPen() { 
    mem_pool = NULL;
    mem_pool_size = 0;
    if(DEBUG>=2)
      for(int c = 0; c < colors; c++)
	deltaX[c] = Pen[c] = OutlierBias[c] = PenBias[c] = NULL;  
  };
} ; /* XPen[k=0..nummaps-1].<field>[m=1..DELTA_X][J=m+1..M+1] */

#if 0
class CXPen {/* HERE : convert to Flat arrays as in RGentigRefScore.h */
 public:
  TFLOAT deltaX;/* X[J] - X[J-m] */
  TFLOAT Pen;/* m*FpPenalty[c] + (X[c][J] - X[c][J-m])*Frate[c] */
  TFLOAT OutlierBias; /* OutlierBias + (X[c][J]-X[c][J-m])*FrateBias[c] + m*FBiasWT[c] */
  TFLOAT PenBias; /* (X[c][J]-X[c][J-m])*FrateBias[c] + m*FBiasWT[c] - Pen */
}; /* XPen[k=0..nummaps-1][c][J=1..M+1][m=1..min(J,DELTA_X)] */
#endif

static long long maxN[MAXCOLOR];
static long long maxM[MAXCOLOR];

/* Gentig Guided Alignment scoring functions */

/* precompute values to speed up likelihood score computation */
static TFLOAT LBias[MAXCOLOR], FBias[MAXCOLOR], EBias[MAXCOLOR], PRbias[MAXCOLOR] = {0.0,0.0};
static TFLOAT FnPenalty[MAXCOLOR], FpPenalty[MAXCOLOR], MinScore[MAXCOLOR];
static TFLOAT OutlierPenalty[MAXCOLOR], OutlierBias, OutlierPenalty2[MAXCOLOR], OutlierEndPenalty, OutlierEndBias, OutlierLambdaInv, FpPenaltyBias[MAXCOLOR], Frate[MAXCOLOR], pTP[MAXCOLOR], log_pTP[MAXCOLOR];
static TFLOAT FrateBias[MAXCOLOR];/* biasWT * (Frate + LBias) */
static TFLOAT FBiasWT[MAXCOLOR];/* biasWT * FBIAS */
static TFLOAT EBiasWT[MAXCOLOR];/* biasWT * EBIAS */
static TFLOAT F2var[MAXCOLOR];/* 2*SF[c]*SF[c] */
static TFLOAT S2var[MAXCOLOR];/* 2*fabs(SD[c])*SD[c] */
static TFLOAT R2var[MAXCOLOR]; /* 2*SR[c]*SR[c] */
static TFLOAT E2var[MAXCOLOR]; /* 2*SE[c]*SE[c] */
static TFLOAT RefPen[MAXCOLOR];/* log(Xtheta[c]) - 0.5*log(M_PI) */

static double IresSD[MAXCOLOR];/* 1.0/(resSD * PixelLen * sqrt(2.0)) */
static double pTPd[MAXCOLOR], log_pTPd[MAXCOLOR], log_qTPd[MAXCOLOR];
static double Ylambda[MAXCOLOR], Xlambda[MAXCOLOR], Xtheta[MAXCOLOR] = {0.0,0.0};
static double resKB[MAXCOLOR], resKB2[MAXCOLOR], Ygm[MAXCOLOR], GsitesPerMap[MAXCOLOR];

static double Ylambda12;/* (Ylambda[0]*Ylambda[0]+Ylambda[1]*Ylambda[1])/(Ylambda[0]+Ylambda[1]) */
static double Xtheta12;/* (Xtheta[0]*Xtheta[0]+Xtheta[1]*Xtheta[1])/(Xtheta[0]+Xtheta[1]) */

static TFLOAT FM2var; /* 2 * SF[2] * SF[2] */
static TFLOAT SM2var; /* 2 * fabs(SD[2]) * SD[2] */
static TFLOAT RefPen12;/* log(Xtheta12) - 0.5*log(M_PI) */

static inline double Pr(double y, int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  double pr = 0.5*(1.0+erf((y - resKB[c])*IresSD[c]));
  return max(1.0e-100, pr);
}

/* Pn(y) == 1-Pr(y) */
static inline double Pn(double y, int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  return 0.5*erfc((y - resKB[c])*IresSD[c]);
}

static inline double log_Pn(double y, int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  return log(Pn(y,c));
}

/* NOTE : return value can be -Inf if Y[I] - Y[I-K] is large */
static inline double Sm(int J, int I, int K, double *Y, int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);

  if(SCORE_APPROX>=2){
    if(!K)
      return log_pTPd[c];
    if(K==1)
      return 2.0 * log_pTPd[c] + log_Pn(Y[I]-Y[I-K],c);
    return 2.0 * log_pTPd[c] + log_Pn(Y[I]-Y[I-K],c) + (K-1.0) * log_qTPd[c];
  } else if(MIS_VITERBI){ /* score = (K+1) * log_pTPd[c] + sum(t = 0 .. K-1) Log_Pn(Y[I-t]-Y[I-t-1],c) */
    double ret = (K+1.0) * log_pTPd[c];
    for(int t = 0; t < K; t++)
      ret += log_Pn(Y[I-t]-Y[I-t-1],c);
    return ret;

  } else {/* The following evaluation is NOT consistent with Viterbi scoring used elsewhere and will give higher score to alignments with misaligned sites, since
	     the Likelihood for all possible ways of resolving misaligned sites are being added up. To avoid this bias use the Viterbi score (see above) */
    double p = pTP[c];
#if RESSD==0 // NEW CODE : only K(K+1)/2 calls to Pn() instead of (K+1)/2 * 2^(K-1) calls AND always exact
    if(K <= 0)
      return log_pTPd[c];

    double *PMK = (double *)alloca(KMAX * sizeof(double));
    double *PowQ = (double *)alloca((KMAX+1) * sizeof(double));
    if(DEBUG>=2) assert(1 <= K && K <= KMAX);
    double q = 1.0 - p;

    PowQ[0] = 1.0;
    for(int k = 0; k < K; k++){
      PowQ[k+1] = q * PowQ[k];      /* PowQ[t = 0..k+1] == pow(q,t) */
      double Yk = Y[I-k-1];
      double sum = Pn(Y[I]-Yk, c) * PowQ[k];
      for(int t = 1; t <= k; t++)
	sum += Pn(Y[I-t] - Yk, c) * PowQ[k-t] * p * PMK[t-1];
      PMK[k] = sum;
    }
    return 2.0 * log_pTPd[c] + log(PMK[K-1]);

#else // OLD CODE

    switch(K){
    case 0:
      return log(p);
    case 1:
      return log(p*p * Pn(Y[I]-Y[I-1],c));
    case 2:
      if(RESSD)
	return log(p*p * p*Pn(Y[I]-Y[I-1],c)*Pn(Y[I-1]-Y[I-2],c));
      else
	return log(p*p * (p*Pn(Y[I]-Y[I-1],c)*Pn(Y[I-1]-Y[I-2],c) + 
			  (1.0-p)*Pn(Y[I]-Y[I-2],c)));
    case 3:
      if(RESSD)
	return log(p*p * p*p * Pn(Y[I]-Y[I-1],c)*Pn(Y[I-1]-Y[I-2],c)*Pn(Y[I-2]-Y[I-3],c)); 
      else
	return log(p*p * (p*p * Pn(Y[I]-Y[I-1],c)*Pn(Y[I-1]-Y[I-2],c)*Pn(Y[I-2]-Y[I-3],c) + 
			  p*(1.0-p) * (Pn(Y[I]-Y[I-2],c)*Pn(Y[I-2]-Y[I-3],c)+Pn(Y[I]-Y[I-1],c)*Pn(Y[I-1]-Y[I-3],c)) +
			  (1.0-p)*(1.0-p) * Pn(Y[I]-Y[I-3],c)));    
    default:  /* general case K>3 is approximate */
      double ret = 0;
      if(RESSD)
	ret = p;
      else {
	for(int j=1; j < K; j++){
	  double den = Pn(Y[I-j+1]-Y[I-j],c)*Pn(Y[I-j]-Y[I-j-1],c);
	  //	if(den > 0.0)
	  ret += Pn(Y[I-j+1]-Y[I-j-1],c)/den;
	}
	ret = p + (1.0 - p)*ret;
      }
      for(int i = 1; i <= K; i++)
	ret *= p*Pn(Y[I-i+1]-Y[I-i],c);
      //    if(ret <= 0.0)
      //      return -1.0e+300;
      return log(ret);
    }
#endif // OLD CODE
  }
}

static void score_free(int refstart, int refend)
{
  if(VERB>=2){
    printf("score_free():refstart=%d,refend=%d\n",refstart,refend);
    fflush(stdout);
  }
  for(int c = 0; c < colors; c++){
    if(PRtab[c]){
      delete [] cprtab_p_pool[c]; cprtab_p_pool[c] = NULL;
      delete [] cprtab_pool[c]; cprtab_pool[c] = NULL;
      delete [] PRtab[c]; PRtab[c] = 0;
    }
    if(YPen[c]){
      delete [] cypen_ppp_pool[c]; cypen_ppp_pool[c] = NULL;
      delete [] cypen_pp_pool[c]; cypen_pp_pool[c] = NULL;
      delete [] cypen_p_pool[c]; cypen_p_pool[c] = NULL;
      delete [] cypen_pool[c]; cypen_pool[c] = NULL;
      delete [] YPen[c]; YPen[c] = 0;
    }
  }
  delete [] OutlierPenaltyBC; OutlierPenaltyBC = NULL;
  score_init_mem = 0;
}

/* return value is 1, if Y[c][] has no intervals above minKB */

static int score_init(Cmap **refmap,int refstart, int refend, int nummaps, int ***RKKmax)
{
  score_free(refstart,refend);

  if(OUTLIER_TYPE != 0){
    printf("RGentigRefScore2.h:OUTLIER_TYPE==1 not yet implemented\n");
    exit(1);
  }
  if(OUTLIER_LTYPE==0 && maptype==0){
    printf("RGentigRefScore2.h: OUTLIER_LTYPE==0 not implemented with maptype==0\n");
    exit(1);
  }
  if(OUTLIER_TYPE==1 && maptype==1){// HERE HERE : remove after implementing for maptype==1
    printf("RGentigRefScore2.h: OUTLIER_TYPE==1 not implemented with maptype==1: Assuming OUTLIER_TYPE==0\n");
    fflush(stdout);
  }
  if(FAST_OUTLIER != 1){
    printf("RGentigRefScore2.h: FAST_OUTLIER != 1 no longer supported\n");
    exit(1);
  }
  if(biasWToutlier != 1.0){
    printf("biaswtOutlier != 1 : not yet supported with colors >= 2\n");
    exit(1);
  }
  if(biasWTend != 1.0){
    printf("biaswtEnd != 1 : not yet supported with colors >= 2\n");
    exit(1);
  }
  if(OUTLIER_TYPE != 0){
    printf("RGentigRefScore2.h: OUTLIER_TYPE==1 not yet supported\n");
    exit(1);
  }

  if(Beta != 1.0){
    printf("Beta != 1.0 no longer supported\n");
    exit(1);
  }

  if(DEBUG)assert(colors == 2);

  for(int c = 0; c < colors; c++){
    resKB[c] = res[c] * PixelLen;
    resKB2[c] = 0.0;
  }

  if(VERB){
    printf("K=%0.2f,KL=%0.2f,KF=%0.2f,G=%0.2f,B=%0.2f\n", KBIAS,KLBIAS,KFBIAS,GAMMA,BETA);
    fflush(stdout);
  }

  int numthreads = 1;
  #ifdef _OPENMP
  numthreads = omp_get_max_threads();
  #endif
  numthreads = min(numthreads,MaxThreads);

  for(int c = 0; c < colors; c++){
    double pFN = FN[c];
    pTP[c] = pTPd[c] = 1.0 - pFN;
    log_pTP[c] = log_pTPd[c] = log(pTP[c]);
    log_qTPd[c] = log(pFN);
    double F = FP[c]*0.01;
    double var = fabs(SD[c])*SD[c];
    double varF = SF[c]*SF[c];
    double varR = SR[c]*SR[c];
    double varE = SE[c]*SE[c];

    Ygm[c] = 0.0;/* to compute Geometric Mean of Y intervals */
    int Ycnt = 0;
    double Ysum = 0.0;

    for(int refid = refstart; refid < refend; refid++){
      Cmap *Ymap = refmap[refid];
      int N = Ymap->numsite[c];
      FLOAT *Y = Ymap->site[c];/* Y[0..N+1] */
      for(int i=1;i<N;i++){
	double y = Y[i+1]-Y[i];
	if(y > resKB[c]){
	  Ygm[c] += log(y);
	  Ysum += y;
	  Ycnt++;
	}
      }
    }
    if(!Ycnt || Ysum <= 0.0){
      printf("color=%d:No intervals over resKB=%0.3f in any of %d..%d reference maps\n",c+1,resKB[c],refstart,refend-1);
      exit(1);
    }

    Ygm[c] /= max(1,Ycnt);
    Ygm[c] = exp(Ygm[c]);

    Ylambda[c] = Ysum / max(1,Ycnt);

    /* compute theta based on X */
    double theta = 0.0;
    int sitecnt = 0;
    maxM[c] = 0;
    for(int i= 0; i < nummaps;i++){
      Cmap *Xmap = Gmap[i];
      int M = Xmap->numsite[c];
      sitecnt += M;
      theta += Xmap->site[c][M+1];
      if(M > maxM[c])
	maxM[c] = M;
    }
    if(DEBUG) assert(0 <= maxM[c] && maxM[c] <= MAXINT);
    
    int Rnummaps = nummaps;
    if(XmapCount[c] > 0){/* use command line values for XmapCount,XmapSites,XmapLength (adjusted for any change in PixelLen) */
      if(VERB>=2){
	printf("Using XmapCount[c]=%d,XmapSites=%lld,XmapLength=%0.8f,XmapGtheta=%0.8f\n",XmapCount[c],XmapSites[c],XmapLength[c],XmapGtheta[c]);
	fflush(stdout);
      }
      Rnummaps = XmapCount[c];
      sitecnt = XmapSites[c];
      theta = XmapLength[c] * PixelLen / origPixelLen;
    }

    if(sitecnt <= 0 || theta <= 0.0){
      fprintf(stderr,"score_init():Insufficient data:nummaps=%d,color=%d,total sites=%d,total length=%0.3f\n",
	      nummaps,c,sitecnt,theta);
      exit(1);
    }
    Xtheta[c] = theta/((double)sitecnt) - resKB2[c] * 0.5;
    double R = RTHETA_FIX ? Rnummaps/theta : 0.0;/* should be reduced by inverse of average linear chromosome size */
    GsitesPerMap[c] = ((double)sitecnt)/Rnummaps;

    if(F*Xtheta[c] >= 1.0){
      printf("WARNING:c=%d:FP=%0.4f (per 100 kb), average site interval Xtheta=%0.3f (kb): FP density exceeds total site density\n",c,F*100.0,Xtheta[c]);
      printf(" Reducing FP[c] to %0.4f (per 100kb)\n", 50.0/Xtheta[c]);
      F = 0.5/Xtheta[c];
      FP[c] = F*100.0;
    }

    Xlambda[c] = pTP[c]/(1.0/(Xtheta[c] - var*0.5) - F);/* HERE : does not include correction for SR[],SE[] */
    if(DEBUG && !SCORE_Y && !(Xlambda[c] > 0.0)){
      printf("c=%d:Xlambda[c]=%0.6e : pTP=%0.6e,Xtheta[c]=%0.6e (theta=%0.6e,resKB2=%0.6e,sitecnt=%d),var=%0.6e,F=%0.6e\n",
	     c,Xlambda[c], pTP[c], Xtheta[c], theta/sitecnt, resKB2[c], sitecnt, var, F);
      fflush(stdout);
      assert(Xlambda[c] > 0.0);
    }

    if(XmapCount[c] && RefRepeats==1){
      Ylambda[c] = Xlambda[c];
      Ygm[c] = XmapGtheta[c] * Xlambda[c] / Xtheta[c];
    } else if(SCORE_Y){ /* make score dependent only on Y (not X) */
      Xlambda[c] = Ylambda[c];
      Xtheta[c] = 1.0/(F + pTP[c]/Xlambda[c]) + var*0.5;
    }

    if(DEBUG){
      assert(FN[c] > 0.0);
      assert(FN[c] < 1.0);
      assert(theta > 0.0);
      assert(Ylambda[c] > 0.0);
      assert(Xlambda[c] > 0.0);
      assert(Xtheta[c] > 0.0);
    }

    Frate[c] = F + R - thetaScale/Xtheta[c];// WAS F - 1.0/Xtheta[c]

    double YgmBpTP = (Ygm[c]/pTP[c]) * (Xlambda[c]/Ylambda[c]);
    double Sbias = log(Xtheta[c]) - 0.5*log(2.0*M_PI*(varF+var*YgmBpTP+varR*YgmBpTP*YgmBpTP));
    FBias[c] = 0.5*KBIAS - KFBIAS*pFN*log(pFN)/pTP[c] - log(pTP[c]) - Sbias; /* Bias per site interval in X (overlapping Y)  = Abias - Sbias in document */
    LBias[c] = -KLBIAS * F*log(F*Xtheta[c]) - F*FBias[c];/* Bias per kb in X (overlapping Y) : In document LBias is larger by Frate = F - 1/Xtheta */
    EBias[c] = 0.5*(Sbias + FBias[c]);/* Bias per aligned end , corresponding to Abias/2 in the document. This is larger by KBIAS/4 compared to the document, so the final
					 score will be KBIAS/2 larger than in the document */
    if(resSD[c] > 0.0) /* PRbias =  Bias per site interval in X for log(Pr(y)) and Sm()-log(pTP) : will be estimated based on average of actual alignments in refalign */
      FBias[c] += PRbias[c];

    /* Note that TotalBias(X) = sum(c) (Len[c] * (Frate[c]+LBias[c]) + (M[c]-1)*FBias[c] + 2*EBias[c]) = sum(c) (Len[c] * (Frate[c]+LBias[c]) + M[c]*FBias[c] + Sbias) (scaled by biasWT)
       where Len[c] is the sum of the intervals of X[c], and M[c] is the number of sites in X[c]. This is larger by KBIAS/2 than in the document due to EBias (see above).

       The score of Nanomap X can be adjusted to determine the true -log(LR) :
       -log(LR(X)) = score(X) - TotalBias(X) + (1.0-BETA)*Sum(c,i) log(Xtheta[c]/sqrt(2*M_PI*(varF[c]+var[c]*y[c][i]))),
       where y[c][i] are the reference intervals aligned with X[c].
    */

    if(VERB){
      printf("   c=%d:pTP=%0.4f,F=%0.6f,Xtheta=%0.6f,Ylambda=%0.6f(Xl=%0.6f),Ygm=%0.6f,Frate=%0.8f,res=%0.6f,resSD=%0.6f,PixelLen= %0.10f\n",
	     c,pTP[c],F,Xtheta[c],Ylambda[c],Xlambda[c],Ygm[c], Frate[c], res[c], resSD[c], PixelLen);
      fflush(stdout);
    }

    FnPenalty[c] = -log(pFN);/* missing cut penalty*/
    FpPenalty[c] = -log(F*Xtheta[c]);/* false cut penalty */

    if(GAMMA != 1.0){
      FnPenalty[c] *= GAMMA;
      FpPenalty[c] *= GAMMA;
    }

    MinScore[c] = -(DELTA_Y-1)*FnPenalty[c]-(DELTA_X-1)*FpPenalty[c];
    MinScore[c] *= 2.0;

    if(Poutlier <= 0.0){
      OutlierPenalty[c] = LARGE_NEGATIVE;
      OutlierBias = 0.0;
    } else {
      if(maptype==0 && !OUTLIER_LTYPE){
	printf("WARNING: maptype=0,OUTLIER_LTYPE=0: this combination not implemented (treating as maptype=0, OUTLIER_LTYPE=1)\n");
	fflush(stdout);
      }
      if(maptype){
	if(outlierLambda >= 1000.0){
	  outlierLambda = 1e+30;
	  OutlierPenalty[c] = log(Poutlier);/* Outlier Likelihood = 1 (NOTE : this has an unbounded integral!) */
	} else {
	  if(OUTLIER_LTYPE == 0){
	    OutlierPenalty[c] = log(Poutlier * Xtheta[c] / outlierLambda);/* Outlier Likelihood = exp(-(x+y)/outlierLambda */
	  } else {
	    OutlierPenalty[c] = log(Poutlier * Xtheta[c] / (2.0 * outlierLambda));/* Outlier Likelihood = exp(-|x-y|/outlierLambda */
	  }
	}
      } else {/* maptype == 0 */
	if(outlierMax >= 1000.0){
	  outlierMax = 1e+30;
	  if(outlierLambda >= 1000.0)
	    OutlierPenalty[c] = log(Poutlier);/* Outlier Likelihood = 1 (NOTE : this has an unbounded integral!) */
	  else
	    OutlierPenalty[c] = log(Poutlier * Xtheta[c] / (2.0 * outlierLambda)); /* Outlier likelihood = exp(-|x-y|/outlierLambda) */
	} else {
	  if(outlierLambda >= 1000.0){
	    outlierLambda = 1e+30;
	    OutlierPenalty[c] = log(Poutlier * Xtheta[c] / (2.0 * outlierMax));/* Outlier Likelihood = 1 IFF |x-y| < outlierMax (0 otherwise) */
	  } else
	    OutlierPenalty[c] = log(Poutlier * Xtheta[c] / ( 2.0 * outlierLambda * (1.0 - exp(-outlierMax/outlierLambda))));/* Outlier Likelihood = exp(-|x-y|/outlierLambda) IFF |x-y| < outlierMax (0 otherwise) */
	}
      }
      OutlierPenalty[c] -= log(max(0.0001,1.0-Poutlier));
      OutlierBias = FAST_OUTLIER ? 0.0 : log(max(0.0001,1.0-Poutlier));
    }
    OutlierLambdaInv = 1.0/outlierLambda;
    OutlierPenalty2[c] = OutlierPenalty[c] + OutlierBias;

    if(PoutlierEnd <= 0.0){
      OutlierEndPenalty = LARGE_NEGATIVE;
      OutlierEndBias = 0.0;
    } else {
      OutlierEndPenalty = log(PoutlierEnd/max(0.0001,1.0-PoutlierEnd));
      OutlierEndBias = log(max(0.0001,1.0-PoutlierEnd));
    }

    if(VERB && (Poutlier > 0.0 || PoutlierEnd > 0.0)){
      printf("Poutlier=%0.6f(Bias=%0.6f,Pen=%0.6e),PoutlierEnd=%0.6f(EndBias=%0.6f,EndPen=%0.6e),Ch=%0.6f,ChFdr=%0.6f\n",
	     Poutlier,OutlierBias,OutlierPenalty[c],PoutlierEnd,OutlierEndBias,OutlierEndPenalty,Ch,ChFdr);
      fflush(stdout);
    }

    if(VERB){
      if(c <= 0)
	printf("       LBias=%0.6f,EBias=%0.4f,FBias=%0.4f(PRbias=%0.4f,Sbias=%0.4f)FnPenalty=%0.6f,FpPenalty=%0.6f,sf=%0.6f,sd=%0.6f,sr=%0.6f,se=%0.6f,minscore=%0.3f\n",
	       LBias[c],EBias[c],FBias[c],PRbias[c],Sbias,FnPenalty[c],FpPenalty[c],SF[c],SD[c],SR[c],SE[c],MinScore[c]);
      else
	printf("       LBias=%0.6f,EBias=%0.4f,FBias=%0.4f(PRbias=%0.4f,Sbias=%0.4f)FnPenalty=%0.6f,FpPenalty=%0.6f,sf=%0.6f,sd=%0.6f,sr=%0.6f,se=%0.6f,minscore=%0.3f,MA:sf[2]=%0.6f,sd[2]=%0.6f,MA_mean=%0.6f\n",
	       LBias[c],EBias[c],FBias[c],PRbias[c],Sbias,FnPenalty[c],FpPenalty[c],SF[c],SD[c],SR[c],SE[c],MinScore[c], SF[2],SD[2],MA_mean);
      fflush(stdout);
    }
    if(DEBUG) assert(isfinite(LBias[c]));
    if(DEBUG) assert(isfinite(EBias[c]));
    if(DEBUG) assert(isfinite(FBias[c]));
    if(DEBUG) assert(isfinite(Sbias));

    //    register int N = refmap->numsite[c];
    //    register FLOAT *Y = refmap->site[c];/* Y[0..N+1] */

    FpPenaltyBias[c] = FpPenalty[c] - FBias[c];
    FrateBias[c] = biasWT * (Frate[c] + LBias[c]);
    FBiasWT[c] = biasWT * FBias[c];
    EBiasWT[c] = biasWT * EBias[c];

    F2var[c] = 2.0 * varF;
    S2var[c] = 2.0 * var;
    R2var[c] = 2.0 * varR;
    E2var[c] = 2.0 * varE;

    RefPen[c] = log(Xtheta[c]) - 0.5*log(M_PI);

    if(DEBUG) assert(resSD[c] > 0.0);
    IresSD[c] = 1.0/(resSD[c] * PixelLen * sqrt(2.0));

    /* allocate PRtab[c][refid= refstart .. refend-1][K= 0..max(max(2,KMAX),DELTA_Y)][I= K+1..N]
       YPen[c][refid= refstart .. refend-1][K= 0..min(N[c]-1,KMAX)][T= 0..KMAX][I= K+1..N[c]][n= 1..DELTA_Y] */
    PRtab[c] = new Cprtab**[refend];
    YPen[c] = new CYPen****[refend];
    score_init_mem += refend*(sizeof(Cprtab**)+sizeof(CYPen****));
    
    int PRsiz = max(max(2,KMAX),DELTA_Y);

    /* First compute total memory required */
    long long YPmemPPP = 0, YPmemPP = 0, YPmemP = 0, YPmem = 0;/* size of total memory of type CYPen***,CYPen**,CYPen*,CYPen respectively */
    long long PRmemP = 0, PRmem = 0;/* size of total memory of type Cprtab*,Cprtab respectively */
    CMemCnt *MemCnt = new CMemCnt[refend-refstart+1];
    MemCnt -= refstart;/* MemCnt[i=refstart..refend-1] will have start location in memory pool(s) for reference contig i */

    for(int refid = refstart; refid < refend; refid++){
      Cmap *Ymap = refmap[refid];
      int N = Ymap->numsite[c];
      int Kmax = min(N-1,KMAX);

      CMemCnt *p = &MemCnt[refid];
      p->PRmemP = PRmemP;
      p->PRmem = PRmem;
      p->YPmemPPP = YPmemPPP;
      p->YPmemPP = YPmemPP;
      p->YPmemP = YPmemP;
      p->YPmem = YPmem;

      PRmemP += PRsiz+1;
      PRmem += (PRsiz+1)*(N+1);

      YPmemPPP += Kmax+1;
      int KK = (KMAX+1)*(Kmax+1);
      YPmemPP += KK;
      YPmemP += KK*(N+1);
      YPmem += KK*N*DELTA_Y;
    }
    if(DEBUG) assert(PRmem >= 0);
    if(DEBUG) assert(YPmem >= 0);
    
    long long memsiz = PRmemP*sizeof(Cprtab*) + PRmem*sizeof(Cprtab);
    if(VERB>=2){
      printf("c=%d:PRmemP=%lld,PRmem=%lld, total memory needed = %lld bytes\n",c,PRmemP,PRmem,memsiz);
      fflush(stdout);
    }
    cprtab_p_pool[c] = new Cprtab*[PRmemP];
    cprtab_pool[c] = new Cprtab[PRmem];
    score_init_mem += memsiz;

    memsiz = YPmemPPP*sizeof(CYPen***)+YPmemPP*sizeof(CYPen**) + YPmemP*sizeof(CYPen*)+YPmem*sizeof(CYPen);
    if(VERB>=2){
      printf("c=%d:YPmemPPP=%lld,YPmemPP=%lld,YPmemP=%lld,YPmem=%lld, total memory needed = %lld bytes\n",
	     c,YPmemPPP,YPmemPP,YPmemP,YPmem,memsiz);
      fflush(stdout);
    }
    cypen_ppp_pool[c] = new CYPen***[YPmemPPP];
    cypen_pp_pool[c] = new CYPen**[YPmemPP];
    cypen_p_pool[c] = new CYPen*[YPmemP];
    cypen_pool[c] = new CYPen[YPmem];
    score_init_mem += memsiz;

    if(VERB>=2){
      printf("c=%d:score_init:%lu bytes succesfully allocated\n",c,score_init_mem);
      fflush(stdout);
    }
    TFLOAT aNaN = nan("NaN");

    if(DEBUG>=2){

      #pragma omp parallel for num_threads(numthreads) schedule(static,256)
      for(long long i = 0; i < PRmem; i++){
	Cprtab *p = &cprtab_pool[c][i];
	p->LogPr = aNaN;
	p->Sm = aNaN;
      }

      #pragma omp parallel for num_threads(numthreads) schedule(static,256)
      for(long long i = 0; i < YPmem; i++){
	CYPen *p = &cypen_pool[c][i];
	p->Sm = aNaN;
	p->Pen = aNaN;
	p->deltaY = aNaN;
	p->Ivar = aNaN;
      }
    }

    int nthreads = min(refend-refstart,numthreads);

    #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1) if(nthreads > 1)
    for(int refid = refstart; refid < refend; refid++){
      Cmap *Ymap = refmap[refid];
      int N = Ymap->numsite[c];
      FLOAT *Y = Ymap->site[c];
      int *Kmax = RKKmax[refid][c];

      CMemCnt *pm = &MemCnt[refid];

      /* allocate PRtab[c][refid][K=0..PRsiz][I=0..N] */
      long long PRcnt = pm->PRmem;
      PRtab[c][refid] = &cprtab_p_pool[c][pm->PRmemP];

      for(int K = 0; K <= PRsiz; K++){
	int n = K+1;
	PRtab[c][refid][K] = &cprtab_pool[c][PRcnt]; PRcnt += N+1;
	Cprtab *pp = &PRtab[c][refid][K][0];

	if(n <= DELTA_Y){
          //      #pragma omp parallel for num_threads(numthreads) schedule(static,16)
	  for(int I = K+1; I <= N;I++){
	    Cprtab *p = &pp[I];
	    FLOAT y = Y[I] - Y[I-n];
	
	    if(DEBUG>=2 && !(isnan(p->LogPr))){
	      #pragma omp critical 
	      {
		printf("c=%d,refid=%d,K=%d,n=%d,I=%d(N=%d):p->LogPr=%0.8e, PRtab[c][refid][K][I].LogPr=%0.8e\n",
		       c,refid,K,n,I,N,p->LogPr, PRtab[c][refid][K][I].LogPr);
		fflush(stdout);
		assert(isnan(p->LogPr));
	      }
	    }
	    p->LogPr = log(Pr(y,c));
	    if(DEBUG) assert(isfinite(p->LogPr));/* could fail if Pr(y) == 0 */
	  }
	}

	if(K <= KMAX){
          //      #pragma omp parallel for num_threads(numthreads) schedule(static,16)
	  for(int I = K+1; I <= N;I++){
	    if(K <= Kmax[I]){
	      Cprtab *p = &pp[I];
			    
	      if(DEBUG>=2) assert(isnan(p->Sm));
	      p->Sm = Sm(0,I,K,Y,c);
	      if(DEBUG>=2) assert(!isnan(p->Sm) && p->Sm < 0.0);
	    }
	  }
	}
      }

      /* allocate YPen[c][refid][K=0..min(N-1,KMAX)][T=0..KMAX][I=1..N][n=1..DELTA_Y] */
      int kmax = min(N-1,KMAX);
      YPen[c][refid] = &cypen_ppp_pool[c][pm->YPmemPPP];
      long long YPcntPP = pm->YPmemPP, YPcntP = pm->YPmemP, YPcnt = pm->YPmem;
      long long YPcntMax = (refid < refend-1) ? pm[1].YPmem : YPmem;
      if(DEBUG) assert(YPcntMax  <= YPmem);

      for(int K = kmax; K >= 0; K--){
	if(DEBUG) assert(N > K);
	YPen[c][refid][K] = &cypen_pp_pool[c][YPcntPP]; YPcntPP += KMAX+1;
	for(int T = KMAX; T >= 0; T--){
	  YPen[c][refid][K][T] = &cypen_p_pool[c][YPcntP]; YPcntP += N+1;
	  CYPen *pmem = &cypen_pool[c][YPcnt]; YPcnt += N*DELTA_Y;
	  if(DEBUG>=2) assert(YPcnt <= YPcntMax);

	  //  #pragma omp parallel for num_threads(numthreads) schedule(static,16)
	  for(int I = N; I > K; I--){
	    CYPen *p = YPen[c][refid][K][T][I] = &pmem[-1];
	    pmem += DELTA_Y;
	    if(K > Kmax[I])
	      continue;
	    FLOAT Yik = Yc(Y,I,K);
	    FLOAT sm = Sm(0,I,K,Y,c);
	    for(int n = min(I-K-T-1,DELTA_Y); n > 0; n--){
	      if(DEBUG>=2) assert(0 < I-K-n);
	      if (T > Kmax[I-K-n])
		continue;
	      FLOAT Ynt = Yc(Y,I-K-n,T);
	      FLOAT y = Yik - Ynt;
	      FLOAT var = F2var[c] + S2var[c]*y;
	      if(QUADRATIC_VARIANCE)
		var += R2var[c] * y * y;
	      if(RES_VARIANCE){
		FLOAT resR = Y[I] - Y[I-K];
		FLOAT resL = Y[I-K-n] - Y[I-K-n-T];
		FLOAT resvar = resL*resL + resR*resR;
		var += E2var[c] * resvar;
	      }
	      if(DEBUG>=2) assert(var > 0.0);
	      if(DEBUG>=2) assert(fabs(PRtab[c][refid][n-1][I-K].LogPr - log(Pr(Y[I-K]-Y[I-K-n],c))) < 1.0e-7);
	      if(DEBUG>=2) assert(&p[n] < pmem);
	      if(DEBUG>=2) assert(&p[n] < &cypen_pool[c][YPmem]);
	      if(DEBUG>=2) assert(isnan(p[n].deltaY));
	      if(DEBUG>=2) assert(isnan(p[n].Pen));
	      if(DEBUG>=2) assert(isnan(p[n].Sm));
	      if(DEBUG>=2) assert(isnan(p[n].Ivar));
	      p[n].deltaY = Yik - Ynt;
	      if(FIX_OUTLIER_MIS){
		p[n].Pen = RefPen[c] - 0.5*mylog(var) - (n-1)*FnPenalty[c] + FpPenalty[c] + PRtab[c][refid][n-1][I-K].LogPr;
		p[n].Sm = sm;
	      } else
		p[n].Pen = RefPen[c] - 0.5*mylog(var) - (n-1)*FnPenalty[c] + FpPenalty[c] + PRtab[c][refid][n-1][I-K].LogPr + sm;
	      p[n].Ivar = 1.0/var;
	    }
	  }
	  if(DEBUG>=2) assert(pmem <= &cypen_pool[c][YPcnt]);
	}
      }

      if(DEBUG){
	if(refid < refend-1){
	  assert(PRcnt <= pm[1].PRmem);
	  assert(YPcntPP <= pm[1].YPmemPP);
	  assert(YPcntP <= pm[1].YPmemP);
	  if(!(YPcnt <= pm[1].YPmem)){
	    printf("c=%d,refid=%d(refstart=%d,refend=%d):YPcnt=%lld,MemCnt[refid].YPmem=%lld,MemCnt[refid+1].YPmem=%lld,N=%d,kmax=%d,KMAX=%d\n",
		   c,refid,refstart,refend,YPcnt,MemCnt[refid].YPmem,MemCnt[refid+1].YPmem,N,kmax,KMAX);
	    fflush(stdout);
	    assert(YPcnt <= pm[1].YPmem);
	  }
	} else {
	  assert(PRcnt <= PRmem);
	  assert(YPcntPP <= YPmemPP);
	  assert(YPcntP <= YPmemP);
	  assert(YPcnt <= YPmem);
	}
      }
    }
    delete [] MemCnt;
  }

  Ylambda12 = (Ylambda[0]*Ylambda[0] + Ylambda[1]*Ylambda[1])/(Ylambda[0]+Ylambda[1]);/* average interval distance between sites of opposite colors in YY */
  Xtheta12 = (Xtheta[0]*Xtheta[0] + Xtheta[1]*Xtheta[1])/(Xtheta[0]+Xtheta[1]);/* average interval distance between sites of opposite colors in XX */
  RefPen12 = BETA*(log(Xtheta12) - 0.5*log(M_PI));
  FM2var = 2.0 * SF[2] * SF[2];
  SM2var = 2.0 * fabs(SD[2]) * SD[2];
  if(VERB){
    printf("    Ylambda12=%0.6f,Xtheta12=%0.6f,RefPen12=%0.8f\n",Ylambda12,Xtheta12,RefPen12);
    fflush(stdout);
  }

  long long maxMN = max(maxN[0],maxM[0]);
  for(int c = 1; c < colors; c++){
    maxMN = max(maxMN,maxN[c]);
    maxMN = max(maxMN,maxM[c]);
  }

  OutlierPenaltyBC = new TFLOAT[maxMN];
  if(!outlierBC)
    memset(OutlierPenaltyBC, 0, maxMN*sizeof(TFLOAT));
  else {
    OutlierPenaltyBC[0] = 0.0;

    int nthreads = min((int)(maxMN/512), numthreads);

    #pragma omp parallel for num_threads(nthreads) schedule(static, 512) if(nthreads > 1)
    for(int i = 1; i < maxMN; i++)
      OutlierPenaltyBC[i] = -log((double)i);
  }
  if(VERB>=2){
    printf("score_init():maxMN=%lld,OutlierPenaltyBC=%p,outlierBC=%d\n",maxMN,OutlierPenaltyBC,outlierBC);
    fflush(stdout);
  }

  return 0;
}

/* allocate XPen and Qmin[],Qmax[] */

static void XPen_alloc(CXPen *XPen, int * &Qmin, int * &Qmax)
{
  if(DEBUG/* HERE >=2 */) assert(colors==2);

  /* allocate XPen-><field>[c=0..colors-1][m=1..DELTA_X][J=m..maxM+1] */
  if(DEBUG>=2 && !(XPen->deltaX[0] == 0)){
    #pragma omp critical
    {
      assert(XPen->deltaX[0] == 0);
    }
  }

  /* compute total memory needed */
  XPen->mem_pool_size = 0;
  int padding = 2*64;/* alignment fudge : padding for _mm512_loadu_ps (used in place of _mm512_maskloadu_ps) */
  for(int c = 0; c < colors; c++){
    int M = maxM[c]+2;
    XPen->mem_pool_size += DELTA_X * 4 * sizeof(TFLOAT*) + DELTA_X * M * 4 * sizeof(TFLOAT) + padding;
  }
  XPen->mem_pool = (char *) malloc(XPen->mem_pool_size);
  size_t offset = 0;

  for(int c = 0; c < colors; c++){
    int M = maxM[c]+2;

    XPen->deltaX[c] = (TFLOAT **) &XPen->mem_pool[offset];
    offset += DELTA_X * 4 * sizeof(TFLOAT*);

    XPen->deltaX[c]--;

    XPen->Pen[c] = &XPen->deltaX[c][DELTA_X];
    XPen->OutlierBias[c] = &XPen->deltaX[c][DELTA_X*2];
    XPen->PenBias[c] = &XPen->deltaX[c][DELTA_X*3];

#pragma ivdep
    for(int m = 1; m <= DELTA_X; m++){
      XPen->deltaX[c][m] = (TFLOAT *)&XPen->mem_pool[offset];
      offset += M * 4 * sizeof(TFLOAT);
      XPen->Pen[c][m] = &XPen->deltaX[c][m][M];
      XPen->OutlierBias[c][m] = &XPen->deltaX[c][m][M*2];
      XPen->PenBias[c][m] = &XPen->deltaX[c][m][M*3];
    }
    if(VERB>=3 && DELTA_X >=1 && c==1){
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num ();
#endif
      #pragma omp critical
      {
	printf("tid=%d,XPen=%p:XPen->deltaX[1]=%p,XPen->deltaX[1][1]=%p\n",tid,XPen,XPen->deltaX[1],XPen->deltaX[1][1]);
	fflush(stdout);
	if(DEBUG) assert(XPen->deltaX[1] != NULL && XPen->deltaX[1][1] != NULL);
      }
    }
    offset += padding;
  }

  if(DEBUG/* HERE >=2 */) assert(offset == XPen->mem_pool_size);

  int maxMM = 1+max(maxM[0],maxM[1]);
  Qmin = new int [maxMM];
  Qmax = new int [maxMM];
}

static void XPen_free(CXPen *XPen, int * &Qmin, int * &Qmax)
{
  free(XPen->mem_pool);
  XPen->mem_pool = NULL;
  XPen->mem_pool_size = 0;
  if(DEBUG>=2)
    for(int c = 0; c < colors; c++)
      XPen->deltaX[c] = XPen->Pen[c] = XPen->OutlierBias[c] = XPen->PenBias[c] = NULL;  

  delete [] Qmin; Qmin = NULL;
  delete [] Qmax; Qmax = NULL;
}

static void XPen_init(FLOAT **XX, int *MM, CXPen *XPen)
{
  int tid = 0;
#ifdef _OPENMP
  tid = omp_get_thread_num ();
#endif
  if(DEBUG>=3 && DELTA_X >= 1 && (!XPen->deltaX[1] || !XPen->deltaX[1][1])){
    #pragma omp critical
    {
      printf("tid=%d,XPen=%p:XPen->deltaX[1]=%p,XPen->deltaX[1][1]=%p:M=%d,%d,X=%p,%p\n",
	     tid,XPen,XPen->deltaX[1],XPen->deltaX[1][1],MM[0],MM[1],XX[0],XX[1]);
      fflush(stdout);
      assert(XPen->deltaX[1] != NULL && XPen->deltaX[1][1] != NULL);
    }
  }

  /* initialize XPen-><field>[c=0..1][m=1..DELTA_X][J=m..MM[c]+1] */
  for(int c = 0; c < colors; c++){
    if(DEBUG>=3 && DELTA_X >= 1) assert(XPen->deltaX[1] != NULL && XPen->deltaX[1][1] != NULL);

    FLOAT *X = XX[c];
    if(DEBUG) assert(XPen->deltaX[c] != NULL && MM[c] <= maxM[c]);

    for(int m = 1; m <= DELTA_X; m++){
      if(DEBUG>=3 && DELTA_X >= 1 && (!XPen->deltaX[c] || !XPen->deltaX[c][m])){
        #pragma omp critical
	{
	  printf("tid=%d,XPen=%p,c=%d,m=%d:XPen->deltaX[c]=%p,XPen->deltaX[c][m]=%p:M=%d,%d,X=%p,%p\n",
		 tid,XPen,c,m,XPen->deltaX[c],XPen->deltaX[c][m],MM[0],MM[1],XX[0],XX[1]);
	  fflush(stdout);
	  assert(XPen->deltaX[c] != NULL && XPen->deltaX[c][m] != NULL);
	}
      }

      TFLOAT *deltaXm = XPen->deltaX[c][m];
      TFLOAT *Penm = XPen->Pen[c][m];
      TFLOAT *OutlierBiasm = XPen->OutlierBias[c][m];
      TFLOAT *PenBiasm = XPen->PenBias[c][m];

      if(DEBUG>=3 && DELTA_X >= 1) assert(XPen->deltaX[1] != NULL && XPen->deltaX[1][1] != NULL);

#pragma ivdep
      for(int J = m; J <= MM[c]+1; J++){
	FLOAT deltaX = X[J] - X[J-m];
	TFLOAT Pen = m * FpPenalty[c] + deltaX * Frate[c];
	TFLOAT Bias = deltaX * FrateBias[c] + m * FBiasWT[c];
	deltaXm[J] = deltaX;
	Penm[J] = Pen;
	OutlierBiasm[J] = OutlierBias + Bias;
	PenBiasm[J] = Bias - Pen;

	if(DEBUG>=3 && DELTA_X >= 1) assert(XPen->deltaX[1] != NULL && XPen->deltaX[1][1] != NULL);
      }

      if(DEBUG>=3 && DELTA_X >= 1) assert(XPen->deltaX[1] != NULL && XPen->deltaX[1][1] != NULL);
    }
  }
}

#if 0 // OLD code

static void XPen_alloc(CXPen ***XPen, int * &Qmin, int * &Qmax)
{
  /* allocate XPen[c=0..colors-1][0..maxM[c]+1][1..DELTA_X] */
  for(int c = 0; c < colors; c++){
    int M = maxM[c];
    if(DEBUG) assert(XPen[c] == 0);
    XPen[c] = new CXPen*[M+2];
    for(int J = 1; J <= M+1; J++){
      XPen[c][J] = new CXPen[DELTA_X];
      XPen[c][J]--;
    }
  }

  int maxMM = 1+max(maxM[0],maxM[1]);
  Qmin = new int [maxMM];
  Qmax = new int [maxMM];
}

static void XPen_free(CXPen ***XPen, int * &Qmin, int * &Qmax)
{
  for(int c = 0; c < colors; c++)
    if(XPen[c]){
      int M = maxM[c];
      for(int J = 1; J <= M+1; J++)
	delete [] ++XPen[c][J];

      delete [] XPen[c];
      XPen[c] = 0;
    }

  delete [] Qmin; Qmin = NULL;
  delete [] Qmax; Qmax = NULL;
}

static void XPen_init(FLOAT **X, int *M, CXPen ***XPen)
{
  /* allocate XPen[c=0..1][J=0..M+1][m=1..DELTA_X] */
  for(int c = 0; c < colors; c++){
    if(DEBUG) assert(XPen[c] && M[c] <= maxM[c]);

    for(int J = 1; J <= M[c]+1; J++){
      int m = min(J,DELTA_X);
      for(CXPen *p = &XPen[c][J][m]; m >= 1; m--, p--){
	FLOAT delX = X[c][J] - X[c][J-m];
	TFLOAT Pen = m * FpPenalty[c] + delX * Frate[c];
	TFLOAT Bias = delX * FrateBias[c] + m * FBiasWT[c];
	p->deltaX = delX;
	p->Pen = Pen;
	p->OutlierBias = OutlierBias + Bias;
	p->PenBias = Bias - Pen;
      }
    }
  }
}
#endif

/* NOTE : this is a slower version of Sint() that does not use tables that depend of res or resSD and always uses double AND 
   forces the interval to be scored as an outlier (or NOT) based on "outlier" 
   Only used for debugging noise parameter estimation with REFDEBUG_STRICT >= 1
 */

static inline double Sint(FLOAT X,
                          FLOAT Y,
			  int m,
			  int n,
			  int J,
			  int I,
			  int K,/* at right end of interval */
			  int T,/* at left end of interval : only used with RES_VARIANCE */
			  FLOAT *Ya,/**< reference map array Ya[0..N+1] */
			  int outlier,/* force interval to be scored as outlier (or non-outlier) */
			  int c)
  
{
  if(DEBUG) assert(REFDEBUG_STRICT >= 1);
  if(outlierLambda < 1000.0 || outlierMax < 1000.0){
    printf("Sint() with REFDEBUG_STRICT not yet implemented with -outlierMax or -outlierLambda\n");
    exit(1);
  }

  if(DEBUG>=2) assert(c==0 || c==1);
  if(DEBUG>=2) assert(n>=1 && n <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2) assert(0 <= K && K <= KMAX);
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2) assert(BETA==1.0);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= KMAX);
    assert(I-K-n-T >= 1);
  }

  FLOAT err = X-Y;
  FLOAT var = F2var[c] + S2var[c] * Y;
  if(QUADRATIC_VARIANCE)
    var += R2var[c] * Y * Y;
  if(RES_VARIANCE){
    FLOAT resR = Ya[I] - Ya[I-K];
    FLOAT resL = Ya[I-K-n] - Ya[I-K-n-T];
    FLOAT resvar = resL*resL + resR*resR;
    var += E2var[c] * resvar;
  }
  if(DEBUG>=2) assert(var > 0.0);

  FLOAT Bias = X * FrateBias[c] + m * FBiasWT[c];
  FLOAT Pen = RefPen[c] - 0.5*mylog(var) - err*err/var - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] - X*Frate[c] + log(Pr(Ya[I-K]-Ya[I-K-n],c));
  FLOAT PenSm = Sm(0,I,K,Ya,c);
  if(outlier){
    FLOAT OutPen = OutlierPenalty[c];
    if(outlierBC)
      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
    if(FIX_OUTLIER_MIS)
      return Bias + OutlierBias + OutPen + PenSm;
    else
      return Bias + OutlierBias + OutPen;
  } else
    return Bias + Pen + PenSm;
}

/* NOTE : this is a slower version of Sint() that does not use tables that depend of res or resSD and always uses FLOAT */

static inline double Sint(FLOAT X,
			  FLOAT Y,
			  int m,
			  int n,
			  int J,
			  int I,
			  int K,/* at right end of interval */
			  int T,/* at left end of interval : only used with RES_VARIANCE */
			  FLOAT *Ya,/**< reference map array Ya[0..N+1] */
			  int c)
  
{
  if(DEBUG>=2) assert(c==0 || c==1);
  if(DEBUG>=2 && !(outlierExtend && !outlierExtendLimit)) assert(n>=1 && n <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2) assert(0 <= K && K <= KMAX);
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2) assert(BETA==1.0);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= KMAX);
    assert(I-K-n-T >= 1);
  }

  double err = X-Y;
  double var = F2var[c] + S2var[c] * Y;
  if(QUADRATIC_VARIANCE)
    var += R2var[c] * Y * Y;
  if(RES_VARIANCE){
    double resR = Ya[I] - Ya[I-K];
    double resL = Ya[I-K-n] - Ya[I-K-n-T];
    double resvar = resL*resL + resR*resR;
    var += E2var[c] * resvar;
  }
  if(DEBUG>=2) assert(var > 0.0);
  double Bias = X * FrateBias[c] + m * FBiasWT[c];
  double Pen = RefPen[c] - 0.5*mylog(var) - err*err/var - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] - X*Frate[c] + log(Pr(Ya[I-K]-Ya[I-K-n],c));
  double sm = Sm(0,I,K,Ya,c);
  if(maptype){
    double OutPen = OutlierPenalty[c];
    if(outlierBC)
      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
    if(OUTLIER_LTYPE==0)
      OutPen -= (X+Y) * OutlierLambdaInv;
    else
      OutPen -= fabs(X-Y) * OutlierLambdaInv;
    if(FIX_OUTLIER_MIS)
      return Bias + OutlierBias + max(Pen, OutPen) + sm;
    else
      return Bias + OutlierBias + max(Pen + sm, OutPen);
  } else if(OUTLIER_DELTA(X-Y)){
    double OutPen = OutlierPenalty[c];
    if(outlierBC)
      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
    OutPen -= fabs(X-Y) * OutlierLambdaInv;
    if(FIX_OUTLIER_MIS)
      return Bias + OutlierBias + max(Pen, OutPen) + sm;
    else
      return Bias + OutlierBias + max(Pen + sm, OutPen);
  } else
    return Bias + Pen + sm;
}


/* Sint() for resSD>0 : includes Sm(J,I,K,c) */

/* outlier penalty terms excluding OutlierPenalty or OutlierBias */
static inline TFLOAT OutPenExtend(TFLOAT X, TFLOAT Y, int m, int n)
{
  TFLOAT OutPen = outlierBC ? (OutlierPenaltyBC[m] + OutlierPenaltyBC[n]) : 0.0;
  if(maptype){
    if(OUTLIER_LTYPE==0)
      OutPen -= (X+Y) * OutlierLambdaInv;
    else
      OutPen -= fabs(X-Y) * OutlierLambdaInv;
  } else {
    OutPen -= fabs(X-Y) * OutlierLambdaInv;
  }
  return OutPen;
}

/* Can be computed faster in two steps:
  TFLOAT Y,Ivar,Sm;
  TFLOAT PenY = SintY(n,I,K,T,Y,Ivar,Sm,YPen[c][refid]);
  TFLOAT score = SintX(X,Y,PenY,Sm,Ivar,XPen[c][m][J].Pen,XPen[c][m][J].OutlierBias,XPen[c][m][J].PenBias,outlier,OutlierPenaltyBC[n] + OutlierPenaltyBC[m]);
*/
static inline TFLOAT Sint(TFLOAT X,
			  TFLOAT Y,
			  int m,
			  int n,
			  int J,
			  int I,
			  int K,/* at right end of interval! */
			  int T, /* at left end of interval : only used with RES_VARIANCE */
			  Cprtab **PRtab, /* actually PRtab[c][refid] */
			  FLOAT *Ya,/* reference map array Ya[c][0..N[c]+1] */
			  int c)
  
{
  if(DEBUG>=2) assert(c==0 || c==1);
  if(DEBUG>=2) assert(n>=1 && n <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2) assert(0 <= K && K <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2) assert(fabs(PRtab[n-1][I-K].LogPr - log(Pr(Ya[I-K]-Ya[I-K-n],c))) <= 1.0e-10);
  if(DEBUG>=2){
    TFLOAT sm = Sm(0,I,K,Ya,c);
    if(!(fabs(PRtab[K][I].Sm - sm) <= 1.0e-10)){
      #pragma omp critical 
      {
	printf("Sint(X=%0.3f,Y=%0.3f,m=%d,n=%d,J=%d,I=%d,K=%d,c=%d):Y[I]=%0.3f,Y[I-K]=%0.3f:PRtab[K][I].Sm=%0.11f,Sm(0,I,K,Ya)=%0.11f(err=%0.11f)\n",
	       X,Y,m,n,J,I,K,c,Ya[I],Ya[I-K],PRtab[K][I].Sm,sm,PRtab[K][I].Sm - sm);
	fflush(stdout);
	assert(fabs(PRtab[K][I].Sm - sm) <= 1.0e-10);
      }
    }
  }
  //  if(DEBUG>=2) assert(fabs(PRtab[K][I].Sm - Sm(0,I,K,Ya,c)) <= 1.0e-10);
  if(DEBUG>=2) assert(BETA==1.0);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= KMAX);
    assert(I-K-n-T >= 1);
  }

  TFLOAT err = X-Y;
  TFLOAT var = F2var[c] + S2var[c] * Y;
  if(QUADRATIC_VARIANCE)
    var += R2var[c] * Y * Y;
  if(RES_VARIANCE){
    FLOAT resR = Ya[I] - Ya[I-K];
    FLOAT resL = Ya[I-K-n] - Ya[I-K-n-T];
    TFLOAT resvar = resL*resL + resR*resR;
    var += E2var[c] * resvar;
  }
  if(DEBUG>=2) assert(var > 0.0);
  TFLOAT Bias = X * FrateBias[c] + m * FBiasWT[c];
  TFLOAT Pen = RefPen[c] - 0.5*mylog(var) - err*err/var - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] - X*Frate[c] + PRtab[n-1][I-K].LogPr;
  TFLOAT sm = PRtab[K][I].Sm;
  if(maptype){
    TFLOAT OutPen = OutlierPenalty[c];
    if(outlierBC)
      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
    if(OUTLIER_LTYPE==0)
      OutPen -= (X+Y) * OutlierLambdaInv;
    else
      OutPen -= fabs(X-Y) * OutlierLambdaInv;
    if(FIX_OUTLIER_MIS)
      return Bias + OutlierBias + max(Pen, OutPen) + sm;
    else
      return Bias + OutlierBias + max(Pen + sm, OutPen);
  } else if(OUTLIER_DELTA(X-Y)){
    TFLOAT OutPen = OutlierPenalty[c];
    if(outlierBC)
      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
    OutPen -= fabs(X-Y) * OutlierLambdaInv;
    if(FIX_OUTLIER_MIS)
      return Bias + OutlierBias + max(Pen, OutPen) + sm;
    else
      return Bias + OutlierBias + max(Pen + sm, OutPen);
  } else
    return Bias + Pen + sm;
}

static inline TFLOAT SintY(int n,
			   int I,
			   int K,/* at right end of interval! */
			   int T,/* at left end of interval */
			   TFLOAT &deltaY,
			   TFLOAT &Ivar,
			   TFLOAT &Sm,
			   CYPen ****YPenR) /* actually YPen[c][refid] */
{
  if(DEBUG>=2) assert(n>=1 && n <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2) assert(0 <= K && K <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2){
    assert(0 <= T && T <= KMAX);
    assert(I-K-n-T >= 1);
  }

  CYPen *p = &YPenR[K][T][I][n];
  deltaY = p->deltaY;
  Ivar = p->Ivar;
  if(FIX_OUTLIER_MIS) Sm = p->Sm;
  return p->Pen;
}

static inline TFLOAT SintX(TFLOAT X,
			   TFLOAT Y,
			   TFLOAT penY,
			   TFLOAT Sm,
			   TFLOAT Ivar,
			   TFLOAT XPen,
			   TFLOAT XOutlierBias,
			   TFLOAT XPenBias,
			   TFLOAT OutPenBC,
			   int c)

{
  if(DEBUG>=2) assert(c==0 || c==1);

  TFLOAT err = X-Y;
  TFLOAT XYpen = penY - err*err*Ivar;

  // HERE HERE : continue update from RGentigRefScore.h to support -maxOutlier
  if(maptype){
    TFLOAT Pen = XYpen - XPen;
    TFLOAT OutPen = OutlierPenalty[c] + OutPenBC;
    if(OUTLIER_LTYPE==0)
      OutPen -= (X+Y) * OutlierLambdaInv;
    else
      OutPen -= fabs(err) * OutlierLambdaInv;
    if(FIX_OUTLIER_MIS)
      return XOutlierBias + max(Pen, OutPen) + Sm;
    else /* in this case Pen includes Sm */
      return XOutlierBias + max(Pen, OutPen);
  } else if(OUTLIER_DELTA(X-Y)){
    TFLOAT Pen = XYpen - XPen;
    TFLOAT OutPen = OutlierPenalty[c] + OutPenBC;
    OutPen -= fabs(X-Y) * OutlierLambdaInv;
    if(FIX_OUTLIER_MIS)
      return XOutlierBias + max(Pen, OutPen) + Sm;
    else /* in this case Pen includes Sm */
      return XOutlierBias + max(Pen, OutPen + OutPenBC);
  }
  if(FIX_OUTLIER_MIS)
     return XPenBias + XYpen + Sm;
  else
    return XPenBias + XYpen;/* in this case Pen includes Sm */
}

/* SintDetail with resSD>0.0 */
static void SintDetail(double X,
		       double Y,
		       int m,
		       int n,
		       int J,
		       int I,
		       int K,
		       int T,/* at left end of interval : only used with RES_VARIANCE */
		       Cprtab **PRtab, /* actually PRtab[c][refid] */
		       FLOAT *Ya,/* reference array Ya[0..N+1] */
		       int c,
		       double &Bias,
		       double &Pen,
		       double &PenSm,
		       int verb)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  if(DEBUG>=2) assert(0 <= K && K <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2){
    TFLOAT sm = Sm(J,I,K,Ya,c);
    if(!(fabs(PRtab[K][I].Sm - sm) <= 1.0e-10)){
      #pragma omp critical
      {
	printf("SintDetail(x=%0.4f,y=%0.4f,m=%d,n=%d,J=%d,I=%d,K=%d,T=%d,c=%d):Ya[I]=%0.6f,Ya[I-K]=%0.6f:PRtab[K][I].Sm=%0.11f,Sm(J,I,K,Ya,c)=%0.11f(err=%0.11f)\n",
	       X,Y,m,n,J,I,K,T,c,Ya[I],Ya[I-K],PRtab[K][I].Sm,sm,PRtab[K][I].Sm -sm);
	fflush(stdout);
	assert(fabs(PRtab[K][I].Sm - sm) <= 1.0e-10);
      }
    }
  }
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= max(max(2,KMAX),DELTA_Y));
    assert(I-K-n-T >= 1);
  }

  PenSm = PRtab[K][I].Sm;
  if(n > DELTA_Y){
    if(DEBUG>=2) assert(outlierExtend > 0);
    if(DEBUG>=2) assert(n >= 1);
    Bias = 0.0;
    Pen = -999999.9999;
  } else {

    if(DEBUG>=2) assert(n>=1 && n <= max(max(2,KMAX),DELTA_Y));
    if(DEBUG>=2) assert(fabs(PRtab[n-1][I-K].LogPr - log(Pr(Ya[I-K]-Ya[I-K-n],c))) <= 1.0e-10);

    double err = X-Y;
    double var = F2var[c] + S2var[c] * Y;
    if(QUADRATIC_VARIANCE)
      var += R2var[c] * Y * Y;
    if(RES_VARIANCE){
      FLOAT resR = Ya[I] - Ya[I-K];
      FLOAT resL = Ya[I-K-n] - Ya[I-K-n-T];
      double resvar = resL*resL + resR*resR;
      var += E2var[c] * resvar;
    }
    if(DEBUG>=2) assert(var > 0.0);

    Bias = X * FrateBias[c] + m * FBiasWT[c];
    Pen = BETA * (RefPen[c] - 0.5*log(var)) - err*err/var - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] - X*Frate[c] + PRtab[n-1][I-K].LogPr;
  
    if(verb){
      printf("    c=%d,I=%d,K=%d,T=%d,J=%d,n=%d,m=%d:PenVar=%0.6f, PenErr=%0.6f, X*Frate=%0.6f, X*FrateBias=%0.6f, m*FBiasWT=%0.6f, PenPr=%0.6f, PenSm=%0.6f,%0.6f (err=%0.6f,var=%0.8f,sf=%0.6f,sd=%0.6f,sr=%0.6f,se=%0.6f)\n", 
	     c,I,K,T,J,n,m, BETA * (RefPen[c] - 0.5*log(var)), -err*err/var, X*Frate[c], X*FrateBias[c], m*FBiasWT[c], PRtab[n-1][I-K].LogPr, PRtab[K][I].Sm, Sm(0,I,K,Ya,c), err, var, SF[c],SD[c],SR[c],SE[c]);
      fflush(stdout);
    }
  }
}

#if USE_RFLOAT==1
static void SintDetail(float X,
		       float Y,
		       int m,
		       int n,
		       int J,
		       int I,
		       int K,
		       int T,/* at left end of interval : only used with RES_VARIANCE */
		       Cprtab **PRtab, /* actually PRtab[c][refid] */
		       FLOAT *Ya,/* reference array Ya[0..N+1] */
		       int c,
		       float &Bias,
		       float &Pen,
		       float &PenSm,
		       int verb)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  if(DEBUG>=2) assert(0 <= K && K <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2){
    TFLOAT sm = Sm(J,I,K,Ya,c);
    if(!(fabs(PRtab[K][I].Sm - sm) <= 1.0e-10)){
      #pragma omp critical
      {
	printf("SintDetail(x=%0.4f,y=%0.4f,m=%d,n=%d,J=%d,I=%d,K=%d,T=%d,c=%d):Ya[I]=%0.6f,Ya[I-K]=%0.6f:PRtab[K][I].Sm=%0.11f,Sm(J,I,K,Ya,c)=%0.11f(err=%0.11f)\n",
	       X,Y,m,n,J,I,K,T,c,Ya[I],Ya[I-K],PRtab[K][I].Sm,sm,PRtab[K][I].Sm -sm);
	fflush(stdout);
	assert(fabs(PRtab[K][I].Sm - sm) <= 1.0e-10);
      }
    }
  }
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= max(max(2,KMAX),DELTA_Y));
    assert(I-K-n-T >= 1);
  }

  PenSm = PRtab[K][I].Sm;
  if(n > DELTA_Y){
    if(DEBUG>=2) assert(outlierExtend > 0);
    if(DEBUG>=2) assert(n >= 1);
    Bias = 0.0;
    Pen = -999999.9999;
  } else {

    if(DEBUG>=2) assert(n>=1 && n <= max(max(2,KMAX),DELTA_Y));
    if(DEBUG>=2) assert(fabs(PRtab[n-1][I-K].LogPr - log(Pr(Ya[I-K]-Ya[I-K-n],c))) <= 1.0e-10);

    float err = X-Y;
    float var = F2var[c] + S2var[c] * Y;
    if(QUADRATIC_VARIANCE)
      var += R2var[c] * Y * Y;
    if(RES_VARIANCE){
      FLOAT resR = Ya[I] - Ya[I-K];
      FLOAT resL = Ya[I-K-n] - Ya[I-K-n-T];
      float resvar = resL*resL + resR*resR;
      var += E2var[c] * resvar;
    }
    if(DEBUG>=2) assert(var > 0.0);

    Bias = X * FrateBias[c] + m * FBiasWT[c];
    Pen = BETA * (RefPen[c] - 0.5*log(var)) - err*err/var - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] - X*Frate[c] + PRtab[n-1][I-K].LogPr;
  
    if(verb){
      printf("c=%d,I=%d,K=%d,T=%d,J=%d,n=%d,m=%d:PenVar = %0.6f, PenErr=%0.6f, X*Frate=%0.6f, X*LBias=%0.6f, m*FBias=%0.6f, PenPr=%0.6f, PenSm=%0.6f,%0.6f (err=%0.6f,sf=%0.6f,sd=%0.6f,sr=%0.6f,se=%0.6f)\n", 
	     c,I,K,T,J,n,m, BETA * (RefPen[c] - 0.5*log(var)), -err*err/var, X*Frate[c], X*LBias[c], m*FBias[c], PRtab[n-1][I-K].LogPr, PRtab[K][I].Sm, Sm(0,I,K,Ya,c), err, SF[c],SD[c],SR[c],SE[c]);
      fflush(stdout);
    }
  }
}
#endif

/* short for Send(X,m,n) + Sint(X,Y,1,1) - Sint((X+Y)/2,(X+Y)/2,1,1)  with resSD > 0.0*/
static inline TFLOAT Sbnd(TFLOAT X,
			  TFLOAT Y,
			  int m,
			  int n,
			  int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  if(DEBUG>=2)assert(n>=1);
  if(DEBUG>=2)assert(m>=1);
  TFLOAT err = X-Y;
  TFLOAT var = F2var[c] + S2var[c] * Y;
  if(QUADRATIC_VARIANCE)
    var += R2var[c] * Y*Y; 

  TFLOAT Bias = X*FrateBias[c] + (m-1)*FBiasWT[c] + EBiasWT[c];
  return Bias - (n-1)*FnPenalty[c] - (m-1)*FpPenalty[c] - X*Frate[c] - err*err/var;
}

static inline TFLOAT Send(TFLOAT X,
			  int m,
			  int n,
			  int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  TFLOAT Bias = X * FrateBias[c] + (m-1)*FBiasWT[c] + EBiasWT[c];
  return Bias - (n-1)*FnPenalty[c] -(m-1)*FpPenalty[c] - X*Frate[c];
}

static inline TFLOAT BiasEnd(TFLOAT X,
			     int m,
			     int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  return X * FrateBias[c] + (m-1)*FBiasWT[c] + EBiasWT[c];
}

static inline TFLOAT BiasEnd2(TFLOAT X,
			      int m,
			      int c)
{
  TFLOAT r = BiasEnd(X,m,c);
  return max((TFLOAT)0.0, r);
}

/* 2-color misalignment penalty (based on the rightmost alignmed sites of each color) */
static inline TFLOAT SMA(FLOAT X, /* difference between rightmost aligned sites of each color on Nanomap */
			 FLOAT Y, /* difference between rightmost aligned sites of each color on Reference */
			 FLOAT len, /* average of rightmost aligned sites on Query or Reference from left end) */
			 int Xorientation)
{
  if(DEBUG>=2 && COLOR_SHIFT) assert(len > 0.0);

  FLOAT err = X - Y + (Xorientation ? MA_mean : -MA_mean);
  FLOAT var = FM2var + SM2var * (COLOR_SHIFT > 0 ? len : fabs(Y));

  if(DEBUG>=2) assert(var > 0.0);
  if(DEBUG>=2) assert(BETA==1.0);

  return RefPen12 - err*err/var - 0.5*log(var);
}

#endif // RGENTIGREFSCORE2_H
