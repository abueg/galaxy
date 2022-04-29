#include <stdlib.h>
#include <stdio.h>

// #include <immintrin.h>

#include "globals.h"
#include "parameters.h"

#ifndef RGENTIGREFSCORE_H
#define RGENTIGREFSCORE_H

static Ident RGentigRefScore_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/RGentigRefScore.h 11384 2020-07-29 20:08:32Z tanantharaman $");

extern double mtime();/* microsecond timer function : see refine.cpp */
extern double wtime();

#define LAMBDA_FIX0 outlierNorm0 /* Don't adjust outlier penalty term normalization factor based on outlierLambda value (if maptype==0)
						    If >= 2 : also don't adjust based on outlierMax */
#define LAMBDA_FIX1 outlierNorm1 /* Don't adjust outlier penalty term normalization factor based on outlierLambda value (if maptype==1) */ 

#define OUTLIER_TYPE OutlierType0 /* WAS 1 : slightly faster if compile time constant */
                       /* OUTLIER type (when maptype==0) :
			  0 : Include both misaligned sites and sizing error in outlier (requires DELTA_X = DELTA_Y >= n to evaluate exactly)
			  1 : Include only sizing error in outlier
			  NOTE : Only applied if maptype==0. If maptype==1 uses OUTLIER_TYPE1 instead 
		       */

#define OUTLIER_TYPE1 OutlierType1 /* OUTLIER type (when maptype==1) :
			  0 : Include both misaligned sites and sizing error in outlier (requires DELTA_X = DELTA_Y >= n to evaluate exactly)
			  1 : Include only sizing error in outlier
			  NOTE : Only applied if maptype==1. If maptype==0 uses OUTLIER_TYPE instead 
		       */


#define OUTLIER_LTYPE 1 /* OUTLIER exponential decay type :
			   0 : outlier likelihood is proportional to exp(-(x+y)/outlierLambda)
			   1 : outlier likelihood is proportional to exp(-|x-y|/outlierLambda) 
			   NOTE : Only implemented with maptype=1 : If maptype==0, always behaves like OUTLIER_LTYPE==1.
			*/

#define OUTLIER_DELTA(err) (!OUTLIER_MAX || fabs(err) < ((RFLOAT)outlierMax)) 

#define OUTLIER(maptype,x,y) (maptype || OUTLIER_DELTA(x-y)) /* only test alignment intervals for outliers if they satisfy this condition : x = Nanomap, y = reference (interval size), maptype = 0(BNX) or 1(CMAP) */

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

static RFLOAT *OutlierPenaltyBC = NULL;

FLOAT **PRsum = NULL;/* PRsum[refid][I] = sum(i=1..I-1) Pr(Y[I-i+1]-Y[I-i]) */

static class Cprtab {
 public:
  FLOAT FnPenalty; /** max(0.0, PRsum[refid][I] - PRsum[refid][I-n] - 1.0) * FnPenalty, where n == K+1 */
  FLOAT LogPr;/**< log(Pr(Y[I]-Y[I-n])), where n == K+1 */
  RFLOAT Sm;/**< Sm(0,I,K,Y) */
} ***PRtab = NULL;/**< PRtab[refid][K][I] has precomputed penalties that depend on refmap[refid]->site[0][I-K ... I] */

static class CYPen {
 public:
  RFLOAT deltaY;/**< Yc(YY[refid],I,K)-Yc(YY[refid],I-K-n,T) */
  RFLOAT Pen;/**< -(FN_ADJUST ? PRtab[refid][n-1][I-K].FnPenalty : (n-1)*FnPenalty) + FpPenalty + PRtab[refid][n-1][I-K].LogPr */
  RFLOAT Gauss;/**< RefPen - 0.5*log(var) */
  RFLOAT Sm;/**< Sm(0,I,K,Y) */
  RFLOAT Ivar;/**< 1.0/var where var = F2var + S2var*deltaY + R2var*deltaY*deltaY + E2var*(resL*resL+resR*resR) */
} *****YPen = NULL;/**< YPen[refid][K=0..min(N-1,KMAX)][T=0..KMAX][I=K+1..N][n=1..DELTA_Y] */

class CMemCnt {
 public:
  long long YPmemPPP, YPmemPP, YPmemP, YPmem;
  long long PRmemP,PRmem;
};

static Cprtab **cprtab_p_pool = NULL;
static Cprtab *cprtab_pool = NULL;
static CYPen **** cypen_ppp_pool = NULL;
static CYPen *** cypen_pp_pool = NULL;
static CYPen **  cypen_p_pool = NULL;
static CYPen *   cypen_pool   = NULL;
static size_t score_init_mem = 0;/* total memory that will be freed by score_free() */

class CXPen {
public:
  RFLOAT **deltaX;/* X[J] - X[J-m] */
  RFLOAT **Pen;/* m*FpPenalty */ // WAS + deltaX * Frate
  RFLOAT **Bias; /* deltaX*FrateBiasWT + m*FBiasWT */
  RFLOAT **PenBias; /**< deltaX*FrateBiasWT + m*FBiasWT - Pen */
  char *mem_pool;
  size_t mem_pool_size;
  CXPen() { 
    mem_pool = NULL;
    mem_pool_size = 0;
    if(DEBUG>=2) deltaX = Pen = Bias = PenBias = NULL;  
  };
} ; /* XPen[k=0..nummaps-1].<field>[m=1..DELTA_X][J=m+1..M+1] */

static int grefend = 0;

static long long maxN = 0;
static long long maxM = 0;

/* Gentig Guided Alignment scoring functions */

/* precompute values to speed up likelihood score computation */
static RFLOAT LBias, FBias, EBias, PRbias = 0.0, FnPenalty, FpPenalty, MinScore;
static RFLOAT OutlierPenalty, OutlierPenalty2, OutlierLambdaInv, pTP, log_pTP;// FpPenaltyBias;
RFLOAT OutlierBias, OutlierEndBias, OutlierEndPenalty, Frate;
static RFLOAT biasWToutlierF; // biasWToutlier

static RFLOAT FrateBiasWT;/* biasWT * (Frate + LBias) */
static RFLOAT FBiasWT;/* biasWT * FBias */
static RFLOAT EBiasWT;/* biasWT * EBias */

static RFLOAT F2var;/* 2*SF[0]*SF[0] */
static RFLOAT S2var;/* 2*fabs(SD[0])*SD[0] */
static RFLOAT R2var; /* 2*SR[0]*SR[0] */
static RFLOAT E2var; /* 2*SE[0]*SE[0] */
static RFLOAT RefPen;/* log(Xtheta) - 0.5*log(M_PI) */

static double IresSD;/* 1.0/(resSD * PixelLen * sqrt(2.0)) */
static double pTPd, log_pTPd, log_qTPd;
static double Ylambda, Xlambda, Xtheta = 0.0, resKB, resKB2, Ygm, GsitesPerMap, theta;
static int Rnummaps;

double raYlambda = -1.0;/* global copy of Ylambda used by refalign */

#if (USE_SSE==1 && USE_AVX==1 && USE_RFLOAT==1)
__m256 v256_Frate, v256_OutlierPenalty, v256_micro_unit, v256_outlierMax, v256_outlierNegMax, v256_OutlierLambdaInv, v256_biasWToutlierF;
__m256 v256_OUTLIER_TYPE, v256_OUTLIER_TYPE1;
__m256 mask8to256[256];/* table to convert 8 bit mask to 8x32 bit mask : bit 0 maps to bits 0-31, bit 1 maps to bits 32..63 etc */
#ifdef __AVX2__
__m256i v256i_76543210;/* constant 32-bit integer values 0,1,2,3,4,5,6,7 in a 256 bit vector (0 in lowest order 32 bits and 7 in highest order 32 bits)*/
#else
__m128i v128i_3210, v128i_7654;/* constants 32-bit integer values (0,1,2,3) and (4,5,6,7) */
#endif
#endif

#if (USE_MIC==1 && USE_RFLOAT==1)
__m512 v512_Frate, v512_OutlierPenalty, v512_micro_unit, v512_outlierMax, v512_outlierNegMax, v512_OutlierLambdaInv, v512_biasWToutlierF, v512_true;
__m512 v512_OUTLIER_TYPE, v512_OUTLIER_TYPE1;
#endif


static inline double Pr(double y)
{
  double pr = 0.5*(1.0+erf((y - resKB)*IresSD));
  return max(1.0e-100, pr);
}

/** Pn(y) == 1-Pr(y) */
static inline double Pn(double y)
{
  return 0.5*erfc((y - resKB)*IresSD);
}

static inline double log_Pn(double y)
{
  return log(Pn(y));
}

/* NOTE : return value can be -Inf if Y[I] - Y[I-K] is large */
double Sm(int J, int I, int K, double *Y)
{
  if(SCORE_APPROX >= 2){
    if(!K)
      return log_pTPd;
    if(K==1)
      return 2.0 * log_pTPd + log_Pn(Y[I]-Y[I-K]);

    return 2.0 * log_pTPd + log_Pn(Y[I]-Y[I-K]) + (K-1.0) * log_qTPd;

  } else if(MIS_VITERBI){/* score = (K+1) * log_pTPd + sum(t = 0 .. K-1) Log_Pn(Y[I-t]-Y[I-t-1]), along with pTPcnt += K+1 */
    double ret = (K+1.0) * log_pTPd;
    for(int t = 0; t < K; t++)
      ret += log_Pn(Y[I-t]-Y[I-t-1]);
    return ret;
    
  } else {/* The following evaluation is NOT consistent with Viterbi scoring used elsewhere and will give higher score to alignments with misresolved sites, since
	     the Likelihood for all possible ways of resolving misresolved sites are being added up. To avoid this bias use the Viterbi score (see above MIS_VITERBI) 
	     NOTE : the downside of Viterbi scoring (above) is that it results in poor scores in high label density regions with many misresolved sites, resulting is spurious endoutliers, especially with biaswt 0.7 -biaswtEnd 0
	  */
    double p = pTPd;

#if RESSD==0 // NEW CODE : only K(K+1)/2 calls to Pn() instead of (K+1)/2 * 2^(K-1) calls AND always exact
    if(K <= 0)
      return log_pTPd;

    double *PMK = (double *)alloca(KMAX * sizeof(double));
    double *PowQ = (double *)alloca((KMAX+1) * sizeof(double));
    if(DEBUG>=2) assert(1 <= K && K <= KMAX);
    double q = 1.0 - p;

    PowQ[0] = 1.0;
    for(int k = 0; k < K; k++){
      PowQ[k+1] = q * PowQ[k];      /* PowQ[t = 0..k+1] == pow(q,t) */
      double Yk = Y[I-k-1];
      double sum = Pn(Y[I]-Yk) * PowQ[k];
      for(int t = 1; t <= k; t++)
	sum += Pn(Y[I-t] - Yk) * PowQ[k-t] * p * PMK[t-1];
      PMK[k] = sum;
    }
    return 2.0 * log_pTPd + log(PMK[K-1]);

#else // OLD CODE
    switch(K){
    case 0:
      return log_pTP;
    case 1:
      return 2*log_pTP+log_Pn(Y[I]-Y[I-1]);
    case 2:
      if(RESSD)
	return 3*log_pTP+log_Pn(Y[I]-Y[I-1])+log_Pn(Y[I-1]-Y[I-2]);
      else
	return 2*log_pTP+log((p*Pn(Y[I]-Y[I-1])*Pn(Y[I-1]-Y[I-2]) + 
			  (1.0-p)*Pn(Y[I]-Y[I-2])));
    case 3:
      if(RESSD)
	return 4*log_pTP+log_Pn(Y[I]-Y[I-1])+log_Pn(Y[I-1]-Y[I-2])+log_Pn(Y[I-2]-Y[I-3]); 
      else
	return 2*log_pTP+log((p*p * Pn(Y[I]-Y[I-1])*Pn(Y[I-1]-Y[I-2])*Pn(Y[I-2]-Y[I-3]) + 
			  p*(1.0-p) * (Pn(Y[I]-Y[I-2])*Pn(Y[I-2]-Y[I-3])+Pn(Y[I]-Y[I-1])*Pn(Y[I-1]-Y[I-3])) +
			  (1.0-p)*(1.0-p) * Pn(Y[I]-Y[I-3])));    
    default:  /* general case K>3 is approximate */
      double ret = 0.0;
      if(RESSD)
	ret = p;
      else {
	for(int j=1; j < K; j++){
	  double den = Pn(Y[I-j+1]-Y[I-j])*Pn(Y[I-j]-Y[I-j-1]);
	  //	if(den > 0.0)
	  ret += Pn(Y[I-j+1]-Y[I-j-1])/den;
	}
	ret = p + (1.0 - p)*ret;
      }
      for(int i = 1; i <= K; i++)
	ret *= p*Pn(Y[I-i+1]-Y[I-i]);
      //    if(ret <= 0.0)
      //      return -1.0e+300;
      return log(ret);
    }
#endif // OLD CODE
  }
}

double SmRA(int J, int I, int K, double *Y)
{
  return Sm(J,I,K,Y);
}

static void score_free(int refstart, int refend)
{
  if(VERB>=2){
    printf("score_free():refstart=%d,refend=%d\n",refstart,refend);
    fflush(stdout);
  }
  if(PRsum){
    if(PRsum[0])
      delete PRsum[0];
    delete [] PRsum;
    PRsum = NULL;
  }

  if(PRtab){
    delete [] cprtab_p_pool; cprtab_p_pool = NULL;
    delete [] cprtab_pool; cprtab_pool = NULL;
    delete [] PRtab;
    PRtab = NULL;
  }
  if(YPen){
    delete [] cypen_ppp_pool; cypen_ppp_pool = NULL;
    delete [] cypen_pp_pool; cypen_pp_pool = NULL;
    delete [] cypen_p_pool; cypen_p_pool = NULL;
    delete [] cypen_pool; cypen_pool = NULL;
    delete [] YPen;
    YPen = NULL;
  }
  delete [] OutlierPenaltyBC; OutlierPenaltyBC = NULL;
  score_init_mem = 0;
}

//static int Verb=0;

/** return value is 1, if Y[] has no intervals above minKB */
/* NOTE : KKmax[refid = refstart .. refend][I=1..N] contains precomputed maximum K values */

static int score_init(Cmap **refmap,int refstart, int refend, int nummaps, int **RKmax)
{
  grefend = refend;

  score_free(refstart,refend);

  if(OUTLIER_TYPE != OutlierType0){
    printf("WARNING: -outlierType0 %0.3f ignored : using fixed value %0.3f instead (see OUTLIER_TYPE in RGentigRefScore.h)\n", OutlierType0, OUTLIER_TYPE);
    fflush(stdout);
  }
  if(OUTLIER_LTYPE==0 && maptype==0){
    printf("RGentigRefScore.h: OUTLIER_LTYPE==0 not implemented with maptype=%d: Assuming OUTLIER_LTYPE==1\n",maptype);
    fflush(stdout);
  }

  if(OUTLIER_TYPE1 && maptype==1){
    if(0 && DEFER_BP && !EXTENDONLY && DELTA_X==4 && DELTA_Y==6 && !outlierExtend && !outlierBC){// remove after implementing with vector code
      printf("RGentigRefScore.h: -outlierType1 1 not implemented with maptype=%d except with -outlierBC OR -outlierExtend OR -deltaX > 4 OR -deltaY > 6\n",maptype);
      fflush(stdout);exit(1);
    }
  }
  if(FAST_OUTLIER != 1){
    printf("RGentigRefScore.h: FAST_OUTLIER != 1 no longer supported\n");
    fflush(stdout);exit(1);
  }
  if(Beta != 1.0){
    printf("Beta != 1.0 no longer supported\n");
    fflush(stdout);exit(1);
  }
  if(BETA != 1.0){
    printf("BETA != 1.0 no longer supported\n");
    fflush(stdout);exit(1);
  }
  if(Refine && maptype != 0){
    printf("WARNING:refinement not supported with maptype=%d: using -maptype 0\n",maptype);
    fflush(stdout);

    maptype = 0;
  }

  if(DEBUG) assert(colors==1);

  resKB = res[0] * PixelLen;
  resKB2 = 0.0;

  double pFN = FN[0];
  pTP = pTPd = 1.0 - pFN;
  log_pTP = log_pTPd = log(pTPd);
  log_qTPd = log(pFN);
  double F = FP[0]*0.01;
  double var = fabs(SD[0])*SD[0];
  double varF = SF[0]*SF[0];
  double varR = SR[0]*SR[0];
  double varE = SE[0]*SE[0];

  Ygm = 0.0;/* to compute Geometric Mean of Y intervals */
  int Ycnt = 0;
  double Ysum = 0.0;

  if(VERB>=2){
    printf("score_init: %d reference contigs (res=%0.4f,resSD=%0.4f)\n", refend-refstart,res[0],resSD[0]);
    fflush(stdout);
  }
  for(int refid = refstart; refid < refend; refid++){
    Cmap *Ymap = refmap[refid];
    int N = Ymap->numsite[0];
    if(DEBUG) assert(N >= 0 && N <= (int)MAXINT);
    FLOAT *Y = Ymap->site[0];/* Y[0..N+1] */
    for(int i=1;i<N;i++){
      double y = Y[i+1]-Y[i];
      if(y > resKB){
	Ygm += log(y);
	Ysum += y;
	Ycnt++;
      }
    }
  }
  if(!Ycnt || Ysum <= 0.0){
    printf("WARNING:No intervals over resKB=%0.3f in any of %d..%d reference maps: Ycnt=%d,Ysum=%0.3f (Defaulting to Ycnt=1,Ysum=6.0kb)\n",resKB,refstart,refend-1,Ycnt,Ysum);
    fflush(stdout);

    Ycnt = 1;
    Ysum = 6.0;
    Ygm = log(6.0);
  }

  Ygm /= max(1,Ycnt);
  Ygm = exp(Ygm);

  Ylambda = Ysum / max(1,Ycnt);

  /* compute theta based on X */
  theta = 0.0;
  long long sitecnt = 0;
  maxM = 0;
  for(int i= 0; i < nummaps;i++){
    Cmap *Xmap = Gmap[i];
    int M = Xmap->numsite[0];
    sitecnt += M;
    theta += Xmap->site[0][M+1];
    if(M > maxM)
      maxM = M;
  }
  if(DEBUG) assert(0 <= maxM && maxM <= MAXINT);

  Rnummaps = nummaps;

  if(XmapCount[0] > 0){/* use command line values for XmapCount,XmapSites,XmapLength (adjusted for any change in PixelLen) */
    if(VERB>=2){
      printf("Using XmapCount[0]=%d,XmapSites=%lld,XmapLength=%0.8f,XmapGtheta=%0.8f, bpp= %0.2f\n",XmapCount[0],XmapSites[0],XmapLength[0],XmapGtheta[0], 500.0 * PixelLen/ origPixelLen);
      fflush(stdout);
    }
    Rnummaps = XmapCount[0];
    sitecnt = XmapSites[0];
    theta = XmapLength[0] * PixelLen / origPixelLen;
  }

  if(sitecnt <= 0 || theta <= 0.0){
    fprintf(stderr,"score_init():Insufficient -i map data:nummaps=%d,total sites=%lld,total length=%0.3f\n",
	    Rnummaps,sitecnt,theta);
    fflush(stdout);exit(1);
  }
  Xtheta = theta/((double)sitecnt) - resKB2*0.5;
  if(VERB>=2 && XmapCount[0]){
    printf("score_init():XmapLength[0]=%0.6f,XmapSites[0]=%lld,PixelLen=%0.12f,origPixelLen=%0.12f,resKB2=%0.6f:theta=%0.10f,sitecnt=%lld,Xtheta=%0.10f\n",
	   XmapLength[0],XmapSites[0],PixelLen,origPixelLen,resKB2,theta,sitecnt,Xtheta);
    fflush(stdout);
  }
  
  double R = RTHETA_FIX ? Rnummaps/theta : 0.0;/* should be reduced by inverse of average linear chromosome size */
  GsitesPerMap = ((double)sitecnt)/Rnummaps;

  if(F*Xtheta >= 1.0){
    printf("WARNING:FP=%0.4f (per 100 kb), average site interval Xtheta=%0.3f (kb): FP density exceeds total site density\n",F*100.0,Xtheta);
    printf("  Reducing FP to %0.4f\n", 50.0/Xtheta);
    F = 0.5/Xtheta;
    FP[0] = F*100.0;
  }

  Xlambda = pTP/(1.0/(Xtheta - var*0.5) - F);/* HERE : does not include correction for SR[],SE[] */
  if(DEBUG && !SCORE_Y && !(Xlambda > 0.0)){
    printf("Xlambda=%0.6e : pTP=%0.6e,Xtheta=%0.6e (theta=%0.6e,resKB2=%0.6e,sitecnt=%lld),var=%0.6e,F=%0.6e\n",
	   Xlambda, pTP, Xtheta, theta/sitecnt, resKB2, sitecnt, var, F);
    fflush(stdout);
    assert(Xlambda > 0.0);
  }

  if(XmapCount[0] /* WAS 169 && RefRepeats==1 */){
    Ylambda = Xlambda;
    Ygm = XmapGtheta[0] * /* NEW12 */ (PixelLen / origPixelLen) * Xlambda/Xtheta;
  } else if(SCORE_Y){ /* make score dependent only on Y (not X) : for debugging only*/
    Xlambda = Ylambda;
    Xtheta = 1.0/(F + pTP/Xlambda) + var*0.5;
  }

  if(DEBUG){
    assert(FN[0] > 0.0);
    assert(FN[0] < 1.0);
    assert(theta > 0.0);
    assert(Ylambda > 0.0);
    assert(Xlambda > 0.0);
    assert(Xtheta > 0.0);
  }

  Frate = F + R - thetaScale/Xtheta; // WAS F - 1.0/Xtheta

  double YgmBpTP = (Ygm/pTPd) * (Xlambda/Ylambda);
  double Sbias = log(Xtheta) - 0.5*log(2.0*M_PI*(varF+var*YgmBpTP+varR*YgmBpTP*YgmBpTP));
  FBias = 0.5*KBIAS - KFBIAS*pFN*log(pFN)/pTPd - log_pTPd - Sbias; /* Bias per site interval in X (overlapping Y)  = Abias - Sbias in document (see FrateBiasWT) */
  LBias = -KLBIAS * F*log(F*Xtheta) - F*FBias;/* Bias per kb in X (overlapping Y) : In document LBias is larger by Frate = F - 1/Xtheta */
  EBias = 0.5*(Sbias + FBias);/* Bias per aligned end , corresponding to Abias/2 in the document. This is larger by KBIAS/4 compared to the document, so the final
			       score will be KBIAS/2 larger than in the document */
  if(resSD[0] > 0.0) /* PRbias =  Bias per site interval in X for log(Pr(y)) and Sm()-log(pTP) : will be estimated based on average of actual alignments in refalign */
    FBias += PRbias;

  /* NOTE: If biasWTend == 1 : TotalBias(X) = Len * (Frate+LBias) + (M-1)*FBias + 2*EBias = Len * (Frate+LBias) + M*FBias + Sbias (scaled by biasWT)
     where Len is the sum of the intervals of X (above resKB2, if RES_BND), and M is the number of sites in X. This is larger by KBIAS/2 than in the document due to EBias (see above).
     The score of Nanomap X can be adjusted to determine the true -log(LR) :
     -log(LR(X)) = score(X) - TotalBias(X) + (1.0-BETA)*Sum(i) log(Xtheta/sqrt(2*M_PI*(varF+var*y[i]+varR*y[i]*y[i]))),
       where y[i] are the reference intervals aligned with X.
 */

  if(VERB){
    printf("K=%0.2f,KL=%0.2f,KF=%0.2f,G=%0.2f,B=%0.2f,pTP=%0.4f,F=%0.8f,R=%0.6f,Xtheta=%0.8f,Ylambda=%0.4f(Xl=%0.4f),Ygm=%0.4f,Frate=%0.8f,PixelLen= %0.8f,outlierMax=%0.4e\n",
	   //    printf("K=%0.2f,KL=%0.2f,KF=%0.2f,G=%0.2f,B=%0.2f,pTP=%0.8f,F=%0.10f,R=%0.10f,Xtheta=%0.8f,Ylambda=%0.8f(Xl=%0.8f),Ygm=%0.8f,Frate=%0.10f\n",
	   KBIAS,KLBIAS,KFBIAS,GAMMA,BETA,pTP,F,R,Xtheta,Ylambda,Xlambda,Ygm, Frate, PixelLen, outlierMax);
    fflush(stdout);
  }
  raYlambda = Ylambda;

  FnPenalty = -log(pFN);/* missing cut penalty*/
  FpPenalty = -log(F*Xtheta);/* false cut penalty */

  if(GAMMA != 1.0){
    FnPenalty *= GAMMA;
    FpPenalty *= GAMMA;
  }

  if(DEBUG) assert(isfinite(FnPenalty));
  if(DEBUG && !isfinite(FpPenalty)){
    printf("FpPenalty -> %0.8f : F= %0.8f, Xtheta= %0.8f, GAMMA=%0.8f\n",FpPenalty, F, Xtheta, GAMMA);
    fflush(stdout);
    assert(isfinite(FpPenalty));
  }

  MinScore = -(DELTA_Y-1)*FnPenalty-(DELTA_X-1)*FpPenalty;
  MinScore *= 2.0;

  int repeatcnt = 0;

  // Loutlier:

  repeatcnt++;
  if(Poutlier <= 0.0){
    OutlierPenalty = MINSCORE;
    OutlierBias = 0.0;
  } else {
    if(maptype==0 && !OUTLIER_LTYPE){
      printf("WARNING: maptype=0,OUTLIER_LTYPE=0: this combination not implemented (treating as maptype=0, OUTLIER_LTYPE=1)\n");
      fflush(stdout);
    }

    if(outlierLambda >= 1000.0)
      outlierLambda = 1e+30;
    if(maptype){
      if(LAMBDA_FIX1 || outlierLambda >= 1000.0){
	OutlierPenalty = log(Poutlier);/* Outlier Likelihood = 1 (NOTE : this has an unbounded integral!) */
      } else {
	if(OUTLIER_LTYPE == 0){
	  OutlierPenalty = log(Poutlier * Xtheta / outlierLambda);/* Outlier Likelihood = exp(-(x+y)/outlierLambda */
	} else {
	  OutlierPenalty = log(Poutlier * Xtheta / (2.0 * outlierLambda));/* Outlier Likelihood = exp(-|x-y|/outlierLambda */
	}
      }
    } else {/* maptype == 0 */
      if(LAMBDA_FIX0 >= 2 || outlierMax >= 1000.0){
	if(LAMBDA_FIX0 || outlierLambda >= 1000.0)
	  OutlierPenalty = log(Poutlier);/* Outlier Likelihood = 1 (NOTE : this has an unbounded integral, but this is OK since it is a Probability density Ratio)  */
	else
	  OutlierPenalty = log(Poutlier * Xtheta / (2.0 * outlierLambda)); /* Outlier likelihood = exp(-|x-y|/outlierLambda) */
      } else {
	if(LAMBDA_FIX0 || outlierLambda >= 1000.0){
	  OutlierPenalty = log(Poutlier * Xtheta / (2.0 * outlierMax));/* Outlier Likelihood = 1 IFF |x-y| < outlierMax (0 otherwise) */
	} else
	  OutlierPenalty = log(Poutlier * Xtheta / ( 2.0 * outlierLambda * (1.0 - exp(-outlierMax/outlierLambda))));/* Outlier Likelihood = exp(-|x-y|/outlierLambda) IFF |x-y| < outlierMax (0 otherwise) */
      }
    }
    OutlierPenalty -= log(max(0.0001,1.0-Poutlier));
    OutlierBias = FAST_OUTLIER ? 0.0 : log(max(0.0001,1.0-Poutlier));
  }/* Poutlier > 0.0 */
  OutlierLambdaInv = 1.0/outlierLambda;
  OutlierPenalty2 = OutlierPenalty + OutlierBias;

  if(PoutlierEnd <= 0.0){
    OutlierEndPenalty = MINSCORE;
    OutlierEndBias = 0.0;
  } else {
    OutlierEndPenalty = log(PoutlierEnd/max(0.0001,1.0-PoutlierEnd));
    OutlierEndBias = log(max(0.0001,1.0-PoutlierEnd));
  }
  
  biasWToutlierF = biasWToutlier;

#if (USE_SSE==1 && USE_AVX==1 && USE_RFLOAT==1)
  v256_Frate = _mm256_broadcast_ss(&Frate);
  v256_OutlierPenalty = _mm256_broadcast_ss(&OutlierPenalty);

  float micro_unit = 1e-6;
  v256_micro_unit = _mm256_broadcast_ss(&micro_unit);

  float outlierMaxf = outlierMax;
  v256_outlierMax = _mm256_broadcast_ss(&outlierMaxf);

  float outlierNegMaxf = -outlierMax;
  v256_outlierNegMax = _mm256_broadcast_ss(&outlierNegMaxf);

  float OutlierLambdaInvf = OutlierLambdaInv;
  v256_OutlierLambdaInv = _mm256_broadcast_ss(&OutlierLambdaInvf);
  if(VERB>=2){
    float vOutlierLambdaInvf[8];
    _mm256_storeu_ps(vOutlierLambdaInvf, v256_OutlierLambdaInv);    
    printf("OutlierLambdaInv = %0.6f, %0.6f vOutlierLambdaInv[0]= %0.6f\n", OutlierLambdaInv, OutlierLambdaInvf, vOutlierLambdaInvf[0]);
    fflush(stdout);
  }
  v256_biasWToutlierF = _mm256_broadcast_ss(&biasWToutlierF);
  
  float OUTLIER_TYPE_F = OUTLIER_TYPE;
  v256_OUTLIER_TYPE = _mm256_broadcast_ss(&OUTLIER_TYPE_F);

  float OUTLIER_TYPE1_F = OUTLIER_TYPE1;
  v256_OUTLIER_TYPE1 = _mm256_broadcast_ss(&OUTLIER_TYPE1_F);

#ifdef __AVX2__
  v256i_76543210 = _mm256_set_epi32(7,6,5,4,3,2,1,0);
#else
  v128i_7654 = _mm_set_epi32(7,6,5,4);
  v128i_3210 = _mm_set_epi32(3,2,1,0);
#endif

  for(unsigned int mask8 = 0; mask8 < 256; mask8++)
    mask8to256[mask8] =_mm256_castsi256_ps(_mm256_set_epi32(((mask8 >> 7) & 0x1) ? 0xffffffff : 0x0,
							    ((mask8 >> 6) & 0x1) ? 0xffffffff : 0x0,
							    ((mask8 >> 5) & 0x1) ? 0xffffffff : 0x0,
							    ((mask8 >> 4) & 0x1) ? 0xffffffff : 0x0,
							    ((mask8 >> 3) & 0x1) ? 0xffffffff : 0x0,
							    ((mask8 >> 2) & 0x1) ? 0xffffffff : 0x0,
							    ((mask8 >> 1) & 0x1) ? 0xffffffff : 0x0,
							    ((mask8 >> 0) & 0x1) ? 0xffffffff : 0x0));
#endif
  
#if (USE_MIC==1 && USE_RFLOAT==1)
  v512_Frate = _mm512_set1_ps(Frate);
  v512_OutlierPenalty = _mm512_set1_ps(OutlierPenalty);

  float tmp = 1e-6;
  v512_micro_unit = _mm512_set1_ps(tmp);

  tmp = outlierMax;
  v512_outlierMax = _mm512_set1_ps(tmp);
  tmp = -outlierMax;
  v512_outlierNegMax = _mm512_set1_ps(tmp);
  
  tmp = OutlierLambdaInv;
  v512_OutlierLambdaInv = _mm512_set1_ps(tmp);

  tmp = biasWToutlierF;
  v512_biasWToutlierF = _mm512_set1_ps(tmp);  

  tmp = 1.0;
  v512_true = _mm512_set1_ps(tmp);

  tmp = OUTLIER_TYPE;
  v512_OUTLIER_TYPE = _mm512_set1_ps(tmp);
  
  tmp = OUTLIER_TYPE1;
  v512_OUTLIER_TYPE1 = _mm512_set1_ps(tmp);
#endif
  
  //  FpPenaltyBias = FpPenalty - FBias;
  FrateBiasWT = biasWT * (Frate + LBias);
  FBiasWT = biasWT * FBias;
  EBiasWT = biasWT * EBias;

  F2var = 2.0 * varF;
  S2var = 2.0 * var;
  R2var = 2.0 * varR;
  E2var = 2.0 * varE;

  RefPen = log(Xtheta) - 0.5*log(M_PI);

  if(VERB){
    printf("LBias=%0.6f,EBias=%0.4f,FBias=%0.4f(PRbias=%0.4f,Sbias=%0.4f)FnPenalty=%0.6f,FpPenalty=%0.6f,sf=%0.6f,sd=%0.6f,sr=%0.6f,se=%0.6f,res=%0.6f,resSD=%0.6f,minscore=%0.3f,RefPen=%0.6f\n",
	   LBias,EBias,FBias,PRbias,Sbias,FnPenalty,FpPenalty,SF[0],SD[0],SR[0],SE[0],res[0],resSD[0],MinScore,RefPen);
    if(Poutlier > 0.0)
      printf("Poutlier=%0.6f(OutlierBias=%0.6f,Pen=%0.6e),PoutlierEnd=%0.6f(Bias=%0.6f,Pen=%0.6e),maptype=%d:cpu time=%0.6f, wall time=%0.6f\n",
	     Poutlier,OutlierBias,OutlierPenalty,PoutlierEnd, OutlierEndBias, OutlierEndPenalty, maptype, mtime(), wtime());
    if((maptype ? OUTLIER_TYPE1 : OUTLIER_TYPE) && fabs(OutlierPenalty + OutlierBias) < FnPenalty + FpPenalty){
      double maxOutlier = ((LAMBDA_FIX1 || outlierLambda >= 1000.0) ? Xtheta : (OUTLIER_LTYPE ? 2.0 : 1.0) * outlierLambda) * pFN * F;

#if 0
      printf("WARNING : -outlier should be at least as stringent as %0.4e to avoid back to back outliers : changing to -outlier %0.4e\n", maxOutlier,maxOutlier * 0.8);
      fflush(stdout);
      
      Poutlier = maxOutlier * 0.8;
      if(DEBUG) assert(repeatcnt <= 1);
      if(repeatcnt <= 1)
	goto Loutlier;/* repeat computation of outlier penalty terms : WARNING may cause matchgroups to break up near ends due to large outlier */
#else
      printf("WARNING : -outlier should be at least as stringent as %0.4e to avoid back to back outliers\n", maxOutlier);
      fflush(stdout);
#endif
    }
    fflush(stdout);
  }


  if(DEBUG) assert(resSD[0] > 0.0);
  IresSD = 1.0/(resSD[0] * PixelLen * sqrt(2.0));

  int revmult = (MultiMatches && MultiMatchesRev) ? 2 : 1;

  /* allocate PRsum[refid= refstart .. refend-1][I=1..N] */
  if(FN_ADJUST){
    PRsum = new double*[refend * revmult];
    long long Nsum = 0;
    for(int refid = refstart; refid < refend; refid++)
      Nsum += refmap[refid]->numsite[0];
    PRsum[0] = new double[Nsum * revmult];
    long long Ncnt = 0;
    for(int refid = refstart; refid < refend; refid++){
      PRsum[refid] = &PRsum[0][Ncnt];
      Ncnt += refmap[refid]->numsite[0];
      if(MultiMatches && MultiMatchesRev){
	PRsum[refid + refend] = &PRsum[0][Ncnt];
	Ncnt += refmap[refid]->numsite[0];
      }
    }
  }

  /* allocate PRtab[refid= refstart .. refend-1][K= 0..max(max(2,KMAX),DELTA_Y)][I= K+1..N]
              YPen[refid= refstart .. refend-1][K= 0..min(N[c]-1,KMAX)][T= 0..KMAX][I= K+1..N[c]][n= 1..DELTA_Y] */
  PRtab = new Cprtab**[refend * revmult];
  YPen = new CYPen****[refend * revmult];
  score_init_mem += refend*revmult*(sizeof(Cprtab**)+sizeof(CYPen****));

  int PRsiz = max(max(2,KMAX),DELTA_Y);

  /* First compute total memory required */
  long long YPmemPPP = 0, YPmemPP = 0, YPmemP = 0, YPmem = 0;/* size of total memory of type CYPen***,CYPen**,CYPen*,CYPen respectively */
  long long PRmemP = 0, PRmem = 0;/* size of total memory of type Cprtab*,Cprtab respectively */
  CMemCnt *MemCnt = new CMemCnt[refend];/* MemCnt[i=refstart..refend-1] will have start location(s) in memory pool(s) for reference contig i (along with its reverse, revmult==2) */

  for(int refid = refstart; refid < refend; refid++){
    Cmap *Ymap = refmap[refid];
    int N = Ymap->numsite[0];
    int Kmax = min(N-1,KMAX);

    CMemCnt *p = &MemCnt[refid];
    p->PRmemP = PRmemP;
    p->PRmem = PRmem;
    p->YPmemPPP = YPmemPPP;
    p->YPmemPP = YPmemPP;
    p->YPmemP = YPmemP;
    p->YPmem = YPmem;

    PRmemP += (PRsiz+1)*revmult;
    PRmem += (PRsiz+1)*(N+1)*revmult;

    YPmemPPP += (Kmax+1)*revmult;
    int KK = (KMAX+1)*(Kmax+1);
    YPmemPP += KK * revmult;
    YPmemP += KK*(N+1) * revmult;
    YPmem += KK*N*DELTA_Y * revmult;
  }
  if(DEBUG) assert(PRmem >= 0);
  if(DEBUG) assert(YPmem >= 0);
  
  long long memsiz = PRmemP*sizeof(Cprtab*) + PRmem*sizeof(Cprtab); 
  if(VERB/* HERE >=2 */){
    printf("PRmemP=%lld,PRmem=%lld, total memory needed = %lld bytes\n",PRmemP,PRmem,memsiz);
    fflush(stdout);
  }
  cprtab_p_pool = new Cprtab*[PRmemP];
  cprtab_pool = new Cprtab[PRmem];
  score_init_mem += memsiz;

  memsiz = YPmemPPP*sizeof(CYPen***)+YPmemPP*sizeof(CYPen**) + YPmemP*sizeof(CYPen*)+YPmem*sizeof(CYPen);
  if(VERB/* HERE >=2 */){
    printf("YPmemPPP=%lld,YPmemPP=%lld,YPmemP=%lld,YPmem=%lld, total memory needed = %lld bytes (refid=%d..%d,revmult=%d)\n",
	   YPmemPPP,YPmemPP,YPmemP,YPmem,memsiz, refstart,refend-1,revmult);
    fflush(stdout);
  }
  cypen_ppp_pool = new CYPen***[YPmemPPP];
  cypen_pp_pool = new CYPen**[YPmemPP];
  cypen_p_pool = new CYPen*[YPmemP];
  cypen_pool = new CYPen[YPmem];
  score_init_mem += memsiz;

  if(VERB/* HERE >=2 */){
    printf("score_init:Total %lu bytes for score tables succesfully allocated (cypen_pool= x%p)\n",score_init_mem, cypen_pool);
    if(VERB/* HERE HERE >=2 */){
      static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
      getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
      printf("\t VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb\n",  VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9);
      fflush(stdout);
    }
    fflush(stdout);
  }
  RFLOAT aNaN = nan("NaN");
  if(DEBUG) assert(IsNan(aNaN));

  int numthreads = 1;
  #ifdef _OPENMP
  numthreads = omp_get_max_threads();
  #endif
  numthreads = min(numthreads,MaxThreads);

  if(DEBUG>=2){

    #pragma omp parallel for num_threads(numthreads) schedule(static,256)
    for(long long i = 0; i < PRmem; i++){
      Cprtab *p = &cprtab_pool[i];
      p->LogPr = aNaN;
      p->Sm = aNaN;
    }

    #pragma omp parallel for num_threads(numthreads) schedule(static,256)
    for(long long i = 0; i < YPmem; i++){
      CYPen *p = &cypen_pool[i];
      p->Sm = aNaN;
      p->Pen = aNaN;
      p->Gauss = aNaN;
      p->deltaY = aNaN;
      //      assert(IsNan(p->deltaY));
      p->Ivar = aNaN;
    }
  }

  int nthreads = min(refend-refstart,numthreads);

  #pragma omp parallel num_threads(nthreads) if(nthreads > 1)
  {
    FLOAT *Yrev = 0;
    if(MultiMatches && MultiMatchesRev)
      Yrev = new FLOAT[maxN + 2];
  
    #pragma omp for schedule(dynamic,1)
    for(int refid = refstart; refid < refend; refid++){
      Cmap *Ymap = refmap[refid];
      int N = Ymap->numsite[0];
      FLOAT *Y = Ymap->site[0];
      if(MultiMatches && MultiMatchesRev)      
	for(int I=0;I <= N+1;I++)
	  Yrev[I] = Y[N+1] - Y[N+1-I];

      int *Kmax = RKmax[refid];

      double *pPRsum = NULL;
      if(FN_ADJUST){
	/* compute PRsum[refid][1..N] */
	pPRsum = PRsum[refid];

	double PRcum = 0.0;
	pPRsum[1] = PRcum;
    
	for(int I = 2; I <= N; I++){
	  FLOAT y = Y[I] - Y[I-1];
	  pPRsum[I] = PRcum += Pr(y);
	}
      }

      CMemCnt *pm = &MemCnt[refid];
      long long PRcnt = pm->PRmem;

      /* allocate PRtab[refid][K=0..PRsiz][I=0..N] */
      PRtab[refid] = &cprtab_p_pool[pm->PRmemP];

      for(int K = 0; K <= PRsiz; K++){
	int n = K+1;
	PRtab[refid][K] = &cprtab_pool[PRcnt]; PRcnt += N+1;
	Cprtab *pp = &PRtab[refid][K][0];

	if(n <= DELTA_Y){
	  for(int I = K+1; I <= N;I++){
	    Cprtab *p = &pp[I];
	  
	    if(FN_ADJUST){
	      double prsum = pPRsum[I] - pPRsum[I-n] - 1.0;
	      p->FnPenalty = max(0.0, prsum) * FnPenalty;
	    }

	    FLOAT y = Y[I] - Y[I-n];
	
	    if(DEBUG>=2 && !(IsNan(p->LogPr))){
#pragma omp critical 
	      {
		printf("refid=%d,K=%d,n=%d,I=%d(N=%d):p->LogPr=%0.8e, PRtab[refid][K][I].LogPr=%0.8e\n",
		       refid,K,n,I,N,p->LogPr, PRtab[refid][K][I].LogPr);
		fflush(stdout);
		assert(IsNan(p->LogPr));
	      }
	    }
	    p->LogPr = log(Pr(y));
	    if(DEBUG) assert(isfinite(p->LogPr));/* could fail if Pr(y) == 0 */
	  }
	}

	if(K <= KMAX){
	  for(int I = K+1; I <= N;I++){
	    if(K <= Kmax[I]){
	      Cprtab *p = &pp[I];
			    
	      if(DEBUG>=2) assert(IsNan(p->Sm));
	      p->Sm = Sm(0,I,K,Y);
	      if(DEBUG>=2) assert(!IsNan(p->Sm) && p->Sm < 0.0);
	    }
	  }
	}
      }

      if(MultiMatches && MultiMatchesRev){
        int *Kmax = RKmax[refid + refend];

        if(FN_ADJUST){ /* compute PRsum[refid+refend][1..N] */
          pPRsum = PRsum[refid + refend];
	  double PRcum = 0.0;
	  pPRsum[1] = PRcum;
	  for(int I = 2; I <= N; I++){
            FLOAT y = Yrev[I] - Yrev[I-1];
	    pPRsum[I] = PRcum += Pr(y);
          }
        }

	/* allocate PRtab[refid+refidend][K=0..PRsiz][I=0..N] */
	PRtab[refid + refend] = &cprtab_p_pool[pm->PRmemP + (PRsiz+1)];
	
	for(int K = 0; K <= PRsiz; K++){
          int n = K+1;
	  PRtab[refid+refend][K] = &cprtab_pool[PRcnt]; PRcnt += N+1;
	  Cprtab *pp = &PRtab[refid+refend][K][0];

	  if(n <= DELTA_Y){
            for(int I = K+1; I <= N;I++){
              Cprtab *p = &pp[I];
	  
	      if(FN_ADJUST){
                double prsum = pPRsum[I] - pPRsum[I-n] - 1.0;
		p->FnPenalty = max(0.0, prsum) * FnPenalty;
              }

	      FLOAT y = Yrev[I] - Yrev[I-n];
	
	      if(DEBUG>=2 && !(IsNan(p->LogPr))){
#pragma omp critical 
                {
                  printf("refid=%d,K=%d,n=%d,I=%d(N=%d):p->LogPr=%0.8e, PRtab[refid+refend][K][I].LogPr=%0.8e\n",
                    refid,K,n,I,N,p->LogPr, PRtab[refid+refend][K][I].LogPr);
		  fflush(stdout);
		  assert(IsNan(p->LogPr));
                }
              }   
	      p->LogPr = log(Pr(y));
	      if(DEBUG) assert(isfinite(p->LogPr));/* could fail if Pr(y) == 0 */
            }
          }

	  if(K <= KMAX){
            for(int I = K+1; I <= N;I++){
              if(K <= Kmax[I]){
                Cprtab *p = &pp[I];
			    
		if(DEBUG>=2) assert(IsNan(p->Sm));
		p->Sm = Sm(0,I,K,Yrev);
		if(DEBUG>=2) assert(!IsNan(p->Sm) && p->Sm < 0.0);
              }
            }
          }
        }
      }

      int kmax = min(N-1,KMAX);
      long long YPcntPP = pm->YPmemPP, YPcntP = pm->YPmemP, YPcnt = pm->YPmem;
      long long YPcntMax = (refid < refend-1) ? pm[1].YPmem : YPmem;
      if(DEBUG) assert(YPcntMax  <= YPmem);

      /* allocate YPen[refid][K=0..min(N-1,KMAX)][T=0..KMAX][I=1..N][n=1..DELTA_Y] */
      YPen[refid] = &cypen_ppp_pool[pm->YPmemPPP];

      for(int K = kmax; K >= 0; K--){
	if(DEBUG) assert(N > K);
	YPen[refid][K] = &cypen_pp_pool[YPcntPP]; YPcntPP += KMAX+1;
	for(int T = KMAX; T >= 0; T--){
	  YPen[refid][K][T] = &cypen_p_pool[YPcntP]; YPcntP += N+1;
	  CYPen *pmem = &cypen_pool[YPcnt]; YPcnt += N*DELTA_Y;
	  if(DEBUG>=1+RELEASE) assert(YPcnt <= YPcntMax);

	  for(int I = N; I > K; I--){
	    CYPen *p = YPen[refid][K][T][I] = &pmem[-1]; pmem += DELTA_Y;
	    if(K > Kmax[I])
	      continue;
	    FLOAT Yik = Yc(Y,I,K);
	    FLOAT sm = Sm(0,I,K,Y);
	    for(int n = min(I-K-T-1,DELTA_Y); n > 0; n--){
              if(DEBUG>=1+RELEASE) assert(0 < I-K-n);
	      if (T > Kmax[I-K-n])
		continue;
	      FLOAT Ynt = Yc(Y,I-K-n,T);
	      FLOAT y = Yik - Ynt;
	      FLOAT var = F2var + S2var*y;
	      if(QUADRATIC_VARIANCE)
		var += R2var * y * y;
	      if(RES_VARIANCE){
		FLOAT resR = Y[I] - Y[I-K];
		FLOAT resL = Y[I-K-n] - Y[I-K-n-T];
		FLOAT resvar = resL*resL + resR*resR;
		var += E2var*resvar;
	      }
	      if(DEBUG>=1+RELEASE && !(var > 0.0)){
#pragma omp critical
		{
		  printf("refid=%d,N=%d:K=%d,T=%d,I=%d,n=%d:Yik=%0.6f,Ynt=%0.6f,y=%0.6f,var=%0.8e(F2var=%0.6e,S2var=%0.6e,R2var=%0.6e,E2var=%0.6e)\n",
			 refid,N,K,T,I,n,Yik,Ynt,y,var*0.5,F2var,S2var,R2var,E2var);
		  fflush(stdout);
		  assert(var > 0.0);
		}
	      }
	      if(DEBUG>=1+RELEASE) assert(fabs(PRtab[refid][n-1][I-K].LogPr - log(Pr(Y[I-K]-Y[I-K-n]))) < 1.0e-7);
	      if(DEBUG>=1+RELEASE) assert(&p[n] < pmem);
	      if(DEBUG>=1+RELEASE) assert(cypen_pool <= &p[n] && &p[n] < &cypen_pool[YPmem]);
	      if(DEBUG>=2 && !IsNan(p[n].deltaY)){
                #pragma omp critical
                {
		  printf("refid=%d,N=%d:K=%d,T=%d,I=%d,n=%d:Yik=%0.6f,Ynt=%0.6f,y=%0.6f,var=%0.8e(F2var=%0.6e,S2var=%0.6e,R2var=%0.6e,E2var=%0.6e)\n",
			 refid,N,K,T,I,n,Yik,Ynt,y,var*0.5,F2var,S2var,R2var,E2var);
		  printf("p[n].deltaY= %0.6e, aNaN= %0.6e\n", p[n].deltaY, aNaN);
		  fflush(stdout);
                  assert(IsNan(p[n].deltaY));
                }
              }
	      if(DEBUG>=2) assert(IsNan(p[n].Pen));
	      if(DEBUG>=2) assert(IsNan(p[n].Sm));
	      if(DEBUG>=2) assert(IsNan(p[n].Ivar));
	      p[n].deltaY = Yik - Ynt;
	      Cprtab *pr = &PRtab[refid][n-1][I-K];
	      p[n].Pen = pr->LogPr - (FN_ADJUST ? pr->FnPenalty : (n-1)*FnPenalty) + FpPenalty;
	      p[n].Gauss = RefPen - 0.5*mylog(var);
	      p[n].Sm = sm;
	      p[n].Ivar = 1.0/var;

	      if(VERB>=2 && refid==0 && I==27151 && K==0 && T==0 && n==1){
#pragma omp critical
		{
		  FLOAT resR = Y[I] - Y[I-K];
		  FLOAT resL = Y[I-K-n] - Y[I-K-n-T];
		  FLOAT resvar = resL*resL + resR*resR;

		  printf("YPen[refid=%d][K=%d][T=%d][I=%d][n=%d]:Pen=%0.6f,Gauss=%0.6f,sm=%0.6f,var=%0.8e,Ivar=%0.8e(F2var=%0.8e,S2var=%0.8e,R2var=%0.8e,E2var=%0.8e,y=%0.8f,resR=%0.6f,resL=%0.6f,resvar=%0.8e)\n",
			 refid,K,T,I,n,p[n].Pen,p[n].Gauss,p[n].Sm,var*0.5,p[n].Ivar,F2var,S2var,R2var,E2var,y,resR,resL,resvar);
		  fflush(stdout);
		}
	      }
	    }
	  }
	  if(DEBUG>=1+RELEASE) assert(pmem <= &cypen_pool[YPcnt]);
	}
      }

      if(MultiMatches && MultiMatchesRev){

      int *Kmax = RKmax[refid + refend];

      /* allocate YPen[refid+refend][K=0..min(N-1,KMAX)][T=0..KMAX][I=1..N][n=1..DELTA_Y] */
      YPen[refid+refend] = &cypen_ppp_pool[pm->YPmemPPP + (kmax+1)];

      for(int K = kmax; K >= 0; K--){
	if(DEBUG) assert(N > K);
	YPen[refid+refend][K] = &cypen_pp_pool[YPcntPP]; YPcntPP += KMAX+1;
	for(int T = KMAX; T >= 0; T--){
	  YPen[refid+refend][K][T] = &cypen_p_pool[YPcntP]; YPcntP += N+1;
	  CYPen *pmem = &cypen_pool[YPcnt]; YPcnt += N*DELTA_Y;
	  if(DEBUG>=1+RELEASE) assert(YPcnt <= YPcntMax);

	  for(int I = N; I > K; I--){
	    CYPen *p = YPen[refid+refend][K][T][I] = &pmem[-1];   pmem += DELTA_Y;
	    if(K > Kmax[I])
	      continue;
	    FLOAT Yik = Yc(Yrev,I,K);
	    FLOAT sm = Sm(0,I,K,Yrev);
	    for(int n = min(I-K-T-1,DELTA_Y); n > 0; n--){
              if(DEBUG>=1+RELEASE) assert(0 < I-K-n);
	      if (T > Kmax[I-K-n])
		continue;
	      FLOAT Ynt = Yc(Yrev,I-K-n,T);
	      FLOAT y = Yik - Ynt;
	      FLOAT var = F2var + S2var*y;
	      if(QUADRATIC_VARIANCE)
		var += R2var * y * y;
	      if(RES_VARIANCE){
		FLOAT resR = Yrev[I] - Yrev[I-K];
		FLOAT resL = Yrev[I-K-n] - Yrev[I-K-n-T];
		FLOAT resvar = resL*resL + resR*resR;
		var += E2var*resvar;
	      }
	      if(DEBUG>=1+RELEASE && !(var > 0.0)){
#pragma omp critical
		{
		  printf("refid+refend=%d,N=%d:K=%d,T=%d,I=%d,n=%d:Yik=%0.6f,Ynt=%0.6f,y=%0.6f,var=%0.8e(F2var=%0.6e,S2var=%0.6e,R2var=%0.6e,E2var=%0.6e)\n",
			 refid+refend,N,K,T,I,n,Yik,Ynt,y,var*0.5,F2var,S2var,R2var,E2var);
		  fflush(stdout);
		  assert(var > 0.0);
		}
	      }
	      if(DEBUG>=1+RELEASE) assert(fabs(PRtab[refid+refend][n-1][I-K].LogPr - log(Pr(Yrev[I-K]-Yrev[I-K-n]))) < 1.0e-7);
	      if(DEBUG>=1+RELEASE) assert(&p[n] < pmem);
	      if(DEBUG>=1+RELEASE) assert(&p[n] < &cypen_pool[YPmem]);
	      if(DEBUG>=2) assert(IsNan(p[n].deltaY));
	      if(DEBUG>=2) assert(IsNan(p[n].Pen));
	      if(DEBUG>=2) assert(IsNan(p[n].Sm));
	      if(DEBUG>=2) assert(IsNan(p[n].Ivar));
	      p[n].deltaY = Yik - Ynt;
	      Cprtab *pr = &PRtab[refid+refend][n-1][I-K];
	      p[n].Pen = pr->LogPr - (FN_ADJUST ? pr->FnPenalty : (n-1)*FnPenalty) + FpPenalty;
	      p[n].Gauss = RefPen - 0.5*mylog(var);
	      p[n].Sm = sm;
	      p[n].Ivar = 1.0/var;

	      if(VERB>=2 && refid==0 && I==27151 && K==0 && T==0 && n==1){
#pragma omp critical
		{
		  FLOAT resR = Yrev[I] - Yrev[I-K];
		  FLOAT resL = Yrev[I-K-n] - Yrev[I-K-n-T];
		  FLOAT resvar = resL*resL + resR*resR;

		  printf("YPen[refid=%d+refend][K=%d][T=%d][I=%d][n=%d]:Pen=%0.6f,Gauss=%0.6f,sm=%0.6f,var=%0.8e,Ivar=%0.8e(F2var=%0.8e,S2var=%0.8e,R2var=%0.8e,E2var=%0.8e,y=%0.8f,resR=%0.6f,resL=%0.6f,resvar=%0.8e)\n",
			 refid,K,T,I,n,p[n].Pen,p[n].Gauss,p[n].Sm,var*0.5,p[n].Ivar,F2var,S2var,R2var,E2var,y,resR,resL,resvar);
		  fflush(stdout);
		}
	      }
	    }
	  }
	  if(DEBUG>=1+RELEASE) assert(pmem <= &cypen_pool[YPcnt]);
	}
      }

      }
      if(DEBUG){
	if(refid < refend-1){
	  assert(PRcnt <= pm[1].PRmem);
	  assert(YPcntPP <= pm[1].YPmemPP);
	  assert(YPcntP <= pm[1].YPmemP);
	  if(!(YPcnt <= pm[1].YPmem)){
	    printf("refid=%d(refstart=%d,refend=%d):YPcnt=%lld,MemCnt[refid].YPmem=%lld,MemCnt[refid+1].YPmem=%lld,N=%d,kmax=%d,KMAX=%d\n",
		   refid,refstart,refend,YPcnt,MemCnt[refid].YPmem,MemCnt[refid+1].YPmem,N,kmax,KMAX);
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
    }// for int refid = refstart .. refend - 1

    delete [] Yrev;
  }// omp parallel

  delete [] MemCnt;

  long long maxMN = max(maxN,maxM);
  OutlierPenaltyBC = new RFLOAT[maxMN];
  if(!outlierBC)
    memset(OutlierPenaltyBC, 0, maxMN*sizeof(RFLOAT));
  else {
    OutlierPenaltyBC[0] = 0.0;

    int nthreads = min((int)(maxMN/512), numthreads);

    #pragma omp parallel for num_threads(nthreads) schedule(static, 512) if(nthreads > 1)
    for(int i = 1; i < maxMN; i++)
      OutlierPenaltyBC[i] = -log((double)i);
  }
  if(VERB>=2){
    printf("score_init():maxM=%lld,maxN=%lld,OutlierPenaltyBC=%p,outlierBC=%d\n",maxM,maxN,OutlierPenaltyBC,outlierBC);
    static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
    getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
    printf("\t VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb:wall= %0.6f secs\n",  VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9,wtime());
    fflush(stdout);
  }
  return 0;
}

/* perform score term updates due to changes in FP[0] and FN[0] : only update those terms used by logLR(). Assumes SCORE_Y==0 */
static void score_update(int SNRadjust = 0)
{
  resKB = res[0] * PixelLen;
  resKB2 = 0.0;

  double pFN = FN[0];
  pTP = pTPd = 1.0 - pFN;
  log_pTP = log_pTPd = log(pTPd);
  log_qTPd = log(pFN);
  double F = FP[0]*0.01;
  double var = fabs(SD[0])*SD[0];
  double varF = SF[0]*SF[0];
  double varR = SR[0]*SR[0];

  Xlambda = pTP/(1.0/(Xtheta - var*0.5) - F);/* HERE : does not include correction for SR[],SE[] */
  if(XmapCount[0] /* WAS169 && RefRepeats==1 */ && !SNRadjust){
    Ylambda = Xlambda;
    Ygm = XmapGtheta[0] * /* NEW12 */ (PixelLen / origPixelLen) * Xlambda/Xtheta;
  } else if(SCORE_Y || SNRadjust){ /* make score dependent only on Y (not X) : for debugging only or when adjusting minSNR so scoring function does not depend on map label density */
    Xlambda = Ylambda;
    Xtheta = 1.0/(F + pTP/Xlambda) + var*0.5;
  }

  if(DEBUG){
    assert(FN[0] > 0.0);
    assert(FN[0] < 1.0);
    assert(Ylambda > 0.0);
    assert(Xlambda > 0.0);
  }

  double R = RTHETA_FIX ? Rnummaps/theta : 0.0;/* should be reduced by inverse of average linear chromosome size */

  Frate = F + R - thetaScale/Xtheta; // WAS F - 1.0/Xtheta

  double YgmBpTP = (Ygm/pTP) * (Xlambda/Ylambda);
  double Sbias = log(Xtheta) - 0.5*log(2.0*M_PI*(varF+var*YgmBpTP+varR*YgmBpTP*YgmBpTP));
  FBias = 0.5*KBIAS - KFBIAS*pFN*log(pFN)/pTP - log_pTP - Sbias; /* Bias per site interval in X (overlapping Y)  = Abias - Sbias in document */
  LBias = -KLBIAS * F*log(F*Xtheta) - F*FBias;/* Bias per kb in X (overlapping Y) : In document LBias is larger by Frate = F - 1/Xtheta */
  EBias = 0.5*(Sbias + FBias);/* Bias per aligned end , corresponding to Abias/2 in the document. This is larger by KBIAS/4 compared to the document, so the final
				 score will be KBIAS/2 larger than in the document */
  if(resSD[0] > 0.0) /* PRbias =  Bias per site interval in X for log(Pr(y)) and Sm()-log(pTP) : will be estimated based on average of actual alignments in refalign */
    FBias += PRbias;
    
  FnPenalty = -log(pFN);/* missing cut penalty*/
  FpPenalty = -log(F*Xtheta);/* false cut penalty */
  if(GAMMA != 1.0){
    FnPenalty *= GAMMA;
    FpPenalty *= GAMMA;
  }

  if(DEBUG) assert(isfinite(FnPenalty));
  if(DEBUG && !isfinite(FpPenalty)){
    printf("FpPenalty -> %0.8f : F= %0.8f, Xtheta= %0.8f, GAMMA=%0.8f\n",FpPenalty, F, Xtheta, GAMMA);
    fflush(stdout);
    assert(isfinite(FpPenalty));
  }
}

static void XPen_alloc(CXPen *XPen)
{
  /* allocate XPen-><field>[m=1..DELTA_X][J=m..maxM+1] */
  if(DEBUG>=2 && !(XPen->deltaX == 0)){
    #pragma omp critical
    {
      assert(XPen->deltaX == 0);
    }
  }
  int M = maxM+2;
  if(PADDING > sizeof(RFLOAT))
    M = (M + (PADDING/sizeof(RFLOAT))-1) & (~((PADDING/sizeof(RFLOAT))-1));/* round up M, so M*sizeof(RFLOAT) is a multiple of PADDING */
  size_t siz = DELTA_X * 4 * M * sizeof(RFLOAT) + 3*PADDING;
  XPen->mem_pool_size = DELTA_X * 5 * sizeof(RFLOAT*) + siz;
  if(LOAD_ALIGNED && PADDING > sizeof(RFLOAT)){
    POSIX_MEMALIGN(XPen->mem_pool, PADDING, XPen->mem_pool_size);
    if(DEBUG) assert(((size_t)XPen->mem_pool & (size_t)(PADDING-1)) == 0);
  } else
    XPen->mem_pool = (char *) malloc(XPen->mem_pool_size);


  XPen->deltaX = (RFLOAT **) &XPen->mem_pool[siz];
  XPen->deltaX--;// so 1st index runs m=1 .. DELTA_X (+1 to support prefetch)
  XPen->Pen = &XPen->deltaX[DELTA_X+1];
  XPen->Bias = &XPen->deltaX[(DELTA_X+1)*2];
  XPen->PenBias = &XPen->deltaX[(DELTA_X+1)*3];

  size_t offset = PADDING;

#pragma ivdep
  for(int m = 1; m <= DELTA_X; m++){
    XPen->deltaX[m] = (RFLOAT *)&XPen->mem_pool[offset];
    XPen->Pen[m] = &XPen->deltaX[m][M];
    XPen->Bias[m] = &XPen->deltaX[m][M*2];
    XPen->PenBias[m] = &XPen->deltaX[m][M*3];
    if(DEBUG && LOAD_ALIGNED && PADDING > sizeof(RFLOAT)){
      if(DEBUG) assert((offset & (size_t)(PADDING-1)) == 0);
      if(DEBUG) assert(((size_t)XPen->deltaX[m] & (size_t)(PADDING-1)) == 0);
      if(DEBUG) assert(((size_t)XPen->Pen[m] & (size_t)(PADDING-1)) == 0);
      if(DEBUG) assert(((size_t)XPen->Bias[m] & (size_t)(PADDING-1)) == 0);
      if(DEBUG) assert(((size_t)XPen->PenBias[m] & (size_t)(PADDING-1)) == 0);
    }
    offset += M * 4 * sizeof(RFLOAT);
    if(m==DELTA_X){/* to support prefetch */
      XPen->deltaX[m+1] = XPen->deltaX[m];
      XPen->Pen[m+1] = XPen->Pen[m];
      XPen->Bias[m+1] = XPen->Bias[m];
      XPen->PenBias[m+1] = XPen->PenBias[m];
    }
  }
  offset += 2*PADDING + DELTA_X * 5 * sizeof(RFLOAT*);

  if(DEBUG/* HERE >=2 */) assert(offset == XPen->mem_pool_size);
}

static void XPen_free(CXPen *XPen)
{
  free(XPen->mem_pool);
  XPen->mem_pool = NULL;
  XPen->mem_pool_size = 0;
  if(DEBUG>=2) { XPen->deltaX = XPen->Pen = XPen->Bias = XPen->PenBias = 0;}
}

static void XPen_init(FLOAT *X, int M, CXPen *XPen)
{
  if(DEBUG) assert(XPen->deltaX != 0 && M <= maxM);

  /* initialize XPen-><field>[m=1..DELTA_X][J=m..M+1] */
  for(int m = 1; m <= DELTA_X; m++){
    RFLOAT *deltaXm = XPen->deltaX[m];
    RFLOAT *Penm = XPen->Pen[m];
    RFLOAT *Biasm = XPen->Bias[m];
    RFLOAT *PenBiasm = XPen->PenBias[m];

#pragma ivdep
    for(int J = m; J <= M+1; J++){
      FLOAT deltaX = X[J] - X[J-m];
      RFLOAT Pen = m * FpPenalty;// WAS + deltaX * Frate;
      RFLOAT Bias = deltaX * FrateBiasWT + m * FBiasWT;
      deltaXm[J] = deltaX;
      Penm[J] = Pen;
      Biasm[J] = Bias;
      PenBiasm[J] = Bias - Pen;
    }
  }
}


/* NOTE : this is a slower version of Sint() that does not use tables that depend of res or resSD and always uses FLOAT AND 
   forces the interval to be scored as an outlier (or NOT) based on "outlier" AND assumes maptype==0
   Only used for debugging noise parameter estimation with REFDEBUG_STRICT >= 1 
 */

static inline double Sint(FLOAT X,
                          FLOAT Y,
			  int m,
                          int n,
			  int J,
			  int I,
			  int K,/**< right end of interval */
			  int T, /**< left end of interval : only used with RES_VARIANCE */
			  FLOAT *Ya,/**< reference map array Ya[0..N+1] */
			  int outlier)/* force interval to be scored as outlier (or NOT) */
{
  if(DEBUG) assert(REFDEBUG_STRICT >= 1);
  if(DEBUG>=2) assert(maptype==0);

  if(outlierLambda < 1000.0 || outlierMax < 1000.0){
    printf("Sint() with REFDEBUG_STRICT not yet implemented with -outlierMax or -outlierLambda\n");
    fflush(stdout);exit(1);
  }

  if(DEBUG>=2 && !(outlierExtend && !outlierExtendLimit) && !(n>=1 && n <= max(max(2,KMAX), outlierExtend ? outlierExtendLimit : DELTA_Y))){
    printf("Sint:X=%0.4f,Y=%0.4f,m=%d,n=%d,J=%d,I=%d,K=%d,T=%d,KMAX=%d,DELTA_Y=%d,outlierExtend=%d,outlierExtendLimit=%d\n",X,Y,m,n,J,I,K,T,KMAX,DELTA_Y,outlierExtend,outlierExtendLimit);
    fflush(stdout);
    assert(n>=1 && n <= max(max(2,KMAX),outlierExtend ? outlierExtendLimit : DELTA_Y));
  }
  if(DEBUG>=2) assert(0 <= K && K <= KMAX);
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= KMAX);
    assert(I-K-n-T >= 1);
  }

  FLOAT err = X-Y;
  FLOAT var;
  if(RES_VARIANCE){
    FLOAT resR = Ya[I] - Ya[I-K];
    FLOAT resL = Ya[I-K-n] - Ya[I-K-n-T];
    FLOAT resvar = resL*resL + resR*resR;
    var = QUADRATIC_VARIANCE ? (F2var + S2var*Y + R2var*Y*Y + E2var*resvar) : (F2var + S2var*Y + E2var*resvar);
  } else
    var = QUADRATIC_VARIANCE ? (F2var + S2var*Y + R2var*Y*Y) : (F2var + S2var*Y);
  if(DEBUG>=2) assert(var > 0.0);
  FLOAT Bias = X * FrateBiasWT + m * FBiasWT;
  FLOAT Gauss = RefPen - 0.5*mylog(var) - err*err/var;
  FLOAT Pen = log(Pr(Ya[I-K]-Ya[I-K-n])) - (m-1)*FpPenalty; // WAS : - X*Frate;
  if(FRATE_FIX0)
    Gauss -= X*Frate;
  else
    Pen -= X*Frate;
  if(FN_ADJUST){
    double sum = 0.0;
    for(int i = 1; i <= n; i++)
      sum += Pr(Ya[I-i+1]-Ya[I-i]);
    Pen -= max(0.0,sum-1.0)*FnPenalty;
  } else
    Pen -= (n-1)*FnPenalty;

  FLOAT sm = Sm(0,I,K,Ya);
  if(outlier){
    FLOAT OutPen = OutlierPenalty;
    if(outlierBC)
      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];

#if 0 // OLD CODE
    if(OUTLIER_TYPE == 0)
      return OutlierBias + Bias * biasWToutlier + OutPen + sm;
    else
      return OutlierBias + Pen + Bias * biasWToutlier + OutPen + sm;
#else // NEW CODE
    return OutlierBias + Pen * OUTLIER_TYPE + Bias * biasWToutlier + OutPen + sm;
#endif
  } else
    return Bias + Pen + Gauss + sm;
}

/* NOTE : this is a slower version of Sint() that does not use tables that depend of res or resSD and always uses FLOAT */

static inline double Sint(FLOAT X,
			  FLOAT Y,
			  int m,
			  int n,
			  int J,
			  int I,
			  int K,/**< right end of interval */
			  int T, /**< left end of interval : only used with RES_VARIANCE */
			  FLOAT *Ya/**< reference map array Ya[0..N+1] */
                          /* , bool verb = false*/)
  
{
  if(DEBUG>=2 && !(outlierExtend && !outlierExtendLimit) && !(n>=1 && n <= max(max(2,KMAX), outlierExtend ? outlierExtendLimit : DELTA_Y))){
    printf("Sint:X=%0.4f,Y=%0.4f,m=%d,n=%d,J=%d,I=%d,K=%d,T=%d,KMAX=%d,DELTA_Y=%d,outlierExtend=%d,outlierExtendLimit=%d\n",X,Y,m,n,J,I,K,T,KMAX,DELTA_Y,outlierExtend,outlierExtendLimit);
    fflush(stdout);
    assert(n>=1 && n <= max(max(2,KMAX),outlierExtend ? outlierExtendLimit : DELTA_Y));
  }
  if(DEBUG>=2) assert(0 <= K && K <= KMAX);
  if(DEBUG>=2 && !(I-K-n >= 1)){
    printf("Sint:X=%0.4f,Y=%0.4f,m=%d,n=%d,J=%d,I=%d,K=%d,T=%d,KMAX=%d\n",X,Y,m,n,J,I,K,T,KMAX);
    fflush(stdout);
    assert(I-K-n >= 1);
  }
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= KMAX);
    assert(I-K-n-T >= 1);
  }

  FLOAT err = X-Y;
  FLOAT var = F2var + S2var*Y;
  if(QUADRATIC_VARIANCE)
    var += R2var * Y * Y;
  if(RES_VARIANCE){
    FLOAT resR = Ya[I] - Ya[I-K];
    FLOAT resL = Ya[I-K-n] - Ya[I-K-n-T];
    FLOAT resvar = resL*resL + resR*resR;
    var += E2var*resvar;
  }
  if(DEBUG>=2) assert(var > 0.0);

  FLOAT Bias = X * FrateBiasWT + m * FBiasWT;
  FLOAT Gauss = RefPen - 0.5*mylog(var) - err*err/var;
  FLOAT Pen =  log(Pr(Ya[I-K]-Ya[I-K-n])) - (m-1)*FpPenalty;// WAS : - X*Frate
  if(FN_ADJUST){
    double sum = 0.0;
    for(int i = 1; i <= n; i++)
      sum += Pr(Ya[I-i+1]-Ya[I-i]);
    Pen -= max(0.0,sum-1.0)*FnPenalty;
  } else
    Pen -= (n-1)*FnPenalty;

  FLOAT PenSm = Sm(0,I,K,Ya);

  if(maptype){
    if(FRATE_FIX1)
      Gauss -= X*Frate;
    else
      Pen -= X*Frate;
    FLOAT OutPen = OutlierPenalty;
    if(outlierBC)
      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
    if(OUTLIER_LTYPE==0)
      OutPen -= (X+Y) * OutlierLambdaInv;
    else
      OutPen -= fabs(err) * OutlierLambdaInv;
#if 0 // OLD CODE
    if(OUTLIER_TYPE1 == 0)
      return OutlierBias + max(Bias + Pen + Gauss, Bias * biasWToutlier + OutPen) + PenSm;
    else
      return OutlierBias + Pen + max(Bias + Gauss, Bias * biasWToutlier + OutPen) + PenSm;
#else // NEW CODE
      return OutlierBias + max(Bias + Pen + Gauss, Bias * biasWToutlier + Pen*OUTLIER_TYPE1 + OutPen) + PenSm;
#endif
  } else {// maptype == 0
    if(FRATE_FIX0)
      Gauss -= X*Frate;
    else
      Pen -= X*Frate;
    if(OUTLIER_DELTA(err)){
      FLOAT OutPen = OutlierPenalty;
      if(outlierBC)
	OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
      OutPen -= fabs(err) * OutlierLambdaInv;
#if 0 // OLD CODE
      if(OUTLIER_TYPE == 0)
	return OutlierBias + max(Bias + Pen + Gauss, Bias * biasWToutlier + OutPen) + PenSm;
      else
	return OutlierBias + Pen + max(Bias + Gauss, Bias * biasWToutlier + OutPen) + PenSm;
#else // NEW CODE
	return OutlierBias + max(Bias + Pen + Gauss, Bias * biasWToutlier + Pen*OUTLIER_TYPE + OutPen) + PenSm;
#endif
    } else {
      /*      if(verb){
	printf("sint():Bias=%0.6f,Pen=%0.6f,Gauss=%0.6f,X*Frate=%0.6f,PenSm=%0.6f:err=%0.6f,var=%0.8f\n",Bias,Pen,Gauss,X*Frate,PenSm,err,var);
	fflush(stdout);
	} */
      return Bias + Pen + Gauss + PenSm;
    }
  }
}

/* outlier penalty terms excluding OutlierPenalty or OutlierBias or PenSm : assumes interval is an outlier */
static inline RFLOAT OutPenExtend(RFLOAT X, RFLOAT Y, int m, int n, int I, int K, Cprtab **PRtab, int refid)
{
  RFLOAT OutPen = outlierBC ? (OutlierPenaltyBC[m] + OutlierPenaltyBC[n]) : 0.0;
  RFLOAT Bias = X * FrateBiasWT + m * FBiasWT;
  if(maptype){
    if(OUTLIER_LTYPE==0)
      OutPen -= (X+Y) * OutlierLambdaInv;
    else
      OutPen -= fabs(X-Y) * OutlierLambdaInv;
    RFLOAT Pen =  (n > DELTA_Y ? 0.0 /* log(Pr(Ya[I-K]-Ya[I-K-n])) */ : PRtab[n-1][I-K].LogPr) - (FN_ADJUST ? max(0.0,PRsum[refid][I-K]-PRsum[refid][I-K-n]-1.0) : (n-1)) * FnPenalty - (m-1)*FpPenalty;
    if(!FRATE_FIX1)
      Pen -= X * Frate;
    return Pen*OUTLIER_TYPE1 + Bias * biasWToutlierF + OutPen;
  } else {
    OutPen -= fabs(X-Y) * OutlierLambdaInv;
    RFLOAT Pen =  (n > DELTA_Y ? 0.0 /* log(Pr(Ya[I-K]-Ya[I-K-n])) */ : PRtab[n-1][I-K].LogPr) - (FN_ADJUST ? max(0.0,PRsum[refid][I-K]-PRsum[refid][I-K-n]-1.0) : (n-1)) * FnPenalty - (m-1)*FpPenalty;
    if(!FRATE_FIX0)
      Pen -= X * Frate;
    return Pen*OUTLIER_TYPE + Bias * biasWToutlierF + OutPen;
  }
}

/** Score of internal aligned interval with resSD>0 : includes Sm(J,I,K) */

/** Sint() can be computed faster in two steps:
  RFLOAT Y,Ivar,GaussY,Sm;
  RFLOAT PenY = SintY(n,I,K,Y,Ivar,GaussY,Sm,YPen[refid]);
  RFLOAT score = SintX(X,Y,PenY,GaussY,Sm,Ivar,XPen[mapid].Pen[m][J],XPen[mapid].Bias[m][J],XPen[mapid].PenBias[m][J],OutlierPenaltyBC[m] + OutlierPenaltyBC[n]);
*/

static inline RFLOAT Sint(RFLOAT X,
			  RFLOAT Y,
			  int m,
			  int n,
			  int J,
			  int I,
			  int K,/**< right end of interval! */
			  int T,/**< left end of interval : only used with RES_VARIANCE */
			  Cprtab **PRtab,
			  FLOAT *Ya)/**< reference map array Ya[0..N+1] */
  
{
  if(DEBUG>=2 && !(n>=1 && n <= max(max(2,KMAX),DELTA_Y))){
    printf("Sint(X=%0.3f,Y=%0.3f,m=%d,n=%d,J=%d,I=%d,K=%d):Y[I]=%0.3f,Y[I-K]=%0.3f:DELTA_Y=%d\n",
	   X,Y,m,n,J,I,K,Ya[I],Ya[I-K],DELTA_Y);
    fflush(stdout);

    assert(n>=1 && n <= max(max(2,KMAX),DELTA_Y));
  }
  if(DEBUG>=2) assert(0 <= K && K <= KMAX);
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= KMAX);
    assert(I-K-n-T >= 1);
  }
  if(DEBUG>=2) assert(fabs(PRtab[n-1][I-K].LogPr - log(Pr(Ya[I-K]-Ya[I-K-n]))) <= 1.0e-10);
  if(DEBUG>=2){
    RFLOAT sm = Sm(0,I,K,Ya);
    if(!(fabs(PRtab[K][I].Sm - sm) <= 1.0e-10)){
      printf("Sint(X=%0.3f,Y=%0.3f,m=%d,n=%d,J=%d,I=%d,K=%d):Y[I]=%0.3f,Y[I-K]=%0.3f:PRtab[K][I].Sm=%0.11f,Sm(0,I,K,Ya)=%0.11f(err=%0.11f)\n",
	     X,Y,m,n,J,I,K,Ya[I],Ya[I-K],PRtab[K][I].Sm,sm,PRtab[K][I].Sm - sm);
      fflush(stdout);
      assert(fabs(PRtab[K][I].Sm - sm) <= 1.0e-10);
    }
  }

  RFLOAT err = X-Y;
  RFLOAT var = F2var + S2var * Y;;
  if(QUADRATIC_VARIANCE)
    var += R2var * Y * Y;
  if(RES_VARIANCE){
    FLOAT resR = Ya[I] - Ya[I-K];
    FLOAT resL = Ya[I-K-n] - Ya[I-K-n-T];
    RFLOAT resvar = resL*resL + resR*resR;
    var += E2var*resvar;
  }
  if(DEBUG>=2) assert(var > 0.0);
  RFLOAT Bias = X * FrateBiasWT + m * FBiasWT;
  RFLOAT Gauss = RefPen - 0.5*mylog(var) - err*err/var;
  Cprtab *pr = &PRtab[n-1][I-K];
  RFLOAT Pen =  pr->LogPr - (FN_ADJUST ? pr->FnPenalty : (n-1)*FnPenalty) - (m-1)*FpPenalty;// WAS - X*Frate;
  RFLOAT Sm = PRtab[K][I].Sm;
  if(maptype){
    if(FRATE_FIX1)
      Gauss -= X*Frate;
    else
      Pen -= X*Frate;
    RFLOAT OutPen = OutlierPenalty;
    if(outlierBC)
      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
    if(OUTLIER_LTYPE==0)
      OutPen -= (X+Y) * OutlierLambdaInv;
    else
      OutPen -= fabs(X-Y) * OutlierLambdaInv;
    return /* OutlierBias + */ max(Bias + Pen + Gauss, Bias * biasWToutlierF + Pen*OUTLIER_TYPE1 + OutPen) + Sm;
  } else if(OUTLIER_DELTA(X-Y)){
    if(FRATE_FIX0)
      Gauss -= X*Frate;
    else
      Pen -= X*Frate;
    RFLOAT OutPen = OutlierPenalty;
    if(outlierBC)
      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
    OutPen -= fabs(X-Y) * OutlierLambdaInv;
    return /* OutlierBias + */ max(Bias + Pen + Gauss, Bias * biasWToutlierF + Pen*OUTLIER_TYPE + OutPen) + Sm;
  } else 
    return Bias + Pen + Gauss + Sm - X*Frate;
}

static inline RFLOAT SintY(int n,
                           int I,
			   int K,/**< at right end of interval! */
			   int T,/**< at left end of interval */
			   RFLOAT &deltaY,
			   RFLOAT &Ivar,
			   RFLOAT &Gauss,
			   RFLOAT &Sm,
                           CYPen ****YPenR,
                           int verb = 0)
{
  if(DEBUG>=2) assert(n>=1 && n <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2) assert(0 <= K && K <= KMAX);
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2){
    assert(0 <= T && T <= KMAX);
    assert(I-K-n-T >= 1);
  }

  CYPen *p = &YPenR[K][T][I][n];
  deltaY = p->deltaY;
  Ivar = p->Ivar;
  Gauss = p->Gauss;
  Sm = p->Sm;
  return p->Pen;
}

static inline RFLOAT SintX(RFLOAT X,
			   RFLOAT Y,
			   RFLOAT penY,// -(n-1)*FnPenalty + FpPenalty + PRtab[refid][n-1][I-K].LogPr
			   RFLOAT Gauss,
			   RFLOAT Sm,
			   RFLOAT Ivar,
                           RFLOAT XPen, // m*FpPenalty
			   RFLOAT XBias,
                           RFLOAT XPenBias,// XBias - XPen
			   RFLOAT OutPenBC)
{
  RFLOAT err = X-Y;
  Gauss -= err*err*Ivar;

  if(maptype){
    if(FRATE_FIX1)
      Gauss -= X*Frate;
    else
      penY -= X*Frate;
    RFLOAT OutPen = OutlierPenalty + OutPenBC;
    if(OUTLIER_LTYPE==0)
      OutPen -= (X+Y) * OutlierLambdaInv;
    else
      OutPen -= fabs(err) * OutlierLambdaInv;
    return /* OutlierBias + */ max(XPenBias + penY + Gauss, XBias * biasWToutlierF + (penY-XPen)*OUTLIER_TYPE1 + OutPen) + Sm;
  } else if(OUTLIER_DELTA(X-Y)){
    if(FRATE_FIX0)
      Gauss -= X * Frate;
    else
      penY -= X * Frate;
    RFLOAT OutPen = OutlierPenalty + OutPenBC;
    OutPen -= fabs(X-Y) * OutlierLambdaInv;
    return /* OutlierBias + */ max(XPenBias + penY + Gauss, XBias * biasWToutlierF + (penY-XPen)*OUTLIER_TYPE + OutPen) + Sm;
  } else
    return XPenBias + penY + Gauss + Sm - X * Frate;
}

/* same as SintX() with maptype==0 : Does NOT support outlierBC */
static inline RFLOAT SintX_maptype0(RFLOAT X,
				    RFLOAT Y,
				    RFLOAT penY,
				    RFLOAT GaussY,
				    RFLOAT Sm,
				    RFLOAT Ivar,
				    RFLOAT XPen,
				    RFLOAT XBias,
				    RFLOAT XPenBias,
				    int verb)
{
  RFLOAT err = X-Y;
  RFLOAT Gauss = GaussY - err*err*Ivar;
  RFLOAT OutPen = OutlierPenalty - fabs(X-Y) * OutlierLambdaInv;

  if(FRATE_FIX0)
    Gauss -= X*Frate;
  else
    penY -= X*Frate;

  RFLOAT Res =  OUTLIER_DELTA(err) ? 
    max(XPenBias + penY + Gauss, XBias * biasWToutlierF + (penY - XPen) * OUTLIER_TYPE + OutPen) : 
    XPenBias + penY + Gauss;
  if(verb){
    printf("SintX_maptype0(X=%0.6f,Y=%0.6f,penY=%0.6f,GaussY=%0.6f,Sm=%0.6f,var=%0.8e,ErrPen=%0.6f,XPen=%0.6f,XBias=%0.6f,XPenBias=%0.6f):Gauss=%0.6f,OutPen=%0.6f,Res=%0.6f,OUTLIER_DELTA(err)=%d\n",
	   X,Y,penY,GaussY,Sm,1.0/Ivar,-err*err*Ivar,XPen,XBias,XPenBias,Gauss,OutPen,Res,OUTLIER_DELTA(err) ? 1 : 0);
    fflush(stdout);
  }

  return Res + Sm;
}

/* same as SintX() with maptype==1 : Does NOT support outlierBC */
static inline RFLOAT SintX_maptype1(RFLOAT X,
				    RFLOAT Y,
				    RFLOAT penY,
				    RFLOAT Gauss,
				    RFLOAT Sm,
				    RFLOAT Ivar,
				    RFLOAT XPen,
				    RFLOAT XBias,
				    RFLOAT XPenBias)
{
  RFLOAT err = X-Y;
  Gauss -= err*err*Ivar;

  if(FRATE_FIX1)
    Gauss -= X*Frate;
  else
    penY -= X*Frate;

  RFLOAT OutPen = OutlierPenalty;
  if(OUTLIER_LTYPE==0)
    OutPen -= (X+Y) * OutlierLambdaInv;
  else
    OutPen -= fabs(err) * OutlierLambdaInv;

  RFLOAT Res = /* OutlierBias + */ max(XPenBias + penY + Gauss, XBias * biasWToutlierF + (penY - XPen)*OUTLIER_TYPE1 + OutPen);
  return Res + Sm;
}

#define VS_deltaX 0
#define VS_Pen 4
#define VS_Bias 8
#define VS_PenBias 12

static inline RFLOAT SintX_BLOCK4_maptype0(RFLOAT score, CYPen *p, RFLOAT *var_scores, RFLOAT * AscoreGT, int J, int Hmin, int Hmax, int G, int *JMIN, int *JMAX)
{
  RFLOAT deltaY = p->deltaY;
  RFLOAT Ivar = p->Ivar;
  RFLOAT penY = p->Pen;
  RFLOAT Gauss = p->Gauss;
  RFLOAT Sm = p->Sm;

  #define BLOCK(m) { \
    int H = J - (m);   \
    if(H >= Hmin && (!DIAGONAL || H <= Hmax)) {				\
      if(DEBUG>=2) assert(JMIN[G] <= H); \
      if(DEBUG>=2) assert(H <= JMAX[G]); \
      RFLOAT newscore = AscoreGT[H] + SintX_maptype0(var_scores[4-m+VS_deltaX],deltaY,penY,Gauss,Sm,Ivar,var_scores[4-m+VS_Pen],var_scores[4-m+VS_Bias], var_scores[4-m+VS_PenBias],0); \
      score = max(score,newscore); \
     } \
  }
  BLOCK(4)
  BLOCK(3)
  BLOCK(2)
  BLOCK(1)
  #undef BLOCK

  return score;
}


#if (USE_MIC==1 && USE_RFLOAT==1)

/* MIC vectorized version of SintX_maptype0() */
static inline __m512 SintX_mm512_maptype0(__m512 X,
					  __m512 Y,
					  __m512 penY,
					  __m512 Gauss,
					  __m512 Sm,
					  __m512 Ivar,
					  __m512 XPenJM_Bias,
					  __m512 XPenJM_Pen,
					  __m512 XPenJM_PenBias,
					  int verb, float *deltaX)
{
  if(DEBUG>=2) assert(OUTLIER_MAX>=1);

  __m512 err = _mm512_sub_ps(X, Y);
  Gauss = _mm512_sub_ps(Gauss, _mm512_mul_ps(err, _mm512_mul_ps(err, Ivar)));

  if(FRATE_FIX0)
    Gauss = _mm512_sub_ps(Gauss, _mm512_mul_ps(X, v512_Frate));
  else
    penY = _mm512_sub_ps(penY, _mm512_mul_ps(X, v512_Frate));

  __m512 OutPenLambda = _mm512_mul_ps(err, v512_OutlierLambdaInv);
  __m512 OutPen = _mm512_mask_sub_ps(_mm512_add_ps(v512_OutlierPenalty, OutPenLambda),
				     _mm512_cmp_ps_mask(X, Y, _MM_CMPINT_GT),
				     v512_OutlierPenalty, OutPenLambda);// (X > Y) ? OutlierPenalty - (X-Y) * OutlierLambdaInv : OutlierPenalty + (X-Y) * OutlierLambdaInv;

  __m512 XYpen = _mm512_add_ps(penY, Gauss);// penY + Gauss
  __m512 Res2 = _mm512_add_ps(XPenJM_PenBias, XYpen);// result without outlier : XPenBias + penY + Gauss
  __m512 OutPenBias = _mm512_add_ps(OutPen, _mm512_mul_ps(XPenJM_Bias, v512_biasWToutlierF));// OutPen + XBias * biasWToutlierF
  __m512 Res1;
  __m512 Pen = _mm512_sub_ps(penY, XPenJM_Pen);// penY - XPen
  Res1 = _mm512_gmax_ps(Res2, _mm512_add_ps(OutPenBias, _mm512_mul_ps(Pen, v512_OUTLIER_TYPE)));

  if(verb){
    __mmask16 mask1 = _mm512_cmp_ps_mask(err, v512_outlierNegMax, _MM_CMPINT_GT);
    __mmask16 mask2 =  _mm512_cmp_ps_mask(err, v512_outlierMax, _MM_CMPINT_LT);
    __mmask16 mask = mask1 & mask2; // OUTLIER_DELTA(err)

    float er[16],gauss[16],outpen[16],xypen[16],res2[16],penbias[16],outpenbias[16],res1[16],pen[16],outlierMax[16],outlierNegMax[16];
    _mm512_storeu_ps(er, err);
    _mm512_storeu_ps(gauss, Gauss);
    _mm512_storeu_ps(outpen,OutPen);
    _mm512_storeu_ps(xypen, XYpen);
    _mm512_storeu_ps(res2,Res2);
    _mm512_storeu_ps(outpenbias,OutPenBias);
    _mm512_storeu_ps(res1,Res1);
    _mm512_storeu_ps(pen, Pen);
    _mm512_storeu_ps(outlierMax, v512_outlierMax);
    _mm512_storeu_ps(outlierNegMax, v512_outlierNegMax);

    #pragma omp critical
    {
      for(int j = 0; j < 16; j++)
	if(deltaX[j] > 0.0)
	  if(OUTLIER_TYPE==0)
	    printf("j=%d:X=%0.4f,err=%0.6f,Gauss=%0.6f,OutPen=%0.6f,XYpen=%0.6f,Res2=%0.6f,OutPenBias=%0.6f,Res1=%0.6f,msk=%d(msk1=%d,msk2=%d,outMax=%0.2f,%0.2f)\n",
		   j, deltaX[j], er[j],gauss[j],outpen[j],xypen[j],res2[j],outpenbias[j],res1[j], (mask >> j) & 1, (mask1 >> j) & 1, (mask2 >> j) & 1, outlierMax[j], outlierNegMax[j]);
	  else
	    printf("j=%d:X=%0.4f,err=%0.6f,Gauss=%0.6f,OutPen=%0.6f,XYpen=%0.6f,Res2=%0.6f,OutPenBias=%0.6f,Res1=%0.6f,Pen=%0.6f,msk=%d(msk1=%d,msk2=%d,outMax=%0.2f,%0.2f)\n",
		   j, deltaX[j], er[j],gauss[j],outpen[j],xypen[j],res2[j],outpenbias[j],res1[j],pen[j], (mask >> j) & 1, (mask1 >> j) & 1, (mask2 >> j) & 1, outlierMax[j], outlierNegMax[j]);
      fflush(stdout);
    }
  }

  __m512 Res = _mm512_mask_mov_ps(Res2, 
				  _mm512_cmp_ps_mask(err, v512_outlierNegMax, _MM_CMPINT_GT) & _mm512_cmp_ps_mask(err, v512_outlierMax, _MM_CMPINT_LT), // OUTLIER_DELTA(err)
				  Res1);
  return _mm512_add_ps(Res, Sm);
}

static inline __m512 SintX_mm512_maptype1(__m512 X,
					  __m512 Y,
					  __m512 penY,
					  __m512 Gauss,
					  __m512 Sm,
					  __m512 Ivar,
					  __m512 XPenJM_Bias,
					  __m512 XPenJM_Pen,
					  __m512 XPenJM_PenBias
 					)
{
  __m512 err = _mm512_sub_ps(X, Y);
  Gauss = _mm512_sub_ps(Gauss, _mm512_mul_ps(err, _mm512_mul_ps(err, Ivar)));

  if(FRATE_FIX1)
    Gauss = _mm512_sub_ps(Gauss, _mm512_mul_ps(X, v512_Frate));
  else
    penY = _mm512_sub_ps(penY, _mm512_mul_ps(X, v512_Frate));

  __m512 OutPen;
  if(OUTLIER_LTYPE==0)
    OutPen = _mm512_sub_ps(v512_OutlierPenalty, _mm512_mul_ps(_mm512_add_ps(X, Y), v512_OutlierLambdaInv));
  else {
    __m512 OutPenLambda = _mm512_mul_ps(err, v512_OutlierLambdaInv);
    OutPen = _mm512_mask_sub_ps(_mm512_add_ps(v512_OutlierPenalty, OutPenLambda),
				_mm512_cmp_ps_mask(X, Y, _MM_CMPINT_GT),
				v512_OutlierPenalty, OutPenLambda);
  }
 
  __m512 XYPen  = _mm512_add_ps(penY, Gauss);// penY + Gauss
  __m512 OutPenBias = _mm512_add_ps(OutPen, _mm512_mul_ps(XPenJM_Bias, v512_biasWToutlierF));// XPenBias * biasWToutlierF + OutPen
  __m512 Res2 = _mm512_add_ps(XPenJM_PenBias, XYPen);/* result without outlier : XPenBias + penY + Guass */
  __m512 Res1;
  __m512 Pen = _mm512_sub_ps(penY, XPenJM_Pen);
  Res1 = _mm512_max_ps(Res2, _mm512_add_ps(OutPenBias, _mm512_mul_ps(Pen, v512_OUTLIER_TYPE1)));
    
  return _mm512_add_ps(Res1, Sm);
}

typedef struct {
  __m512 scoreGT;
  
  __m512 deltaX;
  __m512 deltaY;
  __m512 penY;
  __m512 Gauss;
  __m512 Sm;
  __m512 Ivar;
  __m512 Bias;
  __m512 Pen;
  __m512 PenBias;
  
  int usage;
  
  __mmask16 mask;
} SintX_scratch_mm512;
	
static inline __m512 SintX_scratch_collapse_maptype0(SintX_scratch_mm512 &s, __m512 score)
{
  __m512 newscore = _mm512_mask_gmax_ps(score, s.mask, score,
					_mm512_add_ps(
						      s.scoreGT,
						      SintX_mm512_maptype0(
									   s.deltaX,
									   s.deltaY,
									   s.penY,
									   s.Gauss,
									   s.Sm,
									   s.Ivar,
									   s.Bias,
									   s.Pen,
									   s.PenBias,
									   0,NULL
									   )
						      ));
  s.mask = _mm512_int2mask(0);
  s.usage = 0;
  return newscore;
}

#if 0
static inline void SintX_load_scratch(SintX_scratch_mm512 &s, const int USAGE, __mmask16 mask, __m512 vs, __m512 vAscoreGT, RFLOAT deltaY, RFLOAT penY, RFLOAT Ivar)
{
  s.scoreGT = _mm512_mask_alignr_ps(s.scoreGT, mask, vAscoreGT, vAscoreGT, (16-USAGE) & 0xf); 
}
#endif

/* NOTE : The number of score values to copy is J - m_min .. J - m_max which is a subset of the range J-1 .. J-4 (ie m_min >= 1 && m_max <= 4) */
static inline __m512 SintX_scratch_feed_maptype0(SintX_scratch_mm512 &s, __m512 score, CYPen *p, RFLOAT * AscoreGT, int J, int m_min, int m_max, int verb, int T)
{
  RFLOAT deltaY = p->deltaY;
  RFLOAT Ivar = p->Ivar;
  RFLOAT penY = p->Pen;
  RFLOAT Gauss = p->Gauss;
  RFLOAT Sm = p->Sm;

  __m512 newscore;
  if(s.usage > 12) {
    if(verb){
      __m512 new_score = score;
      SintX_scratch_mm512 new_scratch = s;
      __mmask16 mask = s.mask;
      new_score = SintX_scratch_collapse_maptype0(new_scratch, new_score);
      float origscore[16], newscore[16], scoreGT[16], deltaX[16], deltaY[16], penY[16], Gauss[16], Sm[16], Ivar[16], Bias[16], Pen[16], PenBias[16];
      _mm512_maskstoreu_ps(origscore, mask, score);
      _mm512_maskstoreu_ps(newscore, mask, new_score);
      _mm512_maskstoreu_ps(scoreGT, mask, s.scoreGT);
      _mm512_maskstoreu_ps(deltaX, mask, s.deltaX);
      _mm512_maskstoreu_ps(deltaY, mask, s.deltaY);
      _mm512_maskstoreu_ps(penY, mask, s.penY);
      _mm512_maskstoreu_ps(Sm, mask, s.Gauss);
      _mm512_maskstoreu_ps(Sm, mask, s.Sm);
      _mm512_maskstoreu_ps(Ivar, mask, s.Ivar);
      _mm512_maskstoreu_ps(Bias, mask, s.Bias);
      _mm512_maskstoreu_ps(Pen, mask, s.Pen);
      _mm512_maskstoreu_ps(PenBias, mask, s.PenBias);

      #pragma omp critical
      {
	printf("SintX_scratch_collapse_maptype0() called in SintX_scratch_feed_maptype0:s.usage=%d,J=%d,T=%d\n",s.usage,J,T);
	for(int j = 0; j < 16; j++)
	  printf("ht=%d:deltaX=%0.4f,deltaY=%0.4f,penY=%0.6f,Gauss=%0.6f,Sm=%0.6f,Ivar=%0.8f,Bias=%0.6f,Pen=%0.6f,PenBias=%0.6f:score=%0.6f -> %0.6f\n",
		 j,deltaX[j],deltaY[j],penY[j],Gauss[j],Sm[j],Ivar[j],Bias[j],Pen[j],PenBias[j],origscore[j],newscore[j]);
	fflush(stdout);
      }
    }
    newscore = SintX_scratch_collapse_maptype0(s, score);
  } else 
    newscore = score;

  if(verb || (DEBUG>=2 && !(1 <= m_min && m_min <= m_max && m_max <= 4))){
    #pragma omp critical
    {
      printf("SintX_scratch_feed_maptype0:usage=%d->%d,J=%d,T=%d,H=%d..%d:deltaY=%0.4f,penY=%0.6f,Gauss=%0.6f,Sm=%0.6f,Ivar=%0.8f:\n",
	     s.usage, s.usage+4, J, T, J - m_max, J - m_min, deltaY, penY, Gauss, Sm, Ivar);
      for(int H = J - m_max; H <= J - m_min; H++)
	printf("  AscoreGT[H=%d] = %0.6f\n",H, AscoreGT[H]);
      fflush(stdout);

      assert(1 <= m_min && m_min <= m_max && m_max <= 4);
    }
  }
  __mmask16 mask = ((0xf >> (m_min-1)) & (0xf << (4-m_max))) << s.usage;
  s.mask |= mask;

  s.deltaY = _mm512_mask_mov_ps(s.deltaY, mask, _mm512_set1_ps(deltaY));
  s.penY = _mm512_mask_mov_ps(s.penY, mask, _mm512_set1_ps(penY));
  s.Gauss = _mm512_mask_mov_ps(s.Gauss, mask, _mm512_set1_ps(Gauss));
  s.Sm = _mm512_mask_mov_ps(s.Sm, mask, _mm512_set1_ps(Sm));
  s.Ivar = _mm512_mask_mov_ps(s.Ivar, mask, _mm512_set1_ps(Ivar));

  s.scoreGT = _mm512_mask_mov_ps(s.scoreGT, mask, _mm512_loadu_ps(&AscoreGT[J-4-s.usage]));

  s.usage += 4;
  return newscore;
}

#endif // USE_MIC==1 && USE_RFLOAT==1 


#if (USE_SSE==1 && USE_AVX==1 && USE_RFLOAT==1)

/* same as SintX() with maptype==0 : Does NOT support outlierBC */
static inline __m256 SintX_mm256_maptype0(__m256 X,
					  __m256 Y,
					  __m256 penY,
					  __m256 Gauss,
					  __m256 Sm,
					  __m256 Ivar,
					  __m256 XPenJM_Bias,
					  __m256 XPenJM_Pen,
					  __m256 XPenJM_PenBias
					  )
{
  if(DEBUG>=2) assert(OUTLIER_MAX >= 1);

  __m256 err = _mm256_sub_ps(X, Y);
  Gauss = _mm256_sub_ps(Gauss, _mm256_mul_ps(err, _mm256_mul_ps(err, Ivar)));

  if(FRATE_FIX0)
    Gauss = _mm256_sub_ps(Gauss, _mm256_mul_ps(X, v256_Frate));
  else
    penY = _mm256_sub_ps(penY, _mm256_mul_ps(X, v256_Frate));

  __m256 OutPenLambda = _mm256_mul_ps(err, v256_OutlierLambdaInv);
  __m256 OutPen = _mm256_blendv_ps(_mm256_add_ps(v256_OutlierPenalty, OutPenLambda),
				   _mm256_sub_ps(v256_OutlierPenalty, OutPenLambda),
				   _mm256_cmp_ps(X, Y, _CMP_GT_OQ)); // (X > Y) ? OutlierPenalty - (X-Y) * OutlierLambdaInv : OutlierPenalty + (X-Y) * OutlierLambdaInv;

  __m256 XYPen = _mm256_add_ps(penY, Gauss);
  __m256 Res2 = _mm256_add_ps(XPenJM_PenBias, XYPen);
  __m256 OutPenBias = _mm256_add_ps(_mm256_mul_ps(XPenJM_Bias, v256_biasWToutlierF), OutPen);
  __m256 Res1 = _mm256_max_ps(Res2, _mm256_add_ps(OutPenBias, _mm256_mul_ps(_mm256_sub_ps(penY, XPenJM_Pen), v256_OUTLIER_TYPE)));

  __m256 Res = _mm256_blendv_ps(Res2, Res1,
				_mm256_and_ps(_mm256_cmp_ps(err, v256_outlierNegMax, _CMP_GT_OQ), _mm256_cmp_ps(err, v256_outlierMax, _CMP_LT_OQ))  /* OUTLIER_DELTA(err) */
				);

  return _mm256_add_ps(Res, Sm);
}

/* same as SintX() with maptype==1 : Does NOT support outlierBC */
static inline __m256 SintX_mm256_maptype1(__m256 X,
					  __m256 Y,
					  __m256 penY,
					  __m256 Gauss,
					  __m256 Sm,
					  __m256 Ivar,
					  __m256 XPenJM_Bias,
					  __m256 XPenJM_Pen,
					  __m256 XPenJM_PenBias
 					)
{
  __m256 err = X - Y;
  Gauss -= err * err * Ivar;
 
  if(FRATE_FIX1)
    Gauss -= X * v256_Frate;
  else
    penY -= X * v256_Frate;

  __m256 OutPen;
  if(OUTLIER_LTYPE==0)
    OutPen = v256_OutlierPenalty - (X+Y) * v256_OutlierLambdaInv;
  else {
    __m256 OutPenLambda = err * v256_OutlierLambdaInv;
    OutPen = _mm256_blendv_ps(v256_OutlierPenalty + OutPenLambda,
			      v256_OutlierPenalty - OutPenLambda,
			      _mm256_cmp_ps(X, Y, _CMP_GT_OQ)); // (X > Y) ? OutlierPenalty - (X-Y) * OutlierLambdaInv : OutlierPenalty + (X-Y) * OutlierLambdaInv;
  }

  __m256 PenBias = XPenJM_PenBias + penY + Gauss;
  __m256 OutPenBias = XPenJM_Bias * v256_biasWToutlierF + OutPen;// XBias * biasWToutlierF + OutPen
  __m256 Res1 = _mm256_max_ps(PenBias, OutPenBias + (penY - XPenJM_Pen) * v256_OUTLIER_TYPE1);// max(XPenBias + penY + Gauss, OutPen+XBias*biasWToutlierF + (penY-XPen)*OUTLIER_TYPE1)

  return Res1 + Sm;
}

/* AVX256 version of SintX() with maptype as template parameter AND supports outlierBC */
template<int maptype>
static inline __m256 SintX_mm256(__m256 X,
				 __m256 Y,
				 __m256 penY,
				 __m256 Gauss,
				 __m256 Sm,
				 __m256 Ivar,
				 __m256 XPenJM_Bias,
				 __m256 XPenJM_Pen,
				 __m256 XPenJM_PenBias,
				 RFLOAT outlierBC
				 )
{
  if(DEBUG>=2) assert(OUTLIER_MAX >= 1);

  __m256 err = X - Y;
  Gauss -= err * err * Ivar;

  if(maptype ? FRATE_FIX1 : FRATE_FIX0)
    Gauss -= X * v256_Frate;
  else
    penY -= X * v256_Frate;

  __m256 v_OutlierPenalty = _mm256_set1_ps(OutlierPenalty + outlierBC);

  __m256 OutPen;
  if(OUTLIER_LTYPE==0 && maptype==1)
    OutPen = v_OutlierPenalty - (X+Y) * v256_OutlierLambdaInv;
  else {
    __m256 OutPenLambda = err * v256_OutlierLambdaInv;
    OutPen = _mm256_blendv_ps(v_OutlierPenalty + OutPenLambda,   v_OutlierPenalty - OutPenLambda,
			      _mm256_cmp_ps(X, Y, _CMP_GT_OQ)); // (X > Y) ? OutlierPenalty - (X-Y) * OutlierLambdaInv : OutlierPenalty + (X-Y) * OutlierLambdaInv;
  }

  __m256 PenBias = XPenJM_PenBias + penY + Gauss;
  __m256 OutPenBias = XPenJM_Bias * v256_biasWToutlierF + OutPen;
  __m256 Res1 = _mm256_max_ps(PenBias, OutPenBias + (penY - XPenJM_Pen) * v256_OUTLIER_TYPE);
  if(maptype == 0){
    __m256 Res = _mm256_blendv_ps(PenBias, Res1,
				  _mm256_and_ps(_mm256_cmp_ps(err, v256_outlierNegMax, _CMP_GT_OQ), _mm256_cmp_ps(err, v256_outlierMax, _CMP_LT_OQ)));  /* !OUTLIER_DELTA(err) ? PenBias : Res1 */
    return Res + Sm;
  } else
    return Res1 + Sm;
}

#endif // USE_SSE==1 && USE_AVX==1 && USE_RFLOAT==1

// This version of SintDetail does not use Tables and also computes the iscore that includes the Outlier terms
void SintDetail(double X,
		  double Y,
		  int m,
		  int n,
		  int J,
		  int I,
		  int K,/**< right end of interval! */
		  int T,/**< left end of interval : only used with RES_VARIANCE */
		  FLOAT *Ya,/**< reference array Ya[0..N+1] */
		  double &Bias,
		  double &Pen,
		  double &Gauss,
		  double &PenSm,
		  double &OutPen,
		  double &iscore,
		  int verb)
{
  if(DEBUG>=2) assert(0 <= K && K <= KMAX);
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= KMAX);
    assert(I-K-n-T >= 1);
  }

  double err = X-Y;
  double var = F2var + S2var*Y;
  if(QUADRATIC_VARIANCE)
    var += R2var * Y * Y;
  if(RES_VARIANCE){
    double resR = Ya[I] - Ya[I-K];
    double resL = Ya[I-K-n] - Ya[I-K-n-T];
    double resvar = resL*resL + resR*resR;
    var += E2var * resvar;
  }
  if(DEBUG>=2) assert(var > 0.0);
  
  Bias = X * FrateBiasWT + m * FBiasWT;
  Gauss = RefPen - 0.5*mylog(var) - err*err/var;
  Pen = log(Pr(Ya[I-K]-Ya[I-K-n])) - (m-1)*FpPenalty;// WAS - X*Frate;
  double FnPen = 0.0;
  if(FN_ADJUST){
    double sum = 0.0;
    for(int i = 1; i <= n; i++)
      sum += Pr(Ya[I-i+1]-Ya[I-i]);
    FnPen = max(0.0,sum-1.0)*FnPenalty;
  } else
    FnPen = (n-1)*FnPenalty;
  Pen -= FnPen;

  PenSm = Sm(0,I,K,Ya);

  if(maptype){
    if(FRATE_FIX1)
      Gauss -= X*Frate;
    else 
      Pen -= X*Frate;
    OutPen = OutlierPenalty;
    if(outlierBC)
      OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
    if(OUTLIER_LTYPE==0)
      OutPen -= (X+Y) * OutlierLambdaInv;
    else
      OutPen -= fabs(err) * OutlierLambdaInv;
    iscore = OutlierBias + max(Bias + Pen + Gauss, Bias * biasWToutlier + Pen * OUTLIER_TYPE1 + OutPen) + PenSm;
  } else {// maptype == 0
    if(FRATE_FIX0)
      Gauss -= X*Frate;
    else
      Pen -= X*Frate;
    if(OUTLIER_DELTA(err)){
      OutPen = OutlierPenalty;
      if(outlierBC)
	OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];
      OutPen -= fabs(err) * OutlierLambdaInv;
      iscore = OutlierBias + max(Bias + Pen + Gauss, Bias * biasWToutlier + Pen * OUTLIER_TYPE + OutPen) + PenSm;
    } else {
      OutPen = 0.0;
      iscore = Bias + Pen + Gauss + PenSm;
    }
  }

  if(verb){
    printf("\n\t\tI=%d,K=%d,J=%d,n=%d,m=%d,X=%0.4f,Y=%0.4f:Gterm=%0.6f,ErrPen=%0.6f,X*Frate=%0.6f,FnPen=%0.6f,FpPen=%0.6f,X*LBias=%0.6f,m*FBias=%0.6f,PenPr=%0.6f,PenSm=%0.6f,err=%0.6f,iscore= %0.6f\n\t\t", 
	   I,K,J,n,m, X, Y, RefPen - 0.5*mylog(var), -err*err/var, -X*Frate, -FnPen, -(m-1)*FpPenalty, X*LBias, m*FBias, log(Pr(Ya[I-K]-Ya[I-K-n])), Sm(0,I,K,Ya), err, iscore);
    fflush(stdout);
  }
}

/** SintDetail with resSD>0.0 */
static void SintDetail(double X,
		       double Y,
		       int m,
		       int n,
		       int J,
		       int I,
		       int K,/**< right end of interval! */
		       int T,/**< left end of interval : only used with RES_VARIANCE */
		       Cprtab **PRtab,
		       FLOAT *Ya,/**< reference array Ya[0..N+1] */
		       double &Bias,
		       double &Pen,
		       double &Gauss,
		       double &PenSm,
		       int verb)
{
  if(DEBUG>=2) assert(0 <= K && K <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2) assert(fabs(PRtab[K][I].Sm - (RFLOAT)Sm(0,I,K,Ya)) <= 1.0e-10);
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= max(max(2,KMAX),DELTA_Y));
    assert(I-K-n-T >= 1);
  }
  PenSm = PRtab[K][I].Sm;

  double err = X-Y;
  double var = F2var + S2var*Y;
  if(QUADRATIC_VARIANCE)
    var += R2var * Y * Y;
  if(RES_VARIANCE){
    double resR = Ya[I] - Ya[I-K];
    double resL = Ya[I-K-n] - Ya[I-K-n-T];
    double resvar = resL*resL + resR*resR;
    var += E2var * resvar;
  }
  if(DEBUG>=2) assert(var > 0.0);

  Bias = X * FrateBiasWT + m * FBiasWT;
  Gauss = RefPen - 0.5*mylog(var) - err*err/var;

  if(n > DELTA_Y){
    if(DEBUG>=2) assert(outlierExtend > 0);
    if(DEBUG>=2) assert(n >= 1);

    Pen = /* log(Pr(Ya[I-K]-Ya[I-K-n])) */ - (m-1)*FpPenalty; // WAS - X*Frate;
    if(FN_ADJUST){
      double sum = 0.0;
      for(int i = 1; i <= n; i++)
	sum += Pr(Ya[I-i+1]-Ya[I-i]);
      Pen -= max(0.0,sum-1.0)*FnPenalty;
    } else
      Pen -= (n-1)*FnPenalty;

  } else {

    if(DEBUG>=2) assert(n >= 1 && n <= max(max(2,KMAX),DELTA_Y));
    if(DEBUG>=2) assert(fabs(PRtab[n-1][I-K].LogPr - log(Pr(Ya[I-K]-Ya[I-K-n]))) <= 1.0e-10);

    Cprtab *pr = &PRtab[n-1][I-K];
    Pen = pr->LogPr - (FN_ADJUST ? pr->FnPenalty : (n-1)*FnPenalty) - (m-1)*FpPenalty;// WAS - X*Frate;

    if(verb){
      printf("\n\t\tI=%d,K=%d,J=%d,n=%d,m=%d:x=%0.4f,y=%0.4f,Gterm=%0.6f,ErrPen=%0.6f,X*Frate=%0.6f,FnPen=%0.6f,FpPen=%0.6f,X*LBias=%0.6f,m*FBias=%0.6f,PenPr=%0.6f,PenSm=%0.6f,%0.6f (err=%0.6f,sf=%0.4f,sd=%0.6f,sr=%0.6f,se=%0.6f,var2=%0.8f)\n\t\t", 
	     I,K,J,n,m, X,Y,RefPen - 0.5*mylog(var), -err*err/var, -X*Frate, -(FN_ADJUST ? pr->FnPenalty : (n-1)*FnPenalty), -(m-1)*FpPenalty, X*LBias, m*FBias, pr->LogPr, PRtab[K][I].Sm, Sm(0,I,K,Ya), fabs(err), SF[0],SD[0],SR[0],SE[0],var);
      printf("\t\t Ya[I-K]= %0.6f,Ya[I-K-n]= %0.6f  :delta= %0.6e,log(Pr(delta))= %0.6e), PRtab= %p\n",Ya[I-K],Ya[I-K-n],Ya[I-K]-Ya[I-K-n],log(Pr(Ya[I-K]-Ya[I-K-n])), PRtab);
      fflush(stdout);
    }
  }
  if(maptype){
    if(FRATE_FIX1)
      Gauss -= X*Frate;
    else
      Pen -= X*Frate;
  } else {
    if(FRATE_FIX0)
      Gauss -= X*Frate;
    else
      Pen -= X*Frate;
  }
}

#if USE_RFLOAT==1
static void SintDetail(float X,
		       float Y,
		       int m,
		       int n,
		       int J,
		       int I,
		       int K,/**< right end of interval! */
		       int T,/**< left end of interval : only used with RES_VARIANCE */
		       Cprtab **PRtab,
		       FLOAT *Ya,/**< reference array Ya[0..N+1] */
		       float &Bias,
		       float &Pen,
		       float &Gauss,
		       float &PenSm,
		       int verb)
{
  if(DEBUG>=2) assert(0 <= K && K <= max(max(2,KMAX),DELTA_Y));
  if(DEBUG>=2) assert(fabs(PRtab[K][I].Sm - (RFLOAT)Sm(J,I,K,Ya)) <= 1.0e-10);
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= max(max(2,KMAX),DELTA_Y));
    assert(I-K-n-T >= 1);
  }
  PenSm = PRtab[K][I].Sm;
  
  float err = X-Y;
  FLOAT var, y = Y;
  if(RES_VARIANCE){
    FLOAT resR = Ya[I] - Ya[I-K];
    FLOAT resL = Ya[I-K-n] - Ya[I-K-n-T];
    FLOAT resvar = resL*resL + resR*resR;
    var = QUADRATIC_VARIANCE ? (F2var + S2var*y + R2var*y*y + E2var*resvar) : (F2var + S2var*y + E2var*resvar);
  } else
    var = QUADRATIC_VARIANCE ? (F2var + S2var*y + R2var*y*y) : (F2var + S2var*y);

  float Ivar = 1.0/var;

  if(DEBUG>=2) assert(var > 0.0);

  Bias = X * FrateBiasWT + m * FBiasWT;
  Gauss = RefPen - 0.5f * mylog((float)var) - err*err * Ivar;

  if(n > DELTA_Y){
    if(DEBUG>=2) assert(outlierExtend > 0 || RepeatRec > 0);
    if(DEBUG>=2) assert(n >= 1);

    if(FN_ADJUST){
      double sum = 0.0;
      for(int i = 1; i <= n; i++)
	sum += Pr(Ya[I-i+1]-Ya[I-i]);
      Pen = /* log(Pr(Ya[I-K]-Ya[I-K-n])) */ - max(0.0,sum-1.0)*FnPenalty - (m-1)*FpPenalty;// WAS - X*Frate;
    } else
      Pen = -(n-1)*FnPenalty - (m-1)*FpPenalty;// WAS -X*Frate

    if(verb){
      printf("\n\t I=%d,K=%d,J=%d,n=%d,m=%d,D=%d:X=%0.4f,Y=%0.4f:Gterm=%0.6f(var=%0.8e,Ivar=%0.8e,Xtheta=%0.6f),ErrPen=%0.6f,X*Frate=%0.6f,FnPen=%0.6f,FpPen=%0.6f,X*LBias=%0.6f,m*FBias=%0.6f,PenSm=%0.6f,%0.6f (err=%0.6f,SF=%0.4f,SD=%0.6f,SR=%0.6f,var=%0.8f)\n\t\t", 
	     I,K,J,n,m,T,X,Y,RefPen - 0.5*mylog(var), var*0.5, Ivar, Xtheta, -err*err*Ivar, -X*Frate, Pen + (m-1)*FpPenalty, -(m-1)*FpPenalty, X*LBias, m*FBias, PRtab[K][I].Sm, Sm(0,I,K,Ya), fabs(err), SF[0],SD[0],SR[0],var);
      FLOAT resR = Ya[I] - Ya[I-K];
      FLOAT resL = Ya[I-K-n] - Ya[I-K-n-T];
      FLOAT resvar = resL*resL + resR*resR;
      printf("\t  F2var=%0.8e,S2var=%0.8e,R2var=%0.8e,E2var=%0.8e,Y=%0.8f,resL=%0.6f,resR=%0.6f,resvar=%0.8e\n",F2var,S2var,R2var,E2var,Y,resL,resR,resvar);
      fflush(stdout);
    }
  } else {
    if(DEBUG>=2) assert(n >= 1 && n <= max(max(2,KMAX),DELTA_Y));
    if(DEBUG>=2 && !(fabs(PRtab[n-1][I-K].LogPr - log(Pr(Ya[I-K]-Ya[I-K-n]))) <= 1.0e-10)){
      #pragma omp critical
      {
	printf("SintDetail:X=%0.4f,Y=%0.4f,m=%d,n=%d,I=%d,K=%d,J=%d,T=%d:PRtab[n-1][I-K].LogPr=%0.10e,Ya[I-K]=%0.4f,Ya[I-K-n]=%0.4f,log(Pr(Ya[I-K]-Ya[I-K-n]))=%0.10e\n",
	       X,Y,m,n,I,K,J,T,PRtab[n-1][I-K].LogPr,Ya[I-K],Ya[I-K-n],log(Pr(Ya[I-K]-Ya[I-K-n])));
	fflush(stdout);
	assert(fabs(PRtab[n-1][I-K].LogPr - log(Pr(Ya[I-K]-Ya[I-K-n]))) <= 1.0e-10);
      }
    }

    Cprtab *pr = &PRtab[n-1][I-K];
    Pen = pr->LogPr - (FN_ADJUST ? pr->FnPenalty : (n-1)*FnPenalty) - (m-1)*FpPenalty;// WAS - X*Frate;
  
    if(verb){
      printf("\n\t I=%d,K=%d,J=%d,n=%d,m=%d,D=%d:X=%0.4f,Y=%0.4f:Gterm=%0.6f(var=%0.8e,Ivar=%0.8e,Xtheta=%0.6f),ErrPen=%0.6f,X*Frate=%0.6f,FnPen=%0.6f,FpPen=%0.6f,X*LBias=%0.6f,m*FBias=%0.6f,PenPr=%0.6f,PenSm=%0.6f,%0.6f (err=%0.6f,SF=%0.4f,SD=%0.6f,SR=%0.6f,var2=%0.8f)\n\t\t", 
	     I,K,J,n,m,T,X,Y,RefPen - 0.5*mylog(var), var*0.5, Ivar, Xtheta, -err*err*Ivar, -X*Frate, -(FN_ADJUST ? pr->FnPenalty : (n-1)*FnPenalty), -(m-1)*FpPenalty, X*LBias, m*FBias, pr->LogPr, PRtab[K][I].Sm, Sm(0,I,K,Ya), fabs(err), SF[0],SD[0],SR[0],var);
      FLOAT resR = Ya[I] - Ya[I-K];
      FLOAT resL = Ya[I-K-n] - Ya[I-K-n-T];
      FLOAT resvar = resL*resL + resR*resR;
      printf("\t  F2var=%0.8e,S2var=%0.8e,R2var=%0.8e,E2var=%0.8e,Y=%0.8f,resL=%0.6f,resR=%0.6f,resvar=%0.8e\n",F2var,S2var,R2var,E2var,Y,resL,resR,resvar);
      fflush(stdout);
    }
  }
  if(maptype){
    if(FRATE_FIX1)
      Gauss -= X*Frate;
    else
      Pen -= X*Frate;
  } else {
    if(FRATE_FIX0)
      Gauss -= X*Frate;
    else
      Pen -= X*Frate;
  }
}
#endif

/** short for Send(X,m,n) + Sint(X,Y,1,1) - Sint((X+Y)/2,(X+Y)/2,1,1)  with resSD > 0.0 */
/* NOTE : ignores RES_VARIANCE, since information about mis-resolved interval is hard to pass */
static inline RFLOAT Sbnd(RFLOAT X,
			  RFLOAT Y,
			  int m,
			  int n,
			  int I,/* right end of labels 1..n, excluding end */
			  int refid)
{
  if(DEBUG>=2)assert(n>=1);
  if(DEBUG>=2)assert(m>=1);
  RFLOAT err = X-Y;
  RFLOAT var = QUADRATIC_VARIANCE ? (F2var + S2var*Y + R2var*Y*Y) : (F2var + S2var*Y);
  if(DEBUG>=2) assert(var > 0.0);

  RFLOAT Bias = X * FrateBiasWT + (m-1)*FBiasWT +  EBiasWT;
  return Bias - (FN_ADJUST ? max(0.0, PRsum[refid][I] - PRsum[refid][I-n+1]) : (n-1)) * FnPenalty - (m-1)*FpPenalty - X*Frate - err*err/var;
}

static inline RFLOAT Send(RFLOAT X,
			  int m,
			  int n,
			  int I,/* rightmost of FN labels 1..n, excluding end */
			  int refid,
			  int verb)
{
  RFLOAT Bias = X * FrateBiasWT + (m-1)*FBiasWT + EBiasWT;
  if(verb){
    printf("Send(X=%0.3f,m=%d,n=%d,I=%d,refid=%d):Bias=%0.6f,FnPen= %0.6f, PfPen= %0.6f, FratePen= %0.6f\n",
	   X,m,n,I,refid,Bias,- (FN_ADJUST ? max(0.0, PRsum[refid][I] - PRsum[refid][I-n+1]) : (n-1)) * FnPenalty,-(m-1)*FpPenalty,- X*Frate);
    fflush(stdout);
  }
  return Bias - (FN_ADJUST ? max(0.0, PRsum[refid][I] - PRsum[refid][I-n+1]) : (n-1)) * FnPenalty -(m-1)*FpPenalty - X*Frate;
}

static inline RFLOAT BiasEnd(RFLOAT X, int m)
{
  if(DEBUG>=2) assert(m >= 1);
  return biasWTend * (X * FrateBiasWT + (m-1)*FBiasWT + EBiasWT);
}

/* NOTE : difference between BiasEnd2() and BiasEnd() is no longer needed since biasWTend can be used instead (use -biaswtEnd 0.0 to get effect of old BiasEnd2) */
static inline RFLOAT BiasEnd2(RFLOAT X,
			      int m)
{
  return BiasEnd(X,m);
  //  RFLOAT r = BiasEnd(X,m);
  //  return max((RFLOAT)0.0f, r);
}

/* external version of SintDetail : does not use tables */
void SintDetailRA(double X,
		  double Y,
		  int m,
		  int n,
		  int J,
		  int I,
		  int K,/**< right end of interval! */
		  int T,/**< left end of interval : only used with RES_VARIANCE */
		  FLOAT *Ya,/**< reference array Ya[0..N+1] */
		  double &Bias,
		  double &Pen,
		  double &Gauss,
		  double &PenSm,
		  int verb)
{
  if(DEBUG>=2) assert(0 <= K && K <= KMAX);
  if(DEBUG>=2) assert(I-K-n >= 1);
  if(DEBUG>=2 && RES_VARIANCE){
    assert(0 <= T && T <= KMAX);
    assert(I-K-n-T >= 1);
  }

  double err = X-Y;
  double var = F2var + S2var*Y;
  if(QUADRATIC_VARIANCE)
    var += R2var * Y * Y;
  if(RES_VARIANCE){
    double resR = Ya[I] - Ya[I-K];
    double resL = Ya[I-K-n] - Ya[I-K-n-T];
    double resvar = resL*resL + resR*resR;
    var += E2var * resvar;
  }
  if(DEBUG>=2) assert(var > 0.0);
  
  Bias = X * FrateBiasWT + m * FBiasWT;
  Gauss = RefPen - 0.5*mylog(var) - err*err/var;
  Pen = log(Pr(Ya[I-K]-Ya[I-K-n])) - (m-1)*FpPenalty;// WAS - X*Frate;
  double FnPen = 0.0;
  if(FN_ADJUST){
    double sum = 0.0;
    for(int i = 1; i <= n; i++)
      sum += Pr(Ya[I-i+1]-Ya[I-i]);
    FnPen = max(0.0,sum-1.0)*FnPenalty;
  } else
    FnPen = (n-1)*FnPenalty;
  Pen -= FnPen;

  PenSm = Sm(0,I,K,Ya);

  if(maptype){
    if(FRATE_FIX1)
      Gauss -= X*Frate;
    else 
      Pen -= X*Frate;
  } else {// maptype == 0
    if(FRATE_FIX0)
      Gauss -= X*Frate;
    else
      Pen -= X*Frate;
  }

  if(verb){
    printf("\n\t\tI=%d,K=%d,J=%d,n=%d,m=%d,T=%d,X=%0.4f,Y=%0.4f:SD=%0.6f,Gterm=%0.6f,ErrPen=%0.6f,X*Frate=%0.6f,FnPen=%0.6f,FpPen=%0.6f,X*LBias=%0.6f,m*FBias=%0.6f,PenPr=%0.6f,PenSm=%0.6f (err=%0.6f,sf=%0.4f,sd=%0.6f,sr=%0.6f,se=%0.6f)\n", 
	   I,K,J,n,m, T, X, Y, sqrt(var), RefPen - 0.5*mylog(var), -err*err/var, -X*Frate, -FnPen, -(m-1)*FpPenalty, X*LBias, m*FBias, log(Pr(Ya[I-K]-Ya[I-K-n])), Sm(0,I,K,Ya), fabs(err), SF[0],SD[0],SR[0],SE[0]);
    fflush(stdout);
  }
}

#endif // RGENTIGREFSCORE_H
