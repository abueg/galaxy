#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

#include <omp.h>

#include "globals.h"
#include "Assembler_parameters.h"

namespace probeval {

#ifndef LSCORE_H
#define LSCORE_H

static Ident RGentigScore_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/RGentigScore.h 10830 2020-04-05 22:24:03Z tanantharaman $");

#define LAMBDA_FIX outlierNormRef /* Don't adjust outlier penalty term normalization factor based on outlierLambda value
				     If >= 2 : don't adjust based on outlierMax either */

#define OUTLIER_LTYPE 1 /* OUTLIER exponential decay type :
			   0 : outlier likelihood is proportional to exp(-(x+y)/outlierLambda)
			   1 : outlier likelihood is proportional to exp(-|x-y|/outlierLambda)
			   NOTE : Currently only implemented with OUTLIER_LTYPE==1 and maptype==0
			*/

MFLOAT outlierMaxBC;/* OUTLIER_MAXBC ? min(OUTLIER_MAXBC, outlierMax) : outlierMax */

#define SCORE_Y 0 /* make score independent of X (less accurate Theta estimate and depends on Y, but result does not change when mapset is filtered) (Mprobeval not supported) */

#define RES_BND 0 /* Apply LBias only to interval size above resKB x 2 */

/* Gentig Refinement Likelihood scoring functions */

#define KBIAS KbiasRef /* to scale Gaussian penalty bias */
#define KLBIAS KLbiasRef /* to scale False cut penalty bias */
#define KFBIAS KFbiasRef /* to scale Missing cut penalty bias */
#define BETA 1.0 /* scale log(G/py) term */

int Extend = 20;/* default EXTEND_MAX : will be increased up to EXTEND, if needed, at the expense of run time (see refine.h) */

int *Kmax = 0;/* Kmax[1..N] */
int kmax = 0;
int KMAX_APPROX;
double maxdist;
double AlignResMax;

MFLOAT *VARtab = 0;/* VARtab[(K = 0..kmax) * N + (I=K+1..N)] = (Y[I]-Y[I-K])^2 */

Cprtab **PRtab = 0;/* PRtab[K=0..PRsiz][I=K+1..N] has precomputed penalties that depend on Y[I-K-1 .. I] */
double *PrBiasWT = 0;/* PrBiasWT[0..N] is precomputed value of Pr bias: PrBiasWT[I] = -PRbiasWT * biasWT * sum(i=0..I) Pr(yi) * log(Pr(yi)) + (1-Pr(yi)) * log(1-Pr(yi)), where yi == Y[i+1]-Y[i] */
// NOTE : using PrBiasWT breaks mprobeval's assumption that deleting a label I changes FAJ for an interval that includes I only by 1/FN[0] : 
//        Consider using the refalign() instead to avoid spurious endoutliers in high label density regions

/* precompute values to speed up likelihood score computation */
double minKB;
static double maxKB;

static double LBias, FBias/* , GBias */, EBias, FnPenalty, FpPenalty, Ygm, Xlambda, resKB2;
double Ylambda;
double Xtheta = 0.0;
double resKB;
float resKBf;
float IresSDf;
MFLOAT Frate;/* F + R - 1/Xtheta */
double Rtheta;/* sqrt(R * Xtheta) */
MFLOAT RthetaF;
double InvRtheta;/* 1/sqrt(R * Xtheta) */
MFLOAT InvRthetaF;
MFLOAT InvRtheta2;/* InvRtheta/Rtheta */
double pTP;/* 1-FN[0] */

MFLOAT FrateBiasWT;/* biasWT * (Frate + LBias) */
MFLOAT FBiasWT;/* biasWT * FBias */
// static MFLOAT GBiasWT;/* biasWT * GBias */
MFLOAT EBiasWT;/* biasWT * EBias */
MFLOAT OutlierLambdaInv;/* (outlierLambda < 1000.0) ? 1.0/outlierLambda : 0.0 */

MFLOAT F2var;/* 2*SF[0]*SF[0] */
MFLOAT S2var;/* 2*fabs(SD[0])*SD[0] */
MFLOAT R2var;/* 2*SR[0]*SR[0] */
MFLOAT E2var;/* 2*SE[0]*SE[0] */
MFLOAT RefPen;/* Xtheta/sqrt(PI) */
double IresSD;/* 1.0/(resSD * PixelLen * sqrt(2.0)) */
MFLOAT PoutlierM1;/* 1.0 - Poutlier */

MFLOAT PoutlierXthetaBYlambda; /* Poutlier                                                                             IF outlierLambda >= 1000.0 || LAMBDA_FIX
				  Poutlier * Xtheta / (2 * outlierLambda)                                              IF outlierLambda < 1000.0  && !LAMBDA_FIX */
MFLOAT PoutlierXthetaBYmax;/* Poutlier                                                                                 IF LAMBDA_FIX >= 2
			      Poutlier * Xtheta/(2 * outlierMaxBC)                                                     IF outlierLambda >= 1000.0 || LAMBDA_FIX == 1
			      Poutlier * Xtheta / (2 * outlierLambda * (1.0 - exp(-outlierMaxBC/outlierLambda)))       IF outlierLambda < 1000.0  && !LAMBDA_FIX */

MFLOAT PoutlierPIbXtheta;/* PoutlierXthetaBYlambda / RefPen */
MFLOAT PoutlierPIBYmax;/* PoutlierXthetaBYmax / RefPen */

MFLOAT PoutlierEndM1;/* 1.0 - PoutlierEnd */
MFLOAT biasWToutlierM1;/* 1.0 - biasWToutlier */
MFLOAT biasWToutlierF;/* biasWToutlier */

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT 
 __m256 ZeroHigh_mm256[DELTA_MAX + EXTEND + 1];

 __m256 v256_F2var, v256_S2var, v256_R2var, v256_E2var, v256_RefPen, v256_zero, v256_0p5, v256_1p0, v256_1p5, v256_2p0;
 __m256 v256_biasWToutlierF, v256_biasWToutlierM1;
 __m256 v256_Frate, v256_PoutlierM1, v256_FrateBiasWT;
 __m256 v256_OutlierLambdaInv, v256_PoutlierXthetaBYlambda, v256_PoutlierXthetaBYmax, v256_PoutlierPIbXtheta, v256_PoutlierPIBYmax;
 __m256 v256_outlierNegMaxBC;

 __m256 v256_logMaxScale, v256_logMaxUR, v256_logMaxFAG, v256_logMaxFAG2;
 __m256 v256_RthetaF, v256_InvRthetaF, v256_InvRtheta2;
 
 __m256d v256d_0p5, v256d_1p0;

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT
 __mmask16 MaskLow_mask16[DELTA_MAX + EXTEND + 1];

 __m512 v512_F2var, v512_S2var, v512_R2var, v512_E2var, v512_RefPen, v512_1p0, v512_0p5, v512_1p5, v512_2p0;
 __m512d v512d_0p5, v512d_1p0;
 __m512 v512_biasWToutlierF, v512_biasWToutlierM1;
 __m512 v512_Frate, v512_PoutlierM1, v512_FrateBiasWT;
 __m512 v512_OutlierLambdaInv, v512_PoutlierXthetaBYlambda, v512_PoutlierXthetaBYmax, v512_PoutlierPIbXtheta, v512_PoutlierPIBYmax;
 __m512 v512_outlierNegMaxBC, v512_outlierMaxBC;

 __m512 v512_logMaxScale, v512_logMaxUR, v512_logMaxFAG, v512_logMaxFAG2;
 __m512 v512_LOG2E, v512_MIN_EXP2ARG, v512_MAX_EXP2ARG;

 __m512 v512_RthetaF, v512_InvRthetaF;
#endif

int PRsiz = 0;
int PRNsiz = 0;

int gN = 0;/* value of N when score_init() was called */

double *FnPow=0;/* FnPow[i=0..N] = pow(1-p,i) */
MFLOAT *FnPowf=0;/* FnPowf[i=0..N] = pow(1-p,i) */
MFLOAT *LogFnPow=0;/* LogFnPow[i=0..N] = i * log(1-p) */

double FpPow[DELTA_MAX + EXTEND + 1];/* FpPow[i=0..max(DELTA_X,deltaExtXRef) + EXTEND_MAX] == pow(F*Xtheta,i) */
MFLOAT LogFpPow[(DELTA_MAX + EXTEND + 1 + VPAD)*2 + VPAD];/* LogFpPow[i=0..max(DELTA_X,deltaExtXRef) + EXTEND_MAX] == i * log(F*Xtheta) */
MFLOAT FBiasWTm[(DELTA_MAX + EXTEND + 2 + VPAD)*2 + VPAD];/* FBiasWTm[m=0..max(DELTA_X,deltaExtXRef) + EXTEND_MAX + 1] == m * FBiasWT */

MFLOAT *LogFpPowRev=0; /* LogFpPow[i=0..-max(DELTA_X,deltaExtXRef) + EXTEND_MAX] == -i * log(F*Xtheta) */
MFLOAT *FBiasWTmRev=0; /* FBiasWTm[m=0..-max(DELTA_X,deltaExtXRef) + EXTEND_MAX + 1] == -m * FBiasWT */

static void score_free()
{
  if(PRtab){
    for(register int K = 0; K <= PRsiz; K++)
      delete [] PRtab[K];
    delete [] PRtab;
    PRtab = 0;
    delete [] PrBiasWT;
    PrBiasWT = 0;
  }
  if(FnPow){
    delete [] FnPow;
    FnPow = 0;
  }
  if(FnPowf){
    delete [] FnPowf;
    FnPowf = 0;
  }
  if(LogFnPow){
    delete [] LogFnPow;
    LogFnPow = 0;
  }
  if(Kmax){
    delete [] Kmax;
    Kmax = 0;
  }
  if(VARtab){
    delete [] &VARtab[1];
    VARtab = 0;
  }
}

/* score_init()  updates caused by change in Poutlier, PoutlierEnd OR biasWT only (if PRbiasWT==0) */
void score_initOutlier()
{
  FrateBiasWT = biasWT * (Frate + LBias);
  FBiasWT = biasWT * FBias;
  //  GBiasWT = biasWT * GBias;
  EBiasWT = biasWT * EBias;

#if 1
  Poutlier = min(0.999999,Poutlier);
  PoutlierEnd = min(0.999999,PoutlierEnd);
#endif

  PoutlierM1 = 1.0 - Poutlier;
  PoutlierEndM1 = 1.0 - PoutlierEnd;
  
  biasWToutlierM1 = 1.0 - biasWToutlier;
  biasWToutlierF = biasWToutlier;

  outlierMaxBC = (OUTLIER_MAXBC ? min(OUTLIER_MAXBC, (double)outlierMax) : outlierMax);

  if(outlierLambda >= 1000.0)
    outlierLambda = 1e+30;
  if(LAMBDA_FIX || outlierLambda >= 1000.0){
    PoutlierXthetaBYlambda = Poutlier;
    if(LAMBDA_FIX >= 2)
      PoutlierXthetaBYmax = Poutlier;
    else
      PoutlierXthetaBYmax = Poutlier * Xtheta/(2.0 * outlierMaxBC);
  } else {
    PoutlierXthetaBYlambda = Poutlier * Xtheta / ( 2.0 * outlierLambda);
    PoutlierXthetaBYmax = Poutlier * Xtheta/(2.0 * outlierLambda * (1.0 - exp(-outlierMaxBC/outlierLambda)));
  }
  OutlierLambdaInv = 1.0/outlierLambda;
  if(DEBUG && !(isfinite(PoutlierXthetaBYlambda))){
    printf("score_init:outlierLambda= %0.6e, Poutlier=%0.6e, Xtheta=%0.6e:PoutlierXthetaBYlambda=%0.6e\n",
	   outlierLambda,Poutlier,Xtheta,PoutlierXthetaBYlambda);
    fflush(stdout);
    assert(isfinite(PoutlierXthetaBYlambda));
  }
  if(DEBUG) assert(isfinite(PoutlierXthetaBYmax));

  PoutlierPIbXtheta = PoutlierXthetaBYlambda / RefPen;
  if(DEBUG) assert(isfinite(PoutlierPIbXtheta));
  PoutlierPIBYmax = PoutlierXthetaBYmax / RefPen;
  if(DEBUG) assert(isfinite(PoutlierPIBYmax));

#if (VECTOR_AVX || VECTOR2_AVX) && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT
  v256_biasWToutlierF = _mm256_broadcast_ss(&biasWToutlierF);
  v256_biasWToutlierM1 = _mm256_broadcast_ss(&biasWToutlierM1);

  v256_FrateBiasWT = _mm256_broadcast_ss(&FrateBiasWT);
  v256_PoutlierM1 = _mm256_broadcast_ss(&PoutlierM1);

  v256_OutlierLambdaInv = _mm256_broadcast_ss(&OutlierLambdaInv);
  v256_PoutlierXthetaBYlambda = _mm256_broadcast_ss(&PoutlierXthetaBYlambda);
  v256_PoutlierXthetaBYmax = _mm256_broadcast_ss(&PoutlierXthetaBYmax);
  v256_PoutlierPIbXtheta = _mm256_broadcast_ss(&PoutlierPIbXtheta);
  v256_PoutlierPIBYmax = _mm256_broadcast_ss(&PoutlierPIBYmax);

  float tmp = outlierMaxBC;
  v256_outlierNegMaxBC = _mm256_set1_ps(-tmp);

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

  v512_biasWToutlierF = _mm512_set1_ps(biasWToutlierF);
  v512_biasWToutlierM1 = _mm512_set1_ps(biasWToutlierM1);
  v512_FrateBiasWT = _mm512_set1_ps(FrateBiasWT);
  v512_PoutlierM1 = _mm512_set1_ps(PoutlierM1);

  v512_OutlierLambdaInv = _mm512_set1_ps(OutlierLambdaInv);
  v512_PoutlierXthetaBYlambda = _mm512_set1_ps(PoutlierXthetaBYlambda);
  v512_PoutlierXthetaBYmax = _mm512_set1_ps(PoutlierXthetaBYmax);
  v512_PoutlierPIbXtheta = _mm512_set1_ps(PoutlierPIbXtheta);
  v512_PoutlierPIBYmax = _mm512_set1_ps(PoutlierPIBYmax);

  float tmp = outlierMaxBC;
  v512_outlierMaxBC = _mm512_set1_ps(tmp);
  v512_outlierNegMaxBC = _mm512_set1_ps(-tmp);
#endif
  
}

/* must be called initially and then whenever Y, X[0..MD-1] or error parameters change */
/* HERE : compute special case when only Y changed since last call, since many tables and constants do NOT depend on Y */
static void score_init(int N, double *Y, int MD, int* MX, double **X)
{
  if(colors != 1){
    printf("score_init(): colors=%d not supported\n",colors);
    exit(1);
  }

  if(Beta != 1.0){
    printf("Beta != 1.0 no longer supported\n");
    exit(1);
  }
  if(BETA != 1.0){
    printf("BETA != 1.0 no longer supported\n");
    exit(1);
  }
  if(maptype != 0){
    printf("WARNING:refinement not supported with maptype=%d: using -maptype 0\n",maptype);
    fflush(stdout);
    maptype = 0;
  }

  if(OUTLIER_TYPE != OutlierType){
    printf("WARNING: -outlierType %d ignored : using fixed value %d instead (see OUTLIER_TYPE in RGentigRefScore.h)\n", OutlierType, OUTLIER_TYPE);
    fflush(stdout);
  }

  if(!(OUTLIER_LTYPE==1 && maptype==0)){
    printf("RGentigScore.h: Refinement not implemented except IF maptype==0 AND OUTLIER_LTYPE==1\n");
    exit(1);
  }
  if(outlierBC){
    printf("WARNING:-outlierBC ignored during refinement\n");
    fflush(stdout);
  }
  if(outlierMax >= 1000.0 && OUTLIER_MAXBC > 0.0){
    printf("-outlierMax must be <= %0.1f kb (change OUTLIER_MAXBC in probeval.h to support larger values)\n",OUTLIER_MAXBC);
    exit(1);
  }
  /*  if(outlierLambda < 1000.0){
    printf("-outlierLambda %0.3f not yet implemented during refinement\n",outlierLambda);
    exit(1);
    }*/

  resKBf = resKB = res[0] * PixelLen;
  resKB2 = (RES_BND ? 2.0 * resKB : 0.0);

  double pFN = FN[0];
  pTP = 1.0-FN[0];
  double F = FP[0]*0.01;
  double var = fabs(SD[0])*SD[0];
  double varF = max(SF[0]*SF[0], 1e-10);
  double varR = SR[0]*SR[0];
  double varE = SE[0]*SE[0];

  if(!PRtab){ /* compute Ylambda and Ygm (mean and geometric mean of reference intervals, excluding intervals below resKB)) */
    Ygm = 0.0;/* to compute Geometric Mean of Y intervals */
    int Ycnt = 0;
    double Ysum = 0.0;
    for(int i= 1; i < N; i++){
      double y = Y[i+1]-Y[i];
      Ysum += y;
      if(y > resKB){
	Ygm += log(y);
	Ycnt++;
      }
    }
    Ygm /= max(1,Ycnt);
    Ygm = exp(Ygm);

    Ylambda = Ysum / max(1,Ycnt);
    //  Ylambda = (Y[N] - Y[1])/max(1,Ycnt);
    //  Ylambda = (Y[N] - Y[1])/((N>1) ? N-1 : 1);
    if(DEBUG && !XmapCount[0] && !(Ylambda > 0.0)){
      printf("Ylambda=%0.3f: Ycnt=%d,Ysum=%0.3f,N=%d,resKB=%0.3f:\n",Ylambda,Ycnt,Ysum,N,resKB);
      for(int i = 1; i <= N+1; i++)
	printf("Y[%d]=%0.3f\n",i,Y[i]);
      fflush(stdout);
      assert(Ylambda > 0.0);
    }
  }

  /* compute theta based on X */
  double theta = 0.0;
  int sitecnt = 0;

  if(XmapCount[0] > 0){/* use command line values for XmapCount,XmapSites,XmapLength */
    MD = XmapCount[0];
    sitecnt = XmapSites[0];
    theta = XmapLength[0];
  } else {
    for(int m = 0; m < MD; m++){
      sitecnt += MX[m];
      theta += X[m][MX[m]+1];
    }
  }
  if(sitecnt <= 0 || theta <= 0.0 || MD <= 0){
    if(MD <= 0 && Refine <= 1){/* If Refine <= 1, no alignments will be performed with MD == 0 */
      printf("score_init(): Insufficcient maps:MD=%d,total sites=%d, total length=%0.3f: faking it with theta= 10kb and sitecnt = 1\n",
	     MD,sitecnt,theta);
      fflush(stdout);

      theta = 10.0;
      sitecnt = 1;
    } else {
      printf("score_init(): Insufficcient maps:MD=%d,total sites=%d, total length=%0.3f\n",
	     MD,sitecnt,theta);
      fflush(stdout);
      exit(1);
    }
  }
  double R = RTHETA_FIX ? (MD / theta) : 0.0;/* NOTE : MD/theta is inverse of the average molecule size. This should be reduced by inverse of average linear chromosome size */

  Xtheta = theta / sitecnt - resKB2*0.5;
  if(VERB>=2){
    printf("%d maps:total sitecnt=%d, len=%0.6f kb, resKB2=%0.6f, Xtheta=%0.6f\n",MD,sitecnt,theta,resKB2, Xtheta);
    fflush(stdout);
  }
  Xlambda = pTP/(1.0/(Xtheta - var*0.5) - F);/* HERE : does not include correction for SR[] */

  if(XmapCount[0]){
    Ylambda = Xlambda;
    Ygm = XmapGtheta[0] * Xlambda/Xtheta;
  } else if(SCORE_Y || MD <= 0){ /* make score dependent only on Y (not X) */
    Xlambda = Ylambda;
    Xtheta = 1.0/(F + pTP/Xlambda) + var*0.5;
  }

  if(DEBUG){
    assert(FN[0] > 0.0);
    if(!(FN[0] < 1.0)){
      printf("score_init:FN[0]=%0.6f\n",FN[0]);
      fflush(stdout);
      assert(FN[0] < 1.0);
    }
    assert(theta > 0.0);
    assert(Ylambda > 0.0);
    assert(Xlambda > 0.0);
    assert(Xtheta > 0.0);
  }

  Frate = F + R - thetaScaleRef/Xtheta;
  RthetaF = Rtheta = RTHETA_FIX ? sqrt(R * Xtheta) : 1.0;
  InvRthetaF = InvRtheta = RTHETA_FIX ? 1.0/Rtheta : 1.0/(R * Xtheta);// WAS98 (RTHETA_FIX ? 1.0/(R*Xtheta) : 1.0/Rtheta)
  InvRtheta2 = InvRtheta/Rtheta;

  if(F*Xtheta >= 1.0){
    printf("FP density=%0.4f, Xtheta=%0.3f: FP is too large for meaningful alignment scoring\n",F,Xtheta);
    printf("var=%0.6f,pTP=%0.6f,Ylambda=%0.6f,resKB2=%0.6f,Xtheta=%0.6f,sitecnt=%d\n",var,pTP,Ylambda,resKB2,theta/sitecnt,sitecnt);
    exit(1);
  }

  double YgmBpTP = (Ygm/pTP) * (Xlambda/Ylambda);
  double Sbias = log(Xtheta) - 0.5*log(2.0*M_PI*(varF+var*YgmBpTP +varR*YgmBpTP*YgmBpTP));
  FBias = 0.5*KBIAS - KFBIAS*pFN*log(pFN)/pTP - log(pTP) - Sbias; /* Bias per site in X (overlapping Y)  = Abias - Sbias in document */
  LBias = -KLBIAS * F*log(F*Xtheta) - F*FBias;/* Bias per kb in X (overlapping Y) : In document LBias is larger by Frate = F - 1/Xtheta (see FrateBiasWT) */
  EBias = 0.5*(Sbias + FBias);/* Bias per aligned end , corresponding to Abias/2 in the document. This is larger by KBIAS/4 compared to the document, so the final
			       score will be KBIAS/2 larger than in the document */
  //  GBias = 0.5*KBIAS - Sbias;/* Bias per aligned interval : part of FBias */

  /* Note : Unlike in RGentigRefScore.h, Bias is not applied to unaligned EndOutliers, to help detect chimeric regions (if biasWT > 0.0). Also biasWT is based on -biaswtRefine instead of -biaswt */

  if(VERB && (VERB>=2 || rverb || !PRtab || ERRPLOT)){
    printf("K=%0.2f,KL=%0.2f,KF=%0.2f,B=%0.2f,pTP=%0.4f,F=%0.8f,R=%0.8f,Xtheta=%0.6f,Ylambda=%0.6f(Xl=%0.6f),Ygm=%0.6f,Frate=%0.8f,N=%d,Sbias=%0.6f\n",
	   KBIAS,KLBIAS,KFBIAS,BETA,pTP,F,R,Xtheta,Ylambda,Xlambda,Ygm,Frate,N,Sbias);
    fflush(stdout);
  }

  FnPenalty = -log(pFN);/* missing cut penalty*/
  FpPenalty = -log(F*Xtheta);/* false cut penalty */

  minKB = ZERO_MINKB ? ZERO_MINKB /* WAS45 0.001 * PixelLen */ : rres * 0.500;
  maxKB = (res[0] + resSD[0] * SDrange)*PixelLen;

  F2var = 2.0 * varF;
  S2var = 2.0 * var;
  R2var = 2.0 * varR;
  E2var = 2.0 * varE;
  RefPen = Xtheta/sqrt(M_PI);

  /* set parameters used for rescaling by qprobeval and mprobeval */
  logMaxScale = log(MaxScale);
  logMaxUR = log(MaxUR);
  double Vmin = SF[0]*SF[0] + SD[0]*SD[0]*SD[0]*fabs(SD[0])/(4.0*SR[0]*SR[0]);
  if(DEBUG) assert(Vmin > 0.0);
  Smin = sqrt(Vmin);
  Smin = min(SF[0],Smin);
  MaxFAG = FN[0] * MaxFA * sqrt(2.0) /(Smin * RefPen);
  logMaxFAG = log(MaxFAG);
  MaxFAG2 = FN[0] * MaxFA / RefPen;
  logMaxFAG2 = log(MaxFAG2);
  if(VERB && (VERB>=2 || !PRtab || rverb)){
    printf("\tSF=%0.8f,SD=%0.8f,SR=%0.8f: Vmin= %0.10f,Smin= %0.10f:MaxFAG=%0.8e,MaxFAG2=%0.8e\n",SF[0],SD[0],SR[0],Vmin,Smin,MaxFAG,MaxFAG2);
    fflush(stdout);
  }

#if (VECTOR_AVX || VECTOR2_AVX) && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT
  v256_F2var = _mm256_broadcast_ss(&F2var);
  v256_S2var = _mm256_broadcast_ss(&S2var);
  v256_R2var = _mm256_broadcast_ss(&R2var);
  v256_E2var = _mm256_broadcast_ss(&E2var);
  v256_RefPen = _mm256_broadcast_ss(&RefPen);
  float tmp = 0.0f;
  v256_zero = _mm256_broadcast_ss(&tmp);
  tmp = 0.5f;
  v256_0p5 = _mm256_broadcast_ss(&tmp);
  tmp = 1.0f;
  v256_1p0 = _mm256_broadcast_ss(&tmp);
  tmp = 1.5f;
  v256_1p5 = _mm256_broadcast_ss(&tmp);
  tmp = 2.0f;
  v256_2p0 = _mm256_broadcast_ss(&tmp);

  v256_Frate = _mm256_broadcast_ss(&Frate);

  v256_logMaxScale = _mm256_broadcast_ss(&logMaxScale);
  v256_logMaxUR = _mm256_broadcast_ss(&logMaxUR);
  v256_logMaxFAG = _mm256_broadcast_ss(&logMaxFAG);
  v256_logMaxFAG2 = _mm256_broadcast_ss(&logMaxFAG2);

  v256_RthetaF = _mm256_set1_ps(RthetaF);
  v256_InvRthetaF = _mm256_set1_ps(InvRthetaF);
  v256_InvRtheta2 = _mm256_set1_ps(InvRtheta2);

  v256d_0p5 = _mm256_set1_pd(0.5);
  v256d_1p0 = _mm256_set1_pd(1.0);

  // Initialize ZeroHigh_mm256[n = 0 .. DELTA_MAX+EXTEND] : float vector with first min(n,8) values 1.0f, all other high order values 0.0f */
  float tmp8[8] __attribute__ ((aligned(32)));
  for(int i = 0; i < 8; i++)
    tmp8[i] = 0.0f;
  for(int n = 0; n < 8; n++){
    ZeroHigh_mm256[n] = _mm256_load_ps(tmp8);
    tmp8[n] = 1.0f;
  }
  for(int n = 8; n <= DELTA_MAX+EXTEND; n++)
    ZeroHigh_mm256[n] = _mm256_load_ps(tmp8);

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

  v512_F2var = _mm512_set1_ps(F2var);
  v512_S2var = _mm512_set1_ps(S2var);
  v512_R2var = _mm512_set1_ps(R2var);
  v512_E2var = _mm512_set1_ps(E2var);
  v512_RefPen = _mm512_set1_ps(RefPen);
  float tmp = 0.5f;
  v512_0p5 = _mm512_set1_ps(tmp);
  tmp = 1.0f;
  v512_1p0 = _mm512_set1_ps(tmp);
  tmp = 1.5f;
  v512_1p5 = _mm512_set1_ps(tmp);
  tmp = 2.0f;
  v512_2p0 = _mm512_set1_ps(tmp);
  
  double tmpd = 0.5;
  v512d_0p5 = _mm512_set1_pd(tmpd);

  tmpd = 1.0;
  v512d_1p0 = _mm512_set1_pd(tmpd);

  v512_Frate = _mm512_set1_ps(Frate);

  v512_logMaxScale = _mm512_set1_ps(logMaxScale);
  v512_logMaxUR = _mm512_set1_ps(logMaxUR);
  v512_logMaxFAG = _mm512_set1_ps(logMaxFAG);
  v512_logMaxFAG2 = _mm512_set1_ps(logMaxFAG2);

  v512_LOG2E = _mm512_set1_ps((float)M_LOG2E);
  //  v512_MIN_EXP2ARG = _mm512_set1_ps(MIN_EXP2ARG);
  //  v512_MAX_EXP2ARG = _mm512_set1_ps(MAX_EXP2ARG);

  v512_RthetaF = _mm512_set1_ps(RthetaF);
  v512_InvRthetaF = _mm512_set1_ps(InvRthetaF);

  // Initialize MaskLow_mask16[n = 0 .. DELTA_MA+EXTEND] : 16 bit mask with lowest min(n,16) bits 0, all other high order bits 1 */
  for(int n = 0; n < 16; n++)
    MaskLow_mask16[n] = _mm512_int2mask(((~(0ul)) << n) & 0xffff);
  for(int n = 16; n <= DELTA_MAX + EXTEND; n++)
    MaskLow_mask16[n] = _mm512_int2mask(0);

#endif

  score_initOutlier();

  if(VERB && (VERB>=2 || rverb || !PRtab || ERRPLOT)){
    printf("LBias=%0.8f,EBias=%0.6f,FBias=%0.6f,FnPenalty=%0.6f,FpPenalty=%0.6f,sf=%0.6f,sd=%0.6f,sr=%0.6f,Poutlier=%0.4e,%0.4e(init=%d),outlierLambda=%0.3e(Type=%d),resKB=%0.6f,minKB=%0.6f,maxKB=%0.6f,PRbiasWT=%0.2f\n",
	   LBias,EBias,FBias,FnPenalty,FpPenalty,SF[0],SD[0],SR[0],Poutlier,PoutlierEnd,PoutlierEnd_init,outlierLambda,OutlierType,resKB,minKB,maxKB,PRbiasWT);
    fflush(stdout);
  }

  IresSDf = IresSD = 1.0/(resSD[0] * PixelLen * sqrt(2.0));

  if(VERB>=2){
    printf("score_init():N=%d,MD=%d,Xtheta=%0.6f,outlierMax=%0.6e,Poutlier=%0.6e,PoutlierXthetaBYmax=%0.6e,PoutlierPIBYmax=%0.6e\n",
	   N,MD,Xtheta,outlierMaxBC,Poutlier,PoutlierXthetaBYmax,PoutlierPIBYmax);
    fflush(stdout);
  }

  score_free();

  /* initialize Kmax[1..N] */
  Kmax = new int[N+1];
  KMAX_APPROX = (SCORE_APPROX >= (AlignRes ? 2 : 1) || AlignResMult <= 1.0) ? 1 : 0;
  maxdist = (res[0] + SDrange * resSD[0]) * PixelLen;
  AlignResMax = maxdist * AlignResMult;
  for(int I = 1; I <= N; I++){
    int K = 0;
    if(KMAX_APPROX){/* determine Kmax[I] based on total distance Y[I] - Y[I-K] < maxdist */
      while(I-K > 1){
	if(Y[I] - Y[I-K-1] >= maxdist)
	  break;
	K++;
      }
    } else {
      while(I-K > 1){
	if((!KMAX_FIX && Y[I-K] - Y[I-K-1] >= maxdist) || Y[I]-Y[I-K-1] >= AlignResMax)
	  break;
	K++;
      }
    }
    Kmax[I] = min(K,KMAX);
  }
  kmax = 0;
  for(int I = 1; I <= N; I++)
    kmax = max(kmax,Kmax[I]);

  if(RES_VARIANCE){
    VARtab = new MFLOAT[(kmax+1) * N];
    VARtab--;
    for(int K = 0; K <= kmax; K++){
      MFLOAT *p = &VARtab[K*N];
      for(int I = 1+K; I <= N; I++){
	double res = Y[I] - Y[I-K];
	p[I] = res * res;
      }
    }
  }
  

  /* allocate PRtab[K=0..max(max(2,KMAX),max(DELTA_Y,deltaExtYRef) + (DELTA_MIN ? DELTA_X-DELTA_MIN:0)][I=K+1..N] */
  PRsiz = max(max(2,kmax), max(DELTA_Y,deltaExtYRef) + (DELTA_MIN ? DELTA_LIM : 0)+1/*for deleted site*/ -1);
  PRNsiz = N+1/*for added site*/;/* for debugging */
  PRtab = new Cprtab*[PRsiz+1];
  if(VERB>=2){
    printf("score_init:PRsiz=%d,N=%d\n",PRsiz,N);
    fflush(stdout);
  }

  int numthreads = 1;
  #ifdef _OPENMP
  numthreads = omp_get_max_threads();
  if(MaxThreads > 0)
    numthreads = min(numthreads,MaxThreads);
  numthreads = min(numthreads,max(1,N/16));
  #endif

  if(VERB>=3){
    printf("score_init: computing PRtab[K=0..%d][I=K+1..%d] = Pr(Y[I] - Y[I-K-1])\n",PRsiz,N+1);
    fflush(stdout);
  }

  delete [] PrBiasWT;
  PrBiasWT = new double[N+2];
  double cum = 0.0;
  PrBiasWT[0] = cum;
  double PRbiaswt = PRbiasWT * biasWT;

  #pragma omp parallel for num_threads(numthreads) schedule(static,16)
  for(int I = 1; I <= N+1; I++){
    double y = Y[I] - Y[I-1];
    double pr = Pr(y);
    if(DEBUG>=2) assert(isfinite(pr) && 0.0 <= pr && pr <= 1.0);
    double prbias = (pr <= 0.0 || pr >= 1.0) ? 0.0 : pr * log(pr) + (1.0 - pr) * log(1.0 - pr);
    if(DEBUG>=2 && !(isfinite(prbias) && prbias <= 0.0)){
      printf("I=%d,N=%d:Y[I]=%0.6f,Y[I+1]=%0.6f,y=%0.6f,pr=%0.8e,prbias=%0.8e\n",I,N,Y[I],Y[I+1],y,pr,prbias);
      fflush(stdout);
      assert(isfinite(prbias) && prbias <= 0.0);
    }
    PrBiasWT[I] = cum -= PRbiaswt * prbias;
    if(DEBUG>=2) assert(isfinite(cum) && cum >= 0.0);
  }

  for(int K = 0; K <= kmax; K++){
    int n = K+1;
    PRtab[K] = new Cprtab[N+1+/*for added site*/1];

    #pragma omp parallel for num_threads(numthreads) schedule(static,16)
    for(int I = K+1; I <= N+/*for added site*/1;I++){
      Cprtab *p = &PRtab[K][I];
      double y = Y[I] - Y[I-n];
      p->Pr = Pr(y);
      if(DEBUG>=2) assert(isfinite(p->Pr) && 0.0 <= p->Pr && p->Pr <= 1.0);

      if(VERB>=2 && I==3199 && n==1){
	#pragma omp critical
	{
	  printf("I=%d,N=%d,K=%d,n=%d:Y[I]-Y[I-n]= y= %0.10f, Pr(y)= %0.10e, p->Pr= %0.10e\n",I,N,K,n,y,Pr(y),p->Pr);
	  fflush(stdout);
	}
      }

      p->PM = PM(0,I,K,Y);
      if(DEBUG>=2 && !(isfinite(p->PM) && p->PM >= 0.0)){
        printf("I=%d,K=%d,N=%d,Y[I]=%0.3f,Y[I-K]=%0.3f:PM(0,I,K,Y)=%0.8e,p->PM=%0.8e\n",
          I,K,N,Y[I],Y[I-K],PM(0,I,K,Y),p->PM);
        for(register int i = I-K;i <= I; i++)
          printf("Y[%d]=%0.6f\n",i,Y[i]);
	fflush(stdout);
	assert(isfinite(p->PM) && p->PM >= 0.0);
      }
    }
  }

  for(int K = kmax+1; K <= PRsiz; K++){
    int n = K+1;
    PRtab[K] = new Cprtab[N+1+/*for added site*/1];

    #pragma omp parallel for num_threads(numthreads) schedule(static,16)
    for(int I = K+1; I <= N+/*for added site*/1;I++){
      Cprtab *p = &PRtab[K][I];
      double y = Y[I] - Y[I-n];
      p->Pr = Pr(y);
      if(DEBUG>=2) assert(isfinite(p->Pr) && 0.0 <= p->Pr && p->Pr <= 1.0);
    }
  }

  /* precompute powers of (1-p) and F*Xtheta */
  gN = N;
  FnPow = new double[N+1];
  FnPowf = new MFLOAT[N+1];
  LogFnPow = new MFLOAT[N+1];
  double LR = 1.0;
  double LogFN = log(pFN);

  for(int i= 0; i <= N; i++){// NOTE : cannot be multithreaded due to sequential dependency
    FnPowf[i] = FnPow[i] = LR;
    LR *= pFN;
    LogFnPow[i] = i * LogFN;
  }

  LR = 1.0;
  double FXtheta = F*Xtheta;
  double LogFXtheta = log(FXtheta);
  double FBiasWTd = FBias * biasWT;
  
  LogFpPowRev = &LogFpPow[(DELTA_MAX + EXTEND + 1 + VPAD)*2-1];
  FBiasWTmRev = &FBiasWTm[(DELTA_MAX + EXTEND + 2 + VPAD)*2-1];

  for(int i = 0; i <= max(DELTA_X,deltaExtXRef) + EXTEND_MAX + VPAD; i++){
    FpPow[i] = LR;
    LR *= FXtheta;
    LogFpPowRev[-i] = LogFpPow[i] = i * LogFXtheta;
    FBiasWTmRev[-i] = FBiasWTm[i] = i * FBiasWTd;
  }
  memset(LogFpPowRev, 0, VPAD * sizeof(MFLOAT));
  memset(FBiasWTmRev, 0, VPAD * sizeof(MFLOAT));
}

#endif

} // namespace probeval
