#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

// #include <immintrin.h>

#include "constants.h"
#include "parameters.h"
#include "globals.h"

#ifndef NGENTIGPAIRSCORE_H
#define NGENTIGPAIRSCORE_H

static Ident NGentigPairScore_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/NGentigPairScore.h 11381 2020-07-28 23:26:39Z tanantharaman $");

#define OUTLIER_TYPE 0 /* OUTLIER type :
			  0 : Include both misaligned sites and sizing error in outlier (requires DELTA_X = DELTA_Y >= n to evaluate exactly)
			  1 : Include only sizing error in outlier (HERE : NOT YET IMPLEMENTED : see RGentigScore.h for corresponding implementation) */

#define OUTLIER_DELTA(err) (!OUTLIER_MAX || fabs(err) < ((RFLOAT)outlierMax)) // WAS (!OUTLIER_MAX || (((RFLOAT)-outlierMax) < err && err < ((RFLOAT)outlierMax)))

#define OUTLIER(x,y) OUTLIER_DELTA(x-y)

#define SDFIX 1 /* WAS (PairsplitRef >= 2 ? 1 : 0) */ /* correctly use -ve sd values */

#define RES_BND 1 /* Apply LBias only to interval size above resKB x 2 */

#define HBIAS 1.0 /* 0.0 is original R based score, 1.0 corresponds more closely to GentigPairScore.h */

#define KBIAS (Kbias * biasWT)
#define KLBIAS (KLbias * biasWT)
#define KFBIAS (KFbias * biasWT)

#define BETA 0.0 /* Beta */ /* scale 0.5*log(varF+var*(X+Y)+varR*(X^2+Y^2)) and log(X+Y) terms and corresponding Bias term */

/* Gentig Pairwise Alignment scoring functions */

/* precompute values to speed up likelihood score computation */
static double MatchScore, Ylambda, Ytheta, Glambda, MinScore;
static double OutlierBias, OutlierEndPenalty, OutlierEndBias, MatchBiasD;
static PFLOAT MatchBias;
static PFLOAT biasWToutlierF; // biasWToutlier

static PFLOAT EndThresh, EndThresh2;/* Right End thresholds to avoid calling Send(),Sbnd() */
static FLOAT LBiasD, LBiasD2, FBiasD, FBiasD_2, MissPenaltyD, varD, varFD, varRD, resKB2D, EndBiasD;
static PFLOAT LBias, LBias2, FBias, MissPenalty, Frate, OutlierPenalty, FLBias, var, varF, varR, resKB2, EndBias;

static PFLOAT MatchBiasOutlier;/* MatchBias + OutlierBias */
static FLOAT varFD2,varD2,varRD2;/* double of varFD,varD,varRD */
static PFLOAT varF2,var2,varR2;/* double of varF,var, varR */
static PFLOAT MatchScore2;/* MatchScore + 2*MissPenalty */

static long long maxM,maxN;

static PFLOAT *OutlierPenaltyBC = NULL;

#if USE_SSE && USE_PFLOAT
static __m128 v_resKB2, v_MissPenalty, v_Frate, v_FLBias, v_var2, v_varF2, v_varR2, v_MatchScore2, v_MatchBiasOutlier, v_OutlierPenalty, v_outlierMax, v_outlierMaxneg, v_biasWToutlier, v_MatchBias, v_OutlierBias; // v_FBias;
#endif

#if USE_AVX && USE_PFLOAT
static __m256 

v256_resKB2, v256_MissPenalty, v256_FBias, v256_EndBias, v256_Frate, v256_FLBias, v256_var2, v256_varF2, v256_varR2, v256_MatchScore2, v256_MatchBiasOutlier, v256_OutlierPenalty, v256_outlierMax, v256_outlierMaxneg, v256_PZERO, v256_biasWToutlier, v256_MatchBias, v256_OutlierBias;
		
#endif

		
#if USE_MIC && USE_PFLOAT
static __m512 	v512_resKB2, v512_MissPenalty, v512_FBias, v512_EndBias, v512_FBiasD_2, v512_Frate, v512_LBias, v512_LBiasD2, v512_FLBias, v512_var, v512_var2, v512_varF, v512_varF2, v512_varR, v512_varR2, v512_MatchScore2, v512_MatchBiasOutlier, v512_OutlierPenalty, v512_outlierMax, v512_outlierMaxneg, v512_PZERO, v512_biasWToutlier, v512_MatchBias, v512_OutlierBias;


static __m512 v512_ascending0123;
static __m512i v512i_ascending0123;

#endif // USE_MIC && USE_PFLOAT

static void score_init(int mapstart, int mapend)
{
  if(mapend < mapstart)
    return;

  if(OUTLIER_TYPE != 0){
    printf("NGentigPairScore.h:OUTLIER_TYPE==1 not yet implemented\n");
    fflush(stdout);exit(1);
  }
  if(outlierLambda <= 1000.0){
    printf("-outlierLambda %0.3f not implemented for pairwise alignment\n",outlierLambda);
    fflush(stdout);exit(1);
  }
  if(RES_BND != 1){
    printf("RES_BND=%d no longer supported\n",RES_BND);
    fflush(stdout);exit(1);
  }

  if(Beta != 0.0){
    printf("Beta=%0.3f: -Beta != 0.0 not supported\n",Beta);
    fflush(stdout);exit(1);
  }

  if(DEBUG) assert(colors == 1);

  biasWToutlierF = biasWToutlier;

  double pFN = FN[0];
  double pTP = 1.0 - pFN;
  if(DEBUG) assert(pFN < 1.0);
  if(DEBUG) assert(pTP > 0.0);
  double F = FP[0]*0.01;
  var = varD = (SDFIX ? fabs(SD[0]) : SD[0]) * SD[0];
  varF = varFD = SF[0]*SF[0]*2.0;
  varR = varRD = SR[0]*SR[0];

  /* compute theta, the average observed interval size */
  double Gtheta = 0.0;
  double theta = 0.0;
  int sitecnt = 0;
  maxM = 0;
  for(int i=mapstart; i <= mapend;i++){
    Cmap *Xmap = Gmap[i];
    int M = Xmap->numsite[0];
    sitecnt += M;
    theta += Xmap->site[0][M+1];
    if(M > maxM)
      maxM = M;
  }
  if(XmapCount[0] > 0){/* use command line values for XmapCount,XmapSites,XmapLength,Gtheta */
    //    int nummaps = XmapCount[0];
    sitecnt = XmapSites[0];
    theta = XmapLength[0];
    Gtheta = XmapGtheta[0];
  } else {
    /* compute Gtheta, the geometric mean of the observed interval size */
    for(int i=mapstart; i <= mapend;i++){
      Cmap *Xmap = Gmap[i];
      int M = Xmap->numsite[0];
      FLOAT *X = Xmap->site[0];
      // END fragments are ignored since they are often truncated to 0 
      //    Gtheta += log(X[1]+X[M+1]-X[M]);/* treat end frags as one sample */
      for(int i=1;i<M;i++)
	Gtheta += log(X[i+1]-X[i]);
    }
    Gtheta /= (sitecnt - (mapend-mapstart+1));
    Gtheta = exp(Gtheta);
  }
  
  //  maxM = (maxM+0xf) & ~0xf;
   
  maxN = maxM;
  if(DEBUG) assert(maxN <= MAXINT);
  if(VERB>=2){
    printf("score_init:maxN=maxM = %lld\n",maxM);
    fflush(stdout);
  }

  if(sitecnt <= 0 || theta <= 0.0){
    printf("score_init():Insufficient data:mapstart=%d,mapend=%d,sitecnt=%d\n",
	   mapstart,mapend,sitecnt);
    fflush(stdout);exit(1);
  }

  theta /= sitecnt;

  double resKB = (RES_BND ? mres * PixelLen : 0.0);
  resKB2 = resKB2D = 2.0 * resKB;

  double lambda = pTP/(1.0/(theta - resKB - varD*0.5) - F);/* HERE : does not include correction for SR[] */
  if(VERB>=2){
    printf("pTP=%0.17f,theta=%0.17f,resKB=%0.17f,varD=%0.17f(sd=%0.10f),F=%0.17f:lambda=%0.10f\n",
	   pTP,theta,resKB,varD,SD[0],F,lambda);
    fflush(stdout);
  }
  if(lambda <= 0.0){
    printf("ERROR:theta=%0.6f,pTP=%0.5f,resKB=%0.6f,varD=%0.6f,F=%0.6f:lambda=%0.6f!\n",
	   theta,pTP,resKB,varD,F,lambda);
    printf("      mres=%0.6f,PixelLen=%0.6f\n",mres,PixelLen);
    if(1.0/theta < F)
      printf("      False Positive rate exceeds average site density\n");
    printf("Error parameters inconsistent with average site density\n");
    fflush(stdout);exit(1);
  }
  Glambda = (Gtheta/theta)*lambda;

  Ylambda = lambda;
  Ytheta = theta;

  if(F*theta >= 1.0){
    printf("FP density=%0.4f, theta=%0.3f: FP is too large for meaningful alignment scoring\n",F,theta);
    fflush(stdout);exit(1);
  }

  if(DEBUG)assert(Poutlier >= 0.0 && Poutlier <= 1.0);
  if(Poutlier <= 0.0){
    OutlierPenalty = LARGE_NEGATIVE;
    OutlierBias = 0.0;
  } else {
    OutlierPenalty = log(Poutlier/max(0.0001,1.0-Poutlier));
    OutlierBias = log(max(0.0001,1.0-Poutlier));
  }
  if(PoutlierEnd <= 0.0){
    OutlierEndPenalty = LARGE_NEGATIVE;
    OutlierEndBias = 0.0;
  } else {
    OutlierEndPenalty = log(PoutlierEnd/max(0.0001,1.0-PoutlierEnd));
    OutlierEndBias = log(max(0.0001,1.0-PoutlierEnd));
  }

  if(VERB){
    printf("Kbias=%0.2f,KLbias=%0.2f,KFbias=%0.2f(biasWT=%0.4f,biasWToutlier=%0.4f),Beta=%0.2f,HBIAS=%0.2f,pTP=%0.6f,F=%0.6f,theta=%0.6f,lambda=%0.6f\n",
	   KBIAS,KLBIAS,KFBIAS,biasWT,biasWToutlier,0.0,HBIAS,pTP,F,theta,lambda);
    if(Poutlier>0.0)
      printf("Poutlier=%0.6f(Bias=%0.6f,Pen=%0.6e),PoutlierEnd=%0.6f\n",
	     Poutlier,OutlierBias,OutlierPenalty,PoutlierEnd);
    fflush(stdout);
  }

  double R = F+pTP*pFN/lambda;

  MissPenalty = MissPenaltyD = -log((theta - resKB)*R);
  if(GAMMA != 1.0)
    MissPenalty = MissPenaltyD *= GAMMA;

  EndThresh = EndScoreThresh * MissPenaltyD;
  EndThresh2 = EndScoreThresh2 * MissPenaltyD;

  LBias = LBiasD = KLBIAS* (R - HBIAS*pTP*pFN/lambda) * MissPenaltyD;
  FBias = FBiasD = KFBIAS* HBIAS * (2.0 * pFN/pTP) * MissPenaltyD;
  
  LBias2 = LBiasD2 = LBiasD*2.0;
  FBiasD_2 = FBiasD*0.5;
  
  EndBias = EndBiasD = FBiasD * 0.5 + MissPenaltyD * 2.0;

  Frate = F + R - thetaScale/theta;// WAS107 R - 1.0/theta;
  FLBias = LBias + Frate * biasWT;

  if(PairSplit && PairsplitRef>=2){/* assume one of the two maps is a in-silico reference */
    var = varD *= 0.5;
    varF = varFD *= 0.5;
    varR = varRD *= 0.5;
    //    MissPenalty = MissPenaltyD *= 0.5;
    Frate *= 0.5;
  }

  //  MatchBias = FBias + KBIAS*0.5 + (BETA ? Beta : 0.0)*(log((2.0*Glambda)/(pTP*pTP)) + 0.5*log(varF+2.0*Glambda*var/(pTP*pTP) + 2.0*Glambda*Glambda*varR/(pTP*pTP*pTP*pTP)));
  //  MatchScore = 0.0;
  MatchBias = MatchBiasD = FBiasD + KBIAS * 0.5 + (BETA ? Beta : 0.0)*(log(Glambda/(theta*pTP*pTP))-log(theta)+ 0.5*log(2.0*M_PI) + 0.5*log(varFD+2.0*Glambda*varD/(pTP*pTP) + 2.0*Glambda*Glambda*varR/(pTP*pTP*pTP*pTP)));
  MatchScore = (BETA ? Beta : 0.0)*(log(2.0*theta)+log(theta)-0.5*log(2*M_PI));

  /* BetaPen adjustment only makes a difference if Poutlier > 0.0 or biastWT < 1*/
  //  double BetaPen = (BETA ? 1.0-Beta : 1.0) * (log(theta) - 0.5*log(2.0*M_PI) + 0.5*log(varF+2.0*Glambda*varD/(pTP*pTP) + 2.0*Glambda*Glambda*varR/(pTP*pTP*pTP*pTP)));
  double BetaPen = (BETA ? 1.0-Beta : 1.0) * (log(theta) - 0.5*log(2.0*M_PI) + (OUTLIER_FIXES ? -1.0 : 1.0) * 0.5*log(varF+2.0*Glambda*varD/(pTP*pTP) + 2.0*Glambda*Glambda*varR/(pTP*pTP*pTP*pTP)));

  MatchBias -= BetaPen * biasWT;
  MatchScore += BetaPen;

  MinScore = -2*(DELTA-1)*MissPenalty;
  MinScore *= 2.0;

  MatchBiasOutlier = MatchBias + OutlierBias;

  varF2 = 2.0f * varF;
  var2 = 2.0f * var;
  varR2 = 2.0f * varR;
  varFD2 = 2.0 * varFD;
  varD2 = 2.0 * varD;
  varRD2 = 2.0 * varRD;

  MatchScore2 = MatchScore + 2.0f * MissPenalty;

  int maxMN = max(maxN,maxM);
  OutlierPenaltyBC = new PFLOAT[maxMN];
  if(!outlierBC)
    memset(OutlierPenaltyBC, 0, maxMN*sizeof(PFLOAT));
  else {
    OutlierPenaltyBC[0] = 0.0;

    int numthreads = 1;
    #ifdef _OPENMP
    numthreads = omp_get_max_threads();
    numthreads = min(numthreads,MaxThreads);
    #endif
    int nthreads = min(maxMN/512, numthreads);

    #pragma omp parallel for num_threads(nthreads) schedule(static, 512) if(nthreads > 1)
    for(int i = 1; i < maxMN; i++)
      OutlierPenaltyBC[i] = -log((double)i);
  }

  if(VERB){
    printf("LBias=%0.6f,FBias=%0.4f,MatchBias=%0.4f,MatchScore=%0.4f,MissPenalty=%0.5f,FLBias=%0.5f,EndBias=%0.5f,sf=%0.4f,sd=%0.4f,sr=%0.6f,Frate=%0.6f,MinScore=%0.3f\n",
	   LBias,FBias,MatchBias,MatchScore,MissPenalty,FLBias,EndBias,SF[0],SD[0],SR[0],Frate,MinScore);
    printf("\t    MatchScore=%0.10e,MissPenalty=%0.10e,MatchScore2=%0.10e,theta=%0.10e,resKB=%0.10e,lambda=%0.10e,F=%0.10e,pFN=%0.10e\n",MatchScore,MissPenalty,MatchScore2,theta,resKB,lambda,F,pFN);
    fflush(stdout);
  }
  
  /* allocate and initialize deltaX[][][] */
  deltaXfreeBlock(Gmap,mapstart,mapend);
  deltaXinitBlock(Gmap,mapstart,mapend);
  /*  for(int i = mapstart; i <= mapend; i++){
    Cmap *Xmap = Gmap[i];
    Xmap->deltaXfree();
    Xmap->deltaXinit();
    }*/

  biasWToutlierF = biasWToutlier;

#if USE_SSE==1
  {
    float tmp;

    tmp = resKB2; v_resKB2=_mm_load_ps1(&tmp);
    tmp = MissPenalty; v_MissPenalty=_mm_load_ps1(&tmp);
    //    tmp = FBias; v_FBias=_mm_load_ps1(&tmp);
    tmp = Frate; v_Frate=_mm_load_ps1(&tmp);
    tmp = var2; v_var2 = _mm_load_ps1(&tmp);
    tmp = varF2; v_varF2 = _mm_load_ps1(&tmp);
    tmp = varR2; v_varR2 = _mm_load_ps1(&tmp);
    tmp = FLBias; v_FLBias=_mm_load_ps1(&tmp);
    tmp = MatchScore2; v_MatchScore2=_mm_load_ps1(&tmp);
    tmp = MatchBiasOutlier; v_MatchBiasOutlier = _mm_load_ps1(&tmp);
    tmp = OutlierPenalty; v_OutlierPenalty = _mm_load_ps1(&tmp);
    tmp = outlierMax; v_outlierMax = _mm_load_ps1(&tmp);
    tmp = -outlierMax; v_outlierMaxneg = _mm_load_ps1(&tmp);
    tmp = biasWToutlier; v_biasWToutlier = _mm_load_ps1(&tmp);
    tmp = MatchBias; v_MatchBias = _mm_load_ps1(&tmp);
    tmp = OutlierBias; v_OutlierBias = _mm_load_ps1(&tmp);

    if(VERB/* HERE >=2 */){
      printf("OutlierPenalty=v_OutlierPenalty=%0.4f\n",OutlierPenalty);
      fflush(stdout);
    }
  }
  
#if USE_AVX==1

  {
    float tmp;
    tmp = resKB2; v256_resKB2=_mm256_broadcast_ss(&tmp);
    tmp = MissPenalty; v256_MissPenalty=_mm256_broadcast_ss(&tmp);
    tmp = FBias; v256_FBias=_mm256_broadcast_ss(&tmp);
    tmp = EndBias; v256_EndBias = _mm256_broadcast_ss(&tmp);
    tmp = Frate; v256_Frate=_mm256_broadcast_ss(&tmp);
    tmp = var2; v256_var2=_mm256_broadcast_ss(&tmp);
    tmp = varF2; v256_varF2=_mm256_broadcast_ss(&tmp);
    tmp = varR2; v256_varR2=_mm256_broadcast_ss(&tmp);
    tmp = FLBias; v256_FLBias=_mm256_broadcast_ss(&tmp);
    tmp = MatchScore2; v256_MatchScore2=_mm256_broadcast_ss(&tmp);
    tmp = MatchBiasOutlier; v256_MatchBiasOutlier=_mm256_broadcast_ss(&tmp);
    tmp = OutlierPenalty; v256_OutlierPenalty = _mm256_broadcast_ss(&tmp);
    tmp = outlierMax; v256_outlierMax = _mm256_broadcast_ss(&tmp);
    tmp = -outlierMax; v256_outlierMaxneg = _mm256_broadcast_ss(&tmp);
    tmp = 0; v256_PZERO = _mm256_broadcast_ss(&tmp);
    v256_biasWToutlier = _mm256_broadcast_ss(&biasWToutlierF);
    tmp = MatchBias; v256_MatchBias = _mm256_broadcast_ss(&tmp);
    tmp = OutlierBias; v256_OutlierBias = _mm256_broadcast_ss(&tmp);
  }

#endif // USE_AVX==1
#endif // USE_SSE==1

#if USE_MIC
  #if USE_PFLOAT==0
    printf("USE_MIC==1 not implemented with USE_PFLOAT==0\n");
    fflush(stdout);exit(1);
  #endif

    {
      float tmp;
      tmp = resKB2; v512_resKB2=_mm512_set1_ps(tmp);
      tmp = MissPenalty; v512_MissPenalty=_mm512_set1_ps(tmp);
      tmp = FBias; v512_FBias=_mm512_set1_ps(tmp);
      tmp = FBiasD_2; v512_FBiasD_2=_mm512_set1_ps(tmp);
      tmp = EndBias; v512_EndBias = _mm512_set1_ps(tmp);
      tmp = Frate; v512_Frate=_mm512_set1_ps(tmp);
      tmp = var; v512_var=_mm512_set1_ps(tmp);
      tmp = var2; v512_var2=_mm512_set1_ps(tmp);
      tmp = varF; v512_varF=_mm512_set1_ps(tmp);
      tmp = varF2; v512_varF2=_mm512_set1_ps(tmp);
      tmp = varR; v512_varR=_mm512_set1_ps(tmp);
      tmp = varR2; v512_varR2=_mm512_set1_ps(tmp);
      tmp = LBias; v512_LBias=_mm512_set1_ps(tmp);
      tmp = LBiasD2; v512_LBiasD2=_mm512_set1_ps(tmp);
      tmp = FLBias; v512_FLBias=_mm512_set1_ps(tmp);
      tmp = MatchScore2; v512_MatchScore2=_mm512_set1_ps(tmp);
      tmp = MatchBiasOutlier; v512_MatchBiasOutlier=_mm512_set1_ps(tmp);
      tmp = OutlierPenalty; v512_OutlierPenalty = _mm512_set1_ps(tmp);
      tmp = outlierMax; v512_outlierMax = _mm512_set1_ps(tmp);
      tmp = -outlierMax; v512_outlierMaxneg = _mm512_set1_ps(tmp);
      tmp = 0; v512_PZERO=_mm512_set1_ps(tmp);
      tmp = biasWToutlier; v512_biasWToutlier = _mm512_set1_ps(tmp);
      tmp = MatchBias; v512_MatchBias = _mm512_set1_ps(tmp);
      tmp = OutlierBias; v512_OutlierBias = _mm512_set1_ps(tmp);
    }
   
    {
      int tmp16[16];
      for(int i= 0; i < 16; i++) 
	tmp16[i] = i;
      v512i_ascending0123 = _mm512_load_epi32(tmp16);
    }

    {
      float tmp16[16];
      for(int i= 0; i < 16; i++) 
	tmp16[i] = i;
      v512_ascending0123 = _mm512_load_ps(tmp16);
    }
#endif // USE_MIC

  if(DEBUG)assert(isfinite(LBias));
  if(DEBUG)assert(isfinite(FBias));
  if(DEBUG)assert(isfinite(MatchBias));
  if(DEBUG)assert(isfinite(MatchScore));
  if(DEBUG)assert(isfinite(MissPenalty));
  if(DEBUG)assert(isfinite(MinScore));
}

static void score_free(Cmap **XXmap,int numX)
{
  deltaXfreeBlock(XXmap,0,numX-1);
  /*  for(int i = 0; i < numX; i++)
      XXmap[i]->deltaXfree();*/
  delete [] OutlierPenaltyBC; OutlierPenaltyBC = NULL;
}

/* Sint() includes Sm(J,I) */
static inline PFLOAT Sint(PFLOAT X,
			  PFLOAT Y,
			  int m,
			  int n)
{
  PFLOAT XYsum = X+Y;
  PFLOAT err = X-Y;
  PFLOAT XYlen = PMAX(PZERO, X-resKB2) + PMAX(PZERO, Y-resKB2);
  PFLOAT err2 = err*err;
  PFLOAT varXY = varF2 + var2 * XYsum;
  if(QUADRATIC_VARIANCE)
    varXY += varR2 * (X*X + Y*Y);
#if USE_MIC && USE_PFLOAT
  PFLOAT Pen = (MatchScore2 - err2*reciprocal_mic(varXY)) - (MissPenalty*(m+n) + XYlen*Frate);      
#else
  PFLOAT Pen = (MatchScore2 - err2/varXY) - (MissPenalty*((PFLOAT)(m+n)) + XYlen*Frate);      
#endif

  PFLOAT Bias = XYlen * FLBias + MatchBias;
  return OutlierBias + (OUTLIER_DELTA(err) ? PMAX(Pen + Bias, OutlierPenalty + Bias * biasWToutlierF) : Pen + Bias);
}

/* Sint() with outlierBC support */
static inline PFLOAT Sint(PFLOAT X,
			  PFLOAT Y,
			  int m,
			  int n,
			  int outlierBC)
{
  PFLOAT XYsum = X+Y;
  PFLOAT err = X-Y;
  PFLOAT XYlen = PMAX(PZERO, X-resKB2) + PMAX(PZERO, Y-resKB2);
  PFLOAT err2 = err*err;
  PFLOAT varXY = varF2 + var2 * XYsum;
  if(QUADRATIC_VARIANCE)
    varXY += varR2 * (X*X + Y*Y);
#if USE_MIC && USE_PFLOAT
  PFLOAT Pen = (MatchScore2 - err2*reciprocal_mic(varXY)) - (MissPenalty*(m+n) + XYlen*Frate);      
#else
  PFLOAT Pen = (MatchScore2 - err2/varXY) - (MissPenalty*((PFLOAT)(m+n)) + XYlen*Frate);      
#endif

  PFLOAT Bias = XYlen*FLBias + MatchBias; 

  if(!OUTLIER_DELTA(err))
    return OutlierBias + Bias + Pen;

  PFLOAT OutPen = OutlierPenalty;
  if(outlierBC)
    OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];

  return OutlierBias + PMAX(Pen + Bias, OutPen + Bias * biasWToutlierF);
}

static inline PFLOAT SintOutlier(PFLOAT X, PFLOAT Y, PFLOAT OutPen)
{
  PFLOAT XYlen = PMAX(PZERO, X-resKB2) + PMAX(PZERO, Y-resKB2);
  PFLOAT Bias = XYlen * FLBias + MatchBias;

  return Bias * biasWToutlierF + OutPen;
  // WAS return (PMAX(PZERO, X-resKB2) + PMAX(PZERO, Y-resKB2)) * FLBias + OutPen;
}

static inline PFLOAT Sint(PFLOAT X,
			  PFLOAT Y,
			  PFLOAT m_n)
  
{
  PFLOAT XYsum = X + Y;
  PFLOAT err = X - Y;
  PFLOAT XYlen = PMAX(PZERO, X - resKB2) + PMAX(PZERO, Y - resKB2);
  PFLOAT err2 = err * err;
  PFLOAT varXY = varF2 + var2 * XYsum;
  if(QUADRATIC_VARIANCE)
    varXY += varR2 * (X*X + Y*Y);
#if USE_MIC && USE_PFLOAT
  PFLOAT Pen = (MatchScore2 - err2*reciprocal_mic(varXY)) - (MissPenalty*(m_n) + XYlen*Frate);      
#else
  PFLOAT Pen = (MatchScore2 - err2/varXY) - (MissPenalty*(m_n) + XYlen*Frate);      
#endif

  PFLOAT Bias = XYlen * FLBias + MatchBias;

  return OutlierBias + (OUTLIER_DELTA(err) ? PMAX(Pen + Bias, OutlierPenalty + Bias * biasWToutlierF) : Pen + Bias);
}

#if USE_SSE && USE_AVX && USE_PFLOAT
static inline __m128 Sint128( __m128 X,
			      __m128 Y,
			      __m128 m,
			      __m128 n)
  
{
   __m128 XYsum = _mm_add_ps(X, Y);
   __m128 err = _mm_sub_ps(X, Y);
   __m128 XYlen = _mm_add_ps(_mm_max_ps(ZERO128, _mm_sub_ps(X,v_resKB2)), _mm_max_ps(ZERO128, _mm_sub_ps(Y,v_resKB2)));
   __m128 err2 = _mm_mul_ps(err, err);
#if QUADRATIC_VARIANCE
   __m128 XYsqsum = _mm_add_ps(_mm_mul_ps(X,X),_mm_mul_ps(Y,Y));
   __m128 varXY = _mm_add_ps(v_varF2, _mm_add_ps(_mm_mul_ps(v_var2, XYsum), _mm_mul_ps(v_varR2, XYsqsum)));
#else
   __m128 varXY = _mm_add_ps(v_varF2, _mm_mul_ps(v_var2, XYsum));
#endif

   __m128 Pen = _mm_sub_ps(_mm_sub_ps(v_MatchScore2, _mm_mul_ps(XYlen, v_Frate)), _mm_add_ps(_mm_div_ps(err2, varXY), _mm_mul_ps(v_MissPenalty, _mm_add_ps(m, n))));     
   __m128 Bias = _mm_add_ps(_mm_mul_ps(XYlen, v_FLBias), v_MatchBias);
   __m128 PenBias = _mm_add_ps(Pen, Bias);
   __m128 PenBiasOutlier = _mm_max_ps(PenBias, _mm_add_ps(v_OutlierPenalty, _mm_mul_ps(Bias, v_biasWToutlier)));
#if OUTLIER_MAX>=1
   PenBiasOutlier = _mm_blendv_ps(PenBias, PenBiasOutlier,
				  _mm_and_ps(_mm_cmp_ps(err, v_outlierMax, _CMP_LT_OQ), _mm_cmp_ps(err, v_outlierMaxneg, _CMP_GT_OQ)) /* OUTLIER_DELTA(err) */);
#endif
  return _mm_add_ps(v_OutlierBias, PenBiasOutlier);
}

static inline __m256 Sint256(__m256 X,
			     __m256 Y,
			     int m,
			     int n)
  
{
  __m256 XYsum = _mm256_add_ps(X, Y);
  __m256 err = _mm256_sub_ps(X, Y);
  __m256 XYlen = _mm256_add_ps(_mm256_max_ps(v256_PZERO, _mm256_sub_ps(X, v256_resKB2)), _mm256_max_ps(v256_PZERO, _mm256_sub_ps(Y, v256_resKB2)));
  __m256 err2 = _mm256_mul_ps(err, err);

#if QUADRATIC_VARIANCE
  __m256 XYsqsum = _mm256_add_ps(_mm256_mul_ps(X, X), _mm256_mul_ps(Y, Y));
  __m256 varXY = _mm256_add_ps(v256_varF2, _mm256_add_ps(_mm256_mul_ps(v256_var2, XYsum), _mm256_mul_ps(v256_varR2, XYsqsum)));
#else
  __m256 varXY = _mm256_add_ps(v256_varF2, _mm256_mul_ps(v256_var2, XYsum));
#endif

  __m256 Pen = _mm256_sub_ps(_mm256_sub_ps(v256_MatchScore2, _mm256_div_ps(err2, varXY)), _mm256_add_ps(_mm256_mul_ps(v256_MissPenalty, _mm256_set1_ps((float)(m+n))), _mm256_mul_ps(XYlen, v256_Frate)));      

  __m256 Bias = _mm256_add_ps(_mm256_mul_ps(XYlen, v256_FLBias), v256_MatchBias);
  __m256 PenBias = _mm256_add_ps(Pen, Bias);
  __m256 PenBiasOutlier = _mm256_max_ps(PenBias, _mm256_add_ps(v256_OutlierPenalty, _mm256_mul_ps(Bias, v256_biasWToutlier)));
#if OUTLIER_MAX>=1
  PenBiasOutlier = _mm256_blendv_ps(PenBias, PenBiasOutlier,
			     _mm256_and_ps(_mm256_cmp_ps(err, v256_outlierMax, _CMP_LT_OQ), _mm256_cmp_ps(err, v256_outlierMaxneg, _CMP_GT_OQ)) /* OUTLIER_DELTA(err) */
			     );
#endif

  return _mm256_add_ps(v256_OutlierBias, PenBiasOutlier);
}

#endif // USE_SSE && USE_AVX && USE_PFLOAT

#if USE_MIC && USE_PFLOAT

static inline void Sint512_prefetch(void)
{
	_mm_prefetch((char*)&v512_resKB2, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_varF2, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_var2, _MM_HINT_T0);
#if QUADRATIC_VARIANCE
	_mm_prefetch((char*)&v512_varR2, _MM_HINT_T0);
#endif
	_mm_prefetch((char*)&v512_MatchScore2, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_MissPenalty, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_Frate, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_FLBias, _MM_HINT_T0);
	//	_mm_prefetch((char*)&v512_MatchBiasOutlier, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_MatchBias, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_OutlierBias, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_OutlierPenalty, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_outlierMax, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_outlierMaxneg, _MM_HINT_T0);
	_mm_prefetch((char*)&v512_biasWToutlier, _MM_HINT_T0);
}

static  inline __m512 Sint512b(__m512 X,
			       __m512 Y,
			       __m512 m_n)
  
{
  __m512 XYsum = _mm512_add_ps(X, Y);
  __m512 err = _mm512_sub_ps(X, Y);
  __m512 XYlen = _mm512_add_ps(_mm512_gmax_ps(v512_PZERO, _mm512_sub_ps(X, v512_resKB2)), _mm512_gmax_ps(v512_PZERO, _mm512_sub_ps(Y, v512_resKB2)));
  __m512 err2 = _mm512_mul_ps(err, err);
#if QUADRATIC_VARIANCE
  __m512 XYsqsum = _mm512_add_ps(_mm512_mul_ps(X,X), _mm512_mul_ps(Y,Y));
  __m512 varXY = _mm512_add_ps(v512_varF2, _mm512_add_ps(_mm512_mul_ps(v512_var2, XYsum), _mm512_mul_ps(v512_varR2, XYsqsum)));
#else
  __m512 varXY = _mm512_add_ps(v512_varF2, _mm512_mul_ps(v512_var2, XYsum));
#endif
  __m512 Pen = _mm512_sub_ps(_mm512_sub_ps(v512_MatchScore2, _mm512_mul_ps(err2, _mm512_rcp23_ps(varXY))), 
			     _mm512_add_ps(_mm512_mul_ps(v512_MissPenalty, m_n), _mm512_mul_ps(XYlen, v512_Frate)));      
  __m512 Bias = _mm512_add_ps(_mm512_mul_ps(XYlen, v512_FLBias), v512_MatchBias);
  __m512 PenBias = _mm512_add_ps(Pen, Bias);
#if OUTLIER_MAX>=1
  __m512 PenBiasOutlier = _mm512_mask_gmax_ps(PenBias,  _mm512_cmp_ps_mask(err, v512_outlierMaxneg, _MM_CMPINT_GT) & _mm512_cmp_ps_mask(err, v512_outlierMax, _MM_CMPINT_LT) /* OUTLIER_DELTA(err) */, 
					      PenBias,  _mm512_add_ps(v512_OutlierPenalty, _mm512_mul_ps(Bias, v512_biasWToutlier)));
#else
  __m512 PenBiasOutlier = _mm512_max_ps(PenBias, _mm512_add_ps(v512_OutlierPenalty, _mm512_mul_ps(Bias, v512_biasWToutlier)));
#endif
  return _mm512_add_ps(v512_OutlierBias, PenBiasOutlier);
}

#endif // USE_MIC && USE_PFLOAT

/* convenience function for dealing with outliers */
static inline void SintDetail(PFLOAT X,
			      PFLOAT Y,
			      int m,
			      int n,
			      int J,
			      int I,
			      PFLOAT &Bias,
			      PFLOAT &Pen /* , int verb=0 */)
{
  PFLOAT XYsum = X+Y;
  PFLOAT err = X-Y;
  PFLOAT XYlen = PMAX(PZERO,X-resKB2) + PMAX(PZERO, Y-resKB2);
  PFLOAT err2 = err*err;
  PFLOAT varXY = varF2 + var2 * XYsum;
  if(QUADRATIC_VARIANCE)
    varXY += varR2 * (X*X + Y*Y);
#if USE_MIC && USE_PFLOAT
  Pen = (MatchScore2 - err2*reciprocal_mic(varXY)) - (MissPenalty*(m+n) + XYlen*Frate);      
#else
  Pen = (MatchScore2 - err2/varXY) - (MissPenalty*((PFLOAT)(m+n)) + XYlen*Frate);      
#endif

  Bias = XYlen * FLBias + MatchBias;

  /*  if(verb){
    printf("\n\t XYsum=%0.8f,err=%0.8f,XYlen=%0.8f,err2=%0.10e,varXY=%0.10e,MatchScore2=%0.8f,err2/varXY=%0.8f,MissPenalty=%0.10e,Frate=%0.10e,m=%d,n=%d:Pen=%0.10e\n\t",
	   XYsum,err,XYlen,err2,varXY,MatchScore2,err2/varXY,MissPenalty,Frate,m,n,Pen);
    fflush(stdout);
    }*/
}

/* external version of SintDetail() */
void SintDetailPW(PFLOAT X,
		  PFLOAT Y,
		  int m,
		  int n,
		  int J,
		  int I,
		  PFLOAT &Bias,
		  PFLOAT &Pen)
{
  SintDetail(X,Y,m,n,J,I,Bias,Pen);
}


/* short for Send((X+Y)/2,m,n) + S(X,Y,1,1) - S((X+Y)/2,(X+Y)/2,1,1) */
static inline PFLOAT Sbnd(PFLOAT X,
			  PFLOAT Y,
			  PFLOAT m_n)
{
  PFLOAT err = X-Y;
  PFLOAT XYsum = X+Y;

#if USE_MIC && USE_PFLOAT
  return (max(PZERO,X-resKB2)+max(PZERO,Y-resKB2))*LBias + EndBias - err*err*reciprocal_mic(QUADRATIC_VARIANCE ? varF2+var2*XYsum+varR2*(X*X+Y*Y) : varF2+var2*XYsum) - MissPenalty*m_n;
#else
  return (max(PZERO,X-resKB2)+max(PZERO,Y-resKB2))*LBias + EndBias - err*err/(QUADRATIC_VARIANCE ? varF2+var2*XYsum+varR2*(X*X+Y*Y) : varF2+var2*XYsum) - MissPenalty*m_n;
#endif
}

static inline FLOAT Sbnd(FLOAT X,
			 FLOAT Y,
			 int m,
			 int n)
{
  return Sbnd((PFLOAT)X,(PFLOAT)Y,(PFLOAT)(m+n));
}


static inline PFLOAT Send(PFLOAT X,
			  PFLOAT m_n)
{
  return max(PZERO, X - resKB2) * LBias2 + EndBias - MissPenalty * m_n;
}

static inline FLOAT Send(FLOAT X,
			 int m,
			 int n)
{
  return Send((PFLOAT)X, (PFLOAT)(m+n));
}

#if USE_MIC

static inline __m512 Send512(__m512 X,
			     __m512 m_n)
{
  return _mm512_add_ps(_mm512_mul_ps(_mm512_gmax_ps(_mm512_setzero_ps(), _mm512_sub_ps(X, v512_resKB2)), v512_LBiasD2),
		       _mm512_sub_ps(v512_FBiasD_2, _mm512_mul_ps(v512_MissPenalty, _mm512_sub_ps(m_n, _mm512_set1_ps(2.0)))));
}

#endif // USE_MIC

static inline PFLOAT Sm(int J, int I)
{
  return PZERO;
}

#endif // NGENTIGPAIRSCORE_H
