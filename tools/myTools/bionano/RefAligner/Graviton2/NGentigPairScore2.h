#include <stdlib.h>
#include <stdio.h>
//#include <math.h>
//#include <immintrin.h>

#include "globals.h"
#include "parameters.h"

#ifndef NGENTIGPAIRSCORE2_H
#define NGENTIGPAIRSCORE2_H

static Ident NGentigPairScore2_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/NGentigPairScore2.h 11381 2020-07-28 23:26:39Z tanantharaman $");

#define COLOR_SHIFT 1 /* 0 : assume color shift (at right end) has a fix variance 
			 1 : assume color shift (at right end) depends on Molecule location of rightmost site (approximate Molecule length) */

#define OUTLIER_DELTA(err) (!OUTLIER_MAX || (((RFLOAT)-outlierMax) < err && err < ((RFLOAT)outlierMax)))
#define OUTLIER(x,y) OUTLIER_DELTA(x-y)

#define SDFIX (PairsplitRef>=2 ? 1 : 0) /* HERE HERE 1 */ /* correctly use -ve sd values */

#define RES_BND 1 /* Apply LBias only to interval size above resKB x 2 */

#define HBIAS 1.0 /* 0.0 is original R based score, 1.0 corresponds more closely to GentigPairScore.h */

#define KBIAS (Kbias * biasWT)
#define KLBIAS (KLbias * biasWT)
#define KFBIAS (KFbias * biasWT)

#define BETA 0.0 /* Beta */ /* scale 0.5*log(varF+var*(X+Y)+varR*(X^2+Y^2)) and log(X+Y) terms and corresponding Bias term */

/* Gentig Pairwise Alignment scoring functions */

/* precompute values to speed up likelihood score computation */
static double MatchScore[MAXCOLOR], Ylambda[MAXCOLOR], Ytheta[MAXCOLOR], Glambda[MAXCOLOR], MinScore[MAXCOLOR];
static double OutlierBias, OutlierEndPenalty, OutlierEndBias, MatchBiasD[MAXCOLOR];
static PFLOAT MatchBias[MAXCOLOR];

static PFLOAT EndThresh[MAXCOLOR], EndThresh2[MAXCOLOR];/* Right End thresholds to avoid calling Send(),Sbnd() */
static FLOAT LBiasD[MAXCOLOR], LBiasD2[MAXCOLOR], FBiasD[MAXCOLOR], FBiasD_2[MAXCOLOR], MissPenaltyD[MAXCOLOR], varD[MAXCOLOR], varFD[MAXCOLOR], varRD[MAXCOLOR], resKB2D[MAXCOLOR], EndBiasD[MAXCOLOR];
static PFLOAT LBias[MAXCOLOR], LBias2[MAXCOLOR], FBias[MAXCOLOR], MissPenalty[MAXCOLOR], Frate[MAXCOLOR], OutlierPenalty, FLBias[MAXCOLOR], var[MAXCOLOR], varF[MAXCOLOR], varR[MAXCOLOR], resKB2[MAXCOLOR], EndBias[MAXCOLOR];

static PFLOAT MatchBiasOutlier[MAXCOLOR];/* MatchBias + OutlierBias */
static FLOAT varFD2[MAXCOLOR],varD2[MAXCOLOR],varRD2[MAXCOLOR];/* double of varFD,varD,varRD */
static PFLOAT varF2[MAXCOLOR],var2[MAXCOLOR],varR2[MAXCOLOR];/* double of varF,var, varR */
static PFLOAT MatchScore2[MAXCOLOR];/* MatchScore + 2*MissPenalty */

static double Ytheta12;/* (Ytheta[0]*Ytheta[0]+Ytheta[1]*Ytheta[1])/(Ytheta[0]+Ytheta[1]) */

static PFLOAT MA_mean2;/* 2.0 * MA_mean */
static PFLOAT FM2var; /* 2 * SF[2] * SF[2] */
static PFLOAT SM2var; /* 2 * fabs(SD[2]) * SD[2] */
static PFLOAT RefPen12;/* log(Xtheta12) - 0.5*log(M_PI) */

static long long maxM[MAXCOLOR],maxN[MAXCOLOR];

static PFLOAT *OutlierPenaltyBC = NULL;

#if USE_SSE && USE_PFLOAT
static __m128 

v_resKB2[MAXCOLOR], v_MissPenalty[MAXCOLOR], v_Frate[MAXCOLOR], v_FLBias[MAXCOLOR], v_var2[MAXCOLOR], v_varF2[MAXCOLOR], v_varR2[MAXCOLOR], 
  v_MatchScore2[MAXCOLOR], v_MatchBiasOutlier[MAXCOLOR], v_OutlierPenalty, v_outlierMax, v_outlierMaxneg; // v_FBias[MAXCOLOR];
#endif

#if USE_AVX && USE_PFLOAT
static __m256 

v256_resKB2[MAXCOLOR], v256_MissPenalty[MAXCOLOR], v256_FBias[MAXCOLOR], v256_EndBias[MAXCOLOR], v256_Frate[MAXCOLOR], v256_FLBias[MAXCOLOR], v256_var2[MAXCOLOR], v256_varF2[MAXCOLOR], v256_varR2[MAXCOLOR], 
  v256_MatchScore2[MAXCOLOR], v256_MatchBiasOutlier[MAXCOLOR], v256_OutlierPenalty, v256_outlierMax, v256_outlierMaxneg, v256_PZERO;
		
#endif

		
#if USE_MIC && USE_PFLOAT
static __m512 	

v512_resKB2[MAXCOLOR], v512_MissPenalty[MAXCOLOR], v512_FBias[MAXCOLOR], v512_EndBias[MAXCOLOR], v512_FBiasD_2[MAXCOLOR], v512_Frate[MAXCOLOR], v512_LBias[MAXCOLOR], v512_LBiasD2[MAXCOLOR], v512_FLBias[MAXCOLOR],
  v512_var[MAXCOLOR], v512_var2[MAXCOLOR], v512_varF[MAXCOLOR], v512_varF2[MAXCOLOR], v512_varR[MAXCOLOR], v512_varR2[MAXCOLOR], 
  v512_MatchScore2[MAXCOLOR], v512_MatchBiasOutlier[MAXCOLOR], v512_OutlierPenalty, v512_outlierMax, v512_outlierMaxneg, v512_PZERO;


static __m512 v512_ascending0123;
static __m512i v512i_ascending0123;

#endif // USE_MIC && USE_PFLOAT

static void score_init(int mapstart, int mapend)
{
  if(mapend < mapstart)
    return;

  if(RES_BND != 1){
    printf("RES_BND=%d no longer supported\n",RES_BND);
    exit(1);
  }

  if(Beta != 0.0){
    printf("Beta=%0.3f: -Beta != 0.0 not supported\n",Beta);
    exit(1);
  }
  if(biasWT > 0.0 && biasWToutlier < 1.0){
    printf("-biaswtOutlier %0.3f not yet implemented for 2 colors\n", biasWToutlier);
    exit(1);
  }

  if(DEBUG) assert(colors==2);

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
  if(VERB)
    printf("Poutlier=%0.6f(Bias=%0.6f,Pen=%0.6e),PoutlierEnd=%0.6f\n",
	   Poutlier,OutlierBias,OutlierPenalty,PoutlierEnd);

  for(int c = 0; c < colors; c++){
    double pFN = FN[c];
    double pTP = 1.0 - pFN;
    if(DEBUG) assert(pFN < 1.0);
    if(DEBUG) assert(pTP > 0.0);
    double F = FP[c]*0.01;
    var[c] = varD[c] = (SDFIX ? fabs(SD[c]) : SD[c]) * SD[c];
    varF[c] = varFD[c] = SF[c]*SF[c]*2.0;
    varR[c] = varRD[c] = SR[c]*SR[c];

    /* compute theta, the average observed interval size */
    double Gtheta = 0.0;
    double theta = 0.0;
    int sitecnt = 0;
    maxM[c] = 0;
    for(int i = mapstart; i <= mapend;i++){
      Cmap *Xmap = Gmap[i];
      int M = Xmap->numsite[c];
      sitecnt += M;
      theta += Xmap->site[c][M+1];
      if(M > maxM[c])
	maxM[c] = M;
    }
    // HERE HERE : continue
    if(XmapCount[c] > 0){/* use command line values for XmapCount,XmapSites,XmapLength,Gtheta */
      //    int nummaps = XmapCount[c];
      sitecnt = XmapSites[c];
      theta = XmapLength[c];
      Gtheta = XmapGtheta[c];
    } else {
      /* compute Gtheta, the geometric mean of the observed interval size */
      for(int i= mapstart; i <= mapend;i++){
	Cmap *Xmap = Gmap[i];
	int M = Xmap->numsite[c];
	FLOAT *X = Xmap->site[c];
	// END fragments are ignored since they are often truncated to 0 
	//    Gtheta += log(X[1]+X[M+1]-X[M]);/* treat end frags as one sample */
	for(int i= 1; i < M;i++)
	  Gtheta += log(X[i+1]-X[i]);
      }
      Gtheta /= (sitecnt - (mapend-mapstart+1));
      Gtheta = exp(Gtheta);
    }
  
    maxM[c] = (maxM[c] + 0xf) & ~0xf;
   
    maxN[c] = maxM[c];
    if(DEBUG) assert(maxN[c] <= MAXINT);

    if(sitecnt <= 0 || theta <= 0.0){
      fprintf(stderr,"score_init():Insufficient data for c=%d:mapstart=%d,mapend=%d,sitecnt=%d\n",
	      c,mapstart,mapend,sitecnt);
      exit(1);
    }

    theta /= sitecnt;

    double resKB = (RES_BND ? mres * PixelLen : 0.0);
    resKB2[c] = resKB2D[c] = 2.0 * resKB;

    double lambda = pTP/(1.0/(theta - resKB - varD[c]*0.5) - F);/* HERE : does not include correction for SR[] */
    if(VERB>=2){
      printf("c=%d:pTP=%0.17f,theta=%0.17f,resKB=%0.17f,varD=%0.17f(sd=%0.10f),F=%0.17f:lambda=%0.10f\n",
	     c,pTP,theta,resKB,varD[c],SD[c],F,lambda);
      fflush(stdout);
    }
    if(lambda <= 0.0){
      printf("ERROR:c=%d,theta=%0.6f,pTP=%0.5f,resKB=%0.6f,varD=%0.6f,F=%0.6f:lambda=%0.6f!\n",
	     c,theta,pTP,resKB,varD[c],F,lambda);
      printf("      mres=%0.6f,PixelLen=%0.6f\n",mres,PixelLen);
      if(1.0/theta < F)
	printf("      False Positive rate exceeds average site density for c=%d\n",c);
      printf("Error parameters inconsistent with average site density\n");
      exit(1);
    }
    Glambda[c] = (Gtheta/theta)*lambda;

    Ylambda[c] = lambda;
    Ytheta[c] = theta;

    if(F*theta >= 1.0){
      printf("c=%d:FP density=%0.4f, theta=%0.3f: FP is too large for meaningful alignment scoring\n",c,F,theta);
      exit(1);
    }

    if(VERB){
      printf("c=%d:Kbias=%0.2f,KLbias=%0.2f,KFbias=%0.2f(biasWT=%0.2f),Beta=%0.2f,HBIAS=%0.2f,pTP=%0.6f,F=%0.6f,theta=%0.6f,lambda=%0.6f\n",
	     c,KBIAS,KLBIAS,KFBIAS,biasWT,0.0,HBIAS,pTP,F,theta,lambda);
      fflush(stdout);
    }

    double R = F + pTP*pFN/lambda;

    MissPenalty[c] = MissPenaltyD[c] = -log((theta - resKB)*R);
    if(GAMMA != 1.0)
      MissPenalty[c] = MissPenaltyD[c] *= GAMMA;

    EndThresh[c] = EndScoreThresh * MissPenaltyD[c];
    EndThresh2[c] = EndScoreThresh2 * MissPenaltyD[c];

    LBias[c] = LBiasD[c] = KLBIAS* (R - HBIAS*pTP*pFN/lambda) * MissPenaltyD[c];
    FBias[c] = FBiasD[c] = KFBIAS* HBIAS * (2.0 * pFN/pTP) * MissPenaltyD[c];
  
    LBias2[c] = LBiasD2[c] = LBiasD[c] * 2.0;
    FBiasD_2[c] = FBiasD[c] * 0.5;
  
    EndBias[c] = EndBiasD[c] = FBiasD[c] * 0.5 + MissPenaltyD[c] * 2.0;

    Frate[c] = F+ R - thetaScale/theta;
    FLBias[c] = LBias[c] + Frate[c] * biasWT;

    if(PairSplit && PairsplitRef>=2){/* assume one of the two maps is a in-silico reference */
      var[c] = varD[c] *= 0.5;
      varF[c] = varFD[c] *= 0.5;
      varR[c] = varRD[c] *= 0.5;
      //    MissPenalty[c] = MissPenaltyD[c] *= 0.5;
      Frate[c] *= 0.5;
    }

    //  MatchBias = FBias + KBIAS*0.5 + (BETA ? Beta : 0.0)*(log((2.0*Glambda)/(pTP*pTP)) + 0.5*log(varF+2.0*Glambda*var/(pTP*pTP) + 2.0*Glambda*Glambda*varR/(pTP*pTP*pTP*pTP)));
    //  MatchScore = 0.0;
    MatchBias[c] = MatchBiasD[c] = FBiasD[c] + KBIAS * 0.5 + (BETA ? Beta : 0.0)*(log(Glambda[c]/(theta*pTP*pTP))-log(theta)+ 0.5*log(2.0*M_PI) + 0.5*log(varFD[c]+2.0*Glambda[c]*varD[c]/(pTP*pTP) + 
																			  2.0*Glambda[c]*Glambda[c]*varR[c]/(pTP*pTP*pTP*pTP)));
    MatchScore[c] = (BETA ? Beta : 0.0)*(log(2.0*theta)+log(theta)-0.5*log(2*M_PI));

    /* NOTE : BetaPen adjustment only makes a difference if Poutlier > 0.0 or biastWT < 1*/
    //  double BetaPen = (BETA ? 1.0-Beta : 1.0) * (log(theta) - 0.5*log(2.0*M_PI) + 0.5*log(varF+2.0*Glambda*varD/(pTP*pTP) + 2.0*Glambda*Glambda*varR/(pTP*pTP*pTP*pTP)));
    double BetaPen = (BETA ? 1.0-Beta : 1.0) * (log(theta) - 0.5*log(2.0*M_PI) + (OUTLIER_FIXES ? -1.0 : 1.0) * 0.5*log(varF[c]+2.0*Glambda[c]*varD[c]/(pTP*pTP) + 2.0*Glambda[c]*Glambda[c]*varR[c]/(pTP*pTP*pTP*pTP)));

    MatchBias[c] -= BetaPen * biasWT;
    MatchScore[c] += BetaPen;

    MinScore[c] = -2*(DELTA-1)*MissPenalty[c];
    MinScore[c] *= 2.0;

    MatchBiasOutlier[c] = MatchBias[c] + OutlierBias;

    varF2[c] = 2.0f * varF[c];
    var2[c] = 2.0f * var[c];
    varR2[c] = 2.0f * varR[c];
    varFD2[c] = 2.0 * varFD[c];
    varD2[c] = 2.0 * varD[c];
    varRD2[c] = 2.0 * varRD[c];

    MatchScore2[c] = MatchScore[c] + 2.0f * MissPenalty[c];

    if(VERB){
      printf("c=%d:LBias=%0.6f,FBias=%0.4f,MatchBias=%0.4f,MatchScore=%0.4f,MissPenalty=%0.5f,FLBias=%0.5f,sf=%0.4f,sd=%0.4f,sr=%0.6f,Frate=%0.6f,MinScore=%0.3f\n",
	     c,LBias[c],FBias[c],MatchBias[c],MatchScore[c],MissPenalty[c],FLBias[c],SF[c],SD[c],SR[c],Frate[c],MinScore[c]);
      fflush(stdout);
    }

    if(DEBUG)assert(isfinite(LBias[c]));
    if(DEBUG)assert(isfinite(FBias[c]));
    if(DEBUG)assert(isfinite(MatchBias[c]));
    if(DEBUG)assert(isfinite(MatchScore[c]));
    if(DEBUG)assert(isfinite(MissPenalty[c]));
    if(DEBUG)assert(isfinite(MinScore[c]));

  } // c = 0..colors-1

  Ytheta12 = (Ytheta[0]*Ytheta[0] + Ytheta[1]*Ytheta[1])/(Ytheta[0]+Ytheta[1]);/* average interval distance between sites of opposite colors in XX */
  RefPen12 = log(Ytheta12) - 0.5*log(M_PI);
  MA_mean2 = MA_mean * 2.0;
  FM2var = 4.0 * SF[2] * SF[2];
  SM2var = 2.0 * fabs(SD[2]) * SD[2];
  if(VERB){
    printf("    Ytheta12=%0.6f,RefPen12=%0.8f\n",Ytheta12,RefPen12);
    fflush(stdout);
  }

  int maxMN = max(maxN[0],maxM[0]);
  for(int c = 1; c < colors; c++){
    maxMN = max(maxMN,maxN[c]);
    maxMN = max(maxMN,maxM[c]);
  }

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

  /* allocate and initialize deltaX[][][] */
  deltaXfreeBlock(Gmap,mapstart,mapend);
  deltaXinitBlock(Gmap,mapstart,mapend);

#if USE_SSE==1
  #if USE_PFLOAT==0
  printf("USE_SSE==1 or USE_AVX==1 not implemented with USE_PFLOAT==0\n");
  exit(1);
  #endif

  for(int c = 0; c < colors; c++){
    float tmp;

    tmp = resKB2[c]; v_resKB2[c] = _mm_load_ps1(&tmp);
    tmp = MissPenalty[c]; v_MissPenalty[c] = _mm_load_ps1(&tmp);
    //    tmp = FBias[c]; v_FBias[c] = _mm_load_ps1(&tmp);
    tmp = Frate[c]; v_Frate[c] = _mm_load_ps1(&tmp);
    tmp = var2[c]; v_var2[c] = _mm_load_ps1(&tmp);
    tmp = varF2[c]; v_varF2[c] = _mm_load_ps1(&tmp);
    tmp = varR2[c]; v_varR2[c] = _mm_load_ps1(&tmp);
    tmp = FLBias[c]; v_FLBias[c] = _mm_load_ps1(&tmp);
    tmp = MatchScore2[c]; v_MatchScore2[c] = _mm_load_ps1(&tmp);
    tmp = MatchBiasOutlier[c]; v_MatchBiasOutlier[c] = _mm_load_ps1(&tmp);
  }
  {
    float tmp = OutlierPenalty; v_OutlierPenalty = _mm_load_ps1(&tmp);
    tmp = outlierMax; v_outlierMax = _mm_load_ps1(&tmp);
    tmp = -outlierMax; v_outlierMaxneg = _mm_load_ps1(&tmp);

    if(VERB>=2){
      printf("OutlierPenalty = v_OutlierPenalty=%0.4f\n",OutlierPenalty);
      fflush(stdout);
    }
  } 
  
#if USE_AVX==1

  for(int c = 0; c < colors; c++){
    float tmp;

    tmp = resKB2[c]; v256_resKB2[c] = _mm256_broadcast_ss(&tmp);
    tmp = MissPenalty[c]; v256_MissPenalty[c] = _mm256_broadcast_ss(&tmp);
    tmp = FBias[c]; v256_FBias[c] = _mm256_broadcast_ss(&tmp);
    tmp = EndBias[c]; v256_EndBias[c] = _mm256_broadcast_ss(&tmp);
    tmp = Frate[c]; v256_Frate[c] = _mm256_broadcast_ss(&tmp);
    tmp = var2[c]; v256_var2[c] = _mm256_broadcast_ss(&tmp);
    tmp = varF2[c]; v256_varF2[c] = _mm256_broadcast_ss(&tmp);
    tmp = varR2[c]; v256_varR2[c] = _mm256_broadcast_ss(&tmp);
    tmp = FLBias[c]; v256_FLBias[c] = _mm256_broadcast_ss(&tmp);
    tmp = MatchScore2[c]; v256_MatchScore2[c] = _mm256_broadcast_ss(&tmp);
    tmp = MatchBiasOutlier[c]; v256_MatchBiasOutlier[c] = _mm256_broadcast_ss(&tmp);
  }
  {
    float tmp;
    tmp = OutlierPenalty; v256_OutlierPenalty = _mm256_broadcast_ss(&tmp);
    tmp = outlierMax; v256_outlierMax = _mm256_broadcast_ss(&tmp);
    tmp = -outlierMax; v256_outlierMaxneg = _mm256_broadcast_ss(&tmp);
    tmp = 0; v256_PZERO=_mm256_broadcast_ss(&tmp);
  }

#endif // USE_AVX==1
#endif // USE_SSE==1

#if USE_MIC
  if(USE_PFLOAT==0){
    printf("USE_MIC==1 not implemented with USE_PFLOAT==0\n");
    exit(1);
  }

  for(int c = 0; c < colors; c++){
    float tmp;
    tmp = resKB2[c]; v512_resKB2[c] = _mm512_set1_ps(tmp);
    tmp = MissPenalty[c]; v512_MissPenalty[c] = _mm512_set1_ps(tmp);
    tmp = FBias[c]; v512_FBias[c] = _mm512_set1_ps(tmp);
    tmp = FBiasD_2[c]; v512_FBiasD_2[c] = _mm512_set1_ps(tmp);
    tmp = EndBias[c]; v512_EndBias[c] = _mm512_set1_ps(tmp);
    tmp = Frate[c]; v512_Frate[c] = _mm512_set1_ps(tmp);
    tmp = var[c]; v512_var[c] = _mm512_set1_ps(tmp);
    tmp = var2[c]; v512_var2[c] = _mm512_set1_ps(tmp);
    tmp = varF[c]; v512_varF[c] = _mm512_set1_ps(tmp);
    tmp = varF2[c]; v512_varF2[c] = _mm512_set1_ps(tmp);
    tmp = varR[c]; v512_varR[c] = _mm512_set1_ps(tmp);
    tmp = varR2[c]; v512_varR2[c] = _mm512_set1_ps(tmp);
    tmp = LBias[c]; v512_LBias[c] = _mm512_set1_ps(tmp);
    tmp = LBiasD2[c]; v512_LBiasD2[c] = _mm512_set1_ps(tmp);
    tmp = FLBias[c]; v512_FLBias[c] = _mm512_set1_ps(tmp);
    tmp = MatchScore2[c]; v512_MatchScore2[c] = _mm512_set1_ps(tmp);
    tmp = MatchBiasOutlier[c]; v512_MatchBiasOutlier[c] = _mm512_set1_ps(tmp);
  }
  {
    float tmp;
    tmp = OutlierPenalty; v512_OutlierPenalty = _mm512_set1_ps(tmp);
    tmp = outlierMax; v512_outlierMax = _mm512_set1_ps(tmp);
    tmp = -outlierMax; v512_outlierMaxneg = _mm512_set1_ps(tmp);
    tmp = 0; v512_PZERO = _mm512_set1_ps(tmp);
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

}

static void score_free(Cmap **XXmap,int numX)
{
  deltaXfreeBlock(XXmap,0,numX-1);
  delete [] OutlierPenaltyBC; OutlierPenaltyBC = NULL;
}

/* Sint() includes Sm(J,I) */
static inline PFLOAT Sint(PFLOAT X,
			  PFLOAT Y,
			  int m,
			  int n,
			  int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  PFLOAT XYsum = X + Y;
  PFLOAT err = X - Y;
  PFLOAT XYlen = PMAX(PZERO, X-resKB2[c]) + PMAX(PZERO, Y-resKB2[c]);
  PFLOAT err2 = err * err;
  PFLOAT varXY = varF2[c] + var2[c] * XYsum;
  if(QUADRATIC_VARIANCE)
    varXY += varR2[c] * (X*X + Y*Y);

#if USE_MIC && USE_PFLOAT
  PFLOAT Pen = (MatchScore2[c] - err2*reciprocal_mic(varXY)) - (MissPenalty[c]*(m+n) + XYlen*Frate[c]);      
#else
  PFLOAT Pen = (MatchScore2[c] - err2/varXY) - (MissPenalty[c]*((PFLOAT)(m+n)) + XYlen*Frate[c]);      
#endif

  return XYlen*FLBias[c] + MatchBiasOutlier[c] + (OUTLIER_DELTA(err) ? PMAX(Pen,OutlierPenalty) : Pen);
}

/* Sint() with outlierBC support */
static inline PFLOAT Sint(PFLOAT X,
			  PFLOAT Y,
			  int m,
			  int n,
			  int outlierBC,
			  int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  PFLOAT XYsum = X + Y;
  PFLOAT err = X - Y;
  PFLOAT XYlen = PMAX(PZERO, X - resKB2[c]) + PMAX(PZERO, Y - resKB2[c]);
  PFLOAT err2 = err * err;
  PFLOAT varXY = varF2[c] + var2[c] * XYsum;
  if(QUADRATIC_VARIANCE)
    varXY += varR2[c] * (X*X + Y*Y);
#if USE_MIC && USE_PFLOAT
  PFLOAT Pen = (MatchScore2[c] - err2*reciprocal_mic(varXY)) - (MissPenalty[c]*(m+n) + XYlen*Frate[c]);      
#else
  PFLOAT Pen = (MatchScore2[c] - err2/varXY) - (MissPenalty[c]*((PFLOAT)(m+n)) + XYlen*Frate[c]);      
#endif

  PFLOAT Bias = XYlen*FLBias[c] + MatchBiasOutlier[c]; 

  if(!OUTLIER_DELTA(err))
    return Bias + Pen;

  PFLOAT OutPen = OutlierPenalty;
  if(outlierBC)
    OutPen += OutlierPenaltyBC[m] + OutlierPenaltyBC[n];

  return Bias + PMAX(Pen,OutPen);
}

static inline PFLOAT SintOutlier(PFLOAT X, PFLOAT Y, PFLOAT OutPen, int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  return (PMAX(PZERO, X - resKB2[c]) + PMAX(PZERO, Y - resKB2[c])) * FLBias[c] + OutPen;
}

static inline PFLOAT Sint(PFLOAT X,
			  PFLOAT Y,
			  PFLOAT m_n,
			  int c)
  
{
  if(DEBUG>=2) assert(c==0 || c==1);
  PFLOAT XYsum = X + Y;
  PFLOAT err = X - Y;
  PFLOAT XYlen = PMAX(PZERO, X - resKB2[c]) + PMAX(PZERO, Y - resKB2[c]);
  PFLOAT err2 = err * err;
  PFLOAT varXY = varF2[c] + var2[c] * XYsum;
  if(QUADRATIC_VARIANCE)
    varXY += varR2[c] * (X*X + Y*Y);
#if USE_MIC && USE_PFLOAT
  PFLOAT Pen = (MatchScore2[c] - err2*reciprocal_mic(varXY)) - (MissPenalty[c]*(m_n) + XYlen*Frate[c]);      
#else
  PFLOAT Pen = (MatchScore2[c] - err2/varXY) - (MissPenalty[c]*(m_n) + XYlen*Frate[c]);      
#endif

  return XYlen * FLBias[c] + MatchBiasOutlier[c] + (OUTLIER_DELTA(err) ? PMAX(Pen,OutlierPenalty) : Pen);
}

#if USE_SSE && USE_AVX && USE_PFLOAT
static inline __m128 Sint128( __m128 X,
			      __m128 Y,
			      __m128 m,
			      __m128 n,
			      int c)
  
{
  if(DEBUG>=2) assert(c==0 || c==1);
   __m128 XYsum = _mm_add_ps(X, Y);
   __m128 XYlen = _mm_add_ps(_mm_max_ps(ZERO128, _mm_sub_ps(X, v_resKB2[c])), _mm_max_ps(ZERO128, _mm_sub_ps(Y,v_resKB2[c])));
   __m128 err = _mm_sub_ps(X, Y);
#if QUADRATIC_VARIANCE
   __m128 XYsqsum = _mm_add_ps(_mm_mul_ps(X,X),_mm_mul_ps(Y,Y));
   __m128 varXY = _mm_add_ps(v_varF2[c], _mm_add_ps(_mm_mul_ps(v_var2[c], XYsum), _mm_mul_ps(v_varR2[c], XYsqsum)));
#else
   __m128 varXY = _mm_add_ps(v_varF2[c], _mm_mul_ps(v_var2[c], XYsum));
#endif
   
   __m128 v_a = _mm_mul_ps(XYlen, v_Frate[c]);
   __m128 v_err_sq = _mm_mul_ps(err, err);
   
   __m128 Pen = _mm_sub_ps(_mm_sub_ps(v_MatchScore2[c], v_a), _mm_add_ps(_mm_div_ps(v_err_sq, varXY), _mm_mul_ps(v_MissPenalty[c], _mm_add_ps(m, n))));     
   
   __m128 PenBias = _mm_max_ps(Pen, v_OutlierPenalty);
#if OUTLIER_MAX>=1
   PenBias = _mm_blendv_ps(Pen, PenBias,
			      _mm_and_ps(_mm_cmp_ps(err, v_outlierMax, _CMP_LT_OQ), _mm_cmp_ps(err, v_outlierMaxneg, _CMP_GT_OQ)) /* OUTLIER_DELTA(err) */
			      );
#endif
  return _mm_add_ps(_mm_mul_ps(XYlen, v_FLBias[c]), _mm_add_ps(v_MatchBiasOutlier[c], PenBias));
}

static inline __m256 Sint256(__m256 X,
			     __m256 Y,
			     int m,
			     int n,
			     int c)
  
{
  if(DEBUG>=2) assert(c==0 || c==1);
  __m256 XYsum = _mm256_add_ps(X, Y);
  __m256 XYlen = _mm256_add_ps(_mm256_max_ps(v256_PZERO, _mm256_sub_ps(X, v256_resKB2[c])), _mm256_max_ps(v256_PZERO, _mm256_sub_ps(Y, v256_resKB2[c])));
  __m256 err = _mm256_sub_ps(X, Y);
#if QUADRATIC_VARIANCE
  __m256 XYsqsum = _mm256_add_ps(_mm256_mul_ps(X, X), _mm256_mul_ps(Y, Y));
  __m256 varXY = _mm256_add_ps(v256_varF2[c], _mm256_add_ps(_mm256_mul_ps(v256_var2[c], XYsum), _mm256_mul_ps(v256_varR2[c], XYsqsum)));
#else
  __m256 varXY = _mm256_add_ps(v256_varF2[c], _mm256_mul_ps(v256_var2[c], XYsum));
#endif

  __m256 Pen = _mm256_sub_ps(_mm256_sub_ps(v256_MatchScore2[c], _mm256_div_ps(_mm256_mul_ps(err, err), varXY)), _mm256_add_ps(_mm256_mul_ps(v256_MissPenalty[c], _mm256_set1_ps((float)(m+n))), _mm256_mul_ps(XYlen, v256_Frate[c])));      

  __m256 PenBias = _mm256_max_ps(Pen, v256_OutlierPenalty);
#if OUTLIER_MAX>=1
  PenBias = _mm256_blendv_ps(Pen, PenBias,
			     _mm256_and_ps(_mm256_cmp_ps(err, v256_outlierMax, _CMP_LT_OQ), _mm256_cmp_ps(err, v256_outlierMaxneg, _CMP_GT_OQ)) /* OUTLIER_DELTA(err) */
			     );
#endif

  return _mm256_add_ps(_mm256_mul_ps(XYlen, v256_FLBias[c]), _mm256_add_ps(v256_MatchBiasOutlier[c], PenBias));
}

#endif // USE_SSE && USE_AVX && USE_PFLOAT

#if USE_MIC && USE_PFLOAT

static inline void Sint512_prefetch(int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  _mm_prefetch((char*)&v512_resKB2[c], _MM_HINT_T0);
  _mm_prefetch((char*)&v512_varF2[c], _MM_HINT_T0);
  _mm_prefetch((char*)&v512_var2[c], _MM_HINT_T0);
#if QUADRATIC_VARIANCE
  _mm_prefetch((char*)&v512_varR2[c], _MM_HINT_T0);
#endif
  _mm_prefetch((char*)&v512_MatchScore2[c], _MM_HINT_T0);
  _mm_prefetch((char*)&v512_MissPenalty[c], _MM_HINT_T0);
  _mm_prefetch((char*)&v512_Frate[c], _MM_HINT_T0);
  _mm_prefetch((char*)&v512_FLBias[c], _MM_HINT_T0);
  _mm_prefetch((char*)&v512_MatchBiasOutlier[c], _MM_HINT_T0);
  _mm_prefetch((char*)&v512_OutlierPenalty, _MM_HINT_T0);
}

static inline __m512 Sint512(__m512 X,
			     __m512 Y,
			     int m_n,
			     int c)
  
{
  if(DEBUG>=2) assert(c==0 || c==1);
  __m512 XYsum = _mm512_add_ps(X, Y);
//  __m512 XYlen = _mm512_add_ps(_mm512_gmax_ps(v512_PZERO, _mm512_sub_ps(X, v512_resKB2[c])), _mm512_gmax_ps(v512_PZERO, _mm512_sub_ps(Y, v512_resKB2[c])));
  __m512 XYlen = _mm512_add_ps(
	_mm512_gmax_ps(
		_mm512_setzero_ps(), 
		_mm512_sub_ps(X, v512_resKB2[c])
		), 
	_mm512_gmax_ps(
		_mm512_setzero_ps(), 
		_mm512_sub_ps(Y, v512_resKB2[c])
	));
  __m512 err = _mm512_sub_ps(X, Y);
#if QUADRATIC_VARIANCE
  __m512 XYsqsum = _mm512_add_ps(_mm512_mul_ps(X, X), _mm512_mul_ps(Y, Y));
  __m512 varXY = _mm512_add_ps(v512_varF2[c], _mm512_add_ps(_mm512_mul_ps(v512_var2[c], XYsum), _mm512_mul_ps(v512_varR2[c], XYsqsum)));
#else
  __m512 varXY = _mm512_add_ps(v512_varF2[c], _mm512_mul_ps(v512_var2[c], XYsum));
#endif
  __m512 Pen = _mm512_sub_ps(
	_mm512_sub_ps(v512_MatchScore2[c], 
			_mm512_mul_ps(
				_mm512_mul_ps(err, err), 
				_mm512_rcp23_ps(varXY))
		     ), 
	_mm512_add_ps(
		_mm512_mul_ps(v512_MissPenalty[c], 
			      _mm512_set1_ps(m_n)), 
		_mm512_mul_ps(XYlen, v512_Frate[c])));      

#if OUTLIER_MAX>=1
  __m512 PenBias = _mm512_mask_gmax_ps(Pen, 
				       _mm512_cmp_ps_mask(err, v512_outlierMaxneg, _MM_CMPINT_GT) & _mm512_cmp_ps_mask(err, v512_outlierMax, _MM_CMPINT_LT) /* OUTLIER_DELTA(err) */, 
				       Pen,  v512_OutlierPenalty);
#else
  __m512 PenBias = _mm512_gmax_ps(Pen, v512_OutlierPenalty);
#endif

  return _mm512_add_ps(
		       _mm512_add_ps(v512_MatchBiasOutlier[c], 
				     _mm512_mul_ps(XYlen, v512_FLBias[c])),
		       PenBias);
}

static  inline __m512 Sint512b(__m512 X,
			       __m512 Y,
			       __m512 m_n,
			       int c)
  
{
  if(DEBUG>=2) assert(c==0 || c==1);
  __m512 XYsum = _mm512_add_ps(X, Y);
  __m512 XYlen = _mm512_add_ps(
			       _mm512_gmax_ps(v512_PZERO, 
					      _mm512_sub_ps(X, v512_resKB2[c])
					      ), 
			       _mm512_gmax_ps(v512_PZERO, 
					      _mm512_sub_ps(Y, v512_resKB2[c])));
  __m512 err = _mm512_sub_ps(X, Y);
#if QUADRATIC_VARIANCE
  __m512 XYsqsum = _mm512_add_ps(_mm512_mul_ps(X,X), _mm512_mul_ps(Y,Y));
  __m512 varXY = _mm512_add_ps(v512_varF2[c], _mm512_add_ps(_mm512_mul_ps(v512_var2[c], XYsum), _mm512_mul_ps(v512_varR2[c], XYsqsum)));
#else
  __m512 varXY = _mm512_add_ps(v512_varF2[c], _mm512_mul_ps(v512_var2[c], XYsum));
#endif
  __m512 Pen = _mm512_sub_ps(
			     _mm512_sub_ps(v512_MatchScore2[c], 
// 			_mm512_div_ps(
// 				_mm512_mul_ps(err, err), 
					   _mm512_mul_ps(
							 _mm512_mul_ps(err, err), 
							 _mm512_rcp23_ps(varXY))), 
			     _mm512_add_ps(
					   _mm512_mul_ps(v512_MissPenalty[c], m_n), 
					   _mm512_mul_ps(XYlen, v512_Frate[c])));      

#if OUTLIER_MAX>=1
  __m512 PenBias = _mm512_mask_gmax_ps(Pen, _mm512_cmp_ps_mask(err, v512_outlierMaxneg, _MM_CMPINT_GT) & _mm512_cmp_ps_mask(err, v512_outlierMax, _MM_CMPINT_LT), Pen,  v512_OutlierPenalty);
#else
  __m512 PenBias = _mm512_gmax_ps(Pen, v512_OutlierPenalty);
#endif

  return _mm512_add_ps(
		_mm512_add_ps(v512_MatchBiasOutlier[c], 
			      _mm512_mul_ps(XYlen, v512_FLBias[c])), 
		PenBias);
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
			      PFLOAT &Pen,
			      int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  PFLOAT XYsum = X + Y;
  PFLOAT err = X - Y;
  PFLOAT XYlen = PMAX(PZERO, X - resKB2[c]) + PMAX(PZERO, Y - resKB2[c]);
  PFLOAT err2 = err * err;
  PFLOAT varXY = varF2[c] + var2[c] * XYsum;
  if(QUADRATIC_VARIANCE)
    varXY += varR2[c] * (X*X + Y*Y);

#if USE_MIC && USE_PFLOAT
  Pen = (MatchScore2[c] - err2*reciprocal_mic(varXY)) - (MissPenalty[c]*(m+n) + XYlen*Frate[c]);      
#else
  Pen = (MatchScore2[c] - err2/varXY) - (MissPenalty[c]*((PFLOAT)(m+n)) + XYlen*Frate[c]);      
#endif

  Bias = XYlen * FLBias[c] + MatchBias[c];
}


/* short for Send((X+Y)/2,m,n) + S(X,Y,1,1) - S((X+Y)/2,(X+Y)/2,1,1) */
static inline PFLOAT Sbnd(PFLOAT X,
			  PFLOAT Y,
			  PFLOAT m_n,
			  int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  PFLOAT err = X - Y; 
  PFLOAT XYsum = X + Y;
  PFLOAT varXY = varF2[c] + var2[c] * XYsum;
  if(QUADRATIC_VARIANCE)
    varXY += varR2[c] * (X*X + Y*Y);

#if USE_MIC && USE_PFLOAT
  return (max(PZERO, X - resKB2[c]) + max(PZERO, Y - resKB2[c])) * LBias[c] + EndBias[c] - err*err*reciprocal_mic(varXY) - MissPenalty[c] * m_n;
#else
  return (max(PZERO, X - resKB2[c]) + max(PZERO, Y - resKB2[c])) * LBias[c] + EndBias[c] - err*err/(varXY) - MissPenalty[c] * m_n;
#endif
}

static inline FLOAT Sbnd(FLOAT X,
			 FLOAT Y,
			 int m,
			 int n,
			 int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  return Sbnd((PFLOAT)X, (PFLOAT)Y, (PFLOAT)(m+n), c);
}


static inline PFLOAT Send(PFLOAT X,
			  PFLOAT m_n,
			  int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  return max(PZERO, X - resKB2[c]) * LBias2[c] + EndBias[c] - MissPenalty[c] * m_n;
}

static inline FLOAT Send(FLOAT X,
			 int m,
			 int n,
			 int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  return Send((PFLOAT)X, (PFLOAT)(m+n), c);
}

#if USE_MIC

static inline __m512 Send512(__m512 X,
			     __m512 m_n,
			     int c)
{
  if(DEBUG>=2) assert(c==0 || c==1);
  return _mm512_add_ps(_mm512_mul_ps(_mm512_gmax_ps(_mm512_setzero_ps(), _mm512_sub_ps(X, v512_resKB2[c])), v512_LBiasD2[c]),
		       _mm512_sub_ps(v512_FBiasD_2[c], _mm512_mul_ps(v512_MissPenalty[c], _mm512_sub_ps(m_n, _mm512_set1_ps(2.0)))));
}
#endif // USE_MIC

static inline PFLOAT Sm(int J, int I)
{
  return PZERO;
}

/* 2-color misalignment penalty (based on the rightmost alignmed sites of each color) */
static inline PFLOAT SMA(PFLOAT X, /* difference between rightmost aligned sites of each color on Xmap */
			 PFLOAT Y, /* difference between rightmost aligned sites of each color on Ymap */
			 PFLOAT len, /* average of rightmost aligned sites on Query & Reference from left end) */
			 int Xorientation)
{
  if(DEBUG>=2 && COLOR_SHIFT) assert(len > 0.0);

  PFLOAT err = X - Y + (Xorientation ? MA_mean2 : 0.0);/* NOTE : if orientation is same, the misalignment bias cancels out */
  PFLOAT var = FM2var;
  if(COLOR_SHIFT)
    var += SM2var * len;

  if(DEBUG>=2) assert(var > 0.0);

  return RefPen12 - err*err/var - 0.5*log(var);
}

#endif // NGENTIGPAIRSCORE2_H
