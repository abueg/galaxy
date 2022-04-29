#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <omp.h>


#ifdef VALGRIND
#include <valgrind/memcheck.h>
#endif

#include "parameters.h"
#include "Ccontig.h"
#include "timers.h"
#include "probeval.h"

#define TIME_VERB 0 /* more timing statistics */
#define CNT_VERB 0 /* more count statistics */

#undef DIAGONAL
#define DIAGONAL 1 // Reduce virtual memory used with FLAT_MEMORY

#undef MIN_MEM
#define MIN_MEM QP_MIN_MEM // Minimize virtual memory used with FLAT_MEMORY (larger values reduce memory use but call malloc/free more often

#undef CACHE
#define CACHE 1 // improve cache usage of memory by avoiding changing start of arrays (slight increase in total memory use)
// CACHE uses ALIGN and ALIGN2 defined in probeval.h

#define TWO_LOOPS 0 /* Split main "omp for" into 2 loops, handling slow iterations with setlimit() or alignment in 2nd loop only */ /* NOTE : does not seem to help with OMP_WAIT_POLICY PASSIVE */

#define MIC_DEBUG 0 // debug qprobeval3D_GTH_mm512() by comparing results with non-intrinsic version

#define MISSCNT_FIX 0 /* Make Lij,Rij follow convention in refalign, so they correctly call actual overlap endpoint (rather than Viterbi coverage endpoint) */

#define TRACE_ID -1LL

#define QDEBUG 0 /* debug thread order independed computation (by comparing with result from older code) */
 
static int qprobevalcnt = 0;

static Ident qprobeval_cpp_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/qprobeval.cpp 10733 2020-03-11 17:58:12Z tanantharaman $");

#define URDEBUG ((DEBUG>=2 || TRACE_ID > 0) ? 1 : 0) /* better tracking of UR values (slower) */

namespace probeval {

class QA {
public:
  int Kstride;/**< distance between K index increments <= N * M */
  int Istride;/**< distance between I index increments <= M */
  int kmax;/**< largest value of Kmax[IMIN..IMAX] for current molecule m */
  /* NOTE : the J index stride is always 1 */

#define ARRAY(type, field) type *field##_base; inline type& field(int I, int K, int J) { return field##_base[I*Istride + K*Kstride + J];}

  ARRAY(MFLOAT1, LR); /**< alignment likelihood ratio for all alignments between Y[0..I] & X[0..J] with Y[I-K..I] aligned with X[J] */
  ARRAY(MFLOAT1, bLR);/**< best alignment contribution to LR */
  ARRAY(MFLOAT1, UL);/** only used with DEFER_BP : initial LR score based on Unaligned left end only, need to redo the computation of LR */
#if URDEBUG
  ARRAY(MFLOAT1, UR);/**< Used for debugging only */
#endif
  ARRAY(int, G);/**< previous site in Y in best alignment (-1 if none) */
  ARRAY(int, T);/**< sites Y[G-T .. G] were merged to align with X[H] */
  ARRAY(int, H);/**< previous site in X in best alignment (-1 if none) */
  ARRAY(int, Lij);/**< rightmost site/end in Y that overlaps X when Y[I-K..I] is aligned with X[J] (equals I-K if left end score (UL) was an outlier) */
  ARRAY(int, Rij);/**< rightmost site/end in Y that overlaps X when Y[I-K..I] is aligned with X[J] (equals I if right end score (UR) was an outlier)*/
#undef ARRAY

 QA() {
    LR_base = bLR_base = NULL;
    UL_base = NULL;
#if URDEBUG
    UR_base = NULL;
#endif
    G_base = T_base = H_base = Lij_base = Rij_base = NULL;
  };
};

static long long DPsizF = 3 + (URDEBUG ? 1 : 0);/* number of MFLOAT fields in RecA */
static long long DPsizI = 5;/* number of int fields in RecA */

#if USE_MIC || USE_AVX512
#define _mm512_LOADU_PS _mm512_loadu_ps
//  static inline __m512 _mm512_LOADU_PS(void const *a)
//  {
//    return _mm512_loadu_ps(a);
//  }
#endif

int verb_numthreads=1;

static void LRtrace(int m, int I, int K, int J, int M, int N, int Gmin, int Gmin2, int IMIN, int Hmin, int Hmin2, int *Jmin, double Yik, FLOAT *Xm, FLOAT *Y, MFLOAT1 LR, MFLOAT1 Lscale)
{
  if((VERB>=3 && m==0) || ((LTRACE /*|| 1*/) && MDEBUG_I >= 0 && MDEBUG_M>=0 && N==MDEBUG_N3 && m==MDEBUG_M && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)){
    printf("m=%d,I=%d,K=%d,J=%d,M=%d:Gmin=%d,%d(IMIN=%d),Hmin=%d,%d(Jmin[Gmin]=%d,%d):Yik-Y[Gmin]=%0.3f,X[J]-X[Hmin]=%0.3f,LR=%0.8e,Lscale=%0.8e\n",
	   m,I,K,J,M,Gmin,Gmin2,IMIN,Hmin,Hmin2,Jmin[Gmin],Jmin[Gmin2],Yik-Y[Gmin],Xm[J]-Xm[Hmin],LR,Lscale);
    fflush(stdout);
  }
}

static void LRtrace2(int m, int I, int K, int J, int G, int M, int N, int Gmin, int Hmin, int *Kmax, int hmin, int hmax, int *Jmin, int *Jmax, MFLOAT1 Scale)
{
  if(VERB>=3 || ((LTRACE /* || 1*/) && MDEBUG_I >= 0 && MDEBUG_M>=0 && N==MDEBUG_N3 && m==MDEBUG_M && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)){
    printf("m=%d,I=%d,K=%d,J=%d:Gmin=%d,Hmin=%d:G=%d,T=%d..%d,hmin=%d,hmax=%d,Jmin[G]=%d,Jmax[G]=%d,Scale=%0.8e\n",m,I,K,J,Gmin,Hmin,G,0,Kmax[G],hmin,hmax,Jmin[G],Jmax[G],Scale);
    fflush(stdout);
  }
}

static void LRtrace3(QA &A, int m, int MD, int I, int K, int J, int Gmin, int N, int *Kmax, MFLOAT1 LR, double *LogScale)
{
  if((DEBUG && !(isfinite(LR) && LR < MaxFloat)) || /* (rverb && m==102 && N==70 && I==64 && K==0 && J==28) || */
     (LTRACE && MDEBUG_I >= 0 && MDEBUG_M>=0 && N==MDEBUG_N3 && m==MDEBUG_M && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)  ||
     (HDEBUG && rverb && m==311 && MD==2005 && ((I==78 && K==0 && (J==24 || J==23)) || (I==79 && K==0 && J==25) || (I==77 && K==0 && J==23)))){
#pragma omp critical
    {
      printf("m=%d:I=%d,K=%d(%d..%d),J=%d,G=%d..%d:LR=%0.10e,A.UL(I,K,J)=%0.10e,LogScale=%0.6f\n",
	     m,I,K,0,Kmax[I],J,Gmin,I-K-1,LR,A.UL(I,K,J),LogScale[I]);
      fflush(stdout);
      assert(isfinite(LR) && LR < MaxFloat);
    }
  }
}

/* debug code to check accuracy of FA() == FAIK() * FAJ() computation */
static void FAcheck(int m,int I,int K,int J,int G,int T,int H, int N,FLOAT *Xm, FLOAT *Y, double Yik, MFLOAT deltaX,MFLOAT deltaY,MFLOAT LRik,MFLOAT Ivar,MFLOAT FAj, const int OUT_TYPE, const int OUTLIER_OLD)
{
  MFLOAT FAijk = FA(deltaX,deltaY,J-H,I-K-G,J,I,K,T,N,VARtab[K*N + I], VARtab[T*N + G], FBiasWTm[J-H],LogFpPow[J-H-1],LogFnPow[I-K-G-1], PrBiasWT[I-K] - PrBiasWT[G], PRtab[I-K-G-1][I-K].Pr,Y,1, OUT_TYPE, OUTLIER_OLD);

  if(DEBUG && !(fabs(FAijk - FAj) <= FAijk * (USE_MFLOAT ? 1e-4 : 2e-9) + (USE_MFLOAT ? 1e-32 : 1e-305))){
    #pragma omp critical
    {
      sleep(5);
      rverb = 1;
      MFLOAT FAj1 = FAJ(deltaX,deltaY,J-H,I-K-G,FBiasWTm[J-H],LogFpPow[J-H-1],LogFnPow[I-K-G-1], PrBiasWT[I-K] - PrBiasWT[G], LRik,Ivar,1, OUT_TYPE, OUTLIER_OLD/* ,0*/);
      MFLOAT FAijk1 = FA(deltaX,deltaY,J-H,I-K-G,J,I,K,T,N,VARtab[K*N + I], VARtab[T*N + G], FBiasWTm[J-H],LogFpPow[J-H-1],LogFnPow[I-K-G-1],  PrBiasWT[I-K] - PrBiasWT[G], PRtab[I-K-G-1][I-K].Pr, Y,1, OUT_TYPE, OUTLIER_OLD);

      MFLOAT Ivar2;
      MFLOAT LRik2 = FAIK(deltaY,PRtab[I-K-G-1][I-K].Pr,I-K-G,I,K,T,N,VARtab[K*N+I],VARtab[T*N+G],Ivar2,Y);
      MFLOAT FAj2 = FAJ(deltaX,deltaY,J-H,I-K-G,FBiasWTm[J-H],LogFpPow[J-H-1],LogFnPow[I-K-G-1],PrBiasWT[I-K] - PrBiasWT[G], LRik2,Ivar2,1, OUT_TYPE, OUTLIER_OLD/* ,0*/);
      rverb = 0;
      printf("m=%d,I=%d,K=%d,J=%d,G=%d,T=%d,H=%d:LRik=%0.10e,FAj=%0.10e,FAijk=%0.10e,FAj1=%0.10e,FAijk1=%0.10e,Ivar2=%0.8e,FAik2=%0.8e,FAj2=%0.10e\n",
	     m,I,K,J,G,T,H,LRik,FAj,FAijk,FAj1,FAijk1,Ivar2,LRik2,FAj2);
      printf("Yik=%0.6f,Ygt=%0.6f,X[J]=%0.6f,X[H]=%0.6f,deltaX=%0.6f,deltaY=%0.6f,Ivar=%0.10e,Xtheta=%0.10e\n",
        Yik,Yc(Y,G,T),Xm[J],Xm[H],deltaX,deltaY,Ivar,Xtheta);
      double err = (deltaX-deltaY)*Ivar;
      printf("exp(-deltaX*Frate)=%0.10e,FpPow[J-H-1]=%0.10e,err=%0.10e,exp(-err*err)=%0.10e\n",
        exp(-deltaX*Frate),FpPow[J-H-1],err,exp(-err*err));
      printf("RefPen=%0.10e,FnPow[I-K-G-1]=%0.10e,PrTab[I-K-G-1][I-K].Pr=%0.10e,sqrt(1/(F2var+S2var*deltaY+R2var*deltaY*deltaY))=%0.10e,FAIK()=%0.10e\n",
	     RefPen,FnPow[I-K-G-1],PRtab[I-K-G-1][I-K].Pr,sqrt(1.0/(F2var+S2var*deltaY+R2var*deltaY*deltaY)),FAIK(deltaY,PRtab[I-K-G-1][I-K].Pr,I-K-G,I,K,T,N,VARtab[K*N+I], VARtab[T*N+G],Ivar,Y));
      printf("maptype=%d,outlierMax=%0.6e,OUTLIER(deltaX,deltaY)=%d,expF(-err*err)=%0.10e,Poutlier=%0.6e,PoutlierM1=%0.10e,PoutlierPIbXtheta=%0.10e,PoutlierPIBYmax=%0.10e\n",
        maptype,outlierMax,OUTLIER(deltaX,deltaY),expF(-err*err),Poutlier,PoutlierM1,PoutlierPIbXtheta,PoutlierPIBYmax);
      fflush(stdout);

      assert(fabs(FAijk - FAj) <= FAijk * (USE_MFLOAT ? 1e-4 : 2e-9) + (USE_MFLOAT ? 1e-32 : 1e-305));
    }
  }
}

static void LRcheck(QA &A, int m, int MD, int I, int K, int J, int G, int T, int H, int N, FLOAT *Xm, FLOAT *Y, MFLOAT deltaX, MFLOAT deltaY, MFLOAT LRik, MFLOAT Ivar, MFLOAT FAj, 
		     MFLOAT PMik, MFLOAT newLR, MFLOAT1 LR, MFLOAT1 LRg, MFLOAT1 Scale, const int OUT_TYPE, const int OUTLIER_OLD)
{
  if((DEBUG>=2 && !(isfinite(newLR) && (Poutlier <= 0.0 || !OUTLIER(deltaX,deltaY) || PMik <= 0.0 ||  FAj >  0.0))) ||
     /* (rverb && m==102 && N==70 && I==64 && K==0 && J==28) ||*/
     ((LTRACE /*||1*/) && MDEBUG_I >= 0 && MDEBUG_M>=0 && m==MDEBUG_M && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J && N==MDEBUG_N3 && newLR >= LR*1e-10) ||
     /* (rverb && m==MDEBUG_M && newLR >= LR*1e-10 && ((I==32 && K==0 && J==7) || (I==155 && K==0 && J==25) || (I==599 && K==0 && J==37) )) || */
     (HDEBUG && rverb && m==1 && MD==2005 && ((I==MDEBUG_I && K==0 && (J==24 || J==23)) || (I==79 && K==0 && J==25) || (I==77 && K==0 && J==23)))){
    //		    rverb = 2;
    MFLOAT fa = FA(deltaX,deltaY,J-H,I-K-G,J,I,K,T,N,VARtab[K*N + I], VARtab[T*N + G], FBiasWTm[J-H],LogFpPow[J-H-1],LogFnPow[I-K-G-1], PrBiasWT[I-K] - PrBiasWT[G], PRtab[I-K-G-1][I-K].Pr, Y,1,OUT_TYPE,OUTLIER_OLD);
    printf("m=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d:PMik=%0.10f,FAj=%0.10e,FA=%0.10e(LRik=%0.10e,Ivar=%0.10e,FpPow[%d]=%0.6f,x=%0.6f,y=%0.6f,err=%0.6f)\n\tprev AL=%0.10e(ScaleI[G]=%0.10e),AL=%0.10e(LRg=%0.10e,LRgth=%0.10e)\n\tUL=%0.10e,Ivar=%0.10e,Y[I-K..I]=%0.6f,N=%d,OUTLIER(deltaX,deltaY)=%d\n",
	   m,I,K,J,G,T,H,PMik,FAj,fa, LRik, Ivar,J-H-1,FpPow[J-H-1],deltaX,deltaY,(deltaX-deltaY)*Ivar,A.LR(G,T,H), Scale,LR + LRg*Scale,LRg*Scale,newLR*Scale,A.LR(I,K,J),
	   Ivar,Y[I]-Y[I-K], N, OUTLIER(deltaX,deltaY)?1:0);
    fflush(stdout);
    assert(isfinite(newLR) && (Poutlier <=0.0 || !OUTLIER(deltaX,deltaY) || PMik <= 0.0 || FAj > 0.0));
  }
}

static size_t V3D_totcnt = 0, V3D_rcnt = 0;

// scalar version of qprobeval3D_GTH

template<int OUT_TYPE, int OUT_OLD>
static double inline qprobeval3D_GTH(QA &A, 
				     int m, int MD, int *Jmin, int *Jmax, int M, int N, 
				     double *LogScale, double *ScaleI, double *Xm, double *Y,
				     int I, int K, double PMik, double Yik, int Gmin, int Gmin2, 
				     int J, int Hmin, int Hmin2, MFLOAT *delX, Ccontig *pcontig/* for debugging only*/)
{
  MFLOAT *pFBiasWTm = &FBiasWTmRev[-J];
  MFLOAT *pLogFpPow = &LogFpPowRev[1-J];

  MFLOAT vecR = VARtab[K*N+I];
  int IK = I-K;

  double LR = A.UL(I,K,J);

  for(int G= IK; --G >= Gmin;){
    if(DEBUG>=2) assert(fabs(ScaleI[G] - exp(LogScale[G]-LogScale[I])) < fabs(ScaleI[G]) * 1e-5 + 1e-6);
    if(DEBUG>=2 && SCALE_MONOTONIC) assert(ScaleI[G] <= 1.0);
    if(DEBUG>=2 && !SCALE_MONOTONIC) assert(ScaleI[G] <= DBL_MAX);
    int hmin = (G >= Gmin2) ? max(Hmin,Jmin[G]) : max(Hmin2, Jmin[G]);
    int hmax = min(J-1,Jmax[G]);

    //    double origLR = LR;

    if(hmin <= hmax){// NEW10
      if(VERB>=3 || (LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0)) LRtrace2(m,I,K,J,G,M,N,Gmin,Hmin,Kmax,hmin,hmax,Jmin,Jmax,ScaleI[G]*PMik);
      int n = IK - G;
      MFLOAT LogFnPowIKG = LogFnPow[n-1];
      MFLOAT PRikn = PRtab[n-1][IK].Pr;
      MFLOAT PrBiasWTikg = PrBiasWT[I-K] - PrBiasWT[G];
      if(DEBUG>=2) assert(PRikn <= 1.0);

#if DEBUG>=2
      MFLOAT LRmax = 0.0, bFAj;
      int bT,bH;
#endif

      MFLOAT1 LRg = 0.0;

      for(int T = 0; T <= Kmax[G]; T++){
	MFLOAT Ivar;
	MFLOAT deltaY = Yik - Yc(Y,G,T);
	MFLOAT LRik = FAIK(deltaY,PRikn, n,I,K,T,N,vecR, VARtab[T*N+G],Ivar,Y);
	if(DEBUG>=2) assert(LRik <= RefPen * 1.00001);
	MFLOAT1 *LRgt = &A.LR(G,T,0);

#if DEBUG < 2
#pragma omp simd reduction(+:LRg)
#endif
	for(int H= hmin; H <= hmax; H++){/* critical loop : may have to be combined with G,T loops to accumulate 16 iterations at a time */
	  if(DEBUG>=2) assert(G-T >= 1 && G <= N && H >= 1 && H <= M);
	  if(DEBUG>=2) assert(fabs((Xm[J] - Xm[H]) - delX[H]) < Xm[J]*((MFLOAT)1e-6));
	  MFLOAT deltaX = delX[H]; // Xm[J] - Xm[H];
	  if(DEBUG>=2) assert(pFBiasWTm[H] == FBiasWTm[J-H]);
	  if(DEBUG>=2) assert(pLogFpPow[H] == LogFpPow[J-H-1]);

	  MFLOAT FAj = FAJ(deltaX, deltaY, J-H, n, pFBiasWTm[H], pLogFpPow[H], LogFnPowIKG, PrBiasWTikg, LRik, Ivar, 1, OUT_TYPE, OUT_OLD);
	  if(DEBUG>=2) FAcheck(m,I,K,J,G,T,H,N,Xm,Y,Yik,deltaX,deltaY,LRik,Ivar,FAj, OUT_TYPE, OUT_OLD);
	  MFLOAT1 newLR = LRgt[H] * FAj;
	
	  if(VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && I == 14997 && K == 0 && J == 138){
	    printf("m=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d:x=%0.4f,y=%0.4f:FAj=%0.8e(PMik=%0.8e,LRik=%0.8e,Ivar=%0.8e),prev LR=%0.8e(Scale(G,I)=%0.8e,LS[I]=%0.8f),LRg=%0.8e,cum=%0.8e(delta=%0.8e)\n",
		   m,I,K,J,G,T,H,deltaX,deltaY,FAj,PMik,LRik,Ivar,A.LR(G,T,H), ScaleI[G], LogScale[I], LRg * ScaleI[G] * PMik, (LR + (LRg + newLR) * ScaleI[G]) * PMik, newLR * PMik * ScaleI[G]);
	    fflush(stdout);
	  }

#if DEBUG >= 2
	  if(DEBUG>=2 && !(isfinite(newLR) && newLR < MaxFloat)){
#pragma omp critical
	    {
	      printf("m=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d:newLR=%0.8e,A.LR(G,T,H)=%0.8e,FAj=%0.8e,deltaX=%0.6f,deltaY=%0.6f,LRik=%0.8e,Ivar=%0.8f\n",
		     m,I,K,J,G,T,H,newLR,A.LR(G,T,H),FAj,deltaX,deltaY,LRik,Ivar);
	      fflush(stdout);
	      sleep(10);
	      rverb = 1;
	      MFLOAT FAj2 = FAJ(deltaX, deltaY, J-H, n, FBiasWTm[J-H], LogFpPow[J-H-1], LogFnPowIKG, LRik, Ivar, 1, OUT_TYPE, OUT_OLD);
	      rverb = 0;
	      printf("\t FAj2= %0.8e\n",FAj2);
	      fflush(stdout);
	      assert(isfinite(newLR) && newLR < MaxFloat);
	    }
	  }
#endif

#if DEBUG>=2
	  if(newLR > LRmax){
	    LRmax = newLR;
	    bFAj = FAj;
	    bT = T;
	    bH = H;
	  }
#endif
	  LRg += newLR;
	  if(DEBUG>=2) assert(isfinite(LRg) && LRg < MaxFloat);
	  if(DEBUG>=2 || (LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0) || HDEBUG) 
	    LRcheck(A,m,MD,I,K,J,G,T,H,N,Xm,Y,deltaX,deltaY,LRik,Ivar,FAj,PMik,newLR,LR,LRg,ScaleI[G]*PMik, OUT_TYPE, OUT_OLD);
	} // for H = hmax .. hmin
      } // for T = 0 .. Kmax[G]

      if(DEBUG>=2) assert(isfinite(LRg) && LRg < MaxFloat);
      LR += LRg * ScaleI[G];/* ScaleI[G] == exp(LogScale[G] - LogScale[I]) */
    }// if(hmin <= hmax)
    
    //    if(DEBUG>=2 && hmin > hmax) assert(LR == origLR);

  } // for G = I-K-1 .. Gmin
  
  if(DEBUG>=2 && J <= 1) assert(LR == A.UL(I,K,J));

  return LR;
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT
// AVX256 version of qprobeval3D_GTH

template<int OUT_TYPE, int OUT_OLD>
static double inline qprobeval3D_GTH_mm256(QA &A, 
					   int m, int MD, int *Jmin, int *Jmax, int M, int N, 
					   double *LogScale, double *ScaleI, double *Xm, double *Y,
					   int I, int K, double PMik, double Yik, int Gmin, int Gmin2, 
					   int J, int Hmin, int Hmin2, MFLOAT *delX)
{
  MFLOAT *pFBiasWTm = &FBiasWTmRev[-J];
  MFLOAT *pLogFpPow = &LogFpPowRev[1-J];

  MFLOAT vecR = VARtab[K*N+I];
  int IK = I-K;

  double LR = A.UL(I,K,J);
  if(DEBUG>=2) assert(!IsNan(LR) && LR >= 0.0);

  for(int G= IK; --G >= Gmin;){
    int hmin = (G >= Gmin2) ? max(Hmin,Jmin[G]) : max(Hmin2, Jmin[G]);
    int hmax = min(J-1,Jmax[G]);
    if(DEBUG/* HERE >=2 */) assert(hmax-hmin+1 <= DELTA_MAX + EXTEND);
    int n = IK - G;
    float LogFnPowIKG = LogFnPow[n-1];
    __m256 v_LogFnPowIKG = _mm256_set1_ps(LogFnPowIKG);

    float PrBiasWTikg = PrBiasWT[I-K] - PrBiasWT[G];
    __m256 v_PrBiasWTikg = _mm256_set1_ps(PrBiasWTikg);

    MFLOAT PRikn = PRtab[n-1][IK].Pr;

    MFLOAT LRg = 0.0;

    __m256 v_LRg = _mm256_setzero_ps();

    for(int T = 0; T <= Kmax[G]; T++){
      float Ivar;
      float deltaY = Yik - Yc(Y,G,T);
      float LRik = FAIK(deltaY,PRikn, n,I,K,T,N,vecR, VARtab[T*N+G],Ivar,Y);
      float *LRgt = &A.LR(G,T,0);
      
      __m256 v_Ivar = _mm256_set1_ps(Ivar);
      __m256 v_deltaY = _mm256_set1_ps(deltaY);
      __m256 v_LRik = _mm256_set1_ps(LRik);

      int hmaxV = hmax - VEC_MSIZQ + 1;
      int H = hmin;
      for(; H <= hmaxV; H += VEC_SIZ){/* critical loop : vectorized */
	// NOTE : delX[],pFBiasWTm[] & pLogFpPow[],LRgt[] arrays are padded by at least 7 elements initialized to some finite value
	__m256 v_deltaX = _mm256_loadu_ps(&delX[H]); // delX[H] == Xm[J] - Xm[H]
	__m256 v_FAj = FAJ_mm256<1,OUT_TYPE,OUT_OLD>(v_deltaX, v_deltaY, _mm256_loadu_ps(&pFBiasWTm[H]), _mm256_loadu_ps(&pLogFpPow[H]), v_LogFnPowIKG, v_PrBiasWTikg, v_LRik, v_Ivar);
	
	/* only (hmax - H + 1) low order values in vector v_FAj are valid, so force the rest to zero :
	   ZeroHigh_mm256[n = 0.. DELTA_MAX + EXTEND] is a precomputed _mm256 float vector with low order min(n,8) values 1.0 all others 0.0 */
	v_LRg += _mm256_loadu_ps(&LRgt[H]) * v_FAj * ZeroHigh_mm256[hmax - H + 1];

#if DEBUG >= 2
	float LRgV[8]; _mm256_storeu_ps(LRgV, v_LRg);
	float FAjV[8]; _mm256_storeu_ps(FAjV, v_FAj);
	float ZeroHighV[8]; _mm256_storeu_ps(ZeroHighV, ZeroHigh_mm256[hmax - H + 1]);
	float deltaXV[8]; _mm256_storeu_ps(deltaXV, v_deltaX);
	for(int i = 0; i < 8; i++){
	  if(DEBUG && !(!IsNan(LRgV[i]) && LRgV[i] >= 0.0f)){
	    #pragma omp critical
	    {
	      printf("qprobeval3D_GTH_mm256: m=%d,MD=%d,M=%d,N=%d,I=%d,K=%d,J=%d,G=%d,T=%d,H=%d(hmax=%d): i=%d, LRgV[i]= %0.8e (ZeroHigh_mm256[%d/%d][i]= %0.8e, LRgt[H+i]= %0.8e, FAj[i]= %0.8e, deltaX[i]= %0.6e)\n",
		     m,MD,M,N,I,K,J,G,T,H,hmax,i,LRgV[i], hmax - H + 1, DELTA_MAX+EXTEND, ZeroHighV[i], LRgt[H+i], FAjV[i], deltaXV[i]);
	      printf("\t Ivar= %0.8e, deltaY= %0.6e, LRik= %0.8e, pFBiasWTm[H+i]= %0.8e, pLogFpPow[H+i]= %0.8e, LogFnPowIKG= %0.8e, delX[H+i]= %0.8f, &delX[H+i]= x%p\n",
		     Ivar, deltaY, LRik, pFBiasWTm[H+i], pLogFpPow[H+i], LogFnPowIKG,delX[H+i], &delX[H+i]);
	      fflush(stdout);
	      assert(!IsNan(LRgV[i]) && LRgV[i] >= 0.0f);
	    }
	  }
	}
#endif
      } // for H = hmin .. hmaxV

#if VEC_MSIZQ > 1
      for(int h = H; h <= hmax; h++){/* tail end of critical loop : not vectorized */
	MFLOAT deltaX = delX[h]; // Xm[J] - Xm[h];
	MFLOAT FAj = FAJ(deltaX, deltaY, J-H, n, pFBiasWTm[h], pLogFpPow[h], LogFnPowIKG, LRik, Ivar, 1, OUT_TYPE, OUT_OLD);

	LRg += LRgt[h] * FAj;

	if(DEBUG>=2) assert(!IsNan(LRg) && LRg >= 0.0);
      } // for h = H .. hmax
#endif
    } // for T = 0 .. Kmax[G]

    LRg += horizontal_add_mm256(v_LRg);
    if(DEBUG>=2) assert(!IsNan(LRg) && LRg >= 0.0);

    LR += LRg * ScaleI[G];
    if(DEBUG>=2) assert(!IsNan(LR) && LR >= 0.0);
  } // for G = I-K-1 .. Gmin
  
#if DEBUG>=2
  assert(!IsNan(LR) && LR >= 0.0);

  double sLR = qprobeval3D_GTH<OUT_TYPE,OUT_OLD>(A,m,MD,Jmin,Jmax,M,N,LogScale,ScaleI,Xm,Y,I,K,PMik,Yik,Gmin,Gmin2,J,Hmin,Hmin2,delX,pcontig);
  if(DEBUG && !(fabs(sLR - LR) < fabs(sLR) * 1e-5 + 1e-28)){
    #pragma omp critical
    {
      printf("m=%d/%d:I=%d,K=%d,J=%d,M=%d,N=%d,Yik= %0.4f,Gmin=%d,Gmin2=%d,Hmin=%d,Hmin2=%d:LR= %0.8e, sLR= %0.8e\n",
	     m,MD,I,K,J,M,N,Yik,Gmin,Gmin2,Hmin,Hmin2,LR,sLR);
      fflush(stdout);
      assert(fabs(sLR - LR) < fabs(sLR) * 1e-5 + 1e-28) ;
    }
  }
#endif

  return LR;
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT
// MIC version of qprobeval3D_GTH


template<int OUT_TYPE, int OUT_OLD>
static double inline qprobeval3D_GTH_mm512(QA &A, 
					   int m, int MD, int *Jmin, int *Jmax, int M, int N, 
					   double *LogScale, double *ScaleI, double *Xm, double *Y,
					   int I, int K, double PMik, double Yik, int Gmin, int Gmin2, 
					   int J, int Hmin, int Hmin2, MFLOAT *delX, Ccontig *pcontig/* for debugging only*/)
{
  MFLOAT *pFBiasWTm = &FBiasWTmRev[-J];
  MFLOAT *pLogFpPow = &LogFpPowRev[1-J];

  MFLOAT vecR = VARtab[K*N+I];
  int IK = I-K;

  double LR = A.UL(I,K,J);
#if MIC_DEBUG
  double dLR = LR;
#endif
  if(DEBUG>=2 || MIC_DEBUG) assert(!IsNan(LR) && LR >= 0.0);
  //   __m512 v_LR = _mm512_setzero_ps();

  for(int G= IK; --G >= Gmin;){
    int hmin = (G >= Gmin2) ? max(Hmin,Jmin[G]) : max(Hmin2, Jmin[G]);
    int hmax = min(J-1,Jmax[G]);
    if(DEBUG/* HERE >=2 */) assert(hmax-hmin+1 <= DELTA_MAX + EXTEND);
    int n = IK - G;
    float LogFnPowIKG = LogFnPow[n-1];
    __m512 v_LogFnPowIKG = _mm512_set1_ps(LogFnPowIKG);

    float PrBiasWTikg = PrBiasWT[IK] - PrBiasWT[G];
    __m512 v_PrBiasWTikg = _mm512_set1_ps(PrBiasWTikg);

    MFLOAT PRikn = PRtab[n-1][IK].Pr;

    MFLOAT LRg = 0.0;
#if MIC_DEBUG
    MFLOAT dLRg = 0.0;
#endif    
    __m512 v_LRg = _mm512_setzero_ps();

    for(int T = 0; T <= Kmax[G]; T++){
      float Ivar;
      float deltaY = Yik - Yc(Y,G,T);
      float LRik = FAIK(deltaY,PRikn, n,I,K,T,N,vecR, VARtab[T*N+G],Ivar,Y);
      float *LRgt = &A.LR(G,T,0);
      
      __m512 v_Ivar = _mm512_set1_ps(Ivar);
      __m512 v_deltaY = _mm512_set1_ps(deltaY);
      __m512 v_LRik = _mm512_set1_ps(LRik);

      int hmaxV = hmax - VEC_MSIZQ + 1;
      int H = hmin;
      for(; H <= hmaxV; H += VEC_SIZ){/* critical loop : vectorized */
	// NOTE : delX[],pFBiasWTm[] & pLogFpPow[],LRgt[] arrays are padded by at least 15 elements initialized to some finite value
	__m512 v_deltaX = _mm512_LOADU_PS(&delX[H]); // delX[H] == Xm[J] - Xm[H]
	__m512 v_FAj = FAJ_mm512<1,OUT_TYPE,OUT_OLD>(v_deltaX, v_deltaY, _mm512_LOADU_PS(&pFBiasWTm[H]), _mm512_LOADU_PS(&pLogFpPow[H]), v_LogFnPowIKG, v_PrBiasWTikg, v_LRik, v_Ivar);

        __m512 v_LRg_inc = _mm512_mul_ps(_mm512_LOADU_PS(&LRgt[H]), v_FAj);

        /* only (hmax - H + 1) low order values in vector v_FAj or v_LRg_inc are valid, so force the rest to zero */
	__m512 v_LRg_masked = _mm512_mask_blend_ps(MaskLow_mask16[hmax-H+1], v_LRg_inc, _mm512_setzero_ps());

#if MIC_DEBUG
	__m512 v_LRg_orig = v_LRg;
	float dLRg_orig = dLRg;
#endif
	v_LRg = _mm512_add_ps(v_LRg, v_LRg_masked);

#if MIC_DEBUG
	float FAj[16], FAj_inc[16], FAjT[16], dLRg_cum[16];
	for(int h = H; h <= hmax && h < H + 16; h++){
	  MFLOAT deltaX = delX[h]; // Xm[J] - Xm[h];
	  FAj[h-H] = FAJ<1,OUT_TYPE,OUT_OLD>(deltaX, deltaY, pFBiasWTm[h], pLogFpPow[h], LogFnPowIKG, PrBiasWTikg, LRik, Ivar);
	  if(!(isfinite(FAj[h-H]))){
	    #pragma omp critical
	    {
	      printf("m=%d/%d:M=%d,N=%d,I=%d,K=%d,J=%d,G=%d(%d..%d),T=%d/%d,H=%d(%d..%d):h=%d:FAj[h-H]= %0.8e\n",
		     m,MD,M,N,I,K,J,G,IK-1,Gmin,T,Kmax[G],H,hmin,hmax,h,FAj[h-H]);
	      printf("deltaX= %0.6f, deltaY=%0.6f, pFBiasWTm[h]=%0.8e, pLogFpPow[h]= %0.8e, LogFnPowIKG= %0.8e, LRik= %0.8e, Ivar= %0.8e\n",
		     deltaX, deltaY, pFBiasWTm[h], pLogFpPow[h], LogFnPowIKG, LRik, Ivar);
	      fflush(stdout);
	      assert(isfinite(FAj[h-H]));
	      assert(0);
	    }
	  }

	  FAj_inc[h-H] = LRgt[h] * FAj[h-H];
	  dLRg += FAj_inc[h-H];
	  dLRg_cum[h-H] = dLRg;
	} // for h = H .. hmax

	float FAjV[16] __attribute__((aligned(64))); _mm512_store_ps(FAjV,v_FAj);
        float LRg_incV[16] __attribute__((aligned(64))); _mm512_store_ps(LRg_incV, v_LRg_inc);
	float LRg_maskedV[16] __attribute__((aligned(64))); _mm512_store_ps(LRg_maskedV, v_LRg_masked);
	for(int i = 0; i <= hmax - H && i < 16; i++){
	  if(DEBUG && !(fabs(FAj[i] - FAjV[i]) < fabs(FAj[i]) * 1e-5f + 1e-28f && LRg_incV[i] == LRg_maskedV[i])){
	    #pragma omp critical
	    {
	      printf("m=%d/%d:M=%d,N=%d,I=%d,K=%d,J=%d,G=%d(%d..%d),T=%d/%d,H=%d(%d..%d):i=%d,FAj[i]= %0.8e (should be %0.8e)\n",
		     m,MD,M,N,I,K,J,G,IK-1,Gmin,T,Kmax[G],H,hmin,hmax,i,FAjV[i], FAj[i]);
	      printf("\t x= %0.4f, y= %0.4f, FBiasWTm[J-H-i]= %0.6e, LogFpPow[J-1-H-i]= %0.6e, LogFnPowIKG= %0.6e, LRik= %0.6e, Ivar= %0.6e, LRgt[H+i]= %0.6e\n",
	         delX[H+i], deltaY, pFBiasWTm[H+i], pLogFpPow[H+i], LogFnPowIKG, LRik,Ivar, LRgt[H+i]);
              printf("\t LRg_inc[i]= %0.8e, LRg_masked[i]= %0.8e\n", LRg_incV[i], LRg_maskedV[i]);
	      fflush(stdout);
	      assert(fabs(FAj[i] - FAjV[i]) < fabs(FAj[i]) * 1e-5f + 1e-28f && LRg_incV[i] == LRg_maskedV[i]);
	    }
	  }
	}
	for(int i = hmax - H + 1; i < 16; i++){
	  if(DEBUG && LRg_maskedV[i] != 0.0f){
	    #pragma omp critical
	    {
	      printf("m=%d/%d:M=%d,N=%d,I=%d,K=%d,J=%d,G=%d(%d..%d),T=%d/%d,H=%d(%d..%d):i=%d,LRg_masked[i]= %0.8e (should be 0), mask= %0x\n",
		     m,MD,M,N,I,K,J,G,IK-1,Gmin,T,Kmax[G],H,hmin,hmax,i,LRg_maskedV[i], (unsigned int)MaskLow_mask16[hmax-H+1]);
	      fflush(stdout);
	      assert(LRg_maskedV[i] == 0.0f);
	    }
	  }
	}

	/* check if cumulative values are correct : dLRg should match LRg + reduce_add_ps(v_LRg) */
	float sLRg = LRg + _mm512_reduce_add_ps(v_LRg);
	if(DEBUG && !(fabs(sLRg - dLRg) < fabs(dLRg) * 1e-5f + 1e-37f)){
	  float LRgVorig[16] __attribute__((aligned(64))); _mm512_store_ps(LRgVorig, v_LRg_orig);
	  float LRgV[16] __attribute__((aligned(64))); _mm512_store_ps(LRgV, v_LRg);
	  #pragma omp critical
	  {
	    float sLRg_orig = LRg + _mm512_reduce_add_ps(v_LRg_orig);
	    printf("m=%d/%d:M=%d,N=%d,I=%d,K=%d,J=%d,G=%d(%d..%d),T=%d/%d,H=%d(%d..%d):LRg + hsum(v_LRg)= %0.8e -> %0.8e, dLRg= %0.8e -> %0.8e\n",
		   m,MD,M,N,I,K,J,G,IK-1,Gmin,T,Kmax[G],H,hmin,hmax,sLRg_orig, sLRg, dLRg_orig, dLRg);
	    for(int i = 0; i <= hmax - H && i < 16; i++){
	      FAjT[i] = FAJ<1,OUT_TYPE,OUT_OLD>(delX[H + i], deltaY, pFBiasWTm[H+i], pLogFpPow[H+i], LogFnPowIKG, PrBiasWTikg,LRik, Ivar);
	      printf("\t i=%d: FAj[i]= %0.8e(%0.8e,%0.8e), LRgt[H+i]= %0.8e, v_LRg_inc[i]= %0.8e(%0.8e,cum=%0.8e), v_LRg[i]= %0.8e -> %0.8e\n",
		     i, FAjV[i], FAj[i], FAjT[i], LRgt[H+i], LRg_incV[i], FAj_inc[i], dLRg_cum[i], LRgVorig[i], LRgV[i]);
	    }
	    fflush(stdout);
	    assert(fabs(sLRg - dLRg) < fabs(dLRg) * 1e-5f + 1e-37f);
	  }
	}
#endif	 // DEBUG
      } // for H = hmin .. hmaxV

#if VEC_MSIZQ > 1
      for(int h = H; h <= hmax; h++){/* tail end of critical loop : not vectorized */
	MFLOAT deltaX = delX[h]; // Xm[J] - Xm[h];
	MFLOAT FAj = FAJ(deltaX, deltaY, J-H, n, pFBiasWTm[h], pLogFpPow[h], LogFnPowIKG, LRik, Ivar, 1, OUT_TYPE, OUT_OLD);

	LRg += LRgt[h] * FAj;

#if MIC_DEBUG
	dLRg += LRgt[h] * FAj;
#endif
	if(DEBUG>=2) assert(!IsNan(LRg) && LRg >= 0.0);
      } // for h = H .. hmax
#endif
    } // for T = 0 .. Kmax[G]

#if MIC_DEBUG
    double pLR = LR;
    double pdLR = dLR;
    if(DEBUG && !(fabs(pLR - pdLR) < fabs(pdLR) * 1e-5f + 1e-16f)){
      #pragma omp critical
      {
        printf("m=%d/%d:M=%d,N=%d,I=%d,K=%d,J=%d,G=%d(%d..%d):PRE:LR= %0.8e, dLR= %0.8e\n",
	  m,MD,M,N,I,K,J,G,IK-1,Gmin, pLR, pdLR);
	fflush(stdout);
        assert(fabs(pLR - pdLR) < fabs(pdLR) * 1e-5f + 1e-16f);
      }
    }

    float sLRg = _mm512_reduce_add_ps(v_LRg) + LRg;
    if(DEBUG && !(fabs(sLRg - dLRg) < fabs(dLRg) * 1e-5f + 1e-37f)){
      #pragma omp critical
      {
        printf("m=%d/%d:M=%d,N=%d,I=%d,K=%d,J=%d,G=%d(%d..%d):INC:LRg + hsum(v_LRg)= %0.8e, dLRg= %0.8e\n",
	  m,MD,M,N,I,K,J,G,IK-1,Gmin, sLRg, dLRg);
	fflush(stdout);
        assert(fabs(sLRg - dLRg) < fabs(dLRg) * 1e-5f + 1e-37f);
      }
    }
#endif

    LRg += _mm512_reduce_add_ps(v_LRg);// HERE : check if doing a float to double conversion, following by reduce_add outside the G loop is faster
    LR += LRg * ScaleI[G];

#if MIC_DEBUG
    dLR += dLRg * ScaleI[G];
    if(DEBUG && !(fabs(LR - dLR) < fabs(dLR) * 1e-5f + 1e-16f)){
      #pragma omp critical
      {
        printf("m=%d/%d:M=%d,N=%d,I=%d,K=%d,J=%d,G=%d(%d..%d):\n\t POST:LR= %0.8e, dLR= %0.8e\n",
	  m,MD,M,N,I,K,J,G,IK-1,Gmin, LR, dLR);
	printf("\t PRE: LR= %0.8e, dLR= %0.8e\n",pLR,pdLR);
	printf("\t INC: LRg + hsum(v_LRg)= %0.8e, dLRg= %0.8e\n", sLRg, dLRg);
	printf("\t ScaleI[G]= %0.8e\n",ScaleI[G]);
	fflush(stdout);
        assert(fabs(LR - dLR) < fabs(dLR) * 1e-5f + 1e-16f);
      }
    }
#endif
  } // for G = I-K-1 .. Gmin
  
#if MIC_DEBUG
  assert(!IsNan(LR) && LR >= 0.0);

  double sLR = qprobeval3D_GTH<OUT_TYPE,OUT_OLD>(A,m,MD,Jmin,Jmax,M,N,LogScale,ScaleI,Xm,Y,I,K,PMik,Yik,Gmin,Gmin2,J,Hmin,Hmin2,delX,pcontig);
  if(DEBUG && !(fabs(sLR - LR) < fabs(sLR) * 1e-4f + 1e-16f)){
    #pragma omp critical
    {
      printf("m=%d/%d:I=%d,K=%d,J=%d,M=%d,N=%d,Yik= %0.4f,Gmin=%d,Gmin2=%d,Hmin=%d,Hmin2=%d:LR= %0.8e, sLR= %0.8e (dLR= %0.8e)\n",
	     m,MD,I,K,J,M,N,Yik,Gmin,Gmin2,Hmin,Hmin2,LR,sLR, dLR);
      fflush(stdout);
      assert(fabs(sLR - LR) < fabs(sLR) * 1e-4f + 1e-16f) ;
    }
  }
#endif

  return LR;
}

#endif	  // ! (VECTOR_AVX || VECTOR_MIC)

template<int OUT_TYPE, int OUT_OLD>
static void __attribute__ ((noinline))
qprobeval3D(QA &A, int m, int MD, int IMIN, int IMAX, int *Kmax, int *Jmin, int *Jmax, int M, int N, double *LogScale, double *ScaleI, double *Xm, MFLOAT **delXm, double *Y, Ccontig *pcontig)
{
  LogScale[IMIN-1] = 0.0;
  ScaleI[IMIN-1] = 1.0;

  for(int I= IMIN; I <= IMAX; I++){
    LogScale[I] = LogScale[I-1];
    ScaleI[I] = 1.0;/* ScaleI[G=IMIN..I] = exp(LogScale[G]-LogScale[I]) */
    if(DEBUG) assert(LogScale[I] >= 0.0);
    MFLOAT1 Lscale = exp(-LogScale[I]);
    if(DEBUG) assert(isfinite(Lscale) && Lscale <= 1.0);
    MFLOAT1 maxLR = 0.0;
    for(int K = 0; K <= Kmax[I]; K++){
      MFLOAT PMik = PRtab[K][I].PM;
      if(DEBUG>=2) assert(isfinite(PMik) && PMik >= 0.0);
      double Yik = Yc(Y,I,K);
      int Gmin = max(IMIN, I-K - deltaExtYRef);
      int Gmin2 = max(IMIN, I-K - DELTA_Y);
      for(int J= Jmin[I]; J <= Jmax[I]; J++){
	A.bLR(I,K,J) *= Lscale;/* original value is unscaled, so must be multiplied by Lscale */
	A.UL(I,K,J) = A.LR(I,K,J) * Lscale; // remember original scaled LR value, so we can rerun the code in this iteration later to determine bestG,bestH,bestT
	if(DEBUG>=2) assert(isfinite(A.UL(I,K,J)) && A.UL(I,K,J) < MaxFloat);
	int Hmin = max(1, J - deltaExtXRef);
	int Hmin2 = max(1, J - DELTA_X);
	if(VERB>=3 || (LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0)) LRtrace(m,I,K,J,M,N,Gmin,Gmin2,IMIN,Hmin,Hmin2,Jmin,Yik,Xm,Y,(double)A.UL(I,K,J),Lscale);


	double LR;
	if(J <= 1)// NEW10 : delXm[1] is NOT defined and J==1 has no valid values of H <= J-1
	  LR = A.UL(I,K,J);
	else {
	  MFLOAT *delX = &delXm[J][deltaExtXRef-J];

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT
	  LR = qprobeval3D_GTH_mm256<OUT_TYPE,OUT_OLD>(A, m, MD, Jmin, Jmax, M, N, LogScale, ScaleI, Xm, Y, I, K, PMik, Yik, Gmin, Gmin2, J, Hmin, Hmin2, delX);
#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT
	  LR = qprobeval3D_GTH_mm512<OUT_TYPE,OUT_OLD>(A, m, MD, Jmin, Jmax, M, N, LogScale, ScaleI, Xm, Y, I, K, PMik, Yik, Gmin, Gmin2, J, Hmin, Hmin2, delX, pcontig);
#else
	  LR = qprobeval3D_GTH<OUT_TYPE,OUT_OLD>(A, m, MD, Jmin, Jmax, M, N, LogScale, ScaleI, Xm, Y, I, K, PMik, Yik, Gmin, Gmin2, J, Hmin, Hmin2, delX, pcontig);
#endif
	  if(DEBUG/* HERE >=2 */ && !(!IsNan(LR) && LR >= 0.0)){
#pragma omp critical
	    {
	      printf("m=%d:I=%d,K=%d,J=%d,G=%d..%d:LR= %0.10e, LogScale[I]= %0.10f, isfinite(LR)=%d\n",m,I,K,J,Gmin,I-K-1,LR, LogScale[I],isfinite(LR) ? 1 : 0);
	      fflush(stdout);
	      assert(!IsNan(LR) && LR >= 0.0);
	    }
	  }

	  if(DEBUG>=2 && J <= 1 && !(LR == A.UL(I,K,J))) {
	    #pragma omp critical
	    {
	      printf("m=%d:I=%d,K=%d,J=%d,G=%d..%d: LR= %0.8e, A.UL(I,K,J)= %0.8e\n",m,I,K,J,Gmin,I-K-1,LR, A.UL(I,K,J));
	      fflush(stdout);
	      assert(LR == A.UL(I,K,J));
	    }
	  }

	  if(!SCALE_MONOTONIC && LR >= MaxFloat){/* early re-normalization for I */
	    double rescale = 0.99999 * MaxLR/LR;
	    if(DEBUG) assert(rescale < 1.0);
	    for(int k= 0; k <= Kmax[I]; k++)
	      for(int j = Jmin[I]; j <= Jmax[I]; j++){
		A.LR(I,k,j) *= rescale;
		A.bLR(I,k,j) *= rescale;
		A.UL(I,k,j) *= rescale;
	      }
	    LogScale[I] -= log(rescale);
	    if(DEBUG/* HERE >=2 */ && !isfinite(LogScale[I])){
	      printf("m=%d,I=%d,K=%d,J=%d:LR= %0.8e, rescale= %0.8e, LogScale[I]= %0.8e\n",m,I,K,J,LR,rescale,LogScale[I]);
	      fflush(stdout);
	      assert(isfinite(LogScale[I]));
	    }
	    for(int g = max(IMIN,I - A.kmax - deltaExtYRef); g < I; g++){
	      ScaleI[g] *= rescale;
	      if(DEBUG) assert(isfinite(ScaleI[g]) && ScaleI[g] <= DBL_MAX);
	    }

	    maxLR *= rescale;
	    Lscale = exp(-LogScale[I]);// NEW
	    if(DEBUG) assert(isfinite(Lscale) && Lscale <= 1.0);
	  
	    /* need to redo G,T,H loop */
#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT
	    LR = qprobeval3D_GTH_mm256<OUT_TYPE,OUT_OLD>(A, m, MD, Jmin, Jmax, M, N, LogScale, ScaleI, Xm, Y, I, K, PMik, Yik, Gmin, Gmin2, J, Hmin, Hmin2, delX);
#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT
	    LR = qprobeval3D_GTH_mm512<OUT_TYPE,OUT_OLD>(A, m, MD, Jmin, Jmax, M, N, LogScale, ScaleI, Xm, Y, I, K, PMik, Yik, Gmin, Gmin2, J, Hmin, Hmin2, delX, pcontig);
#else
	    LR = qprobeval3D_GTH<OUT_TYPE,OUT_OLD>(A, m, MD, Jmin, Jmax, M, N, LogScale, ScaleI, Xm, Y, I, K, PMik, Yik, Gmin, Gmin2, J, Hmin, Hmin2, delX, pcontig);
#endif

#if DEBUG>=2
	    if(DEBUG && !(isfinite(LR) && LR < MaxFloat)){
#pragma omp critical
	      {
		printf("m=%d:I=%d,K=%d,J=%d,G=%d..%d:LR= %0.10e, PMik= %0.10f, LogScale[I]= %0.10f, rescale= %0.10e\n",m,I,K,J,Gmin,I-K-1,LR, PMik,LogScale[I],rescale);
		fflush(stdout);
		assert(isfinite(LR) && LR < MaxFloat);
	      }
	    }
#endif
	  }

#if DEBUG>=2
	  if(DEBUG && !(isfinite(LR) && LR < MaxFloat)){
#pragma omp critical
	    {
	      printf("m=%d(mapid=%d,id=%lld):I=%d,K=%d,J=%d,G=%d..%d:LR= %0.10e, PMik= %0.10f, LogScale[I]= %0.10f, IsNan(LR)=%d\n",
		     m,pcontig->contig[m].mapid, gmap[pcontig->contig[m].mapid]->id,I,K,J,Gmin,I-K-1,LR, PMik,LogScale[I],IsNan(LR) ? 1 : 0);
	      fflush(stdout);
	      assert(isfinite(LR) && LR < MaxFloat);
	    }
	  }
#endif

#if DEBUG>=2
	  if(DEBUG && !(isfinite(LR) && LR < MaxFloat)){
	    printf("m=%d:I=%d,K=%d,J=%d,G=%d..%d:LR= %0.10e, PMik= %0.10f, LogScale[I]= %0.10f\n",m,I,K,J,Gmin,I-K-1,LR, PMik,LogScale[I]);
	    fflush(stdout);
	    assert(isfinite(LR) && LR < MaxFloat);
	  }
#endif

	  if((DEBUG && !(isfinite(LR) && LR < MaxFloat)) || (LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0) || HDEBUG) LRtrace3(A,m,MD,I,K,J,Gmin,N,Kmax,LR,LogScale);
	} // J >= 2

	if(DEBUG>=2 && J <= 1) assert(LR== A.UL(I,K,J));

	if(VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && I == 14997 && K == 0 && J == 138){
	  printf("m=%d:I=%d,K=%d,J=%d,G=%d..%d:A.AL(I,K,J)= %10e -> %0.10e, UL= %0.10e, Lscale= %0.10e, LR= %0.10e, maxLR= %0.10e, LogScale[I]= %0.10f\n",
		 m,I,K,J,Gmin,I-K-1,A.LR(I,K,J), LR*PMik, A.UL(I,K,J), Lscale, LR, max((double)maxLR, LR), LogScale[I]);
	  fflush(stdout);
	}

	MFLOAT1 nLR = A.LR(I,K,J) = LR * PMik;
	maxLR = max(maxLR, nLR);

	if(VERB>=3 && m==38){
	  printf("m=%d,I=%d,K=%d,J=%d:A.LR(I,K,J)= %0.10e, maxLR= %0.10e, V3D_rcnt=%lu/%lu\n",m,I,K,J,A.LR(I,K,J),maxLR, V3D_rcnt, V3D_totcnt);
	  fflush(stdout);
	}
      } // J = Jmin[I] .. Jmax[I]
    } // K = 0 .. Kmax[I]

    if(DEBUG>=2){
      for(int K = 0; K <= Kmax[I]; K++){
	for(int J= Jmin[I]; J <= Jmax[I]; J++){
	  MFLOAT origLR = A.LR(I,K,J);
	  if(DEBUG && !(origLR <= maxLR)){
	    #pragma omp critical
	    {
	      printf("m=%d,I=%d:Kmax[I]=%d,Jmin[I]=%d,Jmax[I]=%d,maxLR=%0.10e:K=%d,J=%d,A.LR(I,K,J)= %0.10e, V3D_rcnt=%lu/%lu\n",
		     m,I,Kmax[I],Jmin[I],Jmax[I],maxLR,K,J,A.LR(I,K,J), V3D_rcnt, V3D_totcnt);
	      fflush(stdout); 
	      assert(origLR <= maxLR);
	    }
	  }
	}
      }
    }

    /* maxLR may be 0.0 if all UL values were 0 and no previous alignments are present (eg I=1,J=7 with DELTA=6) */
    if(DEBUG) assert(isfinite(maxLR) && maxLR >= ((MFLOAT)0.0) && maxLR < MaxFloat);

    if(!SCALE_MONOTONIC || maxLR > MaxLR){ /* need to rescale A[I][*][*].LR so largest value (maxLR) is <= MaxLR and update LogScale[I] and Scale[G=IMIN..I-1] */
      MFLOAT1 rescale = ((maxLR > 1.0 || maxLR > (MaxLR/MaxFloat)) ? max(1.00001f * MinFloat, 0.99999f * MaxLR/maxLR) : 0.99999f * MaxFloat);/* multiplier for A.LR(I,*,*)  and ScaleI[1..I-1] */
      if(DEBUG && !(MinFloat < rescale && rescale < MaxFloat)){
	#pragma omp critical
	{
	  printf("\nm=%d:I=%d,maxLR=%0.10e,MaxLR=%0.4e,LogScale[I]=%0.10f,rescale=%0.10e,MinFloat=%0.4e,MaxFloat=%0.4e, MaxLR/maxLR=%0.10e\n", 
		 m,I,maxLR,MaxLR,LogScale[I],rescale,MinFloat,MaxFloat, (MFLOAT1)(MaxLR/maxLR));
	  fflush(stdout);
	  assert(MinFloat < rescale && rescale < MaxFloat);
	}
      }
      double MaxScaleI = DBL_MAX;
      if(!SCALE_MONOTONIC && rescale > 1.0){
	for(int G= max(IMIN, I - A.kmax - deltaExtYRef); G < I; G++){/* avoid overflow of ScaleLI[] */
	  if(rescale <= 1.0)
	    break;
	  if(ScaleI[G] >= MaxScaleI/rescale)
	    rescale = 0.99999 * MaxScaleI/ScaleI[G];
	}
	if(rescale > 1.0){/* make sure LogScale[I] does not drop below 0.0 */
	  if(LogScale[I] - log(rescale) < 0.0)
	    rescale = exp(LogScale[I]);
	}
	if(DEBUG) assert(MinFloat < rescale && rescale < MaxFloat);
      }
      if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792*/){
	printf("m=%d:I=%d,maxLR=%0.10e,MaxLR=%0.4e,rescale=%0.10e,LogScale[I]=%0.10f (log(maxLR)+LogScale[I]=%0.10f)\n", 
	       m,I,maxLR,MaxLR,rescale,LogScale[I], (maxLR > 0.0 ? log(maxLR) : -1000.0) + LogScale[I]);
	fflush(stdout);
      }
      for(int K= 0; K <= Kmax[I]; K++)
	for(int J= Jmin[I]; J <= Jmax[I]; J++){
	  MFLOAT origLR = A.LR(I,K,J);
	  if(DEBUG>=2) assert(origLR <= maxLR);
	  A.LR(I,K,J) = origLR * rescale;
	  if(DEBUG>=2) assert(isfinite(A.LR(I,K,J)) && A.LR(I,K,J) <= MaxLR);
	  A.bLR(I,K,J) *= rescale;
	  A.UL(I,K,J) *= rescale;
	}
      double origLogScale = LogScale[I];
      LogScale[I] -= log(rescale);
      if(DEBUG) assert(isfinite(LogScale[I]) && LogScale[I] >= -1e-5);
      LogScale[I] = max(0.0, LogScale[I]);
      if(DEBUG>=2) assert(isfinite(LogScale[I]) && LogScale[I] >= 0.0);
      for(int G= max(IMIN, I - A.kmax - deltaExtYRef); G < I; G++){
	ScaleI[G] *= rescale;
	if(DEBUG) assert(isfinite(ScaleI[G]) && ScaleI[G] <= (SCALE_MONOTONIC ? 1.0 : MaxScaleI));
      }
      if(VERB>=3 && rverb && m==MDEBUG_M){
	printf("I=%d,LogScale[I]= %0.10f -> %0.10f, maxLR= %0.10e -> %0.10e, A.kmax=%d\n",I,origLogScale,LogScale[I],maxLR,maxLR * rescale, A.kmax);
	fflush(stdout);
      }
      //      maxLR *= rescale;// redundant since maxLR is about to go out of scope
    } // if(maxLR > MAXLR)

    if(VERB>=3 && rverb && m==MDEBUG_M){
      for(int K= 0; K <= Kmax[I]; K++)
	for(int J= Jmin[I]; J <= Jmax[I]; J++)
          printf("\t m=%d:I=%d,K=%d,J=%d:AL(I,K,J)= %0.10e, UL(I,K,J)= %0.10e, LogScale[I]= %0.10f\n",m,I,K,J,A.LR(I,K,J),A.UL(I,K,J), LogScale[I]);
      fflush(stdout);
    }
  } // I = IMIN .. IMAX

  LogScale[IMAX+1] = LogScale[IMAX];
  if(DEBUG) assert(isfinite(LogScale[IMAX+1]));
}

/** return log likelihood (based on additive recurrance) and new best alignment of each Map X[k] to consensus Y[] */
__attribute__ ((noinline))
double qprobeval(int n, double *Hcuts,/**< Hcuts[i=0..n+1] : complete consensus map */
		 int N, double *Y,/**< Y[I=0..N+1] : current consensus map hypothesis */
		 int MD, /**< number of Maps */
		 int *MX, double **X, /**< X[m=0..MD-1][j=0..MX[m]+1] : Map sites */
		 int lc,/**< If left end of Y is a linear chromosome end */
		 int rc,/**< If right end of Y is a linear chromosome end */
		 int **map, /**< map[m=0..MD-1][j=0..MX[m]+1] is index in Hcuts[0..n+1] of X[m][j], 
					       map[MD][j=0..n+1] is index in Y[0..N+1] of Hcuts[j] */
		 int **mapK, /**< mapK[m=0..MD-1][j=0..MX[m]+1] is K index offset in Hcuts[0..n+1] of X[m][j] */ 
		 Csetlimit *limit, /* setlimit() cache : a function of map[], mapK[] and Hcuts[] and slightly dependent on nmap[MD],nmapK[MD] */
		 int **nmap, /**< If nmapK ! = 0 : nmap[m=0..MD-1][j=0..MX[m]+1] will be new index in Y[0..N+1] of X[m][j], 
					 nmap[MD][I=0..N+1] is index in Hcuts[0..n+1] of Y[I] and always defined */
		 int **nmapK, /**< If nmapK != 0 : nmapK[m=0..M-1][j=1..MX[m]] will be new K index offset in Y[0..N+1] of X[m][j] 
					  However nmapK[0] == Lij, mapK[MX[m]+1] == Rij (-ve if ends are outliers, unless OUTLIER_FIX >= 1) */
		 double *TBmapWT, /**< If != 0 : TBmapWT[m] is a the -TB based weight (multiplied by pcontig->mapWT[m] due to -Erefine or -BestRefWT) 
				     for map m to be applied when summing up the total log(LR+LRbias) value */
		 double *newLP,/**< return value of Log Likelihood for each Map */
		 double *mapWT,/**< return value of Map weights = LR/(LR + LRbias) (mapWT[m] is NOT updated IF (newLP[m] is copied from oldLP[m]) OR (nmapK, LogPV, poutlier AND TBmapWT[m] are all 0)) */
		 double *oldLP,/**< If != 0 : previous estimate of Log Likelihood for each Map */
		 double Ymin, double Ymax, /**< location range in Y[] that changed since oldLP was computed : only valid if oldLP != 0 */
		 double *LogPV, /**< If != 0 : log10(Pvalue)  of each map m will be returned in LogPV[m] */
		 int **poutlier, /**< If != 0 : outlier[m][j=0..MX[m]] returns outlier flag for interval X[m][j..j+1] */
		 Ccontig *pcontig, /**< point to complete contig information for debugging */
		 int *forcemap /* If != 0 : forcemap[m=0..MD-1] forces/disables use of oldLP : 0 means force use of oldLP[m], 1 means never use oldLP[m], -1 means use setlimit() to decide and record decision in forcemap[m] */
		 )
{
  START(qprobeval);

  qprobevalcnt++;

  DELTA_YU = FIX_U ? deltaExtYRef + FIX_U - 1 : DELTA_Y;// NEW46
  DELTA_XU = FIX_U ? deltaExtXRef + FIX_U - 1 : DELTA_X;// NEW46
  if(DEBUG) assert(DELTA_YU >= DELTA_Y);
  if(DEBUG) assert(DELTA_XU >= DELTA_X);

  if(VERB>=2){
    printf("qprobeval:start: cum wall time= %0.6f:qprobevalcnt=%d\n", wtime(),qprobevalcnt);
    fflush(stdout);
  }

  /* check that overflow and underflow is not likely to happen : see refine.cpp definitions of MaxLR etc for details */
  if(DEBUG) assert(MaxLR * MaxFA < MaxFloat);// to avoid overflow during 3-D recurrance 
  if(DEBUG) assert(MaxLR * MaxUR * MaxScale < DBL_MAX);// to avoid overflow during LRcum summation
  if(DEBUG && SCALE_MONOTONIC) assert((MaxLR * MinLR) * MinUR > DBL_MIN);// to avoid underflow during LRcum summation
  if(DEBUG && !SCALE_MONOTONIC) assert((MaxFloat * MinLR) * MinUR > DBL_MIN);// to avoid underflow during LRcum summation
  if(DEBUG) assert(MinUR  > MinFloat);
  if(DEBUG) assert(LRCUM_MAX > DBL_MIN);// to avoid underflow during 1st LRcum rescaling with LS
  if(DEBUG) assert(MaxLR * MaxScale > 1.0);// to avoid underflow during 2nd LRcum rescaling with LS
  if(DEBUG) assert(MaxLR >= 1.0);
  if(DEBUG) assert(MaxLR * MaxUR < DBL_MAX);
  if(PoutlierEnd < MinUR){
    printf("qprobeval : PoutlierEnd = %0.3e is less than MinUR = %0.3e. To avoid floating point underflow please raise -endoutlierRef and -endoutlierSwitch values to at least MinUR\n",
	   PoutlierEnd, MinUR);
    fflush(stdout);exit(1);
  }

  double LLRsum = 0.0;

  if(DEBUG) assert(map != 0);
  if(DEBUG) assert(mapK != 0);
  if(DEBUG) assert(nmap != 0);

  if(!ISLINEAR && (rc || lc)){
    printf("qprobeval:rc=%d,lc=%d not supported (set ISLINEAR in constants.h)\n",rc,lc);
    exit(1);
  }

  if(DEBUG) assert(colors==1);
  if(DEBUG) assert(!END_VITERBI);
  if(DEBUG) assert(resSD[0] > 0.0);
  if(DEBUG) assert(LRbias >= 0.0);

  /* NOTE : If oldLP, not all values of nmapK will be updated 
  if(oldLP && nmapK){
    printf("qprobeval: Cannot use both oldLP and update nmap\n");
    fflush(stdout);
    assert(0);
    }*/ 

  if(VERB>=2 && rverb && (TRACE_ID > 0 || MDEBUG_M >= 0)){
    for(int i = 0; i <= N+1; i++)
      printf("Y[%d]=%0.10f\n",i,Y[i]);
    fflush(stdout);
  }

  if(DEBUG>=2){  /* make sure no real cuts are closer than minKB */
    for(int I = 1; I < N; I++)
      if(DEBUG && !(Y[I+1]-Y[I] >= minKB - 2e-9)){
	printf("\nqprobeval:I=%d:Y[I]=%0.10f,Y[I+1]=%0.10f, delta=%0.10f, minKB=%0.6f\n",
	       I,Y[I],Y[I+1], Y[I+1]-Y[I], minKB);
	for(int i = 1; i <= N+1; i++)
	  printf("Y[%d]=%0.10f\n",i,Y[i]);
	fflush(stdout);
	assert(Y[I+1]-Y[I] >= minKB - 2e-9);
      }
  }
  if(DEBUG>=2){/* make sure map[] and mapK[] are consistent */
    for(int m = 0; m < MD; m++){
      int M = MX[m];
      int lastJ = -1, lasti = -1, lastk = -1, i, k;
      for(int J = 1; J <= M; J++){
	if((i=map[m][J]) > 0){
	  k = mapK[m][J];
	  assert(0 <= k && k < i);
	  if(lasti >= 0 && !(i-k > lasti)){
	    printf("m=%d,J=%d,i=%d,k=%d:lastJ=%d,lasti=%d,lastk=%d\n",m,J,i,k,lastJ,lasti,lastk);
	    fflush(stdout);
	    assert(i-k > lasti);
	  }
	  lasti = i;
	  lastk = k;
	  lastJ = J;
	}
      }
    }
  }

  /* Note: I-Kmax[I] must be monotonically non-decreasing functions of I (unless Kmax[I] < 0) */
  if(DEBUG>=2){
    int Tmax = 1;
    for(int I= 1; I <= N; I++){
      assert(I-Kmax[I] >= 1);
      assert(I-Kmax[I] >= Tmax);
      if(I-Kmax[I] > Tmax)
	Tmax = I-Kmax[I];
    }
  }

  int maxMX = 0;
  double AvMX = 0.0;
  for(int m = 0; m < MD; m++){
    maxMX = max(maxMX,MX[m]);
    AvMX += MX[m];
  }
  if(MD > 0)
    AvMX /= MD;
  if(DEBUG) assert(maxMX <= (int)MAXINT);

  long long double_pool_size = (4LL * N + 5LL + (maxMX+1) * 2LL)*sizeof(double);/* 6 allocations */
  long long int_pool_size = (3LL * (N+1) + 5LL * (maxMX+1))*sizeof(int);/* 8 allocations */

  long long IKcnt = 0;
  for(int I=1; I <= N; I++)
    IKcnt += Kmax[I] + 1;

  kmax = 0;
  double AvKmax = 0.0;
  for(int I = 1; I <= N; I++){
    kmax = max(kmax,Kmax[I]);
    AvKmax += Kmax[I];
  }
  AvKmax /= N;

  if(DEBUG && CACHE && ((long long)N) * ((long long)maxMX) > (long long)MASK_LL(31)){
    printf("WARNING: N=%d,maxMX=%d : N*maxMX=%llu is too large to fit in 32 bit int\n",
	   N,maxMX,((long long)N) * ((long long)maxMX));
    fflush(stdout);exit(1);
  }
  int KstrideMax = N * maxMX;
  if(CACHE > 0 && ALIGN > 0)
    KstrideMax = (KstrideMax + ALIGN - 1) & (~(ALIGN-1));
  long long StrideMax = (1LL + kmax) * KstrideMax + VPAD;
  if(CACHE > 0 && ALIGN2 > 0)
    StrideMax = (StrideMax + ALIGN2 - 1) & (~(ALIGN2-1));

  /* HERE : If DIAGONAL && MIN_MEM, Fmem[] and Imem[] are not allocated via heap->alloc() but instead using malloc/free to minimize memory usage
     For predicted virtual memory usage for Fmem[] and Imem[] replace StrideMax with  (1 + kmax) * (IMAX - IMIN + 1) * deltaJ ~ (1 + kmax) * (maxMX + 2*RANGE) * maxMX
     For predicted real memory usage for Fmem[] and Imem[] replace IKcnt * maxMX with min(IKcnt, (1 + av(Kmax[]) * (av(MX[]) + 2 * RANGE)) * av(MX[]) */

  long long strideEstimate = (DIAGONAL && MIN_MEM) ? min(StrideMax, (1LL + kmax) * (maxMX + 2LL * RANGE) * maxMX) : StrideMax;
  long long pool_Fmem_siz = strideEstimate * DPsizF * sizeof(MFLOAT1);
  long long pool_Imem_siz = strideEstimate * DPsizI * sizeof(int);

  long long pool_Vec_siz = 0;
  int num_allocations = 16;

  long long virtperthread_size = pool_Vec_siz + pool_Fmem_siz + pool_Imem_siz + double_pool_size + int_pool_size + /* alignment fudge */num_allocations*(4*CACHE_LINE_GUARD);
  long long perthread_size = (long long)((1 + AvKmax) * (AvMX + 2 * RANGE) * AvMX) /* WAS IKcnt * maxMX */ 
    * (DPsizF * sizeof(MFLOAT1) + DPsizI * sizeof(int)) + double_pool_size + int_pool_size + /* alignment fudge */num_allocations*(4*CACHE_LINE_GUARD);

  if(DEBUG) assert(perthread_size >= 0);
  double virtmemperthread = (double)(virtperthread_size + OMP_stack_size);
  double memperthread = (double)perthread_size;
  double maxmem = (MaxMem > 0.0) ? MaxMem * 1024.0 * 1024.0 * 1024.0 : 1e+13;
  double maxvirtmem = (MaxVirtMem > 0.0) ? MaxVirtMem * 1024.0 * 1024.0 * 1024.0 : 1e+13;
  
  /* determine number of threads (numthreads) to use */
  int numthreads = 1;
#ifdef _OPENMP
  numthreads = MaxThreads; // WAS omp_get_max_threads();
  if(RefineThreads > 0){
    if(BoostThreads && RefineThreads < BoostThreads && wtime() > BoostThreadsTime){
      int origRefineThreads = RefineThreads;
      MaxThreads = omp_get_max_threads();
      numthreads = MaxThreads = RefineThreads = BoostThreads = min(MaxThreads, BoostThreads);
      if(TotalThreads)
	MaxVirtMem *= ((double) RefineThreads)/origRefineThreads;
      if(VERB){
	printf("Increasing number of threads for Refinment from %d to %d (BoostThreads= %d)(wall time=%0.6f secs)\n",
	       min(MaxThreads,origRefineThreads),min(numthreads,RefineThreads), BoostThreads, wtime());
	fflush(stdout);
      }
    }
    numthreads = min(numthreads,RefineThreads);
  }
  if(MaxMem > 0 && perthread_size > 0){
    delete heap;    /* free lightweight heap : this should release virtual memory space with glibc */
    heap = new lightweight_heap(0,0);      
    if(TotalThreads){
      //      getmem(VmSize,VmRSS, VmSwap);
      maxvirtmem -= min(VmSize, VmRSS + VmSwap);
      // maxmem -= VmRSS;// NOTE : this memory can be allowed to swap out
    }
    int limitthreads = max(1,(int)floor(min(maxmem/memperthread, maxvirtmem/virtmemperthread) + 0.5));
    if(DEBUG) assert(memperthread > 0.0 && maxmem > 0.0 && limitthreads >= 1);
    limitthreads = min(numthreads,limitthreads);
    if(limitthreads < numthreads){
      if(TotalThreads){/* update VmSize,VmRSS,VmSwap and recompute */
	maxvirtmem = (MaxVirtMem > 0.0) ? MaxVirtMem * 1024.0 * 1024.0 * 1024.0 : 1e+40;
	getmem(VmSize,VmRSS, VmSwap);
	maxvirtmem -= min(VmSize, VmRSS + VmSwap);
	limitthreads = max(1,(int)floor(min(maxmem/memperthread, maxvirtmem/virtmemperthread) + 0.5));
	if(DEBUG) assert(memperthread > 0.0 && maxmem > 0.0 && limitthreads >= 1);
	limitthreads = min(numthreads,limitthreads);
      }
      if(VERB && (VERB/* HERE >=2 */ || DEBUG>=2 || verb_numthreads)){
	double InvGB = 1e-9;
	printf("threads reduced from %d to %d due to MaxVirtMem=%0.4f(%0.4f) Gb, MaxMem=%0.4f(%0.4f) Gb, per thread virt= %0.4f Gb, real= %0.4f Gb\n\t VmRSS=%0.4f,VmSwap=%0.4f,VmSize=%0.4f Gb,Stack=%0.2f Mb:N=%d,IKcnt=%lld,kmax=%d,maxMX=%d, AvKmax=%0.3f,AvMX=%0.1f,RANGE=%d,pools:double=%0.4f,int=%0.4f Gb\n",
	       numthreads,limitthreads, MaxVirtMem, maxvirtmem*InvGB, MaxMem, maxmem*InvGB, virtmemperthread*InvGB, memperthread*InvGB,VmRSS*InvGB, VmSwap*InvGB,VmSize*InvGB, OMP_stack_size*InvGB*1024,
	       N, IKcnt, kmax,maxMX,AvKmax,AvMX,RANGE,double_pool_size*InvGB, int_pool_size*InvGB);
	fflush(stdout);
      }
      if(DEBUG) assert(limitthreads >= 1);
      numthreads = limitthreads;
      extern int mverb_numthreads;
      mverb_numthreads = verb_numthreads = 1;
    }
  }
  if(numthreads > MD)
    numthreads = max(1,MD);
  if(VERB && (VERB>=2 || verb_numthreads)){
    double InvGB = 1e-9;
    printf("using %d threads:MaxVirtMem=%0.4f(%0.4f) Gb, MaxMem=%0.4f(%0.4f) Gb, per thread virt= %0.4f Gb, real= %0.4f Gb\n\t VmRSS=%0.4f,VmSwap=%0.4f,VmSize=%0.4f Gb,Stack=%0.2f Mb:N=%d,IKcnt=%lld,kmax=%d,maxMX=%d, AvKmax=%0.3f,AvMX=%0.1f,RANGE=%d,pools:double=%0.4f,int=%0.4f Gb\n",
	   numthreads, MaxVirtMem, maxvirtmem*InvGB, MaxMem, maxmem*InvGB, virtmemperthread*InvGB, memperthread*InvGB,VmRSS*InvGB, VmSwap*InvGB,VmSize*InvGB, OMP_stack_size*InvGB*1024.0,
	   N, IKcnt, kmax,maxMX,AvKmax,AvMX,RANGE,double_pool_size*InvGB, int_pool_size*InvGB);
    if(VERB/* HERE >=2 */)
      printf("\t StrideMax=%lld,DPsizF=%lld,DPsizI=%lld:Fmem=%lld,Imem=%lld,fudge=%d(%d):total per thread=%lld,%lld\n",
 	   StrideMax,DPsizF,DPsizI,pool_Fmem_siz,pool_Imem_siz, num_allocations*(4*CACHE_LINE_GUARD),num_allocations, virtperthread_size,perthread_size);
    fflush(stdout);
    verb_numthreads = 0;
  }
  if(DEBUG) assert(numthreads >= 1);

  #endif // _OPENMP

  if(DIAGONAL && MIN_MEM){/* Fmem[] and Imem[] will be allocated inside each thread */
    virtperthread_size -= pool_Fmem_siz + pool_Imem_siz;
    pool_Fmem_siz = pool_Imem_siz = 0;
  }

  heap->reset_lightweight_heap(virtperthread_size, numthreads);
    
  if((DEBUG && !(Ymax >= Ymin && Ymax <= Y[N+1])) || (VERB>=2)){
    printf("qprobeval:N=%d:Ymin=%0.6f,Ymax=%0.6f,Y[N+1]=%0.6f:cum wall time= %0.6f\n",
	   N,Ymin,Ymax,Y[N+1], wtime());
    fflush(stdout);
    assert(Ymax>=Ymin && Ymax <= Y[N+1]);
  }

  int skipcnt;/* how often the alignment was skipped */
  int setcnt;/* how often setlimit was called */
  int qsetcnt;/* how often setlimit() used SETLIMIT_CACHE */
  int globalcnt;/* how often setlimit returned a global alignment limit (IMIN==1,IMAX==N) */
  if(TIME_VERB || CNT_VERB)
    skipcnt = setcnt = qsetcnt = globalcnt = 0;
  
  bool TwoLoops = (TWO_LOOPS && SETLIMIT_CACHE && (oldLP || (!nmapK && !LogPV && !poutlier)));

  int *mlist = NULL;
  int nindex = MD;
  if(TWO_LOOPS && TwoLoops){
    mlist = new int[MD];
    nindex = 0;
  }

  double cumwt1,cumwt2,cumwt3,cumwt4,cumwt5,cumwt5b,cumwt6,cumwt7,cumwt8;
  if(TIME_VERB)
    cumwt1 = cumwt2 = cumwt3 = cumwt4 = cumwt5 = cumwt5b = cumwt6 = cumwt7 = cumwt8 = 0.0;

  int LoopCnt = 0;
  START(qprobeval_parallel_wait);
  #pragma omp parallel num_threads(numthreads) if (rverb==0 && numthreads > 1)
  {
    START(qprobeval_parallel);
    int myLoopCnt = 0;

    int myskipcnt;
    int mysetcnt;
    int myqsetcnt;
    int myglobalcnt;
    if(CNT_VERB || TIME_VERB)
      myskipcnt = mysetcnt = myqsetcnt = myglobalcnt = 0;

    double wt, wt0;
    if(TIME_VERB)
      wt0 = wt = wtime();

    double mywt1,mywt2,mywt3,mywt4,mywt5,mywt5b,mywt6,mywt7,mywt8;
    if(TIME_VERB)
      mywt1 = mywt2 = mywt3 = mywt4 = mywt5 = mywt5b = mywt6 = mywt7 = mywt8 = 0.0;

    /* Obtain thread id */
    int tid=0;
    #ifdef _OPENMP
    tid = omp_get_thread_num ();
    #endif

    /* allocate all heap memory that is local to each thread */
    QA A;
    const int OUT_TYPE = OUTLIER_TYPE;
    const int OUT_OLD = OUTLIER_OLD(outlierMaxBC);

    long long myStrideMax = 0;// Used with MIN_MEM : currently Allocated memory can handle stride values up to myStrideMax
    int myKstrideMax = 0;// Used with CACHE
    MFLOAT1 *Fmem = NULL;
    int *Imem = NULL;
    if(!(DIAGONAL && MIN_MEM)){ /* allocate Fmem[],Imem[] */
      if(sizeof(MFLOAT1) == sizeof(MFLOAT))
        Fmem = (MFLOAT1 *)heap->alloc_MFLOAT(StrideMax * DPsizF, tid);
      else
	Fmem = (MFLOAT1 *)heap->alloc_double(StrideMax * DPsizF, tid);
      Imem = heap->alloc_int(StrideMax * DPsizI, tid);
      myStrideMax = StrideMax;
      myKstrideMax = KstrideMax;
#ifdef VALGRIND
      VALGRIND_MAKE_MEM_UNDEFINED(Fmem,StrideMax * DPsizF * sizeof(MFLOAT1));
      VALGRIND_MAKE_MEM_UNDEFINED(Imem,StrideMax * DPsizI * sizeof(int));
#endif
    }

    double *LogScale = heap->alloc_double(N+2, tid);
    double *ScaleI = heap->alloc_double(N+1, tid); /* ScaleI[G=1..I] = exp(LogScale[G]-logScale[I]) */
    double *JXmin = heap->alloc_double(N+1, tid);
    double *JXmax = heap->alloc_double(N+1, tid);

    int IMIN = 1, IMAX = N;
    int *Jmin = heap->alloc_int(N+1, tid);
    int *Jmax = heap->alloc_int(N+1, tid);
    int *Imin = heap->alloc_int(maxMX+1, tid); 
    int *Imax = heap->alloc_int(maxMX+1, tid); 
    int IMIN1,IMAX1;/* value of IMIN,IMAX before extrapolating from the current alignment */

    int *Jmid = heap->alloc_int(N+1, tid);
    
    int *Ilist = (int *)heap->alloc_int(maxMX+1,tid);
    int *Jlist = (int *)heap->alloc_int(maxMX+1,tid);
    int *Klist = (int *)heap->alloc_int(maxMX+1,tid);
    FLOAT *Xfrag = (FLOAT *)heap->alloc_double(maxMX+1,tid);
    FLOAT *Yfrag = (FLOAT *)heap->alloc_double(maxMX+1,tid);

    START(qprobeval_for_wait);

    if(TWO_LOOPS && TwoLoops){/* first do a fast static loop that don't need setlimit() or alignment() and save the remaining loop m values in mlist[0..nindex-1] */

      #pragma omp for schedule(static,8)
      for(int m = 0; m < MD; m++) {
        START(qprobeval_for);
	if(OMP_DEBUG) myLoopCnt++;
	int M = MX[m];
	double *Xm = X[m];
	double TBwtm = TBmapWT ? TBmapWT[m] : 1.0;
	if(DEBUG>=1+RELEASE/* HERE >=2 */ && TBmapWT && !isfinite(TBmapWT[m])){
          printf("qprobeval:m=%d/%d:TBmapWT[m]= %0.8e\n",m,MD,TBmapWT[m]);
	  fflush(stdout);
	  assert(isfinite(TBmapWT[m]));
        }

	/* set default values for LogPV[m] and poutlier[m][0..M] */
	if(LogPV)
	  LogPV[m] = 0.0;
	if(poutlier)
	  for(int j = 0; j <= M; j++)
	    poutlier[m][j] = 0;

	int noskip = 1;

	if(VERB/* HERE >=2 */ && SETLIMIT_CACHE && oldLP && Ymin >= 0.0 && (VERB>=2 || (rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID  || m == MDEBUG_M)))){
          Csetlimit *plimit = &limit[m];
	  if(plimit->ileft >= 0){
            double Yleft = Hcuts[plimit->ileft] - plimit->Xleft;
	    double Yright = Hcuts[plimit->iright] + plimit->Xright;

	    printf("m=%d/%d:oldLP[m]= %0.6f,TBwtm=%0.6e, Ymin=%0.4f,Ymax=%0.6f,Yleft=%0.4f,Yright=%0.4f (Hcuts[ileft=%d,iright=%d]=%0.4f,%0.4f,Xleft=%0.4f,Xright=%0.4f):IMIN=%d,IMAX=%d,overlap=%d,noskip=%d,nmapK=%d,LogPV=%d,poutlier=%d,forcemap[m]=%d\n",
           	m,MD,oldLP[m],TBwtm,Ymin,Ymax,Yleft,Yright,plimit->ileft,plimit->iright,Hcuts[plimit->ileft],Hcuts[plimit->iright],plimit->Xleft,plimit->Xright,IMIN,IMAX,
	        (max(Yleft,Ymin) > min(Yright,Ymax) + 1e-6) ? 0 : 1,noskip,nmapK ? 1 : 0, LogPV ? 1 : 0, poutlier ? 1 : 0, forcemap ? forcemap[m] : -2);
          } else
	    printf("m=%d/%d:oldLP[m]= %0.6f,TBwtm=%0.6e, Ymin=%0.4f,Ymax=%0.6f, ileft=%d:IMIN=%d,IMAX=%d,noskip=%d,nmapK=%d,LogPV=%d,poutlier=%d,forcemap[m]=%d\n",
	         m,MD,oldLP[m],TBwtm,Ymin,Ymax,plimit->ileft,IMIN,IMAX,noskip,nmapK ? 1 : 0, LogPV ? 1 : 0, poutlier ? 1 : 0, forcemap ? forcemap[m] : -2);
	  fflush(stdout);
        }

	/* check if this map can be skipped by re-using previous alignment score oldLP[m] (or 0.0 if TBwtm <= 0.0) */
	if(oldLP || (!nmapK && !LogPV && !poutlier && TBwtm <= 0.0)){/* HERE HERE : If nmapK != 0 it should still be OK to skip based on TBwtm <= 0.0 and rely on previous values of map/mapK being correct, provided UpdateMap() also skips updating map/mapK for maps with TBwtm <= 0.0 AND reposition() keeps map/mapK updated for all maps */

	  if(TBwtm > 0.0){
	    int force = 0;
	    if(SETLIMIT_CACHE && Ymin >= 0.0 && limit[m].ileft >= 0){
	      Csetlimit *plimit = &limit[m];
	      double Yleft = Hcuts[plimit->ileft] - plimit->Xleft;
	      double Yright = Hcuts[plimit->iright] + plimit->Xright;
	      if(max(Yleft,Ymin) > min(Yright,Ymax) + 1e-6){
		if(CNT_VERB || TIME_VERB)
		  myqsetcnt++;
		noskip = 0;
	      } else if(SETLIMIT_CACHE_FORCE)
		force = 1;
	    }

	    if(forcemap && forcemap[m] >= 0){
	      if(forcemap[m] == 1)
		noskip = force = 1;// forces use of alignment
	      else
		noskip = force = 0;// forces use of oldLP
	    }

	    if(noskip){
	      int v;

#pragma omp atomic capture
	      v = nindex++;

	      mlist[v] = m;
	      continue;// postpone setlimit()
	    }
	  }// if(TBwtm > 0.0)

	  if(forcemap){
	    int newforce = (TBwtm <= 0.0 || !noskip) ? 0 : 1; 

	    if(DEBUG && forcemap[m] >= 0) assert(newforce == forcemap[m]);

	    if(forcemap[m] < 0)
	      forcemap[m] = newforce;
	  }

	  /* Using cache to hold last estimate of IMIN,IMAX,Y[IMIN-LRANGE],Y[IMAX+LRANGE] that is reset only when map[],mapK[] are updated with UpdateMap  */
          if(TBwtm <= 0.0 || !noskip) { /* Xm does not overlap part of Y[] (Ymin..Ymax) changed since oldLP[] was computed */
	    if(DEBUG>=2 && TBwtm > 0.0) assert(oldLP);
	    newLP[m] = oldLP ? oldLP[m] : 0.0;

	    if(VERB && (VERB>=2 || QDEBUG || rverb || (MDEBUG_M >= 0 && m==MDEBUG_M && MTRACE && N==MDEBUG_N))){
              #pragma omp critical
	      {
	        LLRsum += newLP[m] * TBwtm;
	      
		if(VERB && (VERB>=2 || rverb>=3 || (MDEBUG_M >= 0 && m==MDEBUG_M && MTRACE && N==MDEBUG_N))){
	          if(TBwtm > 0.0){
	            double Yleft = Y[max(0,IMIN-LRANGE)] - Xm[Jmax[IMIN]];
		    double Yright = Y[min(N+1,IMAX+LRANGE)] + Xm[M+1] - Xm[Jmin[IMAX]];
		    printf("m=%d/%d:mapid=%d(id=%lld):LLRsum=%0.10f(LLR=%0.10f): Xm overlaps Y at %0.3f..%0.3f, Ychange=%0.3f..%0.3f(IMIN=%d,IMAX=%d,LRANGE=%d,X[Jmax[IMIN]=%d]=%03f,X[Jmin[IMAX]=%d]=%0.3f,X[M+1]=%0.3f)\n",
		             m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,LLRsum,newLP[m],Yleft,Yright,Ymin,Ymax,IMIN,IMAX,
		             LRANGE,Jmax[IMIN],Xm[Jmax[IMIN]],Jmin[IMAX],Xm[Jmin[IMAX]],Xm[M+1]);
	          } else
	            printf("m=%d/%d:mapid=%d(id=%lld):LLRsum=%0.10f(LLR=%0.10f): TBwtm= %0.6e (skipping)\n",
		         m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,LLRsum,newLP[m],TBwtm);
		  fflush(stdout);
	        }
	        if(DEBUG) assert(isfinite(LLRsum));
	      }
	    }

	    if(TIME_VERB || CNT_VERB) myskipcnt++;

	    if(DEBUG>=2 && !(isfinite(newLP[m]) && isfinite(TBwtm))){
	      printf("qprobeval:m=%d/%d:newLP[m]= %0.8e, TBwtm= %0.8e\n",m,MD,newLP[m],TBwtm);
	      fflush(stdout);
	      assert(isfinite(newLP[m]));
	      assert(isfinite(TBwtm));
	    }

	    STOP(qprobeval_for);
	    continue;
	  } // if(TBwtm <= 0.0 || !noskip)
        } 

        int v;

        #pragma omp atomic capture
        v = nindex++;

        mlist[v] = m;
        continue;// postpone setlimit & alignment
        STOP(qprobeval_for);
      } // omp for(int m = 0; m < MD; m++)
    } // if(TwoLoops)

    #pragma omp for schedule(dynamic,1)
    for(int v = 0; v < nindex; v++) {
      START(qprobeval_for);
      int m = (TWO_LOOPS && TwoLoops) ? mlist[v] : v;
      int M = MX[m];
      double *Xm = X[m];
      double TBwtm = TBmapWT ? TBmapWT[m] : 1.0;
      if(DEBUG>=1+RELEASE/* HERE >=2 */ && TBmapWT && !isfinite(TBmapWT[m])){
	printf("qprobeval:m=%d/%d:TBmapWT[m]= %0.8e\n",m,MD,TBmapWT[m]);
	fflush(stdout);
	assert(isfinite(TBmapWT[m]));
      }

      double wt1;
      if(TIME_VERB)
        wt1 = wtime();	       

      /* set default values for LogPV[m] and poutlier[m][0..M] */
      if(LogPV)
	LogPV[m] = 0.0;
      if(poutlier)
	for(int j = 0; j <= M; j++)
	  poutlier[m][j] = 0;

      int noskip = 1;

      if(VERB/* HERE >=2 */ && SETLIMIT_CACHE && oldLP && Ymin >= 0.0 && (VERB>=2 || (rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID  || m == MDEBUG_M)))){
	Csetlimit *plimit = &limit[m];
	if(plimit->ileft >= 0){
	  double Yleft = Hcuts[plimit->ileft] - plimit->Xleft;
	  double Yright = Hcuts[plimit->iright] + plimit->Xright;

	  printf("m=%d/%d:oldLP[m]= %0.6f,TBwtm=%0.6e, Ymin=%0.4f,Ymax=%0.6f,Yleft=%0.4f,Yright=%0.4f (Hcuts[ileft=%d,iright=%d]=%0.4f,%0.4f,Xleft=%0.4f,Xright=%0.4f):IMIN=%d,IMAX=%d,overlap=%d,noskip=%d,nmapK=%d,LogPV=%d,poutlier=%d,forcemap[m]=%d\n",
		 m,MD,oldLP[m],TBwtm,Ymin,Ymax,Yleft,Yright,plimit->ileft,plimit->iright,Hcuts[plimit->ileft],Hcuts[plimit->iright],plimit->Xleft,plimit->Xright,IMIN,IMAX,(max(Yleft,Ymin) > min(Yright,Ymax) + 1e-6) ? 0 : 1,noskip,
		 nmapK ? 1 : 0, LogPV ? 1 : 0, poutlier ? 1 : 0, forcemap ? forcemap[m] : -2);
	} else
	  printf("m=%d/%d:oldLP[m]= %0.6f,TBwtm=%0.6e, Ymin=%0.4f,Ymax=%0.6f, ileft=%d:IMIN=%d,IMAX=%d,noskip=%d,nmapK=%d,LogPV=%d,poutlier=%d,forcemap[m]=%d\n",
		 m,MD,oldLP[m],TBwtm,Ymin,Ymax,plimit->ileft,IMIN,IMAX,noskip,nmapK ? 1 : 0, LogPV ? 1 : 0, poutlier ? 1 : 0, forcemap ? forcemap[m] : -2);
	fflush(stdout);
      }

      /* check if this map can be skipped by re-using previous alignment score oldLP[m] (or 0.0 if TBwtm <= 0.0) */
      if(oldLP || (!nmapK && !LogPV && !poutlier && TBwtm <= 0.0)){/* HERE HERE : If nmapK != 0 it should still be OK to skip based on TBwtm <= 0.0 and rely on previous values of map/mapK being correct, provided UpdateMap() also skips updating map/mapK for maps with TBwtm <= 0.0 AND reposition() keeps map/mapK updated for all maps */

	bool setlimit_called = false;

	if(TBwtm > 0.0){
	  int force = 0;
	  if(SETLIMIT_CACHE && Ymin >= 0.0 && limit[m].ileft >= 0){
	    Csetlimit *plimit = &limit[m];
	    double Yleft = Hcuts[plimit->ileft] - plimit->Xleft;
	    double Yright = Hcuts[plimit->iright] + plimit->Xright;
	    if(max(Yleft,Ymin) > min(Yright,Ymax) + 1e-6){
	      if(CNT_VERB || TIME_VERB) 
		myqsetcnt++;
	      noskip = 0;
	    } else if(SETLIMIT_CACHE_FORCE)
	      force = 1;
	  }

	  if(forcemap && forcemap[m] >= 0){
	    if(forcemap[m] == 1)
	      noskip = force = 1;// forces use of alignment
	    else
	      noskip = force = 0;// forces use of oldLP
	  }

	  if(noskip){
	    if(CNT_VERB || TIME_VERB)
	      mysetcnt++;
	    setlimit_called = true;

	    noskip = setlimit(m,n,Hcuts,N,Y,M,Xm,map[m],mapK[m],&limit[m],map[MD],nmap[MD],Ymin,Ymax,IMIN,IMAX,Jmin,Jmax,Imin,Imax,JXmin,JXmax,Jmid,IMIN1,IMAX1,pcontig,force);
	  }
	}// if(TBwtm > 0.0)

	if(forcemap){
	  int newforce = (TBwtm <= 0.0 || !noskip) ? 0 : 1; 

	  if(DEBUG && forcemap[m] >= 0) assert(newforce == forcemap[m]);

	  if(forcemap[m] < 0)
	    forcemap[m] = newforce;
	}

	/* Using cache to hold last estimate of IMIN,IMAX,Y[IMIN-LRANGE],Y[IMAX+LRANGE] that is reset only when map[],mapK[] are updated with UpdateMap  */
	if(TBwtm <= 0.0 || !noskip) {
          if(DEBUG>=1+RELEASE && TWO_LOOPS && TwoLoops && !setlimit_called){
	    printf("v=%d/%d:m=%d: this iteration should have been skipped previously: TBwtm= %0.6f, noskip=%d,qprobevalcnt=%d\n",v,nindex,m,TBwtm,noskip,qprobevalcnt);
	    fflush(stdout); 
	    assert(!(TBwtm <= 0.0 || !noskip));
          }

          if(DEBUG>=2 && TBwtm > 0.0) assert(oldLP);
	  newLP[m] = oldLP ? oldLP[m] : 0.0;

	  if(VERB && (VERB>=2 || QDEBUG || rverb || (MDEBUG_M >= 0 && m==MDEBUG_M && MTRACE && N==MDEBUG_N))){
            #pragma omp critical
	    {
	      LLRsum += newLP[m] * TBwtm;
	      
	      if(VERB && (VERB>=2 || rverb>=3 || (MDEBUG_M >= 0 && m==MDEBUG_M && MTRACE && N==MDEBUG_N))){
		if(TBwtm > 0.0){
		  double Yleft = Y[max(0,IMIN-LRANGE)] - Xm[Jmax[IMIN]];
		  double Yright = Y[min(N+1,IMAX+LRANGE)] + Xm[M+1] - Xm[Jmin[IMAX]];
		  printf("m=%d/%d:mapid=%d(id=%lld):LLRsum=%0.10f(LLR=%0.10f): Xm overlaps Y at %0.3f..%0.3f, Ychange=%0.3f..%0.3f(IMIN=%d,IMAX=%d,LRANGE=%d,X[Jmax[IMIN]=%d]=%03f,X[Jmin[IMAX]=%d]=%0.3f,X[M+1]=%0.3f)\n",
			 m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,LLRsum,newLP[m],Yleft,Yright,Ymin,Ymax,IMIN,IMAX,
			 LRANGE,Jmax[IMIN],Xm[Jmax[IMIN]],Jmin[IMAX],Xm[Jmin[IMAX]],Xm[M+1]);
		} else
		  printf("m=%d/%d:mapid=%d(id=%lld):LLRsum=%0.10f(LLR=%0.10f): TBwtm= %0.6e (skipping)\n",
		         m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,LLRsum,newLP[m],TBwtm);
		fflush(stdout);
	      }
	      if(DEBUG) assert(isfinite(LLRsum));
	    }
	  }

	  if(TIME_VERB || CNT_VERB) myskipcnt++;

	  if((DEBUG>=2 && !(isfinite(newLP[m]) && isfinite(TBwtm))) || (VERB>=3 && rverb)){
	    printf("qprobeval:m=%d/%d:newLP[m]= %0.8e, TBwtm= %0.8e\n",m,MD,newLP[m],TBwtm);
	    fflush(stdout);
	    assert(isfinite(newLP[m]));
	    assert(isfinite(TBwtm));
	  }
	  
	  STOP(qprobeval_for);
	  continue;
        }
      } else {       /* Limit DP to near the previous alignment defined by Xm[J] to Y[I=map[MD][map[m][J]]] (I= -1 if Xm[J] is not aligned) */
	if(CNT_VERB || TIME_VERB)
	  mysetcnt++;

	if(DEBUG && forcemap){
	  assert(forcemap[m] != 0);
	  forcemap[m] = 1;
	}

	(void) setlimit(m,
			n,
			Hcuts,
			N,
			Y,
			M,
			Xm,
			map[m],
			mapK[m],
			&limit[m],
			map[MD],
			nmap[MD],
			-1.0,
			-1.0,
			IMIN,
			IMAX,
			Jmin,
			Jmax,
			Imin,
			Imax,
			JXmin,
			JXmax,
			Jmid,
			IMIN1,
			IMAX1,
			pcontig,
			1);/* sets IMIN,IMAX,Jmin[IMIN..IMAX],Jmax[IMIN..IMAX] */
      }

      if(CNT_VERB || TIME_VERB)
	myglobalcnt += (N > 0 && IMIN==1 && IMAX==N) ? 1 : 0;     

      if(DEBUG && forcemap) assert(forcemap[m] >= 0);
      
      if(VERB>=3 && (N > 0 && IMIN==1 && IMAX==N)){
	printf("\t m=%d/%d:setlimit() returned global alignment, IMIN=%d,IMAX=%d,N=%d!\n",m,MD,IMIN,IMAX,N);
	fflush(stdout);
	assert(!(N > 0 && IMIN==1 && IMAX==N));
      }

      if(VERB/* HERE >=2 */ && SETLIMIT_CACHE && oldLP && Ymin >= 0.0 && (VERB>=2 || (rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID  || m == MDEBUG_M)))){
	Csetlimit *plimit = &limit[m];
	if(plimit->ileft >= 0){
	  double Yleft = Hcuts[plimit->ileft] - plimit->Xleft;
	  double Yright = Hcuts[plimit->iright] + plimit->Xright;

	  printf("m=%d/%d:oldLP[m]= %0.6f,TBwtm=%0.6e, Ymin=%0.4f,Ymax=%0.6f,Yleft=%0.4f,Yright=%0.4f (Hcuts[ileft=%d,iright=%d]=%0.4f,%0.4f,Xleft=%0.4f,Xright=%0.4f):IMIN=%d,IMAX=%d,overlap=%d,noskip=%d,nmapK=%d,LogPV=%d,poutlier=%d\n",
		 m,MD,oldLP[m],TBwtm,Ymin,Ymax,Yleft,Yright,plimit->ileft,plimit->iright,Hcuts[plimit->ileft],Hcuts[plimit->iright],plimit->Xleft,plimit->Xright,IMIN,IMAX,(max(Yleft,Ymin) > min(Yright,Ymax) + 1e-6) ? 0 : 1,noskip,
		 nmapK ? 1 : 0, LogPV ? 1 : 0, poutlier ? 1 : 0);
	} else
	  printf("m=%d/%d:oldLP[m]= %0.6f,TBwtm=%0.6e, Ymin=%0.4f,Ymax=%0.6f, ileft=%d:IMIN=%d,IMAX=%d,noskip=%d,nmapK=%d,LogPV=%d,poutlier=%d\n",
		 m,MD,oldLP[m],TBwtm,Ymin,Ymax,plimit->ileft,IMIN,IMAX,noskip,nmapK ? 1 : 0, LogPV ? 1 : 0, poutlier ? 1 : 0);
	fflush(stdout);
      }

      if(IMIN > IMAX){/* setlimit found no valid previous alignment and globafallback==0 */
	if(DEBUG) assert(!globalfallback);

	double LRcum = max(1e-300,LRbias) + PoutlierEnd * PoutlierEnd;
	newLP[m] = log(LRcum);
	mapWT[m] = PoutlierEnd * PoutlierEnd / LRcum;// NEW2
	if(DEBUG>=1+RELEASE/* HERE HERE >=2 */ && !(isfinite(newLP[m]) && isfinite(mapWT[m]))){
#pragma omp critical
	  {
	    printf("qprobeval:m=%d/%d:mapid=%d(id=%lld):newLP[m]= %0.8e, mapWT[m]= %0.8e, LRbias= %0.8e, PoutlierEnd= %0.8e,IMIN=%d,IMAX=%d,N=%d,M=%d\n",
		   m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,newLP[m],mapWT[m],LRbias,PoutlierEnd,IMIN,IMAX,N,M);
	    fflush(stdout);
	    assert(isfinite(newLP[m]));
	    assert(isfinite(mapWT[m]));
	  }
	}

	if(VERB && (VERB>=2 || QDEBUG || rverb || (MDEBUG_M >= 0 && m==MDEBUG_M && MTRACE && N==MDEBUG_N))){
#pragma omp critical
	  {
	    LLRsum += newLP[m] * TBwtm;
	   
	    if(VERB && (VERB>=2 || rverb>=3 || (MDEBUG && TBwtm > 0.0) || (MDEBUG_M >= 0 && m==MDEBUG_M && MTRACE && N==MDEBUG_N))){
	      printf("m=%d/%d:mapid=%d(id=%lld):LLRsum=%0.10f(LLR=%0.10f, TBmapWT[m]= %0.10f): setlimit() found no valid previous alignment\n",
		     m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,LLRsum,newLP[m], TBwtm);
	      fflush(stdout);
	    }
	   
	    if(DEBUG) assert(isfinite(LLRsum));
	  }
	}

	if(nmapK)
	  for(int J = 0; J <= M+1; J++)
	    nmap[m][J] = -1;
       
	if(TIME_VERB || CNT_VERB) myskipcnt++;
       
	STOP(qprobeval_for);
	continue;
      } // if(IMIN > IMAX)

      if((VERB>=2 && tverb>=2 && (m==728 || m==9244 || tverb >= 3)) /* || (DEBUG>=1+RELEASE && IMAX-IMIN >= 10000)*/){
#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After setlimit():id=%lld:M=%d,Xm[M+1]= %0.3f kb, Y[IMIN=%d,IMAX=%d,N=%d]= %0.3f, %0.3f, %0.3f:wall time= %0.6f\n", 
		 m, tid, gmap[pcontig->contig[m].mapid]->id, M, Xm[M+1], IMIN,IMAX, N, Y[IMIN], Y[IMAX],Y[N], wtime() - wt);
	  if(stage >= 2){
	    for(int J = 1; J <= M; J++){
	      int i = map[m][J];
	      if(i < 0)
		continue;
	      int k = mapK[m][J];
	      int I = map[MD][i];
	      int K = map[MD][i-k];
	      printf("\t\t\t J=%d:map[m][J]=i=%d, mapK[m][J]=k= %d, n=%d : map[MD][i]=I=%d, map[MD][i-k]=K=%d: Xm[J]=%0.3f: Hcuts[i-k,i]= %0.3f, %0.3f, Y[K,I]= %0.3f,%0.3f\n", 
		     J,i,k,n,I,K,Xm[J],Hcuts[i-k], Hcuts[i], Y[max(0,K)],Y[max(0,I)]);
	    }
	    for(int I = IMIN; I <= IMAX; I++)
	      printf("\t\t\t I=%d:Kmax[I]=%d,Jmin[I]=%d,Jmax[I]=%d:Y[I]= %0.3f, Xm[Jmin[I]]= %0.3f, Xm[Jmax[I]]= %0.3f\n",I,Kmax[I],Jmin[I],Jmax[I],Y[I],Xm[Jmin[I]], Xm[Jmax[I]]);
	  }
	  fflush(stdout);
	  if(DEBUG>=1+RELEASE) assert(IMAX-IMIN < 10000);
	}
      }

      double wt2;
      if(TIME_VERB)
        wt2 = wtime();	       

      if(MDEBUG_M >= 0 && rverb && m==MDEBUG_M){
	printf("m=%d/%d:global=%d,IMIN=%d,IMAX=%d:\n",m,MD,global,IMIN,IMAX);
	for(int I=IMIN; I <= IMAX; I++)
	  printf("  I=%d:Kmax[I]=%d,Jmin[I]=%d,Jmax[I]=%d\n",I,Kmax[I],Jmin[I],Jmax[I]);
	fflush(stdout);
      }

      int NeedBP = (nmapK || ((VERB && (VERB>=2 || rverb || (MTRACE && MDEBUG_M>=0 && N==MDEBUG_N && m==MDEBUG_M))) || LogPV || poutlier)) ? 1 : 0;/* If we need best alignment(bLR,G,T,H,Rij) or bestLR */

      /* re-allocate memory for A.field(I= IMIN..IMAX, K=0..Kmax[I], J= Jmin[I]..Jmax[I]) */
      A.kmax = Kmax[IMIN];
      for(int I = IMIN+1; I <= IMAX; I++)
	A.kmax = max(A.kmax,Kmax[I]);
      if(DEBUG) assert(A.kmax <= kmax);

      long long origStrideMax = myStrideMax;
      int origKstrideMax = myKstrideMax;
      long long stride, firstindex;
      if(DIAGONAL){
	// reduce real memory usage by reallocating array within each thread to take advantage of IMIN,IMAX,Jmin[],Jmax[] (set by setlimit()) : c.f. refalign() hashdelta based memory allocation
	int Jmin1 = Jmin[IMIN];
	int Jmax1 = Jmax[IMIN];
	for(int I = IMIN+1; I <= IMAX; I++){
	  Jmin1 = min(Jmin1,Jmin[I]-(I-IMIN));
	  Jmax1 = max(Jmax1,Jmax[I]-(I-IMIN));
	}
	/* NOTE: Allocated J range for each I=IMIN..IMAX is : Jmin1+(I-IMIN) .. Jmax1+(I-IMIN) and involves a slight overallocation if there are indels or label density differences in the alignment */
	int deltaJ = Jmax1 - Jmin1 + 1;
	if(deltaJ >= M){/* revert to non-diagonal method : will still use less memory since only I = IMIN..IMAX, J = 1..M is allocated */
	  deltaJ = M;
	  Jmin1 = 1;
	  Jmax1 = M;
	  A.Istride = M;
	} else
	  A.Istride = deltaJ - 1;

	if(DEBUG) assert(IMAX - IMIN + 1LL <= N);

	if(DEBUG && sizeof(A.Istride) <= 4 && ((long long)A.Istride * (long long)IMAX > (long long)MASK_LL(31))){// NEW7
	  printf("ERROR in qprobeval.cpp: IMIN=%d,IMAX=%d,N=%d,Jmin1=%d,Jmax1=%d,deltaJ=%d,M=%d,Istride=%d : Istride * IMAX = %llu is too large to fit in 32 bit int\n",
		 IMIN,IMAX,N,Jmin1,Jmax1,deltaJ,M,A.Istride, ((long long)A.Istride) * ((long long)IMAX));
	  fflush(stdout);exit(1);
	}

	if(DEBUG && sizeof(A.Kstride) <= 4 && ((long long)(IMAX-IMIN + 1LL)) * ((long long)deltaJ) > (long long)MASK_LL(31)){
	  printf("ERROR in qprobeval.cpp: IMIN=%d,IMAX=%d,N=%d,Jmin1=%d,Jmax1=%d,deltaJ=%d,M=%d : (IMAX-IMIN+1)*deltaJ=%llu is too large to fit in 32 bit int\n",
		 IMIN,IMAX,N,Jmin1,Jmax1,deltaJ,M,((long long)(IMAX-IMIN + 1LL)) * ((long long)deltaJ));
	  fflush(stdout);exit(1);
	}

	A.Kstride = (IMAX - IMIN + 1LL) * deltaJ;

	if(DEBUG && sizeof(A.Kstride) <= 4 && ((long long)A.Kstride * (long long)A.kmax > (long long)MASK_LL(31))){// NEW7
	  printf("ERROR in qprobeval.cpp: IMIN=%d,IMAX=%d,N=%d,Jmin1=%d,Jmax1=%d,deltaJ=%d,M=%d,Istride=%d,Kstride=%d,kmax=%d : Kstride * kmax = %llu is too large to fit in 32 bit int\n",
		 IMIN,IMAX,N,Jmin1,Jmax1,deltaJ,M,A.Istride,A.Kstride,A.kmax, ((long long)A.Kstride) * ((long long)A.kmax));
	  fflush(stdout);exit(1);
	}

	if(CACHE > 0 && ALIGN > 0)
	  A.Kstride = (A.Kstride + ALIGN - 1) & (~(ALIGN-1));
	stride = (1LL + (CACHE ? kmax : A.kmax)) * A.Kstride;
	if(CACHE > 0 && ALIGN2 > 0)
	  stride = (stride + ALIGN2 - 1) & (~(ALIGN2-1));
	firstindex = IMIN * A.Istride + Jmin1;

	if(DEBUG && CACHE && (ALIGN > 0 || ALIGN2 > 0) && stride <= myStrideMax && !(A.Kstride <= myKstrideMax)){
	  #pragma omp critical
	  {
	    printf("m=%d:IMIN=%d,IMAX=%d,N=%d,deltaJ=%d:A.kmax=%d(kmax=%d),A.Kstride=%d(myKstrideMax=%d),stride=%lld(myStrideMax=%lld)\n",
		   m,IMIN,IMAX,N,deltaJ,A.kmax,kmax,A.Kstride,myKstrideMax,stride,myStrideMax);
	    fflush(stdout);
	    assert(A.Kstride <= myKstrideMax);
	  }
	}

	int MemBins = max(2,MIN_MEM);/* Allow for some extra memory of 1 part per MemBins and reduce memory when need drops below (MemBins-1)/(MemBins+1) */

	if(MIN_MEM && (stride > myStrideMax || MIN_MEM >= 8 || (MIN_MEM>=2 && oldLP&&/*NEW*/ stride * (MemBins+1) < myStrideMax * (MemBins-1)))){

	  /* minimize virtual memory and swapping by reallocating Fmem[],Imem[] using either free()/malloc() or mmap()/madvise(DONT_NEED) */
#if USE_MMAP
	  if(Fmem) free_huge_pages(Fmem);
#else
	  if(Fmem) free(Fmem);
#endif

	  if(DEBUG) assert(stride <= StrideMax);
	  if(DEBUG && CACHE) assert(A.Kstride <= KstrideMax);

          myKstrideMax = min(A.Kstride * (MIN_MEM >= 8 ? MemBins : MemBins+1) / MemBins, KstrideMax);
	  if(CACHE > 0 && ALIGN > 0){
	    int origmyKstrideMax = myKstrideMax;
	    myKstrideMax = (myKstrideMax + ALIGN - 1) & (~(ALIGN-1));
	    if(DEBUG>=2 && !(myKstrideMax >= origmyKstrideMax)){
	      #pragma omp critical
	      {
		printf("origmyKstrideMax=0x%x(%d),ALIGN=0x%x(%d),myKstrideMax=0x%x(%d)\n",
		       origmyKstrideMax,origmyKstrideMax,ALIGN,ALIGN,myKstrideMax,myKstrideMax);
		fflush(stdout);
		assert(myKstrideMax >= origmyKstrideMax);
	      }
	    }
	  }
	  if(DEBUG) assert(A.Kstride <= myKstrideMax && myKstrideMax <= KstrideMax);

	  myStrideMax = (1LL + (CACHE ? kmax : A.kmax)) * myKstrideMax;
	  if(CACHE > 0 && ALIGN2 > 0){
	    long long origmyStrideMax = myStrideMax;
	    myStrideMax = (myStrideMax + ALIGN2 - 1) & (~(ALIGN2 - 1));
	    if(DEBUG>=2) assert(myStrideMax >= origmyStrideMax);
	  }
	  if(DEBUG) assert(stride <= myStrideMax && myStrideMax <= StrideMax);

	  if(CACHE > 0 && ALIGN2 > 0){/* increase myKstrideMax to largest possible value <= myStrideMax / (1LL + kmax) */
	    int origmyKstrideMax = myKstrideMax;
	    myKstrideMax = (myStrideMax / (1LL + kmax)) & (~(ALIGN - 1));
	    myKstrideMax = min(myKstrideMax, KstrideMax);
	    if(DEBUG>=2) assert(myKstrideMax >= origmyKstrideMax);
	    if(DEBUG>=2){
	      long long newStrideMax = (1LL + kmax) * myKstrideMax;
	      newStrideMax = (newStrideMax + ALIGN2 - 1) & (~(ALIGN2 - 1));	      
	      assert(newStrideMax == myStrideMax);
	    }
	  }

	  size_t memsiz = myStrideMax * (DPsizF * sizeof(MFLOAT1) + DPsizI*sizeof(int)) + PADDING;

#if USE_MMAP
	  Fmem = (MFLOAT1 *) malloc_huge_pages(memsiz);
#else
	  if(CACHE && (ALIGN > 0 || ALIGN2 > 0))
	    POSIX_MEMALIGN(Fmem, max(ALIGN,ALIGN2)*sizeof(int), memsiz);
	  else
	    Fmem = (MFLOAT1 *)malloc(memsiz);
#endif
	  if(!Fmem){
	    printf("qprobeval:tid=%d:malloc(%llu) failed : stride=%lld,myStrideMax=%lld(StrideMax=%lld),DPsizF=%lld,DPsizI=%lld\n",
		   tid,(unsigned long long)memsiz,  stride,myStrideMax,StrideMax,DPsizF,DPsizI);
	    fflush(stdout); assert(Fmem);
	  }
	  Imem = (int *) &Fmem[DPsizF * myStrideMax];
#ifdef VALGRIND
	  VALGRIND_MAKE_MEM_UNDEFINED(Fmem,memsiz);
#endif

#if VECTOR_AVX && USE_AVX && !USE_AVX512 && USE_MFLOAT
	  memset(Fmem, 0, memsiz);// make sure there are no NaN or Inf values lurking in Padded regions of memory : HERE HERE reduce overhead by only initializing unused regions
#endif
	}
	if(CACHE){
	  if(DEBUG) assert(A.kmax <= kmax);
	  if(DEBUG) assert(stride <= myStrideMax);
	  if(DEBUG) assert(A.Kstride <= myKstrideMax);
	  A.Kstride = myKstrideMax;
	  stride = myStrideMax;
	  if(DEBUG && !ALIGN2) assert(stride == (1LL + (CACHE ? kmax : A.kmax)) * A.Kstride);
	  if(DEBUG && ALIGN2 > 0) assert(stride >= (1LL + (CACHE ? kmax : A.kmax)) * A.Kstride);
	}
      } else {// DIAGONAL == 0
	if(DEBUG && ((long long)N) * ((long long)M) > (long long)MASK_LL(31)){
	  printf("WARNING: N=%d,M=%d : N*M=%llu is too large to fit in 32 bit int\n",N,M,((long long)N)*((long long)M));
	  fflush(stdout);exit(1);
	}
	A.Kstride = N * M;
	A.Istride = M;
	stride = (1LL + A.kmax) * A.Kstride;
	firstindex = A.Istride + 1LL;
      }

      if(DEBUG) assert(stride <= StrideMax);
      if(DEBUG && CACHE && DIAGONAL){
	assert(A.Kstride == myKstrideMax);
	assert(myKstrideMax <= KstrideMax);
	if(myStrideMax == origStrideMax)
	  assert(myKstrideMax == origKstrideMax);
	assert(stride == myStrideMax);
	if(!ALIGN2) assert(stride == (1LL + (CACHE ? kmax : A.kmax)) * A.Kstride);
	if(ALIGN2) assert(stride >= (1LL + (CACHE ? kmax : A.kmax)) * A.Kstride);
      }

      if(DEBUG) assert(stride <= myStrideMax);

      long long mcnt = 0;
#define FBLOCK(field) { A.field##_base = &Fmem[mcnt - firstindex]; mcnt += stride;}
      FBLOCK(LR);
      FBLOCK(bLR);
      FBLOCK(UL);
#if URDEBUG
      FBLOCK(UR);
#endif
#undef FBLOCK
      if(DEBUG) assert(mcnt <= stride * DPsizF);
      
      mcnt = 0;
#define BLOCK(field) { A.field##_base = &Imem[mcnt - firstindex]; mcnt += stride;}
      BLOCK(G);
      BLOCK(T);
      BLOCK(H);
      BLOCK(Lij);
      BLOCK(Rij);
#undef BLOCK
      if(DEBUG) assert(mcnt <= stride * DPsizI);

      if(LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M){
	#pragma omp critical
	{
	  printf("m=%d:IMIN=%d,IMAX=%d\n",m,IMIN,IMAX);
	  fflush(stdout);
	}
      }

      double wt3;
      if(TIME_VERB)
        wt3 = wtime();	       

      /* initialize A[I= IMIN..IMAX][K=0..Kmax[I]][J= Jmin[I]..Jmax[I]].LR  = PoutlierEnd */
      MFLOAT PoutlierEndf = PoutlierEnd;
      for(int I= IMIN; I <= IMAX; I++){
	if(LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M && I==MDEBUG_I){
	  #pragma omp critical
	  {
	    printf("m=%d:I=%d,Kmax[I]==%d,Jmin[I]=%d,Jmax[I]=%d\n",m,I,Kmax[I],Jmin[I],Jmax[I]);
	    fflush(stdout);
	  }
	}
	for(int K = 0; K <= Kmax[I]; K++){
	  int jmin = Jmin[I];
	  int jmax = Jmax[I];
	  MFLOAT1 LRend = PoutlierEndf;// WAS * PRtab[K][I].PM;
	  MFLOAT1 *pLR = &A.LR(I,K,0);
	  MFLOAT1 *pbLR = &A.bLR(I,K,0);
	  int *pG = &A.G(I,K,0);
	  int *pLij = &A.Lij(I,K,0);
	  int *pRij = &A.Rij(I,K,0);

	  for(int J= jmin; J <= jmax; J++){
	    if((LTRACE && MDEBUG_I >= 0 && MDEBUG_M>=0 && N==MDEBUG_N3 && m==MDEBUG_M /* && A.LR(I,K,J) >= MDEBUG_MIN */ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)
	       || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && I == 2 && K == 0 && J == 1)){
	      #pragma omp critical
	      {
		printf("m=%d:I=%d,K=%d,J=%d:UL=%0.10e(PMik=%0.10e,PoutlierEnd=%0.10e)\n",m,I,K,J,LRend,PRtab[K][I].PM,PoutlierEnd);
		fflush(stdout);
	      }
	    }
	    pLR[J] = pbLR[J] = LRend; /* should never be the best score if PoutlierEnd == 0 */
	    pG[J] = -2;/* outlier (local alignment) termination */
	    pLij[J] = I-K;
	    pRij[J] = I;
	  }
	}
      }

      double wt4;
      if(TIME_VERB)
        wt4 = wtime();	       

      /* initialize A[I=Imin[J]..Imax[J]][K=0..Kmax[I]][J=1..M].Lij as leftmost site in Y[] that overlaps Xm with Y[I-K..I] aligned with Xm[J] */
      /* initialize corresponding LR score with left ends unaligned */
      /* only cases where J <= DELTA_XU (OR EXTEND && I-K <= DELTA_YU && J <= DELTA_XU + EXTEND_MAX) are initialized  */
      int jmax = EXTEND ? min(M, DELTA_XU + EXTEND_MAX) : min(M, DELTA_XU);
      if((VERB>=2 && MDEBUG_M >= 0 && N==MDEBUG_N && m==MDEBUG_M)
	 || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M))){
	printf("m=%d:jmax=%d,M=%d,DELTA_XU=%d\n",m,jmax,M,DELTA_XU);
	fflush(stdout);
      }

      for(int J = 1; J <= jmax; J++){
	int Lij = 0;
	double IvarX;
	double FLj = FLJ(J,Xm[J],IvarX);// HERE : try MFLOAT FLJf()
	int ImKmax = (!EXTEND || J <= DELTA_XU) ? Imax[J] : min(DELTA_YU, Imax[J]);
	if((VERB>=2 && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M && J==MDEBUG_J)
	   || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && J == 138)){
	  printf("m=%d,J=%d:ImKmax=%d,Imin[J]=%d,Imax[J]=%d,DELTA_XU=%d,DELTA_YU=%d\n",m,J,ImKmax,Imin[J],Imax[J],DELTA_XU,DELTA_YU);
	  fflush(stdout);
	}
	for(int ImK = max(1,Imin[J] - KMAX); ImK <= ImKmax; ImK++){
	  int Lijk = Lij;
	  int I;
	  if((VERB>=2 && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M && J==MDEBUG_J && ImK == (MDEBUG_I-MDEBUG_K))
	     || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && ImK == 14997 && J == 138)){
	    printf("m=%d,J=%d,ImK=%d:Imin[J]=%d,Imax[J]=%d\n",m,J,ImK,Imin[J],Imax[J]);
	    fflush(stdout);
	  }
	  for(int K = max(0,Imin[J]-ImK); (I=ImK+K) <= Imax[J] && K <= Kmax[I]; K++){
	    if((VERB>=2 && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M && J==MDEBUG_J && I == MDEBUG_I && K == MDEBUG_K)
	       || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && I == 14997 && K == 0 && J == 138)){
	      printf("m=%d,J=%d,I=%d,K=%d:Jmin[I]=%d,Jmax[I]=%d,Kmax[I]=%d\n",m,J,I,K,Jmin[I],Jmax[I],Kmax[I]);
	      fflush(stdout);
	    }
	    if(DEBUG>=2) assert(!(J < Jmin[I] || J > Jmax[I]));
	    if(DEBUG>=2) assert(J <= DELTA_XU || (EXTEND && I-K <= DELTA_YU && J <= DELTA_XU + EXTEND_MAX));
	    if(DEBUG>=2) assert(I >= Imin[J]);

	    double y = Yc(Y,I,K);
	    MFLOAT PMik = PRtab[K][I].PM;

	    /* update Lijk by incrementing it until Y[Lijk] overlaps Xm */
	    while(y-Y[Lijk] > Xm[J])
	      Lijk++;
	    if(DEBUG>=2) assert(Lijk <= I && !(y-Y[Lijk] > Xm[J]) && (Lijk==0 || y-Y[Lijk-1] > Xm[J]));
	    A.Lij(I,K,J) = min(Lijk,I-K);
	    if(K == 0)
	      Lij = Lijk;
	  
	    double NotOutlier = FLIK_FIX ? PoutlierEndM1 : FLj * PoutlierEndM1;// WAS * PMik;

	    MFLOAT1 UL,bUL;
	    int Umin,Umax,Ub;
	    Ub = Umin = Umax = max(1,A.Lij(I,K,J)) - 1;
	    if(VERB>=3 && rverb && m==MDEBUG_M && I==33 && K==0 && J==4)
	      rverb = 2;
	    bUL = UL = FLIK<0>(I,K,J,Umin,Umin+1,y,Xm[J],y-Y[Umin],Y[Umin+1]-Y[Umin],Y,IvarX, FLj, END_FIX ? lc : 1);
	    if(VERB>=3 && rverb>=2 && m==MDEBUG_M && I==33 && K==0 && J==4)
	      rverb = 1;
	    if(DEBUG>=2){
	      MFLOAT FLikj = FL<0>(I,K,J,Umin,Umin+1,y,Xm[J],y-Y[Umin],Y[Umin+1]-Y[Umin],Y, END_FIX ? lc : 1);
	      if(DEBUG && !(fabs(FLikj - (FLIK_FIX ? UL : UL * FLj)) <= FLikj * (USE_MFLOAT ? 1e-4 : 2e-9) + (USE_MFLOAT ? 1e-32 : 1e-305))){
		#pragma omp critical
		{
		  printf("m=%d,I=%d,K=%d,J=%d,U=%d,N=%d:FLik=%0.8e,FLj=%0.8e,FLikj=%0.8e,FLik*FLj=%0.8e\n",
			 m,I,K,J,Umin,N,UL,FLj,FLikj,(FLIK_FIX ? UL : UL*FLj));
		  printf("\t Xm[J]=%0.8e,IvarX=%0.8e,Yik=%0.8e,Yiu=%0.8e,Yvu=%0.8e\n",Xm[J],IvarX,y,y-Y[Umin],Y[Umin+1]-Y[Umin]);
		  printf("tid=%d:sleeping for 60 seconds\n",tid);
		  fflush(stdout);
		  
		  sleep(60);// allow other threads to terminate, so rverb does not produce too much output from other threads

		  rverb = 1;
		  MFLOAT FLj2 = FLJ(J,Xm[J],IvarX);
		  MFLOAT UL2 = FLIK<0>(I,K,J,Umin,Umin+1,y,Xm[J],y-Y[Umin],Y[Umin+1]-Y[Umin],Y,IvarX, FLj2, END_FIX ? lc : 1);
		  rverb = 0;

		  printf("\t FLj2=%0.8e,FLik2=%0.8e\n",FLj2,UL2);
		  fflush(stdout);
		  assert(fabs(FLikj - (FLIK_FIX ? UL : UL * FLj)) <= FLikj * (USE_MFLOAT ? 1e-4 : 2e-9) + (USE_MFLOAT ? 1e-32 : 1e-305));
		}
	      }
	    }
	    if((LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M /* && p->UL+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)
	       || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && I == 2 && K == 0 && J == 1)){
	      #pragma omp critical
	      {
		printf("m=%d:I=%d,K=%d,J=%d,Lij=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e: FLj=%0.10e,PoutlierEndM1=%0.6e,PMik=%0.6e,FLIK=%0.10e),X[J]=%0.3f,Y[I]=%0.3f,Yc(I,K)=%0.4f,Y[U]=%0.3f,Y[U+1]=%0.3f,Rtheta=%0.10e,Frate=%0.10e,IvarX=%0.6e\n",
		       m,I,K,J,A.Lij(I,K,J),Umin,N,A.LR(I,K,J) + UL*NotOutlier,UL*NotOutlier, FLj, PoutlierEndM1,PMik, UL, Xm[J], Y[I], y, Y[Umin], Y[Umin+1],Rtheta,Frate,IvarX);
		fflush(stdout);
	      }
	    }
	    if(DEBUG>=2 && !(isfinite(UL) && UL >= 0.0 && (!END_FIX || USE_MFLOAT || UL > 0.0))){
	      printf("m=%d:I=%d,K=%d,J=%d,U=%d:UL=bUL=%0.8e,Lij=%d,Umin=%d,Umax=%d,X[J]=%0.3f,y=%0.3f,Y[U]=%0.3f,Y[U+1]=%0.3f,Rtheta=%0.8e,IvarX=%0.8e\n",
		     m,I,K,J,Umin,UL, A.Lij(I,K,J), Umin,Umax, Xm[J], y, Y[Umin],Y[Umin+1], Rtheta, 2.0*(F2var + S2var * Xm[J] + R2var * Xm[J] * Xm[J]));
	      fflush(stdout);
	      assert(isfinite(UL) && UL >= 0.0 && (!END_FIX || USE_MFLOAT || UL > 0.0));
	    }

	    for(int U = Umax; ++U < I-K;){
	      if(VERB>=3 && rverb && m==MDEBUG_M && I==33 && K==0 && J==4)
		rverb = 2;
	      MFLOAT LR = ERFC_FIX ? FLIK<1>(I,K,J,U,U+1,y,Xm[J],y - Y[U], Y[U+1] - Y[U], Y,IvarX, FLj, END_FIX ? lc : 0)
	                           : FLIK<0>(I,K,J,U,U+1,y,Xm[J],y - Y[U], Y[U+1] - Y[U], Y,IvarX, FLj, END_FIX ? lc : 0);
	      if(VERB>=3 && rverb>=2 && m==MDEBUG_M && I==33 && K==0 && J==4)
		rverb = 1;
	      if(URANGE_QFIX ? (LR <= UL * EPS) : (LR < UL * EPS)){
		if((LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M /* && p->UL+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)
		   || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && I == 2 && K == 0 && J == 1)){
		  #pragma omp critical
		  {
		    printf("m=%d:I=%d,K=%d,J=%d,U=%d,Umin=%d,N=%d:UL=%0.10e:delta=%0.10e(ignored):FLIK=%0.10e,X[J]=%0.3f,Y[I]=%0.3f,Yc(I,K)=%0.4f,Y[U]=%0.3f,Y[U+1]=%0.3f,Rtheta=%0.10e,Frate=%0.10e,IvarX=%0.6e\n",
			   m,I,K,J,U,Umin,N,A.LR(I,K,J) + UL*NotOutlier,LR*NotOutlier,LR, Xm[J], Y[I], y, Y[U], Y[U+1], Rtheta, Frate, IvarX);
		    fflush(stdout);
		  }
		}
		//		U = I;
		break;
	      }
	      Umax = U;
	      if(LR > bUL){
		if(URANGE_FIX) {
		  Umin = U;
		  UL = 0.0;
		}
		bUL = LR;
		if(!MISSCNT_FIX)
		  A.Lij(I,K,J) = min(I-K,U+1);
	      }
	      if((LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M /* && p->UL+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)
		 || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && I == 2 && K == 0 && J == 1)){
		#pragma omp critical
		{
		  printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e):FLIK=%0.10e,X[J]=%0.3f,Y[I]=%0.3f,Yc(I,K)=%0.4f,Y[U]=%0.3f,Y[U+1]=%0.3f,Rtheta=%0.10e,Frate=%0.10e,IvarX=%0.6e\n",
			 m,I,K,J,U,N,A.LR(I,K,J) + (UL+LR)*NotOutlier,LR*NotOutlier,LR, Xm[J], Y[I], y, Y[U], Y[U+1], Rtheta, Frate, IvarX);
		  fflush(stdout);
		}
	      }
	      UL += LR;
	    }

	    for(int U = Umin; --U >= 0;){
	      if(VERB>=3 && rverb && m==MDEBUG_M && I==33 && K==0 && J==4)
		rverb = 2;
	      MFLOAT LR = (ERFC_FIX && U > Ub) ? FLIK<1>(I,K,J,U,U+1, y, Xm[J] , y - Y[U], Y[U+1] - Y[U], Y, IvarX, FLj, END_FIX ? lc : 1)
		: FLIK<0>(I,K,J,U,U+1, y, Xm[J] , y - Y[U], Y[U+1] - Y[U], Y, IvarX, FLj, END_FIX ? lc : 1);
	      if(VERB>=3 && rverb>=2 && m==MDEBUG_M && I==33 && K==0 && J==4)
		rverb = 1;
	      if(URANGE_QFIX ? (LR <= UL * EPS) : (LR < UL * EPS)){
		if((LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M /* && A.UL(I,K,J)+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)
		   || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && I == 2 && K == 0 && J == 1)){
		  #pragma omp critical
		  {
		    printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e, delta = %10e(ignored):FLIK=%0.10e,X[J]=%0.3f,Y[I]=%0.3f,Yc(I,K)=%0.4f,Y[U]=%0.3f,Y[U+1]=%0.3f,Rtheta=%0.10e,Frate=%0.10e,IvarX=%0.6e\n",
			   m,I,K,J,U,N,A.LR(I,K,J) + UL*NotOutlier,LR*NotOutlier, LR, Xm[J], Y[I], y, Y[U], Y[U+1], Rtheta, Frate, IvarX);
		    fflush(stdout);
		  }
		}
		break;
	      }
	      Umin = U;
	      if(LR > bUL){
		bUL = LR;
		if(!MISSCNT_FIX)
		  A.Lij(I,K,J) = U+1;/* NOTE : refalign() later uses Lij to compute misaligned sites. Here Lij is only used to compute coverage endpoint and only if NeedBP == 1 (normally not true) */
	      }
	      if((LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M /* && A.UL(I,K,J)+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)
		 || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && I == 2 && K == 0 && J == 1)){
                #pragma omp critical
		{
		  printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e:FLIK=%0.10e),X[J]=%0.3f,Y[I]=%0.3f,Yc(I,K)=%0.4f,Y[U]=%0.3f,Y[U+1]=%0.3f,Rtheta=%0.10e,Frate=%0.10e,IvarX=%0.6e\n",
			 m,I,K,J,U,N,A.LR(I,K,J) + (UL+LR)*NotOutlier,LR*NotOutlier,LR, Xm[J],Y[I],y, Y[U], Y[U+1], Rtheta,Frate,IvarX);
		  fflush(stdout);
		}
	      }
	      UL += LR;
	    }

	    if(ISLINEAR && lc){
	      double LR = FLEIK(I,K,J,0,1,y,Xm[J], FLj);
	      if(LR > bUL){
		bUL = LR;
		if(!MISSCNT_FIX)
		  A.Lij(I,K,J) = 0;
		if(bUL * NotOutlier > A.bLR(I,K,J))
		  A.G(I,K,J) = A.H(I,K,J) = A.T(I,K,J) = 0;
	      }
	      UL += LR;
	    }
	    if(DEBUG>=2)assert(isfinite(UL) && (!END_FIX || UL > 0.0));

	    UL *= NotOutlier;
	    bUL *= NotOutlier;
	    if(bUL <= A.bLR(I,K,J)){/* p->bLR corresponds to local/outlier score */
	      A.Lij(I,K,J) = OUTLIER_FIX ? I-K : -1 - A.Lij(I,K,J);
	      if(DEBUG) assert(A.G(I,K,J) == -2);
	    } else {
	      if(DEBUG>=2 && !lc && PoutlierEnd >= 1.0){
		printf("m=%d:I=%d,K=%d,J=%d:UL=%0.8e(bUL=%0.8e),PMik=%0.8e,p->LR=%0.8e(bLR=%0.8e),Lij=%d,Umin=%d,Umax=%d,X[J]=%0.3f,y=%0.3f,Y[%d]=%0.3f,FLj=%0.8e,PoutlierEndM1=%0.8e,NotOutlier=%0.8e\n",
		       m,I,K,J,UL,bUL,PRtab[K][I].PM,A.LR(I,K,J),A.bLR(I,K,J),A.Lij(I,K,J),Umin,Umax,Xm[J],y,Umax,Y[Umax], FLj, PoutlierEndM1, NotOutlier);
		fflush(stdout);
		assert(bUL <= A.bLR(I,K,J));
	      }
	      if(A.G(I,K,J) == -2)
		A.G(I,K,J) = -1;/* no futher aligned sites, but no local alignment */
	      A.bLR(I,K,J) = bUL;
	    }
	    if((LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0 && N==MDEBUG_N3 && m==MDEBUG_M /* && p->LR+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J) 
	       || (VERB>=3 && rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M) && I == 2 && K == 0 && J == 1)){
	      #pragma omp critical
	      {
		printf("m=%d:I=%d,K=%d,J=%d:ULmin=%d,ULmax=%d,N=%d:UL= %0.10e -> %0.10e(delta=%0.10e,FLj=%0.10e,PoutlierEndM1=%0.10e,PMik=%0.10e)\n",
		       m,I,K,J,Umin,Umax,N,A.LR(I,K,J),A.LR(I,K,J) + UL,UL,FLj,PoutlierEndM1,PMik);
		fflush(stdout);
	      }
	    }

	    A.LR(I,K,J) += UL;

	    if(DEBUG>=2 && !(isfinite(A.LR(I,K,J)) && (!END_FIX || PMik <= 0.0 || PoutlierEnd <= 0.0 || A.LR(I,K,J) > 0.0))){
	      printf("m=%d:I=%d,K=%d,J=%d:UL=%0.8e(bUL=%0.8e),PMik=%0.8e,FLj=%0.8e,PoutlierEndM1=%0.8e,p->LR=%0.8e,Lij=%d,Umin=%d,Umax=%d,X[J]=%0.3f,y=%0.3f,Y[%d]=%0.3f\n",
		     m,I,K,J,UL,bUL,PMik, FLj,PoutlierEndM1, A.LR(I,K,J), A.Lij(I,K,J),Umin,Umax,Xm[J],y,Umax,Y[Umax]);
	      printf("NotOutlier=%0.8e\n",NotOutlier);
	      fflush(stdout);
	      assert(isfinite(A.LR(I,K,J)) && (!END_FIX || PMik <= 0.0 || PoutlierEnd <= 0.0 || A.LR(I,K,J) > 0.0));
	    }
	  } // K = 
	} // ImK = 
      } // J = 

      double wt5;
      if(TIME_VERB)
        wt5 = wtime();	       

      /* initialize A[I=Imin[J]..Imax[J]][K=0..Kmax[I]][J= M+1-DELTA_XU to M].Rij as rightmost site/end in Y that overlaps Xm when
	 Y[I-K..I] is aligned with Xm[J] */
      /* only cases where J >= M+1 - DELTA_XU (OR EXTEND && I >= N+1-DELTA_YU && J >= M+1 - DELTA_X - EXTEND_MAX) */
      /* NOTE : A.Rij is needed later to compute UR subject to EPS, so score matches what mprobeval computes */
      int jmin = EXTEND ? max(1,M+1-DELTA_XU - EXTEND_MAX) : max(1,M+1 - DELTA_XU);
      for(int J= jmin; J <= M; J++){
	int Rij = N+1;
	if(DEBUG>=2 && Imin[J] <= Imax[J]) assert(IMIN <= Imin[J] && Imin[J] <= Imax[J] && Imax[J] <= IMAX);
	int imin = (J >= M+1-DELTA_XU) ? Imin[J] : max(Imin[J], N+1-DELTA_YU);
	for(int I= Imax[J]; I >= imin; I--){
	  if(DEBUG>=2) assert(1 <= Jmin[I] && Jmin[I] <= J && J <= Jmax[I] && Jmax[I] <= M);
	  int Rijk = Rij;
	  for(int K = 0; K <= Kmax[I]; K++){
	    double y = Yc(Y,I,K);
	    /* update Rijk by decrementing it until Y[Rijk] overlaps X */
	    while(Y[Rijk]-y > Xm[M+1]-Xm[J])
	      Rijk--;
	    if(DEBUG>=2) assert(Rijk >= I-K && !(Y[Rijk]-y > Xm[M+1]-Xm[J]) && (Rijk > N || Y[Rijk+1]-y > Xm[M+1]-Xm[J]));
	    if(VERB>=2 && rverb>=2 && Rijk > N){
	      printf("m=%d/%d:mapid=%d,id=%lld:N=%d,M=%d,I=%d,K=%d,J=%d,Rijk -> %d,Yik=%0.3f,Y[Rijk]=%0.3f,X[M+1]=%0.3f,X[J]=%0.3f\n",
		     m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,N,M,I,K,J,Rijk,y,Y[Rijk],Xm[M+1],Xm[J]);
	      fflush(stdout);
	    }
	    A.Rij(I,K,J) = max(Rijk,I);
	    if(K==0)
	      Rij = Rijk;
	  }
	}
      }

      double wt5b;
      if(TIME_VERB)
        wt5b = wtime();	       

      if(DEBUG) assert(!DELTA_MIN);
      if(DEBUG) assert(PM_SPLIT);

      if(VERB && (VERB>=2 || (rverb && (gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m == MDEBUG_M))) && !(MDEBUG_M >= 0 && rverb && m==MDEBUG_M)){
	printf("m=%d/%d:id=%lld:global=%d,IMIN=%d,IMAX=%d:\n",m,MD,gmap[pcontig->contig[m].mapid]->id,global,IMIN,IMAX);
	for(int I=IMIN; I <= IMAX; I++)
	  printf("  I=%d:Kmax[I]=%d,Jmin[I]=%d,Jmax[I]=%d,Y[I]=%0.3f\n",I,Kmax[I],Jmin[I],Jmax[I],Y[I]);
	fflush(stdout);
      }

      MFLOAT **delXm = pcontig->delX[0][m];

      /* Dynamic programming 3-D recurrance */
      if(OUT_TYPE){
	if(OUT_OLD)
	  qprobeval3D<1,1>(A,m,MD,IMIN,IMAX,Kmax,Jmin,Jmax,M,N,LogScale,ScaleI,Xm,delXm,Y,pcontig);	  
	else
	  qprobeval3D<1,0>(A,m,MD,IMIN,IMAX,Kmax,Jmin,Jmax,M,N,LogScale,ScaleI,Xm,delXm,Y,pcontig);	  
      } else {
	if(OUT_OLD)
	  qprobeval3D<0,1>(A,m,MD,IMIN,IMAX,Kmax,Jmin,Jmax,M,N,LogScale,ScaleI,Xm,delXm,Y,pcontig);
	else
	  qprobeval3D<0,0>(A,m,MD,IMIN,IMAX,Kmax,Jmin,Jmax,M,N,LogScale,ScaleI,Xm,delXm,Y,pcontig);
      }

      double wt6;
      if(TIME_VERB)
        wt6 = wtime();	       

      /* Compute LRcum (total LR) and locate best alignment result : only check results where J >= M+1-DELTA_XU (OR EXTEND && I >= N+1-DELTA_YU && J >= M+1 - DELTA_XU - EXTEND_MAX)
	 unless PoutlierEnd > 0.0, in which case check end outliers for all cases (within the bounds of IMIN,IMAX,Jmin[],Jmax[]) */
      int bestI = -1, bestK = 0, bestJ = -1, bestR = -2;
      double bestLR = 0.0, LRcum = 0.0, LS = LogScale[IMIN], bestUR = 0.0;

      for(int I = IMIN; I <= IMAX; I++){
	if(DEBUG>=2) assert(isfinite(LogScale[I]));
	double Scale = exp(LogScale[I] - LS);
	if(DEBUG && !(isfinite(Scale) && Scale >= 0.0 && Scale < MaxScale * ((MFLOAT)1.01))){
	  #pragma omp critical
	  {
	    printf("m=%d,I=%d:Scale=%0.6e,LogScale[I]=%0.8e,LS=%0.8e,IMIN=%d,IMAX=%d,MaxScale=%0.3e\n",m,I,Scale,LogScale[I],LS,IMIN,IMAX,MaxScale * ((MFLOAT)1.01));
	    for(int i = IMIN; i <= IMAX; i++)
	      printf("LogScale[%d]=%0.8e\n", i, LogScale[i]);
	    fflush(stdout);
	    assert(isfinite(Scale) && Scale >= 0.0 && Scale < MaxScale * ((MFLOAT)1.01));
	  }
	}
	for(int K = 0; K <= Kmax[I]; K++){
	  int jmin = Jmin[I];
	  int jmax = Jmax[I];
	  MFLOAT1 *pLR = &A.LR(I,K,0);
	  int *pRij = &A.Rij(I,K,0);

	  if(PoutlierEnd > 0.0){/* consider right end outlier (local alignment) for all cases (not just near right end of X) */
	    double endLR = Scale * PoutlierEnd;
	    if(DEBUG>=2) assert(isfinite(endLR));

	    for(int J = jmin; J <= jmax; J++){
	      // RecA *p = &A[I][K][J];
	      if(DEBUG>=2) assert(isfinite(pLR[J]));
	      double LR = pLR[J] * endLR;
	      if(DEBUG && ((DEBUG>=2 && !isfinite(LR)) || (MTRACE && MDEBUG_M>=0 && N==MDEBUG_N && m==MDEBUG_M && LR >= MDEBUG_MIN) ||
			   (LTRACE && MDEBUG_I >= 0 && MDEBUG_M>=0 && N==MDEBUG_N3 && m==MDEBUG_M && LR >= MDEBUG_MIN && MDEBUG_I == I && MDEBUG_K==K && MDEBUG_J == J)
			   || (VERB>=3 && rverb && gmap[pcontig->contig[m].mapid]->id == TRACE_ID && ((K==0 && J==54 && I == 9786)/* ||
														       (K==0 && J==93 && I == 5187) ||
														       (K==0 && J==92 && I == 5186) ||
														       (K==0 && J==91 && I == 5184) ||
														       (K==0 && J==90 && I == 5183) ||
														       (K==0 && J==89 && I == 5182) ||
														       (K==0 && J==88 && I == 5181)*/) /* && m==MDEBUG_M && LR > LRcum * 1e-6 */))){/* Note : LR can underflow and be 0 */
		printf("m=%d,I=%d,K=%d,J=%d:AL=%0.10e,UR=PoutlierEnd=%0.10e,LogScale[I]=%0.6e,LS=%0.6e:LRcum = %0.10e(delta=%0.10e)\n",
		  m,I,K,J,A.LR(I,K,J) * Scale,PoutlierEnd,LogScale[I],LS,LRcum+LR,LR);
		fflush(stdout);
		assert(isfinite(LR));
	      }
#if URDEBUG
	      A.UR(I,K,J) = PoutlierEnd;
#endif
	      LRcum += LR;
	      if(DEBUG) assert(isfinite(LRcum) && LRcum < DBL_MAX);

	      //  NOTE : to defer computing bestI,bestK,bestJ,bestR we need to save A.UR
	      if(LR > bestLR){
		bestLR = LR;
		bestI = I;
		bestK = K;
		bestJ = J;
		bestR = -2;
		if(URDEBUG)
		  bestUR = PoutlierEnd;
	      }
	    }
	  } else {
#if URDEBUG
	    for(int J = jmin; J <= jmax; J++)
	      A.UR(I,K,J) = 0.0;
#endif
	  }

	  double Yik = Yc(Y,I,K);
	  jmin = (EXTEND && I >= N+1 - DELTA_YU) ? max(Jmin[I],M+1 - DELTA_XU - EXTEND_MAX) : max(Jmin[I],M+1-DELTA_XU);
	  if(VERB>=3 && rverb && gmap[pcontig->contig[m].mapid]->id == TRACE_ID && K==0 && I == 9786){
	    printf("m=%d,I=%d,K=%d:jmin=%d,Jmin[I]=%d,Jmax[I]=%d:N=%d,M=%d,DELTA_YU=%d,DELTA_XU=%d\n",
		   m,I,K,jmin,Jmin[I],Jmax[I],N,M,DELTA_YU,DELTA_XU);
	    fflush(stdout);
	  }

	  for(int J = jmin; J <= Jmax[I]; J++){
	    MFLOAT delX = Xm[M+1]-Xm[J];
	    int Rij = pRij[J];
	    MFLOAT LRijk = pLR[J];

	    if(DEBUG>=2 && !(LRijk <= MaxLR)){
              #pragma omp critical
	      {
	        printf("m=%d:M=%d,N=%d:I=%d,K=%d,J=%d:LRijk=%0.10e, MaxLR=%0.10e\n",m,M,N,I,K,J,LRijk,MaxLR);
		fflush(stdout);
		assert(LRijk <= MaxLR);
	      }
	    }

	    if((EXTEND || Rij <= N) &&  LRijk > 0.0){
	      MFLOAT UR,bUR;
	      int linear = -1;/* end type : -1 = normal, 0 = chr end, -2 = local/outlier */
	      int Umin,Umax,Ub;
	      Ub = Umin = Umax = min(N,Rij)+1;
	      bUR = UR = FR<0>(I,K,J,Umin,Umin-1,Yik,delX,Y[Umin] - Yik, Y[Umin] - Y[Umin-1], Xm,Y,M,N,END_FIX ? rc : 1);

	      if((DEBUG>=2 && !(isfinite(UR) && UR >= 0.0 && UR <= MaxUR ))
		 || (VERB>=3 && rverb && gmap[pcontig->contig[m].mapid]->id == TRACE_ID && K==0 && J==54 && I == 9786)){
		#pragma omp critical
		{
		  printf("m=%d,I=%d,K=%d,J=%d,U=Ub=%d,M=%d,N=%d,Rij=%d:UR=FR=%0.10e:x=%0.6f,y=%0.6f,b=%0.6f\n",
			 m,I,K,J,Ub,M,N,Rij,UR,Xm[M+1]-Xm[J],Y[Umin]-Yik,Y[Umin] - ((END_FIX>=2 && Umin-1==I) ? Yik : Y[Umin-1]));
		  //		  rverb = 1;
		  MFLOAT FR2 = FR<0>(I,K,J,Umin,Umin-1,Yik,Xm[M+1]-Xm[J],Y[Umin]-Yik,Y[Umin]-Y[Umin-1],Xm,Y,M,N,END_FIX ? rc : 1);
		  //		  rverb = 0;
		  printf("\t FR2=%0.8e\n",FR2);
		  fflush(stdout);
		  assert(isfinite(UR) && UR >= 0.0 && UR <= MaxUR);
		}
	      }

	      for(int U = Umin; --U > I; ){
		MFLOAT FRu = ERFC_FIX ? FR<1>(I,K,J,U,U-1,Yik,delX,Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc :1)
		  : FR<0>(I,K,J,U,U-1,Yik,delX,Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc :1);
		if(DEBUG>=2) assert(isfinite(FRu) && FRu >= 0.0 && FRu <= MaxUR);
		if(URANGE_QFIX ? (FRu <= UR * EPS) : (FRu < UR * EPS)){
		  if(VERB>=3 && rverb && (m==MDEBUG_M || gmap[pcontig->contig[m].mapid]->id == TRACE_ID) && K==0 && J==54 && I == 9786){
		    printf("m=%d,I=%d,K=%d,J=%d,U=%d:UR=%0.10e(FRu=%0.10e),bUR=%0.10e,Rij=%d,Umax=%d:x=%0.3f,y=%0.3f,b=%0.3f: ignoring FRu\n",
			   m,I,K,J,Umin,UR,FRu,bUR,pRij[J],Umax,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U] - ((END_FIX>=2 && U-1==I) ? Yik : Y[U-1])); 
		    fflush(stdout);
		  }
		  break;
		}
		Umin = U;
		if(FRu > bUR){
		  if(URANGE_FIX){
		    Umax = U;
		    UR = 0.0;
		  }
		  bUR = FRu;
		  if(!MISSCNT_FIX)
		    pRij[J] = U-1;/* NOTE : In refalign.cpp Rij is used to compute misaligned sites. Here Rij will only be used to compute coverage endpoint */
		}
		UR += FRu;
		if(DEBUG>=2) assert(isfinite(UR) && UR >= 0.0 && UR <= MaxUR);
		if(VERB>=3 && rverb && (m==MDEBUG_M || gmap[pcontig->contig[m].mapid]->id == TRACE_ID) && K==0 && J==54 && I == 9786){
		  printf("m=%d,I=%d,K=%d,J=%d,U=%d:UR=%0.10e(FRu=%0.10e),bUR=%0.10e,Rij=%d,Umax=%d:x=%0.3f,y=%0.3f,b=%0.3f\n",
			 m,I,K,J,Umin,UR,FRu,bUR,pRij[J],Umax,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U] - ((END_FIX>=2 && U-1==I) ? Yik : Y[U-1])); 
		  fflush(stdout);
		}
	      }

	      for(int U = Umax; ++U <= N+1; ){
		MFLOAT FRu = (ERFC_FIX && U < Ub) ? FR<1>(I,K,J,U,U-1,Yik,delX,Y[U] - Yik, Y[U] - Y[U-1], Xm,Y,M,N,END_FIX ? rc : 1)
		  : FR<0>(I,K,J,U,U-1,Yik,delX,Y[U] - Yik, Y[U] - Y[U-1], Xm,Y,M,N,END_FIX ? rc : 1);
		if(DEBUG>=2) assert(isfinite(FRu) && FRu >= 0.0 && FRu <= MaxUR);
		if(URANGE_QFIX ? (FRu <= UR * EPS) : (FRu < UR * EPS)){
		  if(VERB>=3 && rverb && (m==MDEBUG_M || gmap[pcontig->contig[m].mapid]->id == TRACE_ID) && K==0 && J==54 && I == 9786){
		    printf("m=%d,I=%d,K=%d,J=%d,U=%d:UR=%0.10e(FRu=%0.10e),bUR=%0.10e,Rij=%d,Umax=%d:x=%0.3f,y=%0.3f,b=%0.3f: ignoring FRu\n",
			   m,I,K,J,Umin,UR,FRu,bUR,pRij[J],Umax,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U] - ((END_FIX>=2 && U-1==I) ? Yik : Y[U-1]));
		    fflush(stdout);
		  }

		  break;
		}
		Umax = U;
		UR += FRu;
		if(DEBUG>=2) assert(isfinite(UR) && UR >= 0.0 && UR <= MaxUR);
		if(FRu > bUR){
		  bUR = FRu;
		  if(!MISSCNT_FIX)
		    pRij[J] = U-1;
		}
		if(VERB>=3 && rverb && (m==MDEBUG_M || gmap[pcontig->contig[m].mapid]->id == TRACE_ID) && K==0 && J==54 && I == 9786){
		  printf("m=%d,I=%d,K=%d,J=%d,U=%d:UR=%0.10e(FRu=%0.10e),bUR=%0.10e,Rij=%d,Umax=%d:x=%0.3f,y=%0.3f,b=%0.3f\n",
			 m,I,K,J,Umin,UR,FRu,bUR,pRij[J],Umax,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U] - ((END_FIX>=2 && U-1==I) ? Yik : Y[U-1]));
		  fflush(stdout);
		}
	      }
	      if(ISLINEAR && rc){
		MFLOAT FRu = FRE(I,K,J,N+1,N,(MFLOAT)(Xm[M+1]-Xm[J]),(MFLOAT)(Y[N+1]-Yik),M);
		if(DEBUG>=2) assert(isfinite(FRu) && FRu >= 0.0 && FRu <= MaxUR);
		if(FRu > bUR){
		  bUR = FRu;
		  if(!MISSCNT_FIX)
		    pRij[J] = N+1;
		  linear = 0;
		}
		UR += FRu;
		if(DEBUG>=2) assert(isfinite(UR) && UR >= 0.0 && UR <= MaxUR);
	      }

	      /* LRcum += A.LR (I,K,J) * Scale * PoutlierEnd has already been processed previously */
	      UR = UR * PoutlierEndM1;
	      bUR = bUR * PoutlierEndM1;
	      if(bUR <= PoutlierEnd){
		if(DEBUG>=2) assert(A.Rij(I,K,J) > 0);
		pRij[J] = OUTLIER_FIX ? I : - pRij[J];
		linear = -2;/* local alignment OR endoutlier */
		//	      bUR = PoutlierEnd;
	      }
	      if(DEBUG>=2 && !(isfinite(UR) && (!END_FIX || PoutlierEndM1 <= 0.0 || ((0 && EXTEND_MAX) ? (UR > 0.0) : (UR >= 0.0))))){
		#pragma omp critical
		{
		  printf("m=%d:M=%d,N=%d:I=%d,K=%d,J=%d:UR=%0.8e,bUR=%0.8e,Umin=%d,Umax=%d,Rij=%d,A.Rij=%d,Yik=%0.3f(%0.3f),Y[Rij]=%0.3f,Y[Rij+1]=%0.3f,Y[N+1]=%0.3f,Xm[J]=%0.3f,Xm[M+1]=%0.3f,EXTEND_MAX=%d,PoutlierEnd=%0.8e,PoutlierEndM1=%0.8e\n",
			 m,M,N,I,K,J,UR,bUR,Umin,Umax,Rij,A.Rij(I,K,J),Yik,Yc(Y,I,K),Y[Rij],Y[Rij+1],Y[N+1],Xm[J],Xm[M+1],EXTEND_MAX,PoutlierEnd,PoutlierEndM1);
		  fflush(stdout);
		  assert(isfinite(UR) && (!END_FIX || PoutlierEndM1 <= 0.0 ||  ((!USE_MFLOAT && EXTEND_MAX) ? (UR > 0.0) : (UR >= 0.0))));
		}
	      }

#if URDEBUG
	      A.UR(I,K,J) += UR; 
#endif
   	      if(DEBUG>=2 && !(UR <= MaxUR)){
		#pragma omp critical
		{
		  printf("m=%d:M=%d,N=%d:I=%d,K=%d,J=%d:UR=%0.10e, MaxUR=%0.4ee\n",m,M,N,I,K,J,UR,MaxUR);
		  fflush(stdout);
		  assert(UR <= MaxUR);
		}
	      }

	      double LR = (double)UR * (double)LRijk * Scale;/* Scale = exp(LogScale[I]-LS) */
	      double bLR = (double)bUR * (double)LRijk * Scale;/* Scale = exp(LogScale[I]-LS) */
	      if((DEBUG && !(isfinite(LR) && LR < DBL_MAX)) || (MTRACE && MDEBUG_M>=0 && N==MDEBUG_N && m==MDEBUG_M && LR >= MDEBUG_MIN) 
		 || (VERB>=3 && rverb && gmap[pcontig->contig[m].mapid]->id == TRACE_ID && K==0 && J==54 && I == 9786 /* m==MDEBUG_M && LR > LRcum * 1e-6*/)){/* Note : LR can underflow and be 0 */
		#pragma omp critical
		{
		  printf("m=%d,I=%d,K=%d,J=%d:A.AL(I,K,J)=%0.10e,Scale=%0.10e,UR=%0.10e,bUR=%0.10e", m,I,K,J,A.LR(I,K,J),Scale,UR,bUR);
#if URDEBUG
		  printf("(A.UR(I,K,J)=%0.10e)",A.UR(I,K,J));
#endif
		  printf("LogScale[I]=%0.6e,LS=%0.6e: LRcum=%0.10e,LR=%0.10e,bLR=%0.10e:", LogScale[I],LS,LRcum,LR,bLR);
		  printf("(best:LR=%0.6e,I=%d,K=%d,J=%d,UR=%0.6f)\n", bestLR,bestI,bestK,bestJ,bestUR);
		  fflush(stdout);
		  assert(isfinite(LR) && LR < DBL_MAX);
		}
	      }
	      LRcum += LR;

	      // NOTE : to defer computing bestI,bestK,bestJ,bestR we need to save A.UR
		 //	      int test = (rverb && gmap[gcontig->contig[m].mapid]->id == TRACE_ID && K==0 && J==94) ? 1 : 0;
	      if(bLR > bestLR /* test ? (G==5187 && T==0 && H==93) : bLR > bestLR */){
#if URDEBUG
		if(DEBUG>=2 && linear == -2){
		  printf("m=%d,I=%d,K=%d,J=%d,R=%d:UR=%0.6e(bUR=%0.6e,tot=%0.6e),LogScale[I]=%0.6e,LS=%0.6e,Scale=%0.6e,AL=%0.6e,LR=%0.6e,bLR=%0.6e",
			 m,I,K,J,linear,UR,bUR,A.UR(I,K,J),LogScale[I],LS,Scale,A.LR(I,K,J),LR,bLR);
		  printf("(LRcum=%0.6e,best:LR=%0.6e,I=%d,K=%d,J=%d,R=%d,UR=%0.6e)\n", LRcum,bestLR,bestI,bestK,bestJ,bestR,bestUR);
		  fflush(stdout);
		  assert(linear != -2);
		}
#endif
		bestLR = LR;
		bestI = I;
		bestK = K;
		bestJ = J;
		bestR = linear;
		if(URDEBUG)
		  bestUR = bUR;
	      }
	      if(VERB >= 3  || (DEBUG && !isfinite(LR)) || (DEBUG && !isfinite(LRcum)) 
		 || (HDEBUG && rverb && m==311 && MD==2005 && 
		     ((I==78 && K==0 && (J==24 || J==23)) || (I==79 && K==0 && J==25) || (I==77 && K==0 && J==23)))){
		printf("m=%d:I=%d,K=%d,J=%d:UR=%0.10e(LS=%0.8f),AL=%0.10e(Scale=%0.10e,LogScale[I]=%0.8f),LR=%0.10e,LRcum=%0.10e,bestLR=%0.10e(I=%d,K=%d,J=%d,UR=%0.6e)\n", 
		       m,I,K,J, UR, LS, A.LR(I,K,J), Scale, LogScale[I], LR, LRcum, bestLR,bestI,bestK,bestJ,bestUR);
		fflush(stdout);
	      }
	      if(DEBUG)assert(isfinite(LR));
	      if(DEBUG)assert(isfinite(LRcum));
	    }
	  } // for J = jmin .. Jmax[I]
	}// for K = 0 .. Kmax[I]
	if(1/* WAS LRcum > 1.0e+100 */){  
	  double origLS = LS, origLRcum = LRcum, origbestLR = bestLR, origLogScaleI1 = LogScale[I+1];
	  if(LRcum > LRCUM_MAX){/* update scaling factor LS -> LS+log(LRcum) */
	    bestLR /= LRcum;
	    if(VERB>=3){
	      printf("m=%d,I=%d:LRcum=%0.6e,bestLR=%0.6e,LS=%0.6f->%0.6f(1)\n",m,I,LRcum,bestLR,LS,LS+log(LRcum));
	      fflush(stdout);
	    }
	    LS += log(LRcum);
	    if(DEBUG>=2) assert(isfinite(LS));
	    LRcum = 1.0;
	  } 
	  if(LS < LogScale[I+1] - logMaxScale + 0.001){
	    if(VERB>=3){
	      printf("m=%d,I=%d:LRcum=%0.6e,bestLR=%0.6e,LS=%0.6f->%0.6f(2)\n",m,I,LRcum,bestLR,LS,LogScale[I+1] - logMaxScale + 0.001);
	      fflush(stdout);
	    }
	    double Scale = exp(LS - LogScale[I+1] + logMaxScale - 0.001);
	    LS = LogScale[I+1] - logMaxScale + 0.001;
	    if(DEBUG>=2) assert(isfinite(LS));
	    LRcum *= Scale;
	    bestLR *= Scale;
	  }
	  if(LRCUM_FIX && LRcum > 0.0 && LRcum <= MinFloat && LS > LogScale[I+1] - logMaxScale + 0.001){
	    double pLS = LS;
	    LS = LogScale[I+1] - logMaxScale + 0.001;
	    double Scale = exp(pLS - LS);
	    if(LRcum * Scale > 1.0){
	      Scale = 1.0 / LRcum;
	      LS = pLS - log(Scale);
	    }
	    LRcum *= Scale;
	    bestLR *= Scale;
	  }
	  if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 */){
	    printf("m=%d,I=%d:j=%d..%d:LRcum=%0.8e->%0.8e(bestLR=%0.8e),Scale=%0.8e,LS=%0.8f->%0.8f,log(LRcum)+LS= %0.8f,LogScale[I+1]=%0.8f\n",
		   m,I,Jmin[I],Jmax[I],origLRcum,LRcum,bestLR,Scale,origLS,LS,LRcum > 0.0 ? log(LRcum)+LS : 0.0, LogScale[I+1]);
	    fflush(stdout);
	  }

	  if(DEBUG && !(isfinite(LS) && isfinite(LRcum) && LRcum <= LRCUM_MAX && (LRcum > 0.0 || origLRcum <= 0.0))){
	    #pragma omp critical
	    {
	      printf("\nm=%d,I=%d(%d..%d),N=%d:LRcum=%0.6e,bestLR=%0.6e,LS=%0.6f,LogScale[I+1]=%0.6f\n",m,I,IMIN,IMAX,N,LRcum,bestLR,LS,LogScale[I+1]);
	      printf("   original values: LRcum=%0.6e,bestLR=%0.6e,LS=%0.6f, LogScale[I+1]= %0.6f,logMaxScale= %0.6f\n",origLRcum,origbestLR,origLS,origLogScaleI1,logMaxScale);
	      fflush(stdout);
	      if(DEBUG)assert(isfinite(LS));
	      if(DEBUG)assert(isfinite(LRcum) && LRcum <= LRCUM_MAX);
	    }
	  }
	  if(DEBUG && !(LRcum > 0.0 || origLRcum <= 0.0)){
	    #pragma omp critical
	    {
	      printf("WARNING: Best Alignment score for molecule m=%d underflowed!\n",m);
	      fflush(stdout);
	      //	      assert((LRcum > 0.0 || origLRcum <= 0.0));
	    }
	  } else if(DEBUG && !(isfinite(bestLR) && (bestLR > 0.0 || origbestLR <= 0.0))){
	    #pragma omp critical
	    {
	      printf("\nm=%d,I=%d:LRcum=%0.6e,bestLR=%0.6e,LS=%0.6f,LogScale[I+1]=%0.6f\n",m,I,LRcum,bestLR,LS,LogScale[I+1]);
	      printf("original values: LRcum=%0.6e,bestLR=%0.6e,LS=%0.6f\n",origLRcum,origbestLR,origLS);
	      fflush(stdout);
	      assert(isfinite(bestLR) && (bestLR > 0.0 || origbestLR <= 0.0));
	    }
	  }
	} // if(1)

        if(VERB>=3 && rverb && m==MDEBUG_M){
          #pragma omp critical
	  {
	    printf("m=%d/%d:mapid=%d,id=%lld,flip=%d:I=%d:LRcum=%0.10e, LS=%0.10f(B)\n",m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,pcontig->flip[m],I,LRcum, LS);
	    fflush(stdout);
	  }
	}
      }// for I = IMIN .. IMAX

      double wt7;
      if(TIME_VERB)
        wt7 = wtime();	       

      if(NeedBP){ // NOTE : If bestI,bestJ,bestK was defered, we need to compute them first

	/* compute A.bLR,A.G,A.T.A.H along best path */
	int I = bestI;
	int K = bestK;
	int J = bestJ;
	
	while(I > 0){ /* duplicate 3-D code : should compute same A.LR */
	  MFLOAT PMik = PRtab[K][I].PM;
	  if(DEBUG>=2) assert(isfinite(PMik) && PMik >= 0.0);
	  double Yik = Yc(Y,I,K);

	  MFLOAT1 LR = A.UL(I,K,J);
	  double bLR = A.bLR(I,K,J);
	  int bestG = A.G(I,K,J);
	  int bestH = A.H(I,K,J);
	  int bestT = A.T(I,K,J);
	  int Gmin = max(IMIN, I-K - deltaExtYRef);
	  int Gmin2 = max(IMIN,I-K - DELTA_Y);
	  int Hmin = max(1, J - deltaExtXRef);
	  int Hmin2 = max(1, J - DELTA_X);
	  if(VERB>=3 || (LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0)) LRtrace(m,I,K,J,M,N,Gmin,Gmin2,IMIN,Hmin,Hmin2,Jmin,Yik,Xm,Y,LR,exp(-LogScale[I]));
	  for(int G= I-K; --G >= Gmin;){
	    double Scale = exp(LogScale[G] - LogScale[I]);// WAS * PMik;
	    if(DEBUG>=2 && SCALE_MONOTONIC) assert(isfinite(Scale) && Scale <= 1.0);
	    if(DEBUG>=2 && !SCALE_MONOTONIC) assert(isfinite(Scale) && Scale <= DBL_MAX);
	    MFLOAT1 LRg = 0.0;
	    int hmin = (G >= Gmin2) ? max(Hmin,Jmin[G]) : max(Hmin2, Jmin[G]);
	    int hmax = min(J-1,Jmax[G]);
	    if(hmin <= hmax){// NEW10 : avoid segfault from attempt to stack allocate an array of -ve size : newLRs[hmax-hmin+1]
	      if(VERB>=3 || (LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0)) LRtrace2(m,I,K,J,G,M,N,Gmin,Hmin,Kmax,hmin,hmax,Jmin,Jmax,Scale*PMik);
	      for(int T = 0; T <= Kmax[G]; T++){
		MFLOAT Ivar;
		MFLOAT deltaY = Yik - Yc(Y,G,T);
		MFLOAT LRik = FAIK(deltaY,PRtab[I-K-G-1][I-K].Pr,I-K-G,I,K,T,N,VARtab[K*N+I],VARtab[T*N+G],Ivar,Y);
		MFLOAT1 *LRgt = &A.LR(G,T,0);
		double newLRs[hmax-hmin+1];

	        #if (USE_EPOW) && (DEBUG < 2)
	        #pragma omp simd reduction(+:LRg)
	        #endif
	        for(int H= hmax; H >= hmin;H--){/* critical loop */
		  if(DEBUG>=2) assert(G-T >= 1 && G <= N && H >= 1 && H <= M);
		  MFLOAT deltaX = Xm[J] - Xm[H];
		  MFLOAT FAj = FAJ(deltaX, deltaY, J-H, I-K-G, FBiasWTm[J-H], LogFpPow[J-H-1], LogFnPow[I-K-G-1], PrBiasWT[I-K]-PrBiasWT[G], LRik, Ivar, 1, OUT_TYPE, OUT_OLD);
		  if(DEBUG>=2) FAcheck(m,I,K,J,G,T,H,N,Xm,Y,Yik,deltaX,deltaY,LRik,Ivar,FAj, OUT_TYPE, OUT_OLD);
		  MFLOAT1 newLR = LRgt[H] * FAj;
		  if(DEBUG>=2) assert(isfinite(newLR) && newLR < MaxFloat);
		  LRg += newLR;
		  if(DEBUG>=2) assert(isfinite(LRg) && LRg < MaxFloat);
		  newLRs[hmax-H] = newLR * Scale;/* Scale = exp(LogScale[G]-LogScale[I]) */ /* WAS * PMik */
		}
		for (int H = hmax; H >= hmin; --H) {
		  const double newLR = newLRs[hmax-H];
		  if(newLR > bLR){
		    bLR = newLR;
		    bestG = G;
		    bestT = T;
		    bestH = H;
		  }
		} // H = hmax .. hmin
	      } // T = 0 .. Kmax[g]
	    } // if(hmin <= hmax)
	    if(DEBUG>=2 && hmin > hmax) assert(LRg == 0.0 && isfinite(LR) && LR < MaxFloat);
            if(DEBUG>=2) assert(isfinite(LRg) && LRg < MaxFloat);
	    //	    MFLOAT1 origLR = LR;
	    LR += LRg * Scale;

	    if(DEBUG>=2 && !(isfinite(LR) && LR < MaxFloat)){// NOTE : if this overflows, apply same fix as in original 3-D loop
	      #pragma omp critical
	      {
	        printf("m=%d:I=%d,K=%d,J=%d,G=%d,T=0..%d,H=%d..%d:LR= %0.8e : LRg= %0.8e, Scale= %0.8e,LogScale[G,I]=%0.8f,%0.8f\n",
		  m,I,K,J,G,Kmax[G],hmin,hmax,/* origLR, */LR,LRg,Scale,LogScale[G],LogScale[I]);
		fflush(stdout);
	        assert(isfinite(LR) && LR < MaxFloat);
	      }
	    }
	  } // for G = I-K-1 .. Gmin

	  if(DEBUG>=2 && !(isfinite(LR) && LR < MaxFloat)){
	    printf("m=%d:I=%d,K=%d,J=%d:LR= %0.8e\n", m,I,K,J,LR);
	    fflush(stdout);
	    assert(isfinite(LR) && LR < MaxFloat);
	  }

	  if((DEBUG && !(isfinite(LR) && LR < MaxFloat)) || (LTRACE && MDEBUG_I >= 0 && MDEBUG_M >= 0) || HDEBUG) LRtrace3(A,m,MD,I,K,J,Gmin,N,Kmax,LR,LogScale);
	  A.bLR(I,K,J) = bLR * PMik;
	  A.G(I,K,J) = bestG;
	  A.T(I,K,J) = bestT;
	  A.H(I,K,J) = bestH;

	  I = bestG;
	  K = bestT;
	  J = bestH;
	}
      }

      if(DEBUG) assert(LS >= 0.0);

      double LPm = -99999.0;
      double WTm = 0.0;
      double LRbiasLS = LRbias * exp(-LS);
      if(LRcum > 0.0 || LRbiasLS > 0.0){
	double LRtot = LRcum + LRbiasLS;
	if(DEBUG) assert(LRtot > 0.0);
	LPm = LS+log(LRtot);
	WTm = LRcum/LRtot;
	if(DEBUG>=1+RELEASE/* HERE HERE >=2 */ && !(isfinite(LPm) && isfinite(WTm))){
	  printf("qprobeval:m=%d/%d: LPm= %0.8e, WTm= %0.8e, LRcum= %0.8e, LRbiasLS= %0.8e, LRtot= %0.8e, LS= %0.8e\n", m,MD,LPm,WTm,LRcum,LRbiasLS,LRtot,LS);
	  fflush(stdout);
	  assert(isfinite(LPm));
	  assert(isfinite(WTm));
	}
      }
      if((MTRACE || MDEBUG) && MDEBUG_M >= 0 && m==MDEBUG_M){
	printf("m=%d:LPm=%0.10f:LRcum=%0.10e,LRbias=%0.8e,LS=%0.8f(N=%d)\n",m,LPm,LRcum,LRbias,LS,N);
	fflush(stdout);
      }
    
      if(QDEBUG || rverb){
        #pragma omp critical
        {
	  LLRsum += LPm * TBwtm;
	}
      }

      if(newLP) newLP[m] = LPm;
      if(mapWT) mapWT[m] = WTm;
      if(pcontig->LP) pcontig->LP[m] = LPm;
      if(pcontig->mapWT) pcontig->mapWT[m] = WTm * TBwtm;
      if(LogPV) LogPV[m] = 0.0;/* default if no alignement was found */

      if((VERB && (VERB>=2 || rverb || (MTRACE && MDEBUG_M>=0 && N==MDEBUG_N && m==MDEBUG_M))) || LogPV || poutlier){
	/* display scores and (if LogPV) compute LogPV[m] and outlier locations */
        if(VERB && (VERB>=2 || rverb || (MTRACE && MDEBUG_M>=0 && N==MDEBUG_N && m==MDEBUG_M))){ 
          #pragma omp critical
          {
	    printf("m=%d/%d:mapid=%d,id=%lld,flip=%d:LLRsum=%0.10f(delta=%14.10f,TBwt=%0.6f,WT=%0.6f,LS=%0.6f,LRcum=%0.10e,LRbias=%0.6e),bestLR/LRcum=%0.10f,M=%d(%d),N=%d,best I=%d,K=%d,J=%d,Y[I]=%0.3f,X[J]=%0.3f,X[M+1]=%0.3f,TID=%lld\n",
		   m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,pcontig->flip[m],LLRsum,LPm,TBwtm,WTm,LS,LRcum,LRbias,bestLR/LRcum,M,MX[m],N,bestI,bestK,bestJ,
		bestI>=1 ? Y[bestI] : -1.0, bestI>=1 ? Xm[bestJ] : -1.0, Xm[M+1],TRACE_ID);
	    fflush(stdout);
	  }
	}
	if(bestI>=1){
	  /* allocate space to store alignment in reverse order */
	  if(DEBUG) assert(M <= maxMX);
	  int U = 0;

	  int I = bestI;
	  int K = bestK;
	  int J = bestJ;

	  while(I > 0){
	    if(DEBUG){
	      if(DEBUG && !(I-K >= 1 && I <= N)){
		printf("m=%d:U=%d:I=%d,K=%d,J=%d,N=%d,M=%d\n",
		       m,U,I,K,J,N,M);
		fflush(stdout);
		assert(I-K >= 1 && I <= N);
	      }
	      assert(J >= 1 && J <= M);
	    }
	    Ilist[U] = I;
	    Klist[U] = K;
	    Jlist[U] = J;
	    U++;

	    int G = A.G(I,K,J);
	    int T = A.T(I,K,J);
	    int H = A.H(I,K,J);

	    if(DEBUG && !(G < 0 || T <= KMAX)){
	      printf("m=%d:U=%d:I=%d,K=%d,J=%d,N=%d,M=%d,T=%d,KMAX=%d,A.kmax=%d,kmax=%d\n",
		     m,U,I,K,J,N,M,T,KMAX,A.kmax,kmax);
	      printf("bestI=%d,bestK=%d,bestJ=%d\n",bestI,bestK,bestJ);
	      fflush(stdout);
	      assert(G < 0 || T <= KMAX);
	    }
	    
	    I = G;
	    K = T;
	    J = H;
	  }
	  if(DEBUG) assert(U>0);
	  if(DEBUG) assert(U<=M);
	  if(DEBUG) assert(I <= (lc ? 0 : -1));
	  if(DEBUG) assert(bestR <= (rc ? 0 : -1));

	  int lverb = 0; // (U>=2 && Yc(Y,Ilist[U-1],Klist[U-1]) <= 1416.0 && Yc(Y,Ilist[0],Klist[0]) >= 1433.0) ? 1 : 0;

	  if(VERB && (VERB>=2 || (rverb && (rverb>=2 || lverb || gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m == MDEBUG_M)) || (MTRACE && MDEBUG_M>=0 && N==MDEBUG_N && m==MDEBUG_M))){ 
            #pragma omp critical
	    {
	      printf("\t G=%d,T=%d,H=%d:Y[G-T]=%0.3f,X[H]=%0.3f,scale=%0.8f,log(LRcum)=%0.10f,UR=%0.10e,R=%d,L=%d\n",
		     Ilist[U-1],Klist[U-1],Jlist[U-1],Y[Ilist[U-1]-Klist[U-1]],Xm[Jlist[U-1]],
		     (U <= 1) ? 1.0 : (Yc(Y,bestI,bestK)-Yc(Y,Ilist[U-1],Klist[U-1]))/(Xm[bestJ]-Xm[Jlist[U-1]]), log(LRcum) + LS,bestUR,bestR,I);
	      fflush(stdout);
	    }
	  }
	  if(1){/* display segment by segment scores and compute logPV[m] */
	    int misscnt = 0;
	    int falsecnt = 0;
	    int OutlierCnt = 0;/* number of internal outliers */
	    double OutlierMN = 1.0;/* product over all outlier of m x n (the number of intervals in the internal outlier aligned segment) */
	    int EndOutlierCnt = (bestR <= -2 ? 1 : 0) + (I <= -2 ? 1 : 0);
	    
	    int F = 0;/* segment (aligned interval) count */

	    int G= I, D= K,H= J;
	    /*	    if(poutlier)
	      for(int j = 0; j <= M; j++)
	      poutlier[m][j] = 0;*/

	    for(int T= 0; T < U; T++, G = I, H = J, D = K){
	      I = Ilist[U-1-T];
	      J = Jlist[U-1-T];
	      K = Klist[U-1-T];
	      if(T > 0){
		double y = Yc(Y,I,K)-Yc(Y,G,D);
		double x = Xm[J]-Xm[H];

		if(VERB>=3 && rverb && gmap[pcontig->contig[m].mapid]->id == TRACE_ID && I==121 && K==0 && J==27 && G==119 && D==0 && H==25){
		  double y1 = Yc(Y,I,K) - Yc(Y,I-1,K);
		  double x1 = Xm[J] - Xm[J-1];

		  int origrverb = rverb;
		  rverb = 2;
		  double FAijk1 = FA(x1,y1,1,1,J,I,K,0,N,VARtab[K*N+I],VARtab[0*N+I-1],FBiasWTm[1],LogFpPow[0],LogFnPow[0], PrBiasWT[I-K]-PrBiasWT[I-1], PRtab[0][I-K].Pr,Y,1, OUT_TYPE, OUT_OLD);

		  printf("\tT=%d(A):y1=%0.6f,x1=%0.6f:FAijk1= %0.10e\n",T,y1,x1,FAijk1);
		  fflush(stdout);

		  double y2 = Yc(Y,I-1,K) - Yc(Y,G,D);
		  double x2 = Xm[J-1] - Xm[H];

		  double FAijk2 = FA(x2,y2,1,1,J-1,I-1,0,D,N,VARtab[0*N+I-1],VARtab[D*N+G],FBiasWTm[1],LogFpPow[0],LogFnPow[0], PrBiasWT[I-1-0]-PrBiasWT[G], PRtab[0][I-1-0].Pr,Y,1, OUT_TYPE, OUT_OLD);
		  printf("\tT=%d(A):y2=%0.6f,x2=%0.6f:FAijk2= %0.10e\n",T,y2,x2,FAijk2);
		  fflush(stdout);

		  rverb = origrverb;
		}

		//		int origrverb = rverb;
		//		if(rverb && gmap[pcontig->contig[m].mapid]->id == TRACE_ID && I==6 && K==0 && J==4) rverb = 2;
		double FAijk = FA(x,y,J-H,I-K-G,J,I,K,D,N,VARtab[K*N+I],VARtab[D*N+G],FBiasWTm[J-H],LogFpPow[J-H-1],LogFnPow[I-K-G-1], PrBiasWT[I-K]-PrBiasWT[G], PRtab[I-K-G-1][I-K].Pr,Y,1, OUT_TYPE, OUT_OLD);
		//		rverb = origrverb;

		//		if(rverb && gmap[pcontig->contig[m].mapid]->id == TRACE_ID && I==6 && K==0 && J==4) rverb = 2;
		double FAijk0 = FA(x,y,J-H,I-K-G,J,I,K,D,N,VARtab[K*N+I],VARtab[D*N+G],FBiasWTm[J-H],LogFpPow[J-H-1],LogFnPow[I-K-G-1],PrBiasWT[I-K]-PrBiasWT[G],PRtab[I-K-G-1][I-K].Pr,Y,0, OUT_TYPE, OUT_OLD);
		//		rverb = origrverb;

		int outlier = 0;
		if(FAijk0 < FAijk * 0.5){/* outlier score */
		  outlier = 1;
		  if(poutlier)
		    for(int j = H; j < J; j++)
		      poutlier[m][j] |= 1;
		  if(OUTLIER_PV){
		    OutlierCnt++;
		    int m = J-H;
		    int n = I-K-G;
		    OutlierMN *= m*n;
		  }
		} else {
		  Yfrag[F] = y;
		  Xfrag[F] = x;
		  F++;
		  misscnt += I-K-G-1;
		  falsecnt += J-H-1;
		}
		double FAijkBW = FAijk;
		if(VERB>=2 && rverb && biasWT < 1.0){// WARNING : this code may introduce a race condition by changing global score values : should be safe since rverb disables multithreading
		  #pragma omp critical
		  {
		    double origbiasWT = biasWT;
		    biasWT = 1.0;
		    score_initOutlier();

		    int origrverb = rverb;
		    rverb = 2;
		    FAijkBW = FA(x,y,J-H,I-K-G,J,I,K,D,N,VARtab[K*N+I],VARtab[D*N+G],FBiasWTm[J-H],LogFpPow[J-H-1],LogFnPow[I-K-G-1],PrBiasWT[I-K]-PrBiasWT[G],PRtab[I-K-G-1][I-K].Pr,Y,1, OUT_TYPE, OUT_OLD);
		    rverb = origrverb;

		    biasWT = origbiasWT;
		    score_initOutlier();
		  }
		}
		
		if(T==U-1) {	  /* include UR of right end segment */
		  int Rij = (bestR <= -2) ? I : min(N,A.Rij(I,K,J));
		  int Rijx = (bestR <= -2) ? J : M;
		  if(VERB && (VERB>=2 || (rverb && (rverb>=2 || lverb || gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M)) || (MTRACE && MDEBUG_M>=0 && N==MDEBUG_N && m==MDEBUG_M))){ 
		    #pragma omp critical 
		    {
#if URDEBUG
		      printf("   T=%d/%d:I=%d,K=%d,J=%d:Y[I-K]=%0.3f,Y[I]=%0.3f,Yc=%0.3f,X[J]=%0.3f,y=%0.6f,x=%0.6f,PMik=%0.10f,FA=%0.10e(FA0=%0.10e,outlier=%d,FAbw=%0.10e)\n\t log(A.LR(I,K,J))=%0.10f,bLR/LR=%0.10f, A.LR(I,K,J)= %0.10e, LogScale[I]= %0.10f, UR=%0.10e,log(UR)=%0.10f,Rij=%d,Rijx=%d,bestR=%d\n",
			     T,U,I,K,J,Y[I-K],Y[I],Yc(Y,I,K),Xm[J],y,x,PRtab[K][I].PM,FAijk,FAijk0, outlier, FAijkBW, log(A.LR(I,K,J)) + LogScale[I], A.bLR(I,K,J)/A.LR(I,K,J), A.LR(I,K,J), LogScale[I], A.UR(I,K,J),log(A.UR(I,K,J)), Rij, Rijx, bestR);
#else
		      printf("   T=%d/%d:I=%d,K=%d,J=%d:Y[I-K]=%0.3f,Y[I]=%0.3f,Yc=%0.3f,X[J]=%0.3f,y=%0.6f,x=%0.6f,PMik=%0.10f,FA=%0.10e(FA0=%0.10e,outlier=%d,FAbw=%0.10e)\n\t log(A[I][K][J].LR)=%0.10f,bLR/LR=%0.10f\n",
			     T,U,I,K,J,Y[I-K],Y[I],Yc(Y,I,K),Xm[J],y,x,PRtab[K][I].PM,FAijk,FAijk0, outlier, FAijkBW, log(A.LR(I,K,J)) + LogScale[I], A.bLR(I,K,J)/A.LR(I,K,J));
#endif
		      if(VERB>=3 && I >= 3)
			printf("       Y[I-3..I+3]=%0.6f,%0.6f,%0.6f,%0.6f,%0.6f,%0.6f,%0.6f,minKB=%0.6f\n",Y[I-3],Y[I-2],Y[I-1],Y[I],Y[I+1],Y[I+2],Y[I+3],minKB);
		      fflush(stdout);
		    }
		  }
		  if(I < Rij)
		    misscnt += Rij-I;
		  if(J < Rijx)
		    falsecnt += Rijx-J;
		} else {
		  if(VERB && (VERB>=2 || (rverb && (rverb>=2 || lverb || gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M)) || (MTRACE && MDEBUG_M>=0 && N==MDEBUG_N && m==MDEBUG_M))){ 
		    #pragma omp critical 
		    {
		      printf("   T=%d/%d:I=%d,K=%d,J=%d:Y[I-K]=%0.3f,Y[I]=%0.3f,Yc=%0.3f,X[J]=%0.3f,y=%0.6f,x=%0.6f,PMik=%0.10f,FA=%0.10e(FA0=%0.10e,outlier=%d,FAbw=%0.10e)\n\t log(A.LR(I,K,J))=%0.10f,bLR/LR=%0.10f, A.LR(I,K,J)= %0.10e, LogScale[I]= %0.10f\n",
			     T,U,I,K,J,Y[I-K],Y[I],Yc(Y,I,K),Xm[J],Yc(Y,I,K)-Yc(Y,G,D),Xm[J]-Xm[H],PRtab[K][I].PM,FAijk,FAijk0,outlier,FAijkBW,log(A.LR(I,K,J)) + LogScale[I], A.bLR(I,K,J)/A.LR(I,K,J), A.LR(I,K,J), LogScale[I]);
		      fflush(stdout);
		    }
		  }
		}
	      } else if(T==0) { /* left end segment */
		int Lij = (G <= -2) ? I-K : max(1,A.Lij(I,K,J));
		int Lijx = (G <= -2) ? J : 1;/* since X should never extends beyond Y */
		if(VERB && (VERB>=2 || (rverb && (rverb>=2 || lverb || gmap[pcontig->contig[m].mapid]->id == TRACE_ID || m==MDEBUG_M)) || (MTRACE && MDEBUG_M>=0 && N==MDEBUG_N && m==MDEBUG_M))){ 
		  #pragma omp critical 
		  {
		    printf("   T=%d/%d:I=%d,K=%d,J=%d,N=%d,M=%d:", T,U,I,K,J,N,M);
		    printf("Y[I-K]=%0.3f,",Y[I-K]);
		    printf("Y[I]=%0.3f,",Y[I]);
		    printf("Yc=%0.6f,", Yc(Y,I,K));
		    printf("X[J]=%0.6f,",Xm[J]);
		    printf("PMik=%0.6f,",PRtab[K][I].PM);
		    printf("A.LR(I,K,J)=%0.6e(G=%d,Lij=%d,Lijx=%d),",A.LR(I,K,J),G,Lij,Lijx);
		    printf("LogScale[I]=%0.6f,",LogScale[I]);
		    printf("A.bLR(I,K,J)=%0.6e\n",A.bLR(I,K,J));
		    printf("\t log(A.LR(I,K,J))=%0.10f,bLR/LR=%0.10f, A.LR(I,K,J)= %0.10e, LogScale[I]= %0.10f\n", log(A.LR(I,K,J))+LogScale[I], A.bLR(I,K,J)/A.LR(I,K,J), A.LR(I,K,J),LogScale[I]);
		    fflush(stdout);
		  }
		}
		if(I-K > Lij)
		  misscnt += I-K-Lij;
		if(J > Lijx)
		  falsecnt += J-Lijx;
	      }
	    }

	    double logOutlierMN = 0.0;
	    if(OUTLIER_PV && OutlierMN > 1e+300){
	      for(int T = 0; T < U; T++, G = I, H = J, D = K){
		I = Ilist[U-1-T];
		J = Jlist[U-1-T];
		K = Klist[U-1-T];
		if(T > 0){
		  double y = Yc(Y,I,K)-Yc(Y,G,D);
		  double x = Xm[J]-Xm[H];
		  double FAijk = FA(x,y,J-H,I-K-G,J,I,K,D,N,VARtab[K*N+I],VARtab[D*N+G],FBiasWTm[J-H],LogFpPow[J-H-1],LogFnPow[I-K-G-1],PrBiasWT[I-K]-PrBiasWT[G],PRtab[I-K-G-1][I-K].Pr,Y,1,OUT_TYPE, OUT_OLD);
		  double FAijk0 = FA(x,y,J-H,I-K-G,J,I,K,D,N,VARtab[K*N+I],VARtab[D*N+G],FBiasWTm[J-H],LogFpPow[J-H-1],LogFnPow[I-K-G-1],PrBiasWT[I-K]-PrBiasWT[G],PRtab[I-K-G-1][I-K].Pr,Y,0,OUT_TYPE, OUT_OLD);
		  if(FAijk0 < FAijk * 0.5){/* outlier score */
		    int m = J-H;
		    int n = I-K-G;
		    logOutlierMN += log((double)(m*n));
		  }
		}
	      }
	    } else
	      logOutlierMN = log(OutlierMN);
	    int EndOutlierNL = 1;
	    if(OUTLIER_PV && EndOutlierPV && EndOutlierCnt >= 1){
	      if(EndOutlierCnt >= 2){
		int L = min(Ilist[0] - Ilist[U-1], Jlist[0] - Jlist[U-1]);
		int T = min(N, M);
		EndOutlierNL = T - L + 1;
		if(DEBUG) assert(EndOutlierNL > 0);
	      } else
		EndOutlierNL = 2;
	    }

	    double logPV = pvalue(Xfrag,Yfrag,F, misscnt,falsecnt,OutlierCnt,logOutlierMN,EndOutlierNL,Ylambda,res[0]*PixelLen, (rverb>=2 ? 1 : 0), 0, pcontig->contig[m].mapid, 0);
	    if(LogPV) LogPV[m] = logPV;
	    if(VERB && (VERB>=2 || rverb || (MTRACE && MDEBUG_M>=0 && N==MDEBUG_N && m==MDEBUG_M))){ 
              #pragma omp critical 
	      {
		//		printf("   N=%d,Y[1]=%0.3f,Y[N]=%0.3f,Y[N+1]=%0.3f,M=%d,X[M+1]=%0.3f\n",N,Y[1],Y[N],Y[N+1],M,Xm[M+1]);
		if(gmap[pcontig->contig[m].mapid]->align->mapid1 >= 0)
		  printf("   logLR=%0.6f(wt=%0.6f), LogFP = %0.4f (original logLR score=%0.6f, LogFP= %0.4f)\n",
			 log(max(1e-300,LRcum+LRbiasLS))+LS-2.0*log(Rtheta),WTm,logPV, gmap[pcontig->contig[m].mapid]->align->score, gmap[pcontig->contig[m].mapid]->align->logPV);
		else
		  printf("   logLR=%0.6f(wt=%0.6f), LogFP = %0.4f (LRcum=%0.10e,LRbiasLS=%0.10f,LS=%0.10f,log(Rtheta)=%0.10f,LPm=%0.10f)\n",
			 log(max(1e-300,LRcum+LRbiasLS))+LS-2.0*log(Rtheta),WTm,logPV, LRcum, LRbiasLS, LS, log(Rtheta), LPm);
		fflush(stdout);
	      }
	    }
	  }
	} /* if(bestI >= 1) */
      } // if(VERB... || logPV || poutlier)
      
      if(DEBUG)assert(isfinite(LRcum));
      if(DEBUG) assert(isfinite(bestLR) && (LRcum <= 0.0 || bestLR > 0.0));
      if(DEBUG && !(isfinite(LLRsum))){
	printf("m=%d/%d:mapid=%d,id=%lld,flip=%d:LLR=%0.6f(delta=%0.6f,WT=%0.6f,LS=%0.6f,LRcum=%0.6e,LRbias=%0.6e),bestLR/LRcum=%0.8f,M=%d(%d),N=%d,best I=%d,K=%d,J=%d,Y[I]=%0.3f,X[J]=%0.3f",
	       m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,pcontig->flip[m],LLRsum,LPm,WTm,LS,LRcum,LRbias,bestLR/LRcum,M,MX[m],N,bestI,bestK,bestJ,
	       bestI>=1 ? Y[bestI] : -1.0, bestI>=1 ? Xm[bestJ] : -1.0);
	printf("bestLR=%0.10e,LRcum=%0.10e\n",bestLR,LRcum);
	fflush(stdout);
	assert(isfinite(LLRsum));
      }

      if(nmapK)	{ /* backtrack through array starting from bestI,bestK,bestJ to compute nmap[m] */
	if(DEBUG) assert(nmapK && nmap[m] && nmapK[m]);

	for(register int J=0; J <= M+1; J++)
	  nmap[m][J] = -1;/* default -1 :  unaligned */
	nmapK[m][0] = nmapK[m][M+1] = 0;/* default for unaligned : Lij = Rij = 0 */

	if(bestI>=1){
	  /* allocate space to store alignment in reverse order */
	  register int *Ilist = new int[3*M];
	  register int *Klist = &Ilist[M];
	  register int *Jlist = &Ilist[2*M];
	  register int U = 0;

	  register int I = bestI;
	  register int K = bestK;
	  register int J = bestJ;

	  while(I > 0){
	    if(DEBUG){
	      if(DEBUG && !(I-K >= 1 && I <= N)){
		printf("m=%d:U=%d:I=%d,K=%d,J=%d,N=%d,M=%d\n",
		       m,U,I,K,J,N,M);
		fflush(stdout);
		assert(I-K >= 1 && I <= N);
	      }
	      if(!(ISLINEAR && lc)) assert(J>=1);
	      if(!(ISLINEAR && rc)) assert(J<=M);
	    }
	    Ilist[U] = I;
	    Klist[U] = K;
	    Jlist[U] = J;
	    /*	  if(m==1){
		  printf("m=%d:U=%d:I=%d,K=%d,J=%d:LR=%0.8e,G=%d,T=%d,H=%d\n",m,U,I,K,J,p->LR,p->G,p->T,p->H);
		  fflush(stdout);
		  }*/
	    U++;
	    // RecA *p = &A[I][K][J];
	    int G = A.G(I,K,J);
	    int T = A.T(I,K,J);
	    int H = A.H(I,K,J);
	    if(DEBUG && !(G < 0 || T <= A.kmax)){
	      printf("m=%d:U=%d:I=%d,K=%d,J=%d,N=%d,M=%d,T=%d,KMAX=%d,A.kmax=%d,kmax=%d\n",
		     m,U,I,K,J,N,M,T,KMAX,A.kmax,kmax);
	      printf("bestI=%d,bestK=%d,bestJ=%d\n",bestI,bestK,bestJ);
	      fflush(stdout);
	      assert(G < 0 || T <= A.kmax);
	    }

	    I = G;
	    K = T;
	    J = H;
	  }
	  if(DEBUG) assert(U > 0);
	  if(DEBUG) assert(I <= (lc ? 0 : -1));
	  if(DEBUG) assert(bestR <= (rc ? 0 : -1));

	  if(I==0) {/* left chromosome end */
	    if(DEBUG) assert(lc && J==0);
	    nmap[m][0] = nmapK[m][0] = 0;
	  } else { /* compute nmapK[m][0] = Lij */
	    nmapK[m][0] = A.Lij(Ilist[U-1],Klist[U-1],Jlist[U-1]);
	    
	    if(poutlier && I <= -2 && /* WAS nmapK[m][0] == Ilist[U-1]-Klist[U-1] && */ Ilist[U-1]-Klist[U-1] > 1 && Xm[Jlist[U-1]] > Yc(Y,Ilist[U-1],Klist[U-1]) - Y[Ilist[U-1]-Klist[U-1]-1])
	      poutlier[m][0] |= 2;
	  }
	  if(VERB/* HERE >=2 */ && rverb && (rverb>=2 || gmap[pcontig->contig[m].mapid]->id == TRACE_ID)){
	    printf("m=%d/%d:mapid=%d,id=%lld:N=%d,M=%d,I=%d,K=%d,J=%d:Lij=%d,nmapK[m][0]=%d\n",
		   m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,N,M,Ilist[U-1],Klist[U-1],Jlist[U-1],A.Lij(Ilist[U-1],Klist[U-1],Jlist[U-1]),nmapK[m][0]);
	    fflush(stdout);
	  }

	  int D=0,H=0,G=0;
	  for(int T = 0; T < U; T++, D = K, H = I, G = J){
	    I = Ilist[U-1-T];
	    J = Jlist[U-1-T];
	    K = Klist[U-1-T];
	    if(DEBUG && !lc && !(J>0 && I>0)){
	      printf("m=%d:T=%d/%d:I=%d,K=%d,J=%d,Y[I-K..I]=Hcuts[%d..%d]=%0.3f..%0.3f(%0.3f),X[J]=%0.3f:",
		     m, T, U, I,K,J,nmap[MD][I-K],nmap[MD][I],Y[I-K],Y[I],Yc(Y,I,K),Xm[J]);
	      if(T<=0)
		printf("\n");
	      else
		printf("X[J]-X[G]=%0.3f,Y[I,K]-Y[H,D]=%0.3f\n",
		       Xm[J]-Xm[G],Yc(Y,I,K)-Yc(Y,H,D));
	      fflush(stdout);
	      assert(J>0 && I>0);
	    }
	    if(DEBUG && !(ISLINEAR && rc)) assert(J<=M && I <= N);
	    nmap[m][J] = I;
	    nmapK[m][J] = K;
	    if(VERB && rverb && (rverb >= 2 || gmap[pcontig->contig[m].mapid]->id == TRACE_ID)){
	      printf("m=%d:T=%d/%d:I=%d,K=%d,J=%d,Y[I-K..I]=Hcuts[%d..%d]=%0.3f..%0.3f(%0.3f),X[J]=%0.3f:",
		     m, T, U, I,K,J,nmap[MD][I-K],nmap[MD][I],Y[I-K],Y[I],Yc(Y,I,K),Xm[J]);
	      if(T<=0)
		printf("\n");
	      else
		printf("X[J]-X[G]=%0.3f,Y[I,K]-Y[H,D]=%0.3f\n",
		       Xm[J]-Xm[G],Yc(Y,I,K)-Yc(Y,H,D));
	    }
	  }
	  if(bestR==0)
	    nmap[m][M+1] = nmapK[m][M+1] = N+1;
	  else {/* compute nmapK[m][M+1] = Rij */
	    nmapK[m][M+1] = A.Rij(bestI,bestK,bestJ);
	    if(poutlier && bestR <= -2 && /* WAS nmapK[m][M+1] == bestI && */ bestI < N && Xm[M+1] - Xm[bestJ] > Y[bestI+1] - Yc(Y,bestI,bestK))
	      poutlier[m][M] |= 2;
	  }
	  if(VERB/* HERE >=2 */ && rverb && (rverb>=2 || gmap[pcontig->contig[m].mapid]->id == TRACE_ID)){
	    printf("m=%d/%d:mapid=%d,id=%lld:N=%d,M=%d,bestI=%d,bestK=%d,bestJ=%d,bestR=%d:Rij=%d,nmapK[m][M+1]=%d\n",
		   m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,N,M,bestI,bestK,bestJ,bestR,A.Rij(bestI,bestK,bestJ),nmapK[m][M+1]);
	    fflush(stdout);
	  }

	  if(DEBUG){
	    register int Lij = OUTLIER_FIX ? nmapK[m][0] : (nmapK[m][0] < 0 ? -1-nmapK[m][0] : nmapK[m][0]);
	    register int Rij = OUTLIER_FIX ? nmapK[m][M+1] : abs(nmapK[m][M+1]);
	    assert(0 <= Lij && Lij <= Rij && Rij <= N+1);
	  }
	  if(VERB>=2)
	    fflush(stdout);
	  delete [] Ilist;
	}
      } // if (nmapK)

      if(TIME_VERB && tverb>=2){
	double wt8 = wtime();
	#pragma omp critical
	{
	  printf("\t m=%d/%d:tid=%d,lp=%d,sk=%d,set=%d/%d,gl=%d:LPm=%0.4f,WT=%0.3f(LS=%0.4f,LRcum=%0.3e,LRbias=%0.3e),TBmapWT=%0.6f,M=%d,best I=%d,K=%d,J=%d,IMIN=%d,IMAX=%d,N=%d,stride=%lld/%lld->%lld\n",
		 m,MD,tid,myLoopCnt,myskipcnt,myqsetcnt,mysetcnt,myglobalcnt,LPm,WTm,LS,LRcum,LRbias,TBmapWT[m],M,bestI,bestK,bestJ,IMIN,IMAX,N,stride,origStrideMax, myStrideMax);
	  if(TIME_VERB)
	    printf("\t\t wall= %0.6f + %0.6f + %0.6f (wt3= %0.6f, wt4= %0.6f, wt5= %0.6f, wt5b= %0.6f, wt6= %0.6f, wt7= %0.6f, wt8= %0.6f), cum= %0.6f\n", 
		  wt1-wt,wt2-wt1,wt8-wt2, wt3-wt2, wt4-wt3, wt5-wt4, wt5b-wt5, wt6-wt5b, wt7-wt6, wt8-wt7, wt8);
	  else
	    printf("\t\t wall= %0.6f + %0.6f + %0.6f, cum= %0.6f\n", wt1-wt,wt2-wt1,wt8-wt2, wt8);
	  fflush(stdout);
	}
	if(TIME_VERB){
	  mywt1 += wt1-wt;
	  mywt2 += wt2-wt1;
	  mywt3 += wt3-wt2;
	  mywt4 += wt4-wt3;
	  mywt5 += wt5-wt4;
	  mywt5b += wt5b-wt5;
	  mywt6 += wt6 - wt5b;
	  mywt7 += wt7 - wt6;
 	  mywt8 += wt8 - wt7;
	}

	wt = wt8;
      }

      STOP(qprobeval_for); 
    } // omp for m = 0 .. MD-1
    STOP(qprobeval_for_wait);

    /* free up memory that is local to each thread (excluding rs_heap allocations) */
    if(DIAGONAL && MIN_MEM && Fmem){
#if USE_MMAP
      free_huge_pages(Fmem);
#else
      free(Fmem);
#endif // USE_MMAP==0
    }

    #pragma omp critical
    {
      if(OMP_DEBUG) LoopCnt += myLoopCnt;
      if(TIME_VERB || CNT_VERB){
	setcnt += mysetcnt;
	qsetcnt += myqsetcnt;
	globalcnt += myglobalcnt;
	skipcnt += myskipcnt;
      }
      if(TIME_VERB && tverb>=2){
        double wt1 = wtime();
	printf("\t tid=%d:loop=%d,skip=%d,set=%d/%d,gl=%d:wall=%0.6f secs (cum CPU= %0.6f, wall= %0.6f)\n",tid,myLoopCnt,myskipcnt,myqsetcnt,mysetcnt,myglobalcnt,wt1-wt0,mtime(),wt1);
	printf("\t\t wall= %0.6f + %0.6f + %0.6f (wt3= %0.6f, wt4= %0.6f, wt5= %0.6f, wt5b= %0.6f, wt6= %0.6f, wt7= %0.6f, wt8= %0.6f)\n",
          mywt1,mywt2,mywt3+mywt4+mywt5+mywt5b+mywt6+mywt7+mywt8,mywt3,mywt4,mywt5,mywt5b,mywt6,mywt7,mywt8);
	fflush(stdout);

        cumwt1 += mywt1;
        cumwt2 += mywt2;
	cumwt3 += mywt3;
	cumwt4 += mywt4;
	cumwt5 += mywt5;
	cumwt5b += mywt5b;
	cumwt6 += mywt6;
	cumwt7 += mywt7;
        cumwt8 += mywt8;
      }
    }
    STOP(qprobeval_parallel);
  } // pragma omp parallel
  STOP(qprobeval_parallel_wait);
  if(OMP_DEBUG) assert(LoopCnt == MD);// verify OMP loop worked

  if(TWO_LOOPS && TwoLoops) delete [] mlist;

  if(newLP){/* recompute LLRsum = sum(m=0..MD-1) newLP[m] * TBmapWT[m] so that is does not depend on thread ordering */
    double origLLRsum = LLRsum;
    LLRsum = 0.0;
    if(TBmapWT)
      for(int m = 0; m < MD; m++){
	LLRsum += newLP[m] * TBmapWT[m];
	if(DEBUG>=1+RELEASE/* HERE HERE >=2 */ && !isfinite(LLRsum)){
	  printf("qprobeval:m=%d/%d:newLP[m]= %0.8e, TBmapWT[m]= %0.8e, LLRsum= %0.8e\n",m,MD,newLP[m],TBmapWT[m],LLRsum);
	  fflush(stdout);
	  assert(isfinite(LLRsum));
	}
      }
    else
      for(int m = 0; m < MD; m++){
	LLRsum += newLP[m];
	if(DEBUG>=1+RELEASE/* HERE HERE >=2 */ && !isfinite(LLRsum)){
	  printf("qprobeval:m=%d/%d:newLP[m]= %0.8e, TBmapWT[m]= %0.8e, LLRsum= %0.8e\n",m,MD,newLP[m],TBmapWT[m],LLRsum);
	  fflush(stdout);
	  assert(isfinite(LLRsum));
	}
      }
    if((QDEBUG || rverb) && !(fabs(origLLRsum - LLRsum) < 1e-7)){
      printf("WARNING:origLLRsum=%0.10f,LLRsum=%0.10f(err=%0.10f)\n",origLLRsum,LLRsum,LLRsum-origLLRsum);
      fflush(stdout);
      if(QDEBUG) assert(fabs(origLLRsum - LLRsum) < 1e-6);
    }
  }

  LLRsum -= N*SITE_PEN;

  if((TIME_VERB || CNT_VERB) && VERB>=2){
    printf("qprobeval:skip=%d/%d/%d/%d,gl=%d\n",MD-setcnt-qsetcnt,MD-setcnt,skipcnt,MD,globalcnt);
    fflush(stdout);
  }

  if((TIME_VERB || CNT_VERB) && (VERB>=2 || tverb || rverb || (DEBUG && !isfinite(LLRsum)))){
    if(CNT_VERB)
      printf("qprobeval:LR=%0.8f: skip=%d/%d/%d/%d,gl=%d, V3D_cnt=%lu/%lu, isfinite(LR)=%d): cum wall time= %0.6f\n",
	     LLRsum,MD-setcnt-qsetcnt,MD-setcnt,skipcnt,MD,globalcnt,V3D_rcnt,V3D_totcnt,isfinite(LLRsum), wtime());
    else
      printf("qprobeval:LR=%0.8f: MD=%d, V3D_cnt=%lu/%lu, isfinite(LR)=%d): cum wall time= %0.6f\n",
	     LLRsum,MD,V3D_rcnt,V3D_totcnt,isfinite(LLRsum), wtime());

    if(TIME_VERB && tverb)
      printf("\t\t wall= %0.6f + %0.6f + %0.6f (wt3= %0.6f, wt4= %0.6f, wt5= %0.6f, wt5b= %0.6f, wt6= %0.6f, wt7= %0.6f, wt8= %0.6f)\n",
         cumwt1,cumwt2,cumwt3+cumwt4+cumwt5+cumwt5b+cumwt6+cumwt7+cumwt8,cumwt3,cumwt4,cumwt5,cumwt5b,cumwt6,cumwt7,cumwt8);
    fflush(stdout);
  }

  if(DEBUG) assert(isfinite(LLRsum));
 
  STOP(qprobeval);
  return LLRsum;
}

} // namespace probeval
