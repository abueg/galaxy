#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <float.h>

#include "constants.h"
#include "parameters.h"
#include "Ccontig.h"
#include "globals.h"
#include "timers.h"
#include "probeval.h"

static Ident mprobeval_cpp_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/mprobeval.cpp 11397 2020-08-03 00:37:11Z tanantharaman $");

// #undef DEBUG
// #define DEBUG 2

#define MISSCNT_FIX 1 /* Make Lij,Rij follow convention in refalign, so they correctly call actual overlap endpoint (rather than Viterbi coverage endpoint) */

#undef DIAGONAL
#define DIAGONAL 1 // Reduce memory used with FLAT_MEMORY

#undef MIN_MEM
#define MIN_MEM MP_MIN_MEM // Minimize virtual memory used with FLAT_MEMORY (larger values reduce memory use but call malloc/free more often

#undef CACHE
#define CACHE 1  // improve cache usage of memory by avoiding changing start of arrays (slight increase in total memory use)
// CACHE uses ALIGN and ALIGN2 defined in probeval.h

#define LDEBUG ((MDEBUG && MDEBUG_M>=0) ? 1 : 0) /* 1 : debugging code (slower, but better validation when using MDEBUG) (2: forces arrays sizes to 0..N+1)*/
                                                        // NOTE : LDEBUG == 2 appears to be broken and triggers false alarms

#define LRSITE_ERR 1e-4 // 1e-8 /* relative error in LRsite[1..S]  relative to LRcum that triggers an assertion */
#define LRFRAG_ERR 1e-8 /* relative error in LRfrag[1..N+1]  relative to LRcum that triggers an assertion */

#define FIXLRDEL 1 /* force computation of LRdel[0] (required since LRdel[0]==LRdel[1] uses a reduced deltaY) */

#define FIXADD /* WAS 1  */ (addcnt ? 1 : 0) /* fix alignment span to account for site added */

#define TDEBUG 0 /* debug thread order independent LP values by comparing with original method (requires additional critical sections) */

#ifdef WIN32
    #ifndef nan
        static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
        #define nan() (*(const double *) __nan)
    #endif
  #define FL128 long double
  #define FL128b long double
#else // NOT WIN32
#ifdef __INTEL_COMPILER
  #define FL128 _Quad
  #define FL128b long double /* 80 bit float */
#else
#ifdef CLANG // Clang compiler does not implement 128bit float
  // HERE HERE HERE : consider using an extended precision library like f107_o.cpp,f107_o.h or qd-2.3.22
  #define FL128 long double /* 80 bit float */
  #define FL128b long double /* 80 bit float */
#else
#ifdef __ARM_ARCH // gcc-arm does not implement float128 or 80-bit long double
  // HERE HERE HERE : consider using an extended precision library like f107_o.cpp,f107_o.h or qd-2.3.22
  #define FL128 long double
  #define FL128b long double
#else // GCC based X86
  #include <quadmath.h>
  #define FL128 __float128
  #define FL128b long double /* 80 bit float */
#endif // NOT __ARM_ARCH
#endif // NOT CLANG
#endif // NOT __INTEL_COMPILER
#endif // NOT WIN32

#if !defined(WIN32) && !defined(CLANG) && !defined(__ARM_ARCH)
  static bool isfinite(FL128 x)
  {
    return isfinite((long double)x);
  }
#endif

#define VERB_EPS2 0 // verbose display of savings due to EPS2
#define CRITLOOP_MERGE 1 /* merge 2 critical loops (when OutlierType==1) into 1 loop */

#define CUMERR (sizeof(FL128)==8 ? 2e-9 : 1e-14)
#define PQ_BND 0 // Enable bounding of P and Q during 3-D recurrance based on EPS3 (Only implemented with FLAT_MEMORY==1)
#define EPS_ADJUST 1 // Adjust final LP values by drop caused by EPS2 and EPS3
// HERE : undo EPS2 & EPS3 range reduction by 1 or 2 labels to account of changes in interval sizes

namespace probeval {

class QB {
public:
  int Kstride;/**< distance between K index increments == maxN*maxM */
  int Istride;/**< distance between I index increments <= M */
  int kmax;/**< largest value of Kmax[1..N] */
  /* NOTE : the J index stride is always 1 */

#define ARRAY(type, field) type *field##_base; inline type& field(int I, int K, int J) { return field##_base[I*Istride + K*Kstride + J];}

  ARRAY(MFLOAT2, UL);
  ARRAY(MFLOAT2, UR);
  ARRAY(MFLOAT2, bUL);
  ARRAY(MFLOAT2, bUR);
  
  ARRAY(MFLOAT2, AL);
  ARRAY(MFLOAT2, AR);

#if PQ_BND
  // Pmin,Pmax,Qmin,Qmax needed to keep track of next site range in best alignment (computed during 3-D recurrance of AR), to better handle large DELTA_X,DELTA_Y values
  ARRAY(int, Pmin);
  ARRAY(int, Pmax);
  ARRAY(int, Qmin);
  ARRAY(int, Qmax);
#endif

  ARRAY(int,ULmin);
  ARRAY(int,ULmax);
  ARRAY(int,URmin);
  ARRAY(int,URmax);
  ARRAY(int,Lij);
  ARRAY(int,Rij);
#undef ARRAY
  
  QB() {
    UL_base = UR_base = bUL_base = bUR_base = AL_base = AR_base = NULL;
#if PQ_BND >= 1
    Pmin_base = Pmax_base = Qmin_base = Qmax_base = NULL;
#endif
    ULmin_base = ULmax_base = URmin_base = URmax_base = Lij_base = Rij_base = NULL;
  };
};

static void LPdeltaTrace(QB &A, int m, int I, int K, int J, int P, int T, int QminP, int QmaxP, int N, int S, int d, int *Jmin, FLOAT *Xm, double ALscale, double *LogScaleL, double *LogScaleR, double LS, double del, MFLOAT deltaY, MFLOAT delY, MFLOAT Ivar, MFLOAT FAikD, double LRiptD128, int Amin, int Amax, double **LRdeltaD, double LRdeltaDSd, FL128 *LRdelta, const int OUT_TYPE, const int OUT_OLD)
{
  if(MDEBUG_M >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S == MDEBUG_TS && d == MDEBUG_TT && LRiptD128 >= MDEBUG_MIN){
    FL128 sum = LRdelta[0];
    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
      sum += LRdelta[i];
    FL128 sum3 = sum + LRdeltaD[MDEBUG_TS][MDEBUG_TT] + LRdeltaDSd;
    printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d:S=%d,d=%d(%0.3f),deltaY=%0.3f->%0.3f:FAikD=%0.10e,sum(AL*FA*AR)=%0.10e,LRdelta[0..S]=%0.10e,LRdeltaD[S][d]=%0.10e(sum=%0.10e)\n",
	   m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxP,S,d,del,deltaY,delY,FAikD,(double)LRiptD128,(double)sum,(double)LRdeltaD[S][d], (double)sum3);
    double tsum = 0.0;
    for(int Q = QminP; Q <= QmaxP; Q++){
      double FAikjD =   FAJ(Xm[Q]-Xm[J], delY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], FAikD,Ivar,1, OUT_TYPE, OUT_OLD);
      tsum += FAikjD * A.AR(P,T,Q);
      printf("    Q=%d:FAikjD=%0.10e,AR=%0.10e,ALscale=%0.10e(LS=%0.10f), sum(AL*FA*AR)=%0.10e(delta= %0.10e)\n",
	     Q, FAikjD, A.AR(P,T,Q),(double)ALscale,LogScaleL[I]+LogScaleR[P]-LS,(double)(tsum*ALscale), (double)(FAikjD * A.AR(P,T,Q) * ALscale));
      }
  }
}

template<int OUT_TYPE, int OUT_OLD>
static __attribute__ ((noinline)) void mprobeval_delta(QB &A, 
						       int m, int I, int K, int J, int N, int *D,
						       MFLOAT **delta, double **LRdeltaD,
						       int IMAX, int *Kmax, int *Jmin, int *Jmax, int Pmax, int Pmax2, int Qmin, int Qmax, int Qmax2, 
						       int Amin, int Amax,
						       double *Y, double *Xm, double Yik, 
						       double *LogScaleL, double *LogScaleR, double LS,
						       double AL, 
						       FL128 *LRdelta)
{
  // New loop order S,d,P,T,Q : allows larger vectors to be filled

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT
  __m256 v_resL2;
  if(RES_VARIANCE)
    v_resL2 = _mm256_broadcast_ss(&VARtab[K*N + I]);
#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT
  __m512 v_resL2;
  if(RES_VARIANCE)
    v_resL2 = _mm512_set1_ps(VARtab[K*N + I]);
#endif

  for(int S = I; S < Pmax; S++){
    int d = 0;
#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT
    /* Use AVX256 intrinsincs (manual optimization) */

    for(; d < D[S] - VEC_MSIZ; d += VEC_SIZ){/* NOTE : vectorization is useful even if only 3 out of 8 entries of each vector will be used */
      /* NOTE : delta[S] and LRdeltaD[S] are padded, if needed, to a multiple of 8 entries : padded entries are initialized to zero */
      __m256 v_del = _mm256_load_ps(&delta[S][d]);
      __m256d v_LRdeltaDS1 = _mm256_load_pd(&LRdeltaD[S][d]);
      __m256d v_LRdeltaDS2 = _mm256_load_pd(&LRdeltaD[S][d+4]);
      
      for(int PT = S+1; PT <= Pmax; PT++){// NOTE : PT == P - T
	for(int P, T = 0; (P = PT+T) <= IMAX && T <= Kmax[P]; T++){
	  MFLOAT Ypti = Y[P-T] - Y[I];
	  __m256 v_Ypti = _mm256_set1_ps(Ypti);
	  double ALscale = AL * exp(LogScaleL[I] + LogScaleR[P] - LS);
	  int QminP = max(Qmin, Jmin[P]);
	  int QmaxP = (PT <= Pmax2) ? min(Qmax, Jmax[P]) : min(Qmax2, Jmax[P]);
	  MFLOAT deltaY = Yc(Y,P,T) - Yik;
	  __m256 v_deltaY = _mm256_set1_ps(deltaY);// WAS broadcast_ss(&deltaY);
	  __m256 v_delY = _mm256_add_ps(v_deltaY, v_del);
	  v_Ypti = _mm256_add_ps(v_Ypti, v_del);
	  __m256 v_Ivar;
	  __m256 v_FAikD = FAIKD_mm256(v_delY, VARtab[T*N + P], v_resL2, v_Ypti, v_Ivar);
	  __m256 v_LogFnPown1 = _mm256_broadcast_ss(&LogFnPow[P-T-I-1]);

	  float PrBiasWTpti = PrBiasWT[P-T] - PrBiasWT[I];
	  __m256 v_PrBiasWT = _mm256_set1_ps(PrBiasWTpti);

	  MFLOAT2 *ARpt = &A.AR(P,T,0);
	  __m256 v_LRipD = _mm256_setzero_ps();

	  for(int Q = QminP; Q <= QmaxP; Q++){
	    __m256 v_ARptq = _mm256_broadcast_ss(&ARpt[Q]);
	    __m256 v_Faj = FAJ_mm256<1, OUT_TYPE, OUT_OLD>(Xm[Q]-Xm[J], v_delY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], v_LogFnPown1, v_PrBiasWT, v_FAikD, v_Ivar);
	    v_LRipD = _mm256_add_ps(v_LRipD, _mm256_mul_ps(v_Faj, v_ARptq));
	  }
	  
	  /* LRdeltaDSd += LRipD * ALscale */
	  __m256d v_ALscale = _mm256_set1_pd(ALscale);
	  __m256d v_LRipD1 = _mm256_cvtps_pd(_mm256_extractf128_ps(v_LRipD,0));/* low order 4 floats to double */
	  __m256d v_LRipD2 = _mm256_cvtps_pd(_mm256_extractf128_ps(v_LRipD,1));/* high order 4 floats to double */

#if VECTOR_DEBUG
	  float LRipDV[8]; _mm256_storeu_ps(LRipDV, v_LRipD);
	  double LRipD1V[4]; _mm256_storeu_pd(LRipD1V, v_LRipD1);
	  double LRipD2V[4]; _mm256_storeu_pd(LRipD2V, v_LRipD2);
	  for(int i = 0; i < 4; i++){
	    if(!(fabs(LRipD1V[i] - LRipDV[i]) <= 1e-8 * LRipDV[i])){
	      printf("i=%d:LRipD1V[i]= %0.8e, LRipDV[i]= %0.8e\n",i, LRipD1V[i], LRipDV[i]);
	      fflush(stdout);
	      assert(fabs(LRipD1V[i] - LRipDV[i]) <= 1e-8 * LRipDV[i]);
	    }
	    if(!(fabs(LRipD2V[i] - LRipDV[i+4]) <= 1e-8 * LRipDV[i+4])){
	      printf("i=%d:LRipD2V[i]= %0.8e, LRipDV[i+4]= %0.8e\n",i, LRipD2V[i], LRipDV[i+4]);
	      fflush(stdout);
	      assert(fabs(LRipD2V[i] - LRipDV[i+4]) <= 1e-8 * LRipDV[i+4]);
	    }
	  }
#endif

	  v_LRdeltaDS1 = _mm256_add_pd(v_LRdeltaDS1, _mm256_mul_pd(v_LRipD1, v_ALscale));
	  v_LRdeltaDS2 = _mm256_add_pd(v_LRdeltaDS2, _mm256_mul_pd(v_LRipD2, v_ALscale));
	}// T = 0 .. Kmax[P]
      }// PT = S+1 .. Pmax

      _mm256_store_pd(&LRdeltaD[S][d], v_LRdeltaDS1);
      _mm256_store_pd(&LRdeltaD[S][d+4], v_LRdeltaDS2);
    }

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT==1 && !(DEBUG>=2 || (MDEBUG_M >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0))
    /* Use MIC intrinsincs (manual optimization) */

    for(; d < D[S] - VEC_MSIZ; d += VEC_SIZ){/* NOTE : vectorization is useful even if only 3 out of 16 entries of each vector will be used */
      /* NOTE : delta[S] and LRdeltaD[S] are padded, if needed, to a multiple of 16 entries : padded entries are initialized to zero */
      __m512 v_del = _mm512_load_ps(&delta[S][d]);
      __m512d v_LRdeltaDS1 = _mm512_load_pd(&LRdeltaD[S][d]);
      __m512d v_LRdeltaDS2 = _mm512_load_pd(&LRdeltaD[S][d+8]);
      
      for(int PT = S+1; PT <= Pmax; PT++){// NOTE : PT == P - T
	for(int P, T = 0; (P = PT+T) <= IMAX && T <= Kmax[P]; T++){
	  MFLOAT Ypti = Y[P-T] - Y[I];
	  __m512 v_Ypti = _mm512_set1_ps(Ypti);
	  double ALscale = AL * exp(LogScaleL[I] + LogScaleR[P] - LS);
	  int QminP = max(Qmin, Jmin[P]);
	  int QmaxP = (PT <= Pmax2) ? min(Qmax, Jmax[P]) : min(Qmax2, Jmax[P]);
	  MFLOAT deltaY = Yc(Y,P,T) - Yik;
	  __m512 v_deltaY = _mm512_set1_ps(deltaY);
	  __m512 v_delY = _mm512_add_ps(v_deltaY, v_del);
	  v_Ypti = _mm512_add_ps(v_Ypti, v_del);
	  __m512 v_Ivar;
	  __m512 v_FAikD = FAIKD_mm512(v_delY, VARtab[T*N + P], v_resL2, v_Ypti, v_Ivar);

	  __m512 v_LogFnPown1 = _mm512_set1_ps(LogFnPow[P-T-I-1]);

	  MFLOAT PrBiasWTpti = PrBiasWT[P-T] - PrBiasWT[I];
	  __m512 v_PrBiasWTpti = _mm512_set1_ps(PrBiasWTpti);

	  MFLOAT2 *ARpt = &A.AR(P,T,0);
	  __m512 v_LRipD = _mm512_setzero_ps();

	  for(int Q = QminP; Q <= QmaxP; Q++){
	    __m512 v_ARptq = _mm512_set1_ps(ARpt[Q]);
	    __m512 v_Faj = FAJ_mm512<1, OUT_TYPE, OUT_OLD>(Xm[Q]-Xm[J], v_delY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], v_LogFnPown1, v_PrBiasWTpti, v_FAikD, v_Ivar);
	    v_LRipD = _mm512_add_ps(v_LRipD, _mm512_mul_ps(v_Faj, v_ARptq));
	  }
	  
	  /* LRdeltaDSd += LRipD * ALscale */
	  __m512d v_ALscale = _mm512_set1_pd(ALscale);
#if USE_MIC
	  __m512d v_LRipD1 = _mm512_cvtpslo_pd(v_LRipD);/* low order 8 floats to double */
	  __m512d v_LRipD2 = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_LRipD, 0x2 | (0x3 << 2)));/* high order 8 floats to double */
#elif USE_AVX512
	  __m512d v_LRipD1 = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_castsi512_ps(_mm512_castsi256_si512(_mm512_extracti64x4_epi64(_mm512_castps_si512(v_LRipD),0)))));/* low order 8 floats to double */
	  __m512d v_LRipD2 = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_castsi512_ps(_mm512_castsi256_si512(_mm512_extracti64x4_epi64(_mm512_castps_si512(v_LRipD),1)))));/* high order 8 floats to double */
#endif

#if VECTOR_DEBUG
	  float LRipDV[16]; _mm512_storeu_ps(LRipDV, v_LRipD);
	  double LRipD1V[8]; _mm512_storeu_pd(LRipD1V, v_LRipD1);
	  double LRipD2V[8]; _mm512_storeu_pd(LRipD2V, v_LRipD2);
	  for(int i = 0; i < 8; i++){
	    if(!(fabs(LRipD1V[i] - LRipDV[i]) <= 1e-8 * LRipDV[i])){
	      printf("i=%d:LRipD1V[i]= %0.8e, LRipDV[i]= %0.8e\n",i, LRipD1V[i], LRipDV[i]);
	      fflush(stdout);
	      assert(fabs(LRipD1V[i] - LRipDV[i]) <= 1e-8 * LRipDV[i]);
	    }
	    if(!(fabs(LRipD2V[i] - LRipDV[i+8]) <= 1e-8 * LRipDV[i+8])){
	      printf("i=%d:LRipD2V[i]= %0.8e, LRipDV[i+8]= %0.8e\n",i, LRipD2V[i], LRipDV[i+8]);
	      fflush(stdout);
	      assert(fabs(LRipD2V[i] - LRipDV[i+8]) <= 1e-8 * LRipDV[i+8]);
	    }
	  }
#endif

	  v_LRdeltaDS1 = _mm512_add_pd(v_LRdeltaDS1, _mm512_mul_pd(v_LRipD1, v_ALscale));
	  v_LRdeltaDS2 = _mm512_add_pd(v_LRdeltaDS2, _mm512_mul_pd(v_LRipD2, v_ALscale));
	}// T = 0 .. Kmax[P]
      }// PT = S+1 .. Pmax

      _mm512_store_pd(&LRdeltaD[S][d], v_LRdeltaDS1);
      _mm512_store_pd(&LRdeltaD[S][d+8], v_LRdeltaDS2);
    }

#else // neither VECTOR_AVX nor VECTOR_MIC

    if(0 && D[S] >= (USE_MIC ? 16 : 8)){/* rely on omp simd vectorization over d loop */
      MFLOAT *deltaS = delta[S];
      double *LRdeltaDS = LRdeltaD[S];

#if USE_MIC
      #pragma omp simd aligned(deltaS, LRdeltaDS : 64)
#else
      #pragma omp simd aligned(deltaS, LRdeltaDS : 32)
#endif
      for(int t= 0; t < D[S]; t++){
	MFLOAT del = deltaS[t];
	double LRdeltaDSd = 0.0;

	for(int PT = S+1; PT <= Pmax; PT++){// NOTE : PT == P - T
	  for(int P, T = 0; (P = PT+T) <= IMAX && T <= Kmax[P]; T++){
	    MFLOAT Ypti = Y[P-T] - Y[I];
	    double ALscale = AL * exp(LogScaleL[I] + LogScaleR[P] - LS);
	    if(DEBUG>=2) assert(isfinite(ALscale));
	    int QminP = max(Qmin, Jmin[P]);
	    int QmaxP = (PT <= Pmax2) ? min(Qmax, Jmax[P]) : min(Qmax2, Jmax[P]);
	    MFLOAT deltaY = Yc(Y,P,T) - Yik;
	    MFLOAT delY = deltaY + del;
	    if(DEBUG>=2) assert(delY >= minKB - (USE_MFLOAT ? 1e-6 : 2e-9));
	    MFLOAT Ivar;
	    MFLOAT FAikD = FAIKD(delY,P-T-I,P,T,K,Ypti + del,Ivar,Y);
	    MFLOAT2 *ARpt = &A.AR(P,T,0);
	    MFLOAT2 LRiptD = 0.0;

	    for(int Q = QminP; Q <= QmaxP; Q++)
	      LRiptD += FAJ<1,OUT_TYPE,OUT_OLD>(Xm[Q]-Xm[J], delY, FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], FAikD, Ivar) * ARpt[Q];

	    double LRiptD128 = LRiptD * ALscale;
	    if(DEBUG>=2) assert(Amin <= S && S <= min(N,Amax));
	    LRdeltaDSd += LRiptD128;
	    if((MDEBUG_M >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0) || (DEBUG>=2 && !(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0))) 
	      LPdeltaTrace(A,m,I,K,J,P,T,QminP,QmaxP,N,S,d,Jmin,Xm,ALscale,LogScaleL,LogScaleR,LS,del,deltaY,delY,Ivar,FAikD,LRiptD128,Amin,Amax,LRdeltaD,LRdeltaDSd,LRdelta,OUT_TYPE,OUT_OLD);
	    if(DEBUG>=2) assert(isfinite(LRiptD128));
	    if(DEBUG>=2) assert(isfinite(LRdeltaDSd) && LRdeltaDSd >= 0.0);
	  }// T = 0 .. Kmax[P]
	}// PT = S+1 .. Pmax
	LRdeltaDS[t] += LRdeltaDSd;
	if(DEBUG>=2) assert(isfinite(LRdeltaDS[t]) && LRdeltaDS[t] >= 0.0);
      }// d = 0 .. D[S]-1

      d = D[S];
    }
#endif // neither VECTOR_AVX or VECTOR_MIC

    /* Default : rely on omp simd vectorization over innermost Q loop */
    for(; d < D[S]; d++){// original code : used to process remaining iterations of d
      MFLOAT del = delta[S][d];
      double LRdeltaDSd = 0.0;

      for(int PT = S+1; PT <= Pmax; PT++){// NOTE : PT == P - T
	for(int P, T = 0; (P = PT+T) <= IMAX && T <= Kmax[P]; T++){
	  MFLOAT Ypti = Y[P-T] - Y[I];
	  double ALscale = AL * exp(LogScaleL[I] + LogScaleR[P] - LS);
	  if(DEBUG>=2) assert(isfinite(ALscale));
	  int QminP = max(Qmin, Jmin[P]);
	  int QmaxP = (PT <= Pmax2) ? min(Qmax, Jmax[P]) : min(Qmax2, Jmax[P]);
	  MFLOAT deltaY = Yc(Y,P,T) - Yik;
	  MFLOAT delY = deltaY + del;
	  if(DEBUG>=2) assert(delY >= minKB - (USE_MFLOAT ? 1e-6 : 2e-9));
	  MFLOAT Ivar;
	  MFLOAT FAikD = FAIKD(delY,P-T-I,P,T,K,Ypti + del,Ivar,Y);
	  MFLOAT2 *ARpt = &A.AR(P,T,0);
	  MFLOAT2 LRiptD = 0.0;

	  #if !(((VECTOR_MIC && USE_MIC) || (VECTOR_AVX && USE_AVX && USE_SSE)) && USE_MFLOAT)
          #pragma omp simd reduction(+:LRiptD)
	  #endif
	  for(int Q = QminP; Q <= QmaxP; Q++)
	    LRiptD += FAJ(Xm[Q]-Xm[J],delY,Q-J,P-T-I,FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], FAikD,Ivar,1, OUT_TYPE, OUT_OLD) * ARpt[Q];

	  double LRiptD128 = LRiptD * ALscale;
	  if(DEBUG>=2) assert(Amin <= S && S <= min(N,Amax));
	  LRdeltaDSd += LRiptD128;
	  if((MDEBUG_M >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0) || (DEBUG>=2 && !(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0))) 
	    LPdeltaTrace(A,m,I,K,J,P,T,QminP,QmaxP,N,S,d,Jmin,Xm,ALscale,LogScaleL,LogScaleR,LS,del,deltaY,delY,Ivar,FAikD,LRiptD128,Amin,Amax,LRdeltaD,LRdeltaDSd,LRdelta,OUT_TYPE,OUT_OLD);
	  if(DEBUG>=2) assert(isfinite(LRiptD128));
	  if(DEBUG>=2) assert(isfinite(LRdeltaDSd) && LRdeltaDSd >= 0.0);
	}// T = 0 .. Kmax[P]
      }// PT = S+1 .. Pmax
      LRdeltaD[S][d] += LRdeltaDSd;
      if(DEBUG>=2) assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
    }// d = 0 .. D[S]-1
  }// S = I .. Pmax-1
}

static long long mDPsizF = 6;
static long long mDPsizI = 6 + 4*PQ_BND;

/** compute AL[I==S+K+1][K][J] (see line 6 of Equation 10 : does not include PM(I,K)) */
/** NOTE : When called to compute ALT[K][J]: Y[I] is new site and can be Y[N+1] */
static double AL(int m, /* for debugging only */
		 int I, 
		 int K,
		 int J,
		 int U, /**< initialized to Umin */
		 int Umax,/**< If Umax < 0, U range is unknown: IFF Umax== -1, Update A[I][K][J].ULmin,ULmax */
		 double Yik,
		 double *  Xm,
		 double *  Y,
		 int M,
		 int N,
		 int lc,
		 double Scale,/* exp(-LogScaleL[I]) */
		 double *  ScaleLI, /* ScaleLI[P] = exp(LogScaleL[min(max(P,IMIN-1),IMAX+1)] - LogScaleL[I]) */
		 double * LogScaleL,/* for debugging only */		 
		 QB &  A,
		 int IMIN,
		 int IMAX,
		 int *  Kmax,
		 int *  Jmin,
		 int *  Jmax,
		 const int OUT_TYPE,
		 const int OUT_OLD)
{
  int origU,origUmax;
  if(DEBUG>=2){
    origU = U;
    origUmax = Umax;
  }

  double UL = PoutlierEnd;
  int JMAX = (EXTEND && I-K <= DELTA_YU) ? (DELTA_XU + EXTEND_MAX) : DELTA_XU;
  if(DEBUG>=2) assert(J > 0);

  if(J <= JMAX){
    double IvarX;
    double FLj = FLJ(J,Xm[J],IvarX);
    double NotOutlier = FLIK_FIX ? PoutlierEndM1 : FLj * PoutlierEndM1;

    MFLOAT ULu;
    if(DEBUG>=2 && !(isfinite(FLj) && FLj >= 0.0 && FLj <= MaxUR * 3.0)){
      printf("AL(m=%d,I=%d,K=%d,J=%d,U=%d,%d):UL=%0.10e,IMIN=%d,IMAX=%d\n",
	     m,I,K,J,U,Umax,(double)UL,IMIN,IMAX);
      printf("FLJ(J,Xm[J],IvarX)=%0.10e:IvarX=%0.10e\n",FLj,IvarX);
      fflush(stdout);
      assert(isfinite(FLj) && FLj >= 0.0 && FLj <= MaxUR * 3.0);
    }

    if(DEBUG/* HERE >=2 */) assert(U > Umax);// NOTE : Ub code below is NOT correct, see ALq() to fix it if assertion triggers
    if(U <= Umax){
      if(ERFC_FIX && DEBUG /* HERE >=2 */) assert(IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]);
      if(DEBUG>=2) assert(U >= 0);
      if(DEBUG>=2) assert(Umax < I-K);
      int Umin = U;
      int Ub = max(1,A.Lij(I,K,J)) - 1;

      U = Ub;
      ULu = FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);

      for(int U = Ub; ++U <= Umax;)
	ULu += ERFC_FIX ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
	  : FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);

      for(int U = Ub; --U >= Umin;)
	ULu += FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);

      if(lc && Umin <= 0)
	ULu += FLEIK(I,K,J,0,1,Yik,Xm[J], FLj);
    } else {/* need to determine U range */
      int save = (Umax == -1) ? 1 : 0;
      int Lijk = (IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]) ? max(1,A.Lij(I,K,J)) : 0;
      if(DEBUG>=2)assert(Lijk <= I-K);
      /* increment Lijk until Y[Lijk] overlap Xm */
      while(Lijk < I-K && Yik - Y[Lijk] > Xm[J])
	Lijk++;
      /* decrement Lijk until Y[Lijk-1] no longer overlaps Xm */
      while(Lijk > 0 && Yik - Y[Lijk-1] <= Xm[J])
	Lijk--;
      int Umin, Ub;
      Ub = Umin = Umax = max(1,Lijk) - 1;
      if(DEBUG>=2 && !(0 <= Umin && Umin < I-K)){
	printf("AL:I=%d,K=%d,J=%d:Umin=%d,Lijk=%d,IMIN=%d,IMAX=%d",I,K,J,Umin,Lijk,IMIN,IMAX);
	if(IMIN <= I && I <= IMAX){
	  printf(",Kmax[I]=%d,Jmin[I]=%d,Jmax[I]=%d",Kmax[I],Jmin[I],Jmax[I]);
	  if(0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I])
	    printf(",A[I][K][J].Lij=%d",A.Lij(I,K,J));
	}
	printf("\n");
	fflush(stdout);
	assert(0 <= Umin && Umin < I-K);
      }
      ULu = FLIK<0>(I,K,J,Umin,Umin+1,Yik,Xm[J],Yik-Y[Umin],Y[Umin+1]-Y[Umin],Y,IvarX, FLj, END_FIX ? lc : 1);
      MFLOAT bUL = ULu;
      if(DEBUG>=2)assert(isfinite(ULu) && ULu >= 0.0);
      if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N+1 && N+1 <= MDEBUG_N2 && m==MDEBUG_M && K==MDEBUG_K && J==MDEBUG_J){
	printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e,FLIK=%0.10e,FLj=%0.10e,PoutlierEndM1=%0.10e,PM=%0.10e):Yik=%0.4f,X[J]=%0.3f,Y[I-K]=%0.3f,Y[U]=%0.3f,Y[U+1]=%0.3f\n",
	       m,I,K,J,Umin,N,(double)(UL+ULu*NotOutlier),ULu*NotOutlier,ULu,FLj,PoutlierEndM1,PM(0,I,K,Y),Yik,Xm[J],Y[I-K],Y[Umin],Y[Umin+1]);
	fflush(stdout);
      }
      for(U = Umax; ++U < I-K;){
	MFLOAT LR = ERFC_FIX ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
	  : FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);

	if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N+1 && N+1 <= MDEBUG_N2 && m==MDEBUG_M && K==MDEBUG_K && J==MDEBUG_J){
	  printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e,FLIK=%0.10e,FLj=%0.10e)\n",
		 m,I,K,J,U,N,(double)(UL+(ULu+LR)*NotOutlier),LR*NotOutlier, LR, FLj);
	  fflush(stdout);
	}
        if(LR < ULu * EPS)
	  break;
        if(URANGE_FIX && LR > bUL){
	  ULu = 0.0;
	  Umin = U;
	  bUL = LR;
        }
	Umax = U;
	ULu += LR;
	if(DEBUG>=2)assert(isfinite(ULu) && ULu >= 0.0);
      }
      for(U = Umin ;--U >= 0;){
	MFLOAT LR = (ERFC_FIX && U > Ub) ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
	  : FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);
	if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N+1 && N+1 <= MDEBUG_N2 && m==MDEBUG_M && K==MDEBUG_K && J==MDEBUG_J){
	  printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e,FLIK=%0.10e,FLj=%0.10e)\n",
		 m,I,K,J,U,N,(double)(UL+(ULu+LR)*NotOutlier),LR*NotOutlier,LR,FLj);
	  fflush(stdout);
	}
	if(LR < ULu * EPS)
	  break;
	Umin = U;
	ULu += LR;
	if(DEBUG>=2)assert(isfinite(ULu) && ULu >= 0.0);
      }

      if(ISLINEAR && lc){
	MFLOAT LR = FLEIK(I,K,J,0,1,Yik,Xm[J],FLj);
	if(LR >= ULu * EPS){
	  Umin = 0;
	  ULu += LR;
	  if(DEBUG>=2)assert(isfinite(ULu) && ULu >= 0.0);
	}
      }
      if(DEBUG>=2)assert(Umax < I-K);
      if(save){
	if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]);
	A.ULmin(I,K,J) = Umin;
	A.ULmax(I,K,J) = Umax;
      }
    }

    double origUL;
    if(DEBUG>=2) origUL = UL;
    UL += ULu * NotOutlier;
    if(DEBUG>=2 && !(isfinite(UL) && UL >= 0.0)){
      printf("AL(m=%d,I=%d,K=%d,J=%d,U=%d,%d):UL= %0.10e -> %0.10e,IMIN=%d,IMAX=%d\n",
	     m,I,K,J,U,Umax,(double)origUL,(double)UL,IMIN,IMAX);
      printf("PoutlierEndM1=%0.10e,ULu=%0.10e,FLj=%0.10e\n",
	     PoutlierEndM1,ULu,FLj);
      fflush(stdout);
      assert(isfinite(UL) && UL >= 0.0);
    }
  }
  if(DEBUG && !(isfinite(UL) && UL >= 0.0)){
    printf("AL(m=%d,I=%d,K=%d,J=%d,U=%d,%d):UL=%0.10e,IMIN=%d,IMAX=%d\n",
	   m,I,K,J,U,Umax,(double)UL,IMIN,IMAX);
    fflush(stdout);
    assert(isfinite(UL) && UL >= 0.0);
  }
  double AL = UL * Scale;
  if(DEBUG) assert(isfinite(AL) && AL >= 0.0);
  int Pmin = max(IMIN, I-K - deltaExtYRef);
  int Pmin2 = max(IMIN, I-K - DELTA_Y);
  int Qmin = max(1, J - deltaExtXRef);
  int Qmin2 = max(1, J - DELTA_X);

  if(DEBUG>=2 && !(I-K-1 <= IMAX)){
    printf("\n m=%d,I=%d,K=%d,J=%d:IMIN=%d,IMAX=%d,N=%d\n",m,I,K,J,IMIN,IMAX,N);
    fflush(stdout);
    assert(I-K-1 <= IMAX);
  }

  double PrBiasWTik;/* value of PrBiasWT[I-K] is not valid (if K==0), since Y[I] is a new label, but differentials PrBiasWT[t=0..I-1] remain valid */
  if(K)
    PrBiasWTik = PrBiasWT[I-K];
  else {/* compute correct value of PrBiasWT[I] */
    double y = Y[I] - Y[I-1];
    double pr = Pr(y);
    if(DEBUG>=2) assert(isfinite(pr) && 0.0 <= pr && pr <= 1.0);
    double prbias = (pr <= 0.0 || pr >= 1.0) ? 0.0 : pr * log(pr) + (1.0 - pr) * log(1.0 - pr);    
    if(DEBUG>=2)  assert(isfinite(prbias) && prbias <= 0.0);
    PrBiasWTik = PrBiasWT[I-1] - PRbiasWT * prbias;
  }

  for(int P = I-K; --P >= Pmin;){
    MFLOAT2 LRp = 0.0;
    MFLOAT PrBiasWTikp = PrBiasWTik - PrBiasWT[P];
    for(int T = 0; T <= Kmax[P]; T++){
      MFLOAT Ivar;
      MFLOAT deltaY = Yik - Yc(Y,P,T);

      MFLOAT LRik = FAIKP(deltaY,I-K-P,I,K,T,Ivar,Y);
      if(DEBUG>=2) assert(isfinite(LRik));

      int QminP = (P >= Pmin2) ? max(Qmin,Jmin[P]) : max(Qmin2, Jmin[P]);
      MFLOAT2 *ALpt = &A.AL(P,T,0);

      // Critical Loop for vectorization
      MFLOAT2 LRpq = 0.0;

      #pragma omp simd reduction(+:LRpq)
      for(int Q= min(J-1,Jmax[P]); Q >= QminP; Q--){
#if DEBUG>=2
	if(DEBUG>=2) assert(isfinite(ALpt[Q]) && ALpt[Q] >= 0.0);
#endif
        LRpq += ALpt[Q] * FAJ(Xm[J]-Xm[Q], deltaY, J-Q, I-K-P, FBiasWTm[J-Q], LogFpPow[J-Q-1], LogFnPow[I-K-P-1], PrBiasWTikp, LRik, Ivar,1, OUT_TYPE, OUT_OLD);
#if DEBUG>=2
	if(DEBUG>=2) assert(isfinite(LRp));
#endif
      }
      LRp += LRpq;
    }

    if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
    if(isfinite(ScaleLI[P]) || LRp != 0.0){// NEW144 : otherwise assume ScaleLI[P] * LRp == 0
      if(DEBUG>=2) assert(isfinite(ScaleLI[P]) && fabs(ScaleLI[P] - exp(LogScaleL[P] - LogScaleL[I])) < fabs(ScaleLI[P]) * 1e-5 + 1e-6);
      AL += LRp * ScaleLI[P];//  WAS exp(LogScaleL[min(max(P,IMIN-1),IMAX+1)] - LogScaleL[I]);
      if(DEBUG>=2)      assert(isfinite(AL) && AL >= 0.0);
    }
  }
  if(DEBUG && !(isfinite(AL) && AL >= 0.0)){
    if(DEBUG>=2)
      printf("\nm=%d,I=%d,K=%d,J=%d:AL=%0.10e (U=%d,Umax=%d,M=%d,N=%d,IMIN=%d,IMAX=%d)\n",m,I,K,J,(double)AL, origU, origUmax,M,N,IMIN,IMAX);
    else
      printf("\nm=%d,I=%d,K=%d,J=%d:AL=%0.10e (M=%d,N=%d,IMIN=%d,IMAX=%d)\n",m,I,K,J,(double)AL, M,N,IMIN,IMAX);
    fflush(stdout);
    assert(isfinite(AL) && AL >= 0.0);
  }
  return AL;
}

#define TRACE 0  // (m==17 && I==25 && K==1 && J==28 && origU==0 && origUmax == -2)
#define YTRACE 844.6893
#define HTRACE 10


/* specialized version of AL() that never updates ULmin,ULmax (Umax != -1) */
template<int OUT_TYPE, int OUT_OLD>
static inline double ALq(int m, /* for debugging only */
			 int I, 
			 int K,
			 int J,
			 int U, /* initialized to Umin */
			 int Ub,/* estimated peak of range U .. Umax (only used with ERFC_FIX if U <= Umax and can be outside the range U .. Umax) */
			 int Umax,/* If Umax < U, U range is unknown */
			 double Yik,
			 double Yt,/* New value of Y[t] */
			 int t, /* NOTE : t is either I or I-K */
			 double *  Xm,
			 double *  Y,// NOTE : ALq depends on Y[t] (see FAIKP), but should use new value Yt instead
			 int M,
			 int N,
			 int lc,
			 double Scale,/* exp(-LogScaleL[I]) */
			 double *  ScaleLI, /* ScaleLI[P] = exp(LogScaleL[min(max(P,IMIN-1),IMAX+1)] - LogScaleL[I]) */
			 double * LogScaleL,/* for debugging only */
			 QB &  A,
			 int IMIN,
			 int IMAX,
			 int *  Kmax,
			 int *  Jmin,
			 int *  Jmax)
{
#if VECTOR2_DEBUG
  int origU = U, origUmax = Umax, origUb = Ub, finalUmin = -1, finalUmax = -1, finalUb = Ub;
#endif

  double UL = PoutlierEnd;
  int JMAX = (EXTEND && I-K <= DELTA_YU) ? (DELTA_XU + EXTEND_MAX) : DELTA_XU;
  if(DEBUG>=2) assert(J > 0);

  if(J <= JMAX){
    double IvarX;
    double FLj = FLJ(J,Xm[J],IvarX);
    if(DEBUG>=2 && !(isfinite(FLj) && FLj >= 0.0 && FLj <= MaxUR * 3.0)){
      printf("AL(m=%d,I=%d,K=%d,J=%d,U=%d,%d):UL=%0.10e,IMIN=%d,IMAX=%d\n",
	     m,I,K,J,U,Umax,(double)UL,IMIN,IMAX);
      printf("FLJ(J,Xm[J],IvarX)=%0.10e:IvarX=%0.10e\n",FLj,IvarX);
      fflush(stdout);
      assert(isfinite(FLj) && FLj >= 0.0 && FLj <= MaxUR * 3.0);
    }

    double NotOutlier = FLIK_FIX ? PoutlierEndM1 : FLj * PoutlierEndM1;
    MFLOAT ULu;

    if(U <= Umax){
      if(DEBUG>=2) assert(U >= 0);
      if(DEBUG>=2) assert(Umax < I-K);

      ULu = 0.0f;
      if(lc && U <= 0){
	MFLOAT FLu = FLEIK(I,K,J,0,1,Yik,Xm[J],FLj);
	ULu += FLu;
#if VECTOR2_DEBUG
        if(TRACE && fabs(Yt - YTRACE) < 0.001){
	  printf("m=%d,I=%d,K=%d,J=%d,Yt=%0.4f,U=%d: FLEIK= %0.8e, ULu= %0.8e\n",m,I,K,J,Yt,U, FLu, ULu);
	  fflush(stdout);
	}
#endif
	if(DEBUG>=2)assert(isfinite(ULu) && ULu >= 0.0);
      }
      for(; U <= Umax; U++){
	MFLOAT FLu = (ERFC_FIX && U > Ub) ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
	  : FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);
	ULu += FLu;

#if VECTOR2_DEBUG
        if(TRACE && fabs(Yt - YTRACE) < 0.001){
	  printf("m=%d,I=%d,K=%d,J=%d,Yt=%0.4f,U=%d: FLu= %0.8e, ULu= %0.8e\n",m,I,K,J,Yt,U, FLu, ULu);
	  fflush(stdout);
	}
#endif
	if(DEBUG>=2)assert(isfinite(ULu) && ULu >= 0.0);
      }


    } else {/* need to determine U range */
      int Lijk = (IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]) ? max(1,A.Lij(I,K,J)) : 0;
      if(DEBUG>=2)assert(Lijk <= I-K);
      /* increment Lijk until Y[Lijk] overlap Xm */
      while(Lijk < I-K && Yik - Y[Lijk] > Xm[J])
	Lijk++;
      /* decrement Lijk until Y[Lijk-1] no longer overlaps Xm */
      while(Lijk > 0 && Yik - Y[Lijk-1] <= Xm[J])
	Lijk--;
      int Umin;
      Ub = Umin = Umax = max(1,Lijk) - 1;
      if(DEBUG>=2 && !(0 <= Umin && Umin < I-K)){
	printf("AL:I=%d,K=%d,J=%d:Umin=%d,Lijk=%d,IMIN=%d,IMAX=%d",I,K,J,Umin,Lijk,IMIN,IMAX);
	if(IMIN <= I && I <= IMAX){
	  printf(",Kmax[I]=%d,Jmin[I]=%d,Jmax[I]=%d",Kmax[I],Jmin[I],Jmax[I]);
	  if(0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I])
	    printf(",A[I][K][J].Lij=%d",A.Lij(I,K,J));
	}
	printf("\n");
	fflush(stdout);
	assert(0 <= Umin && Umin < I-K);
      }

      ULu = FLIK<0>(I,K,J,Umin,Umin+1,Yik,Xm[J],Yik-Y[Umin],Y[Umin+1]-Y[Umin],Y,IvarX, FLj, END_FIX ? lc : 1);
      if(DEBUG>=2)assert(isfinite(ULu) && ULu >= 0.0);
      MFLOAT bUL = ULu;
#if VECTOR2_DEBUG
      if(TRACE && fabs(Yt - YTRACE) < 0.001){
	float Yiu = Yik - Y[Umin], Yvu = Y[Umin+1] - Y[Umin];
	printf("m=%d,I=%d,K=%d,J=%d,Yt=%0.6f,Umin=%d: ULu= %0.8e (Xj=%0.6f,Yiu=%0.6f,Yvu=%0.6f,IvarX= %0.6e)\n",
	       m,I,K,J,Yt,Umin, ULu, Xm[J], Yiu, Yvu, IvarX);
	fflush(stdout);
      }
#endif
      if(DEBUG>=2)assert(isfinite(ULu) && ULu >= 0.0);
      if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N+1 && N+1 <= MDEBUG_N2 && m==MDEBUG_M && K==MDEBUG_K && J==MDEBUG_J){
	printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e,FLIK=%0.10e,FLj=%0.10e,PoutlierEndM1=%0.10e,PM=%0.10e):Yik=%0.4f,X[J]=%0.3f,Y[I-K]=%0.3f,Y[U]=%0.3f,Y[U+1]=%0.3f\n",
	       m,I,K,J,Umin,N,(double)(UL+ULu*NotOutlier),ULu*NotOutlier,ULu,FLj,PoutlierEndM1,PM(0,I,K,Y),Yik,Xm[J],Y[I-K],Y[Umin],Y[Umin+1]);
	fflush(stdout);
      }

      for(int U = Umax; ++U < I-K;){
	MFLOAT LR = ERFC_FIX ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
	  : FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);
	if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N+1 && N+1 <= MDEBUG_N2 && m==MDEBUG_M && K==MDEBUG_K && J==MDEBUG_J){
	  printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e,FLIK=%0.10e,FLj=%0.10e)\n",
		 m,I,K,J,U,N,(double)(UL+(ULu+LR)*NotOutlier),LR*NotOutlier, LR, FLj);
	  fflush(stdout);
	}
	if(DEBUG>=2)assert(isfinite(LR) && LR >= 0.0);

	if(LR < ULu * EPS)
	  break;
        if(URANGE_FIX && LR > bUL){
	  Umin = U;
	  ULu = 0.0;
	  bUL = LR;
        }
	Umax = U;
	ULu += LR;
	if(DEBUG>=2)assert(isfinite(ULu) && ULu >= 0.0);
      }

      for(int U = Umin ;--U >= 0;){
	MFLOAT LR = (ERFC_FIX && U > Ub) ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
	  : FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);
	if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N+1 && N+1 <= MDEBUG_N2 && m==MDEBUG_M && K==MDEBUG_K && J==MDEBUG_J){
	  printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e,FLIK=%0.10e,FLj=%0.10e)\n",
		 m,I,K,J,U,N,(double)(UL+(ULu+LR)*NotOutlier),LR*NotOutlier,LR,FLj);
	  fflush(stdout);
	}
	if(DEBUG>=2) assert(isfinite(LR) && LR >= 0.0);
	if(LR < ULu * EPS)
	  break;
	Umin = U;
	ULu += LR;
	if(DEBUG>=2) assert(isfinite(ULu) && ULu >= 0.0);
      }

      if(ISLINEAR && lc){
	MFLOAT LR = FLEIK(I,K,J,0,1,Yik,Xm[J],FLj);
	if(LR >= ULu * EPS){
	  Umin = 0;
	  ULu += LR;
	  if(DEBUG>=2) assert(isfinite(ULu) && ULu >= 0.0);
	}
      }
#if VECTOR2_DEBUG
      finalUmin = Umin;
      finalUmax = Umax;
      finalUb = Ub;
#endif
      if(DEBUG>=2)assert(Umax < I-K);
    }

    UL += ULu * NotOutlier;
    if(DEBUG>=2) assert(isfinite(UL) && UL >= 0.0);
  }
  if(DEBUG && !(isfinite(UL) && UL >= 0.0)){
    printf("\nAL(m=%d,I=%d,K=%d,J=%d,U=%d,Umax=%d):UL=%0.10e,IMIN=%d,IMAX=%d\n",
	   m,I,K,J,U,Umax,(double)UL,IMIN,IMAX);
    //    printf("\t origU=%d, origUmax= %d\n",origU,origUmax);
    fflush(stdout);
    assert(isfinite(UL) && UL >= 0.0);
  }
  double AL = UL * Scale;
#if VECTOR2_DEBUG
  if(TRACE && fabs(Yt - YTRACE) < 0.001){
    printf("ALq:m=%d,I=%d,K=%d,J=%d,U=%d..%d,Ub=%d (final U=%d..%d,Ub=%d), Yt= %0.4f: UL= %0.8e, Scale= %0.8e, UL * Scale= %0.8e\n",m,I,K,J,origU,origUmax,origUb,finalUmin,finalUmax,finalUb,Yt,UL,Scale,AL);
    fflush(stdout);
  }
#endif  

  if(DEBUG) assert(isfinite(AL) && AL >= 0.0);
  int Pmin = max(IMIN, I-K - deltaExtYRef);
  int Pmin2 = max(IMIN, I-K - DELTA_Y);
  int Qmin = max(1, J - deltaExtXRef);
  int Qmin2 = max(1, J - DELTA_X);

  if(DEBUG>=2 && !(I-K-1 <= IMAX)){
    printf("\n m=%d,I=%d,K=%d,J=%d:IMIN=%d,IMAX=%d,N=%d\n",m,I,K,J,IMIN,IMAX,N);
    fflush(stdout);
    assert(I-K-1 <= IMAX);
  }

  double Yi, Yimk;
  if(t == I){
    Yi = Yt;
    Yimk = K ? Y[I-K] : Yt;
  } else {
    if(DEBUG/* HERE >=2 */) assert(K > 0 && t == I-K /* NEW3 */ && I-K > 0);
    Yi = Y[I];
    Yimk = Yt;
  }

  double PrBiasWTik;/* value of PrBiasWT[I-K] may not be valid, since Y[t..I] are new labels, but differentials PrBias[0..t-1] remain valid */
  if(t==I && K > 0)
    PrBiasWTik = PrBiasWT[I-K];
  else { /* compute correct value of PrBiasWT[I-K] */
    double y = Yt - Y[t-1];
    double pr = Pr(y);
    if(DEBUG>=2) assert(isfinite(pr) && 0.0 <= pr && pr <= 1.0);
    double prbias = (pr <= 0.0 || pr >= 1.0) ? 0.0 : pr * log(pr) + (1.0 - pr) * log(1.0 - pr);    
    if(DEBUG>=2)  assert(isfinite(prbias) && prbias <= 0.0);
    PrBiasWTik = PrBiasWT[t-1] - PRbiasWT * prbias;
  }

  MFLOAT resR = Yi - Yimk;

  if(VERB>=3 && MDEBUG && m==MDEBUG_M && N==MDEBUG_N && I==2 && K==0 && J==19 && fabs(Yt - 763.509) < 0.01){
    printf("m=%d,I=%d,K=%d,J=%d,Pmin=%d,%d,Qmin=%d,%d:Yt= %0.4f: UL= %0.8e, Scale= %0.8e: AL= %0.8e\n",m,I,K,J,Pmin,Pmin2,Qmin,Qmin2,Yt,UL, Scale, AL);
    fflush(stdout);
  }

  for(int P = I-K; --P >= Pmin;){
    MFLOAT2 LRp = 0.0;
    MFLOAT PrBiasWTikp = PrBiasWTik - PrBiasWT[P];

    for(int T = 0; T <= Kmax[P]; T++){
      MFLOAT Ivar;
      MFLOAT deltaY = Yik - Yc(Y,P,T);
      MFLOAT resL = Y[P] - Y[P-T];
      MFLOAT Yint = Yimk - Y[P];
      
      MFLOAT LRik = FAIKP(deltaY, resR, resL, Yint, I-K-P,I,K,T,Ivar,Y);/* NOTE : does NOT index Y[] */
      if(DEBUG>=2) assert(isfinite(LRik));


      int QminP = (P >= Pmin2) ? max(Qmin,Jmin[P]) : max(Qmin2, Jmin[P]);
      MFLOAT2 *ALpt = &A.AL(P,T,0);

      // Critical Loop for vectorization
      MFLOAT2 LRpq = 0.0;

      #pragma omp simd reduction(+:LRpq)
      for(int Q= min(J-1,Jmax[P]); Q >= QminP; Q--){
	if(DEBUG>=2) assert(isfinite(ALpt[Q]) && ALpt[Q] >= 0.0);

	MFLOAT FAj = FAJ(Xm[J]-Xm[Q], deltaY, J-Q, I-K-P, FBiasWTm[J-Q], LogFpPow[J-Q-1], LogFnPow[I-K-P-1], PrBiasWTikp, LRik,Ivar,1,OUT_TYPE,OUT_OLD);
	LRpq += ALpt[Q] * FAj;

	if(VERB>=3 && MDEBUG && m==MDEBUG_M && N==MDEBUG_N && I==2 && K==0 && J==19 && fabs(Yt - 763.509) < 0.01){
	  printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:x=%0.4f,y=%0.4f:LRik= %0.8e,Ivar= %0.8e,FAj= %0.8e,ALpt[Q]= %0.8e, LRp= %0.8e, ScaleLI[P]= %0.8e, AL= %0.8e (delta= %0.8e)\n",
		 m,I,K,J,P,T,Q,Xm[J]-Xm[Q],deltaY,LRik,Ivar,FAj,ALpt[Q], LRpq+LRp, ScaleLI[P], AL + (LRpq + LRp) * ScaleLI[P], ALpt[Q] * FAj * ScaleLI[P]);
	  fflush(stdout);
	}

#if VECTOR2_DEBUG
	if(TRACE && fabs(Yt - YTRACE) < 0.001){
          printf("\t ALq:P=%d,T=%d,Q=%d:ALpt[Q]= %0.6e, x=%0.4f,y=%0.4f,LRik=%0.6e,Ivar=%0.6e:LRp= %0.8e,ScaleLI[P]= %0.8e,AL= %0.8e (delta= %0.8e)\n", 
		 P, T, Q, ALpt[Q], Xm[J]-Xm[Q], deltaY, LRik, Ivar, LRp + LRpq, ScaleLI[P], (LRp + LRpq) * ScaleLI[P] + AL, ALpt[Q] * FAj * ScaleLI[P]);
	  fflush(stdout);
        }
#endif

	if(DEBUG>=2) assert(isfinite(LRp));
      }
      LRp += LRpq;
    }

    if(DEBUG>=2) assert(IMIN <= P && P <= IMAX); // WAS assert(IMIN-1 <= P && P <= IMAX+1);
    if(isfinite(ScaleLI[P]) || LRp != 0.0){// NEW144 : otherwise assume ScaleLI[P] * LRp == 0
      if(DEBUG>=2) assert(isfinite(ScaleLI[P]) && fabs(ScaleLI[P] - exp(LogScaleL[P] - LogScaleL[I])) < fabs(ScaleLI[P]) * 1e-5 + 1e-6);
      AL += LRp * ScaleLI[P];//  WAS exp(LogScaleL[min(max(P,IMIN-1),IMAX+1)] - LogScaleL[I]);
      if(DEBUG>=2)      assert(isfinite(AL) && AL >= 0.0);
    }
  }
  if(DEBUG && !(isfinite(AL) && AL >= 0.0)){
    printf("\nm=%d,I=%d,K=%d,J=%d:AL=%0.10e (M=%d,N=%d,IMIN=%d,IMAX=%d)\n",m,I,K,J,(double)AL, M,N,IMIN,IMAX);
    fflush(stdout);
    assert(isfinite(AL) && AL >= 0.0);
  }
  return AL;
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT

/* AVX256 version of ALq() */
template<int OUT_TYPE, int OUT_OLD>
static inline void ALq_mm256(double *AL,/* return values AL[0..7] */
			     int m, /* for debugging only */
			     int I, 
			     int K,
			     int J,
			     int U, /* initialized to Umin */
                   	     int Ub,/* estimated peak of range U .. Umax (only used with ERFC_FIX if U <= Umax and can be outside the range U .. Umax) */
			     int Umax,/* If Umax < U, U range is unknown */
			     __m256d v_Yik1, __m256d v_Yik2, /* 8 New values for Yik */
			     __m256d v_Yt1, __m256d v_Yt2, /* 8 New values for Y[t] */
			     int t, /* NOTE : t is either I or I-K */
			     double *  Xm,
			     double *  Y,// NOTE : ALq depends on Y[t] (see FAIKP), but should use never use Y[t]
			     int M,
			     int N,
			     int lc,
			     double Scale,/* exp(-LogScaleL[I]) */
			     double *  ScaleLI, /* ScaleLI[P] = exp(LogScaleL[min(max(P,IMIN-1),IMAX+1)] - LogScaleL[I]) */
			     double * LogScaleL,/* for debugging only */
			     QB &  A,
			     int IMIN,
			     int IMAX,
			     int *  Kmax,
			     int *  Jmin,
			     int *  Jmax)
{
#if VECTOR2_DEBUG
  int origU = U, origUmax = Umax, origUb = Ub, finalUmin = -1, finalUmax = -1, finalUb = -1;
  double ALqh[8];
  double Yik[8];
  double Yt[8];
  _mm256_storeu_pd(Yik, v_Yik1);  _mm256_storeu_pd(&Yik[4], v_Yik2);
  _mm256_storeu_pd(Yt, v_Yt1); _mm256_storeu_pd(&Yt[4], v_Yt2);
  for(int h = 0; h < 8; h++){
    if(TRACE){
      printf("Calling ALq:m=%d,I=%d,K=%d,J=%d,U=%d..%d,Ub=%d,h=%d,Yik=%0.4f,Yt=%0.4f,t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d\n",m,I,K,J,U,Umax,Ub,h,Yik[h],Yt[h],t,M,N,IMIN,IMAX);
      fflush(stdout);
    }
    ALqh[h] = ALq<OUT_TYPE,OUT_OLD>(m,I,K,J,U,Ub,Umax,Yik[h],Yt[h],t,Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax);
  }
#endif  

  __m256d v_UL1 = _mm256_set1_pd(PoutlierEnd);// UL = PoutlierEnd
  __m256d v_UL2 = v_UL1;

  int JMAX = (EXTEND && I-K <= DELTA_YU) ? (DELTA_XU + EXTEND_MAX) : DELTA_XU;

  if(J <= JMAX){
    double IvarX;
    double FLj = FLJ(J,Xm[J],IvarX);
    if(DEBUG>=2) assert(isfinite(FLj) && FLj >= 0.0 && FLj <= MaxUR * 3.0);
    double NotOutlier = FLIK_FIX ? PoutlierEndM1 : FLj * PoutlierEndM1;

    // NOTE : U range is U .. Umax, with Ub corresponding to Lij-1, the estimated peak value (can be outside the range U..Umax)

    if(U > Umax){/* need to determine U range : Use first and last Yik values to compute upper bound of range */

      double Yik8[8] __attribute__((aligned(32)));
      _mm256_store_pd(Yik8, v_Yik1);
      _mm256_store_pd(&Yik8[4], v_Yik2);
      int Ulb = I, Uub = 0, Umb = -1;/* U range will be computed Ulb .. Uub, with Umb corresponding to Lij-1, the estimate peak value (can lie outside the range Ulb .. Uub */

      for(int T = 0; T <= 1; T++){
	double Yik = T ? Yik8[7] : Yik8[0];

	int Lijk = (IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]) ? max(1,A.Lij(I,K,J)) : 0;

	/* increment Lijk until Y[Lijk] overlap Xm */
	while(Lijk < I-K && Yik - Y[Lijk] > Xm[J])
	  Lijk++;
	/* decrement Lijk until Y[Lijk-1] no longer overlaps Xm */
	while(Lijk > 0 && Yik - Y[Lijk-1] <= Xm[J])
	  Lijk--;

	int Umid, Umin;
	Umid = Umin = Umax = max(1,Lijk) - 1;
	MFLOAT ULu = FLIK<0>(I,K,J,Umin,Umin+1,Yik,Xm[J],Yik-Y[Umin],Y[Umin+1]-Y[Umin],Y,IvarX, FLj, END_FIX ? lc : 1);
	MFLOAT bUL = ULu;

	for(int U = Umax; ++U < I-K;){
	  MFLOAT LR = ERFC_FIX ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
	    : FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);
	  if(LR < ULu * EPS)
	    break;
	  if(URANGE_FIX && LR > bUL){
	    Umin = U;
	    ULu = 0.0;
	    bUL = LR;
	  }
	  Umax = U;
	  ULu += LR;
	}

	for(int U = Umin ;--U >= 0;){
	  MFLOAT LR = (ERFC_FIX && U > Umid) ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
	    : FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);
	  if(LR < ULu * EPS)
	    break;
	  Umin = U;
	  ULu += LR;
	}

	if(ISLINEAR && lc){
	  MFLOAT LR = FLEIK(I,K,J,0,1,Yik,Xm[J],FLj);
	  if(LR >= ULu * EPS){
	    Umin = 0;
	    ULu += LR;
	  }
	}

	if(DEBUG>=2 && !(0 <= Umin && Umin <= Umax && Umax < I-K)){
	  printf("ALq_mm256:m=%d,I=%d,K=%d,J=%d:T=%d:Umin=%d,Umax=%d\n",m,I,K,J,T,Umin,Umax);
	  fflush(stdout);
	  assert(0 <= Umin && Umin <= Umax && Umax < I-K);
	}

	Ulb = min(Ulb, Umin);
	Uub = max(Uub, Umax);
	Umb = T ? (Umb + Umid)/2 : Umid;
#if VECTOR2_DEBUG
	if(TRACE){
	  printf("ALq_mm256:T=%d,Yik= %0.6f,Umin=%d,Umax=%d,Umid=%d: Ulb=%d,Uub=%d,Umb=%d\n",	 T,Yik,Umin,Umax,Umid,Ulb,Uub,Umb);
	  fflush(stdout);
	}
#endif
      } // T = 0..1
      
      if(DEBUG) assert(Ulb <= Uub);

      U = Ulb;
      Ub = Umb;
      Umax = Uub;

#if VECTOR2_DEBUG
      finalUmin = U;
      finalUmax = Umax;
      finalUb = Ub;
#endif
    } // U > Umax

    if(DEBUG>=2 && !(0 <= U && U <= Umax && Umax < I-K)){
      printf("ALq_mm256:m=%d,I=%d,K=%d,J=%d:U=%d,Umax=%d\n",m,I,K,J,U,Umax);
      fflush(stdout);
      assert(0 <= U && U <= Umax && Umax < I-K);
    }

    __m256 v_Xj = _mm256_set1_ps((float)Xm[J]);
    __m256 v_IvarX = _mm256_set1_ps((float)IvarX);
    __m256 v_FLj = _mm256_set1_ps((float)FLj);
							   
    __m256 ULu = _mm256_setzero_ps();

    if(lc && U <= 0) {
      __m128 v_YikF1 = _mm256_cvtpd_ps(v_Yik1);
      __m128 v_YikF2 = _mm256_cvtpd_ps(v_Yik2);
      ULu += FLEIK_mm256(I,K,0,1, _mm256_set_m128(v_YikF2, v_YikF1), v_Xj, v_FLj);
    }
    for(; U <= Umax; U++){
      __m256d Yu = _mm256_set1_pd(Y[U]);
      __m128 v_Yiu1 = _mm256_cvtpd_ps(v_Yik1 - Yu);
      __m128 v_Yiu2 = _mm256_cvtpd_ps(v_Yik2 - Yu);
      __m256 v_Yiu = _mm256_set_m128(v_Yiu2, v_Yiu1);
      __m256 v_Yvu = _mm256_set1_ps((float)(Y[U+1]-Y[U]));      

      ULu += (ERFC_FIX && U > Ub) ? FLIK_mm256<1>(I,K,U,U+1, v_Xj, v_Yiu, v_Yvu, v_IvarX, v_FLj, END_FIX ? lc : 1)
	: FLIK_mm256<0>(I,K,U,U+1, v_Xj, v_Yiu, v_Yvu, v_IvarX, v_FLj, END_FIX ? lc : 1);
    }

    /* convert ULu to double precision */
    __m256d v_ULu1 = _mm256_cvtps_pd(_mm256_extractf128_ps(ULu,0));/* low order 4 floats to double */
    __m256d v_ULu2 = _mm256_cvtps_pd(_mm256_extractf128_ps(ULu,1));/* high order 4 floats to double */
    __m256d v_NotOutlier = _mm256_set1_pd(NotOutlier);
    v_UL1 += v_ULu1 * v_NotOutlier;
    v_UL2 += v_ULu2 * v_NotOutlier;
  }

#if DEBUG>=2
    double UL8[8];
    _mm256_store_pd(UL8, v_UL1);
    _mm256_store_pd(&UL8[4], v_UL2);
    for(int h = 0; h < 8; h++)
      assert(isfinite(UL8[h]) && UL8[h] >= 0.0);
#endif

  __m256d v_Scale = _mm256_set1_pd(Scale);
  __m256d v_AL1 = v_UL1 * v_Scale;
  __m256d v_AL2 = v_UL2 * v_Scale;

#if DEBUG>=2
  double AL8[8];
  _mm256_store_pd(AL8, v_AL1);
  _mm256_store_pd(&AL8[4], v_AL2);
  for(int h = 0; h < 8; h++)
    assert(isfinite(AL8[h]) && AL8[h] >= 0.0);
#if VECTOR2_DEBUG
  if(TRACE){
    printf("ALq_mm256:m=%d,I=%d,K=%d,J=%d,U=%d..%d,Ub=%d (final U=%d..%d,Ub=%d),h=%d: UL[h]= %0.8e, Scale= %0.8e, UL[h]*Scale = %0.8e\n",
	   m,I,K,J,origU,origUmax,origUb,finalUmin,finalUmax,finalUb, HTRACE,UL8[HTRACE],Scale,AL8[HTRACE]);
    fflush(stdout);
  }
#endif    
#endif

  int Pmin = max(IMIN, I-K - deltaExtYRef);
  int Pmin2 = max(IMIN, I-K - DELTA_Y);
  int Qmin = max(1, J - deltaExtXRef);
  int Qmin2 = max(1, J - DELTA_X);

  __m256d v_Yi1, v_Yi2, v_Yimk1, v_Yimk2;
  if(t == I){
    v_Yi1 = v_Yt1; v_Yi2 = v_Yt2;
    if(K){
      v_Yimk1 = _mm256_set1_pd(Y[I-K]); v_Yimk2 = v_Yimk1;
    } else {
      v_Yimk1 = v_Yt1; v_Yimk2 = v_Yt2;
    }
  } else {
    if(DEBUG/* HERE >=2 */) assert(K > 0 && t == I-K /* NEW3 */&& I-K > 0);
    v_Yi1 = _mm256_set1_pd(Y[I]); v_Yi2 = v_Yi1;// Yi = Y[I]
    v_Yimk1 = v_Yt1; v_Yimk2 = v_Yt2; // Yimk = Yt
  }

  __m256d v_PrBiasWTik1, v_PrBiasWTik2;/* values of PrBiasWT[I-K] may not be valid, since Y[t..I] are new labels, but differentials PrBias[0..t-1] remain valid */
  if(t==I && K > 0)
    v_PrBiasWTik1 = v_PrBiasWTik2 = _mm256_set1_pd(PrBiasWT[I-K]);
  else { /* compute correct value of PrBiasWT[I-K] as a function of Yt */
    __m256d ytm1 = _mm256_set1_pd(Y[t-1]);
    __m256d v_y1 = v_Yt1 - ytm1, v_y2 = v_Yt2 - ytm1;// y = Yt - Y[t-1]
      
    __m256d v_pr1 = Pr_mm256d(v_y1), v_pr2 = Pr_mm256d(v_y2);// pr = Pr(y)
    __m256d v_1Mpr1 = v256d_1p0 - v_pr1, v_1Mpr2 = v256d_1p0 - v_pr2;// 1.0 - pr
    // prbias = pr * log(max(pr,1e-300)) + (1.0 - pr) * log(max(1.0 - pr, 1e-300))
    __m256d v_eps = _mm256_set1_pd(1e-300);
    __m256d v_prbias1 = v_pr1 * log_mm256d(_mm256_max_pd(v_pr1, v_eps)) + v_1Mpr1 * log_mm256d(_mm256_max_pd(v_1Mpr1, v_eps));
    __m256d v_prbias2 = v_pr2 * log_mm256d(_mm256_max_pd(v_pr2, v_eps)) + v_1Mpr2 * log_mm256d(_mm256_max_pd(v_1Mpr2, v_eps));
    __m256d v_PrBiasWTtm1 = _mm256_set1_pd(PrBiasWT[t-1]);
    __m256d v_PRbiasWT = _mm256_set1_pd(PRbiasWT);
    v_PrBiasWTik1 = v_PrBiasWTtm1 - v_PRbiasWT * v_prbias1; v_PrBiasWTik2 = v_PrBiasWTtm1 - v_PRbiasWT * v_prbias2;// PRbiasWTik = PrBiasWT[t-1] - PRbiasWT * prbias
  }

  __m256 v_resR = _mm256_set_m128(_mm256_cvtpd_ps(v_Yi2 - v_Yimk2), _mm256_cvtpd_ps(v_Yi1 - v_Yimk1));// resR = Yi - Yimk

  for(int P = I-K; --P >= Pmin;){
    __m256 v_LRp = _mm256_setzero_ps(); // LRp = 0.0

    for(int T = 0; T <= Kmax[P]; T++){
      __m256 v_Ivar;
      __m256d Ypt = _mm256_set1_pd(Yc(Y,P,T));
      __m256 v_deltaY = _mm256_set_m128(_mm256_cvtpd_ps(v_Yik2 - Ypt), _mm256_cvtpd_ps(v_Yik1 - Ypt));// deltaY = Yik - Yc(Y,P,T)
      __m256 v_resL = _mm256_set1_ps((float)(Y[P] - Y[P-T])); // resL = Y[P] - Y[P-T]
      __m256d Yp = _mm256_set1_pd(Y[P]);
      __m256 v_Yint =  _mm256_set_m128(_mm256_cvtpd_ps(v_Yimk2 - Yp), _mm256_cvtpd_ps(v_Yimk1 - Yp));// Yint = Yimk - Y[P]
      
      __m256 v_LRik = FAIKP_mm256(v_deltaY, v_resR, v_resL, v_Yint, v_Ivar);

      int QminP = (P >= Pmin2) ? max(Qmin,Jmin[P]) : max(Qmin2, Jmin[P]);
      MFLOAT2 *ALpt = &A.AL(P,T,0);

      __m256 v_LogFnPowIKP = _mm256_broadcast_ss(&LogFnPow[I-K-P-1]);

      __m256d PrBiasWTp = _mm256_set1_pd(PrBiasWT[P]);
      __m256 v_PrBiasWTikp = _mm256_set_m128(_mm256_cvtpd_ps(v_PrBiasWTik2 - PrBiasWTp), _mm256_cvtpd_ps(v_PrBiasWTik1 - PrBiasWTp));// PrBiasWTip = PrBiasWTi - PrBiasWT[P]

      // Critical Loop
      for(int Q= min(J-1,Jmax[P]); Q >= QminP; Q--){
	v_LRp += _mm256_set1_ps(ALpt[Q]) * FAJ_mm256<1,OUT_TYPE,OUT_OLD>((float)(Xm[J]-Xm[Q]), v_deltaY, J-Q, I-K-P, FBiasWTm[J-Q], LogFpPow[J-Q-1], v_LogFnPowIKP, v_PrBiasWTikp, v_LRik, v_Ivar);
#if VECTOR2_DEBUG
	if(TRACE){
	  float LRp[8], deltaY[8], LRik[8], Ivar[8];
	  _mm256_store_ps(LRp, v_LRp);
	  _mm256_store_pd(AL, v_AL1);
	  _mm256_store_pd(&AL[4], v_AL2);
	  _mm256_store_ps(deltaY, v_deltaY);
	  _mm256_store_ps(LRik, v_LRik);
	  _mm256_store_ps(Ivar, v_Ivar);

	  printf("\t ALq_mm256:P=%d,T=%d,Q=%d,h=%d:ALpt[Q]= %0.6e,x=%0.4f,y=%0.4f,LRik=%0.6e,Ivar=%0.6e,:LRp= %0.6e, ScaleLI[P]= %0.6e (AL= %0.6e)\n", 
		 P, T, Q, HTRACE, ALpt[Q], Xm[J]-Xm[Q], deltaY[HTRACE], LRik[HTRACE], Ivar[HTRACE], LRp[HTRACE], ScaleLI[P], LRp[HTRACE] * ScaleLI[P] + AL[HTRACE]);
	  fflush(stdout);
	}
#endif	
      }
    }

    /* convert v_LRp to double precision */
    __m256d v_LRp1 = _mm256_cvtps_pd(_mm256_extractf128_ps(v_LRp, 0));/* low order 4 floats to double */
    __m256d v_LRp2 = _mm256_cvtps_pd(_mm256_extractf128_ps(v_LRp, 1));/* high order 4 floats to double */
    if(DEBUG>=2) assert(IMIN <= P && P <= IMAX); // WAS assert(IMIN-1 <= P && P <= IMAX+1);
    if(!isfinite(ScaleLI[P])){/* NEW144 : check if all LRp[] values are 0 */
      double LRp[8];
      _mm256_store_pd(LRp, v_LRp1);
      _mm256_store_pd(&LRp[4], v_LRp2);
      double maxLRp = 0.0;
      for(int h = 0; h < 8; h++)
	maxLRp = max(maxLRp, LRp[h]);
      if(maxLRp > 0.0){
	printf("ALq_mm256:m=%d,I=%d,K=%d,J=%d,P=%d:t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d:\n",m,I,K,J,P,t,M,N,IMIN,IMAX);
	printf("\t ScaleLI[P]= %0.17e, LogScaleL[P]= %0.17f, LogScaleL[I]= %0.17f, maxLRp= %0.17e\n",ScaleLI[P], LogScaleL[P], LogScaleL[I],maxLRp);
	fflush(stdout);
	assert(isfinite(ScaleLI[P]));
      }
    } else {
      if(DEBUG>=2) assert(fabs(ScaleLI[P] - exp(LogScaleL[P] - LogScaleL[I])) < fabs(ScaleLI[P]) * 1e-5 + 1e-6);
      __m256d ScaleLIP = _mm256_set1_pd(ScaleLI[P]);
      v_AL1 += v_LRp1 * ScaleLIP; v_AL2 += v_LRp2 * ScaleLIP;// AL += LRp * ScaleLI[P]
    }
  }

  /* return results */
  _mm256_store_pd(AL, v_AL1);
  _mm256_store_pd(&AL[4], v_AL2);

  if(DEBUG/* HERE >=2 */){
    for(int h = 0; h < 8; h++){
#if VECTOR2_DEBUG      
      if(DEBUG && !(origUmax < 0 ? AL[h] - ALqh[h] > -fabs(ALqh[h]) * 1e-4 - 1e-9 
		    : fabs(ALqh[h] - AL[h]) < fabs(ALqh[h]) * 1e-4 + 1e-9)){
	#pragma omp critical
	{
	  printf("ALq_mm256:m=%d,I=%d,K=%d,J=%d,U=%d..%d,Ub=%d,Yik[h]=%0.4f,Yt[h]=%0.4f,t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d:\n",m,I,K,J,origU,origUmax,origUb,Yik[h],Yt[h],t,M,N,IMIN,IMAX);
	  printf("\t h=%d: AL[h]= %0.8e, ALqh[h]= %0.8e (final U=%d..%d,Ub=%d)\n",h,AL[h],ALqh[h],finalUmin,finalUmax,finalUb);
	  fflush(stdout);
	  if(origUmax < 0)
	    assert(AL[h] - ALqh[h] > -(fabs(ALqh[h]) * 1e-4 - 1e-9));
	  else
	    assert(fabs(ALqh[h] - AL[h]) < fabs(ALqh[h]) * 1e-4 + 1e-9);
	}
      }
#endif
      assert(isfinite(AL[h]) && AL[h] >= 0.0);
    }
  }
}

#elif VECTOR_MIC && ( USE_MIC || USE_AVX512 ) && USE_MFLOAT && DEBUG < 2

/* MIC AVX512 version of ALq() */
template<int OUT_TYPE, int OUT_OLD>
static inline void ALq_mm512(double *AL,/* return values AL[0..7] */
			     int m, /* for debugging only */
			     int I, 
			     int K,
			     int J,
			     int U, /**< initialized to Umin */
                   	     int Ub,/* estimated peak of range U .. Umax (only used with ERFC_FIX if U <= Umax and can be outside the range U .. Umax) */
			     int Umax,/**< If Umax < U, U range is unknown */
			     __m512d v_Yik1, __m512d v_Yik2, /* 16 New values for Yik */
			     __m512d v_Yt1, __m512d v_Yt2, /* 16 New values for Y[t] */
			     int t, /* NOTE : t is either I or I-K */
			     double *Xm,
			     double *Y,// NOTE : ALq depends on Y[t] (see FAIKP), but should use never use Y[t]
			     int M,
			     int N,
			     int lc,
			     double Scale,/* exp(-LogScaleL[I]) */
			     double *ScaleLI, /* ScaleLI[P] = exp(LogScaleL[min(max(P,IMIN-1),IMAX+1)] - LogScaleL[I]) */
			     double * LogScaleL,/* for debugging only */
			     QB &A,
			     int IMIN,
			     int IMAX,
			     int *Kmax,
			     int *Jmin,
			     int *Jmax)
{
#if VECTOR2_DEBUG
  int origU = U, origUmax = Umax, origUb = Ub, finalUmin = -1, finalUmax = -1, finalUb = -1;
  double ALqh[16];
  double Yik[16];
  double Yt[16];
  _mm512_storeu_pd(Yik, v_Yik1);  _mm512_storeu_pd(&Yik[8], v_Yik2);
  _mm512_storeu_pd(Yt, v_Yt1); _mm512_storeu_pd(&Yt[8], v_Yt2);
  for(int h = 0; h < 16; h++){
    if(TRACE){
      printf("Calling ALq:m=%d,I=%d,K=%d,J=%d,U=%d..%d,Ub=%d,h=%d,Yik=%0.4f,Yt=%0.4f,t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d\n",m,I,K,J,U,Umax,Ub,h,Yik[h],Yt[h],t,M,N,IMIN,IMAX);
      fflush(stdout);
    }
    ALqh[h] = ALq<OUT_TYPE,OUT_OLD>(m,I,K,J,U,Ub,Umax,Yik[h],Yt[h],t,Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax);
  }
#endif  

  __m512d v_UL1 = _mm512_set1_pd(PoutlierEnd);// UL = PoutlierEnd
  __m512d v_UL2 = v_UL1;

  if(DEBUG>=2){
    double UL16[16];
    _mm512_storeu_pd(UL16, v_UL1);
    _mm512_storeu_pd(&UL16[8], v_UL2);
    for(int t = 0; t < 16; t++)
      assert(isfinite(UL16[t]) && UL16[t] >= 0.0);
  }

  int JMAX = (EXTEND && I-K <= DELTA_YU) ? (DELTA_XU + EXTEND_MAX) : DELTA_XU;

  if(J <= JMAX){
    double IvarX;
    double FLj = FLJ(J,Xm[J],IvarX);
    if(DEBUG>=2) assert(isfinite(FLj) && FLj >= 0.0 && FLj <= MaxUR * 3.0);
    double NotOutlier = FLIK_FIX ? PoutlierEndM1 : FLj * PoutlierEndM1;

    // NOTE : U range is U .. Umax, with Ub corresponding to Lij-1, the estimated peak value (can be outside the range U..Umax)

    if(U > Umax){/* need to determine U range : Use first and last Yik values to compute upper bound of range */
      double YikV[16] __attribute__((aligned(64)));
      _mm512_store_pd(YikV, v_Yik1);
      _mm512_store_pd(&YikV[8], v_Yik2);
      int Ulb = I, Uub = 0, Umb = -1;/* U range will be computed Ulb .. Umb .. Uub, with Umb corresponding to Lij-1, the estimate peak value */

      for(int T = 0; T <= 1; T++){
	double Yik = T ? YikV[15] : YikV[0];

	int Lijk = (IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]) ? max(1,A.Lij(I,K,J)) : 0;

	/* increment Lijk until Y[Lijk] overlap Xm */
	while(Lijk < I-K && Yik - Y[Lijk] > Xm[J])
	  Lijk++;
	/* decrement Lijk until Y[Lijk-1] no longer overlaps Xm */
	while(Lijk > 0 && Yik - Y[Lijk-1] <= Xm[J])
	  Lijk--;

	int Umid, Umin;
	Umid = Umin = Umax = max(1,Lijk) - 1;
	MFLOAT ULu = FLIK<0>(I,K,J,Umin,Umin+1,Yik,Xm[J],Yik-Y[Umin],Y[Umin+1]-Y[Umin],Y,IvarX, FLj,END_FIX ? lc : 1);
	MFLOAT bUL = ULu;
	
	for(int U = Umax; ++U < I-K;){
	  MFLOAT LR = ERFC_FIX ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
	    : FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);
	  if(LR < ULu * EPS)
	    break;
	  if(URANGE_FIX && LR > bUL){
	    Umin = U;
	    ULu = 0.0;
	    bUL = LR;
	  }
	  Umax = U;
	  ULu += LR;
	}

	for(int U = Umin ;--U >= 0;){
	  MFLOAT LR = (ERFC_FIX && U > Umid) ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
	    : FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);

	  if(LR < ULu * EPS)
	    break;
	  Umin = U;
	  ULu += LR;
	}

	if(ISLINEAR && lc){
	  MFLOAT LR = FLEIK(I,K,J,0,1,Yik,Xm[J],FLj);
	  if(LR >= ULu * EPS){
	    Umin = 0;
	    ULu += LR;
	  }
	}

	if(DEBUG>=2 && !(0 <= Umin && Umin <= Umax && Umax < I-K)){
	  printf("ALq_mm512:m=%d,I=%d,K=%d,J=%d:T=%d:Umin=%d,Umax=%d\n",m,I,K,J,T,Umin,Umax);
	  fflush(stdout);
	  assert(0 <= Umin && Umin <= Umax && Umax < I-K);
	}

	Ulb = min(Ulb, Umin);
	Uub = max(Uub, Umax);
	Umb = T ? (Umb + Umid)/2 : Umid;
#if VECTOR2_DEBUG
	if(TRACE){
	  printf("ALq_mm512:T=%d,Yik= %0.6f,Umin=%d,Umax=%d,Umid=%d: Ulb=%d,Uub=%d,Umb=%d\n", T,Yik,Umin,Umax,Umid,Ulb,Uub,Umb);
	  fflush(stdout);
	}
#endif
      } // T = 0..1
      
      if(DEBUG) assert(Ulb <= Uub);
      U = Ulb;
      Ub = Umb;
      Umax = Uub;
#if VECTOR2_DEBUG
      finalUmin = U;
      finalUmax = Umax;
      finalUb = Ub;
#endif
    } // U > Umax

    if(DEBUG>=2 && !(0 <= U && U <= Umax && Umax < I-K)){
      printf("ALq_mm512:m=%d,I=%d,K=%d,J=%d:U=%d,Umax=%d\n",m,I,K,J,U,Umax);
      fflush(stdout);
      assert(0 <= U && U <= Umax && Umax < I-K);
    }

    __m512 v_Xj = _mm512_set1_ps((float)Xm[J]);
    __m512 v_IvarX = _mm512_set1_ps((float)IvarX);
    __m512 v_FLj = _mm512_set1_ps((float)FLj);
							   
    __m512 v_ULu = _mm512_setzero_ps();

#if 0 // HERE : remove this when Intel Compiler bug is fixed
    if(lc && U <= 0) 
      v_ULu = _mm512_add_ps(v_ULu, 
			    FLEIK_mm512(I,K,0,1, DoubleToFloat_mm512(v_Yik1, v_Yik2), v_Xj, v_FLj));

#endif // HERE

    for(; U <= Umax; U++){
      // Yiu = (float)(Y[U] - Yik)
      __m512d Yu = _mm512_set1_pd(Y[U]);
      __m512 v_Yiu = DoubleToFloat_mm512(_mm512_sub_pd(v_Yik1, Yu),_mm512_sub_pd(v_Yik2, Yu));

      __m512 v_Yvu = _mm512_set1_ps((float)(Y[U+1]-Y[U]));// Yvu = Y[U+1]-Y[U]

      v_ULu = _mm512_add_ps(v_ULu, 
			    (ERFC_FIX && U > Ub) ? FLIK_mm512<1>(I,K,U,U+1, v_Xj, v_Yiu, v_Yvu, v_IvarX, v_FLj, END_FIX ? lc : 1)
			    : FLIK_mm512<0>(I,K,U,U+1, v_Xj, v_Yiu, v_Yvu, v_IvarX, v_FLj, END_FIX ? lc : 1));
    }

    /* convert v_ULu to double precision */
#if USE_MIC
    __m512d v_ULu1 = _mm512_cvtpslo_pd(v_ULu);/* low order 8 floats to double */
    __m512d v_ULu2 = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_ULu, 0x2 | (0x3 << 2)));/* high order 8 floats to double */
#elif USE_AVX512
    __m512d v_ULu1 = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_castsi512_ps(_mm512_castsi256_si512(_mm512_extracti64x4_epi64(_mm512_castps_si512(v_ULu),0)))));/* low order 8 floats to double */
    __m512d v_ULu2 = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_castsi512_ps(_mm512_castsi256_si512(_mm512_extracti64x4_epi64(_mm512_castps_si512(v_ULu),1)))));/* high order 8 floats to double */
#endif

    // UL += ULu * FLj
    __m512d v_NotOutlier = _mm512_set1_pd(NotOutlier);
    v_UL1 = _mm512_add_pd(v_UL1, _mm512_mul_pd(v_ULu1 , v_NotOutlier));
    v_UL2 = _mm512_add_pd(v_UL2, _mm512_mul_pd(v_ULu2 , v_NotOutlier));
  }

#if DEBUG>=2 || VECTOR2_DEBUG
  double UL16[16];
  _mm512_storeu_pd(UL16, v_UL1);
  _mm512_storeu_pd(&UL16[8], v_UL2);
  for(int t = 0; t < 16; t++)
    assert(isfinite(UL16[t]) && UL16[t] >= 0.0);
#endif

  __m512d v_Scale = _mm512_set1_pd(Scale);
  __m512d v_AL1 = _mm512_mul_pd(v_UL1, v_Scale);
  __m512d v_AL2 = _mm512_mul_pd(v_UL2, v_Scale);

#if DEBUG>=2 || VECTOR2_DEBUG
  double AL16[16];
  _mm512_storeu_pd(AL16, v_AL1);
  _mm512_storeu_pd(&AL16[8], v_AL2);
  for(int t = 0; t < 16; t++)
    assert(isfinite(AL16[t]) && AL16[t] >= 0.0);
  if(TRACE){
    printf("m=%d,I=%d,K=%d,J=%d,U=%d..%d,Ub=%d (final U=%d..%d,Ub=%d),h=%d: UL[h[= %0.8e, Scale= %0.8e, UL[h]*Scale = %0.8e\n",
	   m,I,K,J,origU,origUmax,origUb,finalUmin,finalUmax,finalUb,HTRACE, UL16[HTRACE], Scale, AL16[HTRACE]);
    fflush(stdout);
  }
#endif    

  int Pmin = max(IMIN, I-K - deltaExtYRef);
  int Pmin2 = max(IMIN, I-K - DELTA_Y);
  int Qmin = max(1, J - deltaExtXRef);
  int Qmin2 = max(1, J - DELTA_X);


  __m512d v_Yi1, v_Yi2, v_Yimk1, v_Yimk2;
  if(t == I){
    v_Yi1 = v_Yt1; v_Yi2 = v_Yt2;
    if(K){
      v_Yimk1 = _mm512_set1_pd(Y[I-K]); v_Yimk2 = v_Yimk1;
    } else {
      v_Yimk1 = v_Yt1; v_Yimk2 = v_Yt2;
    }
  } else {
    if(DEBUG/* HERE >=2 */) assert(K > 0 && t == I-K);
    v_Yi1 = _mm512_set1_pd(Y[I]); v_Yi2 = v_Yi1;// Yi = Y[I]
    v_Yimk1 = v_Yt1; v_Yimk2 = v_Yt2; // Yimk = Yt
  }

  __m512d v_PrBiasWTik1, v_PrBiasWTik2;/* values of PrBiasWT[I-K] is not valid, since Y[t..I] are new labels, but differentials PrBias[0..t-1] remain valid */
  if(t==I && K > 0)
    v_PrBiasWTik1 = v_PrBiasWTik2 = _mm512_set1_pd(PrBiasWT[I-K]);
  else {/* compute correct value of PrBiasWT[I-K] as a function of Yt */
    __m512d v_1p0 = _mm512_set1_pd(1.0);
    __m512d v_eps = _mm512_set1_pd(1e-300);

    __m512d ytm1 = _mm512_set1_pd(Y[t-1]);
    __m512d v_y1 = _mm512_sub_pd(v_Yt1, ytm1), v_y2 = _mm512_sub_pd(v_Yt2, ytm1);// y = Yt - Y[t-1]

    __m512d v_pr1 = Pr_mm512d(v_y1), v_pr2 = Pr_mm512d(v_y2);// pr = Pr(y)
    __m512d v_1Mpr1 = _mm512_sub_pd(v_1p0, v_pr1), v_1Mpr2 = _mm512_sub_pd(v_1p0, v_pr2);// 1.0 - pr
    // prbias = pr * log(max(pr,1e-300)) + (1.0 - pr) * log(max(1.0 - pr, 1e-300))
    __m512d v_prbias1 = _mm512_add_pd(_mm512_mul_pd(v_pr1, log_mm512d(_mm512_max_pd(v_pr1, v_eps))), _mm512_mul_pd(v_1Mpr1, log_mm512d(_mm512_max_pd(v_1Mpr1, v_eps))));
    __m512d v_prbias2 = _mm512_add_pd(_mm512_mul_pd(v_pr2, log_mm512d(_mm512_max_pd(v_pr2, v_eps))), _mm512_mul_pd(v_1Mpr2, log_mm512d(_mm512_max_pd(v_1Mpr2, v_eps))));

    __m512d v_PrBiasWTtm1 = _mm512_set1_pd(PrBiasWT[t-1]);
    __m512d v_PRbiasWT = _mm512_set1_pd(PRbiasWT);
      v_PrBiasWTik1 = _mm512_sub_pd(v_PrBiasWTtm1, _mm512_mul_pd(v_PRbiasWT, v_prbias1));  v_PrBiasWTik2 = _mm512_sub_pd(v_PrBiasWTtm1, _mm512_mul_pd(v_PRbiasWT, v_prbias2));// PrBiasWTik = PrBiasWT[t-1] - PrBiasWT * prbias
  }

  __m512 v_resR = DoubleToFloat_mm512(_mm512_sub_pd(v_Yi1, v_Yimk1), _mm512_sub_pd(v_Yi2, v_Yimk2));// resR = (float)(Yi - Yimk)

  for(int P = I-K; --P >= Pmin;){
    __m512 v_LRp = _mm512_setzero_ps(); // LRp = 0.0

    for(int T = 0; T <= Kmax[P]; T++){
      __m512 v_Ivar;
      __m512d Ypt = _mm512_set1_pd(Yc(Y,P,T));
      __m512 v_deltaY = DoubleToFloat_mm512(_mm512_sub_pd(v_Yik1, Ypt), _mm512_sub_pd(v_Yik2, Ypt));// deltaY = (float)(Yik - Yc(Y,P,T))
      __m512 v_resL = _mm512_set1_ps((float)(Y[P] - Y[P-T])); // resL = Y[P] - Y[P-T]
      __m512d Yp = _mm512_set1_pd(Y[P]);
      __m512 v_Yint =  DoubleToFloat_mm512(_mm512_sub_pd(v_Yimk1, Yp), _mm512_sub_pd(v_Yimk2, Yp));// Yint = (float)(Yimk - Y[P])
      
      __m512 v_LRik = FAIKP_mm512(v_deltaY, v_resR, v_resL, v_Yint, v_Ivar);

      int QminP = (P >= Pmin2) ? max(Qmin,Jmin[P]) : max(Qmin2, Jmin[P]);
      MFLOAT2 *ALpt = &A.AL(P,T,0);

      __m512 LogFnPowIKP = _mm512_set1_ps(LogFnPow[I-K-P-1]);

      __m512d PrBiasWTp = _mm512_set1_pd(PrBiasWT[P]);
      __m512 v_PrBiasWTikp = DoubleToFloat_mm512(_mm512_sub_pd(v_PrBiasWTik1, PrBiasWTp), _mm512_sub_pd(v_PrBiasWTik2, PrBiasWTp));// PrBiasWTikp = PrBiasWTik - PrBiasWT[P]

      // Critical Loop
      for(int Q= min(J-1,Jmax[P]); Q >= QminP; Q--){
	__m512 v_LRpq = _mm512_mul_ps(_mm512_set1_ps(ALpt[Q]), 
				      FAJ_mm512<1,OUT_TYPE,OUT_OLD>((float)(Xm[J]-Xm[Q]), v_deltaY, J-Q, I-K-P, FBiasWTm[J-Q], LogFpPow[J-Q-1], LogFnPowIKP, v_PrBiasWTikp, v_LRik, v_Ivar));
	v_LRp = _mm512_add_ps(v_LRp, v_LRpq);

#if VECTOR2_DEBUG
	if(TRACE){
	  float LRpq[16], LRp[16], deltaY[16], LRik[16], Ivar[16];
	  _mm512_storeu_ps(LRpq, v_LRpq);
	  _mm512_storeu_ps(LRp, v_LRp);
	  _mm512_storeu_pd(AL, v_AL1);
	  _mm512_storeu_pd(&AL[8], v_AL2);
	  _mm512_storeu_ps(deltaY, v_deltaY);
	  _mm512_storeu_ps(LRik, v_LRik);
	  _mm512_storeu_ps(Ivar, v_Ivar);

	  printf("\t P=%d,T=%d,Q=%d,h=%d:ALpt[Q]= %0.6e,x=%0.4f,y=%0.4f,LRik=%0.6e,Ivar=%0.6e:LRp= %0.8e, ScaleLI[P]= %0.8e, AL= %0.8e (delta= %0.8e)\n", 
		 P, T, Q, HTRACE, ALpt[Q], Xm[J]-Xm[Q], deltaY[HTRACE], LRik[HTRACE], Ivar[HTRACE], LRp[HTRACE], ScaleLI[P], LRp[HTRACE] * ScaleLI[P] + AL[HTRACE], LRpq[HTRACE] * ScaleLI[P]);
	  fflush(stdout);
	}
#endif	
      }
    }

    /* convert v_LRp to double precision */
#if USE_MIC
    __m512d v_LRp1 = _mm512_cvtpslo_pd(v_LRp);/* low order 8 floats to double */
    __m512d v_LRp2 = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_LRp, 0x2 | (0x3 << 2)));/* high order 8 floats to double */
#elif USE_AVX512
    __m512d v_LRp1 = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_castsi512_ps(_mm512_castsi256_si512(_mm512_extracti64x4_epi64(_mm512_castps_si512(v_LRp),0)))));/* low order 8 floats to double */
    __m512d v_LRp2 = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_castsi512_ps(_mm512_castsi256_si512(_mm512_extracti64x4_epi64(_mm512_castps_si512(v_LRp),1)))));/* high order 8 floats to double */
#endif

    // AL += LRp * ScaleLI[P]
    if(DEBUG>=2) assert(IMIN <= P && P <= IMAX); // WAS assert(IMIN-1 <= P && P <= IMAX+1);
    if(!isfinite(ScaleLI[P])){/* NEW144 : check if all LRp[] values are 0 */
      double LRp[16];
      _mm512_store_pd(LRp, v_LRp1);
      _mm512_store_pd(&LRp[8], v_LRp2);
      double maxLRp = 0.0;
      for(int h = 0; h < 16; h++)
	maxLRp = max(maxLRp, LRp[h]);
      if(maxLRp > 0.0){
	printf("ALq_mm512:m=%d,I=%d,K=%d,J=%d,P=%d:t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d:\n",m,I,K,J,P,t,M,N,IMIN,IMAX);
	printf("\t ScaleLI[P]= %0.17e, LogScaleL[P]= %0.17f, LogScaleL[I]= %0.17f, maxLRp= %0.17e\n",ScaleLI[P], LogScaleL[P], LogScaleL[I],maxLRp);
	fflush(stdout);
	assert(isfinite(ScaleLI[P]));
      }
    } else {
      if(DEBUG>=2) assert(fabs(ScaleLI[P] - exp(LogScaleL[P] - LogScaleL[I])) < fabs(ScaleLI[P]) * 1e-5 + 1e-6);
      __m512d ScaleLIP = _mm512_set1_pd(ScaleLI[P]);
      v_AL1 = _mm512_add_pd(v_AL1, _mm512_mul_pd(v_LRp1, ScaleLIP));
      v_AL2 = _mm512_add_pd(v_AL2, _mm512_mul_pd(v_LRp2, ScaleLIP));
    }
  }

  /* return results */
  _mm512_store_pd(AL, v_AL1);
  _mm512_store_pd(&AL[8], v_AL2);

#if DEBUG // HERE >= 2
  for(int h = 0; h < 16; h++){
#if VECTOR2_DEBUG      
    if(DEBUG && !(origUmax < 0 ? AL[h] - ALqh[h] > -fabs(ALqh[h]) * 1e-4 - 1e-9 
		  : fabs(ALqh[h] - AL[h]) < fabs(ALqh[h]) * 1e-4 + 1e-9)){
#pragma omp critical
      {
	printf("ALq_mm512(m=%d,I=%d,K=%d,J=%d,U=%d..%d,Ub=%d,Yik[h]=%0.4f,Yt[h]=%0.4f,t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d:\n",m,I,K,J,origU,origUmax,origUb,Yik[h],Yt[h],t,M,N,IMIN,IMAX);
	printf("\t h=%d: AL[h]= %0.8e, ALqh[h]= %0.8e (final U=%d..%d,Ub=%d)\n",h,AL[h],ALqh[h],finalUmin,finalUmax,finalUb);
	fflush(stdout);
	if(origUmax < 0)
	  assert(AL[h] - ALqh[h] > -(fabs(ALqh[h]) * 1e-4 - 1e-9));
	else
	  assert(fabs(ALqh[h] - AL[h]) < fabs(ALqh[h]) * 1e-4 + 1e-9);
      }
    }
#endif
    assert(isfinite(AL[h]) && AL[h] >= 0.0);
  }
#endif
}

#endif // VECTOR_MIC || VECTOR_AVX

/** compute AR[I][K][J] (see Equation 1 : does not include PM(I,K) */
/** NOTE : When called to compute ART[K][J]: Y[I] is new site and can be Y[0] */
static double AR(int m, /* for debugging only */
		 int I,
		 int K,
		 int J,
		 int U, /**< initialized to Umin */
		 int Umax,/**< If Umax < U, U range is unknown: IFF Umax== -1, Update A[I][K][J].URmin,URmax */
		 double Yik,
		 double *  Xm,
		 double *  Y,
		 int M,
		 int N,
		 int rc,
		 double Scale,/* exp(-LogScaleR[I]) */
		 double *  ScaleRI,/* scaleRI[P] = exp(LogScaleR[P] - LogScaleR[I]) */
		 double *LogScaleR, /* only used for debugging */
		 QB &  A,
		 int IMIN,
		 int IMAX,
		 int *  Kmax,
		 int *  Jmin,
		 int *  Jmax,
		 const int OUT_TYPE,
		 const int OUT_OLD)
{
  double UR = PoutlierEnd;
  int JMIN = (EXTEND && I >= N+1 - DELTA_YU) ? (M+1 - DELTA_XU - EXTEND_MAX) : (M+1 - DELTA_XU);
  if(J >= JMIN){
    MFLOAT URu;
    if(DEBUG/* HERE >=2 */) assert(U > Umax);// NOTE : Ub code below is missing, see ARq() to fix it if assertion triggers
    if(U <= Umax){
      if(DEBUG>=2) assert(U > I);
      if(DEBUG>=2) assert(Umax <= N+1);
      URu = FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik, Y[U]-Y[U-1], Xm,Y,M,N,END_FIX ? rc : 1);
      while(++U <= Umax){
	URu += FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik, Y[U]-Y[U-1], Xm,Y,M,N, END_FIX ? rc : 1);
	if(rc && U==N+1)
	  URu += FRE(I,K,J,N+1,N,Xm[M+1]-Xm[J],Y[U]-Yik,M);
      }
    } else {/* need to determine U range */
      int save = (Umax == -1) ? 1 : 0;
      int Rijk = (IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]) ? A.Rij(I,K,J) : N+1;
      if(DEBUG>=2)assert(Rijk >= I);
      /* decrement Rijk until Y[Rijk] overlaps X */
      while(Rijk > I && Y[Rijk] - Yik > Xm[M+1]-Xm[J])
	Rijk--;
      /* increment Rijk until Y[Rijk+1] no longer overlaps X */
      while(Rijk <= N && Y[Rijk+1] - Yik <= Xm[M+1]-Xm[J])
	Rijk++;
      int Umin,Ub;
      Ub = Umin = Umax = min(N,Rijk)+1;
      URu = FR<0>(I,K,J,Umin,Umin-1,Yik,Xm[M+1]-Xm[J],Y[Umin]-Yik,Y[Umin]-Y[Umin-1],Xm,Y,M,N,END_FIX ? rc : 1);
      MFLOAT bUR = URu;

      for(U=Umin; --U > I;){
	MFLOAT FRu = ERFC_FIX ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)
	  : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
	if(FRu < URu * EPS)
	  break;
	if(URANGE_FIX && FRu > bUR){
	  Umax = U;
	  URu = 0.0;
	  bUR = FRu;
	}
	Umin = U;
	URu += FRu;
      }

      for(U = Umax; ++U <= N+1;){
	MFLOAT FRu = (ERFC_FIX && U < Ub) ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)
	  : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
	if(FRu < URu * EPS)
	  break;
	Umax = U;
	URu += FRu;
      }

      if(ISLINEAR && rc){
	MFLOAT FRu = FRE(I,K,J,N+1,N,Xm[M+1]-Xm[J],Y[N+1]-Yik,M);
	if(FRu >= URu * EPS){
	  Umax = N+1;
	  URu += FRu;
	}
      }

      if(save){
	if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]);
	A.URmin(I,K,J) = Umin;
	A.URmax(I,K,J) = Umax;
      }
    } /* unknown U range */
    UR += PoutlierEndM1 * URu;
  } /* J <= JMIN */

  if(DEBUG>=2) assert(isfinite(UR) && UR >= 0.0);
  double AR = UR * Scale;
  if(DEBUG>=2) assert(isfinite(AR) && AR >= 0.0);
  int Pmax = min(IMAX, I + deltaExtYRef);
  int Pmax2 = min(IMAX, I + DELTA_Y);
  int Qmax = min(M, J + deltaExtXRef);
  int Qmax2 = min(M, J + DELTA_X);
  if(DEBUG>=2 && !(I+1 >= IMIN)){
    #pragma omp critical
    {
      printf("\n m=%d,I=%d,K=%d,J=%d:IMIN=%d,IMAX=%d,N=%d\n",m,I,K,J,IMIN,IMAX,N);
      fflush(stdout);
      assert(I+1 >= IMIN);
    }
  }

  double PrBiasWTi;/* value of PrBiasWT[I] is not valid (if K==0), since Y[I] is a new label, but differentials for PrBiasWT[I+1..N+1] remain valid */
  if(1){/* compute correct value of PrBiasWT[I] so differential to PrBiasWT[I+1] is correct */
    double y = Y[I+1] - Y[I];
    double pr = Pr(y);
    if(DEBUG>=2) assert(isfinite(pr) && 0.0 <= pr && pr <= 1.0);
    double prbias = (pr <= 0.0 || pr >= 1.0) ? 0.0 : pr * log(pr) + (1.0 - pr) * log(1.0 - pr);    
    if(DEBUG>=2)  assert(isfinite(prbias) && prbias <= 0.0);
    PrBiasWTi = PrBiasWT[I+1] + PRbiasWT * prbias;
  }

  for(int P= I; ++P <= Pmax + KMAX && P <= IMAX;){
    MFLOAT2 LRp = 0.0;
    int Tmax = min(P-I-1,Kmax[P]);
    for(int T = max(0,P-Pmax); T <= Tmax; T++){
      MFLOAT Ivar;
      MFLOAT deltaY = Yc(Y,P,T) - Yik;

      if(DEBUG>=2) assert(I >= 0);
      if(DEBUG>=2) assert(P <= PRNsiz);
      MFLOAT LRik = FAIKP(deltaY,P-T-I,P,T,K,Ivar,Y);
      if(DEBUG>=2) assert(isfinite(LRik));
      MFLOAT PrBiasWTpti = PrBiasWT[P-T] - PrBiasWTi;

      int QmaxP = (P - T <= Pmax2) ? min(Qmax,Jmax[P]) : min(Qmax2, Jmax[P]);
      MFLOAT2 *ARpt = &A.AR(P,T,0);

      // Critical Loop for vectorization
      MFLOAT2 LRpq = 0.0;

      #pragma omp simd reduction(+:LRpq)
      for(int Q= max(J+1,Jmin[P]); Q <= QmaxP; Q++){
#if DEBUG>=2
        if(DEBUG>=2) assert(A.AR(P,T,Q) >= 0.0);
#endif
        LRpq += ARpt[Q] * FAJ(Xm[Q]-Xm[J], deltaY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWTpti, LRik, Ivar,1,OUT_TYPE,OUT_OLD);
#if DEBUG>=2
	if(DEBUG>=2 && !(isfinite(LRpq) && LRpq >= 0.0)){
	  printf("I=%d,K=%d,J=%d,P=%d,Q=%d:LRpq=%0.10e:q->AR=%0.10e,LRik=%0.10e,FAJ=%0.10e,x=%0.3f,y=%0.3f,Ivar=%0.10e\n",
		   I,K,J,P,Q,LRpq,A.AR(P,T,Q),LRik,FAJ(Xm[Q]-Xm[J], deltaY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWTpti, LRik, Ivar, 1, OUT_TYPE, OUT_OLD), Xm[Q]-Xm[J],deltaY,Ivar);
	  fflush(stdout);
	  assert(isfinite(LRpq) && LRpq >= 0.0);
	}
#endif
      }
      LRp += LRpq;
    }

    if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
    if(DEBUG>=2) assert(isfinite(LRp) && LRp >= 0.0);
    if(LRp <= 0.0)// NEW100 : ScaleRI[P] may be Inf in this case
      continue;

    if(isfinite(ScaleRI[P]) || LRp != 0.0){// NEW144 : otherwise assume ScaleRI[P] * LRp == 0
      if(DEBUG>=2) assert(isfinite(ScaleRI[P]) && fabs(ScaleRI[P] - exp(LogScaleR[P]-LogScaleR[I])) < fabs(ScaleRI[P]) * 1e-5 + 1e-6);
      AR += LRp * ScaleRI[P]; // WAS exp(LogScaleR[P] - LogScaleR[I]);
      if(DEBUG>=2) assert(isfinite(AR) && AR >= 0.0);
    }
  }
  if(DEBUG && !(isfinite(AR) && AR >= 0.0)){
    printf("\nI=%d,K=%d,J=%d:AR=%0.10e\n",I,K,J,(double)AR);
    fflush(stdout);
    assert(isfinite(AR) && AR >= 0.0);
  }
  return AR;
}

/* specialized version of AR() that never updates URmin,URmax (Umax != -1) */
template<int OUT_TYPE, int OUT_OLD>
static inline double ARq(int m, /* for debugging only */
			 int I,
			 int K,
			 int J,
			 int U, /* initialized to Umin */
			 int Ub,/* estimated peak of range U .. Umax (only used with ERFC_FIX if U <= Umax and can be outside the range U .. Umax) */
			 int Umax,/* If Umax < U, U range is unknown */
			 double Yik,
			 double Yt,/* new value of Y[t] */
			 int t,/* t is either I or I-K */
			 double *  Xm,
			 double *  Y,// NOTE : ARq depends on Y[t] (see FAIKP), but should use new value Yt instead
			 int M,
			 int N,
			 int rc,
			 double Scale,/* exp(-LogScaleR[I]) */
			 double *  ScaleRI,/* scaleRI[P] = exp(LogScaleR[P] - LogScaleR[I]) */
			 double *LogScaleR, /* only used for debugging */
			 QB &  A,
			 int IMIN,
			 int IMAX,
			 int *  Kmax,
			 int *  Jmin,
			 int *  Jmax)
{
#if VECTOR2_DEBUG
  int origU = U, origUmax = Umax, origUb = Ub, finalUmin = -1, finalUmax = -1, finalUb = -1;
#endif
  double UR = PoutlierEnd;
  int JMIN = (EXTEND && I >= N+1 - DELTA_YU) ? (M+1 -DELTA_XU  - EXTEND_MAX) : (M+1 - DELTA_XU);
  if(J >= JMIN){
    MFLOAT URu;
    if(U <= Umax){
      if(DEBUG>=2) assert(U > I);
      if(DEBUG>=2) assert(Umax <= N+1);
      URu = 0.0f;
      if(rc && Umax==N+1)
	URu += FRE(I,K,J,N+1,N,Xm[M+1]-Xm[J],Y[Umax]-Yik,M);
      for(; U <= Umax; U++)
	URu += (ERFC_FIX && U < Ub) ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik, Y[U]-Y[U-1], Xm,Y,M,N, END_FIX ? rc : 1)
	  : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik, Y[U]-Y[U-1], Xm,Y,M,N, END_FIX ? rc : 1);
    } else {/* need to determine U range */
      int Rijk = (IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]) ? A.Rij(I,K,J) : N+1;
      if(DEBUG>=2)assert(Rijk >= I);
      /* decrement Rijk until Y[Rijk] overlaps X */
      while(Rijk > I && Y[Rijk] - Yik > Xm[M+1]-Xm[J])
	Rijk--;
      /* increment Rijk until Y[Rijk+1] no longer overlaps X */
      while(Rijk <= N && Y[Rijk+1] - Yik <= Xm[M+1]-Xm[J])
	Rijk++;
      int Umin;
      Ub = Umin = Umax = min(N,Rijk)+1;

      URu = FR<0>(I,K,J,Umin,Umin-1,Yik,Xm[M+1]-Xm[J],Y[Umin]-Yik,Y[Umin]-Y[Umin-1],Xm,Y,M,N,END_FIX ? rc : 1);
      MFLOAT bUR = URu;
#if VECTOR2_DEBUG
	if(TRACE){
	  printf("\t ARq:m=%d,I=%d,K=%d,J=%d,U=%d(Umin=%d,Umax=%d,Ub=%d,Rijk=%d):URu= %0.8e\n",
	      m,I,K,J,U,Umin,Umax,Ub,Rijk,URu);
	  fflush(stdout);
        }
#endif

      for(int U=Umin; --U > I;){
	MFLOAT FRu = ERFC_FIX ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)
	  : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
#if VECTOR2_DEBUG
	if(TRACE){
	  printf("\t ARq:m=%d,I=%d,K=%d,J=%d,U=%d(Umin=%d,Umax=%d,Ub=%d):FRu= %0.8e, URu= %0.8e\n",
	    m,I,K,J,U,Umin,Umax,Ub,FRu,URu);
	  fflush(stdout);
        }
#endif
	if(FRu < URu * EPS)
	  break;
	if(URANGE_FIX && FRu > bUR){
	  Umax = U;
	  URu = 0.0;
	  bUR = FRu;
	}
	Umin = U;
	URu += FRu;
      }
      for(int U = Umax; ++U <= N+1;){
	MFLOAT FRu = (ERFC_FIX && U < Ub) ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)
	  : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
#if VECTOR2_DEBUG
	if(TRACE){
	  printf("\t ARq:m=%d,I=%d,K=%d,J=%d,U=%d(Umin=%d,Umax=%d,Ub=%d):FRu= %0.8e, URu= %0.8e\n",
	    m,I,K,J,U,Umin,Umax,Ub,FRu,URu);
	  fflush(stdout);
        }
#endif
	if(FRu < URu * EPS)
	  break;
	Umax = U;
	URu += FRu;
      }
      if(ISLINEAR && rc){
	MFLOAT FRu = FRE(I,K,J,N+1,N,Xm[M+1]-Xm[J],Y[N+1]-Yik,M);
	if(FRu >= URu * EPS){
	  Umax = N+1;
	  URu += FRu;
	}
      }
#if VECTOR2_DEBUG
      finalUmin = Umin;
      finalUmax = Umax;
      finalUb = Ub;
#endif
    } /* unknown U range */
    UR += PoutlierEndM1 * URu;
  } /* J <= JMIN */

  if(DEBUG) assert(isfinite(UR) && UR >= 0.0);
  double AR = UR * Scale;

#if VECTOR2_DEBUG
  if(TRACE && fabs(Yt - YTRACE) < 0.001){
    printf("\nARq:m=%d,I=%d,K=%d,J=%d,JMIN=%d,U=%d..%d, Ub=%d(final:U=%d..%d,Ub=%d), Yik= %0.6f, Yt= %0.6f, t=%d: UR= %0.8e, Scale= %0.8e, UR * Scale= %0.8e\n",
	   m,I,K,J,JMIN,origU,origUmax,origUb,finalUmin,finalUmax,finalUb,Yik,Yt,t, UR, Scale, AR);
    fflush(stdout);
  }
#endif  

  if(DEBUG) assert(isfinite(AR) && AR >= 0.0);
  int Pmax = min(IMAX, I + deltaExtYRef);
  int Pmax2 = min(IMAX, I + DELTA_Y);
  int Qmax = min(M, J + deltaExtXRef);
  int Qmax2 = min(M, J + DELTA_X);
  if(DEBUG>=2 && !(I+1 >= IMIN)){
    #pragma omp critical
    {
      printf("\n m=%d,I=%d,K=%d,J=%d:IMIN=%d,IMAX=%d,N=%d\n",m,I,K,J,IMIN,IMAX,N);
      fflush(stdout);
      assert(I+1 >= IMIN);
    }
  }

  double Yi,Yimk;
  if(t == I){
    Yi = Yt;
    Yimk = K ? Y[I-K] : Yt;
  } else {
    if(DEBUG/* HERE >=2 */) assert(K > 0 && t == I - K);
    Yi = Y[I];
    Yimk = Yt;
  }

  double PrBiasWTi;/* value of PrBiasWT[I] may not be valid, since Y[t] is a new label, but differentials for PrBiasWT[t+1..N+1] remain valid */
  if(t==I){/* compute correct value of PrBiasWT[I] (relative to PrBiasWT[t+1..N+1]) */
    double y = Y[t+1] - Yt;
    double pr = Pr(y);
    if(DEBUG>=2) assert(isfinite(pr) && 0.0 <= pr && pr <= 1.0);
    double prbias = (pr <= 0.0 || pr >= 1.0) ? 0.0 : pr * log(pr) + (1.0 - pr) * log(1.0 - pr);    
    if(DEBUG>=2)  assert(isfinite(prbias) && prbias <= 0.0);
    PrBiasWTi = PrBiasWT[t+1] + PRbiasWT * prbias;
  } else
    PrBiasWTi = PrBiasWT[I];

  MFLOAT resL = Yi - Yimk;

#if VECTOR2_DEBUG
  if(TRACE && fabs(Yt - YTRACE) < 0.001){
    printf("ARq:m=%d,I=%d,K=%d,J=%d,Pmax=%d,%d,Qmax=%d,%d\n",m,I,K,J,Pmax,Pmax2,Qmax,Qmax2);
    fflush(stdout);
  }
#endif

  for(int P= I; ++P <= Pmax + KMAX && P <= IMAX;){
    MFLOAT2 LRp = 0.0;
    int Tmax = min(P-I-1,Kmax[P]);
    for(int T = max(0,P-Pmax); T <= Tmax; T++){
      MFLOAT Ivar;
      MFLOAT deltaY = Yc(Y,P,T) - Yik;
      MFLOAT resR = Y[P] - Y[P-T];
      MFLOAT Yint = Y[P-T] - Yi;

      if(DEBUG>=2) assert(I >= 0);
      if(DEBUG>=2) assert(P <= PRNsiz);
      MFLOAT LRik = FAIKP(deltaY, resR, resL, Yint, P-T-I, P, T, K, Ivar, Y);/* NOTE : does NOT index Y[] */
      if(DEBUG>=2) assert(isfinite(LRik));

      MFLOAT PrBiasWTpti = PrBiasWT[P-T] - PrBiasWTi;

      int QmaxP = (P - T <= Pmax2) ? min(Qmax,Jmax[P]) : min(Qmax2, Jmax[P]);
      MFLOAT2 *ARpt = &A.AR(P,T,0);

#if VECTOR2_DEBUG
      if(TRACE && fabs(Yt - YTRACE) < 0.001){
	printf("\t ARq:m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d:Jmin[P]=%d,Jmax[P]=%d\n",m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxP,Jmin[P],Jmax[P]);
	fflush(stdout);
      }
#endif
      // Critical Loop for vectorization
      MFLOAT2 LRpq = 0.0;

      #pragma omp simd reduction(+:LRpq)
      for(int Q= max(J+1,Jmin[P]); Q <= QmaxP; Q++){
#if DEBUG>=2
        if(DEBUG>=2) assert(A.AR(P,T,Q) >= 0.0);
#endif
	MFLOAT FAj = FAJ(Xm[Q]-Xm[J], deltaY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWTpti, LRik, Ivar, 1, OUT_TYPE, OUT_OLD);
        LRpq += ARpt[Q] * FAj;

	if(VERB>=2){
	  printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:x=%0.4f,y=%0.4f:LRik= %0.8e,Ivar= %0.8e,FAj= %0.8e,ARpt[Q]= %0.8e, LRp= %0.8e, ScaleRI[P]= %0.8e, AR= %0.8e (delta= %0.8e)\n",
		 m,I,K,J,P,T,Q,Xm[Q]-Xm[J],deltaY,LRik,Ivar,FAj,ARpt[Q], LRpq+LRp, ScaleRI[P], AR + (LRpq + LRp) * ScaleRI[P], ARpt[Q] * FAj * ScaleRI[P]);
	  fflush(stdout);
	}

#if VECTOR2_DEBUG
	if(TRACE && fabs(Yt - YTRACE) < 0.001){
          printf("\t ARq:P=%d,T=%d,Q=%d: LRp= %0.8e, ScaleRI[P]= %0.8e (AR= %0.8e)\n", P, T, Q, LRp + LRpq, ScaleRI[P], (LRp + LRpq) * ScaleRI[P] + AR);
	  fflush(stdout);
        }
#endif

#if DEBUG>=2
	if(DEBUG>=2 && !(isfinite(LRpq) && LRpq >= 0.0)){
	  printf("I=%d,K=%d,J=%d,P=%d,Q=%d:LRpq=%0.10e:q->AR=%0.10e,LRik=%0.10e,FAJ=%0.10e,x=%0.3f,y=%0.3f,Ivar=%0.10e\n",
		   I,K,J,P,Q,LRpq,A.AR(P,T,Q),LRik,FAJ(Xm[Q]-Xm[J], deltaY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], LRik, PrBiasWTpti, Ivar, 1, OUT_TYPE, OUT_OLD), Xm[Q]-Xm[J],deltaY,Ivar);
	  fflush(stdout);
	  assert(isfinite(LRpq) && LRpq >= 0.0);
	}
#endif
      }
      LRp += LRpq;
    }

    if(DEBUG>=2) assert(IMIN <= P && P <= IMAX); // WAS assert(IMIN-1 <= P && P <= IMAX+1);
    if(isfinite(ScaleRI[P]) || LRp != 0.0){// NEW144 : otherwise assume ScaleRI[P] * LRp == 0
      if(DEBUG>=2 && !(isfinite(ScaleRI[P]) && fabs(ScaleRI[P] - exp(LogScaleR[P]-LogScaleR[I])) < fabs(ScaleRI[P]) * 1e-5 + 1e-6)){
        #pragma omp critical
        {
	  printf("ARq:m=%d,I=%d,K=%d,J=%d,P=%d:IMIN=%d,IMAX=%d,ScaleRI[P]= %0.10e, LogScaleR[P]= %0.10f, LogScaleR[I]= %0.10f, exp(LogScaleR[P] - LogScaleR[I])= %0.10e\n",
	       m,I,K,J,P,IMIN,IMAX,ScaleRI[P], LogScaleR[P], LogScaleR[I], exp(LogScaleR[P] - LogScaleR[I]));
	  for(int i = IMAX; i >= IMIN; i--)
	    printf("\t LogScaleR[%d]= %0.10f\n", i, LogScaleR[i]);
	  fflush(stdout);
	  assert(isfinite(ScaleRI[P]) && fabs(ScaleRI[P] - exp(LogScaleR[P]-LogScaleR[I])) < fabs(ScaleRI[P]) * 1e-5 + 1e-6);
        }
      }
      AR += LRp * ScaleRI[P]; // WAS exp(LogScaleR[P] - LogScaleR[I]);
      if(DEBUG>=2 && !(isfinite(AR) && AR >= 0.0)){
	printf("ARq:m=%d,I=%d,K=%d,J=%d,P=%d:IMIN=%d,IMAX=%d,ScaleRI[P]= %0.10e, LogScaleR[P]= %0.10f, LogScaleR[I]= %0.10f, exp(LogScaleR[P] - LogScaleR[I])= %0.10e\n",
	  m,I,K,J,P,IMIN,IMAX,ScaleRI[P], LogScaleR[P], LogScaleR[I], exp(LogScaleR[P] - LogScaleR[I]));
	printf("\t LRp= %0.17e, ScaleRI[P]= %0.17e, AR= %0.17e\n",LRp,ScaleRI[P],AR);
	fflush(stdout);
	assert(isfinite(AR) && AR >= 0.0);
      }
    }
  }
  if(DEBUG && !(isfinite(AR) && AR >= 0.0)){
    printf("\nI=%d,K=%d,J=%d:AR=%0.10e\n",I,K,J,(double)AR);
    fflush(stdout);
    assert(isfinite(AR) && AR >= 0.0);
  }
  return AR;
}

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT

/* AVX256 version of ARq */
template<int OUT_TYPE, int OUT_OLD>
static inline void ARq_mm256(double *AR,/* return values AR[0..7] */
			     int m, /* for debugging only */
			     int I,
			     int K,
			     int J,
			     int U, /* initialized to Umin */
			     int Ub,/* estimated peak of range U .. Umax (only used with ERFC_FIX if U <= Umax and can be outside the range U .. Umax) */
			     int Umax,/* If Umax < U, U range is unknown */
			     __m256d v_Yik1, __m256d v_Yik2, /* 8 New values for Yik */
			     __m256d v_Yt1, __m256d v_Yt2, /* 8 New values for Y[t] */
			     int t,/* t is either I or I-K */
			     double *  Xm,
			     double *  Y,// NOTE : ARq depends on Y[t] (see FAIKP), but should never use Y[t] directly
			     int M,
			     int N,
			     int rc,
			     double Scale,/* exp(-LogScaleR[I]) */
			     double *  ScaleRI,/* scaleRI[P] = exp(LogScaleR[P] - LogScaleR[I]) */
			     double * LogScaleR, /* for debugging only */
			     QB &  A,
			     int IMIN,
			     int IMAX,
			     int *  Kmax,
			     int *  Jmin,
			     int *  Jmax)
{
#if DEBUG>=2 || VECTOR2_DEBUG
  int origU = U, origUmax = Umax, origUb = Ub, finalUmin = -1, finalUmax = -1, finalUb = -1 ;
  double Yik[8];
  double Yt[8];
  _mm256_storeu_pd(Yik, v_Yik1);  _mm256_storeu_pd(&Yik[4], v_Yik2);
  _mm256_storeu_pd(Yt, v_Yt1); _mm256_storeu_pd(&Yt[4], v_Yt2);
#if VECTOR2_DEBUG
  double ARqh[8];
  for(int h = 0; h < 8; h++)
    ARqh[h] = ARq<OUT_TYPE,OUT_OLD>(m,I,K,J,U,Ub,Umax,Yik[h],Yt[h],t,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax);
#endif
#endif  

  __m256d v_UR1 = _mm256_set1_pd(PoutlierEnd);// UR = PoutlierEnd
  __m256d v_UR2 = v_UR1;

  int JMIN = (EXTEND && I >= N+1-DELTA_Y) ? (M+1 -DELTA_X-EXTEND_MAX) : (M+1 -DELTA_X);

  if(J >= JMIN){
    // NOTE : U range is U .. Umax, with Ub corresponding to Rij+1, the estimated peak value (can be outside the range U..Umax)

    if(U > Umax){/* need to determine U range : Use first and last Yik values to compute upper bound of range */
      double Yik8[VEC_SIZ] __attribute__((aligned(32)));
      _mm256_store_pd(Yik8, v_Yik1);
      _mm256_store_pd(&Yik8[4], v_Yik2);
      int Ulb = N+1, Uub = I, Umb = -1 ;/* U range will be computed Ulb .. Uub, with Umb corresponding to Rij+1, the estimated peak value (can lie outside the range Ulb .. Uub) */
      
      for(int T = 0; T <= 1; T++){
	double Yik = T ? Yik8[7] : Yik8[0];

	int Rijk = (IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]) ? A.Rij(I,K,J) : N+1;

	/* decrement Rijk until Y[Rijk] overlaps X */
	while(Rijk > I && Y[Rijk] - Yik > Xm[M+1]-Xm[J])
	  Rijk--;
	/* increment Rijk until Y[Rijk+1] no longer overlaps X */
	while(Rijk <= N && Y[Rijk+1] - Yik <= Xm[M+1]-Xm[J])
	  Rijk++;
	int Umid, Umin;
	Umid = Umin = Umax = min(N,Rijk)+1;

	MFLOAT URu = FR<0>(I,K,J,Umin,Umin-1,Yik,Xm[M+1]-Xm[J],Y[Umin]-Yik,Y[Umin]-Y[Umin-1],Xm,Y,M,N,END_FIX ? rc : 1);
	MFLOAT bUR = URu;
#if VECTOR2_DEBUG
	if(TRACE){
	  printf("\t ARq_mm256:m=%d,I=%d,K=%d,J=%d,T=%d,Yik=%0.6f,U=%d(Umin=%d,Umax=%d,Ub=%d,Rijk=%d):URu= %0.8e\n",
		 m,I,K,J,T,Yik,U,Umin,Umax,Umid,Rijk,URu);
	  fflush(stdout);
        }
#endif

	for(U=Umin; --U > I;){
	  MFLOAT FRu = ERFC_FIX ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)
	    : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
#if VECTOR2_DEBUG
	  if(TRACE){
	    printf("\t ARq_mm256:m=%d,I=%d,K=%d,J=%d,T=%d,Yik=%0.6f,U=%d(Umin=%d,Umax=%d,Ub=%d):FRu= %0.8e,URu= %0.8e\n",
		   m,I,K,J,T,Yik,U,Umin,Umax,Ub,FRu,URu);
	    fflush(stdout);
	  }
#endif
	  if(FRu < URu * EPS)
	    break;
	  if(URANGE_FIX && FRu > bUR){
	    Umax = U;
	    URu = 0.0;
	    bUR = FRu;
	  }
	  Umin = U;
	  URu += FRu;
	}

	for(U = Umax; ++U <= N+1;){
	  MFLOAT FRu = (ERFC_FIX && U < Umid) ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)
	    : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
#if VECTOR2_DEBUG
	  if(TRACE){
	    printf("\t ARq_mm256:m=%d,I=%d,K=%d,J=%d,T=%d,Yik=%0.6f,U=%d(Umin=%d,Umax=%d,Ub=%d):FRu= %0.8e,URu= %0.8e\n",
		   m,I,K,J,T,Yik,U,Umin,Umax,Ub,FRu,URu);
	    fflush(stdout);
	  }
#endif
	  if(FRu < URu * EPS)
	    break;
	  Umax = U;
	  URu += FRu;
	}

	if(ISLINEAR && rc){
	  MFLOAT FRu = FRE(I,K,J,N+1,N,Xm[M+1]-Xm[J],Y[N+1]-Yik,M);
	  if(FRu >= URu * EPS){
	    Umax = N+1;
	    URu += FRu;
	  }
	}

	if(DEBUG>=2 && !(I < Umin && Umin <= Umax && Umax <= N+1)){
	  printf("ARq_mm256:m=%d,I=%d,K=%d,J=%d,N=%d:T=%d, Yik= %0.4f, Rijk= %d, Umin=%d,Umax=%d\n",
		 m,I,K,J,N,T,Yik, Rijk, Umin, Umax);
	  fflush(stdout);
	  assert(I < Umin && Umin <= Umax && Umax <= N+1);
	}
	Ulb = min(Ulb, Umin);
	Uub = max(Uub, Umax);
	Umb = T ? (Umb + Umid)/2 : Umid;
#if VECTOR2_DEBUG
	if(TRACE){
	  printf("ARq_mm256:T=%d,Yik= %0.6f,Umin=%d,Umax=%d,Umid=%d: Ulb=%d,Uub=%d,Umb=%d\n",	 T,Yik,Umin,Umax,Umid,Ulb,Uub,Umb);
	  fflush(stdout);
	}
#endif
      } // T = 0..1
      
      if(DEBUG) assert(Ulb <= Uub);
      U = Ulb;
      Ub = Umb;
      Umax = Uub;
#if DEBUG>=2 || VECTOR2_DEBUG
      finalUmin = U;
      finalUmax = Umax;
      finalUb = Ub;
#endif
    } // U > Umax

    if(DEBUG>=2 && !(I < U && U <= Umax && Umax <= N+1)){
      printf("ARq_mm256:m=%d,I=%d,K=%d,J=%d,N=%d: U=%d,Umax=%d,Ub=%d\n",
	     m,I,K,J,N, U, Umax,Ub);
      fflush(stdout);
      assert(I < U && U <= Umax && Umax <= N+1);
    }

    __m256 Xj = _mm256_set1_ps((float)(Xm[M+1]-Xm[J]));
    __m256 v_URu = _mm256_setzero_ps();

    for(;U <= Umax; U++){
      __m256d Yu = _mm256_set1_pd(Y[U]);
      __m128 v_Yui1 = _mm256_cvtpd_ps(Yu - v_Yik1);
      __m128 v_Yui2 = _mm256_cvtpd_ps(Yu - v_Yik2);
      __m256 v_Yui = _mm256_set_m128(v_Yui2, v_Yui1); // Yui = Y[U] - Yik

      __m256 v_Yuv = _mm256_set1_ps((float)(Y[U] - Y[U-1])); // Yuv = Y[U] - Y[U-1]

      v_URu += (ERFC_FIX && U < Ub) ? FR_mm256<1>(I,K,J,U,U-1, Xj, v_Yui, v_Yuv, M, N, END_FIX ? rc : 1)
	: FR_mm256<0>(I,K,J,U,U-1, Xj, v_Yui, v_Yuv, M, N, END_FIX ? rc : 1);

      if(rc && U==N+1)
	v_URu += FRE_mm256(I,K,J,N+1,N, Xj, v_Yui, M);
    }

    /* convert v_URu to double precision */
    __m256d v_URu1 = _mm256_cvtps_pd(_mm256_extractf128_ps(v_URu,0));/* low order 4 floats to double */
    __m256d v_URu2 = _mm256_cvtps_pd(_mm256_extractf128_ps(v_URu,1));/* high order 4 floats to double */
    __m256d v_PoutlierEndM1 = _mm256_set1_pd(PoutlierEndM1);    
    v_UR1 += v_URu1 * v_PoutlierEndM1;
    v_UR2 += v_URu2 * v_PoutlierEndM1;// UR += URu * PoutlierEndM1
  } /* J <= JMIN */

#if DEBUG>=2 || VECTOR2_DEBUG
  double UR8[8];
  _mm256_storeu_pd(UR8, v_UR1);
  _mm256_storeu_pd(&UR8[4], v_UR2);
  for(int h = 0; h < 8; h++)
    assert(isfinite(UR8[h]) && UR8[h] >= 0.0);
#endif

  __m256d v_Scale = _mm256_set1_pd(Scale);
  __m256d v_AR1 = v_UR1 * v_Scale;
  __m256d v_AR2 = v_UR2 * v_Scale;// AR = UR * Scale

#if DEBUG>=2 || VECTOR2_DEBUG
    double AR8[8];
    _mm256_storeu_pd(AR8, v_AR1);
    _mm256_storeu_pd(&AR8[4], v_AR2);
    for(int h = 0; h < 8; h++){
      if(!(isfinite(AR8[h]) && AR8[h] >= 0.0)){
	#pragma omp critical
	{
	  printf("\nARq_mm256:m=%d,I=%d,K=%d,J=%d,JMIN=%d,U=%d..%d,Ub=%d(final:U=%d..%d,Ub=%d),h=%d: Yik[h]= %0.6f, Yt[h]= %0.6f, t= %d : UR[h]= %0.8e, Scale= %0.8e, AR[h] = %0.8e\n",
		 m,I,K,J,JMIN,origU,origUmax,origUb,finalUmin,finalUmax,finalUb,h,Yik[h],Yt[h],t,UR8[h], Scale, AR8[h]);
	  fflush(stdout);
	  assert(isfinite(AR8[h]) && AR8[h] >= 0.0);
	}
      }
    }

    if(TRACE){
      printf("ARq_mm256:m=%d,I=%d,K=%d,J=%d,JMIN=%d,U=%d..%d,Ub=%d(final:U=%d..%d,Ub=%d),h=%d: Yik[h]= %0.6f, Yt[h]= %0.6f, t= %d : UR[h]= %0.8e, Scale= %0.8e, UR[h]*Scale = %0.8e\n",
	     m,I,K,J,JMIN,origU,origUmax,origUb,finalUmin,finalUmax,finalUb,HTRACE,Yik[HTRACE],Yt[HTRACE],t,UR8[HTRACE], Scale, AR8[HTRACE]);
      fflush(stdout);
    }
#endif

  int Pmax = min(IMAX, I + deltaExtYRef);
  int Pmax2 = min(IMAX, I + DELTA_Y);
  int Qmax = min(M, J + deltaExtXRef);
  int Qmax2 = min(M, J + DELTA_X);

  __m256d v_Yi1, v_Yi2, v_Yimk1, v_Yimk2;
  if(t == I){
    v_Yi1 = v_Yt1; v_Yi2 = v_Yt2; // Yi = Yt
    if(K){    //    Yimk = K ? Y[I-K] : Yt;
      v_Yimk1 = _mm256_set1_pd(Y[I-K]); v_Yimk2 = v_Yimk1;
    } else {
      v_Yimk1 = v_Yt1; v_Yimk2 = v_Yt2;
    }
  } else {
    if(DEBUG/* HERE >=2 */) assert(K > 0 && t == I - K);
    v_Yi1 = _mm256_set1_pd(Y[I]); v_Yi2 = v_Yi1;// Yi = Y[I]
    v_Yimk1 = v_Yt1; v_Yimk2 = v_Yt2; // Yimk = Yt
  }

  __m256d v_PrBiasWTi1, v_PrBiasWTi2;/* value of PrBiasWT[I] may not be valid, since Y[t] is a new label, but differentials for PrBiasWT[t+1..N+1] remain valid */
  if(t==I){/* compute correct value of PrBiasWT[I] (relative to PrBiasWT[t+1..N+1]) */
    __m256d ytp1 = _mm256_set1_pd(Y[t+1]);
    __m256d v_y1 = ytp1 - v_Yt1, v_y2 = ytp1 - v_Yt2;// y = Y[t+1] - Yt

    __m256d v_pr1 = Pr_mm256d(v_y1), v_pr2 = Pr_mm256d(v_y2);// pr = Pr(y)
    __m256d v_1Mpr1 = v256d_1p0 - v_pr1, v_1Mpr2 = v256d_1p0 - v_pr2;// 1.0 - pr
    // prbias = pr * log(max(pr,1e-300)) + (1.0 - pr) * log(max(1.0 - pr, 1e-300))
    __m256d v_eps = _mm256_set1_pd(1e-300);
    __m256d v_prbias1 = v_pr1 * log_mm256d(_mm256_max_pd(v_pr1, v_eps)) + v_1Mpr1 * log_mm256d(_mm256_max_pd(v_1Mpr1, v_eps));
    __m256d v_prbias2 = v_pr2 * log_mm256d(_mm256_max_pd(v_pr2, v_eps)) + v_1Mpr2 * log_mm256d(_mm256_max_pd(v_1Mpr2, v_eps));

    __m256d v_PrBiasWTtp1 = _mm256_set1_pd(PrBiasWT[t+1]);
    __m256d v_PRbiasWT = _mm256_set1_pd(PRbiasWT);

    v_PrBiasWTi1 = v_PrBiasWTtp1 + v_PRbiasWT * v_prbias1; v_PrBiasWTi2 = v_PrBiasWTtp1 + v_PRbiasWT * v_prbias2;// PrBiasWTi = PrBiasWT[t+1] + PRbiasWT * prbias;
    if(DEBUG>=2){/* check that PrBiasWT[I] <= PrBiasWT[I+1] */
      double PrBiasWTi[8]; _mm256_storeu_pd(PrBiasWTi, v_PrBiasWTi1); _mm256_storeu_pd(&PrBiasWTi[4], v_PrBiasWTi2);
      for(int k = 0; k < 8; k++){
	if(DEBUG && !(PrBiasWTi[k] <= PrBiasWT[t+1])){
	  double Yt[8]; _mm256_storeu_pd(Yt, v_Yt1); _mm256_storeu_pd(&Yt[4], v_Yt2);
	  double Y[8]; _mm256_storeu_pd(Y, v_y1); _mm256_storeu_pd(&Y[4], v_y2);	  
	  double Pr[8]; _mm256_storeu_pd(Pr, v_pr1); _mm256_storeu_pd(&Pr[4], v_pr2);	  
	  double PrM1[8]; _mm256_storeu_pd(PrM1, v_1Mpr1); _mm256_storeu_pd(&PrM1[4], v_1Mpr2);	  
	  double PrBias[8]; _mm256_storeu_pd(PrBias, v_prbias1); _mm256_storeu_pd(&PrBias[4], v_prbias2);	  
	  #pragma omp critical
	  {
	    printf("ARq_mm256:I=%d,K=%d,J=%d,m=%d:t=%d,k=%d:Y[t+1]=%0.6f,Yt[k]=%0.6f,Y[k]=%0.6f,Pr[k]=%0.6f,PrM1[k]=%0.6f,PrBias[k]=%0.6f,PrBiasWTi[k]=%0.6f,PrBiasWT[t+1]=%0.6f,PRbiasWT=%0.3f\n",
		   I,K,J,m,t,k,Y[t+1],Yt[k],Y[k],Pr[k],PrM1[k],PrBias[k],PrBiasWTi[k],PrBiasWT[t+1],PRbiasWT);
	    fflush(stdout);
	    assert(PrBiasWTi[k] <= PrBiasWT[t+1]);
	  }
	}
      }
    }
  } else
    v_PrBiasWTi1 = v_PrBiasWTi2 = _mm256_set1_pd(PrBiasWT[I]);

  __m256 v_resL = _mm256_set_m128(_mm256_cvtpd_ps(v_Yi2 - v_Yimk2), _mm256_cvtpd_ps(v_Yi1 - v_Yimk1));// resL = Yi - Yimk

#if VECTOR2_DEBUG
  if(TRACE){
    printf("ARq_mm256:m=%d,I=%d,K=%d,J=%d,Pmax=%d,%d,Qmax=%d,%d\n",m,I,K,J,Pmax,Pmax2,Qmax,Qmax2);
    fflush(stdout);
  }
#endif

  for(int P= I; ++P <= Pmax + KMAX && P <= IMAX;){
    __m256 v_LRp = _mm256_setzero_ps(); // LRp = 0.0

    int Tmax = min(P-I-1,Kmax[P]);
    for(int T = max(0,P-Pmax); T <= Tmax; T++){
      __m256 v_Ivar;
      __m256d Ypt = _mm256_set1_pd(Yc(Y,P,T));
      __m256 v_deltaY = _mm256_set_m128(_mm256_cvtpd_ps(Ypt - v_Yik2), _mm256_cvtpd_ps(Ypt - v_Yik1));// deltaY = Yc(Y,P,T) - Yik
      __m256 v_resR = _mm256_set1_ps((float)(Y[P] - Y[P-T])); // resR = Y[P] - Y[P-T]
      __m256d Ypmt = _mm256_set1_pd(Y[P-T]);
      __m256 v_Yint =  _mm256_set_m128(_mm256_cvtpd_ps(Ypmt - v_Yi2), _mm256_cvtpd_ps(Ypmt - v_Yi1));// Yint = Y[P-T] - Yi

      __m256 v_LRik = FAIKP_mm256(v_deltaY, v_resR, v_resL, v_Yint, v_Ivar);

      int QmaxP = (P - T <= Pmax2) ? min(Qmax,Jmax[P]) : min(Qmax2, Jmax[P]);
      MFLOAT2 *ARpt = &A.AR(P,T,0);

#if VECTOR2_DEBUG
      if(TRACE){
	printf("\t ARq_mm256:m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d:Jmin[P]=%d,Jmax[P]=%d\n",m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxP,Jmin[P],Jmax[P]);
	fflush(stdout);
      }
#endif
      __m256 v_LogFnPowPTI = _mm256_broadcast_ss(&LogFnPow[P-T-I-1]);

      __m256d v_PrBiasWTpt = _mm256_set1_pd(PrBiasWT[P-T]);
      __m256 v_PrBiasWTpti = _mm256_set_m128(_mm256_cvtpd_ps(v_PrBiasWTpt - v_PrBiasWTi2), _mm256_cvtpd_ps(v_PrBiasWTpt - v_PrBiasWTi1));// PrBiasWTpti = PrBiasWT[P-T] - PrBiasWTi
      if(DEBUG>=2){/* check that PrBiasWTpti >= 0.0 */
	float PrBiasWTpti[8]; _mm256_storeu_ps(PrBiasWTpti, v_PrBiasWTpti);
	for(int k = 0; k < 8; k++){
	  if(DEBUG && !(PrBiasWTpti[k] >= 0.0)){
	    double PrBiasWTi[8]; _mm256_storeu_pd(PrBiasWTi, v_PrBiasWTi1); _mm256_storeu_pd(&PrBiasWTi[4], v_PrBiasWTi2);
	    double PrBiasWTpt[4]; _mm256_storeu_pd(PrBiasWTpt,v_PrBiasWTpt);
	    #pragma omp critical
	    {
	      printf("ARq_mm256:m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d:Jmin[P]=%d,Jmax[P]=%d\n",m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxP,Jmin[P],Jmax[P]);
	      printf("\t k=%d:PrBiasWTi[k]= %0.8f, PrBiasWT[P-T]= %0.8f (PrBiasWTpt[k%%4]= %0.8f), PrBiasWTpti[k]= %0.8f\n",k, PrBiasWTi[k], PrBiasWT[P-T], PrBiasWTpt[k%4], PrBiasWTpti[k]);
	      fflush(stdout);
	      assert(PrBiasWTpti[k] >= 0.0);
	    }
	  }
	}
      }

      // Critical Loop
      for(int Q= max(J+1,Jmin[P]); Q <= QmaxP; Q++){
        v_LRp += _mm256_set1_ps(ARpt[Q]) * FAJ_mm256<1,OUT_TYPE,OUT_OLD>((float)(Xm[Q]-Xm[J]), v_deltaY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], v_LogFnPowPTI, v_PrBiasWTpti, v_LRik, v_Ivar);
#if VECTOR2_DEBUG
	if(TRACE){
	  float LRp[8];
	  _mm256_storeu_ps(LRp, v_LRp);
	  _mm256_storeu_pd(AR, v_AR1);
	  _mm256_storeu_pd(&AR[4], v_AR2);
	  
	  printf("\t ARq_mm256:P=%d,T=%d,Q=%d,t=0: LRp= %0.8e, ScaleRI[P]= %0.8e (AR= %0.8e)\n", P, T, Q, LRp[0], ScaleRI[P], LRp[0] * ScaleRI[P] + AR[0]);
	  fflush(stdout);
	}
#endif	
      }
    }

    /* convert v_LRp to double precision */
    __m256d v_LRp1 = _mm256_cvtps_pd(_mm256_extractf128_ps(v_LRp, 0));/* low order 4 floats to double */
    __m256d v_LRp2 = _mm256_cvtps_pd(_mm256_extractf128_ps(v_LRp, 1));/* high order 4 floats to double */

    /*if(isfinite(ScaleRI[P]) || !_mm256_testz_si256(_mm256_castpd_si256(v_LRp1),_mm256_castpd_si256(v_LRp1)) 
      || !_mm256_testz_si256(_mm256_castpd_si256(v_LRp2),_mm256_castpd_si256(v_LRp2))) */
    if(!isfinite(ScaleRI[P])){/* NEW144 : check if all LRp[] values are 0 */
      double LRp[8];
      _mm256_store_pd(LRp, v_LRp1);
      _mm256_store_pd(&LRp[4], v_LRp2);
      double maxLRp = 0.0;
      for(int h = 0; h < 8; h++)
	maxLRp = max(maxLRp, LRp[h]);
      if(maxLRp > 0.0){
	printf("ARq_mm256:m=%d,I=%d,K=%d,J=%d,P=%d(Pmax=%d,KMAX=%d),t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d:\n",m,I,K,J,P,Pmax,KMAX,t,M,N,IMIN,IMAX);
	printf("\t ScaleRI[P]= %0.17e, LogScaleR[P]= %0.17f, LogScaleR[I]= %0.17f,maxLRp= %0.17e\n",ScaleRI[P], LogScaleR[P], LogScaleR[I],maxLRp);
	fflush(stdout);
	assert(isfinite(ScaleRI[P]));
      }
    } else {
      if(DEBUG>=2) assert(isfinite(ScaleRI[P]) && fabs(ScaleRI[P] - exp(LogScaleR[P]-LogScaleR[I])) < fabs(ScaleRI[P]) * 1e-5 + 1e-6);
      __m256d ScaleRIP = _mm256_set1_pd(ScaleRI[P]);
      v_AR1 += v_LRp1 * ScaleRIP; v_AR2 += v_LRp2 * ScaleRIP;// AR += LRp * ScaleRI[P]
    }
  }

  /* return results */
  _mm256_store_pd(AR, v_AR1);
  _mm256_store_pd(&AR[4], v_AR2);

  if(DEBUG/* HERE >=2 */){
    for(int h = 0; h < 8; h++){
#if VECTOR2_DEBUG
      if(DEBUG && !(origUmax < 0 ? AR[h] - ARqh[h] > -fabs(ARqh[h]) * 1e-4 - 1e-9 
		    : fabs(ARqh[h] - AR[h]) < fabs(ARqh[h]) * 1e-4 + 1e-9)){
	#pragma omp critical
	{
	  printf("\nARq_mm256:m=%d,I=%d,K=%d,J=%d,U=%d..%d,Ub=%d,Yik[h]=%0.4f,Yt[h]=%0.4f,t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d:\n",m,I,K,J,origU,origUmax,origUb,Yik[h],Yt[h],t,M,N,IMIN,IMAX);
	  printf("\t h=%d: AR[h]= %0.8e, ARqh[h]= %0.8e (final U=%d..%d,Ub=%d)\n",h,AR[h],ARqh[h],finalUmin,finalUmax,finalUb);
	  for(int t = 0; t < 8; t++)
	    printf("\t t=%d: Yik= %0.6f, Yt= %0.6f: AR[t]= %0.8e, ARqh[h]= %0.8e\n", t, Yik[t], Yt[t], AR[t], ARqh[t]);
	  fflush(stdout);
	  if(origUmax < 0)
	    assert(AR[h] - ARqh[h] > -(fabs(ARqh[h]) * 1e-4 - 1e-9));
	  else
	    assert(fabs(ARqh[h] - AR[h]) < fabs(ARqh[h]) * 1e-4 + 1e-9);
	}
      }
#endif      
      if(DEBUG && !(isfinite(AR[h]) && AR[h] >= 0.0)){
	#pragma omp critical
	{
	  printf("ARq_mm256:m=%d,I=%d,K=%d,J=%d,t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d:AR[%d]= %0.8e\n",m,I,K,J,t,M,N,IMIN,IMAX,h,AR[h]);
	  assert(isfinite(AR[h]) && AR[h] >= 0.0);
	}
      }
    }
  }
}

#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT

template<int OUT_TYPE, int OUT_OLD>
static inline void ARq_mm512(double *AR,/* return values AR[0..7] */
			     int m, /* for debugging only */
			     int I,
			     int K,
			     int J,
			     int U, /* initialized to Umin */
			     int Ub,/* estimated peak of range U .. Umax (only used with ERFC_FIX if U <= Umax and can be outside the range U .. Umax) */
			     int Umax,/* If Umax < U, U range is unknown */
			     __m512d v_Yik1, __m512d v_Yik2, /* 8 New values for Yik */
			     __m512d v_Yt1, __m512d v_Yt2, /* 8 New values for Y[t] */
			     int t,/* t is either I or I-K */
			     double *  Xm,
			     double *  Y,// NOTE : ARq depends on Y[t] (see FAIKP), but should never use Y[t] directly
			     int M,
			     int N,
			     int rc,
			     double Scale,/* exp(-LogScaleR[I]) */
			     double *  ScaleRI,/* scaleRI[P] = exp(LogScaleR[P] - LogScaleR[I]) */
			     double * LogScaleR, /* for debugging only */
			     QB &  A,
			     int IMIN,
			     int IMAX,
			     int *  Kmax,
			     int *  Jmin,
			     int *  Jmax)
{
#if DEBUG>=2 || VECTOR2_DEBUG
  int origU = U, origUmax = Umax, origUb = Ub, finalUmin = -1, finalUmax = -1, finalUb = -1;
  double ARqh[16];
  double Yik[16];
  double Yt[16];
  _mm512_storeu_pd(Yik, v_Yik1);  _mm512_storeu_pd(&Yik[8], v_Yik2);
  _mm512_storeu_pd(Yt, v_Yt1); _mm512_storeu_pd(&Yt[8], v_Yt2);
  for(int h = 0; h < 16; h++)
    ARqh[h] = ARq<OUT_TYPE,OUT_OLD>(m,I,K,J,U,Ub,Umax,Yik[h],Yt[h],t,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax);
#endif  

  __m512d v_UR1 = _mm512_set1_pd(PoutlierEnd);// UR = PoutlierEnd
  __m512d v_UR2 = v_UR1;

  int JMIN = (EXTEND && I >= N+1 - DELTA_YU) ? (M+1 -DELTA_XU - EXTEND_MAX) : (M+1 - DELTA_XU);

  if(J >= JMIN){
    // NOTE : U range is U .. Umax, with Ub corresponding to Rij+1, the estimated peak value (can be outside the range U..Umax)

    if(U > Umax){/* need to determine U range : Use first and last Yik values to compute upper bound of range */
      double Yik16[VEC_SIZ] __attribute__((aligned(64)));
      _mm512_store_pd(Yik16, v_Yik1);
      _mm512_store_pd(&Yik16[8], v_Yik2);
      int Ulb = N+1, Uub = I, Umb = -1;/* U range will be computed Ulb .. Uub, with Umb corresponding to Rij+1, the estimated peak value (can lie outside the range Ulb .. Uub) */
      
      for(int T = 0; T <= 1; T++){
	double Yik = T ? Yik16[15] : Yik16[0];

	int Rijk = (IMIN <= I && I <= IMAX && 0 <= K && K <= Kmax[I] && Jmin[I] <= J && J <= Jmax[I]) ? A.Rij(I,K,J) : N+1;

	/* decrement Rijk until Y[Rijk] overlaps X */
	while(Rijk > I && Y[Rijk] - Yik > Xm[M+1]-Xm[J])
	  Rijk--;
	/* increment Rijk until Y[Rijk+1] no longer overlaps X */
	while(Rijk <= N && Y[Rijk+1] - Yik <= Xm[M+1]-Xm[J])
	  Rijk++;
	int Umid, Umin;
	Umid = Umin = Umax = min(N,Rijk)+1;

	MFLOAT URu = FR<0>(I,K,J,Umin,Umin-1,Yik,Xm[M+1]-Xm[J],Y[Umin]-Yik,Y[Umin]-Y[Umin-1],Xm,Y,M,N,END_FIX ? rc : 1);
	MFLOAT bUR = URu;

	for(U=Umin; --U > I;){
	  MFLOAT FRu = ERFC_FIX ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)	  
	    : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
	  if(FRu < URu * EPS)
	    break;
	  if(URANGE_FIX && FRu > bUR){
	    Umax = U;
	    URu = 0.0;
	    bUR = FRu;
	  }
	  Umin = U;
	  URu += FRu;
	}

	for(U = Umax; ++U <= N+1;){
	  MFLOAT FRu = (ERFC_FIX && U < Umid) ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)
	    : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
	  if(FRu < URu * EPS)
	    break;
	  Umax = U;
	  URu += FRu;
	}

	if(ISLINEAR && rc){
	  MFLOAT FRu = FRE(I,K,J,N+1,N,Xm[M+1]-Xm[J],Y[N+1]-Yik,M);
	  if(FRu >= URu * EPS){
	    Umax = N+1;
	    URu += FRu;
	  }
	}

	if(DEBUG) assert(Umin <= Umax);
	Ulb = min(Ulb, Umin);
	Uub = max(Uub, Umax);
	Umb = T ? (Umb + Umid)/2 : Umid;
      } // T = 0..1
      
      if(DEBUG) assert(Ulb <= Uub);
      U = Ulb;
      Ub = Umb;
      Umax = Uub;
#if DEBUG>=2 || VECTOR2_DEBUG
      finalUmin = U;
      finalUmax = Umax;
      finalUb = Ub;
#endif
    } // U > Umax

    __m512 Xj = _mm512_set1_ps((float)(Xm[M+1]-Xm[J]));
    __m512 v_URu = _mm512_setzero_ps();

    for(;U <= Umax;U++){
      // Yui = (float)(Y[U] - Yik)
      __m512d Yu = _mm512_set1_pd(Y[U]);
      __m512 v_Yui = DoubleToFloat_mm512(_mm512_sub_pd(Yu, v_Yik1),_mm512_sub_pd(Yu, v_Yik2));

      __m512 v_Yuv = _mm512_set1_ps((float)(Y[U] - Y[U-1])); // Yuv = Y[U] - Y[U-1]

      v_URu = _mm512_add_ps(v_URu, 
			    (ERFC_FIX && U < Ub) ? FR_mm512<1>(I,K,J,U,U-1, Xj, v_Yui, v_Yuv, M, N, END_FIX ? rc : 1)
			    : FR_mm512<0>(I,K,J,U,U-1, Xj, v_Yui, v_Yuv, M, N, END_FIX ? rc : 1));

#if 0 // HERE : remove this when Intel compiler bug is fixed
      if(rc && U==N+1)
	v_URu = _mm512_add_ps(v_URu, 
			      FRE_mm512(I,K,J,N+1,N, Xj, v_Yui, M));
#endif // HERE
    }

    /* convert v_URu to double precision */
#if USE_MIC
    __m512d v_URu1 = _mm512_cvtpslo_pd(v_URu);/* low order 8 floats to double */
    __m512d v_URu2 = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_URu, 0x2 | (0x3 << 2)));/* high order 8 floats to double */
#elif USE_AVX512
    __m512d v_URu1 = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_castsi512_ps(_mm512_castsi256_si512(_mm512_extracti64x4_epi64(_mm512_castps_si512(v_URu),0)))));/* low order 8 floats to double */
    __m512d v_URu2 = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_castsi512_ps(_mm512_castsi256_si512(_mm512_extracti64x4_epi64(_mm512_castps_si512(v_URu),1)))));/* high order 8 floats to double */
#endif
    __m512d v_PoutlierEndM1 = _mm512_set1_pd(PoutlierEndM1);    

    // UR += URu * PoutlierEndM1
    v_UR1 = _mm512_add_pd(v_UR1, _mm512_mul_pd(v_URu1, v_PoutlierEndM1));
    v_UR2 = _mm512_add_pd(v_UR2, _mm512_mul_pd(v_URu2, v_PoutlierEndM1));
  } /* J <= JMIN */

#if DEBUG>=2 || VECTOR2_DEBUG
  double UR16[16];
  _mm512_storeu_pd(UR16, v_UR1);
  _mm512_storeu_pd(&UR16[8], v_UR2);
  for(int t = 0; t < 16; t++)
    assert(isfinite(UR16[t]) && UR16[t] >= 0.0);
#endif

  // AR = UR * Scale
  __m512d v_Scale = _mm512_set1_pd(Scale);
  __m512d v_AR1 = _mm512_mul_pd(v_UR1, v_Scale);
  __m512d v_AR2 = _mm512_mul_pd(v_UR2, v_Scale);

#if DEBUG>=2 || VECTOR2_DEBUG
  double AR16[16];
  _mm512_storeu_pd(AR16, v_AR1);
  _mm512_storeu_pd(&AR16[8], v_AR2);
  for(int t = 0; t < 16; t++)
    assert(isfinite(AR16[t]) && AR16[t] >= 0.0);
  if(TRACE){
    printf("m=%d,I=%d,K=%d,J=%d,U=%d..%d,Ub=%d(final U=%d..%d,Ub=%d),h=%d: UR[h]= %0.8e, Scale = %0.8e, AR[h]= %0.8e\n",
	m,I,K,J,origU,origUmax,origUb,finalUmin,finalUmax,finalUb,HTRACE,UR16[HTRACE], Scale, AR16[HTRACE]);
    fflush(stdout);
  }
#endif    

  int Pmax = min(IMAX, I + deltaExtYRef);
  int Pmax2 = min(IMAX, I + DELTA_Y);
  int Qmax = min(M, J + deltaExtXRef);
  int Qmax2 = min(M, J + DELTA_X);

  __m512d v_Yi1, v_Yi2, v_Yimk1, v_Yimk2;
  if(t == I){
    v_Yi1 = v_Yt1; v_Yi2 = v_Yt2; // Yi = Yt
    if(K){    //    Yimk = K ? Y[I-K] : Yt;
      v_Yimk1 = _mm512_set1_pd(Y[I-K]); v_Yimk2 = v_Yimk1;
    } else {
      v_Yimk1 = v_Yt1; v_Yimk2 = v_Yt2;
    }
  } else {
    if(DEBUG/* HERE >=2 */) assert(K > 0 && t == I - K);
    v_Yi1 = _mm512_set1_pd(Y[I]); v_Yi2 = v_Yi1;// Yi = Y[I]
    v_Yimk1 = v_Yt1; v_Yimk2 = v_Yt2; // Yimk = Yt
  }

  __m512d v_PrBiasWTi1, v_PrBiasWTi2;/* value of PrBiasWT[I] may not be valid, since Y[t] is a new label, but differentials for PrBiasWT[t+1..N+1] remain valid */
  if(t==I){/* compute correct value of PrBiasWT[I] (relative to PrBiasWT[t+1..N+1]) */
    __m512d v_1p0 = _mm512_set1_pd(1.0);
    __m512d v_eps = _mm512_set1_pd(1e-300);

    __m512d ytp1 = _mm512_set1_pd(Y[t+1]);
    __m512d v_y1 = _mm512_sub_pd(ytp1, v_Yt1), v_y2 = _mm512_sub_pd(ytp1, v_Yt2);// y = Y[t+1] - Yt

    __m512d v_pr1 = Pr_mm512d(v_y1), v_pr2 = Pr_mm512d(v_y2);// pr = Pr(y)
    __m512d v_1Mpr1 = _mm512_sub_pd(v_1p0, v_pr1), v_1Mpr2 = _mm512_sub_pd(v_1p0, v_pr2);// 1.0 - pr
    // prbias = pr * log(max(pr,1e-300)) + (1.0 - pr) * log(max(1.0 - pr, 1e-300))
    __m512d v_prbias1 = _mm512_add_pd(_mm512_mul_pd(v_pr1, log_mm512d(_mm512_max_pd(v_pr1, v_eps))), _mm512_mul_pd(v_1Mpr1, log_mm512d(_mm512_max_pd(v_1Mpr1, v_eps))));
    __m512d v_prbias2 = _mm512_add_pd(_mm512_mul_pd(v_pr2, log_mm512d(_mm512_max_pd(v_pr2, v_eps))), _mm512_mul_pd(v_1Mpr2, log_mm512d(_mm512_max_pd(v_1Mpr2, v_eps))));

    __m512d v_PrBiasWTtp1 = _mm512_set1_pd(PrBiasWT[t+1]);
    __m512d v_PRbiasWT = _mm512_set1_pd(PRbiasWT);
    v_PrBiasWTi1 = _mm512_add_pd(v_PrBiasWTtp1, _mm512_mul_pd(v_PRbiasWT, v_prbias1));  v_PrBiasWTi2 = _mm512_add_pd(v_PrBiasWTtp1, _mm512_mul_pd(v_PRbiasWT, v_prbias2));// PrBiasWTi = PrBiasWT[t+1] + PrBiasWT * prbias
  } else
    v_PrBiasWTi1 = v_PrBiasWTi2 = _mm512_set1_pd(PrBiasWT[I]);

  __m512 v_resL = DoubleToFloat_mm512(_mm512_sub_pd(v_Yi1, v_Yimk1), _mm512_sub_pd(v_Yi2, v_Yimk2));// resL = (float)(Yi - Yimk)

  for(int P= I; ++P <= Pmax + KMAX && P <= IMAX;){
    __m512 v_LRp = _mm512_setzero_ps(); // LRp = 0.0

    int Tmax = min(P-I-1,Kmax[P]);
    for(int T = max(0,P-Pmax); T <= Tmax; T++){
      __m512 v_Ivar;
      __m512d Ypt = _mm512_set1_pd(Yc(Y,P,T));
      __m512 v_deltaY = DoubleToFloat_mm512(_mm512_sub_pd(Ypt, v_Yik1), _mm512_sub_pd(Ypt, v_Yik2));// deltaY = (float)(Yc(Y,P,T) - Yik)

      __m512 v_resR = _mm512_set1_ps((float)(Y[P] - Y[P-T])); // resR = Y[P] - Y[P-T]

      __m512d Ypmt = _mm512_set1_pd(Y[P-T]);
      __m512 v_Yint =  DoubleToFloat_mm512(_mm512_sub_pd(Ypmt, v_Yi1), _mm512_sub_pd(Ypmt, v_Yi2));// Yint = (float)(Y[P-T] - Yi)

      __m512 v_LRik = FAIKP_mm512(v_deltaY, v_resR, v_resL, v_Yint, v_Ivar);

      int QmaxP = (P - T <= Pmax2) ? min(Qmax,Jmax[P]) : min(Qmax2, Jmax[P]);
      MFLOAT2 *ARpt = &A.AR(P,T,0);

      __m512 LogFnPowPTI = _mm512_set1_ps(LogFnPow[P-T-I-1]);

      __m512d v_PrBiasWTpt = _mm512_set1_pd(PrBiasWT[P-T]);
      __m512 v_PrBiasWTpti = DoubleToFloat_mm512(_mm512_sub_pd(v_PrBiasWTpt, v_PrBiasWTi1), _mm512_sub_pd(v_PrBiasWTpt, v_PrBiasWTi2));// PrBiasWTpti = PrBiasWT[P-T] - PrBiasWTi

      if(DEBUG>=2){
	float PrBiasWTpti[16]; _mm512_storeu_ps(PrBiasWTpti, v_PrBiasWTpti);
	for(int k = 0; k < 16; k++){
	  if(DEBUG && !(PrBiasWTpti[k] >= 0.0)){
	    double PrBiasWTpt[8]; _mm512_storeu_pd(PrBiasWTpt, v_PrBiasWTpt);	    
	    double PrBiasWTi[16]; _mm512_storeu_pd(PrBiasWTi, v_PrBiasWTi1); _mm512_storeu_pd(&PrBiasWTi[8], v_PrBiasWTi2);	    
	    #pragma omp critical
	    {
	      printf("m=%d:I=%d,K=%d,J=%d,t=%d,P=%d,T=%d:k=%d:PrBiasWT[P-T]= %0.8e, PrBiasWTpt[k %% 8]= %0.8e, PrBiasWTi[k]= %0.8e, PrBiasWTpti[k]= %0.8e\n",
		     m,I,K,J,t,P,T,k,PrBiasWT[P-T],PrBiasWTpt[k % 8], PrBiasWTi[k], PrBiasWTpti[k]);
	      fflush(stdout);
	      assert(PrBiasWTpti[k] >= 0.0);
	    }
	  }
	}
      }

      // Critical Loop
      for(int Q= max(J+1,Jmin[P]); Q <= QmaxP; Q++){
        v_LRp = _mm512_add_ps(v_LRp, _mm512_mul_ps(_mm512_set1_ps(ARpt[Q]), 
		   FAJ_mm512<1,OUT_TYPE,OUT_OLD>((float)(Xm[Q]-Xm[J]), v_deltaY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPowPTI, v_PrBiasWTpti, v_LRik, v_Ivar)));

#if VECTOR2_DEBUG
	if(TRACE){
	  float LRp[16];
	  _mm512_storeu_ps(LRp, v_LRp);
	  _mm512_storeu_pd(AR, v_AR1);
	  _mm512_storeu_pd(&AR[8], v_AR2);
	  
	  printf("\t P=%d,T=%d,Q=%d,h=%d: LRp= %0.8e, ScaleRI[P]= %0.8e (AR= %0.8e)\n", P, T, Q, HTRACE, LRp[HTRACE], ScaleRI[P], LRp[HTRACE] * ScaleRI[P] + AR[HTRACE]);
	  fflush(stdout);
	}
#endif	
      }
    }

    /* convert v_LRp to double precision */
#if USE_MIC
    __m512d v_LRp1 = _mm512_cvtpslo_pd(v_LRp);/* low order 8 floats to double */
    __m512d v_LRp2 = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_LRp, 0x2 | (0x3 << 2)));/* high order 8 floats to double */
#elif USE_AVX512
    __m512d v_LRp1 = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_castsi512_ps(_mm512_castsi256_si512(_mm512_extracti64x4_epi64(_mm512_castps_si512(v_LRp),0)))));/* low order 8 floats to double */
    __m512d v_LRp2 = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_castsi512_ps(_mm512_castsi256_si512(_mm512_extracti64x4_epi64(_mm512_castps_si512(v_LRp),1)))));/* high order 8 floats to double */
#endif

    if(!isfinite(ScaleRI[P])){/* NEW144 : check if all LRp[] values are 0 */
      double LRp[16];
      _mm512_store_pd(LRp, v_LRp1);
      _mm512_store_pd(&LRp[8], v_LRp2);
      double maxLRp = 0.0;
      for(int h = 0; h < 16; h++)
	maxLRp = max(maxLRp, LRp[h]);
      if(maxLRp > 0.0){
	printf("ARq_mm256:m=%d,I=%d,K=%d,J=%d,P=%d(Pmax=%d,KMAX=%d),t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d:\n",m,I,K,J,P,Pmax,KMAX,t,M,N,IMIN,IMAX);
	printf("\t ScaleRI[P]= %0.17e, LogScaleR[P]= %0.17f, LogScaleR[I]= %0.17f,maxLRp= %0.17e\n",ScaleRI[P], LogScaleR[P], LogScaleR[I],maxLRp);
	fflush(stdout);
	assert(isfinite(ScaleRI[P]));
      }
    } else {
      // AR += LRp * ScaleRI[P]
      if(DEBUG>=2) assert(isfinite(ScaleRI[P]) && fabs(ScaleRI[P] - exp(LogScaleR[P]-LogScaleR[I])) < fabs(ScaleRI[P]) * 1e-5 + 1e-6);
      __m512d ScaleRIP = _mm512_set1_pd(ScaleRI[P]);
      v_AR1 = _mm512_add_pd(v_AR1, _mm512_mul_pd(v_LRp1, ScaleRIP));
      v_AR2 = _mm512_add_pd(v_AR2, _mm512_mul_pd(v_LRp2, ScaleRIP));
    }
  }

  /* return results */
  _mm512_store_pd(AR, v_AR1);
  _mm512_store_pd(&AR[8], v_AR2);

  if(DEBUG/* HERE >=2 */){
    for(int h = 0; h < 16; h++){
#if VECTOR2_DEBUG
      if(DEBUG && !(origUmax < 0 ? AR[h] - ARqh[h] > -fabs(ARqh[h]) * 1e-4 - 1e-9 
		    : fabs(ARqh[h] - AR[h]) < fabs(ARqh[h]) * 1e-4 + 1e-9)){
	#pragma omp critical
	{
	  printf("ARq_mm512(m=%d,I=%d,K=%d,J=%d,U=%d..%d,Ub=%d,Yik[h]=%0.4f,Yt[h]=%0.4f,t=%d,M=%d,N=%d,IMIN=%d,IMAX=%d:\n",m,I,K,J,origU,origUmax,origUb,Yik[h],Yt[h],t,M,N,IMIN,IMAX);
	  printf("\t h=%d: AR[h]= %0.8e, ARqh[h]= %0.8e (final U=%d..%d,Ub=%d)\n",h,AR[h],ARqh[h],finalUmin,finalUmax,finalUb);
	  fflush(stdout);
	  if(origUmax < 0)
	    assert(AR[h] - ARqh[h] > -(fabs(ARqh[h]) * 1e-4 - 1e-9));
	  else
	    assert(fabs(ARqh[h] - AR[h]) < fabs(ARqh[h]) * 1e-4 - 1e-9);
	}
      }
#endif      
      assert(isfinite(AR[h]) && AR[h] >= 0.0);
    }
  }
}

#endif // VECTOR_AVX || VECTOR_MIC

	/* Vectorize for IMIN <= S && S < IMAX1 : 
	   1. jmin = Jmin[S]; jmax = Jmax[S+1]
	   2. Compute Kt and Ks for each Ht (sequentially) : Kt will be monotonically decreasing, Ks monotonically increasing with t (but mostly 0)
	   3. Approximate ULmin,ULmax,URmin,URmax as bounds for current S and each K= 0..Kt, J= jmin..jmax value:
	       ULmin[K,J] = (J <= Jmax[S]) ? A.ULmin(S,max(0,K-1),J) : "re-compute since A.ULmin(S,*,J) is out of bounds"
	       ULmax[K,J] = (J >= Jmin[S+1]) ? A.ULmax(S+1,K,J) : "re-compute since A.ULmax(S+1,*,J) is out of bounds"
	       URmin[K,J] = (J <= Jmax[S]) ? A.URmin(S,max(0,K-1),J) : "re-compute since A.URmin(S,*,J) is out of bounds"
	       URmax[K,J] = (J >= Jmin[S+1]) ? A.URmax(S+1,K,J) : "re-compute since A.URmax(S+1,*,J) is out of bounds"
	   4. Approximate ULmin,ULmax,URmin,URmax as bounds for current S and each K= 1..Ks, J= Jmin[S+K]..Jmax[S+K]  value:
	       ULmin[K,J] = (K <= Kmax[S+K]) ? A.ULmin(S+K,K,J)  : "re-compute since A.ULmin(S+K,K,*) is out of bounds"
	       ULmax[K,J] = (K-1 <= Kmax[S+K]) ? A.ULmax(S+K,K-1,J) : "re-compute since A.ULmax(S+K,K-1,*) is out of bounds"
	       URmin[K,J] = (K <= Kmax[S+K]) ? A.URmin(S+K,K,J) : "re-compute since A.URmin(S+K,K,*) is out of bounds"
	       URmax[K,J] = (K-1 <= Kmax[S+K]) ? A.URmax(S+K,K-1,J) : "re-compute since A.URmax(S+K,K-1,*) is out of bounds"
	   5. Vectorize just the case for K = 0 (and possibly K == 1 if there are enough Kt or Ks values >= 1 to form a vector), which should cover most of the run time:
	      A. Vectorize across Ht values, calling vectorized version of AL() and AR().
	         In AL,AR, when computing Umin, use smallest value of Ht and when computing Umax, use largest value of Ht.
	   
	 */

template<int OUT_TYPE, int OUT_OLD>
static __attribute__ ((noinline)) void mprobeval_add(QB &A, int m, int M, int N, int lc, int rc,
				 int IMIN, int IMAX, int Amin, int Amax, int *Kmax, int *Jmin, int *Jmax,
				 double *Y, double *Xm, double *origY,
				 double *LogScaleL, double *LogScaleR, double *ScaleLI, double *ScaleRI, double LS,
				 int *KtVec, int *KsVec,
				 double **PMTH, double ***ALTH, double ***ARTH,
			         FL128 *LRadd, FL128 *LRadd2, // debug only
				 double **addlocI, int *addcntI, int *addcntIcum,
				 double **LRaddT)
{
  if(DEBUG) assert(SCORE_APPROX);
  if(DEBUG) assert(Y_APPROX);

  for(int S= IMIN; S < IMAX; S++){
    if(addcntI[S] <= 0)
      continue;

    //    FL128b ScaleS = exp(LogScaleL[min(max(S+1,IMIN-1),IMAX+1)]+LogScaleR[min(max(S,IMIN-1),IMAX+1)]-LS);
    FL128b ScaleS = exp(LogScaleL[S+1] + LogScaleR[S] - LS);

    if(DEBUG>=2) assert(isfinite(ScaleS) && ScaleS >= 0.0);
    if(DEBUG>=2) assert(Amin <= S && S <= Amax);
    if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && VERB){
      printf("Looking for updates to LRaddT[%d][%d]:S=%d,addcntIcum[S]=%d,addcntI[S]=%d (IMIN=%d,IMAX=%d)\n",
	     MDEBUG_TS,MDEBUG_TT,S,addcntIcum[S],addcntI[S],IMIN,IMAX);
      fflush(stdout);
    }
	
    int KtEnd[KMAX+1];/* KtEnd[k=0..KMAX] is last (largest) t value for which KtVec[t] >= k. KtEnd[k] == -1 if KtVec[t] is never >= k. */
    //    int KsStart[KMAX+1];/* KsStart[k=0..KMAX] is first (smallest) t value for which KsVec[t] >= k. KsStart[k] == -1 if KsVec[t] is always >= k. KsStart[k] == addcntI[S] if KsVec[t] is never >= k. */

    for(int t = 0; t < addcntI[S]; t++){
      double Ht = addlocI[S][t];
      if(DEBUG>=2) assert(Y[S] < Ht && Ht < Y[S+1]);
      if(DEBUG>=2 && t > 0 && !(Ht >= addlocI[S][t-1])){
	#pragma omp critical
	{
	  printf("m=%d,S=%d,N=%d: addcntI[S]=%d, t=%d : addlocI[S][t-1]= %0.8f, addlocI[S][t]= %0.8f, skip_dist= %0.8f\n",
		 m,S,N,addcntI[S],t, addlocI[S][t-1], addlocI[S][t], skip_dist);
	  fflush(stdout);
	  assert(Ht > addlocI[S][t-1]);
	}
      }

      /* compute largest K so that Ht and Y[S+1-Kt] may fail to resolve */
      /* determine Kt based on total distance Ht - Y[S+1-Kt] < maxdist (if KMAX_APPROX)*/
      int Kt = 0;
      while (S+1-Kt > 1){
	if((KMAX_APPROX || (!KMAX_FIX && Kt==0)) ? (Ht - Y[S-Kt] >= maxdist) : ((!KMAX_FIX && Y[S+1-Kt] - Y[S-Kt] >= maxdist) || Ht - Y[S-Kt] >= AlignResMax))
	  break;
	Kt++;
      }
      if(Kt > KMAX)
	Kt = KMAX;

      /* compute largest K so that Ht and Y[S+K] may fail to resolve */
      /* determine Ks based on total distance Y[S+Ks] - Ht < maxdist (if KMAX_APPROX) */
      int Ks = 0;
      while (S + Ks < N) {
	if((KMAX_APPROX || (!KMAX_FIX && Ks==0)) ? (Y[S+Ks+1] - Ht >= maxdist) : ((!KMAX_FIX && Y[S+Ks+1] - Y[S+Ks] >= maxdist) || Y[S+Ks+1] - Ht >= AlignResMax))
	  break;
	Ks++;
      }
	  
      if(Ks > KMAX)
	Ks = KMAX;
      if(S+Ks > IMAX)
	Ks = IMAX-S;

      if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT){
	printf("m=%d:S=%d,t=%d(T=%d):Kt=%d,Ks=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+Ks]=%0.3f\n",m,S,t,addcntIcum[S]+t,Kt,Ks,Ht,Y[S],Y[S+Ks]);
	fflush(stdout);
      }

      KtVec[t] = Kt;
      KsVec[t] = Ks;
      if(DEBUG && t > 0){
	assert(Kt <= KtVec[t-1]);
	assert(Ks >= KsVec[t-1]);
      }
      if(!t){
	for(int k = Kt + 1; k <= KMAX; k++)
	  KtEnd[k] = -1;

	/*	for(int k = 0; k < Ks; k++)
	  KsStart[k] = -1;
	  KsStart[Ks] = t;*/
      } else {
	if(Kt < KtVec[t-1])
	  for(int k = KtVec[t-1]; k > Kt; k--)
	    KtEnd[k] = t-1;
	/*	if(Ks > KsVec[t-1])
	  for(int k = KsVec[t-1] + 1; k <= Ks; k++)
	  KsStart[k] = t;*/
      }
    }// t = 0 .. addcntI[S] - 1 
    int t = addcntI[S] - 1;
    for(int k = KtVec[t]; k >= 0; k--)
      KtEnd[k] = t;
    /*    for(int k = KsVec[t] + 1; k <= KMAX; k++)
	  KsStart[k] = t + 1;*/

    /* compute Jmin[] and Jmax[] for new site T */
    /* NOTE : this will occasionally be a larger range than computed by setlimit() after site addition */
    int jmin = Jmin[S];
    int jmax = Jmax[S+1];
    if(DEBUG && !(jmin <= jmax)){
      #pragma omp critical
      {
	printf("m=%d,M=%d,N=%d:IMIN=%d,IMAX=%d,S=%d:Jmin[S,S+1]=%d,%d,Jmax[S,S+1]=%d,%d, jmin=%d,jmax=%d\n",
	       m,M,N,IMIN,IMAX,S,Jmin[S],Jmin[S+1],Jmax[S+1],Jmax[S+1],jmin,jmax);
	fflush(stdout);
	assert(jmin <= jmax);
      }
    }
    if(DEBUG>=2) assert(0 < jmin && jmax <= M);

    for(int T = 0; T < addcntI[S]; T += VEC_SIZ){

      int I = S+1;
      /* precompute ScaleLI[P= I-1 .. max(IMIN, I - KtVec[T] - deltaExtYRef)]  = exp(LogScaleL[P] - LogScaleL[I]) */
      int Pmin = max(IMIN, I - KtVec[0] - deltaExtYRef);
      double LogScaleLI = LogScaleL[I];
      for(int P = Pmin; P < I; P++){
	if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
	ScaleLI[P] = exp(LogScaleL[P] - LogScaleLI);
	if(DEBUG>=2) assert(isfinite(ScaleLI[P]));
      }

      /* precompute ScaleRI[P= S+1 .. min(IMAX, S + KMAX + deltaExtYRef)] = exp(LogScaleR[P] - LogScaleR[S]) */
      int Pmax = min(IMAX, S + KMAX /* WAS KtVec[T] */ + deltaExtYRef);
      double LogScaleRI = LogScaleR[S];
      for(int P = max(IMIN,S); P <= Pmax; P++ ){
	if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
	ScaleRI[P] = exp(LogScaleR[P] - LogScaleRI);
	if(DEBUG>=2) assert(isfinite(ScaleRI[P]));		
      }

      /* compute Line 5 of Equation 10:
	 ALT[K= 0..Kt][J=jmin..jmax][h = 0 .. VEC_SIZ -1] == A[S+1][K][J].AL 
	 PMT[K][h] == PM(S+1,K)
	 with Y[S+1] = Ht, t = T + h */

      for(int K = 0; K <= KtVec[T]; K++){
	int vlen = min(VEC_SIZ, min(KtEnd[K]+1,addcntI[S]) - T);
	int h = 0;

#if VECTOR_AVX && USE_SSE && USE_AVX && !USE_AVX512 && USE_MFLOAT
	if(vlen >= VEC_MSIZ2){	/* Use AVX256 intrinsincs (manual optimization) */
	  /* NOTE : addlocI[S][] is aligned and padded, if needed, to a multiple of 8 entries : padded entries repeat the last actual value.
	     Results PMTH[K][], ALTH[K][J][], ARTH[K][J][] are also aligned and padded to multiples of 8 entries */
	  int t = T;
	  __m256d v_Ht1 = _mm256_load_pd(&addlocI[S][t]);
	  __m256d v_Ht2 = _mm256_load_pd(&addlocI[S][t+4]);

	  int I = S+1;
	  double Scale = exp(-LogScaleL[I]);

	  __m256d v_Yimk1, v_Yimk2;
	  if(K){
	    v_Yimk1 = _mm256_set1_pd(Y[I-K]); v_Yimk2 = v_Yimk1;
	  } else {
	    v_Yimk1 = v_Ht1; v_Yimk2 = v_Ht2;
	  }
	  __m256d v_Yik1 = v256d_0p5 * (v_Ht1 + v_Yimk1); __m256d v_Yik2 = v256d_0p5 * (v_Ht2 + v_Yimk2);// Yik = 0.5 * (Ht + Y[I-K])

	  for(int H = h; H < vlen; H++){/* hard to vectorize PM() : not worth it */
	    int t = T + H;
	    double Ht = addlocI[S][t];
	    Y[I] = Ht;
	    PMTH[K][H] = PM(0,I,K,Y);
	  }

	  Y[I] = origY[I];

	  int JMAX = (EXTEND && I-K <= DELTA_YU) ? (DELTA_XU + EXTEND_MAX) : DELTA_XU;

	  for(int J = jmin; J <= jmax; J++){
	    if(Jmin[S+1] <= J && J <= Jmax[S] && K-1 <= Kmax[S] && K <= Kmax[S+1]){
	      int Umin = I, Umax = 0, Ub = -1;
	      if(J <= JMAX){/* otherwise Umin,Umax are ignored by ALq() */
		if(DEBUG>=2) assert(IMIN <= S+1 && S+1 <= IMAX && 0 <= K && K <= Kmax[S+1] && Jmin[S+1] <= J && J <= Jmax[S+1]);
		int ULmin = A.ULmin(S+1,K,J), ULmax = A.ULmax(S+1,K,J);
		int ULb = max(1, A.Lij(S+1,K,J)) - 1;
		if(DEBUG>=2) assert(ULmin <= ULmax);
		if(ULmin <= ULmax){
		  if(DEBUG>=2) assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		  Umin = min(Umin, ULmin);
		  Umax = max(Umax, ULmax);
		  Ub = (Ub < 0) ? ULb : (Ub + ULb)/2;
		}

		if(DEBUG>=2) assert(IMIN <= S && S <= IMAX && 0 <= max(0,K-1) && max(0,K-1) <= Kmax[S] && Jmin[S] <= J && J <= Jmax[S]);
		ULmin = A.ULmin(S,max(0,K-1),J); ULmax = A.ULmax(S,max(0,K-1),J);
		ULb = max(1, A.Lij(S,max(0,K-1),J)) - 1;
		if(DEBUG>=2) assert(ULmin <= ULmax);
		if(ULmin <= ULmax){
		  if(DEBUG>=2) assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		  Umin = min(Umin, ULmin);
		  Umax = max(Umax, ULmax);
		  Ub = (Ub < 0) ? ULb : (Ub + ULb)/2;
		}

		if(DEBUG>=2) assert(Umin <= Umax);
		if(Umin > Umax)
		  Ub = Umin = Umax = I-K-1;
	      }

	      ALq_mm256<OUT_TYPE,OUT_OLD>(ALTH[K][J], m,I,K,J, Umin, Ub, Umax, v_Yik1, v_Yik2, v_Ht1, v_Ht2, I, Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	    } else
	      ALq_mm256<OUT_TYPE,OUT_OLD>(ALTH[K][J], m,I,K,J, 0, 0, -2, v_Yik1, v_Yik2, v_Ht1, v_Ht2, I, Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	  }
	  /* ALTH[K][J][h] does not include PMik : matches AL[S+1][K][J] on line 5 Equation 10 */

	  /* compute ARTH[K=0..Kt][J=jmin..jmax][h] == A[S][K][J].AR with Y[r-1] = origY[r] (for r=S+1-K ..S), Y[S] = Ht */
	  I = S;
	  Scale = exp(-LogScaleR[I]);

	  /* modify Y[] */
	  for(int r = S+1-K; r <= S; r++)
	    Y[r-1] = origY[r];

	  int JMIN = (EXTEND && I >= N+1 - DELTA_YU) ? (M+1 -DELTA_XU - EXTEND_MAX) : (M+1 -DELTA_XU);

	  for(int J = jmin; J <= jmax; J++){
	    if(Jmin[S+1] <= J && J <= Jmax[S] && K-1 <= Kmax[S] && K <= Kmax[S+1]){
	      int Umin = N+1, Umax = I, Ub = -1;
	      if(J >= JMIN){/* otherwise Umin,Umax are ignored by ARq() */
		if(DEBUG>=2) assert(IMIN <= S+1 && S+1 <= IMAX && 0 <= K && K <= Kmax[S+1] && Jmin[S+1] <= J && J <= Jmax[S+1]);
		int URmin = A.URmin(S+1,K,J), URmax = A.URmax(S+1,K,J);
		int URb = min(N, A.Rij(S+1,K,J)) + 1;
		if(DEBUG>=2) assert(URmin <= URmax);
		if(URmin <= URmax){
		  if(DEBUG>=2) assert(I < URmin && URmin <= URmax && URmax <= N+1);
		  Umin = min(Umin, URmin);
		  Umax = max(Umax, URmax);
		  Ub = (Ub < 0) ? URb : (Ub + URb)/2;
		}

		if(DEBUG>=2) assert(IMIN <= S && S <= IMAX && 0 <= max(0,K-1) && max(0,K-1) <= Kmax[S] && Jmin[S] <= J && J <= Jmax[S]);
		URmin = A.URmin(S,max(0,K-1),J); URmax = A.URmax(S,max(0,K-1),J);
		URb = min(N, A.Rij(S,max(0,K-1),J)) + 1;
		if(DEBUG>=2) assert(URmin <= URmax);
		if(URmin <= URmax){
		  if(DEBUG>=2) assert(I < URmin && URmin <= URmax && URmax <= N+1);
		  Umin = min(Umin, URmin);
		  Umax = max(Umax, URmax);
		  Ub = (Ub < 0) ? URb : (Ub + URb)/2;
		}

		if(DEBUG>=2) assert(Umin <= Umax);
		if(Umin > Umax)
		  Ub = Umin = Umax = I+1;
	      }
		  
	      ARq_mm256<OUT_TYPE,OUT_OLD>(ARTH[K][J], m,I,K,J,Umin,Ub,Umax,v_Yik1,v_Yik2,v_Ht1,v_Ht2, I,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	    } else
	      ARq_mm256<OUT_TYPE,OUT_OLD>(ARTH[K][J], m,I,K,J,0,0,-2,v_Yik1,v_Yik2,v_Ht1,v_Ht2 ,I,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	  }
	  /* No PMik : matches AR[S][K][J] on line 5 Equation 10 */

	  /* restore Y[S-K .. S-1] */
	  for(int r = S-K; r < S; r++)
	    Y[r] = origY[r];

	  h += VEC_SIZ;
	}
#elif VECTOR_MIC && (USE_MIC || USE_AVX512) && USE_MFLOAT
	if(DEBUG) assert(lc==0 && rc==0);/* HERE : remove this when Intel Compiler bug is fixed */

	if(vlen >= VEC_MSIZ2){	/* Use MIC intrinsincs (manual optimization) */
	  /* NOTE : addlocI[S][] is aligned and padded, if needed, to a multiple of 8 entries : padded entries repeat the last actual value.
	     Results PMTH[K][], ALTH[K][J][], ARTH[K][J][] are also aligned and padded to multiples of 8 entries */
	  int t = T;
	  __m512d v_Ht1 = _mm512_load_pd(&addlocI[S][t]);
	  __m512d v_Ht2 = _mm512_load_pd(&addlocI[S][t+8]);

	  int I = S+1;
	  double Scale = exp(-LogScaleL[I]);

	  __m512d v_Yimk1, v_Yimk2;
	  if(K){
	    v_Yimk1 = _mm512_set1_pd(Y[I-K]); v_Yimk2 = v_Yimk1;
	  } else {
	    v_Yimk1 = v_Ht1; v_Yimk2 = v_Ht2;
	  }

	  // Yik = 0.5 * (Ht + Y[I-K])
	  __m512d v_Yik1 = _mm512_mul_pd(v512d_0p5, _mm512_add_pd(v_Ht1, v_Yimk1));
	  __m512d v_Yik2 = _mm512_mul_pd(v512d_0p5, _mm512_add_pd(v_Ht2, v_Yimk2));

	  for(int H = h; H < vlen; H++){/* hard to vectorize PM() : not worth it */
	    int t = T + H;
	    double Ht = addlocI[S][t];
	    Y[I] = Ht;
	    PMTH[K][H] = PM(0,I,K,Y);
	  }

	  Y[I] = origY[I];

	  int JMAX = (EXTEND && I-K <= DELTA_YU) ? (DELTA_XU + EXTEND_MAX) : DELTA_XU;

	  for(int J = jmin; J <= jmax; J++){
	    if(Jmin[S+1] <= J && J <= Jmax[S] && K-1 <= Kmax[S] && K <= Kmax[S+1]){
	      int Umin = I, Umax = 0, Ub = -1;
	      if(J <= JMAX){/* otherwise Umin,Umax are ignored by ALq() */
		if(DEBUG>=2) assert(IMIN <= S+1 && S+1 <= IMAX && 0 <= K && K <= Kmax[S+1] && Jmin[S+1] <= J && J <= Jmax[S+1]);
		int ULmin = A.ULmin(S+1,K,J), ULmax = A.ULmax(S+1,K,J);
		int ULb = max(1, A.Lij(S+1,K,J)) - 1;
		if(DEBUG>=2) assert(ULmin <= ULmax);
		if(ULmin <= ULmax){
		  if(DEBUG>=2) assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		  Umin = min(Umin, ULmin);
		  Umax = max(Umax, ULmax);
		  Ub = (Ub < 0) ? ULb : (Ub + ULb)/2;
		}

		if(DEBUG>=2) assert(IMIN <= S && S <= IMAX && 0 <= max(0,K-1) && max(0,K-1) <= Kmax[S] && Jmin[S] <= J && J <= Jmax[S]);
		ULmin = A.ULmin(S,max(0,K-1),J); ULmax = A.ULmax(S,max(0,K-1),J);
		ULb = max(1, A.Lij(S,max(0,K-1),J)) - 1;
		if(DEBUG>=2) assert(ULmin <= ULmax);
		if(ULmin <= ULmax){
		  if(DEBUG>=2) assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		  Umin = min(Umin, ULmin);
		  Umax = max(Umax, ULmax);
		  Ub = (Ub < 0) ? ULb : (Ub + ULb)/2;
		}
		if(DEBUG>=2) assert(Umin <= Umax);
		if(Umin > Umax)
		  Ub = Umin = Umax = I-K-1;
	      }

	      ALq_mm512<OUT_TYPE,OUT_OLD>(ALTH[K][J], m,I,K,J, Umin, Ub, Umax, v_Yik1, v_Yik2, v_Ht1, v_Ht2, I, Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	    } else
	      ALq_mm512<OUT_TYPE,OUT_OLD>(ALTH[K][J], m,I,K,J, 0, 0, -2, v_Yik1, v_Yik2, v_Ht1, v_Ht2, I, Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	  }
	  /* ALTH[K][J][h] does not include PMik : matches AL[S+1][K][J] on line 5 Equation 10 */

	  /* compute ARTH[K=0..Kt][J=jmin..jmax][h] == A[S][K][J].AR 
	     with Y[r-1] = origY[r] (for r=S+1-K ..S), Y[S] = Ht, t = T+h */
	  I = S;
	  Scale = exp(-LogScaleR[I]);

	  /* modify Y[] */
	  for(int r = S+1-K; r <= S; r++)
	    Y[r-1] = origY[r];

	  int JMIN = (EXTEND && I >= N+1 - DELTA_YU) ? (M+1 -DELTA_XU - EXTEND_MAX) : (M+1 - DELTA_XU);

	  for(int J = jmin; J <= jmax; J++){
	    if(Jmin[S+1] <= J && J <= Jmax[S] && K-1 <= Kmax[S] && K <= Kmax[S+1]){
	      int Umin = N+1, Umax = I, Ub = -1;
	      if(J >= JMIN){/* otherwise Umin,Umax are ignored by ARq() */
		if(DEBUG>=2) assert(IMIN <= S+1 && S+1 <= IMAX && 0 <= K && K <= Kmax[S+1] && Jmin[S+1] <= J && J <= Jmax[S+1]);
		int URmin = A.URmin(S+1,K,J), URmax = A.URmax(S+1,K,J);
		int URb = min(N, A.Rij(S+1,K,J)) + 1;
		if(DEBUG>=2) assert(URmin <= URmax);
		if(URmin <= URmax){
		  if(DEBUG>=2) assert(I < URmin && URmin <= URmax && URmax <= N+1);
		  Umin = min(Umin, URmin);
		  Umax = max(Umax, URmax);
		  Ub = (Ub < 0) ? URb : (Ub + URb)/2;
		}

		if(DEBUG>=2) assert(IMIN <= S && S <= IMAX && 0 <= max(0,K-1) && max(0,K-1) <= Kmax[S] && Jmin[S] <= J && J <= Jmax[S]);
		URmin = A.URmin(S,max(0,K-1),J); URmax = A.URmax(S,max(0,K-1),J);
		URb = min(N, A.Rij(S,max(0,K-1),J)) + 1;
		if(DEBUG>=2) assert(URmin <= URmax);
		if(URmin <= URmax){
		  if(DEBUG>=2) assert(I < URmin && URmin <= URmax && URmax <= N+1);
		  Umin = min(Umin, URmin);
		  Umax = max(Umax, URmax);
		  Ub = (Ub < 0) ? URb : (Ub + URb)/2;
		}
		if(DEBUG>=2) assert(Umin <= Umax);
		if(Umin > Umax)
		  Ub = Umin = Umax = I+1;
	      }

	      ARq_mm512<OUT_TYPE, OUT_OLD>(ARTH[K][J], m,I,K,J,Umin,Ub,Umax,v_Yik1,v_Yik2,v_Ht1,v_Ht2, I,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	    } else
	      ARq_mm512<OUT_TYPE, OUT_OLD>(ARTH[K][J], m,I,K,J,0,0,-2,v_Yik1,v_Yik2,v_Ht1,v_Ht2 ,I,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	  }
	  /* No PMik : matches AR[S][K][J] on line 5 Equation 10 */

	  /* restore Y[S-K .. S-1] */
	  for(int r = S-K; r < S; r++)
	    Y[r] = origY[r];

	  h += VEC_SIZ;
	}
#endif // VECTOR_MIC

	for(; h < vlen; h++){
	  int t = T + h;
	  double Ht = addlocI[S][t];
	  int Kt = KtVec[t];

	  int I = S+1;
	  if(DEBUG>=2) assert(I >= Kt+1);

	  /* NOTE: AL(),AR() will NOT use PRtab[].Pr */

	  if(DEBUG>=2) assert(I <= N+1);
	  if(DEBUG>=2) assert(!(I > IMAX+1 || I < IMIN-1));
	  double Scale = exp(-LogScaleL[I]);
	  if(DEBUG>=2) assert(isfinite(Scale));
	  if(DEBUG>=2) assert(S-K <= IMAX);

	  Y[I] = Ht;
	  double Yik = 0.5*(Ht + Y[I-K]);
	  PMTH[K][h] = PM(0,S+1,K,Y);
	  Y[I] = origY[I];

	  if(DEBUG>=2) assert(I-K-1 <= IMAX);
		
	  int JMAX = (EXTEND && I-K <= DELTA_YU) ? (DELTA_XU + EXTEND_MAX) : DELTA_XU;

	  for(int J = jmin; J <= jmax; J++){
	    if(((Jmin[S+1] <= J && J <= Jmax[S] && K-1 <= Kmax[S] && K <= Kmax[S+1]) /*NEW*/|| J > JMAX)){
	      int Umin = I, Umax = 0, Ub = -1;

	      if(J <= JMAX){/* otherwise Umin,Umax are ignored by ALq() */
		if(DEBUG>=2) assert(IMIN <= S+1 && S+1 <= IMAX && 0 <= K && K <= Kmax[S+1] && Jmin[S+1] <= J && J <= Jmax[S+1]);
		int ULmin = A.ULmin(S+1,K,J), ULmax = A.ULmax(S+1,K,J);
		int ULb = max(1, A.Lij(S+1,K,J)) - 1;
		if(DEBUG>=2) assert(ULmin <= ULmax);
		if(ULmin <= ULmax){
		  if(DEBUG>=2) assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		  Umin = min(Umin, ULmin);
		  Umax = max(Umax, ULmax);
		  Ub = (Ub < 0) ? ULb : (Ub + ULb)/2;
		}

		if(DEBUG>=2) assert(IMIN <= S && S <= IMAX && 0 <= max(0,K-1) && max(0,K-1) <= Kmax[S] && Jmin[S] <= J && J <= Jmax[S]);
		ULmin = A.ULmin(S,max(0,K-1),J); ULmax = A.ULmax(S,max(0,K-1),J);
		ULb = max(1, A.Lij(S,max(0,K-1),J)) - 1;
		if(DEBUG>=2) assert(ULmin <= ULmax);
		if(ULmin <= ULmax){
		  if(DEBUG>=2) assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		  Umin = min(Umin, ULmin);
		  Umax = max(Umax, ULmax);
		  Ub = (Ub < 0) ? ULb : (Ub + ULb)/2;
		}

		if(DEBUG>=2) assert(Umin <= Umax);
		if(Umin > Umax)
		  Ub = Umin = Umax = I-K-1;
	      }

	      ALTH[K][J][h] = ALq<OUT_TYPE,OUT_OLD>(m,I,K,J, Umin, Ub, Umax, Yik, Ht, I, Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	    } else
	      ALTH[K][J][h] = ALq<OUT_TYPE,OUT_OLD>(m,I,K,J, 0, 0, -2, Yik, Ht, I, Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	  }
	  /* ALTH[K][J][h] does not include PMik : matches AL[S+1][K][J] on line 5 Equation 10 */

	  //	  Y[I] = origY[I];

	  /* compute ARTH[K=0..Kt][J=jmin..jmax][h] == A[S][K][J].AR with Y[r-1] = origY[r] (for r=S+1-K ..S), Y[S] = Ht */
	  if(DEBUG>=2) assert(S-Kt >= 0);
	  I = S;

	  if(DEBUG>=2) assert(0 <= I  && I <= N+1);
	  if(DEBUG>=2) assert(!(I < IMIN-1 || I > IMAX+1));
	  Scale = exp(-LogScaleR[I]);
	  if(DEBUG>=2) assert(isfinite(Scale));
	  if(DEBUG>=2) assert(S-K <= IMAX);
	  if(DEBUG>=2) assert(IMIN <= I+1);

	  /* modify Y[] */
	  Y[S] = Ht;
	  for(int r = S+1-K; r <= S; r++)
	    Y[r-1] = origY[r];

	  Yik = Yc(Y,I,K);
	  Y[S] = origY[S];

	  int JMIN = (EXTEND && I >= N+1-DELTA_YU) ? (M+1 - DELTA_XU - EXTEND_MAX) : (M+1 -DELTA_XU);

	  for(int J = jmin; J <= jmax; J++){
	    if(((Jmin[S+1] <= J && J <= Jmax[S] && K-1 <= Kmax[S] && K <= Kmax[S+1]) /*NEW*/|| J < JMIN)){
	      int Umin = N+1, Umax = I, Ub = -1;
	      if(J >= JMIN){/* otherwise Umin,Umax are ignored by ARq() */
		if(DEBUG>=2) assert(IMIN <= S+1 && S+1 <= IMAX && 0 <= K && K <= Kmax[S+1] && Jmin[S+1] <= J && J <= Jmax[S+1]);
		int URmin = A.URmin(S+1,K,J), URmax = A.URmax(S+1,K,J);
		int URb = min(N, A.Rij(S+1,K,J)) + 1;
		if(DEBUG>=2) assert(URmin <= URmax);
		if(URmin <= URmax){
		  if(DEBUG>=2 && !(I < URmin && URmin <= URmax && URmax <= N+1)){
		    printf("m=%d,S=%d,I=%d,K=%d,J=%d,N=%d: A.URmin(S+1,K,J)= %d, A.URmax(S+1,K,J)= %d, Kmax[S+1]= %d, Jmin[S+1]= %d, Jmax[S,S+1]= %d,%d\n",
			   m,S,I,K,J, N, A.URmin(S+1,K,J), A.URmax(S+1,K,J), Kmax[S+1], Jmin[S+1], Jmax[S], Jmax[S+1]);
		    fflush(stdout);
		    assert(I < URmin && URmin <= URmax && URmax <= N+1);
		  }
		  Umin = min(Umin, URmin);
		  Umax = max(Umax, URmax);
		  Ub = (Ub < 0) ? URb : (Ub + URb)/2;
		}

		if(DEBUG>=2) assert(IMIN <= S && S <= IMAX && 0 <= max(0,K-1) && max(0,K-1) <= Kmax[S] && Jmin[S] <= J && J <= Jmax[S]);
		URmin = A.URmin(S,max(0,K-1),J); URmax = A.URmax(S,max(0,K-1),J);
		URb = min(N, A.Rij(S,max(0,K-1),J)) + 1;
		if(DEBUG>=2) assert(URmin <= URmax);
		if(URmin <= URmax){
		  if(DEBUG>=2) assert(I < URmin && URmin <= URmax && URmax <= N+1);
		  Umin = min(Umin, URmin);
		  Umax = max(Umax, URmax);
		  Ub = (Ub < 0) ? URb : (Ub + URb)/2;
		}
		if(DEBUG>=2) assert(Umin <= Umax);
		if(Umin > Umax)
		  Ub = Umin = Umax = I+1;
	      }
		  
	      ARTH[K][J][h] = ARq<OUT_TYPE,OUT_OLD>(m,I,K,J,Umin,Ub,Umax,Yik,Ht,I,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	    } else
	      ARTH[K][J][h] = ARq<OUT_TYPE,OUT_OLD>(m,I,K,J,0,0,-2,Yik,Ht,I,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	  }
	  /* No PMik : matches AR[S][K][J] on line 5 Equation 10 */

	  /* restore Y[S-K .. S-1] */
	  for(int r = S-K; r < S; r++)
	    Y[r] = origY[r];
	  //	  Y[S] = origY[S];
	} /* h = 0 .. VEC_SIZ - 1 */
      }/* K = 0 .. KtVec[T] */

      for(int h = 0; h < VEC_SIZ; h++){
	int t = T + h;
	if(t >= addcntI[S])
	  break;
	int Kt = KtVec[t];

	/* add to LRaddT[S][t] : sum(J=jmin..jmax,K=0..Kt) ALT[K][J][h] * PMT[K][h] * ART[K][J][h] */
	for(int K = 0; K <= Kt; K++){
	  if(DEBUG>=2) assert(S-K <= IMAX);

	  double LR = 0.0;
	  for(int J = jmin; J <= jmax; J++)
	    LR += ALTH[K][J][h] * ARTH[K][J][h];

	  double LR2 = LR * PMTH[K][h] * ScaleS; // ScaleS == exp(LogScaleL[S+1] + LogScaleR[S] - LS)
	      
	  if(DEBUG>=2) assert(isfinite(LR2) && LR2>=0.0);

	  LRaddT[S][t] += LR2;

	  if(DEBUG>=2) assert(isfinite(LRaddT[S][t]) && LRaddT[S][t] >= 0.0);
	  if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT && LR != 0.0){
	    FL128 sum = LRadd[0],sum2 = LRadd2[0];
	    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
	      sum += LRadd[i];
	      sum2 += LRadd2[i];
	    }
	    FL128 sum3 = sum;
	    if(Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
	      sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
	      printf("m=%d:S=%d,t=%d(T=%d):I=%d,K=%d,J=%d..%d:sum(ALT[K][J]*PMT[K]*ART[K][J])=%0.10e(ScaleS=%0.10e):LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
		     m,S,t,addcntIcum[S]+t,S+1,K,jmin,jmax,(double)LR,(double)ScaleS,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
	      printf("\t LogScaleL[S+1]= %0.10f, LogScaleR[S]= %0.10f, LS= %0.10f\n",LogScaleL[S+1], LogScaleR[S], LS);

	      double LR3 = 0.0;
	      for(int J = jmin; J <= jmax; J++){
		LR3 += ALTH[K][J][h] * ARTH[K][J][h] * PMTH[K][h] * ScaleS;
		printf("    J=%d:ALT[K][J]*PM=%0.10e,ART[K][J]*PM=%0.10e,PMT[K]=%0.10f:sum(ALT*ART*PMT)=%0.10e(delta=%0.10e)\n",
		       J,(double)(ALTH[K][J][h]*PMTH[K][h]), (double)(ARTH[K][J][h]*PMTH[K][h]), (double)PMTH[K][h], (double)LR3, (double)(ALTH[K][J][h] * ARTH[K][J][h] * PMTH[K][h] * ScaleS));
	      }
	    } else
	      printf("m=%d:S=%d,t=%d(T=%d):K=%d,J=%d..%d:sum(ALT[K][J]*PMT[K]*ART[K][J])=%0.10e(ScaleS=%0.10e):LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
		     m,S,t,addcntIcum[S]+t,K,jmin,jmax,(double)LR,(double)ScaleS,(double)sum,(double)sum2,(double)(sum-sum2));
	    fflush(stdout);
	  }
	}
      }

      /* Now work on line 6 Equation 10 */
      // HERE HERE : not yet vectorized : accounts for only 10-20% of run time compared to code above, since K >= 1

      for(int h = 0; h < VEC_SIZ; h++){
	int t = T + h;
	if(t >= addcntI[S])
	  break;
	double Ht = addlocI[S][t];
	int Ks = KsVec[t];

	if(Ks <= 0)
	  continue;/* Line 6 of Equation 10 vanishes for current S, Ht */
	  
	int I = S+1;
	Y[I] = Ht;

	/* compute ALTH[K=1..Ks][J= jmin..Jmax[S+K]][h] == A[S+K+1][K][J].AL with Y[r+1] = origY[r], r = S+K ... S+1, Y[S+1] = Ht (line 6 of Equation 10) */
	/* NOTE: If _OPENMP, ALq(),ARq() will NOT use PRtab[].Pr */
	for(int K = 1; K <= Ks; K++){
	  I = S+K+1;
	  for(int r = S+K; r > S; r--)
	    Y[r+1] = origY[r];

	  if(DEBUG>=2) assert(I > IMIN);
	  
	  double Scale = exp(-LogScaleL[I]);
	  double Yik = Yc(Y,I,K);

	  if(DEBUG>=2) assert(IMIN <= S+K && S+K <= IMAX);
	  if(DEBUG>=2) assert(0 < Jmin[S+K] && Jmax[S+K] <= M);
	  if(DEBUG>=2) assert(I-K-1 <= IMAX);

	  /* precompute ScaleLI[P=I-K-1 ..max(IMIN, I-Kt-deltaExtYRef)]  = exp(LogScaleL[P] - LogScaleL[I]) */
	  int Pmin = max(IMIN, I-K-deltaExtYRef);
	  double LogScaleLI = LogScaleL[I];
	  for(int P = Pmin; P < I-K; P++){
	    if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
	    ScaleLI[P] = exp(LogScaleL[P] - LogScaleLI);
	    if(DEBUG>=2) assert(isfinite(ScaleLI[P]));
	  }

	  int JMAX = (EXTEND && I-K <= DELTA_YU) ? (DELTA_XU + EXTEND_MAX) : DELTA_XU;

	  if(DEBUG/* HERE >=2 */) assert(0 < jmin && jmin <= Jmax[S+K] && Jmax[S+K] <= M);

	  for(int J = jmin /* WAS Jmin[S+K] */; J <= Jmax[S+K]; J++){
	    if((Jmin[S+K] <= J && K <= Kmax[S+K]) || J > JMAX){ // NEW
	      int Umin = I, Umax = 0, Ub = -1;
	      if(J <= JMAX){/* otherwise Umin,Umax are ignored by ALq() */
		if(DEBUG>=2) assert(IMIN <= S+K && S+K <= IMAX && 0 <= K-1 && K-1 <= Kmax[S+K] && Jmin[S+K] <= J && J <= Jmax[S+K]);
		int ULmin = A.ULmin(S+K,K-1,J), ULmax = A.ULmax(S+K,K-1,J);
		int ULb = max(1, A.Lij(S+K,K-1,J)) - 1;
		if(DEBUG>=2) assert(ULmin <= ULmax);
		if(ULmin <= ULmax){
		  if(DEBUG>=2) assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		  Umin = min(Umin, ULmin);
		  Umax = max(Umax, ULmax);
		  Ub = (Ub < 0) ? ULb : (Ub + ULb)/2;
		}

		if(DEBUG>=2) assert(IMIN <= S+K && S+K <= IMAX && 0 <= K && K <= Kmax[S+K] && Jmin[S+K] <= J && J <= Jmax[S+K]);
		ULmin = A.ULmin(S+K,K,J); ULmax = A.ULmax(S+K,K,J);
		ULb = max(1, A.Lij(S+K,K,J)) - 1;
		if(DEBUG>=2) assert(ULmin <= ULmax);
		if(ULmin <= ULmax){
		  if(DEBUG>=2) assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		  Umin = min(Umin, ULmin);
		  Umax = max(Umax, ULmax);
		  Ub = (Ub < 0) ? ULb : (Ub + ULb)/2;
		}
		if(DEBUG>=2) assert(Umin <= Umax);
		if(Umin > Umax)
		  Ub = Umin = Umax = I-K-1;
	      }

	      ALTH[K][J][h] = ALq<OUT_TYPE,OUT_OLD>(m,I,K,J,Umin,Ub,Umax,Yik,Ht,I-K,Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	    } else
	      ALTH[K][J][h] = ALq<OUT_TYPE,OUT_OLD>(m,I,K,J,0,0,-2,Yik,Ht,I-K,Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	  }
	  /* No PMik : matches AL[S+K+1][K][J] on line 6 Equation 10 */

	  /* restore Y[S+K+1 .. S+2] */
	  for(int r = S+K; r > S; r--)
	    Y[r+1] = origY[r+1];
	} /* K = 1 .. Ks */

	/* restore Y[S+1] */
	Y[S+1] = origY[S+1];

	I = S+1;

	/* compute ARTH[K=1..Ks][J= jmin..Jmax[S+K]][h] == A[S+K][K][J].AR AND PMT[K] == PM(S+K,K) with Y[S] = Ht (Line 6 of Equation 10) */	
	Y[I=S] = Ht;
	for(int K = 1; K <= Ks; K++){
	  I = S+K;
	  if(DEBUG>=2) assert(I >= IMIN);

	  /* precompute ScaleRI[P= I+1 .. min(IMAX, I + KMAX + deltaExtYRef)] = exp(LogScaleR[P] - LogScaleR[I]) */
	  int Pmax = min(IMAX, I + KMAX + deltaExtYRef);
	  double LogScaleRI = LogScaleR[I];
	  for(int P = I+1; P <= Pmax; P++){
	    if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
	    ScaleRI[P] = exp(LogScaleR[P] - LogScaleRI);
	    if(DEBUG>=2) assert(isfinite(ScaleRI[P]));		
	  }

	  if(DEBUG>=2) assert(IMIN <= I && I <= IMAX);
	  double Scale = exp(-LogScaleR[I]);	  
	  double Yik = Yc(Y,I,K);

	  int JMIN = (EXTEND && I >= N+1-DELTA_YU) ? (M+1 -DELTA_XU - EXTEND_MAX) : (M+1 -DELTA_XU);

	  for(int J = jmin /* WAS Jmin[S+K] */; J <= Jmax[S+K]; J++){
	    if((Jmin[S+K] <= J && K <= Kmax[S+K]) || J < JMIN){// NEW
	      int Umin = N+1, Umax = I, Ub = -1;
	      
	      if(J >= JMIN){/* otherwise Umin,Umax are ignored by ARq() */
		if(DEBUG>=2) assert(IMIN <= S+K && S+K <= IMAX && 0 <= K-1 && K-1 <= Kmax[S+K] && Jmin[S+K] <= J && J <= Jmax[S+K]);
		int URmin = A.URmin(S+K,K-1,J), URmax = A.URmax(S+K,K-1,J);
		int URb = min(N, A.Rij(S+K,K-1,J)) + 1;
		if(DEBUG>=2) assert(URmin <= URmax);
		if(URmin <= URmax){
		  if(DEBUG>=2) assert(I < URmin && URmin <= URmax && URmax <= N+1);
		  Umin = min(Umin, URmin);
		  Umax = max(Umax, URmax);
		  Ub = (Ub < 0) ? URb : (Ub + URb)/2;
		}

		if(DEBUG>=2) assert(IMIN <= S+K && S+K <= IMAX && 0 <= K && K <= Kmax[S+K] && Jmin[S+K] <= J && J <= Jmax[S+K]);
		URmin = A.URmin(S+K,K,J), URmax = A.URmax(S+K,K,J);
		URb = min(N, A.Rij(S+K,K,J)) + 1;
		if(DEBUG>=2) assert(URmin <= URmax);
		if(URmin <= URmax){
		  if(DEBUG>=2) assert(I < URmin && URmin <= URmax && URmax <= N+1);
		  Umin = min(Umin, URmin);
		  Umax = max(Umax, URmax);
		  Ub = (Ub < 0) ? URb : (Ub + URb)/2;
		}

		if(DEBUG>=2) assert(Umin <= Umax);
		if(Umin > Umax)
		  Ub = Umin = Umax = I+1;
	      }

	      ARTH[K][J][h] = ARq<OUT_TYPE,OUT_OLD>(m,I,K,J,Umin,Ub,Umax,Yik,Ht,I-K,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	    } else
	      ARTH[K][J][h] = ARq<OUT_TYPE,OUT_OLD>(m,I,K,J,0,0,-2,Yik,Ht,I-K,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax);
	  }
	  /* No PMik : matches AR[S+K][K][J] on line 6 Equation 10 */

	  PMTH[K][h] = PM(0,S+K,K,Y);
	}/* K = 1 .. Ks */

	/* restore Y[S] */
	Y[S] = origY[S];
      }

      for(int h = 0; h < VEC_SIZ; h++){
	int t = T + h;
	if(t >= addcntI[S])
	  break;
	int Ks = KsVec[t];

	/* add to LRaddT[S][t] : sum(J=jmin..Jmax[S+K],K=1..Ks) ALTH[K][J][h] * PMTH[K][h] * ARTH[K][J][h] */

	for(int K = 1; K <= Ks; K++){
	  if(DEBUG>=2) assert(S+K >= IMIN);

	  double LR = 0.0;
	  for(int J = jmin; J <= Jmax[S+K]; J++)
	    LR += ALTH[K][J][h] * ARTH[K][J][h];

	  if(DEBUG>=2) assert(IMIN <= S+K && S+K <= IMAX);

	  LR *= PMTH[K][h] * exp(LogScaleL[S+K+1]+LogScaleR[S+K]-LS);

	  if(DEBUG>=2) assert(isfinite(LR) && LR >= 0.0);
	  if(DEBUG>=2) assert(Amin <= S && S <= min(N,Amax));

	  LRaddT[S][t] += LR;

	  if(DEBUG>=2) assert(isfinite(LRaddT[S][t]) && LRaddT[S][t] >= 0.0);
	  if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT && LR != 0.0){
	    double ScaleS = exp(LogScaleL[S+K+1] + LogScaleR[S+K] - LS);
	    FL128 sum = LRadd[0],sum2 = LRadd2[0];
	    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
	      sum += LRadd[i];
	      sum2 += LRadd2[i];
	    }
	    FL128 sum3 = sum;
	    if(Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
	      sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
	      printf("m=%d:S=%d,t=%d(T=%d):I=%d,K=%d,J=%d..%d:sum(ALT[K][J]*PMT[K]*ART[K][J])=%0.10e(ScaleS=%0.10e):LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
		     m,S,t,addcntIcum[S]+t,S+K+1,K,jmin,Jmax[S+K],(double)LR,(double)ScaleS,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
	      double LR2 = 0.0;
	      for(int J = jmin; J <= Jmax[S+K]; J++){
		LR2 += ALTH[K][J][h] * ARTH[K][J][h] * PMTH[K][h] * ScaleS;
		printf("    J=%d:ALT[K][J]*PM=%0.10e(LS=%0.10f),ART[K][J]*PM=%0.10e,PMT[K]=%0.10e:sum(ALT*ART*PMT)=%0.10e(delta=%0.10e)\n",
		       J,(double)(ALTH[K][J][h]*PMTH[K][h]),LogScaleL[S+K+1]+LogScaleR[S+K]-LS, (double)(ARTH[K][J][h]*PMTH[K][h]), (double)PMTH[K][h], 
		       (double)LR2, (double)(ALTH[K][J][h] * ARTH[K][J][h] * PMTH[K][h] * ScaleS));
	      }
	    } else
	      printf("m=%d:S=%d,t=%d(T=%d):I=%d,K=%d,J=%d..%d:sum(ALT[K][J]*PMT[K]*ART[K][J])=%0.10e(ScaleS=%0.10e):LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
		     m,S,t,addcntIcum[S]+t,S+K+1,K,jmin,Jmax[S+K],(double)LR,(double)ScaleS,(double)sum,(double)sum2,(double)(sum-sum2));
	    fflush(stdout);
	  }	    
	}
      }/* t = T .. min(T + VEC_SIZ, addcntI[S])-1 */
    } /* for(T = 0; T < addcntI[S]; T += VEC_SIZ) */
  } // S = IMIN .. IMAX-1 
}

int mverb_numthreads=1;

static int mprobeval_cnt = 0;

// NOTE : class CdoubleBlock and CdoubleAllocation are declared in probeval.h

std::vector<double *> DoubleAllocations;
std::vector<CdoubleAllocation *> SwappedAllocations;/* array of thread specific swapped allocation lists */
CdoubleBlock **blockList = 0;
int blockListLen = 0, blockListMax = 0;
int blockListSwapped = 0;

size_t firstblocksiz = 0;/* excess size of (first) larger block(s) (in doubles) */

/* lightweight block allocator for small double arrays:
   p : pointer to current list of memory blocks used by current thread.
   RefLoc : Left end location of on Reference (in kb) of map requesting allocation.
   ppBlock : If != 0 : save pointer to memory block that contains the returns array. The caller can use this pointer to check ppBlock[0]->swapped : If 0, the block is still in memory, otherwise caller can read
       in the missing data of ppBlock[0]->siz bytes from ppBlock[0]->fp at offset ppBlock[0]->offset to ppBlock[0]->block[0..siz-1], then reset ppBlock[0]->swapped to 0.
 */

static double lastVmRSS = 0.0;

static inline double *thread_alloc(size_t siz, CdoubleAllocation *p, double RefLoc, CdoubleBlock** ppBlock, int tid )
{
  if ( siz > p->tsiz){/* start a new memory block : the remainder of the last block is wasted */
    if(DEBUG) assert(p->tid >= 0 && p->tid == tid);

    if(MPROBEVAL_SWAPOUT && !MDEBUG && p->cnt > 1 && p->Allele && p->filename[0]){ /* If needed, write out the previous memory block to disk and call madvice(DONTNEED) on it to free up memory */
      double wt = wtime();
      if(wt > lastVmRSS + 1.0){
	#pragma omp critical
	{
	  getmem(VmSize,VmRSS,VmSwap);
	  lastVmRSS = wt;
	}
      }
      if(VmRSS + VmSwap > MaxMem * (double) (1e7 * VMRSS_MIN)){/* VmRSS memory use exceeds threshold */
	CdoubleBlock *block = &p->blocks[p->cnt-1];
	if(DEBUG>=2) assert(p->threadmem == &block->block[block->siz - p->tsiz]);

	if(!p->fp){
	  p->fileopen();
	  if(DEBUG) assert(p->fp);
	}

	block->writeout(p->tsiz,p->fp,p->offset,RefLoc,p->filename,p->tid);

	p->offset += block->siz;

	size_t msiz = (block->siz + p->tsiz) * sizeof(double);
	if(madvise(block->block, msiz, MADV_DONTNEED)){
	  int eno = errno;
	  char *err = strerror(eno);
	
#pragma omp critical
	  {
	    printf("madvise(%p,%lu,MADV_DONTNEED) failed:filename=%s, errno=%d:%s\n",block->block,msiz,p->filename,eno,err);
	    fflush(stdout);exit(1);
	  } 
	} else if(VERB >= 3){
	  printf("tid=%d:madvise(%p,%lu,MADV_DONTNEED) succeeded: filename=%s\n",p->tid,block->block,msiz,p->filename);
	  fflush(stdout);
	}
      }
    }

    p->tsiz = VMEM_MINIMIZE * 1024;
    if(siz > p->tsiz)
      p->tsiz = ((siz*sizeof(double) + (PAGE-1)) & ~(PAGE-1)) / sizeof(double);
    
    char *ptr = (char *)mmap(NULL, p->tsiz * sizeof(double), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if(ptr == MAP_FAILED){
      int eno = errno;
      char *err = strerror(eno);
      
      #pragma omp critical
      {
	printf("thread_alloc:tid=%d: mmap of %lu bytes failed:errno=%d:%s\n",p->tid, p->tsiz * sizeof(double),eno,err);
	dumpmemmap();
	fflush(stdout);exit(1);
      }
    }
    p->threadmem = (double *)ptr;

    if(p->cnt >= p->mcnt){/* reallocate p->blocks[] */
      unsigned int mcnt = p->mcnt * 2 + 1;
      if(DEBUG && mcnt < p->mcnt){
        #pragma omp critical
	{
	  printf("thread_alloc:tid=%d:p->mcnt= %u -> %u (overflow)\n",p->tid,p->mcnt,mcnt);
	  fflush(stdout);exit(1);
	}
      }
      CdoubleBlock *blocks = new CdoubleBlock[mcnt];
      memcpy(blocks,p->blocks,p->cnt * sizeof(CdoubleBlock));
      delete [] p->blocks;
      p->blocks = blocks;
      p->mcnt = mcnt;
    }
    CdoubleBlock *block = &p->blocks[p->cnt++];
    block->block = p->threadmem;
    block->siz = p->tsiz;
    block->swapout = 0;
  }
   
  double *ret = p->threadmem;
  p->threadmem += siz;
  p->tsiz -= siz;
  if(ppBlock)
    ppBlock[0] = &p->blocks[p->cnt - 1];
  
  return ret;
}


/* qsort compare function to sort pointers to memory blocks that have been swapped out in Descending order of RefLoc value */
static inline int RefLocDec(CdoubleBlock **pp1, CdoubleBlock **pp2)
{
  return (int)pp2[0]->RefLoc - (int)pp1[0]->RefLoc;
}

/* create a sorted List of memory blocks that have been swapped out to disk */
void SwappedAllocations_ListSort()
{
  blockListMax = 0;
  for(std::vector<CdoubleAllocation *>::iterator it = SwappedAllocations.begin(); it != SwappedAllocations.end(); ++it){
    CdoubleAllocation *p = *it;

    if(DEBUG) assert(p->cnt >= 1);// the first block should always be present

    blockListMax += max(0, p->cnt - ((p->tid >= 0) ? 2 : 1));
  }
  
  if(VERB>=3 && blockListMax > 0){
    printf("SwappedAllocations_ListSort(): found %d memory blocks that may have been swapped out to disk\n",blockListMax);
    fflush(stdout);
  }

  delete [] blockList;
  blockList = new CdoubleBlock*[blockListMax];

  blockListLen = 0;
  for(std::vector<CdoubleAllocation *>::iterator it = SwappedAllocations.begin(); it != SwappedAllocations.end(); ++it){
    CdoubleAllocation *p = *it;
    
    int imax = max(1, p->cnt) - 1;
    for(int i = (p->tid >= 0) ? 1 : 0; i < imax; i++){
      CdoubleBlock *pblock = &p->blocks[i];
      if(pblock->swapout){
	if(DEBUG>=2) assert(blockListLen < blockListMax);
	blockList[blockListLen++] = pblock;
      }
    }
    if(DEBUG) assert(blockListLen <= blockListMax);
  }

  if(VERB>=3 && blockListLen > 0){
    printf("SwappedAllocations_ListSort(): created combined list of %d/%d memory blocks that are swapped out to disk\n",blockListLen,blockListMax);
    fflush(stdout);
  }

  qsort(blockList, blockListLen, sizeof(CdoubleBlock *), (intcmp*)RefLocDec);

  if(VERB/* HERE >=2 */ && blockListLen > 0){
    printf("\nSwappedAllocations_ListSort(): sorted list of %d block pointers by RefLoc (blockListSwapped = blockListLen = %d): cum wall= %0.6f secs\n",blockListLen,blockListLen,wtime());
    fflush(stdout);
  }

  blockListSwapped = blockListLen;
  if(DEBUG>=2) assert(blockListSwapped == SwappedAllocations_swapcnt());
}

// read in any double blocks that have been written out to disk
void SwappedAllocations_restore()
{
  size_t siz = 0;

  #pragma omp parallel num_threads(MaxThreads)
  {
    size_t mysiz = 0;
    size_t mycnt = 0;

    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num ();
#endif

    for(int i = 0; i < blockListLen; i++){
      CdoubleBlock *p = blockList[i];
#ifdef _OPENMP
      if(p->tid != tid)// Each thread will only read memory blocks that it wrote out, avoiding shared access to the file descriptor
	continue;
#endif // _OPENMP
      p->readin();
      mysiz += p->siz;
      mycnt++;
    }

    #pragma omp atomic
    blockListSwapped -= mycnt;

    #pragma omp atomic
    siz += mysiz;
  } // parallel
  if(DEBUG) assert(blockListSwapped == 0);
  if(DEBUG>=2) assert(blockListSwapped == SwappedAllocations_swapcnt());

  if(VERB/* HERE >=2 */ && siz > 0){
    printf("Read in %u blocks totalling %0.4f Gb from swap files using %d threads: cum wall= %0.6f secs\n", blockListLen,siz*sizeof(double)*1e-9,MaxThreads,wtime());
    fflush(stdout);
  }
}

// return count of memory blocks that are swapped out and NOT in memory (swapout == 1)
int SwappedAllocations_swapcnt()
{
  int cnt = 0;
  for(int i = 0; i < blockListLen; i++)
    if(blockList[i]->swapout > 0)
      cnt++;
  return cnt;
}

/* apply madvice(DONTNEED) to memory blocks blockList[blockListLen-1 .. T] that are backed up on disk, until the number of swapped out blocks (blockListSwapped) is at least MinSwapCnt */
void SwappedAllocations_DONTNEED(int MinSwapCnt)
{
  double wt1 = wtime();

  if(DEBUG && !(0 <= blockListSwapped && blockListSwapped < MinSwapCnt && MinSwapCnt <= blockListLen && blockListLen <= blockListMax)){
    printf("SwappedAllocations_DONTNEED(MinSwapCnt= %d): blockListSwapped= %d, blockListLen= %d, blockListMax= %d, LAZYRESTORE_MIN= %d,MAX= %d\n",
	   MinSwapCnt,blockListSwapped,blockListLen,blockListMax, LAZYRESTORE_MIN, LAZYRESTORE_MAX);
    fflush(stdout);
    assert(0 <= blockListSwapped && blockListSwapped < MinSwapCnt && MinSwapCnt <= blockListLen && blockListLen <= blockListMax);
  }

  int origblockListSwapped = blockListSwapped;

  size_t totsiz = 0;

  for(int i = blockListLen; --i >= 0;){
    CdoubleBlock *p = blockList[i];
    if(DEBUG/* HERE >=2 */) assert(p->swapout);
    if(p->swapout < 0){
      size_t msiz = (p->siz * sizeof(double) + (PAGE-1)) & ~(PAGE-1);
      if(madvise(p->block, msiz, MADV_DONTNEED)){
	int eno = errno;
	char *err = strerror(eno);

	#pragma omp critical
	{
	  printf("madvise(%p,%lu,MADV_DONTNEED) failed: fd=%d: errno=%d:%s\n",p->block,msiz,fileno(p->fp),eno,err);
	  fflush(stdout);exit(1);
	}
      }
      totsiz += msiz;

      if(DEBUG>=2) assert(blockListSwapped == SwappedAllocations_swapcnt());
      p->swapout = 1;
      if(++blockListSwapped >= MinSwapCnt)
	break;
    }
  }
  if(DEBUG) assert(MinSwapCnt <= blockListSwapped && blockListSwapped <= blockListLen);
  if(DEBUG>=2) assert(blockListSwapped == SwappedAllocations_swapcnt());

  if(VERB>=1+RELEASE && (totsiz > 0 || blockListSwapped > origblockListSwapped)){
    getmem(VmSize,VmRSS,VmSwap);
    double wt = wtime();
    printf("Swapped out memory blocks = %u -> %u/%u (bytes removed= %0.3f G):VmRSS=%0.3f,VmSwap=%0.3f,VmSize=%0.3f, MaxMem= %0.3f Gb:wall=%0.6f(cum=%0.6f)\n",
      origblockListSwapped,blockListSwapped,blockListLen,totsiz*1e-9,VmRSS*1e-9,VmSwap*1e-9,VmSize*1e-9,MaxMem,wt-wt1,wt);
    fflush(stdout);
  }
}

void DoubleAllocations_free()
{
  if(VERB>=2){
    long long memsiz = 0, cnt = 0;

    if(DoubleAllocations.size() > 0){
      memsiz = firstblocksiz + ((long long)DoubleAllocations.size()) * VMEM_MINIMIZE * 1024;
      cnt = DoubleAllocations.size();
    }

    if(SwappedAllocations.size() > 0){
      for(std::vector<CdoubleAllocation *>::iterator it = SwappedAllocations.begin(); it != SwappedAllocations.end(); ++it){
        CdoubleAllocation *p = *it;
	for(unsigned int i = (p->tid >= 0) ? 1 : 0; i < p->cnt; i++, cnt++)
	  memsiz += p->blocks[i].siz;
      }
    }

    if(cnt > 0 || memsiz > 0){
      printf("mprobeval:Freeing %lld blocks with total size of %llu bytes, firstblocksiz=%lu bytes\n",cnt, memsiz * sizeof(double), firstblocksiz * sizeof(double));
      fflush(stdout);
    }
  }

  for(std::vector<double *>::iterator it = DoubleAllocations.begin(); it != DoubleAllocations.end(); ++it)
    (void)free(*it);

  DoubleAllocations.clear();
  firstblocksiz = 0;

  for(std::vector<CdoubleAllocation *>::iterator it = SwappedAllocations.begin(); it != SwappedAllocations.end(); ++it){
    CdoubleAllocation *p = *it;
    delete p;
  }

  SwappedAllocations.clear();

  delete [] blockList; blockList = 0;
  blockListSwapped = blockListLen = blockListMax = 0;
}


static Csetlimit *glimit = NULL;// global copy of mprobeval() arg "limit", for use by qsort

/* qsort compare function to order map index by location of map in contig */
 static inline int RefSiteInc(int *p1, int *p2)
 {
    int m1 = *p1;
    int m2 = *p2;
    Csetlimit *plimit1 = &glimit[m1];
    Csetlimit *plimit2 = &glimit[m2];
    
    return (plimit1->ileft > plimit2->ileft) ? 1 : (plimit1->ileft < plimit2->ileft) ? -1 :
      (plimit1->ileft < 0 || plimit2->ileft < 0) ? (m1-m2) : /* preserve original order IF glimit[m] has not been computed for either map */
      (plimit1->iright > plimit2->iright) ? 1 : (plimit1->iright < plimit2->iright) ? -1 :
      m1 - m2;/* preserve original order IF tied */
 }

extern int localINITIAL_DELTA_RANGE;/* see refine.h */


/** return LP (log likelihood) of current map Y[] and single changes to Y[] */
/** This assumes mprobeval and qprobeval are not called recursively */

// NOTE : If Allele > 0, some of the memory of newLPdelta[] or newLPa[] may be written out to disk and needs to be read back in by calling SwappedAllocations_restore() after mprobeval() call(s) for all Alleles

// NOTE : compilation of mprobeval() template is too slow with Intel compiler. Instead the template specialization is done at critical points only

#if USE_MIC==0
template<int OUT_TYPE, int OUT_OLD>
static __attribute__ ((noinline))
#endif
double mprobeval(int n, double *Hcuts,/**< Hcuts[i=0..n+1] : complete consensus map : only used to set alignment bounds */
		 int N, double *origY,/**< Y[I=0..N+1] : current consensus map hypothesis */
		 int MD, /**< number of maps */
		 int *MX, double **X, /**< X[m=0..MD-1][j=0..MX[m]+1] : Map sites */
		 int lc,/**< If left end of Y is a linear chromosome end */
		 int rc,/**< If right end of Y is a linear chromosome end */
		 int **map, /**< map[m=0..MD-1][j=0..MX[m]+1] is index in Hcuts[0..n+1] of X[m][j], 
					map[MD][j=0..n+1] is index in Y[0..N+1] of Hcuts[j] */
		 int **mapK, /**< mapK[m=0..MD-1][j=0..MX[m]+1] is K index offset in Hcuts[0..n+1] of X[m][j] */ 
		 Csetlimit *limit, /* setlimit() cache : a function of map[], mapK[] and Hcuts[] and slightly dependent on nmap[MD],nmapK[MD] : Hcuts[imin,imax] correspond to Y[IMIN-LRANGE,IMAX+LRANGE] */
		 int *nmapMD, /**< nmapMD[I=0..N+1] is index in Hcuts[0..n+1] of Y[I] and always defined */
		 double *TBmapWT, /**< If != 0 : TBmapWT[m] is a the -TB based weight (multiplied by pcontig->mapWT[m] due to -Erefine or -BestRefWT) 
				     for map m to be applied when summing up the total log(LR+LRbias) value(s) */
		 double *newLP,/**< return value of Log Likelihood for each Map */
		 double *oldLP,/**< If != 0 : previous estimate of Log Likelihood for each map */
		 double Ymin, double Ymax, /**< location range in Y[] that changed since oldLP was computed (or WILL change due to addloc[] or delta[][])  : only valid if oldLP != 0 && Ymin >= 0.0 */
		 int addcnt, double *addloc,/**< addloc[t=0..addcnt-1] are the locations at which a site will be added to Y[] (one at a time) : Y[0] < addloc[t] < Y[N+1] */
		 int *D, double **delta,/**< If LPdelta != 0 : delta[I=0..N][d=0..D[I]-1] is the amount by which the interval Y[I+1]-Y[I] is incremented (NOTE : D[I] can be 0 or -1, meaning no values) */
		 double *LPdelG,/**< If != 0 : LPdelG[I=1..N] is the LP value when the site Y[I] is deleted from Y[0..N+1] (NOTE: LPdelG[I] for Y[I] outside of Ymin..Ymax are NOT computed) */
		 double *LPaddG,/**< If addcnt != 0 : LPaddG[t=0..addcnt-1] is the LP value when a site is added at addloc[t] to Y[0..N+1] */
		 double *LPgrad[2],/**< If != 0 : LPgrad[d=0,1][I=0..N] is (d+1)'th gradient of LP vs Y[I+1]-Y[I] (Not yet implemented) */
		 double **LPdelta, /**< If != 0 : LPdelta[I=0..N][d=0..D[I]-1] is the LP value when Y[I+1]-Y[I] is incremented by delta[I][d] */
		 double **newLPd, /**< If != 0 : newLPd[m=0..MD-1][I=DminM[m]..DmaxM[m]] is the LP value of map m when the site Y[I] is deleted from Y[0..N+1]
				     Values newLPd[m][I] for 1 <= I < DminM[m] are not defined but should be assumed equal to newLPd[m][0] 
				     Values newLPd[m][I] for DmaxM[m] < I <= N are not defined but should be assumed equal to newLPd[m][DmaxM[m]] (NOTE: newLPd[m][I] for Y[I] outside of Ymin..Ymax are NOT computed
				     Also if SKIPMAP_ADD && addcnt > 0 && Ymin < 0.0 AND no addloc[t] value falls within Y[I-1]..Y[I+1], then newLPd[m][I] is forced to newLPd[m][0]).
				     NOTE: If VMEM_MINIMIZE, the 2nd index of newLPd[][] has not been allocated and newLPd[m][0] is replaced by newLPd0[m] */
                 double *newLPd0,/* IFF VMEM_MINIMIZE then newLPd0[m] replaces newLPd[m][0] (which is no longer defined, unless DminM[m] <= 0) */
		 int **newLPdPr,  /**< If MDEBUG : newLPdPr[m=0..MD-1][I=DminM[m]..DmaxM[m]] is 1 IFF newLPd[I][m] is low precision (due to drop in LP) */
		 int *DminM,  /**< DminM[m=0..MD-1] : see newLPd[] */ 
		 int *DmaxM,  /**< DmaxM[m=0..MD-1] : see newLPd[] */
		 double **newLPa, /**< If != 0 : newLPa[m=0..MD-1][t=TminM[m]..TmaxM[m]] is the LP value of map m when a site is added at addloc[t] to Y[0..N+1]
				     Values newLPa[m][t] for 0 <= t < TminM[m] are not defined and should be assumed equal to newLPa[m][TminM[m]]
				     Values newLPa[m][t] for TmaxM[m] < t < addcnt are not defined and should be assumed equal to newLPa[m][TmaxM[m]]
				     NOTE: If VMEM_MINIMIZE, the 2nd index of newLPa[][] has not been allocated. The allocation list is stored in DoubleAllocation[] */
		 int **newLPaPr,  /**< If MDEBUG : newLPaPr[m=0..MD-1][t=TminM[m]..TmaxM[m]] is 1 IFF newLPa[t][m] is low precision (due to drop in LP)   */
		 int *TminM,  /**< TminM[m=0..MD-1] : see newLPa[] */ 
		 int *TmaxM,  /**< TmaxM[m=0..MD-1] : see newLPa[] */
		 double **newLPdelta, /**< If != 0 : newLPdelta[m=0..MD-1][Dcum[I=AminM[m]..min(N-1,AmaxM[m])] + (d=0..D[I]-1) - Doffset] is the LP value of map m when Y[I+1]-Y[I] is incremented by delta[I][d] 
					 Values newLPdelta[m][] for I < AminM[m] or I > min(N-1,Amax[m]) are not defined and assumed to match newLP[m] regardless of delta[I][d]  (presumably molecule m does not overlap Y[I..I+1])
					 To save real memory, the offset Doffset == Dcum[AminM[m]] is used to shift the used memory to the left end of each m's memory block.
					 NOTE1: LPdelta != 0 implies newLPdelta != 0, but it is possible that LPdelta == 0 && newLPdelta != 0
					 NOTE2: If VMEM_MINIMIZE, only the first index has been allocated and the 2nd index is allocated for each m by this function. The allocation list is stored in DoubleAllocation[] */
		 int *Dcum,   /* Dcum[I=0..N] = sum(t = 0..I-1) max(0,D[I]) */
		 /* saves real memory by using flattened newLPdelta[m=0..MD-1][d = 0 .. (sum(I=0..N)max(0,D[I]))-1], which saves most of the 2nd index memory that is never de-referenced */
		 int *AminM, /**< AminM[m=0..MD-1] : see newLPdelta[] */
		 int *AmaxM, /**< AmaxM[m=0..MD-1] : see newLPdelta[] */
		 Ccontig *pcontig, /**< point to complete contig information for debugging */
		 char *skipmap, /** If != 0 : skipmap[m=0..MD-1] == 1 IFF map m is already known not to overlap intervals Y[I..I+1] with D[I] > 0. skipmap[] is updated with each call */
		 int Allele, /* If != 0, identifies the Allele that Hcuts[] and origY[] correspond to (1 or 2 if called from HaploTypeR). Result arrays newLPdelta[], newLPa[] etc may be swapped out and
			       caller must call SwappedAllocations_restore() after mprobeval() for both Alleles have been called */
		 CdoubleBlock **pBlock /* If != 0, zero'd memory pBlock[0..MD-1], that is used to save a pointer to the CdoubleBlock used by newLPdelta[m][] or newLPa[m][], so it can be swapped out and restored on demand
					   by the caller, instead of all at once by calling SwappedAllocations_restore() */
		 )
{
  START(mprobeval);

  DELTA_YU = FIX_U ? deltaExtYRef + FIX_U - 1 : DELTA_Y;// NEW46
  DELTA_XU = FIX_U ? deltaExtXRef + FIX_U - 1 : DELTA_X;// NEW46
  if(DEBUG) assert(DELTA_YU >= DELTA_Y);
  if(DEBUG) assert(DELTA_XU >= DELTA_X);

  if(DEBUG && LPdelta != 0) assert(newLPdelta != 0);

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
  if(DEBUG) assert(MaxLR * MaxLR <= MaxFloat);

  if(PoutlierEnd < MinUR){
    printf("mprobeval : PoutlierEnd = %0.3e is less than MinUR = %0.3e. To avoid floating point underflow please raise -endoutlierRef and -endoutlierSwitch values to at least MinUR\n",
	   PoutlierEnd, MinUR);
    fflush(stdout);exit(1);
  }

  if(DEBUG) assert(deltaExtYRef >= DELTA_Y);
  if(DEBUG) assert(deltaExtXRef >= DELTA_X);

  if(EPS_ADJUST && PQ_BND && EPS3 > 0.0 && EPS2 <= 0.0){
    printf("mprobeval() : EPS2=0, EPS3= %0.3e : EPS_ADJUST will not work with EPS2=0 and EPS3 > 0. Please use a non-zero value of EPS2\n",EPS3);
    fflush(stdout);exit(1);
  }
  assert(0 <= PQ_BND && PQ_BND <= 1);
  if(PQ_BND && VERB>=2){
    printf("EPS= %0.3e, EPS2= %0.3e, EPS3= %0.3e\n",EPS,EPS2,EPS3);
    fflush(stdout);
  }

  mprobeval_cnt++;
  if(VERB>=2){
    printf("mprobeval:mprobeval_cnt=%d,n=%d,N=%d,addcnt=%d\n",mprobeval_cnt,n,N,addcnt);
    fflush(stdout);
  }
  if(DEBUG) assert(N > 0);

  if(MDEBUG && MDEBUG_M>=0){
    printf("\nMDEBUG_M=%d(MD=%d),_S=%d,_F=%d,_N=%d(N=%d),_T=%d(addcnt=%d),_TS=%d,_TT=%d\n",
	   MDEBUG_M,MD,MDEBUG_S,MDEBUG_F,MDEBUG_N,N,MDEBUG_T,addcnt,MDEBUG_TS,MDEBUG_TT);
    fflush(stdout);
  }
  if(!ISLINEAR && (rc || lc)){
    printf("mprobeval:rc=%d,lc=%d not supported (set ISLINEAR in constants.h)\n",rc,lc);
    exit(1);
  }
  if(DEBUG && !LPdelG) assert(!addcnt);

  int *maporder = NULL;
  if(MPROBEVAL_MAPORDER && SETLIMIT_CACHE && MP_MIN_TIME > 0 && limit){

    maporder = new int[MD];
    for(int m = 0; m < MD; m++)
      maporder[m] = m;

    glimit = limit;

    qsort(maporder, MD, sizeof(int), (intcmp*)RefSiteInc);
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
	    printf("\nm=%d,J=%d,i=%d,k=%d:lastJ=%d,lasti=%d,lastk=%d\n",m,J,i,k,lastJ,lasti,lastk);
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
  if(DEBUG>=2 && oldLP && Ymin >= 0.0){/* make sure Ymin..Ymax includes all changes addloc[] and delta[][] */
    assert(0.0 <= Ymin && Ymin <= Ymax && Ymax <= origY[N+1]);
    for(int t = 0; t < addcnt;t++)
      assert(Ymin <= addloc[t] && addloc[t] <= Ymax);
    if(newLPdelta)
      for(int I = 0; I <= N; I++)
	if(D[I] > 0)
	  assert(Ymin <= origY[I] && origY[I+1] <= Ymax);
  }

  // NOTE : duplicates can rarely occur on one allele when an interval is completely erased, even if the original map had no duplicates */
  // NOTE : It is even possible that the values in addloc[] are not in ascending order : can be fixed by reordering and keeping track of the index changes so results can be mapped back to the original order
  if(DEBUG && SKIP_DIST_FIX){ /* check that addloc[] values are in ascending order */
    for(int t = 1; t < addcnt;t++){
      if(!(addloc[t] >= addloc[t-1])){
	printf("t=%d/%d:addloc[t-1]= %0.8f, addloct[t]= %0.8f, skip_dist= %0.8f\n",t,addcnt,addloc[t-1],addloc[t], skip_dist);
	fflush(stdout);
	assert(addloc[t] >= addloc[t-1]);
      }
    }
  }

  /*  if(DEBUG>=2 && newLPdelta){
    for(int I = 0; I < N; I++){
      assert(Dcum[I+1]-Dcum[I] == max(0,D[I]));
      assert(D[I] <= (2*localINITIAL_DELTA_RANGE + MAX_INDEL));
    }
    }*/

//  int *Kmax = new int[(N+1)+2*(N+2)];

  double **addlocI = new double*[N+1];/* addlocI[I=0..N][t=0..addcntI[I]-1] are the addloc[] sites in the interval Y[I] .. Y[I+1] */
  int *addcntI = new int[2*(N+2)];/* addcntI[N+1] == 0 */
  int *addcntIcum = &addcntI[N+2];/* addcntIcum[I=0..N+1] == sum(addcntI[0..I-1]) */
  double *addlocImem = NULL;// Only used if PADDING > 0 
  int addcntImax = 0;/* max(addcntI[I=0..N]) */

  FL128 LPaddI;/* LPaddI is the contribution to LPadd[t] that do not depend on the value of addloc[t] (same as LRcum with DELTA_Y reduced by 1) */
  FL128 *LPadd = new FL128[addcnt+N+1]; /* higher precision version of LPaddG[] */
  FL128 **LPaddTI = new FL128*[N+1];/* LPaddTI[I=0..N][t=0..addcntI[I]-1] is a different way of indexing LPadd[0..addcnt-1] (to match addlocI[I][t]) */
  FL128 *LPdel = &LPadd[addcnt];/* higher precision version of LPdelG[] */
  
  if(DEBUG && MDEBUG && !(MDEBUG_F < 0 || MDEBUG_TS < 0)){
    printf("mprobeval:MDEBUG_F=%d,MDEBUG_TS=%d: at least one of these values must be -1\n",MDEBUG_F,MDEBUG_TS);
    exit(1);
  }

  if(addcnt){
    int T = 0;
    int cntIcum = 0;
    int cntIcumPad = 0;
    double *Y = origY;

    if(DEBUG && addloc[addcnt-1] >= Y[N+1]){// NEW2
      if(DEBUG) assert(addloc[addcnt-1] < Y[N+1] + 1e-6);
      printf("\nWARNING: addloc[addcnt-1] = %0.10f, Y[N+1]= %0.10f, Hcuts[n+1]= %0.10f : Increasing Y[N+1] & Hcuts[n+1] to %0.10f\n",addloc[addcnt-1],Y[N+1],Hcuts[n+1], addloc[addcnt-1]+1e-7);
      fflush(stdout);
      Y[N+1] = Hcuts[n+1] = addloc[addcnt-1] + 1e-7;
    }

    for(int I = 0; I <= N; I++){
      if(DEBUG) assert(T==cntIcum);
      addlocI[I] = &addloc[T];
      LPaddTI[I] = &LPadd[T];

      int cntI = 0;
      for(; T < addcnt && addloc[T] < Y[I+1]; T++)
	cntI++;

      addcntIcum[I] = cntIcum;
      addcntI[I] = cntI;
      addcntImax = max(addcntImax,cntI);
      cntIcum += cntI;
      cntIcumPad += (cntI + VPAD) & (~VPAD);
      if(MDEBUG && MDEBUG_M>=0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && MDEBUG_N==N){
	printf("I=%d,N=%d:addcntI[I]=%d,addcntIcum[I]=%d..%d,T=%d,Y[I]=%0.3f,Y[I+1]=%0.3f\n",
	       I,N,addcntI[I],addcntIcum[I],cntIcum,T,Y[I],Y[I+1]);
	if(I==MDEBUG_TS)
	  for(int t = addcntIcum[I]; t < addcntIcum[I] + addcntI[I]; t++)
	    printf("   t=%d:addloc[t]=%0.3f\n",t,addloc[t]);
	fflush(stdout);
      }
    }
    if(DEBUG) assert(cntIcum == T);
    if(DEBUG && !(cntIcum == addcnt)){
      printf("addcnt=%d,cntIcum=%d,T=%d, Y[N+1] = %0.6f\n",addcnt,cntIcum,T, Y[N+1]);
      for(int t = T; T < addcnt; T++)
	printf("addloc[%d] = %0.6f\n",t,addloc[t]);
      fflush(stdout);
      assert(cntIcum == addcnt);
    }

    if(PADDING > 0 && (USE_MIC ? VECTOR_MIC : VECTOR_AVX)){
      POSIX_MEMALIGN(addlocImem,  PADDING * 2, cntIcumPad * sizeof(double));
      int cntIcum = 0;
      for(int I = 0; I <= N; I++){
	int cntI = addcntI[I];
	int cntIpad = (cntI + (PADDING/sizeof(float)) - 1) & ~((PADDING/sizeof(float)) - 1);

	memcpy(&addlocImem[cntIcum], addlocI[I], cntI * sizeof(double));

	double *p = addlocI[I] = &addlocImem[cntIcum];

	if(cntI < cntIpad){
	  if(DEBUG) assert(cntI > 0);
	  double lastloc = p[cntI-1];
	  for(int t = cntI; t < cntIpad; t++)
	    p[t] = lastloc;
	}

	cntIcum += cntIpad;
      }
      if(DEBUG) assert(cntIcum == cntIcumPad);
    }

    addcntIcum[N+1] = addcnt;
    addcntI[N+1] = 0;

    if(DEBUG/* HERE >=2 */){
      for(int S= 0; S <= N; S++){
	for(int t = 0; t < addcntI[S]; t++){
	  if(DEBUG && !(Y[S] < addlocI[S][t] && addlocI[S][t] < Y[S+1])){
	    printf("S=%d,t=%d:addlocI[S][t]=%0.6f,Y[S]=%0.6f,Y[S+1]=%0.6f\n",S,t,addlocI[S][t],Y[S],Y[S+1]);
	    fflush(stdout);
	    assert(Y[S] < addlocI[S][t] && addlocI[S][t] < Y[S+1]);
	  }
	}
      }
    }
  }

  if(DEBUG) assert(resSD[0] > 0.0);
  if(DEBUG) assert(!END_VITERBI);
  if(DEBUG){  /* make sure no real cuts are closer than minKB */
    double *Y = origY;
    for(int I = 1; I < N; I++)
      if(DEBUG && !(Y[I+1]-Y[I] >= minKB - 2e-9)){
	printf("\nqprobeval:I=%d:Y[I]=%0.6f,Y[I+1]=%0.6f, delta=%0.6f, minKB=%0.6f\n",
	       I,Y[I],Y[I+1], Y[I+1]-Y[I], minKB);
	fflush(stdout);
	assert(Y[I+1]-Y[I] >= minKB - 2e-9);
      }
  }


  /* Note: I-Kmax[I] must be monotonically non-decreasing functions of I (unless Kmax[I] < 0)*/
  if(DEBUG>=2){
    int Tmax = 1;
    for(int I= 1; I <= N; I++){
      assert(I-Kmax[I] >= 1);
      assert(I-Kmax[I] >= Tmax);
      if(I-Kmax[I] > Tmax)
	Tmax = I-Kmax[I];
    }
  }

  double LLRsum = 0.0;

  /* Initialize LPdel[0..N] */
  if(LPdelG)
    for(int S = 0; S <= N; S++)
      LPdel[S] = 0.0;

  /* Initialize LPaddI, LPadd[0..addcnt-1] */
  if(addcnt){
    LPaddI = 0.0;
    for(int t = 0; t < addcnt; t++)
      LPadd[t] = 0.0;
  }

  int *mindelta = NULL;

  if(newLPdelta){ 
    mindelta = new int[N+1]; /* mindelta[I] is index of smallest value of delta[I][0..D[I]-1] */

    for(int S = 0; S <= N; S++){
      if(LPdelta){ /* Initialize LPdelta[S=0..N][d=0..D[S]-1] */
	if(DEBUG && D[S] > 0)assert(delta[S] != 0);
	if(D[S] > 0) memset(LPdelta[S], 0, D[S]*sizeof(*LPdelta[S]));
      }

      mindelta[S] = 0;
      if(D[S] <= 0)
	continue;
      if(HFIX){
	double mdelta = delta[S][0];
	for(int d = 0; d < D[S]; d++){
	  if(delta[S][d] < mdelta){
	    mdelta = delta[S][d];
	    mindelta[S] = d;
	  }
	}
      }
    }
    if(DEBUG>=2){
      for(int S = 0; S <= N; S++){
	for(int d = 0; d < D[S]; d++){
	  double del = delta[S][d];
	  if(!(origY[S+1]-origY[S] + del >= 0.0)){
	    printf("S=%d,d=%d,D[S]=%d:Y[S]=%0.6f,Y[S+1]=%0.6f,delta[S][d]=%0.6f=del\n",
		   S,d,D[S],origY[S],origY[S+1],delta[S][d]);
	    fflush(stdout);
	    assert(origY[S+1]-origY[S]+del >= 0.0);
	  }
	  if(0 && d > 0 && !(delta[S][d-1] < del)){
	    printf("S=%d,d=%d,D[S]=%d:delta[S][d-1]=%0.6f,delta[S][d]=%0.6f=del,Y[S]=%0.6f,Y[S+1]=%0.6f(delta=%0.6f)\n",
		   S,d,D[S],delta[S][d-1],delta[S][d],origY[S],origY[S+1],origY[S+1]-origY[S]);
	    for(int t = 0; t < D[S]; t++)
	      printf("    delta[S][%d]= %0.6f\n",t,delta[S][t]);
	    fflush(stdout);
	    assert(delta[S][d-1] < del);
	  }
	}
      }
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

  long long double_pool_size = (5LL * (N+2) + 2LL * (N+1) + (2LL * (KMAX+2))*(maxMX+1) + (EPS2 > 0.0 ? 3*maxMX*N + (N+1) : 0)) * sizeof(double);
  double_pool_size += (KMAX+1) * (2LL * sizeof(double*) + sizeof(double));// NEW91
  if(PQ_BND)
    double_pool_size += deltaExtYRef * deltaExtXRef * sizeof(MFLOAT2) * 2;
  long long int_pool_size = ((EPS2 > 0.0 ? 5LL : 3LL)*(N+1) + 2LL*(maxMX+1))*sizeof(int);
#if VECTOR2_ADD
  double_pool_size += (KMAX+1)*(2*(maxMX+1) * sizeof(double *) + VEC_SIZ * (1 + 2*(maxMX+1)) * sizeof(double)) + PADDING;/* only 1 PADDING allowance since these are allocated back to back in multiples of PADDING */
  int_pool_size += 2LL * addcntImax * sizeof(int *);
  double_pool_size += (KMAX+1) * (2LL * sizeof(double**) + sizeof(double*));
#endif

  long long FL128_addcnt_pool_size = 0;
  if(addcnt > 0){
    for(int I = 0; I <= N; I++){
      if(addcntI[I] <= 0)
	continue;
      FL128_addcnt_pool_size += addcntI[I];
    }
  }
  long long FL128_lpdelta_pool_size = 0;
  if(newLPdelta){
    for(int I = 0; I <= N; I++){
      if(D[I] <= 0)
	continue;
      FL128_lpdelta_pool_size += (D[I] + VPAD) & (~VPAD);
    }
  }
  long long FL128_pool_size = 4*(N+2)*sizeof(FL128) + 
    FL128_addcnt_pool_size * sizeof(double) + (N+1)*sizeof(double *)
    /*    + FL128_lpdelta_pool_size * sizeof(double) + (N+1) * (USE_MIC ? 64 : USE_AVX ? 32 : 8) - 8) */ + (N+1) * sizeof(double *) 
    /*    + FL128_lpdelta_pool_size * sizeof(MFLOAT) + (N+1) * (sizeof(MFLOAT *) + (USE_MIC ? 64 : USE_AVX ? 32 : 8) - 8) */;

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
    KstrideMax = (KstrideMax + ALIGN - 1) & (~(ALIGN - 1));
  long long StrideMax = (1LL + kmax) * KstrideMax;
  if(CACHE > 0 && ALIGN2 > 0)
    StrideMax = (StrideMax + ALIGN2 - 1) & (~(ALIGN2 - 1));
  
  /* HERE : If DIAGONAL && MIN_MEM, Fmem[] and Imem[] are not allocated via heap->alloc() but instead using malloc/free to minimize memory usage
     For predicted virtual memory usage for Fmem[] and Imem[] replace StrideMax with  (1 + kmax) * (IMAX - IMIN + 1) * deltaJ ~ (1 + kmax) * (maxMX + 2*RANGE) * maxMX
     For predicted real memory usage for Fmem[] and Imem[] replace IKcnt * maxMX with (1 + av(Kmax[]) * (av(MX[]) + 2 * RANGE) * av(MX[]) (May use biased av() function, eg mean + sd)
     Test case : see /home/users3/tanantharaman/data/ *24sd5_r57??/output/contigs/exp_refineFinal1/ *_id1.stdout (run with -maxthreads 240 -maxmem 3.6 -maxvirtmem 10.8)
     
     Also use hugepage mmap instead of malloc for lightweight_heap and Fmem/Imem on MIC
  */

  long long strideEstimate = (DIAGONAL && MIN_MEM) ? min(StrideMax, (1LL + kmax) * (maxMX + 2LL * RANGE) * maxMX) : StrideMax;
  long long pool_Fmem_siz = strideEstimate /* WAS StrideMax */ * mDPsizF * sizeof(MFLOAT2);
  long long pool_Imem_siz = strideEstimate /* WAS StrideMax */ * mDPsizI * sizeof(int);
  long long pool_LRdeltaD_siz = FL128_lpdelta_pool_size * sizeof(double);

  int num_allocations = 22 + 2*(N+1) + 2*(KMAX+1) + (EPS2 > 0.0 ? 6 : 0) + (PQ_BND ? 2 : 0);
  long long virtperthread_size = pool_Fmem_siz + pool_Imem_siz + pool_LRdeltaD_siz + double_pool_size + int_pool_size + FL128_pool_size + /* alignment fudge */ num_allocations*(4*CACHE_LINE_GUARD);
  long long perthread_size = (long long)((1.0 + AvKmax) * (AvMX + 2.0 * RANGE) * AvMX) /* WAS IKcnt * maxMX */  * (mDPsizF * sizeof(MFLOAT2) + mDPsizI * sizeof(int))
			  + double_pool_size + int_pool_size + (FL128_pool_size + pool_LRdeltaD_siz) * min(N+1.0, AvMX + 2.0 * RANGE)/(N+1.0) + /* alignment fudge */ num_allocations*(4*CACHE_LINE_GUARD);

  if(DEBUG) assert(perthread_size >= 0);
  double virtmemperthread = (double)(virtperthread_size + OMP_stack_size + 2 * VMEM_MINIMIZE * 1024 * sizeof(double));  /* * 1.1 ? see reset_lightweight_heap */;
  double memperthread = (double)(perthread_size + VMEM_MINIMIZE * 1024 * sizeof(double));  /* * 1.1 ? see reset_lightweight_heap */;
  double maxmem = (MaxMem > 0) ? MaxMem * 1024.0 * 1024.0 * 1024.0 : 1e+13;
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
  if(MaxMem > 0 && (perthread_size > 0 || VMEM_MINIMIZE)){
    /* free lightweight heap : this should release virtual memory space with glibc */
    delete heap;
    heap = new lightweight_heap(0,0);      
    if(TotalThreads){
      getmem(VmSize,VmRSS,VmSwap);
      maxvirtmem -= min(VmSize,VmRSS+VmSwap);
      //      maxmem -= VmRSS; // NOTE : this memory can be allowed to swap out
    }
    int limitthreads = max(1,(int)floor(min(maxmem/memperthread, maxvirtmem/virtmemperthread) + 0.5));
    if(DEBUG) assert(memperthread > 0.0 && maxmem > 0.0 && limitthreads >= 1);
    limitthreads = min(numthreads,limitthreads);
    if(/* 1||*/ limitthreads < numthreads){
      if(VERB && (VERB>=1+RELEASE || DEBUG>=2 || mverb_numthreads)){
	double InvGB = 1e-9;
	printf(" threads reduced from %d to %d due to MaxVirtMem=%0.4f(%0.4f) Gb, MaxMem=%0.4f(%0.4f) Gb, per thread = %0.4f Gb, real = %0.4f Gb\n\t VmRSS=%0.4f,VmSwap=%0.4f,VmSize=%0.4f Gb,Stack=%0.2f Mb:N=%d,IKcnt=%lld,kmax=%d,maxMX=%d,AvKmax=%0.3f,AvMX=%0.1f,RANGE=%d,pools:double=%0.4f,int=%0.4f,FL128=%0.4f(add=%0.4f,delta=%0.4f),res=%0.4f Gb\n",
               numthreads,limitthreads, MaxVirtMem, maxvirtmem*InvGB, MaxMem, maxmem*InvGB, virtmemperthread*InvGB,memperthread*InvGB,VmRSS*InvGB,VmSwap*InvGB,VmSize*InvGB,OMP_stack_size*InvGB*1024.0,
	       N,IKcnt,kmax,maxMX,AvKmax,AvMX,RANGE,double_pool_size*InvGB, int_pool_size*InvGB,FL128_pool_size*InvGB, FL128_addcnt_pool_size*sizeof(double)*InvGB, FL128_lpdelta_pool_size*(sizeof(double)+sizeof(MFLOAT))*InvGB,
	       VMEM_MINIMIZE * 1024 * sizeof(double) * InvGB);
	fflush(stdout);
      }
      if(DEBUG) assert(limitthreads >= 1);
      numthreads = limitthreads;
      extern int verb_numthreads;
      mverb_numthreads = verb_numthreads = 1;
   }
  }
  if(numthreads > MD)
    numthreads = MD;
  if(VERB && (VERB>=2 || mverb_numthreads /* ||LPdelta!=0 */)){
    double InvGB = 1e-9;
    printf(" using %d threads: MaxVirtMem=%0.4f(%0.4f) Gb, MaxMem=%0.4f(%0.4f) Gb, per thread memory = %0.4f Gb, real memory = %0.4f G\n\t VmRSS=%0.4f,VmSwap=%0.4f,VmSize=%0.4f Gb,Stack=%0.2f Mb:N=%d,IKcnt=%lld,kmax=%d,maxMX=%d,AvKmax=%0.3f,AvMX=%0.1f,RANGE=%d,pools:double=%0.4f,int=%0.4f,FL128=%0.4f,res=%0.4f Gb\n",
      numthreads, MaxVirtMem, maxvirtmem*InvGB, MaxMem, maxmem*InvGB, virtmemperthread*InvGB,memperthread*InvGB,VmRSS*InvGB,VmSwap*InvGB,VmSize*InvGB,OMP_stack_size*InvGB*1024.0, 
	   N,IKcnt,kmax,maxMX,AvKmax,AvMX,RANGE,double_pool_size*InvGB,int_pool_size*InvGB,FL128_pool_size*InvGB,VMEM_MINIMIZE*1024*sizeof(double)*InvGB);
    if(VERB/* HERE >=2 */)
      printf("\t StrideMax=%lld,mDPsizF=%lld,mDPsizI=%lld:fudge=%d(%d):total per thread=%lld,%lld,wtime= %0.6f secs\n",
	     StrideMax,mDPsizF,mDPsizI,  num_allocations*(4*CACHE_LINE_GUARD),num_allocations, virtperthread_size,perthread_size,wtime());
    fflush(stdout);
    mverb_numthreads = 0;
  }
  if(DEBUG) assert(numthreads >= 1);

  #endif // _OPENMP

  if(/* (MDEBUG && MDEBUG_M>=0) || */ (VERB>=3 && N==9375 && addcnt==92764 && MD==24210))
    numthreads = 1;

  if(DIAGONAL && MIN_MEM){/* Fmem[], Imem[] and LRdeltaD[] will be allocated inside each thread */
    virtperthread_size -= pool_Fmem_siz + pool_Imem_siz + pool_LRdeltaD_siz;
    pool_Fmem_siz = pool_Imem_siz = pool_LRdeltaD_siz = 0;
  }

  heap->reset_lightweight_heap(virtperthread_size,numthreads);
    
  double *threadmeminit = 0;/* initialize threadmem[] allocation for all threads */
  if(VMEM_MINIMIZE){
    size_t siz = VMEM_MINIMIZE * 1024 * numthreads;
    POSIX_MEMALIGN(threadmeminit, PAGE, siz * sizeof(double));

    DoubleAllocations.push_back(threadmeminit);
    firstblocksiz += siz - VMEM_MINIMIZE * 1024;
  }

  long long origIJtot,IJtot;
  double LLRtotEPS2;
  if(EPS2 > 0.0){
    origIJtot = IJtot = 0;
    LLRtotEPS2 = 0.0;
  }

  if(VERB>= (HAPSTOP_DEBUG>=2 ? 1 : 3) && newLPdelta && DELTA_STOP){
    printf("\nN=%d,Y[N+1]=%0.4f,oldLP=%p: Locations of interval changes:\n",N,origY[N+1],oldLP);
    for(int I = 0; I <= N; I++){
      if(D[I] > 0){
        double mindelta = origY[N+1], maxdelta = -origY[N+1];
	for(int d = 0; d < D[I]; d++){
          double del = delta[I][d];
	  mindelta = min(del,mindelta);
	  maxdelta = max(del,maxdelta);
        }

	printf("I=%d:D[I]=%d,Y[I..I+1]=%0.4f..%0.4f,delta= %0.4f .. %0.4f (min=%0.4f,max=%0.4f)\n",I,D[I],origY[I],origY[I+1], delta[I][0], delta[I][D[I]-1], mindelta,maxdelta);
      }
    }
    fflush(stdout);
  }

  if(VERB>=3){
    printf("mprobeval : numthreads=%d, before entering multithread loop:\n",numthreads);
    fflush(stdout);
  }
  
  MFLOAT **deltaF = 0, *deltaFmem = 0;
  if(newLPdelta){
    int memsiz = FL128_lpdelta_pool_size;
    /*    for(int I = 0; I <= N; I++){
      if(D[I] <= 0)
	continue;
      memsiz += ((D[I] + VPAD) & (~VPAD));
      }
    if(DEBUG) assert(memsiz >= 0);// in case 32-bit overflow occurs
    */
    deltaF = (MFLOAT **) new MFLOAT*[N+1];
    POSIX_MEMALIGN(deltaFmem, PADDING * 2, memsiz * sizeof(MFLOAT));
    
    if(VERB>=2){
      printf("deltaFmem=%p, memsiz=%d floats, alignment=%d bytes, &deltaFmem[memsiz]= %p\n",deltaFmem,memsiz,PADDING*2, &deltaFmem[memsiz]);
      fflush(stdout);
    }

    int cnt = 0, inc;
    for(int I = 0; I <= N; I++){
      if(D[I] <= 0)
	continue;
      deltaF[I] = &deltaFmem[cnt];
      cnt += inc = (D[I] + VPAD) & (~VPAD);
      
      for(int d = 0; d < D[I]; d++)
	deltaF[I][d] = delta[I][d];
      for(int d = D[I]; d < inc; d++)
	deltaF[I][d] = 0.0f;
    }
    if(DEBUG) assert(cnt == memsiz);
  }


  int skipcnt = 0;/* how often the alignment was skipped */
  int setcnt = 0;/* how often setlimit was called */
  int qsetcnt = 0;/* how often setlimit call was skipped due to SETLIMIT_CACHE */
  int globalcnt = 0;/* how often setlimit returned a global alignment limit (IMIN==1,IMAX==N) */
  int LoopCnt = 0;

  START(mprobeval_parallel_wait);
  #pragma omp parallel num_threads(numthreads) if (/* !MDEBUG && *//* !VECTOR2_DEBUG && */ /* HERE rverb==0 && */ numthreads > 1) 
  {
    START(mprobeval_parallel);
    int myskipcnt = 0;
    int mysetcnt = 0;
    int myqsetcnt = 0;
    int myglobalcnt = 0;
    int myLoopCnt = 0;

    double wt, wt0;
    if(VERB>=2)
      wt0 = wt = wtime();
    double lastwt = 0.0;
    if(MP_MIN_TIME > 0)
      lastwt = wtime();

    /* Obtain thread id */
    int tid = 0;
    #ifdef _OPENMP
    tid = omp_get_thread_num ();
    #endif

    /* allocate all heap memory that is local to each thread */
    double *Y =  heap->alloc_double(N+2, tid);// Need a local copy in each thread since we are temporarily modifying the map (see mprobeval_add())
    double *LogScaleL =  heap->alloc_double(N+2, tid);
    double *ScaleLI = heap->alloc_double(N+2, tid);/* ScaleLI[G=1..I] = exp(LogScaleL[G]-LogScaleL[I]), where I takes on values IMIN ... IMAX */

    double *LogScaleR =  heap->alloc_double(N+2, tid);
    double *ScaleRI = heap->alloc_double(N+2, tid); /* ScaleRI[G=1..I] = exp(LogScaleR[G]-LogScaleR[I]), where I takes on values IMAX ... IMIN */
    double *JXmin = heap->alloc_double(N+1, tid);
    double *JXmax = heap->alloc_double(N+1, tid);

    for(int I = 0; I <= N+1; I++)
      Y[I] = origY[I];

    /* allocate DP array A[I=1..N][K=0..Kmax[I]][J==1..MX[]] */
    QB A;

#if USE_MIC
    const int OUT_TYPE = OUTLIER_TYPE;
    const int OUT_OLD = OUTLIER_OLD(outlierMaxBC);
#endif

#if PQ_BND >= 1
    double *LRmem = heap->alloc_double(deltaExtYRef * deltaExtXRef, tid);
    double *LRmem2 = heap->alloc_double(deltaExtYRef * deltaExtXRef, tid);
#endif

    long long myStrideMax = 0;// Used with MIN_MEM : currently Allocated memory can handle stride values up to myStrideMax
    int myKstrideMax = 0;// Used with CACHE
    MFLOAT2 *Fmem = NULL;
    int *Imem = NULL;
    size_t memsiz = 0;// current memory block size allocated at Fmem[] by mmap
    if(!(DIAGONAL && MIN_MEM)){ /* allocate Fmem[],Imem[] */
      if(sizeof(MFLOAT2) == sizeof(MFLOAT))
	Fmem = (MFLOAT2 *) heap->alloc_MFLOAT(StrideMax * mDPsizF, tid);
      else
	Fmem = (MFLOAT2 *) heap->alloc_double(StrideMax * mDPsizF, tid);
      Imem = heap->alloc_int(StrideMax * mDPsizI, tid);
      myStrideMax = StrideMax;
      myKstrideMax = KstrideMax;
#ifdef VALGRIND
      VALGRIND_MAKE_MEM_UNDEFINED(Fmem,StrideMax * mDPsizF * sizeof(MFLOAT2));
      VALGRIND_MAKE_MEM_UNDEFINED(Imem,StrideMax * mDPsizI * sizeof(int));
#endif
    }

    if(DEBUG >= 2){
      for(int I = 0; I <= N+1; I++){
	LogScaleL[I] = aNaN;
	LogScaleR[I] = aNaN;
      }
      for(int I = 0; I <= N+1; I++){
	ScaleLI[I] = aNaN;
	ScaleRI[I] = aNaN;
      }
    }

    int IMIN=1,IMAX=N;
    int *Jmin = heap->alloc_int(N+1, tid);
    int *Jmax = heap->alloc_int(N+1, tid);
    int *Imin = heap->alloc_int(maxMX+1, tid); 
    int *Imax = heap->alloc_int(maxMX+1, tid); 
    double *MSumIJ = EPS2 > 0.0 ? heap->alloc_double(maxMX*N, tid) : NULL; 
    double *MSumIJ2 = EPS2 > 0.0 ? heap->alloc_double(maxMX*N, tid) : NULL; 
    double *MSumIJ3 = EPS2 > 0.0 ? heap->alloc_double(maxMX*N, tid) : NULL; 
    double *E2sum = EPS2 > 0.0 ? heap->alloc_double(N+1, tid) : NULL;
    int *origJmin = EPS2 > 0.0 ? heap->alloc_int(N+1,tid) : NULL;
    int *origJmax = EPS2 > 0.0 ? heap->alloc_int(N+1,tid) : NULL;
    //    int origIMIN=1, origIMAX = N;

    /* statistics for EPS2 */
    int origIJcnt,IJcnt;/* count how many DP cells were used before and after applying EPS2 */
    double LLRsumEPS2;
    if(VERB_EPS2 && EPS2 > 0.0){
      origIJcnt = IJcnt = 0;
      LLRsumEPS2 = 0.0;
    }

    int IMIN1,IMAX1;/* value of IMIN,IMAX before extrapolating from the current alignment */

    int *Jmid = heap->alloc_int(N+1, tid);
    
    /* Allocate ALT[K=0..KMAX][J=1..M], ART[K=0..KMAX][J=1..M] */
    double **ALT = (double **)heap->alloc_ptr(KMAX+1, tid);// WAS91   (double **)alloca((KMAX+1)*sizeof(double *));
    double **ART = (double **)heap->alloc_ptr(KMAX+1, tid);// WAS91 (double **)alloca((KMAX+1)*sizeof(double *));
    double *PMT = (double *)heap->alloc_double(KMAX+1, tid);// WAS91 (double *)alloca((KMAX+1)*sizeof(double));
    for(int K = 0; K <= KMAX; K++){
      ALT[K] = heap->alloc_double(maxMX+1, tid);
      ART[K] = heap->alloc_double(maxMX+1, tid);
    }

#if VECTOR2_ADD
    // Allocate KtVec[t = 0..addcntImax-1], KsVec[t = 0.. addcntImax - 1]
    int *KtVec = (int *)heap->alloc_ptr(addcntImax, tid); // WAS91 (int *)alloca(addcntImax * sizeof(int));
    int *KsVec = (int *)heap->alloc_ptr(addcntImax, tid); // WAS91 (int *)alloca(addcntImax * sizeof(int));

    if(VERB>=2){
      printf("mprobeval:KtVec= %p, KsVec= %p, addcntImax=%d,tid=%d\n",KtVec,KsVec,addcntImax,tid);
      fflush(stdout);
    }
    if(DEBUG>=2 && addcntImax > 0){
      for(int t = 0; t < addcntImax; t++){
	KtVec[t] = 0;
	KsVec[t] = 0;
      }
    }

    /* Allocate ALTH[K=0..KMAX][J=1..M][H= 0..(USE_MIC ? 15 : 7)], ARTH[K][J][H], PMTH[K=0..KMAX][H=0..(USE_MIC ? 15 : 7)]  */
    double ***ALTH = (double ***)heap->alloc_ptr(KMAX+1, tid); // WAS91 (double ***)alloca((KMAX+1)*sizeof(double **));
    double ***ARTH = (double ***)heap->alloc_ptr(KMAX+1, tid); // WAS91 (double ***)alloca((KMAX+1)*sizeof(double **));
    double **PMTH = (double **)heap->alloc_ptr(KMAX+1, tid); // WAS91 (double **)alloca((KMAX+1)*sizeof(double*));
    for(int K = 0; K <= KMAX; K++){
      ALTH[K] = (double **)heap->alloc_ptr(maxMX+1, tid);
      ARTH[K] = (double **)heap->alloc_ptr(maxMX+1, tid);
    }

    for(int K = 0; K <= KMAX; K++){
      PMTH[K] = (double *)heap->alloc(VEC_SIZ*sizeof(double), tid, max(PADDING,8));
      memset(PMTH[K], 0, VEC_SIZ * sizeof(double));// avoid unitialized values in Padded regions 
      for(int J = 1; J <= maxMX; J++){
	ALTH[K][J] = (double *)heap->alloc(VEC_SIZ*sizeof(double), tid, max(PADDING,8));
	ARTH[K][J] = (double *)heap->alloc(VEC_SIZ*sizeof(double), tid, max(PADDING,8));
      }
    }
#endif  

    double *ARD = heap->alloc_double(maxMX+1, tid);
    double *ALD = heap->alloc_double(maxMX+1, tid);

    /* allocate LR arrays for current map */
    FL128 *LRdel = 0;
    if(LPdelG)
      LRdel = (FL128 *)heap->alloc((N+2)*sizeof(FL128), tid, sizeof(FL128)); /* difference array to compute LPdel[] */

    // HERE : try to use FL128b for LRadd & LRadd2

    FL128* LRadd = (FL128 *)heap->alloc((N+2)*sizeof(FL128), tid, sizeof(FL128)); /* difference array to compute LPadd[] (also used by LPdelta[]) */
    FL128* LRadd2 = (FL128 *)heap->alloc((N+2)*sizeof(FL128), tid, sizeof(FL128)); /* negative term array (also used by LPdelta[]) */
    double **LRaddT = 0;
    if(addcnt > 0){
      LRaddT = (double **)heap->alloc_ptr((N+1), tid); /* LRaddT[I][t] is contribution to LPadd[] from terms specific to site added at addlocI[I][t] */
      for(int I = 0; I <= N; I++){
	if(addcntI[I] <= 0)
	  continue;
	LRaddT[I] = (double *)heap->alloc(addcntI[I]*sizeof(double), tid, sizeof(double));
	if(DEBUG>=1+RELEASE /* HERE HERE >=2 */)
	  for(int t = 0; t < addcntI[I]; t++)
	    LRaddT[I][t] = aNaN;
      }
    }

    //    MFLOAT **deltaF = 0;
    double **LRdeltaD = 0;
    FL128 *LRdelta = (FL128 *)heap->alloc((N+2)*sizeof(FL128), tid, sizeof(FL128)); /* difference array to compute LPdelta[I] from terms that do NOT depend on specific site interval change */
    if(newLPdelta){
      //      deltaF = (MFLOAT **)heap->alloc_ptr((N+1), tid);
      LRdeltaD = (double **)heap->alloc_ptr((N+1), tid); /* LRdeltaD[I][d] is contribution to LPdelta[I][d] from terms that depend on site interval change delta[I][d] */
      for(int I = 0; I <= N; I++){
	if(D[I] <= 0)
	  continue;
	/*	deltaF[I] = (MFLOAT *)heap->alloc(((D[I] + VPAD) & (~VPAD)) * sizeof(MFLOAT), tid, (USE_MIC ? 64 : USE_AVX ? 32 : 8));
	for(int d = 0; d < D[I]; d++)
	  deltaF[I][d] = delta[I][d];
	for(int d = D[I]; d < ((D[I] + VPAD) & (~VPAD)); d++)
	deltaF[I][d] = 0.0f;*/

	/*	
	LRdeltaD[I] = (double *)heap->alloc(((D[I] + VPAD) & (~VPAD)) * sizeof(double), tid, (USE_MIC ? 64 : USE_AVX ? 32 : sizeof(double)));
	if(DEBUG>=2)
	  for(int d = 0; d < D[I]; d++)
	  LRdeltaD[I][d] = -1.0e+37; */
	/*	for(int d = D[I]; d < ((D[I] + VPAD) & (~VPAD)); d++)
		LRdeltaD[I][d] = 0.0;*/
      }
    }

    FL128 *LRsite,*LRfrag,*LRfrag2;
    if(LDEBUG){
      LRsite = new FL128[N+2];/* multiple-way LR computation by site decomposition for debugging */
      LRfrag = new FL128[N+2];/* multiple-way LR computation by frag (site interval) decomposition for debugging */
      LRfrag2 = new FL128[N+2];/* negative term array */
    }

    /* initialize perthread memory pool, only used if VMEM_MINIMIZE */
    CdoubleAllocation *doubleAlloc = NULL;
    if(VMEM_MINIMIZE){
      unsigned int threadmemsiz = VMEM_MINIMIZE * 1024;
      double *threadmem = &threadmeminit[threadmemsiz * tid];
      doubleAlloc = new CdoubleAllocation(tid, Allele, threadmem, threadmemsiz,output_prefix, CMapID ? CMapID : pcontig->id);
      
      #pragma omp critical(SwappedAllocations)
      {
	SwappedAllocations.push_back(doubleAlloc);
      }
    }

    START(mprobeval_for_wait);
    #pragma omp for schedule(dynamic,1)
    for(int mt = 0; mt < MD; mt++){
      START(mprobeval_for);
      if(OMP_DEBUG) myLoopCnt++;

      int m = (MPROBEVAL_MAPORDER && SETLIMIT_CACHE && maporder) ? maporder[mt] : mt;
      int M = MX[m];
      double *Xm = X[m];
      FL128 TBwtm = TBmapWT ? TBmapWT[m] : 1.0;

      CdoubleBlock **ppBlock = pBlock ? &pBlock[m] : 0;

      int myaddcnt = addcnt;/* will either equal addcnt OR 0 (to skip LPaddT[] and addloc[] code for current molecule m) */

      if(VERB/* HERE HERE >=3 */ && MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M){
        if(oldLP)
          printf("\t mt=%d/%d, m=%d(id=%lld): M=%d,Ymin= %0.4f, Ymax= %0.4f, oldLP[m]= %0.8e, TBwtm= %0.8e, skipmap[m]= %d: Starting mprobeval: wall=%0.6f\n",
		 mt,MD,m,gmap[pcontig->contig[m].mapid]->id,M,Ymin,Ymax,oldLP[m],(double)TBwtm, skipmap ? skipmap[m] : 0, wtime());
        else
          printf("\t mt=%d/%d, m=%d(id=%lld): M=%d: Starting mprobeval: wall=%0.6f\n",mt,MD,m,gmap[pcontig->contig[m].mapid]->id,M,wtime());

        fflush(stdout);
      }

      /* Limit DP to near the previous alignment defined by Xm[J] to Y[I=map[MD][map[m][J]]] (I= -1 if Xm[J] is not aligned) */
      int noskip = -1;

      if(oldLP || TBwtm <= 0.0){
        double LPm = oldLP ? oldLP[m] : 0.0;

	noskip = ((skipmap && skipmap[m]) ? 0 : 1);
	if(noskip && TBwtm > 0.0){
	  Csetlimit *plimit = &limit[m];
          if(SETLIMIT_CACHE && Ymin >= 0.0 && plimit->ileft >= 0){
	    double Yleft = Hcuts[plimit->ileft] - plimit->Xleft;
	    double Yright = Hcuts[plimit->iright] + plimit->Xright;
	    if(max(Yleft,Ymin) > min(Yright,Ymax) + 1e-6){
	      myqsetcnt++;
	      noskip = 0;
	    }
          }
	  if(SETLIMIT_CACHE && plimit->ileft >= 0 && (newLPdelta || (SKIPMAP_ADD /* WAS && addcnt > 0 */ && Ymin < 0.0))){// quicker check to see if region Y[IMIN-LRANGE..IMAX+LRANGE] has D[I] > 0 or addloc[0..addcnt-1] 
	    int Imin = plimit->Ileft;
	    int Imax = plimit->Iright;
	    if(!(Imax <= N+1 && nmapMD[Imin] == plimit->ileft && nmapMD[Imax] == plimit->iright))
	      plimit->ileft = -1;// need to recompute plimit->Ileft & plimit->Iright since labels have been added/deleted 
	    else {
	      int overlap = newLPdelta ? 1 : (SKIPMAP_ADD /* WAS && addcnt > 0 */ && Ymin < 0.0) ? 1 : 0;
	      if(newLPdelta){
		int I = Imin;
		for(; I < Imax; I++)
		  if(D[I] > 0)
		    break;
		if(!(I < Imax))
		  overlap--;
                if(VERB>=2 && DELTA_STOP && (rverb || (MDEBUG && m==MDEBUG_M))){
                  #pragma omp critical
		  {
		    if(overlap)
		      printf("m=%d:Imin=%d,IMIN=%d,IMIN1=%d,IMAX1=%d,IMAX=%d,Imax=%d,Y[Imin,Imax]=%0.4f..%0.4f,I=%d,D[I]=%d: molecule not skipped\n",
			     m,Imin,IMIN,IMIN1,IMAX1,IMAX,Imax,Y[Imin],Y[Imax],I,D[I]);
		    else
		      printf("m=%d:Imin=%d,IMIN=%d,IMIN1=%d,IMAX1=%d,IMAX=%d,Imax=%d,Y[Imin,Imax]=%0.4f..%0.4f,I=%d: molecule skipped\n",
			     m,Imin,IMIN,IMIN1,IMAX1,IMAX,Imax,Y[Imin],Y[Imax],I);
		    fflush(stdout);
		  }
		}
	      } else if(SKIPMAP_ADD /* WAS && addcnt > 0 */ && Ymin < 0.0){/* check if any of the values addcntI[Imin .. Imax-1] are > 0 */
		if(addcnt <= 0 || addcntIcum[Imax] <= addcntIcum[Imin]){
		  if(SKIPMAP_ADD>=2)
		    overlap--;
		  else {
		    if(newLPa){
		      TminM[m] = TmaxM[m] = 0;
		      if(VMEM_MINIMIZE) newLPa[m] = thread_alloc((size_t)1, doubleAlloc, 0.0, ppBlock,tid);// threadmem, threadmemsiz);
		      newLPa[m][0] = LPm;
		    } else
		      LPaddI += LPm;
		    myaddcnt = 0;
		  }
		}
                if((VERB>=2 && (rverb || (MDEBUG && m==MDEBUG_M))) || (DEBUG && !(0 <= Imin && Imin <= Imax && Imax <= N+1))){
		  #pragma omp critical
		  {
		    printf("m=%d: Imin=%d,Imax=%d,N=%d,overlap=%d:addcntIcum[Imin]=%d,addcntIcum[Imax]=%d,addcnt=%d, SKIPMAP_ADD=%d\n",
			   m,Imin,Imax,N,overlap,addcntIcum[Imin],addcntIcum[Imax], addcnt, SKIPMAP_ADD);
		    fflush(stdout);
		    if(DEBUG) assert(0 <= Imin && Imin <= Imax && Imax <= N+1);
		  }
		}
	      }

	      if(overlap <= 0){
                if(VERB/* HERE >=2 */ && MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M){
		  #pragma omp critical
		  {
		    printf("m=%d:overlap=%d:SKIPMAP_ADD=%d,addcnt=%d,Ymin=%0.3f,Ymax=%0.3f\n",m,overlap,SKIPMAP_ADD,addcnt,Ymin,Ymax);
		    fflush(stdout); exit(1);
		  }
		}
		myqsetcnt++;
		noskip = 0;
		if(skipmap)
		  skipmap[m] = 1;
	      }
	    }
	  }

          /* uses cache limit[m] to hold last estimate of Y[IMIN-LRANGE],Y[IMAX+LRANGE] that is reset only when map[m],mapK[m] are updated with UpdateMap  */
	  if(noskip){
	    mysetcnt++;
	    noskip = setlimit(m,n,Hcuts,N,Y,M,Xm,map[m],mapK[m],&limit[m],map[MD],nmapMD,Ymin,Ymax,IMIN,IMAX,Jmin,Jmax,Imin,Imax,JXmin,JXmax,Jmid,IMIN1,IMAX1,pcontig, 0);
            
            if(VERB/* HERE >=2 */ && MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M){
              printf("m=%d: setlimit() returned noskip = %d\n",m,noskip);
              fflush(stdout);
            }
	  }
	  /* If noskip==0 : 
	     Xm does not overlap part of Y[] (Ymin..Ymax) that changed (or WILL change due to addloc[] or delta[]) since oldLP[] was computed. Specifically (if LRANGE_FIX >=2) :
	     The range Ymin..Ymax does not overlap Y[IMIN-LRANGE..IMAX+LRANGE] (Y[IMIN..IMAX] is the maximum alignment range of Xm) */

	}
	if(noskip && TBwtm > 0.0 && (newLPdelta || (SKIPMAP_ADD && addcnt > 0 && Ymin < 0.0))){	/* Also check if Xm overlaps any region of Y[I..I+1] with D[I] > 0 or addloc[0..addcnt-1] (ignores LPdelG)*/
          int Imin = max(0,IMIN - LRANGE);
	  int Imax = min(N+1, IMAX + LRANGE);

	  int overlap = (newLPdelta ? 1 : 0) + ((SKIPMAP_ADD && addcnt > 0 && Ymin < 0.0) ? 1 : 0);
	  if(newLPdelta){
            int I = Imin;
            for(; I < Imax; I++)
	      if(D[I] > 0)
		break;
	    if(!(I < Imax))
	      overlap--;
	    if(VERB>=2 && DELTA_STOP && rverb){
              #pragma omp critical
	      {
  	        if(overlap)
		  printf("m=%d:Imin=%d,IMIN=%d,IMIN1=%d,IMAX1=%d,IMAX=%d,Imax=%d,Y[Imin,Imax]=%0.4f..%0.4f,I=%d,D[I]=%d: molecule not skipped\n",
		       m,Imin,IMIN,IMIN1,IMAX1,IMAX,Imax,Y[Imin],Y[Imax],I,D[I]);
		else
		  printf("m=%d:Imin=%d,IMIN=%d,IMIN1=%d,IMAX1=%d,IMAX=%d,Imax=%d,Y[Imin,Imax]=%0.4f..%0.4f,I=%d: molecule skipped\n",
		       m,Imin,IMIN,IMIN1,IMAX1,IMAX,Imax,Y[Imin],Y[Imax],I);
		fflush(stdout);
	      }
	    }
          }
	  if(SKIPMAP_ADD && addcnt > 0 && Ymin < 0.0){/* check if any of the values addcntI[Imin .. Imax-1] are > 0 */
	    if(addcntIcum[Imax] <= addcntIcum[Imin]){
	      if(SKIPMAP_ADD>=2)
		overlap--;
	      else {
		if(newLPa){
		  TminM[m] = TmaxM[m] = 0;
		  if(VMEM_MINIMIZE) newLPa[m] = thread_alloc((size_t)1, doubleAlloc, 0.0, ppBlock,tid);// threadmem, threadmemsiz);
		  newLPa[m][0] = LPm;
		} else
		  LPaddI += LPm;
		skipcnt++;// partial skip, but with SKIPMAP_ADD==1, full skip is not possible
		myaddcnt = 0;// disables LPaddT[] and addloc[] related code for current molecule m
	      }
	    }
	    if(VERB>=2 && rverb){
              #pragma omp critical
	      {
		printf("m=%d: IMIN=%d,IMAX=%d,Imin=%d,Imax=%d,N=%d,overlap=%d:addcntIcum[Imin]=%d,addcntIcum[Imax]=%d,addcnt=%d,myaddcnt=%d\n",
		       m,IMIN,IMAX,Imin,Imax,N,overlap,addcntIcum[Imin],addcntIcum[Imax], addcnt,myaddcnt);
		fflush(stdout);
	      }
	    }
	  }
	  if(overlap <= 0){
            noskip = 0;
	    if(skipmap)
	      skipmap[m] = 1;
	  }
        }

	if(!noskip || TBwtm <= 0.0){
          newLP[m] = LPm;// NEW2
	  if(DEBUG && !isfinite(newLP[m])){
            printf("m=%d:newLP[m]= %0.8e (oldLP[m]= %0.8e)\n",m,newLP[m],oldLP ? oldLP[m] : 0.0);
	    fflush(stdout);
	    assert(isfinite(newLP[m]));
          }      

	  if(VERB_EPS2 && EPS2 > 0)
	    LLRsumEPS2 += LPm * TBwtm;

	  if(myaddcnt > 0 && newLPa){
	    TminM[m] = TmaxM[m] = 0;
	    if(VMEM_MINIMIZE) newLPa[m] = thread_alloc((size_t)1, doubleAlloc, 0.0, ppBlock,tid);
	    newLPa[m][0] = LPm;
	  }
	  if(newLPd){
	    DminM[m] = 0 /* WAS8 1 */;
	    DmaxM[m] = 0 /* WAS 1 */;
	    if(VMEM_MINIMIZE){/* NEW8 */
	      newLPd[m] = thread_alloc((size_t)1, doubleAlloc, 0.0, ppBlock, tid);
	      newLPd0[m] = LPm;
	    }
	    newLPd[m][0] = LPm;
	  }
	  if(newLPdelta){
	    AminM[m] = 1;
	    AmaxM[m] = 0;
	  }
	  
          if(TDEBUG || MDEBUG || !newLPd || (myaddcnt > 0 && !newLPa) || (VERB && (VERB>=2 || rverb || (MDEBUG_M >= 0 && m==MDEBUG_M && MTRACE)))){
            #pragma omp critical
            {// NOTE : this critical section should not normally be triggered, except for debugging
	      LLRsum += LPm *= TBwtm;
	      if(LPdelG)
		LPdel[0] += LPm;
	      if(myaddcnt > 0)
		LPaddI += LPm;
	      /* LPdelta[] array is difference from LLRsum, so no update required for map m, if values are the same */

              if(VERB && (VERB>=2 || rverb>=2 || (MDEBUG_M >= 0 && m==MDEBUG_M && (MTRACE || MDEBUG)))){
		printf("m=%d/%d:mapid=%d(id=%lld):LLR=%0.10f(delta=%0.10f from oldLP[]):tid=%d:wall=%0.6f\n",
		       m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,LLRsum,(double)LPm,tid,wtime());
		fflush(stdout);
	      }
	    } // pragma omp critical
	  }

	  if(VERB) myskipcnt++;

          if(VERB/* HERE >=2 */ && MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M){
printf("\t m=%d/%d: Skipping mprobeval, noskip=%d, TBwtm= %0.8e: wall=%0.6f\n",m,MD,noskip, (double)TBwtm, wtime());
            fflush(stdout);
          }
	  STOP(mprobeval_for);
	  continue;
        }
      } else {
        (void) setlimit(m,n,Hcuts,N,Y,M,Xm,map[m],mapK[m],&limit[m],map[MD],nmapMD,-1.0,-1.0,IMIN,IMAX,Jmin,Jmax,Imin,Imax,JXmin,JXmax,Jmid,IMIN1,IMAX1,pcontig,1);/* sets IMIN,IMAX,Jmin[IMIN..IMAX],Jmax[IMIN..IMAX] */
	mysetcnt++;
      }

      if(DEBUG && globalfallback && !(IMIN <= IMIN1 && IMIN1 <= IMAX1 && IMAX1 <= IMAX)){
	#pragma omp critical
	{
	  printf("\n m=%d/%d:mapid=%d(id=%lld),oldLP=%p:IMIN=%d,IMIN1=%d,IMAX1=%d,IMAX=%d,N=%d,LRANGE=%d,noskip=%d\n",
		 m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,oldLP,IMIN,IMIN1,IMAX1,IMAX,N,LRANGE,noskip);
	  fflush(stdout);
	  assert(IMIN <= IMIN1 && IMIN1 <= IMAX1 && IMAX1 <= IMAX);
	}
      }

      myglobalcnt += (N > 0 && IMIN==1 && IMAX==N) ? 1 : 0;
      if(VERB>=3 && noskip && (N > 0 && IMIN==1 && IMAX==N)){
	printf("\t m=%d/%d:setlimit() returned global alignment, IMIN=%d,IMAX=%d,N=%d!\n",m,MD,IMIN,IMAX,N);
	fflush(stdout);
	assert(!(noskip && N > 0 && IMIN==1 && IMAX==N));
      }

      if(IMIN > IMAX){/* setlimit found no valid previous alignment and globafallback==0 */
	double LPm = newLP[m] = log(max(1e-300,LRbias) + PoutlierEnd * PoutlierEnd);
	if(DEBUG && !isfinite(newLP[m])){
	  printf("m=%d:newLP[m]= %0.8e (oldLP[m]= %0.8e)\n",m,newLP[m],oldLP ? oldLP[m] : 0.0);
	  fflush(stdout);
	  assert(isfinite(newLP[m]));
	}      
	
	if(VERB_EPS2 && EPS2 > 0)
	  LLRsumEPS2 += LPm * TBwtm;
	
	if(newLPa && myaddcnt > 0){
	  TminM[m] = TmaxM[m] = 0;
	  if(VMEM_MINIMIZE) newLPa[m] = thread_alloc((size_t)1, doubleAlloc, 0.0, ppBlock,tid);// threadmem, threadmemsiz);
	  newLPa[m][0] = LPm;
	}
	if(newLPd){
	  DminM[m] = 0 /* WAS8 1 */;
	  DmaxM[m] = 0 /* WAS 1 */;
	  if(VMEM_MINIMIZE){/* NEW8 */
	    newLPd[m] = thread_alloc((size_t)1, doubleAlloc, 0.0, ppBlock, tid);
	    newLPd0[m] = LPm;
	  }
	  newLPd[m][0] = LPm;
	}
	if(newLPdelta){
	  AminM[m] = 1;
	  AmaxM[m] = 0;
	}

	if(TDEBUG || MDEBUG || !newLPd || (myaddcnt > 0 && !newLPa) || (VERB && (VERB>=2 || rverb || (MDEBUG_M >= 0 && m==MDEBUG_M && MTRACE)))){
          #pragma omp critical
	  {// NOTE : this critical section is only triggered for debugging
	    LLRsum += LPm *= TBwtm;
	    if(LPdelG)
	      LPdel[0] += LPm;
	    if(myaddcnt > 0)
	      LPaddI += LPm;
	    /* LPdelta[] array is difference from LLRsum, so no update required for map m, if values are the same */

            if(VERB && (VERB>=2 || rverb>=2 || ((MTRACE || MDEBUG) && MDEBUG_M >= 0 && m==MDEBUG_M))){
	      printf("m=%d/%d:mapid=%d(id=%lld):LLR=%0.10f(delta=%0.10f from oldLP[])\n",
		     m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,LLRsum,(double)LPm);
	      fflush(stdout);
	    }
	  } // pragma omp critical
	}

	if(VERB/* HERE >=2 */) myskipcnt++;

        if(VERB/* HERE >=2 */ && MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M){
          printf("\t m=%d/%d: Skipping mprobeval(B) : wall=%0.6f\n",m,MD,wtime());
          fflush(stdout);
        }

	STOP(mprobeval_for);
	continue;
      }

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244)){
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
	      printf("\t\t\t J=%d:map[m][J]=i=%d, mapK[m][J]=k= %d, n=%d : map[MD][i]=I=%d, map[MD][i-k]=K=%d: Xm[J]=%0.3f: Y[K,I]= %0.3f,%0.3f\n", J,i,k,n,I,K,Xm[J],Y[K],Y[I]);
	    }
	    for(int I = IMIN; I <= IMAX; I++)
	      printf("\t\t\t I=%d:Kmax[I]=%d,Jmin[I]=%d,Jmax[I]=%d:Y[I]= %0.3f, Xm[Jmin[I]]= %0.3f, Xm[Jmax[I]]= %0.3f\n",I,Kmax[I],Jmin[I],Jmax[I],Y[I],Xm[Jmin[I]], Xm[Jmax[I]]);
	  }
	  fflush(stdout);
        }
      }

      if(VERB && rverb && MDEBUG_M >= 0 && m==MDEBUG_M){
	printf("m=%d/%d:global=%d,IMIN=%d,IMAX=%d:\n",m,MD,global,IMIN,IMAX);
	for(int I=IMIN; I <= IMAX; I++)
	  printf("  I=%d:Kmax[I]=%d,Jmin[I]=%d,Jmax[I]=%d\n",I,Kmax[I],Jmin[I],Jmax[I]);
	fflush(stdout);
      }

      if(DEBUG>=2){
	assert(1 <= IMIN && IMIN <= IMAX && IMAX <= N);
	for(int I = IMIN; I <= IMAX; I++){
          if(!(1 <= Jmin[I] && Jmin[I] <= Jmax[I] && Jmax[I] <= M)){
            #pragma omp critical
            {
              printf("\nm=%d,M=%d,N=%d,IMIN=%d,IMAX=%d:I=%d,Jmin[I]=%d,Jmax[I]=%d\n",m,M,N,IMIN,IMAX,I,Jmin[I],Jmax[I]);
              fflush(stdout);
	      assert(1 <= Jmin[I] && Jmin[I] <= Jmax[I] && Jmax[I] <= M);
            }
          }
        }
	for(int J = Jmin[IMIN]; J <= Jmax[IMAX]; J++){
	  if(DEBUG && Imin[J] <= Imax[J] && !(1 <= Imin[J] && Imax[J] <= N)){
	    #pragma omp critical
	    {
	      printf("\nm=%d,M=%d,N=%d,IMIN=%d,IMAX=%d\n",m,M,N,IMIN,IMAX);
	      for(int I = IMIN; I <= IMAX; I++)
		printf("I=%d:Jmin[I]=%d,Jmax[I]=%d\n",I,Jmin[I],Jmax[I]);
	      for(int j = 1; j <= M; j++)
		printf("J=%d:Imin[J]=%d,Imax[J]=%d\n",j,Imin[j],Imax[j]);
	      fflush(stdout);
	      assert(1 <= Imin[J] && Imin[J] <= Imax[J] && Imax[J] <= N);
	    }
	  }
	}
      }

      /* set up parameters from setlimit(), used when extrapolating alignment limits for new sites */
      double Ltheta = LOCAL_RANGE ? max(1.0,Xm[M]-Xm[1])/max(1,M-1) : Xtheta;
      double Xrange = RANGE * Ltheta;/* used when FIX_RANGE_Y && RANGE_Y==0 as maximum alignment adjustment in kb */
      double Xrange2 = Xrange; /* used when FIX_RANGE_Y && RANGE_Y==0 as maximum alignment adjustment in kb during alignment extrapolation : should reflect both RANGE and max(DELTA_X,deltaExtXRef) */
      if(RefineXRange)
	Xrange2 += max(DELTA_X,deltaExtXRef) * Ltheta;

      /* re-allocate memory for A.field(I= IMIN..IMAX, K= 0..Kmax[I], J= Jmin[I]..Jmax[I]) */
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
	  printf("ERROR in mprobeval.cpp: IMIN=%d,IMAX=%d,N=%d,Jmin1=%d,Jmax1=%d,deltaJ=%d,M=%d,Istride=%d : Istride * IMAX = %llu is too large to fit in 32 bit int\n",
		 IMIN,IMAX,N,Jmin1,Jmax1,deltaJ,M,A.Istride, ((long long)A.Istride) * ((long long)IMAX));
	  fflush(stdout);exit(1);
	}

        if(DEBUG && sizeof(A.Kstride) <= 4 && ((long long)(IMAX-IMIN + 1LL)) * ((long long)deltaJ) > (long long)MASK_LL(31)){
	  printf("ERROR in mprobeval.cpp: IMIN=%d,IMAX=%d,N=%d,Jmin1=%d,Jmax1=%d,deltaJ=%d,M=%d : (IMAX-IMIN+1)*deltaJ=%llu is too large to fit in 32 bit int\n",
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
	  A.Kstride = (A.Kstride + ALIGN - 1) & (~(ALIGN - 1));
	stride = (1LL + (CACHE ? kmax : A.kmax)) * A.Kstride;
	if(CACHE > 0&& ALIGN2 > 0)
	  stride = (stride + ALIGN2 - 1) & (~(ALIGN2 - 1));
	firstindex = IMIN * A.Istride + Jmin1;

	if(DEBUG && CACHE > 0 && (ALIGN > 0 || ALIGN2 > 0) && stride <= myStrideMax && !(A.Kstride <= myKstrideMax)){
	  #pragma omp critical
	  {
	    printf("m=%d:IMIN=%d,IMAX=%d,N=%d,deltaJ=%d:A.kmax=%d(kmax=%d),A.Kstride=%d(myKstrideMax=%d),stride=%lld(myStrideMax=%lld)\n",
		   m,IMIN,IMAX,N,deltaJ,A.kmax,kmax,A.Kstride,myKstrideMax,stride,myStrideMax);
	    fflush(stdout);
	    assert(A.Kstride <= myKstrideMax);
	  }
	}

	int MemBins = max(2,MIN_MEM);/* Allow for some extra memory of 1 part per MemBins and reduce memory when need drops below (MemBins-1)/(MemBins+1) (unless MP_MIN_TIME && maporder) */

	if(MIN_MEM && (stride > myStrideMax || MIN_MEM >= 8 || (MIN_MEM>=2 && MP_MIN_TIME <= 0 && stride * (MemBins+1) < myStrideMax * (MemBins-1)))){
	  /* minimize virtual memory and swapping by reallocating Fmem[],Imem[] using either free()/malloc() or mmap()/madvise(DONT_NEED) */
#if USE_MMAP
	  if(Fmem) free_huge_pages(Fmem);
#else
	  //	  if(Fmem) free(Fmem);
	  if(Fmem && munmap(Fmem, memsiz)){
	    int eno = errno;
	    char *err = strerror(eno);
	    printf("tid=%d:munmap(%p,%lu) failed: errno=%d:%s\n",tid,Fmem,memsiz,eno,err);
	    dumpmemmap();
	    fflush(stdout);exit(1);
	  }
#endif

	  if(DEBUG) assert(stride <= StrideMax);
	  if(DEBUG && CACHE) assert(A.Kstride <= KstrideMax);

          myKstrideMax = min(A.Kstride * ((MIN_MEM >= 8 /* || MP_MIN_TIME */) ? MemBins : MemBins+1) / MemBins, KstrideMax);
	  if(CACHE > 0 && ALIGN > 0){
	    int origmyKstrideMax = myKstrideMax;
	    myKstrideMax = (myKstrideMax + ALIGN - 1) & (~(ALIGN - 1));
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

	  memsiz = myStrideMax * (mDPsizF * sizeof(MFLOAT2) + mDPsizI*sizeof(int)) + PADDING;
	  if(newLPdelta)
	    memsiz += FL128_lpdelta_pool_size * sizeof(double);

	  if(VERB>=3 && newLPdelta && !origStrideMax /* && LPdelta */){
	    #pragma omp critical
	    {
	      double InvMB = 1.0/(1024.0*1024.0);
	      printf("tid=%d:m=%d/%d: Allocating %0.2f Mb for Fmem[],Imem[] & LRdeltaD[I][]:wtime()= %0.6f\n",tid,m,MD,memsiz*InvMB,wtime());
	      fflush(stdout);
	    }
	  }

	  char *mem = 0;
	  size_t cnt = 0;

#if USE_MMAP
	  mem = (char *) malloc_huge_pages(memsiz);
#else
	  if(CACHE > 0 && (ALIGN > 0 || ALIGN2 > 0 || PADDING > 0)){
	    mem = (char *) mmap(NULL, memsiz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	    if(mem == MAP_FAILED){
	      int eno = errno;
	      char *err = strerror(eno);
	      #pragma omp critical
	      {
		printf("mprobeval:tid=%d: mmap of %lu bytes failed: errno=%d: %s\n",tid,memsiz,eno,err);
		dumpmemmap();
		fflush(stdout); exit(1);
	      }
	    }
	  } else {
	    mem = (char *)malloc(memsiz);
	    if(!mem){
	      #pragma omp critical
	      {
		printf("mprobeval:tid=%d:malloc(%llu) failed : stride=%lld,myStrideMax=%lld(StrideMax=%lld),mDPsizF=%lld,mDPsizI=%lld\n",
		       tid,(unsigned long long)memsiz,  stride,myStrideMax,StrideMax,mDPsizF,mDPsizI);
		fflush(stdout); exit(1);
	      }
	    }
	  }
#endif

	  lastwt = wtime();

#ifdef VALGRIND
	  VALGRIND_MAKE_MEM_UNDEFINED(Fmem, memsiz);
#endif

	  Fmem = (MFLOAT2 *)&mem[cnt];
	  cnt += mDPsizF * myStrideMax * sizeof(MFLOAT2);

	  Imem = (int *) &mem[cnt];
	  cnt += mDPsizI * myStrideMax * sizeof(int);

	  if(newLPdelta){
	    if(PADDING > 0)
	      cnt = (cnt + (PADDING-1)) & (~(PADDING-1));
	    for(int I = 0; I <= N; I++){
	      if(D[I] <= 0)
		continue;

	      LRdeltaD[I] = (double *) &mem[cnt];
	      cnt += ((D[I] + VPAD) & (~VPAD)) * sizeof(double);

	      if(DEBUG>=2)
		for(int d = 0; d < D[I]; d++)
		  LRdeltaD[I][d] = -1.0e+37;
	    }
	  }
	  if(DEBUG) assert(cnt <= memsiz);
        } // Fmem[],Imem[],LRdeltaD[] reallocated
	else if(MP_MIN_TIME > 0 /* WAS MPROBEVAL_MAPORDER && maporder */){
	  double wt = wtime();
	  if(MP_MIN_TIME > 0 && wt - lastwt >= MP_MIN_TIME){
	    if(madvise(Fmem, memsiz, MADV_DONTNEED)){
	      int eno = errno;
	      char *err = strerror(eno);
	      
	      #pragma omp critical
	      {
		printf("tid=%d,lastwt=%0.6f,wt=%0.6f: madvise(%p,%lu,MADV_DONTNEED) failed: errno=%d: %s\n", tid, lastwt, wt, Fmem, memsiz,eno,err);
		fflush(stdout);exit(1);
	      }
	    }
	    lastwt = wt;
	  }
	}

	if(CACHE){
	  if(DEBUG) assert(A.Kstride <= myKstrideMax);
	  if(DEBUG) assert(stride <= myStrideMax);
	  if(DEBUG) assert(A.kmax <= kmax);
	  A.Kstride = myKstrideMax;
	  stride = myStrideMax;
	  if(DEBUG && !ALIGN2) assert(stride == (1LL + (CACHE ? kmax : A.kmax)) * A.Kstride);
	  if(DEBUG && ALIGN2 > 0) assert(stride >= (1LL + (CACHE ? kmax : A.kmax)) * A.Kstride);
	}
      } else { // DIAGONAL==0
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

      long long mcnt = 0;
      if(DEBUG) assert(stride <= myStrideMax);
#define FBLOCK(field) { A.field##_base = &Fmem[mcnt - firstindex]; mcnt += stride;}
      FBLOCK(UL);
      FBLOCK(UR);
      FBLOCK(bUL);
      FBLOCK(bUR);
      FBLOCK(AL);
      FBLOCK(AR);
#undef FBLOCK
      if(DEBUG && !(mcnt <= stride * mDPsizF)){
        #pragma omp critical
        {
          printf("tid=%d,m=%d/%d:mcnt=%lld,stride=%lld,mDPsizF=%lld\n",tid,m,MD,mcnt,stride,mDPsizF);
	  fflush(stdout);
          assert(mcnt <= stride * mDPsizF);
        }
      }
      
      mcnt = 0;
#define BLOCK(field) { A.field##_base = &Imem[mcnt - firstindex]; mcnt += stride;}
#if PQ_BND
      BLOCK(Pmin);
      BLOCK(Pmax);
      BLOCK(Qmin);
      BLOCK(Qmax);
#endif
      BLOCK(ULmin);
      BLOCK(ULmax);
      BLOCK(URmin);
      BLOCK(URmax);
      BLOCK(Lij);
      BLOCK(Rij);
#undef BLOCK
      if(DEBUG) assert(mcnt <= stride * mDPsizI);

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) /* && LPdelta*/){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After memory allocation:IMIN=%d,IMAX=%d,A.kmax=%d,A.Kstride=%d,A.Istride=%d:wall time= %0.6f\n", m,tid,IMIN,IMAX,A.kmax,A.Kstride,A.Istride,wtime() - wt);
	  fflush(stdout);
        }
      }

      /* initialize A[I=IMIN..IMAX][K=0..Kmax[I]][J=Jmin[I]..Jmax[I]].UL,UR, Lij,Rij  based on outlier termination */
      MFLOAT PoutlierEndf = PoutlierEnd;
      for(int I= IMIN; I <= IMAX; I++){
	for(int K = 0; K <= Kmax[I]; K++){
	  int J = Jmin[I];
	  for(; J <= Jmax[I]; J++){
	    A.UL(I,K,J) = A.bUL(I,K,J) = PoutlierEndf; /* should never be the best score if PoutlierEnd == 0 */
	    if((MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M /* && p->LR+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)
	       || (VERB>=3 && rverb && m==MDEBUG_M && I==33 && K==0 && J==4 && PoutlierEndf > 0.0)){
	      printf("m=%d:I=%d,K=%d,J=%d:UR=UL=%0.10e(PMik=%0.10e,PoutlierEnd=%0.10e)\n",m,I,K,J,A.UL(I,K,J),PRtab[K][I].PM,PoutlierEnd);
	      fflush(stdout);
	    }
	    A.Lij(I,K,J) = I-K;
	    A.Rij(I,K,J) = I;
	    A.UR(I,K,J) = A.bUR(I,K,J) = PoutlierEndf;
	    if(DEBUG>=2){/* If these values are ever used it is a bug */
	      A.AL(I,K,J) = A.AR(I,K,J) = -1.0;
	      A.ULmin(I,K,J) = A.ULmax(I,K,J) = A.URmin(I,K,J) = A.URmax(I,K,J) = -1;
	    }
	  }
	  if(DEBUG>=2 && !DIAGONAL){
	    for(J = Jmax[I]+1; J <= M; J++){/* if these values are ever used it is a bug */
	      A.UL(I,K,J) = A.bUL(I,K,J) = aNaN;
	      A.UR(I,K,J) = A.bUR(I,K,J) = aNaN;
	      A.Lij(I,K,J) = A.Rij(I,K,J) = -1;
	      A.AL(I,K,J) = A.AR(I,K,J) = aNaN;
	      A.ULmin(I,K,J) = A.ULmax(I,K,J) = A.URmin(I,K,J) = A.URmax(I,K,J) = -1;
	    }
	    for(J = 1; J < Jmin[I]; J++){/* if these values are ever used it is a bug */
	      A.UL(I,K,J) = A.bUL(I,K,J) = aNaN;
	      A.UR(I,K,J) = A.bUR(I,K,J) = aNaN;
	      A.Lij(I,K,J) = A.Rij(I,K,J) = -1;
	      A.AL(I,K,J) = A.AR(I,K,J) = aNaN;
	      A.ULmin(I,K,J) = A.ULmax(I,K,J) = A.URmin(I,K,J) = A.URmax(I,K,J) = -1;
	    }
	  }
	}
      }
    
      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After initializing A.UL,UR,Lij,Rij based on endoutliers:wall time= %0.6f\n", m,tid,wtime() - wt);
	  fflush(stdout);
        }
      }

      /* initialize A[I = Imin[J]..Imax[J]][K=0..Kmax[I]][J=1..M].Lij as leftmost site in Y[] that overlaps Xm with Y[I-K..I] aligned with Xm[J] */
      /* Initialize corresponding UL score with left ends unaligned */
      /* only compute cases where J <= DELTA_XU (OR EXTEND && I-K <= DELTA_YU && J <= DELTA_XU+EXTEND_MAX) */
      int jmax = EXTEND ? min(M,DELTA_XU+EXTEND_MAX) : min(M,DELTA_XU);
      if((0 && MDEBUG_I>=0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M)||
	 (VERB>=3 && m==19 && addcnt==279)){
	printf("m=%d:jmax=%d,M=%d\n",m,jmax,M);
	fflush(stdout);
      }
      for(int J = 1; J <= jmax; J++){
	int Lij = 0;
	double IvarX;
	double FLj = FLJ(J,Xm[J],IvarX);
	int ImKmax = (!EXTEND || J <= DELTA_XU) ? Imax[J] : min(DELTA_YU + 1/*for deleted site*/,Imax[J]);
	if((0 && MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && J== MDEBUG_J) ||
	   (VERB>=3 && m==94 && newLPdelta && J==6)){
	  printf("m=%d,J=%d,IMIN=%d,IMAX=%d:ImKmax=%d,Imin[J]=%d,Imax[J]=%d\n",m,J,IMIN,IMAX,ImKmax,Imin[J],Imax[J]);
	  fflush(stdout);
	}
	for(int ImK = max(1,Imin[J] - A.kmax); ImK <= ImKmax; ImK++){
	  int Lijk = Lij;
          if(DEBUG>=2)assert(Lijk <= ImK);
	  int I;
	  if((0 && MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && J==MDEBUG_J && ImK == (MDEBUG_I-MDEBUG_K)) ||
	     (VERB>=3 && m==94 && newLPdelta && ImK == 106 && J==6)){
	    printf(" m=%d,J=%d,ImK=%d,IMIN=%d,IMAX=%d:Imin[J]=%d,Imax[J]=%d,Jmin[I-K]=%d,Lij=%d,Lijk=%d\n",m,J,ImK,IMIN,IMAX,Imin[J],Imax[J],Jmin[ImK],Lij,Lijk);
	    fflush(stdout);
	  }
	  for(int K = max(0,Imin[J]-ImK); (I = ImK+K) <= Imax[J] && K <= Kmax[I]; K++){
	    if((0 && MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && J==MDEBUG_J && I == MDEBUG_I && K == MDEBUG_K)||
	       (VERB>=3 && m==94 && newLPdelta && ImK == 106 && J==6)){
	      printf("m=%d,J=%d,I=%d,K=%d,IMIN=%d,IMAX=%d:Jmin[I]=%d,Jmax[I]=%d,Kmax[I]=%d,Lij=%d,Lijk=%d\n",m,J,I,K,IMIN,IMAX,Jmin[I],Jmax[I],Kmax[I],Lij,Lijk);
	      fflush(stdout);
	    }
	    if(DEBUG>=2) assert(IMIN <= I && I <= IMAX);
	    if(DEBUG>=2) assert(Jmin[I] <= J && J <= Jmax[I]);
	    if(DEBUG>=2) assert(0 <= K && K <= Kmax[I] && Kmax[I] <= A.kmax);
	    if(DEBUG>=2) assert(J <= DELTA_XU || (EXTEND && I-K <= DELTA_YU + 1/*for deleted site*/ && J <= DELTA_XU + EXTEND_MAX));

	    if(DEBUG>=2) assert(I >= Imin[J]);
	    if(VERB>=3 && newLPdelta && m==94 && I==106 && K==0 && J==6){
	      printf("   m=%d,I=%d,K=%d,J=%d,IMIN=%d,IMAX=%d:A.ULmin=%d,A.ULmax=%d,A.Lij=%d,Lij=%d,Lijk=%d\n",m,I,K,J,IMIN,IMAX,A.ULmin(I,K,J),A.ULmax(I,K,J),A.Lij(I,K,J),Lij,Lijk);
	      fflush(stdout);
	    }
	    double Yik = Yc(Y,I,K);

	    /* update Lijk by incrementing it until Y[Lijk] overlaps Xm */
	    while(Yik-Y[Lijk] > Xm[J])
	      Lijk++;
	    if((DEBUG>=2 && !(Lijk <= I && !(Yik-Y[Lijk] > Xm[J]) && (Lijk==0 || Yik-Y[Lijk-1] > Xm[J]))) ||
	       (VERB>=3 && newLPdelta && m==94 && I==106 && K==0 && J==6)){
	      #pragma omp critical
	      {
		printf("m=%d,J=%d,I=%d,K=%d:Jmin[I]=%d,Jmax[I]=%d,Kmax[I]=%d,Yik=%0.6f,Lijk=%d,Y[Lijk]=%0.6f,Y[Lijk-1]=%0.6f,Xm[J]=%0.6f\n",
		       m,J,I,K,Jmin[I],Jmax[I],Kmax[I],Yik,Lijk,Y[Lijk],Y[Lijk-1],Xm[J]);	      
		fflush(stdout);
		assert(Lijk <= I && !(Yik-Y[Lijk] > Xm[J]) && (Lijk==0 || Yik-Y[Lijk-1] > Xm[J]));
	      }
	    }
	    A.Lij(I,K,J) = min(Lijk, I-K);
	    if(K == 0){
              if(DEBUG>=2) assert(Lijk <= I-K);
	      Lij = Lijk;
	    }

	    double NotOutlier = FLIK_FIX ? PoutlierEndM1 : FLj * PoutlierEndM1;

	    MFLOAT2 UL,bUL;
	    int Umin,Ub,Umax;
	    Ub = Umin = Umax = max(1,A.Lij(I,K,J)) - 1;
	    if(DEBUG>=2 && !(0 <= Umin && Umax < I-K)){
	      printf("I=%d,K=%d,J=%d:p->Lij=%d,Umax=%d\n",I,K,J,A.Lij(I,K,J),Umax);
	      fflush(stdout);
	      assert(0 <= Umin && Umax < I-K);
	    }
	    if(DEBUG>=2) assert(isfinite(FLj) && 0.0 <= FLj && FLj <= MaxUR * 3.0);
	    bUL = UL = FLIK<0>(I,K,J,Umin,Umin+1,Yik,Xm[J],Yik-Y[Umin],Y[Umin+1]-Y[Umin], Y,IvarX, FLj, END_FIX ? lc : 1);
	    if((MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M /* && A.UL(I,K,J)+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J) ||
	      (VERB>=3 && rverb && m==MDEBUG_M && I==33 && K==0 && J==4)){
	      printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d,A.Lij=%d,Lij=%d,Lijk=%d:UL=%0.10e(delta=%0.10e,FLIK=%0.10e,FLj=%0.10e,PoutlierEndM1=%0.10e,NotOutlier=%0.10e,PM=%0.10e):\n\t Yik=%0.4f,X[J]=%0.3f,Y[I-K]=%0.3f,Y[U]=%0.3f,Y[U+1]=%0.3f\n",
		     m,I,K,J,Umin,N,A.Lij(I,K,J), Lij,Lijk,A.UL(I,K,J)+UL*NotOutlier,UL*NotOutlier,UL,FLj,PoutlierEndM1,NotOutlier,PRtab[K][I].PM,Yik,Xm[J],Y[I-K],Y[Umin],Y[Umin+1]);
	      fflush(stdout);
	    }

	    for(int U = Umax; ++U < I-K;){
	      double LR = ERFC_FIX ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
		: FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);
	      if((MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M /* && A.UL(I,K,J)+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J) ||
		(VERB>=3 && rverb && m==MDEBUG_M && I==33 && K==0 && J==4)){
		printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e,FLIK=%0.10e,FLj=%0.10e)\n",
		       m,I,K,J,U,N,A.UL(I,K,J) + (UL+LR)*NotOutlier,LR*NotOutlier, LR, FLj);
		fflush(stdout);
	      }
	      if(LR <= UL * EPS) /* NOTE : must use <= since UL might be 0 */
		break;
	      Umax = U;
	      if(LR > bUL){
		if(URANGE_FIX){
		  Umin = U;
		  UL = 0.0;
		}
		bUL = LR;
		if(!MISSCNT_FIX)
		  A.Lij(I,K,J) = min(I-K,U+1);
	      }
	      UL += LR;
	    }

	    for(int U = Umin; --U >= 0;){
	      double LR = (ERFC_FIX && U > Ub) ? FLIK<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
		: FLIK<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);

	      if((MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M /* && A.UL(I,K,J)+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J) ||
		 (VERB>=3 && rverb && m==MDEBUG_M && I==33 && K==0 && J==4)){
		printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UL=%0.10e(delta=%0.10e,FLIK=%0.10e,FLj=%0.10e),A.Lij=%d,Lij=%d,Lijk=%d\n",
		       m,I,K,J,U,N,A.UL(I,K,J)+(UL+LR)*NotOutlier,LR*NotOutlier,LR,FLj,A.Lij(I,K,J),Lij,Lijk);
		fflush(stdout);
	      }
	      if(LR <= UL * EPS)/* NOTE : must use <= since UL might be 0 */
		break;
	      Umin = U;
	      if(LR > bUL){
		bUL = LR;
		if(!MISSCNT_FIX)
		  A.Lij(I,K,J) = U+1;/* NOTE : refalign() later uses Lij to compute misaligned sites. Here Lij is only used to compute coverage endpoint */
	      }
	      UL += LR;
	    }
	    if(DEBUG>=2) assert(Umin >= 0 && Umin < I-K);
	    if((DEBUG>=1+RELEASE && !(Umin > A.Lij(I,K,J) - 95 || Y[Umin] > Y[A.Lij(I,K,J)] - 950.0)) || (VERB>=3 && m==94 && I==106 && K==0 && J==6)){
	      printf("WARNING: m=%d:I=%d,K=%d,J=%d:IMIN=%d,IMAX=%d:A.Lij(I,K,J)= %d, Umin=%d, Y[Lij]=%0.4f, Y[Umin]=%0.4f\n",m,I,K,J,IMIN,IMAX,A.Lij(I,K,J), Umin,Y[A.Lij(I,K,J)],Y[Umin]);
	      fflush(stdout);
	      assert(Umin > A.Lij(I,K,J) - 95 || Y[Umin] > Y[A.Lij(I,K,J)] - 950.0);
	    }

	    A.ULmin(I,K,J) = Umin;
	    A.ULmax(I,K,J) = Umax;
	    if(DEBUG>=2 && !(A.ULmax(I,K,J) < I-K)){
	      printf("I=%d,K=%d,J=%d:A.ULmax(I,K,J)=%d\n",I,K,J,A.ULmax(I,K,J));
	      fflush(stdout);
	      assert(A.ULmax(I,K,J) < I-K);
	    }
	    if(DEBUG>=2) assert(0 <= Umin && Umin <= Umax);

	    if(ISLINEAR && lc){
	      double LR = FLEIK(I,K,J,0,1,Yik,Xm[J], FLj);
	      if(LR >= UL * EPS){
		A.ULmin(I,K,J) = 0;
		if(LR > bUL){
		  bUL = LR;
		  if(!MISSCNT_FIX)
		    A.Lij(I,K,J) = 0;
		}
		UL += LR;
	      }
	    }
	    if(DEBUG>=2)assert(isfinite(UL));

	    if(!EXTEND || J <= DELTA_XU || I-K <= DELTA_YU){/* don't update additional UL for deleted site so LP is same as for qprobeval */
	      UL *= NotOutlier;
	      bUL *= NotOutlier;
	      if(bUL <= A.bUL(I,K,J)) {/* A.bUL(I,K,J) corresponds to local/outlier score */
		if(DEBUG>=2) assert(A.Lij(I,K,J) >= 0);
		if(bUL <= A.bUL(I,K,J) * EPS){
		  /* NOTE : Since likelihoods for all U values are negligible, we can set ULmin == ULmax to any single value <= I-K-1 : 
		     1. For LRadd, values closer to I reduce computation cost, so I-K-1 would be optimal
		     2. For mprobeval_add, an estimate based on end of molecule (Rij + 1) would be optimal since it allows rapid estimation of U range for site added between existing sites
		            by Union of ULmin..ULmax ranges for nearest existing sites.
		     3. For LRdelta, a value corresponding to end of molecules (Rij + 1) is required, since U ranges are updated when interval sizes change in the right direction.
		  */
		  if(DEBUG>=2) assert(I-K >= 1);
		  //		  A.ULmin(I,K,J) = I - K;
		  //		  A.ULmax(I,K,J) = I - K - 1;
		  A.ULmin(I,K,J) = A.ULmax(I,K,J) = max(1,A.Lij(I,K,J)) - 1;
		}
		// WAS		A.Lij(I,K,J) = OUTLIER_FIX ? I-K : -1 - A.Lij(I,K,J);
	      } else
		A.bUL(I,K,J) = bUL;

	      if((MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M /* && A.UL(I,K,J)+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J) ||
		 (VERB>=3 && rverb && m==MDEBUG_M && I==33 && K==0 && J==4)){
		printf("m=%d:I=%d,K=%d,J=%d:ULmin=%d,ULmax=%d,N=%d:UL=%0.10e->%0.10e(delta=%0.10e,bUL=%0.10e,FLj=%0.10e,PoutlierEndM1=%0.10e, PMik=%0.10e)\n",
		       m,I,K,J,A.ULmin(I,K,J),A.ULmax(I,K,J),N,A.UL(I,K,J),(A.UL(I,K,J)+UL),UL,bUL,FLj,PoutlierEndM1,PRtab[K][I].PM);
		fflush(stdout);
	      }
	      if(DEBUG>=2) assert(A.UL(I,K,J) >= 0.0);
	      if(DEBUG>=2) assert(UL <= MaxUR * 1.001);
	      A.UL(I,K,J) += UL;
	      if(DEBUG>=2) assert(isfinite(A.UL(I,K,J)) && (!END_FIX || PoutlierEnd <= 0.0 || A.UL(I,K,J) > 0.0));
	    }
	  }
	}
      }

      if(DEBUG>=3){
        for(int I = IMIN; I <= IMAX; I++){
          for(int K = 0; K <= Kmax[I]; K++){
	    for(int J = Jmin[I]; J <= Jmax[I]; J++){
	      if(J <= DELTA_XU || (EXTEND && I-K <= DELTA_YU && J <= DELTA_XU + EXTEND_MAX)){
		int ULmin = A.ULmin(I,K,J);
		int ULmax = A.ULmax(I,K,J);
		if(DEBUG && !(ULmin == I-K-1 && ULmax == I-K) && !(0 <= ULmin && ULmin <= ULmax && ULmax < I-K)){
		  printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:ULmin=%d,ULmax=%d\n",m,I,K,J,N,M,ULmin,ULmax);
		  fflush(stdout);
		  assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		}
	      }
	    }
          }
        }
      }

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After initializing A.UL,Lij : wall time= %0.6f\n", m, tid, wtime() - wt);
	  fflush(stdout);
        }
      }

      /* initialize A[I = Imin[J]..Imax[J]][K=0..Kmax[I]][J= M+1-DELTA_X to M].Rij as rightmost site/end in Y that overlaps Xm with Y[I-K..I] is aligned with Xm[J] */
      /* Initialize corresponding UR score with right ends unaligned */
      /* only compute cases where J >= M+1 - DELTA_XU (OR EXTEND && I >= N+1 - DELTA_YU && J >= M+1 - DELTA_XU - EXTEND_MAX) : all other cases are assumed neglible */
      int jmin = EXTEND ? max(1,M+1 - DELTA_XU - EXTEND_MAX) : max(1,M+1 - DELTA_XU);
      for(int J= jmin; J <= M; J++){
	int Rij = N+1;
	int imin = (J >= M+1-DELTA_XU) ? Imin[J] : max(Imin[J], N+1-DELTA_YU - 1/* for deleted site*/);
	for(int I= Imax[J]; I >= imin; I--){
	  if(DEBUG>=2) assert(IMIN <= I && I <= IMAX);
	  if(DEBUG>=2) assert(Jmin[I] <= J && J <= Jmax[I]);
	  int Rijk = Rij;
	  if(VERB>=3 && m==43070 && I==43589 && J==36){
	    printf("\nm=%d,I=%d,K=0..%d,J=%d,M=%d,N=%d:Rij=%d\n",m,I,Kmax[I],J,M,N,Rij);
	    fflush(stdout);
	  }
	  for(int K = 0; K <= Kmax[I]; K++){
	    double Yik = Yc(Y,I,K);
	    MFLOAT delX = Xm[M+1]-Xm[J];

	    /* update Rijk by decrementing it until Y[Rijk] overlaps X */
	    while(Y[Rijk] - Yik > Xm[M+1]-Xm[J])
	      Rijk--;
	    A.Rij(I,K,J) = max(Rijk,I);
	    if(K==0)
	      Rij = Rijk;

	    if(VERB>=3 && m==43070 && I==43589 && K==0 && J==36){
	      printf("\nm=%d,I=%d,K=%d,J=%d,M=%d,N=%d:A.Rij(I,K,J) = %d, Rijk=%d, Yik= %0.4f, Y[I]= %0.4f, Y[Rijk]= %0.4f, Xm[M+1]= %0.4f, Xm[J]= %0.4f\n",
		     m,I,K,J,M,N,A.Rij(I,K,J), Rijk, Yik, Y[I], Y[Rijk], Xm[M+1], Xm[J]);
	      for(int i = I; i <= Rijk; i++)
		printf("\n\t Y[%d]= %0.6f\n",i, Y[i]);
	      fflush(stdout);
	    }

	    /* update UR score with right ends unaligned */
	    MFLOAT NotOutlier = PoutlierEndM1;
	    MFLOAT UR,bUR;
	    int Umin,Umax,Ub;
	    Ub = Umin = Umax = min(N,A.Rij(I,K,J)) + 1;
	    bUR = UR = FR<0>(I,K,J,Umin,Umin-1,Yik,delX,Y[Umin]-Yik,Y[Umin]-Y[Umin-1],Xm,Y,M,N,END_FIX ? rc : 1);
	    if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M /* && A.UL(I,K,J)+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J){
	      printf("m=%d:I=%d,K=%d,J=%d,Ub=U=%d,N=%d:UR*PM=%0.10e(delta=%0.10e)\n",
		     m,I,K,J,Umin,N,(A.UR(I,K,J)+UR*NotOutlier)*PRtab[K][I].PM,UR*NotOutlier*PRtab[K][I].PM);
	      fflush(stdout);
	    }
	    for(int U=Umin; --U > I; ){
	      MFLOAT FRu = ERFC_FIX ? FR<1>(I,K,J,U,U-1,Yik,delX,Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)
			    : FR<0>(I,K,J,U,U-1,Yik,delX,Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
	      if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M /* && A.UL(I,K,J)+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J){
		printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UR*PM=%0.10e(delta=%0.10e)\n",m,I,K,J,U,N,(A.UR(I,K,J)+(UR+FRu)*NotOutlier)*PRtab[K][I].PM,FRu*NotOutlier*PRtab[K][I].PM);
		fflush(stdout);
	      }
	      if(FRu <= UR * EPS) // NOTE : must use <= since UR may be 0
		break;
	      Umin = U;
	      if(FRu > bUR){
		if(URANGE_FIX){
		  Umax = U;
		  UR = 0.0;
		}
		bUR = FRu;
		if(!MISSCNT_FIX)
		  A.Rij(I,K,J) = U-1;/* NOTE : refalign uses Rij later to compute misaligned sites. Here Rij will only be used to compute coverage endpoint */
	      }
	      UR += FRu;
	    }
	    for(int U = Umax; ++U <= N+1; ){
	      MFLOAT FRu = (ERFC_FIX && U < Ub) ? FR<1>(I,K,J,U,U-1,Yik,delX,Y[U]-Yik, Y[U] - Y[U-1], Xm,Y,M,N,END_FIX ? rc : 1)
		: FR<0>(I,K,J,U,U-1,Yik,delX,Y[U]-Yik, Y[U] - Y[U-1], Xm,Y,M,N,END_FIX ? rc : 1);
	      if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M /* && A.UL(I,K,J)+UL >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J){
		printf("m=%d:I=%d,K=%d,J=%d,U=%d,N=%d:UR*PM=%0.10e(delta=%0.10e)\n",m,I,K,J,U,N,(A.UR(I,K,J)+(UR+FRu)*NotOutlier)*PRtab[K][I].PM,FRu*NotOutlier*PRtab[K][I].PM);
		fflush(stdout);
	      }
	      if(FRu <= UR * EPS) // NOTE : must use <= since UR may be 0
		break;
	      Umax = U;
	      UR += FRu;
	      if(FRu > bUR){
		bUR = FRu;
		if(!MISSCNT_FIX)
		  A.Rij(I,K,J) = U-1;
	      }
	    }
	    if(DEBUG>=2) assert(Umin > I);
	    if(DEBUG>=2) assert(Umin <= Umax);
	    if(DEBUG>=2) assert(Umax <= N+1);
	    if(DEBUG>=1+RELEASE && !(Umax < A.Rij(I,K,J) + 95 || Y[Umax] < Y[A.Rij(I,K,J)] + 950.0)){
	      printf("WARNING: m=%d:I=%d,K=%d,J=%d:IMIN=%d,IMAX=%d:A.Rij(I,K,J)= %d, Umax=%d,N=%d\n",m,I,K,J,IMIN,IMAX,A.Rij(I,K,J), Umax, N);
	      fflush(stdout);
	      if(DEBUG) assert(Umax < A.Rij(I,K,J) + 95 || Y[Umax] < Y[A.Rij(I,K,J)] + 950.0);
	    }
	    A.URmin(I,K,J) = Umin;
	    A.URmax(I,K,J) = Umax;

	    if(ISLINEAR && rc){
	      MFLOAT FRu = FRE(I,K,J,N+1,N,Xm[M+1]-Xm[J],Y[N+1]-Yik,M);
	      if(FRu >= UR * EPS){
		A.URmax(I,K,J) = N+1;
		if(FRu > bUR){
		  bUR = FRu;
		  if(!MISSCNT_FIX)
		    A.Rij(I,K,J) = N+1;
		}
		UR += FRu;
	      }
	    }

	    if(!EXTEND || J >= M+1-DELTA_XU || I >= N+1-DELTA_YU){/* don't update additional UR for deleted site so LP is same as for qprobeval */
	      UR = UR * NotOutlier;
	      bUR = bUR * NotOutlier;
	      if(DEBUG>=2 && A.bUR(I,K,J) != PoutlierEndf){
		#pragma omp critical
		{
		  printf("m=%d:I=%d,K=%d,J=%d:A.bUR(I,K,J)=%0.10e,PoutlierEndf=%0.10e\n",m,I,K,J,A.bUR(I,K,J),PoutlierEndf);
		  fflush(stdout);
		  assert(A.bUR(I,K,J) == PoutlierEndf);
		}
	      }
	      if(bUR <= A.bUR(I,K,J)){
		if(DEBUG>=2) assert(A.Rij(I,K,J) > 0);
		if(bUR <= A.bUR(I,K,J) * EPS){
		  /* NOTE : Since likelihoods for all U values are negligible, we can set URmin == URmax to any single value > I : 
		     1. For LRadd,LRaddT, values closer to I reduce computation cost due to Amin..Amax range, so I+1 would be optimal
		     2. For mprobeval_add, an estimate based on end of molecule (Rij + 1) would be optimal since it allows rapid estimation of U range for site added between existing sites
		            by Union of URmin..URmax ranges for nearest existing sites.
		     3. For LRdelta, a value corresponding to end of molecules (Rij + 1) is required, since U ranges are updated when interval sizes change in the right direction.
		  */
		  //		  A.URmin(I,K,J) = I+1;
		  //		  A.URmax(I,K,J) = I;
		  A.URmin(I,K,J) = A.URmax(I,K,J) = min(N,A.Rij(I,K,J)) + 1;
		}
		// WAS		A.Rij(I,K,J) = OUTLIER_FIX ? I : -A.Rij(I,K,J);
	      } else
		A.bUR(I,K,J) = bUR;

	      if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M /* && A.UR(I,K,J)+UR >= MDEBUG_MIN*/ && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J){
		printf("m=%d:I=%d,K=%d,J=%d:URmin=%d,URmax=%d,N=%d:UR*PM=%0.10e->%0.10e(delta=%0.10e,PoutlierEndM1=%0.10e)\n",
		       m,I,K,J,A.URmin(I,K,J),A.URmax(I,K,J),N,A.UR(I,K,J)*PRtab[K][I].PM,(A.UR(I,K,J)+UR)*PRtab[K][I].PM,UR*PRtab[K][I].PM,PoutlierEndM1);
		fflush(stdout);
	      }
	      if(DEBUG>=2) assert(A.UR(I,K,J) >= 0.0);
   	      if(DEBUG>=2 && !(UR <= MaxUR * 1.001)){
		#pragma omp critical
		{
		  printf("m=%d:M=%d,N=%d:I=%d,K=%d,J=%d:UR=%0.10e, MaxUR=%0.4e\n",m,M,N,I,K,J,UR,MaxUR);
		  fflush(stdout);
		  assert(UR <= MaxUR * 1.001);
		}
	      }
	      A.UR(I,K,J) += UR;
	      if(DEBUG>=2)  assert(isfinite(A.UR(I,K,J)) && (!END_FIX || A.UR(I,K,J) > 0.0));
	    }
	  }
	}
      }

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After initializing A.Rij, UR:wall time= %0.6f\n", m, tid, wtime() - wt);
	  fflush(stdout);
        }
      }

      if(DEBUG) assert(!DELTA_MIN);
      if(DEBUG) assert(PM_SPLIT);

      // HERE : duplicate vectorization code from qprobeval.cpp

      MFLOAT **delXm = pcontig->delX[0][m];

      /* Dynamic programming 3-D recurrance to compute AL */
      LogScaleL[IMIN-1] = 0.0;
      ScaleLI[IMIN-1] = 1.0;

      for(int I=IMIN; I <= IMAX; I++){
	LogScaleL[I] = LogScaleL[I-1];
	ScaleLI[I] = 1.0;/* ScaleLI[G= IMIN .. I] = exp(LogScaleL[G]-LogScaleL[I]) */
	MFLOAT2 Lscale = exp(-LogScaleL[I]);/* Lscale < 1.0 */
	if(DEBUG) assert(isfinite(Lscale) && Lscale <= 1.0);
	MFLOAT2 maxLR = 0.0;
	for(int K= 0; K <= Kmax[I]; K++){
	  MFLOAT PMik = PRtab[K][I].PM;
	  if(DEBUG>=2 && !(PMik > 0.0)){
	    printf("m=%d,I=%d,K=%d:PMik=%0.10e,PM(0,I,K,Y)=%0.10e,Y[I]=%0.3f,Y[I-K]=%0.3f\n",
		   m,I,K,PMik,PM(0,I,K,Y),Y[I],Y[I-K]);
	    printf("PRsiz=%d,PRNsiz=%d,N=%d,PRtab[K][I].PM=%0.10e\n",PRsiz,PRNsiz,N,PRtab[K][I].PM);
	  
	    fflush(stdout);
	    assert(PMik > 0.0);
	  }
	  if(DEBUG>=2) assert(fabs(log(PMik) - log(PM(0,I,K,Y))) <= (USE_MFLOAT ? 1e-5 : 2e-9));
	  double Yik = Yc(Y,I,K);
	  for(int J = Jmin[I]; J <= Jmax[I]; J++){
	    //	    MFLOAT *pFBiasWTm = &FBiasWTmRev[-J];
	    //	    MFLOAT *pLogFpPow = &LogFpPowRev[1-J];

	    if(DEBUG>=2) assert(isfinite(A.UL(I,K,J)) && A.UL(I,K,J) >= 0.0);
	    double AL = A.UL(I,K,J) * Lscale;/* UL value is unscaled, so must be multiplied by Lscale == exp(-LogScaleL[I]) */
	    if(DEBUG>=2) assert(isfinite(AL) && AL >= 0.0);
	    
	    if(J >= 2){// NEW10 : delXm[1] is NOT define AND with J = 1 the Q range (1 <= Q < J) is empty, hence no updates to AL are possible 
	      MFLOAT *delX = &delXm[J][deltaExtXRef-J];
	      int Pmin = max(IMIN, I-K - deltaExtYRef);
	      int Pmin2 = max(IMIN, I-K - DELTA_Y);
	      int Qmin = max(1, J - deltaExtXRef);
	      int Qmin2 = max(1, J - DELTA_X);
	      if(MDEBUG && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J){
		printf("m=%d,I=%d,K=%d,J=%d,P=%d..%d,Qmin=%d,Pmin2=%d,Qmin2=%d,Jmin[%d]=%d:UL=%0.10e,Lscale=%0.10e,PMik=%0.10e,AL=%0.10e\n",
		       m,I,K,J,Pmin,I-K-1,Qmin,Pmin2,Qmin2,Pmin,Jmin[Pmin],A.UL(I,K,J),Lscale,PMik,AL*PMik);
		fflush(stdout);
	      }

	      if(VERB>=3 && MDEBUG && m==MDEBUG_M && N==MDEBUG_N+1 && I==2 && K==0 && J==19){
		printf("m=%d,I=%d,K=%d,J=%d,Pmin=%d,%d,Qmin=%d,%d: UL= %0.8e, Lscale= %0.8e, LogScaleL[I]= %0.8f: AL= %0.8e\n",m,I,K,J,Pmin,Pmin2,Qmin,Qmin2,A.UL(I,K,J),Lscale,LogScaleL[I],AL);
		fflush(stdout);
	      }

	      for(int P= I-K; --P >= Pmin;){
		MFLOAT PRikn = PRtab[I-K-P-1][I-K].Pr;
		MFLOAT LogFnPowIKP = LogFnPow[I-K-P-1];
		MFLOAT PrBiasWTikp = PrBiasWT[I-K] - PrBiasWT[P];
		MFLOAT2 LRp = 0.0;
		int QminP = (P >= Pmin2) ? max(Qmin,Jmin[P]) : max(Qmin2, Jmin[P]);
		int QmaxP = min(J-1,Jmax[P]);
		if(MDEBUG && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J){
		  printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d..%d,Q=%d..%d:Jmin[P]=%d,Jmax[P]=%d: AL= %0.8e\n",m,I,K,J,P,0,Kmax[P],QminP,QmaxP,Jmin[P],Jmax[P], AL * PMik);
		  fflush(stdout);
		}

#if DEBUG>=2
		MFLOAT LRmax = 0.0, bFAj;
		int bT,bQ;
#endif

		for(int T = 0; T <= Kmax[P]; T++){
		  MFLOAT Ivar;
		  MFLOAT deltaY = Yik - Yc(Y,P,T);
		  MFLOAT2 *ALpt = &A.AL(P,T,0);

		  if(DEBUG>=2) assert(deltaY >= minKB - 2e-9);
		  MFLOAT LRik = FAIK(deltaY,PRikn,I-K-P,I,K,T,N,VARtab[K*N+I],VARtab[T*N+P],Ivar,Y);
		  if(DEBUG>=2) assert(isfinite(Ivar));

		  if(DEBUG>=2 && J <= 1) assert(QmaxP < QminP);

		  // critical loop for vectorization
#if DEBUG < 2 && MDEBUG_I < 0
#pragma omp simd reduction(+:LRp)
#endif
		  for(int Q= QminP; Q <= QmaxP; Q++){// NOTE : qprobeval iterates in ascending order (may be more accurate since it adds LR values in ascending order)
		    if(DEBUG>=2) assert(fabs((Xm[J] - Xm[Q]) - delX[Q]) < Xm[J]*((MFLOAT)1e-6));
		    //		  if(DEBUG>=2) assert(pFBiasWTm[Q] == FBiasWTm[J-Q]);
		    //		  if(DEBUG>=2) assert(pLogFpPow[Q] == LogFpPow[J-Q-1]);

		    MFLOAT FAj = FAJ(delX[Q]/*Xm[J]-Xm[Q]*/, deltaY, J-Q, I-K-P,FBiasWTm[J-Q], LogFpPow[J-Q-1], LogFnPowIKP /* LogFnPow[I-K-P-1]*/ , PrBiasWTikp, LRik, Ivar, 1, OUT_TYPE, OUT_OLD);
		    //		  MFLOAT FAj = FAJ(delX[Q]/*Xm[J]-Xm[Q]*/, deltaY, J-Q, I-K-P, pFBiasWTm[Q]/* FBiasWTm[J-Q] */, pLogFpPow[Q] /* LogFpPow[J-Q-1] */, LogFnPowIKP /* LogFnPow[I-K-P-1]*/ , LRik, Ivar, 1, OUT_TYPE, OUT_OLD);
#if DEBUG>=2
		    MFLOAT FAijk = FA(Xm[J]-Xm[Q],deltaY,J-Q,I-K-P,J,I,K,T,N,VARtab[K*N+I],VARtab[T*N+P],FBiasWTm[J-Q],LogFpPow[J-Q-1],LogFnPow[I-K-P-1],PrBiasWT[I-K]-PrBiasWT[P],PRtab[I-K-P-1][I-K].Pr,Y,1, OUT_TYPE, OUT_OLD);
		    if(DEBUG>=2) assert(isfinite(FAj) && FAj >= 0.0);
		    if(DEBUG>=2) assert(isfinite(FAijk));
		    if(DEBUG>=2 && !(fabs(FAijk - FAj) <= FAijk * (USE_MFLOAT ? 1e-4 : 2e-9) + (USE_MFLOAT ? 1e-32 : 1e-280))){
#pragma omp critical
		      {
			sleep(15);
			rverb = 1;
			MFLOAT LRik2 = FAIK(deltaY,PRikn,I-K-P,I,K,T,N,VARtab[K*N+I],VARtab[T*N+P],Ivar,Y);
			MFLOAT FAj2 = FAJ(Xm[J]-Xm[Q], deltaY, J-Q, I-K-P, FBiasWTm[J-Q], LogFpPow[J-Q-1], LogFnPow[I-K-P-1], PrBiasWTikp, LRik2, Ivar, 1, OUT_TYPE, OUT_OLD);
			MFLOAT FAijk2 = FA(Xm[J]-Xm[Q],deltaY,J-Q,I-K-P,J,I,K,T,N,VARtab[K*N+I],VARtab[T*N+P],FBiasWTm[J-Q],LogFpPow[J-Q-1],LogFnPow[I-K-P-1],PrBiasWT[I-K]-PrBiasWT[P],PRtab[I-K-P-1][I-K].Pr,Y,1, OUT_TYPE, OUT_OLD);
			printf("\nm=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:FAj=%0.10e(%0.10e),FAijk=%0.10e(%0.10e)\n",
			       m,I,K,J,P,T,Q,FAj,FAj2,FAijk,FAijk2);
			fflush(stdout);
			assert(fabs(FAijk - FAj) <= FAijk * (USE_MFLOAT ? 1e-4 : 2e-9) + (USE_MFLOAT ? 1e-32 : 1e-280));
		      }
		    }
		    if(DEBUG>=2) assert(isfinite(A.AL(P,T,Q)) && A.AL(P,T,Q) >= 0.0);
#endif // DEBUG >= 2

		    MFLOAT2 newAL = ALpt[Q] * FAj;

		    if((MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && newAL * ScaleLI[P] >= (AL + LRp*ScaleLI[P]) * 1e-8 &&
			I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J) || (VERB>=3 && rverb && m==MDEBUG_M && I==33 && K==0 && J==4)){
		      printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:FAj=%0.8e(PMik=%0.8e,LRik=%0.8e,Ivar=%0.8e),prev AL=%0.8e(Scale(P,I)=%0.8e,LSL[I]=%0.8f),LRp=%0.8e,AL=%0.8e(delta=%0.8e)\n",
			     m,I,K,J,P,T,Q,FAj,PMik,LRik,Ivar,A.AL(P,T,Q),ScaleLI[P]*PMik,LogScaleL[I], LRp * ScaleLI[P] * PMik, (AL + (LRp + newAL) * ScaleLI[P]) * PMik, newAL * PMik * ScaleLI[P]);
		      fflush(stdout);
		    }

		    LRp += newAL;

		    if(VERB>=3 && MDEBUG && m==MDEBUG_M && N==MDEBUG_N+1 && I==2 && K==0 && J==19){
		      printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:x=%0.4f,y=%0.4f,LRik= %0.8e, Ivar= %0.8e, FAj= %0.8e, ALpt[Q]= %0.8e, LRp= %0.8e, ScaleLI[P]= %0.8e, AL= %0.8e (delta= %0.8e)\n",
			     m,I,K,J,P,T,Q,delX[Q],deltaY,LRik,Ivar,FAj,ALpt[Q], LRp, ScaleLI[P], AL + LRp * ScaleLI[P], FAj * ALpt[Q] * ScaleLI[P]);
		      fflush(stdout);
		    }

#if DEBUG>=2
		    if(A.AL(P,T,Q) * FAj > LRmax){
		      LRmax = A.AL(P,T,Q) * FAj;
		      bFAj = FAj;
		      bT = T;
		      bQ = Q;
		    }
		    if(DEBUG>=2) assert(isfinite(LRp) && LRp >= 0.0 && LRp < MaxFloat);
#endif
		  } // Q = QminP .. QmaxP
		} // T = 0 .. Kmax[P]

		if(DEBUG>=2 && J <= 1) assert(LRp == 0.0 && AL <= MaxFloat);

		if(DEBUG>=2) assert(isfinite(LRp) && LRp >= 0.0 && LRp < MaxFloat);
		if(DEBUG>=2) assert(fabs(ScaleLI[P] - exp(LogScaleL[P]-LogScaleL[I])) < fabs(ScaleLI[P]) * 1e-5 + 1e-6);
		double sLR = AL + ((double)LRp) * ((double) ScaleLI[P]);
		if(!SCALE_MONOTONIC && sLR > MaxFloat){/* early re-normalization for A.AL(I,*,*) */
		  double rescale = 0.99999 * MaxLR/sLR; 
		  if(DEBUG) assert(rescale < 1.0);
		  //		double origALR = A.AL(48,0,61);
		  for(int k= 0; k <= Kmax[I]; k++)
		    for(int J = Jmin[I]; J <= Jmax[I]; J++)
		      A.AL(I,k,J) *= rescale;
		  //		if(VERB>=3 && rverb && m==MDEBUG_M && I==48 && J <= 61)
		  //		  printf("m=%d,I=%d,K=%d,J=%d(%d..%d):A.AL(I,0,61) = %0.8e -> %0.8e, rescale= %0.8e\n",m,I,K,J,Jmin[I],Jmax[I],origALR,A.AL(I,0,61),rescale);

		  double origLogScale = LogScaleL[I];
		  double origScale = ScaleLI[P];
		  LogScaleL[I] -= log(rescale);
		  for(int G= max(IMIN,I - A.kmax - deltaExtYRef); G < I; G++)
		    ScaleLI[G] *= rescale;

		  if(VERB/* HERE >=2 */ && MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J){
		    printf("I=%d,K=%d,J=%d,P=%d:LogScaleL[I]=%0.8f->%0.8f,maxLR=%0.8e->%0.8e,PMik=%0.8e,LRp=%0.8e,AL=%0.8e->%0.8e,ScaleLI[P]=%0.8e->%0.8e\n",
			   I,K,J,P,origLogScale,LogScaleL[I],maxLR,maxLR * rescale, PMik, LRp * ScaleLI[P] * PMik, AL * PMik, AL*rescale * PMik, origScale,ScaleLI[P]);
#if DEBUG>=2
		    printf("\t best T=%d,Q=%d : LR = %0.8e, A.AL(P,T,Q) = %0.8e, FAj = %0.8e\n",bT,bQ, LRmax, A.AL(P,bT,bQ), bFAj);
#endif
		    for(int g = max(IMIN, I - A.kmax - deltaExtYRef); g < I; g++)
		      printf("\t ScaleLI[%d] -> %0.8e\n", g, ScaleLI[g]);
		    fflush(stdout);
		  }
		  AL *= rescale;
		  maxLR *= rescale;
		  Lscale = exp(-LogScaleL[I]);
		  if(DEBUG) assert(isfinite(Lscale) && Lscale <= 1.0);
		}
		if((DEBUG>=2 && LRp != 0.0 && !(isfinite(ScaleLI[P]) && ScaleLI[P] >= 0.0 && ScaleLI[P] <= DBL_MAX)) ||
		   (MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J) ||
		   (VERB>=3 && rverb && m==MDEBUG_M && I==102 && K==0 && J==31)){
		  printf("m=%d:I=%d,K=%d,J=%d,P=%d:AL=%0.8e,LRp=%0.8e(ScaleLI[P]=%0.8e,PMik=%0.8e):AL*PMik=%0.8e -> %0.8e,LogScaleL[P,I]=%0.8f,%0.8f\n",
			 m,I,K,J,P,AL,LRp * ScaleLI[P] * PMik, ScaleLI[P],PMik, AL*PMik, (AL + LRp * ScaleLI[P]) * PMik, LogScaleL[P],LogScaleL[I]);
		  fflush(stdout);
		  if(LRp != 0.0) assert(isfinite(ScaleLI[P]) && ScaleLI[P] >= 0.0 && ScaleLI[P] <= DBL_MAX);
		}
		if(isfinite(ScaleLI[P]) || LRp != 0.0){// NEW144 : otherwise assume ScaleLI[P] * LRp == 0		
		  if(DEBUG>=2) assert(isfinite(ScaleLI[P]) && fabs(ScaleLI[P] - exp(LogScaleL[P]-LogScaleL[I])) < fabs(ScaleLI[P]) * 1e-5 + 1e-6);
		  AL += LRp * ScaleLI[P];/* ScaleLI[P] == exp(LogScaleL[P]-LogScaleL[I]) */
		  if(DEBUG>=2) assert(isfinite(AL) && AL >= 0.0 && AL < MaxFloat);
		}
	      } // P = I-K-1 .. Pmin
	    }// if(J >= 2)

	    if(DEBUG>=2) assert(isfinite(AL) && AL >= 0.0 && AL < MaxFloat);
	    if(DEBUG>=2) assert(isfinite(PMik) && PMik >= 0.0);
	    if((MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)
	       || (VERB>=3 && rverb && m==MDEBUG_M && I==33 && K==0 && J==4/* && MD==792 && addcnt==468*/)){
	      printf("m=%d:I=%d,K=%d,J=%d:AL=%0.10e,UL=%0.10e,Lscale=%0.10e,PMik=%0.10e,LogScaleL[I]=%0.6e\n",
		     m,I,K,J,AL*PMik,A.UL(I,K,J),Lscale,PMik,LogScaleL[I]);
	      fflush(stdout);
	    }
	    MFLOAT2 nAL = A.AL(I,K,J) = AL * PMik;
	    if(DEBUG && !(isfinite(nAL) && nAL >= 0.0 && nAL < MaxFloat)){
              #pragma omp critical
	      {
		printf("\nmprobeval:thread=%d:",tid);
		printf("m=%d:I=%d,K=%d,J=%d:AL(I,K,J)->%0.10e,UL=%0.10e,AL=%0.10e,PMik=%0.10e,LogScaleL[I]=%0.10f,Lscale=%0.10e\n",
		       m,I,K,J,A.AL(I,K,J),A.UL(I,K,J),AL, PMik,LogScaleL[I],Lscale);
		fflush(stdout);
		assert(isfinite(nAL) && nAL >= 0.0 && nAL < MaxFloat);
	      }
	    }
	    maxLR = max(nAL, maxLR);
	  }
	}

	if(DEBUG) assert(isfinite(maxLR) && maxLR >= ((MFLOAT)0.0) && maxLR < MaxFloat);
	if(!SCALE_MONOTONIC || maxLR > MaxLR){ /* need to rescale A[I][*][*].AL so largest value (maxLR) is <= MaxLR and update LogScaleL[I] and ScaleLI[G=max(IMIN,I-DELTA_Y)..I-1] */
	  MFLOAT2 rescale = ((maxLR > 1.0 || maxLR > (MaxLR/MaxFloat)) ? max(1.00001f * MinFloat, 0.99999f * MaxLR/maxLR) : 0.99999f * MaxFloat);/* multiplier for A[I][*][*].AL  and ScaleLI[IMIN..I-1] */
	  if(DEBUG) assert(MinFloat < rescale && rescale < MaxFloat);
	  double MaxScaleI = DBL_MAX;
	  if(!SCALE_MONOTONIC && rescale > 1.0){
	    MFLOAT2 origrescale = rescale;

	    for(int G= max(IMIN, I - A.kmax - deltaExtYRef); G < I; G++){/* avoid overflow of ScaleLI[] */
	      if(rescale <= 1.0)
		break;
	      if(ScaleLI[G] >= MaxScaleI/rescale){
		MFLOAT2 rescale2 = rescale;
		//		rescale = (long double)(0.99999 * MaxScaleI) / (long double) ScaleLI[G];
		rescale = 0.99999 * MaxScaleI / ScaleLI[G];
		if(DEBUG>=1+RELEASE && !(MinFloat < rescale && rescale < MaxFloat)){
		  printf("\nI=%d,IMAX=%d,A.kmax=%d,deltaExtYRef=%d,G=%d,maxLR=%0.17e,MaxLR=%0.17e,MinFloat=%0.6e,MaxFloat=%0.6e,MaxScaleI=%0.6e,ScaleLI[G]=%0.17e:rescale= %0.17e -> %0.17e -> %0.17e\n",
			 I,IMAX,G,A.kmax,deltaExtYRef,maxLR,MaxLR,MinFloat,MaxFloat,MaxScaleI,ScaleLI[G],origrescale,rescale2,rescale);
		  fflush(stdout);
		  assert(MinFloat < rescale && rescale < MaxFloat);
		}
	      }
	    }
	    if(rescale > 1.0){/* make sure LogScaleL[I] does not drop below 0.0 */
	      if(LogScaleL[I] - log(rescale) < 0.0){
		MFLOAT2 rescale2 = rescale;
		rescale = exp(LogScaleL[I]);
		if(DEBUG>=1+RELEASE && !(MinFloat < rescale && rescale < MaxFloat)){
		  printf("\nI=%d,maxLR=%0.17e,MaxLR=%0.17e,MinFloat=%0.6e,MaxFloat=%0.6e,LogScaleL[I]=%0.17f:rescale= %0.17e -> %0.17e -> %0.17e\n",
			 I,maxLR,MaxLR,MinFloat,MaxFloat,LogScaleR[I],origrescale,rescale2,rescale);
		  fflush(stdout);
		  assert(MinFloat < rescale && rescale < MaxFloat);
		}
	      }
	    }
	    if(DEBUG && !(MinFloat < rescale && rescale < MaxFloat)){
	      printf("\nI=%d,maxLR=%0.17e,MaxLR=%0.17e,MinFloat=%0.6e,MaxFloat=%0.6e:rescale= %0.17e -> %0.17e\n",
		     I,maxLR,MaxLR,MinFloat,MaxFloat,origrescale,rescale);
	      fflush(stdout);
	      assert(MinFloat < rescale && rescale < MaxFloat);
	    }
	  }
	  if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792*/){
	    printf("m=%d:I=%d,maxLR=%0.10e,MaxLR=%0.4e,rescale=%0.10e,LogScaleL[I]=%0.10f (log(maxLR)+LogScaleL[I]=%0.10f)\n", 
		  m,I,maxLR,MaxLR,rescale,LogScaleL[I], (maxLR > 0.0 ? log(maxLR) : -1000.0) + LogScaleL[I]);
	    fflush(stdout);
	  }
	  for(int K= 0; K <= Kmax[I]; K++)
	    for(int J = Jmin[I]; J <= Jmax[I]; J++){
	      if(DEBUG>=2) assert(A.AL(I,K,J) <= maxLR);
	      A.AL(I,K,J) *= rescale;
	      if(DEBUG>=2) assert(isfinite(A.AL(I,K,J)) && A.AL(I,K,J) >= 0.0);
	    }
	  double origLogScale = LogScaleL[I];
	  LogScaleL[I] -= log(rescale);
	  if(DEBUG) assert(LogScaleL[I] >= -1e-5);
	  LogScaleL[I] = max(0.0, LogScaleL[I]);
	  for(int G= max(IMIN, I - A.kmax - deltaExtYRef); G < I; G++){
	    ScaleLI[G] *= rescale;
	    if(DEBUG) assert(isfinite(ScaleLI[G]) && ScaleLI[G] <= (SCALE_MONOTONIC ? 1.0 : MaxScaleI));
	  }
	  if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 && origIMIN==430 && origIMAX==538 && addcnt==468 */){
	    printf("I=%d,LogScaleL[I]= %0.8f -> %0.8f, maxLR= %0.8e -> %0.8e, A.kmax=%d\n",I,origLogScale,LogScaleL[I],maxLR,maxLR * rescale, A.kmax);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(ScaleLI[I] == 1.0);
	  // maxLR *= rescale;// redundant since maxLR is about to go out of scope
	}

	if(VERB>=3 && MDEBUG && m==MDEBUG_M && (N==MDEBUG_N+1 || N==MDEBUG_N)){
	  printf("m=%d:I=%d:LogScaleL[I]= %0.10f\n",m,I,LogScaleL[I]);
	  if(I<=4)
	    for(int K= 0; K <= Kmax[I]; K++)
	      for(int J = Jmin[I]; J <= Jmax[I]; J++)
		printf("\t m=%d:I=%d,K=%d,J=%d: AL(I,K,J)= %0.10e, UR(I,K,J)= %0.10e\n", m,I,K,J, A.AL(I,K,J), A.UR(I,K,J));
	}

	if(VERB>=3 && rverb && m==MDEBUG_M){
	  for(int K= 0; K <= Kmax[I]; K++)
	    for(int J = Jmin[I]; J <= Jmax[I]; J++)
	      printf("\t m=%d:I=%d,K=%d,J=%d: AL(I,K,J)= %0.10e, UR(I,K,J)= %0.10e, LogScaleL[I]= %0.10f\n", m,I,K,J, A.AL(I,K,J), A.UR(I,K,J), LogScaleL[I]);
	  fflush(stdout);
	}
      }
      LogScaleL[IMAX+1] = LogScaleL[IMAX];

      /* Dynamic programming 3-D recurrance to compute AR */
      LogScaleR[IMAX+1] = 0.0;
      ScaleRI[IMAX+1] = 1.0;

      for(int I = IMAX; I >= IMIN; I--){
	LogScaleR[I] = LogScaleR[I+1];
	ScaleRI[I] = 1.0;/* ScaleRI[G= I .. IMAX] = exp(LogScaleR[G]-LogScaleR[I]) */
	MFLOAT2 Rscale = exp(-LogScaleR[I]);
	if(DEBUG) assert(isfinite(Rscale) && Rscale <= 1.0);
	MFLOAT2 maxLR = 0.0;
	for(int K= 0; K <= Kmax[I]; K++){
	  MFLOAT PMik = PRtab[K][I].PM;
	  if(DEBUG>=2) assert(PMik > 0.0);
	  if(DEBUG>=2) assert(fabs(log(PMik) - log(PM(0,I,K,Y))) <= (USE_MFLOAT ? 1e-5 : 2e-9));
	  double Yik = Yc(Y,I,K);
	  for(int J = Jmax[I]; J >= Jmin[I]; J--){
	    if(DEBUG>=2) assert(isfinite(A.UR(I,K,J)) && A.UR(I,K,J) >= 0.0);
	    MFLOAT2 AR = A.UR(I,K,J) * Rscale;/* UR value is unscaled, so must be multiplied by exp(-LogScaleR[I]) */
	    if(VERB>=2 && MDEBUG_ST >= 0 && m==0 && I==53 && K==0 && J==1 && AR){
	      printf("m=%d,I=%d,K=%d,J=%d:AR[I][K][J](init)=%0.10e:UR=%0.10e,PMik=%0.10e\n",
		     m,I,K,J,AR*PMik*exp(LogScaleR[I]),A.UR(I,K,J),PMik);
	      fflush(stdout);
	    }

#if PQ_BND
	    memset(LRmem, 0, deltaExtYRef * deltaExtXRef * sizeof(double));
#endif	    

	    int Pmax = min(IMAX, I + deltaExtYRef);
	    int Pmax2 = min(IMAX, I + DELTA_Y);
	    int Qmax = min(M, J + deltaExtXRef);
	    int Qmax2 = min(M, J + DELTA_X);

	    if(VERB>=3 && MDEBUG && m==MDEBUG_M && N==MDEBUG_N+1 && I==2 && K==0 && J==19){
	      printf("m=%d,I=%d,K=%d,J=%d,Pmax=%d,%d,Qmax=%d,%d: UR= %0.8e, Rscale= %0.8e, LogScaleR[I]= %0.8f: AR= %0.8e\n",m,I,K,J,Pmax,Pmax2,Qmax,Qmax2,A.UR(I,K,J),Rscale,LogScaleR[I],AR);
	      fflush(stdout);
	    }

	    for(int P= I; ++P <= IMAX && P <= Pmax + A.kmax /* WAS KMAX */;){
	      MFLOAT2 LRp = 0.0;
	      int Tmax = min(P-I-1,Kmax[P]);
	      int QminP = max(J+1,Jmin[P]);

	      for(int T = max(0,P-Pmax); T <= Tmax; T++){
		if(DEBUG>=2) assert(T <= A.kmax && P <= Pmax + A.kmax && P <= Pmax + Kmax[P]);
		int QmaxP = (P - T <= Pmax2) ? min(Qmax,Jmax[P]) : min(Qmax2, Jmax[P]);
		if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && (I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J)){
		  printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d:Jmin[P]=%d,Jmax[P]=%d,QminP=%d,QmaxP=%d\n", m,I,K,J,P,T,Jmin[P],Jmax[P],QminP,QmaxP);
		  fflush(stdout);
		}
		MFLOAT Ivar;
		MFLOAT deltaY = Yc(Y,P,T) - Yik;
		MFLOAT2 *ARpt = &A.AR(P,T,0);
		MFLOAT LRik = FAIK(deltaY,PRtab[P-T-I-1][P-T].Pr, P-T-I,P,T,K,N,VARtab[T*N+P],VARtab[K*N+I],Ivar,Y);
		if(DEBUG>=2) assert(isfinite(Ivar));
#if PQ_BND
		double *const __restrict__ LRptq = &LRmem2[(P-T-I-1) * deltaExtXRef - J  - 1];
#endif

		// critical loop for vectorization
		#if DEBUG < 2
		#pragma omp simd reduction(+:LRp)
		#endif
		for(int Q= QminP; Q <= QmaxP; Q++){
		  MFLOAT FAj = FAJ(Xm[Q]-Xm[J], deltaY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], LRik,Ivar,1, OUT_TYPE, OUT_OLD);
#if DEBUG>=2
		  MFLOAT FAijk = FA(Xm[Q]-Xm[J],deltaY,Q-J,P-T-I,Q,P,T,K,N,VARtab[T*N+P],VARtab[K*N+I],FBiasWTm[Q-J],LogFpPow[Q-J-1],LogFnPow[P-T-I-1],PrBiasWT[P-T]-PrBiasWT[I],PRtab[P-T-I-1][P-T].Pr,Y,1, OUT_TYPE, OUT_OLD);
		  if(DEBUG && !(isfinite(FAj) && FAj >= 0.0)){
                    #pragma omp critical
		    {
		      printf("\nmprobeval:thread=%d:",tid);
		      printf("m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:FAj=%0.10e,FAijk=%0.10e\n",
		        m,I,K,J,P,T,Q,FAj,FAijk);
		      fflush(stdout);
		      assert(isfinite(FAj) && FAj >= 0.0);
		    }
		  }
		  if(DEBUG>=2) assert(isfinite(FAijk) && FAijk >= 0.0);
		  if(DEBUG>=2 && !(fabs(FAijk - FAj) <= FAijk * (USE_MFLOAT ? 1e-4 : 2e-9) + (USE_MFLOAT ? 1e-32 : 1e-280))){
		    printf("WARNING:m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:FAj=%0.10e,FAijk=%0.10e\n",
		      m,I,K,J,P,T,Q,FAj,FAijk);
		    fflush(stdout);
		    assert(fabs(FAijk - FAj) <= FAijk * (USE_MFLOAT ? 1e-3 : 2e-9) + (USE_MFLOAT ? 1e-32 : 1e-280));
		  }
		  if(DEBUG>=2) assert(isfinite(A.AR(P,T,Q)) && A.AR(P,T,Q) >= 0.0);
#endif // DEBUG >= 2

		  if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && A.AR(P,T,Q) * FAj*ScaleRI[P] >= (AR+LRp*ScaleRI[P]) * 1e-8 && (I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J) && MD==585 && addcnt==407){
		    printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:FAj=%0.10e(PMik=%0.10e,LRik=%0.10e),prev AR=%0.10e(ScaleRI[P]=%0.10e,LSR[I]=%0.8f,LSR[P]=%0.8f),AR=%0.10e(delta=%0.10e)\n",
			   m,I,K,J,P,T,Q,FAj*PMik,PMik,LRik,A.AR(P,T,Q)*ScaleRI[P],ScaleRI[P],LogScaleR[I],LogScaleR[P],(AR + (LRp + A.AR(P,T,Q) * FAj) * ScaleRI[P])*PMik, A.AR(P,T,Q)*FAj*PMik*ScaleRI[P]);
		    fflush(stdout);
		  }

		  MFLOAT2 LR = ARpt[Q] * FAj;
		  if(DEBUG>=2) assert(isfinite(LR));

		  LRp += LR;

		  if(VERB>=3 && MDEBUG && m==MDEBUG_M && N==MDEBUG_N+1 && I==2 && K==0 && J==19){
		    printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:x=%0.4f,y=%0.4f,LRik= %0.8e, Ivar= %0.8e, FAj= %0.8e, ARpt[Q]= %0.8e, LRp= %0.8e, ScaleRI[P]= %0.8e, AR= %0.8e (delta= %0.8e)\n",
			   m,I,K,J,P,T,Q,Xm[Q]-Xm[J],deltaY,LRik,Ivar,FAj,ARpt[Q],LRp, ScaleRI[P], AR + LRp * ScaleRI[P], FAj * ARpt[Q] * ScaleRI[P]);
		    fflush(stdout);
		  }
#if PQ_BND
		  LRptq[Q] = LR;// save AR sub-totals for P-T,Q for later application of EPS3
#endif

		  if((DEBUG>=2 && !(isfinite(LRp) && LRp >= 0.0)) || (VERB>=2 && MDEBUG_ST >= 0 && m==0 && I==53 && K==0 && J==1)){
		    double FAikj = FAJ(Xm[Q]-Xm[J], deltaY, Q-J, P-T-I, FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T]-PrBiasWT[I], LRik,Ivar,1, OUT_TYPE, OUT_OLD);
		    if(1 || A.AR(P,T,Q) * FAikj*ScaleRI[P]*PMik > MDEBUG_MIN * 0.000001)
		      printf("m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d:AR[I][K][J]=%0.10e(delta=%0.10e):AR[P][T][Q]=%0.10e,FA=%0.10e,PMik=%0.10e\n",
			     m,I,K,J,P,T,Q,(AR+LRp*ScaleRI[P])*PMik*exp(LogScaleR[I]), A.AR(P,T,Q) * FAikj * ScaleRI[P] * PMik * exp(LogScaleR[I]), A.AR(P,T,Q) * exp(LogScaleR[P]), FAikj, PMik);
		    assert(isfinite(LRp) && LRp >= 0.0);
		  }
		}// for Q = QminP .. QmaxP
	      }// T = max(0,P-Pmax) .. Tmax

	      if(DEBUG>=2) assert(isfinite(LRp) && LRp >= 0.0 && LRp < MaxFloat);
	      if(DEBUG>=2) assert(fabs(ScaleRI[P] - exp(LogScaleR[P]-LogScaleR[I])) < fabs(ScaleRI[P]) * 1e-5 + 1e-6);

	      double LR = AR + ((double)LRp) * ((double) ScaleRI[P]);
	      if(!SCALE_MONOTONIC && LR > MaxFloat){/* early renormalization for A.AR(I,*,*) */
		double rescale = 0.99999 * MaxLR/LR;

		double origLogScale = LogScaleR[I];
		double origAR = AR;
		double origScaleRIP = ScaleRI[P];
		if(DEBUG) assert(rescale < 1.0);
		for(int k= 0; k <= Kmax[I]; k++)
		  for(int J = Jmin[I]; J <= Jmax[I]; J++)
		    A.AR(I,k,J) *= rescale;
		LogScaleR[I] -= log(rescale);
		AR *= rescale;
		
		// NOTE : no need to rescale LRp OR LRmem2[] since it will be rescaled by ScaleRI[P] later
#if 0 
		for(int T = max(0,P-Pmax); T <= Tmax; T++){
		  int QmaxP = (P - T <= Pmax2) ? min(Qmax,Jmax[P]) : min(Qmax2, Jmax[P]);
		  double *LRptq = &LRmem2[(P-T-I-1) * deltaExtXRef - J  - 1];
		  for(int Q = QminP; Q <= QmaxP; Q++)
		    LRptq[Q] *= rescale;
	        }
#endif
		for(int G= min(IMAX,I + A.kmax + deltaExtYRef); G > I; G--)
		  ScaleRI[G] *= rescale;

		if(VERB>=3 && ((MDEBUG && m==MDEBUG_M && N==MDEBUG_N+1 && I==2 && K==0 && J==19) ||
		    ((m==130 && I==1 && K==0 && J==30)||(m==168 && I==1 && K==0 && J==2))) && addcnt==44){
		  printf("I=%d,K=%d,J=%d,P=%d:LR= %0.8e, LRp= %0.8e, rescale= %0.8e: AR= %0.8e -> %0.8e, LogScaleR[I]=%0.8f -> %0.8f, maxLR= %0.8e -> %0.8e, ScaleRI[P]= %0.8e -> %0.8e\n",
		    I,K,J,P,LR,LRp, rescale, origAR, AR, origLogScale,LogScaleR[I],maxLR,maxLR * rescale, origScaleRIP, ScaleRI[P]);
		  fflush(stdout);
		}
		maxLR *= rescale;
		Rscale = exp(-LogScaleR[I]);// NEW
		if(DEBUG) assert(isfinite(Rscale) && Rscale <= 1.0);
	      }

              if(isfinite(ScaleRI[P]) || LRp != 0.0){// NEW144 : otherwise assume ScaleRI[P] * LRp == 0
	        if(DEBUG>=2) assert(isfinite(ScaleRI[P]) && ScaleRI[P] >= 0.0 && ScaleRI[P] <= (SCALE_MONOTONIC ? 1.0 : DBL_MAX));
     	        if(DEBUG>=2) assert(fabs(ScaleRI[P] - exp(LogScaleR[P]-LogScaleR[I])) < fabs(ScaleRI[P]) * 1e-5 + 1e-6);
	        AR += LRp * ScaleRI[P];/* ScaleRI[P] == exp(LogScaleR[P]-LogScaleR[I]) */
	        if(DEBUG>=2) assert(isfinite(AR) && AR >= 0.0 && AR < MaxFloat);
              }

#if PQ_BND
	      for(int T = max(0,P-Pmax); T <= Tmax; T++){
	        int QmaxP = (P - T <= Pmax2) ? min(Qmax,Jmax[P]) : min(Qmax2, Jmax[P]);
		double *LRpq = &LRmem[(P-T-I-1) * deltaExtXRef - J  - 1];
		double *LRptq = &LRmem2[(P-T-I-1) * deltaExtXRef - J  - 1];
		for(int Q= QminP; Q <= QmaxP; Q++){
	          LRpq[Q] += LRptq[Q] * ScaleRI[P];
		  if(DEBUG>=2) assert(isfinite(LRpq[Q]));
	        }
	      }
#endif
	    } // P loop

	    if(DEBUG>=2) assert(isfinite(AR) && AR >= 0.0 && AR < MaxFloat);
	    if(DEBUG>=2) assert(isfinite(PMik) && PMik >= 0.0);
	    if(MDEBUG_I >= 0 && MDEBUG_M >= 0 && MDEBUG_N1 <= N && N <= MDEBUG_N2 && m==MDEBUG_M && ((I==MDEBUG_I && K==MDEBUG_K && J==MDEBUG_J) || (0 && I==1 && K==0 && J==8))){
	      printf("m=%d:I=%d,K=%d,J=%d:AR=%0.10e,UR=%0.10e,PMik=%0.10e,LogScaleR[I]=%0.10f\n",
		     m,I,K,J,AR*PMik,A.UR(I,K,J),PMik,LogScaleR[I]);
	      fflush(stdout);
	    }

#if PQ_BND
	    // compute A.Pmin(I,K,J),A.Pmax(I,K,J),A.Qmin(I,K,J),A.Qmax(I,K,J), by looping over all AR terms again and checking which ones are >= AR * E1 and use them to create a bounding
	    // box P=Pmin..Pmax and Q=Qmin..Qmax that is used to intersect with the DELTA_X,DELTA_Y based box. Note that E1 < EPS3, since many such terms outside the bounding box are being discarded.
	    // Approximately, E1 * (Pmax-Pmin+Qmax-Qmin+2)*2  ~ EPS3, hence E1 ~ EPS3_scaled
	    // NOTE : Pmin1..Pmax1 is range of P-T
	    // HERE : add statistics on savings

	    if(DEBUG>=2){
	      if(VERB>=3 && MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M && MDEBUG_S >= 0 && N==MDEBUG_N && I==1 && J == 2 && K== 0){
		printf("m=%d:I=%d,K=%d,J=%d:UL=%0.8e,AR=%0.8e,Rscale=%0.8e:\n", m,I,K,J,A.UL(I,K,J), AR, Rscale);
		fflush(stdout);
	      }
	      MFLOAT2 AR2 = A.UR(I,K,J) * Rscale;
	      for(int P= 0; P < deltaExtYRef; P++){
	        double *LRpq = &LRmem[P * deltaExtXRef];
		int Qmax3 = (P < DELTA_Y) ? deltaExtXRef : DELTA_X;
		for(int Q = 0; Q < Qmax3; Q++){
		  AR2 += LRpq[Q];
		  if(VERB>=3 && MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M && MDEBUG_S >= 0 && N==MDEBUG_N && I==1 && J == 2 && K== 0){
		    MFLOAT2 x = X[J+Q+1] - X[J];
		    MFLOAT2 y = Y[I+P+1] - Yc(Y,I,K);
		    printf("\t\t LRpq[%d][%d] = %0.8e, LRpq/AR = %0.8e, x = %0.4f, y <= %0.4f\n",P,Q,LRpq[Q], LRpq[Q]/AR, x, y);
		  }
		}
	      }
	      if(!(fabs(AR2-AR) <= AR * (USE_MFLOAT ? 1e-4 : 1e-8) + (USE_MFLOAT ? 1e-32 : 1e-300))){
                #pragma omp critical
		{
		  printf("\nm=%d,I=%d,K=%d,J=%d:AR=%0.10e,AR2=%0.10e,(err=%0.10e),PMik=%0.10f\n",m,I,K,J,AR,AR2,AR2-AR,PMik);
		  fflush(stdout);
		  assert(fabs(AR2-AR) <= AR * ( USE_MFLOAT ? 1e-4 : 1e-8) + (USE_MFLOAT ? 1e-32 : 1e-300));
	        }
	      }
	    }
	    int Pmin1 = 0;
	    int Pmax1 = deltaExtYRef - 1;
	    int Qmin1 = 0;
	    int Qmax1 = deltaExtXRef - 1;;
	    if(EPS3 > 0.0){
	      double ARloss = AR * EPS3;/* LR loss allowance remaining */
	      while(Pmin1 <= Pmax1 && Qmin1 <= Qmax1){
		MFLOAT deltaPmin = 0.0, deltaPmax = 0.0, deltaQmin = 0.0, deltaQmax = 0.0;

		/* accumulate LR values that would be lost if Pmin1 is incremented by 1 */
		double *LRpq = &LRmem[Pmin1 * deltaExtXRef];
		int QmaxPmin = (Pmin1 < DELTA_Y) ?  Qmax1 : min(Qmax1, DELTA_X - 1);
		for(int Q = Qmin1; Q <= QmaxPmin; Q++)
		  deltaPmin += LRpq[Q];

		MFLOAT Loss = deltaPmin, LossRatio = deltaPmin / max(1, QmaxPmin - Qmin1 + 1);/* record best ratio of lost LR values to number of P,Q index values saved */
		int BestLoss = 0;/* 0 means Pmin1, 1 means Pmax1, 2 means Qmin1, 3 means Qmax1 */

		/* accumulate LR values that would be lost if Pmax1 is decremented by 1 */
		LRpq = &LRmem[Pmax1 * deltaExtXRef];
		int QmaxPmax = (Pmax1 < DELTA_Y) ? Qmax1 : min(Qmax1, DELTA_X - 1);
		for(int Q = Qmin1; Q <= QmaxPmax; Q++)
		  deltaPmax += LRpq[Q];
		if(deltaPmax < (QmaxPmax - Qmin1 + 1) * LossRatio && deltaPmax < ARloss){
		  Loss = deltaPmax;
		  LossRatio = deltaPmax / max(1,QmaxPmax - Qmin1 + 1);
		  BestLoss = 1;
		}

		/* accumulate LR values that would be lost if Qmin1 is incremented by 1 */
		int PmaxQmin = (Qmin1 < DELTA_X) ? Pmax1 : min(Pmax1, DELTA_Y - 1);
		LRpq = &LRmem[Pmin1*deltaExtXRef + Qmin1];
		for(int P = Pmin1; P <= PmaxQmin; P++){
		  deltaQmin += *LRpq;
		  LRpq += deltaExtXRef;
		}
		if(deltaQmin < (PmaxQmin - Pmin1 + 1) * LossRatio && deltaQmin < ARloss){
		  Loss = deltaQmin;
		  LossRatio = deltaQmin / max(1, PmaxQmin - Pmin1 + 1);
		  BestLoss = 2;
		}

		/* accumulate LR values that would be lost if Qmax1 is decremented by 1 */
		int PmaxQmax = (Qmax1 < DELTA_X) ? Pmax1 : min(Pmax1, DELTA_Y - 1);
		LRpq = &LRmem[Pmin1*deltaExtXRef + Qmax1];
		for(int P = Pmin1; P <= PmaxQmax; P++){
		  deltaQmax += *LRpq;
		  LRpq += deltaExtXRef;
		}
		if(deltaQmax < (PmaxQmax - Pmin1 + 1) * LossRatio && deltaQmax < ARloss){
		  Loss = deltaQmax;
		  LossRatio = deltaQmax / max(1, PmaxQmax - Pmin1 + 1);
		  BestLoss = 3;
		}

		if(VERB>=3 && MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M && MDEBUG_S >= 0 && N==MDEBUG_N && I==1 && J == 2 && K== 0 && LPdelG){
		  #pragma omp critical
		  {
		    printf("\t ARloss/AR=%0.8e:Pmin1=%d,Pmax1=%d,Qmin1=%d,Qmax1=%d:BestLoss=%d(Loss/AR=%0.8e,LossRatio/AR=%0.8e),dPmin/AR=%0.8e,dPmax/AR=%0.8e,dQmin/AR=%0.8e,dQmax/AR=%0.8e\n",
			   ARloss/AR,Pmin1,Pmax1,Qmin1,Qmax1,BestLoss,Loss/AR,LossRatio/AR, deltaPmin/AR, deltaPmax/AR,deltaQmin/AR,deltaQmax/AR);
		    fflush(stdout);
		  }
		}

		if(Loss > ARloss)
 		  break;

		ARloss -= Loss;

		switch(BestLoss){
		case 0:
		  Pmin1++;
		  break;
		case 1:
		  Pmax1--;
		  break;
		case 2:
		  Qmin1++;
		  break;
		case 3:
		  Qmax1--;
		  break;
		}
	      }

	      if(EPS_ADJUST){	      /* adjust AR to reflect reduced P,Q ranges */
		AR *= (1.0 - EPS3);
		AR += ARloss;
	      }

	      if(VERB>=3 && MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M && MDEBUG_S >= 0 && N==MDEBUG_N && I==17 && J == 25 && K== 0 && LPdelG){
                #pragma omp critical
		{
		  printf("m=%d,I=%d,K=%d,J=%d:UL=%0.8e,AR=%0.8e,ARloss/AR=%0.8e:Pmin1=%d,Pmax1=%d,Qmin1=%d,Qmax1=%d\n",
			 m,I,K,J,A.UL(I,K,J),AR,ARloss/AR,Pmin1,Pmax1,Qmin1,Qmax1);
		  fflush(stdout);
		}
	      }
	    }

	    A.Pmin(I,K,J) = I+1 + Pmin1;
	    A.Pmax(I,K,J) = I+1 + Pmax1;
	    A.Qmin(I,K,J) = J+1 + Qmin1;
	    A.Qmax(I,K,J) = J+1 + Qmax1;
#endif

	    A.AR(I,K,J) = AR *= PMik;
	    if(DEBUG && !(isfinite(AR) && AR>=0.0 && AR < MaxFloat)){
              #pragma omp critical
	      {
		printf("mprobeval:thread=%d:",tid);
		printf("m=%d:I=%d,K=%d,J=%d:AR=%0.10e,UR=%0.10e,PMik=%0.10e,LogScaleR[I]=%0.10f,Rscale=%0.10e\n",
		       m,I,K,J,AR*PMik,A.UR(I,K,J),PMik,LogScaleR[I],Rscale);
		fflush(stdout);
		assert(isfinite(AR) && AR>=0.0 && AR < MaxFloat);
	      }
	    }

	    maxLR = max(AR,maxLR);
	  }
	}

	if(DEBUG) assert(isfinite(maxLR) && maxLR >= ((MFLOAT)0.0) && maxLR < MaxFloat);
	if(!SCALE_MONOTONIC || maxLR > MaxLR){ /* need to rescale A[I][*][*].AR so largest value (maxLR) is <= MaxLR and update LogScaleR[I] and ScaleRI[G=I+1 .. min(IMAX,I+DELTA_Y)] */
	  MFLOAT2 rescale = ((maxLR > 1.0 || maxLR > (MaxLR/MaxFloat)) ? max(1.00001f * MinFloat, 0.99999f * MaxLR/maxLR) : 0.99999f * MaxFloat);/* multiplier for A[I][*][*].AR  and ScaleRI[I+1..IMAX] */
	  if(DEBUG) assert(MinFloat < rescale && rescale < MaxFloat);
	  double MaxScaleI = DBL_MAX;
	  if(!SCALE_MONOTONIC && rescale > 1.0){
	    MFLOAT2 origrescale = rescale;

	    for(int G= min(IMAX,I + A.kmax + deltaExtYRef); G > I; G--){
	      if(rescale <= 1.0)
		break;
	      if(ScaleRI[G] >= MaxScaleI/rescale){
		MFLOAT2 rescale2 = rescale;
		//		rescale = (long double)(0.99999 * MaxScaleI) / (long double) ScaleRI[G];
		rescale = 0.99999 * MaxScaleI / ScaleRI[G];
		if(DEBUG>=1+RELEASE && !(MinFloat < rescale && rescale < MaxFloat)){
		  printf("\nI=%d,IMAX=%d,A.kmax=%d,deltaExtYRef=%d,G=%d,maxLR=%0.17e,MaxLR=%0.17e,MinFloat=%0.6e,MaxFloat=%0.6e,MaxScaleI=%0.6e,ScaleRI[G]=%0.17e:rescale= %0.17e -> %0.17e -> %0.17e\n",
			 I,IMAX,G,A.kmax,deltaExtYRef,maxLR,MaxLR,MinFloat,MaxFloat,MaxScaleI,ScaleRI[G],origrescale,rescale2,rescale);
		  fflush(stdout);
		  assert(MinFloat < rescale && rescale < MaxFloat);
		}
	      }
	    }
	    if(rescale > 1.0){/* make sure LogScaleR[I] does not drop below 0.0 */
	      if(LogScaleR[I] - log(rescale) < 0.0){
		MFLOAT2 rescale2 = rescale;
		rescale = exp(LogScaleR[I]);
		if(DEBUG>=1+RELEASE && !(MinFloat < rescale && rescale < MaxFloat)){
		  printf("\nI=%d,maxLR=%0.17e,MaxLR=%0.17e,MinFloat=%0.6e,MaxFloat=%0.6e,LogScaleR[I]=%0.17f:rescale= %0.17e -> %0.17e -> %0.17e\n",
			 I,maxLR,MaxLR,MinFloat,MaxFloat,LogScaleR[I],origrescale,rescale2,rescale);
		  fflush(stdout);
		  assert(MinFloat < rescale && rescale < MaxFloat);
		}
	      }
	    }
	    if(DEBUG && !(MinFloat < rescale && rescale < MaxFloat)){
	      printf("\nI=%d,maxLR=%0.17e,MaxLR=%0.17e,MinFloat=%0.6e,MaxFloat=%0.6e:rescale= %0.17e -> %0.17e\n",
		     I,maxLR,MaxLR,MinFloat,MaxFloat,origrescale,rescale);
	      fflush(stdout);
	      assert(MinFloat < rescale && rescale < MaxFloat);
	    }
	  }
	  for(int K= 0; K <= Kmax[I]; K++)
	    for(int J = Jmin[I]; J <= Jmax[I]; J++)
	      A.AR(I,K,J) *= rescale;
	  double origLogScale = LogScaleR[I];
	  LogScaleR[I] -= log(rescale);
	  if(DEBUG) assert(LogScaleR[I] >= -1e-5);
	  LogScaleR[I] = max(0.0, LogScaleR[I]);
	  for(int G= min(IMAX,I + A.kmax + deltaExtYRef); G > I; G--){
	    ScaleRI[G] *= rescale;
	    if(DEBUG) assert(isfinite(ScaleRI[G]) && ScaleRI[G] <= (SCALE_MONOTONIC ? 1.0 : MaxScaleI));
	  }
	  if(VERB>=3 && m==54 && MD==357 && IMIN==1 && IMAX==147 && addcnt==407){
	    printf("I=%d,LogScaleR[I]= %0.8f -> %0.8f, maxLR=%0.8e->%0.8e\n",I,origLogScale,LogScaleR[I],maxLR,MaxLR);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(ScaleRI[I] == 1.0);
	}

	if(VERB>=3 && MDEBUG && m==MDEBUG_M && (N==MDEBUG_N+1 || N==MDEBUG_N)){
	  printf("m=%d:I=%d:LogScaleR[I]= %0.10f\n",m,I,LogScaleR[I]);
	  if(I<=4)
	    for(int K= 0; K <= Kmax[I]; K++)
	      for(int J = Jmin[I]; J <= Jmax[I]; J++)
		printf("\t I=%d,K=%d,J=%d: AR(I,K,J)= %0.10e, UL(I,K,J)= %0.10e\n", I,K,J, A.AR(I,K,J), A.UL(I,K,J));
	}

	if(VERB>=3 && rverb && m==553 && MD==585 /* && origIMIN==430 && origIMAX==538 */ && addcnt==37){
	  for(int K= 0; K <= Kmax[I]; K++)
	    for(int J = Jmin[I]; J <= Jmax[I]; J++)
	      printf("\t I=%d,K=%d,J=%d: AR(I,K,J)= %0.10e, UL(I,K,J)= %0.10e, LogScaleR[I]= %0.10f\n", I,K,J, A.AR(I,K,J), A.UL(I,K,J), LogScaleR[I]);
	  fflush(stdout);
	}
      }

      LogScaleR[IMIN-1] = LogScaleR[IMIN];

      // If EPS2 > 0 : Narrow Jmin[I]..Jmax[I] ranges based on AL,AR values that are negligible compared to the largest value
      // Adjust so for each I, Sum(J=Jmin[I]..Jmax[I]) Sum(K= 0..Kmax[I]) AL(I,K,J)*PM(I,K)*AR(I,K,J) is not reduced by more than EPS2 * original value :
      // A. For each I = IMIN .. IMAX :
      //     1. Precompute SumIJ[J = Jmin[I] .. Jmax[I]] = Sum(K= 0..Kmax[I]) AL(I,K,J) * AR(I,K,J)
      //                   SumIJ2[J = Jmin[I] .. Jmax[I]] = ScaleL(I) * Sum(K= 0..Kmax[I]) AL(I,K,J) * UR(I,K,J), where ScaleL(I) == exp(LogScaleL[I] - LS)
      //                   SumIJ3[J = Jmin[I] .. Jmax[I]] = ScaleR(I) * Sum(K= 0..Kmax[I]) AR(I,K,J) * UL(I,K,J), where ScaleR(I) == exp(LogScaleR[I] - LS)
      //      1b. Compute LRcumInvLI = ScaleL(I)/LRcum, LRcumInvRI = ScaleR(I)/LRcum
      //           NOTE : In the code all SumIJ[J], SumIJ2[J] values have NOT been scaled by ScaleL(I) and hence LRcumInvLI is set to ScaleL(I)/LRcum
      //                  Similarly SumIJ3[J] values have NOT been scaled by ScaleR(I) and hence LRcumInvRI is set to ScaleR(I)/LRcum
      //     2a. Compute sumI = sum(J=Jmin[I]..Jmax[I])  SumIJ[J]
      //     3. Keep trying to reduce J = Jmax[I] or increase J = Jmin[I], as long as the corresponding term SumIJ[J]/SumI < EPS2 - Sum(SumIJ[J] values discarded previously)/SumI
      //                                                                                               AND   SumIJ2[J] * LRcumInvLI < EPS2 - Sum(SumIJ2[J]*LRcumInvLI, J discarded previously)
      //                                                                                               AND   SumIJ3[J] * LRcumInvRI < EPS2 - Sum(SumIJ3[J]*LRcumInvRI, J discarded previously)
      // B. Finally re-enforce monotonicity of Jmin[I] .. Jmax[I] as a function of I, by expanding some ranges Jmin[I]..Jmax[I]
      //
      // Step B means that step A3 should not make the choice (either reduce Jmax[I] or increase Jmin[I]) in a greedy fashion, one choice may be undone in step B and provide no benefit 
      // A semi-optimal way would be to pick the choice with the sure benefit over the one with the unclear benefit even if a larger share of the EPS error allowance is used up
      // In some cases the benefit may be more than one.
      // To reduce the unclear benefit cases, the entire algorithm can be applied in two passes using EPS*0.1 then EPS2 as the error allowance.
      // The total relative error may be larger than EPS2, since the errors introduced at each I may be partly additive, but typically will not be.

      double LRcum = 0.0, LS = max(0.0,LogScaleL[IMIN] - logMaxScale);

      double startLRcum= 0.0, startLS= max(0.0, LogScaleL[IMIN] - logMaxScale);/* If (EPS2 > 0) startLRcum,startLS will be LRcum,LS before applying EPS2,EPS3 */

      int origIMIN = IMIN;
      int origIMAX = IMAX;

      if(EPS2 > 0.0){
	memcpy(&origJmin[IMIN],&Jmin[IMIN],(IMAX-IMIN+1)*sizeof(int));
	memcpy(&origJmax[IMIN],&Jmax[IMIN],(IMAX-IMIN+1)*sizeof(int));
      }

      if(EPS2 > 0.0){
	if(DEBUG>=2){
	  for(int I = IMIN; I <= IMAX; I++){
	    if(!(1 <= origJmin[I] && origJmin[I] <= origJmax[I] && origJmax[I] <= M)){
              #pragma omp critical
	      {
	        printf("m=%d:I=%d:origJmin[I]=%d,origJmax[I]=%d,M=%d,origIMIN=%d,origIMAX=%d\n",
		  m,I,origJmin[I],origJmax[I],M,origIMIN,origIMAX);
		for(int i = origIMIN; i <= origIMAX; i++)
		  printf("\t i=%d:origJmin[i]=%d,origJmax[i]=%d\n",i,origJmin[i],origJmax[i]);
		fflush(stdout);
		assert(1 <= origJmin[I] && origJmin[I] <= origJmax[I] && origJmax[I] <= M);
	      }
	    }
	  }
	}

	/* compute original LRcum without EPS2 */
	/* locate total alignment LR based on AL & UR : only check results where J >= M+1-DELTA_XU (OR EXTEND && I >= N+1-DELTA_YU && J >= M+1 - DELTA_XU - EXTEND_MAX)
	   unless PoutlierEnd > 0.0, in which case check end outliers for all cases (within the bounds of IMIN,IMAX,Jmin[],Jmax[]) */
	LS = max(0.0,LogScaleL[IMIN] - logMaxScale);
	LRcum = 0.0;

	for(int I = origIMIN; I <= origIMAX; I++){
	  double *SumIJ = &MSumIJ[(I-origIMIN) * M - 1];
	  double *SumIJ2 = &MSumIJ2[(I-origIMIN) * M - 1];
	  double *SumIJ3 = &MSumIJ3[(I-origIMIN) * M - 1];

	  double ScaleL = exp(LogScaleL[I] - LS);
	  if(DEBUG) assert(isfinite(ScaleL) && ScaleL >= 0.0 && ScaleL < MaxScale * ((MFLOAT)1.01));

	  double LRcumI = 0.0;
	  int jmin = (PoutlierEnd > 0.0) ? Jmin[I] : (EXTEND && I >= N+1 - DELTA_YU) ? max(Jmin[I],M+1 - DELTA_XU - EXTEND_MAX) : max(Jmin[I],M+1 - DELTA_XU);
	  int jmax = Jmax[I];
	  for(int J = Jmin[I]; J <= jmax; J++)
	    SumIJ3[J] = SumIJ2[J] = SumIJ[J] = 0.0;
	  for(int J = jmin; J <= jmax; J++){
	    for(int K = 0; K <= Kmax[I]; K++){
	      // NOTE : must use double, since AL*AR or AL*UR may underflow without ScaleL
	      double AL = A.AL(I,K,J);
	      double AR = A.AR(I,K,J);
	      double UR = A.UR(I,K,J);
	      double UL = A.UL(I,K,J);
	      double PMik = PRtab[K][I].PM;
	      if(DEBUG>=2) assert(isfinite(AL) && AL >= 0.0);
	      if(DEBUG>=2) assert(isfinite(AR) && AR >= 0.0);
	      if(DEBUG>=2) assert(isfinite(UR) && UR >= 0.0);
	      if(DEBUG>=2) assert(isfinite(UL) && UL >= 0.0);
	      if(DEBUG>=2 && !(isfinite(PMik) && PMik > 0.0)){
		printf("m=%d,M=%d,N=%d:I=%d,K=%d,J=%d:AL=%0.10e,AR=%0.10e,PMik=%0.10e (Y[I]=%0.8f,Y[I-K]=%0.8f)\n",
		       m,M,N,I,K,J,AL,AR,PMik,Y[I],Y[I-K]);
		fflush(stdout);
		assert(isfinite(PMik) && PMik > 0.0);
	      }
	      SumIJ[J] += AL * AR / PMik; /* Note that AL,AR already include PM(J,I,K,Y) */
	      SumIJ2[J] += AL * UR; 
	      SumIJ3[J] += AR * UL;
	    }
	    LRcumI += SumIJ2[J];
	    if(DEBUG) assert(isfinite(LRcumI));
	  }
	  LRcum += LRcumI * ScaleL;
	  if(DEBUG) assert(isfinite(LRcum));
	  
	  double origLRcum = LRcum;
	  double origLS = LS;
	  double origScaleL = ScaleL;
	  if(LRcum > LRCUM_MAX /* WAS && LS+log(LRcum) > LogScaleL[I+1] */){/* update scaling factor LS = LS+log(LRcum) */
	    LS += log(LRcum);
	    LRcum = 1.0;
	  }
	  if(LS < LogScaleL[I+1] - logMaxScale + 0.001){
	    double pLS = LS;
	    LS = LogScaleL[I+1] - logMaxScale + 0.001;
	    double Scale = exp(pLS - LS);
	    if(DEBUG && LRcum > DBL_MIN && !(Scale > 0.0 && LRcum * Scale > DBL_MIN)){
	      printf("I=%d:pLS=%0.8f,LS=%0.8f,LogScaleL[I]=%0.8f,LogScaleL[I+1]=%0.8f,logMaxScale=%0.8f:Scale = %0.8e, LRcum=%0.8e\n",
		     I,pLS,LS,LogScaleL[I],LogScaleL[I+1],logMaxScale,Scale,LRcum);
	      fflush(stdout);
	      assert(Scale > 0.0 && LRcum * Scale > DBL_MIN);
	    }
	    LRcum *= Scale;
	  }
	  if(LRCUM_FIX && LRcum > 0.0 && LRcum <= MinFloat && LS > LogScaleL[I+1] - logMaxScale + 0.001){
	    double pLS = LS;
	    LS = LogScaleL[I+1] - logMaxScale + 0.001;
	    double Scale = exp(pLS - LS);
	    if(LRcum * Scale > 1.0){
	      Scale = 1.0 / LRcum;
	      LS = pLS - log(Scale);
	    }
	    if(DEBUG && !(Scale >= 1.0 && LRcum * Scale > MinFloat && LRcum * Scale <= MaxFloat)){
	      printf("I=%d:pLS=%0.8f,LS=%0.8f,LogScaleL[I]=%0.8f,LogScaleL[I+1]=%0.8f,logMaxScale=%0.8f:Scale = %0.8e, LRcum=%0.8e\n",
		     I,pLS,LS,LogScaleL[I],LogScaleL[I+1],logMaxScale,Scale,LRcum);
	      fflush(stdout);
	      assert(Scale >= 1.0 && LRcum * Scale > MinFloat && LRcum * Scale <= MaxFloat);
	    }
	    LRcum *= Scale;
	  }
	  if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792&& origIMIN==430 && origIMAX==538 && addcnt==14149 */){
	    printf("m=%d,I=%d:j=%d..%d:LRcumI=%0.8e,LRcum=%0.8e->%0.8e,ScaleL=%0.8e,LS=%0.8f->%0.8f,log(LRcum)+LS=%0.8f,LogScaleL[I+1]=%0.8f,LogMaxScale=%0.8f(A)\n",
		  m,I,jmin,jmax,LRcumI,origLRcum,LRcum,origScaleL,origLS,LS,(LRcum > 0.0) ? log(LRcum)+LS : 0.0, LogScaleL[I+1], logMaxScale);
	    fflush(stdout);
	  }

	  if(DEBUG && !(isfinite(LS) && isfinite(LRcum))) assert((LRcum > MinFloat || origLRcum <= DBL_MIN));
	}

	int error = 0;

	/* Also adjust LS so ScaleR(I) == exp(LogScaleR[I] - LS) < MaxScale */
	for(int I = IMIN; I <= IMAX; I++){
	  if(LS < LogScaleR[I] - logMaxScale + 0.001){
	    double pLRcum = LRcum;
	    double pLS = LS;
	    LS = LogScaleR[I] - logMaxScale + 0.001;
	    double ScaleR = exp(pLS - LS);

	    if(DEBUG && !SCALE_MONOTONIC && LRcum > DBL_MIN && !(ScaleR > 0.0 && LRcum * ScaleR > DBL_MIN)){
	      printf("\n m=%d:IMIN=%d,%d,IMAX=%d,%d,I=%d:LogScaleR[I]=%0.8f,logMaxScale=%0.8f,pLS=%0.8f,LS=%0.8f,ScaleR=%0.8e,LRcum=%0.8e,MaxFloat=%0.8e,PoutlierEnd=%0.4e\n",
		     m,origIMIN,IMIN,IMAX,origIMAX,I,LogScaleR[I],logMaxScale,pLS,LS,ScaleR,LRcum,MaxFloat,PoutlierEnd);
	      for(int i = origIMIN; i <= origIMAX; i++)
		printf("I=%d:LogScaleL[I]=%0.8f,LogScaleR[I]=%0.8f\n",i,LogScaleL[i],LogScaleR[i]);
	      if(DEBUG>=1+RELEASE){// Extra debug code : compute LRcum using AR and UL : this LRcum should not have underflowed
		double nLRcum = LRcum, nLS = LS;
		LS = max(0.0,LogScaleR[origIMAX] - logMaxScale);
		LRcum = 0.0;

		for(int I = origIMAX; I >= origIMIN; I--){
		  double ScaleR = exp(LogScaleR[I] - LS);
		  double LRcumI = 0.0;
		  for(int K = 0; K <= Kmax[I]; K++){
		    int jmin = origJmin[I];
		    int jmax = (PoutlierEnd > 0.0) ? origJmax[I] : (EXTEND && I-K <= DELTA_YU) ? min(origJmax[I], DELTA_XU + EXTEND_MAX) : min(origJmax[I], DELTA_XU);
		    for(int J = jmin; J <= jmax; J++){
		      LRcumI += (double)A.AR(I,K,J) * (double)A.UL(I,K,J); /* Note that AR already includes PM */
		      printf("m=%d,I=%d,K=%d,J=%d: AR(I,K,J)= %0.10e, UL(I,K,J)= %0.10e, AL(I,K,J)= %0.10e, UR(I,K,J)= %0.10e : LRcumI=%0.10e\n",
			     m,I,K,J,A.AR(I,K,J),A.UL(I,K,J), A.AL(I,K,J), A.UR(I,K,J), LRcumI);
		    }
		  }
		  LRcum += LRcumI * ScaleR;

		  if(1){  
		    double origLS = LS;
		    double origLRcum = LRcum;
		    if(LRcum > LRCUM_MAX){/* update scaling factor LS = LS+log(LRcum) */
		      LS += log(LRcum);
		      LRcum = 1.0;
		    }
		    if(LS < LogScaleR[I-1] - logMaxScale + 0.001){
		      double Scale = exp(LS - LogScaleR[I-1] + logMaxScale - 0.001);
		      LS = LogScaleR[I-1] - logMaxScale + 0.001;
		      LRcum *= Scale;
		    }
		    printf("m=%d,I=%d:j=%d..%d:LRcumI=%0.8e,LRcum=%0.8e->%0.8e,ScaleR=%0.8e,LS=%0.8f->%0.8f,LogScaleR[I-1]=%0.8f,logMaxScale=%0.8f(D)\n",
			   m,I,jmin,jmax,LRcumI,origLRcum,LRcum,ScaleR,origLS,LS,LogScaleR[I-1],logMaxScale);
		  }
		}

		LRcum = nLRcum;
		LS = nLS;
	      }// extra debug code

	      printf("\nWARNING: m=%d(id=%lld):IMIN=%d,%d,IMAX=%d,%d,I=%d:LogScaleR[I]=%0.8f,logMaxScale=%0.8f,pLS=%0.8f,LS=%0.8f,ScaleR=%0.8e,LRcum=%0.8e,MaxFloat=%0.8e,PoutlierEnd=%0.4e\n",
		     m,gmap[pcontig->contig[m].mapid]->id,origIMIN,IMIN,IMAX,origIMAX,I,LogScaleR[I],logMaxScale,pLS,LS,ScaleR,LRcum,MaxFloat,PoutlierEnd);
	      fflush(stdout);
	      error = 1;
	      assert(ScaleR > 0.0 && LRcum * ScaleR > DBL_MIN);// NOTE this means AL * UR based LRcum underflows compared to AR * UL based LRcum
	      // Since LRcum is double precision and subsequent mprobeval evaluation is based on AR,UL, this is not a fatal error, but may result result in a later assertion failure
	      // since qprobeval is based on AL * UR based LRcum and hence qprobeval() may not agree with mprobeval
	    }
	    LRcum *= ScaleR;

	    if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 && origIMIN==430 && origIMAX==538 && addcnt==14149 */){
	      printf("m=%d,I=%d:j=%d..%d:LRcum=%0.8e->%0.8e,ScaleR=%0.8e,LS=%0.8f->%0.8f,LogScaleR[I]=%0.8f,LogMaxScale=%0.8f(A)\n",
		     m,I,jmin,jmax,pLRcum,LRcum,ScaleR,pLS,LS,LogScaleR[I],logMaxScale);
	      fflush(stdout);
	    }
	  }
	}

	startLRcum = LRcum;
	startLS = LS;

	if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 && origIMIN==430 && origIMAX==538 && addcnt==468 */){
	  printf("m=%d:LRcum=%0.10e,LS=%0.10f(A),error=%d\n",m,LRcum,LS,error);
	  fflush(stdout);
	}

	double LRbiasLS = LRbias * exp(-LS);
	double LRcumInv = (LRcum > 0.0 || LRbiasLS > 0.0) ? 1.0/(LRcum + LRbiasLS) : 0.0;

	if(VERB_EPS2){/* collect some statistics */
	  double LPm = (LRcum > 0.0 || LRbiasLS > 0.0) ? LS+log(LRcum + LRbiasLS) : -99999.0;	  
	  LLRsumEPS2 += LPm * TBwtm;

	  if(VERB>=3 && m==25 && MD==440 && origIMIN==1 && origIMAX==103 && addcnt==882){
            #pragma omp critical
	    {
	      printf("m=%d: original LPm=%0.10f\n",m, LPm);
	      fflush(stdout);
	    }
	  }

	  for(int I = IMIN; I <= IMAX; I++)
	    origIJcnt += Jmax[I]-Jmin[I]+1;
	}

	for(int pass = 0; pass < 2; pass++){
	  double E2 = pass ? EPS2 : EPS2 * 0.1;

	  if(VERB>=3 && m==MDEBUG_M /* && MD==440 && origIMIN==1 && origIMAX==103 && addcnt==882*/){
	    printf("m=%d: Before pass p=%d : IMIN->%d,IMAX=%d\n",m,pass,IMIN,IMAX);
	    for(int I=IMIN; I <= IMAX; I++)
	      printf("\t I=%d:Jmin[I]=%d,Jmax[I]=%d\n",I,Jmin[I],Jmax[I]);
	    fflush(stdout);
	  }

	  for(int I = IMIN; I <= IMAX; I++)
	    E2sum[I] = 0.0;/* keep track of fractional error for LR computed via index I using sum(K) AL * AR / PM sum */
	  double E2sumIJ = 0.0;/* same fractional error for LR computed using sum(I,J,K) (AL*UR) */
          double E3sumIJ = 0.0;/* same fractional error for LR computed using sum(I,J,K) (AR*UL) */

	  double pLS = LS;
	  double pLRcum = LRcum;
	  LRcum = 0.0;
	  for(int I = origIMIN; I <= origIMAX; I++){
	    double *SumIJ = &MSumIJ[(I-origIMIN) * M - 1];
	    double *SumIJ2 = &MSumIJ2[(I-origIMIN) * M - 1];
	    double *SumIJ3 = &MSumIJ3[(I-origIMIN) * M - 1];

	    double ScaleL = exp(LogScaleL[I] - LS);
	    if(DEBUG>=2 && !(isfinite(ScaleL) && ScaleL >= 0.0 && ScaleL < MaxScale * ((MFLOAT)1.01))){
              #pragma omp critical
	      {
		printf("m=%d,I=%d:LogScaleL[I]=%0.8f,LS=%0.8f,ScaleL=%0.8e\n",m,I,LogScaleL[I],LS,ScaleL);
		for(int i = IMIN; i <= IMAX; i++)
		  if(LogScaleL[i] != 0.0)
		    printf("\t LogScaleL[%d]= %0.8f\n",i,LogScaleL[i]);
		fflush(stdout);
		assert(isfinite(ScaleL) && ScaleL >= 0.0 && ScaleL < MaxScale * ((MFLOAT)1.01));
	      }
	    }
	    double ScaleR = exp(LogScaleR[I] - LS);
	    if(DEBUG>=2 && !(isfinite(ScaleR) && ScaleR >= 0.0 && ScaleR < MaxScale * ((MFLOAT)1.01))){
              #pragma omp critical
	      {
		printf("m=%d,I=%d:LogScaleR[I]=%0.8f,LS=%0.8f,ScaleR=%0.8e,pLRcum=%0.8e\n",m,I,LogScaleR[I],LS,ScaleR,pLRcum);
		for(int i = IMIN; i <= IMAX; i++)
		  if(LogScaleR[i] != 0.0)
		    printf("\t LogScaleR[%d]= %0.8f\n",i,LogScaleR[i]);
		fflush(stdout);
		assert(isfinite(ScaleR) && ScaleR >= 0.0 && ScaleR < MaxScale * ((MFLOAT)1.01));
	      }
	    }

	    double LRcumInvLI = ScaleL * LRcumInv;
	    double LRcumInvRI = ScaleR * LRcumInv;
	    double SumI = 0.0;

	    double LRcumI = 0.0;
	    int jmin = (PoutlierEnd > 0.0) ? origJmin[I] : (EXTEND && I >= N+1 - DELTA_YU) ? max(origJmin[I],M+1 - DELTA_XU - EXTEND_MAX) : max(origJmin[I],M+1 - DELTA_XU);
	    int jmax = origJmax[I];
	    for(int J = jmin; J <= jmax; J++){
	      LRcumI += SumIJ2[J];
	      SumI += SumIJ[J];
	    }
	    LRcum += LRcumI * ScaleL;

	    double origLRcum = LRcum;
	    double origLS = LS;
	    if(DEBUG>=2) assert(LRcum <= LRCUM_MAX * 1.0001);
	    if(LRcum > LRCUM_MAX /* WAS && LS+log(LRcum) > LogScaleL[I+1] */){/* update scaling factor LS = LS+log(LRcum) */
	      LS += log(LRcum);
	      LRcum = 1.0;
	    } 
	    if(LS < LogScaleL[I+1] - logMaxScale + 0.001){
	      ScaleL = exp(LS - LogScaleL[I+1] + logMaxScale - 0.001);
	      LS = LogScaleL[I+1] - logMaxScale + 0.001;
	      LRcum *= ScaleL;
	    }
	    if(LRCUM_FIX && LRcum > 0.0 && LRcum <= MinFloat && LS > LogScaleL[I+1] - logMaxScale + 0.001){
	      double pLS = LS;
	      LS = LogScaleL[I+1] - logMaxScale + 0.001;
	      double Scale = exp(pLS - LS);
	      if(LRcum * Scale > 1.0){
		Scale = 1.0 / LRcum;
		LS = pLS - log(Scale);
	      }
	      LRcum *= Scale;
	    }
	    if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 && origIMIN==430 && origIMAX==538 && addcnt==14149 */){
	      printf("m=%d,pass=%d,I=%d:j=%d..%d:LRcumI=%0.8e,SumI=%0.8e,LRcum=%0.8e->%0.8e,ScaleL=%0.8e,LS=%0.8f->%0.8f,log(LRcum)+LS=%0.8f,LogScaleL[I+1]=%0.8f\n",
		  m,pass,I,jmin,jmax,LRcumI,SumI,origLRcum,LRcum,ScaleL,origLS,LS,(LRcum > 0.0) ? log(LRcum)+LS : 0.0, LogScaleL[I+1]);
	      fflush(stdout);
	    }
	    if(DEBUG && !(isfinite(LS) && isfinite(LRcum))) assert((LRcum > MinFloat || origLRcum <= MinFloat));

	    if(pass > 0){ /* update E2sumIJ,E3sumIJ based on J values discarded in previous pass */
	      for(int J = origJmin[I]; J < Jmin[I]; J++){
		E2sumIJ += SumIJ2[J] * LRcumInvLI;
		E3sumIJ += SumIJ3[J] * LRcumInvRI;
	      }
	      for(int J = Jmax[I] + 1; J <= origJmax[I]; J++){
		E2sumIJ += SumIJ2[J] * LRcumInvLI;
		E3sumIJ += SumIJ3[J] * LRcumInvRI;
	      }
	    }

	    if(VERB>=3 && m==MDEBUG_M /* && MD==440 && origIMIN==1 && origIMAX==103 && addcnt==882 */){
	      printf("m=%d,I=%d,origJmin[I]=%d,Jmin[I]=%d,Jmax[I]=%d,origJmax[I]=%d:SumI=%0.8e, E2sumIJ=%0.8e, E3sumIJ=%0.8e,E2=%0.2e:\n",
		     m,I,origJmin[I],Jmin[I],Jmax[I],origJmax[I],SumI,E2sumIJ,E3sumIJ,E2);
	      for(int J = origJmin[I]; J <= origJmax[I]; J++)
		printf("\t J=%d:SumIJ[J]=%0.8e (Ratio= %0.8e), SumIJ2[J]=%0.8e (Ratio= %0.8e), SumIJ3[J]=%0.8e (Ratio= %0.8e)\n",
		       J,SumIJ[J], (SumI > 0.0) ? SumIJ[J]/SumI : 0.0, SumIJ2[J], SumIJ2[J] * LRcumInvLI, SumIJ3[J], SumIJ3[J] * LRcumInvRI);
	      fflush(stdout);
	    }

	    if(SumI <= 0.0){/* trivial case */
	      double origE2sumIJ = E2sumIJ;
	      double origE3sumIJ = E3sumIJ;
	      for(int J = Jmin[I]; J <= Jmax[I]; J++){
		E2sumIJ += SumIJ2[J] * LRcumInvLI;
		E3sumIJ += SumIJ3[J] * LRcumInvRI;
	      }
	      if(E2sumIJ < E2 && E3sumIJ < E2){
		Jmin[I] = M;
		Jmax[I] = 1;
		continue;
	      }

	      /* try again more cautiously */
	      E2sumIJ = origE2sumIJ;
	      E3sumIJ = origE3sumIJ;
	      while(Jmin[I] <= Jmax[I]){
		double minRatio2 = SumIJ2[Jmin[I]] * LRcumInvLI;
		double maxRatio2 = SumIJ2[Jmax[I]] * LRcumInvLI;
		double minRatio3 = SumIJ3[Jmin[I]] * LRcumInvRI;
		double maxRatio3 = SumIJ3[Jmax[I]] * LRcumInvRI;
		
		double minBenefit = ((I < IMAX && Jmin[I] >= Jmin[I+1]) || (I+1 < IMAX && Jmin[I] >= Jmin[I+2])) ? 0.1 :
		  (I > IMIN && Jmin[I] < Jmin[I-1]) ? (I-1 > IMIN && Jmin[I] < Jmin[I-2] ? 3.0 : 2.0) : 1.0;
		double maxBenefit = ((I > IMIN && Jmax[I] <= Jmax[I-1]) || (I-1 > IMIN && Jmax[I] <= Jmax[I-2])) ? 0.1 :
		  (I < IMAX && Jmax[I] > Jmax[I+1]) ? (I+1 < IMAX && Jmax[I] > Jmax[I+2] ? 3.0 : 2.0) : 1.0;

		if(VERB>=3 && m==MDEBUG_M /* && MD==440 && origIMIN==1 && origIMAX==103 && addcnt==882 */){
		  printf("m=%d,p=%d,I=%d:Jmin[I]=%d,Jmax[I]=%d;minR=%0.6e,%0.6e,maxR=%0.6e,%0.6e,E2=%0.3e,E2sumIJ=%0.6e,E3sumIJ=%0.6e:B=%0.1f,%0.1f\n",
			 m,pass,I,Jmin[I],Jmax[I],minRatio2,minRatio3,maxRatio2,maxRatio3,E2,E2sumIJ,E3sumIJ,minBenefit,maxBenefit);
		  fflush(stdout);
		}

		if(pass==0 ? (maxRatio2 > minRatio2) : (maxRatio2 * minBenefit > minRatio2 * maxBenefit)){
		  if(minRatio2 + E2sumIJ <= E2 && minRatio3 + E3sumIJ <= E2){
		    Jmin[I]++;
		    E2sumIJ += minRatio2;
		    E3sumIJ += minRatio3;
		    continue;
		  }
		  if(maxRatio2 + E2sumIJ <= E2 && maxRatio3 + E3sumIJ <= E2){
		    Jmax[I]--;
		    E2sumIJ += maxRatio2;
		    E3sumIJ += maxRatio3;
		    continue;
		  }
		} else {
		  if(maxRatio2 + E2sumIJ <= E2 && maxRatio3 + E3sumIJ <= E2){
		    Jmax[I]--;
		    E2sumIJ += maxRatio2;
		    E3sumIJ += maxRatio3;
		    continue;
		  }
		  if(minRatio2 + E2sumIJ <= E2 && minRatio3 + E3sumIJ <= E2){
		    Jmin[I]++;
		    E2sumIJ += minRatio2;
		    E3sumIJ += minRatio3;
		    continue;
		  }
		}
		break;
	      }

	      continue;
	    }
	  
	    double SumInv = 1.0/SumI;
	    if(pass > 0){ /* update E2sum[I] based on J values discarded in previous pass */
	      for(int J = origJmin[I]; J < Jmin[I]; J++)
		E2sum[I] += SumIJ[J] * SumInv;
	      for(int J = Jmax[I] + 1; J <= origJmax[I]; J++)
		E2sum[I] += SumIJ[J] * SumInv;
	    }
	    if(VERB>=3 && m==MDEBUG_M /* && MD==440 && origIMIN==1 && origIMAX==103 && addcnt==882 */){
	      printf("m=%d,I=%d,origJmin[I]=%d,Jmin[I]=%d,Jmax[I]=%d,origJmax[I]=%d:SumI=%0.8e, E2sumIJ=%0.8e, E2sum[I]=%0.8e:\n",
		     m,I,origJmin[I],Jmin[I],Jmax[I],origJmax[I],SumI,E2sumIJ,E2sum[I]);	      
	      fflush(stdout);
	    }

	    while(Jmin[I] <= Jmax[I]){
	      double minRatio = SumIJ[Jmin[I]] * SumInv;
	      double maxRatio = SumIJ[Jmax[I]] * SumInv;
	      double minRatio2 = SumIJ2[Jmin[I]] * LRcumInvLI;
	      double maxRatio2 = SumIJ2[Jmax[I]] * LRcumInvLI;
	      double minRatio3 = SumIJ3[Jmin[I]] * LRcumInvRI;
	      double maxRatio3 = SumIJ3[Jmax[I]] * LRcumInvRI;

	      double minBenefit = ((I < IMAX && Jmin[I] >= Jmin[I+1]) || (I+1 < IMAX && Jmin[I] >= Jmin[I+2])) ? 0.1 :
		(I > IMIN && Jmin[I] < Jmin[I-1]) ? (I-1 > IMIN && Jmin[I] < Jmin[I-2] ? 3.0 : 2.0) : 1.0;
	      double maxBenefit = ((I > IMIN && Jmax[I] <= Jmax[I-1]) || (I-1 > IMIN && Jmax[I] <= Jmax[I-2])) ? 0.1 :
		(I < IMAX && Jmax[I] > Jmax[I+1]) ? (I+1 < IMAX && Jmax[I] > Jmax[I+2] ? 3.0 : 2.0) : 1.0;

	      if(VERB>=3 && m==MDEBUG_M /* && MD==440 && origIMIN==1 && origIMAX==103 && addcnt==882 */){
		printf("m=%d,p=%d,I=%d:SumIJ[Jmin[I]=%d]=%0.6e,SumIJ[Jmax[I]=%d]=%0.6e;minR=%0.6e(%0.6e,%0.6e),maxR=%0.6e(%0.6e,%0.6e),E2=%0.6e,E2sum[I]=%0.6e(%0.6e,%0.6e):B=%0.1f,%0.1f\n",
		       m,pass,I,Jmin[I],SumIJ[Jmin[I]],Jmax[I],SumIJ[Jmax[I]],minRatio,minRatio2,minRatio3,maxRatio,maxRatio2,maxRatio3,E2,E2sum[I],E2sumIJ,E3sumIJ,minBenefit,maxBenefit);
		fflush(stdout);
	      }

	      if(pass==0 ? (maxRatio > minRatio) : (maxRatio * minBenefit > minRatio * maxBenefit)){
		if(minRatio + E2sum[I] < E2 && minRatio2 + E2sumIJ <= E2 && minRatio3 + E2sumIJ <= E2){
		  Jmin[I]++;
		  E2sum[I] += minRatio;
		  E2sumIJ += minRatio2;
		  E3sumIJ += minRatio3;
		  continue;
		}
		if(maxRatio + E2sum[I] < E2 && maxRatio2 + E2sumIJ <= E2 && maxRatio3 + E3sumIJ <= E2){
		  Jmax[I]--;
		  E2sum[I] += maxRatio;
		  E2sumIJ += maxRatio2;
		  E3sumIJ += maxRatio3;
		  continue;
		}
	      } else {
		if(maxRatio + E2sum[I] < E2 && maxRatio2 + E2sumIJ <= E2 && maxRatio3 + E3sumIJ <= E2){
		  Jmax[I]--;
		  E2sum[I] += maxRatio;
		  E2sumIJ += maxRatio2;
		  E3sumIJ += maxRatio3;
		  continue;
		}
		if(minRatio + E2sum[I] < E2 && minRatio2 + E2sumIJ <= E2 && minRatio3 + E3sumIJ <= E2){
		  Jmin[I]++;
		  E2sum[I] += minRatio;
		  E2sumIJ += minRatio2;
		  E3sumIJ += minRatio3;
		  continue;
		}
	      }
	      break;
	    }
	    if(VERB>=3 && m==MDEBUG_M /* && MD==440 && origIMIN==1 && origIMAX==103 && addcnt==882 */){
	      printf("m=%d,p=%d,I=%d:Jmin[I]=%d -> %d,Jmax[I]=%d -> %d, E2sum[I]= %0.6e\n",m,pass,I,origJmin[I],Jmin[I],origJmax[I],Jmax[I],E2sum[I]);
	      fflush(stdout);
	    }
	  }// I = origIMIN .. origIMAX
	  if(DEBUG && pLRcum > 0.0) assert(LRcum > 0.0);
	  if(DEBUG && pLRcum > 0.0 && !(fabs(log(LRcum) + LS - log(pLRcum) - pLS) <= fabs(log(LRcum)+LS) * 1e-8 + 1e-8)){
	    printf("m=%d,p=%d:LRcum= %0.8e, pLRcum = %0.8e, LS= %0.8f, pLS= %0.8f\n",m,pass,LRcum, pLRcum, LS, pLS);
	    fflush(stdout);
	    assert(fabs(log(LRcum) + LS - log(pLRcum) - pLS) <= fabs(log(LRcum)+LS) * 1e-8 + 1e-8);
	  //	    assert(fabs(LRcum - pLRcum) <= 1e-8 * fabs(LRcum));
	  }
	  LRcum = pLRcum;
	  LS = pLS;

	  /* re-enforce monotonicity of Jmin[I] .. Jmax[I] and Jmin[I] <= Jmax[I]
	     1. Jmin[I=IMAX..IMIN] is a non-increasing function as I decreases.
	     2. Jmax[I=IMIN..IMAX] is a non-decreasing function as I increases.
	  */
	  while(IMIN < IMAX){
	    if(Jmin[IMIN] > Jmax[IMIN]){
	      if(VERB>=3 && m==25 && MD==440 && origIMIN==1 && origIMAX==103 && addcnt==882){
		printf("m=%d,IMIN=%d:Jmin[IMIN]=%d,Jmax[IMIN]=%d:IMIN=%d -> %d\n",m,IMIN,Jmin[IMIN],Jmax[IMIN],IMIN,IMIN+1);
		fflush(stdout);
	      }
	      IMIN++;
	      continue;
	    }
	    if(Jmin[IMAX] > Jmax[IMAX]){
	      if(VERB>=3 && m==25 && MD==440 && origIMIN==1 && origIMAX==103 && addcnt==882){
		printf("m=%d,IMAX=%d:Jmin[IMAX]=%d,Jmax[IMAX]=%d:IMAX=%d -> %d\n",m,IMAX,Jmin[IMAX],Jmax[IMAX],IMAX,IMAX+1);
		fflush(stdout);
	      }
	      IMAX--;
	      continue;
	    }
	    break;
	  }
	  if(DEBUG) assert(IMIN <= IMAX);

	  for(int I = IMIN; I <= IMAX; I++){
	    if(Jmin[I] > origJmax[I])
	      Jmin[I] = origJmax[I];
	    if(Jmax[I] < origJmin[I])
	      Jmax[I] = origJmin[I];
	    if(Jmin[I] > Jmax[I])
	      Jmax[I] = Jmin[I] = (Jmax[I] + Jmin[I])/2;
	    if(DEBUG>=2) assert(origJmin[I] <= Jmin[I] && Jmin[I] <= Jmax[I] && Jmax[I] <= origJmax[I]);
	  }

	  for(int I = IMAX; --I >= IMIN;)
	    Jmin[I] = min(Jmin[I],Jmin[I+1]);
	  for(int I = IMIN; ++I <= IMAX;)
	    Jmax[I] = max(Jmax[I], Jmax[I-1]);

	  if(VERB>=3 && m==25 && MD==440 && origIMIN==1 && origIMAX==103 && addcnt==882){
	    printf("m=%d,p=%d: After enforcing monotonicity : IMIN->%d,IMAX=%d\n",m,pass,IMIN,IMAX);
	    for(int I=IMIN; I <= IMAX; I++)
	      printf("\t I=%d:Jmin[I]=%d,Jmax[I]=%d\n",I,Jmin[I],Jmax[I]);
	    fflush(stdout);
	  }

	  if(DEBUG>=2 && !(IMIN >= IMAX))
	    for(int I = IMIN; I <= IMAX; I++){
	      if(!(1 <= origJmin[I] && origJmin[I] <= Jmin[I] && (/* WAS I > IMIN || I < IMAX || */Jmin[I] <= Jmax[I]) && Jmax[I] <= origJmax[I] && origJmax[I] <= M)){
                #pragma omp critical
		{
		  printf("m=%d,pass=%d:I=%d:origJmin[I]=%d,Jmin[I]=%d,Jmax[I]=%d,origJmax[I]=%d,M=%d,origIMIN=%d,IMIN=%d,IMAX=%d,origIMAX-=%d\n",
			 m,pass,I,origJmin[I],Jmin[I],Jmax[I],origJmax[I],M,origIMIN,IMIN,IMAX,origIMAX);
		  for(int i = origIMIN; i <= origIMAX; i++)
		    printf("\t i=%d:origJmin[i]=%d,Jmin[i]=%d,Jmax[i]=%d,origJmax[i]=%d\n",i,origJmin[i],Jmin[i],Jmax[i],origJmax[i]);
		  fflush(stdout);
		  assert(1 <= origJmin[I] && origJmin[I] <= Jmin[I] && (/* WAS I > IMIN || I < IMAX || */ Jmin[I] <= Jmax[I]) && Jmax[I] <= origJmax[I] && origJmax[I] <= M);
		}
	      }
	    }

	  /* re-compute Imin[1..M],Imax[1..M] */
	  for(int J = 1; J <= M; J++){
	    Imin[J] = N+1;
	    Imax[J] = 0;
	  }
	  for(int I = IMIN; I <= IMAX; I++){
	    for(int J = Jmin[I]; J <= Jmax[I]; J++){
	      Imin[J] = min(Imin[J],I);
	      Imax[J] = max(Imax[J],I);
	    }
	  }    
	  if(DEBUG>=2)
	    for(int J = 1; J <= M; J++)
	      for(int I = Imin[J]; I <= Imax[J]; I++)
		assert(Jmin[I] <= J && J <= Jmax[I]);
	}

	if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 && origIMIN==430 && origIMAX==538 && addcnt==468 */){
	  printf("m=%d:LRcum=%0.10e,LS=%0.10f(A),error=%d\n",m,LRcum,LS,error);
	  fflush(stdout);
	}

	if(VERB_EPS2){/* collect some statistics */
	  int origIJcnt = IJcnt;
	  for(int I = IMIN; I <= IMAX; I++)
	    IJcnt += Jmax[I]-Jmin[I]+1;
	  if(VERB>=2){
	    printf("m=%d:IJcnt=%d\n",m,IJcnt-origIJcnt);
	    fflush(stdout);
	  }
	}
      } // EPS2 > 0.0

      if(DEBUG>=2){
        for(int I = IMIN; I <= IMAX; I++){
          for(int K = 0; K <= Kmax[I]; K++){
	    for(int J = Jmin[I]; J <= Jmax[I]; J++){
	      if(J <= DELTA_XU || (EXTEND && I-K <= DELTA_YU && J <= DELTA_XU + EXTEND_MAX)){
		int ULmin = A.ULmin(I,K,J);
		int ULmax = A.ULmax(I,K,J);
		if(DEBUG && !(ULmin == I-K && ULmax == I-K-1) && !(0 <= ULmin && ULmin <= ULmax && ULmax < I-K)){
		  printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:ULmin=%d,ULmax=%d\n",m,I,K,J,N,M,ULmin,ULmax);
		  fflush(stdout);
		  assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		}
	      }
	      if(J >= M+1-DELTA_XU || (EXTEND && I >= N+1 - DELTA_YU && J >= M+1-DELTA_XU - EXTEND_MAX)){
		int URmin = A.URmin(I,K,J);
		int URmax = A.URmax(I,K,J);
		if(DEBUG && !(URmin == I+1 && URmax == I) && !(I < URmin && URmin <= URmax && URmax <= N+1)){
		  printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:URmin=%d,URmax=%d\n",m,I,K,J,N,M,URmin,URmax);
		  fflush(stdout);
		  assert(I < URmin && URmin <= URmax && URmax <= N+1);
		}
	      }
	    }
          }
        }
      }

      if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 && origIMIN==430 && origIMAX==538 && addcnt==468*/){
	printf("m=%d: LRcum= %0.10e, LS= %0.10f (startLRcum= %0.10e, startLS= %0.10f)\n", m, LRcum, LS, startLRcum, startLS);
	fflush(stdout);
      }      

      /* locate total alignment LR based on AL & UR : only check results where J >= M+1-DELTA_XU (OR EXTEND && I >= N+1-DELTA_YU && J >= M+1 -DELTA_XU - EXTEND_MAX)
	 unless PoutlierEnd > 0.0, in which case check end outliers for all cases (within the bounds of IMIN,IMAX,Jmin[],Jmax[]) */
      //      double LRcum = 0.0, LS = max(0.0,LogScaleL[IMIN] - logMaxScale);

      if(DEBUG>=2 || !PQ_BND){
	double origLRcum = LRcum, origLS = LS;

	LRcum = 0.0;
	LS = max(0.0, LogScaleL[IMIN] - logMaxScale);
	for(int I = IMIN; I <= IMAX; I++){
	  double ScaleL = exp(LogScaleL[I] - LS);
	  if(DEBUG) assert(isfinite(ScaleL) && ScaleL >= 0.0 && ScaleL < MaxScale * ((MFLOAT)1.01));

	  double LRcumI = 0.0;
	  for(int K = 0; K <= Kmax[I]; K++){
	    int jmin = (PoutlierEnd > 0.0) ? Jmin[I] : (EXTEND && I >= N+1-DELTA_YU) ? max(Jmin[I],M+1 - DELTA_XU - EXTEND_MAX) : max(Jmin[I],M+1-DELTA_XU);
	    int jmax = Jmax[I];
	    for(int J = jmin; J <= jmax; J++){
	      if(DEBUG>=2) assert(A.AL(I,K,J) >= 0.0);
	      if((MTRACE && MDEBUG_M >= 0 && N==MDEBUG_N && m==MDEBUG_M && A.AL(I,K,J) * A.UR(I,K,J) * ScaleL > MDEBUG_MIN) ||
		 (VERB>=3 && rverb && m==MDEBUG_M /* && (double)A.AL(I,K,J) * (double)A.UR(I,K,J) > LRcumI * 1e-6 */)){
		printf("m=%d:I=%d,K=%d,J=%d:AL=%0.10e,UR=%0.10e,LogScaleL[I]=%0.6e,LS=%0.6e:LRcum = %0.10e(delta=%0.10e)\n",
		       m,I,K,J,A.AL(I,K,J), A.UR(I,K,J), LogScaleL[I], LS, LRcum + (LRcumI + (double)A.AL(I,K,J) * (double)A.UR(I,K,J))*ScaleL, (double)A.AL(I,K,J) * (double)A.UR(I,K,J) * ScaleL);
		fflush(stdout);
	      }
	      LRcumI += (double)A.AL(I,K,J) * (double)A.UR(I,K,J); /* Note that AL already includes PM */
	    }
	    if(DEBUG) assert(isfinite(LRcumI));
	  }
	  LRcum += LRcumI * ScaleL;

	  if(1 /* WAS LRcum > 1.0e+100 */){  
	    double origLS = LS, origLRcum = LRcum;
	    if(LRcum > LRCUM_MAX /* WAS && LS+log(LRcum) > LogScaleL[I+1] */){/* update scaling factor LS = LS+log(LRcum) */
	      LS += log(LRcum);
	      LRcum = 1.0;
	    }
	    if(LS < LogScaleL[I+1] - logMaxScale + 0.001){
	      double Scale = exp(LS - LogScaleL[I+1] + logMaxScale - 0.001);
	      LS = LogScaleL[I+1] - logMaxScale + 0.001;
	      LRcum *= Scale;
	    }
	    if(LRCUM_FIX && LRcum > 0.0 && LRcum <= MinFloat && LS > LogScaleL[I+1] - logMaxScale + 0.001){
	      double pLS = LS;
	      LS = LogScaleL[I+1] - logMaxScale + 0.001;
	      double Scale = exp(pLS - LS);
	      if(LRcum * Scale > 1.0){
		Scale = 1.0 / LRcum;
		LS = pLS - log(Scale);
	      }
	      LRcum *= Scale;
	    }
	    if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 && origIMIN==430 && origIMAX==538 && addcnt==468 */){
	      printf("m=%d,I=%d:j=%d..%d:LRcumI=%0.8e,LRcum=%0.8e->%0.8e,ScaleL=%0.8e,LS=%0.8f->%0.8f,LogScaleL[I+1]=%0.8f(B)\n",
		     m,I,Jmin[I],Jmax[I],LRcumI,origLRcum,LRcum,ScaleL,origLS,LS,LogScaleL[I+1]);
	      fflush(stdout);
	    }

	    if(DEBUG && !(isfinite(LS) && isfinite(LRcum) && (LRcum > DBL_MIN || origLRcum <= DBL_MIN))){
              #pragma omp critical
	      {
		printf("\nm=%d,I=%d,IMIN=%d,IMAX=%d:LRcum= %0.6e -> %0.6e,LS= %0.6e -> %0.6e,LogScaleL[I]=%0.6e,LogScaleL[I+1]=%0.6e,LRcumI=%0.6e,ScaleL=%0.6e\n",
		       m,I,IMIN,IMAX,origLRcum,LRcum,origLS,LS,LogScaleL[I],LogScaleL[I+1],LRcumI,ScaleL);
		fflush(stdout);
		assert(isfinite(LS));
		assert(isfinite(LRcum));
		if(DEBUG && !(LRcum > DBL_MIN || origLRcum <= DBL_MIN)){
		  printf("WARNING: Best Alignment score for molecule m=%d underflowed!\n",m);
		  fflush(stdout);
		  //		assert((LRcum > MinFloat || origLRcum <= MinFloat));
		}
	      }
	    }
	  }
	}// for I = IMIN .. IMAX
    
        if(VERB>=3 && rverb && m==MDEBUG_M){
          #pragma omp critical
	  {
	    printf("m=%d/%d:mapid=%d,id=%lld,flip=%d:LRcum=%0.10e, LS=%0.10f (origLRcum=%0.10e, origLS=%0.10f)\n",
		 m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,pcontig->flip[m],LRcum, LS, origLRcum, origLS);
	    fflush(stdout);
	  }
	}

        if(DEBUG>=1+RELEASE && origLRcum > 0.0 && LRcum > 0.0 && !(fabs(log(LRcum)+LS - log(origLRcum) - origLS) <= fabs(log(origLRcum)+origLS) * 1e-4 + 1e-2)){
	  #pragma omp critical
	  {
	    printf("\nWARNING:m=%d/%d:mapid=%d(id=%lld):origLRcum = %0.10e, LRcum=%0.10e, origLS=%0.10f, LS=%0.10f, oldLP[m]= %0.10f,origIMIN=%d,IMIN=%d,IMAX=%d,origIMAX=%d\n",
		   m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,origLRcum,LRcum, origLS, LS, oldLP ? oldLP[m] : -999.999, origIMIN,IMIN,IMAX,origIMAX);
	    printf("\t log(LRcum)+LS = %0.10f, log(origLRcum) + origLS = %0.10f (err= %0.10f)\n", log(LRcum) + LS, log(origLRcum) + origLS, log(LRcum) + LS - log(origLRcum) - origLS);
	    fflush(stdout);
	  }
	}
	if(DEBUG) assert(LS >= 0.0);
      }

      if(DEBUG>=2){
        for(int I = IMIN; I <= IMAX; I++){
          for(int K = 0; K <= Kmax[I]; K++){
	    for(int J = Jmin[I]; J <= Jmax[I]; J++){
	      if(J <= DELTA_XU || (EXTEND && I-K <= DELTA_YU && J <= DELTA_XU + EXTEND_MAX)){
		int ULmin = A.ULmin(I,K,J);
		int ULmax = A.ULmax(I,K,J);
		if(DEBUG && !(ULmin == I-K && ULmax == I-K-1) && !(0 <= ULmin && ULmin <= ULmax && ULmax < I-K)){
		  printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:ULmin=%d,ULmax=%d\n",m,I,K,J,N,M,ULmin,ULmax);
		  fflush(stdout);
		  assert(0 <= ULmin && ULmin <= ULmax && ULmax < I-K);
		}
	      }
	      if(J >= M+1-DELTA_XU || (EXTEND && I >= N+1-DELTA_YU && J >= M+1-DELTA_XU - EXTEND_MAX)){
		int URmin = A.URmin(I,K,J);
		int URmax = A.URmax(I,K,J);
		if(DEBUG && !(URmin == I+1 && URmax == I) && !(I < URmin && URmin <= URmax && URmax <= N+1)){
		  printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:URmin=%d,URmax=%d\n",m,I,K,J,N,M,URmin,URmax);
		  fflush(stdout);
		  assert(I < URmin && URmin <= URmax && URmax <= N+1);
		}
	      }
	    }
          }
        }
      }

      if(DEBUG>=2 || PQ_BND){/* compute LRcum using AR and UL : only check results where J <= DELTA_XU (OR EXTEND && I-K <= DELTA_YU && J <= DELTA_XU + EXTEND_MAX)
				unless PoutlierEnd > 0.0, in which case check end outliers for all cases (within the bounds of IMIN,IMAX,Jmin[],Jmax[]) */
      
	double origLRcum = LRcum;
	double origLS = LS;

	LRcum = 0.0;
	LS = max(0.0,LogScaleR[IMAX] - logMaxScale);

	for(int I = IMAX; I >= IMIN; I--){
	  double ScaleR = exp(LogScaleR[I] - LS);
	  if((VERB>=3 && m==0) || (DEBUG && !(isfinite(ScaleR) && ScaleR >= 0.0 && ScaleR < MaxScale * ((MFLOAT)1.01)))){
	    printf("m=%d,I=%d:LRcum=%0.8e,ScaleR=%0.8e,LS=%0.8f,LogScaleR[I]=%0.8f,LogScaleR[I-1]=%0.8f,IMAX=%d,IMIN=%d,MaxScale1=%0.8e,logMaxScale=%0.8f(C1)\n",
		   m,I,LRcum,ScaleR,LS,LogScaleR[I],LogScaleR[I-1],IMAX,IMIN,logMaxScale,MaxScale * ((MFLOAT)1.01));
	    for(int i = IMAX; i >= I; i--)
	      printf("LogScaleR[%d] = %0.8f\n",i,LogScaleR[i]);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(isfinite(ScaleR) && ScaleR >= 0.0 && ScaleR < MaxScale * ((MFLOAT)1.01));

	  double LRcumI = 0.0;
	  for(int K = 0; K <= Kmax[I]; K++){
	    int jmin = Jmin[I];
	    int jmax = (PoutlierEnd > 0.0) ? Jmax[I] : (EXTEND && I-K <= DELTA_YU) ? min(Jmax[I],DELTA_XU + EXTEND_MAX) : min(Jmax[I],DELTA_XU);
	    for(int J = jmin; J <= jmax; J++){
	      if(MTRACE && MDEBUG_M >= 0 && N==MDEBUG_N && m==MDEBUG_M && A.AR(I,K,J) * A.UL(I,K,J) * ScaleR > MDEBUG_MIN){
		printf("m=%d:I=%d,K=%d,J=%d:AR=%0.10e,UL=%0.10e,LogScaleR[I]=%0.6e,LS=%0.6e:LRcum = %0.10e(delta=%0.10e)\n",
		       m,I,K,J, A.AR(I,K,J), A.UL(I,K,J), LogScaleR[I], LS, LRcum + (LRcumI + A.AR(I,K,J) * A.UL(I,K,J))*ScaleR, A.AR(I,K,J)*A.UL(I,K,J)*ScaleR);
		fflush(stdout);
	      }
	      LRcumI += (double)A.AR(I,K,J) * (double)A.UL(I,K,J); /* Note that AR already includes PM */
	      if(DEBUG>=2 && !(A.AR(I,K,J) >= 0.0)){
		printf("m=%d/%d:mapid=%d,I=%d(%d..%d),K=%d(%d..%d),J=%d(%d..%d):LRcumI=%0.10e,AR=%0.10e,UL=%0.10e\n",
		       m,MD,pcontig->contig[m].mapid,I,IMIN,IMAX,K,0,Kmax[I],J,Jmin[I],Jmax[I], LRcumI, A.AR(I,K,J), A.UL(I,K,J));
		fflush(stdout);
		assert(A.AR(I,K,J) >= 0.0);
	      }
	      if(DEBUG && !isfinite(LRcumI)){
		printf("\nm=%d/%d:mapid=%d,I=%d(%d..%d),K=%d(%d..%d),J=%d(%d..%d):LRcumI=%0.10e,AR=%0.10e,UL=%0.10e\n",
		       m,MD,pcontig->contig[m].mapid,I,IMIN,IMAX,K,0,Kmax[I],J,Jmin[I],Jmax[I], LRcumI, A.AR(I,K,J), A.UL(I,K,J));
		fflush(stdout);
		assert(isfinite(LRcumI));
	      }
	    }
	    if(DEBUG) assert(isfinite(LRcumI) && LRcumI >= 0.0 && LRcumI < DBL_MAX);
	  }
	  LRcum += LRcumI * ScaleR;
	  if(DEBUG) assert(isfinite(LRcum) && LRcum >= 0.0 && LRcum < DBL_MAX);

	  if(1/* WAS LRcum > 1e+100*/){  
	    double origLRcum = LRcum;
	    double origLS = LS;
	    double origScaleR = ScaleR;
	    if(LRcum > LRCUM_MAX){/* update scaling factor LS = LS+log(LRcum) */
	      LS += log(LRcum);
	      LRcum = 1.0;
	    }
	    if(LS < LogScaleR[I-1] - logMaxScale + 0.001){
	      double Scale = exp(LS - LogScaleR[I-1] + logMaxScale - 0.001);
	      LS = LogScaleR[I-1] - logMaxScale + 0.001;
	      LRcum *= Scale;
	    }
	    if(LRCUM_FIX && LRcum > 0.0 && LRcum <= MinFloat && LS > LogScaleR[I-1] - logMaxScale + 0.001){
	      double pLS = LS;
	      LS = LogScaleR[I-1] - logMaxScale + 0.001;
	      double Scale = exp(pLS - LS);
	      if(LRcum * Scale > 1.0){
		Scale = 1.0 / LRcum;
		LS = pLS - log(Scale);
	      }
	      LRcum *= Scale;
	    }
	    if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 && origIMIN==430 && origIMAX==538 && addcnt==468 */){
	      printf("m=%d,I=%d:j=%d..%d:LRcumI=%0.8e,LRcum=%0.8e->%0.8e,ScaleR=%0.8e,LS=%0.8f->%0.8f,LogScaleR[I-1]=%0.8f,logMaxScale=%0.8f(C)\n",
		     m,I,Jmin[I],Jmax[I],LRcumI,origLRcum,LRcum,origScaleR,origLS,LS,LogScaleR[I-1],logMaxScale);
	      fflush(stdout);
	    }
	    if(DEBUG){
	      assert(isfinite(LS));
	      assert(isfinite(LRcum) && LRcum <= LRCUM_MAX);
	      assert(LRcum > 0.0 || origLRcum <= 0.0);
	    }
	  }
	  if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 && origIMIN==430 && origIMAX==538 && addcnt==468*/){
	    printf("m=%d,I=%d:j=%d..%d:LRcum=%0.8e,ScaleR=%0.8e,LS -> %0.8f,LogScaleR[I-1]=%0.8f(C2)\n",
		   m,I,Jmin[I],Jmax[I],LRcum,ScaleR,LS,LogScaleR[I-1]);
	    fflush(stdout);
	  }
	}

	if(VERB>=3 && rverb && m==MDEBUG_M /* && MD==792 && origIMIN==430 && origIMAX==538 && addcnt==468*/){
	  printf("m=%d: origLRcum= %0.10e, origLS= %0.10f, LRcum= %0.10e, LS= %0.10f\n",m,origLRcum,origLS,LRcum,LS);
	  fflush(stdout);
	}

	/* NOTE : If the following condition triggers, LRcum differs between summation of AL*UR and AR*UL, so probably float precision is not sufficient. To recover use version with MFLOAT=0 */

        if(DEBUG>=1+RELEASE && EPS2 > 0.0 && origLRcum > 0.0 && LRcum > 0.0 && !(fabs(log(LRcum)+LS - log(origLRcum) - origLS) <= fabs(log(origLRcum)+origLS) * 1e-4 + 1e-2) && TBwtm > 0.01){
	  #pragma omp critical
	  {
	    printf("\nWARNING:m=%d/%d:mapid=%d(id=%lld),TBwt= %0.6f:origLRcum = %0.10e, LRcum=%0.10e, origLS=%0.10f, LS=%0.10f, oldLP[m]= %0.10f,origIMIN=%d,IMIN=%d,IMAX=%d,origIMAX=%d\n",
	  m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,(double)TBwtm,origLRcum,LRcum, origLS, LS, oldLP ? oldLP[m] : -999.999, origIMIN,IMIN,IMAX,origIMAX);
#if DEBUG>=2 // Extra debug code	    
	    for(int I = origIMIN; I <= origIMAX; I++){
	      printf("I=%d,K=0..%d,J=%d,%d..%d,%d:LogScaleL[I]= %0.10f, LogScaleR[I]= %0.10f\n",I,Kmax[I],origJmin[I],Jmin[I],Jmax[I],origJmax[I],LogScaleL[I],LogScaleR[I]);
	      for(int K = 0; K <= Kmax[I]; K++){
		for(int J = origJmin[I];J <= origJmax[I]; J++)
		  printf("m=%d,I=%d,K=%d,J=%d: AR(I,K,J)= %0.10e, UL(I,K,J)= %0.10e, AL(I,K,J)= %0.10e, UR(I,K,J)= %0.10e, PM(I,K)= %0.10e\n",
			 m,I,K,J,A.AR(I,K,J),A.UL(I,K,J), A.AL(I,K,J), A.UR(I,K,J), PM(0,I,K,Y));
	      }
	    }
	    double nLRcum = LRcum, nLS = LS;
	    LS = max(0.0,LogScaleR[origIMAX] - logMaxScale);
	    LRcum = 0.0;

	    for(int I = origIMAX; I >= origIMIN; I--){
	      double ScaleR = exp(LogScaleR[I] - LS);
	      double LRcumI = 0.0;
	      for(int K = 0; K <= Kmax[I]; K++){
		int jmin = origJmin[I];
		int jmax = (PoutlierEnd > 0.0) ? origJmax[I] : (EXTEND && I-K <= DELTA_YU) ? min(origJmax[I],DELTA_XU + EXTEND_MAX) : min(origJmax[I],DELTA_XU);
		for(int J = jmin; J <= jmax; J++)
		  LRcumI += (double)A.AR(I,K,J) * (double)A.UL(I,K,J); /* Note that AR already includes PM */
	      }
	      LRcum += LRcumI * ScaleR;

	      if(1){  
		double origLS = LS;
		double origLRcum = LRcum;
		if(LRcum > LRCUM_MAX){/* update scaling factor LS = LS+log(LRcum) */
		  LS += log(LRcum);
		  LRcum = 1.0;
		}
		if(LS < LogScaleR[I-1] - logMaxScale + 0.001){
		  double Scale = exp(LS - LogScaleR[I-1] + logMaxScale - 0.001);
		  LS = LogScaleR[I-1] - logMaxScale + 0.001;
		  LRcum *= Scale;
		}
		if(LRCUM_FIX && LRcum > 0.0 && LRcum <= MinFloat && LS > LogScaleR[I-1] - logMaxScale + 0.001){
		  double pLS = LS;
		  LS = LogScaleR[I-1] - logMaxScale + 0.001;
		  double Scale = exp(pLS - LS);
		  if(LRcum * Scale > 1.0){
		    Scale = 1.0 / LRcum;
		    LS = pLS - log(Scale);
		  }
		  LRcum *= Scale;
		}
		if(VERB>=3 && rverb && m==553 && MD==585 /* && origIMIN==430 && origIMAX==538 */ && addcnt==37){
		  printf("m=%d,I=%d:j=%d..%d:LRcumI=%0.8e,LRcum=%0.8e->%0.8e,ScaleR=%0.8e,LS=%0.8f->%0.8f,LogScaleR[I-1]=%0.8f,logMaxScale=%0.8f(D)\n",
			 m,I,jmin,jmax,LRcumI,origLRcum,LRcum,ScaleR,origLS,LS,LogScaleR[I-1],logMaxScale);
		  fflush(stdout);
		}
	      }
	      if(VERB>=3 && rverb && m==553 && MD==585 /* && origIMIN==430 && origIMAX==538 */ && addcnt==37){
		printf("m=%d,I=%d:j=%d..%d:LRcum=%0.8e,ScaleR=%0.8e,LS=%0.8f,LogScaleR[I-1]=%0.8f(D2)\n",
		       m,I,jmin,jmax,LRcum,ScaleR,LS,LogScaleR[I-1]);
		fflush(stdout);
	      }
	    }

	    LRcum = nLRcum;
	    LS = nLS;
	    printf("\nSummary:m=%d/%d:mapid=%d(id=%lld):origLRcum = %0.10e, LRcum=%0.10e, origLS=%0.10f, LS=%0.10f, oldLP[m]= %0.10f,origIMIN=%d,IMIN=%d,IMAX=%d,origIMAX=%d\n",
		   m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,origLRcum,LRcum, origLS, LS, oldLP ? oldLP[m] : -999.999, origIMIN,IMIN,IMAX,origIMAX);
#endif // Extra debug code
	    printf("\t log(LRcum)+LS = %0.10f, log(origLRcum) + origLS = %0.10f (err= %0.10f)\n", log(LRcum) + LS, log(origLRcum) + origLS, log(LRcum) + LS - log(origLRcum) - origLS);
	    fflush(stdout);

	    if(!RELEASE && (USE_MFLOAT==0 || DEBUG>=2))
	      assert(fabs(log(LRcum) + LS - log(origLRcum) - origLS) <= fabs(log(origLRcum) + origLS)*1e-4 + 1e-4);
	  }
	}
	if(!PQ_BND){
          LRcum = origLRcum;
	  LS = origLS;
	}
      }

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After 3-D recurrance to compute AL,AR:wall time= %0.6f\n", m, tid, wtime() - wt);
	  fflush(stdout);
        }
      }

      FL128b LRratio = 1.0;/* the drop in exp(LS + log(LRcum + LRbiasLS)) due to EPS2 & EPS3 */
      if(EPS_ADJUST && EPS2 > 0.0){
	if(startLS != LS){/* adjust startLS to match LS */
	  startLRcum *= exp(startLS - LS);
	  if(DEBUG) assert(isfinite(startLRcum));
	  if(DEBUG && LRcum > 0.0) assert(startLRcum > 0.0);
	}
	/* NOTE1 : the predicted LR values for current molecule m will be scaled up by startLRcum/LRcum to compenstate for EPS2 & EPS3.
	   This correction will be applied when computing newLP[m], newLPd[m][S], newLPa[m][t] and newLPdelta[m][d] */
	/* NOTE2 : if PQ_BND && USE_MFLOAT==0, it is possible that LRcum (derived from AR & UL) differs from origLRcum (derived from AL * UR after EPS & EPS3) : this adjustment
	   will also correct for that difference so the result matches the AL*UR based value computed by qprobeval */

	if(LRcum > 0.0){
	  LRratio = startLRcum/LRcum;
	  if(DEBUG) assert(isfinite(LRratio));
	}
      } else
	startLRcum = LRcum;
      startLS = LS;

      /* Compute LPdel[0 ... N] contribution of current map m using difference array LRdel[S= Dmin .. Dmax] (along with LRdel[0] if Dmin > 0) */
      /*       LPdel[0] += log(LRdel[0]+LRbias),  LPdel[S=Dmin..Dmax] += log(LRdel[0] + sum(LRdel[Dmin..S]) + LRbias) - log(LRdel[0] + sum(LRdel[Dmin..S-1]) + LRbias) */
      /* NOTE : the differential method of accumulation LPdel[] is faster when range Dmin..Dmax is much smaller than 0..N */
      /* NOTE : LRdel[0] == LRcum is only computed if (FIXLRDEL || LDEBUG>=1), LRdel[N+1] == LRcum will only be computed if (IMAX >= N or LDEBUG>=2) */
      int Dmin = (LDEBUG>=2 ? 0 : IMIN), Dmax = (LDEBUG>=2 ? N : IMAX) + 1;
      if(LPdelG){
	if(KMAX_APPROX){
	  if(IMIN - Kmax[IMIN] < Dmin)
	    Dmin = IMIN - Kmax[IMIN];
	} else {
	  for(int I= IMIN; I <= min(IMAX,IMIN+KMAX-1); I++)
	    if(I - Kmax[I] < Dmin)
	      Dmin = I - Kmax[I];/* need to accomodate smallest value of I-K */
	}
	for(int i = Dmin; i <= Dmax; i++)
	  LRdel[i] = 0.0;
	if(LDEBUG||FIXLRDEL) LRdel[0] = 0.0;
      }

      /* Compute LPaddI & LPadd[0..addcnt-1] contribution of current map m using difference array LRadd[S= 0,Amin .. Amax] and LRaddT[S= Amin..Amax][t= 0 .. addcntI[S]-1] for added site location addlocI[S][t] */
      /* Compute LPdelta[S=0..N][d=0..D[S]-1]] contribution of current map m using difference array LRdelta[S=0,Amin .. Amax] and LRdeltaD[S=Amin..min(N,Amax)][d=0..D[S]-1] for change delta[I][d] in interval Y[S] .. Y[S+1] */
      int Amin = (LDEBUG>=2 ? 0 : IMIN-1), Amax = (LDEBUG>=2 ? N : IMAX) + 1, Doffset = 0;
      if(SCORE_APPROX >= 2){
	if(IMIN - Kmax[IMIN] < Amin)
	  Amin = IMIN - Kmax[IMIN];
      } else {
	for(int I=IMIN; I <= min(IMAX,IMIN+KMAX-1); I++)
	  if(I - Kmax[I] < Amin)
	    Amin = I - Kmax[I];/* need to accomodate smallest value of I-K */
      }

      if(VERB>=2 && m==5430 && newLPdelta){
	printf("m=%d/%d: Initial Amin=%d,Amax=%d\n",m,MD,Amin,Amax);
	fflush(stdout);
      }

      LRadd[0] = LRadd2[0] = LRdelta[0] = 0.0;
      for(int i = Amin; i <= Amax; i++){
	LRadd[i] = LRadd2[i] = LRdelta[i] = 0.0;
	if(i <= N){
	  if(myaddcnt > 0)
	    for(int t = 0; t < addcntI[i]; t++)
	      LRaddT[i][t] = 0.0;
	  if(newLPdelta)
	    for(int d = 0; d < ((max(0,D[i])/* NEW10*/ + VPAD) & (~VPAD)); d++)
	      LRdeltaD[i][d] = 0.0;
	}
      }
      if(LDEBUG)
	for(int i = 0; i <= N+1; i++)
	  LRsite[i] = LRfrag[i] = LRfrag2[i] = 0.0;

      if(VERB>=2){
	printf("mprobeval:m=%d,N=%d,MD=%d,MDEBUG_M=%d,MDEBUG_N=%d,MDEBUG_S=%d,MDEBUG_F=%d\n",m,N,MD,MDEBUG_M,MDEBUG_N,MDEBUG_S,MDEBUG_F);
	fflush(stdout);
      }

      if(VERB>=2 && newLPdelta && rverb && LPdelta){
	printf("mprobeval:m=%d,Amin=%d,Amax=%d:Before any updates:time=%0.6f\n",m,Amin,Amax,wtime());
	fflush(stdout);
      }

      if(DEBUG>=2 && newLPdelta){
	if(VERB>=3){
	  if(m==7 && Amin <= 12 && 12 <= Amax && D[12] >= 1 /* && !isfinite(LRdeltaD[12][0])*/)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[12]=%d,delta[12][0]= %0.4f,LRdeltaD[12][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[12],delta[12][0],LRdeltaD[12][0]);
	  if(m==56 && Amin <= 31 && 31  <= Amax && D[31] >= 1 /* && !isfinite(LRdeltaD[31][0] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[31]=%d,delta[31][0]= %0.4f,LRdeltaD[31][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[31],delta[31][0],LRdeltaD[31][0]);
	  if(m==20 && Amin <= 33 && 33 <= Amax && D[33] >= 11 /* && !isfinite(LRdeltaD[33][10] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[33]=%d,delta[33][10]= %0.4f,LRdeltaD[33][10]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[33],delta[33][10],LRdeltaD[33][10]);
	  fflush(stdout);
	}
	for(int S = Amin;S <= min(N,Amax); S++)
	  for(int d = 0; d < D[S]; d++)
	    assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
      }

      if(MDEBUG_M >= 0){
	if(MDEBUG_S >= 0 && m==MDEBUG_M && N==MDEBUG_N && LPdelG){
	  FL128 sum = (LDEBUG||FIXLRDEL) ? LRdel[0] : LRcum;
	  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
	    sum += LRdel[i];
	  printf("m=%d:Before any updates:LRdel[0..%d]=%0.10e(LS=%0.10f),Dmin=%d,Dmax=%d,N=%d,IMIN=%d,IMAX=%d\n",m,MDEBUG_S,(double)sum,LS,Dmin,Dmax,N,IMIN,IMAX);
	  if(DEBUG && sum != 0.0){
	    for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
	      printf("LRdel[%d]=%0.10e\n",i,(double)LRdel[i]);
	    assert(sum != 0.0);
	  }
	}
	if(LDEBUG && MDEBUG_F >= 0 && m==MDEBUG_M && N==MDEBUG_N-1){
	  FL128 sum = 0.0, sum2 = 0.0;
	  for(int i = 0; i <= MDEBUG_F; i++){
	    sum += LRfrag[i];
	    sum2 += LRfrag2[i];
	  }
	  printf("m=%d:Before any updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f),IMIN=%d,IMAX=%d,N=%d\n",m,MDEBUG_F,(double)sum,(double)sum2,(double)(sum-sum2),LS,IMIN,IMAX,N);
	}
	if((MDEBUG_T >= 0 || newLPdelta) && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  FL128 sum = LRadd[0], sum2 = LRadd2[0];
	  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
	    sum += LRadd[i];
	    sum2 += LRadd2[i];
	  }
	  FL128 sum3 = sum;
	  if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
	    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:Before any updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
	  } else
	    printf("m=%d:Before any updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e\n", m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2));
	  if(LDEBUG){
	    sum = 0.0,sum2=0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:Before any updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2));
	  }
	  for(int I = IMIN; I <= IMAX; I++)
	    printf("\t I=%d:Jmin[I]=%d,Jmax[I]=%d,Kmax[I]=%d\n",I,Jmin[I],Jmax[I],Kmax[I]);
	}
	if(MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N+1){
	  FL128 sum = 0.0;
	  for(int i = 1; i <= MDEBUG_TS+1; i++)
	    sum += LRsite[i]; 
	  printf("m=%d:Before any updates:LRsite[0..%d]=%0.10e\n",m,MDEBUG_TS+1,(double)sum);
	  for(int I = IMIN; I <= IMAX; I++)
	    printf("\t I=%d:Jmin[I]=%d,Jmax[I]=%d,Kmax[I]=%d\n",I,Jmin[I],Jmax[I],Kmax[I]);
	}
	if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  if(newLPdelta && MDEBUG_TT < D[min(N,MDEBUG_TS)]){
	    FL128 sum = LRdelta[0];
	    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
	      sum += LRdelta[i];
	    FL128 sum3 = sum;
	    if(Amin <= MDEBUG_TS && MDEBUG_TS <= min(N,Amax) && MDEBUG_TT < D[MDEBUG_TS])
	      sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:Before any updates:S=%d,D=%d,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3);
	  }
	  if(LDEBUG){
	    FL128 sum = 0.0,sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:Before any updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2));
	  }
	  for(int I = IMIN; I <= IMAX; I++)
	    printf("\t I=%d:Jmin[I]=%d,Jmax[I]=%d,Kmax[I]=%d\n",I,Jmin[I],Jmax[I],Kmax[I]);
	}
	fflush(stdout);
      }
      if(LDEBUG){/* check that cumulative value of LRsite[1..N+1] approximates 0 */
	FL128 LRsitecum = 0.0;
	for(int S = 1; S <= N+1; S++){
	  LRsitecum += LRsite[S];
	  if(VERB && (VERB>=2 || (DEBUG && !(S <= N || fabs((double)LRsitecum) <= LRcum * 1e-8)))){
	    printf("m=%d:S=%d,N=%d:LRsite[S]=%0.10e,LRsitecum=%0.10e,LRcum=%0.10e,LRfrag[0]=%0.10e,LRsitecum-LRcum=%0.10e,IMIN=%d,IMAX=%d\n",
		   m,S,N,(double)LRsite[S],(double)LRsitecum,LRcum,(double)LRfrag[0],(double)(LRsitecum-LRcum),IMIN,IMAX);
	    for(int s = 1; s <= N+1; s++)
	      printf("LRsite[%d]=%0.10e,LRfrag[%d]=%0.10e\n",s,(double)LRsite[s],s,(double)LRfrag[s]);
	    fflush(stdout);
	  }
	  if(DEBUG && S == N+1) assert(fabs((double)LRsitecum) <= LRcum * 1e-8);
	}
      }

      if(DEBUG>=2 && newLPdelta){
	if(VERB>=3){
	  if(m==7 && Amin <= 12 && 12 <= Amax && D[12] >= 1 /* && !isfinite(LRdeltaD[12][0])*/)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[12]=%d,delta[12][0]= %0.4f,LRdeltaD[12][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[12],delta[12][0],LRdeltaD[12][0]);
	  if(m==56 && Amin <= 31 && 31  <= Amax && D[31] >= 1 /* && !isfinite(LRdeltaD[31][0] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[31]=%d,delta[31][0]= %0.4f,LRdeltaD[31][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[31],delta[31][0],LRdeltaD[31][0]);
	  if(m==20 && Amin <= 33 && 33 <= Amax && D[33] >= 11 /* && !isfinite(LRdeltaD[33][10] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[33]=%d,delta[33][10]= %0.4f,LRdeltaD[33][10]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[33],delta[33][10],LRdeltaD[33][10]);
	  fflush(stdout);
	}
	for(int S = Amin;S <= min(N,Amax); S++)
	  for(int d = 0; d < D[S]; d++)
	    assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
      }

      /* Updates based on PoutlierEnd (Ploc term only on Line 1 & 2 of Equation 7 for LRdel and Equation 10 for LRadd) */
      /* Update of LRdelta[] (Line 1 & 2 of Equation 21) is same as for LRadd and will just be copied over */
      for(int I= IMIN; I <= IMAX; I++){
	FL128b ScaleL = PoutlierEnd * exp(LogScaleL[I]-LS);
	FL128b ScaleR = PoutlierEnd * exp(LogScaleR[I]-LS);
	for(int K = 0; K <= Kmax[I]; K++){
          for(int J = Jmin[I]; J <= Jmax[I]; J++){
	    if(DEBUG>=2) assert(A.AL(I,K,J) >= 0.0);
	    FL128b AL = A.AL(I,K,J) * ScaleL;
	    if(DEBUG>=2) assert(isfinite(AL) && AL >= 0.0);

	    if(LPdelG){
	      /* LRdel[] S range for AL * Ploc is I+1 .. N  */
	      if(DEBUG>=2 && !(Dmin <= I+1  && I+1 <= Dmax)){
		printf("m=%d:I=%d,K=%d,J=%d:Dmin=%d,Dmax=%d\n",m,I,K,J,Dmin,Dmax);
		fflush(stdout);
		assert(Dmin <= I+1  && I+1 <= Dmax);
	      }
	      FL128 origLRdelI1 = LRdel[I+1];
	      if(DEBUG>=2 && !isfinite(LRdel[I+1])){	    
		double origLR1 = origLRdelI1;
		FL128 origLR2 = origLRdelI1 - origLR1;
		printf("m=%d:I=%d,K=%d,J=%d:orig LRdel[I+1]=%0.17le+%0.17le\n",
		       m,I,K,J,origLR1,(double)origLR2);
		fflush(stdout);
		assert(isfinite(LRdel[I+1]));
	      }
	      LRdel[I+1] += AL;
	      if(DEBUG>=2 && !isfinite(LRdel[I+1])){	    
		double origLR1 = origLRdelI1;
		FL128 origLR2 = origLRdelI1 - origLR1;
		double LR1 = LRdel[I+1];
		FL128 LR2 = LRdel[I+1]-LR1;
		printf("m=%d:I=%d,K=%d,J=%d:p->AL=%0.17le,ScaleL=%0.24Le,AL=%0.24Le:\n\t LRdel[I+1]=%0.17le+%0.17le->%0.17le+%0.17le\n",
		       m,I,K,J, A.AL(I,K,J),(long double)ScaleL,(long double)AL,origLR1,(double)origLR2,LR1,(double)LR2);
		fflush(stdout);
		assert(isfinite(LRdel[I+1]));	    
	      }
	      if(LDEBUG>=2) LRdel[N+1] -= AL;
	      if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && I+1 <= MDEBUG_S && MDEBUG_S <= N && AL >= MDEBUG_MIN){
		FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		  sum += LRdel[i];
		printf("m=%d:I=%d,K=%d,J=%d:AL*Ploc=%0.10e,LRdel[0..%d]=%0.10e\n",m,I,K,J,(double)AL,MDEBUG_S,(double)sum);
		fflush(stdout);
	      }
	    }

	    if(LDEBUG){/* LRsite[S]: S range for AL*Ploc is I+1 .. N */
	      LRsite[I+1] += AL;
	      LRsite[N+1] -= AL;
	      if(LDEBUG && DEBUG>=2){/* verify that LRsite[1..N+1] adds up to approximately 0.0 */
		FL128 sum = 0.0;
		for(int S = 1; S <= N+1; S++)
		  sum += LRsite[S];
		if((DEBUG && !(fabs((double)sum) <= LRcum * 1e-8)) ||(0 &&m==MDEBUG_M && AL >= MDEBUG_MIN && (I+1==STRACE || N+1==STRACE))){
		  printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:AL=%0.10e,LRsite[1..N+1]=%0.10e,LRsite[%d]+=%0.10e,LRsite[%d]-=%0.10e\n",
			 m,I,K,J,N,M,(double)AL,(double)sum,I+1,(double)LRsite[I+1],N+1,(double)LRsite[N+1]);
		  fflush(stdout);
		  assert(fabs((double)sum) <= LRcum * 1e-8);
		}
		if((MDEBUG_M < 0 || m==MDEBUG_M) && ((MDEBUG_ST >= 0 && I+1 <= MDEBUG_ST && MDEBUG_ST <= N) ||
						     (MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && I+1 <= MDEBUG_TS+1 && MDEBUG_TS+1 <= N && N == MDEBUG_N+1)) && AL >= MDEBUG_MIN){
		  FL128 sum = 0.0;
		  for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++)
		    sum += LRsite[i];
		  if(sizeof(FL128) > sizeof(double)){
		    double AL_double = AL;
		    double sum_double = sum;
		    int t = max(MDEBUG_TS+1,MDEBUG_ST);
		    printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:AL*Ploc=%0.17e,%0.17e,LRsite[1..%d]=%0.17e,%0.17e (LRsite[%d]=%0.10e,LRsite[%d]=%0.10e)\n",
			   m,I,K,J,N,M,AL_double,(double)(AL-AL_double),t,sum_double,(double)(sum-sum_double),t,(double)LRsite[t],t+1,(double)LRsite[t+1]);
		    /* sum = 0.0;
		       for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++){
		       sum += LRsite[i];
		       sum_double = sum;
		       double LR_double = LRsite[i];
		       printf("   LRsite[%d]=%0.17e,%0.17e,sum=%0.17e,%0.17e\n",i,LR_double,(double)(LRsite[i]-LR_double),sum_double,(double)(sum-sum_double));
		       }*/
		  } else
		    printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:AL*Ploc=%0.10e,LRsite[1..%d]=%0.10e\n",
			   m,I,K,J,N,M,(double)AL,max(MDEBUG_TS+1,MDEBUG_ST),(double)sum);
		}
	      }
	    }

	    /* LRadd[S]: S range for AL * Ploc is I .. N */
	    if(DEBUG>=2) assert(Amin <= I && I <= Amax && I <= N);
	    LRadd[I] += AL;
	    if(LDEBUG>=2) LRadd2[N+1] += AL;
	    if(MDEBUG_M >= 0 && (MDEBUG_T >= 0 || newLPdelta) && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && I <= MDEBUG_TS && MDEBUG_TS <= N && AL >= MDEBUG_MIN){
	      FL128 sum = LRadd[0],sum2 = LRadd2[0];
	      for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
		sum += LRadd[i];
		sum2 += LRadd2[i];
	      }
	      FL128 sum3 = sum;
	      if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
		sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
		printf("m=%d:I=%d,K=%d,J=%d:AL*Ploc=%0.10e,S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
		       m,I,K,J,(double)AL,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
	      } else { 
		printf("m=%d:I=%d,K=%d,J=%d:AL*Ploc=%0.10e,S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
		       m,I,K,J,(double)AL,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2));
		if(0 && I==41 && K==0 && J==2){
		  printf("LRadd[0]=%0.10e-%0.10e = %0.10e\n",(double)LRadd[0],(double)LRadd2[0],(double)(LRadd[0]-LRadd2[0]));
		  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
		    printf("LRadd[%d]=%0.10e-%0.10e = %0.10e\n",i,(double)LRadd[i],(double)LRadd2[i],(double)(LRadd[i]-LRadd2[i]));
		}
	      }
	      fflush(stdout);
	    }
	    if(LDEBUG){
	      LRfrag[I] += AL;
	      LRfrag2[N+1] += AL;
	      if((MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N-1 && I <= MDEBUG_F && MDEBUG_F <= N && AL >= MDEBUG_MIN) ||
		 (MDEBUG_M >= 0 /* && MDEBUG_T >= 0 */ && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && I <= MDEBUG_TS && MDEBUG_TS <= N && AL >= MDEBUG_MIN)){
		FL128 sum = 0.0,sum2=0.0;
		for(int i = 0; i <= max(MDEBUG_F,MDEBUG_TS); i++){
		  sum += LRfrag[i];
		  sum2 += LRfrag2[i];
		}
		printf("m=%d:I=%d,K=%d,J=%d:AL*Ploc=%0.10e,LRfrag[0..%d]=%0.10e-%0.10e=%0.10e\n",m,I,K,J,(double)AL,max(MDEBUG_F,MDEBUG_TS),(double)sum,(double)sum2,(double)(sum-sum2));
		fflush(stdout);
	      }
	    }

	    if(DEBUG>=2 && !(isfinite(A.AR(I,K,J)) && A.AR(I,K,J) >= 0.0)){
	      printf("m=%d:I=%d,K=%d,J=%d:ARik[J]=%0.17e\n",m,I,K,J,A.AR(I,K,J));
	      fflush(stdout);
	      assert(isfinite(A.AR(I,K,J)) && A.AR(I,K,J) >= 0.0);
	    }
	    FL128b AR = A.AR(I,K,J) * ScaleR;
	    if(DEBUG>=2 && !(isfinite(AR))){
	      printf("m=%d:I=%d,K=%d,J=%d:p->AR=%0.17e,Scale=%0.24Le,AR=%0.24Le\n",m,I,K,J,A.AR(I,K,J),(long double)ScaleR,(long double)AR);
	      fflush(stdout);
	      
	      assert(isfinite(AR));
	    }
	    if(LPdelG){
	      /* LRdel[] S range for AR * Ploc is 1 .. I-K-1 */
	      if(LDEBUG||FIXLRDEL){
		LRdel[0] += AR;
		if(DEBUG>=2) assert(isfinite(LRdel[0]));
	      }
	      if(DEBUG>=2)  assert(Dmin <= I-K && I-K <= Dmax);
	      LRdel[I-K] -= AR;
	      if(DEBUG>=2) assert(isfinite(LRdel[I-K]));
	      if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && 1 <= MDEBUG_S && MDEBUG_S <= I-K-1 && AR >= MDEBUG_MIN){
		FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		  sum += LRdel[i];
		printf("m=%d:I=%d,K=%d,J=%d:AR*Ploc=%0.10e,LRdel[0..%d]=%0.10e\n",m,I,K,J,(double)AR,MDEBUG_S,(double)sum);
		fflush(stdout);
	      }
	    }
	    if(LDEBUG){/* LRsite[S]: S range for AR * Ploc is 1 .. I-K-1 */
	      LRsite[1] += AR;
	      LRsite[I-K] -= AR;
	      if(LDEBUG && DEBUG>=2){/* verify that LRsite[1..N+1] adds up to approximately 0.0 */
		FL128 sum = 0.0;
		for(int S = 1; S <= N+1; S++)
		  sum += LRsite[S];
		if((DEBUG && !(fabs((double)sum) <= LRcum * 1e-8)) ||(0 &&m==MDEBUG_M && AR >= MDEBUG_MIN && (1==STRACE || I-K==STRACE))){
		  printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:AR=%0.10e,LRsite[1..N+1]=%0.10e,LRsite[%d]+=%0.10e,LRsite[%d]-=%0.10e\n",
			 m,I,K,J,N,M,(double)AR,(double)sum,1,(double)LRsite[1],I-K,(double)LRsite[I-K]);
		  fflush(stdout);
		}
	      }

	      if((MDEBUG_M < 0 || m==MDEBUG_M) && ((MDEBUG_ST >= 0 && 1 <= MDEBUG_ST && MDEBUG_ST <= I-K-1) ||
						   (MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && 1 <= MDEBUG_TS+1 && N == MDEBUG_N+1)) && AR >= MDEBUG_MIN){
		FL128 sum = 0.0;
		for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++)
		  sum += LRsite[i];
		if(sizeof(FL128) > sizeof(double)){
		  double AR_double = AR;
		  double sum_double = sum;
		  int t = max(MDEBUG_TS+1,MDEBUG_ST);
		  printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:AR*Ploc=%0.17e,%0.17e,LRsite[1..%d]=%0.17e,%0.17e(LRsite[%d]=%0.10e,LRsite[%d]=%0.10e)\n",
			 m,I,K,J,N,M,AR_double,(double)(AR-AR_double),t,sum_double,(double)(sum-sum_double),t,(double)LRsite[t],t+1,(double)LRsite[t+1]);
		  /*		  sum = 0.0;
				  for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++){
				  sum += LRsite[i];
				  sum_double = sum;
				  double LR_double = LRsite[i];
				  printf("   LRsite[%d]=%0.17e,%0.17e,sum=%0.17e,%0.17e\n",i,LR_double,(double)(LRsite[i]-LR_double),sum_double,(double)(sum-sum_double));
				  }*/
		} else {
		  printf("m=%d,I=%d,K=%d,J=%d,N=%d,M=%d:AR*Ploc=%0.17e,LRsite[1..%d]=%0.17e\n",m,I,K,J,N,M,(double)AR,max(MDEBUG_TS+1,MDEBUG_ST),(double)sum);
		}
	      }
	    }

	    /* LRadd[S]: S range for AR * Ploc = 0 .. I-K-1 */
	    if(DEBUG>=2) assert(Amin <= I-K && I-K <= Amax && I-K >= 1);
	    LRadd[0] += AR;
	    LRadd2[I-K] += AR;
	    if(MDEBUG_M >= 0 && (MDEBUG_T >= 0 || newLPdelta) && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && 0 <= MDEBUG_TS && MDEBUG_TS <= I-K-1 && AR >= MDEBUG_MIN){
	      FL128 sum = LRadd[0], sum2 = LRadd2[0];
	      for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
		sum += LRadd[i];
		sum2 += LRadd2[i];
	      }
	      FL128 sum3 = sum;
	      if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
		sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
		printf("m=%d:I=%d,K=%d,J=%d:AR*Ploc=%0.10e,S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
		       m,I,K,J,(double)AR,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
	      } else
		printf("m=%d:I=%d,K=%d,J=%d:AR*Ploc=%0.10e,S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
		       m,I,K,J,(double)AR,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2));
	      fflush(stdout);
	    }
	    if(LDEBUG){/* LRfrag[S]: S range for AR * Ploc = 0 .. I-K-1 */
	      LRfrag[0] += AR;
	      LRfrag2[I-K] += AR;
	      if((MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N-1 && 0 <= MDEBUG_F && MDEBUG_F <= I-K-1 && AR >= MDEBUG_MIN) ||
		 (MDEBUG_M >= 0 /* && MDEBUG_T >= 0 */ && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && 0 <= MDEBUG_TS && MDEBUG_TS <= I-K-1 && AR >= MDEBUG_MIN)){
		FL128 sum = 0.0,sum2 = 0.0;
		for(int i = 0; i <= max(MDEBUG_F,MDEBUG_TS); i++){
		  sum += LRfrag[i];
		  sum2 += LRfrag2[i];
		}
		printf("m=%d:I=%d,K=%d,J=%d:AR*Ploc=%0.10e,LRfrag[0..%d]=%0.10e-%0.10e=%0.10e\n",m,I,K,J,(double)AR,max(MDEBUG_F,MDEBUG_TS),(double)sum,(double)sum2,(double)(sum-sum2));
		fflush(stdout);
	      }
	    }
	  }
	}
      }
      if(newLPdelta){/* copy LRadd[] - LRadd2[] to LRdelta[] */
	LRdelta[0] = LRadd[0] - LRadd2[0];
	for(int i = Amin; i <= Amax; i++)
	  LRdelta[i] = LRadd[i] - LRadd2[i];
      }

      if(VERB>=2 && newLPdelta && rverb && LPdelta){
	printf("mprobeval:m=%d,Amin=%d,Amax=%d:After PoutlierEnd updates:time=%0.6f\n",m,Amin,Amax,wtime());
	fflush(stdout);
      }

      if(MDEBUG_M >= 0){
	if(MDEBUG_S >= 0 && m==MDEBUG_M && N==MDEBUG_N && LPdelG){
	  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
	    sum += LRdel[i];
	  printf("m=%d:After PoutlierEnd updates:LRdel[0..%d]=%0.10e(LS=%0.10f, LL=%0.10f)\n",m,MDEBUG_S,(double)sum,LS, LS + log((double)sum));
	}
	if(LDEBUG && MDEBUG_F >= 0 && m==MDEBUG_M && N==MDEBUG_N-1){
	  FL128 sum = 0.0,sum2= 0.0;
	  for(int i = 0; i <= MDEBUG_F; i++){
	    sum += LRfrag[i];
	    sum2 += LRfrag2[i];
	  }
	  printf("m=%d:After PoutlierEnd updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f,LL=%0.10f)\n",m,MDEBUG_F,(double)sum,(double)sum2,(double)(sum-sum2),LS,LS+log((double)(sum-sum2)));
	}
	if(MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  FL128 sum = LRadd[0],sum2=LRadd2[0];
	  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
	    sum += LRadd[i];
	    sum2 += LRadd2[i];
	  }
	  FL128 sum3 = sum;
	  if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
	    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After PoutlierEnd updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e,LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2),LS);
	  } else 
	    printf("m=%d:After PoutlierEnd updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e(LS=%0.10f),Amin=%d,Amax=%d,IMIN=%d,IMAX=%d\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),LS,Amin,Amax,IMIN,IMAX);
	  if(LDEBUG){
	    sum = sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After PoutlierEnd updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS);
	  }
	}
	if((MDEBUG_ST >= 0 && m==MDEBUG_M) || (MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N+1)){
	  FL128 sum = 0.0;
	  for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++)
	    sum += LRsite[i]; 
	  printf("m=%d:After PoutlierEnd updates:LRsite[0..%d]=%0.10e(LS=%0.10f)\n",m,max(MDEBUG_TS+1,MDEBUG_ST),(double)sum,LS);
	}
	if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  if(newLPdelta && MDEBUG_TT < D[min(N,MDEBUG_TS)]){
	    FL128 sum = LRdelta[0];
	    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
	      sum += LRdelta[i];
	    FL128 sum3 = sum;
	    if(Amin <= MDEBUG_TS && MDEBUG_TS <= min(N,Amax) && MDEBUG_TT < D[MDEBUG_TS])
	      sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After PoutlierEnd updates:S=%d,D=%d,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e,LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3, LS);
	  }
	  if(LDEBUG){
	    FL128 sum = 0.0,sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After PoutlierEnd updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS);
	  }
	}
	fflush(stdout);
      }
      if(0 && MDEBUG && LDEBUG && !FIXADD){
	FL128 addsum = LRadd[0],addsum2 = LRadd2[0];
	for(int i = 1; i <= 52;i++){
	  addsum += LRadd[i];
	  addsum2 += LRadd2[i];
	}
	addsum -= addsum2;
	FL128 fragsum = LRfrag[0],fragsum2=LRfrag2[0];
	for(int i = 1; i <= 52; i++){
	  fragsum += LRfrag[i];
	  fragsum2 += LRfrag2[i];
	}
	fragsum -= fragsum2;
	if(DEBUG && (!(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190)||(VERB>=2&&m==1))){
	  printf("m=%d:After PoutlierEnd updates:sum(LRadd[0..52])=%0.17e,sum(LRfrag[0..52]=%0.17e\n",
		 m,(double)addsum,(double)fragsum);
	  fflush(stdout);
	  assert(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190);
	}
      }
      if(LDEBUG){/* check that cumulative value of LRsite[1..N+1] approximates 0 */
	FL128 LRsitecum = 0.0;
	for(int S = 1; S <= N+1; S++){
	  LRsitecum += LRsite[S];
	  if(VERB && (VERB>=2 || (DEBUG && !(S <= N || fabs((double)LRsitecum) <= LRcum * 1e-8)))){
	    printf("m=%d:S=%d,N=%d:LRsite[S]=%0.10e,LRsitecum=%0.10e,LRcum=%0.10e,LRfrag[0]=%0.10e,LRsitecum-LRcum=%0.10e,IMIN=%d,IMAX=%d\n",
		   m,S,N,(double)LRsite[S],(double)LRsitecum,LRcum,(double)LRfrag[0],(double)(LRsitecum-LRcum),IMIN,IMAX);
	    for(int s = 1; s <= N+1; s++)
	      printf("LRsite[%d]=%0.10e,LRfrag[%d]=%0.10e\n",s,(double)LRsite[s],s,(double)LRfrag[s]);
	    fflush(stdout);
	  }
	  if(DEBUG && S == N+1) assert(fabs((double)LRsitecum) <= LRcum * 1e-8);
	}
      }

      FL128b Fn = FN[0];
      FL128b InvFN = 1.0/Fn;

      if(DEBUG>=2 && newLPdelta){
	if(VERB>=3){
	  if(m==7 && Amin <= 12 && 12 <= Amax && D[12] >= 1 /* && !isfinite(LRdeltaD[12][0])*/)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[12]=%d,delta[12][0]= %0.4f,LRdeltaD[12][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[12],delta[12][0],LRdeltaD[12][0]);
	  if(m==56 && Amin <= 31 && 31  <= Amax && D[31] >= 1 /* && !isfinite(LRdeltaD[31][0] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[31]=%d,delta[31][0]= %0.4f,LRdeltaD[31][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[31],delta[31][0],LRdeltaD[31][0]);
	  if(m==20 && Amin <= 33 && 33 <= Amax && D[33] >= 11 /* && !isfinite(LRdeltaD[33][10] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[33]=%d,delta[33][10]= %0.4f,LRdeltaD[33][10]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[33],delta[33][10],LRdeltaD[33][10]);
	  fflush(stdout);
	}
	for(int S = Amin;S <= min(N,Amax); S++)
	  for(int d = 0; d < D[S]; d++)
	    assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
      }

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After Updates based on PoutlierEnd:wall time= %0.6f\n", m, tid, wtime() - wt);
	  fflush(stdout);
        }
      }

      /* Updates based on FR term (Line 1 of Equation 7 for LRdel) (Line 1 of Equation 10 for LRadd) (Line 1 of Equation 21 for LRdelta) */
      /* As with UR only compute cases where J >= M + 1-DELTA_XU (OR EXTEND && I >= N+1-DELTA_YU && J >= M+1 - DELTA_XU - EXTEND_MAX) */
      for(int I = IMIN; I <= IMAX; I++){
	FL128b ScaleL = PoutlierEndM1 * exp(LogScaleL[I]-LS);
	for(int K = 0; K <= Kmax[I]; K++){
	  double Yik = Yc(Y,I,K);
	  //	  double *ALik = &A.AL(I,K,0);
	  int J = (EXTEND && I >= N+1-DELTA_YU - 1/*for deleted site*/) ? max(Jmin[I],M+1 - DELTA_XU - EXTEND_MAX) : max(Jmin[I],M+1-DELTA_XU);
	  for(; J <= Jmax[I]; J++){
	    if(DEBUG>=2) assert(A.AL(I,K,J) >= 0.0);
	    FL128b AL = A.AL(I,K,J) * ScaleL;
	    if(DEBUG>=1+RELEASE/* HERE HERE >= 2 */ && !isfinite(AL)){
	      printf("\nm=%d,I=%d,K=%d,J=%d:ScaleL= %0.24Le, A.AL(I,K,J)= %0.17e, AL= %0.24Le\n",m,I,K,J,(long double)ScaleL, (double)A.AL(I,K,J), (long double)AL);
	      fflush(stdout);
	      assert(isfinite(AL));
	    }
	    int Umin = A.URmin(I,K,J);/* will be reduced further if delta[S][d] > 0 */
	    int URmax = A.URmax(I,K,J);

	    if(DEBUG>=2) assert(!EXTEND || J >= M+1-DELTA_XU || I >= N+1-DELTA_YU - 1/* for deleted site*/);
	    int Ub = min(N,A.Rij(I,K,J)) + 1;

	    int Umax = min(N,URmax) + 1;/* increase by one so deletion of site min(N,URmax) is handled */

	    if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	      printf("m=%d:I=%d,K=%d,J=%d:Umin=%d,URmax=%d,Umax=%d,Ub=%d,Rij=%d\n",m,I,K,J,Umin,URmax,Umax,Ub,A.Rij(I,K,J));
	      for(int i = Umin; i >= I; i--)
		printf("  Y[%d]=%0.3f\n",i,Y[i]);
	      fflush(stdout);
	    }

	    if(DEBUG>=2) assert(Umin > I);
	    if(DEBUG>=2) assert(Umin <= Umax);
	    if(DEBUG>=2) assert(URmax <= Umax);
	    if(DEBUG>=2) assert(URmax <= N+1);

	    if(LPdelG){	    /* expand LRdel[Dmin..Dmax] to include LRdel[I..Umax+1] */
	      if(DEBUG>=2) assert(Dmin <= I);
	      for(int i = Dmax+1; i <= min(N,Umax)+1; i++)
		LRdel[Dmax = i] = 0.0;
	    }

	    /* expand LRadd[0,Amin..Amax] to include LRadd[I..Umax] (also LRdelta[I..Umax])*/
	    if(DEBUG>=2) assert(Amin <= I);
	    int origAmax = Amax;
	    for(int i = Amax+1; i <= Umax; i++){
	      LRadd[Amax = i] = LRadd2[i] = LRdelta[i] = 0.0;
	      if(i <= N){
		if(myaddcnt > 0)
		  for(int t = 0; t < addcntI[i]; t++)
		    LRaddT[i][t] = 0.0;
		if(newLPdelta)
		  for(int d = 0; d < ((D[i] + VPAD) & (~VPAD)); d++)
		    LRdeltaD[i][d] = 0.0;
	      }
	    }

	    if(DEBUG>=1+RELEASE && Amax > origAmax + 100 && !(Amax < IMAX + 100 || Y[Amax] < Y[IMAX] + 500.0)){
	      printf("\nWARNING:m=%d/%d:I=%d,K=%d,J=%d,Umin=%d,Umax=%d,N=%d, IMIN=%d,IMAX=%d: Amin=%d,Amax=%d->%d, Y[IMIN]=%0.4f,Y[IMAX]=%0.4f,Y[Amin]=%0.4f,Y[Amax]=%0.4f,A.URmin=%d,A.URmax=%d,A.Rij=%d\n",
		     m,MD,I,K,J,Umin,Umax,N,IMIN,IMAX,Amin,origAmax,Amax,Y[IMIN],Y[IMAX],Y[Amin],Y[Amax], A.URmin(I,K,J),A.URmax(I,K,J),A.Rij(I,K,J));
	      printf("\t M= %d, X[J]= %0.4f, X[M]= %0.4f, X[M+1]= %0.4f\n", M,Xm[J],Xm[M],Xm[M+1]);
	      fflush(stdout);
	      //	      if(m==43070 && I==43589 && K==0 && J==36)	assert((Amin > IMIN - 100 || Y[Amin] > Y[IMIN] - 500.0) && (Amax < IMAX + 100 || Y[Amax] < Y[IMAX] + 500.0));
	    }

	    double Xj = Xm[M+1] - Xm[J];

	    /* The following loop is seperated for LRdeltaD[S][d] updates since the range is different */
	    if(newLPdelta && (!EXTEND || J >= M+1-DELTA_XU || I >= N+1-DELTA_YU)){/* FR/GE term */

	      /* U = N+1 .. URmax+1 */
	      for(int S = I; S < URmax; S++){
		/* LRdelta[S][d] : S range for FRA(Y[U-1]+delta[S][d],Y[U]+delta[S][d]) is I .. U-2 */
		for(int d = 0; d < D[S]; d++){
		  double del = delta[S][d];
		  if(del >= 0.0)
		    continue;
		  for(int U = URmax; ++U <= N+1;){
		    if(URANGE_FIX<2/* HERE : tighten this up */ && Y[URmax]-Y[U-1] <= del)
		      break;
		    double Yiu = Y[U]+del - Yik;
		    double FRu1 = (ERFC_FIX && Yiu < Xj /* U < Ub */) ? FRA<1>(I,K,J,U,U-1,Yik,Xm,Y,Y[U-1]+del,Y[U]+del,M,N,END_FIX ? rc : 1) 
		      : FRA<0>(I,K,J,U,U-1,Yik,Xm,Y,Y[U-1]+del,Y[U]+del,M,N,END_FIX ? rc : 1); 
		    FRu1 *= AL;
		    if(DEBUG>=2) assert(isfinite(FRu1));
		    if(DEBUG>=2) assert(Amin <= S && S <= min(N,Amax));
		    LRdeltaD[S][d] += FRu1;
		    if(DEBUG>=2) assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
		    if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && d == MDEBUG_TT && FRu1 >= MDEBUG_MIN){
		      FL128 sum = LRdelta[0];
		      for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
			sum += LRdelta[i];
		      FL128 sum3 = sum;
		      if(Amin <= MDEBUG_TS && MDEBUG_TS <= min(N,Amax) && MDEBUG_TT < D[MDEBUG_TS])
			sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
		      printf("m=%d:I=%d,K=%d,J=%d,U=%d:FRu1(Y[U-1,U]+=del)=%0.10e:S=%d,D=%d,del=%0.3f,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e)\n",
			     m,I,K,J,U,(double)FRu1,MDEBUG_TS,MDEBUG_TT,del,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3);
		      fflush(stdout);
		    }
		  }
		}
	      }

	      /* U = URmax .. I+1 */
	      for(int U = URmax; U > I; U--){
		/* LRdelta[S][d] : S range for FRA(Y[U-1]+delta[S][d],Y[U]+delta[S][d]) is I .. U-2 */
		for(int S = I; S < U-1; S++){
		  for(int d = 0; d < D[S]; d++){
		    double del = delta[S][d];
		    if((del > 0.0) ? (URANGE_FIX>=2/* HERE : tighten this up */ || Y[Umin-1]-Y[U] < del) : (U >= Umin)){
		      double Yiu = Y[U]+del - Yik;
		      double FRu1 = (ERFC_FIX && Yiu < Xj /* U < Ub */) ? FRA<1>(I,K,J,U,U-1,Yik,Xm,Y,Y[U-1]+del,Y[U]+del,M,N,END_FIX ? rc : 1)
			: FRA<0>(I,K,J,U,U-1,Yik,Xm,Y,Y[U-1]+del,Y[U]+del,M,N,END_FIX ? rc : 1); 
		      FRu1 *= AL;
		      if(DEBUG>=2) assert(isfinite(FRu1));
		      if(DEBUG>=2) assert(Amin <= S && S <= min(N,Amax));
		      LRdeltaD[S][d] += FRu1;
		      if(DEBUG>=2) assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
		      if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && d == MDEBUG_TT && FRu1 >= MDEBUG_MIN){
			FL128 sum = LRdelta[0];
			for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
			  sum += LRdelta[i];
			FL128 sum3 = sum;
			if(Amin <= MDEBUG_TS && MDEBUG_TS <= min(N,Amax) && MDEBUG_TT < D[MDEBUG_TS])
			  sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
			printf("m=%d:I=%d,K=%d,J=%d,U=%d:FRu1(Y[U-1,U]+=del)=%0.10e:S=%d,D=%d,del=%0.3f,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e)\n",
			       m,I,K,J,U,(double)FRu1,MDEBUG_TS,MDEBUG_TT,del,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3);
			fflush(stdout);
		      }
		    }
		  }
		}

		/* LRdelta[S][d] : S range for FRA(Y[U-1],Y[U]+delta) is U-1 */
		for(int d = 0; d < D[U-1]; d++){
		  double del = delta[U-1][d];
		  if((del > 0.0) ? (URANGE_FIX>=2/* HERE : tighten this up */ || Y[Umin-1]-Y[U] < del) : (U >= Umin)){
		    double Yiu = Y[U]+del - Yik;
		    double FRu1 = (ERFC_FIX && Yiu < Xj /* U < Ub */) ? FRA<1>(I,K,J,U,U-1,Yik,Xm,Y,Y[U-1],Y[U]+del,M,N,END_FIX ? rc : 1)
		      : FRA<0>(I,K,J,U,U-1,Yik,Xm,Y,Y[U-1],Y[U]+del,M,N,END_FIX ? rc : 1);
		    FRu1 *= AL;
		    if(DEBUG>=2) assert(isfinite(FRu1));
		    if(DEBUG>=2) assert(Amin <= U-1 && U-1 <= min(N,Amax));
		    LRdeltaD[U-1][d] += FRu1;
		    if(DEBUG>=2) assert(isfinite(LRdeltaD[U-1][d]) && LRdeltaD[U-1][d] >= 0.0);
		    if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && U-1==MDEBUG_TS && d == MDEBUG_TT && FRu1 >= MDEBUG_MIN){
		      FL128 sum = LRdelta[0];
		      for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
			sum += LRdelta[i];
		      FL128 sum3 = sum;
		      if(Amin <= MDEBUG_TS && MDEBUG_TS <= min(N,Amax) && MDEBUG_TT < D[MDEBUG_TS])
			sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
		      printf("m=%d:I=%d,K=%d,J=%d,U=%d:FRu1(Y[U]+=del)=%0.10e:S=%d,D=%d,del=%0.3f,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e)\n",
			     m,I,K,J,U,(double)FRu1,MDEBUG_TS,MDEBUG_TT,del,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3);
		      fflush(stdout);
		    }
		  }
		}
	      }/* U = Umax .. I+1 */

	      if(ISLINEAR && rc){ /* FRE/G term  : U == N+1 */
		int U = N+1;
		/* LRdelta[S][d] : S range for FRE(Y[N+1] += delta[S][d]) is I .. U-1 = I .. N */
		/* Note : If delta values for different S match each other, a faster algorithm is possible */
		double origYend = Y[N+1];
		for(int S = I; S <= N; S++){
		  for(int d = 0; d < D[S]; d++){
		    double del = delta[S][d];
		    if((del > 0.0) ? (U <= URmax) : (Y[URmax]-Y[U-1] > del)){
		      if(DEBUG>=2) assert(S < URmax);//NOTE: Hence this loop can be speeded up
		      Y[N+1] = origYend + del;
		      if(DEBUG>=2) assert(Amin <= S && S <= min(N,Amax));
		      LRdeltaD[S][d] += AL * FRE(I,K,J,N+1,N,Yik,Xm,Y,M);
		      if(DEBUG>=2) assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
		    }
		  }
		}
		Y[N+1] = origYend;
	      }/* rc && U==N+1 */
	    }/* newLPdelta */

	    for(int U = Umax; U >= Umin; U--){/* FR/GE terms (except LRdeltaD[S][d]) */
	      MFLOAT FRikj = (ERFC_FIX && U < Ub) ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)
			      : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
	      FL128b FRu1 = AL * FRikj;
	      if(DEBUG>=2 && !(isfinite(FRu1) && FRu1 >= 0.0)){
		printf("AL=%0.10e,FR=%0.10e,FRu1=%0.10e,x=%0.6f,y=%0.6f,b=%0.6f\n",(double)AL,FRikj,(double)FRu1,Xm[M+1]-Xm[J], Y[U] - Yik, (I==U-1) ? Y[U] - Yik : Y[U]-Y[U-1]);
		fflush(stdout);
		assert(isfinite(FRu1) && FRu1 >= 0.0);
	      }

	      if(DEBUG>=2) assert(isfinite(FRu1));

	      if(newLPdelta && (!EXTEND || J >= M+1-DELTA_XU || I >= N+1-DELTA_YU) && U <= N && U <= URmax){ /* LRdelta[S]: S range for FRu1 is U .. N */
		if(DEBUG>=2) assert(U==0 || (Amin <= U && U <= Amax));
		LRdelta[U] += FRu1;
		if(LDEBUG>=2) LRdelta[N+1] -= FRu1;

		if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && U<=MDEBUG_TS && FRu1 >= MDEBUG_MIN){
		  FL128 sum = LRdelta[0];
		  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
		    sum += LRdelta[i];
		  FL128 sum3 = sum;
		  if(Amin <= MDEBUG_TS && MDEBUG_TS <= min(N,Amax) && MDEBUG_TT < D[MDEBUG_TS])
		    sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
		  printf("m=%d:I=%d,K=%d,J=%d,U=%d:FRu1=%0.10e:S=%d,D=%d,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e)\n",
			 m,I,K,J,U,(double)FRu1,MDEBUG_TS,MDEBUG_TT,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3);
		  fflush(stdout);
		}
	      }

	      if(LPdelG){      /* LRdel[S]:  S range is min(N,U)+1 .. N  */
		if(DEBUG>=2) assert(Dmin <= min(N,U)+1  && min(N,U)+1 <= Dmax);
		LRdel[min(N,U)+1] += FRu1;
		if(DEBUG>=2) assert(isfinite(LRdel[min(N,U)+1]));
		if(LDEBUG>=2) LRdel[N+1] -= FRu1;
		if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && min(N,U)+1 <= MDEBUG_S && MDEBUG_S <= N && FRu1 >= MDEBUG_MIN){
		  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		    sum += LRdel[i];
		  printf("m=%d:I=%d,K=%d,J=%d,U=%d:FRu1=%0.10e(x=%0.4f,y=%0.4f,b=%0.4f,AL=%0.10e,ScaleL=%0.10e),LRdel[0..%d]=%0.10e\n",
			 m,I,K,J,U,(double)FRu1,Xm[M+1]-Xm[J],Y[U]-Yik, Y[U]-Y[U-1], A.AL(I,K,J),(double)ScaleL,MDEBUG_S,(double)sum);
		  fflush(stdout);
		}
	      }

	      if(LDEBUG && (!EXTEND || J >= M+1-DELTA_XU || I >= N+1-DELTA_YU) && U <= URmax){
		/* LRsite[S]: S range for FRu1 is I+1..N */
		LRsite[I+1] += FRu1;
		LRsite[N+1] -= FRu1;
		if(LDEBUG && DEBUG>=2){/* verify that LRsite[1..N+1] adds up to approximately 0.0 */
		  FL128 sum = 0.0;
		  for(int S = 1; S <= N+1; S++)
		    sum += LRsite[S];
		  if((DEBUG && !(fabs((double)sum) <= LRcum * 1e-8)) || (0 && m==MDEBUG_M && FRu1 >= MDEBUG_MIN && (I+1==STRACE || N+1==STRACE))){
		    printf("m=%d,I=%d,K=%d,J=%d,U=%d,N=%d,M=%d:FRu1=%0.10e,LRsite[1..N+1]=%0.10e,LRsite[%d]+=%0.10e,LRsite[%d]-=%0.10e\n",
			   m,I,K,J,U,N,M,(double)FRu1,(double)sum,I+1,(double)LRsite[I+1],N+1,(double)LRsite[N+1]);
		    fflush(stdout);
		    assert(fabs((double)sum) <= LRcum * 1e-8);
		  }
		}
		if((MDEBUG_M < 0 || m==MDEBUG_M) && ((MDEBUG_ST >= 0 && I+1 <= MDEBUG_ST && MDEBUG_ST <= N) ||
						     (MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && I+1 <= MDEBUG_TS+1 && MDEBUG_TS+1 <= N && N == MDEBUG_N+1)) && AL >= MDEBUG_MIN){
		  FL128 sum = 0.0;
		  for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++)
		    sum += LRsite[i];
		  int t = max(MDEBUG_TS+1,MDEBUG_ST);
		  printf("m=%d,I=%d,K=%d,J=%d,U=%d,N=%d,M=%d:FRu1=%0.10e(AL=%0.10e,FRikj=%0.10e),LRsite[1..%d]=%0.10e(LRsite[%d]=%0.10e,LRsite[%d]=%0.10e)\n",
			 m,I,K,J,U,N,M,(double)FRu1,(double)AL,FRikj,t,(double)sum,t,(double)LRsite[t],t+1,(double)LRsite[t+1]);
		}

		/* LRfrag[S] S range is I..N */
		LRfrag[I] += FRu1;
		LRfrag2[N+1] += FRu1;
		if(MDEBUG_M >= 0 && m==MDEBUG_M && ((N==MDEBUG_N-1 && I <= MDEBUG_F && MDEBUG_F <= N) || 
						    (N==MDEBUG_N && MDEBUG_T < 0 && MDEBUG_TT >= 0 && I <= MDEBUG_TS && MDEBUG_TS <= N)) && FRu1 >= MDEBUG_MIN){
		  FL128 sum = 0.0, sum2 = 0.0;
		  for(int i = 0; i <= max(MDEBUG_F,MDEBUG_TS); i++){
		    sum += LRfrag[i];
		    sum2 += LRfrag2[i];
		  }
		  printf("m=%d:I=%d,K=%d,J=%d,U=%d(Ub=%d):FRu1=%0.10e(x=%0.4f,y=%0.4f,b=%0.4f,AL=%0.10e,ScaleL=%0.10e),LRfrag[0..%d]=%0.10e\n",
			 m,I,K,J,U,Ub,(double)FRu1,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],A.AL(I,K,J),(double)ScaleL,max(MDEBUG_F,MDEBUG_TS),(double)(sum-sum2));
		  fflush(stdout);
		}
	      }

	      if(LPdelG && U >= I+2){
		if(U > I+2){
		  FL128 inc = FRu1 * InvFN;	    
		  /* LRdel[S]: S range = I+1 .. U-2 */
		  if(DEBUG>=2) assert(isfinite(inc));
		  if(DEBUG>=2) assert(Dmin <= I+1 && I+1 < U-1 && U-1 <= Dmax);
		  LRdel[I+1] += inc;
		  LRdel[U-1] -= inc;
		  if(DEBUG>=2) assert(isfinite(LRdel[U-1]));
		  if(DEBUG>=2) assert(isfinite(LRdel[I+1]));
		  if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && I+1 <= MDEBUG_S && MDEBUG_S <= U-2 && inc >= MDEBUG_MIN){
		    FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		    for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		      sum += LRdel[i];
		    MFLOAT FRikj = (ERFC_FIX && U < Ub) ? FR<1>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1)
			      : FR<0>(I,K,J,U,U-1,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-1],Xm,Y,M,N,END_FIX ? rc : 1);
		    printf("m=%d:I=%d,K=%d,J=%d,U=%d(Ub=%d):FRu1*InvFN=%0.10e(x=%0.4f,y=%0.4f,b=%0.4f,AL=%0.10e,FRikj=%0.10e,ScaleL=%0.10e),LRdel[0..%d]=%0.10e\n",
			   m,I,K,J,U,Ub,(double)inc,Xm[M+1]-Xm[J], Y[U]-Yik, Y[U]-Y[U-1], A.AL(I,K,J), FRikj, (double)ScaleL, MDEBUG_S,(double)sum);
		    fflush(stdout);
		  }
		}

		FL128b FRu2 = AL * ((ERFC_FIX && U < Ub) ? FR<1>(I,K,J,U,U-2,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-2],Xm,Y,M,N,END_FIX ? rc : 1)
				    : FR<0>(I,K,J,U,U-2,Yik,Xm[M+1]-Xm[J],Y[U]-Yik,Y[U]-Y[U-2],Xm,Y,M,N,END_FIX ? rc : 1));
		if(DEBUG>=2) assert(isfinite(FRu2) && FRu2 >= 0.0);
		/* LRdel[S]: S range is U-1 */
		if(DEBUG>=2) assert(Dmin <= U-1  && U <= Dmax);
		LRdel[U-1] += FRu2;
		LRdel[U] -= FRu2;
		if(DEBUG>=2) assert(isfinite(LRdel[U-1]));
		if(DEBUG>=2) assert(isfinite(LRdel[U]));
		if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && U-1 == MDEBUG_S && FRu2 >= MDEBUG_MIN && LPdelG){
		  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		    sum += LRdel[i];
		  printf("m=%d:I=%d,K=%d,J=%d,U=%d:FRu1=%0.10e,LRdel[0..%d]=%0.10e\n",m,I,K,J,U,(double)FRu1,MDEBUG_S,(double)sum);
		  fflush(stdout);
		}
	      }/* if(LPdelG .. ) */

	      if((!EXTEND || J >= M+1-DELTA_XU || I >= N+1-DELTA_YU + FIXADD) && U <= URmax){
		/* LRadd[S]: S range for FRu1 is U .. N */
		if(U <= N){
		  if(DEBUG>=2) assert(Amin <= U && U <= Amax);
		  LRadd[U] += FRu1;
		  if(LDEBUG>=2) LRadd2[N+1] += FRu1;
		}

		/* LRadd[S]: S range for FRu1*FN is I .. U-1 */
		if(DEBUG>=2) assert(U > I && U >= 1);
		FL128b inc = FRu1 * Fn;
		if(DEBUG>=2) assert(Amin <= I && U <= Amax);
		LRadd[I] += inc;
		LRadd2[U] += inc;

		if(myaddcnt > 0){ /* LRaddT[S][] : S range is U-1 */
		  double incf = inc;
		  for(int t = 0; t < addcntI[U-1]; t++){
		    if(DEBUG>=2) assert(Y[U-1] < addlocI[U-1][t] && addlocI[U-1][t] < Y[U]);
		    /*		    double FRu = (ERFC_FIX && U < Ub) ? 
		      (FRA<1>(I,K,J,U,U-1,Yik,Xm,Y,Y[U-1],addlocI[U-1][t],M,N,1) + Fn * FRA<1>(I,K,J,U,U-1,Yik,Xm,Y,addlocI[U-1][t],Y[U],M,N,END_FIX ? rc : 1))
		      : (FRA<0>(I,K,J,U,U-1,Yik,Xm,Y,Y[U-1],addlocI[U-1][t],M,N,1) + Fn * FRA<0>(I,K,J,U,U-1,Yik,Xm,Y,addlocI[U-1][t],Y[U],M,N,END_FIX ? rc : 1));*/
		    double FRuA = (ERFC_FIX && U < Ub) ? FRA<1>(I,K,J,U,U-1,Yik,Xm,Y,Y[U-1],addlocI[U-1][t],M,N,1) : FRA<0>(I,K,J,U,U-1,Yik,Xm,Y,Y[U-1],addlocI[U-1][t],M,N,1);
		    double FRuB = (ERFC_FIX && U < Ub) ? FRA<1>(I,K,J,U,U-1,Yik,Xm,Y,addlocI[U-1][t],Y[U],M,N,END_FIX ? rc : 1) : FRA<0>(I,K,J,U,U-1,Yik,Xm,Y,addlocI[U-1][t],Y[U],M,N,END_FIX ? rc : 1);
		    double FRu = FRuA + FN[0] * FRuB;
		    double FRu3 = AL * FRu;
		    if(DEBUG>=1+RELEASE /* HERE HERE >=2 */ && !(isfinite(FRu3) && FRu3 - incf >= -incf * (USE_MFLOAT ? 1e-4 : 2e-5) - (USE_MFLOAT ? 1e-16 : 1e-30))){
		      #pragma omp critical
		      {
			printf("\nm=%d,I=%d,K=%d,J=%d,U=%d(%d..%d),Ub=%d,N=%d,M=%d:Y[U-1]=%0.4f,Y[U]=%0.4f,t=%d,Ht=%0.4f:FRu1=%0.10e,FN[0]=%0.6f,FRu1*FN[0]=inc=%0.24Le,AL=%0.24Le,FRu=%0.17e,FRu3=%0.17e\n",
			       m,I,K,J,U,Umax,Umin,Ub,N,M,Y[U-1],Y[U],t,addlocI[U-1][t],(double)FRu1,FN[0],(long double)inc,(long double)AL,(double)FRu,(double)FRu3);
			printf("\t FRuA= %0.17e, FRuB= %0.17e, Yik= %0.8f, incf= %0.17e, Fn= %0.6Le, FRu3-incf= %0.17e\n", FRuA, FRuB, Yik, incf, (long double)Fn, FRu3-incf);
			fflush(stdout);
			assert(isfinite(FRu3) && isfinite(FRu3 - incf));
			//			assert(isfinite(FRu3) && FRu3 - inc >= -inc * (USE_MFLOAT ? 1e-4 : 2e-5) - (USE_MFLOAT ? 1e-16 : 1e-30));
			//			assert(0);
		      }
		    }
		    if(DEBUG>=1+RELEASE/* HERE HERE >=2 */) assert(Amin <= U-1 && U-1 <= min(N,Amax));
		    LRaddT[U-1][t] += max(((double)0.0),FRu3 - incf);
		    if(DEBUG>=1+RELEASE/* HERE HERE >=2 */) assert(isfinite(LRaddT[U-1][t]) && LRaddT[U-1][t] >= 0.0);
		  }
		}
	      }

	      if(rc && U==N+1){ /* FRE/G term  : U == N+1 */
		FL128b FRu1 = AL * FRE(I,K,J,N+1,N,Yik,Xm,Y,M);		
		if(DEBUG>=2) assert(isfinite(FRu1) && FRu1 >= 0.0);
		/* LRdel[S]: S range = U .. N (empty) */

		if(LDEBUG && (!EXTEND || J >= M+1-DELTA_XU || I >= N+1-DELTA_YU) && U <= URmax){
		  /* LRsite[S]: S range is I+1..N */
		  LRsite[I+1] += FRu1;
		  LRsite[N+1] -= FRu1;
		  /* LRfrag[] S range is I..N */
		  LRfrag[I] += FRu1;
		  LRfrag2[N+1] += FRu1;
		}

		/* LRadd[S]: S range for FREu1 is U .. N (empty) */

		if(LPdelG && U >= I+2){
		  if(DEBUG>=2) assert(Dmax >= N+1);
		  if(U > I+2){
		    FL128b inc = FRu1 * InvFN;		
		    /* LRdel[] S range = I+1 .. U-2 = I+1 .. N-1 */
		    if(DEBUG>=2) assert(isfinite(inc));
		    if(DEBUG>=2) assert(Dmin <= I+1 && I+1 < N && N <= Dmax);
		    LRdel[I+1] += inc;
		    LRdel[N] -= inc;
		    if(DEBUG>=2) assert(isfinite(LRdel[I+1]));
		    if(DEBUG>=2) assert(isfinite(LRdel[N]));
		  }

		  FL128b FRu2 = AL * FRE(I,K,J,N+1,N-1,Yik,Xm,Y,M);
		  if(DEBUG>=2) assert(isfinite(FRu2) && FRu2 >= 0.0);
		  /* LRdel[] S range = U-1 = N */
		  if(DEBUG>=2) assert(Dmin <= N  && N+1 <= Dmax);
		  LRdel[N] += FRu2;
		  LRdel[N+1] -= FRu2;
		  if(DEBUG>=2) assert(isfinite(LRdel[N]));
		  if(DEBUG>=2) assert(isfinite(LRdel[N+1]));
		} /* if (LPdelG && U >= I+1) */

		if((!EXTEND || J >= M+1-DELTA_XU || I >= N+1-DELTA_YU + FIXADD) && U <= URmax){
		  FL128b inc = FRu1 * Fn;
		  /* LRadd[S]: S range for FREu1 * FN is I .. U-1 = I .. N (includes LRaddT[S] terms for S range = U-1 since it does not depend on site location) */
		  if(DEBUG>=2) assert(Amin <= I && I <= N && N+1 <= Amax);
		  LRadd[I] += inc;
		  LRadd2[N+1] += inc;
		}
		/* LRaddT[S][] : S range is U-1 : term is included in LRadd[] above since it does not depend on new site location */
	      }/* rc && U==N+1 */
	    } /* U = Umax .. Umin */
	  }
	}
      }

      if(VERB>=2 && newLPdelta && rverb){
	printf("mprobeval:m=%d,Amin=%d,Amax=%d:After FR updates:time=%0.6f\n",m,Amin,Amax,wtime());
	fflush(stdout);
      }

      if(MDEBUG_M >= 0){
	if(MDEBUG_S >= 0 && m==MDEBUG_M && N==MDEBUG_N && LPdelG){
	  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
	    sum += LRdel[i];
	  printf("m=%d:After FR updates:LRdel[0..%d]=%0.10e(LS=%0.10f,LL=%0.10f)\n",m,MDEBUG_S,(double)sum,LS,LS+log((double)sum));
	  if(VERB>=2){
	    printf("LRdel[0]=%0.10e\n",(double)LRdel[0]);
	    sum = LRdel[0];
	    for(int i = max(1,Dmin); i <= Dmax;i++){
	      sum += LRdel[i];
	      printf("LRdel[%d]=%0.10e,LRdel[0..%d]=%0.10e\n",i,(double)LRdel[i],i,(double)sum);
	    }
	    fflush(stdout);
	  }
	}
	if(LDEBUG && MDEBUG_F >= 0 && m==MDEBUG_M && N==MDEBUG_N-1){
	  FL128 sum = 0.0,sum2=0.0;
	  for(int i = 0; i <= MDEBUG_F; i++){
	    sum += LRfrag[i];
	    sum2 += LRfrag2[i];
	  }
	  printf("m=%d:After FR updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f,LL=%0.10f)\n",m,MDEBUG_F,(double)sum,(double)sum2,(double)(sum-sum2),LS,LS+log((double)(sum-sum2)));
	}
	if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  FL128 sum = LRadd[0], sum2 = LRadd2[0];
	  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
	    sum += LRadd[i];
	    sum2 += LRadd2[i];
	  }
	  FL128 sum3 = sum;
	  if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
	    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After FR updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e,LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2),LS);
	  } else
	    printf("m=%d:After FR updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e(LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),LS);
	  if(LDEBUG){
	    sum = sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After FR updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f,LL=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS,LS+log((double)(sum-sum2)));
	  }
	}
	if((MDEBUG_ST >= 0 && m==MDEBUG_M) || (MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N+1)){
	  FL128 sum = 0.0;
	  for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++)
	    sum += LRsite[i]; 
	  printf("m=%d:After FR updates:LRsite[0..%d]=%0.10e(LS=%0.10f)\n",m,max(MDEBUG_TS+1,MDEBUG_ST),(double)sum,LS);
	}
	if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  if(newLPdelta && MDEBUG_TT < D[min(N,MDEBUG_TS)]){
	    FL128 sum = LRdelta[0];
	    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
	      sum += LRdelta[i];
	    FL128 sum3 = sum;
	    if(Amin <= MDEBUG_TS && MDEBUG_TS <= min(N,Amax) && MDEBUG_TT < D[MDEBUG_TS])
	      sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After FR updates:S=%d,D=%d,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e,LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3,LS);
	  }
	  if(LDEBUG){
	    FL128 sum = 0.0,sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After FR updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e,LS=%0.10f\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS);
	  }
	}
	fflush(stdout);
      }
      if(0 && MDEBUG && LDEBUG && !FIXADD){
	FL128 addsum = LRadd[0],addsum2 = LRadd2[0];
	for(int i = 1; i <= 52;i++){
	  addsum += LRadd[i];
	  addsum2 += LRadd2[i];
	}
	addsum -= addsum2;
	FL128 fragsum = LRfrag[0],fragsum2 = LRfrag2[0];
	for(int i = 1; i <= 52; i++){
	  fragsum += LRfrag[i];
	  fragsum2 += LRfrag2[i];
	}
	fragsum -= fragsum2;
	if(DEBUG && (!(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190) || (VERB>=2&&m==1))){
	  printf("m=%d:After FR updates:sum(LRadd[0..52])=%0.17e,sum(LRfrag[0..52]=%0.17e\n",
		 m,(double)addsum,(double)fragsum);
	  fflush(stdout);
	  assert(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190);
	}
      }
      if(LDEBUG){/* check that cumulative value of LRsite[1..N+1] approximates 0 */
	FL128 LRsitecum = 0.0;
	for(int S = 1; S <= N+1; S++){
	  LRsitecum += LRsite[S];
	  if(VERB && (VERB>=2 || (DEBUG && !(S <= N || fabs((double)LRsitecum) <= LRcum * 1e-8)))){
	    printf("m=%d:S=%d,N=%d:LRsite[S]=%0.10e,LRsitecum=%0.10e,LRcum=%0.10e,LRfrag[0]=%0.10e,LRsitecum-LRcum=%0.10e,IMIN=%d,IMAX=%d\n",
		   m,S,N,(double)LRsite[S],(double)LRsitecum,LRcum,(double)LRfrag[0],(double)(LRsitecum-LRcum),IMIN,IMAX);
	    for(int s = 1; s <= N+1; s++)
	      printf("LRsite[%d]=%0.10e,LRfrag[%d]=%0.10e\n",s,(double)LRsite[s],s,(double)LRfrag[s]);
	    fflush(stdout);
	  }
	  if(DEBUG && S == N+1) assert(fabs((double)LRsitecum) <= LRcum * 1e-8);
	}
      }

      if(0 && MDEBUG_S >= 0 && m==MDEBUG_M && N==MDEBUG_N && LPdelG){
	FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
	  sum += LRdel[i];
	printf("m=%d:Before FL updates:LRdel[0..%d]=%0.10e,Dmin=%d,Dmax=%d\n",m,MDEBUG_S,(double)sum,Dmin,Dmax);
	if(VERB>=2){
	  printf("LRdel[0]=%0.10e\n",(double)LRdel[0]);
	  sum = LRdel[0];
	  for(int i = max(1,Dmin); i <= Dmax;i++){
	    sum += LRdel[i];
	    printf("LRdel[%d]=%0.10e,LRdel[0..%d]=%0.10e\n",i,(double)LRdel[i],i,(double)sum);
	  }
	  fflush(stdout);
	}
      }

      if(DEBUG>=2 && newLPdelta){
	if(VERB>=3){
	  if(m==7 && Amin <= 12 && 12 <= Amax && D[12] >= 1 /* && !isfinite(LRdeltaD[12][0])*/)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[12]=%d,delta[12][0]= %0.4f,LRdeltaD[12][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[12],delta[12][0],LRdeltaD[12][0]);
	  if(m==56 && Amin <= 31 && 31  <= Amax && D[31] >= 1 /* && !isfinite(LRdeltaD[31][0] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[31]=%d,delta[31][0]= %0.4f,LRdeltaD[31][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[31],delta[31][0],LRdeltaD[31][0]);
	  if(m==20 && Amin <= 33 && 33 <= Amax && D[33] >= 11 /* && !isfinite(LRdeltaD[33][10] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[33]=%d,delta[33][10]= %0.4f,LRdeltaD[33][10]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[33],delta[33][10],LRdeltaD[33][10]);
	  fflush(stdout);
	}
	for(int S = Amin;S <= min(N,Amax); S++)
	  for(int d = 0; d < D[S]; d++)
	    assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
      }

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After Updates based on FR term:wall time= %0.6f\n", m, tid, wtime() - wt);
	  fflush(stdout);
        }
      }

      /* Updates based on FL term (Line 2 of Equation 7 for LRdel and Equation 10 for LRadd and Equation ?? for LRdelta) */
      /* As with UL only compute cases where J <= DELTA_XU (OR EXTEND && I-K <= DELTA_YU && J <= DELTA_XU + EXTEND_MAX) */
      for(int I = IMIN; I <= IMAX; I++){
	FL128b ScaleR = PoutlierEndM1 * exp(LogScaleR[I]-LS);
	for(int K = 0; K <= Kmax[I]; K++){
	  double Yik = Yc(Y,I,K);
	  int jmax = (EXTEND && I-K <= DELTA_YU + 1/*For deleted site*/) ? min(Jmax[I],DELTA_XU + EXTEND_MAX) : min(Jmax[I],DELTA_XU);
	  int J = Jmin[I];
	  for(; J <= jmax; J++){
	    if(DEBUG>=2) assert(A.AR(I,K,J) >= 0.0);
	    FL128b AR = A.AR(I,K,J) * ScaleR;
	    if(DEBUG>=2) assert(isfinite(AR));
	    int ULmin = A.ULmin(I,K,J);
	    int Umin = max(1,ULmin) - 1;/* decrease by one so deletion of site ULmin is handled */
	    int Umax = A.ULmax(I,K,J);/* will be increased, when needed, to handle size interval increases */

	    if(DEBUG>=2) assert(!EXTEND || J <= DELTA_XU || I-K <= DELTA_YU + 1/* for deleted site*/);
	    int Ub = max(1,A.Lij(I,K,J)) - 1;

	    if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && 0 <= MDEBUG_TS && (MDEBUG_TS < I-K || MDEBUG_TS < Umax)){
	      printf("m=%d,I=%d,K=%d,J=%d:Umin=%d,ULmin=%d,Umax=%d\n",m,I,K,J,Umin,ULmin,Umax);
	      for(int i = Umax; i <= I-K; i++)
		printf("Y[%d]=%0.3f\n",i,Y[i]);
	      fflush(stdout);
	    }


	    if(DEBUG>=2 && !(A.ULmin(I,K,J) >= 0 && A.ULmin(I,K,J) < I-K)){
	      printf("m=%d,I=%d(%d..%d),K=%d(%d..%d),J=%d(%d..%d):A.ULmin(I,K,J)=%d,ULmax=%d,\n",
		     m,I,IMIN,IMAX,K,0,Kmax[I],J,Jmin[I],jmax,A.ULmin(I,K,J),A.ULmax(I,K,J));
	      fflush(stdout);
	      assert(A.ULmin(I,K,J) >= 0 && A.ULmin(I,K,J) < I-K);
	    }
	    if(DEBUG>=2) assert(A.ULmin(I,K,J) <= A.ULmax(I,K,J));
	    if(DEBUG>=2 && !(A.ULmax(I,K,J) >= 0 && A.ULmax(I,K,J) < I-K)){
	      printf("I=%d,K=%d,J=%d:A.ULmax(I,K,J)=%d,Kmax[I]=%d,Jmin[I]=%d,Jmax[I]=%d\n",I,K,J,A.ULmax(I,K,J),Kmax[I],Jmin[I],Jmax[I]);
	      fflush(stdout);
	      assert(A.ULmax(I,K,J) >= 0 && A.ULmax(I,K,J) < I-K);
	    }

	    if(LPdelG){ /* expand LRdel[Dmin..Dmax] to include LRdel[Umin .. I-K] */
	      if(DEBUG>=2) assert(Dmax >= I-K);
	      for(int i = Dmin-1; i >= max(1,Umin); i--)
		LRdel[Dmin = i] = 0.0;
	    }

	    /* expand LRadd[0,Amin..Amax] to include LRadd[Umin .. I-K] (also LRdelta[Umin .. I-K]) */
	    if(DEBUG>=2)assert(Amax >= I-K);
	    int origAmin = Amin;
	    for(int i = Amin-1; i >= Umin; i--){
	      if((Amin = i) > 0){
		LRadd[Amin] = LRadd2[Amin] = 0.0;
		LRdelta[Amin] = 0.0;
	      }
	      if(DEBUG>=2) assert(i <= N);
	      if(myaddcnt > 0)
		for(int t = 0; t < addcntI[i]; t++)
		  LRaddT[i][t] = 0.0;
	      if(newLPdelta)
		for(int d = 0; d < ((D[i] + VPAD) & (~VPAD)); d++)
		  LRdeltaD[i][d] = 0.0;
	    }
	    if(DEBUG>=1+RELEASE && Amin < origAmin - 100 && !(Amin > IMIN - 100 || Y[Amin] > Y[IMIN] - 500.0)){
	      printf("WARNING:m=%d/%d: I=%d,K=%d,J=%d,IMIN=%d,IMAX=%d: A.Lij(I,K,J)=%d,A.ULmin(I,K,J)=%d,A.ULmax(I,K,J)=%d,Umin=%d Amin=%d->%d,Amax=%d\n",
		     m,MD,I,K,J,IMIN,IMAX,A.Lij(I,K,J),A.ULmin(I,K,J),A.ULmax(I,K,J),Umin,origAmin,Amin,Amax);
	      fflush(stdout);
	      //	      assert((Amin > IMIN - 100 || Y[Amin] > Y[IMIN] - 500.0) && (Amax < IMAX + 100 || Y[Amax] < Y[IMAX] + 500.0));
	    }

	    double IvarX;
	    double FLj = FLJ(J,Xm[J],IvarX) * AR;
	    if(DEBUG>=2) assert(isfinite(FLj) && FLj >= 0.0);

	    if(DEBUG>=2) assert(Umin <= ULmin && ULmin <= Umax);

	    /* The following loop is seperated for LRdeltaD[S][d] updates since the S range is different */
	    if(newLPdelta && (!EXTEND || J <= DELTA_XU || I-K <= DELTA_YU)){/* FL/GE terms for LRdeltaD[S][d] */
	      MFLOAT Xj = Xm[J];


	      /* U = 0 .. ULmin-1 */
	      for(int S = ULmin; S < I-K; S++){
		/* Case 2 in last line of Equation 21 with I-K > S > U :
		   LRdelta[S][d] : S range for FLj * FLIKA(Y[U] -= delta[S][d],Y[U+1] -= delta[S][d]) is U+1 .. I-K-1 */
		for(int d = 0; d < D[S]; d++){
		  double del = delta[S][d];
		  if(del >= 0.0)
		    continue;
		  for(int U = ULmin; --U >= 0;){
		    if(URANGE_FIX<2/* HERE : tighten this up */ && Y[U+1]-Y[ULmin] <= del)
		      break;
		    if(DEBUG>=2) assert(0 <= Y[U+1]-Y[U] && Y[U+1]-Y[U] <= Yik-Y[U]+del);
		    MFLOAT2 Yiu = Yik - Y[U] + del;
		    double FLu1 = (ERFC_FIX && Yiu < Xj /* U > Ub */) ? FLIKA<1>(I,K,J,U,U+1,Yik,Xj,Y,Y[U]-del,Y[U+1]-del, Yiu, Y[U+1]-Y[U], IvarX, FLj, END_FIX ? lc : 1)
		      : FLIKA<0>(I,K,J,U,U+1,Yik,Xj,Y,Y[U]-del,Y[U+1]-del, Yiu, Y[U+1]-Y[U], IvarX, FLj, END_FIX ? lc : 1);
		    if(!FLIK_FIX)
		      FLu1 *= FLj;

		    if(DEBUG>=2) assert(isfinite(FLu1));
		    if(DEBUG>=2) assert(Amin <= S && S <= min(N,Amax));
		    LRdeltaD[S][d] += FLu1;
		    if(DEBUG>=2) assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
		    if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S == MDEBUG_TS && d == MDEBUG_TT && FLu1 >= MDEBUG_MIN){
		      FL128 sum = LRdelta[0];
		      for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
			sum += LRdelta[i];
		      FL128 sum3 = sum + LRdeltaD[MDEBUG_TS][MDEBUG_TT];
		      printf("m=%d,I=%d,K=%d,J=%d,U=%d:FLu1=%0.10e,S=%d,d=%d:LRdelta[0..S]=%0.10e,LRdeltaD[S][d]=%0.10e(sum=%0.10e)\n",
			     m,I,K,J,U,(double)FLu1,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)LRdeltaD[MDEBUG_TS][MDEBUG_TT],(double)sum3);
		    }
		  }
		}
		
		/* NOTE : Case 1 in last line of Equation 21 with S == U : ignored since original terms with U < ULmin were ignored (negligible) */
	      }

	      for(int U = ULmin; U < I-K; U++){
		/* Case 2 in last line of Equation 21 with I-K > S > U :
		   LRdelta[S][d] : S range for FLj * FLIKA(Y[U] -= delta[S][d],Y[U+1] -= delta[S][d]) is U+1 .. I-K-1 */
		for(int S = U+1; S < I-K; S++){
		  for(int d = 0; d < D[S]; d++){
		    double del = delta[S][d];
		    if((del > 0.0) ? (URANGE_FIX>=2/* HERE : tighten this up */ || Y[U]-Y[Umax+1] < del) : (U <= Umax)){
		      if(DEBUG>=2) assert(0 <= Y[U+1]-Y[U] && Y[U+1]-Y[U] <= Yik-Y[U]+del);
		      MFLOAT2 Yiu = Yik - Y[U] + del;
		      double FLu1 = (ERFC_FIX && Yiu < Xj /*U > Ub */) ? FLIKA<1>(I,K,J,U,U+1,Yik,Xj,Y,Y[U]-del,Y[U+1]-del, Yik-Y[U]+del, Y[U+1]-Y[U], IvarX, FLj, END_FIX ? lc : 1)
			: FLIKA<0>(I,K,J,U,U+1,Yik,Xj,Y,Y[U]-del,Y[U+1]-del, Yik-Y[U]+del, Y[U+1]-Y[U], IvarX, FLj, END_FIX ? lc : 1);
		      if(!FLIK_FIX)
			FLu1 *= FLj;

		      if(DEBUG>=2) assert(Amin <= S && S <= min(N,Amax));
		      if(DEBUG>=2) assert(isfinite(FLu1));
		      LRdeltaD[S][d] += FLu1;
		      if(DEBUG>=2) assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
		      if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S == MDEBUG_TS && d == MDEBUG_TT && FLu1 >= MDEBUG_MIN){
			FL128 sum = LRdelta[0];
			for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
			  sum += LRdelta[i];
			FL128 sum3 = sum + LRdeltaD[MDEBUG_TS][MDEBUG_TT];
			printf("m=%d,I=%d,K=%d,J=%d,U=%d:FLu1=%0.10e,S=%d,d=%d:LRdelta[0..S]=%0.10e,LRdeltaD[S][d]=%0.10e(sum=%0.10e)\n",
			       m,I,K,J,U,(double)FLu1,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)LRdeltaD[MDEBUG_TS][MDEBUG_TT],(double)sum3);
		      }
		    }
		  }
		}

		/* Case 1 in last line of Equation 21 with S == U :
		   LRdelta[S][d] : S range for FLj * FLIKA(Y[U] -= delta[S][d],Y[U+1]) is U */
		for(int d = 0; d < D[U]; d++){
		  double del = delta[U][d];
		  if((del > 0.0) ? (URANGE_FIX>=2/* HERE : tighten this up */ || Y[U]-Y[Umax+1] < del) : (U <= Umax)){
		    if(DEBUG>=2) assert(0 <= Y[U+1]-Y[U]+del && Y[U+1]-Y[U]+del <= Yik-Y[U]+del);
		    MFLOAT2 Yiu = Yik-Y[U]+del;
		    double FLu1 = (ERFC_FIX && Yiu < Xj /*U > Ub*/) ? FLIKA<1>(I,K,J,U,U+1,Yik,Xj,Y,Y[U]-del,Y[U+1], Yiu,Y[U+1]-Y[U]+del,IvarX, FLj, END_FIX ? lc : 1)
		      : FLIKA<0>(I,K,J,U,U+1,Yik,Xj,Y,Y[U]-del,Y[U+1], Yiu, Y[U+1]-Y[U]+del,IvarX, FLj, END_FIX ? lc : 1);
		    if(!FLIK_FIX)
		      FLu1 *= FLj;

		    if(DEBUG>=2) assert(Amin <= U && U <= min(N,Amax));
		    if(DEBUG>=2) assert(isfinite(FLu1));
		    LRdeltaD[U][d] += FLu1;
		    if(DEBUG>=2) assert(isfinite(LRdeltaD[U][d]) && LRdeltaD[U][d] >= 0.0);
		    if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && U == MDEBUG_TS && d == MDEBUG_TT && FLu1 >= MDEBUG_MIN){
		      FL128 sum = LRdelta[0];
		      for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
			sum += LRdelta[i];
		      FL128 sum3 = sum + LRdeltaD[MDEBUG_TS][MDEBUG_TT];
		      printf("m=%d,I=%d,K=%d,J=%d,U=%d:FLu1=%0.10e,S=%d,d=%d:LRdelta[0..S]=%0.10e,LRdeltaD[S][d]=%0.10e(sum=%0.10e)\n",
			     m,I,K,J,U,(double)FLu1,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)LRdeltaD[MDEBUG_TS][MDEBUG_TT],(double)sum3);
		    }
		  }
		}
	      }/* U = ULmin .. I-K-1 */

	      if(ISLINEAR && lc){ /* FLE/G term : U == 0 */
		int U = 0;
		/* LRdelta[S][d] : S range for FLj * FLEIK(Y[0] = -delta[S][d]) is U .. I-K-1 = 0 .. I-K-1 */
		for(int S = 0; S < I-K; S++){
		  for(int d = 0; d < D[S]; d++){
		    double del = delta[S][d];
		    if((del > 0.0) ? (U >= ULmin) : (Y[1]-Y[ULmin] > del)){
		      if(DEBUG>=2) assert(S >= ULmin);// HERE : possibility to speed up this loop
		      Y[0] = -delta[S][d];
		      if(DEBUG>=2) assert(Amin <= S && S <= min(N,Amax));
		      MFLOAT FLu1 = FLEIK(I,K,J,0,1,Yik,Xm[J], FLj);
		      if(!FLIK_FIX)
			FLu1 *= FLj;
		      LRdeltaD[S][d] += FLu1;
		      if(DEBUG>=2) assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
		    }
		  }
		}
		Y[0] = 0.0;
	      } /* lc && U==0 */
	    } // FL/GE terms for LRdeltaD[S][d]

	    for(int U = Umin; U <= Umax; U++){/* FL/GE terms (except LRdeltaD[S][d]) */
	      FL128b FLu1 = (ERFC_FIX && U > Ub) ? FLIKd<1>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
		: FLIKd<0>(I,K,J,U,U+1,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);
	      if(!FLIK_FIX)
		FLu1 *= FLj;

	      if(DEBUG>=2 && !(isfinite(FLu1) && FLu1 >= 0.0)){
 		printf("\nm=%d,I=%d,K=%d,J=%d,U=%d(%d..%d):FLu1= %0.10Le, Ub=%d, Yik=%0.4f,Xm[J]=%0.4f,Yik-Y[U]=%0.6f,Y[U+1]-Y[U]=%0.4f,FLj=%0.10e(IvarX=%0.10e,AR(I,K,J)=%0.10Le,ScaleR=%0.10Le)\n",
		       m,I,K,J,U,Umin,Umax,(long double)FLu1,Ub,Yik,Xm[J],Yik-Y[U],Y[U+1]-Y[U],FLj,IvarX,(long double)A.AR(I,K,J),(long double)ScaleR);
		fflush(stdout);
		assert(isfinite(FLu1) && FLu1 >= 0.0);
	      }

	      if(newLPdelta && (!EXTEND || J <= DELTA_XU || I-K <= DELTA_YU) && U > 0 && U >= ULmin){
		/* LRdelta[S]: S range for FLu1 is 0 .. U-1 */
		LRdelta[0] += FLu1;
		if(DEBUG>=2) assert(Amin <= U && U <= min(N,Amax));
		LRdelta[U] -= FLu1;
		if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && 0 <= MDEBUG_TS && MDEBUG_TS <= U-1 && FLu1 >= MDEBUG_MIN){
		  FL128 sum = LRdelta[0];
		  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
		    sum += LRdelta[i];
		  FL128 sum3 = sum + LRdeltaD[MDEBUG_TS][MDEBUG_TT];
		  printf("m=%d,I=%d,K=%d,J=%d,U=%d:FLu1=%0.10e,S=%d,d=%d:LRdelta[0..S]=%0.10e,LRdeltaD[S][d]=%0.10e(sum=%0.10e)\n",
			 m,I,K,J,U,(double)FLu1,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)LRdeltaD[MDEBUG_TS][MDEBUG_TT],(double)sum3);
		}
	      }

	      if(LPdelG){ /* LRdel[S]: S range for FLu1 is 1 .. max(1,U)-1 */
		if(LDEBUG||FIXLRDEL){
		  LRdel[0] += FLu1;
		  if(DEBUG>=2) assert(isfinite(LRdel[0]));
		}
		if(DEBUG>=2) assert(Dmin <= max(1,U) && max(1,U) <= Dmax);
		LRdel[max(1,U)] -= FLu1;
		if(DEBUG>=2) assert(isfinite(LRdel[max(1,U)]));
		if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && 1 <= MDEBUG_S && MDEBUG_S <= max(1,U)-1 && FLu1 >= MDEBUG_MIN){
		  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		    sum += LRdel[i];
		  printf("m=%d:I=%d,K=%d,J=%d,U=%d:FLu1=%0.10e,LRdel[0..%d]=%0.10e,Dmin=%d,Dmax=%d\n",m,I,K,J,U,(double)FLu1,MDEBUG_S,(double)sum,Dmin,Dmax);
		  if(DEBUG && sum < 0.0){
		    printf("LRdel[0]=%0.10e\n",(double)LRdel[0]);
		    FL128 sum2 = LRdel[0];
		    for(int i = max(1,Dmin); i <= Dmax;i++){
		      sum2 += LRdel[i];
		      printf("LRdel[%d]=%0.10e,LRdel[0..%d]=%0.10e\n",i,(double)LRdel[i],i,(double)sum2);
		    }
		    printf("MDEBUG_M=%d,MDEBUG_N=%d,MDEBUG_S=%d,MDEBUG_F=%d\n",MDEBUG_M,MDEBUG_N,MDEBUG_S,MDEBUG_F);
		    fflush(stdout);
		  
		    assert(sum >= 0.0);
		  }
		  fflush(stdout);
		}
	      } /* if(LPdelG) */

	      if(LDEBUG && (!EXTEND || J <= DELTA_XU || I-K <= DELTA_YU) && U >= ULmin){
		/* LRsite[S]: S range is 1 .. I-K-1 */
		LRsite[1] += FLu1;
		LRsite[I-K] -= FLu1;
		if(LDEBUG && DEBUG>=2){/* verify that LRsite[1..N+1] adds up to approximately 0.0 */
		  FL128 sum = 0.0;
		  for(int S = 1; S <= N+1; S++)
		    sum += LRsite[S];
		  if((DEBUG && !(fabs((double)sum) <= LRcum * 1e-8)) ||(0 && m==MDEBUG_M && FLu1 > MDEBUG_MIN && (1==STRACE || I-K==STRACE))){
		    printf("m=%d,I=%d,K=%d,J=%d,U=%d,N=%d,M=%d:FLu1=%0.10e,LRsite[1..N+1]=%0.10e,LRsite[%d]+=%0.10e,LRsite[%d]-=%0.10e\n",
			   m,I,K,J,U,N,M,(double)FLu1,(double)sum,1,(double)LRsite[1],I-K,(double)LRsite[I-K]);
		    fflush(stdout);
		    assert(fabs((double)sum) <= LRcum * 1e-8);
		  }
		}
		if((MDEBUG_M < 0 || m==MDEBUG_M) && ((MDEBUG_ST>=0 && 1 <= MDEBUG_ST && MDEBUG_ST <= I-K-1) ||
						     (MDEBUG_TS>=0 && MDEBUG_TT>=0 && m==MDEBUG_M && 1 <= MDEBUG_TS+1 && MDEBUG_TS+1 <= I-K-1 && N==MDEBUG_N+1)) && FLu1 >= MDEBUG_MIN){
		  FL128 sum = 0.0;
		  for(int i = 1; i <= max(MDEBUG_ST,MDEBUG_TS+1); i++)
		    sum += LRsite[i];
		  int t = max(MDEBUG_ST,MDEBUG_TS+1);
		  printf("m=%d,I=%d,K=%d,J=%d,U=%d(%d..%d,%d),N=%d,M=%d:FLu1=%0.10e,LRsite[1..%d]=%0.10e(LRsite[%d]=%0.10e,LRsite[%d]=%0.10e)\n",
			 m,I,K,J,U,Umin,Umax,Ub,N,M,(double)FLu1,t,(double)sum,t,(double)LRsite[t],t+1,(double)LRsite[t+1]);
		  //		    printf("\t FLj=%0.10e,IvarX=%0.10e,Yik=%0.8f,Xm[J]=%0.8f(AR=%0.10e,ScaleR=%0.10e)\n",(double)FLj,IvarX,Yik,Xm[J],A.AR(I,K,J),(double)ScaleR);
		}

		/* LRfrag[S]: S range is 0 .. I-K-1 */
		LRfrag[0] += FLu1;
		LRfrag2[I-K] += FLu1;
		if(MDEBUG_M >= 0 && m==MDEBUG_M && ((N==MDEBUG_N-1 && 0 <= MDEBUG_F && MDEBUG_F <= I-K-1) ||
						    (N==MDEBUG_N && MDEBUG_T<0 && MDEBUG_TT >= 0 && 0 <= MDEBUG_TS && MDEBUG_TS <= I-K-1)) && FLu1 >= MDEBUG_MIN){
		  FL128 sum = 0.0, sum2 = 0.0;
		  for(int i = 0; i <= max(MDEBUG_F,MDEBUG_TS); i++){
		    sum += LRfrag[i];
		    sum2 += LRfrag2[i];
		  }
		  printf("m=%d:I=%d,K=%d,J=%d,U=%d:FLu1=%0.10e,LRfrag[0..%d]=%0.10e-%0.10e=%0.10e\n",m,I,K,J,U,(double)FLu1,max(MDEBUG_F,MDEBUG_TS),(double)sum,(double)sum2,(double)(sum-sum2));
		  fflush(stdout);
		}
	      }

	      if(LPdelG && U <= I-K-2){
		if(U < I-K-2){
		  FL128 inc = FLu1 * InvFN;
		  if(DEBUG>=2) assert(isfinite(inc) && inc >= 0.0 && U+2 <= I-K);
		  /* LRdel[S] : S range for FLEu1 * FN is U+2 .. I-K-1 */
		  LRdel[U+2] += inc;
		  LRdel[I-K] -= inc;
		  if(DEBUG>=2) assert(isfinite(LRdel[U+2]));
		  if(DEBUG>=2) assert(isfinite(LRdel[I-K]));
		  if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && U+2 <= MDEBUG_S && MDEBUG_S <= I-K-1 && inc >= MDEBUG_MIN){
		    FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		    for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		      sum += LRdel[i];
		    printf("m=%d:I=%d,K=%d,J=%d,U=%d:FLu1*InvFN=%0.10e,LRdel[0..%d]=%0.10e\n",m,I,K,J,U,(double)inc,MDEBUG_S,(double)sum);
		    fflush(stdout);
		  }
		}

		FL128b FLu2 = (ERFC_FIX && U > Ub) ? FLIKd<1>(I,K,J,U,U+2,Yik,Xm[J],Yik-Y[U],Y[U+2]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1)
		  : FLIKd<0>(I,K,J,U,U+2,Yik,Xm[J],Yik-Y[U],Y[U+2]-Y[U],Y,IvarX, FLj, END_FIX ? lc : 1);
		if(!FLIK_FIX)
		  FLu2 *= FLj;

		if(DEBUG>=2) assert(isfinite(FLu2) && FLu2 >= 0.0);
		/* LRdel[S]: S range = U+1 */
		if(DEBUG>=2) assert(Dmin <= U+1 && U+2 <= Dmax);
		LRdel[U+1] += FLu2;
		LRdel[U+2] -= FLu2;
		if(DEBUG>=2) assert(isfinite(LRdel[U+1]));
		if(DEBUG>=2) assert(isfinite(LRdel[U+2]));

		if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && U+1 == MDEBUG_S && FLu2 >= MDEBUG_MIN){
		  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		    sum += LRdel[i];
		  printf("m=%d:I=%d,K=%d,J=%d,U=%d:FLu2=%0.10e,LRdel[0..%d]=%0.10e\n",m,I,K,J,U,(double)FLu2,MDEBUG_S,(double)sum);
		  fflush(stdout);
		}
	      } /* if(LPdelG ...) */

	      if((!EXTEND || J <= DELTA_XU || I-K <= DELTA_YU - FIXADD) && U >= ULmin){
		/* LRadd[S]: S range =  0 .. U-1 */
		if(U > 0){
		  LRadd[0] += FLu1;
		  if(DEBUG>=2) assert(Amin <= U && U <= Amax);
		  LRadd2[U] += FLu1;
		  if(MDEBUG_M >= 0 && MDEBUG_T >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && 0 <= MDEBUG_TS && MDEBUG_TS <= U-1 && FLu1 >= MDEBUG_MIN){
		    FL128 sum = LRadd[0],sum2 = LRadd2[0];
		    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
		      sum += LRadd[i];
		      sum2 += LRadd2[i];
		    }
		    FL128 sum3 = sum;
		    if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
		      sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
		      printf("m=%d:I=%d,K=%d,J=%d,U=%d(%d..%d,%d):FLu1=%0.10e,S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
			     m,I,K,J,U,Umin,Umax,Ub,(double)FLu1,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
		      //		      printf("\t FLj=%0.10e,IvarX=%0.10e,Yik=%0.8f,Xm[J]=%0.8f(AR=%0.10e,ScaleR=%0.10e)\n",(double)FLj,IvarX,Yik,Xm[J],A.AR(I,K,J),(double)ScaleR);
		    } else printf("m=%d:I=%d,K=%d,J=%d,U=%d(%d..%d,%d):FLu1=%0.10e,S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
				  m,I,K,J,U,Umin,Umax,Ub,(double)FLu1,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2));
		    fflush(stdout);
		  }
		}

		/* LRadd[S]: S range = U .. I-K-1 */
		FL128 inc = FLu1 * Fn;
		if(DEBUG>=2) assert((U==0 || Amin <= U) && Amin <= I-K && 1 <= I-K && U < I-K && I-K <= Amax);
		LRadd[U] += inc;
		LRadd2[I-K] += inc;
		if(MDEBUG_M >= 0 && MDEBUG_T >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && U <= MDEBUG_TS && MDEBUG_TS <= I-K-1 && inc >= MDEBUG_MIN){
		  FL128 sum = LRadd[0],sum2 = LRadd2[0];
		  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
		    sum += LRadd[i];
		    sum2 += LRadd2[i];
		  }
		  FL128 sum3 = sum;
		  if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
		    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
		    printf("m=%d:I=%d,K=%d,J=%d,U=%d(%d..%d,%d):FLu1*FN=%0.10e,S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
			   m,I,K,J,U,ULmin,Umax,Ub,(double)inc,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
		  } else printf("m=%d:I=%d,K=%d,J=%d,U=%d(%d..%d,%d):FLu1*FN=%0.10e,S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
				m,I,K,J,U,Umin,Umax,Ub,(double)inc,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2));
		  fflush(stdout);
		}
		if(0 && MDEBUG && LDEBUG && !FIXADD && m==1){
		  FL128 addsum = LRadd[0], addsum2 = LRadd2[0];
		  for(int i = 1; i <= 52;i++){
		    addsum += LRadd[i];
		    addsum2 += LRadd2[i];
		  }
		  addsum -= addsum2;
		  FL128 fragsum = LRfrag[0], fragsum2 = LRfrag2[0];;
		  for(int i = 1; i <= 52; i++){
		    fragsum += LRfrag[i];
		    fragsum2 += LRfrag2[i];
		  }
		  fragsum -= fragsum2;
		  if(DEBUG && (!(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190) || (VERB>=2&&m==1&&I==44&&K==0&&J==1))){
		    printf("m=%d:I=%d,K=%d,J=%d,U=%d:sum(LRadd[0..52])=%0.17e,sum(LRfrag[0..52]=%0.17e\n",
			   m,I,K,J,U,(double)addsum,(double)fragsum);
		    printf("    FLu1=%0.17e(0..U-1),FLu1*FN=%0.17e(U..I-K-1)\n",(double)FLu1,(double)inc);
		    for(int i=0;i <= 52; i++)
		      printf("        LRadd[%d]=%0.17e,LRfrag[%d]=%0.17e\n",i,(double)LRadd[i],i,(double)LRfrag[i]);
		    fflush(stdout);
		    assert(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190);
		  }
		}

		if(myaddcnt > 0){ /* LRaddT[S][] : S range is U */
		  double incf = inc;
		  for(int t = 0; t < addcntI[U]; t++){
		    double Yut = addlocI[U][t];
		    if(DEBUG>=2) assert(Y[U] < Yut && Yut < Y[U+1]);	      
		    if(DEBUG>=2) assert(0 <= Y[U+1] - Yut && Y[U+1] - Yut <= Yik - Yut);
		    if(DEBUG>=2) assert(0 <= Yut - Y[U] && Yut - Y[U] <= Yik - Y[U]);
		    /*		    double FLu3 = (ERFC_FIX && U > Ub) ? 
		      (FLIKA<1>(I,K,J,U,U+1,Yik,Xm[J],Y, Yut, Y[U+1], Yik - Yut, Y[U+1] - Yut, IvarX, FLj, 1) + Fn * FLIKA<1>(I,K,J,U,U+1,Yik,Xm[J], Y, Y[U], Yut, Yik - Y[U], Yut - Y[U], IvarX, FLj, END_FIX ? lc: 1))
		      : (FLIKA<0>(I,K,J,U,U+1,Yik,Xm[J],Y, Yut, Y[U+1], Yik - Yut, Y[U+1] - Yut, IvarX, FLj, 1) + Fn * FLIKA<0>(I,K,J,U,U+1,Yik,Xm[J], Y, Y[U], Yut, Yik - Y[U], Yut - Y[U], IvarX, FLj, END_FIX ? lc: 1));*/
		    double FLuA = (ERFC_FIX && U > Ub) ? 
		      FLIKA<1>(I,K,J,U,U+1,Yik,Xm[J],Y, Yut, Y[U+1], Yik - Yut, Y[U+1] - Yut, IvarX, FLj, 1) : 
		      FLIKA<0>(I,K,J,U,U+1,Yik,Xm[J],Y, Yut, Y[U+1], Yik - Yut, Y[U+1] - Yut, IvarX, FLj, 1);
		    double FLuB = (ERFC_FIX && U > Ub) ? 
		      FLIKA<1>(I,K,J,U,U+1,Yik,Xm[J], Y, Y[U], Yut, Yik - Y[U], Yut - Y[U], IvarX, FLj, END_FIX ? lc: 1) : 
		      FLIKA<0>(I,K,J,U,U+1,Yik,Xm[J], Y, Y[U], Yut, Yik - Y[U], Yut - Y[U], IvarX, FLj, END_FIX ? lc: 1);
		    double FLu3 = FLuA + FN[0] * FLuB;
		    if(!FLIK_FIX)
		      FLu3 *= FLj;
 		    if(DEBUG>=1+RELEASE /* HERE HERE >=2 */ && !(isfinite(FLu3) && FLu3 - inc >= -inc * (USE_MFLOAT ? 1e-4 : 2e-9) - (USE_MFLOAT ? 1e-16 : 1e-280))){
                      #pragma omp critical
		      {
			printf("m=%d:I=%d,K=%d,J=%d,U=%d(%d..%d,%d):FLu3=%0.10e,incf=%0.10e(FLu1=%0.10e,FLuA=%0.10e,FLuB=%0.10e,FN[0]=%0.10e)\n",
			       m,I,K,J,U,Umin,Umax,Ub,(double)FLu3,incf,(double)FLu1,FLuA,FLuB,FN[0]);
			fflush(stdout);
			//			assert(isfinite(FLu3) && FLu3 - inc >= -inc * (USE_MFLOAT ? 1e-4 : 2e-9) - (USE_MFLOAT ? 1e-16 : 1e-280));
			assert(isfinite(FLu3) && isfinite(FLu3 - incf));
		      }
		    }
		    if(DEBUG>=1+RELEASE) assert(Amin <= U && U <= min(N,Amax));		
		    LRaddT[U][t] += max((double)0.0, FLu3 - incf);
		    if(DEBUG>=1+RELEASE && !(isfinite(LRaddT[U][t]) && LRaddT[U][t] >= 0.0)){
		      #pragma omp critical
		      {
			printf("m=%d:I=%d,K=%d,J=%d,U=%d(%d..%d,%d):FLu3=%0.10e,incf=%0.10e(FLu1=%0.10e,Fn=%0.10e),LRaddT[U][t]= %0.8e\n",
			       m,I,K,J,U,Umin,Umax,Ub,(double)FLu3,incf,(double)FLu1,(double)Fn, LRaddT[U][t]);
			fflush(stdout);
			assert(isfinite(LRaddT[U][t]) && LRaddT[U][t] >= 0.0);
		      }
		    }
		    if(MDEBUG_M >= 0 && MDEBUG_T >= 0 && MDEBUG_TT == t && m==MDEBUG_M && N==MDEBUG_N && U == MDEBUG_TS && fabs((double)(FLu3-inc)) >= MDEBUG_MIN){
		      FL128 sum = LRadd[0],sum2 = LRadd2[0];
		      for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
			sum += LRadd[i];
			sum2 += LRadd2[i];
		      }
		      FL128 sum3 = sum;
		      sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
		      printf("m=%d:I=%d,K=%d,J=%d,U=%d(%d..%d,%d):FLu3-FLu1*FN=%0.10e,S=%d,T=%d(%0.3f):LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
			     m,I,K,J,U,Umin,Umax,Ub,(double)(FLu3-inc),U,t,addlocI[U][t],(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
		      fflush(stdout);
		    }
		  }
		}
	      }

	      if(lc && U == 0){ /* FLE/G term : U == 0 */
		FL128b FLu1 = FLEIK(I,K,J,0,1,Yik,Xm[J], FLj);
		if(!FLIK_FIX)
		  FLu1 *= FLj;
		if(DEBUG>=2) assert(isfinite(FLu1) && FLu1 >= 0.0);

		/* LRdel[S]: S range for FLEu1 is 1 .. max(1,U)-1 (empty) */

		if(LDEBUG && (!EXTEND || J <= DELTA_XU || I-K <= DELTA_YU) && U >= ULmin){
		  /* LRsite[S]: S range is 1..I-K-1 */
		  LRsite[1] += FLu1;
		  LRsite[I-K] -= FLu1;
		  if(LDEBUG && DEBUG>=2){/* verify that LRsite[1..N+1] adds up to approximately 0.0 */
		    FL128 sum = 0.0;
		    for(int S = 1; S <= N+1; S++)
		      sum += LRsite[S];
		    if((DEBUG && !(fabs((double)sum) <= LRcum * 1e-8)) || (0 &&m==MDEBUG_M && FLu1 >= MDEBUG_MIN && (1==STRACE || I-K==STRACE))){
		      printf("m=%d,I=%d,K=%d,J=%d,U=%d,N=%d,M=%d:FLu1=%0.10e,LRsite[1..N+1]=%0.10e,LRsite[%d]+=%0.10e,LRsite[%d]-=%0.10e\n",
			     m,I,K,J,U,N,M,(double)FLu1,(double)sum,1,(double)LRsite[1],I-K,(double)LRsite[I-K]);
		      fflush(stdout);
		      assert(fabs((double)sum) <= LRcum * 1e-8);
		    }
		  }

		  /* LRfrag[S]: S range is 0 .. I-K-1 */
		  LRfrag[0] += FLu1;
		  LRfrag2[I-K] += FLu1;
		}
		/* LRadd[S]: S range = 0 .. U-1 (empty) */

		if(LPdelG && I-K >= 2){
		  if(DEBUG>=2) assert(Dmin <= 1);
	      
		  if(I-K > 2){
		    FL128b inc = FLu1 * InvFN;
		    if(DEBUG>=2) assert(isfinite(inc) && inc >= 0.0);
		    /* LRdel[S] : S range for FLEu1 * FN is U+2 .. I-K-1 = 2 .. I-K-1 */
		    if(DEBUG>=2) assert(Dmin <= 2 && I-K <= Dmax);
		    LRdel[2] += inc;
		    LRdel[I-K] -= inc;
		    if(DEBUG>=2)assert(isfinite(LRdel[2]));
		    if(DEBUG>=2)assert(isfinite(LRdel[I-K]));
		  }

		  FL128b FLu2 = FLEIK(I,K,J,0,2,Yik,Xm[J], FLj);
		  if(!FLIK_FIX)
		    FLu2 *= FLj;
		  if(DEBUG>=2) assert(isfinite(FLu2) && FLu2 >= 0.0);
		  /* LRdel[S]: S range = U+1 = 1 */
		  if(DEBUG>=2) assert(Dmin <= 1 && 2 <= Dmax);
		  LRdel[1] += FLu2;
		  LRdel[2] -= FLu2;
		  if(DEBUG>=2) assert(isfinite(LRdel[1]));
		  if(DEBUG>=2) assert(isfinite(LRdel[2]));
		}/* if(LPdelG...) */

		if((!EXTEND || J <= DELTA_XU || I-K <= DELTA_YU - FIXADD) && U >= ULmin){
		  FL128b inc = FLu1 * Fn;
		  /* LRadd[S]: S range = U .. I-K-1 = 0 .. I-K-1 : include U from LRaddT[] since it does not depend on new site location */
		  if(DEBUG>=2) assert(Amin <= I-K && I-K <= Amax);
		  LRadd[0] += inc;
		  LRadd2[I-K] += inc;
		}
		/* LRaddT[S][] : S range is U : included in LRadd[S] above since it does not depend on new site location */
	      } /* lc && U==0 */
	    }/* U = Umin .. Umax */
	  }
	}
      }

      if(VERB>=2 && newLPdelta && rverb){
	printf("mprobeval:m=%d,Amin=%d,Amax=%d:After FL updates:time=%0.6f\n",m,Amin,Amax,wtime());
	fflush(stdout);
      }

      if(MDEBUG_M >= 0){
	if(MDEBUG_S >= 0 && m==MDEBUG_M && N==MDEBUG_N && LPdelG){
	  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
	    sum += LRdel[i];
	  printf("m=%d:After FL updates:LRdel[0..%d]=%0.10e(LS=%0.10f,LL=%0.10f)\n",m,MDEBUG_S,(double)sum,LS,LS+log((double)sum));
	  if(DEBUG && sum < 0.0){
	    printf("LRdel[0]=%0.10e\n",(double)LRdel[0]);
	    sum = LRdel[0];
	    for(int i = max(1,Dmin); i <= Dmax;i++){
	      sum += LRdel[i];
	      printf("LRdel[%d]=%0.10e,LRdel[0..%d]=%0.10e\n",i,(double)LRdel[i],i,(double)sum);
	    }
	    printf("MDEBUG_M=%d,MDEBUG_N=%d,MDEBUG_S=%d,MDEBUG_F=%d\n",MDEBUG_M,MDEBUG_N,MDEBUG_S,MDEBUG_F);
	    fflush(stdout);

	    assert(sum >= 0.0);
	  }
	}
	if(MDEBUG_F >= 0 && m==MDEBUG_M && N==MDEBUG_N-1){
	  FL128 sum = 0.0, sum2 = 0.0;
	  for(int i = 0; i <= MDEBUG_F; i++){
	    sum += LRfrag[i];
	    sum2 += LRfrag2[i];
	  }
	  printf("m=%d:After FL updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f,LL=%0.10f)\n",m,MDEBUG_F,(double)sum,(double)sum2,(double)(sum-sum2),LS,LS+log((double)(sum-sum2)));
	}
	if(MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  FL128 sum = LRadd[0], sum2 = LRadd2[0];
	  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
	    sum += LRadd[i];
	    sum2 += LRadd2[i];
	  }
	  FL128 sum3 = sum;
	  if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
	    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After FL updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e,LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2),LS);
	  } else
	    printf("m=%d:After FL updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e(LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),LS);
	  if(LDEBUG){
	    sum = sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After FL updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f,LL=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS,LS+log((double)(sum-sum2)));
	  }
	}
	if((MDEBUG_ST >= 0 && m==MDEBUG_M) || (MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N+1)){
	  FL128 sum = 0.0;
	  for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++)
	    sum += LRsite[i]; 
	  printf("m=%d:After FL updates:LRsite[0..%d]=%0.10e(LS=%0.10f)\n",m,max(MDEBUG_TS+1,MDEBUG_ST),(double)sum,LS);
	}
	if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  if(newLPdelta && MDEBUG_TT < D[min(N,MDEBUG_TS)]){
	    FL128 sum = LRdelta[0];
	    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
	      sum += LRdelta[i];
	    FL128 sum3 = sum;
	    if(Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < D[MDEBUG_TS])
	      sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After FL updates:S=%d,D=%d,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e,LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3,LS);
	  }
	  if(LDEBUG){
	    FL128 sum = 0.0,sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After FL updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f,LL=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS,LS+log((double)(sum-sum2)));
	  }
	}
	fflush(stdout);
      }
      if(0 && MDEBUG && LDEBUG && !FIXADD){
	FL128 addsum = LRadd[0],addsum2 = LRadd2[0];
	for(int i = 1; i <= 52;i++){
	  addsum += LRadd[i];
	  addsum2 += LRadd2[i];
	}
	addsum -= addsum2;
	FL128 fragsum = LRfrag[0], fragsum2 = LRfrag2[0];
	for(int i = 1; i <= 52; i++){
	  fragsum += LRfrag[i];
	  fragsum2 += LRfrag2[i];
	}
	fragsum -= fragsum2;
	if(DEBUG && (!(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190) || (VERB>=2&&m==1))){
	  printf("m=%d:After FL updates:sum(LRadd[0..52])=%0.17e,sum(LRfrag[0..52]=%0.17e\n", m,(double)addsum,(double)fragsum);
	  fflush(stdout);
	  assert(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190);
	}
      }
      if(LDEBUG){/* check that cumulative value of LRsite[1..N+1] approximates 0 */
	FL128 LRsitecum = 0.0;
	for(int S = 1; S <= N+1; S++){
	  LRsitecum += LRsite[S];
	  if(DEBUG && (VERB>=2 || (DEBUG && !(S <= N || fabs((double)LRsitecum) <= LRcum * 1e-8)))){
	    #pragma omp critical
	    {
	      printf("m=%d:S=%d,N=%d:LRsite[S]=%0.10e,LRsitecum=%0.10e,LRcum=%0.10e,LRfrag[0]=%0.10e,LRsitecum-LRcum=%0.10e,IMIN=%d,IMAX=%d\n",
		     m,S,N,(double)LRsite[S],(double)LRsitecum,LRcum,(double)LRfrag[0],(double)(LRsitecum-LRcum),IMIN,IMAX);
	      FL128 LRsitecum2 = 0.0;
	      for(int s = 1; s <= N+1; s++){
		LRsitecum2 += LRsite[s];
		printf("LRsite[%d]= %0.10e, LRfrag[%d]=%0.10e, LRsitecum= %0.10e\n",s,(double)LRsite[s],s,(double)LRfrag[s],(double)LRsitecum2);
	      }
	      fflush(stdout);
	      if(DEBUG && S == N+1) assert(fabs((double)LRsitecum) <= LRcum * 1e-8);
	    }
	  }
	}
      }

      if(DEBUG>=2 && newLPdelta){
	if(VERB>=3){
	  if(m==7 && Amin <= 12 && 12 <= Amax && D[12] >= 1 /* && !isfinite(LRdeltaD[12][0])*/)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[12]=%d,delta[12][0]= %0.4f,LRdeltaD[12][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[12],delta[12][0],LRdeltaD[12][0]);
	  if(m==56 && Amin <= 31 && 31  <= Amax && D[31] >= 1 /* && !isfinite(LRdeltaD[31][0] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[31]=%d,delta[31][0]= %0.4f,LRdeltaD[31][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[31],delta[31][0],LRdeltaD[31][0]);
	  if(m==20 && Amin <= 33 && 33 <= Amax && D[33] >= 11 /* && !isfinite(LRdeltaD[33][10] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[33]=%d,delta[33][10]= %0.4f,LRdeltaD[33][10]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[33],delta[33][10],LRdeltaD[33][10]);
	  fflush(stdout);
	}
	for(int S = Amin;S <= min(N,Amax); S++)
	  for(int d = 0; d < D[S]; d++)
	    assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
      }

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After Updates based on FL term:wall time= %0.6f\n", m, tid, wtime() - wt);
	  fflush(stdout);
        }
      }

      /* Updates based on FA term (Line 3 of Equation 7 for LRdel and Equation 10 for LRadd and Equation ?? for LRdelta) */
      long long Jcnt,Kcnt;
      if(VERB>=2)
	Jcnt = Kcnt = 0;

      for(int I = IMIN; I <= IMAX; I++){
	for(int K = 0; K <= Kmax[I]; K++){
	  double Yik = Yc(Y,I,K);
	  if(MDEBUG_M >= 0 && m==MDEBUG_M && ((MDEBUG_S >= 0 && N==MDEBUG_N && I+1 <= MDEBUG_S) ||
					      (MDEBUG_F >= 0 && N==MDEBUG_N-1 && I+1 <= MDEBUG_S)/*  ||
					      (MDEBUG_TS >=0 && MDEBUG_TT >= 0 && I <= MDEBUG_TS && N==MDEBUG_N+1) ||
					      (MDEBUG_T >= 0 && MDEBUG_TT >= 0 && I <= MDEBUG_TS && N==MDEBUG_N) */) && LPdelG){
	    FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	    for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
	      sum += LRdel[i];
	    printf("m=%d,I=%d,K=%d,J=%d..%d,N=%d,M=%d,LRdel[0..%d]=%0.10e\n",m,I,K,Jmin[I],Jmax[I], N, M, MDEBUG_S,(double)sum);
	    if(VERB>=2 && sum > 0.0){
	      sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	      double tmp = sum;
	      printf("    LRdel[0]=%0.17e+%0.17e\n",tmp,(double)(LRdel[0]-tmp));
	      for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax);i++){
		sum += LRdel[i];
		tmp = sum;
		printf("    LRdel[%d]=%0.17e (cum=%0.17e+%0.17e)\n",i,(double)LRdel[i],tmp,(double)(sum-tmp));
	      }
	    }
	    fflush(stdout);
	  }
	  for(int J = Jmin[I]; J <= Jmax[I]; J++){
	    double AL = A.AL(I,K,J);
	    if(DEBUG>=2) assert(A.AL(I,K,J) >= 0.0);
	    int Qmin = J+1;
	    int Qmax = min(M, J + deltaExtXRef);
	    int Qmax2 = min(M, J + DELTA_X);

	    int Pmin = I+1;
	    int Pmax = min(IMAX, I + deltaExtYRef);
	    int PmaxD = min(IMAX, I + deltaExtYRef + 1);/* for site deletion */
	    int PmaxA = min(IMAX, I + deltaExtYRef - FIXADD);/* for site addition */
	    int Pmax2 = min(IMAX, I + DELTA_Y);
	    int PmaxD2 = min(IMAX, I + DELTA_Y + 1);/* for site deletion */
	    int PmaxA2 = min(IMAX, I + DELTA_Y - FIXADD);/* for site addition */
	    int QmaxD = Qmax;/* for site deletion : only differs from Qmax if PQ_BND > 0 */
	    int QmaxD2 = Qmax2;/* for site deletion : only differs from Qmax2 if PQ_BNX > 0 */

#if PQ_BND  
	    if(PQ_BND){// narrow range P = Pmin .. Pmax* and Q = max(Qmin,Jmin[P]) .. min(Qmax*,Jmax[P]) based on EPS3 (NOTE : NOT applied for site deletion)
	      int Pmin3 = A.Pmin(I,K,J);
	      int Pmax3 = A.Pmax(I,K,J);
	      int Qmin3 = A.Qmin(I,K,J);
	      int Qmax3 = A.Qmax(I,K,J);
	      if(VERB>=3 || (MDEBUG_M >= 0 && m==MDEBUG_M && ((MDEBUG_S >= 0 && N==MDEBUG_N && I+1 <= MDEBUG_S && MDEBUG_S <= PmaxD-1) ||
							      (MDEBUG_F >= 0 && N==MDEBUG_N-1 && I+1 <= MDEBUG_S && MDEBUG_S <= Pmax-1)) && LPdelG)){
		printf("m=%d,I=%d,K=%d,J=%d:A.Pmin=%d,A.Pmax=%d,A.Qmin=%d,A.Qmax=%d:Pmin=%d,Pmax=%d,PmaxD=%d,PmaxA=%d,Pmax2=%d,PmaxD2=%d,PmaxA2=%d,Qmin=%d,Qmax=%d,Qmax2=%d\n",
		       m,I,K,J,Pmin3,Pmax3,Qmin3,Qmax3,Pmin,Pmax,PmaxD,PmaxA,Pmax2,PmaxD2,PmaxA2,Qmin,Qmax,Qmax2);
		fflush(stdout);
	      }

	      if(DEBUG>=2 && EPS3 <= 0.0){
		assert(Pmin3 == I+1);
		assert(Pmax3 == I + deltaExtYRef);
		assert(Qmin3 == J+1);
		assert(Qmax3 == J + deltaExtXRef);

		assert(Pmin3 <= Pmin);
		assert(Pmax3 >= Pmax);
		assert(Pmax3 + 1 >= PmaxD);
		assert(Pmax3 - FIXADD >= PmaxA);
		assert(Pmax3 >= Pmax2);
		assert(Pmax3 + 1 >= PmaxD2);
		assert(Pmax3 - FIXADD >= PmaxA2);

		assert(Qmin3 <= Qmin);
		assert(Qmax3 >= Qmax);
		assert(Qmax3 >= Qmax2);
	      }

	      Pmin = max(Pmin,Pmin3);
	      Pmax = min(Pmax,Pmax3);
	      PmaxD = min(PmaxD,Pmax3 + 1);// HERE : should this remain unchanged since PQ_BNX did not see upto terms at PmaxD and PmaxD2 ?
	      PmaxA = min(PmaxA,Pmax3/* WAS - FIXADD */);
	      Pmax2 = min(Pmax2,Pmax3);
	      PmaxD2 = min(PmaxD2,Pmax3 + 1);// HERE : should this remain unchanged since PQ_BNX did not see upto terms at PmaxD and PmaxD2 ?
	      PmaxA2 = min(PmaxA2,Pmax3/* WAS - FIXADD */);

	      Qmin = max(Qmin,Qmin3);
	      Qmax = min(Qmax,Qmax3);
	      Qmax2 = min(Qmax2,Qmax3);
	    }
#endif
	    if(MDEBUG_M >= 0 && m==MDEBUG_M && ((MDEBUG_S >= 0 && N==MDEBUG_N && I+1 <= MDEBUG_S && MDEBUG_S <= PmaxD-1) ||
						(MDEBUG_F >= 0 && N==MDEBUG_N-1 && I+1 <= MDEBUG_S && MDEBUG_S <= Pmax-1) /* ||
						(MDEBUG_TS >=0 && MDEBUG_TT >= 0 && I <= MDEBUG_TS && MDEBUG_TS+1 <= Pmax-1 && N==MDEBUG_N+1) ||
						(MDEBUG_T >= 0 && MDEBUG_TT >= 0 && I <= MDEBUG_TS && MDEBUG_TS <= PmaxA-1 && N==MDEBUG_N) */) && LPdelG){
	      FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	      for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		sum += LRdel[i];
	      printf("m=%d,I=%d,K=%d,J=%d,Pmax=%d,%d,Qmax=%d,%d,PmaxD=%d,%d,QmaxD=%d,%d,PmaxA=%d,%d,N=%d,M=%d,LRdel[0..%d]=%0.10e\n",
		     m,I,K,J,Pmax,Pmax2,Qmax,Qmax2,PmaxD,PmaxD2,QmaxD,QmaxD2,PmaxA,PmaxA2,N, M, MDEBUG_S, (double)sum);
	      fflush(stdout);
	    }

	    for(int P = Pmin; P <= IMAX && P <= PmaxD + A.kmax /* WAS KMAX */;P++){
	      int QminP = max(Qmin, Jmin[P]);
	      int Tmin = max(0,P-PmaxD);
	      int Tmax = min(Kmax[P],P-I-1);
	      if(Tmax < Tmin)
		continue;

	      if(DEBUG>=2) assert(Tmax <= A.kmax && P <= PmaxD + A.kmax && P <= PmaxD + Kmax[P]);	      

	      double ALscale = AL * exp(LogScaleL[I] + LogScaleR[P] - LS);
	      if(DEBUG>=2) assert(isfinite(ALscale));

	      /* expand LRdel[Dmin .. Dmax] to include LRdel[I+1 .. P - Tmin] */
	      if(DEBUG>=2 && LPdelG) assert(Dmin <= I+1);
	      if(DEBUG>=2 && LPdelG) assert(Dmax >= P-Tmin);	    
	      /* expand LRadd[0,Amin .. Amax] to include LRadd[I .. PmaxA] (also LRdelta[]) */
	      if(DEBUG>=2) assert(Amin <= I);
	      if(DEBUG>=2) assert(Amax >= PmaxA);

	      for(int T = Tmin; T <= Tmax; T++){
		if(DEBUG>=2) assert(P-T >= I+1 && P-T <= PmaxD);

		MFLOAT Ivar;
		MFLOAT deltaY = Yc(Y,P,T) - Yik;
		MFLOAT FAik = FAIK(deltaY,PRtab[P-T-I-1][P-T].Pr,P-T-I,P,T,K,N,VARtab[T*N+P],VARtab[K*N+I],Ivar,Y);
		MFLOAT LRik = FAik * InvFN;/* term modified for site deletion */
		int QmaxPD = (P-T <= PmaxD2) ? min(QmaxD, Jmax[P]) : min(QmaxD2, Jmax[P]);
		int QmaxPA = (P-T <= PmaxA2) ? min(Qmax, Jmax[P]) : min(Qmax2, Jmax[P]);
		if(DEBUG>=2) assert(QmaxPA <= QmaxPD);
		MFLOAT2 *ARpt = &A.AR(P,T,0);
		FL128 LRipt128 = 0.0, LRiptD128 = 0.0, LRiptA128 = 0.0;
		MFLOAT2 LRipt = 0.0;
		int Q = QminP;

		if(CRITLOOP_MERGE && OutlierType){ // Critical Loop for vectorization : used by both LPdel[] and LPadd[],LPadd2[] (if OutlierType==1) */

		  #pragma omp simd reduction(+:LRipt)
		  for(int q = Q; q <= QmaxPA; q++){
#if DEBUG>=2
		    if(DEBUG>=2)assert(A.AR(P,T,q) >= 0.0);
#endif
		    LRipt += FAJ(Xm[q]-Xm[J], deltaY, q-J, P-T-I, FBiasWTm[q-J], LogFpPow[q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], LRik, Ivar, 1, OUT_TYPE, OUT_OLD) * ARpt[q];
		  }

		  Q = max(Q, QmaxPA + 1);
		  LRiptA128 = LRipt * ALscale * FN[0] * FN[0];
		  if(DEBUG>=2) assert(isfinite(LRiptA128) && LRiptA128 >= 0.0);
		  if(DEBUG>=2 && max(J+1, Jmin[P]) > QmaxPA) assert(LRipt == 0.0 && LRiptA128 == 0.0);
		}

		if(LPdelG){/* handle site deletion */
		  if(MDEBUG_M >= 0 && m==MDEBUG_M && ((MDEBUG_S >= 0 && N==MDEBUG_N && I+1 <= MDEBUG_S && MDEBUG_S <= P-T-1) ||
						      (MDEBUG_F >= 0 && N==MDEBUG_N-1 && I+1 <= MDEBUG_S && MDEBUG_S <= P-T-1) /* ||
						      (MDEBUG_TS>= 0 && MDEBUG_TT >= 0 && I+1 <= MDEBUG_TS+1 && MDEBUG_TS+1 <= P-T-1 && N==MDEBUG_N+1) ||
						      (MDEBUG_T >= 0 && MDEBUG_TT >= 0 && N==MDEBUG_N && I <= MDEBUG_TS && MDEBUG_TS <= P-T-1 && N==MDEBUG_N)*/)){
		    FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		    for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		      sum += LRdel[i];
		    printf("m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d(Jmin[P]=%d,Jmax[P]=%d,QmaxPA=%d),LRdel[0..%d]=%0.10Le\n",m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxPD,Jmin[P],Jmax[P],QmaxPA,MDEBUG_S, (long double)sum);
		    fflush(stdout);
		  }
		
		  //		  #pragma omp simd reduction(+:LRipt)
		  for(int q = Q; q <= QmaxPD; q++){// If OutlierType==1 : short extension of critical loop above
#if DEBUG>=2
		    if(DEBUG>=2)assert(A.AR(P,T,q) >= 0.0);
#endif
		    LRipt += FAJ(Xm[q]-Xm[J],deltaY,q-J,P-T-I,FBiasWTm[q-J], LogFpPow[q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I],LRik, Ivar,1, OUT_TYPE, OUT_OLD) * ARpt[q];
		  }
		  // Q = max(Q, QmaxPD + 1);

		  LRiptD128 = LRipt * ALscale;
		  if(DEBUG>=2) assert(isfinite(LRiptD128) && LRiptD128 >= 0.0);
		  if(DEBUG>=2 && max(J+1,Jmin[P]) > QmaxPD && !(LRipt == 0.0 && LRiptD128 == 0.0)){
		    #pragma omp critical
		    {
		      printf("m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d(Jmin[P]=%d,Jmax[P]=%d,QmaxPA=%d),LRipt=%0.24Le,LRiptD128=%0.24Le\n",
			     m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxPD,Jmin[P],Jmax[P],QmaxPA,(long double)LRipt, (long double)LRiptD128);
		      fflush(stdout);
		      assert(LRipt == 0.0 && LRiptD128 == 0.0);
		    }
		  }
		  
		  /* LRdel[S]: S range is I+1 .. P-T-1 */
		  if(DEBUG>=2) assert(Dmin <= I+1 && I+1 <= P-T && P-T <= Dmax);
		  LRdel[I+1] += LRiptD128;
		  LRdel[P-T] -= LRiptD128;
		  if(DEBUG>=2) assert(isfinite(LRdel[I+1]));
		  if(DEBUG>=2) assert(isfinite(LRdel[P-T]));

		  if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && I+1 <= MDEBUG_S && MDEBUG_S <= P-T-1 && LRiptD128 >= MDEBUG_MIN){
		    FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		    for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		      sum += LRdel[i];
		    printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d:AL(FA/FN)*AR=%0.10e,LRdel[0..%d]=%0.10e:Yc(I,K)=%0.3f,Yc(P,T)=%0.3f,X[J]=%0.3f,X[%d]=%0.3f,FAik=%0.10e,FAik/FN=%0.10e\n",
			   m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxPD,(double)LRiptD128,MDEBUG_S,(double)sum,Yc(Y,I,K),Yc(Y,P,T),Xm[J],QmaxPD,Xm[QmaxPD],FAik,LRik);
		    FL128b LRipt2 = 0.0;
		    for(int Q = QminP; Q <= QmaxPD; Q++){
		      //		      int verb = (MDEBUG && m==MDEBUG_M && N==MDEBUG_N && I==2 && K==0 && J==7 && P==4 && T==0 && 8 <= Q && Q <= 8) ? 1 : 0;
		      double FAijkd = FAJ(Xm[Q]-Xm[J],deltaY,Q-J,P-T-I,FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], LRik,Ivar,1, OUT_TYPE, OUT_OLD/*, verb*/);
		      //		        rverb = 0;

		      LRipt2 += FAijkd * A.AR(P,T,Q);
		      printf("    Q=%d:AL*Scale=%0.10e,FAijk/FN=%0.10e(x=%0.4f,y=%0.4f),AR=%0.10e, Scale=%0.10e: sum(AL(FA/FN)*AR)=%0.10e(delta=%0.10e)\n",
			     Q, (double)ALscale,FAijkd,Xm[Q]-Xm[J],deltaY,A.AR(P,T,Q), exp(LogScaleL[I] + LogScaleR[P] -LS), (double)(LRipt2 * ALscale), (double)(ALscale * FAijkd * A.AR(P,T,Q)));
		    }
		    if(0 && MDEBUG && m==MDEBUG_M && N==MDEBUG_N && I==30 && K==0 && J==1 && P==33 && T==1){
		      printf("LogScaleL[I]= %0.10f, LogScaleR[P]= %0.10f, LS= %0.10f, FAik=%0.10e, LRik= %0.10e, Ivar=%0.10e, FN= %0.10e, LogFnPow[P-T-I-1]=%0.10e\n", 
			     LogScaleL[I], LogScaleR[P], LS, FAik, LRik, Ivar,FN[0], LogFnPow[P-T-I-1]);
		      for(int Q = QminP; Q <= QmaxPD; Q++)
			printf(" Q=%d:FBiasWTm[Q-J]=%0.10e,LogFpPow[Q-J-1]=%0.10e\n",Q,FBiasWTm[Q-J],LogFpPow[Q-J-1]);
		    }
		    fflush(stdout);
		  }
		} /* if(LPdelG) */

		/* remaining terms below do not involve any site deletion */
		if(P-T > Pmax)
		  continue;
		int QmaxP = (P - T <= Pmax2) ? min(Qmax, Jmax[P]) : min(Qmax2, Jmax[P]);
		if(DEBUG>=2) assert(QmaxPA <= QmaxP && QmaxP <= QmaxPD);
		if(DEBUG>=2 && !FIXADD) assert(QmaxPA == QmaxP);
		
		if(LDEBUG){	      /* use un-modified term */
		  FL128b LRipt = 0.0;
		  for(int Q = QminP; Q <= QmaxP; Q++){
		    if(DEBUG>=2)assert(A.AR(P,T,Q) >= 0.0);
		    LRipt += FAJ(Xm[Q]-Xm[J],deltaY,Q-J,P-T-I,FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], FAik,Ivar,1, OUT_TYPE, OUT_OLD) * A.AR(P,T,Q);
		  }
		  LRipt128 = LRipt * ALscale;
		  if(DEBUG>=2) assert(isfinite(LRipt128) && LRipt128 >= 0.0);

		  /* LRsite[S] : S range is I+1 .. P-T-1 */
		  LRsite[I+1] += LRipt128;
		  LRsite[P-T] -= LRipt128;
		  if(LDEBUG && DEBUG>=2){/* verify that LRsite[1..N+1] adds up to approximately 0.0 */
		    FL128 sum = 0.0;
		    for(int S = 1; S <= N+1; S++)
		      sum += LRsite[S];
		    if((DEBUG && !(fabs((double)sum) <= LRcum * 1e-8)) ||(0 &&m==MDEBUG_M && LRipt128 >= MDEBUG_MIN && (I+1==STRACE || P-T==STRACE))){
		      printf("m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d,N=%d,M=%d:sum(AL*FA*AR)=%0.10e,LRsite[1..N+1]=%0.10e,LRsite[%d]+=%0.10e,LRsite[%d]-=%0.10e\n",
			     m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxP,N,M,(double)LRipt128,(double)sum,I+1,(double)LRsite[I+1],P-T,(double)LRsite[P-T]);
		      fflush(stdout);
		      assert(fabs((double)sum) <= LRcum * 1e-8);
		    }
                  }
		  if((MDEBUG_M < 0 || m==MDEBUG_M) && ((MDEBUG_ST >= 0 && I+1 <= MDEBUG_ST && MDEBUG_ST <= P-T-1)||
						       (MDEBUG_TS>=0 && MDEBUG_TT>=0 && I+1 <= MDEBUG_TS+1 && MDEBUG_TS+1 <= P-T-1 && m==MDEBUG_M && N==MDEBUG_N+1)) && LRipt128 >= MDEBUG_MIN){
		    FL128 sum = 0.0;
		    for(int i = 1; i <= max(MDEBUG_ST,MDEBUG_TS+1); i++)
		      sum += LRsite[i];
		    int t = max(MDEBUG_ST,MDEBUG_TS+1);
		    printf("m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d,N=%d,M=%d:sum(AL*FA*AR)=%0.10e,LRsite[1..%d]=%0.10e(LRsite[%d]=%0.10e,LRsite[%d]=%0.10e)\n",
			   m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxP,N,M,(double)LRipt128,t,(double)sum,t,(double)LRsite[t],t+1,(double)LRsite[t+1]);
		    for(int Q = QminP/* WAS max(J+1,Jmin[P]) */; Q <= QmaxP; Q++)
		      printf("    Q=%d:AL*FA*AR=%0.10e (AL=%0.10e,FA=%0.10e,AR=%0.10e):Y[I]=%0.3f,Y[P-T]=%0.3f,Xm[J]=%0.3f,Xm[Q]=%0.3f\n",
			     Q, (double)(ALscale) * FAJ(Xm[Q]-Xm[J],deltaY,Q-J,P-T-I,FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], FAik,Ivar,1, OUT_TYPE, OUT_OLD) * A.AR(P,T,Q), A.AL(I,K,J) * exp(LogScaleL[I]), 
			     FAJ(Xm[Q]-Xm[J],deltaY,Q-J,P-T-I,FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], FAik,Ivar,1, OUT_TYPE, OUT_OLD)*exp(-LS), A.AR(P,T,Q) * exp(LogScaleR[P]), Y[I],Y[P-T],Xm[J],Xm[Q]);
		    fflush(stdout);
		  }

		  /* LRfrag[S] : S range is I .. P-T-1 */
		  LRfrag[I] += LRipt128;
		  LRfrag2[P-T] += LRipt128;
		  if((MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N-1 && I <= MDEBUG_F && MDEBUG_F <= P-T-1 && LRipt128 >= MDEBUG_MIN) ||
		     (MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && MDEBUG_TT >= 0 && I <= MDEBUG_TS && MDEBUG_TS <= P-T-1 && LRipt128 >= MDEBUG_MIN /* && AL >= MDEBUG_MIN*/)){
		    FL128 sum = 0.0, sum2 = 0.0;
		    for(int i = 0; i <= max(MDEBUG_F,MDEBUG_TS); i++){
		      sum += LRfrag[i];
		      sum2 += LRfrag2[i];
		    }
		    printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d:deltaY=%0.3f,FAik=%0.10e,sum(AL*FA*AR)=%0.10e:Yc(I,K)=%0.3f,Yc(P,T)=%0.3f,X[J]=%0.3f,X[%d]=%0.3f,LRfrag[0..%d]=%0.10e-%0.10e=%0.10e\n",
			   m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxP,deltaY,FAik,(double)LRipt128,Yc(Y,I,K),Yc(Y,P,T),Xm[J],QmaxP,Xm[QmaxP],max(MDEBUG_F,MDEBUG_TS),(double)sum,(double)sum2,(double)(sum-sum2));
		    double tsum = 0.0;
		    for(int Q = QminP /* WAS max(J+1,Jmin[P]) */; Q <= QmaxP; Q++){
		      //		      int verb = (MDEBUG && m==MDEBUG_M && N==MDEBUG_N - 1 && I==2 && K==0 && J==7 && P==3 && T==0 && 8 <= Q && Q <= 8) ? 1 : 0;
		      double FAikj = FAJ(Xm[Q]-Xm[J],deltaY,Q-J,P-T-I,FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], FAik,Ivar,1, OUT_TYPE, OUT_OLD/*, verb */);
		      //		        rverb = 0;

		      tsum += FAikj * A.AR(P,T,Q);
		      printf("    Q=%d:AL*Scale=%0.10e,FAikj=%0.10e(x=%0.4f,y=%0.4f),AR=%0.10e,Scale=%0.10e:sum(AL*FA*AR)=%0.10e (delta=%0.10e)\n",
			     Q,(double)ALscale,FAikj, Xm[Q]-Xm[J],deltaY,A.AR(P,T,Q),exp(LogScaleL[I]+LogScaleR[P]-LS), (double)(tsum*ALscale),(double)(ALscale * FAikj * A.AR(P,T,Q)));
		    }
		    if(0 && MDEBUG && m==MDEBUG_M && N==MDEBUG_N - 1 && I==30 && K==0 && J==1 && P==32 && T==1){
		      printf("LogScaleL[I]= %0.10f, LogScaleR[P]= %0.10f, LS= %0.10f,FAik=%0.10e,Ivar=%0.10e,FN= %0.10e,LogFnPow[P-T-I-1]=%0.10e\n", LogScaleL[I], LogScaleR[P], LS,FAik,Ivar,FN[0],LogFnPow[P-T-I-1]);
		      for(int Q = QminP; Q <= QmaxP; Q++)
			printf(" Q=%d:FBiasWTm[Q-J]=%0.10e,LogFpPow[Q-J-1]=%0.10e\n",Q,FBiasWTm[Q-J],LogFpPow[Q-J-1]);
		    }
		  }
		}

		// LPdelta moved into its own loop below

		/* remaining terms below involve site addition only */
		if(P-T > PmaxA)
		  continue;
	      
		MFLOAT LRikA = FAik * FN[0];/* term modified for site addition */
		if(DEBUG>=3 || !(CRITLOOP_MERGE && OutlierType)){// NOTE : following loop is a subrange (QmaxPA <= QmaxPD) of the LRipt summation critical loop above
		  MFLOAT2 LRiptA = 0.0;
		  for(int Q = QminP /* WAS max(J+1, Jmin[P]) */; Q <= QmaxPA; Q++){
		    if(DEBUG>=2)assert(A.AR(P,T,Q) >= 0.0);
		    LRiptA += FAJ(Xm[Q]-Xm[J],deltaY,Q-J,P-T-I,FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], LRikA, Ivar,1, OUT_TYPE, OUT_OLD) * ARpt[Q];
		  }
		  if(DEBUG>=3 && CRITLOOP_MERGE && OutlierType && !(fabsl((long double)LRiptA128 - (LRiptA * ALscale)) <= fabsl((long double)LRiptA128) * 1e-8)){
		    #pragma omp critical
		    {
		      printf("m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d:LRiptA128=%0.24Le,LRiptA*ALscale=%0.24Le(err=%0.24Le)\n",
			     m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxPA,(long double)LRiptA128,(long double)(LRiptA * ALscale), (long double)(LRiptA128 - LRiptA*ALscale));
		      fflush(stdout);
		      assert(fabsl((long double)LRiptA128 - (LRiptA * ALscale)) <= fabsl((long double)LRiptA128) * 1e-8);
		    }
		  }
	          LRiptA128 = LRiptA * ALscale;
		  if(DEBUG>=2) assert(isfinite(LRiptA128) && LRiptA128 >= 0.0);
		} // 2nd copy of crit. loop : should not be needed unless OutlierType == 0

		if(LDEBUG && DEBUG>=2 && !((QmaxPA < QmaxP || LRiptA128 >= LRipt128*Fn*(1.0 - (USE_MFLOAT ? 1e-4 : 1e-6)) - max(1.0,ALscale)*(USE_MFLOAT ? 1e-20 : 1e-190)) 
		    && LRiptA128 <= LRipt128*Fn*(1.0 + (USE_MFLOAT ? 1e-4 : 1e-6)) + max(1.0,ALscale)*(USE_MFLOAT ? 1e-20 : 1e-190))){
		  printf("m=%d,I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d:FAik*FN=%0.10e,FAik=%0.10e,LRiptD128=%0.10e,LRipt128=%0.10e,LRiptA128=%0.10e,LRipt128*FN=%0.10e,ALscale=%0.10e,LSL[I]=%0.6f,LSR[P]=%0.6f,LS=%0.6f\n",
			 m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxP,LRikA,FAik,(double)LRiptD128,(double)LRipt128,(double)LRiptA128,(double)(LRipt128*Fn),(double)ALscale, LogScaleL[I],LogScaleR[P],LS);
		  FL128b LRipt2 = 0.0, LRiptA2 = 0.0;
		  for(int Q = QminP /* WAS max(J+1,Jmin[P]) */; Q <= QmaxP; Q++){
		    double FAikj = FAJ(Xm[Q]-Xm[J],deltaY,Q-J,P-T-I,FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], FAik,Ivar,1, OUT_TYPE, OUT_OLD);
		    LRipt2 += FAikj * A.AR(P,T,Q);
		    double FAikjA = FAJ(Xm[Q]-Xm[J],deltaY,Q-J,P-T-I,FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], LRikA,Ivar,1, OUT_TYPE, OUT_OLD);
		    LRiptA2 += FAikjA * A.AR(P,T,Q);
		    printf("  Q=%d:FAikj=%0.10e,FAikjA=%0.10e,FAikj*FN=%0.10e, AL*sum(FAikj*AR)=%0.10e,AL*sum(FAikjA*AR)=%0.10e,AL*sum(FAikj*AR)*FN=%0.10e,AR=%0.10e\n",
			   Q, (double)FAikj,(double)FAikjA,(double)(FAikj*Fn), (double)(ALscale*LRipt2), (double)(ALscale*LRiptA2), (double)(ALscale*LRipt2*Fn), A.AR(P,T,Q));
		  }
		  fflush(stdout);
		  if(QmaxPA == QmaxP)
		    assert(LRiptA128 >= LRipt128*Fn*(1.0 - (USE_MFLOAT ? 1e-4 : 1e-6)) - max(1.0,ALscale)*(USE_MFLOAT ? 1e-20 : 1e-190));
		  assert(LRiptA128 <= LRipt128*Fn*(1.0 + (USE_MFLOAT ? 1e-4 : 1e-6)) + max(1.0,ALscale)*(USE_MFLOAT ? 1e-20 : 1e-190));
		  /* NOTE for USE_MFLOAT=0 : using 1e-190 since AR can be as large as 1e+100, hence with FAikj = 1e-307, underflow may have occured even with values up to 1e-190*ALscale */
		}

		/* LRadd[S]: S range is I .. P-T-1 */
		if(DEBUG>=2) assert(1 <= I && Amin <= I && I < P-T && P-T <= Amax);
		LRadd[I] += LRiptA128;
		LRadd2[P-T] += LRiptA128;
		if(MDEBUG_M >= 0 && MDEBUG_T >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && I <= MDEBUG_TS && MDEBUG_TS <= P-T-1 && LRiptA128 >= MDEBUG_MIN){
		  FL128 sum = LRadd[0],sum2 = LRadd2[0];
		  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
		    sum += LRadd[i];
		    sum2 += LRadd2[i];
		  }
		  FL128 sum3 = sum;
		  if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
		    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
		    printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d:FAik*FN=%0.10e,LRiptA=%0.10e:S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
			   m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxP,LRikA,(double)LRiptA128,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
		  } else
		    printf("m=%d:I=%d,K=%d,J=%d,P=%d,T=%d,Q=%d..%d:FAik*FN=%0.10e,LRiptA=%0.10e:S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
			   m,I,K,J,P,T,max(J+1,Jmin[P]),QmaxP,LRikA,(double)LRiptA128,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2));
		  FL128 LRiptA2 = 0.0;
		  for(int Q = QminP /* WAS max(J+1,Jmin[P]) */; Q <= QmaxP; Q++){
		    double FAikjA = FAJ(Xm[Q]-Xm[J],deltaY,Q-J,P-T-I,FBiasWTm[Q-J], LogFpPow[Q-J-1], LogFnPow[P-T-I-1], PrBiasWT[P-T] - PrBiasWT[I], LRikA, Ivar,1, OUT_TYPE, OUT_OLD);
		    LRiptA2 += FAikjA * A.AR(P,T,Q);
		    printf("    Q=%d:FAikjA=%0.10e, AL=%0.10e(LS=%0.10f),AR=%0.10e:Y[I]=%0.3f,Y[P-T]=%0.3f,Xm[J]=%0.3f,Xm[Q]=%0.3f,sum=%0.10e(delta=%0.10e)\n",
			   Q,(double)FAikjA, (double)ALscale, LogScaleL[I]+LogScaleR[P]-LS, A.AR(P,T,Q), Y[I], Y[P-T], Xm[J], Xm[Q], (double)(LRiptA2*ALscale), (double)(FAikjA * A.AR(P,T,Q) * ALscale));
		  }
		  
		  fflush(stdout);
		}
		if(0 && MDEBUG && LDEBUG && DEBUG>=2 && !FIXADD){
		  FL128 addsum = LRadd[0], addsum2 = LRadd2[0];
		  for(int i = max(1,Amin); i <= I; i++){
		    addsum += LRadd[i];
		    addsum2 += LRadd2[i];
		  }
		  addsum -= addsum2;
		  FL128 fragsum = LRfrag[0], fragsum2 = LRfrag2[0];
		  for(int i = max(1,Amin); i <= I; i++){
		    fragsum += LRfrag[i];
		    fragsum2 += LRfrag2[i];
		  }
		  fragsum -= fragsum2;
		  if(DEBUG && !(addsum <= fragsum*FN[0]*(1.0+(USE_MFLOAT ? 1e-4 : 1e-6)) + (USE_MFLOAT ? 1e-20 : 1e-190) &&
		                addsum >= fragsum*FN[0]*(1.0-(USE_MFLOAT ? 1e-4 : 1e-6)) - (USE_MFLOAT ? 1e-20 : 1e-190))){
 		    printf("m=%d,I=%d,K=%d,J=%d,P=%d,T=%d: sum(LRadd[0..I)= %0.17e, sum(LRfrag[0..I]= %0.17e, sum(LRfrag[0..I])*Fn= %0.17e\n",
			   m,I,K,J,P,T,(double)addsum,(double)fragsum, (double)(fragsum * FN[0]));
		    printf("LRiptA128= %0.24Le, LRipt128= %0.24Le, LRipt128*Fn= %0.24Le\n",(long double)LRiptA128,(long double)LRipt128,(long double)(LRipt128 * FN[0]));
		    printf("LRadd[I]= %0.24Le, LRfrag[I]= %0.24Le, LRfrag[I]*Fn= %0.24Le\n",(long double)LRadd[I],(long double)LRfrag[I],(long double)(LRfrag[I] * FN[0]));
		    fflush(stdout);
		    assert(addsum <= fragsum*FN[0]*(1.0+(USE_MFLOAT ? 1e-4 : 1e-6)) + (USE_MFLOAT ? 1e-20 : 1e-190));
		    assert(addsum >= fragsum*FN[0]*(1.0-(USE_MFLOAT ? 1e-4 : 1e-6)) - (USE_MFLOAT ? 1e-20 : 1e-190));
		  }
		}
	      } // T = Tmin .. Tmax
	    } // P = Pmin .. min(IMAX, PmaxD + Kmax[P])

	    if(newLPdelta){  /* Repeat Inner Loops with S,d,P,T,Q to handle LPdelta (moved into seperate inline function mprobeval_delta(), see top of this file) */
#if USE_MIC==0
	      mprobeval_delta<OUT_TYPE,OUT_OLD>(A,m,I,K,J,N,D,deltaF,LRdeltaD,IMAX,Kmax,Jmin,Jmax,Pmax,Pmax2,Qmin,Qmax,Qmax2,Amin,Amax,Y,Xm,Yik,LogScaleL,LogScaleR,LS,AL,LRdelta);
#else
	      if(OUT_TYPE){
		if(OUT_OLD)
		  mprobeval_delta<1,1>(A,m,I,K,J,N,D,deltaF,LRdeltaD,IMAX,Kmax,Jmin,Jmax,Pmax,Pmax2,Qmin,Qmax,Qmax2,Amin,Amax,Y,Xm,Yik,LogScaleL,LogScaleR,LS,AL,LRdelta);
		else
		  mprobeval_delta<1,0>(A,m,I,K,J,N,D,deltaF,LRdeltaD,IMAX,Kmax,Jmin,Jmax,Pmax,Pmax2,Qmin,Qmax,Qmax2,Amin,Amax,Y,Xm,Yik,LogScaleL,LogScaleR,LS,AL,LRdelta);
	      } else {
		if(OUT_OLD)
		  mprobeval_delta<0,1>(A,m,I,K,J,N,D,deltaF,LRdeltaD,IMAX,Kmax,Jmin,Jmax,Pmax,Pmax2,Qmin,Qmax,Qmax2,Amin,Amax,Y,Xm,Yik,LogScaleL,LogScaleR,LS,AL,LRdelta);
		else
		  mprobeval_delta<0,0>(A,m,I,K,J,N,D,deltaF,LRdeltaD,IMAX,Kmax,Jmin,Jmax,Pmax,Pmax2,Qmin,Qmax,Qmax2,Amin,Amax,Y,Xm,Yik,LogScaleL,LogScaleR,LS,AL,LRdelta);
	      }
#endif
            }

	  }// J = Jmin[I] .. Jmax[I]
	  if(VERB>=2){
	    Jcnt += Jmax[I] - Jmin[I] + 1;
	    Kcnt++;
	  }
	} // K = 0 .. Kmax[I]
      } // I = IMIN .. IMAX

      if(VERB>=2 && newLPdelta && rverb){
	printf("mprobeval:m=%d,Amin=%d,Amax=%d:After FA updates:time=%0.6f\n",m,Amin,Amax,wtime());
	fflush(stdout);
      }

      if(MDEBUG_M >= 0){
	if(MDEBUG_S >= 0 && m==MDEBUG_M && N==MDEBUG_N && LPdelG){
	  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
	    sum += LRdel[i];
	  printf("m=%d:After FA updates:LRdel[0..%d]=%0.10e(LS=%0.10f,LL=%0.10f),IMIN=%d,IMAX=%d,N=%d\n",m,MDEBUG_S,(double)sum,LS,LS+log((double)sum),IMIN,IMAX,N);
	}
	if(LDEBUG && MDEBUG_F >= 0 && m==MDEBUG_M && N==MDEBUG_N-1){
	  FL128 sum = 0.0, sum2 = 0.0;
	  for(int i = 0; i <= MDEBUG_F; i++){
	    sum += LRfrag[i];
	    sum2 += LRfrag2[i];
	  }
	  printf("m=%d:After FA updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f,LL=%0.10f),IMIN=%d,IMAX=%d,N=%d\n",
		 m,MDEBUG_F,(double)sum,(double)sum2,(double)(sum-sum2),LS,LS+log((double)(sum-sum2)),IMIN,IMAX,N);
	}
	if(MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  FL128 sum = LRadd[0], sum2 = LRadd2[0];
	  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
	    sum += LRadd[i];
	    sum2 += LRadd2[i];
	  }
	  FL128 sum3 = sum;
	  if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
	    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After FA updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e,LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2),LS);
	  } else
	    printf("m=%d:After FA updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e(LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),LS);
	  if(LDEBUG){
	    sum = sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After FA updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f,LL=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS,LS+log((double)(sum-sum2)));
	  }
	}
	if((MDEBUG_ST >= 0 && m==MDEBUG_M) || (MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N+1)){
	  FL128 sum = 0.0;
	  int t = max(MDEBUG_ST,MDEBUG_TS+1);
	  for(int i = 1; i <= t; i++)
	    sum += LRsite[i]; 
	  printf("m=%d:After FA updates:LRsite[0..%d]=%0.10e(LS=%0.10f)\n",m,t,(double)sum,LS);
	}
	if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  if(newLPdelta && MDEBUG_TT < D[min(N,MDEBUG_TS)]){
	    FL128 sum = LRdelta[0];
	    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
	      sum += LRdelta[i];
	    FL128 sum3 = sum;
	    if(Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < D[MDEBUG_TS])
	      sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After FA updates:S=%d,D=%d,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e,LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3,LS);
	  }
	  if(LDEBUG){
	    FL128 sum = 0.0,sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After FA updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS);
	  }
	}
	fflush(stdout);
      }
      if(0 && MDEBUG && LDEBUG && !FIXADD){
	FL128 addsum = LRadd[0], addsum2 = LRadd2[0];;
	for(int i = 1; i <= 52;i++){
	  addsum += LRadd[i];
	  addsum2 += LRadd2[i];
	}
	addsum -= addsum2;
	FL128 fragsum = LRfrag[0], fragsum2 = LRfrag2[0];
	for(int i = 1; i <= 52; i++){
	  fragsum += LRfrag[i];
	  fragsum2 += LRfrag2[i];
	}
	fragsum -= fragsum2;
	if(DEBUG && (!(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190) || (VERB>=2&&m==1))){
	  printf("m=%d:After FA updates:sum(LRadd[0..52])=%0.17e,sum(LRfrag[0..52]=%0.17e\n",
		 m,(double)addsum,(double)fragsum);
	  fflush(stdout);
	  assert(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190);
	}
      }
      if(LDEBUG){/* check that cumulative value of LRsite[1..N+1] approximates 0 */
	FL128 LRsitecum = 0.0;
	for(int S = 1; S <= N+1; S++){
	  LRsitecum += LRsite[S];
	  if(VERB && (VERB>=2 || (DEBUG && !(S <= N || fabs((double)LRsitecum) <= LRcum * 1e-8)))){
	    printf("m=%d:S=%d,N=%d:LRsite[S]=%0.10e,LRsitecum=%0.10e,LRcum=%0.10e,LRfrag[0]=%0.10e,LRsitecum-LRcum=%0.10e,IMIN=%d,IMAX=%d\n",
		   m,S,N,(double)LRsite[S],(double)LRsitecum,LRcum,(double)LRfrag[0],(double)(LRsitecum-LRcum),IMIN,IMAX);
	    for(int s = 1; s <= N+1; s++)
	      printf("LRsite[%d]=%0.10e,LRfrag[%d]=%0.10e\n",s,(double)LRsite[s],s,(double)LRfrag[s]);
	    fflush(stdout);
	  }
	  if(DEBUG && S == N+1) assert(fabs((double)LRsitecum) <= LRcum * 1e-8);
	}
      }
    
      if(DEBUG>=2 && newLPdelta){
	if(VERB>=3){
	  if(m==7 && Amin <= 12 && 12 <= Amax && D[12] >= 1 /* && !isfinite(LRdeltaD[12][0])*/)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[12]=%d,delta[12][0]= %0.4f,LRdeltaD[12][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[12],delta[12][0],LRdeltaD[12][0]);
	  if(m==56 && Amin <= 31 && 31  <= Amax && D[31] >= 1 /* && !isfinite(LRdeltaD[31][0] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[31]=%d,delta[31][0]= %0.4f,LRdeltaD[31][0]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[31],delta[31][0],LRdeltaD[31][0]);
	  if(m==20 && Amin <= 33 && 33 <= Amax && D[33] >= 11 /* && !isfinite(LRdeltaD[33][10] */)
	    printf("\nm=%d:Amin=%d,Amax=%d,IMIN=%d,IMAX=%d,N=%d:D[33]=%d,delta[33][10]= %0.4f,LRdeltaD[33][10]= %0.8e\n",
		   m,Amin,Amax,IMIN,IMAX,N,D[33],delta[33][10],LRdeltaD[33][10]);
	  fflush(stdout);
	}
	for(int S = Amin;S <= min(N,Amax); S++)
	  for(int d = 0; d < D[S]; d++)
	    assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
      }

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:I=%d..%d,Kcnt=%lld,Jcnt=%lld, After Updates based on FA term:wall time= %0.6f\n", m, tid, IMIN,IMAX,Kcnt,Jcnt,wtime() - wt);
	  fflush(stdout);
        }
      }

      /* LRdel[S] or LRsite[S] Updates based on PM term with I-K < S < I (Line 4 of Equation 7) */
      /* LRaddT[S][] or LRfrag[S] Updates based on PM term with I-K <= S < I (Line 4 of Equation 10) */
      /* LRdeltaD[S][d] Updates based on PM term with I-K <= S < I (Line 4 of Equation 21) */
      for(int I= IMIN; I <= IMAX; I++){
	double ScaleL = exp(-LogScaleL[I]);
	double ScaleR = exp(-LogScaleR[I]);
	FL128b ScaleI = exp(LogScaleL[I]+LogScaleR[I]-LS);
	if(DEBUG)assert(isfinite(ScaleI));

	/* precompute ScaleLI[P=I-1 .. Pmin] = exp(LogScaleL[min(max(P,IMIN-1),IMAX+1)] - LogScaleL[I]) */
	int Pmin = max(IMIN-1, I - KMAX - deltaExtYRef);
	for(int P = I; --P >= Pmin; ){
	  if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
	  ScaleLI[P] = exp(LogScaleL[P] - LogScaleL[I]);
	  if(DEBUG>=2) assert(isfinite(ScaleLI[P]));
	}

	/* precompute ScaleRI[P= I+1 .. min(IMAX+1, I + KMAX + deltaExtYRef)] = exp(LogScaleR[P] - LogScaleR[I]) */
	int Pmax = min(IMAX+1, I + KMAX + deltaExtYRef);
	double LogScaleRI = LogScaleR[I];
	for(int P = I+1; P <= Pmax; P++ ){
	  if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
	  ScaleRI[P] = exp(LogScaleR[P] - LogScaleRI);
	  if(DEBUG>=2) assert(isfinite(ScaleRI[P]));		
	}

	for(int K = 1; K <= KMAX && K < I/*Kmax[I]*/; K++){// NOTE : for site addition, K can be larger than Kmax[I] or A.kmax
	  double Yik = Yc(Y,I,K);
	  int J;
	  if(newLPdelta){
	    /* LRdelta[S][d] : S range is I-K .. I-1 for modified sum(AL * AR) / PM(I,K) with Y[I] += delta[S][d] */
	    double origYI = Y[I];
	    double origYIK = Y[I-K];
	    for(int S = I-K; S < I;S++){
	      if(D[S] <= 0)// WAS10 if(!D[S])
		continue;
	      if(!KMAX_APPROX){
		if(origY[I] - origY[I-K] + delta[S][mindelta[S]] >= AlignResMax)
		  continue;/* assumes delta[S][d] values are in ascending order */

		if(!KMAX_FIX){	/* check that all intervals in Y[I-K .. I] (except Y[S..S+1]) are < maxdist */
		  int i = I-K;
		  for(; i < I; i++)
		    if(i != S && origY[i+1]-origY[i] >= maxdist)
		      break;
		  if(i < I)
		    continue;
		}
	      } else if(origY[I] - origY[I-K] + delta[S][mindelta[S]] >= maxdist)
		continue;
		
	      for(int d = 0; d < D[S]; d++){
		double del = delta[S][d];
		Y[I] = origYI + del;
		if(KMAX_APPROX ? (Y[I] - Y[I-K] < maxdist) : ((KMAX_FIX || del+origY[S+1]-origY[S] < maxdist) && Y[I]-Y[I-K] < AlignResMax)){/* compute Line 4 of Equation 21 */
		  /* expand Amin so Amin <= S */
		  while(Amin > S){
		    Amin--;
		    LRadd[Amin] = LRadd2[Amin] = 0.0;
		    LRdelta[Amin] = 0.0;
		    if(myaddcnt > 0)
		      for(int t = 0; t < addcntI[Amin]; t++)
			LRaddT[Amin][t] = 0.0;
		    if(newLPdelta)
		      for(int d = 0; d < D[Amin]; d++)
			LRdeltaD[Amin][d] = 0.0;
		  }
		  if(DEBUG>=2 && !(Amin <= S && S <= min(N,Amax))){
		    printf("m=%d:I=%d,K=%d:S=%d:Amin=%d,Amax=%d,N=%d\n",
			   m,I,K,S,Amin,Amax,N);
		    fflush(stdout);
		    assert(Amin <= S && S <= min(N,Amax));
		  }
		  if(SCORE_APPROX <= 1)
		    for(int t = S+1; t < I; t++)
		      Y[t] = origY[t] + del;
		  FL128b PMikD = PM(0,I,K,Y);
		  if(SCORE_APPROX <= 1)
		    for(int t = S+1; t < I; t++)
		      Y[t] = origY[t];

		  /* NOTE : If _OPENMP, AL(),AR() will NOT use PRtab[].Pr */

		  //		  for(RecB *p = &A[I][K][J=Jmin[I]]; J <= Jmax[I]; J++, p++)
		  for(J = Jmin[I]; J <= Jmax[I]; J++)
		    /* compute ALD[J] == AL[I][K][J] with modified Y[I] (line 3 of Equation 22) : does not include PMik */
		    ALD[J] = AL(m,I,K,J,0,-2,Yik+0.5*del,Xm,Y,M,N,lc,ScaleL,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax,OUT_TYPE, OUT_OLD);

		  Y[I] = origYI;
		  Y[I-K] = origYIK - del;
		  FL128b LRikD = 0.0;
		  for(J = Jmin[I]; J <= Jmax[I]; J++){
		    /* compute ARD[J] == AR[I][K][J] with modified Y[I-K] (line 4 of Equation 22) : does not include PMik */
		    ARD[J] = AR(m,I,K,J,0,-2,Yik-0.5*del,Xm,Y,M,N,rc,ScaleR,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax,OUT_TYPE,OUT_OLD);

		    /* sum ALD[J] * ARD[J] */
		    LRikD += ALD[J] * ARD[J];
		  }
		  Y[I-K] = origYIK;

		  FL128b PMmulD = ScaleI * LRikD * PMikD;
		  if(DEBUG>=2) assert(isfinite(PMmulD));
		  double origLRdeltaD = LRdeltaD[S][d];
		  LRdeltaD[S][d] += PMmulD;
		  if(DEBUG>=2 && !(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0) ){
		    #pragma omp critical
		    {
		      printf("\nm=%d:I=%d,K=%d,Amin=%d,Amax=%d,N=%d:S=%d,d=%d/%d(%0.3f):sum(ALD*ARD)=%0.10e,PMikD=%0.10e,PMmulD=%0.10e,LRdeltaD[S][d]=%0.10e -> %0.10e\n",
			     m,I,K,Amin,Amax,N,S,d,D[S],delta[S][d],(double)(LRikD*ScaleI),(double)PMikD,(double)PMmulD,(double)origLRdeltaD, LRdeltaD[S][d]);
		      fflush(stdout);
		      assert(isfinite(LRdeltaD[S][d]) && LRdeltaD[S][d] >= 0.0);
		    }
		  }
		  if(MDEBUG_M >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S == MDEBUG_TS && d == MDEBUG_TT /* && PMmulD >= MDEBUG_MIN*/){
		    FL128 sum = LRdelta[0];
		    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
		      sum += LRdelta[i];
		    FL128 sum3 = sum + LRdeltaD[MDEBUG_TS][MDEBUG_TT];
		    printf("m=%d:I=%d,K=%d:S=%d,d=%d(%0.3f):sum(ALD*ARD)=%0.10e,PMikD=%0.10e,PMmulD=%0.10e,LRdelta[0..S]=%0.10e,LRdeltaD[S][d]=%0.10e(sum=%0.10e)\n",
			   m,I,K,S,d,delta[S][d],(double)(LRikD*ScaleI),(double)PMikD,(double)PMmulD,(double)sum,(double)LRdeltaD[S][d],(double)sum3);
		  }
		}
	      } /* d = 0 .. D[S] -1 */
	    }/* S = I-K .. I-1 */
	    Y[I] = origYI;
	    Y[I-K] = origYIK;
	  }/* newLPdelta */

	  if(K > Kmax[I]){
	    if(SCORE_APPROX >= 2 || KMAX_APPROX || myaddcnt <= 0 || K >= KMAX/* WAS KMAX - 1 */ || Y[I] - Y[I-K] >= AlignResMax)
	      continue;

	    if(!KMAX_FIX){    /* KMAX_APPROX==0 : make sure all intervals (i=I-K..I-1) satisfy Y[i+1]-Y[i] < 2*maxdist && (Y[i+1]-Y[i] < maxdist || addcntI[i] > 0) */
	      int i = I-K;
	      if(DEBUG>=2) assert(i >= 1 && I <= N);
	      int bigcnt = 0;/* number of intervals >= maxdist */
	      for(;i < I; i++){
		double Fi = Y[i+1]-Y[i];
		if(Fi >= 2.0*maxdist || (Fi >= maxdist && (addcntI[i] <= 0 || addlocI[i][0]-Y[i] >= maxdist || Y[i+1]-addlocI[i][addcntI[i]-1] >= maxdist)))
		  break;
		if(Fi >= maxdist)
		  bigcnt++;
	      }
	      if(i < I || bigcnt >= 2)
		continue;
	    }
	  }

	  if(SCORE_APPROX <= 1 && K > Kmax[I]){/* special case to handle KMAX_APPROX==0 where added site Ht is inside interval Y[S..S+1] >= maxdist, but Ht-Y[S] < maxdist && Y[S+1]-Ht < maxdist */
	    if(DEBUG>=2) assert(SCORE_APPROX==1 && KMAX_APPROX==0 && addcnt > 0);

	    /* Cannot use A[I][K][J].AL or A[I][K][J].AR since they are only valid for K <= Kmax[I] */
	    /* Cannot use PRtab[K][I].PM since it may underflow for K > Kmax[I] */
	    
	    if(DEBUG>=2) assert(IMIN <= I && I <= IMAX);
	    double LRikPM = 0.0;
	    for(J = Jmin[max(IMIN,I-K)]/* WAS Jmin[I] */; J <= Jmax[I]; J++){
	      double ALijk = AL(m,I,K,J,0,-2,Yik,Xm,Y,M,N,lc,ScaleL,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax,OUT_TYPE,OUT_OLD);/* AL() does NOT include PMik */
	      double ARijk = AR(m,I,K,J,0,-2,Yik,Xm,Y,M,N,rc,ScaleR,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax,OUT_TYPE,OUT_OLD);/* AR() does NOT include PMik */
	      LRikPM += ALijk * ARijk;
	    }
	    LRikPM *= ScaleI;
	    if(DEBUG>=2) assert(isfinite(LRikPM));

	    /* expand LRadd[0,Amin..Amax] to include LRadd[I-K .. I] (similarly for LRdelta[]) */
	    if(DEBUG>=2) assert(Amax >= I);
	    int origAmin = Amin;
	    for(int i = Amin-1; i >= I-K; i--){
	      if(DEBUG>=2) assert(i > 0 && i <= N);
	      Amin = i;
	      LRadd[Amin] = LRadd2[Amin] = LRdelta[Amin] = 0.0;
	      if(myaddcnt > 0)
		for(int t = 0; t < addcntI[i]; t++)
		  LRaddT[i][t] = 0.0;
	      if(newLPdelta)
		for(int d = 0; d < ((D[i] + VPAD) & (~VPAD)); d++)
		  LRdeltaD[i][d] = 0.0;
	    }

	    if(DEBUG>=1+RELEASE && Amin != origAmin && !(Amin > IMIN - 100 || Y[Amin] > Y[IMIN] - 500.0)){
	      printf("WARNING:m=%d/%d: I=%d,K=%d,IMIN=%d,IMAX=%d: Amin -> %d,Amax=%d\n",m,MD,I,K,IMIN,IMAX,Amin,Amax);
	      fflush(stdout);
	      //	      assert((Amin > IMIN - 100 || Y[Amin] > Y[IMIN] - 500.0) && (Amax < IMAX + 100 || Y[Amax] < Y[IMAX] + 500.0));
	    }

	    /* No LRdel[] updates, since K > Kmax[I], hence 
	       If KMAX_FIX==0 : Some interval Y[I-K..I] must be >= maxdist, which remains true when a site S (I-K < S < I) is deleted
	       If KMAX_FIX==1 : The interval Y[I-K..I] must be >= AlignResMax, which remains true when a site S (I-K < S < I) is deleted */

	    if(/* addcnt && */ K < KMAX){ /* LRaddT[S][t=0..addcntI[S]-1] or LRfrag[S] Updates based on PM term with I-K <= S < I (Line 4 of Equation 10) */
	      FL128b PMik = PM(0,I,K,Y);
	      if(0 && MDEBUG_M >= 0 && MDEBUG_T >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && I-K <= MDEBUG_TS && MDEBUG_TS < I){
		printf("m=%d:I=%d,K=%d:sum(AL*AR)/(PMik*PMik)=%0.10e,PMik=%0.10e:S=%d,t=%d..%d\n",m,I,K,(double)LRikPM,(double)PMik,MDEBUG_TS,0,addcntI[MDEBUG_TS]-1);
		fflush(stdout);
	      }
	      for(int S = I-K; S < I; S++){
		if(KMAX_FIX==0){/* make sure all intervals in Y[I-K..I] except Y[S..S+1] are < maxdist */
		  int i = I-K;
		  for(;i < I; i++)
		    if(i != S && origY[i+1]-origY[i] >= maxdist)
		      break;
		  if(i < I)
		    continue;
		} else if(DEBUG>=2) assert(Y[I] - Y[I-K] <= AlignResMax);

		double PMmul = LRikPM * PMik;
		if(LDEBUG){
		  /* LRfrag[S] */
		  LRfrag[S] += PMmul;
		  LRfrag2[S+1] += PMmul;
		  if((MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N-1 && I-K <= MDEBUG_F && MDEBUG_F <= S && PMmul >= MDEBUG_MIN) || 
		     (MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && MDEBUG_TT >= 0 && S <= MDEBUG_TS && MDEBUG_TS <= S && PMmul >= MDEBUG_MIN)){
		    FL128 sum = 0.0, sum2 = 0.0;
		    for(int i = 0; i <= max(MDEBUG_F,MDEBUG_TS); i++){
		      sum += LRfrag[i];
		      sum2 += LRfrag2[i];
		    }
		    printf("m=%d:I=%d,K=%d:sum(AL*AR)/(PMik*PMik)=%0.10e,PMik=%0.10e,PMmul=%0.10e,LRfrag[0..%d]=%0.10e-%0.10e=%0.10e: Y[I]=%0.3f,Y[I-K]=%0.3f\n",
			   m,I,K,(double)LRikPM,(double)PMik,(double)PMmul,max(MDEBUG_F,MDEBUG_TS),(double)sum,(double)sum2,(double)(sum-sum2),Y[I],Y[I-K]);
		    fflush(stdout);
		  }
		}
		
		if(DEBUG>=2) assert(0 < S && Amin <= S && K >= 1 && S+1 <= Amax);

		for(int t = 0; t < addcntI[S]; t++){
		  double Ht = addlocI[S][t];
		  if(!KMAX_FIX && (Y[S+1]-Ht >= maxdist || Ht-Y[S] >= maxdist))
		    continue;
		  for(int r = I; r > S; r--)
		    Y[r+1] = origY[r];
		  Y[S+1] = Ht;
		  double PMikt = PM(0,I+1,K+1,Y);
		  PMmul = LRikPM * PMikt;
		  for(int r = I+1; r > S; r--)
		    Y[r] = origY[r];

		  /* LRaddT[S][t] */
		  LRaddT[S][t] += PMmul;
		  if(DEBUG>=1+RELEASE /* HERE HERE >=2 */ && !(isfinite(LRaddT[S][t]) && LRaddT[S][t] >= 0.0)){
		    #pragma omp critical
		    {
		      printf("\t m=%d:S=%d,t=%d(T=%d)Ht=%0.3f:I=%d,K=%d:LRaddT[S][t]=>%0.10e,PMmul=%0.10e\n",
			     m,S,t,addcntIcum[S]+t,Ht,I,K,(double)LRaddT[S][t],(double)PMmul);
		      fflush(stdout);
		      assert(isfinite(LRaddT[S][t]) && LRaddT[S][t] >= 0.0);
		    }
		  }
		  if(MDEBUG_M >= 0 && MDEBUG_T >= 0 && MDEBUG_TT == t && m==MDEBUG_M && N==MDEBUG_N && S <= MDEBUG_TS && MDEBUG_TS <= S && PMmul >= MDEBUG_MIN){
		    FL128 sum = LRadd[0],sum2 = LRadd2[0];
		    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
		      sum += LRadd[i];
		      sum2 += LRadd2[i];
		    }
		    FL128 sum3 = sum;
		    if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
		      sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
		      printf("m=%d:I=%d,K=%d,J=%d..%d:S=%d,T=%d:sum(AL*AR)/(PMik*PMik)=%0.10e(ScaleI=%0.10e),PMik=%0.10e,PMikt=%0.10e,PMmul=%0.10e:LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
			     m,I,K,Jmin[I-K],Jmax[I],MDEBUG_TS,t,(double)LRikPM,(double)ScaleI,(double)PMik,PMikt,(double)PMmul,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
		    } else
		      printf("m=%d:I=%d,K=%d,J=%d..%d:S=%d,T=%d:sum(AL*AR)/(PMik*PMik)=%0.10e(ScaleI=%0.10e),PMik=%0.10e,PMikt=%0.10e,PMmul=%0.10e:LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
			     m,I,K,Jmin[I-K],Jmax[I],MDEBUG_TS,t,(double)LRikPM,(double)ScaleI,(double)PMik,PMikt,(double)PMmul,(double)sum,(double)sum2,(double)(sum-sum2));
		  }/* if(MDEBUG_M>=0 ...) */
		}/* t = 0 .. addcntI[S]-1 */
	      }/* S = I-K .. I-1 */
	    }/* if(myaddcnt > 0 && K < KMAX) */
	  } else {/* K <= Kmax[I] */
	    if(DEBUG>=2) assert(K <= Kmax[I]);
	    double PMik = PRtab[K][I].PM;
	    if(DEBUG>=2) assert(PMik > 0.0);
	    if(DEBUG>=2) assert(fabs(log(PMik) - log(PM(0,I,K,Y))) <= (USE_MFLOAT ? 1e-5 : 2e-9));
	    FL128b LRik = 0.0;
	    for(J = Jmin[I]; J <= Jmax[I]; J++){
	      if(DEBUG>=2)assert(A.AR(I,K,J) >= 0.0);
	      if(DEBUG>=2)assert(A.AL(I,K,J) >= 0.0);
	      if(VERB>=2 && MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && I-K+1 <= MDEBUG_S && MDEBUG_S <= I-1){
		printf("m=%d:I=%d,K=%d,J=%d:AL=%0.10e,AR=%0.10e,PMik=%0.10e,LS=%0.6e,LogScaleL[I]=%0.6e,LogScaleR[I]=%0.6e\n", 
		       m,I,K,J,A.AL(I,K,J)*exp(LogScaleL[I]),A.AR(I,K,J)*exp(LogScaleR[I]),PMik,LS,LogScaleL[I],LogScaleR[I]);
		fflush(stdout);
	      }

	      LRik += A.AL(I,K,J) * A.AR(I,K,J);
	    }
	    if(DEBUG>=2 && !isfinite(LRik)){
	      #pragma omp critical
	      {
		printf("m=%d,I=%d,K=%d:LRik= %0.10e, ScaleI=%0.10e\n",m,I,K,(double)LRik,(double)ScaleI);
		fflush(stdout);
		assert(isfinite(LRik));
	      }
	    }
	    //	    FL128b origLRik = LRik;
	    LRik *= ScaleI;
	    if(DEBUG>=2 && !(isfinite(LRik))){
	      #pragma omp critical
	      {
		printf("m=%d,I=%d,K=%d:LRik= %0.10e, ScaleI=%0.10e\n",m,I,K,(double)LRik,(double)ScaleI);
		//		printf("\t origLRik= %0.10e\n",(double)origLRik);
		fflush(stdout);
		assert(isfinite(LRik));
	      }
	    }

	    if(LPdelG){ /* expand LRdel[Dmin..Dmax] to include LRdel[I-K+1..I] */
	      if(DEBUG>=2) assert(Dmax >= I);
	      if(DEBUG>=2) assert(I-K+1 > 0);
	      for(int i = Dmin-1; i >= I-K+1;i--)
		LRdel[Dmin=i] = 0.0;
	    }

	    /* expand LRadd[0,Amin..Amax] to include LRadd[I-K .. I] (similarly for LRdelta[]) */
	    if(DEBUG>=2) assert(Amax >= I);
	    int origAmin = Amin;
	    for(int i = Amin-1; i >= I-K; i--){
	      if(DEBUG>=2) assert(i > 0 && i <= N);
	      Amin = i;
	      LRadd[Amin] = LRadd2[Amin] = LRdelta[Amin] = 0.0;
	      if(myaddcnt > 0 && K < KMAX)
		for(int t = 0; t < addcntI[i]; t++)
		  LRaddT[i][t] = 0.0;
	      if(newLPdelta)
		for(int d = 0; d < ((D[i] + VPAD) & (~VPAD)); d++)
		  LRdeltaD[i][d] = 0.0;
	    }
	    if(DEBUG>=1+RELEASE && Amin < origAmin - 100 && !(Amin > IMIN - 100 || Y[Amin] > Y[IMIN] - 500.0)){
	      printf("WARNING:m=%d/%d: I=%d,K=%d,IMIN=%d,IMAX=%d: Amin=%d,Amax=%d(B)\n",m,MD,I,K,IMIN,IMAX,Amin,Amax);
	      fflush(stdout);
	      //	      assert((Amin > IMIN - 100 || Y[Amin] > Y[IMIN] - 500.0) && (Amax < IMAX + 100 || Y[Amax] < Y[IMAX] + 500.0));
	    }

	    if(SCORE_APPROX>=2){
	      FL128b PMmul = LRik / PMik;
	      if(LDEBUG){
		/* LRsite[S] : S range = I-K+1 .. I-1 */
		LRsite[I-K+1] += PMmul;
		LRsite[I] -= PMmul;
		if(LDEBUG && DEBUG>=2){/* verify that LRsite[1..N+1] adds up to approximately 0.0 */
		  FL128 sum = 0.0;
		  for(int S = 1; S <= N+1; S++)
		    sum += LRsite[S];
		  if((DEBUG && !(fabs((double)sum) <= LRcum * 1e-8)) ||(STRACE && m==MDEBUG_M && PMmul >= MDEBUG_MIN && (I-K+1==STRACE || I==STRACE))){
		    printf("m=%d,I=%d,K=%d,N=%d,M=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMmul=%0.10e,LRsite[1..N+1]=%0.10e,LRsite[%d]+=%0.10e,LRsite[%d]-=%0.10e\n",
			   m,I,K,N,M,(double)LRik,PMik,(double)PMmul,(double)sum,I-K+1,(double)LRsite[I-K+1],I,(double)LRsite[I]);
		    fflush(stdout);
		    assert(fabs((double)sum) <= LRcum * 1e-8);
		  }
		}
		if((MDEBUG_M < 0 || m==MDEBUG_M) && ((MDEBUG_ST >= 0 && I-K+1 <= MDEBUG_ST && MDEBUG_ST <= I-1)||
						     (MDEBUG_TS>=0 && MDEBUG_TT>=0 && I-K+1 <= MDEBUG_TS+1 && MDEBUG_TS+1 <= I-1 && m==MDEBUG_M && N==MDEBUG_N+1)) && PMmul >= MDEBUG_MIN){
		  FL128 sum = 0.0;
		  for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++)
		    sum += LRsite[i];
		  printf("m=%d,I=%d,K=%d,N=%d,M=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMmul=%0.10e,LRsite[1..%d]=%0.10e\n",
			 m,I,K,N,M,(double)LRik,PMik,(double)PMmul,max(MDEBUG_TS+1,MDEBUG_ST),(double)sum);
		}

		/* LRfrag[S] : S range = I-K .. I-1 */
		LRfrag[I-K] += PMmul;
		LRfrag2[I] += PMmul;
		if((MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N-1 && I-K <= MDEBUG_F && MDEBUG_F <= I-1 && PMmul >= MDEBUG_MIN) || 
		   (MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && MDEBUG_TT >= 0 && I-K <= MDEBUG_TS && MDEBUG_TS <= I-1 && PMmul >= MDEBUG_MIN)){
		  FL128 sum = 0.0, sum2 = 0.0;
		  for(int i = 0; i <= max(MDEBUG_F,MDEBUG_TS); i++){
		    sum += LRfrag[i];
		    sum2 += LRfrag2[i];
		  }
		  printf("m=%d:I=%d,K=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMmul=%0.10e,LRfrag[0..%d]=%0.10e-%0.10e=%0.10e: Y[I]=%0.3f,Y[I-K]=%0.3f\n",
			 m,I,K,(double)LRik,PMik,(double)PMmul,max(MDEBUG_F,MDEBUG_TS),(double)sum,(double)sum2,(double)(sum-sum2),Y[I],Y[I-K]);
		  fflush(stdout);
		}
	      }

	      if(LPdelG){ /* LRdel[S] : S range = I-K+1 .. I-1 */
		if(DEBUG>=2) assert(Dmin <= I-K+1 && I <= Dmax);
		if(DEBUG>=2) assert(isfinite(PMmul));
		LRdel[I-K+1] += PMmul;
		LRdel[I] -= PMmul;
		if(DEBUG>=2) assert(isfinite(LRdel[I-K+1]));
		if(DEBUG>=2) assert(isfinite(LRdel[I]));
		if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && I-K+1 <= MDEBUG_S && MDEBUG_S <= I-1 && PMmul >= MDEBUG_MIN){
		  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
		    sum += LRdel[i];
		  printf("m=%d:I=%d,K=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMd=%0.10e,LRdel[0..%d]=%0.10e: Y[I]=%0.3f,Y[I-K]=%0.3f\n",
			 m,I,K,(double)LRik,PMik, (double)PMmul, MDEBUG_S, (double)sum, Y[I], Y[I-K]);
		  fflush(stdout);
		}
	      }

	      if(myaddcnt > 0 && K < KMAX){
		/* LRadd[S] : S range is I-K .. I-1 */
		if(DEBUG>=2) assert(0 < I-K && Amin <= I-K && K >= 1 && I <= Amax);
		LRadd[I-K] += PMmul;
		LRadd2[I] += PMmul;
		if(MDEBUG_M >= 0 && MDEBUG_T >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && I-K <= MDEBUG_TS && MDEBUG_TS <= I && PMmul >= MDEBUG_MIN){
		  FL128 sum = LRadd[0],sum2 = LRadd2[0];
		  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
		    sum += LRadd[i];
		    sum2 += LRadd2[i];
		  }
		  FL128 sum3 = sum;
		  if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
		    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
		    printf("m=%d:I=%d,K=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMmul=%0.10e:S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
			   m,I,K,(double)LRik,PMik,(double)PMmul,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
		  } else
		    printf("m=%d:I=%d,K=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMmul=%0.10e:S=%d,T=%d:LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
			   m,I,K,(double)LRik,PMik,(double)PMmul,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2));
		}
	      }
	    } else {/* SCORE_APPROX <= 1 */
	      FL128b LRikPM = LRik / (PMik * PMik);

	      if(LPdelG || LDEBUG){/* LRdel[S] or LRsite[S] Updates based on PM term with I-K < S < I (Line 4 of Equation 7) */
		for(int S = I-K+1; S < I; S++){
		  for(int r = S; r < I; r++)
		    Y[r] = origY[r+1];
		  FL128b PMikD = PM(0,I-1,K-1,Y);
		  FL128b PMmulD = LRikPM * PMikD;
		  for(int r = S; r < I; r++)
		    Y[r] = origY[r];

		  if(DEBUG>=2) assert(isfinite(PMmulD));
		  if(LDEBUG){
		    FL128b PMmul = LRikPM * PMik;

		    /* LRsite[S] : single S value */
		    LRsite[S] += PMmul;
		    LRsite[S+1] -= PMmul;
		    if(LDEBUG && DEBUG>=2){/* verify that LRsite[1..N+1] adds up to approximately 0.0 */
		      FL128 sum = 0.0;
		      for(int S = 1; S <= N+1; S++)
			sum += LRsite[S];
		      if((DEBUG && !(fabs((double)sum) <= LRcum * 1e-8)) ||(STRACE && m==MDEBUG_M && PMmul >= MDEBUG_MIN && (I-K+1==STRACE || I==STRACE))){
			printf("m=%d,I=%d,K=%d,N=%d,M=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMmul=%0.10e,LRsite[1..N+1]=%0.10e,LRsite[%d]+=%0.10e,LRsite[%d]-=%0.10e\n",
			       m,I,K,N,M,(double)LRik,PMik,(double)PMmul,(double)sum,I-K+1,(double)LRsite[I-K+1],I,(double)LRsite[I]);
			fflush(stdout);
			assert(fabs((double)sum) <= LRcum * 1e-8);
		      }
		    }
		    if((MDEBUG_M < 0 || m==MDEBUG_M) && ((MDEBUG_ST >= 0 && S == MDEBUG_ST)||
							 (MDEBUG_TS>=0 && MDEBUG_TT>=0 && S == MDEBUG_TS+1 && m==MDEBUG_M && N==MDEBUG_N+1)) && PMmul >= MDEBUG_MIN){
		      FL128 sum = 0.0;
		      for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++)
			sum += LRsite[i];
		      int t = max(MDEBUG_TS+1,MDEBUG_ST);
		      printf("m=%d,I=%d,K=%d,J=%d..%d:S=%d,N=%d,M=%d:sum(AL*AR)=%0.10e(ScaleI=%0.10e),PMik=%0.10e,sum(AL*AR)/PMik=%0.10e,LRsite[1..%d]=%0.10e(LRsite[%d]=%0.10e,LRsite[%d]=%0.10e)\n",
			     m,I,K,Jmin[I],Jmax[I],S,N,M,(double)LRik,(double)ScaleI,(double)PMik,(double)PMmul,t,(double)sum,t,(double)LRsite[t],t+1,(double)LRsite[t+1]);
		    }
		  }

		  if(LPdelG){  /* LRdel[S] */
		    if(KMAX_APPROX==0 && origY[S+1]-origY[S-1] >= maxdist)
		      continue;/* PM(0,S+1,2,Y) will not be computed after deleting site S */

		    if(DEBUG>=2) assert(Dmin <= S && S <= Dmax);
		    if(DEBUG>=2) assert(isfinite(PMmulD));
		    LRdel[S] += PMmulD;
		    LRdel[S+1] -= PMmulD;
		    if(DEBUG>=2) assert(isfinite(LRdel[S]));
		    if(DEBUG>=2) assert(isfinite(LRdel[S+1]));

		    if(MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && S <= MDEBUG_S && MDEBUG_S <= S && PMmulD >= MDEBUG_MIN){
		      FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
		      for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
			sum += LRdel[i];
		      printf("m=%d:I=%d,K=%d,S=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMikD=%0.10e,PMmulD=%0.10e,LRdel[0..%d]=%0.10e: Y[I]=%0.3f,Y[I-K]=%0.3f\n",
			     m,I,K,S,(double)LRik, PMik, (double)PMikD, (double)PMmulD, MDEBUG_S, (double)sum, Y[I], Y[I-K]);
		      fflush(stdout);
		    }
		  }
		}/* S = I-K+1 .. I-1 */
	      }/* if(LPdelG) */

	      if(myaddcnt > 0 && K < KMAX){ /* LRaddT[S][t=0..addcntI[S]-1] or LRfrag[S] Updates based on PM term with I-K <= S < I (Line 4 of Equation 10) */
		if(0 && MDEBUG_M >= 0 && MDEBUG_T >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && I-K <= MDEBUG_TS && MDEBUG_TS < I){
		  printf("m=%d:I=%d,K=%d:sum(AL*AR)=%0.10e,PMik=%0.10e:S=%d,t=%d..%d\n",m,I,K,(double)LRik,PMik,MDEBUG_TS,0,addcntI[MDEBUG_TS]-1);
		  fflush(stdout);
		}
		double PMmul = LRikPM * PMik;

		for(int S = I-K; S < I; S++){
		  if(LDEBUG){
		    /* LRfrag[S] */
		    LRfrag[S] += PMmul;
		    LRfrag2[S+1] += PMmul;
		    if((MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N-1 && I-K <= MDEBUG_F && MDEBUG_F <= S && PMmul >= MDEBUG_MIN) || 
		       (MDEBUG_M >= 0 && m==MDEBUG_M && N==MDEBUG_N && MDEBUG_TT >= 0 && S == MDEBUG_TS && PMmul >= MDEBUG_MIN)){
		      FL128 sum = 0.0, sum2 = 0.0;
		      for(int i = 0; i <= max(MDEBUG_F,MDEBUG_TS); i++){
			sum += LRfrag[i];
			sum2 += LRfrag2[i];
		      }
		      printf("m=%d:I=%d,K=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMmul=%0.10e,LRfrag[0..%d]=%0.10e-%0.10e=%0.10e: Y[I]=%0.3f,Y[I-K]=%0.3f\n",
			     m,I,K,(double)LRik,PMik,(double)PMmul,max(MDEBUG_F,MDEBUG_TS),(double)sum,(double)sum2,(double)(sum-sum2),Y[I],Y[I-K]);
		      fflush(stdout);
		    }
		  }
		
		  if(DEBUG>=2) assert(0 < S && Amin <= S && K >= 1 && S+1 <= Amax);

		  for(int t = 0; t < addcntI[S]; t++){
		    for(int r = I; r > S; r--)
		      Y[r+1] = origY[r];
		    Y[S+1] = addlocI[S][t];
		    double PMikt = PM(0,I+1,K+1,Y);
		    if(DEBUG>=1+RELEASE && !(isfinite(PMikt) && PMikt >= 0.0)){
		      printf("\n m=%d:S=%d,t=%d(T=%d):I=%d,K=%d,N=%d:PMikt= %0.10e, Y[I+1]= %0.8f, Y[I-K]= %0.8f,Y[N+1]= %0.4f, pTP= %0.8f\n",
			     m,S,t,addcntIcum[S]+t,I,K,N,PMikt,Y[I+1],Y[I-K],Y[N+1],pTP);
		      double PMTk = PM(0,I+1,K+1,Y);
		      double PMT2k = PM2(0,I+1,K+1,Y);
		      printf("\t PM(0,I+1,K+1,Y)= %0.10e, PM2(0,I+1,K+1,Y)= %0.10e\n",PMTk,PMT2k);
		      fflush(stdout);
		      assert(isfinite(PMikt) && PMikt >= 0.0);
		    }
		    double PMmul = LRikPM * PMikt;
		    for(int r = I+1; r > S; r--)
		      Y[r] = origY[r];

		    /* LRaddT[S][t] */
		    LRaddT[S][t] += PMmul;
		    if(DEBUG>=1+RELEASE /* HERE HERE >=2 */ && !(isfinite(LRaddT[S][t]) && LRaddT[S][t] >= 0.0)){
		      #pragma omp critical
		      {
			printf("\n m=%d:S=%d,t=%d(T=%d):I=%d,K=%d:LRaddT[S][t]=>%0.10e,PMmul=%0.10e ( PMikt=%0.10e, LRikPM= %0.10e, LRik= %0.10e, PMik=%0.10e)\n",
			       m,S,t,addcntIcum[S]+t,I,K,(double)LRaddT[S][t],(double)PMmul, PMikt, (double)LRikPM, (double)LRik, (double)PMik);
			fflush(stdout);
			assert(isfinite(LRaddT[S][t]) && LRaddT[S][t] >= 0.0);
		      }
		    }
		    if(MDEBUG_M >= 0 && MDEBUG_T >= 0 && MDEBUG_TT == t && m==MDEBUG_M && N==MDEBUG_N && S <= MDEBUG_TS && MDEBUG_TS <= S && PMmul >= MDEBUG_MIN){
		      FL128 sum = LRadd[0],sum2 = LRadd2[0];
		      for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
			sum += LRadd[i];
			sum2 += LRadd2[i];
		      }
		      FL128 sum3 = sum;
		      if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
			sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
			printf("m=%d:I=%d,K=%d,S=%d,T=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMikt=%0.10e,PMmul=%0.10e:LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
			       m,I,K,MDEBUG_TS,t,(double)LRik,PMik,PMikt,(double)PMmul,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
		      } else
			printf("m=%d:I=%d,K=%d,S=%d,T=%d:sum(AL*AR)=%0.10e,PMik=%0.10e,PMikt=%0.10e,PMmul=%0.10e:LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
			       m,I,K,MDEBUG_TS,t,(double)LRik,PMik,PMikt,(double)PMmul,(double)sum,(double)sum2,(double)(sum-sum2));
		    }
		  }
		}/* S = I-K .. I-1 */
	      }/* if(myaddcnt > 0) */
	    }/* SCORE_APPROX <= 1 */
	  }/* K <= Kmax[I] */

	  if(0 && MDEBUG && LDEBUG && !FIXADD){
	    FL128 addsum = LRadd[0], addsum2 = LRadd2[0];
	    for(int i = 1; i <= 52;i++){
	      addsum += LRadd[i];
	      addsum2 += LRadd2[i];
	    }
	    addsum -= addsum2;
	    FL128 fragsum = LRfrag[0], fragsum2 = LRfrag2[0];
	    for(int i = 1; i <= 52; i++){
	      fragsum += LRfrag[i];
	      fragsum2 += LRfrag2[i];
	    }
	    fragsum -= fragsum2;
	    if(DEBUG && !(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190)){
	      printf("m=%d,I=%d,K=%d:PM internal:sum(LRadd[0..52])=%0.17e,sum(LRfrag[0..52]=%0.17e\n",
		     m,I,K,(double)addsum,(double)fragsum);
	      fflush(stdout);
	      assert(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190);
	    }
	  }
	}/* K = 1 .. Kmax[I] */
      } /* I = IMIN .. IMAX */

      if(VERB>=2 && newLPdelta && rverb){
	printf("mprobeval:m=%d,Amin=%d,Amax=%d:After PM internal updates:time=%0.6f\n",m,Amin,Amax,wtime());
	fflush(stdout);
      }

      if(MDEBUG_M >= 0){
	if(MDEBUG_S >= 0 && m==MDEBUG_M && N==MDEBUG_N && LPdelG){
	  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
	    sum += LRdel[i];
	  printf("m=%d:After PM internal updates:LRdel[0..%d]=%0.10e(LS=%0.10f),IMIN=%d,IMAX=%d,N=%d\n",m,MDEBUG_S,(double)sum,LS,IMIN,IMAX,N);
	}
	if(LDEBUG && MDEBUG_F >= 0 && m==MDEBUG_M && N==MDEBUG_N-1){
	  FL128 sum = 0.0, sum2 = 0.0;
	  for(int i = 0; i <= MDEBUG_F; i++){
	    sum += LRfrag[i];
	    sum2 += LRfrag2[i];
	  }
	  printf("m=%d:After PM internal updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f),IMIN=%d,IMAX=%d,N=%d\n",m,MDEBUG_F,(double)sum,(double)sum2,(double)(sum-sum2),LS,IMIN,IMAX,N);
	}
	if(MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  FL128 sum = LRadd[0], sum2 = LRadd2[0];
	  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
	    sum += LRadd[i];
	    sum2 += LRadd2[i];
	  }
	  FL128 sum3 = sum;
	  if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
	    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After PM internal updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e,LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2),LS);
	  } else
	    printf("m=%d:After PM internal updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e(LS=%0.10f),Amin=%d,Amax=%d,IMIN=%d,IMAX=%d\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),LS,Amin,Amax,IMIN,IMAX);
	  if(LDEBUG){
	    sum = sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After PM internal updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS);
	  }
	}
	if((MDEBUG_ST >= 0 && m==MDEBUG_M) || (MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N+1)){
	  FL128 sum = 0.0;
	  for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++)
	    sum += LRsite[i]; 
	  printf("m=%d:After PM internal updates:LRsite[0..%d]=%0.10e(LS=%0.10f)\n",m,max(MDEBUG_TS+1,MDEBUG_ST),(double)sum,LS);
	}
	if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  if(newLPdelta && MDEBUG_TT < D[min(N,MDEBUG_TS)]){
	    FL128 sum = LRdelta[0];
	    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
	      sum += LRdelta[i];
	    FL128 sum3 = sum;
	    if(Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < D[MDEBUG_TS])
	      sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After PM internal updates:S=%d,D=%d,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e,LS=%0.10f)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3, LS);
	  }
	  if(LDEBUG){
	    FL128 sum = 0.0,sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After PM internal updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2), LS);
	  }
	}
	fflush(stdout);
      }
      if(0 && MDEBUG && LDEBUG && !FIXADD){
	FL128 addsum = LRadd[0], addsum2 = LRadd2[0];
	for(int i = 1; i <= 52;i++){
	  addsum += LRadd[i];
	  addsum2 += LRadd2[i];
	}
	addsum -= addsum2;
	FL128 fragsum = LRfrag[0], fragsum2 = LRfrag2[0];
	for(int i = 1; i <= 52; i++){
	  fragsum += LRfrag[i];
	  fragsum2 += LRfrag2[i];
	}
	fragsum -= fragsum2;
	if(DEBUG && (!(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190) || (VERB>=2&&m==1))){
	  printf("m=%d:After PM internal updates:sum(LRadd[0..52])=%0.17e,sum(LRfrag[0..52]=%0.17e\n", m,(double)addsum,(double)fragsum);
	  fflush(stdout);
	  assert(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190);
	}
      }
      if(LDEBUG){/* check that cumulative value of LRsite[1..N+1] approximates 0 */
	FL128 LRsitecum = 0.0;
	for(int S = 1; S <= N+1; S++){
	  LRsitecum += LRsite[S];
	  if(VERB && (VERB>=2 || (DEBUG && !(S <= N || fabs((double)LRsitecum) <= LRcum * 1e-8)))){
	    printf("m=%d:S=%d,N=%d:LRsite[S]=%0.10e,LRsitecum=%0.10e,LRcum=%0.10e,LRfrag[0]=%0.10e,LRsitecum-LRcum=%0.10e,IMIN=%d,IMAX=%d\n",
		   m,S,N,(double)LRsite[S],(double)LRsitecum,LRcum,(double)LRfrag[0],(double)(LRsitecum-LRcum),IMIN,IMAX);
	    for(int s = 1; s <= N+1; s++)
	      printf("LRsite[%d]=%0.10e,LRfrag[%d]=%0.10e\n",s,(double)LRsite[s],s,(double)LRfrag[s]);
	    fflush(stdout);
	  }
	  if(DEBUG && S == N+1) assert(fabs((double)LRsitecum) <= LRcum * 1e-8);
	}
      }

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After PM internal updates:wall time= %0.6f\n", m, tid, wtime() - wt);
	  fflush(stdout);
        }
      }

      /* Update LRsite[S] based on PM term with S == I or S== I-K (Line 5 of Equation 6), no updates for LRdel[] */
      if(LDEBUG){
	for(int I= IMIN; I <= IMAX; I++){
	  FL128b ScaleI = exp(LogScaleL[I]+LogScaleR[I]-LS);
	  if(DEBUG)assert(isfinite(ScaleI) && ScaleI >= 0.0);
	  for(int K = 0; K <= Kmax[I]; K++){
	    FL128b LRik = 0.0;
	    int J = Jmin[I];
	    for(; J <= Jmax[I]; J++){
	      if(DEBUG>=2)assert(A.AR(I,K,J) >= 0.0);
	      if(DEBUG>=2)assert(A.AL(I,K,J) >= 0.0);
	      LRik += A.AL(I,K,J) * A.AR(I,K,J);
	    }
	    double PMik = PRtab[K][I].PM;
	    if(DEBUG>=2) assert(PMik > 0.0);
	    if(DEBUG>=2) assert(fabs(log(PMik) - log(PM(0,I,K,Y))) <= (USE_MFLOAT ? 1e-5 : 2e-9));
	    LRik *= ScaleI;
	    FL128b PMmul = LRik / PMik;
	    if(DEBUG>=2)assert(isfinite(PMmul) && PMmul >= 0.0);
	    LRsite[I-K] += PMmul;
	    LRsite[I+1] -= PMmul;
	    if(LDEBUG && DEBUG>=2){/* verify that LRsite[1..N+1] adds up to approximately 0.0 */
	      FL128 LRsitecum = 0.0;
	      for(int S = 1; S <= N+1; S++)
		LRsitecum += LRsite[S];
	      if((DEBUG && !(fabs((double)LRsitecum) <= LRcum * 1e-8)) ||(STRACE && m==MDEBUG_M && PMmul >= MDEBUG_MIN && (I-K==STRACE || I+1==STRACE))){
		printf("m=%d,I=%d,K=%d,J=%d..%d,N=%d,M=%d:sum(AR*AR)=%0.10e,PMik=%0.10e,PMmul=%0.10e,LRsite[1..N+1]=%0.10e,LRcum=%0.10e:LRsite[%d]+=%0.10e,LRsite[%d]-=%0.10e\n",
		       m,I,K,Jmin[I],Jmax[I],N,M,(double)LRik,PMik,(double)PMmul,(double)LRsitecum,(double)LRcum, I-K, (double)LRsite[I-K], I+1, (double)LRsite[I+1]);
		fflush(stdout);
		assert(fabs((double)LRsitecum) <= LRcum * 1e-8);
	      }
	    }
	    
	    if(K>=2){
	      LRsite[I-K+1] -= PMmul;
	      LRsite[I] += PMmul;
	      if(LDEBUG && DEBUG>=2){/* verify that LRsite[1..N+1] adds up to approximately 0.0 */
		FL128 LRsitecum = 0.0;
		for(int S = 1; S <= N+1; S++)
		  LRsitecum += LRsite[S];
		if((DEBUG && !(fabs((double)LRsitecum) <= LRcum * 1e-8)) ||(STRACE && m==MDEBUG_M && PMmul >= MDEBUG_MIN && (I-K+1==STRACE || I==STRACE))){
		  printf("m=%d,I=%d,K=%d,J=%d..%d,N=%d,M=%d:sum(AR*AR)=%0.10e,PMik=%0.10e,PMmul=%0.10e,LRsite[1..N+1]=%0.10e,LRcum=%0.10e:LRsite[%d]-=%0.10e, LRsite[%d]+=%0.10e\n",
			 m,I,K,Jmin[I],Jmax[I],N,M,(double)LRik,PMik,(double)PMmul,(double)LRsitecum,(double)LRcum, I-K+1, (double)LRsite[I-K+1],I, (double)LRsite[I]);
		  fflush(stdout);
		  assert(fabs((double)LRsitecum) <= LRcum * 1e-8);
		}
	      }
	    }

	    if((MDEBUG_M < 0 || m==MDEBUG_M) && ((I-K == MDEBUG_ST || MDEBUG_ST == I) || 
						 (m==MDEBUG_M && N==MDEBUG_N+1 && (I-K == MDEBUG_TS+1 || MDEBUG_TS+1 == I))) && PMmul >= MDEBUG_MIN){
	      FL128 sum = 0.0;
	      for(int i = 1; i <= max(MDEBUG_ST,MDEBUG_TS+1); i++)
		sum += LRsite[i];
	      int t = max(MDEBUG_ST,MDEBUG_TS+1);
	      printf("m=%d,I=%d,K=%d,J=%d..%d,N=%d,M=%d:Y[I]=%0.3f,Y[I-K]=%0.3f,sum(AL*AR)=%0.10e(ScaleI=%0.10e),PMik=%0.10e,sum(AL*AR)/PMik=%0.10e,LRsite[1..%d]=%0.10e(LRsite[%d]=%0.10e,LRsite[%d]=%0.10e)\n",
		     m,I,K,Jmin[I],Jmax[I],N,M,Y[I],Y[I-K],(double)LRik,(double)ScaleI,PMik,(double)PMmul,t,(double)sum,t,(double)LRsite[t],t+1,(double)LRsite[t+1]);
	      printf("\t LogScaleL[I]= %0.10f, LogScaleR[I]= %0.10f, LS= %0.10f\n", LogScaleL[I], LogScaleR[I], LS);

	      FL128b LRik2 = 0.0;
	      for(J=Jmin[I]; J <= Jmax[I]; J++){
		LRik2 += A.AL(I,K,J) * A.AR(I,K,J) * ScaleI / PMik;
		printf("    J=%d:AL=%0.10e(LS=%0.10f),AR=%0.10e,PMik=%0.10f:sum(AL*AR)/PMik=%0.10e(delta=%0.10e)\n",
		       J, A.AL(I,K,J), LogScaleL[I]+LogScaleR[I]-LS, A.AR(I,K,J), PMik,(double)LRik2, (double)(A.AL(I,K,J) * A.AR(I,K,J) * ScaleI / PMik));
	      }
	    }
	  }
	}
      }

      /* Update LRaddT[S=0..N][t=0..addcntI[S]-1] based on PM term with either end at new site (Line 5 & 6 of Equation 10), no updates for LRfrag[] or LRdeltaD[][] */
      if(myaddcnt > 0){
	int Smin = IMIN;
	int Smax = IMAX;
	if(FIX_RANGE_Y || RANGE_Y){
	  while(Smin > 1 && (KMAX_APPROX ? (Y[IMIN] - Y[Smin-1] < maxdist) : ((KMAX_FIX || Y[Smin]-Y[Smin-1] < maxdist) && Y[IMIN] - Y[Smin-1] < AlignResMax)))
	    Smin--;
	  while(Smax < N && (KMAX_APPROX ? (Y[Smax+1] - Y[IMAX] < maxdist) : ((KMAX_FIX || Y[Smax+1]-Y[Smax] < maxdist) && Y[Smax+1] - Y[IMAX] < AlignResMax)))
	    Smax++;
	}
	/* expand LRadd[0,Amin..Amax] to include Smin-1 .. Smax (also LRdelta[]) */
	int origAmin = Amin, origAmax = Amax;
	for(int i = Amin-1; i >= Smin-1; i--){
	  if((Amin = i) > 0){
	    LRadd[Amin] = LRadd2[Amin] = 0.0;
	    LRdelta[Amin] = 0.0;
	  }
	  if(DEBUG>=2) assert(i <= N);
	  if(1/* myaddcnt > 0 */)
	    for(int t = 0; t < addcntI[i]; t++)
	      LRaddT[i][t] = 0.0;
	  if(newLPdelta)
	    for(int d = 0; d < ((D[i] + VPAD) & (~VPAD)); d++)
	      LRdeltaD[i][d] = 0.0;
	}
	for(int i = Amax+1; i <= Smax; i++){
	  LRadd[Amax = i] = LRadd2[i] = LRdelta[i] = 0.0;
	  if(i <= N){
	    if(1 /* addcnt > 0 */)
	      for(int t = 0; t < addcntI[i]; t++)
		LRaddT[i][t] = 0.0;
	    if(newLPdelta)
	      for(int d = 0; d < ((D[i] + VPAD) & (~VPAD)); d++)
		LRdeltaD[i][d] = 0.0;
	  }
	}
	if(DEBUG>=1+RELEASE && !((Amin == origAmin || Amin > IMIN - 100 || Y[Amin] > Y[IMIN] - 500.0) && (Amax == origAmax || Amax < IMAX + 100 || Y[Amax] < Y[IMAX] + 500.0))){
          printf("WARNING:m=%d/%d:IMIN=%d,IMAX=%d,Smin=%d, Smax=%d, N=%d: Amin = %d,Amax = %d\n",m,MD,IMIN,IMAX,Smin,Smax,N,Amin,Amax);
	  fflush(stdout);
	  //	  assert((Amin > IMIN - 100 || Y[Amin] > Y[IMIN] - 500.0) && (Amax < IMAX + 100 || Y[Amax] < Y[IMAX] + 500.0));
	}

	if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && VERB>=2){
	  printf("Looking for updates to LRaddT[%d][%d]:IMIN=%d,IMAX=%d,Smin=%d,Smax=%d\n",MDEBUG_TS,MDEBUG_TT,IMIN,IMAX,Smin,Smax);
	  fflush(stdout);
	}

#if VECTOR2_ADD
#if USE_MIC==0
	mprobeval_add<OUT_TYPE, OUT_OLD>(A,m,M,N,lc,rc,IMIN,IMAX,Amin,Amax,Kmax,Jmin,Jmax,Y,Xm,origY,LogScaleL,LogScaleR,ScaleLI,ScaleRI,LS,
					 KtVec,KsVec,PMTH,ALTH,ARTH,LRadd,LRadd2,addlocI,addcntI,addcntIcum,LRaddT);
#else // USE_MIC==1
        if(OUT_TYPE){
	  if(OUT_OLD)
	    mprobeval_add<1,1>(A,m,M,N,lc,rc,IMIN,IMAX,Amin,Amax,Kmax,Jmin,Jmax,Y,Xm,origY,LogScaleL,LogScaleR,ScaleLI,ScaleRI,LS,
			       KtVec,KsVec,PMTH,ALTH,ARTH,LRadd,LRadd2,addlocI,addcntI,addcntIcum,LRaddT);
          else
	    mprobeval_add<1,0>(A,m,M,N,lc,rc,IMIN,IMAX,Amin,Amax,Kmax,Jmin,Jmax,Y,Xm,origY,LogScaleL,LogScaleR,ScaleLI,ScaleRI,LS,
			       KtVec,KsVec,PMTH,ALTH,ARTH,LRadd,LRadd2,addlocI,addcntI,addcntIcum,LRaddT);
	} else {
	  if(OUT_OLD)
	    mprobeval_add<0,1>(A,m,M,N,lc,rc,IMIN,IMAX,Amin,Amax,Kmax,Jmin,Jmax,Y,Xm,origY,LogScaleL,LogScaleR,ScaleLI,ScaleRI,LS,
			       KtVec,KsVec,PMTH,ALTH,ARTH,LRadd,LRadd2,addlocI,addcntI,addcntIcum,LRaddT);
          else
	    mprobeval_add<0,0>(A,m,M,N,lc,rc,IMIN,IMAX,Amin,Amax,Kmax,Jmin,Jmax,Y,Xm,origY,LogScaleL,LogScaleR,ScaleLI,ScaleRI,LS,
			       KtVec,KsVec,PMTH,ALTH,ARTH,LRadd,LRadd2,addlocI,addcntI,addcntIcum,LRaddT);
        }
#endif // USE_MIC==1
#endif // VECTOR2_ADD

	for(int S= Smin-1; S <= Smax; S++){
	  if(VECTOR2_ADD && IMIN <= S && S < IMAX)
	    continue;

	  FL128b ScaleS = exp(LogScaleL[min(max(S+1,IMIN-1),IMAX+1)]+LogScaleR[min(max(S,IMIN-1),IMAX+1)]-LS);
	  if(DEBUG>=2) assert(isfinite(ScaleS) && ScaleS >= 0.0);
	  if(DEBUG && !(Amin <= S && S <= Amax)){
	    printf("m=%d:S=%d,Amin=%d,Amax=%d,Smin=%d,IMIN=%d,IMAX=%d,Smax=%d\n",m,S,Amin,Amax,Smin,IMIN,IMAX,Smax);
	    fflush(stdout);
	    assert(Amin <= S && S <= Amax);
	  }
	  if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && VERB){
	    printf("Looking for updates to LRaddT[%d][%d]:S=%d,addcntIcum[S]=%d,addcntI[S]=%d (IMIN=%d,IMIN1=%d,IMAX1=%d,IMAX=%d,Smin=%d,Smax=%d)\n",
		   MDEBUG_TS,MDEBUG_TT,S,addcntIcum[S],addcntI[S],IMIN,IMIN1,IMAX1,IMAX,Smin,Smax);
	    fflush(stdout);
	  }
	
	  for(int t = 0; t < addcntI[S]; t++){
	    double Ht = addlocI[S][t];
	    if(DEBUG>=2) assert(Y[S] < Ht && Ht < Y[S+1]);

	    /* compute largest K so that Ht and Y[S+1-Kt] may fail to resolve */
	    if(DEBUG) assert(SCORE_APPROX);
	    /* determine Kt based on total distance Ht - Y[S+1-Kt] < maxdist (if KMAX_APPROX)*/
	    int Kt = 0;
	    while (S+1-Kt > 1){
	      if((KMAX_APPROX || (!KMAX_FIX && Kt==0)) ? (Ht - Y[S-Kt] >= maxdist) : ((!KMAX_FIX && Y[S+1-Kt] - Y[S-Kt] >= maxdist) || Ht - Y[S-Kt] >= AlignResMax))
		break;
	      Kt++;
	    }
	    if(Kt > KMAX){
	      if(VERB>=2){
		printf("WARNING:m=%d,S=%d,t=%d/%d:Ht=%0.3f,Kt=%d reduced to KMAX=%d,A.kmax=%d\n",m,S,t,addcntI[S],Ht,Kt,KMAX,A.kmax);
		fflush(stdout);
	      }
	      Kt = KMAX;
	    }

	    /* compute Jmin[] and Jmax[] for new site T */
	    /* NOTE : this will occasionally be a larger range than computed by setlimit() after site addition, except if S < IMIN or S+1 > IMAX */
	    int I1 = min(S,S+1-Kt);
	    int I2 = S+1;
	    int jmin = M+1, jmax = 0;
	    if((MDEBUG_SETLIMIT>=2 && MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT) ||(VERB>=2 && rverb && m==14)){
	      printf(" m=%d:S=%d,I1=%d,I2=%d,IMIN1=%d,IMAX1=%d,Jmid[IMIN1]=%d,Jmid[IMAX1]=%d,M=%d\n",m,S,I1,I2,IMIN1,IMAX1,Jmid[IMIN1],Jmid[IMAX1],M);
	      fflush(stdout);
	    }
	    if(GLOBAL || Jmid[IMIN1] == -1){
	      if(!GLOBAL && FIX_RANGE_Y && RANGE_Y==0){/* use code for IMIN > IMAX from setlimit() */
		int JR = M+1,JL = -1,iR= -1,iL= -1;		

		/* check for nearest site X[JL] left of X[JR] that is aligned with Hcuts[iL] <= Ht */
		for(JL = JR; --JL >= 1; )
		  if((iL=map[m][JL]) >= 0){
		    if(Hcuts[iL] <= Ht)
		      break;
		    JR = JL;
		    iR = iL;/* ignores mapK[JL] > 0 */
		  }
		double XI;
		if(JL <= 0){/* extrapolate left of Hcuts[iR] and X[JR] */
		  if(JR > M)/* revert to global alignment */
		    goto Lnoalignment;
		  XI = Xm[JR] - (Hcuts[iR] - Ht);
		} else if(JR > M){/* extrapolate right of Hcuts[iL] and X[JL] */
		  if(DEBUG) assert(JL >= 1);
		  XI = Xm[JL] + (Ht - Hcuts[iL]);
		} else {/* interpolate location of Y[I] between Hcuts[iL] and Hcuts[iR] */
		  if(DEBUG) assert(Hcuts[iL] <= Ht && Ht <= Hcuts[iR] && Hcuts[iL] < Hcuts[iR]);
		  XI = Xm[JL] + (Ht - Hcuts[iL])*(Xm[JR]-Xm[JL])/(Hcuts[iR]-Hcuts[iL]);
		}
		double Xleft = XI - Xrange;
		double Xright = XI + Xrange;

		/* convert Xleft,Xright to jmin,jmax values */
		if(Xleft > Xm[M+1] || Xright < 0.0){/* no alignment from Y[I] */
		  jmin = M+1;
		  jmax = 0;
		} else {
		  for(jmin = M; jmin > 1; jmin --)
		    if(Xm[jmin] <= Xleft)
		      break;
		  for(jmax = 1; jmax < M; jmax++)
		    if(Xm[jmax] >= Xright)
		      break;
		  if((MDEBUG_SETLIMIT>=2 && MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT) || (VERB>=2 && rverb && m==14)){
		    printf("m=%d:S=%d,Kt=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f:JL=%d,JR=%d,iL=%d,iR=%d,XI=%0.3f,Xleft=%0.3f,Xright=%0.3f,jmin=%d,jmax=%d,X[jmin]=%0.3f,X[jmax]=%0.3f,M=%d,X[M+1]=%0.3f\n",
			   m,S,Kt,Ht,Y[S],Y[S+1],JL,JR,iL,iR,XI,Xleft,Xright,jmin,jmax,Xm[jmin],Xm[jmax],M,Xm[M+1]);
		    fflush(stdout);
		  }
		}
	      } else {/* revert to global alignment */
	      Lnoalignment:
		jmin = 1;
		jmax = M;
	      }
	    } else { /* NOT global */
	      if(S < IMIN1){/* duplicate code from setlimit to compute Jmin[I1], Jmax[I1]  */
		int I = S;
		if(!FIX_RANGE_Y && !RANGE_Y){
		  if(DEBUG && !(I+1 >= IMIN)){
		    printf("m=%d,S=%d(%0.3f..%0.3f),t=%d(%0.3f),Kt=%d:I1=%d,I2=%d,IMIN=%d,IMIN1=%d,IMAX1=%d,IMAX=%d\n",
			   m,S,Y[S],Y[S+1],t,Ht,Kt,I1,I2,IMIN,IMIN1,IMAX1,IMAX);
		    fflush(stdout);
		    assert(I+1 >= IMIN);
		  }
		}
		double Xleft = -1e100, Xright = -1e100;

		/* first check if any site Xm[J] is already aligned with new site Ht */
		int J;
		for(J = 1; J <= M; J++){
		  int right = map[m][J];
		  if(right < 0)
		    continue;
		  if(DEBUG) assert(mapK[m][J] >= 0);
		  int left = right - mapK[m][J];
		  if(DEBUG) assert(left <= right);
		  /*		  if(right > 0)
		    while(right < n && Hcuts[right+1] <= Hcuts[right] + MAP_RES)
		    right++;*/
		  if(Hcuts[right] < Ht - 0.1*MAP_RES)
		    continue;

		  /*		  if(left <= n)
		    while(left > 1 && Hcuts[left-1] >= Hcuts[left] - MAP_RES)
		    left--;*/

		  if(Hcuts[left] > Ht + 0.1*MAP_RES){
		    J = M+1;
		    break;
		  }
		
		  if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT){
		    printf("    I1=%d:JL=%d,M=%d:RANGE=%d,Ht=%0.3f,map[m][JL]=%d,mapK[m][JL]=%d,Hcuts[%d..%d]=%0.3f..%0.3f,X[JL]=%0.3f\n",
			   I1,J,M,RANGE,Ht,map[m][J],mapK[m][J],left,right,Hcuts[left],Hcuts[right],Xm[J]);
		    fflush(stdout);
		  }
		  break;
		}
		if((MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT) || (VERB>=2 && rverb && m==14)){
		  printf("    I1=%d:J=%d,M=%d,Jmin[IMIN1]=%d,I=%d,RANGE_Y=%d,Y[S+1]=%0.3f,Xrange=%0.4f,Xrange2=%0.4f\n",
			 I1,J,M,Jmin[IMIN1],I,RANGE_Y,Y[S+1],Xrange,Xrange2);
		  if(J > M)
		    printf("    Y[IMIN1]=%0.3f,Ht=%0.3f,Xm[Jmin[IMIN1]]=%0.3f,Xm[1]=%0.3f,Jmin[IMIN1]=%d\n",
			   Y[IMIN1],Ht,Xm[Jmin[IMIN1]],Xm[1],Jmin[IMIN1]);
		  fflush(stdout);
		}

		if(J <= M){/* site J is already aligned with Ht */
		  if(FIX_RANGE_Y && RANGE_Y == 0){ /* only update Xleft,Xright (same as JXmin[I],JXmax[I] in setlimit) */
		    if(RefineXRange >= 2){
		      Xleft = Xm[J] - Xrange2;
		      Xright = min(JXmax[IMIN1] - 0.001, Xm[J] + Xrange2);
		    } else {
		      Xleft = Xm[J] - Xrange;
		      Xright = Xm[J] + Xrange;
		    }
		    if((MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT) || (VERB>=2 && rverb && m==14)){
		      printf("    Xm[J]=%0.3f,Xleft=%0.3f,Xright=%0.3f\n",Xm[J],Xleft,Xright);
		      fflush(stdout);
		    }
		  } else {
		    jmin = max(1,J - RANGE);
		    jmax = min(M, J + RANGE);
		  }
		} else if(!FIX_RANGE_Y){
		  if(DEBUG) assert(!EXTRAPOLATE);
		  /* use basic extrapolation from Jmin[IMIN1],Jmax[IMIN1] to compute jmin,jmax (same as Jmin[I],Jmax[I] in setlimit) */
		  if(RANGE_Y){
		    jmin = max(1, Jmin[IMIN1] - (IMIN1-I));
		    jmax = min(Jmax[IMIN1], Jmax[IMIN1] - (IMIN1-I) + max(0,RANGE_Y-RANGE));
		  } else {
		    for(jmin = Jmin[I+1]; jmin > 1; jmin--)
		      if(Y[IMIN1]-Ht < Xm[Jmin[IMIN1]] - Xm[jmin])
			break;
		    for(jmax = 0; jmax < Jmax[I+1]; jmax++)
		      if(Y[IMIN1]-Ht > Xm[Jmax[IMIN1]] - Xm[jmax])
			break;
		  }
		} else {/* only update Xleft,Xright (same as JXmin[I],JXmax[I] in setlimit) */
		  int JR = -1,JL = -1,iR= -1,kR = -1, iL= -1, kL = -1;
		  JR = Jmid[IMIN1];/* If Jmid[IMIN1] < 0, need to extrapolate in both directions, see case IMIN > IMAX in setlimit(), handled above */
		  if(DEBUG) assert(JR >= 1 && JR <= M);
		  iR = map[m][JR];//nmapMD[IMIN1];
		  if(DEBUG) assert(iR >= 1 && iR <= n);
		  kR = mapK[m][JR];
		  if(DEBUG) assert(0 <= kR && kR < iR);

		  double XI;
		  if(Hcuts[iR] <= Ht){/* can happen if Ht is very close to Y[IMIN1] */
		    XI = Xm[JR];
		  } else {

		    /* check for nearest site X[JL] left of X[JR] that is aligned with Hcuts[iL-kL .. iL] <= Ht */
		    for(JL = JR; --JL >= 1; ){
		      if((iL = map[m][JL]) >= 0){
			kL = mapK[m][JL];
			if(DEBUG) assert(0 <= kL && kL < iL);

			if(Hcuts[iL] <= Ht)
			  break;

			JR = JL;
			iR = iL;
			kR = kL;
		      }
		    }
		    if(JL >= 1){/* interpolate location of Ht between Hcuts[iL] and Hcuts[iR] */
		      if(DEBUG) assert(Hcuts[iL] <= Ht && Ht <= Hcuts[iR] && Hcuts[iL] < Hcuts[iR]);
		      XI = Xm[JL] + (Ht - Hcuts[iL])*(Xm[JR]-Xm[JL])/(Hcuts[iR] - Hcuts[iL]);
		    } else /* extraplate left of Hcuts[iR] and X[JR] */
		      XI = Xm[JR] - (Hcuts[iR] - Ht);
		  }

		  Xright = min(JXmax[IMIN1] - 0.001, XI + Xrange2);// NOTE : Xrange2 only applies in extension direction to accomodate large outlier
		  Xleft = XI - Xrange2;

		  if((DEBUG && (FIX_RANGE_Y && RANGE_Y==0) && !(Xright < JXmax[IMIN1])) || (VERB>=2 && rverb && m==14) ||
		     (MDEBUG_SETLIMIT>=2 && MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT)){
		    printf("  JL=%d,JR=%d,iL=%d,iR=%d:Ht=%0.3f,XI=%0.3f,Xleft=%0.3f,Xright=%0.3f,M=%d,n=%d,IMIN=%d,IMIN1=%d,IMAX1=%d,IMAX=%d,Jmid[IMIN1]=%d\n",
			   JL,JR,iL,iR, Ht,XI,Xleft,Xright,M,n,IMIN,IMIN1,IMAX1,IMAX,Jmid[IMIN1]);
		    fflush(stdout);
		    printf("X[JL]=%0.3f,X[JR]=%0.3f,Hcuts[iL]=%0.3f,Hcuts[iR]=%0.3f,X[Jmid[IMIN1]]=%0.3f\n",
			   Xm[JL],Xm[JR],Hcuts[max(0,iL)],Hcuts[iR], Xm[Jmid[IMIN1]]);
		    printf("m=%d:S=%d,IMIN1=%d,N=%d,Kt=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f:JXmin[IMIN1]=%0.3f,JXmax[IMIN1]=%0.3f,M=%d,X[M+1]=%0.3f\n",
			   m,S,IMIN1,N,Kt,Ht,Y[S],Y[S+1],JXmin[IMIN1],JXmax[IMIN1],M,Xm[M+1]);
		    fflush(stdout);
		    assert(Xright < JXmax[IMIN1]);
		  }
		}

		if(FIX_RANGE_Y && RANGE_Y==0){  /* convert Xleft,Xright to jmin,jmax values */
		  if(DEBUG) assert(Xleft > LARGE_NEGATIVE && Xleft <= Xright);
		  if(Xright < 0.0){
		    jmin = M+1;
		    jmax = 0;
		  } else {
		    for(jmin = M; jmin > 1; jmin--)
		      if(Xm[jmin] <= Xleft)
			break;
		    for(jmax = 1; jmax < M; jmax++)
		      if(Xm[jmax] >= Xright)
			break;
		    if((DEBUG>=2 && !(jmin <= jmax)) || (VERB>=2 && rverb && m==14) ||
		       (MDEBUG_SETLIMIT>=2 && MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT)){
		      printf("m=%d:S=%d,IMIN1=%d,N=%d,Kt=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f:JXmin[IMIN1]=%0.3f,JXmax[IMIN1]=%0.3f,Xleft=%0.3f,Xright=%0.3f,Jmid[IMIN1]=%d,jmin=%d,jmax=%d,M=%d,X[M+1]=%0.3f\n",
			     m,S,IMIN1,N,Kt,Ht,Y[S],Y[S+1],JXmin[IMIN1],JXmax[IMIN1],Xleft,Xright,Jmid[IMIN1],jmin,jmax,M,Xm[M+1]);
		      fflush(stdout);
		      assert(jmin <= jmax);
		    }
		  }
		}

		if(DEBUG && !(FIX_RANGE_Y && RANGE_Y==0) && !((1 <= jmin && jmin <= jmax && jmax <= M) || jmax <= 0)){
		  printf("    I1=%d:J=%d,M=%d,Jmin[IMIN1]=%d,I=%d,RANGE_Y=%d\n",
			 I1,J,M,Jmin[IMIN1],I,RANGE_Y);
		  if(J > M)
		    printf("Y[IMIN1]=%0.3f,Ht=%0.3f,Xm[Jmin[IMIN1]]=%0.3f,Xm[jmin]=%0.3f\n",
			   Y[IMIN1],Ht,Xm[Jmin[IMIN1]],Xm[jmin]);
		  printf("    I1=%d:JR=%d,M=%d,Jmax[IMIN1]=%d,I=%d,RANGE_Y=%d,RANGE=%d\n",
			 I1,J,M,Jmax[IMIN1],I,RANGE_Y,RANGE);
		  if(J <= 0)
		    printf("Y[IMIN1]=%0.3f,Ht=%0.3f,Xm[Jmax[IMIN1]]=%0.3f,Xm[jmax]=%0.3f\n",
			   Y[IMIN1],Ht,Xm[Jmax[IMIN1]],Xm[jmax]);
		  printf("m=%d,S=%d,Kt=%d,I1=%d,I2=%d:jmin=%d,jmax=%d:IMIN1=%d,IMAX1=%d\n",
			 m,S,Kt,I1,I2,jmin,jmax,IMIN1,IMAX1);
		  fflush(stdout);
		  assert((1 <= jmin && jmin <= jmax && jmax <= M) || jmax <= 0);
		}
	      } else if(S+1 > IMAX1){/* duplicate code from setlimit to compute Jmin[I2],Jmax[I2] */
		int I = S+1;
		if(DEBUG && !FIX_RANGE_Y && !RANGE_Y) assert(I-1 <= IMAX);
		double Xleft = -1e100,Xright = -1e100;

		/* first check if any site Xm[J] is already aligned with new site Ht */
		int J;
		for(J = 1; J <= M; J++){
		  int right = map[m][J];
		  if(right < 0)
		    continue;
		  if(DEBUG) assert(mapK[m][J] >= 0);
		  int left = right - mapK[m][J];
		  if(DEBUG) assert(left <= right);
		  /*		  if(right > 0)
		    while(right < n && Hcuts[right+1] <= Hcuts[right] + MAP_RES)
		    right++;*/
		  if(Hcuts[right] < Ht - 0.1*MAP_RES)
		    continue;

		  /*		  if(left <= n)
		    while(left > 1 && Hcuts[left-1] >= Hcuts[left] - MAP_RES)
		    left--;*/

		  if(Hcuts[left] > Ht + 0.1*MAP_RES){
		    J = M+1;
		    break;
		  }

		  break;
		}
		if((MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT) || (VERB>=2 && rverb && m==14)){
		  printf("    I2=%d:J=%d,M=%d,Jmin[IMAX1]=%d,I=%d,RANGE_Y=%d,RANGE=%d,Xrange=%0.4f,Xrange2=%0.4f\n",
			 I2,J,M,Jmin[IMAX1],I,RANGE_Y,RANGE,Xrange,Xrange2);
		  if(J > M)
		    printf("Y[IMAX1]=%0.3f,Ht=%0.3f,Xm[Jmin[IMAX1]]=%0.3f,Xm[M+1]=%0.3f,Jmin[IMAX1]=%d\n",
			   Y[IMAX1],Ht,Xm[Jmin[IMAX1]],Xm[M+1],Jmin[IMAX1]);
		}

		if(J <= M){/* site J is already aligned with Ht */
		  if(FIX_RANGE_Y && RANGE_Y==0){/* only update Xleft,Xright */
		    if(RefineXRange >= 2){/* duplicate Xrange expansion in setlimi() */
		      Xleft = max(JXmin[IMAX1] + 0.001, Xm[J] - Xrange2);
		      Xright = Xm[J] + Xrange2;
		    } else {
		      Xleft = Xm[J] - Xrange;
		      Xright = Xm[J] + Xrange;
		    }
		  } else {
		    jmin = max(1,J - RANGE);
		    jmax = min(M, J + RANGE);
		  }
		} else {
		  if(DEBUG) assert(!EXTRAPOLATE);// EXTRAPOLATE not implemented

		  /* use basic extrapolation from Jmin[IMAX1],Jmax[IMAX1] */
		  if(RANGE_Y){
		    jmin = max(Jmin[IMAX1],Jmin[IMAX1] + (I-IMAX1) - max(0,RANGE_Y-RANGE));
		    jmax = min(M, Jmax[IMAX1] + (I-IMAX1));
		  } else if(!FIX_RANGE_Y) {
		    for(jmin = M+1; jmin > Jmin[I-1]; jmin--)
		      if(Ht-Y[IMAX1] > Xm[jmin] - Xm[Jmin[IMAX1]])
			break;
		    for(jmax = Jmax[I-1]; jmax < M; jmax++)
		      if(Ht-Y[IMAX1] < Xm[jmax] - Xm[Jmax[IMAX1]])
			break;
		  } else {/* only update Xleft,Xright (same as JXmin[I],JXmax[I] in setlimit) */
		    int JR = -1,JL = -1,iR= -1,iL= -1;
		    JL = Jmid[IMAX1];/* if Jmid[IMAX1] < 0, need to extrapolate in both directions, see IMIN > IMAX in setlimit(), handled above */
		    if(DEBUG) assert(JL >= 1 && JL <= M);
		    iL = map[m][JL];//nmapMD[IMAX1];
		    if(DEBUG && !(iL >= 1 && iL <= n)){
		      printf("\n m=%d:S=%d,IMAX1=%d,N=%d,Kt=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f:JL=Jmid[IMAX1]=%d,iL=map[m][JL]=%d,n=%d\n",
			     m,S,IMAX1,N,Kt,Ht,Y[S],Y[S+1],JL,iL,n);
		      fflush(stdout);
		      assert(iL >= 1 && iL <= n);
		    }

		    double XI;
		    if(Hcuts[iL] >= Ht){/* can happen if Ht is very close to Y[IMAX1] */
		      XI = Xm[JL];
		    } else {
		      /* check for first site X[JR] right of X[JL] that is aligned with Hcuts[iR] >= Ht */
		      for(JR = JL; ++JR <= M;)
			if((iR=map[m][JR]) >= 0){
			  if(Hcuts[iR] >= Ht)
			    break;
			  /* ignore mapK[JL] > 0 */
			  JL = JR;
			  iL = iR;
			}

		      if(JR <= M){/* interpolate location of Ht between Hcuts[iL] and Hcuts[iR] */
			if(DEBUG && !(Hcuts[iL] <= Ht && Ht <= Hcuts[iR] && Hcuts[iL] < Hcuts[iR])){
			  printf("\n m=%d:S=%d,IMAX1=%d,N=%d,Kt=%d,Y[S]=%0.4f,Y[S+1]=%0.4f:Jmid[IMAX1]=%d,JL=%d,iL=%d,n=%d,M=%d,JR=%d,iR=%d\n",
				 m,S,IMAX1,N,Kt,Y[S],Y[S+1],Jmid[IMAX1],JL,iL,n,M,JR,iR);
			  printf("\t Hcuts[iL]= %0.6f, Ht= %0.6f, Hcuts[iR]= %0.6f, Xm[JL]= %0.6f, Xm[JR]= %0.6f\n",Hcuts[iL],Ht,Hcuts[iR],Xm[JL],Xm[JR]);
			  fflush(stdout);
			  assert(Hcuts[iL] <= Ht && Ht <= Hcuts[iR] && Hcuts[iL] < Hcuts[iR]);
			}
			XI = Xm[JL] + (Ht - Hcuts[iL])*(Xm[JR]-Xm[JL])/(Hcuts[iR]-Hcuts[iL]);
		      } else /* extrapolate right of Hcuts[iL] and X[JL] */
			XI = Xm[JL] + (Ht - Hcuts[iL]);
		    }

		    Xleft = max(JXmin[IMAX1] + 0.001, XI - Xrange2);// NOTE : Xrange2 only applies in extension direction to accomodate large outlier
		    Xright = XI + Xrange2;
		  }
		}	      

		if(FIX_RANGE_Y && RANGE_Y==0){/* convert Xleft,Xright to jmin,jmax values */
		  if(Xleft > Xm[M+1]){
		    jmin = M+1;
		    jmax = 0;
		  } else {
		    for(jmin = M; jmin > 1; jmin --)
		      if(Xm[jmin] <= Xleft)
			break;
		    for(jmax = 1; jmax < M; jmax++)
		      if(Xm[jmax] >= Xright)
			break;
		    if((DEBUG>=2 && !(jmin <= jmax)) || (VERB>=2 && rverb && m==14) ||
		       (MDEBUG_SETLIMIT>=2 && MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT)){
		      printf("m=%d:S=%d,IMAX1=%d,N=%d,Kt=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f:JXmin[IMAX1]=%0.3f,JXmax[IMAX1]=%0.3f,Xleft=%0.3f,Xright=%0.3f,jmin=%d,jmax=%d,M=%d,X[M+1]=%0.3f\n",
			     m,S,IMAX1,N,Kt,Ht,Y[S],Y[S+1],JXmin[IMAX1],JXmax[IMAX1], Xleft,Xright,jmin,jmax,M,Xm[M+1]);
		      fflush(stdout);
		      assert(1 <= jmin && jmin <= jmax && jmax <= M);
		    }
		  }
		}
		if(DEBUG && !(FIX_RANGE_Y && RANGE_Y==0) && !((1 <= jmin && jmin <= jmax && jmax <= M) || jmin > M)){
		  printf("m=%d,S=%d,Kt=%d,I1=%d,I2=%d:jmin=%d,jmax=%d:IMIN1=%d,IMAX1=%d,IMIN=%d,IMAX=%d\n",
			 m,S,Kt,I1,I2,jmin,jmax,IMIN1,IMAX1,IMIN,IMAX);
		  fflush(stdout);
		  assert((1 <= jmin && jmin <= jmax && jmax <= M) || jmin > M);
		}
	      } else {/* internal interval : duplicate code from setlimit */
		if(DEBUG) assert(!EXTRAPOLATE);// EXTRAPOLATE not implemented
		int L = S,R = S+1;
		if((MDEBUG_SETLIMIT>=2 && MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT && VERB>=2) ||
		   (VERB>=2 && rverb && m==14)){
		  printf(" m=%d:S=%d,L=%d,R=%d,RANGE_Y=%d,IMIN1=%d,IMAX1=%d\n",m,S,L,R,RANGE_Y,IMIN1,IMAX1);
		  fflush(stdout);
		}
		if(RANGE_Y){
		  jmin = min(Jmin[L] + 1 - max(0,RANGE_Y-RANGE), Jmin[R] - 1);
		  jmax = max(Jmax[L] + 1, Jmax[R] - 1 + max(0,RANGE_Y-RANGE));
		  jmin = max(jmin,Jmin[L]);
		  jmax = min(jmax,Jmax[R]);
		} else if(!FIX_RANGE_Y){
		  for(jmin = Jmin[R]; jmin > 1; jmin--)
		    if(!(Ht - Y[L] <= Xm[jmin] - Xm[Jmin[L]] || Y[R]- Ht >= Xm[Jmin[R]] - Xm[jmin]))
		      break;
		  for(jmax = Jmax[L]; jmax < M; jmax++)
		    if(!(Ht-Y[L] >= Xm[jmax] - Xm[Jmax[L]] || Y[R]-Ht <= Xm[Jmax[R]] - Xm[jmax]))
		      break;
		  jmin = max(jmin,Jmin[L]);
		  jmax = min(jmax,Jmax[R]);
		} else {
		  /* First Locate nearest left alignment L <= S s.t. Jmid[L] > 0,
		     and nearest right alignment R > S s.t. Jmid[R] > 0 */
                  int JL = 0, JR = M+1;
		  for(; L >= IMIN1; L--)
		    if((JL = Jmid[L]) > 0)
		      break;
		  if(DEBUG) assert(L >= IMIN1 && Jmid[L] > 0);
		  if(DEBUG) assert(JL > 0 && JL <= M);
	      
		  for(; R <= IMAX1; R++)
		    if((JR = Jmid[R]) > 0)
		      break;
		  if(DEBUG) assert(R <= IMAX1 && Jmid[R] > 0);
		  if(DEBUG) assert(JR > 0 && JR <= M);

		  /* next check if any site Xm[J] (JL <= J <= JR) is already aligned with new site Ht */
		  int J;
		  for(J = JL; J <= JR; J++){
		    int right = map[m][J];
		    if(right < 0)
		      continue;
		    if(DEBUG) assert(mapK[m][J] >= 0);
		    int left = right - mapK[m][J];
		    if(DEBUG) assert(left <= right);
		    if(Hcuts[right] < Ht - 0.1*MAP_RES)
		      continue;

		    if(Hcuts[left] > Ht + 0.1*MAP_RES){
		      J = M+1;
		      break;
		    }
		    break;
		  }

		  if((DEBUG && JL >= JR && J != JL) || (VERB>=2 && rverb && m==14) ||
		     (MDEBUG_SETLIMIT>=2 && MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT)){
		    #pragma omp critical 
		    {
		      printf("\n m=%d,S=%d,Kt=%d,Ht=%0.3f:Y[S]=%0.3f,Y[S+1]=%0.3f:L=%d,R=%d,JL=%d,JR=%d,J=%d,M=%d,N=%d,IMIN1=%d,IMAX1=%d,Jmid[L]=%d,Jmid[R]=%d\n",
			     m,S,Kt,Ht,Y[S],Y[S+1],L,R,JL,JR,J,M,N,IMIN1,IMAX1,Jmid[L],Jmid[R]);
		      printf(" map[m][JL]=%d,mapK[m][JL]=%d, Hcuts[%d]=%0.3f,Hcuts[%d]=%0.3f,Ht=%0.3f,Xm[J]=%0.3f,Xrange=%0.3f\n",
			     map[m][JL],mapK[m][JL], map[m][JL],Hcuts[map[m][JL]],map[m][JL]-mapK[m][JL], Hcuts[map[m][JL]-mapK[m][JL]], Ht, Xm[J],Xrange);
		      if(VERB>=2)
			for(int i = IMIN1; i <= IMAX1; i++)
			  printf("\t Jmid[%d]= %d\n",i,Jmid[i]);
		      fflush(stdout);

		      if(DEBUG && JL >= JR && J != JL){
			for(int k = map[m][JL]; ++k <= n;){
			  printf("Hcuts[%d]=%0.3f\n",k,Hcuts[k]);
			  if(Hcuts[k] > Y[S+1])
			    break;
			}

			for(int i = IMIN1; i <= IMAX1; i++)
			  if(Jmid[i] > 0){
			    int j = Jmid[i];
			    if(j > MX[m]){
			      printf("Jmid[%d]=%d,MX[m]=%d ???\n",i,j,MX[m]);
			      continue;
			    }
			    int right = map[m][j];
			    int left = right - mapK[m][j];
			    printf("Jmid[%d]=%d=j,MX[m]=%d:X[m][j]=%0.3f,map[m][j]=%d,mapK[m][j]=%d:\n",i,j,MX[m],Xm[j],map[m][j],mapK[m][j]);
			    for(int t = left; t <= right; t++)
			      if(map[MD][t]>=0)
				printf("\t t=%d:map[MD][t]=%d, Hcuts[t]=%0.3f,Y[map[MD][t]]=%0.3f\n",
				       t,map[MD][t],Hcuts[t],Y[map[MD][t]]);
			  }
			fflush(stdout);
			assert(J==JL);
		      }
		    }
		  }

		  double Xleft,Xright;/* update Xleft,Xright based on distances (same as JXmin[],JXmax[] in setlimit) */

		  if(J <= JR){/* site J is already aligned with Ht */
		    Xleft = Xm[J] - Xrange;
		    Xright = Xm[J] + Xrange;
		  } else { 
		    /* locate smallest range Hcuts[iL..iR] that includes Ht and is aligned with X[JL..JR] */
		    if(DEBUG && !(1 <= JL && JL <= JR && JR <= M)){
		      #pragma omp critical
		      {
			printf("\nm=%d,S=%d,Kt=%d,Ht=%0.3f:Y[S]=%0.3f,Y[S+1]=%0.3f:JL=%d,JR=%d,M=%d,N=%d\n",
			       m,S,Kt,Ht,Y[S],Y[S+1],JL,JR,M,N);
			fflush(stdout);
			assert(1 <= JL && JL <= JR && JR <= M);
		      }
		    }
		    int iL = map[m][JL] , iR = map[m][JR];
		    if(DEBUG) assert(1 <= iL && iL <= iR && iR <= n);
		    for(int J = JL; ++J < JR;){
		      int i = map[m][J];
		      if(i >= 0){
			if(Hcuts[i] > Ht)
			  break;
			JL = J;
			iL = i;
		      }
		    }
		    for(int J = JR; --J > JL;){
		      int i = map[m][J];
		      if(i >= 0){
			if(Hcuts[i] < Ht)
			  break;
			JR = J;
			iR = i;
		      }
		    }
		    if(DEBUG && !(JL < JR && iL < iR && Hcuts[iL] <= Ht && Ht <= Hcuts[iR])){
		      printf("\n    S=%d,t=%d,Y[S]=%0.3f,Y[S+1]=%0.3f:L=%d,R=%d,JL=%d,JR=%d,iL=%d,iR=%d,X[JL]=%0.3f,X[JR]=%0.3f,Hcuts[iL]=%0.3f,Hcuts[iR]=%0.3f,Ht=%0.3f,M=%d,n=%d\n",
			     S,t,Y[S],Y[S+1],L,R,JL,JR,iL,iR,Xm[JL],Xm[JR],Hcuts[iL],Hcuts[iR],Ht,M,n);
		      printf("    Jmid[S]=%d,X[Jmid[S]]=%0.3f,Jmid[S+1]=%d,X[Jmid[S+1]=%0.3f\n",
			     Jmid[S],Xm[Jmid[S]],Jmid[S+1],Xm[Jmid[S+1]]);
		      printf("m=%d:S=%d,Kt=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f,JXmin[S]=%0.3f,JXmax[S]=%0.3f,JXmin[S+1]=%0.3f,JXmax[S+1]=%0.3f\n",
			     m,S,Kt,Ht,Y[L],Y[R],JXmin[L],JXmax[L],JXmin[R],JXmax[R]);
		      printf("    map[m][JL]=%d,map[m][JR]=%d\n",map[m][JL],map[m][JR]);
		      fflush(stdout);
		      assert(JL < JR && iL < iR && Hcuts[iL] <= Ht && Ht <= Hcuts[iR]);
		    }

		    /* interpolate between Hcuts[iL] and Hcuts[iR] */
		    double XI = Xm[JL] + (Ht - Hcuts[iL])*(Xm[JR] - Xm[JL])/(Hcuts[iR]-Hcuts[iL]);
		    Xleft = max(JXmin[L] + 0.001, XI - Xrange2);
		    Xright = min(JXmax[R] - 0.001, XI + Xrange2);
		    if((DEBUG && !(Xleft <= Xright)) || 
		       (MDEBUG_SETLIMIT>=2 && MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT)){
		      printf("m=%d:S=%d,Kt=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f,JXmin[S]=%0.3f,JXmax[S]=%0.3f,JXmin[S+1]=%0.3f,JXmax[S+1]=%0.3f,Xleft=%0.3f,Xright=%0.3f\n",
			     m,S,Kt,Ht,Y[L],Y[R],JXmin[L],JXmax[L],JXmin[R],JXmax[R],Xleft,Xright);
		      fflush(stdout);
		      assert(Xleft <= Xright);
		    }
		  }

		  /* convert Xleft,Xright to jmin,jmax values */
		  if(Xleft > Xm[M+1] || Xright < 0.0){/* no alignment from Y[I] */
		    jmin = M+1;
		    jmax = 0;
		  } else {
		    for(jmin = M; jmin > 1; jmin --)
		      if(Xm[jmin] <= Xleft)
			break;
		    for(jmax = 1; jmax < M; jmax++)
		      if(Xm[jmax] >= Xright)
			break;
		    if((DEBUG>=2 && !(jmin <= jmax)) ||
		       (MDEBUG_SETLIMIT>=2 && MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT)){
		      printf("m=%d:S=%d,Kt=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f:Xleft=%0.3f,Xright=%0.3f,jmin=%d,jmax=%d,X[jmin]=%0.3f,X[jmax]=%0.3f,M=%d,X[M+1]=%0.3f\n",
			     m,S,Kt,Ht,Y[S],Y[S+1],Xleft,Xright,jmin,jmax,Xm[jmin],Xm[jmax],M,Xm[M+1]);
		      fflush(stdout);
		      assert(jmin <= jmax);
		    }
		  }
		}
		if(DEBUG && !(jmin <= jmax)){
		  printf("m=%d,S=%d,Kt=%d,I1=%d,I2=%d:jmin=%d,jmax=%d:IMIN1=%d,IMAX1=%d,IMIN=%d,IMAX=%d\n",
			 m,S,Kt,I1,I2,jmin,jmax,IMIN1,IMAX1,IMIN,IMAX);
		  fflush(stdout);
		  assert(jmin <= jmax);
		}
	      }
	    }

	    if(MDEBUG_SETLIMIT>=2 && MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT){
	      printf("Looking for updates to LRaddT[%d][%d]:S=%d,t=%d(T=%d):Ht=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f,Kt=%d,I1=%d,I2=%d,IMIN1=%d,IMAX1=%d,IMIN=%d,IMAX=%d,jmin=%d,jmax=%d,X[jmin]=%0.3f,X[jmax]=%0.3f\n",
		     MDEBUG_TS,MDEBUG_TT,S,t,addcntIcum[S]+t,Ht,Y[S],Y[S+1],Kt,I1,I2,IMIN1,IMAX1,IMIN,IMAX,jmin,jmax,Xm[jmin],Xm[jmax]); 
	      printf("    X[Jmin[IMIN1]]=%0.3f,X[Jmax[IMIN1]]=%0.3f,X[Jmin[IMAX1]]=%0.3f,X[Jmax[IMAX1]]=%0.3f,X[M+1]=%0.3f\n",
		     Xm[Jmin[IMIN1]],Xm[Jmax[IMIN1]],Xm[Jmin[IMAX1]],Xm[Jmax[IMAX1]],Xm[M+1]);
	      fflush(stdout);
	    }
	  
	    /* NOTE 1:	if (jmax < jmin) Line 5 of Equation 10 vanish for current S, Ht
	       NOTE 2 : if (S < IMIN-1) Line 5 of Equation 10 also vanishes, but line 6 (below) may still be valid */
	    if(jmax >= jmin && S >= IMIN-1){
	      if(DEBUG>=2) assert(0 < jmin && jmax <= M);

	      /* compute ALT[K=0..Kt][J=jmin..jmax] == A[S+1][K][J].AL AND PMT[K] == PM(S+1,K) with Y[S+1] = Ht (Line 5 of Equation 10) */
	      int I = S+1;
	      if(DEBUG>=2) assert(I >= Kt+1);
	      Y[I] = Ht;

	      /* NOTE: AL(),AR() will NOT use PRtab[].Pr */

	      if(DEBUG>=2) assert(I <= N+1);
	      if(I > IMAX+1 || I < IMIN-1)
		LogScaleL[I] = LogScaleL[min(max(I,IMIN-1),IMAX+1)];
	      double Scale = exp(-LogScaleL[I]);

	      /* precompute ScaleLI[P=I-1 ..max(IMIN, I-Kt-deltaExtYRef)]  = exp(LogScaleL[min(max(P,IMIN-1),IMAX+1)] - LogScaleL[I]) */
	      int Pmin = max(IMIN, I-Kt-deltaExtYRef);
	      double LogScaleLI = LogScaleL[I];
	      for(int P = Pmin; P < I; P++){
		if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
		ScaleLI[P] = exp(LogScaleL[P] - LogScaleLI);
		if(DEBUG>=2) assert(isfinite(ScaleLI[P]));
	      }
	
	      for(int K = 0; K <= Kt; K++){
		if(S-K > IMAX)
		  continue;// S can be greater than IMAX to handle all cases with S-K <= IMAX

		double Yik = Yc(Y,I,K);
		PMT[K] = PM(0,S+1,K,Y);
		if(DEBUG>=1+RELEASE/* HERE HERE >=2 */ && !(isfinite(PMT[K]) && PMT[K] >= 0.0)){
		  #pragma omp critical
		  {
		    printf("m=%d:S=%d,t=%d(T=%d):Kt=%d,Ht=%0.3f, Y[S]= %0.8f, Y[S+1]= %0.8f\n",m,S,t,addcntIcum[S]+t,Kt,Ht,Y[S],origY[I]);
		    printf("\t I=%d,K=%d,N=%d,IMIN=%d,IMAX=%d,Smin=%d,Smax=%d\n",I,K,N,IMIN,IMAX,Smin,Smax);
		    printf("\t Y[S+1-K]= %0.8f, PMT[K]= %0.10e\n",Y[S+1-K], PMT[K]);
		    fflush(stdout);

		    double PMTk = PM(0,S+1,K,Y);
		    double PMT2k = PM2(0,S+1,K,Y);

		    printf("\t PM(0,S+1,K,Y)= %0.10e, PM2(0,S+1,K,Y)= %0.10e\n", PMTk,PMT2k);
		    fflush(stdout);

		    assert(isfinite(PMT[K]) && PMT[K] >= 0.0);
		  }
		}

		if(DEBUG>=2 && !(I-K-1 <= IMAX)){
		  #pragma omp critical
		  {
		    printf("m=%d:S=%d,t=%d(T=%d):Kt=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f\n",m,S,t,addcntIcum[S]+t,Kt,Ht,Y[S],origY[I]);
		    printf("\t I=%d,K=%d,IMIN=%d,IMAX=%d,Smin=%d,Smax=%d\n",I,K,IMIN,IMAX,Smin,Smax);
		    fflush(stdout);
		    assert(I-K-1 <= IMAX);
		  }
		}
		
		for(int J = jmin; J <= jmax; J++)
		  ALT[K][J] = AL(m,I,K,J,0,-2,Yik,Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax,OUT_TYPE,OUT_OLD);
		/* ALT[K][J] does not include PMik : matches AL[S+1][K][J] on line 5 Equation 10 */
	      }
	      Y[I] = origY[I];

	      if(DEBUG>=2)
		for(int I = 0; I <= N+1; I++)
		  assert(fabs(Y[I]-origY[I]) <= 2e-9);

	      /* compute ART[K=0..Kt][J=jmin..jmax] == A[S][K][J].AR with Y[r-1] = origY[r] (for r=S+1-K ..S), Y[S] = Ht */
	      if(DEBUG>=2) assert(S-Kt >= 0);
	      Y[I=S] = Ht;

	      if(DEBUG>=2) assert(0 <= I  && I <= N+1);
	      if(I < IMIN-1 || I > IMAX+1)
		LogScaleR[I] = LogScaleR[min(max(I,IMIN-1),IMAX+1)];
	      Scale = exp(-LogScaleR[I]);

	      /* precompute ScaleRI[P= I+1 .. min(IMAX, I + KMAX + deltaExtYRef)] = exp(LogScaleR[P] - LogScaleR[I]) */
	      int Pmax = min(IMAX, I + KMAX + deltaExtYRef);
	      double LogScaleRI = LogScaleR[I];
	      for(int P = I+1; P <= Pmax; P++ ){
		if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
		ScaleRI[P] = exp(LogScaleR[P] - LogScaleRI);
		if(DEBUG>=2) assert(isfinite(ScaleRI[P]));		
	      }

	      if(DEBUG>=2) assert(isfinite(Scale));
	      for(int K = 0; K <= Kt; K++){
		if(K)/* modify Y[] incrementally for each K : no further updates to PRtab[][].Pr needed in AR() */
		  Y[S-K] = origY[S+1-K];

		if(S-K > IMAX)
		  continue;// S can be greater than IMAX to handle all cases with S-K <= IMAX

		if(DEBUG>=2 && !(I+1 >= IMIN)){
		  #pragma omp critical
		  {
		    printf("m=%d:S=%d,I=%d,K=%d,J=%d..%d:IMIN=%d,IMAX=%d,N=%d,Smin=%d,Smax=%d\n",m,S,I,K,jmin,jmax,IMIN,IMAX,N,Smin,Smax);
		    fflush(stdout);
		    assert(IMIN <= I+1);
		  }
		}

		double Yik = Yc(Y,I,K);
		for(int J = jmin; J <= jmax; J++)
		  ART[K][J] = AR(m,I,K,J,0,-2,Yik,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax,OUT_TYPE,OUT_OLD);
		/* No PMik : matches AR[S][K][J] on line 5 Equation 10 */
	      }
	      /* restore Y[S-Kt .. S] */
	      for(int K = 0; K <= Kt; K++)
		Y[S-K] = origY[S-K];

	      if(DEBUG>=2)
		for(int I = 0; I <= N+1; I++)
		  assert(fabs(Y[I]-origY[I]) <= 2e-9);

	      /* add to LRaddT[S][t] : sum(J=jmin..jmax,K=0..Kt) ALT[K][J]*PMT[K]*ART[K][J] */
	      for(int K = 0; K <= Kt; K++){
		if(S-K > IMAX)
		  continue; // S can be greater than IMAX to handle all  cases with S-K <= IMAX 
		double LR = 0.0;
		for(int J = jmin; J <= jmax; J++)
		  LR += ALT[K][J] * ART[K][J];
		double LR2 = LR * PMT[K] * ScaleS;
		if(DEBUG>=1+RELEASE/* HERE HERE >=2 */ && !(isfinite(LR2) && LR2>=0.0)){
		  printf("m=%d:S=%d,t=%d(T=%d):I=%d,K=%d,J=%d..%d:LRaddT[S][t]=%0.10e, LR2=%0.10e(LR=%0.10e,PMT[K]=%0.10e,ScaleS=%0.10e)\n",
			 m,S,t,addcntIcum[S]+t,S,K,jmin,jmax,(double)LRaddT[S][t],(double)LR2,(double)LR,(double)PMT[K],(double)ScaleS);
		  double sum = 0.0;
		  for(int J = jmin; J <= jmax; J++){
		    sum += ALT[K][J] * ART[K][J];
		    printf("  J=%d:ALT[K][J]=%0.10e,ART[K][J]=%0.10e, sum(ALT*ART)=%0.10e\n",
			   J,ALT[K][J],ART[K][J],(double)sum);
		  }
		  fflush(stdout);
		  assert(isfinite(LR2) && LR2>=0.0);
		}

		LRaddT[S][t] += LR2;
		if(DEBUG>=1+RELEASE /* HERE HERE >=2 */ && !(isfinite(LRaddT[S][t]) && LRaddT[S][t] >= 0.0)){
		  #pragma omp critical
		  {
		    printf("m=%d:S=%d,t=%d(T=%d):I=%d,K=%d,J=%d..%d:LRaddT[S][t]->%0.10e,delta(LR2)=%0.10e:LR=%0.10e,PMT[K]=%0.10e,ScaleS=%0.10e\n",
			   m,S,t,addcntIcum[S]+t,S,K,jmin,jmax,(double)LRaddT[S][t],(double)LR2,(double)LR,(double)PMT[K],(double)ScaleS);
		    double sum = 0.0;
		    for(int J = jmin; J <= jmax; J++){
		      sum += ALT[K][J] * ART[K][J];
		      printf("  J=%d:ALT[K][J]=%0.10e,ART[K][J]=%0.10e, sum(ALT*ART)=%0.10e\n",
			     J,ALT[K][J],ART[K][J],(double)sum);
		    }
		    fflush(stdout);
		    assert(isfinite(LRaddT[S][t]) && LRaddT[S][t] >= 0.0);
		  }
		}
		if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT && LR != 0.0){
		  FL128 sum = LRadd[0],sum2 = LRadd2[0];
		  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
		    sum += LRadd[i];
		    sum2 += LRadd2[i];
		  }
		  FL128 sum3 = sum;
		  if(Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
		    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
		    printf("m=%d:S=%d,t=%d(T=%d):I=%d,K=%d,J=%d..%d:sum(ALT[K][J]*PMT[K]*ART[K][J])=%0.10e(ScaleS=%0.10e):LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
			   m,S,t,addcntIcum[S]+t,I,K,jmin,jmax,(double)LR,(double)ScaleS,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
		    double LR3 = 0.0;
		    for(int J = jmin; J <= jmax; J++){
		      LR3 += ALT[K][J] * ART[K][J] * PMT[K] * ScaleS;
		      printf("    J=%d:ALT[K][J]*PM=%0.10e(LS=%0.10f),ART[K][J]*PM=%0.10e,PMT[K]=%0.10f:sum(ALT*ART*PMT)=%0.10e(delta=%0.10e)\n",
			     J,(double)(ALT[K][J]*PMT[K]),LogScaleL[S+1]+LogScaleR[S]-LS, (double)(ART[K][J]*PMT[K]), (double)PMT[K], (double)LR3, (double)(ALT[K][J] * ART[K][J] * PMT[K] * ScaleS));
		    }
		  } else
		    printf("m=%d:S=%d,t=%d(T=%d):I=%d,K=%d,J=%d..%d:sum(ALT[K][J]*PMT[K]*ART[K][J])=%0.10e(ScaleS=%0.10e):LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
			   m,S,t,addcntIcum[S]+t,I,K,jmin,jmax,(double)LR,(double)ScaleS,(double)sum,(double)sum2,(double)(sum-sum2));
		  fflush(stdout);
		}
	      }
	    }/* if(jmax >= jmin) */

	    /* Now work on line 6 Equation 10 */

	    /* compute largest K so that Ht and Y[S+K] may fail to resolve */
	    if(DEBUG) assert(SCORE_APPROX);
	    /* determine Kt based on total distance Y[S+Kt] - Ht < maxdist (if KMAX_APPROX) */
	    int Ks = 0;
	    while (S + Ks < N) {
	      if((KMAX_APPROX || (!KMAX_FIX && Ks==0)) ? (Y[S+Ks+1] - Ht >= maxdist) : ((!KMAX_FIX && Y[S+Ks+1] - Y[S+Ks] >= maxdist) || Y[S+Ks+1] - Ht >= AlignResMax))
		break;
	      Ks++;
	    }
	  
	    if(Ks > KMAX){
	      if(VERB>=2){
		printf("WARNING:m=%d,S=%d,t=%d/%d:Ht=%0.3f,Ks=%d reduced to KMAX=%d (N=%d)\n",m,S,t,addcntI[S],Ht,Ks,KMAX,N);
		fflush(stdout);
	      }
	      Ks = KMAX;
	    }
	    if(S+Ks > IMAX)
	      Ks = IMAX-S;

	    if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT){
	      printf("m=%d:S=%d,t=%d(T=%d):Ks=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+Ks]=%0.3f\n",m,S,t,addcntIcum[S]+t,Ks,Ht,Y[S],Y[S+Ks]);
	      fflush(stdout);
	    }

	    if(Ks <= 0)
	      continue;/* Line 6 of Equation 10 vanishes for current S, Ht */
	  
	    /* compute ALT[K=1..Ks][J=Jmin[S+K]..Jmax[S+K]] == A[S+K+1][K][J].AL with Y[r+1] = origY[r], r = S+K ... S+1, Y[S+1] = Ht (Line 6 of Equation 10) */
	    int I = S+1;
	    Y[I] = Ht;

	    /* NOTE: If _OPENMP, AL(),AR() will NOT use PRtab[].Pr */
	    for(int K = 1; K <= Ks; K++){
	      I = S+K+1;
	      Y[I] = origY[S+K]; /* modify Y[] incrementally for each K : updates to PRtab[][].Pr not needed in AL() */

	      if(I <= IMIN) /* S can be less than IMIN-1 to handle all cases with S+K >= IMIN */
		continue;
	  
	      double Scale = exp(-LogScaleL[I]);
	      double Yik = Yc(Y,I,K);
	      if(DEBUG>=2 && !(IMIN <= S+K && S+K <= IMAX)){
		printf("m=%d:S=%d,t=%d(T=%d):Ks=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+Ks]=%0.3f\n",m,S,t,addcntIcum[S]+t,Ks,Ht,Y[S],Y[S+Ks]);
		printf("\t I=%d,K=%d,S+K=%d,IMIN=%d,IMAX=%d,Smin=%d,Smax=%d\n",I,K,S+K,IMIN,IMAX,Smin,Smax);
		fflush(stdout);
		assert(IMIN <= S+K && S+K <= IMAX);
	      }
	      if(DEBUG>=2) assert(0 < Jmin[S+K] && Jmax[S+K] <= M);
	      if(DEBUG>=2 && !(I-K-1 <= IMAX)){
		printf("m=%d:S=%d,t=%d(T=%d):Ks=%d,Ht=%0.3f,Y[S]=%0.3f,Y[S+Ks]=%0.3f\n",m,S,t,addcntIcum[S]+t,Ks,Ht,Y[S],Y[S+Ks]);
		printf("\t I=%d,K=%d,S+K=%d,IMIN=%d,IMAX=%d,Smin=%d,Smax=%d\n",I,K,S+K,IMIN,IMAX,Smin,Smax);
		fflush(stdout);
		assert(I-K-1 <= IMAX);
	      }

	      /* precompute ScaleLI[P=I-K-1 ..max(IMIN, I-Kt-deltaExtYRef)]  = exp(LogScaleL[min(max(P,IMIN-1),IMAX+1)] - LogScaleL[I]) */
	      int Pmin = max(IMIN, I-K-deltaExtYRef);
	      double LogScaleLI = LogScaleL[I];
	      for(int P = Pmin; P < I-K; P++){
		if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
		ScaleLI[P] = exp(LogScaleL[P] - LogScaleLI);
		if(DEBUG>=2) assert(isfinite(ScaleLI[P]));
	      }

	      for(int J = Jmin[S+K]; J <= Jmax[S+K]; J++)
		ALT[K][J] = AL(m,I,K,J,0,-2,Yik,Xm,Y,M,N,lc,Scale,ScaleLI,LogScaleL,A,IMIN,IMAX,Kmax,Jmin,Jmax,OUT_TYPE,OUT_OLD);
	      /* No PMik : matches AL[S+K+1][K][J] on line 6 Equation 10 */
	    } /* K = 1 .. Ks */

	    /* restore Y[] */
	    for(int r = S+Ks+1; r > S; r--)
	      Y[r] = origY[r];

	    I = S+1;

	    if(DEBUG>=2)
	      for(int I = 0; I <= N+1; I++)
		assert(fabs(Y[I]-origY[I]) <= 2e-9);

	    /* compute ART[K=1..Ks][J=Jmin[S+K]..Jmax[S+K]] == A[S+K][K][J].AR AND PMT[K] == PM(S+K,K) with Y[S] = Ht (Line 6 of Equation 10) */	
	    Y[I=S] = Ht;
	    for(int K = 1; K <= Ks; K++){
	      I = S+K;
	      if(I < IMIN) /* S can be less than IMIN-1 to handle all cases with S+K >= IMIN */
		continue;

	      /* precompute ScaleRI[P= I+1 .. min(IMAX, I + KMAX + deltaExtYRef)] = exp(LogScaleR[P] - LogScaleR[I]) */
	      int Pmax = min(IMAX, I + KMAX + deltaExtYRef);
	      double LogScaleRI = LogScaleR[I];
	      for(int P = I+1; P <= Pmax; P++){
		if(DEBUG>=2) assert(IMIN-1 <= P && P <= IMAX+1);
		ScaleRI[P] = exp(LogScaleR[P] - LogScaleRI);
		if(DEBUG>=2) assert(isfinite(ScaleRI[P]));		
	      }

	      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX);
	      double Scale = exp(-LogScaleR[I]);	  
	      double Yik = Yc(Y,I,K);
	      for(int J = Jmin[S+K]; J <= Jmax[S+K]; J++)
		ART[K][J] = AR(m,I,K,J,0,-2,Yik,Xm,Y,M,N,rc,Scale,ScaleRI,LogScaleR,A,IMIN,IMAX,Kmax,Jmin,Jmax,OUT_TYPE,OUT_OLD);
	      /* No PMik : matches AR[S+K][K][J] on line 6 Equation 10 */

	      PMT[K] = PM(0,S+K,K,Y);
	    }/* K = 1 .. Ks */

	    /* restore Y[] */
	    Y[I=S] = origY[S];

	    if(DEBUG>=2)
	      for(int I = 0; I <= N+1; I++)
		assert(fabs(Y[I]-origY[I]) <= 2e-9);

	    /* add to LRaddT[S][t] : sum(J=jmin..jmax,K=1..Ks) ALT[K][J]*PMT[K]*ART[K][J] */
	    for(int K = 1; K <= Ks; K++){
	      if(S+K < IMIN) /* S can be less than IMIN-1 to handle all cases with S+K >= IMIN */
		continue;
	      double LR = 0.0;
	      for(int J = Jmin[S+K]; J <= Jmax[S+K]; J++)
		LR += ALT[K][J] * ART[K][J];
	      if(DEBUG>=2) assert(IMIN <= S+K && S+K <= IMAX);
	      LR *= PMT[K] * exp(LogScaleL[S+K+1]+LogScaleR[S+K]-LS);
	      if(DEBUG>=2) assert(isfinite(LR) && LR >= 0.0);
	      if(DEBUG>=2) assert(Amin <= S && S <= min(N,Amax));
	      LRaddT[S][t] += LR;
	      if(DEBUG>=1+RELEASE /* HERE HERE >= 2 */ && !(isfinite(LRaddT[S][t]) && LRaddT[S][t] >= 0.0)){
		#pragma omp critical
		{
		  printf("\n m=%d/%d:S=%d(Smin=%d,Smax=%d),t=%d (addcntI[S]=%d), K=%d,Ks=%d, IMIN=%d,IMAX=%d, Jmin[S+K]=%d,Jmax[S+K]=%d:LR= %0.8e, LRaddT[S][t]= %0.8e\n",
			 m,MD,S,Smin,Smax,t,addcntI[S],K,Ks,IMIN,IMAX,Jmin[S+K],Jmax[S+K],LR,LRaddT[S][t]);
		  double cumLR = 0.0;
		  for(int J = Jmin[S+K]; J <= Jmax[S+K]; J++){
		    cumLR += ALT[K][J] * ART[K][J];
		    printf("\t J=%d: ALT[K][J]= %0.8e, ART[K][J]= %0.8e : cum(ALT * ART)=  %0.8e\n",J,ALT[K][J],ART[K][J],cumLR);
		  }
		  printf("\t\t PMT[K]= %0.8e, LogScaleL[S+K+1]= %0.8f, LogScaleR[S+K]= %0.8f, LS= %0.8f, PMT[K] * exp(%0.8f) = %0.8e\n",
			 PMT[K],LogScaleL[S+K+1],LogScaleR[S+K],LS, LogScaleL[S+K+1]+LogScaleR[S+K]-LS, PMT[K] * exp(LogScaleL[S+K+1]+LogScaleR[S+K]-LS));
		  fflush(stdout);

		  assert(isfinite(LRaddT[S][t]) && LRaddT[S][t] >= 0.0);
		}
	      }
	      if(MDEBUG_T >= 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N && S==MDEBUG_TS && t==MDEBUG_TT && LR != 0.0){
		double ScaleS = exp(LogScaleL[S+K+1] + LogScaleR[S+K] - LS);
		FL128 sum = LRadd[0],sum2 = LRadd2[0];
		for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
		  sum += LRadd[i];
		  sum2 += LRadd2[i];
		}
		FL128 sum3 = sum;
		if(Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
		  sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
		  printf("m=%d:S=%d,t=%d(T=%d):I=%d,K=%d,J=%d..%d:sum(ALT[K][J]*PMT[K]*ART[K][J])=%0.10e(ScaleS=%0.10e):LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e)\n",
			 m,S,t,addcntIcum[S]+t,S+K,K,Jmin[S+K],Jmax[S+K],(double)LR,(double)ScaleS,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2));
		  double LR2 = 0.0;
		  for(int J = Jmin[S+K]; J <= Jmax[S+K]; J++){
		    LR2 += ALT[K][J] * ART[K][J] * PMT[K] * ScaleS;
		    printf("    J=%d:ALT[K][J]*PM=%0.10e(LS=%0.10f),ART[K][J]*PM=%0.10e,PMT[K]=%0.10e:sum(ALT*ART*PMT)=%0.10e(delta=%0.10e)\n",
			   J,(double)(ALT[K][J]*PMT[K]),LogScaleL[S+K+1]+LogScaleR[S+K]-LS, (double)(ART[K][J]*PMT[K]), (double)PMT[K], (double)LR2, (double)(ALT[K][J] * ART[K][J] * PMT[K] * ScaleS));
		  }
		} else
		  printf("m=%d:S=%d,t=%d(T=%d):I=%d,K=%d,J=%d..%d:sum(ALT[K][J]*PMT[K]*ART[K][J])=%0.10e(ScaleS=%0.10e):LRadd[0..S]=%0.10e-%0.10e=%0.10e\n",
			 m,S,t,addcntIcum[S]+t,S+K,K,jmin,jmax,(double)LR,(double)ScaleS,(double)sum,(double)sum2,(double)(sum-sum2));
		fflush(stdout);
	      }	    
	    }
	  } /* t = 0 .. addcntI[S] -1 */
	} /* S = IMIN-1 .. IMAX */
      } /* if(myaddcnt > 0) */

      if(VERB>=2 && newLPdelta && rverb){
	printf("mprobeval:m=%d,Amin=%d,Amax=%d:After PM boundary updates:time=%0.6f\n",m,Amin,Amax,wtime());
	fflush(stdout);
      }

      if(MDEBUG_M >= 0){
	if(MDEBUG_S >= 0 && m==MDEBUG_M && N==MDEBUG_N && LPdelG){
	  FL128 sum = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	  for(int i = max(1,Dmin); i <= min(MDEBUG_S,Dmax); i++)
	    sum += LRdel[i];
	  printf("m=%d:After PM boundary updates:LRdel[0..%d]=%0.10e(LS=%0.10f),IMIN=%d,IMAX=%d,Dmin=%d,Dmax=%d,N=%d\n",m,MDEBUG_S,(double)sum,LS,IMIN,IMAX,Dmin,Dmax,N);
	}
	if(LDEBUG && MDEBUG_F >= 0 && m==MDEBUG_M && N==MDEBUG_N-1){
	  FL128 sum = 0.0, sum2 = 0.0;
	  for(int i = 0; i <= MDEBUG_F; i++){
	    sum += LRfrag[i];
	    sum2 += LRfrag2[i];
	  }
	  printf("m=%d:After PM boundary updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f),IMIN=%d,IMAX=%d,N=%d\n",m,MDEBUG_F,(double)sum,(double)sum2,(double)(sum-sum2),LS,IMIN,IMAX,N);
	}
	if(MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  FL128 sum = LRadd[0], sum2 = LRadd2[0];
	  for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++){
	    sum += LRadd[i];
	    sum2 += LRadd2[i];
	  }
	  FL128 sum3 = sum;
	  if(myaddcnt > 0 && Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < addcntI[MDEBUG_TS]){
	    sum3 += LRaddT[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After PM boundary updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e,LRaddT[S][T]=%0.10e(sum=%0.10e,LS=%0.10f))\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),(double)LRaddT[MDEBUG_TS][MDEBUG_TT],(double)(sum3-sum2),LS);
	  } else
	    printf("m=%d:After PM boundary updates:S=%d,T=%d,LRadd[0..S]=%0.10e-%0.10e=%0.10e(LS=%0.10f),Amin=%d,Amax=%d,addcntI[%d]=%d\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum,(double)sum2,(double)(sum-sum2),LS,Amin,Amax,MDEBUG_TS,addcntI[MDEBUG_TS]);
	  if(LDEBUG){
	    sum = sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After PM boundary updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS);
	  }
	}
	if((MDEBUG_ST >= 0 && m==MDEBUG_M) || (MDEBUG_T >= 0 && MDEBUG_TS>=0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N+1)){
	  FL128 sum = 0.0;
	  for(int i = 1; i <= max(MDEBUG_TS+1,MDEBUG_ST); i++)
	    sum += LRsite[i]; 
	  printf("m=%d:After PM boundary updates:LRsite[0..%d]=%0.10e,LRcum=%0.10e(LS=%0.10f)\n",m,max(MDEBUG_TS+1,MDEBUG_ST),(double)sum,(double)LRcum,LS);
	}
	if(MDEBUG_T < 0 && MDEBUG_TS >= 0 && MDEBUG_TT >= 0 && m==MDEBUG_M && N==MDEBUG_N){
	  if(newLPdelta && MDEBUG_TT < D[min(N,MDEBUG_TS)]){
	    FL128 sum = LRdelta[0];
	    for(int i = max(1,Amin); i <= min(Amax,MDEBUG_TS); i++)
	      sum += LRdelta[i];
	    FL128 sum3 = sum;
	    if(Amin <= MDEBUG_TS && Amax >= MDEBUG_TS && MDEBUG_TT < D[MDEBUG_TS])
	      sum3 += LRdeltaD[MDEBUG_TS][MDEBUG_TT];
	    printf("m=%d:After PM boundary updates:S=%d,D=%d,LRdelta[0..S]=%0.10e,LRdeltaD[S][D]=%0.10e(sum=%0.10e,LS=%0.10f)(LRdelta[0]=%0.10e)\n",
		   m,MDEBUG_TS,MDEBUG_TT,(double)sum, (double)LRdeltaD[MDEBUG_TS][MDEBUG_TT], (double)sum3, LS, (double)LRdelta[0]);
	  }
	  if(LDEBUG){
	    FL128 sum = 0.0, sum2 = 0.0;
	    for(int i = 0;i <= MDEBUG_TS; i++){
	      sum += LRfrag[i];
	      sum2 += LRfrag2[i];
	    }
	    printf("m=%d:After PM boundary updates:LRfrag[0..%d]=%0.10e-%0.10e=%0.10e(LS=%0.10f)\n",m,MDEBUG_TS,(double)sum,(double)sum2,(double)(sum-sum2),LS);
	  }
	}
	fflush(stdout);
      }
      if(0 && MDEBUG && LDEBUG && !FIXADD){
	FL128 addsum = LRadd[0], addsum2 = LRadd2[0];
	for(int i = 1; i <= 52;i++){
	  addsum += LRadd[i];
	  addsum2 += LRadd2[i];
	}
	addsum -= addsum2;
	FL128 fragsum = LRfrag[0], fragsum2 = LRfrag2[0];
	for(int i = 1; i <= 52; i++){
	  fragsum += LRfrag[i];
	  fragsum2 += LRfrag2[i];
	}
	fragsum -= fragsum2;
	if(DEBUG && (!(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190) || (VERB>=2&&m==1))){
	  printf("m=%d:After PM boundary updates:sum(LRadd[0..52])=%0.17e,sum(LRfrag[0..52]=%0.17e\n", m,(double)addsum,(double)fragsum);
	  fflush(stdout);
	  assert(addsum <= fragsum*1.000001+1e-190 && addsum >= fragsum*FN[0]*0.999999-1e-190);
	}
      }
      if(VERB && (LVERB >=2 || (0 && MDEBUG && m==MDEBUG_M))){
	printf("m=%d:After PM (boundry) term updates:\n",m);
	FL128 LRsitecum = 0.0;
	for(int s = 1; s <= (LDEBUG ? N+1 : N); s++){
	  LRsitecum += LRsite[s];
	  printf("m=%d,s=%d:LRsite[s]=%0.10e,cum=%0.10e,LRcum=%0.10e\n",m,s,(double)LRsite[s],(double)LRsitecum,(double)LRcum);
	}
	fflush(stdout);
      }

      if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	#pragma omp critical
	{
	  printf("\t\t m=%d:tid=%d:After PM boundary updates:wall time= %0.6f\n", m, tid, wtime() - wt);
	  fflush(stdout);
        }
      }

      if(DEBUG) assert(isfinite(LRcum));
      if(DEBUG) assert(LS >= 0.0);
      double LRbiasLS = LRbias * exp(-LS);
      double LPm = (startLRcum > 0.0 || LRbiasLS > 0.0) ? LS+log(startLRcum + LRbiasLS) : -99999.0;
      if(DEBUG) assert(isfinite(LPm));
      if(newLP) newLP[m] = LPm;
      if(VERB>=2 && EPS2 > 0.0 && m==MDEBUG_M && rverb){
        #pragma omp critical
	{
	  printf("m=%d: EPS2 adjusted LPm= %0.10f\n",m, LPm);
	  fflush(stdout);
	}
      }

      if(LDEBUG && (MDEBUG_M < 0 || m==MDEBUG_M)){/* check that cumulative value of LRsite[1..N] matches LRcum and LRsite[1..N+1] approximates 0 */
	if(VERB>=2){
	  printf("m=%d:IMIN=%d,IMAX=%d,N=%d\n",m,IMIN,IMAX,N);
	  fflush(stdout);
	}
	FL128 LRsitecum = 0.0;
	int S = 1;
	for(; S <= N+1; S++){
	  LRsitecum += LRsite[S];
	  if(VERB && (VERB>=2 || (DEBUG && ((S <= N) ? !(fabs((double)(LRsitecum)-LRcum) <= (LRcum+LRbiasLS) * LRSITE_ERR) : !(fabs((double)LRsitecum) <= (LRcum+LRbiasLS) * LRSITE_ERR))))){
	    printf("m=%d:S=%d,N=%d:LRsite[S]=%0.10e,LRsitecum=%0.10e,LRcum=%0.10e,LRfrag[0]=%0.10e,LRsitecum-LRcum=%0.10e\n",
		   m,S,N,(double)LRsite[S],(double)LRsitecum,LRcum,(double)LRfrag[0],(double)(LRsitecum-LRcum));
	    FL128 LRsitecum2 = 0.0;
	    for(int s = 1; s <= N+1; s++){
	      LRsitecum2 += LRsite[s];
	      printf("  LRsite[%d]=%0.10e,LRfrag[%d]=%0.10e:cum=%0.10e,cum-LRcum=%0.10e\n",s,(double)LRsite[s],s,(double)LRfrag[s],(double)LRsitecum2,(double)(LRsitecum2-LRcum));
	    }
	    printf("m=%d:S=%d,N=%d:LRsite[S]=%0.10e,LRsitecum=%0.10e,LRcum=%0.10e,LRbiasLS=%0.10e,LRfrag[0]=%0.10e,LRsitecum-LRcum=%0.10e,lc=%d,rc=%d\n",
		   m,S,N,(double)LRsite[S],(double)LRsitecum,LRcum,LRbiasLS,(double)LRfrag[0],(double)(LRsitecum-LRcum),lc,rc);
	    fflush(stdout);
	    if(DEBUG && S <= N) assert(fabs((double)LRsitecum-LRcum) <= (LRcum+LRbiasLS) * LRSITE_ERR);
	    if(DEBUG && S == N+1) assert(fabs((double)LRsitecum) <= (LRcum+LRbiasLS) * LRSITE_ERR);
	  }
	}
      }

      if(LDEBUG>=2){/* check that cumulative value of LRfrag[0..N+1] matches LRcum */
	FL128 LRfragcum = 0.0, LRfragcum2 = 0.0;;
	for(int S = 0; S <= N+1; S++){
	  LRfragcum += LRfrag[S];
	  LRfragcum2 += LRfrag2[S];
	  if(VERB && (VERB>=2 || (DEBUG && (S <= N ? !(fabs((double)(LRfragcum-LRfragcum2)-LRcum) <= (LRcum+LRbiasLS) * LRFRAG_ERR) : !(fabs((double)(LRfragcum-LRfragcum2)) <= (LRcum+LRbiasLS) * LRFRAG_ERR ))))){
	    printf("m=%d:S=%d,N=%d:LRfrag[S]=%0.10e,LRfrag2[S]=%0.10e,diff=%0.10e(cum=%0.10e,%0.10e,diff= %0.10e),LRcum=%0.10e,LRbiasLS=%0.10e\n",
		   m,S,N,(double)LRfrag[S],(double)LRfrag2[S],(double)(LRfrag[S]-LRfrag2[S]),(double)LRfragcum,(double)LRfragcum2,(double)(LRfragcum-LRfragcum2),LRcum,LRbiasLS);
	    FL128 cum = 0.0, cum2 = 0.0;
	    for(int s = 0; s <= N+1; s++){
	      cum += LRfrag[s];
	      cum2 += LRfrag2[s];
	      printf("  s=%d:LRfrag[s]=%0.10e,LRfrag2[s]=%0.10e,cum=%0.10e,cum2=%0.10e,diff=%0.10e:cum-cum2-LRcum=%0.10e\n",
		     s,(double)LRfrag[s],(double)LRfrag2[s],(double)(LRfrag[s]-LRfrag2[s]),(double)cum,(double)cum2,(double)(cum-cum2-LRcum));
	    }
	    printf("m=%d:S=%d,N=%d:LRfrag[S]=%0.10e,LRfrag2[S]=%0.10e,diff=%0.10e(cum=%0.10e,%0.10e,diff= %0.10e),LRcum=%0.10e,LRbiasLS=%0.10e\n",
		   m,S,N,(double)LRfrag[S],(double)LRfrag2[S],(double)(LRfrag[S]-LRfrag2[S]),(double)LRfragcum,(double)LRfragcum2,(double)(LRfragcum-LRfragcum2),LRcum,LRbiasLS);
	    fflush(stdout);
	    if(DEBUG && S <= N) assert(fabs((double)(LRfragcum-LRfragcum2)-LRcum) <= (LRcum + LRbiasLS) * LRFRAG_ERR);
	    if(DEBUG && S == N+1) assert(fabs((double)(LRfragcum-LRfragcum2)) <= (LRcum + LRbiasLS) * LRFRAG_ERR);
	  }
	}
      }

      if(TDEBUG){
        #pragma omp critical
        {
	  if(VERB>=2 && m==1){
	    printf("m=%d:LPm=%0.10f,TBwtm=%0.10f:Amin=%d,Amax=%d,N=%d\n",m,LPm,(double)TBwtm,Amin,Amax,N);
	    fflush(stdout);
	  }
	  LLRsum += LPm * TBwtm;
	  if(VERB && (VERB>=2 || rverb || (MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M))){
	    printf("m=%d:LPm=%0.10f:LRcum=%0.10e,LRbias=%0.8e,LS=%0.8f(sum=%0.10f)\n",m,LPm,LRcum,LRbias,LS,LLRsum);
	    fflush(stdout);
	  }
	  if(VERB && (VERB >=2 || rverb)){
	    if(newLPdelta)
	      printf("m=%d/%d:mapid=%d(id=%lld):LLR=%0.6f(delta=%0.6f,LS=%0.6e,LRcum=%0.6e,LRbias=%0.6e),M=%d(%d),N=%d,Amin=%d,Amax=%d,time=%0.6f\n",
		     m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,LLRsum,LPm,LS,LRcum,LRbias,M,MX[m],N,Amin,Amax,wtime());
	    else
	      printf("m=%d/%d:mapid=%d(id=%lld):LLR=%0.6f(delta=%0.6f,LS=%0.6e,LRcum=%0.6e,LRbias=%0.6e),M=%d(%d),N=%d\n",
		     m,MD,pcontig->contig[m].mapid,gmap[pcontig->contig[m].mapid]->id,LLRsum,LPm,LS,LRcum,LRbias,M,MX[m],N);
	    fflush(stdout);
	  }
	  if(DEBUG) assert(isfinite(LLRsum));
	} // critical
      }
      
      FL128 LR, LR2;
      double LP;

      if(LPdelG){ /* update differential array LPdel[1..N] from LRdel[Dmin..Dmax] and LRdel[0] 
		     LPdel[0] += log(LRdel[0]+LRbias),  LPdel[S=Dmin..Dmax] += log(LRdel[0] + sum(LRdel[Dmin..S]) + LRbias) - log(LRdel[0] + sum(LRdel[Dmin..S-1]) + LRbias) */
	FL128b LR = (LDEBUG || FIXLRDEL) ? LRdel[0] : LRcum;
	FL128b LRd = LR * LRratio;/* to avoid log(0.0) when LR < 1.0e-306 */
	LP = (LRd > 0.0 || LRbiasLS > 0.0) ? LS+logl(max((FL128b)0.0,LRd) + LRbiasLS) : -99999.0;// WAS LPm;

	if(TDEBUG || MDEBUG || !newLPd){
          #pragma omp critical
	  {
	    LPdel[0] += LP * TBwtm;
	    if(DEBUG && !LDEBUG && Dmin <= 0) assert(LRdel[0] == 0.0);
	  }
	}
	if(newLPd){
	  if(MDEBUG) newLPdPr[m][0] = 0;
	  DminM[m] = max(1,Dmin);
	  DmaxM[m] = min(N,Dmax);
	  if(VMEM_MINIMIZE){
	    newLPd[m] = thread_alloc(DmaxM[m] - DminM[m] + 1, doubleAlloc, origY[IMIN], ppBlock, tid);
	    newLPd[m] -= DminM[m];
	    newLPd0[m] = LP;
	  } else
	    newLPd[m][0] = LP;
	}
	for(int S = max(1,Dmin); S <= Dmax; S++){
	  if(DEBUG && !(isfinite(LRdel[S]))){
            #pragma omp critical
	    {
	      printf("\nm=%d/%d:S=%d(Dmin=%d,Dmax=%d),N=%d:LRdel[S]=%Le,LRdel[0]=%Le,IMIN=%d,IMAX=%d\n",m,MD,S,Dmin,Dmax,N,(long double)LRdel[S],(long double)LRdel[0],IMIN,IMAX);
	      fflush(stdout);
	      assert(isfinite(LRdel[S]));
	    }
	  }
	  LR += LRdel[S];
	  if(S >= N+1)
	    break;
	  FL128b LRd = LR * LRratio;/* to avoid log(0.0) when LR < 1.0e-306 */
	  double LPd = (LRd > 0.0 || LRbiasLS > 0.0) ? LS+logl(max((FL128b)0.0,LRd) + LRbiasLS) : -99999.0;
	  if(DEBUG && !isfinite(LPd)){
	    printf("m=%d:S=%d,N=%d:LPdel[S] += %13.10f (TBwtm=%0.10f,LP=%0.10f,LPd=%0.10f,LR=%0.24Le,LRbias=%0.8e,LS=%0.8f,LRbiasLS=%0.8e),LRcum=%0.10e,LRdel[S]=%0.24Le,LRdel[0]=%0.24Le\n",
		   m,S,N, (double)(TBwtm * (LPd - LP)), (double)TBwtm, LP,LPd,(long double)LRd,LRbias,LS,LRbiasLS,LRcum,(long double)LRdel[S],(long double)LRdel[0]);
	    assert(isfinite(LPd));
	  }
	  if(MDEBUG && MDEBUG_M >= 0 && MDEBUG_S >= 0 && m==MDEBUG_M && S==MDEBUG_S && N==MDEBUG_N){
	    printf("m=%d:S=%d,N=%d:LPdel[S] += %13.10f (LR=%0.24Le,LRbias=%0.8e,LS=%0.8f),LRcum=%0.10e,LRdel[S]=%0.24Le,LRdel[0]=%0.24Le\n",
		   m,S,N,LPd,(long double)LRd,LRbias,LS,LRcum,(long double)LRdel[S],(long double)LRdel[0]);
	    fflush(stdout);
	  }
	  if(TDEBUG || MDEBUG || !newLPd){
            #pragma omp critical
	    {
	      LPdel[S] += TBwtm * (LPd - LP);
	    }
	  }
	  LP = LPd;
	  if(newLPd){
	    if(MDEBUG) newLPdPr[m][S] = (LR < LRcum * ((sizeof(FL128)==8) ? 1e-7 : 1e-11)) ? 1 : 0;/* excessive loss of precision */
	    if(VMEM_MINIMIZE && DEBUG>=2) assert(DminM[m] <= S && S <= DmaxM[m]);
	    newLPd[m][S] = LPd;
	  }
	}

	if(newLPd)
	  if(DEBUG) assert(fabs(LP - newLPd[m][DmaxM[m]]) < 1e-6);

	if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	  #pragma omp critical
	  {
	    printf("\t\t m=%d:tid=%d:After LPdelG[] updates:wall time= %0.6f\n", m, tid, wtime() - wt);
	    fflush(stdout);
	  }
	}

      } // if(LPdelG)

      size_t LPaddSiz = 0;

      if(myaddcnt){
	/* update LPaddI and differential array LPadd[t=0..addcnt-1] from LRadd[] and LRaddT[][] of current map m */
	/* NOTE : LRadd[0]-LRadd2[0] is approximately the same as LRcum (but computed with DELTA_Y reduced by 1) */
	/* LPaddI += log(LRadd[0]+LRbias), 
	   LPadd[t=0..addcnt-1] += log(LRadd[0] + sum(LRadd[max(1,Amin)..S(t)]   + LRaddT[S(t)][T(t)] + LRbias) - log(LRadd[0] + sum(LRadd[max(1,Amin)..S(t-1)] + sum(LRaddT[S(t-1)][T(t-1)]) + LRbias) */
	if(DEBUG) assert(isfinite(LRadd[0]) && LRadd[0] >= 0);
	if(DEBUG && !(isfinite(LRadd2[0]) && LRadd2[0] == 0)){
	  printf("m=%d:LRadd[0]=%0.10e,LRadd2[0]=%0.10e\n",m,(double)LRadd[0],(double)LRadd2[0]);
	  fflush(stdout);
	  assert(isfinite(LRadd2[0]) && LRadd2[0] == 0);
	}
	LR = LRadd[0]-LRadd2[0];
	if(DEBUG)assert(isfinite(LR));
	FL128b LRd = LR * LRratio;/* to avoid log(0.0) when LR < 1.0e-306 */
	LP = (LRd > 0.0 || LRbiasLS > 0.0) ? LS+logl(max((FL128b)0.0,LRd) + LRbiasLS) : -99999.0;
	if(DEBUG) assert(isfinite(LP));
	if(TDEBUG || MDEBUG || !newLPa){
          #pragma omp critical
	  {
	    LPaddI += LP * TBwtm;
	    if(DEBUG) assert(isfinite(LPaddI));
	  }
	}

	if(newLPa){
	  TminM[m] = max(0, addcntIcum[Amin]-1);
	  TmaxM[m] = min(addcnt-1, addcntIcum[min(N,Amax)+1]);
	  if(DEBUG>=2 && !(TmaxM[m] >= TminM[m])){
	    #pragma omp critical
	    {
	      printf("m=%d:Amin=%d,Amax=%d,N=%d,addcnt=%d,addcntIcum[Amin]=%d,addcntIcum[min(N,Amax)+1]=%d,TminM[m]=%d,TmaxM[m]=%d\n",
		     m,Amin,Amax,N,addcnt,addcntIcum[Amin],addcntIcum[min(N,Amax)+1],TminM[m],TmaxM[m]);
	      fflush(stdout);
	      assert(TmaxM[m] >= TminM[m] && (TmaxM[m] > TminM[m] || TminM[m] > 0));// NEW6
	    }
	  }
	  if(VMEM_MINIMIZE){
	    LPaddSiz = TmaxM[m] - TminM[m] + 1;
	    if(DEBUG>=2) assert(LPaddSiz > 0);// NEW6
	    newLPa[m] = thread_alloc(LPaddSiz, doubleAlloc, origY[IMIN], ppBlock,tid);// threadmem, threadmemsiz);
	    newLPa[m] -= TminM[m];
	  }
	  if(VERB>=2 && rverb){
	    #pragma omp critical
	    {
	      printf("m=%d:IMIN=%d,IMAX=%d,N=%d,Amin=%d,Amax=%d,TminM[m]=%d(addloc=%0.3f),TmaxM[m]=%d(addloc=%0.3f),addcnt=%d,Y[N+1]=%0.3f\n",
		     m,IMIN,IMAX,N,Amin,Amax,TminM[m],addloc[TminM[m]],TmaxM[m],addloc[TmaxM[m]],addcnt,Y[N+1]);
	      fflush(stdout);
	    }
	  }
	  if(Amin > 0 && addcntIcum[Amin] > 0){
	    int t = addcntIcum[Amin]-1;
	    if(MDEBUG) newLPaPr[m][t] = (LR < LRcum*1e-7) ? 1 : 0;/* loss of precision */
	    if(DEBUG>=2 && !VMEM_MINIMIZE) assert(newLPa[m][t] <= LARGE_NEGATIVE*0.99);
	    if(DEBUG) assert(TminM[m] <= t && t <= TmaxM[m]);
	    newLPa[m][t] = LP;
#ifdef _OPENMP
	    if(MDEBUG && MDEBUG_M >= 0 && MDEBUG_T >= 0 && t==MDEBUG_T && m==MDEBUG_M && N==MDEBUG_N){
	      printf("m=%d:N=%d,t=%d,addcnt=%d:LPadd[t] += (newLPa[m][t]=%13.10f) (LR=%0.10e,LRbias=%0.8e,LS=%0.8f), Amin=%d,Amax=%d:addloc[t]=%0.3f\n",
		     m, N, t, addcnt, newLPa[m][t], (double)LR, LRbias, LS, Amin,Amax, addloc[t]);
	      fflush(stdout);
	    }
#else
	    if(MDEBUG || DEBUG>=2){ /* verify that LPaddI + sum(LPadd[0..t]) = sum(newLPa[i=0..m][t] * TBmapWT[i] ) */
	      FL128 sum1 = LPaddI;
	      for(int i = 0; i <= t; i++)
		sum1 += LPadd[i];
	      FL128 sum2 = 0.0;
	      if(TBmapWT){
		for(int i = 0; i <= m; i++){
		  int T = (t < TminM[i]) ? TminM[i] : (i < m && t > TmaxM[i]) ? TmaxM[i] : t;
		  if(DEBUG>=2 && !VMEM_MINIMIZE) assert(newLPa[i][T] > LARGE_NEGATIVE*0.99);
		  sum2 += newLPa[i][T] * TBmapWT[i];
		}
	      } else {
		for(int i = 0; i <= m; i++){
		  int T = (t < TminM[i]) ? TminM[i] : (i < m && t > TmaxM[i]) ? TmaxM[i] : t;
		  if(DEBUG>=2 && !VMEM_MINIMIZE) assert(newLPa[i][T] > LARGE_NEGATIVE*0.99);
		  sum2 += newLPa[i][T];
		}
	      }
	      if(DEBUG && (VERB>=2 || (MDEBUG_T >= 0 && t==MDEBUG_T && m==MDEBUG_M) || !(fabs(sum1-sum2) < max(1.0,fabs(sum1)) * CUMERR))){
		printf("m=%d:N=%d,t=%d,addcnt=%d:LPadd[0..t]=%0.10f,newLPa[0..m][t]=%0.10f (error = %0.14f), LPadd[t] += newLPa[m][t]=%13.10f (LR=%0.10e,LRbias=%0.8e,LS=%0.8f), addloc[t]=%0.3f\n",
		       m, N, t, addcnt, (double)sum1, (double)sum2, (double)(sum2-sum1), newLPa[m][t], (double)LR, LRbias, LS, addloc[t]);
		printf("  TBmapWT[m]=%0.10f,TBwtm=%0.10f,LP=%0.10f\n", TBmapWT[m],(double)TBwtm,LP);
		fflush(stdout);
		assert(fabs(sum1-sum2) < max(1.0,fabs(sum1)) * CUMERR);
	      }
	    }
#endif
	  }
	}// if(newLPa)

	for(int S = Amin; S <= min(N,Amax); S++){
	  LR = LRadd[0];
	  LR2 = LRadd2[0];
	  FL128 fragsum = (LDEBUG ? LRfrag[0] : 0.0), fragsum2 = (LDEBUG ? LRfrag2[0] : 0.0);
	  for(int s = max(1,Amin);s <= S; s++){
	    LR += LRadd[s];
	    LR2 += LRadd2[s];
	    if(LDEBUG){
	      fragsum += LRfrag[s];
	      fragsum2 += LRfrag2[s];
	    }
	  }
	  LR -= LR2;
	  if(DEBUG) assert(isfinite(LR));
	  if(LDEBUG) fragsum -= fragsum2;
	  int t = addcntIcum[S];
	  for(int T = 0; T < addcntI[S]; T++, t++){
	    if(DEBUG && !(isfinite(LRaddT[S][T]) && LRaddT[S][T] >= 0.0)){
	      #pragma omp critical
	      {
		printf("\n m=%d:S=%d(Amin=%d,Amax=%d,N=%d),T=%d,addcntIcum[S]=%d,addcnt=%d:LRaddT[S][T]=%0.10e\n",
		       m,S,Amin,Amax,N,T,addcntIcum[S],addcnt,(double)LRaddT[S][T]);
		fflush(stdout);
		assert(isfinite(LRaddT[S][T]) && LRaddT[S][T] >= 0.0);
	      }
	    }
	    FL128b LRa = (LR + LRaddT[S][T]) * LRratio;
	    if(DEBUG && !(isfinite(LRa) && (LRa <= 0.0 || LRa >= -LRcum*1e-8))){
	      printf("m=%d:S=%d,N=%d,T=%d,t=%d,addcnt=%d:LRa=%0.24Le:LRaddT[S][T]=%0.24Le,LR=%0.24Le,LRbias=%0.8e,LS=%0.8f,LRcum=%0.8e,Amin=%d,Amax=%d:addloc[t]=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f\n",
		     m,S,N,T,t,addcnt,(long double)LRa,(long double)LRaddT[S][T],(long double)LR, LRbias,LS,LRcum,Amin,Amax,addloc[t],Y[S],Y[S+1]);
	      fflush(stdout);
	      assert(isfinite(LRa) && (LRa <= 0.0 || LRa >= -LRcum*1e-8));
	    }
	    double LPa = (LRa > 0.0 || LRbiasLS > 0.0) ? LS + logl(max((FL128b)0.0,LRa) + LRbiasLS) : -99999.0;// HERE : vectorize call to logl()
	    if((DEBUG && !isfinite(LPa)) || (MDEBUG && MDEBUG_M >= 0 && MDEBUG_T >= 0 && m==MDEBUG_M && t==MDEBUG_T && N==MDEBUG_N)){
	      printf("m=%d:S=%d,N=%d,T=%d,t=%d,addcnt=%d:LPadd[t] += newLPa[m][t]=%13.10f (LR=%0.10e,LRbias=%0.8e,LS=%0.8f),Amin=%d,Amax=%d:addloc[t]=%0.3f,Y[S]=%0.3f,Y[S+1]=%0.3f\n",
		     m,S,N,T,t,addcnt,LPa,(double)LRa,LRbias,LS,Amin,Amax,addloc[t],Y[S],Y[S+1]);
	      printf("sum(LRadd[0..S]-LRadd2[0..S])=%0.10e,sum(LRadd2[0..S])=%0.10e,LRaddT[S][T]=%0.10e\n",(double)LR,(double)LR2,(double)LRaddT[S][T]);
	      fflush(stdout);
	      assert(isfinite(LPa));
	    }
	    if(DEBUG) assert(t < addcnt);

	    if(TDEBUG || MDEBUG || !newLPa){
              #pragma omp critical
	      {
		LPadd[t] += (LPa - LP) * TBwtm;
		if(DEBUG)assert(isfinite(LPadd[t]));
	      }
	    }
	    if(newLPa){
	      if(MDEBUG) newLPaPr[m][t] = (LRa < LRcum*1e-7) ? 1 : 0;/* loss of precision */
	      if(DEBUG>=2 && !VMEM_MINIMIZE) assert(newLPa[m][t] <= LARGE_NEGATIVE*0.99);
	      if(DEBUG) assert(TminM[m] <= t && t <= TmaxM[m]);
	      newLPa[m][t] = LPa;
	    }
	    LP = LPa;
	  } /* T, t loop */
	} /* S = Amin .. min(Amax,N) */
	if(Amax < N){
	  int t = addcntIcum[Amax+1];
	  if(t < addcnt){
	    FL128b LRd = LR * LRratio;/* to avoid log(0.0) when LR < 1.0e-306 */
	    double LPa = (LRd > 0.0 || LRbiasLS > 0.0) ? LS + logl(max((FL128b)0.0,LRd) + LRbiasLS) : -99999.0;// HERE : vectorize call to logl()
	    if(TDEBUG || MDEBUG || !newLPa){
              #pragma omp critical
	      {
		LPadd[t] += (LPa - LP) * TBwtm;
	      }
	    }
	    if(newLPa){
	      if(MDEBUG) newLPaPr[m][t] = (LR < LRcum*1e-7) ? 1 : 0; /* loss of precision */
	      if(DEBUG>=2 && !VMEM_MINIMIZE) assert(newLPa[m][t] <= LARGE_NEGATIVE*0.99);
	      if(DEBUG) assert(TminM[m] <= t && t <= TmaxM[m]);
	      newLPa[m][t] = LPa;
#ifdef _OPENMP
	      if(MDEBUG && MDEBUG_M >= 0 && MDEBUG_T >= 0 && t==MDEBUG_T && m==MDEBUG_M && N==MDEBUG_N){
		printf("m=%d:N=%d,t=%d,addcnt=%d:LPadd[t] += newLPa[m][t]=%13.10f (LR=%0.24Le,LRbias=%0.8e,LS=%0.8f), Amin=%d,Amax=%d:addloc[t]=%0.3f\n",
		       m, N, t, addcnt, LPa, (long double)LRd, LRbias, LS, Amin,Amax, addloc[t]);
		fflush(stdout);
	      }
#else
	      if(MDEBUG || DEBUG>=2){ /* verify that LPaddI+sum(LPadd[0..t]) = sum(newLPa[i=0..m][t]*TBmapWT[i]) */
		FL128 sum1 = LPaddI;
		for(int i = 0; i <= t; i++)
		  sum1 += LPadd[i];
		FL128 sum2 = 0.0;
		if(TBmapWT){
		  for(int i = 0; i <= m; i++){
		    int T = (t < TminM[i]) ? TminM[i] : (i < m && t > TmaxM[i]) ? TmaxM[i] : t;
		    if(DEBUG>=2 && !VMEM_MINIMIZE) assert(newLPa[i][T] > LARGE_NEGATIVE*0.99);
		    sum2 += newLPa[i][T] * TBmapWT[i];
		  }
		} else {
		  for(int i = 0; i <= m; i++){
		    int T = (t < TminM[i]) ? TminM[i] : (i < m && t > TmaxM[i]) ? TmaxM[i] : t;
		    if(DEBUG>=2 && !VMEM_MINIMIZE) assert(newLPa[i][T] > LARGE_NEGATIVE*0.99);
		    sum2 += newLPa[i][T];
		  }
		}
		if(DEBUG && (VERB>=2 || (MDEBUG_T >= 0 && t==MDEBUG_T && m==MDEBUG_M) || !(fabs(sum1-sum2) < max(1.0,fabs(sum1)) * CUMERR))){
		  int T = (t < TminM[m]) ? TminM[m] : t;
		  printf("m=%d:N=%d,t=%d,addcnt=%d:LPadd[0..t]=%0.10f,newLPa[0..m][t]=%0.10f (error = %0.14f),LPadd[t] += newLPa[m][t]=%13.10f (LR=%0.10e,LRbias=%0.8e,LS=%0.8f), addloc[t]=%0.3f\n",
			 m, N, t, addcnt, (double)sum1, (double)sum2, (double)(sum2-sum1), newLPa[m][T], (double)LR, LRbias,LS, addloc[t]);
		  fflush(stdout);
		  if(DEBUG && !(fabs(sum1-sum2) < max(1.0,fabs(sum1))*CUMERR)){
		    for(int i = 0; i <= t; i++)
		      printf("LPadd[i=%d]=%0.10f, TBmapWT[i]=%0.6f\n",i,(double)LPadd[i], (TBmapWT && i > 0) ? TBmapWT[i] : 1.0);
		    printf("  TBmapWT[m]=%0.10f,TBwtm=%0.10f,LPa=%0.10f,LP=%0.10f\n", TBmapWT[m],(double)TBwtm,LPa,LP);
		    fflush(stdout);
		    assert(fabs(sum1-sum2) < max(1.0,fabs(sum1))*CUMERR);
		  }
		}
	      }
#endif
	    }
	    LP = LPa;
	  } /* if(t < addcnt) */
	} /* if (Amax < N) */

	if(DEBUG>=2 && !VMEM_MINIMIZE){
	  for(int t = TminM[m]; t <= TmaxM[m]; t++){
	    if(!(newLPa[m][t] > LARGE_NEGATIVE*0.99)){
	      printf("mprobeval_cnt=%d:m=%d,addcnt=%d,Amin=%d,Amax=%d,N=%d:t=%d,newLPa[m][t]=%0.10e\n",mprobeval_cnt,m,addcnt,Amin,Amax,N,t,newLPa[m][t]);
	      fflush(stdout);
	      assert(newLPa[m][t] > LARGE_NEGATIVE*0.99);
	    }
	  }
	}

	if(VERB>=2 && tverb>=2 && (m==728 || m==9244) && LPdelta){
	  #pragma omp critical
	  {
	    printf("\t\t m=%d:tid=%d:After LPaddG[] updates:wall time= %0.6f\n", m, tid, wtime() - wt);
	    fflush(stdout);
	  }
	}
      } /* if(myaddcnt) */

      size_t LPdeltaSiz = 0;

      if(newLPdelta){

	/* update array LPdelta[S=0..N][d=0..D[S]-1] (difference from LPm) from LRdelta[S=0,Amin..Amax] and LRdeltaD[S=Amin..min(N,Amax)][d=0..D[S]-1] of current map m */
	/* LPdelta[S=Amin..min(N,Amax)][d=0..D[S]-1] += (log(LRdelta[0] + sum(LRdelta[max(1,Amin)..S]) + LRdeltaD[S][d] + LRbias) - LPm) * TBmapWT[m] */
	FL128b LR = LRdelta[0] * LRratio;
	if(DEBUG)assert(isfinite(LR));
	LP = (LR > 0.0 || LRbiasLS > 0.0) ? LS + logl(max((FL128b)0.0,(FL128b)LR) + LRbiasLS) : -99999.0;
	if(DEBUG && !isfinite(LP)){
	  printf("LR=%0.24Le,LRbiasLS=%0.8e,LS=%0.8f:LP=%0.8e\n",(long double)LR,LRbiasLS,LS,LP);
	  fflush(stdout);
	  assert(isfinite(LP));
	}
	if(1/*newLPdelta*/){
	  AminM[m] = Amin;
	  AmaxM[m] = Amax;
	  Doffset = Dcum[Amin];
	  /* allocation runs newLPdelta[m][(Dcum[Amin]..Dcum[min(N-1,Amax)+1]-1) - Doffset] */
	  if(VMEM_MINIMIZE){
	    LPdeltaSiz = Dcum[min(N-1,Amax)+1] - Dcum[Amin];
	    /*	    if(DEBUG && !(LPdeltaSiz < (size_t)(min(N-1,Amax)+1-Amin) * (2*localINITIAL_DELTA_RANGE + MAX_INDEL))){
	      printf("m=%d/%d:IMIN=%d,IMAX=%d,Amin=%d,Amax=%d,N=%d:localINITIAL_DELTA_RANGE=%d,MAX_INDEL=%d,LPdeltaSiz=%lu(Dcum[Amin]=%d,Dcum[min(N-1,Amax)+1]=%d)\n",
		     m,MD,IMIN,IMAX,Amin,Amax,N,localINITIAL_DELTA_RANGE,MAX_INDEL,LPdeltaSiz,Dcum[Amin],Dcum[min(N-1,Amax)+1]);
	      for(int A = Amin; A <= min(N-1,Amax); A++)
		printf("\t A=%d:Dcum[A]= %d, Dcum[A+1]-Dcum[A]= %d, D[A]= %d %c\n",A, Dcum[A], Dcum[A+1]-Dcum[A], D[A], Dcum[A+1]-Dcum[A] > (2*localINITIAL_DELTA_RANGE + MAX_INDEL) ? '!' : ' ');
	      fflush(stdout);
	      assert(LPdeltaSiz < (size_t)(min(N-1,Amax)+1-Amin) * (2*localINITIAL_DELTA_RANGE + MAX_INDEL));
	      }*/
	    if(DEBUG) assert(Dcum[min(N-1,Amax)+1] >= Dcum[Amin]);
	    newLPdelta[m] = thread_alloc(LPdeltaSiz, doubleAlloc, origY[IMIN], ppBlock,tid);// threadmem, threadmemsiz);
	  }
	} // newLPdelta

	for(int S = Amin; S <= min(N,Amax); S++){
	  if(S > 0)
	    LR += LRdelta[S];
	  if(DEBUG) assert(isfinite(LR));
	  for(int d = 0; d < D[S]; d++){
	    FL128b LRdel = (LR + LRdeltaD[S][d]) * LRratio;
	    double LPdel = (LRdel > 0.0 || LRbiasLS > 0.0) ? LS+logl(max((FL128b)0.0,LRdel) + LRbiasLS) : -99999.0;// HERE : vectorize call to logl() ???
	    if((DEBUG && !isfinite(LPdel)) || (MDEBUG && MDEBUG_M >= 0 && m==MDEBUG_M && S==MDEBUG_TS && d==MDEBUG_TT && N==MDEBUG_N)){
	      printf("\nm=%d:S=%d(N=%d),d=%d/%d:LPdelta[S][d] += LPdel = %0.10e (LRdel=%0.24Le,LRbiasLS=%0.8e,LS=%0.8f,LPm=%0.10f,LR=%0.24Le,LRdeltaD[S][d]=%0.24Le),Amin=%d,Amax=%d\n",
		     m,S,N,d,D[S],LPdel,(long double)LRdel,LRbiasLS,LS,LPm,(long double)LR,(long double)LRdeltaD[S][d],Amin,Amax);
	      fflush(stdout);
	      assert(isfinite(LPdel));
	    }
	    if((TDEBUG || MDEBUG /* || !newLPdelta */) && LPdelta){
              #pragma omp critical
	      {
		if(DEBUG>=2 && !(isfinite(LPdelta[S][d]))){
		  printf("\nm=%d:S=%d(N=%d),d=%d:LPdel = %13.10f (LR=%0.10e,LRbias=%0.8e,LS=%0.8f,LPm=%0.10f),Amin=%d,Amax=%d\n",
			 m,S,N,d,LPdel,(double)LRdel,LRbias,LS,LPm,Amin,Amax);
		  printf("LPdelta[S][d]=%0.10e,D[S]=%d\n",LPdelta[S][d],D[S]);
		  fflush(stdout);
		  assert(isfinite(LPdelta[S][d]));
		}

		LPdelta[S][d] += (LPdel - LPm) * TBwtm;/* LPm will be added back later */
		if(DEBUG && !(isfinite(LPdelta[S][d]))){
		  printf("m=%d:S=%d(N=%d),d=%d:LPdelta[S][d] += %13.10f - LPm (LR=%0.10e,LRbias=%0.8e,LS=%0.8f,LPm=%0.10f),Amin=%d,Amax=%d\n",
			 m,S,N,d,LPdel,(double)LRdel,LRbias,LS,LPm,Amin,Amax);
		  printf("LPdelta[S][d]=%0.10e,D[S]=%d\n",LPdelta[S][d],D[S]);
		  fflush(stdout);
		  assert(isfinite(LPdelta[S][d]));
		}
	      }
	    }// LPdelta
	    if(1/* newLPdelta */){
	      if(DEBUG) assert(Dcum[Amin] <= Dcum[S]+d && Dcum[S]+d < Dcum[min(N-1,Amax)+1]);
	      newLPdelta[m][Dcum[S]+d-Doffset] = LPdel;
	      if(VERB>=2 && S==1 && d==0){
		#pragma omp critical
		{
		  printf("m=%d,S=%d,d=%d:newLPdelta[m][%d]= %0.6f,Dcum[S]=%d,Doffset=%d\n", m,S,d,Dcum[S]+d-Doffset, LPdel, Dcum[S], Doffset);
		  fflush(stdout);
		}
	      }
	    }
	  }
	}/* S = Amin .. min(Amax,N) */

	if(VERB>=2 && tverb>=2 && (m==728 || m==9244)){
	  #pragma omp critical
	  {
	    printf("\t\t m=%d:tid=%d:After LPdelta[] updates:wall time= %0.6f\n", m, tid, wtime() - wt);
	    fflush(stdout);
	  }
	}
      }// if(newLPdelta)

      if(VERB>=3 && rverb){
	long long stride = (1LL + (CACHE ? kmax : A.kmax)) * A.Kstride;
	if(CACHE > 0 && ALIGN2 > 0)
	  stride = (stride + ALIGN2 - 1) & (~(ALIGN2 - 1));
        printf("\t m=%d/%d,tid=%d(cnt=%d),M=%d,len=%0.3f:A.Istride=%d,A.Kstride=%d,kmax=%d(stride=%lld),IMIN=%d,IMAX=%d(IMIN1=%d,IMAX1=%d),N=%d;LPaddSiz=%lu,LPdeltaSiz=%lu(Amin=%d,Amax=%d): Completed mprobeval: CPU=%0.6f,wall=%0.6f\n",
	       m,MD,tid,myLoopCnt,M,Xm[M+1],A.Istride,A.Kstride,(CACHE ? kmax : A.kmax), stride, IMIN,IMAX,IMIN1,IMAX1,N,LPaddSiz,LPdeltaSiz,Amin,Amax,mtime(),wtime());
        fflush(stdout);
      }
      if(VERB>=2 && tverb>=2 && LPdelta){
	#pragma omp critical
	{
	  double wt1 = wtime();
	  printf("\t m=%d:tid=%d,loopcnt=%d,wall=%0.6f secs(cum wall= %0.6f)\n",m,tid,myLoopCnt - myskipcnt,wt1-wt, wt1);
	  fflush(stdout);
	  wt = wt1;
	}
      }

      STOP(mprobeval_for);
    }// parallel for(int m = 0; m < MD; m++)
    STOP(mprobeval_for_wait);

    #pragma omp critical
    {
      if(OMP_DEBUG) LoopCnt += myLoopCnt;
      if(VERB_EPS2 && EPS2 > 0.0){
        origIJtot += origIJcnt;
	IJtot += IJcnt;
	LLRtotEPS2 += LLRsumEPS2;
      }
      setcnt += mysetcnt;
      qsetcnt += myqsetcnt;
      globalcnt += myglobalcnt;
      skipcnt += myskipcnt;
      if(VERB>=2 && (tverb>=2 || rverb)){
	double wt1 = wtime();
	printf("\t tid=%d(loopcnt=%d):wall=%0.6f sec (cum CPU= %0.6f, wall= %0.6f)\n",tid,myLoopCnt-myskipcnt,wt1-wt0,mtime(),wt1);
	fflush(stdout);
      }
    } 

    /* free up thread local memory, excluding lightweight_heap allocations */
    if(LDEBUG){
      delete [] LRfrag2;
      delete [] LRfrag;
      delete [] LRsite;
    }
    if(DIAGONAL && MIN_MEM && Fmem){

#if USE_MMAP
      free_huge_pages(Fmem);
#else
      //      free(Fmem);
      if(Fmem && munmap(Fmem, memsiz)){
	int eno = errno;
	char *err = strerror(eno);
	printf("tid=%d:munmap(%p,%lu) failed: errno=%d:%s\n",tid,Fmem,memsiz,eno,err);
	dumpmemmap();
	fflush(stdout);exit(1);
      }
      Fmem = 0;
#endif // USE_MMAP==0
    }

    STOP(mprobeval_parallel);
  } // pragma omp parallel
  STOP(mprobeval_parallel_wait);
  if(OMP_DEBUG) assert(LoopCnt == MD);// verify OMP loop worked
  
  if(newLPdelta){
    free(deltaFmem);deltaFmem = 0;
    delete [] deltaF; deltaF = 0;
  }

  if(VERB/* HERE >=2 */ && (rverb || tverb) /* && LPdelta*/){
    printf("\t parallel loop end : CPU=%0.3f,wall=%0.3f sec\n",mtime(),wtime());
    fflush(stdout);
  }

  if(TDEBUG || MDEBUG || !newLPd || !newLPa){
    if(LPdelG){  /* accumulate LPdel[S=0..N] */
      for(int S = 1; S <= N; S++)
	LPdel[S] += LPdel[S-1];
      for(int S = 0; S <= N; S++){
	LPdelG[S] = LPdel[S];
	if(DEBUG) assert(isfinite(LPdelG[S]));
      }
    }

    if(addcnt){    /* accumulate LPadd[t=0..addcnt-1] */
      LPadd[0] += LPaddI;
      for(int t = 1; t < addcnt; t++)
	LPadd[t] += LPadd[t-1];
      for(int t = 0; t < addcnt; t++){
	LPaddG[t] = LPadd[t];
	if(DEBUG) assert(isfinite(LPaddG[t]));
      }
    }
  }

  if(newLP){/* recompute LLRsum = sum(m=0..MD-1) newLP[m] * TBmapWT[m] so that is does not depend on thread ordering */
    double origLLRsum = LLRsum;
    LLRsum = 0.0;
    if(TBmapWT){
      for(int m = 0; m < MD; m++){
	LLRsum += newLP[m] * TBmapWT[m];
	if(DEBUG && !isfinite(LLRsum)){
	  printf("LLRsum= %0.8e (m=%d,newLP[m]=%0.8e,TBmapWT[m]=%0.8e)\n",LLRsum, m,newLP[m],TBmapWT[m]);
	  fflush(stdout);
	  assert(isfinite(LLRsum));
	}    
      }
    } else {
      for(int m = 0; m < MD; m++){
	LLRsum += newLP[m];
	if(DEBUG && !isfinite(LLRsum)){
	  printf("LLRsum= %0.8e (m=%d,newLP[m]=%0.8e)\n",LLRsum, m,newLP[m]);
	  fflush(stdout);
	  assert(isfinite(LLRsum));
	}    
      }
    }
    if(TDEBUG && !(fabs(origLLRsum - LLRsum) < 1e-7)){
      printf("origLLRsum=%0.10f,LLRsum=%0.10f(err=%0.10f)\n",origLLRsum,LLRsum,LLRsum - origLLRsum);
      fflush(stdout);
      assert(fabs(origLLRsum - LLRsum) < 1e-6);
    }

    if(VERB/* HERE >=2 */ && (rverb || tverb)){
      printf("\t After recomputing LLRsum= %0.10f -> %0.10f : CPU=%0.3f,wall=%0.3f sec\n",origLLRsum, LLRsum, mtime(),wtime());
      fflush(stdout);
    }
  }

  if(VERB_EPS2 && EPS2 > 0.0 && origIJtot > 0){
    double D = origIJtot;
    printf("\n\t EPS2 statistics: origIJtot = %lld, IJtot = %lld (%0.2f%% reduction), LLRsum = %0.6f -> %0.6f (delta= %0.6e, EPS2*MD=%0.6e)\n",
	   origIJtot, IJtot, (D - IJtot)*100.0/D, LLRtotEPS2, LLRsum, LLRtotEPS2 - LLRsum, EPS2 * MD);
    fflush(stdout);
    //    if(DEBUG>=2) assert(LLRtotEPS2 - LLRsum < EPS2*MD);
    //    exit(1);
  }

  if(LPdelG && newLPd){/* recompute LPdelG[S=1..N] = sum(newLPd[m=0..MD-1][S] * TBmapWT[m] (LPdelG[0] is not used)*/
    int nthreads = max(1,min(numthreads,N/8));

    #pragma omp parallel for schedule(static,8) num_threads(nthreads) if(nthreads>1)
    for(int S = 1; S <= N; S++){
      if(origY[S] < Ymin || Ymax < origY[S])
	continue;

      double sum = 0.0;
      if(TBmapWT){
	for(int m = 0; m < MD; m++){
          int s = (!VMEM_MINIMIZE && S < DminM[m]) ? 0 : min(S,DmaxM[m]);// NEW8
	  sum += ((VMEM_MINIMIZE && S < DminM[m]) ? newLPd0[m] : newLPd[m][s]) * TBmapWT[m];// NEW8
	}
      } else {
	for(int m = 0; m < MD; m++){
          int s = (!VMEM_MINIMIZE && S < DminM[m]) ? 0 : min(S,DmaxM[m]);// NEW8
	  sum += (VMEM_MINIMIZE && S < DminM[m]) ? newLPd0[m] : newLPd[m][s];// NEW8
	}
      }
      if(TDEBUG && !(fabs(sum - LPdelG[S]) < 1e-7)){
	printf("WARNING:LPdelG[%d]=%0.10f -> %0.10f (err=%0.10f)\n",S,LPdelG[S],sum,sum-LPdelG[S]);
	fflush(stdout);
	assert(fabs(sum - LPdelG[S]) < 1e-6);
      }
      LPdelG[S] = sum;
    }

    if(VERB/* HERE >=2 */ && (rverb || tverb)){
      printf("\t After recomputing LPdelG: CPU=%0.3f,wall=%0.3f sec\n",mtime(),wtime());
      fflush(stdout);
    }
  }

  if(addcnt && newLPa){/* recompute LPaddG[t = 0..addcnt-1] = sum(newLPa[m=0..MD-1][t]*TBmapWT[m]) */
    /* HERE : speed up by using differentials */

    int nthreads = max(1,min(numthreads,addcnt/64));

    #pragma omp parallel for schedule(static,64) num_threads(nthreads) if(nthreads>1)
    for(int t = 0; t < addcnt; t++){
      double sum = 0.0;
      if(TBmapWT){
	for(int m = 0; m < MD; m++){
	  int T = (t < TminM[m]) ? TminM[m] : (t > TmaxM[m]) ? TmaxM[m] : t;
	  if(DEBUG>=2 && !VMEM_MINIMIZE) assert(newLPa[m][T] > LARGE_NEGATIVE*0.99);
	  sum += newLPa[m][T] * TBmapWT[m];
	}
      } else {
	for(int m = 0; m < MD; m++){
	  int T = (t < TminM[m]) ? TminM[m] : (t > TmaxM[m]) ? TmaxM[m] : t;
	  if(DEBUG>=2 && !VMEM_MINIMIZE) assert(newLPa[m][T] > LARGE_NEGATIVE*0.99);
	  sum += newLPa[m][T];
	}
      }
      if(TDEBUG && !(fabs(sum - LPaddG[t]) < 1e-7)){
	printf("WARNING:mprobeval:t=%d(addloc[t]=%0.3f):LPaddG[t]=%0.10f,sum(newLPa[0..MD-1][t]*TBmapWT[m])=%0.10f (err=%0.14f)\n", t, addloc[t],LPaddG[t],sum, sum-LPaddG[t]);
	double sum = 0.0, sumLP= 0.0;
	for(int m = 0; m < MD; m++){
	  int T = (t < TminM[m]) ? TminM[m] : (t > TmaxM[m]) ? TmaxM[m] : t;
	  sum += newLPa[m][T];
	  sumLP += newLP[m];
	  printf("m=%d:TminM[m]=%d,TmaxM[m]=%d:newLP[m]=%0.10f(cum=%0.10f),newLPa[m][%d]=%0.10f(cum=%0.10f)\n",
		 m,TminM[m],TmaxM[m],newLP[m],sumLP,T,newLPa[m][T],sum);
	}
	fflush(stdout);
	assert(fabs(sum - LPaddG[t]) < 1e-6);
      }
      LPaddG[t] = sum;
    }
    if(VERB/* HERE >=2 */ && (rverb || tverb)){
      printf("\t After recomputing LPaddG: CPU=%0.3f,wall=%0.3f sec\n",mtime(),wtime());
      fflush(stdout);
    }
  }

  if(LPdelta){    /* NOTE: Can speed up by using differentials (relative to newLR[m]) */
    if(DEBUG) assert(newLPdelta);

    int nthreads = max(1,min(numthreads,N/16));

    double mt1 = 0.0,wt1 = 0.0;
    if(VERB/* HERE >=2 */ && (rverb || tverb)){
      mt1 = mtime();
      wt1 = wtime();
    }

    #pragma omp parallel for schedule(static,16) num_threads(nthreads) if(nthreads>1)
    for(int S = 0; S <= N; S++){
      for(int d = 0; d < D[S]; d++){/* recompute LPdelta[S][d] = sum(newLPdelta[m=0..MD-1][Dcum[S]-Dcum[Amin[m]]+d]*TBmapWT[i]) - LLRsum */
        int SDoffset = Dcum[S]+d;
	double sum = -LLRsum;
	if(TBmapWT)
	  for(int m = 0; m < MD; m++)
	    sum += ((AminM[m] <= S && S <= AmaxM[m]) ? newLPdelta[m][SDoffset - Dcum[AminM[m]]] : newLP[m]) * TBmapWT[m];
	else
	  for(int m = 0; m < MD; m++){
            double inc = (AminM[m] <= S && S <= AmaxM[m]) ? newLPdelta[m][SDoffset - Dcum[AminM[m]]] : newLP[m];
	    if(VERB>=2 && S==1 && d==0){
              printf("S=%d,d=%d,m=%d:Amin=%d,Amax=%d,Dcum[Amin]=%d,Dcum[S]=%d:inc=%0.6f,sum=%0.6f\n",
                 S,d,m,AminM[m],AmaxM[m],Dcum[AminM[m]],Dcum[S],inc,sum);
	      fflush(stdout);
            }
	    sum += inc;
          }
	if(TDEBUG && !(fabs(sum - LPdelta[S][d]) < 1e-7)){
	  printf("WARNING:mprobeval:S=%d/%d,d=%d/%d:LPdelta[S][d]=%0.10f -> %0.10f (err=%0.10f)\n",S,N,d,D[S],LPdelta[S][d],sum,sum-LPdelta[S][d]);
	  fflush(stdout);
	  assert(fabs(sum - LPdelta[S][d]) < 1e-6);
	}
	LPdelta[S][d] = sum;
      }
    } 

    if(VERB>=2){
      for(int S = 0; S <= N; S++)
	for(int d = 0; d < D[S]; d++)
	  printf("S=%d,d=%d:LPdelta[S][d]=%0.6f,%0.6f\n",S,d,LPdelta[S][d],LPdelta[S][d]+LLRsum);
      fflush(stdout);
    }
    if(VERB/* HERE >=2 */ && (rverb || tverb)){
      double mt2 = mtime(), wt2 = wtime();
      printf("\t After recomputing LPdelta: CPU=%0.6f,wall=%0.6f (cum= %0.6f, %0.6f) sec\n",mt2-mt1,wt2-wt1,mt2,wt2);
      fflush(stdout);
    }
  }// LPdelta

  if(SITE_PEN > 0){/* Bayesian prior penalty */
    LLRsum -= N*SITE_PEN;
    if(DEBUG && !isfinite(LLRsum)){
      printf("LLRsum= %0.8e (N=%d, SITE_PEN=%0.8e)\n",LLRsum, N, SITE_PEN);
      fflush(stdout);
      assert(isfinite(LLRsum));
    }    

    if(LPdelG){
      for(int S = 1; S <= N; S++){
        LPdelG[S] -= (N-1)*SITE_PEN;
        if(VERB>=2)
	  printf("I=%d:LPdel[I]=%0.6f, Y[I]=%0.3f (LP=%0.6f)\n",S, (double)LPdel[S],origY[S], LLRsum);
      }
    }

    for(int t = 0; t < addcnt; t++){
      LPaddG[t] -= (N+1)*SITE_PEN;
      if(VERB>=2){
	printf("t=%d,addcnt=%d:LPadd[t]=%0.10e, addloc[t]=%0.3f (LP=%0.10e, N=%d)\n",
	       t,addcnt,(double)LPadd[t], addloc[t], LLRsum, N);
	fflush(stdout);
      }
    }

    if(VERB/* HERE >=2 */ && (rverb || tverb)){
      printf("\t After adjusting LPdelG(N=%d),LPaddG(addcnt=%d) for SITE_PEN: CPU=%0.3f,wall=%0.3f sec\n",LPdelG ? N : 0, addcnt, mtime(),wtime());
      fflush(stdout);
    }
  }

  if(LPdelta) {/* add back LLRsum to LPdelta[S=0..N][d=0..D[S]] */
    int nthreads = max(1,min(numthreads,N/64));

    #pragma omp parallel for schedule(static,64) num_threads(nthreads) if(nthreads>1)
    for(int S = 0; S <= N; S++){
      for(int d = 0; d < D[S]; d++)
	LPdelta[S][d] += LLRsum;
    }
  }
  if(newLPdelta)
    delete [] mindelta;

  if(LPdelG && N <= 2)/* avoid deleting sites, since score_init() will fail if N <= 1 */
    for(int S = 1; S <= N; S++)
      LPdelG[S] = LLRsum - 1000.0;

  delete [] addcntI;
  delete [] addlocI;
  delete [] LPadd;
  delete [] LPaddTI;
  if(PADDING && addlocImem){
    free(addlocImem);
    addlocImem = NULL;
  }

  if(HEAP_MINIMIZE){/* free lightweight heap to minimize real memory usage for calls to qprobeval() */
    delete heap;
    heap = new lightweight_heap(0,0);      
  }

  if(MPROBEVAL_MAPORDER && SETLIMIT_CACHE)
    delete [] maporder;

  if(VERB/* HERE >=2 */){
    printf("skip=%d/%d/%d/%d,gl=%d,",MD-setcnt-qsetcnt,MD-setcnt,skipcnt,MD,globalcnt);
    fflush(stdout);
  }
  if(DEBUG && !isfinite(LLRsum)){
    printf("LLRsum= %0.8e\n",LLRsum);
    fflush(stdout);
    assert(isfinite(LLRsum));
  }
  STOP(mprobeval);
  return LLRsum;
}

#if USE_MIC==0
__attribute__ ((noinline))
double mprobeval(int n, double *Hcuts,/**< Hcuts[i=0..n+1] : complete consensus map : only used to set alignment bounds */
		 int N, double *origY,/**< Y[I=0..N+1] : current consensus map hypothesis */
		 int MD, /**< number of maps */
		 int *MX, double **X, /**< X[m=0..MD-1][j=0..MX[m]+1] : Map sites */
		 int lc,/**< If left end of Y is a linear chromosome end */
		 int rc,/**< If right end of Y is a linear chromosome end */
		 int **map, /**< map[m=0..MD-1][j=0..MX[m]+1] is index in Hcuts[0..n+1] of X[m][j], 
					map[MD][j=0..n+1] is index in Y[0..N+1] of Hcuts[j] */
		 int **mapK, /**< mapK[m=0..MD-1][j=0..MX[m]+1] is K index offset in Hcuts[0..n+1] of X[m][j] */ 
		 Csetlimit *limit,
		 int *nmapMD, /**< nmapMD[j=0..N+1] is index in Hcuts[0..n+1] of Y[j] and always defined */
		 double *TBmapWT, /**< If != 0 : TBmapWT[m] is a the -TB based weight (multiplied by pcontig->mapWT[m] due to -Erefine or -BestRefWT) 
				     for map m to be applied when summing up the total log(LR+LRbias) value(s) */
		 double *newLP,/**< return value of Log Likelihood for each Map */
		 double *oldLP,/**< If != 0 : previous estimate of Log Likelihood for each map */
		 double Ymin, double Ymax, /**< location range in Y[] that changed since oldLP was computed (or WILL change due to addloc[] or delta[][])  : only valid if oldLP != 0 && Ymin >= 0.0 */
		 int addcnt, double *addloc,/**< addloc[t=0..addcnt-1] are the locations at which a site will be added to Y[] (one at a time) : Y[0] < addloc[t] < Y[N+1] */
		 int *D, double **delta,/**< If LPdelta != 0 : delta[I=0..N][d=0..D[I]-1] is the amount by which the interval Y[I+1]-Y[I] is incremented */
		 double *LPdelG,/**< If != 0 : LPdelG[I=1..N] is the LP value when the site Y[I] is deleted from Y[0..N+1] (NOTE: LPdelG[I] for Y[I] outside of Ymin..Ymax are NOT computed) */
		 double *LPaddG,/**< If addcnt != 0 : LPadd[t=0..addcnt-1] is the LP value when a site is added at addloc[t] to Y[0..N+1] */
		 double *LPgrad[2],/**< If != 0 : LPgrad[d=0,1][I=0..N] is (d+1)'th gradient of LP vs Y[I+1]-Y[I] (Not yet implemented) */
		 double **LPdelta, /**< If != 0 : LPdelta[I=0..N][d=0..D[I]-1] is the LP value when Y[I+1]-Y[I] is incremented by delta[I][d] */
		 double **newLPd, /**< If != 0 : newLPd[m=0..MD-1][I=DminM[m]..DmaxM[m]] is the LP value of map m when the site Y[I] is deleted from Y[0..N+1]
				     Values newLPd[m][I] for 1 <= I < DminM[m] are not defined but should be assumed equal to newLPd[m][0] 
				     Values newLPd[m][I] for DmaxM[m] < I <= N are not defined but should be assumed equal to newLPd[m][DmaxM[m]] (NOTE: newLPd[m][I] for Y[I] outside of Ymin..Ymax are NOT computed
				     Also if SKIPMAP_ADD && addcnt > 0 && Ymin < 0.0 AND no addloc[t] value falls within Y[I-1]..Y[I+1], then newLPd[m][I] is forced to newLPd[m][0]).
				     NOTE: If VMEM_MINIMIZE, the 2nd index of newLPd[][] has not been allocated and newLPd[m][0] is replaced by newLPd0[m] */
                 double *newLPd0,/* IFF VMEM_MINIMIZE then newLPd0[m] replaces newLPd[m][0] (which is no longer defined, unless DminM[m] <= 0) */
		 int **newLPdPr,  /**< If MDEBUG : newLPdPr[m=0..MD-1][I=DminM[m]..DmaxM[m]] is 1 IFF newLPd[I][m] is low precision (due to drop in LP) */
		 int *DminM,  /**< DminM[m=0..MD-1] : see newLPd[] */ 
		 int *DmaxM,  /**< DmaxM[m=0..MD-1] : see newLPd[] */
		 double **newLPa, /**< If != 0 : newLPa[m=0..MD-1][t=TminM[m]..TmaxM[m]] is the LP value of map m when a site is added at addloc[t] to Y[0..N+1]
				     Values newLPa[m][t] for 0 <= t < TminM[m] are not defined and should be assumed equal to newLPa[m][TminM[m]]
				     Values newLPa[m][t] for TmaxM[m] < t < addcnt are not defined and should be assumed equal to newLPa[m][TmaxM[m]]
				     NOTE: If VMEM_MINIMIZE, the 2nd index of newLPa[][] has not been allocated. The allocation list is stored in DoubleAllocation[] */
		 int **newLPaPr,  /**< If MDEBUG : newLPaPr[m=0..MD-1][t=TminM[m]..TmaxM[m]] is 1 IFF newLPa[t][m] is low precision (due to drop in LP)   */
		 int *TminM,  /**< TminM[m=0..MD-1] : see newLPa[] */ 
		 int *TmaxM,  /**< TmaxM[m=0..MD-1] : see newLPa[] */
		 double **newLPdelta, /**< If != 0 : newLPdelta[m=0..MD-1][Dcum[I=AminM[m]..min(N-1,AmaxM[m])] + (d=0..D[I]-1) - Doffset] is the LP value of map m when Y[I+1]-Y[I] is incremented by delta[I][d] 
					 To save real memory, the offset Doffset == Dcum[AminM[m]] is used to shift the used memory to the left end of each m's memory block.
					 NOTE: If VMEM_MINIMIZE, only the first index has been allocated and the 2nd index is allocated for each m by this function. The allocation list is stored in DoubleAllocation[] */
		 int *Dcum,   /* Dcum[I=0..N] = sum(t = 0..I-1) max(0,D[I]) */
		 /* save real memory by using flattened newLPdelta[m=0..MD-1][d = 0 .. (sum(I=0..N)D[I])-1], which saves most of the 2nd index memory that is never de-referenced */
		 int *AminM, /**< AminM[m=0..MD-1] : see newLPdelta[] */
		 int *AmaxM, /**< AmaxM[m=0..MD-1] : see newLPdelta[] */
		 Ccontig *pcontig, /**< point to complete contig information for debugging */
		 char *skipmap, /** If != 0 : skipmap[m=0..MD-1] == 1 IFF map m is already known not to overlap intervals Y[I..I+1] with D[I] > 0. skipmap[] is updated with each call */
		 int Allele, /* If != 0, identifies the Allele that Hcuts[] and origY[] correspond to (1 or 2 if called from HaploTypeR) */
		 CdoubleBlock **pBlock /* If != 0, preallocated memory pBlock[0..MD-1], that is used to save a pointer to the CdoubleBlock used by newLPdelta[m][] or newLPa[m][], so it can be swapped out and restored on demand
					   by the caller, instead of all at once by calling SwappedAllocations_restore() */
		 )
{
  int OUT_TYPE = OUTLIER_TYPE;
  int OUT_OLD = OUTLIER_OLD(outlierMaxBC);

  if(OUT_TYPE){
    if(OUT_OLD)
      return mprobeval<1,1>(n,Hcuts,N,origY,MD,MX,X,lc,rc,map,mapK,limit,nmapMD,TBmapWT,newLP,oldLP,Ymin,Ymax,addcnt,addloc,D,delta,LPdelG,LPaddG,LPgrad,LPdelta,newLPd,newLPd0,newLPdPr,DminM,DmaxM,newLPa,newLPaPr,TminM,TmaxM,newLPdelta,Dcum,AminM,AmaxM,pcontig,skipmap,Allele,pBlock);
    else
      return mprobeval<1,0>(n,Hcuts,N,origY,MD,MX,X,lc,rc,map,mapK,limit,nmapMD,TBmapWT,newLP,oldLP,Ymin,Ymax,addcnt,addloc,D,delta,LPdelG,LPaddG,LPgrad,LPdelta,newLPd,newLPd0,newLPdPr,DminM,DmaxM,newLPa,newLPaPr,TminM,TmaxM,newLPdelta,Dcum,AminM,AmaxM,pcontig,skipmap,Allele,pBlock);
  } else {
    if(OUT_OLD)
      return mprobeval<0,1>(n,Hcuts,N,origY,MD,MX,X,lc,rc,map,mapK,limit,nmapMD,TBmapWT,newLP,oldLP,Ymin,Ymax,addcnt,addloc,D,delta,LPdelG,LPaddG,LPgrad,LPdelta,newLPd,newLPd0,newLPdPr,DminM,DmaxM,newLPa,newLPaPr,TminM,TmaxM,newLPdelta,Dcum,AminM,AmaxM,pcontig,skipmap,Allele,pBlock);
    else
      return mprobeval<0,0>(n,Hcuts,N,origY,MD,MX,X,lc,rc,map,mapK,limit,nmapMD,TBmapWT,newLP,oldLP,Ymin,Ymax,addcnt,addloc,D,delta,LPdelG,LPaddG,LPgrad,LPdelta,newLPd,newLPd0,newLPdPr,DminM,DmaxM,newLPa,newLPaPr,TminM,TmaxM,newLPdelta,Dcum,AminM,AmaxM,pcontig,skipmap,Allele,pBlock);
  }
}
#endif // USE_MIC==0

} // namespace probeval
