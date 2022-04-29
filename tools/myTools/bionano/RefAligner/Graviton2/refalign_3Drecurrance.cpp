#ifndef REFALIGN_3DRECURRANCE
#define REFALIGN_3DRECURRANCE

static Ident refalign_3Drecurrance_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/refalign_3Drecurrance.cpp 11381 2020-07-28 23:26:39Z tanantharaman $");

#define PREFETCH 1 /* 0 : no prefetch instructions
		      1 : minimal prefetch instructions (1.5% speedup) */

/* refalign_3Drecurrance<maptype> handles new AVX256 vectorized code and general case */

template<int maptype>
static void refalign_3Drecurrance(FLOAT *Y, int N, FLOAT *X, int M,
				  int IMIN, int IMAX, int *JMIN, int *JMAX, int *Kmax, 
				  AL *A, int refid, int mapid, int orientation, int scaleID, int rev,
				  Cmap *rmap, Cmap *nanomap, CYPen ****YPenR, CXPen *XPen, Cprtab **PRtabY,
				  int tid, int numthreads)
{
  if(1){
    if(VERB && (TSAN || !vector_warning || PVERB>=2)){
      #pragma omp critical(vector_warning)
      {
	if(!vector_warning){
	  if(USE_SSE && (USE_AVX || USE_AVX512) && USE_RFLOAT && !EXTENDONLY)
	    printf("Using vectorized code: MultiMatches=%d,deltaX=%d,deltaY=%d,outlierBC=%d,outlierExtend=%d,%d\n",
		   MultiMatches,DELTA_X,DELTA_Y,outlierBC,outlierExtend,outlierExtendLimit);
	  else
	    printf("WARNING: Not using vectorized code: MultiMatches=%d,deltaX=%d,deltaY=%d,outlierBC=%d,outlierExtend=%d,%d\n",
		   MultiMatches,DELTA_X,DELTA_Y,outlierBC,outlierExtend,outlierExtendLimit);
	  fflush(stdout);
	  vector_warning = 1;
	}
      }
    }

#if !(USE_SSE && (USE_AVX || USE_AVX512) && USE_RFLOAT && EXTENDONLY==0) 
    RFLOAT ExtLen = (extendonly && Refine==2) ? (2.0*Ylambda + EndLen) : 0.0;
#endif

    if(PVERB){
      printf("refid=%d(%lld),mapid=%d(%lld),or=%d,rev=%d,sc=%d:IMIN=%d,IMAX=%d:general case\n",
	refid,rmap->id, mapid, nanomap->id, orientation,rev,scaleID,IMIN,IMAX);
      fflush(stdout);
    }

    for(int I=IMIN; I <= IMAX; I++){
      int jmin = JMIN[I];
      int jmax = JMAX[I];
      if((PVERB && I == I_TRACE)){
	#pragma omp critical
	{
	  printf("refid=%d(%lld),mapid=%d(%lld),or=%d:I=%d:JMIN[I]=%d,JMAX[I]=%d,Kmax[I]=%d\n",
		 refid,rmap->id, mapid, nanomap->id, orientation,I,JMIN[I],JMAX[I],Kmax[I]);
	  fflush(stdout);
	}
      }
      for(int K= Kmax[I]; K >= 0; K--){
        if(DEFER_BP && !outlierExtend){/* faster case without updating BackPointer A->G,H,T */
          if(DEBUG>=2) assert(!MultiMatches);

#if USE_SSE && (USE_AVX || USE_AVX512) && USE_RFLOAT && EXTENDONLY==0 // vectorized AVX256 code
	  int Gmin = max(IMIN,I - K - DELTA_Y);
	  RFLOAT *AscoreIK = &A->score(I,K,0);

	  for(int G= I - K; --G >= Gmin;){
	    int HminG = JMIN[G];
	    int HmaxG = JMAX[G];
	    int jminG = max(jmin,HminG+1);
	    int jmaxG = min(jmax,HmaxG+DELTA_X);

	    for(int T = Kmax[G]; T >= 0; T--){
	      if(DEBUG>=2) assert(G-T >= 1);
	      RFLOAT *AscoreGT = &A->score(G,T,0);
	      if(PREFETCH){
		int mmax = min(DELTA_X, min(jminG+7,jmaxG) - HminG);
		_mm_prefetch((const char *)&AscoreGT[jminG-mmax], _MM_HINT_T0);
		_mm_prefetch((const char *)&AscoreGT[jminG-mmax+8], _MM_HINT_T0);
	      }
	      RFLOAT deltaY,Ivar,GaussY,Sm;
	      RFLOAT penY = SintY(I-K-G,I,K,T,deltaY,Ivar,GaussY,Sm,YPenR);
	      RFLOAT OutPenBC = OutlierPenaltyBC[I-K-G];

	      __m256 v_deltaY = _mm256_set1_ps(deltaY);
	      __m256 v_penY = _mm256_set1_ps(penY);
	      __m256 v_GaussY = _mm256_set1_ps(GaussY);
	      __m256 v_Sm = _mm256_set1_ps(Sm);
	      __m256 v_Ivar = _mm256_set1_ps(Ivar);

	      for(int J = jminG; J <= jmaxG; J += 8){
		int mmax = min(DELTA_X, min(J+7, jmaxG) - HminG);
		if(PREFETCH) _mm_prefetch((const char *)&AscoreGT[J+16-mmax], _MM_HINT_T0);
		int maskJ = (1u << min(8,jmaxG-J+1))-1;// MASK(min(8,jmaxG-J+1))
		if(DEBUG>=2) assert(0 < maskJ && maskJ <= 255);
		__m256 v_score = _mm256_loadu_ps(&AscoreIK[J]);

		for(int m = 1; m <= mmax; m++){// NOTE : vectorizes iterations j = J..J+7 (m == j-H). Only save results from iterations with j <= jmaxG && (j-HmaxG <= m && m <= j-HminG)
		  // Encode condition (j <= jmaxG && j-HmaxG <= m && m <= j-HminG) as 8 bit mask : bit 0 for j=J, bit 1 for j=J+1 ... bit 7 for j=J+7
		  int maskJm = maskJ & ((1u << max(0,min(8,m+HmaxG-J+1)))-1) & ~((1u << max(0,m+HminG-J))-1);// maskJ & MASK(max(0,min(8,m+HmaxG-J+1))) & ~MASK(max(0,m+HminG-J))
		  // if(!maskJm) continue; // faster to run this loop even if none of the 8 sub-iterations will be used!
		  if(DEBUG>=2) assert(maskJm <= 255);
		  __m256 maskJm256 = mask8to256[maskJm];/* convert 8 bit mask to 8x32 bit mask by table lookup */
		  if(DEBUG>=2){/* verify accuracy of maskJm & maskJm256 */
		    unsigned int mask[8];
		    _mm256_storeu_si256((__m256i*)mask,_mm256_castps_si256(maskJm256));
		    for(int t = 0; t < 8; t++){
		      assert(((maskJm >> t) & 1) == (J+t <= jmaxG && J+t-HmaxG <= m && m <= J+t-HminG));
		      assert(mask[t] == (J+t <= jmaxG && J+t-HmaxG <= m && m <= J+t-HminG ? 0xffffffff : 0));
		    }
		  }

		  __m256 v_pscore = _mm256_loadu_ps(&AscoreGT[J-m]);
		  __m256 v_deltaX = _mm256_loadu_ps(&XPen->deltaX[m][J]);
		  __m256 XPenJM_Bias = _mm256_loadu_ps(&XPen->Bias[m][J]);
		  __m256 XPenJM_Pen = _mm256_loadu_ps(&XPen->Pen[m][J]);
		  __m256 XPenJM_PenBias = _mm256_loadu_ps(&XPen->PenBias[m][J]);
		  __m256 v_newscore = v_pscore + SintX_mm256<maptype>(v_deltaX,v_deltaY,v_penY,v_GaussY,v_Sm,v_Ivar,XPenJM_Bias,XPenJM_Pen,XPenJM_PenBias, OutPenBC + OutlierPenaltyBC[m]);
		  if(DEBUG>=2){/* verify correctness of v_newscore */
		    float newscore[8];
		    _mm256_storeu_ps(newscore,v_newscore);
		    for(int t = 0; t < 8; t++){
		      if(!(J+t <= jmaxG && J+t-HmaxG <= m && m <= J+t-HminG))
			continue;/* no need to check since newscore[t] will not be used */
		      RFLOAT newscoret = AscoreGT[J+t-m] + SintX(XPen->deltaX[m][J+t],deltaY,penY,GaussY,Sm,Ivar,XPen->Bias[m][J+t],XPen->Pen[m][J+t],XPen->PenBias[m][J+t], OutPenBC + OutlierPenaltyBC[m]);
		      if(fabs(newscoret - newscore[t]) >= fabs(newscoret)*1e-5 + 1e-5){
			printf("refid=%d(%lld),mapid=%d(%lld),or=%d:I=%d,K=%d,G=%d,T=%d,J=%d,m=%d:t=%d:newscoret= %0.6f, newscore[t]= %0.6f\n",
			       refid,rmap->id,mapid,nanomap->id,orientation,I,K,G,T,J,m,t,newscoret,newscore[t]);
			fflush(stdout);
			assert(fabs(newscoret - newscore[t]) < fabs(newscoret)*1e-5 + 1e-5);
		      }
		    }
		  }

		  __m256 scoreGT = _mm256_and_ps(maskJm256, _mm256_cmp_ps(v_newscore, v_score, _CMP_GT_OS));// if (newscore > score  && J <= jmaxG && J-HmaxG <= m && m <= J-HminG) 
		  v_score = _mm256_blendv_ps(v_score, v_newscore, scoreGT);
		}/* m = 1 .. DELTA_X */

		_mm256_storeu_ps(&AscoreIK[J], v_score);
	      } // for(int J = jminG; J <= jmaxG; J += 8)
	    } // T = Kmax[G] .. 0
	  } // G = I-K .. Gmin
#else // original non-vectorized code
	  FLOAT Yik = Yc(Y,I,K);

	  for(int J= jmin; J <= jmax; J++){
            if(EXTENDONLY && ExtLen && Yik - X[J] > ExtLen && Y[N+1]-Yik - (X[M+1]-X[J]) > ExtLen)
	      continue; /* skip if map is not within ExtLen = (2*Ylambda + EndLen) kb of extending beyond the end of the reference */

	    RFLOAT score = A->score(I,K,J);
	    if((DEBUG>=2 && !isfinite(score)) || (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE)){
	      #pragma omp critical
	      {
	        if(DEFER_BP)
		  printf("refid=%d,mapid=%d:I=%d,K=%d,J=%d:A->score(I,K,J)=%0.8f\n",refid,mapid,I,K,J,score);
		else
		  printf("refid=%d,mapid=%d:I=%d,K=%d,J=%d:A->score(I,K,J)=%0.8f, G=%d,T=%d,H=%d\n",refid,mapid,I,K,J,score, A->G(I,K,J),A->T(I,K,J),A->H(I,K,J));
		fflush(stdout);
		assert(isfinite(score));
	      }
	    }
	    int Gmin = max(IMIN, I - K - DELTA_Y);
	    for(int G= I-K; --G >= Gmin;){
	      int Hmin = max(JMIN[G],J-DELTA_X);
	      int Hmax = min(JMAX[G],J-1);
	      for(int T = Kmax[G]; T >= 0; T--){
		if(DEBUG>=2) assert(G-T >= 1);
		RFLOAT *AscoreGT = &A->score(G,T,0);

		RFLOAT deltaY,Ivar,GaussY,Sm;
		RFLOAT penY = SintY(I-K-G,I,K,T,deltaY,Ivar,GaussY,Sm,YPenR);
		RFLOAT OutPenBC = OutlierPenaltyBC[I-K-G];
		for(int H = Hmax; H >= Hmin;H--){
		  if(DEBUG>=2) assert(H >= (max(1,JMIN[G])));
		  if(DEBUG>=2) assert(H <= (min(J,JMAX[G])));
		  int m = J-H;
		  RFLOAT deltaX = XPen->deltaX[m][J];
		  RFLOAT newscore = AscoreGT[H] + SintX(deltaX,deltaY,penY,GaussY,Sm,Ivar,XPen->Pen[m][J],XPen->Bias[m][J],XPen->PenBias[m][J], OutPenBC + OutlierPenaltyBC[m]);
		  if((DEBUG>=2 && !isfinite(newscore)) || (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE) /* || (PVERB && I==129 && K==0 && J==56) || (PVERB && I==130 && K==0 && J==57) || (PVERB && I==131 && K==0 && J==59)*/){
                    #pragma omp critical 
		    {
		      printf("E:refid=%d,mapid=%d:I=%d,K=%d/%d,J=%d,G=%d,T=%d/%d,H=%d,IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d:AscoreGT[H]=%0.8f,newscore=%0.8f,score=%0.8f\n",
			     refid,mapid,I,K,Kmax[I],J,G,T,Kmax[G],H,IMIN,IMAX,JMIN[I],JMAX[I],JMIN[G],JMAX[G],A->score(G,T,H),newscore,max(score,newscore));
		      fflush(stdout);
		      assert(isfinite(newscore));
		    }
		  }
		  score = max(score,newscore);
		}/* H = Hmax .. Hmin */
	      } /* T = Kmax[G] .. 0 */
	    } /* G = I-K-1 ... Gmin */
	    if((DEBUG>=2 && !isfinite(score)) || (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE) /* || (PVERB && I==129 && K==0 && J==56) || (PVERB && I==130 && K==0 && J==57) || (PVERB && I==131 && K==0 && J==59)*/){
	      #pragma omp critical
	      {
		printf("F:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d:score-> %0.8f\n",refid,mapid,orientation,I,K,J,score);
		fflush(stdout);
		assert(isfinite(score));
	      }
	    }
	    A->score(I,K,J) = score;
	  } // J = jmin .. jmax
#endif // non-vectorized code
        } else {/* !DEFER_BP || outlierExtend */
#if USE_SSE && (USE_AVX || USE_AVX512) && USE_RFLOAT && EXTENDONLY==0 // vectorized AVX256 code
	  int Gmin = max(IMIN, I - K - DELTA_Y);
	  RFLOAT *AscoreIK = &A->score(I,K,0);
	  int *AG_IK = &A->G(I,K,0);
	  int *AH_IK = &A->H(I,K,0);
	  int *AT_IK = &A->T(I,K,0);

	  for(int G= I - K; --G >= Gmin;){
	    int HminG = JMIN[G];
	    int HmaxG = JMAX[G];
	    int jminG = max(jmin,HminG+1);
	    int jmaxG = min(jmax,HmaxG+DELTA_X);

	    __m256i v_G = _mm256_set1_epi32(G);

	    for(int T = Kmax[G]; T >= 0; T--){
	      RFLOAT *AscoreGT = &A->score(G,T,0);
	      if(PREFETCH){
		int mmax = min(DELTA_X, min(jminG+7,jmaxG) - HminG);
		_mm_prefetch((const char *)&AscoreGT[jminG-mmax], _MM_HINT_T0);
		_mm_prefetch((const char *)&AscoreGT[jminG-mmax+8], _MM_HINT_T0);
	      }

	      RFLOAT deltaY,Ivar,GaussY,Sm;
	      RFLOAT penY = SintY(I-K-G,I,K,T,deltaY,Ivar,GaussY,Sm,YPenR);
	      RFLOAT OutPenBC = OutlierPenaltyBC[I-K-G];

	      __m256i v_T = _mm256_set1_epi32(T);
	      __m256 v_deltaY = _mm256_set1_ps(deltaY);
	      __m256 v_penY = _mm256_set1_ps(penY);
	      __m256 v_GaussY = _mm256_set1_ps(GaussY);
	      __m256 v_Sm = _mm256_set1_ps(Sm);
	      __m256 v_Ivar = _mm256_set1_ps(Ivar);

	      for(int J = jminG; J <= jmaxG; J += 8){
		int maskJ = (1u << min(8,jmaxG-J+1))-1;// MASK(min(8,jmaxG-J+1))
		if(DEBUG>=2) assert(0 < maskJ && maskJ <= 255);
		int mmax = min(DELTA_X, min(J+7, jmaxG) - HminG);// avoid m values with maskJm==0 (or that might access AscoreGT[J-m] before AscoreGT[HminG-8])
		__m256 v_score = _mm256_loadu_ps(&AscoreIK[J]);
		__m256i v_bestG = _mm256_loadu_si256((const __m256i*)&AG_IK[J]);
		__m256i v_bestH = _mm256_loadu_si256((const __m256i*)&AH_IK[J]);
		__m256i v_bestT = _mm256_loadu_si256((const __m256i*)&AT_IK[J]);
		if(PREFETCH) _mm_prefetch((const char *)&AscoreGT[J+16-mmax], _MM_HINT_T0);

		for(int m = 1; m <= mmax; m++){// NOTE : Vectorizes iterations j = J..J+7 with m==j-H. Only saves results from iterations with j <= jmaxG && (j-HmaxG <= m && m <= j-HminG)
		  // Encode condition (j <= jmaxG && j-HmaxG <= m && m <= j-HminG) as 8 bit mask : bit 0 for j=J, bit 1 for j=J+1 ... bit 7 for j=J+7
		  int maskJm = maskJ & ((1u << max(0,min(8,m+HmaxG-J+1)))-1) & ~((1u << max(0,m+HminG-J))-1);// maskJ & MASK(max(0,min(8,m+HmaxG-J+1))) & ~MASK(max(0,m+HminG-J))
		  // if(!maskJm) continue; // faster to run this loop even if none of the 8 iterations will be used!
		  if(DEBUG>=2) assert(maskJm <= 255);
		  __m256 maskJm256 = mask8to256[maskJm];/* convert 8 bit mask to 8x32 bit mask by table lookup */
		  if(DEBUG>=2){/* verify accuracy of maskJm & maskJm256 */
		    unsigned int mask[8];
		    _mm256_storeu_si256((__m256i*)mask,_mm256_castps_si256(maskJm256));
		    for(int t = 0; t < 8; t++){
		      if(!(((maskJm >> t) & 0x1) == ((J+t <= jmaxG && J+t-HmaxG <= m && m <= J+t-HminG) ? 1 : 0))){
                        #pragma omp critical
			{
			  printf("refid=%d(%lld),mapid=%d(%lld),or=%d:I=%d,K=%d,G=%d,T=%d,J=%d,m=%d:maskJm= 0x%x,t=%d:maskJm[t]=%d,jmaxG=%d,HminG=%d,HmaxG=%d,H=%d\n",
				 refid,rmap->id,mapid,nanomap->id,orientation,I,K,G,T,J,m,(unsigned int)maskJm,t,((maskJm >> t) & 0x1),jmaxG,HminG,HmaxG,J+t-m);
			  fflush(stdout);
			  assert(((maskJm >> t) & 0x1) == ((J+t <= jmaxG && J+t-HmaxG <= m && m <= J+t-HminG) ? 1 : 0));
			}
		      }
		      assert(mask[t] == ((J+t <= jmaxG && J+t-HmaxG <= m && m <= J+t-HminG) ? 0xffffffff : 0x0));
		    }
		  }

		  __m256 v_pscore = _mm256_loadu_ps(&AscoreGT[J-m]);
#if LOAD_ALIGNED
		  __m256 v_deltaX = _mm256_load_ps(&XPen->deltaX[m][J]);
		  __m256 XPenJM_Bias = _mm256_load_ps(&XPen->Bias[m][J]);
		  __m256 XPenJM_Pen = _mm256_load_ps(&XPen->Pen[m][J]);
		  __m256 XPenJM_PenBias = _mm256_load_ps(&XPen->PenBias[m][J]);
#else
		  __m256 v_deltaX = _mm256_loadu_ps(&XPen->deltaX[m][J]);
		  __m256 XPenJM_Bias = _mm256_loadu_ps(&XPen->Bias[m][J]);
		  __m256 XPenJM_Pen = _mm256_loadu_ps(&XPen->Pen[m][J]);
		  __m256 XPenJM_PenBias = _mm256_loadu_ps(&XPen->PenBias[m][J]);
#endif

		  __m256 v_newscore = v_pscore + SintX_mm256<maptype>(v_deltaX,v_deltaY,v_penY,v_GaussY,v_Sm,v_Ivar,XPenJM_Bias,XPenJM_Pen,XPenJM_PenBias, OutPenBC + OutlierPenaltyBC[m]);

		  if(DEBUG>=2 || PVERB){/* verify correctness of v_newscore (or display values) */
		    float newscore[8],score[8];
		    _mm256_storeu_ps(score,v_score);
		    _mm256_storeu_ps(newscore,v_newscore);
		    for(int t = 0; t < 8; t++){
		      if(!(J+t <= jmaxG && J+t-HmaxG <= m && m <= J+t-HminG))
			continue;/* no need to check since newscore[t] will not be used */
		      int H = J+t - m;
		      RFLOAT deltaX = XPen->deltaX[m][J+t];
		      RFLOAT newscoret = AscoreGT[H] + SintX(deltaX,deltaY,penY,GaussY,Sm,Ivar,XPen->Bias[m][J+t],XPen->Pen[m][J+t],XPen->PenBias[m][J+t], OutPenBC + OutlierPenaltyBC[m]);
		      if(DEBUG && fabs(newscoret - newscore[t]) >= fabs(newscoret)*1e-5 + 1e-5){
                        #pragma omp critical
			{
                          printf("refid=%d(%lld),mapid=%d(%lld),or=%d:I=%d,K=%d,G=%d,T=%d,J=%d,m=%d:t=%d:newscoret= %0.6f, newscore[t]= %0.6f\n",
                     	    refid,rmap->id,mapid,nanomap->id,orientation,I,K,G,T,J,m,t,newscoret,newscore[t]);
			  fflush(stdout);
			  assert(fabs(newscoret - newscore[t]) < fabs(newscoret)*1e-5 + 1e-5);
                        }
                      }

		      if((DEBUG && !isfinite(newscore[t])) || (PVERB && ((I== I_TRACE && K== K_TRACE && J+t== J_TRACE)) /*  && newscore[t] > score[t])*/)){
                        #pragma omp critical 
			{
			  int bestG[8],bestH[8],bestT[8];
			  _mm256_storeu_si256((__m256i *)bestG,v_bestG);		      
			  _mm256_storeu_si256((__m256i *)bestT,v_bestT);		      
			  _mm256_storeu_si256((__m256i *)bestH,v_bestH);		      

			  RFLOAT Bias,Pen,Gauss,PenSm,iscore;
			  SintDetail(deltaX,deltaY, J+t-H, I-K-G, J+t,I,K,T, PRtabY, Y, Bias,Pen,Gauss,PenSm, 0);
			  RFLOAT OutPen = OutlierPenalty;
			  if(outlierBC)
			    OutPen += OutlierPenaltyBC[J+t-H] + OutlierPenaltyBC[I-K-G];
			  if(maptype){
			    if(OUTLIER_LTYPE==0)
			      OutPen -= (deltaX+deltaY) * OutlierLambdaInv;
			    else
			      OutPen -= fabs(deltaX-deltaY) * OutlierLambdaInv;
			    iscore = OutlierBias + max(Bias + Pen + Gauss, Bias * biasWToutlierF + Pen * OUTLIER_TYPE1 + OutPen) + PenSm;
		          } else if(OUTLIER_DELTA(deltaX-deltaY)){
			    OutPen -= fabs(deltaX-deltaY) * OutlierLambdaInv;
			    iscore = OutlierBias + max(Bias + Pen + Gauss, Bias * biasWToutlierF + Pen*OUTLIER_TYPE + OutPen) + PenSm;
		          } else
			    iscore = Bias + Pen + Gauss + PenSm;

			  RFLOAT sint = SintX(deltaX,deltaY,penY,GaussY,Sm,Ivar,XPen->Pen[m][J+t],XPen->Bias[m][J+t],XPen->PenBias[m][J+t], OutPenBC + OutlierPenaltyBC[m]);

			  double var = SF[0]*SF[0] + fabs(SD[0])*SD[0]*deltaY;
			  if(QUADRATIC_VARIANCE)
			    var += SR[0]*SR[0]*deltaY*deltaY;
			  if(RES_VARIANCE){
			    FLOAT resR = Y[I] - Y[I-K];
			    FLOAT resL = Y[G] - Y[G-T];
			    FLOAT resvar = resL*resL + resR*resR;
			    var += SE[0]*SE[0]*resvar;
		          }
			  if(MultiMatches)
			    printf("  I=%d,K=%d,J=%d:G=%d,T=%d,H=%d(%d..%d)(IL=%d,KL=%d,JL=%d):deltaX=%0.3f,deltaY=%0.3f,norm=%0.4f,AscoreGT[H]=%0.6f,newscore=%0.6f,score=%0.6f (SintX=%0.6f,iscore=%0.6f,Bias=%0.6f,Pen=%0.6f,Gauss=%0.6f,PenSm=%0.6f,OutBias=%0.6f,OutPen=%0.6f),best G=%d,T=%d,H=%d\n",
			       I,K,J+t,G,T,H,HminG,HmaxG,A->IL(G,T,H),A->KL(G,T,H),A->JL(G,T,H),deltaX,deltaY,(deltaX-deltaY)*(deltaX-deltaY)/var,AscoreGT[H],newscore[t],score[t],sint,iscore,Bias,Pen,Gauss,PenSm,Bias*biasWToutlierF,OutPen,bestG[t],bestT[t],bestH[t]);
			  else
			    printf("  I=%d,K=%d,J=%d:G=%d,T=%d,H=%d(%d..%d):deltaX=%0.3f(SNR=%0.1f),deltaY=%0.3f,norm=%0.4f,AscoreGT[H]=%0.6f,newscore=%0.6f,score=%0.6f (SintX=%0.6f,iscore=%0.6f,Bias=%0.6f,Pen=%0.6f,Gauss=%0.6f,PenSm=%0.6f),best G=%d,T=%d,H=%d\n",
			       I,K,J+t,G,T,H,HminG,HmaxG,deltaX, MapSNR ? Gmap[mapid]->SNR[0][orientation ? M+1-J-t : M] : -1.0, deltaY,(deltaX-deltaY)*(deltaX-deltaY)/var,AscoreGT[H],newscore[t],score[t],sint,iscore,Bias,Pen,Gauss,PenSm,bestG[t],bestT[t],bestH[t]);
			  fflush(stdout);
			  assert(isfinite(newscore[t]));
		        }// pragma omp critical
		      }//if((DEBUG..) || (PVERB ...))
	            }// t = 0 .. 7
		  }// if(DEBUG>=2 || PVERB)

		  __m256 scoreGT = _mm256_and_ps(maskJm256, _mm256_cmp_ps(v_newscore, v_score, _CMP_GT_OS));// scoreGT[j=0..7] == (newscore > score && J+j <= jmaxG && J+j-HmaxG <= m && m <= J+j-HminG) ? ~0:0
		  v_score = _mm256_blendv_ps(v_score, v_newscore, scoreGT);
#ifdef __AVX2__
		  __m256i v_H = _mm256_add_epi32(_mm256_set1_epi32(J-m), v256i_76543210);// convert J-m into a vector J-m,J-m+1,J-m+2,J-m+3..,J-m+7 (low to high order 32 bit ints)
#else // AVX only
		  __m128i v_H1 = _mm_add_epi32(_mm_set1_epi32(J-m), v128i_3210);
		  __m128i v_H2 = _mm_add_epi32(_mm_set1_epi32(J-m), v128i_7654);
		  __m256i v_H = _mm256_set_m128i(v_H2, v_H1);
#endif

		  v_bestG = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(v_bestG), _mm256_castsi256_ps(v_G), scoreGT));
		  v_bestH = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(v_bestH), _mm256_castsi256_ps(v_H), scoreGT));
		  v_bestT = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(v_bestT), _mm256_castsi256_ps(v_T), scoreGT));		  
	        }/* m = 1 .. min(DELTA_X, jmaxG - HminG) */

	        // use unmasked store : should be faster and will just write back the same values for some extra iterations of J
		_mm256_storeu_ps(&AscoreIK[J], v_score);
		_mm256_storeu_si256((__m256i *)&AG_IK[J], v_bestG);
		_mm256_storeu_si256((__m256i *)&AT_IK[J], v_bestT);
		_mm256_storeu_si256((__m256i *)&AH_IK[J], v_bestH);
	      } // for(int J = jminG; J <= jmaxG; J += 8)
	    } // T = Kmax[G] .. 0
	  } // G = I-K .. Gmin

#else // original non-vectorized code

	  for(int J= jmin; J <= jmax; J++){
	    RFLOAT score = A->score(I,K,J);
	    int bestG = A->G(I,K,J);
	    int bestH = A->H(I,K,J);
	    int bestT = A->T(I,K,J);

	    int Gmin = max(IMIN,I - K - DELTA_Y);

	    if(PVERB>=2 && I_TRACE >= 0 && I== I_TRACE && K== K_TRACE && J== J_TRACE){
	      if(MultiMatches)
		printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:I=%d,K=%d,J=%d(IL=%d,KL=%d,JL=%d):IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,G=%d..%d\n",
		       refid,rmap->id,mapid,nanomap->id,orientation,I,K,J,A->IL(I,K,J),A->KL(I,K,J),A->JL(I,K,J),IMIN,IMAX,JMIN[I],JMAX[I],I-K-1,Gmin);
	      else
		printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:I=%d,K=%d,J=%d:IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,G=%d..%d\n",
		       refid,rmap->id,mapid,nanomap->id,orientation,I,K,J,IMIN,IMAX,JMIN[I],JMAX[I],I-K-1,Gmin);
	      fflush(stdout);
	    }

	    for(int G= I - K; --G >= Gmin;){
	      int Hmin = max(JMIN[G], J - DELTA_X);
	      int Hmax = min(JMAX[G], J - 1);
	      if(PVERB>=2 && I_TRACE >= 0 && I== I_TRACE && K== K_TRACE && J== J_TRACE /* && G==385 && T==0 && H==11031 */){
		if(MultiMatches)
		  printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:I=%d,K=%d,J=%d(IL=%d,KL=%d,JL=%d):IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,G=%d,T=0..%d,H=%d..%d\n",
			 refid,rmap->id,mapid,nanomap->id,orientation,I,K,J,A->IL(I,K,J),A->KL(I,K,J),A->JL(I,K,J),IMIN,IMAX,JMIN[I],JMAX[I],G,Kmax[G],Hmin,Hmax);
		else
		  printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:I=%d,K=%d,J=%d:IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,G=%d,T=0..%d,H=%d..%d\n",
			 refid,rmap->id,mapid,nanomap->id,orientation,I,K,J,IMIN,IMAX,JMIN[I],JMAX[I],G,Kmax[G],Hmin,Hmax);
		fflush(stdout);
	      }
	      for(int T = Kmax[G]; T >= 0; T--){
		if(DEBUG>=2) assert(G-T >= 1);
		RFLOAT deltaY,Ivar,GaussY,Sm;
		RFLOAT penY = SintY(I-K-G,I,K,T,deltaY,Ivar,GaussY,Sm,YPenR);
		RFLOAT OutPenBC = OutlierPenaltyBC[I-K-G];
		RFLOAT *AscoreGT = &A->score(G,T,0);
		for(int H = Hmax; H >= Hmin;H--){
		  if(DEBUG>=2) assert(H >= (max(1,JMIN[G])));
		  if(DEBUG>=2) assert(H <= (min(J,JMAX[G])));
		  int m = J-H;
		  RFLOAT deltaX = XPen->deltaX[m][J];
		  RFLOAT newscore = AscoreGT[H] + SintX(deltaX,deltaY,penY,GaussY,Sm,Ivar,XPen->Pen[m][J],XPen->Bias[m][J],XPen->PenBias[m][J], OutPenBC + OutlierPenaltyBC[m]);
		  if((DEBUG>=2 && !isfinite(newscore)) || (PVERB /* && outlierExtend */ && ((I== I_TRACE && K== K_TRACE && J== J_TRACE /* && ((G >= I - 2 || H >= J - 2) || newscore > score)*/)))){
                    #pragma omp critical 
		    {
		      RFLOAT Bias,Pen,Gauss,PenSm,iscore;
		      SintDetail(deltaX,deltaY, J-H, I-K-G, J,I,K,T, PRtabY, Y, Bias,Pen,Gauss,PenSm, 0);
		      RFLOAT OutPen = OutlierPenalty;
		      if(outlierBC)
			OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-K-G];
		      if(maptype){
			if(OUTLIER_LTYPE==0)
			  OutPen -= (deltaX+deltaY) * OutlierLambdaInv;
			else
			  OutPen -= fabs(deltaX-deltaY) * OutlierLambdaInv;
			iscore = OutlierBias + max(Bias + Pen + Gauss, Bias * biasWToutlierF + Pen * OUTLIER_TYPE1 + OutPen) + PenSm;
		      } else if(OUTLIER_DELTA(deltaX-deltaY)){
			OutPen -= fabs(deltaX-deltaY) * OutlierLambdaInv;
			iscore = OutlierBias + max(Bias + Pen + Gauss, Bias * biasWToutlierF + Pen*OUTLIER_TYPE + OutPen) + PenSm;
		      } else
			iscore = Bias + Pen + Gauss + PenSm;

		      RFLOAT sint = SintX(deltaX,deltaY,penY,GaussY,Sm,Ivar,XPen->Pen[m][J],XPen->Bias[m][J],XPen->PenBias[m][J], OutPenBC + OutlierPenaltyBC[m]);

		      double var = SF[0]*SF[0] + fabs(SD[0])*SD[0]*deltaY;
		      if(QUADRATIC_VARIANCE)
			var += SR[0]*SR[0]*deltaY*deltaY;
		      if(RES_VARIANCE){
			FLOAT resR = Y[I] - Y[I-K];
			FLOAT resL = Y[G] - Y[G-T];
			FLOAT resvar = resL*resL + resR*resR;
			var += SE[0]*SE[0]*resvar;
		      }
		      if(MultiMatches)
			printf("  I=%d,K=%d,J=%d:G=%d,T=%d,H=%d(%d..%d)(IL=%d,KL=%d,JL=%d):deltaX=%0.3f,deltaY=%0.3f,norm=%0.4f,AscoreGT[H]=%0.6f,newscore=%0.6f,score=%0.6f (SintX=%0.6f,iscore=%0.6f,Bias=%0.6f,Pen=%0.6f,Gauss=%0.6f,PenSm=%0.6f,OutBias=%0.6f,OutPen=%0.6f),best G=%d,T=%d,H=%d\n",
			       I,K,J,G,T,H,Hmin,Hmax,A->IL(G,T,H),A->KL(G,T,H),A->JL(G,T,H),deltaX,deltaY,(deltaX-deltaY)*(deltaX-deltaY)/var,AscoreGT[H],newscore,score,sint,iscore,Bias,Pen,Gauss,PenSm,Bias*biasWToutlierF,OutPen,bestG,bestT,bestH);
		      else
			printf("  I=%d,K=%d,J=%d:G=%d,T=%d,H=%d(%d..%d):deltaX=%0.3f(SNR=%0.1f),deltaY=%0.3f,norm=%0.4f,AscoreGT[H]=%0.6f,newscore=%0.6f,score=%0.6f (SintX=%0.6f,iscore=%0.6f,Bias=%0.6f,Pen=%0.6f,Gauss=%0.6f,PenSm=%0.6f),best G=%d,T=%d,H=%d\n",
			       I,K,J,G,T,H,Hmin,Hmax,deltaX, MapSNR ? Gmap[mapid]->SNR[0][orientation ? M+1-J : M] : -1.0, deltaY,(deltaX-deltaY)*(deltaX-deltaY)/var,AscoreGT[H],newscore,score,sint,iscore,Bias,Pen,Gauss,PenSm,bestG,bestT,bestH);
		      fflush(stdout);
		      assert(isfinite(newscore));
		    }
		  }
		  if(newscore > score){
		    score = newscore;
		    bestG = G;
		    bestH = H;
		    bestT = T;
		  }
		}/* H = Hmax .. Hmin */
	      } /* T = Kmax[G] .. 0 */
	    } /* G = I-K-1 ... Gmin */

	    if(PVERB && I== I_TRACE && K== K_TRACE && J== J_TRACE){
	      #pragma omp critical
	      {
		if(bestG <= 0){
		  printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:I=%d,K=%d,J=%d:best G=%d,T=%d,H=%d,score=%0.6f(Sm=%0.6f) (DELTA_X=%d,DELTA_Y=%d,outlierExtend=%d,%d)\n",
			 refid,rmap->id, mapid, nanomap->id, orientation, I,K,J,bestG,bestT,bestH,score, Sm(0,I,K,Y), DELTA_X,DELTA_Y,outlierExtend,outlierExtendLimit);
		} else {
		  printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:I=%d,K=%d,J=%d:best G=%d,T=%d,H=%d,score=%0.6f(Sm=%0.6f),y=%0.3f,x=%0.3f (DELTA_X=%d,DELTA_Y=%d,outlierExtend=%d,%d)\n",
			 refid,rmap->id, mapid, nanomap->id, orientation, I,K,J,bestG,bestT,bestH,score, Sm(0,I,K,Y), 
			 Yc(Y,I,K) - Yc(Y,bestG,bestT), X[J]-X[bestH],DELTA_X,DELTA_Y,outlierExtend,outlierExtendLimit);
		  for(int k = K; k > 0; k--)
		    printf("   k=%d:Y[I-k]=%0.3f,Y[I-k+1]=%0.3f:y=%0.3f,Pn(y)=%0.6f,log(p*Pn(y))=%0.6f\n",k,Y[I-k],Y[I-k+1],Y[I-k+1]-Y[I-k],Pn(Y[I-k+1]-Y[I-k]),log(pTP*Pn(Y[I-k+1]-Y[I-k])));
		  for(int pI = bestG, pK = bestT, pJ = bestH; pI > 0;){
		    int pG = A->G(pI,pK,pJ);
		    int pT = A->T(pI,pK,pJ);
		    int pH = A->H(pI,pK,pJ);
		    RFLOAT pscore = A->score(pI,pK,pJ);
		    
		    if(pG > 0)
		      printf("      I=%d,T=%d,J=%d (G=%d,T=%d,H=%d) score=%0.6f(Sm=%0.6f),y=%0.3f,x=%0.3f\n",pI,pK,pJ,pG,pT,pH,pscore, Sm(0,pI,pK,Y), Yc(Y,pI,pK)-Yc(Y,pG,pT), X[pJ]-X[pH]);
		    else
		      printf("      I=%d,T=%d,J=%d (G=%d,T=%d,H=%d) score=%0.6f(Sm=%0.6f)\n",pI,pK,pJ,pG,pT,pH,pscore, Sm(0,pI,pK,Y));
		    pI = pG;
		    pK = pT;
		    pJ = pH;
		  }
		}
		fflush(stdout);
	      }
	    }

	    A->score(I,K,J) = score;
	    A->G(I,K,J) = bestG;
	    A->T(I,K,J) = bestT;
	    A->H(I,K,J) = bestH;
          } // J = jmin .. jmax
#endif // non-vectorized code

	  if(outlierExtend || MultiMatches){ // this code is NOT vectorized
	    for(int J= jmin; J <= jmax; J++){
	      RFLOAT score = A->score(I,K,J);
	      int bestG = A->G(I,K,J);
	      int bestT = A->T(I,K,J);
	      int bestH = A->H(I,K,J);

	      if(outlierExtend){
		/* try to merge with any previous outlier interval that ends at I-K or J */
		// NOTE1 : actual score (as computed by Sint() with larger deltaX) may be better than computed here (with outlierExtend) if merger of 
		//         two outlier intervals is not an outlier : can only happen with (maptype ? OUTLIER_TYPE1 : OUTLIER_TYPE0) with back to back "insertion" + "deletion" of similar size

		RFLOAT OutlierPenalty3 = OutlierPenalty2;
		if(DEBUG>=2) assert(K <= KMAX && I-K >= 1);
		OutlierPenalty3 += PRtabY[K][I].Sm;// Sm(0,I,K,Y);

		/* first merge any previous interval that ends at G==I-K && H < J, score the combined interval as an outlier */
		RFLOAT Yik = Yc(Y,I,K);
		int Gmin2 = max(IMIN+1, I - K - DELTA_Y/* HERE outlierExtend */);
		int G = I-K;
		if(G >= Gmin2){
		  int Hmin = max(JMIN[G], J - DELTA_X/* HERE outlierExtend */);
		  int Hmax = min(JMAX[G], J - 1);
		  for(int T = Kmax[G]; T >= 0; T--){
		    if(DEBUG>=2) assert(G-T >= 1);
		    for(int H = Hmax; H >= Hmin; H--){
		      if(DEBUG>=2) assert(H >= (max(1,JMIN[G])));
		      if(DEBUG>=2) assert(H <= (min(J,JMAX[G])));
		      int G2 = A->G(G,T,H);
		      if(G2 > 0){
			int T2 = A->T(G,T,H);
			int H2 = A->H(G,T,H);
			if(DEBUG>=2) assert(G2-T2 >= 1 );		
			if(DEBUG>=2) assert(H2 >= (max(1,JMIN[G2])));
			if(DEBUG>=2) assert(H2 <= (min(H,JMAX[G2])));
			if(!outlierExtendLimit || max(I-K-G2,J-H2) <= outlierExtendLimit){
			  RFLOAT x = X[J] - X[H2], y = Yik - Yc(Y,G2,T2);
			  if(OUTLIER_DELTA(x-y)){ // OutlierPenalty3 == OutlierPenalty + OutlierBias + PRtabY[K][I].Sm
			    RFLOAT deltascore = OutlierPenalty3 + OutPenExtend(x,y, J-H2, I-K-G2, I, K, PRtabY, refid);
			    if(DEBUG>=2){
			      RFLOAT Bias,Pen,Gauss,PenSm;
			      SintDetail(x,y, J-H2,I-K-G2,J,I,K,T2,PRtabY,Y,Bias,Pen,Gauss,PenSm,0);

			      RFLOAT OutPen = OutlierPenalty;
			      if(outlierBC)
				OutPen += OutlierPenaltyBC[J-H2] + OutlierPenaltyBC[I-K-G2];
			      if(OUTLIER_LTYPE==0)
				OutPen -= (x+y) * OutlierLambdaInv;
			      else
				OutPen -= fabs(x-y) * OutlierLambdaInv;
			      RFLOAT iscore2 = OutlierBias + Pen*OUTLIER_TYPE1 + Bias * biasWToutlierF + OutPen + PenSm;
			      
			      if(!(fabs(iscore2 - deltascore) < SCORE_MARGIN * fabs(iscore2) + ASCORE_MARGIN)){
				printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d,rev=%d:I=%d,K=%d,J=%d, G2=%d,T2=%d,H2=%d:x= %0.4f, y= %0.4f (DELTA_X=%d,DELTA_Y=%d,outlierExtend=%d,%d)\n",
				       refid,rmap->id, mapid, nanomap->id, orientation, rev, I,K,J,G2,T2,H2,x,y, DELTA_X,DELTA_Y,outlierExtend,outlierExtendLimit);
				printf("\t deltascore= %0.8f: OutlierPenalty3= %0.8f, OutPenExtend(x,y,J-H2,I-K-G2,I,K)= %0.8f\n",
				       deltascore,OutlierPenalty3, OutPenExtend(x,y,J-H2,I-K-G2,I,K,PRtabY,refid));
				printf("\t iscore2= %0.8f: Bias=%0.8f,Pen=%0.8f,Gauss=%0.8f,PenSm=%0.8f; OutPen=%0.8f,OutlierBias=%0.8f,biasWToutlierF=%0.8f\n",
				       iscore2,Bias,Pen,Gauss,PenSm,OutPen,OutlierBias,biasWToutlierF);
				fflush(stdout);
				assert(fabs(iscore2 - deltascore) < SCORE_MARGIN * fabs(iscore2) + ASCORE_MARGIN);
			      }
			    }
			    RFLOAT newscore = A->score(G2,T2,H2) + deltascore;
			    if(PVERB && I== I_TRACE && K== K_TRACE && J== J_TRACE && newscore > score){
			      #pragma omp critical
			      {
				printf("   A: I=%d,K=%d,J=%d: G=%d,T=%d,H=%d: G2=%d,T2=%d,H2=%d,score=%0.6f(%0.6f + %0.6f)%c best=%0.6f\n",I,K,J,G,T,H,G2,T2,H2,newscore,A->score(G2,T2,H2),
				       OutlierPenalty3 + OutPenExtend(X[J]-X[H2],Yc(Y,I,K)-Yc(Y,G2,T2), J-H2, I-K-G2, I, K, PRtabY, refid), (newscore > score) ? '!' : ' ',score);
				fflush(stdout);
			      }
			    }
			    if(newscore > score){
			      score = newscore;
			      bestG = G2;
			      bestH = H2;
			      bestT = T2;
			    }
			  }
			}
		      }
		    }
		  }
		}

		/* try merging previous intervals that end at G < I-K && H == J and score the combined interval as an outlier */
		int H = J;
		if(H >= 2){
		  for(int G= I - K; --G >= Gmin2;){
		    if(JMIN[G] <= H && H <= JMAX[G]){
		      for(int T = Kmax[G]; T >= 0; T--){
			if(DEBUG>=2) assert(G-T >= 1);
			int G2 = A->G(G,T,H);
			if(G2 > 0){
			  if(DEBUG>=2) assert(IMIN <= G2 && G2 < G);
			  int T2 = A->T(G,T,H);
			  int H2 = A->H(G,T,H);
			  if(DEBUG>=2) assert(G2-T2 >= 1 );		
			  if(DEBUG>=2) assert(H2 >= (max(1,JMIN[G2])));
			  if(DEBUG>=2) assert(H2 <= (min(H,JMAX[G2])));
			  if(!outlierExtendLimit || max(I-K-G2, J-H2) <= outlierExtendLimit){
			    RFLOAT x = X[J] - X[H2], y = Yik - Yc(Y,G2,T2);
			    if(OUTLIER_DELTA(x-y)){ // OutlierPenalty3 == OutlierPenalty + OutlierBias + PRtabY[K][I].Sm
			      RFLOAT deltascore = OutlierPenalty3 + OutPenExtend(x, y, J-H2, I-K-G2, I, K, PRtabY, refid);// OutlierPenaltyBC[I-K-G2] + OutlierPenaltyBC[J-H2];
			      if(DEBUG>=2){
				RFLOAT Bias,Pen,Gauss,PenSm;
				SintDetail(x,y, J-H2,I-K-G2,J,I,K,T2,PRtabY,Y,Bias,Pen,Gauss,PenSm,0);
				
				RFLOAT OutPen = OutlierPenalty;
				if(outlierBC)
				  OutPen += OutlierPenaltyBC[J-H2] + OutlierPenaltyBC[I-K-G2];
				RFLOAT iscore2;
				if(maptype){
				  if(OUTLIER_LTYPE==0)
				    OutPen -= (x+y) * OutlierLambdaInv;
				  else
				    OutPen -= fabs(x-y) * OutlierLambdaInv;
				  iscore2 = OutlierBias + Pen*OUTLIER_TYPE1 + Bias * biasWToutlierF + OutPen + PenSm;
				} else {
				  OutPen -= fabs(x-y) * OutlierLambdaInv;				
				  iscore2 = OutlierBias + Pen*OUTLIER_TYPE + Bias * biasWToutlierF + OutPen + PenSm;
				}
				
				if(!(fabs(iscore2 - deltascore) < SCORE_MARGIN * fabs(iscore2) + ASCORE_MARGIN)){
				  printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d,rev=%d:I=%d,K=%d,J=%d, G2=%d,T2=%d,H2=%d:x= %0.4f, y= %0.4f (DELTA_X=%d,DELTA_Y=%d,outlierExtend=%d,%d)\n",
					 refid,rmap->id, mapid, nanomap->id, orientation, rev, I,K,J,G2,T2,H2,x,y, DELTA_X,DELTA_Y,outlierExtend,outlierExtendLimit);
				  printf("\t deltascore= %0.8f: OutlierPenalty3= %0.8f, OutPenExtend(x,y,J-H2,I-K-G2,I,K)= %0.8f\n",
					 deltascore,OutlierPenalty3, OutPenExtend(x,y,J-H2,I-K-G2,I,K,PRtabY,refid));
				  printf("\t iscore2= %0.8f: Bias=%0.8f,Pen=%0.8f,Gauss=%0.8f,PenSm=%0.8f; OutPen=%0.8f,OutlierBias=%0.8f,biasWToutlierF=%0.8f\n",
					 iscore2,Bias,Pen,Gauss,PenSm,OutPen,OutlierBias,biasWToutlierF);
				  fflush(stdout);
				  assert(fabs(iscore2 - deltascore) < SCORE_MARGIN * fabs(iscore2) + ASCORE_MARGIN);
				}
			      }
			      RFLOAT newscore = A->score(G2,T2,H2) + deltascore;
			      if(PVERB && orientation == OR_TRACE && I== I_TRACE && K== K_TRACE && J== J_TRACE && newscore > score){
#pragma omp critical
				{
				  printf("   B: I=%d,K=%d,J=%d: G=%d,T=%d,H=%d: G2=%d,T2=%d,H2=%d,score=%0.6f(%0.6f + %0.6f)%c best=%0.6f\n",I,K,J,G,T,H,G2,T2,H2,newscore,A->score(G2,T2,H2),
				    OutlierPenalty3 + OutPenExtend(X[J]-X[H2], Yc(Y,I,K) - Yc(Y,G2,T2), J-H2, I-K-G2, I, K, PRtabY, refid), (newscore > score) ? '!' : ' ',score);
				  fflush(stdout);
			        }
			      }
			      if(newscore > score){
				score = newscore;
				bestG = G2;
				bestH = H2;
				bestT = T2;
			      }
			    }
			  }
		        }
		      }
		    }
		  }
	        }

		if(bestG >= 1){/* check if last interval is an outlier */
		  if(score <= A->score(bestG,bestT,bestH) + OutlierPenalty3 + OutlierPenaltyBC[J-bestH]+OutlierPenaltyBC[I-K-bestG] + 1e-9){
		    /* try to extend outlier interval to the left by reducing bestG,bestH by up to outlierExtend sites */
		    if(PVERB && orientation == OR_TRACE && I== I_TRACE && K== K_TRACE && J== J_TRACE){
		      #pragma omp critical
		      {
		        printf("    I=%d,K=%d,J=%d:score=%0.6f : A->score(%d,%d,%d)=%0.6f,OutlierPenalty3=%0.6f: trying to extend outlier\n",
			     I,K,J,score, bestG,bestT,bestH,A->score(bestG,bestT,bestH),OutlierPenalty3 + OutlierPenaltyBC[J-bestH]+OutlierPenaltyBC[I-K-bestG]);
			fflush(stdout);
		      }
		    }
		    for(int cnt = 0; cnt < 100; cnt++){/* try multiple rounds of extension, each by up to outlierExtend labels, but not more than outlierExtendLimit total labels */
		      int origbestG = bestG, origbestH = bestH;
		      RFLOAT origscore = score;

		      int Gmin2 = max(IMIN, origbestG - outlierExtend);
		      if(outlierExtendLimit)
			Gmin2 = max(I-K - outlierExtendLimit, Gmin2);

		      /* NOTE : not clear if we need to consider cases G > origbestG or H > origbestH (but not both at the same time) */
		      for(int G = /* WAS origbestG+1 */ I-K; --G >= Gmin2; ){
			int Hmin2 = max(JMIN[G], origbestH - outlierExtend);
			int Hmax2 = /* WAS min(JMAX[G], origbestH) */ min(JMAX[G],(G <= origbestG) ? J-1 : origbestH);
			if(outlierExtendLimit)
			  Hmin2 = max(J - outlierExtendLimit, Hmin2);
			for(int T = Kmax[G]; T >= 0; T--){
			  if(DEBUG>=2) assert(G-T >= 1 );
			  RFLOAT *AscoreGT = &A->score(G,T,0);
			  for(int H = Hmax2; H >= Hmin2; H--){
			    if(DEBUG>=2) assert(H >= (max(1,JMIN[G])));
			    if(DEBUG>=2) assert(H <= (min(J,JMAX[G])));
			    if(DEBUG>=2 && outlierExtendLimit) assert(max(I-K-G, J-H) <= outlierExtendLimit);

			    if(G == origbestG && H == origbestH) 
			      continue;// NEW110

			    RFLOAT x = X[J]-X[H], y = Yik - Yc(Y,G,T);
			    if(OUTLIER_DELTA(x-y)){ // OutlierPenalty3 == OutlierPenalty + OutlierBias + PRtabY[K][I].Sm
			      RFLOAT deltascore = OutlierPenalty3 + OutPenExtend(x, y, J-H, I-K-G, I , K, PRtabY, refid); // OutlierPenaltyBC[I-K-G] + OutlierPenaltyBC[J-H];
			      if(DEBUG>=2){
				RFLOAT Bias,Pen,Gauss,PenSm;
				SintDetail(x,y, J-H,I-K-G,J,I,K,T,PRtabY,Y,Bias,Pen,Gauss,PenSm,0);

				RFLOAT OutPen = OutlierPenalty;
				if(outlierBC)
				  OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-K-G];
				if(OUTLIER_LTYPE==0)
				  OutPen -= (x+y) * OutlierLambdaInv;
				else
				  OutPen -= fabs(x-y) * OutlierLambdaInv;
				RFLOAT iscore2 = OutlierBias + Pen*OUTLIER_TYPE1 + Bias * biasWToutlierF + OutPen + PenSm;

				if(!(fabs(iscore2 - deltascore) < SCORE_MARGIN * fabs(iscore2) + ASCORE_MARGIN)){
				  printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d,rev=%d:I=%d,K=%d,J=%d, G=%d,T=%d,H=%d:x= %0.4f, y= %0.4f (DELTA_X=%d,DELTA_Y=%d,outlierExtend=%d,%d)\n",
	                               refid,rmap->id, mapid, nanomap->id, orientation, rev, I,K,J,G,T,H,x,y, DELTA_X,DELTA_Y,outlierExtend,outlierExtendLimit);
				  printf("\t deltascore= %0.8f: OutlierPenalty3= %0.8f, OutPenExtend(x,y,J-H,I-K-G,I,K)= %0.8f\n",
				       deltascore,OutlierPenalty3, OutPenExtend(x,y,J-H,I-K-G,I,K,PRtabY,refid));
				  printf("\t iscore2= %0.8f: Bias=%0.8f,Pen=%0.8f,Gauss=%0.8f,PenSm=%0.8f; OutPen=%0.8f,OutlierBias=%0.8f,biasWToutlierF=%0.8f\n",
				       iscore2,Bias,Pen,Gauss,PenSm,OutPen,OutlierBias,biasWToutlierF);
				  fflush(stdout);
				  assert(fabs(iscore2 - deltascore) < SCORE_MARGIN * fabs(iscore2) + ASCORE_MARGIN);
			        }
			      }
			      RFLOAT newscore = AscoreGT[H] + deltascore;
			      if(PVERB && orientation == OR_TRACE && I== I_TRACE && K== K_TRACE && J== J_TRACE && newscore > score){
                                #pragma omp critical
				{
				  printf("   testing G=%d,T=%d,H=%d,newscore=%0.6f(%0.6f + %0.6f)%c best=%0.6f\n",G,T,H,newscore,AscoreGT[H],
				       OutlierPenalty3 + OutPenExtend(X[J]-X[H],Yc(Y,I,K)-Yc(Y,G,T), J-H, I-K-G, I, K, PRtabY, refid), (newscore > score) ? '!' : ' ',score);
				  fflush(stdout);
			        }
			      }
			      if(newscore > score){
				score = newscore;
				bestG = G;
				bestH = H;
				bestT = T;
			      }
			    }
			  }
		        }
		      }
		      if(DEBUG>=2 && !((bestG <= origbestG || bestH <= origbestH) && score >= origscore)){
		        #pragma omp critical
			{
			  printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d:I=%d,K=%d,J=%d:origbestG=%d,origbestH=%d,bestG=%d,bestH=%d,score= %0.6f->%0.6f\n",
			       refid,rmap->id,mapid,nanomap->id,orientation,I,K,J,origbestG,origbestH,bestG,bestH,origscore,score);
			  fflush(stdout);
			  assert((bestG <= origbestG || bestH <= origbestH) && score >= origscore);
		        }
		      }

		      /* also check if best interval before (G,T,H) can be merged with this outlier interval : include checking G == I-K (if H <= min(J-1,origbestH)) and H==J (if G <= min(I-K-1,origbestG)) */
		      for(int G = I-K; G >= Gmin2; G--){
			if(G < I-K) G = min(origbestG,G);
			int Hmin2 = max(JMIN[G], origbestH - outlierExtend);
			int Hmax2 = min(JMAX[G],(G <= (min(I-K-1,origbestG))) ? J : min(J-1,origbestH));
			for(int T = Kmax[G]; T >= 0; T--){
			  if(DEBUG>=2) assert(G-T >= 1 );
			  for(int H = Hmax2; H >= Hmin2; H--){
			    if(H < J) H = min(origbestH,H);
			    if(DEBUG>=2) assert(H >= (max(1,JMIN[G])));
			    if(DEBUG>=2) assert(H <= (min(J,JMAX[G])));
			    if(DEBUG>=2) assert(G <= (min(I-K-1,origbestG)) || H <= (min(J-1,origbestH)));
			    int G2 = A->G(G,T,H);
			    if(G2 > 0){
			      int T2 = A->T(G,T,H);
			      int H2 = A->H(G,T,H);
			      if(DEBUG>=2) assert(G2-T2 >= 1 );		
			      if(DEBUG>=2) assert(H2 >= (max(1,JMIN[G2])));
			      if(DEBUG>=2) assert(H2 <= (min(H,JMAX[G2])));
			      if(!outlierExtendLimit || (max(I-K-G2,J-H2)) <= outlierExtendLimit){
				RFLOAT x = X[J] - X[H2], y = Yik - Yc(Y,G2,T2);
				if(OUTLIER_DELTA(x-y)){ // OutlierPenalty3 == OutlierPenalty + OutlierBias + PRtabY[K][I].Sm
				  RFLOAT deltascore = OutlierPenalty3 + OutPenExtend(x, y, J-H2, I-K-G2, I, K, PRtabY, refid) ; // OutlierPenaltyBC[I-K-G2] + OutlierPenaltyBC[J-H2];
				  if(DEBUG>=2){
				    RFLOAT Bias,Pen,Gauss,PenSm;
				    SintDetail(x,y, J-H2,I-K-G2,J,I,K,T2,PRtabY,Y,Bias,Pen,Gauss,PenSm,0);

				    RFLOAT OutPen = OutlierPenalty;
				    if(outlierBC)
				      OutPen += OutlierPenaltyBC[J-H2] + OutlierPenaltyBC[I-K-G2];
				    if(OUTLIER_LTYPE==0)
				      OutPen -= (x+y) * OutlierLambdaInv;
				    else
				      OutPen -= fabs(x-y) * OutlierLambdaInv;
				    RFLOAT iscore2 = OutlierBias + Pen*OUTLIER_TYPE1 + Bias * biasWToutlierF + OutPen + PenSm;
				    if(!(fabs(iscore2 - deltascore) < SCORE_MARGIN * fabs(iscore2) + ASCORE_MARGIN)){
				      printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d,rev=%d:I=%d,K=%d,J=%d, G2=%d,T2=%d,H2=%d:x= %0.4f, y= %0.4f (DELTA_X=%d,DELTA_Y=%d,outlierExtend=%d,%d)\n",
					     refid,rmap->id, mapid, nanomap->id, orientation, rev, I,K,J,G2,T2,H2,x,y, DELTA_X,DELTA_Y,outlierExtend,outlierExtendLimit);
				      printf("\t deltascore= %0.8f: OutlierPenalty3= %0.8f, OutPenExtend(x,y,J-H2,I-K-G2,I,K)= %0.8f\n",
					     deltascore,OutlierPenalty3, OutPenExtend(x,y,J-H2,I-K-G2,I,K,PRtabY,refid));
				      printf("\t iscore2= %0.8f: Bias=%0.8f,Pen=%0.8f,Gauss=%0.8f,PenSm=%0.8f; OutPen=%0.8f,OutlierBias=%0.8f,biasWToutlierF=%0.8f\n",
					     iscore2,Bias,Pen,Gauss,PenSm,OutPen,OutlierBias,biasWToutlierF);
				      fflush(stdout);
				      assert(fabs(iscore2 - deltascore) < SCORE_MARGIN * fabs(iscore2) + ASCORE_MARGIN);
				    }
			          }
				  RFLOAT newscore = A->score(G2,T2,H2) + deltascore;
				  if(PVERB && orientation == OR_TRACE && I== I_TRACE && K== K_TRACE && J== J_TRACE && newscore > score + 1e-6){
				    #pragma omp critical
				    {
				      printf("   testing G=%d,T=%d,H=%d: G2=%d,T2=%d,H2=%d,score=%0.6f(%0.6f + %0.6f)%c best=%0.6f\n",G,T,H,G2,T2,H2,newscore,A->score(G2,T2,H2),
					   OutlierPenalty3 + OutPenExtend(X[J]-X[H2],Yc(Y,I,K)-Yc(Y,G2,T2), J-H2, I-K-G2, I, K, PRtabY, refid), (newscore > score) ? '!' : ' ',score);
				      fflush(stdout);
				    }
				  }
				  if(newscore > score + 1e-6){
				    score = newscore;
				    bestG = G2;
				    bestH = H2;
				    bestT = T2;
				  }
			        }
			      }
			    }
			  }
		        }
		      }

		      if(origbestG == bestG && origbestH == bestH)
			break;

		      if(DEBUG>=2) assert(bestG <= origbestG || bestH <= origbestH);

		      if(PVERB && outlierExtend && I== I_TRACE && K== K_TRACE && J== J_TRACE){
		        #pragma omp critical
			{
			  printf("   cnt=%d:I=%d,K=%d,J=%d:best G=%d->%d,H=%d->%d,T=%d,score=%0.6f\n",
			       cnt, I,K,J,origbestG,bestG,origbestH, bestH,bestT,score);
			  fflush(stdout);
			}
		      }
		    }/* for(int cnt = 0; cnt < 100; cnt++) */
		  } else {/* try to merge non-outlier interval with previous interval (in case previous interval is an outlier) */
		    int G = A->G(bestG,bestT,bestH);
		    if(G > 0){
		      int T = A->T(bestG,bestT,bestH);
		      int H = A->H(bestG,bestT,bestH);
		      if(DEBUG>=2) assert(G-T >= 1 );		
		      if(DEBUG>=2) assert(H >= (max(1,JMIN[G])));
		      if(DEBUG>=2) assert(H <= (min(bestH,JMAX[G])));
		      if(!outlierExtendLimit || (max(I-K-G,J-H)) <= outlierExtendLimit){
			RFLOAT x = X[J] - X[H], y = Yik - Yc(Y,G,T);
			if(OUTLIER_DELTA(x-y)){
			  RFLOAT deltascore = OutlierPenalty3 + OutPenExtend(x, y, J-H, I-K-G, I, K, PRtabY, refid);// OutlierPenaltyBC[I-K-G] + OutlierPenaltyBC[J-H];
			  if(DEBUG>=2){
			    RFLOAT Bias,Pen,Gauss,PenSm;
			    SintDetail(x,y, J-H,I-K-G,J,I,K,T,PRtabY,Y,Bias,Pen,Gauss,PenSm,0);

			    RFLOAT OutPen = OutlierPenalty;
			    if(outlierBC)
			      OutPen += OutlierPenaltyBC[J-H] + OutlierPenaltyBC[I-K-G];
			    if(OUTLIER_LTYPE==0)
			      OutPen -= (x+y) * OutlierLambdaInv;
			    else
			      OutPen -= fabs(x-y) * OutlierLambdaInv;
			    RFLOAT iscore2 = OutlierBias + Pen*OUTLIER_TYPE1 + Bias * biasWToutlierF + OutPen + PenSm;
			    if(!(fabs(iscore2 - deltascore) < SCORE_MARGIN * fabs(iscore2) + ASCORE_MARGIN)){
			      printf("refid=%d(id=%lld),mapid=%d(id=%lld),or=%d,rev=%d:I=%d,K=%d,J=%d, G=%d,T=%d,H=%d:x= %0.4f, y= %0.4f (DELTA_X=%d,DELTA_Y=%d,outlierExtend=%d,%d)\n",
				  refid,rmap->id, mapid, nanomap->id, orientation, rev, I,K,J,G,T,H,x,y, DELTA_X,DELTA_Y,outlierExtend,outlierExtendLimit);
			      printf("\t deltascore= %0.8f: OutlierPenalty3= %0.8f, OutPenExtend(x,y,J-H,I-K-G,I,K)= %0.8f\n",
				     deltascore,OutlierPenalty3, OutPenExtend(x,y,J-H,I-K-G,I,K,PRtabY,refid));
			      printf("\t iscore2= %0.8f: Bias=%0.8f,Pen=%0.8f,Gauss=%0.8f,PenSm=%0.8f; OutPen=%0.8f,OutlierBias=%0.8f,biasWToutlierF=%0.8f\n",
				     iscore2,Bias,Pen,Gauss,PenSm,OutPen,OutlierBias,biasWToutlierF);
			      fflush(stdout);
			      assert(fabs(iscore2 - deltascore) < SCORE_MARGIN * fabs(iscore2) + ASCORE_MARGIN);
			    }
			  }
			  RFLOAT newscore = A->score(G,T,H) + deltascore;
			  if(newscore > score){
			    if(PVERB && orientation == OR_TRACE && I== I_TRACE && K== K_TRACE && J== J_TRACE){
			      #pragma omp critical
			      {
				printf("   best G=%d,T=%d,H=%d: G2=%d,T2=%d,H2=%d,score=%0.6f(%0.6f + %0.6f)\n",bestG,bestT,bestH,G,T,H,newscore,A->score(G,T,H),
				     OutlierPenalty3 + OutPenExtend(X[J]-X[H], Yc(Y,I,K)-Yc(Y,G,T), J-H, I-K-G, I, K, PRtabY, refid));
				fflush(stdout);
			      }
			    }
			    score = newscore;
			    bestG = G;
			    bestH = H;
			    bestT = T;
			  }
		        }
		      }
		    }
		  }
	        }// if(bestG >= 1)

		if(PVERB && (I==I_TRACE && K== K_TRACE && J== J_TRACE)){
	          #pragma omp critical
		  {
		    if(MultiMatches){
		      if(bestG > 0)
		        printf("refid=%d(%lld),mapid=%d(%lld),or=%d:I=%d,K=%d,J=%d(IL=%d->%d,KL=%d->%d,JL=%d->%d):best G=%d,H=%d,T=%d, A->score(I,K,J)=%0.6f -> %0.6f\n",
			     refid,rmap->id, mapid, nanomap->id, orientation, I,K,J,A->IL(I,K,J),A->IL(bestG,bestT,bestH),A->KL(I,K,J),A->KL(bestG,bestT,bestH),A->JL(I,K,J),A->JL(bestG,bestT,bestH),bestG,bestH,bestT,A->score(I,K,J),score);
		      else
			printf("refid=%d(%lld),mapid=%d(%lld),or=%d:I=%d,K=%d,J=%d(IL=%d,KL=%d,JL=%d):best G=%d,H=%d,T=%d, A->score(I,K,J)=%0.6f -> %0.6f\n",
			     refid,rmap->id, mapid, nanomap->id, orientation, I,K,J,A->IL(I,K,J),A->KL(I,K,J),A->JL(I,K,J),bestG,bestH,bestT,A->score(I,K,J),score);
		    } else
		      printf("refid=%d(%lld),mapid=%d(%lld),or=%d:I=%d,K=%d,J=%d:best G=%d,H=%d,T=%d, A->score(I,K,J)=%0.6f -> %0.6f\n",
				  refid,rmap->id, mapid, nanomap->id, orientation, I,K,J,bestG,bestH,bestT,A->score(I,K,J),score);
		    fflush(stdout);
		  }
	        }
	      
		if(DEBUG>=2 && !isfinite(score)){
		  printf("refid=%d,mapid=%d:I=%d,K=%d,J=%d:score=%e\n",refid,mapid,I,K,J,score);
		  fflush(stdout);
		  assert(isfinite(score));
		}

		A->score(I,K,J) = score;
		A->G(I,K,J) = bestG;
		A->T(I,K,J) = bestT;
		A->H(I,K,J) = bestH;
	      }// if(outlierExtend)

	      if(MultiMatches && bestG > 0){/* propagate IL,KL,JL values */
		if(DEBUG>=2) assert(0 <= bestT && bestT <= Kmax[bestG]);
		if(DEBUG>=2) assert(JMIN[bestG] <= bestH && bestH <= JMAX[bestG]);
		A->IL(I,K,J) = A->IL(bestG,bestT,bestH);
		A->KL(I,K,J) = A->KL(bestG,bestT,bestH);
		A->JL(I,K,J) = A->JL(bestG,bestT,bestH);
	      }
	    }//  J = jmin .. jmax
	  }// if(outlierExtend || MultiMatches)
        } // DEFER_BP && !outlierExtend
      }// K = Kmax[I] .. 0
    } // I = IMIN .. IMAX
  } // if(1)
}

#if USE_MIC
/* Vladimir's vectorized code of 3D recurrance */
static void refalign_3Drecurrance(FLOAT *Y, int N, FLOAT *X, int M,
				  int IMIN, int IMAX, int *JMIN, int *JMAX, int *Kmax, 
				  AL *A, int refid, int mapid, int orientation, int scaleID, int rev,
				  Cmap *rmap, Cmap *nanomap, CYPen ****YPenR, CXPen *XPen, Cprtab **PRtabY,
				  int tid, int numthreads)
{
  if(1){
    if(PVERB){
      printf("refid=%d(%lld),mapid=%d(%lld),or=%d,rev=%d,sc=%d:IMIN=%d,IMAX=%d: vectorized case\n",
	     refid,rmap->id, mapid, nanomap->id, orientation,rev,scaleID,IMIN,IMAX);
      fflush(stdout);
    }
    int I = IMIN;
    int imax = min(IMAX,IMIN + KMAX + DELTA_Y - 1);/* first few iterations of I */
    
    if(maptype) {
#define VAR_BLOCK(m)					\
      RFLOAT *deltaX##m = XPen->deltaX[m];		\
      RFLOAT *Pen##m = XPen->Pen[m];			\
      RFLOAT *Bias##m = XPen->Bias[m];	\
      RFLOAT *PenBias##m = XPen->PenBias[m];
    
      if(DEBUG) assert(DELTA_X == 4);
      VAR_BLOCK(1);
      VAR_BLOCK(2);
      VAR_BLOCK(3);
      VAR_BLOCK(4);

#undef VAR_BLOCK
    
      for(; I <= imax; I++){/* first few iterations of I */
	int jmin = JMIN[I];
	int jmax = JMAX[I];
	for(int K= Kmax[I]; K >= 0; K--){
	  RFLOAT *AscoreIK = &A->score(I,K,0);
	  int Gmin = max(IMIN,I-K- DELTA_Y);
	  for(int G= I-K; --G >= Gmin;){
	    for(int T = Kmax[G]; T >= 0; T--){
	      if(DEBUG>=2) assert(G-T >= 1);
	      RFLOAT deltaY,Ivar,GaussY,Sm;
	      RFLOAT penY = SintY(I-K-G,I,K,T,deltaY,Ivar,GaussY,Sm,YPenR);
	      RFLOAT *AscoreGT = &A->score(G,T,0);
#pragma ivdep
	      for(int J= jmin; J <= jmax; J++){
		RFLOAT score = AscoreIK[J];
		int Hmin = max(JMIN[G],J - DELTA_X);
		int Hmax = min(JMAX[G],J - 1);
		if (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE){
		  #pragma omp critical
		  {
		    printf("A1:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d..%d,IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d\n",
			   refid,mapid,orientation,I,K,J,G,T,Hmin,Hmax,IMIN,IMAX,JMIN[I],JMAX[I],JMIN[G],JMAX[G]);
		    fflush(stdout);
		  }
		}

#define BLOCK(m) {\
		  int H = J-m;						\
		  if(H >= Hmin && (!DIAGONAL || H <= Hmax)) {		\
		    if(DEBUG>=2) assert(H >= JMIN[G] && H <= JMAX[G]);	\
		    RFLOAT newscore = AscoreGT[H] + SintX_maptype1(deltaX##m[J],deltaY,penY,GaussY,Sm,Ivar, Pen##m[J], Bias##m[J], PenBias##m[J]); \
		    if (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE){ \
		      printf("A1:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d,IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d:AscoreGT[H]=%0.8f,newscore=%0.8f,score=%0.8f -> %0.8f\n", \
			     refid,mapid,orientation,I,K,J,G,T,H,IMIN,IMAX,JMIN[I],JMAX[I],JMIN[G],JMAX[G],AscoreGT[H],newscore,score,max(score,newscore)); \
		      fflush(stdout);					\
		    }							\
		    score = max(score,newscore);			\
		  }							\
		}
		if(DEBUG>=2) assert(DELTA_X == 4);
		BLOCK(4);
		BLOCK(3);
		BLOCK(2);
		BLOCK(1);
#undef BLOCK
		AscoreIK[J] = score;
	      }
	    }
	  }
	}
      }
    } else {// maptype == 0
 
#define VAR_BLOCK(m) \
      RFLOAT *deltaX##m = XPen->deltaX[m];	\
      RFLOAT *Pen##m = XPen->Pen[m];		   \
      RFLOAT *Bias##m = XPen->Bias[m];	\
      RFLOAT *PenBias##m = XPen->PenBias[m];
    
      VAR_BLOCK(1);
      VAR_BLOCK(2);
      VAR_BLOCK(3);
      VAR_BLOCK(4);

#undef VAR_BLOCK
    
      for(; I <= imax; I++){/* first few iterations of I */
	int jmin = max(2,JMIN[I]);
	int jmax = min(JMIN[I]+3, JMAX[I]);/* first few iterations of J = jmin..jmax */
	int jmax3 = JMAX[I];/* last few iteration J = jmax2+1 .. jmax3 : only needed if DIAGONAL */
	for(int K= Kmax[I]; K >= 0; K--){
	  RFLOAT *AscoreIK = &A->score(I,K,0);
	  int Gmin = max(IMIN,I - K - DELTA_Y);

	  for(int G= I-K; --G >= Gmin;){
#if (USE_MIC && USE_RFLOAT)
	    int jmax2 = min(jmax3,1+JMAX[G]);/* iterations in vector mode J = jmax+1..jmax2 (if USE_MIC) */
#endif
	    for(int T = Kmax[G]; T >= 0; T--){
	      if(DEBUG>=2) assert(G-T >= 1);
	      RFLOAT deltaY,Ivar,GaussY,Sm;
	      RFLOAT penY = SintY(I-K-G,I,K,T,deltaY,Ivar,GaussY,Sm,YPenR);
	      RFLOAT *AscoreGT = &A->score(G,T,0);
	      int J= jmin;

#pragma ivdep
	      for(; J <= jmax; J++){
		RFLOAT score = AscoreIK[J];
		int Hmin = max(JMIN[G],J-4);
		int Hmax = min(JMAX[G],J-1);

#if 0 // original code
		for(int H = Hmax; H >= Hmin; H-){
		  int m = J-H;
		  RFLOAT deltaX = XPen->deltaX[m][J];
		  RFLOAT newscore = AscoreGT[H] + SintX_maptype0(deltaX,deltaY,penY,GaussY,Sm,Ivar,XPen->Pen[m][J],XPen->OutlierBias[m][J],XPen->PenBias[m][J],0);
		  score = max(score,newscore);
		}
#else // new unrolled code

#define BLOCK(m) {\
		  int H = J-m;			\
		  if(H >= Hmin && (!DIAGONAL || H <= Hmax)) {		\
		    RFLOAT newscore = AscoreGT[H] + SintX_maptype0(deltaX##m[J],deltaY,penY,GaussY,Sm,Ivar, Pen##m[J], Bias##m[J], PenBias##m[J],0); \
		    if (PVERB>=2 && I==I_TRACE && K==K_TRACE && J==J_TRACE){ \
		      printf("A:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d,IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d:AscoreGT[H]=%0.8f,newscore=%0.8f,score=%0.8f -> %0.8f\n", \
			     refid,mapid,orientation,I,K,J,G,T,H,IMIN,IMAX,JMIN[I],JMAX[I],JMIN[G],JMAX[G],AscoreGT[H],newscore,score,max(score,newscore)); \
		      fflush(stdout);					\
		    }							\
		    score = max(score,newscore);			\
		  }							\
		}
	        if(DIAGONAL) BLOCK(4); // NOTE if(DIAGONAL==0) : J <= 4, hence H <= J-1..J-3 so BLOCK(4) would never execute
		BLOCK(3);
		BLOCK(2);
		BLOCK(1);
#undef BLOCK
#endif // new unrolled code

		if(PVERB>=2 && I==I_TRACE && K==K_TRACE && J==J_TRACE){
		  printf("8: A->score(I=%d,K=%d,J=%d)=%0.6f -> %0.6f\n",I,K,J,A->score(I,K,J),score);
		  fflush(stdout);
		}
		AscoreIK[J] = score;
	      }

	      /* iterations in vector mode J = jmax+1..jmax2 */

#if (USE_MIC && USE_RFLOAT) // A 
	      __m512 v_deltaY = _mm512_set1_ps(deltaY);
	      __m512 v_Ivar = _mm512_set1_ps(Ivar);
	      __m512 v_penY = _mm512_set1_ps(penY);
	      __m512 v_GaussY = _mm512_set1_ps(GaussY);
	      __m512 v_Sm = _mm512_set1_ps(Sm);

	      for(; J <= jmax2; J += 16){/* critical J loop (with maptype==0, first few iterations of I) */
		__mmask16 mask = (J <= jmax2-15) ? 0xffff : ((1u << (jmax2-J+1))-1);
		__m512 v_score = _mm512_loadu_ps(&AscoreIK[J], rs_heap->thread_block_start(tid), rs_heap->thread_block_siz(), tid, 1, I, J, K);
#define	BLOCK(m){							\
		  if(DEBUG>=2) assert(JMIN[G] <= J-m && min(J+15,jmax2)-m <= JMAX[G]); \
		  __m512 v_prev_score = _mm512_loadu_ps(&AscoreGT[J-m], rs_heap->thread_block_start(tid), rs_heap->thread_block_siz(), tid, 2); \
		  __m512 v_deltaX = _mm512_loadu_ps(&deltaX##m[J], XPen->mem_pool, XPen->mem_pool_size, tid, 3); \
		  __m512 v_XPenJM_Bias = _mm512_loadu_ps(&Bias##m[J], XPen->mem_pool, XPen->mem_pool_size, tid, 4); \
		  __m512 v_XPenJM_Pen = _mm512_loadu_ps(&Pen##m[J], XPen->mem_pool, XPen->mem_pool_size, tid, 5); \
		  __m512 v_XPenJM_PenBias = _mm512_loadu_ps(&PenBias##m[J], XPen->mem_pool, XPen->mem_pool_size, tid, 6); \
		  __m512 v_score_orig = v_score;	\
		  __m512 v_newscore = _mm512_add_ps(v_prev_score, SintX_mm512_maptype0(v_deltaX,v_deltaY, v_penY, v_GaussY, v_Sm, v_Ivar, v_XPenJM_Bias, v_XPenJM_Pen, v_XPenJM_PenBias,0,NULL)); \
		  v_score = _mm512_mask_gmax_ps(v_score, mask, v_score, v_newscore); \
		  if(DEBUG>=2 || (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE)){ \
		    RFLOAT nscore[16], origscore[16], newscore[16];	\
		    _mm512_maskstoreu_ps(origscore, mask, v_score_orig); \
		    _mm512_maskstoreu_ps(newscore, mask, v_newscore);	\
		    _mm512_maskstoreu_ps(nscore, mask, v_score);	\
		    for(int j = J; j <= min(jmax2,J+15); j++)		\
		      if ((DEBUG && !(isfinite(nscore[j-J]) && isfinite(origscore[j-J]))) || (TRACE && rmap->id==REF_TRACE && nanomap->id==MAP_TRACE && orientation==OR_TRACE && I==I_TRACE && K==K_TRACE && j==J_TRACE)){ \
			printf("P:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d,IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d,jmax2=%d,mask=x%04x:AscoreGT[H]=%0.8f,newscore=%0.8f,score=%0.8e -> %0.8f\n", \
			       refid,mapid,orientation,I,K,j,G,T,j-m,IMIN,IMAX,JMIN[I],JMAX[I],JMIN[G],JMAX[G],jmax2,mask,AscoreGT[j-m],newscore[j-J],origscore[j-J],nscore[j-J]); \
			fflush(stdout);					\
			assert(isfinite(nscore[j-J]) && isfinite(origscore[j-J])); \
		      }							\
		  }							\
		}
		BLOCK(1);
		BLOCK(2);
		BLOCK(3);
		BLOCK(4);
#undef BLOCK
		if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && K <= A->kmax);
		_mm512_maskstoreu_ps(&AscoreIK[J], mask, v_score);
	      }
	      if(DIAGONAL) 
		J = max(jmax,jmax2)+1;
	      else if(DEBUG>=2) assert(max(jmax,jmax2) == jmax3);

#endif // if (USE_MIC && USE_RFLOAT)
	    
#pragma ivdep
	      for(; J <= jmax3; J++){
		int m_min = J - JMAX[G];
		RFLOAT score = AscoreIK[J];
#define BLOCK(m){							\
		   if(!DIAGONAL || m >= m_min){ \
		     if(DEBUG>=2) assert(J-m <= JMAX[G] && J-m >= JMIN[G]);	\
		     RFLOAT newscore = AscoreGT[J-m] + SintX_maptype0(deltaX##m[J],deltaY,penY,GaussY,Sm,Ivar, Pen##m[J], Bias##m[J], PenBias##m[J],0); \
		     if (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE){ \
		       printf("B:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d,IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d:AscoreGT[H]=%0.8f,newscore=%0.8f,score=%0.8f -> %0.8f\n", \
			      refid,mapid,orientation,I,K,J,G,T,J-m,IMIN,IMAX,JMIN[I],JMAX[I],JMIN[G],JMAX[G],AscoreGT[J-m],newscore,score,max(score,newscore)); \
		       fflush(stdout);					\
		     }							\
		     score = max(score,newscore);			\
		   } \
		}
		BLOCK(4);
		BLOCK(3);
		BLOCK(2);
		BLOCK(1);
#undef BLOCK
		if(DEBUG>=2 && !isfinite(score)){
		  #pragma omp critical
		  {
		    printf("refid=%lld,mapid=%lld,orientation=%d:N=%d,M=%d:I=%d,K=%d,J=%d,JMIN[I]=%d,JMAX[I]=%d,G=%d,T=%d,JMIN[G]=%d,JMAX[G]=%d:A->score(I,K,J)=%0.8e,score=%0.8e\n",
			   rmap->id,nanomap->id,orientation,N,M,I,K,J,JMIN[I],JMAX[I],G,T,JMIN[G],JMAX[G],A->score(I,K,J),score);
		    fflush(stdout);
		    assert(isfinite(score));
		  }
		}
		if(PVERB>=2 && I==I_TRACE && K==K_TRACE && J==J_TRACE){
		  printf("9: A->score(I=%d,K=%d,J=%d)=%0.6f -> %0.6f\n",I,K,J,A->score(I,K,J),score);
		  fflush(stdout);
		}
		AscoreIK[J] = score;
	      }
	    }
	  }
	}
      }
    } // maptype==0

    if(PVERB>=2 && IMIN <= I_TRACE && I_TRACE <= IMAX && JMIN[I_TRACE] <= J_TRACE && J_TRACE <= JMAX[I_TRACE]){
      printf("10: refid=%d(%lld),mapid=%d(%lld),or=%d:A->score(I=%d,K=%d,J=%d)=%0.6f\n",refid,rmap->id,mapid,nanomap->id,orientation,I_TRACE,K_TRACE,J_TRACE,A->score(I_TRACE,K_TRACE,J_TRACE));
      fflush(stdout);
    }

#define VAR_BLOCK(m) \
    RFLOAT *deltaX##m = XPen->deltaX[m];	\
    RFLOAT *Pen##m = XPen->Pen[m];		   \
    RFLOAT *Bias##m = XPen->Bias[m]; \
    RFLOAT *PenBias##m = XPen->PenBias[m];
    
    VAR_BLOCK(1);
    VAR_BLOCK(2);
    VAR_BLOCK(3);
    VAR_BLOCK(4);

#undef VAR_BLOCK

    /* most iterations of I : */
    if(maptype){/* maptype == 1 : no conditional code in inner loop */
      for(; I <= IMAX; I++){
	int jmin = max(2,JMIN[I]);
	int jmax = min(JMIN[I]+3,JMAX[I]);/* first few iterations of J = jmin..jmax */
	for(int K= Kmax[I]; K >= 0; K--){
          CYPen *p_base = &YPenR[K][0][I][I-K];
	  RFLOAT *AscoreIK = &A->score(I,K,0);

	  int Gmin = max(IMIN,I - K - DELTA_Y);
	  for(int G= I-K; --G >= Gmin;){
	    // int Hmin = JMIN[G]; // NOTE : long burried bug when Hmin exceeds J - DELTA_X : moved into J loop to fix
	    int J = jmin;
	    for(; J <= jmax; J++){/* first few iterations of J (maptype==1, most iterations of I) */
	      int Hmin = max(JMIN[G], J - DELTA_X);
	      int Hmax = min(JMAX[G], J-1);
	      if(Hmin > Hmax)
		continue;

	      if(DEBUG>=2 && !(Hmin <= Hmax && Hmax-Hmin < DELTA_X && DELTA_X == 4)){
		#pragma omp critical
		{
		  printf("refid=%d(%lld),mapid=%d(%lld),or=%d:I=%d,K=%d,J=%d,G=%d:Hmin=%d,Hmax=%d,DELTA_X=%d(JMIN[G]=%d,JMAX[G]=%d)\n",
			 refid,rmap->id,mapid,nanomap->id,orientation,I,K,J,G,Hmin,Hmax,DELTA_X,JMIN[G],JMAX[G]);
		  fflush(stdout);
		  assert(Hmin <= Hmax && Hmax-Hmin < DELTA_X && DELTA_X == 4);
		}
	      }

	      RFLOAT score = AscoreIK[J];

	      for(int T = Kmax[G]; T >= 0; T--){
		RFLOAT * AscoreGT = &A->score(G,T,0);
		//		CYPen *p = &YPenR[K][T][I][I-K-G];
		CYPen *p= &p_base[-T*N*DELTA_Y-G];
		if(DEBUG>=2 && ! (p == &YPenR[K][T][I][I-K-G])){
		  printf("*** INTERNAL ERROR: Pointers don't match p=%p p_expected=%p p_base=%p K=%d T=%d I=%d G=%d N=%d DELTA_Y=%d\n", p, &YPenR[K][T][I][I-K-G], p_base, K, T, I, G, N, DELTA_Y);
		  fflush(stdout);
		  assert(p == &YPenR[K][T][I][I-K-G]);
		}
		_mm_prefetch((const char *)&AscoreGT[J-4], _MM_HINT_T0);
		RFLOAT deltaY = p->deltaY;
		RFLOAT Ivar = p->Ivar;
		RFLOAT penY = p->Pen;
		RFLOAT GaussY = p->Gauss;
		RFLOAT Sm = p->Sm;

#define BLOCK(m) {	   \
		  int H = J-(m);					\
		  if(H >= Hmin && (!DIAGONAL || H <= Hmax)) {		\
		    if(DEBUG>=2) assert(JMIN[G] <= H && H <= JMAX[G]);	\
		    RFLOAT newscore = AscoreGT[H] + SintX_maptype1(deltaX##m[J],deltaY,penY,GaussY,Sm,Ivar,Pen##m[J],Bias##m[J], PenBias##m[J]); \
		    if (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE){ \
		      printf("B1:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d,IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d:AscoreGT[H]=%0.8f,newscore=%0.8f,score=%0.8f -> %0.8f\n", \
			     refid,mapid,orientation,I,K,J,G,T,H,IMIN,IMAX,JMIN[I],JMAX[I],JMIN[G],JMAX[G],AscoreGT[H],newscore,score,max(score,newscore)); \
		      fflush(stdout);					\
		    }							\
		    score = max(score,newscore);			\
		  }							\
		}
		BLOCK(4);
		BLOCK(3);
		BLOCK(2);
		BLOCK(1);
#undef BLOCK

	      }
	      if(DEBUG>=2) assert(isfinite(score));
	      if(PVERB>=2 && I==I_TRACE && K==K_TRACE && J==J_TRACE){
		printf("11: A->score(I=%d,K=%d,J=%d)=%0.6f -> %0.6f\n",I,K,J,A->score(I,K,J),score);
		fflush(stdout);
	      }
	      AscoreIK[J] = score;
	    } /* first few iterations of J (maptype==1, most iterations of I) */

#if (USE_SSE==1 && USE_AVX==1 && USE_RFLOAT==1)
	    for(int jmax3 = min(JMAX[I],JMAX[G]+1) - 7; J <= jmax3; J += 8){/* critical J loop (with maptype==1, most iterations of I) with USE_SSE && USE_AVX && USE_RFLOAT */
	      __m256 v_score = _mm256_loadu_ps(&(AscoreIK[J]));
	      for(int T = Kmax[G]; T >= 0; T--){
		//		CYPen *p = &YPenR[K][T][I][I-K-G];
		CYPen *p= &p_base[-T*N*DELTA_Y-G];
		if(DEBUG>=2) assert(p == &YPenR[K][T][I][I-K-G]);

		__m256 v_deltaY = _mm256_set1_ps(p->deltaY);
		__m256 v_Ivar = _mm256_broadcast_ss(&p->Ivar);
		__m256 v_penY = _mm256_broadcast_ss(&p->Pen);
		__m256 v_GaussY = _mm256_broadcast_ss(&p->Gauss);
		__m256 v_Sm = _mm256_broadcast_ss(&p->Sm);
		RFLOAT *AscoreGT = &A->score(G,T,J);
		
		for(int m = 1; m <= 4;m++){/* NOTE : J-m corresponds to the H index in the original code and assumes DELTA_X == 4 */
		  if(DEBUG>=2) assert(J-m >= JMIN[G]);
		  if(DEBUG>=2) assert(J+7-m <= JMAX[G]);
		  __m256 v_deltaX = _mm256_loadu_ps(&XPen->deltaX[m][J]);
		  __m256 v_prev_score = _mm256_loadu_ps(&AscoreGT[-m]);
		  __m256 v_XPenJM_Bias = _mm256_loadu_ps(&XPen->Bias[m][J]);
		  __m256 v_XPenJM_Pen = _mm256_loadu_ps(&XPen->Pen[m][J]);
		  __m256 v_XPenJM_PenBias = _mm256_loadu_ps(&XPen->PenBias[m][J]);
		  __m256 v_newscore = _mm256_add_ps(v_prev_score, SintX_mm256_maptype1(v_deltaX,v_deltaY, v_penY, v_GaussY, v_Sm, v_Ivar, v_XPenJM_Bias, v_XPenJM_Pen, v_XPenJM_PenBias));
		  v_score = _mm256_max_ps(v_score, v_newscore);
		  if(PVERB && I==I_TRACE && K==K_TRACE && J <= J_TRACE && J_TRACE <= J+7){
		    int H = J-m;
		    float newscore[8], prevscore[8],score[8], deltaX[8], deltaY[8], OutlierBias[8], Pen[8], penY[8], GaussY[8], Sm[8], Ivar[8], OutlierPenalty[8], vOutlierLambdaInv[8], SintX[8];
		    __m256 v_SintX = SintX_mm256_maptype1(v_deltaX,v_deltaY, v_penY, v_GaussY, v_Sm, v_Ivar, v_XPenJM_Bias, v_XPenJM_Pen, v_XPenJM_PenBias);

		    _mm256_storeu_ps(newscore, v_newscore);
		    _mm256_storeu_ps(prevscore, v_prev_score);
		    _mm256_storeu_ps(score, v_score);
		    _mm256_storeu_ps(deltaY, v_deltaY);
		    _mm256_storeu_ps(penY, v_penY);
		    _mm256_storeu_ps(GaussY, v_GaussY);
		    _mm256_storeu_ps(Sm, v_Sm);
		    _mm256_storeu_ps(Ivar, v_Ivar);
		    _mm256_storeu_ps(deltaX, v_deltaX);
		    _mm256_storeu_ps(OutlierBias, v_XPenJM_Bias);
		    _mm256_storeu_ps(Pen, v_XPenJM_Pen);
		    _mm256_storeu_ps(OutlierPenalty, v256_OutlierPenalty);
		    _mm256_storeu_ps(vOutlierLambdaInv, v256_OutlierLambdaInv);
		    _mm256_storeu_ps(SintX, v_SintX);

		    for(int j = 0; j < 8; j++)
		      if(J+j == J_TRACE){
			printf("B2:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d,IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d:j=%d:AscoreIK[j]=%0.6f,prevscore[j]=%0.6f,newscore[j]=%0.6f(delta=%0.6f,SintX=%0.6f),score->%0.6f\n",
			       refid,mapid,orientation,I,K,J+j,G,T,H+j,IMIN,IMAX,JMIN[I],JMAX[I],JMIN[G],JMAX[G],j,AscoreIK[j],prevscore[j],newscore[j],newscore[j]-prevscore[j],SintX[j],score[j]);
			printf("    x=%0.4f,y=%0.4f,OutlierBias=%0.6f,Pen=%0.6f,GaussY=%0.6f\n", deltaX[j],deltaY[j],OutlierBias[j],Pen[j],GaussY[j]);
			printf("    penY=%0.6f,Ivar=%0.8f,Sm=%0.6f,OutlierPenalty=%0.6f,OutlierLambdaInv=%0.6f,%0.6f\n", penY[j], Ivar[j], Sm[j], OutlierPenalty[j], vOutlierLambdaInv[j],OutlierLambdaInv);
		      }		    
		    fflush(stdout);
		  }
		}
	      } /* T loop */
	      _mm256_storeu_ps(&AscoreIK[J], v_score);
	      if(DEBUG>=2) 
		for(int j = J; j < J+8; j++)
		  assert(isfinite(A->score(I,K,j)));
	    } /* critical J loop (with maptype==1, most iterations of I) with USE_SSE && USE_AVX && USE_RFLOAT */

#endif // USE_SSE==1 && USE_AVX==1 && USE_RFLOAT==1

#if (USE_MIC && USE_RFLOAT==1) // B
	    int jmax3 = min(JMAX[I],JMAX[G]+1);
	    for(; J <= jmax3; J += 16){/* critical J loop (with maptype==1, most iterations of I) with USE_MIC && USE_RFLOAT */
	      __mmask16 mask = (J <= jmax3-15) ? 0xffff : ((1u<<(jmax3-J+1))-1);
	      __m512 v_score = _mm512_maskloadu_ps(&AscoreIK[J], mask, rs_heap->thread_block_start(tid), rs_heap->thread_block_siz(), tid, 13);
	      for(int T = Kmax[G]; T >= 0; T--){
		//		CYPen *p = &YPenR[K][T][I][I-K-G];
		CYPen *p= &p_base[-T*N*DELTA_Y-G];
		if(DEBUG>=2 && ! (p == &YPenR[K][T][I][I-K-G])){
		  printf("*** INTERNAL ERROR: Pointers don't match p=%p p_expected=%p p_base=%p K=%d T=%d I=%d G=%d N=%d DELTA_Y=%d\n", p, &YPenR[K][T][I][I-K-G], p_base, K, T, I, G, N, DELTA_Y);
		  fflush(stdout);
		  assert(p == &YPenR[K][T][I][I-K-G]);
		}
		__m512 v_deltaY = _mm512_set1_ps(p->deltaY);
		__m512 v_Ivar = _mm512_set1_ps(p->Ivar);
		__m512 v_penY = _mm512_set1_ps(p->Pen);
		__m512 v_GaussY = _mm512_set1_ps(p->Gauss);
		__m512 v_Sm = _mm512_set1_ps(p->Sm);
		RFLOAT *AscoreGT = &A->score(G,T,J);
		
		for(int H = 1; H <= 4;H++){/* NOTE : J-H corresponds to H in the original code */
		  if(DEBUG>=2) assert(J-H >= JMIN[G]);
		  if(DEBUG>=2) assert(min(jmax3,J+15)-H <= JMAX[G]);
		  __m512 v_prev_score = _mm512_maskloadu_ps(&AscoreGT[-H], mask, rs_heap->thread_block_start(tid), rs_heap->thread_block_siz(), tid, 14);
		  __m512 v_deltaX = _mm512_maskloadu_ps(&XPen->deltaX[H][J], mask, XPen->mem_pool, XPen->mem_pool_size, tid, 15);
		  __m512 v_XPenJM_Bias = _mm512_maskloadu_ps(&XPen->Bias[H][J], mask, XPen->mem_pool, XPen->mem_pool_size, tid, 16);
		  __m512 v_XPenJM_Pen = _mm512_maskloadu_ps(&XPen->Pen[H][J], mask, XPen->mem_pool, XPen->mem_pool_size, tid, 17);
		  __m512 v_XPenJM_PenBias = _mm512_maskloadu_ps(&XPen->PenBias[H][J], mask, XPen->mem_pool, XPen->mem_pool_size, tid, 18);
		  v_score = _mm512_gmax_ps(v_score, _mm512_add_ps(v_prev_score, SintX_mm512_maptype1(v_deltaX,v_deltaY, v_penY, v_GaussY, v_Sm, v_Ivar, v_XPenJM_Bias, v_XPenJM_Pen, v_XPenJM_PenBias)));
		}
	      } /* T loop */
	      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && K <= A->kmax);
	      _mm512_maskstoreu_ps(&AscoreIK[J], mask, v_score);
	      if(PVERB>=2 && I==I_TRACE && K==K_TRACE && J <= J_TRACE && J_TRACE < J+16){
		printf("12: A->score(I=%d,K=%d,J=%d) -> %0.6f\n",I,K,J,A->score(I,K,J));
		fflush(stdout);
	      }

	      if(DEBUG>=2) 
		for(int j = J; j < J+16; j++)
		  assert(isfinite(A->score(I,K,j)));
	    } /* critical loop (with maptype==1, most iterations of I) with USE_MIC && USE_RFLOAT */
	    if(DIAGONAL)
	      J = max(jmax,jmax3)+1;
	    else if(DEBUG>=2) assert(max(jmax,jmax3) == JMAX[I]);
#endif // USE_MIC && USE_RFLOAT==1

	    for(int jmax3 = JMAX[I]; J <= jmax3; J++){/* remaining iterations of J (with maptype==1, most iterations of I) all use cases */
	      int m_min = max(1,J - JMAX[G]);
	      RFLOAT score = AscoreIK[J];
	      for(int T = Kmax[G]; T >= 0; T--){
		CYPen *p= &p_base[-T*N*DELTA_Y-G];
		if(DEBUG>=2 && ! (p == &YPenR[K][T][I][I-K-G])){
		  printf("*** INTERNAL ERROR: Pointers don't match p=%p p_expected=%p p_base=%p K=%d T=%d I=%d G=%d N=%d DELTA_Y=%d\n", p, &YPenR[K][T][I][I-K-G], p_base, K, T, I, G, N, DELTA_Y);
		  fflush(stdout);
		  assert(p == &YPenR[K][T][I][I-K-G]);
		}
		RFLOAT deltaY = p->deltaY;
		RFLOAT Ivar = p->Ivar;
		RFLOAT penY = p->Pen;
		RFLOAT GaussY = p->Gauss;
		RFLOAT Sm = p->Sm;

		RFLOAT *AscoreGT = &A->score(G,T,0);

#define BLOCK(m){							\
		  if(!DIAGONAL || m >= m_min) { \
		    if(DEBUG>=2) assert(J-m >= JMIN[G]);	\
		    if(DEBUG>=2) assert(J-m <= JMAX[G]);	\
		    RFLOAT newscore = AscoreGT[J-m] + SintX_maptype1(deltaX##m[J],deltaY,penY,GaussY,Sm,Ivar,Pen##m[J],Bias##m[J], PenBias##m[J]); \
		    if (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE){ \
		      printf("C1:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d,IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d:AscoreGT[H]=%0.8f,newscore=%0.8f,score=%0.8f -> %0.8f\n", \
			     refid,mapid,orientation,I,K,J,G,T,J-m,IMIN,IMAX,JMIN[I],JMAX[I],JMIN[G],JMAX[G],AscoreGT[J-m],newscore,score,max(score,newscore)); \
		      fflush(stdout);					\
		    }							\
		    score = max(score,newscore);			\
	          } \
		}
		BLOCK(4);
		BLOCK(3);
		BLOCK(2);
		BLOCK(1);
#undef BLOCK
	      } /* T loop */
	      if(DEBUG>=2) assert(isfinite(score));
	      if(PVERB>=2 && I==I_TRACE && K==K_TRACE && J <= J_TRACE && J_TRACE < J+16){
		printf("13: A->score(I=%d,K=%d,J=%d)=%0.6f -> %0.6f\n",I,K,J,A->score(I,K,J),score);
		fflush(stdout);
	      }
	      AscoreIK[J] = score;
	    } /* remaining iterations of J (with maptype==1, most iterations of I) all use cases */
	  } /* G loop */
	} /* K loop */
      } /* I loop */
    } else {/* maptype == 0 : handles internal outliers with conditional code in inner loop : OUTLIER(0,deltaX,deltaY) */

#define VAR_BLOCK(m) \
      RFLOAT *deltaX##m = XPen->deltaX[m];	\
      RFLOAT *Pen##m = XPen->Pen[m];		   \
      RFLOAT *Bias##m = XPen->Bias[m];	\
      RFLOAT *PenBias##m = XPen->PenBias[m];
      
      VAR_BLOCK(1);
      VAR_BLOCK(2);
      VAR_BLOCK(3);
      VAR_BLOCK(4);
   
#undef VAR_BLOCK
    
      RFLOAT *var_scores = (RFLOAT *)alloca(16*sizeof(RFLOAT));

#if (USE_MIC==1 && USE_RFLOAT==1)
      SintX_scratch_mm512 scratch;
      scratch.usage = 0;
      scratch.mask = 0;
#endif
      
      for(; I <= IMAX; I++){/* most iterations of I */
	int jmin = max(2,JMIN[I]);
	int jmax = min(JMIN[I]+3, JMAX[I]);/* first few iterations of J = jmin..jmax */
	for(int K= Kmax[I]; K >= 0; K--){
          CYPen *p_base = &YPenR[K][0][I][I-K];
	  RFLOAT *AscoreIK = &A->score(I,K,0);

	  int Gmin = max(IMIN,I - K - DELTA_Y);// WAS I-K-6
	  for(int G= Gmin; G < I-K; G++){
	    if(PVERB && I==I_TRACE && K==K_TRACE){
	      #pragma omp critical
	      {
		printf("refid=%d,mapid=%d,or=%d,M=%d,N=%d:I=%d,K=%d,Gmin=%d,G=%d,jmin=%d,jmax=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d (top of G loop)\n",
		       refid,mapid,orientation,M,N,I,K,Gmin,G,jmin,jmax,JMIN[I],JMAX[I],JMIN[G],JMAX[G]);
		fflush(stdout);
	      }
	    }
	    // int Hmin = JMIN[G]; // NOTE : long burried bug when Hmin exceeds J - DELTA_X : moved into J loop to fix
	    int J = jmin;
	    for(; J <= jmax; J++){/* first few iterations of J (maptype==0, most iterations of I) */
	      int Hmin = max(JMIN[G], J - DELTA_X);
	      int Hmax = min(JMAX[G], J-1);
	      if(Hmin > Hmax)
		continue;

	      if(DEBUG>=2 && !(Hmin <= Hmax && Hmax-Hmin < DELTA_X && DELTA_X == 4)){
		#pragma omp critical
		{
		  printf("refid=%d(%lld),mapid=%d(%lld),or=%d:I=%d,K=%d,J=%d,G=%d:Hmin=%d,Hmax=%d,DELTA_X=%d(JMIN[G]=%d,JMAX[G]=%d)\n",
			 refid,rmap->id,mapid,nanomap->id,orientation,I,K,J,G,Hmin,Hmax,DELTA_X,JMIN[G],JMAX[G]);
		  fflush(stdout);
		  assert(Hmin <= Hmax && Hmax-Hmin < DELTA_X && DELTA_X == 4);
		}
	      }

	      RFLOAT score = AscoreIK[J];
	    
#define VARSS_BLOCK(m) \
	      var_scores[4-m+VS_deltaX] = deltaX##m[J];	      \
	      var_scores[4-m+VS_Pen] = Pen##m[J];		      \
	      var_scores[4-m+VS_Bias] = Bias##m[J]; \
	      var_scores[4-m+VS_PenBias] = PenBias##m[J];
	    
	      VARSS_BLOCK(4);
	      VARSS_BLOCK(3);
	      VARSS_BLOCK(2);
	      VARSS_BLOCK(1);
#undef VARSS_BLOCK
	    
#if (USE_MIC && USE_RFLOAT==1) // C1 

#define SCRATCH_BLOCK(var) \
	      scratch.var = _mm512_setr4_ps(var ##4[J], var ##3[J], var ##2[J], var ##1[J]); // why the space after var ???

	      SCRATCH_BLOCK(deltaX);
	      SCRATCH_BLOCK(Pen);
	      SCRATCH_BLOCK(Bias);
	      SCRATCH_BLOCK(PenBias);

	      __m512 v_score = _mm512_set1_ps(score);
#endif // USE_MIC && USE_RFLOAT==1
	    
	      //for(int T = Kmax[G]; T >= 0; T--){
	      int T=0;
#pragma nounroll
	      do {
		RFLOAT * AscoreGT = &A->score(G,T,0);
		CYPen *p= &p_base[-T*N*DELTA_Y-G];
		if(DEBUG>=2) {
		  if(p != &YPenR[K][T][I][I-K-G]){
		    printf("*** INTERNAL ERROR: Pointers don't match p=%p p_expected=%p p_base=%p K=%d T=%d I=%d G=%d N=%d DELTA_Y=%d\n", p, &YPenR[K][T][I][I-K-G], p_base, K, T, I, G, N, DELTA_Y);
		    fflush(stdout);
		    assert(p == &YPenR[K][T][I][I-K-G]);
		  }
		}
		_mm_prefetch((const char *)&AscoreGT[Hmin-1], _MM_HINT_T0);
		
#if (USE_MIC==1 && USE_RFLOAT==1) // C2
		if(DEBUG>=2 && !(1 <= J - Hmax && J - Hmin <= 4)){
		  #pragma omp critical
		  {
		    printf("refid=%d(%lld),mapid=%d(%lld),or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d..%d,DELTA_X=%d\n",refid,rmap->id,mapid,nanomap->id,orientation,I,K,J,G,T,Hmin,Hmax,DELTA_X);
		    fflush(stdout);
		    assert(1 <= J - Hmax && J - Hmin <= 4);
		  }
		}
		v_score = SintX_scratch_feed_maptype0(scratch, v_score, p, AscoreGT, J, J - Hmax, J - Hmin, 
						      (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE) ? 1 : 0/* verbose*/, T);
#else // !(USE_MIC==1 && USE_RFLOAT==1)
		score = SintX_BLOCK4_maptype0(score, p, var_scores, AscoreGT, J, Hmin, Hmax, G, JMIN, JMAX);
#endif
		T++;
	      } while(T<=Kmax[G]);

#if (USE_MIC==1 && USE_RFLOAT==1) // C3
	      if(PVERB && I==I_TRACE && K == K_TRACE && J==J_TRACE){
		__m512 new_v_score = v_score;
		SintX_scratch_mm512 s = scratch;
		__mmask16 mask = scratch.mask;
		new_v_score = SintX_scratch_collapse_maptype0(s, new_v_score);
		float origscore[16], newscore[16], scoreGT[16], deltaX[16], deltaY[16], penY[16], Gauss[16],Sm[16], Ivar[16], Bias[16], Pen[16], PenBias[16], Sint[16], OutPen[16];
		for(int j = 0; j < 16; j++)
		  deltaX[j] = -1.0;
		_mm512_maskstoreu_ps(deltaX, mask, s.deltaX);

		__m512 v_Sint = SintX_mm512_maptype0(s.deltaX,s.deltaY,s.penY,s.Gauss,s.Sm,s.Ivar,s.Bias,s.Pen,s.PenBias,1,deltaX);

		_mm512_maskstoreu_ps(origscore, mask, v_score);
		_mm512_maskstoreu_ps(newscore, mask, new_v_score);
		_mm512_maskstoreu_ps(scoreGT, mask, scratch.scoreGT);
		_mm512_maskstoreu_ps(deltaY, mask, scratch.deltaY);
		_mm512_maskstoreu_ps(penY, mask, scratch.penY);
		_mm512_maskstoreu_ps(Gauss, mask, scratch.Gauss);
		_mm512_maskstoreu_ps(Sm, mask, scratch.Sm);
		_mm512_maskstoreu_ps(Ivar, mask, scratch.Ivar);
		_mm512_maskstoreu_ps(Bias, mask, scratch.Bias);
		_mm512_maskstoreu_ps(Pen, mask, scratch.Pen);
		_mm512_maskstoreu_ps(PenBias, mask, scratch.PenBias);
		_mm512_maskstoreu_ps(Sint, mask, v_Sint);
		_mm512_maskstoreu_ps(OutPen, mask, v512_OutlierPenalty);

		#pragma omp critical
		{
		  printf("refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d:G=%d,T=%d..%d,H=%d..%d:AscoreIK[J]=%0.6f\n",refid,mapid,orientation,I,K,J,G,0,Kmax[G],Hmin,Hmax,AscoreIK[J]);
		  for(int j = 0; j < 16; j++)
		    if(deltaX[j] > 0.0)
		      printf("ht=%d:deltaX=%0.4f,deltaY=%0.4f,penY=%0.6f,Sm=%0.6f,Ivar=%0.8f,Bias=%0.6f,Pen=%0.6f,GaussY=%0.6f,PenBias=%0.6f,OutPen=%0.6f,Sint=%0.6f:AscoreGT[H]=%0.6f,score=%0.6f -> %0.6f\n",
			     j,deltaX[j],deltaY[j],penY[j],Sm[j],Ivar[j],Bias[j],Pen[j],Gauss[j],PenBias[j],OutPen[j],Sint[j],scoreGT[j],origscore[j],newscore[j]);
		  fflush(stdout);
		}
	      }

	      v_score = SintX_scratch_collapse_maptype0(scratch, v_score);
	      score = _mm512_reduce_gmax_ps(v_score);
#endif

	      if(PVERB>=3 && I==I_TRACE && K==K_TRACE && J==J_TRACE){
		#pragma omp critical
		{
		  printf("refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d:G=%d,T=%d..%d,H=%d..%d:A->score(I,K,J)=%0.6f -> %0.6f\n",refid,mapid,orientation,I,K,J,G,0,Kmax[G],Hmin,Hmax,AscoreIK[J],score);
		  fflush(stdout);
		}
	      }
	      if(DEBUG>=2) assert(isfinite(score));
	      if(PVERB>=2 && I==I_TRACE && K==K_TRACE && J==J_TRACE){
		printf("11: A->score(I=%d,K=%d,J=%d)=%0.6f -> %0.6f\n",I,K,J,A->score(I,K,J),score);
		fflush(stdout);
	      }
	      AscoreIK[J] = score;
	    } /* first few iterations of J (maptype==0, most iterations of I) */
	    if(DEBUG>=2) assert(J == jmax+1);
	  
#if (USE_SSE==1 && USE_AVX && USE_RFLOAT==1)
	    for(int jmax3 = min(JMAX[I],JMAX[G]+1) - 7; J <= jmax3; J += 8){/* critical J loop (with maptype==0, most iterations of I) with USE_SSE && USE_AVX && USE_RFLOAT */
	      __m256 v_score = _mm256_loadu_ps(&(AscoreIK[J]));
	      __m256 v_deltaX, v_prev_score, v_XPenJM_Bias, v_XPenJM_Pen, v_XPenJM_PenBias, v_deltaY, v_Ivar, v_penY, v_GaussY, v_Sm;
	      for(int T = Kmax[G]; T >= 0; T--){
		//		CYPen *p = &YPenR[K][T][I][I-K-G];
		CYPen *p= &p_base[-T*N*DELTA_Y-G];
		if(DEBUG>=2 && !(p == &YPenR[K][T][I][I-K-G])){
		  printf("*** INTERNAL ERROR: Pointers don't` match p=%p p_expected=%p p_base=%p K=%d T=%d I=%d G=%d N=%d DELTA_Y=%d\n", p, &YPenR[K][T][I][I-K-G], p_base, K, T, I, G, N, DELTA_Y);
		  fflush(stdout);
		  assert(p == &YPenR[K][T][I][I-K-G]);
		}

		v_deltaY = _mm256_set1_ps(p->deltaY);
		v_Ivar = _mm256_broadcast_ss(&p->Ivar);
		v_penY = _mm256_broadcast_ss(&p->Pen);
		v_GaussY = _mm256_broadcast_ss(&p->Gauss);
		v_Sm = _mm256_broadcast_ss(&p->Sm);
		RFLOAT *AscoreGT = &A->score(G,T,J);
		
		for(int m = 1; m <= 4;m++){
		  if(DEBUG>=2) assert(J-m >= JMIN[G]);
		  if(DEBUG>=2) assert(J-m+7 <= JMAX[G]);
		  if(DEBUG>=2) 
		    for(int j = 0; j < 8; j++)
		      assert(isfinite(AscoreGT[j-m]));
		  v_prev_score = _mm256_loadu_ps(&AscoreGT[-m]);
		  v_deltaX = _mm256_loadu_ps(&XPen->deltaX[m][J]);
		  v_XPenJM_Bias = _mm256_loadu_ps(&XPen->Bias[m][J]);
		  v_XPenJM_Pen = _mm256_loadu_ps(&XPen->Pen[m][J]);
		  v_XPenJM_PenBias = _mm256_loadu_ps(&XPen->PenBias[m][J]);
		  v_score = _mm256_max_ps(v_score, _mm256_add_ps(v_prev_score, SintX_mm256_maptype0(v_deltaX,v_deltaY, v_penY, v_GaussY, v_Sm, v_Ivar, v_XPenJM_Bias, v_XPenJM_Pen, v_XPenJM_PenBias)));
		  if(DEBUG>=2){
		    RFLOAT score[8];
		    _mm256_storeu_ps(score,v_score);
		    for(int j = 0; j < 8; j++)
		      if(!isfinite(score[j])){
			#pragma omp critical
			{
			  printf("refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d+j,G=%d,T=%d,H=%d:JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d:j=%d,prevscore[j]=%0.8f,newscore[j]=%0.8f\n",
				 refid,mapid,orientation,I,K,J,G,T,J-m,JMIN[I],JMAX[I],JMIN[G],JMAX[G],j,AscoreGT[j-m],score[j]);
			  fflush(stdout);
			  assert(isfinite(score[j]));
			}
		      }
		  }
		}
	      }
	      _mm256_storeu_ps(&(AscoreIK[J]), v_score);
	      if(DEBUG>=2) 
		for(int j = J; j < J+8; j++){
		  if(!isfinite(A->score(I,K,j))){
                    #pragma omp critical
		    {
		      printf("refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d(JMIN[I]=%d,JMAX[I]=%d):j=%d,AscoreIK[j]=%0.8f\n",
			     refid,mapid,orientation,I,K,J,JMIN[I],JMAX[I],j,AscoreIK[j]);
		      fflush(stdout);
		      assert(isfinite(A->score(I,K,j)));
		    }
		  }
		}
	    } /* critical J loop (with maptype==0, most iterations of I) with USE_SSE && USE_AVX && USE_RFLOAT */
#endif // (USE_SSE==1 && USE_AVX==1 && USE_RFLOAT==1)

#if (USE_MIC==1 && USE_RFLOAT==1) // D1  
	    int jmax3 = min(JMAX[I],JMAX[G]+1);
	    for(; J <= jmax3; J += 16){/* critical J loop (with maptype==0, most iterations of I) with USE_MIC && USE_RFLOAT */
	      if(DEBUG>=2 && !DIAGONAL) assert(jmax3==M);
	      __mmask16 mask = (J <= jmax3-15) ? 0xffff : ((1u<<(jmax3-J+1))-1);
	      __m512 v_score = _mm512_maskloadu_ps(&AscoreIK[J], mask, rs_heap->thread_block_start(tid), rs_heap->thread_block_siz(), tid, 25);
	      if(PVERB>=3 && I == I_TRACE && K== K_TRACE && J >= J_TRACE && J_TRACE < J+16){ \
		float origscore[16];
		_mm512_maskstoreu_ps(origscore, mask, v_score);

		#pragma omp critical
		{
		  printf("refid=%d,mapid=%d,or=%d:I=%d,K=%d,G=%d,J=%d,jmax3=%d,M=%d,N=%d:mask=0x%04x\n",
			 refid,mapid,orientation,I,K,G,J,jmax3,M,N,mask);
		  for(int j = 0; j < 16; j++)
		    printf("v_score[%d]=%0.6e\n",j,origscore[j]);
		  fflush(stdout);
		}
	      }

	      //for(int T = Kmax[G]; T >= 0; T--){
	      int T=0;
#pragma unroll(0)
              do {
		//		CYPen *p = &YPenR[K][T][I][I-K-G];
		CYPen *p= &p_base[-T*N*DELTA_Y-G];
		if(DEBUG>=2 && ! (p == &YPenR[K][T][I][I-K-G])){
		  printf("*** INTERNAL ERROR: Pointers don't match p=%p p_expected=%p p_base=%p K=%d T=%d I=%d G=%d N=%d DELTA_Y=%d\n", p, &YPenR[K][T][I][I-K-G], p_base, K, T, I, G, N, DELTA_Y);
		  fflush(stdout);
		  assert(p == &YPenR[K][T][I][I-K-G]);
		}
		__m512 v_deltaY = _mm512_set1_ps(p->deltaY);
		__m512 v_Ivar = _mm512_set1_ps(p->Ivar);
		__m512 v_penY = _mm512_set1_ps(p->Pen);
		__m512 v_GaussY = _mm512_set1_ps(p->Gauss);
		__m512 v_Sm = _mm512_set1_ps(p->Sm);

		float *AscoreGT = &A->score(G,T,J);

		// new unrolled vector code

                #define BLOCK(H) {		  /* NOTE : J >= JMIN[I] + 4 && JMIN[I] >= JMIN[G] */ \
		  if(DEBUG>=2) assert(J-H >= JMIN[G]);		                                      \
		  if(DEBUG>=2) assert(min(jmax3,J+15)-H <= JMAX[G]);	                              \
		  __m512 v_prev_score = _mm512_maskloadu_ps(&AscoreGT[-H], mask, rs_heap->thread_block_start(tid), rs_heap->thread_block_siz(), tid, 31, H); \
		  __m512 v_deltaX = _mm512_maskloadu_ps(&deltaX##H[J],mask, XPen->mem_pool, XPen->mem_pool_size, tid, 32, H); \
		  __m512 v_XPenJM_Bias = _mm512_maskloadu_ps(&Bias##H[J],mask, XPen->mem_pool, XPen->mem_pool_size, tid, 33, H); \
		  __m512 v_XPenJM_Pen = _mm512_maskloadu_ps(&Pen##H[J],mask, XPen->mem_pool, XPen->mem_pool_size, tid, 34, H); \
		  __m512 v_XPenJM_PenBias = _mm512_maskloadu_ps(&PenBias##H[J],mask, XPen->mem_pool, XPen->mem_pool_size, tid, 35, H, J, M); \
		  __m512 v_SintX = SintX_mm512_maptype0(v_deltaX,v_deltaY, v_penY, v_GaussY, v_Sm, v_Ivar, v_XPenJM_Bias, v_XPenJM_Pen, v_XPenJM_PenBias,0,NULL); \
		  __m512 v_newscore = _mm512_add_ps(v_prev_score, v_SintX);\
		  __m512 v_score_orig; if(DEBUG>=2 || TRACE) v_score_orig = v_score;		\
		  v_score = _mm512_gmax_ps(v_score, v_newscore);                                        \
		  if(DEBUG>=2 || (PVERB && I==I_TRACE && J <= J_TRACE && J_TRACE <= min(jmax3,J+15))){ \
		    RFLOAT nscore[16], origscore[16], newscore[16], deltaX[16], deltaY[16], prevscore[16], SintX[16], penY[16],GaussY[16],Sm[16],Ivar[16],XPen[16],XBias[16],XPenBias[16]; \
		    for(int j = 0; j < 16; j++) \
		      deltaX[j] = -1.0;         \
		    _mm512_maskstoreu_ps(deltaX, mask, v_deltaX); \
		    _mm512_maskstoreu_ps(deltaY, mask, v_deltaY); \
		    _mm512_maskstoreu_ps(origscore, mask, v_score_orig);		        \
		    _mm512_maskstoreu_ps(newscore, mask, v_newscore);			\
		    _mm512_maskstoreu_ps(nscore, mask, v_score);			        \
		    _mm512_maskstoreu_ps(prevscore, mask, v_prev_score);			        \
		    _mm512_maskstoreu_ps(SintX, mask, v_SintX); \
		    _mm512_maskstoreu_ps(penY, mask, v_penY); \
		    _mm512_maskstoreu_ps(GaussY, mask, v_GaussY); \
		    _mm512_maskstoreu_ps(Sm, mask, v_Sm); \
		    _mm512_maskstoreu_ps(Ivar, mask, v_Ivar); \
		    _mm512_maskstoreu_ps(XPen, mask, v_XPenJM_Pen); \
		    _mm512_maskstoreu_ps(XBias, mask, v_XPenJM_Bias); \
		    _mm512_maskstoreu_ps(XPenBias, mask, v_XPenJM_PenBias); \
		    for(int j = J; j <= min(jmax3,J+15); j++){				\
		      if (deltaX[j] > 0.0 && ((DEBUG && !(isfinite(nscore[j-J]) && isfinite(origscore[j-J] && isfinite(AscoreGT[j-J-H])))) || (PVERB && I==I_TRACE && K==K_TRACE && j==J_TRACE))){ \
		        printf("S:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d,jmax3=%d:AscoreGT[H]=%0.8f,SintX=%0.8f,newscore=%0.8f,score=%0.8e -> %0.8f\n", \
			       refid,mapid,orientation,I,K,j,G,T,j-H,jmax3,AscoreGT[j-J-H],SintX[j-J],newscore[j-J],origscore[j-J],nscore[j-J]); \
			printf("\t deltaX=%0.4f,deltaY=%0.4f,penY=%0.6f,Sm=%0.6f,Ivar=%0.8f,GaussY=%0.6f,XBias=%0.6f,XPen=%0.6f,XPenBias=%0.6f\n",\
			       deltaX[j-J],deltaY[j-J],penY[j-J],Sm[j-J],Ivar[j-J],GaussY[j-J],XBias[j-J],XPen[j-J],XPenBias[j-J]); \
		        fflush(stdout);							\
			assert(isfinite(nscore[j-J])); \
			assert(isfinite(origscore[j-J]));		\
			assert(isfinite(AscoreGT[j-J-H]));                \
			assert(fabs(prevscore[j-J] - AscoreGT[j-J-H]) < 1e-10); \
		      }                                            			\
	            }				                                        \
		  } 									\
                }

	        BLOCK(1);
		BLOCK(2);
		BLOCK(3);
		BLOCK(4);
#undef BLOCK

		T++;
	      } while(T <= Kmax[G]);

	      if(DEBUG>=2) assert(IMIN <= I && I <= IMAX && K <= A->kmax);
	      _mm512_maskstoreu_ps(&AscoreIK[J], mask, v_score);
	      if(DEBUG>=2) 
		for(int j = J; j <= min(jmax3,J+15); j++)
		  assert(isfinite(A->score(I,K,j)));
	    } /* critical J loop (with maptype==0, most iterations of I) with USE_MIC && USE_RFLOAT */
	    if(DIAGONAL) 
	      J = max(jmax,jmax3)+1;
	    else if(DEBUG>=2) assert(max(jmax,jmax3) == JMAX[I]);

#endif // USE_MIC && USE_RFLOAT==1

	    for(int jmax3 = JMAX[I]; J <= jmax3; J++){/* remaining iterations of J (with maptype==1, most iterations of I) all use cases */
	      RFLOAT score = AscoreIK[J];
	      int m_min = max(1,J - JMAX[G]);
	      for(int T = Kmax[G]; T >= 0; T--){
		//		CYPen *p = &YPenR[K][T][I][I-K-G];
		CYPen *p= &p_base[-T*N*DELTA_Y-G];
		if(DEBUG>=2 && ! (p == &YPenR[K][T][I][I-K-G])){
		  printf("*** INTERNAL ERROR: Pointers don't match p=%p p_expected=%p p_base=%p K=%d T=%d I=%d G=%d N=%d DELTA_Y=%d\n", p, &YPenR[K][T][I][I-K-G], p_base, K, T, I, G, N, DELTA_Y);
		  fflush(stdout);
		  assert(p == &YPenR[K][T][I][I-K-G]);
		}
		RFLOAT deltaY = p->deltaY;
		RFLOAT Ivar = p->Ivar;
		RFLOAT penY = p->Pen;
		RFLOAT GaussY = p->Gauss;
		RFLOAT Sm = p->Sm;
		
		RFLOAT *AscoreGT = &A->score(G,T,0);
#define BLOCK(m) {							\
		  if(!DIAGONAL || m >= m_min){\
		    if(DEBUG>=2) assert(J-m <= JMAX[G] && J-m >= JMIN[G]);	\
		    RFLOAT newscore = AscoreGT[J-m] + SintX_maptype0(deltaX##m[J],deltaY,penY,GaussY,Sm,Ivar,Pen##m[J],Bias##m[J], PenBias##m[J],0); \
		    if ((DEBUG>=2 && !isfinite(score)) || (PVERB && I==I_TRACE && K==K_TRACE && J==J_TRACE)){ \
		      printf("D:refid=%d,mapid=%d,or=%d:I=%d,K=%d,J=%d,G=%d,T=%d,H=%d,IMIN=%d,IMAX=%d,JMIN[I]=%d,JMAX[I]=%d,JMIN[G]=%d,JMAX[G]=%d:AscoreGT[H]=%0.8f,newscore=%0.8f,score=%0.8e -> %0.8f\n", \
			     refid,mapid,orientation,I,K,J,G,T,J-m,IMIN,IMAX,JMIN[I],JMAX[I],JMIN[G],JMAX[G],AscoreGT[J-m],newscore,score, max(score,newscore)); \
		      fflush(stdout);					\
		      assert(isfinite(score));	\
		    }							\
		    score = max(score,newscore);			\
		  }							\
		}
		BLOCK(4);
		BLOCK(3);
		BLOCK(2);
		BLOCK(1);
#undef BLOCK
	      }/* T loop */
	      if(DEBUG>=2) assert(isfinite(score));
	      if(PVERB>=2 && I==I_TRACE && K==K_TRACE && J==J_TRACE){
		printf("12: A->score(I=%d,K=%d,J=%d)=%0.6f -> %0.6f\n",I,K,J,A->score(I,K,J),score);
		fflush(stdout);
	      }
	      AscoreIK[J] = score;
	    } /* remaining iterations of J (with maptype==0, most iterations of I) all use cases */

	  } /* G loop */
	} /* K loop */
      } /* I loop */
    } /* maptype == 0 */
    
    if(PVERB>=2 && IMIN <= I_TRACE && I_TRACE <= IMAX && JMIN[I_TRACE] <= J_TRACE && J_TRACE <= JMAX[I_TRACE]){
      printf("20: refid=%d(%lld),mapid=%d(%lld),or=%d:A->score(I=%d,K=%d,J=%d)=%0.6f\n",refid,rmap->id,mapid,nanomap->id,orientation,I_TRACE,K_TRACE,J_TRACE,A->score(I_TRACE,K_TRACE,J_TRACE));
      fflush(stdout);
    }
  } // if(1)
}
#endif // USE_MIC

#endif // REFALIGN_3Drecurrance
