#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "parameters.h"

#include "Calign.h"
#include "Calign_brief.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Calign.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

// see Calign_brief.h
void Calign_brief::update(const Calign &align, int M)
{
  score = align.score;
  numpairs = align.numpairs;
  orientation = align.orientation;
  mapid1 = align.mapid1;
  mapid2 = align.mapid2;
  scaleID = align.scaleID;
  logPV = align.logPV;
  if(GOODMAPS && numpairs <= 0){
    sites1_0 = 0;
    sites1_M = M;
  } else {
    sites1_0 = align.sites1[0];
    sites1_M = align.sites1[align.numpairs-1];
  }
}

/** increase size of array alignment[] to at least [0..num-1] and allocate each element */
void maxalignalloc(size_t num, Calign **&alignment, size_t &numaligns, size_t &maxaligns)
{
  size_t newaligns = num;
  if(newaligns <= maxaligns)
    return;

  if(newaligns < 1024)
    newaligns = 1024;
  if(newaligns*4 < maxaligns*5)
    newaligns = (maxaligns*5)/4;

  if(VERB/* HERE HERE HERE >=2 */){
    printf("maxalignalloc(%llu):maxaligns=%llu->%llu,sizeof(Calign)=%lu bytes\n",
	   (unsigned long long)num, (unsigned long long)maxaligns,(unsigned long long)newaligns,sizeof(Calign));
    fflush(stdout);
  }

  if(maxaligns <= 0){
    Calign **newalignment = new Calign*[newaligns];
    for(register size_t i= 0; i < newaligns; i++)
      newalignment[i] = new Calign[colors==1 ? 1 : colors+1];
    alignment = newalignment;
  } else {
    Calign **origalignment = alignment;
    Calign **newalignment = new Calign*[newaligns];
    for(register size_t i=0;i<maxaligns;i++)
      newalignment[i] = origalignment[i];
    if(colors==1) {
      for(register size_t i=maxaligns;i<newaligns;i++)
        newalignment[i] = new Calign[1];
    } else {
      for(register size_t i=maxaligns;i<newaligns;i++)
        newalignment[i] = new Calign[colors+1];
    }

    alignment = newalignment;
    delete [] origalignment;
  }
  maxaligns = newaligns;
}

/** same as maxalignalloc, but each new element is initialized to NULL */
void maxalignallocNULL(size_t num, Calign **&alignment, size_t &numaligns, size_t &maxaligns, int *alignment_blockid)
{
  size_t newaligns = num;
  if(newaligns <= maxaligns)
    return;

  if(newaligns < 1024)
    newaligns = 1024;
  if(newaligns*4 < maxaligns*5)
    newaligns = (maxaligns*5)/4;

  if(VERB/* HERE HERE HERE >=2 */){
    if(alignment_blockid != NULL)
      printf("maxalignallocNULL(%lu):maxaligns=%lu->%lu, sizeof(Calign)= %lu bytes, alignment_blockid=%d\n",num,maxaligns,newaligns,sizeof(Calign),*alignment_blockid);
    else
      printf("maxalignallocNULL(%lu):maxaligns=%lu->%lu, sizeof(Calign)= %lu bytes\n",num,maxaligns,newaligns,sizeof(Calign));
    fflush(stdout);
  }

  if(maxaligns <= 0){
    Calign **newalignment = new Calign*[newaligns];
    for(register size_t i=0;i < newaligns; i++)
      newalignment[i] = NULL;
    alignment = newalignment;
    if(alignment_blockid != NULL)
      *alignment_blockid = -1;
  } else {
    Calign **origalignment = alignment;
    Calign **newalignment = new Calign*[newaligns];
    for(register size_t i=0;i<maxaligns;i++)
      newalignment[i] = origalignment[i];
    for(register size_t i=maxaligns;i<newaligns;i++)
      newalignment[i] = NULL;
    alignment = newalignment;
    if(alignment_blockid != NULL){
      if(*alignment_blockid < 0)
	delete [] origalignment;
      *alignment_blockid = -1;
    } else 
      delete [] origalignment;
  }
  maxaligns = newaligns;
}

extern double Lgamma(double x);// see pvalue.cpp 

/** compute Probability of Chi-Square distribution of r degrees of freedom exceeding value ChiSq */

double ChiSqPvalue(int r, double ChiSq)
{
  /* Based on computing gamma_p(r/2,ChiSq/2), which has a finite expansion with floor(r/2) terms */

  if(DEBUG && !(r > 0)){
    printf("ChiSqPvalue(r=%d,ChiSq=%0.4f)\n",r,ChiSq);
    fflush(stdout);
    assert(r > 0);
  }
  if(DEBUG) assert(ChiSq >= 0.0);

  register double x = ChiSq*0.5;
  register double n = r/2;

  register double sum;

  if((r%2)){ /* odd value of r */
    if(DEBUG>=2 && !(fabs(2.0*n+1.0 - r) <= 1e-8)){
      printf("r=%d,ChiSq=%0.4f: n=%0.4f,x=%0.4f\n", r, ChiSq, n, x);
      fflush(stdout);
      assert(fabs(2.0*n + 0.5 - r) <= 0.001);
    }
    if(r == 1)
      return erfc(sqrt(x));
    if(x < 1e-200){
      double ret = 2.0 * sqrt(x/M_PI) * exp(-x) + erfc(sqrt(x));
      if(DEBUG/* HERE >=2 */ && !(isfinite(ret) && ret >= 0.0)){
        #pragma omp critical
	{
	  printf("ChiSqPvalue(%d,%0.8e)=%0.8e: n=%0.1f,x=%0.8e\n",
		 r,ChiSq,ret,n,x);
	  fflush(stdout);
	  assert(isfinite(ret) && ret >= 0.0);
	}
      }
      return ret;
    }

    double ret = 0.0;
    double logx = log(x);
    double m = (r - 1.0)*0.5;
    if(m * logx < 700.0){
      register double term = 2.0 * x;
      sum = term;/* num/den */
      for(register double k = 1.5; k < n; k++){
	sum += term *= x/k;
	if(DEBUG>=2 && !(isfinite(sum) && sum >= 0.0)){
	  printf("ChiSqPvalue(r=%d,ChiSq=%0.8e):k=%0.1f,n=%0.1f,term=%0.8e,sum=%0.8e\n",
		 r,ChiSq,k,n,term,sum);
	  fflush(stdout);
	  assert(isfinite(sum) && sum >= 0.0); 
	}
      }
      ret = exp(-x + log(sum))/sqrt(M_PI * x) + erfc(sqrt(x));

      if(DEBUG/* HERE >=2 */ && !(isfinite(ret) && ret >= 0.0)){
        #pragma omp critical
	{
	  printf("ChiSqPvalue(%d,%0.8e)=%0.8e: n=%0.1f,x=%0.8e,sum=%0.8e\n",
		 r,ChiSq,ret,n,x,sum);
	  fflush(stdout);
	  sum = term = 2.0*x;
	  for(double k = 1.5; k < n; k++){
	    sum += term *= x/k;
	    printf("  k=%0.1f:term=%0.8e,sum=%0.8e\n",k,term,sum);
	  }
	  fflush(stdout);
	  assert(isfinite(ret) && ret >= 0.0);
	}
      }
    } else {
      register double term = 1.0;
      sum = term;
      for(register double k = 1; k < m; k++)
	sum += term *= (n-k)/x;
      double log2 = log(2.0);
      ret = exp(-x + log(sum) + 2.0*n*log2 + m*logx - Lgamma(2.0*m+1.0) - Lgamma(m+1.0))/sqrt(M_PI * x) + erfc(sqrt(x));
      if(DEBUG/* HERE >=2 */ && !(isfinite(ret) && ret >= 0.0)){
        #pragma omp critical
	{
	  printf("ChiSqPvalue(%d,%0.8e)=%0.8e: n=%0.1f,x=%0.8e,logx=%0.8f,m=%0.1f,sum=%0.8e\n",
		 r,ChiSq,ret,n,x,logx,m,sum);
	  sum = term = 1.0;
	  for(register double k = 1; k < m; k++){
	    sum += term *= (n-k)/x;
	    printf("  k=%0.1f: term=%0.8e,sum=%0.8e\n",
		   k,term,sum);
	  }
	  fflush(stdout);
	  assert(isfinite(ret) && ret >= 0.0);
	}
      }
    }
    return ret;
  } else {/* even value of r */
    if(DEBUG>=2) assert(fabs(2.0*n - r) <= 1e-8);
    if(r == 2)
      return exp(-x);
    if(x < 1e-300){
      double ret = 1.0;
      if(DEBUG/* HERE >=2 */ && !(isfinite(ret) && ret >= 0.0)){
        #pragma omp critical
	{
	  printf("ChiSqPvalue(%d,%0.8e)=%0.8e: n=%0.1f,x=%0.8e\n",
		 r,ChiSq,ret,n,x);
	  fflush(stdout);
	  assert(isfinite(ret) && ret >= 0.0);
	}
      }
      return ret;
    }

    double ret = 0.0;
    double logx = log(x);
    register double term = 1.0;
    sum = term;
    if((n-1.0)*logx < 700.0){
      for(register double k = 1.0; k < n; k++)
	sum += term *= x/k;
      ret = exp(-x + log(sum));

      if(DEBUG/* HERE >=2 */ && !(isfinite(ret) && ret >= 0.0)){
        #pragma omp critical
	{
	  printf("ChiSqPvalue(%d,%0.8e)=%0.8e: n=%0.1f,x=%0.8e,sum=%0.8e\n", r,ChiSq,ret,n,x,sum);
	  fflush(stdout);
	  assert(isfinite(ret) && ret >= 0.0);
        }
      }
    } else {/* last term x^(n-1) overflows : use alternate summation in reverse order */
      for(register double k = 1.0; k < n; k++)
	sum += term *= (n-k)/x;
      ret = exp(-x + (n-1.0)*logx - Lgamma(n) + log(sum));

      if(DEBUG/* HERE >=2 */ && !(isfinite(ret) && ret >= 0.0)){
        #pragma omp critical
	{
	  printf("ChiSqPvalue(%d,%0.8e)=%0.8e: n=%0.1f,x=%0.8e,logx=%0.8f,Lgamma(n)=%0.8f,sum=%0.8e\n", r,ChiSq,ret,n,x,logx,Lgamma(n),sum);
	  fflush(stdout);
	  assert(isfinite(ret) && ret >= 0.0);
        }
      }
    }

    return ret;
  }
}

extern Cmap **YYmap, **XXmap;

/** check if entire repeat region can be shifted by 1 or 2 (up to RepeatMaxShift) sites without increasing the mean-square sizing error significantly */
/** return shift value that revealed repeat */

/** depends on command line parameters RepeatMaxShift and RepeatPvalueRatio */
int Calign::IsRepeatRegion()
{
  if(((spots_filename[0] || (svcheck && numrefmaps > 0)) && resSD[0] > 0.0) || PairSplit || PairMerge || RMS_MIN < 0.0 /* HERE HERE*/||colors==2){ /* NOTE : If RMS_MIN >= 0, pairwise alignments are checked for duplicate maps */
    if(RepeatMaxShift <= 0)
      return 0;
    if(repeat >= 0)/* repeat already checked using -RepeatRec (superior method) */
      return repeat;
    if(DEBUG) assert(!RepeatRec);
  }

  if(VERB>=2 && ((spots_filename[0] || (svcheck && numrefmaps > 0)) && resSD[0] > 0.0)){
    printf("IsRepeatRegion() start: mapid1=%d,mapid2=%d,or=%d:score=%0.6f,logPV=%0.2f,numpairs=%d,repeat=%d\n",mapid1,mapid2,orientation,score,logPV,numpairs,repeat);
    fflush(stdout);
  }
  if(colors > 1){
    fprintf(stderr,"IsRepeatRegion(): not implemented for -colors %d\n",colors);
    exit(1);
  }

  int K = numpairs-1;
  if(K <= 1)
    return 0;

  int result = 0;

  Cmap *Ymap = spots_filename[0] ? refmap[mapid1] : YYmap[mapid1];
  Cmap *Xmap = spots_filename[0] ? Gmap[mapid2] : XXmap[mapid2];
  FLOAT *Y = Ymap->site[0];
  FLOAT *X = Xmap->site[0];
  int N = Ymap->numsite[0];
  int M = Xmap->numsite[0];
  
  int maxshift = min(RepeatMaxShift,K-1);

  /* compute the interval sizes of the current alignment */
  double *Yfrag = new double[K*2+maxshift+1];
  double *Xfrag = &Yfrag[K];
  double *ChiSquare = &Yfrag[2*K];
  int *outlier = new int[K + maxshift + 1];
  int *ChiSqDegree = &outlier[K];

  if(spots_filename[0] && resSD[0] > 0.0){/* Guided (-ref) Alignment with resolution modeling */
    double Xscale = scaleID ? ScaleFactor[scaleID] : 1.0;
    double y,py = Yc(Y,sites1[0],sitesK1[0]);
    double x,px = orientation ? X[M+1] - X[M+1-sites2[0]] : X[sites2[0]];
    for(int i = 1; i <= K; i++, py = y, px = x){
      y = Yc(Y,sites1[i],sitesK1[i]);
      x = orientation ? X[M+1] - X[M+1-sites2[i]] : X[sites2[i]];
      Yfrag[i-1] = y - py;
      Xfrag[i-1] = (x - px)*Xscale;
      outlier[i-1] = (outscore[i] + (FLOAT)0.01 /* WAS14 1e-6 */ < iscore[i] || (outlierExtend && (sites2[i-1] >= sites2[i] || sites1[i-1] >= sites1[i]-sitesK1[i]))) ? 1 : 0;
    }
  } else {/* pairwise Alignment OR Guided Alignment without resolution modeling */
    double y,py = Y[sites1[0]];
    double x,px = orientation ? X[M+1] - X[M+1-sites2[0]] : X[sites2[0]];
    for(int i = 1; i <= K; i++, py = y, px = x){
      y = Y[sites1[i]];
      x = orientation ? X[M+1] - X[M+1-sites2[i]] : X[sites2[i]];
      Yfrag[i-1] = y - py;
      Xfrag[i-1] = x - px;
      outlier[i-1] = (outscore[i] + (PFLOAT)0.01 /* WAS14 1e-6 */ < iscore[i] || (outlierExtend && (sites2[i-1] >= sites2[i] || sites1[i-1] >= sites1[i]))) ? 1 : 0;
    }
  }

  /* compute the Chi-square Pvalue (for MSE greater than the observed value) for the current alignment and for each shifted alignment */
  double varF = spots_filename[0] ? SF[0]*SF[0] : 2.0*SF[0]*SF[0];
  double var = fabs(SD[0])*SD[0];
  double varR = SR[0]*SR[0];

  if(spots_filename[0]){

    /* compute unshifted ChiSquare */
    double sum = 0.0;
    int cnt = 0;
#pragma vector unaligned
    for(int i = 0; i < K; i++){
      if(outlier[i])
	continue;
      double y = Yfrag[i];
      double x = Xfrag[i];
      double vary = varF + var*y;
      if(QUADRATIC_VARIANCE)
	vary += varR * (x*x + y*y);

      double err = x-y;
      if(DEBUG) assert(vary > 0.0);
      sum += err*err/vary;
      cnt++;
    }
    if(DEBUG && !(cnt > 0)){
      printf("mapid1=%d(id=%lld),mapid2=%d(id=%lld),or=%d:numpairs=%d,score=%0.2f,logPV=%0.2f,resSD=%0.3f\n",
	     mapid1,Ymap->id,mapid2,Xmap->id,orientation,numpairs,score,logPV,resSD[0]);
      for(int i = 0; i < K; i++)
	printf("i=%d:outlier[i]=%d(outscore[i+1]=%0.6f,iscore[i+1]=%0.6f),Yfrag[i]=%0.3f,Xfrag[i]=%0.3f\n",i,outlier[i],outscore[i],iscore[i],Yfrag[i],Xfrag[i]);
      fflush(stdout);
      assert(cnt > 0);
    }
    ChiSquare[0] = sum;
    ChiSqDegree[0] = cnt;
    
    double Pvalue1 = ChiSqPvalue(cnt,sum);

    for(int shift = 1; shift <= maxshift; shift++){
      /* first try incrementing Y index (shifting Y left) */
      sum = 0.0;
      cnt = 0;
#pragma vector unaligned
      for(int i = shift; i < K; i++){
	if(outlier[i] || outlier[i-shift])
	  continue;
	double y = Yfrag[i];
	double x = Xfrag[i-shift];
	double vary = varF + var*y;
	if(QUADRATIC_VARIANCE)
	  vary += varR * (x*x + y*y);
	double err = x-y;
	if(DEBUG) assert(vary > 0.0);
	sum += err*err/vary;
	cnt++;
      }
      if(!cnt)/* cannot check for repeat without at least one non-outlier aligned interval */
	break;

      ChiSquare[shift] = sum;
      ChiSqDegree[shift] = cnt;

      /* next try decrementing Y index (shifting Y right) */
      sum = 0.0;
      cnt = 0;
#pragma vector unaligned
      for(register int i = shift; i < K; i++){
	if(outlier[i] || outlier[i-shift])
	  continue;
	register double y = Yfrag[i-shift];
	register double x = Xfrag[i];
	register double vary = varF + var*y;
	register double err = x-y;
	if(QUADRATIC_VARIANCE)
	  vary += varR * (x*x + y*y);
	if(DEBUG) assert(vary > 0.0);
	sum += err*err/vary;
	cnt++;
      }

      /* to be classified as a repeat both forward and reverse shifted alignments should confirm repeat pattern by having small mean-square error (relative to original alignment) */
      if(DEBUG) assert(cnt == ChiSqDegree[shift]);
      if(sum > ChiSquare[shift])
	ChiSquare[shift] = sum;

      register double Pvalue2 = ChiSqPvalue(cnt,ChiSquare[shift]);
      if(Pvalue2 > Pvalue1 * RepeatPvalueRatio){/* probably a repeat */
	if(VERB>=2){
	  printf("refid=%d(%lld),mapid=%d(%lld),or=%d,score=%0.4f,logPV=%0.4f:K=%d:ChiSq[0]=%0.4f(DF=%d),Pvalue1=%0.8e:ChiSq[%d]=%0.4f(DF=%d),Pvalue2=%0.8e: Repeat detected\n",
		 mapid1,Ymap->id,mapid2,Xmap->id,orientation,score,logPV,K,ChiSquare[0],ChiSqDegree[0],Pvalue1,shift,ChiSquare[shift],ChiSqDegree[shift],Pvalue2);
	  for(register int i = 0; i < K; i++)
	    printf("i=%d:Xfrag[i]=%0.3f,Yfrag[i]=%0.3f,outlier[i]=%d,outscore[i+1]=%0.4f,iscore[i+1]=%0.4f\n",
		   i,Xfrag[i],Yfrag[i],outlier[i],outscore[i+1],iscore[i+1]);
	  fflush(stdout);
	}
	if(DEBUG) assert(shift > 0);
	result = shift;
	if(DEBUG) maxshift = shift;
	break;
      }
      if(VERB>=2){
	printf("refid=%d(%lld),mapid=%d(%lld),or=%d,score=%0.4f,logPV=%0.4f:K=%d:ChiSq[0]=%0.4f(DF=%d),Pvalue1=%0.8e:ChiSq[%d]=%0.4f(DF=%d),Pvalue2=%0.8e: No Repeat detected\n",
	       mapid1,Ymap->id,mapid2,Xmap->id,orientation,score,logPV,K,ChiSquare[0],ChiSqDegree[0],Pvalue1,shift,ChiSquare[shift],ChiSqDegree[shift],Pvalue2);
	fflush(stdout);
      }
    }
    if((DEBUG && !(Pvalue1 > 1e-16)) || VERB>=2){/* This can happen for False Positive alignments : very rare if stringent -T -A -S thresholds are used */
      #pragma omp critical
      {
	if(VERB<=1)
	  for(int shift = 1; shift <= maxshift; shift++){
	    assert(ChiSqDegree[shift] > 0);
	    printf("refid=%d(%lld),mapid=%d(%lld),or=%d,score=%0.4f,logPV=%0.4f:K=%d:ChiSq[0]=%0.4f(DF=%d),Pvalue1=%0.8e:ChiSq[%d]=%0.4f(DF=%d),Pvalue2=%0.8e: No Repeat detected\n",
		   mapid1,Ymap->id,mapid2,Xmap->id,orientation,score,logPV,K,ChiSquare[0],ChiSqDegree[0],ChiSqPvalue(ChiSqDegree[0],ChiSquare[0]),shift,ChiSquare[shift],ChiSqDegree[shift],ChiSqPvalue(ChiSqDegree[shift],ChiSquare[shift]));
	  }
	double rsum = 0.0;
	for(int i = 0; i < K; i++){
	  if(outlier[i]){
	    printf("i=%d:Xfrag[i]=%0.3f,Yfrag[i]=%0.3f,outlier[i]=%d,outscore[i+1]=%0.6f,iscore[i+1]=%0.6f\n",
		   i,Xfrag[i],Yfrag[i],outlier[i],outscore[i+1],iscore[i+1]);
	    continue;
	  }
	  double y = Yfrag[i];
	  double x = Xfrag[i];
	  double vary = varF + var*y;
	  if(QUADRATIC_VARIANCE)
	    vary += varR * (x*x + y*y);

	  double err = x-y;
	  rsum += err*err/vary;
	  printf("i=%d:Xfrag[i]=%0.3f,Yfrag[i]=%0.3f,outlier[i]=%d,outscore[i+1]=%0.6f,iscore[i+1]=%0.6f,rms=%0.4f,cumrms=%0.4f\n",
		 i,Xfrag[i],Yfrag[i],outlier[i],outscore[i+1],iscore[i+1],err*err/vary,rsum);
	}
	printf("varF=%0.6f,var=%0.6f\n",varF,var);
	if(!(Pvalue1 > 1e-16))
	  printf("WARNING: Pvalue1=%0.8e\n",Pvalue1);
	fflush(stdout);

	//	assert(Pvalue1 > 1e-16);
	if(DEBUG) assert(isfinite(Pvalue1));
      }
    }
    if(DEBUG>=2) assert(isfinite(Pvalue1));
  } else {/* NOTE : only difference is that vary = varF + var*(x+y) instead of vary = varF + var*y */
    /* compute unshifted ChiSquare */
    double sum = 0.0;
    int cnt = 0;
#pragma vector unaligned
    for(int i = 0; i < K; i++){
      if(outlier[i])
	continue;
      double y = Yfrag[i];
      double x = Xfrag[i];
      double vary = varF + var*(x+y);
      if(QUADRATIC_VARIANCE)
	vary += varR * (x*x + y*y);

      double err = x-y;
      if(DEBUG) assert(vary > 0.0);
      sum += err*err/vary;
      cnt++;
    }
    if(DEBUG && !(cnt > 0)){
      #pragma omp critical
      {
	printf("WARNING:IsRepeatRegion():mapid1=%d(id=%lld),mapid2=%d(id=%lld),or=%d:score=%0.6f,logPV=%0.2f,numpairs=%d,repeat=%d:K=%d,cnt=%d:\n",
	       mapid1,Ymap->id,mapid2,Xmap->id,orientation,score,logPV,numpairs,repeat,K,cnt);
	for(int i = 0; i < K; i++)
	  printf("i=%d:outlier[i]=%d,Xfrag[i]=%0.4f,Yfrag[i]=%0.4f,outscore[i+1]=%0.6f,iscore[i+1]=%0.6f\n",i,outlier[i],Xfrag[i],Yfrag[i],outscore[i],iscore[i]);
	printf("No non-outlier intervals!\n");
	fflush(stdout);
	assert(cnt > 0);
      }
    }

    ChiSquare[0] = sum;
    ChiSqDegree[0] = cnt;

    double Pvalue1 = ChiSqPvalue(cnt,sum);

    if((this->norm = (cnt>0) ? sum/cnt : 1.0) < RMS_MIN && cnt >= RMS_SITES && abs(M-N) <= RMS_MIS && numpairs >= min(N,M) - RMS_MIS && !orientation && Pvalue1 < RMS_PVALUE){
      if(VERB){
	printf("WARNING:mapid1=%d(id=%lld from %s),mapid2=%d(id=%lld from %s),or=%d,N=%d,M=%d:numpairs=%d,score=%0.2f,logPV=%0.2f,RMS norm=%0.6f(DF=%d,Pvalue=%0.6e): Maps are almost identical\n",
	       mapid1,Ymap->id,vfixx_filename[Ymap->fileid],mapid2,Xmap->id,vfixx_filename[Xmap->fileid],orientation,N,M,numpairs,score,logPV,norm,cnt, 1.0-Pvalue1);
	fflush(stdout);
      }
    }

    if(repeat >= 0 || RepeatMaxShift <= 0){/* repeat already checked using -RepeatRec (superior method) OR repeat checking disabled */
      delete [] Yfrag;
      delete [] outlier;
      return max(0,repeat);
    }
    if(DEBUG) assert(!RepeatRec);

    for(int shift = 1; shift <= maxshift; shift++){
      /* first try incrementing Y index (shifting Y left) */
      sum = 0.0;
      cnt = 0;
#pragma vector unaligned
      for(int i = shift; i < K; i++){
	if(outlier[i] || outlier[i-shift])
	  continue;
	double y = Yfrag[i];
	double x = Xfrag[i-shift];
	double vary = varF + var*(x+y);
	if(QUADRATIC_VARIANCE)
	  vary += varR * (x*x + y*y);

	double err = x-y;
	if(DEBUG) assert(vary > 0.0);
	sum += err*err/vary;
	cnt++;
      }
      if(!cnt)/* cannot check for repeat without at least one non-outlier aligned interval */
	break;

      ChiSquare[shift] = sum;
      ChiSqDegree[shift] = cnt;

      /* next try decrementing Y index (shifting Y right) */
      sum = 0.0;
      cnt = 0;
#pragma vector unaligned
      for(int i = shift; i < K; i++){
	if(outlier[i] || outlier[i-shift])
	  continue;
	double y = Yfrag[i-shift];
	double x = Xfrag[i];
	double vary = varF + var*(x+y);
	if(QUADRATIC_VARIANCE)
	  vary += varR * (x*x + y*y);

	double err = x-y;
	if(DEBUG) assert(vary > 0.0);
	sum += err*err/vary;
	cnt++;
      }

      /* to be classified as a repeat both forward and reverse shifted alignments should confirm repeat pattern by having small mean-square error (relative to original alignment) */
      if(DEBUG) assert(cnt > 0);
      if(DEBUG) assert(cnt == ChiSqDegree[shift]);
      if(sum > ChiSquare[shift])
	ChiSquare[shift] = sum;

      register double Pvalue2 = ChiSqPvalue(cnt,ChiSquare[shift]);
      if(Pvalue2 > Pvalue1 * RepeatPvalueRatio){/* probably a repeat */
	if(VERB>=2){
	  #pragma omp critical
	  {
	    printf("Yid=%d(id=%lld),Xid=%d(id=%lld),M=%d,N=%d,or=%d,score=%0.4f,logPV=%0.4f,sites1[0]=%d,sites2[0]=%d:K=%d:ChiSq[0]=%0.4f(DF=%d),Pvalue1=%0.8e:ChiSq[%d]=%0.4f(DF=%d),Pvalue2=%0.8e: Repeat detected\n",
		   mapid1,Ymap->id,mapid2,Xmap->id,M,N,orientation,score,logPV,sites1[0],sites2[0],K,ChiSquare[0],ChiSqDegree[0],Pvalue1,shift,ChiSquare[shift],ChiSqDegree[shift],Pvalue2);
	    for(register int i = 0; i < K; i++)
	      printf("i=%d:Xfrag[i]=%0.3f,Yfrag[i]=%0.3f,sites1[i+1]=%d,sites2[i+1]=%d,outlier[i]=%d,outscore[i+1]=%0.4f,iscore[i+1]=%0.4f\n",
		     i,Xfrag[i],Yfrag[i],sites1[i+1],sites2[i+1],outlier[i],outscore[i+1],iscore[i+1]);
	    fflush(stdout);
	  }
	}
	if(DEBUG) assert(shift > 0);
	result = shift;
	if(DEBUG) maxshift = shift;
	break;
      }
      if(VERB>=2){
	printf("id=%lld,%lld,or=%d,score=%0.4f,logPV=%0.4f:K=%d:ChiSq[0]=%0.4f(DF=%d),Pvalue1=%0.8e:ChiSq[%d]=%0.4f(DF=%d),Pvalue2=%0.8e: No Repeat detected\n",
	       Ymap->id,Xmap->id,orientation,score,logPV,K,ChiSquare[0],ChiSqDegree[0],Pvalue1,shift,ChiSquare[shift],ChiSqDegree[shift],Pvalue2);
	fflush(stdout);
      }
    }
    if((DEBUG && !(Pvalue1 > 1e-16)) || VERB>=2){/* This can happen for False Positive pairwise alignments : very rare if stringent -T -A -S thresholds are used */
      #pragma omp critical
      {
	if(VERB <= 1)
	  for(int shift = 1; shift <= maxshift; shift++){
	    assert(ChiSqDegree[shift] > 0);
	    printf("id=%lld,%lld,or=%d,score=%0.4f,logPV=%0.4f:K=%d:ChiSq[0]=%0.4f(DF=%d),Pvalue1=%0.8e:ChiSq[%d]=%0.4f(DF=%d),Pvalue2=%0.8e\n",
		   Ymap->id,Xmap->id,orientation,score,logPV,K,ChiSquare[0],ChiSqDegree[0],ChiSqPvalue(ChiSqDegree[0],ChiSquare[0]),shift,ChiSquare[shift],ChiSqDegree[shift],ChiSqPvalue(ChiSqDegree[shift],ChiSquare[shift]));
	  }

	double rsum = 0.0;
	for(int i = 0; i < K; i++){
	  if(outlier[i]){
	    printf("i=%d:Xfrag[i]=%0.3f,Yfrag[i]=%0.3f,outlier[i]=%d,outscore[i+1]=%0.4f,iscore[i+1]=%0.4f\n",
		   i,Xfrag[i],Yfrag[i],outlier[i],outscore[i+1],iscore[i+1]);
	    continue;
	  }
	  double y = Yfrag[i];
	  double x = Xfrag[i];
	  double vary = varF + var*(x+y);
	  if(QUADRATIC_VARIANCE)
	    vary += varR * (x*x + y*y);

	  double err = x-y;
	  if(DEBUG) assert(vary > 0.0);
	  rsum += err*err/vary;

	  printf("i=%d:Xfrag[i]=%0.3f,Yfrag[i]=%0.3f,outlier[i]=%d,outscore[i+1]=%0.6f,iscore[i+1]=%0.6f,rms=%0.4f,cumrms=%0.4f\n",
		 i,Xfrag[i],Yfrag[i],outlier[i],outscore[i+1],iscore[i+1],err*err/vary,rsum);
	}
	printf("varF=%0.6f,var=%0.6f\n",varF,var);
	if(!(Pvalue1 > 1e-16))
	  printf("WARNING: Pvalue1=%0.8e\n",Pvalue1);
	fflush(stdout);

	if(DEBUG) assert(isfinite(Pvalue1));
	//      assert(Pvalue1 > 1e-16);
      }
    }
    if(DEBUG>=2) assert(isfinite(Pvalue1));
  }
  
  delete [] Yfrag;
  delete [] outlier;

  if(VERB>=2 && result && !spots_filename[0]){
    #pragma omp critical
    {
      printf("Discarding pairwise alignment:id=%lld,%lld:repeat period=%d,A=%d,S=%0.2f,T=%0.2f\n",Gmap[mapid1]->id,Gmap[mapid2]->id,result,numpairs,score,logPV);
      fflush(stdout);
    }
  }

  return result;
}

void copy(Calign *dest, Calign *pinit, int score, int deepcopy){
  *dest = *pinit;/* first copy all scalars */

  dest->align2 = NULL;
  dest->Malign = NULL;
  dest->multicnt = dest->multimax = 0;
  dest->suppress = NULL;

  if(!deepcopy){
    dest->block_id = -2;
    return;
  }

  dest->sites1 = NULL;
  dest->sitesK1 = NULL;
  dest->sites2 = NULL;
  dest->iscore = NULL;
  dest->outscore = NULL;
    
  dest->block_id = -1;

  dest->allocated_size = pinit->numpairs;
  if(pinit->numpairs > 0){
    dest->numpairs = pinit->numpairs;
    dest->sites1 = new int[pinit->numpairs];
    dest->sites2 = new int[pinit->numpairs];
    for(int k= pinit->numpairs;--k>=0;){
      dest->sites1[k] = pinit->sites1[k];
      dest->sites2[k] = pinit->sites2[k];
    }
    if(pinit->sitesK1){
      dest->sitesK1 = new int[pinit->numpairs];
      for(int k=pinit->numpairs;--k>=0;)
	dest->sitesK1[k] = pinit->sitesK1[k];
    }
    if(score){
      if(pinit->iscore){
	dest->iscore = new FLOAT[pinit->numpairs+1];
	for(int k=pinit->numpairs+1;--k>=0;)
	  dest->iscore[k] = pinit->iscore[k];
      }
      if(pinit->outscore){
	dest->outscore = new FLOAT[pinit->numpairs+1];
	for(register int k=pinit->numpairs+1;--k>=0;)
	  dest->outscore[k] = pinit->outscore[k];
      }
    }
  }
}

void copyrev(Calign *dest, Calign *pinit, int score){
  *dest = *pinit;/* first copy all scalars */

  dest->align2 = NULL;
  dest->Malign = NULL;
  dest->multicnt = dest->multimax = 0;
  dest->suppress = NULL;
    
  dest->sites1 = NULL;
  dest->sitesK1 = NULL;
  dest->sites2 = NULL;
  dest->iscore = NULL;
  dest->outscore = NULL;

  dest->block_id = -1;

  int U = dest->allocated_size = pinit->numpairs;
  if(U > 0){
    Cmap *Ymap = YYmap[pinit->mapid1];
    Cmap *Xmap = XXmap[pinit->mapid2];
    int N = Ymap->numsite[0];
    int M = Xmap->numsite[0];
    //    FLOAT *Y = Ymap->site[0];
    //    FLOAT *X = Xmap->site[0];

    // correct Lij1,Rij1,Lij2,Rij2 and flip Lend,Rend
    dest->Lij1 = N+1 - pinit->Rij1;
    dest->Rij1 = N+1 - pinit->Lij1;
    dest->Lij2 = M+1 - pinit->Rij2;
    dest->Rij2 = M+1 - pinit->Lij2;
#if CALIGN_END==0
    dest->LijY = N+1 - pinit->RijY;
    dest->RijY = N+1 - pinit->LijY;
    dest->LijX = M+1 - pinit->RijX;
    dest->RijX = M+1 - pinit->LijX;

#endif
    dest->Lend = pinit->Rend;
    dest->Rend = pinit->Lend;

    dest->numpairs = U;
    dest->sites1 = new int[U];
    dest->sites2 = new int[U];

    for(int k= U; --k >= 0;){
      dest->sites1[k] = N + 1 - pinit->sites1[U-1-k];
    }
    if(pinit->sitesK1){
      dest->sitesK1 = new int[U];
      for(int k=U; --k >= 0;){
	dest->sites2[k] = M + 1 - pinit->sites2[U-1-k];
	int I = pinit->sites1[U-1-k];
	int K = pinit->sitesK1[U-1-k];
	dest->sites1[k] = N+1 - (I-K);
	dest->sitesK1[k] = K;
      }
    } else {
      for(int k= U; --k >= 0;){
	dest->sites1[k] = N + 1 - pinit->sites1[U-1-k];
	dest->sites2[k] = M + 1 - pinit->sites2[U-1-k];
      }
    }

    if(score){
      if(pinit->iscore){
	dest->iscore = new FLOAT[U+1];
	for(int k=U+1; --k >= 0;)
	  dest->iscore[k] = pinit->iscore[U-k];
      }
      if(pinit->outscore){
	dest->outscore = new FLOAT[U+1];
	for(register int k=U+1; --k >= 0;)
	  dest->outscore[k] = pinit->outscore[U-k];
      }
    }
  }

  dest->rev = 1 - dest->rev;
  if(RefSplitStitch>=2 && dest->stitch)
    dest->stitch = ((pinit->stitch & 1) ? 2 : 0) | ((pinit->stitch & 2) ? 1 : 0);
}


int num_align_blocks = 0;/* number of Calign_block(s) allocated */
Calign_block *align_block[MAX_ALIGNBLOCKS];

/* allocate a Calign_block to hold "num" alignments with a total of "sum" aligned sites. Append the pointer to align_blocks[0..num_align_blocks-1] and return it. */
Calign_block *align_blockalloc(size_t num, size_t sum, size_t Msiz, int AlignScore, int sitesK, int &block_id)
{
  if(VERB>=2){
    printf("align_blockalloc:time=%0.6f(wall time=%0.6f)\n",mtime(),wtime());
    fflush(stdout);
  }

  if(num_align_blocks >= MAX_ALIGNBLOCKS){
    printf("Exceeded maximum of %d alignment blocks\n",MAX_ALIGNBLOCKS);
    exit(1);
  }
  Calign_block *p = align_block[num_align_blocks] = new Calign_block(num, sum, Msiz, AlignScore, sitesK,num_align_blocks);
  block_id = num_align_blocks++;
  return p;
}

/* compact a Calign_block align_block[block_index] to only include alignment[block->start .. block->end-1] != NULL, all of which must point to memory in align_block[block_index].
   numpairs[] and numpairs_cum[] are pre-allocated arrays at least as large as (block->start - block->end)
   relocate the remaining "num" entries to alignment[start .. start + num - 1], where start <= block->start and alignment[start .. block->start - 1] are all NULL
   align_block[block_index] will be replaced with a reduced size memory block, and the original memory block deleted.

   returns num, the compacted number of alignments at alignment[start .. start + num - 1]
 */
size_t align_blockcompact(int block_index, size_t start, int *numpairs, size_t *numpairs_cum, int AlignScore, int numthreads, Calign **alignment)
{
  Calign_block *origblock = align_block[block_index];
  size_t origstart = origblock->start;
  size_t origend = origblock->end;

  /* remove alignments from alignment[start .. end-1] that are NULL*/
  size_t j = start, k = origstart;
  if(start == origstart){
    for(; k < origend; k++){
      if(!alignment[k]){
	k++;
	break;
      }
      j++;
    }
  }
  for(; k < origend; k++){
    if(!alignment[k])
      continue;
    if(DEBUG>=2) assert(j < k);
    if(DEBUG>=2) assert(alignment[j] == NULL);
    alignment[j] = alignment[k];
    alignment[k] = NULL;
    j++;
  }
  if(VERB>=2){
    printf("align_block[%d]: numaligns reduced from %lu to %lu:time=%0.6f(wall time=%0.6f)\n",block_index, origend - origstart, j-start, mtime(),wtime());
    fflush(stdout);
  }
  size_t end = j;

  size_t num = end - start;
  size_t numpairs_sum = 0;

  for(size_t i = 0; i < num; i++){
    Calign *p = alignment[start + i];
    numpairs[i] = p->numpairs;

    numpairs_cum[i] = numpairs_sum;
    numpairs_sum += numpairs[i];
  }
  numpairs_cum[num] = numpairs_sum;

  Calign_block *block = align_block[block_index] = new Calign_block(num, numpairs_sum, AlignScore, 0, block_index,0);
  block->start = start;
  block->end = end;

  long long numLL = num;

  size_t clen = (colors==1) ? 1 : colors + 1;

#pragma omp parallel for num_threads(numthreads) schedule(static,512)
  for(long long i = 0; i < numLL; i++){
    Calign *p = &block->block[i*clen];

    size_t nump = numpairs[i];
    size_t offset = 2 * numpairs_cum[i];
    p->block_id = block_index;
    p->allocated_size = nump;
    p->sites1 = &block->intblock[offset];
    p->sites2 = &block->intblock[offset+nump];
    p->sitesK1 = NULL;
    if(AlignScore){
      p->iscore = &block->scoreblock[offset];
      p->outscore = &block->scoreblock[offset+nump];
    } else
      p->iscore = p->outscore = NULL;
  }
  
  block->nextint = 2 * numpairs_cum[num];
  if(DEBUG && !(block->nextint <= block->intsiz)){
    printf("num=%lu,numpairs_cum[num]=%lu, block->nextint=%lu, block->intsiz=%lu\n",num,numpairs_cum[num],block->nextint, block->intsiz);
    fflush(stdout);
    assert(block->nextint <= block->intsiz);
  }
  if(AlignScore){
    block->nextscore = 2 * numpairs_cum[num];
    if(DEBUG && !(block->nextscore <= block->scoresiz)){
      printf("num=%lu,numpairs_cum[num]=%lu, block->nextscore=%lu, block->scoresiz=%lu\n", num, numpairs_cum[num], block->nextscore, block->scoresiz);
      printf("    block->nextint=%lu, block->intsiz=%lu\n", block->nextint, block->intsiz);
      fflush(stdout);
      assert(block->nextscore <= block->scoresiz);
    }
  }
  
#pragma omp parallel for num_threads(numthreads) schedule(static,512)  
  for(long long i = 0; i < numLL; i++){/* copy data from alignment[start + i] to block->block[i*clen] */
    Calign *q = alignment[start+i];
    Calign *p = &block->block[i*clen];

    p->mapid1 = q->mapid1;
    p->mapid2 = q->mapid2;
    p->score = q->score;
    p->logPV = q->logPV;
    p->trueoffset = q->trueoffset;
    p->trueoverlap = q->trueoverlap;
    p->orientation = q->orientation;
    p->scaleID = q->scaleID;
    p->trueTP = q->trueTP;
    p->numpairs = numpairs[i];// == q->numpairs

    p->Lend = q->Lend;// NEW10
    p->Rend = q->Rend;// NEW10
    p->fileid = q->fileid;// NEW10

    if(DEBUG>=2) assert(p->numpairs == q->numpairs);
    if(DEBUG>=2) assert(&p->sites2[0] == &p->sites1[p->numpairs]);
    if(DEBUG>=2) assert(&q->sites2[0] == &q->sites1[p->numpairs]);
    memcpy(&p->sites1[0], &q->sites1[0], sizeof(int) * p->numpairs * 2);
    //    memcpy(&p->sites2[0], &q->sites2[0], sizeof(int) * p->numpairs);
    if(AlignScore){
      if(DEBUG>=2) assert(&p->outscore[0] == &p->iscore[p->numpairs]);
      if(DEBUG>=2) assert(&q->outscore[0] == &q->iscore[p->numpairs]);
      memcpy(&p->iscore[0], &q->iscore[0], sizeof(FLOAT) * p->numpairs * 2);
      //      memcpy(&p->outscore[0], &q->outscore[0], sizeof(FLOAT) * p->numpairs);
    }
    
    alignment[start+i] = p;
  }

  delete origblock;

  return num;
}

/* free all Calign_block(s) */
void align_blockfree()
{
  for(int i = 0; i < num_align_blocks; i++)
    delete align_block[i];
  num_align_blocks = 0;
}

void alignfree()
{
  if(maxaligns > 0){
    if(alignment){
      for(register size_t i = 0; i < maxaligns; i++){
	if(VERB>=2){
	  printf("i=%lu/%lu:alignment[i]=%p",
		 i,maxaligns,alignment[i]);
	  fflush(stdout);
	  if(alignment[i]){
	    Calign *p = alignment[i];
	    printf(",alignment[i][0]:block_id=%d,mapid1=%d(id=%lld),mapid2=%d(id=%lld),or=%d:Malign=%p,multicnt=%d,multimax=%d\n",
		   p->block_id,p->mapid1,YYmap[p->mapid1]->id,p->mapid2,XXmap[p->mapid2]->id,p->orientation,p->Malign,p->Malign?p->multicnt:0,p->Malign?p->multimax:0);
	  }
	  printf("\n");
	  fflush(stdout);
	}
	if(alignment[i] && alignment[i]->block_id < 0){
	  delete [] alignment[i];
	  alignment[i] = 0;
	}
      }
      delete [] alignment;
      alignment = 0;
    }
    maxaligns = 0;
  }
  numaligns = 0;

  if(VERB/* HERE HERE HERE */){
    static long long VmSize = 0, VmRSS = 0, VmSwap = 0, VmHWM = 0;
    getmem(VmSize,VmRSS,VmSwap,&VmHWM);	    
    printf("Before align_blockfree():VmSize= %0.4f, VmRSS= %0.4f, VmHWM= %0.4f, VmSwap= %0.4f Gb : CPU time= %0.6f, wall time= %0.6f\n",VmSize * 1e-9, VmRSS * 1e-9, VmHWM * 1e-9, VmSwap * 1e-9, mtime(), wtime());
    fflush(stdout);
  }

  align_blockfree();
}

void expand_arrays(Calign *p, int num, int AlignScore, Calign_block *block, int block_index){
  if(num <= p->allocated_size) return;
  if(p->allocated_size > 0 && p->block_id < 0){
    delete [] p->sites1;
    delete [] p->sitesK1;
    delete [] p->sites2;
    delete [] p->iscore;
    delete [] p->outscore;
  }
  p->block_id = block_index;
  p->allocated_size = num;
  p->sites1 = block->int_alloc(num);
  p->sites2 = block->int_alloc(num);
  p->sitesK1 = NULL;
  if(AlignScore){
    p->iscore = block->FLOAT_alloc(num);
    p->outscore = block->FLOAT_alloc(num);
  } else
    p->iscore = p->outscore = NULL;
}

