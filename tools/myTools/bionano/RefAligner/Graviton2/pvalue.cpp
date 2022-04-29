#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

#include "parameters.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/pvalue.cpp 10584 2020-02-12 02:25:46Z tanantharaman $");

#define PV_SCALE 1 /* try scaling X vs Y to improve Pvalue */
#define PV_RES PVres /* correct Pvalue for resolution limit (res value) */

// Define lgamma() since it is not supported under windows
// natural log of the gamma function 
#define G 7
double P[G+2] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
		 771.32342877765313, -176.61502916214059, 12.507343278686905,
		 -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

double Lgamma(register double x)
{
  if(x < 0.5)
    return log(M_PI/sin(M_PI*x))-Lgamma(1.0-x);

  // Based on Lanczo's approximation for the Gamma function http://en.wikipedia.org/wiki/Lanczos_approximation 
  // Accuracy better than 13 decimal places vs Linux lgamma()(0.001 to 100.000), 
  // (theoretically 15 decimal places absolute accuracy)
  x -= 1.0;
  register double v = P[0];
  for(register int i = 1; i <= G+1; i++)
    v += P[i]/(x+i);
  register double t = x + (G + 0.5);
  return 0.5*log(2.0*M_PI) + (x+0.5)*log(t) - t + log(v);
}

/* compute -log10(Pvalue) for given alignment:

Pvalue = choose(2n+r+2,r) * pow(Rsq * e * PI/2 , n/2) / sqrt(n * PI),

where:

n = number of aligned intervals X[0..n-1], Y[0..n-1]
Rsq = "weighted average squared relative error" = (sum(i=0..n-1) (X[i]+Y[i]) * ((X[i]-Y[i])/(X[i]+Y[i]))^2 / (n*G)
G = "geometric mean of X[i]+Y[i]" = pow(prod(i=0..n-1) (X[i]+Y[i]),1/n) 
r = misaligned sites (either fp or fn)
e = base of natural logarithms = 2.71...
PI = 3.141592...
choose(N,K) = "N choose K" = gamma(N+1)/(gamma(K+1)*gamma(N-K+1))

Correction for rounding to nearest pixel : force error |X[i]-Y[i]| to be at least 0.12*(res*PixelLen)

Correction for limited resolution (PV_RES && res > 0.0) : increase Rsq to 1-(1-Rsq)exp(-2*res*PixelLen*Rsq/((Lambda+rres*PixelLen) * (1-Rsq)))
    Also the number mis-resolved sites in excess of expected amount is added to r (the number of misaligned sites) (if PV_RES >= 2)

Correction for internal outliers in alignment : If there are OutlierCnt outliers, compute
          OutlierMN = prod(i=0..OutlierCnt) M[i]*N[i]    (where M[i],N[i] are the number of site intervals in the ith outlier interval)
	  OutlierProd = prod(i=0..Outliercnt) (n + Outliercnt - 1 - i)
	  Multiply Pvalue by Correction Factor = 1 + OutlierProd * OutlierMN

Correction for 1 endoutlier : Multiply by 2
Correction for 2 endoutliers, with aligned length L sites, total length of smaller map N sites (not counting sites closer than res) : Multiply by (N - L + 1)

 */

double pvalue(FLOAT *Y,/* interval sizes for map Y */
	      FLOAT *X,/* interval sizes for map X */
	      int n,/* number of intervals in X and Y (excluding outliers) */
	      int fn, /* number of misaligned fn sites */
	      double fp, /* number of misaligned fp sites (plus a fraction of the mis-resolved sites if PVres >= 2) */
	      int OutlierCnt,/* number of internal outlier aligned segments */
	      double logOutlierMN, /* log of product of m x n over all outlier aligned segments */
	      int EndOutlierNL, /* multiplicative correction factor for endoutliers */
	      double Ylambda,
	      double resKB,
	      int verbose, /* display sub-results */
	      int Yid, int Xid, int orientation /* for debugging */
	      )
{
  if(n <=0 )
    return 0.0;

  if(DEBUG>=1+RELEASE && !(fn >= 0)){
    printf("pvalue:n=%d,fn=%d,fp=%0.1f\n",n,fn,fp);
    fflush(stdout);
    assert(fn >= 0);
  }
  if(DEBUG>=1+RELEASE)assert(fp >= 0.0);

  register double r = max(0,fn) + max(0.0, fp) /* WAS fn + fp */;
  register double Rsq = 0.0;
  register double logG = 0.0;
  //  double resKB = res * PixelLen;
  register double minvar = 0.12*resKB;
  if(VERB>=2 && verbose)
    printf("minerr=%0.4f,minvar=%0.6f\n",minvar,minvar*minvar);
  minvar *= minvar;

  for(register int i = 0; i < n; i++){
    register double sum = X[i]+Y[i];
    register double delta = X[i]-Y[i];
    delta *= delta;
    if(delta < minvar)
      delta = minvar;
    Rsq += delta/sum;
    logG += log(sum);
    if((VERB>=2 && verbose) || (DEBUG/* HERE >=2 */ && !isfinite(logG)) ){
      double nInv = 1.0/(double)(i+1);
      double logR = log(Rsq*nInv * M_PI * M_E * 0.5) - logG * nInv;
      printf("i=%d/%d:X[i]=%0.3f,Y[i]=%0.3f,sum=%0.3f,var=%0.6f:Rsq=%0.6f,logG=%0.6f(G=%0.6e,R=%0.6e)\n",i,n,X[i],Y[i],sum,delta,Rsq,logG,
	     exp(logG*nInv),exp(logR*0.5));
      fflush(stdout);
      if(DEBUG) assert(isfinite(logG));
    }
  }

  register double nInv = 1.0/(double)n;
  double logR = log(Rsq*nInv * M_PI * M_E * 0.5) - logG*nInv;
  double minKB = rres * 0.500;
  if(PV_RES && resKB > 0.0 && logR < 0.0){
    double R = exp(logR);
    if(DEBUG) assert(R < 1.0);
    double Rres = 1.0 - (1.0 - R)*exp(-2.0*resKB*R/((Ylambda+minKB)*(1.0-R)));
    if((VERB>=2 && verbose) || !(Rres <= 1.0)){
      printf("  logR=%0.8f -> %0.8f,R=%0.10e,resKB=%0.3f,Ylambda=%0.3f,minKB=%0.3f,Rres=%0.10e\n",
	     logR,log(Rres),sqrt(R),resKB,Ylambda,minKB,sqrt(Rres));
      fflush(stdout);
    }
    if(DEBUG) assert(Rres <= 1.0);
    logR = log(max(R,Rres));
  }
  double logP = 0.5*(n*(logR) - log(M_PI * n));
  register double logC = Lgamma(2*n+r+3.0)-Lgamma(r+1.0)-Lgamma(2*n+3.0);
  register double logI = 0.0;
  if(OUTLIER_PV && (OutlierCnt > 0 || EndOutlierNL > 1) ){/* add correction for outliers */
    double OutlierMN = exp(logOutlierMN);
    register double prod = OutlierMN, nlogI = 0.0;
    double m = n + OutlierCnt;
    if(prod <= 1e+300){
      for(int i = 0; i < OutlierCnt; i++)
	prod *= PVALUE_FIX ? m-i : m - i - 1.0;
      if(DEBUG && PVALUE_FIX) assert(prod >= 1.0);
      nlogI = log(PVALUE_FIX ? prod : 1.0 + prod);
    } 
    if(!isfinite(prod) || prod > 1e+300) {
      nlogI = logOutlierMN;
      for(int i = 0; i < OutlierCnt; i++)
	nlogI += log(PVALUE_FIX ? m-i : m-i-1.0);
    }
    if((VERB>=2 && verbose) || (DEBUG && !(isfinite(nlogI)))){      
      printf("WARNING:pvalue(Yid=%d,Xid=%d,orientation=%d):n=%d:fn=%d,fp=%0.1f,Outliers=%d:EndOutlierNL=%d,logMN=%0.6f,prod=%0.6e,logI=%0.6f,nLogI=%0.6f\n",
	     Yid,Xid,orientation,n,fn,fp,OutlierCnt,EndOutlierNL,logOutlierMN,prod,logI,nlogI);
      if(VERB>=2)
	for(register int i = 0; i < n; i++)
	  printf("i=%d/%d:X[i]=%0.3f,Y[i]=%0.3f\n",i,n,X[i],Y[i]);
      fflush(stdout);
    }
    if(DEBUG) assert(isfinite(nlogI));
    if(DEBUG && !(EndOutlierNL > 0)){
      printf("pvalue:EndOutlierNL=%d\n",EndOutlierNL);
      fflush(stdout);
      assert(EndOutlierNL > 0);
    }
    logI = nlogI + log((double)EndOutlierNL);
    if(DEBUG) assert(isfinite(logI));
  }
  if(PV_SCALE){/* compute optimal scaling of X vs Y */
    double var= 0.0, delvar = 0.0;
    for(int i = 0; i < n; i++){
      double sum = X[i] + Y[i];
      var += sum*sum;
      delvar += Y[i]*Y[i] - X[i]*X[i];
    }
    double c = delvar / var;
    double Rsq2 =0.0 ,logG2 = 0.0;
    for(int i = 0; i < n; i++){
      double x = X[i]*(1.0+c);
      double y = Y[i]*(1.0-c);
      double sum = x+y;
      double delta = x-y;
      delta *= delta;
      if(delta < minvar)
	delta = minvar;
      Rsq2 += delta/sum;
      logG2 += log(sum);
      if(DEBUG/* HERE >=2 */) assert(isfinite(logG2));
    }
    double logR2 = log(Rsq2*nInv * M_PI * M_E * 0.5) - logG2*nInv;
    if(PV_RES && resKB > 0.0 && logR2 < 0.0){
      double R = exp(logR2);
      if(DEBUG) assert(R < 1.0);
      double Rres = 1.0 - (1.0 - R)*exp(-2.0*resKB*R/((Ylambda+minKB)*(1.0-R)));
      if((VERB>=2 && verbose) || !(Rres <= 1.0)){
	printf("  logR2=%0.8f -> %0.8f,R=%0.8e,resKB=%0.8f,Ylambda=%0.8f,minKB=%0.8f,Rres=%0.8e\n",
	       logR2, log(Rres), sqrt(R),resKB,Ylambda,minKB,sqrt(Rres));
	fflush(stdout);
      }
      if(DEBUG) assert(Rres <= 1.0);
      logR2 = log(max(R,Rres));
    }
    double logP2 = 0.5*(n*(logR2) - log(M_PI * n)); // 0.5*(n*log(Rsq2 * M_PI * M_E * 0.5*nInv) - logG2 - log(M_PI * n));    
    if(VERB>=2 && verbose){
      printf(" n=%d:fp=%0.2f,fn=%d,Outliers=%d,logMN=%0.8f,LogPV=%0.8f:c=%0.8f,LogPV2=%0.8f",
	     n, fp,fn,OutlierCnt,logOutlierMN, -(logP+logC+logI)/log(10.0), c, -(logP2+logC+logI)/log(10.0));
      if(logP2 < logP)
	printf(":logR->%0.8f,logP->%0.8f,Rsq->%0.8f,logG->%0.8f\n",logR2,logP2,Rsq2,logG2);
      else
	printf("\n");
      fflush(stdout);
    }
    if(logP2 < logP){
      logR = logR2;
      logP = logP2;
      Rsq = Rsq2;
      logG = logG2;
    }
  }

  if((VERB && verbose) || (DEBUG && !(isfinite(logP) && isfinite(logC) && isfinite(logI)))){
    register double R = sqrt(Rsq * M_PI * M_E * 0.5 * nInv/exp(logG*nInv));
    printf("    n=%d:fp=%0.2f,fn=%d,resKB=%0.3f,Outliers=%d,logMN=%0.6f,EndNL=%d,logG=%0.4f(G=%0.4e),Rsq=%0.6e,R=%0.6e(%0.6e),LogP=%0.4f,LogC=%0.4f,LogI=%0.4f,LogPV=%0.4f\n",
	   n, fp,fn,resKB,OutlierCnt, logOutlierMN, EndOutlierNL, logG,exp(logG*nInv), Rsq, R,exp(logR*0.5),-logP/log(10.0),-logC/log(10.0),-logI/log(10.0),-(logP+logC+logI)/log(10.0));
    if(VERB/* HERE >=2 */){
      register double NRsq = 0.0;
      register double NlogG = 0.0;
      for(register int i = 0; i < n; i++){
	register double sum = X[i]+Y[i];
	register double delta = X[i]-Y[i];
	delta *= delta;
	if(delta < minvar)
	  delta = minvar;
	NRsq += delta/sum;
	NlogG += log(sum);
	printf("i=%d/%d:X[i]=%0.6f,Y[i]=%0.6f,sum=%0.6f,delta=%0.6f:NRsq=%0.8e,NlogG=%0.8e\n",
	       i,n,X[i],Y[i],sum,delta,NRsq,NlogG);
      }
    }
    fflush(stdout);
    if(DEBUG) assert(isfinite(logG));
    if(DEBUG) assert(isfinite(logP));
    if(DEBUG) assert(isfinite(logC));
    if(DEBUG) assert(isfinite(logI));
  }
  return max(0.0, -(logP+logC+logI)/log(10.0));
}

#define KDEBUG (DEBUG>=2) // Debug KSPvalue() calls

#define KS_FIX 1 /* reduce significance of KSPvalue() by reducing sample sizes used to the smaller of the two samples sizes times this value (if needed) : 
		    only used after computing the KS statistic */

/* Kolmogorov-Smirnov Two-Sample log(Pvalue) (approximate). NOTE : samples S1[0..n1-1] and S2[0..n2-1] are in ascending order */
double KSPvalue(int n1, double *S1, int n2, double *S2, int verb)
{
  /* first compute the Kolmogorov smirnov statistic */
  double D12 = 0.0;
  double F1 = 0.0,F2 = 0.0;
  int i1 = 0,i2 = 0;
  double Inv1 = 1.0/n1, Inv2 = 1.0/n2;

  if(KDEBUG){
    for(int i = 1; i < n1; i++){
      if(!(S1[i] >= S1[i-1])){
	printf("i=%d,n1=%d,S1[i-1]=%0.8e,S1[i]=%0.8e\n",i,n1,S1[i-1],S1[i]);
	fflush(stdout);
	assert(S1[i] >= S1[i-1]);
      }
    }
    for(int i = 1; i < n2; i++)
      assert(S2[i] >= S2[i-1]);
    assert(n1 > 0 && n2 > 0);
  }

  while(i1 < n1 && i2 < n2){

    if(S1[i1] < S2[i2])/* process point S1[i] next */
      F1 = (double)(++i1) * Inv1;
    else /* process point S2[i2] next */
      F2 = (double)(++i2) * Inv2;
    if(KDEBUG) assert(0 <= F1 && F1 <= 1.0);
    if(KDEBUG) assert(0 <= F2 && F2 <= 1.0);
    double val = fabs(F1-F2);
    if(KDEBUG) assert(0 <= val && val <= 1.0);
    if(val > D12)
      D12 = val;
  }
  if(KDEBUG) assert(i1 >= n1 || i2 >= n2);

  if(KDEBUG){
    if(i1 >= n1){
      while(i2 < n2){
	F2 = (double)(++i2) * Inv2;
	if(DEBUG) assert(0 <= F1 && F1 <= 1.0);
	if(DEBUG) assert(0 <= F2 && F2 <= 1.0);
	double val = fabs(F1-F2);
	if(DEBUG) assert(0 <= val && val <= 1.0);
	if(DEBUG) assert(val <= D12);
	if(val > D12)
	  D12 = val;
      }
    } else {
      while(i1 < n1){
	F1 = (double)(++i1) * Inv1;
	if(DEBUG) assert(0 <= F1 && F1 <= 1.0);
	if(DEBUG) assert(0 <= F2 && F2 <= 1.0);
	double val = fabs(F1-F2);
	if(DEBUG) assert(0 <= val && val <= 1.0);
	if(DEBUG) assert(val <= D12);
	if(val > D12)
	  D12 = val;
      }
    }
  }/* if(KDEBUG) */

  if(KS_FIX){
    int n = min(n1,n2) * KS_FIX;
    n1 = min(n,n1);
    n2 = min(n,n2);
  }

  /* next compute the Pvalue based on D12*sqrt(n1*n2/(n1+n2)) being approximately a Kolmogorov-Smirnov distribution */
  double s = D12 * sqrt((double)n1*n2/(n1+n2));
  double ret = 0.0;

  if(s < 1.25){/* compute 1 - sqrt(2*PI)/s * sum(k=1...) exp(-(2k-1)^2*PI^2/8s^2) */
    double u = -M_PI*M_PI/(8.0*s*s);
    double z = exp(u);
    double z2 = z*z;
    double z4 = z2*z2;
    double z8 = z4*z4;
    double z16 = z8*z8;
    double z24 = z16*z8;
    double z32 = z16*z16;
    double sum = z * (1.0 + (z8 * (1.0 + (z16 *(1.0 + z24 * (1.0 + z32))))));
    if(KDEBUG){
      double rsum = exp(u) + (exp(9.0*u) + (exp(25.0*u) + (exp(49.0*u) + exp(81.0*u))));
      assert(fabs(sum-rsum) < 1e-8);
    }
    sum *= sqrt(2.0*M_PI)/s;
    ret = (sum < 1e-8) ? -sum : log(1.0 - sum);
  } else {
    /* compute 2 * sum(k=1 ... ) (-1)^(k-1) exp(-2*(k*s)^2) */
    double u = -2.0*s*s;
    double z = exp(u);
    double z2 = z*z;
    double z3 = z2*z;
    double z4 = z2*z2;
    double z5 = z*z4;
    double z7 = z3*z4;
    double z9 = z5*z4;
    double sum = 1.0 - (z3 * (1.0 - z5 * (1.0 - z7 * (1.0 - z9))));
    if(KDEBUG){
      double rsum = 1.0 - (exp(3.0*u) - (exp(8.0*u) - (exp(15.0*u) - exp(24.0*u))));
      assert(fabs(sum-rsum) < 1e-8);
    }
    ret = u + log(2.0*sum);
  }
  if(VERB>=2 && verb && n1 > 0 && n2 > 0){
    printf("KSPvalue:n1=%d(%0.1f",n1,S1[0]);
    for(int i = 1; i < n1; i++)
      printf(",%0.1f", S1[i]);
    printf("),n2=%d(%0.1f",n2,S2[0]);
    for(int i = 1; i < n2; i++)
      printf(",%0.1f", S2[i]);
    printf("):D12=%0.6f,s=%0.6f,-log10pvalue=%0.6f\n",D12,s,-ret/log(10.0));
  }

  return ret;
}
