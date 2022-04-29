#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>

const char *SVN_ID =  (char *)"$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/test_generator.cpp 7559 2018-04-09 20:16:55Z tanantharaman $";

#ifndef DEBUG
#define DEBUG 1
#endif

#define CMAP_REF 1 /* Output CMAP_REF reference with -A option, instead of spotsfile */

#define SF_FIX 1 /* apply half the -sf variance to each site as a local peturbation (more physically realistic but estimated sf no longer exactly matches simulated value) */

#define VERB 1

#define MASK(n) (~((~(0U))<<(n))) // rightmost n bits set to 1 

int PIXELFRACTION = 1; /* allow fractional pixel site value */
int LINEAR = 1;/* assume linear genome : sample molecules that are shorter due to ends */

#define MAXCOLORS 2 /* maximum number of simultaneous probes */
#define MAXENZYMES 3 /* maximum number of simultaneous enzymes per probe color */
#define MAXPROBESET 64 /* maximum base pairs in probe */
#define MAXPROBES 1024 /* maximum probes per map (initial allocation : can increase) */

#define CHR_OFFSET 1000000000LL /* bp location offset as multiple of chromosome number */

#define CHR_WIDTH 100LL /* decimal field width of map chromosome number (in modulo notation) */
#define LOC_WIDTH 1000000000LL /* decimal field width of map Location (in modulo notation) */
#define LEN_WIDTH 10000000LL   /* decimal field width of map size (in modulo notation) */
#define FLAG_WIDTH 10LL /* decimal field width of map flags (3 boolean flags per decimal) */


/* compute location of condensed site formed from the K+1 consecutive sites Y[I-K ... I] */
static inline double Yc(register double *Y,
			register int I,
			register int K)
{
  return 0.5*(Y[I] + Y[I-K]);
}

static double Pr(register double y, register double resKB, register double IresSD)
{
  return 0.5*(1.0+erf((y - resKB)*IresSD));
}

static inline double min(register double A,register double B)
{
  return (A <= B) ? A : B;
}

static inline double max(register double A,register double B)
{
  return (A >= B) ? A : B;
}

static char **argv0;
static int argc0;

/* write a randome genome data file */
#define TABSIZE 1024
#define TABMASK 1023

static int
random_genome (double genome_length, double gc_fract, unsigned int seed)
{
  int gctable[TABSIZE];

  genome_length = floor ((genome_length + 2.0) / 4.0) * 4.0;

  if (isnan (gc_fract)) gc_fract = 0.0;
  if (gc_fract < 0.0) gc_fract = 0.0;
  if (gc_fract > 1.0) gc_fract = 1.0;
  if (genome_length < 0.0) genome_length = 0.0;

  fprintf (stderr, "Genome size (kb) -G%.3f\n", genome_length/1000.0);
  fprintf (stderr, "GC fraction -g%f\n", gc_fract);
  fprintf (stderr, "random number seed -R%u\n", seed);

  srandom (seed);

  {
    int i, a, c, g, t;
    double at, gc;
    gc = floor (gc_fract * TABSIZE);
    at = floor (TABSIZE - gc);
    g = gc / 2.0;
    c = gc - g;
    a = at / 2.0;
    t = at - a;
    fprintf (stderr, "a = %d c = %d g = %d t=%d\n", a, c, g, t);
    c += a;
    g += c;
    t += g;
    fprintf (stderr, "a = %d c = %d g = %d t=%d\n", a, c, g, t);
    if (t != TABSIZE) {
      fprintf (stderr, "t = %d\n", t);
      exit (1);
    }
    for (i = 0; i < a; i++) {
      gctable[i] = 0;
    }
    for ( ; i < c; i++) {
      gctable[i] = 1;
    }
    for ( ; i < g; i++) {
      gctable[i] = 2;
    }
    for ( ; i < t; i++) {
      gctable[i] = 3;
    }
  }

  while ( (genome_length -= 4.0) >= 0.0) {
    int byte;
    unsigned char uc;
    byte = ((gctable[random()&TABMASK]<<6) +
	    (gctable[random()&TABMASK]<<4) +
	    (gctable[random()&TABMASK]<<2) +
	    (gctable[random()&TABMASK]));
    uc = byte;
    fwrite (&uc, 1, 1, stdout);
  }
  return 0;
}
#undef TABSIZE
#undef TABMASK

/*
From statlib:

     REAL FUNCTION INVNOR(P)
                    NORmal distribution INVerse


                              Function


     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
     infinity to X of (1/SQRT(2*PI)) EXP(-U*U) dU is P


                              Arguments


     P --> The probability whose normal deviate is sought.
                    P is REAL


                              Method


     The  rational   function   on  page 95    of Kennedy  and  Gentle,
     Statistical Computing, Marcel Dekker, NY , 1980.


                              Note


     If P .lt. 1.0e-20 then INVNOR returns -10.0.
     If P .ge. 1.0 then INVNOR returns 10.0.

 */

static const double xnum[] = {
  -0.322232431088e0,
  -1.000000000000e0,
  -0.342242088547e0,
  -0.204231210245e-1,
  -0.453642210148e-4
};

static const double xden[] = {
  0.993484626060e-1,
  0.588581570495e0,
  0.531103462366e0,
  0.103537752850e0,
  0.38560700634e-2
};

static double
poly (const double *coef, int n, double x)
{
  int i;
  double term;
  term = coef[n-1];
  for (i = n-2; i >= 0; i--) {
    term = coef[i] + term * x;
  }
  return term;
}

static double 
inv_cumul_normal (double prob)
{
  double sign, y, z, result;
  if (prob < 1.0e-20) return -10.0;
  if (prob >= 1.0) return 10.0;
  if (prob <= 0.5) {
    sign = -1.0;
    z = prob;
  } else {
    sign = 1.0;
    z = 1.0 - prob;
  }
  y = sqrt (-2.0 * log(z));
  result = y + poly(xnum, 5, y) / poly (xden, 5, y);
  result *= sign;
  return result;
}

/* A uniformly distributed random number on [0,1) */
/* 0.0 <= urandom() < 1.0 */
static double urandom(void)
{
  register double r = random();
  return r / 2147483648.0; /* 2^31 */
}

/* force to lower case */
static char *lowercase(register const char *seq)
{
  char *result = (char *)malloc(strlen(seq)+1);
  register char *pt = result;
  for(;*seq;seq++)
    *pt++ = tolower(*seq);
  *pt = '\0';
  return result;
} 

/* force to upper case */
static char *uppercase(register const char *seq)
{
  char *result = (char *)malloc(strlen(seq)+1);
  if(!result){
    fprintf(stderr,"uppercase:malloc(%lu) failed\n",(unsigned long)(strlen(seq)+1));
    exit(1);
  }
  register char *pt = result;
  for(;*seq;seq++)
    *pt++ = toupper(*seq);
  *pt = '\0';
  return result;
} 

/* seqenzmatch[seqchar][enzchar] will be 1 IFF enzyme character
   "enzchar" will always match sequence character "seqchar" */
static int seqenzmatch[256][256];

static void nextcutinit(void)
{
  register int i,j;
  for(i=0;i<256;i++)
    for(j=0;j<256;j++)
      seqenzmatch[i][j] = 0;

  seqenzmatch[(int)'a'][(int)'a'] = 1;
  seqenzmatch[(int)'c'][(int)'c'] = 1;
  seqenzmatch[(int)'g'][(int)'g'] = 1;
  seqenzmatch[(int)'t'][(int)'t'] = 1;

  seqenzmatch[(int)'r'][(int)'r'] = 1;
  seqenzmatch[(int)'g'][(int)'r'] = 1;
  seqenzmatch[(int)'a'][(int)'r'] = 1;

  seqenzmatch[(int)'y'][(int)'y'] = 1;
  seqenzmatch[(int)'t'][(int)'y'] = 1;
  seqenzmatch[(int)'c'][(int)'y'] = 1;
  
  seqenzmatch[(int)'k'][(int)'k'] = 1;
  seqenzmatch[(int)'t'][(int)'k'] = 1;
  seqenzmatch[(int)'g'][(int)'k'] = 1;

  seqenzmatch[(int)'m'][(int)'m'] = 1;
  seqenzmatch[(int)'c'][(int)'m'] = 1;
  seqenzmatch[(int)'a'][(int)'m'] = 1;

  seqenzmatch[(int)'s'][(int)'s'] = 1;
  seqenzmatch[(int)'c'][(int)'s'] = 1;
  seqenzmatch[(int)'g'][(int)'s'] = 1;

  seqenzmatch[(int)'w'][(int)'w'] = 1;
  seqenzmatch[(int)'t'][(int)'w'] = 1;
  seqenzmatch[(int)'a'][(int)'w'] = 1;

  seqenzmatch[(int)'b'][(int)'b'] = 1;
  seqenzmatch[(int)'g'][(int)'b'] = 1;
  seqenzmatch[(int)'t'][(int)'b'] = 1;
  seqenzmatch[(int)'c'][(int)'b'] = 1;
  seqenzmatch[(int)'k'][(int)'b'] = 1;
  seqenzmatch[(int)'y'][(int)'b'] = 1;
  seqenzmatch[(int)'s'][(int)'b'] = 1;

  seqenzmatch[(int)'d'][(int)'d'] = 1;
  seqenzmatch[(int)'g'][(int)'d'] = 1;
  seqenzmatch[(int)'a'][(int)'d'] = 1;
  seqenzmatch[(int)'t'][(int)'d'] = 1;
  seqenzmatch[(int)'r'][(int)'d'] = 1;
  seqenzmatch[(int)'w'][(int)'d'] = 1;
  seqenzmatch[(int)'k'][(int)'d'] = 1;

  seqenzmatch[(int)'h'][(int)'h'] = 1;
  seqenzmatch[(int)'a'][(int)'h'] = 1;
  seqenzmatch[(int)'c'][(int)'h'] = 1;
  seqenzmatch[(int)'t'][(int)'h'] = 1;
  seqenzmatch[(int)'m'][(int)'h'] = 1;
  seqenzmatch[(int)'y'][(int)'h'] = 1;
  seqenzmatch[(int)'w'][(int)'h'] = 1;

  seqenzmatch[(int)'v'][(int)'v'] = 1;
  seqenzmatch[(int)'g'][(int)'v'] = 1;
  seqenzmatch[(int)'c'][(int)'v'] = 1;
  seqenzmatch[(int)'a'][(int)'v'] = 1;
  seqenzmatch[(int)'s'][(int)'v'] = 1;
  seqenzmatch[(int)'m'][(int)'v'] = 1;
  seqenzmatch[(int)'r'][(int)'v'] = 1;

  seqenzmatch[(int)'n'][(int)'n'] = 1;
  seqenzmatch[(int)'a'][(int)'n'] = 1;
  seqenzmatch[(int)'c'][(int)'n'] = 1;
  seqenzmatch[(int)'g'][(int)'n'] = 1;
  seqenzmatch[(int)'t'][(int)'n'] = 1;
  seqenzmatch[(int)'r'][(int)'n'] = 1;
  seqenzmatch[(int)'y'][(int)'n'] = 1;
  seqenzmatch[(int)'k'][(int)'n'] = 1;
  seqenzmatch[(int)'m'][(int)'n'] = 1;
  seqenzmatch[(int)'s'][(int)'n'] = 1;
  seqenzmatch[(int)'w'][(int)'n'] = 1;
  seqenzmatch[(int)'b'][(int)'n'] = 1;
  seqenzmatch[(int)'d'][(int)'n'] = 1;
  seqenzmatch[(int)'h'][(int)'n'] = 1;
  seqenzmatch[(int)'v'][(int)'n'] = 1;
}

/* return first char of first match of enzyme in sequence (or NULL if no match) */
static inline char *nextcut(char *seq, const char *enz)
{
  if(!seq)
    return NULL;

  return strstr(seq,enz);
}

static char comp[256];/* complement of "acgt" */

/* Invert (reverse) and complement a DNA sequence stored as string
   over alphabet "acgtrykmswbdhvn"
 */
static char *invcomp (const char *seq)
{
  char *result = (char *)malloc (strlen (seq) + 1);
  const char *scan = seq;
  char *rscan = result;

  while (*scan) 
    scan++;

  while (scan != seq) {
    switch(*--scan){
    case 'a': *rscan = 't';break;
    case 'c': *rscan = 'g';break;
    case 'g': *rscan = 'c';break;
    case 't': *rscan = 'a';break;
    case 'r': *rscan = 'y';break;
    case 'y': *rscan = 'r';break;
    case 'k': *rscan = 'm';break;
    case 'm': *rscan = 'k';break;
    case 'b': *rscan = 'v';break;
    case 'v': *rscan = 'b';break;
    case 'd': *rscan = 'h';break;
    case 'h': *rscan = 'd';break;
    case 's':
    case 'w':
    case 'n':
    default : *rscan = *scan;break;
    }
    rscan++;
  }
  *rscan = 0;
  return result;
}

typedef int intcmp (const void *,const void *);

static int randcutinc(int *p1, int *p2)
{
  return (*p1 > *p2) ? 1 : (*p1 < *p2) ? -1 : 0;
}

static int doubleinc(double *p1, double *p2)
{
  return (*p1 > *p2) ? 1 : (*p1 < *p2) ? -1 : 0;
}

#define EVENTTABSIZE 1000

static int random_poisson (double lambda)
{
  register int n;
  register double x, y;
  
  /* n is a random integer drawn from the Poisson distribution
     with mean lambda (and variance lambda). */

  /* 
     y is one term of Taylor series for exp (lambda).  n will be the
     number of terms of the Taylor series required to exceed x.
     (This amounts to a derivation of the Poisson PDF).

     Or, y is a factor of the Poisson PDF at n.  We draw a random
     number, and then sum the PDF until the sum exceeds the random
     number.

     PDF = exp(-lambda) lambda^n/ n!
     target = urandom ()

     y = lambda^n / n!
     x = urandom () * exp (lambda)

     x becomes negative when the cumulated PDF exceeds the random
     number.

     This fails for lambda > 709.7 because exp(709.8) is infinity.

  */
  if (lambda < 500.0) {
    x = urandom () * exp (lambda);
    for (y = 1.0, n = 0; ; y *= lambda / (++n)) {
      if ( (x -= y) <= 0.0)
	break;
      if ((!finite (x)) || (!finite (y)) || (n > 1000)) {
	fprintf (stderr, "random_poisson error n = %d x = %f y = %f lambda = %f\n", n, x, y, lambda);
	exit (1);
      }
    }
  } else {
    /* Can't draw a random poisson value if exp(lambda) is infinity. */
    /* So approximate it as a normal distribution with mean lambda */
    /* and standard deviation sqrt(lambda). */
    x = urandom ();
    n = floor (lambda + inv_cumul_normal (x) * sqrt(lambda));
  }

  return n;
}

/* Decide how many events given length of fragment */
/* and the rate per kb, generate that many uniformly */
/* distributed random events, and sort them. */

static int make_poisson_events (int len_bp, double events_per_kb, int *eventloc)
{
  int i, nevents; /* number of events */
  double lambda = len_bp * events_per_kb / 1000.0;

  nevents = random_poisson (lambda);

  if (nevents > EVENTTABSIZE) {
    fprintf (stderr, "Event table is too small, need %d, have %d.\n", nevents, EVENTTABSIZE);
    nevents = EVENTTABSIZE;
    exit (1);
  }
  for (i = 0; i < nevents; i++) {
    eventloc[i] = floor (urandom () * len_bp);
  }
  qsort (eventloc, (unsigned int)nevents, sizeof(int),(intcmp*) randcutinc);

  /*  if(fabs(events_per_kb-0.005)<=0.0001 ){
    printf("\nlen=%0.2f,F=%0.4f,lambda=%0.2f,nevents=%d\n",
	   len_bp*0.001,events_per_kb,lambda,nevents);
    printf("cuts=");
    for(i=0;i<nevents;i++)
      printf(" %0.2f",eventloc[i]*0.001);
    printf("\n");
    fflush(stdout);
    }*/

  return nevents;
}

/* 

   Here are two alternative Poisson random generators.  The first is
   slower than ours, the second is much faster.

   An alternate generator:
   http://ikpe1101.ikp.kfa-juelich.de/briefbook_data_analysis/node208.html
   
   A simple generator for random numbers taken from a Poisson
   distribution is obtained using this simple recipe: if $x_1, x_2,
   ...$ is a sequence of random numbers with uniform distribution
   between zero and one, $k$ is the first integer for which the product
   $x_1 x_2 ... x_{k+1} < e^{-\lambda}$.

   The alternate method draws many more random numbers than the method
   used here.

Here is the code from Mathematica 2.2.  It is much more efficient than
our method when the answer n is large, because it uses log(n) steps to
discover a value for n which gives a CDF which exceeds the required
probability, and then a binary search (another log(n) steps) to find
the exact answer.  It uses the Gamma function to directly calculate
the CDF rather than summing the PDF, but if n is say 2^20, it evaluates
the Gamma function about 40 times in lieu of evaluating the PDF 1000000
times.

PoissonDistribution/: PDF[PoissonDistribution[mu_], x_] :=
    If[!Negative[x], If[IntegerQ[x], Exp[-mu] mu^x / x!, 0], 0] /;
        PoissonPQ[mu]

PoissonDistribution/: CDF[PoissonDistribution[mu_], x_] :=
    Which[
      Negative[x], 0,
      x == Infinity, 1,
      True,
	GammaRegularized[Floor[x] + 1, mu, Infinity]
    ]	/; PoissonPQ[mu]

PoissonDistribution/: Quantile[PoissonDistribution[mu_], q_] :=
    With[{result = iPoissonQuantile[mu, q]},
        result /; result =!= Fail
    ] /; PoissonPQ[mu] && QuantileQ[q]
 
iPoissonQuantile[mu_, q_] := Module[{high, low, mid},
    If[q == 1, Return[Infinity]];
    If[N[Exp[-mu]] < q,
        high = 1;
        While[N[CDF[PoissonDistribution[mu], high]] < q, high *= 2];
        low = high/2;
        While[high - low > 1,
            mid = (high+low)/2;
            If[N[CDF[PoissonDistribution[mu], mid]] < q,
                low = mid,
                high = mid,
                low = high
            ]
        ]; high,
    (* else *)
        0,
    (* indeterminate *)
        Fail
    ]
]

PoissonDistribution/: Random[PoissonDistribution[mu_]] :=
       Quantile[PoissonDistribution[mu],Random[]]

*/

static void test_poisson (int argc, char **argv, const char *usage)
{
  int i, j, n, r, print, ix;
  double lambda;
  double x,sumx, sumxsqd, sqdsumx, numsamp, averagecoef, variancecoef;
  double average, variance, stddev;
  double meanmean, meanvar;

  if (argc != 4) {
    fprintf (stderr, usage, argv0[0]);
    fprintf (stderr, "missing parameters for %s\n", argv0[0]);
    exit (1);
  }
  n = atoi (argv[1]);
  r = atoi (argv[2]); /* number of repetitions */
  print = 1;
  if (n < 0) {
    n = -n;
    print = 0;
  }
  lambda = atof (argv[3]);
  
  fprintf (stderr, "Testing random Poisson generator.\n");
  fprintf (stderr, "Draw n = %d integers per sample and calculate mean and variance.\n", n);
  fprintf (stderr, "Draw r = %d samples.\n", r);
  if (print) {
    fprintf (stderr, "Print the integers as well as the statistics.\n");
  } else {
    fprintf (stderr, "Do not print the integers.\n");
  }
  fprintf (stderr, "Poisson rate lambda = %f.\n", lambda);
  meanmean = 0.0;
  meanvar = 0.0;
  for (j = 0; j < r; j++) {
    sumx = 0.0;
    sumxsqd = 0.0;
    sqdsumx = 0.0;
    numsamp = n;
    for (i = 0; i < n; i++) {
      x = ix = random_poisson (lambda);
      if (print) fprintf (stdout, "%d ", ix);
      sumx += x;
      x *= x;
      sumxsqd += x;
    }
    if (print) fprintf (stdout, "\n");
    sqdsumx = sumx * sumx;
    averagecoef = 1.0 / numsamp;
    variancecoef = 1.0 / (numsamp * (numsamp - 1.0));
    average = sumx * averagecoef;
    variance = (numsamp * sumxsqd - sqdsumx) * variancecoef;
    stddev = sqrt(variance);
    fprintf (stderr, "mean = %f ", average);
    fprintf (stderr, "standard deviation = %f ", stddev);
    fprintf (stderr, "variance = %f\n", variance);
    meanmean += average;
    meanvar += variance;
  }
  fprintf (stderr, "mean of means = %f\n", meanmean / r);
  fprintf (stderr, "mean of variances = %f\n", meanvar / r);
  exit (0);
}

class Cspot {
public:
  double site;/* observed bp location of nick site */
  double deres;/* range of sites that resolved into this site */
  long long refsite;/* true bp location of nick site (or -1 for false positive site) */
  long long refsiteL;/* If multiple nick sites contributed to this site, leftmost true bp location */
  long long color;/* color(Nickase) id of nick site (<= 0 for fragile site break) */
};

class Cmap {
public:
  int *numsite;/* numsite[c] : number of sites for color c in this map */
  Cspot **spots;/* spots[c][i=0..numsite[c]-1] : the information of site i in color c */

  long long startloc;/* start location in Reference Genome (in kb) For chimeric maps this is of the 1st fragment */
  long long endloc;/* end location in Reference Genome (in kb) For chimeric maps this is of the 1st part */
  long long len;/* length of map in Reference Genome (in kb) For chimeric maps this is the 1st part */
  long long startloc2;/* start location in Reference Genome (in kb) of 2nd chimeric fragment */
  long long endloc2;/* end location in Reference Genome (in kb) of 2nd chimeric fragment */
  long long len2;/* length in Reference Genome (in kb) of 2nd chimeric fragment */
  char *name;/* name of map */
  char chr;/* chromosome number */
  char flip;/* if map was flipped (for chimeric maps this happens after the two parts are combined, with the 2nd part flipped IF chimflip) */
  char chim;/* if map is chimeric */
  char chimflip;/* if 2nd part of chimeric map is flipped relative to the first part */
  Cmap(){
    flip = chim = chimflip = 0;
  }
  long long id(){/* encodes chromosome(2 digits <= 92),startloc(9 digits),(endloc-startloc)(7 digits),chipflip+chim+flip(3 bits as 1 digit) as a single 19 decimal number */
    long long val = chr;
    val = (val*LOC_WIDTH) + startloc % LOC_WIDTH;
    val = (val*LEN_WIDTH) + (endloc-startloc) % LEN_WIDTH;
    val = (val * 10LL) + (chimflip?4:0) + (chim?2:0) + (flip?1:0);
    if(DEBUG) assert(val >= 0);
    return val;
  }
};

static int spotcmp(Cspot *p1, Cspot *p2)
{
  return (p1->site > p2->site) ? 1 : (p1->site < p2->site) ? -1 : 0;
}

int main (int argc, char *argv [])
{
  argv0 = argv;
  argc0 = argc;

  double chimerism_probability = 0.0;
  /* false cut location table */
  int fcutloc[EVENTTABSIZE+1], nfcuts, fcut;
  char state1[128],state2[128],state3[128],state4[128],state5[128],state6[128];
  
  long long size = 0;/* space allocated for dnabin[]
			space for dna[] is 4*size+1 */
  long long ilen,ioffset;
  const char *acgt = "acgt";
  int acgt2bin[256];

  char *dnabin = 0;/* full genome (compressed) */
  char *dnafull = 0;/* full genome (uncompressed) */

  char *dna = 0;/* current molecule (1 or 2 segments of full genome) */
#if 0
  char *rdna = 0;/* reverse strand of current molecule */
#endif

  char *found1[MAXENZYMES] = {0,0,0};
  char *found2[MAXENZYMES] = {0,0,0};
  long long lastcut, thiscut, thiscut1[MAXENZYMES] = {0LL,0LL,0LL},thiscut2[MAXENZYMES] = {0LL,0LL,0LL};// thiscutf;
  register long long i, j, k; 
  double r, r1,r2,r3;
  double len, len2, offset, offset2;
  long long tlen,tlen2, toffset,toffset2;
  double noise;
  double rubberscale = 1.0;
  double noise_stddev;
  double PixelLen = 500.0;
  int Stealth = 0;/* If 1 : do not output true locations of nanomaps */
  char namebuf[1024];

  const char *usage = 
    "Usage: %s -G<genomesizekb>[_<offsetkb>] -N<nummol> -c<coverage> -M<meanmolkb> -m<stddevmolkb> \\\n"
    "-minlen<kb> -chr<int>\\\n"
    "-Eacgtacgt -FP<fp/100kb> -FN<fn rate> -R<randomseed> -res<res> -resSD<resSD> -A<reference>\\\n"
    "-C<chimerism rate> -sf<sf> -sd<sd> -sr<sr> -se<se> genome_file [... genome_fileN]\n"
    "  -G Genome size and optional offset in kb, selects specific region from genome file\n"
    "  -N Number of molecules\n"
    "  -c Genome Coverage achieved by output molecules\n"
    "  -M mean Molecule size (kb)\n"
    "  -m standard deviation Molecule size (kb) assuming Gaussian length distribution\n"
    "  -minlen Minimum Molecule size (kb) assuming exponential length distribution\n"
    "  -chr Chromosome number (added to map id and offset)\n"
    "  -colors Number of nicking enzymes (must be 1 or 2)\n"
    "  -C Chimerism probability.  If a random number is less than this\n"
    "     then the molecule will be made by joining sequence from two\n"
    "     parts of the genome.\n"
    "  -p pixel size in bp\n"
    "  -E Enzyme sequence recognition pattern(s), multiple Nickase patterns separated by comma (one option for each color probe used)\n"
    "  -FP False Positive density per 100 kb\n"
    "  -FN False Negative rate (probability)\n"
    "  -res site resolution in pixels\n"
    "  -resSD site resolution standard deviation in pixels\n"
    "  Interval sizing stddev is specified as a function of interval size x in kb:\n"
    "    variance of interval sizing = sf^2 + sd^2 x + sr^2 x^2\n"
    "  -sf specifies a sizing error in kb independent of interval size\n"
    "  -sd specifies a sizing error in kb^1/2 which scales variance by interval size (can be negative)\n"
    "  -sr specifies a relative sizing error\n"
    "  -se specifies an additional sizing error as fraction of unresolved interval\n"
    "       NOTE : If there is more than one Nickase, then the -FP through -se options will have multiple values separated by underscore\n"
    "  Options for rubberbanding type sizing error: (Example)\n"
    "    -r0.05 Scale all interval sizes by normal(1,0.05) to reflect error in scaling estimate\n"
    "  -rnd round site to nearest pixel\n"
    "  -R random number seed\n"
    "  -A reference filename (outputs reference as .cmap file)\n" // "  -A reference filename (outputs reference as .spots file)\n"
    "  -stealth do not output true locations of maps\n"
    "  -D specifies deterministic Fragile site interval size in kb (along with Fragile rate range over 10 bins, if 2 additional values are specified seperated by spaces)\n"
    "  -FR<mean>_<sd> specifies Fragile sites probability as a function of opposite strand nicking distance : Uses a standard normal function with specified <mean> and <sd>\n"
    "  -KN<Prob>_<mean>_<SD> specifies Knots between sites with probability <Prob> and size distributed as log-Normal with specified <mean> and <SD>\n"
    "  -ST<Prob>_<SD> specifies stitching error between sites with probability <Prob> and size distributed as Normal with mean zero and specied <SD>\n"
    "Specifying nonzero -g option causes program to write random genome file:\n"
    "Write random genome:  -Ggenomesizekb -ggcfrac -Rrandomseed > random_genome\n"
    "  -G Genome size (kb)\n"
    "  -g gc fraction e.g. -g0.6\n"
    "  -R random number seed\n"
    "Specifying -fasta causes a fasta sequence file to read in and written out in compact binary format\n"
    "Test poisson distribution function: -poisson draws repeats lambda\n"
    ;

  /* common parameters */
  double genome_length = 1000 * 1000 * 1000.0; /* -G4000000 (1GB) */
  double genome_offset = 0.0;

  double FP[MAXCOLORS] = {0.0, 0.0}; /* -FP0.0 */
  double FN[MAXCOLORS] = {0.2, 0.2}; /* -FN0.2_0.2 */
  double sf[MAXCOLORS] = {0.0, 0.0}; /* -sf0.0 */
  double sd[MAXCOLORS] = {0.0, 0.0}; /* -sd0.0 */
  double sr[MAXCOLORS] = {0.0, 0.0}; /* -sr0.0 */
  double se[MAXCOLORS] = {0.0, 0.0}; /* -se0.0 */
  double res[MAXCOLORS] = {0.0,0.0};/* -res0.0 : resolution in pixels */
  double resSD[MAXCOLORS] = {0.0,0.0};/* -resSD0.0 : resolution SD in pixels */
  double rubberband = 0.0;/* -r0.0 */

  double DetFragile = 0.0;/* create a deterministic fragile site every this many bases */
  double DetFragileLevel1 = 0.50;/* fraction of molecules that are fragile range from DetFragileLevel1 to DetFragileLevel2 (unformly distributed based on 10 bins) */
  double DetFragileLevel2 = 1.00;/* Bin is determined by (location % (DetFragile*10))/DetFragile : Bin 0 uses DetFragileLevel1 and Bin 9 uses DetFragileLevel1 + (9/10)*(DetFragileLevel2-DetFragileLevel1) */

  double FragileMean = 0.0;/* If > 0 : Use a fragile site probability as a function of opposite strand nicking distance : Uses standard normal function with specified Mean and SD */
  double FragileSD = 0.0;

  double KnotProb = 0.0;/* knot outlier probability (per site interval). Knot is treated as deletion in middle of site interval */
  double KnotMean = 25.0;/* mean knot size in kb */
  double KnotSD = 10.0;/* SD of knot size in kb */

  double StitchProb = 0.0;/* Stitch error probability (per site interval). Stitch is treated as duplication (+ve) or overlap (-ve) in middle of site interval */
  double StitchSD = 0.0;/* SD of stitch error in kb */

  /* optical map dataset parameters */
  int numcolors = 0;
  int numenzymes[MAXCOLORS] = {1,1};
  char *restrict[MAXCOLORS][MAXENZYMES];
  char *invtrict[MAXCOLORS][MAXENZYMES];
  int nmols = 0; /* -N0 */
  double coverage = 0.0;/* -c0.0 */
  double mean = 4000*1000.0; /* -M4000 */
  double stddev = 0.0; /* -m0 */
  double minlen = 0.0; /* -minlen 100.0 */
  int chr = 0;/* chromosome id */

  /* random genome parameters */
  double gc_fract = 0.0;/* -g options */

  int seed = 1; /* -R1 */
  char *spotsfile = 0; /* -A */
  FILE *fp=NULL;
  int fastafile = 0;

  if (argc == 1) {
    fprintf (stderr, usage, argv0[0]);
    fprintf (stderr, "at least one option is required.\n");
    exit (1);
  }

  int plen;

  char buf[PATH_MAX];

  argc--; argv++;
  while (argc > 0 && argv[0][0]=='-') {
    char *pt, *qt;

    switch (argv[0][1]) {
    case '-':
      fprintf (stderr, usage, argv0[0]);
      fprintf (stderr, "unknown option: %s\n", argv[0]);
      exit (1);
    case 'A':
      if(argv[0][2]){
	sprintf(buf,"%s.%s", &argv[0][2], CMAP_REF ? "cmap" : "spots");
	spotsfile = strdup(buf);
      }
      break;
    case 'c':
      plen = strlen("-chr");
      if(!strncmp(argv[0],"-chr",plen)){
	if(!argv[0][plen] || !isdigit(argv[0][plen])){
	  fprintf (stderr, usage, argv0[0]);
	  fprintf (stderr, "unknown option: %s\n", argv[0]);
	  exit (1);
	}
	chr = atoi(&argv[0][plen]);
	if(chr > 92 || chr < 0){
	  fprintf(stderr, "-chr value must be between 0 and 92\n");
	  exit(1);
	}
	break;
      }

      if(!argv[0][2] || !isdigit(argv[0][2])){
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      coverage = atof(&argv[0][2]);
      break;
    case 'D':
      DetFragile = argv[0][2] ? atof(&argv[0][2])*1000.0 : DetFragile;
      DetFragileLevel1 = DetFragileLevel2 = 1.0;
      if(argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-'){
	DetFragileLevel1 = atof(&argv[1][0]);
	DetFragileLevel2 = atof(&argv[2][0]);
	argv += 2;
	argc -= 2;
      }
      break;
    case 'g':
      gc_fract = argv[0][2] ? atof(&argv[0][2]) : gc_fract;
      break;
    case 'G':
      genome_length = 1000 * 1000 * 1000.0;
      genome_offset = 0.0;
      if(argv[0][2]){
	if((pt = strchr(&argv[0][2],'_')))
	  genome_offset = atof(&pt[1])*1000.0;
	genome_length = atof(&argv[0][2])*1000.0;
	genome_offset = floor(genome_offset*0.25)*4.0;
	genome_length = floor(genome_length*0.25)*4.0;
      } 
      break;
    case 'K':
      switch(argv[0][2]){
      case 'N':
	// HERE HERE : continue KN implementation (similar to -FR)

	break;
      default:
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit(1);
      }
      break;
    case 'N':
      if(!argv[0][2] || !isdigit(argv[0][2])){
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      nmols = atof(&argv[0][2]);
      break;
    case 'M':
      if(!argv[0][2] || !isdigit(argv[0][2])){
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      mean = atof(&argv[0][2])*1000.0;
      break;
    case 'm':
      plen = strlen("-minlen");
      if(!strncmp(argv[0],"-minlen",plen)){
	if(!argv[0][plen] || !isdigit(argv[0][plen])){
	  fprintf(stderr, usage, argv0[0]);
	  fprintf(stderr, "unknown option: %s\n",argv[0]);
	  exit(1);
	}
	minlen = atof(&argv[0][plen]) * 1000.0;
	break;
      }
      if(!argv[0][2] || !isdigit(argv[0][2])){
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      stddev = atof(&argv[0][2])*1000.0;
      break;
    case 'E':
      if(numcolors >= MAXCOLORS){
	fprintf(stderr,"Exceeded maximum probes (colors)=%d\n",MAXCOLORS);
	exit(1);
      }
      if(!argv[0][2]){
	fprintf(stderr,"-E option must be followed (without spaces) by Nickase sequence recognition pattern(s), multiple patterns seperated by commas\n");
	exit(1);
      }
      pt = &argv[0][2];
      {
	int numenz = 0;
	while((qt = strchr(pt,','))){
	  *qt++ = '\0';
	  restrict[numcolors][numenz] = lowercase(pt);
	  invtrict[numcolors][numenz] = invcomp(restrict[numcolors][numenz]);
	  numenz++;
	  pt = qt;
	}
	restrict[numcolors][numenz] = lowercase(pt);
	invtrict[numcolors][numenz] = invcomp(restrict[numcolors][numenz]);
	if((numenzymes[numcolors++] = ++numenz) > MAXENZYMES){
	  fprintf(stderr,"Exceeded maximum enzymes (%d) for probe color %d\n",MAXENZYMES, numcolors);
	  exit(1);
	}
      }
      break;
    case 'F':
      if(strlen(argv[0]) <= 4 && !isdigit(argv[0][3])){
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      switch(argv[0][2]){
      case 'R':
        pt = &argv[0][3];
        FragileMean = strtod(pt,&qt);
	if(qt == pt || !(*qt=='_') || FragileMean < 0.0){
	  printf("-FR first value %s is not valid (must be >= 0.0 and be followed by underscore)\n",pt);
	  exit(1);
	}
	pt = qt+1;
	
	FragileSD = strtod(pt,&qt);
	if(qt == pt || !(*qt == 0 || isspace(*qt)) || FragileSD <= 0.0){
	  printf("-FR 2nd value %s i not valid (must be > 0.0)\n",pt);
	  exit(1);
	}
        break;        
      case 'P':
	pt = &argv[0][3];
	for(int i = 0; i < numcolors; i++){
	  FP[i] = strtod(pt,&qt) * 0.01;
	  if(qt == pt || !(*qt == 0 || *qt=='_') || FP[i] < 0.0){
	    printf("FP[%d] value %s is not valid\n",i,qt);
	    exit(1);
	  }
	  pt = qt;
	  if(*pt == '_')
	    pt++;
	}
	break;
      case 'N':
	pt = &argv[0][3];
	for(int i = 0; i < numcolors; i++){
	  FN[i] = strtod(pt,&qt);
	  if(qt == pt || !(*qt == 0 || *qt=='_') || FN[i] < 0.0){
	    printf("FN[%d] value %s is not valid\n",i,qt);
	    exit(1);
	  }
	  pt = qt;
	  if(*pt == '_')
	    pt++;
	}
	break;
      default:
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      break;
    case 'C':
      if(!argv[0][2] || !isdigit(argv[0][2])){
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      chimerism_probability = atof(&argv[0][2]);
      break;
    case 's':
      if(!strcmp(argv[0], "-stealth")){
	Stealth = 1;
	break;
      }
      if(strlen(argv[0]) <= 4 && !isdigit(argv[0][3])){
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      switch (argv[0][2]) {
      case 'f':
	pt = &argv[0][3];
	for(int i = 0; i < numcolors; i++){
	  sf[i] = strtod(pt,&qt);
	  if(qt == pt || !(*qt == 0 || *qt=='_') || sf[i] < 0.0){
	    printf("sf[%d] value %s is not valid\n",i,qt);
	    exit(1);
	  }
	  pt = qt;
	  if(*pt == '_')
	    pt++;
	}
	break;
      case 'd':
	pt = &argv[0][3];
	for(int i = 0; i < numcolors; i++){
	  sd[i] = strtod(pt,&qt);
	  if(qt == pt || !(*qt == 0 || *qt=='_')){
	    printf("sd[%d] value %s is not valid\n",i,qt);
	    exit(1);
	  }
	  pt = qt;
	  if(*pt == '_')
	    pt++;
	}
	break;
      case 'r':
	pt = &argv[0][3];
	for(int i = 0; i < numcolors; i++){
	  sr[i] = strtod(pt,&qt);
	  if(qt == pt || !(*qt == 0 || *qt=='_') || sr[i] < 0.0){
	    printf("sr[%d] value %s is not valid\n",i,qt);
	    exit(1);
	  }
	  pt = qt;
	  if(*pt == '_')
	    pt++;
	}
	break;
      case 'e':
	pt = &argv[0][3];
	for(int i = 0; i < numcolors; i++){
	  se[i] = strtod(pt,&qt);
	  if(qt == pt || !(*qt == 0 || *qt=='_') || se[i] < 0.0){
	    printf("se[%d] value %s is not valid\n",i,qt);
	    exit(1);
	  }
	  pt = qt;
	  if(*pt == '_')
	    pt++;
	}
	break;
      default:
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      break;
    case 'R':
      if(!argv[0][2] || !isdigit(argv[0][2])){
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      seed = atoi(&argv[0][2]);
      break;
    case 'f':
      if (!strncmp(argv[0],"-fasta",6)){
	fastafile = 1;

	/*	fprintf(stderr,"fastafile=%d,argc=%d\n",fastafile,argc);
	  fflush(stderr);*/

	break;
      }
      fprintf (stderr, usage, argv0[0]);
      fprintf (stderr, "unknown option: %s\n", argv[0]);
      exit (1);
    case 'p':
      if(!strcmp(argv[0],"-poisson")){
	test_poisson(argc,argv,usage);
	exit(0);
      }
      if(!argv[0][2] || !isdigit(argv[0][2])){
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      PixelLen = atof(&argv[0][2]);
      break;
    case 'r':
      if(!strncmp(argv[0],"-resSD",6)){
	pt = &argv[0][6];
	for(int i = 0; i < numcolors; i++){
	  resSD[i] = strtod(pt,&qt);
	  if(qt == pt || !(*qt == 0 || *qt=='_') || resSD[i] < 0.0){
	    printf("resSD[%d] value %s is not valid\n",i,qt);
	    exit(1);
	  }
	  pt = qt;
	  if(*pt == '_')
	    pt++;
	}
	break;
      }
      if(!strncmp(argv[0],"-res",4)){
	pt = &argv[0][4];
	for(int i = 0; i < numcolors; i++){
	  res[i] = strtod(pt,&qt);
	  if(qt == pt || !(*qt == 0 || *qt=='_') || res[i] < 0.0){
	    printf("res[%d] value %s is not valid\n",i,qt);
	    exit(1);
	  }
	  pt = qt;
	  if(*pt == '_')
	    pt++;
	}
	break;
      }
      if(!strncmp(argv[0],"-rnd",4)){
	PIXELFRACTION = 0;
	break;
      }
      if(!argv[0][2] || !isdigit(argv[0][2])){
	fprintf (stderr, usage, argv0[0]);
	fprintf (stderr, "unknown option: %s\n", argv[0]);
	exit (1);
      }
      rubberband = atof(&argv[0][2]);
      break;
    default:
      fprintf (stderr, usage, argv0[0]);
      fprintf (stderr, "unknown option: %s\n", argv[0]);
      exit (1);
    }
    argv++;argc--;
  }

  if(SF_FIX){
    for(int i = 0; i < numcolors; i++)
      if(sd[i] < 0.0){
	printf("i=%d:sd[i]= %0.6f: must be >= 0.0\n",i,sd[i]);
	exit(1);
      }
  } else {
    for(int i = 0; i < numcolors; i++)
      if(sd[i] < -sqrt(2.0*sf[i]*sr[i])){
	printf("i=%d:sd[i]= %0.6f, sf[i]= %0.6f, sr[i]= %0.6f : sd[i] must be at least -sqrt(2*sfF[i]*sr[i]) = %0.6f\n",
	       i,sd[i],sf[i],sr[i],-sqrt(2.0*sf[i]*sr[i]));
	exit(1);
      }
  }

  if(DetFragile && FragileMean > 0.0){
    fprintf(stderr,"Cannot use both -D and -FR options : ignoring -D option\n");
    fflush(stderr);
    DetFragile = 0;
  }

  if (gc_fract != 0.0)
    return random_genome (genome_length, gc_fract, seed);

  if(minlen > 0.0 && stddev > 0.0){
    fprintf(stderr, "Cannot specify both -minlen<molkb> and -m<stddevmolkb>\n");
    exit(1);
  }
  if(coverage > 0.0 && nmols > 0){
    fprintf(stderr, "Cannot specify both -N<nmols> and -c<coverage>\n");
    exit(1);
  }

  nextcutinit();

  /* replace stdin by input file argv[0] */
  if( argc != 1 ){
    fprintf(stderr,usage,argv0[0]);
    fprintf(stderr,"missing input file (genome)\n");
    exit(1);
  }
  if((fp=fopen(argv[0],"r"))==NULL){
    fprintf(stderr,"cannot open file %s for input\n",argv[0]);
    exit(1);
  }
  (void)fclose(stdin);
  *stdin = *fp;

  /* initialize comp[0..255] so that comp['a'] == 't' etc */
  for(j = 0; j < 256; j++)
    comp[j] = j;
  comp[(unsigned int)'a'] = 't';
  comp[(unsigned int)'c'] = 'g';
  comp[(unsigned int)'g'] = 'c';
  comp[(unsigned int)'t'] = 'a';
  comp[(unsigned int)'r'] = 'y';
  comp[(unsigned int)'y'] = 'r';
  comp[(unsigned int)'k'] = 'm';
  comp[(unsigned int)'m'] = 'k';
  comp[(unsigned int)'b'] = 'v';
  comp[(unsigned int)'v'] = 'b';
  comp[(unsigned int)'d'] = 'h';
  comp[(unsigned int)'h'] = 'd';

  if(fastafile){
    long long tlen = 0;
    char buf[BUFSIZ];

    fprintf (stderr, "Reading fasta file %s & converting to compact binary format\n",argv[0]);
    fflush(stderr);

    if(fseek(stdin,0,SEEK_SET)){
      fprintf(stderr,"fseek of stdin failed\n");
      exit(1);
    }
    
    while(fgets(buf,BUFSIZ,stdin)!=NULL){
      size_t linesiz;

      if(buf[0]=='>')/* comment line */
	continue;

      linesiz = strlen(buf)-1;
      if(DEBUG) assert(buf[linesiz]=='\n');
      buf[linesiz] = '\0';
      while(linesiz>0 && isspace(buf[linesiz-1]))
	buf[--linesiz] = '\0';

      if((tlen += linesiz) >= size){
	size *= 2;
	if(size < ((int)((tlen+3)/4))*4)
	  size = ((int)((tlen+3)/4))*4;
	dna = (char *)realloc(dna,size+1);
	if(!dna){
	  fprintf(stderr,"realloc of %lld bytes for dna[] failed\n", size+1);
	  exit(1);
	}
      }
      memcpy(&dna[tlen-linesiz],buf,linesiz+1);
    }

    ilen = (tlen+3)/4;
    dnabin = (char *)realloc(dnabin,ilen+1);    
    if(!dnabin){
      fprintf(stderr,"realloc of %lld bytes for dnabin[] failed\n",ilen+1);
      exit(1);
    }
    
    for(i=tlen;i<ilen*4;i++)
      dna[i] = 'A';
    tlen = ilen*4;
    dna[ilen*4] = '\0';

    /* convert dna[0..tlen-1] to dnabin[0..ilen-1] */
    for(j=ilen;--j>=0;)
      dnabin[j] = 0;
    for(j=0;j<256;j++)
      acgt2bin[j] = -1;
    for(j=0;j<4;j++)
      acgt2bin[(int)acgt[j]] = j;
    acgt2bin[(int)'n'] = 0;/* treat all N as A */
    
    for(j=0;j<tlen;j++){
      k = acgt2bin[tolower((int)dna[j])];
      if(!(k>=0 && k<=3)){
	fprintf(stderr,"position %lld:base not recognized:\n%c\n", j,dna[j]);
	fflush(stderr);
	if(DEBUG) assert(k>=0 && k<=3);
      }
      dnabin[j/4] |= k << (6-(j%4)*2); 
    }

    /* write out dnabin[0..ilen-1] to stdout */
    fprintf(stderr,"Writing compacted sequence of %lld bases\n",ilen*4);
    fflush(stderr);
    if(ilen != (int)fwrite(dnabin,1,ilen,stdout)){
      fprintf(stderr,"write of %lld bytes to stdout failed\n",ilen);
      exit(1);
    }
    fflush(stdout);
    exit(0);
  }

  genome_length = floor ((genome_length + 2.0) / 4.0) * 4.0;
  if (genome_length < 0.0) genome_length = 0.0;

  if(genome_length > 2000*1000*1000.0){
    printf("Genome length larger than 2GB not supported : please break up into smaller contigs\n");
    exit(1);
  }

  int maxsites = MAXPROBES;
  Cspot *sites = new Cspot[maxsites];
  double *Y = new double[maxsites];
  int numsites;

  /* compute average interval size for specified enzyme over genome_length */
  if(nmols > 0 || coverage > 0.0)
    fprintf (stderr, "Reading genome from %s\n",argv[0]);

  if((ilen = (len = genome_length) / 4.0) > size){
    size = ilen;
    if(!(dnabin = (char *)realloc(dnabin,size))){
      fprintf(stderr,"realloc(dnabin,%lld) failed\n",size);
      exit(1);
    }
    if(!(dnafull = (char *)realloc(dnafull,4*size+1))){
      fprintf(stderr,"realloc(dnafull,%lld) failed\n",4*size+1);
      exit(1);
    }
    if(!(dna = (char *)realloc(dna,4*size+1))){
      fprintf(stderr,"realloc(dna,%lld) failed\n",4*size+1);
      exit(1);
    }
#if 0
    if(!(rdna = (char *)realloc(rdna,4*size+1))){
      fprintf(stderr,"realloc(rdna,%lld) failed\n",4*size+1);
      exit(1);
    }
#endif
  }
  if(fseek(stdin,0,SEEK_SET)){
    fprintf(stderr,"fseek of stdin failed\n");
    exit(1);
  }
  if(DEBUG) assert(size >= 0);
  if(size != (i = fread(dnabin,1,size,stdin))){
    if(i<=0){
      fprintf(stderr,"fread of %lld bytes failed(i=%lld)\n",size,i);
      exit(1);
    } else {
      if(nmols > 0 || coverage > 0.0)
	fprintf(stderr,"genome size read in is %lld bases\n",i*4);
      ilen = size = i;
      len = genome_length = size*4;
    }
  } else if(nmols>0 || coverage > 0.0)
    fprintf(stderr,"read in first %lld bases of genome\n",i*4);

  long long *ref= (long long*) malloc((4*size+1)*sizeof(long long));/* ref[i] is index in Genome of each base dna[i] in current fragment */
  long long *refrev = (long long*)malloc((4*size+1)*sizeof(long long));
  if(!ref || !refrev){
    fprintf(stderr,"malloc of ref[] failed(%lu bytes)\n",(unsigned long)((4*size+1)*sizeof(long long)));
    exit(1);
  }

  for(j= 0; j < len; j++)
    dnafull[j] = acgt[(dnabin[j/4] >> (6-(j%4)*2)) & 3];
  dnafull[j] = '\0';

  if(spotsfile){
    FILE *fp = NULL;

    fprintf(stderr,"Writing reference map to %s\n",spotsfile);
    fflush(stderr);

    if((fp=fopen(spotsfile,"w"))==NULL){
      fprintf(stderr,"Cannot write reference to %s\n",spotsfile);
      exit(1);
    }
      
    /* display original command line as comment */
    fprintf(fp,"#\t");
    for(int i=0;i<argc0;i++)
      fprintf(fp,"%s ",argv0[i]);
    fprintf(fp,"\n");

    if(CMAP_REF){  /* output key header lines for CMAP_REF */
      fprintf(fp,"# CMAP_REF File Version:\t0.1\n");
      fprintf(fp,"# Label Channels:\t%d\n", numcolors);
    } else
      fprintf(fp,"#\tN Colors: %d\n", numcolors);
    for(int c = 0; c < numcolors;c++){
      char *p = uppercase(restrict[c][0]);
      if(CMAP_REF)
	fprintf(fp,"# Nickase Recognition Site %d:\t%s", c+1, p);
      else
	fprintf(fp,"#\tNickase%d: %s", c+1, p);
      (void)free(p);
      for(int t = 1; t < numenzymes[c]; t++){
	p = uppercase(restrict[c][t]);
	fprintf(fp,",%s", p);
	(void)free(p);
      }
      fprintf(fp,"\n");
    }

    /* display Genome Name */
    fprintf(fp,"#\tReference Name: %s\n",argv[0]);
    fprintf(fp,"#\tReference Size (Bp): %d\n",(int)(genome_length));

    if(CMAP_REF){
      fprintf(fp,"#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\n");
      fprintf(fp,"#f int\tfloat\tint\tint\tint\tfloat\tfloat\tfloat\tfloat\n");
    } else
      fprintf(fp,"NickID Color Location\n");

    numsites = 0;

    for(int c = 0; c < numcolors; c++){
      int numenz = numenzymes[c];
      int enzlen[MAXENZYMES];
      for(int e = 0; e < numenz; e++){
	enzlen[e] = strlen(invtrict[c][e])/2;
	found1[e] = found2[e] = dnafull;
      }
      for(j=0,thiscut=lastcut=0;thiscut!=len;lastcut=thiscut){
	thiscut = len;
	for(int e = 0; e < numenz; e++){
	  found1[e] = nextcut(found1[e],restrict[c][e]);
	  found2[e] = nextcut(found2[e],invtrict[c][e]);

	  thiscut1[e] = found1[e] ? (found1[e] - dnafull) + enzlen[e] : len;
	  thiscut2[e] = found2[e] ? (found2[e] - dnafull) + enzlen[e] : len;

	  thiscut = min(thiscut,min(thiscut1[e],thiscut2[e]));
	}

	/* skip over cut used */
	for(int e = 0; e < numenz; e++){
	  if (found1[e] && thiscut == thiscut1[e]) found1[e]++;
	  if (found2[e] && thiscut == thiscut2[e]) found2[e]++;
	}

	if(thiscut == len)
	  break;

	/* add thiscut to list of sites */
	if(numsites >= maxsites){/* reallocate sites */
	  int newmaxsites = maxsites*2;
	  Cspot *newsites = new Cspot[newmaxsites];
	  for(int i=0;i < numsites;i++)
	    newsites[i] = sites[i];
	  maxsites = newmaxsites;
	  delete [] sites;
	  sites = newsites;
	  delete [] Y;
	  Y = new double[maxsites];
	}
	
	Cspot *p = &sites[numsites++];
	p->color = c+1;
	p->site = thiscut;
      }
    }

    /*if(numcolors > 1)*/ /* sort all sites in order of site location */
    qsort(sites,numsites,sizeof(Cspot), (intcmp*) spotcmp);
    
    for(int i = 0; i < numsites; i++){
      if(CMAP_REF)
	fprintf(fp,"%lld\t%0.1f\t%d\t%d\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\n",
		(long long)(chr ? chr : 1), genome_length, numsites, i+1, sites[i].color, sites[i].site, 0.0, 1.0, 1.0);
      else
	fprintf(fp,"%d %lld %d\n",i,sites[i].color,(int)floor(sites[i].site+0.5));
    }

     /* add right end */
    if(CMAP_REF)
      fprintf(fp,"%lld\t%0.1f\t%d\t%d\t%lld\t%0.1f\t%0.1f\t%0.1f\t%0.1f\n",
	      (long long)(chr ? chr : 1), genome_length, numsites, numsites+1, 0LL, genome_length, 0.0, 1.0, 1.0);
    else
      fprintf(fp,"%d %lld %d\n",numsites+1, 0LL,(int)floor(genome_length+0.5));

    fflush(fp);
    (void)fclose(fp);

    fprintf(stderr,"Genome has %d sites\n", numsites);
  }

  fprintf (stderr, "Genome size = %0.3f (kb)\n", genome_length/1000.0);
  if(!(nmols > 0 || coverage > 0.0))
    exit(0);

  if(nmols > 0)
    fprintf (stderr, "Number of Maps -N%d\n", nmols);
  if(coverage > 0.0)
    fprintf (stderr, "Coverage of output -c%0.2f\n", coverage);
  fprintf (stderr, "mean Map size (kb) -M%.3f\n", mean/1000.0);
  if(stddev > 0.0)
    fprintf (stderr, "standard deviation of Map size (kb) -m%.3f\n", stddev/1000.0);
  else 
    fprintf (stderr, "minimum Map size (kb) -minlen%.3f\n", minlen/1000.0);
  for(int c= 0; c < numcolors; c++){
    fprintf (stderr, "Nickase%d = %s",c+1,restrict[c][0]);
    for(int e = 1; e < numenzymes[c]; e++)
      fprintf(stderr, ",%s", restrict[c][e]);
    fprintf(stderr," compliment = %s", invtrict[c][0]);
    for(int e = 1; e < numenzymes[c]; e++)
      fprintf(stderr,",%s", invtrict[c][e]);
    fprintf(stderr,"\n");
  }
  if(numcolors <= 1){
    fprintf (stderr, "False Postive density per 100 kb -FP%0.4f\n", FP[0]*100.0);
    fprintf (stderr, "False Negative fraction -FN%0.4f\n", FN[0]);
    fprintf (stderr, "sizing errors: -sf%0.6f -sd%0.6f -sr%0.6f -se%0.6f\n", 
	     sf[0],sd[0],sr[0],se[0]);
    if(res[0] > 0.0)
      fprintf (stderr, "Resolution in pixels : -res%0.4f  -resSD%0.4f\n", res[0], resSD[0]);
  } else {
    if(DEBUG) assert(numcolors==2);
    fprintf (stderr, "False Postive density per 100 kb -FP%0.4f_%0.4f\n", FP[0]*100.0,FP[1]*100.0);
    fprintf (stderr, "False Negative fraction -FN%0.4f_%0.4f\n", FN[0],FN[1]);
    fprintf (stderr, "sizing errors: -sf%0.6f_%0.6f -sd%0.6f_%0.6f -sr%0.6f_%0.6f -se%0.6f_%0.6f\n", 
	     sf[0],sf[1],sd[0],sd[1],sr[0],sr[1],se[0],sd[1]);
    if(res[0] > 0.0 || res[1] > 0.0)
      fprintf (stderr, "Resolution in pixels : -res%0.4f_%0.4f -resSD%0.4f_%0.4f\n", res[0], res[1], resSD[0], resSD[1]);
  }
  fprintf (stderr, "Rubber-banding error rate=%.3f\n",rubberband);

  fprintf (stderr, "random number seed -R%u\n", seed);
  fprintf (stderr, "chimerism probability -C%f\n", chimerism_probability);
  fprintf (stderr, "Pixel size = %0.3f bp\n", PixelLen);
  if(DetFragile)
    fprintf(stderr,"Deterministic Fragile sites every %0.3f kb (Fragile Level ranging from %0.1f%% to %0.1f%%)\n", DetFragile*0.001, DetFragileLevel1*100.0,(DetFragileLevel1+0.9*(DetFragileLevel2-DetFragileLevel1))*100.0);
  if(FragileMean > 0.0)
    fprintf(stderr,"Fragile sites based on Standard Normal Function with Mean= %0.1f bp and SD= %0.1f bp\n",FragileMean,FragileSD);

  if(!PIXELFRACTION)
    fprintf (stderr, "site locations rounded to nearest pixel\n");
  if(Stealth)
    fprintf (stderr, "true site locations will not be revealed\n");
  fprintf (stderr, "Generating sampled maps\n");

  if(coverage > 0.0){
    if(stddev > 0.0)
      nmols = coverage * genome_length / mean;
    else
      nmols = coverage * genome_length / (mean + minlen);
    fprintf (stderr, "Number of Maps -N%d\n", nmols);
  }
  if(nmols < 0){
    if(coverage > 0.0){
      long long nmolsLL = coverage * genome_length / (stddev > 0.0 ? mean : mean+minlen);
      printf("nmols=%d : possible 32 bit overflow : (true value = %lld)\n", nmols, nmolsLL);
    } else
      printf("nmols=%d : possible 32 bit overflow ???\n", nmols);
    exit(1);
  }
  fflush(stderr);

  if(nmols > 0 || coverage > 0.0){
    /* display original command line as comment */
    printf("#\t");
    for(i=0;i<argc0;i++)
      printf("%s ",argv0[i]);
    printf("\n");
    
    /* VFIXX file header */
    printf("#\tReference Name: %s\n",argv[0]);
    printf("#\tReference Size (Bp): %d\n", (int)floor(genome_length+0.5));
    printf("#\tN Colors:%d\n",numcolors);
    printf("#\tFragment Length (mean bp): %d\n",(int)floor(mean+0.5));
    if(stddev > 0.0)
      printf("#\tFragment Length standard deviation (bp): %d\n", (int)floor(stddev+0.5));
    else
      printf("#\tFragment Length Exponential distribution Lower Bound cutoff (bp): %d\n", (int)floor(minlen+0.5));
    for(int c= 0; c < numcolors; c++){
      char *p = uppercase(restrict[c][0]);
      printf("#\tNickase %d: %s",c+1, p);
      (void)free(p);
      for(int e = 1; e < numenzymes[c]; e++){
	char *p = uppercase(restrict[c][0]);
	printf(",%s",p);
	(void)free(p);
      }
      printf("\n");
    }
    printf("#\tN Fragments: %d\n", nmols);
    printf("#\tBases per Pixel: %0.6f\n",PixelLen);
    if(numcolors <= 1){
      printf("#\tFalse Positives/100kb: %0.4f\n", FP[0]*100.0);
      printf("#\tFalse Negative fraction: %0.4f\n", FN[0]);
      printf("#\tResolution (Px): %0.4f\n", res[0]);
      printf("#\tResolution SD (Px): %0.4f\n", resSD[0]);
      if(sf[0] > 0.0)
	printf("#\tScaling Variation sf(kb): %0.6f\n", sf[0]);
      if(sd[0] != 0.0)
	printf("#\tScaling Variation sd (kb^1/2): %0.6f\n", sd[0]);
      if(sr[0] > 0.0)
	printf("#\tScaling Variation sr (ratio): %0.6f\n", sr[0]);
      if(se[0] > 0.0)
	printf("#\tResolution Variation se (ratio): %0.6f\n", se[0]);
    } else {
      printf("#\tFalse Positives/100kb: %0.4f %0.4f\n", FP[0]*100.0,FP[1]*100.0);
      printf("#\tFalse Negative fraction: %0.4f %0.4f\n", FN[0],FN[1]);
      printf("#\tResolution (Px): %0.4f %0.4f\n", res[0],res[1]);
      printf("#\tResolution SD (Px): %0.4f %0.4f\n", resSD[0],resSD[0]);
      if(sf[0] > 0.0 || sf[1] > 0.0)
	printf("#\tScaling Variation sf(kb): %0.6f %0.6f\n", sf[0],sf[1]);
      if(sd[0] != 0.0 || sd[1] != 0.0)
	printf("#\tScaling Variation sd (kb^1/2): %0.6f %0.6f\n", sd[0],sd[1]);
      if(sr[0] > 0.0 || sr[1] > 0.0)
	printf("#\tScaling Variation sr (ratio): %0.6f %0.6f\n", sr[0], sr[1]);
      if(se[0] > 0.0 || se[1] > 0.0)
	printf("#\tResolution Variation se (ratio): %0.6f %0.6f\n", se[0], se[1]);
    }
    if(DetFragile)
      printf("#\tDeterministic Fragile sites every %0.3f kb (Fragile Level ranging from %0.1f%% to %0.1f%%)\n", DetFragile*0.001, DetFragileLevel1*100.0,(DetFragileLevel1+0.9*(DetFragileLevel2-DetFragileLevel1))*100.0);
    if(!PIXELFRACTION)
      printf("#\tPixel locations rounded to nearest integer\n");
    else
      printf("#\tPixel locations can be fractional\n");
  }
  fflush(stdout);

  for(int i = 0; i < numcolors; i++){
    sf[i] *= sf[i];
    sd[i] *= fabs(sd[i]);
    sr[i] *= sr[i];
    se[i] *= se[i];
  }

  /* initialize array of maps */
  int maxmols = (FragileMean > 0.0) ? 2*nmols : nmols;
  Cmap *map = new Cmap[maxmols];
  int *allnumsite = new int[maxmols*numcolors];
  Cspot **allspots = new Cspot*[maxmols*numcolors];
  for(i = 0; i < maxmols; i++){
    map[i].numsite = &allnumsite[i*numcolors];
    map[i].spots = &allspots[i*numcolors];
    map[i].chr = chr;
  }

  /* initialize 6 independent random number generator state arrays */
  /* so that theoretical contig will not vary with parameter settings */
  initstate (seed+5, (char *)state6, 128);
  initstate (seed+4, (char *)state5, 128);
  initstate (seed+3, (char *)state4, 128);
  initstate (seed+2, (char *)state3, 128);
  initstate (seed+1, (char *)state2, 128);
  initstate (seed, (char *)state1, 128);

  int sitecnt[MAXCOLORS], fpcnt[MAXCOLORS], tpcnt[MAXCOLORS];/* statistics of actual fn and fp rate in simulated data */
  double lensum = 0.0;
  for(int c = 0; c < numcolors; c++)
    sitecnt[c] = fpcnt[c] = tpcnt[c] = 0;

  if(DetFragile && DetFragile < 2 * mean){
    printf("ERROR: -D fragile site interval = %0.3f kb is smaller than 2x mean molecule size of %0.3f kb\n",
	   DetFragile*0.001, mean*0.001);
    exit(1);
  }

  double cumnormsq = 0.0;
  int normcnt = 0;

  int chim_cnt = 0;
  int nummols = nmols;/* fragile sites may break molecules and additional fragments are appended at map[nmols .. nummols-1] */
  for (int i = 0; i < nmols; i++) {
    long long totlen = 0;

    /* first pick the location(s) and length(s) of molecule i (2 locations/lengths for chimeric molecules), then copy the base pairs into dna[0..totlen-1] */

    setstate ((char *)state1); /* use independent random number generator */
    r = urandom ();
    if (r < chimerism_probability) {
      /* chimeric molecule */
      chim_cnt++;
      do {
	r = urandom();
	if(stddev > 0.0)
	  len = floor ((inv_cumul_normal (r) * stddev + mean));
	else
	  len = -mean * log(1.0 - r);
	r = urandom();
	if(LINEAR){ /* offset can range from -len to genome_length */
	  offset = floor ((genome_length + len) * r) - len;
	  /* if left end of sequence is before beginning of genome, chop it off */
	  if (offset < 0) {
	    len += offset; /* reduce the length */
	    offset = 0;
	  }
	  /* if right end of sequence is after end of genome, chop off right end */
	  if (len + offset > genome_length) {
	    len = genome_length - offset;
	  }
	} else /* offset can range from 0 to genome_length-len */
	  offset = floor((genome_length - len)*r);

	if(DEBUG) assert(0 <= offset && 0 <= len && offset+len <= genome_length);

	if(DetFragile){/* If region crosses fragile site, chop off the shorter piece */
	  double origoffset = offset;
	  double origlen = len;

	  /* NOTE : assumes DetFragile is larger than molecule size, so only one fragile size can be present */
	  long long DF = floor(DetFragile);
	  long long right = floor(len+offset);
	  long long left = floor(offset);
	  long long FragileSite = right - (right % DF);
	  if(FragileSite > left){
	    int Bin = (FragileSite % (DF*10))/DF;
	    double FragileLevel = DetFragileLevel1 + Bin*0.1*(DetFragileLevel2-DetFragileLevel1);
	    double r = urandom();
	    if(r < FragileLevel){
	      if(DEBUG) assert(FragileSite <= right);
	      if(FragileSite > (left + right)/2){/* chop off right of Fragile Site */
		len -= right - FragileSite;
	      } else {/* chop off left of Fragile Site */
		len -= FragileSite - left;
		offset = FragileSite;
	      }
	      if(DEBUG) assert(len <= origlen && offset >= origoffset);
	      if(DEBUG) assert(0 <= offset && 0 <= len && offset+len <= genome_length);
	    }
	  }
	}

	r = urandom();
	if(stddev > 0.0)
	  len2 = floor ((inv_cumul_normal (r) * stddev + mean));
	else
	  len2 = -mean * log(1.0-r);
	r = urandom();
	if(LINEAR){
	  /* offset can range from -len to genome_length */
	  offset2 = floor ((genome_length + len2) * r) - len2;
	  /* if left end of sequence is before beginning of genome, chop it off */
	  if (offset2 < 0) {
	    len2 += offset2; /* reduce the length */
	    offset2 = 0;
	  }
	  /* if right end of sequence is after end of genome, chop off right end */
	  if (len2 + offset2 > genome_length) {
	    len2 = genome_length - offset2;
	  }
	} else /* offset can range from 0 to genome_length-len */
	  offset2 = floor((genome_length - len2)*r);

	if(DEBUG) assert(0 <= offset2 && 0 <= len2 && offset2+len2 <= genome_length);

	if(DetFragile){/* If region crosses fragile site, chop off the shorter piece */
	  double origoffset2 = offset2;
	  double origlen2 = len2;

	  /* NOTE : assumes DetFragile is larger than molecule size, so only one fragile size can be present */
	  long long DF = floor(DetFragile);
	  long long right = floor(len2 + offset2);
	  long long left = floor(offset2);
	  long long FragileSite = right - (right % DF);
	  if(FragileSite > left){
	    int Bin = (FragileSite % (DF*10))/DF;
	    double FragileLevel = DetFragileLevel1 + Bin*0.1*(DetFragileLevel2-DetFragileLevel1);
	    double r = urandom();
	    if(r < FragileLevel){
	      if(DEBUG) assert(FragileSite <= right);
	      if(FragileSite > (left + right)/2){/* chop off right of Fragile Site */
		len2 -= right - FragileSite;
	      } else {/* chop off left of Fragile Site */
		len2 -= FragileSite - left;
		offset2 = FragileSite;
	      }
	      if(DEBUG) assert(len2 <= origlen2 && offset2 >= origoffset2);
	      if(DEBUG) assert(0 <= offset2 && 0 <= len2 && offset2+len2 <= genome_length);
	    }
	  }
	}
      } while (len <= 0.0 || len2 <= 0.0 || len + len2 < (stddev > 0.0 ? floor(mean/10.0)*2.0 : minlen));

      tlen = len;
      tlen2 = len2;
      toffset = offset;
      toffset2 = offset2;

      /* check allocation of dna[], rdna[], ref[] */
      ilen = (tlen+tlen2)/4LL + 2LL;
      if (ilen > size) {
	size = ilen;
	if(!(dna = (char *)realloc (dna, 4*size+1))){
	  fprintf(stderr,"realloc(dna,%lld) failed\n", 4*size+1);
	  exit(1);
	}
	if(!(ref = (long long *)realloc(ref, (4*size+1)*sizeof(long long)))){
	  fprintf(stderr,"realloc(ref,%lu) failed\n", (unsigned long)((4*size+1)*sizeof(long long)));
	  exit(1);
	}
	if(!(refrev = (long long *)realloc(refrev, (4*size+1)*sizeof(long long)))){
	  fprintf(stderr,"realloc(refrev,%lu) failed\n", (unsigned long)((4*size+1)*sizeof(long long)));
	  exit(1);
	}
      }

      /* copy the first of two molecules */
      for (j = 0; j < tlen; j++){
	dna[j] = dnafull[toffset + j];
	if(DEBUG/* HERE >=2 */) assert(dna[j]!=0);
	ref[j] = toffset + j;
      }
      totlen = j;

      map[i].startloc = toffset;
      map[i].endloc = toffset+tlen;
      map[i].len = tlen;

      map[i].startloc2 = toffset2;
      map[i].endloc2 = toffset2+tlen2;
      map[i].len2 = tlen2;

      map[i].chim = 1;

      if ((r=urandom()) < 0.5) {
	/* chimeric with 2nd part in same direction */
	/* save ">0" line */
	sprintf (namebuf, "%s_%lld_%lldC_%lld_%lld",
		 argv[0],toffset,tlen,toffset2,tlen2);
	map[i].name = strdup(namebuf);

	/* copy in forward direction */
	long long c = j;
	for (j = 0; j < tlen2; j++){
	  dna[c+j] = dnafull[toffset2 + j];
	  if(DEBUG/* HERE >=2 */) assert(dna[c+j]!=0);
	  ref[c+j] = toffset2 + j;
	}
	dna[totlen = c+j] = 0;
	if(DEBUG)assert(totlen == tlen+tlen2);
      } else {
	map[i].chimflip = 1;
	/* chimeric with reversed 2nd part */
	/* save ">0" line */
	sprintf(namebuf,"%s_%lld_%lldR_%lld_%lld",
		argv[0],toffset,tlen,toffset2,tlen2);
	map[i].name = strdup(namebuf);

	/* copy reverse-complement strand */
	long long c = tlen + tlen2 - 1;
	for (j = 0; j < tlen2; j++){
	  char tmp = dnafull[toffset2+j];
	  if(DEBUG>=2) assert(tmp != 0);
	  if(DEBUG>=2) assert(comp[(int)comp[(int)tmp]] == tmp);
	  dna[c-j] = comp[(int)tmp];
	  if(DEBUG/* HERE >=2 */) assert(dna[c-j] != 0);
	  ref[c-j] = toffset2 + j;
	}
	dna[totlen = tlen + tlen2] = 0;
      }

      tlen = len = len + len2;
    } else {/* not chimeric */
      do {
	r = urandom();
	if(stddev > 0.0)
	  len = floor (inv_cumul_normal (r) * stddev + mean);
	else /* HERE : sample using minimum r value of 1 - exp(-minlen/mean), or max 1-r value of exp(-minlen/mean), so len >= minlen */
	  len = -mean * log(1.0-r);
	r = urandom();
	if(LINEAR){ /* offset can range from -len to genome_length */
	  if(mean >= genome_length)
	    offset = 0;
	  else
	    offset = floor ((genome_length + len) * r) - len;
	  /* if left end of sequence is before beginning of genome, chop it off */
	  if (offset < 0) {
	    len += offset; /* reduce the length */
	    offset = 0;
	  }
	  /* if right end of sequence is after end of genome, chop off right end */
	  if (len + offset > genome_length) {
	    len = genome_length - offset;
	  }
	} else /* offset can range from 0 to genome_length - len */
	  offset = floor((genome_length - len)*r);

	if(DetFragile && len >= minlen){/* If region crosses fragile site, chop off the shorter piece */
	  if(DEBUG) assert(0 <= offset && 0 <= len && offset+len <= genome_length);
	  double origoffset = offset;
	  double origlen = len;

	  /* NOTE : assumes DetFragile is larger than molecule size, so only one fragile size can be present */
	  long long DF = floor(DetFragile);
	  long long right = floor(len+offset);
	  long long left = floor(offset);
	  long long FragileSite = right - (right % DF);
	  if(0){
	    fprintf(stderr,"offset=%0.1f,len=%0.1f:DF=%lld,left=%lld,right=%lld,FragileSite=%lld,chr=%d\n",
		    offset,len,DF,left,right,FragileSite,chr);
	    fflush(stderr);
	  }

	  if(FragileSite > left){
	    int Bin = (FragileSite % (DF*10))/DF;
	    double FragileLevel = DetFragileLevel1 + Bin*0.1*(DetFragileLevel2-DetFragileLevel1);
	    double r = urandom();
	    if(r < FragileLevel){
	      if(DEBUG) assert(FragileSite <= right);
	      if(FragileSite > (left + right)/2){/* chop off right of Fragile Site */
		len -= right - FragileSite;
	      } else {/* chop off left of Fragile Site */
		len -= FragileSite - left;
		offset = FragileSite;
	      }
	      if(0){
		fprintf(stderr,"     offset->%0.1f,len->%0.1f\n",offset,len);
		fflush(stderr);
	      }
	      if(DEBUG) assert(len <= origlen && offset >= origoffset);
	      if(DEBUG) assert(0 <= offset && 0 <= len && offset+len <= genome_length);
	    }
	  }
	}
      } while (len < (stddev > 0.0 ? floor(mean / 10.0) : minlen));

      /*      fprintf(stderr,"offset=%d,len=%d\n",(int)floor(offset+0.5),(int)floor(len+0.5));
	      fflush(stderr);*/

      tlen = len = floor(len);
      toffset = offset = floor(offset);
      if(tlen > MASK(31)){
	printf("Molecule length larger than %d not supported:tlen=%lld\n",MASK(31),tlen);
	exit(1);
      }

      /* save output ">0" information */
      map[i].startloc = toffset;
      map[i].endloc = toffset+tlen;
      map[i].len = tlen;
      map[i].startloc2 = map[i].endloc2 = map[i].len2 = 0;
      sprintf(namebuf,"%s-%lld-%lld",argv[0],toffset,tlen);
      map[i].name = strdup(namebuf);

      ioffset = toffset/4;
      ilen = ceil((toffset+tlen)/4.0)-ioffset;
      if(fseek (stdin, ioffset, SEEK_SET)){
	fprintf(stderr,"fseek of stdin failed\n");
	exit(1);
      }
      if (ilen > size) {
	size = ilen;
	if(!(dnabin = (char *)realloc (dnabin, size))){
	  fprintf(stderr,"realloc(dnabin,%lld) failed\n",size);
	  exit(1);
	}
	if(!(dna = (char *)realloc (dna, 4*size+1))){
	  fprintf(stderr,"realloc(dna,%lld) failed\n",4*size+1);
	  exit(1);
	}
#if 0
	if(!(rdna = (char *)realloc (rdna, 4*size+1))){
	  fprintf(stderr,"realloc(rdna,%lld) failed\n",4*size+1);
	  exit(1);
	}
#endif
	if(!(ref = (long long *)realloc(ref, (4*size+1)*sizeof(long long)))){
	  fprintf(stderr,"realloc(ref,%lu) failed\n",(unsigned long)((4*size+1)*sizeof(long long)));
	  exit(1);
	}
	if(!(refrev = (long long *)realloc(refrev, (4*size+1)*sizeof(long long)))){
	  fprintf(stderr,"realloc(refrev,%lu) failed\n", (unsigned long)((4*size+1)*sizeof(long long)));
	  exit(1);
	}
      }
      if(DEBUG) assert(ilen >= 0);
      k = toffset - ioffset*4;
      size_t t2;
      if(ilen != (t2 = fread (dnabin, 1, ilen, stdin))){
	fprintf(stderr,"fread of %lld bytes failed(t2=%lu)\n",ilen,t2);
	exit(1);
      }
      for (j = 0; j < tlen; j++){
	dna[j] = acgt[(dnabin[(j+k)/4] >> (6-((j+k)%4)*2)) & 3];
	if(DEBUG/* HERE >=2 */) assert(dna[j] != 0);
	ref[j] = offset + (j+k);
      }
      dna[totlen = j] = 0;
      if(DEBUG) assert(totlen == tlen);
    }/* not chimeric */

    /* randomly flip orientation of molecule */
    setstate ((char *)state6);
    if((r=urandom()) < 0.5){/* flip entire molecule and add "R" to end of name */
      map[i].flip = 1;
      sprintf(namebuf,"%sR",map[i].name);
      free(map[i].name);
      map[i].name = strdup(namebuf);
      if(DEBUG)assert(dna[totlen]==0);
      if(DEBUG)assert((size_t)totlen == strlen(dna));
      char *reversed = invcomp(dna);
      if(DEBUG)assert(reversed[totlen]==0);
      for(j = 0; j < totlen; j++){
	dna[j] = reversed[j];
	if(DEBUG>=2) assert(dna[j]!=0);
      }
      free(reversed);

      /* flip ref[0..totlen-1] */
      for(j = 0; j < totlen; j++)
	refrev[totlen-1-j] = ref[j];
      for(j = 0; j < totlen; j++)
	ref[j] = refrev[j];
    }

    if(VERB>=2){
      fprintf(stderr,"%s:startloc=%lld,endloc=%lld,tlen=%lld,id=%lld\n",namebuf,map[i].startloc,map[i].endloc,tlen,map[i].id());
      fflush(stderr);
    }

    lensum += len * 0.001;

    double maxlen = 0.0;

    if(rubberband){
      setstate ((char *)state5);/* generator for rubber_banding of current molecule (same for all colors) */
      rubberscale = /* WAS8 1.0/ */ (1.0+inv_cumul_normal(urandom())*rubberband);
    }

    for(int c= 0; c < numcolors; c++){
      int numenz = numenzymes[c];
      int enzlen[MAXENZYMES];
      for(int e = 0; e < numenz; e++){
	enzlen[e] = strlen(invtrict[c][e])/2;
	found1[e] = found2[e] = dna;
      }

      numsites = 0;

      /* generate locations of false cuts */
      setstate ((char *)state2); /* use independent random number generator for false positives */
      nfcuts = make_poisson_events ((int)len, FP[c], fcutloc);
      fcutloc[nfcuts] = len;

      fpcnt[c] += nfcuts;

      if(VERB>=2){
	fprintf(stderr,"i=%d/%d,c=%d:false cuts=%d:",i,nmols,c,nfcuts);
	for(int k= 0; k < nfcuts; k++)
	  fprintf(stderr,"  %0.3f",fcutloc[k]*0.001);
	fprintf(stderr," (kb)\n");
      }
	
      setstate ((char *)state3); /* use independent random number generator for false negative */

      double cumsite = 0.0;/* cumulative sum of intervals */
      thiscut = lastcut = fcut = 0;
      while (thiscut != tlen) {
	/* search for either orientation of restriction enzyme */
	long long mincut1, mincut2;
	mincut1 = mincut2 = len;
	for(int e = 0; e < numenz; e++){
	  found1[e] = nextcut(found1[e],restrict[c][e]);
	  found2[e] = nextcut(found2[e],invtrict[c][e]);

	  thiscut1[e] = found1[e] ? (found1[e] - dna) + enzlen[e] : len;
	  thiscut2[e] = found2[e] ? (found2[e] - dna) + enzlen[e] : len;

	  mincut1 = min(mincut1,thiscut1[e]);
	  mincut2 = min(mincut2,thiscut2[e]);
	}

	/* pick smallest of mincut1, mincut2 */
	thiscut = min(mincut1,mincut2);

	/* Skip over all cuts which are used this time. */
	for(int e = 0; e < numenz; e++){
	  if (found1[e] && thiscut == thiscut1[e]) found1[e]++;
	  if (found2[e] && thiscut == thiscut2[e]) found2[e]++;
	}

	if(thiscut != tlen)
	  sitecnt[c]++;

	/* decide on cutting efficiency but end of molecule or falsecut always cuts. */
	r = r1 = r2 = urandom();      
	if (r >= FN[c] || thiscut == tlen) {
	  if(thiscut != tlen)
	    tpcnt[c]++;

	  /* NOTE : add noise to intervals later after applying res/resSD, FP */

	  double f = thiscut - lastcut;
	  cumsite += f;

	  /* save current site in sites[numsites] */
	  if(numsites >= maxsites){/* reallocate sites */
	    int newmaxsites = maxsites*2;
	    Cspot *newsites = new Cspot[newmaxsites];
	    for(int j=0;j < numsites;j++)
	      newsites[j] = sites[j];
	    maxsites = newmaxsites;
	    delete [] sites;
	    sites = newsites;
	    delete [] Y;
	    Y = new double[maxsites];
	  }

	  Cspot *p = &sites[numsites];
	  p->color = c+1;
	  p->site = cumsite;
	  p->refsite = ref[thiscut];
	  lastcut = thiscut;	    
	  if(thiscut != tlen)
	    numsites++;
	}
      
	/* if mincut1 and mincut2 are close enough introduce a fragile site break : save as special site (color <= 0) and break molecule later on */
	if(FragileMean > 0.0 && max(mincut1,mincut2) < tlen){
	  double delta = fabs((double)(mincut1 - mincut2));
	  if(fabs(delta) <= FragileMean + 10.0 * FragileSD){
	    long long fsite = 0;
	    if(FragileSD > 0.0){
	      double InvFragileSD = 1.0/FragileSD;
	      double IfragileSD = InvFragileSD * (1.0/sqrt(2.0));
	      double LogDelta = (1.0 + log(delta * InvFragileSD)) * FragileSD;/* log-transform of delta */
	      double p = Pr(LogDelta, FragileMean, IfragileSD);
	      double r = urandom();
	      if(r > p)
		fsite = (mincut1 + mincut2)/2;
	    } else {
	      if(delta < FragileMean)
		fsite = (mincut1 + mincut2)/2;
	    }	    
	    if(fsite){
	      if(VERB>=2){
		fprintf(stderr,"%s:map[%d],c=%d:Found fragile site break at fsite=%lld,ref=%lld (mincut1=%lld,mincut2=%lld,thiscut=%lld,tlen=%lld):FragileMean=%0.1f,FragileSD=%0.1f\n",
			namebuf,i,c,fsite,ref[fsite],mincut1,mincut2,thiscut,tlen,FragileMean,FragileSD);
		fflush(stdout);
	      }
	      if(DEBUG) assert(fsite >= thiscut);
	      if(DEBUG) assert(fsite < tlen);

	      /* save fragile site midpoint in sites[numsites] */
	      if(numsites >= maxsites){/* reallocate sites */
		int newmaxsites = maxsites*2;
		Cspot *newsites = new Cspot[newmaxsites];
		for(int j=0;j < numsites;j++)
		  newsites[j] = sites[j];
		maxsites = newmaxsites;
		delete [] sites;
		sites = newsites;
		delete [] Y;
		Y = new double[maxsites];
	      }

	      Cspot *p = &sites[numsites++];
	      p->color = -max(mincut1,mincut2);// fragile site
	      p->site = fsite;
	      p->refsite = ref[fsite];
	    }
	  }
	}
      }

      if(VERB>=2){
	fprintf(stderr,"%s:Before FP,res,var:c=%d,numsites=%d,len=%0.4f kb\n",namebuf,c+1,numsites,sites[numsites].site*0.001);
	for(int I = 0; I <= numsites; I++)
	  fprintf(stderr,"    site[%d]=%0.4f kb,c=%lld,ref=%0.4f kb\n",I,sites[I].site*0.001, sites[I].color, sites[I].refsite*0.001);
	fflush(stderr);
      }
      if(DEBUG) assert(fabs(sites[numsites].site - tlen) < 1.0);

      /* add in the false positive cuts */
      double finallen = sites[numsites].site;
      for(int fc = 0; fc < nfcuts; fc++){
	cumsite = fcutloc[fc] * (finallen / len);
	if(DEBUG) assert(cumsite < finallen);

	/* save current site in sites[numsites] */
	if(numsites >= maxsites-1){/* reallocate sites */
	  int newmaxsites = maxsites*2 + 1;
	  Cspot *newsites = new Cspot[newmaxsites];
	  for(int i=0;i < numsites;i++)
	    newsites[i] = sites[i];
	  maxsites = newmaxsites;
	  delete [] sites;
	  sites = newsites;
	  delete [] Y;
	  Y = new double[maxsites];
	}

	Cspot *p = &sites[++numsites];
	p->color = c+1;
	p->site = cumsite;
	p->refsite = -1;
      }

      /* sort sites in ascending order */
      qsort(sites,numsites+1,sizeof(Cspot), (intcmp*) spotcmp);

      if(VERB>=2){
	fprintf(stderr,"%s:Before res,var:c=%d,numsites=%d,len=%0.4f kb\n",namebuf,c+1,numsites,sites[numsites].site*0.001);
	for(int I = 0; I <= numsites; I++)
	  fprintf(stderr,"    site[%d]=%0.4f kb,ref=%0.4f kb\n",I,sites[I].site*0.001, sites[I].refsite*0.001);
	fflush(stderr);
      }
      if(DEBUG) assert(fabs(sites[numsites].site - tlen) < 1.0);

      /* apply res,resSD to sites[0..numsites-1] */
      for(int I = 0; I <= numsites; I++){
	Y[I] = sites[I].site;
	sites[I].deres = 0.0;
	sites[I].refsiteL = sites[I].refsite;
      }
      double resKB = res[c] * PixelLen;
      double IresSD = (resSD[c] > 0.0) ? 1.0/(resSD[c]*PixelLen * sqrt(2.0)) : 1.0;
      for(int I = 0; I < numsites-1; I++){
	if(sites[I].color <= 0)// fragile site
	  continue;
	int J,K;
	for(J = I+1; J < numsites; J++){
	  if(sites[J].color <= 0)// fragile site
	    break;
	  if(resSD[c] > 0.0){
	    double r = urandom();
	    double p = Pr(Y[J] - Y[J-1],resKB,IresSD);
	    if(r < p)
	      break;
	  } else if(Y[J] - Y[J-1] > resKB)
	    break;
	}
	if(--J > I){/* merge sites I..J */
	  sites[I].deres = Y[J] - Y[I];
	  sites[I].site = Y[I] = Yc(Y,J,J-I);

	  /* scan refsite values (excluding false positive sites) and record range as refsiteL .. refsite */
	  long long refsite = -1, refsiteL = -1;
	  int k = I;
	  for(; k <= J; k++){
	    if(sites[k].refsite < 0)
	      continue;
	    refsite = refsiteL = sites[k].refsite;
	    break;
	  }
	  for(k++; k <= J; k++){
	    if(sites[k].refsite < 0)
	      continue;
	    refsite = sites[k].refsite;
	  }
	  sites[I].refsite = refsite;
	  sites[I].refsiteL = refsiteL;
	}
	/* shift remaining sites */
	for(K= I+1; ++J <= numsites;K++){
	  Y[K] = Y[J];
	  sites[K] = sites[J];
	}
	numsites = K-1;
      }

      if(VERB>=2){
	fprintf(stderr,"%s:Before var:c=%d,numsites=%d,len=%0.4f kb\n",namebuf,c+1,numsites,sites[numsites].site*0.001);
	for(int I = 0; I <= numsites; I++)
	  fprintf(stderr,"    site[%d]=%0.4f kb,ref=%0.4f,%0.4f kb, deres=%0.4f kb\n",I,sites[I].site*0.001, sites[I].refsiteL*0.001, sites[I].refsite*0.001, sites[I].deres * 0.001);
	fflush(stderr);
      }
      if(DEBUG) assert(fabs(sites[numsites].site - tlen) < 1.0);

      /* add sizing error to all intervals Y[0..numsites] */
      Y[0] = sites[0].site;
      for(int j = 1; j <= numsites;j++)
	Y[j] = sites[j].site - sites[j-1].site;

      double lastsite = 0.0;
      int lastj = -1;
      for(int j = 0; j <= numsites;j++){
	if(j < numsites && (sites[j].refsite < 0 || sites[j].color <= 0))
	  continue;/* skip false positive or fragile sites : they will be rescaled linearly by interpolation */

	/* process intervals lastsite .. sites[j].site by rescaling Y[lastj+1 .. j] */
	double x = (sites[j].site - lastsite) * 0.001;/* interval size in kb */
	if(DEBUG && !(x > 0.0)){
	  fprintf(stderr,"i=%d/%d,nummols=%d,c=%d:j=%d/%d:lastsite=%0.3f bp,sites[j].site=%0.3f bp, x=%6f kb\n",i,nmols,nummols,c,j,numsites,lastsite,sites[j].site,x);
	  for(int I = 0; I <= numsites; I++)
	    fprintf(stderr,"    site[I=%d]=%0.4f kb,ref=%0.4f,%0.4f kb, deres=%0.4f kb,color=%lld:Y[I]=%0.4f kb\n", 
		    I, sites[I].site*0.001, sites[I].refsiteL*0.001, sites[I].refsite*0.001, sites[I].deres * 0.001, sites[I].color, Y[I] * 0.001);

	  fflush(stderr);
	  assert (x > 0.0);
	}
	double var = (SF_FIX ? 0.0 : sf[c]) + x * (sd[c] + x * sr[c]);
	if(se[c] > 0.0){
	  double resvar = 0.0;
	  if(j > 0) 
	    resvar += sites[j-1].deres * sites[j-1].deres;
	  if(j < numsites)
	    resvar += sites[j].deres * sites[j].deres;
	  var += resvar * se[c] * 1e-6;
	}
	if(DEBUG) assert(var >= 0.0);
	noise_stddev = sqrt(var);

	double delta = noise_stddev*noise_stddev/x;
	double plambda = x/delta;
	if(plambda  < 0.25 /* WAS 1.0 */){/* uniform distribution from -x .. + x */
	  r = r3 = urandom();
	  noise = (2.0*r-1.0)*x;
	} else if(0 && plambda < 16.0){
	  /* use poisson noise with linear interpolation */
	  r = r3 = random_poisson(plambda-0.5)+urandom();
	  noise = r*delta-x;
	} else {/* use gaussian noise truncated at -x and +x*/
	  double gauss = 0.0;
	  if(noise_stddev > 0.0){
	    do { 
	      r = r3 = urandom();
	      gauss = inv_cumul_normal(r);
	      noise = gauss * noise_stddev;
	    } while (noise <= -x || noise >= x);
	  } else
	    noise = 0.0;
	  
	  if(VERB>=2 && x > 0.001) {
	    double norm = noise / noise_stddev;
	    cumnormsq += norm *= norm;
	    normcnt++;
	    fprintf(stderr,"    x=%0.6f kb, y=%0.6f kb, var=%0.8f, r=%0.8f, gauss= %0.8f,norm= %0.8f, err=%0.6f(cumnorm=%0.8f/%d=%0.8f\n",
		    max(0.0,x + noise), x, var, r , gauss, norm, noise, cumnormsq,normcnt,cumnormsq/normcnt);
	  }
	}

	double y = x + noise;
	if(y < 0.0)
	  y = 0.0;
	y /= x;

	for(int J = lastj; ++J <= j;)
	  Y[J] *= y;

	lastsite = sites[j].site;
	lastj = j;
      }

      /* accumulate interval sizes Y[0..numsites] */
      double x = 0.0;
      for(int j = 0; j <= numsites; j++)
	sites[j].site = x += Y[j];

      if(SF_FIX && sf[c] > 0.0){ /* apply sf error as local site peturbation */
	if(VERB>=2){
	  fprintf(stderr,"%s:Before sf:c=%d,numsites=%d,len=%0.4f kb\n",namebuf,c+1,numsites,sites[numsites].site*0.001);
	  for(int I = 0; I <= numsites; I++)
	    fprintf(stderr,"    site[%d]=%0.4f kb,ref=%0.4f,%0.4f kb\n",I,sites[I].site*0.001, sites[I].refsiteL*0.001,sites[I].refsite*0.001);
	  fflush(stderr);
	}

	double var = sf[c]*0.5;
	noise_stddev = sqrt(var) * 1000.0;
	double minsite = sites[numsites].site;
	double maxsite = 0.0;
	for(int j = 0; j < numsites; j++){
	  double r = urandom();
	  double gauss = inv_cumul_normal(r);
	  double noise = gauss * noise_stddev;
	  sites[j].site += noise;
	  minsite = min(minsite,sites[j].site);
	  maxsite = max(maxsite,sites[j].site);
	}
	sites[maxsites].site = max(maxsite,sites[maxsites].site);

	if(minsite < 0.0){/* shift all sites right by this amount */
	  double shift = -minsite;
	  for(int j = 0; j <= numsites; j++)
	    sites[j].site += shift;
	}
	/* sort sites in ascending order */
	qsort(sites,numsites+1,sizeof(Cspot), (intcmp*) spotcmp);
      }

      if(rubberscale != 1.0){
	if(VERB>=2){
	  fprintf(stderr,"%s:Before rubberscale=%0.6f:c=%d,numsites=%d,len=%0.4f kb\n",namebuf,rubberscale,c+1,numsites,sites[numsites].site*0.001);
	  for(int I = 0; I <= numsites; I++)
	    fprintf(stderr,"    site[%d]=%0.4f kb,ref=%0.4f,%0.4f kb\n",I,sites[I].site*0.001, sites[I].refsiteL*0.001,sites[I].refsite*0.001);
	  fflush(stderr);
	}
	for(int j = 0; j <= numsites; j++){
	  sites[j].site *= rubberscale;
	  if(DEBUG && !((j > 0 ? sites[j-1].site : 0.0) <= sites[j].site)){
	    fprintf(stderr,"j=%d:sites[j].site=%0.3f bp, sites[j-1].site=%0.3f bp\n",
		   j,sites[j].site,(j > 0) ? sites[j-1].site : 0.0);
	    fflush(stderr);
	    assert((j > 0 ? sites[j-1].site : 0.0) <= sites[j].site);
	  }
	}
      }

      maxlen = max(sites[numsites].site,maxlen);
      
      if(VERB>=2){
	fprintf(stderr,"i=%d,c=%d:numsites=%d,len=%0.4f kb,maxlen=%0.4f kb,PixelLen=%0.1f bp\n",i,c,numsites,sites[numsites].site*0.001,maxlen*0.001,PixelLen);
	for(int I = 0; I <= numsites; I++)
	  fprintf(stderr,"    site[%d]=%0.4f kb,ref=%0.4f,%0.4f kb\n",I,sites[I].site*0.001, sites[I].refsiteL*0.001,sites[I].refsite*0.001);
	fflush(stderr);
      }

      /* save Reference and Data map for current color c */
      map[i].numsite[c] = numsites;
      map[i].spots[c] = new Cspot[numsites+1];
      for(int j= 0; j < numsites; j++)
	map[i].spots[c][j] = sites[j];
    }// for c = 0 .. numcolors-1
    
    for(int c = 0; c < numcolors; c++){
      map[i].spots[c][map[i].numsite[c]].site = maxlen;
      map[i].spots[c][map[i].numsite[c]].color = 0;
      map[i].spots[c][map[i].numsite[c]].refsiteL = -2;
      map[i].spots[c][map[i].numsite[c]].refsite = -2;
    }

    if(FragileMean > 0.0){   // break molecules at fragile sites with color <= 0 and append to map[nummols .. maxmols-1]
      double frsite[256];/* list of fragile site breaks */
      int frcnt = 0;
      for(int c = 0; c < numcolors; c++)
	for(int j = 0; j < map[i].numsite[c]; j++)
	  if(map[i].spots[c][j].color <= 0 && frcnt < 255)
	    frsite[frcnt++] = map[i].spots[c][j].site;

      if(frcnt > 0){
	/* sort fragile sites in ascending order */
	qsort(frsite,frcnt,sizeof(double),(intcmp*) doubleinc);
	frsite[frcnt] = maxlen;
      
	if(VERB>=2){
	  fprintf(stderr,"map[%d]:found %d fragile site breaks:\n",i,frcnt);
	  for(int c = 0; c < numcolors; c++)
	    for(int j = 0; j <= map[i].numsite[c]; j++)
	      fprintf(stderr,"   site[%d][%d]= %0.4f kb, ref=%0.4f,%0.4f kb, color=%lld\n",
		      c,j,map[i].spots[c][j].site*0.001, map[i].spots[c][j].refsiteL*0.001, map[i].spots[c][j].refsite*0.001, map[i].spots[c][j].color);
	  fflush(stderr);
	}
	
	for(int k = 0; k < frcnt; k++){/* split off molecule from frsite[k] to frsite[k+1] */
	  if(nummols >= maxmols)
	    break;
	  double left = frsite[k], right = frsite[k+1];
	  if(right - left < minlen)
	    continue;
	  for(int c = 0; c < numcolors; c++){
	    int j = 0;
	    for(; j < map[i].numsite[c];  j++)
	      if(map[i].spots[c][j].site >= left)
		break;
	    int jleft = j;
	    int numsites = 0;
	    for(; j < map[i].numsite[c]; j++){
	      double site = map[i].spots[c][j].site;
	      if(site > right)
		break;
	      if(map[i].spots[c][j].color > 0)
		numsites++;
	    }
	    int jright = j;
	    map[nummols].numsite[c] = numsites;
	    map[nummols].spots[c] = new Cspot[numsites+1];
	    int k = 0;
	    for(j = jleft; j < jright; j++){
	      if(map[i].spots[c][j].color <= 0)
		continue;
	      if(DEBUG) assert(k < numsites);
	      map[nummols].spots[c][k] = map[i].spots[c][j];
	      map[nummols].spots[c][k].site -= left;
	      k++;
	    }
	    if(DEBUG) assert(k == numsites);
	    Cspot *p = &map[nummols].spots[c][numsites];
	    p->site = right - left;
	    p->color = 0;
	    p->refsiteL = p->refsite = -2;
	  }
	  
	  if(VERB>=2){
	    fprintf(stderr,"map[%d]: created fragile site fragment from %0.4f .. %0.4f kb as map[%d]:\n",i,left,right,nummols);
	    for(int c = 0; c < numcolors; c++)
	      for(int j = 0; j <= map[nummols].numsite[c]; j++)
		fprintf(stderr,"   site[%d][%d]= %0.4f kb, ref=%0.4f,%0.4f kb, color=%lld\n",
			c,j,map[nummols].spots[c][j].site*0.001, map[nummols].spots[c][j].refsiteL*0.001, map[nummols].spots[c][j].refsite*0.001, map[nummols].spots[c][j].color);
	    fflush(stderr);
	  }

	  nummols++;
	}

	/* truncate map[i] at length = frsite[0] */
	double left = frsite[0];
	if(left < minlen){/* discard map[i] by setting len = len2 = 0 (will be filtered out later) */
	  if(VERB>=2){
	    fprintf(stderr,"Discarding leftmost fragile fragment of map[%d] : length=%0.4f kb\n",i,left * 0.001);
	    fflush(stderr);
	  }
	  map[i].len = map[i].len2 = 0.0;
	} else {
	  double origlen = map[i].spots[0][map[i].numsite[0]].site;  
	  for(int c = 0; c < numcolors; c++)
	    for(int j = 0; j < map[i].numsite[c]; j++)
	      if(map[i].spots[c][j].site >= left){
		Cspot *p = &map[i].spots[c][j];
		p->site = left;
		p->color = 0;
		p->refsite = -2;
		p->refsiteL = -2;
		map[i].numsite[c] = j-1;
		break;
	      }
	  if(VERB>=2){
	    fprintf(stderr,"Truncated map[%d] from length= %0.4f to %0.4f kb\n",i, origlen * 0.001, left * 0.001);
	    for(int c = 0; c < numcolors; c++)
	      for(int j = 0; j <= map[i].numsite[c]; j++)
		fprintf(stderr,"   site[%d][%d]= %0.4f..%0.4f kb, ref=%0.4f kb, color=%lld\n",c,j,map[i].spots[c][j].site*0.001, 
			map[i].spots[c][j].refsiteL*0.0001, map[i].spots[c][j].refsite*0.001, map[i].spots[c][j].color);
	    fflush(stderr);
	  }
	}

	if(VERB>=2){
	  fprintf(stderr,"Debug exit\n");
	  fflush(stderr);
	  exit(1);
	}
      }
    }
  } // i = 0 .. nmols -1

  (void)fclose(fp);/* same as stdin */

  /* output remainder of VFIXX header information with actual error statistics */
  if(nmols > 0){
    for(register int c= 0; c < numcolors; c++){
      printf("#\tObserved False Negative Fraction for color %d: %0.6f\n", 
	     c+1, 1.0 - ((double)tpcnt[c])/sitecnt[c]);
      printf("#\tObserved False Positive Density for color %d (per 100kb): %0.5f\n", 
	     c+1, (fpcnt[c]*100.0/lensum));
    }
    printf("#\tAssigned Chimerism rate: %0.4f\n", chimerism_probability);
    printf("#\tObserved Chimerism rate: %0.4f\n", ((double)chim_cnt)/nmols);
    printf("#\t>0\tMolID\tStart(Bp)\tEnd(Bp)\tLength(Bp)\tStart2(Bp)\tEnd2(Bp)\tLength2(Bp)\tFlipped\tChimflip\tRefString\n");
    for(register int c=0; c < numcolors; c++){
      printf("#\t>%d\tNick%d Locations(Bp) ... \n",c+1,c+1);
      printf("#\tFragLen(Px)\tNick%d Locations(Px) ...\n",c+1);
    }
  }

  /* output all map information */
  Cmap *p = map;
  for(int i = 0; i < nummols; i++, p++){
    if(p->len + p->len2 < minlen)
      continue;

    /* output ">0" line */
    if(Stealth)
      printf (">0\t%lld\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d%s\n",
	      CHR_OFFSET*chr + i + 1LL, 0, 0, 0, 0, 0, 0, 0, 0, p->name);
    else
      printf (">0\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%d\t%d\t%s\n",
	      p->id(), chr*CHR_OFFSET + p->startloc, chr*CHR_OFFSET + p->endloc, p->len, 
	      p->startloc2 ? chr*CHR_OFFSET + p->startloc2 : 0, p->endloc2 ? chr*CHR_OFFSET + p->endloc2 : 0, p->len2, (int)p->flip, p->chim ? (int)p->chimflip : 0, p->name);

    for(int c= 0; c < numcolors; c++){
      int numsites = p->numsite[c];

      /* output Reference and Data map for current color c */
      printf(">%d",c+1);
      for(int j= 0; j < numsites; j++)
	printf("\t%lld\t%lld", Stealth ? 0 : p->spots[c][j].refsiteL < 0 ? -1 : p->spots[c][j].refsiteL + chr*CHR_OFFSET, 
	       Stealth ? 0 : p->spots[c][j].refsite < 0 ? -1 : p->spots[c][j].refsite + chr*CHR_OFFSET);
      printf("\n");
      
      if(PIXELFRACTION)
	printf("%0.6f", p->spots[c][numsites].site/PixelLen);
      else
	printf("%d", (int)floor((p->spots[c][numsites].site/PixelLen)+0.5));
      for(int j= 0; j < numsites; j++)
	if(PIXELFRACTION)
	  printf("\t%0.6f", p->spots[c][j].site/PixelLen);
	else
	  printf("\t%d", (int)floor((p->spots[c][j].site/PixelLen)+0.5));
      printf("\n");

      delete [] p->spots[c];
    }
    (void)free(p->name);
  }

  if(ref){
    (void)free(ref);
    ref = 0;
  }
  if(refrev){
    (void)free(refrev);
    refrev = 0;
  }
  if(dnabin){
    free(dnabin);
    dnabin = 0;
  }
  if(dna){
    free(dna);
    dna = 0;
  }
#if 0
  if(rdna){
    free(rdna);
    rdna = 0;
  }
#endif
  if(spotsfile){
    (void)free(spotsfile);
    spotsfile = 0;
  }
  delete [] allnumsite;
  delete [] allspots;
  delete [] map;
  delete [] sites;
  delete [] Y;
  for(int c = 0; c < numcolors; c++){
    for(int e = 0; e < numenzymes[c]; e++){
      (void)free(restrict[c][e]);
      (void)free(invtrict[c][e]);
    }
  }
  return 0;
}

