#ifndef VEC_EXP_H
#define VEC_EXP_H

static Ident vec_exp_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/vec_exp.h 5675 2016-12-07 02:48:17Z tanantharaman $");

#include <math.h>

//#ifndef FP_FAST_FMAF
//#define fmaf(x,y,z) (x*y + z)
//#endif

/* vectorizable single precision exp and erfc functions */

/* Compute exponential base e. Maximum ulp error = 0.87161 */

#define MAXARG 88.722839f
#define EPOWF_SHIFT 23
#define EPOWF_IX 127

// #pragma omp declare simd
static inline float vec_expf (float a)
{
    float c, f, r;
    int i;

    // exp(a) = 2**i * exp(f) where i = rint (a / log(2)); f = a - i * log(2)
    c = 0x1.800000p+23f; // 1.25829120e+7
    r = fmaf (0x1.715476p+0f, a, c) - c; // 1.44269502e+0 = 1.0/log(2)
    f = fmaf (r, -0x1.62e400p-01f, a); // -6.93145752e-1 // log_2_hi 
    f = fmaf (r, -0x1.7f7d1cp-20f, f); // -1.42860677e-6 // log_2_lo
    i = (int)r;

    // approximate r = exp(f) on interval [-log(2)/2,+log(2)/2]
    r =             0x1.6a98dap-10f;  // 1.38319808e-3
    r = fmaf (r, f, 0x1.1272cap-07f); // 8.37550033e-3
    r = fmaf (r, f, 0x1.555a20p-05f); // 4.16689515e-2
    r = fmaf (r, f, 0x1.55542ep-03f); // 1.66664466e-1
    r = fmaf (r, f, 0x1.fffff6p-02f); // 4.99999851e-1
    r = fmaf (r, f, 0x1.000000p+00f); // 1.00000000e+0
    r = fmaf (r, f, 0x1.000000p+00f); // 1.00000000e+0

    // exp(a) = 2**i * exp(f);
#if 0 // original code : more portable, may not vectorize for gcc
    r = ldexpf (r, i); 
#else
    const int ie = ((i + EPOWF_IX) << EPOWF_SHIFT);
    const float &de = reinterpret_cast<const float &>(ie);
    r = r * de;
#endif

    // handle special cases
#if 0 // original code
    if (!(fabsf (a) < 104.0f)) {
        r = a + a; // handle NaNs
        if (a < 0.0f) r = 0.0f;
        if (a > 0.0f) r = 1e38f * 1e38f; // + INF
    }
#else
    // handle only overflow,underflow
    r = (i < -EPOWF_IX) ? 0.0f : (a > MAXARG) ? (1e38f * 1e38f) : r;
#endif

    if(DEBUG>=2 && !isfinite(r) && !(r > 1e38f)){
      printf("vec_exp(%0.8f): i= %d, f= %0.8f, r = %0.8e (MAXARG= %0.8f)\n",a,i,f,r,MAXARG);
      fflush(stdout);
    }

    return r;
}

/*  
 * Based on: M. M. Shepherd and J. G. Laframboise, "Chebyshev Approximation of 
 * (1+2x)exp(x^2)erfc x in 0 <= x < INF", Mathematics of Computation, Vol. 36,
 * No. 153, January 1981, pp. 249-253.  
 */  
// #pragma omp declare simd
static inline float vec_erfcf (float x)
{
    float a, d, e, m, p, q, r, s, t;

    a = fabsf (x); 

    /* Compute q = (a-2)/(a+2) accurately. [0, 10.0546875] -> [-1, 0.66818] */
    m = a - 2.0f;
    p = a + 2.0f;
    r = 1.0f / p;
    q = m * r;
    t = fmaf (q + 1.0f, -2.0f, a); 
    e = fmaf (q, -a, t); 
    q = fmaf (r, e, q); 

    /* Approximate (1+2*a)*exp(a*a)*erfc(a) as p(q)+1 for q in [-1, 0.66818] */
    p =             -0x1.a48024p-12f;  // -4.01020574e-4
    p = fmaf (p, q, -0x1.42a172p-10f); // -1.23073824e-3
    p = fmaf (p, q,  0x1.585784p-10f); //  1.31355994e-3
    p = fmaf (p, q,  0x1.1ade24p-07f); //  8.63243826e-3
    p = fmaf (p, q, -0x1.081b72p-07f); // -8.05991236e-3
    p = fmaf (p, q, -0x1.bc0b94p-05f); // -5.42047396e-2
    p = fmaf (p, q,  0x1.4ffc40p-03f); //  1.64055347e-1
    p = fmaf (p, q, -0x1.540840p-03f); // -1.66031361e-1
    p = fmaf (p, q, -0x1.7bf612p-04f); // -9.27639678e-2
    p = fmaf (p, q,  0x1.1ba03ap-02f); //  2.76978403e-1

    /* Divide (1+p) by (1+2*a) ==> exp(a*a)*erfc(a) */
    t = a + a;
    d = t + 1.0f;
    r = 1.0f / d;
    q = fmaf (p, r, r); // q = (p+1)/(1+2*a)
    e = (p - q) + fmaf (q, -t, 1.0f); // (p+1) - q*(1+2*a)
    r = fmaf (e, r, q);

    /* Multiply by exp(-a*a) ==> erfc(a) */
    s = a * a; 
    e = vec_expf (-s);  
    t = fmaf (a, -a, s);
    r = fmaf (r, e, r * e * t);

    /* Handle NaN arguments to erfc() */
    if (!(a <= 0x1.fffffep127f)) r = x + x;

    /* Clamp result for large arguments */
    if (a > 10.0546875f) r = 0.0f;

    /* Handle negative arguments to erfc() */
    if (x < 0.0f) r = 2.0f - r; 

    return r;
}

#endif
