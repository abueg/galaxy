#ifndef EPOW_H
#define EPOW_H

static Ident epow_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/epow.h 6184 2017-03-28 06:50:46Z tanantharaman $");

#define EPOWD_SHIFT 52
#define EPOWD_ADD 1023l
#define EPOWD_IXMAX 1022l
#define EPOWD_IXMIN -1021l

#define EPOWD_LOG2E 1.442695040888963407359924681001892137426645954152985934135

#define EPOWD_C1 0.693147180559945309417232121458176568075500134360255254120
#define EPOWD_C2 0.240226506959100712333551263163332485865276475797272793433
#define EPOWD_C3 0.055504108664821579953142263768621757359354414224767775398
#define EPOWD_C4 0.009618129107628477161979071573658865479929845523272117657
#define EPOWD_C5 0.001333355814642844342341222198799617473078225125878224529
#define EPOWD_C6 0.000154035303933816099544370973327423479615309708979558496
#define EPOWD_C7 0.000015252733804059840280025439012009638169614541237242699
#define EPOWD_C8 1.3215486790144309488403758228288360755748282365107953e-6
#define EPOWD_C9 1.0178086009239699727490007597744629254947068155642169e-7
#define EPOWD_C10 7.0549116208011233298753921815507382589380121206357913e-9
#define EPOWD_C11 4.445538271870811497596408558888112052997018680469243e-10

#pragma omp declare simd
static inline double epow(const double x)
{
  const double x2 = x*EPOWD_LOG2E;

  const double offset = (x < 0.0) ? -0.5 : 0.5;
  const long ix = long(x2 + offset);
  const double d = x2 - double(ix);

  const long ixc = ix < EPOWD_IXMAX ? ix : EPOWD_IXMAX;
//  const long ie = ((ixc+EPOWD_ADD) << EPOWD_SHIFT);
//  const double &de = reinterpret_cast<const double &>(ie);
  const union { long ie; double de;} u = { ((ixc+EPOWD_ADD) << EPOWD_SHIFT) };
  const double dec = (ixc < EPOWD_IXMIN) ? 0.0 : u.de;
  const double d2 = 1.0+d*(EPOWD_C1+d*(EPOWD_C2+d*(EPOWD_C3+d*(EPOWD_C4+d*(EPOWD_C5+d*(EPOWD_C6+d*(EPOWD_C7+d*(EPOWD_C8+d*(EPOWD_C9+d*(EPOWD_C10+d*EPOWD_C11))))))))));

  return d2 * dec;
}

#define EPOWF_SHIFT 23
#define EPOWF_IX 127
#define EPOWF_IXMAX 126
#define EPOWF_IXMIN -127 // WAS -125

#define EPOWF_LOG2E 1.442695040888963407359924681001892137426645954152985934135f

#define EPOWF_C1 0.693147180559945309417232121458176568075500134360255254120f
#define EPOWF_C2 0.240226506959100712333551263163332485865276475797272793433f
#define EPOWF_C3 0.055504108664821579953142263768621757359354414224767775398f
#define EPOWF_C4 0.009618129107628477161979071573658865479929845523272117657f
#define EPOWF_C5 0.001333355814642844342341222198799617473078225125878224529f

#pragma omp declare simd
static inline float epowF(const float x)
{
  const float x2 = x * EPOWF_LOG2E;

#if 1 // original code 
  const float offset = (x < 0.0f) ? -0.5f : 0.5f;
  const int ix = int(x2 + offset);
  const float d = x2 - float(ix);
#else // avoid int to float conversion (is slower !)
  const float c = 0x1.800000p+23f;
  const float rx2 = fmaf(1.0f, x2, c) - c;// rounds x2 to nearest integer value rx2
  const int ix = int(rx2);
  const float d = x2 - rx2;
#endif

#if 0 // original code
  const int ixc = ix < EPOWF_IXMAX ? ix : EPOWF_IXMAX;
  const int ie = ((ixc + EPOWF_IX) << EPOWF_SHIFT);
  const float &de = reinterpret_cast<const float &>(ie);
  const float dec = (ixc < EPOWF_IXMIN) ? 0.0f : de;
  const float d2 = 1.0f + d*(EPOWF_C1 + d*(EPOWF_C2 + d*(EPOWF_C3 + d*(EPOWF_C4 + d*EPOWF_C5))));/* polynomial for pow(2, d), where -0.5 <= d <= 0.5 */
  return d2*dec;
#else // preserve overflow to FLT_MAX by delaying bounds check (is faster !)
  const union { int ie; float de;} u = { ((ix + EPOWF_IX) << EPOWF_SHIFT) };
  const float d2 = 1.0f + d*(EPOWF_C1 + d*(EPOWF_C2 + d*(EPOWF_C3 + d*(EPOWF_C4 + d*EPOWF_C5))));/* polynomial for pow(2, d), where -0.5 <= d <= 0.5 */
  return (ix >= EPOWF_IXMAX) ? FLT_MAX : (ix < EPOWF_IXMIN) ? 0.0f : d2 * u.de;
#endif
}

#endif
