#include "structuralVariation.h"
#include "Calign.h"

// #undef DEBUG
// #define DEBUG 2

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/trunk/structuralVariation.cpp 5337 2016-09-18 09:02:09Z tanantharaman $");

//globals used only in output_smap.cpp and structuralVariation.h

const char* transintra_name           = "translocation_intrachr";
const char* transinter_name           = "translocation_interchr";
const char* transintra_overlap_name   = "trans_intrachr_overlap";
const char* transinter_overlap_name   = "trans_interchr_overlap";
const char* transintra_repeat_name    = "trans_intrachr_repeat";
const char* transinter_repeat_name    = "trans_interchr_repeat";
const char* transintra_duplicate_name = "trans_intrachr_duplicate";
const char* transinter_duplicate_name = "trans_interchr_duplicate";
const char* transintra_common_name    = "trans_intrachr_common";
const char* transinter_common_name    = "trans_interchr_common";
const char* transintra_segdup_name    = "trans_intrachr_segdupe";
const char* transinter_segdup_name    = "trans_interchr_segdupe";


const bool smap_xmapid_col = true; //include the xmap entry ids in smap; could be moved to parameters.cpp if promote to command line arg
const bool smap_linkid_col = true; //add column after XmapIDs for LinkID
const bool smap_labelind_col = true; //add four columns for label indices after XmapIDs and LinkID
const bool smap_xmap_allpairs = true; //don't just take neighboring xmap entries, take all (unique) pairs
const float sv_overlap_frac = 0.8; //fraction two svs must overlap to be called 'the same'
//int smap_genotype_ngroup = 0; //count number of genotypeGroups
const double translocation_overlap_buf = 50e3; //buffer around translocation breapoint for overlap

//need next three vars for tranlocation confidence: mean vector, 
//                               Occurrence.mean      Occurrence.max   Fragile.max.max    Occurrence.mean.bp
const double trans_mean_vec[] = {8.475774e-01,        1.163477e+00,    1.297754e+01,      8.117338e-01,
				 //Occurrence.max.bp  ChimQuality.max  nMG12              gapNlab      
				 8.583354e-01,        1.339591e+01,    2.035769e+00,      8.816215e-01,
				 //BPvrlpSize         nMG              gapSize 
				 3.212623e+02,        2.403168e+00,    1.346999e+04};
//covariance matrix for computing D^2
double cov_inv[] = {
2.344681e+02, -9.187475e+01, 4.671873e-01, -9.978937e+01, -6.022090e+00, 1.325726e-01, 1.277356e+00, -8.820539e-01, 2.006067e-04, -1.250822e+00, 9.053457e-05,
-9.187475e+01, 5.365974e+01, -3.707881e-01, 4.251637e+01, -2.071599e+01, -1.476644e-01, -6.556261e-01, 3.365950e-01, -9.327286e-05, 3.748683e-01, -8.692991e-07,
4.671873e-01, -3.707881e-01, 2.612783e-02, 2.501215e-01, -1.621776e-01, -2.165902e-03, -2.838879e-02, 1.112197e-02, -3.029364e-06, -1.514543e-03, -7.599997e-07,
-9.978937e+01, 4.251637e+01, 2.501215e-01, 5.623775e+02, -4.918152e+02, 1.995284e-01, -5.301193e+00, 1.616752e+01, -9.929470e-04, 4.281066e-03, 3.723877e-05,
-6.022090e+00, -2.071599e+01, -1.621776e-01, -4.918152e+02, 5.001455e+02, -1.229658e-01, 3.893722e+00, -1.475109e+01, 8.199814e-04, 5.611489e-01, -1.280411e-04,
1.325726e-01, -1.476644e-01, -2.165902e-03, 1.995284e-01, -1.229658e-01, 2.559985e-02, -2.224769e-02, 7.786831e-03, -2.199259e-06, 4.113321e-03, -3.237864e-07,
1.277356e+00, -6.556261e-01, -2.838879e-02, -5.301193e+00, 3.893722e+00, -2.224769e-02, 1.270530e+01, -1.623433e+00, -4.884736e-05, -2.344787e+00, 5.575866e-07,
-8.820539e-01, 3.365950e-01, 1.112197e-02, 1.616752e+01, -1.475109e+01, 7.786831e-03, -1.623433e+00, 4.903512e+00, -2.992718e-04, -9.226783e-02, -8.848676e-05,
2.006067e-04, -9.327286e-05, -3.029364e-06, -9.929470e-04, 8.199814e-04, -2.199259e-06, -4.884736e-05, -2.992718e-04, 1.565327e-07, -1.661547e-05, 1.020272e-08,
-1.250822e+00, 3.748683e-01, -1.514543e-03, 4.281066e-03, 5.611489e-01, 4.113321e-03, -2.344787e+00, -9.226783e-02, -1.661547e-05, 2.085475e+00, -5.044990e-07,
9.053457e-05, -8.692991e-07, -7.599997e-07, 3.723877e-05, -1.280411e-04, -3.237864e-07, 5.575866e-07, -8.848676e-05, 1.020272e-08, -5.044990e-07, 8.341426e-09};

//cdf for D^2 -> probability (confidence)
double ecdfx[] = {1.327981 , 1.734056 , 1.923069 , 2.069802 , 2.182068 , 2.295027 , 2.390234 , 2.485691 , 2.572938 , 2.661399 , 2.736775 , 2.824159 , 2.91665 , 2.996994 , 3.102753 , 3.17878 , 3.255857 , 3.342606 , 3.439689 , 3.516775 , 3.590741 , 3.672423 , 3.749547 , 3.828869 , 3.904792 , 3.978178 , 4.060434 , 4.13618 , 4.206156 , 4.291036 , 4.360105 , 4.426624 , 4.513181 , 4.586971 , 4.671488 , 4.762922 , 4.832319 , 4.9303 , 5.038754 , 5.122038 , 5.207412 , 5.293966 , 5.421972 , 5.537831 , 5.639657 , 5.765992 , 5.875916 , 5.979721 , 6.083787 , 6.228431 , 6.399593 , 6.520387 , 6.678645 , 6.818251 , 6.966424 , 7.093409 , 7.280555 , 7.441616 , 7.582972 , 7.746432 , 7.879245 , 8.046074 , 8.207415 , 8.383305 , 8.529455 , 8.707731 , 8.884509 , 9.123491 , 9.292371 , 9.491697 , 9.728372 , 10.01255 , 10.23273 , 10.49377 , 10.75922 , 11.01175 , 11.36122 , 11.70735 , 11.99954 , 12.36517 , 12.75183 , 13.09071 , 13.55096 , 14.03642 , 14.64201 , 15.24801 , 15.95631 , 16.72448 , 17.38185 , 18.08515 , 19.10594 , 20.33532 , 21.96433 , 23.92873 , 26.36668 , 29.65647 , 34.13201 , 41.5654 , 57.50552 , 92.1187};

//inversion confidence
//                               Occurrence.max    ChimNorm.max Fragile.max.max     nMG      nMG12        vrlpsz12      BPvrlpSize
const double inver_mean_vec[] = {9.664869e-01, 3.524110e+01, 1.354197e+01, 3.588671e+00, 2.139650e+00, 1.956001e+05, 1.369105e+03,
				 //gapNlab         gapSize         invSize         invNlab      minRawConf  nFarLab
				 2.280330e+00, 1.758491e+04, 3.543253e+05, 4.026097e+01, 5.098298e+01, 2.008829e+04};

//covariance matrix for inversions (above is translocations)
double cov_inv_inv[] = {
7.410747e+00, -3.156344e-03, 1.934357e-02, 1.276613e-03, -4.432703e-01, -6.177425e-06, -8.080355e-06, 7.894713e-02, -1.619922e-06, 3.230675e-07, -1.612661e-02, 6.848401e-03, -4.479051e-06,
-3.156344e-03, 1.370309e-03, -1.110220e-03, -3.032518e-03, 2.614407e-03, 3.163748e-08, -1.294696e-08, -9.112113e-04, 8.463249e-08, -3.472834e-08, -5.058310e-05, 6.380762e-05, 3.165679e-08,
1.934357e-02, -1.110220e-03, 2.235730e-02, 4.709375e-03, -1.817094e-02, -4.817148e-07, -6.052807e-07, -4.389719e-03, 7.234730e-07, -4.426235e-08, 5.208391e-04, -3.306058e-05, -1.273217e-06,
1.276613e-03, -3.032518e-03, 4.709375e-03, 2.451645e-01, -3.243969e-01, -6.750663e-07, -4.838933e-07, 1.719500e-03, 8.590444e-07, -4.576877e-07, 2.600072e-03, -2.556454e-03, -7.971123e-09,
-4.432703e-01, 2.614407e-03, -1.817094e-02, -3.243969e-01, 2.575869e+00, -1.451643e-06, -2.582531e-05, -3.749573e-02, -6.815943e-07, 1.606314e-06, -7.855433e-04, -3.110269e-05, -5.104618e-06,
-6.177425e-06, 3.163748e-08, -4.817148e-07, -6.750663e-07, -1.451643e-06, 4.894739e-10, 8.132374e-11, -3.889921e-06, -8.960324e-11, -1.340819e-10, 1.245581e-06, -1.548637e-07, -1.656530e-10,
-8.080355e-06, -1.294696e-08, -6.052807e-07, -4.838933e-07, -2.582531e-05, 8.132374e-11, 1.970097e-09, -1.117080e-05, 9.178936e-10, -4.885048e-11, 4.756278e-07, -2.070649e-07, -8.274574e-10,
7.894713e-02, -9.112113e-04, -4.389719e-03, 1.719500e-03, -3.749573e-02, -3.889921e-06, -1.117080e-05, 6.273107e-01, -5.521343e-05, 2.620798e-06, -2.829185e-02, 5.525098e-03, 3.419240e-05,
-1.619922e-06, 8.463249e-08, 7.234730e-07, 8.590444e-07, -6.815943e-07, -8.960324e-11, 9.178936e-10, -5.521343e-05, 5.475091e-09, -1.772821e-10, 1.385788e-06, 1.841162e-09, -3.152527e-09,
3.230675e-07, -3.472834e-08, -4.426235e-08, -4.576877e-07, 1.606314e-06, -1.340819e-10, -4.885048e-11, 2.620798e-06, -1.772821e-10, 2.186753e-10, -1.677580e-06, -1.083235e-07, -9.011601e-11,
-1.612661e-02, -5.058310e-05, 5.208391e-04, 2.600072e-03, -7.855433e-04, 1.245581e-06, 4.756278e-07, -2.829185e-02, 1.385788e-06, -1.677580e-06, 1.841924e-02, -2.727464e-03, 9.397678e-07,
6.848401e-03, 6.380762e-05, -3.306058e-05, -2.556454e-03, -3.110269e-05, -1.548637e-07, -2.070649e-07, 5.525098e-03, 1.841162e-09, -1.083235e-07, -2.727464e-03, 2.990654e-03, 4.813759e-08,
-4.479051e-06, 3.165679e-08, -1.273217e-06, -7.971123e-09, -5.104618e-06, -1.656530e-10, -8.274574e-10, 3.419240e-05, -3.152527e-09, -9.011601e-11, 9.397678e-07, 4.813759e-08, 9.423796e-09};

//cdf values are different for inversions, too (algo is same)
double ecdfx_inv[] = {0.000000, 1.179212, 1.556351, 1.759510, 1.931891, 2.053569, 2.177434, 2.296571, 2.417840, 2.530012, 2.668030, 2.755910, 2.867986, 2.970680, 3.077678, 3.182430, 3.255534, 3.333596, 3.404182, 3.490746, 3.564960, 3.660836, 3.757613, 3.849160, 3.936571, 4.006723, 4.080853, 4.141076, 4.211544, 4.309972, 4.401351, 4.489169, 4.580247, 4.668514, 4.754075, 4.838443, 4.934211, 5.031668, 5.112287, 5.187984, 5.260692, 5.372681, 5.462173, 5.577066, 5.670361, 5.767454, 5.873470, 5.990387, 6.115443, 6.234778, 6.333973, 6.430796, 6.557156, 6.680984, 6.846087, 6.975087, 7.086531, 7.167944, 7.303815, 7.441668, 7.575359, 7.765069, 7.942601, 8.122321, 8.305435, 8.490242, 8.782865, 9.005400, 9.178808, 9.359080, 9.614720, 9.774810, 10.048189, 10.284763, 10.537549, 10.795655, 11.138676, 11.375123, 11.637166, 11.926350, 12.322332, 12.762691, 13.168662, 13.612391, 14.044004, 14.492835, 15.070250, 15.691794, 16.324981, 17.242256, 18.246238, 19.355080, 20.566574, 22.302625, 24.018374, 26.320804, 29.002722, 33.830106, 40.192365, 52.612010, 104.005368};


extern Cmap **YYmap,**XXmap;
extern int numY,numX;

//write all data to first argument--a file handle
void writeToSmap(FILE *fp, int nsv, long long qrycontigid, long long refcontigid1, long long refcontigid2, double qrystart, double qrystop, double refstart, double refstop, const char *svtype, int xmapid1, int xmapid2, int linkid, int qsi, int qei, int rsi, int rei, double conf, double confl, double confr, double confc, const char* zyg, int gen, int ggrp, float confs, double SVsize, double SVfreq, double SVcov, double SVtotcov, char *orientation) {
  assert(qrycontigid > 0); //using int for qrycontigid causes a long long to wrap to negative
  //qrycontigid and refcontigid don't change from xmap
  fprintf(fp, "%d\t%lld\t%lld\t%lld\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%s", nsv, //normal
  //fprintf(fp, "%-2d\t%4i\t%2i\t%2i\t%9.1f\t%9.1f\t%10.1f\t%10.1f\t%.2f\t%24s", nsv, //better viewing
	  qrycontigid,
	  refcontigid1,
	  refcontigid2,
	  qrystart, //usexmapentries[i]->qrystartpos,
	  qrystop, //usexmapentries[i]->qryendpos,
	  refstart, //usexmapentries[i-1]->refendpos,
	  refstop, //usexmapentries[i]->refstartpos,
	  //'+', //(usexmapentries[i]->orientforward ? '+' : '-'),
	  confs, //usexmapentries[i]->confidence
	  svtype
	  );
  if( smap_xmapid_col ) {
    fprintf(fp, "\t%i", xmapid1 > 0 ? xmapid1 : -1);
    fprintf(fp, "\t%i", xmapid2 > 0 ? xmapid2 : -1);
  }
  if( smap_linkid_col )
    fprintf(fp, "\t%i", linkid);
  if( smap_labelind_col )
    fprintf(fp, "\t%i\t%i\t%i\t%i", qsi, qei, rsi, rei);
  if( smap_zygosity ) { //include genotype with zygosity
    fprintf(fp, "\t%s", zyg != 0 ? zyg : "unknown");
    fprintf(fp, "\t%i\t%i", gen, ggrp);
  }
  fprintf(fp, "\t%.2f\t%.2f\t%.2f\t%.2f", conf, confl, confr, confc);
  if(SV_SIZE || SV_FREQ || SV_ORIENT)
    fprintf(fp, "\t%0.1f", SV_SIZE ? SVsize : -1.0); 
  if(SV_FREQ || SV_ORIENT){
    fprintf(fp, "\t%0.3f", SV_FREQ ? SVfreq : -1.0);
    if(SV_FREQ >= 2)
      fprintf(fp,"\t%2.2f\t%2.2f", SVcov, SVtotcov);
  }
  if(SV_ORIENT){
    if(orientation != NULL)
      fprintf(fp,"\t%3s", orientation);
    else
      fprintf(fp,"\tNA ");
  }
  fprintf(fp, "\n");
}

//for qsort'ing structuralVariation objects
int compare_svs(const void *a, const void* b) {
  structuralVariation *x = (structuralVariation*)a;
  structuralVariation *y = (structuralVariation*)b;
  if( x->xmapid1 < y->xmapid1 ) //first compare xmapID1
    return -1;
  if( x->xmapid1 > y->xmapid1 )
    return 1;
  if( x->xmapid2 < y->xmapid2 ) //then compare xmapID2
    return -1;
  if( x->xmapid2 > y->xmapid2 )
    return 1;

  /* Need to first break ties in xmapid1 by algo_enum, so paired sv_indel type (with same xmapid1,xmapid2) stay together, in case one of the other SVs is using the pair of xmapid's */

  /* sort algo_enum in descending order so sv_indel (= 1) comes before sv_ps (= 0) */
  if(x->algo_enum > y->algo_enum)
    return -1;
  if(x->algo_enum < y->algo_enum)
    return 1;

  //all PS entries should be sorted by above because there are no duplicates.
  // however, indel entries will be sorted by ref (remaining ifs)
  if( x->refcontigid1 < y->refcontigid1 ) //next compare ref id 1
    return -1;
  if( x->refcontigid1 > y->refcontigid1 )
    return 1;
  if( x->refcontigid2 < y->refcontigid2 ) //next compare ref id 2
    return -1;
  if( x->refcontigid2 > y->refcontigid2 )
    return 1;

  if( x->refstart < y->refstart ) //if equal, compare refstart
    return -1;
  if( x->refstart > y->refstart )
    return 1;
  if( x->refstop < y->refstop ) //if equal, compare refstop
    return -1;
  if( x->refstop > y->refstop )
    return 1;

  /* a match on all the above fields should never happen (except for paired Smalled Inversions ?), but if it does maintain original order */
  if( x < y)
    return -1;
  if( x > y)
    return 1;// WAS29 -1
  return 0; // should never happen
}

double makeConfidence(double a, double b, double c)
{
  //use max 0 bc of extreme case where all three confs are large
  double MinValue = min(a,min(b,c));
  double Pvalue = pow(10.0, MinValue - a) + pow(10.0, MinValue - b) + pow(10.0, MinValue - c);
  if(DEBUG) assert(isfinite(Pvalue) && Pvalue >= 1.0);

  double Confidence = MinValue - log10(Pvalue);
  if(DEBUG && !(isfinite(Confidence))){
    printf("makeConfidence(%0.6f,%0.6f,%0.6f):Confidence=%0.6e\n",a,b,c,Confidence);
    fflush(stdout);

    assert(isfinite(Confidence));
  }

  return Confidence;
}


bool overlapPoint(double min1, double max1, double p) {
  return min1 < p && max1 > p;
}

//utility fn to get overlap length if there is overlap
double overlapLength(double min1, double max1, double min2, double max2) {
  return min(max1,max2) - max(min1,min2); //<= 0 if no overlap. Negative values are size of gap 
}


//analogous to alignmentOverlap, but for translocations
void alignmentOverlapTranslocation(xmapEntry **XMAPentries, int numXMAPentries) {// NOTE : XMAPentries[] could be usexmapentries[] or xmapentries[]

  if( !smap_zygosity )
    return;

#if LATE_ZYGOSITY==0
  assert(XMAPentries == usexmapentries);
#else
  assert(XMAPentries == xmapentries);
#endif

  for(int indp=0; indp < numsv; indp++) {
    structuralVariation &sv = sv_array[indp];

    if( (sv.type_enum != intrachr_translocation &&
	 sv.type_enum != interchr_translocation)
	|| sv.duplicate ) {
      //printf("Skipping sv %i\n", indp+1); //debug
      continue;
    }

    for(int ie=0; ie < 2; ie++) { //loop over 2 xmap entries which make the translocation
      
      xmapEntry *entry1 = (sv.algo_enum == sv_ps) ? usexmapentries[sv.indices[ie]] : xmapentries[sv.indices[ie]];
      double p = (ie == 0 ? sv.refstart : sv.refstop); //translocation breakpoint on entry1
      if(0 && DEBUG && !(p == entry1->refstartpos || p == entry1->refendpos)){/* NOTE : this condition is not true for translocation called by one MG overlapping insertion in other MG */
	printf("indp=%d,numsv=%d: sv.algo_enum= %d, sv.refstart= %0.3f, sv.refstop= %0.3f kb, sv.indices[ie=%d] = %d\n",
	       indp,numsv,sv.algo_enum,sv.refstart*1e-3,sv.refstop*1e-3,ie,sv.indices[ie]);
	fflush(stdout);
	assert(p == entry1->refstartpos || p == entry1->refendpos);
      }

      for(int i= 0; i < numXMAPentries; i++) {
	xmapEntry *pxmap = XMAPentries[i];
	if(!pxmap->output)
	  continue;// NEW148: skip XMAPs that will NOT be output

	if( pxmap->refcontigid != entry1->refcontigid || //must overlap on reference
	    pxmap->qrycontigid == entry1->qrycontigid ) { //a contig cannot overlap itself
	  //printf("Skipping xmapEntry %i\n", pxmap->xmapid); //debug
	  continue;
	}

	//The condition below is too loose: nearby alignments shouldn't count as overlap because, eg, the other side
	// of the translocation is close enough to cause it to be called overlap in this case when it should not
	//instead, require entire region of p+-buf to be spanned
	if( pxmap->refstartpos > p - translocation_overlap_buf ||
	    pxmap->refendpos   < p + translocation_overlap_buf ) 
	//if( !overlapPoint(pxmap->refstartpos-translocation_overlap_buf,
	//		  pxmap->refendpos+translocation_overlap_buf, p) )
	  continue;

	// HERE HERE : mappings from usexmapentries[xi] to xmapentries[yi] should be precomputed for better performance

#if LATE_ZYGOSITY == 0
	int xi = i;
#else
	/* locate xi s.t. XMAPentries[i] == usexmapentries[xi] (may not be present) */
	int xi = -1;
	for(int t = 0; t < smapsize; t++)
	  if(usexmapentries[t] == pxmap){
	    xi = t;
	    break;
	  }
#endif
	if(VERB>=2){
	  if(xi >= 0)
	    printf("translocationOverlap: indp=%4i ie=%i: qryid= %4lld, refid= %2lld: ref= %0.3f .. %0.3f kb (p= %0.3f kb): i=%4i, xi=%d (qryid= %4lld ref= %0.3f .. %0.3f kb)\n", 
		   indp, ie, entry1->qrycontigid, entry1->refcontigid, entry1->refstartpos * 1e-3, entry1->refendpos* 1e-3, p * 1e-3, i, xi, pxmap->qrycontigid, pxmap->refstartpos*1e-3, pxmap->refendpos*1e-3);
	  else
	    printf("translocationOverlap: indp=%4i ie=%i: qryid= %4lld, refid= %2lld: ref= %0.3f .. %0.3f kb (p= %0.3f kb): i=%4i (qryid= %4lld ref= %0.3f .. %0.3f kb)\n", 
		   indp, ie, entry1->qrycontigid, entry1->refcontigid, entry1->refstartpos * 1e-3, entry1->refendpos* 1e-3, p * 1e-3, i, pxmap->qrycontigid, pxmap->refstartpos*1e-3, pxmap->refendpos*1e-3);
	  fflush(stdout);
	}

	if(xi < 0){/* XMAPentries[i] cannot be part of a translocation call */
	  sv.alignment_overlap = true;
	  break; //once one overlap is found, no need to check for more (bool) -- this is only if no indel overlap
	}

	//now must check if this xmap entry is part of a translocation because
	// if it has a similar breakpoint, then the alignment overlap doesn't count
	bool haveoverlap = true;

	for(int indr= 0; indr < numsv; indr++) { //same as outer loop
	  structuralVariation &sv2 = sv_array[indr];

	  //only interested in xmap entry xi: must be translocation
	  if( (sv2.type_enum != intrachr_translocation &&
	       sv2.type_enum != interchr_translocation) ||
	      sv2.duplicate || indp == indr || //sv can't overlap itself
	      (sv2.indices[0] != xi && sv2.indices[1] != xi) )
	    //pxmap->refcontigid != sv2.refcontigid1 ) //must be on same ref
	    continue;
	  
	  xmapEntry *entry2 = 0; //not really necessary (only used in assertion)
	  double p2 = 0;
	  if( sv2.indices[0] == xi ) {
	    entry2 = usexmapentries[sv2.indices[0]];
	    p2 = sv2.refstart;
	  } else if( sv2.indices[1] == xi ) {
	    entry2 = usexmapentries[sv2.indices[1]];
	    p2 = sv2.refstop;
	  } else
	    assert(0); //dummy check

	  if(DEBUG && pxmap->xmapid != entry2->xmapid){
	    printf("translocationOverlap: indp=%4i ie=%i: qryid= %4lld, refid= %2lld: ref= %0.3f .. %0.3f kb (p= %0.3f kb): i=%4i, xi=%d (qryid= %4lld ref= %0.3f .. %0.3f kb)\n", 
		   indp, ie, entry1->qrycontigid, entry1->refcontigid, entry1->refstartpos * 1e-3, entry1->refendpos* 1e-3, p * 1e-3, i, xi, pxmap->qrycontigid, pxmap->refstartpos*1e-3, pxmap->refendpos*1e-3);
	    printf("\t pxmap->xmapid= %d, entry2->xmapid= %d, pxmap= %p, entry2= %p\n",pxmap->xmapid,entry2->xmapid, pxmap, entry2);
	    fflush(stdout);
	    assert( pxmap->xmapid == entry2->xmapid ); //dummy check: index checked above
	  }

	  //must require that the breakpoint from indr (p2) is near the breakpoint from indp (p)
	  //otherwise, the alignment overlap still counts because the two translocation breakpoints are not similar
	  //within 2*overlapbuf is a bit loose, but better to make it too large rather than too small
	  bool ovlp = overlap( p-translocation_overlap_buf, p+translocation_overlap_buf,
			       p2-translocation_overlap_buf, p2+translocation_overlap_buf );
	  if( ovlp ) {
	    haveoverlap = false;
	    printf("translocationOverlap: sv2: %i: gm %4lld chr %2lld, %.1f; chr %2lld, %.1f\n", indr+1, entry2->qrycontigid, sv2.refcontigid1, sv2.refstart, sv2.refcontigid2, sv2.refstop);
	    break; //since haveoverlap is false, nothing else can change
	  }
	} //end loop on sv2

	if( haveoverlap ) {
	  sv.alignment_overlap = haveoverlap; 
	  break; //once one overlap is found, no need to check for more (bool) -- this is only if no indel overlap
	}

      } //end for all xmap entries
    } //end for translocation breakpoints
  } //end for svs
} //end alignmentOverlapTranslocation


//analogous to alignmentOverlap, but for inversions
void alignmentOverlapInversion(xmapEntry **XMAPentries, int numXMAPentries) {// NOTE : XMAPentries[] could be usexmapentries[] or xmapentries[]

  if( !smap_zygosity )
    return;

#if LATE_ZYGOSITY==0
  assert(XMAPentries == usexmapentries);
#else
  assert(XMAPentries == xmapentries);
#endif

  for(int indp=0; indp < numsv; indp++) {
    structuralVariation &sv = sv_array[indp];

    if( (sv.type_enum != inversion_type &&
	 sv.type_enum != inversion_small_type)
	|| sv.duplicate ) {
      //printf("Skipping sv %i\n", indp+1); //debug
      continue;
    }

    int i1 = (sv.link_index == sv.indices[0] ? sv.indices[1] : sv.indices[0]); //the non inverted mg
    xmapEntry *entry1 = (sv.algo_enum == sv_ps) ? usexmapentries[i1] : xmapentries[i1];
    //for inversions, require overlap of both breakpoints separately
    bool ovlp_start = false, ovlp_end = false;

    for(int i= 0; i < numXMAPentries; i++) {
      xmapEntry *pxmap = XMAPentries[i];
      if(!pxmap->output)
	continue;// NEW148: skip XMAPs that will NOT be output

      if( pxmap->refcontigid != entry1->refcontigid || //must overlap on reference
	  pxmap->qrycontigid == entry1->qrycontigid ) { //a contig cannot overlap itself
	//printf("Skipping xmapEntry %i\n", pxmap->xmapid); //debug
	continue;
      }
      //currently not actually using both, just as an ||
      bool thisovlp_start = (pxmap->refstartpos < sv.refstart - translocation_overlap_buf &&
			     pxmap->refendpos   > sv.refstart + translocation_overlap_buf );
      bool thisovlp_end   = (pxmap->refstartpos < sv.refstop - translocation_overlap_buf &&
			     pxmap->refendpos   > sv.refstop + translocation_overlap_buf );

      if( !thisovlp_start && !thisovlp_end )
	continue;

#if LATE_ZYGOSITY == 0
      /* locate yi s.t. XMAPentries[i] == xmapentries[yi] (may not be present) */
      int xi = i;
      int yi = -1;
      for(int t = 0; t < numxmapentries; t++)
	if(xmapentries[t] == pxmap){
	  yi = t;
	  break;
	}
#else
      /* locate xi s.t. XMAPentries[i] == usexmapentries[xi] (may not be present) */
      int yi = i;
      int xi = -1;
      for(int t = 0; t < smapsize; t++)
	if(usexmapentries[t] == pxmap){
	  xi = t;
	  break;
	}
#endif

      if(VERB >= 2){
	printf("inversionOverlap : indp=%4i(xid=%d,%d): qryid=%4lld refid=%2lld: ref=%0.3f .. %0.3f kb: i=%4i,xi=%d,yi=%d: qryid=%4lld ref=%0.3f ... %0.3f kb\n", 
	       indp, sv.indices[0], sv.indices[1], entry1->qrycontigid, entry1->refcontigid, entry1->refstartpos*1e-3, entry1->refendpos*1e-3, 
	       i, xi, yi, pxmap->qrycontigid, pxmap->refstartpos*1e-3, pxmap->refendpos*1e-3);
	fflush(stdout);
      }

      //if this xmapentry is part of an inversion which is 'similar', then the overlap doesn't count (exclude)
      int indr=0;
      structuralVariation *sv2 = NULL;
      for(; indr < numsv; indr++) { //same as outer loop
	structuralVariation &svr = sv_array[indr];
	//only interested in xmap entry xi (or yi): must be inversion which is 'close' to whichever of thisovlp_start/end
	// in order for overlap to not count (because that will be counted as inversion overlap inversion)
	if( (svr.type_enum != inversion_type &&
	     svr.type_enum != inversion_small_type) ||
	    svr.duplicate || indp == indr ) //sv can't overlap itself
	  continue;

	sv2 = &svr;
	xmapEntry *entry2 = (sv2->algo_enum == sv_ps) ? usexmapentries[sv2->indices[0]] : xmapentries[sv2->indices[0]];
	//though checking xi is correct, a necessary condition is the GM must be different
	if( entry2->qrycontigid == entry1->qrycontigid ) 
	  continue;
	//normal inversion:
	if( sv2->type_enum == inversion_type ) {// use xi index
	  if( sv2->indices[0] != xi && sv2->indices[1] != xi ) //this sv must be on xi in order to exclude
	    continue;
	}
	if(sv2->type_enum == inversion_small_type){// use yi index
	  if( sv2->indices[0] != yi && sv2->indices[1] != yi ) //this sv must be on yi in order to exclude
	    continue;
        }

	//check if it is 'close': any overlap from start-stop +- 50 kb is 'close' enough
	if( overlap(min(sv.refstart, sv.refstop)-translocation_overlap_buf,
		    max(sv.refstart, sv.refstop)+translocation_overlap_buf,
		    min(sv2->refstart, sv2->refstop)-translocation_overlap_buf,
		    max(sv2->refstart, sv2->refstop)+translocation_overlap_buf) ) {
	  thisovlp_start = thisovlp_end = false;
	  break;
	}
	
      } //inner sv loop
      
      if(sv2 == NULL)/* NEW72 : this happens if there is only one inversion SV */
	break;

      if(DEBUG>=1+RELEASE) assert(sv2 != NULL);

      if( !thisovlp_start && !thisovlp_end ){ //this means overlap was excluded
	if(VERB>=2){
	  printf("inversionOverlap2: %4i: chr %2lld: %10.1f %10.1f idx: %i %i\n", indr+1, sv2->refcontigid1, sv2->refstart, sv2->refstop, sv2->indices[0], sv2->indices[1]);
	  fflush(stdout);
	}
      } else {
	ovlp_start = (ovlp_start || thisovlp_start);
	ovlp_end   = (ovlp_end   || thisovlp_end);
      }

    } //xmapentries

    if( ovlp_start && ovlp_end )
      sv.alignment_overlap = true;

  } //outer sv loop
} //end alignmentOverlapInversion

// Check all svs (sv_array) for overlap with all alignments in usexmapentries
// For the purpose of zygosity calling, a single contig can't have two alleles
// usexmapentries must be the overlap-filtered array created in output_smap
// Currently only check indels; if overlap found, check that there is not
// also an indel at the same location (to avoid redundancy with svOverlap)

void alignmentOverlap(xmapEntry **XMAPentries, int numXMAPentries) {// NOTE : XMAPentries[] could be usexmapentries[] or xmapentries[]

  if( !smap_zygosity )
    return;

#if LATE_ZYGOSITY==0
  assert(XMAPentries == usexmapentries);
#else
  assert(XMAPentries == xmapentries);
#endif

  for(int indp = 0; indp < numsv; indp++) {
    structuralVariation &sv = sv_array[indp];// Thomas : use reference, rather than deep copy

    //only check indels, also, do not bother for duplicates
    if( (sv.type_enum != insertion && sv.type_enum != deletion && sv.type_enum != compression) ||
	sv.duplicate ) {
      //printf("Skipping sv %i\n", indp+1); //debug
      continue;
    }

    if(DEBUG>=2 && !(sv.refstart <= sv.refstop) ) {    // refstart can be == refstop, otherwise must be < -- not used until overlapLength call
      printf("ERROR: sv : %i, %.1f - %.1f; xmapids %i %i\n", indp+1, sv.refstart, sv.refstop, sv.xmapid1, sv.xmapid2);
      fflush(stdout);
      assert(sv.refstart <= sv.refstop);
    }

    //if equal, add 1.0 bp so that overlapLength returns 1.0; otherwise it would return 0.0 and these calls would be misclassified
    double stop = sv.refstop + (sv.refstart < sv.refstop ? 0.0 : 1.0);

    xmapEntry *entry1 = (sv.algo_enum == sv_ps) ? usexmapentries[sv.indices[0]] : xmapentries[sv.indices[0]];
    if(VERB>=2 && entry1->refcontigid == 21 && entry1->qrycontigid == 11){
      printf("sv %d, type=%d, ref= %0.3f .. %0.3f(%0.3f) kb, xmapid=%d,%d, refid= %lld, qryid= %lld: checking for overlapping xmapids\n",
	     indp+1,sv.type_enum,sv.refstart*1e-3,sv.refstop*1e-3,stop*1e-3,sv.xmapid1,sv.xmapid2, entry1->refcontigid, entry1->qrycontigid);
      fflush(stdout);
    }

    for(int xi= 0; xi < numXMAPentries; xi++) {
      xmapEntry *pxmap =  XMAPentries[xi];
      if(!pxmap->output)
	continue;// NEW148 : skip matchgroups that will not be output (for whatever reason)

      if( pxmap->refcontigid != sv.refcontigid1 || //must overlap on reference
	  pxmap->qrycontigid == entry1->qrycontigid ) { //a contig cannot overlap itself
	if(VERB>=2 && entry1->refcontigid == 21 && entry1->qrycontigid == 11 && pxmap->qrycontigid == 35 && pxmap->refcontigid == 21){
	  printf("\t xi=%d/%d: refid=%lld,qryid=%lld, ref= %0.3f .. %0.3f kb: skipping due to mismatched contig ids\n", 
		 xi, numXMAPentries, pxmap->refcontigid,pxmap->qrycontigid,pxmap->refstartpos*1e-3,pxmap->refendpos*1e-3);	  
	  fflush(stdout);
	}

	continue;
      }

      double ovlp = overlapLength( min(pxmap->refstartpos, pxmap->refendpos),
				   max(pxmap->refstartpos, pxmap->refendpos),
				   sv.refstart, stop);

      if( ovlp <= 0.0 ){
	if(VERB>=2 && entry1->refcontigid == 21 && entry1->qrycontigid == 11 && pxmap->qrycontigid == 35){
	  printf("\t xi=%d: refid=%lld,qryid=%lld, ref= %0.3f .. %0.3f kb, overlap= %0.3f kb : skipping\n", xi, pxmap->refcontigid,pxmap->qrycontigid,pxmap->refstartpos*1e-3,pxmap->refendpos*1e-3,ovlp*1e-3);
	  fflush(stdout);
	}
	continue;
      }
      
      double denom = fabs(stop - sv.refstart); //overlap is a fraction of the sv's size on reference

      //postpone print until after sv overlap check (?) -- just print both separately

      //can move this above printf immediately above if don't want to see partial overlaps
      if( ovlp < denom * smap_Indel_Ref_Alignment_OverlapFrac /* WAS28 sv_overlap_frac */ ) {
	if(VERB>=2 && entry1->refcontigid == 21 && entry1->qrycontigid == 11){
	  printf("\t alignmentOverlap: sv= %i, type=%d, ref=%.3f - %.3f jb; xmapid1= %i, ref=%.3f - %.3f kb; ovlp: %.3f, %0.3f kb, ratio= %.6f, thresh= %0.3f: failed, skipping\n", 
		 indp+1, sv.type_enum,sv.refstart*1e-3, sv.refstop*1e-3, pxmap->xmapid, pxmap->refstartpos*1e-3, pxmap->refendpos*1e-3, ovlp*1e-3, denom*1e-3, ovlp/denom, smap_Indel_Ref_Alignment_OverlapFrac);
	  fflush(stdout);
	}

	continue;
      }

      if(VERB>=2 && entry1->refcontigid == 21 && entry1->qrycontigid == 11){
	printf("\t alignmentOverlap: sv= %i, type=%d, ref=%.3f - %.3f jb; xmapid1= %i, ref=%.3f - %.3f kb; ovlp: %.3f, %0.3f kb, ratio= %.6f, thresh= %0.3f: passed, checking for sv overlap\n", 
	       indp+1, sv.type_enum,sv.refstart*1e-3, sv.refstop*1e-3, pxmap->xmapid, pxmap->refstartpos*1e-3, pxmap->refendpos*1e-3, ovlp*1e-3, denom*1e-3, ovlp/denom, smap_Indel_Ref_Alignment_OverlapFrac);
	fflush(stdout);
      }
      bool haveoverlap = true;

      //need to check for -indel overlap on this alignment bc if there is one,
      // it's an sv overlap, not an indel overlap--don't want to double count
      //only -indel can be in the middle of a matchgroup, so only these are possible
      for(int indr = 0; indr < numsv; indr++) { //same as outer loop
	structuralVariation &sv2 = sv_array[indr];

	//the type_enum and duplicate checking are not necessary since these are always true and
	// false respectively for -indel calls (just check algo_enum instead)
	//if( svr.algo_enum != sv_indel
	//change this: for the case of both split indel and alignment overlap on the same GM, do not
	// count the alignment overlap--only the sv overlap.


	if( (sv2.type_enum != insertion && sv2.type_enum != deletion && sv2.type_enum != compression) ||
	    sv2.duplicate || indp == indr || //sv can't overlap itself
	    pxmap->refcontigid != sv2.refcontigid1 ) //must be on same ref
	  continue;

	xmapEntry *entry2 = (sv2.algo_enum == sv_ps) ? usexmapentries[sv2.indices[1]] : xmapentries[sv2.indices[1]];

	if( pxmap->qrycontigid != entry2->qrycontigid ) //must be same query contig
	  continue;
	
	//check if these two svs are candidates for 'the same' (and therefore the alignment overlap is redundant):
	ovlp = overlapLength( sv.refstart, stop, sv2.refstart, sv2.refstop );
	
	if( ovlp > smap_Indel_Ref_SV_Overlap /* WAS28 0 */ ) {
	  haveoverlap = false;
	  if(VERB>=2 && entry1->refcontigid == 21 && entry1->qrycontigid == 11){
	    printf("\t SVoverlap: sv2=%i,type=%d, ref= %.3f .. %.3f kb; ovlp= %0.4f kb\n", indr+1, sv2.type_enum,sv2.refstart*1e-3, sv2.refstop*1e-3, ovlp*1e-3);
	    fflush(stdout);
	  }
	  break; //since haveoverlap is false, nothing else can change
	}
      } //end loop on sv2

      if( haveoverlap ) {
	sv.alignment_overlap = haveoverlap; 
	if(VERB>=2 && entry1->refcontigid == 21 && entry1->qrycontigid == 11){
	  printf("\t alignmentOverlap confirmed: sv= %d, type= %d, ref= %0.3f .. %0.3f kb\n",indp+1,sv.type_enum,sv.refstart*1e-3,sv.refstop*1e-3);
	  fflush(stdout);
	}
	break; //once one overlap is found, no need to check for more (bool) -- this is only if no indel overlap
      }
    } //end loop on usexmapentries
  } //end loop on svs
} //end alignmentOverlap


int oppositeGenotype(int g1) {
  if( g1 == 1 )
    return 2;
  else if( g1 == 2 )
    return 1;
  printf("Invalid genotype: %i\n", g1);
  assert(false);
}


// call zygosity using overlap of SVs with each other which are called on
// different contigs--originally indels only, add translocations
// here overlap is in reference coordinates, like alignmentOverlap
// Indels are only compared if RefOverlap exceeds smap_Indel_Ref_SV_Overlap (defaults to -20kb)
// Result of comparision are stored in indel_overlap data member of sv objects
//   There are two types of results: 'same', ie, An indel which appears
//     to be the same indel on a different contig, and 'different': overlapping an indel which is not the same.
//   Similarity of Indels is based on having the same type and an Indel sizes ratio (smaller/larger) exceeding smap_Indel_Size_MatchRation (defaults to 0.8)
//
// Interpretation of overlap:
// 0: no overlap
// 1: overlap same and not different 
// 2: overlap different and not same 
// 3: overlap both same and different OR multiple different
//see structuralVariation.getZygosity for how this is interpreted
//also add genotype (-1, 1, or 2) and genotypeGroup for overlapping indels
void svOverlap(xmapEntry **usexmapentries, int numxmapent) {

  if( !smap_zygosity )
    return;

  int smap_genotype_ngroup = 0; //count number of genotypeGroups

  for(int indp=0; indp < numsv; indp++) {
    structuralVariation &sv = sv_array[indp];

    if( DEBUG>=1+RELEASE && (sv.type_enum == insertion || sv.type_enum == deletion || sv.type_enum == compression) &&
	!sv.duplicate && sv.confidence < 0 ) { //dummy check that confidence was set--may be redundant?
      printf("WARNING in svOverlap: sv %i, %i: SVconf=%f: %2lld %10.1f %10.1f %10.1f\n", indp+1, sv.algo_enum, sv.confidence, sv.refcontigid1, sv.refstart, sv.refstop, sv.querysize);
      fflush(stdout);
      //assert(sv.confidence >= 0);
    }

    //now checking translocations and inversions, too; do not bother for duplicates or overlaps (exclude)
    //if this indel is already overlapping both types, nothing more to check
    //also do not consider indels which do not pass the confidence threshold
    if( sv.duplicate || sv.exclude || sv.indel_overlap == 3 )
      continue;

    bool isindel = (sv.type_enum == insertion || sv.type_enum == deletion || sv.type_enum == compression);
    bool istrans = (sv.type_enum == intrachr_translocation ||
		    sv.type_enum == interchr_translocation );
    bool isinv   = (sv.type_enum == inversion_type ||
		    sv.type_enum == inversion_small_type);

    if( isindel && sv.confidence < smap_min_conf )
      continue;

    if( !isindel && !istrans && !isinv )
      continue;

    xmapEntry *entry1 = (sv.algo_enum == sv_ps) ? usexmapentries[sv.indices[0]] : xmapentries[sv.indices[0]];
    
    for(int indr = indp+1; indr < numsv; indr++) { //compare to subsequent; cannot skip based on existing overlap
      if( (isindel && sv_array[indr].type_enum != insertion && sv_array[indr].type_enum != deletion && sv_array[indr].type_enum != compression) ||
	  (istrans && sv_array[indr].type_enum != intrachr_translocation && sv_array[indr].type_enum != interchr_translocation ) ||
	  (isinv   && sv_array[indr].type_enum != inversion_type && sv_array[indr].type_enum != inversion_small_type) ||
	  sv_array[indr].duplicate || sv_array[indr].exclude ||
	  //careful with refcontigids: for indels and inversions, 1==2, so just check 1; for trans, either matching is enough
	  ((isindel || isinv) && sv.refcontigid1 != sv_array[indr].refcontigid1) ||
	  (istrans && sv.refcontigid1 != sv_array[indr].refcontigid1 &&
	   sv.refcontigid2 != sv_array[indr].refcontigid2) ||
	  (isindel && sv_array[indr].confidence < smap_min_conf) )
	continue;

      structuralVariation &sv2 = sv_array[indr];
      xmapEntry *entry2 = (sv2.algo_enum == sv_ps) ? usexmapentries[sv2.indices[0]] : xmapentries[sv2.indices[0]];
      if( entry1->qrycontigid == entry2->qrycontigid ) //no overlap allowed on same contig
	continue;

      bool same = false;
      double ovlp=0, IndelRatio=0;
      bool ol1=0, ol2=0, d11=0, d21=0, d12=0, d22=0;
      if( isindel ) {
	//same as in alignmentOverlap: these two svs are candidates for 'the same' if overlap length exceeds smap_Indel_Ref_SV_Overlap (kb).
	ovlp = overlapLength( sv.refstart, sv.refstop, sv2.refstart, sv2.refstop );
	if( ovlp <= smap_Indel_Ref_SV_Overlap /* WAS28 0. */ ) //nothing to do if insufficient overlap
	  continue;

	// NOTE: refsize and querysize are in bp for *internal outliers* and in kb for *split indels* (convert to kb here)
	double ref_size1 = sv.refsize * (sv.algo_enum == sv_ps ? 1.0 : 1.0e-3);
	double ref_size2 = sv2.refsize * (sv.algo_enum == sv_ps ? 1.0 : 1.0e-3);
	double qry_size1 = sv.querysize  * (sv.algo_enum  == sv_ps ? 1.0 : 1.0e-3);
	double qry_size2 = sv2.querysize * (sv2.algo_enum == sv_ps ? 1.0 : 1.0e-3);
	double Indel_size1 = fabs(ref_size1 - qry_size1);
	double Indel_size2 = fabs(ref_size2 - qry_size2);

	IndelRatio = min(Indel_size1,Indel_size2)/max(0.001,max(Indel_size1, Indel_size2));

	//'same' means that Indel size ratio must be at least 80% and same type : Insertion of same size as deletion is NOT the same
	same = (sv.type_enum == sv2.type_enum && 
		ovlp > smap_Indel_Ref_SV_Overlap + ((smap_Del_Ref_SV_OverlapIncrease && (sv.type_enum == deletion || sv.type_enum == compression)) ? 0.5*(Indel_size1 + Indel_size2) : 0.0)/* WAS28 ovlp/denom > 0.8 */ &&
		IndelRatio > smap_Indel_Size_MatchRatio /* WAS28 qryratio > 0.8 */); 
      } else if( istrans ) {
	//translocations: breakpoints must be similar AND in the same direction (ie, not reciprocal) to be same
	//do nothing if at least one breakpoint is not in common: match criteria in alignmentOverlapTranslocation,
	// namely, that breakpoints are within twice translocation_overlap_buf (refcontigid checked above)
	double p11 = sv.refstart, p12 = sv.refstop, p21 = sv2.refstart, p22 = sv2.refstop;
	ol1 = overlap( p11-translocation_overlap_buf, p11+translocation_overlap_buf,
		       p21-translocation_overlap_buf, p21+translocation_overlap_buf );
	//printf("svOverlap: sv0 %4i: %10.1f - %10.1f; %10.1f - %10.1f : %i\n",
	//       p11-translocation_overlap_buf, p11+translocation_overlap_buf,
	//       p21-translocation_overlap_buf, p21+translocation_overlap_buf, ol1 );
	ol2 = overlap( p12-translocation_overlap_buf, p12+translocation_overlap_buf,
		       p22-translocation_overlap_buf, p22+translocation_overlap_buf );
	//don't forget to also require the same chromosome (above is just one in common)
	ol1 = ol1 && sv.refcontigid1 == sv2.refcontigid1;
	ol2 = ol2 && sv.refcontigid2 == sv2.refcontigid2;
	if( !ol1 && !ol2 )
	  continue;
	if( ol1 && ol2 ) { //must satisfy this to be same
	  //get directions: each mg is either right or left of each bp
	  d11 = ( entry1->refstartpos == p11 ); //first mg on first trans
	  d21 = ( usexmapentries[sv.indices[1]]->refstartpos == p12 ); //second mg on first trans
	  d12 = ( entry2->refstartpos == p21 ); //first mg on second trans
	  d22 = ( usexmapentries[sv2.indices[1]]->refstartpos == p22 ); //second mg on second trans
	  same = (d11 == d12 && d21 == d22);
	}
	//for translocations, use indel_overlap_id to prevent a pair of translocations on one GM from making
	// its partner GM translocations into unknowns (they should all be hom if no third GM overlaps):
	// if same, record it, and if not same but id matches, call it same: should be logged below
	if( same ) {
	  sv_array[indp].indel_overlap_id = entry2->qrycontigid;
	  sv_array[indr].indel_overlap_id = entry1->qrycontigid;
	} else if( sv_array[indp].indel_overlap_id == entry2->qrycontigid ) {
	  same = true;
	}
      } else { //isinv
	if( sv.pair_index == sv2.pair_index ) {
	  if(VERB >= 2){
	    printf("svOverlap: same pair: %i %i\n", indp, indr);
	    fflush(stdout);
	  }
	}
	//'dumb' algo: both breakpoints must be within X to be 'same', otherwise different.
	//But if no overlap at all, then it's neither same nor different, then they don't affect each other
	// Just use translocation_overlap_buf for X (actually that is X/2).
	double bp1_1 = min(sv.refstart, sv.refstop); //not sure if start/stop is always same, so sort again
	double bp1_2 = max(sv.refstart, sv.refstop);
	double bp2_1 = min(sv2.refstart, sv2.refstop);
	double bp2_2 = max(sv2.refstart, sv2.refstop);
	if( !overlap(bp1_1, bp1_2, bp2_1, bp2_2) )
	  continue;
	ol1 = overlapPoint( bp1_1-translocation_overlap_buf, bp1_1+translocation_overlap_buf, bp2_1 );
	ol2 = overlapPoint( bp1_2-translocation_overlap_buf, bp1_2+translocation_overlap_buf, bp2_2 );
	same = (ol1 && ol2);
      }

      //'same': set indel_overlap accordingly for both indels
      if( same ) { 
	if( sv_array[indp].indel_overlap == 0 )
	  sv_array[indp].indel_overlap = 1;
	else if( sv_array[indp].indel_overlap == 2 )
	  sv_array[indp].indel_overlap = 3;

	if( sv_array[indr].indel_overlap == 0 )
	  sv_array[indr].indel_overlap = 1;
	else if( sv_array[indr].indel_overlap == 2 )
	  sv_array[indr].indel_overlap = 3;

	//if same, both get same genotype, but must check if already assigned
	if( sv_array[indr].genotype == -1 && sv_array[indp].genotype == -1 )
	  sv_array[indr].genotype = sv_array[indp].genotype = 1;
	else if( sv_array[indr].genotype != -1 )
	  sv_array[indp].genotype = sv_array[indr].genotype;
	else if( sv_array[indp].genotype != -1 )
	  sv_array[indr].genotype = sv_array[indp].genotype;
      } else { //'different' but overlap: also change indel_overlap
	if( sv_array[indp].indel_overlap == 0 ) {
	  sv_array[indp].indel_overlap = 2;
	  sv_array[indp].indel_overlap_id = entry2->qrycontigid; //multiple 'different' on same contig are ignored: this is the first
	}
	//both 1 and 2 go to 3 now, and there's no harm in assigning 3 if it was already 3
	else if( sv_array[indp].indel_overlap_id != entry2->qrycontigid ) //already have one on this contig
	  sv_array[indp].indel_overlap = 3;

	if( sv_array[indr].indel_overlap == 0 ) {
	  sv_array[indr].indel_overlap = 2;
	  sv_array[indr].indel_overlap_id = entry1->qrycontigid;
	}
	else if( sv_array[indr].indel_overlap_id != entry1->qrycontigid )
	  sv_array[indr].indel_overlap = 3;

	//here, genotype must be different: similar to above, but instead of assigning same, assign opposite
	if( sv_array[indr].genotype == -1 && sv_array[indp].genotype == -1 ) {
	  sv_array[indp].genotype = 1;
	  sv_array[indr].genotype = 2;
	}
	else if( sv_array[indr].genotype != -1 ) 
	  sv_array[indp].genotype = oppositeGenotype( sv_array[indr].genotype );
	else if( sv_array[indp].genotype != -1 )
	  sv_array[indr].genotype = oppositeGenotype( sv_array[indp].genotype );
      }

      //genotype group: if either is set, assign to other, otherwise use smap_genotype_ngroup
      if( sv_array[indp].genotypeGroup != -1 ) //-1 is default -- must check the sv_array because the sv object is not modified
	sv_array[indr].genotypeGroup = sv.genotypeGroup;
      else if( sv_array[indr].genotypeGroup != -1 ) //-1 is default -- same as previous
	sv_array[indp].genotypeGroup = sv2.genotypeGroup;
      else {
	smap_genotype_ngroup++;
	sv_array[indr].genotypeGroup = smap_genotype_ngroup;
	sv_array[indp].genotypeGroup = smap_genotype_ngroup;
      }

      const char* t = ( isindel ? "ind" : "tra" );
      t = ( isinv ? "inv" : t ); //note this logic breaks if more types are added

      if(VERB>=2){
	printf("svOverlap: %s: sv1 %4i: %2lld %10.1f - %2lld %10.1f; sv2 %4i: %2lld %10.1f - %2lld %10.1f; ", t, indp+1, sv.refcontigid1, sv.refstart, sv.refcontigid2, sv.refstop, indr+1, sv2.refcontigid1, sv2.refstart, sv2.refcontigid2, sv2.refstop);
	if( isindel )
	  printf("ref ovlp= %7.1f Indel Size ratio: %.3f;", ovlp, IndelRatio);
	else
	  printf("bp ovlp: %i %i, dir: %i %i %i %i, %i;", ol1, ol2, d11, d12, d21, d22, same);
	printf(" ids %4lld %4lld; gen: %i %i, %i %i\n", sv_array[indp].indel_overlap_id, sv_array[indr].indel_overlap_id, sv_array[indp].genotypeGroup, sv_array[indr].genotypeGroup, sv_array[indp].genotype, sv_array[indr].genotype);
	fflush(stdout);
      }
    } //end inner loop on svs

    //if( sv_array[indp].genotypeGroup == -1 ) { //assign group to homozygous (no overlap) indels -- if not in a group, no reason to use this field
    //  smap_genotype_ngroup++;
    //  sv_array[indp].genotypeGroup = smap_genotype_ngroup;
    //}
    //but make all genotypes of indels either 1 or 2: assign 1
    if( sv_array[indp].genotype == -1 )
      sv_array[indp].genotype = 1;
  } //end outter loop on svs

  //summarize results of alignmentOverlap and svOverlap (debug)
  // for(int indp=0; indp<numsv; indp++) {
  //   printf("Summary: %2i (%9.1f - %9.1f) alignment_overlap: %i svOverlap: %i\n", indp+1, sv_array[indp].refstart, sv_array[indp].refstop, sv_array[indp].alignment_overlap, sv_array[indp].indel_overlap);    
  //}
} //end fn svOverlap

void structuralVariation :: setIndelConfidence(xmapEntry **usexmapentries, int verbose = 1)
{
  if( algo_enum != sv_ps ) //-indel case is handeled by Thomas in output_xmap.cpp
    return;

  //dummy checks: size is in kb, while positions are in bp
  if(DEBUG && !(fabs(querysize - (querystop-querystart)*1e-3) < 1e-10 && fabs(refsize - (refstop-refstart)*1e-3) < 1e-10)){
    xmapEntry *entry1 = ((algo_enum == sv_ps) ? usexmapentries : xmapentries)[indices[0]];
    xmapEntry *entry2 = ((algo_enum == sv_ps) ? usexmapentries : xmapentries)[indices[1]];
    
    printf("setIndelConfidence:this=%p: sv_ps=%d,type_enum=%d: qry= %0.6f..%0.6f(siz=%0.6f), ref= %0.6f .. %0.6f(siz=%0.6f)\n", 
	   this, algo_enum == sv_ps ? 1 : 0, type_enum,querystart*1e-3,querystop*1e-3,querysize,refstart*1e-3,refstop*1e-3,refsize);
    printf("\t refcontigid1=%lld,refcontigid2=%lld,qryid=%lld,%lld: xmapid1=%d,xmapid2=%d\n",refcontigid1, refcontigid2,entry1->qrycontigid, entry2->qrycontigid, xmapid1, xmapid2);
    fflush(stdout);

    assert(fabs(querysize - (querystop-querystart)*1e-3) < 1e-6);
    assert(fabs(refsize - (refstop-refstart)*1e-3) < 1e-6);
  }

  //3rd and 4th args are number of intervals, ie, one more than number of misaligned sites
  if(DEBUG) assert(querystopidx >= querystartidx);// NOTE : if this is not true, use absolute value of m
  int m = querystopidx-querystartidx;
  int n = refstopidx-refstartidx;

  double outscore = 0.0;
  int bc       = 0;
  double iscore   = 0.0;
  double IndelConfidence = 0.0;

  if(PairSplit) {
    // declared in NGentigPairScore.h
    PFLOAT Bias=0, Pen=0; //Pen is passed by reference--result is score
    extern void SintDetailPW(PFLOAT X, PFLOAT Y, int m, int n, int J, int I, PFLOAT &Bias, PFLOAT &Pen);

    // if the start and end label is the same, this fn cannot handle this case. use single hypothetical small interval of size 0.1 kb
    SintDetailPW(max(0.1,querysize), max(0.1,refsize), max(1,m), max(1,n), 0, 0, Bias, Pen); 

    //note in RefAlign.cpp: if Psplit == 0, it is set to a default, if Poutlier not assigned, also given default; if Poutlier < Psplit, Poutlier set == to Psplit

    outscore = Pen + Bias;
    bc       = (outlierBC ? max(1,m)*max(1,n) : 1); //similar to above: cannot have 0
    iscore   = log(Psplit/bc); 

    // Confidence is degree to which actual score is worse than minimum outlier score (+ Psplit threshold value)
    IndelConfidence = (iscore - outscore - log(Psplit))/log(10.0);
  } else { // RefSplit
    // declared in RGentigRefScore.h
    extern void SintDetailRA(double X, double Y, int m, int n, int J, int I, int K /* right */, int T /* left */, FLOAT *Yref, double &Bias, double &Pen, double &Gauss, double &PenSm, int verb);
    extern RFLOAT Frate;

    Calign *p1 = usexmapentries[indices[0]]->align;
    Calign *p2 = usexmapentries[indices[1]]->align;
    if(DEBUG) assert(p1->mapid1 == p2->mapid1);
    if(DEBUG) assert(p1->mapid2 == p2->mapid2);
    if(DEBUG) assert(p1->orientation == p2->orientation);

    Cmap *Ymap = refmap[p1->mapid1];
    if(DEBUG) assert(Ymap->origmap == 0);
    double *Y1 = Ymap->site[0];// reference
    int N = Ymap->numsite[0];// number of labels in reference
    Cmap *Xmap = Gmap[p1->mapid2];
    if(DEBUG) assert(Xmap->origmap == 0);
    double *X1 = Xmap->site[0];// query
    int M = Xmap->numsite[0];// number of labels in query

    int U1 = p1->numpairs;
    int U2 = p2->numpairs;
    int L1 = p1->sites1[0];
    int LK1 = p1->sitesK1[0];
    int LQ1 = p1->sites2[0];
    int R1 = p1->sites1[U1-1];//rightmost aligned site in reference for matchgroup1
    int RK1 = p1->sitesK1[U1-1];
    int RQ1 = p1->sites2[U1-1];//rightmost aligned site on query for matchgroup1

    int L2 = p2->sites1[0]; //leftmost aligned site on reference for matchgroup2
    int LK2 = p2->sitesK1[0];
    int LQ2 = p2->sites2[0];// leftmost aligned site on query for matchgroup2
    int R2 = p2->sites1[U2-1];
    int RK2 = p2->sitesK1[U2-1];
    int RQ2 = p2->sites2[U2-1];
    //extern Cmap **YYmap,**XXmap; //moved above
      
    if(0 && DEBUG && !(R1 <= L2 || R2 <= L1)){
      if(VERB){
	printf("setIndelConfidence:xmapid=%d,%d:mapid1=%d(id=%lld),mapid2=%d,%d(id=%lld,%lld),or=%d,%d:U1=%d,U2=%d,L1=%d,R1=%d,L2=%d,R2=%d,align->score=%0.6f,%0.6f,logPV=%0.2f,%0.2f\n",
	       indices[0],indices[1],p1->mapid1,YYmap[p1->mapid1]->id,p1->mapid2,p2->mapid2,XXmap[p1->mapid2]->id,XXmap[p2->mapid2]->id,p1->orientation,p2->orientation,
	       U1,U2,L1,R1,L2,R2,p1->score,p2->score,p1->logPV,p2->logPV);
	fflush(stdout);
      }
      assert(R1 <= L2 || R2 <= L1);
    }

    double Bias=0, Pen=0, Gauss=0, PenSm=0, deltaX = 0.0;
    if(R1 <= L2-LK2 && RQ1 <= LQ2){/* interval (L1..R1) is left of interval (L2..R2) */
      int T = p1->sitesK1[U1-1];
      int K = p2->sitesK1[0];
      double delY = Yc(Y1,L2,K) - Yc(Y1,R1,T);
	
      if(DEBUG) assert(R1-T > 1);
      //if(DEBUG && L2-K-R1 <= 0) assert(m > 0); //was this: suppress assertion because can be triggered in rare cases by RefSplit
      if(DEBUG && L2-K-R1 <= 0 && m <= 0) {
	if(VERB>=2){
	  printf("Warning: setIndelConfidence: L2=%i K=%i R1=%i m=%i; xmapids=%i %i\n", L2, K, R1, m, usexmapentries[indices[0]]->xmapid, usexmapentries[indices[1]]->xmapid);
	  fflush(stdout);
	}
	//assert(m > 0);
      } /* WAS26 else */
      deltaX = max(0.1,querysize);
      SintDetailRA(deltaX, max(0.1,delY), max(1,m), max(1,L2-K-R1), p2->sites2[0], L2, K, T, Y1, Bias, Pen, Gauss, PenSm, 0);

      if(DEBUG && Bias + Pen + Gauss > 0.0){
	if(VERB>=2){
	  printf("WARNING:setIndelConfidence:xmapid=%d,%d:mapid1=%d(id=%lld),mapid2=%d,%d(id=%lld,%lld),or=%d,%d:U1=%d,U2=%d,L1=%d,R1=%d,L2=%d,R2=%d,align->score=%0.6f,%0.6f,logPV=%0.2f,%0.2f\n",
		 indices[0],indices[1],p1->mapid1,YYmap[p1->mapid1]->id,p1->mapid2,p2->mapid2,XXmap[p1->mapid2]->id,XXmap[p2->mapid2]->id,p1->orientation,p2->orientation,
		 U1,U2,L1,R1,L2,R2,p1->score,p2->score,p1->logPV,p2->logPV);
	  printf("    T=%d,K=%d:querysize=%0.4f,refsize=%0.4f,delY=%0.4f,m=%d,n=%d,J1=%d,J2=%d:Bias=%0.6f,Pen=%0.6f,Gauss=%0.6f,PenSm=%0.6f\n",
		 T,K,querysize,refsize,delY,m,n,p1->sites2[0],p2->sites2[0],Bias,Pen,Gauss,PenSm);
	  fflush(stdout);
	}
	// WAS	  assert(Bias + Pen + Gauss <= 0.0);
      }

    } else if(R2 <= L1-LK1 && RQ2 <= LQ1) {/* interval (L1..R1) is right of interval (L2..R2) */
      int T = p2->sitesK1[U2-1];
      int K = p1->sitesK1[0];
      double delY = Yc(Y1,L1,K) - Yc(Y1,R2,T);

      if(DEBUG) assert(R2-T > 1);
      if(DEBUG && L1-K-R2 <= 0) assert(m > 0);

      deltaX = max(0.1,querysize);
      SintDetailRA(deltaX, max(0.1,delY), max(1,m), max(1,L1-K-R2), p1->sites2[0], L1, K, T, Y1, Bias, Pen, Gauss, PenSm, 0);
	
      if(DEBUG && Bias + Pen + Gauss > 0.0){
	if(VERB>=2){
	  printf("WARNING:setIndelConfidence:xmapid=%d,%d:mapid1=%d(id=%lld),mapid2=%d,%d(id=%lld,%lld),or=%d,%d:U1=%d,U2=%d,L1=%d,R1=%d,L2=%d,R2=%d,align->score=%0.6f,%0.6f,logPV=%0.2f,%0.2f\n",
		 indices[0],indices[1],p1->mapid1,YYmap[p1->mapid1]->id,p1->mapid2,p2->mapid2,XXmap[p1->mapid2]->id,XXmap[p2->mapid2]->id,p1->orientation,p2->orientation,
		 U1,U2,L1,R1,L2,R2,p1->score,p2->score,p1->logPV,p2->logPV);
	  printf("    T=%d,K=%d:querysize=%0.4f,refsize=%0.4f,delY=%0.4f,m=%d,n=%d,J1=%d,J2=%d:Bias=%0.6f,Pen=%0.6f,Gauss=%0.6f,PenSm=%0.6f\n",
		 T,K,querysize,refsize,delY,m,n,p1->sites2[0],p2->sites2[0],Bias,Pen,Gauss,PenSm);
	  fflush(stdout);
	}
	// WAS	  assert(Bias + Pen + Gauss <= 0.0);
      }
    } else if(L1 <= L2-LK2 && LQ1 <= LQ2 && (L2-LK2 - L1 + LQ2 - LQ1 > 0) &&
	      R1 <= R2-RK2 && RQ1 <= RQ2 && (R2-RK2 - R1 + RQ2 - RQ1 > 0)){/* overlapped matchgroups : matchgroup1 left of matchgroup2 */
      /* right end of Indel is (L2,LK2,LQ2) */
      /* find left end of Indel by scaning matchgroup1 */
      int t = U1-1;/* index into matchgroup1 */
      for(;t >= 0; t--){
	int r = p1->sites1[t];
	int s = p1->sitesK1[t];// NEW
	int q = p1->sites2[t];
	if(r <= L2 && r-s <= L2-LK2 && q <= LQ2)
	  break;
      }

      if(DEBUG) assert(t >= 0);
      int R = p1->sites1[t];
      int K = p1->sitesK1[t];
      int Q = p1->sites2[t];
      /* now left end of Indel is (R,K,Q) */
	
      int m = LQ2 - Q;
      int n = L2-LK2 - R;
      double rsize = Yc(Y1, L2, LK2) - Yc(Y1, R, K);
      int Qend   = p2->orientation ? M+1-LQ2 : LQ2;
      int Qstart = p1->orientation ? M+1-Q : Q;
      double qsize = fabs(X1[Qend] - X1[Qstart]);

      if(verbose){
	printf("overlap indel1: xid=%d,%d, refid=%lld, qryid=%lld, fwd=%d: qsize=%.4f, querysize=%.4f; rsize=%.4f, refsize=%.4f, m=%d, n=%d; refstartidx=%i refendidx=%i; qrystartidx=%i qryendidx=%i\n", 
	       xmapid1,xmapid2,refcontigid1, xmapentries[indices[0]]->qrycontigid, xmapentries[indices[0]]->orientforward ? 1 : 0, 
	       qsize, querysize, rsize, refsize, m, n, R-K, L2, Qstart, Qend);
	fflush(stdout);
      }
      deltaX = max(0.1,qsize);
      SintDetailRA(deltaX, max(0.1,rsize), max(1,m), max(1,n), LQ2, L2, LK2, K, Y1, Bias, Pen, Gauss, PenSm, 0);

      if(DEBUG && !(Bias + Pen + Gauss <= 0.0)){
	printf("\tBias= %0.6f, Pen= %0.6f, Gauss=%0.6f,PenSm=%0.6f\n",Bias, Pen, Gauss, PenSm);
	fflush(stdout);
	//	assert(Bias + Pen + Gauss <= 0.0);// this can fail when gap size is 0
      }

      //note: Qend and Qstart depend on orientation: just take min/max for assign to data members
    } else if(L1-LK1 >= L2 && LQ1 >= LQ2 && (L2-L1+LK1 + LQ2-LQ1 < 0) &&
      R1-RK1 >= R2 && RQ1 >= RQ2 && (R2-R1+RK1 + RQ2-RQ1 < 0)){/* overlapped matchgroup : matchgroup2 left of matchgroup1 */
      /* right end of Indel is (L1,LK1,LQ1) */
      int t = U2-1;/* index into matchgroup2 */
      for(;t >= 0; t--){
	int r = p2->sites1[t];
	int s = p2->sitesK1[t];// NEW
	int q = p2->sites2[t];
	if(r <= L1 && r-s <= L1-LK1 && q <= LQ1)
	  break;
      }
      if(DEBUG) assert(t >= 0);
      int R = p2->sites1[t];
      int K = p2->sitesK1[t];
      int Q = p2->sites2[t];
      /* now left end of Indel is (R,K,Q) */
	
      int m = LQ1 - Q;
      int n = L1-LK1 - R;
      double rsize = Yc(Y1, L1, LK1) - Yc(Y1, R, K);
      int Qend   = p2->orientation ? M+1-LQ1 : LQ1;
      int Qstart = p1->orientation ? M+1-Q : Q;
      double qsize = fabs(X1[Qend] - X1[Qstart]);

      deltaX = max(0.1,qsize);
      SintDetailRA(deltaX, max(0.1,rsize), max(1,m), max(1,n), LQ1, L1, LK1, K, Y1, Bias, Pen, Gauss, PenSm, 0);
      if(verbose)
	printf("overlap indel2: qsize=%.4f, qrysize=%.4f; rsize=%.4f, refsize=%.4f, m=%d, n=%d; refstartidx=%i refendidx=%i; qrystartidx=%i qryendidx=%i\n", qsize, querysize, rsize, refsize, m, n, R-K, L1, Qstart, Qend);
      //      if(DEBUG) assert(Bias + Pen + Gauss <= 0.0); // this can fail when gap size is 0

    } else {/* should never happen */
      printf("ERROR : unhandled case (crossed matchgroups): xmapid=%d,%d:mapid1=%d(id=%lld),N=%d,mapid2=%d(id=%lld),M=%d,or=%d,%d:U1=%d,U2=%d:L1=%d,R1=%d,L2=%d,R2=%d,LQ1=%d,RQ1=%d,LQ2=%d,RQ2=%d,LK1=%d,RK1=%d,LK2=%d,RK2=%d\n",
	     usexmapentries[indices[0]]->xmapid,usexmapentries[indices[1]]->xmapid,
	     p1->mapid1,YYmap[p1->mapid1]->id,N,p1->mapid2,XXmap[p1->mapid2]->id,M,
	     p1->orientation,p2->orientation,U1,U2,L1,R1,L2,R2,LQ1,RQ1,LQ2,RQ2,LK1,RK1,LK2,RK2);
      fflush(stdout);
      exit(1);
    }

    outscore = smap_IndelConfBySize ? Gauss + (FRATE_FIX1 ? deltaX * Frate : 0.0) : Bias + Pen + Gauss; // NOTE : excludes PenSm which is the misresolved penalty for the right end of the interval
    //      bc       = (outlierBC ? max(1,m)*max(1,L2-K-R1) : 1); //similar to above: cannot have 0
    //      iscore = log(Psplit/bc);
      
    // Confidence is negative raw score (scaled to log10) */
    IndelConfidence =  -outscore/log(10.0);
    // WAS      if(DEBUG) assert(IndelConfidence > 0.0);
  } // RefSplit

    //algo_enum checked above
  double left_conf  = usexmapentries[indices[0]]->confidence;
  double right_conf = usexmapentries[indices[1]]->confidence;
  if(VERB>=2 && verbose){
    printf("querysize=%.2f, refsize=%.2f, m=%i, n=%i, outscore=%.2f; left_conf=%.2f, right_conf=%.2f\n", querysize, refsize, m, n, outscore, left_conf, right_conf);
    fflush(stdout); //note: if disable above printf, then use if below instead
  }

  //double finalconf = min(min(IndelConfidence, left_conf), right_conf); //this is the approximate version
  if( IndelConfidence > 0 ) //confidence < 0 is not an outlier: set to 0 and filter in writeSVToSmap
    confidence = makeConfidence(IndelConfidence, left_conf, right_conf);
  else
    confidence = 0;

  //assign to data members
  confidence_left   = left_conf;
  confidence_right  = right_conf;
  confidence_center = (IndelConfidence > 0 ? IndelConfidence : 0);

  //if(PairSplit && PairsplitMerge && !(outscore <= iscore)){
  if(PairSplit && !(outscore <= iscore)){
    //this assertion is only valid if pairsplitMerge is enabled; otherwise, just assign 0 becasue it's a bad split
    //however, I don't want to crash, so only enable assertion if DEBUG 2, otherwise make it a warning
    if( DEBUG >= 2 && PairsplitMerge )
      assert(outscore <= iscore); //otherwise shouldn't have split
    else {
      printf("WARNING: outscore > iscore; setting confidence to 0 (try using -pairsplitMerge)\n");
      fflush(stdout);
      confidence = 0;
    }
  }

} //end setIndelConfidence


//compare argument for qsort
static int compareFloat(const void *p1, const void *p2) {
  //return (*(float*)p1 < *(float*)p2) ? 1 : (*(float*)p1 > *(float*)p2) ? -1 : 0; //lovely as this line is, below is even moreso
  return (*(const float*)p1 < *(const float*)p2) - (*(const float*)p1 > *(const float*)p2);
}


//needed for Luna's predictors: median occurence of all labels in the sample
//since Thomas stores occurrence in float, use float here, then convert to double when used in predictors
//also need mean ChimNorm: make args occ: true for occurrence, false for ChimNorm, and median (mean if false)
float getQueryMapStat(bool occ=true, bool median=true) {
  int nsites = 0; //first loop on all maps to get total nsites

  for(int i=0; i < numX; i++) { //numX is length of XXmap
    nsites += XXmap[i]->numsite[0]+1; //add 1 bc Luna's calculation included end of map
  }
  float* occs = new float[nsites];
  int gi = 0; //index in occs
  for(int i=0; i < numX; i++) { //numX is length of XXmap
    for(int j=1; j <= XXmap[i]->numsite[0]; j++) { //point to XXmap[i]? not necessary
      occs[gi] = occ ?
	XXmap[i]->sitecnt[0][j] : //sitecnt is the occurrence
	XXmap[i]->ChimNorm[0][j]; //ChimNorm if not occurrence
      gi++;
    }
    occs[gi] = 0; //add end of map
    gi++;
    assert(gi <= nsites); //dummy check
  }
  assert(gi == nsites); //make sure all the entries in occs are filled, otherwise can be used uninitialized 
  if( !median ) { //didn't actually need array of values, but simpler this way
    float sum = 0;
    for(int i=0; i < nsites; i++)
      sum += occs[i];
    return sum/nsites;
  }
  qsort(occs, nsites, sizeof(float), compareFloat);
  float med; //don't return yet bc need to delete occs
  if( nsites % 2 == 0 )
    med = (occs[nsites/2] + occs[nsites/2-1])/2.;
  else
    med = occs[nsites/2];
  delete [] occs;
  return med;
}


//predictors which are common to both translocations and inversions:
//overlapsize, mg, mg12 : should be exactly the same (if 11 == 12)
//overlapsize12 : inversions only, but in same loop as mg/12
//gapsize : translocations only but same conditional as overlapsize
int getCommonPredictors(structuralVariation* sv, xmapEntry** usexmapentries, int lbuf, double medocc, double& overlapsize, double& overlapsize12, double& nmg, double& nmg11, double& gapsize, double& gapNlab, int& startlab, int& stoplab, double& occ_mean, double& occ_max, double& fragile_max, double& chimq_max, double& occ_mean_bp, double& occ_max_bp, double& chimn_max) {
  //const int lbuf = 11; //label buffer -- promote to argument
  const int endbuf = 7; //exclude first/last labels on GM

  xmapEntry* xe1 = (sv->algo_enum == sv_ps) ? usexmapentries[sv->indices[0]] : xmapentries[sv->indices[0]];
  xmapEntry* xe2 = (sv->algo_enum == sv_ps) ? usexmapentries[sv->indices[1]] : xmapentries[sv->indices[1]];

  if(DEBUG) assert(xe1->qrycontigid == xe2->qrycontigid);

  Calign *p1 = xe1->align;
  Calign *p2 = xe2->align;
  if(DEBUG) assert(p1->mapid2 == p2->mapid2); //mapid1 is reference: these are not the same for translocations
  Cmap *Xmap = Gmap[p1->mapid2]; //query map
  if(DEBUG) assert(Xmap->origmap == 0);
  if( Xmap->FragileEndL[0] == 0 || Xmap->ChimQuality[0] == 0 ) { //the arrays are not allocated
    printf("Missing fragile/chimQuality in map id %lld: disabling confidence\n", xe1->qrycontigid);
    fflush(stdout);
    return 1;
  }
  
  int M = Xmap->numsite[0];// number of labels in query
  double meanchimn = 0;
  for(int i=1; i <= M; i++)
    meanchimn += Xmap->ChimNorm[0][i];
  meanchimn /= M;
  gapNlab = double(abs(sv->querystartidx - sv->querystopidx));
  //gapsize is distance between breakpoints (bp) if no overlap and 0 if overlap
  //overlapsize is opposite: overlapsize (bp) if overlap otherwise 0
  if( overlap(min(xe1->qrystartpos, xe1->qryendpos), max(xe1->qrystartpos, xe1->qryendpos),
	      min(xe2->qrystartpos, xe2->qryendpos), max(xe2->qrystartpos, xe2->qryendpos)) )
    overlapsize = fabs(sv->querystart - sv->querystop);    
  else
    gapsize = fabs(sv->querystart - sv->querystop);
  //nmg is total num matchgroups on this query; nmg11 is num matchgroups which start/end within
  // 11 labels to left of left breakpoint and 11 labels to right of right breakpoint
  //these are really ints, but since they're going to end up in a matrix, they must be double
  startlab = min(sv->querystartidx, sv->querystopidx)-lbuf-1;
  stoplab  = max(sv->querystartidx, sv->querystopidx)+lbuf+1;
  for(register int i=0; i < numxmapentries; i++){
    if( xmapentries[i]->qrycontigid != xe1->qrycontigid ||
	xmapentries[i]->confidence < LogPvThreshold )
      continue;
    nmg++;
    int mgstartlab = min(xmapentries[i]->qrystartidx, xmapentries[i]->qryendidx);
    int mgstoplab  = max(xmapentries[i]->qrystartidx, xmapentries[i]->qryendidx);
    //overlap use </> NOT <=/>=, so add/subtract lbuf+1 to include stop/start on next label from bp
    if( overlap(startlab, stoplab, mgstartlab, mgstoplab) ) {
      nmg11++;
      //Cmap.h: FLOAT *site[MAXCOLOR];/**< site[c=0..colors-1][i=0..numsite[c]+1] : the location of site i in color c (in kb from the left end site[0]) **/
      int i1 = max(startlab, mgstartlab);
      int i2 = min(stoplab , mgstoplab);
      overlapsize12 += (Xmap->site[0][i2] - Xmap->site[0][i1])*1e3; //bp
      if(VERB>=2){
	printf("  gCP: gm %5lld, mcn %f, ovlp xid %4i, conf %6.2f, ovlp lab: %2i %2i, pos: %8.1f %8.1f, diff %8.1f, ovlp %8.1f\n",
	       xe1->qrycontigid, meanchimn, xmapentries[i]->xmapid, xmapentries[i]->confidence, i1, i2, Xmap->site[0][i2]*1e3, Xmap->site[0][i1]*1e3, (Xmap->site[0][i2] - Xmap->site[0][i1])*1e3, overlapsize12); //debug
	fflush(stdout);
      }
    }
  }

  startlab = min(sv->querystartidx, sv->querystopidx);
  stoplab  = max(sv->querystartidx, sv->querystopidx);
  if(DEBUG) assert(1 <= startlab && startlab <= stoplab && stoplab <= M);

  float oc1 = Xmap->sitecnt[0][startlab], oc2 = Xmap->sitecnt[0][stoplab]; //sitecnt is the occurrence
  occ_mean_bp = (oc1 + oc2)/2./medocc;
  occ_max_bp  = max(oc1, oc2); //need to normalize by medocc
  //last set are range from left bp - lbuf to right bp + rbuf, but excluding first/last endbuf labels on GM
  //indices are 1 to N, so exclude endbuf means start at endbuf+1. 
  int i1 = /* WAS370 max(endbuf+1, startlab - lbuf) */ max(min(startlab,endbuf+1), startlab - lbuf);
  int i2 = /* WAS370 min(M-endbuf, stoplab + lbuf)*/ min(max(stoplab,M-endbuf), stoplab + lbuf);
  for(int i=i1; i <= i2; i++) {
    occ_mean   += Xmap->sitecnt[0][i];
    if(DEBUG) assert(isfinite(occ_mean));
    occ_max     = max(occ_max    , (double)Xmap->sitecnt[0][i]);
    fragile_max = max(fragile_max, (double)max(Xmap->FragileEndL[0][i], Xmap->FragileEndR[0][i]));
    chimq_max   = max(chimq_max  , (double)(100.-Xmap->ChimQuality[0][i])); //convert to same scale as FragileEndL/R
    chimn_max   = max(chimn_max  , (double)100.*(1.-Xmap->ChimNorm[0][i]/meanchimn));
  }
  occ_mean /= i2-i1+1; //both i1 and i2 are included, so +1
  if(DEBUG && !isfinite(occ_mean)){
    printf("gCP: occ_mean=%0.3f, i1=%d,i2=%d: lbuf=%d,endbuf=%d,M=%d,startlab=%d,stoplab=%d\n",occ_mean,i1,i2,lbuf,endbuf,M,startlab,stoplab);
    fflush(stdout);
    assert(isfinite(occ_mean));
  }

  return 0;
} //end getCommonPredictors


//get all Luna's predictors for translocation confidence for this sv
void getTransPredictors(structuralVariation* sv, xmapEntry **usexmapentries, double medocc, double* pred, int npred) {
  for( int i=0; i < npred; i++ )
    pred[0] = 0.0; //null all pred elements in case return early

  if(medocc <= 0.0){
    printf("Median occurance is zero in map id %lld: disabling confidence\n", usexmapentries[sv->indices[0]]->qrycontigid);
    fflush(stdout);
    return;
  }

  //gapNlab is number of labels between two breakpoints (regardless of overlap)
  //double gapNlab = double(abs(sv->querystartidx - sv->querystopidx));
  xmapEntry* xe1 = usexmapentries[sv->indices[0]];
  //xmapEntry* xe2 = usexmapentries[sv->indices[1]];
  double gapNlab = 0, gapsize = 0, overlapsize = 0, overlapsize12 = 0;
  double nmg = 0, nmg11 = 0;
  //int startlab = min(sv->querystartidx, sv->querystopidx);
  //int stoplab  = max(sv->querystartidx, sv->querystopidx);
  int startlab = 0, stoplab = 0;

  //label-based (quality score) predictors:
  double occ_mean = 0, occ_max = 0, fragile_max = 0, chimq_max = 0, chimn_max = 0; //last predictors
  double occ_mean_bp = 0, occ_max_bp = 0;
  //third arg is int lbuf
  int err = getCommonPredictors(sv, usexmapentries, 11, medocc, overlapsize, overlapsize12, nmg, nmg11, gapsize, gapNlab, startlab, stoplab, occ_mean, occ_max, fragile_max, chimq_max, occ_mean_bp, occ_max_bp, chimn_max);
  if(err)
    return;

  if(0){
    //printf("getTransPredictors: qry id %5lld: sv qry idx: %3i %3i, occ_mean_bp: %.3f (%.3f %.3f), occ_max_bp: %.3f (%.3f), nmg/11: %.0f %.0f, gapnlab/size: %.0f %.3f, ovlpsz: %.3f;\n  range %3i - %3i (max %3i): occ_mean: %.3f (%.3f), occ_max: %.3f (%.3f), fragile_max: %.3f, chimq_max: %.3f\n", //old
    //xe1->qrycontigid, startlab, stoplab, occ_mean_bp, oc1, oc2, occ_max_bp, occ_max_bp/medocc, nmg, nmg11, gapNlab, gapsize, overlapsize, //old
    printf("getTransPredictors: qry id %5lld: sv qry idx: %3i %3i, medocc=%0.3f,occ_mean_bp: %.3f, occ_max_bp: %.3f (%.3f), nmg/11: %.0f %.0f, gapnlab/size: %.0f %.3f, ovlpsz: %.3f;\n  occ_mean: %.3f (%.3f), occ_max: %.3f (%.3f), fragile_max: %.3f, chimq_max: %.3f\n",
	   xe1->qrycontigid, startlab, stoplab, medocc,occ_mean_bp, occ_max_bp, occ_max_bp/medocc, nmg, nmg11, gapNlab, gapsize, overlapsize,
	   occ_mean, occ_mean/medocc, occ_max, occ_max/medocc, fragile_max, chimq_max);
    fflush(stdout);
  }

  //need in the order Luna specified; note: occ_mean, occ_max, and occ_max_bp still need to be normalized by medocc
  pred[ 0] = occ_mean/medocc;
  pred[ 1] = occ_max/medocc;
  pred[ 2] = fragile_max;
  pred[ 3] = occ_mean_bp;
  pred[ 4] = occ_max_bp/medocc;
  pred[ 5] = chimq_max;
  pred[ 6] = nmg11;
  pred[ 7] = gapNlab;
  pred[ 8] = overlapsize;
  pred[ 9] = nmg;
  pred[10] = gapsize;

  if(DEBUG){
    for(int t = 0; t <= 10; t++){
      if(!isfinite(pred[t])){
	printf("GetTransPredictors: pred[t=%d]= %0.3e\n",t,pred[t]);
	printf("getTransPredictors: qry id %5lld: sv qry idx: %3i %3i, medocc=%0.3f,occ_mean_bp: %.3f, occ_max_bp: %.3f (%.3f), nmg/11: %.0f %.0f, gapnlab/size: %.0f %.3f, ovlpsz: %.3f;\n  occ_mean: %.3f (%.3f), occ_max: %.3f (%.3f), fragile_max: %.3f, chimq_max: %.3f\n",
	       xe1->qrycontigid, startlab, stoplab, medocc,occ_mean_bp, occ_max_bp, occ_max_bp/medocc, nmg, nmg11, gapNlab, gapsize, overlapsize,
	       occ_mean, occ_mean/medocc, occ_max, occ_max/medocc, fragile_max, chimq_max);
	fflush(stdout);
	assert(isfinite(pred[t]));
      }
    }
  }
} //end getTransPredictors

//get all Luna's predictors for inversion confidence for this sv
void getInvPredictors(structuralVariation* sv, xmapEntry **usexmapentries, double medocc, double* pred, int npred) {
  //const int endbuf = 7; //exclude first/last labels on GM
  for( int i=0; i<npred; i++ )
    pred[0] = 0; //null all pred elements in case return early
  //don't forget to pick the right xmapEntry array, and get inverted
  //xmapEntry* xe1 = usexmapentries[sv->indices[0]];
  //xmapEntry* xe2 = usexmapentries[sv->indices[1]];
  int i1 = (sv->link_index == sv->indices[0] ? sv->indices[1] : sv->indices[0]); //the non-inverted mg
  int i2 = (sv->link_index == sv->indices[0] ? sv->indices[0] : sv->indices[1]); //the inverted mg
  xmapEntry *xes = (sv->algo_enum == sv_ps) ? usexmapentries[i1] : xmapentries[i1]; //non-inverted (straight)
  xmapEntry *xei = (sv->algo_enum == sv_ps) ? usexmapentries[i2] : xmapentries[i2]; //inverted

  Calign *p1 = xei->align;
  //Calign *p2 = xes->align;
  //if(DEBUG) assert(p1->mapid2 == p2->mapid2); //mapid1 is reference: these are not the same for translocations
  if(DEBUG) assert(p1->mapid2 == xes->align->mapid2);
  Cmap *Xmap = Gmap[p1->mapid2]; //query map
  if(DEBUG) assert(Xmap->origmap == 0);
  if( Xmap->FragileEndL[0] == 0 || Xmap->ChimQuality[0] == 0 ) { //the arrays are not allocated
    printf("Missing fragile/chimQuality in map id %lld: disabling inversion confidence\n", xei->qrycontigid);
    fflush(stdout);
    return;
  }

  double gapNlab = 0, gapsize = 0, overlapsize = 0, overlapsize12 = 0;
  double nmg = 0, nmg11 = 0;
  //int startlab = min(sv->querystartidx, sv->querystopidx);
  //int stoplab  = max(sv->querystartidx, sv->querystopidx);
  int startlab = 0, stoplab = 0;
  double occ_mean = 0, occ_max = 0, fragile_max = 0, chimq_max = 0, chimn_max = 0; //last predictors
  double occ_mean_bp = 0, occ_max_bp = 0;
  //third arg is int lbuf
  getCommonPredictors(sv, usexmapentries, 10, medocc, overlapsize, overlapsize12, nmg, nmg11, gapsize, gapNlab, startlab, stoplab, occ_mean, occ_max, fragile_max, chimq_max, occ_mean_bp, occ_max_bp, chimn_max);
  //the correct condition to return is ALL the predictors are 0, however, this one should never be 0, so assume that if it is, then all the rest are, too
  if( nmg == 0 )
    return;

  //for all predictors in bp, add 1

  double minRawConf = min(sv->confidence_left, sv->confidence_right); //!= MG conf for inversion_small

  //gap predictors:
  //double gapSizeRef = fabs(sv->refstart - sv->refstop);
  //double gapNlabRef = double(abs(sv->refstartidx - sv->refstopidx));
  double gapSizeQry = fabs(sv->querystart - sv->querystop);
  double gapNlabQry = double(fabs(sv->querystartidx - sv->querystopidx))+1;
  /*
  bool refoverlap = overlap(min(xes->refstartpos, xes->refendpos), max(xes->refstartpos, xes->refendpos),
			    min(xei->refstartpos, xei->refendpos), max(xei->refstartpos, xei->refendpos));
  double invSizeRef = fabs(xei->refstartpos - xei->refendpos);
  double invNlabRef = double(abs(xei->refstartidx - xei->refendidx));
  if( !refoverlap ) { //if there is no overlap, then add gap (ie, breakpoint) size
    invSizeRef += gapSizeRef;
    invNlabRef += gapNlabRef;
  }
  */
  bool qryoverlap = overlap(min(xes->qrystartpos, xes->qryendpos), max(xes->qrystartpos, xes->qryendpos),
			    min(xei->qrystartpos, xei->qryendpos), max(xei->qrystartpos, xei->qryendpos));  
  double invSizeQry = fabs(xei->qrystartpos - xei->qryendpos);
  double invNlabQry = double(abs(xei->qrystartidx - xei->qryendidx));
  if( !qryoverlap ) { //same for query
    invSizeQry += gapSizeQry;
    invNlabQry += gapNlabQry;
  }
  //max partial interval: this is what Luna calls nFarLab
  double maxPint = max(  fabs(Xmap->site[0][sv->querystartidx] - Xmap->site[0][sv->querystartidx-1]),
		         fabs(Xmap->site[0][sv->querystartidx] - Xmap->site[0][sv->querystartidx+1]));
  maxPint = max(maxPint, fabs(Xmap->site[0][sv->querystopidx ] - Xmap->site[0][sv->querystopidx -1]));
  maxPint = max(maxPint, fabs(Xmap->site[0][sv->querystopidx ] - Xmap->site[0][sv->querystopidx +1])) * 1e3; //bp
  
  //removed: invType vrlpprcntg12 vrlpprcntg
  //printf(" invpred: nMG12=%.0f vrlpsz12=%.1f nMG=%.1f minRawConf=%.1f vrlpsz=%.1f BPvrlpSize=%.1f gapSize=%.1f gapSizeRef=%.1f gapNlab=%.1f gapNlabRef=%.1f invSizeRef=%.1f invNlabRef=%.1f invSize=%.1f invNlab=%.1f nFarLab=%.1f\n", nmg11, overlapsize12, nmg, minRawConf, overlapsize12, overlapsize, gapSizeQry, gapSizeRef, gapNlabQry, gapNlabRef, invSizeRef, invNlabRef, invSizeQry, invNlabQry, maxPint); //debug
  printf(" invpred: occ_max=%.4f (%.4f) chimn_max=%.4f frag.max=%.4f nMG=%.1f nMG12=%.0f vrlpsz12=%.1f BPvrlpSize=%.1f gapNlab=%.1f gapSize=%.1f invSize=%.1f invNlab=%.1f minRawConf=%.2f nFarLab=%.1f\n", occ_max/medocc, occ_max, chimn_max, fragile_max, nmg, nmg11, overlapsize12, overlapsize, gapNlabQry, gapSizeQry, invSizeQry, invNlabQry, minRawConf, maxPint); //debug
  
  // link_coord == true means xei->qrystart/refstop is the partial coord (so false is qryend/refstart) (don't need)
  //int pidx = ( sv->link_coord ? xei->qrystartidx : xei->qryendidx );

  //store in pred
  pred[0 ] = occ_max/medocc;
  pred[1 ] = chimq_max;
  pred[2 ] = fragile_max;
  pred[3 ] = nmg;
  pred[4 ] = nmg11;
  pred[5 ] = overlapsize12;
  pred[6 ] = overlapsize;
  pred[7 ] = gapNlabQry;
  pred[8 ] = gapSizeQry;
  pred[9 ] = invSizeQry;
  pred[10] = invNlabQry;
  pred[11] = minRawConf;
  pred[12] = maxPint;
} //end getInvPredictors


//allocate memory for copy of matrix and return it -- not used; use below instead
//double* copyMatrix(double* matrix, int nrow, int ncol) {
//  double* copy = new double[ nrow * ncol ];
//  for( int i=0; i<nrow; i++ ) {
//    for( int j=0; j<ncol; j++ ) {
//      copy[ i*ncol + j ] = matrix[ i*ncol + j ];
//    }
//  }
//  return(copy);
//} //end copyMatrix


//allocate memory for copy of matrix and return its _transpose_
double* copyTranspose(double* matrix, size_t nrow, int ncol) {
  double* copy = new double[ nrow * ncol ];
  for( size_t i=0; i<nrow; i++ ) {
    for( int j=0; j<ncol; j++ ) {
      copy[ j*nrow + i ] = matrix[ i*ncol + j ];
    }
  }
  return(copy);
} //end copyTranspose


void printMatrix(double* matrix, size_t nrow, int ncol) {
  for( size_t i=0; i < nrow; i++ ) {
    for( int j=0; j < ncol; j++ ) {
      printf("%e,  ", matrix[ i*ncol + j ] );
    }
    printf("\n");
  }
} //end printMatrix


//multiply a*b, where a has r1 rows and c1 columns, and b has c1 rows and c2 columns, 
// and store result in r (which, obviously, must be pre-allocated and initialized to 0)
void multiplyMatrix(double* a, double* b, size_t r1, size_t c1, size_t c2, double* r) {
  for(size_t i = 0; i < r1; ++i)
    for(size_t j = 0; j < c2; ++j)
      for(size_t k = 0; k < c1; ++k){
	if(DEBUG>=2) assert(isfinite(a[i*c1 + k]));
	if(DEBUG>=2) assert(isfinite(b[k*c2 + j]));
	r[i*c2 + j] += a[i*c1 + k] * b[k*c2 + j];
	if(DEBUG>=2) assert(isfinite(r[i*c2 + j]));
      }
}


void doTransInvConf(xmapEntry **usexmapentries, bool dotrans=true) {
  if(!XXmap[0]->sitecnt[0])/* BNX maps : no coverage information available : cannot compute confidence_scaled (leave at original value of -1.0) */
    return;

  int verbose = 0;
  const int npred = (dotrans ? 11 : 13); //hardcode number of predictors here: this is size of array returned by getTransPredictors

  double medocc    = getQueryMapStat(); //default is median occurrence
  //double meanchimn = getQueryMapStat(false, false); //mean(ChimNorm) -- NOT USED
  if( dotrans ) {
    printf("Median occurrence = %f\n", medocc);
    //printf("Mean ChimNorm = %f\n", meanchimn);
    fflush(stdout);
  }
  if(medocc <= 0.0){
    if(VERB>= 2 && !dotrans){
      printf("doTransInvConf: Median occurrence = %f : skipping confidence computation\n",medocc);
      fflush(stdout);
    }
    return;
  }

  size_t ntrans = 0; //first count number of translocations (or inversions) : Also handle small inversions here (ntrans * ntrans can overflow 31 bit int)
  for(int indp=0; indp < numsv; indp++) {
    structuralVariation* sv = &sv_array[indp];
    if((dotrans && sv->type_enum != interchr_translocation && sv->type_enum != intrachr_translocation) || sv->duplicate)
      continue;

    if( !dotrans && sv->type_enum == inversion_small_type ) {    // small inversions are assigned confidence based on number of aligned labels in smaller MG
      if(DEBUG) assert(sv->algo_enum == sv_indel);
      xmapEntry *entry2 = xmapentries[sv->indices[1]];/* indices[1] is always the smaller MG for small inversions */
      Calign *p = entry2->align;
      int numpairs = p->numpairs;/* number of aligned labels in small MG */
      Cmap *Ymap = YYmap[p->mapid1];
      FLOAT *Y = Ymap->site[0];
      //      int N = Ymap->numsite[0];
      if(DEBUG) assert(YYmap[p->mapid1]->id == sv->refcontigid1);

      // should exclude MG overlapping breakpoint intervals (due to MG trimming). Also exclude aligned intervals under smap_Inv_ConfMinIntervalSize (kb) (on reference)

      // Only count intervals that do NOT overlap refstart..refstop breakpoint intervals on either  sv
      if(DEBUG) assert(numpairs >= 2);

      double refstart1 = sv->refstart * 1e-3, refstart2 = sv->refstart2 * 1e-3;
      double refstop1 = sv->refstop * 1e-3, refstop2 = sv->refstop2 * 1e-3;
      if(VERB>=2){
	printf("small inversion: qryid=%lld,refid=%lld:ref= %0.3f .. %0.3f AND %0.3f .. %0.3f, np=%d:\n",
	       XXmap[p->mapid2]->id,YYmap[p->mapid1]->id,refstart1,refstop1,refstart2,refstop2,numpairs);
	fflush(stdout);
      }

      int cnt = 0;
      int lastI = p->sites1[0],I;
      int lastK = p->sitesK1[0],K;
      double lastYik = Yc(Y,lastI,lastK),Yik;

      for(int i = 1; i < numpairs; i++, lastI = I, lastK = K, lastYik = Yik){
	I = p->sites1[i];
	K = p->sitesK1[i];
	Yik = Yc(Y,I,K);
	if(DEBUG) assert(I-K >= lastI);
	if(DEBUG) assert(Yik >= lastYik);
	if(Yik - lastYik < smap_Inv_ConfMinIntervalSize)
	  continue;
	if(overlap(refstart1,refstop1, lastYik, Yik))
	  continue;
	if(overlap(refstart2,refstop2, lastYik, Yik))
	  continue;

	cnt++;

	if(VERB>=2){
	  printf("\t Interval %d: Y= %0.3f .. %0.3f, I=%d,K=%d,lastI=%d,lastK=%d\n",cnt,lastYik,Yik,I,K,lastI,lastK);
	  fflush(stdout);
	}
      }

      cnt = min(cnt, INV_PAIRED_LEN - 1);
      sv->confidence_scaled = inv_paired_vec[cnt];
      if(cnt < smap_Inv_ConfMinNumIntervals)// NEW151
	sv->duplicate = true;
      continue;
    }

    if( !dotrans && sv->type_enum != inversion_type)
      continue;
    ntrans++;
  }
  if( ntrans <= 0 ){ //nothing to do
    if(VERB>=2 && !dotrans){
      printf("doTransInvConf: No inversions found : skipping confidence computation\n");
      fflush(stdout);
    }

    return;
  }

  double* pred = new double[ npred ]; //temp
  //matrix of predictors: rows are translocations (ntrans), columns are predictors (npred)
  double* preds = new double[ ntrans * npred ];
  double* temp1 = new double[ ntrans * npred ]; //temp for preds * cov_inv
  double* temp2 = new double[ ntrans * ntrans ]; //temp for predsT * cov_inv * pred
  for( int j=0; j < npred; j++ ) { //initialize
    pred[j] = 0.;
    for( size_t i= 0; i < ntrans; i++ ) {
      preds[ i*npred + j ] = 0.0;
      temp1[ i*npred + j ] = 0.0;
    }
  }  
  for( size_t j= 0; j < ntrans; j++ )
    for( size_t i= 0; i < ntrans; i++ ) 
      temp2[ i*ntrans + j ] = 0.0;

  //now get predictors
  size_t itr = 0;
  for(int indp= 0; indp < numsv; indp++) {
    structuralVariation* sv = &sv_array[indp];
    //this condition must match above loop (note: nicer to put this in a method of structuralVariation class)
    if((dotrans && sv->type_enum != interchr_translocation && sv->type_enum != intrachr_translocation) || sv->duplicate)
      continue;
    //small inversions are handled later
    if( !dotrans && sv->type_enum != inversion_type ) 
      continue;
    dotrans ?
      getTransPredictors(sv, usexmapentries, medocc, pred, npred) : 
      getInvPredictors  (sv, usexmapentries, medocc, pred, npred); //store result in pred
    bool allzero = true;
    for( int j= 0; j < npred; j++ ) {
      if( pred[j] != 0.0 )
	allzero = false;
      preds[ itr*npred + j ] = pred[j]; //row is itr (index of translocation), column is j
      if(DEBUG>=2 && !(isfinite(preds[itr * npred + j]))){
	printf("itr= %lu, npred=%d,ntrans=%lu,j=%d:itr*npred+j=%lu,pred[j]= %0.3e,preds[itr*npred+j]= %0.3e,dotrans=%d\n",
	       itr,npred,ntrans,j,itr*npred+j,pred[j],preds[itr*npred+j],dotrans);
	fflush(stdout);
	assert(isfinite(preds[itr * npred + j]));
      }
    }
    itr++;
    if( allzero ) {
      printf("Terminating doTransInvConf\n");
      fflush(stdout);
      return;
    }
  }
  delete [] pred; //now done with pred
  if( verbose > 1 ) { //DEBUG: print all predictors again
    printf("All predictors:\n");
    printMatrix(preds, ntrans, npred);
  }
  //now subtract trans_mean_vec from each row to get 'x-m'
  for( size_t i= 0; i < ntrans; i++ ) {
    for( int j= 0; j < npred; j++ ) {
      double orig = preds[i*npred + j];
      preds[ i*npred + j ] -= (dotrans ? trans_mean_vec[j] : inver_mean_vec[j]);
      if(DEBUG /* HERE >=2 */ && !(isfinite(preds[i * npred + j]))){
	printf("i=%lu,npred=%d,j=%d:i*npred+j=%lu,preds[i*npred+j]= %0.3e -> %0.3e,dotrans=%d,trans_mean_vec[j]= %0.3e,inver_mean_vec[j]= %0.3e\n",
	       i,npred,j,i*npred+j, orig, preds[i*npred+j], dotrans ? 1 : 0, trans_mean_vec[j], inver_mean_vec[j]);
	fflush(stdout);
	assert(isfinite(preds[i * npred + j]));
      }
    }
  }
  if( verbose > 1 ) { //DEBUG: print all predictors again
    printf("Subtracted predictors:\n");
    printMatrix(preds, ntrans, npred);
  }
  double* predsT = copyTranspose(preds, ntrans, npred); //allocates this array (free below)
  //printf("predsT:\n"); //DEBUG
  //printMatrix(predsT, npred, ntrans); //DEBUG
  dotrans ? 
    multiplyMatrix(preds, cov_inv    , ntrans, npred, npred, temp1) : //temp1 = preds * cov_inv
    multiplyMatrix(preds, cov_inv_inv, ntrans, npred, npred, temp1);
  //printf("temp1:\n"); //DEBUG
  //printMatrix(temp1, ntrans, npred); //DEBUG
  multiplyMatrix(temp1, predsT, ntrans, npred, ntrans, temp2); //this is D^2
  //printf("temp2:\n"); //DEBUG
  //printMatrix(temp2, ntrans, ntrans); //DEBUG
  delete [] preds;
  delete [] predsT;
  delete [] temp1;
  //last step is to get diagonals of D^2 (temp2), use those in ecdfx, and then assign to confidence data member
  //note left/right cofidence are set in getTransPredictors. There is no center confidence.
  itr = 0;
  for(int indp=0; indp < numsv; indp++) {
    structuralVariation* sv = &sv_array[indp];
    //this condition must match above loop
    if((dotrans && sv->type_enum != interchr_translocation && sv->type_enum != intrachr_translocation) || sv->duplicate)
      continue;
    if( !dotrans && sv->type_enum != inversion_type )
      continue;
    double d = temp2[ itr*ntrans + itr ];
    if(!RELEASE) assert(isfinite(d) && d > 0.0); //doesn't make sense, I think...
    int i=0;
    for( ; i < 100; i++ ) { //sizeof(ecdfx)/sizeof(double) == 100
      if( d < (dotrans ? ecdfx[i] : ecdfx_inv[i]) ) //note < is higher precedence than ?
	break;
    }
    //writeSVToSmap filters based on confidence: use confidence_scaled
    sv->confidence_scaled = double(100-i)/100.;
    printf("%sConf: GM %5lld: itr %lu: d %.2f i %2i conf %.2f\n", dotrans ? "trans" : "inv", usexmapentries[sv->indices[0]]->qrycontigid, itr, d, i, sv->confidence_scaled); //debug
    itr++;
  }
  delete [] temp2;
} //end doTransInvConf


//pair inversion breakpoints: arguments are inversion SV indices (1 < 2)
//if paired, store pairing data in sv1 and sv2 objects: (where?)
void pairInversionBP(xmapEntry **usexmapentries, int svind1, int svind2, bool verbose) {

  structuralVariation &sv1 = sv_array[svind1];
  structuralVariation &sv2 = sv_array[svind2];
  //note that refcontigid1 == refcontigid2, but check both anyway
  if( sv1.refcontigid1 != sv2.refcontigid1 || sv1.refcontigid2 != sv2.refcontigid2 )
    return;

  int i1 = (sv1.link_index == sv1.indices[0] ? sv1.indices[1] : sv1.indices[0]); //the non inverted on sv1
  int i2 = (sv2.link_index == sv2.indices[0] ? sv2.indices[1] : sv2.indices[0]); //the non inverted on sv2
  xmapEntry* xe1 = (sv1.algo_enum == sv_ps) ? usexmapentries[i1] : xmapentries[i1]; //the non inverted mg on sv1
  xmapEntry* xe2 = (sv2.algo_enum == sv_ps) ? usexmapentries[i2] : xmapentries[i2]; //the non inverted mg on sv2
  
  if (verbose)
    printf("pairInversionBP: chr %2lld: xids %3i-%3i %3i-%3i: ", sv1.refcontigid1, xe1->xmapid, sv1.link_index, xe2->xmapid, sv2.link_index);

  //case 1: same inverted mg
  if( sv1.link_index == sv2.link_index ) {
    //the only possibility for NOT pairing this case is if both the other MGs are on the same side of the inverted mg
    // these should be duplicates, but maybe not always? Check.
    xmapEntry* xei = (sv1.algo_enum == sv_ps) ? usexmapentries[sv1.link_index] : xmapentries[sv1.link_index]; //the inverted mg (both svs)
    assert( xe1->qrycontigid == xe2->qrycontigid && xe1->qrycontigid == xei->qrycontigid ); //dummy check

    if(verbose)
      printf("samemg ");
    if( (xe1->refstartpos < xei->refstartpos && xe2->refstartpos < xei->refstartpos) ||
	(xe1->refendpos > xei->refendpos && xe2->refendpos > xei->refendpos) ) {
      if(verbose)
	printf("conflict  :");
    } else {
      //store pair data (need to change sv_array object, sv1/2 are copies)
      sv_array[svind1].pair_index = svind2;
      sv_array[svind2].pair_index = svind1;
      if(verbose)
	printf("pair      :");
    }
    if(verbose)
      printf(" xe1: %10.1f - %10.1f xe2: %10.1f - %10.1f xei: %10.1f - %10.1f\n",
	     xe1->refstartpos, xe1->refendpos, xe2->refstartpos, xe2->refendpos, xei->refstartpos, xei->refendpos);
    return; //end of case1: do not check other cases
  }

  //common case which is never paired is when the two non-inverted MGs are 'the same.' May be identical coordinates, but better to allow some non-overlap. This is a sub-case of the 'xe1/2i outside xe2/1' cases below, but it may be easier to separate.
  //just do fixed 20 kb window (+-10) for start and stop
  double buf = 10e3; //bp, same as refstart/endpos
  //overlapPoint(double min1, double max1, double p) {
  if( overlapPoint(xe1->refstartpos-buf, xe1->refstartpos+buf, xe2->refstartpos) &&
      overlapPoint(xe1->refendpos  -buf, xe1->refendpos  +buf, xe2->refendpos) ) {
    if(verbose)
      printf("fail same mg: xe1: %10.1f - %10.1f xe2: %10.1f - %10.1f\n", xe1->refstartpos, xe1->refendpos, xe2->refstartpos, xe2->refendpos);
    return;
  }

  //if( !(xe1->refstartpos < xe2->refstartpos) ) { //note: this is NOT guaranteed by sorting because xe1/2 are NOT indices[0]: they are the non-inverted MGs which can be either indices
  //printf("fail assertion: xe1: %10.1f - %10.1f xe2: %10.1f - %10.1f\n", xe1->refstartpos, xe1->refendpos, xe2->refstartpos, xe2->refendpos); fflush(stdout);

  //distance requirement: maximum separation on reference must be < some value. For now, use smap_sv_sizekb (parameters.h/cpp), which is 5000 (kb): may want to lower this, although that means won't pair inversions larger than this
  //if( fabs(xe2->refstartpos - xe1->refendpos)/1e3 > smap_sv_sizekb ) { //pos in bp: convert to kb
  //LHS is > 0 for no overlap, < 0 for overlap (hence the -)
  if( -(min(xe1->refendpos,xe2->refendpos) - max(xe1->refstartpos,xe2->refstartpos))/1e3 > smap_sv_sizekb ) {
    if(verbose)
      printf("fail max sep: xe1: %10.1f - %10.1f xe2: %10.1f - %10.1f\n", xe1->refstartpos, xe1->refendpos, xe2->refstartpos, xe2->refendpos); //turn off this print since there are lots of them (?)
    return;
  }

  xmapEntry* xe1i = (sv1.algo_enum == sv_ps) ? usexmapentries[ sv1.link_index ] : xmapentries[ sv1.link_index ]; //inverted mg on sv1
  xmapEntry* xe2i = (sv2.algo_enum == sv_ps) ? usexmapentries[ sv2.link_index ] : xmapentries[ sv2.link_index ]; //inverted mg on sv2
  //bool pass = true;
  //xe1i must be between xe2 and xe2i: allow buffer of sv_indel_maxoverlap (50 kb)
  if( (xe1i->refstartpos < min(xe2->refstartpos, xe2i->refstartpos)/1e3 - sv_indel_maxoverlap) ||
      (xe1i->refendpos   > max(xe2->refendpos  , xe2i->refendpos  )/1e3 - sv_indel_maxoverlap) ) {
    //pass = false;
    if(verbose)
      printf("xe1i outside xe2s: ");
  }
  //and flip 1-2: xe2i must be between xe1 and xe1i
  else if( (xe2i->refstartpos < min(xe1->refstartpos, xe1i->refstartpos)/1e3 - sv_indel_maxoverlap) ||
	   (xe2i->refendpos   > max(xe1->refendpos  , xe1i->refendpos  )/1e3 - sv_indel_maxoverlap) ) {
    //pass = false;
    if(verbose)
      printf("xe2i outside xe1s: ");
  }
  else {
    if(verbose)
      printf("compatible xes   : ");
    sv_array[svind1].pair_index = svind2;
    sv_array[svind2].pair_index = svind1;
  }
  if(verbose)
    printf("xe1: %10.1f - %10.1f xe1i: %10.1f - %10.1f; xe2: %10.1f - %10.1f xe2i: %10.1f - %10.1f\n",
	   xe1->refstartpos, xe1->refendpos, xe1i->refstartpos, xe1i->refendpos,
	   xe2->refstartpos, xe2->refendpos, xe2i->refstartpos, xe2i->refendpos);
} //end pairInversionBP

//wrapper to pairInversionBP, also fix the pair_index for the inversion_small_types:
// it is broken by the qsort of sv_arry in output_smap which is before this is called
void doPairInversionBP(xmapEntry **usexmapentries, structuralVariation* sv_array, int numsvIndel, int numsv, bool verbose) {
  for(int indp= numsvIndel; indp < numsv; indp++) {
    if( !sv_array[indp].doInversionPartial() )
      continue;
    for(int indr= indp+1; indr < numsv; indr++) { //compare to subsequent
      if( !sv_array[indr].doInversionPartial() )
	continue;
      pairInversionBP(usexmapentries, indp, indr, verbose);
    }
  }
  if(verbose)
    fflush(stdout);
}
