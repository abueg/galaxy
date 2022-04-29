#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "parameters.h"

#include "CalignA.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/CalignA.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

/** same as maxalignalloc, but each new element is initialized to NULL */
void maxalignallocNULL(size_t num, CalignA **&alignmentA, size_t &numaligns, size_t &maxaligns)
{
  size_t newaligns = num;
  if(newaligns <= maxaligns)
    return;

  if(newaligns < 1024)
    newaligns = 1024;
  if(newaligns*4 < maxaligns*5)
    newaligns = (maxaligns*5)/4;

  if(VERB>=2){
    printf("maxalignallocNULL(%lu):maxaligns=%lu->%lu\n",num,maxaligns,newaligns);
    fflush(stdout);
  }

  if(maxaligns <= 0){
    CalignA **newalignment = new CalignA*[newaligns];
    for(register size_t i=0;i < newaligns; i++)
      newalignment[i] = NULL;
    alignmentA = newalignment;
  } else {
    CalignA **origalignment = alignmentA;
    CalignA **newalignment = new CalignA*[newaligns];
    for(register size_t i=0;i<maxaligns;i++)
      newalignment[i] = origalignment[i];
    for(register size_t i=maxaligns;i<newaligns;i++)
      newalignment[i] = NULL;
    alignmentA = newalignment;
    delete [] origalignment;
  }
  maxaligns = newaligns;
}


extern double ChiSqPvalue(int r, double ChiSq);// see Calign.cpp

/** depends on command line parameters RepeatMaxShift and RepeatPvalueRatio */
int CalignA::IsRepeatRegion()
{
  int repeat = 0;// If repeat had been > 0, it would already have been filtered out before being read in by Assembler_input()

  if(spots_filename[0] || (svcheck && numrefmaps > 0)){
    printf("ERROR: CalignA::IsRepeatRegion() : Only implemented for pairwise alignment\n");
    fflush(stdout);exit(1);
  }

  if(PairSplit || PairMerge || RMS_MIN < 0.0 /* HERE HERE*/||colors==2){ /* NOTE : If RMS_MIN >= 0, pairwise alignments are checked for duplicate maps */
    if(RepeatMaxShift <= 0)
      return 0;
    if(repeat >= 0)/* repeat already checked using -RepeatRec (superior method) */
      return repeat;
    if(DEBUG) assert(!RepeatRec);
  }

  if(colors > 1){
    fprintf(stderr,"IsRepeatRegion(): not implemented for -colors %d\n",colors);
    exit(1);
  }

  int K = numpairs-1;
  if(K <= 1)
    return 0;

  int result = 0;

  Cmap *Ymap = YYmap[mapid1];
  Cmap *Xmap = XXmap[mapid2];
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

  if(1){  /* pairwise Alignment OR Guided Alignment without resolution modeling */
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

  if(1){// !spot_filename[0]:  only difference is that vary = varF + var*(x+y) instead of vary = varF + var*y */
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
	printf("WARNING:IsRepeatRegion():mapid1=%d(id=%lld),mapid2=%d(id=%lld),or=%d:score=%0.6f,logPV=%0.2f,numpairs=%d:K=%d,cnt=%d:\n",
	       mapid1,Ymap->id,mapid2,Xmap->id,orientation,score,logPV,numpairs,K,cnt);
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
    double norm;

    if((norm = (cnt>0) ? sum/cnt : 1.0) < RMS_MIN && cnt >= RMS_SITES && abs(M-N) <= RMS_MIS && numpairs >= min(N,M) - RMS_MIS && !orientation && Pvalue1 < RMS_PVALUE){
      if(VERB){
	printf("WARNING:mapid1=%d(id=%lld from %s),mapid2=%d(id=%lld from %s),or=%d,N=%d,M=%d:numpairs=%d,score=%0.2f,logPV=%0.2f,RMS norm=%0.6f(DF=%d,Pvalue=%0.6e): Maps are almost identical\n",
	       mapid1,Ymap->id,vfixx_filename[Ymap->fileid],mapid2,Xmap->id,vfixx_filename[Xmap->fileid],orientation,N,M,numpairs,score,logPV,norm,cnt, 1.0-Pvalue1);
	fflush(stdout);
      }
    }

    if(RepeatMaxShift <= 0){/* repeat already checked using -RepeatRec (superior method) OR repeat checking disabled */
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

void copy(CalignA *dest, CalignA *pinit, int score, int deepcopy){
  *dest = *pinit;/* first copy all scalars */

  if(!deepcopy){
    dest->block_id = -1;
    return;
  }

  dest->sites1 = NULL;
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


int num_alignA_blocks = 0;/* number of CalignA_block(s) allocated */
CalignA_block *alignA_block[MAX_ALIGNBLOCKS];

/* allocate a CalignA_block to hold "num" alignments with a total of "sum" aligned sites. Append the pointer to align_blocks[0..num_align_blocks-1] and return it. */
CalignA_block *alignA_blockalloc(size_t num, size_t sum, size_t Msiz, int AlignScore, int sitesK, int &block_id)
{
  if(VERB>=2){
    printf("alignA_blockalloc:time=%0.6f(wall time=%0.6f)\n",mtime(),wtime());
    fflush(stdout);
  }

  if(num_alignA_blocks >= MAX_ALIGNBLOCKS){
    printf("Exceeded maximum of %d alignmentA blocks\n",MAX_ALIGNBLOCKS);
    exit(1);
  }
  CalignA_block *p = alignA_block[num_alignA_blocks] = new CalignA_block(num, sum, Msiz, AlignScore, sitesK,num_alignA_blocks);
  block_id = num_alignA_blocks++;
  return p;
}

/* compact a CalignA_block alignA_block[block_index] to only include alignmentA[block->start .. block->end-1] != NULL, all of which must point to memory in alignA_block[block_index].
   numpairs[] and numpairs_cum[] are pre-allocated arrays at least as large as (block->start - block->end)
   relocate the remaining "num" entries to alignmentA[start .. start + num - 1], where start <= block->start and alignmentA[start .. block->start - 1] are all NULL
   alignA_block[block_index] will be replaced with a reduced size memory block, and the original memory block deleted.

   returns num, the compacted number of alignments at alignmentA[start .. start + num - 1]
 */
size_t align_blockcompact(int block_index, size_t start, int *numpairs, size_t *numpairs_cum, int AlignScore, int numthreads, CalignA **alignmentA)
{
  CalignA_block *origblock = alignA_block[block_index];
  size_t origstart = origblock->start;
  size_t origend = origblock->end;

  /* remove alignments from alignmentA[start .. end-1] that are NULL*/
  size_t j = start, k = origstart;
  if(start == origstart){
    for(; k < origend; k++){
      if(!alignmentA[k]){
	k++;
	break;
      }
      j++;
    }
  }
  for(; k < origend; k++){
    if(!alignmentA[k])
      continue;
    if(DEBUG>=2) assert(j < k);
    if(DEBUG>=2) assert(alignmentA[j] == NULL);
    alignmentA[j] = alignmentA[k];
    alignmentA[k] = NULL;
    j++;
  }
  if(VERB>=2){
    printf("alignA_block[%d]: numaligns reduced from %lu to %lu:time=%0.6f(wall time=%0.6f)\n",block_index, origend - origstart, j-start, mtime(),wtime());
    fflush(stdout);
  }
  size_t end = j;

  size_t num = end - start;
  size_t numpairs_sum = 0;

  for(size_t i = 0; i < num; i++){
    CalignA *p = alignmentA[start + i];
    numpairs[i] = p->numpairs;

    numpairs_cum[i] = numpairs_sum;
    numpairs_sum += numpairs[i];
  }
  numpairs_cum[num] = numpairs_sum;

  CalignA_block *block = alignA_block[block_index] = new CalignA_block(num, numpairs_sum, 0ul, AlignScore, 0, block_index);
  block->start = start;
  block->end = end;

  long long numLL = num;

  size_t clen = (colors==1) ? 1 : colors + 1;

#pragma omp parallel for num_threads(numthreads) schedule(static,512)
  for(long long i = 0; i < numLL; i++){
    CalignA *p = &block->block[i*clen];

    size_t nump = numpairs[i];
    size_t offset = 2 * numpairs_cum[i];
    p->block_id = block_index;
    p->allocated_size = nump;
    p->sites1 = &block->intblock[offset];
    p->sites2 = &block->intblock[offset+nump];
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
  for(long long i = 0; i < numLL; i++){/* copy data from alignmentA[start + i] to block->block[i*clen] */
    CalignA *q = alignmentA[start+i];
    CalignA *p = &block->block[i*clen];

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
    
    alignmentA[start+i] = p;
  }

  delete origblock;

  return num;
}

/* free all CalignA_block(s) */
void alignA_blockfree()
{
  for(int i = 0; i < num_alignA_blocks; i++)
    delete alignA_block[i];
  num_alignA_blocks = 0;
}


void expand_arrays(CalignA *p, int num, int AlignScore, CalignA_block *block, int block_index){
  if(num <= p->allocated_size) return;
  if(p->allocated_size > 0 && p->block_id < 0){
    delete [] p->sites1;
    delete [] p->sites2;
    delete [] p->iscore;
    delete [] p->outscore;
  }
  p->block_id = block_index;
  p->allocated_size = num;
  p->sites1 = block->int_alloc(num);
  p->sites2 = block->int_alloc(num);
  if(AlignScore){
    p->iscore = block->FLOAT_alloc(num);
    p->outscore = block->FLOAT_alloc(num);
  } else
    p->iscore = p->outscore = NULL;
}

