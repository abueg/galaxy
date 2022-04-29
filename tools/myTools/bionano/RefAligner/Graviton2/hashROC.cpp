#include <stdlib.h>
#include <stdio.h>
#ifndef WIN32
#include <unistd.h>
#endif
#include <omp.h>

#include "parameters.h"
#include "Assembler_parameters.h"
#include "hash.h"
#include "Calign.h"

#undef TRACE

#define ROCSIZ 32 /* maximum hash score up to which ROC curve is computed */

#define FP_FIX 0 /* Fix scaling of FP scale by 2x : currently not turned on to keep previous ROC curves comparable */

#define RTRACE 0 /* trace specified IDs in input file :
		   1 : Only trace specific ID pair
		   2 : trace either ID in combination with any other map id 
		 */
static const long long RTraceID2 = 2727000461LL, RTraceID1 = 1060000525LL;/* Note: In the output  ID1 ID2 are reversed : make sure RTraceID1 < RTraceID2 */

extern double mtime();/* microsecond timer function : see refine.cpp */
extern double wtime();

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/hashROC.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

extern Cid2mapid *id2mapid;

static int CalignIdInc(Calign **p1, Calign **p2)
{
  long long p1_id1 = Gmap[p1[0]->mapid1]->id;
  long long p1_id2 = Gmap[p1[0]->mapid2]->id;
  long long p2_id1 = Gmap[p2[0]->mapid1]->id;
  long long p2_id2 = Gmap[p2[0]->mapid2]->id;

  /* Cannot use p1[0]->id - p2[0]->id since that may be out of range for integer */
  return (p1_id1 > p2_id1) ? 1 : (p1_id1 < p2_id1) ? -1 :
    (p1_id2 > p2_id2) ? 1 : (p1_id2 < p2_id2) ? -1 : (p1[0]->orientation - p2[0]->orientation);
}

void hashROC(char *hash_filename)
{
  char filename[PATH_MAX];/* for output filename(s) */

  if(VERB){
    printf("Completed map and alignment input: CPU time=%0.6f, wall time=%0.6f secs\n", mtime(), wtime());
    fflush(stdout);
  }

  if(MaxCoverageDump){
    /* sort alignment[0..numaligns-1] in order id1,id2,orientation */
    qsort(alignment,numaligns,sizeof(Calign*),(intcmp*)CalignIdInc);
    if(DEBUG) assert(!Cfirst);
    output_align(hash_filename,0);
    exit(0);
  }

  /* Note : The alignments have theirs id's translated into mapids[] based on the complete Gmap[0..nummaps-1] array, so to translate in the reverse direction
     there is no need to split the input array Gmap[] into XXmap[] and YYmap[] if Cfirst > 0 */

  /* Check alignment[0..numaligns-1] and make sure id1 < id2 (by swapping mapid1,mapid2 if needed : this invalidates alignment information, but this is never used except for numpairs) */
  for(register size_t i = 0; i < numaligns;i++){
    register Calign *palign = alignment[i];
    long long id1 = Gmap[palign->mapid1]->id;
    long long id2 = Gmap[palign->mapid2]->id;
    if(id1 > id2){/* swap mapid1,mapid2 */
      int mapid = palign->mapid1;
      palign->mapid1 = palign->mapid2;
      palign->mapid2 = mapid;
    }
    if(RTRACE && ((id1 == RTraceID2 && id2 == RTraceID1) || (id1 == RTraceID1 && id2 == RTraceID2))){
      printf("alignment[%lu]:id1=%lld,id2=%lld,or=%d:score=%0.4f,logPV=%0.2f,numpairs=%d\n",
	     i, Gmap[palign->mapid1]->id,Gmap[palign->mapid2]->id,palign->orientation,palign->score,palign->logPV, palign->numpairs);
      fflush(stdout);
    }
  }

  /* sort alignment[0..numaligns-1] in order id1,id2,orientation to match order in matches[0..matchcnt-1] */
  qsort(alignment,numaligns,sizeof(Calign*),(intcmp*)CalignIdInc);

  if(VERB){
    printf("Completed sorting of alignments: CPU time=%0.6f, wall time=%0.6f secs\n", mtime(), wtime());
    fflush(stdout);
  }

  if(RTRACE){/* locate alignment location */
    for(register size_t i = 0; i < numaligns; i++){
      register Calign *palign = alignment[i];
      long long id1 = Gmap[palign->mapid1]->id;
      long long id2 = Gmap[palign->mapid2]->id;
      if(RTRACE && ((id1 == RTraceID2 && id2 == RTraceID1) || (id1 == RTraceID1 && id2 == RTraceID2))){
	printf("alignment[%lu]:id1=%lld,id2=%lld,or=%d:score=%0.4f,logPV=%0.2f,numpairs=%d\n",
	       i, Gmap[palign->mapid1]->id,Gmap[palign->mapid2]->id,palign->orientation,palign->score,palign->logPV, palign->numpairs);
	fflush(stdout);
      }
    }
  }

  if(VERB){
    printf("Streaming in %s ...\n",hash_filename);
    fflush(stdout);
  }

  /* open buffered output stream (if needed for -hashtxt) */
  FILE *fptxt = NULL;
  if(HashTxt){
    sprintf(filename,"%s.hash",output_prefix);
    if(VERB){
      printf("Writing %s ...\n", filename);
      fflush(stdout);
    }
    
    checkFile(filename);
    if((fptxt = fopen(filename,"w"))==NULL){
      int eno = errno;
      char *err = strerror(errno);
      printf("failed to open file %s for writing:errno=%d:%s\n",filename,eno,err);
      exit(1);
    }
    
    /* write out commandline */
    printversion(fptxt);

    fprintf(fptxt,"# ID1 ID2 Orientation2 Score Offset(Left end of ID2 relative to Left end of ID1 in kb)\n");
  }

  /* open buffered input stream */
  FILE *fp;
  if((fp = fopen(hash_filename,"r"))==NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("failed to open file %s for reading:errno=%d:%s\n", hash_filename,eno,err);
    exit(1);
  }
  size_t L = strlen(hash_filename);
  int txt = 1;
  if(L >= strlen(".hashbin") && !strcmp(&hash_filename[L-strlen(".hashbin")],".hashbin"))/* binary file input */
    txt = 0;
  int linecnt = 0;/* only used of txt==1 */

  /* allocate input buffer */
  size_t matchcnt = 0;
  size_t matchmax = 64*1024;/*  1.5 Mbytes */
  CHashMatch *matches = new CHashMatch[matchmax];

  /* set file buffer to 4x the input buffer (should be enough even if file is in txt mode) */
  size_t fbufsiz = 4*matchmax*sizeof(CHashMatch);
  char *fbuf = new char[fbufsiz];
  setbuffer(fp, fbuf, fbufsiz);

  /* classify matches[0..matchcnt-1] as True Positive based on alignment[] being present */
  long long totalTPcnt = 0;
  int maxscore = 0;
  
  /* Also locate FN alignment with highest logPV */
  double FNmaxPV = 0.0;
  int FNindex = -1;

  /* Also Locate FP hashtable entry with highest score */
  int FPmaxscore = 0;
  CHashMatch FPmatch;

  /* Also update ROC curve (assume score <= ROCSIZ) : update TotalCnt[] and TPcnt[] */
  long long *TotalCnt = new long long[ROCSIZ+1];
  long long *TPcnt = new long long[ROCSIZ+1];
  for(int i = 0; i <= ROCSIZ;i++)
    TotalCnt[i] = TPcnt[i] = 0;

  size_t totalmatchcnt = 0;
  size_t j = 0;/* index into next alignment[] entry */

  CHashMatch lastmatch;/* for debugging */
  lastmatch.id1 = -1;
  lastmatch.id2 = -1;
  lastmatch.orientation = 0;
  lastmatch.hashscore = 0;

  while((matchcnt = hash_read(fp, txt, matches, matchmax, linecnt)) > 0){
    if(RTRACE){/* check for presence of id==RTraceID1 or id==RTraceID2 */
      for(size_t i = 0; i < matchcnt; i++){
	CHashMatch *p = &matches[i];
	if(RTRACE>=2){
	  if(p->id1 == RTraceID1 || p->id1 == RTraceID2 || p->id2 == RTraceID1 || p->id2 == RTraceID2)
	    printf("id1=%lld,id2=%lld,or=%d:score=%d,offset=%d(Kb)\n",(long long)p->id1,(long long)p->id2,p->orientation,p->hashscore,p->offset);
	} else if(p->id1 == RTraceID1 && p->id2 == RTraceID2)
	  printf("id1=%lld,id2=%lld,or=%d:score=%d,offset=%d(Kb)\n",(long long)p->id1,(long long)p->id2,p->orientation,p->hashscore,p->offset);
      }
      fflush(stdout);
    }
    if(DEBUG>=3){//  check for entires with id1 > id2 and check if the corresponding reversed pair is also present and has the same score 
      printf("Checking for errors in %lu matches\n", matchcnt);
      fflush(stdout);

      long long found = 0, err = 0;
      for(size_t i = 0; i < matchcnt; i++){
	if(i > 0 && !(i%1000000)){
	  printf("Checked %lu Matches, including %lld id reversals, %lld errors\n",i,found,err);
	  fflush(stdout);
	}
	CHashMatch *p = &matches[i];
	if(p->id1 > p->id2){
	  found++;
	  for(size_t k = 0; k < matchcnt; k++){
	    CHashMatch *q = &matches[i];
	    if(q->id2 == p->id1 && q->id1 == p->id2){
	      if(q->hashscore != p->hashscore){
		err++;
		printf("Inconsistent Matches:(id1=%lld,id2=%lld,score=%d,offset=%d),(id1=%lld,id2=%lld,score=%d,offset=%d)\n",
 		       (long long)p->id1, (long long)p->id2, p->hashscore, p->offset, (long long)q->id1,(long long)q->id2, q->hashscore, q->offset);
		fflush(stdout);
	      }
	      break;
	    }
	  }
	}
      }
      if(found){
	printf("Found %lld Matches with id1 > id2 and %lld score errors\n",found,err);
	exit(1);
      }
    }
    if(DEBUG){/* check ordering of matches[0] with lastmatch */
      assert(lastmatch.id1 <= matches[0].id1);
      if(lastmatch.id1 == matches[0].id1){
	assert(lastmatch.id2 <= matches[0].id2);
	if(lastmatch.id2 == matches[0].id2)
	  assert(lastmatch.orientation < matches[0].orientation);
      }
    }
    for(size_t i = 0; i < matchcnt; i++){
      if(DEBUG && i > 0){/* check for out of order entries */
	assert(matches[i-1].id1 <= matches[i].id1);
	if(matches[i-1].id1 == matches[i].id1){
	  if(!(matches[i-1].id2 <= matches[i].id2)){
	    printf("The following matches from %s are not in ascending order of id1,id2,orientation:\n",hash_filename);
	    for(size_t k = i-1; k <= i; k++)
	      printf("matches[%lu]:id1=%lld,id2=%lld,orientation=%d,score=%d,offset=%d\n",
		     k+totalmatchcnt,(long long)matches[k].id1,(long long)matches[k].id2,matches[k].orientation,matches[k].hashscore,matches[k].offset);
	    fflush(stdout);
	    assert(matches[i-1].id2 <= matches[i].id2);
	  }
	  if(matches[i-1].id2 == matches[i].id2){
	    if(!(matches[i-1].orientation < matches[i].orientation)){
	      printf("The following matches from %s are not in ascending order of id1,id2,orientation:\n",hash_filename);
	      for(size_t k = i-1; k <= i; k++)
		printf("matches[%lu]:id1=%lld,id2=%lld,orientation=%d,score=%d,offset=%d\n",
		       k+totalmatchcnt,(long long)matches[k].id1,(long long)matches[k].id2,matches[k].orientation,matches[k].hashscore,matches[k].offset);
	      fflush(stdout);
	      assert(matches[i-1].orientation < matches[i].orientation);
	    }
	  }
	}
      }

      CHashMatch *p = &matches[i];
      int score = p->hashscore;
      if(score > maxscore)
	maxscore = score;
      int bndscore = min(ROCSIZ,score);
      TotalCnt[bndscore]++;

      for(;j < numaligns;j++){
	Calign *palign = alignment[j];
	long long id1 = Gmap[palign->mapid1]->id;
	long long id2 = Gmap[palign->mapid2]->id;
	if(DEBUG) assert(id1 < id2);
	if(p->id1 > id1){/* False Negative for HashTable */
	  if(palign->logPV > FNmaxPV){
	    if(RTRACE && id1 == RTraceID1 && id2 == RTraceID2){
	      printf("False Negative for id1=%lld,id2=%lld,or=%d,score=%0.4f,logPV=%0.2f,numpairs=%d\n",
		     id1,id2,palign->orientation,palign->score,palign->logPV,palign->numpairs);
	      if(i > 0){
		printf("Previous 16 HashTable entries:\n");
		if(DEBUG && totalmatchcnt > 0 && i > 16)
		  printf("Match[%lu]:id1=%lld,id2=%lld,or=%d,score=%d,offset=%d\n", 
			 totalmatchcnt-1,(long long)lastmatch.id1,(long long)lastmatch.id2,lastmatch.orientation,lastmatch.hashscore,lastmatch.offset);
		for(size_t k = max((size_t)0,i-16); k < i; k++)
		  printf("Match[%lu]:id1=%lld,id2=%lld,or=%d,score=%d,offset=%d\n", 
			 k+totalmatchcnt,(long long)matches[k].id1,(long long)matches[k].id2,matches[k].orientation,matches[k].hashscore,matches[k].offset);
	      }
	      printf("Next 16 HashTable entries:\n");
	      for(size_t k = i; k < i+16 && k < matchcnt; k++)
		printf("Match[%lu]:id1=%lld,id2=%lld,or=%d,score=%d,offset=%d\n", 
		       k+totalmatchcnt,(long long)matches[k].id1,(long long)matches[k].id2,matches[k].orientation,matches[k].hashscore,matches[k].offset);
	      fflush(stdout);
	    }
	    FNmaxPV = palign->logPV;
	    FNindex = j;
	  }
	  continue;
	}
	if(p->id1 < id1){/* False Positive for Hashtable */
	  if(p->hashscore > FPmaxscore && !(j > 0 && Gmap[alignment[j-1]->mapid1]->id == p->id1 && Gmap[alignment[j-1]->mapid2]->id == p->id2)){
	    if(RTRACE && p->id1 == RTraceID1 && p->id2 == RTraceID2){
	      printf("False Positive for id1=%lld,id2=%lld,or=%d,score=%d,offset=%d(alignment[j=%lu]:id1=%lld,id2=%lld,or=%d)\n",(long long)p->id1,(long long)p->id2,p->orientation,p->hashscore,p->offset,j,id1,id2,palign->orientation);
	      if(j > 0){
		printf("Previous %lu Alignments: \n", min(j,(size_t)4));
		for(size_t k = max((size_t)0,j-4); k < j; k++){
		  Calign *qalign = alignment[k];
		  printf("alignment[%lu]:id1=%lld,id2=%lld,or=%d,score=%0.4f,logPV=%0.2f,numpairs=%d\n",
			 k,Gmap[qalign->mapid1]->id,Gmap[qalign->mapid2]->id,qalign->orientation,qalign->score,qalign->logPV,qalign->numpairs);
		}
	      }
	      if(j+1 < numaligns){
		printf("Next %lu Alignments:\n", min(numaligns-j-1,(size_t)4));
		for(size_t k = j; k < min(j+4,numaligns); k++){
		  Calign *qalign = alignment[k];
		  printf("alignment[%lu]:id1=%lld,id2=%lld,or=%d,score=%0.4f,logPV=%0.2f,numpairs=%d\n",
			 k,Gmap[qalign->mapid1]->id,Gmap[qalign->mapid2]->id,qalign->orientation,qalign->score,qalign->logPV,qalign->numpairs);
		}
	      }
	    }
	    FPmaxscore = p->hashscore;
	    FPmatch = p[0];
	  }
	  break;
	}

	if(p->id2 > id2){/* False Negative for HashTable */
	  if(palign->logPV > FNmaxPV){
	    if(RTRACE && id1 == RTraceID1 && id2 == RTraceID2){
	      printf("False Negative for id1=%lld,id2=%lld,or=%d,score=%0.4f,logPV=%0.2f,numpairs=%d\n",
		     id1,id2,palign->orientation,palign->score,palign->logPV,palign->numpairs);
	      if(i > 0)
		printf("Previous HashTable entry: id1=%lld,id2=%lld,or=%d,score=%d,offset=%d\n", (long long)matches[i-1].id1,(long long)matches[i-1].id2,matches[i-1].orientation,matches[i-1].hashscore,matches[i-1].offset);
	      printf("Next HashTable entry: id1=%lld,id2=%lld,or=%d,score=%d,offset=%d\n", (long long)p->id1,(long long)p->id2,p->orientation,p->hashscore,p->offset);
	    }
	    FNmaxPV = palign->logPV;
	    FNindex = j;
	  }
	  continue;
	}
	if(p->id2 < id2){/* False Positive for HashTable */
	  if(p->hashscore > FPmaxscore && !(j > 0 && Gmap[alignment[j-1]->mapid1]->id == p->id1 && Gmap[alignment[j-1]->mapid2]->id == p->id2)){
	    if(RTRACE && p->id1 == RTraceID1 && p->id2 == RTraceID2){
	      printf("False Positive for id1=%lld,id2=%lld,or=%d,score=%d,offset=%d(alignment[j=%lu]:id1=%lld,id2=%lld,or=%d)\n",(long long)p->id1,(long long)p->id2,p->orientation,p->hashscore,p->offset,j,id1,id2,palign->orientation);
	      if(j > 0){
		printf("Previous %lu Alignments: \n", min(j,(size_t)4));
		for(size_t k = max((size_t)0,j-4); k < j; k++){
		  Calign *qalign = alignment[k];
		  printf("alignment[%lu]:id1=%lld,id2=%lld,or=%d,score=%0.4f,logPV=%0.2f,numpairs=%d\n",
			 k,Gmap[qalign->mapid1]->id,Gmap[qalign->mapid2]->id,qalign->orientation,qalign->score,qalign->logPV,qalign->numpairs);
		}
	      }
	      if(j+1 < numaligns){
		printf("Next %lu Alignments:\n", min(numaligns-j-1,(size_t)4));
		for(size_t k = j; k < min(j+4,numaligns); k++){
		  Calign *qalign = alignment[k];
		  printf("alignment[%lu]:id1=%lld,id2=%lld,or=%d,score=%0.4f,logPV=%0.2f,numpairs=%d\n",
			 k,Gmap[qalign->mapid1]->id,Gmap[qalign->mapid2]->id,qalign->orientation,qalign->score,qalign->logPV,qalign->numpairs);
		}
	      }
	      /*	      if(j > 0)
			      printf("Gmap[alignment[j-1]->mapid1]->id=%lld,Gmap[alignment[j-1]->mapid2]->id=%lld,alignment[j-1]->or=%d\n",Gmap[alignment[j-1]->mapid1]->id,Gmap[alignment[j-1]->mapid2]->id,alignment[j-1]->orientation);*/
	    }
	    FPmaxscore = p->hashscore;
	    FPmatch = *p;
	  }
	  break;
	}
	if(p->orientation > palign->orientation){/* False Negative for HashTable */
	  if(palign->logPV > FNmaxPV){
	    if(RTRACE && id1 == RTraceID1 && id2 == RTraceID2){
	      printf("False Negative for id1=%lld,id2=%lld,or=%d,score=%0.4f,logPV=%0.2f,numpairs=%d\n",
		     id1,id2,palign->orientation,palign->score,palign->logPV,palign->numpairs);
	      if(i > 0)
		printf("Previous HashTable entry: id1=%lld,id2=%lld,or=%d,score=%d,offset=%d\n", (long long)matches[i-1].id1,(long long)matches[i-1].id2,matches[i-1].orientation,matches[i-1].hashscore,matches[i-1].offset);
	      printf("Next HashTable entry: id1=%lld,id2=%lld,or=%d,score=%d,offset=%d\n", (long long)p->id1,(long long)p->id2,p->orientation,p->hashscore,p->offset);
	    }
	    FNmaxPV = palign->logPV;
	    FNindex = j;
	  }
	  continue;
	}
	if(p->orientation < palign->orientation){/* False Positive for HashTable */
	  if(p->hashscore > FPmaxscore && !(i+1 < matchcnt && p[1].id1==p->id1 && p[1].id2==p->id2 && p[1].orientation == palign->orientation)){
	    if(RTRACE && p->id1 == RTraceID1 && p->id2 == RTraceID2){
	      printf("False Positive for id1=%lld,id2=%lld,or=%d,score=%d,offset=%d(i=%lu,matchcnt=%lu,alignment[j=%lu]:id1=%lld,id2=%lld,or=%d)\n",
		     (long long)p->id1,(long long)p->id2,p->orientation,p->hashscore,p->offset,i,matchcnt,j,id1,id2,palign->orientation);
	      if(j > 0){
		printf("Previous %lu Alignments: \n", min(j,(size_t)4));
		for(size_t k = max((size_t)0,j-4); k < j; k++){
		  Calign *qalign = alignment[k];
		  printf("alignment[%lu]:id1=%lld,id2=%lld,or=%d,score=%0.4f,logPV=%0.2f,numpairs=%d\n",
			 k,Gmap[qalign->mapid1]->id,Gmap[qalign->mapid2]->id,qalign->orientation,qalign->score,qalign->logPV,qalign->numpairs);
		}
	      }
	      if(j+1 < numaligns){
		printf("Next %lu Alignments:\n", min(numaligns-j-1,(size_t)4));
		for(size_t k = j; k < min(j+4,numaligns); k++){
		  Calign *qalign = alignment[k];
		  printf("alignment[%lu]:id1=%lld,id2=%lld,or=%d,score=%0.4f,logPV=%0.2f,numpairs=%d\n",
			 k,Gmap[qalign->mapid1]->id,Gmap[qalign->mapid2]->id,qalign->orientation,qalign->score,qalign->logPV,qalign->numpairs);
		}
	      }
	      if(i > 0)
		printf("Previous HashTable entry: id1=%lld,id2=%lld,or=%d,score=%d,offset=%d\n", (long long)p[-1].id1,(long long)p[-1].id2,p[-1].orientation,p[-1].hashscore,p[-1].offset);
	      else if(totalmatchcnt > 0)
		printf("Previous HashTable entry: id1=%lld,id2=%lld,or=%d,score=%d,offset=%d\n", (long long)lastmatch.id1,(long long)lastmatch.id2,lastmatch.orientation,lastmatch.hashscore,lastmatch.offset);
	      if(i+1 < matchcnt)
		printf("Next HashTable entry: id1=%lld,id2=%lld,or=%d,score=%d,offset=%d\n", (long long)p[+1].id1,(long long)p[+1].id2,p[+1].orientation,p[+1].hashscore,p[+1].offset);
	    }
	    FPmaxscore = p->hashscore;
	    FPmatch = *p;
	  }
	  break;
	}

	/* found a match */
	TPcnt[bndscore]++;
	j++;
	totalTPcnt++;
	break;
      }
    }

    if(HashTxt){/* output hashtable in text format */
      for(size_t i = 0; i < matchcnt; i++){
	CHashMatch *phash = &matches[i];
	fprintf(fptxt, "%lld %lld %d %d %d\n", (long long)phash->id1, (long long)phash->id2, phash->orientation ? -1 : 1, phash->hashscore, phash->offset);
      }
    }

    totalmatchcnt += matchcnt;
    if(VERB && !(totalmatchcnt % (1000000000))){
      printf("Read %lu Billion matches: wall time=%0.6f\n", totalmatchcnt/(1000000000), wtime());
      fflush(stdout);
    }
    if(matchcnt < matchmax)/* last read reached EOF */
      break;
    if(DEBUG) lastmatch = matches[matchcnt-1];/* to check ordering with next input match */
  }
  fclose(fp);
  delete [] fbuf;
  delete [] matches;
  if(HashTxt)
    FILEclose(fptxt);

  if(VERB){
    printf("Read %lu Matches from %s: %lld TP entries, maxscore=%d: CPU time=%0.6f, wall time=%0.6f\n", totalmatchcnt, hash_filename, totalTPcnt, maxscore,mtime(),wtime());
    fflush(stdout);
  }
  
  int bmaxscore = min(ROCSIZ,maxscore);/* limit ROC curve range to 1..min(ROCSIZ,maxscore) */

  /* accumulate the counts to make it an ROC curve */
  for(int j = bmaxscore; --j >= 1;){
    TotalCnt[j] += TotalCnt[j+1];
    TPcnt[j] += TPcnt[j+1];
  }

  if(FPmaxscore > 0){
    if(DEBUG) assert(FPmaxscore == FPmatch.hashscore);
#if ID_HASH==0
    int index1 = findid(FPmatch.id1,id2mapid,nummaps);
    int index2 = findid(FPmatch.id2,id2mapid,nummaps);
#else
    int index1 = HashFindID(FPmatch.id1);
    int index2 = HashFindID(FPmatch.id2);
#endif
    
    if(index1 >= 0 && index2 >= 0){
      Cmap *p1 = Gmap[index1];
      Cmap *p2 = Gmap[index2];
      printf("Highest scoring FP in hashtable: id1=%lld(len=%0.3f,sites=%d) id2=%lld(len=%0.3f,sites=%d) orientation=%d hashscore=%d offset=%d(Kb): \n", 
	     (long long)FPmatch.id1, p1->site[0][p1->numsite[0]+1], p1->numsite[0], (long long)FPmatch.id2, p2->site[0][p2->numsite[0]+1], p2->numsite[0], FPmatch.orientation, FPmaxscore, FPmatch.offset);
    } else
      printf("Highest scoring FP in hashtable: id1=%lld id2=%lld orientation=%d score=%d offset=%d(Kb): \n", 
	     (long long)FPmatch.id1, (long long)FPmatch.id2, FPmatch.orientation, FPmaxscore, FPmatch.offset);
  }
  
  if(FNindex >= 0){
    Calign *palign = alignment[FNindex];
    Cmap *p1 = Gmap[palign->mapid1];
    Cmap *p2 = Gmap[palign->mapid2];

    printf("Highest scoring alignment not in hashtable: id1=%lld(len=%0.3f,sites=%d),id2=%lld(len=%0.3f,sites=%d),orientation=%d,score=%0.3f,logPV=%0.2f,numpairs=%d\n",
	   Gmap[palign->mapid1]->id,p1->site[0][p1->numsite[0]+1],p1->numsite[0],Gmap[palign->mapid2]->id, p2->site[0][p2->numsite[0]+1],p2->numsite[0], palign->orientation, palign->score, palign->logPV, palign->numpairs);
  }
  fflush(stdout);

  /* now output the ROC curve */
  sprintf(filename,"%s.ROC",output_prefix);
  
  if(VERB){
    printf("Writing %s : CPU time=%0.6f, wall time=%0.6f secs\n", filename, mtime(), wtime());
    fflush(stdout);
  }

  checkFile(filename);

  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(errno);
    printf("failed to open file %s for writing .cmap file of %d maps:errno=%d:%s\n",filename, nummaps,eno,err);
    exit(1);
  }

  /* write out commandline */
  char cwd[PATH_MAX];
  fprintf(fp,"# $ cd %s; %s", getcwd(cwd,PATH_MAX), gargv[0]);
  for(register int i = 1; i < gargc; i++)
    fprintf(fp," %s", gargv[i]);
  fprintf(fp,"\n");
  printversion(fp);

  for(int j = 1; j <= bmaxscore; j++)
    fprintf(fp,"%10lld ",TPcnt[j]);
  fprintf(fp,"\n");
  for(int j = 1; j <= bmaxscore; j++)
    fprintf(fp,"%10lld ",TotalCnt[j] - TPcnt[j]);
  fprintf(fp,"\n");
  fprintf(fp,"%10lu\n", numaligns);

  /* NOTE : totalFN should  not have /2 but we need to keep ROC curves comparable to past ROC curves which had the same error. The true FP rate is half the reported rate! */
  long long CfirstLL = Cfirst;
  long long totalFN = ((CfirstLL > 0) ? CfirstLL * (nummaps-CfirstLL) * 2LL : nummaps*(nummaps-1LL))/(FP_FIX ? 1LL : 2LL) - totalTPcnt;
  fprintf(fp,"%10lld\n", totalFN);

  /* compute normalization factors */
  double TPnorm = numaligns;
  double FPnorm = totalFN;
  double TPnormInv = 1.0/TPnorm;
  double FPnormInv = 1.0/FPnorm;

  fprintf(fp,"# score TPnorm FPnorm\n");
  for(int j = 1; j <= bmaxscore; j++)
    fprintf(fp,"%d %0.6e %0.6e\n", j, TPcnt[j] * TPnormInv, (TotalCnt[j] - TPcnt[j]) * FPnormInv);
  fprintf(fp,"\n");

  (void)FILEclose(fp);

  delete [] TPcnt;
  delete [] TotalCnt;
}
