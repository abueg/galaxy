#include "globals.h" //globals.h includes xmapEntry.h
#include "parameters.h"
#include "structuralVariation.h"

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/simpleRepeat.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

//if oldsize > 0, double the size of the argument array, return the new size
//otherwise just make single new object
int growRepeatArray(repeat*& ina, int oldsize) {
  repeat* tmp = ina; //make pointer to start of old array
  int newsize = (oldsize > 0 ? oldsize*2 : 1); //new size is twice old size as long as old size > 0, else 1
  ina = new repeat[newsize]; //allocate new size in old object ina
  for(int i=0; i<oldsize; i++) //copy eles of old into new
    ina[i] = tmp[i];
  if( oldsize ) //if no elements, can't delete
    delete [] tmp; //tmp is old array--container no longer necessary
  return newsize;
}


//templated version of above -- moved to structuralVariation.h
/*
template<class T> 
int growArray(T*& ina, int oldsize) {
  T* tmp = ina;
  int newsize = (oldsize > 0 ? oldsize*2 : 1); //new size is twice old size as long as old size > 0, else 1
  ina = new T[newsize];
  for(int i=0; i<oldsize; i++) //copy eles of old into new
    ina[i] = tmp[i];
  if( oldsize ) //if no elements, can't delete
    delete [] tmp; //tmp is old array--container no longer necessary
  return newsize;
}
*/

//version of above for array of pointers
template<class T> 
int growPointerArray(T**& ina, int oldsize) {
  T** tmp = ina;
  int newsize = (oldsize > 0 ? oldsize*2 : 1); //new size is twice old size as long as old size > 0, else 1
  ina = new T*[newsize];
  for(int i=0; i<oldsize; i++) //copy eles of old into new
    ina[i] = tmp[i];
  if( oldsize ) //if no elements, can't delete
    delete [] tmp; //tmp is old array--container no longer necessary
  return newsize;
} //end growArray


double total_map_length_kb; //global to this file, sum in findRepeatSimpleMaps, and output in output_repeatStats

//call findRepeatSimple on all -i arguments -- data stored in global repeat_array
void findRepeatSimpleMaps(int filtertype, float tolerance, int minele, bool verbose) {

  printf("Starting simple repeat detection on -i inputs\n"); //always print this
  total_map_length_kb = 0;

  Cmap **partmaps = 0;
  int numpartmaps = 0; //number of elements in partmaps
  int maxpartmaps = 0; //size of partmaps array
  if( filtertype == 1 || filtertype == 2 ) {
    maxpartmaps = max(1,int(nummaps/10.)); //just start with 10% and grow in this case
    partmaps = new Cmap *[maxpartmaps];
  }
  else if( filtertype == 3 ) { //here, all maps are output, so there's no need for a new array--point to maps
    partmaps = Gmap;
    numpartmaps = nummaps;
  }

  for( int i=0; i<nummaps; i++ ) {
    //note: for site and numsite, have color index--should always be 0 because of usecolor
    //also note that the units of site are kb, but since findRepeatSimple works in relative distance, it doesn't matter
    //more importantly, zeroth element is zero, so skip that, and numsite'th element is last (not one less that, as normal)

    if( verbose )
      printf("Starting map %lli\n", Gmap[i]->id); 

    //for( int j=1; j<=Gmap[i]->numsite[0]+1; j++) printf("%i: %.4f\n", j, Gmap[i]->site[0][j]); //debug
    //int nrepeats = findRepeatSimple(Gmap[i]->site[0]+1, Gmap[i]->numsite[0], tolerance, minele); //OLD, no FP/FN
    int nrepeats = findSimpleRepeatFPFN(Gmap[i]->site[0]+1, Gmap[i]->numsite[0], tolerance, minele );
    //if( nrepeats == -1 ) { printf("Error in map %lli\n", Gmap[i]->id); return; } //debug

    double maplen = Gmap[i]->site[0][Gmap[i]->numsite[0]+1]; //see Cmap.h
    total_map_length_kb += maplen;

    if( filtertype > 0 ) {
      if( (nrepeats == 0 && filtertype == 1) ||
	  (nrepeats >  0 && filtertype == 2) ) {
	if( numpartmaps == maxpartmaps )
	  maxpartmaps = growPointerArray(partmaps, numpartmaps);
	partmaps[numpartmaps] = Gmap[i];
	numpartmaps += 1;
      }
      else if( filtertype == 3 ) 
	maskRepeatLabels(Gmap[i], numrepeats-nrepeats, numrepeats, verbose);
    }

    //these fields are not populated in findRepeatSimple, here instead
    for( int j=numrepeats-nrepeats; j<numrepeats; j++ ) {
      repeat_array[j].id     = Gmap[i]->id;
      repeat_array[j].length = maplen;
    }
  } //end loop on maps
  
  for( int i=0; verbose && i<numrepeats; i++ )
    printf("repeat %2i: mapID %2lli: %.4f  %.4f  %2i  %.1f\n", i+1, repeat_array[i].id, repeat_array[i].startpos, repeat_array[i].endpos, repeat_array[i].nele, repeat_array[i].avgsize);

  if( verbose )
    printf("Done with simple repeat detection\n");

  output_repeatAll(filtertype, partmaps, numpartmaps); //see below, simpleRepeat.h
}
//end findRepeatSimpleMaps


//given a map and indices of repeat objects in repeat_arry, remove labels in those repeat regions
// remove from repeat_start to repeat_end - 1, ie, don't include repeat_end (it may not exist),
// and don't remove the first and last labels in the repeat--these are the border labels
// except if these are the first or last labels in the map, because then they can't align and will
// penalize the alignments
//Since the map is edited in place, there's no need to return anything--just pass the Cmap* by
// reference.
//If repeat_start == repeat_end, there are no repeats, so nothing to do.
void maskRepeatLabels(Cmap*& inmap, int repeat_start, int repeat_end, bool verbose) {
  if( repeat_start == repeat_end )
    return;

  if( verbose )
    printf("Removing %i repeat regions from molecule %lli\n", repeat_end-repeat_start, inmap->id);

  int icol = 0; //color index
  //for( int i=1; i<=inmap->numsite[icol]+1; i++) printf("%i: %.4f\n", i, inmap->site[icol][i]);
  //Loop on repeat_array, and for each one, edit the sites--shift array by simply assigning values
  // and changing length (extra values at end will be assigned 0 to be safe)
  //note that number of sites is not checked, ie, stopidx - startidx >= minele (argument to findRepeatSimple)
  for( int irep=repeat_start; irep<repeat_end; irep++ ) {
    //need to get index in the map of the start and end and shift accordingly
    int startidx = 0, stopidx = 0;
    for(int isit=1; isit <= inmap->numsite[icol]; isit++) { //array is 0 to N+1, but 0 and N+1 are not sites (1 through N are)
      if( repeat_array[irep].startpos == inmap->site[icol][isit] ) {
	assert( startidx == 0 ); //this means that two sites have the same position, which shouldn't be possible
	if( isit == 1 ) //if first site is start of repeat, remove it
	  startidx = isit;
	else //otherwise, remove starting at next site
	  startidx = isit+1; //this will fail only if isit is the last entry in the array
      }
      else if( repeat_array[irep].endpos == inmap->site[icol][isit] ) {
	assert( stopidx == 0 ); //same here
	if( isit == inmap->numsite[icol] ) //if last site is end of repeat, remove it
	  stopidx = isit;
	else //if not, keep it, remove previous
	  stopidx = isit-1;
	break; //once you found the end of this repeat, you're done
      } 
    } //end for sites in inmap
    //now the indices are in startidx and stopidx--remove them
    if( verbose )
      printf("repeat %i: %.4f %.4f\n", irep+1, repeat_array[irep].startpos, repeat_array[irep].endpos);
    //printf("mol %lli, nsite: %i  repeat: %.4f %.4f  remove indices: %i %i\n", inmap->id, inmap->numsite[icol], repeat_array[irep].startpos, repeat_array[irep].endpos, startidx, stopidx); //debug -- note: indices are after shifting, so don't print them
    if( startidx == 0 || stopidx == 0 ) { 
      printf("ERROR in maskRepeatLabels: startidx (%i) or stopidx (%i) not found in molecule %lli (%i sites, last site %.4f)--skipping: repeat start=%.4f end=%.4f\n", startidx, stopidx, inmap->id, inmap->numsite[icol], inmap->site[icol][inmap->numsite[icol]], repeat_array[irep].startpos, repeat_array[irep].endpos);
      return; 
    }
    int delta = stopidx-startidx+1; //quantity to shift by, also n sites removed--since stopidx is removed, add 1
    for( int i=startidx; i <= inmap->numsite[icol]+1; i++ ) { //numsite+1 is length--need to copy this also
      if( i <= inmap->numsite[icol]+1-delta ) //copy all sites plus length
	inmap->site[icol][i] = inmap->site[icol][i+delta];
      else
	inmap->site[icol][i] = 0; //null elements past the end
    }
    inmap->numsite[icol] -= delta; //numsites reduced
    //for( int i=1; i<=inmap->numsite[icol]; i++) printf("%i: %.4f\n", i, inmap->site[icol][i]);
  } //end for repeat_array elements
} //end maskRepeatLabels


/* Simple repeat detection: if neighboring intervals are similar, they are repeat.
   Similar means that abs(interval1-interval2) <= tolerance * sqrt(SF*SF + SR*SR*X*X), where X = (interval1+interval2)/2

   Shortest possible number of repeat intervals is 2 (minele), use more to find longer repeats.
   inlist must be sorted.
   Return number of repeats found : for each repeat :
   first ele is start position, second is end position, third is number of repeat elements, fourth is average size of elements.
*/
int findRepeatSimple(double* inlist, int inlist_len, double tolerance, int minele, bool verbose) {
  //interval i starts at inlist[i] and ends at inlist[i+1] (obviously one fewer element)
  double *intervals = (double *)malloc((inlist_len-1)*sizeof(double));
  if(!intervals){
    printf("findRepeatSimple:malloc of %lu bytes failed: inlist_len=%d\n",(inlist_len-1)*sizeof(double),inlist_len);
    fflush(stdout);exit(1);
  }

  //  double intervals[inlist_len-1];
  for(int i=0; i < inlist_len-1; i++) {
    intervals[i] = inlist[i+1] - inlist[i];
    if(DEBUG>=2) assert(0 <= inlist[i] && inlist[i] <= inlist[i+1]);
  }
  //repeats i are in interval[i] and interval[i+1], which means inlist[i] to inlist[i+2] (no tolerance applied yet)
  double *repeats = (double *)malloc((inlist_len-2)*sizeof(double));
  if(!repeats){
    printf("findRepeatSimple:malloc of %lu bytes failed:inlist_len=%d\n",(inlist_len-2)*sizeof(double),inlist_len);
    fflush(stdout);exit(1);
  }

  //  double repeats[inlist_len-2];
  bool isrepeat = false;
  for(int i=0; i < inlist_len-2; i++) {
    // WAS136 : repeats[i] = fabs(intervals[i+1] - intervals[i])/min(intervals[i+1],intervals[i]);
    repeats[i] = fabs(intervals[i+1] - intervals[i]);
    double avg = (intervals[i+1] + intervals[i]) * 0.5 * 0.001;
    double sd = sqrt(SF[0]*SF[0] + SR[0]*SR[0]*avg*avg) * 1000.0;
    if(verbose){
      printf("findRepeatSimple: i=%d: intervals[i]= %0.1f,intervals[i+1]= %0.1f, repeats[i]= %0.6f (sd= %0.1f, tolerance=%0.1f)\n",i,intervals[i],intervals[i+1],repeats[i],sd, tolerance * sd);
      fflush(stdout);
    }
    if (repeats[i] > tolerance * sd) 
      repeats[i] = -repeats[i];
    else
      isrepeat = true;
  }
  if( !isrepeat ) {
    free(repeats);
    free(intervals);
    return 0;
  }

  //in order to get number of eles, need number of sequential repeats--if it's unbroken, it's the same element
  minele = max(1,minele-1); //take one away from argument (minimum 1) because it's applied to number of eles in repeats, which is number of intervals less one
  isrepeat = false; //reuse for different purpose
  double startpos = -1.0;
  int nele = 0, startidx = 0;
  int nrepeat = 0; //number of new repeat objects created -- return this
  int i=0;
  for( ; i < inlist_len-2; i++ ) {
    if(verbose){
      printf("i=%d/%d: inlist[i]= %0.4f, intervals[i,i+1]= %0.4f, %0.4f : repeats= %0.4f (nele=%d, isrepeat=%d, minele=%d)\n", 
	     i,inlist_len-2, inlist[i], intervals[i], intervals[i+1], repeats[i], nele,isrepeat,minele);
      fflush(stdout);
    }
    //repeats is -1 for non-repeat intervals; >= is very important, since exact match is 0
    if( repeats[i] >= 0.0 ) {
      nele++;
      if( !isrepeat ) {
	startidx = i;
	startpos = inlist[i];
      }
      isrepeat = true;
    }
    else if( isrepeat ) { //rep < 0 : store current repeat, reset isrepeat, startpos, nele
      assert( startpos >= 0.0); //must be >= 0
      if( nele > minele ) {
	double avgint = 0; //get average element size
	for( int j=startidx; j < i+1; j++ )
	  avgint += intervals[j];
	avgint /= (i+1-startidx);
	//create new repeat object
	if( numrepeats == maxrepeats ) 
	  maxrepeats = growRepeatArray(repeat_array, numrepeats);
	//nele is n elements in repeat, but n intervals is actually one more
	//and since this is the interval after the repeat ended, we already added 1--only add one more for inlist
	repeat_array[numrepeats] = repeat(startpos, inlist[i+1], nele+1, avgint);
	numrepeats++;
	nrepeat++;
      }
      isrepeat = false;
      startpos = -1.0;
      startidx = nele = 0;
    }
  } //end for eles in repeats

  //if the end of the list is a repeat element, the else if isrepeat isn't hit--got to catch this
  if( isrepeat && nele > minele ) {
    assert( startpos >= 0.0); //dummy check
    double avgint = 0; //get average element size
    for( int j=startidx; j<i+1; j++ )
      avgint += intervals[j];
    avgint /= (i+1-startidx);
    //create new repeat object -- copy logic in loop above
    if( numrepeats == maxrepeats ) 
      maxrepeats = growRepeatArray(repeat_array, numrepeats);
    repeat_array[numrepeats] = repeat(startpos, inlist[i+1], nele+1, avgint);
    numrepeats++;
    nrepeat++;
  }
  free(repeats);
  free(intervals);
  return nrepeat;
}
//end findRepeatSimple


//this is the similarity criteria for finding a simple repeat
bool similarInterval(double a, double b, float tolerance) {
  return fabs(a-b)/(a+b) < tolerance/2;
}


//given an array, it's size, n elements used, and element to add,
// check if new element will fit, and if so, add it
// if not, grow the size, then add it
// return the maxsize, new size if changed, passed maxsize if not
template<class T>
int addToArray(T*& array, int maxsize, int nele, T toadd) {
  assert(nele <= maxsize); //arguments are wrong
  if( nele < maxsize ) {
    array[nele] = toadd;
    return maxsize;
  }
  maxsize = growArray(array, maxsize); //see above
  array[nele] = toadd;
  return maxsize;
} //end addToArray


//store repeat data in repeat object if number of repeat elements
// (from repUnitLengths_num, not repFreq, because that can be > the former)
// is at least minele
//needs the start and end pos,
// repUnitLengths array, repFreq (length of repUnitLengths), FP/FN arrays.
//also need to null all of these, so this function modifies all of its
//  parameters passed by reference, but does not deallocate the memory
//need to return whether element added so that findSimpleRepeatFPFN knows
int endRepeat(double& startpos, double& endpos, double*& repUnitLengths, int& repUnitLengths_num, int*& fpIndices, int& fpIndices_num, int*& fnIndices, int& fnIndices_num, int minele) {
  //printf("  endRepeat: %i  %.4f %.4f, %i %i\n", repUnitLengths_num+fnIndices_num, startpos, endpos, fpIndices_num, fnIndices_num); //debug
  int add = 0;
  if( repUnitLengths_num >= minele ){
    if( numrepeats == maxrepeats ) 
      maxrepeats = growRepeatArray(repeat_array, numrepeats);
    //  repeat(double start, double end, double* intervals, int nintervals, int* fpind, int nfpind, int* fnind, int nfnind) {
    repeat_array[numrepeats] = repeat(startpos, endpos, repUnitLengths, repUnitLengths_num, fpIndices, fpIndices_num, fnIndices, fnIndices_num);
    numrepeats++;
    add++;
  }
  startpos = 0; //endpos = 0; //endpos is NOT copied, so cannot null it (in one case it is, but it's temporary anyway)
  for(int i=0; i<repUnitLengths_num; i++)
    repUnitLengths[i] = 0;
  repUnitLengths_num = 0;
  for(int i=0; i<fpIndices_num; i++)
    fpIndices[i] = 0;
  fpIndices_num = 0;
  for(int i=0; i<fnIndices_num; i++)
    fnIndices[i] = 0;
  fnIndices_num = 0;
  return add;
}


int findSimpleRepeatFPFN(double* labelPosList, int labelPosList_len, float tolerance, int minele ) {
  bool debug = false; //true
  //interval i starts at inlist[i] and ends at inlist[i+1] (obviously one fewer element)
  double *distList = (double *)malloc((labelPosList_len-1)*sizeof(double));
  for(int i=0; i < labelPosList_len-1; i++) 
    distList[i] = labelPosList[i+1] - labelPosList[i];

  int repUnitLengths_len = 10; //size of repUnitLengths array
  int repUnitLengths_num = 0; //number elements of repUnitLengths array used
  //number of repeat elements == repUnitLengths_num + fnIndices_num because FN are two intervals, but only one get stored in repUnitLengths
  double *repUnitLengths = (double *)malloc(repUnitLengths_len*sizeof(double));
  int nrepeat = 0; //return number of repeats added to repeat_array

  int  fpIndices_len = 1; //change back to 5
  int  fpIndices_num = 0;
  int* fpIndices = (int*)malloc(fpIndices_len*sizeof(double));
  int  fnIndices_len = 1; //change back to 5
  int  fnIndices_num = 0;
  int* fnIndices = (int*)malloc(fnIndices_len*sizeof(double));

  bool insideRepeat = false, falsePosFlag = false, falseNegFlag = false;
  double repeatStart = 0; //like startpos in findRepeatSimple
  double L = distList[0]; //left interval
  double M = distList[1]; //middle interval
  double R = 0; //distList[2]; //assigned inside loop below -- right interval
  double falsePosDistHolder = 0.; //stores distance to the left of an FP

  for( int d=2; d<labelPosList_len-1; d++ ) {
    R = distList[d];
    if( falsePosFlag )
      R += falsePosDistHolder;
    if(debug) printf("d=%2i: L=%.4f M=%.4f R=%.4f; fpflag=%i fnflag=%i\n", d, L, M, R, falsePosFlag, falseNegFlag);
    
    if( !insideRepeat ) { //check for start of repeat
      if( similarInterval(L, M, tolerance) && similarInterval(M, R, tolerance) ) {
	if(debug) printf("  !insideRepeat, similar interval\n");
	//store the first three intervals in the repUnitLengths array
	repUnitLengths_len = addToArray(repUnitLengths, repUnitLengths_len, repUnitLengths_num, L);
	repUnitLengths_num++;
	repUnitLengths_len = addToArray(repUnitLengths, repUnitLengths_len, repUnitLengths_num, M);
	repUnitLengths_num++;
	repUnitLengths_len = addToArray(repUnitLengths, repUnitLengths_len, repUnitLengths_num, R);
	repUnitLengths_num++;
	insideRepeat = true;
	repeatStart = labelPosList[d-2];
      }
    }
    else { //inside Repeat
      if( similarInterval(M, R, tolerance) ) { //continue repeat
	if(debug) printf("  insideRepeat, similar interval\n");
	if( falseNegFlag ) {
	  //if( falsePosFlag ) return -1; //debug
	  assert(!falsePosFlag); //if FN flag, then FP flag must be false
	  //store the element: M is already half the true distance of this interval, so store that
	  repUnitLengths_len = addToArray(repUnitLengths, repUnitLengths_len, repUnitLengths_num, M);
	  repUnitLengths_num++;
	  fnIndices_len = addToArray(fnIndices, fnIndices_len, fnIndices_num, d); //also store FN (prev label: -1+1=0)
	  fnIndices_num++;
	  falseNegFlag = false;
	}
	repUnitLengths_len = addToArray(repUnitLengths, repUnitLengths_len, repUnitLengths_num, R);
	repUnitLengths_num++;
	if( falsePosFlag ) { //need to store it
	  fpIndices_len = addToArray(fpIndices, fpIndices_len, fpIndices_num, d+1); //+1 bc start at 1, not 0
	  fpIndices_num++;
	  falsePosFlag = false; //if M and R are similar, this cannot be a FP interval (and above case has assert !this)
	}
      }
      else if( M > R ) {
	if( falsePosFlag ) { //reset R for next interval
	  if( debug ) {
	    printf("  reset R: old R=%.4f, new R=%.4f\n", R, distList[d]);
	    printf("  reset M: old M=%.4f, new M=%.4f\n", M, distList[d-1]);
	  }
	  R = distList[d];
	  M = distList[d-1]; //assign this so that L becomes the interval which is too short, so must skip this before next start
	}
	if( falsePosFlag || falseNegFlag ) { //two elements in a row are not similar--end repeat
	  nrepeat += endRepeat(repeatStart, labelPosList[d-1], repUnitLengths, repUnitLengths_num, fpIndices, fpIndices_num, fnIndices, fnIndices_num, minele);
	  insideRepeat = falsePosFlag = falseNegFlag = false;
	}
	else {
	  falsePosFlag = true;
	  falsePosDistHolder = R;
	}
      }
      else if( M < R ) { //this condition is redundant, but more readable
	//this is the condition for a FN (if already FP or FN, end)
	if( similarInterval(R/2, M, tolerance) && !falsePosFlag && !falseNegFlag ) { 
	  falseNegFlag = true;
	  R /= 2; //very clever
	}
	else {
	  if( falsePosFlag ) { //reset R for next interval
	    if( debug ) {
	      printf("  reset R: old R=%.4f, new R=%.4f\n", R, distList[d]);
	      printf("  reset M: old M=%.4f, new M=%.4f\n", M, distList[d-1]);
	    }
	    R = distList[d];
	    M = distList[d-1]; //assign this so that L becomes the interval which is too short, so must skip this before next start
	  }
	  double endpos = ( falsePosFlag || falseNegFlag ? labelPosList[d-1] : labelPosList[d] ); 
	  nrepeat += endRepeat(repeatStart, endpos, repUnitLengths, repUnitLengths_num, fpIndices, fpIndices_num, fnIndices, fnIndices_num, minele);
	  insideRepeat = falsePosFlag = falseNegFlag = false;
	}
      }
    } //end inside repeat
    if( !falsePosFlag ) { //if FP, L and M remain, and R gets FP dist added on next iteration
      L = M;
      M = R;
    }
  } //end loop over intervals

  //need to catch case where last label is end of repeat
  nrepeat += endRepeat(repeatStart, labelPosList[labelPosList_len-1], repUnitLengths, repUnitLengths_num, fpIndices, fpIndices_num, fnIndices, fnIndices_num, minele);

  free(distList);
  free(repUnitLengths);
  return nrepeat;
} //findSimpleRepeatFPFN


//call fn declared in simpleRepeat.h: findRepeatSimple
/* command line values from -svRepeat :
   tolerance * sqrt(sf*sf + sr*sr*X*X) is max abs error (X = average of 2 intervals) 
   minele is minimum number of repeats
   confPenalty : If > 0 compute adjusted confidence and verify it is below -T threshold. */

void xmapEntry::callRepeat(bool verbose, double tolerance, int minele, double confPenalty) {
  if( refsites != NULL && numrefsites != 0 && refendidx-refstartidx >= 1){
    if(DEBUG) assert(numrefsites >= refendidx-refstartidx+1);
    int nrepeats = findRepeatSimple(refsites, refendidx-refstartidx+1/* WAS130 numrefsites*/, tolerance, minele, /* (qrycontigid==2830 && refcontigid==10)*/ false);
    if( nrepeats > 0 ) { 
      have_ref_repeat = true;
      //print the repeats found -- the elements of repeat_array to print are last nrepeats
      if(/* qrycontigid==51 && refcontigid==5 */  verbose ){
	for( int i=numrepeats-nrepeats; i < numrepeats; i++ )
	  printf("ref repeat %2i: xmapid=%2i:  refid=%lld,qryid=%lld  startpos=%.1f  endpos=%.1f  nele=%2i  avgsize=%.1f\n", 
		 i+1, xmapid, refcontigid, qrycontigid, repeat_array[i].startpos, repeat_array[i].endpos, repeat_array[i].nele, repeat_array[i].avgsize);
	fflush(stdout);
      }
      if(confPenalty > 0.0){/* compute adjusted confidence of matchgroup and make sure it drops below LogPvThreshold */
	int np = align->numpairs;
	int repeatIntervals = 0;
	for(int i = numrepeats - nrepeats; i < numrepeats; i++)
	  repeatIntervals += repeat_array[i].nele;

	double adjconf = confidence * max(0.0,np - repeatIntervals - 1.0)/(np-1.0)   - confPenalty * nrepeats;
	if(adjconf >= LogPvThreshold){
	  if(VERB>=2 && qrycontigid==51 && refcontigid==5){
	    printf("Invalid ref repeat found: xmapid=%d:refid=%lld,qryid=%lld: numpairs=%d,repeats=%d,repeatIntervals=%d,confidence = %0.2f -> %0.2f\n",
		   xmapid,refcontigid,qrycontigid,np,nrepeats,repeatIntervals,confidence,adjconf);
	    fflush(stdout);
	  }

	  have_ref_repeat = false;
	}
      }
    } else if(VERB>=2 && qrycontigid==51 && refcontigid==5){
      printf("No ref repeat found: xmapid=%2i: refid=%lld,qryid=%lld\n",xmapid,refcontigid,qrycontigid);
      fflush(stdout);
    }
  }

  if(QRY_REPEAT && numqrysites != 0 && qrysites != NULL && abs(qryendidx-qrystartidx) >= 1 && !have_ref_repeat){
    if(DEBUG) assert(numqrysites >= abs(qryendidx-qrystartidx)+1);
    int nrepeats = findRepeatSimple(qrysites, abs(qryendidx-qrystartidx)+1, tolerance, minele, /* (qrycontigid==2830 && refcontigid==10) */ false);
    if( nrepeats > 0 ) { 
      have_ref_repeat = true;
      //print the repeats found -- the elements of repeat_array to print are last nrepeats
      if(/* qrycontigid==51 && refcontigid==5*/verbose){
	for( int i=numrepeats-nrepeats; i < numrepeats; i++ )
	  printf("qry repeat %2i: xmapid=%2i:  refid=%lld,qryid=%lld  startpos=%.1f  endpos=%.1f  nele=%2i  avgsize=%.1f\n", 
		 i+1, xmapid, refcontigid, qrycontigid, repeat_array[i].startpos, repeat_array[i].endpos, repeat_array[i].nele, repeat_array[i].avgsize);
	fflush(stdout);
      }
      if(confPenalty > 0.0){/* compute adjusted confidence of matchgroup and make sure it drops below LogPvThreshold */
	int np = align->numpairs;
	int repeatIntervals = 0;
	for(int i = numrepeats - nrepeats; i < numrepeats; i++)
	  repeatIntervals += repeat_array[i].nele;

	double adjconf = confidence * max(0.0,np - repeatIntervals - 1.0)/(np-1.0)   - confPenalty * nrepeats;
	if(VERB>=2 && qrycontigid==51 && refcontigid==5){
	  printf("%s qry repeat found: xmapid=%d:refid=%lld,qryid=%lld: numpairs=%d,repeats=%d,repeatIntervals=%d,confidence = %0.2f -> %0.2f (T=%0.2f)\n",
		 adjconf >= LogPvThreshold ? "Invalid" : "Valid", xmapid,refcontigid,qrycontigid,np,nrepeats,repeatIntervals,confidence,adjconf,LogPvThreshold);
	  fflush(stdout);
	}
	if(adjconf >= LogPvThreshold)
	  have_ref_repeat = false;
      }
    } else if(VERB>=2 && qrycontigid==51 && refcontigid==5){
      printf("No qry repeat found: xmapid=%2i: refid=%lld,qryid=%lld\n",xmapid,refcontigid,qrycontigid);
      fflush(stdout);
    }
  }
}


//output repeat_array data -- basename will be output_prefix, which is global
void output_repeat(char *basename)
{
  if(strstr(basename,"/dev/null"))
    return;

  char filename[PATH_MAX];
  strcpy(filename,basename);
  int i = strlen(filename);
  strcpy(&filename[i],".rmap");

  if(checkFile(filename))
    return;

  if(VERB){
    printf("Generating %s\n",filename);
    fflush(stdout);
  }

  FILE *fp;
  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("failed to open file %s for writing smap file:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  printversion(fp); /* write out commandline */
  fprintf(fp,"# RMAP File Version:\t0.1\n"); 
  //header
  fprintf(fp,"#h Index\tMapID\tLength\tColor\tRepeatStart\tRepeatEnd\tAvgRepeatLength\tNRepeatUnits\tStandardDeviation\tConfidence\tRepeatString\tFPstring\tFNstring\n");
  fprintf(fp,"#f int \tint \tfloat \tint \tfloat \tfloat      \tfloat    \tint      \tfloat   \tfloat     \tstring  \tstring  \tstring \n");

  int color = (usecolor==2 ? 2 : 1); //usecolor should handle this

  for( int i=0; i<numrepeats; i++ ) {
    //note: length, startpos, endpos, size, stddev are still in kb here
    fprintf(fp, "%i\t%lli\t%.1f\t%i\t%.1f\t%.1f\t%.1f\t%i\t%.1f\t%.1f\t(%.1f;%.1f)\t(", i+1, repeat_array[i].id, repeat_array[i].length*1000., color, repeat_array[i].startpos*1000., repeat_array[i].endpos*1000., repeat_array[i].avgsize*1000., repeat_array[i].nele, repeat_array[i].stddev*1000., -1., repeat_array[i].avgsize, repeat_array[i].avgsize);
    for( int j=0; j<repeat_array[i].fpIndices_len; j++ ) {
      fprintf(fp, "%i", repeat_array[i].fpIndices[j]);
      if( j != repeat_array[i].fpIndices_len-1 )
	fprintf(fp, ",");	
    }
    fprintf(fp, ")\t(");
    for( int j=0; j<repeat_array[i].fnIndices_len; j++ ) {
      fprintf(fp, "%i", repeat_array[i].fnIndices[j]);
      if( j != repeat_array[i].fnIndices_len-1 )
	fprintf(fp, ",");	
    }
    fprintf(fp, ")\n");
  }

  FILEclose(fp);
} //end void output_repeat(char *basename)


//output statistics file -- return n mol w repeats for use in output_repeatAll
void output_repeatStats(char *basename) {
  if(strstr(basename,"/dev/null"))
    return;

  char filename[PATH_MAX];
  strcpy(filename,basename);
  int i = strlen(filename);
  strcpy(&filename[i],"_repeatStats.txt");

  if(checkFile(filename))
    return;

  if(VERB){
    printf("Generating %s\n",filename);
    fflush(stdout);
  }

  FILE *fp;
  if((fp = fopen(filename,"w"))==NULL){
    int eno = errno;
    char *err = strerror(eno);
    printf("failed to open file %s for writing smap file:errno=%d:%s\n",filename,eno,err);
    exit(1);
  }

  printversion(fp); /* write out commandline */
  //no column description header

  long long previd = 0;
  int nrepeatmol = 0; //need to filter mols with multiple repeats
  double total_repeat_length_kb = 0;
  for( int i=0; i<numrepeats; i++ ) {
    if( repeat_array[i].id != previd )
      nrepeatmol++;
    total_repeat_length_kb += repeat_array[i].endpos - repeat_array[i].startpos; //startpos and endpos are still in kb here
    previd = repeat_array[i].id;
  }

  fprintf(fp,"N molecules: %i\n", nummaps);
  fprintf(fp,"N molecules with repeats: %i\n", nrepeatmol);
  fprintf(fp,"%% molecules with repeats: %.3f\n", nrepeatmol/double(nummaps)*100.); //*100 for %
  fprintf(fp,"Total N repeats: %i\n", numrepeats);
  fprintf(fp,"Avg %% repeats/molecule: %.3f\n", numrepeats/double(nummaps)*100.); //*100 for %
  fprintf(fp,"Total bases (Mb): %.3f\n", total_map_length_kb/1000.); //convert kb to Mb
  fprintf(fp,"Total bases in repeats (Mb): %.3f\n", total_repeat_length_kb/1000.); //convert kb to Mb
  fprintf(fp,"%% bases in repeats: %.3f\n", total_repeat_length_kb/total_map_length_kb*100.); //*100 for %
  FILEclose(fp);
} //end void output_repeatStats(char *basename)


//call above output functions, also, optionally take array of maps and its size if outputting filtered bnx
void output_repeatAll(int filtertype, Cmap **partmaps, int numpartmaps) {

  output_repeat(output_prefix);
  output_repeatStats(output_prefix); //don't need return here

  //0 is no filter output, 1 is with repeat, 2 is no repeat, 3 is mask
  if( filtertype ) { 
    float filtpercent = numpartmaps/float(nummaps)*100;
    if( filtertype == 1 ) //always print a message about filter
      printf("Creating bnx with no repeat maps: %i maps kept (%.3f%%)\n", numpartmaps, filtpercent);
    else if( filtertype == 2 )
      printf("Creating bnx with only repeat maps: %i maps kept (%.3f%%)\n", numpartmaps, filtpercent);
    else if( filtertype == 3 )
      printf("Creating bnx with masked repeats: %i maps kept (%.3f%%)\n", numpartmaps, filtpercent);
    else { //this means that there is a bug in RefAligner.cpp bc this should be checked at command line input
      printf("ERROR: invalid value for filter type: %i\n", filtertype);
      return;
    }
    //for( int i=0; i<numpartmaps; i++ ) //debug
    //printf("map %lli, len=%.4f\n", partmaps[i]->id, partmaps[i]->site[0][partmaps[i]->numsite[0]+1]);

    //for map output, see this block from refalign.cpp around line 9641 -- it takes array of pointers to Cmap objects
    //PartMappedIdPrefix --this is just the output_prefix + _repeat or _norepeat as appropriate -- get from filtertype
    char PartMappedIdPrefix[PATH_MAX];
    const char* suf;
    if( filtertype == 1 )
      suf = "nonrepeat";
    else if( filtertype == 2 )
      suf = "repeat";
    else if( filtertype == 3 )
      suf = "repeatmask";
    else {
      printf("Unknown filtertype=%d\n",filtertype);
      exit(1);
    }
    sprintf(PartMappedIdPrefix, "%s_%s", output_prefix, suf);
    if(usecolor)
      colors = 2;
    if(usecolor==2){
      for(int i = 0; i < numpartmaps; i++)
	partmaps[i]->colorswap(usecolor);
      char *tmp = Nickase[0]; Nickase[0] = Nickase[1]; Nickase[1] = tmp;
    }
    //printf("output_repeatAll: maptype=%i MapType=%i CmapMergeBnx=%i\n", maptype, MapType, CmapMergeBnx);
    //this logic is from the -mapped argument--defaults to input argument
    //logic for -merge is to only use CmapMergeBnx, not maptype -- just leave it this way
    if(maptype==0 || MapType >= 0 || CmapMergeBnx)
      output_bnx(PartMappedIdPrefix,partmaps,0,numpartmaps-1, 1, NULL, NULL,-1);
    else
      output_cmap(PartMappedIdPrefix,partmaps,0,numpartmaps-1);
    if(usecolor==2){
      for(int i = 0; i < numpartmaps; i++)
	partmaps[i]->colorswap(usecolor);
      char *tmp = Nickase[0]; Nickase[0] = Nickase[1]; Nickase[1] = tmp;
    }
    if(usecolor)
      colors = 1;
  }
} //end output_repeatAll


float stdDev(double* inlist, int nin, double mean) {
  double stddev=0;
  for( int i=0; i<nin; i++ )
    stddev += (inlist[i]-mean)*(inlist[i]-mean);
  return sqrt(stddev/nin);
}
