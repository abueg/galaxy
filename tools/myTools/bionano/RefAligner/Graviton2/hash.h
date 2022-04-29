#ifndef HASH_H
#define HASH_H

#include <omp.h>

#include <sys/types.h>

#include "constants.h"

static Ident Hash_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/hash.h 11310 2020-07-14 01:00:37Z tanantharaman $");


#define OFFSET_SCALE_FIX 0 /* WAS402 2 */ /* 0 : With -hashScaleDelta, offset2 is based in original unscaled probe map (scaling is only used for sizing errors)
			      > 0 : With -hashScaleDelta, offset2 is based on rescaled probe map and output offset2 is reverse scaled to reflect that scaling is applied to qry (mapid1) in output
			      1 : Original implementation does not handle -hashMultiMatch with -hashScaleDelta
			      2 : New implementation handles -hashMultiMatch with -hashScaleDelta
			   */

#define HASH_STREAM 3 /* Used by pairalign.cpp to support streamed hashtable IO (can avoid the IO associated with writing and reading an actual hashtable file)
			     2 : Use double buffering to hide disk IO in the parallel section of pairalign.cpp
			     3 : Use double buffering AND streamed merge sort to avoid disk IO and hide merge sort in the parallel section of pairalign.cpp or refalign.cpp */

#define OFFSET_SPREAD 2 /* Maximum offset mismatch allowance (as multiple of OffsetKB) : OffsetSpread (based on -hashoffset) is the actual value used */

#define HASHBUF_SIZ (16*1024) /* size of input buffer per thread (as multiple of CHashMatch) (256 Kbytes x 2) */

extern int GroupedScore;/* see -hashGrouped */

extern size_t maxMatchesPerRefid, maxMapsPerRefid;
extern short *besthashscore;/* If -hashbest, besthashscore[i=0..numrmaps-1] is the best hashscore for any match with mapid1 (output id2) == i. If pairwise also any match with mapid2 (output id1) == i*/

/** Match between two maps accumulated from hash hits of corresponding offset. This is also the binary format used for file output */
class CHashMatch {
 public:
  int id1;/* probe map mapid : same as mapid2 */
  int id2;/* insert map mapid : same as mapid1 */
  int offset;/**< Location of left end of map id2 relative to left end of map id1 (approximated in kb) */
  short hashscore;
  char orientation;
  char scaleID;/* NOTE : this value represents scaling of id2 (query) vs id1 (ref) unlike elsewhere in hashtable, where it is the other way */
  float MinLoc2;/* leftmost location (in kb, in correct orientation) of map id2 involved in match */
  float MaxLoc2;/* rightmost location (in kb, in correct orientation) of map id2 involved in match */
};

/* comparison function to sort CHashMatch array in increasing order of id1,id2,orientation,(-hashscore),offset */
static inline int CHashMatchIncId1(CHashMatch *p1, CHashMatch *p2)
{
  return (p1->id1 > p2->id1) ? 1 : (p1->id1 < p2->id1) ? -1 :
    (p1->id2 > p2->id2) ? 1 : (p1->id2 < p2->id2) ? -1 : 
    (p1->orientation > p2->orientation) ? 1 : (p1->orientation < p2->orientation) ? -1 :
    (p1->hashscore < p2->hashscore) ? 1 : (p1->hashscore > p2->hashscore) ? -1 :
    (p1->offset - p2->offset);
}


extern void hash_generate(Cmap *rmap[], int numrmaps, Cmap *qmap[], int numqmaps, int pairwise, char *prefix);
/*  NOTE : if rmap[] and qmap[] are not the same, then all ids in qmap[] must strictly smaller than all ids in rmap[] */

extern void hash_input(char *filename, size_t &matchcnt, CHashMatch* &matches);/* reads the entire file (not recommended for files over 100GB) */

extern size_t hash_open(FILE* &HashFP, char *hash_filename);
extern int hash_eof(FILE *HashFP);
extern size_t hash_read(FILE *fp, int txt, CHashMatch* matches, size_t matchmax, int &linecnt);/* reads next section of file (can be used for very large files) */
extern void hash_close(FILE *HashFP, int final);

#endif
