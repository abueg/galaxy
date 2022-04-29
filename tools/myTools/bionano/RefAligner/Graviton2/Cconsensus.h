/** data structure to hold a Contig/map and the alignment 
   of each component Nanomap(or another Contig consensus) with the consensus map */

class Cconsensus {
  int numsite; ///< total number of sites (interleaved if there are multiple colors) 
  double *site; ///< site[i=0..numsite-1] is the location of site i 
  char *color; ///< color[i=0..numsite-1] is the color id of site i 
  char *siteFP; ///< siteFP[i=0..numsite-1] is a flag that is set if site i is a False Positive 
  int mapid; /**< If >=0, the index in nummaps[] of the map of which this contig is a copy 
	      * NOTE : Either (mapid>=0 AND nummaps == 0) OR (mapid < 0 AND nummaps > 0) */

  int *sitecnt; ///< sitecnt[i=0..numsite-1] is a count of the number of aligned maps aligned to site i (If siteTP[i]==0, then site[i] is the end of a map and NOT an actual site) 
  int *siteFN; ///< siteFN[i=0..numsite-1] is a count of the number of aligned maps spanning site i, but NOT aligned to site i 
  int *fragcnt; ///< fragcnt[i=0..numsite] is the count of the number of aligned map spanning the consensus interval site[i .. i+1] 
  double *fragsum; ///< fragsum[i=0..numsite] is the sum of the aligned map intervals that span the interval site[i .. i+1] 

  int nummaps; ///< total number of maps/contigs that make up this contig (0 if this contig is a map) 
   
   /* NOTE : The following fields are not defined if nummaps == 0 */

  Cconsensus *contig; ///< contig[i=0..nummaps-1] is the i'th component of this contig 
  int *flip; ///< flip[i=0..nummaps-1] is the orientation of the i'th component of this contig 
  int **sitemap; /**< sitemap[i=0..nummaps-1][j=0..contig[i].numsite+1] maps site j of map i to a site of this contig 
		  * all values lie in the range 0 ... numsite+1, where 0 refers to the contig's left end and 
		  * numsite+1 to the contig's right end (value can be -1 if a site is not present in the consensus : currently not utilized) */
};
