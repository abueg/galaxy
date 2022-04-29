#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "globals.h"
#include "parameters.h"
#include "Cmap.h"
#include "Ccontig.h"

#define DELTAX_BLOCK 1 /* new block based memory allocation for deltaX[],deltaY[],deltaR[] */

#define BIASCORRECT_FIX 1 /* limit shift of each site so it does get closer to its other neighbor than the original neighbor */

static Ident Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/Cmap.cpp 10779 2020-03-25 17:43:32Z tanantharaman $");

void Cmap::allfree()
{

    if(VERB >= 2 || (LEAKDEBUG>=2  && id >= 0))
      printf("allfree():this=x%p:id=%lld,colors=%d,blockmem=%d,site[0]=%p\n",this,id,colors,blockmem,site[0]);

    if(refname){
      free(refname);
      refname = 0;
    }
    if(name){
      free(name);
      name = 0;
    }
    for(int c = 0; c < colors; c++){
      if(SNRcnt[c]){
	if(SNRdist[c]){
	  for(int i = 0; i <= numsite[c]+1; i++){
	    if(SNRdist[c][i]){
	      //	      if(!blockmem)
	      if(blockmem >= 0){
		if(VERB>=3){
		  printf("Calling delete[] Gmap[%d](%p)->SNRdist[%d][%d]=%p(SNRcnt=%d) blockmem=%d\n",mapid,this,c,i,SNRdist[c][i],SNRcnt[c][i],blockmem);
		  fflush(stdout);
		}
		delete [] SNRdist[c][i];
	      } else if(VERB>=3){
		printf("Skipped calling delete[] Gmap[%d](%p)->SNRdist[%d][%d]=%p(SNRcnt=%d) due to blockmem=%d\n",mapid,this,c,i,SNRdist[c][i],SNRcnt[c][i],blockmem);
		fflush(stdout);
	      }
	      SNRdist[c][i] = 0;
	    }
	  }
	  //	      if(!blockmem)
	  if(blockmem >= 0)
	    delete [] SNRdist[c];
	  SNRdist[c] = 0;
	}
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] SNRcnt[c];
	SNRcnt[c] = 0;
	if(SNRgmean[c]){
	  //	      if(!blockmem)
	  if(blockmem >= 0)
	    delete [] SNRgmean[c];
	  SNRgmean[c] = 0;
	}
	if(lnSNRsd[c]){
	  //	      if(!blockmem)
	  if(blockmem >= 0)
	    delete [] lnSNRsd[c];
	  lnSNRsd[c] = 0;
	}
      }
      if(remap[c]){
	if(VERB >= 2)
	  printf("allfree():id=%lld,c=%d,deleting remap[c]\n",id,c);
	delete [] remap[c];
	remap[c] = 0;
      }
      if(NickCnt[c]){
	if(VERB >= 2)
	  printf("allfree():id=%lld,c=%d,deleting NickCnt[c]\n",id,c);
	delete [] NickCnt[c];
	NickCnt[c] = 0;
      }
      if(site[c]){
	if(!blockmem){
	  if(VERB>=2 || LEAKDEBUG>=2){
	    printf("allfree():id=%lld,c=%d,blockmem=%d:deleting site[%d]=%p\n",id,c,blockmem,c,site[c]);
	    fflush(stdout);
	  }
	  delete [] site[c];
	}
	site[c] = 0;
      }
      if(SNR[c]){
	if(!blockmem)
	  delete [] SNR[c];
	SNR[c] = NULL;
      }
      if(Intensity[c]){
	if(!blockmem)
	  delete [] Intensity[c];
	Intensity[c] = NULL;
      }
      if(stitch[c]){
	if(!blockmem)
	  delete [] stitch[c];
	stitch[c] = NULL;
      }
      if(stitchLocation[c]){
	delete [] stitchLocation[c];
	stitchLocation[c] = NULL;
      }
      if(PSFWidth[c]){
	if(!blockmem)
	  delete [] PSFWidth[c];
	PSFWidth[c] = NULL;
      }
      if(truesiteL[c]){
	if(!blockmem)
	  delete [] truesiteL[c];
	truesiteL[c] = NULL;
      }
      if(truesiteR[c]){
	if(!blockmem)
	  delete [] truesiteR[c];
	truesiteR[c] = NULL;
      }
      if(ImageFOV[c]){
	if(!blockmem)
	  delete [] ImageFOV[c];
	ImageFOV[c] = NULL;
      }
      if(ImageX[c]){
	if(!blockmem)
	  delete [] ImageX[c];
	ImageX[c] = NULL;
      }
      if(ImageY[c]){
	if(!blockmem)
	  delete [] ImageY[c];
	ImageY[c] = NULL;
      }

      if(rawsite[c]){
	if(!rblockmem)
	  delete [] rawsite[c];
	rawsite[c] = NULL;
      }

      if(siteSD[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] siteSD[c];
	siteSD[c] = NULL;
      }
      if(sitecov[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] sitecov[c];
	sitecov[c] = NULL;
      }
      if(sitecnt[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] sitecnt[c];
	sitecnt[c] = NULL;
      }
      if(ChimQuality[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] ChimQuality[c];
	ChimQuality[c] = NULL;
      }
      if(ChimNorm[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] ChimNorm[c];
	ChimNorm[c] = NULL;
      }
      if(SegDupL[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] SegDupL[c];
	SegDupL[c] = NULL;
      }
      if(SegDupR[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] SegDupR[c];
	SegDupR[c] = NULL;
      }
      if(FragileEndL[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] FragileEndL[c];
	FragileEndL[c] = NULL;
      }
      if(FragileEndR[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] FragileEndR[c];
	FragileEndR[c] = NULL;
      }
      if(OutlierFrac[c]){
	if(blockmem >= 0)
	  delete [] OutlierFrac[c];
	OutlierFrac[c] = NULL;
      }
      if(FragSd[c]){
	if(blockmem >= 0)
	  delete [] FragSd[c];
	FragSd[c] = NULL;
      }
      if(ExpSd[c]){
	if(blockmem >= 0)
	  delete [] ExpSd[c];
	ExpSd[c] = NULL;
      }
      if(FragCov[c]){
	if(blockmem >= 0)
	  delete [] FragCov[c];
	FragCov[c] = NULL;
      }
      if(FragChiSq[c]){
	if(blockmem >= 0)
	  delete [] FragChiSq[c];
	FragChiSq[c] = NULL;
      }

      if(Mask[c]){
	if(blockmem >= 0)
	  delete [] Mask[c];
	Mask[c] = NULL;
      }

      // NOTE : orig* arrays are usually block allocated in mapcorrect(), but not always (see input_vfixx.cpp and with PairMerge in output.cpp)
      if(origSNRcnt[c]){
	if(origSNRdist[c]){
	  //	  for(int i = 0; i <= orignumsite[c]+1; i++){
	  //	    if(origSNRdist[c][i]){
	  //	      delete [] origSNRdist[c][i];
	  //	      origSNRdist[c][i] = 0;
	  //	    }
	  //	  }
	  //	  delete [] origSNRdist[c];
	  origSNRdist[c] = 0;
	}
	//	delete [] origSNRcnt[c];
	origSNRcnt[c] = 0;

	if(origSNRgmean[c]){
	  //	  delete [] origSNRgmean[c];
	  origSNRgmean[c] = 0;
	}
	if(origlnSNRsd[c]){
	  //	  delete [] origlnSNRsd[c];
	  origlnSNRsd[c] = 0;
	}
      }
      if(origsite[c]){
	//	delete [] origsite[c];
	origsite[c] = 0;
      }
      if(origsiteSD[c]){
	//	delete [] origsiteSD[c];
	origsiteSD[c] = 0;
      }
      if(origsitecov[c]){
	//	delete [] origsitecov[c];
	origsitecov[c] = 0;
      }
      if(origsitecnt[c]){
	//	delete [] origsitecnt[c];
	origsitecnt[c] = 0;
      }
      if(origChimQuality[c]){
	//	delete [] origChimQuality[c];
	origChimQuality[c] = 0;
      }
      if(origChimNorm[c]){
	//	delete [] origChimNorm[c];
	origChimNorm[c] = 0;
      }
      if(origSegDupL[c]){
	//	delete [] origSegDupL[c];
	origSegDupL[c] = 0;
      }
      if(origSegDupR[c]){
	//	delete [] origSegDupR[c];
	origSegDupR[c] = 0;
      }
      if(origFragileEndL[c]){
	//	delete [] origFragileEndL[c];
	origFragileEndL[c] = 0;
      }
      if(origFragileEndR[c]){
	//	delete [] origFragileEndR[c];
	origFragileEndR[c] = 0;
      }
      if(origOutlierFrac[c]){
	//	delete [] origOutlierFrac[c];
	origOutlierFrac[c] = 0;
      }
      if(origFragSd[c]){
	origFragSd[c] = 0;
      }
      if(origExpSd[c]){
	origExpSd[c] = 0;
      }
      if(origFragCov[c]){
	origFragCov[c] = 0;
      }
      if(origFragChiSq[c]){
	origFragChiSq[c] = 0;
      }
      
      if(origMask[c]){
	//	delete [] origMask[c];
	origMask[c] = 0;
      }

      if(origtruesiteL[c]){
	origtruesiteL[c] = NULL;
      }
      if(origtruesiteR[c]){
	origtruesiteR[c] = NULL;
      }
      if(origSNR[c]){
	origSNR[c] = NULL;
      }
      if(origIntensity[c]){
	origIntensity[c] = NULL;
      }
      if(origstitch[c]){
	origstitch[c] = NULL;
      }
      if(origPSFWidth[c]){
	origPSFWidth[c] = NULL;
      }
      if(origImageFOV[c]){
	origImageFOV[c] = NULL;
      }
      if(origImageX[c]){
	origImageX[c] = NULL;
      }
      if(origImageY[c]){
	origImageY[c] = NULL;
      }

      if(Nickase[c]){
	/*	printf("Cmap=%p:c=%d,Nickase[c]=%p\n",this,c,Nickase[c]);
		fflush(stdout);*/

	//	free(Nickase[c]);
	Nickase[c] = 0;
      }
      if(contig){
	delete contig;
	contig = NULL;
      }

      numsite[c]=0;
      orignumsite[c] = 0;
    }

    if(DEBUG) RunIndex = -1;
    /* NOTE : cannot reset mapid since it may be set in a critical section before calling allfree() */

    blockmem = rblockmem = 0;

    flipped = 0;
    origmap = NULL;
    OriginalMoleculeId = -1;
    ScanDirection = -1;
    UniqueScanId = -1;
    if(DEBUG>=2)
      ScanNumber = -1;
    startloc = endloc = 0.0;
}

/* free memory of all except 1st color */
void Cmap::allfree2()
{
    if(VERB >= 2)
      printf("allfree2():id=%lld,colors=%d\n",id,colors);

    for(int c = 1; c < colors; c++){
      if(SNRcnt[c]){
	if(SNRdist[c]){
	  for(int i = 0; i <= numsite[c]+1; i++){
	    if(SNRdist[c][i]){
	      delete [] SNRdist[c][i];
	      SNRdist[c][i] = 0;
	    }
	  }
	  if(blockmem >= 0)
	    delete [] SNRdist[c];
	  SNRdist[c] = 0;
	}
	if(blockmem >= 0)
	  delete [] SNRcnt[c];
	SNRcnt[c] = 0;
	if(SNRgmean[c]){
	  if(blockmem >= 0)
	    delete [] SNRgmean[c];
	  SNRgmean[c] = 0;
	}
	if(lnSNRsd[c]){
	  if(blockmem >= 0)
	    delete [] lnSNRsd[c];
	  lnSNRsd[c] = 0;
	}
      }
      numsite[c]=0;
      if(remap[c]){
	if(VERB >= 2)
	  printf("allfree2():id=%lld,c=%d,deleting remap[c]\n",id,c);
	delete [] remap[c];
	remap[c] = 0;
      }
      if(NickCnt[c]){
	if(VERB >= 2)
	  printf("allfree2():id=%lld,c=%d,deleting NickCnt[c]\n",id,c);
	delete [] NickCnt[c];
	NickCnt[c] = 0;
      }
      if(site[c]){
	if(!blockmem)
	  delete [] site[c];
	site[c] = 0;
      }
      if(SNR[c]){
	if(!blockmem)
	  delete [] SNR[c];
	SNR[c] = NULL;
      }
      if(Intensity[c]){
	if(!blockmem)
	  delete [] Intensity[c];
	Intensity[c] = NULL;
      }
      if(stitch[c]){
	if(!blockmem)
	  delete [] stitch[c];
	stitch[c] = NULL;
      }
      if(stitchLocation[c]){
	delete [] stitchLocation[c];
	stitchLocation[c] = NULL;
      }
      if(PSFWidth[c]){
	if(!blockmem)
	  delete [] PSFWidth[c];
	PSFWidth[c] = NULL;
      }
      if(truesiteL[c]){
	if(!blockmem)
	  delete [] truesiteL[c];
	truesiteL[c] = NULL;
      }
      if(truesiteR[c]){
	if(!blockmem)
	  delete [] truesiteR[c];
	truesiteR[c] = NULL;
      }
      if(ImageFOV[c]){
	if(!blockmem)
	  delete [] ImageFOV[c];
	ImageFOV[c] = NULL;
      }
      if(ImageX[c]){
	if(!blockmem)
	  delete [] ImageX[c];
	ImageX[c] = NULL;
      }
      if(ImageY[c]){
	if(!blockmem)
	  delete [] ImageY[c];
	ImageY[c] = NULL;
      }

      if(rawsite[c]){
	if(!rblockmem)
	  delete [] rawsite[c];
	rawsite[c] = NULL;
      }

      if(siteSD[c]){
	if(blockmem >= 0)
	  delete [] siteSD[c];
	siteSD[c] = 0;
      }
      if(sitecov[c]){
	if(blockmem >= 0)
	  delete [] sitecov[c];
	sitecov[c] = 0;
      }
      if(sitecnt[c]){
	if(blockmem >= 0)
	  delete [] sitecnt[c];
	sitecnt[c] = 0;
      }
      if(ChimQuality[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] ChimQuality[c];
	ChimQuality[c] = NULL;
      }
      if(ChimNorm[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] ChimNorm[c];
	ChimNorm[c] = NULL;
      }
      if(SegDupL[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] SegDupL[c];
	SegDupL[c] = NULL;
      }
      if(SegDupR[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] SegDupR[c];
	SegDupR[c] = NULL;
      }
      if(FragileEndL[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] FragileEndL[c];
	FragileEndL[c] = NULL;
      }
      if(FragileEndR[c]){
	//	      if(!blockmem)
	if(blockmem >= 0)
	  delete [] FragileEndR[c];
	FragileEndR[c] = NULL;
      }

      if(OutlierFrac[c]){
	if(blockmem >= 0)
	  delete [] OutlierFrac[c];
	OutlierFrac[c] = NULL;
      }
      if(FragSd[c]){
	if(blockmem >= 0)
	  delete [] FragSd[c];
	FragSd[c] = NULL;
      }
      if(ExpSd[c]){
	if(blockmem >= 0)
	  delete [] ExpSd[c];
	ExpSd[c] = NULL;
      }
      if(FragCov[c]){
	if(blockmem >= 0)
	  delete [] FragCov[c];
	FragCov[c] = NULL;
      }
      if(FragChiSq[c]){
	if(blockmem >= 0)
	  delete [] FragChiSq[c];
	FragChiSq[c] = NULL;
      }

      if(Mask[c]){
	delete [] Mask[c];
	Mask[c] = NULL;
      }

      // NOTE : orig* arrays are usually block allocated in mapcorrect(), but not always (see input_vfixx.cpp and with PairMerge in output.cpp)
      if(origSNRcnt[c]){
	if(origSNRdist[c]){
	  for(int i = 0; i <= orignumsite[c]+1; i++){
	    if(origSNRdist[c][i]){
	      //	      delete [] origSNRdist[c][i];
	      origSNRdist[c][i] = 0;
	    }
	  }
	  //	  delete [] origSNRdist[c];
	  origSNRdist[c] = 0;
	}
	//	delete [] origSNRcnt[c];
	origSNRcnt[c] = 0;

	if(origSNRgmean[c]){
	  //	  delete [] origSNRgmean[c];
	  origSNRgmean[c] = 0;
	}
	if(origlnSNRsd[c]){
	  //	  delete [] origlnSNRsd[c];
	  origlnSNRsd[c] = 0;
	}
      }
      orignumsite[c] = 0;
      if(origsite[c]){
	//	delete [] origsite[c];
	origsite[c] = 0;
      }
      if(origsiteSD[c]){
	//	delete [] origsiteSD[c];
	origsiteSD[c] = 0;
      }
      if(origsitecov[c]){
	//	delete [] origsitecov[c];
	origsitecov[c] = 0;
      }
      if(origsitecnt[c]){
	//	delete [] origsitecnt[c];
	origsitecnt[c] = 0;
      }

      if(origChimQuality[c]){
	//	delete [] origChimQuality[c];
	origChimQuality[c] = 0;
      }
      if(origChimNorm[c]){
	//	delete [] origChimNorm[c];
	origChimNorm[c] = 0;
      }
      if(origSegDupL[c]){
	//	delete [] origSegDupL[c];
	origSegDupL[c] = 0;
      }
      if(origSegDupR[c]){
	//	delete [] origSegDupR[c];
	origSegDupR[c] = 0;
      }
      if(origFragileEndL[c]){
	//	delete [] origFragileEndL[c];
	origFragileEndL[c] = 0;
      }
      if(origFragileEndR[c]){
	//	delete [] origFragileEndR[c];
	origFragileEndR[c] = 0;
      }
      if(origOutlierFrac[c]){
	origOutlierFrac[c] = 0;
      }
      if(origFragSd[c]){
	origFragSd[c] = 0;
      }
      if(origExpSd[c]){
	origExpSd[c] = 0;
      }
      if(origFragCov[c]){
	origFragCov[c] = 0;
      }
      if(origFragChiSq[c]){
	origFragChiSq[c] = 0;
      }
      
      if(origMask[c]){
	delete [] origMask[c];
	origMask[c] = 0;
      }

      if(origtruesiteL[c]){
	origtruesiteL[c] = NULL;
      }
      if(origtruesiteR[c]){
	origtruesiteR[c] = NULL;
      }
      if(origSNR[c]){
	origSNR[c] = NULL;
      }
      if(origIntensity[c]){
	origIntensity[c] = NULL;
      }
      if(origstitch[c]){
	origstitch[c] = NULL;
      }
      if(origPSFWidth[c]){
	origPSFWidth[c] = NULL;
      }
      if(origImageFOV[c]){
	origImageFOV[c] = NULL;
      }
      if(origImageX[c]){
	origImageX[c] = NULL;
      }
      if(origImageY[c]){
	origImageY[c] = NULL;
      }

      if(Nickase[c]){
	/*	printf("Cmap=%p:c=%d,Nickase[c]=%p\n",this,c,Nickase[c]);
		fflush(stdout);*/

	//	free(Nickase[c]);
	Nickase[c] = 0;
      }
    }
}


Cmap::Cmap(Cmap *pmap)
{
  if(DEBUG && !(pmap->siteSD[0] && pmap->sitecov[0] && pmap->sitecnt[0])){
    printf("Contig copy constructor is only implemented for Contig maps (from CMAP input): this map (id=%lld) is from %s\n",
	   pmap->id, vfixx_filename[pmap->fileid]);
    fflush(stdout);exit(1);
  }

  init();/* default initialization of map */

  mapid = pmap->mapid;
  id = pmap->id;
  fileid = pmap->fileid;
  ScanNumber = pmap->ScanNumber;
  UniqueScanId = pmap->UniqueScanId;

  for(int c = 0; c < colors; c++){
    if(DEBUG>=2) assert(pmap->siteSD[c] && pmap->sitecov[c] && pmap->sitecnt[c]);

    numsite[c] = pmap->numsite[c];
    site[c] = new FLOAT[numsite[c] + 2];
    blockmem = 0;
    site[c][0] = 0.0;
    for(int i = 1; i <= numsite[c] + 1; i++)
      site[c][i] = pmap->site[c][i];

    if(pmap->siteSD[c]){
      siteSD[c] = new double[numsite[c]+2];
      sitecov[c] = new float[numsite[c]+2];
      sitecnt[c] = new float[numsite[c]+2];
      for(int i = 1; i <= numsite[c]+1; i++){
	siteSD[c][i] = pmap->siteSD[c][i];
	sitecov[c][i] = pmap->sitecov[c][i];
	sitecnt[c][i] = pmap->sitecnt[c][i];
      }
    }
    if(pmap->SNRcnt[c]){
      SNRcnt[c] = new int[numsite[c]+2];
      SNRgmean[c] = new double[numsite[c]+2];
      lnSNRsd[c] = new double[numsite[c]+2];
      SNRdist[c] = new double*[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++){
	SNRcnt[c][J] = pmap->SNRcnt[c][J];
	SNRgmean[c][J] = pmap->SNRgmean[c][J];
	lnSNRsd[c][J] = pmap->lnSNRsd[c][J];
	SNRdist[c][J] = NULL;
	if(SNRcnt[c][J] > 0){
	  SNRdist[c][J] = new double[SNRcnt[c][J]];
	  if(DEBUG) assert(pmap->SNRdist[c][J] != NULL);
	  for(int t = 0; t < SNRcnt[c][J]; t++)
	    SNRdist[c][J][t] = pmap->SNRdist[c][J][t];
	}
      }
    }
    if(pmap->ChimQuality[c]){
      ChimQuality[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        ChimQuality[c][J] = pmap->ChimQuality[c][J];
    }
    if(pmap->ChimNorm[c]){
      ChimNorm[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        ChimNorm[c][J] = pmap->ChimNorm[c][J];
    }
    if(pmap->SegDupL[c]){
      SegDupL[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        SegDupL[c][J] = pmap->SegDupL[c][J];
    }
    if(pmap->SegDupR[c]){
      SegDupR[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        SegDupR[c][J] = pmap->SegDupR[c][J];
    }
    if(pmap->FragileEndL[c]){
      FragileEndL[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        FragileEndL[c][J] = pmap->FragileEndL[c][J];
    }
    if(pmap->FragileEndR[c]){
      FragileEndR[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        FragileEndR[c][J] = pmap->FragileEndR[c][J];
    }
    if(pmap->OutlierFrac[c]){
      OutlierFrac[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        OutlierFrac[c][J] = pmap->OutlierFrac[c][J];
    }
    if(pmap->FragSd[c]){
      FragSd[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        FragSd[c][J] = pmap->FragSd[c][J];
    }
    if(pmap->ExpSd[c]){
      ExpSd[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        ExpSd[c][J] = pmap->ExpSd[c][J];
    }
    if(pmap->FragCov[c]){
      FragCov[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        FragCov[c][J] = pmap->FragCov[c][J];
    }
    if(pmap->FragChiSq[c]){
      FragChiSq[c] = new float[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        FragChiSq[c][J] = pmap->FragChiSq[c][J];
    }

    if(pmap->Mask[c]){
      Mask[c] = new size_t[numsite[c]+2];
      for(int J = 0; J <= numsite[c]+1; J++)
        Mask[c][J] = pmap->Mask[c][J];
    }

    if(DEBUG>=2) assert(siteSD[c] && sitecov[c] && sitecnt[c]);
  }
}

Cmap::Cmap(Ccontig *pcontig)
{
  init();/* default initialization of map */
  
  if(colors > 1){
    printf("Cmap::Cmap(Ccontig*): more than 1 color not supported\n");
    exit(1);
  }

  if(pcontig->mapid >= 0){/* just copy Gmap[mapid] from site trimL to trimR */
    Cmap *pmap = Gmap[pcontig->mapid];
    int M = pmap->numsite[0];
    int left = pcontig->trimL[0];
    int right = pcontig->trimR[0];
    if(DEBUG) assert(left >= 0 && right <= M+1);
    numsite[0] = right - left -1;
    site[0] = new FLOAT[numsite[0] + 2];
    blockmem = 0;
    site[0][0] = 0.0;
    for(int i = 1; i <= numsite[0]; i++)
      site[0][i] = pmap->site[0][i+left] - pmap->site[0][left];
    site[0][numsite[0]+1] = pmap->site[0][numsite[0]+1+left] - pmap->site[0][left];

    return;
  }

  /* copy the top level consensus map */
  numsite[0] = pcontig->numsite[0];
  site[0] = new double[numsite[0]+2];
  blockmem = 0;
  for(int i = 0; i <= numsite[0]+1; i++)
    site[0][i] = pcontig->site[0][i];
}

// HERE HERE : add support for Haplotype contig : already implemented for simple case with contig.nummaps == -1 (as should be the case right after hap map input if using -merge)

void Cmap::trim(double RangeLeft, double RangeRight)
{
  if(!(siteSD[0] && sitecov[0] && sitecnt[0])){
    printf("trim() only implemented for Contig maps (from CMAP input with StdDev, Coverage and Occurence): this map (id=%lld) is from %s\n",
	   id, vfixx_filename[fileid]);
    printf("\t siteSD[0]=%p, sitecov[0]=%p, sitecnt[0]=%p\n",siteSD[0],sitecov[0],sitecnt[0]);
    fflush(stdout);exit(1);
  }

  double Len = site[0][numsite[0]+1];

  RangeLeft = max(ZERO, RangeLeft);
  RangeRight = min(Len, RangeRight);

  int left[MAXCOLOR],right[MAXCOLOR];

  for(int c = 0; c < colors; c++){
    left[c] = 0, right[c] = numsite[c]+1;
    while(left[c] < right[c] && site[c][left[c]] < RangeLeft - 0.0001)
      left[c]++;
    while(left[c] < right[c] && site[c][right[c]] > RangeRight + 0.0001)
      right[c]--;
    RangeLeft = min(RangeLeft, max(ZERO, site[c][left[c]] - ((FLOAT)MININTERVAL)));
    RangeRight = max(RangeRight, min(Len, site[c][right[c]] + ((FLOAT)MININTERVAL)));
  }

  for(int c = 0; c < colors; c++){
    if(left[c] > 0 || right[c] <= numsite[c]){
      left[c] = max(1,left[c]);
      right[c] = min(right[c],numsite[c]);
      for(int i = left[c]; i <= right[c]; i++){
	site[c][i-left[c]+1] = site[c][i] - RangeLeft;
	if(left[c] == 1)
	  continue;
	if(siteSD[c]){
	  siteSD[c][i-left[c]+1] = siteSD[c][i];
	  sitecov[c][i-left[c]+1] = sitecov[c][i];
	  sitecnt[c][i-left[c]+1] = sitecnt[c][i];
	}
	if(Mask[c])
	  Mask[c][i-left[c]+1] = Mask[c][i];
	if(ChimQuality[c])
	  ChimQuality[c][i-left[c]+1] = ChimQuality[c][i];
	if(ChimNorm[c])
	  ChimNorm[c][i-left[c]+1] = ChimNorm[c][i];
	if(SegDupL[c])
	  SegDupL[c][i-left[c]+1] = SegDupL[c][i];
	if(SegDupR[c])
	  SegDupR[c][i-left[c]+1] = SegDupR[c][i];
	if(FragileEndL[c])
	  FragileEndL[c][i-left[c]+1] = FragileEndL[c][i];
	if(FragileEndR[c])
	  FragileEndR[c][i-left[c]+1] = FragileEndR[c][i];
	if(OutlierFrac[c])
	  OutlierFrac[c][i-left[c]+1] = OutlierFrac[c][i];
	if(FragSd[c])
	  FragSd[c][i-left[c]+1] = FragSd[c][i];
	if(ExpSd[c])
	  ExpSd[c][i-left[c]+1] = ExpSd[c][i];
	if(FragCov[c])
	  FragCov[c][i-left[c]+1] = FragCov[c][i];
	if(FragChiSq[c])
	  FragChiSq[c][i-left[c]+1] = FragChiSq[c][i];
      }
      if(contig && contig->nummaps == -1 && contig->HapSite[c] && contig->HapDelta[c]){/* update Haplotype Map information */
	for(int i = left[c]; i < right[c]; i++){
	  contig->HapSite[c][i-left[c]+1] = contig->HapSite[c][i];
	  contig->HapDelta[c][i-left[c]+1] = contig->HapDelta[c][i];
	}
	contig->HapSite[c][right[c]-left[c]+1] = contig->HapSite[c][right[c]];
      }
      if(SNRcnt[c] && left[c] != 1){
	for(int i = left[c]; i <= right[c]; i++){

	  SNRcnt[c][i-left[c]+1] = SNRcnt[c][i];
	  SNRgmean[c][i-left[c]+1] = SNRgmean[c][i];
	  lnSNRsd[c][i-left[c]+1] = lnSNRsd[c][i];

	  if(SNRdist[c][i-left[c]+1])
	    delete [] SNRdist[c][i-left[c]+1];
	  SNRdist[c][i-left[c]+1] = SNRdist[c][i];
	  SNRdist[c][i] = NULL;// to mark SNRdist[c][i] as invalid, so it is not deleted later
	}
	SNRcnt[c][right[c]+1-left[c]+1] = 0;
      }
      if(Mask[c] && right[c] >= numsite[c] && RangeRight >= site[c][numsite[c]+1] - 1e-6)
	Mask[c][right[c]-left[c]+2] = Mask[c][numsite[c]+1];// right end
      numsite[c] = right[c]-left[c]+1;
      site[c][numsite[c]+1] = RangeRight - RangeLeft;
      if(siteSD[c]){
	siteSD[c][numsite[c]+1] = 0.0;
	sitecov[c][numsite[c]+1] = sitecov[c][right[c]+1];
	sitecnt[c][numsite[c]+1] = sitecnt[c][right[c]+1];
      }
      if(ChimQuality[c])
	ChimQuality[c][numsite[c]+1] = ChimQuality[c][right[c]+1];
      if(ChimNorm[c])
	ChimNorm[c][numsite[c]+1] = ChimNorm[c][right[c]+1];
      if(SegDupL[c])
	SegDupL[c][numsite[c]+1] = SegDupL[c][right[c]+1];
      if(SegDupR[c])
	SegDupR[c][numsite[c]+1] = SegDupR[c][right[c]+1];
      if(FragileEndL[c])
	FragileEndL[c][numsite[c]+1] = FragileEndL[c][right[c]+1];
      if(FragileEndR[c])
	FragileEndR[c][numsite[c]+1] = FragileEndR[c][right[c]+1];
      if(OutlierFrac[c])
	OutlierFrac[c][numsite[c]+1] = OutlierFrac[c][right[c]+1];
      if(FragSd[c])
	FragSd[c][numsite[c]+1] = FragSd[c][right[c]+1];
      if(ExpSd[c])
	ExpSd[c][numsite[c]+1] = ExpSd[c][right[c]+1];
      if(FragCov[c])
	FragCov[c][numsite[c]+1] = FragCov[c][right[c]+1];
      if(FragChiSq[c])
	FragChiSq[c][numsite[c]+1] = FragChiSq[c][right[c]+1];
    }
  }
}

void Cmap::inversion(double RangeLeft, double RangeRight)
{
  if(!(siteSD[0] && sitecov[0] && sitecnt[0])){
    printf("Cmap copy constructor is only implemented for Contig maps (from CMAP input): this map (id=%lld) is from %s\n",
	   id, vfixx_filename[fileid]);
    fflush(stdout);exit(1);
  }

  if(colors != 1){
    printf("Cmap::inversion():multiple color channels not implemented\n");
    exit(1);
  }
  int left = 0, right = numsite[0]+1;
  while(left < right && site[0][left] < RangeLeft - 0.0001)
    left++;
  while(left < right && site[0][right] > RangeRight + 0.0001)
    right--;
  if(site[0][left] < RangeLeft)
    RangeLeft = max(ZERO,site[0][left] - ((FLOAT)MININTERVAL));
  if(site[0][right] > RangeRight)
    RangeRight = min(site[0][numsite[0]+1],site[0][right] + ((FLOAT)MININTERVAL));
  
  if(left < right){
    /* first make a copy of the region to be inverted (site[0][left .. right]) */
    double *nsite = new double[right-left+1];
    double *nsiteSD = 0;
    float *nsitecov = 0,*nsitecnt = 0, *nChimQuality = 0, *nChimNorm = 0, *nSegDupL = 0, *nSegDupR = 0, *nFragileEndL = 0, *nFragileEndR = 0, *nOutlierFrac = 0, *nFragSd = 0, *nExpSd = 0, *nFragCov = 0, *nFragChiSq = 0;
    for(int i = left; i <= right; i++)
      nsite[i-left] = site[0][i];
    if(siteSD[0]){
      nsiteSD = new double[right-left+1];
      nsitecov = new float[right-left+1];
      nsitecnt = new float[right-left+1];

      for(int i = left; i <= right; i++){
	nsiteSD[i-left] = siteSD[0][i];
	nsitecov[i-left] = sitecov[0][i];
	nsitecnt[i-left] = sitecnt[0][i];
      }
      if(ChimQuality[0]){
	nChimQuality = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nChimQuality[i-left] = ChimQuality[0][i];
      }
      if(ChimNorm[0]){
	nChimNorm = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nChimNorm[i-left] = ChimNorm[0][i];
      }
      if(SegDupL[0]){
	nSegDupL = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nSegDupL[i-left] = SegDupL[0][i];
      }
      if(SegDupR[0]){
	nSegDupR = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nSegDupR[i-left] = SegDupR[0][i];
      }
      if(FragileEndL[0]){
	nFragileEndL = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nFragileEndL[i-left] = FragileEndL[0][i];
      }
      if(FragileEndR[0]){
	nFragileEndR = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nFragileEndR[i-left] = FragileEndR[0][i];
      }
      if(OutlierFrac[0]){
	nOutlierFrac = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nOutlierFrac[i-left] = OutlierFrac[0][i];
      }
      if(FragSd[0]){
	nFragSd = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nFragSd[i-left] = FragSd[0][i];
      }
      if(ExpSd[0]){
	nExpSd = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nExpSd[i-left] = ExpSd[0][i];
      }
      if(FragCov[0]){
	nFragCov = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nFragCov[i-left] = FragCov[0][i];
      }
      if(FragChiSq[0]){
	nFragChiSq = new float[right-left+1];
	for(int i = left; i <= right; i++)
	  nFragChiSq[i-left] = FragChiSq[0][i];
      }
    }

    /* copy the region in reverse order */
    for(int i = left; i <= right; i++){
      site[0][i] = RangeLeft + RangeRight - nsite[right-i];
      if(siteSD[0]){
	siteSD[0][i] = nsiteSD[right-i];
	sitecov[0][i] = nsitecov[right-i];
	sitecnt[0][i] = nsitecnt[right-i];
      }
      if(ChimQuality[0])
	ChimQuality[0][i] = nChimQuality[right-i];
      if(ChimNorm[0])
	ChimNorm[0][i] = nChimNorm[right-i];
      if(SegDupL[0])
	SegDupL[0][i] = nSegDupL[right-i];
      if(SegDupR[0])
	SegDupR[0][i] = nSegDupR[right-i];
      if(FragileEndL[0])
	FragileEndL[0][i] = nFragileEndL[right-i];
      if(FragileEndR[0])
	FragileEndR[0][i] = nFragileEndR[right-i];
      if(OutlierFrac[0])
	OutlierFrac[0][i] = nOutlierFrac[right-i];
      if(FragSd[0])
	FragSd[0][i] = nFragSd[right-i];
      if(ExpSd[0])
	ExpSd[0][i] = nExpSd[right-i];
      if(FragCov[0])
	FragCov[0][i] = nFragCov[right-i];
      if(FragChiSq[0])
	FragChiSq[0][i] = nFragChiSq[right-i];
    }
    delete [] nsite;
    if(siteSD[0]){
      delete [] nsiteSD;
      delete [] nsitecov;
      delete [] nsitecnt;
    }
    if(ChimQuality[0])
      delete [] nChimQuality;
    if(ChimNorm[0])
      delete [] nChimNorm;
    if(SegDupL[0])
      delete [] nSegDupL;
    if(SegDupR[0])
      delete [] nSegDupR;
    if(FragileEndL[0])
      delete [] nFragileEndL;
    if(FragileEndR[0])
      delete [] nFragileEndR;
    if(OutlierFrac[0])
      delete [] nOutlierFrac;
    if(FragSd[0])
      delete [] nFragSd;
    if(ExpSd[0])
      delete [] nExpSd;
    if(FragCov[0])
      delete [] nFragCov;
    if(FragChiSq[0])
      delete [] nFragChiSq;
  }
}

void maxmap_free()
{
  for(int i = 0; i < num_blocks; i++){
    if(VERB>=2 || LEAKDEBUG>=2){
      printf("maxmap_free():deleting Cmap_blocks[%d]=%p\n",i,Cmap_blocks[i]);
      fflush(stdout);
    }
    delete [] Cmap_blocks[i];
    Cmap_blocks[i] = NULL;
    if(site_blocks[i]){
      if(VERB>=2){
	printf("freeing site_blocks[%d]=%p\n",i,site_blocks[i]);
	fflush(stdout);
      }
      free(site_blocks[i]);
      site_blocks[i] = NULL;
    }
  }
  num_blocks = 0;
  
  if(MAX_BLOCKS > 0){
    free(Cmap_blocks); Cmap_blocks = NULL;
    if(VERB>=2 || LEAKDEBUG>=2){
      printf("MAX_BLOCKS = %d -> 0\n",MAX_BLOCKS);
      fflush(stdout);
    }
    free(site_blocks); site_blocks = NULL;
    free(id_blocks); id_blocks = NULL;
    free(mapcnt_blocks); mapcnt_blocks = NULL;

    MAX_BLOCKS = 0;
  }
  
  if(ChipIds){
    for(int i = 0; i < numChipIds; i++)
      free(ChipIds[i]);
    free(ChipIds);
    ChipIds = NULL;
    numChipIds = maxChipIds = 0;
  }
}

#define MINMAPS (64*1024)

static FLOAT aNaN = nan("NaN");

/* truncate any excess map allocations (excess block allocated memory is wasted). Used when the block allocation size has changed */
void maxmaptruncate(int nummaps, int &maxmaps, Cmap **&map)
{
  for(int i = nummaps; i < maxmaps; i++){
    Gmap[i]->allfree(); 
    Gmap[i] = NULL;
  }

  maxmaps = nummaps;
}

/** increase size of array Map[0..maxmaps-1] to at least Map[0..num-1] */
/** If minsites > 0, pre-allocate site[0..colors-1] etc (except stitch locations) */
void maxmapalloc(int num, int &maxmaps, Cmap **&Map, int minsites, int numthreads)
{
  int newmaps = num;
  if(newmaps > (int)(MASK(31)) || newmaps < 0){
    printf("maxmapalloc:num=%d is out of range for 32 bit signed number\n", newmaps);
    exit(1);
  }
  if(newmaps <= maxmaps)
    return;
  if(LEAKDEBUG <= 1){
    if(newmaps < MINMAPS)
      newmaps = MINMAPS;
    if(newmaps*2 < maxmaps*3)
      newmaps = (maxmaps*3)/2;
  }

  if(VERB>=2 || LEAKDEBUG>=2){
    printf("maxmapalloc(num=%d,maxmaps=%d(maxrefmaps=%d),minsites=%d,Map=%p(Gmap=%p,refmap=%p)):num_blocks=%d,newmaps=%d\n",
	   num,maxmaps,maxrefmaps,minsites,Map,Gmap,refmap,num_blocks,newmaps);
    fflush(stdout);
  }

  size_t mapsiz = (minsites+2 + (MapSNR ? minsites+1 : 0) + (MapIntensity ? minsites+1 : 0) + (MapPSFWidth ? minsites+1 : 0)) * colors * sizeof(double);
  if(MapTrueSite || Tracksites)
    mapsiz += (minsites+1) * 2 * colors * sizeof(long long);
  if(MapImageLoc){
    mapsiz += (minsites + 2 ) * colors * 2 * sizeof(double);
    mapsiz += (((minsites + 2 ) + 1) & ~1) * colors * sizeof(int);
  }
  if(MapStitched)
    mapsiz += (((minsites + 2 ) + 1) & ~1) * colors * sizeof(int);

  if(maxmaps <= 0){
    Cmap **newmap = new Cmap*[newmaps];
    Cmap *pmap = new Cmap[newmaps];
    maxblockalloc(num_blocks);
    mapcnt_blocks[num_blocks] = newmaps;

    // #pragma omp for num_threads(numthreads) schedule(static, 2048)
    for(int i = 0; i < newmaps; i++)
      newmap[i] = &pmap[i];

    site_blocks[num_blocks] = 0;
    if(minsites > 0){
      size_t blocksiz = mapsiz * newmaps;
      char *pmemblock = (char *)malloc(blocksiz);
      if(VERB>=2 || LEAKDEBUG>=2){
	printf("maxmapalloc: allocated map block %d with %d maps and %llu site bytes (site_blocks[%d]=%p)\n", 
	       num_blocks, newmaps, (unsigned long long)blocksiz, num_blocks, pmemblock);
	fflush(stdout);
      }
      if(!pmemblock){
	printf("maxmapalloc: malloc(%llu) failed\n", (unsigned long long)blocksiz);
	exit(1);
      }
      site_blocks[num_blocks] = pmemblock;

      #pragma omp parallel for num_threads(numthreads) schedule(static, 512) if(numthreads > 1)
      for(int i = 0; i < newmaps; i++){
	Cmap *p = newmap[i];
	char *pmem = &pmemblock[i * mapsiz];

	//#pragma unroll(0)
#pragma nounroll
	for(register int c = 0; c < colors; c++){
	  p->site[c] = (FLOAT *)pmem;
	  if(DEBUG>=2) 
	    for(int i = 0; i < minsites+2;i++)
	      p->site[c][i] = aNaN;
	  pmem += (minsites+2) * sizeof(FLOAT);

	  p->blockmem = minsites+2;
	  p->numsite[c] = minsites;

	  if(MapTrueSite || Tracksites){
	    p->truesiteL[c] = (long long *)pmem;
	    pmem += (minsites+1)*sizeof(long long);
	    p->truesiteR[c] = (long long *)pmem;
	    pmem += (minsites+1)*sizeof(long long);
	    if(DEBUG>=2) 
	      for(int i = 0; i <= minsites;i++)
		p->truesiteL[c][i] = p->truesiteR[c][i] = -1;
	  }
	  if(MapSNR){
	    p->SNR[c] = (double *)pmem;
	    if(DEBUG>=2) 
	      for(int i = 0; i < minsites+1;i++)
		p->SNR[c][i] = aNaN;
	    pmem += (minsites+1) * sizeof(double);
	  }
	  if(MapIntensity){
	    p->Intensity[c] = (double *)pmem;
	    if(DEBUG>=2) 
	      for(int i = 0; i < minsites+1;i++)
		p->Intensity[c][i] = aNaN;
	    pmem += (minsites+1) * sizeof(double);
	  }
	  if(MapStitched){
	    p->stitch[c] = (int *)pmem;
	    pmem += (((minsites+2) + 1) & ~1) * sizeof(int);
	  }
	  if(MapPSFWidth){
	    p->PSFWidth[c] = (double *)pmem;
	    if(DEBUG>=2)
	      for(int i = 1; i <= minsites; i++)
		p->PSFWidth[c][i] = aNaN;
	    pmem += (minsites+1) * sizeof(double);
	  }
	  if(MapImageLoc){
	    p->ImageFOV[c] = (int *)pmem;
	    pmem += (((minsites+2) + 1) & ~1) * sizeof(int);
	    p->ImageX[c] = (double *)pmem;
	    pmem += (minsites + 2) * sizeof(double);
	    p->ImageY[c] = (double *)pmem;
	    pmem += (minsites + 2) * sizeof(double);
	    if(DEBUG>=2)
	      for(int i = 0; i <= minsites+1; i++){
		p->ImageFOV[c][i] = -1;
		p->ImageX[c][i] = p->ImageY[c][i] = aNaN;
	      }
	  }
	}
	if(DEBUG>=2) assert(pmem <= &pmemblock[(i+1)*mapsiz]);
      }
    }
    id_blocks[num_blocks] = &Map;
    if(VERB>=2 || LEAKDEBUG>=2){
      printf("maxmapalloc:Cmap_blocks[%d] -> %p\n",num_blocks,pmap);
      fflush(stdout);
    }
    Cmap_blocks[num_blocks++] = pmap;

    Map = newmap;/* atomic update of global map[] pointer */

  } else {// maxmaps > 0

    Cmap **origmap = Map;
    Cmap **newmap = new Cmap*[newmaps];

    // #pragma omp for num_threads(numthreads) schedule(static, 2048)
    for(int i = 0; i < maxmaps; i++)
      newmap[i] = origmap[i];

    maxblockalloc(num_blocks);

    Cmap *pmap = new Cmap[newmaps - maxmaps];
    mapcnt_blocks[num_blocks] = newmaps - maxmaps;

    // #pragma omp for num_threads(numthreads) schedule(static, 2048)
    for(int i = maxmaps; i < newmaps; i++)
      newmap[i] = &pmap[i-maxmaps];

    site_blocks[num_blocks] = 0;
    if(minsites > 0){
      size_t blocksiz = mapsiz * (newmaps - maxmaps);
      char *pmemblock = (char *)malloc(blocksiz);
      // char *pmemblock = (char *)calloc(blocksiz,1);
      if(VERB>=2 || LEAKDEBUG>=2){
	printf("maxmapalloc: allocated map block %d with maps=%d->%d and %llu site bytes (site_blocks[%d]=%p)\n", 
	       num_blocks, maxmaps, newmaps, (unsigned long long)blocksiz, num_blocks, pmemblock);
	fflush(stdout);
      }
      if(!pmemblock){
	printf("maxmapalloc: malloc(%llu) failed\n", (unsigned long long)blocksiz);
	exit(1);
      }
      site_blocks[num_blocks] = pmemblock;

      #pragma omp parallel for num_threads(numthreads) schedule(static, 512) if(numthreads > 1)
      for(int i = maxmaps; i < newmaps; i++){
	Cmap *p = newmap[i];
	char *pmem = &pmemblock[(i-maxmaps) * mapsiz];

	for(int c = 0; c < colors; c++){
	  p->site[c] = (FLOAT *)pmem;
	  if(DEBUG>=2) 
	    for(int i = 0; i < minsites+2;i++)
	      p->site[c][i] = aNaN;
	  pmem += (minsites+2) * sizeof(FLOAT);

	  p->blockmem = minsites+2;
	  p->numsite[c] = minsites;
	  
	  if(MapTrueSite || Tracksites){
	    p->truesiteL[c] = (long long *)pmem;
	    pmem += (minsites+1)*sizeof(long long);
	    p->truesiteR[c] = (long long *)pmem;
	    pmem += (minsites+1)*sizeof(long long);
	    if(DEBUG>=2) 
	      for(int i = 0; i <= minsites;i++)
		p->truesiteL[c][i] = p->truesiteR[c][i] = -1;
	  }
	  if(MapSNR){
	    p->SNR[c] = (FLOAT *)pmem;
	    if(DEBUG>=2) 
	      for(int i = 0; i < minsites+1;i++)
		p->SNR[c][i] = aNaN;
	    pmem += (minsites+1) * sizeof(FLOAT);
	  }
	  if(MapIntensity){
	    p->Intensity[c] = (FLOAT *)pmem;
	    if(DEBUG>=2) 
	      for(int i = 0; i < minsites+1;i++)
		p->Intensity[c][i] = aNaN;
	    pmem += (minsites+1) * sizeof(FLOAT);
	  }
	  if(MapStitched){
	    p->stitch[c] = (int *)pmem;
	    pmem += (((minsites+2) + 1) & ~1) * sizeof(int);
	  }
	  if(MapPSFWidth){
	    p->PSFWidth[c] = (double *)pmem;
	    if(DEBUG>=2)
	      for(int i = 1; i <= minsites; i++)
		p->PSFWidth[c][i] = aNaN;
	    pmem += (minsites+1) * sizeof(double);
	  }
	  if(MapImageLoc){
	    p->ImageFOV[c] = (int *)pmem;
	    pmem += (((minsites+2) + 1) & ~1) * sizeof(int);
	    p->ImageX[c] = (double *)pmem;
	    pmem += (minsites + 2) * sizeof(double);
	    p->ImageY[c] = (double *)pmem;
	    pmem += (minsites + 2) * sizeof(double);
	    if(DEBUG>=2)
	      for(int i = 0; i <= minsites+1; i++){
		p->ImageFOV[c][i] = -1;
		p->ImageX[c][i] = p->ImageY[c][i] = aNaN;
	      }
	  }
	}
	if(DEBUG>=2) assert(pmem <= &pmemblock[(i+1-maxmaps)*mapsiz]);
      }
    }
    id_blocks[num_blocks] = &Map;
    if(VERB>=2 || LEAKDEBUG>=2){
      printf("maxmapalloc:Cmap_blocks[%d] -> %p\n",num_blocks,pmap);
      fflush(stdout);
    }
    Cmap_blocks[num_blocks++] = pmap;

    Map = newmap;/* atomic update of global Map[] (either Gmap[] or refmap[]) */
    delete [] origmap;
  }
  maxmaps = newmaps;
}

extern double wtime();

/** campact all maps in Map[0..nmaps-1] that are using block memory to use minimum needed memory */
void mapcompact(int nmaps, int &maxmaps, Cmap ** &Map)
{
  if(VERB>=2){
    printf("mapcompact: Estimating site memory: nmaps=%d,maxmaps=%d, time=%0.6f(wall time=%0.6f)\n", nmaps,maxmaps,mtime(), wtime());
    fflush(stdout);
  }

  int numthreads = 1;
#ifdef _OPENMP
  numthreads = omp_get_max_threads();
  numthreads = min(numthreads,MaxThreads);
  numthreads = max(1,min(numthreads, nmaps/2048));
#endif

  if(VERB>=2){
    int blockcnt = 0;
    #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
    {
      int myblockcnt = 0;

      #pragma omp for schedule(static,2048)
      for(int i = 0; i < nmaps; i++)
	myblockcnt += (Map[i]->blockmem ? 1 : 0);
      
      #pragma omp atomic
      blockcnt += myblockcnt;
    }

    printf("mapcompact: nmaps=%d of which %d had their block memory re-allocated:time=%0.6f(wall time=%0.6f)\n",nmaps,nmaps-blockcnt,mtime(),wtime());
    fflush(stdout);
  }

  /* count how many total sites we need (and previously used) */
  size_t sitecnt = 0, origcnt = 0, sitecntI = 0;

  #pragma omp parallel num_threads(numthreads) if(numthreads > 1)
  {
    size_t mysitecnt = 0, myorigcnt = 0, mysitecntI = 0;

    #pragma omp for schedule(static,2048)
    for(int i = 0; i < nmaps; i++){
      Cmap *p = Map[i];
      if(p->blockmem)
        myorigcnt += (p->blockmem - 2) * colors;
      else
        for(int c = 0; c < colors; c++)
	  myorigcnt += p->numsite[c];

      for(int c = 0; c < colors; c++){
        mysitecnt += p->numsite[c];
	mysitecntI += (p->numsite[c] + 1) & ~1;
      }
    }

    #pragma omp atomic
    sitecnt += mysitecnt;

    #pragma omp atomic
    origcnt += myorigcnt;

    #pragma omp atomic
    sitecntI += mysitecntI;
  }

  if(sitecnt <= 0)
    return;

  /* allocate a new block */
  maxblockalloc(num_blocks);

  size_t memsiz = (sitecnt + 2 * colors * nmaps) * sizeof(FLOAT);
  if(MapTrueSite || Tracksites)
    memsiz += (sitecnt + colors * nmaps) * 2 * sizeof(long long);
  if(MapSNR)
    memsiz += (sitecnt + colors * nmaps) * sizeof(double);
  if(MapIntensity)
    memsiz += (sitecnt + colors * nmaps) * sizeof(double);
  if(MapStitched)
    memsiz += (sitecntI + 2 * colors * nmaps) * sizeof(int);
  if(MapPSFWidth)
    memsiz += (sitecnt + colors * nmaps) * sizeof(double);
  if(MapImageLoc){
    memsiz += (sitecntI + 2 * colors * nmaps) * sizeof(int);
    memsiz += (sitecnt + 2 * colors * nmaps) * 2 * sizeof(double);
  }
  size_t offset = (nmaps+1) * sizeof(size_t);
  memsiz += offset;

  if(VERB){
    printf("reducing site memory from %lu sites to %lu sites (%lu bytes) for %d maps (nummaps=%d): time=%0.6f(wall time=%0.6f)\n",origcnt, sitecnt, memsiz, nmaps, nummaps, mtime(), wtime());
    fflush(stdout);
  }
  if(DEBUG) assert(nmaps >= 0);

  id_blocks[num_blocks] = NULL;
  char *mem = site_blocks[num_blocks] = (char *)malloc(memsiz);
  if(!mem){
    printf("mapcompact:malloc(%llu) failed\n", (unsigned long long)memsiz);
    exit(1);
  }
  if(VERB>=2){
    printf("mapcompact:site_blocks[%d] -> %p, memsiz= %lu, &mem[memsiz]= %p\n",num_blocks, mem, memsiz, &mem[memsiz]);
    fflush(stdout);
  }

  /* compute start of each maps memory as offset into mem */
  size_t *mapoffset = (size_t *)mem;
  mapoffset[0] = offset;

  if(VERB>=2){
    printf("nmaps=%d,nummaps=%d,maxmaps=%d,numthreads=%d\n",nmaps,nummaps,maxmaps,numthreads);
    fflush(stdout);
  }
  if(DEBUG) assert(nmaps >= 0);

  /* The following vectorized loop segfaults in MIC due to a Intel compiler bug, hence has been disabled when USE_MIC==1  */
  #pragma omp parallel for schedule(static, 2048) num_threads(numthreads) if(USE_MIC==0 && numthreads > 1)
  #pragma novector
  for(int i = 0; i < nmaps; i++){
    Cmap *p = Map[i];
    size_t siz = 0;
    
    for(int c = 0; c < colors; c++){
      int numsite = p->numsite[c];
      siz += (numsite + 2) * sizeof(FLOAT);
      if(MapTrueSite || Tracksites)
	siz += (numsite + 1) * 2 * sizeof(long long);
      if(MapSNR)
	siz += (numsite + 1) * sizeof(double);
      if(MapIntensity)
	siz += (numsite + 1) * sizeof(double);
      if(MapStitched)
	siz += (((numsite + 2) + 1) & ~1) * sizeof(int);
      if(MapPSFWidth)
	siz += (numsite + 1) * sizeof(double);
      if(MapImageLoc){
	siz += (numsite + 2) * 2 * sizeof(double);
	siz += (((numsite + 2) + 1) & ~1) * sizeof(int);
      }
    }
    mapoffset[i+1] = siz;
  }

  /* prefix sum computation : not easy to multithread */
  for(int i = 1; i <= nmaps; i++){
    offset += mapoffset[i];
    mapoffset[i] = offset;
  }

  if(VERB>=2){
    printf("offset=%lu,memsiz=%lu\n",offset,memsiz);
    fflush(stdout);
  }

  if(DEBUG) assert(offset == memsiz);

  if(VERB>=2){
    printf("Allocated compact site memory: time=%0.6f(wall time=%0.6f)\n", mtime(), wtime());
    fflush(stdout);
  }

  #pragma omp parallel for num_threads(numthreads) schedule(dynamic,1024) if(numthreads > 1)
  for(int i = 0; i < nmaps; i++){
    Cmap *p = Map[i];
    char *pmem = &mem[mapoffset[i]];
    for(int c = 0; c < colors; c++){
      int numsite = p->numsite[c];

      FLOAT *origmem = p->site[c];
      p->site[c] = (FLOAT *)pmem;
      size_t siz = (numsite + 2) * sizeof(FLOAT);
      pmem += siz;
      memcpy(p->site[c],origmem,siz);
      if(!p->blockmem)
	delete [] origmem;

      if(MapTrueSite || Tracksites){
	long long *origLLmem = p->truesiteL[c];
	p->truesiteL[c] = (long long *)pmem;
	siz = (numsite + 1) * sizeof(long long);
	pmem += siz;
	memcpy(p->truesiteL[c],origLLmem,siz);
	if(!p->blockmem)
	  delete [] origLLmem;
	
	origLLmem = p->truesiteR[c];
	p->truesiteR[c] = (long long *)pmem;
	pmem += siz;
	memcpy(p->truesiteR[c],origLLmem,siz);
	if(!p->blockmem)
	  delete [] origLLmem;
      }

      origmem = p->SNR[c];
      if(MapSNR){
	p->SNR[c] = (double *)pmem;
	siz = (numsite + 1) * sizeof(double);
	pmem += siz;
	memcpy(p->SNR[c],origmem,siz);
      }
      if(!p->blockmem && origmem)
	delete [] origmem;

      origmem = p->Intensity[c];
      if(MapIntensity){
	p->Intensity[c] = (double *)pmem;
	siz = (numsite + 1) * sizeof(double);
	pmem += siz;
	memcpy(p->Intensity[c], origmem, siz);
      }
      if(!p->blockmem && origmem)
	delete [] origmem;

      int *origstitch = p->stitch[c];
      if(MapStitched){
	p->stitch[c] = (int *)pmem;
	siz = (numsite + 2) * sizeof(int);
	pmem += (siz + 0x7) & ~0x7;
	memcpy(p->stitch[c], origstitch, siz);
      }
      if(!p->blockmem && origstitch)
	delete [] origstitch;

      origmem = p->PSFWidth[c];
      if(MapPSFWidth){
	p->PSFWidth[c] = (double *)pmem;
	siz = (numsite + 1) * sizeof(double);
	pmem += siz;
	memcpy(p->PSFWidth[c], origmem, siz);
      }
      if(!p->blockmem)
	delete [] origmem;

      int *origFOV = p->ImageFOV[c];
      double *origX = p->ImageX[c];
      double *origY = p->ImageY[c];
      if(MapImageLoc){
        p->ImageFOV[c] = (int *)pmem;
	siz = (numsite + 2) * sizeof(int);
	pmem += (siz + 0x7) & ~0x7;
	memcpy(p->ImageFOV[c], origFOV, siz);

	p->ImageX[c] = (double *)pmem;
	siz = (numsite + 2) * sizeof(double);
	pmem += siz;
	memcpy(p->ImageX[c], origX, siz);

	p->ImageY[c] = (double *)pmem;
	siz = (numsite + 2) * sizeof(double);
	pmem += siz;
	memcpy(p->ImageY[c], origY, siz);
      }
      if(!p->blockmem){
        delete [] origFOV;
	delete [] origX;
	delete [] origY;
      }
    }
    if(DEBUG) assert(pmem <= &mem[memsiz]);

    p->blockmem = p->numsite[0] + 2;
    for(int c = 1; c < colors; c++)
      p->blockmem = max(p->blockmem,p->numsite[c]+2);
  }

  if(VERB>=2){
    printf("Copied site memory: time=%0.6f(wall time=%0.6f)\n", mtime(), wtime());
    fflush(stdout);
  }

  int fcnt = 0;
  for(int i = 0; i < num_blocks; i++)
    if(id_blocks[i] == &Map){
      if(site_blocks[i]){
       if(VERB>=2){
	  printf("freeing site_blocks[%d]=%p\n",i,site_blocks[i]);
	  fflush(stdout);
	}
	fcnt++;
	free(site_blocks[i]);
	site_blocks[i] = 0;
      }
    }

  size_t mcnt = 0;
  for(int i = 0; i < num_blocks; i++)
    if(id_blocks[i] == &Map)    
      mcnt += mapcnt_blocks[i];

  if(VERB){
    if(VERB>=2)
      printf("Freed %d/%d site blocks:",fcnt,num_blocks);
    printf("reducing Cmap memory from %lu to %lu bytes (%lu to %d maps): time=%0.6f(wall time=%0.6f)\n", mcnt*sizeof(Cmap), nmaps*sizeof(Cmap), mcnt, nmaps, mtime(),wtime());
    fflush(stdout);
  }

  if(VERB>=2 || LEAKDEBUG>=2){
    printf("mapcompact:Cmap_blocks[%d] -> 0\n",num_blocks);
    fflush(stdout);
  }

  Cmap_blocks[num_blocks] = 0;

  int newmaps = min(maxmaps,nmaps+16);

  Cmap *pmap = new Cmap[newmaps];
  if(VERB>=2 || LEAKDEBUG>=2){
    printf("mapcompact:Cmap_blocks[%d]=%p -> %p\n",num_blocks,Cmap_blocks[num_blocks],pmap);
    fflush(stdout);
  }
  Cmap_blocks[num_blocks] = pmap;

  if(VERB>=2){
    printf("Starting copy of Cmap memory: time=%0.6f(wall time=%0.6f)\n", mtime(), wtime());
    fflush(stdout);
  }

  #pragma omp parallel for num_threads(numthreads) schedule(dynamic, 1024) if(numthreads > 1)
  for(int i = 0; i < nmaps /* WAS31 newmaps */; i++){
    Cmap *origmap = Map[i];
    Map[i] = &pmap[i];
    memcpy(Map[i],origmap,sizeof(Cmap));
    memset(origmap,0,sizeof(Cmap));/* to keep Cmap site arrays from being deallocated when original Cmap's are deallocated */
  }

  for(int i = nmaps; i < newmaps; i++)// NEW33
    Map[i] = &pmap[i];

  if(VERB>=2){
    printf("Copied Cmap memory: time=%0.6f(wall time=%0.6f)\n", mtime(),wtime());
    fflush(stdout);
  }

  /* NOTE : the Map[] pointer array is not compacted since that saves little memory */

  maxmaps = newmaps;

  fcnt = 0;
  for(int i = 0; i < num_blocks; i++)
    if(id_blocks[i] == &Map){
      if(Cmap_blocks[i]){
	fcnt++;
	if(VERB>=2 || LEAKDEBUG>=2){
          printf("mapcompact:deleting Cmap_blocks[%d]=%p -> 0\n",i,Cmap_blocks[i]);
	  fflush(stdout);
        }

	delete [] Cmap_blocks[i];
	Cmap_blocks[i] = 0;
      }
    }

  if(VERB>=2){
    printf("Free %d/%d Cmap blocks:time=%0.6f(wall time=%0.6f)\n", fcnt, num_blocks, mtime(), wtime());
    fflush(stdout);
  }

  num_blocks++;
}

#define DELTA_BLOCK 16
static PFLOAT *deltaXmem[DELTA_BLOCK];
static PFLOAT **deltaXpmem[DELTA_BLOCK];
static int numblock = 0;

void deltaXinitBlock(Cmap **rmap, int mapstart, int mapend)
{
  if(VERB>=2){
    printf("deltaXinitBlocK:mapstart=%d,mapend=%d,DELTA=%d,numblock=%d (starting):wall time=%0.6f\n",mapstart,mapend,DELTA,numblock,wtime());
    fflush(stdout);
  }

  long long *fmem_start = new long long[(mapend+2)*2];
  long long *pmem_start = &fmem_start[mapend+2];

  long long fmem = 0, pmem = 0 ;
  for(int i = mapstart; i <= mapend; i++){
    Cmap *pmap = rmap[i];
    fmem_start[i] = fmem;
    pmem_start[i] = pmem;
    for(int c = 0; c < colors; c++){
      int M = pmap->numsite[c];
      if(((USE_SSE && USE_AVX) || USE_MIC) && USE_PFLOAT){
	int M_stride = COMPUTE_STRIDE4(M-1);
	if(DEBUG) assert(M_stride >= M-1);
	fmem += M_stride*DELTA*3LL;
	pmem += 2LL*(DELTA+1) + M+1;
      } else {
	fmem += max(0,M-1)*2LL*DELTA;
	pmem += 2LL*(M+1);
      }
    }    
  }

  if(DEBUG){
    fmem_start[mapend+1] = fmem;
    pmem_start[mapend+1] = pmem;
    if(DEBUG>=2){
      for(int i = mapstart + 1; i <= mapend; i++){
	assert(fmem_start[i] >= fmem_start[i-1]);
	assert(pmem_start[i] >= pmem_start[i-1]);
      }
    }
  }

  if(numblock >= DELTA_BLOCK){
    printf("deltaXinitBlock:numblock=%d : increase DELTA_BLOCK in Cmap.cpp\n",numblock);
    fflush(stdout);
  }
  PFLOAT *fmemblock = deltaXmem[numblock] = new PFLOAT[fmem];
  PFLOAT **pmemblock = deltaXpmem[numblock] = new PFLOAT*[pmem];

  int nthreads = 1;
  #ifdef _OPENMP
  nthreads = omp_get_max_threads();
  nthreads = min(nthreads,MaxThreads);
  #endif

  #pragma omp parallel for num_threads(nthreads) schedule(dynamic,64) if(nthreads > 1)
  for(int i = mapstart; i <= mapend; i++){
    Cmap *pmap = rmap[i];
    long long fcnt = fmem_start[i];
    long long pcnt = pmem_start[i];

    for(int c = 0; c < colors; c++){
      int M = pmap->numsite[c];
      FLOAT *X = pmap->site[c];
      if(((USE_SSE && USE_AVX) || USE_MIC) && USE_PFLOAT){
	int M_stride = COMPUTE_STRIDE4(M-1);
	PFLOAT **pdeltaXc = pmap->deltaX[c] = &pmemblock[pcnt]; pcnt += DELTA+1;
	PFLOAT **pdeltaRc = pmap->deltaR[c] = &pmemblock[pcnt]; pcnt += DELTA+1;
	for(int n = 1; n <= DELTA; n++){
	  PFLOAT *pdeltaXcn = pdeltaXc[n]  = &fmemblock[fcnt-2]; fcnt += M_stride;
	  for(int J = n+1; J <= M; J++)
	    pdeltaXcn[J] = X[J] - X[J-n];
	}
	for(int n = 1; n <= DELTA; n++){
	  PFLOAT *pdeltaRcn = pdeltaRc[n] = &fmemblock[fcnt-2]; fcnt += M_stride;
	  for(int J = n+1; J <= M; J++)
	    pdeltaRcn[J] = X[M+1-J+n] - X[M+1-J];
	}
	PFLOAT **pdeltaYc = pmap->deltaY[c] = &pmemblock[pcnt]; pcnt += M+1;
	for(int J = 2; J <= M; J++){
	  PFLOAT *pdeltaYcJ = pdeltaYc[J] = &fmemblock[fcnt]; fcnt += DELTA;
	  FLOAT Xj = X[J];
	  FLOAT *pX = &X[J-DELTA];
	  int s = 0;
	  if(J < DELTA+1)
	    s = DELTA+1-J;
	  for(int n = s; n < DELTA; n++)
	    pdeltaYcJ[n] = Xj - pX[n];
	}
      } else {
	PFLOAT **pdeltaXc = pmap->deltaY[c] = pmap->deltaX[c] = &pmemblock[pcnt]; pcnt += M+1;
	PFLOAT **pdeltaRc = pmap->deltaR[c] = &pmemblock[pcnt]; pcnt += M+1;
	for(int J = 2; J <= M; J++){
	  PFLOAT *pdeltaXcJ = pdeltaXc[J] = &fmemblock[fcnt]; fcnt += DELTA;
	  FLOAT Xj = X[J];
	  FLOAT *pX = &X[J-DELTA];
	  int s = 0;
	  if(J < DELTA+1)
	    s = DELTA+1-J;
	  for(int n = s; n < DELTA; n++)
	    pdeltaXcJ[n] = Xj - pX[n];
	  PFLOAT *pdeltaRcJ = pdeltaRc[J] = &fmemblock[fcnt]; fcnt += DELTA;
	  Xj = X[M+1-J];
	  pX = &X[M+1+DELTA-J];
	  for(int n = s; n < DELTA; n++)
	    pdeltaRcJ[n] = pX[-n] - Xj;
	}
      }
    }/* c = 0 .. colors-1 */
    if(DEBUG){
      if(!(fcnt <= fmem_start[i+1])){
	printf("i=%d(%d..%d):fmem_start[i+1]=%lld,fmem_start[i]=%lld,fcnt=%lld,DELTA=%d,colors=%d,M=%d\n",i,mapstart,mapend,fmem_start[i+1],fmem_start[i],fcnt,DELTA,colors,pmap->numsite[0]);
	fflush(stdout);
	assert(fcnt <= fmem_start[i+1]);
      }
      assert(pcnt <= pmem_start[i+1]);
    }
  }/* i = mapstart .. mapend */

  numblock++;

  delete [] fmem_start;

  if(VERB>=2){
    printf("deltaXinitBlocK:mapstart=%d,mapend=%d,DELTA=%d,numblock=%d (completed):wall time=%0.6f\n",mapstart,mapend,DELTA,numblock,wtime());
    fflush(stdout);
  }
}

void deltaXfreeBlock(Cmap **rmap, int mapstart, int mapend)
{
  if(VERB>=2 || LEAKDEBUG>=2){
    printf("deltaXfreeBlocK:mapstart=%d,mapend=%d,DELTA=%d,numblock=%d (starting):wall time=%0.6f\n",mapstart,mapend,DELTA,numblock,wtime());
    fflush(stdout);
  }

#if DELTAX_BLOCK==0 // OLD
  for(int i = mapstart; i <= mapend; i++)
    rmap[i]->deltaXfree();
#else // NEW

  /* free up deltaX for split maps */ 
  for(int i = mapstart; i <= mapend; i++){
    if(DEBUG) assert(rmap[i]);
    if(rmap[i]->origmap)
      rmap[i]->deltaXfree();/* HERE causes segfault with LEAKDEBUG (check with Valgrind) */
  }

  while(--numblock >= 0){
    delete [] deltaXmem[numblock];
    delete [] deltaXpmem[numblock];
  }
  numblock = 0;

  for(int i = mapstart; i <= mapend; i++){
    if(DEBUG) assert(rmap[i]);
    Cmap *pmap = rmap[i];
    for(int c = 0;c < colors;c++)
      pmap->deltaX[c] = pmap->deltaR[c] = 0;
  }
#endif

  if(VERB>=2 || LEAKDEBUG>=2){
    printf("deltaXfreeBlocK:mapstart=%d,mapend=%d,DELTA=%d,numblock=%d (completed):wall time=%0.6f\n",mapstart,mapend,DELTA,numblock,wtime());
    fflush(stdout);
  }
}

static int rawsiteblock = -1;/* If >= 0, memory block last used to allocate rawsites with start==0 was rawsiteblock */
static int rawsiteblockend = -1;/* value of end arg when rawsitealloc() was last called with start==0 */
/* NOTE : if -nosplit 1 or 0 is combined with -minSNRestimate, a memory leak will occur, since only the first block with start==0 will be 
   tracked and deallocated when minSNR is reapplied globally */

/* uses block allocation to reduce new/delete calls */
void rawsitealloc(Cmap **Map, int start, int end)
{
  if(VERB>=2){
    printf("rawsitealloc:start=%d,end=%d,rawsiteblock=%d,Map=%p, map=%p:time=%0.6f(wall time=%0.6f)\n",start,end,rawsiteblock,Map,Gmap,mtime(),wtime());
    fflush(stdout);
  }

  if(start >= end)
    return;

  int nthreads = max(1,min(MaxThreads, nummaps/2048));
  int origcolors = colors;
  if(usecolor)
    colors = 2;

  /* compute total memory needed for rawsite[] */
  int sitecnt = 0;
  #pragma omp parallel num_threads(nthreads) 
  {
    int mysitecnt = 0;
    
    #pragma omp for schedule(static,2048)
    for(int mapid = start; mapid < end; mapid++){
      Cmap *p = Map[mapid];
      for(int c = 0; c < colors; c++)
	mysitecnt += p->numsite[c] + 2;
    }

    #pragma omp atomic
    sitecnt += mysitecnt;
  }

  /* allocate new block */
  if(num_blocks >= MAX_BLOCKS){
    printf("rawsitealloc:num_blocks=%d:increase MAX_BLOCKS in Cmap.cpp\n",num_blocks);
    exit(1);
  }

  Cmap_blocks[num_blocks] = NULL;
  if(VERB>=2 || LEAKDEBUG>=2){
    printf("rawsitealloc:Cmap_blocks[%d] -> 0\n", num_blocks);
    fflush(stdout);
  }
  mapcnt_blocks[num_blocks] = 0;
  id_blocks[num_blocks] = &Map;
  size_t memsiz = sitecnt * sizeof(FLOAT);
  memsiz += (end-start+2)*sizeof(size_t);
  char *mem = site_blocks[num_blocks] = (char *)malloc(memsiz);
  if(!mem){
    printf("rawsitealloc:malloc(%llu) failed\n", (unsigned long long)memsiz);
    exit(1);
  }
  if(VERB>=2){
    printf("rawsitealloc:site_blocks[%d] -> %p (memsiz= %lu, sitecnt= %d, start= %d, end=%d)\n",num_blocks,mem, memsiz, sitecnt, start, end);
    fflush(stdout);
  }

  size_t *mapoffset = (size_t *)mem;
  mapoffset -= start;
  size_t offset = (end-start+2) * sizeof(size_t);

  /* compute mapoffset[start..end] as byte index into mem[] */

  #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1024)
  for(int i = start; i < end; i++){
    Cmap *p = Map[i];
    size_t siz = 0;
    for(int c = 0; c < colors; c++)
      siz += (p->numsite[c] + 2) * sizeof(FLOAT);
    mapoffset[i+1] = siz;
  }

  for(int i = start; i < end; i++){
    mapoffset[i] = offset;
    offset += mapoffset[i+1];
  }
  if(DEBUG && !(offset == memsiz)){
    printf("offset = %lu, memsiz= %lu\n",offset,memsiz);
    int totsite = 0;
    for(int i = start; i < end; i++){
      Cmap *p = Map[i];
      for(int c = 0; c < colors; c++)
	totsite += p->numsite[c] + 2;
    }
    printf("\t sitecnt= %d, start=%d, end=%d : totsite= %d\n",sitecnt, start, end, totsite);
    fflush(stdout);
    assert(offset == memsiz);
  }

  if(VERB>=2){
    printf("rawsitealloc:memory allocated and block offsets computed:time=%0.6f(wall time=%0.6f)\n",mtime(),wtime());
    fflush(stdout);
  }
  
  /* use Map[]->rawsite[][] to save current site locations */

  #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1024)
  for(int mapid = start; mapid < end; mapid++){
    Cmap *pmap = Map[mapid];
    FLOAT *pmem = (FLOAT *)&mem[mapoffset[mapid]];
    for(int c = 0; c < colors; c++){
      if(pmap->rawsite[c] && rawsiteblock < 0 && start==0){
#pragma omp critical
	{
	  printf("WARNING in rawsitealloc for mapid=%d(id=%lld),pmap->rblockmem=%d,color=%d:rawsite[%d]=%p already allocated!\n", mapid,pmap->id,pmap->rblockmem,c+1,c,pmap->rawsite[c]);
	  fflush(stdout);
	}
      }
      if(!pmap->rblockmem)
	delete [] pmap->rawsite[c];
      pmap->rawsite[c] = pmem;
      int inc = pmap->numsite[c]+2;
      pmem += inc;
      memcpy(pmap->rawsite[c],pmap->site[c],inc*sizeof(FLOAT));
      if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
        #pragma omp critical
	{
	  printf("rawsitealloc:pmap->id=%lld,c=%d,inc=%d,numsite[c]=%d:\n",pmap->id,c,inc,pmap->numsite[c]);
	  for(int I = 0; I <= pmap->numsite[c] + 1; I++)
	    if(MapSNR)
	      printf("  rawsite[c][%d]= %0.6f, SNR= %0.3f\n",I,pmap->rawsite[c][I], (1 <= I && I <= pmap->numsite[c]) ? pmap->SNR[c][I] : -1.0);
	    else
	      printf("  rawsite[c][%d]= %0.6f\n",I,pmap->rawsite[c][I]);
	  fflush(stdout);
        }
      }
    }
    pmap->rblockmem = pmap->numsite[0]+2;
    for(int c = 1; c < colors; c++)
      pmap->rblockmem = max(pmap->blockmem,pmap->numsite[c]+2);
    if(DEBUG) assert(pmap->rblockmem > 0);
  }

  if(start == 0){/* try to free up main block previously allocated with start==0 */
    if(rawsiteblock >= 0){
      if(VERB>=2){
	printf("rawsitealloc:free memory block site_blocks[%d]=%p\n",rawsiteblock,site_blocks[rawsiteblock]);
	fflush(stdout);
      }
      free(site_blocks[rawsiteblock]);
      site_blocks[rawsiteblock] = NULL;

      if(NoSplit != 2 && rawsiteblockend > 0 && end > rawsiteblockend){
	printf("WARNING in rawsitealloc(): using -minSNRestimate with -nosplit 1 or 0 will cause some memory leak:start=%d,end=%d,rawsiteblockend=%d\n",
	       start,end,rawsiteblockend);
	fflush(stdout);
      }
    }

    rawsiteblock = num_blocks++;
    rawsiteblockend = end;
  }

  if(usecolor)
    colors = origcolors;

  if(VERB>=2){
    printf("rawsitealloc completed: time=%0.6f (wall time=%0.6f)\n",mtime(),wtime());
    fflush(stdout);
  }
}

void BiasCorrectSlope(double **Offset, double **Slope)
{
  for(int c = 0; c < colors; c++){
    Slope[c][-1] = Slope[c][ResBins[c]] = 0.0;
    Offset[c][-1] = resbias[c][0];
    Offset[c][ResBins[c]] = resbias[c][ResBins[c]];

    for(int Bin = 0; Bin < ResBins[c]; Bin++){
      if(DEBUG) assert(resbiasX[c][Bin+1] > resbiasX[c][Bin]);
#if SIMPLE_BIAS>=2
      Slope[c][Bin] = (resbias[c][Bin+1] - resbias[c][Bin])/(resbiasX[c][Bin+1] - resbiasX[c][Bin]);
      Offset[c][Bin] = resbias[c][Bin] - Slope[c][Bin] * resbiasX[c][Bin];
#else
      Slope[c][Bin] = 0.0;
      Offset[c][Bin] = resbias[c][Bin];
#endif
      if(VERB>=2 || (DEBUG/* HERE >=2 */ && !(isfinite(Slope[c][Bin]) && isfinite(Offset[c][Bin])))){
	printf("c=%d,Bin=%d:range=%0.4f..%0.4f, bias=%0.4f..%0.4f, Slope[c][Bin]=%0.8f,Offset[c][Bin]=%0.8f\n",c,Bin,
	       resbiasX[c][Bin],resbiasX[c][Bin+1],resbias[c][Bin],resbias[c][Bin+1],Slope[c][Bin],Offset[c][Bin]);
	fflush(stdout);
	if(DEBUG) assert(isfinite(Slope[c][Bin]));
	if(DEBUG) assert(isfinite(Offset[c][Bin]));
      }
    }
  }
} 

void BiasCorrect(Cmap **Map, int orignummaps, int nummaps, int rawalloc)
{
  if(DEBUG && rawalloc && !(orignummaps == 0 && ResBins[0] > 0)){
    printf("BiasCorrect:orignummaps=%d,nummaps=%d,ResBins[0]=%d,parametersfile=%s\n",orignummaps,nummaps, ResBins[0],parametersfile);
    fflush(stdout);
    assert(orignummaps == 0 && ResBins[0] > 0);
  }
  if(VERB>=2){
    printf("BiasCorrect:orignummaps=%d,nummaps=%d,rawalloc=%d\n",orignummaps,nummaps,rawalloc);
    for(int c = 0; c < colors; c++){
      printf("BiasCorrect:c=%d,maxresbias=%0.10f,ResBins[c]=%d:\n",c,maxresbias,ResBins[c]);
      for(int Bin = 0; Bin <= ResBins[c]; Bin++)
	printf("Bin=%d:resbiasX=%0.10f,resbias=%0.10f\n",Bin,resbiasX[c][Bin],resbias[c][Bin]);
    }
    fflush(stdout);
  }

  /* apply bias correction to all small intervals */
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
  nthreads = min(nthreads,MaxThreads);
  nthreads = max(1,min(nthreads, nummaps/2048));
#endif  

  if(rawalloc) /* allocate rawsite[] and copy site[] values to it */
    rawsitealloc(Map,orignummaps,nummaps);

  /* compute Offset[c][-1..ResBins], Slope[c][-1..ResBins] */
  double OffsetMem[MAXCOLOR][RESBINS+2], SlopeMem[MAXCOLOR][RESBINS+2];
  double *Offset[MAXCOLOR], *Slope[MAXCOLOR];
  for(int c = 0; c < colors; c++){
    Offset[c] = &OffsetMem[c][1];
    Slope[c] = &SlopeMem[c][1];
  }
  BiasCorrectSlope(Offset,Slope);

  /* create table to quickly determine the Bin of any interval size up to maxresbias : should match code in refalign.cpp */
  int Imaxresbias = (int)floor(maxresbias*1000.0+0.5);
  int *SizeToBin[MAXCOLOR] = {0,0};
  for(int c = 0; c < colors; c++)
    SizeToBin[c] = new int[Imaxresbias + 1];

  for(int c = 0; c < colors; c++){
    if(DEBUG && !(resbiasX[c][ResBins[c]] <= maxresbias)){
      printf("BiasCorrect:c=%d,colors=%d:ResBins[c]=%d,resbiasX[c][ResBins[c]]=%0.6f,maxresbias=%0.6f\n",
	     c,colors,ResBins[c],resbiasX[c][ResBins[c]],maxresbias);
      fflush(stdout);
      assert(resbiasX[c][ResBins[c]] <= maxresbias);
    }

    int Ibot = (int)floor(resbiasX[c][0]*1000.0 + 0.5), Itop = 0;
    for(int i = 0; i < Ibot; i++)
      SizeToBin[c][i] = 0;
    for(int Bin = 0; Bin < ResBins[c]; Bin++, Ibot = Itop){
      Itop = (int)floor(resbiasX[c][Bin+1]*1000.0 + 0.5);
      if(DEBUG/* HERE >=2 */ && Bin) assert(Ibot >= 0);
      if(DEBUG/* HERE >=2 */) assert(Ibot <= Itop && Itop <= Imaxresbias);
      for(int i = Ibot; i < Itop; i++)
	SizeToBin[c][i] = Bin;
      //	  SizeToBin[c][Itop] = Bin + 1;
    }
    if(DEBUG) assert(Itop <= Imaxresbias);
    for(int i = Ibot /* Itop+1*/ ; i <= Imaxresbias; i++)
      SizeToBin[c][i] = ResBins[c];
  }

  BiasCorrect2(Map, orignummaps, nummaps, nthreads, SizeToBin, Imaxresbias, Offset, Slope);

  for(int c = 0; c < colors; c++){
    delete [] SizeToBin[c];
   
    if(DEBUG>=2 && ResBins[c] > 0 && !(maxresbias >= resbiasX[c][ResBins[c]] - 1e-6)){
      printf("BiasCorrect:c=%d:ResBins[c]=%d,maxresbias=%0.6f,resbiasX[c][ResBins[c]]=%0.6f\n",
	     c,ResBins[c],maxresbias,resbiasX[c][ResBins[c]]);
      fflush(stdout);
      assert(maxresbias >= resbiasX[c][ResBins[c]] - 1e-6);
    }
  }
}

void BiasCorrect2(Cmap **Map, int orignummaps, int nummaps, int nthreads, int **SizeToBin, int Imaxresbias, double **Offset, double **Slope)
{

  int origcolors = colors;

  if(colors==1 && usecolor){/* create reasonable defaults for 2nd color */
    ResBins[1] = ResBins[0];
    for(int bin = 0; bin <= ResBins[0]; bin++){
      resbias[1][bin] = resbias[0][bin];
      resbiasX[1][bin] = resbiasX[0][bin];
    }

    Offset[1] = Offset[0];
    Slope[1] = Slope[0];
    SizeToBin[1] = SizeToBin[0];

    colors = 2;
  }

  /* Apply full bias removal to all sites in Gmap[orignummaps..nummaps-1] */

  double minInterval = 1.0e10;
  int minMapid = -1;

  #pragma omp parallel num_threads(nthreads)
  {
   double myMinInterval = 1.0e10;
   int myMapid = -1;

   #pragma omp for schedule(dynamic,256)
   for(int mapid = orignummaps; mapid < nummaps; mapid++){
    Cmap *pmap = Map[mapid];
    if(DEBUG>=2 && colors>=2) assert(fabs(pmap->site[0][pmap->numsite[0]+1] - pmap->site[1][pmap->numsite[1]+1]) < 1e-3);
    
    for(int c = 0; c < colors; c++){
      int M = pmap->numsite[c];
      if(M <= 1)
	continue;

      FLOAT *X = pmap->site[c];
      FLOAT *RawX = pmap->rawsite[c];

      if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
        #pragma omp critical
        {
          printf("mapid=%d(id=%lld),M=%d\n",pmap->mapid,pmap->id, M);
          fflush(stdout);
        }
      }

      for(int I = 1; I <= M; I++){

	double site = RawX[I], interval;
	double shift = 0.0;

	if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
	  #pragma omp critical
          {
            printf("  I=%d,rawX[I]=%0.10f,old X[I]=%0.10f\n",I,RawX[I],X[I]);
	    fflush(stdout);
          }
        }
	if(I > 1 && RawX[I] - RawX[I-1] < myMinInterval){
	  myMinInterval = RawX[I]-RawX[I-1];
	  myMapid = pmap->mapid;
	}

	/* first handle left neighbors */
	for(int J = I; --J >= 1;){
	  if((interval = RawX[I] - RawX[J]) > maxresbias)
	    break;
	  int isize = (int)floor(interval * 1000.0 + 0.5);
	  if(DEBUG) assert(isize <= Imaxresbias);
	  int Bin = SizeToBin[c][isize];
	  if(DEBUG) assert(0 <= Bin && Bin <= ResBins[c]);

	  /* correct Bin for boundary cases */
	  while(Bin > 0 && interval < resbiasX[c][Bin])
	    Bin--;
	  if(DEBUG) assert(0 <= Bin && Bin <= ResBins[c]);
	  if(DEBUG && Bin > 0 && !(resbiasX[c][Bin]  <= interval)){
	    printf("c=%d, interval = %0.8f, isize = %d, SizeToBin[c][isize]=%d, Bin= %d, resbiasX[c][%d]=%0.8f, resbiasX[c][%d]=%0.8f,resbiasX[c][%d]=%0.8f\n",
		   c, interval, isize, SizeToBin[c][isize], Bin, max(0,Bin-1), resbiasX[c][max(0,Bin-1)], Bin, resbiasX[c][Bin], min(Bin+1,ResBins[c]),resbiasX[c][min(Bin+1,ResBins[c])]);
	    fflush(stdout);
	    assert(resbiasX[c][Bin] <= interval);
	  }
	  if(DEBUG && Bin < ResBins[c]) assert(interval < resbiasX[c][Bin+1]);

	  double bias = Offset[c][Bin] + Slope[c][Bin] * interval;
	  shift -= bias;/* NOTE : bias < 0 */
	  if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
	    #pragma omp critical
	    {
              printf("    c=%d,I=%d,rawX[I]=%0.10f:rawX[J=%d]=%0.10f,bias=%0.10f,X[I]-> %0.10f(shift=%0.10f)\n", c,I, RawX[I], J, RawX[J], bias, site + shift*0.5, shift*0.5);
	      fflush(stdout);
	    }
          }
	}

	/* next handle right neigbors */
	for(int J = I; ++J <= M; ){
	  if((interval = RawX[J] - RawX[I]) > maxresbias)
	    break;
	  int isize = (int)floor(interval * 1000.0 + 0.5);
	  if(DEBUG) assert(isize <= Imaxresbias);
	  int Bin = SizeToBin[c][isize];
	  if(DEBUG) assert(0 <= Bin && Bin <= ResBins[c]);

	  /* correct Bin for boundary cases */
	  while(Bin > 0 && interval < resbiasX[c][Bin])
	    Bin--;
	  if(DEBUG) assert(0 <= Bin && Bin <= ResBins[c]);
	  if(DEBUG && Bin > 0 && !(resbiasX[c][Bin]  <= interval)){
	    printf("interval = %0.8f, isize = %d, SizeToBin[c][isize]=%d, Bin= %d, resbiasX[c][%d]=%0.8f, resbiasX[c][%d]=%0.8f,resbiasX[c][%d]=%0.8f\n",
		   interval, isize, SizeToBin[c][isize], Bin, max(0,Bin-1), resbiasX[c][max(0,Bin-1)], Bin, resbiasX[c][Bin], min(Bin+1,ResBins[c]),resbiasX[c][min(Bin+1,ResBins[c])]);
	    fflush(stdout);
	    assert(resbiasX[c][Bin] <= interval);
	  }
	  if(DEBUG && Bin < ResBins[c]) assert(interval < resbiasX[c][Bin+1]);

	  double bias = Offset[c][Bin] + Slope[c][Bin] * interval;
	  shift += bias;

	  if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
	    #pragma omp critical
	    {
	      printf("    c=%d,I=%d,rawX[I]=%0.10f:rawX[J=%d]=%0.10f,bias=%0.10f,X[I]-> %0.10f(shift=%0.10f)\n", c, I, RawX[I], J, RawX[J], bias, site + shift*0.5,shift*0.5);
	      fflush(stdout);
	    }
	  }
	}
	if(BIASCORRECT_FIX && I > 1 && I < M){
	  double maxshift = (RawX[I-1] + RawX[I+1])*0.5 - site;
	  if(maxshift > 0.0){
	    shift = min(shift,maxshift);
	    shift = max(shift, 0.0);
	  } else {
	    shift = min(shift, 0.0);
	    shift = max(shift, maxshift);
	  }
	}
	site += shift * 0.5;

	if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
          #pragma omp critical
	  {
	    printf("mapid=%d:c=%d,I=%d,rawX[I]=%0.10f,X[I]=%0.10f->%0.10f\n",mapid,c,I,RawX[I],X[I],site);
	    fflush(stdout);
	  }
	}

	X[I] = site;
      }
	  
      if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
        #pragma omp critical
	{
	  printf("mapid=%d:c=%d,I=M+1=%d,rawX[I]=%0.10f,X[I]=%0.10f->%0.10f\n",mapid,c,M+1,RawX[M+1],X[M+1],max(X[M]+MININTERVAL,RawX[M+1]));
	  fflush(stdout);
	}
      }
      X[M+1] = max(X[M]+MININTERVAL,RawX[M+1]);
      if(X[1] < MININTERVAL){
	FLOAT shift = MININTERVAL - X[1];

	for(int I = 1; I <= M+1; I++)
	  X[I] += shift;

	if(BIAS_TRACE && pmap->id == BIAS_TRACE_ID){
          #pragma omp critical
	  {
	    printf("mapid=%d:c=%d:shifted all X[I] by %0.10f due to first site being <= MININTERVAL=%0.6f\n",mapid,c,shift,MININTERVAL);
	    fflush(stdout);
	  }
	}
	if(colors==2){/* also shift the other color */
          int c2 = (c ? 1 : 0);
	  int M2 = pmap->numsite[c2];
	  FLOAT *X2 = pmap->site[c2];
	  for(int I = 1; I <= M2+1; I++)
	    X2[I] += shift;
        }
      }
    }
    if(colors >= 2){  /* equalize lengths of 2 colors */
      FLOAT len = 0.0;
      for(int c = 0; c < colors; c++){
        int M = pmap->numsite[c];
        len = max(len, pmap->site[c][M+1]);
      }
      for(int c = 0; c < colors; c++){
        int M = pmap->numsite[c];
	pmap->site[c][M+1] = len;
      }
    }
    if(DEBUG/* HERE >=2 */){/* check that sites are still monotonic */
      for(int c = 0; c < colors; c++){
	int M = pmap->numsite[c];
	for(int I = 1; I <= M+1; I++)
	  if(DEBUG && !(pmap->site[c][I] >= pmap->site[c][I-1])){
            #pragma omp critical
	    {
	      printf("After bias removal:pmap->id=%lld,c=%d,I=%d,M=%d:pmap->site[c][I-1]=%0.6f,pmap->site[c][I]=%0.6f\n",
		     pmap->id,c,I,M,pmap->site[c][I-1],pmap->site[c][I]);
	      for(int t = 1; t <= M; t++)
		printf("t=%d:rawsite[c][t]=%0.6f, site[c][t]= %0.6f\n",t,pmap->rawsite[c][t],pmap->site[c][t]);
	      fflush(stdout);
	      assert(pmap->site[c][I] >= pmap->site[c][I-1]);
	    }
	  }
      }
    }
   }/* mapid = orignumaps .. mapid */

   #pragma omp critical
   {
     if(myMinInterval < minInterval){
       minInterval = myMinInterval;
       minMapid = myMapid;
     }
   }
  } // omp parallel 

  if(VERB>=2){
    printf("Minimum interval = %0.8f in mapid=%d,id=%lld\n",minInterval, minMapid, Gmap[minMapid]->id);
    fflush(stdout);
  }

  colors = origcolors;
}

/* create pairmerge map from Ymap,Xmap in orientation Yor,Xor with crossover point at labels I,J respectively
   assumes .cmap input : only site[],siteSD[],sitecov[],sitecnt[],ChimQuality[](if present) are updated */
Cmap::Cmap(Cmap *Ymap, Cmap *Xmap, int Yor, int Xor, int I, int J)
{
  if(!(Ymap->siteSD[0] && Ymap->sitecov[0] && Ymap->sitecnt[0])){
	printf("Cmap pairmerge constructor is only implemented for Contig maps (from CMAP input): this map (id=%lld) is from %s\n",
	   Ymap->id, vfixx_filename[Ymap->fileid]);
    fflush(stdout);exit(1);
  }
  if(!(Xmap->siteSD[0] && Xmap->sitecov[0] && Xmap->sitecnt[0])){
	printf("Cmap pairmerge constructor is only implemented for Contig maps (from CMAP input): this map (id=%lld) is from %s\n",
	   Xmap->id, vfixx_filename[Xmap->fileid]);
    fflush(stdout);exit(1);
  }

  if(colors != 1){
    printf("Cmap(Ymap,Xmap,...) not implemented for colors=%d\n",colors);
    fflush(stdout);exit(1);
  }

  int N = Ymap->numsite[0];
  int M = Xmap->numsite[0];
  double Ylen = Ymap->site[0][N+1];
  //  double Xlen = Xmap->site[0][M+1];

  if(Yor){/* flip Ymap */
    I = N+1 - I;/* code from pairmerge() in pairalign.cpp (used below) assumes LabelID (I) is from left end as oriented, while input specifies label ID in normal orientation */
    Ymap->inversion(0.0, Ylen);
    Yor = 0;
  }

  if(Xor)
    J = M+1-J;/* convert from actual label ID J to label count from left end as oriented, which is the convention used in pairmerge() in pairalign.cpp */

  init();/* default initialization of map */

  FLOAT *Y = Ymap->site[0];
  FLOAT *X = Xmap->site[0];

  /* compute number of sites and allocate memory */
  int NumSites = I + (Xor ? M+1-J-1 : M - J);
  
  for(int c = 0; c < colors; c++){
    if( site[c]){
      if(! blockmem) delete []  site[c];
      if( siteSD[c]) delete []  siteSD[c];
      if( sitecov[c]) delete []  sitecov[c];
      if( sitecnt[c]) delete []  sitecnt[c];
      if(CmapChimQuality >= 1){
	if(ChimQuality[c]) delete []  ChimQuality[c];
	if(CmapChimQuality >= 2){
	  if(ChimNorm[c]) delete []  ChimNorm[c];
	  if(SegDupL[c]) delete []  SegDupL[c];
	  if(SegDupR[c]) delete []  SegDupR[c];
	  if(FragileEndL[c]) delete []  FragileEndL[c];
	  if(FragileEndR[c]) delete []  FragileEndR[c];
	  if(OutlierFrac[c]) delete []  OutlierFrac[c];
	  if(CmapChimQuality >= 3){
	    if(FragSd[c]) delete []  FragSd[c];
	    if(ExpSd[c]) delete []  ExpSd[c];
	    if(FragCov[c]) delete []  FragCov[c];
	    if(FragChiSq[c]) delete []  FragChiSq[c];
	  }
	}
      }
    }
    if(SNRcnt[c]){
      if(SNRgmean[c]) delete []  SNRgmean[c];
      if(lnSNRsd[c]) delete []  lnSNRsd[c];
      if(SNRdist[c]){
	for(int i = 0; i <=  numsite[0]+1; i++)
	  if( SNRdist[c][i]) delete []  SNRdist[c][i];
	delete []  SNRdist[c];
      }
      delete []  SNRcnt[c];
    }
     numsite[c] = NumSites;
     site[c] = new FLOAT[NumSites+2];
     siteSD[c] = new double[NumSites+2];
     sitecov[c] = new float[NumSites+2];
     sitecnt[c] = new float[NumSites+2];
     if(CmapChimQuality >= 1){
       ChimQuality[c] = new float[NumSites+2];
       if(CmapChimQuality >= 2){
	 ChimNorm[c] = new float[NumSites+2];
	 SegDupL[c] = new float[NumSites+2];
	 SegDupR[c] = new float[NumSites+2];
	 FragileEndL[c] = new float[NumSites+2];
	 FragileEndR[c] = new float[NumSites+2];
	 OutlierFrac[c] = new float[NumSites+2];
	 if(CmapChimQuality >= 3){
	   FragSd[c] = new float[NumSites+2];
	   ExpSd[c] = new float[NumSites+2];
	   FragCov[c] = new float[NumSites+2];
	   FragChiSq[c] = new float[NumSites+2];
	 }
       }
     }
     SNRcnt[c] = 0;
     if(Ymap->SNRcnt[c] || Xmap->SNRcnt[c]){
       SNRcnt[c] = new int[NumSites+2];
       SNRgmean[c] = new double[NumSites+2];
       lnSNRsd[c] = new double[NumSites+2];
       SNRdist[c] = new double*[NumSites+2];
     }

    /* left end */
     site[c][0] = 0.0;
     siteSD[c][0] = 0.0;
     sitecov[c][0] = 1;
     sitecnt[c][0] = 1;
     if(CmapChimQuality >= 1){
       ChimQuality[c][0] = END_CHIMQUAL;
       if(CmapChimQuality >= 2){
	 ChimNorm[c][0] = END_CHIMQUAL;
	 SegDupL[c][0] = END_CHIMQUAL;
	 SegDupR[c][0] = END_CHIMQUAL;
	 FragileEndL[c][0] = 0.0;
	 FragileEndR[c][0] = 0.0;
	 OutlierFrac[c][0] = 0.0;
	 if(CmapChimQuality >= 3){
	   FragSd[c][0] = 0.0;
	   ExpSd[c][0] = 0.0;
	   FragCov[c][0] = 0.0;
	   FragChiSq[c][0] = 1.0;
	 }
       }
     }
     if( SNRcnt[c]){
       SNRcnt[c][0] = 0;
       SNRgmean[c][0] = 0.0;
       lnSNRsd[c][0] = 0.0;
       SNRdist[c][0] = NULL;
     }

     /* right end */
     siteSD[c][NumSites+1] = 0.0;
     sitecov[c][NumSites+1] = 1;
     sitecnt[c][NumSites+1] = 1;
     if(CmapChimQuality >= 1){
       ChimQuality[c][NumSites+1] = END_CHIMQUAL;
       if(CmapChimQuality >= 2){
	 ChimNorm[c][NumSites+1] = END_CHIMQUAL;
	 SegDupL[c][NumSites+1] = END_CHIMQUAL;
	 SegDupR[c][NumSites+1] = END_CHIMQUAL;
	 FragileEndL[c][NumSites+1] = 0.0;
	 FragileEndR[c][NumSites+1] = 0.0;
	 OutlierFrac[c][NumSites+1] = 0.0;
	 if(CmapChimQuality >= 3){
	   FragSd[c][NumSites+1] = 0.0f;
	   ExpSd[c][NumSites+1] = 0.0f;
	   FragCov[c][NumSites+1] = 0.0f;
	   FragChiSq[c][NumSites+1] = 1.0f;
	 }
       }
     }
    if(SNRcnt[c]){
      SNRcnt[c][NumSites+1] = 0;
      SNRgmean[c][NumSites+1] = 0.0;
      lnSNRsd[c][NumSites+1] = 0.0;
      SNRdist[c][NumSites+1] = NULL;
    }
  }/* c = 0 .. colors-1 */

  blockmem = 0;
  register FLOAT *Z = site[0];

  /* copy left end of Ymap from site 0 to site I */
  register int i = 0;
  for(; i <= I; i++){
    Z[i] = Y[i];
    if(DEBUG && i > 0) assert(Z[i] >= Z[i-1]);
    siteSD[0][i] = Ymap->siteSD[0][i];
    sitecov[0][i] = Ymap->sitecov[0][i];
    sitecnt[0][i] = Ymap->sitecnt[0][i];
    if(CmapChimQuality >= 1){
      ChimQuality[0][i] = Ymap->ChimQuality[0] ? Ymap->ChimQuality[0][i] : END_CHIMQUAL;
      if(CmapChimQuality >= 2){
	ChimNorm[0][i] = Ymap->ChimNorm[0] ? Ymap->ChimNorm[0][i] : END_CHIMQUAL;	
	SegDupL[0][i] = Ymap->SegDupL[0] ? Ymap->SegDupL[0][i] : END_CHIMQUAL;	
	SegDupR[0][i] = Ymap->SegDupR[0] ? Ymap->SegDupR[0][i] : END_CHIMQUAL;	
	FragileEndL[0][i] = Ymap->FragileEndL[0] ? Ymap->FragileEndL[0][i] : 0.0;
	FragileEndR[0][i] = Ymap->FragileEndR[0] ? Ymap->FragileEndR[0][i] : 0.0;
	OutlierFrac[0][i] = Ymap->OutlierFrac[0] ? Ymap->OutlierFrac[0][i] : 0.0;
	if(CmapChimQuality >= 3){
	  FragSd[0][i] = Ymap->FragSd[0] ? Ymap->FragSd[0][i] : 0.0;
	  ExpSd[0][i] = Ymap->ExpSd[0] ? Ymap->ExpSd[0][i] : 0.0;
	  FragCov[0][i] = Ymap->FragCov[0] ? Ymap->FragCov[0][i] : 0.0;
	  FragChiSq[0][i] = Ymap->FragChiSq[0] ? Ymap->FragChiSq[0][i] : 1.0;
	}
      }
    }
    if(Ymap->SNRcnt[0]){
      SNRcnt[0][i] = Ymap->SNRcnt[0][i];
      SNRgmean[0][i] = Ymap->SNRgmean[0][i];
      lnSNRsd[0][i] = Ymap->lnSNRsd[0][i];
      SNRdist[0][i] = NULL;
      if(SNRcnt[0][i]){
	SNRdist[0][i] = new double[SNRcnt[0][i]];
	for(int t = 0; t < SNRcnt[0][i]; t++)
	  SNRdist[0][i][t] = Ymap->SNRdist[0][i][t];
      }
    } else if(SNRcnt[0]){
      SNRcnt[0][i] = 0;
      SNRgmean[0][i] = 0.0;
      lnSNRsd[0][i] = 0.0;
      SNRdist[0][i] = NULL;
    }
    if(VERB>=3){
      printf("Copying site Y[i=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
	     i,N, Y[i],Ymap->id, i, numsite[0], Z[i], sitecov[0][i], sitecnt[0][i]);
      fflush(stdout);
    }
  }

  if(Xor){/* copy right end of Xmap from M+1-J downto 0 */
    register int k = M+1-J;
    register double Xstart = Z[i-1] + X[k];
    for(register int j = k; --j >= 0; i++){
      Z[i] = Xstart - X[j];
      if(DEBUG) assert(Z[i] >= Z[i-1]);
      siteSD[0][i-1] = Xmap->siteSD[0][j];/* see definition of siteSD[] */
      sitecov[0][i] = Xmap->sitecov[0][j];
      sitecnt[0][i] = Xmap->sitecnt[0][j];
      if(CmapChimQuality >= 1){
	ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][j] : END_CHIMQUAL;
	if(CmapChimQuality >= 2){
	  ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][j] : END_CHIMQUAL;
	  SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][j] : END_CHIMQUAL;
	  SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][j] : END_CHIMQUAL;
	  FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][j] : 0.0;
	  FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][j] : 0.0;
	  OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][j] : 0.0;
	  if(CmapChimQuality >= 3){
	    FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][j] : 0.0;
	    ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][j] : 0.0;
	    FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][j] : 0.0;
	    FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][j] : 1.0;
	  }
	}
      }
      if(Xmap->SNRcnt[0]){
	SNRcnt[0][i] = Xmap->SNRcnt[0][j];
	SNRgmean[0][i] = Xmap->SNRgmean[0][j];
	lnSNRsd[0][i] = Xmap->lnSNRsd[0][j];
	SNRdist[0][i] = NULL;
	if(SNRcnt[0][i]){
	  SNRdist[0][i] = new double[SNRcnt[0][i]];
	  for(int t = 0; t < SNRcnt[0][i]; t++)
	    SNRdist[0][i][t] = Xmap->SNRdist[0][j][t];
	}
      } else if(SNRcnt[0]){
	SNRcnt[0][i] = 0;
	SNRgmean[0][i] = 0.0;
	lnSNRsd[0][i] = 0.0;
	SNRdist[0][i] = NULL;
      }
      if(VERB>=3){
	printf("Copying site X[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
	       j, M, X[j],Xmap->id, i, numsite[0], Z[i], sitecov[0][i], sitecnt[0][i]);
	fflush(stdout);
      }
    }
    if(DEBUG && !(i-1 == NumSites+1)){
      printf("N=%d,M=%d,N=%d,M=%d,I=%d,J=%d,i=%d,NumSites=%d,orientation=%d\n",
	     N,M,N,M,I,J,i,NumSites,Xor);
      fflush(stdout);
      assert(i-1 == NumSites+1);
    }
  } else {/* copy right end of Xmap from  J to M+1 */
    int k = J;
    double Xstart = Z[i-1] - X[k];
    for(int j = k; ++j <= M+1; i++){
      Z[i] = Xstart + X[j];
      if(DEBUG) assert(Z[i] >= Z[i-1]);
      siteSD[0][i-1] = Xmap->siteSD[0][j-1];/* see definition of siteSD[] */
      sitecov[0][i] = Xmap->sitecov[0][j];
      sitecnt[0][i] = Xmap->sitecnt[0][j];
      if(CmapChimQuality >= 1){
	ChimQuality[0][i] = Xmap->ChimQuality[0] ? Xmap->ChimQuality[0][j] : END_CHIMQUAL;
	if(CmapChimQuality >= 2){
	  ChimNorm[0][i] = Xmap->ChimNorm[0] ? Xmap->ChimNorm[0][j] : END_CHIMQUAL;
	  SegDupL[0][i] = Xmap->SegDupL[0] ? Xmap->SegDupL[0][j] : END_CHIMQUAL;
	  SegDupR[0][i] = Xmap->SegDupR[0] ? Xmap->SegDupR[0][j] : END_CHIMQUAL;
	  FragileEndL[0][i] = Xmap->FragileEndL[0] ? Xmap->FragileEndL[0][j] : 0.0;
	  FragileEndR[0][i] = Xmap->FragileEndR[0] ? Xmap->FragileEndR[0][j] : 0.0;
	  OutlierFrac[0][i] = Xmap->OutlierFrac[0] ? Xmap->OutlierFrac[0][j] : 0.0;
	  if(CmapChimQuality >= 3){
	    FragSd[0][i] = Xmap->FragSd[0] ? Xmap->FragSd[0][j] : 0.0;
	    ExpSd[0][i] = Xmap->ExpSd[0] ? Xmap->ExpSd[0][j] : 0.0;
	    FragCov[0][i] = Xmap->FragCov[0] ? Xmap->FragCov[0][j] : 0.0;
	    FragChiSq[0][i] = Xmap->FragChiSq[0] ? Xmap->FragChiSq[0][j] : 1.0;
	  }
	}
      }
      if(Xmap->SNRcnt[0]){
	SNRcnt[0][i] = Xmap->SNRcnt[0][j];
	SNRgmean[0][i] = Xmap->SNRgmean[0][j];
	lnSNRsd[0][i] = Xmap->lnSNRsd[0][j];
	SNRdist[0][i] = NULL;
	if(SNRcnt[0][i]){
	  SNRdist[0][i] = new double[SNRcnt[0][i]];
	  for(int t = 0; t < SNRcnt[0][i]; t++)
	    SNRdist[0][i][t] = Xmap->SNRdist[0][j][t];
	}
      } else if(SNRcnt[0]){
	SNRcnt[0][i] = 0;
	SNRgmean[0][i] = 0.0;
	lnSNRsd[0][i] = 0.0;
	SNRdist[0][i] = NULL;
      }
      if(VERB >=3){
	printf("Copying site X[j=%d/%d]=%0.3f (id=%lld) to site Z[i=%d/%d]=%0.3f:sitecov=%0.1f,sitecnt=%0.1f\n",
	       j, M, X[j],Xmap->id, i, numsite[0], Z[i], sitecov[0][i], sitecnt[0][i]);
	fflush(stdout);
      }
    }
    if(DEBUG) assert(i-1 == NumSites+1);
  }
}

