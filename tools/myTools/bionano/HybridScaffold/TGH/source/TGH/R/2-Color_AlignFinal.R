#This script create final alignments of NGS sequence to the output assemblies from sandwich
#library(data.table)
#source("~/workspace/codebase/MoleculeSimulator/IO/readmaps.R")
#source("~/workspace/codebase/MoleculeSimulator/IO/exportMaps.r")
#source("~/workspace/codebase/MoleculeSimulator/HybridScaffoldTwo/sandwichScript/ConflictResolveTwo.R")

#merge alignment from multiple sources, using confidence score to pick alignment from each source
merge.align <- function(align1, align2, use.confidence=T){
  if(is.null(align1) || is.null(align2)){
    stop("Input alignment to cannot be null")
  }
  if(nrow(align1) == 0 || nrow(align2) == 0){
    return(rbind(align1, align2))
  }
  align1$Confidence <- align1$Confidence + 0.001 #breaking tie in confidence score, in subsequent step of picking the best alignments
  align.merge <- rbind(align1, align2)
  setkeyv(align.merge, c("QryContigID", "Confidence")) #need to sort it in order of the selector below
  best.score <- align.merge[,Confidence == max(Confidence), by=QryContigID] #picking best alignment
  align.merge <- align.merge[best.score$V1==T]
  return(align.merge)
}

#this function filter align.final from single-enzyme hybrid in order to get a
#set of high-quality alignment for sandwich-scaffolding
filterAlignFinal1 <- function(align.file1, align.file2, cmap.file1, cmap.file2, output.path){
  min.align.overlap <- 50000
  #align.file1 <- '~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/BSPQI/align_final_M1/NoBestRef/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_soap_result_scafSeq_fa_NGScontigs_HYBRID_SCAFFOLD.xmap'
  #align.file2 <- '~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/BSSSI/align_final_M1/NoBestRef/BSSSI_EXP_REFINEFINAL1_bppAdjust_cmap_soap_result_scafSeq_fa_NGScontigs_HYBRID_SCAFFOLD.xmap'
  align1 <- as.data.table(readxmap(align.file1))
  align2 <- as.data.table(readxmap(align.file2))
  maps1 <- as.data.table(readcmap(cmap.file1))
  maps2 <- as.data.table(readcmap(cmap.file2))

  #filtering by score and uniqueness
  align1.filtered <- remove.duplicated.align(align1, delta.score = 3, align.len = 0.00001, bestRef = F)
  align2.filtered <- remove.duplicated.align(align2, delta.score = 3, align.len = 0.00001, bestRef = F)
  align1.keep <- align1.filtered[,length(Confidence) == 1, by=QryContigID]
  align2.keep <- align2.filtered[,length(Confidence) == 1, by=QryContigID]
  print(paste0("Filtering out alignment from 1st-enzyme due to alignments to more than one hybrids: ", sum(!align1.keep$V1)))
  print(paste0("Filtering out alignment from 2nd-enzyme due to alignments to more than one hybrids: ", sum(!align2.keep$V1)))
  align1.filtered <- align1.filtered[QryContigID %in% align1.keep[V1==TRUE, QryContigID]]
  align2.filtered <- align2.filtered[QryContigID %in% align2.keep[V1==TRUE, QryContigID]]

  #filtering alignment with sticky-end
  align1.filtered2 <- filterAlignWithEndOultier(align1.filtered)
  align2.filtered2 <- filterAlignWithEndOultier(align2.filtered)
  print(paste0("Filtering out alignment from 1st-enzyme with endOutlier: ", nrow(align1.filtered) - nrow(align1.filtered2)))
  print(paste0("Filtering out alignment from 2nd-enzyme with endOutlier: ", nrow(align2.filtered) - nrow(align2.filtered2)))

  #filtering alignment in segdup regions
  align1.filtered3 <- filterAlignInSegDup(align1.filtered2, maps1)
  align2.filtered3 <- filterAlignInSegDup(align2.filtered2, maps2)
  print(paste0("Filtering out alignment from 1st-enzyme overlapping with segDups: ", nrow(align1.filtered2) - nrow(align1.filtered3)))
  print(paste0("Filtering out alignment from 2nd-enzyme overlapping with segDups: ", nrow(align2.filtered2) - nrow(align2.filtered3)))

  align.merged <- merge(align1.filtered3, align2.filtered3, by='QryContigID')
  align.overlap <- align.merged[,pCompareInterval(QryStartPos.x, QryEndPos.x, QryStartPos.y, QryEndPos.y)]
  print(paste0("Filtering out alignment pairs due to  non-overlap of aligned region: ", sum(align.overlap$dist > -1*min.align.overlap)))
  align.merged <- cbind(align.merged, align.overlap)
  align.merged <- align.merged[align.overlap$dist < -1*min.align.overlap]
  align1.filtered <- align1.filtered[QryContigID %in% align.merged$QryContigID]
  align2.filtered <- align2.filtered[QryContigID %in% align.merged$QryContigID]

  #filtering out connected pairs
  aligns.filter <- filterSparsePairs(align1.filtered, align2.filtered)
  print(paste0("Total pairs ", nrow(aligns.filter$pairstats), " filter out ", nrow(aligns.filter$pairstats[fr < 0.1]), " due to weak sparse glue alignment"))

  align1.filtered <- aligns.filter[[1]]
  align2.filtered <- aligns.filter[[2]]
  #exportXmap(align1.filtered, "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/BSPQI/hybrid_scaffolds_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_soap_result_scafSeq_fa_NGScontigs_HYBRID_SCAFFOLD.xmap")
  #exportXmap(align2.filtered, "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/BSSSI/hybrid_scaffolds_M1/BSSSI_EXP_REFINEFINAL1_bppAdjust_cmap_soap_result_scafSeq_fa_NGScontigs_HYBRID_SCAFFOLD.xmap")
  return(list(align1.filtered, align2.filtered))
}

#filtering alignment with big end outliers
filterAlignWithEndOultier <- function(align, max.end.len=50e3, max.frac=0.1){
  EndLen1 <- align[,pmin(pmin(QryStartPos, QryEndPos), RefStartPos)]
  EndLen2 <- align[,pmin(QryLen - pmax(QryStartPos, QryEndPos), RefLen - RefEndPos)]
  align.filtered <- align[!(EndLen1 > max.end.len & EndLen1 / QryLen > max.frac |  EndLen2 > max.end.len & EndLen2/QryLen > max.frac),]
  return(align.filtered)
}


#this filters out alignment that falls into segDup regions masked by refAligner
filterAlignInSegDup <- function(align, cmap, map.type='Ref', min.non.segdup=80e3){
  if(map.type=='Qry'){
    align_dt <- data.table(XmapID = align$XmapEntryID, start=align$QryStartPos, end=align$QryEndPos, mapId=align$QryContigID)
  }else if(map.type=='Ref'){
    align_dt <- data.table(XmapID = align$XmapEntryID, start=align$RefStartPos, end=align$RefEndPos, mapId=align$RefcontigID)
  }else{
    stop(paste0("Unrecognize map type: ", map.type, ".", " Use Ref/Qry as map-tyep only"))
  }

  segDupRegions <- getSegDupRegions(cmap)

  is_segdup <- unlist(lapply(1:nrow(align_dt),
                        function(i){
                            segDup <- segDupRegions[CMapId==align_dt$mapId[[i]]]
                            if(nrow(segDup) > 0){
                              overlap <- pCompareInterval(segDup$startPos, segDup$endPos, align_dt$start[[i]], align_dt$end[[i]])
                              if(any(overlap$dist < 0)){
                                return(abs(align_dt$start[i] - align_dt$end[i]) - abs(min(overlap$dist)) < min.non.segdup)
                              }
                            }
                            return(FALSE)
                        }
                      ))
  #print(paste0("Number of alignments overlap with segdup: ", sum(is_segdup)))
  align_filter <- align[!is_segdup,]
  return(align_filter)
}


getSegDupRegions <- function(cmap){
   if(is.null(cmap$Mask)){
     return(data.table(RefContigID=numeric(), start=numeric(), end=numeric(), startPos=numeric(), endPos=numeric(), CMapId=numeric(), SegDupId=numeric()))
   }
   cmap <- as.data.table(cmap)
   segdup <- cmap[Mask != 0,]
   same_contig <- segdup[,c(FALSE, diff(as.numeric(CMapId)) == 0)]
   is_consecutive <- segdup[,c(FALSE, diff(SiteID)==1)]
   is_start <- !(same_contig & is_consecutive)
   segdup$SegDupId <- cumsum(is_start)
   segdup_regions <- segdup[, .(begin=min(SiteID), end=max(SiteID), startPos=min(Position), endPos=max(Position), CMapId=CMapId[[1]]), by=SegDupId]
   return(segdup_regions)
}

#for each pair connected maps check the length  of the glue-alignment with repect to the length of the ref maps
#if only small fraction of the overlap regions are aligned, it indicates this pairing is not reliable
filterSparsePairs <- function(align1, align2, min.frac = 0.1){
  indata <- merge(align1, align2, by='QryContigID')
  indata$id_str <- paste0(indata$RefcontigID.x, "_", indata$RefcontigID.y)
  pos_offsets <- computeOffset(indata)
  stats <- computeRelativeAlignRegions(pos_offsets, indata, valid.region = T)
  L1 <- stats$RelEnd1 - stats$RelStart1
  L2 <- stats$RelEnd2 - stats$RelStart2
  stats$OverlapL <- ifelse(L1 < L2, L1, L2)
  AlignL <- unlist(lapply(1:nrow(stats), function(i){
        pair <- stats$id_str[[i]]
        aligned.len <- indata[id_str==pair, mean(sum(abs(RefEndPos.x - RefStartPos.x)), sum(abs(RefEndPos.y - RefStartPos.y)))]
  }))
  stats$AlignL <- AlignL
  stats$fr <- stats$AlignL/stats$OverlapL
  #hist(stats$fr)
  aligns <- indata[id_str %in% stats[fr > min.frac, id_str]]
  align1.filter <- align1[QryContigID %in% aligns$QryContigID]
  align2.filter <- align2[QryContigID %in% aligns$QryContigID]
  return(list(align1=align1.filter, align2=align2.filter, pairstats=stats))
}

computeOffset <- function(indata){
  computeOffset_helper <- function(data){
    final_slope <- ifelse(data[,Orientation.x][[1]] == data[,Orientation.y][[1]], 1, -1)
    dir <- ifelse(data[,Orientation.x]=='+', 1, -1)
    intercept <- mean(c(data[,RefStartPos.y] - final_slope*(data[,RefStartPos.x] - (data[,QryStartPos.x] - data[,QryStartPos.y])*dir),
                      data[,RefEndPos.y] - final_slope*(data[,RefEndPos.x] - (data[,QryEndPos.x] - data[,QryEndPos.y])*dir)))
    list(slope=final_slope, intercept=intercept)
  }
  offsets <- indata[, computeOffset_helper(.SD), by=id_str]
  return(offsets)
}

computeRelativeAlignRegions <- function(pairstats, align.merged, valid.region=F){
  pairstats$slope <- as.numeric(pairstats$slope)
  pairstats$intercept <- as.numeric(pairstats$intercept)
  id_str <- unique(align.merged$id_str)
  setkey(pairstats, id_str)
  pairstats <- pairstats[id_str]
  pairstats <- pairstats[!is.na(intercept)]
  setkey(align.merged, id_str)
  pairstats$RefLen.x <- align.merged[pairstats$id_str, RefLen.x, mult="first"]
  pairstats$RefLen.y <- align.merged[pairstats$id_str, RefLen.y, mult="first"]
  pairstats$RefcontigID.x <- align.merged[pairstats$id_str, RefcontigID.x, mult="first"]
  pairstats$RefcontigID.y <- align.merged[pairstats$id_str, RefcontigID.y, mult="first"]
  pairstats$RelStart1 <- pairstats[,(1-intercept)/slope]
  pairstats$RelEnd1 <- pairstats[,(RefLen.y-intercept)/slope]
  pairstats$RelStart2 <- pairstats[,slope*(1)+intercept]
  pairstats$RelEnd2 <- pairstats[,slope*(RefLen.x)+intercept]

  pairstats.ordered <- flip.columns(pairstats, 'RelStart1', 'RelEnd1', `>`)
  pairstats.ordered <- flip.columns(pairstats.ordered, 'RelStart2', 'RelEnd2', `>`)
  stats <- pairstats.ordered
  if(valid.region){
    stats$RelStart1 <- ifelse(stats$RelStart1 <= 0, 1, stats$RelStart1)
    stats$RelStart2 <- ifelse(stats$RelStart2 <= 0, 1, stats$RelStart2)
    stats$RelEnd1 <- ifelse(stats$RelEnd1 <= 0, 1, stats$RelEnd1)
    stats$RelEnd2 <- ifelse(stats$RelEnd2 <= 0, 1, stats$RelEnd2)

    stats$RelStart1 <- ifelse(stats$RelStart1 > stats$RefLen.x, stats$RefLen.x, stats$RelStart1)
    stats$RelStart2 <- ifelse(stats$RelStart2 > stats$RefLen.y, stats$RefLen.y, stats$RelStart2)
    stats$RelEnd1 <- ifelse(stats$RelEnd1 > stats$RefLen.x, stats$RefLen.x, stats$RelEnd1)
    stats$RelEnd2 <- ifelse(stats$RelEnd2 > stats$RefLen.y, stats$RefLen.y, stats$RelEnd2)
  }
  return(stats)
}



#this function covert a two-color alignment to a single entry for exporting the scaffold
#only the alignment boundary stats and scores etc are valid, the actual alignment string is only
#for single-color
getExportTwoColorAlign <- function(xmap.twocolor){
  xmap1 <- xmap.twocolor[LabelChannel == 1 & QryStartPos != QryEndPos]
  xmap2 <- xmap.twocolor[LabelChannel == 2 & QryStartPos != QryEndPos & !(QryContigID %in% xmap1$QryContigID)]
  #xmap1 <- xmap.twocolor[LabelChannel==1 ]
  #xmap2 <- xmap.twocolor[LabelChannel == 2]
  #browser()
  #xmap1[,QryStartPos := ifelse(QryStartPos < xmap2$QryStartPos, QryStartPos, xmap2$QryStartPos)]
  #xmap1[,QryEndPos := ifelse(QryEndPos > xmap2$QryEndPos, QryEndPos, xmap2$QryEndPos)]
  #xmap1[,RefStartPos := ifelse(RefStartPos < xmap2$RefStartPos, RefStartPos, xmap2$RefStartPos)]
  #xmap1[,RefEndPos := ifelse(RefEndPos > xmap2$RefEndPos, RefEndPos, xmap2$RefEndPos)]
  #xmap1[,QryStartPos := ifelse(nchar(Alignment) > nchar(xmap2$Alignment), QryStartPos, xmap2$QryStartPos)]
  #xmap1[,QryEndPos := ifelse(nchar(Alignment) > nchar(xmap2$Alignment), QryEndPos, xmap2$QryEndPos)]
  #xmap1[,RefStartPos := ifelse(nchar(Alignment) > nchar(xmap2$Alignment), RefStartPos, xmap2$RefStartPos)]
  #xmap1[,RefEndPos := ifelse(nchar(Alignment) > nchar(xmap2$Alignment), RefEndPos, xmap2$RefEndPos)]
  #A1 <- xmap1[,QryStartPos!=QryEndPos]
  #A2 <- xmap2[,QryStartPos != QryEndPos]
  #xmap1[, QryStartPos := (QryStartPos + xmap2$QryStartPos)/2]
  #xmap1[, QryEndPos :=   (QryEndPos + xmap2$QryEndPos)/2]
  #xmap1[, RefStartPos := (RefStartPos + xmap2$RefStartPos)/2]
  #xmap1[, RefEndPos := (RefEndPos + xmap2$RefEndPos)/2]
  print(paste0("number of align from channels 1 : ", nrow(xmap1), " 2: ", nrow(xmap2)))
  xmap1 <- rbind(xmap1, xmap2)
  #print(paste0("number of alignemnts ", nrow(xmap1)))
  return(xmap1)
}

combine_twocolor_aligns <- function(aligns_twocolor1, aligns_twocolor2){
  aligns_twocolor1[, QryStartPos := min(QryStartPos), by=XmapEntryID]
  aligns_twocolor1[, QryEndPos := max(QryEndPos), by=XmapEntryID]

  aligns_twocolor2[, QryStartPos := min(QryStartPos), by=XmapEntryID]
  aligns_twocolor2[, QryEndPos := max(QryEndPos), by=XmapEntryID]

  aligns_twocolor <- rbind(aligns_twocolor1[!(RefcontigID %in% aligns_twocolor2$RefcontigID)], aligns_twocolor2)
  aligns_twocolor <- aligns_twocolor[order(Confidence, decreasing =  T)]
  aligns_twocolor[, Conf_delta := get_delta_score(Confidence), by=RefcontigID]
  aligns_twocolor[, maxScore := max(Confidence), by=RefcontigID]
  aligns_twocolor$is_overlap <- NA
  return(aligns_twocolor)
}

#this function filter the two-color alignments base on single-color alignments as well as
#the two-color alignments
xmap_filter_twocolor <- function(aligns_single1, aligns_single2, aligns_twocolor1, aligns_twocolor2, two.color.file, two.color.file2, Confidence.T=10){
  if(nrow(aligns_twocolor1) == 0 && nrow(aligns_twocolor2) == 0){
    return(aligns_twocolor1)
  }
  #filtering one color alignment
  aligns_single1[, maxScore := max(Confidence), by=RefcontigID]
  aligns_single1 <- get_ref_align_xmap_filter(aligns_single1, cleanUp = F)
  aligns_single2[, maxScore := max(Confidence), by=RefcontigID]
  aligns_single2 <- get_ref_align_xmap_filter(aligns_single2, cleanUp = F)
  aligns_single1$Colors <- 1
  aligns_single2$Colors <- 1
  #combining twocolor alignment
  aligns_twocolor <- combine_twocolor_aligns(aligns_twocolor1, aligns_twocolor2)
  aligns_twocolor1 <- aligns_twocolor[LabelChannel == 1]
  aligns_twocolor1$Colors <- 2

  aligns_onecolor <- rbind(aligns_single1, aligns_single2[!(RefcontigID %in% aligns_single1$RefcontigID)])
  aligns_onecolor <- aligns_onecolor[order(Confidence, decreasing = T)]

  aligns_combined <- rbind(aligns_onecolor, aligns_twocolor1)
  aligns_combined[, is_overlap := check_aligns_overlap(QryStartPos, QryEndPos, QryContigID, RefcontigID), by=QryContigID]
  aligns_combined <- get_ref_align_xmap_filter(aligns_combined)

  print(paste0("Total in combined ", nrow(aligns_combined)))

  aligns_twocolor_filter <- aligns_twocolor[RefcontigID %in% aligns_combined[Colors==2, RefcontigID]]
  aligns_twocolor_filter <- flip.xmap(aligns_twocolor_filter)
  aligns_twocolor_filter$is_overlap <- NULL
  aligns_twocolor_filter$maxScore <- NULL
  aligns_twocolor_filter$Conf_delta <- NULL
  export_twocolor_aligns(aligns_twocolor_filter, two.color.file, two.color.file2, Confidence.T)

  return(aligns_twocolor_filter)
}

export_twocolor_aligns <- function(aligns_twocolor, aligns_twocolor1, aligns_twocolor2, Confidence.T){
  #export alignment
  prefix <- gsub("_hash1.xmap", "", aligns_twocolor1)
  ngs_maps_out <- readcmap_w_headers(gsub("\\.xmap", "_r.cmap", aligns_twocolor1))
  ngs_maps_out$cmap <- get.merged.cmap(
    list(gsub("\\.xmap", "_r.cmap", aligns_twocolor1),
         gsub("\\.xmap", "_r.cmap", aligns_twocolor2))
  )
  hybrid_maps_out <- readcmap_w_headers(gsub("\\.xmap", "_q.cmap", aligns_twocolor1))
  hybrid_maps_out$cmap <- get.merged.cmap(
    list(gsub("\\.xmap", "_q.cmap", aligns_twocolor1),
         gsub("\\.xmap", "_q.cmap", aligns_twocolor2))
  )
  exportCMap(ngs_maps_out$cmap, "", paste0(prefix,  '_q.cmap'), headers = ngs_maps_out$header, printColName = F)
  exportCMap(hybrid_maps_out$cmap, "", paste0(prefix, '_r.cmap'), headers = hybrid_maps_out$header, printColName = F)

  headers <- list(paste0("# XMAP File Version:\t0.2\n",
                         "# Reference Maps From:\t", prefix, '_r.cmap', "\n",
                         "# Query Maps From:\t", prefix, '_q.cmap'))

  exportXmap(aligns_twocolor, paste0(prefix, '_not_filtered.xmap'), headers = headers)
  aligns_twocolor <- aligns_twocolor[Confidence > Confidence.T]
  exportXmap(aligns_twocolor, paste0(prefix, '.xmap'), headers = headers, printColName = T)
}

#Align final consist of six alignment total
#two.color.file and file2 are two-color/enzyme alignment using different channel as hash to be comprehensive
#align1.file and alig2.file are single-color alignment to cover single-color regions in the hybrid-scaffolds
#finally single.align.file and single.align.file2 are alignment to the single-enzyme hybrid-scaffold so account for
#hybrids scaffolds that did not make it into the sandwich scaffolds (i.e. those that only has NGS aligned in one enzyme
#but not the others)
getAlignFinalStats <- function(align1.file, align2.file, two.color.file, two.color.file2, single.align1.file, single.align2.file,
                               Confidence.T=11, Confidence2.T=13, IDprefix1='100000', IDprefix2='200000'){
  xmap1 <- as.data.table(readxmap(align1.file))
  xmap2 <- as.data.table(readxmap(align2.file))

  empty.align <- as.data.table(matrix(data=0, nrow = 0, ncol = ncol(xmap1)))
  setnames(empty.align, colnames(xmap1))
  xmap.twocolor1 <- if(length(two.color.file) > 0){as.data.table(readxmap(two.color.file))}else{empty.align}
  xmap.twocolor2 <- if(length(two.color.file2) > 0){as.data.table(readxmap(two.color.file2))}else{empty.align}
  #xmap.twocolor <- merge.align(xmap.twocolor, xmap.twocolor2)
  xmap.twocolor <- xmap_filter_twocolor(as.data.table(rbind(readxmap(paste0(stripExtension(align1.file), '_1st_pass.xmap')),
                                                            readxmap(paste0(stripExtension(align1.file), '_2nd_pass.xmap')))),
                                        as.data.table(rbind(readxmap(paste0(stripExtension(align2.file), '_1st_pass.xmap')),
                                                            readxmap(paste0(stripExtension(align2.file), '_2nd_pass.xmap')))),
                                        xmap.twocolor1, xmap.twocolor2, two.color.file, two.color.file2, Confidence2.T)

  if(file.exists(single.align1.file) && file.exists(single.align2.file)){
    xmap.single1 <- as.data.table(readxmap(single.align1.file))
    xmap.single2 <- as.data.table(readxmap(single.align2.file))
  }else{
    xmap.single1 <- xmap1[grepl(IDprefix1, RefcontigID)]
    xmap.single2 <- xmap2[grepl(IDprefix2, RefcontigID)]
    xmap1 <- xmap1[!grepl(IDprefix1, RefcontigID)]
    xmap2 <- xmap2[!grepl(IDprefix2, RefcontigID)]
  }

  #twocolor.out <- gsub('hash1.xmap', "", two.color.file)
  #exportXmap(xmap.twocolor, paste0(twocolor.out, "_not_filterd.xmap"))

  xmap.twocolor <- remove.duplicated.align(xmap.twocolor, delta.score = 0.5, align.len = 0.0)

  xmap1 <- remove.duplicated.align(xmap1, delta.score = 0.5, align.len = 0.0)
  xmap2 <- remove.duplicated.align(xmap2, delta.score = 0.5, align.len =  0.0)

  xmap.twocolor <- xmap.twocolor[Confidence >= Confidence2.T]
  #exportXmap(xmap.twocolor, paste0(twocolor.out, ".xmap"))

  xmap1 <- xmap1[Confidence >= Confidence.T]
  xmap2 <- xmap2[Confidence >= Confidence.T]

  xmap.single1 <- xmap.single1[Confidence >= Confidence.T]
  xmap.single2 <- xmap.single2[Confidence >= Confidence.T]

  xmap.single.sandwich <- rbind(xmap1, xmap2[!(QryContigID %in% xmap1$QryContigID)])
  xmap.single.sandwich <- remove.duplicated.align(xmap.single.sandwich, delta.score = 0.0005, align.len = 0.0)

  xmap.single.hybrid <- rbind(xmap.single1, xmap.single2)
  xmap.single.hybrid <- remove.duplicated.align(xmap.single.hybrid, delta.score = 0.0005, align.len =  0.0)

  #browser()

  #alignment from sandwich scaffolds
  ids.two.sandwich <- unique(xmap.twocolor$QryContigID)
  ids1.sandwich <- unique(xmap1$QryContigID)
  ids2.sandwich <- unique(xmap2$QryContigID)
  #alignment from single-enzyme hybrids
  ids.single1 <- unique(xmap.single1$QryContigID)
  ids.single2 <- unique(xmap.single2$QryContigID)

  ids1.only.sandwich <- setdiff(ids1.sandwich, ids.two.sandwich)
  ids2.only.sandwich <- setdiff(ids2.sandwich, ids.two.sandwich)
  ids12.only.overlap <- intersect(ids1.only.sandwich, ids2.only.sandwich)


  ids.single.sandwich <- unique(c(ids1.sandwich, ids2.sandwich))
  ids.single.only.sandwich <- setdiff(ids.single.sandwich, ids.two.sandwich)
  ids.all.sandwich <- unique(c(xmap.twocolor$QryContigID, xmap1$QryContigID, xmap2$QryContigID))

  ids.single1.only <- setdiff(ids.single1, ids.all.sandwich)
  ids.single2.only <- setdiff(ids.single2, ids.all.sandwich)
  ids.single.only <- setdiff(c(ids.single1.only, ids.single2.only), ids.all.sandwich)
  ids.single12.overalp <- intersect(ids.single1.only, ids.single2.only)

  xmap.single.only.sandwich <- xmap.single.sandwich[QryContigID %in% ids.single.only.sandwich]
  xmap.single.only <- xmap.single.hybrid[QryContigID %in% ids.single.only]

  xmap.export <- getExportTwoColorAlign(xmap.twocolor)

  #browser()
  #xmap.export <- rbind(xmap.export, xmap.single.sandwich)
  #xmap.export <- remove.duplicated.align(xmap.export, delta.score = 0.01, align.len = 0.0, bestRef = T)
  xmap.export <- rbind(xmap.export, xmap.single.only.sandwich)
  xmap.export <- rbind(xmap.export, xmap.single.only)

  #removing embedded seq contigs from statistics
  embed.ids <- get.embed.aligns(xmap.export)
  xmap.embed <- xmap.export[QryContigID %in% embed.ids]

  ids.two.sandwich <- ids.two.sandwich[!(ids.two.sandwich %in% embed.ids)]
  ids1.sandwich <- ids1.sandwich[!(ids1.sandwich %in% embed.ids)]
  ids2.sandwich <- ids2.sandwich[!(ids2.sandwich %in% embed.ids)]
  ids12.only.overlap <- ids12.only.overlap[!(ids12.only.overlap %in% embed.ids)]
  ids1.only.sandwich <- ids1.only.sandwich[!(ids1.only.sandwich %in% embed.ids)]
  ids2.only.sandwich <- ids2.only.sandwich[!(ids2.only.sandwich %in% embed.ids)]
  ids.single.only.sandwich <- ids.single.only.sandwich[!(ids.single.only.sandwich %in% embed.ids)]
  ids.all.sandwich <- ids.all.sandwich[!(ids.all.sandwich %in% embed.ids)]
  ids.single1 <- ids.single1[!(ids.single1 %in% embed.ids)]
  ids.single2 <- ids.single2[!(ids.single2 %in% embed.ids)]
  ids.single1.only <- ids.single1.only[!(ids.single1.only %in% embed.ids)]
  ids.single2.only <- ids.single2.only[!(ids.single2.only %in% embed.ids)]
  ids.single12.overalp <- ids.single12.overalp[!(ids.single12.overalp %in% embed.ids)]
  ids.single.only <- ids.single.only[!(ids.single.only %in% embed.ids)]
  ids.all.sandwich <- ids.all.sandwich[!(ids.all.sandwich %in% embed.ids)]


  xmap.export <- xmap.export[!(QryContigID %in% embed.ids)]
  xmap.twocolor <- xmap.twocolor[!(QryContigID %in% embed.ids)]
  xmap1 <- xmap1[!(QryContigID %in% embed.ids)]
  xmap2<- xmap2[!(QryContigID %in% embed.ids)]
  xmap.single1 <- xmap.single1[!(QryContigID %in% embed.ids)]
  xmap.single2 <- xmap.single2[!(QryContigID %in% embed.ids)]
  xmap.single.only <- xmap.single.only[!(QryContigID %in% embed.ids)]
  xmap.single.only.sandwich <- xmap.single.only.sandwich[!(QryContigID %in% embed.ids)]


  #analyzing the alignment difference
  sandwich.ngs.align <- data.frame(
               sandwich_twoEnzy=length(ids.two.sandwich),
               sandwich_enzy1=length(ids1.sandwich), sandwich_enzy2=length(ids2.sandwich),
               sandwich_enzy1_only=length(ids1.only.sandwich), sandwich_enzy2_only=length(ids2.only.sandwich),
               sandwich_enzy1_enzy2_overlap=length(ids12.only.overlap),
               sandwich_total=(length(ids.single.only.sandwich) + length(ids.two.sandwich)))
  sandwich.ngs.Len <- data.frame(
              sandwich_twoEnzy=sum(xmap.twocolor[,QryLen/2]),
              sandwich_enzy1=sum(xmap1[,QryLen]), sandwich_enzy2=sum(xmap2[,QryLen]),
              sandwich_enzy1_only=sum(xmap1[QryContigID %in% ids1.only.sandwich, QryLen]),
              sandwich_enzy2_only=sum(xmap2[QryContigID %in% ids2.only.sandwich, QryLen]),
              sandwich_total=(sum(xmap.twocolor[,QryLen/2]) + sum(xmap.single.only.sandwich[,QryLen])))

  hybrid.ngs.align <- data.frame(
              hybrid_twoEnzy=length(ids.all.sandwich),
              hybrid_enzy1=length(ids.single1), hybrid_enzy2=length(ids.single2),
              hybrid_enzy1_only=length(ids.single1.only), hybrid_enzy2_only=length(ids.single2.only),
              hybrid_enzy1_enzy2_overlap=length(ids.single12.overalp),
              hybrid_total=(length(ids.single.only) + length(ids.all.sandwich)),
              embedded = length(embed.ids))

  hybrid.ngs.Len <- data.frame(
              hybrid_twoEnzy=(sum(xmap.twocolor[,QryLen/2]) + sum(xmap.single.only.sandwich[,QryLen])),
              hybrid_enzy1=sum(xmap.single1[,QryLen]), hybrid_enzy2=sum(xmap.single2[,QryLen]),
              hybrid_enzy1_only=sum(xmap.single1[QryContigID %in% ids.single1.only, QryLen]),
              hybrid_enzy2_only=sum(xmap.single2[QryContigID %in% ids.single2.only, QryLen]),
              hybrid_total=(sum(xmap.twocolor[,QryLen/2]) + sum(xmap.single.only.sandwich[,QryLen]) + sum(xmap.single.only[,QryLen])),
              embedded = xmap.embed[,sum(QryLen)])

  summary.stats <- list(sandwich_ngs_align=sandwich.ngs.align, sandiwch_ngs_len=sandwich.ngs.Len,
                  hybrid_ngs_align = hybrid.ngs.align, hybrid_ngs_Len = hybrid.ngs.Len)

  print(summary.stats)
  #return(list(ExportAlign=xmap.export, TwoColorAlign=xmap.twocolor,
  #            SingleOnlyAlign=xmap.single.only.sandwich, SingleHybridAlign=xmap.single.only,
  #            SummaryStats=summary.stats))
  return(list(ExportAlign=xmap.export, TwoColorAlign=xmap.twocolor,
              SingleOnlyAlign=xmap.single.only.sandwich, SingleHybridAlign=xmap.single.only,
              SummaryStats=summary.stats))
}

testGetAlignFinalStat <- function(){
  two.color.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/RunWithT11Len80/alignfinal/mar3_NA12878.scf.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash1_errEst.xmap"
  two.color.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/RunWithT11Len80/alignfinal/mar3_NA12878.scf.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash2_errEst.xmap"
  align1.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/RunWithT11Len80/alignfinal/mar3_NA12878.scf.fasta.cutted_BSPQI_BSSSI_0kb_0labels_1_errEst.xmap"
  align2.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/RunWithT11Len80/alignfinal/mar3_NA12878.scf.fasta.cutted_BSPQI_BSSSI_0kb_0labels_2_errEst.xmap"

  two.color.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoAndyDebug/alignfinal_old_params/chm13.quivered.fa.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash1_errEst.xmap"
  two.color.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoAndyDebug/alignfinal_old_params/chm13.quivered.fa.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash2_errEst.xmap"
  align1.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoAndyDebug/alignfinal_manual/chm13.quivered.fa.cutted_BSPQI_0kb_0labels_1_errEst.xmap"
  align2.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoAndyDebug/alignfinal_manual/chm13.quivered.fa.cutted_BSSSI_0kb_0labels_1_errEst.xmap"
  align.single1.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunOak/alignfinal_old/AlignToSingleHybrid/Qlob.v1.0.fasta.cutted_BSPQI_BSSSI_0kb_0labels_1_errEst.xmap"
  align.single2.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunOak/alignfinal_old/AlignToSingleHybrid/Qlob.v1.0.fasta.cutted_BSPQI_BSSSI_0kb_0labels_2_errEst.xmap"

  aligns <- getAlignFinalStats(align1.file, align2.file, two.color.file, two.color.file2, "", "",  9, 12)
  export.align.path <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoAndyDebug/alignfinal_manual/chm13.quivered.fa.cutted_BSPQI_BSSSI_0kb_0labels_1_aligns.RData"
  saveRDS(aligns, export.align.path)
  export.align.path <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoAndyDebug/alignfinal_manual/chm13.quivered.fa.cutted_BSPQI_BSSSI_0kb_0labels_1_Export.xmap"
  exportXmap(aligns$ExportAlign, export.align.path)

  #test with new greedy-alg of align final
  two.color.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TestRunWithTwoPass/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_E_CTTAAG_Q_NGS_A_HYBRID_hash1.xmap"
  two.color.file2 <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TestRunWithTwoPass/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_E_CTTAAG_Q_NGS_A_HYBRID_hash2.xmap"
  align1.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TestRunWithTwoPass/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID.xmap"
  align2.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TestRunWithTwoPass/two_enzyme_hybrid_scaffold_M1/alignfinal/E_CTTAAG_Q_NGScontigs_A_HYBRID.xmap"
  aligns <- getAlignFinalStats(align1.file, align2.file, two.color.file, two.color.file2, "", "", 9, 12)

  two.color.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/TestRun/SBJ/PB40x/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_E_BSSSI_Q_NGS_A_HYBRID_hash1.xmap"
  two.color.file2 <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/TestRun/SBJ/PB40x/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_E_BSSSI_Q_NGS_A_HYBRID_hash2.xmap"
  align1.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/TestRun/SBJ/PB40x/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID.xmap"
  align2.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/TestRun/SBJ/PB40x/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSSSI_Q_NGScontigs_A_HYBRID.xmap"
  aligns <- getAlignFinalStats(align1.file, align2.file, two.color.file, two.color.file2, "", "", 9, 12)

  two.color.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/E_BSPQI_E_CTTAAG_Q_NGS_A_HYBRID_test_multimatch_hash1.xmap"
  two.color.file2 <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/E_BSPQI_E_CTTAAG_Q_NGS_A_HYBRID_test_multimatches_hash2.xmap"
  align1.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID_multimatches.xmap"
  align2.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/E_CTTAAG_Q_NGScontigs_A_HYBRID_multimatches.xmap"

  process_align_final(align1.file)
  process_align_final(align2.file)
  aligns <- getAlignFinalStats(align1.file, align2.file, two.color.file, two.color.file2, "", "", 9, 12)

  two.color.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_E_BSSSI_Q_NGS_A_HYBRID_hash1.xmap"
  two.color.file2 <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_E_BSSSI_Q_NGS_A_HYBRID_hash2.xmap"
  align1.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID.xmap"
  align2.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSSSI_Q_NGScontigs_A_HYBRID.xmap"

  aligns <- getAlignFinalStats(align1.file, align2.file, two.color.file, two.color.file2, "", "", 10, 12)


  two.color.file <- numeric(0)
  two.color.file2 <- numeric(0)
  align1.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRun_MaskInd_Debug_UNC_Pfal/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSMI_Q_NGScontigs_A_HYBRID.xmap"
  align2.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRun_MaskInd_Debug_UNC_Pfal/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSMI_Q_NGScontigs_A_HYBRID.xmap"

  aligns <- getAlignFinalStats(align1.file, align2.file, two.color.file, two.color.file2, "", "", 10, 12)

}

testGetNGSAlignFract <- function(){
  ngs.hu.sop <- as.data.table(readcmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestAutoSOP/BSPQI/fa2cmap/a.lines_BSPQI_0kb_0labels.cmap'))
  align.hu.sop <- as.data.table(readxmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestAutoSOP/BSPQI/hybrid_scaffolds_M1/EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap'))
  ngs.oak <- as.data.table((readcmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunOak/BSPQI/fa2cmap/Qlob.v1.0_BSPQI_0kb_0labels.cmap")))
  align.oak <- as.data.table((readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunOak/BSPQI/hybrid_scaffolds_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_Qlob_v1_0_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap")))
  ngs.hu.pb <- as.data.table(readcmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/BSPQI/assignAlignType/cut_conflicts_M1/ngs_cut_pre_exclude.cmap'))
  align.hu.pb <- as.data.table(readxmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/BSPQI/align_final_M1/EXP_REFINEFINAL1_bppAdjust_cmap_mar3_NA12878_scf_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap'))
  ngs.hu.soap <- as.data.table(readcmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/BSPQI/assignAlignType/cut_conflicts_M1/ngs_cut_pre_exclude.cmap'))
  ngs.hu.soap <- ngs.hu.soap[ContigLength > 5000]
  align.hu.soap <- as.data.table(readxmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/BSPQI/hybrid_scaffolds_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_soap_result_scafSeq_fa_NGScontigs_HYBRID_SCAFFOLD.xmap'))
  ngs.hu.pb2 <- as.data.table(readcmap('~/RemoteServer/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_pipeline_020215/v1/output_results/fa2cmap/consensus_diploid_for_andy_forncbi_filtered2b_BspQI_0Kb_0labels.cmap'))
  align.hu.pb2 <- as.data.table(readxmap('~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_pipeline_020215/v1/output_results/align_final/EXP_REFINEFINAL1_bppAdjust_cmap_consensus_diploid_for_andy_forncbi_filtered2b_fa_NGScontigs_HYBRID_SCAFFOLD.xmap'))

  g <- getAlignmentRate(align.hu.sop, ngs.hu.sop, label='Hu-Dis')
  g <- getAlignmentRate(align.oak, ngs.oak, label='Oak', g1 = g[[1]], g2 = g[[2]])
  g <- getAlignmentRate(align.hu.pb, ngs.hu.pb, label='Hu-PacBio', g1 = g[[1]], g2=g[[2]], numlabels.bins = seq(0,1000, 5))
  g <- getAlignmentRate(align.hu.soap, ngs.hu.soap, label='Hu-SOAPScaf', g1=g[[1]], g2=g[[2]], numlabels.bins = seq(0,1000, 5))
  g <- getAlignmentRate(align.hu.pb2, ngs.hu.pb2, label='PacBio-old-run', g1=g[[1]], g2=g[[2]], numlabels.bins=seq(0,1000,5), xlim1=c(4.5,6))
}

#this function compute embedded alignments from
#an xmap file, it assumes each query alignment is has
#one entry in the xmap, if this is not true, it will pick
#the best scoring one
get.embed.aligns <- function(align){
  if(nrow(align) < 2){
    return(list())
  }
  align <- remove.duplicated.align(align, delta.score = 0.0001, align.len = 0.0, bestRef = T)
  align <- as.data.table(align)
  align <- align[order(RefcontigID, RefStartPos, -RefEndPos)]
  align <- getProjAlignPos(align)
  align[,GapLength := c(0, but.first(RefStartPos) - but.last(RefEndPos) - 1), by=RefcontigID]
  align[,AdjGapLength := c(0, but.first(ProjRefStart) - but.last(ProjRefEnd) - 1), by=RefcontigID]
  align[,IsEmbedded := -1*AdjGapLength > QryLen]
  prev_i <- 1
  gap.lens <- align$GapLength
  adj.gap.lens <- align$AdjGapLength
  embed <- align$IsEmbedded
  for(i in 3:nrow(align)){
      if(embed[i-1] && align$RefcontigID[i] == align$RefcontigID[i-1]){
        gap.lens[i] <- align$RefStartPos[i] - align$RefEndPos[prev_i] - 1
        adj.gap.lens[i] <- align$ProjRefStart[i] - align$ProjRefEnd[prev_i] - 1
        if(-1*adj.gap.lens[i] > align$QryLen[i]){
           embed[i] <- TRUE
        }
      }else{
          prev_i <- i - 1
      }
  }
  align$GapLength <- gap.lens
  align$AdjGapLength <- adj.gap.lens
  align$IsEmbedded <- embed
  return(align[IsEmbedded==TRUE, QryContigID])
}


#analyze the two-color and single-color scaffolding for hb
analyzeAlignFinalHB <- function(){
  two.color.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/alignfinal/p_consensus.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash1_errEst.xmap"
  two.color.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/alignfinal/p_consensus.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash2_errEst.xmap"

  align1.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/alignfinal/p_consensus.fasta.cutted_BSPQI_BSSSI_0kb_0labels_1_errEst.xmap"
  align2.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/alignfinal/p_consensus.fasta.cutted_BSPQI_BSSSI_0kb_0labels_2_errEst.xmap"

  single.align1.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/RunWithT11Len80/alignfinal/SingleColorReRun/p_consensus.fasta.cutted_BSPQI_BSSSI_0kb_0labels_1_errEst.xmap"
  single.align2.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/RunWithT11Len80/alignfinal/SingleColorReRun/p_consensus.fasta.cutted_BSPQI_BSSSI_0kb_0labels_2_errEst.xmap"


  aligned.final.ids <- getAlignFinalStats(align1.file, align2.file, two.color.file, two.color.file2, 11, 13)

  xmap1 <- as.data.table(readxmap(single.align1.file))
  xmap2 <- as.data.table(readxmap(single.align2.file))
  xmap1 <- xmap1[Confidence > 11]
  xmap2 <- xmap2[Confidence > 11]

  xmap.intersect <- merge(xmap1, xmap2, by='QryContigID')

  qcmap <- as.data.table(readcmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/RunWithT11Len80/fa2cmap2color/p_consensus.fasta.cutted_BSPQI_BSSSI_0kb_0labels.cmap"))

  bspqi.only <- setdiff(unique(xmap1$QryContigID), unique(xmap2$QryContigID))

  setkey(qcmap, 'CMapId')
  qstats.bspqi.only <- qcmap[bspqi.only][, list(Len=ContigLength[1], Site1=sum(LabelChannel==1), Site2=sum(LabelChannel==2)), by=CMapId]

}


#this function is for merge all single-color and two-color alignments to a single alignfinal
mergeAlignFinal <- function(){
  setkey(xmap.twocolor, 'QryContigID')
  setkey(xmap1, 'QryContigID')
  setkey(xmap2, 'QryContigID')

  intersect1 <- merge(xmap.twocolor, xmap1, by='QryContigID')
  intersect2 <- merge(xmap.twocolor, xmap2, by='QryContigID')
  intersect3 <- merge(xmap1, xmap2, by='QryContigID')

  mid.dif1 <- intersect1[seq(1,nrow(intersect1),2)][RefcontigID.x == RefcontigID.y,abs((RefStartPos.x + RefEndPos.x)/2-(RefStartPos.y+RefEndPos.y)/2)]
  mid.dif2 <- intersect2[seq(2,nrow(intersect1),2)][RefcontigID.x == RefcontigID.y,abs((RefStartPos.x + RefEndPos.x)/2-(RefStartPos.y+RefEndPos.y)/2)]
  mid.dif3 <- intersect3[RefcontigID.x == RefcontigID.y,(RefStartPos.x + RefEndPos.x)/2-(RefStartPos.y+RefEndPos.y)/2]

  d <- data.frame(Label=c(rep("Two-Enz vs Enz1", length(mid.dif1)), rep("Two-Enzy vs Enz2", length(mid.dif2)), rep("Enzy1 vs Enzy2", length(mid.dif3))),
                MidPointDifference = c(mid.dif1, mid.dif2, mid.dif3))
  g <- ggplot(d) + geom_freqpoly(aes(x=MidPointDifference, color = Label, stat_bin=50)) + scale_x_log10(breaks=10^c(0,1,2,3,4,5,6,7))
  plot(g)

  #merge aligned ngs from all three alignment
  xmap1.only <- xmap1[QryContigID %in% ids1.only]
  xmap2.only <- xmap2[QryContigID %in% ids2.only]

  overlap.align1.only <- lapply(1:nrow(xmap1.only), function(i){
                                  align.overlap <- xmap.twocolor[RefcontigID == xmap1.only$RefcontigID[i]
                                                               & (RefStartPos >= xmap1.only$RefStartPos[i] & RefStartPos <= xmap1.only$RefEndPos[i]
                                                                 | RefEndPos >= xmap1.only$RefStartPos[i] & RefEndPos <= xmap1.only$RefEndPos[i])]
                              })
  print(paste0("Length of NGS aligned: 2-color: ", sum(xmap.twocolor[,QryLen/2]), " enzyme1: ", sum(xmap1[,QryLen]), " enzyme2: ", sum(xmap2[,QryLen]),
               " enzyme1-only: ", sum(xmap1[QryContigID %in% ids1.only, QryLen]), " enzyme2-only: ", sum(xmap2[QryContigID %in% ids2.only, QryLen]),
               " total: ", (sum(xmap.twocolor[,QryLen/2]) + sum(xmap1[QryContigID %in% ids1.only, QryLen]) + sum(xmap2[QryContigID %in% ids2.only, QryLen]))))
  overlap.align2.only <- lapply(1:nrow(xmap2.only), function(i){
    #browser()
    align.overlap <- xmap.twocolor[RefcontigID == xmap2.only$RefcontigID[i]
                                   & ((RefStartPos >= xmap2.only$RefStartPos[i] & RefStartPos <= xmap2.only$RefEndPos[i])
                                      | (RefEndPos >= xmap2.only$RefStartPos[i] & RefEndPos <= xmap2.only$RefEndPos[i]))]
    return(align.overlap)
  })
}



#merge two two-color alignments (presumably each use a different hashcolor)
#to generate a merged two-color alignment
#use.color allow only outputting the alignment from one channel
get.merged.twocolor.align <- function(align.file1, align.file2, use.color=3){
  align1 <- as.data.table(readxmap(align.file1))
  align2 <- as.data.table(readxmap(align.file2))
  align <- merge.align(align1, align2)
  if(use.color==1){
    align <- align[c(TRUE,FALSE)]
  }else if(use.color==2){
    align <- align[c(FALSE, TRUE)]
  }
  return(align)
}


#get consecutive ID pairs of a ordered set of ids
get.pairs <- function(ID.list, dist=1){
  true.pairs <- paste0(ID.list, '_', c(ID.list[-(1:dist)], rep(0, dist)))
  rev.true.pairs <- paste0(rev(ID.list), '_', c(rev(ID.list)[-(1:dist)], rep(0, dist)))
  return(c(true.pairs, rev.true.pairs))
}

#based on a two-color alignment, calculate a single-coordinate the represent the location of the
#qry contig on the reference
getContigCoordinate <- function(xmap.two.color){
  xmap.two.color <- xmap.two.color[,Location := mean(0.5*(RefStartPos + RefEndPos)), by=QryContigID]
  return(xmap.two.color)
}

#getting the non-scaffold contigs in sandwich assembly
test.get.nonscaffold.contigs <- function(){

}


#check accuracy of NGS anchoring using two-color
checkNGSAlign <- function(){
  ref1.file <- '~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestTwoColor/OneStep_alignFinalParams/Twocolor/Hash1/a.lines_BSPQI_BSSSI_0kb_0labels_twocolor_errEst.xmap'
  ref2.file <- '~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestTwoColor/OneStep_alignFinalParams/Twocolor/Hash2/a.lines_BSPQI_BSSSI_0kb_0labels_twocolor_errEst.xmap'
  ref.bspqi.file <- '~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestTwoColor/OneStep_alignFinalParams/BSPQI/a.lines_BSPQI_0kb_0labels_1_errEst.xmap'
  ref.bsssi.file <- '~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestTwoColor/OneStep_alignFinalParams/BSSSI/a.lines_BSSSI_0kb_0labels_1_errEst.xmap'
  ref <- get.merged.twocolor.align(ref1.file, ref2.file)
  ref <- ref[Confidence > 12]

  ref.bspqi <- as.data.table(readxmap(ref.bspqi.file))
  ref.bspqi <- ref.bspqi[Confidence > 11]
  ref.bsssi <- as.data.table(readxmap(ref.bsssi.file))
  ref.bsssi <- ref.bsssi[Confidence > 11]

  #analyze consistency of answer itself
  setkey(ref, 'QryContigID')
  setkey(ref.bspqi, 'QryContigID')
  setkey(ref.bsssi, 'QryContigID')

  ref.merged <- merge(ref.bspqi, ref.bsssi, by = 'QryContigID', all.x = T, all.y= T)
  ref.merged <- merge(ref.merged, ref, by='QryContigID', all.x = T, all.y = T)
  ref.merged[is.na(Confidence), Confidence := -1]
  ref.merged[is.na(Confidence.x), Confidence.x := -1]
  ref.merged[is.na(Confidence.y), Confidence.y := -1]

  ref.merged[,sum(Confidence > Confidence.x | Confidence.y > Confidence)]

  #some sanity check
  stat <- quote(list(Confidence.x, Confidence.y, Confidence, QryContigID, RefcontigID, RefcontigID.x, RefcontigID.y, RefStartPos, RefEndPos))
  #check any two-color alignment in conflict with a consensus single-color alignment
  nrow(ref.merged[RefcontigID != RefcontigID.x & RefcontigID.x == RefcontigID.y & !is.na(RefcontigID.x)])

  bspqi.only.count <- nrow(ref.merged[is.na(RefcontigID) & !is.na(RefcontigID.x)])

  bsssi.only.count <- nrow(ref.merged[is.na(RefcontigID) & !is.na(RefcontigID.y)])

  #check NGS anchoring order accuracy from scaffolding methods
  two.color.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies_withM2CutContigs/NGSAlignFinal/Hash1/a.lines_BSPQI_BSSSI_0kb_0labels_twocolor_errEst.xmap"
  two.color.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies_withM2CutContigs/NGSAlignFinal/Hash2/a.lines_BSPQI_BSSSI_0kb_0labels_twocolor_errEst.xmap"

  two.color.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestTwoColor/OneStep_alignFinalParams/Twocolor_withupdateparams/a.lines_BSPQI_BSSSI_0kb_0labels_twocolor_hash1_errEst.xmap"
  two.color.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestTwoColor/OneStep_alignFinalParams/Twocolor_withupdateparams/a.lines_BSPQI_BSSSI_0kb_0labels_twocolor_hash2_errEst.xmap"
  align1.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies_withHybridM3/AlignFinal/BSPQI/a.lines_BSPQI_0kb_0labels_1_errEst.xmap"
  align2.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies_withHybridM3/AlignFinal/BSSSI/a.lines_BSSSI_0kb_0labels_1_errEst.xmap"


  align.final.twocolor <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestTw"))
  align.final.twocolor <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/"))

  align.final.twocolor <- get.merged.twocolor.align(two.color.file, two.color.file2, use.color = 1)
  align.final.twocolor <- align.final.twocolor[Confidence > 12]
  align.final.twocolor <- align.final.twocolor[,Confidence := 100*Confidence] #taking precedent over twocolor alignment

  align.bspqi <- as.data.table(readxmap(align1.file))
  align.bsssi <- as.data.table(readxmap(align2.file))

  align.bspqi <- align.bspqi[Confidence > 11]
  align.bsssi <- align.bsssi[Confidence > 11]

  align.final.twocolor <- merge.align(align.final.twocolor, align.bspqi)
  align.final.twocolor <- merge.align(align.final.twocolor, align.bsssi)

  #align.final.twocolor <- align.final.twocolor[abs(QryStartPos-QryEndPos)/QryLen > 0.9]
  align.final.twocolor <- align.final.twocolor[Confidence > 11]
  align.final.twocolor <- align.final.twocolor[(QryContigID %in% ref$QryContigID)==TRUE]
  align.final.twocolor <- getContigCoordinate(align.final.twocolor)
  setkeyv(align.final.twocolor, c("RefcontigID", "Location"))
  align.final.ids <- align.final.twocolor[, list(RefcontigID,QryContigID)]  #select every other row for qrycontigID
  scaffolded.ngs.pairs <- align.final.ids[,list(QryContigID[-length(QryContigID)],  QryContigID[-1]), by=RefcontigID]


  ref.select <- ref[(QryContigID %in% align.final.twocolor$QryContigID)==T]
  ref.select <- getContigCoordinate(ref.select)
  ref.select <- ref.select[c(TRUE, FALSE)==TRUE]

  setkeyv(ref.select, c("RefcontigID", "Location"))

  ref.order <- ref.select[,list(QryContigID)]
  ref.order <- ref.order[,Indice:=1:nrow(ref.order)]
  setkey(ref.order, QryContigID)

  order.dist <- scaffolded.ngs.pairs[, dist:=ref.order[V2, Indice] - ref.order[V1, Indice]]
  nrow(order.dist) - order.dist[,sum(abs(dist)==1)]
}

mergeByTwoColorAlign <- function(){
  two.color1.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestAutoSOP/alignfinal/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash1_errEst.xmap"
  two.color2.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestAutoSOP/alignfinal/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash2_errEst.xmap"
  two.color.align1 <- as.data.table(readxmap(two.color1.file))
  two.color.align2 <- as.data.table(readxmap(two.color2.file))
  two.color.align <- merge.align(two.color.align1, two.color.align2)
  two.color.align1 <- subset(two.color.align1, Confidence > 11)
  mi <- generateMappingInfo(two.color.align1[seq(1,nrow(two.color.align1),2),])

}

#this function implements the greedy algorithm for incorporation of NGS alignment
#to hybrid-scaffold
get_ref_align_xmap_filter <- function(aligns, T0=13, T1 = 10, T2 = 9, T_delta = 2, cleanUp=TRUE){
  #aligns[, maxScore:=max(Confidence), by=QryContigID]
  #aligns <- aligns[order(maxScore, QryContigID, Confidence, decreasing = T)]
  aligns <- aligns[order(Confidence, decreasing = T)]
  aligns[, is_overlap := check_aligns_overlap(QryStartPos, QryEndPos), by=QryContigID]
  aligns_filter1 <- aligns[is_overlap==FALSE | Confidence > T0]
  aligns_filter1[, Conf_delta := get_delta_score(Confidence), by=RefcontigID]
  aligns_filter2 <- aligns_filter1[Confidence > T1 | (Confidence > T2 & Conf_delta > T_delta)]
  aligns_filter2[, maxScore := max(Confidence), by=RefcontigID]
  aligns_filter3 <- aligns_filter2[Confidence==maxScore]
  #we produce auxilliary columns to perform filtering, we remove them after filtering
  #they can be kept for further processing
  if(cleanUp){
    aligns_filter3$is_overlap <- NULL
    aligns_filter3$Conf_delta <- NULL
    aligns_filter3$maxScore <- NULL
  }
  return(aligns_filter3)
}

get_delta_score <- function(scores){
  if(length(scores)==1){
    return(1000)
  }else{
    c(-1*(diff(scores)), 0)
  }
}

check_aligns_overlap <- function(startPos, endPos, IDs=1, QIDs=1, buffer=20e3){
  ind <- startPos > endPos
  temp <- startPos[ind]
  startPos[ind] <- endPos[ind]
  endPos[ind] <- temp

  qryStart <- startPos + buffer
  qryEnd <- endPos - buffer
  qryEnd <- ifelse(qryEnd <= qryStart, qryStart+1, qryEnd)

  overlap <- intervals_overlap_wrapper(data.frame(qryStart, qryEnd),
                                       data.frame(startPos, endPos))

  is_overlap <- unlist(lapply(1:length(overlap),
                              function(i){
                                currInd <- as.integer(ls(overlap[i]))
                                any(overlap[[i]] < currInd)
                              }))
  return(is_overlap)
}

#this filter the align_final alignments to determine the best alignment for each NGS contigs
#and output the alignment using the hybrid as the reference and NGS contigs as query for visualization
process_align_final <- function(align.final.file, aligns=NULL){
  if(is.null(aligns)){
    aligns <- as.data.table(readxmap(align.final.file))
  }
  aligns_out <- get_ref_align_xmap_filter(aligns)
  exportXmap(aligns, paste0(align.final.file, '.tmp'))
  aligns_out <- flip.xmap(as.data.frame(aligns_out))
  aligns_out <- aligns_out[order(aligns_out$XmapEntryID),]
  exportXmap(aligns_out, align.final.file)
}

#this output a list of unaligned contigs from a previous-run of alignment
get.unalign.contigs.single <- function(align.file, cmap.file, min.conf=13, min.lab=5, Id.col){
  align <- as.data.table(readxmap(align.file))
  cmap <- as.data.table(readcmap(cmap.file))
  if(nrow(align) == 0){
    warning(paste0("There is no alignment in file ", align.file))
    return(cmap)
  }
  if(Id.col == 'Qry'){
    align.ids <-unique(c(align[Confidence > min.conf, QryContigID]))
  }else if(Id.col== 'Ref'){
    align.ids <- unique(c(align[Confidence > min.conf, RefcontigID]))
  }else{
    stop("Unknow ID columsn, please specify either Qry or Ref")
  }
  print(paste0("Number of NGS contig anchored: ", length(align.ids)))
  unaligned_cmap <- cmap[!(CMapId %in% align.ids) & NumSites >= min.lab]
  print(paste0("Total unalign contigs with minimum labels ", nrow(unaligned_cmap[SiteID==1])))
  return(unaligned_cmap)
}

#this function perform a two-pass alignment using one set of parameter to align the large contigs and another
#set of parameters to align the short contigs, this ensure we have both high-sensitivity and speed when
#performing alignfinal
run_alignfinal_two_pass <- function(hybrid, ngs, prefix, run_flag=F, ra=NULL, bestRef=1, filter=F, args1=NULL, args2=NULL){
  #running first pass
  first_pass_out <- paste0(prefix, '_1st_pass')
  r <- run_ra_alignref_final(hybrid, ngs, first_pass_out, run_flag, ra, bestRef, filter=F, args1)
  unaligned_ngs <- get.unalign.contigs.single(paste0(first_pass_out, '.xmap'),  ngs, min.conf=1, min.lab=3, Id.col='Ref') #note here we use NGS as reference
  unaligned_ngs_file <- paste0(prefix, '_1st_pass_unaligned.cmap')
  exportCMap(unaligned_ngs, "", unaligned_ngs_file)

  #running second pass
  second_pass_out <- paste0(prefix, '_2nd_pass')
  if(nrow(unaligned_ngs) > 0){
    r <- run_ra_alignref_final(hybrid, unaligned_ngs_file, second_pass_out, run_flag, ra, bestRef, filter=F, args2)
  }else{
    #if there is no sequence left we output empty result file in second pass
    align1 <- as.data.table(readxmap(paste0(first_pass_out, '.xmap')))
    empty <- align1[QryContigID != QryContigID]
    exportXmap(empty, paste0(second_pass_out, '.xmap'))

    cmap <- as.data.table(readcmap(paste0(first_pass_out, '_r.cmap')))
    empty <- cmap[CMapId != CMapId]
    exportCMap(empty, "", paste0(second_pass_out, '_r.cmap'))

    cmap <- as.data.table(readcmap(paste0(first_pass_out, '_q.cmap')))
    empty <- cmap[CMapId != CMapId]
    exportCMap(empty, "", paste0(second_pass_out, '_q.cmap'))
  }

  #combining alignments
  align1 <- as.data.table(readxmap(paste0(first_pass_out, '.xmap')))
  align2 <- as.data.table(readxmap(paste0(second_pass_out, '.xmap')))
  aligns <- rbind(align1, align2)
  aligns_out_all <- flip.xmap(aligns)
  #filtering alignment
  aligns_filter <- get_ref_align_xmap_filter(aligns)

  #combining aligned NGS and hybrid cmap
  ngs_maps <- readcmap_w_headers(ngs)
  #ngs_maps$cmap <- as.data.table(ngs_maps$cmap)
  #ngs_maps_out <- ngs_maps$cmap[CMapId %in% aligns_filter$RefcontigID]
  maps <- get.merged.cmap(list(
    paste0(first_pass_out, '_r.cmap'),
    paste0(second_pass_out, '_r.cmap')
  ))
  ngs_maps_out <- ngs_maps
  ngs_maps_out$cmap <- as.data.frame(maps)

  hybrid_maps <- readcmap_w_headers(hybrid)
  #hybrid_maps$cmap <- as.data.table(hybrid_maps$cmap)
  #hybrid_maps_out <- hybrid_maps$cmap[CMapId %in% aligns_filter$QryContigID]
  maps <- get.merged.cmap(list(
    paste0(first_pass_out, '_q.cmap'),
    paste0(second_pass_out, '_q.cmap')
  ))
  hybrid_maps_out <- hybrid_maps
  hybrid_maps_out$cmap <- as.data.frame(maps)

  #flip align columns and output alignment
  aligns_out <- flip.xmap(aligns_filter)
  exportCMap(ngs_maps_out$cmap, "", paste0(prefix,  '_q.cmap'), headers = ngs_maps$header, printColName = F)
  exportCMap(hybrid_maps_out$cmap, "", paste0(prefix, '_r.cmap'), headers = hybrid_maps$header, printColName = F)

  headers <- list(paste0("# XMAP File Version:\t0.2\n",
                         "# Reference Maps From:\t", prefix, '_r.cmap', "\n",
                         "# Query Maps From:\t", prefix, '_q.cmap'))

  exportXmap(aligns_out_all, paste0(prefix, '_not_filtered.xmap'), headers = headers)
  exportXmap(aligns_out, paste0(prefix, '.xmap'), headers = headers, printColName = T)
}

test_two_pass_align_final <- function(){
  tgh.xml <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/TGH_two_passes.xml"
  tgh.xml <- "~/RemoteServer//home/users/jwang/workspace/codebase/HybridScaffoldThree_trunk/trunk/TGH/hybridScaffold_two_enzymes_DLE1.xml"
  params <- xmlParse(tgh.xml)
  params <- xmlRoot(params)

  align_final1 <- list(params[['hybridScaffold1']][['align_final_1st_pass']],
                       params[['hybridScaffold1']][['global']])
  align_final2 <- list(params[['hybridScaffold1']][['align_final_2nd_pass']],
                       params[['hybridScaffold1']][['global']] )

  hybrid <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/Joyce_Alignment/two_enzyme_hybrid_scaffold_M1/Sandwich2/test_2col_c2.cmap"
  ngs <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/Joyce_Alignment/two_enzyme_hybrid_scaffold_M1/fa2cmap/p_and_a_ctg.fasta.cut_CTTAAG_0kb_0labels.cmap"
  prefix <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/Joyce_Alignment/two_enzyme_hybrid_scaffold_M1/alignfinal/E_CTTAAG_Q_NGScontigs_A_HYBRID"
  r <- run_alignfinal_two_pass(hybrid, ngs, prefix, args1=align_final1, args2=align_final2, ra = '/home/users3/tanantharaman/tools/7238/RefAligner')

  #testing empty second pass
  hybrid <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/TestRun/SBJ/PB40x/two/Sandwich2/test_2col_c1.cmap"
  ngs <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/TestRun/SBJ/PB40x/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID_1st_pass_r.cmap"
  prefix <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/TestRun/SBJ/PB40x/TGH_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID"
  r <- run_alignfinal_two_pass(hybrid, ngs, prefix, args1=align_final1, args2=align_final2, ra = '/home/users3/tanantharaman/tools/7238/RefAligner')
}

test_greedy_align_final <- function(){
  align_final_file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID_multimatches.xmap"
  process_align_final(align_final_file)
}

test_get_unalign_contigs_single <- function(){
  align_final_file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID.xmap"
  qcmap <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/fa2cmap/a.lines.fasta.cut_BSPQI_0kb_0labels.cmap"
  a <- get.unalign.contigs.single(align_final_file, qcmap, Id.col = 'Qry')
  #testing empty alignment
  align_final_file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID.xmap1"
  qcmap <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/fa2cmap/a.lines.fasta.cut_BSPQI_0kb_0labels.cmap"
  a <- get.unalign.contigs.single(align_final_file, qcmap, Id.col = 'Qry')

}

tet_two_color_xmap_filter <- function(){
  align1 <- as.data.table(rbind(readxmap("~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/HummingBird/Illumn/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID_1st_pass.xmap"),
                                readxmap("~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/HummingBird/Illumn/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID_2nd_pass.xmap")))
  align2 <- as.data.table(rbind(readxmap("~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/HummingBird/Illumn/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSSSI_Q_NGScontigs_A_HYBRID_1st_pass.xmap"),
                                 readxmap("~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/HummingBird/Illumn/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSSSI_Q_NGScontigs_A_HYBRID_2nd_pass.xmap")))
  two.color.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/HummingBird/Illumn/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_E_BSSSI_Q_NGS_A_HYBRID_hash1.xmap"
  two.color.file2 <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/HummingBird/Illumn/two_enzyme_hybrid_scaffold_M1/alignfinal/E_BSPQI_E_BSSSI_Q_NGS_A_HYBRID_hash2.xmap"
  align_twocolor1 <- as.data.table(readxmap(two.color.file))
  align_twocolor2 <- as.data.table(readxmap(two.color.file2))
  xmap_filter_twocolor(align1, align2, align_twocolor1, align_twocolor2, two.color.file, two.color.file2, Confidence.T = 11)
}


test_filterAlign_SegDup <- function(){
  align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/CTTAAG/hybrid_scaffolds_M1/CTTAAG_EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_BNGcontigs_NGScontigs.xmap"
  cmap.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/CTTAAG/hybrid_scaffolds_M1/CTTAAG_EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_BNGcontigs_NGScontigs_q.cmap"
  align <- as.data.table(readxmap(align.file))
  maps <- as.data.table(readcmap(cmap.file))

  align_filter <- filterAlignInSegDup(align, maps, map.type = 'Qry')
  align[!(XmapEntryID %in% align_filter$XmapEntryID), sum(abs(QryStartPos - QryEndPos))]

  #BSPQI
  align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/BSPQI_non_SegDup/align1/align1.xmap"
  cmap.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/BSPQI_non_SegDup/align1/align1_q.cmap"
  align_BSPQI <- as.data.table(readxmap(align.file))
  maps_BSPQI <- as.data.table(readcmap(cmap.file))

  align_filter_BSPQI <- filterAlignInSegDup(align_BSPQI, maps_BSPQI, map.type = 'Qry')
  align_BSPQI[!(XmapEntryID %in% align_filter_BSPQI$XmapEntryID), sum(abs(QryStartPos - QryEndPos))]

  #compare segdup regions
  segDup_DLE1 <- getSegDupRegions(maps)
  segDup_BSPQI <- getSegDupRegions(maps_BSPQI)

  #testing case without mask
  align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_two/BSPQI/align1/align1.xmap"
  cmap.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_two/BSPQI/align1/align1_q.cmap"
  align_BSPQI <- as.data.table(readxmap(align.file))
  maps_BSPQI <- as.data.table(readcmap(cmap.file))

  align_filter_BSPQI <- filterAlignInSegDup(align_BSPQI, maps_BSPQI, map.type = 'Qry')
  align_BSPQI[!(XmapEntryID %in% align_filter_BSPQI$XmapEntryID), sum(abs(QryStartPos - QryEndPos))]

  #filtering ref
  align.file <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/hybrid_scaffolds/hybrid_scaffold_020518/output/align_final/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_not_filtered.xmap"
  cmap.file <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/hybrid_scaffolds/hybrid_scaffold_020518/output/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_r.cmap"
  align <- as.data.table(readxmap(align.file))
  maps <- as.data.table(readcmap(cmap.file))

  align_filter <- filterAlignInSegDup(align, maps, map.type = 'Ref')
  align[!(XmapEntryID %in% align_filter$XmapEntryID), sum(abs(QryStartPos - QryEndPos))]


  align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/BSPQI_non_SegDup/align_final/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_not_filtered.xmap"
  cmap.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/BSPQI_non_SegDup/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_r.cmap"
  align_BSPQI <- as.data.table(readxmap(align.file))
  maps_BSPQI <- as.data.table(readcmap(cmap.file))

  align_filter_BSPQI <- filterAlignInSegDup(align_BSPQI, maps_BSPQI)
  align_BSPQI[!(XmapEntryID %in% align_filter_BSPQI$XmapEntryID), sum(abs(QryStartPos - QryEndPos))]


}


test_run_sandwich_withSegDupMask <- function(){
  align1.file <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/hybrid_scaffolds/hybrid_scaffold_020518/output/align_final/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_not_filtered.xmap"
  align2.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/BSPQI_non_SegDup/align_final/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_not_filtered.xmap"
  cmap1.file <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/hybrid_scaffolds/hybrid_scaffold_020518/output/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_r.cmap"
  cmap2.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/BSPQI_non_SegDup/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_r.cmap"

  a <- filterAlignFinal1(align1.file, align2.file, cmap1.file, cmap2.file, "")
  exportXmap(a[[1]], "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/two_enzyme_hybrid_scaffold_M1/Sandwich2_nonSegdup/test_1_fullAln.xmap")
  exportXmap(a[[2]], "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/two_enzyme_hybrid_scaffold_M1/Sandwich2_nonSegdup/test_2_fullAln.xmap")

  hybrid.alignfinal1 <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/two_enzyme_hybrid_scaffold_M1/Sandwich2_nonSegdup/test_1_fullAln.xmap"
  hybrid.alignfinal2 <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/two_enzyme_hybrid_scaffold_M1/Sandwich2_nonSegdup/test_2_fullAln.xmap"

  sandwich2.out <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_SegDupMask/two_enzyme_hybrid_scaffold_M1/Sandwich2_nonSegdup/"

  run.sandwich(hybrid.alignfinal1, hybrid.alignfinal2, cmap1.file, cmap2.file, sandwich2.out, NULL, NULL, score.T = 13)

  align1.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_FullSegDupMask/RunWithAsm_7420/BSPQI/align_final_M1/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_not_filtered.xmap"
  align2.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_FullSegDupMask/RunWithAsm_7420/CTTAAG/align_final_M1/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_not_filtered.xmap"
  cmap1.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_FullSegDupMask/RunWithAsm_7420/BSPQI/align_final_M1/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_r.cmap"
  cmap2.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_FullSegDupMask/RunWithAsm_7420/CTTAAG/align_final_M1/EXP_REFINEFINAL1_bppAdjust_cmap_GCA_002077035_2_NA12878_prelim_2_1_genomic_fna_NGScontigs_HYBRID_SCAFFOLD_r.cmap"

  align1.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x/test/BSPQI/align_final_M1/test_cmap_bppAdjust_cmap_SBJ_40x_p_contigsGr40kb_fa_NGScontigs_HYBRID_SCAFFOLD_not_filtered.xmap"
  align2.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x/test/BSSSI/align_final_M1/EXP_REFINEFINAL1_bppAdjust_cmap_SBJ_40x_p_contigsGr40kb_fa_NGScontigs_HYBRID_SCAFFOLD_not_filtered.xmap"
  cmap1.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x/test/BSPQI/align_final_M1/test_cmap_bppAdjust_cmap_SBJ_40x_p_contigsGr40kb_fa_NGScontigs_HYBRID_SCAFFOLD_r.cmap"
  cmap2.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x/test/BSSSI/align_final_M1/EXP_REFINEFINAL1_bppAdjust_cmap_SBJ_40x_p_contigsGr40kb_fa_NGScontigs_HYBRID_SCAFFOLD_r.cmap"

  a <- filterAlignFinal1(align1.file, align2.file, cmap1.file, cmap2.file, "")

}

test_align_fraction_filter <- function(){
  align1.file <- '~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_FullSegDupMask/RunWithNewBNX/two_enzyme_hybrid_scaffold_M1/Sandwich2/test_1_fullAln.xmap'
  align2.file <- '~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_FullSegDupMask/RunWithNewBNX/two_enzyme_hybrid_scaffold_M1/Sandwich2/test_2_fullAln.xmap'

  align1.file <- '~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TestRunWith_FullSegDupFix_Asm/WithNewBnx/two_enzyme_hybrid_scaffold_M1/Sandwich2/test_1_fullAln.xmap'
  align2.file <- '~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TestRunWith_FullSegDupFix_Asm/WithNewBnx/two_enzyme_hybrid_scaffold_M1/Sandwich2/test_2_fullAln.xmap'

  align1.file <- '~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_Illmn2/TestRunWith_FullSegDupMask/two_enzyme_hybrid_scaffold_M1/Sandwich2/test_1_fullAln.xmap'
  align2.file <- '~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_Illmn2/TestRunWith_FullSegDupMask/two_enzyme_hybrid_scaffold_M1/Sandwich2/test_2_fullAln.xmap'

  align1 <- as.data.table(readxmap(align1.file))
  align2 <- as.data.table(readxmap(align2.file))

  align.filter <- filterSparsePairs(align1, align2)
  hist(align.filter$pairstats[,fr])
  print(paste0("Total pairs ", nrow(align.filter$pairstats), " Filter out ", nrow(align.filter$pairstats[fr < 0.1])))
  align.filter$pairstats[fr < 0.1]

}






