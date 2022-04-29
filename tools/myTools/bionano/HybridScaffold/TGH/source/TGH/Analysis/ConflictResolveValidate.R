#This script contain function that use simulation to validate/benchmark conflict resolution in hybrid-scaffold
library(data.table)
#library(VennDiagram)
#source('~/workspace/codebase/MoleculeSimulator/IO/readmaps.R')
#source('~/workspace/codebase/MoleculeSimulator/Utility/Utility.R')
#source('~/workspace/codebase/MoleculeSimulator/SingleMoleculeSimulator.R')

#This is an updated addChimera function with finer control of chimera cutting-points
addChimera <- function(molCmap1, molCmap2, chimID, minLeng=150e3){
  addChimera.helper <- function(molCmap1, molCmap2){
    cut.point1 <- 0
    cut.point2 <- molCmap2$ContigLength[1]
    maxAttempts <- 1000
    attempts <-1
    while(cut.point1 + (molCmap2$ContigLength[1] - cut.point2) < minLeng && attempts < maxAttempts){
      cut.point1 <- runif(1, min=molCmap1$Position[[1]]+1, max=tail(molCmap1$Position,1)-1)
      cut.point2 <- runif(1, min=molCmap2$Position[[1]]+1, max=tail(molCmap2$Position,1)-1)
      attempts <- attempts + 1
    }
    labelInd1 <- which(molCmap1$Position < cut.point1)
    labelInd2 <- which(molCmap2$Position > cut.point2)
    Chimera <- rep(NA, length(labelInd1) + length(labelInd2)) #adding a column that keep track of where the chimera mols comes from
    #computing new positions for chimera molecules
    newPositions <- molCmap2$Position[labelInd2] - cut.point2 + cut.point1
    molCmap2$Position[labelInd2] <- newPositions
    #if(ncol(molCmap1) > 12 | ncol(molCmap2) > 12){
    #  print(head(molCmap1))
    #  print(head(molCmap2))
    #}
    chimera.cmap <- rbind(molCmap1[labelInd1,], molCmap2[labelInd2,])
    chimera.cmap$Chimera <- c(rep(molCmap1$CMapId[1], length(labelInd1)),
                              rep(molCmap2$CMapId[1], length(labelInd2)))
    #chimera.cmap$ContigLength <- cut.point1 + (molCmap2$ContigLength[1]-cut.point2)
    chimera.cmap$ContigLength <- max(chimera.cmap$Position)
    chimera.cmap$CMapId <- chimID
    chimera.cmap <- bookeepMolecules(chimera.cmap)
    return(chimera.cmap)
  }
  #there is a 50-50 chance of which molecules we use as the first one
  if(runif(1) > 0.5){
    return(addChimera.helper(molCmap2, molCmap1))
  }else{
    return(addChimera.helper(molCmap1, molCmap2))
  }
}

run.conflictResolve.validate <- function(){
  #seq.file <- '~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/fa2cmap2color/soap.result.scafSeq.fa.cutted_BSPQI_BSSSI_0kb_0labels.cmap'

  #simulated NGS chimera
  seq.file <- '~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/fa2cmap2color/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels.cmap'
  output.dir <- '~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/DisCovar_NGS/'
  chim.cmap.sep <- createNGSChimera(seq.file, output.dir)

  conflicts.detect <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/DisCovar_NGS/BSPQI_conflict_detect/"
  conflicts.detect2 <-"~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/DisCovar_NGS/BSSSI_conflict_detect/"
  chim.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/DisCovar_NGS/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels.cmap_chim_cmap.RData"

  eval <- eval.conflict.detect(conflicts.detect, conflicts.detect2, chim.file)
  chim.cmap <- readRDS(chim.file)
  chim.cmap.sep <- getSeparateChannel(chim.cmap)

  aligned.boundary.dist <- get.Align.boundary(align1, chim.cmap.sep$C1, eval[[1]])
  aligned.boundary.dist2 <- get.Align.boundary(align1.2, chim.cmap.sep$C2, eval[[2]])

  #simulated BNG chimera
  bng.file <- '~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/BSPQI/assignAlignType/cut_conflicts_M1/bn_cut_pre_exclude.cmap'
  output.dir <- '~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/DisCovar_NGS/'
  chim.bng <- createBNGChimera(bng.file, output.dir)

}

test.ngs.len.dist <- function(){
  ngs.seq.file1 <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/fa2cmap2color/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels.cmap"
  ngs.seq.file2 <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/fa2cmap2color/soap.result.scafSeq.fa.cutted_BSPQI_BSSSI_0kb_0labels.cmap"
  seq.cmap1 <- as.data.table(readcmap(ngs.seq.file1))
  seq.cmap2 <- as.data.table(readcmap(ngs.seq.file2))

  Nbreaks <- nrow(seq.cmap2[SiteID==1 & ContigLength > 10e3])
  lens2 <- seq.cmap2[SiteID==1, ContigLength]
  break.prob <- seq.cmap2[SiteID==1, .(CMapId, ContigLength/sum(ContigLength)*Nbreaks)]
  break.prob$nbreaks <- as.integer(round(break.prob$V2))
  breaks.loc <- break.prob[,{p<-sort(runif(nbreaks)); c(head(p,1), diff(p), 1-tail(p,1))}, by=CMapId]
  contig.len <- seq.cmap2[SiteID==1, .(CMapId, ContigLength)]
  setkey(contig.len, 'CMapId')
  breaks.loc$Len <- breaks.loc$V1*contig.len[breaks.loc$CMapId, ContigLength]
  lens2.split <- c(breaks.loc$Len, seq.cmap2[SiteID==1][break.prob$nbreaks == 0 & ContigLength > 10000, ContigLength])

  lens2.split <- rep(seq.cmap2[SiteID==1, ContigLength], nbreaks.per.contig)
  lens2.split <- lens2.split*runif(length(lens2.split))
  lens2 <- c(lens2.split, seq.cmap[SiteID==1][nbreaks.per.contig== 0, ContigLength])

  lens2 <- seq.cmap2[SiteID==1 & ContigLength > 10e3, {r <- runif(length(ContigLength)); c(r*ContigLength, (1-r)*ContigLength)}]

  lens2 <- seq.cmap2[SiteID==1 & ContigLength > 10000, ContigLength]
  lens1 <- seq.cmap1[SiteID==1 & ContigLength > 10000, ContigLength]

  lens2.2 <- lens2*runif(length(lens2))
  g <- ggCompareDensity(log10(lens1[lens1 > 10]), log10(lens2.2[lens2.2 >10]), breaks = seq(1, 7, 0.1))

  g <- ggCompareDensity(log10(lens1[lens1 > 10]), log10(lens2.split[lens2.split > 10]), breaks = seq(1,7,0.1))

}


createNGSChimera <- function(seq.file, output.dir, min.len=100e3){
  seq.cmap <- as.data.table(readcmap(seq.file))
  seq.cmap.large <- seq.cmap[ContigLength > 100e3 & NumSites > 0]
  chim.cmap <- getChimera(seq.cmap.large, subset=1.0)

  chim.cmap.c1 <- chim.cmap[LabelChannel==1 | LabelChannel==0]
  chim.cmap.c1 <- bookeepMolecules(chim.cmap.c1)
  chim.cmap.c2 <- chim.cmap[LabelChannel==2 | LabelChannel==0]
  chim.cmap.c2[LabelChannel > 0]$LabelChannel <- 1
  chim.cmap.c2 <- bookeepMolecules(chim.cmap.c2)

  chim.cmap[SiteID==1, hist(NumSites, seq(0,max(NumSites),1))]
  plot(ecdf(chim.cmap[SiteID==1, NumSites]), xlim=(c(0,50)))

  filename <- basename(seq.file)
  saveRDS(chim.cmap, paste0(output.dir, "/", filename, "_chim_cmap.RData"))
  exportCMap(chim.cmap.c1, paste0(output.dir, "/"), paste0(filename, "_c1.cmap"))
  exportCMap(chim.cmap.c2, paste0(output.dir, "/"), paste0(filename, "_c2.cmap"))
  return(list(chim.cmap.c1, chim.cmap.c2))
}

createBNGChimera <- function(cmap.file, output.dir, min.len=100e3){
  bng.cmap <- as.data.table(readcmap(cmap.file))
  chim.list <- lapply(seq(1:5), function(x){
    chim.cmap <- getChimera(bng.cmap, subset=1.0)
  })
  maxId <- 0
  for(i in 1:(length(chim.list))){
    chim <- chim.list[[i]]
    chim$CMapId <- as.character(as.numeric(chim$CMapId) + maxId)
    chim.list[[i]] <- chim
    maxId <- max(as.numeric(chim$CMapId))
    print(maxId)
  }
  chim.cmap <- rbindlist(chim.list)
  #hist(chim.cmap[SiteID==1, NumSites], seq(0,700,1))
  plot(ecdf(chim.cmap[SiteID==1, NumSites]), xlim=(c(0,50)))

  filename <- basename(cmap.file)
  saveRDS(chim.cmap, paste0(output.dir, "/", filename, "_chim_cmap.RData"))
  exportCMap(chim.cmap, paste0(output.dir, "/"), paste0(filename, "_chim.cmap"))
  return(list(chim.cmap))
}


eval.conflict.detect <- function(conflict.detect.path, conflict.detect.path2, chim.rdata.file){
  #browser()
  conflicts.table.file <- paste0(conflict.detect.path, "/assignAlignType/conflicts.txt")
  conflicts.table.file2 <- paste0(conflict.detect.path2, "/assignAlignType/conflicts.txt")
  align1.file  <- paste0(conflict.detect.path, "/align1/align1.xmap")
  align1.file2 <- paste0(conflict.detect.path2, "/align1/align1.xmap")
  chim.cmap <- readRDS(chim.rdata.file)
  chim.cmap.sep <- getSeparateChannel(chim.cmap)

  chim.loc <- chim.cmap[, Position[Chimera != Chimera[[1]]][[1]], by=CMapId]

  check1 <- check.chim.brkpt(align1.file, chim.loc, span.region = 50000)
  check2 <- check.chim.brkpt(align1.file2, chim.loc, span.region = 50000)

  #setnames(chim.loc, c("RefcontigID", "ChimLoc"))
  #align1 <- merge(align1, chim.loc, by='RefcontigID')
  #align1.span.chim <- align1[RefStartPos < ChimLoc - 20000 & RefEndPos > ChimLoc + 20000]

  detected1 <- getDetectionTable(conflicts.table.file, align1.file, chim.cmap.sep$C1, type='NGS')
  detected2 <- getDetectionTable(conflicts.table.file2, align1.file2, chim.cmap.sep$C2, type='NGS')

  detected1 <- merge(detected1, check1, by='CMapId')
  detected2 <- merge(detected2, check2, by='CMapId')

  detected <- merge(detected1, detected2, by='CMapId', all = T)
  detected[is.na(V1.x), V1.x := V1.y]
  detected[is.na(V2.x), V2.x := V2.y]
  detected[is.na(detected)] <- FALSE

  bins <- c(seq(1,20,1), seq(20, 40, 5),500)
  detected$Bins <- findInterval(detected$V1.x, bins)
  detected$Bins2 <- findInterval(detected$V1.y, bins)
  detected.rate <- detected[,list(sum(Detected.x | Detected.y)/sum(Aligned.x | Aligned.y), sum(Detected.x)/sum(Aligned.x)), by=Bins]
  detected.rate2 <- detected[,list(sum(Detected.x | Detected.y)/sum(Aligned.x | Aligned.y), sum(Detected.y)/sum(Aligned.y)), by=Bins2]
  setnames(detected.rate2, "Bins2", "Bins")
  detected.rate2 <- detected.rate2[Bins > 0]
  #detected.rate <- merge(detected.rate, detected.rate2, by="Bins")
  detected.rate <- detected.rate[order(Bins)]
  detected.rate2 <- detected.rate2[order(Bins)]

  v1 <- detected[V1.x > 0, draw.pairwise.venn(sum(Detected.x), sum(Detected.y), sum(Detected.x & Detected.y),
                                        category = c("BSPQI-Detected", "BSSSI-Detected"))]

  v2 <- detected[V1.x > 0, draw.pairwise.venn(sum(Aligned.x), sum(Aligned.y), sum(Aligned.x & Aligned.y),
                                        category = c("BSPQI-Aligned", "BSSSI-Aligned"))]
  grid.newpage()
  grid.draw(v1)
  grid.newpage()
  grid.draw(v2)

  N <- nrow(detected.rate)
  N2 <- nrow(detected.rate2)
  bin.labels <- bins
#  df <- data.frame(Labels=c(rep(bin.labels[detected.rate$Bins], 2), rep(bin.labels[detected.rate$Bins2], 2)),
#                   DetectedRate=c(detected.rate$V1, detected.rate$V2, detected.rate2$V1, detected.rate2$V2),
#                   Enzyme=c(rep("Two-enzyme (by BSPQI label)", N), rep("BSPQI", N), rep("Two-enzyme (by BSSSI label)", N2), rep("BSSSI", N2)))
  df <- data.frame(Labels=c(rep(bin.labels[detected.rate$Bins], 2)),
                   DetectedRate=c(detected.rate$V1,  detected.rate$V2),
                   Enzyme=c(rep("Two-enzyme", N), rep("BSPQI", N)))


  g <- ggplot(df)
  g <- g + geom_point(aes(x=Labels, y=DetectedRate, colour=Enzyme))
  g <- g + geom_line(aes(x=Labels, y=DetectedRate, colour=Enzyme), data=df[grepl('Two-enzy', df$Enzyme),])
  g <- g + geom_line(aes(x=Labels, y=DetectedRate, colour=Enzyme), data=df[!grepl('Two-enzy', df$Enzyme),], linetype=2)
  g <- g +  xlab("Number of labels in conflict regions") + ylab("Fractions of conflicts detected")
  g <- g + scale_x_continuous(limit=c(1,20)) + theme(text=element_text(size=18))
  plot(g)
  return(list(detected1, detected2))

#  df <- rbind(
#              data.frame(Labels=detected[Detected.x == TRUE, V1.x], Enzyme="BSPQI"),
#              data.frame(Labels=detected[Detected.y==TRUE, V1.y], Enzyme="BSSSI"),
#              data.frame(Labels=detected[(Detected.x | Detected.y) == TRUE, V1.x], Enzyme="Two-Enzyme"))
#  g <- ggplot(df)
#  g <- g + geom_histogram(aes(x=Labels, fill=Enzyme), position = 'dodge')
#  plot(g)
}

check.chim.brkpt <- function(xmap.file , chim.loc, span.region=50000){
  align1 <- as.data.table(readxmap(xmap.file))
  setkey(align1, 'RefcontigID')
  setkey(chim.loc, 'CMapId')
  align1$ChimLoc <- chim.loc[align1$RefcontigID]$V1
  translated.pos <- translatePositions(align1$ChimLoc, align1, 1)
  align1$QrySpan <- align1[,QryLen - translated.pos > span.region & translated.pos > span.region]
  span.count <- align1[,sum(QrySpan), by=RefcontigID]
  setnames(span.count, c("CMapId", "SpanCount"))
  return(span.count)
}

#evaluate the accuracy of alignment boundary (i.e. how close is the end of alignment match the true breakpointo f a chimeric contigs)
get.Align.boundary <- function(xmap, chim, detected.table){
  chim <- removeDresLabels(chim)
  chim <- removeNonChim(chim)
  chim <- bookeepMolecules(chim)
  chim.ind <- chim[, SiteID[Chimera != Chimera[[1]]][[1]], by=CMapId]
  setkey(chim.ind, 'CMapId')
  chim.ind <- chim.ind[detected.table[SpanCount >0, CMapId]]
  setkey(chim.ind, 'CMapId')
  xmap <- processXmapAlignments(xmap)
  aligned.dist <- lapply(1:nrow(xmap),
                         function(i){
                           ind <- chim.ind[xmap$RefcontigID[[i]], V1]
                           if(is.na(ind)){
                             return(-1000)
                           }
                           dist <- xmap$alignRefInd[[i]] - ind
                           left.aligned.cnt <- sum(xmap$alignRefInd[[i]] < ind)
                           right.aligned.cnt <- sum(xmap$alignRefInd[[i]] > ind)
                           if(left.aligned.cnt > right.aligned.cnt){
                             return(max(xmap$alignRefInd[[i]]) - (ind-1)) #if aligned from left contig should end at one label before chim-labels
                           }else{
                             return(ind - min(xmap$alignRefInd[[i]]))
                           }

                        })
  xmap$aligned.dist <- unlist(aligned.dist)
  aligned.dist <- xmap[,aligned.dist[which.min(abs(aligned.dist))], by=RefcontigID]
  aligned.dist <- aligned.dist[V1 > -1000]
  #browser()
  h <- aligned.dist[,hist(V1, (seq(min(V1)-1, max(V1),1)+0.5))]
  df <- data.frame(mids=h$mids, density=h$density)
  g <- ggplot()
  g <- g + geom_bar(aes(x=mids, y=density), stat='identity', data=df) + scale_x_continuous(limit=c(-10,10))
  df$labels <- paste0(format(h$density*100, digits = 1), "%")
  g <- g + geom_text(data=df[df$density > 0.0001,], aes(x=mids, y=density, label=labels, nudge_y=0.5), vjust=0)
  g <- g + xlab('Alignment Boundary Errors') + ylab('Density')
  plot(g)
  #compute  FP errors
  df$errors <- cumsum(h$density)
  df$ErrorLabels <- paste0(format(df$errors*100, digits=2), "%")
  #browser()
  g <- ggplot(df[df$mids < 0,])
  g <- g + geom_bar(aes(x=mids, y=errors), stat='identity') + scale_x_continuous(limit=c(-20, 0))
  g <- g + geom_text(data=df[df$mids < 0,], aes(x=mids, y=errors, label=ErrorLabels, nudge_y=0.5), vjust=0)
  g <- g + xlab("Align boundary errors") + ylab("Accumulative errors")
  plot(g)
  return(aligned.dist)
}

#selecting only chim that has at least one label from the secondary molecules
removeNonChim <- function(chim.cmap){
  chim.count <- chim.cmap[, length(unique(Chimera)), by=CMapId]
  chim.count <- chim.count[V1 > 1]
  return(chim.cmap[CMapId %in% chim.count$CMapId])
}

removeDresLabels <- function(cmap, res=0.9){
  dist <- cmap[,c(5e6, diff(Position)), by=CMapId]
  cmap.res0.9 <- cmap[dist$V1 > res*500 | LabelChannel == 0] #we do not remove end labels, which is not a real lables
  return(cmap.res0.9)
}

#generate a detection table to check which chimeric contigs/maps are detected by conflict detection
getDetectionTable <- function(conflict.table.file, align.file, chim.cmap, type){
  chim.cmap.res0.9 <- removeDresLabels(chim.cmap)
  label.count <- chim.cmap.res0.9[LabelChannel > 0, list(min(table(Chimera)), max(table(Chimera))), by=CMapId]
  setkey(label.count, 'CMapId')
  if(length(chim.cmap$Chimera)){
    chim.count <- chim.cmap.res0.9[LabelChannel >0, length(unique(Chimera)), by=CMapId] #selecting only chim that has at least one label from the secondary molecules
    chim.count <- chim.count[V1 > 1]
    label.count <- label.count[chim.count$CMapId]
  }
  #browser()
  conflicts.table <- read.table(conflict.table.file, stringsAsFactors = F, header = T)
  align1 <- as.data.table(readxmap(align.file))
  if(type=='NGS'){
    ngs.aligned.ids <- data.table(CMapId=as.character(unique(align1$RefcontigID)), Aligned=T)
    ngs.detected.ids <- data.table(CMapId=as.character(unique(conflicts.table$refId)), Detected=T)
  }else if(type=='BNG'){
    ngs.aligned.ids <- data.table(CMapId=as.character(unique(align1$QryContigID)), Aligned=T)
    ngs.detected.ids <- data.table(CMapId=as.character(unique(conflicts.table$qryId)), Detected=T)
  }else{
    stop('Unknown type: type can either be BNG or NGS')
  }
  detected <- merge(label.count, ngs.aligned.ids, by='CMapId')
  detected <- merge(detected, ngs.detected.ids, all.x = T, by='CMapId')
  detected[is.na(Detected)]$Detected <- FALSE
  return(detected)
}

getDetectRate <- function(detected.table, bins, label){
  #bins <- c(seq(1,20,1), seq(21, 50, 5))
  detected.table$Bins <- findInterval(detected.table$V1, bins)
  detected.rate <- detected.table[,sum(Detected)/sum(Aligned), by=Bins]
  detected.rate <- detected.rate[order(Bins)]

  N <- nrow(detected.rate)
  bin.labels <- (bins[-1] + head(bins, length(bins)-1) )/2 - 0.5
  df <- data.frame(Labels=rep(bin.labels[detected.rate$Bins], 1),
                   DetectedRate=c(detected.rate$V1),
                   Labels=c(rep(label, N)))
  return(df)
}

eval.bngconflict.detect <- function(){
  conflicts.table.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/PacBio_NGS/assignAlignType/conflicts.txt"
  conflicts.table.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/TwoEnzymeConflictResolve/PacBioNGS/CombinedCutConflicts/Combined_cut_conflicts_status_1.txt"
  align1.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/PacBio_NGS/align1/align1.xmap"
  chim.cmap <- readRDS('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/DisCovar_NGS/Chims5000/bn_cut_pre_exclude.cmap_chim_cmap.RData')
  #chim.loc <- chim.cmap.res0.9[, Position[Chimera != Chimera[[1]]][[1]], by=CMapId]
  #setnames(chim.loc, c("RefcontigID", "ChimLoc"))
  #align1 <- merge(align1, chim.loc, by='RefcontigID')
  #align1.span.chim <- align1[RefStartPos < ChimLoc - 20000 & RefEndPos > ChimLoc + 20000]

  detected <- getDetectionTable(conflicts.table.file, align1.file, chim.cmap, type='BNG')
  detected2 <- getDetectionTable(conflicts.table.file2, align1.file, chim.cmap, type='BNG')
  detected[is.na(detected)] <- FALSE
  detected2[is.na(detected)] <- FALSE

  bins <- c(seq(1,60,2), seq(61, 100, 10))
  df <- getDetectRate(detected, bins=bins, label = "BSPQI-only")
  df2 <- getDetectRate(detected2, bins = bins, label = "Two-enzyme")

  g <- ggplot(df)
  g <- g + geom_point(aes(x=Labels, y=DetectedRate, colour=Labels.1))
  g <- g + geom_line(aes(x=Labels, y=DetectedRate, colour=Labels.1))
  g <- g + geom_point(aes(x=Labels, y=DetectedRate, colour=Labels.1), data = df2)
  g <- g + geom_line(aes(x=Labels, y=DetectedRate, colour=Labels.1), data = df2)

  g <- g + scale_y_continuous(limit=c(0,1))
  plot(g)
}

eval.twoEnzyme.bngConflictDetect <- function(){
  conflicts.table.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/DisCovar_NGS/BSPQI_GM_conflict_Detect/assignAlignType/conflicts.txt"
  align1.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/DisCovar_NGS/BSPQI_GM_conflict_Detect/align1/align1.xmap"
  conflicts.table.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/BSSSI/assignAlignType/conflicts.txt"
  align1.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/BSSSI/align1/align1.xmap"


}



analyzeAlign1 <- function(){
  align1 <- readxmap("~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/BSPQI_conflict_detect/align1/align1.xmap")
  align1 <- as.data.table(align1)

}

analyzeConflictFromHapHybrid <- function(){
  NGS.align.final1 <- as.data.table(readxmap("~/RemoteServer/home/users/jlee/data/20160727_Aaegypti_BspQI_55C72C_hybrid_seq_spec_config_T11/hybrid_scaffolds/EXP_REFINEFINAL1_Aedes_aegypti_PacBio_assembly_primaryContigs_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"))
  NGS.align.final2 <- as.data.table(readxmap("~/RemoteServer/home/users/jlee/data/20160727_Aaegypti_BspQI_55C72C_hybrid_p_plus_a_seq_spec_config_T11//hybrid_scaffolds/EXP_REFINEFINAL1_Aedes_aegypti_PacBio_p_plus_a_assembly_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"))
  align2.only <- setdiff(NGS.align.final2$QryContigID, NGS.align.final1$QryContigID)

  align1 <- as.data.table(readxmap("~/RemoteServer/home/users/jlee/data/20160727_Aaegypti_BspQI_55C72C_hybrid_seq_spec_config_T11//align1/align1.xmap"))
  align21 <- as.data.table(readxmap("~/RemoteServer/home/users/jlee/data/20160727_Aaegypti_BspQI_55C72C_hybrid_p_plus_a_seq_spec_config_T11/align1/align1.xmap"))

  align2.only <- setdiff(align21$RefcontigID, align1$RefcontigID)

}

#this function analyze the conflicts identified across runs
conflict.overlap.across.runs <- function(){
  conflict1.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/PB_randcut/RandCut_EmpSBJ40x/CombinedConflictsCut/Combined_conflicts_cut_status_BSPQI.txt"
  conflict2.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/ReRun_12082016/CombinedConflictsCut/Combined_conflicts_cut_status_BSPQI.txt"

  conflict1.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/PB30x/CombinedConflictsCut/Combined_conflicts_cut_status_BSPQI.txt"
  conflict2.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/PB60x/CombinedConflictsCut/Combined_conflicts_cut_status_BSPQI.txt"

  conflicts1 <- as.data.table(read.hybrid.conflictfile(conflict1.file))
  conflicts2 <- as.data.table(read.hybrid.conflictfile(conflict2.file))

  conflicts1 <- conflicts1[qry_leftBkpt_toCut == 'cut' | qry_rightBkpt_toCut == 'cut' | ref_leftBkpt_toCut == 'cut' | ref_rightBkpt_toCut == 'cut']
  conflicts2 <- conflicts2[qry_leftBkpt_toCut == 'cut' | qry_rightBkpt_toCut == 'cut' | ref_leftBkpt_toCut == 'cut' | ref_rightBkpt_toCut == 'cut']

  overlap.cnt <- overlap.conflicts.qry(conflicts1, conflicts2)
  overlap.cnt2 <- overlap.conflicts.qry(conflicts2, conflicts1)


  paste0("Total conflicts in conflict1: ", nrow(conflicts1))
  paste0("Total conflicts in conflict2: ", nrow(conflicts2))
  paste0("Total conflicts in conflict1: ", nrow(conflicts1[refId != -1]))
  paste0("Total conflicts in conflict2: ", nrow(conflicts2[refId != -1]))

  paste0("Overlap within 10kb: ", sum(overlap.cnt > 0))
  paste0("Overlap within 10kb: ", sum(overlap.cnt2 > 0))

  overlap.cnt <- overlap.conflicts.qry(conflicts1, conflicts2, 500e3)
  overlap.cnt2 <- overlap.conflicts.qry(conflicts2, conflicts1, 500e3)
  paste0("Overlap within 100kb: ", sum(overlap.cnt > 0))
  paste0("Overlap within 100kb: ", sum(overlap.cnt2 > 0))

}

#this function compute overlap between two sets of conflicts, using qry as the reference point
overlap.conflicts.qry <- function(conflict1, conflict2, max.dist = 10e3){
  #browser()
  conflict.intersect <- lapply(1:nrow(conflict1),
                               function(i){
                                 conflict <- conflict1[i]
                                 overlap.left <- conflict2[qryId == conflict$qryId
                                                           & (conflict$leftQryBkpt > 0
                                                              & (abs(conflict$leftQryBkpt - leftQryBkpt) < max.dist
                                                                 | abs(conflict$leftQryBkpt - rightQryBkpt) < max.dist)), which=T]
                                 overlap.right <- conflict2[qryId == conflict$qryId
                                                            & (conflict$rightQryBkpt > 0
                                                               &  (abs(conflict$rightQryBkpt - leftQryBkpt) < max.dist)
                                                               | abs(conflict$rightQryBkpt - rightQryBkpt) < max.dist), which=T]
                                 overlap <- c(overlap.left, overlap.right)
                               })
  intersect.count <- unlist(lapply(conflict.intersect, function(inter){length(inter)}))
  return(intersect.count)
}




























