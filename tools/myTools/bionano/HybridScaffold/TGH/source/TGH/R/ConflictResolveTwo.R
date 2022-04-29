#This script provide functions that perform a combined conflict resolution on two single-enzyme hybrid scaffolds
#library(igraph)
#library(data.table)

#source("/home/jian/workspace/codebase/MoleculeSimulator/Utility/Utility.R")
#source("/home/jian/workspace/codebase/MoleculeSimulator/IO/readmaps.R")

#analyze the intersection of two conflict table from hybrid-scaffolds
conflict.intersect <- function(hybrid.dir1, hybrid.dir2, xmap1.file1=NULL, xmap2.file=NULL, print.summary=TRUE){
  #browser()
  conflict1.file <- paste0(hybrid.dir1, "/assignAlignType/cut_conflicts/conflicts_cut_status.txt")
  conflict2.file <- paste0(hybrid.dir2, "/assignAlignType/cut_conflicts/conflicts_cut_status.txt")

  conflict1 <- as.data.table(read.hybrid.conflictfile(conflict1.file))
  conflict2 <- as.data.table(read.hybrid.conflictfile(conflict2.file))


  xmap1.file <- paste0(hybrid.dir1, "/align1/align1.xmap")
  xmap2.file <- paste0(hybrid.dir2, "/align1/align1.xmap")

  #BNG cuts comparison
  qry_cut1 <- which(conflict1$qry_leftBkpt_toCut == 'cut'
                   | conflict1$qry_rightBkpt_toCut == 'cut')
  qry_cut2 <- which(conflict2$qry_leftBkpt_toCut == 'cut'
                    | conflict2$qry_rightBkpt_toCut == 'cut')
  qry.cut.overlap <- overlap.conflicts(conflict1[qry_cut1], conflict2[qry_cut2])

  #NGS cut comparison
  ref_cut1 <- which(conflict1$ref_leftBkpt_toCut == 'cut'
                    | conflict1$ref_rightBkpt_toCut == 'cut')
  ref_cut2 <- which(conflict2$ref_leftBkpt_toCut == 'cut'
                    | conflict2$ref_rightBkpt_toCut == 'cut')
  #should add another overlap stats base on location also
  ref.cut.overlap <- overlap.conflicts(conflict1[ref_cut1], conflict2[ref_cut2])

  xmap1 <- as.data.table(readxmap(xmap1.file))
  xmap2 <- as.data.table(readxmap(xmap2.file))
  xmap1$RefcontigID <- as.character(xmap1$RefcontigID)
  xmap2$RefcontigID <- as.character(xmap2$RefcontigID)

  #browser()
  overlap <- data.frame(BNG_cut_overlap_allDist=conflict1[qry_cut1][,sum(refId %in% conflict2[qry_cut2, refId])],
                  BNG_cut_overlap=sum(qry.cut.overlap > 0),
                  NGS_cut_overlap_allDist=conflict1[ref_cut1][,sum(refId %in% conflict2[ref_cut2, refId])],
                  NGS_cut_overlap=sum(ref.cut.overlap > 0),
                  BNG1_NGS2_cut_overlap=sum(overlap.conflicts(conflict1[qry_cut1], conflict2[ref_cut2]) > 0),
                  BNG2_NGS1_cut_overlap=sum(overlap.conflicts(conflict2[qry_cut2], conflict1[ref_cut1]) > 0))

  #checking if second-enzyme map overlap break points
  map.overlap.check1 <- check.brkpt.overlap(xmap2, conflict1)
  map.overlap.check2 <- check.brkpt.overlap(xmap1, conflict2)

  conflict1 <- cbind(conflict1, data.table(leftBrkptOverlap = map.overlap.check1$V1))
  conflict1 <- cbind(conflict1, data.table(rightBrkptOverlap = map.overlap.check1$V2))
  conflict2 <- cbind(conflict2, data.table(leftBrkptOverlap = map.overlap.check2$V1))
  conflict2 <- cbind(conflict2, data.table(rightBrkptOverlap = map.overlap.check2$V2))

  #checking if second-enzyme map cover breaks points, note this diff from above as it doens't need to cross and span the region
  map.cov.check1 <- check.brkpt.overlap(xmap2, conflict1, span.region = -50000)
  map.cov.check2 <- check.brkpt.overlap(xmap1, conflict2, span.region = -50000)

  conflict1 <- cbind(conflict1, data.table(leftBrkptCov = map.cov.check1$V1))
  conflict1 <- cbind(conflict1, data.table(rightBrkptCov = map.cov.check1$V2))
  conflict2 <- cbind(conflict2, data.table(leftBrkptCov = map.cov.check2$V1))
  conflict2 <- cbind(conflict2, data.table(rightBrkptCov = map.cov.check2$V2))

  #checking overhang alignments of maps at break point
  self.check1 <- check.brkpt.supportAln(xmap1, conflict1)
  self.check2 <- check.brkpt.supportAln(xmap2, conflict2)
  #browser()
  conflict1 <- cbind(conflict1, self.check1)
  conflict2 <- cbind(conflict2, self.check2)
  second.map.check1 <- check.brkpt.supportAln(xmap2, conflict1)
  setnames(second.map.check1, c("leftBkOverhang2", "rightBkOverhang2"))
  second.map.check2 <- check.brkpt.supportAln(xmap1, conflict2)
  setnames(second.map.check2, c("leftBkOverhang2", "rightBkOverhang2"))

  conflict1 <- cbind(conflict1, second.map.check1)
  conflict2 <- cbind(conflict2, second.map.check2)


  N_qry_cut1 <- sum(conflict1$qry_leftBkpt_toCut == 'cut') + sum(conflict1$qry_rightBkpt_toCut == 'cut')
  N_qry_cut2 <- sum(conflict2$qry_leftBkpt_toCut == 'cut') + sum(conflict2$qry_rightBkpt_toCut == 'cut')
  N_ref_cut1 <- sum(conflict1$ref_leftBkpt_toCut == 'cut') + sum(conflict1$ref_rightBkpt_toCut == 'cut')
  N_ref_cut2 <- sum(conflict2$ref_leftBkpt_toCut == 'cut') + sum(conflict2$ref_rightBkpt_toCut == 'cut')

  summary <- data.table(totalBNGCut=c(N_qry_cut1, N_qry_cut2),
               BNGCutSpan = c(nrow(conflict1[qry_cut1][leftBrkptOverlap > 0 | rightBrkptOverlap > 0]),
                           nrow(conflict2[qry_cut2][leftBrkptOverlap > 0 | rightBrkptOverlap > 0])),
                           BNGCutCov = c(nrow(conflict1[qry_cut1][leftBrkptCov > 0 | rightBrkptCov > 0]),
                              nrow(conflict2[qry_cut2][leftBrkptCov > 0 | rightBrkptCov > 0])),
               BNGOverhang = c(nrow(conflict1[qry_cut1][leftBkOverhang > 0 | rightBkOverhang > 0]),
                               nrow(conflict2[qry_cut2][leftBkOverhang > 0 | rightBkOverhang > 0])),
               BNGOverhang2 = c(nrow(conflict1[qry_cut1][leftBkOverhang2 > 0 | rightBkOverhang2 > 0])  ,
                                nrow(conflict2[qry_cut2][leftBkOverhang2 > 0 | rightBkOverhang2 > 0])),
               totalNGSCut=c(N_ref_cut1, N_ref_cut2),
               NGSCutSpan = c(nrow(conflict1[ref_cut1][leftBrkptOverlap > 0 | rightBrkptOverlap > 0]),
                              nrow(conflict2[ref_cut2][leftBrkptOverlap > 0 | rightBrkptOverlap > 0])),
               NGSCutCov = c(nrow(conflict1[ref_cut1][leftBrkptCov > 0 | rightBrkptCov > 0]),
                              nrow(conflict2[ref_cut2][leftBrkptCov > 0 | rightBrkptCov > 0])),
               NGSOverhang = c(nrow(conflict1[ref_cut1][leftBkOverhang > 0 | rightBkOverhang > 0]),
                               nrow(conflict2[ref_cut2][leftBkOverhang > 0 | rightBkOverhang > 0])),
               NGSOverhang2 = c(nrow(conflict1[ref_cut1][leftBkOverhang2 > 0 | rightBkOverhang2 > 0]),
                               nrow(conflict2[ref_cut2][leftBkOverhang2 > 0 | rightBkOverhang2 > 0])))
  if(print.summary){
    print(summary)
    print(overlap)
  }
  return(list(summary=summary, conflict1=conflict1, conflict2=conflict2, overlap=overlap))

}

#This function compute conflict overlap between two conflict table
overlap.conflicts <- function(conflict1, conflict2, max.dist = 10e3){
  #browser()
  conflict.intersect <- lapply(1:nrow(conflict1),
                               function(i){
                                 conflict <- conflict1[i]
                                 overlap.left <- conflict2[refId == conflict$refId
                                                      & (conflict$leftRefBkpt > 0
                                                        & (abs(conflict$leftRefBkpt - leftRefBkpt) < max.dist
                                                           | abs(conflict$leftRefBkpt - rightRefBkpt) < max.dist)), which=T]
                                 overlap.right <- conflict2[refId == conflict$refId
                                                      & (conflict$rightRefBkpt > 0
                                                         &  (abs(conflict$rightRefBkpt - leftRefBkpt) < max.dist)
                                                             | abs(conflict$rightRefBkpt - rightRefBkpt) < max.dist), which=T]
                                 overlap <- c(overlap.left, overlap.right)
                               })
  intersect.count <- unlist(lapply(conflict.intersect, function(inter){length(inter)}))
  return(intersect.count)
}


test.conflict.intersect <- function(){
  #conflict1.file <- "~/RemoteServer/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_davidJaffe_031815/hybrid_scaffold/hybrid_scaffold_012416/thomas_haplotype_bsssi_vs_discovar/hybrid_scaffolds/conflicts_cut_status.txt"
  #conflict2.file <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_davidJaffe_031815/hybrid_scaffold/hybrid_scaffold_012416/thomas_haplotype_vs_discovar/hybrid_scaffolds/conflicts_cut_status.txt"
  #xmap1.file <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_davidJaffe_031815/hybrid_scaffold/hybrid_scaffold_012416/thomas_haplotype_bsssi_vs_discovar/align1/align1.xmap"
  #xmap2.file <- "~/RemoteServer/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_davidJaffe_031815/hybrid_scaffold/hybrid_scaffold_012416/thomas_haplotype_vs_discovar/align1/align1.xmap"

  hybrid1 <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x_old/BSPQI/"
  hybrid2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/ConflictDetectEval/DisCovar_NGS/BSSSI_Conflict_Detect/"

  conflict.stats <- conflict.intersect(hybrid1, hybrid2)

}

#merge conflict file from two-hybrid scaffold runs
#using NGS as the common reference, translate the conflict detected in one hybrid run
#to another, so both hybrid will end up making the same consistent set of cut in the final hybrid
merge.conflict.file <- function(hybrid.dir1, hybrid.dir2, span.region = 500){
  conflict1.file <- paste0(hybrid.dir1, "/assignAlignType/cut_conflicts/conflicts_cut_status.txt")
  conflict2.file <- paste0(hybrid.dir2, "/assignAlignType/cut_conflicts/conflicts_cut_status.txt")
  conflict1 <- as.data.table(read.hybrid.conflictfile(conflict1.file))
  conflict2 <- as.data.table(read.hybrid.conflictfile(conflict2.file))
  xmap1.file <- paste0(hybrid.dir1, "/align1/align1.xmap")
  xmap2.file <- paste0(hybrid.dir2, "/align1/align1.xmap")
  xmap1 <- as.data.table(readxmap(xmap1.file))
  xmap2 <- as.data.table(readxmap(xmap2.file))
  cmap1.file <- paste0(hybrid.dir1, "/align1/align1_q.cmap")
  cmap2.file <- paste0(hybrid.dir2, "/align1/align1_q.cmap")
  qcmap1 <- as.data.table(readcmap(cmap1.file))
  qcmap2 <- as.data.table(readcmap(cmap2.file))
  ref1.file <- paste0(hybrid.dir1, "/align1/align1_r.cmap")
  ref2.file <- paste0(hybrid.dir2, "/align1/align1_r.cmap")
  rcmap1 <- as.data.table(readcmap(ref1.file))
  rcmap2 <- as.data.table(readcmap(ref2.file))

  check1 <- check.brkpt.overlap(xmap2, conflict1, span.region =  span.region)
  check2 <- check.brkpt.overlap(xmap1, conflict2, span.region =  span.region)

  #ngs.cut.1 <- conflict1$ref_leftBkpt_toCut == 'cut' | conflict1$ref_rightBkpt_toCut=='cut'
  #ngs.cut.2 <- conflict2$ref_leftBkpt_toCut == 'cut' | conflict2$ref_rightBkpt_toCut=='cut'

  conflict2.from.1 <- getTranslatedConflict(conflict1, xmap2, qcmap2, rcmap2)
  conflict1.from.2 <- getTranslatedConflict(conflict2, xmap1, qcmap1, rcmap1)

  #we should not turn on cut on BNG when doing translation
  #conflict2.from.1 <- conflict2.from.1[qry_leftBkpt_toCut != 'cut' & qry_rightBkpt_toCut != 'cut']
  #conflict1.from.2 <- conflict1.from.2[qry_leftBkpt_toCut != 'cut' & qry_rightBkpt_toCut != 'cut']

  #determining possible conflicting contig pair from hybrid-scaffold conflicts
  #conflict1.pair <- merge(conflict1, conflict2.from.1, by='refId', allow.cartesian = TRUE)
  #conflict1.pairId <- unique(paste0(conflict1.pair$qryId.x, '_', conflict1.pair$qryId.y))

  #conflict2.pair <- merge(conflict2, conflict1.from.2, by='refId', allow.cartesian = TRUE)
  #conflict2.pairId <- unique(paste0(conflict2.pair$qryId.y, '_', conflict1.pair$qryId.x))
  conflict1 <- rbind(conflict1, conflict1.from.2)
  conflict2 <- rbind(conflict2, conflict2.from.1)

  return(list(conflict1, conflict2))
}


#translate the coordiante in the conflict file from one enzyme map to another enzyme
#map using NGS to BNG map alignments as a common reference system
getTranslatedConflict <- function(conflict, xmap.dt, qcmap, rcmap){

  if(nrow(conflict) == 0 || is.null(conflict)){
    return(conflict)
  }

  aligns1 <- getSpanAlns(xmap.dt, conflict$refId, conflict$leftRefBkpt, span.region = 500)
  aligns2 <- getSpanAlns(xmap.dt, conflict$refId, conflict$rightRefBkpt, span.region = 500)
  translate.helper <- function(position, align){
    #browser()
    if(nrow(align) > 0){
      pos <- translatePositionsLocal(position, qcmap, rcmap, align[1], 1)
      return(data.table(qryId=align[1]$QryContigID, qryPos=pos))
    }else{
      return(data.table(qryId='-1', qryPos=-1))
    }
  }
  translated.pos1 <- rbindlist(mapply(translate.helper,  conflict$leftRefBkpt, aligns1, SIMPLIFY = F))
  translated.pos2 <- rbindlist(mapply(translate.helper, conflict$rightRefBkpt, aligns2, SIMPLIFY = F))
  #browser()
  conflict$leftQryBkpt <- translated.pos1$qryPos
  conflict$rightQryBkpt <- translated.pos2$qryPos
  conflict$qryId <- ifelse(as.numeric(translated.pos1$qryId) > 0, translated.pos1$qryId, translated.pos2$qryId)
  #there are cases where NGS is mapped to one map in one case but two in the other case we need to created
  #new cut conflict entries for those
  inds.diff.id <- which(as.numeric(translated.pos1$qryId) * as.numeric(translated.pos2$qryId) > 0
                        & translated.pos1$qryId != translated.pos2$qryId)
  if(length(inds.diff.id) > 0){
    conflict.to.add <- conflict[inds.diff.id]
    conflict.to.add$qryId <- translated.pos2$qryId[inds.diff.id]
    conflict <- clearRightBrkPtEntry(conflict, inds.diff.id)
    conflict.to.add <- clearLeftBrkPtEntry(conflict.to.add, 1:nrow(conflict.to.add))
    conflict <- rbind(conflict, conflict.to.add)
  }
  return(conflict)
}

#test merging conflict file
test.conflict.merge <- function(){
  #hybrid.dir1 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridScaffold_NonHaplotype/ReRun_03242016/BSPQI"
  #hybrid.dir2 <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridScaffold_NonHaplotype/BSSSI/"
  hybrid.dir1 <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x_old/BSPQI/"
  hybrid.dir2 <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/SBJ/PB40x_old/BSSSI/"

  #hybrid.dir1 <- "~/RemoteServer/home/users/jlee/data/20160727_Aaegypti_BspQI_55C72C_hybrid_p_plus_a_seq_spec_config_T11/"
  #hybrid.dir2 <- "~/RemoteServer/home/users/jlee/data/20160727_Aaegypti_BssSI_hybrid_p_plus_a_seq_spec_config_T11/"

  conflicts.merged <- merge.conflict.file(hybrid.dir1, hybrid.dir2)
  conflicts.merged.stats <- conflict.intersect(hybrid.dir1, hybrid.dir2)

  merged.out.file1 <- '~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/CombinedCutConflicts/Combined_cut_conflicts_status_1.txt'
  merged.out.file2 <- '~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/CombinedCutConflicts/Combined_cut_conflicts_status_2.txt'

  #out.df <- data.frame(lapply(conflicts.merged[[1]], as.character), stringsAsFactors=FALSE)
  #write.table(out.df, file=merged.out.file1, quote = F, row.names = F, sep='\t')

  #out.df <- data.frame(lapply(conflicts.merged[[2]], as.character), stringsAsFactors=FALSE)
  #write.table(out.df, file=merged.out.file2, quote = F, row.names = F, sep='\t')
}

#clear the brkpt entries of a conflict file
#i.e setting a coordinate to -1 and cut status to okay
clearRightBrkPtEntry <- function(conflict, inds){
  conflict$rightRefBkpt[inds] <- -1
  conflict$rightQryBkpt[inds] <- -1
  conflict$ref_rightBkpt_toCut[inds] <- 'okay'
  conflict$qry_rightBkpt_toCut[inds] <- 'okay'
  return(conflict)
}

clearLeftBrkPtEntry <- function(conflict, inds){
  conflict$leftRefBkpt[inds] <- -1
  conflict$leftQryBkpt[inds] <- -1
  conflict$ref_leftBkpt_toCut[inds] <- 'okay'
  conflict$qry_leftBkpt_toCut[inds] <- 'okay'
  return(conflict)
}


#check which alignment spans a particular position on the reference
getSpanAlns <- function(xmap.dt, refId, refPos, span.region=30000, topAln=0){
  Pos.list <- data.frame(refId, refPos)
  aligns <- lapply(1:nrow(Pos.list),
                   function(i){
                     align <- xmap.dt[RefcontigID== refId[i] & RefStartPos + span.region <= refPos[i]
                             & RefEndPos - span.region >= refPos[i]]
                     #browser()
                     if(topAln > 0 && nrow(align) > topAln){
                       return(align[1:topAln])
                     }else{
                       return(align)
                     }
                   }
  )
  return(aligns)
}


#Using the NGS as reference, check if a conflict breakpoint is cover
#by alternative map alignment (e.g. another contig from haplotype-assembly or second enzyme).
#A cut in BNG map support by map alignment support this cut, alternatively a cut in NGS coverd
#by map alignment indicate disagreement with the cut
#NOTE may consier output this status anyway even if decide not cut, for comprehensiveness
check.brkpt.overlap <- function(xmap.dt, conflict, span.region=50000){
  if(nrow(conflict)==0){
    return(data.table(V1=list(), V2=list()))
  }
  setkey(xmap.dt, RefcontigID)

  check.helper <- function(conflict){
    if(conflict$qry_leftBkpt_toCut=='cut' | conflict$ref_leftBkpt_toCut=='cut'){
      cover.left <- nrow(xmap.dt[RefcontigID==conflict$refId][RefStartPos + span.region <= conflict$leftRefBkpt
                                                              & RefEndPos - span.region >= conflict$leftRefBkpt])
    }else{
      cover.left <- -1
    }
    if(conflict$qry_rightBkpt_toCut == 'cut' | conflict$ref_rightBkpt_toCut=='cut'){
      cover.right <- nrow(xmap.dt[RefcontigID==conflict$refId][RefStartPos + span.region <= conflict$rightRefBkpt
                                                             & RefEndPos - span.region >= conflict$rightRefBkpt])
    }else{
       cover.right <- -1
    }
    return(list(cover.left, cover.right))
  }

  cover.list <- lapply(1:nrow(conflict),
                       function(i){
                         check.helper(conflict[i,])
                      }
                  )
  return(rbindlist(cover.list))
}


#check for pair of over-hanging alignments supporting/disagreeing the break point
check.brkpt.supportAln <- function(xmap.dt, conflict, overhang.region=50000){
  if(nrow(conflict)==0){
    return(data.table(leftBkOverhang=list(),
                               rightBkOverhang=list()))
  }
  setkey(xmap.dt, RefcontigID)

  check.helper <- function(Bkpt.toCut, refId, brkpt.position){
    #browser()
    if(Bkpt.toCut=='cut'){
      #alignment extends from the left and pass breakpoint
      support.left <- nrow(xmap.dt[RefcontigID==refId][abs(RefEndPos - brkpt.position) < 30000
                                                              & RefStartPos < 20000
                                                              & ((QryLen - QryEndPos > overhang.region & Orientation == '+')
                                                              | (QryEndPos > overhang.region & Orientation == '-'))])
       #alignment extends from the right pass breakpoint
      support.right <- nrow(xmap.dt[RefcontigID==refId][abs(RefStartPos - brkpt.position) < 30000
                                                               & abs(RefEndPos - RefLen) < 20000
                                                               & ((QryLen - QryStartPos > overhang.region & Orientation == '-')
                                                              | (QryStartPos > overhang.region & Orientation == '+'))])
      return(c(support.left, support.right))
    }else{
      return(c(-1,-1))
    }
  }

  cover.list <- lapply(1:nrow(conflict),
                       function(i){
                         left.check <- check.helper(conflict[i,]$ref_leftBkpt_toCut, conflict[i,]$refId, conflict[i,]$leftRefBkpt)
                         right.check <- check.helper(conflict[i,]$ref_rightBkpt_toCut, conflict[i,]$refId, conflict[i,]$rightRefBkpt)
                         return(c(left.check, right.check))
                       }
  )
  #return(rbindlist(cover.list))
  combined <- {
    if(length(cover.list)==1){
      data.table(t(cover.list[[1]]))
    }else{
      data.table(Reduce(function(x,y){rbind(x,y)}, cover.list))
    }
  }

  #setnames(combined, c('leftbk.leftSupport', 'leftbk.rightSupport', 'rightbk.leftSupport', 'rightbk.rightSupport'))
  #browser()
  check.both.end <- function(c1, c2){
    p <- c1*c2
    s <- c1+c2
    p[s < 0] <- -1
    p[p > 0 & s > 0] <- s[p > 0 & s > 0]
    return(p)
  }
  return(data.table(leftBkOverhang=check.both.end(combined$V1, combined$V2),
                    rightBkOverhang=check.both.end(combined$V3, combined$V4)))
}


#find conflicts between two scaffold
find.scaffold.conflict <- function(){
  conflicts.connection <- find.connection.conflict()
  gap.conflicts <- find.gap.conflicts()
}


find.connection.conflict <-function(){
  agp.file1 <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_davidJaffe_031815/hybrid_scaffold/hybrid_scaffold_012416/thomas_haplotype_bsssi_vs_discovar/agp_fasta/EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.agp"
  agp.file2 <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_davidJaffe_031815/hybrid_scaffold/hybrid_scaffold_012416/thomas_haplotype_vs_discovar/agp_fasta/EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.agp"
  xmap1.file1<- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_davidJaffe_031815/hybrid_scaffold/hybrid_scaffold_012416/thomas_haplotype_bsssi_vs_discovar/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_BNGcontigs_HYBRID_SCAFFOLD.xmap"
  xmap2.file <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_davidJaffe_031815/hybrid_scaffold/hybrid_scaffold_012416/thomas_haplotype_vs_discovar/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_BNGcontigs_HYBRID_SCAFFOLD.xmap"
  name.file1 <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_davidJaffe_031815/hybrid_scaffold/hybrid_scaffold_012416/thomas_haplotype_bsssi_vs_discovar/fa2cmap/a.lines_BSSSI_0kb_0labels_key.txt"
  name.file2 <- "~/RemoteServer//mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_davidJaffe_031815/hybrid_scaffold/hybrid_scaffold_012416/thomas_haplotype_vs_discovar/fa2cmap/a.lines_BSPQI_0kb_0labels_key.txt"

  g <- get.graphfrom.agp(agp.file1, directed=F)
  E(g)$color <- 'green'
  g2 <- get.graphfrom.agp(agp.file2, directed=F)
  E(g2)$color <- 'blue'
  g.merged <- union(g, g2, byname = TRUE)
  v.degree <- degree(g.merged)

  triangles <- triangles(g.merged)
  squares <- unlist(find.cycles(g.merged, 4))
  pentagon <- unlist(find.cycles(g.merged, 5))

  v.split <- setdiff(names(v.degree[v.degree==3]), c(names(triangles), names(squares), names(pentagon)))

  c1 <- components(g)
  c2 <- components(g2)
  c <- components(g.merged)

}

find.cycles <- function(graph, k) {
  ring <- graph.ring(k)
  graph.get.subisomorphisms.vf2(graph, ring)
}

#flip the path coordinate so rel_start always less than rel_end
flip.rel.coord <- function(path.coord){
  selected <- path.coord$orientation == -1
  tmp <- path.coord$rel_start[selected]
  path.coord$rel_start[selected] <- path.coord$rel_end[selected]
  path.coord$rel_end[selected] <- tmp
  return(path.coord)
}

#a generlaize version of above function, given a arbirary test function
#we swap the two columns that satisfy the condition given in test.fun
flip.columns <- function(dt, col1, col2, test.fun){
  dt <- as.data.frame(dt)
  selected <- test.fun(dt[,col1], dt[,col2])
  tmp <- dt[selected, col1]
  dt[selected, col1] <- dt[selected, col2]
  dt[selected, col2] <- tmp
  return(as.data.table(dt))
}

find.conflict.alignment <- function(indata, pos_offsets){
  if(nrow(indata) < 2 || nrow(pos_offsets) < 2){
    return(list(data.table(), data.table()))
  }
  align.merged <- as.data.table(indata)
  #align.merged <- align.merged[Confidence.x > 11 & Confidence.y > 11]
  align.merged$RefcontigID.x <- as.character(align.merged$RefcontigID.x)
  align.merged$RefcontigID.y <- as.character(align.merged$RefcontigID.y)
  id_str <- unique(align.merged$id_str)
  connected.pairs <- strsplit(pos_offsets$id_str, "_")

  pairstats <- as.data.table(pos_offsets)
  pairstats$slope <- as.numeric(pairstats$slope)
  pairstats$intercept <- as.numeric(pairstats$intercept)
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

  #finding all single-ton pairs
  single.id.x <- table(pairstats.ordered$RefcontigID.x)
  single.id.x <- names(single.id.x[single.id.x==1])
  single.id.y <- table(pairstats.ordered$RefcontigID.y)
  single.id.y <- names(single.id.y[single.id.y==1])

  pairstats.ordered <- pairstats.ordered[!(RefcontigID.x %in% single.id.x) & !(RefcontigID.y %in% single.id.y)]
  #finding alignment that is conflict
  Npairs <- nrow(pairstats.ordered)
  Npairs <- ifelse(Npairs > 0, Npairs, 1)
  conflict.aligns.list1 <- lapply(1:Npairs,
                                 function(i){
                                   pairstat <- pairstats.ordered[i]
                                   inds <- align.merged[RefcontigID.x == pairstat$RefcontigID.x
                                        & RefStartPos.x >= pairstat$RelStart1 & RefEndPos.x <= pairstat$RelEnd1
                                        & id_str != pairstat$id_str, which=T]
                                   compare <- align.merged[inds][,pCompareInterval(QryStartPos.x, QryEndPos.x, QryStartPos.y, QryEndPos.y)]
                                   inds <- inds[compare$dist < -10e3]
                                   cover.ids <- align.merged[RefcontigID.x == pairstat$RefcontigID.x & RefcontigID.y == pairstat$RefcontigID.y, QryContigID]
                                   ind.removed <- which(align.merged[inds]$QryContigID %in% cover.ids)
                                   if(length(ind.removed) > 0){
                                     inds <- inds[-ind.removed]
                                   }
                                   return(inds)
                                 })
  counts1 <- sapply(conflict.aligns.list1, length)
  conflict.aligns.list2 <- lapply(1:Npairs,
                                  function(i){
                                    pairstat <- pairstats.ordered[i]
                                    inds <- align.merged[RefcontigID.y == pairstat$RefcontigID.y
                                                 & RefStartPos.y >= pairstat$RelStart2 & RefEndPos.y <= pairstat$RelEnd2
                                                 & id_str != pairstat$id_str, which=T]
                                    compare <- align.merged[inds][,pCompareInterval(QryStartPos.x, QryEndPos.x, QryStartPos.y, QryEndPos.y)]
                                    inds <- inds[compare$dist < -10e3]
                                    cover.ids <- align.merged[RefcontigID.x == pairstat$RefcontigID.x & RefcontigID.y == pairstat$RefcontigID.y, QryContigID]
                                    ind.removed <- which(align.merged[inds]$QryContigID %in% cover.ids)
                                    if(length(ind.removed) > 0){
                                      inds <- inds[-ind.removed]
                                    }
                                    return(inds)
                                  })
  counts2 <- sapply(conflict.aligns.list2, length)
  print(paste0('Found conflicts from alignment 1: ', sum(counts1 > 0), ' from alignment 2: ', sum(counts2 > 0)))

  #browser()
  #we need to find conflict region for each conflict pair separately
  #for a pair of conflicting connection of contigs, find out the conflict breakpoint by analyzing NGS to contig alignments
  getConflictRegion <- function(agree.pair.id, conflict.pair.id, align.ids, RefStartPos, RefEndPos){
      conflict.Inds <- align.ids == conflict.pair.id
      agree.Inds <- align.ids == agree.pair.id
      start.conflict <- min(RefStartPos[conflict.Inds])
      end.conflict <- max(RefEndPos[conflict.Inds])
      start.agree <- min(RefStartPos[agree.Inds])
      end.agree <- max(RefEndPos[agree.Inds])
      compare <- compareInterval(start.conflict, end.conflict, start.agree, end.agree)
      #conflict on the right side
      if(compare[2] == 1){
        return(data.table(id_str = agree.pair.id,
                          ConflictStartPos=end.agree, ConflictEndPos=start.conflict, Overlap=compare[1], conflict_pair=conflict.pair.id))
      #conflict on the left side
      }else if(compare[2] == -1){
        return(data.table(id_str = agree.pair.id, ConflictStartPos=end.conflict, ConflictEndPos=start.agree, Overlap=compare[1], conflict_pair=conflict.pair.id))
      }else{
        warning(paste0("conflict ", conflict.pair.id ," is embedded in agreement regions, do not current handle such case"))
        #browser()
        return(data.table(id_str=agree.pair.id, ConflictStartPos=-1, ConflictEndPos=-1, Overlap=compare[1], conflict_pair=conflict.pair.id))
      }
  }

  #This find conflict region in a batch mode for every pair of connected contigs on sandwich
  getConflictRegionBatch <- function(pairstat, conflict.align.inds, align.ids, alignStartPos, alignEndPos){
    conflicted.pairs <- unique(align.merged$id_str[conflict.align.inds])
    conflict.regions <- lapply(conflicted.pairs,
                               function(conflicted.pair){
                                 getConflictRegion(pairstat$id_str, conflicted.pair,
                                                   align.ids, alignStartPos, alignEndPos)
                               })
    return(rbindlist(conflict.regions))
  }

  #find conflict regions from conflict alignments
  conflict.region.list <- lapply(1:length(conflict.aligns.list1),
                             function(i){
                               inds <- conflict.aligns.list1[[i]]
                               if(length(inds) > 0){
                                 region <- getConflictRegionBatch(pairstats.ordered[i], inds, align.merged$id_str, align.merged$RefStartPos.x, align.merged$RefEndPos.x)
                                 region$ngs_align_id <- align.merged[region$id_str, QryContigID, mult='first']
                                 return(region)
                               }
                             })
  conflict.region.combined <- rbindlist(conflict.region.list)

  if(nrow(conflict.region.combined) > 0){
    conflict.region.combined <- conflict.region.combined[Overlap > -20000] #remove potentially haplotype maps
    print(paste0("Remove ", nrow(conflict.region.combined[Overlap < -20000]), " conflicts due to possible haplotype maps"))
  }


  conflict.region.list2 <- lapply(1:length(conflict.aligns.list2),
                                  function(i){
                                    inds <- conflict.aligns.list2[[i]]
                                    if(length(inds)==0){
                                      return(NULL)
                                    }else{
                                      region <- getConflictRegionBatch(pairstats.ordered[i], inds, align.merged$id_str, align.merged$RefStartPos.y, align.merged$RefEndPos.y)
                                      region$ngs_align_id <- align.merged[region$id_str, QryContigID, mult='first']
                                      return(region)
                                    }
                                  })

  conflict.region.combined2 <- rbindlist(conflict.region.list2)
  if(nrow(conflict.region.combined2) > 0){
    conflict.region.combined2 <- conflict.region.combined2[Overlap > -20000]
    print(paste0("Remove ", nrow(conflict.region.combined2[Overlap < -20000]), " conflicts due to possible haplotype maps"))
  }

  setkey(pairstats.ordered, id_str)
  if(nrow(conflict.region.combined) > 0){
    setkey(conflict.region.combined, id_str)
    conflict.region.combined <- pairstats.ordered[conflict.region.combined]
    conflict.region.combined[,ConflictStartPos2 := slope*ConflictStartPos + intercept]
    conflict.region.combined[,ConflictEndPos2 := slope*ConflictEndPos + intercept]
    conflict.region.combined[ConflictStartPos2 > RefLen.y, ConflictStartPos2 := RefLen.y - 5000]
    conflict.region.combined[ConflictEndPos2 > RefLen.y, ConflictEndPos2 := RefLen.y - 5000]
    conflict.region.combined[ConflictStartPos2 < 0, ConflictStartPos2 := 5000]
    conflict.region.combined[ConflictEndPos2 < 0, ConflictEndPos2 := 5000]
  }

  if(nrow(conflict.region.combined2) > 0){
    setkey(conflict.region.combined2, id_str)
    setnames(conflict.region.combined2, c("ConflictStartPos", "ConflictEndPos"), c("ConflictStartPos2", "ConflictEndPos2"))
    conflict.region.combined2 <- pairstats.ordered[conflict.region.combined2]
    conflict.region.combined2[,ConflictStartPos := 1/slope*(ConflictStartPos2 - intercept)]
    conflict.region.combined2[,ConflictEndPos := 1/slope*(ConflictEndPos2 - intercept)]
    conflict.region.combined2[ConflictStartPos > RefLen.x, ConflictStartPos := RefLen.x - 5000]
    conflict.region.combined2[ConflictEndPos > RefLen.x, ConflictEndPos := RefLen.x - 5000]
    conflict.region.combined2[ConflictStartPos < 0, ConflictStartPos := 5000]
    conflict.region.combined2[ConflictEndPos < 0, ConflictEndPos := 5000]
    conflict.region.combined2 <- flip.columns(conflict.region.combined2, 'ConflictStartPos', 'ConflictEndPos', `>`)
  }
  return(list(conflict.region.combined, conflict.region.combined2))
}

#Generate the cut-status file for hybrid scaffold
#type specifies where the cut should be happening, either 'qry' and 'ref'
generate.conflict.table <- function(conflict.region, type='qry'){
  if(type != 'qry' && type != 'ref'){
    stop('Unknown type, must be either qry or ref')
  }
  if(nrow(conflict.region) < 1){
    return(list(data.table(), data.table()))
  }
  conflict.region <- conflict.region[ConflictStartPos > 0 & ConflictEndPos > 0 & ConflictStartPos2 > 0 & ConflictEndPos > 0]
  if(nrow(conflict.region) < 1){
    return(list(data.table(), data.table()))
  }

  generate.conflict.table.helper <- function(Ids, ngs.Ids, Positions1, Positions2, type='qry'){
    if(type=='qry'){
      conflict.table <- data.frame(xMapId = seq(1:nrow(conflict.region)),
                                   refQry ='ref', refId=ngs.Ids, leftRefBkpt='-1', rightRefBkpt='-1', alignmentOrientation='+',
                                   ref_leftBkpt_toCut = 'okay', ref_rightBkpt_toCut = 'okay', ref_toDiscard = 'okay',
                                   refQry = 'qry', qryId =Ids,  leftQryBkpt = Positions1, rightQryBkpt=Positions2,
                                   alignmentOrientation='+', qry_leftBkpt_toCut='cut', qry_rightBkpt_toCut='okay', qry_toDiscard='okay')
      conflict.table <- get.unique.df(conflict.table, 'qryId')
    }
    if(type=='ref'){
      conflict.table <- data.frame(xMapId = seq(1:nrow(conflict.region)),
                                   refQry = 'ref', refId = Ids, leftRefBkpt = Positions1, rightRefBkpt=Positions2,
                                   alignmentOrientation='+', ref_leftBkpt_toCut='cut', ref_rightBkpt_toCut='okay', ref_toDiscard='okay',
                                   refQry ='qry', qryId=ngs.Ids, leftQryBkpt='-1', rightQryBkpt='-1', alignmentOrientation='+',
                                   qry_leftBkpt_toCut = 'okay', qry_rightBkpt_toCut = 'okay', qry_toDiscard = 'okay')
      conflict.table <- get.unique.df(conflict.table, 'refId')
    }
    return(conflict.table)
  }
  #browser()
  conflict.table <- generate.conflict.table.helper(conflict.region$RefcontigID.x, conflict.region$ngs_align_id, conflict.region$ConflictStartPos,
                                                   conflict.region$ConflictEndPos, type)
  conflict.table2 <- generate.conflict.table.helper(conflict.region$RefcontigID.y, conflict.region$ngs_align_id, conflict.region$ConflictStartPos2,
                                                    conflict.region$ConflictEndPos2, type)
  return(list(conflict.table, conflict.table2))
}

test.find.conflict.align <- function(){
  align <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestUWPrimates/Clint_TGH_BNGbased_2/TGH_M1/Sandwich2/test_indata.RData"
  pos.offset <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestUWPrimates/Clint_TGH_BNGbased_2/TGH_M1/Sandwich2/test_pos_offsets.RData"
  load(align)
  load(pos.offset)
  conflict.regions <- find.conflict.alignment(indata, pos_offsets)
  type <- 'qry'
  conflict.tables <- generate.conflict.table(conflict.regions[[1]], type)
  conflict.tables2 <- generate.conflict.table(conflict.regions[[2]], type)

  conflict.tables[[1]] <- rbind(conflict.tables[[1]], conflict.tables2[[1]])
  conflict.tables[[2]] <- rbind(conflict.tables[[2]], conflict.tables2[[2]])

  #clean up header to restore the proper header column for cut_conflict file, data.table do not allow duplicate colnames
  conflict.tables[[1]] <- cleanup.header(conflict.tables[[1]])
  conflict.tables[[2]] <- cleanup.header(conflict.tables[[2]])


}

#gets a unique dataframe
get.unique.df <- function(df, key){
  names <- colnames(df)
  dt <- as.data.table(df)
  setkeyv(dt, key)
  dt <- unique(dt, by=key)
  df <- as.data.frame(dt)
  colnames(df) <- names
  return(df)
}

#this function pick best alignment from a merged alignment possible contain more than one alignment
#per query, it is used in two-color alignFinal to pick an alignment for a NGS contigs
remove.duplicated.align <- function(xmap, delta.score = 2, align.len=0.9, bestRef = T){
  if(nrow(xmap) < 2){
    return(xmap)
  }
  setkey(xmap, 'QryContigID')
  score.filter <- xmap[,max(Confidence) - Confidence < delta.score, by=QryContigID]$V1
  align.len.filter <- xmap[,abs(QryStartPos - QryEndPos)/max(abs(QryStartPos-QryEndPos)), by=QryContigID]
  align.len.filter <- align.len.filter$V1 >= align.len
  xmap.filtered <- xmap[score.filter & align.len.filter]
  #choosing alignment with longest bng contigs
  if(bestRef){
    ref.len.filter <- xmap.filtered[, RefLen==max(RefLen), by=QryContigID]$V1
    xmap.filtered <- xmap.filtered[ref.len.filter]
  }
  return(xmap.filtered)
}

test.remove.dupllicated <- function(){
  xmap.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3/Alignments/HybridM1/BSPQI/a.lines_BSPQI_0kb_0labels_cut_1_fullAln.xmap"
  xmap <- as.data.table(readxmap(xmap.file))
  xmap <- remove.duplicated.align(xmap, delta.score = 3, bestRef = T)
  xmap.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3/Alignments/HybridM1/BSSSI/a.lines_BSSSI_0kb_0labels_cut_1_fullAln.xmap"
  xmap2 <- as.data.table(readxmap(xmap.file2))
  xmap2 <- remove.duplicated.align(xmap2, delta.score = 3, bestRef = T)
  exportXmap(xmap, "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3/Sandwich2_2/test_1_fullAln.xmap")
  exportXmap(xmap2, "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3/Sandwich2_2/test_2_fullAln.xmap")
}

get.chimq.score <- function(conflict.regions, cmap.with.chim1, cmap.with.chim2){
  cmap.with.chim1 <- as.data.table(cmap.with.chim1)
  cmap.with.chim2 <- as.data.table(cmap.with.chim2)
  setkey(cmap.with.chim1, 'CMapId')
  setkey(cmap.with.chim2, 'CMapId')
  conflict.regions <- flip.columns(conflict.regions, 'ConflictStartPos2', 'ConflictEndPos2', `>`)
  buffer <- 10000
  chim.q.scores <- unlist(lapply(1:nrow(conflict.regions),
                          function(i){cmap.with.chim1[conflict.regions$RefcontigID.x[[i]]][Position + buffer >= conflict.regions$ConflictStartPos[[i]] & Position - buffer <= conflict.regions$ConflictEndPos[[i]], min(ChimQuality)]}))
  chim.q.scores2 <- unlist(lapply(1:nrow(conflict.regions),
                          function(i){cmap.with.chim2[conflict.regions$RefcontigID.y[[i]]][Position + buffer >= conflict.regions$ConflictStartPos2[[i]] & Position - buffer <= conflict.regions$ConflictEndPos2[[i]], min(ChimQuality)]}))
  data.table(chim1=chim.q.scores, chim2=chim.q.scores2)
}

test.conflict.resolve <- function(){
  hybrid.dir1 <- "/home/users7/jwang/HybridScaffold/HybridScaffodTwo/TestDebugKeeSan_C076_job612/YingFeng/output/CACGAG"
  hybrid.dir2 <- "/home/users7/jwang/HybridScaffold/HybridScaffodTwo/TestDebugKeeSan_C076_job612/YingFeng/output/GCTCTTC"
  cmap.with.chim1 <- NULL
  cmap.with.chim2 <- NULL
  #hybrid.dir1 <- NULL
  #hybrid.dir2 <- NULL
  alignmerged.rdata.file <- "/home/users7/jwang/HybridScaffold/HybridScaffodTwo/TestDebugKeeSan_C076_job612/YingFeng/output/Sandwich1/test_indata.RData"
  pos.offset.rdata.file <- "/home/users7/jwang/HybridScaffold/HybridScaffodTwo/TestDebugKeeSan_C076_job612/YingFeng/output/Sandwich1/test_pos_offsets.RData"
  out.file1 <- "/home/users7/jwang/HybridScaffold/HybridScaffodTwo/TestDebugKeeSan_C076_job612/YingFeng/output/CombinedConflictsCut/test_BSPQI_merged.txt"
  out.file2 <- "/home/users7/jwang/HybridScaffold/HybridScaffodTwo/TestDebugKeeSan_C076_job612/YingFeng/output/CombinedConflictsCut/test_BSSSI_merged.txt"
  r <- run.conflict.resolve(hybrid.dir1, hybrid.dir2, alignmerged.rdata.file, pos.offset.rdata.file)
}

#run conflict resolution for sandwich assemblies
run.conflict.resolve <- function(hybrid.dir1, hybrid.dir2, align.rdata.file, pos.offset.rdata.file,
                                 out.file1=NULL, out.file2=NULL, type='qry', cmap.with.chim1=NULL, cmap.with.chim2=NULL){
  #Merge cut-status file from hybrid scaffold conflict resolutions
  conflicts.merged1 <- data.table()
  conflicts.merged2 <- data.table()
  if(!is.null(hybrid.dir1) && !is.null(hybrid.dir2)){
    conflicts.merged <- merge.conflict.file(hybrid.dir1, hybrid.dir2)
    conflicts.merged.stats <- conflict.intersect(hybrid.dir1, hybrid.dir2)
    conflicts.merged1 <- data.frame(lapply(conflicts.merged[[1]], as.character), stringsAsFactors=FALSE)
    conflicts.merged2 <- data.frame(lapply(conflicts.merged[[2]], as.character), stringsAsFactors=FALSE)
    #clean up header to restore the proper header column for cut_conflict file, data.table do not allow duplicate colnames
    conflicts.merged1 <- cleanup.header(conflicts.merged1)
    conflicts.merged2 <- cleanup.header(conflicts.merged2)
  }

  #detect and resolve conflicts from the two single-enzyme map
  if(!is.null(align.rdata.file) && !is.null(pos.offset.rdata.file)){
    load(align.rdata.file)
    load(pos.offset.rdata.file)
    ngs.cut.ids <- rbind(subset(conflicts.merged1, ref_leftBkpt_toCut == 'cut' | ref_rightBkpt_toCut == 'cut', 'refId'),
                          subset(conflicts.merged2, ref_leftBkpt_toCut == 'cut' | ref_rightBkpt_toCut == 'cut', 'refId'))
    ngs.cut.ids <- unique(ngs.cut.ids$refId)
    alignment.rdata <- subset(indata, !(QryContigID %in% ngs.cut.ids))
    pos.offset.rdata <- pos_offsets
    conflict.regions <- find.conflict.alignment(alignment.rdata, pos.offset.rdata)
    print(conflict.regions)
    conflict.tables <- generate.conflict.table(conflict.regions[[1]], type)
    conflict.tables2 <- generate.conflict.table(conflict.regions[[2]], type)

    conflict.tables[[1]] <- rbind(conflict.tables[[1]], conflict.tables2[[1]])
    conflict.tables[[2]] <- rbind(conflict.tables[[2]], conflict.tables2[[2]])

    #clean up header to restore the proper header column for cut_conflict file, data.table do not allow duplicate colnames
    conflict.tables[[1]] <- cleanup.header(conflict.tables[[1]])
    conflict.tables[[2]] <- cleanup.header(conflict.tables[[2]])
    conflicts.merged1 <- rbind(conflicts.merged1, conflict.tables[[1]])
    conflicts.merged2 <- rbind(conflicts.merged2, conflict.tables[[2]])

    #converting all coordinate to integer, sub-bp resolution is not really ideal
    conflicts.merged1 <- toIntegerCoord(conflicts.merged1)
    conflicts.merged2 <- toIntegerCoord(conflicts.merged2)
  }

  if(!is.null(out.file1) && !is.null(out.file2)){
    write.conflict.file(conflicts.merged1, out.file1)
    write.conflict.file(conflicts.merged2, out.file2)
  }

  #For chim-score guided conflict resolution
  #cmap.with.chim1 <-readcmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies/ChimQScore/BSPQI/BSPQI_REFINEFINAL1_with_Chim_Score.cmap")
  #cmap.with.chim2 <- readcmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies/ChimQScore/BSSSI/BSSSI_REFINEFINAL1_with_Chim_Score.cmap")
  if(!is.null(cmap.with.chim1) && !is.null(cmap.with.chim2)){
    cmap.with.chim1 <- readcmap(cmap.with.chim1)
    cmap.with.chim2 <- readcmap(cmap.with.chim2)
    chim.scores1 <- get.chimq.score(conflict.regions[[1]], cmap.with.chim1, cmap.with.chim2)
    chim.scores2 <- get.chimq.score(conflict.regions[[2]], cmap.with.chim1, cmap.with.chim2)
  }

  return(list(conflicts.merged1, conflicts.merged2, conflicts.merged.stats))
}

#this function clean up column header names
#note this only remove '.1' from colnames
#and return a data.frame such that duplicated column names are allowed
cleanup.header <- function(dt){
  names <- colnames(dt)
  if(length(names) == 0){
    return(dt)
  }
  names.fixed <- gsub('\\.1', "", names)
  df <- as.data.frame(dt)
  colnames(df) <- names.fixed
  colnames(df)[[1]] <- '# xMapId'
  return(df)
}

toIntegerCoord <- function(conflict){
  conflict$leftRefBkpt <- as.integer(conflict$leftRefBkpt)
  conflict$rightRefBkpt <- as.integer(conflict$rightRefBkpt)
  conflict$leftQryBkpt <- as.integer(conflict$leftQryBkpt)
  conflict$rightQryBkpt <- as.integer(conflict$rightQryBkpt)
  return(conflict)
}


write.conflict.file <- function(conflict.table, conflict.file){
  sink(conflict.file)
  header2 <- "# id/-1\tref\tid/-1\tposition/-1\tposition/-1\t+/-\tokay/cut/-\tokay/cut/-\tokay/exclude/-\tqry\tid/-1\tposition/-1\tposition/-1\t+/-\tokay/cut/-\tokay/cut/-\tokay/exclude/-\n"
  cat(colnames(conflict.table), sep='\t')
  cat('\n')
  cat(header2)
  sink()
  write.table(conflict.table, conflict.file, append=TRUE, sep='\t', col.names = F, quote = F, row.names = F)
}


flip.xmap <- function(align, out.file = NULL, unique.ref=F){
  if(is.character(align) && file.exists(align)){
    #align <- readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunDefault/BSPQI/align1/align1.xmap")
    align <- readxmap(align)
  }
  if(nrow(align)==0){
    return(align)
  }

  align <- as.data.table(align)

  #a helper function that pass through all value in a column
  all.pass.filter <- function(col1, col2){
    return(rep(TRUE, length(col)))
  }

  orientation.filter <- function(col1, col2){
    return(align$Orientation == '-')
  }

  align.flip <- flip.columns(align, "QryContigID", "RefcontigID", test.fun = all.pass.filter )
  align.flip <- flip.columns(align.flip, "QryLen", "RefLen", test.fun = all.pass.filter )
  align.flip <- flip.columns(align.flip, "QryStartPos", "RefStartPos", test.fun = all.pass.filter )
  align.flip <- flip.columns(align.flip, "QryEndPos", "RefEndPos", test.fun = all.pass.filter )
  #make sure position is consistent with orientation
  align.flip <- flip.columns(align.flip, "RefStartPos", "RefEndPos", test.fun = orientation.filter)
  align.flip <- flip.columns(align.flip, "QryStartPos", "QryEndPos", test.fun = orientation.filter)

  align.flip$Alignment <- flip.align.str(align)

  if(unique.ref){
    keep <- align.flip[, Confidence == max(Confidence), by=QryContigID]$V1
    align.flip <- align.flip[keep]
  }

  if(!is.null(out.file)){
    exportXmap(align.flip, out.file)
    #old.name <- colnames(align.flip)[[1]]
    #new.name <- paste0("#h ", old.name)
    #setnames(align.flip, old.name, new.name)
    #write.table(align.flip, out.file,
    #            sep = "\t", row.names = F, quote = F)
    #setnames(align.flip, new.name, old.name)
  }
  return(align.flip)
}

flip.align.str <- function(align){
  if(nrow(align)==0){
    return(align$Alignment)
  }
  align <- processXmapAlignments(align)
  new.align.str <- unlist(lapply(1:nrow(align),
                                 function(i){
                                   if(align$Orientation[[i]]=='+'){
                                    paste("(", align$alignQryInd[[i]], ",", align$alignRefInd[[i]], ")", sep="", collapse="")
                                   }else{
                                     paste("(", rev(align$alignQryInd[[i]]), ",", rev(align$alignRefInd[[i]]), ")", sep="", collapse="")
                                   }
                              }))
  align$alignQryInd <- NULL
  align$alignRefInd <- NULL
  return(new.align.str)
}


test.flip.xmap <- function(){
  xmap1 <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridScaffold_NonHaplotype/ReRun_03242016/BSPQI_with_hash3/align1/align1.xmap"
  xmap1.orig <- read.bng.maps(xmap1)
  xmap1.flip <- flip.xmap(xmap1)
  #xmap1.flip <- remove.duplicated.align(xmap1.flip, delta.score = 5, bestRef = F)
  xmap1.out.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies_with_assignAlignT11_align1Hash3/test_1_fullAln.xmap"
  write.table(xmap1.flip, xmap1.out.file, sep='\t', quote = F, row.names = F, col.names = F)

  xmap2 <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridScaffold_NonHaplotype/ReRun_03242016/BSSSI_with_hash3/align1/align1.xmap"
  xmap2.flip <- flip.xmap(xmap2)
  #xmap2.flip <- remove.duplicated.align(xmap2.flip, delta.score = 5, bestRef = F)
  xmap2.out.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies_with_assignAlignT11_align1Hash3/test_2_fullAln.xmap"
  write.table(xmap2.flip, xmap2.out.file, sep='\t', quote = F, row.names = F, col.names = F)
}

#this operates on the BNG contig graph rather than the NGS graph
#it find contigs whose relative position on the BNG cluster that overlap with each other which signal
#a possible confict of the two-enzyme contigs, this acts as a proxy to evaluate two-enzyme scaffold quality
find.overlap.contigs <- function(){
  #load('~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies_withHybridM3/test_paths.RData')
  load('~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_new_pacbio/TestRunTwoPassNewAsm_two/two_enzyme_hybrid_scaffold_M1/Sandwich2/test_paths.RData')
  component <- outlist$paths[[1]]
  #compute overlap of contigs in an scaffolded component
  getoverlap <- function(component){
    path.coord <- as.data.table(component$path_coordinates)
    path.coord.reorder <- flip.rel.coord(path.coord)
    setkeyv(path.coord.reorder, c('node_type', 'rel_start', 'rel_end'))
    lastInd <- nrow(path.coord.reorder)
    path.coord.reorder$overlap <- path.coord.reorder[,c(0,rel_start[-1]-rel_end[-length(rel_end)]), by=node_type]$V1
    return(path.coord.reorder)
  }
  paths.with.overlap <- lapply(paths, getoverlap)
  paths.with.overlap <- lapply(1:length(paths.with.overlap),
                               function(i){p <- paths.with.overlap[[i]]; p$path_id <- i; p})
  paths.with.overlap <- rbindlist(paths.with.overlap)
}

addKeys <- function(keys, dt=NULL){
  if(is.null(dt) || nrow(dt)==0){
    dt <- data.table(Name=keys, Ind=1:length(keys))
  }else{
    not.found.inds <- which(is.na(getInd(keys, dt)))
    #nothing new to add
    if(length(not.found.inds)==0){
      return(dt)
    }
    #adding new keys to the ind table
    startInd <- max(dt$Ind)+1
    endInd <- startInd + length(not.found.inds)-1
    dt.new <- data.table(Name=keys[not.found.inds], Ind=startInd:endInd)
    dt <- rbind(dt, dt.new)
  }
  setkey(dt, Name)
  return(dt)
}

#this function create a graph that represents the relationship of NGS contigs from hybrid-scaffolding
get.graphfrom.agp <- function(agp.file, vertexTable=NULL, directed=FALSE){
  agp <- read.table(agp.file, header=T, stringsAsFactors = F, sep='\t', comment.char = "", skip=7)
  agp <- as.data.table(agp[grepl('Super', agp$X..Obj_Name),])
  agp.contigs <- agp[Compnt_Type=='W']
  scaffold.g <- make_empty_graph(directed = directed)
  #adding vertex
  vertexTable <- addKeys(agp.contigs$CompntId_GapLength, vertexTable)
  scaffold.g <- scaffold.g + agp.contigs$CompntId_GapLength
  #adding edges for each super-scaffold
  if(nrow(agp.contigs[Obj_Start==1]) != length(unique(agp.contigs$X..Obj_Name))){
    warning("Not all super contigs in agp start at position 1")
  }
  inds <- which(agp.contigs$Obj_Start == 1)
  #agp.contigs[,{scaffold.g <<- scaffold.g + path(CompntId_GapLength);1}, by=X..Obj_Name]
  for(i in 2:length(inds)){
    if(length(inds[i-1]:(inds[i]-1)) == 2){
      scaffold.g <- scaffold.g  + edge(agp.contigs$CompntId_GapLength[inds[i-1]], agp.contigs$CompntId_GapLength[inds[i]-1])
    }else{
      scaffold.g <- scaffold.g + path(agp.contigs$CompntId_GapLength[inds[i-1]:(inds[i]-1)])
    }
    d <- degree(scaffold.g)
    #if(d['flattened_line_15286'] > 1){
    #  browser()
    #}
    #print(paste0("degree is ", d['flattened_line_15286']))
  }
  return(scaffold.g)
}


overlap.sv.with.conflicts <- function(conflict.file, smap){
  conflict.file <- ""
  smap.file <- ""
  conflicts <- read.hybrid.conflictfile(conflict.file)
  sv <- read.bng.maps(smap.file)
}


getInd <- function(key, dt){
  return(dt[key, Ind])
}




find.gap.conflict <- function(){

}








