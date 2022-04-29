#This script evaluate the quality or accuracy of hybrid-scaffold
library(BNGBase)

#check accuracy of NGS anchoring using two-color
checkNGSAlign <- function(){
  #check NGS anchoring order accuracy from scaffolding methods
  ref1.file <- '~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/NGSAlignValidate/AlignToRefTwoColor/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash1_errEst.xmap'
  ref2.file <- '~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/NGSAlignValidate/AlignToRefTwoColor/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash2_errEst.xmap'

  ref1.file <- '~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestTwoColor/OneStep_alignFinalParams/Twocolor/Hash1/a.lines_BSPQI_BSSSI_0kb_0labels_twocolor_errEst.xmap'
  ref2.file <- '~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestTwoColor/OneStep_alignFinalParams/Twocolor/Hash2/a.lines_BSPQI_BSSSI_0kb_0labels_twocolor_errEst.xmap'

  two.color.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoV2Asm/alignfinal/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash1_errEst.xmap"
  two.color.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoV2Asm/alignfinal/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash2_errEst.xmap"
  single.color.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSPQI_0kb_0labels_1_errEst.xmap"
  single.color.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSSSI_0kb_0labels_1_errEst.xmap"


  two.color.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash1_errEst.xmap"
  two.color.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash2_errEst.xmap"
  single.color.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSPQI_0kb_0labels_1_errEst.xmap"
  single.color.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSSSI_0kb_0labels_1_errEst.xmap"

  single.align.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSPQI/align_final_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"
  single.align.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSSSI/align_final_M1/BSSSI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"

  ref <- get.merged.twocolor.align(ref1.file, ref2.file)

  align.final.twocolor <- get.merged.twocolor.align(two.color.file, two.color.file2)
  align.single.color <- as.data.table(readxmap(single.color.file))
  align.single.color2 <- as.data.table(readxmap(single.color.file2))

  align.single <- as.data.table(readxmap(single.align.file))
  align.single2 <- as.data.table(readxmap(single.align.file2))

   aligns <- getAlignFinalStats(single.color.file, single.color.file2, two.color.file , two.color.file2, "", "", 10, 10)


   ngs.align.out <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/NGSAlignValidate/NGS_cutted_align.xmap"
   ngs.name.table <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/fa2cmap2color/a.lines.fasta.cutted_BSPQI_0kb_0labels_key.txt"

   ngs.ref.align <- get.ngs.ref.align(ngs.align.out, ngs.name.table)
   ngs.ref.align2 <- get.ngs.ref.align("", "")
   #order.dist <- evalScaffoldOrder(ref[Confidence > 8], ngs.ref.align[Confidence > 11])
   order.dist <- evalScaffoldOrder(align.final.twocolor[Confidence > 12 & QryStartPos != QryEndPos], ref[Confidence > 9])
   order.dist <- evalScaffoldOrder(align.final.twocolor[Confidence > 12], ngs.ref.align[Confidence >11])

   order.dist.combined <- evalScaffoldOrder(aligns$ExportAlign, ngs.ref.align[Confidence > 11])
   order.dist.twocolor <- evalScaffoldOrder(aligns$TwoColorAlign, ngs.ref.align[Confidence >11])
   order.dist.singlecolor <- evalScaffoldOrder(aligns$SingleHybridAlign, ngs.ref.align[Confidence >11])

   order.dist <- evalScaffoldOrder(align.single.color[Confidence > 12], ngs.ref.align[Confidence >11])
   order.dist <- evalScaffoldOrder(align.single.color2[Confidence > 12], ngs.ref.align[Confidence >11])
   order.dist <- evalScaffoldOrder(align.single[Confidence > 1], ngs.ref.align2[Confidence >11])
   order.dist <- evalScaffoldOrder(align.single2[Confidence > 1], ngs.ref.align2[Confidence >11])

   ids <- setdiff(order.dist.combined[abs(dist) > 1, V1], order.dist[abs(dist) > 1, V1])
   order.dist.add <- order.dist.combined[V1 %in% ids]

}


#evaluate scaffold quality for TGH base on a sequence alignment of input NGS to reference as
#the ground truth
eval.twocolor.scaffold.accuracy <- function(scaffol.out, ngs.align.out){
  align.final.dir <- "/alignfinal"
  two.color.file <- dir(path = paste0(scaffold.out, align.final.dir), pattern = '*hash1.*.xmap', full.names = T)
  two.color.file2 <- dir(path = paste0(scaffold.out, align.final.dir), pattern = '*hash2.*.xmap', full.names = T)
  single.color.file <- dir(path = paste0(scaffold.out,  align.final.dir), pattern = '*BSPQI_0kb_0labels_1_errEst.xmap', full.names = T)
  single.color.file2 <- dir(path = paste0(scaffold.out, align.final.dir), pattern = '*cut.*_BSSSI_0kb_0labels_1_errEst.xmap', full.names = T)


  single.align.file <- dir(path = paste0(scaffold.out, "/BSPQI/align_final_M1/"), pattern = '*NGS.*HYBRID.*.xmap', full.names = T)
  single.align.file2 <- dir(path = paste0(scaffold.out, "/BSSSI/align_final_M1/"), pattern = '*NGS.*HYBRID.*.xmap', full.names = T)

  aligns <- getAlignFinalStats(single.color.file, single.color.file2, two.color.file , two.color.file2, "", "", Confidence.T = 10, Confidence2.T = 12)
  align.single1 <- remove.duplicated.align(as.data.table(readxmap(single.align.file[[1]])), 0.1, 0.00001, bestRef = T)
  align.single2 <- remove.duplicated.align(as.data.table(readxmap(single.align.file2[[1]])), 0.1, 0.000001, bestRef = T)


  #ngs.align.out <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/NGS_cutted_out.mum.mgap.xmap"
  ngs.name.table <- dir(path = paste0(scaffold.out, "/fa2cmap2color/"), pattern = '*_BSPQI_BSSSI_0kb_0labels_key*', full.names = T)
  ngs.name.table.single <- dir(path = paste0(scaffold.out, "/BSPQI/fa2cmap/"), pattern = '*0kb_0labels_key*', full.names = T)
  #ngs.name.table <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/SOAP_denovo/soap.result.scafSeq.fa.cutted.fasta.mum_ReRun_BSPQI_BSSSI_key_fixednames.txt"
  #ngs.name.table.single <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/SOAP_denovo/soap.result.scafSeq.fa.cutted.fasta.mum_ReRun_BSPQI_key_fixednames.txt"

  ngs.align.out <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/NGS_cutted_align.xmap"
  ngs.name.table <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/fa2cmap2color/a.lines.fasta.cutted_BSPQI_0kb_0labels_key.txt"
  ngs.name.table.single <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSPQI/agp_fasta_M1/a.lines_BSPQI_0kb_0labels_key.txt.cutted.txt"


  scaffold.out <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunScaffWithSOAP/ReRun10042016/"
  two.color.file <- dir(path = paste0(scaffold.out, "alignfinal"), pattern = '*hash1.*.xmap', full.names = T)
  two.color.file2 <- dir(path = paste0(scaffold.out, "alignfinal"), pattern = '*hash2.*.xmap', full.names = T)
  single.color.file <- dir(path = paste0(scaffold.out, "alignfinal"), pattern = '*BSPQI_0kb_0labels_1_errEst.xmap', full.names = T)
  single.color.file2 <- dir(path = paste0(scaffold.out, "alignfinal"), pattern = '*cut_BSSSI_0kb_0labels_1_errEst.xmap', full.names = T)

  single.align.file <- dir(path = paste0(scaffold.out, "/BSPQI/align_final_M1/"), pattern = '*NGS.*HYBRID.*.xmap', full.names = T)
  single.align.file2 <- dir(path = paste0(scaffold.out, "/BSSSI/align_final_M1/"), pattern = '*NGS.*HYBRID.*.xmap', full.names = T)


  aligns <- getAlignFinalStats(single.color.file, single.color.file2, two.color.file , two.color.file2, "", "", Confidence.T = 9, Confidence2.T = 11)
  align.single1 <- remove.duplicated.align(as.data.table(readxmap(single.align.file[[1]])), 0.1, 0.00001, bestRef = T)
  align.single2 <- remove.duplicated.align(as.data.table(readxmap(single.align.file2[[1]])), 0.1, 0.000001, bestRef = T)

  ngs.align.out <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/SOAP_denovo/soap.result.scafSeq.fa.cutted.fixednames.fasta.mum.mgap.xmap"
  ngs.name.table <- dir(path = paste0(scaffold.out, "/fa2cmap2color/"), pattern = '*_BSPQI_BSSSI_0kb_0labels_key*', full.names = T)
  ngs.name.table.single <- dir(path = paste0(scaffold.out, "/BSPQI/fa2cmap/"), pattern = '*0kb_0labels_key*', full.names = T)
  ngs.name.table <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/SOAP_denovo/soap.result.scafSeq.fa.cutted.fasta.mum_ReRun_BSPQI_BSSSI_key_fixednames.txt"
  ngs.name.table.single <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/SOAP_denovo/soap.result.scafSeq.fa.cutted.fasta.mum_ReRun_BSPQI_key_fixednames.txt"

  ngs.ref.align <- get.ngs.ref.align(ngs.align.out, ngs.name.table)
  ngs.ref.align2 <- get.ngs.ref.align(ngs.align.out, ngs.name.table.single)

  print(paste0("Scaffolding accuracy for TGH alignfinal: "))
  order.dist.combined <- evalScaffoldOrder(aligns$ExportAlign[!(grepl('100000', RefcontigID) | grepl('200000', RefcontigID))], ngs.ref.align[Confidence > 11])
  order.dist.twocolor <- evalScaffoldOrder(aligns$TwoColorAlign, ngs.ref.align[Confidence >11])
  order.dist.singlecolor <- evalScaffoldOrder(aligns$SingleHybridAlign, ngs.ref.align[Confidence >11])
  order.dist.single <- evalScaffoldOrder(aligns$SingleOnlyAlign, ngs.ref.align[Confidence >11])

  print(paste0("Scaffolding accuracy for single-enzyme hybrid alignfinal: "))
  order.dist.single1 <- evalScaffoldOrder(align.single1, ngs.ref.align2[Confidence > 11])
  order.dist.single2 <- evalScaffoldOrder(align.single2, ngs.ref.align2[Confidence > 11])

  print("Gap sizing errors for TGH: ")
  gap.dist.combined <- evalScaffoldDistance(aligns$ExportAlign, ngs.ref.align[Confidence > 11])
  print("Gap sizing errors for single-enzyme Hybrid1 (BSPQI): ")
  gap.dist.single1 <- evalScaffoldDistance(align.single1, ngs.ref.align2[Confidence > 11])
  print("Gap sizing errors for single-enzyme Hybrid2 (BSSSI): ")
  gap.dist.single2 <- evalScaffoldDistance(align.single2, ngs.ref.align2[Confidence > 11])

  #analyzeCV(gap.dist.combined)
  #analyzeCV(gap.dist.single1)
  return(list(order.dist.combined, ngs.ref.align, order.dist.single1, order.dist.single2, ngs.ref.align2))
}

eval.twocolor.scaffold.accuracy <- function(scaffol.out, ngs.align.out){
  #browser()
  align.final.dir <- "/TGH_M1/alignfinal"
  two.color.file <- dir(path = paste0(scaffold.out, align.final.dir), pattern = '*hash1.*.xmap', full.names = T)
  two.color.file2 <- dir(path = paste0(scaffold.out, align.final.dir), pattern = '*hash2.*.xmap', full.names = T)
  single.color.file <- dir(path = paste0(scaffold.out,  align.final.dir), pattern = 'E_BSPQI_Q_NGScontigs_A_HYBRID.xmap', full.names = T)
  single.color.file2 <- dir(path = paste0(scaffold.out, align.final.dir), pattern = 'E_CTTAAG_Q_NGScontigs_A_HYBRID.xmap', full.names = T)

  single.align.file <- dir(path = paste0(scaffold.out, "/BSPQI/align_final_M1/"), pattern = '*NGS.*HYBRID.*.xmap', full.names = T)
  single.align.file2 <- dir(path = paste0(scaffold.out, "/CTTAAG/align_final_M1/"), pattern = '*NGS.*HYBRID.*.xmap', full.names = T)

  aligns <- getAlignFinalStats(single.color.file, single.color.file2, two.color.file , two.color.file2, "", "", Confidence.T = 10, Confidence2.T = 12)
  align.single1 <- remove.duplicated.align(as.data.table(readxmap(single.align.file[[1]])), 0.1, 0.00001, bestRef = T)
  align.single2 <- remove.duplicated.align(as.data.table(readxmap(single.align.file2[[1]])), 0.1, 0.000001, bestRef = T)

  #ngs.align.out <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/NGS_cutted_out.mum.mgap.xmap"
  ngs.name.table <- dir(path = paste0(scaffold.out, "/TGH_M1/fa2cmap/"), pattern = '*_BSPQI_CTTAAG_0kb_0labels_key*', full.names = T)
  ngs.name.table.single <- dir(path = paste0(scaffold.out, "/BSPQI/fa2cmap/"), pattern = '*0kb_0labels_key*', full.names = T)
  #ngs.name.table <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/SOAP_denovo/soap.result.scafSeq.fa.cutted.fasta.mum_ReRun_BSPQI_BSSSI_key_fixednames.txt"
  #ngs.name.table.single <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/SOAP_denovo/soap.result.scafSeq.fa.cutted.fasta.mum_ReRun_BSPQI_key_fixednames.txt"

  ngs.ref.align <- get.ngs.ref.align(ngs.align.out, ngs.name.table)
  ngs.ref.align2 <- get.ngs.ref.align(ngs.align.out, ngs.name.table.single)

  print(paste0("Scaffolding accuracy for TGH alignfinal: "))
  order.dist.combined <- evalScaffoldOrder(aligns$ExportAlign[!(grepl('100000', RefcontigID) | grepl('200000', RefcontigID))], ngs.ref.align[Confidence > 11])
  order.dist.twocolor <- evalScaffoldOrder(aligns$TwoColorAlign, ngs.ref.align[Confidence >11])
  order.dist.singlecolor <- evalScaffoldOrder(aligns$SingleHybridAlign, ngs.ref.align[Confidence >11])
  order.dist.single <- evalScaffoldOrder(aligns$SingleOnlyAlign, ngs.ref.align[Confidence >11])
  #browser()
  print(paste0("Scaffolding accuracy for single-enzyme hybrid alignfinal: "))
  order.dist.single1 <- evalScaffoldOrder(align.single1, ngs.ref.align2[Confidence > 11])
  order.dist.single2 <- evalScaffoldOrder(align.single2, ngs.ref.align2[Confidence > 11])

  print("Gap sizing errors for TGH: ")
  gap.dist.combined <- evalScaffoldDistance(aligns$ExportAlign, ngs.ref.align[Confidence > 11])
  print("Gap sizing errors for single-enzyme Hybrid1 (BSPQI): ")
  gap.dist.single1 <- evalScaffoldDistance(align.single1, ngs.ref.align2[Confidence > 11])
  print("Gap sizing errors for single-enzyme Hybrid2 (BSSSI): ")
  gap.dist.single2 <- evalScaffoldDistance(align.single2, ngs.ref.align2[Confidence > 11])

  #analyzeCV(gap.dist.combined)
  #analyzeCV(gap.dist.single1)
  return(list(order.dist.combined, ngs.ref.align, order.dist.single1, order.dist.single2, ngs.ref.align2))
}



analyzeCV <- function(gap.dist){
  breaks <- c(seq(0, 5000, 500), seq(5000, 50000, 5000), 10e6)
  gap.dist[, Inds := findInterval(abs(ref.dist), breaks)]
  gap.dist[, PercentDiff :=  difference/(abs(aligndist)+0.1)*100]
  g <- ggplot(gap.dist)
  g <- g + geom_boxplot(aes(x=as.factor(Inds), y=PercentDiff))
  g <- g + coord_cartesian(ylim=c(-20, 20))
  g <- g + scale_x_discrete(lab=breaks)
  g <- g + xlab('Reference distance')
  plot(g)
  return(gap.dist)
}


eval.single.scaffold.accuracy <- function(scaffold.out, ngs.align.out){
  #single.align.file <- dir(path = paste0(scaffold.out, "/align_final/"), pattern = '*NGS.*HYBRID.*.xmap', full.names = T)
  single.align.file <- dir(path = paste0(scaffold.out, "/align_final/"), pattern="EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD_with_T5\\.xmap", full.names=T) #pattern = '*NGS.*HYBRID.*.xmap', full.names = T)
  ngs.name.table.single <- dir(path = paste0(scaffold.out, "/fa2cmap/"), pattern = '*0kb_0labels_key*', full.names = T)
  #single.align.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/HybridRelaxParams/BSPQI/align_final/EXP_REFINEFINAL1_bppAdjust_cmap_mar3_NA12878_scf_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"
  #single.align.file2 <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/HybridRelaxParams//BSSSI/align_final/EXP_REFINEFINAL1_bppAdjust_cmap_mar3_NA12878_scf_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"
  #ngs.align.out <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/PacBio_mar/mar3_NA12878.scf.fasta.cutted.fasta.mum.mgap.xmap"
  #ngs.name.table <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/HybridRelaxParams/BSPQI/agp_fasta/mar3_NA12878.scf_BSPQI_0kb_0labels_key.txt.cutted.txt"
  #browser()
  if(is.null(single.align.file) | is.null(ngs.name.table.single)){
    stop("Cannot find input file")
  }
  ngs.ref.align.single <- get.ngs.ref.align(ngs.align.out, ngs.name.table.single)
  if(nrow(ngs.ref.align.single) < 1){
    stop("Check ngs align or name file, no ngs alignments was parsed")
  }
  align.final1 <- as.data.table(readxmap(single.align.file[[1]]))
  order.dist.single1 <- evalScaffoldOrder(align.final1, ngs.ref.align.single)
  print("Gap sizing error in scaffold: ")
  gap.dist.single1 <-  evalScaffoldDistance(align.final1, ngs.ref.align.single)
}

eval.align.final.strategy <- function(){
  ngs.align.out <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/NGS_cutted_out.mum.mgap.xmap"
  scaffold.out <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/BSPQI/"

  ngs.align.out <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/GIAB_DISCOVAR/GIAB_HG004_discovar_contig.fa.mum.mgap.xmap"
  scaffold.out <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_Illmn2/BSPQI/"

  ngs.name.table.single <- dir(path = paste0(scaffold.out, "/fa2cmap/"), pattern = '*0kb_0labels_key*', full.names = T)
  #browser()
  eval.align.final.strategy.helper <- function(single.align.file, ngs.name.table.single, flip=F){
    if(is.null(single.align.file) | is.null(ngs.name.table.single)){
      stop("Cannot find input file")
    }
    ngs.ref.align.single <- get.ngs.ref.align(ngs.align.out, ngs.name.table.single)
    if(nrow(ngs.ref.align.single) < 1){
      stop("Check ngs align or name file, no ngs alignments was parsed")
    }
    align.final1 <- as.data.table(readxmap(single.align.file[[1]]))
    if(flip){
      align.final1 <- get_flip_align_final_simple(align.final1)
    }
    align.final1[,isBest:=Confidence==max(Confidence), by=QryContigID]
    align.final <- align.final1[isBest==T]
    scores <- 5:15
    ret <- rbindlist(lapply(scores, function(T){
      aligns <- align.final1[align.final1$Confidence > T,]
      order.dist.single1 <- evalScaffoldOrder(aligns, ngs.ref.align.single)
      list(T=T, scaffold_order_accuracy=sum(abs(order.dist.single1$dist)==1)/nrow(order.dist.single1), Total_aligned=nrow(aligns))
    }))
    return(ret)
  }

  single.align.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/E_CTTAAG_Q_NGScontigs_A_HYBRID.xmap"
  r1 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single)
  r1$Align <- "Default"

  single.align.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/E_CTTAAG_Q_NGScontigs_A_HYBRID_test_xmapfilter.xmap"
  r2 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single, flip = T)
  r2$Align <- "Greedy-alg-1color"

  single.align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/align_final_combined.xmap"
  r3 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single)
  r3$Align <- "Greedy-alg-2color"

  single.align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/a.lines.fasta.cut_BSPQI_CTTAAG_0kb_0labels_NGS_contigs_HYBRID_Export.xmap"
  r4 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single)
  r4$Align <- "Default-2color"

  single.align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/TGH_M1/alignfinal/align_final_ref_hg19_combined.xmap"
  r5 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single)
  r5$Align <- "hg19-2color"

  #Illumina 2 dataset
  single.align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_Illmn2/CTTAAG/align_final_M1/EXP_REFINEFINAL1_bppAdjust_cmap_GIAB_HG004_discovar_contig_fa_NGScontigs_HYBRID_SCAFFOLD_AllNGS_T5.xmap"
  r1 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single)
  r1$Align <- "Default"

  single.align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_Illmn2/TGH_M1/alignfinal/CTTAAG_align_final_bestRef.xmap"
  r2 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single, flip=T)
  r2$Align <- "BestRef_reverse"

  single.align.file <-"~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_Illmn2/TGH_M1/alignfinal/CTTAAG_align_final_grdalg.xmap"
  r3 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single, flip = F)
  r3$Align <- "Greedy-alg-1color1"

  single.align.file <-"~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_Illmn2/TGH_M1/alignfinal/CTTAAG_align_final_grdalg1.xmap"
  r4 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single, flip = F)
  r4$Align <- "Greedy-alg-1color"

  #Illumina 3 dataset
  single.align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_Illmn3/CTTAAG/align_final_M1/EXP_REFINEFINAL1_bppAdjust_cmap_GIAB_HG004_abyss1_0_sealer_contig_fa_NGScontigs_HYBRID_SCAFFOLD.xmap"
  r1 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single)
  r1$Align <- "Default"

  single.align.file <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_Illmn3/CTTAAG/align_final_M1/EXP_REFINEFINAL1_bppAdjust_cmap_GIAB_HG004_abyss1_0_sealer_contig_fa_NGScontigs_HYBRID_SCAFFOLD.xmap"
  r2 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single)
  r2$Align <- "All-NGS"

  single.align.file <-"~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI_Illmn3/CTTAAG/hybrid_scaffolds_M1/test_xmapFilter.xmap"
  r3 <- eval.align.final.strategy.helper(single.align.file, ngs.name.table.single, flip = T)
  r3$Align <- "Greedy-alg"

  g <- ggplot(rbind(r1,r2,r3, r4, r5))
  g <- g  + geom_line(aes(x=Total_aligned, y=scaffold_order_accuracy, color=Align))
  g <- g  + geom_point(aes(x=Total_aligned, y=scaffold_order_accuracy, color=Align))
  g <- g + scale_y_continuous(limit=c(0.5,1)) + scale_x_continuous(limit=c(7000,16500))
  g <- g + scale_x_continuous(breaks = seq(7000, 21500, 1000))
  g
}

get_flip_align_final_simple <- function(align_final){
  old_cols <- names(align_final)

  new_cols <- c("XmapEntryID",    "RefcontigID",  "QryContigID",   "RefStartPos",   "RefEndPos",     "QryStartPos",   "QryEndPos",     "Orientation",   "Confidence",    "HitEnum",
                "RefLen",         "QryLen",       "LabelChannel",  "Alignment")

  setnames(align_final, old_cols, new_cols)
  align_final <- align_final[,old_cols, with=F]
  return(align_final)
}



run.scaffold.eval <- function(){
  #human evaluation
  ngs.align.out <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/TestExtendSplit/a.lines.fasta.cut.fasta.mum.mgap.xmap"
  scaffold.out <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunExtenSplit/WithUpdatedRefAligner/Illumina_Discovar/"

  single.align.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSPQI/align_final_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"
  single.align.file2 <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSSSI/align_final_M1/BSSSI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"


  #BNG2 human NA12878
  ngs.align.out <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/TestExtendSplit/a.lines.fasta.cut.fasta.mum.mgap.xmap"
  scaffold.out <- "~/RemoteServer/home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19_with_BSPQI/"

  single.align.file <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSPQI/align_final_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"
  single.align.file2 <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSSSI/align_final_M1/BSSSI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"


  #SBJ evaluation
  ngs.align.out <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/SBJ_PB/SBJ_40x_p_contigsGr40kb.fa.cut.fasta.mum.mgap.xmap"
  scaffold.out <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/ReRun/"

  single.align.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSPQI/align_final_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"
  single.align.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSSSI/align_final_M1/BSSSI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"

  #Bird
  ngs.align.out <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/HummingBird/seq.fa.cut.fasta.mum.mgap.xmap"
  scaffold.out <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/PB_Andy/ReRun_01092017/"

  single.align.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSPQI/align_final_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"
  single.align.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSSSI/align_final_M1/BSSSI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"

  #SoyBean
  ngs.align.out <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/SoyBean/Soybean_Falcon_seed_5.5_60_60_1_Assembly_POSTQUIVER.fasta.cut.fasta.mum.mgap.xmap"
  scaffold.out <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSoyBean/ReRun12122016/"

  single.align.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSPQI/align_final_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"
  single.align.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSSSI/align_final_M1/BSSSI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"


  #ngs.align.out <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/  "
  #ngs.name.table <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/fa2cmap2color/a.lines.fasta.cutted_BSPQI_0kb_0labels_key.txt"
  #ngs.name.table.single <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/BSPQI/agp_fasta_M1/a.lines_BSPQI_0kb_0labels_key.txt.cutted.txt"

  eval.stats <- eval.twocolor.scaffold.accuracy(scaffold.out, ngs.align.out)

  scaffold.out.single.discovar <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/SingleEnzyme/PB30x/"
  scaffold.out.single.discovar <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridAutoTest/TestRun/NA12878/Illumina_Discovar/BSPQI/"
  scaffold.out.single.discovar <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunBNG2/Sample19/CTTAAG/"
  eval.single.scaffold.accuracy(scaffold.out.single.discovar, ngs.align.out)

}

analyze.align.final <- function(){
  align1 <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/ReRun10032016/alignfinal/mar3_NA12878.scf.fasta.cutted_BSPQI_0kb_0labels_1_errEst.xmap"))
  align1.all <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/ReRun10032016/alignfinal/mar3_NA12878.scf.fasta.cutted_BSPQI_0kb_0labels_noBestRef.xmap"))
  align2 <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/ReRun10032016/alignfinal/mar3_NA12878.scf.fasta.cutted_BSSSI_0kb_0labels_1_errEst.xmap"))
  align.single1 <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/ReRun10032016/BSPQI/align_final_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_mar3_NA12878_scf_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"))
  align.single1.rerun <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunAutoWithPacBio/ReRun10032016/BSPQI/align_final_M1/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_mar3_NA12878_scf_fasta_NGScontigs_HYBRID_SCAFFOLD_RERUN2.xmap"))
  align.single1 <- remove.duplicated.align(align.single1, 0.1, 0.00001, bestRef = T)

  ids.overlap <- intersect(align1$QryContigID, align2$QryContigID)
  align.intersect <- merge(align1, align2, by='QryContigID')
  conflict.ids <- align.intersect[RefcontigID.x != RefcontigID.y, QryContigID]
  hist(align.intersect[,((RefEndPos.x - RefStartPos.x) - (RefEndPos.y - RefStartPos.y))/2], 1000)
  conflict.ids <- align.intersect[abs((RefEndPos.x - RefStartPos.x) - (RefEndPos.y - RefStartPos.y))/2 > 40000, QryContigID]
  order.dist1 <- evalScaffoldOrder(align1[Confidence > 10], eval.stats[[2]])
  order.dist2 <- evalScaffoldOrder(align2[Confidence > 10], eval.stats[[2]])
  order.dist.single1 <- evalScaffoldOrder(align.single1.rerun[Confidence > 10], eval.stats[[5]])

  order.dist1.prescaff <- evalScaffoldOrder(align1[Confidence > 10 & QryContigID %in% align.single1$QryContigID], eval.stats[[2]])

  order.dist.single1 <- evalScaffoldOrder(align1[Confidence > 10 & QryContigID %in% ids.overlap], eval.stats[[2]])
  order.dist.single1 <- evalScaffoldOrder(align1[Confidence > 10 & !(QryContigID %in% conflict.ids)], eval.stats[[2]])

  length(intersect(align1$QryContigID, align.single1$QryContigID))
  align.intersect2 <- merge(align1, align.single1, 'QryContigID')

  order.dist1[,pair:=paste0(V1, "_", V2)]
  order.dist2[,pair:=paste0(V1, "_", V2)]
  order.dist.single1[,pair:=paste0(V1, "_", V2)]
  order.dist.single1[,pair2:=paste0(V2, "_", V1)]
  order.dist.single2[,pair:=paste0(V1, "_", V2)]

  order.dist.merge <- merge(order.dist1, order.dist2, 'pair')
  order.dist.merge2 <- merge(order.dist1, order.dist.single1, 'pair')

}




#this function read the mummer alignment in table format and make it into a xmap format
get.ngs.ref.align <- function(ngs.align.out, ngs.name.table){
  #browser()
  ngs.ref.align <- read.table(ngs.align.out, sep='\t', stringsAsFactors = F)
  setnames(ngs.ref.align, c("XmapEntryID", "RefcontigID", "QryContigID",  "RefStartPos", "QryStartPos", "RefEndPos",  "QryEndPos", "Orientation",  "Confidence",   "HitEnum",  "QryLen"       ,"RefLen", "LabelChannel"))
  ngs.ref.align <- as.data.table(ngs.ref.align)
  ngs.name.table <- read.table(ngs.name.table,
                              stringsAsFactors = F, header=T, sep='\t')
  ngs.name.table <- as.data.table(ngs.name.table)
  names <- strsplit(ngs.name.table$CompntName, ' ')
  ngs.name.table$CompntName <- unlist(lapply(names, function(x){x[[1]]}))
  setkey(ngs.name.table, 'CompntName')
  ngs.ref.align$QryContigID <- ngs.name.table[ngs.ref.align$QryContigID, CompntId, mult='first']
  ngs.ref.align$RefcontigID <- substr(ngs.ref.align$RefcontigID, 4, nchar(ngs.ref.align$RefcontigID))
  ngs.ref.align$QryContigID <- as.character(ngs.ref.align$QryContigID)
  ngs.ref.align$RefcontigID <- as.character(ngs.ref.align$RefcontigID)
  ngs.ref.align$RefLen <- ngs.ref.align$RefLen + runif(nrow(ngs.ref.align))
  ngs.ref.align <- remove.duplicated.align(ngs.ref.align, delta.score = 500, align.len = 0.00000001, bestRef = F)
  align.count <- table(ngs.ref.align$QryContigID)
  print(paste0("Removing alignment map to multiple locations: ", sum(align.count > 1)))
  ngs.ref.align <- ngs.ref.align[QryContigID %in% names(align.count[align.count ==1])]
  return(ngs.ref.align)
}


#this function compare the NGS anchoring ordering of two xmap
#the qry is usually anchoring from a scaffold and ref is usually to a reference seq and
#the accuracy of anchoring is calculated
evalScaffoldOrder <- function(qry.xmap, ref.xmap){
  #browser()
  #select alignments that are both in qry and reference alignment
  if(length(unique(qry.xmap$QryContigID)) == nrow(qry.xmap)/2){
    qry.xmap <- qry.xmap[LabelChannel==1]
  }
  if(length(unique(qry.xmap$QryContigID)) == nrow(qry.xmap)/2){
    ref.xmap <- ref.xmap[LabelChannel==1]
  }
  qry.xmap.intersect <- qry.xmap[(QryContigID %in% ref.xmap$QryContigID)==TRUE]
  ref.xmap.intersect <- ref.xmap[(QryContigID %in% qry.xmap$QryContigID)==TRUE]
  #using mid point of the alignment to order NGS contigs
  qry.xmap.intersect <- getContigCoordinate(qry.xmap.intersect)
  ref.xmap.intersect <- getContigCoordinate(ref.xmap.intersect)
  setkeyv(qry.xmap.intersect, c("RefcontigID", "Location"))
  qry.aligned.ids <- qry.xmap.intersect[, list(RefcontigID, QryContigID)]
  scaffolded.ngs.pairs <- qry.aligned.ids[,list(QryContigID[-length(QryContigID)],  QryContigID[-1]), by=RefcontigID]
  print(paste0("Total NGS pairs can be evaluated for order accuracy in scaffolds: ", nrow(scaffolded.ngs.pairs)))
  setkeyv(ref.xmap.intersect, c("RefcontigID", "Location"))
  ref.order <- ref.xmap.intersect[,list(QryContigID, RefcontigID)]
  ref.order <- ref.order[,Indice:=1:nrow(ref.order)]
  ref.order <- ref.order[,isBoundary := (Indice == Indice[[1]] | Indice == tail(Indice,1)), by=RefcontigID]
  setkey(ref.order, QryContigID)
  #browser()
  order.dist <- scaffolded.ngs.pairs[, dist:=ref.order[V2, Indice] - ref.order[V1, Indice]]
  order.dist <- scaffolded.ngs.pairs[, isBoundary:=ref.order[V2, isBoundary] & ref.order[V1, isBoundary]]
  order.dist <- scaffolded.ngs.pairs[, sameRef:=ref.order[V2, RefcontigID] == ref.order[V1, RefcontigID]]
  next.indices <- ref.order[scaffolded.ngs.pairs$V1, Indice+1]
  prev.indices <- ref.order[scaffolded.ngs.pairs$V1, Indice-1]
  next.indices <- ifelse(next.indices > max(ref.order$Indice), max(ref.order$Indice), next.indices)
  prev.indices <- ifelse(prev.indices > 0, prev.indices, 1)
  setkeyv(ref.order, 'Indice')
  order.dist$V3 <- ref.order[next.indices]$QryContigID
  order.dist$V4 <- ref.order[prev.indices]$QryContigID

  n <- nrow(order.dist) - order.dist[,sum(abs(dist)==1)]
  print(paste0("Total NGS aligned in scaffold ", nrow(qry.xmap), " Total NGS aligned in reference ", nrow(ref.xmap)))
  print(paste0("Number of contigs out of order: ", n,
               " Anchoring order accuracy: ", (1-n/nrow(order.dist)),
               " Total evaluated NGS pairs: ", nrow(order.dist)))

  print(paste0("Number of contigs out of order: ", nrow(order.dist[abs(dist) > 1 &  isBoundary ==F]),
              " Intra Contig Anchoring order accuracy: ", (1-nrow(order.dist[abs(dist) > 1 & isBoundary ==F])/nrow(order.dist[isBoundary==F])),
               " Total evaluated NGS pairs: ", nrow(order.dist[sameRef==T])))
  return(order.dist)
}

evalScaffoldDistance <- function(qry.xmap, ref.xmap){
  if(length(unique(qry.xmap$QryContigID)) == nrow(qry.xmap)/2){
    qry.xmap <- qry.xmap[LabelChannel==1]
  }
  if(length(unique(qry.xmap$QryContigID)) == nrow(qry.xmap)/2){
    ref.xmap <- ref.xmap[LabelChannel==1]
  }
  #browser()
  qry.xmap.intersect <- qry.xmap[(QryContigID %in% ref.xmap$QryContigID)==TRUE]
  ref.xmap.intersect <- ref.xmap[(QryContigID %in% qry.xmap$QryContigID)==TRUE]
  qry.xmap.intersect <- getContigCoordinate(qry.xmap.intersect)
  ref.xmap.intersect <- getContigCoordinate(ref.xmap.intersect)
  setkeyv(qry.xmap.intersect, c("RefcontigID", "Location"))
  setkeyv(ref.xmap.intersect, c("QryContigID"))
  qry.xmap.intersect <- getProjAlignPos(qry.xmap.intersect)

  consec.contig.pair <- qry.xmap.intersect[, .(but.last(QryContigID), but.first(QryContigID), but.last(Orientation), but.first(Orientation)), by=RefcontigID]
  setnames(consec.contig.pair, c("RefcontigID", "QryContigID1", "QryContigID2", "QryOrientation1", "QryOrientation2"))
  consec.contig.pair$scaffolddist <- qry.xmap.intersect[, but.first(ProjRefStart) - but.last(ProjRefEnd), by=RefcontigID]$V1
  consec.contig.pair$aligndist <- qry.xmap.intersect[, but.first(RefStartPos) - but.last(RefEndPos), by=RefcontigID]$V1

  ref.dist.forward <- ref.xmap.intersect[consec.contig.pair$QryContigID2, RefStartPos] - ref.xmap.intersect[consec.contig.pair$QryContigID1, RefEndPos]
  ref.dist.rev <- ref.xmap.intersect[consec.contig.pair$QryContigID1, RefStartPos] - ref.xmap.intersect[consec.contig.pair$QryContigID2, RefEndPos]
  consec.contig.pair$sameRef <- ref.xmap.intersect[consec.contig.pair$QryContigID1, RefcontigID] ==  ref.xmap.intersect[consec.contig.pair$QryContigID2, RefcontigID]
  consec.contig.pair$RefOrientation1=ref.xmap.intersect[consec.contig.pair$QryContigID1, Orientation]
  consec.contig.pair$RefOrientation2=ref.xmap.intersect[consec.contig.pair$QryContigID2, Orientation]

  inds.forward <- consec.contig.pair[, QryOrientation1 == RefOrientation1 & QryOrientation2 == RefOrientation2]
  inds.rev <- consec.contig.pair[, QryOrientation1 != RefOrientation1 & QryOrientation2 != RefOrientation2]

  consec.contig.pair$ref.dist <- ifelse(inds.forward, ref.dist.forward, ref.dist.rev)
  consec.contig.pair <- consec.contig.pair[sameRef==T]

  #browser()
  print(paste0("Number of contigs with incorrect orientation: ", sum(!(inds.forward | inds.rev))))
  print(paste0("Number of overlapping contigs: "))

  consec.contig.pair[,difference:=scaffolddist - ref.dist]
  hist(log10(consec.contig.pair$ref.dist))
  h <- consec.contig.pair[abs(difference) < 2e3, hist(difference, 1000)]
  h2 <- consec.contig.pair[, hist(difference/(aligndist), c(-10e6, seq(-1, 1, 0.01), 10e6), xlim=c(-0.20, 0.2))]
  print(paste0("median of difference ", consec.contig.pair[, median((difference))]))
  print(paste0("mad of difference ", consec.contig.pair[, mad((difference))]))
  print("Quantiles of absolute difference: ")
  print(quantile(abs(consec.contig.pair$difference), seq(0,1,0.1)))
  print(paste0("median of percent difference ", consec.contig.pair[, median(difference/abs(aligndist+0.1)*100)]))
  print(paste0("mad of percent difference ", consec.contig.pair[, mad(difference/abs(aligndist+0.1)*100)]))
  print("Quantiles of percent absolute difference: ")
  print(quantile(consec.contig.pair[, abs(difference/(aligndist+0.1))*100], seq(0,1,0.1)))
  return(consec.contig.pair)
}

getProjAlignPos <- function(xmap){
  #browser()
  xmap[, ProjRefStart := ifelse(Orientation=='+', RefStartPos - QryStartPos, RefStartPos - (QryLen - QryStartPos))]
  xmap[, ProjRefEnd := ifelse(Orientation == '+', RefEndPos + (QryLen - QryEndPos), RefEndPos + QryEndPos)]
  return(xmap)
}


eval.block.accuracy <- function(align.file){
  ngs.read.aligns <- read.table(align.file, stringsAsFactors = F)
  setnames(ngs.read.aligns, c('QryContigID', 'RefcontigID', 'Refstrand', 'Qrystrand', 'Confidence',
                              'PercentSimilarity', 'RefStartPos', 'RefEndPos', 'RefLen', 'QryStartPos', 'QryEndPos', 'QryLen', 'NCells'))
  ngs.read.aligns <- as.data.table(ngs.read.aligns)
  ngs.ids <- unique(ngs.read.aligns$QryContigID)
  getIdInd <- function(ids){
    ids.split <- strsplit(ids, '/')
    ids.split <- lapply(ids.split, function(x){strsplit(x[[1]], '_')})
    id.ind <- unlist(lapply(ids.split, function(x){tail(x[[1]],1)}))
    return(id.ind)
  }
  #browser()
  ngs.read.aligns$Confidence <- -1*ngs.read.aligns$Confidence
  read.aligns.filtered <- remove.duplicated.align(ngs.read.aligns, delta.score = 0.01, align.len = 0.0, bestRef = F)
  read.aligns.filtered[,NumAligns:=length(Confidence), by=QryContigID]
  setkey(read.aligns.filtered, 'QryContigID')
  read.aligns.final <- read.aligns.filtered[ngs.ids, mult='first']
  read.aligns.final[,dist:=c(0,diff(RefStartPos))]
  read.aligns.final[, sameChr:=c(0, generalized.diff(RefcontigID, `==`))]
  read.aligns.final[, prevNumAligns:=c(0, but.last(NumAligns))]
  read.aligns.final$IdInd <- as.integer(getIdInd(read.aligns.final$QryContigID))
  read.aligns.final[, IdDist:=c(0,diff(IdInd))==1]
  read.aligns.final[, usePair:=sameChr & prevNumAligns== 1 & NumAligns== 1 & IdDist==1]
  block.dist <- read.aligns.final[sameChr & prevNumAligns== 1 & NumAligns== 1 & IdDist==1, abs(dist)]
  print(paste0("Scaffold block-based accuracy ", sum(abs(block.dist) > 90e3 & abs(block.dist) < 110e3)/length(block.dist)))
  return(block.dist)
}

#this function contain scripts to diagnose the block-accuracy of scaffolds
testBloackAccuray <- function(){
  align1 <- readxmap('~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSPQI_0kb_0labels_1_errEst.xmap')
  align2 <- readxmap('~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSSSI_0kb_0labels_1_errEst.xmap')
  align1 <- as.data.table(align1)
  align2 <- as.data.table(align2)
  align1$Confidence <- align1$Confidence + runif(nrow(align1))*0.1
  align.merge <- rbind(align1, align2)
  align.merge <- remove.duplicated.align(align.merge, 0.000005, 0.0, T)
  exportXmap(align.merge, "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_singleEnzymeAlign_Merged.xmap")

  two.color.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash1_errEst.xmap"
  two.color.file2 <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash2_errEst.xmap"


  align.two.color <- as.data.table(get.merged.twocolor.align(two.color.file, two.color.file2))

  align.two.color1 <- align.two.color[Confidence > 10 & QryStartPos != QryEndPos & LabelChannel==1]
  align.two.color2 <- align.two.color[Confidence > 10 & QryStartPos != QryEndPos & LabelChannel==2 & !QryContigID %in% align.two.color1$QryContigID]
  align.two.color <- rbind(align.two.color1, align.two.color2)
  exportXmap(align.two.color,
             "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_twocolor_reexported.xmap")


  align.two.color.select <- align.two.color[QryContigID %in% align.merge$QryContigID]
  align.two.color.select <- align.two.color.select[LabelChannel == 1]
  exportXmap(align.two.color.select, "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_twocolor_intersect_singleEnzyAln.xmap")

  align.two.color.only <- align.two.color[!(QryContigID %in% align.merge$QryContigID)]

  #running block-accuracy
  align.file <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/NGSAlignValidate/Block_Accuracy/hybrid_Sandwich_with_pacbio_rerun10032016_block_reads_block_reads.fa.aln"
  reads.align.final <- eval.block.accuracy(align.file)

}

