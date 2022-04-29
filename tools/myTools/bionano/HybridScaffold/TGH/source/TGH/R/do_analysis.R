#20151026 EL - main function to drive analysis

run_analysis <- function(input,input_type="rdata",
                         output_dir,
                         alignment_to_ref,
                         selected_chr=NULL,
                         check_ans_flag=F, twocolor.as.ref=F, min.align=1, min.score.T=13) {

  #get alignment data
  save_flag <- T
  indata <- get_input_alignment(input,input_type,save_flag, twocolor.as.ref=twocolor.as.ref, score.T = min.score.T)

  #checking if we have enough "glues" aligned to assembly
  if(is.null(indata) | nrow(indata) < min.align){
     stop("Not sufficient numbers of NGS contigs aligned to both assemblies, stoping sandwich assembly")
  }

  indata[,"id_str"] <- paste0(indata[,"RefcontigID.x"],"_",indata[,"RefcontigID.y"]) #get str representing pair #####SHOULD NOT BE NEEDED ANYMORE

  if(!is.null(alignment_to_ref) && length(alignment_to_ref) > 1){
    xmap1 <- alignment_to_ref[[1]] #alignment of assembly 1 maps to reference
    xmap2 <- alignment_to_ref[[2]] #alignment of assembly 2 maps to reference
  }

  #get answers
  if (check_ans_flag) {
    answers <- check_answers(xmap1,xmap2,parallel_flag=F) #get known pairs based on a certain threshold
  }
  #browser()

  #select alignments from one chromosome (based on map-to-reference alignments)
  #for testing
  if (!is.null(selected_chr)) {
    stopifnot(length(selected_chr)==1)

    print(paste0("selecting alignments from chr",selected_chr))

    xmap1_selected <- select_alignement_chr(xmap1,selected_chr)
    xmap2_selected <- select_alignement_chr(xmap2,selected_chr)

    nrows_initial <- nrow(indata)
    indata <- indata[(indata[,"RefcontigID.x"]%in%xmap1_selected[,"QryContigID"] & indata[,"RefcontigID.y"]%in%xmap2_selected[,"QryContigID"]),]
    nrows_final <- nrow(indata)

    nrows_diff <- nrows_initial-nrows_final

    print(paste0("number of rows reduced by ",nrows_diff))
    print(paste0("number of rows after filtering: ",nrows_final))
  }
  #apply filtering
  indata <- add_offsets(indata)

  indata_filtered <- filter_pairs(indata,1) #filter based on the number of supporting molecules
  indata_filtered <- filter_by_correlation(indata_filtered) #filter based on correlation

  indata_collapsed <- collapse_duplicates(indata_filtered) #collapse pairs

  if (check_ans_flag) {
    indata_annotated <- annotate_pairs(indata_collapsed,answers)
    fraction_true_edges <- count_true_edges(indata_annotated)
    indata_renumbered <- fix_numbering(indata_annotated) #get non-redundant ids
  } else {
    indata_renumbered <- fix_numbering(indata_collapsed) #get non-redundant ids
  }

  #make clusters
  clusters <- get_clusters(indata_renumbered)
  clusters_g <- mclapply(clusters,get_bipartite,mc.cores=64)


  save_flag <- F
  if (save_flag) {
    if (input_type=="prefix") {
      infile <- paste0(input[["prefix"]],"_clusters_g.RData")
      if (!file.exists(infile)) { #save if file does not exist
        print(paste0("saving cluster_g data to ",infile))
        save(list="clusters_g",file=paste0(input[["prefix"]],"_clusters_g.RData"))
      } else { #stop if file already exists
        stop(paste0(infile," already exists."))
      }
    }
  }

  #fit and get offsets
  plot_flag <- F
  parallel_flag <- F
  cores <- 12
  if (plot_flag) {
    if (input_type=="prefix") {
      pdf_file <- paste0(input[["prefix"]],"_indata_analyze_alignments.pdf")
    } else {
      pdf_file <- paste0(output_dir,"indata_analyze_alignments.pdf")
    }
    pdf(pdf_file,width=12,height=8)
  }
  try({
    indata_split <- split(indata_filtered,indata_filtered[,"id_str"])
    pos_offset_list <- plot_alignments(indata_split,plot_flag=plot_flag,parallel_flag=parallel_flag, cores=cores)
  })


  if (plot_flag) {
    dev.off()
  }

  pos_offsets <- do.call(rbind,pos_offset_list)
  pos_offsets <- as.data.frame(pos_offsets,stringsAsFactors=F)
  colnames(pos_offsets) <- c("id_str","intercept","slope")
  #browser()
  save(list="pos_offsets", file=paste0(input[["prefix"]], "_pos_offsets.RData"))

  #analyze each path
  plot_flag <- F
  cores <- 12
  if (plot_flag) {
    if (input_type=="prefix") {
      pdf_file <- paste0(input[["prefix"]],"_indata_paths.pdf")
    } else {
      pdf_file <- paste0(output_dir,"indata_paths.pdf")
    }
    pdf(pdf_file,width=12,height=8)
  }

  try({
     #browser()
     paths <- lapply(clusters_g,get_path1,pos_offsets,indata_filtered,alignment_to_ref,plot_flag=plot_flag)
     #paths <- mclapply(clusters_g,get_path1,pos_offsets,indata_filtered,alignment_to_ref,plot_flag=plot_flag,mc.cores=cores)

#     #testing to make sure the nodes are of the right types
#     asm1 <- "/home/users/elam/data/20151027_NA12878_BspQI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
#     asm2 <- "/home/users/elam/data/20151027_NA12878_BssSI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
#     mode <- 2
#
#     check_source <- check_assembly_source1(paths,asm1,asm2,mode)
    if (input_type=="prefix") {
      save(list="paths",file=paste0(input[["prefix"]],"_paths.RData"))
    }
  })

  if (plot_flag) {
    dev.off()
  }

  return(paths)
}



#@DEPRECATED note this function is replace by run.pre.sandwich() in runHybridScaffoldV3.R
#this function only run the sandwich pipeline up to the point where the pos_offset data.frame is generated
#it it used to generate neccesary information for the bng vs bng contigs conflict resolution
run_pre_analysis <- function(input,input_type="rdata",
                         output_dir,
                         alignment_to_ref,
                         selected_chr=NULL,
                         check_ans_flag=F, twocolor.as.ref=F, min.align=1, min.score.T=13) {

  #get alignment data
  save_flag <- T
  indata <- get_input_alignment(input,input_type,save_flag, twocolor.as.ref=twocolor.as.ref, score.T = min.score.T)

  print(paste0("Number of two-color molecules aligned to both assemblies ", nrow(indata)))
  #checking if we have enough "glues" aligned to assembly
  if(is.null(indata) | nrow(indata) < min.align){
    stop("Not sufficient numbers of two-color molecules aligned to both assemblies, skipping sandwich scaffolding")
  }

  indata[,"id_str"] <- paste0(indata[,"RefcontigID.x"],"_",indata[,"RefcontigID.y"]) #get str representing pair #####SHOULD NOT BE NEEDED ANYMORE

  if(!is.null(alignment_to_ref) && length(alignment_to_ref) > 1){
    xmap1 <- alignment_to_ref[[1]] #alignment of assembly 1 maps to reference
    xmap2 <- alignment_to_ref[[2]] #alignment of assembly 2 maps to reference
  }

  #get answers
  if (check_ans_flag) {
    answers <- check_answers(xmap1,xmap2,parallel_flag=F) #get known pairs based on a certain threshold
  }
  #browser()

  #select alignments from one chromosome (based on map-to-reference alignments)
  #for testing
  if (!is.null(selected_chr)) {
    stopifnot(length(selected_chr)==1)

    print(paste0("selecting alignments from chr",selected_chr))

    xmap1_selected <- select_alignement_chr(xmap1,selected_chr)
    xmap2_selected <- select_alignement_chr(xmap2,selected_chr)

    nrows_initial <- nrow(indata)
    indata <- indata[(indata[,"RefcontigID.x"]%in%xmap1_selected[,"QryContigID"] & indata[,"RefcontigID.y"]%in%xmap2_selected[,"QryContigID"]),]
    nrows_final <- nrow(indata)

    nrows_diff <- nrows_initial-nrows_final

    print(paste0("number of rows reduced by ",nrows_diff))
    print(paste0("number of rows after filtering: ",nrows_final))
  }

  #apply filtering
  indata <- add_offsets(indata)

  indata_filtered <- filter_pairs(indata,1) #filter based on the number of supporting molecules
  indata_filtered <- filter_by_correlation(indata_filtered) #filter based on correlation

  indata_collapsed <- collapse_duplicates(indata_filtered) #collapse pairs

  if (check_ans_flag) {
    indata_annotated <- annotate_pairs(indata_collapsed,answers)
    fraction_true_edges <- count_true_edges(indata_annotated)
    indata_renumbered <- fix_numbering(indata_annotated) #get non-redundant ids
  } else {
    indata_renumbered <- fix_numbering(indata_collapsed) #get non-redundant ids
  }

  #make clusters
  clusters <- get_clusters(indata_renumbered)
  clusters_g <- mclapply(clusters,get_bipartite,mc.cores=64)


  save_flag <- F
  if (save_flag) {
    if (input_type=="prefix") {
      infile <- paste0(input[["prefix"]],"_clusters_g.RData")
      if (!file.exists(infile)) { #save if file does not exist
        print(paste0("saving cluster_g data to ",infile))
        save(list="clusters_g",file=paste0(input[["prefix"]],"_clusters_g.RData"))
      } else { #stop if file already exists
        stop(paste0(infile," already exists."))
      }
    }
  }

  #fit and get offsets
  plot_flag <- F
  parallel_flag <- F
  cores <- 12
  if (plot_flag) {
    if (input_type=="prefix") {
      pdf_file <- paste0(input[["prefix"]],"_indata_analyze_alignments.pdf")
    } else {
      pdf_file <- paste0(output_dir,"indata_analyze_alignments.pdf")
    }
    pdf(pdf_file,width=12,height=8)
  }
  try({
    indata_split <- split(indata_filtered,indata_filtered[,"id_str"])
    pos_offset_list <- plot_alignments(indata_split,plot_flag=plot_flag,parallel_flag=parallel_flag, cores=cores)
  })

  if (plot_flag) {
    dev.off()
  }
  #browser()
  pos_offsets <- do.call(rbind,pos_offset_list)
  pos_offsets <- as.data.frame(pos_offsets,stringsAsFactors=F)
  colnames(pos_offsets) <- c("id_str","intercept","slope")
  save(list="pos_offsets", file=paste0(input[["prefix"]], "_pos_offsets.RData"))

  return()

####################################### NOT RUN ###################################################
  #analyze each path
  plot_flag <- T
  cores <- 12
  if (plot_flag) {
    if (input_type=="prefix") {
      pdf_file <- paste0(input[["prefix"]],"_indata_paths.pdf")
    } else {
      pdf_file <- paste0(output_dir,"indata_paths.pdf")
    }
    pdf(pdf_file,width=12,height=8)
  }

  try({
     #browser()
     paths <- lapply(clusters_g,get_path1,pos_offsets,indata_filtered,alignment_to_ref,plot_flag=plot_flag)
     #paths <- mclapply(clusters_g,get_path1,pos_offsets,indata_filtered,alignment_to_ref,plot_flag=plot_flag,mc.cores=cores)

#     #testing to make sure the nodes are of the right types
#     asm1 <- "/home/users/elam/data/20151027_NA12878_BspQI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
#     asm2 <- "/home/users/elam/data/20151027_NA12878_BssSI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
#     mode <- 2
#
#     check_source <- check_assembly_source1(paths,asm1,asm2,mode)

    if (input_type=="prefix") {
      save(list="paths",file=paste0(input[["prefix"]],"_paths.RData"))
    }
  })

  if (plot_flag) {
    dev.off()
  }

  return(paths)
}
