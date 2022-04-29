#20150710 EL - making summary statistics in a nice format - does path-level summmary and scaffold-level summary
#20151027 EL - modifying from run_2-color_scaffold_info.R

# rm(list=ls())
#
# #general scripts
# scripts_dir <- "/home/users/elam/rscripts/"
# source(paste(scripts_dir,"paste0.R",sep=""))
# source(paste0(scripts_dir,"readmaps.R")) #util functions to read _map formats
# source(paste0(scripts_dir,"calc_n50.R"))
#
# #2-color specific scripts
# scripts_dir <- "/home/users/elam/rscripts/2col_dev/current/"
# source(paste0(scripts_dir,"utils.R"))
# source(paste0(scripts_dir,"do_validation.R"))

merge_output <- function(scaffolds) {

  summaries <- lapply(scaffolds,function(scaffold) {
    return(scaffold[["summary"]])
  })
  summary <- do.call(rbind,summaries)

  paths <- lapply(1:length(scaffolds),function(i) {
    scaffold <- scaffolds[[i]]
    path <- scaffold[["path_coordinates"]]
    path <- cbind("path_idx"=i,"node_idx"=1:nrow(path),path)

    return(path)
  })
  paths <- do.call(rbind,paths)

  outlist <- list("summary"=summary,
                  "paths"=paths)

  return(outlist)
}

check_against_reference <- function(path_df,alignment_to_ref) {
  path_split <- split(path_df,path_df[,"node_type"])

  path_annotated <- lapply(path_split,function(nodes) {
    node_type <- as.character(nodes[1,"node_type"])
#     print(node_type)

    if (node_type=="asm1") {
      xmap <- alignment_to_ref[[1]]
    } else {
      xmap <- alignment_to_ref[[2]]
    }

#     print(xmap)
#     stop()
    if(length(xmap) > 0){
        xmap <- readxmap(xmap)[,c(1:9,12)]
     }else{
        #when no align_to_ref supplied we create a dummy alignment
	#so it does not crash
	xmap <- data.frame(XmapEntryID=1, QryContigID=-100, RefcontigID=-100, QryStartPos=0, QryEndPos=0, RefStartPos=0, RefEndPos=0, Orientation='+', Confidence=0, HitEnum='', QryLen=0, RefLen=0, LabelChannel=1, Alignment='')
     }

    class(nodes[,"node"]) <- "numeric"
    class(xmap[,"QryContigID"]) <- "numeric"

    merged <- merge(nodes,xmap,by.x="node",by.y="QryContigID",sort=F,all.x=T)

#     refid <- merged[,"RefcontigID"]
#     refid[is.na(refid)] <- -1 #if a map is not aligned
#
#     stopifnot(nrow(nodes)==length(refid))
#
#     out_df <- cbind(nodes,"chr"=refid)

    return(merged)
    return(out_df)
  })
  path_annotated <- do.call(rbind,path_annotated)

  path_annotated <- path_annotated[order(path_annotated[,"path_idx"],path_annotated[,"node_idx"]),]

  return(path_annotated)
}

#08/29/2016
#change to return whole dataframe of non-scaffolded contigs instead of just the length
get_non_scaffolded_lengths <- function(path_df,asm1,asm2) {
  ncol1 <- ncol(asm1)
  ncol2 <- ncol(asm2)

  if (ncol2>ncol1) {
    asm1 <- asm1[!duplicated(asm1[,"CMapId"]),]
    asm2 <- asm2[!duplicated(asm2[,"CMapId"]),-((ncol(asm2)-2):ncol(asm2))]
  } else if (ncol2<ncol1) {
    asm1 <- asm1[!duplicated(asm1[,"CMapId"]),-((ncol(asm1)-2):ncol(asm1))]
    asm2 <- asm2[!duplicated(asm2[,"CMapId"]),]
  } else {
    asm1 <- asm1[!duplicated(asm1[,"CMapId"]),]
    asm2 <- asm2[!duplicated(asm2[,"CMapId"]),]
  }
  #browser()
  asm1_nodes <- path_df[path_df[,"node_type"]=="asm1",]
  asm2_nodes <- path_df[path_df[,"node_type"]=="asm2",]
  #   print(head(asm1_nodes))

  asm1_maps <- unique(asm1[,"CMapId"])
  asm2_maps <- unique(asm2[,"CMapId"])
  #   print(head(asm1_maps))

  merged1 <- merge(asm1,asm1_nodes,by.x="CMapId",by.y="node",all.x=T,sort=F)
  merged1$node_type[is.na(merged1$node_type)] <- "asm1"
  merged2 <- merge(asm2,asm2_nodes,by.x="CMapId",by.y="node",all.x=T,sort=F)
  merged2$node_type[is.na(merged2$node_type)] <- "asm2"
  merged <- rbind(merged1,merged2)

  non_scaffolded <- merged[is.na(merged[,"path_idx"]),]

  #merged <- non_scaffolded[,"ContigLength"] #lengths of non-scaffolded maps
  #
  #return(merged)
  #browser()
  return(non_scaffolded)
}

get_path_length <- function(path) { #util function to calculate path length
  min_pos <- min(path[,"rel_start"],path[,"rel_end"])
  max_pos <- max(path[,"rel_start"],path[,"rel_end"])
  total_length <- max_pos - min_pos
  # print(total_length)

  out_list <- list("min"=min_pos,
                   "max"=max_pos,
                   "total"=total_length)

  return(total_length)
  return(out_list)
}

#util function to calculate how much of the path is covered by a particular node type (asm1 or asm2)
#uses the library "intervals"
get_path_coverage <- function(path,node_type) {
  path_length <- get_path_length(path)

  path_selected <- path[path[,"node_type"]==node_type,]

  for (i in 1:nrow(path_selected)) {
    if (path_selected[i,"rel_start"]>path_selected[i,"rel_end"]) {
      tmp <- path_selected[i,"rel_start"]
      path_selected[i,"rel_start"] <- path_selected[i,"rel_end"]
      path_selected[i,"rel_end"] <- tmp
    }
  }

  path_intervals <- Intervals(as.matrix(path_selected[,c("rel_start","rel_end")]))
  path_union <- interval_union(path_intervals)

  union_size <- size(path_union)
  fraction_covered <- sum(union_size)/path_length

  # print(round(fraction_covered,2))
  return(fraction_covered)
  return(path_union)
}

find_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x,ux)))]
}

summarize_path1 <- function(path) {

  n_nodes <- nrow(path)
  path_length <- get_path_length(path)
  cov_asm1 <- get_path_coverage(path,"asm1")
  cov_asm2 <- get_path_coverage(path,"asm2")
#   n_chr <- length(unique(path[,"chr"]))
#   major_chr <- find_mode(path[,"chr"])
#   fraction_chimeric <- 1-sort(table(path[,"chr"]),decreasing=T)[1]/nrow(path)

  n_aligned <- sum(!is.na(path[,"RefcontigID"])) #number of unaligned maps
  percent_aligned <- n_aligned/n_nodes

  uniq_chr <- unique(path[,"RefcontigID"])
  uniq_chr <- uniq_chr[!is.na(uniq_chr)]
  n_chr <- length(uniq_chr)

  chrs <- path[,"RefcontigID"]
  chrs <- chrs[!is.na(chrs)]

  major_chr <- find_mode(chrs)
  n_chim_nodes <- length(chrs)-sort(table(chrs),decreasing=T)[1]
  fraction_chimeric <- n_chim_nodes/length(chrs)

  out_df <- data.frame("n_nodes"=n_nodes,
                       "n_aligned"=n_aligned,
                       "percent_aligned"=percent_aligned,
                       "path_length"=path_length,
                       "cov_asm1"=cov_asm1,
                       "cov_asm2"=cov_asm2,
                       "n_chr"=n_chr,
                       "major_chr"=major_chr,
                       "n_chim_nodes"=n_chim_nodes,
                       "fraction_chimeric"=fraction_chimeric)

  return(out_df)
}

make_summary1 <- function(paths,asm1,asm2,alignment_to_ref) { #optimized for speed
  merged <- merge_output(paths)
  merged_annotated1 <- check_against_reference(merged[[2]],alignment_to_ref) #check where the maps align to on reference
  save(list="merged_annotated1", file=paste0("merged_annotated.RData"))
  paths_split <- split(merged_annotated1,merged_annotated1[,"path_idx"])
  #path-based summary
  paths_annotated <- mclapply(paths_split,function(scaffold) {
    path_annotated <- summarize_path1(scaffold)

    return(path_annotated)
  },mc.cores=8)
  paths_annotated <- do.call(rbind,paths_annotated)
  paths_annotated <- cbind("indx"=(1:nrow(paths_annotated)),paths_annotated)

  summary <- merged[[1]]
  paths_annotated <- cbind(paths_annotated,summary)

  #scaffold-based summary
  asm1 <- readcmap(asm1)
  asm2 <- readcmap(asm2)

  trimmed <- trim_cmap(asm1, asm2)
  asm1 <- trimmed[[1]]
  asm2 <- trimmed[[2]]

  asm1_lens <- get_lengths(asm1)
  asm2_lens <- get_lengths(asm2)

  asm1_n50 <- calc_n50(asm1_lens)
  asm2_n50 <- calc_n50(asm2_lens)
  combined_n50 <- calc_n50(c(asm1_lens,asm2_lens))

  non_scaffolded <- get_non_scaffolded_lengths(merged[[2]],asm1,asm2)
  non_scaffolded_lengths <- non_scaffolded[, "ContigLength"]
  longest_map_len_i <- max(c(asm1_lens, asm2_lens)) #longest map

  path_lengths <- paths_annotated[,"path_length"] #get scaffold lengths
  scaffolded_n50 <- calc_n50(c(path_lengths)) #n50 of scaffolded maps
  assembly_n50 <- calc_n50(c(path_lengths,non_scaffolded_lengths)) #n50 of scaffolded and unscaffold maps
  assembly_mean_len <- mean(c(path_lengths,non_scaffolded_lengths)) #mean of lengths of scaffolded and unscaffold maps

  max_map_len <- max(c(path_lengths,non_scaffolded_lengths)) #longest map after scaffolding

  scaffolded_size <- sum(path_lengths) #total length of scaffolded maps
  non_scaffolded_size <- sum(non_scaffolded_lengths) #total length of scaffolded maps

  n_scaffolded <- length(path_lengths) #number of scaffolds
  n_non_scaffolded <- length(non_scaffolded_lengths) #number of unscaffolded maps

  n_nodes <- nrow(merged[[2]]) #number of nodes in scaffolds

  fold_improvement <- assembly_n50/combined_n50 #fold improvement: before vs after n50

  n_all_maps_i <- length(c(asm1_lens,asm2_lens)) #sanity check to make sure we are not missing any maps
  n_all_maps_f <- n_nodes + n_non_scaffolded

  stopifnot(n_all_maps_i==n_all_maps_f)

  out_df <- data.frame("longest_map_len_i"=longest_map_len_i,
                       "longest_map_len_f"=max_map_len,
                       "n_nodes"=n_nodes,
                       "n_scaffolded"=n_scaffolded,
                       "scaffolded_size"=scaffolded_size,
                       "n_non_scaffolded"=n_non_scaffolded,
                       "non_scafflded_size"=non_scaffolded_size,
                       "asm1_n50"=asm1_n50,
                       "asm2_n50"=asm2_n50,
                       "combined_n50"=combined_n50,
                       "scaffolded_n50"=scaffolded_n50,
                       "overall_n50"=assembly_n50,
                       "overall_mean_len"=assembly_mean_len,
                       "fold_improvement"=fold_improvement)

  outlist <- list("paths_annotated"=paths_annotated,
                  "scaffold_summary"=out_df,
                  "non_scaffolds"=non_scaffolded)

  return(outlist)
  return(paths_annotated)
}

#parse final output
#chimeric rates (fraction of scaffolds that were "chimeric" [mapping to more than one chromosome])
get_chimeric_rates <- function(final_output) {
  chimeric_rates <- lapply(final_output,function(run) {
    path_info <- run[[2]][[1]]

    fraction_chimeric <- sum(path_info[,"fraction_chimeric"]>0)/nrow(path_info)

    return(fraction_chimeric)
  })
  chimeric_rates <- unlist(chimeric_rates)

  return(chimeric_rates)
}

#extract summary information and merged with mapping rates
process_output <- function(final_output,rates) {
  summaries <- lapply(final_output,function(run) {
    outsummary <- run[[2]][[2]]

    return(outsummary)
  })
  summaries <- do.call(rbind,summaries)

  outdf <- cbind(rates,summaries)

  return(outdf)
}

#plot eff cov vs N50-fold improvement
plot_summary <- function(summary,outfile=NULL,plot_flag=F) {

  columns_to_keep <- c("prefixes","chr1_maprate","chr2_maprate","eff_cov","longest_map_len_i","longest_map_len_f","asm1_n50","asm2_n50","combined_n50","overall_n50","fold_improvement")
  summary_selected <- summary[,colnames(summary)%in%columns_to_keep]

  if (plot_flag) {
    if (!is.null(outfile)) {
      pdf(outfile,width=12,height=8)
    }

    xrange <- range(summary_selected["eff_cov"])
    yrange <- range(summary_selected["fold_improvement"])

    try({
      plot(0,type="n",xlim=xrange,ylim=yrange,xlab="Effective coverage",ylab="N50 fold improvement",main="2-color sandwich of single-color assemblies") #make empty plot
      points(summary_selected[,"eff_cov"],summary_selected[,"fold_improvement"],pch=16,col=rgb(0,0,1,0.5))
      text(summary_selected[,"eff_cov"]+1,summary_selected[,"fold_improvement"]+0.01,labels=paste0(round(summary_selected[,"eff_cov"],1),"X"),cex=0.5,col=rgb(0,0,0,0.8))
    })

    if (!is.null(outfile)) {
      dev.off()
    }
  }

  return(summary_selected)
}

# #20151031 EL - debugging
# load("/home/users/elam/20151030_2-col_V3_sandwich/20151030_v3_sub110X_outlist.RData")
#
# #check the path first
# asm1 <- "/home/users/elam/data/20151027_NA12878_BspQI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
# asm2 <- "/home/users/elam/data/20151027_NA12878_BssSI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
#
# output <- check_assembly_source1(outlist[[1]],asm1,asm2)
#
# #check against reference
# paths <- outlist[[1]]
#
# bspqi_xmap <- "/home/users/elam/data/20151027_NA12878_BspQI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap"
# bsssi_xmap <- "/home/users/elam/data/20151027_NA12878_BssSI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap"
#
# alignment_to_ref <- list("asm1"=bspqi_xmap,
#                          "asm2"=bsssi_xmap)
#
# # alignment_to_ref <- list("asm1"=bsssi_xmap,
# #                          "asm2"=bspqi_xmap)
#
# merged <- merge_output(paths)
# merged_annotated1 <- check_against_reference(merged[[2]],alignment_to_ref) #check where the maps align to on reference
#
# paths_split <- split(merged_annotated1,merged_annotated1[,"path_idx"])
#
# #path-based summary
# paths_annotated <- mclapply(paths_split,function(scaffold) {
#   path_annotated <- summarize_path1(scaffold)
#
#   return(path_annotated)
# },mc.cores=8)
# paths_annotated <- do.call(rbind,paths_annotated)
# paths_annotated <- cbind("indx"=(1:nrow(paths_annotated)),paths_annotated)
#
# summary <- merged[[1]]
# paths_annotated <- cbind(paths_annotated,summary)

# #main
# load("/home/users/elam/20150911_2-col_actual_data/sop_and_ml_test/Human_willsample_twocolor_combined_paths.RData")
#
# #where the RData files will be output
# setwd("/home/users3/elam/20150911_2-col_actual_data/sop_and_ml_test") #20151011 - TESTING
#
# align_dir <- paste0(getwd(),"/")
#
# #input parameters
# #prefix input (based on the bnx prefixes)
# prefix <- "Human_willsample_twocolor_combined" #20151011 - TESTING
#
# bspqi_xmap <- "/home/users3/elam/data/NA12891_bspqi_ml_50X_4/output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap" #new 50X assembly (ML)
# bsssi_xmap <- "/home/users3/elam/data/NA12891_bsssi_sop_50X_2/output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap" #new 50X assembly
#
# alignment_to_ref <- list("asm1"=bsssi_xmap,
#                          "asm2"=bspqi_xmap)
#
# asm1 <- paste0(align_dir,prefix,"_1_bppAdj.cmap")
# asm2 <- paste0(align_dir,prefix,"_2_bppAdj.cmap")
#
# path_summary <- make_summary1(paths,asm1,asm2,alignment_to_ref)
