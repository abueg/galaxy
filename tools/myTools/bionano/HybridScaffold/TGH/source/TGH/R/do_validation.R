#20151026 EL - checking "answers" based on alignment of the assemblies to the reference

# rm(list=ls())
#
# #general scripts
# scripts_dir <- "/home/users/elam/rscripts/"
# source(paste(scripts_dir,"paste0.R",sep=""))
# source(paste0(scripts_dir,"readmaps.R")) #util functions to read _map formats

not_in_answer <- function(answer_pairs,sim_data){ #see which observed pairs do not exist in the answer
  pair1 <- paste0(answer_pairs$asm1,"_",answer_pairs$asm2)
  pair2 <- paste0(answer_pairs$asm2,"_",answer_pairs$asm1)
  real_pairs <- c(pair1,pair2)

  return(sim_data[!(sim_data$id_str%in%real_pairs),])
}

annotate_pairs <- function(indata,answers) {
  F <- indata[,"id_str"]%in%answers[,"id_str"]
  #browser()
  indata[,"true_edge"] <- F
  #saveRDS(indata, '/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies/TestConnectivity/true_Pairs.RData')
  return(indata)
}

count_true_edges <- function(indata) {
  fraction <- (sum(indata[,"true_edge"]))/nrow(indata)

  print(paste0("fraction of true edges: ",round(fraction,2)))

  return(fraction)
}

check_overlap <- function(pair) {

  region1 <- pair[1,]
  region2 <- pair[2,]

  q_id1 <- region1[1]
  r_id1 <- region1[2]
  start1 <- region1[3]
  end1 <- region1[4]
  asm_id1 <- region1[5]

  q_id2 <- region2[1]
  r_id2 <- region2[2]
  start2 <- region2[3]
  end2 <- region2[4]
  asm_id2 <- region2[5]

  if (asm_id1==asm_id2) { #if both maps are from the same assembly
    return(-1e5) #arbitrary number
  }
  if (r_id1!=r_id2) { #if the maps are on different chromosomes
    return(-1e5)
  }

  block1 <- c(start1,end1)
  block2 <- c(start2,end2)
  overlap <- calculate_overlap1(block1,block2)

  return(overlap)
}

get_pairs <- function(combo,xmap,thresh=0) {

  pair <- xmap[combo,]
  overlap <- check_overlap(pair)

  if (overlap>thresh) {
    pair <- pair[order(pair[,"asm"]),]
    edge <- c(as.numeric(pair[1,"QryContigID"]),as.numeric(pair[2,"QryContigID"]),overlap)

    return(edge)
  } else {
    return(NULL)
  }
}

do_pairwise <- function(xmap,thresh=0,parallel_flag=T,cores=8) {

  nrows <- nrow(xmap)
  if(nrow(xmap) == 1){
      return(list(NULL))
  }
  combos <- combn(1:nrows,2,simplify=F)
  print(paste0(length(combos)," pairs in chr",xmap[,"RefcontigID"][1]))

  if (parallel_flag) {
    overlaps <- mclapply(combos,get_pairs,xmap,thresh,mc.cores=cores)
  } else {
    overlaps <- lapply(combos,get_pairs,xmap,thresh)
  }

  return(overlaps)
}

do_pairwise_wrapper <- function(xmap) {

  pairwise_output <- do_pairwise(xmap,-100000)
  pairwise_output <- pairwise_output[!sapply(pairwise_output,is.null)]
  pairwise_output <- do.call(rbind,pairwise_output)
  pairwise_output <- as.data.frame(pairwise_output)

  chr <- xmap[1,"RefcontigID"]
  pairwise_output <- as.data.frame(pairwise_output,stringsAsFactors=F)
  if(ncol(pairwise_output) != 3){
      return(data.frame())
  }
  colnames(pairwise_output) <- c("asm1","asm2","overlap")

  id_str <- paste(pairwise_output[,"asm1"],pairwise_output[,"asm2"],sep="_")
  pairwise_output <- cbind("chr"=chr,pairwise_output,"id_str"=id_str)

  return(pairwise_output)
}

check_answers <- function(xmap1,xmap2,parallel_flag=T,cores=8) {

  xmap1 <- readxmap(xmap1) #load data
  xmap2 <- readxmap(xmap2)

  #xmap1 <- xmap1[xmap1$Confidence > 12,]
  #xmap2 <- xmap2[xmap2$Confidence > 12,]
    
  xmap1 <- xmap1[,c(2,3,6,7)] #only keep certain columns
  xmap1[,"asm"] <- 1
  xmap2 <- xmap2[,c(2,3,6,7)]
  xmap2[,"asm"] <- 2

  xmap <- rbind(xmap1,xmap2)
  xmap_by_chr <- split(xmap,xmap[,"RefcontigID"]) #split by aligned chromosome

  #browser()  
  if (parallel_flag) {
    pairwise_output_all <- mclapply(xmap_by_chr,do_pairwise_wrapper,mc.cores=cores)
  } else {
    pairwise_output_all <- lapply(xmap_by_chr,do_pairwise_wrapper)
  }

  pairwise_output_all <- do.call(rbind,pairwise_output_all)

  return(pairwise_output_all)
}

analyze_residuals <- function(residuals) {
  par(mfrow=c(2,1))
  residuals <- residuals[!is.null(residuals)]

  hist(residuals,breaks=500)
  hist(residuals,breaks=500000,xlim=c(-100000,100000))

  quantiles <- quantile(abs(residuals),probs=seq(0.1,1,0.1))

  return(quantiles)
}

#for indata
#quickly check the assembly source of the indata ids
#expectation is that x_ids should correspond to asm1 ids
check_assembly_source <- function(indata,asm1_path,asm2_path) {
  asm1_ids <- unique(readcmap(asm1_path)[,"CMapId"])
  asm2_ids <- unique(readcmap(asm2_path)[,"CMapId"])

  x_ids <- unique(indata[,"RefcontigID.x"])
  y_ids <- unique(indata[,"RefcontigID.y"])

  #fraction of ids are not found
  x_to_asm1 <- sum(!x_ids%in%asm1_ids)/length(x_ids)
  x_to_asm2 <- sum(!x_ids%in%asm2_ids)/length(x_ids)

  y_to_asm1 <- sum(!y_ids%in%asm1_ids)/length(y_ids)
  y_to_asm2 <- sum(!y_ids%in%asm2_ids)/length(y_ids)

  outlist <- list("x_to_asm1"=x_to_asm1,
                  "x_to_asm2"=x_to_asm2,
                  "y_to_asm1"=y_to_asm1,
                  "y_to_asm2"=y_to_asm2)

  return(outlist)
}

# load("/home/users/elam/20151030_2-col_V3_sandwich/20151030_v2_sub30X.RData")
# load("/home/users/elam/20151030_2-col_V3_sandwich/20151030_v3_sub30X.RData")
# asm1 <- "/home/users/elam/data/20151027_NA12878_BspQI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
# asm2 <- "/home/users/elam/data/20151027_NA12878_BssSI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
#
# output <- check_assembly_source(indata,asm1,asm2)

#for path df
#quickly check the assembly source of the path df ids
#expectation is that x_ids should correspond to asm1 ids
check_assembly_source1 <- function(paths,asm1_path,asm2_path,mode=1) {

  paths <- lapply(1:length(paths),function(id) {
    path <- paths[[id]]
    if (mode==1) {
      path <- cbind("id"=id,path[[2]])
    } else if (mode==2) {
      path <- cbind("id"=id,path)
    }

    return(path)
  })
  path_df <- do.call(rbind,paths)

  asm1_ids <- unique(readcmap(asm1_path)[,"CMapId"])
  asm2_ids <- unique(readcmap(asm2_path)[,"CMapId"])

  asm1_nodes <- path_df[path_df[,"node_type"]=="asm1","node"]
  asm2_nodes <- path_df[path_df[,"node_type"]=="asm2","node"]

  #fraction of ids are not found
  asm1_nodes_to_asm1 <- sum(!asm1_nodes%in%asm1_ids)/length(asm1_nodes)
  asm1_nodes_to_asm2 <- sum(!asm1_nodes%in%asm2_ids)/length(asm1_nodes)

  asm2_nodes_to_asm1 <- sum(!asm2_nodes%in%asm1_ids)/length(asm2_nodes)
  asm2_nodes_to_asm2 <- sum(!asm2_nodes%in%asm2_ids)/length(asm2_nodes)

  outlist <- list("asm1_nodes_to_asm1"=asm1_nodes_to_asm1,
                  "asm1_nodes_to_asm2"=asm1_nodes_to_asm2,
                  "asm2_nodes_to_asm1"=asm2_nodes_to_asm1,
                  "asm2_nodes_to_asm2"=asm2_nodes_to_asm2)

  if (asm1_nodes_to_asm1+asm2_nodes_to_asm2==0) {
    return(TRUE)
  } else {
    return(FALSE)
  }

  return(outlist)
}

# # load("/home/users/elam/20151030_2-col_V3_sandwich/20151030_v2_sub110X_outlist.RData")
# load("/home/users/elam/20151030_2-col_V3_sandwich/20151030_v3_sub110X_outlist.RData")
# asm1 <- "/home/users/elam/data/20151027_NA12878_BspQI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
# asm2 <- "/home/users/elam/data/20151027_NA12878_BssSI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
#
# output <- check_assembly_source1(outlist[[1]],asm1,asm2)
#
# load("/home/users/elam/20151030_2-col_V3_sandwich/20151030_v3_sub30X_outlist_TEST.RData")
# asm1 <- "/home/users/elam/data/20151027_NA12878_BspQI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
# asm2 <- "/home/users/elam/data/20151027_NA12878_BssSI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap"
#
# output <- check_assembly_source1(outlist[[1]],asm1,asm2)
