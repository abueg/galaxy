#20151026 EL - analyze alignments

# rm(list=ls())

add_offsets <- function(indata) { #getting offsets
  indata[,"orient.x"] <- (indata[, "Orientation.x"]=="+")*2-1
  indata[,"orient.y"] <- (indata[, "Orientation.y"]=="+")*2-1
  indata[,"rel_orient"] <- indata[,"orient.x"]*indata[,"orient.y"]

  #indata[,"offset.x"]<- 0.5*(indata[,"RefStartPos.x"]+indata[,"RefEndPos.x"]-indata[,"RefLen.x"])-0.5*(indata[,"QryStartPos.x"]+indata[,"QryEndPos.x"]-indata[,"QryLen.x"])*indata[,"orient.x"]
  #indata[,"offset.y"]<- 0.5*(indata[,"RefStartPos.y"]+indata[,"RefEndPos.y"]-indata[,"RefLen.y"])-0.5*(indata[,"QryStartPos.y"]+indata[,"QryEndPos.y"]-indata[,"QryLen.y"])*indata[,"orient.y"]

  indata[,"offset.x"]<- 0.5*(indata[,"RefStartPos.x"]+indata[,"RefEndPos.x"])-0.5*(indata[,"QryStartPos.x"]+indata[,"QryEndPos.x"])
  indata[,"offset.y"]<- 0.5*(indata[,"RefStartPos.y"]+indata[,"RefEndPos.y"])-0.5*(indata[,"QryStartPos.y"]+indata[,"QryEndPos.y"])

  indata[,"rel_offset"]<-indata[,"offset.x"]-indata[,"offset.y"]*indata[,"rel_orient"]

  return(indata)
}

get_num_edges <- function(indata) {
  num_edges <- nrow(indata)

  return(num_edges)
}

get_aligned_fraction <- function(indata,thresh) {
  x_span <- abs(indata[,"QryStartPos.x"]-indata[,"QryEndPos.x"])/indata[,"QryLen.x"]
  y_span <- abs(indata[,"QryStartPos.y"]-indata[,"QryEndPos.y"])/indata[,"QryLen.y"]

  min_span <- pmin(x_span,y_span)
  avg_min_span <- mean(min_span)

  return(avg_min_span)

  bool_span <- min_span>thresh
  kept_fraction <- sum(bool_span)/length(bool_span)

  return(kept_fraction)
}

get_correct_orient_fraction <- function(indata) {
  orientations <- ((indata[,"Orientation.x"]=="+")*2-1)*((indata[,"Orientation.y"]=="+")*2-1)
  majority_orientation <- names(sort(-table(orientations)))[1]

  bool_orient <- orientations==majority_orientation

  kept_fraction <- sum(bool_orient)/length(bool_orient)

  return(kept_fraction)
}

get_columns <- function(indata,selected_columns) {
  indata_filtered <- indata[,selected_columns]

  tmp_names <- colnames(indata_filtered) #fix column names by removing the .x or .y suffix
  tmp_names <- lapply(tmp_names,function(tmp_name) {
    unlist(strsplit(tmp_name,split="\\."))[1]
  })

  colnames(indata_filtered) <- unlist(tmp_names)

  return(indata_filtered)
}

calculate_coverage <- function(indata) {
  span <- sum(abs(indata[,"RefStartPos"]-indata[,"RefEndPos"]))

  return(span)
}

get_refid <- function(indata) {
  return(indata[1,"RefcontigID"])
}

get_reflen <- function(indata) {
  return(indata[1,"RefLen"])
}

get_map_coverage <- function(indata) {
  indata_split_x <- split(indata,indata[,"RefcontigID.x"])
  indata_split_x <- lapply(indata_split_x,get_columns,c(1:9,11:12))
#   print("1")
#   print(head(indata_split_x[[1]]))

  indata_split_y <- split(indata,indata[,"RefcontigID.y"])
  indata_split_y <- lapply(indata_split_y,get_columns,c(1,15:22,24:25))
#   print("2")
#   print(head(indata_split_y[[1]]))

  cov_x <- lapply(indata_split_x,calculate_coverage)
  cov_y <- lapply(indata_split_y,calculate_coverage)

  refid_x <- lapply(indata_split_x,get_refid)
  refid_y <- lapply(indata_split_y,get_refid)

  reflen_x <- lapply(indata_split_x,get_reflen)
  reflen_y <- lapply(indata_split_y,get_reflen)

  avg_cov_x <- unlist(cov_x)/unlist(reflen_x)
  avg_cov_y <- unlist(cov_y)/unlist(reflen_y)

  df_x <- cbind("refid"=unlist(refid_x),"cov"=unlist(avg_cov_x),"type"="asm1")
  df_y <- cbind("refid"=unlist(refid_y),"cov"=unlist(avg_cov_y),"type"="asm2")

  print(df_x)
  print(df_y)

  out_df <- rbind(df_x,df_y)

  return(out_df)
}

get_coverage_range <- function(indata,map_coverage) {
  mapid1 <- indata[1,"RefcontigID.x"]
  mapid2 <- indata[1,"RefcontigID.y"]
  cov1 <- as.numeric(map_coverage[map_coverage[,"type"]=="asm1" & map_coverage[,"refid"]==mapid1,"cov"])
  cov2 <- as.numeric(map_coverage[map_coverage[,"type"]=="asm2" & map_coverage[,"refid"]==mapid2,"cov"])

  cov_range <- range(cov1,cov2)

  out_df <-data.frame("low"=cov_range[1],"high"=cov_range[2],stringsAsFactors=F)

  return(out_df)
}

get_avg_confidence <- function(indata) {
  min_conf <- pmin(indata[,"Confidence.x"],indata[,"Confidence.y"])
  avg_min_conf <- mean(min_conf)

  return(avg_min_conf)
}

get_correlation <- function(indata) {

  pair_correlation <- cor(indata[,"offset.x"],indata[,"offset.y"],method="pearson") #default
  pair_correlation1 <- cor(indata[,"offset.x"],indata[,"offset.y"],method="spearman")
  pair_correlation2 <- cor(indata[,"offset.x"],indata[,"offset.y"],method="kendall")

  r_sq <- pair_correlation^2
  r_sq1 <- pair_correlation1^2
  r_sq2 <- pair_correlation2^2

  r_sq <- max(r_sq,r_sq1,r_sq2)

  if (is.na(r_sq)) { #if NA, turn into 0
    r_sq <- 0
  }

  return(r_sq)
}

check_overlap_ratio <- function(start.a, end.a, start.ref, end.ref) {
  start <- max(min(start.a,end.ref),start.ref)
  end <- max(min(end.a,end.ref),start.ref)

  if ((end-start)==0) { #no overlap
    return(0)
  } else { #overlap ratio
    return((end-start)/(end.ref-start.ref))
  }
}

one_third_coverage_simplified <- function(indata_start_end_len) {
  have_cov <- 0

  molecule_start <- min(indata_start_end_len[,1])
  molecule_end <- max(indata_start_end_len[,2])

  third_index <- round(indata_start_end_len[1,3]/3)

  first_third <- check_overlap_ratio(molecule_start,molecule_end,1,third_index)
  second_third <- check_overlap_ratio(molecule_start,molecule_end,third_index+1,third_index*2)
  third_third <- check_overlap_ratio(molecule_start,molecule_end,(third_index*2)+1,indata_start_end_len[1,3])

  bin_cov <- c(first_third,second_third,third_third)
  have_cov <- as.numeric(bin_cov>0.5)

  if (sum(have_cov)==0) {
    have_cov <- as.numeric(bin_cov>0)
  }

  return(have_cov)
}

get_configuration <- function(indata) {
  x_summary <- one_third_coverage_simplified(indata[,c("RefStartPos.x","RefEndPos.x","RefLen.x")])
  y_summary <- one_third_coverage_simplified(indata[,c("RefStartPos.y","RefEndPos.y","RefLen.y")])

  if (length(which(x_summary==0)>0) && length(which(y_summary==0)>0)) {
    orientation <- median(((indata[,"Orientation.x"]=="+")*2-1)*((indata[,"Orientation.y"]=="+")*2-1))

    if (orientation==-1) {
      y_summary <- rev(y_summary)
    }

    xy_summary <- x_summary + y_summary

    if ((xy_summary[1]==0)||(xy_summary[3]==0)) {
      return(1)
    } else {
      return(0)
    }
  } else {
    return(0)
  }
}

get_slope <- function(indata) {
  slope <- NA

  if (nrow(indata)) {
    return(slope)
  }

  try({
    rlm_out <- rlm(indata[,"offset.y"]~indata[,"offset.x"])
    slope <- coef(summary(rlm_out))[2,1]
  })

  return(slope)
}

analyze_alignment <-function(indata) {
  id_str <- indata[1,"id_str"]

  num_edges <- get_num_edges(indata)

  aligned_fraction <- get_aligned_fraction(indata,0.5)

  correct_orient_fraction <- get_correct_orient_fraction(indata)

  map_coverage <- get_map_coverage(indata)
  map_coverage_range <- get_coverage_range(indata,map_coverage)

  avg_conf <- get_avg_confidence(indata)

  pair_correlation <- get_correlation(indata)

  configuration_notRight <- get_configuration(indata)

  out_df <- data.frame("id_str"=id_str,
                       "num_edges"=num_edges,
                       "aligned_fraction"=aligned_fraction,
                       "correct_orient_fraction"=correct_orient_fraction,
                       "map_coverage_range"=map_coverage_range,
                       "avg_conf"=avg_conf,
                       "pair_correlation"=pair_correlation,
                       "configuration"=configuration_notRight)

  return(out_df)
}

analyze_alignment1 <- function(indata) {
  id_str <- indata[1,"id_str"]

  num_edges <- get_num_edges(indata)
  pair_correlation <- get_correlation(indata)
  slopes <- get_slope(indata)

  out_df <- data.frame("id_str"=id_str,
                       "num_edges"=num_edges,
                       "pair_correlation"=pair_correlation,
                       "slopes"=slopes)

  return(out_df)
}

analyze_alignments <- function(indata,parallel_flag=F,cores=2,mode="simple") {
  indata <- add_offsets(indata)

  indata_split <- split(indata,indata[,"id_str"])

  if (mode=="simple") {
    if (parallel_flag) {
      output <- mclapply(indata_split,analyze_alignment1,mc.cores=cores)
    } else {
      output <- lapply(indata_split,analyze_alignment1)
    }
    output <- do.call(rbind,output)
  } else {
    if (parallel_flag) {
      output <- mclapply(indata_split,analyze_alignment,mc.cores=cores)
    } else {
      output <- lapply(indata_split,analyze_alignment)
    }
    output <- do.call(rbind,output)
  }

  return(output)
}

make_alignment_plots <- function(indata,plot_flag=F) {

  out <- c(0,0,0)

  x <- indata[,c("QryContigID","RefLen.x","RefStartPos.x","RefEndPos.x"),F]
  names(x) <- c("QryContigID","RefLen","RefStartPos","RefEndPos")
  y <- indata[,c("QryContigID","RefLen.y","RefStartPos.y","RefEndPos.y"),F]
  names(y) <- c("QryContigID","RefLen","RefStartPos","RefEndPos")

  #if(nrow(x) == 1){
      #x <- rbind(x,x)n
      #y <- rbind(y,y)
  #}

  id_str <- indata[1,"id_str"]
  ids <- unlist(strsplit(id_str,"_"))

  try({
    #plotting
    cov_x <- plot_coverage(x,paste0(id_str," (",ids[1],")"),plot_flag=plot_flag) #plot 1
    cov_y <- plot_coverage(y,paste0(id_str," (",ids[2],")"),plot_flag=plot_flag) #plot 2
    rlm_out <- plot_correlation1(indata,plot_flag=plot_flag) #plot 3
    plot_alignment1(indata,plot_flag=plot_flag) #plot 4

    out <- c(id_str,rlm_out)
  })

  return(out)
}

#plot alignments and output offsets
plot_alignments <- function(indata_split,plot_flag=F,parallel_flag=F,cores=2) {
  if (plot_flag) { #if plot_flag is on, require that parallel_flag be off
    parallel_flag <- F
  }

  #if plotting, single-core
  if (plot_flag) { #single-thread if plot_flag is true
    par(mfrow=c(2,2))
    offset_list <- lapply(indata_split,make_alignment_plots,plot_flag=plot_flag)
  } else { #if not plotting, either single- or multi-core
    if (parallel_flag) {
      offset_list <- mclapply(indata_split,make_alignment_plots,plot_flag=plot_flag,mc.cores=cores)
    } else {
      offset_list <- lapply(indata_split,make_alignment_plots,plot_flag=plot_flag)
    }
  }
  names(offset_list) <- names(indata_split)

  return(offset_list)
}

# #main
# # full_data <- "/home/users/elam/20150911_2-col_actual_data/sop_and_ml_test/Human_willsample_twocolor_combined_indata.RData"
# # load(full_data)
# #
# # indata <- indata[sample(nrow(indata),nrow(indata)*0.2),]
# partial_data <- "/home/users/elam/20150911_2-col_actual_data/sop_and_ml_test/Human_willsample_twocolor_combined_indata_subset.RData"
# # save(list="indata",file=partial_data)
# load(partial_data)
#
# parallel_flag <- T
# cores <- 64
# output <- analyze_alignments(indata,parallel_flag,cores)
#
# stop()
#
# #TEST
# #testing coverage calculation
# indata_split <- split(indata,indata[,"id_str"])
# test <- indata_split[["1000_128"]]
# test1 <- test[,c("RefStartPos.x","RefEndPos.x","RefStartPos.y","RefEndPos.y")]
# test1[,"diffx"] <- test1[,"RefEndPos.x"]-test1[,"RefStartPos.x"]
# test1[,"diffy"] <- test1[,"RefEndPos.y"]-test1[,"RefStartPos.y"]
#
# c(sum(test1[,"diffx"]),sum(test1[,"diffy"]))/test[1,c("RefLen.x","RefLen.y")]
