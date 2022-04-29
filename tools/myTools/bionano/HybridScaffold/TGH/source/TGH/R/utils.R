#20151026 EL - utility functions

# rm(list=ls())
#
# #general utility scripts
# scripts_dir <- "/home/users/elam/rscripts/"
# source(paste0(scripts_dir,"readmaps.R")) #util functions to read _map formats

#check which number is closer to 1 regardless of signs
get_closer_to_one <- function(a,b) {
  if (is.na(a) && is.na(b)) {
    return(NA)
  } else if (is.na(a)) {
    return(b)
  } else if (is.na(b)) {
    return(a)
  }

  #if (abs(abs(a)-1)>abs(abs(b)-1)) {
  if (abs(abs(a)-1)<abs(abs(b)-1)) {
    return(a)
  } else {
    return(b)
  }
}

#util function to make a value negative
turn_negative <- function(value) {
  return(-1*abs(value))
}

#select alignments from one chromosome
select_alignement_chr <- function(xmap,chr=22){
  xmap <- readxmap(xmap)
  xmap_selected <- xmap[xmap[,"RefcontigID"]==chr,]

  return(xmap_selected)
}

#basic function to calculation mapping rate based on an .err df
get_mapping_rate <- function(err) {
  mapping_rate <- err[,"aligned"]/err[,"maps"]

  return(mapping_rate)
}

#read err files (one for each color) and get rates
analyze_rates <- function(prefix,indir) {

  rate1 <- NA
  rate2 <- NA
  effective_cov <- NA

  try({
    err1 <- paste0(indir,prefix,"_1_fullAln.err") #mapping rates for individual colors
    err2 <- paste0(indir,prefix,"_2_fullAln.err")

    if (file.exists(err1)) {
      err1 <- readerr(err1)
      rate1 <- get_mapping_rate(err1)
    }

    if (file.exists(err2)) {
      err2 <- readerr(err2)
      rate2 <- get_mapping_rate(err2)
    }

    xmap1 <- paste0(indir,prefix,"_1_fullAln.xmap")
    xmap2 <- paste0(indir,prefix,"_2_fullAln.xmap")

    if (file.exists(xmap1) && file.exists(xmap1)) {
      xmaps <- mclapply(c(xmap1,xmap2),readxmap,mc.cores=2)
      xmap1 <- xmaps[[1]]
      xmap2 <- xmaps[[2]]

      merged <- merge(xmap1,xmap2,by="QryContigID",all=FALSE)
      sum_usable <- sum(merged[,"QryLen.x"])
      effective_cov <- sum_usable/3200000000
    }
  })

  out_df <- data.frame("ch1_maprate"=rate1,
                       "ch2_maprate"=rate2,
                       "eff_cov"=effective_cov)

  return(out_df)
}

#get mapping rates for the two colors
#output: df containing the rates
get_alignment_rates <- function(prefixes,indir,repeat_flag=T,save_flag=F,cores=32) {
  infile <- paste0("mapping_rates.RData")

  if (!repeat_flag) { #if not repeating analysis, find RData file
    if (!file.exists(infile)) {
      stop(paste0(infile," not found."))
    }

    load(infile)
  } else { #if repeating analysis, make sure RData not already created
    if (save_flag) {
      if (file.exists(infile)) {
        stop(paste0(infile," already exists."))
      }
    }

    rates <- mclapply(prefixes,analyze_rates,indir,mc.cores=cores) #get mapping rates
    rates <- do.call(rbind,rates)
    rates <- cbind("prefixes"=prefixes,rates)

    if (save_flag) {
      save(list="rates",file=infile)
    }
  }

  return(rates)
}

#givien two vectors, make them the same length
#if two vectors have different lenghts, add dummy numbers to the end of the shorter one
make_same_length <- function(v1,v2,dummy=0) {
  len1 <- length(v1)
  len2 <- length(v2)
  len_diff <- len1-len2

  if (len1>len2) {
    v2 <- c(v2,rep(dummy,abs(len_diff)))
  } else {
    v1 <- c(v1,rep(dummy,abs(len_diff)))
  }

  stopifnot(length(v1)==length(v2))

  outdf <- cbind(v1,v2)

  return(outdf)
}

#normalize two dataframes by keeping only the common columns
normalize_df <- function(df1,df2,run_flag=T) {

  diff_flag <- ncol(df1)==ncol(df2) #checking if they have the same number of columns

  if (run_flag) {
    common_columns <- intersect(colnames(df1),colnames(df2)) #columns common between df input

    df1_trimmed <- df1[,colnames(df1)%in%common_columns]
    df2_trimmed <- df2[,colnames(df2)%in%common_columns]

  } else {
    return(diff_flag)
  }

  outlist <- list("df1_trimmed"=df1_trimmed,
                  "df2_trimmed"=df2_trimmed)

  return(outlist)
}


rlm_internal <- MASS::rlm
rlm<- function(...){
  #make sure there is no singular values for rlm
  e <- environment(...)
  dt <- as.data.table(e$data)
  setkeyv(dt, c('offset.x', 'offset.y'))
  dt <- unique(dt)
  dat.unique <- e$data[e$data$QryContigID %in% dt$QryContigID,]
  if(nrow(dat.unique) != nrow(dt)){
    #browser()
  }
  #e$data <- dat.unique
  rlm_out <- tryCatch(
    {
      rlm_internal(...)
    },
    error=function(err){
      NULL
    })
  if(is.null(rlm_out)){
    x= c(1, 1000)
    d <- data.frame(offset.x=x,
                    offset.y=x*e$data$rel_orient[[1]]+get_offset_single(e$data[1,]))
    #print(d)
    rlm_out <- rlm_internal(d[, 'offset.x'] ~ d[, 'offset.y'])
  }
  return(rlm_out)
}


# #20151103 EL - testing
# bspqi_asm <- "/home/users/elam/data/20151031_NA12878_50X_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap" #V2
# bsssi_asm <- "/home/users/elam/data/20151027_NA12878_BssSI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap" #V3
#
# asm1 <- readcmap(bspqi_asm)
# asm2 <- readcmap(bsssi_asm)
#
# asm_normalized <- normalize_df(asm1,asm2)
