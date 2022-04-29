#20151012 EL - check and get input alignment

#20160420 This function assumes reference is assembly and query is two-color glues
#In hyrid-scaffold, the alignment is done in reverse order, twocolor.as.ref "flip" the alignment
#so it can be handled properly in downstream analysis
get_alignment_data <- function(align_dir,xmap1,xmap2,parallel_flag=T,cores=2, twocolor.as.ref=F,  score.T = 13) { #get xmap data and generated a merged data frame

  xmap1 <- paste0(align_dir,xmap1)
  xmap2 <- paste0(align_dir,xmap2)

  xmap1_exists <- file.exists(xmap1) #check if the files exist
  xmap2_exists <- file.exists(xmap2)

  if (!xmap1_exists) {
    stop("xmap1 is not found")
  } else if (!xmap2_exists) {
    stop("xmap2 is not found")
  }

  if (parallel_flag) { #read xmap in parallel
    xmaps <- mclapply(c(xmap1,xmap2),readxmap,mc.cores=cores)
  } else {
    xmaps <- lapply(c(xmap1,xmap2),readxmap)
  }
  if(twocolor.as.ref){
      print("two.color mols is reference, flipping alignment")
      xmap1 <- flip.xmap(xmaps[[1]])
      xmap2 <- flip.xmap(xmaps[[2]])
  }else{
      xmap1 <- xmaps[[1]]
      xmap2 <- xmaps[[2]]
  }
  xmap1 <- xmap1[xmap1$Confidence > score.T,]
  xmap2 <- xmap2[xmap2$Confidence > score.T,]

  #xmap1 <- rbind(xmap1, xmap1)
  #xmap2 <- rbind(xmap2, xmap2)


  merged <- merge(xmap1,xmap2,by="QryContigID",all=FALSE) #merge to find molecules that map to both assemblies
  if(nrow(merged) > 0){
    merged[,"id_str"] <- paste0(merged[,"RefcontigID.x"],"_",merged[,"RefcontigID.y"]) #get str representing pair
  }
  return(merged)
}

get_input_alignment <- function(input,input_type,save_flag=F, twocolor.as.ref=F, score.T) {
  if (input_type=="rdata") { #expects rdata

    if (class(input)=="list") {
      align_dir <- input[["align_dir"]]
      rdata <- input[["rdata"]]

      input <- paste0(align_dir,rdata)
    }

    if (file.exists(input)) {
      load(input) #has to be named indata

      if (!exists("indata")) {
        stop("indata not found after importing")
      }

    } else {
      stop("rdata file not found")
    }

  } else if (input_type=="xmap") { #expects two xmap files

    align_dir <- input[["align_dir"]]
    xmap1 <- input[["xmap1"]]
    xmap2 <- input[["xmap2"]]

    indata <- get_alignment_data(align_dir=align_dir,xmap1=xmap1,xmap2=xmap2, score.T = score.T) #load alignment data

  } else if (input_type=="prefix") { #expects the bnx/xmap prefix
    align_dir <- input[["align_dir"]]
    prefix <- input[["prefix"]]
    load_previous_flag <- input[["load_previous_flag"]]

    if (is.null(load_previous_flag)) { #if the flag is not specified, enable loading
      load_previous_flag <- T
    }

    if (load_previous_flag) { #load previous data
      load(file=paste0(align_dir,prefix,"_indata.RData"))

      if (!exists("indata")) {
        stop("indata not found after importing")
      }
    } else {
      xmap1 <- paste0(prefix,"_1_fullAln.xmap")
      xmap2 <- paste0(prefix,"_2_fullAln.xmap")

      indata <- get_alignment_data(align_dir=align_dir,xmap1=xmap1,xmap2=xmap2, twocolor.as.ref=twocolor.as.ref, score.T = score.T) #load alignment data

      if (save_flag) {
        rdata_file <- paste0(align_dir,prefix,"_indata.RData")
        if (!file.exists(rdata_file)) {
          save(list="indata",file=rdata_file)
        } else {
          stop(paste0(rdata_file," already exists."))
        }
      }
    }

  } else {
    stop("input type not recognized")
  }

  return(indata)
}
