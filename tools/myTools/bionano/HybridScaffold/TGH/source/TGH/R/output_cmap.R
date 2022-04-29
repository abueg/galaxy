#20150923 EL - some refactoring
#this part deals with preparing the cmap df and outputting cmap files

# rm(list=ls())
#
# #general scripts
# scripts_dir <- "/home/users/elam/rscripts/"
# source(paste(scripts_dir,"paste0.R",sep=""))
# source(paste0(scripts_dir,"readmaps.R")) #util functions to read _map formats
# source(paste0(scripts_dir,"readcmap_w_headers.R"))
#
# #2-color specific scripts
# scripts_dir <- "/home/users/elam/rscripts/2col_dev/current/"
# source(paste0(scripts_dir,"make_plots.R"))
# source(paste0(scripts_dir,"run_ra.R"))
# source(paste0(scripts_dir,"do_validation.R"))
# source(paste0(scripts_dir,"utils.R"))
#
# #global
# ra <- "/home/users/elam/RefAligner/ra_151016/bin/RefAligner"

#get assembly cmap data
get_assemblies <- function(asm1_filename,asm2_filename,trim_columns='both') {

  #read cmap
  asm1_cmap <- readcmap_w_headers(asm1_filename)
  asm2_cmap <- readcmap_w_headers(asm2_filename)

  #trim if necessary
  if (!is.null(trim_columns)) { #if doing any trimming, get common columns
    common_columns <- intersect(colnames(asm1_cmap[[2]]),colnames(asm2_cmap[[2]])) #columns common between the two cmap input
  }

  if (trim_columns=='asm1') { #trim asm1
    asm1_cmap[[2]] <- asm1_cmap[[2]][,colnames(asm1_cmap[[2]])%in%common_columns]
  } else if (trim_columns=='asm2') { #trim asm2
    asm2_cmap[[2]] <- asm2_cmap[[2]][,colnames(asm2_cmap[[2]])%in%common_columns]
  } else if (trim_columns=='both') { #trim both
    asm1_cmap[[2]] <- asm1_cmap[[2]][,colnames(asm1_cmap[[2]])%in%common_columns]
    asm2_cmap[[2]] <- asm2_cmap[[2]][,colnames(asm2_cmap[[2]])%in%common_columns]
  }

  asm1_cmap_split <- split(asm1_cmap[[2]],asm1_cmap[[2]][,"CMapId"]) #split by cmapid
  asm1_cmap_split <- lapply(asm1_cmap_split,head,-1)
  asm2_cmap_split <- split(asm2_cmap[[2]],asm2_cmap[[2]][,"CMapId"]) #split by cmapid
  asm2_cmap_split <- lapply(asm2_cmap_split,head,-1)

  outlist <- list("asm1_filename"=asm1_filename,
                  "asm2_filename"=asm2_filename,
                  "asm1_cmap"=asm1_cmap_split,
                  "asm2_cmap"=asm2_cmap_split)

  return(outlist)
}


#flip dataframe
flip <- function(data) { #util; flip a dataframe 180 degrees
  stopifnot(class(data)=="data.frame") #this function does not work for matrices
  inverted <- data[rev(rownames(data)), ]
  #rownames(inverted) <- NULL

  return(inverted)
}

orient_map <- function(cmap,orientation) { #flip a map based on orientation
  if (orientation==1) { #don't flip if orientation is 1
    return(cmap)
  }

  inverted <- flip(cmap) #rotate dataframe by 180 degrees

  map_len <- inverted[1,"ContigLength"]
  inverted[,"Position"] <- map_len - inverted[,"Position"] #get new positions; the rest (stdev, SNR etc) are kept unchanged

  return(inverted)
}

shift_map <- function(cmap,shift) { #shift map by a certain amount
  if (shift==0) {
    return(cmap)
  }

  cmap[,"Position"] <- cmap[,"Position"] + shift

  return(cmap)
}

cut_map <- function(cmap,left_cut=NULL,right_cut=NULL) {
  if (nrow(cmap)==0) {
    return(cmap)
  }

  if (is.null(left_cut) && is.null(right_cut)) {
    return(cmap)
  }

  if (!is.null(left_cut)) { #if the left needs to be cut; generates cmap_cut for next step
    cmap_cut <- cmap[cmap[,"Position"]>left_cut,]
  } else {
    cmap_cut <- cmap #unmodified
  }

  if (!is.null(right_cut)) { #if the right needs to be cut
    cmap_cut <- cmap_cut[cmap_cut[,"Position"]<right_cut,]
  }

  return(cmap_cut)
}

fix_channel <- function(cmap,channel) { #remove backbone channel and then fix label channel id
  stopifnot(sum(cmap[,"LabelChannel"]==0)==0) #make sure there is no backbone line

  if (nrow(cmap)==0) {
    return(cmap)
  }

  cmap_filtered <- cmap[cmap[,"LabelChannel"]!=0,]
  cmap_filtered[,"LabelChannel"] <- channel

  return(cmap_filtered)
}

fix_cmap <- function(cmap,cmapid=1,id_flag=F) { #fix various things in merged cmap
  stopifnot(sum(cmap[,"LabelChannel"]==0)==0) #make sure there is no backbone line

  #fix cmapid if id_flag on - make cmapid 1 by default
  if (id_flag) {
    cmap[,"CMapId"] <- cmapid
  }

  #fix contig length
  max_pos <- max(cmap[,"Position"])
  cmap[,"ContigLength"] <- max_pos

  #fix num sites
  nrows <- nrow(cmap)
  cmap[,"NumSites"] <- nrows-1

  #add backbond line
  cmap[nrow(cmap),"LabelChannel"] <- 0

  #fix site ids
  cmap[,"SiteID"] <- 1:nrow(cmap)

  #fix the last columns
  if (("GmeanSNR")%in%rownames(cmap)) {
    cmap[nrow(cmap),c("StdDev","Coverage","Occurrence","GmeanSNR","lnSNRsd","SNR")] <- c(0,1,1,0,0,0)
  } else {
    cmap[nrow(cmap),c("StdDev","Coverage","Occurrence")] <- c(0,1,1)
  }

  return(cmap)
}

merge_cmap_df <- function(cmap_list1,cmap_list2) {
  cmap_list1 <- lapply(cmap_list1,fix_channel,1) #fix label channel
  cmap_df1 <- do.call(rbind,cmap_list1)

  cmap_list2 <- lapply(cmap_list2,fix_channel,2)
  cmap_df2 <- do.call(rbind,cmap_list2)

  cmap_merged <- rbind(cmap_df1,cmap_df2)
  cmap_merged <- cmap_merged[order(cmap_merged[,"Position"]),]

  cmap_fixed <- fix_cmap(cmap_merged,id_flag=T)

  return(cmap_fixed)
}

merge_1col_cmap_df <- function(cmap_list) {
  cmap_df <- do.call(rbind,cmap_list)
  stopifnot(sum(cmap_df[,"LabelChannel"]==0)==0) #make sure there is no backbone line

  cmap_df <- cmap_df[order(cmap_df[,"Position"]),]

  cmap_fixed <- fix_cmap(cmap_df,id_flag=T)

  return(cmap_fixed)
}

output_cmap_df <- function(path_coordinates,cmap_split) {
  stopifnot(nrow(path_coordinates)>0)

  starts <- pmin(path_coordinates[,"rel_start"],path_coordinates[,"rel_end"])
  path_coordinates <- path_coordinates[order(starts),]

#   print(path_coordinates)

  #current_path <- matrix(data=NA,nrow=nrow(path_coordinates),ncol=2) #initialize empty matrix
  current_path <- path_coordinates #just make a copy

#   print("initial path")
#   print(current_path)

  last_end <- -1
  cmap_list <- list()

  for (i in 1:nrow(path_coordinates)) {
#     print(i)

    entry <- path_coordinates[i,,drop=F]

#     print("entry")
#     print(entry)

#     if (i==1) { #if first entry, just take current values
#       current_path[1,] <- c(entry[,"rel_start"],entry[,"rel_end"])
#     }

    current_start <- min(entry[,"rel_start"],entry[,"rel_end"])

#     print("current_start")
#     print(current_start)

    cmap_selected <- cmap_split[[as.character(entry[,"node"])]]

#     print("cmap_selected")
#     print(head(cmap_selected,5))

    orientation <- entry[,"orientation"]
    shift <- entry[,"rel_start"]

    if (orientation==1) {
      shift <- entry[,"rel_start"]
    } else {
      shift <- entry[,"rel_end"]
    }

#     print("shift")
#     print(shift)

    cmap_oriented <- orient_map(cmap_selected,orientation) #orient map if necessary
#     print("cmap_oriented")
#     print(head(cmap_oriented,5))
    cmap_shifted <- shift_map(cmap_oriented,shift) #flip map if necessary
#     print("cmap_shifted")
#     print(head(cmap_shifted,5))

    offset <- current_start-last_end

    if (i==1) {
      cmap_cut <- cut_map(cmap_shifted,left_cut=NULL,right_cut=NULL) #don't need cutting for 1st entry
    } else {
      if (offset>0) { #no overlap between maps
#         print("here")
        cmap_cut <- cut_map(cmap_shifted,left_cut=NULL,right_cut=NULL)
#         current_path[i,c("rel_start","rel_end")] <- c(entry[,"rel_start"],entry[,"rel_end"])
      } else { #cut 2nd map if there is overlap
#         print("here1")
        cmap_cut <- cut_map(cmap_shifted,left_cut=last_end,right_cut=NULL)
        if (orientation==1) {
          current_path[i,c("rel_start","rel_end")] <- c(last_end,entry[,"rel_end"])
        } else {
          current_path[i,c("rel_start","rel_end")] <- c(entry[,"rel_start"],last_end)
        }
      }
    }

#     print("cmap_cut")
#     print(head(cmap_cut))

    #     if (nrow(cmap_cut)==0) {
    #       cmap_cut <-NULL
    #     }

    last_end <-  max(entry[,"rel_start"],entry[,"rel_end"])

    cmap_list[[i]] <- cmap_cut
  }

#   print("trimmed path")
#   print(current_path)

#   print(head(lapply(cmap_list,head)))

  outlist <- list("cmap_path"=current_path,
                  "cmap_list"=cmap_list)

  return(outlist)
  return(cmap_list)
}

check_ranges <- function(df) {
  mins <- pmin(df[,1],df[,2])
  maxes <- pmax(df[,1],df[,2])

#   print("mins")
  mins <- mins[-1] #remove first element
#   print(mins)

#   print("maxes")
  maxes <- head(maxes,-1) #all elements except last 1
#   print(maxes)

  stopifnot(length(mins)==length(maxes))

  diffs <- mins-maxes
#   print("diffs")
#   print(diffs)

  return(sum(diffs<=0)==0)
}

#check to make sure the cmap trimming was done correctly
#returns a bool
check_cmap_list <- function(cmap_list) {
  asm1 <- cmap_list[["asm1"]]
  asm2 <- cmap_list[["asm2"]]

  range_asm1 <- do.call(rbind,(lapply(asm1,function(x) range(x[,"Position"]))))
  range_asm2 <- do.call(rbind,(lapply(asm2,function(x) range(x[,"Position"]))))

#   print(range_asm1)
  bool_check_ranges1 <- check_ranges(range_asm1)
#   print(bool_check_ranges1)

#   print(range_asm2)
  bool_check_ranges2 <- check_ranges(range_asm2)
#   print(bool_check_ranges2)

  if (bool_check_ranges1 && bool_check_ranges2) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#write cmap to file (assume cmap has been fixed)
#09/19/2016
#add headers to cmap output
write_cmap <- function(input,mode,renumber_flag=F) {
  if (mode=="single") { #there is only a single cmap to be output
    filename <- input[["filename"]]
    cmap <- input[["cmap"]]

    #write.table(cmap,file=filename,quote=F,sep="\t",row.names=F,col.names=F)
    write.cmap.header(filename, cmap)
    write.table(format(cmap,digits=1, scientific=F),file=filename,quote=F,sep="\t",row.names=F,col.names=F, append=T)

    id_list <- filename
  } else if (mode=="multi") { #there is a list of cmap to be output
    prefix <- input[["prefix"]]
    cmap_list <- input[["cmap_list"]]

    if (renumber_flag) { #if the cmap ids have to be renumbered
      lapply(1:length(cmap_list),function(i) {
        cmap <- cmap_list[[i]]
        filename <- paste0(prefix,i,".cmap")
        #write.table(cmap,file=filename,quote=F,sep="\t",row.names=F,col.names=F)
        write.cmap.header(filename, cmap)
        write.table(format(cmap,digits=1, scientific=F),file=filename,quote=F,sep="\t",row.names=F,col.names=F, append=T)
      })

      id_list <- paste0(prefix,"_id_list.txt")
      write.table(paste0(prefix,1:length(cmap_list),".cmap"), #prepare id list file
                  file=id_list,
                  quote=F,sep="\t",row.names=F,col.names=F)
    } else {
      cmapids <- unlist(lapply(cmap_list,function(cmap) cmap[1,"CMapId"]))

      lapply(1:length(cmap_list),function(i) {
        cmap <- cmap_list[[i]]
        cmapid <- cmapids[[i]]
        filename <- paste0(prefix,cmapid,".cmap")
        #write.table(cmap,file=filename,quote=F,sep="\t",row.names=F,col.names=F)
        write.cmap.header(filename, cmap)
        write.table(format(cmap,digits=1),file=filename,quote=F,sep="\t",row.names=F,col.names=F, append=T)
      })

      id_list <- paste0(prefix,"_id_list.txt")
      write.table(paste0(prefix,cmapids,".cmap"), #prepare id list file
                  file=id_list,
                  quote=F,sep="\t",row.names=F,col.names=F)
    }
  } else {
    stop("write_cmap mode not recognized")
  }
  return(id_list)
}

#add by Jian 09/19/2016
write.cmap.header <- function(file.name, cmap){
    num.maps <- nrow(subset(cmap, SiteID==1))
    labelChannels <- length(unique(subset(cmap, LabelChannel > 0, LabelChannel)[[1]]))
    sink(file.name)
    cat( "# CMAP File Version:\t0.1\n" );
    cat( paste0("# Label Channels:\t", labelChannels, "\n" ));
    cat( "# Nickase Recognition Site 1:\tunknown\n" );
    cat( paste0("# Number of Consensus Nanomaps:\t", num.maps, "\n"));
    header <- colnames(cmap)
    header <- paste("#h", paste(header, collapse = "\t"), collapse="")
    types <- rep('float', ncol(cmap))
    types[unlist(lapply(cmap, is.character))] <- 'string'
    types[unlist(lapply(cmap, is.integer))] <- 'int'
    mask_ind <- which(grepl('Mask', colnames(cmap)))
    if(length(mask_ind) > 0){
      types[mask_ind] <- "Hex"
    }
    header2 <- paste("#f", paste(types, collapse="\t"), collapse="")
    #cat( paste0( "#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\t",
    #    "Position\tStdDev\tCoverage\tOccurrence\tChimQuality\tSegDupL\tSegDupR\tFragileL\tFragileR\tOutlierFrac\tChimNorm\tMask\n" ) );
    #cat( "#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\tfloat\tfloat\tfloat\tfloat\float\tfloat\tfloat\tHex\n" );
    cat(header)
    cat("\n")
    cat(header2)
    cat("\n")
    sink();
}

process_alignment <- function(output_ra_alignref,plot_flag=T) {
  prefix <- output_ra_alignref[["prefix"]]
  xmap <- readxmap(paste0(prefix,".xmap"))

  nrows <- nrow(xmap)

  orientations <- xmap[,"Orientation"]
  n_wrong_orientation <- sum(orientations!="+")

  n_check_ranges <- check_ranges(xmap[,c("RefStartPos","RefEndPos")])

  plot_alignment_schematic(xmap,plot_flag=T)
  plot_alignment_schematic_w_mg(xmap,plot_flag=T)

  outlist <- list("n_alignments"=nrows,
                  "n_check_ranges"=n_check_ranges,
                  "n_wrong_orientation"=n_wrong_orientation)

  return(outlist)
}

test_intermediate_cmap <- function(cmap_list) {
  asm1 <- cmap_list[["asm1"]]
  asm2 <- cmap_list[["asm2"]]

  asm1_fixed <- (lapply(asm1,fix_cmap))
  asm2_fixed <- (lapply(asm2,fix_cmap))

  #prepare cmap files
  #asm1
  prefix_asm1 <- "/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/20151015_test_intermediate_asm1_contig"
  input_asm1 <- list("prefix"=prefix_asm1,
                     "cmap_list"=asm1_fixed)
  id_list_asm1 <- write_cmap(input_asm1,"multi")

  #asm2
  prefix_asm2 <- "/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/20151015_test_intermediate_asm2_contig"
  input_asm2 <- list("prefix"=prefix_asm2,
                     "cmap_list"=asm2_fixed)
  id_list_asm2 <- write_cmap(input_asm2,"multi")

  #run merge using RefAligner
  output_ra_merge1 <- run_ra_merge(id_list_asm1,prefix_asm1,run_flag=F)
  output_ra_merge2 <- run_ra_merge(id_list_asm2,prefix_asm2,run_flag=F)

  #align maps to reference
  #asm1
  ref_asm1 <- "/home/users/csecol/genomes/human/hg19/hg19_BSSSI_0Kb_0labels.cmap" #bsssi reference
  prefix1_asm1 <- paste0(prefix_asm1,"_alignref_final")
  cmap_asm1 <- output_ra_merge1[["cmap"]]
  output_ra_alignref1 <- run_ra_alignref_final(cmap_asm1,ref_asm1,prefix1_asm1,run_flag=F)

  ref_asm2 <- "/home/users/csecol/genomes/human/hg19/hg19_BSPQI_0Kb_0labels.cmap" #bspqi reference
  prefix1_asm2 <- paste0(prefix_asm2,"_alignref_final")
  cmap_asm2 <- output_ra_merge2[["cmap"]]
  output_ra_alignref2 <- run_ra_alignref_final(cmap_asm2,ref_asm2,prefix1_asm2,run_flag=F)

  alignment_analysis <- process_alignment(output_ra_alignref1,plot_flag=T)
#   print(alignment_analysis)

  title_str <- "alignment schematic"
  plot_alignment_schematic_w_mg_2col(output_ra_alignref1,output_ra_alignref2,title_str,plot_flag=T)

  #check one color at a time
  #write cmap
  cmap_df_asm1 <- merge_1col_cmap_df(asm1)
  prefix_asm1 <- "/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/20151016_path1_asm1"
  input_asm1 <- list("filename"=paste0(prefix_asm1,".cmap"),
                     "cmap"=cmap_df_asm1)
  cmap_file1 <- write_cmap(input_asm1,"single")

  cmap_df_asm2 <- merge_1col_cmap_df(asm2)
  prefix_asm2 <- "/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/20151016_path1_asm2"
  input_asm2 <- list("filename"=paste0(prefix_asm2,".cmap"),
                     "cmap"=cmap_df_asm2)
  cmap_file2 <- write_cmap(input_asm2,"single")

  #align to reference
  prefix1_asm1 <- paste0(prefix_asm1,"_alignref_final")
  output_ra_alignref1 <- run_ra_alignref_final(cmap_file1,ref_asm1,prefix1_asm1,run_flag=F)

  prefix1_asm2 <- paste0(prefix_asm2,"_alignref_final")
  output_ra_alignref2 <- run_ra_alignref_final(cmap_file2,ref_asm2,prefix1_asm2,run_flag=F)

  outlist <- list("cmap1"=cmap_df_asm1,
                  "cmap2"=cmap_df_asm2)

  return(outlist)
}

#comparing the single color cmaps with the merged 2 color cmap
compare_cmaps <- function(single_col_cmaps,cmap_merged) {
  cmap1 <- single_col_cmaps[["cmap1"]]
  cmap2 <- single_col_cmaps[["cmap2"]]

  merged_ch1 <- cmap_merged[cmap_merged[,"LabelChannel"]==1,]
  merged_ch2 <- cmap_merged[cmap_merged[,"LabelChannel"]==2,]

  row_diff1 <- nrow(cmap1)-nrow(merged_ch1)
  row_diff2 <- nrow(cmap2)-nrow(merged_ch2)

  if (row_diff1+row_diff2>1) {
    stop("the single color cmaps diff from the merged 2 color cmap by too much in terms of the number of labels")
  }

  #check correlation in positions
  df <- make_same_length(cmap1[,"Position"],merged_ch1[,"Position"])
  plot(df[,1],df[,2],pch=16,main="single-color vs 2-color (channel 1)",xlab="position in single-color",ylab="position in 2-color",col=rgb(0,0,0,0.2),cex=0.5)
  cor1 <- (cor(df[,1],df[,2]))^2
  abline(0,1,col=rgb(1,0,0,0.4))
  legend("bottomright",bty="n",legen=paste0("correlation = ",format(cor1,digits=2)),cex=0.8,text.col=rgb(0,0,0,0.5))

  df <- make_same_length(cmap2[,"Position"],merged_ch2[,"Position"])
  plot(df[,1],df[,2],pch=16,main="single-color vs 2-color (channel 2)",xlab="position in single-color",ylab="position in 2-color",col=rgb(0,0,0,0.2),cex=0.5)
  cor2 <- (cor(df[,1],df[,2]))^2
  abline(0,1,col=rgb(1,0,0,0.4))
  legend("bottomright",bty="n",legen=paste0("correlation = ",format(cor2,digits=2)),cex=0.8,text.col=rgb(0,0,0,0.5))

  outlist <- list("cor1"=cor1,
                  "cor2"=cor2)
#   print(outlist)

  return(outlist)
}

#wrapper for outputting a two-color cmap
#input: the path coordinates and the assembly cmap files
#output: a df formatted like a cmap
make_cmap <- function(path_coordinates,references,plot_flag=F,test_flag=F) {
#   #get assembly data
#   references <- get_assemblies(asm1,asm2)

  #shift coordinates so that the scaffold would start at 0
  coord_range <- range(c(path_coordinates[,"rel_start"],path_coordinates[,"rel_end"]))
  path_coordinates[,"rel_start"] <- path_coordinates[,"rel_start"]-coord_range[1]
  path_coordinates[,"rel_end"] <- path_coordinates[,"rel_end"]-coord_range[1]

  path_coordinates_split <- split(path_coordinates,path_coordinates[,"node_type"])

  #get cmap df
  cmap_list_asm1 <- output_cmap_df(path_coordinates_split[[1]],references[["asm1_cmap"]])
  cmap_list_asm2 <- output_cmap_df(path_coordinates_split[[2]],references[["asm2_cmap"]])

  cmap_path <- rbind(cmap_list_asm1[["cmap_path"]],cmap_list_asm2[["cmap_path"]]) #combine cmap paths

  cmap_list <- list("asm1"=cmap_list_asm1[["cmap_list"]],"asm2"=cmap_list_asm2[["cmap_list"]]) #combine cmap lists
  #save(list="cmap_list",file="/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/Human_willsample_twocolor_cmap_list.RData") #for testing

  #make intermediate cmaps
  if (test_flag) {
    single_col_cmaps <- test_intermediate_cmap(cmap_list)
  }

  if (plot_flag) {
    plot_schematic(cmap_path,order_flag=T,plot_flag=T,title_str='post-trimming')
  }

  #plot cmap starts and ends
  if (plot_flag) {
    title_str <- "cmap post-trimming"
    plot_cmaps(cmap_list,plot_flag=T,title_str)
    bool_cmap_check <- check_cmap_list(cmap_list)
    if (!bool_cmap_check) {
      stop("cmap entries overlap each other")
    }
  }
  cmap_merged <- merge_cmap_df(cmap_list_asm1[["cmap_list"]],cmap_list_asm2[["cmap_list"]])

  if (plot_flag) {
    plot_cmap(cmap_merged,plot_flag=T)
  }

  if (test_flag) {
    cmap_comparison <- compare_cmaps(single_col_cmaps,cmap_merged)
  }

  return(cmap_merged)
}

renumber_mapid <- function(cmap_list,mode="sequential",parallel_flag=F,cores=2) {

  if (mode=="sequential") { #simple renumber the maps based on sequential order
    if (parallel_flag) {
      renumbered <- mclapply(1:length(cmap_list),function(i) {
        cmap <- cmap_list[[i]]
        cmap[,"CMapId"] <- i

        return(cmap)
      },mc.cores=cores)
    } else {
      renumbered <- lapply(1:length(cmap_list),function(i) {
        cmap <- cmap_list[[i]]
        cmap[,"CMapId"] <- i

        return(cmap)
      })
    }
  } else {
    stop("mode not recognized.")
  }

  return(renumbered)
}

#round the numeric columns in a df
round_df <- function(df,digits=0) {

  outlist <- lapply(df,function(y) {
    if (is.numeric(y)) {
      return(round(y,digits))
    } else {
      return(y)
    }
  })

  outdf <- data.frame(outlist)

  return(outdf)
}

#get maps in an assembly
#input is the list of ids and the cmap
#output is a cmap dataframe with the selected maps
get_maps <- function(ids,cmap,invert_flag=F) {

  if (invert_flag) {
    cmap_selected <- cmap[!cmap[,"CMapId"]%in%ids,]
  } else {
    cmap_selected <- cmap[cmap[,"CMapId"]%in%ids,]
  }

  return(cmap_selected)
}

#get maps that did not participate in sandwich
get_leftover_cmap <- function(paths,asm1,asm2) {

  #testing to make sure the nodes are of the right types
  mode <- 1
  check_source <- check_assembly_source1(paths,asm1,asm2,mode)
  stopifnot(check_source)

  paths_df <- lapply(1:length(paths),function(i) {
    path <- paths[[i]][["path_coordinates"]]
    path <- cbind("id"=i,path)
  })
  path_df <- do.call(rbind,paths_df)

  path_df_split <- split(path_df,path_df[,"node_type"])

  nodes_asm1 <- path_df_split[["asm1"]][,"node"]
  nodes_asm2 <- path_df_split[["asm2"]][,"node"]

  asm1 <- readcmap(asm1)
  asm2 <- readcmap(asm2)

  asm_normalized <- normalize_df(asm1,asm2)

  asm1 <- asm_normalized[["df1_trimmed"]]
  asm2 <- asm_normalized[["df2_trimmed"]]

  leftovers_asm1 <- get_maps(nodes_asm1,asm1,invert_flag=T)
  leftovers_asm2 <- get_maps(nodes_asm2,asm2,invert_flag=T)

  outlist <- list("leftovers_asm1"=leftovers_asm1,
                  "leftovers_asm2"=leftovers_asm2)

  return(outlist)
}

check_leftovers <- function(leftovers,asm1,asm2,paths) {
  leftovers_asm1 <- leftovers[["leftovers_asm1"]]
  leftovers_asm2 <- leftovers[["leftovers_asm2"]]

  leftovers_asm1_ids <- unique(leftovers_asm1[,"CMapId"])
  leftovers_asm2_ids <- unique(leftovers_asm2[,"CMapId"])

  asm1 <- readcmap(asm1)
  asm2 <- readcmap(asm2)

  asm1_ids <- unique(asm1[,"CMapId"])
  asm2_ids <- unique(asm2[,"CMapId"])

  #are the ids from leftover maps found in the initial assemblies
  cmap_not_found1 <- sum(!leftovers_asm1_ids%in%asm1_ids)
  cmap_not_found2 <- sum(!leftovers_asm2_ids%in%asm2_ids)

  paths_df <- lapply(1:length(paths),function(i) {
    path <- paths[[i]][["path_coordinates"]]
    path <- cbind("id"=i,path)
  })
  path_df <- do.call(rbind,paths_df)

  path_df_split <- split(path_df,path_df[,"node_type"])

  nodes_asm1 <- path_df_split[["asm1"]][,"node"]
  nodes_asm2 <- path_df_split[["asm2"]][,"node"]

  #are the ids from the leftover maps found in the paths
  nodes_found1 <- sum(leftovers_asm1_ids%in%nodes_asm1)
  nodes_found2 <- sum(leftovers_asm2_ids%in%nodes_asm2)

  outlist <- list("cmap_not_found1"=cmap_not_found1,
                  "cmap_not_found2"=cmap_not_found2,
                  "nodes_found1"=nodes_found1,
                  "nodes_found2"=nodes_found2)

  check_flag <- cmap_not_found1+cmap_not_found2+nodes_found1+nodes_found2

  if (check_flag>0) {
    print(outlist)
  }

  return(outlist)
}

#########################added by Jian 4/3/2016 ################
splitTwocolorCmap <- function(cmap_list){
    splitHelper <- function(cmap, channel=1){
        cmap.dt <- as.data.table(cmap)
        cmap.dt.1 <- cmap.dt[LabelChannel == channel | LabelChannel == 0]
        cmap.dt.1 <- cmap.dt.1[LabelChannel > 0, LabelChannel := 1]
        cmap.dt.1 <- fixcmap(cmap.dt.1)

        return(as.data.frame(cmap.dt.1))
    }
    fixcmap <- function(cmap){
        cmap$NumSites <-sum(cmap$LabelChannel > 0)
        cmap$SiteID <- seq(1:nrow(cmap))
        return(cmap)
    }
    split.cmap.list1 <- lapply(cmap_list, function(cmap){splitHelper(cmap, 1)})
    split.cmap.list2 <- lapply(cmap_list, function(cmap){splitHelper(cmap, 2)})
    return(list(cmap_col1=split.cmap.list1, cmap_col2=split.cmap.list2))
}







# #20151103 EL - outputting leftover maps for analysis
# bspqi_asm <- "/home/users/elam/data/20151027_NA12878_BspQI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap" #V3
# bsssi_asm <- "/home/users/elam/data/20151027_NA12878_BssSI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap" #V3
#
# mode <- "single"
#
# load("/home/users3/elam/20151030_2-col_V3_sandwich/20151030_v3_sub30X_outlist.RData")
# paths <- outlist[["paths"]]
#
# output <- get_leftover_cmap(paths,bspqi_asm,bsssi_asm)
# # output1 <- check_leftovers(output,bspqi_asm,bsssi_asm,paths)
# # print(output1)
#
# input <- list("filename"="/home/users3/elam/20151030_2-col_V3_sandwich/20151030_v3_sub30X_leftovers_BspQI.cmap",
#               "cmap"=output[["leftovers_asm1"]])
#
# out <- write_cmap(input,mode)
#
# load("/home/users3/elam/20151030_2-col_V3_sandwich/20151030_v3_sub130X_outlist.RData")
# paths <- outlist[["paths"]]
#
# output <- get_leftover_cmap(paths,bspqi_asm,bsssi_asm)
# # output1 <- check_leftovers(output,bspqi_asm,bsssi_asm,paths)
# # print(output1)
#
# input <- list("filename"="/home/users3/elam/20151030_2-col_V3_sandwich/20151030_v3_sub130X_leftovers_BspQI.cmap",
#               "cmap"=output[["leftovers_asm1"]])
#
# out1 <- write_cmap(input,mode)

# #20151103 EL - trying to output the leftover maps as well
# load("/home/users3/elam/20151101_2-col_V2_v_V3/20151103_v3_mol_v2_qi_v3_si_sub40X_outlist.RData")
# paths <- outlist[["paths"]]
# bspqi_asm <- "/home/users/elam/data/20151031_NA12878_50X_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap" #V2
# bsssi_asm <- "/home/users/elam/data/20151027_NA12878_BssSI_ML_V3_50x_mapped_AT/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap" #V3
#
# output <- get_leftover_cmap(paths,bspqi_asm,bsssi_asm)
# output1 <- check_leftovers(output,bspqi_asm,bsssi_asm,paths)

# #20151103 EL - testing write_cmap such that the columns are of the right types
# load("/home/users/elam/20151030_2-col_V3_sandwich/20151030_v2_sub110X_cmap_list.RData")
#
# cmap <- cmap_list[[1]]
# cmap <- round_df(cmap,2)
# filename <- paste0("/home/users/elam/20151030_2-col_V3_sandwich/20151103_test.cmap")
# write.table(cmap,file=filename,quote=F,sep="\t",row.names=F,col.names=F)

# #main
# load("/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/Human_willsample_twocolor_combined_paths.RData") #path data
#
# align_dir <- "/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/"
# prefix <- "Human_willsample_twocolor_combined"
#
# asm1 <- paste0(align_dir,prefix,"_1_bppAdj.cmap")
# asm2 <- paste0(align_dir,prefix,"_2_bppAdj.cmap")
#
# #testing path 1
# i <- 4
# path <- paths[[i]]
# path_coordinates <- path[["path_coordinates"]]
# par(mfrow=c(3,1))
#
# title_str <- paste0("path ",i,' pre-trimming')
# plot_flag <- T
# order_flag <-T
#
# plot_schematic(path_coordinates,plot_flag,order_flag,title_str)
#
# plot_flag <- F
# test_flag <- F
#
# output <- make_cmap(path_coordinates,asm1,asm2,plot_flag,test_flag)
#
# stop("done")
#
# #big test
# plot_flag <- T
# order_flag <-T
#
# pdf(file="/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/20151014_cmap_out.pdf",height=8,width=12)
# par(mfrow=c(3,1))
#
# try({
#   cmap_list <- lapply(1:length(paths),function(i) {
#     print(i)
#
#     path <- paths[[i]]
#     path_coordinates <- path[["path_coordinates"]]
#
#     title_str <- paste0("path ",i,' pre-trimming')
#
#     plot_schematic(path_coordinates,plot_flag,order_flag,title_str)
#
#     return(make_cmap(path_coordinates,asm1,asm2,plot_flag=F))
#   })
# })
#
# dev.off()
#
# save(list="cmap_list",file="/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/20151014_cmap_list_testing.RData")
#
# stop("cmap df done")
#
# lapply(1:length(cmap_list),function(i) {
#   cmap <- cmap_list[[i]]
#   filename <- paste0("/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/","20151014_contig",i,".cmap")
#   write.table(cmap,file=filename,quote=F,sep="\t",row.names=F,col.names=F)
# })
#
# write.table(paste0("/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/","20151014_contig",1:length(cmap_list),".cmap"),
#             file="/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/20151014_id_list.txt",
#             quote=F,sep="\t",row.names=F,col.names=F)
#
# stop("cmap generation done")
