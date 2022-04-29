overlap <- function( block1, block2 ) {
    block1 <- range( block1 );
    block2 <- range( block2 );
    return( ( block1[[ 1 ]] < block2[[ 2 ]] ) &
            ( block1[[ 2 ]] > block2[[ 1 ]] ) );
} # overlap

overlaps <- function( block1DataFrame, block2 ) {
    nBlocks <- ncol( block1DataFrame );
    result <- rep( FALSE, nBlocks );
    for ( block1 in 1:nBlocks ) {
        result[[ block1 ]] <- overlap( block1DataFrame[ , block1 ], block2 );
    } # for block1
    return( result );
} # overlaps

#12142014 - check for overlap and calculate amount of overlap
calculate_overlap <- function(block1,block2) {
  block1 <- range(block1)
  block2 <- range(block2)

  flag_overlap <- ((block1[1] < block2[2]) & #check for overlap
                     (block1[2] > block2[1]))

  if (flag_overlap) { #calculate amount of overlap
    diff_merge <- diff(range(block1,block2))
    diff_separate <- diff(block1)+diff(block2)

    overlap <- diff_separate-diff_merge+1
  } else {
    return(-1)
  }

  return(overlap)
}

#12142014 - check for overlap and calculate amount of overlap (based on sort)
calculate_overlap1 <- function(block1,block2) {
  block1 <- range(block1)
  block2 <- range(block2)

  flag_overlap <- ((block1[1] < block2[2]) & #check for overlap
                     (block1[2] > block2[1]))

  sorted <- sort(c(block1,block2))
  if (flag_overlap) { #calculate amount of overlap
    overlap <- sorted[3]-sorted[2]+1
  } else {
    overlap <- -1*(sorted[3]-sorted[2]) #output how far the two intervals are away if no overlap
  }

  return(overlap)
}

#03242015 - comparing an interval vs a point
#for example, a molecule must span a junction on both sides by N bp
# rm(list=ls())
# if (!exists("paste0")) {
#   paste0 <- function(...) { paste(...,sep="") } #load paste0 if undefined
# }
# scripts_dir <- "/home/users/elam/rscripts/"
# source(paste0(scripts_dir,"readmaps.R"))
# xmap <- readxmap("/mnt/bionf_tmp/elam/020215_Broad/hg19_2ver_plus_all_merged_4_BspQ1_v_mol_wMarginT4/hg19_2ver_plus_all_merged_4_BspQ1_v_mol_wMarginT4_contig17.xmap")
# xmap <- head(xmap,5)

overlaps_junction <- function(xmap,junction,pad) { #output counts
  overlap_list <- lapply(1:nrow(xmap),function(i) {
    xmap_row <- xmap[i,]
    start <- xmap_row$RefStartPos
    end <- xmap_row$RefEndPos

    if ((start<junction-pad) & (junction+pad<end)) {
      return(1)
#       print(xmap_row[,1:10])
#       plot(x=c(start,end),y=c(1,1),pch=16)
#       points(x=junction,y=1,col="red")
#       points(c(junction-5000,junction+5000),y=c(1,1),col="blue")
    } else {
      return(0)
    }
  })
  return(sum(unlist(overlap_list)))
}
# test <- overlaps(xmap,150000,5000)

overlaps_junction1 <- function(xmap,junction,pad) { #output xmap entries
  overlap_list <- lapply(1:nrow(xmap),function(i) {
    xmap_row <- xmap[i,]
    start <- xmap_row$RefStartPos
    end <- xmap_row$RefEndPos

    if ((start<junction-pad) & (junction+pad<end)) {
      return(xmap_row)
      #       print(xmap_row[,1:10])
      #       plot(x=c(start,end),y=c(1,1),pch=16)
      #       points(x=junction,y=1,col="red")
      #       points(c(junction-5000,junction+5000),y=c(1,1),col="blue")
    } else {
      return(xmap_row[0,])
    }
  })
  out_xmap <- do.call(rbind,overlap_list)

  return(out_xmap)
}
# test <- overlaps(xmap,150000,94700)
# print(test[,1:9])

#20151201 EL - optimized for comparing a region vs regions in a dataframe
overlaps <- function(region,df,count_flag=T,bool_flag=F) {
  stopifnot(count_flag+bool_flag<=1)

  region_start <- region[1]
  region_end <- region[2]

  #overlap_flag <- (region_start <= df[,2]) & (region_end >= df[,1])
  overlap_flag <- (region_start < df[,2]) & (region_end > df[,1]) #20160930 EL

  if (count_flag) {
    return(sum(overlap_flag))
  }

  if (bool_flag) {
    return(overlap_flag)
  }

  return(df[overlap_flag,])
}
