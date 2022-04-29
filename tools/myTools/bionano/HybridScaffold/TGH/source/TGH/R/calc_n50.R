#contig or scaffold N50 is a weighted median statistic such that 50% of the entire assembly is contained in contigs or scaffolds equal to or larger than this value
#updated 2/6/2017 JW: in previous implementation, the case with no exact contig at 50% of length and need average, the computation is not exactly correct, update code to be
#consistent with single-hybrid
calc_n50 <- function(contigLengths,mode="avg") {
    if(sum(contigLengths) < 0 || any(contigLengths < 0)){
      stop("Contig lengths cannot be negative")
    }else if(sum(contigLengths == 0)){
      return(0)
    }
    sorted <- sort(contigLengths,decreasing=T)
    n50_low <- sorted[cumsum(sorted)==(sum(sorted)/2)][1]
    if(length(n50_low)==0 || is.na(n50_low)){
      n50_low <- sorted[cumsum(sorted) > (sum(sorted)/2)][1]
    }

    if (mode=="avg") {
        n50_high <- tail(sorted[cumsum(sorted) < (sum(sorted)/2)], 1)
        n50 <- mean(c(n50_low,n50_high))

        return(n50)
    } else if (mode=="integer") {
        return(n50_low)
    } else {
        return(NULL)
    }
}

#output a simple vector of lengths
get_lengths <- function(incmap) {
    cmap_single_line <- incmap[!duplicated(incmap[,"CMapId"]),] #get one line per map
    map_lengths <- cmap_single_line[,"ContigLength"]

    return(map_lengths)
}

#outputs a df with cmapids
get_lengths_df <- function(incmap) {
  cmap_single_line <- incmap[!duplicated(incmap[,"CMapId"]),] #get one line per map
  cmap_single_line <- cmap_single_line[,c("CMapId","ContigLength")]

  return(cmap_single_line)
}

get_nummaps <- function(incmap) {
    cmap_single_line <- incmap[!duplicated(incmap[,"CMapId"]),] #get one line per map
    nrows <- nrow(cmap_single_line)

    return(nrows)
}

get_total_size <- function(incmap) {
    cmap_single_line <- incmap[!duplicated(incmap[,"CMapId"]),] #get one line per map
    total_size <- sum(cmap_single_line[,"ContigLength"])

    return(total_size)
}
