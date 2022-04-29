#20150626 EL - read cmap with headers

# readcmap <- function(filename) {
#   cmap <- NULL
#   try( {
#
#     cmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50) #adjust nrows as needed
#
#     if (ncol(cmap) == 9) {
#       colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
#                             "SiteID", "LabelChannel", "Position",
#                             "StdDev", "Coverage", "Occurrence"
#       )
#     } else if (ncol(cmap) == 11) {
#       colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
#                             "SiteID", "LabelChannel", "Position",
#                             "StdDev", "Coverage", "Occurrence",
#                             "GmeanSNR", "lnSNRsd" #new columns
#       )
#     } else if (ncol(cmap) == 12) {
#       colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
#                             "SiteID", "LabelChannel", "Position",
#                             "StdDev", "Coverage", "Occurrence",
#                             "GmeanSNR", "lnSNRsd", #new columns
#                             "SNR" #new columns
#       )
#     }
#
#     class(cmap$CMapId) <- "character"
#
#     cmap <- as.data.frame(scan(filename, what=cmap, comment.char="#", quiet=TRUE),stringsAsFactors=F)
#
#   } ) # try
#   return(cmap)
# }

readcmap_w_headers <- function(filename,sep="") { #output a list with the header lines and the actual content
    header0 <- readLines(filename,n=100) #find header lines
    header <- header0[regexpr("^#",header0)>=0]
    cmap <- readcmap(filename)
    return(list("header"=header,
                "cmap"=cmap))
}


readxmap_w_headers <- function(filename,sep="") { #output a list with the header lines and the actual content
  header0 <- readLines(filename,n=100) #find header lines
  header <- header0[regexpr("^#",header0)>=0]
  xmap <- readxmap(filename)
  return(list("header"=header,
              "xmap"=xmap))
}

#testing
# cmap <- "/mnt/bionf_tmp/elam/20150607_two-color/MR1_sim2/20150521_MR1_bspq1.cmap"
# cmap1 <- readcmap(cmap)
# cmap2 <- readcmap_w_headers(cmap)
