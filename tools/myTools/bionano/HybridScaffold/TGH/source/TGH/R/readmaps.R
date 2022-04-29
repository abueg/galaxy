#utility functions to read standard files
# readCMap <- function(filename) {
#     cmap <- NULL
#     try( {
#
#         #message(filename)
#         cmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50, fill = TRUE) #adjust nrows as needed
#         cmap <- as.data.frame(scan(filename, what=cmap, comment.char="#", fill = TRUE, quiet=TRUE))
#
#         if (ncol(cmap) == 9) {
#             colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
#                                   "SiteID", "LabelChannel", "Position",
#                                   "StdDev", "Coverage", "Occurrence"
#             )
#         } else if (ncol(cmap) == 11) {
#             colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
#                                   "SiteID", "LabelChannel", "Position",
#                                   "StdDev", "Coverage", "Occurrence",
#                                   "GmeanSNR", "lnSNRsd" #new columns
#             )
#         } else if (ncol(cmap) == 12) {
#         colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
#                               "SiteID", "LabelChannel", "Position",
#                               "StdDev", "Coverage", "Occurrence",
#                               "GmeanSNR", "lnSNRsd", #new columns
#                               "SNR" #new columns
#             )
#         }
#     } ) # try
#     return(cmap)
# } #readCMap
#
# readXMap <- function(filename) {
#     xmap <- NULL
#     try( {
#
#         #message(filename)
#         xmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50)
#         xmap <- as.data.frame(scan(filename, what=xmap, comment.char="#", quiet=TRUE))
#
#         colnames(xmap) <- c( "XmapEntryID", "QryContigID", "RefcontigID",
#                              "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
#                              "Orientation", "Confidence", "HitEnum"
#                             )
#     },
#     silent = TRUE
#     ) #try
#     return(xmap)
# } #readXMap
#
readSMap <- function( filename ) {
    smap <- NULL;
    try ( {

        #message(filename)
        smap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=1000000)

        if (ncol(smap) == 12) {
            colnames( smap ) <- c( "SmapEntryID", "QryContigID", "RefcontigID1", "RefcontigID2",
                                   "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                                   "Confidence", #check for column for orientation
                                   "Type", "XmapID1", "XmapID2"
            )
        } else if (ncol(smap) == 13) {
            #01222015 EL
            colnames( smap ) <- c( "SmapEntryID", "QryContigID", "RefcontigID1", "RefcontigID2",
                                   "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                                   "Confidence",
                                   "Type", "XmapID1", "XmapID2","LinkID"
            )

#             colnames( smap ) <- c( "SmapEntryID", "QryContigID", "RefcontigID1", "RefcontigID2",
#                                    "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
#                                    "Orientation", "Confidence", #check for column for orientation
#                                    "Type", "XmapID1", "XmapID2"
#             )
        } else if (ncol(smap)==17) {
            #03242015 EL
            colnames(smap) <- c("SmapEntryID","QryContigID","RefcontigID1","RefcontigID2",
                                "QryStartPos","QryEndPos","RefStartPos","RefEndPos",
                                "Confidence","Type","XmapID1","XmapID2",
                                "LinkID","QryStartIdx","QryEndIdx","RefStartIdx","RefEndIdx")
        }
    },
    silent = TRUE #added 03202014 to suppres read smap error there are no lines
    ) #try
    return(smap)
} #readSMap

#improved functions
#read a .map file
readmap <- function(filename) {
    map <- NULL
    try( {

        con1 <- pipe(paste0("cut -f1-18 ",filename))
        #map <- read.table(con1, skip=8, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, col.names=1:18, nrows=50) #read in partial file
        map <- read.table(con1, skip=6, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, col.names=1:18, nrows=50) #20150324 EL - format change

        colnames(map) <- c("MappedMoleculeId", "MoleculeId", "MoleculeIndex", "ContigId",
                           "Score", "Zscore", "Direction", "StartLocation",
                           "EndLocation", "StartMatchLocation", "EndMatchLocation", "DetectedLabelCount",
                           "TruePositiveLableCount", "FalsePositiveLableCount", "FalseNegativeLabelCount", "Log10Pvalue",
                           "LeftEndChim","RightEndChim")
        class(map$MoleculeId) <- "character" #change MoleculeId to "character" to prevent overflow

        con2 <- pipe(paste0("cut -f1-18 ",filename))
        #map <- as.data.frame(scan(con2, what=map, skip=8, quiet=TRUE),stringsAsFactors=F) #read entire file with scan
        map <- as.data.frame(scan(con2, what=map, skip=6, quiet=TRUE),stringsAsFactors=F) #20150324 EL - format change
        close(con2)

    } ) # try
    return(map)
} #readmap

#read a .map file (11132014 - RStudio Server version)
readmap1 <- function(filename) {
  map <- NULL
  try( {

    con1 <- pipe(paste0("head -n 50 ",filename," | cut -f1-18")) #cut has to finish what it's doing first
    #map <- read.table(con1, skip=8, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, col.names=1:18)
    map <- read.table(con1, skip=6, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, col.names=1:18, nrows=50) #20150324 EL - format change

    colnames(map) <- c("MappedMoleculeId", "MoleculeId", "MoleculeIndex", "ContigId",
                       "Score", "Zscore", "Direction", "StartLocation",
                       "EndLocation", "StartMatchLocation", "EndMatchLocation", "DetectedLabelCount",
                       "TruePositiveLableCount", "FalsePositiveLableCount", "FalseNegativeLabelCount", "Log10Pvalue",
                       "LeftEndChim","RightEndChim")
    class(map$MoleculeId) <- "character"

    con2 <- pipe(paste0("cut -f1-18 ",filename))
    #map <- as.data.frame(scan(con2, what=map, skip=8, quiet=TRUE),stringsAsFactors=F)
    map <- as.data.frame(scan(con2, what=map, skip=6, quiet=TRUE),stringsAsFactors=F) #20150324 EL - format change
    close(con2)

  } ) # try
  return(map)
} #readmap

#read an .xmap file

#112414 - modified to convert type for MoleculeId; assignment of colnames moved up
#deprecated, modfity to version using a general method of reading bng maps
#readxmap <- function(filename) {
#      xmap <- NULL
#      xmap <-  tryCatch({
#                    xmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50, sep="\t")
#               }, error = function(err){ #handling case when xmap has no alignment
#                 if(grepl("no lines available", err$message)){
#                   warning(paste0("xmap file ", filename, " has zero row"))
#                   xmap <- as.data.frame(matrix(0, nrow=0, ncol=10))
#                   return(xmap)
#                 }
#                 stop(paste0("Cannot read xmap file: ", filename))
#              }
#      )
#
#      if (ncol(xmap) == 10) {
#          colnames(xmap) <- c("XmapEntryID", "QryContigID", "RefcontigID",
#                              "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
#                              "Orientation", "Confidence", "HitEnum")
#        } else if (ncol(xmap) == 14) { #121514 - new columns
#          colnames(xmap) <- c("XmapEntryID", "QryContigID", "RefcontigID",
#                              "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
#                              "Orientation", "Confidence", "HitEnum",
#                              "QryLen","RefLen","LabelChannel","Alignment") #121514
#        }
#
#        class(xmap$QryContigID) <- "character" #change QryContigID to "character" to prevent overflow
#        class(xmap$RefcontigID) <- "character" #change QryContigID to "character" to prevent overflow
#        class(xmap$QryStartPos) <- "double"
#        class(xmap$QryEndPos) <- "double"
#        class(xmap$RefStartPos) <- "double"
#        class(xmap$RefEndPos) <- "double"
#        class(xmap$QryLen) <- "double"
#        class(xmap$RefLen) <- "double"
#        xmap <- as.data.frame(scan(filename, what=xmap, comment.char="#", quiet=TRUE),stringsAsFactors=FALSE)
    #} ) #try
#    return(xmap)
#} #readxmap

#01222017 using a central method to read bng map data format and then do file specific processing
readxmap <- function(filename) {
  xmap <-read.bng.maps(filename)
  ind <- which(names(xmap) == 'RefContigID')
  if(length(ind) > 0){
    names(xmap)[[ind]] <- "RefcontigID"
  }
  class(xmap$QryContigID) <- "character" #change QryContigID to "character" to prevent overflow
  class(xmap$RefcontigID) <- "character" #change QryContigID to "character" to prevent overflow
  class(xmap$QryStartPos) <- "double"
  class(xmap$QryEndPos) <- "double"
  class(xmap$RefStartPos) <- "double"
  class(xmap$RefEndPos) <- "double"
  class(xmap$QryLen) <- "double"
  class(xmap$RefLen) <- "double"
  class(xmap$Confidence) <- "double"
  return(xmap)
} #readxmap


#read all xmap from a directory and merge them into a single data structure
read.xmap.batch <- function(dir.path, pattern="*.xmap"){
  xmap.files <- dir(dir.path, pattern=pattern, full.names = T)
  xmap.list <- lapply(xmap.files, function(f){
                                      tryCatch(
                                        readxmap(f),
                                        error=function(e){return(NULL)})})
  xmap.combined <- Reduce(function(x,y){rbind(x,y)}, xmap.list)
  print(paste0("Total number of alignments read: ", nrow(xmap.combined)))
  return(xmap.combined)
}


#read a .cmap file
readcmap <- function(filename) {
    cmap <- NULL
    try( {
        cmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50) #adjust nrows as needed
        header <- read.header(filename)
        if(length(header) < 1){
          warning("Cannot find proper header line, please check file")
        }
        if(length(header) == ncol(cmap)){
            colnames(cmap) <- header
        }else{
          if (ncol(cmap) == 9) {
              colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
                                    "SiteID", "LabelChannel", "Position",
                                    "StdDev", "Coverage", "Occurrence"
              )
          } else if (ncol(cmap) == 11) {
              colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
                                    "SiteID", "LabelChannel", "Position",
                                    "StdDev", "Coverage", "Occurrence",
                                    "GmeanSNR", "lnSNRsd" #new columns
              )
          } else if (ncol(cmap) == 12) {
              colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
                                    "SiteID", "LabelChannel", "Position",
                                    "StdDev", "Coverage", "Occurrence",
                                    "GmeanSNR", "lnSNRsd", #new columns
                                    "SNR" #new columns
              )
          } else{
            warning(paste0("For file ", filename, " Column header do not seem to be consistent with data, trimming header instead"))
            header <- header[1:ncol(cmap)]
            colnames(cmap) <- header
          }
        }
        #R sometime has some inconsistency in determining the type of some columns, this affect downstream processing
        #we standardized them here
        class(cmap$CMapId) <- "character"
        if(!is.null(cmap$Mask)){
          class(cmap$Mask) <- "character"
        }
        cmap$CMapId <- as.character(cmap$CMapId)
        cmap$ContigLength <- as.numeric(cmap$ContigLength)
        cmap$Position <- as.numeric(cmap$Position)
        #browser()
        cmap <- as.data.frame(scan(filename, what=cmap, comment.char="#", quiet=TRUE),stringsAsFactors=F)

    } ) # try
    return(cmap)
} #readcmap

read.header <- function(filename){
  header.lines <- readLines(filename, n=50)
  header.lines <- header.lines[which(grepl("^#h ", header.lines))]
  if(length(header.lines)==0){
    return(data.frame())
  }
  header <- strsplit(header.lines, " ")
  header <- strsplit(header[[1]][[2]], "\t")
  return(header[[1]])
}


#read a generic bng map file that is has the tab-delimited format
read.bng.maps <- function(filename){
  bmap <- data.frame()
  try( {
     header <- read.header(filename)
     bmap <-  tryCatch({
                  bmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50, sep="\t")
                  for( i in 1:ncol(bmap)){
                    if(typeof(bmap[[i]]) == 'integer'){
                      class(bmap[[i]]) <- "double"
                    }
                  }
                  bmap
              }, error = function(err){ #handling case when xmap has no alignment
                    if(grepl("no lines available", err$message)){
                          return(NULL)
                    }
                    stop(paste0("Cannot read xmap file: ", filename))
                }
             )
     if(is.null(bmap)){
       bmap <- as.data.frame(matrix(0, nrow=0, ncol=length(header)))
     }
     if(length(header) == ncol(bmap)){
      colnames(bmap) <- header
    }else{
      warning("header line is not present in file or does not seem to have the same columns as data, not setting header")
    }
    bmap <- as.data.frame(scan(filename, what=bmap, comment.char="#", quiet=TRUE),stringsAsFactors=F)

  } ) # try
  return(bmap)

}



#read a .bed file
readbed <- function(filename) {
    bed <- NULL
    try( {

        bed <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)

        if (ncol(bed)==5) {
            colnames(bed) <- c( "chr", "start", "end", "walkID", "orientation")
        } else if (ncol(bed)==3) {
            colnames(bed) <- c( "chr", "start", "end")
        }

    } ) #try
    return(bed)
} #readbed

#20150522 EL - read an .err file
readerr <- function(filename) {
  err <- NULL
  try({

    err <- read.table(filename,comment.char="#",header=TRUE,stringsAsFactors=FALSE,fill=FALSE)

    if (ncol(err)==19) {
      colnames(err) <- c("iteration","fp","fn","sf","sd","bpp","res","maps","log10LR","aligned",
                         "log10LR.aligned","bppSD","fp.rate","sr","se","label.density","resSD","mres","mresSD")
    }

    use_columns <- c("log10LR","aligned","log10LR.aligned","label.density")

    F <- names(err)%in%use_columns

    err1 <- tail(err,1) #last line
    err2 <- head(tail(err,2),1) #second to last line

    err1[F] <- err2[F]

    err <- err1
    #add a input file name as the first columan
    ExpName <- basename(filename)
    ColNames <- colnames(err)
    err$MoleculesFile <- ExpName
    #re-ordering the columns
    err <- err[c("MoleculesFile", ColNames)]
  }) #try
  return(err)
} #readerr
#err_file <- "/mnt/bionf_tmp/elam/20150519_ESHG/rubin_case.bnx_errEst.err"
#err <- readerr(err_file)

read.err.dir <- function(path){
  err.files <- dir(path, pattern = "*.err$", full.names = T)
  errs.list <- lapply(err.files, readerr)
  return(errs.list)
}



readbnx.header <- function(filename){
  file.pipe <- pipe(paste0("cat ", filename, " | awk '{if($1 ~ /^0/) print}'"))
  mols.header <- read.table(file.pipe, comment.char = c('#'))
  colnames(mols.header) <- c("Channel", "MolID", "Length", "AvgIntensity", "SNR",
                             "NumLabels", "OriginalMolID", "ScanNum", "ScanDirecion",
                             "ChipId", "Flowcell", "RunId", "GlobalScanNum")
  return(mols.header)
}


read.hybrid.conflictfile <- function(filename){
  conflict.header <- read.table(filename, header=T, stringsAsFactors = F, comment.char = '', sep='\t')
  conflict <- tryCatch({read.table(filename, header=F, stringsAsFactors = F, skip = 2)}, error=function(err){NULL})
  #browser()
  if(is.null(conflict)){
    conflict <- as.data.frame(matrix(0, nrow = 0, ncol = ncol(conflict.header)))
  }
  colnames(conflict) <- colnames(conflict.header)
  conflict$refId <- as.character(conflict$refId)
  conflict$qryId <- as.character(conflict$qryId)
  class(conflict$leftRefBkpt) <- "double"
  class(conflict$rightRefBkpt) <- "double"
  class(conflict$leftQryBkpt) <- "double"
  class(conflict$rightQryBkpt) <- "double"
  return(conflict)
}



