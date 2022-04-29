#02122015 EL - optimized for speed

#utility functions
extractNumbers <- function(instr,delimiter="") { #get just numbers from string
  gsub("[^[:digit:]]",delimiter,instr)
} #extractNumbers

extractLetters <- function(instr,delimiter="") { #get just letters from string
  gsub("[[:digit:]]",delimiter,instr)
} #extractLetters

#http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r
#returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

#returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

#returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#parse cigar string
#good for longer cigar strings
parseCigar2 <- function(cigar="1M2D3I13D") {
  numbers <- extractNumbers(cigar,"_") #get just numbers
  numberSplit <- unlist(strsplit(numbers,"_")) #split string

  letters <- extractLetters(cigar,"_") #get just letters
  letterSplit <- unlist(strsplit(letters,"_")) #split string
  letterSplit <- letterSplit[letterSplit!=""] #filter out empty entries

  fullseq <- rep(letterSplit,numberSplit)
  #fullseq <- paste(fullseq,collapse = '') #combine into single string
  #print(fullseq)
} #parseCigar1

enumerateLabels_single <- function(qcmap,rcmap,xmap) {

  result <- rep(NA,nrow(qcmap))

  start_ref <- xmap$RefStartPos
  end_ref <- xmap$RefEndPos
  start_qry <- xmap$QryStartPos
  end_qry <- xmap$QryEndPos

  cigar <- (parseCigar2(xmap$HitEnum))

  if (xmap$Orientation == "+") {
    sites_qry <- which(start_qry <= qcmap$Position & qcmap$Position <= end_qry)
    sites_ref <- which(start_ref <= rcmap$Position & rcmap$Position <= end_ref)
  } else {
    sites_qry <- which(end_qry <= qcmap$Position & qcmap$Position <= start_qry)
    sites_ref <- rev(which(start_ref <= rcmap$Position & rcmap$Position <= end_ref))
    cigar <- rev(cigar)
  }

  cigar_ref <- cigar[cigar!="I"]
  cigar_qry <- cigar[cigar!="D"]

  temp <- rep(NA,length(sites_qry))

  temp[cigar_qry=="M"] <- sites_ref[cigar_ref=="M"]

  result[sites_qry] <- temp

  return(result)
}

enumerateLabels <- function(qcmap,rcmap,xmap,colname="ref_siteID") {

  result <- rep(NA,nrow(qcmap))

  for (xmap_id in (1:nrow(xmap))) {
    xmap_entry <- xmap[xmap_id,]

    molID <- xmap[xmap_id,"QryContigID"]

    labelIndices <- which(as.character(qcmap[,"CMapId"])==as.character(molID))
    qcmap_selected <- qcmap[labelIndices,]

    tmp <- enumerateLabels_single(qcmap_selected,rcmap,xmap_entry)

    result[labelIndices] <- tmp
  }

  result <- cbind(qcmap,result)
  colnames(result)[ncol(result)] <- colname

  return(result)
}

split_alignment <- function(alignment) {
  alignment_split <- unlist(strsplit(alignment,split="\\)\\(")) #split the brackets
  alignment_split <- lapply(alignment_split,function(alignment) {gsub("\\(","",alignment)}) #remove "("
  alignment_split <- lapply(alignment_split,function(alignment) {gsub("\\)","",alignment)}) #remove ")"
  alignment_split <- lapply(alignment_split,function(alignment) {strsplit(alignment,",")}) #split qry and ref
  alignment_split <- lapply(alignment_split,function(alignment) {do.call(rbind,alignment)}) #put in dataframe

  alignment <- do.call(rbind,alignment_split)
  colnames(alignment) <- c("ref","qry")
  # print(alignment)

  return(alignment)
}
