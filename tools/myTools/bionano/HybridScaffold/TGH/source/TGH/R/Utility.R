########################################################################################################
#Section General Utility function
#Strip extension from a file path
stripExtension <- function(filename, ext=""){
  if(ext == ""){
    ind <- tail(gregexpr('\\.', filename)[[1]],1)
  }else{
    ind <- regexpr(ext, filename)

  }
  if(ind < 1){
    warning(paste0("file has not extension", ext))
    return(filename)
  }
  substr(filename,1, ind-1)
}

#get the current path from enviorment
getCurrentPath <- function(args){
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
  return(dirname(script.name))
}

getMemoryUsage <- function(){
  print(sort( sapply(ls(),function(x){object.size(get(x))})))
}


#get appropriate histogram breaks for a sets of data, so their distribution can be ploted on
#the same axis
getHistBreaks <- function(..., breaks){

}

#remove outlier from a distribution
removeOutLier <- function(){

}


#creating hash map functionality using enviroment
createHashMap<-function(keys){
    map <- new.env()
    mapply(function(key, value){map[[key]]<-value}, keys, seq(1,length(keys),1))
    getValues <-function(keys){
        getValue <- function(key){
            return (map[[key]])
        }
        return (sapply(keys, getValue))
    }
    return(getValues)
}

#generic binary search function over a sorted list of values
binarySearch <- function(sorted.pos, qryPos, beginRow, endRow, direction){
    if(sorted.pos[beginRow] > qryPos) {
        #print("base 1")
        if(direction == -1){
            return(-1)
        }else{
            return(beginRow)
        }
    }
    if(sorted.pos[endRow] < qryPos){
        #print("base 2")
        if(direction == -1){
            return(endRow)
        }else{
            return(-1)
        }
    }
    if(endRow == beginRow){
        return(1)
    }
    mid <- as.integer(beginRow + (endRow - beginRow)/2)
    #print(mid)
    if(sorted.pos[mid] < qryPos){
        #print("going right\n")
        rightSearch <- binarySearch(sorted.pos, qryPos, mid+1, endRow, direction)
        if(rightSearch > 0){
            return(rightSearch)
        }else{
            return(mid)
        }
    }else{
        #print("going left\n")
        leftSearch <- binarySearch(sorted.pos, qryPos, beginRow, mid-1, direction)
        if(leftSearch > 0){
            return(leftSearch)
        }else{
            return(mid)
        }
    }

}

#Generate a distribution as a vector of probability or density values
generateExpDensity <- function(targeted.avg, min, max, binsize){
    lengths <- seq(min, max, binsize)
    lambda <- targeted.avg - min
    #d  <- exp(-1*(lengths-targeted.avg)/lambda)
    d <- exp(-1*(lengths-min)/lambda)
    return (data.frame(Length=lengths, p=d))
}

#generate length distribution for each chromosome
generateChrLengthDensity <- function(genome.cmap, lambda, min, max, binsize){
    chr.ids <- unique(genome.cmap$CMapId);
    length.distr <- generateExpDensity(lambda, min, max, binsize)
    length.distr.list <- list()
    #for(id in chr.ids){
    for(id in 1:length(chr.ids)){
        length.distr.list[[id]] <- length.distr
    }
    return (length.distr.list)
}

#Get basic stat summary for data
getStats <- function(dat){
  stat <- data.frame(mean=mean(dat),
               median=median(dat),
               std=var(dat)^0.5,
               mad=mad(dat))
}


################################################################################################################
#Section CMAP Utils
#Utility function for accessing and processing CMAP files

getN50 <- function(contig.lens){
  contig.lens <- as.numeric(sort(contig.lens,decreasing = T))
  sums <- cumsum(contig.lens)
  if(sum(sums) == 0){
    return(0)
  }
  fract.sums <- sums /sum(contig.lens)
  ind_low <- head(which(fract.sums > 0.5),1)
  ind_mid <- head(which(fract.sums == 0.5),1)
  ind_hi <- tail(which(fract.sums < 0.5),1)
  if(length(ind_mid) > 0){
    return(contig.lens[ind_mid])
  }else{
    return(mean(c(contig.lens[ind_low],
                contig.lens[ind_hi])))
  }
}

#get the chromosome start and end row number in a reference cmap
getChrRowMap <- function(cmap.ref){
    numSites <- getNumSites(cmap.ref)
    #note the additional one account for the additonal end sites in cmap file
    startCoords <- Reduce(function(x,y){x+y+1}, numSites, accumulate=T, init=1)
    endCoords <- startCoords[-length(startCoords)] + numSites
    map<-mapply(function(begin, end){list(begin, end)}, startCoords[-length(startCoords)], endCoords)
    colnames(map) <- names(endCoords)
    return(map)
}


#get the number of sites for each cmap contig/molecule, return
#a named vector of the number of sites per molecule
getNumSites <- function(cmap.ref){
    ChrIDs <- unique(cmap.ref$CMapId)
    numSites <- sapply(ChrIDs, function(id){cmap.ref$NumSites[cmap.ref$CMapId == id & cmap.ref$SiteID == 1]})
    return (numSites)
}

#create a empty cmap data frame

#ToDo make this create empy cmap based on a template
#because cmap of different version may have different columns
createEmptyCMap <- function(nRows, templateCols=""){
    dummy <- rep( NA, nRows);
    if(templateCols==""){
      newMoleculesCMap <- data.frame(
            "CMapId"=dummy,
            "ContigLength"=dummy,
            "NumSites"=dummy,
            "SiteID"=dummy,
            "LabelChannel"=dummy,
            "Position"=dummy,
            "StdDev"=dummy,
            "Coverage"=dummy,
            "Occurrence"=dummy,
            "GmeanSNR"=dummy,
            "lnSNRsd"=dummy
            );
        return(newMoleculesCMap)
    }else{
       #not-yet implemented
    }

}

#apply a function to a range of data in a window
rolling.function <- function(data, window.width, f){
  L <- nrow(data)
  values <- unlist(
    lapply(seq(1,L,1),
           function(i){
             left <- i - window.width
             right <- i + window.width
             if(left < 0){
               left <- 1
             }
             if(right > L){
               right <- L
             }
             return(f(data[left:right,]))

           }))
  return(values)
}

getLabelDensity <- function(cmap, windowWidth=10){
  label.density <- rolling.function(cmap, windowWidth,
                                    function(d){
                                      Id <- d$CMapId[order(table(d$CMapId), decreasing=TRUE)[1]]
                                      d <- subset(d, CMapId == Id)
                                      density <- nrow(d)/(d$Position[nrow(d)] - d$Position[1])*100e3
                                      if(density < 0){
                                        return(0)
                                      }
                                      return(density)
                                    })
}

#add snr columns to some older version of cmap that do not have these
#information
addSNRCols <- function(cmap){
    cmap["GmeanSNR"] = rep(1, nrow(cmap))
    cmap["lnSNRsd"] = rep(0, nrow(cmap))
    return(cmap)
}

#someitmes in simulation some positions are duplicated after rounding
#we add a small offset to the position coordinate to make them unique
adjustDuplicatePosition <- function(cmap, offset=0.05){
  diff.Id <- cmap$CMapId[-1] - cmap$CMapId[-nrow(cmap)]
  diff.pos <- cmap$Position[-1] - cmap$Position[-nrow(cmap)]
  duplicate.ind <- which(diff.Id == 0 & diff.pos == 0) +1
  print(paste0("Adjusting total labels with same position: ", length(duplicate.ind)))
  cmap$Position[duplicate.ind] <- cmap$Position[duplicate.ind] + offset
  return(cmap)
}

#check consistensy of molecules
molCmapCheck <- function(mol.cmap, detail=FALSE){
  lengthcheck1 <- which(mol.cmap$ContigLength < mol.cmap$Position)
  lengthcheck2 <- which(mol.cmap$ContigLength < 0)
  lengthcheck <- which(lengthcheck1 & lengthcheck2)
  print(paste0("Contig length not-consistent: ", length(lengthcheck)))
  mol.cmap$diff1 <- c(diff(mol.cmap$Position),0)
  mol.cmap$diff2 <- c(diff(as.numeric(mol.cmap$CMapId)), 0)
  pos.order <- which(mol.cmap$diff1 < 0 & mol.cmap$diff2 == 0)
  print(paste0("Label position not in valid order: ", length(pos.order)))
  mol.cmap$diff3 <- c(diff(mol.cmap$SiteID),1)
  siteID.check <- which(mol.cmap$diff3 != 1 & mol.cmap$diff2 == 0)
  #print(which(mol.cmap$diff3 != 1 & mol.cmap$diff2 == 0))
  print(paste0("SiteID  not in valid order: ", length(siteID.check)))
  return(list(LenCheck=lengthcheck, Poscheck=pos.order, SiteIDCheck=siteID.check))
}


#after modifying a molecules cmap files this function perform some bookeeping
#on numsite and siteID and Position to make sure they of proper values
bookeepMolecules <- function(molecules, sort=FALSE){
    if(sort){
        molecules <- molecules[with(molecules, order(molecules[,"CMapId"], molecules[,"Position"])),]
    }
   # browser()
    numSites <- table(molecules$CMapId);
    #we want to set the class back to basic type, otherwise sometime has weird effects downstream
    molecules$NumSites <-
        as.integer(numSites[molecules$CMapId] - 1);
    #cat( "Assigned NumSites\n" );
    indices <- c( 0, which(molecules$LabelChannel == 0 ) ) + 1;
    numSites <- molecules$NumSites[ head( indices, -1 ) ];
    molecules$SiteID <- unlist( lapply( numSites, function( x ) 1:( x + 1 ) ) );
    #cat( "Assigned SiteID\n" );
    #By convention we leave a little room for the beginning of the molecules
    if(molecules$Position[1] < 20){
      molecules$Position <- molecules$Position + 20
      molecules$ContigLength <- molecules$ContigLength + 20
    }
    return(molecules)
}



# creating a hashtable for a pair of keys in the cmap file
#returning the row number in cmap that correspond to that particular
#pair of keys
createHashCMap<-function(Cmap, key1, key2){
  map<-new.env()
  sep='|'
  keys=paste0(Cmap[[key1]], sep, Cmap[[key2]])
  #browser();
  mapply(function(key, value){map[[key]]<-value}, keys, seq(1,length(keys),1))
  getCoords<-function(molIDs, siteIDs=''){
    molIDs<-as.list(molIDs)
    siteIDs<-as.list(siteIDs)
    getCoord<-function(molID, siteID=''){
      if(length(siteID)>1){
        key<-paste0(molID, sep, siteID)
        #print(paste0("key is", key))
        #if(is.na(Mol_list[key])){
      }else{
        key<-molID
      }
      #browser();
      return(map[[key]])
      #}
    }
    #browser();
    return(unlist(mapply(getCoord, molIDs, siteIDs)))
  }
  return(getCoords)
}

#split genome cmap into a list, which contain one chromosome in each element of the list
#note this is only suitable when number of contigs are small, it is slow for large number
#of contigs/molecules
splitGenome <- function(genome.cmap){
    chr.ids <- unique(genome.cmap$CMapId);
    chr.list <- lapply(chr.ids,
                        function(id){genome.cmap[genome.cmap$CMapId == id,]}
                        )
    return (chr.list)
}


#split cmap into a list
splitMolecules <- function(cmap){
   if(!is.data.table(cmap)){
       cmap <- as.data.table(cmap)
   }
   mol.list <- splitMolecules.dt(cmap)
   #print("Done spliting with DT")
   #for(n in seq_along(mol.list)){
   #    mol.list[[n]] <- as.data.frame(mol.list[[n]])
   #}
   mol.list <- lapply(mol.list, as.data.frame)
   return(mol.list)
}


#split molecules according to their channel
getSeparateChannel <- function(mols){
  ref1 <- subset(mols, LabelChannel == 1 | LabelChannel == 0)
  ref2 <- subset(mols, LabelChannel == 2 | LabelChannel == 0)
  ref1 <- bookeepMolecules(ref1, sort=FALSE)
  ref2 <- bookeepMolecules(ref2, sort=FALSE)
  return(list(C1=ref1, C2=ref2))
}

#split cmap into list using DT
splitMolecules.dt <- function(cmap){
    mol.list <- cmap[,list(list(.SD)), by=CMapId, .SDcols=names(cmap)]$V1
    setattr(mol.list, 'names', unique(mol.list$CMapId))
}




#Given an alignmetn to query, for every position count how many molecules start
#and end at that particular position
getMoleculesStartEndCount <- function(cmap.ref, xmap, cmap){
  if(is.null((cmap.ref$MoleculeBeginCount))){
    cmap.ref$MoleculeBeginCount <- vector(mode='numeric', length=nrow(cmap.ref))
    cmap.ref$MoleculeEndCount <- vector(mode='numeric', length=nrow(cmap.ref))
    cmap.ref$Coverage <- 0
  }
  getInd <- createHashCMap(cmap.ref, "CMapId", "Position")
  beginInds <- paste0(xmap$RefcontigID, "|", xmap$RefStartPos)
  beginCount <- table(beginInds)
  inds <- getInd(names(beginCount))
  if(length(inds )< 1){
    warning("Cannot find corresponding indices in cmap for some positon in xmap ", head(names(beginCount)))
    return(cmap.ref)
  }
  cmap.ref$MoleculeBeginCount[inds] <- cmap.ref$MoleculeBeginCount[inds] + beginCount
  endInds <- paste0(xmap$RefcontigID, "|", xmap$RefEndPos)
  endCount <- table(endInds)
  inds <- getInd(names(endCount))
  cmap.ref$MoleculeEndCount[inds] <- cmap.ref$MoleculeEndCount[inds] + endCount
  cmap.ref$Coverage <- cmap.ref$Coverage + cmap$Coverage
  cmap.ref$Occurrence <- cmap.ref$Occurrence + cmap$Occurrence
  return(cmap.ref)
}

#compute the fraction of molecules that pass through a partricular label position
#for a reference cmap and all the molecuels aligned to it
computeMolPassThrough <- function(cmap, xmap){
  cmap.with.stat <- getMoleculesStartEndCount(cmap, xmap, cmap)
  Fragile.prob <- (cmap.with.stat$MoleculeBeginCount + cmap.with.stat$MoleculeEndCount) / cmap.with.stat$Coverage
  Fragile.prob[is.na(Fragile.prob)] <- 0
  Fragile.prob[Fragile.prob > 1] <- 1 #setting max at one
  cmap.with.stat$PassThrough <- 1 - Fragile.prob
  return(cmap.with.stat)
}


###################################################################################################################
#Section XMAP Utils
#Utility function for accessing and processing XMAP files

#convert the alignment str into a vector of integer where the pair of aligned siteInd
#are subsequent element in the vector
processAlignment <- function(alignmentStr){
  coords <- as.integer(unlist(strsplit(alignmentStr, "[(,)]")))
  coords <- coords[!is.na(coords)]
  return(coords);
}

#parse the xmap alignment string and convert the aligned sites to pair of indices
#in two column
processXmapAlignments <- function(xmap){
  align<-lapply(xmap$Alignment, processAlignment)
  alignRefInd <- lapply(align, function(x){return(x[seq(1,length(x), 2)])})
  alignQryInd <- lapply(align, function(x){return(x[seq(2,length(x), 2)])})
  #assignmennt doesn't seems to work directly when alignRefInd is one element list
  #browser()
  if(nrow(xmap) == 1){
    xmap$alignRefInd <- list(alignRefInd)
    xmap$alignQryInd <- list(alignQryInd)
  }else{
    xmap$alignRefInd <- alignRefInd
    xmap$alignQryInd <- alignQryInd
  }
  return(xmap)
}



#simple way of translate positions from one coordinates to another, only
#shift positions globally, do not really take into account the local distances
#The last argument tells whether one want to translate coordinates from 1: ref to qry
#or 2: qry to ref this is with repected to the alignment in xmap
translatePositions <- function(pos.list, align, direction){
  dists <- pos.list - pos.list[1]
  if(direction == 1){
    scaleFactor <- (align$QryEndPos - align$QryStartPos)/(align$RefEndPos - align$RefStartPos)
    leftOffset <- pos.list[1] - align$RefStartPos
    dists <- scaleFactor*(leftOffset + dists)
    return (dists + align$QryStart)
  }else if(direction == 2){
    scaleFactor <-(align$RefEndPos - align$RefStartPos)/ (align$QryEndPos - align$QryStartPos)
    leftOffset <- pos.list[1] - align$QryStartPos
    dists <- scaleFactor*(leftOffset + dists)
    return (dists + align$RefStart)

  }
}

#same as above but using the individual aligned labels pair to translate the coordinate instead of
#just a global offset. Direction 1: ref to qry 2: qry to ref
translatePositionsLocal <- function(pos.list, qcmap.dt, rcmap.dt, align.dt, direction){
  setkey(qcmap.dt, CMapId)
  setkey(rcmap.dt, CMapId)
  if(length(pos.list) != nrow(align.dt)){
    warning("Number of positions to be translated is not the same as the supplied alignments, returning original position list")
    return(pos.list)
  }
  if(is.null(align.dt$alignQryInd)){
    align.dt <- processXmapAlignments(align.dt)
  }
  translatePositionsHelper <- function(pos.list, qry.pos, ref.pos, align, direction){
    #browser()
    align.qry.Inds <- align$alignQryInd[[1]]
    align.ref.Inds <- align$alignRefInd[[1]]
    if(direction == 1){
      if(align$Orientation=='-'){
        align.qry.Inds <- rev(align$alignQryInd[[1]])
        align.ref.Inds <- rev(align$alignRefInd[[1]])
      }
      return(translatePositionsAlignToRef(pos.list, ref.pos, qry.pos, align$Orientation, align.ref.Inds, align.qry.Inds))
    }else if(direction == 2){
      return(translatePositionsAlignToRef(pos.list, qry.pos, ref.pos, align$Orientation, align$alignQryInd[[1]], align$alignRefInd[[1]]))
    }
  }
  translated.pos <- lapply(1:length(pos.list),
                           function(i){
                                qry.pos.list <- qcmap.dt[align.dt$QryContigID, Position]
                                ref.pos.list <- rcmap.dt[align.dt$RefcontigID, Position]
                                #browser()
                                if(length(qry.pos.list) > 1 & length(ref.pos.list) > 1){
                                  return(translatePositionsHelper(pos.list[[i]], qry.pos.list, ref.pos.list, align.dt[i], direction[[i]]))
                                }else{
                                  return(-1)
                                }
                            }
                          )

}



#translate positions from qry to a ref base on an alignment
translatePositionsAlignToRef <- function(pos.list, qry.pos, ref.pos, orientation, qry.alignInd, ref.alignInd){
  if(length(qry.alignInd) < 1 | length(ref.alignInd) < 1){
    stop("Cannot translate positions without alignment")
  }
  max.qry.pos <- max(qry.pos[qry.alignInd])
  min.qry.pos <- min(qry.pos[qry.alignInd])-0.1  #added a small padding in case there is one aligned inds (which should not really happen)
  max.ref.pos <- max(ref.pos[ref.alignInd])
  min.ref.pos <- min(ref.pos[ref.alignInd])-0.1
  scale <- (max.ref.pos - min.ref.pos)/(max.qry.pos - min.qry.pos)
  #pos.list may have some positions outside the boundary of alignment, compute their projected position boundary
  projected.max.qry <- max(c(max.qry.pos, pos.list)) + 1
  projected.min.qry <- min(c(min.qry.pos, pos.list)) - 1
  if(projected.max.qry > max.qry.pos + 1 | projected.min.qry < min.qry.pos - 1){
    warning(paste0("Qry position fall outside of aligned positions ", " " ,max.qry.pos, " ", min.qry.pos, " ", projected.max.qry, " ", projected.min.qry))
  }
  projected.max.ref <- (projected.max.qry - max.qry.pos)*scale + max.ref.pos
  projected.min.ref <- min.ref.pos - (min.qry.pos - projected.min.qry)*scale

  if(orientation == '+'){
      aligned.qry.pos <- c(projected.min.qry, qry.pos[qry.alignInd], projected.max.qry)
      aligned.ref.pos <- c(projected.min.ref, ref.pos[ref.alignInd], projected.max.ref)
  }else{
      aligned.qry.pos <- c(projected.min.qry, qry.pos[rev(qry.alignInd)], projected.max.qry)
      aligned.ref.pos <- c(projected.max.ref, ref.pos[rev(ref.alignInd)], projected.min.ref)
  }
  inds <- findInterval(pos.list, aligned.qry.pos)
  translated.pos <- translatePosToRef(pos.list, aligned.qry.pos[inds], aligned.qry.pos[inds+1],
                                      aligned.ref.pos[inds], aligned.ref.pos[inds+1])
  return(translated.pos)
}


#translate coordinate from qry positions to reference positions base on a pair of
#anchored positions left and right on qry and ref respectively
translatePosToRef <- function(qry.pos, qry.left, qry.right, ref.left, ref.right){
  scale <- (ref.right - ref.left)/(qry.right-qry.left)
  ref.pos <- (qry.pos - qry.left)*scale + ref.left
  return(ref.pos)
}

#check two alignment consistency
#This check whether molecules align to same region with two separate alignments
compareAlign <- function(dt.xmap1, dt.xmap2){
  if(nrow(dt.xmap1) > nrow(dt.xmap2)){
    dt.ref <- dt.xmap1
    dt.qry <- dt.xmap2
  }else{
    dt.ref <- dt.xmap2
    dt.qry <- dt.xmap1
  }
  numAligns <- nrow(dt.qry)
  setkey(dt.ref, QryContigID)
  check <- rep(TRUE, numAligns)
  for(n in 1:numAligns){
    qry.algn <- dt.qry[n,]
    ref.align <- dt.ref[.(qry.algn$QryContigID)]
    if(ref.align$RefStartPos > qry.algn$RefEndPos || ref.align$RefEndPos < qry.algn$RefStartPos){
      check[[n]] <- FALSE
    }
  }
}

#this function takes in two set of alignment to a common reference and
#for each alignment in xmap1, find the alignment in xmap2 that has the bigest overlap
findAlignOverlap <- function(xmap1.dt, xmap2.dt){
  setkey(xmap1.dt, RefcontigID, RefEndPos, RefStartPos)
  setkey(xmap2.dt, RefcontigID, RefEndPos, RefStartPos)
  findAlignOverlap.helper <- function(align, xmap2){
    overlap <- findOverlapContigsDT(align$RefStartPos, align$RefEndPos, align$RefcontigID, xmap2.dt, 100e3)
  }
  overlaps <- lapply(seq(1,nrow(xmap1.dt),1), function(ind){findAlignOverlap.helper(xmap1.dt[ind,], xmap2.dt)})
  return(overlaps)
}

#this function compuate the alignment rate of molecuels by molecule length and
#number of labels
getAlignmentRate <- function(xmap.dt, qcmap.dt, type='QryContigID', len.bins=50, numlabels.bins = 50, log.len=T, g1=NULL, g2=NULL, label='Sample', xlim1=c(4,6), xlim2=c(0,50), plot=T){
  mol.table <- qcmap.dt[SiteID==1]
  if(log.len){
    mol.table$ContigLength <- log10(mol.table$ContigLength)
  }
  #setkey(xmap.dt, 'QryContigID')
  xmap.dt <- unique(xmap.dt, by=force(type))
  setkeyv(xmap.dt, c(type))
  mol.table$is.aligned <- !is.na(xmap.dt[mol.table$CMapId, QryLen])
  h1.all <- hist(mol.table$ContigLength, len.bins, plot=F)
  h1.aligned <- hist(mol.table[is.aligned==TRUE, ContigLength], h1.all$breaks, plot=F)
  h2.all <- hist(mol.table$NumSites, numlabels.bins, plot=F)
  h2.aligned <- hist(mol.table[is.aligned==TRUE, NumSites], breaks=h2.all$breaks, plot=F)
  #browser()
  if(is.null(g1)){g1 <- ggplot()}
  g1 <- g1 + geom_line(aes(x=len, y=fract, color=label), data=data.frame(len=h1.all$mids, fract=h1.aligned$count/h1.all$count, label=label))
  if(is.null(g2)){g2 <- ggplot()}
  g2 <- g2 + geom_line(aes(x=len, y=fract, color=label), data=data.frame(len=h2.all$mids, fract=h2.aligned$count/h2.all$count, label=label))
  if(log.len){
    g1 <- g1 + xlab('Log10(Molecule length)')
  }else{
    g1 <- g1 + xlab('Molecule length')
  }
  g1 <- g1 + scale_x_continuous(limit=xlim1)
  g2 <- g2 + scale_x_continuous(limit=xlim2)
  g1 <- g1 + ylab('Fraction aligned')
  g2 <- g2 + xlab('Number of labels')
  g2 <- g2 + ylab('Fraction aligned')
  if(plot){
    plot(g1)
    plot(g2)
  }
  return(list(g1, g2, mol.table))
}

#this function takes in an xmap and project the begin and end coordinate on reference
#where the query is aligned
getProjAlignPos <- function(xmap){
  #browser()
  xmap[, ProjRefStart := ifelse(Orientation=='+', RefStartPos - QryStartPos, RefStartPos - (QryLen - QryStartPos))]
  xmap[, ProjRefEnd := ifelse(Orientation == '+', RefEndPos + (QryLen - QryEndPos), RefEndPos + QryEndPos)]
  return(xmap)
}




##############################################################################################################
######
#Visulalization Utils
#plot simulated molecuels along genome
plotSimulatedMolecules <- function(contigData){
    spacing <- 5;
    genomeY <- spacing*nrow(contigData)
    genomeX <- max(contigData$GenomeEnd)
    plot(c(1, genomeX), c(1, genomeY), type="n")
    lines(c(1,genomeX), c(genomeY, genomeY), type="l")
    for(n in 1:nrow(contigData)){
        #print(n)
        lines(c(contigData$GenomeStart[n],
                contigData$GenomeEnd[n]),
              c(genomeY-n*spacing,
                genomeY-n*spacing),
              type="l")
    }

}



##############################################################################################################
######
#Util function for dealing with interval objects
compareInterval <- function(start1, end1, start2, end2){
  if(start1 > end1){
    start1 <- tmp
    start1 <- end1
    end1 <- tmp
  }
  if(start2 > end2){
    start2 <- tmp
    start2 <- end2
    end2 <- tmp
  }
  dist <- 0
  type <- 0
  #intervals I1 and I2
  #I1 left of I2 no overlap
  if(start1 < start2){
    dist <- start2 - end1
    #otherwise I2 is completely embedded in I1
    if(-1*dist < end2 - start2){
        type <- -1
    }
  }else if(start2 < start1){
    dist <- start1 - end2
    if(-1*dist < end1 - start1){
      type <- 1
    }
  }
  return(c(dist, type))
}

#parrellel version that perform many-to-many intervals comparison
pCompareInterval <- function(start1, end1, start2, end2, proper.interval=F){
 #browser()
  if(!proper.interval){
    proper.interval <- getProperInterval(start1, end1)
    proper.interval2 <- getProperInterval(start2, end2)
    return(pCompareInterval(proper.interval$startPos, proper.interval$endPos,
                     proper.interval2$startPos, proper.interval2$endPos, T))
  }
  dist <- ifelse(start1 < start2, start2 - end1, start1 - end2)
  type1 <- ifelse(start1 < start2  & -1*dist < end2 - start2, -1, 0)
  type2 <- ifelse(start1 > start2 & -1*dist < end1 - start1, 1, 0)
  type <- type1 + type2
  embed.dist <- -1*pmin(end1-start1, end2-start2)
  dist <- ifelse(type==0, embed.dist, dist)
  return(data.frame(dist=dist, type=type))
}

#make sure interval is proper in the sense that all ends >= starts
getProperInterval <- function(starts, ends){
  start <- ifelse(starts < ends, starts, ends)
  end <- ifelse(starts > ends, starts, ends)
  return(data.frame(startPos=start, endPos=end))
}

#this is a generalized diff function
#it apply the fun to every consecutive pair of
#element in a list/array and return the value
generalized.diff <- function(x, fun=NULL){
  if(is.null(fun)){
    fun <- `-`
  }
  ele1 <- x[1:(length(x)-1)]
  ele2 <- x[2:(length(x))]
  ret <- fun(ele2, ele1)
  return(ret)
}

but.first <- function(l){
 return(l[-1])
}

but.last <- function(l){
  return(l[-length(l)])
}



#a generlaize version of above function, given a arbirary test function
#we swap the two columns that satisfy the condition given in test.fun
flip.columns <- function(dt, col1, col2, test.fun){
  dt <- as.data.frame(dt)
  selected <- test.fun(dt[,col1], dt[,col2])
  tmp <- dt[selected, col1]
  dt[selected, col1] <- dt[selected, col2]
  dt[selected, col2] <- tmp
  return(as.data.table(dt))
}


#make cmap consistent in columns
trim_cmap <- function(cmap1, cmap2){
  common_columns <- intersect(colnames(cmap1),colnames(cmap2)) #columns common between the two cmap input
  cmap1 <- cmap1[,colnames(cmap1)%in%common_columns]
  cmap2 <- cmap2[,colnames(cmap2)%in%common_columns]
  return(list(cmap1, cmap2))
}


