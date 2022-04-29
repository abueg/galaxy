#20151026 EL - simple filtering of indata

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

#filter pairs with weak molecule support; require a certain number of edges between two maps
filter_pairs <- function(data,min_edges=5) {

  pair_stats <- aggregate(list("num_edges"=rep(1,nrow(data))),list("id_str"=data[,"id_str"]),sum,na.rm=TRUE) #count how many molecules support each pair
  pair_ok <- pair_stats[pair_stats[,"num_edges"]>=min_edges,"id_str"] #filter based on min_edges

  F <- data[,"id_str"] %in% pair_ok
  cat(length(F),"pairs total, ",sum(F),"pairs pass,",sum(!F),"pairs discarded\n")

  data<-data[F,,drop=FALSE]

  return(data)
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
    r_sq <- 1
  }

  return(r_sq)
}

filter_by_correlation <- function(indata,thresh=0.9) {
  indata_split <- split(indata,indata[,"id_str"])

  indata_filtered <- mclapply(indata_split,function(indata) {
    pair_correlation <- get_correlation(indata)
    if (abs(pair_correlation)>thresh) {
      return(indata)
    } else {
      return(NULL)
    }
  },mc.cores=32)

  #indata_filtered <- do.call(rbind,indata_filtered)
  indata_filtered <- rbindlist(indata_filtered) #data.table verseion of rbind
  indata_filtered <- as.data.frame(indata_filtered,stringsAsFactors=F)

  nrows_indata <- nrow(indata)
  nrows_indata_filtered <- nrow(indata_filtered)
  nrows_discarded <- nrows_indata-nrows_indata_filtered

  cat(nrows_indata,"pairs total, ",nrows_indata_filtered,"pairs pass,",nrows_discarded,"pairs discarded\n")

  return(indata_filtered)
}

filter_by_regression_slope <- function(indata,thresh=c(0.9,1.1)) {
  stopifnot(length(thresh)==2)

  indata_split <- split(indata,indata[,"id_str"])

  indata_filtered <- mclapply(indata_split,function(indata) {
    rlm_out <- rlm(indata[,"offset.y"]~indata[,"offset.x"])
    slope <- coef(summary(rlm_out))[2,1]

    if (abs(slope)>thresh[1] & abs(slope)<thresh[2]) {
      return(indata)
    } else {
      return(NULL)
    }
  },mc.cores=32)

  #indata_filtered <- do.call(rbind,indata_filtered)
  indata_filtered <- rbindlist(indata_filtered) #data.table verseion of rbind
  indata_filtered <- as.data.frame(indata_filtered,stringsAsFactors=F)

  nrows_indata <- nrow(indata)
  nrows_indata_filtered <- nrow(indata_filtered)
  nrows_discarded <- nrows_indata-nrows_indata_filtered

  cat(nrows_indata,"pairs total, ",nrows_indata_filtered,"pairs pass,",nrows_discarded,"pairs discarded\n")

  return(indata_filtered)
}

collapse_duplicates <- function(data) { #collapse the single-molecule alignment data and output pairs of maps

  id_str <- data[,"id_str"]
  sorted <- sort(table(id_str),decreasing=T) #get count summary (sorted)

  names_split <- lapply(names(sorted),function(name) {
    name_split <- unlist(strsplit(name,"_"))
  })
  names_split <- do.call(rbind,names_split)

  collapsed <- data.frame("id_str"=names(sorted),"asm1"=names_split[,1],"asm2"=names_split[,2],"num_edges"=as.numeric(sorted),stringsAsFactors=F)

  return(collapsed)
}

fix_numbering <- function(data) { #make sure no duplicates in mapIDs (the two assemblies might have the same mapIDs)

  asm1_uniq <- sort(as.numeric(unique(data[,"asm1"])))
  asm2_uniq <- sort(as.numeric(unique(data[,"asm2"])))

  #plot(1:length(asm1_uniq),asm1_uniq) #show that there is a need for numbering

  lookup1 <- data.frame("asm1"=asm1_uniq,"asm1_renumbered"=1:length(asm1_uniq))
  lookup2 <- data.frame("asm2"=asm2_uniq,"asm2_renumbered"=length(asm1_uniq)+(1:length(asm2_uniq)))

  data1 <- merge(data,lookup1,by="asm1")
  stopifnot(nrow(data1)==nrow(data))

  data2 <- merge(data1,lookup2,by="asm2")
  stopifnot(nrow(data2)==nrow(data1))

  return(data2)
}
