#20151013 EL - plotting utilities

#given the path coordinates, make a quick schematic of the path
#input: the path coordinates
#output: output list - min/max position and the length of path
plot_schematic <- function(path_coordinates,plot_flag=F,order_flag=F,title_str='') { #plot a schematic of what the path looks like

  if (order_flag) { #whether to order by start
    path_coordinates <- path_coordinates[order(path_coordinates[,"rel_start"]),]
  }

  min_pos <- min(path_coordinates[,"rel_start"],path_coordinates[,"rel_end"]) #get xrange
  max_pos <- max(path_coordinates[,"rel_start"],path_coordinates[,"rel_end"])
  total_length <- max_pos - min_pos

  if (plot_flag) { #whether to plot
    nrows <- nrow(path_coordinates)

    plot(0,type="n",xlab="Nucleotide position",ylab="",yaxt='n',bty='n',xlim=c(min_pos,max_pos),ylim=c(0,nrows+1),main=paste0("path schematic"),sub=title_str) #make empty plot
    axis(2,at=1:nrows) #add y-axis

    exit_sig <- lapply(1:nrow(path_coordinates),function(j) { #plot each node
      entry <- path_coordinates[j,,F]
      midpoint <- mean(c(entry[,"rel_start"],entry[,"rel_end"]))

      if (entry[,"node_type"]=="asm1") { #draw map lines
        lines(c(entry[,"rel_start"],entry[,"rel_end"]),c(j,j),col=rgb(1,0,0,0.99),lend=0,lwd=0.5) #plot each match line
      } else {
        lines(c(entry[,"rel_start"],entry[,"rel_end"]),c(j,j),col=rgb(0,0,1,0.99),lend=0,lwd=0.5) #plot each match line
      }

      if (entry[,"rel_start"]<entry[,"rel_end"]) { #draw arrows indicating orientation
        arrows(x0=entry[,"rel_end"]-total_length*0.02,x1=entry[,"rel_end"],length=0.05,angle=25,y0=j)
      } else {
        arrows(x0=entry[,"rel_end"]+total_length*0.02,x1=entry[,"rel_end"],length=0.05,angle=25,y0=j)
      }

      text(midpoint,j+0.1,labels=entry[,"node"],cex=0.75,srt=0)
    })

    legend(max_pos-(total_length*0.2),nrows*0.2,c("asm1","asm2"),lty=c(1,1),lwd=c(1.5,1.5),col=c(rgb(1,0,0,0.99),rgb(0,0,1,0.99)),cex=0.8)
  }

  outlist <- list("min_pos"=min_pos,
                  "max_pos"=max_pos,
                  "total_length"=total_length)

  return(outlist)
}

#given a 2-color cmap (it should handle a one-color cmap as well), make a plot the label distributions (for quick diagnostics)
#input: cmap dataframe
#output: void
plot_cmap <- function(cmap,plot_flag=F) {
  if (!plot_flag) {
    return()
  }

  xrange <- range(cmap[,"Position"])
  yrange <- range(cmap[,"LabelChannel"])
  #   print(xrange)
  #   print(yrange)

  pos1 <- cmap[cmap[,"LabelChannel"]==1,"Position"]
  pos2 <- cmap[cmap[,"LabelChannel"]==2,"Position"]

  plot(0,type="n",xlim=xrange,ylim=yrange)
  try(points(pos1,rep(1,length(pos1)),col="blue",pch="|"))
  try(points(pos2,rep(2,length(pos2)),col="red",pch="|"))

  return()
}

#given a list of cmap dataframes, make a quick schematic
#input: cmap list
#output: void
plot_cmap_individaul <- function(cmap_list,col,yoffset=0,total_length) {
  exit_sig <- lapply((yoffset+1):(length(cmap_list)+yoffset),function(i) {
#     print(i-yoffset)
    entry <- cmap_list[[i-yoffset]]
    midpoint <- mean(entry[["Position"]])

    lines(range(entry[,"Position"]),c(i,i),col=col,lend=0,lwd=0.5) #plot each match line

    if (entry[1,"SiteID"]<tail(entry[,"SiteID"],1)) { #draw arrows indicating orientation
#       print("positive orientation")
      arrows(x0=max(entry[,"Position"])-total_length*0.02,x1=max(entry[,"Position"]),length=0.05,angle=25,y0=i)
    } else {
#       print("negative orientation")
      arrows(x0=min(entry[,"Position"])+total_length*0.02,x1=min(entry[,"Position"]),length=0.05,angle=25,y0=i)
    }

    text(midpoint,i+0.1,labels=entry[1,"CMapId"],cex=0.75,srt=0)
  })

  return()
}

#wrapper for plot_cmap_individaul
#input: cmap list (separated by color)
#output: void
plot_cmaps <- function(cmap_list,plot_flag=F,title_str='') {
  if (!plot_flag) {
    return()
  }

  asm1 <- cmap_list[["asm1"]]
  asm2 <- cmap_list[["asm2"]]

#   print(head(lapply(asm1,head)))
#   print(head(lapply(asm2,head)))

  length_asm1 <- length(asm1)
  length_asm2 <- length(asm2)
  nrows <- length_asm1 + length_asm2

  #get xrange
  min_pos <- min(unlist(lapply(asm1,function(x) min(x[,"Position"]))),unlist(lapply(asm2,function(x) min(x[,"Position"]))))
  max_pos <- max(unlist(lapply(asm1,function(x) max(x[,"Position"]))),unlist(lapply(asm2,function(x) max(x[,"Position"]))))
  total_length <- max_pos-min_pos

#   print(min_pos)
#   print(max_pos)

  plot(0,type="n",xlab="Nucleotide position",ylab="",bty='n',xlim=c(min_pos,max_pos),ylim=c(0,nrows),main=paste0("cmap schematic"),sub=title_str) #make empty plot

  plot_cmap_individaul(asm1,col=rgb(1,0,0,0.99),yoffset=0,total_length=total_length) #plot individual cmap components
  plot_cmap_individaul(asm2,col=rgb(0,0,1,0.99),yoffset=length_asm1,total_length=total_length)

  legend(max_pos-(total_length*0.2),3,c("asm1","asm2"),lty=c(1,1),lwd=c(1.5,1.5),col=c(rgb(1,0,0,0.99),rgb(0,0,1,0.99)),cex=0.8) #add legend
}

# #20151014 EL - testing plot_cmaps
# load(file="/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/Human_willsample_twocolor_cmap_list.RData")
# plot_cmaps(cmap_list,plot_flag=T,title_str='path 2')

#make simple plot of alignments
#input: xmap
#output: void if plot_flag is false; otherwise, xrange and total length of region covered
plot_alignment_schematic <- function(xmap,pad=50000,title_str='',plot_flag=F) {
  if (!plot_flag) {
    return()
  }

#   print(xmap[,c(1:9,11:12)])
  xrange <- range(xmap[,"RefStartPos"],xmap[,"RefEndPos"])
  xrange <- c(xrange[1]-pad,xrange[2]+pad)
#   print(xrange)

  total_length <- diff(xrange)

  nrows <- nrow(xmap)

  plot(0,type="n",xlab="Nucleotide position",ylab="",bty='n',xlim=xrange,ylim=c(0,nrows),main=paste0("alignment schematic"),sub=title_str) #make empty plot

  exit_sig <- lapply(1:nrows,function(i) {
    entry <- xmap[i,,F]
    midpoint <- mean(c(entry[,"RefStartPos"],entry[,"RefEndPos"]))

    lines(c(entry[,"RefStartPos"],entry[,"RefEndPos"]),c(i,i),col=rgb(1,0,0,0.99),lend=0,lwd=0.5) #plot each match line

    if (entry[,"Orientation"]=="+") { #draw arrows indicating orientation
      arrows(x0=entry[,"RefEndPos"]-total_length*0.02,x1=entry[,"RefEndPos"],length=0.05,angle=25,y0=i)
    } else {
      arrows(x0=entry[,"RefStartPos"]+total_length*0.02,x1=entry[,"RefStartPos"],length=0.05,angle=25,y0=i)
    }

    text(midpoint,i+0.1,labels=entry[,"QryContigID"],cex=0.75,srt=0)
  })

  abline(v=xmap[,"RefStartPos"],col="red")
  abline(v=xmap[,"RefEndPos"],col="blue")

  outlist <- list("xrange"=xrange,
                  "total_length"=total_length)

  return(outlist)
}

plot_alignment_schematic_w_mg <- function(xmap,pad=100000,title_str='',plot_flag=F) {
  if (!plot_flag) {
    return()
  }

  print(xmap[,c(1:9,11:12)])
  xrange <- range(xmap[,"RefStartPos"],xmap[,"RefEndPos"])
  xrange <- c(xrange[1]-pad,xrange[2]+pad)
  #   print(xrange)

  total_length <- diff(xrange)

  nrows <- nrow(xmap)

  plot(0,type="n",xlab="Nucleotide position",ylab="",yaxt='n',bty='n',xlim=xrange,ylim=c(0,nrows),main=paste0("alignment schematic"),sub=title_str) #make empty plot
  #axis(2,at=c(0:nrows))
  axis(2,at=c(0:nrows),labels=c("ref",xmap[,"QryContigID"]),las=1)

  abline(h=0)

  exit_sig <- lapply(1:nrows,function(i) {
    entry <- xmap[i,,F]

    q_start <- entry[,"QryStartPos"]
    q_end <- entry[,"QryEndPos"]
    r_start <- entry[,"RefStartPos"]
    r_end <- entry[,"RefEndPos"]
    orientation <- entry[,"Orientation"]

    if (orientation=="+") { #flip positions
      offset <- r_start - q_start
    } else {
      offset <- r_start + q_start
    }

    map_endpoints <- c(0,entry[,"QryLen"])+offset
    lines(map_endpoints,c(i,i),col=rgb(1,0,0,0.99),lend=0,lwd=0.5) #plot map

    lines(c(q_start+offset,r_start),c(i,0)) #plot match lines
    lines(c(q_end+offset,r_end),c(i,0))

    if (orientation=="+") { #draw arrows indicating orientation
      arrows(x0=map_endpoints[2]-total_length*0.02,x1=map_endpoints[2],length=0.05,angle=25,y0=i)
    } else {
      arrows(x0=map_endpoints[1]+total_length*0.02,x1=map_endpoints[1],length=0.05,angle=25,y0=i)
    }

    midpoint <- mean(map_endpoints)
    text(midpoint,i+0.1,labels=entry[,"QryContigID"],cex=0.75,srt=0) #label map ids
  })

  outlist <- list("xrange"=xrange,
                  "total_length"=total_length)

  return(outlist)
}

plot_alignment_schematic_w_mg_2col <- function(output_ra_alignref1,output_ra_alignref2,title_str='',pad=100000,plot_flag=F) {
  if (!plot_flag) {
    return()
  }

  prefix <- output_ra_alignref1[["prefix"]]
  xmap1 <- readxmap(paste0(prefix,".xmap"))
  xmap1 <- flip(xmap1)

  prefix <- output_ra_alignref2[["prefix"]]
  xmap2 <- readxmap(paste0(prefix,".xmap"))

  nrows1 <- nrow(xmap1)
  nrows2 <- nrow(xmap2)

  nrows <- nrows1+nrows2+1 #the extra 1 is for the reference line

  xrange <- range(xmap1[,"RefStartPos"],xmap1[,"RefEndPos"],xmap2[,"RefStartPos"],xmap2[,"RefEndPos"])
  xrange <- c(xrange[1]-pad,xrange[2]+pad)

  total_length <- diff(xrange)

  plot(0,type="n",xlab="Nucleotide position",ylab="",yaxt='n',bty='n',xlim=xrange,ylim=c(1,nrows),main=paste0("alignment schematic"),sub=title_str) #make empty plot
  #axis(2,at=c(0:nrows))
  axis(2,at=c(1:nrows),labels=c(xmap1[,"QryContigID"],"ref",xmap2[,"QryContigID"]),las=1)

  y_ref <- nrows1+1
  abline(h=y_ref)

  xmap1[,"asm"] <- 1
  xmap2[,"asm"] <- 2

  xmap1[,"index"] <- 1:nrows1
  xmap2[,"index"] <- (1:nrows2)+nrows1+1

  xmap <- rbind(xmap1,xmap2)

  exit_sig <- lapply(1:nrow(xmap),function(i) {
    entry <- xmap[i,,F]

    q_start <- entry[,"QryStartPos"]
    q_end <- entry[,"QryEndPos"]
    r_start <- entry[,"RefStartPos"]
    r_end <- entry[,"RefEndPos"]
    orientation <- entry[,"Orientation"]
    asm <- entry[,"asm"]
    index <- entry[,"index"]

    if (orientation=="+") { #calculate offsets
      offset <- r_start - q_start
    } else {
      offset <- r_start + q_start
    }

    if (asm==1) { #decide on color
      color <- rgb(1,0,0,0.8)
    } else {
      color <-rgb(0,1,0,0.8)
    }

    map_endpoints <- c(0,entry[,"QryLen"])+offset
    lines(map_endpoints,c(index,index),col=color,lend=0,lwd=0.5) #plot map

    lines(c(q_start+offset,r_start),c(index,y_ref),col=rgb(0,0,0,0.2)) #plot match lines
    lines(c(q_end+offset,r_end),c(index,y_ref),col=rgb(0,0,0,0.2))

    if (orientation=="+") { #draw arrows indicating orientation
      arrows(x0=map_endpoints[2]-total_length*0.02,x1=map_endpoints[2],length=0.05,angle=25,y0=index)
    } else {
      arrows(x0=map_endpoints[1]+total_length*0.02,x1=map_endpoints[1],length=0.05,angle=25,y0=index)
    }

    midpoint <- mean(map_endpoints)
    text(midpoint,index+0.1,labels=entry[,"QryContigID"],cex=0.75,srt=0) #label map ids
  })

  outlist <- list("xrange"=xrange,
                  "total_length"=total_length)

  return(outlist)
}

get_ref_alignment <- function(path_data,asm_id,real_xmap) {
  #asm_observed_data <- real_xmap[(real_xmap[,"QryContigID"]%in%asm_id),c("QryContigID","RefStartPos","RefEndPos")] #get observed data
  asm_observed_data <- real_xmap[(real_xmap[,"QryContigID"]%in%asm_id),] #get observed data
  asm_predited_data <- path_data[(path_data[,"node"]%in%asm_id),c("node","rel_start","rel_end")] #get path data
  asm_merged <- merge(asm_observed_data,asm_predited_data,by.x = "QryContigID", by.y ="node") #merge

  return(asm_merged)
}

check_predictions <- function(path_data,real_xmap1,real_xmap2,plot_flag=F,title_str=''){
  #   print(real_xmap1)
  #   print(real_xmap2)

  real_xmap1 <- readxmap(real_xmap1) #based on alignment of genomes to reference; this is the "answer"
  real_xmap2 <- readxmap(real_xmap2)

  #   print(head(real_xmap1[,1:9]))
  #   print(head(real_xmap2[,1:9]))
  #browser()
  path_data_asm1 <- subset(path_data, node_type=='asm1')
  path_data_asm2 <- subset(path_data, node_type=='asm2')
  #get data for asm1
  asm1_id <- path_data$node[which(as.character(path_data$node_type)=="asm1")]
  asm1_merged <- get_ref_alignment(path_data_asm1,asm1_id,real_xmap1)

  #   print("asm1")
  #   print(asm1_merged)

  if (nrow(asm1_merged)>0) {
    asm1_merged[,"color"] <- "red"
  }

  #get data for asm2
  asm2_id <- path_data$node[which(as.character(path_data$node_type)=="asm2")]
  asm2_merged <- get_ref_alignment(path_data_asm2,asm2_id,real_xmap2)

  if (nrow(asm2_merged)>0) {
    asm2_merged[,"color"] <- "blue"
  }

  #   print("asm2")
  #   print(asm2_merged)

  if (nrow(asm1_merged)>0 && (nrow(asm2_merged)>0)) {
    coordinates_answer_sim <- rbind(asm1_merged,asm2_merged)
  } else if (nrow(asm1_merged)>0) {
    coordinates_answer_sim <- asm1_merged
  } else {
    coordinates_answer_sim <- asm2_merged
  }

  major_chr <- find_mode(coordinates_answer_sim[,"RefcontigID"]) #making sure all the alignments are on the same chromosome
  coordinates_answer_sim <- coordinates_answer_sim[coordinates_answer_sim[,"RefcontigID"]==major_chr,] #subset

  if (nrow(coordinates_answer_sim)>2) {
    predicted_midpoints <- (coordinates_answer_sim[,"rel_start"]+coordinates_answer_sim[,"rel_end"])/2
    observed_midpoints <- (coordinates_answer_sim[,"RefStartPos"]+coordinates_answer_sim[,"RefEndPos"])/2

    if (plot_flag) {
      plot(predicted_midpoints,observed_midpoints,col=coordinates_answer_sim[,"color"],main=title_str,xlab="Predicted midpoints",ylab="Observed midpoints (alignment-based)",pch=20)
      text(predicted_midpoints,observed_midpoints,labels=coordinates_answer_sim[,"QryContigID"],cex=0.6,pos=3,offset=0.2)
    }
  } else {
    return(NA)
  }

  #correlation
  correlation <- (cor(predicted_midpoints,observed_midpoints))^2

  #linear regression
  lm.out <- lm(observed_midpoints~predicted_midpoints)
  if (plot_flag) {
    try(abline(lm.out,col=rgb(0,0,0,0.3),lty="dotted"))
  }

  #   print(summary(lm.out))
  #   print(coef(summary(lm.out)))

  if (plot_flag) {
    slope <- 0
    if (nrow(coef(summary(lm.out)))==2) {
      slope <- coef(summary(lm.out))[2,1]
    }

    if (slope>0) {
      legend("bottomright",bty="n",legend=paste0("correlation = ",format(correlation,digits=4),"\n",
                                                 "slope = ",format(slope,digits=4)),cex=0.8)
    } else {
      legend("topright",bty="n",legend=paste0("correlation = ",format(correlation,digits=4),"\n",
                                              "slope = ",format(slope,digits=4)),cex=0.8)
    }
  }

  residuals <- summary(lm.out)$residuals
  #print(path_data)
  #adding coordinate accuracy informaton for output
  #coordinates_stat <- data.frame(node=asm1_id, node_type="asm1")
  #coordinates_stat <- rbind(coordinates_stat, data.frame(node=asm2_id, node_type="asm2"))
  #coordinates_stat$residuals <- residuals

  #path_data <- merge(path_data, coordinates_stat, by=c("node", "node_type"))

  #return(list(correlation=correlation, residuals=path_data))
  return(correlation)
  return(residuals)
}

#checking the transformation from start-to-start distance to center-to-center distance
compare_offsets <- function(pairs_merged) {
#   print(pairs_merged)

  #par(mfrow=c(2,2))
  lapply(1:nrow(pairs_merged),function(i) {
    pair <- pairs_merged[i,,F]

    id1 <- pair[,"RefcontigID.x"]
    id2 <- pair[,"RefcontigID.y"]

    len1 <- pair[,"RefLen.x"]
    len2 <- pair[,"RefLen.y"]

    intercept <- as.numeric(pair[,"intercept"])
    slope <- ((as.numeric(pair[,"slope"])>0)*2)-1 #convert to whole numbers

    #simple plot that shows the regression line
    plot(c(0,len1),c(0,len2),type='n',xlab=paste0("Map ",id1," coordinates"),ylab=paste0("Map ",id2," coordinates"))
    abline(intercept,slope)

    legend("bottomright",bty="n",legend=paste0("intercept = ",format(intercept,digits=2,big.mark=","),"\n",
                                               "slope = ",format(slope,digits=2)),cex=0.8)

    #offset schematic
    xstart1 <- 0
    xend1 <- len1

    xstart2 <- -1*intercept*slope
    xend2 <- xstart2+slope*len2

#     print(paste0("xstart1: ",format(xstart1,digits=2),
#                  "; xend1: ",format(xend1,digits=2),
#                  "; xstart2: ",format(xstart2,digits=2),
#                  "; xend2: ",format(xend2,digits=2)))

    xmin <- min(xstart1,xstart2,xend1,xend2)
    xmax <- max(xstart1,xstart2,xend1,xend2)

    plot(c(xmin,xmax),c(0,3),type='n',yaxt='n',xlab="coordinates",ylab=paste0(""))
    axis(2,at=c(1,2),labels=c(paste0("Map",id1),paste0("Map",id2)),las=1) #draw y-axis

    lines(c(xstart1,xend1),c(1,1)) #map x
    lines(c(xstart2,xend2),c(2,2)) #map y

    text(xstart1,1,"S1",cex=0.75)
    text(xstart2,2,"S2",cex=0.75)
    text(xend1,1,"E1",cex=0.75)
    text(xend2,2,"E2",cex=0.75)

    #start-to-start
    lines(c(0,0),c(1,2),lty="dotted")
    if (xstart2>0) {
      lines(c(0,xstart2),c(2,2),lty="dotted")
    }
    points(x=0,y=1,col="blue",pch=16)
    points(x=0,y=2,col="blue",pch=16)

    #center-to-center
    if (slope==1) {
      lines(c(len1/2,(len2/2)+xstart2),c(1,2),lty="dotted")
      points(x=(len2/2)+xstart2,y=2,col="red",pch=16)
    } else {
      lines(c(len1/2,(len2/2)+xend2),c(1,2),lty="dotted")
      points(x=(len2/2)+xend2,y=2,col="red",pch=16)
    }

    points(x=len1/2,y=1,col="red",pch=16)

  })
}

#http://stats.stackexchange.com/questions/83826/is-a-weighted-r2-in-robust-linear-model-meaningful-for-goodness-of-fit-analys
r2 <- function(x){
  SSe <- sum((x$resid)^2)
  observed <- x$resid+x$fitted
  SSt <- sum((observed-mean(observed))^2)
  value <- 1-SSe/SSt
  return(value)
}

r2ww <- function(x){ #weighted r-squared from the robust linear model
  SSe <- sum((x$w*x$resid)^2) #the residual sum of squares is weighted
  observed <- x$resid+x$fitted
  SSt <- sum((x$w*(observed-mean(observed)))^2) #the total sum of squares is weighted
  value <- 1-SSe/SSt
  return(value)
}

#alignment coverage of the map (molecule pileup)
plot_coverage <- function(data,id_str,plot_flag=FALSE) {

  reflen <- data[1,"RefLen"]
  nrows <- nrow(data)

  data <- data[order(data[,"RefStartPos"]),] #order by RefStartPos

  if (nrows>1) {
    nalignment_str <- paste0(nrows," alignments")
  } else {
    nalignment_str <- paste0(nrows," alignment")
  }

  span <- range(data[,"RefStartPos"],data[,"RefEndPos"]) #rough estimate of span
  span_fraction <- diff(span)/reflen*100 #percentage

  data[,"midpoint"] <- (data[,"RefStartPos"]+data[,"RefEndPos"])/2

  if (plot_flag) {
    plot(0,type="n",xlab="Nucleotide position",ylab="",yaxt='n',bty='n',xlim=c(0,reflen),ylim=c(0,nrows+1), #make empty plot
         main=paste0("pair: ",id_str),
         sub=paste0(nalignment_str,"; ","span_fraction: ",round(span_fraction,2),"%"))

    lines(c(0,reflen),c(0.5,0.5),col="red",lend=1,lwd=2) #map line

    exit_sig <- lapply(1:nrow(data),function(i) {
      entry <- data[i,,F]

      lines(c(entry[,"RefStartPos"],entry[,"RefEndPos"]),c(i,i),col="black",lend=0,lwd=1) #plot each alignment
      #text(data[,"midpoint"][i],i+0.5,labels=paste0(entry[,"QryContigID"]),cex=0.1)
    })
  }

  return(span_fraction)
}

#linear regression; center-to-center offset
get_offset_single <- function(data){
    data <- as.data.table(data)
    sign.x <- (data$Orientation.x == '+')*2 - 1
    sign.y <- (data$Orientation.y == '+')*2 - 1
    offset <- data[,sign.y * ((QryStartPos.x - sign.x * RefStartPos.x ) - (QryStartPos.y - sign.y* RefStartPos.y))]
    return(mean(offset))
}

plot_correlation1 <- function(data,plot_flag=FALSE) {
   try({
     if(nrow(data) == 1){
        intercept <- get_offset_single(data)
        slope <- data$rel_orient[[1]]
        #TODO fix this single-case
        return(c(intercept, slope))
    }else{

        rlm_out <- rlm(data[,"offset.y"]~data[,"offset.x"])
        intercept <- coef(summary(rlm_out))[1,1]
        slope <- coef(summary(rlm_out))[2,1]
        rlm_resid <- resid(rlm_out)
    }
    outlier <- which((abs(rlm_resid-median(rlm_resid))/mad(rlm_resid)) > 50)

    if (length(outlier)==0) {
      correlation <- cor(data[,"offset.x"],data[,"offset.y"])
    }else{
      correlation <- cor(data[-outlier,"offset.x"],data[-outlier,"offset.y"])
    }
  })

  if (plot_flag) {
    plot(data[,"offset.x"],data[,"offset.y"],pch=16,main="correlation between x and y",xlab="x coordinates",ylab="y coordinates")
    try(abline(rlm_out,col="red"))
    text(min(data[,"offset.x"]),min(data[,"offset.y"]),label=round(correlation/.02)*.02)

    #plot(data[,"offset.x"],lm_resid,pch=16,main="residual plot",xlab="x coordinates",ylab="residual") #residual plot
    #abline(0, 0)

    #out_x <- data[,"offset.x"][lm_resid %in% outlier]
    #points(out_x,outlier,col="red")
  }

  final_slope <- get_closer_to_one(slope,correlation)

  len_x <- data[1,"RefLen.x"]
  len_y <- data[1,"RefLen.y"]

  slope <- final_slope>0
  slope <- slope*2-1

  #intercept <- intercept+len_y/2-slope*len_x/2

  #using simple offset, not usign slope
  dir <- ifelse(data[,'Orientation.x']=='+', 1, -1)
  final_slope <- ifelse(data[,'Orientation.x'][[1]] == data[,'Orientation.y'][[1]], 1, -1)
  intercept <- mean(c(data[,'RefStartPos.y'] - final_slope*(data[,'RefStartPos.x'] - (data[,'QryStartPos.x'] - data[,'QryStartPos.y'])*dir),
                      data[,'RefEndPos.y'] - final_slope*(data[,'RefEndPos.x'] - (data[,'QryEndPos.x'] - data[,'QryEndPos.y'])*dir)
                                                          ))
  return(c(intercept,final_slope))
}

#alignment between two maps (each alignment line represent a molecule connecting the two consensus maps)
plot_alignment1 <- function(data,plot_flag=F) {
  if (!plot_flag) {
    return()
  }

  len_x <- data[1,"RefLen.x"]
  len_y <- data[1,"RefLen.y"]

  xlim <- c(0,max(len_x,len_y))

  plot(0,type="n",xlab="Nucleotide position",ylab="",yaxt='n',bty='n',xlim=xlim,ylim=c(0,1), #make empty plot
       main=paste0("pair: ",data[,"id_str"][1]))

  lines(c(0,len_x),c(0,0),col="red",lend=1,lwd=2) #map1 line
  lines(c(0,len_y),c(1,1),col="blue",lend=1,lwd=2) #map2 line

  exit_sig <- lapply(1:nrow(data),function(i) {
    entry <- data[i,,F]

    lines(c(entry[,"offset.x"]+len_x/2,entry[,"offset.y"]+len_y/2),c(0,1),col=rgb(1,0,1,0.25),lend=0,lwd=0.25) #plot each match line
  })
}
