#20150630 EL - trying to do a 3-layer view ref-genome map-molecule alignment

# rm(list=ls())
#
# scripts_dir <- "/home/users/elam/rscripts/"
# source(paste(scripts_dir,"paste0.R",sep=""))
# source(paste0(scripts_dir,"readmaps.R")) #util functions to read _map formats
# source(paste0(scripts_dir,"overlap.R"))
# source(paste0(scripts_dir,"cigar1.R"))
#
# options(scipen=999)

#draw reference line based on provided xrange
draw_reference <- function(rcmap,xrange,nrows,pad=50000,region_str='') {

  ref_id <- rcmap[1,"CMapId"]

  #set range (and make sure they are within bounds)
  xrange <- c(xrange[1]-pad,xrange[2]+pad) #apply padding
  xrange[1] <- max(xrange[1],0) #truncate range (start coord cannot be less than 0)
  xrange[2] <- min(xrange[2],rcmap$ContigLength[1]) #truncate range (start coord cannot be more than reference length)

  rcmap_selected <- rcmap[rcmap$Position>=xrange[1] & rcmap$Position<=xrange[2],] #get rcmap entries

  #make empty plot
  #plot(xrange,c(0,nrows),type="n",axes=F,ann=T,xlab="Coordinates",ylab="Molecule indices",main=region_str,sub=paste0("#alignments: ",nrows/2))
  plot(xrange,c(0,nrows),type="n",axes=F,ann=T,xlab="Coordinates",ylab="Molecule indices",main=region_str,sub=paste0("chr",ref_id))

  #draw reference line
  ycoord_ref <- 0.2
  color_ref <- rgb(0,0,0)
  lines(xrange,c(ycoord_ref,ycoord_ref),lwd=2,col=color_ref) #reference line

  #split channels
  rcmap_selected <- split(rcmap_selected,rcmap_selected["LabelChannel"])
  #stopifnot(length(rcmap_selected)==2)

  pos_ref_ch1 <- c()
  pos_ref_ch2 <- c()

  try({
    #plot channel 1
    pos_ref_ch1 <- rcmap_selected[['1']]$Position #prepare reference labels
    nsites_ch1 <- length(pos_ref_ch1)
    color_ref <- rgb(1,0,0,0.5)
    points(pos_ref_ch1,rep(ycoord_ref+0.02,nsites_ch1),xlab="",ylab="",pch="l",col=color_ref,cex=0.75) #reference label positions

  })

  try({
    #plot channel 2
    pos_ref_ch2 <- rcmap_selected[['2']]$Position
    nsites_ch2 <- length(pos_ref_ch2)
    color_ref <- rgb(0,1,0,0.5)
    points(pos_ref_ch2,rep(ycoord_ref-0.02,nsites_ch2),xlab="",ylab="",pch="l",col=color_ref,cex=0.75) #reference label positions

  })

  #   offset_labelID <- 5 #label IDs
  #   text(c(pos_ref[1],tail(pos_ref,1)),rep(ycoord_ref+offset_labelID,2), #first and last label ids
  #        lab=c(rcmap_selected$SiteID[1],tail(rcmap_selected$SiteID,1)),
  #        cex=0.5,col=color_ref)

  mapid <- rcmap$CMapId[1]
  # mtext(paste0("map",mapid," coordinates (bp)"),side=1,line=3) #x-axis

  #add axis
  xlabels <- round(xrange/1000) #x-axis labels
  xsteps <- diff(xlabels)/5 #make intervals
  xlabels <- seq(from=xlabels[1],to=xlabels[2],by=xsteps)*1000
  xtext <- format(xlabels,big.mark=",") #can add formatting (current not doing anything)

  axis(side=1,at=xlabels,labels=xtext,tick=TRUE,line=NA, #axis
       pos=NA,outer=FALSE,font=NA,lty="solid",
       lwd=1,lwd.ticks=1,col=NULL,col.ticks=NULL,
       hadj=NA,padj=NA,
       cex.axis=1)

  arrows(x0=xrange[2]-diff(xrange)/1000,x1=xrange[2],y0=ycoord_ref,angle=20,length=0.12,col="red",lwd=2) #arrow at the end

  outlist <- list("pos_ref_ch1"=pos_ref_ch1,
                  "pos_ref_ch2"=pos_ref_ch2)

  return(outlist)
  return(0)
}

#draw molecules
draw_molecules <- function(qcmap,xmap,stretch_flag=F,colors,yoffset=0,xoffset=NULL) {
  stopifnot(length(unique(xmap[,"RefcontigID"]))==1) #should only align to one reference contig
  stopifnot(nrow(xmap)==1) #expect single alignment for each map

  #print(head(qcmap),1)
  #print(tail(qcmap),1)
  #print(qcmap)
  #print(xmap[,c(1:9,11:13)])

  point_size <- 0.75
  line_width <- 1

  offsets <- lapply(1:nrow(xmap),function(i) { #loop through xmap entries
    entry <- xmap[i,,drop=F]

    mol_start <- entry[,"QryStartPos"]
    map_start <- entry[,"RefStartPos"]
    orientation <- entry[,"Orientation"]
    alignment <- entry[,"Alignment"]
    alignment_split <- split_alignment(alignment)

    mol_id <- entry[,"QryContigID"]
    mol_len <- entry["QryLen"]

    ref_range <- abs(entry[,"RefStartPos"]-entry[,"RefEndPos"])
    qry_range <- abs(entry[,"QryStartPos"]-entry[,"QryEndPos"])
    stretch <- qry_range/ref_range

    #     num_mismatches <- sum(!qcmap[,"SiteID"]%in%alignment_split[,"qry"])
    #     fraction_mismatches <- num_mismatches/(length(qcmap[,"SiteID"])-1)
    #     if  (fraction_mismatches > 0.3) {
    #       return()
    #     }

    mol_positions <- qcmap[qcmap[,"CMapId"]==mol_id,"Position"]
    mol_ids <- qcmap[qcmap[,"CMapId"]==mol_id,"SiteID"]

    if (stretch_flag) {
      mol_positions <- mol_positions/stretch
      mol_start <- mol_start/stretch
    }

    if (orientation=="+") { #flip positions
      if (!is.null(xoffset)) {
        startX <- xoffset
      } else {
        startX <- map_start - mol_start
      }
      moleculeX <- startX + mol_positions
    } else {
      if (!is.null(xoffset)) {
        startX <- xoffset
      } else {
        startX <- map_start + mol_start
      }
      moleculeX <- startX - mol_positions
    }

    #draw molecule backbone
    if (orientation=="+") { #flip positions
      lines(c(moleculeX[1],moleculeX[length(moleculeX)]),c(i,i),col="black",lend=1,lwd=line_width)
    } else {
      lines(c(moleculeX[1],moleculeX[length(moleculeX)]),c(i,i),col="black",lend=1,lwd=line_width)
    }
    # lines(c(entry[,"RefStartPos"],entry[,"RefEndPos"]),c(i,i),col="red",lend=1,lwd=line_width) #highlight aligned portions

    lapply(1:length(moleculeX),function(j) {
      pos <- moleculeX[j]
      id <- mol_ids[j]

      #draw labels
      if (id%in%as.numeric(alignment_split[,"qry"])) {
        points(pos,i+yoffset,pch="|",col=colors[2],cex=point_size)
      } else {
        points(pos,i+yoffset,pch="|",col=colors[1],cex=point_size)
      }
    })

    #moleculeX <- moleculeX[(1:length(moleculeX))%in%alignment_split[,"qry"]] #filter for only labels that were aligned

    return(startX)
    return(moleculeX)

#     #add arrow
#     if (orientation=="+") {
#       arrows(x0=moleculeX[length(moleculeX)]-1000,x1=moleculeX[length(moleculeX)],y0=i,angle=20,length=0.12,col="red",lwd=2) #arrow at the end
#     } else {
#       arrows(x0=moleculeX[length(moleculeX)],x1=moleculeX[length(moleculeX)]-1000,y0=i,angle=20,length=0.12,col="red",lwd=2) #arrow at the beginning
#     }
  })

  return(unlist(offsets))
  return(unlist(molecule_pos))
}

draw_matchgroups <- function(xmap,rcmap,qcmap,color,xoffset=NULL) {
  #print("matchgroups")
  #print(xoffset)

  lapply(1:nrow(xmap),function(i) { #loop through xmap entries
    entry <- xmap[i,,drop=F]

    mol_start <- entry[,"QryStartPos"]
    map_start <- entry[,"RefStartPos"]
    orientation <- entry[,"Orientation"]
    alignment <- entry[,"Alignment"]
    alignment_split <- split_alignment(alignment)

    mol_id <- entry[,"QryContigID"]
    mol_len <- entry["QryLen"]

    ref_id <- entry[,"RefcontigID"]

    ref_range <- abs(entry[,"RefStartPos"]-entry[,"RefEndPos"])
    qry_range <- abs(entry[,"QryStartPos"]-entry[,"QryEndPos"])

#     if (ref_range==0||ref_range==0) {
#       stretch <- 1
#     } else {
#       stretch <- qry_range/ref_range
#     }

    mol_positions <- qcmap[qcmap[,"CMapId"]==mol_id,"Position"]

#     if (stretch_flag) {
#       mol_positions <- mol_positions/stretch
#       mol_start <- mol_start/stretch
#     }

    if (orientation=="+") {
      if (!is.null(xoffset)) {
        startX <- xoffset
        qcmap[,"Position"] <- startX + qcmap[,"Position"]
      }
    } else {
      if (!is.null(xoffset)) {
        startX <- xoffset
        qcmap[,"Position"] <- startX - qcmap[,"Position"]
      }
    }

    tmp_merge <- merge(alignment_split,qcmap,by.x="qry",by.y="SiteID",sort=F)
    #print(head(tmp_merge))
    merged <- merge(tmp_merge,rcmap,by.x="ref",by.y="SiteID",sort=F,suffixes=c(".q",".r"))
    #print(head(merged))

    lapply(1:nrow(merged),function(j) {
      match <- merged[j,,drop=F]
      #print(match)
      lines(c(match[,"Position.r"],match[,"Position.q"]),c(0.2,1),lwd=0.5,lty=1,col=color) #alignment lines
    })
  })
}

####################################

#util to renumber sites
renumber_sites <- function(cmap) {
  cmap[,"SiteID"] <- 1:nrow(cmap)

  return(cmap)
}

display_region <- function(indata,min_cushion) {
  xmap <- indata[["xmap"]]
  rcmap <- indata[["rcmap"]]
  qcmap <- indata[["qcmap"]]

  xrange <- range(xmap[,"RefStartPos"],xmap[,"RefEndPos"])
  pad <- min_cushion
  nrows <- nrow(xmap)
  region_str <- paste0("map",qcmap[1,"CMapId"])

  draw_reference(rcmap,xrange,nrows,pad,region_str)

  xmap <- split(xmap,xmap["LabelChannel"])
  #print(lapply(xmap,head))
  qcmap <- split(qcmap,qcmap["LabelChannel"])
  #qcmap <- lapply(qcmap,renumber_sites)
  #print(lapply(qcmap,head))
  rcmap <- split(rcmap,rcmap["LabelChannel"])

  xoffset <- draw_molecules(qcmap[["1"]],xmap[["1"]],stretch_flag=F,colors=c(rgb(1,0,0,0.1),rgb(1,0,0,0.9)),yoffset=0.015)
#   print(xoffset)
  draw_molecules(qcmap[["2"]],xmap[["2"]],stretch_flag=F,colors=c(rgb(0,1,0,0.1),rgb(0,1,0,0.9)),yoffset=-0.025,xoffset)

  draw_matchgroups(xmap[["1"]],rcmap[["1"]],qcmap[["1"]],rgb(1,0,0,0.75),xoffset)
  draw_matchgroups(xmap[["2"]],rcmap[["2"]],qcmap[["2"]],rgb(0,1,0,0.75),xoffset)

  return(1)
}

################################################################################
################################################################################
################################################################################
################################################################################

#load xmap/cmap
get_alignment <- function(folder_prefix,alignment_folder) {
  prefix <- paste0(alignment_folder,folder_prefix)

  xmap <- paste0(prefix,".xmap")
  xmap <- readxmap(xmap)

  qcmap <- paste0(prefix,"_q.cmap")
  qcmap <- readcmap(qcmap)

  rcmap <- paste0(prefix,"_r.cmap")
  rcmap <- readcmap(rcmap)

  out_list <- list("xmap"=xmap,
                   "qcmap"=qcmap,
                   "rcmap"=rcmap)

  return(out_list)
}

#select relevant entries
get_data <- function(id,alignments) {
  xmap <- alignments[["xmap"]]
  qcmap <- alignments[["qcmap"]]
  rcmap <- alignments[["rcmap"]]

  xmap <- xmap[xmap[,"QryContigID"]==id,]
  stopifnot(length(unique(xmap[,"RefcontigID"]))==1)

  qcmap <- qcmap[qcmap[,"CMapId"]==id,]

  rcmap <- rcmap[rcmap[,"CMapId"]==xmap[,"RefcontigID"],]

  out_list <- list("xmap"=xmap,
                   "qcmap"=qcmap,
                   "rcmap"=rcmap)

  return(out_list)
}

#wrapper
display_regions <- function(folder_prefix,
                            alignment_folder,
                            min_cushion=5e5) {

  alignments <- get_alignment(folder_prefix,alignment_folder)
  #save(list="alignments",file="/home/users/elam/20150911_2-col_actual_data/sop_and_ml/20151014_merged.rdata")

  #save(list="alignments",file="/home/users/elam/20150911_2-col_actual_data/sop_and_ml/test_2-color_align1.rdata")
  #load(file="/home/users/elam/20150911_2-col_actual_data/sop_and_ml/test_2-color_align1.rdata")

  qcmap_ids <- unique(alignments[["qcmap"]][,"CMapId"])#[c(15,29)] #looks good: 15,29

  lapply(qcmap_ids,function(id) {
    indata <- get_data(id,alignments)
    display_region(indata,min_cushion)
  })

  return(alignments)
}

################################################################################
################################################################################
################################################################################
################################################################################

# stop()
#
# alignment_folder <- "/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/" #for output
# alignment_folder <- "/home/users3/elam/20151028_2-col_V3_sandwich/" #20151029 EL
# alignment_folder <- "/home/users3/elam/20151030_2-col_V3_sandwich/cmap_output/" #20151031 EL - preparing input for Andy and Tiff
# alignment_folder <- "/home/users/elam/20151231_2-col_fragileSite_repaired/" #20160103
#
# #folder_prefix <- "test_2-color_align1"
# folder_prefix <- "20151014_merged" #20151014 EL - debugged cmap output
#
# folder_prefix <- "20151014_contig4_align" #20151017 EL - debugging still
# folder_prefix <- "20151014_merged_align" #20151017 EL - full set
# folder_prefix <- "20151014_merged1_align" #20151017 EL - full set (corrected)
#
# folder_prefix <- "sub1_2col" #20151029 EL
# folder_prefix <- "batches_indata_2col_noDelta"
#
# #folder_prefix <- "20151030_v3_sub30X_2col_cut_delta4"
# #folder_prefix <- "20151030_v3_sub30X_2col_cut"
# folder_prefix <- "20151030_v3_sub130X_2col_cut"
#
# folder_prefix <- "20151231_fragRep_sub30X_2col_cut" #20160103
#
# #pdf_output <- "/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/20151009_2-color_alignments.pdf"
# #pdf_output <- "/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/20151014_2-color_alignments.pdf"
# #pdf_output <- "/home/users3/elam/20150911_2-col_actual_data/sop_and_ml/20151017_2-color_alignments.pdf"
# #pdf_output <- "/home/users3/elam/20151028_2-col_V3_sandwich/sub1_2col_2-color_alignments1.pdf"
# #pdf_output <- "/home/users3/elam/20151028_2-col_V3_sandwich/batches_indata_2col_noDelta_2-color_alignments2.pdf"
#
# #pdf_output <- "/home/users3/elam/20151030_2-col_V3_sandwich/cmap_output/20151030_v3_sub30X_2col_cut_delta4_2-color_alignments.pdf"
# #pdf_output <- "/home/users3/elam/20151030_2-col_V3_sandwich/cmap_output/20151030_v3_sub30X_2col_cut_2-color_alignments.pdf"
# pdf_output <- "/home/users3/elam/20151030_2-col_V3_sandwich/cmap_output/20151030_v3_sub130X_2col_cut_2-color_alignments.pdf"
#
# pdf_output <- "/home/users3/elam/20151231_2-col_fragileSite_repaired/20160103_fragRep_sub30X_2col_cut_2-color_alignments.pdf"
#
# pdf(pdf_output,width=12,height=8)
# par(mfrow=c(2,1))
# try({
#   output <- display_regions(folder_prefix,
#                             alignment_folder,
#                             min_cushion=5e5)
# })
# dev.off()
