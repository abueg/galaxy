#20151027 EL - output final paths

# rm(list=ls())
#
# scripts_dir <- "/home/users/elam/rscripts/"
# source(paste(scripts_dir,"paste0.R",sep=""))
# source(paste0(scripts_dir,"readmaps.R")) #util functions to read _map formats
#
# scripts_dir <- "/home/users/elam/rscripts/2col_dev/current/"
# source(paste0(scripts_dir,"plotting.R"))
# source(paste0(scripts_dir,"graph.R"))

#get regression coefficients based on path and fits
get_regression_data <- function(inpath,g,pos_offsets) {

  node_type <- V(g)$type[get.diameter(g,directed=F)] #either 1 or 0
  node_selected <- V(g)$name[V(g)$type]
  node_selected <- node_selected[!is.na(node_selected)] #20150706

  el <- as.data.frame(get.edgelist(g))
  colnames(el) <- c("asm1","asm2")
  id_str <- paste(el[,"asm1"],el[,"asm2"],sep="_")

  stopifnot(sum(node_selected%in%el[,"asm2"])==(length(node_selected))) #type 1 nodes are from asm2

  path_df <- cbind("node"=inpath,node_type) #make dataframe with asm source information
  path_df[,"node_type"][path_df[,"node_type"]==1] <- "asm2"
  path_df[,"node_type"][path_df[,"node_type"]!="asm2"] <- "asm1"
  # print(path_df)

  index <- c(2:length(inpath))
  pair_com1 <- paste(inpath[index-1],inpath[index],sep="_")
  pair_com2 <- paste(inpath[index],inpath[index-1],sep="_")

  names(pair_com1) <- index
  names(pair_com2) <- index

  if (path_df[1,2]=="asm1") {
    pairs <- pair_com1
    pairs[(1:length(pairs))%%2==0] <- pair_com2[(1:length(pairs))%%2==0]
  }
  else {
    pairs <- pair_com2
    pairs[(1:length(pairs))%%2==0] <- pair_com1[(1:length(pairs))%%2==0]
  }
  # print(pairs)

  pairs <- as.data.frame(pairs,stringsAsFactors=F)
  colnames(pairs) <- c("id_str")

  pairs_merged <- merge(pairs,pos_offsets,by="id_str",sort=F)

  stopifnot(nrow(pairs_merged)==nrow(pairs))

  return(list("path_df"=path_df,
              "pairs_merged"=pairs_merged))

  return(pairs_merged)
}

get_positions1 <- function(outpath,pairs,data, node.types) { #get relative positions of the nodes in path
  if(66 %in% outpath){
      #browser()
  }
  data_nr <- data[!duplicated(data["id_str"]),]
  pairs_merged <- merge(pairs,data_nr,by="id_str",sort=F)

  columns_to_keep <- c("id_str","intercept","slope","RefcontigID.x","RefcontigID.y","RefLen.x","RefLen.y")
  pairs_merged <- pairs_merged[,names(pairs_merged)%in%columns_to_keep]

  #add orientation information
  pairs_merged[,"orientation"] <- as.numeric(pairs_merged[,"slope"])>-0.2
  pairs_merged[,"orientation"] <- as.numeric(pairs_merged[,"orientation"])
  pairs_merged[,"orientation"] <- pairs_merged[,"orientation"]*2-1 #turn from {0,1} into {-1,1}

  pairs_merged[,"rel_orientation"] <- cumprod(pairs_merged[,"orientation"])

  #get map ids
  id_str <- lapply(pairs_merged[,"id_str"],function(id_str) {
    return(unlist(strsplit(id_str,split='_')))
  })
  id_str <- do.call(rbind,id_str)
  colnames(id_str) <- c(1,2)

  pairs_merged <- cbind(pairs_merged,id_str)

  #initialize
  outdf <- data.frame(matrix(ncol=4,nrow=length(outpath)))
  colnames(outdf) <- c("node","rel_start","rel_end","orientation")

  #get first node
  node <- outpath[1]
  is.type.I <- (node.types[1,"node_type"] == "asm1") #check which asm the node is from
  xstart1 <- 1


  #if (node==pairs_merged[1,"1"]) {
  if(is.type.I){
    xend1 <- pairs_merged[1,"RefLen.x"]
  } else {
    xend1 <- pairs_merged[1,"RefLen.y"]
  }

  orientation <- 1
  outdf[1,] <- c(node,xstart1,xend1,orientation)


  #get the rest of the nodes
  for (i in 2:length(outpath)) {
    entry <- pairs_merged[i-1,,F]
    node <- outpath[i]
    is.type.I <- (node.types[i, "node_type"] == "asm1")
    #if(node == 66){
    #    browser()
    #}
    previous_node <- outdf[i-1,]

    if (is.type.I) {
      offset <- as.numeric(entry[,"intercept"]) #y=b when x=0
      map_len <- entry[,"RefLen.x"]
    } else {
      offset <- -1*as.numeric(entry[,"intercept"])/entry[,"orientation"] #x=-b/m when y=0
      map_len <- entry[,"RefLen.y"]
    }

    previous_orientation <- previous_node[,"orientation"]
    previous_xstart <- previous_node[,"rel_start"]

    #get map start
    if (previous_orientation==1) {
      xstart <- previous_xstart+offset
    } else {
      xstart <- previous_xstart-offset
    }

    #get map end
    if (entry[,"rel_orientation"]==1) {
      xend <- xstart+map_len
    } else {
      xend <- xstart-map_len
    }

    orientation <- entry[,"rel_orientation"]

    outdf[i,] <- c("node"=node,
                   "rel_start"=xstart,
                   "rel_end"=xend,
                   "orientation"=orientation)
  }

  return(outdf)
}

#find nodes that don't provide scaffolding information but are part of the path
find_embedded_nodes <- function(g) {
  longest_path <- get.diameter(g,directed=F) #get longest path [with nr ids]
  # print(longest_path)

  neighbors <- unique(unlist(neighborhood(g,1,nodes=longest_path))) #all 1-degree neighbors along longest path
  # print(neighbors)

  missed_nodes <- V(g)[V(g)%in%neighbors & !V(g)%in%longest_path] #nodes in neighbors but not longest path
  # print(missed_nodes)

  degree_missed_nodes <- degree(g,missed_nodes) #check degree for each missed node
  # print(degree_missed_nodes)

  missed_nodes_filtered <- missed_nodes[degree_missed_nodes==1] #filter nodes with >1 degree; they likely lead to other paths
  # print(missed_nodes_filtered)

  return(missed_nodes_filtered)
}

get_node_offsets <- function(nodes,g,pos_offsets) { #getting offsets for arbitrary nodes

  node_type <- V(g)$type[nodes] #either 1 or 0

  el <- as.data.frame(get.edgelist(g))
  colnames(el) <- c("asm1","asm2")
  id_str <- paste(el[,"asm1"],el[,"asm2"],sep="_")

  node_df <- cbind("node"=V(g)$name[nodes],node_type) #make dataframe with asm source information
  node_df[,"node_type"][node_df[,"node_type"]==1] <- "asm2"
  node_df[,"node_type"][node_df[,"node_type"]!="asm2"] <- "asm1"
  # print(node_df)

  id_str <- lapply(pos_offsets[,"id_str"],function(id_str) {
    id_str_split <- unlist(strsplit(id_str,"_"))
  })
  id_str <- do.call(rbind,id_str)
  colnames(id_str) <- c("asm1","asm2")

  pos_offsets <- cbind(pos_offsets,id_str)
  # print(pos_offsets)

  node_offset_list <- lapply(1:nrow(node_df),function(i) {
    entry <- node_df[i,,drop=F]
    # print(entry)
    pos_offsets_selected <- pos_offsets[pos_offsets[,entry[,"node_type"]]==entry[,"node"],]

    return(pos_offsets_selected)
  })
  node_offset_list <- do.call(rbind,node_offset_list)

  node_offset_list <- node_offset_list[!duplicated(node_offset_list[,"id_str"]),]
  # print(node_offset_list)

  return(node_offset_list)
}

get_missed_nodes_lengths <- function(pos_offsets,indata) { #getting lengths for embedded nodes
  reflens <- lapply(1:nrow(pos_offsets),function(i) {
    entry <- pos_offsets[i,]
    id_str <- entry[,"id_str"]

    merged <- merge(indata,data.frame("id_str"=id_str),by="id_str")
    merged <- head(merged,1)

    reflen.x <- merged[,"RefLen.x"]
    reflen.y <- merged[,"RefLen.y"]

    #     print(reflen.x)
    #     print(reflen.y)

    out_list <- list("reflen.x"=reflen.x,
                     "reflen.y"=reflen.y)

    out_df <- do.call(cbind,out_list)

    return(out_df)
  })
  reflens <- do.call(rbind,reflens)

  pos_offsets <- cbind(pos_offsets,reflens)
  # print(pos_offsets)

  return(pos_offsets)
}

add_embedded_nodes1 <- function(path_df,pos_offsets) {

  path_df[,"node_str"] <- paste0(path_df[,"node"],"_",path_df[,"node_type"])

  pos_offsets[,"orientation"] <- as.numeric(pos_offsets[,"slope"])>0
  pos_offsets[,"orientation"] <- as.numeric(pos_offsets[,"orientation"])
  pos_offsets[,"orientation"] <- pos_offsets[,"orientation"]*2-1 #turn from {0,1} into {-1,1}

  pos_offsets[,"node_str1"] <- paste0(pos_offsets[,"asm1"],"_","asm1")
  pos_offsets[,"node_str2"] <- paste0(pos_offsets[,"asm2"],"_","asm2")

  new_rows <- lapply(1:nrow(pos_offsets),function(i) {
    entry <- pos_offsets[i,]
    if (entry[,"node_str1"]%in%path_df[,"node_str"]) { #flip depending on if the node is from asm1 or asm2
      offset <- -1*as.numeric(entry[,"intercept"])/entry[,"orientation"] #x=-b/m when y=0
      merged <- merge(entry,path_df,by.x="node_str1",by.y="node_str")
      node <- merged[,"asm2"]
      node_type <- "asm2"
      map_len <- merged[,"reflen.y"]
    } else {
      offset <- as.numeric(entry[,"intercept"]) #y=b when x=0
      merged <- merge(entry,path_df,by.x="node_str2",by.y="node_str")
      node <- merged[,"asm1"]
      node_type <- "asm1"
      map_len <- merged[,"reflen.x"]
    }

    if (nrow(merged)==0) {
      return(NULL)
    }

    anchor_orientation <- merged[,"orientation.y"]
    anchor_xstart <- merged[,"rel_start"]

    current_orientation <- merged[,"orientation.x"]
    rel_orientation <- anchor_orientation*current_orientation

    #get map start
    if (anchor_orientation==1) {
      xstart <- anchor_xstart+offset
    } else {
      xstart <- anchor_xstart-offset
    }

    #get map end
    if (rel_orientation==1) {
      xend <- xstart+map_len
    } else {
      xend <- xstart-map_len
    }

    out_df <- data.frame("node"=node,
                         "rel_start"=xstart,
                         "rel_end"=xend,
                         "orientation"=rel_orientation,
                         "node_type"=node_type)

    return(out_df)
  })
  new_rows <- do.call(rbind,new_rows)

  return(new_rows)
}

add_orientation <- function(path_coordinates) { #add orientation information
  path_coordinates[,"orientation"] <- path_coordinates[,"rel_start"]<path_coordinates[,"rel_end"]
  path_coordinates[,"orientation"] <- path_coordinates[,"orientation"]*2-1

  return(path_coordinates)
}

get_path1 <- function(g,pos_offsets,indata,alignment_to_ref,plot_flag=F) {
  #print(g)
  outpath <- V(g)$name[get.diameter(g,directed=F)] #longest path
  pairs <- get_regression_data(outpath,g,pos_offsets)

#   return(pairs[["path_df"]]) #used for testing node types

  path_coordinates <- get_positions1(outpath,pairs[["pairs_merged"]],indata, pairs$path_df)

#   print(lapply(pairs[["path_df"]],class))
#   print(lapply(path_coordinates,class))

  path_coordinates1 <- cbind(path_coordinates,"node_type"=pairs[["path_df"]][,"node_type"]) #annotate asm type

#   return(path_coordinates1) #used for testing node types

#   plot_schematic(path_coordinates1,T,T)


  xmap1 <- alignment_to_ref[[1]] #alignment of assembly 1 maps to reference
  xmap2 <- alignment_to_ref[[2]] #alignment of assembly 2 maps to reference
  if(!is.null(xmap1) && !is.null(xmap2)){
    correlation <- check_predictions(path_coordinates1,xmap1,xmap2,plot_flag,title_str='')
      #if(!is.na(correlation)){
      #    correlation <- correlation$correlation
      #}
  }
  stopifnot(path_coordinates[,"node"]==pairs[["path_df"]][,"node"])

  missed_nodes <- find_embedded_nodes(g) #get missing mnodes

  if (length(missed_nodes)>0) {
    missed_nodes_offsets <- get_node_offsets(missed_nodes,g,pos_offsets) #get offsets for missing nodes
    #browser()
    missed_nodes_offsets1 <- get_missed_nodes_lengths(missed_nodes_offsets,indata)
    path_coordinates2 <- add_embedded_nodes1(path_coordinates1,missed_nodes_offsets1) #positions for embedded nodes (20151025 EL - new function)
    path_coordinates3 <- rbind(path_coordinates1,path_coordinates2)
  } else{
    path_coordinates1 <- add_orientation(path_coordinates) #adding orientation information
    path_coordinates1 <- cbind(path_coordinates1,"node_type"=pairs[["path_df"]][,"node_type"]) #annotate asm type
    path_coordinates3 <- path_coordinates1
  }

  if (plot_flag) {
    plot_schematic(path_coordinates1,T,T)
    plot_schematic(path_coordinates3,T,T)
  }

  graph_summary <- summarize_graph(g,outpath,missed_nodes)

  xmap1 <- alignment_to_ref[[1]] #alignment of assembly 1 maps to reference
  xmap2 <- alignment_to_ref[[2]] #alignment of assembly 2 maps to reference

  if (plot_flag) {
    correlation <- check_predictions(path_coordinates3,xmap1,xmap2,plot_flag=T,title_str='')
    #if(!is.na(correlation)){
    #   path_coordinates3 <- merge(path_coordinates3, correlation$residuals, by=c("node", "node_type"))
    #}
  }



  outlist <- list("summary"=graph_summary,
                  "path_coordinates"=path_coordinates3)

  return(outlist)
  return(path_coordinates3)
}

# #main
# prefix <- "/home/users/elam/20150911_2-col_actual_data/sop_and_ml_test/Human_willsample_twocolor_combined"
# load(file=paste0(prefix,"_indata_filtered.RData")) #load indata_filtered
#
# load(file=paste0(prefix,"_pos_offsets.RData")) #load pos_offsets
# load(file=paste0(prefix,"_clusters_g.RData")) #load clusters_g
#
# bsssi_xmap <- "/home/users3/elam/data/NA12891_bsssi_sop_50X_2/output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap" #new 50X assembly
# bspqi_xmap <- "/home/users3/elam/data/NA12891_bspqi_ml_50X_4/output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap" #new 50X assembly (ML)
#
# alignment_to_ref <- list("asm1"=bsssi_xmap,
#                          "asm2"=bspqi_xmap)
#
# align_dir <- "/home/users/elam/20150911_2-col_actual_data/sop_and_ml_test/"
# prefix <- "Human_willsample_twocolor_combined"
#
# set.seed(1)
# par(mfrow=c(2,2))
#
# try({
#   cluster_sample <- sample(1:length(clusters_g),length(clusters_g)*1.0)
#   clusters_g <- clusters_g[cluster_sample]
#
#   paths1 <- lapply(clusters_g[2],get_path1,pos_offsets,indata_filtered,alignment_to_ref,plot_flag=T)
# })
