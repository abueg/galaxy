#20151026 EL - graph-related functions

summarize_graph <- function(g,outpath,missed_nodes) {
  n_all_nodes <- sum(!is.na(V(g)$name))
  n_scaffolded_nodes <- length(outpath)
  n_embedded_nodes <- length(missed_nodes)
  n_missed_nodes <- n_all_nodes-n_scaffolded_nodes-n_embedded_nodes

  out_df <- data.frame("n_all_nodes"=n_all_nodes,
                       "n_scaffolded_nodes"=n_scaffolded_nodes,
                       "n_embedded_nodes"=n_embedded_nodes,
                       "n_missed_nodes"=n_missed_nodes)

  return(out_df)
}

plot_bipartite <- function(g,mode="single") { #wrapper for plotting

  g_pruned <- delete_isolates(g) #delete isolates before plotting (vertices deleted earlier)

  try({
    if (mode=="single") {
      plot_graph(g=g_pruned) #make one plot
    } else {
      lapply(graph_types,plot_graph,g_pruned) #make multiple plots
    }
  })

  return(g_pruned)
}

get_clusters <- function(data,thresh=0) { #get connect components

  data <- data[order(data[,"num_edges"],decreasing=TRUE),,drop=FALSE]

  L<-list()

  while(nrow(data)>0) {
    enz1_id <- data[,"asm1_renumbered"]
    enz2_id <- data[,"asm2_renumbered"]

    current_seed <- c(enz1_id[1],enz2_id[1])

    run_flag <- TRUE
    while (run_flag) {
      F <- enz1_id %in% current_seed | enz2_id %in% current_seed
      next_seed <- unique(c(current_seed,enz1_id[F],enz2_id[F]))

      if (!(length(current_seed)<length(next_seed))) {
        run_flag <- FALSE
      }

      current_seed <- next_seed
    }

    if (sum(data[F,"num_edges"])<thresh) { #filter by num_edges
      data <- data[!F,,drop=FALSE]
      next
    }

    #     cat("Component",length(L)+1," with", sum(F), "reduced edges,", sum(data[F,"num_edges"]), "original edges,", length(current_seed), " molecules.\n")

    L[[length(L)+1]] <- data[F,,drop=FALSE]
    data <- data[!F,,drop=FALSE]
  }

  print(paste0("found ",length(L)," clusters"))

  return(L)
}

get_bipartite <- function(data) { #add type and color information
  enz1_max <- max(as.numeric(unique(data[,"asm1_renumbered"])))

  edges <- data[,c("asm1_renumbered","asm2_renumbered")]
  g <- graph.edgelist(as.matrix(edges))
  # plot_graph(g=g)

  V(g)$type <- V(g) > enz1_max

  names1 <- sort(unique(as.numeric(data[,"asm1"]))) #get unique mapIDs
  names2 <- sort(unique(as.numeric(data[,"asm2"]))) #get unique mapIDs
  names <- c(names1,names2)

  non_isolates <- which(degree(g,mode="all")!=0)
  # print(non_isolates)

  stopifnot(length(names)==length(non_isolates))
  V(g)$name[non_isolates] <- names

  color1 <- rgb(1,0,0,0.2)
  color2 <- rgb(0,0,1,0.2)

  V(g)$color<-ifelse(V(g)$type==TRUE,color1,color2)
  # plot_graph(g=g)

  return(g)
}
