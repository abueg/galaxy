#05212015 EL - some igraph utilities

#check if there are cycles in the graph
find_cycles <- function(graph, k) { #functional but not particularly useful
  ring <- graph.ring(k) #only makes cycles at k>2
  graph.get.subisomorphisms.vf2(graph, ring)
}

delete_isolates <- function(g,mode='all') { #delete isolates in graph; will renumber vertices
  isolates <- which(degree(g,mode=mode)==0)
  g_pruned <- delete.vertices(g,isolates) #ids not preserved

  return(g_pruned)
}

graph_types <- c("layout.auto", #possible algorithms to plot graph
                 #"layout.random",
                 #"layout.circle",
                 #"layout.sphere",
                 "layout.fruchterman.reingold",
                 "layout.kamada.kawai",
                 "layout.spring",
                 "layout.reingold.tilford",
                 "layout.fruchterman.reingold.grid",
                 "layout.lgl",
                 "layout.graphopt",
                 "layout.svd")

plot_graph <- function(graph_type="layout.auto",g) { #plot graph with "reasonable" parameters
  plot.igraph(g,
              layout=get(graph_type),
              vertex.size=2,
              #vertex.color=rgb(0,0,0,0.3),
              vertex.frame.color=NULL,
              vertex.label.cex=0.2,
              vertex.label.color="black",
              edge.color=rgb(0,0,0,0.3),
              edge.width=0.5,
              edge.arrow.size=0.2,
              edge.arrow.width=0.75,
              edge.lty=1,
              main="",
              sub=graph_type)
}
