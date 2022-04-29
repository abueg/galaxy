#This script is adopted from Ernest wrapper script
#It loads neccessary scripts and dependency for sandwich

#tested on process2
#R version 3.2.2

#clean workspace first (optional)
#rm(list=ls())    #a little danger of load from middle of codes

options(scipen=999)

#global

library("igraph") #for graph-related functions
library("intervals") #for comparing intervals
library("MASS") #for rlm
library("parallel") #running things in parallel
library("data.table") #quicker dataframe operations

#general utility scripts
#loading depending script for sandwich only need if running outside of package
initSandwich <-function(){
  scripts_dir <- "/home/users/elam/rscripts/"
  source(paste(scripts_dir,"paste0.R",sep=""))
  source(paste0(scripts_dir,"readmaps.R")) #util functions to read _map formats
  source(paste0(scripts_dir,"overlap.R")) #test interval overlaps
  source(paste0(scripts_dir,"readcmap_w_headers.R")) #updated readcmap function to also store header information
  source(paste0(scripts_dir,"calc_n50.R")) #related functions to calculate N50 of cmap
  source(paste0(scripts_dir,"cigar1.R")) #parse cigar string
  source(paste0(scripts_dir,"igraph_utils.R")) #igraph utility functions


  #2-color sandwich specific scripts
  scripts_dir <- "/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/sandwichScript/"
  source(paste0(scripts_dir,"utils.R"))
  source(paste0(scripts_dir,"do_analysis.R"))
  source(paste0(scripts_dir,"get_input_alignment.R"))
  source(paste0(scripts_dir,"do_validation.R"))
  source(paste0(scripts_dir,"analyze_alignment.R"))
  source(paste0(scripts_dir,"analyze_alignment_simple.R"))
  source(paste0(scripts_dir,"make_plots.R"))
  source(paste0(scripts_dir,"make_graph.R"))
  source(paste0(scripts_dir,"output_path.R"))
  source(paste0(scripts_dir,"summarizes_scaffolds.R"))
  source(paste0(scripts_dir,"run_ra.R"))
  source(paste0(scripts_dir,"output_cmap.R"))
  source(paste0(scripts_dir,"visualize_alignment.R"))
  source(paste0(scripts_dir,"readmaps.R"))
  source(paste0(scripts_dir, "run_sandwich.R"))
  source(paste0(scripts_dir, "ConflictResolveTwo.R"))
  source(paste0(scripts_dir, "2-Color_AlignFinal.R"))
  source(paste0(scripts_dir, "readmaps.R"))
  source(paste0(scripts_dir, "Utility.R"))
  source(paste0(scripts_dir, "exportMaps.r"))
}
