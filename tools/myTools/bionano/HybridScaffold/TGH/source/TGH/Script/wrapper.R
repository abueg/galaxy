#20151012 EL - wrapper for 2-color analysis

#clean workspace first (optional)
rm(list=ls())

options(scipen=999)

#global
ra <- "/home/users/elam/RefAligner/ra_151016/bin/RefAligner"

#load libraries
library("igraph") #for graph-related functions
library("intervals") #for comparing intervals
library("MASS") #for rlm
library("parallel") #running things in parallel
library("data.table") #quicker dataframe operations

#general utility scripts
scripts_dir <- "/home/users/elam/rscripts/"
source(paste(scripts_dir,"paste0.R",sep=""))
source(paste0(scripts_dir,"readmaps.R")) #util functions to read _map formats
source(paste0(scripts_dir,"overlap.R"))
source(paste0(scripts_dir,"readcmap_w_headers.R"))
source(paste0(scripts_dir,"calc_n50.R"))
source(paste0(scripts_dir,"cigar1.R"))
source(paste0(scripts_dir,"igraph_utils.R"))

#2-color specific scripts
scripts_dir <- "/home/users/elam/rscripts/2col_dev/current/"
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

#main
#where the RData files will be output
setwd("/home/users3/elam/20150911_2-col_actual_data/sop_and_ml_test1") #20151011 - TESTING

#input parameters
#prefix input (based on the bnx prefixes)
prefixes <- "Human_willsample_twocolor_combined" #20151011 - TESTING

#summarize mapping rates
align_dir <- paste0(getwd(),"/")
repeat_flag <- T
save_flag <- F
cores <- 32 #max should be 32
rates <- get_alignment_rates(prefixes,align_dir,repeat_flag=repeat_flag,save_flag=save_flag,cores=cores)

#main analysis call
final_output <- lapply(prefixes,function(prefix) {
  #variables
  input <- list("align_dir"=align_dir,
                "prefix"=prefix,
                "load_previous_flag"=T)

  input_type <- "prefix"

  output_dir <- "" #generally not used

  #bspqi_xmap <- "/home/users3/elam/data/NA12891_bspqi_sop_50X_1/output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap" #new 50X assembly (SOP)
  bspqi_xmap <- "/home/users3/elam/data/NA12891_bspqi_ml_50X_4/output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap" #new 50X assembly (ML)

  bsssi_xmap <- "/home/users3/elam/data/NA12891_bsssi_sop_50X_2/output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap" #new 50X assembly

  alignment_to_ref <- list("asm1"=bsssi_xmap,
                           "asm2"=bspqi_xmap)

  asm1 <- paste0(align_dir,prefix,"_1_bppAdj.cmap")
  asm2 <- paste0(align_dir,prefix,"_2_bppAdj.cmap")

  #get assembly data
  references <- get_assemblies(asm1,asm2)

  #main call
  selected_chr <- NULL
  check_ans_flag <- F

  paths <- NULL
  path_summary <- NULL

  try({
    paths <- run_analysis(input,input_type,output_dir,alignment_to_ref,selected_chr,check_ans_flag) #get paths
    path_summary <- make_summary1(paths,asm1,asm2,alignment_to_ref) #summarize paths

    cmap_list <- lapply(1:length(paths),function(i) { #prepare cmap df
      path <- paths[[i]]
      path_coordinates <- path[["path_coordinates"]]
      cmap_df <- make_cmap(path_coordinates,references,plot_flag=F,test_flag=F)

      return(cmap_df)
    })

    save(list="cmap_list",file=paste0(prefix,"_cmap_list.RData"))

    #renumber map ids sequentially (otherwise, all maps had an id of 1)
    mode <- "sequential"
    parallel_flag <- T
    cores <- 2

    cmap_list1 <- renumber_mapid(cmap_list,mode,parallel_flag,cores)
    cmap_df <- do.call(rbind,cmap_list1)

    #output a 2-color cmap file
    input <- list("filename"=paste0(prefix,"_2col.cmap"),
                  "cmap"=cmap_df)

    mode <- "single"
    renumber_flag <- F

    id_list <- write_cmap(input,mode,renumber_flag)

    #run 2-color alignment of the scaffolds to reference
    cmap <- id_list
    ref <- "/home/users/csecol/genomes/human/hg19/hg19_BSSSI_BSPQI_0kb_0labels.cmap" #bsssi/bspqi reference
    prefix <- paste0(prefix,"_2col")
    run_flag <- T

    ra_output <- run_ra_2col(cmap,ref,prefix,run_flag)

    #visualize alignment
    alignment_folder <- align_dir #for output
    folder_prefix <- prefix

    pdf_output <- paste0(alignment_folder,folder_prefix,".pdf")
    pdf(pdf_output,width=12,height=8)
    par(mfrow=c(2,1))
    try({
      output <- display_regions(folder_prefix,
                                alignment_folder,
                                min_cushion=5e5)
    })
    dev.off()
  })

  outlist <- list("paths"=paths,
                  "path_summary"=path_summary,
                  "cmap_list"=cmap_list,
                  "output_cmap"=id_list,
                  "ra_output"=ra_output)

  save(list="outlist",file=paste0(prefix,"_outlist.RData"))

  return(outlist)
})

#post-analysis
summary <- process_output(final_output,rates)

outfile <- NULL
plot_flag <- T
plot_summary(summary,outfile,plot_flag)
