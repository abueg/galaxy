#this script runs the full sandwich assemblies

run.sandwich <- function(align1, align2, asm1, asm2, output, bspqi_xmap=NULL, bsssi_xmap=NULL, score.T=13){
    #clean workspace first (optional)
    options(scipen=999)
    #main
    #where the RData and other output files will be output
    curr.wd <- getwd()
    setwd(output)

    #input parameters
    #prefix input (based on the bnx prefixes)
    prefixes <- "test"

    #summarize molecule-to-assembly mapping rates
    align_dir <- paste0(getwd(),"/")
    repeat_flag <- T #whether to refdo this analysis
    save_flag <- T #whether to save df to RData file
    cores <- 32 #max should be 32
    rates <- get_alignment_rates(prefixes,align_dir,repeat_flag=repeat_flag,save_flag=save_flag,cores=cores)

    #main analysis callls
    final_output <- lapply(prefixes,function(prefix) { #process each bnx of interest based on prefix
        #variables
        input <- list("align_dir"=align_dir,
                      "prefix"=prefix,
                      "load_previous_flag"=F)

        input_type <- "prefix"

        output_dir <- "" #generally not used

        #get assembly data
        references <- get_assemblies(asm1,asm2)

        alignment_to_ref <- list("asm1"=bspqi_xmap,
                                 "asm2"=bsssi_xmap)
        #main call
        selected_chr <- NULL
        check_ans_flag <- F

        paths <- NULL
        path_summary <- NULL
        #browser()
        #try({
            paths <- run_analysis(input,input_type,output_dir,alignment_to_ref,selected_chr,check_ans_flag, min.score.T = score.T) #get paths
            path_summary <- make_summary1(paths,asm1,asm2,alignment_to_ref) #summarize paths
            save(list="path_summary", file=paste0(prefix, "_path_summmary.RData"))
            cmap_list <- lapply(1:length(paths),function(i) { #prepare cmap df
                path <- paths[[i]]
                path_coordinates <- path[["path_coordinates"]]
                cmap_df <- make_cmap(path_coordinates,references,plot_flag=F,test_flag=F)
                                        #exit()
                return(cmap_df)
              })

            save(list="cmap_list",file=paste0(prefix,"_cmap_list.RData"))

           #renumber map ids sequentially (otherwise, all maps had an id of 1)
            mode <- "sequential"
            parallel_flag <- T
            cores <- 2

            cmap_list1 <- renumber_mapid(cmap_list,mode,parallel_flag,cores)
            cmap_df <- do.call(rbind,cmap_list1)

            #get the non-scaffold contigs
            #browser()
            non_scaffolds_list <- get_nonscaffolds(references[[1]], references[[2]], path_summary$non_scaffolds, cmap_df)

           #output a 2-color cmap file
            input <- list("filename"=paste0(prefix,"_2col.cmap"),
                          "cmap"=cmap_df)

            mode <- "single"
            renumber_flag <- F

            id_list <- write_cmap(input,mode,renumber_flag)

            #output the two-color map as two-single color map
            split_cmap_list <- splitTwocolorCmap(cmap_list1)
            cmap_df_col1 <- do.call(rbind, split_cmap_list$cmap_col1)
            cmap_df_col1 <- rbind(cmap_df_col1, non_scaffolds_list[[1]])
            cmap_df_col2 <- do.call(rbind, split_cmap_list$cmap_col2)
            cmap_df_col2 <- rbind(cmap_df_col2, non_scaffolds_list[[2]])

            output_file1 <- list("filename"=paste0(prefix,"_2col_c1.cmap"),
                                 "cmap"=cmap_df_col1)

            output_file2 <- list("filename"=paste0(prefix,"_2col_c2.cmap"),
                                 "cmap"=cmap_df_col2)
            id_list1 <- write_cmap(output_file1, mode, renumber_flag)
            id_list2 <- write_cmap(output_file2, mode, renumber_flag)

            ra_output <- NULL

            #trimming cmap
            cmap <- id_list
            outfile <- paste0(prefix,"_2col_cut.cmap")
            system(paste("cut",cmap,"-f 1-9 > ",outfile))

            outlist <- list("paths"=paths,
                            "path_summary"=path_summary,
                            "cmap_list"=cmap_list,
                            "output_cmap"=id_list,
                            "ra_output"=ra_output)

            save(list="outlist",file=paste0(prefix,"_outlist.RData"))
            setwd(curr.wd)
            return(outlist)
        #})
        setwd(curr.wd)
        return()

        #post-analysis
        summary <- process_output(final_output,rates)

        outfile <- NULL
        plot_flag <- T
        plot_summary(summary,outfile,plot_flag)

    })
}

get_nonscaffolds <- function(asm1, asm2, non_scaffolds, scaffold_cmap){
    asm1 <- readcmap(asm1)
    asm2 <- readcmap(asm2)
    trimms <- trim_cmap(asm1, asm2)
    asm1 <- trimms[[1]]
    asm2 <- trimms[[2]]
    non_scaffold_ids1 <- non_scaffolds$CMapId[non_scaffolds$node_type =='asm1']
    non_scaffold_ids2 <- non_scaffolds$CMapId[non_scaffolds$node_type =='asm2']

    non_scaffold_cmap1 <- asm1[asm1$CMapId %in% non_scaffold_ids1,]
    non_scaffold_cmap2 <- asm2[asm2$CMapId %in% non_scaffold_ids2,]

    id_shift <- 100000
    #print(nrow(non_scaffold_cmap1))
    if(nrow(non_scaffold_cmap1) > 0){
        non_scaffold_cmap1$CMapId <- paste(id_shift, non_scaffold_cmap1$CMapId, sep="")
    }

    id_shift2 <- id_shift * 2
    if(nrow(non_scaffold_cmap2) > 0){
        non_scaffold_cmap2$CMapId <- paste(id_shift2, non_scaffold_cmap2$CMapId, sep="")
    }
    return(list(non_scaffold_cmap1, non_scaffold_cmap2))
}
