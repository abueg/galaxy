#This script serves as the main controller script to run the new V3 hybridscaffold pipeline for V3
#scripts_dir <- "/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/sandwichScript/"
#hybrid.script.dir <- "/home/users//jwang/workspace/codebase/HybridScaffold/trunk/scripts/"
#initialize dependent pacakge and dependent script
#source(paste0(scripts_dir,"initialize.R"))

#a wrapper for running system command line
run.system.cmd <- function(cmd.str, run=T, runlog="", wait=T, print=T){
  print(paste0("Running command: ", cmd.str))
  if(run){
    #Piping stderr in stdout as well, so we can capture in output from R command
    out= system(paste0(cmd.str, ' 2>&1'), wait=T, intern=TRUE)
    if(print){
      print(out)
    }
    if(!is.null(attr(out, 'status'))){
        stop(paste0("Command ", cmd.str, " finished with errors.\n", paste(out, collapse = "\n")))
    }
    return(out)
  }
  return("")
}

#a wrapper to use mclapply to run command with error handling
mclapply.with.error <- function(arg.list, fun, mc.cores){
  ret <- mclapply(arg.list, fun, mc.cores=mc.cores)
  error.header <- paste0("Parallel command finished with error: ",
                         "Function: ",
                         paste(as.character(body(fun)), collapse="\n"))

  checks <- lapply(seq(1:length(ret)), function(i){
                            r <- ret[[i]]
                            if(class(r) == "try-error"){
                              error.msg <- paste0("Argument-list index: ", i, "\n",
                                                  "Argument-value: ", arg.list[[i]], "\n",
                                                  "Err msg: ", paste(r, collapse = ""), "\n")
                              return(error.msg)
                            }else{
                              return(NULL)
                            }
                        })
  if(!all(sapply(checks, is.null))){
      stop(paste0(error.header, paste(unlist(checks), collapse = ""), "\n"))
  }
  return(ret)
}

#running single-hybrid pipeline with error capture from its log
run.hybrid.with.error.handle <- function(cmd.list, output.list, run.flag){
  ret <- tryCatch({mclapply.with.error(cmd.list,
                                  function(cmd){run.system.cmd(cmd, run.flag)}, mc.cores=length(cmd.list))},
             error = function(err){
                    #grep the last status update from status.txt file
                    status.list <- lapply(output.list, function(out.dir){
                           try({lines <- readLines(file(paste0(out.dir, "/status.txt"), "rt"))
                           return(tail(lines,1))})
                    })
                    #print("status: ")
                    #print(status.list)
                    cat(err$message)
                    is.err <- grepl("ERROR:", status.list)
                    if(sum(is.err) > 0){
                      stop(paste0("Error running hybridScaffold.pl: ",
                                  paste(status.list[is.err], collapse = " | ")))
                    }else{
                      stop(err$message)
                    }
             })
}

#run the hybridscaffold up to the conflict resolution step
#this step obtain the coordinate translation information for translating
#coordinate in one bng asm to a second bng assembly using ngs alignment as the bridge
run.NGS2BNG.align <- function(run.flag, configs){
  hybridOut1 <- paste0(configs$OutputDir, "/", configs$Enzyme1)
  hybridOut2 <- paste0(configs$OutputDir, "/", configs$Enzyme2)
  runHybrid1 <- paste0("perl ", configs$HybridScaffoldPath, " -r ", configs$RefAlignerPath,
                       " -n ", configs$NGSPath, " -b ", configs$BNGPath1, " -c ", configs$Hybrid.xml1, " -o ", hybridOut1, " -B 2 -N 2 -S")
  runHybrid2 <- paste0("perl ", configs$HybridScaffoldPath, " -r ", configs$RefAlignerPath,
                       " -n ", configs$NGSPath, " -b ", configs$BNGPath2, " -c ", configs$Hybrid.xml2, " -o ", hybridOut2, " -B 2 -N 2 -S")
  #ret <- lapply(list(runHybrid1, runHybrid2),
  #              function(cmd){run.system.cmd(cmd, run.flag)})
  ret <- run.hybrid.with.error.handle(list(runHybrid1, runHybrid2),
                                      list(hybridOut1, hybridOut2), run.flag)

  return(list(hybridOut1, hybridOut2))
}


#run sandwich pipeline up to the point where pos_offset is generated
#we need this to perform the conflict resolution on the two bng assemblies
run.pre.sandwich <- function(asm1, asm2, alignToRef1, alignToRef2, sandwich1.dir, min.score=13){
  options(scipen=999)
  #2-color sandwich specific scripts
  #scripts_dir <- "/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/sandwichScript/"
  #initialize dependent pacakge and dependent script
  #source(paste0(scripts_dir,"initialize.R"))
  #main
  #where the RData and other output files will be output
  curr.wd <- getwd()
  setwd(sandwich1.dir)

  #input parameters
  #prefix input (based on the bnx prefixes)
  prefixes <- "test"

  #summarize molecule-to-assembly mapping rates
  align_dir <- paste0(getwd(),"/")
  repeat_flag <- T #whether to redo this analysis
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
  print(paste0("asm1: " , asm1))
  print(paste0("asm2: " , asm2))

  references <- get_assemblies(asm1,asm2)

  alignment_to_ref <- list("asm1"=alignToRef1,
                           "asm2"=alignToRef2)
  #browser()
  #main call
  selected_chr <- NULL
  check_ans_flag <- F

  paths <- NULL
  path_summary <- NULL
  #print("Hello")

  paths <- run_pre_analysis(input,input_type,output_dir,alignment_to_ref,selected_chr,check_ans_flag, twocolor.as.ref=F, min.score.T = min.score ) #get paths

  setwd(curr.wd)
  return(paths)
})

}

#detect and resolve bng conflicts from two single-enzyme bng assemblies
resolve.merged.conflicts <- function(hybridOut1, hybridOut2, run.flag=T, configs){
  sandwich1.out <- paste0(configs$OutputDir, "/Sandwich1")
  make.sandwich.out <- paste0("mkdir ", sandwich1.out)

  #if maunal cut file is presented we do not perform automatic merged conflict resolution
  if(!is.na(configs$ManualCut1)){
    if(file.exists(configs$ManualCut1) && file.exists(configs$ManualCut2)){
      return(list(configs$ManualCut1, configs$ManualCut2))
    }else{
      stop(paste0("Manual conflict cut file do not exists: \n", configs$ManualCut1, "\n", configs$ManualCut2))
    }
  }

  combined.cut.conflicts <- paste0(configs$OutputDir, "/CombinedConflictsCut")
  combined.conflicts1 <- paste0(combined.cut.conflicts, "/Combined_conflicts_cut_status_", configs$Enzyme1, ".txt")
  combined.conflicts2 <- paste0(combined.cut.conflicts, "/Combined_conflicts_cut_status_", configs$Enzyme2, ".txt")
  configs$ManualCut1 <- combined.conflicts1
  configs$ManualCut2 <- combined.conflicts2
  if(run.flag){
      #setup output directory for sandiwch1
      run.system.cmd(make.sandwich.out, run.flag)
      #copy input alignments to approprirate directory
      cp.align1 <- paste0("cp ", hybridOut1, "/align1/align1.xmap ", sandwich1.out, "/test_1_fullAln.xmap")
      cp.align2 <- paste0("cp ", hybridOut2, "/align1/align1.xmap ", sandwich1.out, "/test_2_fullAln.xmap")
      run.system.cmd(cp.align1, run.flag)
      run.system.cmd(cp.align2, run.flag)
      #browser()
      flip.xmap(paste0(hybridOut1, "/align1/align1.xmap"), paste0(sandwich1.out, "/test_1_fullAln.xmap"))
      flip.xmap(paste0(hybridOut2, "/align1/align1.xmap"), paste0(sandwich1.out, "/test_2_fullAln.xmap"))

      print("Running pre-sandwich assembly")
      run.pre.sandwich(configs$BNGPath1, configs$BNGPath2, configs$bspqi_xmap, configs$bsssi_xmap, sandwich1.out, min.score = 13)#configs$twoEnzymeMergeT)#we should have another parameter for sandwich1 which is used for conflict resoltuion
      #setup output for combined conflict file
      run.system.cmd(paste0("mkdir ", combined.cut.conflicts), run.flag)
      print("Running merged conflict resolution")
      indata.file <- paste0(sandwich1.out, '/test_indata.RData')
      pos.offset.file <- paste0(sandwich1.out, '/test_pos_offsets.RData')
      conflicts <- run.conflict.resolve(hybridOut1, hybridOut2, indata.file, pos.offset.file, combined.conflicts1, combined.conflicts2)
      conflicts.saved <- paste0(combined.cut.conflicts, "/Combined_conflicts.RData")
      saveRDS(conflicts, conflicts.saved)
  }
  return(list(combined.conflicts1, combined.conflicts2))

}

#running hybrid scaffold to completion using the combined cut_conflcits files
run.hybrid.with.combined.conflict <- function(combined.conflict1, combined.conflict2, run.flag=T, configs){
  hybridOut1 <- paste0(configs$OutputDir, "/", configs$Enzyme1)
  hybridOut2 <- paste0(configs$OutputDir, "/", configs$Enzyme2)
  runHybrid1 <- paste0("perl ", configs$HybridScaffoldPath, " -r ", configs$RefAlignerPath, " -n ", configs$NGSPath,
                       " -b ", configs$BNGPath1, " -c ", configs$Hybrid.xml1, " -o ", hybridOut1, " -B 2 -N 2 ", " -M ", combined.conflict1,  " -t -a")
  runHybrid2 <- paste0("perl ", configs$HybridScaffoldPath, " -r " , configs$RefAlignerPath, " -n ", configs$NGSPath,
                       " -b ", configs$BNGPath2, " -c ", configs$Hybrid.xml2, " -o ", hybridOut2, " -B 2 -N 2 ", " -M ", combined.conflict2, " -t -a")
  #ret <- mclapply.with.error(list(runHybrid1, runHybrid2),
  #              function(cmd){run.system.cmd(cmd, run.flag)}, mc.cores=2)
  ret <- run.hybrid.with.error.handle(list(runHybrid1, runHybrid2),
                                      list(hybridOut1, hybridOut2), run.flag)

  return(list(hybridOut1, hybridOut2))
}

#obtain the filename without the file extension from bng file path
getBNG.file.header <- function(path){
    if(length(path) < 1){
        stop(paste0("Path  is empty"))
    }
    filename <- tail(strsplit(path, "/")[[1]],1) #get last entry in the path
    filename <- strsplit(filename, '\\.')[[1]][[1]] #strip the file extension
}

#strip file extension
stripExtension <- function(filename){
  ind <- tail(gregexpr('\\.', filename)[[1]],1)
  if(ind < 1){
     warning("file has not extension")
     return(filename)
  }
  substr(filename,1, ind-1)
}


#run the full sandwich assemblies
#using the ngs sequence as the glue molecules
run.ngs.sandwich <- function(hybridOut1, hybridOut2, run.flag=T, configs, run.ind=1){
    bng1.name <- getBNG.file.header(configs$BNGPath1)
    bng2.name <- getBNG.file.header(configs$BNGPath2)
    #temporarily for development test case
    #bng1.name <- substr(bng1.name, 7, nchar(bng1.name))
    #bng2.name <- substr(bng2.name, 7, nchar(bng2.name))

    ################
    #hybrid.alignfinal1 <- paste0("find ",  hybridOut1, "/hybrid_scaffolds_M", run.ind, "/ -name \"*NGS*HYBRID*.xmap\"")
    #hybrid.alignfinal2 <- paste0("find ",  hybridOut2, "/hybrid_scaffolds_M", run.ind, "/ -name \"*NGS*HYBRID*.xmap\"")
    hybrid.alignfinal1 <- paste0("find ",  hybridOut1, "/align_final_M", run.ind, "/ -name \"*NGS*HYBRID*_not_filtered.xmap\"")
    hybrid.alignfinal2 <- paste0("find ",  hybridOut2, "/align_final_M", run.ind, "/ -name \"*NGS*HYBRID*_not_filtered.xmap\"")
    hybrid1 <- paste0("find ",  hybridOut1, "/hybrid_scaffolds_M", run.ind, "/ -name \"*NGS*HYBRID*_r.cmap\"")
    hybrid2 <- paste0("find ",  hybridOut2, "/hybrid_scaffolds_M", run.ind, "/ -name \"*NGS*HYBRID*_r.cmap\"")

    hybrid.alignfinal1 <- run.system.cmd(hybrid.alignfinal1, run.flag)[[1]]
    hybrid.alignfinal2 <- run.system.cmd(hybrid.alignfinal2, run.flag)[[1]]
    hybrid1 <- run.system.cmd(hybrid1)[[1]]
    hybrid2 <- run.system.cmd(hybrid2)[[1]]

    sandwich2.out <- paste0(configs$TGHOut, "/Sandwich2")
    if(!run.flag){
        return(sandwich2.out)
    }
    run.system.cmd(paste0("mkdir ", sandwich2.out), run.flag)

    #cp.align1 <- paste0("cp ", hybrid.alignfinal1, " ", sandwich2.out, "/test_1_fullAln.xmap")
    #cp.align2 <- paste0("cp ", hybrid.alignfinal2, " ", sandwich2.out, "/test_2_fullAln.xmap")
    #run.system.cmd(cp.align1, run.flag)
    #run.system.cmd(cp.align2, run.flag)
    #we filter alignment base on "delta-alignment score" and overlap align regions
    print(hybrid.alignfinal1)
    align.filtered <- filterAlignFinal1(hybrid.alignfinal1, hybrid.alignfinal2, hybrid1, hybrid2, "")
    exportXmap(align.filtered[[1]], paste0(sandwich2.out, "/test_1_fullAln.xmap"))
    exportXmap(align.filtered[[2]], paste0(sandwich2.out, "/test_2_fullAln.xmap"))

    align2ref1 <- configs$bspqi_xmap #"/mnt/bionf_tmp/jwang/Assembly/V3_BSPQI_04062016_50xmappd_NonhapAsm/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap";
    align2ref2 <- configs$bsssi_xmap #"/home/users/jwang/data/V3_BSSSI_04062016_50xmappd_NonhapAsm/output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.xmap";
    print("Running full sandwich assembly")
    run.sandwich(hybrid.alignfinal1, hybrid.alignfinal2, hybrid1, hybrid2, sandwich2.out, align2ref1, align2ref2, score.T = configs$twoEnzymeMergeT)
    #getting summary stats from sandwich
    load(paste0(sandwich2.out, '/test_path_summmary.RData'))
    print(path_summary$scaffold_summary)
    return(sandwich2.out)
}


#converting the ngs sequence into a 2-color maps for final alignments
fa2cmap.twocolor <- function(hybridOut1, run.flag, configs){
   fa2cmap2color <- paste0(configs$TGHOut, "/fa2cmap")
   fa.cut.file <- paste0(basename(configs$NGSPath), '.cut.fasta')
   if(!dir.exists(fa2cmap2color)){
     run.system.cmd(paste0("mkdir ", fa2cmap2color), run.flag)
     agp.path <- paste0(hybridOut1, "/agp_fasta_M", configs$RunIndex, "/")
     key.file <- dir(path=agp.path, pattern = '*_key.txt.cut.txt', full.names = T)
     create.2color.cmap <- paste0("perl ", hybrid.script.dir, "/fa2cmap_multi_color.pl ", " -i ", agp.path, fa.cut.file,
                                  " -o ", fa2cmap2color, " -e ", configs$Enzyme1, " 1 ", configs$Enzyme2, " 2 ",  " -k ",  key.file)
     create.1color.cmap1 <- paste0("perl ", hybrid.script.dir, "/fa2cmap_multi_color.pl ", " -i ", agp.path, fa.cut.file,
                                   " -o ", fa2cmap2color, " -e ", configs$Enzyme2, " 1 ", " -k ", key.file)

     create.1color.cmap2 <- paste0("perl ", hybrid.script.dir, "/fa2cmap_multi_color.pl ", " -i ", agp.path, fa.cut.file,
                                   " -o ", fa2cmap2color, " -e ", configs$Enzyme1, " 1 ", " -k ", key.file)
     run.system.cmd(create.2color.cmap, run.flag)
     run.system.cmd(create.1color.cmap1, run.flag)
     run.system.cmd(create.1color.cmap2, run.flag)
   }else{
     warning(paste0(fa2cmap2color, " already exisits, skipping fa2cmap step"))
   }
   cmap2color <- paste0(fa2cmap2color, "/", stripExtension(fa.cut.file), "_", configs$Enzyme1, "_", configs$Enzyme2, "_0kb_0labels.cmap")
   cmap.color1 <- paste0(fa2cmap2color, "/", stripExtension(fa.cut.file), "_", configs$Enzyme1, "_0kb_0labels.cmap")
   cmap.color2 <- paste0(fa2cmap2color, "/", stripExtension(fa.cut.file), "_", configs$Enzyme2, "_0kb_0labels.cmap")
   return(list(cmap2color, cmap.color1, cmap.color2))
}


#this function obtain all contigs that are confidently aligned by single-color align-final and
#output the remaining unnaligned contigs for two-color alignmen
get.unalign.contigs <- function(single.align.file1, single.align.file2, qry.two.color, min.conf=13, min.lab=5){
    single.align1 <- as.data.table(readxmap(single.align.file1))
    single.align2 <- as.data.table(readxmap(single.align.file2))
    two.color.cmap <- as.data.table(readcmap(qry.two.color))
    align.ids <-unique(c(single.align1[Confidence > min.conf, QryContigID],
                                           single.align2[Confidence > min.conf, QryContigID]))
    print(paste0("Number of NGS contig single-enzyme alignment anchored: ", length(align.ids)))
    two.color.cmap.unalign <- two.color.cmap[!(CMapId %in% align.ids) & NumSites >= min.lab]
    print(paste0("Total unalign contigs with minimum labels ", nrow(two.color.cmap.unalign[SiteID==1])))
    return(two.color.cmap.unalign)
}

run.align.final <- function(hybridOut1, hybridOut2, sandwich2.out, RefAlignerPath, run.flag, configs){
    cmap.2color <- fa2cmap.twocolor(hybridOut1, run.flag, configs)
    align.final.out <- paste0(configs$TGHOut, "/alignfinal")
    run.system.cmd(paste0("mkdir ", align.final.out), run.flag)
    curr.wd <- getwd()
    setwd(align.final.out)
    ra <- configs$RefAlignerPath
    #aligning NGS to two-enzyme hybrid with single-color alignment
    single.color.qry <- list(cmap.2color[[2]], cmap.2color[[3]])
    single.color.ref <- paste0(sandwich2.out, list("/test_2col_c1.cmap", "/test_2col_c2.cmap"))
    enzymes <- list(configs$Enzyme1, configs$Enzyme2)
    run.singlecolor.aligns <- mclapply.with.error(1:length(single.color.qry),
                                    function(i){
                                        #prefix <- paste0(stripExtension(basename(single.color.qry[[i]])), "_1_errEst")
                                        prefix <- paste0(align.final.out, "/", "E_", enzymes[[i]], "_Q_NGScontigs_A_HYBRID")
                                        ret <- run_alignfinal_two_pass(single.color.ref[[i]], single.color.qry[[i]], prefix,
                                                                     run.flag, ra, bestRef = 0, filter = T,
                                                                     args1=list(configs$XML[['hybridScaffold1']][['align_final_1st_pass']],
                                                                               configs$XML[['hybridScaffold1']][['global']]),
                                                                     args2=list(configs$XML[['hybridScaffold1']][['align_final_2nd_pass']],
                                                                               configs$XML[['hybridScaffold1']][['global']])
                                                                     )
                                        return(paste0(prefix, '.xmap'))
                                    }, mc.cores = 2)

    #aligning BNG to two-enzyme hybrid
    bng.qry <- list(dir(paste0(hybridOut1, '/assignAlignType/cut_conflicts_M', configs$RunIndex, "/"), patter="bn_cut_pre_exclude.cmap", full.names = T),
                    dir(paste0(hybridOut2, '/assignAlignType/cut_conflicts_M', configs$RunIndex, "/"), patter="bn_cut_pre_exclude.cmap", full.names = T))

    run.bng.aligns <- mclapply.with.error(1:length(bng.qry),
                                                  function(i){
                                                    #prefix <- paste0(enzymes[i], "_", stripExtension(basename(bng.qry[[i]])), "_1_errEst")
                                                    prefix <- paste0("E_", enzymes[[i]], "_Q_BNG_A_HYBRID")
                                                    ret <- run_ra_alignref_final(bng.qry[[i]], single.color.ref[[i]], prefix,
                                                                                     run.flag, ra, bestRef=1, filter=F,
                                                                                     args=list(configs$XML[['hybridScaffold1']][['align_final_BNG']],
                                                                                               configs$XML[['hybridScaffold1']][['global']]))
                                                    return(paste0(prefix, '.xmap'))
                                                  }, mc.cores = 2)

    two.color.qry <-cmap.2color[[1]]
    two.color.ref <- paste0(sandwich2.out, "/test_2col_cut.cmap")
    print(run.singlecolor.aligns[[1]])
    print(run.singlecolor.aligns[[2]])
    print(two.color.qry)
    #we find out what contigs are aligned in previous step and get remaining contigs that are not aligned using single-color alignment
    two.color.unalign <- get.unalign.contigs(run.singlecolor.aligns[[1]], run.singlecolor.aligns[[2]], two.color.qry)
    two.color.qry <- paste0(stripExtension(basename(cmap.2color[[1]])), "_single_unaligned.cmap")
    if(run.flag){exportCMap(two.color.unalign, "./", two.color.qry)}
    #print("Start aligninig two-enzyme")
    if(nrow(two.color.unalign) > 0){
      run.twocolor.align <- lapply(1:2,
                                 function(hash.color){
                                     #prefix <- paste0(stripExtension(basename(cmap.2color[[1]])), "_twocolor_hash", hash.color, "_errEst")
                                     prefix <- paste0("E_", enzymes[[1]], "_E_",  enzymes[[2]], "_Q_NGS_A_HYBRID_hash", hash.color)
                                     ret <- run_ra_alignref_final_twocolor(two.color.ref, two.color.qry,
                                                                           prefix, run.flag, hash.color, ra,
                                                                           args = list(configs$XML[['TGH']][['align_final_twoEnzyme']]))
                                 })
    }
    print(getwd())

    #copying bng to ngs alignment to alignfinal folder
    bng2ngs <- list(dir(paste0(hybridOut1, '/align_final_M', configs$RunIndex, "/"), patter="prePairmerge_NGS", full.names = T),
                    dir(paste0(hybridOut2, '/align_final_M', configs$RunIndex, "/"), patter="prePairmerge_NGS", full.names = T))
    r <- lapply(bng2ngs, function(f){file.copy(f, align.final.out)})


    export.path <- get.align.finalstat(cmap.2color, run.flag, sandwich.out = sandwich2.out,
                                       T1 = configs$alignFinal1Color, T2 = configs$alignFinal2Color,
                                       enzyme1 = configs$Enzyme1, enzyme2=configs$Enzyme2)
    setwd(curr.wd)
    return(list(paste0(align.final.out, "/", export.path), cmap.2color))
    #print(getwd())
}

#merge multiple cmaps into a single one. Start with first cmap
#then successively, BNG maps not found in previous files were append to the file
get.merged.cmap <- function(cmap.file.list){
  if(length(cmap.file.list) < 1){
    stop("No cmap file list is specified")
  }
  print("Merging file list: ")
  print(cmap.file.list)
  cmap1 <- as.data.table(readcmap(cmap.file.list[[1]]))
  cmap.combined <- Reduce(function(cmap, f){
                  if(file.exists(f)){
                    cmap2 <- as.data.table(readcmap(f))
                    if(is.null(cmap2) || nrow(cmap2) == 0){
                      return(cmap)
                    }
                    cmap2.unique <- cmap2[!(CMapId %in% cmap1$CMapId)]
                    if(nrow(cmap2.unique)){
                      cmap <-rbind(cmap, cmap2.unique)
                    }
                    return(cmap)
                  }else{
                    warning(paste0("Cmap file ", f, " does not exists, skipping"))
                  }
               },
            cmap.file.list[-1], cmap1
         )
  return(cmap.combined)
}

get.align.finalstat <- function(cmap.2color, run.flag=T, sandwich.out, alignfinal.dir = './', T1=10, T2=12, enzyme1, enzyme2){
    align.file.prefix1 <- stripExtension(basename(cmap.2color[[1]]))
    align.file.prefix2 <- stripExtension(basename(cmap.2color[[2]]))
    align.file.prefix3 <- stripExtension(basename(cmap.2color[[3]]))
    two.color.hash1 <- dir(alignfinal.dir,  pattern="hash1.xmap$", full.names = T)
    two.color.hash2 <- dir(alignfinal.dir,  pattern="hash2.xmap$", full.names=T)
    #single.color <- dir(alignfinal.dir, pattern="Q_NGScontigs_A_HYBRID.xmap", full.names = T)
    single.color.1 <- dir(alignfinal.dir, pattern=paste0(enzyme1, "_Q_NGScontigs_A_HYBRID.xmap$"), full.names = T)
    single.color.2 <- dir(alignfinal.dir, pattern=paste0(enzyme2, "_Q_NGScontigs_A_HYBRID.xmap$"), full.names = T)
    aligns <- getAlignFinalStats(single.color.1, single.color.2, two.color.hash1, two.color.hash2, "", "", T1, T2)
    if(run.flag){
      saveRDS(aligns, paste0(align.file.prefix1, "_alignFinals.RData"))
      hybrid.combined <- get.merged.cmap(paste0(sandwich.out, "/", list("test_2col.cmap", "test_2col_c1.cmap", "test_2col_c2.cmap")))
      exportCMap(hybrid.combined, alignfinal.dir, paste0(align.file.prefix1, "_NGS_contigs_HYBRID_Export_r.cmap"))
      exportXmap(aligns$ExportAlign, paste0(alignfinal.dir, "/", align.file.prefix1, "_NGS_contigs_HYBRID_Export.xmap"))
      exportXmap(aligns$TwoColorAlign, paste0(alignfinal.dir, "/", align.file.prefix1, "_NGS_contigs_HYBRID_TwoEnzyme.xmap"))
      exportXmap(aligns$SingleOnlyAlign, paste0(alignfinal.dir, "/", align.file.prefix1, "_NGS_contigs_HYBRID_SingleEnzymeOnly.xmap"))
    }

    ngs.cmap <- as.data.table(readcmap(cmap.2color[[1]]))
    ngs.cmap.not.scaffold <- ngs.cmap[!(CMapId %in% aligns$ExportAlign$QryContigID)]
    scaffold.len <- aligns$Export[, RefLen[[1]], by=RefcontigID]
    non.scaffold.len <- ngs.cmap.not.scaffold[SiteID==1, ContigLength]
    if(length(non.scaffold.len)==0){
      non.scaffold.len <- c(0)
    }

    seq.summary <- data.frame(NGS_in_scaffold=length(unique(aligns$ExportAlign$QryContigID)),
                              n50_in_scaffold=getN50(aligns$Export$QryLen),
                              Len_in_scaffold=sum(aligns$Export$QryLen),
                              NGS_notIn_scaffold=length(unique(ngs.cmap.not.scaffold$CMapId)),
                              n50_notIn_scaffold=getN50(non.scaffold.len),
                              len_notIn_scaffold=sum(non.scaffold.len),
                              n50_overall_SeqAsm=getN50(c(scaffold.len$V1, non.scaffold.len)))
    print(seq.summary)
    summary.stat <- aligns$SummaryStats
    summary.stat$SeqSummary <- seq.summary
    if(run.flag){saveRDS(summary.stat, paste0(alignfinal.dir, "/summary_stats.RData"))}
    return(paste0(align.file.prefix1, "_NGS_contigs_HYBRID_Export.xmap"))
}


run.agp.export <- function(export.align.path, cmap.2color, hybridOut1, sandwich.out, run.flag, configs){
    agp.out <- paste0(configs$TGHOut, "/AGPExport/")
    align.final <- paste0(configs$TGHOut, "/alignfinal/")
    run.system.cmd(paste0("mkdir ", agp.out), run.flag)
    fasta.file <- run.system.cmd(paste0('find ', hybridOut1, "/agp_fasta_M", configs$RunIndex,  "/ -name *cut*.fasta"))
    cmap.name.map.file <- paste0(stripExtension(cmap.2color[[1]]), "_key.txt")
    #copy input file to export dir, need for final output
    run.system.cmd(paste0("cp ", fasta.file, " ", agp.out))

    merged.hybrid.cmap <- dir(align.final, pattern = "_NGS_contigs_HYBRID_Export_r.cmap", full.names = T)
    file.copy(from=merged.hybrid.cmap, to = agp.out)
    agp.export <- paste0("perl ", hybrid.script.dir, "/ExportAGP_TwoEnzyme.pl ",
                         " -i ", export.align.path, " -m ", cmap.name.map.file,
                         " -c ", merged.hybrid.cmap,
                         " -s ", fasta.file,
                         " -o ", agp.out, ' -e \"', configs$Enzyme1, " 1 ", configs$Enzyme2, " 2 \"")
                         #" -o ", agp.out, ' -e \"N 1 N 2 \"')
    #print(agp.export)
    print(paste0("Trimm option " , configs$trim))
    if(configs$trim){
      agp.export <- paste0(agp.export, " -t ")
    }
    run.system.cmd(agp.export, run = run.flag)

    #generate manifest file for NCBI export
    get.manifest.file(agp.out, configs)

    return(agp.out)
}

#this generates the manifest file for agp export
get.manifest.file <- function(agp.out, configs){
  cut.fasta <- dir(agp.out, pattern = '*cut.fasta$')
  agp.file <- dir(agp.out, pattern = '*.agp$')
  coord.file <- dir(agp.out, pattern ='*.coord$')
  hybrid.cmap <- dir(agp.out, patter='*_r.cmap$')
  conn <- file(paste0(configs$TGHOut, '/ncbi_manifest.txt'), open = 'w')
  writeLines(paste0('./AGPExport/', list(cut.fasta, agp.file, coord.file, hybrid.cmap)), conn)
  flush(conn)
  close(conn)
}

#get the stats from agp/fasta export
get.agp.stat <- function(export.path){
  agp.file <- dir(export.path, pattern = '*.agp', full.names = T)
  if(length(agp.file) < 1){
    stop(paste0("Cannot find the agp file in export directory: ", export.path))
  }
  agp <- as.data.table(read.table(agp.file, sep = '\t', header = F))
  total.hybrids <- agp[grepl('Super-Scaffold_', V1), length(unique(V1))]
  hybrid.lens <- agp[grepl('Super-Scaffold_', V1), max(V3), by='V1']
  unscaffold.lens <- agp[(!grepl('Super-Scaffold_', V1)), max(V3), by='V1']

  total.hybrid.len <- sum(as.numeric(hybrid.lens[[2]]))
  hybrid.n50 <- getN50(as.numeric(hybrid.lens[[2]]))
  all.n50 <- getN50(as.numeric(c(hybrid.lens[[2]], unscaffold.lens[[2]])))
  total.unscaffold.len <- sum(as.numeric(unscaffold.lens[[2]]))

  return(list(total.hybrids=total.hybrids,
              unscaffold.len=total.unscaffold.len,
              hybrid.len=total.hybrid.len,
              hybrid.n50=hybrid.n50,
              all.n50 = all.n50))
}

get.single.enzyme.report <- function(report.file, stats.list=NULL, printHeader=T){
  if(is.null(stats.list)){
    stats.list <- list("Bpp-adjusted BioNano Genome Map statistics:",
                     "Original NGS Genome Map statistics:")
                     #"BNG contigs in hybrid scaffold statistics:",
                     #"NGS contigs in hybrid scaffold statistics:",
                     #"Hybrid scaffold statistics:",
                     #"The statistics of hybrid scaffold plus not scaffolded NGS:")
  }
  report <- read.table(report.file, sep="\t", stringsAsFactors = F)
  report.selected <- unlist(lapply(stats.list, function(stat.header){
                              ind <- which(grepl(stat.header, report$V1))
                              if(length(ind) > 0){
                                if(printHeader){
                                  return(list(report$V1[[ind]],
                                              paste0("  ", report$V1[[ind+1]]),
                                              paste0("  ", report$V1[[ind+5]]),
                                              paste0("  ", report$V1[[ind+7]])))
                                }else{
                                  return(list(
                                              paste0("  ", report$V1[[ind+1]]),
                                              paste0("  ", report$V1[[ind+5]]),
                                              paste0("  ", report$V1[[ind+7]])))
                                }
                              }else{
                                return("")
                              }
                          }))
  return(report.selected)
}

#this function generate a summary output for TGH scaffold
get.scaffold.report <- function(cmap2color, configs, hybridOuts, agp.out){
  conflicts.dir <- paste0(configs$OutputDir, "/CombinedConflictsCut")
  alignfinal.dir <- paste0(configs$TGHOut, "/alignfinal")
  sandwich.dir <- paste0(configs$TGHOut, "/Sandwich2")
  conflicts <- readRDS(paste0(conflicts.dir, "/Combined_conflicts.RData"))
  report.file <- paste0(configs$TGHOut, "/hybrid_scaffold_informatics_report.txt")
  hybrid1.report <- paste0(hybridOuts[[1]], "/hybrid_scaffolds_M", configs$RunIndex, "/hybrid_scaffold_informatics_report.txt")
  hybrid2.report <- paste0(hybridOuts[[2]], "/hybrid_scaffolds_M", configs$RunIndex, "/hybrid_scaffold_informatics_report.txt")
  conn <- file(report.file, open = 'w')

  conflict.summary <- conflicts[[3]]$summary
  conflict.overlap <- conflicts[[3]]$overlap
  writeLines("Automatic Conflicts Detections from BNG-NGS alignment: ", conn)
  writeLines("", conn)
  #writeLines("After merging (or manually edited) conflicts: ", conn)
  combined.conflict <- as.data.table(read.hybrid.conflictfile(configs$ManualCut1))
  combined.conflict2 <- as.data.table(read.hybrid.conflictfile(configs$ManualCut2))
  writeLines(paste0("Number of conflict cuts made to Bionano maps (", configs$Enzyme1, "): ", combined.conflict[, sum(qry_leftBkpt_toCut == 'cut') + sum(qry_rightBkpt_toCut == 'cut')]), conn)
  writeLines(paste0("Number of conflict cuts made to Bionano maps (", configs$Enzyme2, "): ", combined.conflict2[, sum(qry_leftBkpt_toCut == 'cut') + sum(qry_rightBkpt_toCut == 'cut')]), conn)
  writeLines(paste0("Number of conflict cuts made to NGS sequences: ",
                    combined.conflict[, sum(ref_leftBkpt_toCut == 'cut') + sum(ref_rightBkpt_toCut == 'cut')]
                    - conflict.overlap$NGS_cut_overlap), conn)
  writeLines(paste0("Number of BNG maps to be cut (", configs$Enzyme1, "): ", combined.conflict[qry_leftBkpt_toCut == 'cut' | qry_rightBkpt_toCut == 'cut', length(unique(qryId))]), conn)
  writeLines(paste0("Number of BNG maps to be cut (", configs$Enzyme2, "): ", combined.conflict2[qry_leftBkpt_toCut == 'cut' | qry_rightBkpt_toCut == 'cut', length(unique(qryId))]), conn)
  writeLines(paste0("Number of NGS contigs to be cut: ", combined.conflict[ref_leftBkpt_toCut == 'cut' | ref_rightBkpt_toCut == 'cut', length(unique(refId))]), conn)
  writeLines("", conn)

  #writeLines(paste0("Stats for ", configs$Enzyme1), conn)
  single.report1 <- get.single.enzyme.report(hybrid1.report, stats.list = list("Original NGS sequences statistics"), printHeader = F)
  writeLines("Input NGS sequence statistics: ", conn)
  lapply(single.report1, function(l){writeLines(l, conn)})
  writeLines("", conn)


  single.report1 <- get.single.enzyme.report(hybrid1.report, stats.list = list("Original BioNano Genome Map statistics"), printHeader = F)
  writeLines(paste0("Input Bionano Genome map statistics (", configs$Enzyme1, ")"), conn)
  lapply(single.report1, function(l){writeLines(l, conn)})
  writeLines("", conn)

  #writeLines(paste0("Stats for ", configs$Enzyme2), conn)
  single.report2 <- get.single.enzyme.report(hybrid2.report, stats.list = list("Original BioNano Genome Map statistics"), printHeader = F)
  writeLines(paste0("Input Bionano Genome map statistics (", configs$Enzyme2, ")"), conn)
  lapply(single.report2, function(l){writeLines(l, conn)})

  align.final.summary <- readRDS(paste0(configs$TGHOut, "/alignfinal/summary_stats.RData"))
  writeLines("", conn)
  writeLines("Hybrid scaffold stats: ", conn)
  writeLines("", conn)
  load(paste0(sandwich.dir, "/test_path_summmary.RData"))
  #writeLines(paste0("Total number of single-enzyme hybrids: ", path_summary$scaffold_summary$n_nodes), conn)
  #writeLines(paste0("Total number of two-enzyme hybrids: ", path_summary$scaffold_summary$n_scaffolded), conn)
  #writeLines(paste0("Total number of leftover single-enzyme hybrids: ", path_summary$scaffold_summary$n_non_scaffolded), conn)

  writeLines(paste0("Total number of hybrids: ",
                   path_summary$scaffold_summary$n_scaffolded
                 + path_summary$scaffold_summary$n_non_scaffolded), conn)
  #writeLines(paste0("Scaffold N50 of single-enzyme hybrid1 (Mb): ", path_summary$scaffold_summary$asm1_n50/1e6), conn)
  #writeLines(paste0("Scaffold N50 of single-enzyme hybrid2 (Mb): ", path_summary$scaffold_summary$asm2_n50/1e6), conn)
  #writeLines(paste0("Hybrid scaffold N50 (Mbp)  : ", path_summary$scaffold_summary$scaffolded_n50/1e6), conn)
  writeLines(paste0("Two-enzyme hybrid scaffold N50 (Mbp): ", path_summary$scaffold_summary$scaffolded_n50/1e6), conn)
  #writeLines(paste0("Mean length of  final hybrids (two-enzyme + leftover): ", path_summary$scaffold_summary$overall_mean_len/1e6), conn)
  #writeLines(paste0("N50 fold improvement over single-enzyme hybrid: ", path_summary$scaffold_summary$fold_improvement), conn)
  #writeLines(paste0("Total length of two-enzyme hybrids (Mb): ", path_summary$scaffold_summary$scaffolded_size/1e6), conn)
  writeLines(paste0("Scaffold N50 of final hybrids inclusive of single-enzyme hybrids (Mbp): ", path_summary$scaffold_summary$overall_n50/1e6), conn)
  writeLines(paste0("Total length of hybrid scaffolds (Mbp): ",
                    path_summary$scaffold_summary$scaffolded_size/1e6
                    + path_summary$scaffold_summary$non_scafflded_size/1e6), conn)

  writeLines(paste0("Total length of unscaffoled NGS (Mbp): ", align.final.summary$SeqSummary$len_notIn_scaffold/1e6), conn)

  #writeLines(paste0("Total length of left-over hybrids (Mb): ", path_summary$scaffold_summary$non_scafflded_size/1e6), conn)
  #writeLines(paste0("Longest single-enzyme hybrid (Mb): ", path_summary$scaffold_summary$longest_map_len_i/1e6), conn)
  writeLines(paste0("Longest hybrid (Mbp): ", path_summary$scaffold_summary$longest_map_len_f/1e6), conn)

  writeLines(paste0(""), conn)


  #writeLines(paste0("Total Number of NGS anchored in two-enzyme hybrids: ", align.final.summary$sandwich_ngs_align$sandwich_total), conn)
  #writeLines(paste0("Anchored with two-enzyme: ", align.final.summary$sandwich_ngs_align$sandwich_twoEnzy), conn)
  #writeLines(paste0("Anchored with Enzyme1 only: ", align.final.summary$sandwich_ngs_align$sandwich_enzy1_only), conn)
  #writeLines(paste0("Anchored with Enzyme2 only: ", align.final.summary$sandwich_ngs_align$sandwich_enzy2_only), conn)
  #writeLines(paste0("Total length of NGS anchored in two-enzyme hybrids (Mb): ", align.final.summary$sandiwch_ngs_len$sandwich_total/1e6), conn)
  #writeLines(paste0("Anchored with two-enzyme (Mb): ", align.final.summary$sandiwch_ngs_len$sandwich_twoEnzy/1e6), conn)
  #writeLines(paste0("Anchored with Enzyme1 only (Mb): ", align.final.summary$sandiwch_ngs_len$sandwich_enzy1_only/1e6), conn)
  #writeLines(paste0("Anchored with Enzyme2 only (Mb): ", align.final.summary$sandiwch_ngs_len$sandwich_enzy2_only/1e6), conn)
  #writeLines(paste0(""), conn)
  writeLines(paste0("Total Numbers of NGS anchored in final hybrids: ", align.final.summary$hybrid_ngs_align$hybrid_total), conn)
  writeLines(paste0("Anchored in two-enzyme hybrid: ", align.final.summary$hybrid_ngs_align$hybrid_twoEnzy), conn)
  writeLines(paste0("Anchored in single-enzyme hybrid1 only: ", align.final.summary$hybrid_ngs_align$hybrid_enzy1_only), conn)
  writeLines(paste0("Anchored in single-enzyme hybrid2 only: ", align.final.summary$hybrid_ngs_align$hybrid_enzy2_only), conn)
  writeLines(paste0("Total NGS not in scaffold: ", align.final.summary$SeqSummary$NGS_notIn_scaffold), conn)
  writeLines("", conn)
  writeLines(paste0("Total length of NGS anchored in final hybrids (Mbp): ", align.final.summary$hybrid_ngs_Len$hybrid_total/1e6), conn)
  writeLines(paste0("N50 of sequence anchored (Mbp): ", align.final.summary$SeqSummary$n50_in_scaffold/1e6), conn)
  writeLines(paste0("Anchored with two-enzyme (Mbp): ", align.final.summary$hybrid_ngs_Len$hybrid_twoEnzy/1e6), conn)
  writeLines(paste0("Anchored with single-enzyme hybrid1 only (Mbp): ", align.final.summary$hybrid_ngs_Len$hybrid_enzy1_only/1e6), conn)
  writeLines(paste0("Anchored with single-enzyme hybrid2 only (Mbp): ", align.final.summary$hybrid_ngs_Len$hybrid_enzy2_only/1e6), conn)
  writeLines(paste0("Total NGS not in scaffold (Mbp): ", align.final.summary$SeqSummary$len_notIn_scaffold/1e6), conn)
  writeLines(paste0("N50 of sequence not in scaffold (Mbp): ", align.final.summary$SeqSummary$n50_notIn_scaffold/1e6), conn)
  writeLines(paste0("Scaffold N50 of final hybrids inclusive of unscaffolded NGS (Mbp): ", align.final.summary$SeqSummary$n50_overall_SeqAsm/1e6), conn)
  #writeLines(paste0("N50 of Hybrid + Seq not Scaffold: ", align.final.summary$SeqSummary$n50_overall_SeqAsm/1e6), conn)
  writeLines("", conn)
  agp_stats <- get.agp.stat(agp.out)
  writeLines(paste0("Stats calculated from final FASTA file: "), conn)
  writeLines("", conn)
  writeLines(paste0("Total number of hybrid (Mbp): ", agp_stats$total.hybrids), conn)
  writeLines(paste0("Hybrid scaffold n50 (Mbp): ", agp_stats$hybrid.n50/1e6), conn)
  writeLines(paste0("Scaffold N50 of final hybrid inclusive of unscaffolded NGS sequences (Mbp): ", agp_stats$all.n50/1e6), conn)
  writeLines(paste0("Total length of hybrid scaffold (Mbp): ", agp_stats$hybrid.len/1e6), conn)
  writeLines(paste0("Total length of unscaffoled NGS (Mbp): ", agp_stats$unscaffold.len/1e6), conn)

  flush(conn)
  close(conn)

}

writeLinesNow <- function(string, con){
  writeLines(string, con)
  flush(con)
}

#try parse the file to see if it is xml file
isXML <- function(file){
  tryCatch({
    tgh.params <- xmlParse(file)
    return(TRUE)
  }, error=function(e){
     print(paste0("Parsed xml error ", e))
     return(FALSE)
  }
  )
}

#parse parameter file
parseParams <- function(parsed.args){
  params.file <- parsed.args$paramFile
  params <- tryCatch({return(generateHybridXml(parsed.args))},
                      error=function(e){
                         print(e)
                         stop(paste0("Cannot parse parameter XML file: ", params.file, " please check format."))
                      })
  return(params)
}

#checking validity of input parameters
checkParams <- function(params){
  #params$HybridScaffoldPath <- HybridScaffoldPath
  input.files <- list(params$HybridScaffoldPath, params$RefAlignerPath, params$NGSPath, params$BNGPath1, params$BNGPath2, params$OutputDir, params$Hybrid.xml1, params$Hybrid.xml2)
  checks <- unlist(lapply(input.files, file.exists))
  if(!all(checks)){
    stop(paste0("Input file ", input.files[which(!checks)], " cannot be found\n"))
  }
  if(is.null(params$run.flags)){
    params$run.flags <- c(T, T, T, T, T, T)
  }else{
    if(length(params$run.flags) != 6){
      stop("Run flags for some stage of the pipeline is missing, please provide a correct run flag vector or do not use run flags to run the whole pipeline.")
    }
  }
  params <- check.manual.cuts(params)
  return(params)
}

#check manual cut file and set run flag accordingly
check.manual.cuts <- function(params){
  if(is.na(params$ManualCut1)){
    return(params)
  }else{
    if(!file.exists(params$ManualCut1)){
      stop(paste0("Conflict cut file ", params$ManualCut1, " not valid"))
    }
    if(!file.exists(params$ManualCut2)){
      stop(paste0("Conflict cut file ", params$ManualCut2, " not valid"))
    }
  }
  #skipping first two step in pipeline if manual cut file present
  params$run.flags[1] <- F
  params$run.flags[2] <- F
  return(params)
}


#TGH may be re-ran with maunal conflict cut for multiple iteration, find the current iteration
getRunIndex <- function(hybridOut, configs){
  hybrid.out.path <- dir(hybridOut, pattern = 'hybrid_scaffolds_M*')
  if(length(hybrid.out.path)==0){
    return(0)
  }
  split.hybrid.out <- unlist(strsplit(hybrid.out.path, 'hybrid_scaffolds_M'))
  max.ind <- max(as.double(Filter(function(x){x!=""}, split.hybrid.out)))
  return(max.ind)
}

#merge table by their key, if key exist in both table, take value of
#first table. is.valid function check if values from t1 is valid
#if not valid it will use value from t2 instead
merge.table <- function(t1, t2, is.valid = function(x){!is.null(x)}){
  keys <- unique(c(names(t1), names(t2)))
  combined.table <- c()
  combined.table[keys] <- t1[keys]
  key.t1.not.valid <- names(Filter(function(x){is.null(x)}, combined.table))
  key.t1.not.valid <- key.t1.not.valid[key.t1.not.valid != ""]
  combined.table[key.t1.not.valid] <- t2[key.t1.not.valid]
  if(all(sapply(combined.table, is.valid))){
    return(combined.table)
  }else{
   stop(paste0("Cannot find valid values for following key in the table:\n", "t1:\n", toString(t1), "\n", "t2:\n", toString(t2)))
  }
}


copy.cut.coor <- function(params, hybridOuts){
  conflict_dir1 <- paste0(hybridOuts[[1]], '/assignAlignType/cut_conflicts_M', params$RunIndex, "/")
  conflict_dir2 <- paste0(hybridOuts[[2]], '/assignAlignType/cut_conflicts_M', params$RunIndex, "/")
  NGS_coord1 <- dir(conflict_dir1, pattern='auto_cut_NGS_coord_translation')
  NGS_coord2 <- dir(conflict_dir2, pattern='auto_cut_NGS_coord_translation')
  BNG_coord1 <- dir(conflict_dir1, pattern='auto_cut_BN_coord_translation')
  BNG_coord2 <- dir(conflict_dir2, pattern='auto_cut_BN_coord_translation')
  file.copy(paste0(conflict_dir1, "/", NGS_coord1),
            paste0(params$TGHOut, "/CombinedConflictsCut/", params$Enzyme1, "_", NGS_coord1))
  file.copy(paste0(conflict_dir1, "/", BNG_coord1),
            paste0(params$TGHOut, "/CombinedConflictsCut/", params$Enzyme1, "_", BNG_coord1))
  file.copy(paste0(conflict_dir2, "/", NGS_coord2),
            paste0(params$TGHOut, "/CombinedConflictsCut/", params$Enzyme2, "_", NGS_coord2))
  file.copy(paste0(conflict_dir2, "/", BNG_coord2),
            paste0(params$TGHOut, "/CombinedConflictsCut/", params$Enzyme2, "_", BNG_coord2))

}

#this function create a tar ball .of the hybrid-scaffold results for IrysView
create.tar <- function(params, hybridOuts){
  #copy neccessary log file to final result folder
  log.files <- list(params$log.file, paste0(params$log.file, ".errLog"), paste0(params$OutputDir, "/", params$status))
  ret <- lapply(log.files, function(f){
    run.system.cmd(paste0("cp ", f, " ",
                   params$TGHOut, '/'))
   })
  #copy conflict cut status file
  conflicts.dir <- paste0(params$TGHOut, '/CombinedConflictsCut/')
  dir.create(conflicts.dir)
  run.system.cmd(paste0("cp ", params$ManualCut1, " ", conflicts.dir))
  run.system.cmd(paste0("cp ", params$ManualCut2, " ", conflicts.dir))
  tar.name <- strsplit(params$tar, split='\\.')[[1]][[1]]
  #copy conflict cut alignment file
  align1.dir <- paste0(params$TGHOut, '/align1')
  dir.create(align1.dir)
  align1.files1 <- dir(paste0(hybridOuts, '/align1/'), pattern='align1.*map')
  #copy align1 files to output folder
  copy.align1 <- function(path, enzyme){
    align1.files <- dir(paste0(path, '/align1/'), pattern='align1*')
    r <- lapply(align1.files1,
                function(f){
                  new.str <- paste0('E_', enzyme, '_Q_preCut_BNG_A_preCut_NGS')
                  new.f <- sub('align1', new.str, f)
                  #run.system.cmd(paste0("cp " , path, "/align1/", f, " ", align1.dir, "/", new.f))
                  run.system.cmd(paste0("sed -- \'s/align1_/", new.str, "_/g\' ",
                                         path, "/align1/", f,
                                         " > ", align1.dir, "/", new.f))
                })
  }
  r1 <- copy.align1(hybridOuts[[1]], params$Enzyme1)
  r2 <- copy.align1(hybridOuts[[2]], params$Enzyme2)
  copy.cut.coor(params, hybridOuts)
  #align1.file1 <- dir(paste0(hybridOuts[[1]], '/hybrid_scaffolds_M', params$RunIndex, "/"), pattern='*BNGcontigs_NGScontigs*', full.names = T)
  #run.system.cmd(paste0("cp ", align1.file1, " ", align1.dir, "/"))
  #align1.file2 <- dir(paste0(hybridOuts[[2]], '/hybrid_scaffolds_M', params$RunIndex, "/"), pattern='*BNGcontigs_NGScontigs*', full.names = T)
  #run.system.cmd(paste0("cp ", align1.file2, " ", align1.dir, "/"))

  #tar.path <- paste0(params$OutputDir, "/", tar.name)
  #dir.create(tar.path)
  #file.copy(from=params$TGHOut, to=tar.path, recursive = T)
  #run.system.cmd(paste0("cp -r ", params$TGHOut, " ", tar.path))
  con <- file(paste0(params$TGHOut, '/', 'cur_results.txt'), open = 'w')
  writeLinesNow(paste0('{Version:"', params$Version, '"}'), con)
  writeLinesNow(paste0(params$TGHOut, '/alignfinal'), con)
  close(con)
  try({
   #run.system.cmd(paste0("tar -czvf ", params$OutputDir, "/", params$tar, " -C ", params$OutputDir, "/ ",
   #                     'two_enzyme_hybrid_scaffold_M', params$RunIndex))
   
   pwd <- getwd()
   setwd(params$OutputDir)

   run.system.cmd(paste0("zip -r ", params$OutputDir, "/", params$tar, " ", 
                        'two_enzyme_hybrid_scaffold_M', params$RunIndex))
  })
  #unlink(tar.path, recursive = T)

}

#main function
runHybridSandwich <- function(parsed.args){
  params <- parseParams(parsed.args)
  params$log.file <- paste0(params$OutputDir, "/TGH.log")
  errlog <- file(paste0(params$log.file, '.errLog'), open='wt')
  sink(params$log.file)
  sink(errlog, type="message")
  begin.time <- Sys.time()

  print(paste0("Begn timestamp ", begin.time))
  print(params)
  print("from command line")
  print(parsed.args)
  params <- merge.table(params, parsed.args, function(x){TRUE})
  print('after merge')
  print(params)
  params$Version <- "3.0" #specify hybrid-scaffold version to spceified in output
  begin.time <- Sys.time()
  print(paste0("Begin timestamp ", begin.time))

  run.stat <- tryCatch({

    status.file <- file(paste0(params$OutputDir, "/", params$status), open = 'w')
    writeLinesNow("Checking input parameter file", status.file)
    params <- checkParams(params)

    print(paste0("Running stage ", paste(seq(1:6), ": ", params$run.flags)))

    #step 1: run ngs to bng alignment and detect conlfict using hybridscaffold
    run.flag1 <- params$run.flags[[1]]
    stage <- "Aligning BNG to NGS and detecting conflicts"
    writeLinesNow(paste0("{Stage: \"", stage, "\", Status:\"started\"}"), status.file)
    writeLines(paste0("{Stage: \"", stage, "\", Status:"), status.file, sep="")
    hybridOuts <- run.NGS2BNG.align(run.flag1, params)
    writeLinesNow("\"done\"}", status.file)

    #step 2: a) merge cut conflict files from each hybridscaffold pre-run
    #        b) detect and resolve conflicts/inconsistency between the two bng assemblies
    run.flag2 <- params$run.flags[[2]]
    stage <- "Merging conflicts from BNG-NGS alignment and detecting conflicts between two sets of BNG maps"
    writeLinesNow(paste0("{Stage: \"", stage, "\", Status:\"started\"}"), status.file)
    writeLines(paste0("{Stage: \"", stage, "\", Status:"), status.file, sep="")
    combined.conflicts <- resolve.merged.conflicts(hybridOuts[[1]], hybridOuts[[2]], run.flag2, params)
    #updated conflict cut file
    params$ManualCut1 <- combined.conflicts[[1]]
    params$ManualCut2 <- combined.conflicts[[2]]
    writeLinesNow("\"done\"}", status.file)

    #step 3: run hybrid-scaffold with combined conflict cut file
    run.flag3 <- params$run.flags[[3]]
    stage <- "Constructing single-enzyme hybrid scaffold"
    writeLinesNow(paste0("{Stage: \"", stage, "\", Status:\"started\"}"), status.file)
    writeLines(paste0("{Stage: \"", stage, "\", Status:"), status.file, sep="")
    hybridOuts <- run.hybrid.with.combined.conflict(combined.conflicts[[1]], combined.conflicts[[2]], run.flag3, params)
    writeLinesNow("\"done\"}", status.file)

    #create output dir for storing TGH output
    run.ind <- getRunIndex(hybridOuts[[1]], params)
    params$TGHOut <- paste0(params$OutputDir, "/two_enzyme_hybrid_scaffold_M", run.ind)
    params$RunIndex <- run.ind
    dir.create(params$TGHOut)
    print(paste0("Two enzyme output dir: ", params$TGHOut))

    #step 4: run sandwich assemblies using ngs as glue
    run.flag4 <- params$run.flags[[4]]
    stage <- "Merging two single-enzyme scaffolds into two-enzyme scaffolds"
    writeLinesNow(paste0("{Stage: \"", stage, "\", Status:\"started\"}"), status.file)
    writeLines(paste0("{Stage: \"", stage, "\", Status:"), status.file, sep="")
    sandwich.out <- run.ngs.sandwich(hybridOuts[[1]], hybridOuts[[2]], run.flag=run.flag4, params, run.ind)
    writeLinesNow("\"done\"}", status.file)

    #step 5: aligning ngs sequence to two-color maps using both two-color and single-color alignments
    run.flag5 <- params$run.flags[[5]]
    stage <- "Aligning NGS sequence contigs to final scaffold"
    writeLinesNow(paste0("{Stage: \"", stage, "\", Status:\"started\"}"), status.file)
    writeLines(paste0("{Stage: \"", stage, "\", Status:"), status.file, sep="")
    export.align <- run.align.final(hybridOuts[[1]], hybridOuts[[2]], sandwich.out, RefAlignerPath, run.flag5, params)
    writeLinesNow("\"done\"}", status.file)

    #step 6: export output to agp and fasta format
    run.flag6 <- params$run.flags[[6]]
    stage <- "Generating AGP and fasta file for final scaffold"
    writeLinesNow(paste0("{Stage: \"", stage, "\", Status:\"started\"}"), status.file)
    writeLines(paste0("{Stage: \"", stage, "\", Status:"), status.file, sep="")
    agp.out <- run.agp.export(export.align[[1]], export.align[[2]], hybridOuts[[1]], sandwich.out, run.flag6, params)
    writeLinesNow("\"done\"}", status.file)

    #writing out summary report
    report <- get.scaffold.report(export.align[[2]], params, hybridOuts, agp.out)

    #wrapping up
    writeLines("", status.file)
    writeLinesNow("Two-enzyme hybrid scaffold pipeline finished running!", status.file)
    #copying status file to TGH folder so each run/re-run keeps their own logs

    end.time <-  Sys.time()
    print(paste0("End timestamp ", end.time))
    print(paste0("Total running time for two-enzyme hybrid scaffold pipeline: ", (end.time-begin.time)))
    run.stat <- 0
  #handling error from running pipeline
  },  error = function(err){
    writeLines("Hybrid scaffold pipeline encountered errors:", status.file)
    writeLinesNow(paste("ERROR: ", err), status.file)
    print(paste("ERROR: ", err))
    return(1)
  #Ending hybrid-scaffold
  }, finally = {
    sink()
    sink(type="message")
    close(status.file)
    #archiving results in a tar file
  })
  print(paste0("run stat: ", run.stat))
  if(run.stat == 0){
    create.tar(params, hybridOuts)
  }
  return(run.stat)
}

#running the pipeline
run <- function(){
  args <- commandArgs()
  begin.ind <- which(grepl("--args", args))
  if(any(grepl("--help", args)) || length(begin.ind) == 0 || length(args) < begin.ind+1){
      print(paste0("Usage: Rscript runHybridScaffoldV3.R  <ParamFile>"))
  }else{
    param.file <- args[begin.ind + 1]
    runHybridSandwich(param.file)
  }
}
