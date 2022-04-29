#20151028 EL - functions to run RefAligner commands

#merge maps based on a list file
run_ra_merge <- function(cmap_list,prefix,num_colors=1,force_flag=T,run_flag=F) {
#   ra <- "/home/users/elam/RefAligner/ra_150811/bin/RefAligner"
  prefix <- paste0(prefix,"s")

  base_command <- paste(ra,"-if",cmap_list,"-o",prefix,"-merge")

  if (num_colors==1) {
    command <- base_comamnd
  } else if (num_colors==2) {
    command <- paste(base_comamnd,"-colors 2")
  } else {
    stop("the number of colors has to be either 1 or 2.")
  }

  if (force_flag) {
    command <- paste(command,"-f")
  }

  stdout <- NULL
  cmap <- paste0(prefix,".cmap")

  if (run_flag) {
    suppressWarnings({
      stdout <- system(command,ignore.stdout=F,ignore.stderr=F,intern=T)
    })
  } else {
    #     print(command)
  }

  outlist <- list("command"=command,
                  "stdout"=stdout,
                  "cmap"=cmap)

  return(outlist)
  #return(stdout)
}

#run alignref final alignment (single-color)
#10/5/2016: updated alignment parameter to be consistent with the newest hybrid-scaffold align-final parameters
run_ra_alignref_final <- function(cmap,ref,prefix,run_flag=F, ra=NULL, bestRef=1, filter=F, args=NULL) {
    if(is.null(ra)){
        #ra <- "/home/users/elam/RefAligner/ra_150811/bin/RefAligner"
        stop(paste0("Cannot find valid refaligner path: ", ra))
    }

  stdout <- NULL
  #if no argument specified we used the default one
  if(any(sapply(args, is.null))){
    args <- paste(
      "-stdout","-stderr","-f",
      "-maxthreads 64",
      "-minsites 5 -A 5",
      "-res 2.9 -resSD 0.75 -FP 0.6 -FN 0.06 -sf 0.20 -sd 0.0 -sr 0.01 -extend 1 -outlier 1e-4 -endoutlier 1e-3 -PVendoutlier",
      "-deltaX 12 -deltaY 12 -xmapchim 12 -relerr 0.001",
      "-hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 -hashMultiMatch 100 -insertThreads 4",
      "-nosplit 2 -biaswt 0 -T 1e-5 -S -1000 -indel -PVres 2 -MaxSE 0.5 -HSDrange 1.0 -outlierBC",
      "-xmapUnique 12 -AlignRes 2. -outlierExtend 12 24 -Kmax 12 -maxmem 128 -maxvirtmem 0",
      "-BestRef ", 0, " ",
      "-RepeatMask 4 0.01 -RepeatRec 0.7 0.6 1.4 -indel -AlignRes 2 -resEstimate",
      "-mres 0.9 -rres 0.9 -MultiMatches 5"
    )
  }else{
    args <- lapply(args, function(a){get_ref_align_args(a)})
    args <- paste(args, collapse = " ")
  }
  command <- paste(ra,"-i",cmap,"-ref",ref,"-o",prefix, args)

  stdout <- run.system.cmd(command, run=run_flag, print=F)

  if(filter){
    align_out <- paste0(prefix, '.xmap')
    process_align_final(align_out)
  }

  outlist <- list("command"=command,
                  "stdout"=stdout,
                  "prefix"=prefix)

  return(outlist)
  #return(stdout)
}

run_ra_bng_alignref_final <- function(cmap,ref,prefix,run_flag=F, ra=NULL, bestRef=1) {
  if(is.null(ra)){
    #ra <- "/home/users/elam/RefAligner/ra_150811/bin/RefAligner"
    stop(paste0("Cannot find valid refaligner path: ", ra))
  }

  stdout <- NULL
  command <- paste(ra,"-i",cmap,"-ref",ref,"-o",prefix,
                   "-stdout","-stderr","-f",
                   "-maxthreads 64",
                   "-minsites 5 -A 5",
                   "-res 2.9 -resSD 0.75 -FP 0.6 -FN 0.06 -sf 0.20 -sd 0.0 -sr 0.01 -extend 1 -outlier 0.0001 -endoutlier 0.001 -PVendoutlier",
                   "-deltaX 6 -deltaY 6 -xmapchim 12 -relerr 0.001",
                   "-hashgen 5 7 2.4 1.5 0.05 5.0 1 1 3 -hash -hashdelta 26 10 46 -hashMultiMatch 30 10 -insertThreads 4",
                   "-nosplit 2 -biaswt 0 -T 1e-9 -S -1000 -indel -PVres 2 -MaxSE 0.5 -HSDrange 1.0 -outlierBC",
                   "-xmapUnique 12 -AlignRes 2. -outlierExtend 12 24 -Kmax 12",
                   "-BestRef ", bestRef, " ",
                   "-RepeatMask 4 0.01 -RepeatRec 0.7 0.6 1.4 -indel -AlignRes 2 -resEstimate",
                   "-mres 0.9 -rres 0.9 -maxmem 128 -RAmem 3 3 -maxvirtmem 0" )

  if (run_flag) {
      #stdout <- system(command,ignore.stdout=F,ignore.stderr=F,intern=T)
     stdout <- run.system.cmd(command, run = run_flag, print = F)
  } else {
    #     print(command)
  }

  outlist <- list("command"=command,
                  "stdout"=stdout,
                  "prefix"=prefix)

  return(outlist)
  #return(stdout)
}


run_ra_alignref_final_twocolor <- function(cmap,ref,prefix,run_flag=F, hash_color=1, ra=NULL, args=NULL) {
   if(is.null(ra)){
       #ra <- "/home/users/elam/RefAligner/ra_150811/bin/RefAligner"
       stop(paste0("Cannot find valid refaligner path: ", ra))
   }

  if(any(sapply(args, is.null))){
    args <- paste(
      "-stdout","-stderr","-f",
      "-maxthreads 64 -colors 2",
      "-minsites 5 -A 5",
      "-res 2.9 -resSD 0.75 -FP 0.6 -FN 0.06 -sf 0.20 -sd 0.0 -sr 0.01 -extend 1 -outlier 0.0001 -endoutlier 0.001 -PVendoutlier",
      "-deltaX 12 -deltaY 12 -relerr 0.001",
      "-hashgen 5 2 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 -mres 0.9 -hashMultiMatch 100 -insertThreads 4",
      "-nosplit 2 -biaswt 0 -T 1e-9 -S -1000 -indel -PVres 2 -rres 0.9 -MaxSE 0.5 -HSDrange 1.0 -outlierBC",
      "-xmapUnique 12 -AlignRes 2 -Kmax 12 -maxmem 128 -maxvirtmem 0",
      "-indel -AlignRes 2 -resEstimate -MultiMatches 5 -RAmem 3 3"
     )
   }else{
      args <- lapply(args, function(a){get_ref_align_args(a)})
      args <- paste(args, collapse = " ")
   }

   command <- paste(ra,"-ref", ref,"-i", cmap,"-o", prefix, "-hashcolor", hash_color, args)


  if (run_flag) {
      #stdout <- system(command,ignore.stdout=F,ignore.stderr=F,intern=T)
      stdout <- run.system.cmd(command, run = run_flag, print = F)
   } else {
     #     print(command)
   }
   outlist <- list("command"=command,
                  "stdout"=stdout,
                  "prefix"=prefix)

  return(outlist)
  #return(stdout)
}

#get a one-color cmap from a two-color cmap
run_ra_get_color <- function(cmap,prefix,channel,run_flag=F) {
  prefix <- paste0(prefix,"_ch",channel)

  command <- paste(ra,"-i",cmap,"-o",prefix,"-merge","-colors 2","-usecolor",channel,"-f")
  stdout <- NULL

  if (run_flag) {
    suppressWarnings({
      stdout <- system(command,ignore.stdout=F,ignore.stderr=F,intern=T)
    })
  } else {
    #     print(command)
  }

  cmap <- paste0(prefix,".cmap")

  outlist <- list("command"=command,
                  "stdout"=stdout,
                  "cmap"=cmap)

  return(outlist)
}

#run 2-color alignment
run_ra_2col <- function(cmap,ref='',prefix,run_flag=F) {
#   ra <- "/home/users/elam/RefAligner/ra_151016/bin/RefAligner"
#   ref <- "/home/users/csecol/genomes/human/hg19/hg19_BSSSI_BSPQI_0kb_0labels.cmap" #bsssi/bspqi reference

  stdout <- NULL

  command <- paste(ra,"-i",cmap,"-ref",ref,"-o",prefix,
                   "-colors 2 -rres 0.9 -M 1 -biaswt 0 -MaxSR 0.10 -MaxSE 0.5 -BestRef 1 -S -10000 -A 2 -T 1e-12 -BestRefPV 1",
                   "-resEstimate -resbias 4.0 32 -AlignRes 2.0 -Kmax 12 -outlierRate",
                   "-hashgen 5 7 2.4 1.5 0.05 5.0 1 1 3 -hash -hashdelta 10 -insertThreads 16 -hashcolor 2 -hashmaxmem 128",
                   "-nosplit 2 -outlier 1e-3 -endoutlier 1e-3",
                   "-deltaX 6 -deltaY 6",
                   "-FP 0.2 0.2 -FN 0.02 0.02 -sf 0.2 0.2 -sd 0 0 -sr 0.02 0.02 -se 0.1 0.1 -res 2.9 2.9 -resSD 0.7 0.7 -XmapInterleaved 1 -maxmem 1280")

  if (run_flag) {
    suppressWarnings({
      stdout <- system(command,ignore.stdout=F,ignore.stderr=F,intern=T)
      #       print("stdout")
      #       print(stdout)
    })
  } else {
    #     print(command)
  }

  outlist <- list("command"=command,
                  "stdout"=stdout,
                  "prefix"=prefix)

  return(outlist)
  #return(stdout)
}
