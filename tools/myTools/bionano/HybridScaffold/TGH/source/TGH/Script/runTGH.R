#This script serve as the main interface to Two-enzyme hybrid scaffold pipeline

#set up the R libray paths and required dependency
setup.library <- function(script.dir, path.sep, install=FALSE){
    repos <- "http://cran.cnr.berkeley.edu"
    library.path <- paste0(script.dir, path.sep, 'TGH', path.sep, 'library')
    dependent.packages <- list('data.table', 'igraph', 'intervals', 'MASS', 'XML', 'argparser')
    r <- 0
    if(!file.exists(library.path)){
        r <- system(paste0('mkdir ', library.path))
    }
    #setting up library path variable
    library.path <- paste0(script.dir, path.sep, "TGH", path.sep, "library/")
    .libPaths(c(library.path, .libPaths()))

    #check for dependecy and install dependent library if neccesary
    if(r == 0){
        for(package in dependent.packages){
            if(!require(package, character.only=TRUE)){
                install.packages(package, lib=library.path, repos = repos, dependencies=TRUE)
            }
        }
        if(install){
          install.packages(paste0(script.dir, path.sep, "TGH", path.sep, "TGH_0.1.tar.gz"), lib=library.path, dependencies=TRUE)
        }
    }else{
        stop(paste0('Cannot create library path ', library.path))
    }
}

#parse command line arguments and grep the script path
get.script.paths <- function(args){
    #return(list(script.path="/home/users//jwang/workspace/codebase/HybridScaffoldThree/",
    #            param.file="/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRun10102016/runHybridParams.R",
    #            log.file="/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRun10102016/TGH.log"))
    begin.ind <- which(grepl("--args", args))
    dir.sep <- .Platform$file.sep
    script.path.ind <- which(grepl("--file=", args))
    script.path <- strsplit(args[[script.path.ind]], "=")[[1]][[2]]
    script.path <- strsplit(script.path, "runTGH.R")[[1]]
    if(script.path == ""){
        script.path <- paste0('.', dir.sep)
    }

    if(length(args) == begin.ind || length(begin.ind) == 0){
      return(list(script.path=script.path, args=c('-h')))
    }
    return(list(script.path=script.path, args=args[(begin.ind+1):length(args)]))
}

#parsing command line argument for options additional options and argument
parse.cmd.args <- function(args){
  parser <- arg_parser(description = "This script performs two-enzyme hybrid scaffolding of BNG genome maps and NGS sequence",
                       name = "runTGH.R")
  parser <- add_argument(parser, "paramFile", help="parameter file (xml)")
  parser <- add_argument(parser, "--BNGPath1", help="Path to BNG maps for enzyme1", short = '-b1')
  parser <- add_argument(parser, "--BNGPath2", help="Path to BNG maps for enzyme2", short = '-b2')
  parser <- add_argument(parser, "--NGSPath", help="Path to NGS sequence (fasta file)")
  parser <- add_argument(parser, "--OutputDir", help="Output directory", default = './')
  parser <- add_argument(parser, "--RefAlignerPath", help="Path to RefAligner")
  parser <- add_argument(parser, "--RunFlags", help="Specify the stage of scaffold to run/skip")
  parser <- add_argument(parser, "--Enzyme1", help="Enzyme used in BNG maps specify by --bng1", short = '-e1')
  parser <- add_argument(parser, "--Enzyme2", help="Enzyme used in BNG maps specify by --bng2", short = '-e2')
  parser <- add_argument(parser, "--ManualCut1", help="Manual cut file", short = '-m1')
  parser <- add_argument(parser, "--ManualCut2", help="Manual cut file", short = '-m2')
  parser <- add_argument(parser, "--tar", help="Result tar file to be import to IrysView", default = "TGH.tar")
  parser <- add_argument(parser, "--status", help="Path to status file", default="status.txt")
  parser <- add_argument(parser, "--override", help="Override output folder", short='-f', flag = TRUE)
  parser <- add_argument(parser, "--install", help="Install/re-install TGH packages", short='-i', flag = TRUE)
  parsed.args <- parse_args(parser, args)
  return(parsed.args)
}


run.main <- function(){
  args <- commandArgs()
  print("Command line args:")
  print(args)
  install.flag <- any(grepl('--install', args) || any(grepl('-i', args)))
  #get script paths to where other required code reside
  parsed.args <- get.script.paths(args)
  dir.sep <- .Platform$file.sep
  script.path <- parsed.args$script.path
  HybridScaffoldPath <- paste0(script.path, dir.sep,  "hybridScaffold.pl")
  hybrid.script.dir <<- paste0(script.path, dir.sep, "scripts", dir.sep)
  print(paste0("TGH script path ", script.path))
  print(paste0("HybridScaffold path ", HybridScaffoldPath))
  print(paste0("HybridScaffold script path ", hybrid.script.dir))

  #check for required library, if not exists insall them
  print("Setting up libraries")
  setup.library(script.path, dir.sep, install.flag)
  print("Current library paths: ")
  print(.libPaths())
  library(TGH)

  print("Parsing command line arguments")
  parsed.args <- parse.cmd.args(parsed.args$args)

  HybridScaffoldPath <- paste0(script.path, dir.sep,  "hybridScaffold.pl")
  hybrid.script.dir <- paste0(script.path, dir.sep, "scripts", dir.sep)
  Hybrid.xml1 <- paste0(script.path, dir.sep, "TGH", dir.sep, "hybridScaffold_config_aggressive_TGH.xml")
  Hybrid.xml2 <- paste0(script.path, dir.sep, "TGH", dir.sep, "hybridScaffold_config_aggressive_TGH2.xml")
  #we probably should check if xml and script path exists
  #print(paste0("Hybrid xml path ", Hybrid.xml1))
  parsed.args$HybridScaffoldPath <- HybridScaffoldPath
  print("Starting two-enzyme hybrid scaffold")
  r <- runHybridSandwich(parsed.args)
  print("Finish scaffold.")
  quit(save="no", status=r, runLast=FALSE)
}

run.main()
