#this script is a wrapper for running sandwich assemblies
scripts_dir <- "/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/sandwichScript/"
#initialize dependent pacakge and dependent script
source(paste0(scripts_dir,"initialize.R"))

args <- commandArgs()
begin.ind <- which(grepl("--args", args))
if(any(grepl("--help", args)) || length(begin.ind) == 0 || length(args) < begin.ind+5){
    print(paste0("Usage: Rscript run_sandwich_wrapper.R <align1> <align2> <asm1> <asm2> <output>"))
    stop()
}

paths <- lapply(args[(begin.ind+1):(begin.ind+7)], normalizePath) #convert input argument to absolute paths
print(paste0("Input paths: "))
print(paths)
run.sandwich(paths[[1]], paths[[2]], paths[[3]], paths[[4]], paths[[5]], paths[[6]], paths[[7]])
