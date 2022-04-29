#This file contain function that process xml-based documents
#this function read the TGH xml parameter file and generates the corresponding xml file
#for single-enzyme hybrids
generateHybridXml <- function(parsed.args){
  print(parsed.args$paramFile)
  tgh.params <- xmlParse(parsed.args$paramFile)
  tgh.params <- xmlRoot(tgh.params)

  params <- parseTGHParams(tgh.params)
  params <- merge.table(params, parsed.args, function(x){return(TRUE)})
  params$XML <- tgh.params

  #converting path to absolute
  #we need some more generalize way to distinguish path string from regular character argument
  #right now it require the path to have at least one path separator
  params <- lapply(params, function(p){
                            #if(is.character(p) && grepl(.Platform$file.sep, p) && file.exists(p)){
                            if(is.character(p) && file.exists(p)){
                              normalizePath(p)
                            }else{
                              p
                            }
                        })
  #browser()
  #if override we wipe out outptu directory first
  if(parsed.args$override){
    unlink(paste0(params$OutputDir, "/*"), recursive = T)
  }
  #generating single-enzyme hyrid xml parameters file
  if(is.na(params$Enzyme1)){
    params$Enzyme1 <- getChildAttrVal(tgh.params[["hybridScaffold1"]][["fasta2cmap"]], 'enzyme')
  }else{
    fasta2cmap1 <- getChildWithAttr(tgh.params[["hybridScaffold1"]][["fasta2cmap"]], 'enzyme')[[1]]
    xmlAttrs(fasta2cmap1)["val0"] <- params$Enzyme1
  }
  if(is.na(params$Enzyme2)){
    params$Enzyme2 <- getChildAttrVal(tgh.params[["hybridScaffold2"]][["fasta2cmap"]], 'enzyme')
  }else{
    fasta2cmap2 <- getChildWithAttr(tgh.params[["hybridScaffold2"]][["fasta2cmap"]], 'enzyme')[[1]]
    xmlAttrs(fasta2cmap2)["val0"] <- params$Enzyme2
  }

  if(is.null(params$OutputDir)){
    stop(paste0("Cannot find output path in parameter file: ", tgh.xml))
  }

  if(!file.exists(params$OutputDir)){
    system(paste0('mkdir -p ', params$OutputDir))
    params$OutputDir <- normalizePath(params$OutputDir)
  }

  xml.out <- paste0(params$OutputDir, "/Config_file")
  dir.create(xml.out)
  tgh.xml.copy <- paste0(xml.out, "/TGH.xml")
  saveXML(tgh.params, tgh.xml.copy)

  #saving xml for single-enzyme hybrid
  hybrid.xml1 <- paste0(xml.out, "/hybrid_", params$Enzyme1, ".xml")
  hybrid.xml2 <- paste0(xml.out, "/hybrid_", params$Enzyme2, ".xml")
  #changing section name
  hybrid1.params <- tgh.params[["hybridScaffold1"]]
  xmlName(hybrid1.params) <- "hybridScaffold"
  saveXML(hybrid1.params, hybrid.xml1)
  xmlName(hybrid1.params) <- "hybridScaffold1" #restore proper name

  hybrid2.params <- tgh.params[["hybridScaffold2"]]
  xmlName(hybrid2.params) <- "hybridScaffold"
  saveXML(hybrid2.params, hybrid.xml2)
  xmlName(hybrid2.params) <- "hybridScaffold2" #restore proper name

  params <- c(params, TGH.xml = tgh.xml.copy,
              Hybrid.xml1 = hybrid.xml1,
              Hybrid.xml2 = hybrid.xml2)
  #we also store teh root of xml parameters
  return(params)
}

#parse TGH-specific parameters
parseTGHParams <- function(xmlroot){
  tgh.params <- xmlroot
  #parsing TGH required parameters
  tgh.section <- tgh.params[["TGH"]][["common"]]
  TGH.required <- list(
    NGSPath = getChildAttrVal(tgh.section, 'NGSPath'),
    RefAlignerPath = getChildAttrVal(tgh.section, 'RefAlignerPath'),
    OutputDir = getChildAttrVal(tgh.section, 'OutputDir'),
    BNGPath1 = getChildAttrVal(tgh.section, 'BNGPath1'),
    BNGPath2 = getChildAttrVal(tgh.section, 'BNGPath2'),
    Enzyme1 = getChildAttrVal(tgh.section, 'Enzyme1'),
    Enzyme2 = getChildAttrVal(tgh.section, 'Enzyme2')
  )
  #print(paste0("Merge p-value: ", TGH.required$twoEnzymeMergeT))
  #check to make sure all required parameters are there
   params <- lapply(1:length(TGH.required),
                         function(i){
                           if(length(TGH.required[[i]]) < 1){
                             #stop(paste0("Cannot find required parameter: ", names(TGH.required)[[i]],  ", please check xml file: "))
                             return(NULL)
                           }else{
                             return(TGH.required[[i]][[1]])
                           }
                         })
  names(params) <- names(TGH.required)
  #parsing optional parameters

  TGH.optional <- list(
   RunFlags = getChildAttrVal(tgh.section, 'RunFlags'),
   alignFinal2Color = getChildAttrVal(tgh.section, 'alignFinal2Enzyme'),
   twoEnzymeMergeT= getChildAttrVal(tgh.section, 'twoEnzymeMergeT')
  )

  if(length(TGH.optional$RunFlags) > 0){
    RunFlags <- strsplit(TGH.optional$RunFlags[[1]], ",")[[1]]
    RunFlags <- as.logical(RunFlags)
    params <- c(params, run.flags=list(RunFlags))
  }
  if(length(TGH.optional$alignFinal2Color) > 0){
    alignFinal2Color <- -log10(as.double(TGH.optional$alignFinal2Color[[1]]))
    params <- c(params, alignFinal2Color=alignFinal2Color)
  }else{
    params <- c(params, alignFinal2Color=12)
  }
  if(length(TGH.optional$twoEnzymeMergeT) > 0){
    params <- c(params, twoEnzymeMergeT=-log10(as.numeric(TGH.optional$twoEnzymeMergeT)))
  }else{
    params <- c(params, twoEnzymeMergeT=13)
  }
  #getting alignfinal thershold from single-enzyme section
  alignfinal1.T <- getChildAttrVal(tgh.params[["hybridScaffold1"]][["align_final_2nd_pass"]], 'T')
  alignfinal2.T <- getChildAttrVal(tgh.params[["hybridScaffold2"]][["align_final_2nd_pass"]], 'T')
  if(length(alignfinal1.T) > 0 && length(alignfinal2.T) > 0){
      params <- c(params, alignFinal1Color = max(-log10(as.double(alignfinal1.T)),
                                                 -log10(as.double(alignfinal2.T))))
  }else{
    params <- c(params, alignFinal1Color = 11)
  }

  return(params)
}

#this function return the refaligner argument from a particular xml section
get_ref_align_args <- function(align_section){
  children <- xmlChildren(align_section)
  args_list <- lapply(children,
                      function(child){
                        attrs <- xmlAttrs(child)
                        option_name_ind <- which(names(attrs)=='attr')
                        if(length(option_name_ind) > 1){
                          stop(paste0("cannot have more than one option in each argument line: ", attrs))
                        }
                        args <- paste0("-", attrs[option_name_ind])
                        option_values_inds <- which(grepl('val', names(attrs)))
                        option_values_inds <- sort(option_values_inds)
                        option_values <- ""
                        if(length(option_values_inds)){
                          option_values <- paste0(attrs[option_values_inds], collapse = " ")
                        }
                        args <- paste0(args, " ", option_values)
                      })
  args <- paste(args_list, sep="", collapse=" ")
  return(args)
}

#this function search and return the children of a node
#with a specified attribute
getChildWithAttr <- function(xmlnode, attrName){
  children <- xmlChildren(xmlnode)
  children <- lapply(children,
                     function(child){
                       if(xmlAttrs(child)[['attr']] == attrName)
                         {return(child)}})
  children <- Filter(function(x){!is.null(x)}, children)
}

getChildAttrVal <- function(xmlnode, attrName){
  children.nodes <- getChildWithAttr(xmlnode, attrName)
  attribs <- lapply(children.nodes, function(child){
                        attribs <- xmlAttrs(child)
                        if('val0' %in% names(attribs)){
                          return(attribs[['val0']])
                        }else{
                          return(NULL)
                        }})
  attribs <- Filter(function(x){!is.null(x)}, attribs)
}

testProcessXMl <- function(){
  tgh.xml <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/ReRun/TGH.xml"
  single.xml.out <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/Config_Xml/"
  args <- list(paramFile = tgh.xml, OutputDir=single.xml.out)
  params <- generateHybridXml(args)
}

testAlignParamsFromXML <- function(){
  tgh.xml <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/ReRun/Config_file/TGH.xml"
  single.xml.out <- "~/RemoteServer//home/users5/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestGreedyAlignFinal/TGH_two_passes.xml"
  params <- xmlParse(tgh.xml)
  params <- xmlRoot(params)

  p <- parseTGHParams(params)

  align_final <- params[['hybridScaffold1']][['align_final']]
  get_ref_align_args(align_final)

  align0 <- params[['hybridScaffold1']][['align1']]
  get_ref_align_args(align0)
}



