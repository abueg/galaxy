#this file contain function that process xmap alignments
getTranslateIndices <- function(xmap, qcmap, rcmap, qcmap.2color, rcmap.2color, label.channel=1){
  xmap <- processXmapAlignments(xmap)
  xmap.translated <- xmap
  qcmap.2color$GlobalSiteID <- qcmap.2color$SiteID
  rcmap.2color$GlobalSiteID <- rcmap.2color$SiteID
  qcmap.split <- getSeparateChannel(qcmap.2color)
  rcmap.split <- getSeparateChannel(rcmap.2color)

  if(label.channel == 1){
    rcmap.selected <- rcmap.split$C1
    qcmap.selected <- qcmap.split$C1
  }else if(label.channel==2){
    rcmap.selected <- rcmap.split$C2
    qcmap.selected <- qcmap.split$C2
  }else{
    stop("Unknown channel can only be 1 or 2")
  }
  #translating per-channel indices to global indices
  global.ref.inds <- lapply(seq(1:nrow(xmap)),
                              function(i){
                                rcmap.selected[CMapId==xmap$RefcontigID[[i]]][xmap$alignRefInd[[i]], GlobalSiteID]
                              })

  global.qry.inds <- lapply(seq(1:nrow(xmap)),
                            function(i){
                               qcmap.selected[CMapId==xmap$QryContigID[[i]]][xmap$alignQryInd[[i]], GlobalSiteID]
                            })

  xmap$globalRefInds <- global.ref.inds
  xmap$globalQryInds <- global.qry.inds
  #generating new alignment string
  align.new <- mapply(function(x,y){paste0("(", x, ",", y, ")", sep="", collapse="")}, xmap$globalRefInds, xmap$globalQryInds)
  xmap$Alignment <- align.new
  xmap$alignRefInd <- NULL
  xmap$alignQryInd <- NULL
  xmap$globalRefInds <- NULL
  xmap$globalQryInds <- NULL
  return(xmap)
}

test.translate.indice <- function(){
  xmap <- readxmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/TGH_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID.xmap')
  qcmap <- readcmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/TGH_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID_q.cmap')
  rcmap <- readcmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/TGH_M1/alignfinal/E_BSPQI_Q_NGScontigs_A_HYBRID_r.cmap')
  qcmap.2color <- as.data.table(readcmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/TGH_M1/fa2cmap/NGS_2col_mres0.9.cmap'))
  cmap.2color.list <- dir('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/TGH_M1/Sandwich2/', pattern = "test_2col_", full.names = T)
  out.file <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/TGH_M1/alignfinal/test_merged.xmap"
  rcmap.2color <- as.data.table(get.merged.cmap(rev(cmap.2color.list)))
  rcmap.2color <- as.data.table(readcmap('~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/TGH_M1/Sandwich2/hybrids_merged_mre0.9.cmap'))
  exportXmap(xmap[RefcontigID==1], out.file)
  exportCMap(rcmap.2color, "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_PB40x/ReRunTestCmdArgs/TGH_M1/Sandwich2/", 'hybrids_merged.cmap')

}







