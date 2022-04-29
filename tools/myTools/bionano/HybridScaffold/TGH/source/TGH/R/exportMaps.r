################################################################################
################################################################################
################################################################################
################################################################################
exportCMap <- function( reads, folder, fileName, headers=NULL, printColName=T) {
    sink( paste0( folder, fileName ) );
    labelChannels <- max(reads$LabelChannel)
    labelChannels <- ifelse(is.infinite(labelChannels), 1, labelChannels)  #in case it is an empty cmap
    if(is.null(headers)){
      cat( "# CMAP File Version:\t0.1\n" );
      cat( paste0("# Label Channels:\t", labelChannels, "\n" ));
      cat( "# Nickase Recognition Site 1:\tunknown\n" );
      cat( "# Number of Consensus Nanomaps:\t24\n" );
    }else{
      lapply(headers, function(h){cat(h); cat("\n")})
    }
    if(printColName){
      names <- colnames(reads)
      names[[1]] <- paste0("#h ", names[[1]])
      header1 <- paste(names, collapse = "\t")
      cat(paste0(header1,"\n"))
      cat( "#f\tint\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\n");
    }
    sink()
    #formatString <- "%s\t%15.1f\t%15i\t%15i\t%15i\t%15.1f\t%15.1f\t%15i\t%15i\n";
    try( {
      write.table(reads, paste0(folder, fileName), sep="\t", row.names = F, col.names = F, quote=F, append = T)
    } ); # try
} # exportCMap


exportCmap2Bnx <- function(molecules, outputPath){
    print(paste0("Output path: ", outputPath))
    sink( outputPath);
    labelChannels <- max(molecules$LabelChannel)
    cat( "# BNX File Version: 0.1\n")
    cat( paste0("# Label Channels:     ", labelChannels, "\n" ));
    cat( "# Nickase Recognition Site 1: unknown\n" );
    try({
        prevId = "";
        beginInd=1;
        max.channel <- max(molecules$LabelChannel)
        for(n in 1:nrow(molecules)){
            if(molecules$CMapId[n] != prevId){
                if(prevId != ""){
                    cat( paste0("0\t", sprintf("%.0f", as.numeric(molecules$CMapId[n-1])),
                                "\t", sprintf("%.3f", molecules$ContigLength[n-1]), "\n"))
                    for(channel in 1:max.channel){
                        cat(channel)
                        molInds <- beginInd:(n-1)
                        select.channel <- which(molecules$LabelChannel[molInds] == channel
                                                | molecules$LabelChannel[molInds] == 0) #include end label in all each channel
                        #no labels of this specific channel, add in a dummy begin label
                        if(length(select.channel) == 1){
                            cat( paste0("\t20\t",
                                        sprintf("%.3f", molecules$Position[molInds[select.channel]])))
                        }else{
                            cat( paste("\t",
                                       sprintf("%.3f", molecules$Position[molInds[select.channel]])))
                        }
                        cat("\n")
                    }
                }
                prevId <- molecules$CMapId[n]
                beginInd <- n
                #print(paste0("begin is: ", beginInd))
            }
        }
    })
    sink()
}

#export simple xmap
exportXmap <- function(xmap, out.file, headers=NULL, printColName=T){
  if(!is.null(out.file)){
    names <- colnames(xmap)
    names[[1]] <- paste0("#h ", names[[1]])
    header1 <- paste(names, collapse = "\t")
    header2 <- "#f int\tint\tint\tfloat\tfloat\tfloat\tfloat\tstring\tfloat\tstring\tfloat\tfloat\tint\tstring"
    sink(out.file)
    if(!is.null(headers)){
      lapply(headers, function(h){cat(h); cat("\n")})
    }else{
      cat("# XMAP File Version:\t0.2\n")
    }
    if(printColName){
      cat(header1)
      cat('\n')
      cat(paste0(header2, "\n"))
      sink()
    }
    #setnames(xmap, old.name, new.name)
    write.table(xmap, out.file,
                sep = "\t", row.names = F, quote = F, col.names = F, append = T)
    #etnames(xmap, new.name, old.name)
  }
}

