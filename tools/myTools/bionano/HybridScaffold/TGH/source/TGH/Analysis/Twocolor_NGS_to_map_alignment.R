#This script compare the alignment rate of NGS to reference using BspqI, BsssI and two-color alignments
library(data.table)


analyzeNGSAligns <- function(){
#xmap.hybrid.bspqi <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridScaffold_NonHaplotype/ReRun_03242016/BSPQI_with_relaxParams/align_final/BSPQI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"))
#xmap.hybrid.bsssi <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/HybridScaffold_NonHaplotype/ReRun_03242016/BSSSI_with_relaxParams/align_final/BSSSI_EXP_REFINEFINAL1_bppAdjust_cmap_a_lines_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap"))
#xmap.hybrid.twocolor <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies_withHybridM3/AlignFinal/TwoColorHash/Hash1/a.lines_BSPQI_BSSSI_0kb_0labels_twocolor_errEst.xmap"))
#xmap.hybrid.twocolor2 <- as.data.table(readxmap("~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSandwich/NonHaplotype_assemblies_withHybridM3/AlignFinal/TwoColorHash/Hash2/a.lines_BSPQI_BSSSI_0kb_0labels_twocolor_errEst.xmap"))

xmap.hybrid.bspqi <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSPQI_0kb_0labels_1_errEst.xmap"))
xmap.hybrid.bsssi <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSSSI_0kb_0labels_1_errEst.xmap"))
xmap.hybrid.twocolor <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash1_errEst.xmap"))
xmap.hybrid.twocolor2 <- as.data.table(readxmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/alignfinal/a.lines.fasta.cutted_BSPQI_BSSSI_0kb_0labels_twocolor_hash2_errEst.xmap"))


xmap.hybrid.bspqi <- xmap.hybrid.bspqi[xmap.hybrid.bspqi$Confidence > 11,]
xmap.hybrid.bsssi <- xmap.hybrid.bsssi[xmap.hybrid.bsssi$Confidence > 11,]



xmap.hybrid.twocolor$Confidence <- xmap.hybrid.twocolor$Confidence + 0.01
xmap.hybrid.twocolor <- rbind(xmap.hybrid.twocolor, xmap.hybrid.twocolor2)

xmap.hybrid.twocolor <- xmap.hybrid.twocolor[Confidence > 11]
setkeyv(xmap.hybrid.twocolor, c("QryContigID", "Confidence")) #need to sort it in order of the selector below
best.score <- xmap.hybrid.twocolor[,Confidence == max(Confidence), by=QryContigID] #picking best alignment 
xmap.hybrid.twocolor <- xmap.hybrid.twocolor[best.score$V1==T]


qcmap <- readcmap("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunV3Auto/ReRun08232016/fa2cmap2color/a.lines.fasta.cutted_BSPQI_0kb_0labels.cmap")
qcmap <- as.data.table(qcmap)

bspqi.labelcounts <- qcmap[,sum(LabelChannel==1), by=CMapId]
bsssI.labelcounts <- qcmap[,sum(LabelChannel==2), by=CMapId]

mol.length <- qcmap[,log10(ContigLength[[1]]), by=CMapId]


setkey(bspqi.labelcounts, CMapId)
molstats <- bspqi.labelcounts[bsssI.labelcounts]
setkey(molstats, CMapId)
molstats <- molstats[mol.length]
setkey(molstats, CMapId)
setnames(molstats, c("CMapId", "NumSites.bspqi", "NumSites.bsssi", "log10.Length"))


label.counts.bin <- c(seq(1, 100, 1), 500)
mol.length.bin <- c(seq(4, 6, 0.025))

molstats$bspqi.bin.ind <- findInterval(molstats$NumSites.bspqi, label.counts.bin)
molstats$bsssi.bin.ind <- findInterval(molstats$NumSites.bsssi, label.counts.bin)
molstats$twocol.bin.ind <- findInterval(molstats$NumSites.bspqi + molstats$NumSites.bsssi, label.counts.bin)
molstats$length.bin <- findInterval(molstats$log10.Length, mol.length.bin)

h.bspqi <- table(molstats$bspqi.bin.ind)[-1]
h.bsssi <- table(molstats$bsssi.bin.ind)[-1]
h.twocolor <- table(molstats$twocol.bin.ind)[-1]
h.mollen <- table(molstats[log10.Length > 4, length.bin])[-1]
plot(as.numeric(h.mollen))
plot(as.numeric(h.bspqi))
plot(h.bsssi)

setkey(molstats, CMapId)

#getting aligned molecules
aligned.bspqi <- data.table(CMapId = unique(xmap.hybrid.bspqi$QryContigID))
bspqi.aligned.labelcount <- molstats[aligned.bspqi]

aligned.bsssi <- data.table(CMapId = unique(xmap.hybrid.bsssi$QryContigID))
bsssi.aligned.labelcount <- molstats[aligned.bsssi]

aligned.twocolor <- data.table(CMapId = unique(xmap.hybrid.twocolor$QryContigID))
twocol.aligned.labelcount <- molstats[aligned.twocolor]

aligned.union <- data.table(CMapId = unique(c(aligned.bsssi$CMapId, aligned.bspqi$CMapId)))
union.aligned.labelcount <- molstats[aligned.union]

h.bspqi.aligned <- table(bspqi.aligned.labelcount$bspqi.bin.ind)
fract.aligned <- h.bspqi.aligned/h.bspqi[names(h.bspqi.aligned)]
plot(names(fract.aligned), as.numeric(fract.aligned), type='l', col='blue')

h.bsssi.aligned <- table(bsssi.aligned.labelcount$bsssi.bin.ind)
fract.aligned <- h.bsssi.aligned/h.bsssi[names(h.bsssi.aligned)]
lines(names(fract.aligned), as.numeric(fract.aligned), col='red')

h.bspqi.union.aligned <- table(union.aligned.labelcount$bspqi.bin.ind)
h.bsssi.union.aligned <- table(union.aligned.labelcount$bsssi.bin.ind)

#using two-color count
h.twocol.aligned <- table(twocol.aligned.labelcount$twocol.bin.ind)
h.twocol.bspqi.algined <- table(molstats[aligned.bspqi, twocol.bin.ind])
h.twocol.bsssi.aligned <- table(molstats[aligned.bsssi, twocol.bin.ind])
h.twocol.union.aligned <- table(molstats[aligned.union, twocol.bin.ind])

fract.aligned <- h.twocol.aligned/h.twocolor[names(h.twocol.aligned)]
lines(names(fract.aligned), as.numeric(fract.aligned), col='green')

#using bspqi label counts
h.bspqi.twocol.aligned <- table(twocol.aligned.labelcount$bspqi.bin.ind)
fract.aligned <- h.twocol.aligned/h.bspqi[names(h.twocol.aligned)]
lines(names(fract.aligned), as.numeric(fract.aligned), col='purple')

#using bsssI label counts
h.bsssi.twocol.aligned <- table(twocol.aligned.labelcount$bsssi.bin.ind)
fract.aligned <- h.twocol.aligned/h.bsssi[names(h.twocol.aligned)]
lines(names(fract.aligned), as.numeric(fract.aligned), col='orange')

plot(names(h.bspqi), as.numeric(h.bspqi), col='blue', type='l')
lines(names(h.bsssi), as.numeric(h.bsssi), col='red')
lines(names(h.twocolor), as.numeric(h.twocolor), col='green')


plot(names(h.bspqi.twocol.aligned), as.numeric(h.bspqi.twocol.aligned), 
     col='blue', type='l', xlab='Number of BspqI labels', ylab='Number of aligned NGS')
lines(names(h.bspqi.union.aligned), as.numeric(h.bspqi.union.aligned), col='orange')
lines(names(h.bspqi.aligned), as.numeric(h.bspqi.aligned), col='red')

plot(names(h.bsssi.twocol.aligned), as.numeric(h.bsssi.twocol.aligned), 
     col='blue', type='l', xlab='Number of BsssI labels', ylab='Number of aligned NGS')
lines(names(h.bsssi.union.aligned), as.numeric(h.bsssi.union.aligned), col='orange')
lines(names(h.bsssi.aligned), as.numeric(h.bsssi.aligned), col='red')


#plot molecules length vs number of aligned molecules
h.len.two.col.aligned <- table(molstats[aligned.twocolor, length.bin])
h.len.bspqi.aligned <- table(molstats[aligned.bspqi, length.bin])
h.len.bsssi.aligned <- table(molstats[aligned.bsssi, length.bin])
h.len.union.aligned <- table(molstats[aligned.union, length.bin])
h.mollen <- table(molstats[log10.Length > 4, length.bin])

plot(mol.length.bin[as.numeric(names(h.mollen))], as.numeric(h.mollen), 
     xlab='log10(NGS contig length)', ylab='Number of NGS contigs', type='l')
lines(mol.length.bin[as.numeric(names(h.len.two.col.aligned))], as.numeric(h.len.two.col.aligned), col='blue') 
lines(mol.length.bin[as.numeric(names(h.len.bspqi.aligned))], as.numeric(h.len.bspqi.aligned), col='red')
lines(mol.length.bin[as.numeric(names(h.len.bsssi.aligned))], as.numeric(h.len.bsssi.aligned), col='orange')
lines(mol.length.bin[as.numeric(names(h.len.union.aligned))], as.numeric(h.len.union.aligned), col='cyan')
legend(5.4, 450, c("Total", "Two-color aligned", "Union of two enzyme", "BspqI-aligned", "BsssI-aligned"), 
       col=c('black', 'blue', 'cyan', 'red', 'orange'), lwd=1, bty='n', xjust=0, cex=0.8)




#plot molecuels label count vs number of aligned molecules
pdf("~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/Analysis/Label_counts_vs_num_aligned_twocolor_vs_single_zoomIN_updated.pdf")
plot(names(h.bspqi), as.numeric(h.bspqi), 
     xlab='Number of BspqI labels', ylab='Number of NGS', type='l', xlim = c(5, 25), ylim=c(0,1500))
lines(names(h.bspqi.aligned), as.numeric(h.bspqi.aligned), col='red')
lines(names(h.bspqi.twocol.aligned), as.numeric(h.bspqi.twocol.aligned), col='blue')
lines(names(h.bspqi.union.aligned), as.numeric(h.bspqi.union.aligned), col='cyan')
legend(65, 1000, c("Total", "Two-color aligned", "Union of two enzyme", "BspqI-aligned"), 
       col=c('black', 'blue', 'cyan', 'red'), lwd=1, bty='n', xjust=0, cex=0.8)
legend(15, 1000, c("Total", "Two-color aligned", "Union of two enzyme", "BspqI-aligned"), 
       col=c('black', 'blue', 'cyan', 'red'), lwd=1, bty='n', xjust=0, cex=0.8)

plot(names(h.bsssi), as.numeric(h.bsssi), 
     xlab='Number of BssSI labels', ylab='Number of NGS', type='l', xlim=c(5,25), ylim=c(0,1500))
lines(names(h.bsssi.aligned), as.numeric(h.bsssi.aligned), col='orange')
lines(names(h.bsssi.twocol.aligned), as.numeric(h.bsssi.twocol.aligned), col='blue')
lines(names(h.bsssi.union.aligned), as.numeric(h.bsssi.union.aligned), col='cyan')
legend(65, 1000, c("Total", "Two-color aligned", "Union of two enzyme", "BsssI-aligned"), 
       col=c('black', 'blue', 'cyan', 'orange'), lwd=1, bty='n', xjust=0, cex=0.8)
legend(15, 1000, c("Total", "Two-color aligned", "Union of two enzyme", "BsssI-aligned"), 
       col=c('black', 'blue', 'cyan', 'orange'), lwd=1, bty='n', xjust=0, cex=0.8)


plot(names(h.twocolor), as.numeric(h.twocolor), 
     xlab='Number of BspqI + BsssI labels', ylab='Number of NGS', type='l', xlim=c(5,30), ylim=c(0,1500))
lines(names(h.twocol.bspqi.algined), as.numeric(h.twocol.bspqi.algined), col='red')
lines(names(h.twocol.bsssi.aligned), as.numeric(h.twocol.bsssi.aligned), col='orange')
lines(names(h.twocol.aligned), as.numeric(h.twocol.aligned), col='blue')
lines(names(h.twocol.union.aligned), as.numeric(h.twocol.union.aligned), col='cyan')
legend(15, 1000, c("Total", "Two-color aligned", "Union of two enzyme", "BspqI-aligned", "BsssI-aligned"), 
       col=c('black', 'blue', 'cyan', 'red', 'orange'), lwd=1, bty='n', xjust=0, cex=0.8)


dev.off()


#investigate unaligned sequence
molstats2 <- molstats
setnames(molstats2, "CMapId", "QryContigID")

}



