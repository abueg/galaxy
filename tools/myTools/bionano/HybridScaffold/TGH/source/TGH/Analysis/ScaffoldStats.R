#This script analyze and visualize different sandwich-scaffold results
library(ggplot2)

#this plot current benchmark results for TGH on scaffolding BSPQI and BSSSI BNG maps with
#serveral NGS assemblies from NA12878 samples
getTGHScaffoldStats <- function(){
input.ngs.n50 <- c(0.179, 0.554, 0.895)
single.hybrid.n50 <- c(5.029, 7.52, 10.533)
single.ncontigs <- c(893, 596, 479)
single.ngs.aligned <- c(8477, 5498, 3925)
single.ngs.len <- c(2081, 2557, 2665)
sequential.n50 <- c(9.599, 22.365, 35.494)
sequential.ncontigs <- c(725, 215, 201, 90)
sequential.ngs.aligned <- c(8477, 5498, 3925) #double-check
sequential.len <- c(2013.244, 2807.245, 2814.342)
sandwich.hybrid.n50 <- c(18.03, 32.91, 38.127)
sandwich.ncontigs <- c(312, 171, 151)
sandwich.ngs.aligned <- c(11937, 6181, 4387)
sandwich.ngs.len <- c(2319, 2619, 2703)
input.ind <- c(1,2,3)

input.conflicts.ngs1 <- c(4, 20, 2)
input.conflicts.ngs2 <- c(2, 19, 1)
input.conflicts.bng1 <- c(27, 27, 27)
input.conflicts.bng2 <- c(19, 19, 19)
single.conflicts1 <- c(15,10,3)
single.conflicts2 <- c(4,8,6)
sandwich.conflicts1 <- c(6, 5, 3)
sandwich.conflicts2 <- c(3, 7, 6)


#conflicts
df <- data.table(input=input.ind, conflicts=input.conflicts.bng1 + input.conflicts.ngs1, type='Input')
df <- rbind(df, data.table(input=input.ind, conflicts=single.conflicts1, type='Single-enzme'))
df <- rbind(df, data.table(input=input.ind, conflicts=sandwich.conflicts1, type='Two-enzyme'))

df <- data.table(input=input.ind, conflicts=input.conflicts.bng2 + input.conflicts.ngs2, type='Input')
df <- rbind(df, data.table(input=input.ind, conflicts=single.conflicts2, type='Single-enzyme'))
df <- rbind(df, data.table(input=input.ind, conflicts=sandwich.conflicts2, type='Two-enzyme'))

df <- data.table(input=input.ind, conflicts=input.conflicts.bng1 + input.conflicts.ngs1
                 + input.conflicts.bng2 + input.conflicts.ngs2, type='Input (BNG + NGS)')
df <- rbind(df, data.table(input=input.ind, conflicts=single.conflicts1 + single.conflicts2, type='Single-enzyme'))
df <- rbind(df, data.table(input=input.ind, conflicts=sandwich.conflicts1 + sandwich.conflicts2, type='Two-Enzyme'))


g <- ggplot(df)
g <- g + geom_bar(aes(x=input, y=conflicts, fill=type, group=type), stat='identity', position='dodge')
#g <- g + geom_text(aes(x=input, y=scaffold, label=as.character(scaffold)), position=position_dodge(0.9))
g <- g + ylab('Number of chimeric contigs detected') + xlab('')
g <- g + xlim('Illumina-D', 'Illumina-S', 'PacBio')
g <- g + scale_y_continuous(limit=c(0,50))
g <- g+ theme(text=element_text(size=18))
plot(g)


#N50
df <- data.table(input=input.ind, scaffold=input.ngs.n50, Scaffold='Input NGS')
df <- rbind(df, data.table(input=input.ind, scaffold=single.hybrid.n50, Scaffold='Single-enzyme hybrid-scaffold'))
df <- rbind(df, data.table(input=input.ind, scaffold=sequential.n50, Scaffold='tSequential hybrid scaffold'))
df <- rbind(df, data.table(input=input.ind, scaffold=sandwich.hybrid.n50, Scaffold='Two-enzyme hybrid-scaffold'))

g <- ggplot(df)
g <- g + geom_bar(aes(x=input, y=scaffold, fill=Scaffold, group=Scaffold), stat='identity', position='dodge')
#g <- g + geom_text(aes(x=input, y=scaffold, label=as.character(scaffold)), position=position_dodge(0.9))
g <- g + ylab('Hybrid Scaffold N50') + xlab('')
g <- g + xlim('Illumina-D', 'Illumina-S', 'PacBio')
g <- g+ theme(text=element_text(size=18))
plot(g)

#number of scaffold
df <- data.table(input=input.ind, scaffold=c(0,0,0), Scaffold='Input NGS')
df <- rbind(df, data.table(input=input.ind, scaffold=single.ncontigs, Scaffold='Single-enzyme hybrid-scaffold'))
df <- rbind(df, data.table(input=input.ind, scaffold=sequential.ncontigs, Scaffold='tSequential hybrid-scaffold'))
df <- rbind(df, data.table(input=input.ind, scaffold=sandwich.ncontigs, Scaffold='Two-enzyme hybrid-scaffold'))

g <- ggplot(df)
g <- g + geom_bar(aes(x=input, y=scaffold, group=Scaffold, fill=Scaffold), stat='identity', position='dodge')

g <- g + ylab('Number of scaffolds') + xlab('')
g <- g + xlim('Illumina-D', 'Illumina-S', 'PacBio')
g <- g+ theme(text=element_text(size=18))

plot(g)

#NGS anchored
df <- data.table(input=input.ind, scaffold=c(0,0,0), Scaffold='Input NGS')
df <- rbind(df, data.table(input=input.ind, scaffold=single.ngs.aligned, Scaffold='Single-enzyme hybrid-scaffold'))
df <- rbind(df, data.table(input=input.ind, scaffold=sequential.ngs.aligned, Scaffold='tSequential hybrid-scaffold'))
df <- rbind(df, data.table(input=input.ind, scaffold=sandwich.ngs.aligned, Scaffold='Two-enzyme hybrid-scaffold'))

g <- ggplot(df)
g <- g + geom_bar(aes(x=input, y=scaffold, group=Scaffold, fill=Scaffold), stat='identity', position='dodge')
g <- g + ylab('# of NGS anchored in Scaffold') + xlab('')
g <- g + xlim('Illumina-D', 'Illumina-S', 'PacBio')
g <- g+ theme(text=element_text(size=18))

plot(g)

df <- data.table(input=input.ind, scaffold=c(0,0,0), Scaffold='Input NGS')
df <- rbind(df, data.table(input=input.ind, scaffold=single.ngs.len, Scaffold='Single-enzyme hybrid-scaffold'))
df <- rbind(df, data.table(input=input.ind, scaffold=sandwich.ngs.len, Scaffold='Two-enzyme hybrid-scaffold'))

g <- ggplot(df)
g <- g + geom_bar(aes(x=input, y=scaffold, group=Scaffold, fill=Scaffold), stat='identity', position='dodge')
g <- g + ylab('Total NGS length in scaffold') + xlab('')
g <- g + xlim('Illumina-D', 'Illumina-S', 'PacBio')
g <- g+ theme(text=element_text(size=18))

plot(g)

}

getNonHumanTGHStats <- function(){
  file <- '~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/Analysis/TGH_results.csv'
  dt <- as.data.table(read.table(file, sep=',', header=T))
  sbj.dt <- dt[Dataset == 'SugarBeet']
  sbj.dt <- dt[Dataset == 'HummingBird']

  sbj.dt <- sbj.dt[Scaffold_pipeline != 'Sequential' & Scaffold_pipeline !='Single-enzyme BSSSI']
  g <- ggplot(sbj.dt)
  g <- g+ geom_line(aes(x=as.factor(NGS), y=NGS_N50, group="Input", colour="Input"))#, colour = "Input", group="Input"))
  g <- g+ geom_line(aes(x=as.factor(NGS), y=Scaffold_NG50, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  #g <- g+ geom_line(aes(x=as.factor(NGS), y=Overall_Scaffold_N50, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  g <- g+ geom_point(aes(x=as.factor(NGS), y=NGS_N50, colour = "Input"))
  g <- g+ geom_point(aes(x=as.factor(NGS), y=Scaffold_NG50, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  #g <- g+ geom_point(aes(x=as.factor(NGS), y=Overall_Scaffold_N50, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  g <- g + xlab("NGS data") + ylab("Scaffold NG50 (Mb)")
  plot(g)

  g <- ggplot(sbj.dt)
  g <- g+ geom_line(aes(x=as.factor(NGS), y=Genome_Size, colour = "Input", group="Input"))
  g <- g+ geom_point(aes(x=as.factor(NGS), y=Genome_Size, colour = "Input", group="Input"))
  g <- g+ geom_line(aes(x=as.factor(NGS), y=Scaffold_Size, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  g <- g+ geom_point(aes(x=as.factor(NGS), y=Scaffold_Size, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  #g <- g+ scale_y_continuous(limit=c(0,600))
  plot(g)



  g <- ggplot(sbj.dt)
  g <- g+ geom_line(aes(x=as.factor(NGS), y=N_NGS_filtered, group="Input", colour="Input"))#, colour = "Input", group="Input"))
  g <- g+ geom_line(aes(x=as.factor(NGS), y=N_NGS_Aligned, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  g <- g+ geom_point(aes(x=as.factor(NGS), y=N_NGS_filtered, group="Input", colour = "Input"))
  g <- g+ geom_point(aes(x=as.factor(NGS), y=N_NGS_Aligned, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  plot(g)

  g <- ggplot(sbj.dt)
  g <- g+ geom_line(aes(x=as.factor(NGS), y=NGS_Size, colour = "Input", group="Input"))
  g <- g+ geom_point(aes(x=as.factor(NGS), y=NGS_Size, colour = "Input", group="Input"))
  g <- g+ geom_line(aes(x=as.factor(NGS), y=Len_Aligned, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  g <- g+ geom_point(aes(x=as.factor(NGS), y=Len_Aligned, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  #g <- g+ scale_y_continuous(limit=c(0,600))
  plot(g)

  g <- ggplot(sbj.dt)
  g <- g+ geom_line(aes(x=as.factor(NGS), y=Scaffold_order_accuracy, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  g <- g+ geom_point(aes(x=as.factor(NGS), y=Scaffold_order_accuracy, colour = Scaffold_pipeline, group=Scaffold_pipeline))
  g <- g+ scale_y_continuous(limit=c(0.5,1))
  plot(g)
}


#for some scaffold the total size is different in this case we calculate NG50 also
getNG50 <- function(){
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/SingleEnzyme/PB30x/hybrid_scaffolds/"
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/PB30x/BSSSI/hybrid_scaffolds_M1/"
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/SingleEnzyme/PB30x/Step2/hybrid_scaffolds/"

  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/SingleEnzyme/PB40x/hybrid_scaffolds/"
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/PB40x/BSSSI/hybrid_scaffolds_M1/"
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/SingleEnzyme/PB40x/Step2/hybrid_scaffolds/"

  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/SingleEnzyme/PB60x/hybrid_scaffolds/"
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/PB60x/BSSSI/hybrid_scaffolds_M1/"
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/SingleEnzyme/PB60x/Step2/hybrid_scaffolds/"

  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSoyBean/ReRun12122016/SingleEnzyme/hybrid_scaffolds/"
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSoyBean/ReRun12122016/BSSSI/hybrid_scaffolds_M1/"
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSoyBean/ReRun12122016/SingleEnzyme/Step2/hybrid_scaffolds/"

  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSoyBean/NewSeqAsm/SingleEnzyme/hybrid_scaffolds/"
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSoyBean/NewSeqAsm/BSSSI/hybrid_scaffolds_M1/"
  scaffold.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSoyBean/NewSeqAsm/SingleEnzyme/Step2/hybrid_scaffolds/"

  scaffold.dir <- "~/RemoteServer/home/users2/ahastie/data/HybridScaffold_twoEnzyme/Duke/20151212_autocutDuke_hummingbird/20151210BspQI_ILMN_N2B2/hybrid_scaffolds/"
  scaffold.dir <- "~/RemoteServer/home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/PB_Andy/BSSSI/hybrid_scaffolds_M1/"
  scaffold.dir <- "~/RemoteServer/home/users2/ahastie/data/HybridScaffold_twoEnzyme/Duke/20151212_autocutDuke_hummingbird/20151210BssSI_BspQIscaf_ILMN_N2B2//hybrid_scaffolds/"

  scaffold.dir <- "~/RemoteServer/home/users2/ahastie/data/HybridScaffold_twoEnzyme/Duke/20151212_autocutDuke_hummingbird/20151210BssSI_BspQIscaf_N2B2/hybrid_scaffolds/"

  hybrid.file <- dir(scaffold.dir, pattern = "*HYBRID_SCAFFOLD.cmap", full.names =T)
  #hybrid.file <- hybrid.file[[2]]
  single.scaff1 <- as.data.table(readcmap(hybrid.file))
  getN50(single.scaff1[SiteID==1, ContigLength])
  getNG50(single.scaff1[SiteID==1, ContigLength], 1076e6)

  tgh.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/PB30x/TGH_M1/Sandwich2/"
  tgh.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/PB40x/TGH_M1/Sandwich2/"
  tgh.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/PB60x/TGH_M1/Sandwich2/"
  tgh.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSBJ_12092016/FullPB/TGH_M1/Sandwich2/"

  tgh.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSoyBean/ReRun12122016/TGH_M1/Sandwich2/"
  tgh.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestRunSoyBean/NewSeqAsm/TGH_M1/Sandwich2/"
  tgh.dir <- "~/RemoteServer//home/users/jwang/workspace/HybridScaffold/HybridScaffodTwo/TestHummingBird_Auto/ReRun_12082016//TGH_M1/Sandwich2/"

  hybrid.file.list <- dir(tgh.dir, "test_2col_c*", full.names = T)
  combined.cmap <- get.merged.cmap(hybrid.file.list)
  getN50(combined.cmap[SiteID==1, ContigLength])
  getNG50(combined.cmap[SiteID==1, ContigLength], 1076e6)

}


