#!/usr/bin/R
rm(list=ls())
library(gplots)
library(ggplot2)
library(snapCGH)
library(cluster)
library(limma)
library(Vennerable)
library(DESeq)
library(Ringo)
library(tools)
library(plyr)
library(gridExtra)
library(Biostrings)
library(seqLogo)
library(ggseqlogo)
library(RColorBrewer)
source(file.path("/home/anton/data/R-script/R-counts/RUN", "functions.R"))


wde <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/repeat/sample_S1_L001_R1_001_expr/Dump"
wdn <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/repeat/sample_S1_L001_R1_001_norm/Dump"
wdm <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/repeat/sample_S1_L001_R1_001_map/Dump"
#mwd95 <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_BC_Mut_mapping/Dump"
#mwd80 <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_BC_Mut_mapping_80/Dump"
popName <- data.frame("pop1" = c("e18-1-1", "e18-1-2", "n18-1-1", "n18-1-2"), "pop2" = c("e18-2-1", "e18-2-2", "n18-2-1", "n18-2-2"))
mapName <- data.frame("pop1" = c("m18-1-1", "m18-1-2"), "pop2" = c("m18-2-1", "m18-2-2"))

pmi <- c('AGCTC', 'ACGTA', 'CTGCT', 'AGTCA', 'TCAAA', 'TTGAG', 'TCGCT')
names(pmi) <- c("HexA", "Hsp70", "MtnA", "PCNA", "Pyk", "Tbp", "Promoterless")
output <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/repeat/Stat_output"
#setwd(wd)

ShowData <- function() {
	for (i in names(DATA)){print(i);print(head(DATA[[i]]))}
}
RemoveUnvData <- function(df, col1, col2) {
df <- df[!is.na(df[,col1]) & !is.na(df[,col2]),]
df <- df[!is.infinite(df[,col1]) & !is.infinite(df[,col2]),]
return(df)
}

DATA <- list()
for (src in c(wde, wdn)){
  loadFiles <- list.files(path=src, pattern=".*([ATGC]{5}).txt", full.names = T)
  if (exists("frcData")) {
    frcData <- rbind(frcData, read.delim(list.files(path=src, pattern="frc.*.txt", full.names = T)))
  } else {
    frcData <- read.delim(list.files(path=src, pattern="frc.*.txt", full.names = T))  
  }
  for (f in loadFiles){
    fileName <- file_path_sans_ext(basename(f))
    exp <- as.character(strsplit(fileName, "_")[[1]][1])
    subIndex <- strsplit(fileName, "_")[[1]][2]
    popnm <- paste0("pop", sub(".?[0-9]{2}-([0-9]?)-[0-9]?", "\\1", exp))
    if (!(exists(popnm, where = DATA))){
      DATA[[popnm]] <- list()
    }
    if (!(exists(exp, where = DATA[[popnm]]))) {
      DATA[[popnm]][[exp]] <- read.delim(f, stringsAsFactors = F)
      print(paste(popnm, exp, nrow(DATA[[popnm]][[exp]]), "init"))
      next
    }
    DATA[[popnm]][[exp]] <- rbind(DATA[[popnm]][[exp]], read.delim(f, stringsAsFactors = F))
    print(paste(popnm, exp, nrow(DATA[[popnm]][[exp]]), "merge"))
  }
}
for (popnm in names(DATA)){
  for (exp in names(DATA[[popnm]])){
    DATA[[popnm]][[exp]] <- aggregate(value ~ id, data = DATA[[popnm]][[exp]], sum)
    print(paste(popnm, exp, "aggregate", nrow(DATA[[popnm]][[exp]])))
    DATA[[popnm]][[exp]] <- DATA[[popnm]][[exp]][order(DATA[[popnm]][[exp]]$id),]
    DATA[[popnm]][[exp]]$rpm <- (DATA[[popnm]][[exp]]$value / frcData$value[frcData$id == exp]) * 10^6
  }
}

# Merge data between tech. replicates
# mergedMap <- list()
# for (p in names(MAP)){
#     mergedMap[[p]] <- list()
#     for (mpp in names(MAP[[p]])){
#         if (strsplit(mpp, "-")[[1]][3] == 1) {
#             mpp2 <- sub(".$", "2", mpp)
#             for (subIndex in pmi){
#                 mergedMap[[p]][[subIndex]] <- unique(c(MAP[[p]][[mpp]][[subIndex]]$id, MAP[[p]][[mpp2]][[subIndex]]$id))
#             }
#         }
#     }
# }


for (pop in names(DATA)) {
  for (exp in names(DATA[[pop]])){
    if (substr(exp, 1, 1) == "e") {
      norm <- sub("^.", "n", exp)
      DATA[[pop]][[exp]]$rpm_norm <- DATA[[pop]][[norm]]$rpm
      DATA[[pop]][[exp]]$norm <- log2(DATA[[pop]][[exp]]$rpm / DATA[[pop]][[exp]]$rpm_norm)
    }
  }
}

scList <- list()
for (pop in names(DATA)) {
  scList[[pop]] <- list()
  for (exp in names(DATA[[pop]])) {
    if (substr(exp, 1, 1) == "e" & substr(exp, nchar(exp), nchar(exp)) == 1) {
      exp2 <- sub(".$", "2", exp)
      e1df <- DATA[[pop]][[exp]]
      names(e1df)[-1] <- paste(names(e1df)[-1], "e1", sep="_")
      e2df <- DATA[[pop]][[exp2]]
      names(e2df)[-1] <- paste(names(e2df)[-1], "e2", sep="_")
      scList[[pop]] <- cbind(e1df, e2df[, -1])
      scList[[pop]]$mean <- rowMeans(cbind(scList[[pop]]$norm_e1, scList[[pop]]$norm_e2), na.rm=T)
      for (i in 2:(ncol(scList[[pop]]))) {
          nan.index <- is.nan(scList[[pop]][, i])
          inf.index <- is.infinite(scList[[pop]][, i])
          scList[[pop]][nan.index, i] <- NA
          scList[[pop]][inf.index, i] <- NA
          rm(nan.index, inf.index)
      }
    }
  }
}


# Filter scList
# scList_filt <- list()
# for (p in names(scList)) {
#     scList_filt[[p]] <- list()
#     for (subIndex in names(scList[[p]])){
#         scList_filt[[p]][[subIndex]] <- scList[[p]][[subIndex]][ scList[[p]][[subIndex]]$id %in% mergedMap[[p]][[subIndex]], ]
#         print(paste(p, subIndex, "original scList:", nrow(scList[[p]][[subIndex]]), "filtered scList:", nrow(scList_filt[[p]][[subIndex]]) ))
#     }
# }
# scList <- scList_filt


plotList <- list()
for (pop in names(scList)){
    Cor.P <- paste("Pearson.Cor = ", round(with(scList[[pop]], cor(norm_e1, norm_e2, method="pearson", use="pairwise.complete.obs")), digits=2), sep = "")
    Cor.S <- paste("Spearman.Cor = ", round(with(scList[[pop]], cor(norm_e1, norm_e2, method="spearman", use="pairwise.complete.obs")), digits=2), sep = "")
    allValues <- paste("Values with NA = ", nrow(scList[[pop]]), sep = "")
    x.coord <- rep(min(scList[[pop]]$norm_e1, na.rm = T) + 1.5, 3)
    byVal <- -((max(scList[[pop]]$norm_e2, na.rm = T) / 100) * 15)
    y.coord <- seq(max(scList[[pop]]$norm_e2, na.rm = T) * 0.7, by=byVal, length.out=3)
    plotList$scplot <- ggplot(scList[[pop]], aes(norm_e1, norm_e2)) + geom_point(alpha = 5/10, colour = "#65A6D1", size = 3, na.rm = T) + xlab("Tech.norm.exp. 1 (log2)") + ylab("Tech.norm.exp. 2 (log2)") + annotate("text", x = x.coord, y = y.coord, label = c(Cor.P, Cor.S, allValues), size = 5, na.rm=T) + labs(title = paste(pop, "scatter plot tech. expression 1&2")) + theme_bw() + theme(plot.title = element_text(size = 18), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))

    Cor.P0 <- paste("Pearson.Cor = ", round(with(scList[[pop]], cor(rpm_norm_e1, rpm_norm_e2, method="pearson", use="pairwise.complete.obs")), digits=2), sep = "")
    Cor.S0 <- paste("Spearman.Cor = ", round(with(scList[[pop]], cor(rpm_norm_e1, rpm_norm_e2, method="spearman", use="pairwise.complete.obs")), digits=2), sep = "")
    allValues0 <- paste("Values with NA = ", nrow(scList[[pop]]), sep = "")
    x.coord0 <- rep(min(scList[[pop]]$rpm_norm_e1, na.rm = T) + 1.5, 3)
    byVal0 <- -((max(scList[[pop]]$rpm_norm_e2, na.rm = T) / 100) * 15)
    y.coord0 <- seq(max(scList[[pop]]$rpm_norm_e2, na.rm = T) * 0.7, by=byVal0, length.out=3)
    plotList$scplot0 <- ggplot(scList[[pop]], aes(rpm_norm_e1, rpm_norm_e2)) + geom_point(alpha = 5/10, colour = "#65A6D1", size = 3, na.rm = T) + xlab("Tech.rmp.norm. 1 (log2)") + ylab("Tech.rmp.norm. 2 (log2)") + annotate("text", x = x.coord0, y = y.coord0, label = c(Cor.P0, Cor.S, allValues), size = 5, na.rm=T) + labs(title = paste(pop, "scatter plot tech. normalization on RPM 1&2")) + theme_bw() + theme(plot.title = element_text(size = 18), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
    pdf(file=file.path(output, paste(pop, "Scatter plot by expression and normalization data.pdf")), width=14, height=14)
    print(multiplot(plotlist=plotList, cols=2, layout=matrix(c(1:2), ncol = 1, nrow = 2, byrow=T)))
    dev.off()
}


plotList <- list()
for (pop in names(scList)){
  for (subIndex in names(scList[[pop]])) {
    colourVector <- brewer.pal(3, 'Accent')
    names(colourVector) <- c("norm_e1", "norm_e2", "mean")
    piName <- names(pmi)[pmi ==  subIndex]
    Cor.P <- paste("Pearson.Cor = ", round(with(scList[[pop]][[subIndex]], cor(norm_e1, norm_e2, method="pearson", use="pairwise.complete.obs")), digits=2), sep = "")
    Cor.S <- paste("Spearman.Cor = ", round(with(scList[[pop]][[subIndex]], cor(norm_e1, norm_e2, method="spearman", use="pairwise.complete.obs")), digits=2), sep = "")
    allValues <- paste("Values with NA = ", nrow(scList[[pop]][[subIndex]]), sep = "")
    x.coord <- rep(min(scList[[pop]][[subIndex]]$norm_e1, na.rm = T) + 1.5, 3)
    byVal <- -((max(scList[[pop]][[subIndex]]$norm_e2, na.rm = T) / 100) * 15)
    y.coord <- seq(max(scList[[pop]][[subIndex]]$norm_e2, na.rm = T) * 0.7, by=byVal, length.out=3)

    # Density plot for expression replicate 1/2
    dsplotE1E2_main <- ggplot() + geom_density(colour=colourVector['norm_e1'], data=scList[[pop]][[subIndex]], aes(norm_e1)) + geom_density(colour=colourVector['norm_e2'], data=scList[[pop]][[subIndex]], aes(norm_e2)) + theme_bw() + xlab("Norm.exp. log2 value") + ylab("Expression count")
    xee <- rep(-3, 2)
    byee <- -((max(scList[[pop]][[subIndex]]$norm_e1, na.rm = T) / 100) * 5)
    yee <- seq(max(scList[[pop]][[subIndex]]$norm_e1, na.rm=T) * 0.7, by=byee, length.out=2)
    cMean <- c(mean(scList[[pop]][[subIndex]]$norm_e1, na.rm=T), mean(scList[[pop]][[subIndex]]$norm_e2, na.rm=T))
    plotList$dsplotE1E2_main <- dsplotE1E2_main + labs(title = paste(subIndex, piName, pop, "Tech. expression replicates 1&2")) + annotate("text", x = xee, y = yee, label = c(paste("Total barcodes expr#1:", length(scList[[pop]][[subIndex]]$norm_e1[!is.na(scList[[pop]][[subIndex]]$norm_e1)])), paste("Total barcodes expr#2:", length(scList[[pop]][[subIndex]]$norm_e2[!is.na(scList[[pop]][[subIndex]]$norm_e2)])))) + geom_vline(aes(xintercept=cMean), colour=c(colourVector[1], colourVector[2]))

    # Density plot for expression replicate mean
    xm <- 0
    bym <- -((max(scList[[pop]][[subIndex]]$mean, na.rm = T) / 100) * 5)
    ym <- seq(max(scList[[pop]][[subIndex]]$mean, na.rm=T) * 0.7, by=byee, length.out=1)
    dsplotEMean_main <- ggplot() + geom_density(colour=colourVector['mean'], data=scList[[pop]][[subIndex]], aes(mean)) + theme_bw() + xlab("Norm. expression log2 value") + ylab("Expression log2 count")
    plotList$dsplotEMean_main <- dsplotEMean_main + annotate("text", x = xm, y = ym, label = paste("Total barcodes expr. mean:", length(scList[[pop]][[subIndex]]$mean[!is.na(scList[[pop]][[subIndex]]$mean)]))) + labs(title = paste(subIndex, piName, pop, "Mean expression"))
    plotList$scplot <- ggplot(scList[[pop]][[subIndex]], aes(norm_e1, norm_e2)) + geom_point(alpha = 5/10, colour = "#65A6D1", size = 3, na.rm = T) + xlab("Tech.norm.exp. 1 (log2)") + ylab("Tech.norm.exp. 2 (log2)") + annotate("text", x = x.coord, y = y.coord, label = c(Cor.P, Cor.S, allValues), size = 5, na.rm=T) + labs(title = paste(subIndex, piName, pop, "scatter plot tech. expression 1&2")) + theme_bw() + theme(plot.title = element_text(size = 18), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
    
    pdf(file=file.path(output, paste(subIndex, piName, pop, "Density plot by expression and mapping data.pdf")), width=14, height=14)
    print(multiplot(plotlist=plotList, cols=2, layout=matrix(c(1:4), ncol = 2, nrow = 2, byrow=T)))
    dev.off()
  }
}

plotListComp <- list()
for (pop in names(scList)){
  colorComp <- brewer.pal(length(pmi), 'Accent')
  names(colorComp) <- pmi
  if (exists("dsplotEMean_comp")) {rm(dsplotEMean_comp)}
  for (subIndex in pmi) {
    if (exists("dsplotEMean_comp")) {
      dsplotEMean_comp <- dsplotEMean_comp + geom_density(colour=colorComp[[subIndex]], data=scList[[pop]][[subIndex]], aes(mean))
      count <- c(count, length(scList[[pop]][[subIndex]]$mean[!is.na(scList[[pop]][[subIndex]]$mean)]))
      cMean <- c(cMean, round(mean(scList[[pop]][[subIndex]]$mean, na.rm=T), 2))
      Prange <- rbind(Prange, data.frame("min" = round(range(scList[[pop]][[subIndex]]$mean, na.rm=T)[1], 1), "max" = round(range(scList[[pop]][[subIndex]]$mean, na.rm=T)[2], 1)))
    } else {
      dsplotEMean_comp <- ggplot() + geom_density(colour=colorComp[[subIndex]], data=scList[[pop]][[subIndex]], aes(mean))
      count <- length(scList[[pop]][[subIndex]]$mean[!is.na(scList[[pop]][[subIndex]]$mean)])
      cMean <- round(mean(scList[[pop]][[subIndex]]$mean, na.rm=T), 2)
      Prange <- data.frame("min" = round(range(scList[[pop]][[subIndex]]$mean, na.rm=T)[1], 1), "max" = round(range(scList[[pop]][[subIndex]]$mean, na.rm=T)[2], 1))
    }
  }
  xa <- rep(-5, length(pmi))
    ya <- seq(1, by=-0.05, length.out=length(pmi))

    xaa <- rep(-4, length(pmi))
    xaae <- rep(-3.5, length(pmi))

    count <- paste("n:", count)
    xac <- rep(-3, length(pmi))

    xar <- rep(-2, length(pmi))
    Prange <- paste("r:", Prange$min, ":", Prange$max)

    xm <- rep(-0.5, length(pmi))
    cMeantxt <- paste("mean:", cMean)

  plotListComp[[pop]] <- dsplotEMean_comp + theme_bw() + xlab("Norm. expression log2 value") + ylab("Expression count") + labs(title = paste(pop, "Mean expression")) + annotate("text", x = xa, y = ya, label = names(pmi)) + annotate("segment", x = xaa, xend = xaae, y = ya, yend = ya, colour = colorComp) + annotate("text", x = xac, y = ya, label = count) + annotate("text", x = xar, y = ya, label = Prange) + geom_vline(aes(xintercept=cMean), colour=colorComp) + annotate("text", x = xm, y = ya, label = cMeantxt)
}
pdf(file=file.path(output, paste("Summary density plot by expression in all promotor index.pdf")), width=14, height=14)
print(multiplot(plotlist=plotListComp, cols=1))
dev.off()

plotList <- list()
# Density plot for expression replicate 1/2
dsplotE1E2_main <- ggplot() + geom_density(colour="#000000", data=scDF, aes(e1)) + geom_vline(aes(xintercept=ControlDF$e1_norm), colour=c("#99ff00", "#99cc33", "#00ffff", "#006666")) + geom_density(colour="#8c762b", data=scDF, aes(e2)) + geom_vline(aes(xintercept=ControlDF$e2_norm), colour=c("#ff3366", "#990000", "#3333cc", "#9966ff")) + theme_bw() + xlab("Norm. expression value") + ylab("Expression count")
plotList$dsplotE1E2_main <- dsplotE1E2_main + annotate("text", x = c(6, 6), y = c(0.6, 0.5), label = c(paste("Total barcodes expr#1:", nrow(DATA$e1)), paste("Total barcodes expr#2:", nrow(DATA$e2)))) + labs(title = "Individual expression replicates #1/#2")
plotList$dsplotE1E2_trim <- plotList$dsplotE1E2_main + xlim(c(0, 3)) + annotate("text", x = c(1.7, 1.7), y = c(0.6, 0.5), label = c(paste("Total barcodes expr#1:", nrow(DATA$e1)), paste("Total barcodes expr#2:", nrow(DATA$e2)))) + labs(title = "Individual expression replicates #1/#2. Trimmed")

# Density plot for expression replicate mean
dsplotEMean_main <- ggplot() + geom_density(colour="#ff9933", data=scDF, aes(mean)) + geom_vline(aes(xintercept=cMean), colour=c("#336600", "#330066")) + theme_bw() + xlab("Norm. expression value") + ylab("Expression count")

plotList$dsplotEMean_main <- dsplotEMean_main + annotate("text", x = 6, y = c(0.6), label = paste("Total barcodes expr. mean:", nrow(scDF))) + labs(title = "Mean expression")

plotList$dsplotEMean_trim <- plotList$dsplotEMean_main + xlim(c(0, 3)) + annotate("text", x = 1.7, y = c(0.6), label = paste("Total barcodes expr. mean:", nrow(scDF))) + labs(title = "Mean expression. Trimmed")

pdf(file=file.path(output, "Density plot by expression and mapping data.pdf"), width=14, height=14)
print(multiplot(plotlist=plotList, cols=2, layout=matrix(c(1:12), ncol = 2, nrow = 6, byrow=T)))
dev.off()


stringList <- list()
for (i in c("result_mp95_by_mut", "result_mp80_by_mut")) {
	iName <- sub(".*([0-9]{2}).*", "prob_\\1", i)
	vmin <- cMean[1]
	vmax <- cMean[2]
	tmpDF <- mappingList[[i]][order(mappingList[[i]]$mean),]
	tmpDFmin <- tmpDF[tmpDF$mean < vmin,]
	tmpDFmax <- tmpDF[tmpDF$mean > vmax,]
	min <- paste0(iName, "_MIN")
	max <- paste0(iName, "_MAX")
	control <- paste0(iName, "_CONTROL")
	fstpct <- paste0(iName, "_1-st_pct")
	lstpct <- paste0(iName, "_99-th_pct")
	controlpct <- paste0(iName, "_2-nd_to_98-th_pct_CONTROL")
	stringList[[min]] <- tmpDFmin$mut[c(1:round(length(tmpDFmin$mut)*0.1))]
	stringList[[max]] <- tmpDFmax$mut[c(round(length(tmpDFmax$mut)*0.9):length(tmpDFmax$mut))]
	stringList[[control]] <- tmpDF$mut[c(round(length(tmpDFmin$mut)*0.1)+1:which(tmpDF == tmpDFmax$mean[round(length(tmpDFmax$mut)*0.9)], arr.ind=T)[1])]
	stringList[[fstpct]] <- tmpDF$mut[c(1:round(length(tmpDF$mut) * 0.01))]
	stringList[[lstpct]] <- tmpDF$mut[c( round( length( tmpDF$mut ) * 0.99 ) : length(tmpDF$mut) )]
	stringList[[controlpct]] <- tmpDF$mut[c( round(length(tmpDF$mut) * 0.01) + 1 : round( length( tmpDF$mut ) * 0.99 ) - 1)]
}
stringPlotList <- list()
for (ss in names(stringList)){
	for (m in c("bits", "prob")){
		name <- paste0(ss, "_", m)
		stringPlotList[[name]] <- ggplot() + geom_logo(as.character(stringList[[ss]]), method=m, seq_type="dna") + theme_logo() + labs(title=paste0(ss, " by ", m, ". Count: ", length(stringList[[ss]])))
	}
}
pdf(file=file.path(output, "Position weight matrix_95.pdf"), width=14, height=14)
print(multiplot(plotlist=stringPlotList[1:12], cols=2, layout=matrix(c(1:12), ncol = 2, nrow = 6, byrow=T)))
dev.off()
pdf(file=file.path(output, "Position weight matrix_80.pdf"), width=14, height=14)
print(multiplot(plotlist=stringPlotList[13:24], cols=2, layout=matrix(c(1:12), ncol = 2, nrow = 6, byrow=T)))
dev.off()
