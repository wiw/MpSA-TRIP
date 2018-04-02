#!/usr/bin/R
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


rm(list=ls())
wde <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_Genome_Expr_A20-23/Dump"
wdn <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_Genome_Norm_A10-13/Dump"
#mwd95 <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_BC_Mut_mapping/Dump"
#mwd80 <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_BC_Mut_mapping_80/Dump"
popName <- data.frame("pop1" = c("e18-1-1", "e18-1-2", "n18-1-1", "n18-1-2"), "pop2" = c("e18-2-1", "e18-2-2", "n18-2-1", "n18-2-2"))
pmi <- c('AGCTC', 'ACGTA', 'CTGCT', 'AGTCA', 'TTGAG', 'TCGCT')
names(pmi) <- c("HexA", "Hsp70", "MtnA", "PCNA", "Tbp", "Promoterless")
output <- "/home/anton/backup/input/trip/RUN_2017-11-27/results/Stat_output/genome"
#setwd(wd)

ShowData <- function() {
    for (i in names(DATA)){print(i);print(head(DATA[[i]]))}
}
RemoveUnvData <- function(df, col1, col2) {
df <- df[!is.na(df[,col1]) & !is.na(df[,col2]),]
df <- df[!is.infinite(df[,col1]) & !is.infinite(df[,col2]),]
return(df)
}
RemoveUnvDataS <- function(df, col1) {
df <- df[!is.na(df[,col1]),]
df <- df[!is.infinite(df[,col1]),]
return(df)
}

# Load count for PI in replicate
countPIFiles <- list.files(path=wde, pattern=".*_count.txt", full.names=T)
if (exists("countDF")){rm(countDF)}
for (df in countPIFiles) {
    rName <- sub("(.*)_count\\.txt", "\\1", basename(df))
    if (exists("countDF")) {
        tdf <- read.delim(df, stringsAsFactors=F)
        countDF <- cbind(countDF, tdf[, 2])
        names(countDF)[ncol(countDF)] <- rName
        rm(tdf)
        } else {
            countDF <- read.delim(df, stringsAsFactors=F)
            names(countDF) <- c("PI", rName)
        }
}
# Load frc data
for (src in c(wde, wdn)){
  if (exists("frcData")) {
    frcData <- rbind(frcData, read.delim(list.files(path=src, pattern="frc.*.txt", full.names = T)))
  } else {
    frcData <- read.delim(list.files(path=src, pattern="frc.*.txt", full.names = T))  
  }
}
frcData <- data.frame("id" = c("e18-1", "e18-2", "n18-1", "n18-2"), "value"= c(sum(frcData$value[1:2]), sum(frcData$value[3:4]), sum(frcData$value[5:6]), sum(frcData$value[7:8])))
# Load main data
DATA <- list()
for (src in c(wde, wdn)){
  loadFiles <- list.files(path=src, pattern=".*([ATGC]{5}).txt", full.names = T)
  for (f in loadFiles){
    fileName <- file_path_sans_ext(basename(f))
    exp <- as.character(strsplit(fileName, "_")[[1]][1])
    pi <- strsplit(fileName, "_")[[1]][2]
    if (is.list(DATA[[exp]])) {
        DATA[[exp]][[pi]] <- read.delim(f, stringsAsFactors = F)
        DATA[[exp]][[pi]] <- DATA[[exp]][[pi]][order(DATA[[exp]][[pi]]$id),]
        DATA[[exp]][[pi]]$rpm <- (DATA[[exp]][[pi]]$value / frcData$value[frcData$id == exp]) * 10^6
    } else {
        DATA[[exp]] <- list()
        DATA[[exp]][[pi]] <- read.delim(f, stringsAsFactors = F)
        DATA[[exp]][[pi]] <- DATA[[exp]][[pi]][order(DATA[[exp]][[pi]]$id),]
        DATA[[exp]][[pi]]$rpm <- (DATA[[exp]][[pi]]$value / frcData$value[frcData$id == exp]) * 10^6
    }    
  }
}

# loadFiles <- list.files(pattern="*.txt")
# ControlDF <- data.frame("control" = c("wt_bc1", "wt_bc2", "dC_bc3", "dC_bc4"), "n" = c(1633+27, 1654+9, 1641+11, 1632+10), "e1" = c(487+13, 590, 1073+9, 1992+24), "e2" = c(419+11, 503, 916, 1760+9))
# frcData <- read.delim("frc.txt")
# frcData$id <- c("e2", "e1", "n")
# loadFiles <- loadFiles[loadFiles != "frc.txt"]
# DATA <- list()
# # prepare experiment
# for (i in loadFiles){
#   varName <- sub("(.*)\\.txt", "\\1", i)
#   # read from file
#   DATA[[varName]] <- read.delim(i, stringsAsFactors=F)
#   # sorting
#   DATA[[varName]] <- DATA[[varName]][order(DATA[[varName]]$id), ]
#   # rpm
#   DATA[[varName]]$rpm <- (DATA[[varName]]$value / frcData$value[frcData$id == varName]) * 10^6
# }

# prepare control
# for (i in frcData$id){
#   rpm <- paste(i, "rpm", sep="_")
#   norm <- paste(i, "norm", sep="_")
#   ControlDF[[rpm]] <- (ControlDF[[i]] / frcData$value[frcData$id == i]) * 10^6
# }
# normalization experiments&control
# for (expr in c("e1", "e2")){
#   DATA[[expr]]$norm <- DATA[[expr]]$rpm / DATA$n$rpm
  # rpm <- paste(expr, "rpm", sep="_")
  # norm <- paste(expr, "norm", sep="_")
  # ControlDF[[norm]] <- ControlDF[[rpm]] / ControlDF$n_rpm
# }

# Normalization
for (exp in names(DATA)){
    if (substr(exp, 1, 1) == "e") {
        norm <- sub("^.", "n", exp)
        for (pi in names(DATA[[exp]])) {
            DATA[[exp]][[pi]]$norm <- DATA[[exp]][[pi]]$rpm / DATA[[norm]][[pi]]$rpm 
        }
    }
}

# ControlDF$eMean <- rowMeans(cbind(ControlDF$e1_norm, ControlDF$e2_norm))
# cMean <- rowMeans(cbind(ControlDF$eMean[c(1,3)], ControlDF$eMean[c(2,4)]))
# scList <- list()
# for (exp in names(DATA)) {
#     if (substr(exp, 1, 1) == "e") {
#         exp2 <- sub(".$", "2", exp)
#         for (pi in names(DATA[[exp]])){
#             scList[[pi]] <- data.frame("id" = DATA[[exp]][[pi]]$id, "e1" = DATA[[exp]][[pi]]$norm, "e2" = DATA[[exp2]][[pi]]$norm)
#             scList[[pi]] <- RemoveUnvData(scList[[pi]], 2, 3)
#             scList[[pi]]$mean <- rowMeans(cbind(scList[[pi]]$e1, scList[[pi]]$e2))
#         }
#     }
# }

popList <- list()
for (pop in c("e18-1", "e18-2")){
    popList[[pop]] <- DATA[[pop]]
    for (pi in names(popList[[pop]])) {
        popList[[pop]][[pi]] <- RemoveUnvDataS(popList[[pop]][[pi]], 4)
    }
}

# Select data into one dataframe for drawing pictures
# scDF <- data.frame("id" = DATA$e1$id, "e1" = DATA$e1$norm, "e2" = DATA$e2$norm)
# scDF <- RemoveUnvData(scDF, 2, 3)
# scDF$mean <- rowMeans(cbind(scDF$e1, scDF$e2))
# mappingList <- list()
# 
# # load mapping reads 95/80 CutOff
# mlf95 <- list.files(path=mwd95, pattern="m.*\\.txt", full.names=T)
# mlf80 <- list.files(path=mwd80, pattern="m.*\\.txt", full.names=T)
# 
# for (i in c(mlf95, mlf80)) {
#   varName <- sub("(.*)\\.txt", "\\1", basename(i))
#   mappingResult <- paste0("result_", varName)
#   mappingResultUnq <- paste0(mappingResult, "_by_mut")
#   DATA[[varName]] <- read.delim(i, stringsAsFactors=F)
#   names(DATA[[varName]]) <- c("bc", "mut")
#   mappingList[[mappingResult]] <- merge(scDF, DATA[[varName]], by.x = "id", by.y = "bc")
#   mappingList[[mappingResultUnq]] <- data.frame("mean" = mappingList[[mappingResult]]$mean, "mut" = as.factor(mappingList[[mappingResult]]$mut))
#   mappingList[[mappingResultUnq]] <- aggregate(.~mut, mappingList[[mappingResultUnq]], mean)
# }

# for (pi in names(scList)) {
#     plotList <- list()
#     synPI <- names(pmi)[pmi ==  pi]

#     # Density plot for expression replicate 1/2
#     dsplotE1E2_main <- ggplot() + geom_density(colour="#000000", data=scList[[pi]], aes(e1)) + geom_density(colour="#8c762b", data=scList[[pi]], aes(e2)) + theme_bw() + xlab("Norm. expression value") + ylab("Expression count")
#     # plotList$dsplotE1E2_main <- dsplotE1E2_main + annotate("text", x = c(6, 6), y = c(0.6, 0.5), label = c(paste("Total barcodes expr#1:", nrow(DATA$e1)), paste("Total barcodes expr#2:", nrow(DATA$e2)))) + labs(title = "Individual expression replicates #1/#2")
#     # plotList$dsplotE1E2_trim <- plotList$dsplotE1E2_main + xlim(c(0, 3)) + annotate("text", x = c(1.7, 1.7), y = c(0.6, 0.5), label = c(paste("Total barcodes expr#1:", nrow(DATA$e1)), paste("Total barcodes expr#2:", nrow(DATA$e2)))) + labs(title = "Individual expression replicates #1/#2. Trimmed")
#     plotList$dsplotE1E2_main <- dsplotE1E2_main + labs(title = "Individual expression replicates #1/#2")
#     plotList$dsplotE1E2_trim <- plotList$dsplotE1E2_main + xlim(c(0, 3)) + labs(title = "Individual expression replicates #1/#2. Trimmed")

#     # Density plot for expression replicate mean
#     dsplotEMean_main <- ggplot() + geom_density(colour="#ff9933", data=scList[[pi]], aes(mean)) + theme_bw() + xlab("Norm. expression value") + ylab("Expression count")
#     plotList$dsplotEMean_main <- dsplotEMean_main + annotate("text", x = 6, y = c(0.6), label = paste("Total barcodes expr. mean:", nrow(scList[[pi]]))) + labs(title = paste(synPI, pi, "Mean expression"))
#     plotList$dsplotEMean_trim <- plotList$dsplotEMean_main + xlim(c(0, 3)) + annotate("text", x = 1.7, y = c(0.6), label = paste("Total barcodes expr. mean:", nrow(scList[[pi]]))) + labs(title = paste(synPI, pi, "Mean expression. Trimmed"))

#     Cor.P <- paste("Pearson.Cor = ", round(with(scList[[pi]], cor(e1, e2, method="pearson", use="pairwise.complete.obs")), digits=2), sep = "")
#     Cor.S <- paste("Spearman.Cor = ", round(with(scList[[pi]], cor(e1, e2, method="spearman", use="pairwise.complete.obs")), digits=2), sep = "")
#     allValues <- paste("Values = ", nrow(scList[[pi]]), sep = "")
#     x.coord <- rep(max(scList[[pi]]$e1, na.rm = T) * 0.2, 3)
#     byVal <- -((max(scList[[pi]]$e2, na.rm = T) / 100) * 10)
#     y.coord <- seq(max(scList[[pi]]$e2, na.rm = T) * 0.7, by=byVal, length.out=3)
#     plotList$scplot <- ggplot(scList[[pi]], aes(e1, e2)) + geom_point(alpha = 5/10, colour = "#65A6D1", size = 3, na.rm = T) + xlab("Norm. expression replicate 1") + ylab("Norm. expression replicate 2") + annotate("text", x = x.coord, y = y.coord, label = c(Cor.P, Cor.S, allValues), size = 5, na.rm=T) + labs(title = paste(synPI, pi, "scatter plot correlation between expression replicate 1 and 2")) + theme_bw() + theme(plot.title = element_text(size = 18), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
#     pdf(file=file.path(output, paste(synPI, pi, "plots by expression and mapping data.pdf")), width=14, height=14)
#     print(multiplot(plotlist=plotList, cols=2, layout=matrix(c(1:6), ncol = 2, nrow = 3, byrow=T)))
#     dev.off()
# }

plotListComp <- list()
for (pop in names(popList)) {
    colorComp <- brewer.pal(length(pmi), 'Accent')
    names(colorComp) <- names(popList[[pop]])
    if (exists("dsPlotComp")) {rm(dsPlotComp)}
        for (pi in names(popList[[pop]])) {
            if (exists("dsPlotComp")) {
                dsPlotComp <- dsPlotComp + geom_density(colour=colorComp[[pi]], data=popList[[pop]][[pi]], aes(norm))
                count <- c(count, nrow(popList[[pop]][[pi]]))
                Prange <- rbind(Prange, data.frame("min" = round(range(popList[[pop]][[pi]]$norm)[1], 1), "max" = round(range(popList[[pop]][[pi]]$norm)[2], 1)))
            } else {
                dsPlotComp <- ggplot() + geom_density(colour=colorComp[[pi]], data=popList[[pop]][[pi]], aes(norm))
                count <- nrow(popList[[pop]][[pi]])
                Prange <- data.frame("min" = round(range(popList[[pop]][[pi]]$norm)[1], 1), "max" = round(range(popList[[pop]][[pi]]$norm)[2], 1))
            }
    }
    if (pop == "e18-1") {
        popName <- "Population_18-1"
    } else {
        popName <- "Population_18-2"
    }
    plotListComp[[pop]] <- dsPlotComp + theme_bw() + xlab("Norm. expression value") + ylab("Expression count") + labs(title = paste(popName, "expression density"))
    trimmed <- paste0(popName, "trim")

    xa <- rep(1, length(pmi))
    ya <- seq(1, by=-0.05, length.out=length(pmi))

    xaa <- rep(1.37, length(pmi))
    xaae <- rep(1.47, length(pmi))

    count <- paste("n:", count)
    xac <- rep(1.8, length(pmi))

    xar <- rep(2.2, length(pmi))
    Prange <- paste("r:", Prange$min, ":", Prange$max)

    plotListComp[[trimmed]] <- dsPlotComp + xlim(c(0, 3)) + ylim(c(0, 1)) + annotate("text", x = xa, y = ya, label = names(pmi)) + annotate("segment", x = xaa, xend = xaae, y = ya, yend = ya, colour = colorComp) + annotate("text", x = xac, y = ya, label = count) + annotate("text", x = xar, y = ya, label = Prange)
}
pdf(file=file.path(output, paste("Summary density plot by expression data in all promotor index.pdf")), width=14, height=14)
print(multiplot(plotlist=plotListComp, cols=2, layout=matrix(c(1:4), ncol = 2, nrow = 2, byrow=T)))
dev.off()

# for (pop in names(scList)){
#   colorComp <- brewer.pal(length(names(scList)), 'Accent')
#   names(colorComp) <- names(scList)
#   if (exists("dsplotEMean_comp")) {rm(dsplotEMean_comp)}
#   for (pi in names(scList)) {
#     if (exists("dsplotEMean_comp")) {
#       dsplotEMean_comp <- dsplotEMean_comp + geom_density(colour=colorComp[[pi]], data=scList[[pi]], aes(mean))
#     } else {
#       dsplotEMean_comp <- ggplot() + geom_density(colour=colorComp[[pi]], data=scList[[pi]], aes(mean))
#     }
#   }
#   plotListComp <- dsplotEMean_comp + theme_bw() + xlab("Norm. expression value") + ylab("Expression count") + labs(title = paste(pop, "Mean expression"))
#   trimmed <- paste0(pop, "trim")
#   plotListComp[[trimmed]] <- dsplotEMean_comp + xlim(c(0, 3)) + ylim(c(0, 1)) + scale_colour_manual("", breaks = names(colorComp), values = colorComp)
# }
# pdf(file=file.path(output, paste("Summary density plot by expression and mapping data in all promotor index.pdf")), width=14, height=14)
# print(multiplot(plotlist=plotListComp, cols=2))
# dev.off()
# #  Scatter Plot for NormExpr replicate 1/2
# Cor.P <- paste("Pearson.Cor = ", round(with(scDF, cor(e1, e2, method="pearson", use="pairwise.complete.obs")), digits=2), sep = "")
# Cor.S <- paste("Spearman.Cor = ", round(with(scDF, cor(e1, e2, method="spearman", use="pairwise.complete.obs")), digits=2), sep = "")
# allValues <- paste("Values = ", nrow(scDF), sep = "")
# x.coord <- rep(max(scDF$e1, na.rm = T) * 0.2, 3)
# y.coord <- seq(max(scDF$e2, na.rm = T) * 0.7, by=-0.3, length.out=3)
# 
# scplot <- ggplot(scDF, aes(e1, e2)) + geom_point(alpha = 5/10, colour = "#65A6D1", size = 3, na.rm = T) + xlab("Expression replicate 1") + ylab("Expression replicate 2") + annotate("text", x = x.coord, y = y.coord, label = c(Cor.P, Cor.S, allValues), size = 5, na.rm=T) + labs(title = "Scatter plot correlation between expression replicate 1 and 2") + theme_bw() + theme(plot.title = element_text(size = 18),axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
# 
# pdf(file=file.path(output, "Scatter plot correlation between expression replicate 1 and 2.pdf"), width=14, height=14)
# print(scplot)
# dev.off()

plotList <- list()
# Density plot for expression replicate 1/2
dsplotE1E2_main <- ggplot() + geom_density(colour="#000000", data=scDF, aes(e1)) + geom_vline(aes(xintercept=ControlDF$e1_norm), colour=c("#99ff00", "#99cc33", "#00ffff", "#006666")) + geom_density(colour="#8c762b", data=scDF, aes(e2)) + geom_vline(aes(xintercept=ControlDF$e2_norm), colour=c("#ff3366", "#990000", "#3333cc", "#9966ff")) + theme_bw() + xlab("Norm. expression value") + ylab("Expression count")
plotList$dsplotE1E2_main <- dsplotE1E2_main + annotate("text", x = c(6, 6), y = c(0.6, 0.5), label = c(paste("Total barcodes expr#1:", nrow(DATA$e1)), paste("Total barcodes expr#2:", nrow(DATA$e2)))) + labs(title = "Individual expression replicates #1/#2")
plotList$dsplotE1E2_trim <- plotList$dsplotE1E2_main + xlim(c(0, 3)) + annotate("text", x = c(1.7, 1.7), y = c(0.6, 0.5), label = c(paste("Total barcodes expr#1:", nrow(DATA$e1)), paste("Total barcodes expr#2:", nrow(DATA$e2)))) + labs(title = "Individual expression replicates #1/#2. Trimmed")

# Density plot for expression replicate mean
dsplotEMean_main <- ggplot() + geom_density(colour="#ff9933", data=scDF, aes(mean)) + geom_vline(aes(xintercept=cMean), colour=c("#336600", "#330066")) + theme_bw() + xlab("Norm. expression value") + ylab("Expression count")

plotList$dsplotEMean_main <- dsplotEMean_main + annotate("text", x = 6, y = c(0.6), label = paste("Total barcodes expr. mean:", nrow(scDF))) + labs(title = "Mean expression")

plotList$dsplotEMean_trim <- plotList$dsplotEMean_main + xlim(c(0, 3)) + annotate("text", x = 1.7, y = c(0.6), label = paste("Total barcodes expr. mean:", nrow(scDF))) + labs(title = "Mean expression. Trimmed")

# Density plot for mapping expression mean
# for (name in names(mappingList)){
#   main <- paste0(name, "_main")
#   trim <- paste0(name, "_trim")
#   label <- sub(".*([0-9]{2}).*", "\\1", name )
#   mutLabel <- sub(".*(mut)", "\\1", name)
#   tmpPlot <- ggplot() + geom_density(colour="#ff9933", data=mappingList[[name]], aes(mean)) + geom_vline(aes(xintercept=cMean), colour=c("#336600", "#330066")) + theme_bw() + xlab("Norm. expression value") + ylab("Expression count")
#   if (mutLabel != "mut"){
#       plotList[[main]] <- tmpPlot + annotate("text", x = 6, y = c(0.6), label = paste("Total barcodes in", label, "mean:", nrow(mappingList[[name]]))) + labs(title = paste("Mean mapping expression in", label))
#       plotList[[trim]] <- plotList[[main]] + xlim(c(0, 3)) + annotate("text", x = 1.7, y = c(0.6), label = paste("Total barcodes in", label, "mean:", nrow(mappingList[[name]]))) + labs(title = paste("Mean mapping expression in", label, "Trimmed"))
#   } else {
#       plotList[[main]] <- tmpPlot + annotate("text", x = 6, y = c(0.6), label = paste("Total mutations in", label, "mean:", nrow(mappingList[[name]]))) + labs(title = paste("Mean mapping expression in unique mut", label))
#       plotList[[trim]] <- plotList[[main]] + xlim(c(0, 3)) + annotate("text", x = 1.7, y = c(0.6), label = paste("Total mutations in", label, "mean:", nrow(mappingList[[name]]))) + labs(title = paste("Mean mapping expression in unique mut", label, "Trimmed"))
#   }
# }

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
