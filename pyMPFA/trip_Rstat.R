#!/usr/bin/R
library(ggplot2)
args = commandArgs(TRUE)
inputFile = args[1]
outputDir = args[2]
indexLabel = args[3]
bc <- read.delim(inputFile)
bc$BCLength <- nchar(as.character(bc$Barcode))
bc <- bc[,-1]
bc.expanded <- bc[rep(row.names(bc), bc$BCSequenceCount), 2]
bcData <- data.frame("bcLength" = bc.expanded)
variance <- sort(unique(bcData$bcLength))
lengthStat <- vector()
for (i in variance) {
    lengthStat <- append(lengthStat, length(bcData$bcLength[bcData$bcLength == i]))
}
x.coord <- variance
y.coord <- lengthStat + max(lengthStat) * 0.07
bcHist <- ggplot(bcData, aes(x=bcLength)) + geom_histogram(binwidth=1) + labs(title="Distribution of barcodes by length", x="Length of barcode", y="Barcode count") + annotate("text", x = x.coord, y = y.coord, label = lengthStat, size = 6, na.rm=T) + theme(axis.title.x=element_text(size=28), axis.title.y=element_text(size=28), axis.text.x=element_text(size=20), axis.text.y=element_text(size=20), plot.title=element_text(size=30, hjust=0.5), legend.position = "none")
pdf(file=file.path(outputDir, paste(indexLabel, "Barcode_distribution_length.pdf", sep="_")), width=18, height=12)
print(bcHist)
dev.off()