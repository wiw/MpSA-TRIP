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
library(jsonlite)
library(ggseqlogo)
library(RColorBrewer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
options(run.main=TRUE)
source(file.path("/home/anton/data/R-script/R-counts/RUN", "functions.R"))

options <- list(
                "wd" = "/home/anton/backup/input/trip/RUN_2018-11-20/",
                "motif_set" = "/home/anton/backup/input/trip/RUN_2018-11-20/results_12_12_2018/lib 33-40 1pct max expr.csv",
                "output" = "/home/anton/backup/input/trip/RUN_2018-11-20/results_12_12_2018",
                "motif_1" = "TTCGTTGC")

genome <- BSgenome.Hsapiens.UCSC.hg38
genome25names <- names(genome)[1:25]
annotated.genes <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(annotated.genes) <- genome25names # Resctrict area of 25 chromosomes
utr3 <- threeUTRsByTranscript(annotated.genes)

load_motif <- function(options) {
    motif_df <- read.csv(options$motif_set, header = FALSE, stringsAsFactors = FALSE)
    return(motif_df)
}

search_motif_into_genome <- function(motifs, genome, genome25names) {
    motif_list_stat <- list()
    for (mf_item in motifs) {
        CounterPal <- 0
        if (exists("result_df") == T) {
            rm(result_df)
        }
        for (chr in genome25names) {
            MatchPal <- matchPattern(mf_item, genome[[chr]])
            CounterPal <- CounterPal + length(MatchPal)
            if (length(MatchPal@ranges@start) != 0 && length(MatchPal@ranges@width) != 0) {
                if (exists("result_df") == T) {
                    result_df <- rbind(result_df, data.frame("chr" = as.factor(chr), "strand" = "*", "start" = MatchPal@ranges@start, "end" = MatchPal@ranges@start + MatchPal@ranges@width - 1 ))
                } else {
                    result_df <- data.frame("chr" = as.factor(chr), "strand" = "*", "start" = MatchPal@ranges@start, "end" = MatchPal@ranges@start + MatchPal@ranges@width - 1)
                }
            }
        }
        motif_list_stat[[mf_item]] <- result_df
    }
    return(motif_list_stat)
}

genome_length <- function(genome, genome25names) {
    ll <- c()
    for (chr in genome25names) {
        ll <- append(ll, length(genome[[chr]]))
    }
    return(sum(as.numeric(ll)))
}

count_coverage_whole <- function(l_motif, l_genome) {
    output <- c()
    for (mf_item in names(l_motif)) {
        coverage_value <- sum(width(GRanges(l_motif[[mf_item]])))
        output[[mf_item]] <- ( coverage_value / l_genome ) * 100
    }
    return(output)
}

format_utr3 <-  function(utr3) {
    utr3_df <- as.data.frame(utr3)
    utr3_df[utr3_df$strand %in% "+", 4] <- utr3_df[utr3_df$strand %in% "+", 5]
    utr3_df[utr3_df$strand %in% "+", 5] <- utr3_df[utr3_df$strand %in% "+", 5] + 100
    utr3_df[utr3_df$strand %in% "-", 5] <- utr3_df[utr3_df$strand %in% "-", 4]
    utr3_df[utr3_df$strand %in% "-", 4] <- utr3_df[utr3_df$strand %in% "-", 4] - 100
    utr3_gr <- GRanges(utr3_df)
    return(utr3_gr)
}

main <- function(options) {
    motif_df <- load_motif(options)
    p_w_m <- as.data.frame(PWM(motif_df[, 1]))
    motif_2 <- "TTTTTTCC" # handmade by visual examination p_w_m
    motifs <- c(options$motif_1, motif_2)
    lists_of_motif <- search_motif_into_genome(motifs, genome, genome25names)
    length_of_genome <- genome_length(genome, genome25names)
    coverage_whole_genome <- count_coverage_whole(lists_of_motif, length_of_genome)
    utr3f <- format_utr3(utr3)
}


if (getOption('run.main', default=TRUE)) {
   main(options)
}
