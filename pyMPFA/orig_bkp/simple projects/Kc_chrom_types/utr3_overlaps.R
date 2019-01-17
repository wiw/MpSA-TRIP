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
library(rowr)
library(reshape2)
options(run.main=TRUE)
source(file.path("/home/anton/data/R-script/R-counts/RUN", "functions.R"))

options <- list(
                "wd" = "/home/anton/backup/input/trip/RUN_2018-11-20/",
                "motif_set" = list("lib33_40" = "/home/anton/backup/input/trip/RUN_2018-11-20/results_12_12_2018/lib33_40.csv",
                                   "lib29_36" = "/home/anton/backup/input/trip/RUN_2018-11-20/results_12_12_2018/lib29_36.csv",
                                   "lib25_32" = "/home/anton/backup/input/trip/RUN_2018-11-20/results_12_12_2018/lib25_32.csv",
                                   "lib37_44" = "/home/anton/backup/input/trip/RUN_2018-11-20/results_12_12_2018/lib37_44.csv"),
                "output" = "/home/anton/backup/input/trip/RUN_2018-11-20/results_12_12_2018")

genome <- BSgenome.Hsapiens.UCSC.hg38
genome25names <- names(genome)[1:25]
annotated.genes <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(annotated.genes) <- genome25names # Resctrict area of 25 chromosomes
utr3 <- threeUTRsByTranscript(annotated.genes)

load_motif <- function(options) {
    motif_source_list <- list()
    for (motif_csv in names(options$motif_set)) {
        motif_source_list[[motif_csv]] <- read.csv(options$motif_set[[motif_csv]], header = FALSE, stringsAsFactors = FALSE)
    }
    return(motif_source_list)
}

search_motif_into_genome <- function(motif_source_list, genome, genome25names) {
    motif_list_stat <- list()
    for (motif_csv in names(motif_source_list)) {
        print(paste0("Search motifs in ", motif_csv))
        motif_list_stat[[motif_csv]] <- list()
        motif_vector <- motif_source_list[[motif_csv]][, 1]
        for (mf_item in motif_vector) {
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
            motif_list_stat[[motif_csv]][[mf_item]] <- result_df
            print(paste0("Motif ", mf_item, " complete"))
        }
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

count_coverage_whole <- function(lists_of_motif, length_of_genome) {
    output <- list()
    for (motif_csv in names(lists_of_motif)) {
        output[[motif_csv]] <- c()
        for (mf_item in names(lists_of_motif[[motif_csv]])) {
            tmp_data <- lists_of_motif[[motif_csv]][[mf_item]]
            coverage_count <- nrow(tmp_data)
            output[[motif_csv]][[mf_item]] <- (coverage_count / length_of_genome) * 10^6
        }
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

count_coverage_utr3 <- function(utr3f, lists_of_motif, length_of_genome) {
    output <- list()
    for (motif_csv in names(lists_of_motif)) {
        output[[motif_csv]] <- c()
        for (mf_item in names(lists_of_motif[[motif_csv]])) {
            tmp_data <- lists_of_motif[[motif_csv]][[mf_item]]
            motif_gr <- makeGRangesFromDataFrame(tmp_data, keep.extra.columns = T)
            overlap_value <- as.data.frame(findOverlaps(utr3f, motif_gr))
            coverage_count <- nrow(overlap_value)
            output[[motif_csv]][[mf_item]] <- (coverage_count / (length(utr3f) * 100) ) * 10^6
        }
    }
    return(output)
}

make_df_write_result <- function(lists_of_motif, coverage_whole_genome, coverage_utr3_genome, options) {
    df_massive <- list()
    for (motif_csv in names(lists_of_motif)) {
        df_list <- list("Motif" = names(coverage_whole_genome[[motif_csv]]),
                         # "Abundance in the Genome" = "",
                         # "Genome Length" = "",
                         "Frequency.in.the.Genome.per.Mb" = coverage_whole_genome[[motif_csv]],
                         # "Abundance in Gene3 Regions" = "",
                         # "Length of Gene3 Regions" = "",
                         "Frequency.in.Gene3.Regions.per.Mb" = coverage_utr3_genome[[motif_csv]])
        df <- cbind.fill(df_list$Motif, df_list$Frequency.in.the.Genome.per.Mb, df_list$Frequency.in.Gene3.Regions.per.Mb, fill = NA)
        names(df) <- c("Motif", "Frequency.in.the.Genome.per.Mb", "Frequency.in.Gene3.Regions.per.Mb")
        output_filename <- file.path(options$output, paste0(motif_csv, " Genome Coverage Stat.csv"))
        write.csv2(df, output_filename, row.names = FALSE)
        df_massive[[motif_csv]] <- df
    }
    return(df_massive)
}


perform_t.test <- function(df, options) {
    result <- c()
    for (lib in levels(df$library)) {
        result[[lib]] <- wilcox.test(value ~ variable, data = df, paired = TRUE, subset = library %in% lib)$p.value
    }
    result <- paste0("p.value = ", format(result, digits=7, width=11), sep="")
    result <- data.frame("p_value" = result,
                         "x" = 1.5,
                         "y" = Inf,
                         "library" = levels(df$library),
                         stringsAsFactors=FALSE)
    write.csv2(result, file = file.path(options$output, "p_values.csv"), row.names = FALSE)
}

prepare2boxplot <- function(df_list) {
    melted_df <- list()
    for (motif_csv in names(df_list)) {
        melted_df[[motif_csv]] <- melt(df_list[[motif_csv]])
    }
    rebuild_df <- do.call("rbind", melted_df)
    library.vector <- rownames(rebuild_df)
    library.vector <- sub("([_a-z0-9].*)\\..*", "\\1", library.vector)
    rebuild_df$library <- as.factor(library.vector)
    rownames(rebuild_df) <- NULL
    return(rebuild_df)
}

make_boxplot <- function(rebuild_df) {
    graphs_list <- list()
    main_graph <- ggplot(aes(x=variable, y=value, fill=variable), data=rebuild_df) + 
    geom_boxplot(colour = "#595959") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ylab("Frequency, per.Mb") + 
    facet_grid(. ~ library)
    # geom_text(aes(x, y, label = p_value), data = p__values, vjust = 1)
    graphs_list$normal <- main_graph
    return(graphs_list)
}

draw_graph <- function(graphs_list, options, custom.label = "") {
    if (length(graphs_list) <= 4) {
        cols_value <- ceiling( length(graphs_list) / 2 )
    } else {
        cols_value <- 2
    }
    # Draw boxplots
    pdf( file = file.path( options$output, paste0( "Boxplots from ", basename(options$wd), format(Sys.time(), " %Y-%m-%d"), custom.label, ".pdf" ) ), width=14, height=14 )
    print(multiplot(plotlist = graphs_list, cols = cols_value))
    dev.off()
}

main <- function(options) {
    motif_source_list <- load_motif(options)
    # p_w_m <- as.data.frame(PWM(motif_df[, 1]))
    # motif_2 <- "TTTTTTCC" # handmade by visual examination p_w_m
    # motifs <- c(options$motif_1, motif_2)
    lists_of_motif <- search_motif_into_genome(motif_source_list, genome, genome25names)
    length_of_genome <- genome_length(genome, genome25names)
    coverage_whole_genome <- count_coverage_whole(lists_of_motif, length_of_genome)
    utr3f <- format_utr3(utr3)
    coverage_utr3_genome <- count_coverage_utr3(utr3f, lists_of_motif, length_of_genome)
    df_list <- make_df_write_result(lists_of_motif, coverage_whole_genome, coverage_utr3_genome, options)
    rebuild_df <- prepare2boxplot(df_list)
    graphs_list <- make_boxplot(rebuild_df)
    draw_graph(make_boxplot(rebuild_df), options)
    perform_t.test(rebuild_df, options)
}


if (getOption('run.main', default=TRUE)) {
   main(options)
}