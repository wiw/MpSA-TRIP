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
options(run.main=TRUE)
source(file.path("/home/anton/data/R-script/R-counts/RUN", "functions.R"))
args <- commandArgs(TRUE)

options <- list(
                "wd" = "/home/anton/backup/input/trip/RUN_2018-12-14/results/statistics/trip_18_2",
                "output" = "/home/anton/backup/input/trip/RUN_2018-12-14/report_trip18_2018_12",
                "kc_colors" = "/home/anton/backup/input/trip/RUN_2018-06-07/philion/GSE22069_kc167_r6_reformat.tsv",
                "promoter_indexes" = c('AGCTC', 'ACGTA', 'CTGCT', 'AGTCA', 'TCAA', 'TTGAG', 'TCGCT'),
                "promoter_index_names" = c("HexA", "Hsp70", "MtnA", "PCNA", "Pyk", "Tbp", "promoterless"))

kc_load <- function(options){
    kc_colors = read.delim(file = options$kc_colors, header = TRUE, sep = "\t", stringsAsFactors = F, dec = ",")
    kc_colors$chromatin <- as.factor(kc_colors$chromatin)
    return(kc_colors)
}

data_load <- function(options){
    empty_data <- list()
    for (p_index in options$promoter_indexes) {
        results_path <- file.path(options$wd, p_index)
        results_file <- list.files(results_path, pattern = "results.*.csv", full.names = TRUE)
        empty_data[[p_index]] <- read.delim(file = results_file, header = TRUE, sep = ";", stringsAsFactors = F, dec = ",")
        empty_data[[p_index]]$end <- empty_data[[p_index]]$start
        names(empty_data[[p_index]])[names(empty_data[[p_index]]) == "chr"] <- "seqnames"
        empty_data[[p_index]]$seqnames <- as.factor(empty_data[[p_index]]$seqnames)
    }
    return(empty_data)
}

assign_kc_to_data <- function(kc_colors, data){
    merged_data <- list()
    for (p_index in names(data)){
        gr_kc_colors <- makeGRangesFromDataFrame(kc_colors, keep.extra.columns = T)
        gr_data <- makeGRangesFromDataFrame(data[[p_index]], keep.extra.columns = T)
        results <- as.data.frame(findOverlaps(gr_data, gr_kc_colors))
        merged_data[[p_index]] <- cbind(data[[p_index]][results$queryHits, ], data.frame("kc_colors" = kc_colors[results$subjectHits, ]))
    }
    return(merged_data)
}

add_null_values <- function(colour_data) {
    for (p_index in names(colour_data)) {
        uniq_levels <- unique(colour_data[[p_index]]$kc_colors.chromatin)
        all_levels <- levels(colour_data[[p_index]]$kc_colors.chromatin)
        diff_levels <- all_levels[!(all_levels %in% uniq_levels)]
        last_row <- colour_data[[p_index]][nrow(colour_data[[p_index]]),]
        for (column in names(last_row)) {
            if (is.integer(last_row[[column]]) | is.double(last_row[[column]])) {
                last_row[[column]] <- 0
            }
        }
        for (level_items in diff_levels) {
            last_row$kc_colors.chromatin <- factor(level_items)
            colour_data[[p_index]] <- rbind(colour_data[[p_index]], last_row)
        }
    }
    return(colour_data)
}

make_boxplot <- function(colour_data, options, trim = FALSE) {
    get_trim_value <- function(tables) {
        boxplot_stat <- boxplot(expression_mean_norm ~ kc_colors.chromatin, data = tables)
        trim_value <- max(boxplot_stat$conf) * 1.25
        return(trim_value)
    }
    main_boxplot = list()
    for (p_index in names(colour_data)) {
        index_set = options$promoter_indexes
        names(index_set) = options$promoter_index_names
        current_p_index_name = names(index_set[index_set == p_index])
        colours_set = sort(unique(colour_data[[p_index]]$kc_colors.chromatin))
        main_boxplot[[current_p_index_name]] <- ggplot(colour_data[[p_index]], aes(kc_colors.chromatin, expression_mean_norm, fill = kc_colors.chromatin)) + 
        geom_boxplot(colour = "#595959") + 
        scale_fill_manual(breaks = colours_set, values = tolower(colours_set)) + 
        xlab("Kc167 Chromatin domains") +
        ylab("Norm., mean expression") +
        labs(title = paste0("Mean expression by chromatin types in ", current_p_index_name))
        if (trim) {
            trim_value <- get_trim_value(colour_data[[p_index]])
            current_p_index_name_trim <- paste0(current_p_index_name, "_trim")
            main_boxplot[[current_p_index_name]] <- main_boxplot[[current_p_index_name]] + ylim(c(0, trim_value))
        }
    }
    return(main_boxplot)
}

make_boxplot_by_promoter <- function(colour_data, options) {
    summ_boxplot_list <- list()
    current_p_index_name <- basename(options$wd)
    names(colour_data) <- options$promoter_index_names
    reforma_data <- do.call("rbind", colour_data)
    promoter_vector <- rownames(reforma_data)
    promoter_vector <- sub("([a-zA-Z0-9]*)\\.([0-9\\.]*)", "\\1", promoter_vector)
    reforma_data$promoter <- as.factor(promoter_vector)
    rownames(reforma_data) <- NULL
    core_boxplot <- ggplot(aes(x = promoter, y = expression_mean_norm, fill = promoter), data = reforma_data) + 
    geom_boxplot(colour = "#595959") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab("Promoter indexes") +
    ylab("Norm., mean expression") +
    facet_grid(. ~ kc_colors.chromatin)
    summ_boxplot_list$alldata <- core_boxplot + 
    # coord_cartesian(ylim = c(0,8)) + 
    labs(title = paste0("Mean expression by chromatin types & promoter indexes in ", current_p_index_name))
    summ_boxplot_list$trimdata <- core_boxplot + 
    coord_cartesian(ylim = c(0,8)) + 
    labs(title = paste0("Trimmed mean expression by chromatin types & promoter indexes in ", current_p_index_name))
    return(summ_boxplot_list)
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

write_table <- function(df, options) {
    file_name <- file.path(options$output, paste0("boxplot_colour_stat", basename(options$wd), ".csv"))
    write.csv2(df, file_name, row.names = F)
}


main_draw <- function(colour_data, colour_stat, options) {
    boxplot_data <- make_boxplot(colour_data, options)
    boxplot_data_t <- make_boxplot(colour_data, options, trim = T)
    boxplot_data <- c(boxplot_data, boxplot_data_t)
    summ_boxplot_data <- make_boxplot_by_promoter(colour_data, options)
    draw_graph(boxplot_data, options)
    draw_graph(summ_boxplot_data, options, custom.label = "_summary_boxplot_data")
    write_table(colour_stat, options)
}

generate_table_stats <- function(data, kc_colors) {
    counts_by_promoter <- sapply(data, nrow)
    colour_data <- assign_kc_to_data(kc_colors, data)
    output <- data.frame("promoter_indexes" = names(counts_by_promoter),
                         "genuine bc's count" = counts_by_promoter)
    rownames(output) <- NULL
    by_chromatin_count <- lapply(colour_data, function(x) { count(x, "kc_colors.chromatin") })
    tmp_df <- data.frame(by_chromatin_count[1])
    names(tmp_df) <- c("kc_colors.chromatin", names(by_chromatin_count[1]))
    for (p_index in names(by_chromatin_count)[! names(by_chromatin_count) %in% names(by_chromatin_count[1])]) {
        tmp_df <- merge(tmp_df, by_chromatin_count[[p_index]], by = "kc_colors.chromatin", all = TRUE)
        names(tmp_df)[ncol(tmp_df)] <- p_index
    }
    names_df <- as.character(tmp_df[, 1])
    tmp_df <- as.data.frame(t(tmp_df))
    rownames(tmp_df) <- NULL
    tmp_df <- tmp_df[-1, ]
    names(tmp_df) <- names_df
    output <- cbind(output, tmp_df)
    return(output)
}

main <- function(options) {
    kc_colors <- kc_load(options)
    data <- data_load(options)
    colour_data <- add_null_values(assign_kc_to_data(kc_colors, data))
    colour_stat <- generate_table_stats(data, kc_colors)
    main_draw(colour_data, colour_stat, options)
}

if (getOption('run.main', default=TRUE)) {
   main(options)
}