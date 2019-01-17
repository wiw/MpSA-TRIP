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
library(reticulate)
options(run.main=TRUE)
source(file.path("/home/anton/data/R-script/R-counts/RUN", "functions.R"))

options <- list(
                "wd" = "/home/anton/backup/input/trip/RUN_2018-12-14/results/statistics/trip_18_2",
                "output" = "/home/anton/backup/input/trip/RUN_2018-12-14/report_trip18_2018_12",
                "promoter_indexes" = c("HexA" = 'AGCTC',
                                       "Hsp70" = 'ACGTA',
                                       "MtnA" = 'CTGCT', 
                                       "PCNA" = 'AGTCA', 
                                       "Pyk" = 'TCAA', 
                                       "Tbp" = 'TTGAG', 
                                       "promoterless" = 'TCGCT'),
                "target_file" = "map_norm_data.json",
                "rpl_file" = "rpl_count.json",
                "pyScrypt" = "/home/anton/data/simple projects/Kc_chrom_types/python_load_pickle.py",
                "raw_path" = "/home/anton/backup/input/trip/RUN_2018-12-14/results/trip_18_2/sample_S1_L001_R1_001_trip_18_2_mapping/Dump")


load_json <- function(json_file) {
    if (file.exists(json_file) && !dir.exists(json_file)) {
        loaded_data <- jsonlite::fromJSON(json_file, flatten = T)
    }
    return(loaded_data)
}


main_loading <- function(options) {
    output <- list()
    for (p_index in names(options$promoter_indexes)) {
        # load mapping and normalization
        json_file <- file.path(options$wd,
                               options$promoter_indexes[[p_index]],
                               options$target_file)
        json_set <- load_json(json_file)
        output[[p_index]] <- convert_json_set(json_set)
        # load read per million count
        rpl_file <- file.path(options$wd,
                              options$promoter_indexes[[p_index]],
                              options$rpl_file)
        output[[p_index]][["rpl"]] <- as.data.frame(load_json(rpl_file), stringsAsFactors = FALSE)
    }
    return(output)
}


rpm <- function(V1, int) {
    return((V1/int)*10^6)
}


normalization <- function(V1, V2) {
    return(V1/V2)
}

add_raw_mapping_data <- function(options) {
    load_pickle <- function(filepath, options) {
        source_python(file.path(options$pyScrypt))
        pckl.data <- read_pickle_file(filepath)
        return(pckl.data)
    }
    convert_to_list_of_df <- function(pckl.data, options) {
        pickle.list <- list()
        for (promoter_name in names(options$promoter_indexes)) {
            promoter_seq <- options$promoter_indexes[[promoter_name]]
            pickle.list[[promoter_name]] <- data.frame("barcode" = names(pckl.data[[promoter_seq]]))
        }
        return(pickle.list)
    }
    bc.dict.vector <- list.files(path = options$raw_path, pattern = ".*bcDict.*", full.names = TRUE)
    raw.map.list <- list()
    for (replicate in bc.dict.vector) {
        replicate.pckl.data <- load_pickle(replicate, options)
        replicate.name <- strsplit(basename(replicate), "[.]")[[1]][1]
        raw.map.list[[replicate.name]] <- convert_to_list_of_df(replicate.pckl.data, options)
    }
    return(raw.map.list)
}

merge_raw_data <- function(raw.map.list, options) {
    replicate.vector <- names(raw.map.list)
    output.list <- list()
    for (promoter_name in names(options$promoter_indexes)) {
        tmp.df <- rbind(raw.map.list[[replicate.vector[1]]][[promoter_name]], raw.map.list[[replicate.vector[2]]][[promoter_name]])
        output.list[[promoter_name]] <- data.frame("barcode" = unique(tmp.df$barcode))
    }
    return(output.list)
}

map_raw_intersect_norm <- function(raw.map.data, norm.data, options) {
    raw.intersect.count <- c()
    for (promoter_name in names(options$promoter_indexes)) {
        count_value <- nrow(merge(raw.map.data[[promoter_name]], norm.data[[promoter_name]][['normalization']], by = "barcode"))
        raw.intersect.count[[promoter_name]] <- paste0("raw.mapping: ", count_value)
    }
    return(raw.intersect.count)
}

main.raw.data <- function(norm.data, options) {
    raw.map.list <- add_raw_mapping_data(options)
    raw.map.data <- merge_raw_data(raw.map.list, options)
    raw.intersect.count <- map_raw_intersect_norm(raw.map.data, norm.data, options)
    return(raw.intersect.count)
}

count_by_factor <- function(intersected.data, raw.intersect.count) {
    all_counts <- c()
    for(p_index in levels(intersected.data$promoter.index)) {
        tmp_v <- c()
        for (exp in levels(intersected.data$experiment)) {
            count_value <- nrow( intersected.data[ intersected.data$experiment %in% exp & intersected.data$promoter.index %in% p_index, ] )
            tmp_v <- append( tmp_v, paste0(exp, ": ", count_value))
        }
        tmp_v <- append(tmp_v, raw.intersect.count[[p_index]])
        all_counts[[p_index]] <- paste0(tmp_v, collapse = "\n")
    }
    return(all_counts)
}

custom_mean <- function(df, col_vector) {
    return(rowMeans(df[, col_vector]))
}


json2df <- function(list.of.single.string.list) {
    first_column <- names(list.of.single.string.list)
    mtx <- as.data.frame(matrix(unlist(list.of.single.string.list), nrow = length(list.of.single.string.list), byrow=T), stringsAsFactors=F)
    mtxDF <- lapply(mtx, function(x) {x1 <- type.convert(as.character(x))
                                            if(is.factor(x1))
                                                as.character(x1)
                                            else x1})
    converted_df <- cbind(data.frame("V0" = first_column, stringsAsFactors=F), mtxDF, stringsAsFactors=F)
    return(converted_df)
}

convert_json_set <- function(json_set) {
    output <- list()
    for (item in names(json_set)) {
        output[[item]] <- json2df(json_set[[item]])
    }
    return(output)
}

prepare_data <- function(data_set) {
    output <- list()
    for (p_index in names(data_set)) {
        output[[p_index]][["normalization"]] <- data.frame("barcode" = data_set[[p_index]][["normalization-1"]]$V0,
                                        "normalization_1" = data_set[[p_index]][["normalization-1"]]$V1,
                                        "normalization_2" = data_set[[p_index]][["normalization-2"]]$V1)
        for (value in c("normalization_1", "normalization_2")) {
            rpm_name <- paste0(value, ".rpm")
            V1 <- output[[p_index]][["normalization"]][[value]]
            rpm_value <- as.integer(data_set[[p_index]]$rpl[data_set[[p_index]]$rpl$V1 %in% value, 2])
            output[[p_index]][["normalization"]][[rpm_name]] <- rpm(V1, rpm_value)
        }
        output[[p_index]][["normalization"]]$normalization.mean <- custom_mean(output[[p_index]][["normalization"]], c(4, 5))
    }
    return(output)
}


map_norm_intersect <- function(downloaded.data, norm.data) {
    output <- list()
    for (p_index in names(norm.data)) {
        tmp_output <- merge(norm.data[[p_index]]$normalization, downloaded.data[[p_index]]$mapping, by.x = "barcode", by.y = "V0")
        tmp_list <- list("normalization" = norm.data[[p_index]]$normalization$normalization.mean,
                                  "mapping.intersect" = tmp_output$normalization.mean)
        output[[p_index]] <- melt(tmp_list)
        names(output[[p_index]]) <- c("value", "experiment")
        output[[p_index]]$experiment <- as.factor(output[[p_index]]$experiment)
    }
    # Melting output data
    output.melt <- melt(output)
    names(output.melt) <- c("experiment", "undef", "value", "promoter.index")
    output.melt$promoter.index <- as.factor(output.melt$promoter.index)
    return(output.melt)
}


draw_density <- function(intersected.data, norm.data, options) {
    raw.intersect.count <- main.raw.data(norm.data, options)
    label.text <- count_by_factor(intersected.data, raw.intersect.count)
    gplot_item <- ggplot(data = intersected.data, aes(x = log2(value), fill = experiment)) + 
        geom_density(alpha = 0.25) + 
        ylab("Frequency of normalization values") + 
        xlab("log2 from density values") + 
        facet_wrap(vars(promoter.index)) + 
        annotate( "text", 
            x = 6.5, 
            y = 0.2, 
            hjust = 0, 
            label = label.text )
    return(list("data" = gplot_item))
}

draw_graph <- function(graphs_list, options, custom.label = "") {
    if (length(graphs_list) <= 4) {
        cols_value <- ceiling( length(graphs_list) / 2 )
    } else {
        cols_value <- 2
    }
    # Draw Denity
    pdf( file = file.path( options$output, paste0( "Density from ", basename(options$wd), format(Sys.time(), " %Y-%m-%d"), custom.label, ".pdf" ) ), width=14, height=14 )
    if (cols_value == 1) {
        print(graphs_list[[1]])
    } else {
        print(multiplot(plotlist = graphs_list, cols = cols_value))
    }
    dev.off()
}

main <- function(options) {
    downloaded.data <- main_loading(options)
    norm.data <- prepare_data(downloaded.data)
    intersected.data <- map_norm_intersect(downloaded.data, norm.data)
    draw_graph(draw_density(intersected.data, norm.data, options), options)
}


if (getOption('run.main', default=TRUE)) {
    main(options)
}