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
options(run.main=TRUE)
source(file.path("/home/anton/data/R-script/R-counts/RUN", "functions.R"))
args <- commandArgs(TRUE)

# wd <- "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics/Lib_33-40_bcRead-1_bcMutated-0"
# control <- "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics/Lib_33-40_bcRead-1_bcMutated-0/control.json"
# data <- "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics/Lib_33-40_bcRead-1_bcMutated-0/data.json"
# rpl_count <- "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics/Lib_33-40_bcRead-1_bcMutated-0/rpl_count.json"
# arguments <- list("wd"=wd, "control"=control, "data"=data, "rpl"=rpl_count)

parse_args <- function(args) {
    wd <- args[1]
    control <- file.path(args[1], args[2])
    data <- file.path(args[1], args[3])
    rpl_count <- file.path(args[1], args[4])
    return(list("wd"=wd, "control"=control, "data"=data, "rpl"=rpl_count))
}
load_json <- function(json_file) {
    if (file.exists(json_file) && !dir.exists(json_file)) {
        loaded_data <- jsonlite::fromJSON(json_file, flatten = T)
    }
    return(loaded_data)
}
transposition_table <- function(df, column_number_who_transposed) {
    tt_int <- as.integer(paste0(column_number_who_transposed, collapse=""))
    if (is.integer(tt_int) && tt_int > 0 && tt_int <= ncol(df)){
        tt <- as.data.frame(df[, -tt_int])
        rownames(tt) <- df[, tt_int]
        ttt <- as.data.frame(t(tt), stringsAsFactors=F)
        ttt[[names(df)[tt_int]]] <- rownames(ttt)
        rownames(ttt) <- NULL
        ttt <- ttt[, c(ncol(ttt), 1:ncol(ttt)-1)]
    } else {
        print(paste("Checking column number. You are enter:", tt_int))
        break
    }
    return(ttt)
}
convert_data <- function(loaded_jsn, first_item_name) {
    if (first_item_name == "control") {
        synonymous <- list("e" = "expression", "m" = "mapping", "n" = "normalization")
        ab <- sub("control_", "", names(loaded_jsn))
        bc <- sub("-", "_", ab)
        first_column <- c()
        for (element in bc) {
            marker <- unlist(strsplit(element, "_"))
            first_column <- append(first_column, paste(synonymous[[marker[1]]], marker[2], sep="_"))
        }
    } else {
        first_column <- names(loaded_jsn)
    }
    mtx <- as.data.frame(matrix(unlist(loaded_jsn), nrow = length(loaded_jsn), byrow=T), stringsAsFactors=F)
    mtxDF <- cbind(as.data.frame(sapply(mtx[,c(1:ncol(mtx)-1)], as.numeric)), mtx[, ncol(mtx)])
    converted_df <- cbind(data.frame("V0" = first_column, stringsAsFactors=F), mtxDF)
    names(converted_df) <- c(first_item_name, sub("-", "_", names(loaded_jsn[[1]])))
    return(converted_df)
}

main_loading <- function(arguments) {
    output <- list()
    for (item in names(arguments)[!names(arguments) %in% "wd"]) {
        if (item == "rpl") {
            output$rpl <- as.data.frame(load_json(arguments[[item]]), stringsAsFactors=F)
            names(output$rpl) <- c("replicate", "count")
            output$rpl$count <- as.integer(output$rpl$count)
        } else {
            if (item == "control") {
                first_item_name <- item
            } else {
                first_item_name <- "barcode"
            }
            output[[item]] <- convert_data(load_json(arguments[[item]]), first_item_name)
        }
        if (item != "data") {
            output[[item]] <- transposition_table(output[[item]], 1)
        }
    }
    return(output)
}

rpm <- function(V1, int) {
    return((V1/int)*10^6)
}
normalization <- function(V1, V2) {
    return(V1/V2)
}
remove_unv_data <- function(df, col_vector_two) {
    if (length(col_vector_two) == 2) {
        col1 <- col_vector_two[1]
        col2 <- col_vector_two[2]
        df <- df[!is.na(df[,col1]) & !is.na(df[,col2]),]
        df <- df[!is.infinite(df[,col1]) & !is.infinite(df[,col2]),]
    } else {
        print("Use only two column from your data frame. Exit...")
        break
    }
    return(df)
}
custom_mean <- function(df, col_vector) {
    return(rowMeans(df[, col_vector]))
}

prepare_data <- function(DATA) {
    for (set in names(DATA)[!names(DATA) %in% "rpl"]) {
        target_vector <- grep("_", names(DATA[[set]]))
        # Count reads per millon
        for (experiment in names(DATA[[set]])[target_vector]) {
            rpm_column_name <- paste(experiment, "rpm", sep="_")
            DATA[[set]][[rpm_column_name]] <- rpm(DATA[[set]][[experiment]], DATA$rpl[[experiment]])
        }
        # Normalization expression data
        for (rpm in names(DATA[[set]][grep("rpm", names(DATA[[set]]))])){
            if (length(grep("expression", rpm)) != 0) {
                expr <- paste(rpm, "norm", sep="_")
                norm <- sub("expression", "normalization", rpm)
                DATA[[set]][[expr]] <- normalization(DATA[[set]][[rpm]], DATA[[set]][[norm]])
            }
        }
        expr_norm_vector <- grep("expr.*norm", names(DATA[[set]]))
        expr_mean_name <- "expression_mean"
        rudt <- remove_unv_data(DATA[[set]], expr_norm_vector)
        if (ncol(rudt) != ncol(DATA[[set]])) {
            print("Found NA. I don't know why it happened. Take your own. Exit...")
            break
        }
        DATA[[set]][[expr_mean_name]] <- custom_mean(DATA[[set]], expr_norm_vector)
    }
    return(DATA)
}

scatter_plots <- function(DATA) {
    Cor.P <- paste("Pearson.Cor = ", round(with(DATA$data, cor(expression_1_rpm_norm, expression_2_rpm_norm, method="pearson", use="pairwise.complete.obs")), digits=2), sep = "")
    Cor.S <- paste("Spearman.Cor = ", round(with(DATA$data, cor(expression_1_rpm_norm, expression_2_rpm_norm, method="spearman", use="pairwise.complete.obs")), digits=2), sep = "")
    allValues <- paste("Values = ", nrow(DATA$data), sep = "")
    x.coord <- rep(max(DATA$data$expression_1_rpm_norm, na.rm = T) * 0.2, 3)
    y.coord <- seq(max(DATA$data$expression_2_rpm_norm, na.rm = T) * 0.7, by=-0.3, length.out=3)
    scplot <- ggplot(DATA$data, aes(expression_1_rpm_norm, expression_2_rpm_norm)) + geom_point(alpha = 5/10, colour = "#65A6D1", size = 3, na.rm = T) + xlab("Expression replicate 1") + ylab("Expression replicate 2") + annotate("text", x = x.coord, y = y.coord, label = c(Cor.P, Cor.S, allValues), size = 5, na.rm=T) + labs(title = "Scatter plot correlation between expression replicate 1 and 2") + theme_bw() + theme(plot.title = element_text(size = 18),axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
    return(scplot)
}
density_plots <- function(DATA, threshold){
    if (!is.numeric(threshold)) {
        print("Please enter numeric value to threshold! Exit...")
        break
    }
    plotList <- list()
    # Density plot for expression replicate 1/2
    dsplotE1E2_main <- ggplot() + geom_density(colour="#000000", data=DATA$data, aes(expression_1_rpm_norm)) + geom_vline(aes(xintercept=DATA$control$expression_1_rpm_norm), colour=c("#99ff00", "#99cc33", "#00ffff", "#006666")) + geom_density(colour="#8c762b", data=DATA$data, aes(expression_2_rpm_norm)) + geom_vline(aes(xintercept=DATA$control$expression_2_rpm_norm), colour=c("#ff3366", "#990000", "#3333cc", "#9966ff")) + theme_bw() + xlab("Norm. expression value") + ylab("Expression count")
    plotList$dsplotE1E2_main <- dsplotE1E2_main + annotate("text", x = c(6, 6), y = c(0.6, 0.5), label = c(paste("Total barcodes expr#1:", length(DATA$data$expression_1_rpm_norm[!is.na(DATA$data$expression_1_rpm_norm)])), paste("Total barcodes expr#2:", length(DATA$data$expression_2_rpm_norm[!is.na(DATA$data$expression_2_rpm_norm)])))) + labs(title = "Individual expression replicates #1/#2")
    plotList$dsplotE1E2_trim <- plotList$dsplotE1E2_main + xlim(c(0, threshold)) + annotate("text", x = c(1.7, 1.7), y = c(0.6, 0.5), label = c(paste("Total barcodes expr#1:", length(DATA$data$expression_1_rpm_norm[!is.na(DATA$data$expression_1_rpm_norm)])), paste("Total barcodes expr#2:", length(DATA$data$expression_2_rpm_norm[!is.na(DATA$data$expression_2_rpm_norm)])))) + labs(title = "Individual expression replicates #1/#2. Trimmed")
    # Density plot for expression replicate mean
    wt <- grep("wt", DATA$control$control)
    delta <- grep("delta", DATA$control$control)
    expr_mean_name <- grep("mean", names(DATA$control))
    control_mean <- sort(c(mean(DATA$control[[expr_mean_name]][wt]), mean(DATA$control[[expr_mean_name]][delta])))
    dsplotEMean_main <- ggplot() + geom_density(colour="#ff9933", data=DATA$data, aes(expression_mean)) + geom_vline(aes(xintercept=control_mean), colour=c("#336600", "#330066")) + theme_bw() + xlab("Norm. expression value") + ylab("Expression count")
    plotList$dsplotEMean_main <- dsplotEMean_main + annotate("text", x = 6, y = c(0.6), label = paste("Total barcodes expr. mean:", nrow(DATA$data))) + labs(title = "Mean expression")
    plotList$dsplotEMean_trim <- plotList$dsplotEMean_main + xlim(c(0, threshold)) + annotate("text", x = 1.7, y = c(0.6), label = paste("Total barcodes expr. mean:", nrow(DATA$data))) + labs(title = "Mean expression. Trimmed")
    return(plotList)
}

main <- function(args) {
    arguments <- parse_args(args)
    DATA <- main_loading(arguments)
    DATA <- prepare_data(DATA)
    scatter_list <- scatter_plots(DATA)
    density_list <- density_plots(DATA, 4)
    pdf(file=file.path(arguments$wd, paste0("Scatter plot correlation between expression replicate 1 and 2_", format(Sys.time(), "%Y-%m-%d"), ".pdf")), width=14, height=14)
    print(scatter_list)
    dev.off()
    pdf(file=file.path(arguments$wd, paste0("Density plot by expression data_", format(Sys.time(), "%Y-%m-%d"), ".pdf")), width=14, height=14)
    print(multiplot(plotlist=density_list, cols=2, layout=matrix(c(1:4), ncol = 2, nrow = 2, byrow=T)))
    dev.off()
    write.csv2(DATA$data, file.path(arguments$wd, paste0("results_", format(Sys.time(), "%Y-%m-%d"), ".csv")), row.names = FALSE)

}

if (getOption('run.main', default=TRUE)) {
   main(args)
}