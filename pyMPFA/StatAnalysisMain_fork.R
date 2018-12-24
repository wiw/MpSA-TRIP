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
options(run.main=TRUE)
source(file.path("/home/anton/data/R-script/R-counts/RUN", "functions.R"))
args <- commandArgs(TRUE)

####################
# Replace arguments. Use this section only for development! Remove in the production.

# wd <- "/home/anton/backup/input/trip/RUN_2018-06-07/results/statistics/trip6_2"
# control <- "/home/anton/backup/input/trip/RUN_2018-06-07/results/statistics/trip6_2/control.json"
# data <- "/home/anton/backup/input/trip/RUN_2018-06-07/results/statistics/trip6_2/data.json"
# rpl_count <- "/home/anton/backup/input/trip/RUN_2018-06-07/results/statistics/trip6_2/rpl_count.json"

# use_method <- "mpsa"
# wd <- "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics/Lib_33-40"
# control <- "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics/Lib_33-40/control.json"
# data <- "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics/Lib_33-40/data.json"
# rpl_count <- "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics/Lib_33-40/rpl_count.json"
# arguments <- list("use_method"=use_method, "wd"=wd, "control"=control, "data"=data, "rpl"=rpl_count)
###################

parse_args <- function(args) {
    output_arguments = list()
    output_arguments$use_method <- args[1]
    output_arguments$wd <- args[2]
    if (args[1] == "mpsa") {
        output_arguments$control <- args[3]
        output_arguments$data <- args[4]
        output_arguments$rpl <- args[5]
    } else {
        output_arguments$data <- args[3]
        output_arguments$rpl <- args[4]
    }
    return(output_arguments)
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


convert_data <- function(loaded_jsn, first_item_name, use_method) {
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
    mtxDF <- lapply(mtx, function(x) {x1 <- type.convert(as.character(x))
                                            if(is.factor(x1))
                                                as.character(x1)
                                            else x1})
    converted_df <- cbind(data.frame("V0" = first_column, stringsAsFactors=F), mtxDF, stringsAsFactors=F)
    names(converted_df) <- c(first_item_name, sub("-", "_", names(loaded_jsn[[1]])))
    if (use_method == "trip") {
        replacement <- c("chr", "start", "string")
        names(converted_df)[names(converted_df) %in% c("mutation", NA)] <- replacement
    }
    return(converted_df)
}


main_loading <- function(arguments) {
    output <- list()
    for (item in names(arguments)[!names(arguments) %in% c("wd", "use_method")]) {
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
            output[[item]] <- convert_data(load_json(arguments[[item]]), first_item_name, arguments$use_method)
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


remove_unv_data <- function(df, col_vector) {
    all_combinations <- combn(col_vector, 2, simplify = FALSE)
    for (comb in all_combinations) {
        col1 <- comb[1]
        col2 <- comb[2]
        df <- df[!is.na(df[,col1]) & !is.na(df[,col2]),]
        df <- df[!is.infinite(df[,col1]) & !is.infinite(df[,col2]),]
    }
    return(df)
}


custom_mean <- function(df, col_vector) {
    return(rowMeans(df[, col_vector]))
}


vector.is.empty <- function(x) {
    return( length(x) == 0 )
}


prepare_data <- function(DATA) {
    for (set in sort(names(DATA)[!names(DATA) %in% "rpl"])) {
        target_vector <- grep("_", names(DATA[[set]]))
        # Count reads per millon
        for (experiment in names(DATA[[set]])[target_vector]) {
            rpm_column_name <- paste(experiment, "rpm", sep="_")
            DATA[[set]][[rpm_column_name]] <- rpm(DATA[[set]][[experiment]], DATA$rpl[[experiment]])
        }
        # Mean RPM data
        rpm_set <- names(DATA[[set]][grep("rpm", names(DATA[[set]]))])
        rpm_set_unique_names <- unique(sub("([a-z]*)_.*", "\\1", rpm_set))
        for (item in rpm_set_unique_names) {
            mean_vector <- grep(paste0(item, ".*rpm"), names(DATA[[set]]))
            mean_col_name <- paste0(item, "_mean")
            rudt <- remove_unv_data(DATA[[set]], mean_vector)
            if (ncol(rudt) != ncol(DATA[[set]])) {
                print("Founded NA&Inf values. I don't know why it happened. Take your own. Bye...")
                break
            }
            DATA[[set]][[mean_col_name]] <- custom_mean(DATA[[set]], mean_vector)
        }
        # Normalization expression data
        mean_norm_data <- names(DATA[[set]][grep(paste0("norm", ".*mean"), names(DATA[[set]]))])
        mean_expr_data <- names(DATA[[set]][grep(paste0("expr", ".*mean"), names(DATA[[set]]))])
        if (!vector.is.empty(mean_norm_data)) {
            if (!vector.is.empty(mean_expr_data)){
                column_name <- paste0(mean_expr_data, "_norm")
                DATA[[set]][[column_name]] <- normalization(DATA[[set]][[mean_expr_data]], DATA[[set]][[mean_norm_data]])
            } else {
                print("Expression data is empty. Bye...")
                break
            }
        } else {
            print("Normalization data is empty. Bye...")
            break
        }
    }
    return(DATA)
}


normalization_to_control <- function(DATA) {
    # Count mean value for control wt data
    wt_row_vector <- grep("wt", DATA$control$control)
    wt_m_n_expression <- grep("expr.*mean.*norm.*", names(DATA$control))
    wt_values <- DATA$control[wt_row_vector, wt_m_n_expression]
    wt_average <- mean(wt_values)
    # Normalization to average wt control data
    for (set in sort(names(DATA)[!names(DATA) %in% "rpl"])) {
        m_n_expression <- grep("expr.*mean.*norm.*", names( DATA[[ set ]] ))
        names_final_expression <- "expression_mn_control"
        DATA[[ set ]][[ names_final_expression ]] <- log2 ( DATA[[ set ]][[ m_n_expression ]] / wt_average )
        # Cleaning output data
        df <- DATA[[ set ]]
        DATA[[ set ]] <- df[ !is.na( df[ ,names_final_expression ] ) & !is.infinite( df[ ,names_final_expression ] ), ]
    }
    return(DATA)
}


scatter_plots_wrapper <- function(DATA, experiment_set) {
    main_vector <- list()
    i <- 1
    for (experiment in experiment_set) {
        scatter_vector <- list()
        scatter_column_number <- grep( paste0(experiment, ".*rpm"), names(DATA$data) )
        lsc <- length(scatter_column_number)
        if (lsc == 0) {
            print("I don't have any data for scatter plots computing. Bye...")
            break
        } else if (lsc == 1) {
            print("I have only one vector. I can't compare...")
            break
        } else if (lsc >= 3) {
            all_comb <- combn(scatter_column_number, 2, simplify = FALSE)
            for (comb in all_comb) {
                vector_item <- DATA$data[, comb]
                scatter_vector[[ LETTERS[i] ]] <- DATA$data[, comb]
                i <- i + 1
            }
        } else {
            scatter_vector[[ LETTERS[i] ]] <- DATA$data[, scatter_column_number]
            i <- i + 1
        }
        main_vector[[ experiment ]] <- scatter_vector
    }
    return ( scatter_plots(main_vector) )
}


scatter_plots <- function(main_vector) {
    plot.list <- list()
    for ( experiment in names(main_vector) ) {
        for (sets in names( main_vector[[ experiment ]] ) ) {
            names_V <- names( main_vector[[ experiment ]][[ sets ]] )
            V1 <- main_vector[[ experiment ]][[ sets ]][, 1]
            V2 <- main_vector[[ experiment ]][[ sets ]][, 2]
            corr_vector <- c()
            for ( cor_method in c("pearson", "spearman") ) {
                corr_vector <- append( corr_vector, paste0( cor_method, ".cor = ", round( cor( V1, V2, method = cor_method, use = "pairwise.complete.obs" ), digits=2 ) ) )
            }
            if ( length(V1) == length(V2) ) {
                all_values <- paste0( "values = ", length(V1) )
            }
            label_text <- c(sets, corr_vector, all_values)
            x.coord <- rep( max( V1, na.rm = T ) * 0.1, length(label_text) )
            y.coord <- seq( max( V2, na.rm = T ) * 0.8, by = -(3 * max( V2, na.rm = T ) / 100), length.out = length(label_text) )
            scplot <- ggplot( main_vector[[ experiment ]][[ sets ]], aes_string( names_V[1], names_V[2] ) ) + 
                geom_point(alpha = 5/10, colour = "#65A6D1", size = 3, na.rm = T) + 
                xlab( paste0( experiment, " replicate 1" ) ) + 
                ylab( paste0( experiment, " replicate 2" ) ) + 
                annotate( "text", 
                    x = x.coord, 
                    y = y.coord, 
                    label = label_text, 
                    hjust = 0,
                    size = 5, 
                    na.rm=T ) + 
                labs( title = paste0( "Scatter plot correlation\nbetween ", experiment, " replicate 1 and 2" ) ) + 
                theme_bw() + 
                theme(plot.title = element_text(size = 18),
                    axis.title.x = element_text(size = 16), 
                    axis.title.y = element_text(size = 16), 
                    axis.text.x = element_text(size = 14), 
                    axis.text.y = element_text(size = 14)
                    )
            names_of_pltlist_item <- paste(experiment, sets, sep = "_")
            plot.list[[ names_of_pltlist_item ]] <- scplot
        }
    }
    return(plot.list)
}


density_plots <- function(DATA){
    plot.list <- list()
    # Density plot for expression replicate mean
    wt <- grep("wt", DATA$control$control)
    delta_c <- grep("delta", DATA$control$control)
    expr_mean_name <- grep("_mn_control", names(DATA$control))
    control_mean <- sort(c(mean(DATA$control[[expr_mean_name]][wt]), mean(DATA$control[[expr_mean_name]][delta_c])))
    t.coord <- with(DATA$data, density(expression_mn_control))
    x.coord <- quantile(t.coord$x)[names(quantile(t.coord$x)) %in% "0%"]
    y.coord <- quantile(t.coord$x)[names(quantile(t.coord$x)) %in% "50%"]
    e_mean <- ggplot() +
        geom_density(colour="#ff9933", data=DATA$data, aes(expression_mn_control)) +
        geom_vline(aes(xintercept=control_mean), colour=c("#336600", "#330066")) +
        theme_bw() +
        xlab("Norm. expression value ( log2( expr / (norm * control) )") +
        ylab("Expression count (portion)") +
        annotate( "text", 
            x = x.coord, 
            y = y.coord, 
            hjust = 0, 
            label = paste( "Total barcodes\nexpr. mean:", nrow( DATA$data ) ) ) +
        labs(title = "Mean expression")
    plot.list$e_mean <- e_mean
    return(plot.list)
}

trip_density_plots <- function(DATA, to_title){
    plot.list <- list()
    # Density plot for expression replicate mean
    t.coord <- with(DATA$data, density(expression_mean_norm))
    x.coord <- quantile(t.coord$x)[names(quantile(t.coord$x)) %in% "0%"]
    y.coord <- quantile(t.coord$x)[names(quantile(t.coord$x)) %in% "50%"]
    e_mean <- ggplot() +
        geom_density(colour="#ff9933", data=DATA$data, aes(expression_mean_norm)) +
        theme_bw() +
        xlab("Norm. expression value ( expr / norm )") +
        ylab("Expression portion") +
        annotate( "text", 
            x = x.coord, 
            y = y.coord, 
            hjust = 0, 
            label = paste( "Total barcodes\nexpr. mean:", nrow( DATA$data ) ) ) +
        labs(title = paste0("Mean expression ", to_title))
    plot.list$e_mean <- e_mean
    return(plot.list)
}

pwm <- function(DATA, column_name) {
    plot.list <- list()
    for (method in c("bits", "prob")){
        plot.list[[method]] <- ggplot() + geom_logo(as.character(DATA[[column_name]]), method=method, seq_type="dna") + theme_logo() + labs(title=paste0("Position weight matrix by ", method, ". Count: ", length(DATA[[column_name]])))
    }
    return(plot.list)
}


main <- function(args) {
    arguments <- parse_args(args)
    # Warning! Using this commented block only for testing.
    ##########################################
    # DATA <- main_loading( arguments )
    # DATA <- prepare_data( DATA )
    # DATA <- normalization_to_control( DATA )
    ##########################################
    if (arguments$use_method == "mpsa") {
        DATA <- normalization_to_control( prepare_data( main_loading( arguments ) ) )
        filename <- basename(arguments$wd)
        draw_list <- c( scatter_plots_wrapper( DATA, c("expression", "normalization") ), density_plots(DATA), pwm(DATA$data, "mutation") )
    } else {
        DATA <- prepare_data( main_loading( arguments ) )
        print("DATA ... OK")
        filename <- paste(basename(dirname(arguments$wd)), basename(arguments$wd), sep="_")
        draw_list <- c( scatter_plots_wrapper( DATA, c("expression", "normalization") ), trip_density_plots(DATA, filename) )
    }
    if (length(draw_list) <= 4) {
        cols_value <- ceiling( length(draw_list) / 2 )
    } else {
        cols_value <- 2
    }
    # Draw plots
    pdf( file = file.path( arguments$wd, paste0( "Plots from ", filename, format(Sys.time(), " %Y-%m-%d"), ".pdf" ) ), width=14, height=14 )
    print(multiplot(plotlist = draw_list, cols = cols_value))
    dev.off()
    # Write output table to csv file
    write.csv2(DATA$data, file.path(arguments$wd, paste0("results_", format(Sys.time(), "%Y-%m-%d"), ".csv")), row.names = FALSE)
    # Write all environment to RData
    save.image(file=file.path(arguments$wd, "session.RData"))

}


if (getOption('run.main', default=TRUE)) {
   main(args)
}


##########################
# DEVELOPMENT SECTION
# aa <- DATA$data
# unique_mutation <- unique(aa$mutation)
# aa_mut <- data.frame()
# aa_mut_one <- data.frame()
# aa_mut_two <- data.frame()
# for (umt in unique_mutation){
#     selector <- aa[aa$mutation %in% umt, ]
#     len_mut <- nrow(selector)
#     if (len_mut == 2){
#         aa_mut <- rbind(aa_mut, selector)
#         aa_mut_one <- rbind(aa_mut_one, selector[1, ])
#         aa_mut_two <- rbind(aa_mut_two, selector[2, ])
#     }
# }

# scatter_plots <- function(V1, V2) {
#     plot.list <- list()
#     tmp_df <- data.frame("V1"=V1, "V2"=V2)
#     corr_vector <- c()
#     experiment <- "lib 29-36"
#     for ( cor_method in c("pearson", "spearman") ) {
#         corr_vector <- append( corr_vector, paste0( cor_method, ".cor = ", round( cor( V1, V2, method = cor_method, use = "pairwise.complete.obs" ), digits=2 ) ) )
#     }
#     if ( length(V1) == length(V2) ) {
#         all_values <- paste0( "values = ", length(V1) )
#     }
#     label_text <- c(corr_vector, all_values)
#     x.coord <- rep( max( V1, na.rm = T ) * 0.1, length(label_text) )
#     y.coord <- seq( max( V2, na.rm = T ) * 0.8, by = -(3 * max( V2, na.rm = T ) / 100), length.out = length(label_text) )
#     scplot <- ggplot( tmp_df, aes( V1, V2 )) + 
#         geom_point(alpha = 5/10, colour = "#65A6D1", size = 3, na.rm = T) + 
#         xlab( paste0( experiment, " barcode 1" ) ) + 
#         ylab( paste0( experiment, " barcode 2" ) ) + 
#         annotate( "text", 
#             x = x.coord, 
#             y = y.coord, 
#             label = label_text, 
#             hjust = 0,
#             size = 5, 
#             na.rm=T ) + 
#         labs( title = paste0( "Scatter plot correlation\nbetween ", experiment, " replicate 1 and 2" ) ) + 
#         theme_bw() + 
#         theme(plot.title = element_text(size = 18),
#             axis.title.x = element_text(size = 16), 
#             axis.title.y = element_text(size = 16), 
#             axis.text.x = element_text(size = 14), 
#             axis.text.y = element_text(size = 14)
#             )
#     plot.list[[ experiment ]] <- scplot
#     return(plot.list)
# }

# plots <- scatter_plots(aa_mut_one$expression_mn_control, aa_mut_two$expression_mn_control)
# pdf( file = file.path( arguments$wd, paste0( "Experimental plot ", filename, format(Sys.time(), " %Y-%m-%d"), ".pdf" ) ), width=14, height=14 )
# print(plots)
# dev.off()