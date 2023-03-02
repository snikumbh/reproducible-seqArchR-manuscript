## Helper functions used in Rmd file to geenrate plots for analyses of archR
## results from any organism


plot_dend_arch <- function(arch_plot, fname, use_ht = 40,
                           use_cutk = 2, use_cuth = NULL, clusts,
                           clust_assignment = NULL,
                           new_clusts = NULL, rawSeqs = NULL, pos_lab = NULL,
                           plot_png = TRUE, ...
                           ){


    dend_pl2 <- factoextra::fviz_dend(clusts, horiz = TRUE, main = "",
                              k = use_cutk,
                              h = use_cuth,
                              ...
    )

    dend_pl2 <- dend_pl2 + ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(add = c(0.6,0.2))) +
        NULL

    if(is.null(new_clusts)){
        sam_foo2 <- cowplot::plot_grid(dend_pl2, arch_plot,
                                       ncol = 2,
                                       rel_widths = c(0.40,1), rel_heights = c(1,1),
                                       align = "hv")
        ## pdf
        cowplot::ggsave2(paste0(fname, ".pdf"), plot = sam_foo2, width = 30, height = use_ht,
                         units = "cm", dpi = 300)
        ## png
        if(plot_png){
            cowplot::ggsave2(paste0(fname, ".png"), plot = sam_foo2, width = 40, height = use_ht,
                         units = "cm", dpi = 300)
        }
    }else{
        ## new_clusts info is provided
        new_arch_pl <- archR::plot_arch_for_clusters(
            seqs = rawSeqs,
            clust_list = new_clusts, pos_lab = pos_lab,
            xt_freq = 5, set_titles = FALSE, show = FALSE, bits_yax = "auto")

       new_arch_pl2 <- lapply(seq_along(new_arch_pl), function(x){
            pl <- new_arch_pl[[x]] +
                ggtitle(paste0("Cluster C", x ," Obtained by combining: ",
                                       paste(clust_assignment[[x]], collapse= ", "))) +
                ggplot2::theme(axis.text = ggplot2::element_text(size = 0),
                               axis.text.x = ggplot2::element_text(
                                   angle = 0, vjust = 2, hjust = 0.5),
                               axis.text.y = ggplot2::element_text(vjust = 0.5),
                               axis.title.y = ggplot2::element_text(size = 0),
                               axis.ticks.length = ggplot2::unit(0.00, "cm"),
                               plot.title = element_text(margin=margin(1,0,0,0)),
                               plot.margin = ggplot2::unit(c(-0.1,0,-0.4,-0.4), "cm"))
        })
        sam_foo <- cowplot::plot_grid(plotlist = new_arch_pl2, ncol = 1)


        ######
        sam_foo2 <- cowplot::plot_grid(dend_pl2, arch_plot, sam_foo,
                                       ncol = 3,
                                       rel_widths = c(0.40,1,1), rel_heights = c(1,1,1),
                                       align = "hv")
        ## pdf
        cowplot::ggsave2(paste0(fname, ".pdf"), plot = sam_foo2, width = 50, height = use_ht,
                         units = "cm", dpi = 300)
        ## png
        if(plot_png){
            cowplot::ggsave2(paste0(fname, ".png"), plot = sam_foo2, width = 60, height = use_ht,
                             units = "cm", dpi = 300)
        }
    }

    return(sam_foo2)

}
## =============================================================================


## Get scores of how much the sequence matches the motif
simply_scores_mat <- function(seqs, use_unif_bg = TRUE, logodds = TRUE){
    # print(names(seqs[1:10]))
    sam_pfm <- Biostrings::consensusMatrix(seqs)[1:4,] ## select rows A,C,G,T
    # PFM to PWM
    if(use_unif_bg){
        sam_pwm <- PWMEnrich::toPWM(sam_pfm, seq.count = length(seqs))
    }else{
        prior <- getBackgroundFrequencies("dm6", quick = TRUE)
        sam_pwm <- PWMEnrich::toPWM(sam_pfm, seq.count = length(seqs),
                                    prior = prior)
    }
    #
    # sam_scores <- PWMEnrich::motifScores(sequences = seqs,
    #   motifs = sam_pwm, raw.scores = raw)
    sam_scores <- unlist(lapply(seq(length(seqs)), function(x){
        score <- Biostrings::PWMscoreStartingAt(sam_pwm$pwm,
                                                subject=as.character(seqs[x]),
                                                starting.at = 1)
        score <- ifelse(logodds, 2^score, score)
        score
    }))
    sam_scores
}
## =============================================================================


form_str <- function(str){
    return(format(str, trim = TRUE, digits = 3, nsmall = 2, scientific = TRUE))
}
## =============================================================================


get_tissue_spec_scores <- function(df, tissueSpec_values, ensId){#archR_clusts, sn, clust_id){
    # df_obj is the peakAnno as.data.frame

    ##
    idx1 <- which(tissueSpec_values$ensembl_gene_id == ensId)
    subset(tissueSpec_values, ensembl_gene_id == ensId)
    tissueSpec_values$aucRatio[idx1]
}
## =============================================================================


get_motif_scores <- function(result_obj, archR_clusts, sn, clust_id, unif_bg = TRUE, logodds = TRUE){

    # get sequences from archR result object for only the sequences in the cluster
    seqs <- result_obj[[sn]]$rawSeqs[archR_clusts[[sn]][clust_id][[1]]]
    # fetch scores
    motif_scores <- simply_scores_mat(seqs, use_unif_bg = unif_bg, logodds = logodds)

    motif_scores
}
## =============================================================================

## Use this function to get n colors (in sequence) from the specified palette
## -- Can specify n greater than that existing in a palette, in which case
## the colors are recycled and returned
## -- Can also specify a list of colors of choice < n, which are recycled and
## returned. Useful when random colors from a palette are required
get_ncolors <- function(n, palname = "Set1", clrs = NULL){

    ## Recycle color vector from RColorBrewer
    use_colors <- clrs
    if(is.null(clrs)){
        use_colors <- suppressWarnings(
            RColorBrewer::brewer.pal(n = n, name = palname))
    }

    nColor <- length(use_colors)
    if(n <= nColor){
        n_colors <- use_colors[seq_len(n)]
        return(n_colors)
    }

    rep_times <- base::ceiling((n-nColor)/nColor)
    if(n %% nColor == 0) rep_times <- rep_times + 1
    additional <- ((n-nColor) %% nColor)
    col_idx <- c(rep(seq_len(nColor), rep_times), seq_len(additional))
    n_colors <- use_colors[col_idx]
    n_colors
}
## =============================================================================


get_hline_ycoord <- function(res_clusts){
    length(unlist(res_clusts)) - cumsum( unlist(lapply(res_clusts, length)) )
}
## =============================================================================

## check_and_create_dir
check_and_create_dir <- function(dir_path){
    creation_ok <- FALSE
    if(!dir.exists(dir_path)){
        message("Creating directory: ", dir_path)
        creation_ok <- dir.create(dir_path)
    }else{
        message("Directory exists: ", dir_path)
        creation_ok <- TRUE
    }
    ## TRUE: success in creation; FALSE: otherwise
    return(creation_ok)
}
## =============================================================================


get_go_term_dot_plot <- function(peakAnno, useOrgDb, useKeyType,
                                 choose_idx, plot_title_text,
                                 bar_or_dot = "dot", font.size = 10){

    ##
    samarth_go <- tryCatch(enrichGO(gene =
                                as.data.frame(peakAnno)[choose_idx, useKeyType],
                                keyType = useKeyType,
                                OrgDb = useOrgDb, ont = "all"),
                           error = function(e) e)

    if(class(samarth_go)[1] == "simpleError"){
        print("ontology ALL errored, hence using BP")
        samarth_go <- tryCatch(enrichGO(gene =
                                as.data.frame(peakAnno)[choose_idx, useKeyType],
                                keyType = useKeyType,
                                OrgDb = useOrgDb, ont = "BP"),
                               error = function(e) e)
    }
    use_label <- "N/A"
    if(class(samarth_go)[1] == "simpleError"){
        print("ontology ALL and BP both errored, hence N/A")
        samarth_go <- NULL
        use_label <- "All and BP errored, hence N/A"
    }
    ##
    if(!is.null(samarth_go) && nrow(samarth_go) > 0){
        if(bar_or_dot == "dot"){
            p1 <- dotplot(samarth_go,
                          title = plot_title_text,
                          font.size = font.size
            )
        }else{
            p1 <- barplot(samarth_go,
                          title = plot_title_text,
                          font.size = font.size
            )
        }
        # p1 <- p1 + ggplot2::theme(legend.text = ggplot2::element_text(size=15),
        #                           axis.text.x = ggplot2::element_text(angle = 90))
        p1 <- p1 + DOSE::theme_dose(font.size=font.size) +
            ggplot2::scale_y_discrete(labels = scales::label_wrap(30)) +
            ggplot2::theme(title = element_text(size=18),
                           legend.text = element_text(size=14),
                           legend.title = element_text(size=14))
        return(p1)
    }else{
        print("Null/Zero rows")
        p1 <- ggplot() +
            geom_blank() + theme_bw() +
            theme(panel.grid = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  title = element_text(size=18)) +
            geom_text(aes(0,0,label=use_label)) +
            ggplot2::labs(title = plot_title_text)
        return(p1)
    }
}
## =============================================================================

get_fixed_anno_ord <- function(){

    anno_terms_ord <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron",
                paste0(rep("Downstream (", 3), c("<1", "1-2", "2-3"), "kb)"),
                    "Distal Intergenic", "NA")

    anno_terms_ord
}
## =============================================================================


get_named_colors <- function(anno_terms_ord, palname = "Set1"){
    use_colors <- get_ncolors(n=length(anno_terms_ord), palname=palname)
    if(palname=="Set1") use_colors <- use_colors[c(2:length(use_colors),1)]
    names(use_colors) <- anno_terms_ord
    use_colors
}
## =============================================================================


get_go_term_freq_plot <- function(peakAnno, choose_idx, plot_title_text,
                                  palname = "Set1", font.size = 10){
    ##
    anno_terms_ord <- get_fixed_anno_ord()
    use_colors <- get_named_colors(anno_terms_ord, palname)
    # use_colors <- get_ncolors(n=length(anno_terms_ord), palname=palname)
    # if(palname=="Set1") use_colors <- use_colors[c(2:length(use_colors),1)]
    # names(use_colors) <- anno_terms_ord
    ##
    anno_terms1 <- as.data.frame(peakAnno)[choose_idx,]$annotation
    anno_terms_splits <- strsplit(anno_terms1, " ")
    anno_terms <- unlist(lapply(anno_terms_splits, function(x){
        if(is.na(x[1])){ "NA";
        }else if(grepl("Intron", x[1]) || grepl("Exon", x[1])){
            if(x[4] == 1){
                paste("1st", x[1])
            }else{
                paste("Other", x[1])
            }
        }else if(grepl("Promoter", x[1])){
            x[1]
        }else{ paste(x[1], x[2]) }
    }))
    ##
    wordFreqTable <- table(anno_terms)
    ##
    wordFreqDF <- data.frame(word = rownames(wordFreqTable),
                             freq = as.vector(wordFreqTable))
    ## convert to tibble using tidyr::pivot_longer
    wordFreqTib <- tidyr::pivot_longer(wordFreqDF, cols = c("word"),
                                       values_to = "Annotation")
    print(wordFreqTib$Annotation)
    ## ggplot2 frequency barplot
    # breaks <- sort(c(pretty(range(wordFreqDF$freq)), max(wordFreqDF$freq)))
    # max_idx <- which(breaks == max(wordFreqDF$freq))
    # colvec <- rep("black", length(breaks))
    # colvec[max_idx] <- "red"
    p2 <- ggplot(wordFreqTib, aes(y = freq, x = name))+
        ggplot2::theme_bw() +
        ggplot2::geom_col(aes(fill = Annotation), width=0.98) +
        ggplot2::scale_fill_manual(values=use_colors) +
        # ggplot2::scale_y_continuous(breaks = breaks,
        #                    minor_breaks = breaks,
        #                    labels = breaks) + #,
        #labels = scales::percent) +
        ggplot2::theme(#axis.title = element_text(size = 14),
            #axis.text = element_text(size = 12),
            axis.title = element_text(size=font.size),
            axis.text = element_text(size=font.size),
            axis.text.y = element_text(hjust = 1.0),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank(),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14),
            legend.position = "top",
            legend.direction = "vertical") +
        ggplot2::labs(y = "Count") +
        # ggplot2::geom_hline(yintercept = wordFreqDF[which.max(wordFreqDF[,"freq"]),"freq"],
        #            color = "red", linetype = "dashed") +
        # ggplot2::guides(fill = FALSE) +
        NULL

    p2
}
## =============================================================================

get_iqw_ord_plot <- function(iqw = FALSE, tpm = FALSE, phast = FALSE,
                             use_notch = TRUE, y_axis_text = TRUE, text_size = 12,
                             samarth_df, seqs_clust, use_suffix = "X"){

    clust_lens <- unlist(lapply(seqs_clust, length))
    clust_labels <- paste("C", 1:length(seqs_clust), use_suffix,
                          paste0(" (", clust_lens, ")"), sep="")
    clr <- RColorBrewer::brewer.pal(3, "Dark2")
    ##
    if(iqw){
        pl <- ggplot(samarth_df,
                     aes(y=fct_reorder(clust_ID, IQW,
                                       .fun = median,
                                       .desc = TRUE),
                         x=IQW)) +
            ggplot2::geom_boxplot(outlier.size = 1, width = 0.5, notch = use_notch,
                color = "black", fill = clr[1]) +
            xlab("IQW")
        pl <- pl + ggplot2::annotation_logticks(sides="b",
                                                short = unit(0.5, "cm"),
                                                mid = unit(0.6, "cm"),
                                                long = unit(0.7, "cm"))
    }
    ##
    if(tpm){
        pl <- ggplot(samarth_df,
                     aes(y=fct_reorder(clust_ID, IQW,
                                       .fun = median,
                                       .desc = TRUE),
                         x=domTPM)) +
            geom_boxplot(outlier.size = 1, width = 0.5, notch = use_notch,
                         color = "black", fill = clr[2]) +
            # theme_linedraw() +
            scale_x_log10(#trans = 'log10',
                          breaks = scales::trans_format("log10", function(x) 10^x),
                          labels = scales::trans_format("log10", scales::math_format(10^.x))) +
            xlab("TPM") +
            ylab("Clusters") +
            scale_y_discrete(labels = rev(clust_labels))
        pl <- pl + ggplot2::annotation_logticks(sides="b",
                                                short = unit(0.5, "cm"),
                                                mid = unit(0.6, "cm"),
                                                long = unit(0.7, "cm"))
    }
    if(phast){
        pl <- ggplot(samarth_df,
                     aes(y=fct_reorder(clust_ID, IQW,
                                       .fun = median,
                                       .desc = TRUE),
                         x=phast)) +
            geom_boxplot(outlier.size = 1, width = 0.5, notch = use_notch,
                         color = "black", fill = clr[3]) +
            xlab("PhastCons score") +
            ylab("Clusters") +
            scale_y_discrete(labels = rev(clust_labels)) +
            theme_linedraw()
    }
    ##
    if(!phast) pl <- pl + scale_x_log10()
    pl <- pl  +
        ylab("Clusters") +
        scale_y_discrete(labels = rev(clust_labels)) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.title = element_text(colour = "black", size = text_size),
            axis.text.x = element_text(colour = "black", angle = 45,
                                       size = text_size, vjust = 0.5, hjust = 0.5),
            axis.text.y = element_text(colour = "black",
                                       size = text_size, hjust =1)
        )
    ##

    if(!y_axis_text){
        pl <- pl  +
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),)
    }
    return(pl)
    ####
}
## =============================================================================

write_as_track_bed <- function(given_df, seq_ids, track_name, bedFilename){
    track.description <- track_name
    message(bedFilename)
    write(paste('track name="', track_name,'" description="',
                track.description,'" visibility="pack"', ' itemRgb="On"', sep = ''),
          file = bedFilename, append = F)
    use_score <- ifelse("phast" %in% colnames(given_df),
                        given_df$phast, given_df$domTPM)
    write.table(data.frame(given_df$chr,
                           formatC(given_df$start, format = 'f', digits = 0),
                           formatC(given_df$end, format = 'f', digits = 0),
                           seq_ids,
                           score = use_score,
                           given_df$strand
    ),
    file = bedFilename,
    append = T, col.names = F, row.names = F,
    quote = F, sep = '\t')
}
## =============================================================================

get_strand_plot_title <- function(this_id, nclust, clust_names, this_n,
                                  tot_n, strand_val = "+"){
    if(strand_val == "+"){
        title_str <-  paste0("(", this_id , "/", nclust,
                         ") Arch `", clust_names[this_id], "': ",
                         this_n, " (+ strand) /", tot_n)
    }else if(strand_val == "-"){
        title_str <-  paste0("(", this_id , "/", nclust,
                             ") Arch `", clust_names[this_id], "': ",
                             this_n, " (- strand) /", tot_n)
    }

    title_str
}
## =============================================================================

get_strand_specific_indices <- function(df, seq_ids_in_clust, strand_val = "+"){
    return(seq_ids_in_clust[which(df$strand[seq_ids_in_clust] == strand_val)])
}
## =============================================================================




