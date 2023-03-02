add_dm_annotation <- function(clust_label, pl, txt_size = 24){
    ## A list of annotations we wish to add
    ## Identify these by their clust labels, e.g., C7X, C11Y etc.
    ##
    text_annots <- list(
        "C7X" = list(geom = "text", label = "TCT Promoter", x = 59, y = 1.5),
        "C11Y" = list(geom = "text", label = "TCT Promoter", x = 59, y = 1.5),
        "C10Z" = list(geom = "text", label = "TCT Promoter", x = 59, y = 1.5),
        ##
        "C6X" = list(geom = "text", label = "DPE appears", x = 65, y = 1.5),
        "C1Y" = list(geom = "text", label = "DPE", x = 68, y = 1.5),
        "C8Y" = list(geom = "text", label = "DPE", x = 68, y = 1.5),
        "C2Z" = list(geom = "text", label = "Two DPE variants merged", x = 60, y = 1.5),
        "C3Z" = list(geom = "text", label = "TATA-box appears", x = 30, y = 1.2),
        "C4Z" = list(geom = "text", label = "DPE", x = 68, y = 1.2),
        "C15X" = list(geom = "text", label = "Polythymine stretch", x = 68, y = 1.2),
        "C14Y" = list(geom = "text", label = "Polythymine stretch", x = 68, y = 1.2),
        "C15Z" = list(geom = "text", label = "Polythymine stretch", x = 68, y = 1.2),
        "C16X" = list(geom = "text", label = c("Ohler1", "Ohler6"), x = c(59, 30), y = c(1.2)),
        "C16Y" = list(geom = "text", label = c("Ohler1", "Ohler6"), x = c(59, 30), y = c(1.2)),
        "C18Z" = list(geom = "text", label = c("Ohler1", "Ohler6"), x = c(59, 30), y = c(1.2))#,

        ## histone genes
        # "C1X" = list(geom = "text", label = "Histone genes cluster", x = 4, y = 1.7)
    )
    ##

    hist_annots <- c(paste0("C", 1:5, "X"), paste0("C", c(2:4, 6:7), "Y"), paste0("C", c(1, 5:7), "Z"))
    names(hist_annots) <- hist_annots

    ##

    fill_col <- 'white'
    rect_col <- 'black'

        rect_annots <- list(
            "C6X" = list(geom = "rect", xmin = 71, xmax = 77,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col),

            "C7X" = list(geom = "rect", xmin = 41, xmax = 51,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col,
                linewidth = 2),

            "C15X" = list(geom = "rect", xmin = 76, xmax = 91,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col),

            "C16X" = list(geom = "rect", xmin = c(40, 17), xmax = c(52, 27),
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col,
                linewidth = 2),

            "C1Y" = list(geom = "rect", xmin = 71, xmax = 81,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col),

            "C8Y" = list(geom = "rect", xmin = 71, xmax = 77,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col),

            "C11Y" = list(geom = "rect", xmin = 41, xmax = 51,
                ymin = -0.05, ymax = 2.0,
                alpha = .1, col=rect_col, fill=fill_col),

            "C14Y" = list(geom = "rect", xmin = 76, xmax = 91,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col),

            "C16Y" = list(geom = "rect", xmin = c(40, 17), xmax = c(52, 27),
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col,
                linewidth = 2),

            "C2Z" = list(geom = "rect", xmin = 71, xmax = 79,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col),

            "C3Z" = list(geom = "rect", xmin = 14, xmax = 22,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col),

            "C4Z" = list(geom = "rect", xmin = 71, xmax = 79,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col),

            "C10Z" = list(geom = "rect", xmin = 41, xmax = 51,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col),

            "C15Z" = list(geom = "rect", xmin = 76, xmax = 91,
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col),

            "C18Z" = list(geom = "rect", xmin = c(40, 17), xmax = c(52, 27),
                ymin = -0.05, ymax = 2.0,
                alpha = .01, col=rect_col, fill=fill_col,
                linewidth = 2)#,

            # "C1X" = list(geom = "rect", xmin = 1, xmax = 10,
            #     ymin = 1.5, ymax = 2.0,
            #     alpha = 1, col=rect_col, fill=fill_col,
            #     linewidth = 2)


        )

        ## text annotations
        any_match <- match(clust_label, names(text_annots))
        if(!is.na(any_match)){
            message("Added annotation for ", clust_label)
            pl <- pl + ggplot2::annotate(geom = text_annots[[clust_label]]$geom,
                x = text_annots[[clust_label]]$x,
                y = text_annots[[clust_label]]$y,
                label = text_annots[[clust_label]]$label,
                size = 10 + txt_size/.pt)
        }

        ## rect annotations
        any_match <- match(clust_label, names(rect_annots))
        if(!is.na(any_match)){
            message("Added annotation for ", clust_label)
            pl <- pl + ggplot2::annotate(geom = rect_annots[[clust_label]]$geom,
                xmin = rect_annots[[clust_label]]$xmin,
                xmax = rect_annots[[clust_label]]$xmax,
                ymin = rect_annots[[clust_label]]$ymin,
                ymax = rect_annots[[clust_label]]$ymax,
                alpha = rect_annots[[clust_label]]$alpha,
                col = rect_annots[[clust_label]]$col,
                fill = rect_annots[[clust_label]]$fill,
            )
        }
        ## histone genes clusters
        any_match <- match(clust_label, names(hist_annots))
        if(!is.na(any_match)){
            print("******")
            message("Added annotation for ", clust_label)
            pl <- pl + ggplot2::theme(panel.background = element_rect(fill = "grey85"))
        }


        return(pl)

}
## =============================================================================
