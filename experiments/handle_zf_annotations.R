add_zf_annotation <- function(clust_label, pl, txt_size = 24){
    ## A list of annotations we wish to add
    ## Identify these by their clust labels, e.g., C7X, C11Y etc.
    ##
    text_annots <- list(
        "C1X" = list(geom = "text", label = "W-box", x = 33, y = 0.5)
    )
    ##
    fill_col <- NA
    rect_col <- 'black'

    rect_annots <- list(
        "C1X" = list(geom = "rect", xmin = 13, xmax = 24,
            ymin = -0.05, ymax = 0.9,
            alpha = .01, col=rect_col, fill=fill_col,
            linewidth = 2),

        "C2X" = list(geom = "rect", xmin = 13, xmax = 24,
            ymin = -0.05, ymax = 0.6,
            alpha = .01, col=rect_col, fill=fill_col),

        "C3X" = list(geom = "rect", xmin = 13, xmax = 24,
            ymin = -0.05, ymax = 0.6,
            alpha = .01, col=rect_col, fill=fill_col),

        "C4X" = list(geom = "rect", xmin = 13, xmax = 24,
            ymin = -0.05, ymax = 0.3,
            alpha = .1, col=rect_col, fill=fill_col),

        "C5X" = list(geom = "rect", xmin = 13, xmax = 24,
            ymin = -0.05, ymax = 0.4,
            alpha = .01, col=rect_col, fill=fill_col),

        "C6X" = list(geom = "rect", xmin = 13, xmax = 24,
            ymin = -0.05, ymax = 0.4,
            alpha = .01, col=rect_col, fill=fill_col),

        "C10Z" = list(geom = "rect", xmin = 91, xmax = 196,
            ymin = -0.05, ymax = 0.4,
            alpha = .01, col=rect_col, fill=fill_col)

    )

    # segment_annots <- list(
    #     "C10Z" = list(geom = "segment", x = 91, xend=-2, y=0, yend=-10, size=2)
    # )

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

    # ## segment annotations
    # any_match <- match(clust_label, names(segment_annots))
    # if(!is.na(any_match)){
    #     message("Added annotation for ", clust_label)
    #     pl <- pl + ggplot2::annotate(geom = rect_annots[[clust_label]]$geom,
    #         xmin = rect_annots[[clust_label]]$xmin,
    #         xmax = rect_annots[[clust_label]]$xmax,
    #         ymin = rect_annots[[clust_label]]$ymin,
    #         ymax = rect_annots[[clust_label]]$ymax,
    #         alpha = rect_annots[[clust_label]]$alpha,
    #         col = rect_annots[[clust_label]]$col,
    #         fill = rect_annots[[clust_label]]$fill,
    #     )
    # }


    return(pl)

}
## =============================================================================
