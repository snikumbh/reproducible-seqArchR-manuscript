---
author: "Sarvesh Nikumbh"
date: '24 February, 2023'
title: "archR result analysis report for CTSS in Human (inr KO 5_5)"
output: 
  bookdown::html_document2:
    toc: true
    toc_float: true
    number_sections: true
    anchor_sections: true
    number_tables: true
  editor_options:
    chunk_output_type: console
---





This CAGE data was obtained from ENCODE and processed by merging all cell lines 
data.

# Sample: `human_cellGroup_merged`

This file produces additional analyses plots for sample "human_cellGroup_merged" 
for condition Inr knock out [-5,+5] bp.








## Plots for sections marked TRUE are included in this report.





```r

sample_names <- c("human_cellGroup_merged")

archR_human_best_run <- vector("list", length(sample_names))
names(archR_human_best_run) <- sample_names

human_result <- vector("list", length(sample_names))
names(human_result) <- sample_names

human_bed_info <- vector("list", length(sample_names))
names(human_bed_info) <- sample_names

seqs_clusters_as_list <- vector("list", length(sample_names))
names(seqs_clusters_as_list) <- sample_names

seqs_clusters_as_list_ordered <- vector("list", length(sample_names))
names(seqs_clusters_as_list_ordered) <- sample_names

perSample_CAGEobj <- vector("list", length(sample_names))
names(perSample_CAGEobj) <- sample_names

perSample_peakAnno <- vector("list", length(sample_names))
names(perSample_peakAnno) <- sample_names

perSample_tissueSpec_pl <- vector("list", length(sample_names))
names(perSample_tissueSpec_pl) <- sample_names

samarth_df <- vector("list", length(sample_names))
names(samarth_df) <- sample_names
##


##
## lists storing plots
## 
## store iqw+tpm plots
perSample_pl <- vector("list", length(sample_names)) 
names(perSample_pl) <- sample_names


## store result_directory path
result_dir_path <- vector("list", length(sample_names)) 
names(result_dir_path) <- sample_names



## store clusters
perSample_archR_clusts <- vector("list", length(sample_names))
names(perSample_archR_clusts) <- sample_names


################################################################################
## Per section some variables need to be set. These are done here

# txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
## CAGEr object processing


# read TagClusters information from corresponding RDS file for sample
# All samples in one RDS object
cager_obj <- file.path(archR_org_data_path, 
                       paste0("samarth_hsapiens_TC_sample_human_cellGroup_merged_",
                       "minTPM", use_minTPM,
                       "_flank_up500_flank_down500_all70.rds"))

ok_chr_names <- paste0("chr", c(1:21, "X", "Y"))

for(sn in sample_names){
    result_dir_path[[sn]] <- file.path(archR_org_results_path, 
                                       paste0(sn, "_inrKO_5_5_results"))
}

```


```r



for(sn in sample_names){
    archR_human_best_run[[sn]] <- 
    file.path(archR_org_results_path, 
        paste0("archR_result_hsapiens_", sn,
               "_flank_up50_flank_down150__inrKO_5_5",
               "_modSelType_stability_chunkSize_5000", 
               "_bound_1e-06_collate_FTTTF"))
    if(file.exists(archR_human_best_run[[sn]])){
    human_result[[sn]] <- readRDS(file.path(archR_human_best_run[[sn]],
                                          "archRresult.rds"))
    }else{
    stop("Check if file exists. ", archR_human_best_run[[sn]])
    }
    ## Bed Info -- IQW and dominant TPM values 
    bed_fname <- file.path(archR_org_data_path, paste0("samarth_hsapiens_TC_sample_", sn, 
                                                   paste0("_minTPM", use_minTPM, "_all70"),
                                                   ".bed"))
    print(bed_fname)
    
    human_bed_info[[sn]] <- read.delim(file = bed_fname,
                          sep = "\t", header = TRUE)
    colnames(human_bed_info[[sn]]) <- c("chr", "start", "end", 
                                                    "IQW", "domTPM", "strand")
    ##
    message("Size of bed_info: ", nrow(human_bed_info[[sn]]))
}
## [1] "experiments/data/human/samarth_hsapiens_TC_sample_human_cellGroup_merged_minTPM1_all70.bed"
## Size of bed_info: 9523
```

For annotation, TSS region considered is -500 upstream and 100 nt downstream.


```
## SAMPLE: human_cellGroup_merged
## Directory exists: experiments/results/human/human_cellGroup_merged_inrKO_5_5_results
## ℹ Using saved GRanges obj from CAGEr object
## >> preparing features information...		 2023-02-24 14:54:13 
## >> identifying nearest features...		 2023-02-24 14:54:13 
## >> calculating distance from peak to TSS...	 2023-02-24 14:54:14 
## >> assigning genomic annotation...		 2023-02-24 14:54:14 
## >> adding gene annotation...			 2023-02-24 14:54:33
## 'select()' returned 1:many mapping between keys and columns
## >> assigning chromosome lengths			 2023-02-24 14:54:33 
## >> done...					 2023-02-24 14:54:33
```

## Order archR clusters by median IQ width {#prepare-objects}

By default, the clusters are ordered by their _median_ interquantile widths. 
One can choose to order by the _mean_ interquantile widths by setting paramater _iqw_order_by_median_ to FALSE.


### ENCODE CAGE in Hsapiens cellGroup merged


```r


## Sample 1, human_cellGroup_merged
sn <- 1
itr <- 5
use_aggl <- 'ward.D'
use_dist <- 'cor'

iter5_clusts_reord <- archR::collate_archR_result(
  result = human_result[[sn]],
  iter = itr, clust_method = 'hc', aggl_method = use_aggl,
  dist_method = use_dist, regularize = TRUE, topn = 50,
  flag = list(debugFlag = FALSE, verboseFlag = TRUE), collate = FALSE,
  return_order = TRUE)


## these are in default archR ordering
clust_archR_ord_list <- archR::get_seqs_clust_list(
  seqs_clust_lab = human_result[[ sn ]]$seqsClustLabels[[itr]])

## these are now ordered by the hc ordering
clust_hc_ord_list <- lapply(iter5_clusts_reord$order, function(x){
  clust_archR_ord_list[[x]]
})


ordered_arch_pl <- archR::plot_arch_for_clusters(
  seqs = human_result[[ sn ]]$rawSeqs,
  clust_list = clust_hc_ord_list, pos_lab = -50:150,
  xt_freq = 5, set_titles = FALSE, show = FALSE, bits_yax = "auto")

ordered_arch_pl2 <- lapply(rev(ordered_arch_pl), function(pl){
  pl <- pl +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 0),
                       axis.text.x = ggplot2::element_text(
                         angle = 0, vjust = 2, hjust = 0.5),
                       axis.text.y = ggplot2::element_text(vjust = 0.5),
                       axis.title.y = ggplot2::element_text(size = 0),
                       axis.ticks.length = ggplot2::unit(0.00, "cm"),
                       plot.margin = ggplot2::unit(c(-0.1,0,-0.4,-0.4), "cm"))
})
sam_foo <- cowplot::plot_grid(plotlist = ordered_arch_pl2, ncol = 1)

use_cutk <- 16

stopifnot(check_and_create_dir(result_dir_path[[sn]]))
## Directory exists: experiments/results/human/human_cellGroup_merged_inrKO_5_5_results
fname <- file.path(result_dir_path[[sn]], paste0(sample_names[sn], "_dend_arch_list_",
                                           use_dist, "_", use_aggl, "_",
                                           use_cutk, "clusters"))

sam_foo2 <- plot_dend_arch(arch_plot = sam_foo, fname = fname, use_cutk = use_cutk,
              clusts = iter5_clusts_reord, use_ht = 60, plot_png = FALSE,
              lwd = 0.4, repel = TRUE, show_labels = TRUE, labels_track_height = 0.25,
              rect = TRUE, rect_fill = TRUE,
              color_labels_by_k = TRUE)

## IMPORTANT: Seems like no additional manipulation is required. 

# 
# 
temp_clusts <- cutree(iter5_clusts_reord, k = use_cutk)
names(temp_clusts) <- NULL
## Make further few clusters
# nCl <- length(unique(temp_clusts))
# 
# # TATA architectures
# temp_clusts[c(34)] <- nCl + 1
# temp_clusts[c(33)] <- nCl + 2
# temp_clusts[c(35)] <- nCl + 3
# temp_clusts[c(26)] <- nCl + 4
# 
# # TCT architecture
# temp_clusts[c(37)] <- nCl + 5
# 
# temp_clusts[c(31)] <- nCl + 6
# temp_clusts[c(36)] <- nCl + 6
# 
# temp_clusts[c(17)] <- temp_clusts[c(13)]

# Re-plot with proper coloring that shows the manipulations in clusters
use_cutk <- 16

stopifnot(check_and_create_dir(result_dir_path[[sn]]))
## Directory exists: experiments/results/human/human_cellGroup_merged_inrKO_5_5_results
fname2 <- file.path(result_dir_path[[sn]], paste0(sample_names[sn], "_dend_arch_list_",
                                           use_dist, "_", use_aggl, "_",
                                           use_cutk, "clusters_final"))


## Also need to show alongside, how the final clusters' seqlogos look
clust_list <- lapply(unique(temp_clusts), function(x){which(temp_clusts == x)})
seqs_clusters_as_list[[sample_names[sn]]] <- archR::collate_clusters(clust_list, archR::get_seqs_clust_list(human_result[[sn]]$seqsClustLabels[[itr]]))
##

use_color <- scales::hue_pal()(length(unique(temp_clusts)))
sam_foo2 <- plot_dend_arch(arch_plot = sam_foo, fname = fname2, use_ht = 60,
               use_cutk = use_cutk,#length(unique(temp_clusts)), 
               clusts = iter5_clusts_reord, rect = TRUE, rect_fill = TRUE,
               label_cols = use_color[temp_clusts[iter5_clusts_reord$order]],
               k_colors = use_color, 
               clust_assignment = clust_list,
               new_clusts = seqs_clusters_as_list[[sample_names[sn]]],
               rawSeqs = human_result[[sample_names[sn]]]$rawSeqs,
               palette = FALSE, plot_png = FALSE)


cluster_medians_IQW <- unlist(lapply(seqs_clusters_as_list[[sn]], function(x){
    median(human_bed_info[[sn]]$IQW[x])
}))
cluster_means_IQW <- unlist(lapply(seqs_clusters_as_list[[sn]], function(x){
    mean(human_bed_info[[sn]]$IQW[x])
}))

if(iqw_order_by_median){
    ascending_order_IQW <- sort(cluster_medians_IQW, decreasing = FALSE,
                            index.return = TRUE)
}else{
    ascending_order_IQW <- sort(cluster_means_IQW, decreasing = FALSE,
                            index.return = TRUE)
}
##
seqs_clusters_as_list_ordered[[sn]] <-
                      lapply(ascending_order_IQW$ix,
                              function(x){
                                  seqs_clusters_as_list[[sn]][[x]]
                              })

perSample_archR_clusts[[sn]] <- seqs_clusters_as_list_ordered[[sn]]

```




```r

## Using Tau as a measure of tissue specificity
## This measure is defined/used in this paper:
##  https://www.aging-us.com/article/202648/text
## Various tissue specificity measures have been reviewed in this paper:
## https://doi.org/10.1093/bib/bbw008
## The above benchmarking paper finds Tau to be the most robust and sensitive 
## measure of tissue specificity. 
## 
tissueSpecTau <- read_csv(file.path(data_path, "tissueSpecificity-tau", 
                                    "Tau_score", "Tau_gene_V8.csv"), 
                          col_names = TRUE, show_col_types = FALSE)


tau_dist <- ggpubr::gghistogram(tissueSpecTau, x = "tau", title = "Tau distribution/all genes (ensembl ids); N = 56156")
## Warning: Using `bins = 30` by default. Pick better value with the argument `bins`.
ggplot2::ggsave(file.path(result_dir_path[[sn]], "tau_distribution.pdf"), plot = tau_dist, device = "pdf", width = 6, height = 4)
## Warning: Removed 23809 rows containing non-finite values (`stat_bin()`).

for(sn in sample_names) {
    clust_lens <- unlist(lapply(seqs_clusters_as_list_ordered[[sn]], length))
    message("Ordered list lengths: ", paste(clust_lens, collapse = " "))
    message("Original list lengths:", 
      paste(unlist(lapply(seqs_clusters_as_list[[sn]], length)), collapse = " "))
  
  
    ##
    clust_lab <- rep("0", length(human_result[[sn]]$rawSeqs))
    clust_names <- sort(as.character(1:length(seqs_clusters_as_list_ordered[[sn]])))
    for(i in seq_along(seqs_clusters_as_list_ordered[[sn]])){
        clust_lab[seqs_clusters_as_list_ordered[[sn]][[i]] ] <- clust_names[i]
    }
  
    peakAnno_df <- as.data.frame(perSample_peakAnno[[sn]])
    
    
    # tissSpecScores <- rep(-1, nrow(peakAnno_df))
    tissSpecScores <- vapply(seq(nrow(peakAnno_df)), FUN = function(x){
        idx <- which(tissueSpecTau$gene_id == peakAnno_df[x, "geneId"])
        if(length(idx) < 1){
            return(-1)
        }
        return(tissueSpecTau$tau[idx])
        }, FUN.VALUE=numeric(1))
    
    message("Tissue specificity values (tau) collected")
    
    ### Get DF ob BED information
    samarth_df[[sn]] <- data.frame(chr = human_bed_info[[sn]]$chr, 
                             start = human_bed_info[[sn]]$start,
                             end = human_bed_info[[sn]]$end,
                             strand = human_bed_info[[sn]]$strand,
                             IQW = human_bed_info[[sn]]$IQW,
                             domTPM = human_bed_info[[sn]]$domTPM,
                             geneId = peakAnno_df$geneId,
                             tissueSpecScore = tissSpecScores,
                             clust_ID = clust_lab #these are character sorted
                             )
    
    df_fname <- file.path(result_dir_path[[sn]], "samarth_info_df.csv")
    write.csv(samarth_df, file = df_fname)
    
    ## All, includes repeated ones
    fname <- file.path(result_dir_path[[sn]], "tau_existing_all_genes_distribution.pdf")
    p1 <- ggpubr::gghistogram(data.frame(tissSpecScores = tissSpecScores), 
                        x = "tissSpecScores",  y = "..count..")
    ggplot2::ggsave(filename = fname, plot = p1, device = pdf, width = 6, height = 4)
    
    ## Just for the 4th cluster which looks like TATA-box
    fname <- file.path(result_dir_path[[sn]], "tau_clust_4_distribution.pdf")
    p1 <- ggpubr::gghistogram(
        data.frame(tissSpecScores_clust4 = tissSpecScores[perSample_archR_clusts[[1]][[4]]]), 
        x = "tissSpecScores_clust4", y = "..count..")
    ggplot2::ggsave(filename = fname, plot = p1, device = pdf, width = 6, height = 4)
    
}
## Ordered list lengths: 89 302 196 194 257 52 1739 391 507 1094 1335 462 598 835 771 701
## Original list lengths:835 1335 1739 771 701 462 302 89 257 598 196 1094 194 52 507 391
## Tissue specificity values (tau) collected
## Warning: Using `bins = 30` by default. Pick better value with the argument `bins`.
## Warning: Removed 125 rows containing non-finite values (`stat_bin()`).
## Warning: Using `bins = 30` by default. Pick better value with the argument `bins`.
## Warning: Removed 16 rows containing non-finite values (`stat_bin()`).
```



# Session Info


```
## R version 4.2.2 Patched (2022-11-10 r83330)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux bookworm/sid
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.21.so
## 
## locale:
##  [1] LC_CTYPE=en_DK.UTF-8       LC_NUMERIC=C               LC_TIME=en_DK.UTF-8       
##  [4] LC_COLLATE=en_DK.UTF-8     LC_MONETARY=en_DK.UTF-8    LC_MESSAGES=en_DK.UTF-8   
##  [7] LC_PAPER=en_DK.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
## [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_DK.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] dendextend_1.16.0                 DT_0.26                           clusterProfiler_4.6.0            
##  [4] ggeasy_0.1.3                      reshape2_1.4.4                    CAGEr_2.4.0                      
##  [7] MultiAssayExperiment_1.24.0       SummarizedExperiment_1.28.0       MatrixGenerics_1.10.0            
## [10] matrixStats_0.63.0                rmarkdown_2.19                    readr_2.1.3                      
## [13] dplyr_1.0.10                      ggpubr_0.5.0                      forcats_0.5.2                    
## [16] ggsankey_0.0.99999                patchwork_1.1.2                   cowplot_1.1.1                    
## [19] ggplot2_3.4.0                     GenomicFeatures_1.50.3            ChIPseeker_1.34.1                
## [22] org.Hs.eg.db_3.16.0               AnnotationDbi_1.60.0              Biobase_2.58.0                   
## [25] BSgenome.Hsapiens.UCSC.hg19_1.4.3 BSgenome_1.66.2                   rtracklayer_1.58.0               
## [28] Biostrings_2.66.0                 XVector_0.38.0                    seqArchRplus_0.99.8              
## [31] GenomicRanges_1.50.1              GenomeInfoDb_1.34.6               IRanges_2.32.0                   
## [34] S4Vectors_0.36.1                  BiocGenerics_0.44.0               seqArchR_1.2.0                   
## 
## loaded via a namespace (and not attached):
##   [1] storr_1.2.5                             KEGGREST_1.38.0                        
##   [3] lattice_0.20-45                         RMariaDB_1.2.2                         
##   [5] vctrs_0.5.1                             utf8_1.2.3                             
##   [7] blob_1.2.3                              withr_2.5.0                            
##   [9] lifecycle_1.0.3                         stringr_1.5.0                          
##  [11] munsell_0.5.0                           ragg_1.2.5                             
##  [13] codetools_0.2-19                        magick_2.7.3                           
##  [15] fastmatch_1.1-3                         formula.tools_1.7.1                    
##  [17] stringi_1.7.12                          grid_4.2.2                             
##  [19] polyclip_1.10-4                         rhdf5filters_1.10.0                    
##  [21] yulab.utils_0.0.6                       cluster_2.1.4                          
##  [23] ggraph_2.1.0                            ape_5.6-2                              
##  [25] pkgconfig_2.0.3                         prettyunits_1.1.1                      
##  [27] data.table_1.14.6                       lubridate_1.8.0                        
##  [29] sparseMatrixStats_1.10.0                httr_1.4.4                             
##  [31] igraph_1.3.5                            treeio_1.22.0                          
##  [33] progress_1.2.2                          graphlayouts_0.8.3                     
##  [35] ggfun_0.0.9                             gson_0.0.9                             
##  [37] htmltools_0.5.4                         viridisLite_0.4.1                      
##  [39] yaml_2.3.6                              jquerylib_0.1.4                        
##  [41] pillar_1.8.1                            later_1.3.0                            
##  [43] glue_1.6.2                              DBI_1.1.3                              
##  [45] BiocParallel_1.32.5                     plyr_1.8.8                             
##  [47] gtable_0.3.1                            GOSemSim_2.24.0                        
##  [49] caTools_1.18.2                          fastmap_1.1.0                          
##  [51] archR_0.1.8                             crosstalk_1.2.0                        
##  [53] broom_1.0.2                             promises_1.2.0.1                       
##  [55] textshaping_0.3.6                       ggforce_0.4.1                          
##  [57] hms_1.1.2                               png_0.1-8                              
##  [59] ggtree_3.6.2                            lazyeval_0.2.2                         
##  [61] crayon_1.5.2                            boot_1.3-28.1                          
##  [63] tidyselect_1.2.0                        xfun_0.35                              
##  [65] purrr_1.0.1                             splines_4.2.2                          
##  [67] operator.tools_1.6.3                    seqPattern_1.30.0                      
##  [69] rappdirs_0.3.3                          bit64_4.0.5                            
##  [71] factoextra_1.0.7                        ggsignif_0.6.4                         
##  [73] VGAM_1.1-7                              permute_0.9-7                          
##  [75] evd_2.3-6.1                             cachem_1.0.6                           
##  [77] DelayedArray_0.24.0                     gdata_2.18.0.1                         
##  [79] vegan_2.6-4                             abind_1.4-5                            
##  [81] mime_0.12                               systemfonts_1.0.4                      
##  [83] rjson_0.2.21                            aplot_0.1.9                            
##  [85] ggrepel_0.9.2                           rstatix_0.7.1                          
##  [87] processx_3.8.0                          tools_4.2.2                            
##  [89] cli_3.6.0                               magrittr_2.0.3                         
##  [91] tibble_3.1.8                            Matrix_1.5-3                           
##  [93] ggplotify_0.1.0                         DelayedMatrixStats_1.20.0              
##  [95] assertthat_0.2.1                        qvalue_2.30.0                          
##  [97] fgsea_1.24.0                            HDF5Array_1.26.0                       
##  [99] DECIPHER_2.26.0                         BiocFileCache_2.6.0                    
## [101] mgcv_1.8-41                             tweenr_2.0.2                           
## [103] zlibbioc_1.44.0                         restfulr_0.0.15                        
## [105] biomaRt_2.54.0                          shadowtext_0.1.2                       
## [107] tzdb_0.3.0                              ps_1.7.2                               
## [109] fansi_1.0.4                             tidygraph_1.2.2                        
## [111] KernSmooth_2.23-20                      backports_1.4.1                        
## [113] som_0.3-5.1                             hopach_2.58.0                          
## [115] farver_2.1.1                            bit_4.0.5                              
## [117] gplots_3.1.3                            Rsamtools_2.14.0                       
## [119] BiocIO_1.8.0                            DOSE_3.24.2                            
## [121] scatterpie_0.1.8                        sass_0.4.4                             
## [123] downloader_0.4                          viridis_0.6.2                          
## [125] rstudioapi_0.14                         Rhdf5lib_1.20.0                        
## [127] nlme_3.1-161                            stringdist_0.9.10                      
## [129] webshot2_0.1.0                          gtools_3.9.4                           
## [131] bslib_0.4.2                             rhdf5_2.42.0                           
## [133] generics_0.1.3                          colorspace_2.1-0                       
## [135] base64enc_0.1-3                         XML_3.99-0.13                          
## [137] dbplyr_2.2.1                            RColorBrewer_1.1-3                     
## [139] GenomeInfoDbData_1.2.9                  chromote_0.1.1                         
## [141] evaluate_0.19                           memoise_2.0.1                          
## [143] knitr_1.41                              TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
## [145] Rcpp_1.0.9                              BiocManager_1.30.19                    
## [147] seqLogo_1.64.0                          vroom_1.6.0                            
## [149] jsonlite_1.8.4                          digest_0.6.31                          
## [151] bitops_1.0-7                            RSQLite_2.2.19                         
## [153] wesanderson_0.3.6                       slickR_0.6.0                           
## [155] compiler_4.2.2                          carData_3.0-5                          
## [157] gridGraphics_0.5-1                      rlang_1.0.6                            
## [159] ggseqlogo_0.1                           remake_0.3.0                           
## [161] PWMEnrich_4.34.0                        htmlwidgets_1.6.1                      
## [163] websocket_1.4.1                         labeling_0.4.2                         
## [165] curl_5.0.0                              parallel_4.2.2                         
## [167] filelock_1.0.2                          scales_1.2.1                           
## [169] plotrix_3.8-2                           enrichplot_1.18.3                      
## [171] HDO.db_0.99.1                           gridExtra_2.3                          
## [173] RCurl_1.98-1.9                          car_3.1-1                              
## [175] tidyr_1.2.1                             GO.db_3.16.0                           
## [177] MASS_7.3-58.1                           ellipsis_0.3.2                         
## [179] tidytree_0.4.2                          xml2_1.3.3                             
## [181] R6_2.5.1                                GenomicAlignments_1.34.0
```
