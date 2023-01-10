## archR on Homo Sapiens
# Experiments run on server


args <- commandArgs(trailingOnly = TRUE)

conda_env_name <- args[1]
conda_path <- args[2]
ncores <- as.numeric(args[3])

print(ncores)
stopifnot(is.numeric(ncores) && ncores > 0)

## Ensure version/tag archR@b4d9627
repo_name <- "snikumbh/archR@b4d9627"
message("Installing archR: ", repo_name)
remotes::install_github(repo = repo_name, force = TRUE, upgrade = FALSE,
                        quiet = FALSE)


library(archR)
# Setup python with reticulate
reticulate::use_condaenv(condaenv = conda_env_name,
    conda = conda_path,
    required = TRUE)
# Check if OK?
reticulate::import("sklearn")



use_minTPM <- 1
data_path <- file.path("data", "human")
results_path <- file.path("results", "human", "newly_processed", paste0("minTPM_", use_minTPM))
inr_analysis <- TRUE

set.seed(1234)

sample_names <- c("HepG2_cell_rep12", "hMSC-BM_cell_rep12",
                  "SkMC_cell_rep1", "SkMC_cell_rep2",
                  "human_cellGroup_merged")


process_single_sample <- function(sample_name){

    inputFastaFname <- file.path(data_path,
                                 paste0('samarth_hsapiens_TC_sample_',
                                        sample_name, '_minTPM', use_minTPM,
                                        '_flank_up500_flank_down500_all70.fa'))
    message("Using fasta file:", inputFastaFname)

    ## 500bp upstream and downstream of the initiator, excluding the initiator
    ## Initiator at position 501

    tss.seqs_raw <- archR::prepare_data_from_FASTA(fasta_fname = inputFastaFname,
                                                   raw_seq = TRUE)



    use_subseq <- TRUE
    left_flank <- 50
    right_flank <- 150
    nSeqs <- length(tss.seqs_raw)
    tssPosRel500 <- 501

    ## for INR-removed analysis
    if(inr_analysis){
        leftFlankKO <- 5
        rightFlankKO <- 5
        tssPos <- 51
        # seqLen <- Biostrings::width(tss.seqs_raw[1])
        ## This will work if seqLen is odd, to find the right inr position


        endLeft <- tssPosRel500-leftFlankKO-1
        tss.seqs_raw_inr_left <- Biostrings::subseq(tss.seqs_raw,
                                                    start=tssPosRel500 - left_flank,
                                                    end=endLeft)
        ##
        startRight <- tssPosRel500+rightFlankKO+1
        tss.seqs_raw_inr_right <- Biostrings::subseq(tss.seqs_raw,
                                                     start= startRight,
                                                     end=tssPosRel500+right_flank)

        tss.seqs_raw_inr_excluded <- Biostrings::DNAStringSet(paste0(tss.seqs_raw_inr_left,
                                                              tss.seqs_raw_inr_right))
        tss.seqs_raw <- tss.seqs_raw_inr_excluded
        message("INR removed w/ flanks ", -1*leftFlankKO, " and ", rightFlankKO)
        positions <- seq(-left_flank, right_flank)[c(1:(tssPos-leftFlankKO+1),
                                                    (tssPos+rightFlankKO+1):(tssPos+right_flank))]
        print(positions)
        ## This version of archR throws an error if 0 is not included in posiitons
        message("positions altered")
        positions <- seq(-left_flank, right_flank)[c(1:(tssPos-leftFlankKO+1), tssPos,
                                                    (1+(tssPos+rightFlankKO+1)):(tssPos+right_flank))]
        print(positions)
        #c(seq(-50, -5), seq(5,150))
        use_inr_suffix <- paste0("inrKO_", leftFlankKO, "_", rightFlankKO)
    }else{
        use_inr_suffix <- ""
        message("Full analysis/no TSS-flanks removed")
        if(use_subseq){
            tss.seqs_raw <- Biostrings::subseq(tss.seqs_raw,
                                               start = tssPosRel500 - left_flank,
                                               end = tssPosRel500 + right_flank)
        }
        positions <- seq(-left_flank, right_flank)


    }


    tss.seqs_mat <- archR::get_one_hot_encoded_seqs(tss.seqs_raw,
                                                    sinuc_or_dinuc = "dinuc")



    boundValsToRun <- 10^c(-1*c(4,6,8))
    chunkSizesToRun <- c(5000)

    for(chunkSizeVal in chunkSizesToRun){
        for(boundVal in boundValsToRun){
            archRconfig <- archR::archR_set_config(
                parallelize = TRUE,
                n_cores = ncores,
                n_runs = 100,
                k_min = 1,
                k_max = 20,
                mod_sel_type = "stability",
                bound = boundVal,
                min_size = 50,
                result_aggl = "ward.D",
                result_dist = "cor",
                chunk_size = chunkSizeVal,
                flags = list(debug = TRUE, time = TRUE, verbose = TRUE,
                             plot = TRUE)
            )
            perform_iters <- 5
            collationStrategy <- c(FALSE, TRUE, TRUE, TRUE, FALSE)
            collate_str <- paste(unlist(lapply(collationStrategy, function(x){
                ifelse(x, 'T', 'F')
            })), collapse= "")

            path_prefix <- paste0("archR_result_hsapiens_", sample_name,
                                  "_flank_up", left_flank,"_flank_down", right_flank, "_")
            #
            path_suffix_by_config <- paste(
                use_inr_suffix,
                "modSelType", archRconfig$modSelType,
                "chunkSize", archRconfig$chunkSize,
                "bound", format(archRconfig$bound, scientific = TRUE),
                "collate", collate_str,
                sep = "_")
            #
            result_dir_path <- file.path(results_path, paste(path_prefix,
                                                             path_suffix_by_config, sep = "_"))
            #

            message("Processing ", nSeqs, " sequences")
            archRresult <- archR::archR(config = archRconfig,
                                        seqs_ohe_mat = tss.seqs_mat,
                                        seqs_raw = tss.seqs_raw,
                                        total_itr = perform_iters,
                                        seqs_pos = positions,
                                        set_ocollation = collationStrategy,
                                        fresh = TRUE,
                                        o_dir = result_dir_path)
            message("Complete stability run with bound ", boundVal, "for", sample_name)
        }
    }
    message("Sample", sample_name, " completed")
    ##
} ## func process_single_sample ENDS


for(sample_name in sample_names[5]){
    message("archR v0.1.8/b4d9627 for Hsapiens/", sample_name)
    process_single_sample(sample_name)
}


message("All completed!")
stop("SAMARTH")
######
## Function to invisibly capture plot in a temp file and return only computed
## info from logomaker object
invisible_logomaker <- function(...){
    tf <- tempfile()
    png(filename=tf)
    edlogo_info <- Logolas::logomaker(...)
    dev.off()
    unlink(tf)
    edlogo_info
}

# results_path <- file.path("experiments", results_path)

use_leftflank <- 50
use_rightflank <- 5#150
use_prefix <- file.path(results_path, paste0("archR_result_hsapiens_human_cellGroup_merged_flank_up", use_leftflank, "_flank_down", use_rightflank, "__modSelType_stability_chunkSize_"))
use_suffix <- "_collate_FTTTF"
use_chunkSz <- c(5000, 10000)
use_chunkSz <- as.character(use_chunkSz)
use_bound <- c(10^(-1*c(4,6,8)))
use_bound <- as.character(use_bound)
path_comb <- expand.grid(use_prefix, use_chunkSz, "_bound_", use_bound, use_suffix)
final_paths <- apply(path_comb, 1, paste0, collapse="")
print_edlogo <- FALSE
print_ymax_logo <- TRUE


for(final_path in final_paths){
    # final_path <- file.path(results_path, "archR_result_hsapiens_human_cellGroup_merged_modSelType_stability_chunkSize_1000_bound_1e-08_collate_FTFFF")
    result <- readRDS(file.path(final_path,"archRresult.rds"))

    use_labels <- c(-use_leftflank:use_rightflank) #c(-50:150)#c(-50:-6, 6:150)

    for(itr in c(1:5)){
        seq_clusters <- archR::get_seqs_clust_list(result$seqsClustLabels[[itr]])
        ##
        edlogo_info <- vector("list", length(seq_clusters))
        edlogo_mat <- vector("list", length(seq_clusters))
        elogo_mat <- vector("list", length(seq_clusters))
        ## using for loop

        if(print_ymax_logo){
            archR::plot_arch_for_clusters(seqs = result$rawSeqs, clust_list = seq_clusters,
                          bits_yax = "auto", pos_lab = use_labels, set_titles = TRUE,
                          pdf_name = file.path(final_path, paste0("Architectures_ymax_itr", itr, ".pdf")),
                          pdf_width = 15, pdf_height = 2)
        }

        if(print_edlogo){
            for(i in seq_along(seq_clusters)){
                seqs <- seq_clusters[[i]]
                sam <- Biostrings::alphabetFrequency(result$rawSeqs[seqs], as.prob = TRUE)
                ## Remove any sequences that have N
                ## Report how many removed
                if(any(sam[,"N"] > 0)){
                    idx <- which(sam[,"N"] > 0)
                    seqs <- seqs[-c(idx)]
                    sam <- Biostrings::alphabetFrequency(result$rawSeqs[seqs], as.prob = TRUE)
                }
                sam <- sam[, 1:4]
                use_bg <- colSums(sam)/length(seqs)
                #prod(dim(sam))
                ##
                edlogo_info[[i]] <- invisible_logomaker(
                    as.character(result$rawSeqs[seqs]),
                    type="EDLogo",
                    bg = use_bg,
                    colors = c("blue", "red", "darkgreen", "orange"),
                    return_heights = TRUE)
                ##
                edlogo_mat[[i]] <-
                    t(t(edlogo_info[[i]]$table_mat_pos_norm) * edlogo_info[[i]]$pos_ic) +
                    (-1*t(t(edlogo_info[[i]]$table_mat_neg_norm) * edlogo_info[[i]]$neg_ic))

                elogo_mat[[i]] <-
                    t(t(edlogo_info[[i]]$table_mat_pos_norm) * edlogo_info[[i]]$pos_ic)
                # +
                #     (-1*t(t(edlogo_info[[i]]$table_mat_neg_norm) * edlogo_info[[i]]$neg_ic))

            }
            ##
            cumsums <- cumsum(unlist(lapply(seq_clusters, length)))
            # print(cumsums)
            starts <- c(1,cumsums[1:(length(cumsums)-1)])
            ends <- cumsums
            ## print elogo
            pdf(file.path(final_path, paste0("Elogo_Itr", itr, ".pdf")), width = 11, height = 2)
            for(i in seq_along(seq_clusters)){
                sam_elogo <- archR::plot_ggseqlogo(pwm_mat = elogo_mat[[i]],
                                                    method = "custom",
                                                    pos_lab = use_labels)

                sam_elogo <- sam_elogo + ggplot2::scale_x_continuous(
                    # breaks = c(seq(1,45, by=5), 45, seq(46,190,by=5)),
                    # labels = use_labels[c(seq(1,45, by=5), 45, seq(46,190,by=5))],
                    breaks = seq(1, 201, by = 5),
                    labels = seq(-50, 150, by = 5),
                    expand = ggplot2::expansion(mult = c(0, 0))) +
                    ggplot2::ggtitle(label = paste0("Arch ", i, ", ",
                                                    length(seq_clusters[[i]]), " sequences ",
                                                    "(", starts[i], "-", ends[i], ")"))
                print(sam_elogo)
            }

            dev.off()

            pdf(file.path(final_path, paste0("EDlogo_Itr", itr, ".pdf")), width = 11, height = 2)
            for(i in seq_along(seq_clusters)){
                sam_edlogo <- archR::plot_ggseqlogo(pwm_mat = edlogo_mat[[i]],
                                                    method = "custom",
                                                    pos_lab = use_labels)

                sam_edlogo <- sam_edlogo + ggplot2::scale_x_continuous(
                    # breaks = c(seq(1,45, by=5), 45, seq(46,190,by=5)),
                    # labels = use_labels[c(seq(1,45, by=5), 45, seq(46,190,by=5))],
                    breaks = seq(1, 201, by = 5),
                    labels = seq(-50, 150, by = 5),
                    expand = ggplot2::expansion(mult = c(0, 0))) +
                    ggplot2::ggtitle(label = paste0("Arch ", i, ", ",
                                                    length(seq_clusters[[i]]), " sequences ",
                                                    "(", starts[i], "-", ends[i], ")"))
                print(sam_edlogo)
            }

            dev.off()
            }
    }
}
######
# ## using lapply, Logolas logomaker doesn't work
# pdf(file.path(results_path, "EDlogo.pdf"), width = 11, height = 3)
# seq_logos <- lapply(seq_along(seq_clusters), function(x){
#     seqs <- seq_clusters[[x]]
#     sam <- Biostrings::alphabetFrequency(tss.seqs_raw[seqs])
#     sam <- sam[, 1:4]
#     use_bg <- colSums(sam)/(nrow(sam)*Biostrings::width(tss.seqs_raw[1]))
#     Logolas::logomaker(
#         as.character(tss.seqs_raw[seqs]),
#         type="EDLogo",
#         bg = use_bg,
#         colors = c("blue", "red", "darkgreen", "orange"),
#         return_heights = TRUE)
#      })
# print(seq_logos)
# dev.off()

# sam <- Biostrings::alphabetFrequency(x = tss.seqs_raw[seq_clusters[2][[1]]])
# sam <- sam[,1:4]
# (use_bg <- colSums(sam)/nrow(sam))


# sam_edlogo <- archR::plot_ggseqlogo(pwm_mat = sam_icm,
#                       method = "custom",
#                       pos_lab = c(-50:150)
#
# sam_edlogo <- sam_edlogo + ggplot2::scale_x_continuous(
#     breaks = seq(1, 201, by = 5),
#     labels = seq(-50, 150, by = 5),
#     expand = ggplot2::expansion(mult = c(0, 0)))





