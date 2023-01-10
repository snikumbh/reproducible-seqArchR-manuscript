## archR on Zebrafish
# Experiments run on server


args <- commandArgs(trailingOnly = TRUE)

conda_env_name <- args[1]
conda_path <- args[2]
ncores <- as.numeric(args[3])

stopifnot(is.numeric(ncores) && ncores > 0)

## Ensure version/tag archR@b4d9627
repo_name <- "snikumbh/archR@b4d9627"
message("Installing archR: ", repo_name)
remotes::install_github(repo = repo_name, force = TRUE, upgrade = FALSE,
                        quiet = FALSE)

message("archR v0.1.8/b4d9627 for Zebrafish NepalEtAl2013 with bound 1e-06/07/08")
library(archR)
# Setup python with reticulate
reticulate::use_condaenv(condaenv = conda_env_name,
    conda = conda_path,
    required = TRUE)
# Check if OK?
reticulate::import("sklearn")

data_path <- file.path("data", "zebrafish-nepal2013/new")
results_path <- file.path("results", "zebrafish-nepal2013/with_archR_v0.1.8")


set.seed(1234)

sample_names <- c("64_cells","dome_30perc_epiboly", "prim6")


process_single_sample <- function(sample_name){

    inputFastaFname <- file.path("data", paste0('samarth_zf_TC_sample_Drerio_', sample_name, '_flank_500.fa'))
    ## 500bp upstream and downstream of the initiator, excluding the initiator
    ## Initiator at position 501

    tss.seqs_raw <- archR::prepare_data_from_FASTA(fasta_fname = inputFastaFname,
        raw_seq = TRUE)

    left_flank <- 45
    right_flank <- 150
    tss.seqs_raw <- Biostrings::subseq(tss.seqs_raw, start=501-left_flank,
                                    end=501+right_flank)


    tss.seqs_mat <- archR::get_one_hot_encoded_seqs(tss.seqs_raw,
        sinuc_or_dinuc = "dinuc")

    nSeqs <- length(tss.seqs_raw)
    positions <- seq(-left_flank, right_flank)

    boundValsToRun <- 10^c(-1*c(8))
    innerChunkSizesToRun <- c(5000, 10000)

    for(innerChunkSizeVal in innerChunkSizesToRun){
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
                result_dist = "euclid",
                chunk_size = innerChunkSizeVal,
                flags = list(debug = TRUE, time = TRUE, verbose = TRUE,
                    plot = TRUE)
            )
            perform_iters <- 5
            collationStrategy <- c(FALSE, TRUE, TRUE, TRUE, FALSE)
            collate_str <- paste(unlist(lapply(collationStrategy, function(x){
                ifelse(x, 'T', 'F')
            })), collapse= "")

            path_prefix <- paste0("archR_result_zebrafish_nepal2013_",
                sample_name)
            #
            path_suffix_by_config <- paste(
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


for(sample_name in sample_names){
    process_single_sample(sample_name)
}


message("All completed!")



