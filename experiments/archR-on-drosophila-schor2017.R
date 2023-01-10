## archR on drosophila melanogaster Schor et al 2017
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

message("archR v0.1.8/b4d9627 for Schor2017")
library(archR)


data_path <- file.path("data", "drosophila-schor2017")
results_path <- file.path("results", "drosophila-schor2017", "with_archR_v0.1.8")


set.seed(1234)

# sample_names <- paste(rep(paste0('RAL', c(28)), each=3),
#                     rep(paste(c(2,6,10), 'to', c(4,8,12), sep="_"), times=1),
#                     sep="_")
# sample_names[3] <- "RAL28_10_to_12_sample12"

sample_names <- paste(rep(paste0('RAL', c(805, 642, 639)), each=3),
                      rep(paste(c(2,6,10), 'to', c(4,8,12), sep="_"), times=3),
                      sep="_")
#sample_names[3] <- "RAL28_10_to_12_sample12"



process_single_sample <- function(sample_name){

    inputFastaFname <- file.path(data_path,
        paste0('dm6_samarth_schor_et_al_TC_sample_', sample_name,
               '_minTPM1_flank_up500_flank_down500.fa'))

    ## 500bp upstream and downstream of the initiator, excluding the initiator
    ## Initiator at position 501

    tss.seqs_raw <- archR::prepare_data_from_FASTA(fasta_fname = inputFastaFname,
                                            raw_seq = TRUE)


    tss.seqs_raw_45U_45D <- Biostrings::subseq(tss.seqs_raw, start=501-45, end=501+45)


    tss.seqs_mat_45U_45D <- archR::get_one_hot_encoded_seqs(tss.seqs_raw_45U_45D,
        sinuc_or_dinuc = "dinuc")

    nSeqs <- length(tss.seqs_raw_45U_45D)
    positions <- seq(-45,45)

    boundValsToRun <- 10^c(-1*c(6:8))
    innerChunkSizesToRun <- c(5000)

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
                min_size = 25,
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
            print(str(archRconfig))
            path_prefix <- paste0("archR_result_drosophila_schor2017_",
                                sample_name)
            #
            path_suffix_by_config <- paste(
                "modSelType", archRconfig$modSelType,
                "chunkSize", archRconfig$chunkSize,
                "bound", format(archRconfig$bound, scientific = TRUE),
                "aggl", archRconfig$result_aggl,
                "dist", archRconfig$result_dist,
                "collate", collate_str,
                sep = "_")
            #
            result_dir_path <- file.path(results_path, paste(path_prefix,
                path_suffix_by_config, sep = "_"))
            #

            message("Processing ", nSeqs, " sequences")
            archRresult <- archR::archR(config = archRconfig,
                seqs_ohe_mat = tss.seqs_mat_45U_45D,
                seqs_raw = tss.seqs_raw_45U_45D,
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
