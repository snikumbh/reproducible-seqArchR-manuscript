## archR on drosophila melanogaster Chen et al 2014
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

message("archR v0.1.8/b4d9627 for Chen2014 with bound 1e-06/07/08")

library(archR)


data_path <- file.path("data", "drosophila-chen2014")
results_path <- file.path("results", "drosophila-chen2014", "with_archR_v0.1.8")


set.seed(1234)

inputFastaFname <- file.path(data_path,
                    'from_dm3_2bit_Dmel_carcass_300_300_corrected.fasta')
## 300bp upstream and downstream of the initiator, excluding the initiator
## Initiator at position 301

tss.seqs_raw <- archR::prepare_data_from_FASTA(inputFastaFname, raw_seq = TRUE)


tss.seqs_raw_45U_45D <- Biostrings::subseq(tss.seqs_raw, start=301-45, end=301+45)


tss.seqs_mat_45U_45D <- archR::get_one_hot_encoded_seqs(tss.seqs_raw_45U_45D,
                                    sinuc_or_dinuc = "dinuc")


###

nSeqs <- length(tss.seqs_raw_45U_45D)
positions <- seq(-45,45)


boundValsToRun <- 10^c(-1*c(8))
innerChunkSizesToRun <- c(2000)


for(boundVal in boundValsToRun){
    for(innerChunkSizeVal in innerChunkSizesToRun){
        archRconfig <- archR::archR_set_config(
            parallelize = FALSE,
            n_cores = ncores,
            n_runs = 100,
            k_min = 1,
            k_max = 20,
            mod_sel_type =  "stability",
            bound = boundVal,
            min_size = 50,
            result_aggl = "ward.D",
            result_dist = "euclid",
            chunk_size = innerChunkSizeVal,
            flags = list(debug = TRUE, time = TRUE, verbose = TRUE,
                plot = TRUE)
        )
        perform_iters <- 5
        collationStrategy <- c(FALSE, TRUE, FALSE, FALSE, FALSE)
        collate_str <- paste(unlist(lapply(collationStrategy, function(x){
            ifelse(x, 'T', 'F')
        })), collapse= "")

        path_prefix <- "archR_result_drosophila-chen2014_serial_run"
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
                                    seqs_ohe_mat = tss.seqs_mat_45U_45D,
                                    seqs_raw = tss.seqs_raw_45U_45D,
                                    total_itr = perform_iters,
                                    seqs_pos = positions,
                                    set_ocollation = collationStrategy,
                                    fresh = TRUE,
                                    o_dir = result_dir_path)
        message("Complete stability run with bound ", boundVal)
    }
}

message("All completed!")
