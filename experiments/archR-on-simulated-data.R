# archR on synthetic data
# Experiments run on server


args <- commandArgs(trailingOnly = TRUE)

conda_env_name <- args[1]
conda_path <- args[2]
ncores <- as.numeric(args[3])
nruns <- as.numeric(args[4])


stopifnot(is.numeric(ncores) && ncores > 0)


## Ensure version/tag archR@b4d9627
repo_name <- "snikumbh/archR@b4d9627"
message("Installing archR: ", repo_name)
remotes::install_github(repo = repo_name, force = TRUE, upgrade = FALSE,
                        quiet = FALSE)

message("archR v0.1.8/b4d9627 for simdata with bound 1e-08")


# Load libraries
library(archR, lib.loc = "/mnt/biggley/home/sarvesh/R/x86_64-pc-linux-gnu-library/4.1")
# Setup python with reticulate
reticulate::use_condaenv(condaenv = conda_env_name,
             conda = conda_path,
             required = TRUE)
# Check if OK?
reticulate::import("sklearn")


##--------------------------- SETUP TEST CONFIG --------------------------------
# -- Test with n different random seeds for stability check of clustering
nRuns <- nruns

# -- Test with chunkSize 200, 500 and 1000 for nSeqs = 1000
#    For nSeqs = 5000, test chunkSizes 500, 1000, 2000, 5000
#    For nSeqs = 10000, test chunkSizes 500, 1000, 2000, 5000, 10000
chunkSizes_list <- list("n1000" = c(200, 500, 1000),
                   "n5000" = c(500, 1000, 2000, 5000),
                   "n10000" = c(500, 1000, 2000, 5000, 10000),
                   "n25000" = c(2000, 5000, 10000, 15000, 25000),
                   "n50000" = c(5000, 10000, 15000, 25000, 50000))
# use a different collation strategy for chunk_size 1000, because total
# number fo sequences is 1000
collationStrategy <- list(c(TRUE, FALSE))

# -- Test with bound values 10^-8
boundVals <- c(10^-8)

# -- Test speedups so parallelize and serially run experiments
parallelizeVals <- c(TRUE)
##--------------------------- SETUP TEST CONFIG ENDS ---------------------------

##--------------------------- SETUP FILES --------------------------------------
data_path <- file.path("data", "simulated-data")
results_path <- file.path("results", "simulated-data-archR-v0.1.8-bound-1e08", paste0(ncores, "cores"))
if(!dir.exists(results_path)){
    message("Creating dir ", results_path)
    dir.create(results_path)
}else{
    message("Dir exists, not creating a new one ", results_path)
}

# Load fasta files
# -- For n1000, all combinations of mutation rates and mutation positions
# -- For n5000, only mu = 0.1 and p = 1
# -- For n10000, only mu = 0.1 and p = 1
file_names_list <- list(
                      "n1000" = list.files(data_path, pattern = "*n1000.fa$",
                      recursive = FALSE, full.names = TRUE)#,
                      # "n5000" = list.files(data_path, pattern = "*mu0.1_p1_n5000.fa$",
                      #   recursive = FALSE, full.names = TRUE),
                      # "n10000" = list.files(data_path, pattern = "*mu0.1_p1_n10000.fa$",
                      #   recursive = FALSE, full.names = TRUE),
                      # "n25000" = list.files(data_path, pattern = "*mu0.1_p1_n25000.fa$",
                      #   recursive = FALSE, full.names = TRUE),
                      # "n50000" = list.files(data_path, pattern = "*mu0.1_p1_n50000.fa$",
                      #                 recursive = FALSE, full.names = TRUE)
                      )

message("Files to process")
print(file_names_list)


##--------------------------- SETUP FILES ENDS ---------------------------------


##----------------------------RUN EXPS IN A FUNCTION ---------------------------
run_archR_experiment <- function(fname, dir_name, seed_val, results_df,
                                 useNCores = 2){

    set.seed(seed_val)
    message("*****Now processing*****: ", fname)
    param_str <- tail(unlist(strsplit(strsplit(fname, split=".fa")[[1]], "-")),1)
    #
    inputSeqsMat <- archR::prepare_data_from_FASTA(fname,
                                                   sinuc_or_dinuc = "dinuc")
    inputSeqsRaw <- archR::prepare_data_from_FASTA(fname, raw_seq = TRUE)

    nSeqs_this <- length(inputSeqsRaw)
    positions <- seq(1, Biostrings::width(inputSeqsRaw[1]))


    archR::viz_seqs_acgt_mat_from_seqs(seqs = as.character(inputSeqsRaw),
                        pos_lab = positions,
                        save_fname =
                        file.path(dir_name, paste0("inputSeqsRaw_GroundTruth_",
                            param_str, ".png"))
                        )

    # Randomize the sequence order
    changedOrder <- sample.int(nSeqs_this, nSeqs_this, replace = FALSE)
    inputSeqsMat <- inputSeqsMat[ , changedOrder]
    inputSeqsRaw <- inputSeqsRaw[changedOrder]


    chunkSizes <- switch(as.character(nSeqs_this),
           "1000"= chunkSizes_list$n1000,
           "5000"= chunkSizes_list$n5000,
           "10000"= chunkSizes_list$n10000,
           "25000" = chunkSizes_list$n25000,
           "50000" = chunkSizes_list$n50000
           )

    for(parallelize in parallelizeVals){
    for(chunk_size in chunkSizes){
    for(bound_val in boundVals){
            # Set archR configuration
            archRconfig <- archR::archR_set_config(
                parallelize = parallelize,
                n_cores = useNCores,
                n_runs = 100,
                k_min = 1,
                k_max = 20,
                mod_sel_type = "stability",
                bound = bound_val,
                chunk_size = chunk_size,
                result_dist = "euclid",
                result_aggl = "ward.D",
                flags = list(debug = TRUE,
                             time = TRUE,
                             verbose = TRUE,
                             plot = TRUE)
            )
            #
            use_collation <- collationStrategy[[1]]
            #
            path_prefix <- file.path(dir_name,
                                     paste("archR_result_simulateddata",
                                           param_str, sep = "_"))
            path_suffix_by_config <- paste(
                "modSelType", archRconfig$modSelType,
                "chunkSize", archRconfig$chunkSize,
                "bound", format(archRconfig$bound, scientific = TRUE),
                sep = "_")
            #
            result_dir_path <- paste(path_prefix,
                                     path_suffix_by_config, sep = "_")
            #
            archR::viz_seqs_acgt_mat_from_seqs(seqs = as.character(inputSeqsRaw),
                        pos_lab = positions,
                        save_fname =
                        file.path(dir_name, paste0("inputSeqsRaw_Randomized_",
                              param_str, ".png"))
                        )
            # Call/Run archR
            perform_iters <- 2
            if(nSeqs_this == 1000){
                # cat(paste0("Calling archR when nSeqs = ", nSeqs_this, "\n"))
                archRresult <- archR::archR(config = archRconfig,
                                seqs_ohe_mat = inputSeqsMat,
                                seqs_raw = inputSeqsRaw,
                                seqs_pos = positions,
                                total_itr= perform_iters,
                                set_ocollation = use_collation,
                                o_dir = result_dir_path
                                )
            }
            archRresult <- readRDS(file.path(result_dir_path, "archRresult.rds"))
            ##
            check_mutrate <- get_mut_rate_from_param_str(param_str)
            check_mutpos <- get_mut_pos_from_param_str(param_str)
            ##
            relevant_idx <- dplyr::filter(results_df,
                                          nSeqs == nSeqs_this &
                                          parallel == parallelize &
                                          chunkSize == chunk_size &
                                          bound == bound_val &
                                          mutrate == check_mutrate &
                                          mutpos == check_mutpos &
                                          seed == seed_val
                                          )$index
            print(result_dir_path)
            message("Number of rows to fill: ", length(relevant_idx))
            message("Filling row: ", relevant_idx)
            stopifnot(length(relevant_idx) == 1)
            ##
            # ## With v0.1.8, if the overall time taken to complete an iteration
            # ## is less than a minute, it is specified in seconds instead of minutes.
            # ## This creates confusion if the table has some time info in minutes and
            # ## some in seconds. Hence, we use the same strategy as used for NPLB
            # ## to calculate time duration based on files.
            # ##
            # results_df[relevant_idx,
            #            "runTime_m_Iter1"] <- archRresult$timeInfo[[1]]
            # results_df[relevant_idx,
            #            "runTime_m_Iter2"] <- archRresult$timeInfo[[2]]
            ##
            startTime_iter1 <- file.mtime(file.path(result_dir_path, "allSequencesLogo.pdf"))
            startTime_iter2 <- file.mtime(file.path(result_dir_path, "archRresult_checkpoint1.rds"))
            endTime_iter1 <- file.mtime(file.path(result_dir_path, "ClusteringImage_Iteration1.png"))
            endTime_iter2 <- file.mtime(file.path(result_dir_path, "ClusteringImage_Iteration2.png"))

            dur_iter1 <- as.numeric(difftime(endTime_iter1, startTime_iter1, units ="sec"))
            dur_iter2 <- as.numeric(difftime(endTime_iter2, startTime_iter2, units ="sec"))
            results_df[relevant_idx, "runTime_m_Iter1"] <- dur_iter1
            results_df[relevant_idx, "runTime_m_Iter2"] <- dur_iter2
            ##
            trueLabels <- unlist(lapply(names(archRresult$rawSeqs), function(x){
                                        unlist(strsplit(x, "_clust"))[2]
                                    }))
            ##
            ## Get cluster labels at 2 iterations + final clustSol
            predLabels <- vector("list", 3)
            ## Iteration 1 labels
            predLabels[[1]] <- archRresult$seqsClustLabels[[1]]
            ## Iteration 2 labels
            predLabels[[2]] <- archRresult$seqsClustLabels[[2]]
            ##
            predLabels[[3]] <- archRresult$clustSol$seqsClustLabels
            ##
            results_df[relevant_idx, "ARI_Iter1"] <-
                mclust::adjustedRandIndex(predLabels[[1]], trueLabels)
            results_df[relevant_idx, "ARI_Iter2"] <-
                mclust::adjustedRandIndex(predLabels[[2]], trueLabels)
            results_df[relevant_idx, "ARI"] <-
                mclust::adjustedRandIndex(predLabels[[3]], trueLabels)
            ##
            ##
            results_df[relevant_idx, "nClusters_Iter1"] <-
                archRresult$clustBasisVectors[[1]]$nBasisVectors
            results_df[relevant_idx, "nClusters_Iter2"] <-
                archRresult$clustBasisVectors[[2]]$nBasisVectors
            results_df[relevant_idx, "nClusters"] <-
                length(archRresult$clustSol$clusters)
            ##
        } ## parallelize for loop ends
        } ## bound_values for loop ends
    } ## chunkSizes for loop ends
    ##
    ## write updated results_df
    # write.table(results_df, file.path(results_path, "archR_result_simulated_data_summary.tsv"),
    #             sep = "\t", col.names = TRUE, row.names = FALSE)
    ##
    return(results_df)
}
##----------------------------RUN EXPS IN A FUNCTION ENDS-----------------------

get_nSeqs_from_param_str <- function(param_str){
    as.integer(strsplit(unlist(strsplit(param_str, split = "_"))[3], split = "n")[[1]][2])
}

get_mut_rate_from_param_str <- function(param_str){
    as.double(strsplit(unlist(strsplit(param_str, split = "_"))[1], split = "mu")[[1]][2])
}

get_mut_pos_from_param_str <- function(param_str){
    as.integer(strsplit(unlist(strsplit(param_str, split = "_"))[2], split = "p")[[1]][2])
}


set.seed(11992288)
seed_values <- sample.int(.Machine$integer.max, size = nRuns, replace = FALSE)

# -- Create df to store result information
params_sett_list <- lapply(unlist(file_names_list), function(x){
    tail(unlist(strsplit(strsplit(x, split=".fa")[[1]], "-")),1)
})

mutation_rate <-
    as.double(unique(unlist(lapply(params_sett_list, function(x){
    get_mut_rate_from_param_str(x)
}))))
mutation_positions <-
    as.integer(unique(unlist(lapply(params_sett_list, function(x){
    get_mut_pos_from_param_str(x)
}))))
#
results_df_1000 <- purrr::cross_df(list(mutrate = mutation_rate,
                                   mutpos = mutation_positions,
                                   seed = seed_values,
                                   parallel = parallelizeVals,
                                   nSeqs = 1000,
                                   chunkSize = chunkSizes_list$n1000,
                                   bound = boundVals,
                                   runTime_m_Iter1 = -1.0,
                                   runTime_m_Iter2 = -1.0,
                                   ARI_Iter1 = -1.0,
                                   ARI_Iter2 = -1.0,
                                   ARI = -1.0,
                                   nClusters_Iter1 = -1,
                                   nClusters_Iter2 = -1,
                                   nClusters = -1
))
## ^ creates sufficient rows for n1000

## for n5000 and n10000, only one combination of mutrate and mutpos
results_df_5000 <- purrr::cross_df(list(mutrate = c(0.1),
                                   mutpos = c(1),
                                   seed = seed_values,
                                   parallel = parallelizeVals,
                                   nSeqs = 5000,
                                   chunkSize = chunkSizes_list$n5000,
                                   bound = boundVals,
                                   runTime_m_Iter1 = -1.0,
                                   runTime_m_Iter2 = -1.0,
                                   ARI_Iter1 = -1.0,
                                   ARI_Iter2 = -1.0,
                                   ARI = -1.0,
                                   nClusters_Iter1 = -1,
                                   nClusters_Iter2 = -1,
                                   nClusters = -1
))


results_df_10000 <- purrr::cross_df(list(mutrate = c(0.1),
                                        mutpos = c(1),
                                        seed = seed_values,
                                        parallel = parallelizeVals,
                                        nSeqs = 10000,
                                        chunkSize = chunkSizes_list$n10000,
                                        bound = boundVals,
                                        runTime_m_Iter1 = -1.0,
                                        runTime_m_Iter2 = -1.0,
                                        ARI_Iter1 = -1.0,
                                        ARI_Iter2 = -1.0,
                                        ARI = -1.0,
                                        nClusters_Iter1 = -1,
                                        nClusters_Iter2 = -1,
                                        nClusters = -1
))

results_df_25000 <- purrr::cross_df(list(mutrate = c(0.1),
                                        mutpos = c(1),
                                        seed = seed_values,
                                        parallel = parallelizeVals,
                                        nSeqs = 25000,
                                        chunkSize = chunkSizes_list$n25000,
                                        bound = boundVals,
                                        runTime_m_Iter1 = -1.0,
                                        runTime_m_Iter2 = -1.0,
                                        ARI_Iter1 = -1.0,
                                        ARI_Iter2 = -1.0,
                                        ARI = -1.0,
                                        nClusters_Iter1 = -1,
                                        nClusters_Iter2 = -1,
                                        nClusters = -1
))


results_df_50000 <- purrr::cross_df(list(mutrate = c(0.1),
                                         mutpos = c(1),
                                         seed = seed_values,
                                         parallel = parallelizeVals,
                                         nSeqs = 50000,
                                         chunkSize = chunkSizes_list$n50000,
                                         bound = boundVals,
                                         runTime_m_Iter1 = -1.0,
                                         runTime_m_Iter2 = -1.0,
                                         ARI_Iter1 = -1.0,
                                         ARI_Iter2 = -1.0,
                                         ARI = -1.0,
                                         nClusters_Iter1 = -1,
                                         nClusters_Iter2 = -1,
                                         nClusters = -1
))


results_df <- tibble::add_row(results_df_1000, results_df_5000)
results_df <- tibble::add_row(results_df, results_df_10000)
results_df <- tibble::add_row(results_df, results_df_25000)
results_df <- tibble::add_row(results_df, results_df_50000)

results_df <- tibble::add_column(results_df, .before = "mutrate",
                                 index = seq_len(nrow(results_df)))

file_names_list <- as.vector(unlist(file_names_list))


for(seed_val in seed_values){
    #
    dir_name <- file.path(results_path, paste0("seed_", seed_val))
    if(!dir.exists(dir_name)){
        message("Creating directory: ", dir_name)
        dir.create(dir_name)
    }
    #
    for(file_name in file_names_list){
        results_df <- run_archR_experiment(file_name, dir_name, seed_val,
                             results_df, useNCores = ncores)
    }
}

# write.table(results_df, file.path(results_path, "archR_result_simulated_data_summary.tsv"),
#             sep = "\t", col.names = TRUE, row.names = FALSE)





