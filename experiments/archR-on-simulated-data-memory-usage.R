## Benchmark memory usage for seqArchR
##

# archR on synthetic data
# Experiments run on server


args <- commandArgs(trailingOnly = TRUE)

conda_env_name <- args[1]
conda_path <- args[2]
ncores <- as.numeric(args[3])
chooserun <- as.numeric(args[4])


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
# nRuns <- nruns

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


data_path <- file.path("data", "simulated-data")
results_path <- file.path("results", "simulated-data-archR-v0.1.8-bound-1e08", paste0(ncores, "cores"))
if(!dir.exists(results_path)){
    message("Creating dir ", results_path)
    dir.create(results_path)
}else{
    message("Dir exists, not creating a new one ", results_path)
}


file_names_list <- list(
    n1000 = list.files(data_path, pattern = "*n1000.fa$",
        recursive = FALSE, full.names = TRUE),
    "n5000" = list.files(data_path, pattern = "*mu0.1_p1_n5000.fa$",
      recursive = FALSE, full.names = TRUE),
    "n10000" = list.files(data_path, pattern = "*mu0.1_p1_n10000.fa$",
      recursive = FALSE, full.names = TRUE)#,
    # "n25000" = list.files(data_path, pattern = "*mu0.1_p1_n25000.fa$",
    #   recursive = FALSE, full.names = TRUE),
    # "n50000" = list.files(data_path, pattern = "*mu0.1_p1_n50000.fa$",
    #                 recursive = FALSE, full.names = TRUE)
)
file_names_list$n1000 <- file_names_list$n1000[1]
message("Files to process")
print(file_names_list)


##----------------------------RUN EXPS IN A FUNCTION ---------------------------
run_archR_experiment <- function(fname, dir_name, seed_val, results_df,
    useNCores = 2, chunkSizes, tsv_fname = NULL){

    message("Using seed: ", seed_val)
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



    # chunkSizes <- switch(as.character(nSeqs_this),
    #     "1000"= chunkSizes_list$n1000,
    #     "5000"= chunkSizes_list$n5000,
    #     "10000"= chunkSizes_list$n10000,
    #     "25000" = chunkSizes_list$n25000,
    #     "50000" = chunkSizes_list$n50000
    # )

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
                    "mem",
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

                ### Using Rprof ============
                p_fname <- file.path(dir_name, paste0(param_str, "_", chunk_size, "_",
                    "Rprof.out"))
                utils::Rprof(filename = p_fname, append = FALSE,
                    memory.profiling = TRUE)

                # Call/Run archR
                perform_iters <- 2
                # if(nSeqs_this == 1000){
                    # cat(paste0("Calling archR when nSeqs = ", nSeqs_this, "\n"))
                    archRresult <- archR::archR(config = archRconfig,
                        seqs_ohe_mat = inputSeqsMat,
                        seqs_raw = inputSeqsRaw,
                        seqs_pos = positions,
                        total_itr= perform_iters,
                        set_ocollation = use_collation,
                        o_dir = result_dir_path
                    )
                # }
                utils::Rprof(NULL)
                profile <- utils::summaryRprof(filename = p_fname,
                    chunksize = -1L,
                    memory = "tseries", diff = FALSE)
                print(head(profile))
                max_mem <- max(rowSums(profile[,1:3]))/(1024*1024)
                # ^ max mem recorded in MB
                print("====================================")
                message("Maximum memory usage: ", max_mem, " MB")
                print("====================================")
                archRresult <- readRDS(file.path(result_dir_path, "archRresult.rds"))
                ##
                # check_mutrate <- get_mut_rate_from_param_str(param_str)
                # check_mutpos <- get_mut_pos_from_param_str(param_str)
                ##
                relevant_idx <- dplyr::filter(results_df,
                    nSeqs == nSeqs_this &
                        parallel == parallelize &
                        chunkSize == chunk_size &
                        seed == seed_val
                )$index
                print(result_dir_path)
                print(head(results_df))
                message("Number of rows to fill: ", length(relevant_idx))
                message("Filling row: ", relevant_idx)
                stopifnot(length(relevant_idx) == 1)
                ##
                results_df[relevant_idx, "MaxRSS"] <- max_mem
                ##
            } ## parallelize for loop ends
        } ## bound_values for loop ends
    } ## chunkSizes for loop ends
    ##
    ## write updated results_df
    write.table(results_df,
        # file.path(results_path, "archR_result_simulated_data_summary_mem3.tsv"),
        tsv_fname,
        sep = "\t", col.names = TRUE, row.names = FALSE)
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
seed_values <- sample.int(.Machine$integer.max, size = 10, #nRuns,
                            replace = FALSE)

choose_seed <- chooserun
seed_values <- seed_values[choose_seed]

tsv_fname <- file.path(results_path,
    paste0("archR_result_simulated_data_summary_mem", choose_seed, ".tsv"))


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


results_df_1000 <- purrr::cross_df(list(mutrate = c(0.1),
    mutpos = c(1),
    seed = seed_values,
    parallel = parallelizeVals,
    nSeqs = 1000,
    chunkSize = chunkSizes_list$n1000,
    MaxRSS = -1.0
))
## ^ creates sufficient rows for n1000

## for n5000 and n10000, only one combination of mutrate and mutpos
results_df_5000 <- purrr::cross_df(list(mutrate = c(0.1),
    mutpos = c(1),
    seed = seed_values,
    parallel = parallelizeVals,
    nSeqs = 5000,
    chunkSize = chunkSizes_list$n5000,
    MaxRSS = -1.0
))


results_df_10000 <- purrr::cross_df(list(mutrate = c(0.1),
    mutpos = c(1),
    seed = seed_values,
    parallel = parallelizeVals,
    nSeqs = 10000,
    chunkSize = chunkSizes_list$n10000,
    MaxRSS = -1.0
))


# results_df_25000 <- purrr::cross_df(list(mutrate = c(0.1),
#     mutpos = c(1),
#     seed = seed_values,
#     parallel = parallelizeVals,
#     nSeqs = 25000,
#     chunkSize = chunkSizes_list$n25000,
#     MaxRSS = -1.0
# ))


# results_df_50000 <- purrr::cross_df(list(mutrate = c(0.1),
#     mutpos = c(1),
#     seed = seed_values,
#     parallel = parallelizeVals,
#     nSeqs = 50000,
#     chunkSize = chunkSizes_list$n50000,
#     MaxRSS = -1.0
# ))


# results_df <- results_df_1000
results_df <- tibble::add_row(results_df_1000, results_df_5000)
results_df <- tibble::add_row(results_df, results_df_10000)
# results_df <- tibble::add_row(results_df, results_df_25000)
# results_df <- tibble::add_row(results_df, results_df_50000)

results_df <- tibble::add_column(results_df, .before = "mutrate",
    index = seq_len(nrow(results_df)))

file_names_list <- as.vector(unlist(file_names_list))
print(file_names_list)
print(results_df)

for(seed_val in seed_values){
    #
    dir_name <- file.path(results_path, paste0("seed_", seed_val))
    if(!dir.exists(dir_name)){
        message("Creating directory: ", dir_name)
        dir.create(dir_name)
    }
    #
    for(file_name in file_names_list){
        nSeqsStr <- unlist(strsplit(
                        unlist(strsplit(file_name, split = "_"))[3],
                        split = ".fa"))[1]
        use_chunk_sizes <- chunkSizes_list[[nSeqsStr]]
        message("Using chunk sizes: ", paste(use_chunk_sizes))
        for(cs in use_chunk_sizes){
            results_df <- run_archR_experiment(file_name, dir_name, seed_val,
                    results_df, useNCores = ncores, chunkSizes = cs,
                tsv_fname = tsv_fname)
        }
    }
}



write.table(results_df, tsv_fname, sep = "\t", col.names = TRUE, row.names = FALSE)
