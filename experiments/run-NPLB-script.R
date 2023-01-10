
library(readr)

nRuns <- 10

get_nplb_cmd <- function(fasta_file, out_dir){
    NPLB_cmd <- paste0("comparison-approaches/NPLB/NPLB/promoterLearn -f ",
        fasta_file, " -proc $SLURM_CPUS_PER_TASK -o ", out_dir)
    return(NPLB_cmd)
}

write_initial_nplb_script <- function(fname = "NPLB.sh", ncores = 8, memPerCpu = "10M",
                                      partition = "high", useName = "compare"){

    readr::write_lines(c("#!/bin/bash",
    paste0("#SBATCH -c ", ncores),
    "#SBATCH -N 1",
    paste0("#SBATCH -J compare", useName),
    paste0("#SBATCH --mem-per-cpu ", memPerCpu),
    paste0("#SBATCH -p ", partition),
    "### 1 is the fasta file",
    "### 2 is the output dir name",
    "INFILE=$1",
    "OUTDIR=$2",
    "###NPLB/NPLB/promoterLearn -f $INFILE -proc $SLURM_CPUS_PER_TASK -o $OUTDIR"),
    file = fname, sep="\n", append = FALSE)
}

##--------------------------- SETUP FILES --------------------------------------
data_path <- file.path("data", "simulated-data")
results_path <- file.path("results", "comparison-approaches", "simulated-data")

# Load fasta files
# -- For n1000, all combinations of mutation rates and mutation positions
# -- For n5000, only mu = 0.1 and p = 1
# -- For n10000, only mu = 0.1 and p = 1
file_names_list <- list("n1000" = list.files(data_path, pattern = "*n1000.fa$",
    recursive = FALSE, full.names = TRUE),
    "n5000" = list.files(data_path, pattern = "*mu0.1_p1_n5000.fa$",
        recursive = FALSE, full.names = TRUE),
    "n10000" = list.files(data_path, pattern = "*mu0.1_p1_n10000.fa$",
        recursive = FALSE, full.names = TRUE))

message("Files to process")
print(file_names_list)


create_nplb_experiment_script <- function(fname, dir_name, seed_val, nplb_sh){
    ## Set seed, use seed_val
    ## Read FASTA file, use file_name
    ## change to randomize order
    ## write changed order file to disk for reading as input NPLB, use dir_name
    ## Add a line in slurm sbatch script
    ## Run slurm script

    set.seed(seed_val)
    # message("*****Now processing*****: ", fname)
    param_str <- tail(unlist(strsplit(strsplit(fname, split=".fa")[[1]], "-")),1)
    #
    inputSeqsRaw <- Biostrings::readDNAStringSet(fname, format = "fasta",
                                                use.names = TRUE)

    nSeqs_this <- length(inputSeqsRaw)

    # Randomize the sequence order
    changedOrder <- sample.int(nSeqs_this, nSeqs_this, replace = FALSE)
    inputSeqsRaw <- inputSeqsRaw[changedOrder]

    ##
    result_dir <- file.path(dir_name, paste("nplb_result_simulateddata",
                                param_str, sep = "_"))
    fname_randomized <- file.path(dir_name,
        paste0(basename(fname), "_", seed_val, "_randomized"))

    Biostrings::writeXStringSet(inputSeqsRaw, format = "fasta",
        filepath = fname_randomized)


    add_line <- get_nplb_cmd(fname_randomized, result_dir)
    echo_line <- paste("echo ", "*****Now processing*****: ", fname)
    write(paste0(echo_line, "\n", add_line, "\n#\n"), nplb_sh, append = TRUE)
}

summarize_nplb_experiment <- function(fname, dir_name, seed_val, results_df){

    param_str <- tail(unlist(strsplit(strsplit(fname, split=".fa")[[1]], "-")),1)

    ##
    result_dir <- file.path(dir_name, paste("nplb_result_simulateddata",
        param_str, sep = "_"))
    ##
    # print(result_dir)
    nplbArchFile <- file.path(result_dir, "architectureDetails.txt")
    if(!file.exists(nplbArchFile)){
        # message(nplbArchFile, " NOT FOUND")
        return(results_df)
    }
    ##
    nplbArchDetails <- read.delim(nplbArchFile, sep="\t", header=FALSE)
    nSeqs_this <- nrow(nplbArchDetails)
    ##
    check_mutrate <- get_mut_rate_from_param_str(param_str)
    check_mutpos <- get_mut_pos_from_param_str(param_str)
    ##
    relevant_idx <- dplyr::filter(results_df,
            nSeqs == nSeqs_this &
            parallel == TRUE &
            mutrate == check_mutrate &
            mutpos == check_mutpos &
            seed == seed_val
    )$index
    # print(relevant_idx)
    # message("Filling row: ", relevant_idx)

    stopifnot(length(relevant_idx) == 1)
    ##
    startTime <- file.mtime(file.path(result_dir, "settings.txt"))
    endTime <- file.mtime(file.path(result_dir, "architectureDetails.txt"))
    dur <- as.numeric(difftime(endTime, startTime, units ="mins"))
    results_df[relevant_idx, "runTime_m"] <- dur
    ##
    ##
    trueLabels <- unlist(lapply(as.vector(nplbArchDetails[,2]), function(x){
        unlist(strsplit(x, "_clust"))[2]
    }))
    # print(trueLabels)
    # print(table(trueLabels))
    ##
    ## Get cluster labels from NPLB result files
    predLabels <- as.character(nplbArchDetails[,1])
    # print(predLabels)
    # print(table(predLabels))
    # print(mclust::adjustedRandIndex(predLabels, trueLabels))
    # print(mclust::adjustedRandIndex(trueLabels, predLabels))
    # stop("samarth")
    ##
    results_df[relevant_idx, "ARI"] <-
        mclust::adjustedRandIndex(predLabels, trueLabels)
    ##
    ##
    results_df[relevant_idx, "nClusters"] <- length(unique(nplbArchDetails[,1]))

    ## write updated results_df
    write.table(results_df, file.path(results_path, "nplb_result_simulated_data_summary.tsv"),
        sep = "\t", col.names = TRUE, row.names = FALSE)
    ##
    return(results_df)
}



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
params_sett_list <- lapply(file_names_list$n1000, function(x){
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
    parallel = TRUE,
    nSeqs = 1000,
    runTime_m = -1.0,
    ARI = -1.0,
    nClusters = -1
))
## ^ creates sufficient rows for n1000
#
## for n5000 and n10000, only one combination of mutrate and mutpos
results_df_5000 <- purrr::cross_df(list(mutrate = c(0.1),
    mutpos = c(1),
    seed = seed_values,
    parallel = TRUE,
    nSeqs = 5000,
    runTime_m = -1.0,
    ARI = -1.0,
    nClusters = -1
))
#
#
results_df_10000 <- purrr::cross_df(list(mutrate = c(0.1),
    mutpos = c(1),
    seed = seed_values,
    parallel = TRUE,
    nSeqs = 10000,
    runTime_m = -1.0,
    ARI = -1.0,
    nClusters = -1
))
#
results_df <- tibble::add_row(results_df_1000, results_df_5000)
results_df <- tibble::add_row(results_df, results_df_10000)
#
results_df <- tibble::add_column(results_df, .before = "mutrate",
    index = seq_len(nrow(results_df)))

file_names_list <- as.vector(unlist(file_names_list))
print(file_names_list)
# file_names_list <- file_names_list[13:17]

nplb_sh <- "NPLB.sh"
## create a new NPLB.sh script
## First lines will be written
write_initial_nplb_script(fname = nplb_sh, ncores = 16, memPerCpu = "10M",
                          partition = "daniocode", useName = "compare")

for(seed_val in seed_values){
    #
    dir_name <- file.path(results_path, paste0("seed_", seed_val))
    # if(!dir.exists(dir_name)){
    #     dir.create(dir_name)
    # }
    # #
    # for(file_name in file_names_list){
    #     create_nplb_experiment_script(file_name, dir_name, seed_val, nplb_sh)
    # }
}

## results/experiments summary
for(seed_val in seed_values){
    #
    dir_name <- file.path(results_path, paste0("seed_", seed_val))
    # seed dir exists
    #
    for(file_name in file_names_list){
        results_df <-
            summarize_nplb_experiment(file_name, dir_name, seed_val, results_df)
    }
}


write.table(results_df, file.path(results_path, "nplb_result_simulated_data_summary.tsv"),
    sep = "\t", col.names = TRUE, row.names = FALSE)
