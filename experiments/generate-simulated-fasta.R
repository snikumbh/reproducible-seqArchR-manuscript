## Generate simulated data FASTA files




mutation_rates <- seq(0.1, 0.5, by=0.1)
mutate_positions <- 1:3

motif_files <- paste0(file.path("data", "simulated-data", "motifs.list"), seq(1,4))
nseqs <- c(200, 350, 150, 300)

nseqs_mult <- c(1, 5, 10, 25, 50)

ifToMut <- c(FALSE, TRUE, TRUE, TRUE)

py_str <- "python3.9"
prog_name <- "toy-seq-generation/generateData.py"
optionT <- "-T 1" ## motfis given in file
optionF <- "-F"  ## motif filename
optionN <- "-N"  ## number of sequences
optionI <- "-I 100"  ## minimum sequence length
optionX <- "-X 100"  ## maximum sequence length
optionS <- ""  ## start position
optionE <- ""  ## end position
optionM <- "-M"  ## mutation rate
optionP <- "-P"  ## mutation position
optionO <- "-O"  ## output filename
optionA <- "-A"  ## append str


use_file_path <- file.path("data", "simulated-data", "try_script_again_for_25_and_50k")
cmds_filename <- file.path(use_file_path, "cmds-generate-simulated-data_nseqs.sh")
cmds_file <- file(cmds_filename, "w")

print(mutation_rates)
print(mutate_positions)

for(nm in nseqs_mult){

  this_mutpos <- mutate_positions
  if(nm > 1) this_mutpos <- mutate_positions[1]

  this_mutrate <- mutation_rates
  if(nm > 1) this_mutrate <- mutation_rates[1]

  print(nm)
  print("--")
  print(this_mutpos)
  print(this_mutrate)
  print("--")
  print(mutation_rates)
  print(mutate_positions)

  for(mutrate in this_mutrate){
    for(mutpos in this_mutpos){
      filename <- file.path(use_file_path,
                    paste0("simulated-dataset-",
                        "mu", mutrate,
                        "_p", mutpos, "_n", (nm*sum(nseqs)),
                        ".fa"))
      info_str <- "## File generated using automated script. Do not edit by hand."
      cat(info_str, "\n\n", file = cmds_file)
      setting_str <- paste("## Mutation rate", mutrate,
                           "nPositionsMutated", mutpos,
                           "nSeqs", nm*sum(nseqs))
      ##
      cat(setting_str, "\n\n", file = cmds_file, append = TRUE)
      cat(setting_str, "\n\n")
      for(l in 1:length(motif_files)){
        if(!ifToMut[l]){
          use_mutrate <- 0.0
          use_mutpos <- 0
        }else{
          use_mutrate <- mutrate
          use_mutpos <- mutpos
        }
        fname_suffix <- paste0("clust", l, ".fa")
        cmd_str <- c(py_str, prog_name, optionT, optionF, motif_files[l],
                       optionN, nm*nseqs[l], optionI, optionX,
                       optionA, paste0("clust",l),
                       optionM, use_mutrate,
                       optionP, use_mutpos,
                       optionO, fname_suffix)

        cat(cmd_str, "\n", file = cmds_file, append = TRUE)
        cat(cmd_str, "\n")

      }
      cat_cmd_str <- paste("cat",
                           paste0("clust", seq(length(motif_files)), ".fa", collapse = " "),
                           ">", filename)
      cat(cat_cmd_str, "\n\n")
      cat(cat_cmd_str, "\n", file = cmds_file, append = TRUE)
      ## clean up
      rm_cmd_str <- paste("rm", paste0("clust", seq(length(motif_files)), ".fa", collapse = " "))
      cat(rm_cmd_str, "\n\n", file = cmds_file, append = TRUE)
    }
  }
}
close(cmds_file)
## Execute cmds to generate fasta files
system2("bash", args = cmds_filename, stdout = "", stderr = "", wait = TRUE)




