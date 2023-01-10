#!/usr/bin/srun /bin/bash
#SBATCH -c 48
#SBATCH -N 1
#SBATCH -J simdata
#SBATCH --output=simdata_archRv0.1.8_slurm_%j_create_tsv_mem.out
#SBATCH --mem 10G
#SBATCH -p daniocode
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=s.nikumbh@lms.mrc.ac.uk


# archR_outputlog_simdata_slurm_%j.out
# sequence of arguments: for Rscript conda-environment-to-use path-to-python-to-use num-of-cores num-of-runs
# conda-environment: this should be the one where required python modules are installed/available
# 
# path-to-python: the python installation to use. this is important if you want a particular python to be used
#
# num-of-cores: this should be the same as that requested from/allocated by slurm (should match -c option to slurm)
# When not running for memory usage accounting, the fourth arg is
#  num-of-runs: this is the number of repeated runs of experiments to be performed
# Otherwise, it is
#  choose-run/seed (out of 10 runs/seeds used for running timing/performance experiments)

Rscript archR-on-simulated-data-memory-usage.R r-reticulate /mnt/biggley/home/sarvesh/sarvesh_miniconda/bin/python $SLURM_CPUS_PER_TASK 6


