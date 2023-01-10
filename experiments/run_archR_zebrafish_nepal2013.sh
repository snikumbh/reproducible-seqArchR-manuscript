#!/usr/bin/srun /bin/bash
#SBATCH -c 48
#SBATCH -N 1
#SBATCH -J zf
#SBATCH --output=zf_archRv0.1.8_slurm_%j.out
#SBATCH --mem-per-cpu 10M
#SBATCH -p daniocode

# sequence of arguments: conda-environment-to-use path-to-python-to-use num-of-cores
# conda-environment: this should be the one where required python modules are installed/available
#
# path-to-python: the python installation to use. this is important if you want a particular python to be used
#
# num-of-cores: this should be the same as that requested from/allocated by slurm (should match -c option to slurm)


Rscript archR-on-zebrafish-nepal2013.R r-reticulate /mnt/biggley/home/sarvesh/sarvesh_miniconda/bin/python $SLURM_CPUS_PER_TASK

