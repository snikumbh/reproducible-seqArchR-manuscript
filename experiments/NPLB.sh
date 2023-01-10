#!/bin/bash
#SBATCH -c 16
#SBATCH -N 1
#SBATCH -J comparecompare
#SBATCH --mem-per-cpu 10M
#SBATCH -p daniocode
### 1 is the fasta file
### 2 is the output dir name
INFILE=$1
OUTDIR=$2
###NPLB/NPLB/promoterLearn -f $INFILE -proc $SLURM_CPUS_PER_TASK -o $OUTDIR
