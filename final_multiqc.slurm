#!/bin/bash
#SBATCH --job-name=multiqc                # Job name
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=100              # Number of CPU cores per task
#SBATCH --mem=32GB                      # Total memory limit
#SBATCH --time=10:00:00                 # Time limit hrs:min:sec
#SBATCH --output=multiqc_%j.log      # Standard output and error log


####### How to call this Script ##########
#
# sbatch script.slurm
#
##########################################
$OUT_DIR=/mnt/scratch/schavan/Yale/RUN_Meloe1_celllines

source activate rnaseq

# Compile multiqc report from all the fastp, fastqc, and salmon logs--->

     multiqc -f $OUT_DIR ${OUT_DIR}/QUANT/

# Clean up--->
