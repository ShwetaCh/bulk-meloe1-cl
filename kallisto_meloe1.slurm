#!/bin/bash
#SBATCH --job-name=star                # Job name
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=20              # Number of CPU cores per task
#SBATCH --mem=32GB                      # Total memory limit
#SBATCH --time=20:00:00                 # Time limit hrs:min:sec
#SBATCH --output=kallisto.meloe1_%j.log      # Standard output and error log


####### How to call this Script ##########
#
# sbatch script.slurm <SAMPLE> <OUT DIR>
#
##########################################

SAMPLE=$1
OUT_DIR=/mnt/scratch/schavan/Yale/RUN_Meloe1_celllines
COND=NA
file_prefix=$SAMPLE.$COND
num_threads=10

FQ_TR_DIR=/mnt/scratch/schavan/Yale/trimmed_fastqs #dir with/to hold final trimmed fastq files
mkdir -p $FQ_TR_DIR

#Kallisto indexes
KALLISTO_IDX=/mnt/galaxy/home/schavan/projects/bulk_rna_seq/HPVDetection/HPV-Hybrid/Hybrid_Kallisto_CustomMeloe1/transcriptsHPV_Meloe1.idx

source activate rnaseq

#if [ -f ${OUT_DIR}/${SAMPLE}*.fastq.gz ]; then
mkdir -p ${OUT_DIR}/QUANT/${SAMPLE}.${COND}
kallisto quant -i $KALLISTO_IDX -o ${OUT_DIR}/QUANT/${SAMPLE}.${COND} --bias -b 100 -t 8 ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz &> ${OUT_DIR}/QUANT/${SAMPLE}.${COND}/${SAMPLE}.${COND}.kallisto.log   # Align reads to reference
    
    #Clean up--->
#fi
