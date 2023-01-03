#!/bin/bash
#SBATCH --job-name=qc_merge_trim_hervs                # Job name
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=20              # Number of CPU cores per task
#SBATCH --mem=32GB                      # Total memory limit
#SBATCH --time=10:00:00                 # Time limit hrs:min:sec
#SBATCH --output=qc_merge_trim_%j.log      # Standard output and error log


####### How to call this Script ##########
#
# sbatch qc_merge_trim_hervs.slurm <SAMPLE> <OUT DIR>
# ex: sbatch qc_merge_trim_hervs.slurm YUTAS_p /mnt/scratch/schavan/Yale/RUN_herv_090122
#
##########################################

SAMPLE=$1
OUT_DIR=/mnt/scratch/schavan/Yale/RUN_Meloe1_celllines
COND=NA

IN_DIR=/mnt/scratch/schavan/Yale/raw_fastqs #raw fastq dir
FQ_TR_DIR=/mnt/scratch/schavan/Yale/trimmed_fastqs #dir with/to hold final trimmed fastq files
mkdir -p $FQ_TR_DIR

source activate rnaseq

#if [ -f ${IN_DIR}/${SAMPLE}*.fastq.gz ]; then
    
    #  Merge lanes--->

    cat $IN_DIR/$SAMPLE*R1* > ${OUT_DIR}/${SAMPLE}.${COND}.merged.R1.fastq.gz
    cat $IN_DIR/$SAMPLE*R2* > ${OUT_DIR}/${SAMPLE}.${COND}.merged.R2.fastq.gz
    
    # Run fastp--->

    fastp -h ${OUT_DIR}/${SAMPLE}.${COND}.fastp.html -j ${OUT_DIR}/${SAMPLE}.${COND}.fastp.json \
        -i ${OUT_DIR}/${SAMPLE}.${COND}.merged.R1.fastq.gz -I ${OUT_DIR}/${SAMPLE}.${COND}.merged.R2.fastq.gz \
        -o ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz -O ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz \
        -y -3 --cut_tail_window_size 4 -5 --cut_front_window_size 4 -l 40 -n 5 -p

    # Run fastqc--->

    fastqc ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz --outdir $OUT_DIR

    # Clean up--->
#fi
