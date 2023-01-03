#!/bin/bash
##conda activate rnaseq
##~/get_a_worker_node.sh ##Don't forget to get a worker node first!!! Do not run this on the main node


IN_SS=/mnt/scratch/schavan/Yale/RUN_Meloe1_Celllines/SampleSheet.txt
OUT_DIR=/mnt/scratch/schavan/Yale/RUN_Meloe1_Celllines #output dir
IN_DIR=/mnt/scratch/schavan/Yale/raw_fastqs #raw fastq dir
FQ_TR_DIR=/mnt/scratch/schavan/Yale/trimmed_fastqs #dir with/to hold final trimmed fastq files

#For Kallisto
KALLISTO_IDX=/mnt/galaxy/home/schavan/projects/bulk_rna_seq/HPVDetection/HPV-Hybrid/Hybrid_Kallisto_CustomMeloe1/transcriptsHPV_Meloe1.idx

#What all needs to be run? Turn these flags on/off accordingly
DO_MERGE_FQ=1
DO_FASTP=1
DO_FASTQC=1
DO_KALLISTO=1
DO_MULTIQC=1

##Edit block above to point to appropriate locations
##-----------------------------------------------------------------------------------------------------------------------

echo "Input dir: "$IN_DIR
echo "Output dir: "$OUT_DIR
date

while IFS=$'\t' read line
do
    set $line
    SAMPLE=$1
    COND=$2
    echo "Started working on: "$SAMPLE $COND

if [[ $DO_MERGE_FQ -eq 1 ]]
then
    echo "Running merge--->"
    #Merge Lanes
    cat $IN_DIR/$SAMPLE*R1* > ${OUT_DIR}/${SAMPLE}.${COND}.merged.R1.fastq.gz
    cat $IN_DIR/$SAMPLE*R2* > ${OUT_DIR}/${SAMPLE}.${COND}.merged.R2.fastq.gz
fi

if [[ $DO_FASTP -eq 1 ]]
then
  echo "Running fastp--->"
  #Run fastp
    fastp -h ${OUT_DIR}/${SAMPLE}.${COND}.fastp.html -j ${OUT_DIR}/${SAMPLE}.${COND}.fastp.json \
    -i ${OUT_DIR}/${SAMPLE}.${COND}.merged.R1.fastq.gz -I ${OUT_DIR}/${SAMPLE}.${COND}.merged.R2.fastq.gz \
    -o ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz -O ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz \
    -y -3 --cut_tail_window_size 4 -5 --cut_front_window_size 4 -l 40 -n 5 -p
fi

if [[ $DO_FASTQC -eq 1 ]]
then
    echo "Running FASTQC--->"
    #Run fastQC before
    #fastqc $OUT_DIR
    #Run fastQC after https://home.cc.umanitoba.ca/~psgendb/doc/fastqc.help
    #fastqc ${OUT_DIR}/${SAMPLE}.${COND}.merged.R1.fastq.gz ${OUT_DIR}/${SAMPLE}.${COND}.merged.R2.fastq.gz --outdir $OUT_DIR
    fastqc ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz --outdir $OUT_DIR
fi

if [[ $DO_KALLISTO -eq 1 ]]
then
    echo "Running Kallisto--->"
    #Run Kallisto
    mkdir -p ${OUT_DIR}/QUANT/${SAMPLE}.${COND}
    kallisto quant -i $KALLISTO_IDX -o ${OUT_DIR}/QUANT/${SAMPLE}.${COND} --bias -b 100 -t 8 ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R1.fastq.gz ${FQ_TR_DIR}/${SAMPLE}.${COND}.merged.trimmed.R2.fastq.gz &> ${OUT_DIR}/QUANT/${SAMPLE}.${COND}/${SAMPLE}.${COND}.kallisto.log
#--rf-stranded
#--fr-stranded
fi

if [[ $DO_MULTIQC -eq 1 ]]
then
    #Run MultiQC
    echo "Running MultiQC--->"
    multiqc -f $OUT_DIR ${OUT_DIR}/QUANT/
fi

done < <(tail -n +2 $IN_SS)
date

