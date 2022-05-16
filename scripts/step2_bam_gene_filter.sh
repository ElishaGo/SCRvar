#!/bin/bash

# A script to run step 2 in the pipeline:
# prepare and run a bam file with htseq-count, in order to filter out reads that are not on genes.

FILTERED_BAM_PATH=$1
OUTPUT_DIR=$2
ANNOTATION_GTF=$3
FNAME=$4
N=$5
INPUT_DIR=step1_filtered_bam_files

# add chr to chromosome names in bam files (step1)
samtools view -H ${FILTERED_BAM_PATH} | sed  -e '/SN:chr/!s/SN:\([0-9XY]*\)/SN:chr&/' -e '/SN:chrM/!s/SN:MT/SN:chrM&/' | samtools reheader - ${FILTERED_BAM_PATH} > ${INPUT_DIR}/1_${FNAME}_CBfiltered_chr.bam;samtools index ${INPUT_DIR}/1_${FNAME}_CBfiltered_chr.bam

# run htseq to filter non gene reads
htseq-count -f bam -i gene_name -t gene -m union -s yes -o ${OUTPUT_DIR}/2_${FNAME}.gene_filter.sam ${INPUT_DIR}/1_${FNAME}_CBfiltered_chr.bam ${ANNOTATION_GTF} 3>&1 > ${OUTPUT_DIR}/2_${FNAME}_stranded_counts.txt

# add header to the sam file
samtools view -H ${INPUT_DIR}/1_${FNAME}_CBfiltered_chr.bam | cat - ${OUTPUT_DIR}/2_${FNAME}.gene_filter.sam > ${OUTPUT_DIR}/2_${FNAME}.gene_filter_header.sam

# get statistics on file
samtools flagstat -@ ${N} ${OUTPUT_DIR}/2_${FNAME}.gene_filter_header.sam > ${OUTPUT_DIR}/2_flagstat_htseq.tsv

# keep only gene sites
grep -v "__" ${OUTPUT_DIR}/2_${FNAME}.gene_filter_header.sam | samtools view -@ ${N} -Sb - > ${OUTPUT_DIR}/2_${FNAME}.gene_filter_header.bam;samtools sort -@ ${N} ${OUTPUT_DIR}/2_${FNAME}.gene_filter_header.bam -o ${OUTPUT_DIR}/2_${FNAME}.gene_filter_header.bam;samtools index ${OUTPUT_DIR}/2_${FNAME}.gene_filter_header.bam

# remove old bam files and index
rm ${FILTERED_BAM_PATH}*

# remove sam files
rm ${OUTPUT_DIR}/2_${FNAME}.gene_filter_header.sam
rm ${OUTPUT_DIR}/2_${FNAME}.gene_filter.sam
