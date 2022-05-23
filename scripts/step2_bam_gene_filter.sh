#!/bin/bash

# A script to run step 2 in the pipeline:
# prepare and run a bam file with htseq-count, in order to filter out reads that are not on genes.

INPUT_DIR=$1
FILTERED_BAM_PATH=$2
OUTPUT_DIR=$3
ANNOTATION_GTF=$4
EDITING_GTF_INTERSECT=$5
SNP_GTF_INTERSECT=$6
FNAME=$7
N=$8

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

# make intersection between editing A_I sites transcriptomes and filtered bam file
bedtools intersect -u -header -a EDITING_GTF_INTERSECT -b 2_${FNAME}.gene_filter_header.bam > ${OUTPUT_DIR}/2_${FNAME}_editing.bam_intersect.bed6

# make intersection between SNP A_I sites transcriptomes and filtered bam file
bedtools intersect -u -header -a SNP_GTF_INTERSECT -b 2_${FNAME}.gene_filter_header.bam > ${OUTPUT_DIR}/2_${FNAME}_snp.bam_intersect.vcf
