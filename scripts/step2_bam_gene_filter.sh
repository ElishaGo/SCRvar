#!/bin/bash
# A script to run step 2 in the pipeline:
# prepare and run a bam file with htseq-count, in order to filter out reads that are not on genes.

FILTERED_BAM_PATH=$1
ANNOTATION_GTF=$2
FNAME=$3
N=$4


mkdir step2_bam_gene_filter
# add chr to chromosome names in bam files
samtools view -H ${FILTERED_BAM_PATH} | sed  -e '/SN:chr/!s/SN:\([0-9XY]*\)/SN:chr&/' -e '/SN:chrM/!s/SN:MT/SN:chrM&/' | samtools reheader - ${FILTERED_BAM_PATH} > step1_filtered_bam_files/1_${FNAME}_CBfiltered_chr.bam;samtools index step1_filtered_bam_files/1_${FNAME}_CBfiltered_chr.bam

# remove old bam files and index
rm ${FILTERED_BAM_PATH}*

# run htseq to filter non gene reads
htseq-count -f bam -i gene_name -t gene -m union -s yes -o step2_bam_gene_filter/2_${FNAME}.gene_filter.sam step1_filtered_bam_files/1_${FNAME}_CBfiltered_chr.bam ${ANNOTATION_GTF} 3>&1 > step2_bam_gene_filter/2_${FNAME}_stranded_counts.txt

# add header to the sam file
samtools view -H step1_filtered_bam_files/1_${FNAME}_CBfiltered_chr.bam | cat - step2_bam_gene_filter/2_${FNAME}.gene_filter.sam > step2_bam_gene_filter/2_${FNAME}.gene_filter_header.sam

# get statistics on file
samtools flagstat -@ ${N} step2_bam_gene_filter/2_${FNAME}.gene_filter_header.sam > step2_bam_gene_filter/2_flagstat_htseq.tsv

# keep only gene sites
grep -v "__" step2_bam_gene_filter/2_${FNAME}.gene_filter_header.sam | samtools view -@ ${N} -Sb - > step2_bam_gene_filter/2_${FNAME}.gene_filter_header.bam;samtools sort -@ ${N} step2_bam_gene_filter/2_${FNAME}.gene_filter_header.bam -o step2_bam_gene_filter/2_${FNAME}.gene_filter_header.bam;samtools index step2_bam_gene_filter/2_${FNAME}.gene_filter_header.bam

# remove sam files
rm step2_bam_gene_filter/2_${FNAME}.gene_filter_header.sam
rm step2_bam_gene_filter/2_${FNAME}.gene_filter.sam
