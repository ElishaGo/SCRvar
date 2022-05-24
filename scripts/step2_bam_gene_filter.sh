#!/bin/bash

# A script to run step 2 in the pipeline:
# prepare and run a bam file with htseq-count, in order to filter out reads that are not on genes.

FILTERED_BAM_PATH=$1
OUTPUT_DIR=$2  # folder for outputs
ANNOTATION_GTF=$3  # path to genecode gtf file
EDITING_GTF_INTERSECT=$4  # path to intersection of editing DB (only A bases) and gtf files
SNP_GTF_INTERSECT=$5  # path to intersection of SNP and gtf files
FNAME=$6  # sample name
N=$7  # number of threads

LOGFILE=${OUTPUT_DIR}/step2.log

(
    # run htseq to filter non gene reads
    htseq-count -f bam -i gene_name -t gene -m union -s yes -o ${OUTPUT_DIR}/2.${FNAME}.gene_filter.sam ${FILTERED_BAM_PATH} ${ANNOTATION_GTF} 3>&1 > ${OUTPUT_DIR}/2.${FNAME}_stranded_counts.txt

    # add header to the sam file
    samtools view -H ${FILTERED_BAM_PATH} | cat - ${OUTPUT_DIR}/2.${FNAME}.gene_filter.sam > ${OUTPUT_DIR}/2.${FNAME}.gene_filter_header.sam

    # get statistics on file
    samtools flagstat -@ ${N} ${OUTPUT_DIR}/2.${FNAME}.gene_filter_header.sam > ${OUTPUT_DIR}/2.flagstat_htseq.tsv

    # keep only gene sites
    grep -v "__" ${OUTPUT_DIR}/2.${FNAME}.gene_filter_header.sam | samtools view -@ ${N} -Sb - > ${OUTPUT_DIR}/2.${FNAME}.gene_filter_header.bam;samtools sort -@ ${N} ${OUTPUT_DIR}/2.${FNAME}.gene_filter_header.bam -o ${OUTPUT_DIR}/2.${FNAME}.gene_filter_header.bam;samtools index ${OUTPUT_DIR}/2.${FNAME}.gene_filter_header.bam


    # remove sam files
    rm ${OUTPUT_DIR}/2.${FNAME}.gene_filter_header.sam
    rm ${OUTPUT_DIR}/2.${FNAME}.gene_filter.sam

    # make intersection between editing A_I sites transcriptomes and filtered bam file
    bedtools intersect -u -header -a ${EDITING_GTF_INTERSECT} -b ${OUTPUT_DIR}/2.${FNAME}.gene_filter_header.bam > ${OUTPUT_DIR}/2.${FNAME}_editing.bam_intersect.bed

    # make intersection between SNP A_I sites transcriptomes and filtered bam file
    bedtools intersect -u -header -a ${SNP_GTF_INTERSECT} -b ${OUTPUT_DIR}/2.${FNAME}.gene_filter_header.bam > ${OUTPUT_DIR}/2.${FNAME}_snp.bam_intersect.vcf

    echo end_of_log_file 1>&2 # test stderr

) >& $LOGFILE
