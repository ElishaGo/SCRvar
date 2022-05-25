#!/bin/bash

# A script to count number of reads in each barcode in bam file

BAM_FILE=$1
OUT_DIR=$2
SNAME=$3

samtools view ${BAM_FILE} | grep NH:i:1 | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > ${OUT_DIR}/raw_bam_reads_per_barcode_count.${SNAME}.csv