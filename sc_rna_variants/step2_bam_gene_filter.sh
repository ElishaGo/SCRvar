#!/bin/bash
# Retrieve reads only on genes by usage of HTseq-count
# create venn diagrams with overlaps between the data and SNP and editing sites databases if supplied

export FILTERED_BAM_PATH=$1  # bam to filter
export OUTPUT_DIR=$2  # folder for outputs
export ANNOTATION_GTF=$3  # path to genecode gtf file
export EDITING_GTF_INTERSECT=$4  # path to intersection of editing DB (only A bases) and gtf files
export SNP_GTF_INTERSECT=$5  # path to intersection of SNP and gtf files
export EDITING_GTF_BAM_INTERSECT=$6  # path to intersection of editing DB, gtf and bam files
export SNP_GTF_BAM_INTERSECT=$7  # path to intersection of SNP, gtf and bam files
export SNAME=$8  # sample name
export THREADS=$9  # number of threads

export OUTPUT_FILE=${OUTPUT_DIR}/2.${SNAME}.gene_filter

# redirect stdout/stderr to a file
LOGFILE=${OUTPUT_DIR}/step2.log
exec > $LOGFILE 2>&1

#{
# run HTseq to annotate non-gene reads
BAM_DIR=$(dirname $FILTERED_BAM_PATH)
bamtools split -in ${FILTERED_BAM_PATH} -reference  # split bam for parallelism
OUT_NAMES=$(for f in $BAM_DIR/*REF_chr*.bam; do echo -n "-o ${f}_temp_htseq_out.bam "; done)
for f in $BAM_DIR/*REF_chr*.bam; do samtools index "${f}"; done
htseq-count -f bam -i gene_name -t gene -m union -s yes -n ${THREADS} -p bam ${OUT_NAMES} -c ${OUTPUT_DIR}/2.${SNAME}_gene_counts_htseq.tsv $BAM_DIR/*REF_chr*.bam ${ANNOTATION_GTF}

# keep only gene sites reads and save into a bam file
samtools merge -c -@ ${THREADS} -O SAM - $BAM_DIR/*REF*temp_htseq_out.bam | grep -v "__" - | samtools sort -@ ${THREADS} - -o ${OUTPUT_FILE}.bam
samtools index ${OUTPUT_FILE}.bam

rm ${BAM_DIR}/*REF*.bam

# get statistics on filtered file
samtools flagstat -@ ${THREADS} ${OUTPUT_FILE}.bam > ${OUTPUT_DIR}/2.${SNAME}_flagstat.txt
# RUN INTERSECTIONS WITH EDITING DB ANALYSIS
if [ $EDITING_GTF_INTERSECT != "None" ]
then
  # make intersection between editing A_I sites-transcriptomes-filtered bam file
  bedtools intersect -u -header -sorted -a ${EDITING_GTF_INTERSECT} -b ${OUTPUT_FILE}.bam > ${EDITING_GTF_BAM_INTERSECT}

  # make venn diagram for editing and bam file
  bedtools genomecov -max 1 -split -ibam ${OUTPUT_FILE}.bam > ${OUTPUT_DIR}/2.${SNAME}_genomecov.txt
  bam_positions_count=$(tail -n 1 ${OUTPUT_DIR}/2.${SNAME}_genomecov.txt | awk '{print $3}')
  EDITING_GTF_INTERSECT_count=$(grep -v "#" ${EDITING_GTF_INTERSECT} | cut -f 1,2,3 | sort |  uniq | wc -l)
  EDITING_GTF_INTERSECT_output_bam_intersect=$(grep -v "#" ${EDITING_GTF_BAM_INTERSECT} | cut -f 1,2,3 | sort |  uniq | wc -l)
  python -c "import sys;sys.path.append('/home/labs/bioservices/shared/rarevar/scrarevar');import sc_rna_variants.statistic_plots as sp; sp.plot_venn2_diagram({'10':${EDITING_GTF_INTERSECT_count}-${EDITING_GTF_INTERSECT_output_bam_intersect},'01':${bam_positions_count}-${EDITING_GTF_INTERSECT_output_bam_intersect},'11':${EDITING_GTF_INTERSECT_output_bam_intersect}}, ('editing DB positions on genes','aligned positions', ''), '${OUTPUT_DIR}/2.editing_gtf.bam_gtf.venn.png', 'count of positions on genes - ${SNAME}')"
fi

# RUN INTERSECTIONS WITH SNP DB ANALYSIS
if [ $SNP_GTF_INTERSECT != "None" ]
then
  # make intersection between SNP-transcriptomes-filtered bam file
  bedtools intersect -u -header -sorted -a ${SNP_GTF_INTERSECT} -b ${OUTPUT_FILE}.bam > ${SNP_GTF_BAM_INTERSECT}

  # make venn diagram for snp and bam file
  SNP_GTF_INTERSECT_count=$(grep -v "#" ${SNP_GTF_INTERSECT} | cut -f 1,2 | sort | uniq | wc -l)
  SNP_GTF_INTERSECT_output_bam_intersect=$(grep -v "#" ${SNP_GTF_BAM_INTERSECT} | cut -f 1,2 | sort | uniq | wc -l)
  python -c "import sys;sys.path.append('/home/labs/bioservices/shared/rarevar/scrarevar');import sc_rna_variants.statistic_plots as sp; sp.plot_venn2_diagram({'10':${SNP_GTF_INTERSECT_count}-${SNP_GTF_INTERSECT_output_bam_intersect},'01':${bam_positions_count}-${SNP_GTF_INTERSECT_output_bam_intersect},'11':${SNP_GTF_INTERSECT_output_bam_intersect}}, ('SNP DB positions on genes','aligned positions', ''), '${OUTPUT_DIR}/2.snp_gtf.bam_gtf.venn.png', 'count of positions on genes - ${SNAME}')"
fi

echo end_of_log_file 1>&2

#) > $LOGFILE 2>&1