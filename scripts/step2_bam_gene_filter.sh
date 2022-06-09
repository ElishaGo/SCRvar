#!/bin/bash
# A script to run step 2 in the pipeline:
# Run a bam file with htseq-count in order to filter out reads that are not on genes.


FILTERED_BAM_PATH=$1  # bam to filter
OUTPUT_DIR=$2  # folder for outputs
ANNOTATION_GTF=$3  # path to genecode gtf file
EDITING_GTF_INTERSECT=$4  # path to intersection of editing DB (only A bases) and gtf files
SNP_GTF_INTERSECT=$5  # path to intersection of SNP and gtf files
SNAME=$6  # sample name
N=$7  # number of threads

OUTPUT_FILE=${OUTPUT_DIR}/2.${SNAME}.gene_filter
LOGFILE=${OUTPUT_DIR}/step2.log

(
#    # run htseq to filter non gene reads
#    htseq-count -f bam -i gene_name -t gene -m union -s yes -n ${N} -o ${OUTPUT_FILE}.sam ${FILTERED_BAM_PATH} ${ANNOTATION_GTF} 3>&1 > ${OUTPUT_DIR}/2.${SNAME}_read_counts.txt
#
#    # add header to the sam file
#    samtools view -H ${FILTERED_BAM_PATH} | cat - ${OUTPUT_FILE}.sam > ${OUTPUT_FILE}_temp.sam
#
#    # replace reheadered (temp) file with previous file
#    rm ${OUTPUT_FILE}.sam
#    mv ${OUTPUT_FILE}_temp.sam ${OUTPUT_FILE}.sam
#
#    # get statistics on file
#    samtools flagstat -@ ${N} ${OUTPUT_FILE}.sam > ${OUTPUT_DIR}/2.flagstat_htseq.tsv
#
#    # keep only gene sites
#    grep -v "__" ${OUTPUT_FILE}.sam | samtools view -@ ${N} -Sb - > ${OUTPUT_FILE}.bam;samtools sort -@ ${N} ${OUTPUT_FILE}.bam -o ${OUTPUT_FILE}.bam;samtools index ${OUTPUT_FILE}.bam
#
#    # remove sam files
#    rm ${OUTPUT_FILE}.sam
#
#    # make intersection between editing A_I sites-transcriptomes-filtered bam file
#    bedtools intersect -u -header -a ${EDITING_GTF_INTERSECT} -b ${OUTPUT_FILE}.bam > ${OUTPUT_DIR}/2.${SNAME}.editing.genecode.bam_intersect.bed
#
#    # make intersection between SNP-transcriptomes-filtered bam file
#    bedtools intersect -u -header -a ${SNP_GTF_INTERSECT} -b ${OUTPUT_FILE}.bam > ${OUTPUT_DIR}/2.${SNAME}.snp.genecode.bam_intersect.vcf

    # make venn diagram for editing and bam file
    EDITING_GTF_INTERSECT_count=$(wc -l < ${EDITING_GTF_INTERSECT})
    output_bam_count=$(samtools view -c --threads ${N} ${OUTPUT_FILE}.bam)
    EDITING_GTF_INTERSECT_output_bam_intersect=$(wc -l < ${OUTPUT_DIR}/2.${SNAME}.editing.genecode.bam_intersect.bed)
    python -c "import sys;sys.path.append('/home/labs/bioservices/shared/rarevar/scrarevar');import sc_rna_variants.statistic_plots as sp; sp.plot_venn2_diagram((${EDITING_GTF_INTERSECT_count},${output_bam_count},${EDITING_GTF_INTERSECT_output_bam_intersect}), ('editing_gtf_intersect','bam_gtf_intersect', ''), '${OUTPUT_DIR}/2.editing_gtf.bam_gtf.venn.png', '${SNAME}')"

    # make venn diagram for snp and bam file
    SNP_GTF_INTERSECT_count=$(wc -l < ${SNP_GTF_INTERSECT})
    SNP_GTF_INTERSECT_output_bam_intersect=$(wc -l < ${OUTPUT_DIR}/2.${SNAME}.snp.genecode.bam_intersect.vcf)
    python -c "import sys;sys.path.append('/home/labs/bioservices/shared/rarevar/scrarevar');import sc_rna_variants.statistic_plots as sp; sp.plot_venn2_diagram((${SNP_GTF_INTERSECT_count},${output_bam_count},${SNP_GTF_INTERSECT_output_bam_intersect}), ('snp_gtf_intersect','bam_gtf_intersect', ''), '${OUTPUT_DIR}/2.snp_gtf.bam_gtf.venn.png', '${SNAME}')"

    echo end_of_log_file 1>&2 # test stderr

)
#  >& $LOGFILE
