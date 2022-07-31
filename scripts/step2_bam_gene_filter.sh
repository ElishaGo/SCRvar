#!/bin/bash
# A script to run step 2 in the pipeline:
# Run a bam file with htseq-count in order to filter out reads that are not on genes.

export FILTERED_BAM_PATH=$1  # bam to filter
export OUTPUT_DIR=$2  # folder for outputs
ANNOTATION_GTF=$3  # path to genecode gtf file
EDITING_GTF_INTERSECT=$4  # path to intersection of editing DB (only A bases) and gtf files
SNP_GTF_INTERSECT=$5  # path to intersection of SNP and gtf files
EDITING_GTF_BAM_INTERSECT=$6  # path to intersection of editing DB, gtf and bam files
SNP_GTF_BAM_INTERSECT=$7  # path to intersection of SNP, gtf and bam files
export SNAME=$8  # sample name
export N=$9  # number of threads

export OUTPUT_FILE=${OUTPUT_DIR}/2.${SNAME}.gene_filter
LOGFILE=${OUTPUT_DIR}/step2.log

(
    # run htseq to filter non gene reads
    # TODO ask: how to handle if there are REF file already. the OUT_NAMES causes problem
    BAM_DIR=$(dirname $FILTERED_BAM_PATH)
    bamtools split -in ${FILTERED_BAM_PATH} -reference
    OUT_NAMES=$(for f in $BAM_DIR/*REF*.bam; do echo -n "-o ${f}_temp_htseq_out.bam "; done)
    htseq-count -f bam -i gene_name -t gene -m union -s yes -n ${N} ${OUT_NAMES} -c ${OUTPUT_DIR}/2.${SNAME}_gene_counts_htseq.tsv $BAM_DIR/*REF*.bam ${ANNOTATION_GTF}
    samtools merge -@ ${N} -O SAM ${OUTPUT_FILE}.sam $BAM_DIR/*REF*temp_htseq_out.bam
#    rm ${BAM_DIR}/*REF*.bam
#    rm ${OUT_NAMES}

#    htseq-count -f bam -i gene_name -t gene -m union -s yes -n ${N} -o ${OUTPUT_FILE}.sam -c ${OUTPUT_DIR}/2.${SNAME}_gene_counts_htseq.tsv ${FILTERED_BAM_PATH} ${ANNOTATION_GTF}

    # add header to the sam file
#    samtools view -H ${FILTERED_BAM_PATH} | cat - ${OUTPUT_FILE}.sam > ${OUTPUT_FILE}_temp.sam

    # replace reheadered (temp) file with previous file
#    mv ${OUTPUT_FILE}_temp.sam ${OUTPUT_FILE}.sam

    # keep only gene sites and save into a bam file
    grep -v "__" ${OUTPUT_FILE}.sam | samtools sort -@ ${N} - -o ${OUTPUT_FILE}.bam
    samtools index ${OUTPUT_FILE}.bam

    # get statistics on filtered file
    samtools flagstat -@ ${N} ${OUTPUT_FILE}.bam > ${OUTPUT_DIR}/2.${SNAME}_flagstat.txt

    # remove sam files
#    rm ${OUTPUT_FILE}.sam

    #TODO Run without these files
    # RUN INTERSECTIONS WITH EDITING DB ANALYSIS
    if [ ! -z "$EDITING_GTF_INTERSECT" ]
    then
      # make intersection between editing A_I sites-transcriptomes-filtered bam file
      bedtools intersect -u -header -a ${EDITING_GTF_INTERSECT} -b ${OUTPUT_FILE}.bam > ${EDITING_GTF_BAM_INTERSECT}

      # make venn diagram for editing and bam file
      bedtools genomecov -max 1 -split -ibam ${OUTPUT_FILE}.bam > ${OUTPUT_DIR}/2.${SNAME}_genomecov.txt
      bam_positions_count=$(tail -n 1 ${OUTPUT_DIR}/2.${SNAME}_genomecov.txt | awk '{print $3}')
      EDITING_GTF_INTERSECT_count=$(grep -v "#" ${EDITING_GTF_INTERSECT} | cut -f 1,2,3 | sort |  uniq | wc -l)
      EDITING_GTF_INTERSECT_output_bam_intersect=$(grep -v "#" ${EDITING_GTF_BAM_INTERSECT} | cut -f 1,2,3 | sort |  uniq | wc -l)
      python -c "import sys;sys.path.append('/home/labs/bioservices/shared/rarevar/scrarevar');import sc_rna_variants.statistic_plots as sp; sp.plot_venn2_diagram({'10':${EDITING_GTF_INTERSECT_count-EDITING_GTF_INTERSECT_output_bam_intersect},'01':${bam_positions_count-EDITING_GTF_INTERSECT_output_bam_intersect},'11':${EDITING_GTF_INTERSECT_output_bam_intersect}}, ('editing DB positions','aligned positions', ''), '${OUTPUT_DIR}/2.editing_gtf.bam_gtf.venn.png', 'count of positions on gene - ${SNAME}')"
    fi

    # RUN INTERSECTIONS WITH SNP DB ANALYSIS
    if [ ! -z "$SNP_GTF_INTERSECT" ]
    then
      # make intersection between SNP-transcriptomes-filtered bam file
      bedtools intersect -u -header -a ${SNP_GTF_INTERSECT} -b ${OUTPUT_FILE}.bam > ${SNP_GTF_BAM_INTERSECT}

      # make venn diagram for snp and bam file
      SNP_GTF_INTERSECT_count=$(grep -v "#" ${SNP_GTF_INTERSECT} | cut -f 1,2 | sort | uniq | wc -l)
      SNP_GTF_INTERSECT_output_bam_intersect=$(grep -v "#" ${SNP_GTF_BAM_INTERSECT} | cut -f 1,2 | sort | uniq | wc -l)
      python -c "import sys;sys.path.append('/home/labs/bioservices/shared/rarevar/scrarevar');import sc_rna_variants.statistic_plots as sp; sp.plot_venn2_diagram({'10':${SNP_GTF_INTERSECT_count-SNP_GTF_INTERSECT_output_bam_intersect},'01':${bam_positions_count-SNP_GTF_INTERSECT_output_bam_intersect},'11':${SNP_GTF_INTERSECT_output_bam_intersect}}, ('SNP DB positions','aligned positions', ''), '${OUTPUT_DIR}/2.snp_gtf.bam_gtf.venn.png', 'count of positions on gene - ${SNAME}')"
    fi

    echo end_of_log_file 1>&2

) >& $LOGFILE