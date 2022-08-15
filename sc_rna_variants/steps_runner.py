import os
import subprocess
import logging

import sc_rna_variants.analysis_utils
import sc_rna_variants.bio_functions

logger = logging.getLogger(__name__)


def run_step1(input_bam, filtered_barcodes_list, min_mapq, cigar_clipping_allowed, max_gene_length, max_no_basecall,
              tag_for_umi, tag_for_cell_barcode, output_folder, threads):
    # create the filtered bam from which the variants will be counted
    filtered_bam_path = sc_rna_variants.bio_functions.create_filtered_bam(input_bam, filtered_barcodes_list,
                                                                          min_mapq, cigar_clipping_allowed,
                                                                          max_gene_length, max_no_basecall,
                                                                          tag_for_umi, tag_for_cell_barcode,
                                                                          output_folder, threads)

    # add chr to chromosome names in bam files
    subprocess.run([
        f"samtools view -H {filtered_bam_path} | sed -e '/SN:chr/!s/SN:\([0-9XY]*\)/SN:chr&/' -e '/SN:chrM/!s/SN:MT/SN:chrM&/' | samtools reheader - {filtered_bam_path} > {filtered_bam_path}_temp"],
        shell=True)
    os.remove(filtered_bam_path)
    os.rename(f"{filtered_bam_path}_temp", filtered_bam_path)
    subprocess.run(['samtools', 'index', filtered_bam_path])
    return filtered_bam_path


def run_step2(script_path, filtered_bam_path, output_dir, annotation_gtf, editing_gtf_intersect,  snp_gtf_intersect, editing_gtf_bam_intersect, snp_gtf_bam_intersect, sname, threads):
    step2_string = f"sh {script_path} {filtered_bam_path} {output_dir} {annotation_gtf} {editing_gtf_intersect} {snp_gtf_intersect} {editing_gtf_bam_intersect} {snp_gtf_bam_intersect} {sname} {threads}"
    os.system(step2_string)


def run_step3(input_bam, genome_fasta, tag_for_umi, tag_for_cell_barcode, output_folder, threads):
    logger.info('starting to accumulate reads data')
    sc_rna_variants.bio_functions.variants_finder(input_bam, genome_fasta, tag_for_umi, tag_for_cell_barcode, output_folder, threads)


def run_step4(input_dir, output_dir, annotation_gtf, editing_db_path, snp_db_path):
    # load the mutated and unmutated data frames
    df_mutated = sc_rna_variants.analysis_utils.load_tables(os.path.join(input_dir, "3.mismatch_dictionary.bed"), mutated=True)
    df_unmutated = sc_rna_variants.analysis_utils.load_tables(os.path.join(input_dir, "3.no_mismatch_dictionary.bed"), mutated=False)

    # create aggregated file of data
    df_mutated_agg = sc_rna_variants.analysis_utils.aggregate_df(df_mutated)

    # merge aggregated mutated and (aggregated) unmutated table
    df_merged_agg = sc_rna_variants.analysis_utils.merge_dfs(df_mutated_agg, df_unmutated)

    # add columns of percentage of mutated UMIs from mutated umis and from total umis
    df_merged_agg = sc_rna_variants.analysis_utils.add_fractions_of_mutated_umis(df_merged_agg)

    # reorder the aggregated file
    df_merged_agg = sc_rna_variants.analysis_utils.reorder_and_sort_agg_df(df_merged_agg)

    # add gene names
    df_merged_agg = sc_rna_variants.analysis_utils.add_gene_name_from_gtf(df_merged_agg, output_dir, annotation_gtf)

    # save the aggregated file
    sc_rna_variants.analysis_utils.save_df(df_merged_agg, output_dir, "4.aggregated_per_position.bed")

    # find intersection between df and databases
    if editing_db_path != None:
        sc_rna_variants.analysis_utils.add_intersections_with_SNP_and_edit_DB(output_dir, snp_db_path, editing_db_path)
