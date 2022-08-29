import os
import subprocess
import logging

import sc_rna_variants.analysis_utils
import sc_rna_variants.bio_functions

logger = logging.getLogger(__name__)


def run_step1(input_bam, barcodes_cluster_file, min_mapq, cigar_clipping_allowed, max_gene_length, max_no_basecall,
              tag_for_umi, tag_for_cell_barcode, output_folder, threads):
    # create the filtered bam from which the variants will be counted
    filtered_bam_path = sc_rna_variants.bio_functions.create_filtered_bam(input_bam, barcodes_cluster_file,
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


def run_step2(script_path, filtered_bam_path, output_dir, annotation_gtf, editing_gtf_intersect, snp_gtf_intersect,
              editing_gtf_bam_intersect, snp_gtf_bam_intersect, sname, threads):
    step2_string = f"bash {script_path} {filtered_bam_path} {output_dir} {annotation_gtf} {editing_gtf_intersect} {snp_gtf_intersect} {editing_gtf_bam_intersect} {snp_gtf_bam_intersect} {sname} {threads}"
    os.system(step2_string)


def run_step3(input_bam, genome_fasta, tag_for_umi, tag_for_cell_barcode, output_folder, threads):
    logger.info('starting to accumulate reads data')
    sc_rna_variants.bio_functions.variants_finder(input_bam, genome_fasta, tag_for_umi, tag_for_cell_barcode,
                                                  output_folder, threads)


def run_step4(input_dir, output_dir, annotation_gtf, editing_db, snp_db):
    # load the mutated and unmutated data frames
    df_mutated = sc_rna_variants.analysis_utils.load_tables(os.path.join(input_dir, "3.mismatch_dictionary.bed"),
                                                            mutated=True)
    df_unmutated = sc_rna_variants.analysis_utils.load_tables(os.path.join(input_dir, "3.no_mismatch_dictionary.bed"),
                                                              mutated=False)

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
    if editing_db != None:
        sc_rna_variants.analysis_utils.add_intersections_with_SNP_and_edit_DB(output_dir, editing_db, snp_db)


def run_step5(output_dir, SCRvar_aggregated_bed_file, min_cb_per_pos, min_mutation_umis, min_total_umis,
              min_mutation_rate, snp_db, editing_db, sname):
    # get agg position tables
    df_merged_agg, df_merged_agg_filtered = sc_rna_variants.analysis_utils.get_df_and_filtered_df(SCRvar_aggregated_bed_file, min_cb_per_pos,
                                                                   min_mutation_umis, min_total_umis,
                                                                   min_mutation_rate)

    # get open tables and filter them
    df_mut_open, df_unmutated, df_merged_open = sc_rna_variants.analysis_utils.get_open_table(output_dir)
    df_merged_open, df_merged_open_filtered = sc_rna_variants.analysis_utils.combine_data_from_agg_to_open_table(df_merged_open, df_merged_agg,
                                                                                  df_merged_agg_filtered)

    # make plots
    logger.info("started to make plots")
    sc_rna_variants.analysis_utils.get_stat_plots(df_merged_open, df_mut_open, df_unmutated, df_merged_agg, df_merged_open_filtered,
                   df_merged_agg_filtered,
                   output_dir, sname)

    # write statistics to text filescvar
    sc_rna_variants.analysis_utils.write_statistics_numbers(df_merged_open, df_merged_open_filtered, output_dir, min_cb_per_pos, min_mutation_umis,
                             min_total_umis, min_mutation_rate)

    # make intersections with SNP and edit DB
    if editing_db != None:
        logger.info("started to make intersection with Data Bases")
        sc_rna_variants.analysis_utils.run_snp_edit_DB_intersections(df_merged_agg, df_merged_agg_filtered, df_merged_open,
                                      df_merged_open_filtered, output_dir, snp_db, editing_db,
                                      sname)


def run_step6(SCRvar_aggregated_bed_file, output_dir,
              read_per_barcode_raw_bam,
              min_cb_per_pos,
              min_mutation_umis, min_total_umis, min_mutation_rate, reditools_data, annotation_gtf,
              mismatch_dict_bed, barcode_clusters, atacseq_gcoverage_min, atacseq_gfrequency_min,
              sname):
    df, df_filtered = sc_rna_variants.analysis_utils.get_df_and_filtered_df(SCRvar_aggregated_bed_file, min_cb_per_pos, min_mutation_umis,
                                             min_total_umis, min_mutation_rate)
    print("shape of mutation table:", df.shape)
    print("shape of mutation table after filtering:", df_filtered.shape)

    # step 6.1 - add gene names, and fing interactions
    # df_filtered = step6_1_add_gene_name_from_gtf(df_filtered, output_dir, annotation_gtf)  # df_filtered

    # TODO ask what to include if no editing and snp DB are given
    # analysis of editing sites
    if 'is_editing' in df_filtered.columns:
        df_filtered = sc_rna_variants.analysis_utils.drop_editing_and_snp_overlap(df_filtered)
        print("shape of mutation table after drop position with edit and SNP overlaps:", df_filtered.shape)

        df_edit = df_filtered.loc[df_filtered['is_editing'] == 1]
        print("shape of aggregated editing sites table is:", df_edit.shape)

        # get open mismatch table of editing sites
        df_open_mismatches_editing, df_open_mismatches = sc_rna_variants.analysis_utils.get_mismatches_tables(df_edit, mismatch_dict_bed,
                                                                               barcode_clusters)
        print("shape of open mismatch editing sites table is:", df_open_mismatches_editing.shape)

        sc_rna_variants.analysis_utils.exploratory_data_analysis(df_filtered, df_edit, df_open_mismatches_editing, df_open_mismatches,
                                  read_per_barcode_raw_bam, output_dir, sname)

    # Analysis
    # TODO: ask where should the output go
    if barcode_clusters:
        sc_rna_variants.analysis_utils.interactions_anlysis(df_open_mismatches_editing, output_dir, sname)

    # Analysis
    # TODO: ask where should the output go
    if (reditools_data):
        # merge REDItools data to our data
        df_reditools = sc_rna_variants.analysis_utils.get_REDItools_data(reditools_data)
        df_filtered = sc_rna_variants.analysis_utils.merge_REDItools_data(df_filtered, df_reditools)
        # df_filtered = sc_rna_variants.analysis_utils.drop_high_prob_snp(df_filtered, atacseq_gcoverage_min, atacseq_gfrequency_min)
        # print(
        #     f"shape after filtering position with REDItools genomic information: {df_filtered.shape} gCoverage>={atacseq_gcoverage_min}, gFrequency>={atacseq_gfrequency_min}:")
        df_filtered.to_csv(os.path.join(output_dir, "6.aggregated_per_position.filtered.REDItools.bed"), index=False, sep='\t')
