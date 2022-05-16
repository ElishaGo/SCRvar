import os
import logging
from datetime import datetime

import sys # for development environments
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)  # for development environments

from sc_rna_variants.DB_intersections import *
from sc_rna_variants.statistic_plots import *
import sc_rna_variants.utils
import sc_rna_variants.analysis_utils

logging.getLogger('matplotlib').setLevel(logging.CRITICAL)


def print_frequencies(df_merged, df_merged_agg, output_folder):
    with open(os.path.join(output_folder, 'general_numbers.txt'), 'w') as f:
        f.write("General numbers information of table in open mode:\n\n")
        f.write("number of unique (position, CB) in table: %s \n" % str(df_merged.shape[0]))
        f.write("number of unique positions: %s \n" % str(df_merged['position'].nunique()))
        f.write("number of unique cell barcodes: %s \n" % str(df_merged['cell barcode'].nunique()))

        # f.write("\nGeneral numbers information of table in aggregated (by position) mode:\n\n")
        # f.write("number of rows in table: %s \n" %str(df_merged_agg.shape[0]))


def filter_by_cb_count(df_agg, min_mutation_cb_to_filter, min_mutation_umis, min_total_umis, min_mutation_rate):
    """filtering function for aggregated tables and open table by the same positions from aggregated filtered table."""

    #  'true values' - drop positions tations and probably hard to get insights from
    def filter_rare_mut(df, min_mutation_rate):
        df = df[df['percent of non ref from all cells'] > min_mutation_rate]
        return df

    # first condition to filter by
    cond_1 = (df_agg['count of mutated cell barcodes'] >= min_mutation_cb_to_filter)

    # second condition to filter by
    mutation_umi_counts = df_agg['total mutation umi count']
    total_umi_count = mutation_umi_counts + \
                      df_agg['unmutated multi reads'] + \
                      df_agg['unmutated single reads']
    cond_2 = ((mutation_umi_counts >= min_mutation_umis) & (total_umi_count >= min_total_umis))

    # filter the aggregated df
    df_agg_filt = df_agg[cond_1 & cond_2]
    df_agg_filt = filter_rare_mut(df_agg_filt, min_mutation_rate)

    return df_agg_filt


def filter_open_and_agg_tables(df, df_agg, min_mutation_cb_to_filter, min_mutation_umis, min_total_umis, min_mutation_rate):
    # filter aggregated table
    df_agg_filt = filter_by_cb_count(df_agg, min_mutation_cb_to_filter, min_mutation_umis, min_total_umis, min_mutation_rate)

    # filter the open table by the position which were filtered in the aggregated df
    filter_idx = df_agg_filt['position'].values
    df_filt = df[df['position'].isin(filter_idx)]
    return df_filt, df_agg_filt


def get_stat_plots(df_merged_open, df_merged_agg, df_merged_filtered, df_merged_agg_filtered, output_folder, sname):
    """
    Function to manage all plots creation.
    Input - mutation table in open form and aggregated form.
    """
    # plot not grouped data
    plot_cb_occurences_hist(df_merged_open, df_merged_filtered, output_folder, sname)
    plot_umi_per_reference_base(df_merged_open, df_merged_filtered, output_folder, sname)
    plot_heatmap_mutation_per_base(df_merged_open, df_merged_filtered, output_folder, sname)

    # plot data grouped by position
    plot_cb_count_overall(df_merged_agg, df_merged_agg_filtered, output_folder, sname)
    plot_cb_count_per_position(df_merged_agg, df_merged_agg_filtered, output_folder, sname)


def run_snp_edit_DB_intersections(input_dir, output_dir, snp_db_path, editing_db_path, sname, atacseq):
    # find intersection between df and databases
    find_intersections_with_SNP_and_edit_DB(input_dir, output_dir, snp_db_path, editing_db_path)

    # get the df with intersections, before and after filtering
    df_agg_intersect, df_agg_intrsct_filtered = get_df(output_dir)

    # make Venn diagrams of the intersections
    make_venn_diagrams(df_agg_intersect, df_agg_intrsct_filtered, output_dir, snp_db_path, editing_db_path, sname)

    plot_heatmap_mutation_per_base_DB(df_agg_intersect, df_agg_intrsct_filtered, output_dir, sname)

    # if ATACseq data if supplied, remove potential SNP sites
    if (atacseq):
        intersect_with_atacseq(df_agg_intersect, output_dir, atacseq)


def run_step5(args):
    df_merged_agg = sc_rna_variants.analysis_utils.load_df(os.path.join(args.input_dir, "4_test_aggregation_per_position.tsv"))
    # TODO: instead of loading mut_open and unmutated and merge them, look into the plots function and see how we can avoid this
    df_mut_open = sc_rna_variants.analysis_utils.load_tables(
        "/home/eligol/Documents/01_WIS/scrarevar/data/test/step3_mismatch_dictionary/3_mismatch_dictionary.bed6",
        mutated=True)
    df_unmutated = sc_rna_variants.analysis_utils.load_tables(
        "/home/eligol/Documents/01_WIS/scrarevar/data/test/step3_mismatch_dictionary/3_no_mismatch_dictionary.bed6",
        mutated=False)
    # df_mut_open = sc_rna_variants.analysis_utils.load_tables("/home/eligol/Documents/01_WIS/scrarevar/data/outputs/scrarevar_outputs/human_RSS540_RSS451_merged/3_mismatch_dictionary.bed6", mutated=True)
    # df_unmutated = sc_rna_variants.analysis_utils.load_tables("/home/eligol/Documents/01_WIS/scrarevar/data/outputs/scrarevar_outputs/human_RSS540_RSS451_merged/3_no_mismatch_dictionary.bed6", mutated=False)
    df_merged_open = sc_rna_variants.analysis_utils.merge_dfs(df_mut_open, df_unmutated)

    # get filtered data for both aggregated and open aggregated df
    df_merged_filtered, df_merged_agg_filtered = filter_open_and_agg_tables(df_merged_open, df_merged_agg, args.min_cb_per_pos,
                                                                    args.min_mutation_umis, args.min_total_umis,
                                                                            args.min_mutation_rate)
    # make plots
    logger.info("started to make plots")
    get_stat_plots(df_merged_open, df_merged_agg, df_merged_filtered, df_merged_agg_filtered, args.output_dir, args.sname)

    # write data to text file
    logger.info("started to make frequency text file")
    print_frequencies(df_merged_open, df_merged_agg, args.output_dir)

    # make intersections with SNP and edit DB
    logger.info("started to make intersection with SNP and editing data bases")
    run_snp_edit_DB_intersections(args.input_dir, args.output_dir, args.snp_db_path, args.editing_db_path, args.sname, args.atacseq)


##################################################################################################################
def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(formatter_class=sc_rna_variants.utils.ArgparserFormater, description="This script filter positions, and make plots for analysis",)

    # positional arguments
    parser.add_argument('input_dir', type=sc_rna_variants.utils.assert_is_directory,
                        help='folder with raw_stats.tsv and raw_umutated_stats.tsv files from scrnavariants.py')
    parser.add_argument('output_dir', type=sc_rna_variants.utils.assert_is_directory, help='folder for outputs')
    parser.add_argument('snp_db_path', type=sc_rna_variants.utils.assert_is_file, help='path to known SNP sites file')
    parser.add_argument('editing_db_path', type=sc_rna_variants.utils.assert_is_file, help='path to known editing sites file')

    # optional arguments
    parser.add_argument('--min_cb_per_pos', default=5, type=int,
                        help='position with less cell barcodes will be filtered')
    parser.add_argument('--min_mutation_umis', default=10, type=int,
                        help='position with less mutated UMIs will be filtered')
    parser.add_argument('--min_total_umis', default=20, type=int,
                        help='position with less number of mutated + unmutated UMIs will be filtered')
    parser.add_argument('--min_mutation_rate', default=0.1, type=int,
                        help='position with less rate of mutation will be filtered')
    parser.add_argument('--atacseq', type=str, help='path to atacseq file')

    # Meta arguments
    parser.add_argument('--log-file', default=os.path.join(sys.argv[2], "5_filtering_positions_and_snp_editing_DB_intersections.log"),
                        help='a log file for tracking the program\'s progress')
    parser.add_argument('--sname', type=str, help='sample name to add to outputs')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    startTime = datetime.now()
    args = parse_arguments()

    # initialize logger
    sc_rna_variants.config_logging(args.log_file)
    logger = logging.getLogger("positions_filtering_and_plots")
    logger.info('positions_filtering_and_plots started')
    logger.debug('Running with parameters:\n%s' % '\n'.join(
        ['%s: %s' % (key, value) for key, value in vars(args).items()]))

    # run filtering and create plots
    run_step5(args)

    print(datetime.now() - startTime)
    logger.info('Step 5 finished')