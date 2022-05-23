import argparse
import os
import logging
from datetime import datetime
import numpy as np
import pandas as pd
from pandarallel import pandarallel

import sys  # for development environments
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)  # for development environments

import sc_rna_variants.analysis_utils
from sc_rna_variants.utils import assert_is_directory, ArgparserFormater

# logger = logging.getLogger(__name__)
logging.getLogger('matplotlib').setLevel(logging.CRITICAL)


def reorder_and_sort_agg_df(df):
    # reorder columns
    cols = ["#chrom", "chromStart", "chromEnd", 'position', 'percent of non ref from all cells', 'strand',
            'count of unmutated cell barcodes', 'count of mutated cell barcodes',
            'percent of non ref only from mutated cells', 'reference base',
            'same multi reads', 'transition multi reads', 'reverse multi reads', 'transvertion multi reads',
            'same single reads', 'transition single reads', 'reverse single reads', 'transvertion single reads',
            'mixed reads', 'total mutation umi count', 'unmutated single reads',
            'unmutated multi reads', 'bin of 1 UMI per cell - #cell with mut',
            'bin of 1 UMI per cell - #median % non ref umis per barcode',
            'bin of 2 UMI per cell - #cell with mut', 'bin of 2 UMI per cell - #median % non ref umis per barcode'
        , 'bin of 3 UMI per cell - #cell with mut', 'bin of 3 UMI per cell - #median % non ref umis per barcode'
        , 'bin of 4+ UMI per cell - #cell with mut', 'bin of 4+ UMI per cell - #median % non ref umis per barcode'
        , 'aggregated cell barcodes']
    df = df[cols]

    df.sort_values(by=['#chrom', 'chromStart'], inplace=True)  # sort by position?
    return df


def get_non_ref_percent(line, cols_no_unmutated, cols_all_umi_counts):
    """helper function to calculate the fraction of mutated UMIs from UMIs in mutated cells,
    and fraction of mutated UMIs from all UMIs including not mutated cells
    """
    mutation_umi_count = line[
        ~line.index.str.startswith(('same', 'unmutated', 'reference'))].sum()
    umi_count_no_unmutated = line[cols_no_unmutated].sum()
    all_umi_count = line[cols_all_umi_counts].sum()
    percent_of_non_ref_total = round(mutation_umi_count / all_umi_count * 100, ndigits=2)
    percent_of_non_ref_mutated = round(mutation_umi_count / umi_count_no_unmutated * 100, ndigits=2)
    return percent_of_non_ref_mutated, percent_of_non_ref_total


def add_counts_of_umis(df):
    """
    From aggreagted tsv, add a columns with percent of UMIs in cells
    add two columns with counts of UMIs to table
    """
    all_umi_cols = ['same multi reads', 'transition multi reads', 'reverse multi reads', 'transvertion multi reads',
                    'same single reads', 'transition single reads', 'reverse single reads', 'transvertion single reads',
                    'unmutated multi reads', 'unmutated single reads']
    df_umi_cols = df[(['reference base'] + all_umi_cols)]
    cols_no_unmutated = df_umi_cols.columns[
        ~df_umi_cols.columns.str.startswith(('unmutated', 'reference'))]  # not sure if this is faster that way
    cols_all_umi_counts = df_umi_cols.columns[~df_umi_cols.columns.str.startswith('reference')]
    df[['percent of non ref only from mutated cells', 'percent of non ref from all cells']] = \
        df_umi_cols.parallel_apply(get_non_ref_percent, args=(cols_no_unmutated, cols_all_umi_counts),
                                   result_type='expand', axis=1)
    return df


def count_per_position_statistics(df, num_of_umis):
    """Helper function to calculate statistics of the cells withing each position.
    Calculation is made separately for cells with 1,2,3 and 4+ umi counts."""
    if num_of_umis == 4:
        filtered_by_umi_count = df[df['total umi count'] >= num_of_umis]
        num_of_umis = '4+'
    else:
        filtered_by_umi_count = df[df['total umi count'] == num_of_umis]

    # get number of unique cell barcodes within position and with the specific number of UMIs.
    num_cells = filtered_by_umi_count.groupby('position')['cell barcode'].count()

    # calculate the median of the value of fraction of mutated UMIs within position and specific number of UMIs
    if num_of_umis == 1:  # when there is 1 mutated UMI the fraction will always be 100
        fraction_mut_umi_median = pd.Series(np.full(shape=(len(num_cells),), fill_value=100), index=num_cells.index)
    else:
        fraction_mut_umi_median = filtered_by_umi_count.groupby('position')['percent of non ref'].median()

    # return two columns: count of UMIs and median fraction
    return pd.DataFrame([num_cells, fraction_mut_umi_median],
                        index=['bin of {} UMI per cell - #cell with mut'.format(num_of_umis),
                               'bin of {} UMI per cell - #median % non ref umis per barcode'.format(num_of_umis)]).T


def per_position_statistics(df):
    # calculate statistics of UMIs of cells in each position
    umi_stats_1_cell = count_per_position_statistics(df, 1)
    umi_stats_2_cell = count_per_position_statistics(df, 2)
    umi_stats_3_cell = count_per_position_statistics(df, 3)
    umi_stats_4up_cell = count_per_position_statistics(df, 4)

    # merge to one table
    per_position_stats = umi_stats_1_cell.merge(
        umi_stats_2_cell.merge(
            umi_stats_3_cell.merge(
                umi_stats_4up_cell, how='outer', left_index=True, right_index=True),
            how='outer', left_index=True, right_index=True),
        how='outer', left_index=True, right_index=True)

    return per_position_stats


def aggregate_existing_columns(df):
    df_grouped = df.groupby('position')
    df_agg = df_grouped.agg(
        {'#chrom': 'first', 'chromStart': 'first', 'chromEnd': 'first', 'strand': 'first', 'reference base': 'first',
         'same multi reads': 'sum', 'transition multi reads': 'sum', 'reverse multi reads': 'sum',
         'transvertion multi reads': 'sum',
         'same single reads': 'sum', 'transition single reads': 'sum', 'reverse single reads': 'sum',
         'transvertion single reads': 'sum',
         'mixed reads': 'sum', 'total umi count': 'sum', 'cell barcode': lambda x: ','.join(x),
         # 'count of unmutated cell barcodes': 'first',  # keep only first occurrence per position of unmutated data
         # 'unmutated multi reads': 'first', 'unmutated single reads': 'first'
         }
    ).reset_index()

    # rename columns
    df_agg.rename(columns={"cell barcode": "aggregated cell barcodes",
                           'total umi count': 'total mutation umi count'}, inplace=True)
    return df_agg


def add_per_position_stats(df_open, df_agg):
    per_position_stats = per_position_statistics(df_open)

    # merge the UMI statistics with table
    df_agg_merged = df_agg.merge(per_position_stats, left_on='position', right_index=True)
    return df_agg_merged


def aggregate_df(df):
    """
    Aggregate the mutated dataframe on the position column.
    In addition, add count of cells in each position.
    """
    df_agg = aggregate_existing_columns(df)

    # get additional statistics per position
    df_agg = add_per_position_stats(df, df_agg)

    # add column with count of CB in each position . This can maybe be shorter by using groupby and count
    df_agg.loc[:, 'count of mutated cell barcodes'] = df_agg.apply(
        lambda x: len(x['aggregated cell barcodes'].split(',')), axis=1)

    return df_agg


def find_intersections_with_SNP_and_edit_DB(output_dir, snp_db_path, editing_db_path):
    agg_df_path = os.path.join(output_dir, '4_aggregated_per_position.tsv')
    snp_temp_path = os.path.join(output_dir, '4_snp_intersect.tsv')
    df_intersection = os.path.join(output_dir, '4_aggregated_per_position_intersect.tsv')

    # add '#' to header of df_aggregated
    os.system(f"head -c 1 {agg_df_path} | grep -q '#' || sed -i '1s/^/#/' {agg_df_path}")

    # both files must be sorted if you use '-sorted' which reduce memory usage
    # find intersection with snp db
    os.system(f"bedtools intersect -c -header -sorted -a {agg_df_path} -b {snp_db_path} > {snp_temp_path}")

    # add column name 'is_snp'
    os.system(f"sed -i '1 s/.*/&\tis_snp/' {snp_temp_path}")

    # find intersection with editing non rep db
    os.system(f"bedtools intersect -s -c -header -a {snp_temp_path} -b {editing_db_path} > {df_intersection}")

    # add column name 'is_editing_non_rep'
    os.system(f"sed -i '1 s/.*/&\tis_editing/' {df_intersection}")

    # remove temp files
    # os.system(f"rm {snp_temp_path}")


def run_step4(args):

    # load the mutated and unmutated data frames
    logger.info("Loading and preprocessing the data frames")
    df_mutated = sc_rna_variants.analysis_utils.load_tables(os.path.join(args.input_dir, "3_mismatch_dictionary.bed6"),
                                                            mutated=True)
    df_unmutated = sc_rna_variants.analysis_utils.load_tables(
        os.path.join(args.input_dir, "3_no_mismatch_dictionary.bed6"), mutated=False)

    # merge mutated and unmutated files to one file
    # TODO : don't merge the open file. merge only the aggregated table
    # df_merged = merge_dfs(df_mutated, df_unmutated)

    # create aggregated file of data
    logger.info("started to aggregate mismatch_dictionary table")
    df_mutated_agg = aggregate_df(df_mutated)

    # merge aggregated mutated and (aggregated) unmutated table
    logger.info("started to merge the files")
    df_merged_agg = sc_rna_variants.analysis_utils.merge_dfs(df_mutated_agg, df_unmutated)

    # initialize pandarallel for parallel pandas apply. used in the following function
    pandarallel.initialize(nb_workers=args.threads)
    logger.info("started to count UMIs in aggregated file")
    df_merged_agg = add_counts_of_umis(df_merged_agg)

    # reorder and save the aggregated file
    logger.info("reorder and save file")
    df_merged_agg = reorder_and_sort_agg_df(df_merged_agg)
    sc_rna_variants.analysis_utils.save_df(df_merged_agg, args.output_dir,
                                           "4_aggregated_per_position.tsv")

    # find intersection between df and databases
    find_intersections_with_SNP_and_edit_DB(args.output_dir, args.snp_db_path, args.editing_db_path)


def make_output_dir(string):
    assert_is_directory(string)
    output_path = os.path.join(string, "step4_outputs")
    os.makedirs(output_path, exist_ok=True)
    return output_path


def parse_arguments(arguments=None):
    """argument parsing wrapper function
    helper functions and classes are found in sc_rna_variants.utils
    # TODO explain more the 'epilog
    """
    parser = argparse.ArgumentParser(
        formatter_class=ArgparserFormater,
        description="""This script aggregates and learns statistics on the output files from 'step3_mismatch_dictionary.py'.""",
        epilog='''Outputs aggregated tsv, figures and statistics.'''
    )

    # positional arguments
    parser.add_argument('input_dir', type=assert_is_directory,
                        help='folder with raw_stats.tsv and raw_umutated_stats.tsv files from step3_mismatch_dictionary.py')
    parser.add_argument('output_dir', help='folder for step outputs', type=assert_is_directory)
    parser.add_argument('snp_db_path', type=sc_rna_variants.utils.assert_is_file, help='path to known SNP sites file')
    parser.add_argument('editing_db_path', type=sc_rna_variants.utils.assert_is_file,
                        help='path to known editing sites file')

    # optional arguments
    parser.add_argument('--min_cb_per_pos', default=5, type=int,
                        help='position with less cell barcodes will be filtered')
    parser.add_argument('--min_mutation_umis', default=10, type=int,
                        help='position with less mutated UMIs will be filtered')
    parser.add_argument('--min_total_umis', default=20, type=int,
                        help='position with less number of mutated + unmutated UMIs will be filtered')
    parser.add_argument('--min_mutation_rate', default=0.1, type=int,
                        help='position with less rate of mutation will be filtered')

    # Meta arguments
    parser.add_argument('--threads', type=int,
                        help='number of available threads', default=1)
    parser.add_argument('--log-file',
                        default=os.path.join(sys.argv[2], '4_aggregated_per_position_and_statisitcs.log'),
                        help='a log file for tracking the program\'s progress')
    parser.add_argument('--sname', type=str, help='sample name to add to outputs')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    startTime = datetime.now()
    args = parse_arguments()

    # initialize logger
    sc_rna_variants.config_logging(args.log_file)
    logger = logging.getLogger("aggregation_per_position")
    logger.info('Aggregation per position started')
    logger.debug('Running with parameters:\n%s' % '\n'.join(
        ['%s: %s' % (key, value) for key, value in vars(args).items()]))

    # run statistics analysis
    run_step4(args)

    print(datetime.now() - startTime)
    logger.info('Step 4 finished')
    print("Step 4 finished")
