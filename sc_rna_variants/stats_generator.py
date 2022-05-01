import os
import logging

import numpy as np
import pandas as pd
from pandarallel import pandarallel
import sc_rna_variants.analysis_utils

logger = logging.getLogger(__name__)


def reorder_and_sort_agg_df(df):
    # reorder columns
    cols = ["chrom", "chromStart", "chromEnd", 'position', 'percent of non ref from all cells', 'strand',
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

    df.sort_values(by=['chrom', 'chromStart'], inplace=True)  # sort by position?
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
        {'chrom': 'first', 'chromStart': 'first', 'chromEnd': 'first', 'strand': 'first', 'reference base': 'first',
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


def merge_dfs(df_mutated, df_unmutated):
    """Merge open mutation table and aggregated unmutetated table.
    NOTICE - the unmutated data is aggregated, so you should only look at one line per each position"""
    # merge aggregated mutations and non mutations tables
    df_m = df_mutated.merge(df_unmutated.drop(['chrom', 'chromStart', 'chromEnd', 'strand'], axis=1), how='left',
                            on='position')

    # if missing values where added while merging, fill with 0
    na_before_merge = df_mutated.isna().sum().sum() + df_unmutated.isna().sum().sum()
    if (df_m.isna().sum().sum()) != na_before_merge:
        logger.debug('Missing values where found after merging mutated and unmutated files, probably becaues there '
                     'are no matching positions for all the mutations. The missing values are transformed to 0')
        df_m.fillna(0, inplace=True)
    return df_m


def run(args):
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
    df_merged_agg = merge_dfs(df_mutated_agg, df_unmutated)

    # initialize pandarallel for parallel pandas apply. used in the following function
    pandarallel.initialize(nb_workers=args.threads)
    logger.info("started to count UMIs in aggregated file")
    df_merged_agg = add_counts_of_umis(df_merged_agg)

    # reorder and save the aggregated file
    logger.info("reorder and save file")
    df_merged_agg = reorder_and_sort_agg_df(df_merged_agg)
    sc_rna_variants.analysis_utils.save_df(df_merged_agg, args.output_dir, "aggregated_tsv.tsv")

    # TODO: myabe we don't need this
    # save other tables
    # sc_rna_variants.analysis_utils.save_df(df_merged, args.output_dir, "merged_mutated_unmutated_no_agg.tsv")
