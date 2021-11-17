import os
import logging

import numpy as np
import pandas as pd
from pandarallel import pandarallel

import sc_rna_variants.statistic_plots as statistic_plots
# from statistic_plots import *

logger = logging.getLogger(__name__)


def order_and_save_agg(df, out_folder):
    # reorder columns
    cols = ['chromosome', 'start', 'end', 'position', 'percent of non ref from all cells', 'strand',
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
    # sort and save df
    df.sort_values(by=['chromosome', 'start'], inplace=True)
    save_file(df, out_folder, "aggregated_tsv.tsv")
    # df.to_csv(os.path.join(out_folder, "aggregated_tsv.tsv"), index=False, sep='\t')


def add_counts_of_umis(df):
    """
    From aggreagted tsv, add a columns with percent of UMIs in cells
    """
    def get_non_ref_percent(line,cols_no_unmutated, cols_all_umi_counts):
        """helper function to calculate the fraction of mutated UMIs from UMIs in mutated cells,
        and fraction of mutated UMIs from all UMIs including not mutated cells
        # TODO: check if it faster to use explicit column names instead of 'startswith'"""

        mutation_umi_count = line[
            ~line.index.str.startswith(('same', 'unmutated', 'reference'))].sum()
        # mutation_umi_count = line[  #old version 6.9.21
        #     ~line.index.str.startswith(('R->%s' % line['reference base'].upper(), 'unmutated', 'reference'))].sum()
        umi_count_no_unmutated = line[cols_no_unmutated].sum()
        all_umi_count = line[cols_all_umi_counts].sum()
        percent_of_non_ref_total = round(mutation_umi_count / all_umi_count * 100, ndigits=2)
        percent_of_non_ref_mutated = round(mutation_umi_count / umi_count_no_unmutated * 100, ndigits=2)
        return percent_of_non_ref_mutated, percent_of_non_ref_total

    all_umi_cols = ['same multi reads', 'transition multi reads', 'reverse multi reads', 'transvertion multi reads',
                    'same single reads', 'transition single reads', 'reverse single reads', 'transvertion single reads',
                    'unmutated multi reads', 'unmutated single reads']
    # all_umi_cols = [col for col in df.columns if  # old version 9.6.21
    #                 col.startswith(('R->', 'unmutated multi reads', 'unmutated single reads'))]
    df_umi_cols = df[(['reference base'] + all_umi_cols)]
    cols_no_unmutated = df_umi_cols.columns[~df_umi_cols.columns.str.startswith(('unmutated', 'reference'))]  # not sure if this is faster that way
    cols_all_umi_counts = df_umi_cols.columns[~df_umi_cols.columns.str.startswith('reference')]
    df[['percent of non ref only from mutated cells', 'percent of non ref from all cells']] = \
        df_umi_cols.parallel_apply(get_non_ref_percent, args=(cols_no_unmutated, cols_all_umi_counts), result_type='expand', axis=1)
    return df


def print_frequencies(df_merged, df_merged_agg, output_folder):
    with open(os.path.join(output_folder,'general_numbers.txt'), 'w') as f:
        f.write("General numbers information of table in open mode:\n\n")
        f.write("number of rows in table: %s \n" %str(df_merged.shape[0]))
        f.write("number of unique positions: %s \n" %str(df_merged['position'].nunique()))
        f.write("number of unique cell barcodes: %s \n" %str(df_merged['cell barcode'].nunique()))

        f.write("\nGeneral numbers information of table in aggregated (by position) mode:\n\n")
        f.write("number of rows in table: %s \n" %str(df_merged_agg.shape[0]))


def filter_by_cb_count(df, df_agg, min_mutation_cb_to_filter,min_mutation_umis, min_total_umis):
    """filtering function for aggregated tables and open table by the same positions from aggregated filtered table."""

    # first condition to filter by
    cond_1 = (df_agg['count of mutated cell barcodes'] >= min_mutation_cb_to_filter)

    # second condition to filter by
    mutation_umi_counts = df_agg['total mutation umi count']
    total_umi_count = mutation_umi_counts + \
                      df_agg['unmutated multi reads'] +\
                      df_agg['unmutated single reads']
    cond_2 = ((mutation_umi_counts >= min_mutation_umis) & (total_umi_count >= min_total_umis))

    # filter the aggregated df
    df_agg_filt = df_agg[cond_1 & cond_2]

    # filter the open table by the position which were filtered in the aggregated df
    filter_idx = df_agg_filt['position'].values
    df_filt = df[df['position'].isin(filter_idx)]  # TODO: check if you can use directly the filter idx without 'isin'
    return df_filt, df_agg_filt


def get_stat_plots(df_merged_open, df_merged_agg, args):
    """
    Function to manage all plots creation.
    Input - mutation table in open form and aggregated form.
    """
    # get filtered data for both aggregated and open aggregated df
    df_merged_filtered, df_merged_agg_filtered = filter_by_cb_count(df_merged_open, df_merged_agg, args.min_cb_per_pos,
                                                                    args.min_mutation_umis, args.min_total_umis)

    # plot not grouped data
    statistic_plots.plot_cb_occurences_hist(df_merged_open, df_merged_filtered, args.output_folder, args.sname)
    statistic_plots.plot_umi_per_reference_base(df_merged_open, df_merged_filtered, args.output_folder, args.sname)
    statistic_plots.plot_heatmap_mutation_per_base(df_merged_open, df_merged_filtered, args.output_folder, args.sname)

    # plot data grouped by position
    statistic_plots.plot_cb_count_overall(df_merged_agg, df_merged_agg_filtered, args.output_folder, args.sname)
    statistic_plots.plot_cb_count_per_position(df_merged_agg, df_merged_agg_filtered, args.output_folder, args.sname)


def agg_dfs(df):
    """
    Aggregate the mutated dataframe on the positions.
    In addition, add count of cells in each position.
    """
    def get_umis_stats(num_of_umis):
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

        # return df with two columns: count of UMIs and median fraction
        return pd.DataFrame([num_cells, fraction_mut_umi_median],
                            index=['bin of {} UMI per cell - #cell with mut'.format(num_of_umis),
                                   'bin of {} UMI per cell - #median % non ref umis per barcode'.format(num_of_umis)]).T

                            # index = ['total umi counts {} cells'.format(num_of_umis),
                            #          'median percent of non ref umi {} cells'.format(num_of_umis)]).T

    # TODO: check if we can fasten aggregation by replacing the 'first' with some other function
    df_grouped = df.groupby('position')
    df_agg = df_grouped.agg(
        {'chromosome': 'first', 'start': 'first', 'end': 'first', 'strand': 'first', 'reference base': 'first',
         'same multi reads': 'sum', 'transition multi reads': 'sum', 'reverse multi reads': 'sum', 'transvertion multi reads': 'sum',
         'same single reads': 'sum', 'transition single reads': 'sum', 'reverse single reads': 'sum', 'transvertion single reads': 'sum',
         'mixed reads': 'sum', 'total umi count': 'sum', 'cell barcode': lambda x: ','.join(x),
         'count of unmutated cell barcodes': 'first',  # keep only first occurrence per position of unmutated data
         'unmutated multi reads': 'first', 'unmutated single reads': 'first'}
    ).reset_index()

    # old version 9.6.21
    # df_agg = df_grouped.agg(
    #     {'chromosome': 'first', 'start': 'first', 'end': 'first', 'strand': 'first', 'reference base': 'first',
    #     'R->A multi reads': 'sum', 'R->T multi reads': 'sum', 'R->G multi reads': 'sum',
    #      'R->C multi reads': 'sum', 'R->A single reads': 'sum', 'R->T single reads': 'sum',
    #      'R->G single reads': 'sum', 'R->C single reads': 'sum', 'mixed reads': 'sum','total umi count':'sum',
    #      'cell barcode': lambda x: ','.join(x),
    #      'count of unmutated cell barcodes': 'first',  #  keep only first occurrence per position of unmutated
    #      'unmutated multi reads': 'first',
    #          'unmutated single reads': 'first'}).reset_index()

    # rename columns
    df_agg.rename(columns={"cell barcode": "aggregated cell barcodes",
                           'total umi count': 'total mutation umi count'}, inplace=True)

    # add column with count of CB in each position . This can maybe be shorter by using groupby and count
    df_agg.loc[:, 'count of mutated cell barcodes'] = df_agg.apply(
        lambda x: len(x['aggregated cell barcodes'].split(',')), axis=1)

    # calculate statistics of UMIs of cells in each position
    umi_stats_1_cell = get_umis_stats(1)
    umi_stats_2_cell = get_umis_stats(2)
    umi_stats_3_cell = get_umis_stats(3)
    umi_stats_4up_cell = get_umis_stats(4)

    umi_stats = umi_stats_1_cell.merge(  # merge the umi statistics to one table
        umi_stats_2_cell.merge(
            umi_stats_3_cell.merge(
                umi_stats_4up_cell, how='outer', left_index=True, right_index=True),
            how='outer', left_index=True, right_index=True),
        how='outer', left_index=True, right_index=True)

    # merge the UMI statistics with table
    df_agg = df_agg.merge(umi_stats, left_on='position', right_index=True)
    return df_agg


def merge_dfs(df_mutated, df_unmutated):
    """Merge open mutation table and aggregated unmutetated table.
    NOTICE - the unmutated data is aggregated, so you should only look at one line per each position"""
    # merge aggregated mutations and non mutations tables
    df_m = df_mutated.merge(
        df_unmutated[['position', 'count of unmutated cell barcodes',
                      'unmutated single reads', 'unmutated multi reads']], how='left', on='position')

    # check if missing values where added in merged file. If so, transform to 0
    na_before_merge = df_mutated.isna().sum().sum() + df_unmutated.isna().sum().sum()
    if (df_m.isna().sum().sum()) != na_before_merge:
        logger.debug('Missing values where found after merging mutated and unmutated files, probably becaues there '
                     'are no matching positions for all the mutations. The missing values are transformed to 0')
        df_m.fillna(0, inplace=True)
    return df_m


def agg_unmutated(df):
    """
    Checks if there is a need to aggregate the unmutated dataframe on the base position.
    If so, the basic aggregation function is sum.
    TODO: explain when this can happen
    """
    #check if aggregation is needed
    if df.duplicated('position', keep = False).sum() > 0:
        df = df.groupby('position').agg(
            {'count of unmutated cell barcodes': 'sum', 'unmutated multi reads': 'sum',
             'unmutated single reads': 'sum', 'total umi count': 'sum'}).reset_index()
    return df


def add_position(df):
    """
    Add one column with the position notation.
    TODO:
     -change function name
     -change function to be in one line without calling pnadas eac time
     - add documentation including final format string
    """
    # concat 3 first columns into one column of position
    df['position'] = df.chromosome.str.cat(df['start'].astype(str), sep=':')  #change to chromoseme
    df['position'] = df.position.str.cat(df['end'].astype(str), sep='-')
    df['position'] = df.position.str.cat(df['strand'].astype(str), sep = ',')


def load_mutated(path):
    """
    Load and preprocess the data of unmutated cells from raw_stats_unmutated.tsv.
    TODO: Change renames in source funtions
    """
    df = pd.read_csv(path, sep='\t')
    df.rename(columns={'same multi': 'same multi reads',
                         'transition multi': 'transition multi reads', 'reverse multi': 'reverse multi reads',
                         'transvertion multi': 'transvertion multi reads', 'same single': 'same single reads',
                         'transition single': 'transition single reads','reverse single': 'reverse single reads',
                         'transvertion single': 'transvertion single reads'}, inplace=True)

    add_position(df)  # add column of full coordinates
    mutation_cols = ['same multi reads','transition multi reads', 'reverse multi reads', 'transvertion multi reads',
                     'same single reads', 'transition single reads', 'reverse single reads','transvertion single reads']
    df['total umi count'] = df[mutation_cols].sum(axis=1)
    df.sort_values(by=['position'], inplace=True)
    return df


def load_unmutated(path):
    """
    Load and preprocess the data of unmutated cells from raw_stats_unmutated.tsv.
    TODO: Change renames in source funtions
    """
    df = pd.read_csv(path, sep='\t')
    df.rename(
        columns={'direction': 'strand', 'unique cells': 'count of unmutated cell barcodes',
                 'multiples': 'unmutated multi reads', 'singles': 'unmutated single reads'}, inplace=True)
    df['total unmutated umi count'] = df['unmutated multi reads'] + df['unmutated single reads']
    add_position(df)  # add column of full coordination
    df = agg_unmutated(df)
    df.sort_values(by=['position'], inplace=True)
    return df


def get_grouped(df):
    """return a groupby object grouped by 'position' columns"""
    return df.groupby('position')


def save_file(df, output_folder, name):
    """save pandas df to output folder"""
    df.to_csv(os.path.join(output_folder, name), index=False, sep='\t')


def run(args):
    # TODO: replace the 'arguments' variable with explicit arguments.

    #load the mutated and unmutated data frames
    logger.info("Loading and preprocessing the data frames")
    df_mutated = load_mutated(args.input_mutated_file)
    df_unmutated = load_unmutated(args.input_unmutated_file)

    # merge mutated and unmutated files to one file
    logger.info("started to merge the files")
    df_merged = merge_dfs(df_mutated, df_unmutated)

    # create aggregated file of data
    logger.info("started to aggregate the data")
    df_merged_agg = agg_dfs(df_merged)

    # make plots
    # logger.info("started to make plots")
    get_stat_plots(df_merged, df_merged_agg, args)

    # write data to text file
    logger.info("started to make frequency text file")
    print_frequencies(df_merged, df_merged_agg, args.output_folder)

    # add two columns with counts of UMIs to table
    logger.info("started to count UMIs in aggregated file")
    # initialize pandarallel for parallel pandas apply. used in the following function
    pandarallel.initialize(nb_workers=args.threads)
    df_merged_agg = add_counts_of_umis(df_merged_agg)

    # reorder and save the aggregated file
    logger.info("reorder and save file")
    order_and_save_agg(df_merged_agg, args.output_folder)

    # save other tables
    save_file(df_merged, args.output_folder, "merged_mutated_unmutated_no_agg.tsv")

    print("Code finished")
