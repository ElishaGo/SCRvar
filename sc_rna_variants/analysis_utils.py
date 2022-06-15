import os
import logging
import pandas as pd

logger = logging.getLogger(__name__)


def agg_unmutated(df):
    """
    Checks if there is a need to aggregate the unmutated dataframe on the base position.
    If so, the basic aggregation function is sum.
    TODO: explain when this can happen
    """
    # check if aggregation is needed
    if df.duplicated('position', keep=False).sum() > 0:
        df = df.groupby('position').agg(
            {'count of unmutated cell barcodes': 'sum', 'unmutated multi reads': 'sum',
             'unmutated single reads': 'sum', 'total unmutated umi count': 'sum'}).reset_index()
    return df


def add_full_position_notation(df):
    """Add one column with full position notation"""
    return df['#chrom'] + ":" + df.chromStart.astype(str) + "-" + df.chromEnd.astype(str) + "," + df.strand


def load_df(df_path):
    df = pd.read_csv(df_path, sep='\t')
    if 'TABLE1_hg38' in df_path.split(os.sep)[-1]:
        df.rename(columns={df.columns[0]: "#chrom", df.columns[1]: 'chromEnd', 'Strand': 'strand'}, inplace=True)
        df = df.fillna('.')
    return df


def load_tables(path, mutated=True):
    """Load and preprocess the data of unmutated cells from raw_stats_unmutated.tsv.
    TODO: Change renames in source funtions
    """
    df = load_df(path)
    df.rename(
        columns={'chromosome': '#chrom', 'start': 'chromStart', 'end': 'chromEnd', 'same multi': 'same multi reads',
                 'transition multi': 'transition multi reads', 'reverse multi': 'reverse multi reads',
                 'transvertion multi': 'transvertion multi reads', 'same single': 'same single reads',
                 'transition single': 'transition single reads', 'reverse single': 'reverse single reads',
                 'transvertion single': 'transvertion single reads',
                 'direction': 'strand', 'unique cells': 'count of unmutated cell barcodes',
                 'multiples': 'unmutated multi reads', 'singles': 'unmutated single reads'
                 }, inplace=True)

    df['position'] = add_full_position_notation(df)

    if mutated:
        mutation_cols = ['same multi reads', 'transition multi reads', 'reverse multi reads',
                         'transvertion multi reads',
                         'same single reads', 'transition single reads', 'reverse single reads',
                         'transvertion single reads']
        df['total umi count'] = df[mutation_cols].sum(axis=1)
    else:
        # df['total unmutated umi count'] = df['unmutated multi reads'] + df['unmutated single reads']
        df = agg_unmutated(df)
    df.sort_values(by=['position'], inplace=True)
    return df


def get_grouped(df):
    """return a groupby object grouped by 'position' columns"""
    return df.groupby('position')


def save_df(df, output_folder, name):
    """save pandas df to output folder"""
    df.to_csv(os.path.join(output_folder, name), index=False, sep='\t')


def merge_dfs(df_mutated, df_unmutated):
    """Merge open mutation table and aggregated unmutetated table.
    NOTICE - the unmutated data is aggregated, so you should only look at one line per each position"""
    # merge aggregated mutations and non mutations tables
    df_m = df_mutated.merge(df_unmutated.drop(['#chrom', 'chromStart', 'chromEnd', 'strand'], axis=1), how='left',
                            on='position')

    # if missing values where added while merging, fill with 0
    na_before_merge = df_mutated.isna().sum().sum() + df_unmutated.isna().sum().sum()
    if (df_m.isna().sum().sum()) != na_before_merge:
        logger.debug('Missing values where found after merging mutated and unmutated files, probably becaues there '
                     'are no matching positions for all the mutations. The missing values are transformed to 0')
        df_m.fillna(0, inplace=True)
    return df_m


def filter_positions(df_agg, min_mutation_cb_to_filter, min_mutation_umis, min_total_umis, min_mutation_rate):
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


def get_df_and_filtered_df(df_path, min_cb_per_pos, min_mutation_umis, min_total_umis, min_mutation_rate):
    """function to load the df and filter it"""
    def get_edit_intersections(df, colname):
        """helper function to keep only editing sites which intersect witht the DB,
        with 'a' base as reference, or with other references but in the '-' strand"""
        temp = df.loc[df_agg_intrsct[colname] > 0, :]
        temp = temp[(temp['reference base'] == 'a') & (temp['strand'] == '+') |
                    ((temp['reference base'] == 't') & (temp['strand'] == '-'))]

        # set all sites to 0, except those we found as editing sites
        df.loc[:, colname] = 0
        df.loc[temp.index, colname] = 1
        return df

    # load df with intersections notations
    df_agg_intrsct = pd.read_csv(df_path, sep='\t')

    # define intersections to be binary (1 - if any overlap with db occured, 0 otherwise)
    df_agg_intrsct.loc[df_agg_intrsct['is_snp'] > 0, 'is_snp'] = 1
    df_agg_intrsct.loc[df_agg_intrsct['is_editing'] > 0, 'is_editing'] = 1

    # # save df
    # df_agg_intrsct.to_csv(df_intersection_path, sep='\t', index=False)

    # get filtered df
    df_agg_intrsct_filtered = filter_positions(df_agg=df_agg_intrsct,
                                               min_mutation_cb_to_filter=min_cb_per_pos,
                                               min_mutation_umis=min_mutation_umis,
                                               min_total_umis=min_total_umis,
                                               min_mutation_rate=min_mutation_rate)
    return df_agg_intrsct, df_agg_intrsct_filtered


def write_statistics_numbers(df_merged, df_merged_filtered, output_folder):
    def write_stats_to_file(file_handler, df):
        file_handler.write("number of unique (position, CB) in table: %s \n" % str(df.shape[0]))
        file_handler.write("number of unique positions: %s \n" % str(df['position'].nunique()))
        file_handler.write("number of unique cell barcodes: %s \n" % str(df['cell barcode'].nunique()))

    with open(os.path.join(output_folder, '5.statistics_log.txt'), 'w') as f:
        f.write("General numbers information of table in open mode:\n\n")
        f.write("Before filtering:\n")
        write_stats_to_file(f, df_merged)
        f.write("After filtering:\n")
        write_stats_to_file(f, df_merged_filtered)

