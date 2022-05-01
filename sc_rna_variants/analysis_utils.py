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
    return df.chrom + ":" + df.chromStart.astype(str) + "-" + df.chromEnd.astype(str) + "," + df.strand


def load_tables(path, mutated=True):
    """Load and preprocess the data of unmutated cells from raw_stats_unmutated.tsv.
    TODO: Change renames in source funtions
    """
    df = pd.read_csv(path, sep='\t')
    df.rename(
        columns={'chromosome': 'chrom', 'start': 'chromStart', 'end': 'chromEnd', 'same multi': 'same multi reads',
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
    df_m = df_mutated.merge(df_unmutated.drop(['chrom', 'chromStart', 'chromEnd', 'strand'], axis=1), how='left',
                            on='position')

    # if missing values where added while merging, fill with 0
    na_before_merge = df_mutated.isna().sum().sum() + df_unmutated.isna().sum().sum()
    if (df_m.isna().sum().sum()) != na_before_merge:
        logger.debug('Missing values where found after merging mutated and unmutated files, probably becaues there '
                     'are no matching positions for all the mutations. The missing values are transformed to 0')
        df_m.fillna(0, inplace=True)
    return df_m