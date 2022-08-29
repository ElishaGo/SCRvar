import os
import logging

import pandas as pd
import numpy as np
from pathlib import Path

logger = logging.getLogger(__name__)

import sc_rna_variants.statistic_plots
import sc_rna_variants.utils

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
             'unmutated single reads': 'sum'}).reset_index()
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
    """Load and preprocess the data of unmutated cells from raw_stats_unmutated.tsv."""
    df = load_df(path)

    # add column with full coordination
    df['position'] = add_full_position_notation(df)

    if mutated:
        mutation_cols = ['same multi reads', 'transition multi reads', 'reverse multi reads',
                         'transvertion multi reads',
                         'same single reads', 'transition single reads', 'reverse single reads',
                         'transvertion single reads']
        df['total umi count'] = df[mutation_cols].sum(axis=1)
    else:
        df = agg_unmutated(df)

    df.sort_values(by=['#chrom', 'chromStart'], inplace=True)
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
    logger.info("started to merge the files")
    # merge aggregated mutations and non mutations tables
    df_m = df_mutated.merge(df_unmutated.drop(['#chrom', 'chromStart', 'chromEnd', 'strand'], axis=1), how='left',
                            on='position')

    # if missing values where added while merging, fill with 0
    na_before_merge = df_mutated.isna().sum().sum() + df_unmutated.isna().sum().sum()
    if (df_m.isna().sum().sum()) != na_before_merge:
        logger.debug('Missing values where found after merging mutated and unmutated files, probably becaues there '
                     'are no matching positions for all the mutations. The missing values are transformed to 0')
        df_m[['count of unmutated cell barcodes', 'unmutated multi reads', 'unmutated single reads']] = df_m[
            ['count of unmutated cell barcodes', 'unmutated multi reads', 'unmutated single reads']].fillna(0)
    return df_m


def filter_positions(df_agg, min_mutated_cells_per_position, min_mutated_umis_per_position, min_total_umis_per_position, min_mutation_rate):
    """filtering function for aggregated tables and open table by the same positions from aggregated filtered table."""
    # first condition to filter by
    cond_1 = (df_agg['count of mutated cell barcodes'] >= min_mutated_cells_per_position)

    # second condition to filter by
    mutation_umi_counts = df_agg['total mutation umi count']
    total_umi_count = mutation_umi_counts + \
                      df_agg['unmutated multi reads'] + \
                      df_agg['unmutated single reads']
    cond_2 = ((mutation_umi_counts >= min_mutated_umis_per_position) & (total_umi_count >= min_total_umis_per_position))

    # filter the aggregated df
    df_agg_filt = df_agg[cond_1 & cond_2].copy()

    #  'true values' - drop positions tations and probably hard to get insights from
    df_agg_filt = df_agg_filt[df_agg_filt['percent of non ref from all cells'] > min_mutation_rate]

    return df_agg_filt


def get_df_and_filtered_df(df_path, min_mutated_cells_per_position, min_mutation_umis, min_total_umis, min_mutation_rate):
    """function to load the df and filter it"""
    # load df with intersections notations
    df = pd.read_csv(df_path, sep='\t')

    # get filtered df
    df_filtered = filter_positions(df_agg=df,
                                   min_mutated_cells_per_position=min_mutated_cells_per_position,
                                   min_mutated_umis_per_position=min_mutation_umis,
                                   min_total_umis_per_position=min_total_umis,
                                   min_mutation_rate=min_mutation_rate)
    return df, df_filtered


def write_statistics_numbers(df_merged, df_merged_filtered, output_folder, min_mutated_cells_per_position, min_mutated_umis_per_position,
                             min_total_umis_per_position, min_mutation_rate):
    def write_stats_to_file(file_handler, df):
        file_handler.write("number of unique (position, CB) in table: %s \n" % str(df.shape[0]))
        file_handler.write("number of unique positions: %s \n" % str(df['position'].nunique()))
        file_handler.write("number of unique cell barcodes: %s \n" % str(df['cell barcode'].nunique()))

    with open(os.path.join(output_folder, '5.filtering_statistics.txt'), 'w') as f:
        f.write("General numbers information of table in open mode:\n\n")
        f.write("Filtering parameters:\n")
        f.write("min_mutated_cells_per_position: {}\n".format(min_mutated_cells_per_position))
        f.write("min_mutated_umis_per_position: {}\n".format(min_mutated_umis_per_position))
        f.write("min_total_umis_per_position: {}\n".format(min_total_umis_per_position))
        f.write("min_mutation_rate_per_umi: {}\n".format(min_mutation_rate))

        f.write("\nBefore filtering:\n")
        write_stats_to_file(f, df_merged)
        f.write("\nAfter filtering:\n")
        write_stats_to_file(f, df_merged_filtered)


def get_REDItools_data(ATACseq_path):
    # ON LOCAL SYSTEM
    # columns = ['Region', 'Position', 'Strand', 'gCoverage-q20', 'gFrequency']
    # df_atacseq = pd.read_csv(ATACseq_path, sep='\t', usecols=columns, memory_map=True)
    df_atacseq = pd.read_csv(ATACseq_path, sep='\t', memory_map=True)
    df_atacseq.columns = df_atacseq.columns.map(lambda x: "REDItools_" + x)

    # replace missing values with 0
    df_atacseq['REDItools_gCoverage-q20'] = df_atacseq['REDItools_gCoverage-q20'].replace('-', 0).astype('int32')
    df_atacseq['REDItools_gFrequency'] = df_atacseq['REDItools_gFrequency'].replace('-', 0).astype(
        'float32')  # .convert_dtypes()

    # transform column types to reduce memory
    df_atacseq['REDItools_Position'] = df_atacseq['REDItools_Position'].astype('int32')
    df_atacseq['REDItools_Strand'] = df_atacseq['REDItools_Strand'].astype('int8')
    df_atacseq['REDItools_gCoverage-q20'] = df_atacseq['REDItools_gCoverage-q20'].astype('int32')
    df_atacseq['REDItools_MeanQ'] = df_atacseq['REDItools_MeanQ'].astype('float32')
    df_atacseq['REDItools_Frequency'] = df_atacseq['REDItools_Frequency'].astype('float32')

    # rename column names to be the same as the statistics tables
    df_atacseq = df_atacseq.rename(
        columns={'REDItools_Region': '#chrom', 'REDItools_Position': 'chromStart',
                 'REDItools_Strand': 'REDItools_Strand (0:-, 1:+, 2:unknown)'})
    return df_atacseq


def drop_high_prob_snp(df, gcov_min=5, gfreq_min=0.2):
    """remove positions where there is high probabilty to be SNP"""
    idx_to_remove = (df[(df['REDItools_gCoverage-q20'] >= gcov_min)]['REDItools_gFrequency'] >= gfreq_min).index
    return df.drop(idx_to_remove)


# functions for step4
def create_mismatches_gtf_intersections(df_path, path_to_gtf, out_fpath):
    """creates a temporary table with each intersection in A including information from B"""
    os.system(f"bedtools intersect -s -wo -sorted -a {df_path} -b {path_to_gtf} > {out_fpath}")


def add_gene_names(df, genecode_gtf_file):
    """function to add column of gene names.
    input:  -aggregated mismatches by posotion dataframe
            -path to gtf file
    output: aggregated file with aditional column for gene names.
    Notice that the gtf file contains multiple lines per each position. However, the gene name
    should be consistent in all the position duplicate lines. Therfore we drop duplicates in some stage here.
    Positions which didn't apper in the gtf file are removed
    """

    # load only positions column and column with gene name
    # TODO find more robust way to extract positon and gene names.last column of gtf includes all text info including gene name
    gene_names = pd.read_csv(genecode_gtf_file, header=None, sep='\t').iloc[:, [3, -2]]
    gene_names.columns = ['position', 'gene_name']

    # parse the gene name
    gene_names['gene_name'] = gene_names['gene_name'].map(lambda x: x[x.find("gene_name") + 11:])
    gene_names['gene_name'] = gene_names['gene_name'].map(lambda x: x[:x.find("\"")])

    ## NOTE: In early stage of the analysis we check intersection between the bam file and the sampe .gtf to get only reads from genes, by using htseq.
    ## Thus, we would expexct all the positions in our table to be on genes.
    ## However, when we look for intersections between the mutation table and the gtf file, we get some empty records returned.
    ## The reason is the htseq looks on reads, and some reads overlap the gene area and the none gene area.
    ## However now, we look on a position, which is not on the gene, even though the read it came from had some overlap with a gene
    ## We drop the empty records.
    gene_names = gene_names[gene_names['gene_name'] != '']

    # drop duplicates
    gene_names = gene_names.drop_duplicates()

    # merge df with gene names
    return df.merge(gene_names, on='position', how='inner')


def add_gene_name_from_gtf(df, output_dir, annotation_gtf):
    # add gene names
    df_path = os.path.join(output_dir, 'temp_4.df.bed')
    df.to_csv(df_path, index=False, sep='\t')
    intersections_annotation_gtf = os.path.join(output_dir, "temp_4.genecode_intersect.bed")
    create_mismatches_gtf_intersections(df_path=df_path, path_to_gtf=annotation_gtf,
                                        out_fpath=intersections_annotation_gtf)
    df = add_gene_names(df, intersections_annotation_gtf)
    os.remove(df_path)
    os.remove(intersections_annotation_gtf)
    return df


def reorder_and_sort_agg_df(df):
    # reorder columns
    logger.info("reorder and save file")
    cols = ["#chrom", "chromStart", "chromEnd", 'position', 'percent of non ref from all cells', 'strand',
            'count of unmutated cell barcodes', 'count of mutated cell barcodes',
            'percent of non ref only from mutated cells', 'reference base',
            'same multi reads', 'transition multi reads', 'reverse multi reads', 'transvertion multi reads',
            'same single reads', 'transition single reads', 'reverse single reads', 'transvertion single reads',
            'mixed reads', 'total mutation umi count', 'unmutated single reads',
            'unmutated multi reads', 'bin of 1 UMI per cell - #cell with mut',
            'bin of 1 UMI per cell - #median % non ref umis per barcode',
            'bin of 2 UMI per cell - #cell with mut', 'bin of 2 UMI per cell - #median % non ref umis per barcode',
            'bin of 3 UMI per cell - #cell with mut', 'bin of 3 UMI per cell - #median % non ref umis per barcode',
            'bin of 4+ UMI per cell - #cell with mut', 'bin of 4+ UMI per cell - #median % non ref umis per barcode',
            'aggregated cell barcodes']
    df = df[cols]

    df.sort_values(by=['#chrom', 'chromStart'], inplace=True)
    return df


def add_fractions_of_mutated_umis(df):
    """
    calculate the fraction of mutated UMIs from UMIs in mutated cells,
    and the fraction of mutated UMIs from all UMIs (including not mutated cells)
    """
    logger.info("Adding fractions of mutated UMIs")
    all_umi_cols = ['same multi reads', 'transition multi reads', 'reverse multi reads', 'transvertion multi reads',
                    'same single reads', 'transition single reads', 'reverse single reads', 'transvertion single reads',
                    'unmutated multi reads', 'unmutated single reads']
    cols_mutated_umis = [x for x in all_umi_cols if not x.startswith(('same', 'unmutated'))]
    cols_no_unmutated_umis = [x for x in all_umi_cols if not x.startswith('unmutated')]

    mutation_umi_count = df.loc[:, cols_mutated_umis].sum(axis=1)
    no_muatation_umi_count = df.loc[:, cols_no_unmutated_umis].sum(axis=1)
    all_umi_count = df.loc[:, all_umi_cols].sum(axis=1)
    percent_of_non_ref_mutated = round(mutation_umi_count / no_muatation_umi_count * 100, ndigits=2)
    percent_of_non_ref_total = round(mutation_umi_count / all_umi_count * 100, ndigits=2)
    df['percent of non ref only from mutated cells'] = percent_of_non_ref_mutated
    df['percent of non ref from all cells'] = percent_of_non_ref_total

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
    df_agg = df.groupby('position').agg(
        {'#chrom': 'first', 'chromStart': 'first', 'chromEnd': 'first', 'strand': 'first', 'reference base': 'first',
         'same multi reads': 'sum', 'transition multi reads': 'sum', 'reverse multi reads': 'sum',
         'transvertion multi reads': 'sum', 'same single reads': 'sum', 'transition single reads': 'sum',
         'reverse single reads': 'sum', 'transvertion single reads': 'sum', 'mixed reads': 'sum',
         'total umi count': 'sum', 'cell barcode': lambda x: ','.join(x)
         }).reset_index()

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
    logger.info("started to aggregate mismatch_dictionary table")
    df_agg = aggregate_existing_columns(df)

    # get additional statistics per position
    df_agg = add_per_position_stats(df, df_agg)

    # add column with count of CB in each position . This can maybe be shorter by using groupby and count
    df_agg.loc[:, 'count of mutated cell barcodes'] = df_agg.apply(
        lambda x: len(x['aggregated cell barcodes'].split(',')), axis=1)

    return df_agg


# def add_DB_intersection_column(DB_path, mismatch_table_path, out_path, col_name, is_vcf):
#     """find intersection with snp db"""
#     # create temp copy file
#     temp_out_path = mismatch_table_path[:mismatch_table_path.rfind('.')] + '_temp' + mismatch_table_path[mismatch_table_path.rfind('.'):]
#     os.system(f"cp {mismatch_table_path} {temp_out_path}")
#
#     # add intersections column
#     if is_vcf:
#         # both files must be sorted if you use '-sorted' which reduce memory usage
#         os.system(f"bedtools intersect -c -header -sorted -a {mismatch_table_path} -b {DB_path} > {temp_out_path}")
#     else:
#         os.system(f"bedtools intersect -s -c -header -a {mismatch_table_path} -b {DB_path} > {temp_out_path}")
#
#     # add column name 'is_snp'
#     os.system(f"sed -i '1 s/.*/&\t{col_name}/' {temp_out_path}")
#
#     # replace new file and old file
#     os.system(f"mv {temp_out_path} {mismatch_table_path}")
#
# def add_intersections_with_SNP_and_edit_DB(output_dir, snp_db_path, editing_db_path):
#     agg_mismatch_table_path = os.path.join(output_dir, '4.aggregated_per_position.bed')
#     out_path = os.path.join(output_dir, '4.aggregated_per_position_intersect.bed')
#     if snp_db_path:
#         add_DB_intersection_column(snp_db_path, agg_mismatch_table_path, out_path, 'is_snp', is_vcf=True)
#     if editing_db_path:
#         add_DB_intersection_column(editing_db_path, agg_mismatch_table_path, out_path, 'is_editing', is_vcf=False)

def add_intersections_with_SNP_and_edit_DB(output_dir, editing_db_path, snp_db_path):
    """add column for intersection with SNP and editing DB.
    Note, we use here intersect with -c flag which add column of counts, instead of -u flag which only returns the
    intersected entries."""
    agg_df_path = os.path.join(output_dir, '4.aggregated_per_position.bed')
    snp_temp_path = os.path.join(output_dir, 'temp_4.snp_intersect.bed')
    df_intersection = os.path.join(output_dir, '4.aggregated_per_position_intersect.bed')

    # both files must be sorted if you use '-sorted' which reduce memory usage
    # find intersection with snp db
    os.system(f"bedtools intersect -c -header -sorted -a {agg_df_path} -b {snp_db_path} > {snp_temp_path}")

    # add column name 'is_snp'
    # os.system(f"sed -i '1 s/.*/&\tis_snp/' {snp_temp_path}")

    # find intersection with editing non rep db
    os.system(f"bedtools intersect -s -c -header -a {snp_temp_path} -b {editing_db_path} > {df_intersection}")

    # add column name 'is_editing_non_rep'
    # os.system(f"sed -i '1 s/.*/&\tis_editing/' {df_intersection}")

    # remove temp files
    os.system(f"rm {snp_temp_path}")
    os.system(f"mv {df_intersection} {agg_df_path}")

    # define intersections to be binary (1 - if any overlap with db occured, 0 otherwise)
    df = pd.read_csv(agg_df_path, sep='\t')
    # add columns names to new columns
    cols = df.columns.to_list() + ['is_snp', 'is_editing']
    df = df.reset_index()
    df.columns = cols
    df.loc[df['is_snp'] > 0, 'is_snp'] = 1
    df.loc[df['is_editing'] > 0, 'is_editing'] = 1
    df.to_csv(agg_df_path, index=False, sep='\t')


def run_venn(df, df_filtered, column_name, db_total_count, labels, out_dir, sname):
    # DB position set is combination of positions from table, and strings representing non overlaping positions.
    set1 = set(df[df[column_name] != 0].position.to_list() +
               ['not_in_table_position' + str(i) for i in range(db_total_count)])
    set2 = set(df.position)  # aggregated data
    set3 = set(df_filtered.position)  # filtered aggregated data

    sc_rna_variants.statistic_plots.plot_venn3_diagram([set1, set2, set3], labels,
                       output_name=os.path.join(out_dir, f'5.venn_diagram_{labels[0]}.png'),
                       title=f'count of positions on gene - {sname}')
    sc_rna_variants.statistic_plots.plot_mutated_CB_hist(df, column_name, out_dir, labels[0], sname)


def make_venn_diagrams(df_agg_intrsct, df_filtered, output_dir, snp_db_path, editing_db_path, sname):
    """function to get information on intersections of tables with databases"""
    snp_total_count = os.popen("grep -v '#' {} | wc -l".format(snp_db_path)).read()  # count non header lines
    editing_total_count = os.popen("grep -v '#' {} | wc -l".format(editing_db_path)).read()

    # convert to int
    snp_total_count, editing_total_count = int(snp_total_count), int(editing_total_count)

    # make Venn diagrams for snp, editing rep and editing non_rep intersections
    run_venn(df_agg_intrsct, df_filtered, 'is_snp', snp_total_count, ['covered SNP DB positions', 'mismatch positions', 'filtered mismatch positions'],
             output_dir, sname)
    run_venn(df_agg_intrsct, df_filtered, 'is_editing', editing_total_count,
             ['covered editing DB positions', 'mismatch positions', 'filtered mismatch positions'], output_dir, sname)


def combine_data_from_agg_to_open_table(df_merged_open, df_merged_agg, df_merged_agg_filtered):
    # TODO check if df are changed inside function. if so, you dont need to return
    def add_snp_and_editig_notations(df, df_agg_intersect):
        df['is_snp'] = 0
        df.loc[df['position'].isin(df_agg_intersect.loc[df_agg_intersect['is_snp'] == 1, 'position']), 'is_snp'] = 1
        df['is_editing'] = 0
        df.loc[df['position'].isin(
        df_agg_intersect.loc[df_agg_intersect['is_editing'] == 1, 'position']), 'is_editing'] = 1
        return df

    # add snp and editing notations
    df_merged_open = add_snp_and_editig_notations(df_merged_open, df_merged_agg)

    # filter the open table by the position which were filtered in the aggregated df
    filter_idx = df_merged_agg_filtered['position'].values
    df_filt = df_merged_open[df_merged_open['position'].isin(filter_idx)]
    return df_merged_open, df_filt


def get_stat_plots(df_merged_open, df_mut_open, df_unmutated, df_merged_agg, df_merged_filtered, df_merged_agg_filtered,
                   output_folder, sname):
    """
    Function to manage all plots creation.
    Input - mutation table in open form and aggregated form.
    """
    output_folder = os.path.join(output_folder, '5.filtering_effect')
    os.makedirs(output_folder, exist_ok=True)

    # plot not grouped data
    sc_rna_variants.statistic_plots.plot_cb_occurences_hist(df_merged_open, df_merged_filtered, fig_path=os.path.join(output_folder, "5.cb_distribution.png"), sname=sname, is_snp=False)
    sc_rna_variants.statistic_plots.plot_umi_per_reference_base(df_merged_open, df_merged_filtered, output_folder, sname, with_unmut=False, figname="5.umi_per_reference_base")
    sc_rna_variants.statistic_plots.plot_umi_per_reference_base(df_merged_open, df_merged_filtered, output_folder, sname, with_unmut=True, figname="5.umi_per_reference_base_with_unmutated")
    sc_rna_variants.statistic_plots.plot_heatmap_mutation_per_base(df_merged_open, df_merged_filtered, output_folder, sname)  # use nonmut data
    sc_rna_variants.statistic_plots.plot_heatmap_mutation_a_base(df_merged_agg, df_merged_filtered, output_folder, sname)

    # plot data grouped by position
    sc_rna_variants.statistic_plots.plot_cb_count_overall(df_merged_agg, df_merged_agg_filtered, output_folder, sname)
    sc_rna_variants.statistic_plots.plot_cb_count_per_position(df_merged_agg, df_merged_agg_filtered, output_folder, sname, with_unmut=True)
    sc_rna_variants.statistic_plots.plot_cb_count_per_position(df_merged_agg, df_merged_agg_filtered, output_folder, sname, with_unmut=False)


def run_snp_edit_DB_intersections(df_agg_intersect, df_agg_intrsct_filtered, df_merged_open, df_merged_open_filtered, output_folder, snp_db_path, editing_db_path, sname):
    output_folder = os.path.join(output_folder, '5.DB_intersect_effect')
    os.makedirs(output_folder, exist_ok=True)

    # make Venn diagrams of the intersections
    make_venn_diagrams(df_agg_intersect, df_agg_intrsct_filtered, output_folder, snp_db_path, editing_db_path, sname)

    # plot mutations per cell with snp and edit notation
    sc_rna_variants.statistic_plots.plot_cb_occurences_hist(df_merged_open, df_merged_open_filtered,
                            fig_path=os.path.join(output_folder, "5.cb_distribution_snp.png"), sname=sname, is_snp=True)

def get_open_table(dir_path):
    dir3_outputs = os.path.join(str(Path(dir_path).parent), 'step3_mismatch_dictionary')
    df_mut_open = load_tables(os.path.join(dir3_outputs, '3.mismatch_dictionary.bed'), mutated=True)
    df_unmutated = load_tables(os.path.join(dir3_outputs, '3.no_mismatch_dictionary.bed'),mutated=False)
    df_merged_open = merge_dfs(df_mut_open, df_unmutated)
    return df_mut_open, df_unmutated, df_merged_open


def filter_rare_mut(df, min_mutation_rate):
    """drop positions with rare mutations and probably hard to get insights from"""
    try:
        df = df[df['percent of non ref from all cells'] > min_mutation_rate]
    except:
        df = df[df['percent of non ref from all cells'] > min_mutation_rate]
    print("shape after filtering rare mutations (%non ref from all cells>{}):".format(min_mutation_rate), df.shape)
    return df


def get_filtered_mutation_bias(df, min_mutation_cb, min_mutation_umis, min_total_umis, min_mutation_rate):
    """filtering function for aggregated tables"""
    mutation_umi_counts = df['total mutation umi count']
    total_umi_count = mutation_umi_counts + df['unmutated multi reads'] + df['unmutated single reads']

    cond_1 = (df['count of mutated cell barcodes'] >= min_mutation_cb)
    cond_2 = ((mutation_umi_counts >= min_mutation_umis) & (total_umi_count >= min_total_umis))

    df = df[cond_1 & cond_2]
    print("shape after filtering: mutated CB>={}, (mutation UMIs>={} & min total UMIs>={}):".format(min_mutation_cb,
                                                                                                    min_mutation_umis,
                                                                                                    min_total_umis),
          df.shape)
    df = filter_rare_mut(df, min_mutation_rate)

    return df


def get_filtered(df, min_cb_to_filter, min_total_umis):
    """filtering function for aggregated tables"""
    cond_1 = (df['count of mutated cell barcodes'] + df['count of unmutated cell barcodes'] >= min_cb_to_filter)

    total_umi_count = df['total mutation umi count'] + df['unmutated multi reads'] + df['unmutated single reads']
    cond_2 = (total_umi_count >= min_total_umis)
    return df[cond_1 & cond_2]


def drop_ifs(df):
    """remove positions which appear in both editing and snp sites"""
    idx = df[((df['is_editing_non_rep'] == 1) | (df['is_editing_rep'] == 1)) & (df['is_snp'] == 1)].index
    print("number of positin with snp and editing overlap to remove: %d." % len(idx))
    return df.drop(idx)


def drop_editing_and_snp_overlap(df):
    """remove positions which appear in both editing and snp sites"""
    idx_to_drop = df[(df['is_editing'] >= 1) & (df['is_snp'] >= 1)].index
    return df.drop(idx_to_drop)


# def step6_1_add_gene_name_from_gtf(df, output_dir, annotation_gtf):
#     # add gene names
#     df_path = os.path.join(output_dir, 'temp_6.df.bed')
#     df.to_csv(df_path, index=False, sep='\t')
#     intersections_annotation_gtf = os.path.join(output_dir, "temp_6.genecode_intersect.bed")
#     create_mismatches_gtf_intersections(df_path=df_path, path_to_gtf=annotation_gtf, out_fpath=intersections_annotation_gtf)
#     df = add_gene_names(df, intersections_annotation_gtf)
#     os.remove(df_path)
#     os.remove(intersections_annotation_gtf)
#
#     df.to_csv(os.path.join(output_dir, "6.1.aggregated_with_gene_name.bed"), index=False, sep='\t')
#     return df


def add_clusters(df_open, clusters_path):
    """helper function to add clusters notations to open table"""
    # TODO: dont use pandas to read the file. check how the filter barcodes file is read
    cb_clusters = pd.read_csv(clusters_path, sep='\t', names=['cell barcode', 'cluster'])
    cb_clusters['cluster'] = cb_clusters.apply(lambda x: 'c' + str(x['cluster']), axis=1)
    cb_clusters['cluster cell barcode'] = cb_clusters.apply(lambda x: x['cluster'] + '_' + x['cell barcode'], axis=1)

    # add clusters to cells in open table
    return df_open.merge(cb_clusters.loc[:, ['cell barcode', 'cluster', 'cluster cell barcode']], on='cell barcode',
                         how='left')


def load_and_process_mismatch_table(df_edit, mismatches_path, barcodes_clusters_path):
    # load open tables and add 'is_edit column'
    df_mismatches = load_tables(mismatches_path, mutated=True)

    # add 'is_edit' notation to open table
    df_mismatches['is_edit'] = 0
    df_mismatches.loc[df_mismatches['position'].isin(df_edit['position']), 'is_edit'] = 1

    if barcodes_clusters_path:
        df_mismatches = add_clusters(df_mismatches, barcodes_clusters_path)

    return df_mismatches


def get_mismatches_tables(df_edit, mismatch_dict_bed, barcode_clusters):
    """function to get open table of editing sites, and add a column with count of umis per cell"""
    # create open df with only editing sites
    df_open_mismatches = load_and_process_mismatch_table(df_edit, mismatch_dict_bed, barcode_clusters)
    df_open_mismatches_editing = df_open_mismatches[df_open_mismatches['is_edit'] == 1].copy()

    # create one colum with counts of mutated umis
    mut_umi_cols = ['transition multi reads', 'reverse multi reads', 'transvertion multi reads',
                    'transition single reads', 'reverse single reads', 'transvertion single reads']
    df_open_mismatches_editing.loc[:, "mutated umis per cell"] = df_open_mismatches_editing.loc[:, mut_umi_cols].sum(
        axis=1)

    # create one colum with counts of unmutated umis
    unmut_umi_cols = ['same multi reads', 'same single reads']
    df_open_mismatches_editing.loc[:, "unmutated umis per cell"] = df_open_mismatches_editing.loc[:,
                                                                   unmut_umi_cols].sum(axis=1)

    # add gene names to open_edit_df
    if 'gene_name' in df_edit.columns:
        df_open_mismatches_editing = df_open_mismatches_editing.merge(df_edit.loc[:, ['position', 'gene_name']],
                                                                      on='position', how='left')

    return df_open_mismatches_editing, df_open_mismatches


def exploratory_data_analysis(df, df_edit, df_open_edit, df_mismatches, reads_per_barcode_path, output_dir, sname):
    output_dir = os.path.join(output_dir, '6.plots')
    os.makedirs(output_dir, exist_ok=True)

    sc_rna_variants.statistic_plots.non_ref_from_all_cells(df, output_dir)
    sc_rna_variants.statistic_plots.snp_observed_against_mut_UMI_in_position(df, output_dir, sname)
    # sc_rna_variants.statistic_plots.plot_venn2_diagram()

    # editing EDA
    umis_per_cb_editing = df_open_edit.groupby('cell barcode')["mutated umis per cell"].sum()
    sc_rna_variants.statistic_plots.editing_sites_per_chr(df_edit, output_dir, sname)
    sc_rna_variants.statistic_plots.mutated_umis_per_cell(umis_per_cb_editing, output_dir, sname)

    if reads_per_barcode_path:
        more_analysis(df_mismatches, umis_per_cb_editing, reads_per_barcode_path, min_num_of_cells=2, output_dir=output_dir)


def get_pivot_tables(df):
    # create heatmap for Seurat clusters VS. Genes
    clusters_VS_genes_umis_pt = df.groupby(['cluster', 'gene_name']).agg(
        {'mutated umis per cell': 'sum'}).reset_index().pivot(index='cluster', columns='gene_name',
                                                              values="mutated umis per cell").fillna(0)
    clusters_VS_genes_umis_pt = filter_columns_with_less_than_value(clusters_VS_genes_umis_pt, min_counts_in_col=2)

    # create heatmap for Seurat clusters VS. Genes - unmutated umis
    clusters_VS_genes_unmutated_umis_pt = df.groupby(['cluster', 'gene_name', 'position']).agg(
        {'unmutated umis per cell': 'first'}).groupby(['cluster', 'gene_name']).agg(
        {'unmutated umis per cell': 'sum'}).reset_index().pivot(index='cluster', columns='gene_name',
                                                                values="unmutated umis per cell").fillna(0)
    clusters_VS_genes_umis_pt = filter_columns_with_less_than_value(clusters_VS_genes_umis_pt, min_counts_in_col=2)

    # create heatmap for Cells VS. Genes
    cells_VS_genes_umis_pt = df.groupby(['gene_name', 'cluster cell barcode']).agg(
        {'mutated umis per cell': 'sum'}).reset_index().pivot(index='cluster cell barcode', columns='gene_name',
                                                              values="mutated umis per cell").fillna(0)
    cells_VS_genes_umis_pt = get_repetetive(df=cells_VS_genes_umis_pt, min_genes_to_occure_in=2)

    # create heatmap for Cells VS. Genes- unmutated umis
    cells_VS_genes_unmutated_umis_pt = df.groupby(['gene_name', 'cluster cell barcode']).agg(
        {'unmutated umis per cell': 'first'}).reset_index().pivot(index='cluster cell barcode', columns='gene_name',
                                                                  values="unmutated umis per cell").fillna(0)
    cells_VS_genes_unmutated_umis_pt = get_repetetive(df=cells_VS_genes_unmutated_umis_pt, min_genes_to_occure_in=2)

    # # create heatmap for Cells VS. Positions
    # cells_VS_position_mutated_umis_pt = df.pivot_table(index='cluster cell barcode', columns='position',
    #                                                    values="mutated umis per cell",
    #                                                    aggfunc='sum')

    # # create heatmap for Cells VS. Positions - unmutated umis
    # cells_VS_position_unmutated_umis_pt = df.pivot_table(index='cluster cell barcode',
    #                                                      columns='position',
    #                                                      values="unmutated umis per cell",
    #                                                      aggfunc='sum')

    # create a dictionary between chromosomes and genes
    chr_genes_pairs = pd.Series(df['#chrom'].values, index=df.gene_name).to_dict()

    pt_names = ['clusters_VS_genes_umis', 'cells_VS_genes_umis',
                'clusters_VS_genes_unmutated_umis', 'cells_VS_genes_unmutated_umis'
                # 'cells_VS_position_mutated_umis',
                # 'cells_VS_position_unmutated_umis',
                ]

    pts = [clusters_VS_genes_umis_pt, cells_VS_genes_umis_pt,
           clusters_VS_genes_unmutated_umis_pt, cells_VS_genes_unmutated_umis_pt,
           # cells_VS_position_mutated_umis_pt,
           # cells_VS_position_unmutated_umis_pt,
           chr_genes_pairs]
    return pts, pt_names


def get_repetetive(df, min_genes_to_occure_in):
    """function to filter the matrix, such that only mutated cells (rows) which occure in more than X genes are kept.\
    X is a given paramter"""
    # convert all UMI counts to 1, so we count each cell/gene one time
    # get indices of cells which occures im more than X genes
    df_temp = df.copy()
    df_temp[df_temp > 1] = 1
    cells_idx = df_temp[df_temp.sum(axis=1) > min_genes_to_occure_in].index
    df_repetetive_cells = df.loc[cells_idx, :]
    return df_repetetive_cells


def filter_rows_with_less_than_value(df, min_counts_in_row):
    df_temp = df[df.sum(axis=1) > min_counts_in_row]
    print("shape of df is:", df_temp.shape)
    return df_temp


def filter_columns_with_less_than_value(df, min_counts_in_col):
    df_temp = df.loc[:, df.sum(axis=0) > min_counts_in_col]
    print("shape of df is:", df_temp.shape)
    return df_temp


def plot_clustering_heatmaps(pivot_tables, pt_names, sname, output_dir):
    chr_genes_pairs = pivot_tables.pop(-1)

    for pt, name in zip(pivot_tables, pt_names):
        print("shape of {} is: {}".format(name, pt.shape))
        sc_rna_variants.statistic_plots.make_clusters_heatmap(pt, name=name + str(sname),
                                                              output_dir=output_dir,
                                                              chr_genes_pairs=chr_genes_pairs)


def interactions_anlysis(df_open, output_dir, sname):
    output_dir = os.path.join(output_dir, '6.ANALYSIS_interactions_analysis')
    os.makedirs(output_dir, exist_ok=True)

    pts, pt_names = get_pivot_tables(df_open)

    sc_rna_variants.statistic_plots.plot_umis_per_gene(pts[1], output_dir)
    plot_clustering_heatmaps(pts, pt_names, sname, output_dir)


def get_reads_per_cell(rpc_fpath):
    """function to read file with counts of reads per cell barcode after parsing. Return dataframe"""
    with open(rpc_fpath, "r") as f:
        reads_per_cb = f.readlines()

    df_rpb = {rpc.split()[-1].strip('\n'): rpc.split()[0] for rpc in reads_per_cb}
    df_rpb = pd.DataFrame.from_dict(df_rpb, orient='index',
                                    columns=['number of reads per cell']).reset_index().rename(
        columns={'index': 'cell barcode'})
    df_rpb['cell barcode'] = df_rpb['cell barcode'].map(lambda x: x + '-1')
    df_rpb['number of reads per cell'] = df_rpb['number of reads per cell'].astype('int')
    return df_rpb


def more_analysis(df_mismatches, umis_per_cb_editing, reads_per_barcode_path, min_num_of_cells, output_dir):
    """
    the plot needs to be made from the count of all the umis in cells. this information exist in the bam file before any filtering.
    in addition, add the column of 0 editing sites (cell barcodes not in the editing table (with the 191 rows)).
    we want to see if there is correlation between the number of mutated umis per cell and the number of editing sites in the same cell
    If we will see cells which doesn't follow the correlation, will probably have more editing events than usual
    :param df_mismatches:
    :param reads_per_barcode_path:
    :param min_num_of_cells:
    :return:
    """
    # get number of editing sites per cell
    editings_per_cb_stats = df_mismatches.iloc[:, :].groupby('cell barcode')['is_edit'].sum()
    editingsites_per_cb_stats = editings_per_cb_stats[editings_per_cb_stats >= min_num_of_cells]

    # load counts of reads per cell barcode
    df_rpb = get_reads_per_cell(reads_per_barcode_path)
    df_scatter = pd.merge(editingsites_per_cb_stats, df_rpb, on='cell barcode')
    sc_rna_variants.statistic_plots.editing_events_vs_number_of_reads(df_scatter, output_dir)
    sc_rna_variants.statistic_plots.editing_events_vs_number_of_mutated_umis_per_cell(editingsites_per_cb_stats,
                                                                                      umis_per_cb_editing, output_dir)


def merge_REDItools_data(df, df_reditools):
    """left merge the atacseq into the statistics table"""
    df_merged = pd.merge(df, df_reditools, on=['#chrom', 'chromStart'], how='left')

    # replace missing values with 0
    df_merged['REDItools_gCoverage-q20'] = df_merged['REDItools_gCoverage-q20'].replace('-', 0).fillna(0)
    df_merged['REDItools_gFrequency'] = df_merged['REDItools_gFrequency'].replace('-', 0).fillna(0)
    return df_merged