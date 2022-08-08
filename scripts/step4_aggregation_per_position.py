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

import sc_rna_variants
from sc_rna_variants.analysis_utils import load_tables, merge_dfs, save_df
from sc_rna_variants.utils import assert_is_directory, assert_is_file, ArgparserFormater

# logger = logging.getLogger(__name__)
logging.getLogger('matplotlib').setLevel(logging.CRITICAL)


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
    df_path = os.path.join(output_dir, 'temp_6.df.bed')
    df.to_csv(df_path, index=False, sep='\t')
    intersections_annotation_gtf = os.path.join(output_dir, "temp_6.genecode_intersect.bed")
    create_mismatches_gtf_intersections(df_path=df_path, path_to_gtf=annotation_gtf, out_fpath=intersections_annotation_gtf)
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

def add_intersections_with_SNP_and_edit_DB(output_dir, snp_db_path, editing_db_path):
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


def run_step4(input_dir, output_dir, annotation_gtf, snp_db_path, editing_db_path):
    # load the mutated and unmutated data frames
    df_mutated = load_tables(os.path.join(input_dir, "3.mismatch_dictionary.bed"), mutated=True)
    df_unmutated = load_tables(os.path.join(input_dir, "3.no_mismatch_dictionary.bed"), mutated=False)

    # create aggregated file of data
    df_mutated_agg = aggregate_df(df_mutated)

    # merge aggregated mutated and (aggregated) unmutated table
    df_merged_agg = merge_dfs(df_mutated_agg, df_unmutated)

    # add columns of percentage of mutated UMIs from mutated umis and from total umis
    df_merged_agg = add_fractions_of_mutated_umis(df_merged_agg)

    # reorder the aggregated file
    df_merged_agg = reorder_and_sort_agg_df(df_merged_agg)

    # add gene names
    df_merged_agg = add_gene_name_from_gtf(df_merged_agg, output_dir, annotation_gtf)

    # save the aggregated file
    save_df(df_merged_agg, output_dir, "4.aggregated_per_position.bed")

    # find intersection between df and databases
    if editing_db_path != None:
        add_intersections_with_SNP_and_edit_DB(output_dir, snp_db_path, editing_db_path)


def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(ArgparserFormater)

    # positional arguments
    parser.add_argument('input_dir', type=assert_is_directory, help='folder with mismatch_dictionaries (step 3)')
    parser.add_argument('output_dir', help='folder for step4 outputs', type=assert_is_directory)
    parser.add_argument('editing_db_path', type=assert_is_file, help='path to known editing sites file')
    parser.add_argument('snp_db_path', type=assert_is_file, help='path to known SNP sites file')
    parser.add_argument('annotation_gtf', type=assert_is_file, help='path to transcriptome gtf file')

    # optional arguments
    parser.add_argument('--min_cb_per_pos', default=5, type=int,
                        help='position with less cell barcodes will be filtered')
    parser.add_argument('--min_mutation_umis', default=10, type=int,
                        help='position with less mutated UMIs will be filtered')
    parser.add_argument('--min_total_umis', default=20, type=int,
                        help='position with less number of mutated + unmutated UMIs will be filtered')
    parser.add_argument('--min_mutation_rate', default=0.1, type=int,
                        help='position with less rate of mutation will be filtered')
    parser.add_argument('--atacseq_path', type=str, help='path to atacseq file')
    parser.add_argument('--atacseq_gcoverage_min', type=int, default=5)
    parser.add_argument('--atacseq_gfrequency_min', type=float, default=0.2)

    # Meta arguments
    parser.add_argument('--log-file',
                        default=os.path.join(sys.argv[2], '4.aggregated_per_position_and_statisitcs.log'),
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
    run_step4(args.input_dir, args.output_dir, args.annotation_gtf, args.snp_db_path, args.editing_db_path)

    print(datetime.now() - startTime)
    logger.info('Step 4 finished')
    print("Step 4 finished")
