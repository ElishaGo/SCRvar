import os
import logging
import argparse
import pandas as pd
from datetime import datetime

import sys  # for development environments
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)  # for development environments

from sc_rna_variants.analysis_utils import save_df, load_tables, get_df_and_filtered_df, get_REDItools_data, drop_high_prob_snp
import sc_rna_variants.statistic_plots
from sc_rna_variants.utils import ArgparserFormater, assert_is_directory, assert_is_file

pd.set_option('display.max_columns', None)


# def create_mismatches_gtf_intersections(df_path, path_to_gtf, out_fpath):
#     """use left outer join to add genes information to editing table"""
#     os.system(f"bedtools intersect -s -loj -a {df_path} -b {path_to_gtf} > {out_fpath}")
#
#
# def add_gene_names(df, genecode_gtf_file):
#     """function to add column of gene names.
#     input:  -aggregated mismatches by posotion dataframe
#             -path to gtf file
#     output: aggregated file with aditional column for gene names.
#     Notice that the gtf file contains multiple lines per each position. However, the gene name
#     should be consistent in all the position duplicate lines. Therfore we drop duplicates in some stage here.
#     Positions which didn't apper in the gtf file are removed
#     """
#
#     # load df with gene names
#     gene_names = pd.read_csv(genecode_gtf_file, header=None, sep='\t')
#     # extract the gene names
#     # TODO find more robust way to extract positon and gene names. check fot gtf file convention
#     gene_names = gene_names.loc[:, [3, gene_names.shape[1] - 1]]  # get 'postion' and last columnt
#     gene_names.columns = ['position', 'gene_name']
#
#     # parse the gene name
#     gene_names['gene_name'] = gene_names['gene_name'].map(lambda x: x[x.find("gene_name") + 11:])
#     gene_names['gene_name'] = gene_names['gene_name'].map(lambda x: x[:x.find("\"")])
#
#     ## NOTE: In early stage of the analysis we check intersection between the bam file and the sampe .gtf to get only reads from genes, by using htseq.
#     ## Thus, we would expexct all the positions in our table to be on genes.
#     ## However, when we look for intersections between the mutation table and the gtf file, we get some empty records returned.
#     ## The reason is the htseq looks on reads, and some reads overlap the gene area and the none gene area.
#     ## However now, we look on a position, which is not on the gene, even though the read it came from had some overlap with a gene
#     ## We drop the empty records.
#     gene_names['gene_name'].replace('', 'None', inplace=True)
#     gene_names = gene_names[gene_names['gene_name'] != 'None']
#
#     # drop duplicates
#     gene_names = gene_names.drop_duplicates()
#
#     # merge df with gene names
#     df = df.merge(gene_names, on='position', how='inner')
#
#     return df


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


# def step6_1_add_gene_name_from_gtf(df, output_dir, gtf_path):
#     # add gene names
#     df_path = os.path.join(output_dir, 'temp_6.df.bed')
#     df.to_csv(df_path, index=False, sep='\t')
#     intersections_gtf_path = os.path.join(output_dir, "temp_6.genecode_intersect.bed")
#     create_mismatches_gtf_intersections(df_path=df_path, path_to_gtf=gtf_path, out_fpath=intersections_gtf_path)
#     df = add_gene_names(df, intersections_gtf_path)
#     os.remove(df_path)
#     os.remove(intersections_gtf_path)
#
#     df.to_csv(os.path.join(output_dir, "6.1.aggregated_with_gene_name.bed"), index=False, sep='\t')
#     return df


def add_clusters(df_open, clusters_path):
    """helper function to add clusters notations to open table"""
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

    pt_names = ['clusters_VS_genes_umis', 'clusters_VS_genes_unmutated_umis',
                'cells_VS_genes_umis', 'cells_VS_genes_unmutated_umis'
                # 'cells_VS_position_mutated_umis',
                # 'cells_VS_position_unmutated_umis',
                ]

    pts = [clusters_VS_genes_umis_pt, clusters_VS_genes_unmutated_umis_pt,
           cells_VS_genes_umis_pt, cells_VS_genes_unmutated_umis_pt,
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
        print("name")
        print(pt.shape)
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


def run_step6(input_dir, output_dir, read_per_barcode_raw_bam, min_cb_per_pos, min_mutation_umis, min_total_umis,
              min_mutation_rate, reditools_data, gtf_path, mismatch_dict_bed, barcode_clusters, gcoverage_min,
              gfrequency_min, sname):
    stats_agg_path = os.path.join(input_dir, "4.aggregated_per_position.bed")
    df, df_filtered = get_df_and_filtered_df(stats_agg_path, min_cb_per_pos, min_mutation_umis, min_total_umis, min_mutation_rate)
    print("shape of mutation table:", df.shape)
    print("shape of mutation table after filtering:", df_filtered.shape)

    # step 6.1 - add gene names, and fing interactions
    # df_filtered = step6_1_add_gene_name_from_gtf(df_filtered, output_dir, gtf_path)  # df_filtered

    # TODO ask what to include if no editing and snp DB are given
    # analysis of editing sites
    if 'is_editing' in df_filtered.columns:
        df_filtered = drop_editing_and_snp_overlap(df_filtered)
        print("shape of mutation table after drop position with edit and SNP overlaps:", df_filtered.shape)

        df_edit = df_filtered.loc[df_filtered['is_editing'] == 1]
        print("shape of aggregated editing sites table is:", df_edit.shape)

        # get open mismatch table of editing sites
        df_open_mismatches_editing, df_open_mismatches = get_mismatches_tables(df_edit, mismatch_dict_bed, barcode_clusters)
        print("shape of open mismatch editing sites table is:", df_open_mismatches_editing.shape)

        exploratory_data_analysis(df_filtered, df_edit, df_open_mismatches_editing, df_open_mismatches,
                                  read_per_barcode_raw_bam, output_dir, sname)

    # Analysis
    # TODO: ask where should the output go
    if barcode_clusters:
        interactions_anlysis(df_open_mismatches_editing, output_dir, sname)

    # Analysis
    # TODO: ask where should the output go
    if (reditools_data):
        # merge REDItools data to our data
        df_reditools = get_REDItools_data(reditools_data)
        df_filtered = merge_REDItools_data(df_filtered, df_reditools)
        df_filtered = drop_high_prob_snp(df_filtered, gcoverage_min, gfrequency_min)
        df_filtered.to_csv(os.path.join(output_dir, "6_ANALYSIS_with_REDItools_data.bed"), index=False, sep='\t')
        print(
            f"shape after filtering position with REDItools genomic information: {df_filtered.shape} gCoverage>={gcoverage_min}, gFrequency>={gfrequency_min}:")


##################################################################################################################
def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(formatter_class=ArgparserFormater, description="", )

    # positional arguments
    parser.add_argument('input_dir', type=assert_is_directory, help='step 4 output folder')
    parser.add_argument('output_dir', type=assert_is_directory, help='folder for outputs')
    parser.add_argument('mismatch_dict_bed', type=assert_is_file, help='path to 3.mismatch_dictionary.bed')
    parser.add_argument('read_per_barcode_raw_bam', type=assert_is_file, help='count of reads per CB in raw bam file')
    parser.add_argument('gtf_path', type=assert_is_file, help='path to gtf file')

    # optional arguments
    parser.add_argument('--barcode_clusters', type=assert_is_file, help='barcodes and clusters analysed by Seurat')
    parser.add_argument('--min_cb_per_pos', default=5, type=int,
                        help='position with less cell barcodes will be filtered')
    parser.add_argument('--min_mutation_umis', default=10, type=int,
                        help='position with less mutated UMIs will be filtered')
    parser.add_argument('--min_total_umis', default=20, type=int,
                        help='position with less number of mutated + unmutated UMIs will be filtered')
    parser.add_argument('--min_mutation_rate', default=0.1, type=int,
                        help='position with less rate of mutation will be filtered')
    parser.add_argument('--reditools_data', type=str, help='path to atacseq file')
    parser.add_argument('--atacseq_gcoverage_min', type=int, default=5)
    parser.add_argument('--atacseq_gfrequency_min', type=float, default=0.2)

    # Meta arguments
    parser.add_argument('--log-file',
                        default=os.path.join(sys.argv[2], "6.filtering_positions_and_snp_editing_DB_intersections.log"),
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

    # run step
    run_step6(args.input_dir, args.output_dir, args.read_per_barcode_raw_bam, args.min_cb_per_pos,
              args.min_mutation_umis, args.min_total_umis, args.min_mutation_rate, args.reditools_data, args.gtf_path,
              args.mismatch_dict_bed, args.barcode_clusters, args.atacseq_gcoverage_min, args.atacseq_gfrequency_min,
              args.sname)

    print(datetime.now() - startTime)
    logger.info('Step 6 finished')
