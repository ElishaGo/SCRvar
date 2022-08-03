import os
import logging
import argparse
from datetime import datetime
import matplotlib.pyplot as plt

import sys  # for development environments
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)  # for development environments

from sc_rna_variants.analysis_utils import load_tables, merge_dfs, get_df_and_filtered_df, write_statistics_numbers
# from sc_rna_variants.utils import assert_is_file, assert_is_directory, ArgparserFormater
from sc_rna_variants.statistic_plots import *
import sc_rna_variants.utils

pd.set_option('display.max_columns', None)
logging.getLogger('matplotlib').setLevel(logging.CRITICAL)
logger = logging.getLogger(__name__)


def run_venn(df, df_filtered, column_name, db_total_count, labels, input_dir, sname):
    # DB position set is combination of positions from table, and strings representing non overlaping positions.
    set1 = set(df[df[column_name] != 0].position.to_list() +
               ['not_in_table_position' + str(i) for i in range(db_total_count)])
    set2 = set(df.position)  # aggregated data
    set3 = set(df_filtered.position)  # filtered aggregated data

    plot_venn3_diagram([set1, set2, set3], labels,
                       output_name=os.path.join(input_dir, f'5.venn_diagram_{labels[0]}.png'),
                       title=f'count of positions on gene - {sname}')
    plot_mutated_CB_hist(df, column_name, input_dir, labels[0], sname)


def plot_heatmap_mutation_A_base(df_merged, df_merged_filtered, output_dir, sname):
    bases = ['a', 'c', 'g', 't']
    umi_cols = ['same multi reads', 'transition multi reads', 'reverse multi reads', 'transvertion multi reads',
                'same single reads', 'transition single reads', 'reverse single reads', 'transvertion single reads']

    # create matrix with counts of mutations observed
    count_matrices = make_counts_matrix(df_merged, df_merged_filtered, bases, umi_cols)

    # plot only A base mutations
    plt.clf()
    all_a_mutations_mat = np.zeros((4, 5))
    for i, count_matrix_set in enumerate(count_matrices):
        count_matrix = count_matrix_set[0]
        a_mut_data = count_matrix.sum(axis=0)
        all_a_mutations_mat[i, :] = a_mut_data

    vmin, vmax = np.floor(np.log10(all_a_mutations_mat.min())), np.floor(np.log10(all_a_mutations_mat.max()))
    s = sns.heatmap(np.log10(all_a_mutations_mat), linewidth=0.5, annot=np.array(all_a_mutations_mat),
                    cbar_kws={'label': 'log 10'}, cmap='brg', vmin=vmin, vmax=vmax,
                    xticklabels=['same_all', 'same_mut', 'transition', 'reverse', 'transversion'],
                    yticklabels=['a_single_reads', 'a_multi_reads', 'a_single_reads_filtered',
                                 'a_multi_reads_filtered'])
    s.set_yticklabels(s.get_yticklabels(), rotation=360)
    s.set_xticklabels(s.get_xticklabels(), rotation=30, ha='right')
    plt.title("Counts of mutations 'A' reference base - {}".format(sname))
    plt.ylabel("data type")
    plt.xlabel("mutation")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "5.heatmap_A_mutations.png"), bbox_inches='tight')


def make_venn_diagrams(df_agg_intrsct, df_filtered, output_dir, snp_db_path, editing_db_path, sname):
    """function to get information on intersections of tables with databases"""
    snp_total_count = os.popen("grep -v '#' {} | wc -l".format(snp_db_path)).read()  # count non header lines
    editing_total_count = os.popen("grep -v '#' {} | wc -l".format(editing_db_path)).read()

    # convert to int
    snp_total_count, editing_total_count = int(snp_total_count), int(editing_total_count)

    # make Venn diagrams for snp, editing rep and editing non_rep intersections
    run_venn(df_agg_intrsct, df_filtered, 'is_snp', snp_total_count, ['SNP DB positions', 'aligned positions', 'aligned positions filtered'],
             output_dir, sname)
    run_venn(df_agg_intrsct, df_filtered, 'is_editing', editing_total_count,
             ['editing DB positions', 'aligned positions', 'aligned positions filtered'], output_dir, sname)


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
    plot_cb_occurences_hist(df_merged_open, df_merged_filtered, fig_path=os.path.join(output_folder, "5.cb_distribution.png"), sname=sname, is_snp=False)
    plot_umi_per_reference_base(df_merged_open, df_merged_filtered, output_folder, sname, with_unmut=False, figname="5.umi_per_reference_base")
    plot_umi_per_reference_base(df_merged_open, df_merged_filtered, output_folder, sname, with_unmut=True,
                                figname="5.umi_per_reference_base_with_unmutated")
    plot_heatmap_mutation_per_base(df_merged_open, df_merged_filtered, output_folder, sname)  # use nonmut data
    plot_heatmap_mutation_A_base(df_merged_agg, df_merged_filtered, output_folder, sname)

    # plot data grouped by position
    plot_cb_count_overall(df_merged_agg, df_merged_agg_filtered, output_folder, sname)
    plot_cb_count_per_position(df_merged_agg, df_merged_agg_filtered, output_folder, sname, with_unmut=True)
    plot_cb_count_per_position(df_merged_agg, df_merged_agg_filtered, output_folder, sname, with_unmut=False)


def run_snp_edit_DB_intersections(df_agg_intersect, df_agg_intrsct_filtered, df_merged_open, df_merged_open_filtered, output_folder, snp_db_path, editing_db_path, sname):
    output_folder = os.path.join(output_folder, '5.DB_intersect_effect')
    os.makedirs(output_folder, exist_ok=True)

    # make Venn diagrams of the intersections
    make_venn_diagrams(df_agg_intersect, df_agg_intrsct_filtered, output_folder, snp_db_path, editing_db_path, sname)

    # plot mutations per cell with snp and edit notation
    plot_cb_occurences_hist(df_merged_open, df_merged_open_filtered,
                            fig_path=os.path.join(output_folder, "5.cb_distribution_snp.png"), sname=sname, is_snp=True)

def get_open_table(dir_path):
    dir3_outputs = os.path.join(os.path.dirname(dir_path), 'step3_mismatch_dictionary')
    df_mut_open = load_tables(os.path.join(dir3_outputs, '3.mismatch_dictionary.bed'), mutated=True)
    df_unmutated = load_tables(os.path.join(dir3_outputs, '3.no_mismatch_dictionary.bed'),mutated=False)
    df_merged_open = merge_dfs(df_mut_open, df_unmutated)
    return df_mut_open, df_unmutated, df_merged_open


def run_step5(input_dir, output_dir, min_cb_per_pos, min_mutation_umis, min_total_umis, min_mutation_rate, snp_db_path,
              editing_db_path, sname):
    # get agg position tables
    df_merged_agg, df_merged_agg_filtered = get_df_and_filtered_df(os.path.join(input_dir, '4.aggregated_per_position.bed'), min_cb_per_pos, min_mutation_umis, min_total_umis, min_mutation_rate)

    # get open tables and filter them
    df_mut_open, df_unmutated, df_merged_open = get_open_table(output_dir)
    df_merged_open, df_merged_open_filtered = combine_data_from_agg_to_open_table(df_merged_open, df_merged_agg, df_merged_agg_filtered)

    # make plots
    logger.info("started to make plots")
    get_stat_plots(df_merged_open, df_mut_open, df_unmutated, df_merged_agg, df_merged_open_filtered, df_merged_agg_filtered,
                   output_dir, sname)

    # write statistics to text file
    write_statistics_numbers(df_merged_open, df_merged_open_filtered, output_dir, min_cb_per_pos, min_mutation_umis, min_total_umis, min_mutation_rate)

    # make intersections with SNP and edit DB
    if editing_db_path != None:
        logger.info("started to make intersection with Data Bases")
        run_snp_edit_DB_intersections(df_merged_agg, df_merged_agg_filtered, df_merged_open, df_merged_open_filtered, output_dir, snp_db_path, editing_db_path,
                                      sname)



##################################################################################################################
def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(formatter_class=sc_rna_variants.utils.ArgparserFormater, description="", )

    # positional arguments
    parser.add_argument('input_dir', type=sc_rna_variants.utils.assert_is_directory, help='step 4 output folder')
    parser.add_argument('output_dir', type=sc_rna_variants.utils.assert_is_directory, help='folder for outputs')
    parser.add_argument('editing_db_path', help='path to known editing sites file')
    parser.add_argument('snp_db_path', help='path to known SNP sites file')

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
    parser.add_argument('--log-file',
                        default=os.path.join(sys.argv[2], "5.filtering_positions_and_snp_editing_DB_intersections.log"),
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
    # logger.debug('Running with parameters:\n%s' % '\n'.join(
    #     ['%s: %s' % (key, value) for key, value in vars(args).items()]))

    # run filtering and create plots
    run_step5(args.input_dir, args.output_dir, args.min_cb_per_pos, args.min_mutation_umis, args.min_total_umis,
              args.min_mutation_rate, args.snp_db_path, args.editing_db_path, args.sname)

    print(datetime.now() - startTime)
    logger.info('Step 5 finished')
