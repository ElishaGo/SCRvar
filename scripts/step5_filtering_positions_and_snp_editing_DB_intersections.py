import os
import logging
import argparse
from datetime import datetime

import sys # for development environments
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)  # for development environments

from sc_rna_variants.statistic_plots import *
import sc_rna_variants.utils
import sc_rna_variants.analysis_utils
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles
from sc_rna_variants.statistic_plots import get_min_max, make_mut_counts_heatmap


pd.set_option('display.max_columns', None)
logging.getLogger('matplotlib').setLevel(logging.CRITICAL)


def intersect_with_atacseq(df_agg_intersect, output_dir, atacseq_file):
    # load atacseq file
    df_atacseq = pd.read_csv(atacseq_file, sep='\t')
    print(df_atacseq.shape)

    # match column names with the 10X tables
    # TODO: ask what to do with differetn column name
    df_atacseq = df_atacseq.rename(
        columns={'Region': '#chromosome', 'Position': 'start', 'Strand': 'Strand (0:-, 1:+, 2:unknown)'})

    # the atacseq is
    df_merged = pd.merge(df_agg_intersect, df_atacseq, on=['#chromosome', 'start'], how='left')
    df_merged_temp = df_merged[df_merged['gFrequency'].notnull()]

    # replace missing values with 0
    df_merged_temp['gCoverage-q20'] = df_merged_temp['gCoverage-q20'].replace('-', 0).astype(int)
    df_merged_temp['gFrequency'] = df_merged_temp['gFrequency'].replace('-', 0).astype(float)

    df_merged.to_csv(os.path.join(output_dir, '4.aggregated_per_position_intersect.bed'), sep='\t', index=False)


def plot_venn_diagram(df, subset_list, labels, column_name, output_dir, sname):
    plt.title('Intersection of positions - {} and {}'.format(labels[1], labels[0]))

    v = venn3(subsets=subset_list, set_labels=(labels[0], labels[1], labels[2]))
    venn3_circles(subsets=subset_list, color='gray', linewidth=1, linestyle='dashed')

    # get the text from the diagram components
    ids = ['100', '010', '001', '110', '101', '011', '111']
    sets = Counter()
    for id in ids:
        try:
            sets[id] = v.get_label_by_id(id).get_text()
        except:
            continue

    # change the position of the text on the figure
    h, l = [], []
    for i in sets:
        v.get_label_by_id(i).set_text("")  # remove label by setting them to empty string:
        h.append(v.get_patch_by_id(i))  # append patch to handles list
        l.append(sets[i])  # append count to labels list

    # create legend from handles and labels, and save figure
    plt.legend(handles=h, labels=l, title="counts", loc='lower right')  # bbox_to_anchor=(0.95,0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '5.venn_diagram_{}.png'.format(labels[0])), facecolor='white')
    plt.clf()

    # make histogram of mutated CB
    intrsct_list = df[df[column_name] == 1]['count of mutated cell barcodes']

    # histogram on log scale.
    # Use non-equal bin sizes, such that they look equal on log scale.
    hist, bins = np.histogram(intrsct_list, bins=30)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    plt.hist(intrsct_list, bins=logbins)
    plt.hist([i for i in intrsct_list if i >= 5], bins=logbins, alpha=0.6)
    plt.xscale('log', base=10)
    plt.yscale('log', base=10)

    plt.title("Number of mutated cells in intersection positions between table and {} - {}".format(labels[0], sname),
              fontsize=10)
    plt.ylabel("number of positions")
    plt.xlabel("number of mutated cells within position")
    plt.legend(['Intersection positions - not filtered', 'Intersection positions - filtered'])
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '5.CB_intersections_histogram_{}.png'.format(labels[0])), facecolor='white')
    plt.clf()


def run_venn(df, df_filtered, column_name, db_total_count, labels, input_dir, sname):
    # DB position set is combination of positions from table, and strings representing non overlaping positions.
    set1 = set(df[df[column_name] != 0].position.to_list() +
               ['not_in_table_position' + str(i) for i in range(db_total_count)])
    set2 = set(df.position)  # aggregated data
    set3 = set(df_filtered.position)  # filtered aggregated data

    plot_venn_diagram(df, [set1, set2, set3], labels, column_name, input_dir, sname)


def plot_venn2_diagram(subset_list, labels, output_dir, sname):
    # create hisrogram without green
    v = venn2(subsets=subset_list, set_labels=(labels[0], labels[1]))
    venn2_circles(subsets=subset_list, color='gray', linewidth=1, linestyle='dashed')

    sets = Counter()
    sets['10'] = v.get_label_by_id('10').get_text()
    sets['01'] = v.get_label_by_id('01').get_text()
    sets['11'] = v.get_label_by_id('11').get_text()

    h, l = [], []
    for i in sets:
        # remove label by setting them to empty string:
        v.get_label_by_id(i).set_text("")
        # append patch to handles list
        h.append(v.get_patch_by_id(i))
        # append count to labels list
        l.append(sets[i])

    # create legend from handles and labels
    plt.title('Intersection of {} and {} - {}'.format(labels[0], labels[1], sname))
    plt.legend(handles=h, labels=l, title="counts", bbox_to_anchor=(0.95, 0.7))
    plt.savefig(os.path.join(output_dir, 'venn2_diagram_{}.png'.format(labels[0])), facecolor='white')
    plt.clf()


def plot_heatmap_mutation_per_base_DB(df_merged, df_merged_filtered, output_dir, sname):
    # def get_min_max(count_matrices):
    #     """helper function to find the common min and max values for color scaling for all heatmaps in figure"""
    #     vmin, vmax = np.Inf, np.NINF
    #     for mat_to_plot in count_matrices:
    #         if mat_to_plot.min() < vmin:
    #             vmin = mat_to_plot.min()
    #         if mat_to_plot.max() > vmax:
    #             vmax = mat_to_plot.max()
    #
    #     return np.floor(np.log10(vmin)), np.ceil(np.log10(vmax))

    # def make_mut_counts_heatmap(count_matrices, out_folder, sname):
    #     """helper function to plot the heatmap"""
    #     # get min and max values for plotting color scale
    #     vmin, vmax = get_min_max(count_matrices)
    #
    #     fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    #     axs = axs.reshape(-1)
    #     for i, count_matrix in enumerate(count_matrices):
    #         # add first column of unmut data to the mut data
    #         mat_to_plot = count_matrix[0]
    #         axs[i] = sns.heatmap(np.log10(mat_to_plot), linewidth=0.5, annot=np.array(mat_to_plot),
    #                              cbar_kws={'label': 'log 10'}, ax=axs[i], cmap='brg', vmin=vmin, vmax=vmax,
    #                              xticklabels=['same_all', 'same_mut', 'transition', 'reverse', 'transversion'],
    #                              yticklabels=bases)
    #         axs[i].set_yticklabels(axs[i].get_yticklabels(), rotation=360)
    #         axs[i].set_xticklabels(axs[i].get_xticklabels(), rotation=30, ha='right')
    #         axs[i].set_title(
    #             "Counts of mutations per base - {} reads {} - {}".format(count_matrix[1], count_matrix[2], sname))
    #         axs[i].set_ylabel("reference base")
    #         axs[i].set_xlabel("mutation")
    #     plt.tight_layout()
    #     plt.savefig(os.path.join(out_folder, "heatmap_mutation_perbase_intersection.png"), bbox_inches='tight')

    def make_counts_matrix():
        """helper function to create two matices with counts of different umis, one with unmutated data and one with
        unmutated data"""
        umi_cols = ['same multi reads', 'transition multi reads', 'reverse multi reads', 'transvertion multi reads',
                    'same single reads', 'transition single reads', 'reverse single reads', 'transvertion single reads']
        count_matrices = []
        for i, df_tuple in enumerate(zip([df_merged, df_merged_filtered], ['', "- filtered"])):
            df, df_name = df_tuple[0], df_tuple[1]
            for j, read_type in enumerate(['single', 'multi']):
                count_matrix = []
                ref_umi_cols = [col for col in umi_cols if read_type in col]
                for base in bases:
                    idx = df[(df['is_editing'] == 1) & (df['reference base'] == base)].index
                    # idx = df[((df['is_editing_rep'] == 1) | (df['is_editing_non_rep'] == 1)) & (
                    #         df['reference base'] == base)].index
                    df_to_plot = df.loc[idx, ref_umi_cols].sum(axis=0)

                    # add count of 'same' umis in both mutated and un mutated
                    df_by_refbase = df[(df['is_editing'] == 1) & (df['reference base'] == base)]
                    # df_by_refbase = df[((df['is_editing_rep'] == 1) | (df['is_editing_non_rep'] == 1)) & (
                    #         df['reference base'] == base)]
                    unmuteted_read_count = df_by_refbase.drop_duplicates(subset='position')[
                        'unmutated {} reads'.format(read_type)].sum()
                    df_to_plot = pd.concat(
                        [pd.Series(df_to_plot['same {} reads'.format(read_type)] + unmuteted_read_count,
                                   index=['same all single reads']), df_to_plot])

                    count_matrix.append(df_to_plot.values)
                count_matrices.append((np.array(count_matrix), read_type, df_name))
        return count_matrices

    bases = ['a', 'c', 'g', 't']

    # create matrix with counts of mutations observed
    count_matrices = make_counts_matrix()

    # plot and save heatmap
    make_mut_counts_heatmap(count_matrices, output_dir, sname)

    # plot only A base mutations
    plt.clf()
    all_a_mutations_mat = np.zeros((4, 5))
    for i, count_matrix_set in enumerate(count_matrices):
        count_matrix = count_matrix_set[0]
        a_mut_data = count_matrix.sum(axis=0)
        all_a_mutations_mat[i, :] = a_mut_data

    vmin, vmax = all_a_mutations_mat.min(), all_a_mutations_mat.max()
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
    editing_total_count = os.popen("cat {} | wc -l".format(editing_db_path)).read()

    # convert to int
    snp_total_count, editing_total_count = int(snp_total_count), int(editing_total_count)

    # make Venn diagrams for snp, editing rep and editing non_rep intersections
    run_venn(df_agg_intrsct, df_filtered, 'is_snp', snp_total_count, ['SNP_DB', 'Aggregated data', 'Filtered data'],
             output_dir, sname)
    run_venn(df_agg_intrsct, df_filtered, 'is_editing', editing_total_count,
             ['Editing_DB', 'Aggregated data', 'Filtered data'], output_dir, sname)


def print_frequencies(df_merged, df_merged_agg, output_folder):
    with open(os.path.join(output_folder, '5.general_numbers.txt'), 'w') as f:
        f.write("General numbers information of table in open mode:\n\n")
        f.write("number of unique (position, CB) in table: %s \n" % str(df_merged.shape[0]))
        f.write("number of unique positions: %s \n" % str(df_merged['position'].nunique()))
        f.write("number of unique cell barcodes: %s \n" % str(df_merged['cell barcode'].nunique()))


def filter_open_and_agg_tables(df, df_agg, min_mutation_cb_to_filter, min_mutation_umis, min_total_umis, min_mutation_rate):
    # filter aggregated table
    df_agg_filt = sc_rna_variants.analysis_utils.filter_positions(df_agg, min_mutation_cb_to_filter, min_mutation_umis, min_total_umis, min_mutation_rate)

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


def drop_ifs(df):
    """remove positions which appear in both editing and snp sites"""
    idx = df[((df['is_editing_non_rep'] == 1) | (df['is_editing_rep'] == 1)) & (df['is_snp'] == 1)].index
    print("number of positin with snp and editing overlap to remove: %d." % len(idx))
    return df.drop(idx)


def run_snp_edit_DB_intersections(input_dir, output_dir, snp_db_path, editing_db_path, sname, atacseq):
    # get the df with intersections, before and after filtering
    df_agg_intersect, df_agg_intrsct_filtered = sc_rna_variants.analysis_utils.get_df_and_filtered_df(os.path.join(input_dir, '4.aggregated_per_position_intersect.bed'))

    # make Venn diagrams of the intersections
    make_venn_diagrams(df_agg_intersect, df_agg_intrsct_filtered, output_dir, snp_db_path, editing_db_path, sname)

    plot_heatmap_mutation_per_base_DB(df_agg_intersect, df_agg_intrsct_filtered, output_dir, sname)

    # if ATACseq data if supplied, remove potential SNP sites
    if (atacseq):
        intersect_with_atacseq(df_agg_intersect, output_dir, atacseq)


def run_step5(args):
    df_merged_agg = sc_rna_variants.analysis_utils.load_df(os.path.join(args.input_dir, "4.aggregated_per_position.bed"))
    # TODO: instead of loading mut_open and unmutated and merge them, look into the plots function and see how we can avoid this
    df_mut_open = sc_rna_variants.analysis_utils.load_tables(os.path.join(os.path.dirname(args.output_dir), 'step3_mismatch_dictionary', '3.mismatch_dictionary.bed'), mutated=True)
    df_unmutated = sc_rna_variants.analysis_utils.load_tables(os.path.join(os.path.dirname(args.output_dir), 'step3_mismatch_dictionary', '3.no_mismatch_dictionary.bed'), mutated=False)
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

    # TODO: remove snp overlaps with editing sites
    drop_ifs(df)

##################################################################################################################
def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(formatter_class=sc_rna_variants.utils.ArgparserFormater, description="",)

    # positional arguments
    parser.add_argument('input_dir', type=sc_rna_variants.utils.assert_is_directory, help='step 4 output folder')
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
    parser.add_argument('--log-file', default=os.path.join(sys.argv[2], "5.filtering_positions_and_snp_editing_DB_intersections.log"),
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
