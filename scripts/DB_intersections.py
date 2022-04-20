import os
import argparse

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles
from stats.generator import filter_by_cb_count
pd.set_option('display.max_columns', None)


def intersect_with_atacseq(df_agg_intersect, input_dir, atacseq_file):
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

    df_merged.to_csv(os.path.join(args.input_dir, 'aggregated_intersect.tsv'), sep='\t', index=False)


def find_intersections_with_SNP_and_edit_DB(args):
    # define paths
    agg_df_path = os.path.join(args.input_dir, 'aggregated_tsv.tsv')
    snp_temp_path = os.path.join(args.input_dir, 'snp_intersect.tsv')
    edit_temp_path = os.path.join(args.input_dir, 'edit_intersect.tsv')
    df_intersection = os.path.join(args.input_dir, 'aggregated_intersect.tsv')
    snp_db_path = args.snp_db_path
    rep_db_path = args.rep_db_path
    non_rep_db_path = args.non_rep_db_path

    # add '#' to header of df_aggregated
    os.system(f"head -c 1 {agg_df_path} | grep -q '#' || sed -i '1s/^/#/' {agg_df_path}")

    # both files must be sorted if you use '-sorted' which reduce memory usage
    # find intersection with snp db
    os.system(f"bedtools intersect -c -header -sorted -a {agg_df_path} -b {snp_db_path} > {snp_temp_path}")

    # add column name 'is_snp'
    os.system(f"sed -i '1 s/.*/&\tis_snp/' {snp_temp_path}")

    # find intersection with editing rep db
    os.system(f"bedtools intersect -s -c -header -a {snp_temp_path} -b {rep_db_path} > {edit_temp_path}")

    # add column name 'is_editing_rep'
    os.system(f"sed -i '1 s/.*/&\tis_editing_rep/' {edit_temp_path}")

    # find intersection with editing non rep db
    os.system(f"bedtools intersect -s -c -header -a {edit_temp_path} -b {non_rep_db_path} > {df_intersection}")

    # add column name 'is_editing_non_rep'
    os.system("sed -i '1 s/.*/&\tis_editing_non_rep/' {df_intersection}")

    # remove temp files
    os.system(f"rm {snp_temp_path} {edit_temp_path}")


# def make_pathes(input_dir):
#     """make pathes to temporary and output files"""
#     agg_df_path = os.path.join(input_dir, 'aggregated_tsv.tsv')
#     snp_temp = os.path.join(input_dir, 'snp_intersect.tsv')
#     edit_temp = os.path.join(input_dir, 'edit_intersect.tsv')
#     df_intersection = os.path.join(input_dir, 'aggregated_intersect.tsv')
#     return agg_df_path, snp_temp, edit_temp, df_intersection


def plot_venn_diagram(df, subset_list, labels, column_name, input_dir, sname):
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
    plt.savefig(os.path.join(input_dir, 'venn_diagram_{}.png'.format(labels[0])), facecolor='white')
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
    plt.savefig(os.path.join(input_dir, 'CB_intersections_histogram_{}.png'.format(labels[0])), facecolor='white')
    plt.clf()


def run_venn(df, df_filtered, column_name, db_total_count, labels, input_dir, sname):
    # DB position set is combination of positions from table, and strings representing non overlaping positions.
    set1 = set(df[df[column_name] != 0].position.to_list() +
               ['not_in_table_position' + str(i) for i in range(db_total_count)])
    set2 = set(df.position)  # aggregated data
    set3 = set(df_filtered.position)  # filtered aggregated data

    plot_venn_diagram(df, [set1, set2, set3], labels, column_name, input_dir, sname)


def plot_venn2_diagram(subset_list, labels, input_dir, sname):
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
    plt.savefig(os.path.join(input_dir, 'venn2_diagram_{}.png'.format(labels[0])), facecolor='white')
    plt.clf()


def plot_heatmap_mutation_per_base(df_merged, df_merged_filtered, out_folder, sname):
    def get_min_max(count_matrices):
        """helper function to find the common min and max values for color scaling for all heatmaps in figure"""
        vmin, vmax = np.Inf, np.NINF
        for mat_to_plot in count_matrices:
            if mat_to_plot.min() < vmin:
                vmin = mat_to_plot.min()
            if mat_to_plot.max() > vmax:
                vmax = mat_to_plot.max()

        return np.floor(np.log10(vmin)), np.ceil(np.log10(vmax))

    def make_mut_counts_heatmap(count_matrices, out_folder, sname):
        """helper function to plot the heatmap"""
        # get min and max values for plotting color scale
        vmin, vmax = get_min_max(count_matrices)

        fig, axs = plt.subplots(2, 2, figsize=(10, 8))
        axs = axs.reshape(-1)
        for i, count_matrix in enumerate(count_matrices):
            # add first column of unmut data to the mut data
            mat_to_plot = count_matrix[0]
            axs[i] = sns.heatmap(np.log10(mat_to_plot), linewidth=0.5, annot=np.array(mat_to_plot),
                                 cbar_kws={'label': 'log 10'}, ax=axs[i], cmap='brg', vmin=vmin, vmax=vmax,
                                 xticklabels=['same_all', 'same_mut', 'transition', 'reverse', 'transversion'],
                                 yticklabels=bases)
            axs[i].set_yticklabels(axs[i].get_yticklabels(), rotation=360)
            axs[i].set_xticklabels(axs[i].get_xticklabels(), rotation=30, ha='right')
            axs[i].set_title(
                "Counts of mutations per base - {} reads {} - {}".format(count_matrix[1], count_matrix[2], sname))
            axs[i].set_ylabel("reference base")
            axs[i].set_xlabel("mutation")
        plt.tight_layout()
        plt.savefig(os.path.join(out_folder, "heatmap_mutation_perbase_intersection.png"), bbox_inches='tight')

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
                    idx = df[((df['is_editing_rep'] == 1) | (df['is_editing_non_rep'] == 1)) & (
                            df['reference base'] == base)].index
                    df_to_plot = df.loc[idx, ref_umi_cols].sum(axis=0)

                    # add count of 'same' umis in both mutated and un mutated
                    df_by_refbase = df[((df['is_editing_rep'] == 1) | (df['is_editing_non_rep'] == 1)) & (
                            df['reference base'] == base)]
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
    make_mut_counts_heatmap(count_matrices, out_folder, sname)

    plt.clf()
    all_a_mutations_mat = np.zeros((4, 5))
    for i, count_matrix_set in enumerate(count_matrices):
        count_matrix = count_matrix_set[0]
        a_mut_data = count_matrix.sum(axis=0)
        all_a_mutations_mat[i, :] = a_mut_data

    # get min and max values for plotting color scale
    vmin, vmax = get_min_max(all_a_mutations_mat)
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
    plt.savefig(os.path.join(out_folder, "heatmap_A_mutations.png"), bbox_inches='tight')


def count_intersection(df_agg_intrsct, df_filtered, args):
    """function to get information on intersections of tables with databases"""
    snp_total_count = os.popen("grep -v '#' {} | wc -l".format(args.snp_db_path)).read()  # count non header lines
    edit_rep_total_count = os.popen("cat {} | wc -l".format(args.rep_db_path)).read()
    edit_nonrep_total_count = os.popen("cat {} | wc -l".format(args.non_rep_db_path)).read()

    # convert to int
    snp_total_count, edit_rep_total_count, edit_nonrep_total_count = \
        int(snp_total_count), int(edit_rep_total_count), int(edit_nonrep_total_count)

    # make Venn diagram for rep and non_rep Editing sites combined
    # DB position set is combination of positions from table, and strings representing non overlaping positions.
    edit_db_total_count = edit_rep_total_count + edit_nonrep_total_count
    set1 = set(df_agg_intrsct[(df_agg_intrsct['is_editing_rep'] != 0) | (
            df_agg_intrsct['is_editing_non_rep'] != 0)].position.to_list() +
               ['not_in_table_position' + str(i) for i in range(edit_db_total_count)])
    set2 = set(df_filtered.position)
    plot_venn2_diagram([set1, set2], ['Editing_sites_DB', 'Filtered_data'], args.input_dir, args.sname)

    # make Venn diagrams for snp, editing rep and editing non_rep intersections
    run_venn(df_agg_intrsct, df_filtered, 'is_snp', snp_total_count, ['SNP_DB', 'Aggregated data', 'Filtered data'],
             args.input_dir, args.sname)
    run_venn(df_agg_intrsct, df_filtered, 'is_editing_rep', edit_rep_total_count,
             ['Edit_rep_DB', 'Aggregated data', 'Filtered data'], args.input_dir, args.sname)
    run_venn(df_agg_intrsct, df_filtered, 'is_editing_non_rep', edit_nonrep_total_count,
             ['Edit_non_rep_DB', 'Aggregated data', 'Filtered data'], args.input_dir, args.sname)


def get_df(input_dir):
    """function to load the df and filter it"""

    # def get_filtered(df_agg, min_mutation_cb_to_filter, min_mutation_umis, min_total_umis, min_mutation_rate):
    #     #  'true values' - drop positions with rare mutations and probably hard to get insights from
    #     def filter_rare_mut(df, min_mutation_rate):
    #         df = df[df['percent of non ref from all cells'] > min_mutation_rate]
    #         return df
    #
    #     """function to return filtered aggregated table"""
    #     # first condition to filter by
    #     cond_1 = (df_agg['count of mutated cell barcodes'] >= min_mutation_cb_to_filter)
    #
    #     # second condition to filter by
    #     mutation_umi_counts = df_agg['total mutation umi count']
    #     total_umi_count = mutation_umi_counts + \
    #                       df_agg['unmutated multi reads'] + \
    #                       df_agg['unmutated single reads']
    #     cond_2 = ((mutation_umi_counts >= min_mutation_umis) & (total_umi_count >= min_total_umis))
    #
    #     # filter the aggregated df
    #     df_agg_filt = df_agg[cond_1 & cond_2]
    #     df_agg_filt = filter_rare_mut(df_agg_filt, min_mutation_rate)
    #
    #     return df_agg_filt

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

    # extract paths
    df_intersection_path = os.path.join(input_dir, 'aggregated_intersect.tsv')

    # load and df with intersections notations
    df_agg_intrsct = pd.read_csv(df_intersection_path, sep='\t')

    # define intersections to be binary (1 - if any overlap with db occured, 0 otherwise)
    df_agg_intrsct.loc[df_agg_intrsct['is_snp'] > 0, 'is_snp'] = 1
    df_agg_intrsct.loc[df_agg_intrsct['is_editing_rep'] > 0, 'is_editing_rep'] = 1
    df_agg_intrsct.loc[df_agg_intrsct['is_editing_non_rep'] > 0, 'is_editing_non_rep'] = 1
    # df_agg_intrsct = get_edit_intersections(df_agg_intrsct, 'is_editing_rep')  # editing sites need addition filteration
    # df_agg_intrsct = get_edit_intersections(df_agg_intrsct, 'is_editing_non_rep')  # editing sites need addition filteration

    # save df
    df_agg_intrsct.to_csv(df_intersection_path, sep='\t', index=False)

    # get filtered df
    _, df_agg_intrsct_filtered = filter_by_cb_count(df_agg=df_agg_intrsct,
                                           min_mutation_cb_to_filter=5,
                                           min_mutation_umis=10,
                                           min_total_umis=20,
                                           min_mutation_rate=0.1)

    return df_agg_intrsct, df_agg_intrsct_filtered


def run(args):
    # get paths to files we use
    # pathes = make_pathes(args.input_dir)

    # find intersection between df and databases
    find_intersections_with_SNP_and_edit_DB(args)

    # get the df with intersections, before and after filtering
    df_agg_intersect, df_agg_intrsct_filtered = get_df(args.input_dir)

    # make Venn diagrams of the intersections
    count_intersection(df_agg_intersect, df_agg_intrsct_filtered, args)

    plot_heatmap_mutation_per_base(df_agg_intersect, df_agg_intrsct_filtered, args.input_dir, args.sname)

    # if ATACseq data if supplied, remove ppotential SNP sites
    if (args.atacseq):
        intersect_with_atacseq(df_agg_intersect, args.input_dir, args.atacseq)


#####################################################################################################
def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(
        description="""Find intersection between the aggregated file and already known SNP and editing sites.""",
        epilog='''Outputs table with new intersection columns and creates venn diagrams.'''
    )

    # positional arguments
    parser.add_argument('input_dir', help='path to folder with outputs from make_statistics.py')

    parser.add_argument('rep_db_path', help='path to gencode file with known repetitive editing sites')

    parser.add_argument('non_rep_db_path', help='path to gencode file with known non repetitive editing sites')

    parser.add_argument('snp_db_path', help='path to gencode file with known SNP sites')

    # optional arguments
    parser.add_argument('--sname', type=str, help='sample name to add to outputs')

    parser.add_argument('--atacseq', type=str, help='path to atacseq file')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    args = parse_arguments()
    run(args)
    print("Fininshed DB_intersection")
