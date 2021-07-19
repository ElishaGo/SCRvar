import os
import numpy as np
import pandas as pd
from math import ceil
import seaborn as sns
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")


def plot_cb_occurences_hist(df_orig, df_filtered, out_folder):
    """ This function count the number of each Cb in the table, and then plot the distribution of CB occurences.
    This can give an insight regarding the representation of different cell in the table.
    Note that the counts in tables can overlap on positions."""
    fig, axs = plt.subplots(2, 1, figsize=(10, 12), sharex=True, sharey=True)

    def make_plot(ax, cb_counts, df_name):
        ax.hist(cb_counts.values, bins=50)
        ax.set_title("Number of cell barcodes - {} filtering".format(df_name))
        ax.set_ylabel("count of different cell barcodes")
        ax.set_xlabel("number of mutated positions")

    for i, df_tuple in enumerate(zip([df_orig, df_filtered], ['prior', 'post'])):
        df = df_tuple[0]
        df_name = df_tuple[1]
        cb_counts = df['cell barcode'].value_counts()
        make_plot(axs[i], cb_counts, df_name)
    plt.savefig(os.path.join(out_folder, "cb_distribution.png"), bbox_inches='tight')


def plot_umi_per_reference_base(df_merged, df_merged_filtered, out_folder):
    fig, axs = plt.subplots(2, 4, figsize=(15, 8), sharex=True)
    bases = ['a', 'c', 'g', 't']
    ref_umi_cols = ['same multi reads', 'transition multi reads', 'reverse multi reads', 'transvertion multi reads',
                    'same single reads', 'transition single reads', 'reverse single reads', 'transvertion single reads']
    # ref_umi_cols = [col for col in df_merged.columns if col.startswith('R->')]

    for i, df in enumerate([df_merged, df_merged_filtered]):
        for j, base in enumerate(bases):
            idx = df.loc[df['reference base'] == base].index
            df.loc[idx, ref_umi_cols].sum(axis=0).plot(kind='barh', ax=axs[i][j], sharey=True)
            axs[i][j].set_title('refference base: {0}'.format(base))
            axs[i][j].set_xlabel("Umi counts")
            axs[i][j].set_xscale('log', base=10)
            for z in axs[i][j].patches:
                axs[i][j].text(z.get_width() + .09, z.get_y() + .3, str(round((z.get_width()), 1)), fontsize=8)
    fig.text(0, 0.65, 'Before filtering', ha='center', rotation='vertical', style='italic')
    fig.text(0, 0.25, 'After filtering', ha='center', rotation='vertical', style='italic')
    plt.savefig(os.path.join(out_folder, "umi_per_reference_base.png"), bbox_inches='tight')

    ### plot with unmutated data
    fig, axs = plt.subplots(2, 4, figsize=(15, 8), sharex=True)
    for i, df in enumerate([df_merged, df_merged_filtered]):
        for j, base in enumerate(bases):
            idx = df.loc[df['reference base'] == base].index
            df_to_plot = df.loc[idx, ref_umi_cols].sum(axis=0)
            df_by_refbase = df.loc[df['reference base'] == base, :]
            unmuteted_multi_read_count = df_by_refbase.drop_duplicates(subset='position')[
                'unmutated multi reads'].sum()
            unmuteted_single_read_count = df_by_refbase.drop_duplicates(subset='position')[
                'unmutated single reads'].sum()
            df_to_plot['same multi reads'] += unmuteted_multi_read_count
            df_to_plot['same single reads'] += unmuteted_single_read_count
            df_to_plot.plot(kind='barh', ax=axs[i][j], sharey=True)
            axs[i][j].set_title('refference base: {0}'.format(base))
            axs[i][j].set_xscale('log', base=10)
            for z in axs[i][j].patches:
                axs[i][j].text(z.get_width() + .09, z.get_y() + .3, str(round((z.get_width()), 1)), fontsize=8)
    fig.text(0, 0.65, 'Before filtering', ha='center', rotation='vertical', style='italic')
    fig.text(0, 0.25, 'After filtering', ha='center', rotation='vertical', style='italic')
    plt.savefig(os.path.join(out_folder, "umi_per_reference_base_with_unmutated.png"), bbox_inches='tight')


def plot_heatmap_mutation_per_base(df_merged, df_merged_filtered, out_folder):
    def get_min_max(count_matrices):
        """helper function to find the common min and max values for color scaling for all heatmaps in figure"""
        vmin, vmax = np.Inf, np.NINF
        for count_matrix in count_matrices:
            for row in count_matrix[0]:
                if row.min() < vmin:
                    vmin = row.min()
                if row.max() > vmax:
                    vmax = row.max()
        return np.floor(np.log10(vmin)), np.ceil(np.log10(vmax))

    def make_plot(count_matrices):
        """helper function to plot the heatmap"""
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))
        axs = axs.reshape(-1)
        vmin, vmax = get_min_max(count_matrices)
        for i, count_matrix in enumerate(count_matrices):
            axs[i] = sns.heatmap(np.log10(count_matrix[0]), linewidth=0.5, annot=np.array(count_matrix[0]),
                                 cbar_kws={'label': 'log 10'}, ax=axs[i], cmap='coolwarm', vmin=vmin, vmax=vmax,
                                 xticklabels=['same', 'transition', 'reverse','transversion'], yticklabels=bases)
            axs[i].set_yticklabels(axs[i].get_yticklabels(), rotation=360)
            axs[i].set_title("Counts of mutations per base - {} reads {}".format(count_matrix[1], count_matrix[2]))
            axs[i].set_ylabel("reference base")
            axs[i].set_xlabel("mutation")
        plt.tight_layout()

    def create_figure(with_unmutated, suffix=""):
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
                    idx = df[df['reference base'] == base].index
                    df_to_plot = df.loc[idx, ref_umi_cols].sum(axis=0)
                    if with_unmutated:
                        df_by_refbase = df.loc[df['reference base'] == base, :]
                        unmuteted_read_count = df_by_refbase.drop_duplicates(subset='position')[
                            'unmutated {} reads'.format(read_type)].sum()
                        df_to_plot['same {} reads'.format(read_type)] += unmuteted_read_count
                    count_matrix.append(df_to_plot.values)
                count_matrices.append((count_matrix, read_type, df_name))
        # count_matrices = []
        # for i, df_tuple in enumerate(zip([df_merged, df_merged_filtered], ['', "- filtered"])):
        #     df, df_name = df_tuple[0], df_tuple[1]
        #     for j, read_type in enumerate(['single', 'multi']):
        #         count_matrix = []
        #         ref_umi_cols = [col for col in df.columns if (col.startswith('R->') and read_type in col)]
        #         for base in bases:
        #             idx = df[df['reference base'] == base].index
        #             df_to_plot = df.loc[idx, ref_umi_cols].sum(axis=0)
        #             if with_unmutated:
        #                 df_by_refbase = df.loc[df['reference base'] == base, :]
        #                 unmuteted_read_count = df_by_refbase.drop_duplicates(subset='position')['unmutated {} reads'.format(read_type)].sum()
        #                 df_to_plot['R->{} {} reads'.format(base.upper(), read_type)] += unmuteted_read_count
        #             count_matrix.append(df_to_plot.values)
        #         count_matrices.append((count_matrix, read_type, df_name))
        make_plot(count_matrices)
        plt.savefig(os.path.join(out_folder, "heatmap_mutation_perbase{}.png".format(suffix)), bbox_inches='tight')

    bases = ['a', 'c', 'g', 't']
    create_figure(with_unmutated=False)
    create_figure(with_unmutated=True, suffix="_with_unmutated")


def plot_cb_count_overall(df_merged_agg, df_merged_agg_filtered, out_folder):
    """plot count of all CB overall occurences (no unique) per unique postion"""

    def make_plot(ax, cb_counts, df_name, num_of_bins):
        """helper function to plot a histogram"""
        counts, bins = np.histogram(cb_counts, bins=max(cb_counts))
        ax.bar(bins[:num_of_bins], bins[:num_of_bins] * counts[:num_of_bins])
        for i in range(num_of_bins):  # add ticks on top of the bars
            ax.text(bins[i], bins.round()[i] * counts[i], str(int(bins.round()[i])) + '*' + str(counts[i]), rotation=45)
        ax.set_title("count of all cell barcodes (with repetitions) per position - {} filtering".format(df_name))
        ax.set_ylabel("count of positions")
        ax.set_xlabel("number of all Cell Barcodes per position")
        ax.tick_params(axis='x', reset=True)  # show ticks of x axis on both graphs

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(20, 12))
    for i, df_tuple in enumerate(zip([df_merged_agg, df_merged_agg_filtered], ['prior', 'post'])):
        df, df_name = df_tuple[0], df_tuple[1]
        num_of_bins = 30
        cb_counts = df['count of mutated cell barcodes']  # count of CB in each position
        make_plot(axs[i], cb_counts, df_name, num_of_bins)
    plt.savefig(os.path.join(out_folder, "cb_count_not_unique_per_position.png"), bbox_inches='tight')


def plot_cb_count_per_position(df_merged_agg, df_merged_agg_filtered, out_folder):
    def make_plot(ax, cb_counts, df_name):
        """helper function to plot the histogram"""
        bins = [2 ** n for n in list(range(0, ceil(np.log2(max(cb_counts))) + 1))]
        arr = ax.hist(cb_counts.values, bins=bins, weights=np.ones(len(cb_counts.values)) / len(cb_counts.values))
        ax.set_title("{} filtering".format(df_name))
        ax.set_ylabel("count of positions%")
        ax.set_xlabel("number of different Cell Barcodes per position")
        ax.set_xscale('log', base=2)
        ax.set_ylim([0, 1])
        vals = ax.get_yticks()
        ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
        bin_length = len(bins)
        for j in range(bin_length - 1):
            ax.text((arr[1][j] + arr[1][j + 1]) / 2.5, arr[0][j] + 0.01, str(int((arr[0][j] * 100).round())) + '%')

    # plot mutated data
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 12))
    for i, df_tuple in enumerate(zip([df_merged_agg, df_merged_agg_filtered], ['Before', 'After'])):
        df, df_name = df_tuple[0], df_tuple[1]
        # df.to_csv(os.path.join(out_folder, "test1_df.tsv"))
        cb_counts = df['count of mutated cell barcodes']
        make_plot(axs[i], cb_counts, df_name)
    plt.savefig(os.path.join(out_folder, "cb_count_per_position.png"), bbox_inches='tight')

    # plot with unmutated data
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 12))
    for i, df_tuple in enumerate(zip([df_merged_agg, df_merged_agg_filtered], ['Before', 'After'])):
        df, df_name = df_tuple[0], df_tuple[1]
        cb_counts = df['count of mutated cell barcodes']  # mutated counts
        unmutated_cb_counts = df['count of unmutated cell barcodes']  # unmutated counts
        cb_counts = cb_counts.combine(unmutated_cb_counts, np.add, fill_value=0)  # sum mutated and unmutated counts
        make_plot(axs[i], cb_counts, df_name)
    plt.savefig(os.path.join(out_folder, "cb_count_per_position_with_unmutated.png"), bbox_inches='tight')


"""
    data / raw_Stats_files / removed / ADAR / concatenated_tsv / raw_stats_removed_ADAR.tsv
    data / raw_Stats_files / removed / ADAR / concatenated_tsv / raw_stats_unmutated_removed_ADAR.tsv
    --output_folder
    outputs / statistics_output / test
    --log - file
    log_files / log_test.txt
    --threads
    16
"""
