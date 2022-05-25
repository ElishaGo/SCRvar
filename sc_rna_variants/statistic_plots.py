import os
import warnings
import numpy as np
import pandas as pd
from math import ceil
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")


def plot_cb_occurences_hist(df_orig, df_filtered, out_folder, sname):
    """ This function count the number of each Cb in the table, and then plot the distribution of CB occurences.
    This can give an insight regarding the representation of different cell in the table.
    Note that the counts in tables can overlap on positions."""

    def make_plot(ax, cb_counts, prior_post):
        ax.hist(cb_counts.values, bins=50)
        ax.set_title("Number of cell barcodes - {} filtering - {}".format(prior_post, sname))
        ax.set_ylabel("count of different cell barcodes")
        ax.set_xlabel("number of mutated positions")
        ax.tick_params(axis='x', reset=True)  # show ticks of x axis on both graphs

    fig, axs = plt.subplots(2, 1, figsize=(10, 12), sharex=True, sharey=True)

    for i, df_tuple in enumerate(zip([df_orig, df_filtered], ['prior', 'post'])):
        df = df_tuple[0]
        prior_post = df_tuple[1]
        cb_counts = df['cell barcode'].value_counts()
        make_plot(axs[i], cb_counts, prior_post)
    plt.savefig(os.path.join(out_folder, "5.cb_distribution.png"), bbox_inches='tight')


def plot_umi_per_reference_base(df_merged, df_merged_filtered, out_folder, sname):
    fig, axs = plt.subplots(2, 4, figsize=(15, 8), sharex=True)
    bases = ['a', 'c', 'g', 't']
    ref_umi_cols = ['same multi reads', 'transition multi reads', 'reverse multi reads', 'transvertion multi reads',
                    'same single reads', 'transition single reads', 'reverse single reads', 'transvertion single reads']

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
    plt.suptitle(f"UMI per reference base - {sname}")
    plt.savefig(os.path.join(out_folder, "5.umi_per_reference_base.png"), bbox_inches='tight')

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
    plt.suptitle(f"UMI per reference base with unmutated - {sname}")
    plt.savefig(os.path.join(out_folder, "5.umi_per_reference_base_with_unmutated.png"), bbox_inches='tight')


def get_min_max(count_matrices):
    """helper function to find the common min and max values for color scaling for all heatmaps in figure"""
    vmin, vmax = np.Inf, np.NINF
    for count_matrix in count_matrices:
        mat_to_plot = count_matrix[0]
        if mat_to_plot.min() < vmin:
            vmin = mat_to_plot.min()
        if mat_to_plot.max() > vmax:
            vmax = mat_to_plot.max()

    return np.floor(np.log10(vmin)), np.ceil(np.log10(vmax))


def make_mut_counts_heatmap(count_matrices, out_folder, sname):
    """helper function to plot the heatmap"""
    # get min and max values for plotting color scale
    bases = ['a', 'c', 'g', 't']
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
        axs[i].set_title("{} reads {}".format(count_matrix[1], count_matrix[2]))
        axs[i].set_ylabel("reference base")
        axs[i].set_xlabel("mutation")
    plt.suptitle(f"Counts of mutations per base - {sname}")
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, "5.heatmap_mutation_perbase.png"), bbox_inches='tight')


def plot_heatmap_mutation_per_base(df_merged, df_merged_filtered, out_folder, sname):
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
                    idx = df[df['reference base'] == base].index
                    df_to_plot = df.loc[idx, ref_umi_cols].sum(axis=0)

                    # add count of 'same' unis in both mutated and un mutated
                    df_by_refbase = df.loc[df['reference base'] == base, :]
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

    relative_mut_dfs = []
    for count_matrix_set in count_matrices:
        count_matrix = count_matrix_set[0]
        type_name = count_matrix_set[1] + count_matrix_set[2]

        same_col = count_matrix[:, 1]  # get the 'same_mut' column
        relative_mutations = count_matrix[:, 2:] / same_col[:, None]
        relative_mutations_df = pd.DataFrame(relative_mutations,
                                             columns=['transition', 'reverse', 'transversion'],
                                             index=bases).reset_index()
        relative_mutations_long = pd.melt(relative_mutations_df, id_vars=['index'],
                                          value_vars=['transition', 'reverse', 'transversion'])
        relative_mutations_long['type'] = type_name
        relative_mutations_long['variable'] = relative_mutations_long['index'] + '_' + relative_mutations_long[
            'variable']
        relative_mutations_long.drop('index', axis=1, inplace=True)

        relative_mut_dfs.append(relative_mutations_long)

    relative_mut_dfs = pd.concat(relative_mut_dfs)

    # make barplot
    plt.clf()

    s = sns.barplot(data=relative_mut_dfs, x='variable', y='value', hue='type')
    s.set_yscale("log", base=10)  # set y axis in log scale
    _ = s.set(xlabel="mutated UMIs ratio against non mutated UMIs", ylabel="ratio")  # set axis labels
    # reset x labels
    mut_names = ['(a->g)', '(c->t)', '(g->a)', '(t->c)', '(a->t)', '(c->g)', '(g->c)',
                 '(t->a)', '(a->c)', '(c->a)', '(g->t)', '(t->g)']
    s.set_xticklabels([l.get_text() + '\n' + m for l, m in zip(s.get_xticklabels(), mut_names)])
    plt.xticks(rotation=30, ha='right')
    plt.title("ratio of mutation out of non mutated UMIs (same column)")
    plt.savefig(os.path.join(out_folder, "5.relative_mutations_barplot.png"), bbox_inches='tight')


def plot_cb_count_overall(df_merged_agg, df_merged_agg_filtered, out_folder, sname):
    """plot count of all CB overall occurences (no unique) per unique postion"""

    def make_plot(ax, cb_counts, df_name, num_of_bins):
        """helper function to plot a histogram"""
        counts, bins = np.histogram(cb_counts, bins=max(cb_counts))
        ax.bar(bins[:num_of_bins], bins[:num_of_bins] * counts[:num_of_bins])
        for i in range(num_of_bins):  # add ticks on top of the bars
            ax.text(bins[i], bins.round()[i] * counts[i], str(int(bins.round()[i])) + '*' + str(counts[i]), rotation=45)
        ax.set_title(
            "amount of (positions, CB) with a mismatch, binned by the number of cells - {} filtering".format(df_name))
        ax.set_ylabel("count of (positions, CB) pairs")
        ax.set_xlabel("number of all Cell Barcodes per position")
        ax.tick_params(axis='x', reset=True)  # show ticks of x axis on both graphs

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(20, 12))
    for i, df_tuple in enumerate(zip([df_merged_agg, df_merged_agg_filtered], ['prior', 'post'])):
        df, df_name = df_tuple[0], df_tuple[1]
        num_of_bins = 30
        cb_counts = df['count of mutated cell barcodes']  # count of CB in each position
        make_plot(axs[i], cb_counts, df_name, num_of_bins)
    plt.savefig(os.path.join(out_folder, "5.cb_count_not_unique_per_position.png"), bbox_inches='tight')


def plot_cb_count_per_position(df_merged_agg, df_merged_agg_filtered, out_folder, sname):
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
        ax.tick_params(axis='x', reset=True)  # show ticks of x axis on both graphs
        bin_length = len(bins)
        for j in range(bin_length - 1):
            ax.text((arr[1][j] + arr[1][j + 1]) / 2.5, arr[0][j] + 0.01, str(int((arr[0][j] * 100).round())) + '%')

    # plot mutated data
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 12))
    for i, df_tuple in enumerate(zip([df_merged_agg, df_merged_agg_filtered], ['Before', 'After'])):
        df, df_name = df_tuple[0], df_tuple[1]
        cb_counts = df['count of mutated cell barcodes']
        make_plot(axs[i], cb_counts, df_name)
    plt.savefig(os.path.join(out_folder, "5.cb_count_per_position.png"), bbox_inches='tight')

    # plot with unmutated data
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 12))
    for i, df_tuple in enumerate(zip([df_merged_agg, df_merged_agg_filtered], ['Before', 'After'])):
        df, df_name = df_tuple[0], df_tuple[1]
        cb_counts = df['count of mutated cell barcodes']  # mutated counts
        unmutated_cb_counts = df['count of unmutated cell barcodes']  # unmutated counts
        cb_counts = cb_counts.combine(unmutated_cb_counts, np.add, fill_value=0)  # sum mutated and unmutated counts
        make_plot(axs[i], cb_counts, df_name)
    plt.savefig(os.path.join(out_folder, "5.cb_count_per_position_with_unmutated.png"), bbox_inches='tight')


def non_ref_from_all_cells(df, output_dir):
    plt.clf()
    plt.hist(df['percent of non ref from all cells'], bins=21)
    plt.title('percent of non ref from all cells')
    plt.xlabel("percent of mutated UMIs in a position")
    plt.ylabel("number of observed positions")
    plt.savefig(os.path.join(output_dir, "6.mut_UMI_in_positions.png"))


def snp_observed_against_mut_UMI_in_position(df, output_dir, sname):
    plt.clf()
    plt.hist(df.loc[df['is_snp'] >= 1, :]['percent of non ref from all cells'], bins=30)
    plt.hist(df.loc[df['is_snp'] == 0, :]['percent of non ref from all cells'], bins=30, alpha=.7)

    plt.title("Percent of mutated umis from all umis per position, in SNP and non SNP positions \n - {}".format(sname))
    plt.legend(['SNP position', 'No SNP position'])
    plt.xlabel("percent of mutated UMIs in a position")
    plt.ylabel("number of observed positions")
    plt.savefig(os.path.join(output_dir, "6.snp_observed_against_mut_UMI_in_position.png"))


def editing_sites_per_chr(df_edit, output_dir, sname):
    chr_labels = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                  'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',
                  'chrY', 'chrM']

    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=False, figsize=(10, 8))
    chr_series = df_edit['position'].map(lambda x: x.split(':')[0]).value_counts()
    plt.barh(range(0, len(chr_labels)), pd.Series(chr_series, index=chr_labels), tick_label=chr_labels, color='b')
    for i, v in enumerate(pd.Series(chr_series, index=chr_labels).values):
        try:
            ax[0].text(v, i - 0.25, str(int(v)))
        except:
            continue
    plt.title("editing DB")
    plt.ylabel("chromosome")
    plt.xlabel("Number of editing sites")
    plt.savefig(os.path.join(output_dir, "6.editing_sites_per_chr.png"))


def editing_sites_per_cell(df_open_edit, output_dir, sname):
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(6, 4))

    editings_per_cb_stats = df_open_edit.iloc[:, :].groupby('cell barcode')['is_edit'].sum()
    print(
        "{} mutations compattable with known editing sites in \n{} cells out of total number of \n{} cells in sample\n\n".format(
            df_open_edit[df_open_edit['is_edit'] == 1].shape[0], (editings_per_cb_stats > 0).sum(),
            df_open_edit['cell barcode'].nunique()))
    plt.hist(editings_per_cb_stats, bins=30)
    plt.title("Number of editing sites per cell")
    plt.ylabel("number of cells")
    plt.xlabel("number of editing sites")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "6.editing_sites_per_cell.png"))


def mutated_umis_per_cell(umis_per_cb_editing, output_dir, sname):
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(6, 4))

    plt.hist(umis_per_cb_editing, bins=18)
    plt.xlabel("number of mutated UMI's in cell")
    plt.ylabel("observed cells")
    fig.suptitle("number of mutated umis per cell in known editing sites")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "6.umi_per_cell_in_editing_sites.png"))


def make_clusters_heatmap(pivot_table, name, min_umi_per_position, min_umi_per_cell, output_dir, chr_genes_pairs):
    plt.clf()
    # filter rows (positions) and columns (cells)
    df_temp = pivot_table[
        pivot_table.sum(axis=1) > min_umi_per_cell]  # keep cells with more than N mutated umis
    df_temp = df_temp.loc[:,
              df_temp.sum(axis=0) > min_umi_per_position]  # keep positions with more than N mutetaed umis
    print("shape of df is:", df_temp.shape)

    # define clusters and set color map for rows
    clusters = df_temp.index.map(lambda x: x.split('_')[0])  # use Seurat clusters - when rows are cells
    #     clusters = df_temp.index # use Seurat clusters - when rows are clusters
    print("number of colors on X axis:", clusters.nunique())
    gist_rainbow = cm.get_cmap('tab20c', clusters.nunique())
    clsters_lut = dict(zip(clusters.unique(), gist_rainbow(np.linspace(0, 1, clusters.nunique()))))
    cluster_colors = clusters.map(clsters_lut)

    # define gene names and set color map by chromosome
    gist_rainbow = cm.get_cmap('gist_rainbow', len(set(chr_genes_pairs.values())))
    genes_lut = dict(
        zip(set(chr_genes_pairs.values()), gist_rainbow(np.linspace(0, 1, len(set(chr_genes_pairs.values()))))))
    genes = df_temp.columns
    gene_colors = genes.map(chr_genes_pairs).map(genes_lut)  # map genes to chromosomes and chromosomes to colors

    # xticklabels=True, yticklabels=True
    # col_colors=cluster_colors,
    # distance metrics jaccard, euclidean
    cg = sns.clustermap(df_temp.fillna(0).T, metric='jaccard', method='complete', row_cluster=True, col_cluster=True,
                        col_colors=cluster_colors,
                        xticklabels=True, yticklabels=True, cmap='rocket_r', figsize=(30, 30))
    cg.ax_row_dendrogram.set_visible(True)
    cg.ax_col_dendrogram.set_visible(True)

    # color the gene names according to chromosomoes
    for tick, color in zip(cg.ax_heatmap.get_yticklabels(), gene_colors):
        c = genes_lut[chr_genes_pairs[tick.get_text()]]
        plt.setp(tick, color=c, size=18, x=1, ha='left')

    # color the Seurat cluster text
    #     for tick, color in zip(cg.ax_heatmap.get_xticklabels(), cluster_colors):
    #         c = clsters_lut[tick.get_text()]
    #         plt.setp(tick, color=c, y=1, va='baseline')  # 'va': top', 'bottom', 'center', 'baseline', 'center_baseline'

    plt.title(name)
    plt.savefig(os.path.join(output_dir, '6.editing_{}.png'.format(name)))


def plot_umis_per_gene(clusters_VS_genes_umis_pt, output_dir):
    fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True, sharex=True, figsize=(6, 4))

    plt.hist(clusters_VS_genes_umis_pt.sum(), bins=100)
    plt.xlabel("number of mutated UMI's in gene")
    plt.ylabel("observed genes")
    fig.suptitle("number of mutated UMI's in genes in known editing sites")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "6.umi_per_gene_in_editing_sites.png"))


def editing_events_vs_number_of_reads(df_scatter, output_dir):
    plt.clf()
    s = sns.lmplot(data=df_scatter, x="is_edit", y="number of reads per cell", height=7,
                   scatter_kws={'alpha': 0.5, 's': 14})
    s.set(xlim=(1, None))

    # add correlation
    def annotate(data, **kws):
        r, p = pearsonr(df_scatter['is_edit'], df_scatter['number of reads per cell'])
        ax = plt.gca()
        ax.text(.9, .05, 'pearson_r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes, ha='right')
    s.map_dataframe(annotate)

    plt.title("Editing events vs. number of reads, per cell barcode")
    plt.xlabel('number of editing events')
    plt.ylabel('number of reads per cell')
    # plt.yaxis.set_tick_params(labelbottom=True)
    # plt.legend([f"a CB. Number of CB - {len(df_scatter)}"], loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "6.editing_VS_reads_scatter.png"))


def editing_events_vs_number_of_mutated_umis_per_cell(editingsites_per_cb_stats, umis_per_cb_editing, output_dir):
    """
    the plot needs to be made from the count of all the umis in cells. this information exist in the bam file before any filtering.
    in addition, add the column of 0 editing sites (cell barcodes not in the editing table (with the 191 rows)).
    we want to see if there is correlation between the number of mutated umis per cell and the number of editing sites in the same cell
    If we will see cells which doesn't follow the correlation, will probably have more editing events than usual
    :param editingsites_per_cb_stats:
    :param output_dir:
    :return:
    """

    df_scatter = pd.merge(editingsites_per_cb_stats, umis_per_cb_editing, on='cell barcode')
    plt.clf()
    sns.scatterplot(
        data=df_scatter.groupby(['is_edit', 'mutated umis per cell']).size().reset_index(name='count of cells'),
        x="is_edit", y="mutated umis per cell", hue="count of cells", size='count of cells', sizes=(30, 400),
        legend="full", hue_norm=(0, 100))
    plt.xlabel('number of observed editing events (positions)')
    plt.ylabel('number of mutated UMIs per cell')
    plt.legend(title='#cells per position', loc='upper left', ncol=2)

    plt.suptitle("Editing events vs. number of mutated umis per cell barcode")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "editing_events_VS_umi_scatter.png"))
