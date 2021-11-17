import os
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles
pd.set_option('display.max_columns', None)


def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(
        description="""This script produces Venn diagrams to see intersection between the aggregated file and already 
        known editing sites.""",
        epilog='''Outputs venn diagrams.'''
    )

    # positional arguments
    parser.add_argument('input_dir', help='path to folder with outputs from make_statistics.py')

    # positional arguments
    parser.add_argument('rep_db_path', help='path to gencode file with known repetitive editing sites')

    # positional arguments
    parser.add_argument('non_rep_db_path', help='path to gencode file with known non repetitive editing sites')

    # positional arguments
    parser.add_argument('snp_db_path', help='path to gencode file with known SNP sites')

    parser.add_argument('--sname', type=str, help='sample name to add to outputs')
    
    return parser.parse_args(arguments)


def plot_venn_diagram(df, subset_list, labels, column_name, input_dir, sname):

    # fig, axs = plt.subplots(2, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]})
    plt.title('Intersection - {} and {} {}'.format(labels[0], labels[2], labels[1]))
    # axs[0].set_title("Venn diagram \n{a}: {b},   {c}: {d},   {e}: {f} \n{g}: {h},  {i}: {j},  {k}: {l}".format(
    #     a=labels[0], b=subset_list[0],
    #     c=labels[1], d=subset_list[1],
    #     e="intersection", f=subset_list[2],
    #     g=labels[2], h=subset_list[5] + subset_list[3],
    #     i=labels[0] + " & " + labels[2], j=subset_list[4],
    #     k="all_intersect", l=subset_list[6]), fontsize=10, bbox={"facecolor": "gray", "alpha": 0.3})

    v = venn3(subsets=subset_list, set_labels=(labels[0], labels[1], labels[2]))
    venn3_circles(subsets=subset_list, color='gray', linewidth=1, linestyle='dashed')

    # get the text from the diagram components
    ids = ['100', '010', '001', '110', '101','011', '111']
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
    hist, bins = np.histogram(intrsct_list, bins=20)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    plt.hist(intrsct_list, bins=logbins)
    plt.hist([i for i in intrsct_list if i >= 5], bins=logbins, alpha=0.6)
    plt.xscale('log')

    plt.title("Count of mutated cell barcodes in intersection with {} - {}".format(labels[0], sname), fontsize=10)
    plt.xlabel("mutations")
    plt.ylabel("cell barcodes")
    plt.legend(['not filtered', 'fitered'])
    plt.tight_layout()
    plt.savefig(os.path.join(input_dir, 'CB_intersections_histogram_{}.png'.format(labels[0])), facecolor='white')
    plt.clf()


def run_venn(df, df_filtered, column_name, db_total_count, labels, input_dir, sname):
    # DB position set is combination of positions from table, and strings representing non overlaping positions.
    set1 = set(df[df[column_name] != 0].position.to_list() +
               ['not_in_table_position' + str(i) for i in range(db_total_count)])
    set2 = set(df.position)  # aggregated data
    set3 = set(df_filtered.position)  # filtered aggregated data


    # df_pos_count = df.shape[0]
    # filtered_df_intrsct = len(df.index & df_filtered.index)
    # uniq_filtered_count = df_filtered.shape[0] - filtered_df_intrsct
    # db_df_intrsct = (df[column_name] != 0).sum()
    # db_filtered_intrsct = (df_filtered[column_name] != 0).sum()
    # db_filtered_df_intrsct = (df_filtered[column_name] != 0).sum()
    # # {'001':, len(df.position) - ,
    # # '100':,
    # # '010':
    # # '111': db_filtered_df_intrsct}
    # intersection_list = [db_total_count, df_pos_count, db_df_intrsct,
    #                      uniq_filtered_count, db_filtered_intrsct, filtered_df_intrsct,
    #                      db_filtered_df_intrsct]
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


def count_intersection(df_agg_intrsct, df_filtered, args):
    """function to get information on intersections of tables with databases"""
    # count instances in databases
    snp_total_count = os.popen("grep -v '#' {} | wc -l".format(args.snp_db_path)).read()  # count lines without headers
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
    run_venn(df_agg_intrsct, df_filtered, 'is_snp', snp_total_count, ['SNP', 'Aggregated data', 'Filtered data'],
             args.input_dir, args.sname)
    run_venn(df_agg_intrsct, df_filtered, 'is_editing_rep', edit_rep_total_count,
             ['Edit_rep', 'Aggregated data', 'Filtered data'], args.input_dir, args.sname)
    run_venn(df_agg_intrsct, df_filtered, 'is_editing_non_rep', edit_nonrep_total_count,
             ['Edit_non_rep', 'Aggregated data', 'Filtered data'], args.input_dir, args.sname)


def get_df(pathes):
    """function to load the df and filter it"""
    def get_filtered(df_agg, min_mutation_cb_to_filter, min_mutation_umis, min_total_umis):
        """function to return filtered aggregated table"""
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
        return df_agg_filt

    # extract paths
    df_intersection_path = pathes[3]

    # load and change the df with intersections notations to be binary (1 - if any overlap with db occured)
    df_agg_intrsct = pd.read_csv(df_intersection_path, sep='\t')
    df_agg_intrsct.loc[df_agg_intrsct['is_snp'] > 0, 'is_snp'] = 1
    df_agg_intrsct.loc[df_agg_intrsct['is_editing_rep'] > 0, 'is_editing_rep'] = 1
    df_agg_intrsct.loc[df_agg_intrsct['is_editing_non_rep'] > 0, 'is_editing_non_rep'] = 1
    df_agg_intrsct.to_csv(df_intersection_path, sep='\t', index=False)

    # get filtered df
    df_filtered = get_filtered(df_agg_intrsct, 5, 10, 20)

    return df_agg_intrsct, df_filtered


def find_intersections(args, pathes):
    # extract pathes
    agg_df_path = pathes[0]
    snp_temp = pathes[1]
    edit_temp = pathes[2]
    df_intersection = pathes[3]
    snp_db_path = args.snp_db_path
    rep_db_path = args.rep_db_path
    non_rep_db_path = args.non_rep_db_path

    # find intersection with help of bamtools intersect
    os.system("head -c 1 {agg_df} | grep -q '#' || sed -i '1s/^/#/' {agg_df}".format(agg_df = agg_df_path))  # add '#' to header of df_aggregated
    #### both files must be sorted if you use '-sorted' which reduce memoty usage
    os.system("bedtools intersect -a {agg_df} -b {snp_db_path} -c -header -sorted > {snp_temp}".format(
        agg_df=agg_df_path, snp_db_path=snp_db_path, snp_temp=snp_temp))  # find intersection with snp db
        
    os.system("sed -i '1 s/.*/&\tis_snp/' {snp_temp}".format(snp_temp=snp_temp))  # add column name 'is_snp'
    
    os.system("bedtools intersect -a {snp_temp} -b {rep_db_path} -c -header > {edit_temp}".format(snp_temp=snp_temp,rep_db_path=rep_db_path,edit_temp=edit_temp))  # find intersection with editing rep db
    
    os.system("sed -i '1 s/.*/&\tis_editing_rep/' {}".format(edit_temp)) # add column name 'is_editing_rep'
    
    os.system("bedtools intersect -a {edit_temp} -b {non_rep_db_path} -c -header > {df_intersection}".format(
        edit_temp=edit_temp,non_rep_db_path=non_rep_db_path,df_intersection=df_intersection))  # find intersection with editing non rep db
        
    os.system("sed -i '1 s/.*/&\tis_editing_non_rep/' {}".format(df_intersection))  # add column name 'is_editing_non_rep'
    
    os.system("rm {snp_temp} {edit_temp}".format(snp_temp=snp_temp,edit_temp=edit_temp))  #remove temp files


def make_pathes(input_dir):
    """make pathes to temporary and output files"""
    agg_df_path = os.path.join(input_dir, 'aggregated_tsv.tsv')
    snp_temp = os.path.join(input_dir ,'snp_intersect.tsv')
    edit_temp = os.path.join(input_dir ,'edit_intersect.tsv')
    df_intersection = os.path.join(input_dir ,'aggregated_intersect.tsv')
    return agg_df_path, snp_temp, edit_temp, df_intersection


if __name__ == '__main__':
    # parse arguments
    args = parse_arguments()

    # get paths to files we use
    pathes = make_pathes(args.input_dir)

    # find intersection between df and databases
    find_intersections(args, pathes)

    # get the df with intersections, before and after filtering
    df_agg_intersect, df_filtered = get_df(pathes)

    # make Venn diagrams of the intersections
    count_intersection(df_agg_intersect, df_filtered, args)

    print("Fininshed Venn diagrams")
