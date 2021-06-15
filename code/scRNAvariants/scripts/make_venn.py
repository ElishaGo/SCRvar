import pandas as pd
import argparse
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import os
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
    
    return parser.parse_args(arguments)


def make_venn_diagram(df, subset_list, labels, column_name,input_dir):
    fig, axs = plt.subplots(2, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle('Intersection - {} and {} {}'.format(labels[0], labels[2], labels[1]),
                 fontsize=18, bbox={"facecolor": "orange", "alpha": 0.5})
    axs[0].set_title("Venn diagram \n{a}: {b},   {c}: {d},   {e}: {f} \n{g}: {h},  {i}: {j},  {k}: {l}".format(
        a=labels[0], b=subset_list[0],
        c=labels[1], d=subset_list[1],
        e="intersection", f=subset_list[2],
        g=labels[2], h=subset_list[5] + subset_list[3],
        i=labels[0] + " & " + labels[2], j=subset_list[4],
        k="all_intersect", l=subset_list[6]), fontsize=10, bbox={"facecolor": "gray", "alpha": 0.3})

    v = venn3(subsets=subset_list, set_labels=(labels[0], labels[1], labels[2]), ax=axs[0])
    venn3_circles(subsets=subset_list, color='gray', linewidth=1, linestyle='dashed', ax=axs[0])

    # make histogram of mutated CB
    intrsct_list = df[df[column_name] == 1]['count of mutated cell barcodes']
    axs[1].hist(intrsct_list, bins=30)
    axs[1].hist([i for i in intrsct_list if i >= 5], bins=30, alpha=0.7)
    axs[1].set_title("Count of mutated cell barcodes in intersection with {}".format(labels[0]), fontsize=10,
                     bbox={"facecolor": "gray", "alpha": 0.3})
    axs[1].set_xlabel("mutations")
    axs[1].set_ylabel("cell barcodes")
    axs[1].legend(['total CB counts', 'CB counts >= 5'])
    fig.tight_layout()
    plt.savefig(os.path.join(input_dir, 'venn_diagram_{}.png'.format(labels[0])), facecolor='white')


def run_venn(df,df_filtered, column_name, db_total_count, labels,input_dir):
    df_pos_count = df.shape[0]
    filtered_df_intrsct = len(df.index & df_filtered.index)
    filtered_count = df_filtered.shape[0] - filtered_df_intrsct
    db_df_intrsct = df[column_name].sum()
    db_filtered_intrsct = df_filtered[column_name].sum()
    db_filtered_df_intrsct = df_filtered[column_name].sum()
    intersection_list = [db_total_count, df_pos_count, db_df_intrsct,
                         filtered_count, db_filtered_intrsct, filtered_df_intrsct,
                         db_filtered_df_intrsct]
    make_venn_diagram(df, intersection_list, labels, column_name,input_dir)


def make_diagrams(args,pathes):
    def get_filtered(df_agg, min_mutation_cb_to_filter, min_mutation_umis, min_total_umis):
        """function to return filtered aggregated table"""
        cond_1 = (df_agg['count of mutated cell barcodes'] >= min_mutation_cb_to_filter)

        mutation_umi_counts = df_agg['total umi counts 1 cells'] + \
                              df_agg['total umi counts 2 cells'] + \
                              df_agg['total umi counts 3 cells'] + \
                              df_agg['total umi counts 4+ cells']
        total_umi_count = mutation_umi_counts + \
                          df_agg['unmutated multi reads'] + \
                          df_agg['unmutated single reads']
        cood_2 = ((mutation_umi_counts >= min_mutation_umis) & (total_umi_count >= min_total_umis))

        df_agg_filt = df_agg[cond_1 & cood_2]
        return  df_agg_filt

    # extract pathes
    df_intersection = pathes[3]
    input_dir = args.input_dir
    snp_db_path = args.snp_db_path
    rep_db_path = args.rep_db_path
    non_rep_db_path = args.non_rep_db_path
    
    # count lines in databases files
    snp_total_count = os.popen("grep -v '#' {} | wc -l".format(snp_db_path)).read()  # count lines without headers
    edit_rep_total_count = os.popen("cat {} | wc -l".format(rep_db_path)).read()
    edit_nonrep_total_count = os.popen("cat {} | wc -l".format(non_rep_db_path)).read()
    snp_total_count = int(snp_total_count)
    edit_rep_total_count = int(edit_rep_total_count)
    edit_nonrep_total_count = int(edit_nonrep_total_count)

    # load and change the df with intersections notations to be binary (1 - if any overlap with db occured)
    df_agg_intrsct = pd.read_csv(df_intersection, sep='\t')
    df_agg_intrsct.loc[df_agg_intrsct['is_snp'] > 0, 'is_snp'] = 1
    df_agg_intrsct.loc[df_agg_intrsct['is_editing_rep'] > 0, 'is_editing_rep'] = 1
    df_agg_intrsct.loc[df_agg_intrsct['is_editing_non_rep'] > 0, 'is_editing_non_rep'] = 1
    df_filtered = get_filtered(df_agg_intrsct, 5,10, 20)

    # make Venn diagrams for snp, editing rep and editing non_rep intersections
    run_venn(df_agg_intrsct,df_filtered, 'is_snp', snp_total_count, ['SNP', 'Aggregated data', 'Filtered data'],input_dir)
    run_venn(df_agg_intrsct,df_filtered, 'is_editing_rep', edit_rep_total_count, ['Edit_rep', 'Aggregated data', 'Filtered data'],input_dir)
    run_venn(df_agg_intrsct,df_filtered, 'is_editing_non_rep', edit_nonrep_total_count,['Edit_non_rep', 'Aggregated data', 'Filtered data'],input_dir)


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
    pathes = make_pathes(args.input_dir)
    find_intersections(args, pathes)
    make_diagrams(args, pathes)
    print("Fininshed Venn diagrams")
