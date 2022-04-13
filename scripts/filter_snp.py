import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from joblib import load
from matplotlib_venn import venn2, venn2_circles

def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument('input_file', help='tsv file of mutated and unmutetd umis after running "intersection"')

    # positional arguments
    parser.add_argument('model', help='path to trained model to use for prediction. needs to be a joblib file')

    # optional arguments
    parser.add_argument('--output_folder', default=os.path.join(os.getcwd()),
                        help='Path to the output directory')

    parser.add_argument('--filter_known_snp', dest='filter_known_snp', action='store_true')
    parser.add_argument('--no-filter_known_snp', dest='filter_known_snp', action='store_false')
    parser.set_defaults(filter_known_snp=True)

    return parser.parse_args(arguments)


def preprocess(df):
    # preprocess data
    df_numeric_test = df.drop(['#chromosome', 'start', 'end', 'reference base', 'position',
                                    'strand', 'aggregated cell barcodes', 'is_snp', 'is_editing_rep',
                                    'is_editing_non_rep'], axis=1)
    # optianal columns to drop. may change in future
    # df_numeric_test = df_numeric_test.drop(
    #     ['percent of non ref mutated reads', 'unmutated single reads', 'count of unmutated cell barcodes',
    #      'mixed reads'], axis=1)
    # convert ref columns to dummies
    dummies = pd.get_dummies(df['reference base'], prefix='refbase')

    # concat scaled df and dummies columns
    df_numeric_test = pd.concat([df_numeric_test, dummies], join='inner', axis=1)
    return df_numeric_test


def get_filtered(df, min_mutation_cb_to_filter,min_mutation_umis, min_total_umis):
    """filtering function for aggregated tables"""
    mutation_umi_counts = df['total mutation umi count']
    total_umi_count = mutation_umi_counts + df['unmutated multi reads'] + df['unmutated single reads']

    cond_1 = (df['count of mutated cell barcodes'] >= min_mutation_cb_to_filter)
    cond_2 = ((mutation_umi_counts >= min_mutation_umis) & (total_umi_count >= min_total_umis))
    return df[cond_1 & cond_2]


def make_venn(num_pred_pos_snp, num_known_snp_in_table, num_overlaps_clf_db, args):
    venn2(subsets=(num_pred_pos_snp - num_overlaps_clf_db,
                   num_known_snp_in_table - num_overlaps_clf_db,
                   num_overlaps_clf_db),
          set_labels=('predicted SNP by model', 'known SNP in table'))
    plt.savefig(os.path.join(args.output_folder, 'SNP_predictions_venn_diagram.png'), facecolor='white')


if __name__ == '__main__':
    args = parse_arguments()  # parse argument
    df = pd.read_csv(args.input_file, sep= '\t')  # get dataset
    snp_clf = load(args.model)  # load the saved classification model

    # filter the table
    df = get_filtered(df=df, min_mutation_cb_to_filter=5, min_mutation_umis=10, min_total_umis=20)
    df_numeric_test = preprocess(df)  # preprocess the dataframe
    y_pred = snp_clf.predict(df_numeric_test)  # make prediction on data to find snp

    # add prediction to table
    df['predicted_as_snp'] = y_pred

    # get counts of different crosses
    num_known_snp_in_table = df[df['is_snp'] == 1].shape[0]
    num_overlaps_clf_db = df[(df['is_snp'] == 1) & (df['predicted_as_snp'] == 1)].shape[0]
    num_pred_pos_snp = df[df['predicted_as_snp'] == 1].shape[0]
    num_new_snp = df[(df['is_snp'] != 1) & (df['predicted_as_snp'] == 1)].shape[0]
    num_no_snp = df[(df['is_snp'] != 1) & (df['predicted_as_snp'] != 1)].shape[0]

    print("number of known SNPs from database in table:", num_known_snp_in_table)
    print("number of predicted SNPs by the classifier:", num_pred_pos_snp)
    print("number of overlaps between predicted SNP and known SNP:", num_overlaps_clf_db)
    print("number of new assumed SNP by the model:", num_new_snp)
    print("number of rows in table without any SNP:", num_no_snp)

    # make venn diagram
    make_venn(num_pred_pos_snp, num_known_snp_in_table, num_overlaps_clf_db, args)

    # save the table
    if args.filter_known_snp:
        df = df[df['is_snp'] != 1]
        df.to_csv(os.path.join(args.output_folder,'df_with_predicted_snp_excluding_known.tsv'), sep='\t',index=False)
    else:
        df.to_csv(os.path.join(args.output_folder, 'df_with_predicted_snp.tsv'), sep='\t', index=False)