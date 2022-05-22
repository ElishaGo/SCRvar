"""
This scripts processes editing DB file TABLE_hg38.txt downloaded from REDIportal,
and keeps a filtered file with coordinates of only 'a' base as reference.

zero and one base:
fasta file don't have numbers for position, thus depend on the tool reading it.
Bedtools commands uses 'start' with zero-base and 'end' with one-base.
for more information: https://bedtools.readthedocs.io/en/latest/content/overview.html?highlight=1-based#bed-starts-are-zero-based-and-bed-ends-are-one-based
"""
import argparse
import os
import subprocess
from datetime import datetime
import pandas as pd
import numpy as np
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent.parent.absolute()) + os.path.sep)

import sc_rna_variants.analysis_utils
import sc_rna_variants.utils


class BadFileHeaderError(Exception):
    pass


def check_editing_DB_header(editing_file_path):
    with open(editing_file_path) as f:
        first_line = f.readline()
    good_header = first_line[0] == '#'
    return good_header


def load_and_process_fasta_coordinations(editing_DB_bed, fasta_path, editing_DB_path):
    """
    find referenes from fasta and return as dataframe.
    We get the references from the fasta file by using bedtools getfasta which outputs file in the format:
    '>chr1:87157-87158(+)\n',
    't\n', ...
    """
    temp_ref_path = editing_DB_path[:editing_DB_path.rfind(os.sep)] + '/temp_ref_from_fasta_by_bed_file.txt'
    temp_editing_bed_path = editing_DB_path[:editing_DB_path.rfind(os.sep)] + '/editing_DB.bed'
    if not os.path.isfile(temp_ref_path):
        print("running bedtools getfasta")
        editing_DB_bed.to_csv(temp_editing_bed_path, index=False, sep='\t')
        os.system(f"bedtools getfasta -s -fi {fasta_path} -bed {temp_editing_bed_path} > {temp_ref_path}")

    with open(temp_ref_path, 'r') as f:
        coordinates = f.readlines()  # into list

    # process the coordinates to be in one line (e.g 'chr1:87157-87158(+) t')
    coordinates = [''.join(x).replace('\n', '') for x in zip(coordinates[0::2], coordinates[1::2])]

    # parse the coordinates
    coordinates_df = pd.DataFrame()
    coordinates_df[['#chrom', 'chromStart']] = pd.Series(coordinates).str.split(':', expand=True)
    coordinates_df[['chromStart', 'strand']] = coordinates_df['chromStart'].str.split('(', expand=True)
    coordinates_df[['strand', 'Ref']] = coordinates_df['strand'].str.split(')', expand=True)
    coordinates_df[['chromStart', 'chromEnd']] = coordinates_df['chromStart'].str.split('-', expand=True)

    coordinates_df['#chrom'] = coordinates_df['#chrom'].str.replace('>', '')
    coordinates_df[['chromStart', 'chromEnd']] = coordinates_df[['chromStart', 'chromEnd']].astype('int')
    coordinates_df['Ref'] = coordinates_df['Ref'].str.lower()
    return coordinates_df


def sort_and_reorder(df):
    print("start sort")
    df = df.sort_values(by=['#chrom', 'chromStart'])

    # reorder columns
    columns_reorder = ['#chrom', 'chromStart', 'chromEnd', 'Ref', 'score', 'strand', 'Ed', 'db', 'type',
                       'dbsnp', 'repeat',
                       'Func.wgEncodeGencodeBasicV34', 'Gene.wgEncodeGencodeBasicV34',
                       'GeneDetail.wgEncodeGencodeBasicV34', 'ExonicFunc.wgEncodeGencodeBasicV34',
                       'AAChange.wgEncodeGencodeBasicV34', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene',
                       'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene',
                       'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'phastConsElements100way']
    return df[columns_reorder]


def merge_DB_and_fasta_by_coor(editing_DB, fasta_coor):
    return pd.merge(editing_DB.drop_duplicates().reset_index(), fasta_coor.drop_duplicates(), how='inner').set_index(
        'index')


def split_editing_by_A_reference(editing_DB_bed, fasta_path, editing_DB_path):
    fasta_coordinates_df = load_and_process_fasta_coordinations(editing_DB_bed, fasta_path, editing_DB_path)

    # TODO: ask Dena about filter by strand. currently we ignore here strandness
    fasta_A_I = fasta_coordinates_df[fasta_coordinates_df['Ref'] == 'a']
    fasta_other = fasta_coordinates_df[fasta_coordinates_df['Ref'] != 'a']

    editing_other = merge_DB_and_fasta_by_coor(editing_DB_bed.drop(['Ref', 'strand'], axis=1), fasta_other)
    editing_A_I = merge_DB_and_fasta_by_coor(editing_DB_bed.drop(['Ref', 'strand'], axis=1), fasta_A_I)

    return editing_A_I, editing_other


def transform_to_bed(editing_DB_df):
    """
    Since the file is downloaded from REDIportal, and has only one coordinate, it is one-bases. Thus we will use the
    given coordinate as the 'chromEnd' position when we convert the file to bed format, and we add manually the
    'chromStart' column
    :param editing_DB_df:
    :return:
    """
    # add start coordinate
    editing_DB_df['chromStart'] = editing_DB_df['chromEnd'] - 1

    # reorder columns
    columns_reorder = ['#chrom', 'chromStart', 'chromEnd', 'Ref', 'strand', 'Ed', 'db', 'type', 'dbsnp', 'repeat',
                       'Func.wgEncodeGencodeBasicV34', 'Gene.wgEncodeGencodeBasicV34',
                       'GeneDetail.wgEncodeGencodeBasicV34', 'ExonicFunc.wgEncodeGencodeBasicV34',
                       'AAChange.wgEncodeGencodeBasicV34', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene',
                       'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene',
                       'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'phastConsElements100way']
    editing_DB_bed_df = editing_DB_df[columns_reorder]

    editing_DB_bed_df.insert(4, 'score', np.int8(0))  # add numeric column to keep bed format

    return editing_DB_bed_df


def intersect_with_gtf(a_path, gtf_path, output_path):
    subprocess.run(['bedtools', ' intersect', '-u', '-header', '-a', a_path, '-b', gtf_path, '>', output_path])


def process_editing_DB(editing_DB_path, fasta_path, gtf_path):
    # TODO: check how to change the header on the fly
    good_header = check_editing_DB_header(editing_DB_path)
    if not good_header:
        raise BadFileHeaderError(f"Need to add '#' to the beginning of the header in file {editing_DB_path}")

    editing_DB_df = sc_rna_variants.analysis_utils.load_df(editing_DB_path)

    editing_DB_bed = transform_to_bed(editing_DB_df)

    editing_A_I, editing_other = split_editing_by_A_reference(editing_DB_bed, fasta_path,
                                                              editing_DB_path)

    editing_A_I = sort_and_reorder(editing_A_I)
    editing_other = sort_and_reorder(editing_other)

    # save the editing DB file as bed file
    output_dir = editing_DB_path[:editing_DB_path.rfind(os.sep)]
    editing_A_I_path = os.path.join(output_dir, "0_editing_A_I.bed6")

    sc_rna_variants.analysis_utils.save_df(editing_A_I, output_dir, "0_editing_A_I.bed6")
    sc_rna_variants.analysis_utils.save_df(editing_other, output_dir, "0_editing_other.bed6")

    # intersect with gencode.gtf file
    intersect_with_gtf(a_path=editing_A_I_path, gtf_path=gtf_path,
                       output_path=os.path.join(output_dir, '0_editing_A_I.genecode_intersect.bed6'))


def snp_DB_intersections(snp_DB_path, editing_DB_path, gtf_path):
    """
    make intersection of the SNP DB with the entries in the editing files with reference base 'A', and make another
    intersection with the gencode (transcriptome) file
    :param snp_DB_path:
    :param editing_DB_path:
    :param gtf_path:
    :return:
    """
    # intersect with editing_A_I file
    snp_A_output = os.path.join(snp_DB_path[:snp_DB_path.rfind(os.sep)], '0_snp_A.vcf')
    subprocess.run(['bedtools', 'intersect', '-u', '-header', '-a', snp_DB_path, '-b',
                    os.path.join(editing_DB_path, '0_editing_A_I.bed6'), '>', snp_A_output])

    # instersect with gtf file
    snp_A_gtf_output = os.path.join(snp_DB_path[:snp_DB_path.rfind(os.sep)], '0_snp_A.gencode_intersect.vcf')
    subprocess.run(
        ['bedtools', 'intersect', '-u', '-header', '-a', snp_A_output, '-b', gtf_path, '>', snp_A_gtf_output])


def run_step0(args):
    process_editing_DB(args.editing_DB_path, args.fasta_path, args.gtf_path)
    snp_DB_intersections(args.snp_DB_path, args.editing_DB_path, args.gtf_path)


def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser()

    parser.add_argument('editing_DB_path')
    parser.add_argument('fasta_path')
    parser.add_argument('gtf_path')
    parser.add_argument('snp_DB_path')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    startTime = datetime.now()

    args = parse_arguments()
    run_step0(args)
    print(datetime.now() - startTime)
    print("\nfinished to process DB files")
