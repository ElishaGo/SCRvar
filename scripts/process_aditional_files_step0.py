#!/usr/bin/python
"""
This script process the DB we use in the pipeline: editing_DB from REDI portal, SNP vcf file and gtf genecode.
This script should be used once per file, because the output will be saved and used in the other steps.

zero and one base:
fasta file don't have numbers for position, thus depend on the tool reading it.
Bedtools commands uses 'start' with zero-base and 'end' with one-base.
for more information: https://bedtools.readthedocs.io/en/latest/content/overview.html?highlight=1-based#bed-starts-are-zero-based-and-bed-ends-are-one-based
vcf is one based
"""
import argparse
import os
import subprocess
from datetime import datetime
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)

import sc_rna_variants.analysis_utils
import sc_rna_variants.utils


class BadFileHeaderError(Exception):
    pass


def comment_out_header(f_path):
    """check if file header starts with '#'. otherwise, add it"""
    os.system(f"head -c 1 {f_path} | grep -q '#' || sed -i '1s/^/#/' {f_path}")


def load_and_process_fasta_coordinations(editing_DB_bed, fasta_path, editing_db):
    """
    find references from fasta and return as dataframe.
    We get the references from the fasta file by using bedtools getfasta which outputs file in the format:
    '>chr1:87157-87158(+)\n',
    't\n', ...
    """
    temp_ref_path = editing_db[:editing_db.rfind(os.sep)] + '/temp_ref_from_fasta_by_bed_file.txt'
    temp_editing_bed_path = editing_db[:editing_db.rfind(os.sep)] + '/editing_DB.bed'
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
    df = df.sort_values(by=['#chrom', 'chromStart'])

    # reorder columns
    # columns_reorder = ['#chrom', 'chromStart', 'chromEnd', 'Ref', 'score', 'strand', 'Ed', 'db', 'type',
    #                    'dbsnp', 'repeat',
    #                    'Func.wgEncodeGencodeBasicV34', 'Gene.wgEncodeGencodeBasicV34',
    #                    'GeneDetail.wgEncodeGencodeBasicV34', 'ExonicFunc.wgEncodeGencodeBasicV34',
    #                    'AAChange.wgEncodeGencodeBasicV34', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene',
    #                    'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene',
    #                    'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'phastConsElements100way']
    # df = df[columns_reorder]
    return df


def merge_DB_and_fasta_by_coor(editing_DB, fasta_coor):
    return pd.merge(editing_DB.drop_duplicates().reset_index(), fasta_coor.drop_duplicates(), how='inner').set_index(
        'index')


def split_editing_by_A_reference(editing_DB_bed, fasta_path, editing_db):
    print("Find only A base references")
    fasta_coordinates_df = load_and_process_fasta_coordinations(editing_DB_bed, fasta_path, editing_db)

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

    # add numeric column to keep bed format
    editing_DB_bed_df.insert(4, 'score', np.int8(0))

    return editing_DB_bed_df


def bedtools_intersect_u_flag(a_path, b_path, output_path):
    with open(output_path, "w") as outfile:
        subprocess.run(['bedtools', 'intersect', '-u', '-sorted', '-header', '-a', a_path, '-b', b_path], stdout=outfile)


def bedtools_sort(to_sort, f_sorted):
    temp = f_sorted + "temp"
    with open(temp, "w") as outfile:
        subprocess.run(['bedtools', 'sort', '-i', to_sort], stdout=outfile)

    # bedtools sort drops the header, so we add it manually
    with open(f_sorted, "w") as outfile:
        subprocess.run([f"grep '#' {to_sort} | cat - {temp}"], stdout=outfile, shell=True)
    os.remove(temp)


def add_chr_to_vcf(snp_db, f_with_chr):
    """add 'chr' to chromosome notations and sort file"""
    with open(f_with_chr, "w") as outfile:
        subprocess.run(['awk', '{if($0 !~ /^#/) print "chr"$0; else print $0}', snp_db], stdout=outfile)


def process_editing_DB(editing_db, output_dir, fasta_path, annotation_gtf):
    # TODO: what input should we expect?
    print("processing editing DB file")
    comment_out_header(editing_db)

    editing_DB_df = sc_rna_variants.analysis_utils.load_df(editing_db)

    editing_DB_bed = transform_to_bed(editing_DB_df)

    editing_A_I, editing_other = split_editing_by_A_reference(editing_DB_bed, fasta_path, editing_db)

    editing_A_I = sort_and_reorder(editing_A_I)
    editing_other = sort_and_reorder(editing_other)

    # save the editing DB file as bed file
    editing_A_I_path = os.path.join(output_dir, "0.editing_A_I.bed")
    editing_other_path = os.path.join(output_dir, "0.editing_other.bed")

    sc_rna_variants.analysis_utils.save_df(editing_A_I, output_dir, "0.editing_A_I.bed")
    sc_rna_variants.analysis_utils.save_df(editing_other, output_dir, "0.editing_other.bed")

    # intersect with gencode.gtf file
    bedtools_intersect_u_flag(a_path=editing_A_I_path, b_path=annotation_gtf,
                              output_path=os.path.join(output_dir, '0.editing_A_I.genecode_intersect.bed'))

    return editing_A_I_path, editing_other_path


def snp_DB_intersections(snp_db, out_dir, editing_db, annotation_gtf):
    """
    make intersection of the SNP DB with the entries in the editing files with reference base 'A', and make another
    intersection with the gencode (transcriptome) file
    :param snp_db:
    :param editing_db:
    :param annotation_gtf:
    :return:
    """
    # intersect with editing_A_I file
    snp_A_output = os.path.join(out_dir, '0.snp.A_I.vcf')
    bedtools_intersect_u_flag(a_path=snp_db,
                              b_path=os.path.join(editing_db[:editing_db.rfind(os.sep)], '0.editing_A_I.bed'),
                              output_path=snp_A_output)

    # instersect with gtf file
    snp_gtf_output = os.path.join(out_dir, '0.snp.gencode_intersect.vcf')
    bedtools_intersect_u_flag(a_path=snp_db, b_path=annotation_gtf, output_path=snp_gtf_output)

    # instersect with editing A_I and gtf file
    snp_A_gtf_output = os.path.join(out_dir, '0.snp.A_I.gencode_intersect.vcf')
    bedtools_intersect_u_flag(a_path=snp_A_output, b_path=annotation_gtf, output_path=snp_A_gtf_output)

    return snp_A_output, snp_gtf_output, snp_A_gtf_output


def process_gtf(orig_annotation_gtf, out_dir):
    print("processing genecode DB file")
    out_path = os.path.join(out_dir, "0." + os.path.basename(orig_annotation_gtf))
    bedtools_sort(orig_annotation_gtf, out_path)
    return out_path


def sort_and_replace(to_sort):
    f_sorted = to_sort.replace(os.path.basename(to_sort), "temp_" + os.path.basename(to_sort))
    bedtools_sort(to_sort, f_sorted)
    os.remove(to_sort)
    os.rename(f_sorted, f_sorted.replace('temp_', ''))


def process_snp(snp_db, snp_out_dir, editing_A_I_path, processed_annotation_gtf):
    print("processing SNP DB file")
    temp_snp_with_chr = snp_db.replace(os.path.basename(snp_db), "temp_" + os.path.basename(snp_db))
    add_chr_to_vcf(snp_db, temp_snp_with_chr)
    snp_A_output, snp_gtf_output, snp_A_gtf_output = snp_DB_intersections(temp_snp_with_chr, snp_out_dir,
                                                                          editing_A_I_path,
                                                                          processed_annotation_gtf)
    sort_and_replace(snp_A_output)
    sort_and_replace(snp_gtf_output)
    sort_and_replace(snp_A_gtf_output)

    os.remove(temp_snp_with_chr)


def run_step0(out_dir, annotation_gtf, editing_db, fasta_path, snp_db):
    # make output subdirectories
    out_dir = os.path.join(out_dir, 'data_files_processed')
    gtf_out_path = os.path.join(out_dir, "genecode_gtf")
    editing_out_dir = os.path.join(out_dir, "editing")
    snp_out_dir = os.path.join(out_dir, "snp_vcf")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(gtf_out_path, exist_ok=True)
    os.makedirs(editing_out_dir, exist_ok=True)
    os.makedirs(snp_out_dir, exist_ok=True)

    # process the files
    processed_annotation_gtf = process_gtf(annotation_gtf, gtf_out_path)
    if editing_db:
        editing_A_I_path, _ = process_editing_DB(editing_db, editing_out_dir, fasta_path, processed_annotation_gtf)
    if snp_db:
        process_snp(snp_db, snp_out_dir, editing_A_I_path, processed_annotation_gtf)


def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument('out-dir', type=sc_rna_variants.utils.assert_is_directory, help='folder for output files')
    parser.add_argument('annotation-gtf', type=sc_rna_variants.utils.assert_is_file, help='path to annotation gtf file')

    # optional arguments
    parser.add_argument('--editing-db', help='path to editing database file')
    parser.add_argument('--fasta-path', help='needed if editing_db is used. path to genome fasta file')
    parser.add_argument('--snp-db', help='path to SNP database file')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    startTime = datetime.now()
    args = parse_arguments()

    run_step0(args.out_dir, args.annotation_gtf, args.editing_db, args.fasta_path, args.snp_db)
    print(datetime.now() - startTime)
    print("\nfinished to process DB files")
