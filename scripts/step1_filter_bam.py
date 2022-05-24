import argparse
import logging

import sys, os  # for development environments
from pathlib import Path
import subprocess

sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)  # for development environments

import sc_rna_variants.utils
import sc_rna_variants.bio_functions


def parse_arguments(arguments=None):
    """argument parsing wrapper function
    helper functions and classes are found in sc_rna_variants.utils
    """
    parser = argparse.ArgumentParser(
        formatter_class=sc_rna_variants.utils.ArgparserFormater,
        description="""A script to locate cases of RNA modifications in single cell RNAseq data

By default the umi and barcode cells are comptible to bam files from cellranger. 
For other formats you need to change the parameters of cell barcodes and UMI tags.
'  

The script filters out the reads with deletions/insersions/soft and hard clipped.
Reads aligned to genomic segments that don't appear in the genome FASTA are discarded.""",
        epilog='''Outputs a filtered BAM file, and a BED formated file
where each row represents the modifications a certain cell has in a certain position'''
    )

    # positional arguments
    parser.add_argument('input_bam', type=sc_rna_variants.utils.bam_check,
                        help='the BAM file with scRNAseq data')
    parser.add_argument('output_folder', type=sc_rna_variants.utils.assert_is_directory,
                        help='Path to the directory where the output files will be placed')

    # optional arguments
    parser.add_argument('--filtered-barcodes-list',
                        type=sc_rna_variants.utils.filtered_barcodes_processing,
                        # returns a set with the barcodes names
                        help='''Text/tsv file with a list of cell barcodes as first column. Counts only these cells. Please note GEM-well numbers are ignored''')
    parser.add_argument('--min-mapq', help='Minimum quality of the read mapping',
                        type=int, default=255)
    parser.add_argument('--cigar-clipping-allowed', type=int, default=2,
                        help='how many clipped positions are allowed in a read\'s CIGAR')
    parser.add_argument('--max-gene-length',
                        help='Maximum length of the gene. Reads that will be mapped to with longer gaps will be discarded',
                        type=int, default=100000)
    parser.add_argument('--max-no-basecall',
                        help='Maximum number of N bases allowed in a single read',
                        type=int, default=2)
    parser.add_argument('--tag-for-cell-barcode',
                        help='the corrected cell barcode tag in the bam file. Reads without this tag are filtered out',
                        type=str, default='CB')
    parser.add_argument('--tag-for-umi',
                        help='the corrected umi tag in the bam file. Reads without this tag are filtered out',
                        type=str, default='UB')

    # Meta arguments
    parser.add_argument('--threads',
                        help='number of available threads',
                        type=int, default=1)
    # To Do change the path to be the same as output folder
    parser.add_argument('--log-file', default=os.path.join(sys.argv[2], '1.filter_bam.log'),
                        help='a log file for tracking the program\'s progress')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    args = parse_arguments()

    sc_rna_variants.config_logging(args.log_file)
    logger = logging.getLogger("sc_rna_variants")

    logger.info('scRNA variants finding program started')
    logger.debug('Running with parameters:\n%s' % '\n'.join(
        ['%s: %s' % (key, value) for key, value in vars(args).items() if key != 'filtered_barcodes_list'])
                 )

    # create the filtered bam from which the variants will be counted
    filtered_bam_path = sc_rna_variants.bio_functions.create_filtered_bam(args.input_bam, args.filtered_barcodes_list,
                                                      args.min_mapq, args.cigar_clipping_allowed,
                                                      args.max_gene_length, args.max_no_basecall,
                                                      args.tag_for_umi, args.tag_for_cell_barcode,
                                                      args.output_folder, args.threads)

    # add chr to chromosome names in bam files
    subprocess.run([f"samtools view -H {filtered_bam_path} | sed -e '/SN:chr/!s/SN:\([0-9XY]*\)/SN:chr&/' -e '/SN:chrM/!s/SN:MT/SN:chrM&/' | samtools reheader - {filtered_bam_path} > {filtered_bam_path}_temp"],
                   shell=True)
    os.remove(filtered_bam_path)
    os.rename(f"{filtered_bam_path}_temp", filtered_bam_path)
    subprocess.run(['samtools', 'index', filtered_bam_path])
