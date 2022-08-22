import argparse
import logging
import os

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)  # for development environments

import sc_rna_variants.utils
import sc_rna_variants.bio_functions
import sc_rna_variants.steps_runner


def parse_arguments(arguments=None):
    """argument parsing wrapper function
    helper functions and classes are found in sc_rna_variants.utils
    """
    parser = argparse.ArgumentParser(
        formatter_class=sc_rna_variants.utils.ArgparserFormater,
        description=""" The script filters out the reads with deletions/insersions/soft and hard clipped.
Reads aligned to genomic segments that don't appear in the genome FASTA are discarded.""",
        epilog='''Outputs a filtered BAM file'''
    )

    # positional arguments
    parser.add_argument('input_bam', type=sc_rna_variants.utils.bam_check,
                        help='the BAM file with scRNAseq data')
    parser.add_argument('output_folder', type=sc_rna_variants.utils.assert_is_directory,
                        help='Path to the directory where the output files will be placed')

    # optional arguments
    parser.add_argument('--barcodes-cluster-file',
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
    parser.add_argument('--log-file', default=os.path.join(sys.argv[2], '1.filter_bam.log'),
                        help='a log file for tracking the program\'s progress')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    args = parse_arguments()

    sc_rna_variants.config_logging(args.log_file)
    logger = logging.getLogger("sc_rna_variants")

    logger.info('scRNA variants finding program started')
    logger.debug('Running with parameters:\n%s' % '\n'.join(
        ['%s: %s' % (key, value) for key, value in vars(args).items() if key != 'barcodes_cluster_file'])
                 )

    sc_rna_variants.steps_runner.run_step1(args.input_bam, args.barcodes_cluster_file,
                                           args.min_mapq, args.cigar_clipping_allowed,
                                           args.max_gene_length, args.max_no_basecall,
                                           args.tag_for_umi, args.tag_for_cell_barcode,
                                           args.output_folder, args.threads)
