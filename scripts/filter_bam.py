import os
import argparse
import pysam
from subprocess import Popen, PIPE
import argparse
import logging

import sys, os  # for development environments
from pathlib import Path

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
    # parser.add_argument('genome_fasta', type=sc_rna_variants.utils.genome_fasta_check,
    #                     help='a FASTA file of the genome to which the BAM reads are aligned')
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
    # parser.add_argument('--tag-for-sample-barcode',
    #                     help='the sample barcode tag in the bam file',
    #                     type=str, default='BC')  # TBD: Not sure we need this
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
    parser.add_argument('--log-file', default=None,
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
    sc_rna_variants.bio_functions.create_filtered_bam(args.input_bam, args.filtered_barcodes_list,
                                                      args.min_mapq, args.cigar_clipping_allowed,
                                                      args.max_gene_length, args.max_no_basecall,
                                                      args.tag_for_umi, args.tag_for_cell_barcode,
                                                      args.output_folder, args.threads)

# def parse_arguments(arguments=None):
#     parser = argparse.ArgumentParser()
#     # positional arguments
#     parser.add_argument('input_bam', help='bam file to filter.')
#
#     parser.add_argument('filter_list', help='list of cell barcodes to keep in the output bam file')
#
#     parser.add_argument('--name_suffix', type=str, default="", help='Add suffix to output files')
#
#     parser.add_argument('--threads', type=int, default=1)
#
#     return parser.parse_args(arguments)
#
#
# def make_output_folder():
#     """Create a folder to store the filtered bam file"""
#     output_path = os.path.join(os.getcwd(), 'filtered_bam_files')
#     if not os.path.isdir(output_path):
#         print("Creating folder for the filtered bam file")
#         os.mkdir(output_path)
#
#     return output_path
#
#
# def process_filter_list(filter_list_path):
#     """Check for duplicates and optional problematic elements in list"""
#     # load filter list with end signs to object
#     filter_list = os.popen("cat -E {}".format(filter_list_path)).read().split('\n')
#     # filter_list = Popen(["cat", "-E", str(filter_list_path)]).read().split('\n')
#     # remove duplicates
#     filter_list = list(set(filter_list))
#
#     for cb in filter_list:
#         if cb in ['$', '', '\t', '\n']:
#             filter_list.remove(cb)
#     # check if end sign is by it's own, and correct
#     for i, cb in enumerate(filter_list):
#         # remove space from barcodes
#         filter_list[i] = cb.replace(' ', '').replace('$','')
#         # clear empty lines
#         # if cb in ['$', '', '\t', '\n']:
#         #     filter_list.pop(i)
#         #     i =- 1
#         # make sure the barcode are written like in the bam file
#         if not cb.startswith("CB:Z:"):
#             filter_list[i] = "CB:Z:" + filter_list[i]
#
#     # save the clean list
#     clean_list_path = filter_list_path.split("/")
#     clean_list_path[-1] = "temp_" + clean_list_path[-1]
#     clean_list_path = "/".join(clean_list_path)
#     with open(clean_list_path, "w") as f:
#         f.writelines("%s\n" % cb for cb in filter_list)
#
#     return clean_list_path
#
#
# def filter_bam(bam_file, out_folder, args, filter_list):
#     '''Recieve bam file and filter list, and creates a new bam file with only the cell barcodes from the filter list'''
#     suffix = args.name_suffix
#     threads = args.threads
#
#     # os.system("export name={suffix};export BAM_FILE={bam_file};export FILTER_FILE={filter_list}; samtools view -H $BAM_FILE > {out}/SAM_header_$name; samtools view -@ {threads} $BAM_FILE | LC_ALL=C grep -F -f $FILTER_FILE | cat {out}/SAM_header_$name - | samtools view -@ {threads} -b - > {out}/CBfiltered_$name.bam;samtools index {out}/CBfiltered_$name.bam;rm SAM_header_$name;rm $FILTER_FILE".format(bam_file=bam_file, filter_list=filter_list, out= out_folder,suffix=suffix, threads = threads))
#
#     os.system("samtools view -H {bam_file} > {out}/{name}_SAM_header; samtools view -@ {threads} {bam_file} |"
#               " LC_ALL=C grep -F -f {filter_list} | cat {out}/{name}_SAM_header - | "
#               "samtools view -@ {threads} -b - > {out}/{name}_CBfiltered.bam;samtools index {out}/{name}_CBfiltered.bam;"
#               "".format(bam_file=bam_file, filter_list=filter_list,
#                         out=out_folder, name=suffix, threads=threads))
#
#
# def open_bam(input_bam, available_threads=1):
#     """
#     tries to open a bam file, and asserts it is sorted by coordinates
#     returns pysam.libcalignmentfile.AlignmentFile
#     """
#     try:
#         bamfile = pysam.AlignmentFile(input_bam, "rb", threads=available_threads)
#
#         if bamfile.header['HD']['SO'] != "coordinate":
#             msg = "BAM is not sorted by coordinate (missing SO:coordinate from header)"
#             raise AssertionError(msg)
#         return bamfile
#     except Exception as e:
#         msg = "\nError handeling the input bam file %s.\n%s" % (input_bam, str(e))
#         raise AssertionError(msg)
#
#
# def run(args):
#     # create output folder
#     output_path = make_output_folder()
#     # read bam file
#     bamfile = open_bam(args.input_bam)
#     # load and process the filter list
#     filter_list = process_filter_list(args.filter_list)
#     # filter the bam file by filter list
#     filter_bam(bamfile, output_path, args, filter_list)
#
#
# if __name__ == '__main__':
#     os.system("echo Starting program to filter bam file by supplied cell barcodes")
#     args = parse_arguments()
#     run(args)
#     os.system("echo finished to filter bam file by cell barcodes")
