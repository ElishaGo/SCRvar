"""
This script creates a text file in bsub format to run the pipeline on LSF system (like WEXAC),
and execute it as a bsub command
"""
import os
import sys
import argparse
import logging
from datetime import datetime

from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent.parent.absolute()) + os.path.sep)  # for development environments

from sc_rna_variants.utils import assert_is_directory, assert_is_file
from sc_rna_variants import config_logging
from sc_rna_variants.general_utils import save_how_to

from scripts.step5_filtering_positions_and_snp_editing_DB_intersections import run_step5
from scripts.step6_gene_level import run_step6


def makedir(dir_path):
    os.makedirs(dir_path, exist_ok=True)
    return dir_path


def run_steps_5_6(args):
    # run step 5
    run_step5(args.input_dir, args.output_dir, args.min_cb_per_pos, args.min_mutation_umis, args.min_total_umis,
              args.min_mutation_rate, args.snp_db_path, args.editing_db_path, args.atacseq_path, args.sname)

    # run step 6
    run_step6(args.input_dir, args.output_dir, args.read_per_barcode_raw_bam, args.min_cb_per_pos,
              args.min_mutation_umis, args.min_total_umis, args.min_mutation_rate, args.atacseq_path, args.gtf_path,
              args.mismatch_dict_bed, args.barcode_clusters, args.atacseq_gcoverage_min, args.atacseq_gfrequency_min,
              args.sname)


#########################################################################################################
def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(description="""A script to set parameter to a bsub file and send to bsub""", )

    # positional arguments
    parser.add_argument('input_dir', type=assert_is_directory, help='step 4 output folder')
    parser.add_argument('output_dir', type=makedir, help='folder for outputs')
    parser.add_argument('snp_db_path', type=assert_is_file, help='path to known SNP sites file')
    parser.add_argument('editing_db_path', type=assert_is_file, help='path to known editing sites file')
    parser.add_argument('gtf_path', type=assert_is_file, help='path to gtf file')
    parser.add_argument('mismatch_dict_bed', type=assert_is_file,
                        help='path to 3.mismatch_dictionary.bed')
    parser.add_argument('read_per_barcode_raw_bam', type=assert_is_file,
                        help='count of reads per cell barcode in raw bam file')

    # optional arguments
    parser.add_argument('--min_cb_per_pos', default=5, type=int,
                        help='position with less cell barcodes will be filtered')
    parser.add_argument('--min_mutation_umis', default=10, type=int,
                        help='position with less mutated UMIs will be filtered')
    parser.add_argument('--min_total_umis', default=20, type=int,
                        help='position with less number of mutated + unmutated UMIs will be filtered')
    parser.add_argument('--min_mutation_rate', default=0.1, type=int,
                        help='position with less rate of mutation will be filtered')
    parser.add_argument('--atacseq_path', type=str, help='path to atacseq file')
    parser.add_argument('--atacseq_gcoverage_min', type=int, default=5)
    parser.add_argument('--atacseq_gfrequency_min', type=float, default=0.2)
    parser.add_argument('--barcode_clusters', type=assert_is_file,
                        help='table with barcodes and associated clusters analysed by Seurat')

    # Meta arguments
    parser.add_argument('--log-file',
                        default=os.path.join(sys.argv[2], "steps_5_6_customized.log"),
                        help='a log file for tracking the program\'s progress')
    parser.add_argument('--sname', default='', help='sample name to add to outputs')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    startTime = datetime.now()
    args = parse_arguments()

    # initialize logger
    config_logging(args.log_file)
    logger = logging.getLogger("positions_filtering_and_plots")
    logger.info('positions_filtering_and_plots started')
    logger.debug('Running with parameters:\n%s' % '\n'.join(
        ['%s: %s' % (key, value) for key, value in vars(args).items()]))

    run_steps_5_6(args)

    print(datetime.now() - startTime)
    print("finished steps 5 and 6")
