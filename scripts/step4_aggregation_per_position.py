import argparse
import os
import logging
from datetime import datetime

import sys  # for development environments
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)  # for development environments

import sc_rna_variants.steps_runner
import sc_rna_variants.utils

logger = logging.getLogger(__name__)


def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(sc_rna_variants.utils.ArgparserFormater)

    # positional arguments
    parser.add_argument('input_dir', type=sc_rna_variants.utils.assert_is_directory,
                        help='folder with mismatch_dictionaries (step 3)')
    parser.add_argument('output_dir', help='folder for step4 outputs', type=sc_rna_variants.utils.assert_is_directory)
    parser.add_argument('editing_db_path', type=sc_rna_variants.utils.assert_is_file,
                        help='path to known editing sites file')
    parser.add_argument('snp_db_path', type=sc_rna_variants.utils.assert_is_file, help='path to known SNP sites file')
    parser.add_argument('annotation_gtf', type=sc_rna_variants.utils.assert_is_file,
                        help='path to transcriptome gtf file')

    # optional arguments
    parser.add_argument('--barcodes-clusters',
                        type=sc_rna_variants.utils.filtered_barcodes_processing,
                        # returns a set with the barcodes names
                        help='''Text/tsv file with a list of cell barcodes as first column. Counts only these cells. Please note GEM-well numbers are ignored''')
    parser.add_argument('--min_cb_per_pos', default=5, type=int,
                        help='position with less cell barcodes will be filtered')
    parser.add_argument('--min_mutation_umis', default=10, type=int,
                        help='position with less mutated UMIs will be filtered')
    parser.add_argument('--min_total_umis', default=20, type=int,
                        help='position with less number of mutated + unmutated UMIs will be filtered')
    parser.add_argument('--min_mutation_rate', default=0.1, type=int,
                        help='position with less rate of mutation will be filtered')

    # Meta arguments
    parser.add_argument('--log-file',
                        default=os.path.join(sys.argv[2], '4.aggregated_per_position_and_statisitcs.log'),
                        help='a log file for tracking the program\'s progress')
    parser.add_argument('--sname', type=str, help='sample name to add to outputs')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    startTime = datetime.now()
    args = parse_arguments()

    # initialize logger
    sc_rna_variants.config_logging(args.log_file)
    logger = logging.getLogger("aggregation_per_position")
    logger.info('Aggregation per position started')
    logger.debug('Running with parameters:\n%s' % '\n'.join(
        ['%s: %s' % (key, value) for key, value in vars(args).items()]))

    # run statistics analysis
    sc_rna_variants.steps_runner.run_step4(args.input_dir, args.output_dir, args.annotation_gtf, args.barcodes_clusters,
                                           args.editing_db_path, args.snp_db_path)

    print(datetime.now() - startTime)
    logger.info('Step 4 finished')
    print("Step 4 finished")
