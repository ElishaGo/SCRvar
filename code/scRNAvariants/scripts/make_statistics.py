#!/usr/bin/env python
import argparse
import logging
from datetime import datetime

import sys, os  # for development environments
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)  # for development environments

import sc_rna_variants.utils
import sc_rna_variants.stats_generator

logging.getLogger('matplotlib').setLevel(logging.CRITICAL)


def parse_arguments(arguments=None):
    """argument parsing wrapper function
    helper functions and classes are found in sc_rna_variants.utils

    # TODO explain more the 'epilog
    """
    parser = argparse.ArgumentParser(
        formatter_class=sc_rna_variants.utils.ArgparserFormater,
        description="""This script aggregates and learns statistics on the output files from 'scrnavariants.py'.""",
        epilog='''Outputs aggregated tsv, figures and statistics.'''
    )

    # positional arguments
    parser.add_argument('input_mutated_file', type=sc_rna_variants.utils.assert_is_file,
                        help='raw_stats.tsv file of mutated umis as it came out from scrnavariants.py')

    parser.add_argument('input_unmutated_file', type=sc_rna_variants.utils.assert_is_file,
                        help='raw_umutated_stats.tsv file of un-mutated umis as it came out from scrnavariants.py')

    # optional arguments
    parser.add_argument('--output_folder', type=sc_rna_variants.utils.assert_is_directory,
                        default=os.path.join(os.getcwd(), 'statistics_output'),
                        help='Path to the directory where the output files will be placed')

    parser.add_argument('--min_cb_per_pos', default=5, type=int,
                        help='position with less cell barcodes will be filtered')

    parser.add_argument('--min_mutation_umis', default=10, type=int,
                        help='position with less mutated UMIs will be filtered')

    parser.add_argument('--min_total_umis', default=20, type=int,
                        help='position with less number of mutated + unmutated UMIs will be filtered')

    # Meta arguments
    parser.add_argument('--threads', type=int,
                        help='number of available threads', default=1)
    parser.add_argument('--log-file', default=None,
                        help='a log file for tracking the program\'s progress')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    # parse arguments
    args = parse_arguments()

    # initialize logger
    sc_rna_variants.config_logging(args.log_file)
    logger = logging.getLogger("statistics_analysis")
    logger.info('scRNA variants statistics analysis started')
    logger.debug('Running with parameters:\n%s' % '\n'.join(
        ['%s: %s' % (key, value) for key, value in vars(args).items()]))

    startTime = datetime.now()
    # run statistics analysis
    # TODO: replace the 'arguments' variable with explicit arguments.
    sc_rna_variants.stats_generator.run(args)
    print(datetime.now() - startTime)

    logger.info('Program finished')
    logger.info('Time took for analysis: %s' % str(datetime.now() - startTime))
