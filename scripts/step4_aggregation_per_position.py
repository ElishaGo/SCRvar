#!/usr/bin/env python
import argparse
import logging
from datetime import datetime

import sys, os  # for development environments
from pathlib import Path

# sys.path.append(str(Path(__file__).parent.parent.absolute()) + os.path.sep)  # for development environments

from sc_rna_variants.utils import assert_is_directory, ArgparserFormater
import sc_rna_variants.stats_generator

logging.getLogger('matplotlib').setLevel(logging.CRITICAL)


def make_output_dir(string):
    assert_is_directory(string)
    output_path = os.path.join(string, "step4_outputs")
    os.makedirs(output_path, exist_ok=True)
    return output_path


def parse_arguments(arguments=None):
    """argument parsing wrapper function
    helper functions and classes are found in sc_rna_variants.utils
    # TODO explain more the 'epilog
    """
    parser = argparse.ArgumentParser(
        formatter_class=ArgparserFormater,
        description="""This script aggregates and learns statistics on the output files from 'scrnavariants.py'.""",
        epilog='''Outputs aggregated tsv, figures and statistics.'''
    )

    # positional arguments
    parser.add_argument('input_dir', type=assert_is_directory,
                        help='folder with raw_stats.tsv and raw_umutated_stats.tsv files from scrnavariants.py')
    parser.add_argument('output_dir', help='folder for step outputs', type=make_output_dir)

    # optional arguments
    parser.add_argument('--min_cb_per_pos', default=5, type=int,
                        help='position with less cell barcodes will be filtered')
    parser.add_argument('--min_mutation_umis', default=10, type=int,
                        help='position with less mutated UMIs will be filtered')
    parser.add_argument('--min_total_umis', default=20, type=int,
                        help='position with less number of mutated + unmutated UMIs will be filtered')
    parser.add_argument('--min_mutation_rate', default=0.1, type=int,
                        help='position with less rate of mutation will be filtered')

    # Meta arguments
    parser.add_argument('--threads', type=int,
                        help='number of available threads', default=1)
    parser.add_argument('--log-file', default=None,
                        help='a log file for tracking the program\'s progress')
    parser.add_argument('--sname', type=str, help='sample name to add to outputs')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    startTime = datetime.now()
    args = parse_arguments()

    # initialize logger
    sc_rna_variants.config_logging(args.log_file)
    logger = logging.getLogger("agregation_per_position")
    logger.info('Aggregation per position started')
    logger.debug('Running with parameters:\n%s' % '\n'.join(
        ['%s: %s' % (key, value) for key, value in vars(args).items()]))

    # run statistics analysis
    sc_rna_variants.stats_generator.run(args)

    print(datetime.now() - startTime)
    logger.info('Step 4 finished')
    print("Step 4 finished")