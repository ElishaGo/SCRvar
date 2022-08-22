#!/usr/bin/python
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

import sc_rna_variants.utils
import sc_rna_variants.bio_functions
import sc_rna_variants.steps_runner



def run_step_5(output_dir, SCvar_aggregated_bed_file, min_cb_per_pos, min_mutation_umis, min_total_umis,
              min_mutation_rate, snp_db_path, editing_db_path, sname):
    # run step 5
    out_dir5 = os.path.join(output_dir, 'step5_filtering_and_DB_intersections_effects')
    os.makedirs(out_dir5, exist_ok=True)
    sc_rna_variants.steps_runner.run_step5(out_dir5, SCvar_aggregated_bed_file, min_cb_per_pos, min_mutation_umis, min_total_umis,
              min_mutation_rate, snp_db_path, editing_db_path, sname)


#########################################################################################################
def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(formatter_class=sc_rna_variants.utils.ArgparserFormater, description="", )

    # positional arguments
    parser.add_argument('output_dir', type=sc_rna_variants.utils.assert_is_directory, help='folder for outputs')
    parser.add_argument('SCvar_aggregated_bed_file', type=sc_rna_variants.utils.assert_is_file, help='last output of SCvar main program')

    # optional arguments
    parser.add_argument('--editing_db_path', help='path to known editing sites file')
    parser.add_argument('--snp_db_path', help='path to known SNP sites file')
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
                        default=os.path.join(sys.argv[1], "5.filtering_and_analysis.log"),
                        help='a log file for tracking the program\'s progress')
    parser.add_argument('--sname', type=str, help='sample name to add to outputs')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    startTime = datetime.now()
    args = parse_arguments()

    # initialize logger
    sc_rna_variants.config_logging(args.log_file)
    logger = logging.getLogger("positions_filtering_and_plots")
    logger.info('positions_filtering_and_plots started')
    logger.debug('Running with parameters:\n%s' % '\n'.join(
        ['%s: %s' % (key, value) for key, value in vars(args).items()]))

    # run filtering and create plots
    run_step_5(args.output_dir, args.SCvar_aggregated_bed_file, args.min_cb_per_pos, args.min_mutation_umis, args.min_total_umis,
              args.min_mutation_rate, args.snp_db_path, args.editing_db_path, args.sname)

    print(datetime.now() - startTime)
    logger.info('Step 5 finished')
