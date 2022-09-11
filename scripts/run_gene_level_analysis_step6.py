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


def run_step6(SCRvar_aggregated_bed_file, output_dir, read_per_barcode_raw_bam, min_cb_per_pos, min_mutation_umis,
              min_total_umis, min_mutation_rate, reditools_data, annotation_gtf, mismatch_dict_bed, barcode_clusters,
              atacseq_gcoverage_min, atacseq_gfrequency_min, sname):

    out_dir6 = os.path.join(output_dir, 'step6_gene_level')
    os.makedirs(out_dir6, exist_ok=True)
    sc_rna_variants.steps_runner.run_step6(SCRvar_aggregated_bed_file, out_dir6, read_per_barcode_raw_bam,
                                           min_cb_per_pos, min_mutation_umis, min_total_umis, min_mutation_rate,
                                           reditools_data, annotation_gtf, mismatch_dict_bed, barcode_clusters,
                                           atacseq_gcoverage_min, atacseq_gfrequency_min, sname)


#########################################################################################################
def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument('output-dir', type=sc_rna_variants.utils.assert_is_directory, help='folder for outputs')
    parser.add_argument('scrvar-aggregated-bed', type=sc_rna_variants.utils.assert_is_file,
                        help='SCRvar step4 output')
    parser.add_argument('mismatch-dict-bed', type=sc_rna_variants.utils.assert_is_file, help='mismatch_dictionary from SCRvar step3')
    parser.add_argument('annotation-gtf', type=sc_rna_variants.utils.assert_is_file, help='path to gtf file')
    parser.add_argument('--read-per-barcode-raw-bam', type=sc_rna_variants.utils.assert_is_file, help='count of reads per CB in raw bam file')

    # optional arguments
    parser.add_argument('--barcodes-clusters', type=sc_rna_variants.utils.assert_is_file, help='barcodes and clusters analysed by Seurat')
    parser.add_argument('--min-cb-per-pos', default=5, type=int,
                        help='position with fewer cell barcodes will be filtered')
    parser.add_argument('--min-mutation-umis', default=10, type=int,
                        help='position with fewer mutated UMIs will be filtered')
    parser.add_argument('--min-total-umis', default=20, type=int,
                        help='position with fewer mutated + unmutated UMIs will be filtered')
    parser.add_argument('--min-mutation-rate', default=0.1, type=int,
                        help='position with lower mutation rate will be filtered')
    parser.add_argument('--reditools-data', type=str, help='path to atacseq file')
    parser.add_argument('--atacseq-gcoverage-min', type=int, default=5)
    parser.add_argument('--atacseq-gfrequency-min', type=float, default=0.2)

    # Meta arguments
    parser.add_argument('--log-file',
                        default=os.path.join(sys.argv[1], "6.gene_level_analysis.log"),
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

    print(['%s: %s' % (key, value) for key, value in vars(args).items() if key != 'barcodes_clusters'])
    run_step6(args.scrvar_aggregated_bed, args.output_dir,
              args.read_per_barcode_raw_bam,
              args.min_cb_per_pos,
              args.min_mutation_umis, args.min_total_umis, args.min_mutation_rate, args.reditools_data,
              args.annotation_gtf,
              args.mismatch_dict_bed, args.barcodes_clusters, args.atacseq_gcoverage_min, args.atacseq_gfrequency_min,
              args.sname)

    print(datetime.now() - startTime)
    logger.info('Step 6 finished')
