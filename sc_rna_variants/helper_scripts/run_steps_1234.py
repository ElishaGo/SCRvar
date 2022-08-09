"""
This script creates a text file in bsub format to run the pipeline on LSF system (like WEXAC),
and execute it as a bsub command
"""
import os
import sys
import argparse
import logging
import subprocess
from datetime import datetime

from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent.parent.absolute()) + os.path.sep)  # for development environments

import sc_rna_variants.utils
import sc_rna_variants.bio_functions
import sc_rna_variants.steps_runner
from sc_rna_variants import config_logging


def makedir(dir_path):
    os.makedirs(dir_path, exist_ok=True)
    return dir_path


def run_steps_1234(args):
    code_dir = str(Path(__file__).parent.parent.parent.absolute()) + os.path.sep
    logger.info("count number of reads per barcode\n")
    count_reads_exec = f"{code_dir}/sc_rna_variants/count_reads_per_barcode_in_bam.sh {args.input_bam} {args.sample_output_dir} {args.threads} {args.sname}"
    logger.info("running command: " + count_reads_exec)
    # os.system(count_reads_exec)

    # step1 - filter bam file
    step1_output_dir = os.path.join(args.output_folder, 'step1_filtered_bam_files')
    logger.info("# STEP 1 - filter bam file by filter list\n")
    os.makedirs(step1_output_dir, exist_ok=True)
    step1_exec = f"python {code_dir}/scripts/step1_filter_bam.py {args.input_bam} {step1_output_dir} --filtered-barcodes-list {args.filtered_barcodes_list} --threads {args.threads}"
    # logger.info("running command: " + step1_exec)
    # os.system(step1_exec)
    sc_rna_variants.steps_runner.run_step1(args.input_bam, args.filtered_barcodes_list,
              args.min_mapq, args.cigar_clipping_allowed,
              args.max_gene_length, args.max_no_basecall,
              args.tag_for_umi, args.tag_for_cell_barcode,
              step1_output_dir, args.threads)
    
    # get path to filtered bam file
    filtered_bam_path = str(
        os.path.join(args.sample_output_dir, step1_output_dir,
                     "1." + os.path.basename(args.input_bam).replace(".bam", '') + "_filtered.bam"))

    # step 2 - keep only reads from genes with htseq
    step2_output_dir = os.path.join(args.output_folder, 'step2_bam_gene_filter')
    editing_gtf_bam_intersect = None
    snp_gtf_bam_intersect = None
    if args.editing_db_path:
        editing_gtf_bam_intersect = os.path.join(step2_output_dir, f'2.editing.genecode.{args.sname}_intersect.bed')
    if args.snp_db_path:
        snp_gtf_bam_intersect = os.path.join(step2_output_dir, f'2.snp.genecode.{args.sname}_intersect.vcf')
    logger.info("STEP 2 - bam genes filter\n")
    step2_exec = f"sh {code_dir}/scripts/step2_bam_gene_filter.sh {filtered_bam_path} {step2_output_dir} {args.annotation_gtf} {args.editing_db_path} {args.snp_db_path} {editing_gtf_bam_intersect} {snp_gtf_bam_intersect} {args.sname} {args.threads}"
    os.makedirs(step2_output_dir, exist_ok=True)
    logger.info("running command: " + step2_exec)
    subprocess.run([step2_exec])

    # step 3 - create mismatch dictionary
    step3_output_dir = os.path.join(args.output_folder, 'step3_mismatch_dictionary')
    logger.info("STEP 3 - create mismatch dictionary\n")
    os.makedirs(step3_output_dir, exist_ok=True)
    step3_exec = f"python {code_dir}/scripts/step3_mismatch_dictionary.py {step2_output_dir}/2.{args.sname}.gene_filter.bam {args.genome_fasta} {step3_output_dir} --threads {args.threads}"
    logger.info("running command: " + step3_exec)
    # subprocess.run([step3_exec])
    sc_rna_variants.steps_runner.run_step3(args.input_bam, args.genome_fasta, args.tag_for_umi, args.tag_for_cell_barcode, step3_output_dir, args.threads)

    # step 4 - Aggregation per position + statistics
    step4_output_dir = os.path.join(args.output_folder, 'step4_aggregation_per_position_and_statistics')
    logger.info('STEP 4 - aggregation per position + statistics\n')
    os.makedirs(step4_output_dir, exist_ok=True)
    step4_exec = f"python {code_dir}/scripts/step4_aggregation_per_position.py {step3_output_dir} {step4_output_dir} {editing_gtf_bam_intersect} {snp_gtf_bam_intersect} {args.annotation_gtf} --sname {args.sname}"
    logger.info("running command: " + step4_exec)
    # subprocess.run([step4_exec])
    sc_rna_variants.steps_runner.run_step4(args.input_dir, step4_output_dir, args.annotation_gtf, args.snp_db_path, args.editing_db_path)

#########################################################################################################
def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(description="""A script to set parameter to a bsub file and send to bsub""", )

    # positional arguments
    parser.add_argument('sample_output_dir', type=sc_rna_variants.utils.assert_is_directory, help='directory to programs outputs')
    parser.add_argument('input_bam', type=sc_rna_variants.utils.assert_is_file, help='input bam file')
    parser.add_argument('genome_fasta', type=sc_rna_variants.utils.assert_is_file, help='genome reference')
    parser.add_argument('annotation_gtf', type=sc_rna_variants.utils.assert_is_file, help='genecode annotation file (gtf format)')

    # optional arguments
    parser.add_argument('--filtered-barcodes-list',
                        type=sc_rna_variants.utils.filtered_barcodes_processing,
                        # returns a set with the barcodes names
                        help='''Text/tsv file with a list of cell barcodes as first column. Counts only these cells. Please note GEM-well numbers are ignored''')
    parser.add_argument('--editing_db_path', help='Editing sites data base in bed format'
                        # , '/home/labs/bioservices/shared/rarevar/data/DataBases/editing/0.editing_A_I.genecode_intersect.bed'
                        )

    parser.add_argument('--snp_db_path', help='Known SNP sites data base in vcf format'
                        # , '/home/labs/bioservices/shared/rarevar/data/DataBases/snp_vcf/0.snp.gencode_intersect.vcf'
                        )

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
    parser.add_argument('--log-file',
                        default=os.path.join(sys.argv[1], "steps_1234_main.log"),
                        help='a log file for tracking the program\'s progress')

    parser.add_argument('--threads', type=int, default=1, help='Number of cores to use')

    parser.add_argument('--sname', default='', help='sample name to add to outputs')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    startTime = datetime.now()
    os.path.dirname(os.path.realpath(__file__))
    args = parse_arguments()

    # initialize logger
    config_logging(args.log_file)
    logger = logging.getLogger("positions_filtering_and_plots")
    logger.info('positions_filtering_and_plots started')
    logger.debug('Running with parameters:\n%s' % '\n'.join(
        ['%s: %s' % (key, value) for key, value in vars(args).items() if key != 'filtered_barcodes_list'])
                 )

    run_steps_1234(args)

    print(datetime.now() - startTime)
    print("finished")
