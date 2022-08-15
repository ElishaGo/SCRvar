"""
This script creates a text file in bsub format to run the pipeline on LSF system (like WEXAC),
and execute it as a bsub command
"""
import os
import argparse
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parent.parent.parent.absolute()) + os.path.sep)  # for development environments

from sc_rna_variants.utils import assert_is_directory
from sc_rna_variants.utils import save_how_to


def run_SCvar_pipeline(args):
    code_dir = str(Path(__file__).parent.parent.parent.absolute()) + os.path.sep
    print("count number of reads per barcode\n")
    os.system(
        f"{code_dir}/sc_rna_variants/count_reads_per_barcode_in_bam.sh {args.bam_file} {args.sample_output_dir} {args.sname} {args.threads}")

    # step1 - filter bam file
    step1_output_dir = 'step1_filtered_bam_files'
    print("# STEP 1 - filter bam file by filter list\n")
    os.system(f"mkdir {step1_output_dir}")
    os.system(
        f"python {code_dir}/scripts/step1_filter_bam.py {args.bam_file} {step1_output_dir} --filtered-barcodes-list {args.filter_list_bam} --threads {args.threads}")

    # get path to filtered bam file
    filtered_bam_path = str(
        os.path.join(args.sample_output_dir, step1_output_dir,
                     "1." + os.path.basename(args.bam_file).replace(".bam", '') + "_filtered.bam"))

    # step 2 - keep only reads from genes with htseq
    step2_output_dir = 'step2_bam_gene_filter'
    editing_gtf_intersect = os.path.join(args.editing_db_path, '0.editing_A_I.genecode_intersect.bed')
    snp_gtf_intersect = os.path.join(args.snp_db_path, '0.snp.gencode_intersect.vcf')
    editing_gtf_bam_intersect = os.path.join(step2_output_dir, f'2.editing.genecode.{args.sname}_intersect.bed')
    snp_gtf_bam_intersect = os.path.join(step2_output_dir, f'2.snp.genecode.{args.sname}_intersect.vcf')
    print("# STEP 2 - bam genes filter\n")
    os.system(f"mkdir {step2_output_dir}")
    os.system(
        f"sh {code_dir}/scripts/step2_bam_gene_filter.sh {filtered_bam_path} {step2_output_dir} {args.annotation_gtf} {editing_gtf_intersect} {snp_gtf_intersect} {editing_gtf_bam_intersect} {snp_gtf_bam_intersect} {args.sname} {args.threads}")

    # step 3 - create mismatch dictionary
    step3_output_dir = 'step3_mismatch_dictionary'
    print("# STEP 3 - create mismatch dictionary\n")
    os.system(f"mkdir {step3_output_dir}")
    os.system(
        f"python {code_dir}/scripts/step3_mismatch_dictionary.py {step2_output_dir}/2.{args.sname}.gene_filter.bam {args.genome_fasta} {step3_output_dir} --threads {args.threads}")

    # step 4 - Aggregation per position + statistics
    step4_output_dir = 'step4_aggregation_per_position_and_statistics'
    print('# STEP 4 - aggregation per position + statistics\n')
    os.system(f"mkdir {step4_output_dir}")
    os.system(
        f"python {code_dir}/scripts/step4_aggregation_per_position.py {step3_output_dir} {step4_output_dir} {editing_gtf_bam_intersect} {snp_gtf_bam_intersect} --sname {args.sname} --threads {args.threads}")

    # step 5 - filtering positions and SNP/editing DB intersections
    step5_output_dir = 'step5_filtering_positions_and_SNP_editing_DB_intersections'
    print('# STEP 5 - filtering positions and SNP/editing DB intersections\n')
    os.system(f"mkdir {step5_output_dir}")
    os.system(
        f'python {code_dir}/scripts/step5_filtering_positions_and_snp_editing_DB_intersections.py {step4_output_dir} {step5_output_dir} {snp_gtf_bam_intersect} {editing_gtf_bam_intersect} --sname {args.sname}')

    # step 6 - gene level analysis
    step6_output_dir = 'step6_gene_level'
    if args.barcode_clusters:
        barcode_clusters = f"--barcode_clusters {args.barcode_clusters}"
    print('# STEP 6 - gene level\n')
    os.system(f"mkdir {step6_output_dir}")
    os.system(
        f"python {code_dir}/scripts/step6_gene_level.py {step4_output_dir} {step6_output_dir} {os.path.join(step3_output_dir, '3.mismatch_dictionary.bed')} {os.path.join(args.sample_output_dir, f'raw_bam_reads_per_barcode_count.{args.sname}.csv')} {args.annotation_gtf} {barcode_clusters} --sname {args.sname}")


#########################################################################################################
def parse_arguments(arguments=None):
    """argument parsing wrapper function"""
    parser = argparse.ArgumentParser(description="""A script to set parameter to a bsub file and send to bsub""", )

    # positional arguments
    parser.add_argument('sample_output_dir', type=assert_is_directory, help='the directory to run file and get outputs')

    parser.add_argument('bam_file', help='Input bam file')

    parser.add_argument('genome_fasta',
                        default="/home/labs/bioservices/services/expression_references/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
                        help='genome reference')

    # optional arguments
    parser.add_argument('--filter_list_bam', help='List of cell barcodes to use in format as in  the bam file')

    parser.add_argument('--annotation_gtf',
                        default="/home/labs/bioservices/shared/rarevar/data/DataBases/genecode_gtf/0.gencode.v37.annotation.gtf",
                        help='gtf annotation file to find gene sites')

    parser.add_argument('--editing_db_path',
                        default="/home/labs/bioservices/shared/rarevar/data/DataBases/editing",
                        help='Editing sites data base in bed format')

    parser.add_argument('--snp_db_path',
                        default="/home/labs/bioservices/shared/rarevar/data/DataBases/snp_vcf",
                        help='Known SNP sites data base in vcf format')

    parser.add_argument('--barcode_clusters', help='table with barcodes and associated clusters analysed by Seurat')

    parser.add_argument('--sname', default='', help='File name. Will be used in inputs and outputs')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    args = parse_arguments()
    run_SCvar_pipeline(args)
    save_how_to(out_dir=args.sample_output_dir, sub_cmd_str='_run_SCvar_pipeline')
