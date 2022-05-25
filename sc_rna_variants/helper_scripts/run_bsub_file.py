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
from sc_rna_variants.general_utils import  save_how_to

def write_bsub_execution_parameters(f, args):
    """write LSF system parameters"""
    f.write("#BSUB -q {}\n".format(args.q))
    f.write("#BSUB -J {}\n".format(args.J))
    f.write("#BSUB -oo serverlog_files/{}_%J.out\n".format(args.J))
    f.write("#BSUB -eo serverlog_files/{}_%J.err\n".format(args.J))
    f.write('#BSUB -R "rusage[mem={}] span[hosts={}]"\n'.format(args.rusage, args.hosts))
    f.write("#BSUB -n {}\n".format(args.n))
    f.write("#BSUB -cwd {}\n".format(args.sample_output_dir))


def write_params(f, args):
    """write the parameters used in the job"""
    f.write("\n### Parameters used in this job: \n%s" % '\n'.join(
        ['# %s: %s' % (key, value) for key, value in vars(args).items() if
         key not in ['q', 'J', 'rusage', 'hosts', 'n']]))


def write_meta_commands(f, args):
    write_bsub_execution_parameters(f, args)

    write_params(f, args)

    f.write("\n\n\n### Run the commands in serial:\n")

    # load environment and modules
    f.write("# load environment and modules\n")
    f.write(". /home/labs/bioservices/services/miniconda2/etc/profile.d/conda.sh;conda activate rarevar;module load "
            "SAMtools;module load bamtools;module load bedtools\n\n")


def write_pipelines_scripts_execution_commands(f, args):
    # TODO: ask what step to put this
    f.write("# count number of reads per barcode\n")
    f.write(f"{os.getcwd()}/sc_rna_variants/count_reads_per_barcode_in_bam.sh {args.bam_file} {args.sample_output_dir} {args.sname}")

    # step1 - filter bam file
    step1_output_dir = 'step1_filtered_bam_files'
    f.write("# STEP 1 - filter bam file by filter list\n")
    f.write(f"mkdir {step1_output_dir}\n")
    f.write(
        f"python {os.getcwd()}/scripts/step1_filter_bam.py {args.bam_file} {step1_output_dir} --filtered-barcodes-list {args.filter_list_bam} --threads {args.n}\n\n")

    # get path to filtered bam file
    filtered_bam_path = str(
        os.path.join(args.sample_output_dir, step1_output_dir, "1_" + os.path.basename(args.bam_file) + "_filtered.bam"))

    # step 2 - keep only reads from genes with htseq
    step2_output_dir = 'step2_bam_gene_filter'
    editing_gtf_intersect = '/home/labs/bioservices/shared/rarevar/data/DataBases/REDIportal/0_editing_A_I.genecode_intersect.bed'
    snp_gtf_intersect = '/home/labs/bioservices/shared/rarevar/data/DataBases/snp_vcf/0_snp_A.gencode_intersect.vcf'
    f.write("# STEP 2 - bam genes filter\n")
    f.write(f"mkdir {step2_output_dir}\n")
    f.write(
        f"sh {os.getcwd()}/scripts/step2_bam_gene_filter.sh {filtered_bam_path} {step2_output_dir} {args.annotation_gtf} {editing_gtf_intersect} {snp_gtf_intersect} {args.sname} {args.n}\n\n")

    # step 3 - create mismatch dictionary
    step3_output_dir = 'step3_mismatch_dictionary'
    f.write("# STEP 3 - create mismatch dictionary\n")
    f.write(f"mkdir {step3_output_dir}\n")
    f.write(
        f"python {os.getcwd()}/scripts/step3_mismatch_dictionary.py {step2_output_dir}/2_{args.sname}.gene_filter.bam {args.genome_ref} {step3_output_dir} --threads {args.n}\n\n")

    # step 4 - Aggregation per position + statistics
    step4_output_dir = 'step4_aggregation_per_position_and_statistics'
    f.write('# STEP 4 - aggregation per position + statistics\n')
    f.write(f"mkdir {step4_output_dir}\n")
    f.write(
        f"python {os.getcwd()}/scripts/step4_aggregation_per_position.py {step3_output_dir} {step4_output_dir} {args.snp_vcf} {args.editing_DB} --sname {args.sname} --threads {args.n}\n\n")

    # step 5 - filtering positions and SNP/editing DB intersections
    step5_output_dir = 'step5_filtering_positions_and_SNP_editing_DB_intersections'
    f.write('# STEP 5 - filtering positions and SNP/editing DB intersections\n')
    f.write(f"mkdir {step5_output_dir}\n")
    f.write(
        f'python {os.getcwd()}/scripts/step5_filtering_positions_and_snp_editing_DB_intersections.py {step4_output_dir} {step5_output_dir} {args.snp_vcf} {args.editing_DB} --sname {args.sname}\n\n')

    # step 6 - gene level analysis
    step6_output_dir = 'step6_gene_level'
    f.write('# STEP 6 - gene level\n')
    f.write(f"mkdir {step6_output_dir}\n")
    f.write(
        f"python {os.getcwd()}/scripts/step6_gene_level.py {step4_output_dir} {step6_output_dir} {os.path.join(step3_output_dir, '3_mismatch_dictionary.bed')} {os.path.join(args.sample_output_dir, f'raw_bam_reads_per_barcode_count.{args.sname}.csv')} {args.annotation_gtf} --barcode_clusters {args.barcode_clusters} --sname {args.sname}\n\n")


def create_job_file(args):
    # open folder for log files
    os.system("mkdir {}/serverlog_files".format(args.sample_output_dir))

    # write command and pipeline to bsub file
    f = open(os.path.join(args.sample_output_dir, f"bsub_file_SCrarevar_pipline_{args.sname}.txt"), "w")
    write_meta_commands(f, args)
    write_pipelines_scripts_execution_commands(f, args)
    f.close()


#########################################################################################################
def parse_arguments(arguments=None):
    """argument parsing wrapper function"""
    parser = argparse.ArgumentParser(description="""A script to set parameter to a bsub file and send to bsub""", )

    # positional arguments
    parser.add_argument('sample_output_dir', type=assert_is_directory, help='the directory to run file and get outputs')

    parser.add_argument('bam_file', help='Input bam file')

    # optional arguments
    parser.add_argument('--q', default='bio', help='''queue to run on WEXAC''')

    parser.add_argument('--J', default='jobfile_scrarv', help='job name on WEXAC')

    parser.add_argument('--rusage', type=int, default=1024, help='Memory to use in bytes')

    parser.add_argument('--hosts', type=int, default=1, help='Number of hosts to use')

    parser.add_argument('--n', type=int, default=1, help='Number of cores to use')

    parser.add_argument('--filter_list_bam', help='List of cell barcodes to use in format as in  the bam file')

    parser.add_argument('--annotation_gtf',
                        default="/home/labs/bioservices/shared/rarevar/data/DataBases/genecode_gtf/gencode.v37.annotation.gtf",
                        help='gtf annotation file to find gene sites')

    parser.add_argument('--genome_ref',
                        default="/home/labs/bioservices/services/expression_references/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
                        help='genome reference')

    parser.add_argument('--editing_DB',
                        default="/home/labs/bioservices/shared/rarevar/data/DataBases/REDIportal/0.editing_A_I.bed",
                        help='Editing sites data base in bed format')

    parser.add_argument('--snp_vcf',
                        default="/home/labs/bioservices/shared/rarevar/data/DataBases/snp_vcf/snp_chr_sorted.vcf",
                        help='Known SNP sites data base in vcf format')

    parser.add_argument('--barcode_clusters', help='table with barcodes and associated clusters analysed by Seurat')

    parser.add_argument('--sname', default='', help='File name. Will be used in inputs and outputs')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    args = parse_arguments()
    create_job_file(args)
    save_how_to(out_dir=args.sample_output_dir, sub_cmd_str='_run_bsub')
    os.system(f"bsub < {args.sample_output_dir}/bsub_file_SCrarevar_pipline_{args.sname}.txt")
