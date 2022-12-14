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


def write_bsub_execution_parameters(f, args):
    """write LSF system parameters"""
    f.write("#BSUB -q {}\n".format(args.q))
    f.write("#BSUB -J {}\n".format(args.sname))
    f.write("#BSUB -oo serverlog_files/{}.out\n".format(args.J))
    f.write("#BSUB -eo serverlog_files/{}.err\n".format(args.J))
    f.write('#BSUB -R "rusage[mem={}] span[hosts={}]"\n'.format(args.rusage, args.hosts))
    f.write("#BSUB -n {}\n".format(args.threads))
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
    f.write(
        ". /home/labs/bioservices/services/miniconda2/etc/profile.d/conda.sh;conda activate scvar;module load SAMtools;module load bedtools;module load bamtools\n\n")


def write_pipelines_scripts_execution_commands(f, args):
    code_dir = "/home/labs/bioservices/shared/rarevar/scrarevar/"
    f.write("# count number of reads per barcode\n")
    f.write(
        f"#{code_dir}/sc_rna_variants/count_reads_per_barcode_in_bam.sh {args.input_bam} {args.sample_output_dir} {args.threads} {args.sname}\n\n")

    # step1 - filter bam file
    step1_output_dir = 'step1_filtered_bam_files'
    f.write("# STEP 1 - filter bam file by filter list\n")
    f.write(f"mkdir {step1_output_dir}\n")
    f.write(
        f"python {code_dir}/scripts/step1_filter_bam.py {args.input_bam} {step1_output_dir} --barcodes-cluster-file {args.barcodes_clusters} --threads {args.threads}\n\n")

    # get path to filtered bam file
    filtered_bam_path = str(
        os.path.join(args.sample_output_dir, step1_output_dir,
                     "1." + os.path.basename(args.input_bam).replace(".bam", '') + "_filtered.bam"))

    # step 2 - keep only reads from genes with htseq
    step2_output_dir = 'step2_bam_gene_filter'
    editing_gtf_bam_intersect = None
    snp_gtf_bam_intersect = None
    if args.editing_db_path:
        editing_gtf_bam_intersect = os.path.join(step2_output_dir, f'2.editing.genecode.{args.sname}_intersect.bed')
    if args.snp_db_path:
        snp_gtf_bam_intersect = os.path.join(step2_output_dir, f'2.snp.genecode.{args.sname}_intersect.vcf')
    f.write("# STEP 2 - bam genes filter\n")
    f.write(f"mkdir {step2_output_dir}\n")
    f.write(
        f"sh {code_dir}/scripts/step2_bam_gene_filter.sh {filtered_bam_path} {step2_output_dir} {args.annotation_gtf} {args.editing_db_path} {args.snp_db_path} {editing_gtf_bam_intersect} {snp_gtf_bam_intersect} {args.sname} {args.threads}\n\n")

    # step 3 - create mismatch dictionary
    step3_output_dir = 'step3_mismatch_dictionary'
    f.write("# STEP 3 - create mismatch dictionary\n")
    f.write(f"mkdir {step3_output_dir}\n")
    f.write(
        f"python {code_dir}/scripts/step3_mismatch_dictionary.py {step2_output_dir}/2.{args.sname}.gene_filter.bam {args.genome_fasta} {step3_output_dir} --threads {args.threads}\n\n")

    # step 4 - Aggregation per position + statistics
    step4_output_dir = 'step4_aggregation_per_position_and_statistics'
    f.write('# STEP 4 - aggregation per position + statistics\n')
    f.write(f"mkdir {step4_output_dir}\n")
    f.write(
        f"python {code_dir}/scripts/step4_aggregation_per_position.py {step3_output_dir} {step4_output_dir} {args.annotation_gtf} {editing_gtf_bam_intersect} {snp_gtf_bam_intersect} --sname {args.sname}\n\n")

    # step 5 - filtering positions and SNP/editing DB intersections
    step5_output_dir = 'step5_filtering_and_DB_intersections_effects'
    f.write('# STEP 5 - filtering positions and SNP/editing DB intersections\n')
    f.write(f"mkdir {step5_output_dir}\n")
    f.write(
        f'python {code_dir}/scripts/step5_filtering_and_DB_intersections_effects.py {step4_output_dir} {step5_output_dir} {editing_gtf_bam_intersect} {snp_gtf_bam_intersect} --sname {args.sname}\n\n')

    # step 6 - gene level analysis
    step6_output_dir = 'step6_gene_level'
    if args.barcodes_clusters:
        barcode_clusters = f"--barcode_clusters {args.barcode_clusters}"
    f.write('# STEP 6 - gene level\n')
    f.write(f"mkdir {step6_output_dir}\n")
    f.write(
        f"python {code_dir}/scripts/step6_gene_level.py {step4_output_dir} {step6_output_dir} {os.path.join(step3_output_dir, '3.mismatch_dictionary.bed')} {os.path.join(args.sample_output_dir, f'raw_bam_reads_per_barcode_count.{args.sname}.csv')} {args.annotation_gtf} {barcode_clusters} --sname {args.sname}\n\n")


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

    # TODO see how to make this explicit arguments but still mandatory
    # positional arguments
    parser.add_argument('sample_output_dir', type=assert_is_directory, help='the directory to run file and get outputs')

    parser.add_argument('input_bam', help='Input bam file')

    parser.add_argument('genome_fasta',
                        # default="/home/labs/bioservices/services/expression_references/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
                        help='genome reference')

    parser.add_argument('annotation_gtf',
                        # default="/home/labs/bioservices/shared/rarevar/data/DataBases/genecode_gtf/0.gencode.v37.annotation.gtf",
                        help='gtf annotation file to find gene sites')

    # optional arguments
    parser.add_argument('--barcodes-clusters', help='List of cell barcodes to use in format as in  the bam file')

    parser.add_argument('--editing_db_path', help='Editing sites data base in bed format'
                        # , '/home/labs/bioservices/shared/rarevar/data/DataBases/editing/0.editing_A_I.genecode_intersect.bed'
                        )

    parser.add_argument('--snp_db_path', help='Known SNP sites data base in vcf format'
                        # , '/home/labs/bioservices/shared/rarevar/data/DataBases/snp_vcf/0.snp.gencode_intersect.vcf'
                        )

    parser.add_argument('--sname', default='', help='File name. Will be used in inputs and outputs')

    # LSF system arguments
    parser.add_argument('--q', default='bio', help='''queue to run on WEXAC''')

    parser.add_argument('--J', default='jobfile_scrarv', help='job name on WEXAC')

    parser.add_argument('--rusage', type=int, default=1024, help='Memory to use in bytes')

    parser.add_argument('--hosts', type=int, default=1, help='Number of hosts to use')

    parser.add_argument('--threads', type=int, default=1, help='Number of cores to use')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    args = parse_arguments()
    create_job_file(args)
    save_how_to(out_dir=args.sample_output_dir, sub_cmd_str='_run_bsub')
    os.system(f"bsub < {args.sample_output_dir}/bsub_file_SCrarevar_pipline_{args.sname}.txt")
