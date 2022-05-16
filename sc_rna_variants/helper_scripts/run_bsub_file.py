"""
This script creates a text file in bsub format to run the pipeline on LSF system (like WEXAC),
and execute it as a bsub command
"""
import os
import argparse
from sc_rna_variants.utils import assert_is_directory


def create_job_file(args):
    f = open(os.path.join(args.working_dir, "bsub_file_SCrarevar_pipline.txt"), "w")

    # open folder for log files
    os.system("mkdir {}/serverlog_files/".format(args.working_dir))

    # write the LSF system parameters
    f.write("#BSUB -q {}\n".format(args.q))

    f.write("#BSUB -J {}\n".format(args.J))

    f.write("#BSUB -oo serverlog_files/{}_%J.out\n".format(args.J))

    f.write("#BSUB -eo serverlog_files/{}_%J.err\n".format(args.J))

    f.write('#BSUB -R "rusage[mem={}] span[hosts={}]"\n'.format(args.rusage, args.hosts))

    f.write("#BSUB -n {}\n".format(args.n))

    f.write("#BSUB -cwd {}\n".format(args.working_dir))

    # write the parameters used in the job
    f.write("\n### Parameters used in this job: \n%s" % '\n'.join(
        ['# %s: %s' % (key, value) for key, value in vars(args).items() if
         key not in ['q', 'J', 'rusage', 'hosts', 'n']]))

    f.write("\n\n\n### Run the commands in serial:\n")
    
    # load environment and modules
    f.write("# load environment and modules\n")
    f.write(". /home/labs/bioservices/services/miniconda2/etc/profile.d/conda.sh;conda activate rarevar;module load "
            "SAMtools;module load bamtools;module load bedtools\n\n")

    # TODO: ask what step to put this
    f.write("# count number of read per barcode\n")
    f.write(
        f"samtools view {args.bam_file} | grep NH:i:1 | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_barcode\n\n")

    # step1 - filter bam file
    step1_output_dir = 'step1_filtered_bam_files/'
    f.write("# STEP 1 - filter bam file by filter list\n")
    f.write("mkdir step1_filtered_bam_files\n")
    f.write(f"python {os.getcwd}/scripts/step1_filter_bam.py {args.bam_file} {step1_output_dir} --filtered-barcodes-list {args.filter_list_bam} --threads {args.n}\n\n")

    # get path to filtered bam file
    filtered_bam_path = str(
        os.path.join(args.working_dir, step1_output_dir, "1_" + os.path.basename(args.bam_file) + "_filtered.bam"))

    # step 2 - keep only reads from genes with htseq
    step2_output_dir = 'step2_bam_gene_filter/'
    f.write("# STEP 2 - bam genes filter\n")
    f.write(f"mkdir {step2_output_dir}\n")
    f.write(f"sh {os.getcwd}/scripts/step2_bam_gene_filter.sh {filtered_bam_path} {step2_output_dir} {args.annotation_gtf} {args.fname} {args.n}\n\n")
    
    # step 3 - create mismatch dictionary
    step3_output_dir = 'step3_mismatch_dictionary/'
    f.write("# STEP 3 - create mismatch dictionary\n")
    f.write(f"mkdir {step3_output_dir}\n")
    f.write(f"python {os.getcwd}/scripts/scrnavariants.py {step2_output_dir}/2_{args.fname}.gene_filter_header.bam {args.genome_ref} {step3_output_dir} --threads {args.n}\n\n")

    # step 4 - Aggregation per position + statistics
    step4_output_dir = 'step4_aggregation_per_position_and_statistics/'
    f.write('# STEP 4 - aggregation per position + statistics\n')
    f.write(f"mkdir {step4_output_dir}\n")
    f.write(f"python {os.getcwd}/scripts/step4_aggregation_per_position.py {step3_output_dir} {step4_output_dir} --sname {args.fname} --threads {args.n}\n\n")

    # step 5 - filtering positions and SNP/editing DB intersections
    step5_output_dir = 'step5_filtering_positions_and_SNP_editing_DB_intersections/'
    f.write('# STEP 5 - filtering positions and SNP/editing DB intersections\n')
    f.write(f"mkdir {step5_output_dir}\n")
    f.write(f'python {os.getcwd}/scripts/step5_filtering_positions_and_snp_editing_DB_intersections.py {step4_output_dir} {step5_output_dir} {args.editing_DB} {args.snp_vcf} --sname {args.fname}\n\n')

    f.write(f'python {os.getcwd}/scripts/filter_snp.py step5_filtering_positions_and_SNP_editing_DB_intersections/5_aggregated_intersect.tsv {args.snp_clf_weights}')

    f.close()


#########################################################################################################
def parse_arguments(arguments=None):
    """argument parsing wrapper function"""
    parser = argparse.ArgumentParser(description="""A script to set parameter to a bsub file and send to bsub""",)

    # positional arguments
    parser.add_argument('working_dir', type=assert_is_directory, help='the directory to run file and get outputs')

    parser.add_argument('bam_file', type=str, help='Input bam file')

    parser.add_argument('fname', type=str, help='File name. Will be used in inputs and outputs')

    # optional arguments
    parser.add_argument('--q', type=str, default='bio', help='''queue to run on WEXAC''')

    parser.add_argument('--J', default='jobfile_scrarv', help='job name on WEXAC', type=str)

    parser.add_argument('--rusage', type=int, default=1024, help='Memory to use in bytes')

    parser.add_argument('--hosts', type=int, default=1, help='Number of hosts to use')

    parser.add_argument('--n', type=int, default=1, help='Number of cores to use')

    parser.add_argument('--filter_list_bam', type=str,
                        help='List of cell barcodes to use in format as in  the bam file')

    parser.add_argument('--annotation_gtf', type=str,
                        default="/home/labs/bioservices/shared/rarevar/data/DataBases/gencode.v37.annotation.gtf",
                        help='gtf annotation file to find gene cites')

    parser.add_argument('--genome_ref', type=str,
                        default="/home/labs/bioservices/services/expression_references/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
                        help='genome reference')

    parser.add_argument('--editing_DB', type=str,
                        default="/home/labs/bioservices/shared/rarevar/data/DataBases/edit_snp_DB/human/edit_TABLE1_hg38.txt",
                        help='Editing repetitive sites data base in bed format')

    # parser.add_argument('--edit_rep_bed', type=str,
    #                     default="/home/labs/bioservices/shared/rarevar/data/DataBases/edit_snp_DB/human/edit_rep.bed",
    #                     help='Editing repetitive sites data base in bed format')
    #
    # parser.add_argument('--edit_nonrep_bed', type=str,
    #                     default="/home/labs/bioservices/shared/rarevar/data/DataBases/edit_snp_DB/human/edit_nonrep.bed",
    #                     help='Editing non repetitive sites data base in bed format')

    parser.add_argument('--snp_vcf', type=str,
                        default="/home/labs/bioservices/shared/rarevar/data/DataBases/edit_snp_DB/human/snp_chr_sorted.vcf",
                        help='Known SNP sites data base in vcf format')

    parser.add_argument('--snp_clf_weights', type=str,
                        default="/home/labs/bioservices/shared/rarevar/data/trained_models/snp_classifier.joblib",
                        help='SNP classifier pretrained weights')

    return parser.parse_args(arguments)


if __name__ == '__main__':
    args = parse_arguments()
    create_job_file(args)
    os.system(f"bsub < {args.working_dir}/bsub_file_SCrarevar_pipline_{args.fname}.txt")