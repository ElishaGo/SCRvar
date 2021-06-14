import os
import argparse


def assert_is_directory(string):
    """Make sure the given string indicates an OS directory"""
    if not os.path.isdir(string):
        msg="\n%s is not a valid directory" % string
        raise argparse.ArgumentTypeError(msg)
    return string


def parse_arguments(arguments=None):
    """argument parsing wrapper function"""
    parser = argparse.ArgumentParser(
        # formatter_class=sc_rna_variants.utils.ArgparserFormater,
        description="""A script to set parameter to a bsub file and send to bsub""",
    )

    # positional arguments
    parser.add_argument('working_dir', type=assert_is_directory, help='the directory to run file and get outputs')

    parser.add_argument('bam_file', type=str, help='Input bam file')

    parser.add_argument('fname', type=str, help='File name. Will be used in inputs and outputs')

    # optional arguments
    parser.add_argument('--q', type=str,default='bio', help='''queue to run on WEXAC''')

    parser.add_argument('--J', default='jobfile', help='job name on WEXAC', type=str)

    # parser.add_argument('--oo', default='%J.out',  type=str, help='path to output log file')
    #
    # parser.add_argument('--eo', default='%J.err', type=str, help='path to error log file')

    parser.add_argument('--rusage', type=int, default=1024, help='Memory to use in bytes')

    parser.add_argument('--hosts', type=int, default=1, help='Number of hosts to use')

    parser.add_argument('--n', type=int, default=1, help='Number of cores to use')

    parser.add_argument('--cwd', type=str, help='Set working directory')

    parser.add_argument('--filter_list_bam', type=str,
                        help='List of cell barcodes to use in format as in  the bam file')

    parser.add_argument('--filter_list_scrarvar', type=str, help='List of cell barcodes to use in scrarvar program')

    parser.add_argument('--annotation_gtf', type=str,
                        default="/home/labs/bioservices/shared/rarevar/data/gencode.v37.annotation.gtf",
                        help='gtf annotation file to find gene cites')

    parser.add_argument('--genome_ref', type=str,
                        default="/shareDB/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa",
                        help='genome reference')

    parser.add_argument('--edit_rep_bed', type=str,
                        default="/home/labs/bioservices/shared/rarevar/data/edit_snp_DB/human/edit_rep.bed",
                        help='Editing repetitive sites data base in bed format')

    parser.add_argument('--edit_nonrep_bed', type=str,
                        default="/home/labs/bioservices/shared/rarevar/data/edit_snp_DB/human/edit_nonrep.bed",
                        help='Editing non repetitive sites data base in bed format')

    parser.add_argument('--snp_vcf', type=str,
                        default="/home/labs/bioservices/shared/rarevar/data/edit_snp_DB/human/snp_chr_sorted.vcf",
                        help='Known SNP sites data base in vcf format')

    parser.add_argument('--snp_clf_weights', type=str,
                        default="/home/labs/bioservices/shared/rarevar/data/trained_models/snp_classifier.joblib",
                        help='SNP classifier pretrained weights')

    return parser.parse_args(arguments)


def create_job_file(args):
    """function to create file with bsub foramt to run on LSF systm"""
    # open file
    f = open(os.path.join(args.working_dir, "bjobs_file.txt"), "w")

    # write the LSF system parameters
    f.write("#BSUB -q {}\n".format(args.q))

    f.write("#BSUB -J {}\n".format(args.J))

    f.write("#BSUB -oo {}_%J.out\n".format(args.J))

    f.write("#BSUB -eo {}_%J.err\n".format(args.J))

    f.write('#BSUB -R "rusage[mem={}] span[hosts={}]"\n'.format(args.rusage, args.hosts))

    f.write("#BSUB -n {}\n".format(args.n))

    f.write("#BSUB -cwd {}\n".format(args.args.working_dir))


    # write the commands to execute
    f.write("\n")
    f.write("# load conda environment and modules\n")
    f.write(
        ". /home/labs/bioservices/services/miniconda2/etc/profile.d/conda.sh;conda activate rarevar;module load "
        "samtools;module load bamtools;module load bedtools\n\n"
    )

    f.write("# make dirs\n")
    f.write("mkdir filtered_bam_files/ bam_statistics/ /scRarevar_output/ log_files/ statistics_ouputs\n\n")

    f.write("# filter bam file by filter list\n")
    f.write("python /home/labs/bioservices/shared/rarevar/code/scripts/filter_bam.py {bam} {filter_list}"
              " --output_folder ./filtered_bam_files --name_suffix {fname} --threads {n}\n\n".format(bam=args.bam_file,
                                                                                                 filter_list=args.filter_list_bam,
                                                                                                 fname=args.fname,
                                                                                                 n=args.n))
    f.write("# add chr to chromosome names in bam files\n")
    f.write(
        "samtools view -H {fname}_CBfiltered.bam | sed  -e '/SN:chr/!s/SN:\([0-9XY]*\)/SN:chr&/' -e "
        "'/SN:chrM/!s/SN:MT/SN:chrM&/' | samtools reheader - {fname}_CBfiltered.bam > {fname}_CBfiltered_chr.bam;"
        "samtools index {fname}_CBfiltered_chr.bam\n\n".format(
            fname=args.fname))

    f.write("# remove old bam files\n")
    f.write("rm {fname}_CBfiltered.bam {fname}_CBfiltered.bam.bai\n\n".format(fname=args.fname))

    f.write("# run ht-seq to filter non gene cites from bam\n")
    f.write("htseq-count -f bam -i gene_name -t gene -m union -s yes -o {fname}_htseq_gene.sam "
              "{fname}_CBfiltered_chr.bam {gtf} 3>&1 > {fname}_stranded_counts.txt\n\n".format(fname=args.fname,
                                                                                           gtf=args.annotation_gtf))
    f.write("# add header to the bam file\n")
    f.write("samtools view -H {fname}_CBfiltered_chr.bam | "
              "cat - {fname}_htseq_gene.sam > {fname}_htseq_gene_header.sam\n\n".format(fname=args.fname))

    f.write("# remove sam file\n")
    f.write("rm {fname}_htseq_gene.sam\n\n".format(fname=args.fname))

    f.write("Run scrarevar program\n")
    f.write("python /home/labs/bioservices/shared/rarevar/code/scRNAvariants/scripts/scrnavariants.py"
              " {fname}_htseq_gene_header.bam {genome_ref} /scRarevar_output/ --log-file /log_files/log_{fname}.txt "
              "--threads {n}\n\n".format(fname=args.fname, genome_ref=args.genome_ref, n=args.n))

    f.write('Run statistics analysis program\n')
    f.write("python /home/labs/bioservices/shared/rarevar/code/scRNAvariants/scripts/make_statistics.py "
            "scRarevar_output/raw_stats.tsv "
            "scRarevar_output/raw_unmutated_stats.tsv "
            "--output_folder statistics_ouput/ "
            "--log-file scRarevar_log_file/log_statistics_{fname}.txt\n\n".format(fname=args.fname))

    f.write('Find intersections with SNP and edit databases\n')
    f.write('python /home/labs/bioservices/shared/rarevar/code/scRNAvariants/scripts/visul/make_venn.py '
            'statistics_ouput/ {edit_rep_bed} {edit_nonrep_bed} {snp_vcf}\n\n'.format(edit_rep_bed=args.edit_rep_bed,
                                                                                   edit_nonrep_bed=args.edit_nonrep_bed,
                                                                                   snp_vcf=args.snp_vcf))

    f.write('Find intersections with SNP and edit databases\n')
    f.write('python /home/labs/bioservices/shared/rarevar/code/scripts/filter_snp.py '
            'statistics_ouput/aggregated_intersect.tsv '
            '{snp_clf}'.format(snp_clf=args.snp_clf_weights))
    
    f.close()


if __name__ == '__main__':
    args = parse_arguments()

    create_job_file(args)

    # os.system("bsub < {}bjobs_file.txt".format(args.working_dir))