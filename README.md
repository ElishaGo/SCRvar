# Single Cell RNA variants 

A script to locate cases of RNA modifications in single cell RNAseq data.

Given a BAM file with cellRanger format, returns a BED formated file where each
row represents the modifications a certain cell has in a certain position,
categorised by mutation type and UMI counts.

requirements: pysam package

conda environment setup example:
$ conda create -n pysam37
$ conda activate pysam37
$ conda config --env --add channels conda-forge
$ conda config --env --add channels bioconda
$ conda install python=3.7 pysam samtools spyder


~~
# Requirements
- UNIX OS
- python3.7 or higher
- samtools, bedtools, bamtools (For conda users, these can be installed with bioconda channel)

## Setup
It is recommended to first open a virtual environment (e.g python venv or conda).

After activating the environment the scvar package can be installed via PyPi:
~~~
pip install pandas matplotlib-venn (if using test.pypi)
pip install --index-url https://test.pypi.org/simple/ scRNAvariants==2.0.5
~~~

# Mandatory files
SCvar needs external files to run:
- Bam file to analyse
- Genome reference fasta file
- Transcriptome Genecoed file

It is important that the bam file and the other files are compatible.

## Optional files 
SCvar pipeline has some outputs properties available when running with additional databases files.

The file are: 
- Known editing sites (in bed format)
- Known SNP (in vcf format)
- Cell Barcodes and clusters file - a tsv/text file where the first column is cell barcodes of the sample (with no prefixex/suffixes), and the second column is clusters. The file must be without header line.


# Run SCvar pipeline
### 1. Run step0 script:
step 0 script will do some basic preproccessing to the external files (Trancsctiptome gtf, Editing bed file and SNP vcf file if provided).
This script needs to be run *only once* per external file. So each time you use a different file you'll need to run script 0:
~~~
step0_process_editing_and_snp_DB.py OUTPUT_DIR ANNOTATION_GTF --editing_DB_path EDITING_BD_PATH --fasta_path GENOME_FASTA --snp_DB_path SNP_DB_PATH
~~~
### 2. Run steps 1-4
steps 1-4 are the main body of the program. They will process the bam file and output a mismatch table, positions table and mid output files.

A single script execute steps 1-4, *and you need to run it once per bam file*.
~~~
run_steps_1234.py SAMPLE_OUTPUT_DIR INPUT_BAM GENOME_FASTA ANNOTATION_GTF --SNAME YOUR_SAMPLE_NAME OTHER_OPTIONAL_PARAMS
~~~

### 2. Run steps 5, 6 (To come) 
SCvar has it's own basic analysis outputs. 

SCvar can also run with some customized parameters to filter out some of the positions, thus ending up with different analysis outputs.
 ~~~
 step5_filtering_and_DB_intersections_effects.py INPUT_DIR OUTPUT_DIR EDITING_DB_PATH  SNP_DB_PATH --SNAME YOUR_SAMPLE_NAME OTHER_OPTIONAL_PARAMS 
 ~~~


# Reproduce our results
You can download the same files we used in our analysis and run them with the pipeline:
~~~
mkdir reproduce_results_data/

# download the bam file
wget -P reproduce_results_data ftp://dors2.weizmann.ac.il/merged_SRR452_SRR453.bam


# download the transcriptome gtf file
wget -P reproduce_results_data https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
gzip -d reproduce_results_data/gencode.v37.annotation.gtf.gz


# download the reference genome
wget -P reproduce_results_data ASK WHERE THE FILE IS DOWNLOADED FROM


# download the editing sites database from REDIportal
wget -P reproduce_results_data http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz
gzip -d reproduce_results_data/TABLE1_hg38.txt.gz


# download the SNP database
wget -P reproduce_results_data https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz
gzip -d reproduce_results_data/common_all_20180418.vcf.gz

# run the pipeline
step0_process_editing_and_snp_DB.py reproduce_results_data reproduce_results_data/gencode.v37.annotation.gtf --editing_DB_path reproduce_results_data/TABLE1_hg38.txt --fasta_path GENOME_FASTA --snp_DB_path reproduce_results_data/common_all_20180418.vcf
run_steps_1234.py reproduce_results_data/SCvar_output ../scrarevar/data/outputs/test/input_files/merged_chr_2_20_21_22.bam ../scrarevar/data/DB_edit_snp/human/genome/genome.fa reproduce_results_data/data_files_processed/genecode_gtf/0.gencode.v37.annotation.gtf --sname reproduce_scvar --editing_db_path reproduce_results_data/data_files_processed/editing/0.editing_A_I.genecode_intersect.bed --snp_db_path reproduce_results_data/data_files_processed/snp_vcf/0.snp.gencode_intersect.vcf

python sc_rna_variants/helper_scripts/run_bsub_file.py PATH_TO_OUTPUT reproduce_results_data/PATH_TO_BAM reproduce_results_data/GENOME_REF reproduce_results_data/ANNOTAION_GTF --filter_list_bam reproduce_results_data/XXXX --editing_DB_dir reproduce_results_data/XXXX --snp_vcf_dir reproduce_results_data/XXXX --barcode_clusters reproduce_results_data/XXXX --sname reprouce_results
~~~
