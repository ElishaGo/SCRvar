# Single Cell RNA variants (SCvar)

A script to locate cases of RNA modifications in single cell RNAseq data.

Given a BAM file with cellRanger format, returns a BED formatted file where each
row represents the modifications a certain cell has in a certain position,
categorised by mutation type and UMI counts.

~~
# Requirements
- UNIX OS
- python3.7 or higher
- samtools, bedtools, bamtools (For conda users, these can be installed with bioconda channel)
- We recommend on using a virtual environment such as conda or python venv. 

## Create conda virtual environment
Assuming you have Anaconda or miniconda installed, create an environment with the following packages:
~~~
conda create --name scvar_env python=3.7
conda activate scvar_env
conda install -c bioconda samtools bedtools bamtools
~~~

## Setup
After activating the environment, the SCvar package can be installed via PyPi:
~~~
pip install pandas matplotlib-venn (TEMP - while using test.pypi)
pip install --index-url https://test.pypi.org/simple/ scRNAvariants==2.0.25
~~~

## small test sample
We provide a small subset sample for those who first which to run the program and see its outputs.
First download the tar.gz file
~~~
wget PATH_TO_DORS2
TEMP - download the file from box https://weizmann.box.com/s/651v8e8p6z6b0u6ru0p058j9r7mip5kr 
~~~

Then unpack it:
~~~
mkdir sample_data/
tar -xvzf H1_CR_small_sample.tar.gz -C sample_data/
~~~
If you wish to use additional optional files, you can see under 'Reproduce our results' how to download them.


# Mandatory files
SCvar needs external files to run:
- Bam file to analyse
- Genome reference fasta file
- Transcriptome Genecode file

It is important that the bam file, and the other files are compatible.

## Optional files 
SCvar pipeline has some outputs properties available when running with additional databases files.

The file are: 
- Known editing sites (in bed format)
- Known SNP (in vcf format)
- Cell Barcodes and clusters file - a tsv/text file where the first column is cell barcodes of the sample (with no prefixex/suffixes), and the second column is clusters. The file must be without header line.


# Run SCvar pipeline
### 1. Preprocess the additional files:
step 0 script will do some basic preprocessing to the external files (Trancsctiptome gtf, Editing bed file and SNP vcf file if provided).
This script needs to be run *only once* per an external file. So each time you use a different file you'll need to run script 0:
~~~
mkdir external_files_processed/
process_aditional_files(step0).py external_files_processed/ANNOTATION_GTF --editing_DB_path EDITING_BD_PATH_(optional) --fasta_path GENOME_FASTA_(optional) --snp_DB_path SNP_DB_PATH_(optional)
~~~
### 2. Run SCvar main program
steps 1-4 are the main body of the program. They will process the bam file and output a mismatch table, positions table and mid output files.

A single script execute steps 1-4, *and you need to run it once per bam file*.
~~~
mkdir SCvar_sample_output/
run_SCvar_(steps_1234).py SCvar_sample_output/ sample_data/possorted_genome_H1_CR_chr2_20_22.bam GENOME_FASTA external_files_processed/ANNOTATION_GTF --barcodes-cluster-file CELL_BARCODES_AND_CLUSTERS_FILE --editing_db_path EDITING_DB_PATH_(optional) --sname YOUR_SAMPLE_NAME OTHER_OPTIONAL_PARAMS
~~~

### 3. Run more analysis (steps 5, 6 (To come)) 
SCvar has its own basic analysis outputs. 

SCvar can also run with some customized parameters to filter out some positions, thus ending up with different analysis outputs.
 ~~~
run_filtering_and_analysis_(step5).py SCvar_sample_output/ SCvar_sample_output/step4_aggregation_per_position_and_statistics/4.aggregated_per_position.bed --editing_db_path SCvar_sample_output/step2_bam_gene_filter/2.editing.genecode.*_intersect.bed --snp_db_path SCvar_sample_output/step2_bam_gene_filter/2.snp.genecode.*_intersect.bed --sname YOUR_SAMPLE_NAME OTHER_OPTIONAL_PARAMS
 
run_gene_level_analysis_(step6).py SCvar_sample_output/ SCvar_sample_output/step4_aggregation_per_position_and_statistics/4.aggregated_per_position.bed SCvar_sample_output/step3_mismatch_dictionary/3.mismatch_dictionary.bed external_files_processed/ANNOTATION_GTF --barcodes-cluster-file CELL_BARCODES_AND_CLUSTERS_FILE OTHER_OPTIONAL_PARAMS
 ~~~


# Reproduce our results
You can download the same files we used in our analysis and run the steps above with the files:
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
~~~
