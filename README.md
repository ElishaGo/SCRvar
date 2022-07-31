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
- anaconda or miniconda

## Installation
First clone this repository - 
~~~
git clone ...
~~~
You can create a cond environment through:
~~~
conda env create -n RNA_SCvar --file RNA_SCvar.yml 
conda activate RNA_SCvar
~~~

# Mandatory and optioanl files
SCvar needs external files to run.
## Mandatory files
SCvar needs a genome file (fasta file) and a transcriptome file (gtf file).
It is important that the bam file and the genomic files are compatible.
## Optional files 
SCvar pipeline has some outputs properties available when running with some Databases files.
The file are: 
- Known editing sites (bed format)
- Known SNP (vcf format)
- Cell Barcodes filter list
- Table with cell barcodes and associated clustered (e.g. analyzed by Seurat package)

# Usage example data set

In order to reproduce the results from our research, download the files :
~~~
mkdir data/
wget https://weizmann.box.com/shared/static/qglesovsnijhcy2ote1ppk2sw01brcoe -O files.tar.gz
tar -xvzf files.tar.gz -C data/
~~~
(Elisha: tar.gz file created with the command: tar -cvzf RNA_SCvar_files.tar.gz DataBases/editing/0.editing_A_I.genecode_intersect.bed  DataBases/snp_vcf/0.snp.gencode_intersect.vcf DataBases/genecode_gtf/0.gencode.v37.annotation.gtf DataBases/raw_DB_files/genome.fa* human_9week/SRR450_SRR451_concatenated/input_files/raw_bam_merged_SRR450_SRR451.bam* human_9week/SRR452_SRR453_concatenated/input_files/merged_SRR452_SRR453.bam* 
)

# Run SCvar pipeline
In order to run the pipline the following is required:
Bam file to analyse
Genome reference fasta file
Transcriptome Genecoed file

Optional files are:
list of barcodes to use in the bam file
Editing sites file in bed format
SNP notations in vcf format
table with barcodes and their associated clusters, outputed by Seurat package

run the python script end to end:
~~~
python sc_rna_variants/helper_scripts/run_bsub_file.py PATH_TO_OUTPUT PATH_TO_BAM GENOME_REF ANNOTAION_GTF [options] 
~~~

# Reproduce our results
You can download the same files we used in our analysis and run them with the pipeline:
~~~
mkdir reproduce_results_data/
cd data_files/

# transcriptome gtf file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
gzip -d gencode.v37.annotation.gtf.gz

# reference genome
ASK WHERE THE FILE IS DOWNLOADED FROM

# editing sites database
wget http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz
gzip -d TABLE1_hg38.txt.gz

# SNP database
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz
gzip -d common_all_20180418.vcf.gz

python sc_rna_variants/helper_scripts/run_bsub_file.py PATH_TO_OUTPUT reproduce_results_data/PATH_TO_BAM reproduce_results_data/GENOME_REF reproduce_results_data/ANNOTAION_GTF --filter_list_bam reproduce_results_data/XXXX --editing_DB_dir reproduce_results_data/XXXX --snp_vcf_dir reproduce_results_data/XXXX --barcode_clusters reproduce_results_data/XXXX --sname reprouce_results
~~~
