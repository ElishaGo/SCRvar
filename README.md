# Single Cell RNA variants (SCRvar)

A script to detect RNA variations in single cell RNA-seq data.

Given a BAM file produced by cellRanger pipeline,will return a BED formatted file where each
row represents the base variations in a genomic position, and reports mutation base change, UMI counts and UMI counts per cell.

~~
# Requirements
- UNIX OS
- python3.7 or higher
- samtools, bedtools, bamtools (For conda users, these can be installed with bioconda channel)
- We recommend on using a virtual environment such as conda or python venv. 

# Installation and setup

## Create conda virtual environment
Assuming you have Anaconda or miniconda installed, create an environment with the following packages:
~~~
conda create --name scrvar --yes "python>=3.7"
conda install --name scrvar --yes -c bioconda -c conda-forge "samtools>=1.12" pysam bedtools bamtools pandas matplotlib-venn
conda activate scrvar
~~~
To make sure installation went ok, you should see at the beginning of you command line "(scrvar)", and you can check whether the packages are installed by running:
~~~
conda list
~~~
For instance, you should see the package scrvar, samtools, matplotlib-venn are within the list.
## Installation
After activating the environment, the SCRvar package can be installed via PyPi:
~~~
pip install --index-url https://test.pypi.org/simple/ SCRvar==2.0.2
~~~


# Small test sample
We provide a small subset sample to test the installation.

Download the tar.gz file:
~~~
wget ftp://dors2.weizmann.ac.il/scrvar/H1_CR_small_sample.tar.gz
~~~

Then unpack it:
~~~
mkdir sample_data/
tar -xvzf H1_CR_small_sample.tar.gz -C sample_data/
~~~
If you wish to use additional optional files, you can see under 'Reproduce our results' how to download them.


# Mandatory files
SCRvar needs external files to run:
- Bam file to analyse
- Genome reference fasta file (that was used to align the reads)
- Transcriptome annotation file, such as GENCODE gtf file

It is important that the bam file, and the other files are compatible.

See below on "Reproduce our results" section for possible files to use.

## Optional files 
SCRvar pipeline has some outputs properties available when running with additional databases files.

The file are: 
- Known editing sites (in bed format)
- Known SNP (in vcf format)
- Cell Barcodes and clusters file - a tsv/text file where the first column is cell barcodes of the sample (with no prefixex/suffixes), and the second column is clusters. The file must be without header line.


# Run SCRvar pipeline

First, activate the conda environment (If it is not activated yet)
~~~
conda activate scrvar
~~~

### 1. Preprocess the additional files:
step 0 script will do some basic preprocessing to the external files (Trancsctiptome gtf, Editing bed file and SNP vcf file if provided).
This script needs to be run *only once* per an external file. So each time you use a different file you'll need to run script 0:
~~~
mkdir external_files_processed/
process_aditional_files(step0).py external_files_processed/ANNOTATION_GTF --editing_DB_path EDITING_BD_PATH_(optional) --fasta_path GENOME_FASTA_(optional) --snp_DB_path SNP_DB_PATH_(optional)
~~~
### 2. Run SCRvar main program
steps 1-4 are the main body of the program. They will process the bam file and output a mismatch table, positions table and mid output files.

A single script execute steps 1-4, *and you need to run it once per bam file*.
~~~
mkdir SCRvar_sample_output/
run_SCRvar_steps_1234.py SCRvar_sample_output/ sample_data/possorted_genome_H1_CR_chr2_20_22.bam GENOME_FASTA external_files_processed/ANNOTATION_GTF --barcodes-cluster-file CELL_BARCODES_AND_CLUSTERS_FILE --editing_db_path EDITING_DB_PATH_(optional) --sname YOUR_SAMPLE_NAME OTHER_OPTIONAL_PARAMS
~~~

### 3. Run more analysis (steps 5, 6 (To come)) 
SCRvar has its own basic analysis outputs. 

SCRvar can also run with some customized parameters to filter out some positions, thus ending up with different analysis outputs.
 ~~~
run_filtering_and_analysis_step5.py SCRvar_sample_output/ SCRvar_sample_output/step4_aggregation_per_position_and_statistics/4.aggregated_per_position.bed --editing_db_path SCRvar_sample_output/step2_bam_gene_filter/2.editing.genecode.*_intersect.bed --snp_db_path SCRvar_sample_output/step2_bam_gene_filter/2.snp.genecode.*_intersect.vcf --sname YOUR_SAMPLE_NAME OTHER_OPTIONAL_PARAMS
 
run_gene_level_analysis_step6.py SCRvar_sample_output/ SCRvar_sample_output/step4_aggregation_per_position_and_statistics/4.aggregated_per_position.bed SCRvar_sample_output/step3_mismatch_dictionary/3.mismatch_dictionary.bed external_files_processed/ANNOTATION_GTF --barcodes-cluster-file CELL_BARCODES_AND_CLUSTERS_FILE OTHER_OPTIONAL_PARAMS
 ~~~


# Reproduce our results
You can download the same files we used in our analysis and run the steps above with the files:
~~~
mkdir reproduce_results_data/

# download the bam file
wget -P reproduce_results_data ftp://dors2.weizmann.ac.il/scrvar/H1_CR_data.tar.gz


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

# SCRvarOutputs
Step 1 - Filtering reads by alignment quality and cell barcode according to provided list.

Step 1 outputs - 
- 1.\<BAM FILE NAME>_filtered.bam and its index file.

Step 2 - Filtering bam and accessory files to contain reads within gene regions, and computing their intersections.

Step 2 outputs -  
- 2.\<SNAME>.gene_filter.bam and its index, 
- 2.editing.genecode./<SNAME>_intersect.bed, 
- 2.snp.genecode./<SNAME>_intersect.vcf,
- Venn plots of intersections with SNP and editing sites databases after gene filtering.
- 2.\<SNAME>_genomecov.txt - genome coverage output of bedtools genomcov command.
- 2.\<SNAME>_flagstat.txt - output of samtools flagstat command using the bam file as input.
- 2.\<SNAME>_gene_counts_htseq.txt - output of HTseq count command using the bam file as input.

Step 3 - Creating dictionary data structure: {cell_barcode:{(chr,coordinate, direction):{UMI: {'a': 0, 'c': 0 ,'g': 0 ,'t': 0 } } } } for mismatch and no-mismatch bases. Dictionary contains counts for both UMIs with single reads and UMIs with multi reads (more than one read).

Step 3 outputs 
- 3.mismatch_dictionary.bed, 
- 3.no_mismatch_dictionary.bed

Step 4 - Aggregation the dictionary data per genomic position, calculating mutations rate, and intersection with accessory files.

Step 4 outputs -
- 4.aggregated_per_position.bed includes flags for intersection with accessory files.

Step 5 - Filtering positions according to defined parameters and intersection plots with accessoty files.

Step 5 outputs - 
- 5.aggregated_per_position.filtered.bed 
- 5.filtering_effect - Folder with plots demonstrating mutations trends and filtering effect.
- 5.DB_intersect_effect - Folder with plots demonstraing intersection of mutations with accesoy files. 



#
#
#


#TODO add all the optioanl arguments, all parameter with underscore, rempve 'path' from parameters, add explanation on the different steps as much as you cam, in steps 5 and 6 dont need to give the snp and editin paths because you have them already, also save the filtered tables in steps 5 and 6.

