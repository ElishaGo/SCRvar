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