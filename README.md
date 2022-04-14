# Single Cell RNA variants

A script to locate cases of RNA modifications in single cell RNAseq data.

Given a BAM file with cellRanger format, returns a BED formated file where each row represents the modifications a certain cell has in a certain position, categorised by mutation type and UMI counts.

requirements: pysam package


conda environment setup example:
```
conda create -n pysam37
conda activate pysam37
conda config --env --add channels conda-forge
conda config --env --add channels bioconda
conda install python=3.7 pysam samtools
```


## Parameters per step.
## TODO: This section should be based on the output of XXX.py -h, which is automatticaly created by argparser

Input:
-- tsv file
-- filter out positions with less than X barcodes

Output:
-- statistics plots PDF/html - before and after filtering
-- tsv file with aggregated cell barcodes (concatenate the different cell barcodes into one column):
Chromosome
start
end
Count of cell barcode(s)
percent of non ref (recalculate - unique UMIs that are different from the reference divided by unique UMIs)
strand
reference base
R->A multi		(aggregate)
R->T multi		(aggregate)
R->G multi	(aggregate)
R->C multi	(aggregate)
R->A single	(aggregate)
R->T single	(aggregate)
R->G single	(aggregate)
R->C single	(aggregate)
mixed reads	(aggregate)
Cell barcodes ('abc-1','def-1')



Example input:
Chr5	112000	abc	1	1	4
Chr5	112000	def	2	0	1

Example output:
Chr5	112000	2	3	1	5	abc,def
