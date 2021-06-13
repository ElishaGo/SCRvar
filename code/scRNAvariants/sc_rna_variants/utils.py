import os
import logging
import argparse
import re
import random
import string

import pysam

logger = logging.getLogger(__name__)


## Arguments parsing functions
# according to argparse reference, these functions return the values
 # to be stored in the arguments namespace
def assert_is_file(string):
    '''Make sure the string represents an existing file in the OS'''
    if not os.path.isfile(string):
        msg="\n%s is not a valid file" % string
        raise argparse.ArgumentTypeError(msg)    
    return string


def assert_is_directory(string):
    '''Make sure the given string indicates an OS directory'''
    if not os.path.isdir(string):
        msg="\n%s is not a valid directory" % string
        raise argparse.ArgumentTypeError(msg)    
    return string


def bam_check(input_bam):
    '''Make sure the given string represents an existing BAM file in the OS
    if the file isn't valid according to samtools, raises error
    if the file is not indexed, tries to index it (raises error if fails)'''
    assert_is_file(input_bam)
    try:
        pysam.quickcheck(input_bam)
    except Exception as e:
            msg="\nSamtools did not recognise the input bam file %s as valid.\n%s" % (input_bam,str(e))
            raise argparse.ArgumentTypeError(msg)
            
    if not os.path.exists(input_bam+'.bai'):
        try:
            pysam.index(input_bam)
        except Exception as e:
            msg="\nSamtools cannot index the input bam file %s.\n%s" % (input_bam,str(e))
            raise argparse.ArgumentTypeError(msg)
    return input_bam


def genome_fasta_check(input_genome_fasta):
    '''Make sure the given string represents an existing genome FASTA file in the OS
    if the file is not indexed, tries to index it (automatically with pysam.Fastafile)'''
    assert_is_file(input_genome_fasta)
    
    try:
        genome_fasta = pysam.Fastafile(input_genome_fasta)
        number_of_reference_sequences = genome_fasta.nreferences
        if not genome_fasta.closed:
            genome_fasta.close()
    
        if number_of_reference_sequences < 1:
            msg="\nNo reference sequences were found in the genome FASTA file %s" % input_genome_fasta
            raise argparse.ArgumentTypeError(msg)
    except Exception as e:
        msg="\nError processing the genome FASTA file %s.\n%s" % (input_genome_fasta,str(e))
        raise argparse.ArgumentTypeError(msg)

    return input_genome_fasta


def filtered_barcodes_processing(filtered_barcodes_file):
    """Given a path to a file with barcodes as a first column, returns a set of the barcodes"""
    assert_is_file(filtered_barcodes_file)
    try:
        with open(filtered_barcodes_file, 'r') as barcodes_file_handler:
            filtered_barcodes = set()
            for line in barcodes_file_handler:
                split_line = re.split('\W',line.strip())
                assert len(re.findall('[^aAcCgGtT\-\d]',split_line[0])) == 0,\
                    'Barcode containing non-nucleotide letters at line %d' % (1 + len(filtered_barcodes))
                without_GEM_well_number =split_line[0].partition("-")[0]
                filtered_barcodes.add(without_GEM_well_number)
        return filtered_barcodes
    except Exception as e:
        msg="\nCannot process cell barcodes list in %s. This might be a formating issue;\n%s" % (filtered_barcodes_file,str(e))
        raise argparse.ArgumentTypeError(msg) 


class ArgparserFormater(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """
    a helper class used in formating the -h help text. Inherts from the argparse classes:
    RawTextHelpFormatter - can break lines in the help text, but don't print default values
    ArgumentDefaultsHelpFormatter - print default values, but don't break lines in the help text
    """
    pass


# general helper functions
def make_random_names_list(list_length, names_length=12, extension=".txt"):
    """
    Creates a list of random names starting with the process id
    
    Args:
        list_length   (int): amount of names to create
        names_length  (int): length of the names beyond pid (default=12)
        extension     (str): file extension (default=".txt")
    
    Returns:
        chrom_set    (list of str): list of chromosoms
    """
    names = []
    pid = os.getpid()
    for i in range (list_length):
        rand_name = ''.join(random.choices(string.ascii_letters, k=names_length))
        names.append("%d_%s%s" % (pid, rand_name, extension))
    return names


    
    
    
    
