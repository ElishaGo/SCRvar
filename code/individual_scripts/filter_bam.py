import os
import argparse
import pysam

def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser()
    # positional arguments
    parser.add_argument('input_bam', help='bam file to filter.')

    parser.add_argument('filter_list', help='list of cell barcodes to keep in the output bam file')

    # optional arguments
    parser.add_argument('--output_folder', default= os.path.join(os.getcwd(), 'TEST_bam_filter_outputs'), help='Path to the directory where the output files will be placed')

    parser.add_argument('--name_suffix',type = str, default= "", help='Add suffix to output files')
    
    parser.add_argument('--threads', type=int,default=1)
    
    return parser.parse_args(arguments)
    
    
def make_output_folder(output_path):
    '''Make sure the string represents an existing file in the OS.
    If not, create a default file'''
    if not os.path.isdir(output_path):
        print("\n creating a new output folder in: {out}".format(out=output_path))
        os.mkdir(output_path)


def clear_filter_list(filter_list_path):
    '''Check for duplicates and optional problematic elements in list'''
    #save filter list with end signs to object 
    filter_list = os.popen("cat -E {}".format(filter_list_path)).read().split('\n')
    #remove duplicates
    filter_list = list(set(filter_list))
    #check if end sign is by it's own, and correct
    for i, cb in enumerate(filter_list):
        if cb in ['$','',' ','\t','\n']:
            filter_list.pop(i)
        #remove space from barcodes
        filter_list[i] = cb.replace(' ','')

    
    #save the clean list
    clean_list_path = filter_list_path.split("/")
    clean_list_path[-1] = "temp_" + clean_list_path[-1]
    clean_list_path = "/".join(clean_list_path)
    with open(clean_list_path, "w") as f:
        f.writelines("%s\n" % cb for cb in filter_list)
        
    return clean_list_path
        
        
def filter_bam(args, filter_list):
    '''Recieve bam file and filter list, and creates a new bam file with only the cell barcodes from the filter list'''
    bam_file = args.input_bam
    out_folder = args.output_folder
    suffix = args.name_suffix
    threads = args.threads
    
#     os.system("export name={suffix};export BAM_FILE={bam_file};export FILTER_FILE={filter_list}; samtools view -H $BAM_FILE > {out}/SAM_header_$name; samtools view -@ {threads} $BAM_FILE | LC_ALL=C grep -F -f $FILTER_FILE | cat {out}/SAM_header_$name - | samtools view -@ {threads} -b - > {out}/CBfiltered_$name.bam;samtools index {out}/CBfiltered_$name.bam;rm SAM_header_$name;rm $FILTER_FILE".format(bam_file=bam_file, filter_list=filter_list, out= out_folder,suffix=suffix, threads = threads))

    os.system("samtools view -H {bam_file} > {out}/{name}_SAM_header; samtools view -@ {threads} {bam_file} | LC_ALL=C grep -F -f {filter_list} | cat {out}/{name}_SAM_header - | samtools view -@ {threads} -b - > {out}/{name}_CBfiltered.bam;samtools index {out}/{name}_CBfiltered.bam;rm {name}_SAM_header;rm {filter_list}".format(bam_file=bam_file, filter_list=filter_list, out= out_folder,name=suffix, threads = threads))

def open_bam(input_bam, available_threads=1):
    """
    tries to open a bam file, and asserts it is sorted by coordinates
    returns pysam.libcalignmentfile.AlignmentFile
    """
    try:
        bamfile = pysam.AlignmentFile(input_bam, "rb", threads = available_threads)
        
        if bamfile.header['HD']['SO'] != "coordinate":
            msg="BAM is not sorted by coordinate (missing SO:coordinate from header)"
            raise AssertionError(msg)
        return bamfile
    except Exception as e:
        msg="\nError handeling the input bam file %s.\n%s" % (input_bam,str(e))
        raise AssertionError(msg)
        
        
if __name__ == '__main__':
    os.system("echo STARTING filter_bam program")
    args = parse_arguments()
    open_bam(args.input_bam)
    make_output_folder(args.output_folder)
    filter_list = clear_filter_list(args.filter_list)
    filter_bam(args, filter_list)
    os.system("echo program finished")
