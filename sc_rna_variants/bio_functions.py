import logging
import multiprocessing
import subprocess
import copy
import os
import pathlib
import re

import pysam

import sc_rna_variants.utils as utils

# import utils


logger = logging.getLogger(__name__)


class ReadsFilter(object):
    # TODO: filter read with no CB tag as part of general read process where no filter list is suplied
    """ a class to help with reads filtering, 
    
    once the filter is initiated, running .process(pysamAlignmentRead) will
    return the read if its Kosher, or None if it isn't
    """

    def __init__(self, barcodes_set, min_mapq, cigar_clipping_allowed,
                 max_gene_length, max_no_basecall,
                 tag_for_umi, tag_for_cell_barcode,
                 threshold_low_quality=15, max_low_quality_bases=10):
        if barcodes_set is None:
            self.process = self.process_read
        else:
            self.barcodes = barcodes_set.copy()
            self.process = self.process_with_cb_filter
        self.ub_tag = tag_for_umi
        self.cb_tag = tag_for_cell_barcode
        self.mapq = min_mapq
        self.cigar_clip_limit = cigar_clipping_allowed
        self.mgl = max_gene_length
        self.max_n = max_no_basecall
        self.qual_threshold = threshold_low_quality
        self.max_low_quality = max_low_quality_bases
        # if self.mapq >= 255:
        #     self.check_mm = lambda read: read.has_tag('MM')
        # else:
        #     self.check_mm = lambda read: False

    def process_read(self, read):
        """filter reads by the different parameters"""
        if not read.has_tag(self.ub_tag) or not read.has_tag(self.cb_tag):
            return None
        elif read.reference_length > self.mgl:
            return None
        elif read.mapq < self.mapq:  # Mapq 10 and above is uniquely mapped
            return None

        # problematic cigar charachters
        # filter indels
        clipped_positions = sum([cigar_length for cigar_type, cigar_length in read.cigartuples if cigar_type in [4, 5]])
        if clipped_positions > self.cigar_clip_limit or \
                read.cigarstring.count('D') > 0 or read.cigarstring.count('I') > 0:
            return None

        # non called bases
        seq = read.query_alignment_sequence
        if seq.count('n') + seq.count('N') > self.max_n:
            return None

        # quality assurance
        low_quality_bases = sum(
            1 for read_quality in read.query_alignment_qualities if read_quality <= self.qual_threshold)
        if low_quality_bases > self.max_low_quality:
            return None

        # multiple alignment of the read - NH and MM tags filtration
        if read.has_tag('NH') and read.get_tag('NH') > 1:
            return None
        # if self.check_mm(read):   #MM tag appears only when NH > 1
        #     return None

        return read

    def process_with_cb_filter(self, read):
        """Filters reads by cell barcodes list and pipes to the main filter
        the "partition('-')[0]" command is there in case the corrected barcode contains a dash - for GEM well
        """
        if read.has_tag(self.cb_tag) and read.get_tag(self.cb_tag).partition('-')[0] in self.barcodes:
            return self.process_read(read)
        return None


class ChromosomeNamesTranslator(object):
    '''Helper class to quickly translate chromosome names between their BAM and FASTA format
    Change this class if your result neglect to include certain chromosomes'''

    def __init__(self, genome_fasta_file):
        genome_fasta = pysam.Fastafile(genome_fasta_file)
        self.references = genome_fasta.references
        genome_fasta.close()

    def translate_chromosome_name(self, bam_chromosome_name: str):
        '''Given a chromosome/segment name from the BAM file, check for the parallel name in the genomic FASTA'''
        if bam_chromosome_name in self.references:
            return bam_chromosome_name
        elif 'chr' + str(bam_chromosome_name) in self.references:
            return 'chr' + bam_chromosome_name
        if bam_chromosome_name == 'MT' and 'chrM' in self.references:
            return 'chrM'

        return None


def open_bam(input_bam, available_threads=1):
    """
    tries to open a bam file, and asserts it is sorted by coordinates
    returns pysam.libcalignmentfile.AlignmentFile
    """
    try:
        bamfile = pysam.AlignmentFile(input_bam, "rb", threads=available_threads)

        if bamfile.header['HD']['SO'] != "coordinate":
            msg = "BAM is not sorted by coordinate (missing SO:coordinate from header)"
            raise AssertionError(msg)
        return bamfile
    except Exception as e:
        msg = "\nError handeling the input bam file %s.\n%s" % (input_bam, str(e))
        raise AssertionError(msg)


def get_chrom_list(bam_input):
    """
    Creates a dictionary of sequences names (SN) and lengths (LN)
    from the header lines of the bam file
    
    Args:
        bam_input    (AlignmentFile): bamfile to work with
    
    Returns:
        chromosomes    (dictionary of str:int): list of chromosoms
    """

    chromosomes = {}

    for line in bam_input.header['SQ']:
        chromosomes[line['SN']] = line['LN']

    return chromosomes


def divide_chromosome_chunks(bam_input, max_length=100000000):
    """
    Creates a set of chromosome chunks tuple (name, start, end)
    
    Args:
        bam_input    (AlignmentFile): bamfile to work with
        max_length   (int): max length of a chunk
    
    Returns:
        chunks    (set of tuple): set of (chromosoms name, start position, end position)
    """
    chromosomes = get_chrom_list(bam_input)
    chunks = set()

    for chr_name, length in chromosomes.items():
        chunks_borders = [i for i in range(0, length, max_length)]  # devide the length by max_length
        chunks_borders.append(length)

        chunks.update([(chr_name, chunks_borders[i], chunks_borders[i + 1]) for i in range(len(chunks_borders) - 1)])

    return chunks


def divide_chromosome_by_reads(bam_input, max_reads=100000):
    """
    Creates a set of chromosome chunks tuple (name, start, end) where each chunk has a maximum of max_reads reads
    In case where a segment with length < 100 has more than max_reads,
    includes it as-is
    
    Args:
        bam_input    (AlignmentFile): bamfile to work with
        max_reads   (int): max reads per chunk
    
    Returns:
        chunks    (set of tuple): set of 4-tuples:
            (chromosoms name, start position, end position, read_count)
    """
    assert max_reads > 5000
    initial_chunks = divide_chromosome_chunks(bam_input, 200000000)

    good_size_chunks = set()
    counter = 0
    while len(initial_chunks) > 0:
        chunk = initial_chunks.pop()
        read_count = bam_input.count(*chunk)
        if read_count < max_reads or (chunk[2] - chunk[1]) < 500:
            # if we have a good amount of reads, or if we are looking at a small enough segment (super-expressed gene)
            good_size_chunks.add((*chunk, read_count))
        else:
            mid = (chunk[1] + chunk[2]) // 2
            left = (chunk[0], chunk[1], mid)
            initial_chunks.add(left)

            right = (chunk[0], mid + 1, chunk[2])
            initial_chunks.add(right)

        counter += 1
        if counter > 10 ** 6:
            raise ("Error with chromosomal segments devision")

    return good_size_chunks


def filter_chromosome_chunk(
        bamfilename: str,
        chunk_tuple: tuple,
        filtered_bam_path: pathlib.PosixPath,
        filter_object: ReadsFilter,
        sema_helper: multiprocessing.Semaphore,
        bam_header
):
    logger.debug("Working on chunk %s" % str(chunk_tuple))
    bamfile = open_bam(bamfilename)

    filtered_bam = pysam.AlignmentFile(filtered_bam_path, "wb", header=bam_header)

    ##### temporary load reducer ####
    # counter = 0
    ##########

    for read in bamfile.fetch(*chunk_tuple):
        if filter_object.process(read):
            ##### temporary load reducer ####
            # counter += 1
            # if counter%500 < 450:
            #     continue 
            ##########

            filtered_bam.write(read)

    bamfile.close()
    filtered_bam.close()

    logger.debug("Finished working on chunk %s" % str(chunk_tuple))

    sema_helper.release()  # tells multiprocessing we can start a new process


def create_filtered_bam(input_bam, filtered_barcodes_list, min_mapq, cigar_clipping_allowed, max_gene_length,
                        max_no_basecall, tag_for_umi, tag_for_cell_barcode, output_folder, threads):
    """ Given the script's argument, return the path to a filtered version of the bam"""
    # divide the genome to chromosome segments (assert no more than a 1000 exist, otherwise samtools "merge" fails)
    bamfile = open_bam(input_bam, threads)
    for i in range(20):
        chunk_max_length = 100000000 * (2 ** i)
        chunks_set = divide_chromosome_chunks(bamfile, max_length=chunk_max_length)
        if len(chunks_set) < 1000:
            break;

    temp_header_file = pathlib.Path(output_folder) / (os.path.basename(input_bam) + '_temp_header.sam')
    with open(temp_header_file, mode='w') as f:
        f.write(str(bamfile.header))

    # for technical reasons the chunked bam files will receive a lighter version of the header
    skinny_bam_header = bamfile.header.to_dict()
    skinny_bam_header.pop('CO')
    skinny_bam_header.pop('PG')
    skinny_bam_header.pop('RG')
    bamfile.close()

    # choose file names for temporary bam files
    names = utils.make_random_names_list(
        len(chunks_set), names_length=12, extension=".bam")
    temporary_bams_paths = []

    # create multiprocessing pool and helper variables
    sema = multiprocessing.Semaphore(threads)  # a shared variable, indicating number of concurrent processes
    all_processes = []
    basic_fileter = ReadsFilter(
        filtered_barcodes_list, min_mapq,
        cigar_clipping_allowed,
        max_gene_length, max_no_basecall,
        tag_for_umi, tag_for_cell_barcode)
    counter = 0

    logger.info("started creating filtered bam files for chunks of the genome (a total of %s)." % len(chunks_set))
    logger.debug("temporary file names is of the format %s_[random string].bam" % os.getpid())
    # send threads to work
    for chunk in chunks_set:
        filter_instance = copy.deepcopy(basic_fileter)
        filtered_chunk_bam_path = pathlib.Path(output_folder) / names.pop()
        temporary_bams_paths.append(str(filtered_chunk_bam_path))
        p = multiprocessing.Process(
            target=filter_chromosome_chunk,
            args=(input_bam, chunk, filtered_chunk_bam_path,
                  filter_instance, sema, skinny_bam_header)
        )
        all_processes.append(p)
        sema.acquire()  # Blocks until a semaphore is released
        p.start()

        # try and close done processes
        for p in all_processes:
            p.join(0)

        counter += 1
        if counter % 10 == 0:
            logger.info("working on chunk %d out of %d" % (counter, len(chunks_set)))

    for p in all_processes:
        p.join()

    logger.debug("finished creating filtered bam files for different segments of the genome, starting to merge them")

    # create the combined, filtered bam file
    filtered_bam_path = os.path.join(pathlib.Path(output_folder), "1_" + (os.path.basename(input_bam) + "_filtered.bam"))
    pysam.merge("-f", "-h", str(temp_header_file), filtered_bam_path, *temporary_bams_paths, '--threads', str(threads))

    # create index for the filtered file
    subprocess.run(["samtools", "index", "--threads", str(threads), filtered_bam_path], check=True)

    logger.info("finished merging filtered bams into the final, united, filtered bam file. starting cleanup")

    os.remove(temp_header_file)

    for filename in os.listdir(output_folder):
        if re.search('^%s_.*\.bam$' % os.getpid(), filename):
            os.remove(os.path.join(output_folder, filename))

    logger.debug("finished cleanup of temporary files")

    return filtered_bam_path


def calculate_position_parameters(reference_base, umis_dictionary):
    """given the reference base and the umis dictionary, calculates parameters.
    output:
        total_singles_counts,
        total_multiples_counts,
        [same multi, transition multi, reverse multi, transversion multi,
        same single, transition single, reverse single, transversion single,
        mixed, non_reference_to_total_percent]
    # TODO: update documentation
    """

    mixed = 0
    singles = {'a': 0, 'c': 0, 'g': 0, 't': 0}  # base: number_of_umis
    multiples = {'a': 0, 'c': 0, 'g': 0, 't': 0}  # base: number_of_umis

    # sum up things
    for bases_dictionary in umis_dictionary.values():
        # get a list of [(base, appearances)]
        appearances = [item for item in bases_dictionary.items() if item[1] > 0]

        if len(appearances) > 1:  # more than one base for this UMI -> mixed
            mixed += 1
            continue

        read_base, number_of_reads = appearances[0]

        if number_of_reads > 1:
            multiples[read_base] += 1
        else:
            singles[read_base] += 1

    singles_anotated = [None, None, None, None]
    singles_anotated[0] = singles[reference_base]  # same
    singles_anotated[1] = singles['acgt'['gtac'.index(reference_base)]]  # transition
    singles_anotated[2] = singles['acgt'['tgca'.index(reference_base)]]  # reverse_complement
    singles_anotated[3] = singles['acgt'['catg'.index(reference_base)]]  # transversion

    # copy to multiples
    multiples_anotated = [None, None, None, None]
    multiples_anotated[0] = multiples[reference_base]
    multiples_anotated[1] = multiples['acgt'['gtac'.index(reference_base)]]
    multiples_anotated[2] = multiples['acgt'['tgca'.index(reference_base)]]
    multiples_anotated[3] = multiples['acgt'['catg'.index(reference_base)]]

    total_singles_counts = sum(singles_anotated)
    total_multiples_counts = sum(multiples_anotated)
    non_reference_umi_counts = total_singles_counts - singles_anotated[0] + \
                               total_multiples_counts - multiples_anotated[0]

    # if non_reference_umi_counts > 0:    # if we have more than 1 umi different from the reference
    #     result_list = [multiples['a'], multiples['c'], multiples['g'], multiples['t'],
    #             singles['a'], singles['c'], singles['g'], singles['t'],
    #             mixed, int(100*non_reference_umi_counts/len(umis_dictionary))%101]  # seems that %101 is not needed
    #     return total_singles_counts, total_multiples_counts, [str(number) for number in result_list]

    if non_reference_umi_counts > 0:  # if we have more than 1 umi different from the reference
        result_list = [multiples_anotated[0], multiples_anotated[1], multiples_anotated[2], multiples_anotated[3],
                       singles_anotated[0], singles_anotated[1], singles_anotated[2], singles_anotated[3],
                       mixed, int(
                100 * non_reference_umi_counts / len(umis_dictionary)) % 101]  # seems that %101 is not needed
        return total_singles_counts, total_multiples_counts, [str(number) for number in result_list]
    else:
        return total_singles_counts, total_multiples_counts, []


def process_chromosome_chunk(
        bamfilename: str,
        chunk_tuple: tuple,
        fasta_chromosome_name: str,
        fasta_filename: str,
        ub_tag: str,
        cb_tag: str,
        output_csv_path: pathlib.PosixPath,
        min_position_quality: int = 20
):
    '''Process the chromosome chunk, write content to csv formatted file
        
    BED format with the following columns, each row represents a single cell in a single position where at least some of the UMIs were mutated
    1. Chromosome
    2. Start position
    3. End position
    4. Cell-Barcode & "-" 		# do we also need sample barcode?
    5. Percentage of non-reference UMI to total UMIs
    6. Strand			# always . To indicate no strand
    7. Reference nucleotide
    8. Same as reference - UMIs with > 1 read
    9. Point mutation transition - UMIs with > 1 read
    10. Complementary - UMIs with > 1 read
    11. Point mutation transversion - UMIs with > 1 read
    12. Same as reference - UMIs with 1 read
    13. Point mutation transition - UMIs with 1 read
    14. Complementary - UMIs with 1 read
    15. Point mutation transversion - UMIs with 1 read
    16. UMIs with unmatching (mixed) reads

    # 8. Mutations R->A - UMIs with > 1 read
    # 9. Mutations R->C - UMIs with > 1 read
    # 10. Mutations R->G - UMIs with > 1 read
    # 11. Mutations R->T - UMIs with > 1 read
    # 12. Mutations R->A - UMIs with 1 read
    # 13. Mutations R->C - UMIs with 1 read
    # 14. Mutations R->G - UMIs with 1 read
    # 15. Mutations R->T - UMIs with 1 read
    #  TODO: update documentation
    '''

    logger.debug("Working on chunk %s:%d-%d, with %d reads" % chunk_tuple)

    # open bam and start looping over pileup columns
    bamfile = open_bam(bamfilename)

    # create basic data structure:
    # cells_dict is a dictionary with cell barcodes as keys,
    # and values are:
    # positions dictionaries have UMIs as keys, and bases dictionary as values
    cells_dict = {}
    start_position = chunk_tuple[1]
    end_position = chunk_tuple[2]

    for read in bamfile.fetch(chunk_tuple[0], start_position, end_position):
        cell_barcode = read.get_tag(cb_tag)
        umi = read.get_tag(ub_tag)
        alignedRefPositions = read.get_reference_positions()
        readSequence = read.query_alignment_sequence
        direction = '-' if read.is_reverse else '+'

        # make sure the cell has an entry in the general structure
        if cell_barcode not in cells_dict:
            cells_dict[cell_barcode] = {}

        # Go over positions, count only aligned
        for i in range(len(readSequence)):
            position = alignedRefPositions[i]
            # note that at this point the position might be outside chunk

            # checking quality of base
            if (read.query_alignment_qualities[i] < min_position_quality):
                continue
            base = read.query_alignment_sequence[i].lower()
            if base == 'n':
                continue

            if (direction, position) not in cells_dict[cell_barcode]:
                cells_dict[cell_barcode][(direction, position)] = {}

            if umi not in cells_dict[cell_barcode][(direction, position)]:
                cells_dict[cell_barcode][(direction, position)][umi] = {'a': 0, 'c': 0, 'g': 0, 't': 0}

            cells_dict[cell_barcode][(direction, position)][umi][base] += 1

    bamfile.close()

    # get reference nucleotides
    genome_fasta = pysam.FastaFile(fasta_filename)
    reference = genome_fasta.fetch(fasta_chromosome_name, start_position, end_position + 1)
    genome_fasta.close()

    temp_results = []
    has_mutation = set()
    total_umis_for_unmutated_positions = {}

    # process data for position
    for cell in cells_dict:
        for (direction, position) in cells_dict[cell]:
            if position < start_position or position > end_position:
                continue

            ref_base = reference[position - start_position].lower()
            if ref_base == 'n':
                continue

            total_singles_counts, total_multiples_counts, params = calculate_position_parameters(ref_base,
                                                                                                 cells_dict[cell][(
                                                                                                     direction,
                                                                                                     position)])

            if params:  # has mutation
                if direction == '-':  # flip for negative strand
                    # we only flip the reference, the mutation stays the same
                    # params = [params[i] for i in [3,2,1,0,7,6,5,4,8,9]]
                    ref_base = 'acgt'[3 - 'acgt'.index(ref_base)]

                has_mutation.add((direction, position))
                # ouput bed file as positions start (1-based), end (1-based)
                temp_results.append('\t'.join([
                    fasta_chromosome_name, str(position + 1), str(position + 2), cell,  # CHANGED HERE 15/12/2021
                    params[-1], direction, ref_base, *params[:-1]
                ]) + '\n')
            else:  # no mutations
                if (direction, position) not in total_umis_for_unmutated_positions:
                    # singles umis, multiples umis & unique cells
                    total_umis_for_unmutated_positions[(direction, position)] = [0, 0, 0]
                total_umis_for_unmutated_positions[(direction, position)][0] += total_multiples_counts
                total_umis_for_unmutated_positions[(direction, position)][1] += total_singles_counts
                total_umis_for_unmutated_positions[(direction, position)][2] += 1

    # logger.debug("has_mutation length is %d, total_umis_for_unmutated_positions had %d items" % (len(has_mutation), len(total_umis_for_unmutated_positions)))
    with open(output_csv_path.with_suffix('.tsv2'), 'w') as fh:
        fh.writelines(
            ['\t'.join(
                # chromosome, start (1-based), end (1-based), unique cells, direction, multiples, singles
                [fasta_chromosome_name, str(key[1] + 1), str(key[1] + 2), str(value[2]), key[0], str(value[0]),
                 # CHANGED HERE 15/12/2021
                 str(value[1])]
            ) + '\n' for key, value in total_umis_for_unmutated_positions.items() if key in has_mutation]
        )

    with open(output_csv_path, 'w') as fh:
        fh.writelines(temp_results)

    del temp_results
    del cells_dict
    del reference

    logger.debug("finished chunk %s:%d-%d" % chunk_tuple[:3])


def variants_finder(filtered_bam, genome_fasta, tag_for_umi, tag_for_cell_barcode, output_folder, threads):
    """ Given the filtered bam and script's argument, creates the final csv"""
    names_translator = ChromosomeNamesTranslator(genome_fasta)

    # devide the genome to chromosome chunks with limited amount of reads
    bamfile = open_bam(filtered_bam, threads)
    logger.debug("Starting to divide the genome to equaly covered segments, this might take a while")
    # Chunk tuple structure: (chromosoms name, start position, end position, read_count)
    chunks_set = divide_chromosome_by_reads(bamfile, max_reads=200000)
    bamfile.close()

    # assert only chromosomes that appear in the FASTA file are taken
    chunks_set = {chunk for chunk in chunks_set if names_translator.translate_chromosome_name(chunk[0])}

    logger.debug("Genome was devided into %d segments" % len(chunks_set))

    # choose file names for temporary bam files
    names = utils.make_random_names_list(
        len(chunks_set), names_length=12, extension=".tsv")
    temporary_tsvs_paths = []

    # create multiprocessing pool and helper variables
    pool = multiprocessing.Pool(threads, maxtasksperchild=1)
    asyncs = []

    logger.info("starting to process different segments of the genome.")
    logger.debug("temporary file names is of the format %s_[random string].tsv" % os.getpid())

    # send threads to work
    for chunk in chunks_set:
        filtered_chunk_tsv_path = pathlib.Path(output_folder) / names.pop()
        temporary_tsvs_paths.append(str(filtered_chunk_tsv_path))
        ######## DEBUGGING
        # process_chromosome_chunk(
        #         filtered_bam, chunk,
        #         names_translator.translate_chromosome_name(chunk[0]),
        #         arguments.genome_fasta, arguments.tag_for_umi,
        #         arguments.tag_for_cell_barcode, filtered_chunk_tsv_path
        #         )
        #######
        result = pool.apply_async(func=process_chromosome_chunk, args=(
            filtered_bam, chunk,
            names_translator.translate_chromosome_name(chunk[0]),
            genome_fasta, tag_for_umi,
            tag_for_cell_barcode, filtered_chunk_tsv_path
        ))
        asyncs.append(result)

    for i in range(len(asyncs)):
        asyncs[i].get()
        if (i) % 10 == 0:
            logger.debug("completed %d chunks out of %d" % (i + 1, len(asyncs)))
    pool.close()
    pool.join()

    logger.debug("finished creating temporary tsv files for different segments of the genome, starting to merge them")

    # create the combined, filtered bam file
    output_tsv = pathlib.Path(output_folder) / '3_mismatch_dictionary.bed6'
    output_unmutated_tsv = pathlib.Path(output_folder) / '3_no_mismatch_dictionary.bed6'

    header_line = '\t'.join([
        "chrom", "chromStart", "chromEnd", "cell barcode",
        "percent of non ref", "strand", "reference base",
        "same multi reads", "transition multi reads", "reverse multi reads", "transvertion multi reads",
        "same single reads", "transition single reads", "reverse single reads", "transvertion single reads",
        "mixed reads"]) + '\n'

    header_line_unmutated = '\t'.join([
        "chrom", "chromStart", "chromEnd", "count of unmutated cell barcodes", "strand", "unmutated multi reads",
        "unmutated single reads"
    ]) + '\n'

    with open(output_tsv, 'w') as outfile, open(output_unmutated_tsv, 'w') as outfile_unmutated:
        outfile.write(header_line)
        outfile_unmutated.write(header_line_unmutated)

        for fpath in temporary_tsvs_paths:
            with open(fpath) as infile, open(fpath + '2') as infile2:
                outfile.write(infile.read())
                infile.close()
                os.remove(fpath)

                outfile_unmutated.write(infile2.read())
                infile2.close()
                os.remove(fpath + '2')

    logger.debug("finished merging and cleanup of tsv files")

    return
