

import pathlib
import csv
import gzip
import Bio.SeqIO.QualityIO
import Bio.SeqIO.FastaIO as FastaIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import argparse
import collections
import logging
import subprocess
import sys
import re
import itertools
import urllib.request
import urllib.response
import tempfile




'''
Start Global config section
'''

MTAGS_CLASSIFY_MIN_PERCID = 97.0  # the minimum perc id of an alignment to be reported
MTAGS_CLASSIFY_MIN_BITSCORE = 0.0  # the minimum bitscore of an alignment to be considered
MTAGS_CLASSIFY_MAX_ALIGNMENTS_PER_QUERY = 5000  # maximum number of alignments for each query
MTAGS_CLASSIFY_TOP_PERCENT = 95.0  # consider only alignments with that range in percentage of best bitscore of alignment
MTAGS_CLASSIFY_MIN_ALIGNMENT_LENGTH = 100  # Minimum length of alignment. In paired end mode both reads have to be longer then this number combined.
MTAGS_CLASSIFY_QCOV = 0.5
MTAGS_CLASSIFY_MATCHSCORE = 2
MTAGS_CLASSIFY_MISMATCHSCORE = -4

workdir = pathlib.Path(__file__).absolute().parent.parent.joinpath('db')
database = workdir.joinpath('SILVA-138_NR-97_complink_cons.vsearch.udb')
taxmap = workdir.joinpath('SILVA-138_NR-97_complink_cons.taxmap')
db_marker_file = workdir.joinpath('download.done')



R1_SUFFIXES = set(['.1', 'R1', '/1'])
R2_SUFFIXES = set(['.2', 'R2', '/2'])


#ALLOWED_TAXONOMIC_RANKS = ['root', 'domain', 'superkingdom', 'kingdom', 'subkingdom', 'infrakingdom', 'superphylum', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'superorder', 'order', 'suborder', 'superfamily', 'family', 'subfamily', 'genus', 'otu'][::-1]
REPORT_TAXONOMIC_RANKS = ['root', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'otu'][::-1]
REPORT_TAXONOMIC_INDEXES = set([0, 1, 7, 10, 13, 16, 18, 19])


TAXID_UNALIGNED = 'Unaligned'
TAXID_UNASSIGNED = 'Unassigned'
TAXID_ASSIGNED_TO_HIGHER_RANK = 'ASSIGNED_HIGHER_RANK'
RANK_2_NAME_SEPARATOR = '__'

'''

End Global config section
'''

vsearch_aln = collections.namedtuple('vsearch_aln', 'target qcov percid qstart qend rstart rend alnlength matches mismatches gaps, strand')





class AlignedInsert():
    __slots__ = ['insert_name', 'is_paired','reference_2_alignments']

    def __init__(self, insert_name, is_paired):
        self.insert_name = insert_name
        self.reference_2_alignments = collections.defaultdict(list)
        self.is_paired = is_paired
    def add_alignment(self, alignment, orientation=1):
        self.reference_2_alignments[alignment.target].append((orientation, alignment))


    def extract_best_alignments(self):
        '''
        Best alignment is defined by
        1. Having the best score (additive for PE alignments)
        2. Of those, having the shortest alignment length
        3. Of those check if the alignment length is longer then MITAG_CLASSIFY_MIN_ALIGNMENT_LENGTH

        :return: a list with alignments [(reference, score, length), ...] or and empty list of this insert didnt align or was too short
        '''

        # TODO check the orientation of fw and reverse reads and calc common percid. E.g percidfw + percidrev / 2
        alignments_scored = []
        for reference, alignments in self.reference_2_alignments.items():
            if len(alignments) > 2:
                fw_alns = list(filter(lambda aln: aln[0] == '1', alignments))
                rev_alns = list(filter(lambda aln: aln[0] == '2', alignments))
                if len(fw_alns) > 1:
                    max_percid = max(map(lambda aln: aln[1].percid, fw_alns))
                    fw_alns = [list(filter(lambda aln: aln[1].percid == max_percid, fw_alns))[0]]
                if len(rev_alns) > 1:
                    max_percid = max(map(lambda aln: aln[1].percid, rev_alns))
                    rev_alns = [list(filter(lambda aln: aln[1].percid == max_percid, rev_alns))[0]]
                alignments = fw_alns + rev_alns

            assert len(alignments)  <= 2
            score = 0
            alnlength = 0
            for (orientation, alignment) in alignments:
                score = score + MTAGS_CLASSIFY_MATCHSCORE * alignment.matches + MTAGS_CLASSIFY_MISMATCHSCORE * (alignment.mismatches + alignment.gaps)
                alnlength = alnlength + alignment.matches + alignment.mismatches + alignment.gaps
            alignments_scored.append((reference, score, alnlength))

        best_score = max(map(lambda alignment: alignment[1], alignments_scored))
        best_score_alignments = list(filter(lambda alignment: alignment[1] == best_score, alignments_scored))
        shortest_alignment_length = min(map(lambda alignment: alignment[2], best_score_alignments))
        best_score_and_shortest_alignments = list(filter(lambda alignment: alignment[2] == shortest_alignment_length, best_score_alignments))
        if shortest_alignment_length < MTAGS_CLASSIFY_MIN_ALIGNMENT_LENGTH:
            return []
        else:
            return best_score_and_shortest_alignments




class FastA(object):
    '''
    Standard data container for fasta sequences
    '''
    __slots__ = ['header', 'sequence']
    def __init__(self, header: str, sequence: str) -> None:
        self.header = header
        self.sequence = sequence

def revcomp(sequence: str) -> str:
    '''
    Reverse complement a standard nucleotide sequence.
    :param sequence:
    :return:
    '''
    return str(Seq(sequence, generic_dna).reverse_complement())

def stream_fa(sequence_file: pathlib.Path):
    '''
    Read a fastq file either gzipped or not and return it as a stream of tuples
    (Header, Sequence, Quality)
    :param infile:
    :return: Generator[FastA, None, None]
    '''
    sequence_file = str(sequence_file)
    if sequence_file.endswith('fq.gz') or sequence_file.endswith('fastq.gz'):
        with gzip.open(sequence_file, 'rt') as handle:
            for header, sequence, qual in Bio.SeqIO.QualityIO.FastqGeneralIterator(handle):
                yield FastA(header, sequence)
    elif sequence_file.endswith('fq') or sequence_file.endswith('fastq'):
        with open(sequence_file) as handle:
            for header, sequence, qual in Bio.SeqIO.QualityIO.FastqGeneralIterator(handle):
                yield FastA(header, sequence)
    elif sequence_file.endswith('fasta.gz') or sequence_file.endswith('fa.gz'):
        with gzip.open(sequence_file, 'rt') as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                yield FastA(header, sequence)
    elif sequence_file.endswith('fasta') or sequence_file.endswith('fa'):
        with open(sequence_file) as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                yield FastA(header, sequence)
    else:
        raise Exception(f'{sequence_file} not a sequence file.')


def startup():
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)
    logging.info('Starting mTAGs')


def shutdown(status=0):
    logging.info(f'Finishing mTAGs with status:\t{status}')
    sys.exit(status)

def _i_read_blasttab_file(file):
    """
    Read and yield an blasttab formatted file.
    if the file ends with *.gz read as gzip
    :param file:
    :return:
    """
    if not file.exists():
        raise FileNotFoundError('File {} not found'.format(file))
    if str(file).endswith('.gz'):
        with gzip.open(file, "rt") as handle:
            tsvreader = csv.reader(handle, delimiter="\t")
            for line in tsvreader:
                yield line
    else:
        with open(file, 'r', 1024 * 1024 * 64) as handle: # 64MB buffer
            tsvreader = csv.reader(handle, delimiter="\t")
            for line in tsvreader:
                yield line





def load_and_parse_vsearch_tab_file(file):
    """
    Reads and parses a vsearch tab file. Only the specific vsearch tab format that is defined by
    query+target+qcov+id+qilo+qihi+tilo+tihi+alnlen+ids+mism+gaps+qstrand
    can be parsed.
    Additions to the standard format will lead to unpredictable behavior.
    :param file:
    :return: A dict with read 2 alignment mappings
    """
    if not file.exists():
        raise FileNotFoundError('File {} not found'.format(file))
        shutdown(1)
    readname_2_alignments = collections.defaultdict(list)
    for line in _i_read_blasttab_file(file):
        query = line[0]
        target = line[1]
        qcoverage = float(line[2])
        percentid = float(line[3])
        qstart = int(line[4])
        qend = int(line[5])
        refstart = int(line[6])
        refend = int(line[7])
        alignmentlength = int(line[8])
        matches = int(line[9])
        mismatches = int(line[10])
        gaps = int(line[11])
        strand = line[12]
        yield (query, vsearch_aln(target=target, qcov=qcoverage, percid=percentid, qstart=qstart, qend=qend, rstart=refstart, rend=refend, alnlength=alignmentlength, matches=matches, mismatches=mismatches, gaps=gaps, strand=strand))



def read_and_pair_vsearch_tab_files(pe_m8_file, se_m8_file):
    aligned_inserts = {}
    pe_readname_2_alignments = load_and_parse_vsearch_tab_file(pe_m8_file)
    se_readname_2_alignments = load_and_parse_vsearch_tab_file(se_m8_file)
    for readname, se_aln in se_readname_2_alignments:
        aligned_insert = aligned_inserts.get(readname, None)
        if not aligned_insert:
            aligned_insert = AlignedInsert(readname, False)
            aligned_inserts[readname] = aligned_insert

        aligned_inserts[readname].add_alignment(se_aln, '1')
    for readname, pe_aln in pe_readname_2_alignments:
        [insert_name, orientation] = readname.rsplit('.',1)
        #orientation = readname[-1]

        aligned_insert = aligned_inserts.get(insert_name, None)
        if not aligned_insert:
            aligned_insert = AlignedInsert(insert_name, True)
            aligned_inserts[insert_name] = aligned_insert
        aligned_inserts[insert_name].add_alignment(pe_aln, orientation)

    return aligned_inserts








def lca_vsearch(best_alignments, taxid_2_path):
    """
    For every read find the LCA using the taxonomic path in the taxid2map data
    1. If a read maps uniquely to a reference --> This references is assigned to the read
    2. If a read maps to multiple references --> find the LCA of the references using LTREE LCP (longest common prefix) algorithm

    """
    if len(best_alignments) == 0:
        return taxid_2_path[TAXID_UNALIGNED]
    if len(best_alignments) == 1:
        return taxid_2_path[best_alignments[0][0]]

    references = [aln[0] for aln in best_alignments]
    number_of_references = len(references)
    taxid = ([TAXID_UNASSIGNED], False, None)
    taxpath = [taxid_2_path[reference] for reference in references]
    lca_taxpath = []
    for tax_rank in zip(*taxpath):
        if len(set(tax_rank)) == 1:
            lca_taxpath.append(tax_rank[0])
        else:
            break
    taxid = lca_taxpath
    return taxid


def read_id_to_taxpath(taxmap_file):
    file = pathlib.Path(taxmap_file)
    assert file.is_file()
    otu_2_path = {}
    for line in file.read_text().splitlines():
        splits = line.split('\t', 1)
        report_ranks = [rank for cnt, rank in enumerate(splits[1].split(';')) if cnt in REPORT_TAXONOMIC_INDEXES]
        otu_2_path[splits[0]] = report_ranks
    otu_2_path[TAXID_UNALIGNED] = ['root__Root', 'domain__Unaligned']
    rank_2_taxpath = collections.defaultdict(set)
    for _, taxpath in otu_2_path.items():
        for index, taxrank in enumerate(REPORT_TAXONOMIC_RANKS[::-1], 1):
            rank_2_taxpath[taxrank].add(';'.join(taxpath[:index]))

    return otu_2_path, rank_2_taxpath


def check_call(command: str):
    """
    Simple wrapper to execute check_call and catch exceptions
    :param command:
    :return:
    """

    returncode = -1
    try:
        returncode = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error('Command {} failed with message:\t{}'.format(e.cmd, e.stderr))
        shutdown(1)







def mitag_find(input_seq_file: pathlib.Path, output_folder: pathlib.Path, threads: int = 1):


    tmp_files = [] # to be deleted at the end of the program


    logging.info(f'Extracting FastA and revcomp FastA from {input_seq_file}')
    samplename = input_seq_file.name
    # Create forward fasta file in workdir
    fasta_forward = output_folder.joinpath(f'{samplename}_fw.fasta')
    fasta_reverse = output_folder.joinpath(f'{samplename}_rev.fasta')

    tmp_files.append(fasta_forward)
    tmp_files.append(fasta_reverse)

    with open(fasta_forward, 'w') as fw_handle:
        with open(fasta_reverse, 'w') as rev_handle:
            for number_of_sequences, fasta in enumerate(stream_fa(input_seq_file), 1):
                if number_of_sequences % 100000 == 0:
                    logging.info(f'Processed reads:\t{number_of_sequences}')
                header = fasta.header.split()[0]#re.sub('\s+', '_', fasta.header)
                fw_handle.write(f'>{header}\n{fasta.sequence}\n')
                rev_handle.write(f'>{header}\n{revcomp(fasta.sequence)}\n')

    logging.info(f'Finished extracting. Found {number_of_sequences} sequences.')



    logging.info(f'Start detecting rRNA sequences in FastA files')

    read2assignments = collections.defaultdict(list)

    for molecule in ['ssu', 'lsu']:
        logging.info(f'Start detecting rRNA sequences for molecule={molecule}')
        for (tmp_fasta_file, orientation) in [(fasta_forward, 'for'), (fasta_reverse, 'rev')]:
            folder = str(pathlib.Path(__file__).parent.parent.joinpath('data'))
            hmmer_pattern_file = pathlib.Path(folder + f'/{molecule}.hmm').absolute()
            hmmer_out_file = pathlib.Path(str(tmp_fasta_file) + f'_{molecule}.hmmer')
            dom_out_file = pathlib.Path(str(tmp_fasta_file) + f'_{molecule}.dom')

            tmp_files.append(hmmer_out_file)
            tmp_files.append(dom_out_file)

            hmmer_command = f'hmmsearch --cpu {threads} -o {hmmer_out_file} --domtblout {dom_out_file} -E 0.01 {hmmer_pattern_file} {tmp_fasta_file}'

            logging.info(f'Executing:\t{hmmer_command}')
            check_call(hmmer_command)
            logging.info(f'Finished hmmsearch')

            with open(dom_out_file) as handle:
                for line in handle:
                    if not line.startswith('#'):
                        [read, acc, tlen, qname, qaccr, qlen, seq_evalue, seq_score, seq_bias, seq_num, seq_of,
                         dom_c_evalue, dom_i_evalue, dom_score, dom_bias, hmm_start, hmm_end, dom_start, dom_end,
                         env_start, env_end, acc, description] = line.strip().split()

                        read2assignments[read.strip()].append([qname, dom_start, dom_end, read.strip(), float(dom_i_evalue), orientation, hmm_start, hmm_end])

        logging.info(f'Finished detecting rRNA sequences for molecule={molecule}')

    logging.info(f'Found {len(read2assignments)} potential rRNA sequences.')
    logging.info(f'Finished detecting rRNA sequences from FastA files.')

    logging.info(f'Finding best molecule for each read')
    read_2_bestassignment = {}
    for read, assignments in read2assignments.items():
        minevalue = min(map(lambda assignment: assignment[4], assignments))
        bestassignment = list(filter(lambda assignment: assignment[4] <= minevalue, assignments))[0]
        read_2_bestassignment[read] = bestassignment
    logging.info(f'Finished finding best molecule for each read')

    logging.info(f'Start extracting reads/writing output')
    writers = {}
    stats = collections.Counter()
    fasta_iterator = stream_fa(str(fasta_forward))
    for number_of_sequences, fasta in enumerate(fasta_iterator, 1):
        if number_of_sequences % 100000 == 0:
            logging.info(f'Processed reads:\t{number_of_sequences}')
        header = re.sub('\s+', '_', fasta.header)
        best_assignment = read_2_bestassignment.get(header, None)

        if best_assignment:
            stats[best_assignment[0]] += 1
            writer = writers.get(best_assignment[0], None)
            if not writer:
                o_file = output_folder.joinpath(f'{samplename}_{best_assignment[0]}.fasta')
                writer = open(o_file, 'w')
                writers[best_assignment[0]] = writer
            writer.write(f'>{header}\n{fasta.sequence}\n')

    for writer in writers.values():
        writer.close()
    logging.info(f'Finished extracting reads/writing output')
    for stat, count in stats.most_common():
        logging.info(f'{stat}\t{count}')
    for tmpfile in tmp_files:
        tmpfile.unlink()


















def execute_mitag_find(args: list):
    logging.info('Executing command find')
    parser = argparse.ArgumentParser(
        description='Takes a sequence file and extracts sequences that likely come from SSU/LSU genes.')

    parser.add_argument('-i', action='store', dest='input', help='Input fq/fa file. Can be gzipped.', required=True)
    parser.add_argument('-o', action='store', dest='output',help='Output folder.', required=True)
    parser.add_argument('-t', action='store', dest='threads', help='Number of threads for hmmsearch.',
                        default=1, type=int, choices=[1,2,4,8])


    results = parser.parse_args(args)



    if len(args) == 0:
        parser.print_help()
        shutdown(0)

    input_seq_file = pathlib.Path(results.input)
    output_folder = pathlib.Path(results.output)
    threads = results.threads

    if not input_seq_file.is_file():
        logging.error(f'Input {input_seq_file} file does not exist.')
        shutdown(1)
    if output_folder.is_file():
        logging.error('Output folder exist and is file. Select folder')
        shutdown(1)
    output_folder.mkdir(parents=True, exist_ok=True)




    mitag_find(input_seq_file, output_folder, threads=threads)

def test_and_merge_fasta_files(input_seqfiles_r1, input_seqfiles_r2, input_seqfiles_s):
    logging.info('Reading and testing FastA files')
    potentially_paired_sequences = collections.defaultdict(list)
    unpaired_sequences = collections.defaultdict(list)
    for read in itertools.chain([read for seqfile in input_seqfiles_r1 for read in stream_fa(seqfile)], [read for seqfile in input_seqfiles_r2 for read in stream_fa(seqfile)]):
        insert_name = read.header.split()[0]
        potentially_paired_sequences[insert_name].append(read)
    for read in [read for seqfile in input_seqfiles_s for read in stream_fa(seqfile)]:
        unpaired_sequences[read.header] = read

    actually_paired = list(filter(lambda insert: len(insert[1]) == 2, potentially_paired_sequences.items()))
    actually_single = list(filter(lambda insert: len(insert[1]) == 1, potentially_paired_sequences.items()))
    crazy_paired = list(filter(lambda insert: len(insert[1]) > 2, potentially_paired_sequences.items()))
    if len(crazy_paired) != 0:
        logging.error(f'Some inserts inserts have more then 2 reads assigned. This must be an error. Quitting...')
        for cp in crazy_paired:
            logging.error(cp)
        shutdown(1)
    logging.info('Paired files stats:')
    logging.info(f'\tFound {len(potentially_paired_sequences)} inserts in paired files of which {len(actually_paired)} are actually paired and {len(actually_single)} are singletons')
    logging.info('Singleton/Merged files stats:')
    logging.info(f'\tFound {len(unpaired_sequences)} inserts in singleton files.')

    for actually_single_read in actually_single:
        read = actually_single_read[1][0]
        unpaired_sequences[read.header] = read
    paired_sequences = {}
    for paired_read in actually_paired:
        if paired_sequences.get(paired_read[0], None):
            logging.error(f'Read {paired_read[0]} exists in paired and in singleton files.')
            shutdown(1)
        paired_sequences[paired_read[0]] =  paired_read[1]
    logging.info('Total stats:')
    logging.info(f'\tPaired inserts:\t{len(paired_sequences)}')
    logging.info(f'\tSingleton inserts:\t{len(unpaired_sequences)}')
    logging.info('Finished reading and testing FastA files')
    return paired_sequences, unpaired_sequences











def mitag_annotate_vsearch(input_seqfiles_r1, input_seqfiles_r2, input_seqfiles_s, database, taxmap, out_file_pattern, threads=4, enable_binning=True, maxaccepts=1000, maxrejects=1000):
    tmp_files = []
    # 1. Prepare
    # 1.1. Create 1 R1, 1 R2 and 1 S file with the output pattern
    # 1.2. Check if you find unpaired sequences in the R1 files
    #input_seqfile_r1, input_seqfile_r2, input_seqfile_s =

    paired_sequences, unpaired_sequences = test_and_merge_fasta_files(input_seqfiles_r1, input_seqfiles_r2, input_seqfiles_s)
    pe_fasta_file = pathlib.Path(out_file_pattern).with_suffix(f'.pe.fasta')
    se_fasta_file = pathlib.Path(out_file_pattern).with_suffix(f'.se.fasta')
    tmp_files.append(pe_fasta_file)
    tmp_files.append(se_fasta_file)
    logging.info(f'Writing paired sequences to FastA:\t{pe_fasta_file}')
    with open(pe_fasta_file, 'w') as handle:
        for insert_name, reads in paired_sequences.items():
            handle.write(f'>{insert_name}.1\n{reads[0].sequence}\n')
            handle.write(f'>{insert_name}.2\n{reads[1].sequence}\n')
    logging.info(f'Writing finished')
    logging.info(f'Writing singleton sequences to FastA:\t{se_fasta_file}')
    with open(se_fasta_file, 'w') as handle:
        for insert_name, read in unpaired_sequences.items():
            handle.write(f'>{read.header}\n{read.sequence}\n')
    logging.info(f'Writing finished')
    # 2. Map reads against database

    pe_m8_file = pathlib.Path(out_file_pattern).with_suffix(f'.pe.m8')
    se_m8_file = pathlib.Path(out_file_pattern).with_suffix(f'.se.m8')



    aln_command_template = 'vsearch --usearch_global {} --db {} --id 0.97 --maxaccepts {} --maxrejects {} --strand both --userout {} --userfields query+target+qcov+id+qilo+qihi+tilo+tihi+alnlen+ids+mism+gaps+qstrand --query_cov {} --threads {} --output_no_hits'
    pe_aln_command = aln_command_template.format(pe_fasta_file, database, maxaccepts, maxrejects, pe_m8_file, MTAGS_CLASSIFY_QCOV, threads)
    se_aln_command = aln_command_template.format(se_fasta_file, database, maxaccepts, maxrejects, se_m8_file, MTAGS_CLASSIFY_QCOV, threads)
    if True:
        logging.info(f'Aligning paired sequence file:\t{pe_aln_command}')
        if len(paired_sequences) == 0:
            pe_m8_file.touch()
        else:
            check_call(pe_aln_command)
        logging.info(f'Finished alignment')
        logging.info(f'Aligning singleton sequence file:\t{se_aln_command}')
        if len(unpaired_sequences) == 0:
            se_m8_file.touch()
        else:
            check_call(se_aln_command)
        logging.info(f'Finished alignment')
    # 3. Read mapping files
    logging.info(f'Parsing alignment files')
    aligned_inserts = read_and_pair_vsearch_tab_files(pe_m8_file, se_m8_file)
    logging.info(f'Finished parsing alignment files')
    # 4. Filter
    logging.info(f'Filtering alignments by score and length')
    insert_2_best_alignments = {}
    for insert_name, aligned_insert in aligned_inserts.items():
        best_alignments = aligned_insert.extract_best_alignments()
        insert_2_best_alignments[insert_name] = best_alignments
    logging.info(f'Finished filtering')
    logging.info(f'Performing LCA calculations')
    otu_2_taxpath, rank_2_taxpath = read_id_to_taxpath(taxmap)

    insert_2_lca = {}
    for insert_name, best_alignments in insert_2_best_alignments.items():
        lca = lca_vsearch(best_alignments, otu_2_taxpath)
        insert_2_lca[insert_name] = lca
    logging.info(f'Finished LCA calculations')

    logging.info(f'Start writing output files')





    if enable_binning:
        of = pathlib.Path(out_file_pattern).with_suffix(f'.bins')
        logging.info(f'Start writing:\t{of}')
        with open(of, 'w') as binning_file:
            for insert_name, lca in insert_2_lca.items():
                binning_file.write('{}\t{}\t{}\n'.format(insert_name, lca[-1].split(RANK_2_NAME_SEPARATOR)[0], ';'.join(lca)))



    all_counts = {}

    lcas = list(insert_2_lca.values())
    total_insert_count = len(insert_2_lca.values())
    updated_lcas = []
    for taxrank in REPORT_TAXONOMIC_RANKS:

        total_inserts = len(insert_2_lca)
        total_assigned_inserts = 0
        taxon_2_count = collections.Counter()
        for lca in lcas:
            if lca[-1].split(RANK_2_NAME_SEPARATOR)[0] == taxrank:
                total_assigned_inserts += 1
                taxon_2_count[';'.join(lca)] += 1
                updated_lcas.append(lca[:-1])
            else:
                updated_lcas.append(lca)
        taxpath_this_domain = []
        for taxpath in rank_2_taxpath[taxrank]:
            taxpath_this_domain.append(taxpath)
        cnter = collections.Counter()
        for taxpath in sorted(taxpath_this_domain):
            cnter[taxpath] = taxon_2_count.get(taxpath, 0)
        all_counts[taxrank] = cnter
        lcas = updated_lcas
        updated_lcas = []

    for taxrank, counts in all_counts.items():
        if taxrank in REPORT_TAXONOMIC_RANKS:
            of = pathlib.Path(out_file_pattern).with_suffix(f'.{taxrank}.tsv')
            logging.info(f'Start writing:\t{of}')
            with open(of, 'w') as handle:
                total_assigned_inserts = 0
                for taxpath in sorted(counts.keys()):
                    total_assigned_inserts += counts[taxpath]
                    if counts[taxpath] == 0:
                        continue
                    else:
                        handle.write(f'{taxpath}\t{counts[taxpath]}\n')
                if taxrank not in ['domain', 'root']:
                    total_unaligned_inserts = all_counts['domain']['root__Root;domain__Unaligned']
                    handle.write(f'Unaligned\t{total_unaligned_inserts}\n')
                    handle.write(f'Unassigned\t{total_insert_count - total_assigned_inserts - total_unaligned_inserts}\n')
    logging.info(f'Finished writing output files')
















def execute_mitag_annotate(args):
    logging.info('Executing command annotate')
    parser = argparse.ArgumentParser(
        description='Assigns LCA based taxonomy to rRNA ssu reads.')

    parser.add_argument('-i1', action='store', nargs='+', dest='i1', help='Input R1 fasta file(s). Files need to have the same order as in i2. Sequences in paired files will be matched by name. Unmatched sequences will be counted as singletons. Insert names have to be identical after removal of (.1 .2) or (/R1 /R2) or (/1 /2).')
    parser.add_argument('-i2', action='store', nargs='+', dest='i2', help='Input R2 fasta file(s). Files need to have the same order as in i1. Sequences in paired files will be matched by name. Unmatched sequences will be counted as singletons. Insert names have to be identical after removal of (.1 .2) or (/R1 /R2) or (/1 /2).')
    parser.add_argument('-is', action='store', nargs='+', dest='im', help='Input singletons/merged fasta file(s).')

    parser.add_argument('-o', action='store', dest='o',help='Output pattern. This pattern will be prepended to all output files',required=True)
    parser.add_argument('-t', action='store', dest='t', help='Number of threads for vsearch alignment', default=4, type=int)
    #parser.add_argument('-d', action='store', dest='db', help='Path to the vsearch database (udb or fasta)', required=True)
    parser.add_argument('-ma', action='store', dest='ma', help='Maxaccepts, vsearch parameter. Larger numbers increase sensitivity and runtime.', default=1000, type=int)
    parser.add_argument('-mr', action='store', dest='mr', help='Maxrejects, vsearch parameter. Larger numbers increase sensitivity and runtime.', default=1000, type=int)
    #parser.add_argument('-m', action='store', dest='tm', help='Taxmap file that maps tax ids to taxonomic path.', required=True)
    parser.add_argument('-b', dest='bin', help='Enable binning - write taxonomic assignment for each read/insert', action='store_true', default=False)
    results = parser.parse_args(args)
    if len(args) == 0:
        parser.print_help()
        shutdown(0)
    #######
    #Input#
    #######
    input_seqfiles_r1 = results.i1
    input_seqfiles_r2 = results.i2
    input_seqfiles_s = results.im

    #database = results.db
    #taxmap = results.tm

    ########
    #Output#
    ########
    out_file_pattern = results.o


    ########
    #Params#
    ########
    threads = results.t
    enable_binning = results.bin
    maxaccepts = results.ma
    maxrejects = results.mr



    ########
    #Checks#
    ########

    assert threads > 0
    assert maxaccepts > 0
    assert maxrejects > 0

    if not input_seqfiles_r1:
        input_seqfiles_r1 = []
    if not input_seqfiles_r2:
        input_seqfiles_r2 = []
    if not input_seqfiles_s:
        input_seqfiles_s = []


    for input_seqfile in itertools.chain(input_seqfiles_r1, input_seqfiles_r2, input_seqfiles_s):
        if pathlib.Path(input_seqfile).exists() and pathlib.Path(input_seqfile).is_file():
            continue
        else:
            logging.error(f'Input file {input_seqfile} does not exist or is not file.')
            shutdown(1)
    if len(input_seqfiles_r1) != len(input_seqfiles_r2):
        logging.error(f'Different number of R1 and R2 files provided. #R1={len(input_seqfiles_r1)}, #R2={len(input_seqfiles_r2)}')
        shutdown(1)

    if not db_marker_file.exists():
        logging.error(f'The database is not downloaded. Please run `mtags download` first')
        shutdown(1)



    if pathlib.Path(out_file_pattern).is_dir():
        logging.error(f'The output pattern cannot be a directory.')
        shutdown(1)
    pathlib.Path(out_file_pattern).parent.mkdir(parents=True, exist_ok=True)
    mitag_annotate_vsearch(input_seqfiles_r1, input_seqfiles_r2, input_seqfiles_s, database, taxmap, out_file_pattern, threads=threads, enable_binning=enable_binning, maxaccepts=maxaccepts, maxrejects=maxrejects)
    logging.info('Finished annotate command')




def execute_mtags_download():
    database_url = 'https://sunagawalab.ethz.ch/share/MTAGS_DB/silva/s138/SILVA-138_NR-97_complink_cons.vsearch.udb.gz'
    taxmap_url = 'https://sunagawalab.ethz.ch/share/MTAGS_DB/silva/s138/SILVA-138_NR-97_complink_cons.taxmap.gz'
    user = 'mtags'
    password = 'mtags'


    if db_marker_file.exists():
        logging.info('Database is already downloaded. Quitting')
    else:
        logging.info('Start downloading the mTAGs database. ~600MB')
        _execute_download(database_url, database, user, password)
        _execute_download(taxmap_url, taxmap, user, password)

        db_marker_file.touch(exist_ok=True)
        logging.info('Finished downloading the mTAGs database.')


def _execute_download(url, destfile, user, password):
    """Method to download and extract the mTAGs database

    """

    p = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    p.add_password(None, url, user, password);

    auth_handler = urllib.request.HTTPBasicAuthHandler(p)
    opener = urllib.request.build_opener(auth_handler)

    urllib.request.install_opener(opener)
    tfile = str(destfile) + '.gz'
    try:

        req = urllib.request.Request(url, headers={'Content-Type': 'application/json'})
        result = opener.open(req)
        with open(tfile, 'wb') as handle:
            handle.write(result.read())
    except IOError as e:
        logging.error(e)
        shutdown(1)

    try:
        returncode = subprocess.check_call(f'gunzip -c {tfile} > {destfile}',shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(e)
        shutdown(1)
    pathlib.Path(tfile).unlink()



def main():
    parser = argparse.ArgumentParser(description='A toolkit to extract and annotate rRNA reads from metagenomic data',
                                     usage='''mtags <command> [<args>]
    Command options
        download    Download the mTAGs database
        find        Searches for potential rRNA reads in metagenomic samples. (First step)
        annotate    Annotate and quantify rRNA reads. (Second Step)
    ''')

    parser.add_argument('command', help='Subcommand to run: find|annotate|download', choices=['find', 'annotate', 'download'])

    args = parser.parse_args(sys.argv[1:2])
    startup()
    '''
    Paired end read files are identified by name. The names should be identical between the 2 files
    '''

    command = args.command
    if command == 'find':
        execute_mitag_find(sys.argv[2:])
        shutdown(0)
    elif command == 'annotate':
        execute_mitag_annotate(sys.argv[2:])
        shutdown(0)
    elif command == 'download':
        execute_mtags_download()
        shutdown(0)
    else:
        parser.print_help()
        shutdown(0)
    shutdown(0)


if __name__ == '__main__':
    main()

