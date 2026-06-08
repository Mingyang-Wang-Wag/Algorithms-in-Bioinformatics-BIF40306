# !/usr/bin/env python3
"""
Author: Mingyang Wang
Student No.: 1067192

Description:
    Implementation of the SSAHA (Sequence Search and Alignment by Hashing
    Algorithm) algorithm for DNA sequence searching. The program builds a
    hash table of non-overlapping k-mers from reference sequences and searches
    overlapping k-mers from query sequences against this hash table. Matching
    k-mers are clustered according to sequence index and shift, and the
    cluster with the largest number of hits is reported as the best match.

Usage:
    python <script.py> <TAIR10.fasta> <query.fasta>

    where:
        <script.py>: the name of this script.
        <TAIR10.fasta>: a FASTA file containing one or more reference DNA
            sequences.
        <query.fasta>: a FASTA file containing one or more query DNA
            sequences.

Output:
    For a small example dataset, the program prints the hash table, number of
    hits, and best match. For the Arabidopsis reference genome and query
    sequences, it reports the total number of hits and the best hit
    for each query sequence.
"""

from sys import argv

def fasta_reader(filename):
    """Parse a FASTA file into a list of DNA sequences.

    Ambiguous DNA letters such as N, Y, R, W, S, K, M and D are converted
    to A. Header lines starting with ">" are skipped.

    Parameters:
        filename: A string, the name of the input FASTA file.

    Returns:
        seq_ls: A list of strings, where each string is one DNA sequence.
    """
    seq_ls = []
    table = str.maketrans({
        "N": "A",
        "Y": "A",
        "R": "A",
        "W": "A",
        "S": "A",
        "K": "A",
        "M": "A",
        "D": "A",
    })
    with open(filename) as f:
        seq = ''
        for line in f:
            line = line.strip()
            if line == '':
                continue

            if line.startswith('>'):
                if seq != '':
                    seq_ls.append(seq)
                seq = '' #skip the > line
                continue

            if line.isalpha(): #idea from Chatgpt to deal with ambigiou letter
                line = line.upper()
                line = line.translate(table)
                '''
                clean_line = ''
                for let in line:
                    if let in 'ATCG':
            #in the TAIR10 sequence there is Y, which is ambigous dna letter
            #i only stick with ATCG, then it cannot proceed.
            # or if let != 'A' and let != 'C' and let != 'G' and let != 'T':
                        clean_line += let
                    else: #convert all ambigious letter to A
                        non_atcg_num += 1 #keep tracking
                        clean_line += 'A'
                '''
                seq += line
            else:
                raise ValueError('Sequence should only contain letters')

        if seq != '': #the last seq
            seq_ls.append(seq)

    return seq_ls


def hashtable_maker(seq_ls, k):
    """Create a hash table of non-overlapping k-mers from subject sequences.

    Parameters:
        seq_ls: A list of DNA sequences as strings.
        k: An integer, the length of each k-mer.

    Returns:
        genome_kmer_dict: A dictionary where each key is a k-mer string and
            each value is a list of tuples in the form
            {kmer:[(seq_index, position)]}.
            kmer, str, small fragment of sequence, e.g. 'TG'
            seq_index, int, the No. of sequence store in the seq_ls
            position, int, the index of kmer on the sequence
            The seq_index and position are 1-based.
    """
    genome_kmer_dict = {}
    for seq_index, seq in enumerate(seq_ls, start=1):
        if len(seq) < k:
            raise ValueError('The length of sequence '
                             'should not be shorter than the kmer')
        for i in range(0, len(seq) - k + 1, k):
            # if k = 3, the last two cannot be a kmer
            kmer = seq[i:i+k]
            if kmer not in genome_kmer_dict:
                genome_kmer_dict[kmer] = []

            genome_kmer_dict[kmer].append((seq_index, i+1))


    return genome_kmer_dict

def query_hashtable_maker(query_ls, k):
    """Create a hash table of overlapping k-mers from query sequences.

    Parameters:
        query_ls: A list of query DNA sequences as strings.
        k: An integer, the length of each k-mer.

    Returns:
        query_kmer_dict: A nested dictionary in the form
            {query: {kmer: [query_positions]}}.
            Query positions are 0-based.
        query, str, is the input query used for searching
        kmer, str, small fragment of sequence
        query_position, int, 0-based position number at query
    """
    query_kmer_dict = {}
    for query in query_ls:
        if len(query) < k:
            raise ValueError('The length of sequence '
                             'should not be shorter than the kmer')
        query_kmer_dict[query] = {}
        for i in range(0, len(query) - k + 1): #overlapping kmer
            kmer = query[i:i+k]
            if kmer not in query_kmer_dict[query]:
                query_kmer_dict[query][kmer] = []

            query_kmer_dict[query][kmer].append(i)

    return query_kmer_dict


def match_table_maker(query_kmer_dict, genome_kmer_dict):
    """Find k-mer hits between query sequences and subject sequences.

    Parameters:
        query_kmer_dict: A nested dictionary from query_hashtable_maker.
        genome_kmer_dict: A dictionary from hashtable_maker.

    Returns:
        kmer_posi_dict: A nested dictionary in the form
            {query: {kmer: [(seq_index, shift, ref_posi)]}}.
        query, str, is the input DNA query used for searching
        kmer, str, small fragment of sequence
        seq_index, int, the No. of sequence store in the seq_ls,1-based position
        shift, int, is calculated as ref_posi - query_posi.
        ref_position, int, 1-based position number at query

    """
    kmer_posi_dict = {}
    for query, query_dict in query_kmer_dict.items():
        kmer_posi_dict[query] = {}
        for kmer, query_posi_ls in query_dict.items():
            if kmer in genome_kmer_dict:
                for seq_index, ref_posi in genome_kmer_dict[kmer]:
                    for query_posi in query_posi_ls:
                        shift = ref_posi - query_posi
                        if kmer not in kmer_posi_dict[query]:
                            kmer_posi_dict[query][kmer] = []

                        kmer_posi_dict[query][kmer].append((seq_index, shift,
                                                     ref_posi))

                kmer_posi_dict[query][kmer].sort()


    return kmer_posi_dict

def cluster_hits(kmer_posi_dict):
    """Group k-mer hits by query sequence, subject sequence, and shift.

    Parameters:
        kmer_posi_dict: A nested dictionary from match_table_maker.

    Returns:
        hits_dict: A nested dictionary in the form
            {query: {seq_index: {shift: [(ref_posi, kmer)]}}}.

        query, str, is the input DNA query used for searching
        seq_index, int, the No. of sequence store in the seq_ls,1-based position
        shift, int, is calculated as ref_posi - query_posi.
        ref_position, int, 1-based position number at query
        kmer, str, small fragment of sequence
    """
    hits_dict = {}

    for query, query_dict in kmer_posi_dict.items():
        hits_dict[query] = {}

        for kmer, hits_ls in query_dict.items():
            for seq_index, shift, ref_posi in hits_ls:

                if seq_index not in hits_dict[query]:
                    hits_dict[query][seq_index] = {}

                if shift not in hits_dict[query][seq_index]:
                    hits_dict[query][seq_index][shift] = []

                hits_dict[query][seq_index][shift].append((ref_posi, kmer))

    return hits_dict

def best_hit_finder(hits_dict):
    """Find the best hit cluster for each query sequence.

    Parameters:
        hits_dict: A nested dictionary in the form
            {query: {seq_index: {shift: [(ref_posi, kmer)]}}}.
        query, str, is the input DNA query used for searching
        seq_index, int, the No. of sequence store in the seq_ls,1-based position
        shift, int, is calculated as ref_posi - query_posi.
        ref_position, int, 1-based position number at query
        kmer, str, small fragment of sequence

    Returns:
        best_hit_dict: A dictionary in the form
            {query: {"seq_index": seq_index, "shift": shift, "hits": hits}}.
        where query, str, short DNA sequence for searching
        seq_index, int, the No. of sequence store in the seq_ls,1-based position
        shift, int, is calculated as ref_posi - query_posi.
        hits is list of tuple [(sequence position, kmer)]
        sequence position, int, the position of kmer on the sequence
        kmer, str, small fragment of DNA sequence
    """
    best_hit_dict = {}

    for query, seq_dict in hits_dict.items():
        best_seq = None
        best_shift = None
        best_hits = []

        for seq_index, shift_dict in seq_dict.items():
            for shift, hits in shift_dict.items():

                if len(hits) > len(best_hits):
                #len(kmer_posi_ls) == len(query) - k + 1, too strict
                #but
                    best_seq = seq_index
                    best_shift = shift
                    best_hits = hits

        best_hits.sort()



        best_hit_dict[query] = {
            "seq_index": best_seq,
            "shift": best_shift,
            "hits": best_hits
        }

    return best_hit_dict


def matching_printer(best_hit_dict, seq_ls):
    """Print a simple visual alignment of each query against its best hit.

    This function is only intended for short subject sequences. For long
    sequences, it raises a ValueError.

    Parameters:
        best_hit_dict: A dictionary in the form
            {query: {"seq_index": seq_index, "shift": shift, "hits": hits}}.
        where query, str, short DNA sequence for searching
        seq_index, int, the No. of sequence store in the seq_ls,1-based position
        shift, int, is calculated as ref_posi - query_posi.
        hits is list of tuple [(sequence position, kmer)]
        sequence position, int, the position of kmer on the sequence
        kmer, str, small fragment of DNA sequence

        seq_ls: A list of DNA sequences as strings.

    Returns:
        None.
    """
    #relocate the match part into the original sequence
    seq_dict = {}
    for index, seq in enumerate(seq_ls, start=1):
        if len(seq) > 100:
            raise ValueError(f'matching_printer do not fit for long sequence')


    for query, hit_dict in best_hit_dict.items():
        seq_index = hit_dict['seq_index']
        shift = hit_dict['shift']
        space = ' ' * (shift - 1)
        print(f'{space}{query}')
        print(f'\n')
        print(f'{seq_ls[seq_index - 1]}')

#another solution, print the matching part before and after + 5
#also print the relative position


def main():
    """Perform these commands when called from the command line
    """
    if len(argv) != 3: #change the number
        print("Usage: python <script.py> <TAIR10.fasta> <query.fasta>")
        exit()


    k_short = 2
    query_ls_short = ['TGCAACAT']
    
    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seq_ls_short = [s1, s2, s3]


    # Running for short sequence and query
    genome_kmer_dict_short = hashtable_maker(seq_ls_short, k_short)
    print(f'Hash table of short sequences: {genome_kmer_dict_short}\n')
    query_kmer_dict_short = query_hashtable_maker(query_ls_short, k_short)
    kmer_posi_dict_short= match_table_maker(query_kmer_dict_short,
                                      genome_kmer_dict_short)
    total_hits = 0
    for query, kmer_hits in kmer_posi_dict_short.items():
        for key, hits in kmer_hits.items():

            total_hits += len(hits)
        print(f'There are {total_hits} hits in total')
    #my data structure is not a list, the first key of the dictionary is 'TG',
    #the last key is 'AT'.
    print(f'The first hit: {kmer_posi_dict_short['TGCAACAT']['TG'][0]}')
    print(f'The last hit: {kmer_posi_dict_short['TGCAACAT']['AT'][-1]}\n')

    hits_dict_short = cluster_hits(kmer_posi_dict_short)
    best_hit_dict_short = best_hit_finder(hits_dict_short)
    print(f'{best_hit_dict_short}\n')
    matching_printer(best_hit_dict_short, seq_ls_short)
    print(f'\n')


    # Running for TAIR10 sequence and list of queries
    k = 13
    ref_filename = argv[1]
    seq_ls= fasta_reader(ref_filename)

    query_filename = argv[2]
    query_ls= fasta_reader(query_filename)

    genome_kmer_dict = hashtable_maker(seq_ls, k)
    print(len(genome_kmer_dict))

    query_kmer_dict = query_hashtable_maker(query_ls, k)
    kmer_posi_dict= match_table_maker(query_kmer_dict, genome_kmer_dict)
    total_hits = 0
    for query, kmer_hits in kmer_posi_dict.items():
        for key, hits in kmer_hits.items():
            total_hits += len(hits)
    print(f'With k = {k}, the number of hits is {total_hits}')

    hits_dict = cluster_hits(kmer_posi_dict)
    best_hit_dict = best_hit_finder(hits_dict)
    print(best_hit_dict)
    #for query, hits_dict in best_hit_dict.items():
        #print(len(hits_dict['hits']))


if __name__ == "__main__":
    main()

