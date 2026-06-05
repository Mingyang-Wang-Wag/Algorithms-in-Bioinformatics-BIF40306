"""usage: python <script.py> <TAIR10.fasta> <query.fasta> """

# !/usr/bin/env python3
"""
Author:
Student number:
Implementation of the SSAHA algorithm

Hints:
- write a function to generate a hash table of k-mers from a set of sequences
- write a function to return a sorted list of hits in the hash table 
for a query sequence
- write a function to find the best hit
- write a fasta parser to read in the Arabidopsis data

"""
# import statements
from sys import argv

def fasta_reader(filename):
    seq_ls = []
    non_atcg_num = 0
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
                seq += clean_line
            else:
                raise ValueError('Sequence should only contain letters')

        if seq != '': #the last seq
            seq_ls.append(seq)

    return seq_ls, non_atcg_num

def kmer_maker(seq_ls, k):
    """

    Return: genome_kmer_dict, dictionary, {seq1: {kmer1:[positions in seq], }, }
    """
    genome_kmer_dict = {}
    for seq in seq_ls:
        genome_kmer_dict[seq] = {}
        for i in range(len(seq)):
            if i <= len(seq) - k: # if k = 3, the last two cannot be a kmer
                kmer = seq[i:i+k]
                if kmer not in genome_kmer_dict[seq]:
                    genome_kmer_dict[seq][kmer] = []
                    genome_kmer_dict[seq][kmer].append(i) #append position
                else:
                    genome_kmer_dict[seq][kmer].append(i)

    return genome_kmer_dict

#store as (2, 7, 11),
#2 is the No. of sequence, 7 is the shift, 11 is the position at reference

def match_table_maker(query_kmer_dict, genome_kmer_dict):

    kmer_posi_dict = {}
    for index, seq in enumerate(genome_kmer_dict, start=1):
        for kmer in query_kmer_dict: #only search the query kmer in the ref
            if kmer in genome_kmer_dict[seq]:
                if kmer not in kmer_posi_dict:
                    kmer_posi_dict[kmer] = []
                    #no kmer then create a empty list, append posi
                    #else, means there is a list already, append posi

                for ref_posi in genome_kmer_dict[seq][kmer]:
                    for query_posi in query_kmer_dict[kmer]:
                        shift = ref_posi - query_posi
                        kmer_posi_dict[kmer].append((index, shift, ref_posi))


    for ls in kmer_posi_dict.values():
        ls.sort()

    return kmer_posi_dict

def cluster_hits(kmer_posi_dict):
    hits_dict = {}
    for kmer, hits_ls in kmer_posi_dict.items():
        for hit in hits_ls:
            seq_num = hit[0]
            shift = hit[1]
            ref_posi = hit[2]
            if seq_num not in hits_dict:
                hits_dict[seq_num] = {}

            if shift not in hits_dict[seq_num]:
                hits_dict[seq_num][shift] = []
            hits_dict[seq_num][shift].append((ref_posi, kmer))

    return hits_dict

def best_hit_finder(hits_dict, query, k):
    best_hit_dict = {}
    for seq_num, shifts_dict in hits_dict.items():
        for shift, kmer_posi_ls in shifts_dict.items():
            if len(kmer_posi_ls) == len(query) - k + 1:
                if seq_num not in best_hit_dict:
                    best_hit_dict[seq_num] = []

                best_hit_dict[seq_num].append(kmer_posi_ls)


    for best_hit_ls in best_hit_dict.values():
        for best_hit in best_hit_ls:
            best_hit.sort()

    return best_hit_dict

def matching_printer(best_hit_dict, query, seq_ls):
    #relocate the match part into the original sequence
    seq_dict = {}
    for index, seq in enumerate(seq_ls, start=1):
        if len(seq) > 100:
            print(f'matching_printer do not fit for long sequence')
            exit()
        seq_dict[index] = seq

    for seq_num, best_hit_ls in best_hit_dict.items():
        for best_hit in best_hit_ls:
            space_num = best_hit[0][0]
            space = ' ' * space_num
            print(f'{space}{query}')
            print(f'\n')
            seq = seq_dict[seq_num]
            print(f'{seq}')

#another solution, print the matching part before and after + 5
#also print the relative position


def main():
    """Perform these commands when called from the command line
    """
    if len(argv) != 3: #change the number
        print("Usage: python <script.py> <TAIR10.fasta> <query.fasta>")
        exit()
    ref_filename = argv[1]
    seq_ls, non_atcg_num = fasta_reader(ref_filename)

    query_filename = argv[2]
    query_ls, non_atcg_num_query = fasta_reader(query_filename)
    k = 13

    '''######for devlopment
    query_ls = ['TGCAACAT']
    
    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seq_ls = [s1, s2, s3]
    '''

    # Running
    genome_kmer_dict = kmer_maker(seq_ls, k)
    query_kmer_dict = kmer_maker(query_ls, k)

    for query, kmer_dict in query_kmer_dict.items():
    
        kmer_posi_dict = match_table_maker(kmer_dict, genome_kmer_dict)
        hits_dict = cluster_hits(kmer_posi_dict)
        best_hit_dict = best_hit_finder(hits_dict, query, k)
        for seq_num, best_hit in best_hit_dict.items():
            print(f'sequence{seq_num}, have {len(best_hit)} best hits')
        #if len > 100
        #else:
        #matching_printer(best_hit_dict, query, seq_ls)

#with query.fasta and TAIR10.fasta, the memory explode

#even with query.fasta and seq1,2,3 the kmer_posi_dict is a lot.


#Time records: 7:30 running time, checking letter one by one, one short query

#Time records: checking letter one by one, query.fasta


if __name__ == "__main__":
    # the code in your main function should produce the results necessary
    # to answer the questions. In other words, if we run your code, we should
    # see the data that you used to answer the questions
    main()

