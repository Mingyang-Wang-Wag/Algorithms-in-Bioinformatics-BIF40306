#!/usr/bin/env python3
"""
Name: Mingyang Wang
Student Number: 1067192
Description: Parsing the BED file and extract the start and stop coordinates
             and strand for each CDS feature, calculate the total number and
             length of features on positive and negative strand. Print this
             result in the terminal.
             Run bedtools getfasta program to extract the coding sequences
             of whole genome sequences based on BED file, create an output file
             with a name based on the input genome_filename.
Usage: python3 <script_name.py> <genome_filename> <bed_filename>
    where script_name.py: str, is the name of this script (0705_2025exam.py)
        genome_filename: str, the file name of a genome file containing
        DNA sequence. e.g.ecoli_genome.fa
        bed_filename: str, the filename of a bed file, containing the start,
        end, strand etc. info for each CDS in the BED format.
        e.g. ecoli_cds.bed
"""

from pathlib import Path
from sys import argv
import subprocess

def bed_reader(bed_filename):
    """Parse a BED file into a list

        bed_filename: str, the name of the bed file, the bed file contains
        the start, end, the strand etc. info for each CDS in the BED format.
        return(s): a list, this main list contains sub-lists,
                for each sub-list, it is [start, end, strand],
                where start: int, the start coordinate of a CDS feature
                      end: int, the end coordinate of a CDS feature
                      strand: str, on which strand (+ or -) of a CDS feature
    """
    bed_ls = []
    with open(bed_filename) as f:

        for line in f:
            line = line.strip()

            if line == '':
                continue

            parts = line.split()
            if len(parts) == 9:
                start = int(parts[3])
                end = int(parts[4])
                strand = str(parts[6])
                ls = [start, end, strand]
                bed_ls.append(ls)

    return bed_ls

def strand_seperator(bed_ls):
    """Parsing the start, end coordinates based on the strand into a dictionary

        bed_ls: a list, the main list contains sub-lists [start, end, strand]
        return(s): strand_dict, a dictionary, key: str, '+' and '-'
                 value: list of list, the main list contains sub-lists
                 each sub-list is [start, end],
                 where: start: int, the start coordinate of a CDS feature
                      end: int, the end coordinate of a CDS feature
    """
    strand_dict = {}
    for fraction in bed_ls:
        strand = fraction[2]
        start = fraction[0]
        end = fraction[1]

        if strand not in strand_dict.keys():
            strand_dict[strand] = [[start, end]]
        else:
            strand_dict[strand].append([start, end])

    return strand_dict

def length_calculator(strand_dict):
    """Calculate the total number and length for positive and negative strand

        strand_dict: a dictionary, key: str, '+' and '-'; values, list of list,
                     the main list contains sub-list of [start, end].
        return(s): total_neg, total_pos, neg_length, pos_length
                where, total_neg: int, the total numbers of features on negative
                strand
                total_pos: int, the total numbers of features on positive
                strand
                neg_length:int, the total length of features on negative
                strand
                pos_length: int, the total length of features on positive
                strand
    """
    total_neg = len(strand_dict['-'])
    total_pos = len(strand_dict['+'])

    neg_length = 0
    for fraction in strand_dict['-']:
        start = fraction[0]
        end = fraction[1]
        fra_len = end - start
        neg_length += fra_len

    pos_length = 0
    for fraction in strand_dict['+']:
        start = fraction[0]
        end = fraction[1]
        fra_len = end - start
        pos_length += fra_len

    return total_neg, total_pos, neg_length, pos_length

def output(total_neg, total_pos, neg_length, pos_length):
    """Print out the four arguments in the required format

        total_neg, total_pos: all are integers, the total numbers of features on
        negative and positive strand.
        neg_length, pos_length: all are integers, the total length of features
        on negative and positive strand.
        return(s): None
    """
    print(f'{total_pos} features on the + strand (total length: {pos_length})')
    print(f'{total_neg} features on the - strand (total length: {neg_length})')

def run_bedtools(genome_filename, bed_filename):
    """Run the bedtools getfasta program in command line, create an output file

        genome_filename: str, the name of genome file, the file contains the
                        sequence of the whole genome
        bed_filename: str, the name of the bed file, the bed file contains
        the start, end, the strand etc. info for each CDS in the BED format.
        Create an output file with a name: <species>_cds.fa, where species is
        a string taken from the genome filename
        return(s): None
    """
    parts = genome_filename.split('_') #ecoli_genome.fa
    bed_output_name = parts[0] + '_cds.fa'

    bed_output_path = Path(bed_output_name)

    if not bed_output_path.exists():
        subprocess.run(['bedtools',
                        'getfasta',
                        '-fi',
                        f'{genome_filename}',
                        '-fo',
                        f'{bed_output_name}',
                        '-bed',
                        f'{bed_filename}'
                        ], check=True)  # always with check = True

def main():
    """Main function of this script"""
    if len(argv) != 3:
        print("Usage: python3 <script.py> <genome_filename> <bed_filename>")
        exit()

    genome_filename = argv[1]
    bed_filename = argv[2]

    genome_real_path = Path('/home/wang394/resit/data') / genome_filename
    genome_path = Path(genome_filename)

    bed_real_path = Path('/home/wang394/resit/data') / bed_filename
    bed_path = Path(bed_filename)

    if not genome_path.exists():
        genome_path.symlink_to(genome_real_path)

    if not bed_path.exists():
        bed_path.symlink_to(bed_real_path)

    bed_ls = bed_reader(bed_filename)
    strand_dict = strand_seperator(bed_ls)
    total_neg, total_pos, neg_length, pos_length =length_calculator(strand_dict)
    output(total_neg, total_pos, neg_length, pos_length)
    run_bedtools(genome_filename,bed_filename)

if __name__ == "__main__":
    main()
