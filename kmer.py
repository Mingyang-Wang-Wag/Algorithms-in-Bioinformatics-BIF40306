#!/usr/bin/env python3
"""
Name: Mingyang Wang
Student Number: 1067192

Description:
    Read a DNA sequence from a FASTA input file and calculate its 4-mer
    composition.It prints the frequencies of all possible 4-mers in
    lexicographic order on one line.

Usage:
    python3 <script_name.py> <input_file>

    where:
        script_name.py:
            str, Name of this Python script.

        input_file:
            str, name of the input FASTA file containing a DNA sequence.
            The file may contain a FASTA header line starting with '>'.
"""


from sys import argv

def kmer_library(): #what a nice method!!!!
    """Generate all possible DNA 4-mers in lexicographic order.

        Input: None

        Returns:
            A list containing all possible 4-mers in lexicographic order.
    """


    lib = []
    alphabet = 'ACGT'

    for first in alphabet:
        for second in alphabet:
            for third in alphabet:
                for fourth in alphabet:
                    kmer = first + second + third + fourth
                    lib.append(kmer)

    return lib

def reader(filename):
    """Read a DNA sequence from a FASTA file.

        Parameters:
            filename:
                str, name of the FASTA file containing the DNA sequence.

        Returns:
            A string containing the DNA sequence without the FASTA header.

        Effect:
            Empty lines and lines starting with '>' are ignored.
    """
    dna_seq = ''
    with open(filename) as f:
        for line in f:
            line = line.strip()

            if line == '':
                continue

            if line.startswith('>'):
                continue

            dna_seq += line

    return dna_seq

def parser(dna_seq):
    """Count the frequency of each 4-mer in a DNA sequence.

        Parameters:
            dna_seq:
                str, DNA sequence to be analyzed.

        Returns:
            A dictionary where:
                key:
                    A 4-mer string.
                value:
                    The number of times that 4-mer appears in dna_seq.

        Example:
            For dna_seq = 'AAAAA', the 4-mer 'AAAA' appears twice:
                positions 0-3: AAAA
                positions 1-4: AAAA
    """
    kmer_dict = {}
    for i in range(len(dna_seq)):
        if i+4 <= len(dna_seq):
            kmer = dna_seq[i:i+4]
            if kmer not in kmer_dict:
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1

    return kmer_dict

def output(lib, kmer_dict):
    """Print 4-mer frequencies in lexicographic order.

        Parameters:
            lib:
                list, all possible 4-mers in lexicographic order.

            kmer_dict:
                dict, observed 4-mer counts from the DNA sequence.

        Effect:
            Prints the frequency of each 4-mer in lib order on one line.
            If a 4-mer is not found in kmer_dict, its frequency is printed as 0.

        Returns:
            None
    """
    kmer_order_dict = {}
    for kmer in lib:
        if kmer in kmer_dict:
            freq = kmer_dict[kmer]
            kmer_order_dict[kmer] = freq
        else:
            kmer_order_dict[kmer] = 0

    freq_ls = []
    for value in kmer_order_dict.values():
        freq_ls.append(value)

    for value in freq_ls:
        print(value, end=' ')


def main():
    """Main function of this script
    """
    filename = argv[1]
    lib = kmer_library()
    dna_seq = reader(filename)
    kmer_dict = parser(dna_seq)
    output(lib, kmer_dict)

if __name__ == "__main__":
    main()