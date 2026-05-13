#!/usr/bin/env python3
"""
Name: Mingyang Wang
Student Number: 1067192

Description:
    Read a DNA sequence from an input txt file, calculate the sequence length,
    GC content, and print the reverse complement sequence. The GC content is
    calculated as the percentage of G and C bases in the sequence.

Usage:
    python <script.py> <sequence_filename.txt>

Arguments:
    script.py: str
        Name of the script
    sequence_filename.txt: str
        Name of the input txt file containing a DNA sequence.

Output:
    Prints the sequence length, GC content, and reverse complement sequence to
    the terminal.
"""

from sys import argv


def parser(filename):
    """Read a DNA sequence from an input file.

    Args:
        filename: str, name of the input file containing a DNA sequence.

    Returns:
        str: The DNA sequence as a single string.
    """
    seq = ''
    with open(filename) as f:
        for line in f:
            line = line.strip()

            if line == '':
                continue

            seq += line

        return seq

def calculator(seq):
    """Calculate the GC content of a DNA sequence.

        Args:
            seq: str, a DNA sequence.

        Returns:
            float: The GC content of the sequence as a percentage.
    """
    count = 0
    for letter in seq:
        if letter == 'G' or letter == 'g' or letter == 'C' or letter == 'c':
        # equals to if (letter == 'G') or ('g') or ('C') or ('c'):
        #and 'g' will be treated as a non-empty string, which will be always True
            count += 1
    gc_content = (count / len(seq)) * 100
    return gc_content

def reverse_complement(seq):
    """Create the reverse complement sequence of a DNA sequence.

        Args:
            seq: str, a DNA sequence.

        Returns:
            str: The reverse complement sequence, where G is replaced by C,
            C by G, A by T, and T by A, and the final sequence is reversed.
    """
    com_seq = ''
    for letter in seq:
        if letter == 'G' or letter == 'g':
            com_seq += 'C'

        if letter == 'C' or letter == 'c':
            com_seq += 'G'

        if letter == 'A' or letter == 'a':
            com_seq += 'T'

        if letter == 'T' or letter == 't':
            com_seq += 'A'
    re_seq = com_seq[::-1]
    return re_seq

def output(seq, gc_content, re_seq):
    """Print the sequence length, GC content, and reverse complement sequence.

        Args:
            seq: str, the original DNA sequence.
            gc_content: float, the GC content of the sequence as a percentage.
            re_seq: str, the reverse complement sequence.

        Returns:
            None.
    """
    length = len(seq)
    print(f'The length of the sequence is {length} bp')
    print(f'The GC content of the sequence is {gc_content:.2f}%')
    print(f'This is the reverse complement of the sequence:\n{re_seq}')

def main():
    """Main function of this script"""
    if len(argv) != 2: #change the number based on the exam question
        print("Usage: python3 <xxxx.py> <xxxx.txt>")
        exit()

    filename = argv[1]
    seq = parser(filename)
    gc_content = calculator(seq)
    re_seq = reverse_complement(seq)
    output(seq, gc_content, re_seq)

if __name__ == "__main__":
    main()
