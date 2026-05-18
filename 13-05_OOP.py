#!/usr/bin/env python3
"""
Name: Mingyang Wang
Student Number: 1067192

Description:
    Read a DNA sequence from an input file and report several sequence
    properties, including the sequence length, GC content, and reverse
    complement.

    The script stores the DNA sequence in a DNA object, calculates the GC
    content as a percentage, generates the reverse complement sequence, and
    prints the results to the terminal.

Usage:
    python3 <script_name.py> <input_file>

    where:
        script_name.py:
            Name of this Python script.

        input_file:
            str, name of the input file, containing a DNA sequence.
            The file is expected to contain only sequence lines. Empty lines are
            ignored.
"""

from sys import argv


class DNA:
    def __init__(self):
        """Initialize a DNA object with an empty sequence, zero GC content,
           and an empty reverse complement sequence.
        """
        self.seq = ''
        self.gc_content = 0
        self.re_seq = ''

    def parse_input(self, filename):
        """Read a DNA sequence from an input file.

            Parameters:
                filename: str, name of the input file containing sequence

            Effect:
                Stores the sequence in self.seq. Empty lines are ignored.

            Returns:
                None
        """
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if line == '':
                    continue
                self.seq += line

    def calculator(self):
        """Calculate the GC content of the DNA sequence.

            Effect:
                Counts the number of G, C characters in self.seq and stores
                the GC content percentage in self.gc_content.

            Returns:
                None

            Assumption:
                self.seq contains a non-empty DNA sequence.
        """
        count = 0
        if self.seq == '':
            print('seq is empty')
            exit()
        else:
            for letter in self.seq:
                if letter == 'G' or letter == 'g' or letter == 'C' or letter == 'c':
                # equals to if (letter == 'G') or ('g') or ('C') or ('c'):
                # and 'g' will be treated as a non-empty string, which will be always True
                    count += 1

        self.gc_content = (count / len(self.seq)) * 100


    def reverse_complement(self):
        """Generate the reverse complement of the DNA sequence.

            Effect:
                Converts A to T, T to A, G to C, and C to G, then reverses the
                complemented sequence. The result is stored in self.re_seq.

            Returns:
                None

            Assumption:
                self.seq contains only DNA bases: A, T, G, and C.
        """
        com_seq = ''
        for letter in self.seq:
            if letter == 'G' or letter == 'g':
                com_seq += 'C'

            if letter == 'C' or letter == 'c':
                com_seq += 'G'

            if letter == 'A' or letter == 'a':
                com_seq += 'T'

            if letter == 'T' or letter == 't':
                com_seq += 'A'
        self.re_seq = com_seq[::-1]

    def __len__(self):# defines what should happen when you call: len(object)
        """Return the length of the DNA sequence.

        Returns:
            An integer representing the number of bases in self.seq.
        """

        return len(self.seq)

    def __str__(self):
        """Return a string representation of the DNA object.

        Returns:
            A string containing the DNA sequence and its GC content.
        """
        return f'{self.seq} {self.gc_content}'

def main():
    """Main code of the script
    """
    # step 1, get DNA sequence from the file
    dna = DNA()
    dna.parse_input(argv[1])
    dna.calculator()
    dna.reverse_complement()

    # step 2, print length, GC content, reverse complement
    print(f'The length is the sequence is {len(dna)} bp')
    print(f'The GC content of the sequence is {dna.gc_content:.2f}')
    print(f'The reverse complement of the sequence:\n{dna.re_seq}')
    #print(dna)

if __name__ == "__main__":
    main()