#!/usr/bin/env python3
"""
Author:
Script to report a few sequence properties
Usage: python3 <script> <input_file>
    where script is the name of the python script
    and input_file the path to an input file containing a DNA sequence
"""

from sys import argv


class DNA:
    def __init__(self):
        self.seq = ''
        self.gc_content = 0
        self.re_seq = ''

    def parse_input(self, filename):
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if line == '':
                    continue
                self.seq += line


    def calculator(self):
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

    def __len__(self): #t defines what should happen when you call: len(object)
        return len(self.seq)

    def __str__(self):
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
    print("The length is the sequence is {} bp".format(len(dna)))
    print("The GC content of the sequence is {:.2f}%".format(dna.gc_content))
    print("This is the reverse complement of the sequence:\n{}".format(dna.re_seq))
    print(dna)

if __name__ == "__main__":
    main()