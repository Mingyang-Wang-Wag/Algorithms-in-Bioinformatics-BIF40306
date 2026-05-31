#!/usr/bin/env python3
"""
Name: Mingyang Wang
Student Number: 1067192
Description:
Usage:

"""


from sys import argv


def dummy(fasta_file): #######DELETE this part before submit
    """Parse a FASTQ file into xxxx

    fasta_file: a fasta file, xxxxxx
    (if create anything, say it)
    returns: a dictionary, key: str, xxxx; value: str, xxxx
    (if you do not have a return, then do not call this)
    """

    pass

def main():
    """Main function of this script"""
    if len(argv) != 4: #change the number
        print("Usage: python3 xxxx")
        exit()

    filename_input = argv[1]
    other_input = argv[2]
    #parameter_input = argv[3]



if __name__ == "__main__":
    main()
