#!/usr/bin/env python3

"""
Author:

Description: this is a script to ...

Usage:
"""

"""To do:
implement the Needleman-Wunsch algorithm for global sequence alignment
with BlOSUM62 matrix and linear gap penalties.
1. create an initial matrix
2. fill the matrix cells
3. use the 'score' function to look up a value in BLOSUM62
4. a traceback strategy to obtained aligned sequences.
5. a function to calculate the % of identity between two sequence
6.assert statments to test the functionality
7.once the basic function is done, implement the functionality to use different values for end-gaps and regular gaps.

"""
from sys import argv
import unittest


# functions between here and __main__

# matrix source: https://github.com/biopython/biopython/blob/master/
#                Bio/Align/substitution_matrices/data/BLOSUM62
blosum ="""
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""


def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    returns: tuple of order and blosum_matrix, where
    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {} #{'A': 0, 'R': 1, 'N': 2, 'D': 3} #the header of the matrix
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24: #check if it is header?
            for idx, sym in enumerate(parts):
                #print(idx, sym)
                order[sym] = idx
        else: #the matrix part
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int, parts[1:]))) #convert numbers from string to int

    return order, blosum_matrix


BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()


def score(res1, res2):
    """Return score from BLOSUM62 matrix for two residues

    res1: string, an amino acid
    res2: string, an amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1] #coordinate (x)
    lookup2 = BLOSUM62_ORDER[res2] #coordinate (y)
    return BLOSUM62_MATRIX[lookup1][lookup2]
    #the score of two amino acid in the matrix
    #the order of two letter do not really matter because the matrix is symmetry

def initial_matrix_creation(seq1, seq2, gap_penalty):
    """Create and initialize the Needleman-Wunsch score matrix.

    Parameters:
        seq1: str, first amino acid sequence
        seq2: str, second amino acid sequence
        gap_penalty: int, penalty for opening/extending a gap

    Returns:
        initial_matrix, list of lists, initialized score matrix
    """
    rows = len(seq1) + 1 #8+1
    cols = len(seq2) + 1 #9+1
    initial_matrix = []

    for index in range(rows): #what if the length of two sequence is unequal?
                              #range(rows), [0] * cols
                              #number of rows * number of columns
        initial_matrix.append([0] * cols)
        #[0] * 8 = [0, 0, 0, 0, 0, 0, 0, 0]

    for i in range(rows): #the first num in each line is gap penalty
        initial_matrix[i][0] = i * gap_penalty

    for j in range(cols):
        initial_matrix[0][j] = j * gap_penalty

    #for row in initial_matrix:
        #print(row)

    return initial_matrix

def fill_matrix(seq1, seq2, gap_penalty):
    """Fill the Needleman-Wunsch score matrix.

    Parameters:
        seq1: str, first amino acid sequence
        seq2: str, second amino acid sequence
        gap_penalty: int, penalty for opening/extending a gap

    Returns:
        list of lists, filled score matrix
    """
    initial_matrix = initial_matrix_creation(seq1, seq2, gap_penalty)
    filled_matrix = initial_matrix
    for i in range(1, len(filled_matrix)):
        for j in range(1, len(filled_matrix[0])):
            diagonal = filled_matrix[i-1][j-1] + score(seq1[i - 1], seq2[j - 1])
            #why it is i-1 and j-1
            #because the matrix row [0] = 0, seq1[0] = matrix row[1]
            #and because your i is about matrix

            up = filled_matrix[i-1][j] + gap_penalty
            left = filled_matrix[i][j - 1] + gap_penalty

            filled_matrix[i][j] = max(diagonal, up, left)

    #for row in initial_matrix:
        #print(row)
    return filled_matrix

def traceback(seq1, seq2, filled_matrix, gap_penalty):
    #start with empty aligned seq
    aligned1 = ''
    aligned2 = ''

    #start from the end of both sequences
    i = len(seq1) #8
    j = len(seq2) #9

    while i > 0 or j > 0: #continue while there is still something left in seq1 or seq2
                          #stop when i == 0, j == 0
        if i > 0 and j > 0: #need to check if
            current_score = filled_matrix[i][j]
            diagonal_score = filled_matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1]) #i = len(seq1)
                                                                                           #seq1(i - 1) = last letter
            #diagonal_score = the upper diagnoal score + the score of the last letter
            if current_score == diagonal_score: #Did this cell come from the diagonal cell?
                aligned1 = seq1[i-1] + aligned1
                aligned2 = seq2[j-1] + aligned2
                i -= 1
                j -= 1
                continue #after making one move, skip the rest if statement, restart a new while loop after this move

        if i > 0: #if one sequence still have leftover
            current_score = filled_matrix[i][j]
            up_score = filled_matrix[i - 1][j] + gap_penalty

            if current_score == up_score:
                aligned1 = seq1[i - 1] + aligned1
                aligned2 = '-' + aligned2
                i -= 1
                continue

        if j > 0:
            current_score = filled_matrix[i][j]
            left_score = filled_matrix[i][j-1] + gap_penalty

            if current_score == left_score:
                aligned1 = '-' + aligned1 #cannot do like this
                                #traceback is assemble the sequence backward
                                #So every new character must be added to the front
                aligned2 = seq2[j - 1] + aligned2
                j -= 1
                continue

    return aligned1, aligned2

def identity_calculator(aligned1, aligned2):
    if len(aligned1) != len(aligned2):
        print(f'The length of two aligned sequences should be equal')
        exit()
    iden_index = 0
    for i in range(len(aligned1)):
        if aligned1[i] == '-' or aligned2[i] == '-': #excluding the gaps
            continue
        if aligned1[i] == aligned2[i]:
            iden_index += 1


    identity = (iden_index / len(aligned1)) * 100
    return identity


#store the arrow
#now you do the calculation twice but do not store anything


#now it has a preference that it prefer to go diagnoal
#how can i indicate the preference, or i can define which direction i prefer
#store the direction in a list, you can put your preference in the arguments
#there is a default preference, but you can also input yours.


def sequence_dict():
    """Dictionary of amino acid sequence

    Parameter: None

    return(s): a dictionary, {sequential No.1: sequence}
    where, sequential No.1: int, sequence: str
    """
    seq_dict = {}
    seq1 = "THISLINE"
    seq2 = "ISALIGNED"

    # seq3: GPA1_ARATH
    seq3 = ("MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIF"
            "KQIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEK"
            "LSEIGGRLDYPRLTKDIAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIP"
            "TKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAA"
            "ISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSFMLFLNKFDIFEKKVLDVPLNVCEWF"
            "RDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYRTTALDQKLVKKTFKLVDE"
            "TLRRRNLLEA")
    # seq4: GPA1 BRANA
    seq4 = (
        "MGLLCSRSRHHTEDTDENAQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQASS"
        "DKRKIIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDPAKYTLSSEN"
        "MAIGEKLSEIGARLDYPRLTKDLAEGIETLWNDPAIQETCSRGNELQVPDCTKYLMENLK"
        "RLSDVNYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHL"
        "FEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSIMLFLNKFDIF"
        "EKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYR"
        "TTALDQKLVKKTFKLVDETLRRRNLLEAGLL")
    seq_dict = {1:seq1, 2:seq2, 3:seq3, 4:seq4}
    return seq_dict

class TestAlignment(unittest.TestCase):
    def test_score(self):
        seq1 = "THISLINE"
        seq2 = "ISALIGNED"
        gap_penalty = int(argv[1])
        #print(fill_matrix(seq1, seq2, gap_penalty)[-1][-1])
        self.assertEqual(fill_matrix(seq1, seq2, gap_penalty)[-1][-1], 7)


def main():
    """Main function of the script
    """
    if len(argv) != 2:
        print("Usage: python <xxxx.py> <gap_penalty>")
        exit()
    seq_dict = sequence_dict()
    seq1 = seq_dict[1]
    seq2 = seq_dict[2]

    gap_penalty = int(argv[1])
    initial_matrix_creation(seq1, seq2, gap_penalty)
    filled_matrix = fill_matrix(seq1, seq2, gap_penalty)
    aligned1, aligned2 = traceback(seq1, seq2, filled_matrix, gap_penalty)
    identity = identity_calculator(aligned1, aligned2)

    #unit test
    align_test = TestAlignment()
    align_test.test_score()

if __name__ == "__main__":
    main()

