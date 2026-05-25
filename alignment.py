#!/usr/bin/env python3

"""
Author: Mingyang Wang

Description:
    Needleman-Wunsch global alignment using BLOSUM62.
    Supports different penalties for internal gaps and end gaps.

Usage:
    python <alignment.py> <gap_penalty> <end_gap_penalty>
    gap_penalty, int
    end_gap_penalty, int

Example:
    python alignment.py -4 0
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
            blosum_matrix.append(list(map(int, parts[1:])))
            #convert numbers from string to int

    return order, blosum_matrix


BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()


def score(res1, res2):
    """Return score from BLOSUM62 matrix for two residues

    Parameters:
        res1: string, an amino acid
        res2: string, an amino acid

    Return(s): BLOSUM62_MATRIX[lookup1][lookup2], int, substitution score
    """
    lookup1 = BLOSUM62_ORDER[res1] #coordinate (x)
    lookup2 = BLOSUM62_ORDER[res2] #coordinate (y)
    return BLOSUM62_MATRIX[lookup1][lookup2]
    #the score of two amino acid in the matrix
    #the order of two letter do not really matter because the matrix is symmetry

def initial_matrix_creation(seq1, seq2, end_gap_penalty):
    """Create and initialize the Needleman-Wunsch score matrix.

    Parameters:
        seq1: str, first amino acid sequence
        seq2: str, second amino acid sequence
        end_gap_penalty: int, end_penalty for opening/extending a gap

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
        initial_matrix[i][0] = i * end_gap_penalty

    for j in range(cols):
        initial_matrix[0][j] = j * end_gap_penalty

    #for row in initial_matrix:
        #print(row)

    return initial_matrix

def fill_matrix(seq1, seq2, gap_penalty, end_gap_penalty):
    """Fill the Needleman-Wunsch score matrix.

    Parameters:
        seq1: str, first amino acid sequence
        seq2: str, second amino acid sequence
        gap_penalty: int, penalty for opening/extending a gap
        end_gap_penalty: int, penalty for beginning and end gap

    Returns:
        list of lists, filled score matrix
    """

    filled_matrix = initial_matrix_creation(seq1, seq2, end_gap_penalty)

    rows = len(filled_matrix)
    cols = len(filled_matrix[0])

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = filled_matrix[i - 1][j - 1] + score(seq1[i - 1],
                                                           seq2[j - 1])
            #why it is i-1 and j-1
            #because the matrix row [0] = 0, seq1[0] = matrix row[1]
            #and because your i is about matrix

            # if gap reaches the end of seq2
            if j == cols - 1:
                up_gap = end_gap_penalty
            else:
                up_gap = gap_penalty

            # if gap reaches the end of seq1
            if i == rows - 1:
                left_gap = end_gap_penalty
            else:
                left_gap = gap_penalty

            up = filled_matrix[i - 1][j] + up_gap
            left = filled_matrix[i][j - 1] + left_gap

            filled_matrix[i][j] = max(diagonal, up, left)

    #for row in initial_matrix:
        #print(row)
    return filled_matrix

def traceback(seq1, seq2, filled_matrix, gap_penalty, end_gap_penalty):
    """Trace back through the filled score matrix to build the alignment.

    Parameters:
        seq1: str, first amino acid sequence
        seq2: str, second amino acid sequence
        filled_matrix: list of lists, completed Needleman-Wunsch score matrix
        gap_penalty: int, penalty for internal gaps
        end_gap_penalty: int, penalty for gaps at the beginning or end

    Returns:
        tuple of str, the two aligned sequences
    """

    aligned1 = ''
    aligned2 = ''

    i = len(seq1)
    j = len(seq2)

    rows = len(filled_matrix)
    cols = len(filled_matrix[0])

    while i > 0 or j > 0:
        current_score = filled_matrix[i][j]

        if j == 0 and i > 0:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = '-' + aligned2
            i -= 1
            continue

        if i == 0 and j > 0:
            aligned1 = '-' + aligned1
            aligned2 = seq2[j - 1] + aligned2
            j -= 1
            continue

        if i > 0 and j > 0:
            diagonal_score = filled_matrix[i - 1][j - 1] + score(
                seq1[i - 1], seq2[j - 1]
            )

            if current_score == diagonal_score:
                aligned1 = seq1[i - 1] + aligned1
                aligned2 = seq2[j - 1] + aligned2
                i -= 1
                j -= 1
                continue

        if i > 0:
            if j == cols - 1:
                up_gap = end_gap_penalty
            else:
                up_gap = gap_penalty

            up_score = filled_matrix[i - 1][j] + up_gap

            if current_score == up_score:
                aligned1 = seq1[i - 1] + aligned1
                aligned2 = '-' + aligned2
                i -= 1
                continue

        if j > 0:
            if i == rows - 1:
                left_gap = end_gap_penalty
            else:
                left_gap = gap_penalty

            left_score = filled_matrix[i][j - 1] + left_gap

            if current_score == left_score:
                aligned1 = '-' + aligned1
                aligned2 = seq2[j - 1] + aligned2
                j -= 1
                continue

        raise RuntimeError(
            f"Traceback stuck at i={i}, j={j}, current_score={current_score}"
        )

    return aligned1, aligned2

def identity_calculator(aligned1, aligned2):
    """Calculate percentage identity between two aligned sequences.

    Parameters:
        aligned1: str, first aligned sequence
        aligned2: str, second aligned sequence

    Returns:
        identity:float, percentage identity calculated as identical residues /
        alignment length * 100
    """
    if len(aligned1) != len(aligned2):
        print(f'The length of two aligned sequences should be equal')
        exit()
    identical = 0

    for i in range(len(aligned1)):
        if aligned1[i] == "-" or aligned2[i] == "-":
            continue

        if aligned1[i] == aligned2[i]:
            identical += 1

    identity = identical / len(aligned1) * 100
    return identity

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
    def test_regular(self):
        seq1 = "THISLINE"
        seq2 = "ISALIGNED"
        self.assertEqual(fill_matrix(seq1, seq2, -8, -8)[-1][-1], -4)
        self.assertEqual(fill_matrix(seq1, seq2, -4, -4)[-1][-1], 7)
        self.assertEqual(fill_matrix('BXBWANG', 'MINGYANG', -4, -4)
                         [-1][-1], 13)
        filled_matrix = fill_matrix(seq1, seq2, -4, -4)
        aligned1 = 'THIS-LI-NE-'
        aligned2 = '--ISALIGNED'
        self.assertEqual(traceback(seq1, seq2, filled_matrix, -4, -4),
                         (aligned1, aligned2))

def main():
    """Main function of the script
    """
    if len(argv) != 3:
        print("Usage: python <alignment.py> <gap_penalty> <end_gap_penalty>")
        exit()

    seq_dict = sequence_dict()
    seq1 = seq_dict[3]
    seq2 = seq_dict[4]

    gap_penalty = int(argv[1])
    end_gap_penalty = int(argv[2])
    filled_matrix = fill_matrix(seq1, seq2, gap_penalty, end_gap_penalty)
    aligned1, aligned2 = traceback(seq1, seq2, filled_matrix, gap_penalty,
                                   end_gap_penalty)
    print(aligned1)
    print(aligned2)

    identity = identity_calculator(aligned1, aligned2)
    print(f'{identity:.2f}')

    #unit test
    align_test = TestAlignment()
    align_test.test_regular()

if __name__ == "__main__":
    main()

#to do: 1 store the arrow
#2. understand the end_gap_penalty.
#redo without checking others.