#!/usr/bin/env python3

"""
Author:Mingyang Wang
Student No.: 1067192

Description: build a profile HMM trainer from an alignment
"""
# Import statements
import sys
import random

def fasta_reader(filename):
    seq_ls = []
    with open(filename) as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                continue

            if line == '':
                continue

            seq_ls.append(line)
    return seq_ls

def match_state_checker(seq_ls):
    #first checking if all sequence have the same length
    length = len(seq_ls[0])
    for seq in seq_ls:
        if length != len(seq):
            print(f'sequence length should be same')
            exit()
    #check which column is match
    sequence_length = len(seq_ls[0])
    col_num = len(seq_ls)
    match_state = [False] * sequence_length
    for i in range(sequence_length):
        count = 0
        for seq in seq_ls:
            if seq[i] != '-':
                count += 1

        if count >= col_num / 2:
            match_state[i] = True

    return match_state

def reduced_alignment_maker(seq_ls, match_state):
    unmatch_pos = []
    for i, state in enumerate(match_state):
        if state == False:
            unmatch_pos.append(i)

    redu_seq_ls = []
    for seq in seq_ls:
        seq_inls = list(seq)
        for index in sorted(unmatch_pos, reverse = True):
            #remove from the last
            #so the index of former part remained
            del seq_inls[index]

        redu_seq = ''.join(seq_inls)
        redu_seq_ls.append(redu_seq)

    return redu_seq_ls

def transitions_count():
    pass

# Background amino acid probabilities
pa = { 'A':0.074, 'C':0.025, 'D':0.054, 'E':0.054, 'F':0.047, 'G':0.074,
    'H':0.026, 'I':0.068, 'L':0.099, 'K':0.058, 'M':0.025, 'N':0.045,
    'P':0.039, 'Q':0.034, 'R':0.052, 'S':0.057, 'T':0.051, 'V':0.073,
    'W':0.013, 'Y':0.034 }

def sample(events):
    """Return a key from dict based on the probabilities

    :events: dict of {string: probability}, probabilities can also be weights.
    :return: string
    """

    pick = random.choices(list(events.keys()),list(events.values()))[0]
    #the keys, which is the residue letter, is randomly picked based on its probability as weight
    return pick

def transition_counts(alignment,is_match_state):
    """Count the transitions between states in the HMM

    :param alignment: list; protein sequence alignments.
    :param is_match_state: binary list; indication of match state for each
                alignment position
    :return: dictionary; key - string, transitions types. value - list of
                integers, length number of match states + 1. Each entry in the
                list represents the number of transitions between the states.

    This function counts the transitions sequence by sequence. It starts in a
    match state at position 0.
    """

    names = ['MM', 'MI', 'MD', 'IM', 'II', 'ID', 'DM', 'DI', 'DD']
    nm=sum(is_match_state)+1
    trans_dic = {}
    for name in names:
        trans_dic[name] = [0]*nm

    #continue

    return trans_dic

def emission_counts(alignment,is_match_state):

    """
    Count the emissions for each match state

    :param alignment: list; protein sequence alignments.
    :param is_match_state: binary list; indication of match state for each
                alignment position
    :return: list of dictionaries; key - amino acid. value - count.

    This function counts the emissions sequence by sequence.
    """

    nm=sum(is_match_state)+1
    em_dic = []
    for i in range(0,nm):
        em_dic.append({})
        for j in pa.keys():
            em_dic[i][j] = 0

    #continue

    return em_dic

def main():
    infile = 'develop.fasta'
    seq_ls = fasta_reader(infile)
    match_state = match_state_checker(seq_ls)
    print(match_state)
    redu_seq_ls = reduced_alignment_maker(seq_ls, match_state)
    print(redu_seq_ls)
    #the main should produce all the results for the report

if __name__ == "__main__":
    main()
    # implement main code here

    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file)
