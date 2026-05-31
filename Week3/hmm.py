#!/usr/bin/env python3

"""
Author: Mingyang Wang
Student No.: 1067192

Description: build a profile HMM trainer from an alignment
"""
# Import statements
import sys
import random

def fasta_reader(filename):
    seq_ls = [] #not for test_large.fasta
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue

            if line.startswith('>'):
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

def transitions_count(seq_ls, match_state):
    """
    input, the sequence and the match state, reduce sequence

    output, the state count
    emission
    """
    transition_counts = []
    emission_counts = []
    for seq in seq_ls:
        match_count = 0
        for i, state in enumerate(match_state):
            #three section is similar
            if i == 0: #section 1, the beginning and first state
                if seq[i] == '-' and state == False: #no emission in this state
                    trans_state = ('Begin', f'')

                    found = False #this is the flag
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])


                elif seq[i] != '-' and state == False: #Insertion
                    trans_state = ('Begin', f'I{match_count}')
                    #emission adding
                    residue = seq[i]
                    current_state = trans_state[1]
                    emission_found = False
                    for emission_item in emission_counts:
                        if emission_item[0] == (current_state, residue):
                            emission_item[1] += 1
                            emission_found = True
                            break
                    if emission_found == False:
                        emission_counts.append([(current_state, residue), 1])

                    #transition adding
                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

                elif seq[i] == '-' and state == True: #Deletion #no emission
                    match_count += 1
                    trans_state = ('Begin', f'D{match_count}')

                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

                elif seq[i] != '-' and state == True:
                    match_count += 1
                    trans_state = ('Begin', f'M{match_count}')
                    # emission adding
                    residue = seq[i]
                    current_state = trans_state[1]
                    emission_found = False
                    for emission_item in emission_counts:
                        if emission_item[0] == (current_state, residue):
                            emission_item[1] += 1
                            emission_found = True
                            break
                    if emission_found == False:
                        emission_counts.append([(current_state, residue), 1])

                    # transition adding
                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

            elif i < len(match_state) - 1:  # section 2, the state in between
                if seq[i] == '-' and state == False: #no emission
                    pre_state = transition_counts[i - 1][0][1]
                    trans_state = (f'{pre_state}', f'')

                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

                elif seq[i] != '-' and state == False:
                    pre_state = transition_counts[i - 1][0][1]
                    trans_state = (f'{pre_state}', f'I{match_count}')
                    # emission adding
                    residue = seq[i]
                    current_state = trans_state[1]
                    emission_found = False
                    for emission_item in emission_counts:
                        if emission_item[0] == (current_state, residue):
                            emission_item[1] += 1
                            emission_found = True
                            break
                    if emission_found == False:
                        emission_counts.append([(current_state, residue), 1])

                    # transition adding
                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

                elif seq[i] == '-' and state == True: #no emission
                    match_count += 1
                    pre_state = transition_counts[i - 1][0][1]
                    trans_state = (f'{pre_state}', f'D{match_count}')

                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

                elif seq[i] != '-' and state == True:
                    match_count += 1
                    pre_state = transition_counts[i - 1][0][1]
                    trans_state = (f'{pre_state}', f'M{match_count}')
                    # emission adding
                    residue = seq[i]
                    current_state = trans_state[1]
                    emission_found = False
                    for emission_item in emission_counts:
                        if emission_item[0] == (current_state, residue):
                            emission_item[1] += 1
                            emission_found = True
                            break
                    if emission_found == False:
                        emission_counts.append([(current_state, residue), 1])

                    # transition adding
                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

            else: #section 3, the last state and end state
                if seq[i] == '-' and state == False: # no emission
                    pre_state = transition_counts[i - 1][0][1]
                    trans_state = (f'{pre_state}', f'')
                    end_state = (f'{trans_state[1]}', f'End')
                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

                    for item in transition_counts:
                        if item[0] == end_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([end_state, 1])

                elif seq[i] != '-' and state == False:
                    pre_state = transition_counts[i - 1][0][1]
                    trans_state = (f'{pre_state}', f'I{match_count}')
                    end_state = (f'{trans_state[1]}', f'End')
                    # emission adding
                    residue = seq[i]
                    current_state = trans_state[1]
                    emission_found = False
                    for emission_item in emission_counts:
                        if emission_item[0] == (current_state, residue):
                            emission_item[1] += 1
                            emission_found = True
                            break
                    if emission_found == False:
                        emission_counts.append([(current_state, residue), 1])

                    # transition adding
                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

                    for item in transition_counts:
                        if item[0] == end_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([end_state, 1])

                elif seq[i] == '-' and state == True: #no emission
                    match_count += 1
                    pre_state = transition_counts[i - 1][0][1]
                    trans_state = (f'{pre_state}', f'D{match_count}')
                    end_state = (f'{trans_state[1]}', f'End')
                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

                    for item in transition_counts:
                        if item[0] == end_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([end_state, 1])

                elif seq[i] != '-' and state == True:
                    match_count += 1
                    pre_state = transition_counts[i - 1][0][1]
                    trans_state = (f'{pre_state}', f'M{match_count}')
                    end_state = (f'{trans_state[1]}', f'End')
                    # emission adding
                    residue = seq[i]
                    current_state = trans_state[1]
                    emission_found = False
                    for emission_item in emission_counts:
                        if emission_item[0] == (current_state, residue):
                            emission_item[1] += 1
                            emission_found = True
                            break
                    if emission_found == False:
                        emission_counts.append([(current_state, residue), 1])

                    # transition adding
                    found = False
                    for item in transition_counts:
                        if item[0] == trans_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([trans_state, 1])

                    for item in transition_counts:
                        if item[0] == end_state:
                            item[1] += 1
                            found = True
                            break

                    if found == False:
                        transition_counts.append([end_state, 1])

    return transition_counts, emission_counts

def trans_count_normalize(transition_counts, emission_counts):
    #for transition
    trans_total_dict = {}
    for trans_counts in transition_counts:
        trans_name = trans_counts[0][0]
        count = trans_counts[1]
        if trans_name == '':
            continue

        if trans_name not in trans_total_dict:
            trans_total_dict[trans_name] = count
        else:
            trans_total_dict[trans_name] += count

    trans_norm_dict = {}
    for trans_state in trans_total_dict:
        for trans_counts in transition_counts:
            if trans_state == trans_counts[0][0]:
                branch = trans_counts[0]
                count = trans_counts[1]
                total_count = trans_total_dict[trans_state]
                ratio = round((count / total_count), 2)
                trans_norm_dict[branch] = ratio

    #for emission
    emission_total_count = {}
    for counts in emission_counts:
        emission_name = counts[0][0]
        count = counts[1]

        if emission_name not in emission_total_count:
            emission_total_count[emission_name] = count
        else:
            emission_total_count[emission_name] += count

    emission_norm_dict = {}
    for state in emission_total_count:
        for counts in emission_counts:
            if state == counts[0][0]:
                branch = counts[0]
                count = counts[1]
                total_count = emission_total_count[state]
                ratio = round((count / total_count), 2)
                emission_norm_dict[branch] = ratio


    return trans_norm_dict, emission_norm_dict






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

def fake_sequence_generator(trans_norm_dict, emission_norm_dict):
    current_state = 'Begin'


    #look through transition probabilities and collect only transitions that start from current state

    next_possible_state = {} #based on the previous steps, only store the possible steps, empty every loop.
    pass



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
    print(em_dic)
    return em_dic

def main():
    """
    #the main should produce all the results for the report


    """
    infile = 'develop.fasta'
    seq_ls = fasta_reader(infile)
    #print(seq_ls) #['-VFA-A-AGY', '-V---ANVDV', '-VE----VAH', '-VKG-----D', '-VYS--TYES', '-FNA--NIPH', 'AIAGADNGAV']
    match_state = match_state_checker(seq_ls)
    #print(match_state) #[False, True, True, True, False, False, True, True, True, True]
    redu_seq_ls = reduced_alignment_maker(seq_ls, match_state)
    #print(redu_seq_ls) #['VFA-AGY', 'V--NVDV', 'VE--VAH', 'VKG---D', 'VYSTYES', 'FNANIPH', 'IAGNGAV']
    transition_counts, emission_counts = transitions_count(seq_ls, match_state)
    #print(transition_counts)
    #print(emission_counts)
    trans_norm_dict, emission_norm_dict = trans_count_normalize(transition_counts, emission_counts)
    print(trans_norm_dict)
    print(emission_norm_dict)
    fake_sequence_generator(trans_norm_dict, emission_norm_dict)

if __name__ == "__main__":
    main()
    # implement main code here

    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file)
