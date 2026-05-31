#!/usr/bin/env python3
"""
Name: Mingyang Wang
Student Number: 1067192
Description: Unit tests for the hamming_distance function.
Usage:python test_hamm.py
     (python test_hamm.py -v, if you want to check the process)

"""

import unittest
from Week1 import hamm


class TestHammingDistance(unittest.TestCase): #inherrit from class unittest
    def test_regulartests(self):
        self.assertEqual(hamm.hamming_distance('GAGCCTACTAACGGGAT',
                                               'CATCGTAATGACGGCCT'),
                         7)
        self.assertEqual(hamm.hamming_distance('AAAA', 'AAAT'), 1)
        self.assertEqual(hamm.hamming_distance('GATT', 'GAAA'), 2)


    def test_same_string(self):
        s = 'GATATATATATAG'
        t = 'GATATATATATAG'
        result = hamm.hamming_distance(s, t)
        self.assertEqual(0, result)

    def test_complete_diff(self):
        s = 'TTTTTTTTTTTTT'
        t = 'GGGGGGGGGGGGG'
        result = hamm.hamming_distance(s, t)
        self.assertEqual(len(s), result)

    def test_unequal_length(self):
        s = 'TTTT'
        t = 'TTT'
        with self.assertRaises(ValueError): hamm.hamming_distance(s, t)
        #assertRaises is an unittest method used to
        #test whether your code raises an expected error.

    def test_lowercase(self):
        s = 'TTTT'
        t = 'tttt'
        result = hamm.hamming_distance(s, t)
        self.assertEqual(0, result)

    def test_emptystring(self):
        s = ''
        t = ''
        with self.assertRaises(ValueError):
            hamm.hamming_distance(s, t)
        #means, run the hamm.hamming_distance(s,t), i expect a ValueError

if __name__ == "__main__":
    unittest.main() #this is the normal way to call the unittest
                    #python will find all the methods starts with test_


