import hmm
import unittest

class TestHMM(unittest.TestCase):
    def test_regular_tests(self):
        seq_ls = ['-VFA-A-AGY', '-V---ANVDV', '-VE----VAH', '-VKG-----D',
                  '-VYS--TYES', '-FNA--NIPH', 'AIAGADNGAV']
        #the sequences from develop.fasta



        self.assertEqual(hmm.match_state_checker(seq_ls), [False, True,
                                                           True, True,
                                                           False, False, True,
                                                           True, True, True,
                                                    ])

        match_state = hmm.match_state_checker(seq_ls)

        self.assertEqual(hmm.reduced_alignment_maker(seq_ls, match_state)[0],
                         'VFA-AGY')

        transition_counts, emission_counts = (
            hmm.transitions_count(seq_ls,match_state))

        self.assertEqual(hmm.trans_count_normalize
          (transition_counts, emission_counts)[0]['Begin', 'M1'], 6/7)

        self.assertEqual(hmm.trans_count_normalize
            (transition_counts, emission_counts)[1]['M1', 'V'], 5/7)

    def test_unaligned_seq(self):
        unaligned_seq = ['AAAAA', 'AA', 'A', '']
        with self.assertRaises(ValueError):
            hmm.match_state_checker(unaligned_seq)

    def test_all_gaps_column(self):
        seq_ls = ['---', '---', '---']
        self.assertEqual(
            hmm.match_state_checker(seq_ls),
            [False, False, False])

    def test_empty_seq(self):
        empty_seq = ['', '', '']
        with self.assertRaises(ValueError):
            hmm.match_state_checker(empty_seq)

if __name__ == "__main__":
    unittest.main()
