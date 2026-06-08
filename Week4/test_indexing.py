import indexing
import unittest

class Testindexing(unittest.TestCase):
    def test_regulartests(self):
        s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
        s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
        s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
        seq_ls = [s1, s2, s3]
        query_ls = ['TGCAACAT']
        k = 2

        #hashtable_maker function
        self.assertEqual(indexing.hashtable_maker(seq_ls, k)['TG'],
                         [(1,13), (2,7), (3,9)])
        self.assertEqual(indexing.hashtable_maker(seq_ls, k)['AA'],
                         [(2,19)])
        self.assertEqual(indexing.hashtable_maker(seq_ls, k)['AC'],
                         [(1, 9), (2,5), (2,11)])

        #query_kmer_maker function
        self.assertEqual(indexing.query_hashtable_maker
                         (query_ls, k)['TGCAACAT']['CA'], [2, 5])

        #match_table_maker()
        genome_kmer_dict = indexing.hashtable_maker(seq_ls, k)
        query_kmer_dict = indexing.query_hashtable_maker(query_ls, k)
        self.assertEqual(indexing.match_table_maker
                         (query_kmer_dict, genome_kmer_dict)['TGCAACAT']['TG'],
                         [(1,13,13), (2,7,7), (3,9,9)])

        #cluster_hits()
        kmer_posi_dict = (indexing.match_table_maker
                          (query_kmer_dict, genome_kmer_dict))
        self.assertEqual(len(indexing.cluster_hits
                             (kmer_posi_dict)['TGCAACAT'][2][7]), 4)

        #best_hit_dict()
        hits_dict = indexing.cluster_hits(kmer_posi_dict)
        self.assertEqual(indexing.best_hit_finder(hits_dict)
                         ['TGCAACAT']['seq_index'], 2)
        self.assertEqual(indexing.best_hit_finder(hits_dict)
                         ['TGCAACAT']['shift'], 7)
        self.assertEqual(len(indexing.best_hit_finder(hits_dict)
                             ['TGCAACAT']['hits']), 4)


    def test_edgecases(self):

        # Test k larger than sequence
        with self.assertRaises(ValueError):
            indexing.hashtable_maker(['AT'], 3)
            indexing.query_hashtable_maker(['AT'], 3)

        seq = 'A' * 101
        seq_ls = [seq]
        query_ls = ['TGCAACAT']
        k = 2
        genome_kmer_dict = indexing.hashtable_maker(seq_ls, k)
        query_kmer_dict = indexing.query_hashtable_maker(query_ls, k)
        kmer_posi_dict = indexing.match_table_maker(query_kmer_dict,
                                                    genome_kmer_dict)
        hits_dict = indexing.cluster_hits(kmer_posi_dict)
        best_hit_dict = indexing.best_hit_finder(hits_dict)
        with self.assertRaises(ValueError):
            indexing.matching_printer(best_hit_dict, seq_ls)

if __name__ == "__main__":
    unittest.main()